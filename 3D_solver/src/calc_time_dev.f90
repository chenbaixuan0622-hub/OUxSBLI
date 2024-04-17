module calc_time_dev
  use cudafor
  use mod_globals, only : accuracy, id_hybrid, id_muscl, id_visc, id_turbulence, nx, ny, nz, nt, np, blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_steps
  use calc_flux
  use calc_hybrid
  use calc_visc
  use calc_les
  use set
  use print
  implicit none
  interface calc_EFG
    module procedure calc_EFG_basic, calc_EFG_hybrid
  end interface
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th
  end interface
contains
  subroutine calc_EFG_basic(id_hybrid,Q,T,mut,E,F,G)
    integer(kind=2), intent(in) :: id_hybrid
    real(8), intent(in), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), dimension(nx,ny,nz), device :: T, mut
    real(8), intent(out), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p
    integer stat
    call calc_quantities(Q,rho,u,v,w,p,T)

    call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,w,p,E)
    call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,w,p,F)
    call calc_G<<<blocksG,threadsG>>>(id_muscl,rho,u,v,w,p,G)
    !print *, trim(cudaGetErrorString(cudaGetLastError()))
    stat = cudaDeviceSynchronize()

    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads>>>(rho,u,v,w,mut)
      stat = cudaDeviceSynchronize()
      call set_bc_mut(mut)
    endif

    if (id_visc == 1 .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE>>>(u,v,w,T,mut,E)
      call calc_Fv<<<blocksF,threadsF>>>(u,v,w,T,mut,F)
      call calc_Gv<<<blocksG,threadsG>>>(u,v,w,T,mut,G)
    endif
    !print *, trim(cudaGetErrorString(cudaGetLastError()))
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_basic
  
  subroutine calc_EFG_hybrid(id_hybrid,Q,T,mut,E_hybrid,F_hybrid,G_hybrid)
    integer(kind=4), intent(in) :: id_hybrid
    real(8), intent(in), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), dimension(nx,ny,nz), device :: T, mut
    real(8), intent(out), device :: E_hybrid(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device :: F_hybrid(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device :: G_hybrid(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, energy
    real(8), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_keep, E_upwind
    real(8), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_keep, F_upwind
    real(8), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_keep, G_upwind
    integer stat
    integer(kind=2) :: id_muscl1
    call calc_quantities(Q,rho,u,v,w,p,T)

    call calc_E<<<blocksE,threadsE>>>(id_muscl1,rho,u,v,w,p,E_keep)
    call calc_F<<<blocksF,threadsF>>>(id_muscl1,rho,u,v,w,p,F_keep)
    call calc_G<<<blocksG,threadsG>>>(id_muscl1,rho,u,v,w,p,G_keep)
    call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,w,p,E_upwind)
    call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,w,p,F_upwind)
    call calc_G<<<blocksG,threadsG>>>(id_muscl,rho,u,v,w,p,G_upwind)
    stat = cudaDeviceSynchronize()

    energy(:,:,:) = Q(:,:,:,5)
    call calc_E_hybrid<<<blocksE,threadsE>>>(rho,u,v,w,energy,E_keep,E_upwind,E_hybrid)
    call calc_F_hybrid<<<blocksF,threadsF>>>(rho,u,v,w,energy,F_keep,F_upwind,F_hybrid)
    call calc_G_hybrid<<<blocksG,threadsG>>>(rho,u,v,w,energy,G_keep,G_upwind,G_hybrid)
    stat = cudaDeviceSynchronize()
    !print *, trim(cudaGetErrorString(cudaGetLastError()))

    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads>>>(rho,u,v,w,mut)
      stat = cudaDeviceSynchronize()
      call set_bc_mut(mut)
    endif

    if (id_visc == 1 .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE>>>(u,v,w,T,mut,E_hybrid)
      call calc_Fv<<<blocksF,threadsF>>>(u,v,w,T,mut,F_hybrid)
      call calc_Gv<<<blocksG,threadsG>>>(u,v,w,T,mut,G_hybrid)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_hybrid

  subroutine RungeKutta_3rd(id_RungeKutta,T0,Q)
    integer(kind=2), intent(in) :: id_RungeKutta
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    real(8), intent(inout) :: T0(nx,ny,nz)
    integer t1, t2, itr, stat
    integer(kind=2) :: id_muscl
    real(8), dimension(nx,ny,nz) :: mut_cpu = 0.d0
    real(8), dimension(nx,ny,nz,5), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny,nz), device :: T, mut
    real(8), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    ! for plot
    real(8) ke0, entropy0

    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,Q,T0,ke0,entropy0)
    else
      call print_vtk(0,Q,T0,ke0,entropy0,mut_cpu)
    endif
    ! copy on GPU
    Q_d = Q
    T = T0
    mut = mut_cpu
    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(id_hybrid,Q_d,T,mut,E,F,G)
        call calc_step1(E,F,G,Q_d,Q2)
        call set_bc(id_accuracy,Q2,T)

        call calc_EFG(id_hybrid,Q2,T,mut,E,F,G)
        call calc_step2(E,F,G,Q_d,Q2,Q3)
        call set_bc(id_accuracy,Q3,T)

        call calc_EFG(id_hybrid,Q3,T,mut,E,F,G)
        call calc_step3(E,F,G,Q3,Q_d)
        call set_bc(id_accuracy,Q_d,T)
      enddo
      Q = Q_d
      T0 = T
      mut_cpu = mut
      if (id_turbulence == 0) then
        call print_vtk(t2,Q,T0,ke0,entropy0)
      else
        call print_vtk(t2,Q,T0,ke0,entropy0,mut_cpu) 
      endif
    enddo
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,T0,Q)
    integer(kind=4), intent(in) :: id_RungeKutta
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    real(8), intent(inout) :: T0(nx,ny,nz)
    integer t1, t2, itr, stat
    integer(kind=2) :: id_muscl
    real(8), dimension(nx,ny,nz) :: mut_cpu = 0.d0
    real(8), dimension(nx,ny,nz,5), device :: Q_d, Qs
    real(8), dimension(nx,ny,nz), device :: T, mut
    real(8), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), device :: Rs(nx-accuracy,ny-accuracy,nz-accuracy,5)
    ! for plot
    real(8) ke0, entropy0

    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,Q,T0,ke0,entropy0)
    else
      call print_vtk(0,Q,T0,ke0,entropy0,mut_cpu)
    endif
    ! copy on GPU
    Q_d = Q
    T = T0
    mut = mut_cpu
    Rs(:,:,:,:) = 0.d0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(id_hybrid,Q_d,T,mut,E,F,G)
        call calc_step(0.5d0,1.d0,E,F,G,Rs,Q_d,Qs)
        call set_bc(id_accuracy,Qs,T)

        call calc_EFG(id_hybrid,Qs,T,mut,E,F,G)
        call calc_step(0.5d0,2.d0,E,F,G,Rs,Q_d,Qs)
        call set_bc(id_accuracy,Qs,T)

        call calc_EFG(id_hybrid,Qs,T,mut,E,F,G)
        call calc_step(1.d0,2.d0,E,F,G,Rs,Q_d,Qs)
        call set_bc(id_accuracy,Qs,T)

        call calc_EFG(id_hybrid,Qs,T,mut,E,F,G)
        call calc_step4(E,F,G,Rs,Q_d)
        call set_bc(id_accuracy,Q_d,T)
      enddo
      Q = Q_d
      T0 = T
      mut_cpu = mut
      if (id_turbulence == 0) then
        call print_vtk(t2,Q,T0,ke0,entropy0)
      else
        call print_vtk(t2,Q,T0,ke0,entropy0,mut_cpu)
      endif
    enddo
  end subroutine RungeKutta_4th
end module calc_time_dev

