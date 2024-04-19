module calc_time_dev
  use cudafor
  use mod_globals, only : accuracy, id_hybrid, id_muscl, id_visc, nx, ny, nt, np, blocksE, blocksF, threadsE, threadsF, blocks, threads
  use calc_physical_quantities
  use calc_steps
  use calc_flux, only : calc_E, calc_F
  use calc_hybrid
  use calc_visc
  use set
  use print
  implicit none
  interface calc_EF
    module procedure calc_EF_basic, calc_EF_hybrid
  end interface
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th
  end interface
contains
  subroutine calc_EF_basic(id_hybrid,Q,T,E,F)
    integer(kind=2), intent(in) :: id_hybrid
    real(8), intent(in), device :: Q(nx,ny,4)
    real(8), intent(inout), device :: T(nx,ny)
    real(8), intent(out), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(out), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), dimension(nx,ny), device :: rho, u, v, p
    integer stat
    call calc_quantities(Q,rho,u,v,p,T)
    call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,p,E)
    call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,p,F)
    stat = cudaDeviceSynchronize()
        
    if (id_visc == 1) then 
      call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E)
      call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EF_basic

  subroutine calc_EF_hybrid(id_hybrid,Q,T,E_hybrid,F_hybrid)
    integer(kind=4), intent(in) :: id_hybrid
    real(8), intent(in), device :: Q(nx,ny,4)
    real(8), intent(inout), device :: T(nx,ny)
    real(8), intent(out), device :: E_hybrid(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(out), device :: F_hybrid(nx-accuracy,ny-accuracy+1,4)
    real(8), dimension(nx,ny), device :: rho, u, v, p, energy, fd
    real(8), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_keep, E_upwind
    real(8), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_keep, F_upwind
    integer stat
    integer(kind=2) :: id_muscl1
    call calc_quantities(Q,rho,u,v,p,T)
    call calc_E<<<blocksE,threadsE>>>(id_muscl1,rho,u,v,p,E_keep)
    call calc_F<<<blocksF,threadsF>>>(id_muscl1,rho,u,v,p,F_keep)
    call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,p,E_upwind)
    call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,p,F_upwind)
    stat = cudaDeviceSynchronize()

    call calc_Ducros<<<blocks,threads>>>(u,v,fd)
    stat = cudaDeviceSynchronize()
    energy(:,:) = Q(:,:,4)
    call calc_E_hybrid<<<blocksE,threadsE>>>(rho,energy,fd,E_keep,E_upwind,E_hybrid)
    call calc_F_hybrid<<<blocksF,threadsF>>>(rho,energy,fd,F_keep,F_upwind,F_hybrid)
    stat = cudaDeviceSynchronize()

    if (id_visc == 1) then 
      call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E_hybrid)
      call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F_hybrid)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EF_hybrid
  
  subroutine RungeKutta_3rd(id_RungeKutta,T0,Q)
    integer(kind=2), intent(in) :: id_RungeKutta
    real(8), intent(inout) :: Q(nx,ny,4)
    real(8), intent(inout) :: T0(nx,ny)
    integer t1, t2, itr
    real(8), dimension(nx,ny,4), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny), device :: T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,4)

    ! print initial condition
    call print_vtk(0,Q,T0)
    ! copy on GPU
    Q_d = Q
    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_EF(id_hybrid,Q_d,T,E,F)
        call calc_step1(E,F,Q_d,Q2)
        call set_bc(Q2)
        
        call calc_EF(id_hybrid,Q2,T,E,F)
        call calc_step2(E,F,Q_d,Q2,Q3)
        call set_bc(Q3)

        call calc_EF(id_hybrid,Q3,T,E,F)
        call calc_step3(E,F,Q3,Q_d)
        call set_bc(Q_d)
      enddo
      Q = Q_d
      T0 = T
      call print_vtk(t2,Q,T0)
    enddo
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,T0,Q)
    integer(kind=4), intent(in) :: id_RungeKutta
    real(8), intent(inout) :: Q(nx,ny,4)
    real(8), intent(inout) :: T0(nx,ny)
    integer t1, t2, itr
    real(8), dimension(nx,ny,4), device :: Q_d, Qs
    real(8), dimension(nx,ny), device :: T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), device :: Rs(nx-accuracy,ny-accuracy,4)

    ! print initial condition
    call print_vtk(0,Q,T0)
    ! copy on GPU
    Q_d = Q
    T = T0
    Rs(:,:,:) = 0.d0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_EF(id_hybrid,Q_d,T,E,F)
        call calc_step(0.5d0,1.d0,E,F,Rs,Q_d,Qs)
        call set_bc(Qs)
        
        call calc_EF(id_hybrid,Qs,T,E,F)
        call calc_step(0.5d0,2.d0,E,F,Rs,Q_d,Qs)
        call set_bc(Qs)

        call calc_EF(id_hybrid,Qs,T,E,F)
        call calc_step(1.d0,2.d0,E,F,Rs,Q_d,Qs)
        call set_bc(Qs)

        call calc_EF(id_hybrid,Qs,T,E,F)
        call calc_step4(E,F,Rs,Q_d)
        call set_bc(Q_d)
      enddo
      Q = Q_d
      T0 = T
      call print_vtk(t2,Q,T0)
    enddo
  end subroutine RungeKutta_4th
end module calc_time_dev

