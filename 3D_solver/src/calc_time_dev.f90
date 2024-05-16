module calc_time_dev
  use cudafor
  use mod_globals, only : accuracy, id_hybrid, id_muscl, id_visc, id_turbulence, nx, ny, nz, nt, np, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
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
    module procedure RungeKutta_3rd, RungeKutta_4th, RungeKutta_10th
  end interface
contains
  subroutine calc_EFG_basic(id_hybrid,dx,dy,Jacobian,QJ,mut,E,F,G)
    integer(kind=2), intent(in)                         :: id_hybrid
    real(8), intent(in), dimension(nx-1), device        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device        :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device       :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(inout), dimension(nx,ny,nz), device :: mut
    real(8), intent(out), device                        :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                        :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                        :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, T
    integer stat
    call calc_quantities(Jacobian,QJ,rho,u,v,w,p,T)

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
      call calc_Ev<<<blocksE,threadsE>>>(dx,dy,rho,u,v,w,T,p,mut,E)
      call calc_Fv<<<blocksF,threadsF>>>(dy,dx,rho,u,v,w,T,p,mut,F)
      call calc_Gv<<<blocksG,threadsG>>>(dx,dy,rho,u,v,w,T,p,mut,G)
    endif
    !print *, trim(cudaGetErrorString(cudaGetLastError()))
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_basic
  
  subroutine calc_EFG_hybrid(id_hybrid,dx,dy,Jacobian,QJ,mut,E_hybrid,F_hybrid,G_hybrid)
    integer(kind=4), intent(in)                         :: id_hybrid
    real(8), intent(in), dimension(nx-1), device        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device        :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device       :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(inout), dimension(nx,ny,nz), device :: mut
    real(8), intent(out), device                        :: E_hybrid(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                        :: F_hybrid(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                        :: G_hybrid(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device                                :: rho, u, v, w, p, T, energy
    real(8), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_keep, E_upwind
    real(8), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_keep, F_upwind
    real(8), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_keep, G_upwind
    integer stat
    integer(kind=2) :: id_muscl1
    call calc_quantities(Jacobian,QJ,rho,u,v,w,p,T)

    call calc_E<<<blocksE,threadsE>>>(id_muscl1,rho,u,v,w,p,E_keep)
    call calc_F<<<blocksF,threadsF>>>(id_muscl1,rho,u,v,w,p,F_keep)
    call calc_G<<<blocksG,threadsG>>>(id_muscl1,rho,u,v,w,p,G_keep)
    call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,w,p,E_upwind)
    call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,w,p,F_upwind)
    call calc_G<<<blocksG,threadsG>>>(id_muscl,rho,u,v,w,p,G_upwind)
    stat = cudaDeviceSynchronize()

    !!!!not correct !!!!!!!!!!!!!
    energy(:,:,:) = QJ(:,:,:,5)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      call calc_Ev<<<blocksE,threadsE>>>(dx,dy,rho,u,v,w,T,p,mut,E_hybrid)
      call calc_Fv<<<blocksF,threadsF>>>(dy,dx,rho,u,v,w,T,p,mut,F_hybrid)
      call calc_Gv<<<blocksG,threadsG>>>(dx,dy,rho,u,v,w,T,p,mut,G_hybrid)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_hybrid

  subroutine RungeKutta_3rd(id_RungeKutta,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q,Vin)
    integer(kind=2), intent(in) :: id_RungeKutta
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(8), intent(in)         :: Vin(ny,2)
    integer i, j, t1, t2, itr, stat
    real(8), dimension(nx,ny,nz) :: mut_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), dimension(nx,ny,nz,5), device  :: QJ, QJ2, QJ3
    real(8), dimension(nx,ny,nz), device    :: mut
    real(8), device                         :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device                         :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device                         :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx-1), device        :: dx, xix  ! xix = 1 / dx
    real(8), dimension(ny-1), device        :: dy, etay ! etay = 1 / dy
    real(8), dimension(nx,ny), device       :: Jacobian
    ! for plot
    real(8) ke0, entropy0

    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
    enddo;enddo

    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
    else
      call print_vtk(0,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
    endif

    ! copy on GPU
    QJ = Q
    mut = 0.d0
    xix = xix_cpu
    etay = etay_cpu
    dx = dx_cpu
    dy = dy_cpu
    Jacobian = Jacobian_cpu
    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJ,mut,E,F,G)
        call calc_step(1.d0,0.d0,dx,dy,E,F,G,QJ,QJ2)
        call set_bc(id_accuracy,Jacobian,QJ2)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJ2,mut,E,F,G)
        call calc_step2(0.75d0,0.25d0,0.25d0,1.d0,dx,dy,E,F,G,QJ,QJ2,QJ3)
        call set_bc(id_accuracy,Jacobian,QJ3)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJ3,mut,E,F,G)
        call calc_step3(dx,dy,E,F,G,QJ3,QJ)
        call set_bc(id_accuracy,Jacobian,QJ)
      enddo
      Q = QJ
      mut_cpu = mut
      if (id_turbulence == 0) then
        call print_vtk(t2,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
      else
        call print_vtk(t2,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu) 
      endif
    enddo
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q,Vin)
    integer(kind=4), intent(in) :: id_RungeKutta
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(8), intent(in)         :: Vin(ny,2)
    integer i, j, t1, t2, itr, stat
    real(8), dimension(nx,ny,nz) :: mut_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), dimension(nx,ny,nz,5), device  :: QJ, QJs
    real(8), dimension(nx,ny,nz), device    :: mut
    real(8), device                         :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device                         :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device                         :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), device                         :: Rs(nx-accuracy,ny-accuracy,nz-accuracy,5)
    real(8), dimension(nx-1), device        :: dx, xix  ! xix = 1 / dx
    real(8), dimension(ny-1), device        :: dy, etay ! etay = 1 / dy
    real(8), dimension(nx,ny), device       :: Jacobian
    ! for plot
    real(8) ke0, entropy0

    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
    enddo;enddo

    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
    else
      call print_vtk(0,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
    endif

    ! copy on GPU
    QJ = Q
    mut = 0.d0
    Rs(:,:,:,:) = 0.d0
    xix = xix_cpu
    etay = etay_cpu
    dx = dx_cpu
    dy = dy_cpu
    Jacobian = Jacobian_cpu
    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJ,mut,E,F,G)
        call calc_step(0.5d0,1.d0,dx,dy,E,F,G,QJ,QJs,Rs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step(0.5d0,2.d0,dx,dy,E,F,G,QJ,QJs,Rs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step(1.d0,2.d0,dx,dy,E,F,G,QJ,QJs,Rs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step4(dx,dy,E,F,G,Rs,QJ)
        call set_bc(id_accuracy,Jacobian,QJ)
      enddo
      Q = QJ
      mut_cpu = mut
      if (id_turbulence == 0) then
        call print_vtk(t2,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
      else
        call print_vtk(t2,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
      endif
    enddo
  end subroutine RungeKutta_4th

  subroutine RungeKutta_10th(id_RungeKutta,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q,Vin)
    integer(kind=8), intent(in) :: id_RungeKutta  
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(8), intent(in)         :: Vin(ny,2)
    integer i, j, t1, t2, itr, stat
    real(8), dimension(nx,ny,nz) :: mut_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), dimension(nx,ny,nz,5), device  :: QJ, QJs, QJ4
    real(8), dimension(nx,ny,nz), device    :: mut
    real(8), device                         :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device                         :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device                         :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), device                         :: R4(nx-accuracy,ny-accuracy,nz-accuracy,5)
    real(8), dimension(nx-1), device        :: dx, xix  ! xix = 1 / dx
    real(8), dimension(ny-1), device        :: dy, etay ! etay = 1 / dy
    real(8), dimension(nx,ny), device       :: Jacobian
    ! for plot
    real(8) ke0, entropy0

    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
    enddo;enddo

    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
    else
      call print_vtk(0,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
    endif

    ! copy on GPU
    QJ = Q
    mut = 0.d0
    xix = xix_cpu
    etay = etay_cpu
    dx = 1.d0 / dx_cpu
    dy = 1.d0 / dy_cpu
    Jacobian = Jacobian_cpu
    do t2 = 1, np
      do t1 = 1, nt
        QJs = QJ
        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)
        QJ4 = QJs
  
        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step2(9.d0,6.d0,1.d0,15.d0,dx,dy,E,F,G,QJ,QJ4,QJs,R4)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(id_accuracy,Jacobian,QJs)

        call calc_EFG(id_hybrid,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step10(dx,dy,E,F,G,R4,QJ4,QJs,QJ)
        call set_bc(id_accuracy,Jacobian,QJ)
      enddo
      Q = QJ
      mut_cpu = mut
      if (id_turbulence == 0) then
        call print_vtk(t2,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
      else
        call print_vtk(t2,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
      endif
    enddo
  end subroutine RungeKutta_10th
end module calc_time_dev

