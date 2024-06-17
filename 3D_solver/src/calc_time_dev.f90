module calc_time_dev
  use cudafor
  use nvtx
  use mod_globals, only : accuracy, id_hybrid, id_muscl, id_visc, id_turbulence, nt, np, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_steps
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
  subroutine calc_EFG_basic(id_hybrid,nx,ny,nz,dx,dy,Jacobian,QJ,mut,E,F,G)
    use calc_flux
    integer(kind=2), intent(in)                         :: id_hybrid
    integer, intent(in), value                          :: nx, ny, nz
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
    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)

    call calc_E<<<blocksE,threadsE>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,E)
    call calc_F<<<blocksF,threadsF>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,F)
    call calc_G<<<blocksG,threadsG>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,G)
    !print *, trim(cudaGetErrorString(cudaGetLastError()))
    stat = cudaDeviceSynchronize()

    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads>>>(nx,ny,nz,rho,u,v,w,mut)
      stat = cudaDeviceSynchronize()
      call set_bc_mut(nx,ny,nz,mut)
    endif

    if (id_visc == 1 .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,E)
      call calc_Fv<<<blocksF,threadsF>>>(nx,ny,nz,dy,dx,rho,u,v,w,T,p,mut,F)
      call calc_Gv<<<blocksG,threadsG>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,G)
    endif
    !print *, trim(cudaGetErrorString(cudaGetLastError()))
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_basic
  
  subroutine calc_EFG_hybrid(id_hybrid,nx,ny,nz,dx,dy,Jacobian,QJ,mut,E,F,G)
    use calc_flux_hybrid
    integer(kind=4), intent(in)                         :: id_hybrid
    integer, intent(in), value                          :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device        :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device       :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(inout), dimension(nx,ny,nz), device :: mut
    real(8), intent(out), device                        :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                        :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                        :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device                                :: rho, u, v, w, p, T, fd
    !real(8), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_upwind
    !real(8), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_upwind
    !real(8), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_upwind
    integer stat
    !integer(kind=2) :: id_muscl1
    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)

    !call calc_E<<<blocksE,threadsE>>>(id_muscl1,nx,ny,nz,rho,u,v,w,p,E)
    !call calc_F<<<blocksF,threadsF>>>(id_muscl1,nx,ny,nz,rho,u,v,w,p,F)
    !call calc_G<<<blocksG,threadsG>>>(id_muscl1,nx,ny,nz,rho,u,v,w,p,G)
    !call calc_E<<<blocksE,threadsE>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,E_upwind)
    !call calc_F<<<blocksF,threadsF>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,F_upwind)
    !call calc_G<<<blocksG,threadsG>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,G_upwind)
    !stat = cudaDeviceSynchronize()

    call calc_Ducros<<<blocks,threads>>>(nx,ny,nz,dx,dy,u,v,w,rho,p,fd)
    stat = cudaDeviceSynchronize()
    call calc_E_hybrid<<<blocksE,threadsE>>>(nx,ny,nz,rho,u,v,w,p,fd,E)
    call calc_F_hybrid<<<blocksF,threadsF>>>(nx,ny,nz,rho,u,v,w,p,fd,F)
    call calc_G_hybrid<<<blocksG,threadsG>>>(nx,ny,nz,rho,u,v,w,p,fd,G)
    !call calc_E_hybrid<<<blocksE,threadsE>>>(nx,ny,nz,u,v,w,fd,E_upwind,E)
    !call calc_F_hybrid<<<blocksF,threadsF>>>(nx,ny,nz,u,v,w,fd,F_upwind,F)
    !call calc_G_hybrid<<<blocksG,threadsG>>>(nx,ny,nz,u,v,w,fd,G_upwind,G)
    stat = cudaDeviceSynchronize()
    !print *, trim(cudaGetErrorString(cudaGetLastError()))

    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads>>>(nx,ny,nz,rho,u,v,w,mut)
      stat = cudaDeviceSynchronize()
      call set_bc_mut(nx,ny,nz,mut)
    endif

    if (id_visc == 1 .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,E)
      call calc_Fv<<<blocksF,threadsF>>>(nx,ny,nz,dy,dx,rho,u,v,w,T,p,mut,F)
      call calc_Gv<<<blocksG,threadsG>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_hybrid

  subroutine RungeKutta_3rd(id_RungeKutta,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q,Vin)
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(8), intent(in)         :: Vin(ny,2)
    integer i, j, k, t1, t2, itr, stat
    real(8), dimension(nx,ny,nz)    :: mut_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), dimension(nx,ny,nz,5), device  :: QJ, QJ2, QJ3
    real(8), dimension(nx,ny,nz), device    :: mut
    real(8), device                         :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device                         :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device                         :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx-1), device        :: dx, xix  ! xix = 1 / dx
    real(8), dimension(ny-1), device        :: dy, etay ! etay = 1 / dy
    real(8), dimension(nx,ny), device       :: Jacobian
    real(8), dimension(nx,ny,nz,3), device  :: Vmean
    real(8), dimension(nx,ny,nz), device    :: rhomean, Tmean
    ! for plot
    real(8) ke0, entropy0

    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
    enddo;enddo
  
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Vmean_cpu(i,j,k,:) = Q(i,j,k,2:4) / Q(i,j,k,1)
    enddo;enddo;enddo
    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
    else
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
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
        call nvtxStartRange("calc 1step",1)
        call nvtxStartRange("calc flux",2)
        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ,mut,E,F,G)
        call nvtxEndRange
        call nvtxStartRange("calc time dev",3)
        call calc_step(nx,ny,nz,1.d0,0.d0,dx,dy,E,F,G,QJ,QJ2)
        call nvtxEndRange
        call set_bc(nx,ny,nz,Jacobian,QJ2)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ2,mut,E,F,G)
        call calc_step2(nx,ny,nz,0.75d0,0.25d0,0.25d0,1.d0,dx,dy,E,F,G,QJ,QJ2,QJ3)
        call set_bc(nx,ny,nz,Jacobian,QJ3)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ3,mut,E,F,G)
        call calc_step3(nx,ny,nz,dx,dy,E,F,G,QJ3,QJ)
        call set_bc(nx,ny,nz,Jacobian,QJ)
        call calc_mean(t1+(t2-1)*nt,nx,ny,nz,Jacobian,QJ,Vmean,rhomean,Tmean)
        call nvtxEndRange
      enddo
      Q = QJ
      mut_cpu = mut
      Vmean_cpu = Vmean
      if (id_turbulence == 0) then
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
      else
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu) 
      endif
    enddo
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q,Vin)
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(8), intent(in)         :: Vin(ny,2)
    integer i, j, k, t1, t2, itr, stat
    real(8), dimension(nx,ny,nz)    :: mut_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
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
    real(8), dimension(nx,ny,nz,3), device  :: Vmean
    real(8), dimension(nx,ny,nz), device    :: rhomean, Tmean
    ! for plot
    real(8) ke0, entropy0

    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
    enddo;enddo

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Vmean_cpu(i,j,k,:) = Q(i,j,k,2:4) / Q(i,j,k,1)
    enddo;enddo;enddo
    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
    else
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
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
        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ,mut,E,F,G)
        call calc_step(nx,ny,nz,0.5d0,1.d0,dx,dy,E,F,G,QJ,QJs,Rs) ! QJs = Q2
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step(nx,ny,nz,0.5d0,2.d0,dx,dy,E,F,G,QJ,QJs,Rs) ! QJs = Q3
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step(nx,ny,nz,1.0d0,2.d0,dx,dy,E,F,G,QJ,QJs,Rs) ! QJs = Q4
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step4(nx,ny,nz,dx,dy,E,F,G,Rs,QJ)
        call set_bc(nx,ny,nz,Jacobian,QJ)
        call calc_mean(t1+(t2-1)*nt,nx,ny,nz,Jacobian,QJ,Vmean,rhomean,Tmean)
      enddo
      Q = QJ
      mut_cpu = mut
      Vmean_cpu = Vmean
      if (id_turbulence == 0) then
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
      else
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
      endif
    enddo
  end subroutine RungeKutta_4th

  subroutine RungeKutta_10th(id_RungeKutta,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q,Vin)
    integer(kind=8), intent(in) :: id_RungeKutta  
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(8), intent(in)         :: Vin(ny,2)
    integer i, j, k, t1, t2, itr, stat
    real(8), dimension(nx,ny,nz)    :: mut_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
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
    real(8), dimension(nx,ny,nz,3), device  :: Vmean
    real(8), dimension(nx,ny,nz), device    :: rhomean, Tmean
    ! for plot
    real(8) ke0, entropy0

    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
    enddo;enddo

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Vmean_cpu(i,j,k,:) = Q(i,j,k,2:4) / Q(i,j,k,1)
    enddo;enddo;enddo
    ! print initial condition
    if (id_turbulence == 0) then
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
    else
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
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
        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)
        QJ4 = QJs
  
        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step2(nx,ny,nz,9.d0,6.d0,1.d0,15.d0,dx,dy,E,F,G,QJ,QJ4,QJs,R4)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
        call calc_step10(nx,ny,nz,dx,dy,E,F,G,R4,QJ4,QJs,QJ)
        call set_bc(nx,ny,nz,Jacobian,QJ)
        call calc_mean(t1+(t2-1)*nt,nx,ny,nz,Jacobian,QJ,Vmean,rhomean,Tmean)
      enddo
      Q = QJ
      mut_cpu = mut
      Vmean_cpu = Vmean
      if (id_turbulence == 0) then
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0)
      else
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,ke0,entropy0,mut_cpu)
      endif
    enddo
  end subroutine RungeKutta_10th
end module calc_time_dev

