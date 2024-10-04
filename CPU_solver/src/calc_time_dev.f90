module calc_time_dev
  use mod_globals, only : nt, np
  use calc_physical_quantities
  use calc_steps
  use calc_flux
  use set
  use print
  implicit none
contains
  subroutine calc_EFG(nx,ny,nz,dx,dy,dz,Jacobian,QJ,E,F,G)
    integer, intent(in)                         :: nx, ny, nz
    real(8), intent(in), dimension(nx-1)        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1)        :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1)        :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz)    :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5)  :: QJ ! Q / Jacobian
    real(8), intent(out)                        :: E(nx-1,ny-2,nz-2,5)
    real(8), intent(out)                        :: F(nx-2,ny-1,nz-2,5)
    real(8), intent(out)                        :: G(nx-2,ny-2,nz-1,5)
    real(8), allocatable :: rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:), p(:,:,:), T(:,:,:)
    integer stat
    allocate(rho(nx,ny,nz), u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), p(nx,ny,nz), T(nx,ny,nz))

    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)
    
    call calc_E(nx,ny,nz,rho,u,v,w,p,E)
    call calc_F(nx,ny,nz,rho,u,v,w,p,F)
    call calc_G(nx,ny,nz,rho,u,v,w,p,G)
    
    deallocate(rho, u, v, w, p, T)
  end subroutine calc_EFG

  subroutine RungeKutta(nx,ny,nz,x,dx,xix,y,dy,etay,z,dz,zetaz,Jacobian,Q)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(in)    :: x(nx), dx(nx-1), xix(nx-1)
    real(8), intent(in)    :: y(ny), dy(ny-1), etay(ny-1)
    real(8), intent(in)    :: z(nz), dz(nz-1), zetaz(nz-1), Jacobian(nx,ny,nz)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr 
    real(8), allocatable :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    ! for plot
    real(4) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0

    allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),Rs(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))

    ! set Q / Jacobian
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          QJ(i,j,k,:) = Q(i,j,k,:) / Jacobian(i,j,k)
    enddo;enddo;enddo

    ! print initial condition
    call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian),real(QJ),mass0,ke0,entropy0)

    Rs(:,:,:,:) = 0.d0

    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJ,E,F,G)
        call calc_step(nx,ny,nz,0.5d0,1.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q2
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,E,F,G)
        call calc_step(nx,ny,nz,0.5d0,2.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q3
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,E,F,G)
        call calc_step(nx,ny,nz,1.0d0,2.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q4
        call set_bc(nx,ny,nz,Jacobian,QJs)

        call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,E,F,G)
        call calc_step4(nx,ny,nz,dx,dy,dz,E,F,G,Rs,QJ)
        call set_bc(nx,ny,nz,Jacobian,QJ)
      enddo
      call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian),real(QJ),mass0,ke0,entropy0)
    enddo

    deallocate(QJ,QJs,Rs,E,F,G)
  end subroutine RungeKutta
end module calc_time_dev

