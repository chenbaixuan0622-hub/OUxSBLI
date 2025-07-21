module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, rhol => rho0, rhor => rho1, pl => p0, pr => p1
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j
    real(8) dx1, dy1
    dx1 = Lx / dble(nx-1)
    dy1 = Ly / dble(ny-1)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    y(1) = 0.d0
    do j = 1, ny-1
      dy(j) = dy1
      y(j+1) = y(j) + dy(j)
    enddo
  end subroutine set_grid
  
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,4)
    integer i, j
    do i = 1, nx
      if (i < int(0.5*nx)) then
        Q(i,:,1) = rhol
        Q(i,:,2) = 0.d0
        Q(i,:,3) = 0.d0
        Q(i,:,4) = pl / (gamma  - 1.d0)
      else
        Q(i,:,1) = rhor
        Q(i,:,2) = 0.d0
        Q(i,:,3) = 0.d0
        Q(i,:,4) = pr / (gamma - 1.d0)
      endif
    enddo
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, Jacobian, Q)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: Q(nx,ny,4) ! Q / J
    integer i, j, k
  
    ! inlet and outlet
    !$cuf kernel do(1) <<<*,*>>>
    do j = 1, ny
      Q(1,j,1)    = rhol / Jacobian(j)
      Q(1,j,2)    = 0.d0
      Q(1,j,3)    = 0.d0
      Q(1,j,4)    = pl / (gamma - 1.d0) / Jacobian(j)
      Q(2,j,1)    = rhol / Jacobian(j)
      Q(2,j,2)    = 0.d0
      Q(2,j,3)    = 0.d0
      Q(2,j,4)    = pl / (gamma - 1.d0) / Jacobian(j)
      Q(3,j,1)    = rhol / Jacobian(j)
      Q(3,j,2)    = 0.d0
      Q(3,j,3)    = 0.d0
      Q(3,j,4)    = pl / (gamma - 1.d0) / Jacobian(j)
      Q(nx-2,j,1) = rhor / Jacobian(j)
      Q(nx-2,j,2) = 0.d0
      Q(nx-2,j,3) = 0.d0
      Q(nx-2,j,4) = pr / (gamma - 1.d0) / Jacobian(j)
      Q(nx-1,j,1) = rhor / Jacobian(j)
      Q(nx-1,j,2) = 0.d0
      Q(nx-1,j,3) = 0.d0
      Q(nx-1,j,4) = pr / (gamma - 1.d0) / Jacobian(j)
      Q(nx,j,1)   = rhor / Jacobian(j)
      Q(nx,j,2)   = 0.d0
      Q(nx,j,3)   = 0.d0
      Q(nx,j,4)   = pr / (gamma - 1.d0) / Jacobian(j)
    enddo
    
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, 4
      do i = 4, nx-3
        Q(i,1,k)    = Q(i,4,k)
        Q(i,2,k)    = Q(i,4,k)
        Q(i,3,k)    = Q(i,4,k)
        Q(i,ny-2,k) = Q(i,4,k)
        Q(i,ny-1,k) = Q(i,4,k)
        Q(i,ny,k)   = Q(i,4,k)
    enddo;enddo
  end subroutine set_bc

  subroutine calc_forcing(nx, ny, dx, dy, rho, u, v, p, fx, fy)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: rho(nx,ny), u(nx,ny), v(nx,ny), p(nx,ny)
    real(8), intent(out), device :: fx(nx-2,ny-2,4), fy(nx-2,ny-2,4)
  end subroutine calc_forcing
end module set

