module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, rhol => rho0, rhor => rho1, pl => p0, pr => p1
  use set_bc_common
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
    real(8), intent(out) :: Q(4,nx,ny)
    integer i, j
    do i = 1, nx
      if (i < int(0.5*nx)) then
        Q(1,i,:) = rhol
        Q(2,i,:) = 0.d0
        Q(3,i,:) = 0.d0
        Q(4,i,:) = pl / (gamma  - 1.d0)
      else
        Q(1,i,:) = rhor
        Q(2,i,:) = 0.d0
        Q(3,i,:) = 0.d0
        Q(4,i,:) = pr / (gamma - 1.d0)
      endif
    enddo
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, Jacobian, Q)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: Q(4,nx,ny) ! Q / J
    integer i, j, k
    ! inlet and outlet
    !$cuf kernel do(1) <<<*,*>>>
    do j = 1, ny
      Q(1,1,j)    = rhol / Jacobian(j)
      Q(2,1,j)    = 0.d0
      Q(3,1,j)    = 0.d0
      Q(4,1,j)    = pl / (gamma - 1.d0) / Jacobian(j)
      Q(1,2,j)    = rhol / Jacobian(j)
      Q(2,2,j)    = 0.d0
      Q(3,2,j)    = 0.d0
      Q(4,2,j)    = pl / (gamma - 1.d0) / Jacobian(j)
      Q(1,3,j)    = rhol / Jacobian(j)
      Q(2,3,j)    = 0.d0
      Q(3,3,j)    = 0.d0
      Q(4,3,j)    = pl / (gamma - 1.d0) / Jacobian(j)
      Q(1,nx-2,j) = rhor / Jacobian(j)
      Q(2,nx-2,j) = 0.d0
      Q(3,nx-2,j) = 0.d0
      Q(4,nx-2,j) = pr / (gamma - 1.d0) / Jacobian(j)
      Q(1,nx-1,j) = rhor / Jacobian(j)
      Q(2,nx-1,j) = 0.d0
      Q(3,nx-1,j) = 0.d0
      Q(4,nx-1,j) = pr / (gamma - 1.d0) / Jacobian(j)
      Q(1,nx,j)   = rhor / Jacobian(j)
      Q(2,nx,j)   = 0.d0
      Q(3,nx,j)   = 0.d0
      Q(4,nx,j)   = pr / (gamma - 1.d0) / Jacobian(j)
    enddo
    !$cuf kernel do(2)<<<*,*>>>
    do i = 4, nx-3
      do k = 1, 4
        Q(k,i,1)    = Q(k,i,4)
        Q(k,i,2)    = Q(k,i,4)
        Q(k,i,3)    = Q(k,i,4)
        Q(k,i,ny-2) = Q(k,i,4)
        Q(k,i,ny-1) = Q(k,i,4)
        Q(k,i,ny)   = Q(k,i,4)
    enddo;enddo
  end subroutine set_bc

  subroutine calc_forcing(nx, ny, dx, dy, rho, u, v, p, fx, fy)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: rho(nx,ny), u(nx,ny), v(nx,ny), p(nx,ny)
    real(8), intent(out), device :: fx(nx-2,ny-2), fy(nx-2,ny-2)
  end subroutine calc_forcing
end module set

