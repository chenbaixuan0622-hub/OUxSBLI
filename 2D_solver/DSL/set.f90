module set
  use mod_globals, only : nx, ny, nz, gamma, R, dtn
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(1)
    real(8) x(nx+1), y(ny+1)
    integer i, j
    dx(:) = Lx / dble(nx-6)
    dy(:) = Ly / dble(ny-6)

    ! x direction
    do i = 4, nx-2
      x(i) = dx(1) * dble(i-4)
    enddo
    x(1)    = x(4)    - 3.d0 * dx(1)
    x(2)    = x(4)    - 2.d0 * dx(1)
    x(3)    = x(4)    - dx(1)
    x(nx-1) = x(nx-2) + dx(1)
    x(nx)   = x(nx-2) + 2.d0 * dx(1)
    x(nx+1) = x(nx-2) + 3.d0 * dx(1)

    ! y direction
    do j = 4, ny-2
      y(j) = dy(1) * dble(j-4)
    enddo
    y(1)    = y(4)    - 3.d0 * dy(1)
    y(2)    = y(4)    - 2.d0 * dy(1)
    y(3)    = y(4)    - dy(1)
    y(ny-1) = y(ny-2) + dy(1)
    y(ny)   = y(ny-2) + 2.d0 * dy(1)
    y(ny+1) = y(ny-2) + 3.d0 * dy(1)

    ! cell centered
    do i = 1, nx
      xc(i) = 0.5d0 * (x(i) + x(i+1))
    enddo
    do j = 1, ny
      yc(j) = 0.5d0 * (y(j) + y(j+1))
    enddo
  end subroutine set_grid
  
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : pi, M0, rho0, u0, d1, d2
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,4)
    integer i, j
    real(8) :: Cp = R * gamma / (gamma - 1.d0), p = rho0 * u0**2 / (gamma * M0**2)
    do j = 1, ny
      do i = 1, nx
        if (y(j) <= pi) then
          Q(i,j,1) = rho0
          Q(i,j,2) = rho0 * u0 * tanh((y(j) - 0.5d0 * pi) / d1)
          Q(i,j,3) = rho0 * u0 * d2 * sin(x(i))
          Q(i,j,4) = p / (gamma - 1.d0) + 0.5d0 * (Q(i,j,2)**2 + Q(i,j,3)**2) / Q(i,j,1)
        else
          Q(i,j,1) = rho0
          Q(i,j,2) = rho0 * u0 * tanh((1.5d0 * pi - y(j)) / d1)
          Q(i,j,3) = rho0 * u0 * d2 * sin(x(i))
          Q(i,j,4) = p / (gamma - 1.d0) + 0.5d0 * (Q(i,j,2)**2 + Q(i,j,3)**2) / Q(i,j,1)
        endif
    enddo;enddo
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, Jacobian, Q)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, 4
      do j = 4, ny-3
        Q(1,j,k) = Q(nx-5,j,k)
        Q(2,j,k) = Q(nx-4,j,k)
        Q(3,j,k) = Q(nx-3,j,k)
        Q(nx-2,j,k) = Q(4,j,k)
        Q(nx-1,j,k) = Q(5,j,k)
        Q(nx,j,k)   = Q(6,j,k)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, 4
      do i = 4, nx-3
        Q(i,1,k) = Q(i,ny-5,k)
        Q(i,2,k) = Q(i,ny-4,k)
        Q(i,3,k) = Q(i,ny-3,k)
        Q(i,ny-2,k) = Q(i,4,k)
        Q(i,ny-1,k) = Q(i,5,k)
        Q(i,ny,k)   = Q(i,6,k)
    enddo;enddo

    !$cuf kernel do(1)<<<*,*>>>
    do k = 1, 4
      Q(1,1,k) = Q(nx-5,ny-5,k)
      Q(1,2,k) = Q(nx-5,ny-4,k)
      Q(1,3,k) = Q(nx-5,ny-3,k)
      Q(2,1,k) = Q(nx-4,ny-5,k)
      Q(2,2,k) = Q(nx-4,ny-4,k)
      Q(2,3,k) = Q(nx-4,ny-3,k)
      Q(3,1,k) = Q(nx-3,ny-5,k)
      Q(3,2,k) = Q(nx-3,ny-4,k)
      Q(3,3,k) = Q(nx-3,ny-3,k)
      Q(nx-2,1,k) = Q(4,ny-5,k)
      Q(nx-2,2,k) = Q(4,ny-4,k)
      Q(nx-2,3,k) = Q(4,ny-3,k)
      Q(nx-1,1,k) = Q(5,ny-5,k)
      Q(nx-1,2,k) = Q(5,ny-4,k)
      Q(nx-1,3,k) = Q(5,ny-3,k)
      Q(nx,1,k)   = Q(6,ny-5,k)
      Q(nx,2,k)   = Q(6,ny-4,k)
      Q(nx,3,k)   = Q(6,ny-3,k)
      Q(1,ny-2,k) = Q(nx-5,4,k)
      Q(1,ny-1,k) = Q(nx-5,5,k)
      Q(1,ny,k)   = Q(nx-5,6,k)
      Q(2,ny-2,k) = Q(nx-4,4,k)
      Q(2,ny-1,k) = Q(nx-4,5,k)
      Q(2,ny,k)   = Q(nx-4,6,k)
      Q(3,ny-2,k) = Q(nx-3,4,k)
      Q(3,ny-1,k) = Q(nx-3,5,k)
      Q(3,ny,k)   = Q(nx-3,6,k)
      Q(nx-2,ny-2,k) = Q(4,4,k)
      Q(nx-2,ny-1,k) = Q(4,5,k)
      Q(nx-2,ny,k)   = Q(4,6,k)
      Q(nx-1,ny-2,k) = Q(5,4,k)
      Q(nx-1,ny-1,k) = Q(5,5,k)
      Q(nx-1,ny,k)   = Q(5,6,k)
      Q(nx,ny-2,k)   = Q(6,4,k)
      Q(nx,ny-1,k)   = Q(6,5,k)
      Q(nx,ny,k)     = Q(6,6,k)
    enddo
  end subroutine set_bc

  subroutine calc_forcing(nx, ny, dx, dy, rho, u, v, p, fx, fy)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: rho(nx,ny), u(nx,ny), v(nx,ny), p(nx,ny)
    real(8), intent(out), device :: fx(nx-2,ny-2,4), fy(nx-2,ny-2,4)
  end subroutine calc_forcing
end module set

