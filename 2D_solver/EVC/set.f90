module set
  use mod_globals, only : id_accuracy, nx, ny, Lx, Ly, gamma, R, dtn, theta
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(1)
    call set_grid_cyclic(id_accuracy, nx, ny, Lx, Ly, xc, yc, dx, dy)
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : M0, rho0, p0, T0, u0, Rc, beta
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(4,nx,ny)
    integer i, j
    real(8) xr, yr, xc, yc, ex, T, rho, u, v, du, dv, p
    real(8) :: Cp = R * gamma / (gamma - 1.d0)
    xc = x(int(nx/2))
    yc = y(int(ny/2))
    do j = 1, ny
      do i = 1, nx
        ex  = exp(-0.5d0 * ((x(i) - xc)**2 + (y(j) - yc)**2) / (Rc**2))
        T   = T0 - 0.5d0 * (u0 * beta)**2 / Cp * ex**2
        rho = rho0 * (T / T0)**(1.d0  / (gamma - 1.d0))
        du  = -u0 * beta * (y(j) - yc) / Rc * ex
        dv  =  u0 * beta * (x(i) - xc) / Rc * ex
        u   = u0 * cos(theta) + du * cos(theta) - dv * sin(theta)
        v   = u0 * sin(theta) + du * sin(theta) + dv * cos(theta)
        p   = rho * R * T
        Q(1,i,j) = rho
        Q(2,i,j) = rho * u
        Q(3,i,j) = rho * v
        ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
        Q(4,i,j) = p / (gamma - 1.d0) + 0.5d0 * rho * (u**2 + v**2)
    enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, x, y, Jacobian, Q)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: x(nx), y(ny)
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Q(4,nx,ny)
    call set_bc_cyclic(id_accuracy, nx, ny, Q)
  end subroutine set_bc


  attributes(global) subroutine calc_force(nx, ny, x, y, dx, dy, Q, Fout)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    real(8), intent(in), device  :: Q(4,nx,ny)
    real(8), intent(out), device :: Fout(3,nx-2,ny-2)
  end subroutine calc_force
end module set

