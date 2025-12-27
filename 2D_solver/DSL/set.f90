module set
  use mod_globals, only : id_accuracy, nx, ny, gamma, R, dtn
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
    use mod_globals, only : pi, M0, rho0, u0, d1, d2
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(4,nx,ny)
    integer i, j
    real(8) :: p = rho0 * u0**2 / (gamma * M0**2)
    do j = 1, ny
      do i = 1, nx
        if (y(j) <= pi) then
          Q(1,i,j) = rho0
          Q(2,i,j) = rho0 * u0 * tanh((y(j) - 0.5d0 * pi) / d1)
          Q(3,i,j) = rho0 * u0 * d2 * sin(x(i))
          Q(4,i,j) = p / (gamma - 1.d0) + 0.5d0 * (Q(2,i,j)**2 + Q(3,i,j)**2) / Q(1,i,j)
        else
          Q(1,i,j) = rho0
          Q(2,i,j) = rho0 * u0 * tanh((1.5d0 * pi - y(j)) / d1)
          Q(3,i,j) = rho0 * u0 * d2 * sin(x(i))
          Q(4,i,j) = p / (gamma - 1.d0) + 0.5d0 * (Q(2,i,j)**2 + Q(3,i,j)**2) / Q(1,i,j)
        endif
    enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, Q)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: Q(4,nx,ny)
    call set_bc_cyclic(id_accuracy, nx, ny, Q)
  end subroutine set_bc
end module set

