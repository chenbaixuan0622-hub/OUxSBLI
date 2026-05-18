module set
  use mod_globals, only : id_accuracy, nx, ny, gamma, R
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, xc, yc, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: xc(nx), yc(ny), dx(nx-1), dy(ny-1)
    call set_grid_cyclic(id_accuracy, nx, ny, Lx, Ly, xc, yc, dx, dy)
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, x, y, Q)
    use mod_globals, only : pi, M0, rho0, u0, d1, d2
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(out) :: Q(nx,4,ny)
    integer i, j
    real(8) :: p = rho0 * u0**2 / (gamma * M0**2)
    do j = 1, ny
      do i = 1, nx
        if (y(j) <= pi) then
          Q(i,1,j) = rho0
          Q(i,2,j) = rho0 * u0 * tanh((y(j) - 0.5d0 * pi) / d1)
          Q(i,3,j) = rho0 * u0 * d2 * sin(x(i))
          Q(i,4,j) = p / (gamma - 1.d0) + 0.5d0 * (Q(i,2,j)**2 + Q(i,3,j)**2) / Q(i,1,j)
        else
          Q(i,1,j) = rho0
          Q(i,2,j) = rho0 * u0 * tanh((1.5d0 * pi - y(j)) / d1)
          Q(i,3,j) = rho0 * u0 * d2 * sin(x(i))
          Q(i,4,j) = p / (gamma - 1.d0) + 0.5d0 * (Q(i,2,j)**2 + Q(i,3,j)**2) / Q(i,1,j)
        endif
    enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, Q)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Q(nx,4,ny)
    call set_bc_cyclic(id_accuracy, nx, ny, Q)
  end subroutine set_bc
end module set

