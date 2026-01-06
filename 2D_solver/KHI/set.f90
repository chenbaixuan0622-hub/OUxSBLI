module set
  use mod_globals, only : id_accuracy, nx, ny, gamma, R
  use mod_constant, only : Cp
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
    use mod_globals, only : pi, Lx, amp, rho1, u1, rho2, u2, pd
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(4,nx,ny)
    real(8) v
    integer i, j
    do j = 1, ny
      do i = 1, nx
        v = amp * sin(2.d0 * pi * x(i) / Lx)
        if (y(j) > 0.75d0 * Lx .or. y(j) < 0.25d0 * Lx) then
          Q(1,i,j) = rho1
          Q(2,i,j) = rho1 * u1
          Q(3,i,j) = rho1 * v
          Q(4,i,j) = pd / (gamma - 1.d0) + 0.5d0 * rho1 * (u1**2 + v**2)
        else
          Q(1,i,j) = rho2
          Q(2,i,j) = rho2 * u2
          Q(3,i,j) = rho2 * v
          Q(4,i,j) = pd / (gamma - 1.d0) + 0.5d0 * rho2 * (u2**2 + v**2)
        endif
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

