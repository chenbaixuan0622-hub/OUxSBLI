module set
  use mod_globals, only : nx, ny, nz, gamma, R, RHO0, V0, p0, L0, T, dtn
  use set_bc_common
  use set_coordinate
  use set_init_dhit
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    call set_grid_cyclic(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,5,ny,nz)
    call init_spectral_velocity(nx, ny, nz, Q)
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
    use mod_globals, only : id_accuracy
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny,nz)
    real(8), intent(inout), device :: Q(nx,5,ny,nz)
    real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_bc


  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx,ny,nz,mut,qc2)
  end subroutine set_bc_mut
end module set
