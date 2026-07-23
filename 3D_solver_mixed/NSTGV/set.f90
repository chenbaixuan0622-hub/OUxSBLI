module set
  use mod_globals, only : nx, ny, nz, gamma, R, RHO0, V0, p0, L0, T, dtn
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    use mod_constant, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    call set_grid_cyclic(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_constant, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,5,ny,nz)
    integer i, j, k, offset
    real(8) p, RHO
    if (kind(id_accuracy) == 2) then
      offset = 1
    elseif (kind(id_accuracy) == 4) then
      offset = 2
    elseif (kind(id_accuracy) == 8) then
      offset = 3
    endif
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          p = p0+RHO0*(V0**2)*(cos(2.d0*x(i)/L0)+cos(2.d0*y(j)/L0))*(cos(2.d0*z(k)/L0)+2.d0)/16.d0
          RHO = p / (R * T)
          ! rho
          Q(i,1,j,k) = RHO
          ! rho u
          Q(i,2,j,k) = RHO * V0 * sin(x(i)/L0) * cos(y(j)/L0) * cos(z(k)/L0)
          ! rho v
          Q(i,3,j,k) = - RHO * V0 * cos(x(i)/L0) * sin(y(j)/L0) * cos(z(k)/L0)
          ! rho w
          Q(i,4,j,k) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q(i,5,j,k) = p / (gamma - 1.d0) + 0.5d0 * (Q(i,2,j,k) ** 2 + Q(i,3,j,k) ** 2) / Q(i,1,j,k)
    enddo;enddo;enddo
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
    use mod_constant, only : id_accuracy
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
