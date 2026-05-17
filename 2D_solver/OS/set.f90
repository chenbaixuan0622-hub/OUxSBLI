module set
  use cudafor
  use mod_globals, only : rho0, p0, u0, rho2, p2, ux, uy
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1
  implicit none
contains
  ! Uniform Cartesian grid.
  subroutine set_grid(myrank, nx, ny, Lx, Ly, xc, yc, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: xc(nx), yc(ny), dx(nx-1), dy(ny-1)
    real(8) :: dx0, dy0
    integer :: i, j
    dx0 = Lx / dble(nx - 1)
    dy0 = Ly / dble(ny - 1)
    dx  = dx0
    dy  = dy0
    xc(1) = 0.0d0
    do i = 1, nx - 1
      xc(i + 1) = xc(i) + dx0
    end do
    yc(1) = 0.0d0
    do j = 1, ny - 1
      yc(j + 1) = yc(j) + dy0
    end do
  end subroutine set_grid

  subroutine set_init(myrank, nx, ny, x, y, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(out) :: Q(nx,4,ny)
    integer :: i, j
    do j = 1, ny
      do i = 1, nx
        Q(i,1,j) = rho0
        Q(i,2,j) = rho0 * u0
        Q(i,3,j) = 0.d0
        Q(i,4,j) = p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2
    enddo;enddo
    do i = int(0.1d0 * nx), nx
      Q(i,1,ny) = rho2
      Q(i,2,ny) = rho2 * ux
      Q(i,3,ny) = rho2 * uy
      Q(i,4,ny) = p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo
  end subroutine set_init

  subroutine set_bc(myrank, nx, ny, Jacobian, QJ)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny)
    real(8) Jacobian_tmp
    integer :: i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      Jacobian_tmp = 1.d0 / Jacobian(1,ny-1)
      if (i < int(0.1d0 * nx)) then
        QJ(i,1,ny) = rho0 * Jacobian_tmp
        QJ(i,2,ny) = rho0 * u0 * Jacobian_tmp
        QJ(i,3,ny) = 0.d0
        QJ(i,4,ny) = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      else
        QJ(i,1,ny) = rho2 * Jacobian_tmp
        QJ(i,2,ny) = rho2 * ux * Jacobian_tmp
        QJ(i,3,ny) = rho2 * uy * Jacobian_tmp
        QJ(i,4,ny) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
      endif
      ! Slip
      QJ(i,1,1) =  QJ(i,1,2)
      QJ(i,2,1) =  QJ(i,2,2)
      QJ(i,3,1) = -QJ(i,3,2)
      QJ(i,4,1) =  QJ(i,4,2)
    enddo
    !$cuf kernel do(1)<<<*,*>>>
    do j = 1, ny
      Jacobian_tmp = 1.d0 / Jacobian(1,j)
      QJ(1,1,j)  = rho0 * Jacobian_tmp
      QJ(1,2,j)  = rho0 * u0 * Jacobian_tmp
      QJ(1,3,j)  = 0.d0
      QJ(1,4,j)  = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      QJ(nx,1,j) = QJ(nx-1,1,j)
      QJ(nx,2,j) = QJ(nx-1,2,j)
      QJ(nx,3,j) = QJ(nx-1,3,j)
      QJ(nx,4,j) = QJ(nx-1,4,j)
    enddo
  end subroutine set_bc
end module set
