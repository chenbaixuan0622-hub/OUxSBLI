module set_init_common
  use mod_globals, only : gamma, R, Taw, rf
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  use set_compressible_bl
  implicit none
contains
  subroutine set_init_tbl(nx, ny, x, y, blt0, blt, u0, p0, T0, M0, Q)
    integer, intent(in)  :: nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(in)  :: blt0, blt, u0, p0, T0, M0
    real(8), intent(out) :: Q(nx,4,ny)
    integer i, j
    real(8) :: p_wall
    real(8), allocatable :: rho(:), u(:), v(:), T(:)
    allocate(rho(ny), u(ny), v(ny), T(ny))
    call calc_HD_Blasius(ny, y, blt0, u0, T0, p0, M0, rho, u, v, T)
    do j = 1, ny
      do i = 1, nx
        Q(i,1,j) = p0 / (R * T(j))
        Q(i,2,j) = Q(i,1,j) * u(j)
        Q(i,3,j) = Q(i,1,j) * v(j)
        Q(i,4,j) = p0 * over_gamma_1 + 0.5d0 * (Q(i,2,j)**2 + Q(i,3,j)**2) / Q(i,1,j)
    enddo;enddo
    deallocate(rho, u, v, T)
    ! bottom
    Q(:,1,1) = Q(:,1,2)
    Q(:,2,1) = 0.d0
    Q(:,3,1) = 0.d0
    p_wall = gamma_1 * (Q(2,4,2) - 0.5d0 * (Q(2,2,2)**2 + Q(2,3,2)**2) / Q(2,1,2))
    Q(:,4,1) = p_wall * over_gamma_1
  end subroutine set_init_tbl
end module set_init_common

