module set
  use mod_globals, only : nx, Lx, dx, gamma, R, rho0, rho1, p0, p1
  use mod_constant, only : over_gamma_1
  implicit none
contains
  subroutine set_grid(nx, x)
    integer, intent(in)  :: nx
    real(8), intent(out) :: x(nx)
    integer i
    x(1) = 0.d0
    do i = 1, nx-1
      x(i+1) = x(i) + dx
    enddo
  end subroutine set_grid
  
  subroutine set_init(nx,x,Q)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: x(nx)
    real(8), intent(out) :: Q(3,nx)
    integer i
    do i = 1, nx
      if (i < int(0.5 * nx)) then
        Q(1,i) = rho0
        Q(2,i) = 0.d0
        Q(3,i) = p0 * over_gamma_1
      else
        Q(1,i) = rho1
        Q(2,i) = 0.d0
        Q(3,i) = p1 * over_gamma_1
      endif
    enddo
  end subroutine set_init

  subroutine set_bc(nx,Q)
    integer, intent(in), value     :: nx
    real(8), intent(inout), device :: Q(3,nx)
    Q(1,1)  = rho0
    Q(2,1)  = 0.d0
    Q(3,1)  = p0 * over_gamma_1
    Q(1,nx) = rho1
    Q(2,nx) = 0.d0
    Q(3,nx) = p1 * over_gamma_1
  end subroutine set_bc
end module set

