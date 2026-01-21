module set
  use mod_globals, only : nx, Lx, dx, gamma, R, rhod, rhou, ud, uu, pd, pu
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
        Q(1,i) = rhod
        Q(2,i) = rhod * ud
        Q(3,i) = pd * over_gamma_1 + 0.5d0 * rhod * ud**2
      else
        Q(1,i) = rhou
        Q(2,i) = rhou * uu
        Q(3,i) = pu * over_gamma_1 + 0.5d0 * rhou * uu**2
      endif
    enddo
  end subroutine set_init


  subroutine set_bc(nx,Q)
    integer, intent(in), value     :: nx
    real(8), intent(inout), device :: Q(3,nx)
    Q(1,1)  = rhod
    Q(2,1)  = rhod * ud
    Q(3,1)  = pd * over_gamma_1 + 0.5d0 * rhod * ud**2
    Q(1,nx) = rhou
    Q(2,nx) = rhou * uu
    Q(3,nx) = pu * over_gamma_1 + 0.5d0 * rhou * uu**2
  end subroutine set_bc
end module set

