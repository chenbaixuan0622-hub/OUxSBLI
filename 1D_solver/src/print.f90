module print
  use mod_globals, only : gamma, Lx
  implicit none
contains
  subroutine print_1d(step, nx, x, Q)
    integer, intent(in) :: step, nx
    real(4), intent(in) :: x(nx)
    real(4), intent(in) :: Q(nx,3)
    integer i
    real(4), dimension(nx) :: rho, u, p
    character(len=40) filename
    do i = 1, nx
      rho(i) = Q(i,1)
      u(i)   = Q(i,2) / rho(i)
      p(i)   = (real(gamma) - 1.e0) * (Q(i,3) - 0.5e0 * rho(i) * u(i)**2)
    enddo

    ! rho, u, p, T, Mach, mu
    write(filename, "(a, i5.5, a)") "data/Q", int(step), ".d"
    open(10,file=filename)
    do i = 1, nx
      write(10,"(7e12.4)") x(i)/Lx, rho(i), u(i), p(i)
    enddo
    close(10)
  end subroutine print_1d
end module print

