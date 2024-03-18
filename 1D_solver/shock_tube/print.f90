module print
  implicit none
contains
  subroutine print_vtk(step,nx,dx,gamma,Q)
    integer, intent(in) :: step, nx
    real(8), intent(in) :: dx, gamma, Q(nx,3)
    integer i
    real(8), dimension(nx) :: rho, u, p
    character(len=40) filename
    rho = Q(:,1)
    u = Q(:,2) / rho
    p = (gamma - 1.d0) * (Q(:,3) - 0.5d0 * rho * u ** 2)
    
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".d"
    open(10,file=filename)
    write(10,"('# Q.d')")
    write(10,"('# x     y0      y1      y2')")
    do i = 1, nx
      write(10,"(4(f9.4,1x))") (i-1)*dx, rho(i), u(i), p(i)
    enddo
    write(10,"('# end of file')")
    close(10)
  end subroutine print_vtk
end module print

