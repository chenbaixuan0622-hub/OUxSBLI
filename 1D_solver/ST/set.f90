module set
  use cudafor
  use mod_globals, only : dx, rho0, p0, rho1, p1
  use mod_constant, only : over_gamma_1
  implicit none
contains
  subroutine set_grid(myrank, nx, x)
    integer, intent(in)  :: myrank, nx
    real(8), intent(out) :: x(nx)
    integer i
    do i = 1, nx
      x(i) = dble(i-1) * dx
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, x, Q)
    integer, intent(in)  :: myrank, nx
    real(8), intent(in)  :: x(nx)
    real(8), intent(out) :: Q(nx,3)
    integer i
    do i = 1, nx/2
      Q(i,1) = rho0
      Q(i,2) = 0.d0
      Q(i,3) = p0 * over_gamma_1
    enddo
    do i = nx/2+1, nx
      Q(i,1) = rho1
      Q(i,2) = 0.d0
      Q(i,3) = p1 * over_gamma_1
    enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, Q)
    integer, intent(in), value     :: myrank, nx
    real(8), intent(inout), device :: Q(nx,3)
    integer l
    !$cuf kernel do<<<*,*>>>
    do l = 1, 3
      Q(1,l)  = Q(2,l)
      Q(nx,l) = Q(nx-1,l)
    enddo
  end subroutine set_bc
end module set
