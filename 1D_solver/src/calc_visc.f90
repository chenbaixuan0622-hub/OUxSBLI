module calc_visc
  use mod_globals, only : accuracy, offset
  use calc_Sutherland
  implicit none
contains
  subroutine calc_Ev(nx, dX, u, T, E)
    integer, intent(in), value :: nx
    real(8), intent(in), value :: dX
    real(8), intent(in), dimension(nx) :: u, T
    real(8), intent(inout), dimension(nx-accuracy+1,3) :: E
    integer i
    real(8) mu, kappa
    real(8) txx
    do i = 1 + offset - 1, nx - offset
      call calc_mu(T(i),T(i+1),mu)
      txx = 4.d0 * mu * dX * (-u(i) + u(i+1)) / 3.d0
      call calc_kappa(T(i),T(i+1),kappa)
      E(i-offset+1,2) = E(i-offset+1,2) - txx
      E(i-offset+1,3) = E(i-offset+1,3) - txx * 0.5d0 *(u(i) + u(i+1)) - kappa * (-T(i) + T(i+1))
    enddo
  end subroutine calc_Ev
end module calc_visc

