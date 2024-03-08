module calc_KEEP
  use calc_term
  implicit none
  interface calc_E
    module procedure calc_E2, calc_E4
  end interface
contains
  subroutine calc_E2(id, nx, gamma, rho, u, p, E)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx) :: rho, u, p
    real(8), intent(out), dimension(nx-1,3) :: E
    integer i
    real(8) Rho_m, U_m, P_m
    do i = 1, nx-1
      Rho_m = 0.5d0 * (rho(i) + rho(i+1))
      U_m = 0.5d0 * (u(i) + u(i+1))
      P_m = 0.5d0 * (p(i) + p(i+1))
      E(i,1) = Rho_m * U_m
      E(i,2) = E(i,1) * U_m + P_m
      E(i,3) = E(i,1) * P_m / ((gamma - 1.d0) * Rho_m) &
              & + 0.5d0 * E(i,1) * (u(i) * u(i+1)) &
              & + 0.5d0 * (u(i) * p(i+1) + u(i+1) * p(i))
    enddo
  end subroutine calc_E2

  subroutine calc_E4(id, nx, gamma, rho, u, p, E)
    integer(kind=4), intent(in), value :: id
    integer, intent(in), value :: nx
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx) :: rho, u, p
    real(8), intent(out), dimension(nx-3,3) :: E
    integer i
    real(8), dimension(3) :: RhoU, RhoUU_P, IE, Energy
    real(8), dimension(4) :: rhos, us, ps
    do i = 1, nx-3
      rhos(:) = rho(i:i+3)
      us(:) = u(i:i+3)
      ps(:) = p(i:i+3)
      RhoU(:)  = RhoPhi(rhos(:), us(:))
      RhoUU_P(:) = RhoPhiU(RhoU(:), us(:)) + Phi(ps(:))
      IE(:) = p(i:i+3) / ((gamma - 1.d0) * rho(i:i+3))
      Energy(:) = IE(:) + RhoUPhiPhi(RhoU(:), us(:)) + PhiPsi(us(:), ps(:))
      E(i,1) = Flux(RhoU(:))
      E(i,2) = Flux(RhoUU_P(:))
      E(i,3) = Flux(Energy(:))
    enddo
  end subroutine calc_E4
end module calc_KEEP

