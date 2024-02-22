module calc_term
  implicit none
contains
  attributes(device) subroutine calc_Phi(a, Phi)
    real(8), intent(in) :: a(4)
    real(8), intent(out) :: Phi(3)
    Phi(1) = 0.5d0 * (a(2) + a(3))
    Phi(2) = 0.5d0 * (a(2) + a(4))
    Phi(3) = 0.5d0 * (a(1) + a(3))
  end subroutine calc_Phi
  
  attributes(device) subroutine calc_RhoPhi(rho, u, RhoPhi)
    real(8), intent(in), dimension(4) :: rho, u
    real(8), intent(out) :: RhoPhi(3)
    RhoPhi(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))
    RhoPhi(2) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    RhoPhi(3) = 0.25d0 * (rho(1) + rho(3)) * (u(1) + u(3))
  end subroutine calc_RhoPhi
  
  attributes(device) subroutine calc_RhoPhiU(rhou, ph, RhoPhi)
    real(8), intent(in), dimension(3) :: rhou
    real(8), intent(in), dimension(4) :: ph
    real(8), intent(out) :: RhoPhi(3)
    RhoPhi(1) = rhou(1) * 0.5d0 * (ph(2) + ph(3))
    RhoPhi(2) = rhou(2) * 0.5d0 * (ph(2) + ph(4))
    RhoPhi(3) = rhou(3) * 0.5d0 * (ph(1) + ph(3))
  end subroutine calc_RhoPhiU
  
  attributes(device) subroutine calc_RhoUPhi2(rhou, u, v, w, RhoUPhi2)
    real(8), intent(in), dimension(3) :: rhou
    real(8), intent(in), dimension(4) :: u, v, w
    real(8), intent(out) :: RhoUPhi2(3)
    RhoUPhi2(1) = rhou(1) * 0.5d0 * (u(2) * u(3) + v(2) * v(3) + w(2) * w(3))
    RhoUPhi2(2) = rhou(2) * 0.5d0 * (u(2) * u(4) + v(2) * v(4) + w(2) * w(4))
    RhoUPhi2(3) = rhou(3) * 0.5d0 * (u(1) * u(3) + v(1) * v(3) + w(1) * w(3))
  end subroutine calc_RhoUPhi2
  
  attributes(device) subroutine calc_PhiPsi(ph, psi, PhiPsi)
    real(8), intent(in), dimension(4) :: ph, psi
    real(8), intent(out) :: PhiPsi(3)
    PhiPsi(1) = 0.5d0 * (ph(2) * psi(3) + ph(3) * psi(2))
    PhiPsi(2) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    PhiPsi(3) = 0.5d0 * (ph(1) * psi(3) + ph(3) * psi(1))
  end subroutine calc_PhiPsi
  
  attributes(device) subroutine merge_term(ph, flux)
    real(8), intent(in) :: ph(3)
    real(8), intent(out) :: flux
    flux = 2.d0 * ((2.d0/3.d0) * ph(1) - (ph(2) + ph(3)) / 12.d0)
  end subroutine merge_term
end module