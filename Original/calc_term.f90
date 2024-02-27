module calc_term
  implicit none
contains

  function Phi(a) result(ans)
    real(8), intent(in), dimension(4) :: a
    real(8), dimension(3) :: ans
    ans(1) = 0.5d0 * (a(2) + a(3))
    ans(2) = 0.5d0 * (a(2) + a(4))
    ans(3) = 0.5d0 * (a(1) + a(3))
  end function Phi

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function RhoPhi(rho, u) result(ans)
    real(8), intent(in), dimension(4) :: rho, u
    real(8), dimension(3) :: ans
    ans(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))
    ans(2) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    ans(3) = 0.25d0 * (rho(1) + rho(3)) * (u(1) + u(3))
  end function RhoPhi

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function RhoPhiU(rhou, ph) result(ans)
    real(8), intent(in), dimension(3) :: rhou
    real(8), intent(in), dimension(4) :: ph
    real(8), dimension(3) :: ans
    ans(1) = rhou(1) * 0.5d0 * (ph(2) + ph(3))
    ans(2) = rhou(2) * 0.5d0 * (ph(2) + ph(4))
    ans(3) = rhou(3) * 0.5d0 * (ph(1) + ph(3))
  end function RhoPhiU

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function RhoUPhiPhi(rhou, u, v, w) result(ans)
    real(8), intent(in), dimension(3) :: rhou
    real(8), intent(in), dimension(4) :: u, v, w
    real(8), dimension(3) :: ans
    ans(1) = rhou(1) * 0.5d0 * (u(2) * u(3) + v(2) * v(3) + w(2) * w(3))
    ans(2) = rhou(2) * 0.5d0 * (u(2) * u(4) + v(2) * v(4) + w(2) * w(4))
    ans(3) = rhou(3) * 0.5d0 * (u(1) * u(3) + v(1) * v(3) + w(1) * w(3))
  end function RhoUPhiPhi

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function PhiPsi(ph, psi) result(ans)
    real(8), intent(in), dimension(4) :: ph, psi
    real(8), dimension(3) :: ans
    ans(1) = 0.5d0 * (ph(2) * psi(3) + ph(3) * psi(2))
    ans(2) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    ans(3) = 0.5d0 * (ph(1) * psi(3) + ph(3) * psi(1))
  end function PhiPsi

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function Flux(ph) result(ans)
    real(8), intent(in), dimension(3) :: ph
    real(8) :: ans
    ans = 2.d0 * ((2.d0/3.d0) * ph(1) - (ph(2) + ph(3)) / 12.d0)
  end function Flux
end module
