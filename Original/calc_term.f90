module calc_term
  implicit none

  interface Phi
    module procedure Phi2, Phi4
  end interface

  interface RhoPhi
    module procedure RhoPhi2, RhoPhi4
  end interface

  interface RhoPhiU
    module procedure RhoPhiU2, RhoPhiU4
  end interface

  interface RhoUPhiPhi
    module procedure RhoUPhiPhi2, RhoUPhiPhi4
  end interface

  interface PhiPsi
    module procedure PhiPsi2, PhiPsi4
  end interface

  interface Flux
    module procedure Flux2, Flux4
  end interface

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Phi2(a, ans)
    real(8), intent(in), dimension(2) :: a
    real(8), intent(out) :: ans
    ans = 0.5d0 * (a(1) + a(2))
  end subroutine Phi2

  subroutine Phi4(a, ans)
    real(8), intent(in), dimension(4) :: a
    real(8), intent(out), dimension(3) :: ans
    ans(1) = 0.5d0 * (a(2) + a(3))
    ans(2) = 0.5d0 * (a(2) + a(4))
    ans(3) = 0.5d0 * (a(1) + a(3))
  end subroutine Phi4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine RhoPhi2(rho, u, ans)
    real(8), intent(in), dimension(2) :: rho, u
    real(8), intent(out) :: ans
    ans = 0.25d0 * (rho(1) + rho(2)) * (u(1) + u(2))
  end subroutine RhoPhi2

  subroutine RhoPhi4(rho, u, ans)
    real(8), intent(in), dimension(4) :: rho, u
    real(8), intent(out), dimension(3) :: ans
    ans(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))
    ans(2) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    ans(3) = 0.25d0 * (rho(1) + rho(3)) * (u(1) + u(3))
  end subroutine RhoPhi4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function RhoPhiU2(rhou, ph) result(ans)
    real(8), intent(in) :: rhou
    real(8), intent(in), dimension(2) :: ph
    real(8) :: ans
    ans = rhou * 0.5d0 * (ph(1) + ph(2))
  end function RhoPhiU2

  function RhoPhiU4(rhou, ph) result(ans)
    real(8), intent(in), dimension(3) :: rhou
    real(8), intent(in), dimension(4) :: ph
    real(8), dimension(3) :: ans
    ans(1) = rhou(1) * 0.5d0 * (ph(2) + ph(3))
    ans(2) = rhou(2) * 0.5d0 * (ph(2) + ph(4))
    ans(3) = rhou(3) * 0.5d0 * (ph(1) + ph(3))
  end function RhoPhiU4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function RhoUPhiPhi2(rhou, u, v, w) result(ans)
    real(8), intent(in) :: rhou
    real(8), intent(in), dimension(2) :: u, v, w
    real(8) :: ans
    ans = rhou * 0.5d0 * (u(1) * u(2) + v(1) * v(2) + w(1) * w(2))
  end function RhoUPhiPhi2

  function RhoUPhiPhi4(rhou, u, v, w) result(ans)
    real(8), intent(in), dimension(3) :: rhou
    real(8), intent(in), dimension(4) :: u, v, w
    real(8), dimension(3) :: ans
    ans(1) = rhou(1) * 0.5d0 * (u(2) * u(3) + v(2) * v(3) + w(2) * w(3))
    ans(2) = rhou(2) * 0.5d0 * (u(2) * u(4) + v(2) * v(4) + w(2) * w(4))
    ans(3) = rhou(3) * 0.5d0 * (u(1) * u(3) + v(1) * v(3) + w(1) * w(3))
  end function RhoUPhiPhi4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PhiPsi2(ph, psi, ans)
    real(8), intent(in), dimension(2) :: ph, psi
    real(8), intent(out) :: ans
    ans = 0.5d0 * (ph(1) * psi(2) + ph(2) * psi(1))
  end subroutine PhiPsi2

  subroutine PhiPsi4(ph, psi, ans)
    real(8), intent(in), dimension(4) :: ph, psi
    real(8), intent(out), dimension(3) :: ans
    ans(1) = 0.5d0 * (ph(2) * psi(3) + ph(3) * psi(2))
    ans(2) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    ans(3) = 0.5d0 * (ph(1) * psi(3) + ph(3) * psi(1))
  end subroutine PhiPsi4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function Flux2(ph) result(ans)
    real(8), intent(in) :: ph
    real(8) :: ans
    ans = ph
  end function Flux2

  function Flux4(ph) result(ans)
    real(8), intent(in), dimension(3) :: ph
    real(8) :: ans
    ans = 2.d0 * ((2.d0/3.d0) * ph(1) - (ph(2) + ph(3)) / 12.d0)
  end function Flux4
end module