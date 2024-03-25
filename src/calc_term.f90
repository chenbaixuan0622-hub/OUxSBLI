module calc_term
  use mod_globals, only : dimension
  use calc_common_dim
  implicit none
contains
! only 3d version is available

  attributes(device) function Phi(a) result(ans)
    real(8), intent(in), dimension(4), device :: a
    real(8), dimension(3) :: ans
    ans(1) = 0.5d0 * (a(2) + a(3))
    ans(2) = 0.5d0 * (a(2) + a(4))
    ans(3) = 0.5d0 * (a(1) + a(3))
  end function Phi

  attributes(device) function RhoPhi(rho, u) result(ans)
    real(8), intent(in), dimension(4), device :: rho, u
    real(8), dimension(3) :: ans
    ans(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))
    ans(2) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    ans(3) = 0.25d0 * (rho(1) + rho(3)) * (u(1) + u(3))
  end function RhoPhi

  attributes(device) function RhoPhiU(rhou, ph) result(ans)
    real(8), intent(in), dimension(3), device :: rhou
    real(8), intent(in), dimension(4), device :: ph
    real(8), dimension(3) :: ans
    ans(1) = rhou(1) * 0.5d0 * (ph(2) + ph(3))
    ans(2) = rhou(2) * 0.5d0 * (ph(2) + ph(4))
    ans(3) = rhou(3) * 0.5d0 * (ph(1) + ph(3))
  end function RhoPhiU

  attributes(device) function RhoUPhiPhi(rhou, V) result(ans)
    real(8), intent(in), dimension(3), device :: rhou
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), dimension(dimension) :: V1, V2, V3, V4
    real(8), dimension(3) :: ans
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    V3(:) = V(3,:)
    V4(:) = V(4,:)
    ans(1) = rhou(1) * 0.5d0 * vecsum(V2, V3)
    ans(2) = rhou(2) * 0.5d0 * vecsum(V2, V4)
    ans(3) = rhou(3) * 0.5d0 * vecsum(V1, V3)
  end function RhoUPhiPhi

  attributes(device) function PhiPsi(ph, psi) result(ans)
    real(8), intent(in), dimension(4), device :: ph, psi
    real(8), dimension(3) :: ans
    ans(1) = 0.5d0 * (ph(2) * psi(3) + ph(3) * psi(2))
    ans(2) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    ans(3) = 0.5d0 * (ph(1) * psi(3) + ph(3) * psi(1))
  end function PhiPsi

  attributes(device) function Flux(ph) result(ans)
    real(8), intent(in), dimension(3), device :: ph
    real(8) :: ans
    ans = 2.d0 * ((2.d0/3.d0) * ph(1) - (ph(2) + ph(3)) / 12.d0)
  end function Flux
end module

