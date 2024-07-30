module calc_term
  use mod_globals, only : dimension
  use calc_common_dim
  implicit none
contains
  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function Phi4(a) result(ans)
    real(8), intent(in), dimension(4), device :: a
    real(8), dimension(3) :: ans
    ans(1) = 0.5d0 * (a(2) + a(3))
    ans(2) = 0.5d0 * (a(2) + a(4))
    ans(3) = 0.5d0 * (a(1) + a(3))
  end function Phi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function Phi6(a) result(ans)
    real(8), intent(in), dimension(6), device :: a
    real(8), dimension(6) :: ans
    ans(1) = 0.5d0 * (a(3) + a(4))
    ans(2) = 0.5d0 * (a(3) + a(5))
    ans(3) = 0.5d0 * (a(2) + a(4))
    ans(4) = 0.5d0 * (a(3) + a(6))
    ans(5) = 0.5d0 * (a(2) + a(5))
    ans(6) = 0.5d0 * (a(1) + a(4))
  end function Phi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function RhoPhi4(rho, u) result(ans)
    real(8), intent(in), dimension(4), device :: rho, u
    real(8), dimension(3) :: ans
    ans(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))
    ans(2) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    ans(3) = 0.25d0 * (rho(1) + rho(3)) * (u(1) + u(3))
  end function RhoPhi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  attributes(device) function RhoPhi6(rho, u) result(ans)
    real(8), intent(in), dimension(6), device :: rho, u
    real(8), dimension(6) :: ans
    ans(1) = 0.25d0 * (rho(3) + rho(4)) * (u(3) + u(4))
    ans(2) = 0.25d0 * (rho(3) + rho(5)) * (u(3) + u(5))
    ans(3) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    ans(4) = 0.25d0 * (rho(3) + rho(6)) * (u(3) + u(6))
    ans(5) = 0.25d0 * (rho(2) + rho(5)) * (u(2) + u(5))
    ans(6) = 0.25d0 * (rho(1) + rho(4)) * (u(1) + u(4))
  end function RhoPhi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function RhoPhiU4(rhou, ph) result(ans)
    real(8), intent(in), dimension(3), device :: rhou
    real(8), intent(in), dimension(4), device :: ph
    real(8), dimension(3) :: ans
    ans(1) = rhou(1) * 0.5d0 * (ph(2) + ph(3))
    ans(2) = rhou(2) * 0.5d0 * (ph(2) + ph(4))
    ans(3) = rhou(3) * 0.5d0 * (ph(1) + ph(3))
  end function RhoPhiU4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function RhoPhiU6(rhou, ph) result(ans)
    real(8), intent(in), dimension(6), device :: rhou
    real(8), intent(in), dimension(6), device :: ph
    real(8), dimension(6) :: ans
    ans(1) = rhou(1) * 0.5d0 * (ph(3) + ph(4))
    ans(2) = rhou(2) * 0.5d0 * (ph(3) + ph(5))
    ans(3) = rhou(3) * 0.5d0 * (ph(2) + ph(4))
    ans(4) = rhou(4) * 0.5d0 * (ph(3) + ph(6))
    ans(5) = rhou(5) * 0.5d0 * (ph(2) + ph(5))
    ans(6) = rhou(6) * 0.5d0 * (ph(1) + ph(4))
  end function RhoPhiU6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function RhoUPhiPhi4(rhou, V) result(ans)
    real(8), intent(in), dimension(3), device           :: rhou
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
  end function RhoUPhiPhi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function RhoUPhiPhi6(rhou, V) result(ans)
    real(8), intent(in), dimension(6), device           :: rhou
    real(8), intent(in), dimension(6,dimension), device :: V
    real(8), dimension(dimension) :: V1, V2, V3, V4, V5, V6
    real(8), dimension(6) :: ans
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    V3(:) = V(3,:)
    V4(:) = V(4,:)
    V5(:) = V(5,:)
    V6(:) = V(6,:)
    ans(1) = rhou(1) * 0.5d0 * vecsum(V3, V4)
    ans(2) = rhou(2) * 0.5d0 * vecsum(V3, V5)
    ans(3) = rhou(3) * 0.5d0 * vecsum(V2, V4)
    ans(4) = rhou(4) * 0.5d0 * vecsum(V3, V6)
    ans(5) = rhou(5) * 0.5d0 * vecsum(V2, V5)
    ans(6) = rhou(6) * 0.5d0 * vecsum(V1, V4)
  end function RhoUPhiPhi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function PhiPsi4(ph, psi) result(ans)
    real(8), intent(in), dimension(4), device :: ph, psi
    real(8), dimension(3) :: ans
    ans(1) = 0.5d0 * (ph(2) * psi(3) + ph(3) * psi(2))
    ans(2) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    ans(3) = 0.5d0 * (ph(1) * psi(3) + ph(3) * psi(1))
  end function PhiPsi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function PhiPsi6(ph, psi) result(ans)
    real(8), intent(in), dimension(6), device :: ph, psi
    real(8), dimension(6) :: ans
    ans(1) = 0.5d0 * (ph(3) * psi(4) + ph(4) * psi(3))
    ans(2) = 0.5d0 * (ph(3) * psi(5) + ph(5) * psi(3))
    ans(3) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    ans(4) = 0.5d0 * (ph(3) * psi(6) + ph(6) * psi(3))
    ans(5) = 0.5d0 * (ph(2) * psi(5) + ph(5) * psi(2))
    ans(6) = 0.5d0 * (ph(1) * psi(4) + ph(4) * psi(1))
  end function PhiPsi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function Flux4(ph) result(ans)
    real(8), intent(in), dimension(3), device :: ph
    real(8) :: ans
    ans = 2.d0 * ((2.d0/3.d0) * ph(1) - (ph(2) + ph(3)) / 12.d0)
  end function Flux4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  attributes(device) function Flux6(ph) result(ans)
    real(8), intent(in), dimension(6), device :: ph
    real(8) :: ans
    ans = 2.d0 * (0.75d0 * ph(1) - 3.d0 * (ph(2) + ph(3)) / 20.d0 &
          + (ph(4) + ph(5) + ph(6)) / 60.d0)
  end function Flux6
end module

