module calc_keep
  use mod_globals, only : dimension, gamma
  use calc_common_dim
  use calc_term
  implicit none
  interface KEEP
    module procedure KEEP_2nd, KEEP_4th
  end interface
contains
  attributes(device) function KEEP_2nd(id_accuracy,id,rho,p,V,Normal) result(F)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value :: id
    real(8), intent(in), dimension(2), device :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension), device :: Normal
    real(8), dimension(dimension+2) :: F
    real(8) Rho_m, P_m, PRho
    real(8), dimension(dimension) :: V_m, V1, V2
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    Rho_m = 0.5d0 * (rho(1) + rho(2))
    V_m(:) = 0.5d0 * (V1(:) + V2(:))
    P_m = 0.5d0 * (p(1) + p(2))
    PRho = 0.5d0 * (p(1) / rho(1) + p(2) / rho(2))
    F(1) = Rho_m * V_m(id)
    F(2:dimension+1) = F(1) * V_m(:) + P_m * Normal(:)
    F(dimension+2) = F(1) * PRho / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * vecsum(V1, V2) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1))
  end function KEEP_2nd

  attributes(device) function KEEP_4th(id_accuracy,id,rho,p,V,Normal) result(F)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value :: id
    real(8), intent(in), dimension(4), device :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension), device :: Normal
    real(8), dimension(4) :: P_over_Rho
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3) :: RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    integer i
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    ! energy equation
    P_over_Rho(:) = p(:) / rho(:)
    RhoVIE(:) = RhoPhiU(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)
    RhoVKE(:) = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:) = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1) = Flux(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i)
      F(i+1) = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))
  end function KEEP_4th
end module calc_keep

