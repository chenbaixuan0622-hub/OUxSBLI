module calc_keep
  use mod_globals, only : dimension, gamma
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
    real(8), intent(in), dimension(dimension) :: Normal
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), dimension(dimension+2) :: F
    real(8) Rho_m, P_m, PRho
    real(8), dimension(dimension) :: V_m
    Rho_m = 0.5d0 * sum(rho(:))
    V_m(:) = 0.5d0 * sum(V(:,:), dim=1)
    P_m = 0.5d0 * sum(p(:))
    PRho = 0.5d0 * sum(p(:) / rho(:))
    F(1) = Rho_m * V_m(id)
    F(2:dimension+1) = F(1) * V_m(:) + P_m * Normal(:)
    F(dimension+2) = F(1) * PRho / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * sum(V(1,:) * V(2,:)) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1))
  end function KEEP_2nd

  attributes(device) function KEEP_4th(id_accuracy,id,rho,p,V,Normal) result(F)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value :: id
    real(8), intent(in), dimension(4), device :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension), device :: Normal
    real(8), dimension(4) :: PRho
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3) :: RhoV, IE, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    integer i
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i)
    enddo
    PRho(:) = p(:) / rho(:)
    IE(:) = Phi(PRho(:)) / (gamma - 1.d0)
    Energy(:) = IE(:) + RhoUPhiPhi(RhoV(:), V(:,:)) + PhiPsi(V(:,id), p(:))
    F(1) = Flux(RhoV(:))
    do i = 1, dimension
      F(i+1) = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))
  end function KEEP_4th
end module calc_keep

