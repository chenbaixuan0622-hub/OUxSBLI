module calc_keep
  use mod_globals, only : dimension, gamma
  use calc_common_dim
  use calc_term
  use calc_mat
  implicit none
contains
  attributes(device) function KEEP2(id,rho,p,V,Normal) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2) :: F
    real(8)                         :: Rhom, Pm, PRho
    real(8), dimension(dimension)   :: Vm, V1, V2
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    Rhom  = 0.5d0 * (rho(1) + rho(2))
    Vm(:) = 0.5d0 * (V1(:) + V2(:))
    Pm    = 0.5d0 * (p(1) + p(2))
    PRho  = 0.5d0 * (p(1) / rho(1) + p(2) / rho(2))
    F(1)  = Rhom * Vm(id)
    F(2:dimension+1) = F(1) * Vm(:) + Pm * Normal(2:dimension+1)
    !F(dimension+2)   = F(1) * PRho / (gamma - 1.d0) &
    ! KEEPPE
    F(dimension+2)   = Vm(id) * Pm / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * vecsum(V1, V2) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1))
  end function KEEP2

  attributes(device) function KEEPUP(id,rho,p,V,rholr,plr,Vlr,Normal,sensor) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(2), device           :: rholr, plr
    real(8), intent(in), dimension(2,dimension), device :: Vlr
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), intent(in), value                          :: sensor
    real(8), dimension(4)           :: P_over_Rho
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3)           :: Vm, RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    integer i
    !real(8), dimension(dimension+2) :: F
    !real(8)                         :: Rhom, Pm, PRho
    ! upwind
    !real(8) Vroe, dV(dimension)
    ! SLAU
    real(8) press, cl, cr, c, Mp, Mm, M, x, bp, bm, V1(dimension), V2(dimension)
    ! SLAU pressure term
    cl = sqrt(gamma * plr(1) / rholr(1))
    cr = sqrt(gamma * plr(2) / rholr(2))
    c  = 0.5d0 * (cl + cr)
    V1 = Vlr(1,:)
    V2 = Vlr(2,:)
    Mp = V1(id) / c
    Mm = V2(id) / c
    if (abs(Mp) < 1.d0) then
      bp = 0.25d0 * (2.d0 - Mp) * (Mp + 1.d0) ** 2
    else
      bp = 0.5d0 * (1.d0 + sign(1.d0, Mp))
    endif
    if (abs(Mm) < 1.d0) then
      bm = 0.25d0 * (2.d0 + Mm) * (Mm - 1.d0) ** 2
    else
      bm = 0.5d0 * (1.d0 + sign(1.d0, -Mm))
    endif
    M = min(1.d0, sqrt(0.5d0 * q2(V1, V2)) / c)
    x = (1.d0 - M)**2
    !press = 0.5d0 * (plr(1) + plr(2) + (bp - bm) * (plr(1) - plr(2)) + (1.d0 - x) * (bp + bm - 1.d0) * (plr(1) + plr(2)))
    press = 0.5d0 * (plr(1) + plr(2))

    RhoV(:) = RhoPhi(rho(:), V(:,id))
    ! KEEP
    !P_over_Rho(:) = p(:) / rho(:)
    !RhoVIE(:)     = RhoPhiU(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)

    ! KEEP PE
    Vm(:) = Phi(V(:,id))
    RhoVIE(:) = RhoPhiU(Vm(:), p(:)) / (gamma - 1.d0)    

    RhoVKE(:) = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:)     = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1)      = Flux(RhoV(:))
    do i = 1, dimension
      !RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i+1)
      !F(i+1)       = Flux(RhoVV_P(:,i))
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i))
      F(i+1)       = Flux(RhoVV_P(:,i)) + press * Normal(i+1)
    enddo
    F(dimension+2) = Flux(Energy(:))
  end function KEEPUP

  attributes(device) function KEEP4(id,rho,p,V,Normal) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(4)           :: P_over_Rho
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3)           :: Vm, RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    integer i
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    ! KEEP
    !P_over_Rho(:) = p(:) / rho(:)
    !RhoVIE(:)     = RhoPhiU(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)

    ! KEEP PE
    Vm(:) = Phi(V(:,id))
    RhoVIE(:) = RhoPhiU(Vm(:), p(:)) / (gamma - 1.d0)    

    RhoVKE(:) = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:)     = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1)      = Flux(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i+1)
      F(i+1)       = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))
  end function KEEP4

  attributes(device) function KEEPRho(id,rho,p,V,Normal,rho2,p2,V2,sensor) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), intent(in), dimension(2), device           :: rho2, p2
    real(8), intent(in), dimension(2,dimension), device :: V2
    real(8), intent(in), value                          :: sensor
    real(8), dimension(4)           :: P_over_Rho
    real(8), dimension(dimension+2) :: F, Ql, Qr, dQ
    real(8), dimension(3)           :: RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    real(8), dimension(dimension)   :: Vl, Vr, V_ave
    integer i
    real(8) rho_ave, Hl, Hr, H_ave, c_ave, el, er, mat(dimension+2,dimension+2)
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    ! energy equation
    P_over_Rho(:) = p(:) / rho(:)
    RhoVIE(:)     = RhoPhiU(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)
    RhoVKE(:)     = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:)         = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1) = Flux(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i+1)
      F(i+1)       = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))

    Vl = V2(1,:)
    Vr = V2(2,:)
    el = p2(1) / (gamma - 1.d0) + 0.5d0 * rho2(1) * vecsum(Vl(:), Vl(:))
    er = p2(2) / (gamma - 1.d0) + 0.5d0 * rho2(2) * vecsum(Vr(:), Vr(:))
    Hl = (el + p2(1)) / rho2(1)
    Hr = (er + p2(2)) / rho2(2)

    rho_ave  = sqrt(rho(1) * rho(2))
    V_ave(:) = (sqrt(rho(1)) * Vl(:) + sqrt(rho(2)) * Vr(:)) / (sqrt(rho(1)) + sqrt(rho(2)))
    H_ave    = (sqrt(rho(1)) * Hl    + sqrt(rho(2)) * Hr)    / (sqrt(rho(1)) + sqrt(rho(2)))
    c_ave    = sqrt((gamma - 1.d0) * (H_ave - 0.5d0 * vecsum(V_ave(:), V_ave(:))))

    Ql(1)             = rho(1)
    Ql(2:dimension+1) = rho(1) * Vl(:)
    Ql(dimension+2)   = el
    Qr(1)             = rho(2)
    Qr(2:dimension+1) = rho(2) * Vr(:)
    Qr(dimension+2)   = er

    !mat(:,:) = calc_AB(id, rho_ave, H_ave, c_ave, V_ave)
    dQ(:)    = Qr(:) - Ql(:)
    F(:)     = F(:) - 0.5d0 * sensor * cumatmul(mat(:,:), dQ(:))
  end function KEEPRho
end module calc_keep

