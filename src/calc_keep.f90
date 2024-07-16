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
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), dimension(dimension+2) :: F
    real(8)                         :: Rho_m, P_m, PRho
    real(8), dimension(dimension)   :: V_m, V1, V2
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    Rho_m = 0.5d0 * (rho(1) + rho(2))
    ! contravariant velocity
    V_m(:) = 0.5d0 * (V1(:) + V2(:)) * Normal(id)
    P_m    = 0.5d0 * (p(1) + p(2))
    PRho   = 0.5d0 * (p(1) / rho(1) + p(2) / rho(2))
    F(1)   = Rho_m * V_m(id)
    F(2:dimension+1) = F(1) * V_m(:) + P_m * Normal(:)
    F(dimension+2)   = F(1) * PRho / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * vecsum(V1, V2) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1)) * Normal(id)
  end function KEEP2

  attributes(device) function minmod(x,y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = dsign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod

  attributes(device) subroutine MUSCL3(k,a,al,ar)
    real(8), intent(in), value  :: k
    real(8), intent(in), device :: a(4)
    real(8), intent(out)        :: al, ar
    real(8) b, d1, d2, d3
    b = (3.d0 - k) / (1.d0 - k)
    d1 = -a(1) + a(2)
    d2 = -a(2) + a(3)
    d3 = -a(3) + a(4)
    al = a(2) + 0.25d0 * ((1.d0 - k) * d1 + (1.d0 + k) * d2)
    ar = a(3) - 0.25d0 * ((1.d0 - k) * d3 + (1.d0 + k) * d2)
  end subroutine MUSCL3

  attributes(device) function KEEP4(id,rho,p,V,Normal) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), dimension(4)           :: P_over_Rho
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3)           :: Vm, RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    ! MUSCL interpolation
    ! k = 0.d0 2nd-upwind + 1, k = 1.d0 / 3.d0 3rd
    real(8) :: pr, pl, k = 1.d0 / 3.d0, b, d1, d2, d3
    ! Van Albada limiter
    real(8) :: s, eps = 1.d-4
    integer i
    RhoV(:)       = RhoPhi(rho(:), V(:,id))
    ! KEEP
    P_over_Rho(:) = p(:) / rho(:)
    RhoVIE(:)     = RhoPhiU(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)

    ! KEEP PE
    !Vm(:) = Phi(V(:,id))
    !RhoVIE(:) = RhoPhiU(Vm(:), p(:)) / (gamma - 1.d0)    

    ! MUSCL interpolation
    b  = (3.d0 - k) / (1.d0 - k) 
    d1 = -p(1) + p(2)
    d2 = -p(2) + p(3)
    d3 = -p(3) + p(4)
    ! TVD minmod
    !pl = p(2) + 0.25d0 * ((1.d0 - k) * minmod(d1, b * d2) + (1.d0 + k) * minmod(d2, b * d1))
    !pr = p(3) - 0.25d0 * ((1.d0 - k) * minmod(d3, b * d2) + (1.d0 + k) * minmod(d2, b * d3))
    ! TVD van Albada
    !s  = (2.d0 * d1 * d2 + eps) / (d1**2 + d2**2 + eps)
    !pl = p(2) + 0.25d0 * s * ((1.d0 - k * s) * d1 + (1.d0 + k * s) * d2)
    !s  = (2.d0 * d2 * d3 + eps) / (d2**2 + d3**2 + eps)
    !pr = p(3) - 0.25d0 * s * ((1.d0 - k * s) * d3 + (1.d0 + k * s) * d2)
    ! nonTVD
    pl = p(2) + 0.25d0 * ((1.d0 - k) * d1 + (1.d0 + k) * d2)
    pr = p(3) - 0.25d0 * ((1.d0 - k) * d3 + (1.d0 + k) * d2)
    
    RhoVKE(:) = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:)     = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1)      = Flux(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + 0.5d0 * (pl + pr) * Normal(i)!Phi(p(:)) * Normal(i)
      F(i+1)       = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))
  end function KEEP4

  attributes(device) function KEEPRho(id,rho,p,V,Normal,fd) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), intent(in), value                          :: fd
    real(8), dimension(4)           :: P_over_Rho
    real(8), dimension(dimension+2) :: F, Ql, Qr, dQ
    real(8), dimension(3)           :: RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    real(8), dimension(dimension)   :: Vl, Vr, V_ave
    integer i
    real(8) rhol, rhor, rho_ave, pl, pr, Hl, Hr, H_ave, c_ave, el, er, mat(dimension+2,dimension+2)
    real(8) :: k = 1.d0 / 3.d0
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    ! energy equation
    P_over_Rho(:) = p(:) / rho(:)
    RhoVIE(:)     = RhoPhiU(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)
    RhoVKE(:)     = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:)         = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1) = Flux(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i)
      F(i+1)       = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))

    ! Rho viscosity
    call MUSCL3(k,rho,rhol,rhor)
    call MUSCL3(k,p,pl,pr)
    do i = 1, dimension
      call MUSCL3(k,V(:,i),Vl(i),Vr(i))
    enddo
    el = pl / (gamma - 1.d0) + 0.5d0 * rhol * vecsum(Vl(:), Vl(:))
    er = pr / (gamma - 1.d0) + 0.5d0 * rhor * vecsum(Vr(:), Vr(:))
    Hl = (el + pl) / rhol
    Hr = (er + pr) / rhor

    rho_ave  = sqrt(rhol * rhor)
    V_ave(:) = (sqrt(rhol) * Vl(:) + sqrt(rhor) * Vr(:)) / (sqrt(rhol) + sqrt(rhor))
    H_ave    = (sqrt(rhol) * Hl    + sqrt(rhor) * Hr)    / (sqrt(rhol) + sqrt(rhor))
    c_ave    = sqrt((gamma - 1.d0) * (H_ave - 0.5d0 * vecsum(V_ave(:), V_ave(:))))

    Ql(1)             = rhol
    Ql(2:dimension+1) = rhol * Vl(:)
    Ql(dimension+2)   = el
    Qr(1)             = rhor
    Qr(2:dimension+1) = rhor * Vr(:)
    Qr(dimension+2)   = er

    !mat(:,:) = calc_AB(id, rho_ave, H_ave, c_ave, V_ave)
    dQ(:)    = Qr(:) - Ql(:)
    F(:)     = F(:) - 0.5d0 * fd * cumatmul(mat(:,:), dQ(:))
  end function KEEPRho
end module calc_keep

