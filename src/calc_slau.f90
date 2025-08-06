module calc_slau
  use mod_globals, only : dimension, gamma
  use mod_constant, only : over_gamma_1
  use calc_common_dim
  use calc_hybrid
  implicit none
  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU
contains
  attributes(device) subroutine SLAU_common(id,rho,p,V,c,over_c,Mp,Mm,bp,bm,dp,Vtp,Vtm,Vl,Vr,phil,phir)
    integer, intent(in), value                           :: id
    real(8), intent(in), dimension(2), device            :: rho, p
    real(8), intent(in), device                          :: V(2,dimension)
    real(8), intent(out)                                 :: c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm
    real(8), intent(out), dimension(dimension), device   :: Vl, Vr
    real(8), intent(out), dimension(dimension+2), device :: phil, phir
    real(8) cl, cr, g, Vt, el, er, over_rho1, over_rho2
    over_rho1 = 1.d0 / rho(1)
    over_rho2 = 1.d0 / rho(2)
    cl  = sqrt(gamma * p(1) * over_rho1)
    cr  = sqrt(gamma * p(2) * over_rho2)
    c   = 0.5d0 * (cl + cr)
    over_c = 1.d0 / c
    Vl  = V(1,:)
    Vr  = V(2,:)
    Mp  = Vl(id) * over_c
    Mm  = Vr(id) * over_c
    g   = -max(min(Mp, 0.d0), -1.d0) * min(max(Mm, 0.d0), 1.d0)
    Vt  = (rho(1) * abs(Vl(id)) + rho(2) * abs(Vr(id))) / (rho(1) + rho(2))
    Vtp = (1.d0 - g) * Vt + g * abs(Vl(id))
    Vtm = (1.d0 - g) * Vt + g * abs(Vr(id))
    bp = merge(0.25d0 * (2.d0 - Mp) * (Mp + 1.d0)**2, &
               0.5d0 * (1.d0 + sign(1.d0, Mp)), &
               abs(Mp) < 1.d0)
    bm = merge(0.25d0 * (2.d0 + Mm) * (Mm - 1.d0)**2, &
               0.5d0 * (1.d0 + sign(1.d0, -Mm)), &
               abs(Mm) < 1.d0)
    dp = -p(1) + p(2)
    el = p(1) * over_gamma_1 + 0.5d0 * rho(1) * vecsum(Vl(:), Vl(:))
    er = p(2) * over_gamma_1 + 0.5d0 * rho(2) * vecsum(Vr(:), Vr(:))
    phil(1)             = 1.d0
    phil(2:dimension+1) = Vl(:)
    phil(dimension+2)   = (el + p(1)) * over_rho1
    phir(1)             = 1.d0
    phir(2:dimension+1) = Vr(:)
    phir(dimension+2)   = (er + p(2)) * over_rho2
  end subroutine SLAU_common

  attributes(device) function SLAU1(id_slau,id,rho,p,V,Norm,HR,sensor) result(F)
    integer(kind=2), intent(in), value                  :: id_slau
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Norm
    real(8), intent(in), value, optional                :: HR, sensor
    real(8) c, over_c, Vl(dimension), Vr(dimension), Mp, Mm, M, x
    real(8) Vtp, Vtm, dp, bp, bm, mass, pres
    real(8), dimension(dimension+2) :: F, phil, phir
    call SLAU_common(id,rho,p,V,c,over_c,Mp,Mm,bp,bm,dp,Vtp,Vtm,Vl,Vr,phil,phir)
    M = min(1.d0, sqrt(0.5d0 * q2(Vl(:), Vr(:))) * over_c)
    x = (1.d0 - M) ** 2
    mass = 0.5d0 * (rho(1) * (Vl(id) + Vtp) + rho(2) * (Vr(id) - Vtm) - x * dp * over_c)
    pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + (1.d0 - x) * (bp + bm - 1.d0) * (p(1) + p(2)))
    F(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + pres * Norm(:)
  end function SLAU1

  attributes(device) function HRSLAU2(id_slau,id,rho,p,V,Norm,HR,sensor) result(F)
    integer(kind=4), intent(in), value                  :: id_slau
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Norm
    real(8), intent(in), value                          :: HR
    real(8), intent(in), value, optional                :: sensor
    real(8) c, over_c, Vl(dimension), Vr(dimension), Mp, Mm, M, x
    real(8) Vtp, Vtm, dp, bp, bm, mass, pres, V2
    real(8), dimension(dimension+2) :: F, phil, phir
    call SLAU_common(id,rho,p,V,c,over_c,Mp,Mm,bp,bm,dp,Vtp,Vtm,Vl,Vr,phil,phir)
    V2 = sqrt(0.5d0 * q2(Vl(:), Vr(:)))
    M  = min(1.d0, V2 * over_c)
    x  = (1.d0 - M) ** 2
    mass = 0.5d0 * (rho(1) * (Vl(id) + Vtp) + rho(2) * (Vr(id) - Vtm) - x * dp * over_c)
    pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + HR * V2 * (bp + bm - 1.d0) * 0.5d0 * (rho(1) + rho(2)) * c)
    F(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + pres * Norm(:)
  end function HRSLAU2
end module calc_slau

