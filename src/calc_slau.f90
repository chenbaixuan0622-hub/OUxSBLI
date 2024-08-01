module calc_slau
  use mod_globals, only : dimension, gamma
  use calc_common_dim
  use calc_hybrid
  implicit none
  interface SLAU
    module procedure SLAU1, HRSLAU2, VHRSLAU2
  end interface SLAU
contains
  attributes(device) subroutine SLAU_common(id,rho,p,V,c,Mp,Mm,bp,bm,dp,Vtp,Vtm,Vl,Vr,phil,phir)
    integer, intent(in), value                           :: id
    real(8), intent(in), dimension(2), device            :: rho, p
    real(8), intent(in), device                          :: V(2,dimension)
    real(8), intent(out)                                 :: c, Mp, Mm, bp, bm, dp, Vtp, Vtm
    real(8), intent(out), dimension(dimension), device   :: Vl, Vr
    real(8), intent(out), dimension(dimension+2), device :: phil, phir
    real(8) cl, cr, g, Vt, el, er
    cl  = sqrt(gamma * p(1) / rho(1))
    cr  = sqrt(gamma * p(2) / rho(2))
    c   = 0.5d0 * (cl + cr)

    Vl  = V(1,:)
    Vr  = V(2,:)

    Mp  = Vl(id) / c
    Mm  = Vr(id) / c
    
    g   = -max(min(Mp, 0.d0), -1.d0) * min(max(Mm, 0.d0), 1.d0)
    Vt  = (rho(1) * abs(Vl(id)) + rho(2) * abs(Vr(id))) / (rho(1) + rho(2))
    Vtp = (1.d0 - g) * Vt + g * abs(Vl(id))
    Vtm = (1.d0 - g) * Vt + g * abs(Vr(id))

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
    dp = -p(1) + p(2)

    el = p(1) / (gamma - 1.d0) + 0.5d0 * rho(1) * vecsum(Vl(:), Vl(:))
    er = p(2) / (gamma - 1.d0) + 0.5d0 * rho(2) * vecsum(Vr(:), Vr(:))

    phil(1)             = 1.d0
    phil(2:dimension+1) = Vl(:)
    phil(dimension+2)   = (el + p(1)) / rho(1)
    phir(1)             = 1.d0
    phir(2:dimension+1) = Vr(:)
    phir(dimension+2)   = (er + p(2)) / rho(2)
  end subroutine SLAU_common

  attributes(device) function SLAU1(id_slau,id,rho,p,V,Norm,HR,sensor) result(F)
    integer(kind=2), intent(in), value                  :: id_slau
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Norm
    real(8), intent(in), value, optional                :: HR, sensor
    real(8) c, Vl(dimension), Vr(dimension), Mp, Mm, M, x
    real(8) Vtp, Vtm, dp, bp, bm, mass, pres
    real(8), dimension(dimension+2) :: F, phil, phir
    call SLAU_common(id,rho,p,V,c,Mp,Mm,bp,bm,dp,Vtp,Vtm,Vl,Vr,phil,phir)
    M = min(1.d0, sqrt(0.5d0 * q2(Vl(:), Vr(:))) / c)
    x = (1.d0 - M) ** 2

    mass = 0.5d0 * (rho(1) * (Vl(id) + Vtp) + rho(2) * (Vr(id) - Vtm) - x * dp / c)
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
    real(8) c, Vl(dimension), Vr(dimension), Mp, Mm, M, x
    real(8) Vtp, Vtm, dp, bp, bm, mass, pres, V2
    real(8), dimension(dimension+2) :: F, phil, phir
    call SLAU_common(id,rho,p,V,c,Mp,Mm,bp,bm,dp,Vtp,Vtm,Vl,Vr,phil,phir)
    V2 = sqrt(0.5d0 * q2(Vl(:), Vr(:)))
    M  = min(1.d0, V2 / c)
    x  = (1.d0 - M) ** 2

    mass = 0.5d0 * (rho(1) * (Vl(id) + Vtp) + rho(2) * (Vr(id) - Vtm) - x * dp / c)
    pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + HR * V2 * (bp + bm - 1.d0) * 0.5d0 * (rho(1) + rho(2)) * c)
    F(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + pres * Norm(:)
  end function HRSLAU2

  attributes(device) function VHRSLAU2(id_slau,id,rho,p,V,Norm,HR,sensor) result(F)
    integer(kind=8), intent(in), value                  :: id_slau
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Norm
    real(8), intent(in), value                          :: HR, sensor
    real(8) c, Vl(dimension), Vr(dimension), Mp, Mm, M, x
    real(8) Vtp, Vtm, dp, bp, bm, mass, pres, V2
    real(8), dimension(dimension+2) :: F, phil, phir
    call SLAU_common(id,rho,p,V,c,Mp,Mm,bp,bm,dp,Vtp,Vtm,Vl,Vr,phil,phir)
    V2 = sqrt(0.5d0 * q2(Vl(:), Vr(:)))
    M  = min(1.d0, V2 / c)
    x   = (1.d0 - M) ** 2

    mass = 0.5d0 * (rho(1) * (Vl(id) + Vtp) + rho(2) * (Vr(id) - Vtm) - x * dp / c)
    pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + HR * V2 * (bp + bm - 1.d0) * 0.5d0 * (rho(1) + rho(2)) * c)
    F(:) = 0.5d0 * (mass * (phil(:) + phir(:)) - sigmoid(sensor) * abs(mass) * (phir(:) - phil(:))) + pres * Norm(:)
  end function VHRSLAU2
end module calc_slau

