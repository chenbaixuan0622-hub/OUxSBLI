module calc_slau
  use mod_globals, only : dimension, gamma
  use calc_common_dim
  implicit none
contains
  attributes(device) function SLAU(id,rho,p,V,Norm) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Norm
    real(8) cl, cr, c, Vl(dimension), Vr(dimension), Mp, Mm, M, x, g
    real(8) Vt, Vtp, Vtm, dp, bp, bm, el, er, mass, pres
    real(8), dimension(dimension+2) :: F, phil, phir
    cl  = sqrt(gamma * p(1) / rho(1))
    cr  = sqrt(gamma * p(2) / rho(2))
    c   = 0.5d0 * (cl + cr)

    Vl  = V(1,:)
    Vr  = V(2,:)

    Mp  = Vl(id) / c
    Mm  = Vr(id) / c
    M   = min(1.d0, sqrt(0.5d0 * q2(Vl(:), Vr(:))) / c)
    x   = (1.d0 - M) ** 2
    g   = -max(min(Mp, 0.d0), -1.d0) * min(max(Mm, 0.d0), 1.d0)
    Vt  = (rho(1) * abs(Vl(id)) + rho(2) * abs(Vr(id))) / (rho(1) + rho(2))
    Vtp = (1.d0 - g) * Vt + g * abs(Vl(id))
    Vtm = (1.d0 - g) * Vt + g * abs(Vr(id))
    dp  = p(2) - p(1)

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

    el = p(1) / (gamma - 1.d0) + 0.5d0 * rho(1) * vecsum(Vl(:), Vl(:))
    er = p(2) / (gamma - 1.d0) + 0.5d0 * rho(2) * vecsum(Vr(:), Vr(:))

    phil(1)             = 1.d0
    phil(2:dimension+1) = Vl(:)
    phil(dimension+2)   = (el + p(1)) / rho(1)
    phir(1)             = 1.d0
    phir(2:dimension+1) = Vr(:)
    phir(dimension+2)   = (er + p(2)) / rho(2)

    mass = 0.5d0 * (rho(1) * (Vl(id) + Vtp) + rho(2) * (Vr(id) - Vtm) - x * dp / c)
    ! OK
    !mass = 0.5d0 * (rho(1) * (Vl(id) + Vtp) + rho(2) * (Vr(id) - Vtm))
    pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + (1.d0 - x) * (bp + bm - 1.d0) * (p(1) + p(2)))
    ! not OK
    !pres = 0.5d0 * (p(1) + p(2))
    F(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + pres * Norm(:) 
  end function SLAU
end module calc_slau

