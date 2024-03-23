module calc_slau
  use mod_globals, only : dimension
  use calc_common
  use calc_physical_quantities
  implicit none
contains
  attributes(device) subroutine calc_quantities_AUSM(pl,pr,rhol,rhor,Vl,Vr,el,er,Hl,Hr,cl,cr)
    real(8), intent(in), value :: pl, pr, rhol, rhor
    real(8), intent(in), dimension(dimension), device :: Vl, Vr
    real(8), intent(out) :: el, er, Hl, Hr, cl, cr
    el = energy(pl,rhol,Vl(:))
    er = energy(pr,rhor,Vr(:))
    Hl = ENTHALPY(el,pl,rhol)
    Hr = ENTHALPY(er,pr,rhor)
    cl = speed_of_sound(pl,rhol)
    cr = speed_of_sound(pr,rhor)
  end subroutine calc_quantities_AUSM

  attributes(device) subroutine calc_beta(M_p,M_m,beta_p,beta_m)
    real(8), intent(in), value :: M_p, M_m
    real(8), intent(out) :: beta_p, beta_m
    if (abs(M_p) < 1.d0) then
      beta_p = 0.25d0 * (2.d0 - M_p) * (M_p + 1.d0) ** 2
    else
      beta_p = 0.5d0 * (1.d0 + sign(1.d0, M_p))
    endif
    if (abs(M_m) < 1.d0) then
      beta_m = 0.25d0 * (2.d0 + M_m) * (M_m - 1.d0) ** 2
    else
      beta_m = 0.5d0 * (1.d0 + sign(1.d0, -M_m))
    endif
  end subroutine calc_beta

  attributes(device) function flux_AUSM(mass,pressure,Hl,Hr,Vl,Vr,Normal) result(F)
    real(8), intent(in), value :: mass, pressure, Hl, Hr
    real(8), intent(in), dimension(dimension), device :: Vl, Vr
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2) :: F, phil, phir
    phil(1) = 1.d0
    phil(2:dimension+1) = Vl(:)
    phil(dimension+2) = Hl
    phir(1) = 1.d0
    phir(2:dimension+1) = Vr(:)
    phir(dimension+2) = Hr
    F(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + pressure * Normal(:)
  end function flux_AUSM

  attributes(device) function SLAU(id_dim,Ql,Qr,Normal) result(F)
    integer, intent(in), value :: id_dim
    real(8), intent(in), dimension(dimension+2), device :: Ql, Qr, Normal
    real(8) rhol, rhor, pl, pr, el, er, Hl, Hr, cl, cr, c
    real(8) vn, V_p, V_m, V_bar, V_bar_p, V_bar_m, M_p, M_m, M, x, g, dp, mass, beta_p, beta_m, Pressure
    real(8), dimension(dimension) :: Vl, Vr
    real(8), dimension(dimension+2) :: F

    call set_q(Ql,Qr,rhol,rhor,pl,pr,Vl,Vr)
    call calc_quantities_AUSM(pl,pr,rhol,rhor,Vl,Vr,el,er,Hl,Hr,cl,cr)

    c = 0.5d0 * (cl + cr)
    ! What is vn?
    vn = 0.d0
    V_p = Vl(id_dim) - vn
    V_m = Vr(id_dim) - vn
    M_p = V_p / c
    M_m = V_m / c
    M = min(1.d0, sqrt(0.5d0 * (M_p ** 2 + M_m ** 2)))
    x = (1.d0 - M) ** 2
    g = -max(min(M_p, 0.d0), -1.d0) * min(max(M_m, 0.d0), 1.d0)
    V_bar = (rhol * abs(V_p) + rhor * abs(V_m)) / (rhol + rhor)
    V_bar_p = abs((1.d0 - g) * V_bar + g * V_p)
    V_bar_m = abs((1.d0 - g) * V_bar + g * V_m)
    dp = pr - pl
    call calc_beta(M_p,M_m,beta_p,beta_m)
    mass = 0.5d0 * (rhol * (Vl(id_dim) + V_bar_p) + rhor * (Vr(id_dim) - V_bar_m) - x * dp / c)
    Pressure = 0.5d0 * (pl + pr + (beta_p - beta_m) * (pl - pr) + (1.d0 - x) * (beta_p + beta_m - 1.d0) * (pl + pr))
    F(:) = flux_AUSM(mass,Pressure,Hl,Hr,Vl,Vr,Normal)
  end function SLAU
end module calc_slau

