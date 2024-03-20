module calc_slau
  use calc_common
  use calc_physical_quantities
  implicit none
contains
  attributes(device) function SLAU(dim,Ql,Qr,Normal) result(Flux)
    integer, intent(in), value :: dim ! x:1, y:2
    real(8), intent(in), dimension(4), device :: Ql, Qr, Normal
    real(8) :: Flux(4)
    real(8) rhol, rhor, pl, pr, el, er, Hl, Hr, cl, cr, c
    real(8) vn, V_p, V_m, V_bar, V_bar_p, V_bar_m, M_p, M_m, M, x, g, dp, mass, beta_p, beta_m, Pressure
    real(8), dimension(2) :: Vl, Vr
    real(8), dimension(4) :: phil, phir

    call set_q(Ql,Qr,rhol,rhor,pl,pr,Vl,Vr)

    el = energy(pl,rhol,Vl(1),Vl(2))
    er = energy(pr,rhor,Vr(1),Vr(2))
    Hl = ENTHALPY(el,pl,rhol)
    Hr = ENTHALPY(er,pr,rhor)
    cl = speed_of_sound(pl,rhol)
    cr = speed_of_sound(pr,rhor)
    c = 0.5d0 * (cl + cr)
    ! What is vn?
    vn = 0.d0
    V_p = Vl(dim) - vn
    V_m = Vr(dim) - vn
    M_p = V_p / c
    M_m = V_m / c
    M = min(1.d0, sqrt(0.5d0 * (M_p ** 2 + M_m ** 2)))
    x = (1.d0 - M) ** 2
    g = -max(min(M_p, 0.d0), -1.d0) * min(max(M_m, 0.d0), 1.d0)
    V_bar = (rhol * abs(V_p) + rhor * abs(V_m)) / (rhol + rhor)
    V_bar_p = abs((1.d0 - g) * V_bar + g * V_p)
    V_bar_m = abs((1.d0 - g) * V_bar + g * V_m)
    dp = pr - pl
    mass = 0.5d0 * (rhol * (Vl(dim) + V_bar_p) + rhor * (Vr(dim) - V_bar_m) - x * dp / c)
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
    Pressure = 0.5d0 * (pl + pr + (beta_p - beta_m) * (pl - pr) + (1.d0 - x) * (beta_p + beta_m - 1.d0) * (pl + pr))
    phil(:) = (/1.d0, Vl(1), Vr(2), Hl/)
    phir(:) = (/1.d0, Vr(1), Vr(2), Hr/)
    Flux(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + Pressure * Normal(:)
  end function SLAU
end module calc_slau

