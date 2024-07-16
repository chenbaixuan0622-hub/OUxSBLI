module calc_slau1d
  use mod_globals, only : gamma
  implicit none
contains
  attributes(device) function slau(rho,p,V) result(F)
    real(8), intent(in), dimension(2), device :: rho, p, V 
    real(8) cl, cr, c, Mp, Mm, M, x, g, Vt, Vtp, Vtm, dp, bp, bm, etl, etr, Hl, Hr, mass, pres
    real(8) :: phil(3), phir(3), Norm(3) = (/0.d0, 1.d0, 0.d0/), F(3)
    cl  = sqrt(gamma * p(1) / rho(1))
    cr  = sqrt(gamma * p(2) / rho(2))
    c   = 0.5d0 * (cl + cr)
    Mp  = V(1) / c
    Mm  = V(2) / c
    M   = min(1.d0, sqrt(0.5d0 * (V(1)**2 + V(2)**2)) / c)
    x   = (1.d0 - M)**2
    g   = -max(min(Mp, 0.d0), -1.d0) * min(max(Mm, 0.d0), 1.d0)
    Vt  = (rho(1) * abs(V(1)) + rho(2) * abs(V(2))) / (rho(1) + rho(2))
    Vtp = (1.d0 - g) * Vt + g * abs(V(1))
    Vtm = (1.d0 - g) * Vt + g * abs(V(2))
    dp  = p(2) - p(1)
    if (abs(Mp) < 1.d0) then
      bp = 0.25d0 * (2.d0 - Mp) * (Mp + 1.d0)**2
    else
      bp = 0.5d0 * (1.d0 + sign(1.d0, Mp))
    endif
    if (abs(Mm) < 1.d0) then
      bm = 0.25d0 * (2.d0 + Mm) * (Mm - 1.d0)**2
    else
      bm = 0.5d0 * (1.d0 + sign(1.d0, -Mm))
    endif
    etl = p(1) / (gamma - 1.d0) + 0.5d0 * rho(1) * V(1)**2
    etr = p(2) / (gamma - 1.d0) + 0.5d0 * rho(2) * V(2)**2
    Hl = (etl + p(1)) / rho(1)
    Hr = (etr + p(2)) / rho(2)
    phil(:) = (/1.d0, V(1), Hl/)
    phir(:) = (/1.d0, V(2), Hr/)
    mass = 0.5d0 * (rho(1) * (V(1) + Vtp) + rho(2) * (V(2) - Vtm) - x * dp / c)
    pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + (1.d0 - x) * (bp + bm - 1.d0) * (p(1) + p(2)))
    F(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + pres * Norm(:)
  end function slau
end module calc_slau1d

