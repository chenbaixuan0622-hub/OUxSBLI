module calc_slau
  use calc_Qlr
  use calc_common
  implicit none
contains
  subroutine calc_E(nx, gamma, k, b, rho, u, p, F)
    integer, intent(in) :: nx
    real(8), intent(in) :: gamma, k, b
    real(8), intent(in), dimension(nx) :: rho, u, p
    real(8), intent(out), dimension(nx-1,3) :: F
    integer i
    real(8) rhol, rhor, pl, pr, el, er, hl, hr, mass, Pressure, dp, c, cl, cr
    real(8) vn, x, g, M_p, M_m, M, Vl, Vr, V_p, V_m, V_bar, V_bar_p, V_bar_m, beta_p, beta_m
    real(8), dimension(3) :: d1, d2, d3, phil, phir, Normal
    real(8), dimension(nx,3) :: Q
    Normal(1) = 0.d0
    Normal(2) = 1.d0
    Normal(3) = 0.d0
    Q(:,1) = rho
    Q(:,2) = u
    Q(:,3) = p
    do i = 1, nx-1
      if (2 <= i .and. i <= nx-2) then
        call Qlr(k,b,Q(i-1,:),Q(i,:),Q(i+1,:),Q(i+2,:),rhol,rhor,Vl,Vr,pl,pr)
      elseif (i == 1) then
        call Qlr(k,b,0.d0,Q(i,:),Q(i+1,:),Q(i+2,:),rhol,rhor,Vl,Vr,pl,pr)
      else
        call Qlr(k,b,Q(i-1,:),Q(i,:),Q(i+1,:),0.d0,rhol,rhor,Vl,Vr,pl,pr)
      endif

      ! calc SLAU
      el = energy(gamma,pl,rhol,Vl)
      er = energy(gamma,pr,rhor,Vr)
      hl = enthalpy(el,pl,rhol)
      hr = enthalpy(er,pr,rhor)
      phil(1) = 1.d0
      phil(2) = Vl
      phil(3) = hl
      phir(1) = 1.d0
      phir(2) = Vr
      phir(3) = hr
      cl = speed_of_sound(gamma,pl,rhol)
      cr = speed_of_sound(gamma,pr,rhor)
      c = 0.5d0 * (cl + cr)
      ! what is vn?
      vn = 0.d0
      V_p = Vl - vn
      V_m = Vr - vn
      M_p = V_p / c
      M_m = V_m / c
      M = min(1.d0, sqrt(0.5d0 * (M_p ** 2 + M_m ** 2)))
      x = (1.d0 - M) ** 2
      g = -max(min(M_p,0.d0),-1.d0) * min(max(M_m,0.d0), 1.d0)
      V_bar = (rhol * abs(V_p) + rhor * abs(V_m)) / (rhol + rhor)
      V_bar_p = abs((1.d0 - g) * V_bar + g * V_p)
      V_bar_m = abs((1.d0 - g) * V_bar + g * V_m)
      dp = pr - pl
      mass = 0.5d0 * (rhol * (Vl + V_bar_p) + rhor * (Vr - V_bar_m) - x * dp / c)
      if (abs(M_p) < 1.d0) then
        beta_p = 0.25d0 * (2.d0 - M_p) * (M_p + 1.d0) ** 2
      else
        beta_p = 0.5d0 * (1.d0 + sign(1.d0,M_p))
      endif
      if (abs(M_m) < 1.d0) then
        beta_m = 0.25d0 * (2.d0 + M_m) * (M_m - 1.d0) ** 2
      else
        beta_m = 0.5d0 * (1.d0 + sign(1.d0,-M_m))
      endif
      Pressure = 0.5d0 * (pl + pr + (beta_p - beta_m) * (pl - pr) + (1.d0 - x) * (beta_p + beta_m - 1.d0) * (pl + pr))
      F(i,:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + Pressure * Normal(:)
    enddo
  end subroutine
end module calc_slau

