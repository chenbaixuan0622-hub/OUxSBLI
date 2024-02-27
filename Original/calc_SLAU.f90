module calc_SLAU
  implicit none
contains
  function calc_c(gamma,pl,pr,rhol,rhor) result(ans)
    real(8), intent(in), value :: gamma, pl, pr, rhol, rhor
    real(8) cl, cr
    real(8) ans
    cl = sqrt(gamma * pl / rhol)
    cr = sqrt(gamma * pr / rhor)
    ans = 0.5d0 * (cl + cr)
  end function

  function calc_x(M_plus,M_minus) result(ans)
    real(8), intent(in), value :: M_plus, M_minus 
    real(8) M_bar
    real(8) ans
    M_bar = min(1.d0, sqrt(0.5d0 * (M_plus**2 + M_minus**2)))
    ans = (1.d0 - M_bar) ** 2
  end function
  
  function calc_m(rhol,rhor,pl,pr,c,M_plus,M_minus,x,Vl,Vr,Vn) result(ans)
    real(8), intent(in), value :: rhol, rhor, pl, pr, M_plus, M_minus, c, x, Vl, Vr, Vn
    real(8) V_plus, V_minus, V_bar, V_bar_plus, V_bar_minus, g
    real(8) ans
    g = -max(min(M_plus,0.d0), -1.d0) * min(max(M_minus, 0.d0), 1.d0)
    V_plus = abs(Vl - Vn)
    V_minus = abs(Vr - Vn)
    V_bar = (rhol * V_plus + rhor * V_bar_minus) / (rhol + rhor)
    V_bar_plus = (1.d0 - g) * V_bar + g * V_plus
    V_bar_minus = (1.d0 - g) * V_bar + g * V_minus
    ans = 0.5d0*(rhol*(Vl + abs(V_bar_plus)) + rhor*(Vr - abs(V_bar_minus)) - x*(pr - pl)/c)
  end function
  
  function calc_P(pl,pr,M_plus,M_minus,x) result(ans)
    real(8), intent(in), value :: pl, pr, M_plus, M_minus, x
    real(8) p, b_plus, b_minus
    real(8) ans
    p = 0.5d0 * (pl + pr)
    if (abs(M_plus) < 1.d0) then 
      b_plus = 0.25d0 * (2.d0 - M_plus) * (M_plus + 1.d0) ** 2
    else
      b_plus = 0.5d0 * (1.d0 + sign(1.d0, M_plus))
    endif
    if (abs(M_minus) < 1.d0) then
      b_minus = 0.25d0 * (2.d0 + M_minus) * (M_minus - 1.d0) ** 2
    else
      b_minus = 0.5d0 * (1.d0 + sign(1.d0, -M_minus))
    endif
    ans = p + 0.5d0 * (b_plus - b_minus) * (pl - pr) &
    & + (1.d0 - x) * (b_plus + b_minus - 1.d0) * p
  end function
  
  subroutine calc_Flux(nx, ny, nz, gamma, Ql, Qr, N, E)
    integer, intent(in), value :: nx, ny, nz
    integer, intent(in) :: N(3)
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz,5) :: Ql, Qr
    real(8), intent(out), dimension(nx,ny,nz,5) :: E
    integer i, j, k
    real(8) rhol, rhor, ul, ur, vl, vr, wl, wr, el, er, pl, pr, c, Vn, Vnl, Vnr, M_plus, M_minus, x, m, ml, mr, p
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          !Q(rho, u, v, w, p)
          rhol = Ql(i,j,k,1)
          rhor = Qr(i,j,k,1)
          ul = Ql(i,j,k,2)
          ur = Qr(i,j,k,2)
          vl = Ql(i,j,k,3)
          vr = Qr(i,j,k,3)
          wl = Ql(i,j,k,4)
          wr = Qr(i,j,k,4)
          pl = Ql(i,j,k,5)
          pr = Qr(i,j,k,5)
          el = pl / (gamma - 1.d0) + 0.5d0 * rhol * (ul**2 + vl**2 + wl**2)
          er = pr / (gamma - 1.d0) + 0.5d0 * rhor * (ur**2 + vr**2 + wr**2)
          c = calc_c(gamma,pl,pr,rhol,rhor)
          ! What is Vn ?
          Vn = 0.d0
          Vnl = ul * N(1) + vl * N(2) + wl * N(3)
          Vnr = ur * N(1) + vr * N(2) + wr * N(3)
          M_plus = (Vnl - Vn) / c
          M_minus = (Vnr - Vn) / c
          x = calc_x(M_plus,M_minus)
          m = calc_m(rhol,rhor,pl,pr,c,M_plus,M_minus,x,Vl,Vr,Vn)
          ml = 0.5d0 * (m + abs(m))
          mr = 0.5d0 * (m - abs(m))
          p = calc_P(pl,pr,M_plus,M_minus,x)
          E(i,j,k,1) = ml + mr
          E(i,j,k,2) = ml * ul + mr * ur + N(1) * p
          E(i,j,k,3) = ml * vl + mr * vr + N(2) * p
          E(i,j,k,4) = ml * wl + mr * wr + N(3) * p
          E(i,j,k,5) = ml * (el + pl) / rhol + mr * (er + pr) / rhor
        enddo
      enddo
    enddo
  end subroutine calc_Flux
end module calc_SLAU
  