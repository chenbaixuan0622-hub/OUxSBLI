module calc_muscl
  implicit none
  interface minmod
    module procedure minmod2, minmod3
  end interface

  interface al
    module procedure al_3rd, al_4th
  end interface

  interface ar
    module procedure ar_3rd, ar_4th
  end interface
contains

!minmod!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function minmod2(x,y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = dsign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod2

  attributes(device) function minmod3(x,y,z) result(ans)
    real(8), intent(in), value :: x, y, z
    real(8) :: ans, sgn
    sgn = dsign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y, sgn * z), 0.d0)
  end function minmod3

!3rd-order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function al_3rd(a,eps,k,b,d_p,d_m) result(ans)
    real(8), intent(in), value :: a, eps, k, b, d_p, d_m
    real(8) :: ans, delta_p, delta_m
    delta_p = minmod(d_p, b * d_m)
    delta_m = minmod(d_m, b * d_p)
    ans = a + 0.25d0 * eps * ((1.d0 - k) * delta_m + (1.d0 + k) * delta_p)
  end function al_3rd

  attributes(device) function ar_3rd(a,eps,k,b,d_p,d_m) result(ans)
    real(8), intent(in), value :: a, eps, k, b, d_p, d_m
    real(8) :: ans, delta_p, delta_m
    delta_p = minmod(d_p, b * d_m)
    delta_m = minmod(d_m, b * d_p)
    ans = a - 0.25d0 * eps * ((1.d0 - k) * delta_p + (1.d0 + k) * delta_m)
  end function ar_3rd

!4th-order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function d33(sigma,d1,d2,d3) result(ans)
    real(8), intent(in), value :: sigma, d1, d2, d3 
    real(8) :: ans, da, db, dc
    da = minmod(d1, sigma * d2, sigma * d3)
    db = minmod(d2, sigma * d1, sigma * d3)
    dc = minmod(d3, sigma * d1, sigma * d2)
    ans = da - 2.d0 * db + dc
  end function d33

  attributes(device) function al_4th(a,eps,w,sigma,d1,d2,d3,d4) result(ans)
    real(8), intent(in), value :: a, eps, w, sigma, d1, d2, d3, d4
    real(8) :: ans, delta_p, delta_m, dp, dm
    dp = d3 - eps * d33(sigma, d2, d3, d4) / 6.d0
    dm = d2 - eps * d33(sigma, d1, d2, d3) / 6.d0
    delta_p = minmod(dp, w * dm)
    delta_m = minmod(dm, w * dp)
    ans = a + (delta_m + 2.d0 * delta_p) / 6.d0
  end function al_4th

  attributes(device) function ar_4th(a,eps,w,sigma,d1,d2,d3,d4) result(ans)
    real(8), intent(in), value :: a, eps, w, sigma, d1, d2, d3, d4
    real(8) :: ans, delta_p, delta_m, dp, dm
    dp = d3 - eps * d33(sigma, d2, d3, d4) / 6.d0
    dm = d2 - eps * d33(sigma, d1, d2, d3) / 6.d0
    delta_p = minmod(dp, w * dm)
    delta_m = minmod(dm, w * dp)
    ans = a - (delta_p + 2.d0 * delta_m) / 6.d0
  end function ar_4th

!MUSCL interpolation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) subroutine MUSCL(dim,eps_left,eps_right,k,b,Q1,Q2,d1,d2,d3,Ql,Qr)
    integer, intent(in), value :: dim
    real(8), intent(in), value :: eps_left, eps_right, k, b
    real(8), intent(in), dimension(dim), device :: Q1, Q2, d1, d2, d3
    real(8), intent(out), dimension(dim), device :: Ql, Qr
    integer i
    do i = 1, dim
      Ql(i) = al(Q1(i), eps_left, k, b, d2(i), d1(i))
      Qr(i) = ar(Q2(i), eps_right, k, b, d3(i), d2(i))
    enddo
  end subroutine MUSCL

  attributes(device) subroutine MUSCL_4th(dim,eps_left,eps_right,w,sigma,Q1,Q2,d1,d2,d3,d4,d5,Ql,Qr)
    integer, intent(in), value :: dim
    real(8), intent(in), value :: eps_left, eps_right, w, sigma
    real(8), intent(in), dimension(dim), device :: Q1, Q2, d1, d2, d3, d4, d5
    real(8), intent(out), dimension(dim), device :: Ql, Qr
    integer i
    do i = 1, dim
      Ql(i) = al(Q1(i), eps_left, w, sigma, d4(i), d3(i), d2(i), d1(i))
      Qr(i) = ar(Q2(i), eps_right, w, sigma, d5(i), d4(i), d3(i), d2(i))
    enddo
  end subroutine MUSCL_4th
end module calc_muscl

