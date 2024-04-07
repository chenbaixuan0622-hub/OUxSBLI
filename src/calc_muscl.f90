module calc_muscl
  implicit none
contains
  attributes(device) function minmod(x,y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = dsign(1.d0, x)
    ans = sgn * dmax(dmin(abs(x), sgn * y), 0.d0)
  end function minmod

  attributes(device) function al(a,k,b,d_p,d_m) result(ans)
    real(8), intent(in), value :: a, k, b, d_p, d_m
    real(8) :: ans, delta_p, delta_m
    delta_p = minmod(d_p, b * d_m)
    delta_m = minmod(d_m, b * d_p)
    ans = a + 0.25d0 * ((1.d0 - k) * delta_m + (1.d0 + k) * delta_p)
  end function al

  attributes(device) function ar(a,k,b,d_p,d_m) result(ans)
    real(8), intent(in), value :: a, k, b, d_p, d_m
    real(8) :: ans, delta_p, delta_m
    delta_p = minmod(d_p, b * d_m)
    delta_m = minmod(d_m, b * d_p)
    ans = a - 0.25d0 * ((1.d0 - k) * delta_p + (1.d0 + k) * delta_m)
  end function ar

  attributes(device) subroutine MUSCL(dim,k,b,Q1,Q2,d1,d2,d3,Ql,Qr)
    integer, intent(in), value :: dim
    real(8), intent(in), value :: k, b
    real(8), intent(in), dimension(dim), device :: Q1, Q2, d1, d2, d3
    real(8), intent(out), dimension(dim), device :: Ql, Qr
    integer i
    do i = 1, dim
      Ql(i) = al(Q1(i), k, b, d2(i), d1(i))
      Qr(i) = ar(Q2(i), k, b, d3(i), d2(i))
    enddo
  end subroutine MUSCL
end module calc_muscl

