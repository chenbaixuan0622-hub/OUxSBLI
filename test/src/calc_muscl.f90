module calc_muscl
  implicit none
  interface minmod
    module procedure minmod2, minmod3
  end interface
contains
  function minmod2(x, y) result(ans)
    real(8), intent(in) :: x, y
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod2


  function minmod3(x, y, z) result(ans)
    real(8), intent(in) :: x, y, z
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y, sgn * z), 0.d0)
  end function minmod3

  
  function d33(d1,d2,d3) result(ans)
    real(8), intent(in) :: d1, d2, d3 
    real(8) :: ans, da, db, dc
    da = minmod(d1, 2.d0 * d2, 2.d0 * d3)
    db = minmod(d2, 2.d0 * d1, 2.d0 * d3)
    dc = minmod(d3, 2.d0 * d1, 2.d0 * d2)
    ans = da - 2.d0 * db + dc
  end function d33


  function MUSCL_left(a) result(al)
    real(8), intent(in) :: a(3)
    real(8) :: b, d(2), d1, d2, al, k = 1.d0 / 3.d0
    b  = (3.d0 - k) / (1.d0 - k)
    d  = -a(1:2) + a(2:3)
    d1 = minmod(d(1), b * d(2))
    d2 = minmod(d(2), b * d(1))
    al = a(2) + 0.25d0 * ((1.d0 - k) * d1 + (1.d0 + k) * d2)
  end function MUSCL_left


  function MUSCL_right(a) result(ar)
    real(8), intent(in) :: a(3)
    real(8) :: b, d(2), d3, d4, ar, k = 1.d0 / 3.d0
    b  = (3.d0 - k) / (1.d0 - k)
    d  = -a(1:2) + a(2:3)
    d3 = minmod(d(2), b * d(1))
    d4 = minmod(d(1), b * d(2))
    ar = a(2) - 0.25d0 * ((1.d0 - k) * d3 + (1.d0 + k) * d4)
  end function MUSCL_right


  function MUSCL3rd(a) result(alr)
    real(8), intent(in) :: a(4)
    real(8) :: b, d(3), d1, d2, d3, d4, alr(2), k = 1.d0 / 3.d0
    b  = (3.d0 - k) / (1.d0 -k)
    d  = -a(1:3) + a(2:4)
    d1 = minmod(d(1), b * d(2))
    d2 = minmod(d(2), b * d(1))
    d3 = minmod(d(3), b * d(2))
    d4 = minmod(d(2), b * d(3))
    alr(1) = a(2) + 0.25d0 * ((1.d0 - k) * d1 + (1.d0 + K) * d2)
    alr(2) = a(3) - 0.25d0 * ((1.d0 - k) * d3 + (1.d0 + K) * d4)
  end function MUSCL3rd
  

  function MUSCL4th(a) result(alr)
    real(8), intent(in) :: a(6)
    real(8) d(5), delta1, delta2, delta3, dl, dr, alr(2)
    d      = -a(1:5) + a(2:6)
    delta1 = d(2) - d33(d(1), d(2), d(3)) / 6.d0
    delta2 = d(3) - d33(d(2), d(3), d(4)) / 6.d0
    delta3 = d(4) - d33(d(3), d(4), d(5)) / 6.d0
    dl     = minmod(delta1, 4.d0 * delta2)
    dr     = minmod(delta2, 4.d0 * delta1)
    alr(1) = a(3) + (dl + 2.d0 * dr) / 6.d0
    dl     = minmod(delta2, 4.d0 * delta3)
    dr     = minmod(delta3, 4.d0 * delta2)
    alr(2) = a(4) - (dr + 2.d0 * dl) / 6.d0
  end function MUSCL4th
    

  subroutine calc_rho_muscl(nx, rho, rho_m)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: rho(nx)
    real(8), intent(out) :: rho_m(nx-1)
    real(8) rho_lr(2)
    integer i
    do i = 1, nx-1
      if (3 <= i .and. i <= nx-3) then
        rho_lr = MUSCL4th(rho(i-2:i+3))
      elseif (2 <= i .and. i <= nx-2) then
        rho_lr = MUSCL3rd(rho(i-1:i+2))
      elseif (i == 1) then
        rho_lr(1) = 0.5d0 * (rho(1) + rho(2))
        rho_lr(2) = MUSCL_right(rho(i:i+2))
      else
        rho_lr(1) = MUSCL_left(rho(i-1:i+1))
        rho_lr(2) = 0.5d0 * (rho(nx-1) + rho(nx))
      endif
      rho_m(i) = 0.5d0 * (rho_lr(1) + rho_lr(2))
    enddo
  end subroutine calc_rho_muscl
end module calc_muscl

