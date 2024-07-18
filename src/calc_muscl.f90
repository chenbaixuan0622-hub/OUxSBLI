module calc_muscl
  implicit none
  interface minmod
    module procedure minmod2, minmod3
  end interface

  interface MUSCL
    module procedure MUSCL3rdnonTVD, MUSCL3rdMinmod, MUSCL4th
  end interface
contains
  attributes(device) function minmod2(x,y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod2

  attributes(device) function minmod3(x,y,z) result(ans)
    real(8), intent(in), value :: x, y, z
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y, sgn * z), 0.d0)
  end function minmod3

  attributes(device) function d33(sigma,d1,d2,d3) result(ans)
    real(8), intent(in), value :: sigma, d1, d2, d3 
    real(8) :: ans, da, db, dc
    da = minmod(d1, sigma * d2, sigma * d3)
    db = minmod(d2, sigma * d1, sigma * d3)
    dc = minmod(d3, sigma * d1, sigma * d2)
    ans = da - 2.d0 * db + dc
  end function d33

  attributes(device) function MUSCL3rdnonTVD(id_tvd,eps,k,a2,a3,d) result(alr)
    integer(kind=2), intent(in), value :: id_tvd
    real(8), intent(in), value         :: eps, k, a2, a3
    real(8), intent(in), device        :: d(3)
    real(8) :: b, alr(2)
    b      = (3.d0 - k) / (1.d0 -k)
    alr(1) = a2 + 0.25d0 * eps * ((1.d0 - k) * d(1) + (1.d0 + K) * d(2))
    alr(2) = a3 - 0.25d0 * eps * ((1.d0 - k) * d(3) + (1.d0 + K) * d(2))
  end function MUSCL3rdnonTVD

  attributes(device) function MUSCL3rdMinmod(id_tvd,eps,k,a2,a3,d) result(alr)
    integer(kind=4), intent(in), value :: id_tvd
    real(8), intent(in), value         :: eps, k, a2, a3
    real(8), intent(in), device        :: d(3)
    real(8) :: b, dt1, dt2, dt3, dt4, alr(2)
    b   = (3.d0 - k) / (1.d0 -k)
    dt1 = minmod(d(1), b * d(2))
    dt2 = minmod(d(2), b * d(1))
    dt3 = minmod(d(3), b * d(2))
    dt4 = minmod(d(2), b * d(3))
    alr(1) = a2 + 0.25d0 * eps * ((1.d0 - k) * dt1 + (1.d0 + K) * dt2)
    alr(2) = a3 - 0.25d0 * eps * ((1.d0 - k) * dt3 + (1.d0 + K) * dt4)
  end function MUSCL3rdMinmod

  attributes(device) function MUSCL4th(id_tvd,eps,k,a2,a3,d) result(alr)
    integer(kind=8), intent(in), value :: id_tvd
    real(8), intent(in), value         :: eps, k, a2, a3
    real(8), intent(in), device        :: d(5)
    real(8) delta1, delta2, delta3, dpl, dml, dpr, dmr, alr(2)
    real(8) :: sigma = 2.d0, w = 4.d0
    delta1 = d(2) - eps * d33(sigma, d(1), d(2), d(3)) / 6.d0
    delta2 = d(3) - eps * d33(sigma, d(2), d(3), d(4)) / 6.d0
    delta3 = d(4) - eps * d33(sigma, d(3), d(4), d(5)) / 6.d0
    dpl    = minmod(delta1, w * delta2)
    dml    = minmod(delta2, w * delta1)
    alr(1) = a2 + (dml + 2.d0 * dpl) / 6.d0
    dpr    = minmod(delta3, w * delta2)
    dmr    = minmod(delta2, w * delta3)
    alr(2) = a3 - (dpr + 2.d0 * dmr) / 6.d0
  end function MUSCL4th

  attributes(device) function delta4(eps,k,a) result(alr)
    use mod_globals, only : id_tvd
    real(8), intent(in), value  :: eps, k
    real(8), intent(in), device :: a(4)
    real(8) :: alr(2), d(3)
    integer i
    do i = 1, 3
      d(i) = -a(i) + a(i+1)
    enddo
    alr  = MUSCL(id_tvd,eps,k,a(2),a(3),d)
  end function delta4

  attributes(device) function delta6(eps,k,a) result(alr)
    use mod_globals, only : id_tvd
    real(8), intent(in), value  :: eps, k
    real(8), intent(in), device :: a(6)
    real(8) :: alr(2), d(5)
    integer i
    do i = 1, 5
      d(i) = -a(i) + a(i+1)
    enddo
    alr  = MUSCL(id_tvd,eps,k,a(3),a(4),d)
  end function delta6

  attributes(device) subroutine calc_4points(eps1,eps2,eps3,k,rho,p,V,rho2,p2,V2)
    real(8), intent(in), value   :: eps1, eps2, eps3, k
    real(8), intent(in), device  :: rho(4),  p(4),  V(4,3)
    real(8), intent(out), device :: rho2(2), p2(2), V2(2,3)
    real(8) Vtemp(4)
    rho2(:) = delta4(eps1,k,rho)
    p2(:)   = delta4(eps2,k,p)
    Vtemp   = V(:,1)
    V2(:,1) = delta4(eps3,k,Vtemp)
    Vtemp   = V(:,2)
    V2(:,2) = delta4(eps3,k,Vtemp)
    Vtemp   = V(:,3)
    V2(:,3) = delta4(eps3,k,Vtemp)
  end subroutine calc_4points

  attributes(device) subroutine calc_6points(eps1,eps2,eps3,k,rho,p,V,rho2,p2,V2)
    real(8), intent(in), value   :: eps1, eps2, eps3, k
    real(8), intent(in), device  :: rho(6),  p(6),  V(6,3)
    real(8), intent(out), device :: rho2(2), p2(2), V2(2,3)
    real(8) Vtemp(6)
    rho2(:) = delta6(eps1,k,rho)
    p2(:)   = delta6(eps2,k,p)
    Vtemp   = V(:,1)
    V2(:,1) = delta6(eps3,k,Vtemp)
    Vtemp   = V(:,2)
    V2(:,2) = delta6(eps3,k,Vtemp)
    Vtemp   = V(:,3)
    V2(:,3) = delta6(eps3,k,Vtemp)
  end subroutine calc_6points
end module calc_muscl

