module calc_qlr
  use calc_MUSCL
  implicit none
  interface Qlr
    module procedure Qlr1, Qlr2, Qlr3
  end interface
contains
  subroutine common(k,b,Q1,Q2,d1,d2,d3,rhol,rhor,Vl,Vr,pl,pr)
    real(8), intent(in) :: k, b
    real(8), intent(in), dimension(3) :: Q1, Q2, d1, d2, d3
    real(8), intent(out) :: rhol, rhor, Vl, Vr, pl, pr
    real(8), dimension(3) :: Ql, Qr
    call MUSCL(3,k,b,Q1(:),Q2(:),d1(:),d2(:),d3(:),Ql(:),Qr(:))
    rhol = Ql(1)
    rhor = Qr(1)
    Vl = Ql(2)
    Vr = Qr(2)
    pl = Ql(3)
    pr = Qr(3)
  end subroutine common

  subroutine Qlr1(k,b,zero,Q1,Q2,Q3,rhol,rhor,Vl,Vr,pl,pr)
    real(8), intent(in) :: k, b, zero
    real(8), intent(in), dimension(3) :: Q1, Q2, Q3
    real(8), intent(out) :: rhol, rhor, Vl, Vr, pl, pr
    real(8), dimension(3) :: d1, d2, d3
    d1(:) = 0.d0
    d2(:) = -Q1(:) + Q2(:)
    d3(:) = -Q2(:) + Q3(:)
    call common(k,b,Q1,Q2,d1,d2,d3,rhol,rhor,Vl,Vr,pl,pr)
  end subroutine Qlr1

  subroutine Qlr2(k,b,Q1,Q2,Q3,Q4,rhol,rhor,Vl,Vr,pl,pr)
    real(8), intent(in) :: k, b
    real(8), intent(in), dimension(3) :: Q1, Q2, Q3, Q4
    real(8), intent(out) :: rhol, rhor, Vl, Vr, pl, pr
    real(8), dimension(3) :: d1, d2, d3
    d1(:) = -Q1(:) + Q2(:)
    d2(:) = -Q2(:) + Q3(:)
    d3(:) = -Q3(:) + Q4(:)
    call common(k,b,Q2,Q3,d1,d2,d3,rhol,rhor,Vl,Vr,pl,pr)
  end subroutine Qlr2

  subroutine Qlr3(k,b,Q1,Q2,Q3,zero,rhol,rhor,Vl,Vr,pl,pr)
    real(8), intent(in) :: k, b, zero
    real(8), intent(in), dimension(3) :: Q1, Q2, Q3
    real(8), intent(out) :: rhol, rhor, Vl, Vr, pl, pr
    real(8), dimension(3) :: d1, d2, d3
    d1(:) = -Q1(:) + Q2(:)
    d2(:) = -Q2(:) + Q3(:)
    d3(:) = 0.d0
    call common(k,b,Q2,Q3,d1,d2,d3,rhol,rhor,Vl,Vr,pl,pr)
  end subroutine Qlr3
end module calc_qlr

