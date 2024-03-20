module calc_qlr
  use mod_globals, only : dim => dimension, k, b
  use calc_MUSCL
  implicit none
contains
  attributes(device) subroutine Qlr_left(Q1,Q2,Q3,Ql,Qr)
    real(8), intent(in), dimension(dim+2) :: Q1, Q2, Q3 
    real(8), intent(out), dimension(dim+2) :: Ql, Qr
    real(8), dimension(dim+2) :: d1, d2, d3
    d1(:) = 0.d0
    d2(:) = -Q1(:) + Q2(:)
    d3(:) = -Q2(:) + Q3(:)
    call MUSCL(dim+2,k,b,Q1,Q2,d1,d2,d3,Ql,Qr)
  end subroutine Qlr_left
  
  attributes(device) subroutine Qlr_mid(Q1,Q2,Q3,Q4,Ql,Qr)
    real(8), intent(in), dimension(dim+2) :: Q1, Q2, Q3, Q4 
    real(8), intent(out), dimension(dim+2) :: Ql, Qr
    real(8), dimension(dim+2) :: d1, d2, d3
    d1(:) = -Q1(:) + Q2(:)
    d2(:) = -Q2(:) + Q3(:)
    d3(:) = -Q3(:) + Q4(:)
    call MUSCL(dim+2,k,b,Q2,Q3,d1,d2,d3,Ql,Qr)
  end subroutine Qlr_mid
  
  attributes(device) subroutine Qlr_right(Q1,Q2,Q3,Ql,Qr)
    real(8), intent(in), dimension(dim+2) :: Q1, Q2, Q3 
    real(8), intent(out), dimension(dim+2) :: Ql, Qr
    real(8), dimension(dim+2) :: d1, d2, d3
    d1(:) = -Q1(:) + Q2(:)
    d2(:) = -Q2(:) + Q3(:)
    d3(:) = 0.d0
    call MUSCL(dim+2,k,b,Q2,Q3,d1,d2,d3,Ql,Qr)
  end subroutine Qlr_right
end module calc_qlr

