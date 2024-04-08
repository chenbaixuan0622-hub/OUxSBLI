module calc_qlr
  use mod_globals, only : dim => dimension, k, b, omega, sigma, eps
  use calc_MUSCL
  implicit none
  interface Qlr_mid
    module procedure Qlr_4points, Qlr_6points
  end interface Qlr_mid

contains
  attributes(device) subroutine Qlr_left(Q1,Q2,Q3,Ql,Qr)
    real(8), intent(in), dimension(dim+2) :: Q1, Q2, Q3 
    real(8), intent(out), dimension(dim+2) :: Ql, Qr
    real(8), dimension(dim+2) :: d1, d2, d3
    real(8) :: eps_left = 0.d0
    d1(:) = 0.d0
    d2(:) = -Q1(:) + Q2(:)
    d3(:) = -Q2(:) + Q3(:)
    call MUSCL(dim+2,eps_left,eps,k,b,Q1,Q2,d1,d2,d3,Ql,Qr)
  end subroutine Qlr_left
  
!3rd-order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) subroutine Qlr_4points(Q1,Q2,Q3,Q4,Ql,Qr)
    real(8), intent(in), dimension(dim+2) :: Q1, Q2, Q3, Q4 
    real(8), intent(out), dimension(dim+2) :: Ql, Qr
    real(8), dimension(dim+2) :: d1, d2, d3
    d1(:) = -Q1(:) + Q2(:)
    d2(:) = -Q2(:) + Q3(:)
    d3(:) = -Q3(:) + Q4(:)
    call MUSCL(dim+2,eps,eps,k,b,Q2,Q3,d1,d2,d3,Ql,Qr)
  end subroutine Qlr_4points
  
!4th-order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) subroutine Qlr_6points(eps_left,eps_right,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    real(8), intent(in) :: eps_left, eps_right
    real(8), intent(in), dimension(dim+2) :: Q1, Q2, Q3, Q4, Q5, Q6
    real(8), intent(out), dimension(dim+2) :: Ql, Qr
    real(8), dimension(dim+2) :: d1, d2, d3, d4, d5
    d1(:) = -Q1(:) + Q2(:)
    d2(:) = -Q2(:) + Q3(:)
    d3(:) = -Q3(:) + Q4(:)
    d4(:) = -Q4(:) + Q5(:)
    d5(:) = -Q5(:) + Q6(:)
    call MUSCL_4th(dim+2,eps_left,eps_right,omega,sigma,Q3,Q4,d1,d2,d3,d4,d5,Ql,Qr)
  end subroutine Qlr_6points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) subroutine Qlr_right(Q1,Q2,Q3,Ql,Qr)
    real(8), intent(in), dimension(dim+2) :: Q1, Q2, Q3 
    real(8), intent(out), dimension(dim+2) :: Ql, Qr
    real(8), dimension(dim+2) :: d1, d2, d3
    real(8) :: eps_right = 0.d0
    d1(:) = -Q1(:) + Q2(:)
    d2(:) = -Q2(:) + Q3(:)
    d3(:) = 0.d0
    call MUSCL(dim+2,eps,eps_right,k,b,Q2,Q3,d1,d2,d3,Ql,Qr)
  end subroutine Qlr_right
end module calc_qlr

