module calc_steps
  use cudafor
  use mod_constant, only : one_third
  implicit none
contains
  subroutine calc_step1(nx, Q, R, Q2)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: Q(3,nx)
    real(8), intent(in), device  :: R(3,nx-2)
    real(8), intent(out), device :: Q2(3,nx)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 2, nx-1
      do j = 1, 3
        Q2(j,i) = Q(j,i) - R(j,i-1)
    enddo;enddo
  end subroutine calc_step1


  subroutine calc_step2(nx, Q, R, Q2)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q(3,nx)
    real(8), intent(in), device    :: R(3,nx-2)
    real(8), intent(inout), device :: Q2(3,nx)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 2, nx-1
      do j = 1, 3
        Q2(j,i) = 0.25d0 * (3.d0 * Q(j,i) + Q2(j,i) - R(j,i-1))
    enddo;enddo
  end subroutine calc_step2
  

  subroutine calc_step3(nx, Q2, R, Q)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q2(3,nx)
    real(8), intent(in), device    :: R(3,nx-2)
    real(8), intent(inout), device :: Q(3,nx)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 2, nx-1
      do j = 1, 3
        Q(j,i) = (Q(j,i) + 2.d0 * Q2(j,i) - 2.d0 * R(j,i-1)) * one_third
    enddo;enddo
  end subroutine calc_step3
end module calc_steps

