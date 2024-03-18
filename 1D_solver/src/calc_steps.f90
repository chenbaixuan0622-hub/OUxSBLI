module calc_steps
  use mod_globals, only : accuracy, offset
  implicit none
contains
  subroutine calc_step1(nx,dtdx,E,Q,Q2)
    integer, intent(in), value :: nx
    real(8), intent(in), value :: dtdx
    real(8), intent(in) :: E(nx-accuracy+1,3)
    real(8), intent(in), dimension(nx,3) :: Q
    real(8), intent(out), dimension(nx,3) :: Q2
    integer i
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    do i = 1+offset, nx-offset
      Q2(i,:) = Q(i,:) - dtdx * (-E(i-offset,:) + E(i-offset+1,:))
    enddo
  end subroutine calc_step1
  
  subroutine calc_step2(nx,dtdx,E,Q,Q2,Q3)
    integer, intent(in), value :: nx
    real(8), intent(in), value :: dtdx
    real(8), intent(in) :: E(nx-accuracy+1,3)
    real(8), intent(in), dimension(nx,3) :: Q
    real(8), intent(in), dimension(nx,3) :: Q2
    real(8), intent(out), dimension(nx,3) :: Q3
    integer i
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    do i = 1+offset, nx-offset
      Q3(i,:) = 0.75d0 * Q(i,:) + 0.25d0 * Q2(i,:) - 0.25d0 * dtdx * (-E(i-offset,:) + E(i-offset+1,:))
    enddo
  end subroutine calc_step2

  subroutine calc_step3(nx,dtdx,E,Q3,Q)
    integer, intent(in), value :: nx
    real(8), intent(in), value :: dtdx
    real(8), intent(in) :: E(nx-accuracy+1,3)
    real(8), intent(in), dimension(nx,3) :: Q3
    real(8), intent(inout), dimension(nx,3) :: Q
    integer i
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    do i = 1+offset, nx-offset
      Q(i,:) = (Q(i,:) + 2.d0 * Q3(i,:) - 2.d0 * dtdx * (-E(i-offset,:) + E(i-offset+1,:))) / 3.d0
    enddo
  end subroutine calc_step3
end module calc_steps

