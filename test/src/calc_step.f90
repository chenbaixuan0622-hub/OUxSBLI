module calc_step
  implicit none
contains
  subroutine calc_step1(nx, dt, q, RHS, q2)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: dt, q(nx)
    real(8), intent(in)  :: RHS(nx-2)
    real(8), intent(out) :: q2(nx)
    integer i
    do i = 2, nx-1
      q2(i) = q(i) - dt * RHS(i-1)
    enddo
  end subroutine calc_step1


  subroutine calc_step2(nx, dt, q, RHS, q2)
    integer, intent(in)    :: nx
    real(8), intent(in)    :: dt, q(nx)
    real(8), intent(in)    :: RHS(nx-2)
    real(8), intent(inout) :: q2(nx)
    integer i
    do i = 2, nx-1
      q2(i) = 0.25d0 * (3.d0 * q(i) + q2(i) - dt * RHS(i-1))
    enddo
  end subroutine calc_step2
  

  subroutine calc_step3(nx, dt, q2, RHS, q)
    integer, intent(in)    :: nx
    real(8), intent(in)    :: dt, q2(nx)
    real(8), intent(in)    :: RHS(nx-2)
    real(8), intent(inout) :: q(nx)
    integer i
    do i = 2, nx-1
      q(i) = (q(i) + 2.d0 * q2(i) - 2.d0 * dt * RHS(i-1)) / 3.d0
    enddo
  end subroutine calc_step3
end module calc_step

