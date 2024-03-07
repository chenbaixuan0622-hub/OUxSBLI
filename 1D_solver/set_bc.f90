module set_bc
  implicit none
  ! only 2nd-order accuracy is available
contains
  subroutine set_tube_bc(id,nx,gamma,Q)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx
    real(8), intent(in), value :: gamma
    real(8), intent(inout) :: Q(nx,3)
    ! inlet
    Q(1,:) = Q(2,:)
    ! outlet
    Q(nx,:) = Q(nx-1,:)
  end subroutine
end module set_bc

