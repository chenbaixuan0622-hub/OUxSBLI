module set_init
  implicit none
contains
  subroutine shock_tube(nx,gamma,rhol,rhor,pl,pr,Q)
    integer, intent(in) :: nx
    real(8), intent(in) :: gamma, rhol, rhor, pl, pr
    real(8), intent(out), dimension(nx,3) :: Q
    ! left half
    Q(1:int(0.5*nx),1) = rhol
    Q(1:int(0.5*nx),3) = pl / (gamma - 1.d0)
    ! right half
    Q(int(0.5*nx):nx,1) = rhor
    Q(int(0.5*nx):nx,3) = pr / (gamma - 1.d0)
    ! both (velocity = 0)
    Q(:,2) = 0.d0
  end subroutine shock_tube
end module set_init

