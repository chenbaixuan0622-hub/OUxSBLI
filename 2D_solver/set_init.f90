module set_init
  implicit none
contains
  subroutine shock_tube(nx,ny,gamma,rhol,rhor,pl,pr,Q)
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: gamma, rhol, rhor, pl, pr
    real(8), intent(out), dimension(nx,ny,4) :: Q
    ! left half
    Q(1:int(0.5*nx),:,1) = rhol
    Q(1:int(0.5*nx),:,4) = pl / (gamma - 1.d0)
    ! right half
    Q(int(0.5*nx):nx,:,1) = rhor
    Q(int(0.5*nx):nx,:,4) = pr / (gamma - 1.d0)
    ! both (velocity = 0)
    Q(:,:,2) = 0.d0
    Q(:,:,3) = 0.d0
  end subroutine shock_tube

  subroutine wind_tunnel_with_a_step(nx,ny,gamma,Q)
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: gamma
    real(8), parameter :: rho0 = 1.4d0, p0 = 1.d0
    real(8), intent(out), dimension(nx,ny,4) :: Q
    ! wind tunnel
    Q(:,:,1) = rho0
    Q(:,:,2) = 0.d0
    Q(:,:,3) = 0.d0
    Q(:,:,4) = p0 / (gamma - 1.d0)
  end subroutine
end module set_init

