module set
  use mod_globals, only : rho0, p0, u0
  implicit none
contains
  subroutine set_init(nx,ny,gamma,Q)
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: gamma
    real(8), intent(out), dimension(nx,ny,4) :: Q
    ! wind tunnel
    Q(:,:,1) = rho0
    Q(:,:,2) = 0.d0
    Q(:,:,3) = 0.d0
    Q(:,:,4) = p0 / (gamma - 1.d0)
  end subroutine set_init
  
  subroutine set_bc(nx,ny,gamma,Q)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k, nxs, nys
    real(8) p_bottom, p_top, p_step_left, p_step_top
    nxs = int(0.2 * nx)
    nys = int(0.2 * ny)
    
    !$acc kernels deviceptr(Q)
    !$acc loop independent
    do j = 1, ny
      ! inlet
      Q(1,j,1) = rho0
      Q(1,j,2) = rho0 * u0
      Q(1,j,3) = 0.d0
      Q(1,j,4) = p0 / (gamma - 1.d0)
      ! outlet
      Q(nx,j,:) = Q(nx-1,j,:)
    enddo
    !$acc end kernels

    !$acc kernels deviceptr(Q)
    !$acc loop independent private(p_bottom,p_top)
    do i = 2, nx-1
      ! bottom wll
      Q(i,1,1) = Q(i,2,1)
      Q(i,1,2) = 0.d0
      Q(i,1,3) = 0.d0
      p_bottom = (gamma - 1.d0) * (Q(i,2,4) - 0.5d0 * (Q(i,2,2)**2 + Q(i,2,3)**2) / Q(i,2,1))
      Q(i,1,4) = p_bottom / (gamma - 1.d0)
      ! top wall
      Q(i,ny,1) = Q(i,ny-1,2)
      Q(i,ny,2) = 0.d0
      Q(i,ny,3) = 0.d0
      p_top = (gamma - 1.d0) * (Q(i,ny-1,4) - 0.5d0 * (Q(i,ny-1,2)**2 + Q(i,ny-1,3)**2) / Q(i,ny-1,1))
      Q(i,ny,4) = p_top / (gamma - 1.d0)
    enddo
    !$acc end kernels

    ! step
    !$acc kernels deviceptr(Q)
    !$acc loop independent private(p_step_left)
    do j = 2, nys
      ! left
      Q(nxs:nxs+1,j,1) = Q(nxs-1,j,1)
      Q(nxs:nxs+1,j,2) = 0.d0
      Q(nxs:nxs+1,j,3) = 0.d0
      p_step_left = (gamma - 1.d0) * (Q(nxs-1,j,4) - 0.5d0 * (Q(nxs-1,j,2)**2 + Q(nxs-1,j,3)**2) / Q(nxs-1,j,1))
      Q(nxs:nxs+1,j,4) = p_step_left / (gamma - 1.d0)
    enddo
    !$acc end kernels

    !$acc kernels deviceptr(Q)
    !$acc loop independent private(p_step_top)
    do i = nxs+2, nx
      ! top
      Q(i,nys-1:nys,1) = Q(i,nys+1,1)
      Q(i,nys-1:nys,2) = 0.d0
      Q(i,nys-1:nys,3) = 0.d0
      p_step_top = (gamma - 1.d0) * (Q(i,nys+1,4) - 0.5d0 * (Q(i,nys+1,2)**2 + Q(i,nys+1,3)**2) / Q(i,nys+1,1))
      Q(i,nys-1:nys,4) = p_step_top / (gamma - 1.d0)
    enddo
    !$acc end kernels
  end subroutine set_bc
end module set

