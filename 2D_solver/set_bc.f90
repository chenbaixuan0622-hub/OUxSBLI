module set_bc
  implicit none
  ! only 2nd-order accuracy is available
contains
  subroutine set_tube_bc(id,nx,ny,gamma,Q)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(inout), device :: Q(nx,ny,4)
    real(8) p_bottom
    integer i, j, k
    !$acc kernels deviceptr(Q)
    !$acc loop independent
    do k = 1, 4
      !$acc loop independent
      do j = 1, ny
        ! inlet
        Q(1,j,k) = Q(2,j,k)
        ! outlet
        Q(nx,j,k) = Q(nx-1,j,k)
      enddo
    enddo

    !$acc loop independent
    do i = 1, nx
      ! bottom-wall
      Q(i,1,1) = Q(i,2,1)
      Q(i,1,2) = 0.d0
      Q(i,1,3) = 0.d0
      p_bottom = (gamma - 1.d0) * (Q(i,2,4) - 0.5d0 * (Q(i,2,2)**2 + Q(i,2,3)**2) / Q(i,2,1))
      Q(i,1,4) = p_bottom / (gamma - 1.d0)
      ! top-wall (symmetry, ny-1 is the symmetry axis)
      Q(i,ny,:) = Q(i,ny-2,:)
    enddo
    !$acc end kernels
  end subroutine

  subroutine wind_tunnel_with_a_step(nx,ny,gamma,Q)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(in), device :: Q(nx,ny,4)
    integer i, j, k, nxs, nys
    real(8), parameter :: rho0 = 1.4d0, u0 = 3.d0, p0 = 1.d0
    real(8) p_bottom, p_top, p_step_left, p_step_top
    nxs = int(nx * 0.2)
    nys = int(ny * 0.2)
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
    
    !$acc loop independent
    do i = 1, nx
      ! bottom wll
      Q(i,1,1) = Q(i,2,1)
      Q(i,1,2) = 0.d0
      Q(i,1,3) = 0.d0
      p_bottom = (gamma - 1.d0) * (Q(i,2,4) - 0.5d0 * (Q(i,2,2)**2 + Q(i,2,3)**2) / Q(i,2,1))
      Q(i,1,4) = p_bottom / (gamma - 1.d0)
      ! top wall
      Q(i,nx,1) = Q(i,nx-1,2)
      Q(i,nx,2) = 0.d0
      Q(i,nx,3) = 0.d0
      p_top = (gamma - 1.d0) * (Q(i,nx-1,4) - 0.5d0 * (Q(i,nx-1,2)**2 + Q(i,nx-1,3)**2) / Q(i,nx-1,1))
      Q(i,nx,4) = p_top / (gamma - 1.d0)
    enddo

    ! step
    !$acc loop independent
    do j = 1, nys
      ! left
      Q(nxs:nxs+1,j,1) = Q(nxs-1,j,1)
      Q(nxs:nxs+1,j,2) = 0.d0
      Q(nxs:nxs+1,j,3) = 0.d0
      p_step_left = (gamma - 1.d0) * (Q(nxs-1,j,4) - 0.5d0 * (Q(nxs-1,j,2)**2 + Q(nxs-1,j,3)**2) / Q(nxs-1,j,1))
      Q(nxs:nxs+1,j,4) = p_step_left / (gamma - 1.d0)
    enddo

    !$acc loop independent
    do i = nxs, nx
      ! top
      Q(i,nys-1:nys,1) = Q(i,nys+1,1)
      Q(i,nys-1:nys,2) = 0.d0
      Q(i,nys-1:nys,3) = 0.d0
      p_step_top = (gamma - 1.d0) * (Q(i,nys+1,4) - 0.5d0 * (Q(i,nys+1,2)**2 + Q(i,nys+1,3)**2) / Q(i,nys+1,1))
      Q(i,nys-1:nys,4) = p_step_top / (gamma - 1.d0)
    enddo

    !$acc loop independent
    do j = 1, nys-2
      !$acc loop independent
      do i = nxs+2, nx
        Q(i,j,1) = rho0
        Q(i,j,2) = 0.d0
        Q(i,j,3) = 0.d0
        Q(i,j,4) = p0 / (gamma - 1.d0)
      enddo
    enddo
    !$acc end kernels
  end subroutine wind_tunnel_with_a_step
end module set_bc

