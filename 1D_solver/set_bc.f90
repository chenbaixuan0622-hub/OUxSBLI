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
end module set_bc

