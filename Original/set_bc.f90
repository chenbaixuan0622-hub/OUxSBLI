module set_bc
  implicit none
contains
  subroutine set_cyclic_bc(na,nx,ny,nz,Q)
  integer, intent(in) :: na, nx, ny, nz
  real(8), intent(inout) :: Q(nx,ny,nz,5)
  integer i, j, k, offset, offset_1
  ! 2nd-order accuracy
  ! 1 offset
  ! 4th-order accuracy
  ! 2 offset
  offset = na/2           ! 2, 1
  offset_1 = offset - 1   ! 1, 0

  do k = 1+offset, nz-offset
    do j = 1+offset, ny-offset
      ! 2nd
      !Q(1:1,j,k,:) = Q(nx-1:nx-1,j,k,:)
      !Q(nx:nx,j,k,:) = Q(2:2,j,k,:)
      ! 4th
      !Q(1:2,j,k,:) = Q(nx-3:nx-2,j,k,:)
      !Q(nx-1:nx,j,k,:) = Q(3:4,j,k,:)
      ! general
      Q(1:offset,j,k,:) = Q(nx-offset-offset_1:nx-offset,j,k,:)
      Q(nx-offset_1:nx,j,k,:) = Q(1+offset:1+offset+offset_1,j,k,:)
    enddo
  enddo

  do k = 1+offset, nz-offset
    do i = 1+offset, nx-offset
      ! 2nd
      !Q(i,1:1,k,:) = Q(i,ny-1:ny-1,k,:)
      !Q(i,ny:ny,k,:) = Q(i,2:2,k,:)
      ! 4th
      !Q(i,1:2,k,:) = Q(i,ny-3:ny-2,k,:)
      !Q(i,ny-1:ny,k,:) = Q(i,3:4,k,:)
      ! general
      Q(i,1:offset,k,:) = Q(i,ny-offset-offset_1:ny-offset,k,:)
      Q(i,ny-offset_1:ny,k,:) = Q(i,1+offset:1+offset+offset_1,k,:)
    enddo
  enddo

  do k = 1+offset, nz-offset
    ! 2nd
    !Q(1:1,1:1,k,:) = Q(nx-1:nx-1,ny-1:ny-1,k,:)
    !Q(nx:nx,1:1,k,:) = Q(2:2,ny-1:ny-1,k,:)
    !Q(1:1,ny:ny,k,:) = Q(nx-1:nx-1,2:2,k,:)
    !Q(nx:nx,ny:ny,k,:) = Q(2:2,2:2,k,:)
    ! 4th
    !Q(1:2,1:2,k,:) = Q(nx-3:nx-2,ny-3:ny-2,k,:)
    !Q(nx-1:nx,1:2,k,:) = Q(3:4,ny-3:ny-2,k,:)
    !Q(1:2,ny-1:ny,k,:) = Q(nx-3:nx-2,3:4,k,:)
    !Q(nx-1:nx,ny-1:ny,k,:) = Q(3:4,3:4,k,:)
    Q(1:offset,1:offset,k,:) = Q(nx-offset-offset_1:nx-offset,ny-offset-offset_1:ny-offset,k,:)
    Q(nx-offset_1:nx,1:offset,k,:) = Q(1+offset:1+offset+offset_1,ny-offset-offset_1:ny-offset,k,:)
    Q(1:offset,ny-1:ny,k,:) = Q(nx-offset-offset_1:nx-offset,1+offset:1+offset+offset_1,k,:)
    Q(nx-offset_1:nx,ny-1:ny,k,:) = Q(1+offset:1+offset+offset_1,1+offset:1+offset+offset_1,k,:)
  enddo

  do j = 1, ny
    do i = 1, nx
      ! 2nd
      !Q(i,j,1:1,:) = Q(i,j,nz-1:nz-1,:)
      !Q(i,j,nz:nz,:) = Q(i,j,2:2,:)
      ! 4th
      !Q(i,j,1:2,:) = Q(i,j,nz-3:nz-2,:)
      !Q(i,j,nz-1:nz,:) = Q(i,j,3:4,:)
      ! general
      Q(i,j,1:offset,:) = Q(i,j,nz-offset-offset_1:nz-offset,:)
      Q(i,j,nz-offset_1:nz,:) = Q(i,j,1+offset:1+offset+offset_1,:)
    enddo
  enddo
  end subroutine
end module set_bc
