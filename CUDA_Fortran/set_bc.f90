module set_bc
  implicit none

  interface set_cyclic_bc
    module procedure set_cyclic_bc2, set_cyclic_bc4
  end interface

  interface set_cyclic_bc_d
    module procedure set_cyclic_bc_d2, set_cyclic_bc_d4
  end interface
  
contains
  subroutine set_cyclic_bc2(id,nx,ny,nz,Q)
    integer(kind=2), intent(in) :: id
    integer, intent(in) :: nx, ny, nz
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        Q(1,j,k,:) = Q(nx-1,j,k,:)
        Q(nx,j,k,:) = Q(2,j,k,:)
      enddo
    enddo

    do k = 2, nz-1
      do i = 2, nx-1
        Q(i,1,k,:) = Q(i,ny-1,k,:)
        Q(i,ny,k,:) = Q(i,2,k,:)
      enddo
    enddo

    do k = 2, nz-1
      Q(1,1,k,:) = Q(nx-1,ny-1,k,:)
      Q(nx,1,k,:) = Q(2,ny-1,k,:)
      Q(1,ny,k,:) = Q(nx-1,2,k,:)
      Q(nx,ny,k,:) = Q(2,2,k,:)
    enddo

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1,:) = Q(i,j,nz-1,:)
        Q(i,j,nz,:) = Q(i,j,2,:)
      enddo
    enddo
  end subroutine

  subroutine set_cyclic_bc4(id,nx,ny,nz,Q)
    integer(kind=4), intent(in) :: id
    integer, intent(in) :: nx, ny, nz
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k

    do k = 3, nz-2
      do j = 3, ny-2
        Q(1:2,j,k,:) = Q(nx-3:nx-2,j,k,:)
        Q(nx-1:nx,j,k,:) = Q(3:4,j,k,:)
      enddo
    enddo

    do k = 3, nz-2
      do i = 3, nx-2
        Q(i,1:2,k,:) = Q(i,ny-3:ny-2,k,:)
        Q(i,ny-1:ny,k,:) = Q(i,3:4,k,:)
      enddo
    enddo

    do k = 3, nz-2
      Q(1:2,1:2,k,:) = Q(nx-3:nx-2,ny-3:ny-2,k,:)
      Q(nx-1:nx,1:2,k,:) = Q(3:4,ny-3:ny-2,k,:)
      Q(1:2,ny-1:ny,k,:) = Q(nx-3:nx-2,3:4,k,:)
      Q(nx-1:nx,ny-1:ny,k,:) = Q(3:4,3:4,k,:)
    enddo

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1:2,:) = Q(i,j,nz-3:nz-2,:)
        Q(i,j,nz-1:nz,:) = Q(i,j,3:4,:)
      enddo
    enddo
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_cyclic_bc_d2(id,nx,ny,nz,Q)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    integer i, j, k
    !$acc kernels deviceptr(Q)
    !$acc loop independent
    do k = 2, nz-1
      !$acc loop independent
      do j = 2, ny-1
        Q(1,j,k,:) = Q(nx-1,j,k,:)
        Q(nx,j,k,:) = Q(2,j,k,:)
      enddo
    enddo

    !$acc loop independent
    do k = 2, nz-1
      !$acc loop independent
      do i = 2, nx-1
        Q(i,1,k,:) = Q(i,ny-1,k,:)
        Q(i,ny,k,:) = Q(i,2,k,:)
      enddo
    enddo

    !$acc loop independent
    do k = 2, nz-1
      Q(1,1,k,:) = Q(nx-1,ny-1,k,:)
      Q(nx,1,k,:) = Q(2,ny-1,k,:)
      Q(1,ny,k,:) = Q(nx-1,2,k,:)
      Q(nx,ny,k,:) = Q(2,2,k,:)
    enddo

    !$acc loop independent
    do j = 1, ny
      !$acc loop independent
      do i = 1, nx
        Q(i,j,1,:) = Q(i,j,nz-1,:)
        Q(i,j,nz,:) = Q(i,j,2,:)
      enddo
    enddo
    !$acc end kernels
  end subroutine

  subroutine set_cyclic_bc_d4(id,nx,ny,nz,Q)
    integer(kind=4), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    integer i, j, k
  
    !$acc kernels deviceptr(Q)
    !$acc loop independent
    do k = 3, nz-2
      !$acc loop independent
      do j = 3, ny-2
        Q(1:2,j,k,:) = Q(nx-3:nx-2,j,k,:)
        Q(nx-1:nx,j,k,:) = Q(3:4,j,k,:)
      enddo
    enddo
  
    !$acc loop independent
    do k = 3, nz-2
      !$acc loop independent
      do i = 3, nx-2
        Q(i,1:2,k,:) = Q(i,ny-3:ny-2,k,:)
        Q(i,ny-1:ny,k,:) = Q(i,3:4,k,:)
      enddo
    enddo
  
    !$acc loop independent
    do k = 3, nz-2
      Q(1:2,1:2,k,:) = Q(nx-3:nx-2,ny-3:ny-2,k,:)
      Q(nx-1:nx,1:2,k,:) = Q(3:4,ny-3:ny-2,k,:)
      Q(1:2,ny-1:ny,k,:) = Q(nx-3:nx-2,3:4,k,:)
      Q(nx-1:nx,ny-1:ny,k,:) = Q(3:4,3:4,k,:)
    enddo
  
    !$acc loop independent
    do j = 1, ny
      !$acc loop independent
      do i = 1, nx
        Q(i,j,1:2,:) = Q(i,j,nz-3:nz-2,:)
        Q(i,j,nz-1:nz,:) = Q(i,j,3:4,:)
      enddo
    enddo
    !$acc end kernels
  end subroutine
end module set_bc
