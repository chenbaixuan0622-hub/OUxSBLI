module set_bc
  implicit none
  ! calc domain
  !!!!!!!!!!!!!!!!!z
  ! block ! block !
  !   1   !   4   !
  !       !       !
  !!!!!!!!!!!!!!!!!
  ! block ! block !
  !   2   !   3   !
  !       !       !
  !!!!!!!!!!!!!!!!!
  !y
contains
  subroutine set_cyclic_bc(nx,ny,nz,Q)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(inout) :: Q(nx,ny,nz,5,4)
    integer blockIdx, tb, lr, c !tb: top-bottom, lr: left-right, c: corner
    do blockIdx = 1, 4
      ! when blockIdx = 1, c = 3, tb = 2, and lr = 4
      ! when blockIdx = 2, c = 4, tb = 1, and lr = 3
      ! when blockIdx = 3, c = 1, tb = 4, and lr = 2
      ! when blockIdx = 4, c = 2, tb = 3, and lr = 1
      c = mod(blockIdx + 1,4) + 1
      tb = blockIdx - (-1) ** mod(blockIdx,2)
      lr = mod(3 + blockIdx + (-1) ** mod(blockIdx,2),4) + 1

      ! top
      Q(3:nx-2,1:2,3:nz-2,:,blockIdx) = Q(3:nx-2,ny-3:ny-2,3:nz-2,:,tb)
      ! bottom
      Q(3:nx-2,ny-1:ny,3:nz-2,:,blockIdx) = Q(3:nx-2,3:4,3:nz-2,:,tb)
      ! left
      Q(3:nx-2,3:ny-2,1:2,:,blockIdx) = Q(3:nx-2,3:ny-2,nz-3:nz-2,:,lr)
      ! right
      Q(3:nx-2,3:ny-2,nz-1:nz,:,blockIdx) = Q(3:nx-2,3:ny-2,3:4,:,lr)
      ! corner
      Q(3:nx-2,1:2,1:2,:,blockIdx) = Q(3:nx-2,ny-3:ny-2,nz-3:nz-2,:,c)
      Q(3:nx-2,1:2,nz-1:nz,:,blockIdx) = Q(3:nx-2,ny-3:ny-2,3:4,:,c)
      Q(3:nx-2,ny-1:ny,nz-1:nz,:,blockIdx) = Q(3:nx-2,3:4,3:4,:,c)
      Q(3:nx-2,ny-1:ny,1:2,:,blockIdx) = Q(3:nx-2,3:4,nz-3:nz-2,:,c)
      ! x direction
      Q(1:2,:,:,:,blockIdx) = Q(nx-3:nx-2,:,:,:,blockIdx)
      Q(nx-1:nx,:,:,:,blockIdx) = Q(3:4,:,:,:,blockIdx)
    enddo
  end subroutine

  subroutine set_cyclic_bc_d(nx,ny,nz,Q)
    integer, value :: nx, ny, nz
    real(8), device :: Q(:,:,:,:,:)
    integer i, j, k, blockId
    integer tb, lr, c !tb: top-bottom, lr: left-right, c: corner
    ! when blockIdx = 1, c = 3, tb = 2, and lr = 4
    ! when blockIdx = 2, c = 4, tb = 1, and lr = 3
    ! when blockIdx = 3, c = 1, tb = 4, and lr = 2
    ! when blockIdx = 4, c = 2, tb = 3, and lr = 1

    !$acc kernels deviceptr(Q)
    !$acc loop independent gang(4) private(tb)
    do blockId = 1, 4
    !$acc loop independent worker(16) private(tb)
      do k = 3, nz-2
      !$acc loop independent vector(32) private(tb)
        do i = 3, nx-2
          tb = blockId - (-1) ** mod(blockId,2)
          ! top
          Q(i,1:2,k,:,blockId) = Q(i,ny-3:ny-2,k,:,tb)
          ! bottom
          Q(i,ny-1:ny,k,:,blockId) = Q(i,3:4,k,:,tb)
        enddo
      enddo
    enddo
    
    !$acc loop independent gang(4) private(lr)
    do blockId = 1, 4
    !$acc loop independent worker(16) private(lr)
      do j = 3, ny-2
      !$acc loop independent vector(32) private(lr)
        do i = 3, nx-2
          lr = mod(3 + blockId + (-1) ** mod(blockId,2),4) + 1
          ! left
          Q(i,j,1:2,:,blockId) = Q(i,j,nz-3:nz-2,:,lr)
          ! right
          Q(i,j,nz-1:nz,:,blockId) = Q(i,j,3:4,:,lr)
        enddo
      enddo
    enddo
    
    !$acc loop independent gang(4) private(c)
    do blockId = 1, 4
    !$acc loop independent vector(64) private(c)
      do i = 3, nx-2
        c = mod(blockId + 1,4) + 1
        ! corner
        Q(i,1:2,1:2,:,blockId) = Q(i,ny-3:ny-2,nz-3:nz-2,:,c)
        Q(i,1:2,nz-1:nz,:,blockId) = Q(i,ny-3:ny-2,3:4,:,c)
        Q(i,ny-1:ny,nz-1:nz,:,blockId) = Q(i,3:4,3:4,:,c)
        Q(i,ny-1:ny,1:2,:,blockId) = Q(i,3:4,nz-3:nz-2,:,c)
      enddo
    enddo
    
    !$acc loop independent gang(4)
    do blockId = 1, 4
    !$acc loop independent worker(16)
      do k = 1, nz
      !$acc loop independent vector(32)
        do j = 1, ny
          ! x direction
          Q(1:2,j,k,:,blockId) = Q(nx-3:nx-2,j,k,:,blockId)
          Q(nx-1:nx,j,k,:,blockId) = Q(3:4,j,k,:,blockId)
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine
end module set_bc
