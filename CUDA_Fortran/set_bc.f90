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
end module set_bc
