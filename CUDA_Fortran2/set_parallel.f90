module set_parallel
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
    subroutine set_index(ny,nz,nyp,nzp)
      integer, intent(in) :: ny, nz
      integer, intent(out) :: nyp, nzp
      nyp = int(0.5 * ny) + 2
      nzp = int(0.5 * nz) + 2
    end subroutine set_index
  
    subroutine split(nx,ny,nz,nyp,nzp,Q,Qp)
      integer, intent(in) :: nx, ny, nz
      integer, intent(in) :: nyp, nzp
      real(8), intent(in) :: Q(nx,ny,nz,5)
      real(8), intent(out) :: Qp(nx,nyp,nzp,5,4)
      ! split into 4 blocks
      Qp(:,:,:,:,1) = Q(:,1:nyp,1:nzp,:)
      Qp(:,:,:,:,2) = Q(:,ny-nyp+1:ny,1:nzp,:)
      Qp(:,:,:,:,3) = Q(:,ny-nyp+1:ny,nz-nzp+1:nz,:)
      Qp(:,:,:,:,4) = Q(:,1:nyp,nz-nzp+1:nz,:)
    end subroutine split
  
    subroutine merge(nx,ny,nz,nyp,nzp,Qp,Q)
      integer, intent(in) :: nx, ny, nz, nyp, nzp
      real(8), intent(in) :: Qp(nx,nyp,nzp,5,4)
      real(8), intent(out) :: Q(nx,ny,nz,5)
      ! merge blocks
      Q(:,1:nyp,1:nzp,:) = Qp(:,:,:,:,1)
      Q(:,ny-nyp+1:ny,1:nzp,:) = Qp(:,:,:,:,2)
      Q(:,ny-nyp+1:ny,nz-nzp+1:nz,:) = Qp(:,:,:,:,3)
      Q(:,1:nyp,nz-nzp+1:nz,:) = Qp(:,:,:,:,4)
    end subroutine merge
end module set_parallel