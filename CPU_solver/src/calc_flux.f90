module calc_flux
  use mod_globals, only : gamma, id_accuracy
  use calc_keep
  implicit none
contains
  subroutine calc_E(nx, ny, nz, rho, u, v, w, p, E)
    integer, intent(in)                      :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out)                     :: E(nx-1,ny-2,nz-2,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8)                 :: Normal(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 1, nx-1
          if (3 <= i .and. i <= nx-3 .and. 8 == kind(id_accuracy)) then
            rho6(:) = rho(i-2:i+3,j,k)
            p6(:)   =   p(i-2:i+3,j,k)
            V6(:,1) =   u(i-2:i+3,j,k)
            V6(:,2) =   v(i-2:i+3,j,k)
            V6(:,3) =   w(i-2:i+3,j,k)
            E(i,j-1,k-1,:) = KEEP6(1,rho6,p6,V6,Normal)
          elseif (2 <= i .and. i <= nx-2 .and. 4 <= kind(id_accuracy)) then
            rho4(:) = rho(i-1:i+2,j,k)
            p4(:)   =   p(i-1:i+2,j,k)
            V4(:,1) =   u(i-1:i+2,j,k)
            V4(:,2) =   v(i-1:i+2,j,k)
            V4(:,3) =   w(i-1:i+2,j,k)
            E(i,j-1,k-1,:) = KEEP4(1,rho4,p4,V4,Normal)
          else
            rho2(:) = rho(i:i+1,j,k)
            p2(:)   =   p(i:i+1,j,k)
            V2(:,1) =   u(i:i+1,j,k)
            V2(:,2) =   v(i:i+1,j,k)
            V2(:,3) =   w(i:i+1,j,k)
            E(i,j-1,k-1,:) = KEEP2(1,rho2,p2,V2,Normal)
          endif
    enddo;enddo;enddo
  end subroutine calc_E

  subroutine calc_F(nx, ny, nz, rho, u, v, w, p, F)
    integer, intent(in)                      :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out)                     :: F(nx-2,ny-1,nz-2,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8)                 :: Normal(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    do k = 2, nz-1
      do j = 1, ny-1
        do i = 2, nx-1
          if (3 <= j .and. j <= ny-3 .and. 8 == kind(id_accuracy)) then
            rho6(:) = rho(i,j-2:j+3,k)
            p6(:)   =   p(i,j-2:j+3,k)
            V6(:,1) =   u(i,j-2:j+3,k)
            V6(:,2) =   v(i,j-2:j+3,k)
            V6(:,3) =   w(i,j-2:j+3,k)
            F(i-1,j,k-1,:) = KEEP6(2,rho6,p6,V6,Normal)
          elseif (2 <= j .and. j <= ny-2 .and. 4 <= kind(id_accuracy)) then
            rho4(:) = rho(i,j-1:j+2,k)
            p4(:)   =   p(i,j-1:j+2,k)
            V4(:,1) =   u(i,j-1:j+2,k)
            V4(:,2) =   v(i,j-1:j+2,k)
            V4(:,3) =   w(i,j-1:j+2,k)
            F(i-1,j,k-1,:) = KEEP4(2,rho4,p4,V4,Normal)
          else
            rho2(:) = rho(i,j:j+1,k)
            p2(:)   =   p(i,j:j+1,k)
            V2(:,1) =   u(i,j:j+1,k)
            V2(:,2) =   v(i,j:j+1,k)
            V2(:,3) =   w(i,j:j+1,k)
            F(i-1,j,k-1,:) = KEEP2(2,rho2,p2,V2,Normal)
          endif
    enddo;enddo;enddo
  end subroutine calc_F

  subroutine calc_G(nx, ny, nz, rho, u, v, w, p, G)
    integer, intent(in)                      :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out)                     :: G(nx-2,ny-2,nz-1,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8)                 :: Normal(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
    do k = 1, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          if (3 <= k .and. k <= nz-3 .and. 8 == kind(id_accuracy)) then
            rho6(:) = rho(i,j,k-2:k+3)
            p6(:)   =   p(i,j,k-2:k+3)
            V6(:,1) =   u(i,j,k-2:k+3)
            V6(:,2) =   v(i,j,k-2:k+3)
            V6(:,3) =   w(i,j,k-2:k+3)
            G(i-1,j-1,k,:) = KEEP6(3,rho6,p6,V6,Normal)
          elseif (2 <= k .and. k <= nz-2 .and. 4 <= kind(id_accuracy)) then
            rho4(:) = rho(i,j,k-1:k+2)
            p4(:)   =   p(i,j,k-1:k+2)
            V4(:,1) =   u(i,j,k-1:k+2)
            V4(:,2) =   v(i,j,k-1:k+2)
            V4(:,3) =   w(i,j,k-1:k+2)
            G(i-1,j-1,k,:) = KEEP4(3,rho4,p4,V4,Normal)
          else
            rho2(:) = rho(i,j,k:k+1)
            p2(:)   =   p(i,j,k:k+1)
            V2(:,1) =   u(i,j,k:k+1)
            V2(:,2) =   v(i,j,k:k+1)
            V2(:,3) =   w(i,j,k:k+1)
            G(i-1,j-1,k,:) = KEEP2(3,rho2,p2,V2,Normal)
          endif
    enddo;enddo;enddo
  end subroutine calc_G
end module calc_flux

