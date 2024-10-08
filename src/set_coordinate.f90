module set_coordinate
  implicit none
contains
  subroutine set_xix(nx,dx,xix)
    integer, intent(in)   :: nx
    real(8), intent(in)   :: dx(nx-1)
    real(8), intent(out)  :: xix(nx-1)
    integer i
    do i = 1, nx-1
      xix(i) = 1.d0 / dx(i)
    enddo
  end subroutine set_xix

  subroutine set_etay(ny,dy,etay)
    integer, intent(in)   :: ny
    real(8), intent(in)   :: dy(ny-1)
    real(8), intent(out)  :: etay(ny-1)
    integer j
    do j = 1, ny-1
      etay(j) = 1.d0 / dy(j)
    enddo
  end subroutine set_etay

  subroutine set_zetaz(nz,dz,zetaz)
    integer, intent(in)   :: nz
    real(8), intent(in)   :: dz(nz-1)
    real(8), intent(out)  :: zetaz(nz-1)
    integer k
    do k = 1, nz-1
      zetaz(k) = 1.d0 / dz(k)
    enddo
  end subroutine set_zetaz

  subroutine set_Jacobian_y(nx,ny,nz,dx,dy,dz,Jacobian)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(out) :: Jacobian(ny)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          Jacobian(j) = 8.d0 / &
          & ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * (dz(k-1) + dz(k)))
    enddo;enddo;enddo
    do k = 2, nz-1
      do i = 1, nx
        Jacobian(1)  = Jacobian(2)
        Jacobian(ny) = Jacobian(ny-1)
    enddo;enddo
  end subroutine set_Jacobian_y

  subroutine set_Jacobian(nx,ny,nz,dx,dy,dz,Jacobian)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(out) :: Jacobian(nx,ny,nz)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          Jacobian(i,j,k) = 8.d0 / &
          & ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * (dz(k-1) + dz(k)))
    enddo;enddo;enddo
    do k = 2, nz-1
      do j = 2, ny-1
        Jacobian(1,j,k)  = Jacobian(2,j,k)
        Jacobian(nx,j,k) = Jacobian(nx-1,j,k)
    enddo;enddo
    do k = 2, nz-1
      do i = 1, nx
        Jacobian(i,1,k)  = Jacobian(i,2,k)
        Jacobian(i,ny,k) = Jacobian(i,ny-1,k)
    enddo;enddo
    do j = 1, ny
      do i = 1, nx
        Jacobian(i,j,1)  = Jacobian(i,j,2)
        Jacobian(i,j,nz) = Jacobian(i,j,nz-1)
    enddo;enddo
  end subroutine set_Jacobian
end module set_coordinate

