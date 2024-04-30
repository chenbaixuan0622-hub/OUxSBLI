module set_coordinate
  use mod_globals, only : nx, ny, nz, dxi => dx, deta => dy
  implicit none
  interface set_xix
    module procedure set_xix_2D, set_xix_3D
  end interface

  interface set_etay
    module procedure set_etay_2D, set_etay_3D
  end interface

  interface set_Jacobian
    module procedure set_Jacobian_2D, set_Jacobian_3D
  end interface
contains
  subroutine set_xix_2D(x,xix)
    real(8), intent(in) :: x(nx,ny)
    real(8), intent(out) :: xix(nx,ny)
    integer i, j
    xix(1,:) = dxi / (-x(1,:) + x(2,:)) 
    do j = 1, ny
      do i = 2, nx-1
        xix(i,j) = 2.d0 * dxi / (-x(i-1,j) + x(i+1,j))
      enddo
    enddo
    xix(nx,:) = dxi / (-x(nx-1,:) + x(nx,:))
    xix(:,:) = 1.d0
  end subroutine set_xix_2D

  subroutine set_xix_3D(x,xix)
    real(8), intent(in) :: x(nx,ny,nz)
    real(8), intent(out) :: xix(nx,ny,nz)
    integer i, j, k
    xix(1,:,:) = dxi / (-x(1,:,:) + x(2,:,:)) 
    do k = 1, nz
      do j = 1, ny
        do i = 2, nx-1
          xix(i,j,k) = 2.d0 * dxi / (-x(i-1,j,k) + x(i+1,j,k))
        enddo
      enddo
    enddo
    xix(nx,:,:) = dxi / (-x(nx-1,:,:) + x(nx,:,:))
    xix(:,:,:) = 1.d0
  end subroutine set_xix_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_etay_2D(y,etay)
    real(8), intent(in) :: y(nx,ny)
    real(8), intent(out) :: etay(nx,ny)
    integer i, j
    etay(:,1) = deta / (-y(:,1) + y(:,2)) 
    do j = 2, ny-1
      do i = 1, nx
        etay(i,j) = 2.d0 * deta / (-y(i,j-1) + y(i,j+1))
      enddo
    enddo
    etay(:,ny) = deta / (-y(:,ny-1) + y(:,ny))
    etay(:,:) = 1.d0
  end subroutine set_etay_2D

  subroutine set_etay_3D(y,etay)
    real(8), intent(in) :: y(nx,ny,nz)
    real(8), intent(out) :: etay(nx,ny,nz)
    integer i, j, k
    etay(:,1,:) = deta / (-y(:,1,:) + y(:,2,:)) 
    do k = 1, nz
      do j = 2, ny-1
        do i = 1, nx
          etay(i,j,k) = 2.d0 * deta / (-y(i,j-1,k) + y(i,j+1,k))
        enddo
      enddo
    enddo
    etay(:,ny,:) = deta / (-y(:,ny-1,:) + y(:,ny,:))
    etay(:,:,:) = 1.d0
  end subroutine set_etay_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_Jacobian_2D(x,y,xix,etay,Jacobian)
    real(8), intent(in), dimension(nx,ny) :: x, y, xix, etay
    real(8), intent(out) :: Jacobian(nx,ny)
    integer i, j
    do j = 2, ny-1
      Jacobian(1,j) = 0.5d0 * (-x(1,j) + x(2,j)) * (-y(1,j-1) + y(1,j+1)) / (dxi * deta)
      Jacobian(nx,j) = 0.5d0 * (-x(nx-1,j) + x(nx,j)) * (-y(nx,j-1) + y(nx,j+1)) / (dxi * deta)
    enddo
    do i = 2, nx-1
      Jacobian(i,1) = 0.5d0 * (-x(i-1,1) + x(i+1,1)) * (-y(i,1) + y(i,2)) / (dxi * deta)
      Jacobian(i,ny) = 0.5d0 * (-x(i-1,ny) + x(i+1,ny)) * (-y(i,ny-1) + y(i,ny)) / (dxi * deta)
    enddo
    Jacobian(1,1) = (-x(1,1) + x(2,1)) * (-y(1,1) + y(1,2)) / (dxi * deta)
    Jacobian(nx,1) = (-x(nx-1,1) + x(nx,1)) * (-y(nx,1) + y(nx,2)) / (dxi * deta)
    Jacobian(1,ny) = (-x(1,ny) + x(2,ny)) * (-y(1,ny-1) + y(1,ny)) / (dxi * deta)
    Jacobian(nx,ny) = (-x(nx-1,ny) + x(nx,ny)) * (-y(nx,ny-1) + y(nx,ny)) / (dxi * deta)
    do j = 2, ny-1
      do i = 2, nx-1
        Jacobian(i,j) = 0.25d0 * (-x(i-1,j) + x(i+1,j)) * (-y(i,j-1) + y(i,j+1)) / (dxi * deta)
      enddo
    enddo
    Jacobian(:,:) = 1.d0
  end subroutine set_Jacobian_2D

  subroutine set_Jacobian_3D(x,y,xix,etay,Jacobian)
    real(8), intent(in), dimension(nx,ny,nz) :: x, y, xix, etay
    real(8), intent(out) :: Jacobian(nx,ny,nz)
    integer i, j, k
    do k = 1, nz
      do j = 2, ny-1
        Jacobian(1,j,k) = 0.5d0 * (-x(1,j,k) + x(2,j,k)) * (-y(1,j-1,k) + y(1,j+1,k)) / (dxi * deta)
        Jacobian(nx,j,k) = 0.5d0 * (-x(nx-1,j,k) + x(nx,j,k)) * (-y(nx,j-1,k) + y(nx,j+1,k)) / (dxi * deta)
      enddo
    enddo
    do k = 1, nz
      do i = 2, nx-1
        Jacobian(i,1,k) = 0.5d0 * (-x(i-1,1,k) + x(i+1,1,k)) * (-y(i,1,k) + y(i,2,k)) / (dxi * deta)
        Jacobian(i,ny,k) = 0.5d0 * (-x(i-1,ny,k) + x(i+1,ny,k)) * (-y(i,ny-1,k) + y(i,ny,k)) / (dxi * deta)
      enddo
    enddo
    do k = 1, nz
      Jacobian(1,1,k) = (-x(1,1,k) + x(2,1,k)) * (-y(1,1,k) + y(1,2,k)) / (dxi * deta)
      Jacobian(nx,1,k) = (-x(nx-1,1,k) + x(nx,1,k)) * (-y(nx,1,k) + y(nx,2,k)) / (dxi * deta)
      Jacobian(1,ny,k) = (-x(1,ny,k) + x(2,ny,k)) * (-y(1,ny-1,k) + y(1,ny,k)) / (dxi * deta)
      Jacobian(nx,ny,k) = (-x(nx-1,ny,k) + x(nx,ny,k)) * (-y(nx,ny-1,k) + y(nx,ny,k)) / (dxi * deta)
    enddo
    do k = 1, nz
      do j = 2, ny-1
        do i = 2, nx-1
          Jacobian(i,j,k) = 0.25d0 * (-x(i-1,j,k) + x(i+1,j,k)) * (-y(i,j-1,k) + y(i,j+1,k)) / (dxi * deta)
        enddo
      enddo
    enddo
    Jacobian(:,:,:) = 1.d0
  end subroutine set_Jacobian_3D
end module set_coordinate

