module set_coordinate
  use mod_globals, only : nx, ny
  implicit none
contains
  subroutine set_xix(dx,xix)
    real(8), intent(in)   :: dx(nx-1)
    real(8), intent(out)  :: xix(nx-1)
    integer i
    do i = 1, nx-1
      xix(i) = 1.d0 / dx(i)
    enddo
  end subroutine set_xix

  subroutine set_etay(dy,etay)
    real(8), intent(in)   :: dy(ny-1)
    real(8), intent(out)  :: etay(ny-1)
    integer j
    do j = 1, ny-1
      etay(j) = 1.d0 / dy(j)
    enddo
  end subroutine set_etay

  subroutine set_Jacobian(dx,dy,Jacobian)
    real(8), intent(in)  :: dx(nx-1), dy(ny-1)
    real(8), intent(out) :: Jacobian(nx,ny)
    integer i, j
    do j = 2, ny-1
      do i = 2, nx-1
        Jacobian(i,j) = 4.d0 / &
        & ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)))
    enddo;enddo
    do j = 2, ny-1
      Jacobian(1,j) = Jacobian(2,j)
    enddo
    do j = 2, ny-1
      Jacobian(nx,j) = Jacobian(nx-1,j)
    enddo
    do i = 1, nx
      Jacobian(i,1) = Jacobian(i,2)
    enddo
    do i = 1, nx
      Jacobian(i,ny) = Jacobian(i,ny-1)
    enddo
  end subroutine set_Jacobian
end module set_coordinate

