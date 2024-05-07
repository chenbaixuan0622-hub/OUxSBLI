module set_coordinate
  use mod_globals, only : nx, ny, nz, dxi => dx, deta => dy
  implicit none
contains
  subroutine set_xix(dx,xix)
    real(8), intent(in) :: dx(nx)
    real(8), intent(out) :: xix(nx)
    integer i
    do i = 1, nx
      xix(i) = 1.d0!dxi / dx(i)
    enddo
  end subroutine set_xix

  subroutine set_etay(dy,etay)
    real(8), intent(in) :: dy(ny)
    real(8), intent(out) :: etay(ny)
    integer j
    do j = 1, ny
      etay(j) = 1.d0!deta / dy(j)
    enddo
  end subroutine set_etay

  subroutine set_Jacobian(dx,dy,Jacobian)
    real(8), intent(in) :: dx(nx), dy(ny)
    real(8), intent(out) :: Jacobian(nx,ny)
    integer i, j
    do j = 1, ny
      do i = 1, nx
        Jacobian(i,j) = 1.d0!dx(i) * dy(j) / (dxi * deta)
    enddo;enddo
  end subroutine set_Jacobian
end module set_coordinate

