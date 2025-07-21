module set_coordinate
  implicit none
contains
  subroutine set_xix(nx, dx, xix)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: dx(nx-1)
    real(8), intent(out) :: xix(nx-1)
    integer i
    do i = 1, nx-1
      xix(i) = 1.d0 / dx(i)
    enddo
  end subroutine set_xix

  subroutine set_etay(ny, dy, etay)
    integer, intent(in)  :: ny
    real(8), intent(in)  :: dy(ny-1)
    real(8), intent(out) :: etay(ny-1)
    integer j
    do j = 1, ny-1
      etay(j) = 1.d0 / dy(j)
    enddo
  end subroutine set_etay

  subroutine set_zetaz(nz, dz, zetaz)
    integer, intent(in)  :: nz
    real(8), intent(in)  :: dz(nz-1)
    real(8), intent(out) :: zetaz(nz-1)
    integer k
    do k = 1, nz-1
      zetaz(k) = 1.d0 / dz(k)
    enddo
  end subroutine set_zetaz

  subroutine set_Jacobian_y(nx, ny, nz, dx, dy, dz, Jacobian)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(out) :: Jacobian(ny)
    integer i, j, k
    if (2 < nz) then
      do k = 2, nz-1
        do j = 2, ny-1
          do i = 2, nx-1
            Jacobian(j) = 8.d0 / &
            & ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * (dz(k-1) + dz(k)))
      enddo;enddo;enddo
    else
      do j = 2, ny-1
        do i = 2, nx-1
          Jacobian(j) = 4.d0 / &
          & ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)))
      enddo;enddo
    endif
    Jacobian(1)  = Jacobian(2)
    Jacobian(ny) = Jacobian(ny-1)
  end subroutine set_Jacobian_y
end module set_coordinate

