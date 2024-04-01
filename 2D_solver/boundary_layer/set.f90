module set
  use cudafor
  use mod_globals, only : nx, ny, gamma
  implicit none
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,4) :: Q
    Q(:,:,1) = 1.4d0
    Q(:,:,2:3) = 0.d0
    Q(:,:,4) = 1.d0 / (gamma - 1.d0) 
    ! inlet
    Q(1:2,2:ny-1,1) = 1.4d0
    Q(1:2,2:ny-1,2) = 1.4d0 * 3.d0
    Q(1:2,2:ny-1,3) = 0.d0
    Q(1:2,2:ny-1,4) = 1.d0 / (gamma - 1.d0) + 0.5d0 * 1.4d0 * 3.d0**2
  end subroutine set_init

  subroutine set_bc(Q)
    real(8), intent(inout), device :: Q(nx,ny,4)
    real(8), p_left, p_top
    integer i, j
    integer :: N1 = int(0.4 * nx)
    integer :: N2 = int(0.2 * ny)
    !$cuf kernel do <<<*,*>>>
    do j = 2, ny-1
      ! inlet
      Q(1,j,1) = 1.4d0
      Q(1,j,2) = 1.4d0 * 3.d0
      Q(1,j,3) = 0.d0
      Q(1,j,4) = 1.d0 / (gamma - 1.d0) + 0.5d0 * 1.4d0 * 3.d0**2
      ! outlet
      Q(nx,j,1) = Q(nx-1,j,1)
      Q(nx,j,2) = Q(nx-1,j,2)
      Q(nx,j,3) = Q(nx-1,j,3)
      Q(nx,j,4) = Q(nx-1,j,4)
    enddo

    !$cuf kernel do <<<*,*>>>
    do i = 1, nx
      ! bottom
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = 0.d0
      Q(i,2,3) = 0.d0
      Q(i,2,4) = Q(i,3,4)
      ! imaginary
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = Q(i,3,2)
      Q(i,1,3) = -Q(i,3,3)
      Q(i,1,4) = Q(i,3,4)
      ! top
      Q(i,ny,:) = Q(i,ny-1,:)
    enddo
  end subroutine set_bc
end module set

