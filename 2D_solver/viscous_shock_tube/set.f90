module set
  use cudafor
  use mod_globals, only : nx, ny, gamma
  implicit none
contains
  subroutine set_init(Q)
    use mod_globals, only : rhol, rhor, pl, pr
    real(8), intent(out), dimension(nx,ny,4) :: Q
    ! left half
    Q(1:int(0.5*nx),:,1) = rhol
    Q(1:int(0.5*nx),:,4) = pl / (gamma - 1.d0)
    ! right half
    Q(int(0.5*nx):nx,:,1) = rhor
    Q(int(0.5*nx):nx,:,4) = pr / (gamma - 1.d0)
    ! both (velocity = 0)
    Q(:,:,2) = 0.d0
    Q(:,:,3) = 0.d0
  end subroutine set_init

  subroutine set_bc(Q)
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k
    !$cuf kernel do <<<*,*>>>
    do k = 1, 4
      do j = 2, ny-1
        ! inlet
        Q(1,j,k) = Q(2,j,k)
        ! outlet
        Q(nx,j,k) = Q(nx-1,j,k)
      enddo
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
      Q(i,ny-1,1) = Q(i,ny-2,1)
      Q(i,ny-1,2) = 0.d0
      Q(i,ny-1,3) = 0.d0
      Q(i,ny-1,4) = Q(i,ny-2,4)
      ! imaginary
      Q(i,ny,1) = Q(i,ny-2,1)
      Q(i,ny,2) = Q(i,ny-2,2)
      Q(i,ny,3) = -Q(i,ny-2,3)
      Q(i,ny,4) = Q(i,ny-2,4)
    enddo
  end subroutine set_bc
end module set

