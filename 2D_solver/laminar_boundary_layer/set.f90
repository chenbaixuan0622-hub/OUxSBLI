module set
  use cudafor
  use mod_globals, only : nx, ny, gamma, u0, R, p0, rho0
  implicit none
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,4) :: Q
    integer i
    Q(:,:,1) = rho0
    Q(:,:,2) = rho0 * u0
    Q(:,:,3) = 0.d0
    Q(:,:,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    ! slip
    do i = 1, int(0.25 * nx)
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = Q(i,3,2)
      Q(i,2,3) = 0.d0
      Q(i,2,4) = Q(i,3,4)
      ! imaginary
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = Q(i,3,2) 
      Q(i,1,3) = -Q(i,3,3) 
      Q(i,1,4) = Q(i,3,4) 
    enddo
    ! noSlip
    do i = int(0.25 * nx), nx
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = 0.d0
      Q(i,2,3) = 0.d0
      Q(i,2,4) = Q(i,3,4)
      ! imaginray
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = -Q(i,3,2) 
      Q(i,1,3) = -Q(i,3,3) 
      Q(i,1,4) = Q(i,3,4) 
    enddo
  end subroutine set_init

  subroutine set_bc(Q)
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k
    !$cuf kernel do <<<*,*>>>
    do j = 2, ny-1
      ! inlet
      Q(1,j,1) = rho0
      Q(1,j,2) = rho0 * u0
      Q(1,j,3) = 0.d0
      Q(1,j,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    enddo

    !$cuf kernel do <<<*,*>>>
    do k = 1, 4
      do j = 2, ny-1
        ! outlet
        Q(nx,j,k) = Q(nx-1,j,k)
      enddo
    enddo

    !$cuf kernel do <<<*,*>>>
    do k = 1, 4
      do i = 1, nx
        ! top
        Q(i,ny,k) = Q(i,ny-1,k)
      enddo
    enddo
    
    ! slip
    !$cuf kernel do <<<*,*>>>
    do i = 1, int(0.25 * nx)
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = Q(i,3,2)
      Q(i,2,3) = 0.d0
      Q(i,2,4) = Q(i,3,4)
      ! imaginary
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = Q(i,3,2) 
      Q(i,1,3) = -Q(i,3,3) 
      Q(i,1,4) = Q(i,3,4) 
    enddo
    ! noSlip
    !$cuf kernel do <<<*,*>>>
    do i = int(0.25 * nx), nx
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = 0.d0
      Q(i,2,3) = 0.d0
      Q(i,2,4) = Q(i,3,4)
      ! imaginray
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = -Q(i,3,2) 
      Q(i,1,3) = -Q(i,3,3) 
      Q(i,1,4) = Q(i,3,4) 
    enddo
  end subroutine set_bc
end module set

