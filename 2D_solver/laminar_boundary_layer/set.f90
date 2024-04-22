module set
  use cudafor
  use mod_globals, only : nx, ny, dx, dy, gamma, u0, R, p0, rho0
  implicit none
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,4) :: Q
    real(8) x, y, d, u, p_wall
    real(8) :: nu0 = 3.8206d-5
    integer i, j
    Q(:,:,1) = rho0
    Q(:,:,3) = 0.d0
    ! boundary layer
    do j = 3, ny-1
      do i = 1, nx
        x = (i-1) * dx
        y = (j-2) * dy
        d = 5.d0 * sqrt(nu0 * x / u0)
        u = u0 * (1.5d0 * (y / d) - 0.5d0 * (y / d)**3) 
        Q(i,j,2) = min(u0, u) 
        Q(i,j,4) = p0 / (gamma - 1.d0) + 0.5d0 * (Q(i,j,2)**2 + Q(i,j,3)**2) / Q(i,j,1) 
      enddo
    enddo

    ! noSlip
    do i = 1, nx
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = 0.d0
      Q(i,2,3) = 0.d0
      p_wall = (gamma - 1.d0) * (Q(i,3,4) - 0.5d0 * (Q(i,3,2)**2 + Q(i,3,3)**2) / Q(i,3,1))
      Q(i,2,4) = p_wall / (gamma - 1.d0)
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
    real(8) p_wall
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
    
    ! noSlip
    !$cuf kernel do <<<*,*>>>
    do i = 1, nx
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = 0.d0
      Q(i,2,3) = 0.d0
      p_wall = (gamma - 1.d0) * (Q(i,3,4) - 0.5d0 * (Q(i,3,2)**2 + Q(i,3,3)**2) / Q(i,3,1))
      Q(i,2,4) = p_wall / (gamma - 1.d0)
      ! imaginray
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = -Q(i,3,2) 
      Q(i,1,3) = -Q(i,3,3) 
      Q(i,1,4) = Q(i,3,4) 
    enddo
  end subroutine set_bc
end module set

