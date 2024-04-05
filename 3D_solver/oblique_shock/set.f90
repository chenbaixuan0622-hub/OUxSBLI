module set
  use mod_globals, only : nx, ny, nz, dy, gamma, u0, &
  & beta, theta, Ms, Ms2, a1, rho0, rho2, p0, p2, u1, u2, v1, v2, u_magnitude, ux, uy
  implicit none
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    Q(:,:,:,1) = rho0
    Q(:,:,:,2) = rho0 * u0
    Q(:,:,:,3) = 0.d0
    Q(:,:,:,4) = 0.d0
    Q(:,:,:,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
  end subroutine set_init
  
  subroutine set_bc(id_accuracy,Q,T)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), device :: T(nx,ny,nz)
    integer i, j, k, l
    integer :: No = int(0.4 * nx)
    ! inlet
    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do j = 1, ny-1
        Q(1,j,k,1) = rho0
        Q(1,j,k,2) = rho0 * u0
        Q(1,j,k,3) = 0.d0
        Q(1,j,k,4) = 0.d0
        Q(1,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
    enddo

    ! outlet
    !$cuf kernel do<<<*,*>>>
    do l = 1, 5
      do k = 1, nz
        do j = 1, ny-1
          Q(nx,j,k,l) = Q(nx-1,j,k,l)
        enddo
      enddo
    enddo

    ! top
    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do i = 1, No
        Q(i,ny,k,1) = rho0
        Q(i,ny,k,2) = rho0 * u0 
        Q(i,ny,k,3) = 0.d0
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
    enddo

    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do i = No+1, nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * ux
        Q(i,ny,k,3) = rho2 * uy
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
      enddo
    enddo

    ! bottom
    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        ! slip wall
        Q(i,2,k,1) = Q(i,3,k,1)
        Q(i,2,k,2) = Q(i,3,k,2)
        Q(i,2,k,3) = 0.d0
        Q(i,2,k,4) = Q(i,3,k,4)
        Q(i,2,k,5) = Q(i,3,k,5)
        ! imaginary
        Q(i,1,k,1) = Q(i,3,k,1)
        Q(i,1,k,2) = Q(i,3,k,2)
        Q(i,1,k,3) = -Q(i,3,k,3)
        Q(i,1,k,4) = Q(i,3,k,4)
        Q(i,1,k,5) = Q(i,3,k,5)
      enddo
    enddo

    ! span
    !$cuf kernel do<<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q(i,j,1,l) = Q(i,j,2,l)
          Q(i,j,nz,l) = Q(i,j,nz-1,l)
        enddo
      enddo
    enddo
  end subroutine set_bc
end module set

