module set
  use mod_globals, only : nx, ny, nz, dy, gamma, u0, &
  & beta, theta, Ms, Ms2, a1, rho0, rho2, p0, p2, u1, u2, v1, v2, u_magnitude, ux, uy
  implicit none
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    integer i, k
    integer :: No = int(0.4 * nx)
    Q(:,:,:,1) = rho0
    Q(:,:,:,2) = rho0 * u0
    Q(:,:,:,3) = 0.d0
    Q(:,:,:,4) = 0.d0
    Q(:,:,:,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2

    ! top
    do k = 2, nz-1
      do i = 1, No
        Q(i,ny,k,1) = rho0
        Q(i,ny,k,2) = rho0 * u0 
        Q(i,ny,k,3) = 0.d0
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
    enddo

    do k = 2, nz-1
      do i = No+1, nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * ux
        Q(i,ny,k,3) = rho2 * uy
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
      enddo
    enddo

    ! bottom
    Q(:,2,:,1) = Q(:,3,:,1)
    Q(:,2,:,2) = Q(:,3,:,2)
    Q(:,2,:,3) = 0.d0
    Q(:,2,:,4) = Q(:,3,:,4)
    Q(:,2,:,5) = Q(:,3,:,5)
    ! imaginary
    Q(:,1,:,1) = Q(:,3,:,1)
    Q(:,1,:,2) = Q(:,3,:,2)
    Q(:,1,:,3) = -Q(:,3,:,3)
    Q(:,1,:,4) = Q(:,3,:,4)
    Q(:,1,:,5) = Q(:,3,:,5)
  end subroutine set_init
  
  subroutine set_bc(id_accuracy,Q,T)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), device :: T(nx,ny,nz)
    integer i, j, k, l
    integer :: No = int(0.4 * nx)
    !$cuf kernel do<<<*,*>>>
    do k = 2, nz-1
      do j = 3, ny-1
        ! inlet
        Q(1,j,k,1) = rho0
        Q(1,j,k,2) = rho0 * u0
        Q(1,j,k,3) = 0.d0
        Q(1,j,k,4) = 0.d0
        Q(1,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
        ! outlet
        Q(nx,j,k,1) = Q(nx-1,j,k,1)
        Q(nx,j,k,2) = Q(nx-1,j,k,2)
        Q(nx,j,k,3) = Q(nx-1,j,k,3)
        Q(nx,j,k,4) = Q(nx-1,j,k,4)
        Q(nx,j,k,5) = Q(nx-1,j,k,5)
      enddo
    enddo

    ! top
    !$cuf kernel do<<<*,*>>>
    do k = 2, nz-1
      do i = 1, No
        Q(i,ny,k,1) = rho0
        Q(i,ny,k,2) = rho0 * u0 
        Q(i,ny,k,3) = 0.d0
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
    enddo

    !$cuf kernel do<<<*,*>>>
    do k = 2, nz-1
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
    do k = 2, nz-1
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

    ! cyclic
    !$cuf kernel do<<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q(i,j,1,l) = Q(i,j,nz-3,l)
          Q(i,j,2,l) = Q(i,j,nz-2,l)
          Q(i,j,nz-1,l) = Q(i,j,3,l)
          Q(i,j,nz,l) = Q(i,j,4,l)
        enddo
      enddo
    enddo
  end subroutine set_bc
end module set

