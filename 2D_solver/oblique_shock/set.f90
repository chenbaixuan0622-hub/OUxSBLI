module set
  use cudafor
  use mod_globals, only : nx, ny, gamma, u0, R, M0, p0, T, &
  & beta, Ms, Ms2, theta, rho0, rho2, p2, u1, u2, v1, v2, u_magnitude, ux, uy
  implicit none
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,4) :: Q
    integer i
    integer :: No = int(0.4 * nx)
    write(*,"(a,f11.4)") "beta ", beta 
    write(*,"(a,f11.4)") "Ms2  ", Ms2
    write(*,"(a,f11.4)") "theta", theta
    write(*,"(a,f11.4)") "uu1  ", u1
    write(*,"(a,f11.4)") "vv1  ", v1
    write(*,"(a,f11.4)") "uu2  ", u2
    write(*,"(a,f11.4)") "vv2  ", v2
    write(*,"(a,f11.4)") "uuu2 ", u_magnitude
    write(*,"(a,f11.4)") "u2   ", ux
    write(*,"(a,f11.4)") "v2   ", uy
    Q(:,:,1) = rho0
    Q(:,:,2) = rho0 * u0
    Q(:,:,3) = 0.d0
    Q(:,:,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    ! top
    do i = 1, No
      Q(i,ny,1) = rho0 
      Q(i,ny,2) = rho0 * u0
      Q(i,ny,3) = 0.d0
      Q(i,ny,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    enddo
    do i = No+1, nx
      Q(i,ny,1) = rho2
      Q(i,ny,2) = rho2 * ux 
      Q(i,ny,3) = rho2 * uy
      Q(i,ny,4) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo
    
    ! bottom
    Q(:,2,1) = Q(:,3,1)
    Q(:,2,2) = Q(:,3,2)
    Q(:,2,3) = 0.d0
    Q(:,2,4) = Q(:,3,4)
    ! imaginary
    Q(:,1,1) = Q(:,3,1)
    Q(:,1,2) = Q(:,3,2)
    Q(:,1,3) = -Q(:,3,3)
    Q(:,1,4) = Q(:,3,4)
  end subroutine set_init

  subroutine set_bc(Q)
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k
    integer :: No = int(0.4 * nx)
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
    do i = 1, No
      Q(i,ny,1) = rho0 
      Q(i,ny,2) = rho0 * u0
      Q(i,ny,3) = 0.d0
      Q(i,ny,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    enddo

    !$cuf kernel do <<<*,*>>>
    do i = No+1, nx
      Q(i,ny,1) = rho2
      Q(i,ny,2) = rho2 * ux 
      Q(i,ny,3) = rho2 * uy
      Q(i,ny,4) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo

    !$cuf kernel do <<<*,*>>>
    do i = 1, nx
      ! bottom (slip wall)
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
  end subroutine set_bc
end module set

