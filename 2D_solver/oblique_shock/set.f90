module set
  use cudafor
  use mod_globals, only : nx, ny, gamma, u0
  implicit none
  real(8) :: R = 287.03d0
  real(8) :: M0 = 1.9d0
  real(8) :: p0 = 14924.d0
  real(8) :: T0 = 171.31d0
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,4) :: Q
    real(8) rho0
    rho0 = p0 / (R * T0)
    Q(:,:,1) = rho0
    Q(:,:,2) = rho0 * u0
    Q(:,:,3) = 0.d0
    Q(:,:,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
  end subroutine set_init

  subroutine set_bc(Q)
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k
    integer :: No = int(0.4 * nx)
    real(8) beta, theta, Ms, Ms2, a1, rho0, rho2, p2, u1, u2, v1, v2,u_magnitude
    beta = dacos(-1.d0) * 37.2d0 / 180.d0
    Ms = M0 * dsin(beta)
    Ms2 = Ms**2
    theta = datan(2.d0 * (1.d0 / dtan(beta)) * (Ms2 - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
    rho0 = p0 / (R * T0)
    rho2 = rho0 * (gamma + 1.d0) * Ms2 / ((gamma - 1.d0) * Ms2 + 2.d0)
    p2 = p0 * (1.d0 + 2.d0 * gamma * (Ms2 - 1.d0) / (gamma + 1.d0))
    u1 = u0 * dsin(beta)
    v1 = u0 * dcos(beta)
    a1 = u0 / M0
    u2 = u1 - 2.d0 * a1 * (Ms - 1.d0 / Ms) / (gamma + 1.d0)
    v2 = u0 * dcos(beta)
    u_magnitude = sqrt(u2**2 + v2**2)
    write(*,"(a, f9.4)") "beta", beta  
    write(*,"(a, f9.4)") "Ms2", Ms2
    write(*,"(a, f9.4)") "theta", theta  
    write(*,"(a, f9.4)") "a1", a1
    write(*,"(a, f9.4)") "rho2", rho2  
    write(*,"(a, f11.4)") "p1", p0  
    write(*,"(a, f11.4)") "p2", p2  
    write(*,"(a, f9.4)") "u1", u1
    write(*,"(a, f9.4)") "v1", v1  
    write(*,"(a, f9.4)") "u2", u2  
    write(*,"(a, f9.4)") "v2", v2  
    write(*,"(a, f9.4)") "uuu2", u_magnitude
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
      Q(i,ny,2) = rho2 * u_magnitude * dcos(theta) 
      Q(i,ny,3) = - rho2 * u_magnitude * dsin(theta)
      Q(i,ny,4) = p2 / (gamma - 1.d0) + 0.5d0 * (Q(i,ny,2)**2 + Q(i,ny,3)**2) / rho2
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

