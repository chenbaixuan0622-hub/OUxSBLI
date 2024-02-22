module set_init
  implicit none
contains
  function linspace(x1, x2, n) result(x)
    real(8), intent(in) :: x1, x2
    integer, intent(in) :: n
    integer i
    real(8) x(n)
    x = x1 + (x1 + x2) * (/ (dble(i - 1) / dble(n - 1), i = 1, n) /)
  end function linspace
  
  subroutine TaylorGreen(nx,ny,nz,nyp,nzp,gamma,RHO,L,M,Q)
    ! 4-th order accuracy
    use set_parallel
    use set_bc
    integer, intent(in) :: nx, ny, nz, nyp, nzp
    real(8), intent(in) :: gamma, RHO, L, M
    real(8), intent(out), dimension(nx,nyp,nzp,5,4) :: Q
    integer i, j, k
    real(8) pi
    real(8) x(nx-4), y(ny-4), z(nz-4)
    real(8) Q0(nx,ny,nz,5)
    pi = 2.d0 * acos(0.d0)
    x = linspace(0.d0, 2.d0 * pi, nx-4)
    y = linspace(0.d0, 2.d0 * pi, ny-4)
    z = linspace(0.d0, 2.d0 * pi, nz-4)
    do k = 3, nz-2
      do j = 3, ny-2
        do i = 3, nx-2
          ! rho
          Q0(i,j,k,1) = RHO
          ! rho u
          Q0(i,j,k,2) = RHO * M * sin(x(i-2)) * cos(y(j-2)) * cos(z(k-2))
          ! rho v
          Q0(i,j,k,3) = - RHO * M * cos(x(i-2)) * sin(y(j-2)) * cos(z(k-2))
          ! rho w
          Q0(i,j,k,4) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q0(i,j,k,5) = (1.d0/gamma+RHO*M**2*(cos(2.d0*x(i-2)/L)+cos(2.d0*y(j-2)/L))*(cos(2.d0*z(k-2)/L)+2.d0)/16.d0)&
          / (gamma - 1.d0) + 0.5d0 * (Q0(i,j,k,2) ** 2 + Q0(i,j,k,3) ** 2) / Q0(i,j,k,1)
        enddo
      enddo
    enddo
    call split(nx,ny,nz,nyp,nzp,Q0,Q)
    call set_cyclic_bc(nx,nyp,nzp,Q)
  end subroutine TaylorGreen
end module set_init
