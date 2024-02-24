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
  
  subroutine TaylorGreen(na,nx,ny,nz,gamma,RHO,L,M,Q)
    use set_bc
    integer, intent(in) :: na, nx, ny, nz
    real(8), intent(in) :: gamma, RHO, L, M
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    integer i, j, k, offset
    real(8) pi
    real(8) x(nx-na), y(ny-na), z(nz-na)
    pi = 2.d0 * acos(0.d0)
    x = linspace(0.d0, 2.d0 * pi, nx-na)
    y = linspace(0.d0, 2.d0 * pi, ny-na)
    z = linspace(0.d0, 2.d0 * pi, nz-na)
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    offset = na/2
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          ! rho
          Q(i,j,k,1) = RHO
          ! rho u
          Q(i,j,k,2) = RHO * M * sin(x(i-offset)) * cos(y(j-offset)) * cos(z(k-offset))
          ! rho v
          Q(i,j,k,3) = - RHO * M * cos(x(i-offset)) * sin(y(j-offset)) * cos(z(k-offset))
          ! rho w
          Q(i,j,k,4) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q(i,j,k,5) = (1.d0/gamma+RHO*M**2*(cos(2.d0*x(i-offset)/L)+cos(2.d0*y(j-offset)/L))*(cos(2.d0*z(k-offset)/L)+2.d0)/16.d0)&
          / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2) ** 2 + Q(i,j,k,3) ** 2) / Q(i,j,k,1)
        enddo
      enddo
    enddo
    call set_cyclic_bc(na,nx,ny,nz,Q)
  end subroutine TaylorGreen
end module set_init
