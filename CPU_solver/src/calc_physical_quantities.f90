module calc_physical_quantities
  use mod_globals, only : gamma, R
  implicit none
contains
  subroutine calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)
    integer, intent(in)                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz)   :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5) :: QJ ! Q / Jacobian
    real(8), intent(out), dimension(nx,ny,nz)  :: rho, u, v, w, p, T
    integer i, j, k
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Jacobian(i,j,k) * QJ(i,j,k,1)
          u(i,j,k) = QJ(i,j,k,2) / QJ(i,j,k,1)
          v(i,j,k) = QJ(i,j,k,3) / QJ(i,j,k,1)
          w(i,j,k) = QJ(i,j,k,4) / QJ(i,j,k,1)
          p(i,j,k) = (gamma - 1.d0) * (Jacobian(i,j,k) * QJ(i,j,k,5) &
                    - 0.5d0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
          T(i,j,k) = p(i,j,k) / (R * rho(i,j,k))
    enddo;enddo;enddo
  end subroutine calc_quantities
end module calc_physical_quantities

