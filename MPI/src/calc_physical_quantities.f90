module calc_physical_quantities
  use cudafor
  use mod_globals, only : gamma, R
  implicit none
contains
  attributes(device) subroutine set_q(Ql,Qr,rhol,rhor,pl,pr,Vl,Vr)
    real(8), intent(in), dimension(5) :: Ql, Qr
    real(8), intent(out) :: rhol, rhor, pl, pr
    real(8), intent(out), dimension(3) :: Vl, Vr
    rhol = Ql(1)
    rhor = Qr(1)
    pl = Ql(5)
    pr = Qr(5)
    Vl = Ql(2:4)
    Vr = Qr(2:4)
  end subroutine set_q

  subroutine calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)
    integer, intent(in), value                          :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny), device       :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(out), dimension(nx,ny,nz), device   :: rho, u, v, w, p, T
    integer i, j, k
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Jacobian(i,j) * QJ(i,j,k,1)
          u(i,j,k) = Jacobian(i,j) * QJ(i,j,k,2) / rho(i,j,k)
          v(i,j,k) = Jacobian(i,j) * QJ(i,j,k,3) / rho(i,j,k)
          w(i,j,k) = Jacobian(i,j) * QJ(i,j,k,4) / rho(i,j,k)
          p(i,j,k) = (gamma - 1.d0) * (Jacobian(i,j) * QJ(i,j,k,5) - 0.5d0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
          T(i,j,k) = p(i,j,k) / (R * rho(i,j,k))
    enddo;enddo;enddo
  end subroutine calc_quantities
end module calc_physical_quantities

