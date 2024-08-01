module calc_physical_quantities
  use cudafor
  use mod_globals, only : gamma, R, dim => dimension
  implicit none
  interface calc_quantities
    module procedure calc_quantities_2D, calc_quantities_3D
  end interface
contains
  attributes(device) subroutine set_q(Ql,Qr,rhol,rhor,pl,pr,Vl,Vr)
    real(8), intent(in), dimension(dim+2) :: Ql, Qr
    real(8), intent(out) :: rhol, rhor, pl, pr
    real(8), intent(out), dimension(dim) :: Vl, Vr
    rhol = Ql(1)
    rhor = Qr(1)
    pl = Ql(dim+2)
    pr = Qr(dim+2)
    Vl = Ql(2:dim+1)
    Vr = Qr(2:dim+1)
  end subroutine set_q

  subroutine calc_quantities_2D(nx,ny,Q,rho,u,v,p,T)
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny), device  :: rho, u, v, p, T
    integer i, j
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        rho(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2) / rho(i,j)
        v(i,j) = Q(i,j,3) / rho(i,j)
        p(i,j) = (gamma - 1.d0) * (Q(i,j,4)- 0.5d0 * rho(i,j) * (u(i,j)**2 + v(i,j)**2))
        T(i,j) = p(i,j) / (R * rho(i,j))
    enddo;enddo
  end subroutine calc_quantities_2D

  subroutine calc_quantities_3D(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)
    integer, intent(in), value                          :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device    :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(out), dimension(nx,ny,nz), device   :: rho, u, v, w, p, T
    integer i, j, k
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Jacobian(i,j,k) * QJ(i,j,k,1)
          u(i,j,k) = QJ(i,j,k,2) / QJ(i,j,k,1)
          v(i,j,k) = QJ(i,j,k,3) / QJ(i,j,k,1)
          w(i,j,k) = QJ(i,j,k,4) / QJ(i,j,k,1)
          p(i,j,k) = (gamma - 1.d0) * (Jacobian(i,j,k) * QJ(i,j,k,5) - 0.5d0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
          T(i,j,k) = p(i,j,k) / (R * rho(i,j,k))
    enddo;enddo;enddo
  end subroutine calc_quantities_3D
end module calc_physical_quantities

