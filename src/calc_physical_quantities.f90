module calc_physical_quantities
  use cudafor
  use mod_globals, only : gamma, R, dim => dimension
  use mod_constant, only : gamma_1
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
    pl   = Ql(dim+2)
    pr   = Qr(dim+2)
    Vl   = Ql(2:dim+1)
    Vr   = Qr(2:dim+1)
  end subroutine set_q

  subroutine calc_quantities_2D(nx,ny,Jacobian,QJ,rho,u,v,p)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: rho(nx,ny), u(nx,ny), v(nx,ny), p(nx,ny)
    integer i, j
    real(8) :: over_Q1
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        over_Q1  = 1.d0 / QJ(1,i,j)
        rho(i,j) = Jacobian(i,j) * QJ(1,i,j)
        u(i,j)   = QJ(2,i,j) * over_Q1
        v(i,j)   = QJ(3,i,j) * over_Q1
        p(i,j)   = gamma_1 * (Jacobian(i,j) * QJ(4,i,j)- 0.5d0 * rho(i,j) * (u(i,j)*u(i,j) + v(i,j)*v(i,j)))
    enddo;enddo
  end subroutine calc_quantities_2D

  subroutine calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    integer, intent(in), value                          :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny), device       :: Jacobian
    real(8), intent(in), dimension(5,nx,ny,nz), device  :: QJ ! Q / Jacobian
    real(8), intent(out), dimension(5,nx,ny,nz), device :: Q
    integer i, j, k
    real(8) :: over_Q1
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          over_Q1    = 1.d0 / QJ(1,i,j,k)
          Q(1,i,j,k) = Jacobian(i,j) * QJ(1,i,j,k)
          Q(2,i,j,k) = QJ(2,i,j,k) * over_Q1
          Q(3,i,j,k) = QJ(3,i,j,k) * over_Q1
          Q(4,i,j,k) = QJ(4,i,j,k) * over_Q1
          Q(5,i,j,k) = gamma_1 * (Jacobian(i,j) * QJ(5,i,j,k) - 0.5d0 * Q(1,i,j,k) * &
                      (Q(2,i,j,k)*Q(2,i,j,k) + Q(3,i,j,k)*Q(3,i,j,k) + Q(4,i,j,k)*Q(4,i,j,k)))
    enddo;enddo;enddo
  end subroutine calc_quantities_3D
  
  subroutine calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: Q(5,nx,ny,nz)
    real(8), intent(out), device :: T(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          over_Q1    = 1.d0 / QJ(1,i,j,k)
          Q(1,i,j,k) = Jacobian(i,j) * QJ(1,i,j,k)
          Q(2,i,j,k) = QJ(2,i,j,k) * over_Q1
          Q(3,i,j,k) = QJ(3,i,j,k) * over_Q1
          Q(4,i,j,k) = QJ(4,i,j,k) * over_Q1
          Q(5,i,j,k) = gamma_1 * (Jacobian(i,j) * QJ(5,i,j,k) - 0.5d0 * Q(1,i,j,k) * &
                      (Q(2,i,j,k)*Q(2,i,j,k) + Q(3,i,j,k)*Q(3,i,j,k) + Q(4,i,j,k)*Q(4,i,j,k)))
          T(i,j,k)   = Q(5,i,j,k) / (R * Q(1,i,j,k))
    enddo;enddo;enddo
  end subroutine calc_quantities_T_3D
end module calc_physical_quantities

