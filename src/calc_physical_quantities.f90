module calc_physical_quantities
  use cudafor
  use mod_globals, only : gamma, R
  use mod_constant, only : gamma_1, mu0_T0_S_over_T0_2_3
  implicit none
contains
  subroutine calc_quantities_2D(nx, ny, Jacobian, QJ, Q, T)
    integer, intent(in), value               :: nx, ny
    real(8), intent(in), device, contiguous  :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous  :: QJ(nx,4,ny) ! Q / Jacobian
    real(8), intent(out), device, contiguous :: Q(nx,4,ny)
    real(8), intent(out), device, contiguous :: T(nx,ny)
    integer i, j
    real(8) :: over_Q1, rho, u, v, p
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny
      do i = 1, nx
        over_Q1  = 1.d0 / QJ(i,1,j)
        rho      = Jacobian(i,j) * QJ(i,1,j)
        u        = QJ(i,2,j) * over_Q1
        v        = QJ(i,3,j) * over_Q1
        p        = gamma_1 * (Jacobian(i,j) * QJ(i,4,j) - 0.5d0 * rho * (u*u + v*v))
        Q(i,1,j) = rho
        Q(i,2,j) = u
        Q(i,3,j) = v
        Q(i,4,j) = p
        T(i,j)   = p / (R * rho)
    enddo;enddo
  end subroutine calc_quantities_2D
  

  subroutine calc_quantities_T_2D(nx, ny, Jacobian, QJ, Q, T, mu)
    integer, intent(in), value               :: nx, ny
    real(8), intent(in), device, contiguous  :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous  :: QJ(nx,4,ny) ! Q / Jacobian
    real(8), intent(out), device, contiguous :: Q(nx,4,ny)
    real(8), intent(out), device, contiguous :: T(nx,ny)
    real(8), intent(out), device, contiguous :: mu(nx,ny)
    integer i, j
    real(8) :: over_Q1, rho, u, v, p, temp
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny
      do i = 1, nx
        over_Q1  = 1.d0 / QJ(i,1,j)
        rho      = Jacobian(i,j) * QJ(i,1,j)
        u        = QJ(i,2,j) * over_Q1
        v        = QJ(i,3,j) * over_Q1
        p        = gamma_1 * (Jacobian(i,j) * QJ(i,4,j) - 0.5d0 * rho * (u*u + v*v))
        Q(i,1,j) = rho
        Q(i,2,j) = u
        Q(i,3,j) = v
        Q(i,4,j) = p
        temp     = p / (R * rho)
        T(i,j)   = temp
        mu(i,j)  = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
    enddo;enddo
  end subroutine calc_quantities_T_2D


  subroutine calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous    :: QJ(nx,5,ny,nz) ! Q / Jacobian
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: T(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = 1, nx
          over_Q1    = 1.d0 / QJ(i,1,j,k)
          rho        = Jacobian(i,j) * QJ(i,1,j,k)
          u          = QJ(i,2,j,k) * over_Q1
          v          = QJ(i,3,j,k) * over_Q1
          w          = QJ(i,4,j,k) * over_Q1
          p          = gamma_1 * (Jacobian(i,j) * QJ(i,5,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q(i,1,j,k) = rho
          Q(i,2,j,k) = u
          Q(i,3,j,k) = v
          Q(i,4,j,k) = w
          Q(i,5,j,k) = p
          T(i,j,k)   = p / (R * rho)
    enddo;enddo;enddo
  end subroutine calc_quantities_3D


  subroutine calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous    :: QJ(nx,5,ny,nz) ! Q / Jacobian
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mu(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p, temp
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = 1, nx
          over_Q1    = 1.d0 / QJ(i,1,j,k)
          rho        = Jacobian(i,j) * QJ(i,1,j,k)
          u          = QJ(i,2,j,k) * over_Q1
          v          = QJ(i,3,j,k) * over_Q1
          w          = QJ(i,4,j,k) * over_Q1
          p          = gamma_1 * (Jacobian(i,j) * QJ(i,5,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q(i,1,j,k) = rho
          Q(i,2,j,k) = u
          Q(i,3,j,k) = v
          Q(i,4,j,k) = w
          Q(i,5,j,k) = p
          temp       = p / (R * rho)
          T(i,j,k)   = temp
          mu(i,j,k)  = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
    enddo;enddo;enddo
  end subroutine calc_quantities_T_3D
end module calc_physical_quantities

