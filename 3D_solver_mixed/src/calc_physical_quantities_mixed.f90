! MIXED_FORK-only copy of src/calc_physical_quantities.f90's calc_quantities_T_3D,
! specialized for NSTGV (VISC='NS'). Differs from the baseline in two ways:
!   1. `mu` is real(4) (Sutherland's law is still evaluated in real(8)
!      internally; only the final store narrows to real4), matching the
!      MIXED_FORK-gated real(4) mu allocated by 3D_solver/src/preprocess.f90.fypp.
!   2. Persistent real(4) module-level mirrors u_r4/v_r4/w_r4/T_r4 are written
!      from the already-local u/v/w/temp values in the same loop iteration
!      that computes them -- zero extra device memory reads. These are
!      consumed directly by the real(4)-internal viscous kernels
!      (calc_visc_high_internal_mp.f90.fypp / load_smem_visc_cent_mp.f90.fypp),
!      eliminating the redundant per-thread real(8)->real(4) casts those
!      kernels used to perform on Q/T themselves.
! Real(8) Q_1..Q_5 and T are still written exactly as calc_quantities_T_3D
! does today -- the real(8) SLAU convective path still consumes them.
module calc_physical_quantities_mixed
  use cudafor
  use mod_globals, only : gamma, R
  use mod_constant, only : gamma_1, mu0_T0_S_over_T0_2_3
  implicit none
  real(4), allocatable, device, save, public :: u_r4(:,:,:), v_r4(:,:,:), w_r4(:,:,:), T_r4(:,:,:)
contains
  subroutine init_physical_quantities_mixed(nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    allocate(u_r4(nx,ny,nz), v_r4(nx,ny,nz), w_r4(nx,ny,nz), T_r4(nx,ny,nz))
  end subroutine init_physical_quantities_mixed


  subroutine calc_quantities_T_3D_mixed(nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, &
                                         Q_1, Q_2, Q_3, Q_4, Q_5, T, mu, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous    :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device, contiguous    :: QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: T(nx,ny,nz)
    real(4), intent(inout), device, contiguous :: mu(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p, temp
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = 1, nx
          over_Q1  = 1.d0 / QJ_1(i,j,k)
          rho      = Jacobian(i,j) * QJ_1(i,j,k)
          u        = QJ_2(i,j,k) * over_Q1
          v        = QJ_3(i,j,k) * over_Q1
          w        = QJ_4(i,j,k) * over_Q1
          p        = gamma_1 * (Jacobian(i,j) * QJ_5(i,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q_1(i,j,k) = rho
          Q_2(i,j,k) = u
          Q_3(i,j,k) = v
          Q_4(i,j,k) = w
          Q_5(i,j,k) = p
          temp       = p / (R * rho)
          T(i,j,k)   = temp
          mu(i,j,k)  = real(mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp)), 4)
          u_r4(i,j,k) = real(u, 4)
          v_r4(i,j,k) = real(v, 4)
          w_r4(i,j,k) = real(w, 4)
          T_r4(i,j,k) = real(temp, 4)
    enddo;enddo;enddo
  end subroutine calc_quantities_T_3D_mixed
end module calc_physical_quantities_mixed
