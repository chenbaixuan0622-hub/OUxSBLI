! real(4) port of 3D_solver/src/calc_div.f90, 6th-order interior stencil only
! (mixed's NSTGV always uses VISC_ORDER=6, no fypp genericity needed here,
! matching the existing convention in calc_visc_cent_mp.f90.fypp/
! calc_visc_high_internal_mp.f90.fypp). Operates on the persistent real(4)
! mirrors u_r4/v_r4/w_r4 (calc_physical_quantities_mixed) and inv_dx_r4/
! inv_dy_r4/inv_dz_r4 instead of base's real(8) Q_2/Q_3/Q_4 and xix/etay/zetaz,
! producing real(4) ux/vy/wz consumed by calc_visc_high_internal_mp.f90.fypp.
module calc_div_mp
  use cudafor
  implicit none
  real(4), parameter :: one_120_mp = 1._4 / 120._4
contains
  !> du/dx at cell-center (6th-order, interior stencil)
  subroutine calc_div_ux_6_in_mp(nx, ny, nz, inv_dx_r4, u_r4, ux)
    integer, intent(in), value                 :: nx, ny, nz
    real(4), intent(in), device, contiguous    :: inv_dx_r4(nx-1)
    real(4), intent(in), device, contiguous    :: u_r4(nx,ny,nz)
    real(4), intent(out), device, contiguous   :: ux(nx,ny,nz)
    integer, parameter :: io_v = 2
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = 1, ny
        do i = io_v+2, nx-io_v-1
          ux(i,j,k) = (one_120_mp * (-u_r4(i-3,j,k) + u_r4(i+3,j,k)) &
                      + 0.075_4 * (u_r4(i-2,j,k) - u_r4(i+2,j,k)) &
                     + 0.375_4 * (-u_r4(i-1,j,k) + u_r4(i+1,j,k))) * (inv_dx_r4(i-1) + inv_dx_r4(i))
        enddo
      enddo
    enddo
  end subroutine calc_div_ux_6_in_mp

  !> dv/dy at cell-center (6th-order, interior stencil)
  subroutine calc_div_vy_6_in_mp(nx, ny, nz, inv_dy_r4, v_r4, vy)
    integer, intent(in), value                 :: nx, ny, nz
    real(4), intent(in), device, contiguous    :: inv_dy_r4(ny-1)
    real(4), intent(in), device, contiguous    :: v_r4(nx,ny,nz)
    real(4), intent(out), device, contiguous   :: vy(nx,ny,nz)
    integer, parameter :: io_v = 2
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = io_v+2, ny-io_v-1
        do i = 1, nx
          vy(i,j,k) = (one_120_mp * (-v_r4(i,j-3,k) + v_r4(i,j+3,k)) &
                      + 0.075_4 * (v_r4(i,j-2,k) - v_r4(i,j+2,k)) &
                     + 0.375_4 * (-v_r4(i,j-1,k) + v_r4(i,j+1,k))) * (inv_dy_r4(j-1) + inv_dy_r4(j))
        enddo
      enddo
    enddo
  end subroutine calc_div_vy_6_in_mp

  !> dw/dz at cell-center (6th-order, interior stencil)
  subroutine calc_div_wz_6_in_mp(nx, ny, nz, inv_dz_r4, w_r4, wz)
    integer, intent(in), value                 :: nx, ny, nz
    real(4), intent(in), device, contiguous    :: inv_dz_r4(nz-1)
    real(4), intent(in), device, contiguous    :: w_r4(nx,ny,nz)
    real(4), intent(out), device, contiguous   :: wz(nx,ny,nz)
    integer, parameter :: io_v = 2
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = io_v+2, nz-io_v-1
      do j = 1, ny
        do i = 1, nx
          wz(i,j,k) = (one_120_mp * (-w_r4(i,j,k-3) + w_r4(i,j,k+3)) &
                      + 0.075_4 * (w_r4(i,j,k-2) - w_r4(i,j,k+2)) &
                     + 0.375_4 * (-w_r4(i,j,k-1) + w_r4(i,j,k+1))) * (inv_dz_r4(k-1) + inv_dz_r4(k))
        enddo
      enddo
    enddo
  end subroutine calc_div_wz_6_in_mp
end module calc_div_mp
