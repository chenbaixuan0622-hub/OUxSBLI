!> Curvilinear flux dispatcher (Euler only, 2nd-order).
!> Provides calc_EFG_curv with KEEP/SLAU/Hybrid dispatch via id_scheme kind.
module calc_flux_base_curv
  use cudafor
  use mod_globals, only : id_scheme, id_accuracy, id_slau, sp, blocks, threads, &
                          blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_keep_kernel_curv
  use calc_slau_kernel_curv
  use calc_hybrid_kernel_curv
  use calc_hybrid_curv
  implicit none
  private
  public init_sensor_curv, calc_EFG_curv

  real(sp), allocatable, device, save :: sensor(:,:,:)

  interface calc_conv_curv
    module procedure calc_conv_curv_keep, calc_conv_curv_slau, calc_conv_curv_hybrid
  end interface calc_conv_curv

contains
  subroutine init_sensor_curv(nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    allocate(sensor(nx,ny,nz))
  end subroutine init_sensor_curv


  !> KEEP scheme dispatch (id_scheme kind = integer(2)).
  subroutine calc_conv_curv_keep(id_scheme, nx, ny, nz, dz, &
      n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Q, T, E, F, G)
    integer(2), intent(in), value           :: id_scheme
    integer, intent(in), value              :: nx, ny, nz
    real(8), intent(in), value              :: dz
    real(8), intent(in), device, contiguous :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8), intent(in), device, contiguous :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8), intent(in), device, contiguous :: xi_x(nx,ny), xi_y(nx,ny)
    real(8), intent(in), device, contiguous :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous :: T(nx,ny,nz)
    real(8), intent(out),device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out),device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out),device, contiguous :: G(5,nx-2,ny-2,nz-1)
    call calc_keep_xi_curv<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q, T, E)
    call calc_keep_eta_curv<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q, T, F)
    call calc_keep_z_curv<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, T, G)
  end subroutine calc_conv_curv_keep


  !> SLAU scheme dispatch (id_scheme kind = real(2)).
  subroutine calc_conv_curv_slau(id_scheme, nx, ny, nz, dz, &
      n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Q, T, E, F, G)
    real(2), intent(in), value              :: id_scheme
    integer, intent(in), value              :: nx, ny, nz
    real(8), intent(in), value              :: dz
    real(8), intent(in), device, contiguous :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8), intent(in), device, contiguous :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8), intent(in), device, contiguous :: xi_x(nx,ny), xi_y(nx,ny)
    real(8), intent(in), device, contiguous :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous :: T(nx,ny,nz)
    real(8), intent(out),device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out),device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out),device, contiguous :: G(5,nx-2,ny-2,nz-1)
    call calc_Ducros_curv<<<blocks,threads>>>(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, Q, sensor)
    call calc_slau_xi_curv<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q, sensor, E)
    call calc_slau_eta_curv<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q, sensor, F)
    call calc_slau_z_curv<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
  end subroutine calc_conv_curv_slau


  !> Hybrid (KEEP+SLAU) scheme dispatch (id_scheme kind = real(8)).
  subroutine calc_conv_curv_hybrid(id_scheme, nx, ny, nz, dz, &
      n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Q, T, E, F, G)
    real(8), intent(in), value              :: id_scheme
    integer, intent(in), value              :: nx, ny, nz
    real(8), intent(in), value              :: dz
    real(8), intent(in), device, contiguous :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8), intent(in), device, contiguous :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8), intent(in), device, contiguous :: xi_x(nx,ny), xi_y(nx,ny)
    real(8), intent(in), device, contiguous :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous :: T(nx,ny,nz)
    real(8), intent(out),device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out),device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out),device, contiguous :: G(5,nx-2,ny-2,nz-1)
    call calc_Ducros_curv<<<blocks,threads>>>(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, Q, sensor)
    call calc_hybrid_xi_curv<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q, T, sensor, E)
    call calc_hybrid_eta_curv<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q, T, sensor, F)
    call calc_hybrid_z_curv<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, T, sensor, G)
  end subroutine calc_conv_curv_hybrid


  !> Euler curvilinear flux (id_visc kind=2).
  !> Reconstructs primitive Q from QJ, then dispatches to KEEP/SLAU/Hybrid.
  subroutine calc_EFG_curv(id_visc, nx, ny, nz, dz, n_xi_x, n_xi_y, n_eta_x, n_eta_y, &
      xi_x, xi_y, eta_x, eta_y, Jacobian, QJ, Q, T, E, F, G)
    integer(2), intent(in), value            :: id_visc
    integer, intent(in), value               :: nx, ny, nz
    real(8), intent(in), value               :: dz
    real(8), intent(in), device, contiguous  :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8), intent(in), device, contiguous  :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8), intent(in), device, contiguous  :: xi_x(nx,ny), xi_y(nx,ny)
    real(8), intent(in), device, contiguous  :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(in), device, contiguous  :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous  :: QJ(nx,5,ny,nz)
    real(8), intent(out), device, contiguous :: Q(nx,5,ny,nz)
    real(8), intent(out), device, contiguous :: T(nx,ny,nz)
    real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T)
    call calc_conv_curv(id_scheme, nx, ny, nz, dz, &
        n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, &
        Q, T, E, F, G)
  end subroutine calc_EFG_curv
end module calc_flux_base_curv
