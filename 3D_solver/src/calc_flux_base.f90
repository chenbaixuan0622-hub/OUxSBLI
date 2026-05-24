!> Module for computing convective and viscous fluxes
!> Dispatches to different numerical schemes (KEEP, SLAU, Roe, Hybrid)
!> Groups all GPU kernel calls for computing E, F, G flux components
module calc_flux_base
  use mod_globals, only : id_accuracy, id_bc_x, id_bc_y, id_bc_z, sp, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_physical_quantities
  use calc_hybrid
  use calc_keep_kernel
  use calc_keep_kernel_internal
  use calc_slau_kernel
  use calc_slau_kernel_internal
  use calc_roe_kernel
  use calc_roe_kernel_internal
  use calc_hybrid_kernel
  use calc_hybrid_kernel_internal
  use calc_visc2
  use calc_visc4
  use calc_visc4_internal
  use calc_visc4_les_internal
  use calc_les
  use set
  implicit none
  private
  public calc_EFG, calc_EFG_halo, init_sensor
  !> Ghost-cell width: io+1 where io = kind(id_accuracy)/3 (0/1/2 for 2nd/4th/6th order)
  integer, parameter :: io_fb      = kind(id_accuracy) / 3
  integer, parameter :: overlap_fb = io_fb + 1
  !> Persistent Ducros shock sensor on device — allocated once in init_sensor,
  !> reused across all calls to calc_conv_slau / calc_conv_roe / calc_conv_hybrid.
  !> Eliminates repeated device heap alloc/free on every RK stage.
  real(sp), allocatable, device, save :: sensor(:,:,:)
  interface calc_conv
    module procedure calc_conv_keep, calc_conv_slau, calc_conv_roe, calc_conv_hybrid
  end interface calc_conv

  !> Convective E and F only (no G); used in zdec interior phase.
  interface calc_conv_EF
    module procedure calc_conv_EF_keep, calc_conv_EF_slau, calc_conv_EF_roe, calc_conv_EF_hybrid
  end interface calc_conv_EF

  !> Convective G for an explicit [k_lo, k_hi] range; scheme dispatched via id_scheme kind.
  !> Replaces the former calc_conv_G_interior_koff / calc_conv_G_halo_koff pair.
  interface calc_conv_G_koff
    module procedure calc_conv_G_koff_keep, calc_conv_G_koff_slau, calc_conv_G_koff_hybrid
  end interface calc_conv_G_koff

  !> Full domain or zdec interior flux dispatch (Euler/NS/LES via id_visc kind).
  !> k_lo_G=1, k_hi_G=nz-1: full domain; other values: zdec interior restricted G.
  interface calc_EFG
    module procedure calc_EFG_Euler, calc_EFG_visc, calc_EFG_LES
  end interface calc_EFG

  !> Halo-z phase: ghost Q refresh + G halo (lo+hi); called after MPI ghost exchange completes.
  interface calc_EFG_halo
    module procedure calc_EFG_halo_Euler, calc_EFG_halo_visc, calc_EFG_halo_LES
  end interface calc_EFG_halo
contains
  !> Allocate the persistent Ducros sensor array on the device.
  !> Call once from preprocess.f90 before the time loop.
  subroutine init_sensor(nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    allocate(sensor(nx,ny,nz))
  end subroutine init_sensor

  !> Compute convective fluxes using KEEP (energy-preserving) scheme
  !> High-order minimal dissipation scheme for smooth flow regions
  !> Algorithm: F = (H·u) where H = enthalpy, using flux reconstruction via divergence forms
  !> Provides 4th-5th order accuracy by minimizing dispersive errors in smooth regions
  subroutine calc_conv_keep(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    integer(2), intent(in), value            :: id_scheme           !< ID for scheme: int 2 means KEEP
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !< inverse grid spacing y (1/dy)
    real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !< inverse grid spacing z (1/dz)
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)       !< conservative variables Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous  :: T(nx,ny,nz)         !< temperature field
    real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< convective flux in x direction
    real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< convective flux in y direction
    real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !< convective flux in z direction
    if (id_bc_x) then
      call calc_keep_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, T, E)
    else
      call calc_keep_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, T, E)
    endif
    if (id_bc_y) then
      call calc_keep_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, T, F)
    else
      call calc_keep_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, T, F)
    endif
    if (id_bc_z) then
      call calc_keep_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, T, G)
    else
      call calc_keep_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, T, G)
    endif
  end subroutine calc_conv_keep


  !> Compute convective fluxes using SLAU (Simple Low-dissipation Roe-based Upwind) scheme
  !> Low-dissipation scheme with shock-capturing capability via Ducros sensor
  !> Algorithm: F = (F_L + F_R)/2 + |A|(Q_L - Q_R)/2 where A is weighted Jacobian
  !> Dissipation modulated by Ducros sensor: f_d controls blend ratio
  subroutine calc_conv_slau(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(2), intent(in), value               :: id_scheme           !< ID for scheme: real 2 means SLAU
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !< inverse grid spacing y (1/dy)
    real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !< inverse grid spacing z (1/dz)
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous  :: T(nx,ny,nz)         !< temperature field
    real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< convective flux in x direction
    real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< convective flux in y direction
    real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !< convective flux in z direction
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    if (id_bc_x) then
      call calc_slau_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    else
      call calc_slau_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, sensor, E)
    endif
    if (id_bc_y) then
      call calc_slau_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    else
      call calc_slau_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, sensor, F)
    endif
    if (id_bc_z) then
      call calc_slau_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    else
      call calc_slau_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, sensor, G)
    endif
  end subroutine calc_conv_slau


  !> Compute convective fluxes using Roe approximate Riemann solver
  !> Classic approximate Riemann solver with wave decomposition for flux splitting
  !> Algorithm: F = (F_L + F_R)/2 - (1/2)Σ|λ_i|*(p_i·r_i) wave reconstruction
  !> Where λ_i are Roe eigenvalues, p_i are wave strengths, r_i are eigenvectors
  !> Entropy fix via sensor prevents expansion shocks at sonic points
  subroutine calc_conv_roe(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(4), intent(in), value               :: id_scheme           !< ID for scheme: real 4 means Roe
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !< inverse grid spacing y (1/dy)
    real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !< inverse grid spacing z (1/dz)
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous  :: T(nx,ny,nz)         !< temperature field
    real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< convective flux in x direction
    real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< convective flux in y direction
    real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !< convective flux in z direction
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    if (id_bc_x) then
      call calc_roe_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    else
      call calc_roe_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, sensor, E)
    endif
    if (id_bc_y) then
      call calc_roe_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    else
      call calc_roe_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, sensor, F)
    endif
    if (id_bc_z) then
      call calc_roe_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    else
      call calc_roe_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, sensor, G)
    endif
  end subroutine calc_conv_roe


  !> Compute convective fluxes using hybrid KEEP/SLAU scheme  
  !> Automatically blends between KEEP (smooth regions) and SLAU (shock regions) seamlessly
  !> Blending formula: F_hybrid = (1-f_d)·F_keep + f_d·F_slau where f_d ∈ [0,1]
  !> Preserves vortex structures (f_d≈0) and captures shocks accurately (f_d≈1)
  !> Ducros shock sensor: f_d = (∇·u)²/[(∇·u)² + (∇×u)² + ε]
  subroutine calc_conv_hybrid(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(8), intent(in), value               :: id_scheme           !< ID for scheme: real 8 means Hybrid
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous  :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    if (id_bc_x) then
      call calc_hybrid_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    else
      call calc_hybrid_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, T, sensor, E)
    endif
    if (id_bc_y) then
      call calc_hybrid_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    else
      call calc_hybrid_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, T, sensor, F)
    endif
    if (id_bc_z) then
      call calc_hybrid_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, T, sensor, G)
    else
      call calc_hybrid_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, T, sensor, G)
    endif
  end subroutine calc_conv_hybrid

  ! ---------------------------------------------------------------------------
  ! calc_conv_G_halo_koff: private scheme-dispatched wrappers
  ! Each subroutine launches the correct z_in_koff kernel for lo and hi halos.
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! calc_conv_EF: E and F only; G is omitted so the interior phase can later
  ! call calc_conv_G_interior_koff for a restricted k-range.
  ! ---------------------------------------------------------------------------

  subroutine calc_conv_EF_keep(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F)
    use mod_globals, only : id_accuracy
    integer(2), intent(in), value              :: id_scheme
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device, contiguous   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device, contiguous   :: F(5,nx-2,ny-1,nz-2)
    if (id_bc_x) then
      call calc_keep_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, T, E)
    else
      call calc_keep_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, T, E)
    endif
    if (id_bc_y) then
      call calc_keep_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, T, F)
    else
      call calc_keep_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, T, F)
    endif
  end subroutine calc_conv_EF_keep


  subroutine calc_conv_EF_slau(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F)
    use mod_globals, only : id_accuracy
    real(2), intent(in), value                 :: id_scheme
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device, contiguous   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device, contiguous   :: F(5,nx-2,ny-1,nz-2)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    if (id_bc_x) then
      call calc_slau_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    else
      call calc_slau_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, sensor, E)
    endif
    if (id_bc_y) then
      call calc_slau_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    else
      call calc_slau_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, sensor, F)
    endif
  end subroutine calc_conv_EF_slau


  subroutine calc_conv_EF_roe(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F)
    use mod_globals, only : id_accuracy
    real(4), intent(in), value                 :: id_scheme
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device, contiguous   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device, contiguous   :: F(5,nx-2,ny-1,nz-2)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    if (id_bc_x) then
      call calc_roe_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    else
      call calc_roe_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, sensor, E)
    endif
    if (id_bc_y) then
      call calc_roe_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    else
      call calc_roe_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, sensor, F)
    endif
  end subroutine calc_conv_EF_roe


  subroutine calc_conv_EF_hybrid(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F)
    use mod_globals, only : id_accuracy
    real(8), intent(in), value                 :: id_scheme
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device, contiguous   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device, contiguous   :: F(5,nx-2,ny-1,nz-2)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    if (id_bc_x) then
      call calc_hybrid_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    else
      call calc_hybrid_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, T, sensor, E)
    endif
    if (id_bc_y) then
      call calc_hybrid_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    else
      call calc_hybrid_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, T, sensor, F)
    endif
  end subroutine calc_conv_EF_hybrid


  subroutine calc_conv_G_koff_keep(id_scheme, nx, ny, nz, Q, T, G, k_lo, k_hi)
    integer(2), intent(in), value              :: id_scheme
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer :: n_G
    type(dim3) :: bG
    n_G = k_hi - k_lo + 1
    bG = dim3(blocksG%x, blocksG%y, (n_G + threadsG%z - 1) / threadsG%z)
    call calc_keep_z_in_koff<<<bG,threadsG>>>(nx, ny, nz, Q, T, G, k_lo, k_hi)
  end subroutine calc_conv_G_koff_keep


  subroutine calc_conv_G_koff_slau(id_scheme, nx, ny, nz, Q, T, G, k_lo, k_hi)
    real(2), intent(in), value                 :: id_scheme
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer :: n_G
    type(dim3) :: bG
    n_G = k_hi - k_lo + 1
    bG = dim3(blocksG%x, blocksG%y, (n_G + threadsG%z - 1) / threadsG%z)
    call calc_slau_z_in_koff<<<bG,threadsG>>>(nx, ny, nz, Q, sensor, G, k_lo, k_hi)
  end subroutine calc_conv_G_koff_slau


  subroutine calc_conv_G_koff_hybrid(id_scheme, nx, ny, nz, Q, T, G, k_lo, k_hi)
    real(8), intent(in), value                 :: id_scheme
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer :: n_G
    type(dim3) :: bG
    n_G = k_hi - k_lo + 1
    bG = dim3(blocksG%x, blocksG%y, (n_G + threadsG%z - 1) / threadsG%z)
    call calc_hybrid_z_in_koff<<<bG,threadsG>>>(nx, ny, nz, Q, T, sensor, G, k_lo, k_hi)
  end subroutine calc_conv_G_koff_hybrid


  subroutine calc_EFG_Euler(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G, k_lo_G, k_hi_G)
    use mod_globals, only : id_scheme
    integer(2), intent(in), value              :: id_visc
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny), QJ(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mu(1,1,1), mut(1,1,1), qc2(1,1,1)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, intent(in)                        :: k_lo_G, k_hi_G
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T, 1, nz)
    if (k_lo_G == 1 .and. k_hi_G == nz - 1) then
      call calc_conv(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    else
      call calc_conv_EF(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F)
      call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, k_lo_G, k_hi_G)
    end if
  end subroutine calc_EFG_Euler


  subroutine calc_EFG_visc(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G, k_lo_G, k_hi_G)
    use mod_globals, only : id_scheme
    integer(4), intent(in), value              :: id_visc
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny), QJ(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mut(1,1,1), qc2(1,1,1)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, intent(in)                        :: k_lo_G, k_hi_G
    type(dim3) :: blocksGv_int, blocksEv_int, blocksFv_int
    integer :: n_int_Gv, n_int_Ev
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu, 1, nz)
    if (k_lo_G == 1 .and. k_hi_G == nz - 1) then
      call calc_conv(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
      if (id_visc == 2) then
        if (id_bc_x == .false. .and. kind(id_accuracy) == 8) then
          call calc_Ev4_in<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
        else
          call calc_Ev4<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
        endif
        if (id_bc_y == .false. .and. kind(id_accuracy) == 8) then
          call calc_Fv4_in<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
        else
          call calc_Fv4<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
        endif
        if (id_bc_z == .false. .and. kind(id_accuracy) == 8) then
          call calc_Gv4_in<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
        else
          call calc_Gv4<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
        endif
      else
        call calc_Ev2<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
        call calc_Fv2<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
        call calc_Gv2<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
      endif
    else
      call calc_conv_EF(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F)
      call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, k_lo_G, k_hi_G)
      n_int_Ev = nz - 2*overlap_fb - 2
      n_int_Gv = nz - 2*overlap_fb - 1
      blocksEv_int = dim3(blocksEv%x, blocksEv%y, (n_int_Ev + threadsEv%z - 1) / threadsEv%z)
      blocksFv_int = dim3(blocksFv%x, blocksFv%y, (n_int_Ev + threadsFv%z - 1) / threadsFv%z)
      blocksGv_int = dim3(blocksGv%x, blocksGv%y, (n_int_Gv + threadsGv%z - 1) / threadsGv%z)
      if (kind(id_accuracy) >= 4) then
        if (.not. id_bc_x .and. kind(id_accuracy) == 8) then
          call calc_Ev4_in_koff<<<blocksEv_int,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
               overlap_fb+2, nz-overlap_fb-1)
        else
          call calc_Ev4_koff<<<blocksEv_int,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
               overlap_fb+2, nz-overlap_fb-1)
        end if
        if (.not. id_bc_y .and. kind(id_accuracy) == 8) then
          call calc_Fv4_in_koff<<<blocksFv_int,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
               overlap_fb+2, nz-overlap_fb-1)
        else
          call calc_Fv4_koff<<<blocksFv_int,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
               overlap_fb+2, nz-overlap_fb-1)
        end if
        call calc_Gv4_koff<<<blocksGv_int,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G, &
             overlap_fb+1, nz-overlap_fb-1)
      else
        call calc_Ev2_koff<<<blocksEv_int,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
             overlap_fb+2, nz-overlap_fb-1)
        call calc_Fv2_koff<<<blocksFv_int,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
             overlap_fb+2, nz-overlap_fb-1)
        call calc_Gv2_koff<<<blocksGv_int,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G, &
             overlap_fb+1, nz-overlap_fb-1)
      end if
    end if
  end subroutine calc_EFG_visc

 
  subroutine calc_EFG_LES(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G, k_lo_G, k_hi_G)
    use mod_globals, only : id_scheme
    integer(8), intent(in), value              :: id_visc
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny), QJ(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, intent(in)                        :: k_lo_G, k_hi_G
    type(dim3) :: blocksGv_int, blocksEv_int, blocksFv_int
    integer :: n_int_Gv, n_int_Ev
    mut = 0.d0
    qc2 = 0.d0
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu, 1, nz)
    if (k_lo_G == 1 .and. k_hi_G == nz - 1) then
      call calc_conv(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
      call calc_mut<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, mut, qc2)
      call set_bc_mut(nx, ny, nz, mut, qc2)
      if (id_visc == 2) then
        if (id_bc_x == .false. .and. kind(id_accuracy) == 8) then
          call calc_Ev_LES4_in<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E)
        else
          call calc_Ev_LES4<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E)
        endif
        if (id_bc_y == .false. .and. kind(id_accuracy) == 8) then
          call calc_Fv_LES4_in<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F)
        else
          call calc_Fv_LES4<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F)
        endif
        if (id_bc_z == .false. .and. kind(id_accuracy) == 8) then
          call calc_Gv_LES4_in<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
        else
          call calc_Gv_LES4<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
        endif
      else
        call calc_Ev_LES2<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E)
        call calc_Fv_LES2<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F)
        call calc_Gv_LES2<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
      endif
    else
      call calc_conv_EF(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F)
      call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, k_lo_G, k_hi_G)
      call calc_mut<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, mut, qc2)
      call set_bc_mut(nx, ny, nz, mut, qc2)
      n_int_Ev = nz - 2*overlap_fb - 2
      n_int_Gv = nz - 2*overlap_fb - 1
      blocksEv_int = dim3(blocksEv%x, blocksEv%y, (n_int_Ev + threadsEv%z - 1) / threadsEv%z)
      blocksFv_int = dim3(blocksFv%x, blocksFv%y, (n_int_Ev + threadsFv%z - 1) / threadsFv%z)
      blocksGv_int = dim3(blocksGv%x, blocksGv%y, (n_int_Gv + threadsGv%z - 1) / threadsGv%z)
      if (kind(id_accuracy) >= 4) then
        if (.not. id_bc_x .and. kind(id_accuracy) == 8) then
          call calc_Ev_LES4_in_koff<<<blocksEv_int,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
               overlap_fb+2, nz-overlap_fb-1)
        else
          call calc_Ev_LES4_koff<<<blocksEv_int,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
               overlap_fb+2, nz-overlap_fb-1)
        end if
        if (.not. id_bc_y .and. kind(id_accuracy) == 8) then
          call calc_Fv_LES4_in_koff<<<blocksFv_int,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
               overlap_fb+2, nz-overlap_fb-1)
        else
          call calc_Fv_LES4_koff<<<blocksFv_int,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
               overlap_fb+2, nz-overlap_fb-1)
        end if
        call calc_Gv_LES4_koff<<<blocksGv_int,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G, &
             overlap_fb+1, nz-overlap_fb-1)
      else
        call calc_Ev_LES2_koff<<<blocksEv_int,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
             overlap_fb+2, nz-overlap_fb-1)
        call calc_Fv_LES2_koff<<<blocksFv_int,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
             overlap_fb+2, nz-overlap_fb-1)
        call calc_Gv_LES2_koff<<<blocksGv_int,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G, &
             overlap_fb+1, nz-overlap_fb-1)
      end if
    end if
  end subroutine calc_EFG_LES

  ! ---------------------------------------------------------------------------
  ! Halo phase: ghost Q refresh + G halo (lo+hi); called after MPI ghost exchange.
  ! ---------------------------------------------------------------------------

  subroutine calc_EFG_halo_Euler(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, &
                                   Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)
    use mod_globals, only : id_scheme
    integer(2), intent(in), value              :: id_visc
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny), QJ(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mu(1,1,1), mut(1,1,1), qc2(1,1,1)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T, 1, overlap_fb)
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T, nz-overlap_fb+1, nz)
    call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, overlap_fb, 2*overlap_fb-1)
    call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, nz-2*overlap_fb+1, nz-overlap_fb)
  end subroutine calc_EFG_halo_Euler


  subroutine calc_EFG_halo_visc(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, &
                                  Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)
    use mod_globals, only : id_scheme
    integer(4), intent(in), value              :: id_visc
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny), QJ(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mut(1,1,1), qc2(1,1,1)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    type(dim3) :: blocksGv_halo, blocksEv_halo, blocksFv_halo
    integer :: n_halo
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu, 1, overlap_fb)
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu, nz-overlap_fb+1, nz)
    n_halo = overlap_fb
    blocksEv_halo = dim3(blocksEv%x, blocksEv%y, (n_halo + threadsEv%z - 1) / threadsEv%z)
    blocksFv_halo = dim3(blocksFv%x, blocksFv%y, (n_halo + threadsFv%z - 1) / threadsFv%z)
    blocksGv_halo = dim3(blocksGv%x, blocksGv%y, (n_halo + threadsGv%z - 1) / threadsGv%z)
    call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, overlap_fb, 2*overlap_fb-1)
    call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, nz-2*overlap_fb+1, nz-overlap_fb)
    if (kind(id_accuracy) >= 4) then
      if (.not. id_bc_x .and. kind(id_accuracy) == 8) then
        call calc_Ev4_in_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
             2, overlap_fb+1)
        call calc_Ev4_in_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
             nz-overlap_fb, nz-1)
      else
        call calc_Ev4_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
             2, overlap_fb+1)
        call calc_Ev4_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
             nz-overlap_fb, nz-1)
      end if
      if (.not. id_bc_y .and. kind(id_accuracy) == 8) then
        call calc_Fv4_in_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
             2, overlap_fb+1)
        call calc_Fv4_in_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
             nz-overlap_fb, nz-1)
      else
        call calc_Fv4_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
             2, overlap_fb+1)
        call calc_Fv4_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
             nz-overlap_fb, nz-1)
      end if
      call calc_Gv4_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G, &
           1, overlap_fb)
      call calc_Gv4_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G, &
           nz-overlap_fb, nz-1)
    else
      call calc_Ev2_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
           2, overlap_fb+1)
      call calc_Ev2_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E, &
           nz-overlap_fb, nz-1)
      call calc_Fv2_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
           2, overlap_fb+1)
      call calc_Fv2_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F, &
           nz-overlap_fb, nz-1)
      call calc_Gv2_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G, &
           1, overlap_fb)
      call calc_Gv2_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G, &
           nz-overlap_fb, nz-1)
    end if
  end subroutine calc_EFG_halo_visc


  subroutine calc_EFG_halo_LES(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, &
                                 Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)
    use mod_globals, only : id_scheme
    integer(8), intent(in), value              :: id_visc
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: inv_dx(nx-1), inv_dy(ny-1), inv_dz(nz-1)
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny), QJ(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    type(dim3) :: blocksGv_halo, blocksEv_halo, blocksFv_halo
    integer :: n_halo
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu, 1, overlap_fb)
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu, nz-overlap_fb+1, nz)
    n_halo = overlap_fb
    blocksEv_halo = dim3(blocksEv%x, blocksEv%y, (n_halo + threadsEv%z - 1) / threadsEv%z)
    blocksFv_halo = dim3(blocksFv%x, blocksFv%y, (n_halo + threadsFv%z - 1) / threadsFv%z)
    blocksGv_halo = dim3(blocksGv%x, blocksGv%y, (n_halo + threadsGv%z - 1) / threadsGv%z)
    call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, overlap_fb, 2*overlap_fb-1)
    call calc_conv_G_koff(id_scheme, nx, ny, nz, Q, T, G, nz-2*overlap_fb+1, nz-overlap_fb)
    if (kind(id_accuracy) >= 4) then
      if (.not. id_bc_x .and. kind(id_accuracy) == 8) then
        call calc_Ev_LES4_in_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
             2, overlap_fb+1)
        call calc_Ev_LES4_in_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
             nz-overlap_fb, nz-1)
      else
        call calc_Ev_LES4_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
             2, overlap_fb+1)
        call calc_Ev_LES4_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
             nz-overlap_fb, nz-1)
      end if
      if (.not. id_bc_y .and. kind(id_accuracy) == 8) then
        call calc_Fv_LES4_in_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
             2, overlap_fb+1)
        call calc_Fv_LES4_in_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
             nz-overlap_fb, nz-1)
      else
        call calc_Fv_LES4_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
             2, overlap_fb+1)
        call calc_Fv_LES4_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
             nz-overlap_fb, nz-1)
      end if
      call calc_Gv_LES4_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G, &
           1, overlap_fb)
      call calc_Gv_LES4_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G, &
           nz-overlap_fb, nz-1)
    else
      call calc_Ev_LES2_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
           2, overlap_fb+1)
      call calc_Ev_LES2_koff<<<blocksEv_halo,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E, &
           nz-overlap_fb, nz-1)
      call calc_Fv_LES2_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
           2, overlap_fb+1)
      call calc_Fv_LES2_koff<<<blocksFv_halo,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F, &
           nz-overlap_fb, nz-1)
      call calc_Gv_LES2_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G, &
           overlap_fb, 2*overlap_fb-1)
      call calc_Gv_LES2_koff<<<blocksGv_halo,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G, &
           nz-overlap_fb, nz-1)
    end if
  end subroutine calc_EFG_halo_LES
end module calc_flux_base
