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
  !use calc_roe_kernel
  !use calc_roe_kernel_internal
  use calc_hybrid_kernel
  use calc_hybrid_kernel_internal
  ! use calc_visc2
  ! use calc_visc4
  ! use calc_visc4_internal
  ! use calc_visc4_les_internal
  !use calc_les
  use set
  implicit none
  private
  public calc_EFG, init_sensor
  !> Persistent Ducros shock sensor on device — allocated once in init_sensor,
  !> reused across all calls to calc_conv_slau / calc_conv_roe / calc_conv_hybrid.
  !> Eliminates repeated device heap alloc/free on every RK stage.
  real(sp), allocatable, device, save :: sensor(:,:)
  interface calc_conv
    module procedure calc_conv_keep, calc_conv_slau, calc_conv_hybrid ! calc_conv_roe
  end interface calc_conv

  interface calc_EFG
    module procedure calc_EFG_Euler, calc_EFG_visc, calc_EFG_LES
  end interface calc_EFG
contains
  !> Allocate the persistent Ducros sensor array on the device.
  !> Call once from preprocess.f90 before the time loop.
  subroutine init_sensor(nx, ny) !2d
    integer, intent(in) :: nx, ny
    allocate(sensor(nx,ny))
  end subroutine init_sensor

  !> Compute convective fluxes using KEEP (energy-preserving) scheme
  !> High-order minimal dissipation scheme for smooth flow regions
  !> Algorithm: F = (H·u) where H = enthalpy, using flux reconstruction via divergence forms
  !> Provides 4th-5th order accuracy by minimizing dispersive errors in smooth regions
  subroutine calc_conv_keep(id_scheme, nx, ny, inv_dx, inv_dy,  Q, T, E, F)
    use mod_globals, only : id_accuracy
    integer(2), intent(in), value            :: id_scheme           !< ID for scheme: int 2 means KEEP
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    !integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !< inverse grid spacing y (1/dy)
    !real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !< inverse grid spacing z (1/dz)
    real(8), intent(in), device, contiguous  :: Q(nx,4,ny)       !< conservative variables Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous  :: T(nx,ny)         !< temperature field
    real(8), intent(out), device, contiguous :: E(4,nx-1,ny-2) !< convective flux in x direction
    real(8), intent(out), device, contiguous :: F(4,nx-2,ny-1) !< convective flux in y direction
    !real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2) !< convective flux in z direction
    if (id_bc_x) then
      call calc_keep_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny,  Q, T, E)
    else
      call calc_keep_x_in<<<blocksE,threadsE>>>(nx, ny,  Q, T, E)
    endif
    if (id_bc_y) then
      call calc_keep_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, Q, T, F)
    else
      call calc_keep_y_in<<<blocksF,threadsF>>>(nx, ny,  Q, T, F)
    endif
    !if (id_bc_z) then
    !  call calc_keep_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, T, G)
    !else
    !  call calc_keep_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, T, G)
    !endif
  end subroutine calc_conv_keep


  !> Compute convective fluxes using SLAU (Simple Low-dissipation Roe-based Upwind) scheme
  !> Low-dissipation scheme with shock-capturing capability via Ducros sensor
  !> Algorithm: F = (F_L + F_R)/2 + |A|(Q_L - Q_R)/2 where A is weighted Jacobian
  !> Dissipation modulated by Ducros sensor: f_d controls blend ratio
  subroutine calc_conv_slau(id_scheme, nx, ny,  inv_dx, inv_dy, Q, T, E, F)
    use mod_globals, only : id_accuracy
    real(2), intent(in), value               :: id_scheme           !< ID for scheme: real 2 means SLAU
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    !integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !< inverse grid spacing y (1/dy)
    !real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !< inverse grid spacing z (1/dz)
    real(8), intent(in), device, contiguous  :: Q(nx,4,ny)       !< conservative variables
    real(8), intent(in), device, contiguous  :: T(nx,ny)         !< temperature field
    real(8), intent(out), device, contiguous :: E(4,nx-1,ny-2) !< convective flux in x direction
    real(8), intent(out), device, contiguous :: F(4,nx-2,ny-1) !< convective flux in y direction
    !real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2) !< convective flux in z direction
    call calc_Ducros<<<blocks,threads>>>(nx, ny, inv_dx, inv_dy, Q, sensor)
    if (id_bc_x) then
      call calc_slau_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny,  Q, sensor, E)
    else
      call calc_slau_x_in<<<blocksE,threadsE>>>(nx, ny,  Q, sensor, E)
    endif
    if (id_bc_y) then
      call calc_slau_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny,Q, sensor, F)
    else
      call calc_slau_y_in<<<blocksF,threadsF>>>(nx, ny, Q, sensor, F)
    endif
    ! if (id_bc_z) then
    !   call calc_slau_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    ! else
    !   call calc_slau_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, sensor, G)
    ! endif
  end subroutine calc_conv_slau


  !> Compute convective fluxes using Roe approximate Riemann solver
  !> Classic approximate Riemann solver with wave decomposition for flux splitting
  !> Algorithm: F = (F_L + F_R)/2 - (1/2)Σ|λ_i|*(p_i·r_i) wave reconstruction
  !> Where λ_i are Roe eigenvalues, p_i are wave strengths, r_i are eigenvectors
  !> Entropy fix via sensor prevents expansion shocks at sonic points
  ! subroutine calc_conv_roe(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
  !   use mod_globals, only : id_accuracy
  !   real(4), intent(in), value               :: id_scheme           !< ID for scheme: real 4 means Roe
  !   integer, intent(in), value               :: nx                  !< number of grid points in x direction
  !   integer, intent(in), value               :: ny                  !< number of grid points in y direction
  !   integer, intent(in), value               :: nz                  !< number of grid points in z direction
  !   real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
  !   real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !< inverse grid spacing y (1/dy)
  !   real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !< inverse grid spacing z (1/dz)
  !   real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)       !< conservative variables
  !   real(8), intent(in), device, contiguous  :: T(nx,ny,nz)         !< temperature field
  !   real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< convective flux in x direction
  !   real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< convective flux in y direction
  !   real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !< convective flux in z direction
  !   call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
  !   if (id_bc_x) then
  !     call calc_roe_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
  !   else
  !     call calc_roe_x_in<<<blocksE,threadsE>>>(nx, ny, nz, Q, sensor, E)
  !   endif
  !   if (id_bc_y) then
  !     call calc_roe_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
  !   else
  !     call calc_roe_y_in<<<blocksF,threadsF>>>(nx, ny, nz, Q, sensor, F)
  !   endif
  !   if (id_bc_z) then
  !     call calc_roe_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
  !   else
  !     call calc_roe_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, sensor, G)
  !   endif
  ! end subroutine calc_conv_roe


  !> Compute convective fluxes using hybrid KEEP/SLAU scheme  
  !> Automatically blends between KEEP (smooth regions) and SLAU (shock regions) seamlessly
  !> Blending formula: F_hybrid = (1-f_d)·F_keep + f_d·F_slau where f_d ∈ [0,1]
  !> Preserves vortex structures (f_d≈0) and captures shocks accurately (f_d≈1)
  !> Ducros shock sensor: f_d = (∇·u)²/[(∇·u)² + (∇×u)² + ε]
  subroutine calc_conv_hybrid(id_scheme, nx, ny,  inv_dx, inv_dy,  Q, T, E, F)
    use mod_globals, only : id_accuracy
    real(8), intent(in), value               :: id_scheme           !< ID for scheme: real 8 means Hybrid
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    !integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !< inverse grid spacing x (1/dx)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !> 1 / dy
    !real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny)       !> Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous  :: T(nx,ny)         !> temperature
    real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2) !> Flux in x direction
    real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1) !> Flux in y direction
    !real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    call calc_Ducros<<<blocks,threads>>>(nx, ny, inv_dx, inv_dy,Q, sensor)
    if (id_bc_x) then
      call calc_hybrid_x<<<blocksE,threadsE>>>(id_accuracy, nx, ny, Q, T, sensor, E)
    else
      call calc_hybrid_x_in<<<blocksE,threadsE>>>(nx, ny, Q, T, sensor, E)
    endif
    if (id_bc_y) then
      call calc_hybrid_y<<<blocksF,threadsF>>>(id_accuracy, nx, ny, Q, T, sensor, F)
    else
      call calc_hybrid_y_in<<<blocksF,threadsF>>>(nx, ny, Q, T, sensor, F)
    endif
    !if (id_bc_z) then
   !   call calc_hybrid_z<<<blocksG,threadsG>>>(id_accuracy, nx, ny, Q, T, sensor, G)
    !else
     ! call calc_hybrid_z_in<<<blocksG,threadsG>>>(nx, ny, nz, Q, T, sensor, G)
   ! endif
  end subroutine calc_conv_hybrid


  !< calc Flux of Euler equation
  subroutine calc_EFG_Euler(id_visc, nx, ny, inv_dx, inv_dy, Jacobian, QJ, Q, T, mu, mut, qc2, E, F)
    use mod_globals, only : id_scheme
    integer(2), intent(in), value            :: id_visc             !> ID for equation, int 2 means Euler
    integer, intent(in), value               :: nx                  !> number of grid points in x direction
    integer, intent(in), value               :: ny                  !> number of grid points in y direction
    !integer, intent(in), value               :: nz                  !> number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !> 1 / dy
    !real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device, contiguous  :: Jacobian(nx,ny)     !> Jacobian
    real(8), intent(in), device, contiguous  :: QJ(nx,4,ny)      !> Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous :: Q(nx,4,ny)       !> Q(rho, u, v, w, p)
    real(8), intent(out), device, contiguous :: T(nx,ny)         !> temperature
    real(8), intent(out), device, contiguous :: mu(1,1)           !> viscosity, size is (1,1,1) in case of Euler
    real(8), intent(out), device, contiguous :: mut(1,1)          !> SGS viscosity, size is (1,1,1) in case of Euler
    real(8), intent(out), device, contiguous :: qc2(1,1)          !> SGS kinetic energy, size is (1,1,1) in case of Euler
    real(8), intent(out), device, contiguous :: E(4,nx-1,ny-2) !> Flux in x direction
    real(8), intent(out), device, contiguous :: F(4,nx-2,ny-1) !> Flux in y direction
    !real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    integer stat, i, j, k
    call calc_quantities_2D(nx, ny,  Jacobian, QJ, Q, T)
    !call calc_quantities_3D(nx, ny,  Jacobian, QJ, Q, T)
    call calc_conv(id_scheme, nx, ny, inv_dx, inv_dy, Q, T, E, F)
  end subroutine calc_EFG_Euler


  !> Compute fluxes for viscous (Navier-Stokes) flow - convective + viscous components
  !> Computes stress tensor tau_ij = mu*(du_i/dx_j + du_j/dx_i) - (2/3)*mu*delta_ij*(div u)
  !> and heat flux via Fourier's law: q = -k*dT/dx where k depends on Prandtl number
  subroutine calc_EFG_visc(id_visc, nx, ny, inv_dx, inv_dy,  Jacobian, QJ, Q, T, mu, mut, qc2, E, F)
    use mod_globals, only : id_scheme
    integer(4), intent(in), value            :: id_visc             !> ID for equation, int 4 means NS
    integer, intent(in), value               :: nx                  !> number of grid points in x direction
    integer, intent(in), value               :: ny                  !> number of grid points in y direction
    !integer, intent(in), value               :: nz                  !> number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !> 1 / dx (for finite differences)
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !> 1 / dy
    !real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device, contiguous  :: Jacobian(nx,ny)     !> Jacobian determinant for scaling
    real(8), intent(in), device, contiguous  :: QJ(nx,4,ny)      !> Q/Jacobian (scaled conserved variables)
    real(8), intent(out), device, contiguous :: Q(nx,4,ny)       !> Q(rho, u, v, w, p) primitive variables
    real(8), intent(out), device, contiguous :: T(nx,ny)         !> temperature field (for viscosity & heat flux)
    real(8), intent(out), device, contiguous :: mu(nx,ny)        !> molecular viscosity via Sutherland's law
    real(8), intent(out), device, contiguous :: mut(1,1)          !> SGS turbulent viscosity (unused for NS)
    real(8), intent(out), device, contiguous :: qc2(1,1)          !> SGS kinetic energy (unused for NS)
    real(8), intent(out), device, contiguous :: E(4,nx-1,ny-2) !> x-direction flux (convective + viscous)
    real(8), intent(out), device, contiguous :: F(4,nx-2,ny-1) !> y-direction flux (convective + viscous)
    !real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1) !> z-direction flux (convective + viscous)
    integer stat
    ! Step 1: Decode Q and compute T(rho) and mu(T) via Sutherland's formula
    call calc_quantities_T_2D(nx, ny, Jacobian, QJ, Q, T, mu)   !call calc_quantities_T_3D(nx, ny, Jacobian, QJ, Q, T, mu)
    ! Step 2: Compute convective fluxes (KEEP/SLAU/Roe/Hybrid depending on id_scheme)
    call calc_conv(id_scheme, nx, ny, inv_dx, inv_dy, Q, T, E, F)
    ! Step 3: Add viscous fluxes (choose 2nd or 4th-order stencils)
    if (id_visc == 2) then
      ! 4th-order compact finite differences (higher accuracy, larger stencil)
      if (id_bc_x == .false. .and. kind(id_accuracy) == 8) then
        call calc_Ev4_in<<<blocksEv,threadsEv>>>(nx, ny, inv_dx, inv_dy, Q, T, mu, E)
      else
        call calc_Ev4<<<blocksEv,threadsEv>>>(nx, ny, inv_dx, inv_dy, Q, T, mu, E)
      endif
      if (id_bc_y == .false. .and. kind(id_accuracy) == 8) then
        call calc_Fv4_in<<<blocksFv,threadsFv>>>(nx, ny, inv_dy, inv_dx, Q, T, mu, F)
      else
        call calc_Fv4<<<blocksFv,threadsFv>>>(nx, ny, inv_dy, inv_dx, Q, T, mu, F)
      endif
      !if (id_bc_z == .false. .and. kind(id_accuracy) == 8) then
      !  call calc_Gv4_in<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
      !else
        !call calc_Gv4<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
      !endif
    else
      ! 2nd-order centered differences (standard, 3-point stencil)
      call calc_Ev2<<<blocksEv,threadsEv>>>(nx, ny,inv_dx, inv_dy,  Q, T, mu, E)
      call calc_Fv2<<<blocksFv,threadsFv>>>(nx, ny,inv_dy, inv_dx,  Q, T, mu, F)
      !call calc_Gv2<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
    endif
  end subroutine calc_EFG_visc

 
  !> Compute fluxes for Large-Eddy Simulation (LES) with subgrid-scale modeling
  !> Combines molecular viscosity (Navier-Stokes) with turbulent viscosity (from Smagorinsky model)
  !> Filters out subgrid scales: nu_t = (C_s * Delta)^2 * |S_ij| where Delta is grid filter width
  subroutine calc_EFG_LES(id_visc, nx, ny, inv_dx, inv_dy, Jacobian, QJ, Q, T, mu, mut, qc2, E, F)
    use mod_globals, only : id_scheme
    integer(8), intent(in), value            :: id_visc             !> ID for equation, int 8 means LES
    integer, intent(in), value               :: nx                  !> number of grid points in x direction
    integer, intent(in), value               :: ny                  !> number of grid points in y direction
    !integer, intent(in), value               :: nz                  !> number of grid points in z direction
    real(8), intent(in), device, contiguous  :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device, contiguous  :: inv_dy(ny-1)        !> 1 / dy
   ! real(8), intent(in), device, contiguous  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device, contiguous  :: Jacobian(nx,ny)     !> Jacobian
    real(8), intent(in), device, contiguous  :: QJ(nx,4,ny)      !> Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous :: Q(nx,4,ny)       !> Q(rho, u, v, w, p)
    real(8), intent(out), device, contiguous :: T(nx,ny)         !> temperature
    real(8), intent(out), device, contiguous :: mu(nx,ny)        !> viscosity
    real(8), intent(out), device, contiguous :: mut(nx,ny)       !> SGS viscosity
    real(8), intent(out), device, contiguous :: qc2(nx,ny)       !> SGS kinetic energy
    real(8), intent(out), device, contiguous :: E(4,nx-1,ny-2) !> Flux in x direction
    real(8), intent(out), device, contiguous :: F(4,nx-2,ny-1) !> Flux in y direction
    !real(8), intent(out), device, contiguous :: G(4,nx-2,ny-2,) !> Flux in z direction
    integer stat
    mut = 0.d0
    qc2 = 0.d0
    call calc_quantities_T_2D(nx, ny, Jacobian, QJ, Q, T, mu) !    call calc_quantities_T_3D(nx, ny, Jacobian, QJ, Q, T, mu)
    call calc_conv(id_scheme, nx, ny, inv_dx, inv_dy, Q, T, E, F)
    call calc_mut<<<blocks,threads>>>(nx, ny, inv_dx, inv_dy, Q, mut, qc2)
    call set_bc_mut(nx, ny, mut, qc2)
    if (id_visc == 2) then
      if (id_bc_x == .false. .and. kind(id_accuracy) == 8) then
        call calc_Ev_LES4_in<<<blocksEv,threadsEv>>>(nx, ny, inv_dx, inv_dy, Q, T, mu, mut, qc2, E)
      else
        call calc_Ev_LES4<<<blocksEv,threadsEv>>>(nx, ny, inv_dx, inv_dy,  Q, T, mu, mut, qc2, E)
      endif
      if (id_bc_y == .false. .and. kind(id_accuracy) == 8) then
        call calc_Fv_LES4_in<<<blocksFv,threadsFv>>>(nx, ny, inv_dy, inv_dx,  Q, T, mu, mut, qc2, F)
      else
        call calc_Fv_LES4<<<blocksFv,threadsFv>>>(nx, ny, inv_dy, inv_dx, Q, T, mu, mut, qc2, F)
      endif
      !if (id_bc_z == .false. .and. kind(id_accuracy) == 8) then
      !  call calc_Gv_LES4_in<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
      !else
      !  call calc_Gv_LES4<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
      !endif
    else
      call calc_Ev_LES2<<<blocksEv,threadsEv>>>(nx, ny,  inv_dx, inv_dy, Q, T, mu, mut, qc2, E)
      call calc_Fv_LES2<<<blocksFv,threadsFv>>>(nx, ny,  inv_dy, inv_dx, Q, T, mu, mut, qc2, F)
      !call calc_Gv_LES2<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
    endif
  end subroutine calc_EFG_LES
end module calc_flux_base

