module calc_flux_base
  use mod_globals, only : id_accuracy, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_physical_quantities
  use calc_hybrid
  use calc_keep_kernel
  use calc_slau_kernel
  use calc_roe_kernel
  use calc_hybrid_kernel
  use calc_visc2
  use calc_visc4
  use calc_les
  use set
  implicit none
  private
  public calc_EFG
  interface calc_conv
    module procedure calc_conv_keep, calc_conv_slau, calc_conv_roe, calc_conv_hybrid
  end interface calc_conv

  interface calc_EFG
    module procedure calc_EFG_Euler, calc_EFG_visc, calc_EFG_LES
  end interface calc_EFG
contains
  !< calc conv term using KEEP scheme
  subroutine calc_conv_keep(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    integer(2), intent(in), value :: id_scheme           !> ID for scheme, int 2 mesns KEEP
    integer, intent(in), value    :: nx                  !> number of grid points in x direction
    integer, intent(in), value    :: ny                  !> number of grid points in y direction
    integer, intent(in), value    :: nz                  !> number of grid points in x direction
    real(8), intent(in), device   :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device   :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device   :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device   :: Q(5,nx,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(in), device   :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device  :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device  :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device  :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    call calc_keep_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, T, E)
    call calc_keep_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, T, F)
    call calc_keep_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, T, G)
  end subroutine calc_conv_keep


  !< calc conv term using SLAU scheme
  subroutine calc_conv_slau(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(2), intent(in), value   :: id_scheme           !> ID for scheme, real 2 means SLAU
    integer, intent(in), value   :: nx                  !> number of grid points in x direction
    integer, intent(in), value   :: ny                  !> number of grid points in y direction
    integer, intent(in), value   :: nz                  !> number of grid points in x direction
    real(8), intent(in), device  :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device  :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device  :: Q(5,nx,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(in), device  :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    real(8), device :: sensor(nx,ny,nz)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    call calc_slau_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_slau_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_slau_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
  end subroutine calc_conv_slau


  !< calc conv term using Roe scheme
  subroutine calc_conv_roe(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(4), intent(in), value   :: id_scheme           !> ID for scheme, real 4 means Roe
    integer, intent(in), value   :: nx                  !> number of grid points in x direction
    integer, intent(in), value   :: ny                  !> number of grid points in y direction
    integer, intent(in), value   :: nz                  !> number of grid points in x direction
    real(8), intent(in), device  :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device  :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device  :: Q(5,nx,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(in), device  :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    real(8), device :: sensor(nx,ny,nz)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    call calc_roe_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_roe_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_roe_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
  end subroutine calc_conv_roe


  !< calc conv term using KEEP/SLAU Hybrid scheme
  subroutine calc_conv_hybrid(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(8), intent(in), value   :: id_scheme           !> ID for scheme, real 8 means Hybrid
    integer, intent(in), value   :: nx                  !> number of grid points in x direction
    integer, intent(in), value   :: ny                  !> number of grid points in y direction
    integer, intent(in), value   :: nz                  !> number of grid points in x direction
    real(8), intent(in), device  :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device  :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device  :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device  :: Q(5,nx,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(in), device  :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    real(8), device :: sensor(nx,ny,nz)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, sensor)
    call calc_hybrid_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    call calc_hybrid_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    call calc_hybrid_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, T, sensor, G)
  end subroutine calc_conv_hybrid


  !< calc Flux of Euler equation
  subroutine calc_EFG_Euler(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)
    use mod_globals, only : id_scheme
    integer(2), intent(in), value :: id_visc             !> ID for equation, int 2 means Euler
    integer, intent(in), value    :: nx                  !> number of grid points in x direction
    integer, intent(in), value    :: ny                  !> number of grid points in y direction
    integer, intent(in), value    :: nz                  !> number of grid points in z direction
    real(8), intent(in), device   :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device   :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device   :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device   :: Jacobian(nx,ny)     !> Jacobian
    real(8), intent(in), device   :: QJ(5,nx,ny,nz)      !> Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device  :: Q(5,nx,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(out), device  :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device  :: mu(1,1,1)           !> viscosity, size is (1,1,1) in case of Euler
    real(8), intent(out), device  :: mut(1,1,1)          !> SGS viscosity, size is (1,1,1) in case of Euler
    real(8), intent(out), device  :: qc2(1,1,1)          !> SGS kinetic energy, size is (1,1,1) in case of Euler
    real(8), intent(out), device  :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device  :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device  :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    integer stat, i, j, k
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T)
    call calc_conv(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
  end subroutine calc_EFG_Euler


  !< calc FLux of Navier-Stokes equation
  subroutine calc_EFG_visc(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)
    use mod_globals, only : id_scheme
    integer(4), intent(in), value :: id_visc             !> ID for equation, int 4 means NS
    integer, intent(in), value    :: nx                  !> number of grid points in x direction
    integer, intent(in), value    :: ny                  !> number of grid points in y direction
    integer, intent(in), value    :: nz                  !> number of grid points in z direction
    real(8), intent(in), device   :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device   :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device   :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device   :: Jacobian(nx,ny)     !> Jacobian
    real(8), intent(in), device   :: QJ(5,nx,ny,nz)      !> Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device  :: Q(5,nx,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(out), device  :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device  :: mu(nx,ny,nz)        !> viscosity
    real(8), intent(out), device  :: mut(1,1,1)          !> SGS viscosity, size is (1,1,1) in case of NS
    real(8), intent(out), device  :: qc2(1,1,1)          !> SGS kinetic energy, size is (1,1,1) in case of NS
    real(8), intent(out), device  :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device  :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device  :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    integer stat
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu)
    call calc_conv(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    if (id_visc == 2) then
      call calc_Ev4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
      call calc_Fv4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
      call calc_Gv4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
    else
      call calc_Ev2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
      call calc_Fv2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
      call calc_Gv2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc

 
  !< calc FLux of LES
  subroutine calc_EFG_LES(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)
    use mod_globals, only : id_scheme
    integer(8), intent(in), value :: id_visc             !> ID for equation, int 8 means LES
    integer, intent(in), value    :: nx                  !> number of grid points in x direction
    integer, intent(in), value    :: ny                  !> number of grid points in y direction
    integer, intent(in), value    :: nz                  !> number of grid points in z direction
    real(8), intent(in), device   :: inv_dx(nx-1)        !> 1 / dx
    real(8), intent(in), device   :: inv_dy(ny-1)        !> 1 / dy
    real(8), intent(in), device   :: inv_dz(nz-1)        !> 1 / dz
    real(8), intent(in), device   :: Jacobian(nx,ny)     !> Jacobian
    real(8), intent(in), device   :: QJ(5,nx,ny,nz)      !> Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device  :: Q(5,nx,ny,nz)       !> Q(rho, u, v, w, p)
    real(8), intent(out), device  :: T(nx,ny,nz)         !> temperature
    real(8), intent(out), device  :: mu(nx,ny,nz)        !> viscosity
    real(8), intent(out), device  :: mut(nx,ny,nz)       !> SGS viscosity
    real(8), intent(out), device  :: qc2(nx,ny,nz)       !> SGS kinetic energy
    real(8), intent(out), device  :: E(5,nx-1,ny-2,nz-2) !> Flux in x direction
    real(8), intent(out), device  :: F(5,nx-2,ny-1,nz-2) !> Flux in y direction
    real(8), intent(out), device  :: G(5,nx-2,ny-2,nz-1) !> Flux in z direction
    integer stat
    mut = 0.d0
    qc2 = 0.d0
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu)
    call calc_conv(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)
    call calc_mut<<<blocks,threads>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, mut, qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx, ny, nz, mut, qc2)
    if (id_visc == 2) then
      call calc_Ev_LES4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E)
      call calc_Fv_LES4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F)
      call calc_Gv_LES4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
    else
      call calc_Ev_LES2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, E)
      call calc_Fv_LES2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, mut, qc2, F)
      call calc_Gv_LES2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, mut, qc2, G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_LES
end module calc_flux_base

