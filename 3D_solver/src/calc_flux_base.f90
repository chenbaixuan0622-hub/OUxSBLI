module calc_flux_base
  use mod_globals, only : id_accuracy, id_igr, &
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
  use calc_igr
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
  subroutine calc_conv_keep(id_scheme, nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    use mod_globals, only : id_accuracy
    integer(2), intent(in), value :: id_scheme
    integer, intent(in), value    :: nx, ny, nz
    real(8), intent(in), device   :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device   :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device   :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device   :: ruvwp(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device  :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device  :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device  :: G(5,nx-2,ny-2,nz-1)
    call calc_keep_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, T, E)
    call calc_keep_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, T, F)
    call calc_keep_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, T, G)
  end subroutine calc_conv_keep


  subroutine calc_conv_slau(id_scheme, nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(2), intent(in), value   :: id_scheme
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: ruvwp(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(8), device :: sensor(nx,ny,nz)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    call calc_slau_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, E)
    call calc_slau_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, F)
    call calc_slau_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, G)
  end subroutine calc_conv_slau


  subroutine calc_conv_roe(id_scheme, nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(4), intent(in), value   :: id_scheme
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: ruvwp(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(8), device :: sensor(nx,ny,nz)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    call calc_roe_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, E)
    call calc_roe_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, F)
    call calc_roe_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, G)
  end subroutine calc_conv_roe


  subroutine calc_conv_hybrid(id_scheme, nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    use mod_globals, only : id_accuracy
    real(8), intent(in), value   :: id_scheme
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: ruvwp(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(8), device :: sensor(nx,ny,nz)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    call calc_hybrid_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, E)
    call calc_hybrid_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, F)
    call calc_hybrid_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, G)
  end subroutine calc_conv_hybrid

  
  subroutine calc_conv_igr(nx, ny, nz, dx, dy, dz, ruvwp, T, sigma, E, F, G)
    use mod_globals, only : id_accuracy
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: ruvwp(5,nx,ny,nz), T(nx,ny,nz)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    real(8), intent(out), device   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device   :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device   :: G(5,nx-2,ny-2,nz-1)
    real(8), device :: sensor(nx,ny,nz)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    call calc_sigma(nx, ny, nz, dx, dy, dz, sensor, ruvwp, sigma)
    call calc_keep_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, T, E, sigma)
    call calc_keep_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, T, F, sigma)
    call calc_keep_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, T, G, sigma)
  end subroutine calc_conv_igr


  subroutine calc_EFG_Euler(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
    use mod_globals, only : id_scheme
    integer(kind=2), intent(in), value   :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(1,1,1), mut(1,1,1), qc2(1,1,1)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    integer stat, i, j, k
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T)
    if (kind(id_igr) == 4) then
      call calc_conv_igr(nx, ny, nz, dx, dy, dz, ruvwp, T, sigma, E, F, G)
    else
      call calc_conv(id_scheme, nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_Euler

 
  subroutine calc_EFG_visc(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
    use mod_globals, only : id_scheme
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(nx,ny,nz), mut(1,1,1), qc2(1,1,1)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    integer stat
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T, mu)
    if (kind(id_igr) == 4) then
      call calc_conv_igr(nx, ny, nz, dx, dy, dz, ruvwp, T, sigma, E, F, G)
    else
      call calc_conv(id_scheme, nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    endif
    stat = cudaDeviceSynchronize()
    if (present(seed)) then
      if (id_visc == 2) then
        call calc_Ev4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E, seed)
        call calc_Fv4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F, seed)
        call calc_Gv4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G, seed)
      else
        call calc_Ev2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E, seed)
        call calc_Fv2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F, seed)
        call calc_Gv2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G, seed)
      endif
      call update_seed(nx, ny, nz, seed)
    else
      if (id_visc == 2) then
        call calc_Ev4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E)
        call calc_Fv4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F)
        call calc_Gv4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G)
      else
        call calc_Ev2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E)
        call calc_Fv2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F)
        call calc_Gv2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G)
      endif
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc

  
  subroutine calc_EFG_LES(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
    use mod_globals, only : id_scheme
    integer(kind=8), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    integer stat
    mut = 0.d0
    qc2 = 0.d0
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T, mu)
    if (kind(id_igr) == 4) then
      call calc_conv_igr(nx, ny, nz, dx, dy, dz, ruvwp, T, sigma, E, F, G)
    else
      call calc_conv(id_scheme, nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    endif
    call calc_mut<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, mut, qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx, ny, nz, mut, qc2)
    if (id_visc == 2) then
      call calc_Ev_LES4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, E)
      call calc_Fv_LES4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, mut, qc2, F)
      call calc_Gv_LES4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, G)
    else
      call calc_Ev_LES2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, E)
      call calc_Fv_LES2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, mut, qc2, F)
      call calc_Gv_LES2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_LES
end module calc_flux_base

