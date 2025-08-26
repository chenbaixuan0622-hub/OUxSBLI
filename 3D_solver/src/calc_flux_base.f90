module calc_flux_base
  use curand
  use mod_globals, only : accuracy, id_accuracy, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_physical_quantities
  use calc_hybrid
  use calc_flux
  use calc_visc
  use calc_les
  use set
  implicit none
  interface calc_EFG
    module procedure calc_EFG_Euler,         calc_EFG_visc,         calc_EFG_LES, &
                     calc_EFG_Euler_forcing, calc_EFG_visc_forcing, calc_EFG_LES_forcing, calc_EFG_fh
  end interface calc_EFG
contains
  subroutine calc_EFG_Euler(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
    integer(kind=2), intent(in), value   :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    real(8), intent(out), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    real(8), intent(out), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    real(8), dimension(5,nx,ny,nz), device :: Q(5,nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_Euler

  
  subroutine calc_EFG_Euler_forcing(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz)
    integer(kind=2), intent(in), value   :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    real(8), intent(out), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    real(8), intent(out), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    real(8), intent(out), device :: fx(nx-2,ny-2,nz-2), fy(nx-2,ny-2,nz-2), fz(nx-2,ny-2,nz-2)
    real(8), dimension(5,nx,ny,nz), device :: Q(5,nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    stat = cudaDeviceSynchronize()
    call calc_forcing(nx, ny, nz, dx, dy, dz, Q, fx, fy, fz)
  end subroutine calc_EFG_Euler_forcing


  subroutine calc_EFG_visc(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    real(8), intent(out), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    real(8), intent(out), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    real(8), dimension(5,nx,ny,nz), device :: Q(5,nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, Q, E)
    call calc_Fv<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, Q, F)
    call calc_Gv<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, Q, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc

  
  subroutine calc_EFG_visc_forcing(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    real(8), intent(out), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    real(8), intent(out), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    real(8), intent(out), device :: fx(nx-2,ny-2,nz-2), fy(nx-2,ny-2,nz-2), fz(nx-2,ny-2,nz-2)
    real(8), dimension(5,nx,ny,nz), device :: Q(5,nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, Q, E)
    call calc_Fv<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, Q, F)
    call calc_Gv<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, Q, G)
    stat = cudaDeviceSynchronize()
    call calc_forcing(nx, ny, nz, dx, dy, dz, Q, fx, fy, fz)
  end subroutine calc_EFG_visc_forcing


  subroutine calc_EFG_LES(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
    integer(kind=8), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    real(8), intent(out), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    real(8), intent(out), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    real(8), dimension(5,nx,ny,nz), device :: Q(5,nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor, mut, qc2
    integer stat
    mut = 0.d0
    qc2 = 0.d0
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    call calc_mut<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, mut, qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx, ny, nz, mut, qc2)
    call calc_Ev_LES<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, Q, mut, qc2, E)
    call calc_Fv_LES<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, Q, mut, qc2, F)
    call calc_Gv_LES<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, Q, mut, qc2, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_LES
  

  subroutine calc_EFG_LES_forcing(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz)
    integer(kind=8), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    real(8), intent(out), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    real(8), intent(out), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    real(8), intent(out), device :: fx(nx-2,ny-2,nz-2), fy(nx-2,ny-2,nz-2), fz(nx-2,ny-2,nz-2)
    real(8), dimension(5,nx,ny,nz), device :: Q(5,nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor, mut, qc2
    integer stat
    mut = 0.d0
    qc2 = 0.d0
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    call calc_mut<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, mut, qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx, ny, nz, mut, qc2)
    call calc_Ev_LES<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, Q, mut, qc2, E)
    call calc_Fv_LES<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, Q, mut, qc2, F)
    call calc_Gv_LES<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, Q, mut, qc2, G)
    stat = cudaDeviceSynchronize()
    call calc_forcing(nx, ny, nz, dx, dy, dz, Q, fx, fy, fz)
  end subroutine calc_EFG_LES_forcing
  
 
  subroutine calc_EFG_fh(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
    real(kind=2), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    real(8), intent(out), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    real(8), intent(out), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    real(8), dimension(5,nx,ny,nz), device :: Q(5,nx,ny,nz)
    real(8), dimension(4,nx,ny,nz), device :: rand(4,nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor
    type(curandGenerator) :: gen
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, Q, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, Q, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
    stat = cudaDeviceSynchronize()
    stat = curandCreateGenerator(gen, CURAND_RNG_PSEUDO_MT19937)
    stat = curandSetPseudoRandomGeneratorSeed(gen, 12345_8)
    stat = curandGenerateNormalDouble(gen, rand, 4*nx*ny*nz, 0.d0, 1.d0)
    stat = curandDestroyGenerator(gen)
    !call calc_Ev_FH<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, Q, rand, E)
    !call calc_Fv_FH<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, Q, rand, F)
    !call calc_Gv_FH<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, Q, rand, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_fh
end module calc_flux_base

