module calc_flux_base
  use mod_globals, only : accuracy, id_accuracy, &
  & blocks, threads, blocksE, blocksF, threadsE, threadsF, &
  & blocksEv, blocksFv, threadsEv, threadsFv
  use calc_physical_quantities
  use calc_hybrid
  use calc_flux
  use calc_visc
  use set
  implicit none
  interface calc_EFG
    module procedure calc_EFG_Euler,         calc_EFG_visc, &
                     calc_EFG_Euler_forcing, calc_EFG_visc_forcing
  end interface calc_EFG
contains
  subroutine calc_EFG_Euler(id_visc, nx, ny, dx, dy, Jacobian, QJ, E, F)
    integer(kind=2), intent(in), value   :: id_visc
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: Jacobian(ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: E(4,nx-accuracy+1,ny-accuracy)
    real(8), intent(out), device :: F(4,nx-accuracy,ny-accuracy+1)
    real(8), dimension(nx,ny), device :: rho, u, v, p, sensor
    integer stat
    call calc_quantities_2D(nx, ny, Jacobian, QJ, rho, u, v, p)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, dx, dy, u, v, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, rho, u, v, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, rho, u, v, p, sensor, F)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_Euler

  
  subroutine calc_EFG_Euler_forcing(id_visc, nx, ny, dx, dy, Jacobian, QJ, E, F, fx, fy)
    integer(kind=2), intent(in), value   :: id_visc
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: Jacobian(ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: E(4,nx-accuracy+1,ny-accuracy)
    real(8), intent(out), device :: F(4,nx-accuracy,ny-accuracy+1)
    real(8), intent(out), device :: fx(nx-2,ny-2), fy(nx-2,ny-2)
    real(8), dimension(nx,ny), device :: rho, u, v, p, sensor
    integer stat
    call calc_quantities_2D(nx, ny, Jacobian, QJ, rho, u, v, p)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, dx, dy, u, v, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, rho, u, v, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, rho, u, v, p, sensor, F)
    stat = cudaDeviceSynchronize()
    call calc_forcing(nx, ny, dx, dy, rho, u, v, p, fx, fy)
  end subroutine calc_EFG_Euler_forcing


  subroutine calc_EFG_visc(id_visc, nx, ny, dx, dy, Jacobian, QJ, E, F)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: Jacobian(ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: E(4,nx-accuracy+1,ny-accuracy)
    real(8), intent(out), device :: F(4,nx-accuracy,ny-accuracy+1)
    real(8), dimension(nx,ny), device :: rho, u, v, p, sensor
    integer stat
    call calc_quantities_2D(nx, ny, Jacobian, QJ, rho, u, v, p)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, dx, dy, u, v, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, rho, u, v, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, rho, u, v, p, sensor, F)
    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksEv,threadsEv,1>>>(nx, ny, dx, dy, rho, u, v, p, E)
    call calc_Fv<<<blocksFv,threadsFv,2>>>(nx, ny, dy, dx, rho, u, v, p, F)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc


  subroutine calc_EFG_visc_forcing(id_visc, nx, ny, dx, dy, Jacobian, QJ, E, F, fx, fy)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: Jacobian(ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: E(4,nx-accuracy+1,ny-accuracy)
    real(8), intent(out), device :: F(4,nx-accuracy,ny-accuracy+1)
    real(8), intent(out), device :: fx(nx-2,ny-2), fy(nx-2,ny-2)
    real(8), dimension(nx,ny), device :: rho, u, v, p, sensor
    integer stat
    call calc_quantities_2D(nx, ny, Jacobian, QJ, rho, u, v, p)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, dx, dy, u, v, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, rho, u, v, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, rho, u, v, p, sensor, F)
    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksEv,threadsEv,1>>>(nx, ny, dx, dy, rho, u, v, p, E)
    call calc_Fv<<<blocksFv,threadsFv,2>>>(nx, ny, dy, dx, rho, u, v, p, F)
    stat = cudaDeviceSynchronize()
    call calc_forcing(nx, ny, dx, dy, rho, u, v, p, fx, fy)
  end subroutine calc_EFG_visc_forcing
end module calc_flux_base

