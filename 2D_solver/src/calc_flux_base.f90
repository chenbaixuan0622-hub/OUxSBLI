module calc_flux_base
  use mod_globals, only : id_accuracy, id_igr, &
  & blocks, threads, blocksE, blocksF, threadsE, threadsF, &
  & blocksEv, blocksFv, threadsEv, threadsFv
  use calc_physical_quantities
  use calc_hybrid
  use calc_flux
  use calc_visc2
  use calc_visc4
  use set
  use calc_igr
  implicit none
  interface calc_EF
    module procedure calc_EF_Euler, calc_EF_visc
  end interface calc_EF
contains
  subroutine calc_EF_Euler(id_visc, nx, ny, dx, dy, Jacobian, QJ, ruvp, T, mu, E, F, sigma, seed)
    integer(kind=2), intent(in), value   :: id_visc
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: ruvp(4,nx,ny) ! (rho, u, v, p)
    real(8), intent(out), device :: T(nx,ny), mu(1,1) 
    real(8), intent(out), device :: E(4,nx-1,ny-2)
    real(8), intent(out), device :: F(4,nx-2,ny-1)
    real(4), intent(inout), device :: sigma(nx,ny)
    integer(8), intent(inout), device, optional :: seed(nx,ny)
    real(8), dimension(nx,ny), device :: sensor
    integer stat, i, j
    call calc_quantities_2D(nx, ny, Jacobian, QJ, ruvp, T)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, dx, dy, ruvp, sensor)
    if (kind(id_igr) == 4) then
      call calc_sigma(nx, ny, dx, dy, sensor, ruvp, sigma)
      call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, ruvp, T, sensor, E, sigma)
      call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, ruvp, T, sensor, F, sigma)
    else
      call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, ruvp, T, sensor, E)
      call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, ruvp, T, sensor, F)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EF_Euler

 
  subroutine calc_EF_visc(id_visc, nx, ny, dx, dy, Jacobian, QJ, ruvp, T, mu, E, F, sigma, seed)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: ruvp(4,nx,ny) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny), mu(nx,ny)
    real(8), intent(out), device :: E(4,nx-1,ny-2)
    real(8), intent(out), device :: F(4,nx-2,ny-1)
    real(4), intent(inout), device :: sigma(nx,ny)
    integer(8), intent(inout), device, optional :: seed(nx,ny)
    real(8), dimension(nx,ny), device :: sensor
    integer stat
    call calc_quantities_T_2D(nx, ny, Jacobian, QJ, ruvp, T, mu)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, dx, dy, ruvp, sensor)
    if (kind(id_igr) == 4) then
      call calc_sigma(nx, ny, dx, dy, sensor, ruvp, sigma)
      call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, ruvp, T, sensor, E, sigma)
      call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, ruvp, T, sensor, F, sigma)
    else
      call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, ruvp, T, sensor, E)
      call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, ruvp, T, sensor, F)
    endif
    stat = cudaDeviceSynchronize()
    if (present(seed)) then
      if (id_visc == 2) then
        call calc_Ev4<<<blocksEv,threadsEv,1>>>(nx, ny, dx, dy, ruvp, T, mu, E, seed)
        call calc_Fv4<<<blocksFv,threadsFv,2>>>(nx, ny, dy, dx, ruvp, T, mu, F, seed)
      else
        call calc_Ev2<<<blocksEv,threadsEv,1>>>(nx, ny, dx, dy, ruvp, T, mu, E, seed)
        call calc_Fv2<<<blocksFv,threadsFv,2>>>(nx, ny, dy, dx, ruvp, T, mu, F, seed)
      endif
      call update_seed(nx, ny, seed)
    else
      if (id_visc == 2) then
        call calc_Ev4<<<blocksEv,threadsEv,1>>>(nx, ny, dx, dy, ruvp, T, mu, E)
        call calc_Fv4<<<blocksFv,threadsFv,2>>>(nx, ny, dy, dx, ruvp, T, mu, F)
      else
        call calc_Ev2<<<blocksEv,threadsEv,1>>>(nx, ny, dx, dy, ruvp, T, mu, E)
        call calc_Fv2<<<blocksFv,threadsFv,2>>>(nx, ny, dy, dx, ruvp, T, mu, F)
      endif
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EF_visc
end module calc_flux_base

