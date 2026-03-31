module calc_steps
  use cudafor
  use cooperative_groups ! for WarpShuffle
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  implicit none
contains
  !> CUDA Fortran kernel for 1st step of 3-3 TVD Runge-Kutta
  attributes(global) subroutine calc_step1(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2)
    integer, intent(in), value   :: nx                  !< number of grid points in x direction
    integer, intent(in), value   :: ny                  !< number of grid points in y direction
    integer, intent(in), value   :: nz                  !< number of grid points in z direction
    real(8), intent(in), value   :: coef                !< coefficient for Runge-Kutta
    real(8), intent(in), device  :: dx(nx-1)            !< grid size in x direction
    real(8), intent(in), device  :: dy(ny-1)            !< grid size in y direction
    real(8), intent(in), device  :: dz(nz-1)            !< grid size in z direction
    real(8), intent(in), device  :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    real(8), intent(in), device  :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    real(8), intent(in), device  :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    real(8), intent(in), device  :: Q(5,nx,ny,nz)       !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device :: Q2(5,nx,ny,nz)      !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Q2(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - coef * R
    enddo
  end subroutine calc_step1


  !> CUDA Fortran kernel for 1st~3rd step of 4-4 Runge-Kutta
  attributes(global) subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in), value     :: nx                   !< number of grid points in x direction
    integer, intent(in), value     :: ny                   !< number of grid points in y direction
    integer, intent(in), value     :: nz                   !< number of grid points in z direction
    real(8), intent(in), value     :: coef1                !< coefficient for Runge-Kutta
    real(8), intent(in), value     :: coef2                !< coefficient for Runge-Kutta
    real(8), intent(in), device    :: dx(nx-1)             !< grid size in x direction
    real(8), intent(in), device    :: dy(ny-1)             !< grid size in y direction
    real(8), intent(in), device    :: dz(nz-1)             !< grid size in z direction
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)  !< Flux in x direction
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)  !< Flux in y direction
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)  !< Flux in z direction
    real(8), intent(in), device    :: Q(5,nx,ny,nz)        !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device   :: Q2(5,nx,ny,nz)       !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(inout), device :: Rs(5,nx-2,ny-2,nz-2) !< accumulation for 4-4 Runge-Kutta
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Q2(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - coef1 * R
      Rs(l,i,j,k) = Rs(l,i,j,k) + coef2 * R
    enddo
  end subroutine calc_step
 

  !> CUDA Fortran kernel for 2nd & 3rd step of 3-3 TVD Runge-Kutta
  attributes(global) subroutine calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Qin, Qout)
    integer, intent(in), value     :: nx                  !< number of grid points in x direction
    integer, intent(in), value     :: ny                  !< number of grid points in y direction
    integer, intent(in), value     :: nz                  !< number of grid points in z direction
    real(8), intent(in), value     :: coef1               !< coefficient for Runge-Kutta
    real(8), intent(in), value     :: coef2               !< coefficient for Runge-Kutta
    real(8), intent(in), value     :: coef3               !< coefficient for Runge-Kutta
    real(8), intent(in), value     :: coef4               !< coefficient for Runge-Kutta
    real(8), intent(in), device    :: dx(nx-1)            !< grid size in x direction
    real(8), intent(in), device    :: dy(ny-1)            !< grid size in y direction
    real(8), intent(in), device    :: dz(nz-1)            !< grid size in z direction
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    real(8), intent(in), device    :: Qin(5,nx,ny,nz)     !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(inout), device :: Qout(5,nx,ny,nz)    !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Qout(l,i+1,j+1,k+1) = (coef1 * Qin(l,i+1,j+1,k+1) + coef2 * Qout(l,i+1,j+1,k+1) - coef3 * R) / coef4
    enddo
  end subroutine calc_step2_3
  
 
  !> CUDA Fortran kernel for 4th step of 4-4 Runge-Kutta
  attributes(global) subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
    integer, intent(in), value     :: nx                   !< number of grid points in x direction
    integer, intent(in), value     :: ny                   !< number of grid points in y direction
    integer, intent(in), value     :: nz                   !< number of grid points in z direction
    real(8), intent(in), device    :: dx(nx-1)             !< grid size in x direction
    real(8), intent(in), device    :: dy(ny-1)             !< grid size in y direction
    real(8), intent(in), device    :: dz(nz-1)             !< grid size in z direction
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)  !< Flux in x direction
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)  !< Flux in y direction
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)  !< Flux in z direction
    real(8), intent(inout), device :: Rs(5,nx-2,ny-2,nz-2) !< accumulation for 4-4 Runge-Kutta
    real(8), intent(inout), device :: Q(5,nx,ny,nz)        !< next Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Rs(l,i,j,k) = Rs(l,i,j,k) + R
      Q(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - Rs(l,i,j,k) * one_sixth
      Rs(l,i,j,k) = 0.d0
    enddo
  end subroutine calc_step4
end module calc_steps

