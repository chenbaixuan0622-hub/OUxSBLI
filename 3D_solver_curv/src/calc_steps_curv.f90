!> Step kernels for curvilinear grid.
!> Differences from calc_steps.f90:
!>   - E and F fluxes are area-scaled by |S_xi| / |S_eta|; dt_xi and dt_eta are scalars (dt*dz).
!>   - G flux is NOT area-scaled; dt_Szeta(nx-2,ny-2) = dt * J_2D(i,j) varies per cell.
module calc_steps_curv
  use cudafor
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_R_curv(nx, ny, nz, i, j, k, dt_xi, dt_eta, dt_Szeta, E, F, G, R)
    integer, intent(in), value              :: nx, ny, nz, i, j, k
    real(8), intent(in), value              :: dt_xi    !< dt * dz  (uniform; E is area-scaled)
    real(8), intent(in), value              :: dt_eta   !< dt * dz  (uniform; F is area-scaled)
    real(8), intent(in), value              :: dt_Szeta !< dt * J_2D(i,j)
    real(8), intent(in), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(out), contiguous        :: R(5)
    R(1) = dt_xi * (-E(1,i,j,k) + E(1,i+1,j,k)) + dt_eta * (-F(1,i,j,k) + F(1,i,j+1,k)) + dt_Szeta * (-G(1,i,j,k) + G(1,i,j,k+1))
    R(2) = dt_xi * (-E(2,i,j,k) + E(2,i+1,j,k)) + dt_eta * (-F(2,i,j,k) + F(2,i,j+1,k)) + dt_Szeta * (-G(2,i,j,k) + G(2,i,j,k+1))
    R(3) = dt_xi * (-E(3,i,j,k) + E(3,i+1,j,k)) + dt_eta * (-F(3,i,j,k) + F(3,i,j+1,k)) + dt_Szeta * (-G(3,i,j,k) + G(3,i,j,k+1))
    R(4) = dt_xi * (-E(4,i,j,k) + E(4,i+1,j,k)) + dt_eta * (-F(4,i,j,k) + F(4,i,j+1,k)) + dt_Szeta * (-G(4,i,j,k) + G(4,i,j,k+1))
    R(5) = dt_xi * (-E(5,i,j,k) + E(5,i+1,j,k)) + dt_eta * (-F(5,i,j,k) + F(5,i,j+1,k)) + dt_Szeta * (-G(5,i,j,k) + G(5,i,j,k+1))
  end subroutine calc_R_curv


  attributes(global) subroutine calc_step1_curv(nx, ny, nz, coef, dt_xi, dt_eta, dt_Szeta, E, F, G, Q, Q2)
    integer, intent(in), value               :: nx, ny, nz
    real(8), intent(in), value               :: coef, dt_xi, dt_eta
    real(8), intent(in), device, contiguous  :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous  :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device, contiguous  :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device, contiguous  :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)
    real(8), intent(out), device, contiguous :: Q2(nx,5,ny,nz)
    real(8) R(5)
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, coef*dt_xi, coef*dt_eta, coef*dt_Szeta(i,j), E, F, G, R)
    do l = 1, 5
      Q2(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - R(l)
    enddo
  end subroutine calc_step1_curv


  attributes(global) subroutine calc_step2_3_curv(nx, ny, nz, coef1, coef2, coef3, coef4_inv, dt_xi, dt_eta, dt_Szeta, E, F, G, Qin, Qout)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), value                 :: coef1, coef2, coef3, coef4_inv, dt_xi, dt_eta
    real(8), intent(in), device, contiguous    :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device, contiguous    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device, contiguous    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device, contiguous    :: Qin(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Qout(nx,5,ny,nz)
    real(8) R(5)
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, coef3*dt_xi, coef3*dt_eta, coef3*dt_Szeta(i,j), E, F, G, R)
    do l = 1, 5
      Qout(i+1,l,j+1,k+1) = (coef1*Qin(i+1,l,j+1,k+1) + coef2*Qout(i+1,l,j+1,k+1) - R(l)) * coef4_inv
    enddo
  end subroutine calc_step2_3_curv


  attributes(global) subroutine calc_step_curv(nx, ny, nz, coef1, coef2, dt_xi, dt_eta, dt_Szeta, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), value                 :: coef1, coef2, dt_xi, dt_eta
    real(8), intent(in), device, contiguous    :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device, contiguous    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device, contiguous    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(out), device, contiguous   :: Q2(nx,5,ny,nz)
    real(8), intent(inout), device, contiguous :: Rs(nx-2,5,ny-2,nz-2)
    real(8) R(5)
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, dt_xi, dt_eta, dt_Szeta(i,j), E, F, G, R)
    do l = 1, 5
      Q2(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - coef1 * R(l)
      Rs(i,l,j,k) = Rs(i,l,j,k) + coef2 * R(l)
    enddo
  end subroutine calc_step_curv


  attributes(global) subroutine calc_step4_curv(nx, ny, nz, dt_xi, dt_eta, dt_Szeta, E, F, G, Rs, Q)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), value                 :: dt_xi, dt_eta
    real(8), intent(in), device, contiguous    :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device, contiguous    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device, contiguous    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(inout), device, contiguous :: Rs(nx-2,5,ny-2,nz-2)
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz)
    real(8) R(5)
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, dt_xi, dt_eta, dt_Szeta(i,j), E, F, G, R)
    do l = 1, 5
      R(l) = Rs(i,l,j,k) + R(l)
      Q(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - R(l) * one_sixth
      Rs(i,l,j,k) = 0.d0
    enddo
  end subroutine calc_step4_curv
end module calc_steps_curv
