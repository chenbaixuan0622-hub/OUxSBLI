module calc_steps
  use cudafor
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  implicit none
contains
  !> CUDA Fortran kernel for 1st step of 3-3 TVD Runge-Kutta
  attributes(global) subroutine calc_step1(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2)
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), value               :: coef                !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous  :: dx(nx-1)            !< grid size in x direction
    real(8), intent(in), device, contiguous  :: dy(ny-1)            !< grid size in y direction
    real(8), intent(in), device, contiguous  :: dz(nz-1)            !< grid size in z direction
    real(8), intent(in), device, contiguous  :: E(nx-1,5,ny-2,nz-2) !< Flux in x direction
    real(8), intent(in), device, contiguous  :: F(nx-2,5,ny-1,nz-2) !< Flux in y direction
    real(8), intent(in), device, contiguous  :: G(nx-2,5,ny-2,nz-1) !< Flux in z direction
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)       !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous :: Q2(nx,5,ny,nz)      !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== Compute Grid Jacobian Volumes ==========
    ! Cell volume at (i, j, k) to (i+1, j+1, k+1) is product of grid spacings
    ! For flux divergence, use average grid spacing between cell centers
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1)) ! dt*dy*dz x-flux divergence
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1)) ! dt*dz*dx y-flux divergence
    dtdxdy = dt * 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1)) ! dt*dx*dy z-flux divergence
    ! ========== Conservative Update via TVD RK3: Stage 1 ==========
    ! Q^(1) = Q^n - (coef) * dt/vol * (Flux_divergence)
    do l = 1, 5  ! Loop over all conserved variables (rho, rhou, rhov, rhow, E)
      R = dtdydz * (-E(i,l,j,k) + E(i+1,l,j,k)) &
      & + dtdzdx * (-F(i,l,j,k) + F(i,l,j+1,k)) &
      & + dtdxdy * (-G(i,l,j,k) + G(i,l,j,k+1))
      Q2(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - coef * R
    enddo
  end subroutine calc_step1


  !> CUDA Fortran kernel for 1st~3rd step of 4-4 Runge-Kutta
  attributes(global) subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                 :: nx                   !< number of grid points in x direction
    integer, intent(in), value                 :: ny                   !< number of grid points in y direction
    integer, intent(in), value                 :: nz                   !< number of grid points in z direction
    real(8), intent(in), value                 :: coef1                !< coefficient for Runge-Kutta
    real(8), intent(in), value                 :: coef2                !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous    :: dx(nx-1)             !< grid size in x direction
    real(8), intent(in), device, contiguous    :: dy(ny-1)             !< grid size in y direction
    real(8), intent(in), device, contiguous    :: dz(nz-1)             !< grid size in z direction
    real(8), intent(in), device, contiguous    :: E(nx-1,5,ny-2,nz-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(nx-2,5,ny-1,nz-2)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(nx-2,5,ny-2,nz-1)  !< Flux in z direction
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)        !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous   :: Q2(nx,5,ny,nz)       !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(inout), device, contiguous :: Rs(nx-2,5,ny-2,nz-2) !< accumulation for 4-4 Runge-Kutta
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== Compute Grid Jacobian Volumes ==========
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1))
    dtdxdy = dt * 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1))
    ! ========== Conservative Update via 4-4 RK: Stage 1-3 ==========
    ! For stage 1-3: Q^(s) = Q^(s-1) - coef1 * dt/vol * Flux_div + accumulate in Rs
    ! coef2 applies weighting to residual for final 4th stage assembly
    do l = 1, 5
      R = dtdydz * (-E(i,l,j,k) + E(i+1,l,j,k)) &
      & + dtdzdx * (-F(i,l,j,k) + F(i,l,j+1,k)) &
      & + dtdxdy * (-G(i,l,j,k) + G(i,l,j,k+1))
      Q2(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - coef1 * R ! Intermediate Q for next stage
      Rs(i,l,j,k) = Rs(i,l,j,k) + coef2 * R            ! Accumulate weighted residual
    enddo
  end subroutine calc_step
 

  !> CUDA Fortran kernel for 2nd & 3rd step of 3-3 TVD Runge-Kutta
  !> TVD RK3 Stage 2 & 3: Q^(n+1) = (α*Q^n + β*Q^(*) - γ*R)/(α+β)
  attributes(global) subroutine calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Qin, Qout)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), value                 :: coef1               !< α coefficient (weight of original Q^n)
    real(8), intent(in), value                 :: coef2               !< β coefficient (weight of Q^(*))
    real(8), intent(in), value                 :: coef3               !< γ coefficient (weight of flux residual)
    real(8), intent(in), value                 :: coef4               !< 1/(α+β) normalization factor
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< grid size in x direction
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< grid size in y direction
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< grid size in z direction
    real(8), intent(in), device, contiguous    :: E(nx-1,5,ny-2,nz-2) !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(nx-2,5,ny-1,nz-2) !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(nx-2,5,ny-2,nz-1) !< Flux in z direction
    real(8), intent(in), device, contiguous    :: Qin(nx,5,ny,nz)     !< Q^n (original from previous step)
    real(8), intent(inout), device, contiguous :: Qout(nx,5,ny,nz)    !< Q^(*) on input, Q^(n+1) on output
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== Compute Grid Jacobian Volumes ==========
    ! Averaged grid spacing for flux divergence calculation
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1)) ! dt·Δy·Δz for x-flux
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1)) ! dt·Δz·Δx for y-flux
    dtdxdy = dt * 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1)) ! dt·Δx·Δy for z-flux
    ! ========== TVD RK3 Stage 2 & 3 Update ==========
    ! Q^(n+1) = (α·Q^n + β·Q^(*) - γ·dt/vol·∇·F) / (α+β)
    ! Stage 2: α=3/4, β=1/4 (from Q^n and Q^(1))
    ! Stage 3: α=1/3, β=2/3 (from Q^n and Q^(2)), then multiply by 3 (coef4 = 1/3)
    do l = 1, 5  ! All conserved variables
      R = dtdydz * (-E(i,l,j,k) + E(i+1,l,j,k)) &
      & + dtdzdx * (-F(i,l,j,k) + F(i,l,j+1,k)) &
      & + dtdxdy * (-G(i,l,j,k) + G(i,l,j,k+1))
      ! Convex combination: weighted average of Qin and Qout minus scaled residual
      Qout(i+1,l,j+1,k+1) = (coef1 * Qin(i+1,l,j+1,k+1) + coef2 * Qout(i+1,l,j+1,k+1) - coef3 * R) / coef4
    enddo
  end subroutine calc_step2_3
  
 
  !> CUDA Fortran kernel for 4th step of 4-4 Runge-Kutta
  !> Final RK4 Stage: Q^n+1 = Q^n - (1/6)·∑(R_ᵢ) where R_ᵢ indexed over 4 stages
  attributes(global) subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
    integer, intent(in), value                 :: nx                   !< number of grid points in x direction
    integer, intent(in), value                 :: ny                   !< number of grid points in y direction
    integer, intent(in), value                 :: nz                   !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)             !< grid size in x direction
    real(8), intent(in), device, contiguous    :: dy(ny-1)             !< grid size in y direction
    real(8), intent(in), device, contiguous    :: dz(nz-1)             !< grid size in z direction
    real(8), intent(in), device, contiguous    :: E(nx-1,5,ny-2,nz-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(nx-2,5,ny-1,nz-2)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(nx-2,5,ny-2,nz-1)  !< Flux in z direction
    real(8), intent(inout), device, contiguous :: Rs(nx-2,5,ny-2,nz-2) !< accumulated residuals from stages 1-3
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz)        !< Q^n on input, Q^n+1 on output
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== Compute Grid Jacobian Volumes ==========
    ! Cell volume scaling consistent with stages 1-3
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1))
    dtdxdy = dt * 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1))
    ! ========== 4-4 RK Final Assembly ==========
    ! Compute 4th stage residual and accumulate with previous stages
    ! Final update: Q^(n+1) = Q^n - (one_sixth) * (R1 + 2*R2 + 2*R3 + R4)
    ! one_sixth ≈ 1/6 is the standard RK4 weight
    do l = 1, 5 ! All conserved variables
      ! Compute 4th stage flux divergence
      R = dtdydz * (-E(i,l,j,k) + E(i+1,l,j,k)) &
      & + dtdzdx * (-F(i,l,j,k) + F(i,l,j+1,k)) &
      & + dtdxdy * (-G(i,l,j,k) + G(i,l,j,k+1))
      ! Accumulate 4th stage residual (not multiplied by coefficient yet)
      Rs(i,l,j,k) = Rs(i,l,j,k) + R
      ! Apply full RK4 update with (1/6) weighting to final solution
      Q(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - Rs(i,l,j,k) * one_sixth
      ! Clear residual accumulator for next time step
      Rs(i,l,j,k) = 0.d0
    enddo
  end subroutine calc_step4
end module calc_steps
