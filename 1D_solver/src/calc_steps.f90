module calc_steps
  use cudafor
  use mod_globals, only : dt, dx
  use libm
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_R(nx, i, dt_coef, E, R)
    integer, intent(in), value              :: nx
    integer, intent(in), value              :: i
    real(8), intent(in), value              :: dt_coef       !< coef * dt / dx
    real(8), intent(in), device, contiguous :: E(3,nx-1)      !< flux
    real(8), intent(out)                    :: R(3)
    real(8) v0, v1
    v0 = E(1,i); v1 = E(1,i+1); R(1) = dt_coef * (v1 - v0)
    v0 = E(2,i); v1 = E(2,i+1); R(2) = dt_coef * (v1 - v0)
    v0 = E(3,i); v1 = E(3,i+1); R(3) = dt_coef * (v1 - v0)
  end subroutine calc_R


  !> CUDA Fortran kernel for 1st step of 3-3 TVD Runge-Kutta
  attributes(global) subroutine calc_step1(nx, coef, E, Q, Q2)
    integer, intent(in), value               :: nx
    real(8), intent(in), value               :: coef            !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous  :: E(3,nx-1)        !< flux
    real(8), intent(in), device, contiguous  :: Q(nx,3)          !< present Q(rho, rhou, E)
    real(8), intent(out), device, contiguous :: Q2(nx,3)         !< next    Q(rho, rhou, E)
    real(8) R(3), coef_dt_dx
    integer i, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (nx-2 < i) return
    coef_dt_dx = coef * dt / dx
    call calc_R(nx, i, coef_dt_dx, E, R)
    do l = 1, 3
      ! Q2 is write-once here and only consumed by the next kernel launch: __stcs
      call __stcs(Q2(i+1,l), Q(i+1,l) - R(l))
    enddo
  end subroutine calc_step1


  !> CUDA Fortran kernel for 2nd & 3rd step of 3-3 TVD Runge-Kutta
  !> Q^(n+1) = (coef1*Qin + coef2*Qout - R) * coef4_inv
  attributes(global) subroutine calc_step2_3(nx, coef1, coef2, coef3, coef4_inv, E, Qin, Qout)
    integer, intent(in), value                 :: nx
    real(8), intent(in), value                 :: coef1           !< weight of Qin
    real(8), intent(in), value                 :: coef2           !< weight of Qout
    real(8), intent(in), value                 :: coef3           !< weight of flux residual
    real(8), intent(in), value                 :: coef4_inv       !< 1/(coef1+coef2) normalization
    real(8), intent(in), device, contiguous    :: E(3,nx-1)       !< flux
    real(8), intent(in), device, contiguous    :: Qin(nx,3)       !< Q^n (original from previous step)
    real(8), intent(inout), device, contiguous :: Qout(nx,3)      !< Q^(*) on input, Q^(n+1) on output
    real(8) R(3), coef3_dt_dx
    integer i, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (nx-2 < i) return
    coef3_dt_dx = coef3 * dt / dx
    call calc_R(nx, i, coef3_dt_dx, E, R)
    do l = 1, 3
      ! Qout is write-once here, not re-read until the next kernel launch: __stcs
      block
        real(8) qin_val, qout_val, val
        qin_val  = Qin(i+1,l)
        qout_val = Qout(i+1,l)
        val      = fma(coef1, qin_val, fma(coef2, qout_val, -R(l))) * coef4_inv
        call __stcs(Qout(i+1,l), val)
      end block
    enddo
  end subroutine calc_step2_3
end module calc_steps
