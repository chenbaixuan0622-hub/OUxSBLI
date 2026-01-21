module calc_visc
  use cudafor
  use mod_globals, only : dt, dx, over_dx
  use mod_constant, only : one_third, Cp_over_Pr
  implicit none
contains
  attributes(global) subroutine calc_Ev(nx, rho, u, p, T, mu, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx), mu(nx)
    real(8), intent(out), device :: E(3,nx-1)
    real(8) mua, txx, utxx, kTx
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    mua  = 0.5d0 * (mu(i) + mu(i+1))
    txx  = 4.d0 * mua * (-u(i) + u(i+1)) * one_third * over_dx
    utxx = 0.5d0 * (u(i) + u(i+1)) * txx
    kTx  = Cp_over_Pr * mua * (-T(i) + T(i+1)) * over_dx
    E(2,i) = E(2,i) - txx
    E(3,i) = E(3,i) - (utxx + kTx)
  end subroutine calc_Ev


  attributes(global) subroutine calc_Ev_LL(nx, rho, u, p, T, mu, Z, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx), mu(nx)
    real(8), intent(in), device  :: Z(2,nx)  
    real(8), intent(out), device :: E(3,nx-1)
    real(8) :: kb = 1.380650d-23
    real(8) mua, txx, utxx, kappa, kTx, s, q
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    mua   = 0.5d0 * (mu(i) + mu(i+1))
    txx   = 4.d0 * mua * (-u(i) + u(i+1)) * one_third * over_dx
    utxx  = 0.5d0 * (u(i) + u(i+1)) * txx
    kappa = Cp_over_Pr * mua
    kTx   = kappa * (-T(i) + T(i+1)) * over_dx
    s     = sqrt(4.d0 * kb * mua * (T(i) + T(i+1)) / (3.d0 * dt * dx)) * 0.5d0 * (Z(1,i) + Z(1,i+1))
    q     = sqrt(kb * kappa * (T(i)**2 + T(i+1)**2) / (dt * dx)) * 0.5d0 * (Z(2,i) + Z(2,i+1))
    E(2,i) = E(2,i) - (txx + sqrt(2.d0) * s)
    E(3,i) = E(3,i) - (utxx + kTx + sqrt(2.d0) * (q + 0.5d0 * (u(i) + u(i+1)) * s))
  end subroutine calc_Ev_LL
end module calc_visc

