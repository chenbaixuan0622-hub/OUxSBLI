module calc_flux
  use cudafor
  use mod_globals, only : R, gamma
  use mod_constant, only : R_over_gamma_1
  implicit none
contains
  attributes(global) subroutine calc_KEEP(nx, rho, u, p, T, sigma, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx)
    real(4), intent(in), device  :: sigma(nx)
    real(8), intent(out), device :: E(3,nx-1)
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    E(1,i) = 0.25d0 * (rho(i) + rho(i+1)) * (u(i) + u(i+1))
    E(2,i) = E(1,i) * 0.5d0 * (u(i) + u(i+1)) + 0.5d0 * (p(i) + p(i+1))
    E(3,i) = 0.5d0 * E(1,i) * u(i) * u(i+1) &
           + E(1,i) * 0.5d0 * (T(i) + T(i+1)) * R_over_gamma_1 &
           + 0.5d0 * (u(i) * p(i+1) + u(i+1) * p(i))
  end subroutine calc_KEEP


  attributes(global) subroutine calc_KEEP_G(nx, rho, u, p, T, sigma, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx)
    real(4), intent(in), device  :: sigma(nx)
    real(8), intent(out), device :: E(3,nx-1)
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    E(1,i) = sqrt(rho(i) * rho(i+1)) * 0.5d0 * (u(i) + u(i+1))
    E(2,i) = E(1,i) * 0.5d0 * (u(i) + u(i+1)) + 0.5d0 * (p(i) + p(i+1))
    E(3,i) = 0.5d0 * E(1,i) * u(i) * u(i+1) &
           + E(1,i) * sqrt(T(i) * T(i+1)) * R_over_gamma_1 &
           + 0.5d0 * (u(i) * p(i+1) + u(i+1) * p(i))
  end subroutine calc_KEEP_G


  attributes(global) subroutine calc_KEEP_LF(nx, rho, u, p, T, sigma, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx)
    real(4), intent(in), device  :: sigma(nx)
    real(8), intent(out), device :: E(3,nx-1)
    real(8) alpha
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    if (-u(i-1) + u(i+1) < 0.d0) then
      alpha = max(abs(u(i)), abs(u(i+1)))
    else
      alpha = 0.d0
    endif
    E(1,i) = sqrt(rho(i) * rho(i+1)) * 0.5d0 * (u(i) + u(i+1))
    E(2,i) = E(1,i) * 0.5d0 * (u(i) + u(i+1)) + 0.5d0 * (p(i) + p(i+1)) 
    E(3,i) = 0.5d0 * E(1,i) * u(i) * u(i+1) &
           + E(1,i) * sqrt(T(i) * T(i+1)) * R_over_gamma_1 &
           + 0.5d0 * (u(i) * p(i+1) + u(i+1) * p(i) &
           - alpha * dble(-sigma(i) + sigma(i+1))) ! Lax-Friedrichs
  end subroutine calc_KEEP_LF
end module calc_flux

