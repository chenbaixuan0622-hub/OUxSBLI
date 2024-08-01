module calc_visc
  use mod_globals, only : accuracy, offset, gamma, R, Pr
  use calc_Sutherland
  implicit none
contains
  attributes(global) subroutine calc_Ev(nx, dx, rho, u, T, p, E)
    integer, intent(in), value                   :: nx
    real(8), intent(in), dimension(nx-1), device :: dx ! 1 / dx
    real(8), intent(in), dimension(nx), device   :: rho, u, T, p
    real(8), intent(inout), device               :: E(nx-accuracy+1,3)
    integer i
    real(8) :: mu = 0.d0
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) um, ux, T1, T2, Tx, txx, kappa
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1

    if (2 <= i .and. i <= nx-2) then
      T1 = 0.5d0 * (T(i-1) + T(i))
      T2 = 0.5d0 * (T(i+1) + T(i+2))
      call calc_mu(T1,T2,mu)
      txx = 4.d0 * dx(i) * mu * (u(i-1) - 27.d0 * u(i) + 27.d0 * u(i+1) - u(i+2)) / 72.d0
      um = 0.0625d0 * (- u(i-1) + 9.d0 * u(i) + 9.d0 * u(i+1) - u(i+2))
      Tx = (T(i-1) - 27.d0 * T(i) + 27.d0 * T(i+1) - T(i+2)) * dx(i) / 24.d0
    else
      call calc_mu(T(i),T(i+1),mu)
      txx = 4.d0 * dx(i) * mu * (-u(i) + u(i+1)) / 3.d0
      um = 0.5d0 * (u(i) + u(i+1))
      Tx = (-T(i) + T(i+1)) * dx(i)
    endif

    kappa = Cp * mu / Pr
    E(i-offset+1,2) = E(i-offset+1,2) - txx
    E(i-offset+1,3) = E(i-offset+1,3) - (txx * um + kappa * Tx)
  end subroutine calc_Ev
end module calc_visc

