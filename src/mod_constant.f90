module mod_constant
  use cudafor
  use mod_globals, only : R, gamma, Pr
  implicit none
  real(8), parameter :: gamma_1      = gamma - 1.d0
  real(8), parameter :: over_gamma_1 = 1.d0 / (gamma - 1.d0)
  real(8), parameter :: Cp = gamma * R * over_gamma_1
  real(8), parameter :: Cp_over_Pr = Cp / Pr
  real(8), constant  :: Normal_x(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
  real(8), constant  :: Normal_y(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
  real(8), constant  :: Normal_z(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
end module mod_constant

