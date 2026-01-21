module mod_globals
  use cudafor
  integer, parameter    :: dimension     = 1
  integer(2), parameter :: id_LL         = 0
  integer(4), parameter :: id_igr        = 0
  integer(4), parameter :: id_visc       = 0
  integer(2), parameter :: id_scheme     = 0
  integer(8), parameter :: id_accuracy   = 0
  integer(8), parameter :: id_tvd        = 0
  integer(4), parameter :: id_slau       = 0
  integer(4), parameter :: id_rescale    = 0
  integer(2), parameter :: id_gpumpi     = 0
  integer(2), parameter :: id_RungeKutta = 0
  integer(2), parameter :: id_recal      = 0

  ! mesh
  integer, parameter :: nx = 4097
  
  ! GPU
  type(dim3), parameter :: threads = dim3(128,1,1)
  type(dim3)            :: blocks  = dim3((nx-1+threads%x -1)/threads%x, 1, 1)

  ! time
  integer, parameter :: np = 1

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: T0   = 300.d0
  real(8), parameter :: C    = 1.461d-6
  real(8), parameter :: S    = 110.3d0
  real(8), parameter :: mu0  = C * T0**1.5 / (T0 + S)
  real(8), parameter :: Re   = 25000.d0
  real(8), parameter :: rho0 = 1.293d0
  real(8), parameter :: p0   = rho0 * R * T0
  real(8), parameter :: rho1 = 0.125d0 * rho0
  real(8), parameter :: p1   = 0.1d0 * p0
  real(8), parameter :: Lx   = Re * mu0 / sqrt(rho0 * p0)
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lx / (dble(nx-1) * sqrt(p0 / rho0))
  real(8), parameter :: endT = 0.2136d0 * Lx / sqrt(p0 / rho0)
  real(8), parameter :: nt   = endT / (dble(np) * dt)
  real(8), parameter :: dx   = Lx / dble(nx-1)
  real(8), parameter :: over_dx = 1.d0 / dx
  real(8), parameter :: dtdx = dt / dx
end module mod_globals

