module mod_globals
  use cudafor
  implicit none
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
  real(8), parameter :: Lx      = 0.001d0
  integer, parameter :: nx      = 4097
  real(8), parameter :: dx      = Lx / dble(nx-1)
  real(8), parameter :: over_dx = 1.d0 / dx

  ! GPU
  type(dim3), parameter :: threads = dim3(128,1,1)
  type(dim3)            :: blocks  = dim3((nx-1+threads%x -1)/threads%x, 1, 1)

  ! time
  integer, parameter :: np = 1

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: R     = 287.03d0

  ! Re & Ma
  real(8), parameter :: Re   = 250.d0
  real(8), parameter :: T    = 273.2d0 
  real(8), parameter :: S    = 111.d0
  real(8), parameter :: mu0  = 1.716d-5 * (273.2d0 + S) / (T + S) * (T / 273.2d0)**1.5d0
  real(8), parameter :: Ms   = 1.2d0
  ! shock-upstream
  real(8), parameter :: uu_s = Ms * sqrt(gamma * R * T)
  real(8), parameter :: rhou = mu0 * Re / (uu_s * Lx)
  real(8), parameter :: pu   = rhou * R * T
  ! shock-downstream
  real(8), parameter :: pd   =   pu * (1.d0 + 2.d0 * (Ms**2 - 1.d0) * gamma / (gamma + 1.d0))
  real(8), parameter :: ud_s = uu_s * (gamma + 1.d0 + (gamma - 1.d0) * pd / pu) / &
                                      (gamma - 1.d0 + (gamma + 1.d0) * pd / pu)
  real(8), parameter :: rhod = rhou * uu_s / ud_s
  ! Galilean
  real(8), parameter :: uu   = uu_s !+ 0.25d0 * uu_s
  real(8), parameter :: ud   = ud_s !+ 0.25d0 * uu_s

  real(8), parameter :: CFL  = 0.01d0
  real(8), parameter :: dt   = CFL * Lx / (dble(nx-1) * ud)
  real(8), parameter :: endT = 0.1d0 * Lx / ud
  real(8), parameter :: nt   = endT / (dble(np) * dt)
  real(8), parameter :: dtdx = dt / dx
end module mod_globals

