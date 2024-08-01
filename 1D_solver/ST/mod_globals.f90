module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 1
  integer, parameter :: accuracy = 2 
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc   ! 0 no-visc  !
  !           ! 1 visc     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 0  KEEP    !
  !           ! 1  KEP     !
  !           ! 2  KEEPKEP !
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau   ! kind2 slau !
  !           ! kind4 sd   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: id_scheme = 0
  ! mesh
  integer, parameter :: nx = 514

  ! GPU
  type(dim3) :: blocksE  = dim3((nx-accuracy+1)/27,1,1)
  type(dim3) :: blocks   = dim3((nx-accuracy)/128,1,1)
  type(dim3) :: threadsE = dim3(27,1,1)
  type(dim3) :: threads  = dim3(128,1,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !               ! kind=8 ! 10step  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_RungeKutta = 0
  integer, parameter         :: np = 100

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: R     = 287.03d0

  ! MUSCL
  real(8), parameter :: k     = 1.d0 / 3.d0
  real(8), parameter :: b     = (3.d0 - k) / (1.d0 - k)
  real(8), parameter :: eps   = 1.d0
  real(8), parameter :: omega = 4.d0
  real(8), parameter :: sigma = 2.d0

  ! initial condition
  real(8), parameter :: T0   = 300.d0
  real(8), parameter :: C    = 1.456d-6
  real(8), parameter :: S    = 110.4d0
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
end module mod_globals

