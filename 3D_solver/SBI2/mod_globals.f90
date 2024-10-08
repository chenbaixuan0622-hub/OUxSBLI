module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy  = 2
  integer, parameter :: offset    = accuracy / 2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! integer(2)  KEEP    !
  !             ! real(2)     SLAU    !
  !             ! real(4)   Weighted  !
  !             ! real(8)   Threshold !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_sensor   ! 1 Ducros            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_accuracy ! kind2 2nd           !
  !             ! kind4 4th           !
  !             ! kind8 6th           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd      ! kind2 non TVD       !
  !             ! kind4 minmod        !
  !             ! kind8 switch        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_keep     ! kind2 KEEP          !
  !             ! kind4 KEEPPE        !
  !             ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU          !
  !             ! kind4 HR-SLAU2      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! slau_wall   ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(2), parameter         :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.6d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=8), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=4), parameter :: slau_wall   = 0
  integer(kind=4), parameter :: id_rescale  = 0
  integer, parameter         :: start_rescale = 1000

  ! mesh
  real(8), parameter :: Lx = 20d-3 ! * 2 20 delta
  real(8), parameter :: Ly = 8d-3  !  4 delta
  real(8), parameter :: Lz = 4d-3  !  2 delta
  integer, parameter :: nx = 769   ! * 2
  integer, parameter :: ny = 257
  integer, parameter :: nz = 257
  integer, parameter :: rerank = 0
  integer, parameter :: nre = int(nx * 8 / 15)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind2 ! 3rd_TVD !
  !               ! kind4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind2 ! set 0   !
  !               ! kind4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer(kind=2), parameter :: id_recal      = 0
  integer, parameter         :: step_offset   = 0
  real(8), parameter :: endT = 0.5d-3
  integer, parameter :: np   = 20!100
  real(8), parameter :: u0   = 506.8d0
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lx / (dble(nx-1) * u0) ! 7.7e-9
  integer, parameter :: nt   = 1!int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: M0    = 1.9d0
  real(8), parameter :: p0    = 14924.d0
  real(8), parameter :: T0    = 171.31d0
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: beta  = dacos(-1.d0) * 37.2d0 / 180.d0
  real(8), parameter :: Ms    = M0 * dsin(beta)
  real(8), parameter :: Ms2   = Ms**2
  real(8), parameter :: theta = datan(2.d0 * (1.d0 / dtan(beta)) * (Ms2 - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
  real(8), parameter :: T2    = T0 * (1.d0 + 2.d0 * (gamma - 1.d0) * (Ms2 - 1.d0) * (1.d0 + gamma * Ms2) / (Ms2 * (gamma + 1.d0)**2))
  real(8), parameter :: p2    = p0 * (1.d0 + 2.d0 * gamma * (Ms2 - 1.d0) / (gamma + 1.d0))
  real(8), parameter :: rho2  = p2 / (R * T2)
  real(8), parameter :: u1    = u0 * dsin(beta)
  real(8), parameter :: v1    = u0 * dcos(beta)
  real(8), parameter :: a1    = u0 / M0
  real(8), parameter :: u2    = u1 - 2.d0 * a1 * (Ms - 1.d0 / Ms) / (gamma + 1.d0)
  real(8), parameter :: v2    = u0 * dcos(beta)
  real(8), parameter :: u_magnitude = sqrt(u2**2 + v2**2)
  real(8), parameter :: ux    = u_magnitude * dcos(theta)
  real(8), parameter :: uy    = - u_magnitude * dsin(theta)
end module mod_globals

