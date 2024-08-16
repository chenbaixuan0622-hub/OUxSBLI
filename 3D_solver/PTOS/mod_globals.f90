module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy = 2 
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! 0 no-visc               !
  !               ! 1 visc                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar               !
  !               ! 1 Smagorinsky           !
  !               ! 2 selective_mixed_scale !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! 1  KEEP4th          !
  !             ! 2  KEEP MUSCL       !
  !             ! 3  SLAU             !
  !             ! 4  KEEPUP           !
  !             ! 5  Hybrid Weighted  !
  !             ! 6  Hybrid threshold !
  !             ! 7  Hybrid Sigmoid   !
  !             ! 8  KEEP + Rho       !
  !             ! 9  KEEP2nd          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_sensor   ! 1 Ducros            !
  !             ! 2 Albada            !
  !             ! 3 Ducros + Albada   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_accuracy ! kind2 2nd           !
  !             ! kind4 4th           !
  !             ! kind8 6th           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd      ! kind2 non TVD       !
  !             ! kind4 minmod        !
  !             ! kind8 MUSCL4th      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_keep     ! kind2 KEEP          !
  !             ! kind4 KEEPPE        !
  !             ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU          !
  !             ! kind4 HR-SLAU2      !
  !             ! kind8 VHR-SLAU2     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: id_scheme   = 7
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=8), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=4), parameter :: id_rescale  = 0

  ! mesh
  real(8), parameter :: Lx1 = 48d-3
  real(8), parameter :: Ly1 = 12d-3
  real(8), parameter :: Lx2 = Lx1
  real(8), parameter :: Ly2 = 4d-3
  real(8), parameter :: Lz  = 8d-3
  
  integer, parameter :: nx1 = 513
  integer, parameter :: ny1 = 257
  integer, parameter :: nz1 = 65
  integer, parameter :: nx2 = 513
  integer, parameter :: ny2 = 65
  integer, parameter :: nz2 = 65
  integer, parameter :: nre = int(0.3 * nx1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind2 ! set 0   !
  !               ! kind4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_recal = 0
  integer, parameter :: nt  = 1!500
  integer, parameter :: np  = 1!100
  real(8), parameter :: u0  = 506.8d0
  real(8), parameter :: CFL = 0.1d0
  real(8), parameter :: dt  = CFL * Lx1 / (dble(nx1-1) * u0)

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr = 0.71d0
  real(8), parameter :: Prt = 0.9d0
  real(8), parameter :: R = 287.03d0

  ! MUSCL
  real(8), parameter :: k = 1.d0 / 3.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)
  real(8), parameter :: omega = 4.d0
  real(8), parameter :: sigma = 2.d0
  real(8), parameter :: eps = 1.d0

  ! initial condition
  real(8), parameter :: M0 = 1.9d0
  real(8), parameter :: p0 = 14924.d0
  real(8), parameter :: T0 = 171.31d0
  real(8), parameter :: beta = dacos(-1.d0) * 37.2d0 / 180.d0
  real(8), parameter :: Ms = M0 * dsin(beta)
  real(8), parameter :: Ms2 = Ms**2
  real(8), parameter :: theta = datan(2.d0 * (1.d0 / dtan(beta)) * (Ms2 - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
  real(8), parameter :: rho0 = p0 / (R * T0)
  real(8), parameter :: rho2 = rho0 * (gamma + 1.d0) * Ms2 / ((gamma - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
  real(8), parameter :: p2 = p0 * (1.d0 + 2.d0 * gamma * (Ms2 - 1.d0) / (gamma + 1.d0))
  real(8), parameter :: u1 = u0 * dsin(beta)
  real(8), parameter :: v1 = u0 * dcos(beta)
  real(8), parameter :: a1 = u0 / M0
  real(8), parameter :: u2 = u1 - 2.d0 * a1 * (Ms - 1.d0 / Ms) / (gamma + 1.d0)
  real(8), parameter :: v2 = u0 * dcos(beta)
  real(8), parameter :: u_magnitude = sqrt(u2**2 + v2**2)
  real(8), parameter :: ux = u_magnitude * dcos(theta)
  real(8), parameter :: uy = - u_magnitude * dsin(theta)
end module mod_globals

