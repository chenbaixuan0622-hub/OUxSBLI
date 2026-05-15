module mod_globals
  use cudafor
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc     ! kind2 Euler       !
  !             ! kind4 NS          !
  !             ! kind8 LES         !
  !             ! 1 2nd             !
  !             ! 2 4th             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! integer(2) KEEP   !
  !             ! real(2)    SLAU   !
  !             ! real(4)    Roe    !
  !             ! real(8)    Hybrid !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_accuracy ! kind2 2nd         !
  !             ! kind4 4th         !
  !             ! kind8 6th         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd      ! kind2 non TVD     !
  !             ! kind4 minmod      !
  !             ! kind8 MUSCL4th    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU        !
  !             ! kind4 HR-SLAU2    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off         !
  !             ! kind4 on          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: dimension   = 2
  integer(2), parameter      :: id_visc     = 1
  integer(2), parameter      :: id_scheme   = 0
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.4_sp
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=4), parameter :: id_rescale  = 0
  integer(kind=2), parameter :: id_gpumpi   = 0
  real(8), parameter         :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 20.d0 * blt
  real(8), parameter :: Ly = 5.d0 * blt
  ! DNS
  integer, parameter :: nx = 257!513
  integer, parameter :: ny = 129!161
  
  ! boundary condition
  logical, parameter :: id_bc_x = .true.
  logical, parameter :: id_bc_y = .true.

  ! RTX 4090
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsEv = dim3(64,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksEv, blocksFv, blocks

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind2 ! 3rd_TVD !
  !               ! kind4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind2 ! set 0   !
  !               ! kind4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_RungeKutta = 0
  integer(kind=2), parameter :: id_recal      = 0
  integer, parameter         :: step_offset   = 0
  real(8), parameter :: endT  = 0.1d-3
  integer, parameter :: np    = 10
  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.d0
  real(8), parameter :: T0    = 171.31d0
  real(8), parameter :: p0    = 14924.d0
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: dt    = 3.d-9
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = 0.89d0
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
  ! oblique shock
  real(8), parameter :: beta  = dacos(-1.d0) * 40.03d0 / 180.d0
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

