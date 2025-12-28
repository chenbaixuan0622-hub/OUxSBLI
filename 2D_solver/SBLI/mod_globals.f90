module mod_globals
  use cudafor
  implicit none
  integer, parameter    :: dimension = 2
  integer(4), parameter :: id_visc   = 1
  integer(2), parameter :: id_LL     = 0
  integer(2), parameter :: id_igr    = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! kind2 Euler       !
  !               ! kind4 NS          !
  !               ! kind8 LES         !
  !               ! 1 2nd             !
  !               ! 2 4th             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar         !
  !               ! 1 SMS             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_av         ! 0 no              !
  !               ! 1 Neumann         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! integer(2)  KEEP    !
  !             ! real(2)     SLAU    !
  !             ! real(4)   Weighted  !
  !             ! real(8)   Threshold !
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
  ! slau_wall   ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8), parameter         :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=4), parameter :: id_tvd      = 0
  integer(kind=4), parameter :: id_slau     = 0
  real(8), parameter         :: blt         = 1.d-3

  ! shock + boundary layer
  real(8), parameter :: Lx = 100.d0 * blt
  real(8), parameter :: Ly = 10.d0 * blt
  real(8), parameter :: Lz = 0.d0
  integer, parameter :: nx = 2241!1121
  integer, parameter :: ny = 577!321
  integer, parameter :: nz = 1

  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks
  
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
  integer, parameter         :: start_rescale = 10
  real(8), parameter :: endT  = 1.d-3
  integer, parameter :: np    = 20!200
  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.2d0
  real(8), parameter :: T0    = 171.31d0
  real(8), parameter :: p0    = 14924.d0
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: dt    = 2.d-9
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = 0.89d0
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
  ! oblique shock
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
  real(8), parameter :: ux    = sqrt(u2**2 + v2**2) * dcos(theta)
  real(8), parameter :: uy    = - sqrt(u2**2 + v2**2) * dsin(theta)
  ! reflected shock
  real(8), parameter :: a2    = sqrt(gamma * R * T2)
  real(8), parameter :: M2    = sqrt(u2**2 + v2**2) / a2
  real(8), parameter :: beta_r = dacos(-1.d0) * 40.6d0 / 180.d0
  real(8), parameter :: Mr    = M2 * dsin(beta_r)
  real(8), parameter :: Mr2   = Mr**2
  real(8), parameter :: T3    = T2 * (1.d0 + 2.d0 * (gamma - 1.d0) * (Mr2 - 1.d0) * (1.d0 + gamma * Mr2) / (Mr2 * (gamma + 1.d0)**2))
  real(8), parameter :: p3    = p2 * (1.d0 + 2.d0 * gamma * (Mr2 - 1.d0) / (gamma + 1.d0))
  real(8), parameter :: rho3  = p3 / (R * T3)
  real(8), parameter :: un2   = ux * dsin(beta_r) - uy * dcos(beta_r)
  real(8), parameter :: ut2   = ux * dcos(beta_r) + uy * dsin(beta_r)
  real(8), parameter :: un3   = un2 - 2.d0 * a2 * (Mr - 1.d0 / Mr) / (gamma + 1.d0)
  real(8), parameter :: ut3   = ut2
  real(8), parameter :: ux3   = un3 * dsin(beta_r) + ut3 * dcos(beta_r)
  real(8), parameter :: uy3   =-un3 * dcos(beta_r) + ut3 * dsin(beta_r)
end module mod_globals

