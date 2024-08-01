module mod_globals
  use cudafor
  implicit none
  integer, parameter                       :: dimension = 3
  integer, parameter                       :: accuracy = 2 
  integer, parameter                       :: offset = accuracy / 2
  integer, parameter                       :: id_visc = 1 
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
  integer, parameter         :: id_scheme   = 1
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=4), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=4), parameter :: id_rescale  = 0

  ! mesh
  real(8), parameter :: Lx = 20d-3
  real(8), parameter :: Ly = 8d-3
  real(8), parameter :: Lz = 4d-3
  integer, parameter :: nx = 513!257
  integer, parameter :: ny = 257
  integer, parameter :: nz = 129!65

  integer, parameter :: nre = int(0.5 * nx)

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,(nz-accuracy)/1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/1,(ny-accuracy+1)/128,(nz-accuracy)/1)
  type(dim3) :: blocksG = dim3((nx-accuracy)/1,(ny-accuracy)/5,(nz-accuracy+1)/32)
  type(dim3) :: blocks  = dim3((nx-accuracy)/1,(ny-accuracy)/51,(nz-accuracy)/1)
  type(dim3) :: threadsE = dim3(32,5,1)
  type(dim3) :: threadsF = dim3(1,128,1)
  type(dim3) :: threadsG = dim3(1,5,32)
  type(dim3) :: threads  = dim3(1,51,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind2 ! 3rd_TVD !
  !               ! kind4 ! 4th     !
  !               ! kind8 ! 10step  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind2 ! set 0   !
  !               ! kind4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer(kind=2), parameter :: id_recal = 0
  integer, parameter :: nt  = 50
  integer, parameter :: np  = 2000
  real(8), parameter :: u0  = 506.8d0
  real(8), parameter :: CFL = 0.1d0
  real(8), parameter :: dt  = CFL * Lx / (dble(nx-1) * u0)

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

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
  real(8), parameter :: rho0 = p0 / (R * T0)
end module mod_globals

