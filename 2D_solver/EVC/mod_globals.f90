module mod_globals
  use cudafor
  implicit none
  integer, parameter    :: dimension = 2
  integer, parameter    :: accuracy  = 2 
  integer, parameter    :: offset    = accuracy / 2
  integer(2), parameter :: id_visc   = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! kind2 Euler       !
  !               ! kind4 NS          !
  !               ! 1 2nd             !
  !               ! 2 4th             !
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
  integer(2), parameter      :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: slau_wall   = 0

  ! mesh
  real(8), parameter :: Lx = 0.1d0
  real(8), parameter :: Ly = 0.1d0
  real(8), parameter :: Lz = 0.d0
  integer, parameter :: nx = 130!258!66
  integer, parameter :: ny = 130!258!66
  integer, parameter :: nz = 1

  type(dim3) :: blocksE   = dim3((nx-accuracy+1)/1,(ny-accuracy)/32,1)
  type(dim3) :: blocksF   = dim3((nx-accuracy)/32,(ny-accuracy+1)/1,1)
  type(dim3) :: blocksEv  = dim3((nx-accuracy+1)/1,(ny-accuracy)/32,1)
  type(dim3) :: blocksFv  = dim3((nx-accuracy)/32,(ny-accuracy+1)/1,1)
  type(dim3) :: blocks    = dim3((nx-accuracy)/8,(ny-accuracy)/8,1)
  type(dim3) :: threadsE  = dim3(1,32,1)
  type(dim3) :: threadsF  = dim3(32,1,1)
  type(dim3) :: threadsEv = dim3(1,32,1)
  type(dim3) :: threadsFv = dim3(32,1,1)
  type(dim3) :: threads   = dim3(8,8,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter  :: id_recal      = 0
  integer(kind=4), parameter  :: id_RungeKutta = 0
  integer, parameter          :: step_offset   = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.15d0

  ! initial condition
  real(8), parameter :: M0    = 0.05d0
  real(8), parameter :: beta  = 1.d0 / 50.d0
  real(8), parameter :: theta = 0.d0 / 180.d0
  real(8), parameter :: Rc    = 0.005d0
  real(8), parameter :: p0    = 1.d5
  real(8), parameter :: T0    = 300.d0
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: CFL   = 0.05d0
  real(8), parameter :: dt    = CFL * Lx / (dble(nx-1) * u0)
  real(8), parameter :: T     = 1.d0 * Lx / u0
  integer, parameter :: np    = 1
  integer, parameter :: nt    = int(T / (dble(np) * dt))
end module mod_globals

