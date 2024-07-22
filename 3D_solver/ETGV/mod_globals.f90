module mod_globals
  use cudafor
  implicit none
  integer, parameter                        :: dimension = 3
  integer, parameter                        :: accuracy = 2 
  integer(kind=2**(accuracy/2)), parameter  :: id_accuracy = 1
  integer, parameter                        :: offset = accuracy / 2
  integer, parameter                        :: id_visc = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! 0 no-visc               !
  !               ! 1 visc                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar               !
  !               ! 1 Smagorinsky           !
  !               ! 2 selective_mixed_scale !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP             !
  !           ! 2  KEEP MUSCL       !
  !           ! 3  SLAU             !
  !           ! 4  KEEPUP           !
  !           ! 5  Hybrid Weighted  !
  !           ! 6  Hybrid threshold !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_sensor ! 1 Ducros            !
  !           ! 2 Albada            !
  !           ! 3 Ducros + Albada   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd    ! kind2 non TVD       !
  !           ! kind4 minmod        !
  !           ! kind8 MUSCL4th      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_keep   ! kind2 KEEP          !
  !           ! kind4 KEEPPE        !
  !           ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: id_scheme = 3
  integer, parameter         :: id_sensor = 3
  real(8), parameter         :: threshold = 0.4d0
  integer(kind=4), parameter :: id_tvd = 0
  integer(kind=2), parameter :: id_keep = 0

  ! mesh
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  real(8), parameter :: Lz = 2.d0 * pi
  integer, parameter :: nx = 66
  integer, parameter :: ny = 66
  integer, parameter :: nz = 66

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/5,(ny-accuracy)/8,(nz-accuracy)/8)
  type(dim3) :: blocksF = dim3((nx-accuracy)/8,(ny-accuracy+1)/5,(nz-accuracy)/8)
  type(dim3) :: blocksG = dim3((nx-accuracy)/8,(ny-accuracy)/8,(nz-accuracy+1)/5)
  type(dim3) :: blocks  = dim3((nx-accuracy)/8,(ny-accuracy)/8,(nz-accuracy)/8)
  type(dim3) :: threadsE = dim3(5,8,8)
  type(dim3) :: threadsF = dim3(8,5,8)
  type(dim3) :: threadsG = dim3(8,8,5)
  type(dim3) :: threads  = dim3(8,8,8)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !               ! kind=8 ! 10step  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter  :: id_recal = 0
  integer(kind=4), parameter  :: id_RungeKutta = 0
  integer, parameter          :: nt = 200
  integer, parameter          :: np = 200
  real(8), parameter          :: dt = 0.01d0!0.02d0!0.01d0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr = 0.71d0
  real(8), parameter :: Prt = 0.9d0
  real(8), parameter :: R = 287.03d0

  ! MUSCL
  real(8), parameter :: k = 1.d0 / 3.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)
  real(8), parameter :: eps = 1.d0
  real(8), parameter :: omega = 4.d0
  real(8), parameter :: sigma = 2.d0

  ! initial condition
  real(8), parameter :: M0 = 0.4d0
  real(8), parameter :: RHO0 = 1.d0
end module mod_globals

