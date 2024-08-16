module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy = 2 
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! 0 no-visc               !
  !               ! 1 visc2nd               !
  !               ! 2 visc4th               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_dim        ! 1 dimensional           !
  !               ! 2 non-dimensional       !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar               !
  !               ! 1 Smagorinsky           !
  !               ! 2 selective_mixed_scale !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP4th            !
  !           ! 2  KEEP MUSCL         !
  !           ! 3  SLAU               !
  !           ! 4  KEEPUP             !
  !           ! 5  Hybrid Weighted    !
  !           ! 6  Hybrid threshold   !
  !           ! 7  Hybrid Sigmoid     !
  !           ! 8  KEEP + Rho         !
  !           ! 9  KEEP2nd            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_sensor ! 1 Ducros              !
  !           ! 2 Albada              !
  !           ! 3 Ducros + Albada     !
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
  integer(kind=2), parameter :: id_slau     = 0
  integer(kind=2), parameter :: id_rescale  = 0

  ! mesh
  real(8), parameter :: L0 = 1.524d-3
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi * L0
  real(8), parameter :: Ly = 2.d0 * pi * L0
  real(8), parameter :: Lz = 2.d0 * pi * L0
  integer, parameter :: nx = 130!66
  integer, parameter :: ny = 130!66
  integer, parameter :: nz = 130!66

  integer, parameter :: nre = nx-4

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/3,(ny-accuracy)/8,(nz-accuracy)/8)
  type(dim3) :: blocksF = dim3((nx-accuracy)/8,(ny-accuracy+1)/3,(nz-accuracy)/8)
  type(dim3) :: blocksG = dim3((nx-accuracy)/8,(ny-accuracy)/8,(nz-accuracy+1)/3)
  type(dim3) :: blocks  = dim3((nx-accuracy)/8,(ny-accuracy)/8,(nz-accuracy)/8)
  type(dim3) :: threadsE = dim3(3,8,8)
  type(dim3) :: threadsF = dim3(8,3,8)
  type(dim3) :: threadsG = dim3(8,8,3)
  type(dim3) :: threads  = dim3(8,8,8)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_recal = 0
  integer(kind=4), parameter :: id_RungeKutta = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! MUSCL
  real(8), parameter :: k     = 1.d0 / 3.d0
  real(8), parameter :: b     = (3.d0 - k) / (1.d0 - k)
  real(8), parameter :: eps   = 1.d0
  real(8), parameter :: omega = 4.d0
  real(8), parameter :: sigma = 2.d0

  ! initial condition
  real(8), parameter :: Re   = 1600.d0
  real(8), parameter :: M0   = 0.1d0
  real(8), parameter :: T    = 530.d0 * 5.d0 / 9.d0 
  real(8), parameter :: S    = 111.d0
  real(8), parameter :: mu0  = 1.716d-5 * (273.2d0 + S) / (T + S) * (T / 273.2d0)**1.5d0
  real(8), parameter :: V0   = M0 * sqrt(gamma * R * T)
  real(8), parameter :: RHO0 = mu0 * Re / (V0 * L0)
  real(8), parameter :: p0   = RHO0 * R * T

  real(8), parameter :: CFL = 0.03d0
  real(8), parameter :: dt  = CFL * (Lx / dble(nx-1)) / V0
  real(8), parameter :: dtn = V0 * dt / L0
  integer, parameter :: np  = 100
  integer, parameter :: nt  = int(20.d0 / (dble(np) * dtn))
end module mod_globals

