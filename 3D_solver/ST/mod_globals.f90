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
  ! id_dim        ! 1 dimensional           !
  !               ! 2 non-dimensional       !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar               !
  !               ! 1 Smagorinsky           !
  !               ! 2 selective_mixed_scale !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP4th          !
  !           ! 2  KEEP MUSCL       !
  !           ! 3  SLAU             !
  !           ! 4  KEEPUP           !
  !           ! 5  Hybrid Weighted  !
  !           ! 6  Hybrid threshold !
  !           ! 7  Hybrid Sigmoid   !
  !           ! 8  KEEP + Roe       !
  !           ! 9  KEEP2nd          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_sensor ! 1 Ducros            !
  !           ! 2 Albada            !
  !           ! 3 Ducros + Albada   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd    ! kind2 non TVD       !
  !           ! kind4 minmod        !
  !           ! kind8 4thMUSCL      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_keep   ! kind2 KEEP          !
  !           ! kind4 KEEPPE        !
  !           ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau   ! kind2 SLAU          !
  !           ! kind4 HRSLAU2       !
  !           ! kind8 VHRSLAU2      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: id_scheme = 2
  integer, parameter         :: id_sensor = 3
  integer(kind=2), parameter :: id_tvd  = 0
  integer(kind=4), parameter :: id_keep = 0
  integer(kind=8), parameter :: id_slau = 0

  ! mesh
  integer, parameter :: nx = 514
  integer, parameter :: ny = 5
  integer, parameter :: nz = 5

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/514,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/128,(ny-accuracy+1)/1,(nz-accuracy)/1)
  type(dim3) :: blocksG = dim3((nx-accuracy)/128,(ny-accuracy)/1,(nz-accuracy+1)/1)
  type(dim3) :: blocks  = dim3((nx-accuracy)/128,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: threadsE = dim3(514,1,1)
  type(dim3) :: threadsF = dim3(128,1,1)
  type(dim3) :: threadsG = dim3(128,1,1)
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
  integer(kind=2), parameter :: id_recal = 0
  integer(kind=2), parameter :: id_RungeKutta = 0
  integer, parameter         :: np = 100

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: Prt   = 0.9d0
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
  real(8), parameter :: Ly   = 0.1d0 * Lx
  real(8), parameter :: Lz   = 0.1d0 * Lx
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = cfl * lX / (DBLE(NX-1) * SQRT(P0 / RHO0))
  real(8), parameter :: endT = 0.2136d0 * Lx / sqrt(p0 / rho0)
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

