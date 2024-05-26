module mod_globals
  use cudafor
  implicit none
  integer, parameter                       :: dimension = 3
  integer, parameter                       :: accuracy = 2 
  integer(kind=2**(accuracy/2)), parameter :: id_accuracy = 1
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_hybrid ! kind2 off   !
  !           ! kind4 on    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_muscl  ! kind2 no    !
  !           ! kind4 3rd   !
  !           ! kind8 4th   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP     !
  !           ! 2  Roe      !
  !           ! 3  SLAU     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau   ! kind2 slau  !
  !           ! kind4 sd    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_hybrid = 0
  integer(kind=2), parameter :: id_muscl = 0
  integer, parameter         :: id_scheme = 1
  integer(kind=2), parameter :: id_slau = 0
  real(8), parameter         :: dp_max = 0.d0

  ! mesh
  real(8), parameter :: Lx = 20d-3
  real(8), parameter :: Ly = 8d-3
  real(8), parameter :: Lz = 4d-3
  integer, parameter :: nx = 257
  integer, parameter :: ny = 257
  integer, parameter :: nz = 65
  real(8), parameter :: dz = Lz / dble(nz-1)
  real(8), parameter :: dzi = 1.d0 / dz

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,(nz-accuracy)/1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/1,(ny-accuracy+1)/32,(nz-accuracy)/3)
  type(dim3) :: blocksG = dim3((nx-accuracy)/1,(ny-accuracy)/5,(nz-accuracy+1)/32)
  type(dim3) :: blocks = dim3((nx-accuracy)/5,(ny-accuracy)/5,(nz-accuracy)/3)
  type(dim3) :: threadsE = dim3(32,5,1)
  type(dim3) :: threadsF = dim3(1,32,3)
  type(dim3) :: threadsG = dim3(1,5,32)
  type(dim3) :: threads = dim3(5,5,3)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind2 ! 3rd_TVD !
  !               ! kind4 ! 4th     !
  !               ! kind8 ! 10step  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind2 ! set 0   !
  !               ! kind4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_RungeKutta = 0
  integer(kind=2), parameter :: id_recal = 0
  integer, parameter :: nt = 2000
  integer, parameter :: np = 100
  real(8), parameter :: u0 = 506.8d0
  real(8), parameter :: dt = 5d-9!0.025d0 * (Lx / dble(nx-1)) / u0

  real(8), parameter :: dtdz = dt / dz

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
  real(8), parameter :: rho0 = p0 / (R * T0)

  ! variables
  real(8), allocatable :: Q(:,:,:,:)
end module mod_globals

