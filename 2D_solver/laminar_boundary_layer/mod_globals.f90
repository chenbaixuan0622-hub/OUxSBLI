module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 2
  integer, parameter :: accuracy = 2 ! only 2nd-order accuracy is available
  integer(kind=2**(accuracy/2)), parameter :: id_accuracy = 1
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_hybrid ! kind2 off !
  !           ! kind4 on  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_muscl  ! kind2 off !
  !           ! kind4 3rd !
  !           ! kind8 4th !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP   !
  !           ! 2  Roe    !
  !           ! 3  SLAU   !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=4), parameter :: id_hybrid = 0.d0
  integer(kind=8), parameter :: id_muscl = 0.d0
  integer, parameter :: id_scheme = 3

  ! mesh
  real(8), parameter :: Lx = 4.d0
  real(8), parameter :: Ly = 0.05d0
  real(8), parameter :: Lz = 0.d0 
  integer, parameter :: nx = 257 
  integer, parameter :: ny = 65
  integer, parameter :: nz = accuracy+1
  real(8), parameter :: dx = Lx / (nx-1)
  real(8), parameter :: dy = Ly / (ny-1)
  real(8), parameter :: dz = Lz / (nz-1)
  real(8), parameter :: dxi = 1.d0 / dx
  real(8), parameter :: dyi = 1.d0 / dy

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/3,1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,1)
  type(dim3) :: blocks = dim3((nx-accuracy)/5,(ny-accuracy)/3,1)
  type(dim3) :: threadsE = dim3(32,3,1)
  type(dim3) :: threadsF = dim3(5,32,1)
  type(dim3) :: threads = dim3(5,3,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_Rungekutta ! kind=2 ! 3rd-TVD !
  !               ! kind=4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer, parameter :: nt = 10
  integer, parameter :: np = 10 
  real(8), parameter :: u0 = 34.7d0
  real(8), parameter :: dt = 0.01d0 * dx/ u0

  real(8), parameter :: dtdx = dt / dx
  real(8), parameter :: dtdy = dt / dy
  real(8), parameter :: dtdz = 0.d0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0

  ! MUSCL
  real(8), parameter :: k = 1.d0 / 3.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)
  real(8), parameter :: omega = 4.d0
  real(8), parameter :: sigma = 2.d0
  real(8), parameter :: eps = 1.d0

  ! laminar boundary layer
  real(8), parameter :: R = 287.03d0
  real(8), parameter :: p0 = 1013d2
  real(8), parameter :: T = 300.d0
  real(8), parameter :: rho0 = p0 / (R * T)

  ! variables
  real(8), save :: Q(nx,ny,4)
  real(8), save :: T0(nx,ny)
end module mod_globals

