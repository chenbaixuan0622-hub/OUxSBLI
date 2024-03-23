module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 2
  integer, parameter :: accuracy = 2 ! only 2nd-order accuracy is available
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1
  !!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP !
  !           ! 2  Roe  !
  !           ! 3  SLAU !
  !!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_hybrid = 0
  integer, parameter :: id_muscl = 1
  integer, parameter :: id_scheme = 2

  ! mesh
  real(8), parameter :: Lx = 1.d0
  real(8), parameter :: Ly = 1.d0
  real(8), parameter :: Lz = 0.d0
  integer, parameter :: nx = 257
  integer, parameter :: ny = 257
  integer, parameter :: nz = accuracy+1
  real(8), parameter :: dx = Lx / (nx-1)
  real(8), parameter :: dy = Ly / (ny-1)
  real(8), parameter :: dz = Lz / (nz-1)
  real(8), parameter :: dxi = 1.d0 / dx
  real(8), parameter :: dyi = 1.d0 / dy

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,1)
  type(dim3) :: threadsE = dim3(32,5,1)
  type(dim3) :: threadsF = dim3(5,32,1)

  ! time
  integer, parameter :: nt = 200
  integer, parameter :: np = 10
  real(8), parameter :: dt = 0.0001d0

  real(8), parameter :: dtdx = dt / dx
  real(8), parameter :: dtdy = dt / dy
  real(8), parameter :: dtdz = 0.d0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0

  ! MUSCL
  real(8), parameter :: k = -1.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)

  ! initial condition
  real(8), parameter :: T = 300.d0
  real(8), parameter :: rhol = 1.d0
  real(8), parameter :: rhor = 0.125d0
  real(8), parameter :: pl = 1.d0
  real(8), parameter :: pr = 0.1d0

  ! variables
  real(8), save :: Q(nx,ny,4)
  real(8), save :: T0(nx,ny)
end module mod_globals

