module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy = 2 
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1 
  !!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP !
  !           ! 2  Roe  !
  !           ! 3  SLAU !
  !!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_hybrid = 0
  integer, parameter :: id_muscl = 0
  integer, parameter :: id_scheme = 1

  ! mesh
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  real(8), parameter :: Lz = 2.d0 * pi
  integer, parameter :: nx = 257
  integer, parameter :: ny = 257
  integer, parameter :: nz = 257
  real(8), parameter :: dx = Lx / (nx-1)
  real(8), parameter :: dy = Ly / (ny-1)
  real(8), parameter :: dz = Lz / (nz-1)
  real(8), parameter :: dxi = 1.d0 / dx
  real(8), parameter :: dyi = 1.d0 / dy
  real(8), parameter :: dzi = 1.d0 / dz

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,(nz-accuracy)/5)
  type(dim3) :: blocksF = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,(nz-accuracy)/5)
  type(dim3) :: blocksG = dim3((nx-accuracy)/5,(ny-accuracy)/5,(nz-accuracy+1)/32)
  type(dim3) :: threadsE = dim3(32,5,5)
  type(dim3) :: threadsF = dim3(5,32,5)
  type(dim3) :: threadsG = dim3(5,5,32)

  ! time
  integer, parameter :: nt = 1250
  integer, parameter :: np = 80
  real(8), parameter :: dt = 0.01d0

  real(8), parameter :: dtdx = dt / dx
  real(8), parameter :: dtdy = dt / dy
  real(8), parameter :: dtdz = dt / dz

  ! physical properties
  real(8), parameter :: gamma = 1.4d0

  ! MUSCL
  real(8), parameter :: k = -1.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)

  ! initial condition
  real(8), parameter :: T = 300.d0
  real(8), parameter :: RHO0 = 1.d0 
  real(8), parameter :: L0 = 1.d0
  real(8), parameter :: M0 = 0.01d0

  ! variables
  real(8), save :: Q(nx,ny,nz,5)
  real(8), save :: T0(nx,ny,nz)
end module mod_globals

