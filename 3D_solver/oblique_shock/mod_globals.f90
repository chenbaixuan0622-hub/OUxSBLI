module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy = 2 
  integer(kind=2**(accuracy/2)), parameter :: id_accuracy = 1
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 0 
  !!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP !
  !           ! 2  Roe  !
  !           ! 3  SLAU !
  !!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_hybrid = 0
  integer(kind=4), parameter :: id_muscl = 0
  integer, parameter :: id_scheme = 3

  ! mesh
  real(8), parameter :: Lx = 48d-3
  real(8), parameter :: Ly = 16d-3
  real(8), parameter :: Lz = 8d-3
  integer, parameter :: nx = 513
  integer, parameter :: ny = 321
  integer, parameter :: nz = 161
  real(8), parameter :: dx = Lx / (nx-1)
  real(8), parameter :: dy = Ly / (ny-1)
  real(8), parameter :: dz = Lz / (nz-1)
  real(8), parameter :: dxi = 1.d0 / dx
  real(8), parameter :: dyi = 1.d0 / dy
  real(8), parameter :: dzi = 1.d0 / dz

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/1,(nz-accuracy)/3)
  type(dim3) :: blocksF = dim3((nx-accuracy)/7,(ny-accuracy+1)/32,(nz-accuracy)/3)
  type(dim3) :: blocksG = dim3((nx-accuracy)/7,(ny-accuracy)/1,(nz-accuracy+1)/32)
  type(dim3) :: threadsE = dim3(32,1,3)
  type(dim3) :: threadsF = dim3(7,32,3)
  type(dim3) :: threadsG = dim3(7,1,32)

  ! time
  integer, parameter :: nt = 250
  integer, parameter :: np = 40
  real(8), parameter :: dt = 1.5d-7

  real(8), parameter :: dtdx = dt / dx
  real(8), parameter :: dtdy = dt / dy
  real(8), parameter :: dtdz = dt / dz

  ! physical properties
  real(8), parameter :: gamma = 1.4d0

  ! MUSCL
  real(8), parameter :: k = -1.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)

  ! initial condition
  real(8), parameter :: T = 171.31d0

  ! variables
  real(8), save :: Q(nx,ny,nz,5)
  real(8), save :: T0(nx,ny,nz)
end module mod_globals

