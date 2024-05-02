module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy = 2 
  integer(kind=2**(accuracy/2)), parameter :: id_accuracy = 1
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar               !
  !               ! 1 Smagorinsky           !
  !               ! 2 selective_mixed_scale !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP    !
  !           ! 2  Roe     !
  !           ! 3  SLAU    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau   ! kind2 slau !
  !           ! kind4 sd   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_hybrid = 0
  integer(kind=2), parameter :: id_muscl = 0
  integer, parameter :: id_scheme = 1
  integer(kind=2), parameter :: id_slau = 0
  real(8), parameter :: dp_max = 0.d0

  ! mesh
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  real(8), parameter :: Lz = 2.d0 * pi
  integer, parameter :: nx = 66 
  integer, parameter :: ny = 66
  integer, parameter :: nz = 66
  real(8), parameter :: dx = Lx / dble(nx-1)
  real(8), parameter :: dy = Ly / dble(ny-1)
  real(8), parameter :: dz = Lz / dble(nz-1)
  real(8), parameter :: dxi = 1.d0 / dx
  real(8), parameter :: dyi = 1.d0 / dy
  real(8), parameter :: dzi = 1.d0 / dz
  real(8), dimension(nx,ny,nz) :: x, y, z, xix, etay, Jacobian

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/5,(ny-accuracy)/8,(nz-accuracy)/8)
  type(dim3) :: blocksF = dim3((nx-accuracy)/8,(ny-accuracy+1)/5,(nz-accuracy)/8)
  type(dim3) :: blocksG = dim3((nx-accuracy)/8,(ny-accuracy)/8,(nz-accuracy+1)/5)
  type(dim3) :: blocks = dim3((nx-accuracy)/8,(ny-accuracy)/8,(nz-accuracy)/8)
  type(dim3) :: threadsE = dim3(5,8,8)
  type(dim3) :: threadsF = dim3(8,5,8)
  type(dim3) :: threadsG = dim3(8,8,5)
  type(dim3) :: threads = dim3(8,8,8)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer, parameter :: nt = 200
  integer, parameter :: np = 200
  real(8), parameter :: dt = 0.01d0

  real(8), parameter :: dtdx = dt / dx
  real(8), parameter :: dtdy = dt / dy
  real(8), parameter :: dtdz = dt / dz

  ! physical properties
  real(8), parameter :: gamma = 1.4d0

  ! MUSCL
  real(8), parameter :: k = 1.d0 / 3.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)
  real(8), parameter :: eps = 1.d0
  real(8), parameter :: omega = 4.d0
  real(8), parameter :: sigma = 2.d0

  ! initial condition
  real(8), parameter :: R = 287.03d0
  real(8), parameter :: Pr = 0.72d0
  real(8), parameter :: T = 300.d0
  real(8), parameter :: RHO0 = 1.d0 
  real(8), parameter :: L0 = 1.d0
  real(8), parameter :: M0 = 0.4d0

  ! variables
  real(8), save :: Q(nx,ny,nz,5)
  real(8), save :: T0(nx,ny,nz)
end module mod_globals

