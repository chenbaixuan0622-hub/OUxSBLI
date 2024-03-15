module mod_globals
  use cudafor
  implicit none
  integer, parameter :: accuracy = 2 ! only 2nd-order accuracy is available
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1
  !!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP !
  !           ! 2  SLAU !
  !!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_scheme = 2

  ! mesh
  real(8), parameter :: Lx = 3.d0 
  real(8), parameter :: Ly = 1.d0
  integer, parameter :: nx = 769
  integer, parameter :: ny = 257
  
  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,1)
  type(dim3) :: blocks = dim3((nx-accuracy)/5,(ny-accuracy)/5,1)
  type(dim3) :: threadsE = dim3(32,5,1)
  type(dim3) :: threadsF = dim3(13,32,1)
  type(dim3) :: threads = dim3(13,5,1)
 
  ! time
  integer, parameter :: nt = 1000
  integer, parameter :: np = 10
  real(8), parameter :: dt = 0.0001d0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0

  ! initial condition
  real(8), parameter :: T = 300.d0
  real(8), parameter :: rho0 = 1.4d0
  real(8), parameter :: p0 = 1.d0
  real(8), parameter :: u0 = 3.d0
end module mod_globals

