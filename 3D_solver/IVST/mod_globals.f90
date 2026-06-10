module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.4_sp
  real(8), parameter         :: blt         = 2.d-3

  ! mesh
  real(8), parameter :: Lx = 1.d0
  real(8), parameter :: Ly = 0.1d0 * Lx
  real(8), parameter :: Lz = 0.1d0 * Lx
  integer, parameter :: nx = 129
  integer, parameter :: ny = 7
  integer, parameter :: nz = 7
  

  integer, parameter :: nre1 = int(1.d0 * dble(nx) / 7.d0)
  integer, parameter :: nre2 = int(2.d0 * dble(nx) / 7.d0)
  integer, parameter :: rerank = 0

  ! GPU
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(128,1,1)
  type(dim3), parameter :: threadsG  = dim3(128,1,1)
  type(dim3), parameter :: threadsEv = dim3(128,1,1)
  type(dim3), parameter :: threadsFv = dim3(128,1,1)
  type(dim3), parameter :: threadsGv = dim3(128,1,1)
  type(dim3), parameter :: threads   = dim3(128,1,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: start_rescale = 0
  integer, parameter         :: np            = 100

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: rho0 = 1.d0
  real(8), parameter :: p0   = 1.d0
  real(8), parameter :: rho1 = 0.125d0
  real(8), parameter :: p1   = 0.1d00
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lx / (dble(Nx-1) * sqrt(p0 / rho0))
  real(8), parameter :: endT = 0.2d0
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

