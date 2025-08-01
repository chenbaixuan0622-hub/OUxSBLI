module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 2
  integer, parameter :: accuracy  = 2 
  integer, parameter :: offset    = accuracy / 2
  integer(4), parameter :: id_visc   = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! kind2 Euler       !
  !               ! kind4 NS          !
  !               ! 1 2nd             !
  !               ! 2 4th             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! integer(2)  KEEP    !
  !             ! real(2)     SLAU    !
  !             ! real(4)   Weighted  !
  !             ! real(8)   Threshold !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_sensor   ! 1 Ducros            !
  !             ! 2 Albada            !
  !             ! 3 Ducros + Albada   !
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
  !             ! kind4 HRSLAU2       !
  !             ! kind8 VHRSLAU2      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(2), parameter      :: id_scheme   = 1
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=2), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: slau_wall   = 0

  ! mesh
  integer, parameter :: nx = 1026
  integer, parameter :: ny = 7
  integer, parameter :: nz = 1

  ! GPU
  type(dim3) :: blocksE   = dim3((nx-accuracy+1)/205,(ny-accuracy)/1,1)
  type(dim3) :: blocksF   = dim3((nx-accuracy)/128,(ny-accuracy+1)/1,1)
  type(dim3) :: blocksEv  = dim3((nx-accuracy+1)/205,(ny-accuracy)/1,1)
  type(dim3) :: blocksFv  = dim3((nx-accuracy)/128,(ny-accuracy+1)/1,1)
  type(dim3) :: blocks    = dim3((nx-accuracy)/128,(ny-accuracy)/1,1)
  type(dim3) :: threadsE  = dim3(205,1,1)
  type(dim3) :: threadsF  = dim3(128,1,1)
  type(dim3) :: threadsEv = dim3(205,1,1)
  type(dim3) :: threadsFv = dim3(128,1,1)
  type(dim3) :: threads   = dim3(128,1,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !               ! kind=8 ! Gauss   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_recal      = 0
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: np            = 100

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

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
  real(8), parameter :: dt   = CFL * Lx / (dble(Nx-1) * sqrt(p0 / rho0))
  real(8), parameter :: endT = 0.2136d0 * Lx / sqrt(p0 / rho0)
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

