module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy  = 2 
  integer, parameter :: offset    = accuracy / 2
  integer, parameter :: id_visc   = 0
  integer, parameter :: id_turbulence = 0
  integer, parameter :: id_av     = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! 0 no-visc         !
  !               ! 1 visc 2nd        !
  !               ! 2 visc 4th        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_av         ! 0 no              !
  !               ! 1 Neumann         !  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar         !
  !               ! 1 SMS             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! 1  KEEP             !
  !             ! 2  KEEP MUSCL       !
  !             ! 3  SLAU             !
  !             ! 4  KEEPUP           !
  !             ! 5  Hybrid Weighted  !
  !             ! 6  Hybrid threshold !
  !             ! 7  Hybrid Sigmoid   !
  !             ! 8  KEEP + Roe       !
  !             ! 9  KEEP2nd          !
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
  !             ! kind8 4thMUSCL      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_keep     ! kind2 KEEP          !
  !             ! kind4 KEEPPE        !
  !             ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU          !
  !             ! kind4 HRSLAU2       !
  !             ! kind8 VHRSLAU2      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: id_scheme    = 3
  integer, parameter         :: id_sensor    = 1
  real(8), parameter         :: threshold    = 0.9d0
  integer(kind=8), parameter :: id_accuracy  = 0
  integer(kind=4), parameter :: id_tvd       = 0
  integer(kind=2), parameter :: id_keep      = 0
  integer(kind=4), parameter :: id_slau      = 0

  ! mesh
  real(8), parameter :: Lx = 1.d0
  real(8), parameter :: Ly = 0.1d0 * Lx
  real(8), parameter :: Lz = 0.1d0 * Lx
  integer, parameter :: nx = 129
  integer, parameter :: ny = 7
  integer, parameter :: nz = 7

  integer, parameter :: nre = nx-4

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/128,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/127,(ny-accuracy+1)/1,(nz-accuracy)/1)
  type(dim3) :: blocksG = dim3((nx-accuracy)/127,(ny-accuracy)/1,(nz-accuracy+1)/1)
  type(dim3) :: blocks  = dim3((nx-accuracy)/127,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: threadsE = dim3(128,1,1)
  type(dim3) :: threadsF = dim3(127,1,1)
  type(dim3) :: threadsG = dim3(127,1,1)
  type(dim3) :: threads  = dim3(127,1,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !               ! kind=8 ! 10step  !
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
  real(8), parameter :: rho0 = 1.d0
  real(8), parameter :: p0   = 1.d0
  real(8), parameter :: rho1 = 0.125d0
  real(8), parameter :: p1   = 0.1d00
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lx / (dble(Nx-1) * sqrt(p0 / rho0))
  real(8), parameter :: endT = 0.2d0
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

