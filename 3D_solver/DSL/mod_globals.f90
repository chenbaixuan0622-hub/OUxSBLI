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
  !               ! 1 visc            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar         !
  !               ! 1 SMSky           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_av         ! 0 non             !
  !               ! 1 vonNeumann      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! 1  KEEP             !
  !             ! 2  KEEP MUSCL       !
  !             ! 3  SLAU             !
  !             ! 4  KEEPUP           !
  !             ! 5  Hybrid Weighted  !
  !             ! 6  Hybrid threshold !
  !             ! 7  Hybrid Sigmoid   !
  !             ! 8  KEEP + Rho       !
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
  !             ! kind8 MUSCL4th      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_keep     ! kind2 KEEP          !
  !             ! kind4 KEEPPE        !
  !             ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU          !
  !             ! kind4 HR-SLAU2      !
  !             ! kind8 VHR-SLAU2     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: id_scheme   = 3
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.9d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=4), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: id_rescale  = 0

  ! mesh
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  real(8), parameter :: Lz = 0.1 * Lx
  integer, parameter :: nx = 129
  integer, parameter :: ny = 129
  integer, parameter :: nz = 7
  integer, parameter :: nre = nx-4

  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/1,(ny-accuracy+1)/32,(nz-accuracy)/1)
  type(dim3) :: blocksG = dim3((nx-accuracy)/127,(ny-accuracy)/1,(nz-accuracy+1)/1)
  type(dim3) :: blocks  = dim3((nx-accuracy)/127,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: threadsE = dim3(32,1,1)
  type(dim3) :: threadsF = dim3(1,32,1)
  type(dim3) :: threadsG = dim3(127,1,1)
  type(dim3) :: threads  = dim3(127,1,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter  :: id_recal = 0
  integer(kind=4), parameter  :: id_RungeKutta = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! HR-SLAU2 and HR-AUSM+-up towards High Resolution Unsteady Aerodynamic Simulations
  ! Keiichi Kitamura, Atsushi Hashimoto
  ! JAXA-SP-14-010

  ! initial condition
  real(8), parameter :: M0   = 0.1d0
  real(8), parameter :: rho0 = 1.d0
  real(8), parameter :: u0   = 1.d0
  real(8), parameter :: d1   = pi / 15.d0
  real(8), parameter :: d2   = 0.05d0
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = 1.d-3
  real(8), parameter :: endT = 8.d0
  integer, parameter :: np   = 100
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

