module mod_globals
  use cudafor
  implicit none
  integer, parameter    :: dimension = 3
  integer, parameter    :: accuracy  = 2
  integer, parameter    :: offset    = accuracy / 2
  integer(4), parameter :: id_visc   = 2
  integer, parameter    :: id_av     = 0
  integer, parameter    :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! kind2 Euler       !
  !               ! kind4 NS          !
  !               ! kind8 LES         !
  !               ! 1 2nd             !
  !               ! 2 4th             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar         !
  !               ! 1 SMS             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_av         ! 0 no              !
  !               ! 1 Neumann         !
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
  !             ! kind4 HR-SLAU2      !
  !             ! kind8 VHR-SLAU2     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! slau_wall   ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(2), parameter         :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=4), parameter :: slau_wall   = 0
  integer(kind=4), parameter :: id_rescale  = 0
  real(8), parameter         :: blt         = 2.d-3

  ! mesh
  real(8), parameter :: Lx = 20d-3 ! 10 delta
  real(8), parameter :: Ly = 8d-3  !  4 delta
  real(8), parameter :: Lz = 4d-3  !  2 delta
  ! LES
  !integer, parameter :: nx = 257   ! xp = 10   0.04   mm
  !integer, parameter :: ny = 257   ! yp = 0.5, 0.002  mm
  !integer, parameter :: nz = 129   ! zp = 5    0.0156 mm
  ! DNS
  integer, parameter :: nx = 513 ! xp = 10   0.04   mm
  integer, parameter :: ny = 257 ! yp = 0.5, 0.002  mm
  integer, parameter :: nz = 257 ! zp = 5    0.0156 mm

  integer, parameter :: nre = int(0.8 * nx)

 ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,(nz-accuracy)/1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/1,(ny-accuracy+1)/128,(nz-accuracy)/1)
  type(dim3) :: blocksG = dim3((nx-accuracy)/1,(ny-accuracy)/5,(nz-accuracy+1)/32)
  type(dim3) :: blocks  = dim3((nx-accuracy)/1,(ny-accuracy)/51,(nz-accuracy)/1)
  type(dim3) :: threadsE = dim3(32,5,1)
  type(dim3) :: threadsF = dim3(1,128,1)
  type(dim3) :: threadsG = dim3(1,5,32)
  type(dim3) :: threads  = dim3(1,51,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind2 ! 3rd_TVD !
  !               ! kind4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind2 ! set 0   !
  !               ! kind4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer(kind=2), parameter :: id_recal      = 0
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: start_rescale = 80
  real(8), parameter :: endT  = 0.5d-3
  integer, parameter :: np    = 100
  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: T0    = 171.31d0
  real(8), parameter :: u0    = 506.8d0! + sqrt(gamma * R * T0)
  real(8), parameter :: CFL   = 0.1d0
  real(8), parameter :: dt    = CFL * Lx / (dble(nx-1) * u0)
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0

  ! initial condition
  real(8), parameter :: M0   = 1.9d0
  real(8), parameter :: p0   = 14924.d0
  real(8), parameter :: rho0 = p0 / (R * T0)
end module mod_globals

