module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 2
  integer, parameter :: accuracy  = 2 
  integer, parameter :: offset    = accuracy / 2
  integer(4), parameter :: id_visc   = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc     ! kind2 Euler         !
  !             ! kind4 NS            !
  !             ! 1 2nd               !
  !             ! 2 4th               !
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
  integer(2), parameter      :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=2), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: slau_wall   = 0

  integer(kind=2), parameter :: id_forcing  = 0
  ! mesh
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  real(8), parameter :: Lz = 0.d0
  integer, parameter :: nx = 129!257!513
  integer, parameter :: ny = 129!257!513
  integer, parameter :: nz = 1

  type(dim3) :: blocksE   = dim3((nx-accuracy+1)/64,(ny-accuracy)/1,1)
  type(dim3) :: blocksF   = dim3((nx-accuracy)/1,(ny-accuracy+1)/64,1)
  type(dim3) :: blocksEv  = dim3((nx-accuracy+1)/64,(ny-accuracy)/1,1)
  type(dim3) :: blocksFv  = dim3((nx-accuracy)/1,(ny-accuracy+1)/64,1)
  type(dim3) :: blocks    = dim3((nx-accuracy)/73,(ny-accuracy)/1,1)
  type(dim3) :: threadsE  = dim3(64,1,1)
  type(dim3) :: threadsF  = dim3(1,64,1)
  type(dim3) :: threadsEv = dim3(64,1,1)
  type(dim3) :: threadsFv = dim3(1,64,1)
  type(dim3) :: threads   = dim3(73,1,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_recal = 0
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer, parameter         :: step_offset   = 0

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
  real(8), parameter :: Re   = 1d3
  real(8), parameter :: T    = 273.2d0 
  real(8), parameter :: S    = 111.d0
  real(8), parameter :: mu0  = 1.716d-5 * (273.2d0 + S) / (T + S) * (T / 273.2d0)**1.5d0
  real(8), parameter :: u0   = M0 * sqrt(gamma * R * T)
  real(8), parameter :: rho0 = mu0 * Re / (u0 * pi)
  real(8), parameter :: d1   = pi / 15.d0
  real(8), parameter :: d2   = 0.05d0
  real(8), parameter :: dt   = 0.25d0 * 0.25d0 * 1.d-3 / u0
  real(8), parameter :: endT = 8.d0 / u0
  integer, parameter :: np   = 10
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

