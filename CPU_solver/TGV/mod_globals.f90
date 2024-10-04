module mod_globals
  implicit none
  integer, parameter :: dimension = 3
  integer, parameter :: accuracy  = 2 
  integer, parameter :: offset    = accuracy / 2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_keep     ! kind2 KEEP          !
  !             ! kind4 KEEPPE        !
  !             ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_keep     = 0

  ! mesh
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  real(8), parameter :: Lz = 2.d0 * pi
  ! 4th-order accuracy
  integer, parameter :: nx = 66
  integer, parameter :: ny = 66
  integer, parameter :: nz = 66

  integer(kind=2), parameter  :: id_recal = 0
  integer, parameter          :: step_offset = 0
  integer, parameter          :: nt = 200!1
  integer, parameter          :: np = 200
  real(8), parameter          :: dt = 0.01d0!0.02d0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr = 0.71d0
  real(8), parameter :: Prt = 0.9d0
  real(8), parameter :: R = 287.03d0

  ! initial condition
  real(8), parameter :: M0 = 0.4d0
  real(8), parameter :: RHO0 = 1.d0
end module mod_globals

