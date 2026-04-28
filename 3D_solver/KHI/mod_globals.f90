module mod_globals
  use cudafor
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! kind2 Euler     !
  !               ! kind4 NS        !
  !               ! kind8 LES       !
  !               ! 1 2nd           !
  !               ! 2 4th           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! integer(2) KEEP   !
  !             ! real(2)    SLAU   !
  !             ! real(4)    Roe    !
  !             ! real(8)    Hybrid !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_accuracy ! kind2 2nd         !
  !             ! kind4 4th         !
  !             ! kind8 6th         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd      ! kind2 non TVD     !
  !             ! kind4 minmod      !
  !             ! kind8 MUSCL4th    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU        !
  !             ! kind4 HR-SLAU2    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off         !
  !             ! kind4 on          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: dimension   = 3
  integer(4), parameter      :: id_visc     = 2
  integer(2), parameter      :: id_scheme   = 0
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.9_sp
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: id_rescale  = 0
  integer(kind=2), parameter :: id_gpumpi   = 0
  real(8), parameter         :: blt         = 0.d0

  ! mesh
  real(8), parameter :: Lx = 5.d-3
  real(8), parameter :: Ly = 5.d-3
  real(8), parameter :: Lz = 5.d-3
  integer, parameter :: nx = 130
  integer, parameter :: ny = 130
  integer, parameter :: nz = 130
  
  ! boundary condition
  logical, parameter :: id_bc_x = .false.
  logical, parameter :: id_bc_y = .false.
  logical, parameter :: id_bc_z = .false.

  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
  integer, parameter :: rerank = 0

  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

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
  integer, parameter          :: step_offset = 0
  integer, parameter          :: start_rescale = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.15d0

  ! initial condition
  real(8), parameter :: p    = 2.72d3
  ! M2
  real(8), parameter :: T1   = 162.7908d0
  real(8), parameter :: rho1 = p / (R * T1)
  real(8), parameter :: u1   = 2.d0 * sqrt(gamma * R * T1)
  ! M3
  real(8), parameter :: T2   = 104.6303d0
  real(8), parameter :: rho2 = p / (R * T2)
  real(8), parameter :: u2   = 3.d0 * sqrt(gamma * R * T2)
  
  real(8), parameter :: amp  = 0.05d0
  real(8), parameter :: CFL  = 0.05d0
  real(8), parameter :: dt   = CFL * Lx / (dble(nx-1) * abs(u1))
  real(8), parameter :: endT = 10.d0 * Lx / u2
  integer, parameter :: np   = 100
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

