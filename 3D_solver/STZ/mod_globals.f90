module mod_globals
  use cudafor
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc     ! kind2 Euler       !
  !             ! kind4 NS          !
  !             ! kind8 LES         !
  !             ! 1 2nd             !
  !             ! 2 4th             !
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
  integer(2), parameter      :: id_visc     = 2   ! kind=2 → Euler
  real(2), parameter         :: id_scheme   = 0   ! real(2) → SLAU
  integer, parameter         :: sp          = kind(1.d0)
  real(sp),   parameter      :: threshold   = 0.4_sp
  integer(kind=8), parameter :: id_accuracy = 0   ! kind=2 → 2nd order (overlap=1)
  integer(kind=4), parameter :: id_tvd      = 0   ! kind=4 → minmod limiter
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=4), parameter :: id_rescale  = 0   ! kind=4 → selects RungeKutta_4th_zdec
  integer(kind=2), parameter :: id_gpumpi   = 0
  real(8),    parameter      :: blt         = 0.d0

  ! mesh — nz is LOCAL per rank; global nz = nranks*(nz-2)
  real(8), parameter :: Lx = 0.1d0
  real(8), parameter :: Ly = 0.1d0
  real(8), parameter :: Lz = 1.0d0
  integer, parameter :: nx = 6    ! small: 6 interior cells
  integer, parameter :: ny = 6
  integer, parameter :: nz = 513 ! 128 interior + 1 ghost each end

  logical, parameter :: id_bc_x = .false.  ! periodic x (uniform in x)
  logical, parameter :: id_bc_y = .false.  ! periodic y (uniform in y)
  logical, parameter :: id_bc_z = .true.   ! non-periodic z (extrapolation)

  integer, parameter :: nre1   = 1
  integer, parameter :: nre2   = nx
  integer, parameter :: rerank = 0

  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time — RK4 z-decomposition
  ! id_RungeKutta kind=8 → zdec; id_rescale kind=4 → selects RungeKutta_4th_zdec
  integer(kind=2), parameter :: id_recal      = 0
  integer(kind=8), parameter :: id_RungeKutta = 0
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: start_rescale = 0

  ! Physical (dimensionless Sod shock tube)
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0  ! unused for Euler, kept for NS switch

  real(8), parameter :: rho_L = 1.0d0          ! left state density
  real(8), parameter :: p_L   = 1.0d0          ! left state pressure
  real(8), parameter :: rho_R = 0.125d0         ! right state density
  real(8), parameter :: p_R   = 0.1d0           ! right state pressure

  real(8), parameter :: CFL  = 0.01d0
  real(8), parameter :: dt   = CFL * Lz / (dble(Nz-1) * sqrt(p_L / rho_L))
  real(8), parameter :: endT = 0.1d0
  integer, parameter :: np   = 50
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals
