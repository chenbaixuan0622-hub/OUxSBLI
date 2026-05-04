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
  !             ! kind8 Hybrid      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU        !
  !             ! kind4 HR-SLAU2    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off         !
  !             ! kind4 on          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Value is IGNORED; only the kind matters for dispatch
  integer, parameter    :: dimension   = 3
  integer(2), parameter :: id_visc     = 0
  real(8), parameter    :: id_scheme   = 0
  integer, parameter    :: sp          = 4
  real(sp), parameter   :: threshold   = 0.1_sp
  integer(2), parameter :: id_accuracy = 0
  integer(2), parameter :: id_tvd      = 0
  integer(2), parameter :: id_slau     = 0
  integer(2), parameter :: id_rescale  = 0
  integer(2), parameter :: id_gpumpi   = 0
  
  ! mesh
  real(8), parameter :: Lx = 2.0d0
  real(8), parameter :: Ly = 1.0d0
  real(8), parameter :: Lz = 0.1d0
  integer, parameter :: nx = 192
  integer, parameter :: ny = 80
  integer, parameter :: nz = 4
  ! Compression corner geometry
  real(8), parameter :: theta    = 8.d0 * acos(-1.d0) / 180.d0 ! 8° in radians
  real(8), parameter :: x_corner = 0.5d0                       ! corner location

  ! Boundary condition flags
  logical, parameter :: id_bc_x = .true.  ! xi has inlet/outlet walls (non-periodic)
  logical, parameter :: id_bc_y = .true.  ! eta has slip walls at j=1 and j=ny
  logical, parameter :: id_bc_z = .false. ! z-periodic

  ! GPU thread block dimensions (tuned for nx≈192, ny≈80)
  type(dim3), parameter :: threads   = dim3(32, 8, 1)
  type(dim3), parameter :: threadsE  = dim3(32, 4, 2)
  type(dim3), parameter :: threadsFv = dim3(8,  4, 4)
  type(dim3), parameter :: threadsF  = dim3(8,  16, 2)
  type(dim3), parameter :: threadsEv = dim3(16, 4, 2)
  type(dim3), parameter :: threadsG  = dim3(8,  8, 4)
  type(dim3), parameter :: threadsGv = dim3(8,  8, 4)
  type(dim3) :: blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !               ! kind=8 ! Gauss   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(2), parameter :: id_RungeKutta = 0  ! kind=2 → TVD-RK3, kind=4 → RK4
  integer(2), parameter :: id_recal      = 0  ! kind=2 → initialize, kind=4 → restart from file
  integer, parameter    :: step_offset   = 0

  ! Free-stream flow conditions (non-dimensional)
  real(8), parameter :: Ma_inf  = 2.d0    ! Mach number (supersonic)
  real(8), parameter :: gamma   = 1.4d0   ! heat capacity ratio
  real(8), parameter :: R       = 1.d0    ! gas constant (non-dimensional)
  real(8), parameter :: Pr      = 0.72d0  ! Prandtl number

  ! Time stepping parameters
  real(8), parameter :: dt = 5.d-5
  ! Derived time parameters
  integer, parameter :: nt = 10000
  integer, parameter :: np = 10
  integer, parameter :: rerank = -1    ! no re-scaling by default
end module mod_globals
