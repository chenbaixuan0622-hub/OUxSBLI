module mod_globals
  use cudafor
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc     ! kind2 Euler       !
  !             ! kind4 NS          !
  !             ! kind8 LES         !
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
  real(2), parameter    :: id_scheme   = 0
  integer, parameter    :: sp          = 4
  real(sp), parameter   :: threshold   = 0.1_sp
  integer(2), parameter :: id_accuracy = 0
  integer(2), parameter :: id_tvd      = 0
  integer(2), parameter :: id_slau     = 0
  integer(2), parameter :: id_rescale  = 0
  integer(2), parameter :: id_gpumpi   = 0

  ! mesh
  real(8), parameter :: Lx = 2.0d0    ! dummy (main_curv passes it to set_grid; set_grid ignores it)
  real(8), parameter :: Ly = 10.0d0   ! dummy (= far_r, for consistency)
  real(8), parameter :: Lz = 0.5d0    ! quasi-2D spanwise extent
  integer, parameter :: nx = 194      ! 192 interior cells (i=2..193) wrapping airfoil
  integer, parameter :: ny = 80       ! wall-normal cells
  integer, parameter :: nz = 4        ! quasi-2D spanwise

  ! NACA 0012 O-grid geometry
  real(8), parameter :: chord = 1.d0
  real(8), parameter :: aoa   = -acos(-1.d0) * 5.d0 / 180.d0 ! angle of attack [radians]
  real(8), parameter :: far_r = 10.d0                        ! far-field radius [chords]

  real(8), parameter :: beta  = dacos(-1.d0) * 37.2d0 / 180.d0
  ! Boundary condition flags
  logical, parameter :: id_bc_x = .false. ! xi is periodic (O-grid wraps around airfoil)
  logical, parameter :: id_bc_y = .true.  ! eta has airfoil wall (j=1) and far-field (j=ny)
  logical, parameter :: id_bc_z = .false. ! z-periodic

  ! GPU thread block dimensions (tuned for nx≈194, ny≈80)
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
  integer(2), parameter :: id_RungeKutta = 0  ! kind=2 → TVD-RK3
  integer(2), parameter :: id_recal      = 0  ! kind=2 → initialize
  integer, parameter    :: step_offset   = 0

  ! Free-stream flow conditions (non-dimensional)
  real(8), parameter :: Ma_inf  = 0.8d0   ! Mach number (subsonic)
  real(8), parameter :: gamma   = 1.4d0   ! heat capacity ratio
  real(8), parameter :: R       = 1.d0    ! gas constant (non-dimensional)
  real(8), parameter :: Pr      = 0.72d0  ! Prandtl number

  ! Time stepping parameters
  ! dt=1e-3: CFL estimate with stretch=1.05, ny=80, far_r=10 gives Δη_wall≈0.011,
  ! (|u|+a)=1.5 → dt_max≈7e-3; dt=1e-3 is safely below that.
  real(8), parameter :: dt = 2.d-4
  integer, parameter :: nt = 10000
  integer, parameter :: np = 10
  integer, parameter :: rerank = -1
end module mod_globals
