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
  integer(4), parameter      :: id_visc     = 2
  integer(2), parameter      :: id_scheme   = 0
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.4_sp
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=8), parameter :: id_tvd      = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: id_rescale  = 0
  integer(kind=2), parameter :: id_gpumpi   = 0
  real(8), parameter         :: blt         = 0.d0

  ! mesh
  real(8), parameter :: L0 = 1.0d0
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi * L0
  real(8), parameter :: Ly = 2.d0 * pi * L0
  real(8), parameter :: Lz = 2.d0 * pi * L0
  integer, parameter :: nx = 129
  integer, parameter :: ny = 129
  integer, parameter :: nz = 129

  ! boundary condition (all periodic)
  logical, parameter :: id_bc_x = .false.
  logical, parameter :: id_bc_y = .false.
  logical, parameter :: id_bc_z = .false.

  integer, parameter :: nre1  = 1
  integer, parameter :: nre2  = nx
  integer, parameter :: rerank = 0

  ! GPU
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
  integer(kind=2), parameter  :: id_recal      = 0
  integer(kind=2), parameter  :: id_RungeKutta = 0
  integer, parameter          :: step_offset   = 0
  integer, parameter          :: start_rescale = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: Mt   = 0.1d0                         ! turbulent Mach number
  real(8), parameter :: Urms = 100.d0
  real(8), parameter :: T    = Urms**2 / (gamma * R * Mt**2) ! c = sqrt(gamma * R * T)
  real(8), parameter :: S    = 111.d0
  real(8), parameter :: mu0  = 1.716d-5 * (273.2d0 + S) / (T + S) * (T / 273.2d0)**1.5d0
  ! Pope (2000) energy spectrum parameters
  real(8), parameter :: pope_C    = 1.5d0  ! constant
  real(8), parameter :: pope_cL   = 6.78d0 ! constant, large-scale coefficient
  real(8), parameter :: pope_p0   = 2.0d0  ! constant, large-scale exponent (gives E~k^2 at low k)
  real(8), parameter :: pope_beta = 5.2d0  ! constant, dissipation-range coefficient
  real(8), parameter :: pope_ceta = 0.40d0 ! constant, dissipation-range offset constant
  ! integral length scale
  real(8), parameter :: kp        = 4.d0   ! spectrum peak
  real(8), parameter :: pope_L    = sqrt(3.d0 * pope_p0 * pope_cL / 5.d0) / kp
  ! Reynolds numbers (set Re_lambda; Re and pope_eta are derived)
  real(8), parameter :: Re_lambda = 20.d0                           ! target Taylor-scale Re
  real(8), parameter :: Re        = Re_lambda**2 / (15.d0 * pope_L) ! integral Re: ρ₀V₀L₀/μ₀
  real(8), parameter :: pope_eta  = pope_L**0.25d0 * Re**(-0.75d0)  ! Kolmogorov scale
  real(8), parameter :: nu0       = 15.d0 * Urms * pope_L / Re_lambda**2
  real(8), parameter :: RHO0      = mu0 / nu0
  real(8), parameter :: p0        = RHO0 * R * T
  ! Petersen-type solenoidal forcing parameters
  real(8), parameter :: eps_s  = 15.d0 * nu0 * Urms**2 / pope_L**2  ! energy injection rate [m^2/s^3]
  integer, parameter :: kf_min = 1                                  ! forcing band lower bound
  integer, parameter :: kf_max = 2                                  ! forcing band upper bound
  real(8), parameter :: C_T    = 0.1d0                              ! temperature relaxation coefficient
  ! dt is acoustic-limited: c_s = urms / Mt dominates advective speed
  real(8), parameter :: CFL = 0.03d0
  real(8), parameter :: dt  = CFL * (Lx / dble(nx-1)) * Mt / urms
  real(8), parameter :: dtn = urms * dt / L0
  integer, parameter :: np  = 200
  integer, parameter :: nt  = 1000
end module mod_globals
