module mod_globals
  use cudafor
  implicit none
  integer, parameter                       :: dimension = 3
  integer, parameter                       :: accuracy = 2 
  integer(kind=2**(accuracy/2)), parameter :: id_accuracy = 1
  integer, parameter                       :: offset = accuracy / 2
  integer, parameter                       :: id_visc = 1 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! 0 no-visc               !
  !               ! 1 visc                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar               !
  !               ! 1 Smagorinsky           !
  !               ! 2 selective_mixed_scale !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_hybrid ! kind2 off   !
  !           ! kind4 on    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_muscl  ! kind2 no    !
  !           ! kind4 3rd   !
  !           ! kind8 4th   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP     !
  !           ! 2  Roe      !
  !           ! 3  SLAU     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau   ! kind2 slau  !
  !           ! kind4 sd    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=4), parameter :: id_hybrid = 0
  integer(kind=8), parameter :: id_muscl = 0
  integer, parameter         :: id_scheme = 3
  integer(kind=2), parameter :: id_slau = 0
  real(8), parameter         :: dp_max = 0.d0

  ! mesh
  real(8), parameter :: Lx = 48d-3
  real(8), parameter :: Ly = 16d-3
  real(8), parameter :: Lz = 4d-3
  integer, parameter :: nx = 513
  integer, parameter :: ny = 321
  integer, parameter :: nz = 65
  real(8), parameter :: dz = Lz / dble(nz-1)
  real(8), parameter :: dzi = 1.d0 / dz

  ! GPU
  type(dim3) :: blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/11,(nz-accuracy)/1)
  type(dim3) :: blocksF = dim3((nx-accuracy)/1,(ny-accuracy+1)/32,(nz-accuracy)/3)
  type(dim3) :: blocksG = dim3((nx-accuracy)/1,(ny-accuracy)/11,(nz-accuracy+1)/32)
  type(dim3) :: blocks = dim3((nx-accuracy)/1,(ny-accuracy)/11,(nz-accuracy)/3)
  type(dim3) :: threadsE = dim3(32,11,1)
  type(dim3) :: threadsF = dim3(1,32,3)
  type(dim3) :: threadsG = dim3(1,11,32)
  type(dim3) :: threads = dim3(1,11,3)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind2 ! 3rd_TVD !
  !               ! kind4 ! 4th     !
  !               ! kind8 ! 10step  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_RungeKutta = 0
  integer, parameter :: nt = 500
  integer, parameter :: np = 40
  real(8), parameter :: u0 = 506.8d0
  real(8), parameter :: dt = 0.1d0 * (Lx / dble(nx-1)) / u0

  real(8), parameter :: dtdz = dt / dz

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr = 0.71d0
  real(8), parameter :: Prt = 0.9d0

  ! MUSCL
  real(8), parameter :: k = 1.d0 / 3.d0
  real(8), parameter :: b = (3.d0 - k) / (1.d0 - k)
  real(8), parameter :: omega = 4.d0
  real(8), parameter :: sigma = 2.d0
  real(8), parameter :: eps = 1.d0

  ! initial condition
  real(8), parameter :: R = 287.03d0
  real(8), parameter :: M0 = 1.9d0
  real(8), parameter :: p0 = 14924.d0
  real(8), parameter :: T0 = 171.31d0
  real(8), parameter :: beta = dacos(-1.d0) * 37.2d0 / 180.d0
  real(8), parameter :: Ms = M0 * dsin(beta)
  real(8), parameter :: Ms2 = Ms**2
  real(8), parameter :: theta = datan(2.d0 * (1.d0 / dtan(beta)) * (Ms2 - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
  real(8), parameter :: rho0 = p0 / (R * T0)
  real(8), parameter :: rho2 = rho0 * (gamma + 1.d0) * Ms2 / ((gamma - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
  real(8), parameter :: p2 = p0 * (1.d0 + 2.d0 * gamma * (Ms2 - 1.d0) / (gamma + 1.d0))
  real(8), parameter :: u1 = u0 * dsin(beta)
  real(8), parameter :: v1 = u0 * dcos(beta)
  real(8), parameter :: a1 = u0 / M0
  real(8), parameter :: u2 = u1 - 2.d0 * a1 * (Ms - 1.d0 / Ms) / (gamma + 1.d0)
  real(8), parameter :: v2 = u0 * dcos(beta)
  real(8), parameter :: u_magnitude = sqrt(u2**2 + v2**2)
  real(8), parameter :: ux = u_magnitude * dcos(theta)
  real(8), parameter :: uy = - u_magnitude * dsin(theta)

  ! variables
  real(8), allocatable :: Q(:,:,:,:)
end module mod_globals

