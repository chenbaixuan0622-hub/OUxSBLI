!> Interior-only Hybrid kernel variants — skip boundary stencil-order fallback branches.
!> Called when id_bc_x/y/z = .false. (periodic domain) so all points are interior.
!> Eliminates the if/elseif/else order-fallback branches that cause warp divergence
!> in the regular calc_hybrid_kernel.f90 subroutines.
!> Shared memory layout matches the regular kernels (2D shaped, sweep direction first).
!> The KEEP/SLAU physics branch (fdx <= threshold) is preserved — only the stencil-order
!> boundary fallback is removed.
module calc_hybrid_kernel_internal
  use mod_globals, only : id_accuracy, id_slau, gamma, threshold, threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1, R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  use calc_muscl
  use calc_hybrid
  implicit none
  private
  public calc_hybrid_x_in, calc_hybrid_y_in, calc_hybrid_z_in
  !> io = 0, 1, 2 for 2nd, 4th, 6th order — compile-time stencil half-width.
  !> Derived from id_accuracy kind: kind=2 → io=0, kind=4 → io=1, kind=8 → io=2.
  integer, parameter :: io = kind(id_accuracy) / 3
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0

  interface KEEP
    module procedure KEEP2, KEEP4, KEEP6
  end interface KEEP

  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU

  !> Generic delta interpolation dispatched by id_accuracy kind (SLAU path).
  interface delta_r
    module procedure delta_r2, delta_r4, delta_r6
  end interface delta_r
contains
  include 'calc_keep_3d.f90'
  include 'calc_slau_3d.f90'

  !$dir inline
  attributes(device) subroutine delta_r2(id_acc, sensor, a, al, ar)
    integer(2), intent(in), value :: id_acc   !< dispatch key: kind=2 → 1st-order
    real(8), intent(in), value    :: sensor
    real(8), intent(in)           :: a(2)     !< two-point stencil [i, i+1]
    real(8), intent(out)          :: al, ar
    al = a(1); ar = a(2)
  end subroutine delta_r2

  !$dir inline
  attributes(device) subroutine delta_r4(id_acc, sensor, a, al, ar)
    integer(4), intent(in), value :: id_acc   !< dispatch key: kind=4 → MUSCL 3rd
    real(8), intent(in), value    :: sensor
    real(8), intent(in)           :: a(4)     !< four-point stencil [i-1:i+2]
    real(8), intent(out)          :: al, ar
    call delta4(sensor, a, al, ar)
  end subroutine delta_r4

  !$dir inline
  attributes(device) subroutine delta_r6(id_acc, sensor, a, al, ar)
    integer(8), intent(in), value :: id_acc   !< dispatch key: kind=8 → MUSCL 4th
    real(8), intent(in), value    :: sensor
    real(8), intent(in)           :: a(6)     !< six-point stencil [i-2:i+3]
    real(8), intent(out)          :: al, ar
    call delta6(sensor, a, al, ar)
  end subroutine delta_r6


  !> Interior-only Hybrid kernel for x-direction convective flux.
  !> Computes E at interfaces io+1..nx-(io+1) only (no boundary fallback branching).
  !> Shared memory: 2D shaped (x-sweep first), size -(io-1):threadsE%x+io+1.
  attributes(global) subroutine calc_hybrid_x_in(nx, ny, nz, Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer, intent(in), value                :: nx                  !< grid points x
    integer, intent(in), value                :: ny                  !< grid points y
    integer, intent(in), value                :: nz                  !< grid points z
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< temperature (for KEEP path)
    real(8), intent(in), device, contiguous   :: sensor(nx,ny,nz)    !< Ducros shock sensor
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2) !< x-direction flux
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(-(io-1):threadsE%x+io+1, threadsE%y, threadsE%z), shared :: rho, u, v, w, p
    real(8), dimension(threadsE%x, threadsE%y, threadsE%z), shared :: rhor, ur, vr, wr, pr
    real(8) fdx, rhol, ul, vl, wl, pl
    it = threadIdx%x;  jt = threadIdx%y;  kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    ! Phase 1: load shared memory tile (including halo of width io)
    do ii = it-io, threadsE%x+io+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(i,1,j,k);   u(ii,jt,kt) = Q(i,2,j,k)
          v(ii,jt,kt) = Q(i,3,j,k);   w(ii,jt,kt) = Q(i,4,j,k)
          p(ii,jt,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    i = i_base + it
    ! Phase 2: interior-only flux (no stencil-order boundary fallback)
    if (io+1 <= i .and. i <= nx-(io+1) .and. j <= ny-1 .and. k <= nz-1) then
      fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
      if (fdx <= threshold) then
        ! KEEP: entropy/energy preserving flux — reads T from global memory
        block
          real(8) tmp(2*io+2)
          tmp = T(i-io:i+io+1, j, k)
          associate(uu => u)
            E(i,:,j-1,k-1) = KEEP(id_accuracy, &
                                   rho(it-io:it+io+1,jt,kt), u(it-io:it+io+1,jt,kt), &
                                     v(it-io:it+io+1,jt,kt), w(it-io:it+io+1,jt,kt), &
                                    uu(it-io:it+io+1,jt,kt), p(it-io:it+io+1,jt,kt), tmp, Normal_x)
          end associate
        end block
      else
        ! SLAU path: reconstruct left/right states via MUSCL interpolation
        call delta_r(id_accuracy, fdx, rho(it-io:it+io+1,jt,kt), rhol, rhor(it,jt,kt))
        call delta_r(id_accuracy, fdx,   u(it-io:it+io+1,jt,kt),   ul,   ur(it,jt,kt))
        call delta_r(id_accuracy, fdx,   v(it-io:it+io+1,jt,kt),   vl,   vr(it,jt,kt))
        call delta_r(id_accuracy, fdx,   w(it-io:it+io+1,jt,kt),   wl,   wr(it,jt,kt))
        call delta_r(id_accuracy, fdx,   p(it-io:it+io+1,jt,kt),   pl,   pr(it,jt,kt))
      endif
    endif
    ! Sync before left-state overwrite: ensures all delta_r reads of shared memory
    ! are complete before any thread overwrites its shared memory slot with the left state.
    call syncthreads()
    ! Phase 3: overwrite shared memory with left state, then compute SLAU flux
    if (io+1 <= i .and. i <= nx-(io+1) .and. j <= ny-1 .and. k <= nz-1 .and. fdx > threshold) then
      rho(it,jt,kt) = rhol;   u(it,jt,kt) = ul;   v(it,jt,kt) = vl
        w(it,jt,kt) = wl;     p(it,jt,kt) = pl
      associate(un1 => u(it,jt,kt), un2 => ur(it,jt,kt))
        call SLAU(id_slau, rho(it,jt,kt), rhor(it,jt,kt), u(it,jt,kt), ur(it,jt,kt), v(it,jt,kt), vr(it,jt,kt), &
                  w(it,jt,kt), wr(it,jt,kt), un1, un2, p(it,jt,kt), pr(it,jt,kt), Normal_x, 1.d0, &
                  E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1))
      end associate
    endif
  end subroutine calc_hybrid_x_in


  !> Interior-only Hybrid kernel for y-direction convective flux.
  !> Shared memory: 2D shaped (y-sweep first: jj,it,kt).
  attributes(global) subroutine calc_hybrid_y_in(nx, ny, nz, Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer, intent(in), value                :: nx                  !< grid points x
    integer, intent(in), value                :: ny                  !< grid points y
    integer, intent(in), value                :: nz                  !< grid points z
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< temperature (for KEEP path)
    real(8), intent(in), device, contiguous   :: sensor(nx,ny,nz)    !< Ducros shock sensor
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2) !< y-direction flux
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(-(io-1):threadsF%y+io+1, threadsF%x, threadsF%z), shared :: rho, u, v, w, p
    real(8), dimension(threadsF%y, threadsF%x, threadsF%z), shared :: rhor, ur, vr, wr, pr
    real(8) fdy, rhol, ul, vl, wl, pl
    it = threadIdx%x;  jt = threadIdx%y;  kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    ! Phase 1: load shared memory tile
    do jj = jt-io, threadsF%y+io+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(i,1,j,k);   u(jj,it,kt) = Q(i,2,j,k)
          v(jj,it,kt) = Q(i,3,j,k);   w(jj,it,kt) = Q(i,4,j,k)
          p(jj,it,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    j = j_base + jt
    ! Phase 2: interior-only flux
    if (io+1 <= j .and. j <= ny-(io+1) .and. i <= nx-1 .and. k <= nz-1) then
      fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
      if (fdy <= threshold) then
        block
          real(8) tmp(2*io+2)
          tmp = T(i, j-io:j+io+1, k)
          associate(vv => v)
            F(i-1,:,j,k-1) = KEEP(id_accuracy, &
                                   rho(jt-io:jt+io+1,it,kt), u(jt-io:jt+io+1,it,kt), &
                                     v(jt-io:jt+io+1,it,kt), w(jt-io:jt+io+1,it,kt), &
                                    vv(jt-io:jt+io+1,it,kt), p(jt-io:jt+io+1,it,kt), tmp, Normal_y)
          end associate
        end block
      else
        call delta_r(id_accuracy, fdy, rho(jt-io:jt+io+1,it,kt), rhol, rhor(jt,it,kt))
        call delta_r(id_accuracy, fdy,   u(jt-io:jt+io+1,it,kt),   ul,   ur(jt,it,kt))
        call delta_r(id_accuracy, fdy,   v(jt-io:jt+io+1,it,kt),   vl,   vr(jt,it,kt))
        call delta_r(id_accuracy, fdy,   w(jt-io:jt+io+1,it,kt),   wl,   wr(jt,it,kt))
        call delta_r(id_accuracy, fdy,   p(jt-io:jt+io+1,it,kt),   pl,   pr(jt,it,kt))
      endif
    endif
    call syncthreads()
    if (io+1 <= j .and. j <= ny-(io+1) .and. i <= nx-1 .and. k <= nz-1 .and. fdy > threshold) then
      rho(jt,it,kt) = rhol;   u(jt,it,kt) = ul;   v(jt,it,kt) = vl
        w(jt,it,kt) = wl;     p(jt,it,kt) = pl
      associate(un1 => v(jt,it,kt), un2 => vr(jt,it,kt))
        call SLAU(id_slau, rho(jt,it,kt), rhor(jt,it,kt), u(jt,it,kt), ur(jt,it,kt), v(jt,it,kt), vr(jt,it,kt), &
                  w(jt,it,kt), wr(jt,it,kt), un1, un2, p(jt,it,kt), pr(jt,it,kt), Normal_y, 1.d0, &
                  F(i-1,1,j,k-1), F(i-1,2,j,k-1), F(i-1,3,j,k-1), F(i-1,4,j,k-1), F(i-1,5,j,k-1))
      end associate
    endif
  end subroutine calc_hybrid_y_in


  !> Interior-only Hybrid kernel for z-direction convective flux.
  !> Shared memory: 2D shaped (z-sweep first: kk,jt,it).
  attributes(global) subroutine calc_hybrid_z_in(nx, ny, nz, Q, T, sensor, G)
    use mod_constant, only : Normal_z
    integer, intent(in), value                :: nx                  !< grid points x
    integer, intent(in), value                :: ny                  !< grid points y
    integer, intent(in), value                :: nz                  !< grid points z
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< temperature (for KEEP path)
    real(8), intent(in), device, contiguous   :: sensor(nx,ny,nz)    !< Ducros shock sensor
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1) !< z-direction flux
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(-(io-1):threadsG%z+io+1, threadsG%y, threadsG%x), shared :: rho, u, v, w, p
    real(8), dimension(threadsG%z, threadsG%y, threadsG%x), shared :: rhor, ur, vr, wr, pr
    real(8) fdz, rhol, ul, vl, wl, pl
    it = threadIdx%x;  jt = threadIdx%y;  kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    ! Phase 1: load shared memory tile
    do kk = kt-io, threadsG%z+io+1, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(i,1,j,k);   u(kk,jt,it) = Q(i,2,j,k)
          v(kk,jt,it) = Q(i,3,j,k);   w(kk,jt,it) = Q(i,4,j,k)
          p(kk,jt,it) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    k = k_base + kt
    ! Phase 2: interior-only flux
    if (io+1 <= k .and. k <= nz-(io+1) .and. i <= nx-1 .and. j <= ny-1) then
      fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
      if (fdz <= threshold) then
        block
          real(8) tmp(2*io+2)
          tmp = T(i, j, k-io:k+io+1)
          associate(ww => w)
            G(i-1,:,j-1,k) = KEEP(id_accuracy, &
                                   rho(kt-io:kt+io+1,jt,it), u(kt-io:kt+io+1,jt,it), &
                                     v(kt-io:kt+io+1,jt,it), w(kt-io:kt+io+1,jt,it), &
                                    ww(kt-io:kt+io+1,jt,it), p(kt-io:kt+io+1,jt,it), tmp, Normal_z)
          end associate
        end block
      else
        call delta_r(id_accuracy, fdz, rho(kt-io:kt+io+1,jt,it), rhol, rhor(kt,jt,it))
        call delta_r(id_accuracy, fdz,   u(kt-io:kt+io+1,jt,it),   ul,   ur(kt,jt,it))
        call delta_r(id_accuracy, fdz,   v(kt-io:kt+io+1,jt,it),   vl,   vr(kt,jt,it))
        call delta_r(id_accuracy, fdz,   w(kt-io:kt+io+1,jt,it),   wl,   wr(kt,jt,it))
        call delta_r(id_accuracy, fdz,   p(kt-io:kt+io+1,jt,it),   pl,   pr(kt,jt,it))
      endif
    endif
    call syncthreads()
    if (io+1 <= k .and. k <= nz-(io+1) .and. i <= nx-1 .and. j <= ny-1 .and. fdz > threshold) then
      rho(kt,jt,it) = rhol;   u(kt,jt,it) = ul;   v(kt,jt,it) = vl
        w(kt,jt,it) = wl;     p(kt,jt,it) = pl
      associate(un1 => w(kt,jt,it), un2 => wr(kt,jt,it))
        call SLAU(id_slau, rho(kt,jt,it), rhor(kt,jt,it), u(kt,jt,it), ur(kt,jt,it), v(kt,jt,it), vr(kt,jt,it), &
                  w(kt,jt,it), wr(kt,jt,it), un1, un2, p(kt,jt,it), pr(kt,jt,it), Normal_z, 1.d0, &
                  G(i-1,1,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,4,j-1,k), G(i-1,5,j-1,k))
      end associate
    endif
  end subroutine calc_hybrid_z_in
end module calc_hybrid_kernel_internal
