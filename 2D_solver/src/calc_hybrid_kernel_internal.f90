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
    real(sp), intent(in), value   :: sensor
    real(8), intent(in)           :: a(2)     !< two-point stencil [i, i+1]
    real(8), intent(out)          :: al, ar
    al = a(1); ar = a(2)
  end subroutine delta_r2

  !$dir inline
  attributes(device) subroutine delta_r4(id_acc, sensor, a, al, ar)
    integer(4), intent(in), value :: id_acc   !< dispatch key: kind=4 → MUSCL 3rd
    real(sp), intent(in), value   :: sensor
    real(8), intent(in)           :: a(4)     !< four-point stencil [i-1:i+2]
    real(8), intent(out)          :: al, ar
    call delta4(sensor, a, al, ar)
  end subroutine delta_r4

  !$dir inline
  attributes(device) subroutine delta_r6(id_acc, sensor, a, al, ar)
    integer(8), intent(in), value :: id_acc   !< dispatch key: kind=8 → MUSCL 4th
    real(sp), intent(in), value   :: sensor
    real(8), intent(in)           :: a(6)     !< six-point stencil [i-2:i+3]
    real(8), intent(out)          :: al, ar
    call delta6(sensor, a, al, ar)
  end subroutine delta_r6


  !> Interior-only Hybrid kernel for x-direction convective flux.
  !> Computes E at interfaces io+1..nx-(io+1) only (no boundary fallback branching).
  !> Shared memory: 2D shaped (x-sweep first), size -(io-1):threadsE%x+io+1.
  attributes(global) subroutine calc_hybrid_x_in(nx, ny, Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer, intent(in), value                :: nx                  !< grid points x
    integer, intent(in), value                :: ny                  !< grid points y
    !integer, intent(in), value                :: nz                  !< grid points z
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)       !< conservative variables
    real(8), intent(in), device, contiguous   :: T(nx,ny)         !< temperature (for KEEP path)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)    !< Ducros shock sensor
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2) !< x-direction flux
    integer i, j, it, jt, ii, i_base
    real(8), dimension(-(io-1):threadsE%x+io+1, threadsE%y), shared :: rho, u, v, p
    real(8), dimension(threadsE%x, threadsE%y), shared :: rhor, ur, vr, pr
    real(sp) fdx
    real(8) rhol, ul, vl,  pl
    it = threadIdx%x;  jt = threadIdx%y;  !kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    ! Phase 1: load shared memory tile (including halo of width io)
    do ii = it-io, threadsE%x+io+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(ii,jt) = Q(i,1,j);   u(ii,jt) = Q(i,2,j)
          v(ii,jt) = Q(i,3,j);   !w(ii,jt) = Q(i,4,j)
          p(ii,jt) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    i = i_base + it
    ! Phase 2: interior-only flux (no stencil-order boundary fallback)
    if (io+1 <= i .and. i <= nx-(io+1) .and. j <= ny-1 ) then
      fdx = 0.5_sp * (sensor(i,j) + sensor(i+1,j))
      if (fdx <= threshold) then
        ! KEEP: entropy/energy preserving flux — reads T from global memory
        block
          real(8) tmp(2*io+2)
          tmp = T(i-io:i+io+1, j)
          associate(uu => u)
            E(:,i,j-1) = KEEP(id_accuracy, &
                                   rho(it-io:it+io+1,jt), u(it-io:it+io+1,jt), &
                                     v(it-io:it+io+1,jt),  &
                                    uu(it-io:it+io+1,jt), p(it-io:it+io+1,jt), tmp, Normal_x)
          end associate
        end block
      else
        ! SLAU path: reconstruct left/right states via MUSCL interpolation
        call delta_r(id_accuracy, fdx, rho(it-io:it+io+1,jt), rhol, rhor(it,jt))
        call delta_r(id_accuracy, fdx,   u(it-io:it+io+1,jt),   ul,   ur(it,jt))
        call delta_r(id_accuracy, fdx,   v(it-io:it+io+1,jt),   vl,   vr(it,jt))
        !call delta_r(id_accuracy, fdx,   w(it-io:it+io+1,jt,kt),   wl,   wr(it,jt,kt))
        call delta_r(id_accuracy, fdx,   p(it-io:it+io+1,jt),   pl,   pr(it,jt))
      endif
    endif
    ! Sync before left-state overwrite: ensures all delta_r reads of shared memory
    ! are complete before any thread overwrites its shared memory slot with the left state.
    call syncthreads()
    ! Phase 3: overwrite shared memory with left state, then compute SLAU flux
    if (io+1 <= i .and. i <= nx-(io+1) .and. j <= ny-1 .and. fdx > threshold) then
      rho(it,jt) = rhol;   u(it,jt) = ul;   v(it,jt) = vl
           p(it,jt) = pl
      associate(un1 => u(it,jt), un2 => ur(it,jt))
        call SLAU(id_slau, rho(it,jt), rhor(it,jt), u(it,jt), ur(it,jt), v(it,jt), vr(it,jt), &
                   un1, un2, p(it,jt), pr(it,jt), Normal_x, 1.0_sp, &
                  E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))
      end associate
    endif
  end subroutine calc_hybrid_x_in


  !> Interior-only Hybrid kernel for y-direction convective flux.
  !> Shared memory: 2D shaped (y-sweep first: jj,it,kt).
  attributes(global) subroutine calc_hybrid_y_in(nx, ny,  Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer, intent(in), value                :: nx                  !< grid points x
    integer, intent(in), value                :: ny                  !< grid points y
   ! integer, intent(in), value                :: nz                  !< grid points z
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)       !< conservative variables
    real(8), intent(in), device, contiguous   :: T(nx,ny)         !< temperature (for KEEP path)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)    !< Ducros shock sensor
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1) !< y-direction flux
    integer i, j, it, jt, jj, j_base
    real(8), dimension(-(io-1):threadsF%y+io+1, threadsF%x), shared :: rho, u, v, p
    real(8), dimension(threadsF%y, threadsF%x), shared :: rhor, ur, vr, pr
    real(sp) fdy
    real(8) rhol, ul, vl, pl
    it = threadIdx%x;  jt = threadIdx%y; 
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    ! Phase 1: load shared memory tile
    do jj = jt-io, threadsF%y+io+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(jj,it) = Q(i,1,j);   u(jj,it) = Q(i,2,j)
          v(jj,it) = Q(i,3,j);  
          p(jj,it) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    j = j_base + jt
    ! Phase 2: interior-only flux
    if (io+1 <= j .and. j <= ny-(io+1) .and. i <= nx-1 ) then
      fdy = 0.5_sp * (sensor(i,j) + sensor(i,j+1))
      if (fdy <= threshold) then
        block
          real(8) tmp(2*io+2)
          tmp = T(i, j-io:j+io+1)
          associate(vv => v)
            F(:,i-1,j) = KEEP(id_accuracy, &
                                   rho(jt-io:jt+io+1,it), u(jt-io:jt+io+1,it), &
                                     v(jt-io:jt+io+1,it),  &
                                    vv(jt-io:jt+io+1,it), p(jt-io:jt+io+1,it), tmp, Normal_y)
          end associate
        end block
      else
        call delta_r(id_accuracy, fdy, rho(jt-io:jt+io+1,it), rhol, rhor(jt,it))
        call delta_r(id_accuracy, fdy,   u(jt-io:jt+io+1,it),   ul,   ur(jt,it))
        call delta_r(id_accuracy, fdy,   v(jt-io:jt+io+1,it),   vl,   vr(jt,it))
        !call delta_r(id_accuracy, fdy,   w(jt-io:jt+io+1,it,kt),   wl,   wr(jt,it,kt))
        call delta_r(id_accuracy, fdy,   p(jt-io:jt+io+1,it),   pl,   pr(jt,it))
      endif
    endif
    call syncthreads()
    if (io+1 <= j .and. j <= ny-(io+1) .and. i <= nx-1 .and. fdy > threshold) then
      rho(jt,it) = rhol;   u(jt,it) = ul;   v(jt,it) = vl
        p(jt,it) = pl
      associate(un1 => v(jt,it), un2 => vr(jt,it))
        call SLAU(id_slau, rho(jt,it), rhor(jt,it), u(jt,it), ur(jt,it), v(jt,it), vr(jt,it), &
                   un1, un2, p(jt,it), pr(jt,it), Normal_y, 1.0_sp, &
                  F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))
      end associate
    endif
  end subroutine calc_hybrid_y_in


  !> Interior-only Hybrid kernel for z-direction convective flux.
  !> Shared memory: 2D shaped (z-sweep first: kk,jt,it).

end module calc_hybrid_kernel_internal
