module calc_keep_kernel_internal
  use mod_globals, only : id_accuracy, threadsE, threadsF, threadsG
  use mod_constant, only : R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  implicit none
  private
  public calc_keep_x_in, calc_keep_y_in, calc_keep_z_in
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0
  !> io = 0, 1, 2 (2nd, 4th, 6th)
  integer, parameter :: io = kind(id_accuracy) / 3

  interface KEEP
    module procedure KEEP2, KEEP4, KEEP6
  end interface KEEP
contains
  include 'calc_keep_3d.f90'

  !> CUDA Fortran kernel for KEEP scheme in x direction
  attributes(global) subroutine calc_keep_x_in(nx, ny, nz, Q, T, E)
    use mod_constant, only : Normal_x
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(5,nx,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer ii, i_base, idx, offset_yz
    integer, parameter :: sx = threadsE%x + 2*io + 1 !< tile size in x direction
    integer, parameter :: sy = threadsE%y            !< tile size in y direction
    integer, parameter :: sz = threadsE%z            !< tile size in z direction
    real(8), dimension(-(io-1):sx*sy*sz-io), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx + (kt-1) * sx * sy
    do ii = it-io, threadsE%x+io+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(1,i,j,k);   u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k);   w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (io+1 <= i .and. i <= nx-(io+1) .and. j <= ny-1 .and. k <= nz-1) then
      associate(uu => u)
      idx = it + offset_yz
      E(:,i,j-1,k-1) = KEEP(id_accuracy, &
                            rho(idx-io:idx+io+1), u(idx-io:idx+io+1), &
                              v(idx-io:idx+io+1), w(idx-io:idx+io+1), &
                             uu(idx-io:idx+io+1), p(idx-io:idx+io+1), &
                            tmp(idx-io:idx+io+1), Normal_x)
      end associate
    endif
  end subroutine calc_keep_x_in


  !> CUDA Fortran kernel for KEEP scheme in y direction
  attributes(global) subroutine calc_keep_y_in(nx, ny, nz, Q, T, F)
    use mod_constant, only : Normal_y
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(5,nx,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer jj, j_base, idx, offset_xz
    integer, parameter :: sx = threadsF%x            !< tile size in x direction
    integer, parameter :: sy = threadsF%y + 2*io + 1 !< tile size in y direction
    integer, parameter :: sz = threadsF%z            !< tile size in z direction
    real(8), dimension(-(io-1):sx*sy*sz-io), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1) * sy + (kt-1) * sy * sx
    do jj = jt-io, threadsF%y+io+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(1,i,j,k);   u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k);   w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (io+1 <= j .and. j <= ny-(io+1) .and. i <= nx-1 .and. k <= nz-1) then
      associate(vv => v)
      idx = jt + offset_xz
      F(:,i-1,j,k-1) = KEEP(id_accuracy, &
                            rho(idx-io:idx+io+1), u(idx-io:idx+io+1), &
                              v(idx-io:idx+io+1), w(idx-io:idx+io+1), &
                             vv(idx-io:idx+io+1), p(idx-io:idx+io+1), &
                            tmp(idx-io:idx+io+1), Normal_y)
      end associate
    endif
  end subroutine calc_keep_y_in


  !> CUDA Fortran kernel for KEEP scheme in z direction
  attributes(global) subroutine calc_keep_z_in(nx, ny, nz, Q, T, G)
    use mod_constant, only : Normal_z
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(5,nx,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer kk, k_base, idx, offset_xy
    integer, parameter :: sx = threadsG%x            !< tile size in x direction
    integer, parameter :: sy = threadsG%y            !< tile size in y direction
    integer, parameter :: sz = threadsG%z + 2*io + 1 !< tile size in z direction
    real(8), dimension(-(io-1):sx*sy*sz-io), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1) * sz + (it-1) * sz * sy
    do kk = kt-io, threadsG%z+io+1, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(1,i,j,k);   u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k);   w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (io+1 <= k .and. k <= nz-(io+1) .and. i <= nx-1 .and. j <= ny-1) then
      associate(ww => w)
      idx = kt + offset_xy
      G(:,i-1,j-1,k) = KEEP(id_accuracy, &
                            rho(idx-io:idx+io+1), u(idx-io:idx+io+1), &
                              v(idx-io:idx+io+1), w(idx-io:idx+io+1), &
                             ww(idx-io:idx+io+1), p(idx-io:idx+io+1), &
                            tmp(idx-io:idx+io+1), Normal_z)
      end associate
    endif
  end subroutine calc_keep_z_in
end module calc_keep_kernel_internal

