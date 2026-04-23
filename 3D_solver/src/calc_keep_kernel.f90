module calc_keep_kernel
  use mod_globals, only : threadsE, threadsF, threadsG
  use mod_constant, only : R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  implicit none
  private
  public calc_keep_x, calc_keep_y, calc_keep_z
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0
  
  interface calc_keep_x
    module procedure calc_keep_x2, calc_keep_x4, calc_keep_x6
  end interface calc_keep_x

  interface calc_keep_y
    module procedure calc_keep_y2, calc_keep_y4, calc_keep_y6
  end interface calc_keep_y

  interface calc_keep_z
    module procedure calc_keep_z2, calc_keep_z4, calc_keep_z6
  end interface calc_keep_z

  interface KEEP
    module procedure KEEP2, KEEP4, KEEP6
  end interface KEEP
contains
  include 'calc_keep_3d.f90'

  !> CUDA Fortran kernel for 6th-order KEEP scheme in x direction
  attributes(global) subroutine calc_keep_x6(id_accuracy, nx, ny, nz, Q, T, E)
    use mod_constant, only : Normal_x
    integer(8), intent(in), value :: id_accuracy         !< ID for accuracy, 8 means 6th-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2) !< Flux in x direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer ii, i_base, idx, offset_yz
    integer, parameter :: sx = threadsE%x + 5 !< tile size in x direction
    integer, parameter :: sy = threadsE%y     !< tile size in y direction
    integer, parameter :: sz = threadsE%z     !< tile size in z direction
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx + (kt-1) * sx * sy
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j,k);   u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k);   w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = it + offset_yz
    associate(uu => u)
    if (3 <= i .and. i <= nx-3) then
      E(i,:,j-1,k-1) = KEEP(id_accuracy, &
                            rho(idx-2:idx+3), u(idx-2:idx+3), &
                              v(idx-2:idx+3), w(idx-2:idx+3), &
                             uu(idx-2:idx+3), p(idx-2:idx+3), &
                            tmp(idx-2:idx+3), Normal_x)
    elseif (2 <= i .and. i <= nx-2) then
      E(i,:,j-1,k-1) = KEEP(id_accuracy4, & 
                            rho(idx-1:idx+2), u(idx-1:idx+2), &
                              v(idx-1:idx+2), w(idx-1:idx+2), &
                             uu(idx-1:idx+2), p(idx-1:idx+2), &
                            tmp(idx-1:idx+2), Normal_x)
    else
      E(i,:,j-1,k-1) = KEEP(id_accuracy2, &
                            rho(idx:idx+1), u(idx:idx+1), &
                              v(idx:idx+1), w(idx:idx+1), &
                             uu(idx:idx+1), p(idx:idx+1), &
                            tmp(idx:idx+1), Normal_x)
    endif
    end associate
  end subroutine calc_keep_x6


  !> CUDA Fortran kernel for 6th-order KEEP scheme in y direction
  attributes(global) subroutine calc_keep_y6(id_accuracy, nx, ny, nz, Q, T, F)
    use mod_constant, only : Normal_y
    integer(8), intent(in), value :: id_accuracy         !< ID for accuracy, 8 means 6th-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2) !< Flux in y direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer jj, j_base, idx, offset_xz
    integer, parameter :: sx = threadsF%x     !< tile size in x direction
    integer, parameter :: sy = threadsF%y + 5 !< tile size in y direction
    integer, parameter :: sz = threadsF%z     !< tile size in z direction
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1) * sy + (kt-1) * sy * sx
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j,k);   u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k);   w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = jt + offset_xz
    associate(vv => v)
    if (3 <= j .and. j <= ny-3) then
      F(i-1,:,j,k-1) = KEEP(id_accuracy, &
                            rho(idx-2:idx+3), u(idx-2:idx+3), &
                              v(idx-2:idx+3), w(idx-2:idx+3), &
                             vv(idx-2:idx+3), p(idx-2:idx+3), &
                            tmp(idx-2:idx+3), Normal_y)
    elseif (2 <= j .and. j <= ny-2) then
      F(i-1,:,j,k-1) = KEEP(id_accuracy4, &
                            rho(idx-1:idx+2), u(idx-1:idx+2), &
                              v(idx-1:idx+2), w(idx-1:idx+2), &
                             vv(idx-1:idx+2), p(idx-1:idx+2), &
                            tmp(idx-1:idx+2), Normal_y)
    else
      F(i-1,:,j,k-1) = KEEP(id_accuracy2, &
                            rho(idx:idx+1), u(idx:idx+1), &
                              v(idx:idx+1), w(idx:idx+1), &
                             vv(idx:idx+1), p(idx:idx+1), &
                            tmp(idx:idx+1), Normal_y)
    endif
    end associate
  end subroutine calc_keep_y6


  !> CUDA Fortran kernel for 6th-order KEEP scheme in z direction
  attributes(global) subroutine calc_keep_z6(id_accuracy, nx, ny, nz, Q, T, G)
    use mod_constant, only : Normal_z
    integer(8), intent(in), value :: id_accuracy         !< ID for accuracy, 8 means 6th-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1) !< Flux in z direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer kk, k_base, idx, offset_xy
    integer, parameter :: sx = threadsG%x     !< tile size in x direction
    integer, parameter :: sy = threadsG%y     !< tile size in y direction
    integer, parameter :: sz = threadsG%z + 5 !< tile size in z direction
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1) * sz + (it-1) * sz * sy
    do kk = kt-2, threadsG%z+3, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(i,1,j,k);   u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k);   w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = kt + offset_xy
    associate(ww => w)
    if (3 <= k .and. k <= nz-3) then
      G(i-1,:,j-1,k) = KEEP(id_accuracy, &
                            rho(idx-2:idx+3), u(idx-2:idx+3), &
                              v(idx-2:idx+3), w(idx-2:idx+3), &
                             ww(idx-2:idx+3), p(idx-2:idx+3), &
                            tmp(idx-2:idx+3), Normal_z)
    elseif (2 <= k .and. k <= nz-2) then
      G(i-1,:,j-1,k) = KEEP(id_accuracy4, &
                            rho(idx-1:idx+2), u(idx-1:idx+2), &
                              v(idx-1:idx+2), w(idx-1:idx+2), &
                             ww(idx-1:idx+2), p(idx-1:idx+2), &
                            tmp(idx-1:idx+2), Normal_z)
    else
      G(i-1,:,j-1,k) = KEEP(id_accuracy2, &
                            rho(idx:idx+1), u(idx:idx+1), &
                              v(idx:idx+1), w(idx:idx+1), &
                             ww(idx:idx+1), p(idx:idx+1), &
                            tmp(idx:idx+1), Normal_z)
    endif
    end associate
  end subroutine calc_keep_z6


  !> CUDA Fortran kernel for 4th-order KEEP scheme in x direction
  attributes(global) subroutine calc_keep_x4(id_accuracy, nx, ny, nz, Q, T, E)
    use mod_constant, only : Normal_x
    integer(4), intent(in), value :: id_accuracy         !< ID for accuracy, 4 means 4th-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2) !< Flux in x direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer ii, i_base, idx, offset_yz
    integer, parameter :: sx = threadsE%x + 3 !< tile size in x direction
    integer, parameter :: sy = threadsE%y     !< tile size in y direction
    integer, parameter :: sz = threadsE%z     !< tile size in z direction
    real(8), dimension(0:sx*sy*sz-1), shared :: rho, u, v, w, p, tmp
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx + (kt-1) * sx * sy
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j,k);   u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k);   w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = it + offset_yz
    associate(uu => u)
    if (2 <= i .and. i <= nx-2) then
      E(i,:,j-1,k-1) = KEEP(id_accuracy, &
                            rho(idx-1:idx+2), u(idx-1:idx+2), &
                              v(idx-1:idx+2), w(idx-1:idx+2), &
                             uu(idx-1:idx+2), p(idx-1:idx+2), &
                            tmp(idx-1:idx+2), Normal_x)
    else
      E(i,:,j-1,k-1) = KEEP(id_accuracy2, &
                            rho(idx:idx+1), u(idx:idx+1), &
                              v(idx:idx+1), w(idx:idx+1), &
                             uu(idx:idx+1), p(idx:idx+1), &
                            tmp(idx:idx+1), Normal_x)
    endif
    end associate
  end subroutine calc_keep_x4


  !> CUDA Fortran kernel for 4th-order KEEP scheme in y direction
  attributes(global) subroutine calc_keep_y4(id_accuracy, nx, ny, nz, Q, T, F)
    use mod_constant, only : Normal_y
    integer(4), intent(in), value :: id_accuracy         !< ID for accuracy, 4 means 4th-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2) !< Flux in y direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer jj, j_base, idx, offset_xz
    integer, parameter :: sx = threadsF%x     !< tile size in x direction
    integer, parameter :: sy = threadsF%y + 3 !< tile size in y direction
    integer, parameter :: sz = threadsF%z     !< tile size in z direction
    real(8), dimension(0:sx*sy*sz-1), shared :: rho, u, v, w, p, tmp
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1) * sy + (kt-1) * sy * sx
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j,k);   u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k);   w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = jt + offset_xz
    associate(vv => v)
    if (2 <= j .and. j <= ny-2) then
      F(i-1,:,j,k-1) = KEEP(id_accuracy, &
                            rho(idx-1:idx+2), u(idx-1:idx+2), &
                              v(idx-1:idx+2), w(idx-1:idx+2), &
                             vv(idx-1:idx+2), p(idx-1:idx+2), &
                            tmp(idx-1:idx+2), Normal_y)
    else
      F(i-1,:,j,k-1) = KEEP(id_accuracy2, &
                            rho(idx:idx+1), u(idx:idx+1), &
                              v(idx:idx+1), w(idx:idx+1), &
                             vv(idx:idx+1), p(idx:idx+1), &
                            tmp(idx:idx+1), Normal_y)
    endif
    end associate
  end subroutine calc_keep_y4


  !> CUDA Fortran kernel for 4th-order KEEP scheme in z direction
  attributes(global) subroutine calc_keep_z4(id_accuracy, nx, ny, nz, Q, T, G)
    use mod_constant, only : Normal_z
    integer(4), intent(in), value :: id_accuracy         !< ID for accuracy, 4 means 4th-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1) !< Flux in z direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer kk, k_base, idx, offset_xy
    integer, parameter :: sx = threadsG%x     !< tile size in x direction
    integer, parameter :: sy = threadsG%y     !< tile size in y direction
    integer, parameter :: sz = threadsG%z + 3 !< tile size in z direction
    real(8), dimension(0:sx*sy*sz-1), shared :: rho, u, v, w, p, tmp
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1) * sz + (it-1) * sz * sy
    do kk = kt-1, threadsG%z+2, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(i,1,j,k);   u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k);   w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k); tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = kt + offset_xy
    associate(ww => w)
    if (2 <= k .and. k <= nz-2) then
      G(i-1,:,j-1,k) = KEEP(id_accuracy, &
                            rho(idx-1:idx+2), u(idx-1:idx+2), &
                              v(idx-1:idx+2), w(idx-1:idx+2), &
                             ww(idx-1:idx+2), p(idx-1:idx+2), &
                            tmp(idx-1:idx+2), Normal_z)
    else
      G(i-1,:,j-1,k) = KEEP(id_accuracy2, &
                            rho(idx:idx+1), u(idx:idx+1), &
                              v(idx:idx+1), w(idx:idx+1), &
                             ww(idx:idx+1), p(idx:idx+1), &
                            tmp(idx:idx+1), Normal_z)
    endif
    end associate
  end subroutine calc_keep_z4


  !> CUDA Fortran kernel for 2nd-order KEEP scheme in x direction
  attributes(global) subroutine calc_keep_x2(id_accuracy, nx, ny, nz, Q, T, E)
    use mod_constant, only : Normal_x
    integer(2), intent(in), value :: id_accuracy         !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2) !< Flux in x direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    real(8), dimension(2) :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(i:i+1,1,j,k); u   = Q(i:i+1,2,j,k)
    v   = Q(i:i+1,3,j,k); w   = Q(i:i+1,4,j,k)
    p   = Q(i:i+1,5,j,k); tmp = T(i:i+1,j,k)
    E(i,:,j-1,k-1) = KEEP(id_accuracy, rho, u, v, w, u, p, tmp, Normal_x)
  end subroutine calc_keep_x2


  !> CUDA Fortran kernel for 2nd-order KEEP scheme in y direction
  attributes(global) subroutine calc_keep_y2(id_accuracy, nx, ny, nz, Q, T, F)
    use mod_constant, only : Normal_y
    integer(2), intent(in), value :: id_accuracy         !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2) !< Flux in y direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    real(8), dimension(2) :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(i,1,j:j+1,k); u   = Q(i,2,j:j+1,k)
    v   = Q(i,3,j:j+1,k); w   = Q(i,4,j:j+1,k)
    p   = Q(i,5,j:j+1,k); tmp = T(i,j:j+1,k)
    F(i-1,:,j,k-1) = KEEP(id_accuracy, rho, u, v, w, v, p, tmp, Normal_y)
  end subroutine calc_keep_y2


  !> CUDA Fortran kernel for 2nd-order KEEP scheme in z direction
  attributes(global) subroutine calc_keep_z2(id_accuracy, nx, ny, nz, Q, T, G)
    use mod_constant, only : Normal_z
    integer(2), intent(in), value :: id_accuracy         !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(8), intent(in), device, contiguous   :: T(nx,ny,nz)         !< Temperature
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1) !< Flux in z direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    real(8), dimension(2) :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(i,1,j,k:k+1); u   = Q(i,2,j,k:k+1)
    v   = Q(i,3,j,k:k+1); w   = Q(i,4,j,k:k+1)
    p   = Q(i,5,j,k:k+1); tmp = T(i,j,k:k+1)
    G(i-1,:,j-1,k) = KEEP(id_accuracy, rho, u, v, w, w, p, tmp, Normal_z)
  end subroutine calc_keep_z2
end module calc_keep_kernel

