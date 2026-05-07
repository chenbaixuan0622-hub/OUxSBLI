module calc_slau_kernel
  use mod_globals, only : id_slau, gamma, threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1
  use calc_muscl
  use calc_hybrid
  implicit none
  private
  public calc_slau_x, calc_slau_y, calc_slau_z
  
  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU

  interface calc_slau_x
    module procedure calc_slau_x2, calc_slau_x4, calc_slau_x6
  end interface calc_slau_x

  interface calc_slau_y
    module procedure calc_slau_y2, calc_slau_y4, calc_slau_y6
  end interface calc_slau_y

  interface calc_slau_z
    module procedure calc_slau_z2, calc_slau_z4, calc_slau_z6
  end interface calc_slau_z
contains
  include 'calc_slau_3d.f90'

  !> CUDA Fortran kernel for 6 points SLAU scheme in x direction
  attributes(global) subroutine calc_slau_x6(id_accuracy, nx, ny, nz, Q, sensor, E)
    use mod_constant, only : Normal_x
    integer(8), intent(in), value :: id_accuracy         !< ID for accuracy, 8 means 6 ppoints
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer ii, i_base, idx, idx_r, offset_yz, offset_yzr
    integer, parameter :: sx  = threadsE%x + 5 !< tile size in x direction to calc high-order interpolation
    integer, parameter :: sxr = threadsE%x     !< tile size in x direction to store the results
    integer, parameter :: sy  = threadsE%y     !< tile size in y direction
    integer, parameter :: sz  = threadsE%z     !< tile size in z direction
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho,  u,  v,  w,  p
    real(8), dimension(sxr*sy*sz), shared     :: rhor, ur, vr, wr, pr
    real(8) rhol, ul, vl, wl, pl
    real(sp) fdx !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base     = (blockIdx%x-1)*blockDim%x
    offset_yz  = (jt-1) * sx  + (kt-1) * sx  * sy
    offset_yzr = (jt-1) * sxr + (kt-1) * sxr * sy
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx   = it + offset_yz
    idx_r = it + offset_yzr
    fdx   = 0.5_sp * (sensor(i,j,k) + sensor(i+1,j,k))
    if (3 <= i .and. i <= nx-3) then
      call delta6(fdx, rho(idx-2:idx+3), rhol, rhor(idx_r))
      call delta6(fdx,   u(idx-2:idx+3),   ul,   ur(idx_r))
      call delta6(fdx,   v(idx-2:idx+3),   vl,   vr(idx_r))
      call delta6(fdx,   w(idx-2:idx+3),   wl,   wr(idx_r))
      call delta6(fdx,   p(idx-2:idx+3),   pl,   pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    elseif (2 <= i .and. i <= nx-2) then
      call delta4(fdx, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdx,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdx,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdx,   w(idx-1:idx+2),   wl,   wr(idx_r))
      call delta4(fdx,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        wl =   w(idx);   wr(idx_r) =   w(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdx = 1.0_sp
    endif
    associate(un1 => ul, un2 => ur(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                wl, wr(idx_r), un1, un2, pl, pr(idx_r), Normal_x, fdx, &
                E(1,i,j-1,k-1), E(2,i,j-1,k-1), E(3,i,j-1,k-1), E(4,i,j-1,k-1), E(5,i,j-1,k-1))
    end associate
  end subroutine calc_slau_x6


  !> CUDA Fortran kernel for 6 points SLAU scheme in y direction
  attributes(global) subroutine calc_slau_y6(id_accuracy, nx, ny, nz, Q, sensor, F)
    use mod_constant, only : Normal_y
    integer(8), intent(in), value :: id_accuracy         !< ID for accuracy, 8 means 6 ppoints
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer jj, j_base, idx, idx_r, offset_xz, offset_xzr
    integer, parameter :: sx  = threadsF%x     !< tile size in x direction
    integer, parameter :: sy  = threadsF%y + 5 !< tile size in y direction to calc high-order interpolation
    integer, parameter :: syr = threadsF%y     !< tile size in y direction to store the results
    integer, parameter :: sz  = threadsF%z     !< tile size in z direction
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho,  u,  v,  w,  p
    real(8), dimension(sx*syr*sz), shared     :: rhor, ur, vr, wr, pr
    real(8) rhol, ul, vl, wl, pl
    real(sp) fdy !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base     = (blockIdx%y-1)*blockDim%y
    offset_xz  = (it-1) * sy  + (kt-1) * sy  * sx
    offset_xzr = (it-1) * syr + (kt-1) * syr * sx
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx   = jt + offset_xz
    idx_r = jt + offset_xzr
    fdy   = 0.5_sp * (sensor(i,j,k) + sensor(i,j+1,k))
    if (3 <= j .and. j <= ny-3) then
      call delta6(fdy, rho(idx-2:idx+3), rhol, rhor(idx_r))
      call delta6(fdy,   u(idx-2:idx+3),   ul,   ur(idx_r))
      call delta6(fdy,   v(idx-2:idx+3),   vl,   vr(idx_r))
      call delta6(fdy,   w(idx-2:idx+3),   wl,   wr(idx_r))
      call delta6(fdy,   p(idx-2:idx+3),   pl,   pr(idx_r))
      fdy = wiggle_detector(p(idx-1:idx+2))
    elseif (2 <= j .and. j <= ny-2) then
      call delta4(fdy, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdy,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdy,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdy,   w(idx-1:idx+2),   wl,   wr(idx_r))
      call delta4(fdy,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdy = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        wl =   w(idx);   wr(idx_r) =   w(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdy = 1.0_sp
    endif
    associate(un1 => vl, un2 => vr(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                wl, wr(idx_r), un1, un2, pl, pr(idx_r), Normal_y, fdy, &
                F(1,i-1,j,k-1), F(2,i-1,j,k-1), F(3,i-1,j,k-1), F(4,i-1,j,k-1), F(5,i-1,j,k-1))
    end associate
  end subroutine calc_slau_y6


  !> CUDA Fortran kernel for 6 points SLAU scheme in z direction
  attributes(global) subroutine calc_slau_z6(id_accuracy, nx, ny, nz, Q, sensor, G)
    use mod_constant, only : Normal_z
    integer(8), intent(in), value :: id_accuracy         !< ID for accuracy, 8 means 6 ppoints
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer kk, k_base, idx, idx_r, offset_xy, offset_xyr
    integer, parameter :: sx  = threadsG%x     !< tile size in x direction
    integer, parameter :: sy  = threadsG%y     !< tile size in y direction
    integer, parameter :: sz  = threadsG%z + 5 !< tile size in z direction to calc high-order interpolation
    integer, parameter :: szr = threadsG%z     !< tile size in z direction to store the results
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho,  u,  v,  w,  p
    real(8), dimension(sx*sy*szr), shared     :: rhor, ur, vr, wr, pr
    real(8) rhol, ul, vl, wl, pl
    real(sp) fdz !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base     = (blockIdx%z-1)*blockDim%z
    offset_xy  = (jt-1) * sz  + (it-1) * sz  * sy
    offset_xyr = (jt-1) * szr + (it-1) * szr * sy
    do kk = kt-2, threadsG%z+3, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx   = kt + offset_xy
    idx_r = kt + offset_xyr
    fdz   = 0.5_sp * (sensor(i,j,k) + sensor(i,j,k+1))
    if (3 <= k .and. k <= nz-3) then
      call delta6(fdz, rho(idx-2:idx+3), rhol, rhor(idx_r))
      call delta6(fdz,   u(idx-2:idx+3),   ul,   ur(idx_r))
      call delta6(fdz,   v(idx-2:idx+3),   vl,   vr(idx_r))
      call delta6(fdz,   w(idx-2:idx+3),   wl,   wr(idx_r))
      call delta6(fdz,   p(idx-2:idx+3),   pl,   pr(idx_r))
      fdz = wiggle_detector(p(idx-1:idx+2))
    elseif (2 <= k .and. k <= nz-2) then
      call delta4(fdz, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdz,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdz,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdz,   w(idx-1:idx+2),   wl,   wr(idx_r))
      call delta4(fdz,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdz = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        wl =   w(idx);   wr(idx_r) =   w(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdz = 1.0_sp
    endif
    associate(un1 => wl, un2 => wr(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                wl, wr(idx_r), un1, un2, pl, pr(idx_r), Normal_z, fdz, &
                G(1,i-1,j-1,k), G(2,i-1,j-1,k), G(3,i-1,j-1,k), G(4,i-1,j-1,k), G(5,i-1,j-1,k))
    end associate
  end subroutine calc_slau_z6


  !> CUDA Fortran kernel for 4 points SLAU scheme in x direction
  attributes(global) subroutine calc_slau_x4(id_accuracy, nx, ny, nz, Q, sensor, E)
    use mod_constant, only : Normal_x
    integer(4), intent(in), value :: id_accuracy         !< ID for accuracy, 4 means 4 ppoints
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer ii, i_base, idx, idx_r, offset_yz, offset_yzr
    integer, parameter :: sx  = threadsE%x + 3 !< tile size in x direction to calc high-order interpolation
    integer, parameter :: sxr = threadsE%x     !< tile size in x direction to store the results
    integer, parameter :: sy  = threadsE%y     !< tile size in y direction
    integer, parameter :: sz  = threadsE%z     !< tile size in z direction
    real(8), dimension(0:sx*sy*sz-1), shared :: rho,  u,  v,  w,  p  !< smem to calc high-order interpolation
    real(8), dimension(sxr*sy*sz), shared    :: rhor, ur, vr, wr, pr !< smem to store the results
    real(8) rhol, ul, vl, wl, pl
    real(sp) fdx !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base     = (blockIdx%x-1)*blockDim%x
    offset_yz  = (jt-1) * sx  + (kt-1) * sx  * sy
    offset_yzr = (jt-1) * sxr + (kt-1) * sxr * sy
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx   = it + offset_yz
    idx_r = it + offset_yzr
    fdx   = 0.5_sp * (sensor(i,j,k) + sensor(i+1,j,k))
    if (2 <= i .and. i <= nx-2) then
      call delta4(fdx, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdx,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdx,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdx,   w(idx-1:idx+2),   wl,   wr(idx_r))
      call delta4(fdx,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        wl =   w(idx);   wr(idx_r) =   w(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdx = 1.0_sp
    endif
    associate(un1 => ul, un2 => ur(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                wl, wr(idx_r), un1, un2, pl, pr(idx_r), Normal_x, fdx, &
                E(1,i,j-1,k-1), E(2,i,j-1,k-1), E(3,i,j-1,k-1), E(4,i,j-1,k-1), E(5,i,j-1,k-1))
    end associate
  end subroutine calc_slau_x4


  !> CUDA Fortran kernel for 4 points SLAU scheme in y direction
  attributes(global) subroutine calc_slau_y4(id_accuracy, nx, ny, nz, Q, sensor, F)
    use mod_constant, only : Normal_y
    integer(4), intent(in), value :: id_accuracy         !< ID for accuracy, 4 means 4 ppoints
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer jj, j_base, idx, idx_r, offset_xz, offset_xzr
    integer, parameter :: sx  = threadsF%x     !< tile size in x direction
    integer, parameter :: sy  = threadsF%y + 3 !< tile size in y direction to calc high-order interpolation
    integer, parameter :: syr = threadsF%y     !< tile size in y direction to store the results
    integer, parameter :: sz  = threadsF%z     !< tile size in z direction
    real(8), dimension(0:sx*sy*sz-1), shared :: rho,  u,  v,  w,  p  !< smem to calc high-order interpolation
    real(8), dimension(sx*syr*sz), shared    :: rhor, ur, vr, wr, pr !< smem to store the results
    real(8) rhol, ul, vl, wl, pl
    real(sp) fdy !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base     = (blockIdx%y-1)*blockDim%y
    offset_xz  = (it-1) * sy  + (kt-1) * sy  * sx
    offset_xzr = (it-1) * syr + (kt-1) * syr * sx
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx   = jt + offset_xz
    idx_r = jt + offset_xzr
    fdy   = 0.5_sp * (sensor(i,j,k) + sensor(i,j+1,k))
    if (2 <= j .and. j <= ny-2) then
      call delta4(fdy, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdy,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdy,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdy,   w(idx-1:idx+2),   wl,   wr(idx_r))
      call delta4(fdy,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdy = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        wl =   w(idx);   wr(idx_r) =   w(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdy = 1.0_sp
    endif
    associate(un1 => vl, un2 => vr(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                wl, wr(idx_r), un1, un2, pl, pr(idx_r), Normal_y, fdy, &
                F(1,i-1,j,k-1), F(2,i-1,j,k-1), F(3,i-1,j,k-1), F(4,i-1,j,k-1), F(5,i-1,j,k-1))
    end associate
  end subroutine calc_slau_y4


  !> CUDA Fortran kernel for 4 points SLAU scheme in z direction
  attributes(global) subroutine calc_slau_z4(id_accuracy, nx, ny, nz, Q, sensor, G)
    use mod_constant, only : Normal_z
    integer(4), intent(in), value :: id_accuracy         !< ID for accuracy, 4 means 4 ppoints
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer kk, k_base, idx, idx_r, offset_xy, offset_xyr
    integer, parameter :: sx  = threadsG%x     !< tile size in x direction
    integer, parameter :: sy  = threadsG%y     !< tile size in y direction
    integer, parameter :: sz  = threadsG%z + 3 !< tile size in z direction to calc high-order interpolation
    integer, parameter :: szr = threadsG%z     !< tile size in z direction to store the results
    real(8), dimension(0:sx*sy*sz-1), shared :: rho,  u,  v,  w,  p  !< smem to calc high-order interpolation
    real(8), dimension(sx*sy*szr), shared    :: rhor, ur, vr, wr, pr !< smem to store the results
    real(8) rhol, ul, vl, wl, pl
    real(sp) fdz !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base     = (blockIdx%z-1)*blockDim%z
    offset_xy  = (jt-1) * sz  + (it-1) * sz  * sy
    offset_xyr = (jt-1) * szr + (it-1) * szr * sy
    do kk = kt-1, threadsG%z+2, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx   = kt + offset_xy
    idx_r = kt + offset_xyr
    fdz   = 0.5_sp * (sensor(i,j,k) + sensor(i,j,k+1))
    if (2 <= k .and. k <= nz-2) then
      call delta4(fdz, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdz,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdz,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdz,   w(idx-1:idx+2),   wl,   wr(idx_r))
      call delta4(fdz,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdz = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        wl =   w(idx);   wr(idx_r) =   w(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdz = 1.0_sp
    endif
    associate(un1 => wl, un2 => wr(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                wl, wr(idx_r), un1, un2, pl, pr(idx_r), Normal_z, fdz, &
                G(1,i-1,j-1,k), G(2,i-1,j-1,k), G(3,i-1,j-1,k), G(4,i-1,j-1,k), G(5,i-1,j-1,k))
    end associate
  end subroutine calc_slau_z4


  !> CUDA Fortran kernel for 2 points SLAU scheme in x direction
  attributes(global) subroutine calc_slau_x2(id_accuracy, nx, ny, nz, Q, sensor, E)
    use mod_constant, only : Normal_x
    integer(2), intent(in), value :: id_accuracy         !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer ii, idx, i_base, offset_yz
    integer, parameter :: sx = threadsE%x+1 !< tile size in x direction
    integer, parameter :: sy = threadsE%y   !< tile size in y direction
    integer, parameter :: sz = threadsE%z   !< tile size in z direction
    real(8), dimension(sx*sy*sz), shared :: rho, u, v, w, p
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx + (kt-1) * sx * sy
    do ii = it, threadsE%x+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = it + offset_yz
    associate(un1 => u(idx), un2 => u(idx+1))
      call SLAU(id_slau, rho(idx), rho(idx+1), u(idx), u(idx+1), v(idx), v(idx+1), &
                w(idx), w(idx+1), un1, un2, p(idx), p(idx+1), Normal_x, 1.0_sp, &
                E(1,i,j-1,k-1), E(2,i,j-1,k-1), E(3,i,j-1,k-1), E(4,i,j-1,k-1), E(5,i,j-1,k-1))
    end associate
  end subroutine calc_slau_x2


  !> CUDA Fortran kernel for 2 points SLAU scheme in y direction
  attributes(global) subroutine calc_slau_y2(id_accuracy, nx, ny, nz, Q, sensor, F)
    use mod_constant, only : Normal_y
    integer(2), intent(in), value :: id_accuracy         !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer jj, idx, j_base, offset_xz
    integer, parameter :: sx = threadsF%x   !< tile size in x direction
    integer, parameter :: sy = threadsF%y+1 !< tile size in y direction
    integer, parameter :: sz = threadsF%z   !< tile size in z direction
    real(8), dimension(sx*sy*sz), shared :: rho, u, v, w, p
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1) * sy + (kt-1) * sy * sx
    do jj = jt, threadsF%y+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = jt + offset_xz
    associate(un1 => v(idx), un2 => v(idx+1))
      call SLAU(id_slau, rho(idx), rho(idx+1), u(idx), u(idx+1), v(idx), v(idx+1), &
                w(idx), w(idx+1), un1, un2, p(idx), p(idx+1), Normal_y, 1.0_sp, &
                F(1,i-1,j,k-1), F(2,i-1,j,k-1), F(3,i-1,j,k-1), F(4,i-1,j,k-1), F(5,i-1,j,k-1))
    end associate
  end subroutine calc_slau_y2


  !> CUDA Fortran kernel for 2 points SLAU scheme in z direction
  attributes(global) subroutine calc_slau_z2(id_accuracy, nx, ny, nz, Q, sensor, G)
    use mod_constant, only : Normal_z
    integer(2), intent(in), value :: id_accuracy         !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)    !< shock sensor
    real(8), intent(out), device, contiguous  :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    integer i,  j,  k  !< global index in physical space
    integer it, jt, kt !< local index in a block
    integer kk, idx, k_base, offset_xy
    integer, parameter :: sx = threadsG%x   !< tile size in x direction
    integer, parameter :: sy = threadsG%y   !< tile size in y direction
    integer, parameter :: sz = threadsG%z+1 !< tile size in z direction
    real(8), dimension(sx*sy*sz), shared :: rho, u, v, w, p
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1) * sz + (it-1) * sz * sy
    do kk = kt, threadsG%z+1, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(i,1,j,k)
          u(idx) = Q(i,2,j,k)
          v(idx) = Q(i,3,j,k)
          w(idx) = Q(i,4,j,k)
          p(idx) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = kt + offset_xy
    associate(un1 => w(idx), un2 => w(idx+1))
      call SLAU(id_slau, rho(idx), rho(idx+1), u(idx), u(idx+1), v(idx), v(idx+1), &
                w(idx), w(idx+1), un1, un2, p(idx), p(idx+1), Normal_z, 1.0_sp, &
                G(1,i-1,j-1,k), G(2,i-1,j-1,k), G(3,i-1,j-1,k), G(4,i-1,j-1,k), G(5,i-1,j-1,k))
    end associate
  end subroutine calc_slau_z2
end module calc_slau_kernel

