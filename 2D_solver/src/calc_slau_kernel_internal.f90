module calc_slau_kernel_internal
  use mod_globals, only : id_accuracy, id_slau, gamma, threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1
  use calc_muscl
  use calc_hybrid
  implicit none
  private
  public calc_slau_x_in, calc_slau_y_in    !calc_slau_z_in
  !> io = 0, 1, 2 (2nd, 4th, 6th)
  integer, parameter :: io = kind(id_accuracy) / 3

  interface interp
    module procedure interp2, interp4, interp6
  end interface interp

  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU
contains
  include 'calc_slau_3d.f90'

  !$dir inline
  attributes(device) subroutine interp6(id_accuracy, rho,  u,  v, p, &
                                                     rhol, ul, vl, pl, &
                                                     rhor, ur, vr, pr, fd)
    integer(8), intent(in), value :: id_accuracy
    real(8), intent(in), contiguous :: rho(6), u(6), v(6), p(6)
    real(8), intent(out)            :: rhol, ul, vl, pl
    real(8), intent(out)            :: rhor, ur, vr, pr
    real(sp), intent(inout)         :: fd
    call delta6(fd, rho(:), rhol, rhor)
    call delta6(fd,   u(:),   ul,   ur)
    call delta6(fd,   v(:),   vl,   vr)
    !call delta6(fd,   w(:),   wl,   wr)
    call delta6(fd,   p(:),   pl,   pr)
    fd = wiggle_detector(p(2:5))
  end subroutine interp6


  !$dir inline
  attributes(device) subroutine interp4(id_accuracy, rho,  u,  v, p, &
                                                     rhol, ul, vl, pl, &
                                                     rhor, ur, vr, pr, fd)
    integer(4), intent(in), value :: id_accuracy
    real(8), intent(in), contiguous :: rho(4), u(4), v(4), p(4)
    real(8), intent(out)            :: rhol, ul, vl, pl
    real(8), intent(out)            :: rhor, ur, vr, pr
    real(sp), intent(inout)         :: fd
    call delta4(fd, rho(:), rhol, rhor)
    call delta4(fd,   u(:),   ul,   ur)
    call delta4(fd,   v(:),   vl,   vr)
    !call delta4(fd,   w(:),   wl,   wr)
    call delta4(fd,   p(:),   pl,   pr)
    fd = wiggle_detector(p(:))
  end subroutine interp4


  !$dir inline
  attributes(device) subroutine interp2(id_accuracy, rho,  u,  v,  p, &
                                                     rhol, ul, vl, pl, &
                                                     rhor, ur, vr, pr, fd)
    integer(2), intent(in), value :: id_accuracy
    real(8), intent(in), contiguous :: rho(2), u(2), v(2), p(2)
    real(8), intent(out)            :: rhol, ul, vl, pl
    real(8), intent(out)            :: rhor, ur, vr, pr
    real(sp), intent(inout)         :: fd
    rhol = rho(1); rhor = rho(2)
      ul =   u(1);   ur =   u(2)
      vl =   v(1);   vr =   v(2)
      !wl =   w(1);   wr =   w(2)
      pl =   p(1);   pr =   p(2)
      fd = 1.0_sp
  end subroutine interp2


  !> CUDA Fortran kernel for 4 points SLAU scheme in x direction
  attributes(global) subroutine calc_slau_x_in(nx, ny, Q, sensor, E)
    use mod_constant, only : Normal_x
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    !integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)    !< shock sensor
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2) !< Flux in x direction
    integer i,  j !< global index in physical space
    integer it, jt !< local index in a block
    integer ii, i_base, idx, idx_r, offset_yz, offset_yzr, i1, i2
    integer, parameter :: sx  = threadsE%x + 2*io + 1 !< tile size in x direction to calc high-order interpolation
    integer, parameter :: sxr = threadsE%x            !< tile size in x direction to store the results
    integer, parameter :: sy  = threadsE%y            !< tile size in y direction
    !integer, parameter :: sz  = threadsE%z            !< tile size in z direction
    real(8), dimension(-(io-1):sx*sy-io), shared :: rho,  u,  v,  p  !< 1st use: stensils, 2nd use: left
    real(8), dimension(sxr*sy), shared           :: rhor, ur, vr, pr !< right
    real(8) rhol, ul, vl, pl !< they are 40 bytes stack frame
    real(sp) fdx !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base     = (blockIdx%x-1)*blockDim%x
    offset_yz  = (jt-1) * sx  !+ (kt-1) * sx  * sy
    offset_yzr = (jt-1) * sxr !+ (kt-1) * sxr * sy
    do ii = it-io, threadsE%x+io+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j); u(idx) = Q(i,2,j); v(idx) = Q(i,3,j)
          !w(idx) = Q(i,4,j)
        p(idx) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (io+1 <= i .and. i <= nx-(io+1) .and. j <= ny-1) then
      idx   = it + offset_yz
      idx_r = it + offset_yzr
      i1    = idx - io
      i2    = idx + io + 1
      fdx   = 0.50_sp * (sensor(i,j) + sensor(i+1,j))
      call interp(id_accuracy, &
                  rho(i1:i2),  u(i1:i2),  v(i1:i2), p(i1:i2), &
                  rhol,        ul,        vl,       pl, &
                  rhor(idx_r), ur(idx_r), vr(idx_r), pr(idx_r), fdx)
    endif
    call syncthreads()
    if (io+1 <= i .and. i <= nx-(io+1) .and. j <= ny-1) then
      !> reuse smem to avoid using stack frame (rhol, ul, vl, wl, pl)
      rho(idx) = rhol; u(idx) = ul; v(idx) = vl; p(idx) = pl
      associate(un1 => u(idx), un2 => ur(idx_r))
        call SLAU(id_slau, rho(idx), rhor(idx_r), u(idx), ur(idx_r), v(idx), vr(idx_r), &
                  un1, un2, p(idx), pr(idx_r), Normal_x, fdx, &
                  E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))  !E(5,i,j-1)) unnecessary?
      end associate
    endif
  end subroutine calc_slau_x_in


  !> CUDA Fortran kernel for 4 points SLAU scheme in y direction
  attributes(global) subroutine calc_slau_y_in(nx, ny, Q, sensor, F)
    use mod_constant, only : Normal_y
    integer, intent(in), value                :: nx                  !< number of grid points in x direction
    integer, intent(in), value                :: ny                  !< number of grid points in y direction
    !integer, intent(in), value                :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)       !< Q(rho, u, v, w, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)    !< shock sensor
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1) !< Flux in y direction
    integer i,  j  !< global index in physical space
    integer it, jt !< local index in a block
    integer jj, j_base, idx, idx_r, offset_xz, offset_xzr, i1, i2
    integer, parameter :: sx  = threadsF%x            !< tile size in x direction
    integer, parameter :: sy  = threadsF%y + 2*io + 1 !< tile size in y direction to calc high-order interpolation
    integer, parameter :: syr = threadsF%y            !< tile size in y direction to store the results
    !integer, parameter :: sz  = threadsF%z            !< tile size in z direction
    real(8), dimension(-(io-1):sx*sy-io), shared :: rho,  u,  v,  p  !< 1st use: stencils, 2nd use: left
    real(8), dimension(sx*syr), shared           :: rhor, ur, vr, pr !< right
    real(8) rhol, ul, vl, pl !< they are 40 bytes stack frame
    real(sp) fdy !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base     = (blockIdx%y-1)*blockDim%y
    offset_xz  = (it-1) * sy  !+ (kt-1) * sy  * sx
    offset_xzr = (it-1) * syr !+ (kt-1) * syr * sx
    do jj = jt-io, threadsF%y+io+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j); u(idx) = Q(i,2,j); v(idx) = Q(i,3,j)
          !w(idx) = Q(i,4,j,k)
        p(idx) = Q(i,4,j) !!!!!!!!!!!!!!!!!!!!
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (io+1 <= j .and. j <= ny-(io+1) .and. i <= nx-1) then
      idx   = jt + offset_xz
      idx_r = jt + offset_xzr
      i1    = idx - io
      i2    = idx + io + 1
      fdy   = 0.50_sp * (sensor(i,j) + sensor(i,j+1))
      call interp(id_accuracy, &
                  rho(i1:i2),  u(i1:i2),  v(i1:i2), p(i1:i2), &
                  rhol,        ul,        vl,        pl, &
                  rhor(idx_r), ur(idx_r), vr(idx_r), pr(idx_r), fdy)
    endif
    call syncthreads()
    if (io+1 <= j .and. j <= ny-(io+1) .and. i <= nx-1) then
      !> reuse smem to avoid using stack frame (rhol, ul, vl, wl, pl)
      rho(idx) = rhol; u(idx) = ul; v(idx) = vl; p(idx) = pl
      associate(un1 => v(idx), un2 => vr(idx_r))
        call SLAU(id_slau, rho(idx), rhor(idx_r), u(idx), ur(idx_r), v(idx), vr(idx_r), &
                  un1, un2, p(idx), pr(idx_r), Normal_y, fdy, &
                  F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))   !F(5,i-1,j))  unnecessary?
      end associate
    endif
  end subroutine calc_slau_y_in

end module calc_slau_kernel_internal

