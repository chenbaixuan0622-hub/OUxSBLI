module calc_slau_kernel
  use mod_globals, only : gamma, threadsE, threadsF
  use mod_constant, only : over_gamma_1, id_slau
  use calc_muscl
  use calc_hybrid
  implicit none
  private
  public calc_slau_x, calc_slau_y
  
  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU

  interface calc_slau_x
    module procedure calc_slau_x2, calc_slau_x4, calc_slau_x6
  end interface calc_slau_x

  interface calc_slau_y
    module procedure calc_slau_y2, calc_slau_y4, calc_slau_y6
  end interface calc_slau_y

contains
  include 'calc_slau_2d.f90'

  !> CUDA Fortran kernel for 6 points SLAU scheme in x direction
  attributes(global) subroutine calc_slau_x6(id_accuracy, nx, ny, Q, sensor, E)
    use mod_constant, only : Normal_x
    integer(8), intent(in), value             :: id_accuracy    !< ID for accuracy, 8 means 6 ppoints
    integer, intent(in), value                :: nx             !< number of grid points in x direction
    integer, intent(in), value                :: ny             !< number of grid points in y direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)     !< Q(rho, u, v, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)  !< shock sensor
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2) !< Flux in x direction
    integer i,  j  !< global index in physical space
    integer it, jt !< local index in a block
    integer ii, i_base, idx, idx_r, offset_yz, offset_yzr
    integer, parameter :: sx  = threadsE%x + 5 !< tile size in x direction to calc high-order interpolation
    integer, parameter :: sxr = threadsE%x     !< tile size in x direction to store the results
    integer, parameter :: sy  = threadsE%y     !< tile size in y direction
    real(8), dimension(-1:sx*sy-2), shared :: rho,  u,  v,  p
    real(8), dimension(sxr*sy), shared     :: rhor, ur, vr, pr
    real(8) rhol, ul, vl, pl
    real(sp) fdx !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    i_base     = (blockIdx%x-1)*blockDim%x
    offset_yz  = (jt-1) * sx
    offset_yzr = (jt-1) * sxr
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j)
          u(idx) = Q(i,2,j)
          v(idx) = Q(i,3,j)
          p(idx) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j) return
    idx   = it + offset_yz
    idx_r = it + offset_yzr
    fdx   = 0.5_sp * (sensor(i,j) + sensor(i+1,j))
    if (3 <= i .and. i <= nx-3) then
      call delta6(fdx, rho(idx-2:idx+3), rhol, rhor(idx_r))
      call delta6(fdx,   u(idx-2:idx+3),   ul,   ur(idx_r))
      call delta6(fdx,   v(idx-2:idx+3),   vl,   vr(idx_r))
      call delta6(fdx,   p(idx-2:idx+3),   pl,   pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    elseif (2 <= i .and. i <= nx-2) then
      call delta4(fdx, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdx,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdx,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdx,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdx = 1.0_sp
    endif
    associate(un1 => ul, un2 => ur(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                un1, un2, pl, pr(idx_r), Normal_x, fdx, &
                E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))
    end associate 
  end subroutine calc_slau_x6


  !> CUDA Fortran kernel for 6 points SLAU scheme in y direction
  attributes(global) subroutine calc_slau_y6(id_accuracy, nx, ny, Q, sensor, F)
    use mod_constant, only : Normal_y
    integer(8), intent(in), value             :: id_accuracy    !< ID for accuracy, 8 means 6 ppoints
    integer, intent(in), value                :: nx             !< number of grid points in x direction
    integer, intent(in), value                :: ny             !< number of grid points in y direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)     !< Q(rho, u, v, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)  !< shock sensor
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1) !< Flux in y direction
    integer i,  j  !< global index in physical space
    integer it, jt !< local index in a block
    integer jj, j_base, idx, idx_r, offset_xz, offset_xzr
    integer, parameter :: sx  = threadsF%x     !< tile size in x direction
    integer, parameter :: sy  = threadsF%y + 5 !< tile size in y direction to calc high-order interpolation
    integer, parameter :: syr = threadsF%y     !< tile size in y direction to store the results
    real(8), dimension(-1:sx*sy-2), shared :: rho,  u,  v,  p
    real(8), dimension(sx*syr), shared     :: rhor, ur, vr, pr
    real(8) rhol, ul, vl, pl
    real(sp) fdy !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j_base     = (blockIdx%y-1)*blockDim%y
    offset_xz  = (it-1) * sy
    offset_xzr = (it-1) * syr
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j)
          u(idx) = Q(i,2,j)
          v(idx) = Q(i,3,j)
          p(idx) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    idx   = jt + offset_xz
    idx_r = jt + offset_xzr
    fdy   = 0.5_sp * (sensor(i,j) + sensor(i,j+1))
    if (3 <= j .and. j <= ny-3) then
      call delta6(fdy, rho(idx-2:idx+3), rhol, rhor(idx_r))
      call delta6(fdy,   u(idx-2:idx+3),   ul,   ur(idx_r))
      call delta6(fdy,   v(idx-2:idx+3),   vl,   vr(idx_r))
      call delta6(fdy,   p(idx-2:idx+3),   pl,   pr(idx_r))
      fdy = wiggle_detector(p(idx-1:idx+2))
    elseif (2 <= j .and. j <= ny-2) then
      call delta4(fdy, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdy,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdy,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdy,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdy = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdy = 1.0_sp
    endif
    associate(un1 => vl, un2 => vr(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                un1, un2, pl, pr(idx_r), Normal_y, fdy, &
                F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))
    end associate
  end subroutine calc_slau_y6


  !> CUDA Fortran kernel for 4 points SLAU scheme in x direction
  attributes(global) subroutine calc_slau_x4(id_accuracy, nx, ny, Q, sensor, E)
    use mod_constant, only : Normal_x
    integer(4), intent(in), value             :: id_accuracy    !< ID for accuracy, 4 means 4 ppoints
    integer, intent(in), value                :: nx             !< number of grid points in x direction
    integer, intent(in), value                :: ny             !< number of grid points in y direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)     !< Q(rho, u, v, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)  !< shock sensor
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2) !< Flux in x direction
    integer i,  j  !< global index in physical space
    integer it, jt !< local index in a block
    integer ii, i_base, idx, idx_r, offset_yz, offset_yzr
    integer, parameter :: sx  = threadsE%x + 3 !< tile size in x direction to calc high-order interpolation
    integer, parameter :: sxr = threadsE%x     !< tile size in x direction to store the results
    integer, parameter :: sy  = threadsE%y     !< tile size in y direction
    real(8), dimension(0:sx*sy-1), shared :: rho,  u,  v,  p  !< smem to calc high-order interpolation
    real(8), dimension(sxr*sy), shared    :: rhor, ur, vr, pr !< smem to store the results
    real(8) rhol, ul, vl, pl
    real(sp) fdx !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    i_base     = (blockIdx%x-1)*blockDim%x
    offset_yz  = (jt-1) * sx
    offset_yzr = (jt-1) * sxr
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j)
          u(idx) = Q(i,2,j)
          v(idx) = Q(i,3,j)
          p(idx) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j) return
    idx   = it + offset_yz
    idx_r = it + offset_yzr
    fdx   = 0.5_sp * (sensor(i,j) + sensor(i+1,j))
    if (2 <= i .and. i <= nx-2) then
      call delta4(fdx, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdx,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdx,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdx,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdx = 1.0_sp
    endif
    associate(un1 => ul, un2 => ur(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                un1, un2, pl, pr(idx_r), Normal_x, fdx, &
                E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))
    end associate
  end subroutine calc_slau_x4


  !> CUDA Fortran kernel for 4 points SLAU scheme in y direction
  attributes(global) subroutine calc_slau_y4(id_accuracy, nx, ny, Q, sensor, F)
    use mod_constant, only : Normal_y
    integer(4), intent(in), value             :: id_accuracy    !< ID for accuracy, 4 means 4 ppoints
    integer, intent(in), value                :: nx             !< number of grid points in x direction
    integer, intent(in), value                :: ny             !< number of grid points in y direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)     !< Q(rho, u, v, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)  !< shock sensor
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1) !< Flux in y direction
    integer i,  j  !< global index in physical space
    integer it, jt !< local index in a block
    integer jj, j_base, idx, idx_r, offset_xz, offset_xzr
    integer, parameter :: sx  = threadsF%x     !< tile size in x direction
    integer, parameter :: sy  = threadsF%y + 3 !< tile size in y direction to calc high-order interpolation
    integer, parameter :: syr = threadsF%y     !< tile size in y direction to store the results
    real(8), dimension(0:sx*sy-1), shared :: rho,  u,  v,  p  !< smem to calc high-order interpolation
    real(8), dimension(sx*syr), shared    :: rhor, ur, vr, pr !< smem to store the results
    real(8) rhol, ul, vl, pl
    real(sp) fdy !< 1st use: shock sensor, 2nd use: wiggle ditector
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j_base     = (blockIdx%y-1)*blockDim%y
    offset_xz  = (it-1) * sy
    offset_xzr = (it-1) * syr
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j)
          u(idx) = Q(i,2,j)
          v(idx) = Q(i,3,j)
          p(idx) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    idx   = jt + offset_xz
    idx_r = jt + offset_xzr
    fdy   = 0.5_sp * (sensor(i,j) + sensor(i,j+1))
    if (2 <= j .and. j <= ny-2) then
      call delta4(fdy, rho(idx-1:idx+2), rhol, rhor(idx_r))
      call delta4(fdy,   u(idx-1:idx+2),   ul,   ur(idx_r))
      call delta4(fdy,   v(idx-1:idx+2),   vl,   vr(idx_r))
      call delta4(fdy,   p(idx-1:idx+2),   pl,   pr(idx_r))
      fdy = wiggle_detector(p(idx-1:idx+2))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
        ul =   u(idx);   ur(idx_r) =   u(idx+1)
        vl =   v(idx);   vr(idx_r) =   v(idx+1)
        pl =   p(idx);   pr(idx_r) =   p(idx+1)
       fdy = 1.0_sp
    endif
    associate(un1 => vl, un2 => vr(idx_r))
      call SLAU(id_slau, rhol, rhor(idx_r), ul, ur(idx_r), vl, vr(idx_r), &
                un1, un2, pl, pr(idx_r), Normal_y, fdy, &
                F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))
    end associate
  end subroutine calc_slau_y4


  !> CUDA Fortran kernel for 2 points SLAU scheme in x direction
  attributes(global) subroutine calc_slau_x2(id_accuracy, nx, ny, Q, sensor, E)
    use mod_constant, only : Normal_x
    integer(2), intent(in), value             :: id_accuracy    !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx             !< number of grid points in x direction
    integer, intent(in), value                :: ny             !< number of grid points in y direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)     !< Q(rho, u, v, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)  !< shock sensor
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2) !< Flux in x direction
    integer i,  j !< global index in physical space
    integer it, jt !< local index in a block
    integer ii, idx, i_base, offset_yz
    integer, parameter :: sx = threadsE%x+1 !< tile size in x direction
    integer, parameter :: sy = threadsE%y   !< tile size in y direction
    real(8), dimension(sx*sy), shared :: rho, u, v, p
    it = threadIdx%x
    jt = threadIdx%y
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx
    do ii = it, threadsE%x+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = ii + offset_yz
        rho(idx) = Q(i,1,j)
          u(idx) = Q(i,2,j)
          v(idx) = Q(i,3,j)
          p(idx) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j) return
    idx = it + offset_yz
    associate(un1 => u(idx), un2 => u(idx+1))
      call SLAU(id_slau, rho(idx), rho(idx+1), u(idx), u(idx+1), v(idx), v(idx+1), &
                un1, un2, p(idx), p(idx+1), Normal_x, 1.0_sp, &
                E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))
    end associate
  end subroutine calc_slau_x2


  !> CUDA Fortran kernel for 2 points SLAU scheme in y direction
  attributes(global) subroutine calc_slau_y2(id_accuracy, nx, ny, Q, sensor, F)
    use mod_constant, only : Normal_y
    integer(2), intent(in), value             :: id_accuracy    !< ID for accuracy, 2 means 2nd-order
    integer, intent(in), value                :: nx             !< number of grid points in x direction
    integer, intent(in), value                :: ny             !< number of grid points in y direction
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny)     !< Q(rho, u, v, p)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)  !< shock sensor
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1) !< Flux in y direction
    integer i,  j  !< global index in physical space
    integer it, jt !< local index in a block
    integer jj, idx, j_base, offset_xz
    integer, parameter :: sx = threadsF%x   !< tile size in x direction
    integer, parameter :: sy = threadsF%y+1 !< tile size in y direction
    real(8), dimension(sx*sy), shared :: rho, u, v, p
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1) * sy
    do jj = jt, threadsF%y+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        idx = jj + offset_xz
        rho(idx) = Q(i,1,j)
          u(idx) = Q(i,2,j)
          v(idx) = Q(i,3,j)
          p(idx) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    idx = jt + offset_xz
    associate(un1 => v(idx), un2 => v(idx+1))
      call SLAU(id_slau, rho(idx), rho(idx+1), u(idx), u(idx+1), v(idx), v(idx+1), &
                un1, un2, p(idx), p(idx+1), Normal_y, 1.0_sp, &
                F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))
    end associate
  end subroutine calc_slau_y2
end module calc_slau_kernel

