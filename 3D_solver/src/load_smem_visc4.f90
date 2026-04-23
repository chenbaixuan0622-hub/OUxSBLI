module load_smem_visc4
  use mod_globals, only : threadsEv, threadsFv, threadsGv
  use mod_constant, only : two_third, one_twelfth
  implicit none
  private
  public load_smem_visc4_x, load_smem_visc4_y, load_smem_visc4_z
contains
  #if _CUDA_ARCH_ >= 900
  !> TMA version
  attributes(device) subroutine load_smem_visc4_x(it, jt, kt, j, k, &
                                                  nx, ny, nz, inv_dy, inv_dz, Q, u, v, w, uy, vy, uz, wz)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: inv_dy(ny-1)  !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous :: inv_dz(nz-1)  !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz) !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsEv%x + 2*io_v + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: uy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: uz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer i_base, ii, i, idx, offset_yz
    ! TMA variables
    integer(8), shared :: barrier
    integer(8) token
    integer start_i, end_i, count, start_idx, tid_linear
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    ! TMA init
    tid_linear = threadIdx%x + (threadIdx%y-1)*blockDim%x &
                + (threadIdx%z-1)*blockDim%x*blockDim%y
    if (tid_linear == 1) then
      call barrier_init(barrier, blockDim%x * blockDim%y * blockDim%z)
    endif
    call syncthreads()
    ! TMA transfer
    start_i = max(1, i_base + 1 - io_v)
    end_i   = min(nx, i_base + threadsEv%x + io_v +1)
    count   = max(0, end_i - start_i + 1)
    start_idx = (start_i - i_base) + offset_yz
    if (tid_linear == 1) then
      if (count > 0 .and. j <= ny .and. k <= nz) then
        call tma_bulk_load(barrier, Q(start_i,2,j,k), u(start_idx), count)
        call tma_bulk_load(barrier, Q(start_i,3,j,k), v(start_idx), count)
        call tma_bulk_load(barrier, Q(start_i,4,j,k), w(start_idx), count)
      endif
    endif
    ! traditional smem load
    do ii = it-io_v, threadsEv%x+io_v+1, blockDim%x
      i = i_base + ii
      idx = ii + offset_yz
      if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2 .and. k <= nz) then
        uy(idx) = (two_third * (-Q(i,2,j-1,k) + Q(i,2,j+1,k)) - one_twelfth * (-Q(i,2,j-2,k) + Q(i,2,j+2,k))) * inv_dy(j)
        vy(idx) = (two_third * (-Q(i,3,j-1,k) + Q(i,3,j+1,k)) - one_twelfth * (-Q(i,3,j-2,k) + Q(i,3,j+2,k))) * inv_dy(j)
      endif
      if (1 <= i .and. i <= nx .and. j <= ny .and. 3 <= k .and. k <= nz-2) then
        uz(idx) = (two_third * (-Q(i,2,j,k-1) + Q(i,2,j,k+1)) - one_twelfth * (-Q(i,2,j,k-2) + Q(i,2,j,k+2))) * inv_dz(k)
        wz(idx) = (two_third * (-Q(i,4,j,k-1) + Q(i,4,j,k+1)) - one_twelfth * (-Q(i,4,j,k-2) + Q(i,4,j,k+2))) * inv_dz(k)
      endif
    enddo
    ! TMA syncthreads
    call syncthreads()
    token = barrier_arrive(barrier)
    do
      if (barrier_try_wait_sleep(barrier, token, 1000000) .ne. 0) exit
    enddo
  end subroutine load_smem_visc4_x
  #else
  attributes(device) subroutine load_smem_visc4_x(it, jt, kt, j, k, &
                                                  nx, ny, nz, inv_dy, inv_dz, Q, u, v, w, uy, vy, uz, wz)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: inv_dy(ny-1)  !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous :: inv_dz(nz-1)  !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz) !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsEv%x + 2*io_v + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: uy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: uz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer i_base, ii, i, idx, offset_yz
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    do ii = it-io_v, threadsEv%x+io_v+1, blockDim%x
      i = i_base + ii
      idx = ii + offset_yz
      if (1 <= i .and. i <= nx .and. j <= ny .and. k <= nz) then
        u(idx) = Q(i,2,j,k)
        v(idx) = Q(i,3,j,k)
        w(idx) = Q(i,4,j,k)
      endif
      if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2 .and. k <= nz) then
        uy(idx) = (two_third * (-Q(i,2,j-1,k) + Q(i,2,j+1,k)) - one_twelfth * (-Q(i,2,j-2,k) + Q(i,2,j+2,k))) * inv_dy(j)
        vy(idx) = (two_third * (-Q(i,3,j-1,k) + Q(i,3,j+1,k)) - one_twelfth * (-Q(i,3,j-2,k) + Q(i,3,j+2,k))) * inv_dy(j)
      endif
      if (1 <= i .and. i <= nx .and. j <= ny .and. 3 <= k .and. k <= nz-2) then
        uz(idx) = (two_third * (-Q(i,2,j,k-1) + Q(i,2,j,k+1)) - one_twelfth * (-Q(i,2,j,k-2) + Q(i,2,j,k+2))) * inv_dz(k)
        wz(idx) = (two_third * (-Q(i,4,j,k-1) + Q(i,4,j,k+1)) - one_twelfth * (-Q(i,4,j,k-2) + Q(i,4,j,k+2))) * inv_dz(k)
      endif
    enddo
  end subroutine load_smem_visc4_x
  #endif


  attributes(device) subroutine load_smem_visc4_y(it, jt, kt, i, k, &
                                                  nx, ny, nz, inv_dx, inv_dz, Q, u, v, w, ux, vx, vz, wz)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: inv_dx(nx-1)  !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous :: inv_dz(nz-1)  !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz) !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 2*io_v + 1
    integer, parameter :: sz = threadsFv%z
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: ux(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vx(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer j_base, jj, j, idx, offset_xz
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    do jj = jt-io_v, threadsFv%y+io_v+1, blockDim%y
      j = j_base + jj
      idx = jj + offset_xz
      if (i <= nx .and. 1 <= j .and. j <= ny .and. k <= nz) then
        u(idx) = Q(i,2,j,k)
        v(idx) = Q(i,3,j,k)
        w(idx) = Q(i,4,j,k)
      endif
      if (3 <= i .and. i <= nx-2 .and. 1 <= j .and. j <= ny .and. k <= nz) then
        ux(idx) = (two_third * (-Q(i-1,2,j,k) + Q(i+1,2,j,k)) - one_twelfth * (-Q(i-2,2,j,k) + Q(i+2,2,j,k))) * inv_dx(i)
        vx(idx) = (two_third * (-Q(i-1,3,j,k) + Q(i+1,3,j,k)) - one_twelfth * (-Q(i-2,3,j,k) + Q(i+2,3,j,k))) * inv_dx(i)
      endif
      if (i <= nx .and. 1 <= j .and. j <= ny .and. 3 <= k .and. k <= nz-2) then
        vz(idx) = (two_third * (-Q(i,3,j,k-1) + Q(i,3,j,k+1)) - one_twelfth * (-Q(i,3,j,k-2) + Q(i,3,j,k+2))) * inv_dz(k)
        wz(idx) = (two_third * (-Q(i,4,j,k-1) + Q(i,4,j,k+1)) - one_twelfth * (-Q(i,4,j,k-2) + Q(i,4,j,k+2))) * inv_dz(k)
      endif
    enddo
  end subroutine load_smem_visc4_y

  
  attributes(device) subroutine load_smem_visc4_z(it, jt, kt, i, j, &
                                                  nx, ny, nz, inv_dx, inv_dy, Q, u, v, w, ux, wx, vy, wy)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: inv_dx(nx-1)  !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous :: inv_dy(ny-1)  !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz) !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 2*io_v + 1
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: ux(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wx(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer k_base, kk, k, idx, offset_xy
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    do kk = kt-io_v, threadsGv%z+io_v+1, blockDim%z
      k = k_base + kk
      idx = kk + offset_xy
      if (i <= nx .and. j <= ny .and. 1 <= k .and. k <= nz) then
        u(idx) = Q(i,2,j,k)
        v(idx) = Q(i,3,j,k)
        w(idx) = Q(i,4,j,k)
      endif
      if (3 <= i .and. i <= nx-2 .and. j <= ny .and. 1 <= k .and. k <= nz) then
        ux(idx) = (two_third * (-Q(i-1,2,j,k) + Q(i+1,2,j,k)) - one_twelfth * (-Q(i-2,2,j,k) + Q(i+2,2,j,k))) * inv_dx(i)
        wx(idx) = (two_third * (-Q(i-1,4,j,k) + Q(i+1,4,j,k)) - one_twelfth * (-Q(i-2,4,j,k) + Q(i+2,4,j,k))) * inv_dx(i)
      endif
      if (i <= nx .and. 3 <= j .and. j <= ny-2 .and. 1 <= k .and. k <= nz) then
        vy(idx) = (two_third * (-Q(i,3,j-1,k) + Q(i,3,j+1,k)) - one_twelfth * (-Q(i,3,j-2,k) + Q(i,3,j+2,k))) * inv_dy(j)
        wy(idx) = (two_third * (-Q(i,4,j-1,k) + Q(i,4,j+1,k)) - one_twelfth * (-Q(i,4,j-2,k) + Q(i,4,j+2,k))) * inv_dy(j)
      endif
    enddo
  end subroutine load_smem_visc4_z
end module load_smem_visc4

