module load_smem_visc4
  use wmma
  use mod_globals, only : threadsEv, threadsFv
  use mod_constant, only : two_third, one_twelfth
  implicit none
  private
  public load_smem_visc4_x, load_smem_visc4_y
contains
  attributes(device) subroutine load_smem_visc4_x(it, jt, j, &
                                                  nx, ny, inv_dy, Q, u, v, uy, vy)
    integer, intent(in), value              :: it           !< local idx for x direction
    integer, intent(in), value              :: jt           !< local idx for y direction
    integer, intent(in), value              :: j            !< global idx for y direction
    integer, intent(in), value              :: nx           !< number of grid points in x direction
    integer, intent(in), value              :: ny           !< number of grid points in y direction
    real(8), intent(in), device, contiguous :: inv_dy(ny-1) !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous :: Q(nx,4,ny)   !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsEv%x + 2*io_v + 1
    integer, parameter :: sy = threadsEv%y
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(8), intent(inout) :: uy(-(io_v-1):sx*sy-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vy(-(io_v-1):sx*sy-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer i_base, ii, i, idx, offset_y
    i_base = (blockIdx%x-1)*blockDim%x
    offset_y = (jt-1)*sx
    do ii = it-io_v, threadsEv%x+io_v+1, blockDim%x
      i = i_base + ii
      idx = ii + offset_y
      if (1 <= i .and. i <= nx .and. j <= ny) then
        call pipelineMemcpyAsync(u(idx), Q(i,2,j))
        call pipelineMemcpyAsync(v(idx), Q(i,3,j))
      endif
    enddo
    call pipelineCommit()
    do ii = it-io_v, threadsEv%x+io_v+1, blockDim%x
      i = i_base + ii
      idx = ii + offset_y
      if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2) then
        uy(idx) = (two_third * (-Q(i,2,j-1) + Q(i,2,j+1)) - one_twelfth * (-Q(i,2,j-2) + Q(i,2,j+2))) * inv_dy(j)
        vy(idx) = (two_third * (-Q(i,3,j-1) + Q(i,3,j+1)) - one_twelfth * (-Q(i,3,j-2) + Q(i,3,j+2))) * inv_dy(j)
      endif
    enddo
    call pipelineWaitPrior(0)
    call syncthreads()
  end subroutine load_smem_visc4_x


  attributes(device) subroutine load_smem_visc4_y(it, jt, i, &
                                                  nx, ny, inv_dx, Q, u, v, ux, vx)
    integer, intent(in), value              :: it           !< local idx for x direction
    integer, intent(in), value              :: jt           !< local idx for y direction
    integer, intent(in), value              :: i            !< global idx for x direction
    integer, intent(in), value              :: nx           !< number of grid points in x direction
    integer, intent(in), value              :: ny           !< number of grid points in y direction
    real(8), intent(in), device, contiguous :: inv_dx(nx-1) !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous :: Q(nx,4,ny)   !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 2*io_v + 1
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(8), intent(inout) :: ux(-(io_v-1):sx*sy-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vx(-(io_v-1):sx*sy-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer j_base, jj, j, idx, offset_x
    j_base = (blockIdx%y-1)*blockDim%y
    offset_x = (it-1)*sy
    do jj = jt-io_v, threadsFv%y+io_v+1, blockDim%y
      j = j_base + jj
      idx = jj + offset_x
      if (i <= nx .and. 1 <= j .and. j <= ny) then
        call pipelineMemcpyAsync(u(idx), Q(i,2,j))
        call pipelineMemcpyAsync(v(idx), Q(i,3,j))
      endif
    enddo
    call pipelineCommit()
    do jj = jt-io_v, threadsFv%y+io_v+1, blockDim%y
      j = j_base + jj
      idx = jj + offset_x
      if (3 <= i .and. i <= nx-2 .and. 1 <= j .and. j <= ny) then
        ux(idx) = (two_third * (-Q(i-1,2,j) + Q(i+1,2,j)) - one_twelfth * (-Q(i-2,2,j) + Q(i+2,2,j))) * inv_dx(i)
        vx(idx) = (two_third * (-Q(i-1,3,j) + Q(i+1,3,j)) - one_twelfth * (-Q(i-2,3,j) + Q(i+2,3,j))) * inv_dx(i)
      endif
    enddo
    call pipelineWaitPrior(0)
    call syncthreads()
  end subroutine load_smem_visc4_y
end module load_smem_visc4

