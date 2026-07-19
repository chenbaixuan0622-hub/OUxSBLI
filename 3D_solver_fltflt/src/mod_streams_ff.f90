!> Two-stream infrastructure for concurrent convective (fltflt) / viscous (real8)
!> kernel execution: per-kernel occupancy-driven block counts and the
!> cross-phase events used to order the staggered E/F/G buffer reuse in
!> calc_flux_base. Block counts are derived from each kernel's *actual*
!> max resident blocks/SM (queried via cudaOccupancyMaxActiveBlocksPerMultiprocessor
!> in calc_flux_base's init_occupancy, cached there) rather than assumed --
!> a flat "one block per SM" guess badly under-provisions parallelism for
!> these memory/latency-bound stencil kernels (unlike the compute-bound,
!> register-heavy kernel in ~/TMA_practice/practice5/kernel8_fltflt.f90
!> that guess was originally modeled on).
module mod_streams
  use cudafor
  implicit none
  integer(cuda_stream_kind) :: stream_conv, stream_visc
  type(cudaEvent) :: event_E, event_Fv, event_Gv
  integer :: numSM
contains
  subroutine init_streams(mygpu)
    use mod_constant, only : init_fltflt_constants
    use calc_keep_kernel_internal, only : init_keep_ff_constants
    integer, intent(in) :: mygpu
    type(cudaDeviceProp) :: prop
    integer :: istat
    istat = cudaGetDeviceProperties(prop, mygpu)
    numSM = prop%multiProcessorCount
    print '(1x, a, i0)', "SMs: ", numSM
    istat = cudaStreamCreate(stream_conv)
    istat = cudaStreamCreate(stream_visc)
    istat = cudaEventCreate(event_E)
    istat = cudaEventCreate(event_Fv)
    istat = cudaEventCreate(event_Gv)
    call init_fltflt_constants()
    call init_keep_ff_constants()
  end subroutine init_streams

  !> Reserve half of a kernel's true max-resident-blocks-per-SM capacity
  !> (numBlocksPerSM, from cudaOccupancyMaxActiveBlocksPerMultiprocessor)
  !> for this stream, leaving the other half free for the paired stream.
  function query_nblk(numBlocksPerSM) result(nblk)
    integer, intent(in) :: numBlocksPerSM
    integer :: nblk
    nblk = max(1, (numBlocksPerSM * numSM) / 2)
  end function query_nblk
end module mod_streams
