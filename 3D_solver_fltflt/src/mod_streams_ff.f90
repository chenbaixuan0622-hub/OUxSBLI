!> Two-stream infrastructure for concurrent convective (fltflt) / viscous (real8)
!> kernel execution, plus the cross-phase events used to order the staggered
!> E/F/G buffer reuse in calc_flux_base. stream_conv/stream_visc are each
!> bound to their own CUDA Green Context (green_ctx_bindings_ff.f90) so the
!> two streams run on disjoint, hardware-partitioned SM sets -- a spatial
!> guarantee that ordinary CUDA stream concurrency does not provide. Without
!> it, a large enough grid on either stream can opportunistically consume
!> SMs "intended" for the other stream, collapsing conv/visc overlap back
!> toward serial execution; shrinking the grid to compensate just reintroduces
!> the low-occupancy problem this was meant to avoid. Both kernels now launch
!> at their natural full grid size (calc_flux_base_ff.f90.fypp), same as the
!> baseline (non-split) solver -- no persistent-kernel tiling needed once the
!> SM split is enforced by the green context rather than by grid size.
module mod_streams
  use cudafor
  use iso_c_binding
  use green_ctx_bindings
  implicit none
  integer(cuda_stream_kind) :: stream_conv, stream_visc
  type(cudaEvent) :: event_E, event_Fv, event_Gv
  integer :: numSM
  ! Not a flat 50/50: an nsys profile (NSTGV, GH200, one RK3 stage) measured
  ! the convective (fltflt) stream's total sequential kernel time at
  ! ~120ms/stage (E+F+G) vs. the viscous (real8) stream's ~99ms/stage
  ! (Fv+Gv+Ev) -- conv is the critical path each stage, visc finishes with
  ! ~21% slack. Give conv a larger share of the SMs so the two streams'
  ! green-context partitions finish closer together instead of visc idling
  ! while conv catches up. This ratio is a starting point derived from that
  ! one profile, not a universal constant -- re-derive from fresh profiling
  ! if the scheme, order, or grid size change materially enough to shift the
  ! conv/visc time ratio.
  real, parameter :: conv_share = 0.6
contains
  subroutine init_streams(mygpu)
    use mod_constant, only : init_fltflt_constants
    use calc_keep_kernel_internal, only : init_keep_ff_constants
    integer, intent(in) :: mygpu
    integer, parameter :: i_conv = 1, i_visc = 2
    type(cudaDeviceProp) :: prop
    integer :: istat
    integer(c_int) :: cires, dev, nbGroups, nb_conv_sm
    type(CUdevResource) :: whole_sm, sm_parts(2)
    type(c_ptr) :: desc(2), gctx(2), cstream(2)
    integer :: i
    istat = cudaGetDeviceProperties(prop, mygpu)
    numSM = prop%multiProcessorCount
    print '(1x, a, i0)', "SMs: ", numSM

    ! Partition SMs into disjoint conv/visc green contexts. minCount is
    ! rounded to a multiple of 8 (valid on both the cc9.0+ rule -- min 8,
    ! multiple of 8 -- and the coarser cc8.x rule -- min 4, multiple of 2 --
    ! since every multiple of 8 also satisfies the latter); any SMs left
    ! over from that rounding (e.g. 132 SMs -> 80+48, 4 unused) go unused
    ! rather than into either context, per cuDevSmResourceSplitByCount's
    ! documented behavior.
    nb_conv_sm = round_to_multiple8(int(conv_share * numSM), numSM)
    cires = cuInit(0)
    cires = cuDeviceGet(dev, mygpu)
    cires = cuDeviceGetDevResource(dev, whole_sm, CU_DEV_RESOURCE_TYPE_SM)
    nbGroups = 1
    ! "1 group + remaining" pattern (validated in the Stage 1 prototype):
    ! sm_parts(1) gets exactly nb_conv_sm SMs, sm_parts(2) captures whatever
    ! is left. Loop bounds below are the fixed i_conv/i_visc indices, not
    ! nbGroups -- nbGroups is only the *first* split call's group count (1
    ! here, since the second partition comes from `remaining`), and keying a
    ! creation loop off it would silently skip the second green context.
    cires = cuDevSmResourceSplitByCount(sm_parts(i_conv:i_conv), nbGroups, whole_sm, sm_parts(i_visc), 0, nb_conv_sm)
    print '(1x, a, i0, a, i0)', "Green Context SM split -- conv: ", sm_count(sm_parts(i_conv)), &
                                 "  visc: ", sm_count(sm_parts(i_visc))

    do i = 1, 2
      cires = cuDevResourceGenerateDesc(desc(i), sm_parts(i:i), 1)
      cires = cuGreenCtxCreate(gctx(i), desc(i), dev, CU_GREEN_CTX_DEFAULT_STREAM)
      cires = cuGreenCtxStreamCreate(cstream(i), gctx(i), CU_STREAM_NON_BLOCKING, 0)
    end do
    stream_conv = transfer(cstream(i_conv), stream_conv)
    stream_visc = transfer(cstream(i_visc), stream_visc)

    istat = cudaEventCreate(event_E)
    istat = cudaEventCreate(event_Fv)
    istat = cudaEventCreate(event_Gv)
    call init_fltflt_constants()
    call init_keep_ff_constants()
  end subroutine init_streams

  function round_to_multiple8(target, total) result(n)
    integer, intent(in) :: target, total
    integer :: n
    n = max(8, min(((target + 4) / 8) * 8, total - 8))
  end function round_to_multiple8
end module mod_streams
