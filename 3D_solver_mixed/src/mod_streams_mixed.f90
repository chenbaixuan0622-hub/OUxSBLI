!> Two-stream infrastructure for concurrent convective (real8) / viscous
!> (real4-internal) kernel execution, plus the cross-phase events used to
!> order the staggered E/F/G buffer reuse in calc_flux_base_mixed.f90.fypp.
!> stream_conv/stream_visc are each bound to their own CUDA Green Context
!> (green_ctx_bindings_mixed.f90) so the two streams run on disjoint,
!> hardware-partitioned SM sets -- a spatial guarantee ordinary CUDA stream
!> concurrency does not provide. Without it, a large enough grid on either
!> stream can opportunistically consume SMs "intended" for the other stream,
!> collapsing conv/visc overlap back toward serial execution. Adapted from
!> 3D_solver_fltflt/src/mod_streams_ff.f90 -- the SM-split/Green-Context
!> mechanism carries over unchanged; only the fltflt-constant-init calls at
!> the end of init_streams are dropped (nothing here uses the fltflt type,
!> convection stays plain real(8)).
module mod_streams_mixed
  use cudafor
  use iso_c_binding
  use green_ctx_bindings
  implicit none
  integer(cuda_stream_kind) :: stream_conv, stream_visc
  type(cudaEvent) :: event_E, event_Fv, event_Gv
  integer :: numSM
  ! Measured from a real Miyabi GH200 (132 SMs) nsys trace of this fork's
  ! NSTGV kernels at the production 513^3 grid (OVERLAP_GC=False run):
  ! summing full-132-SM per-kernel costs into "SM-ms" work gives
  ! W_conv = (Ducros+E+F+G) * 132 =~ 2710 SM-ms, W_visc = (Ev+Fv+Gv) * 132
  ! =~ 4093 SM-ms -- i.e. the viscous side is the HEAVIER total workload,
  ! the opposite of what the inherited fltflt-fork 0.6 (conv-favoring) split
  ! assumed. Under perfect overlap, wall time = max(W_conv/x, W_visc/(132-x)),
  ! minimized (both terms equal) at x =~ 132*W_conv/(W_conv+W_visc) =~ 53 SMs
  ! to conv / ~79 to visc, i.e. conv_share =~ 53/132 =~ 0.40.
  ! IMPORTANT: this is a work-conservation optimum, not a guarantee of a net
  ! win. At that balance point, best-case wall time = (W_conv+W_visc)/132,
  ! i.e. EXACTLY the sequential (OVERLAP_GC=False) total -- splitting a fixed
  ! SM pool into two concurrent partitions can only reallocate existing
  ! compute capacity, not add to it, and both conv and visc kernels here
  ! scale ~linearly with SM count at this grid size (confirmed in the trace:
  ! per-kernel slowdown when confined to a partition tracks the SM-count
  ! ratio closely). So even at this optimal split, expect OVERLAP_GC=True to
  ! at best TIE plain sequential dispatch, never decisively beat it -- do not
  ! re-tune this expecting a win; a real gain would need a workload where at
  ! least one side is latency/occupancy-bound rather than SM-count-bound.
  real, parameter :: conv_share = 0.4
contains
  subroutine init_streams(mygpu)
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
    ! over from that rounding go unused rather than into either context, per
    ! cuDevSmResourceSplitByCount's documented behavior.
    nb_conv_sm = round_to_multiple8(int(conv_share * numSM), numSM)
    cires = cuInit(0)
    cires = cuDeviceGet(dev, mygpu)
    cires = cuDeviceGetDevResource(dev, whole_sm, CU_DEV_RESOURCE_TYPE_SM)
    nbGroups = 1
    ! "1 group + remaining" pattern: sm_parts(1) gets exactly nb_conv_sm SMs,
    ! sm_parts(2) captures whatever is left. Loop bounds below are the fixed
    ! i_conv/i_visc indices, not nbGroups -- nbGroups is only the *first*
    ! split call's group count (1 here, since the second partition comes
    ! from `remaining`), and keying a creation loop off it would silently
    ! skip the second green context.
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
  end subroutine init_streams

  function round_to_multiple8(target, total) result(n)
    integer, intent(in) :: target, total
    integer :: n
    n = max(8, min(((target + 4) / 8) * 8, total - 8))
  end function round_to_multiple8
end module mod_streams_mixed
