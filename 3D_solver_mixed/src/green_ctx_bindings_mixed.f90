!> Hand-written bind(C) interfaces to the CUDA *driver* API's green-context
!> functions (cuda.h, CUDA 12.4+), used by mod_streams_mixed.f90 to give
!> stream_conv/stream_visc disjoint, hardware-guaranteed SM partitions
!> instead of sharing the whole GPU via a grid-size heuristic. nvfortran's
!> `cudafor` module wraps the *runtime* API only and does not expose these.
!> Verbatim copy of 3D_solver_fltflt/src/green_ctx_bindings_ff.f90 -- this
!> binding layer is precision-agnostic (it only partitions SMs and creates
!> streams), so nothing about it changes for the FP64-conv/FP32-visc split.
!>
!> CUdevResource (cuda.h): a tagged union --
!>   struct { CUdevResourceType type; unsigned char _internal_padding[92];
!>            union { CUdevSmResource sm; unsigned char _oversize[48]; }; }
!> = 4 + 92 + 48 = 144 bytes. Modeled here as a byte buffer since Fortran has
!> no native union; the first 4 bytes of `payload` double as the `type` tag,
!> and (when type == CU_DEV_RESOURCE_TYPE_SM) the union's first 4 bytes
!> (offset 96 within the struct, i.e. payload(1:4)) hold `sm%smCount`.
module green_ctx_bindings
  use iso_c_binding
  implicit none

  integer(c_int), parameter :: CU_DEV_RESOURCE_TYPE_SM       = 1
  integer(c_int), parameter :: CU_GREEN_CTX_DEFAULT_STREAM   = 1  ! 0x1
  integer(c_int), parameter :: CU_STREAM_NON_BLOCKING        = 1  ! 0x1

  type, bind(C) :: CUdevResource
    integer(c_int)    :: restype
    character(c_char) :: reserved_pad(92)
    character(c_char) :: payload(48)
  end type CUdevResource

  interface
    ! CUresult cuInit(unsigned int Flags)
    function cuInit(flags) bind(C, name="cuInit") result(res)
      import :: c_int
      integer(c_int), value :: flags
      integer(c_int) :: res
    end function cuInit

    ! CUresult cuDeviceGet(CUdevice *device, int ordinal)
    function cuDeviceGet(device, ordinal) bind(C, name="cuDeviceGet") result(res)
      import :: c_int
      integer(c_int), intent(out) :: device
      integer(c_int), value :: ordinal
      integer(c_int) :: res
    end function cuDeviceGet

    ! CUresult cuDeviceGetDevResource(CUdevice device, CUdevResource *resource, CUdevResourceType type)
    function cuDeviceGetDevResource(device, resource, restype) &
        bind(C, name="cuDeviceGetDevResource") result(res)
      import :: c_int, CUdevResource
      integer(c_int), value :: device
      type(CUdevResource), intent(out) :: resource
      integer(c_int), value :: restype
      integer(c_int) :: res
    end function cuDeviceGetDevResource

    ! CUresult cuDevSmResourceSplitByCount(CUdevResource *result, unsigned int *nbGroups,
    !     const CUdevResource *input, CUdevResource *remaining, unsigned int useFlags, unsigned int minCount)
    function cuDevSmResourceSplitByCount(result_arr, nbGroups, input, remaining, useFlags, minCount) &
        bind(C, name="cuDevSmResourceSplitByCount") result(res)
      import :: c_int, CUdevResource
      type(CUdevResource), intent(out) :: result_arr(*)
      integer(c_int), intent(inout) :: nbGroups
      type(CUdevResource), intent(in) :: input
      type(CUdevResource), intent(out) :: remaining
      integer(c_int), value :: useFlags
      integer(c_int), value :: minCount
      integer(c_int) :: res
    end function cuDevSmResourceSplitByCount

    ! CUresult cuDevResourceGenerateDesc(CUdevResourceDesc *phDesc, CUdevResource *resources, unsigned int nbResources)
    function cuDevResourceGenerateDesc(phDesc, resources, nbResources) &
        bind(C, name="cuDevResourceGenerateDesc") result(res)
      import :: c_int, c_ptr, CUdevResource
      type(c_ptr), intent(out) :: phDesc
      type(CUdevResource), intent(in) :: resources(*)
      integer(c_int), value :: nbResources
      integer(c_int) :: res
    end function cuDevResourceGenerateDesc

    ! CUresult cuGreenCtxCreate(CUgreenCtx* phCtx, CUdevResourceDesc desc, CUdevice dev, unsigned int flags)
    function cuGreenCtxCreate(phCtx, desc, dev, flags) bind(C, name="cuGreenCtxCreate") result(res)
      import :: c_int, c_ptr
      type(c_ptr), intent(out) :: phCtx
      type(c_ptr), value :: desc
      integer(c_int), value :: dev
      integer(c_int), value :: flags
      integer(c_int) :: res
    end function cuGreenCtxCreate

    ! CUresult cuGreenCtxStreamCreate(CUstream* phStream, CUgreenCtx greenCtx, unsigned int flags, int priority)
    function cuGreenCtxStreamCreate(phStream, greenCtx, flags, priority) &
        bind(C, name="cuGreenCtxStreamCreate") result(res)
      import :: c_int, c_ptr
      type(c_ptr), intent(out) :: phStream
      type(c_ptr), value :: greenCtx
      integer(c_int), value :: flags
      integer(c_int), value :: priority
      integer(c_int) :: res
    end function cuGreenCtxStreamCreate

    ! CUresult cuGreenCtxDestroy(CUgreenCtx hCtx)
    function cuGreenCtxDestroy(hCtx) bind(C, name="cuGreenCtxDestroy") result(res)
      import :: c_int, c_ptr
      type(c_ptr), value :: hCtx
      integer(c_int) :: res
    end function cuGreenCtxDestroy
  end interface

contains

  function sm_count(resource) result(n)
    type(CUdevResource), intent(in) :: resource
    integer(c_int) :: n
    n = transfer(resource%payload(1:4), 0_c_int)
  end function sm_count

end module green_ctx_bindings
