# Optimization Plan: OUxSBLI GPU-Accelerated CFD Solver

## Context

The user wants to optimize the OUxSBLI CUDA Fortran CFD solver for speed. The code solves compressible Navier-Stokes equations using explicit high-order finite-difference schemes on NVIDIA GPUs. The hot path is a TVD-RK3 loop: each of the 3 stages calls `calc_EFG` (flux kernels) → `calc_step*` (state update) → `set_bc` (boundary conditions).

A thorough survey of all source files was performed. Findings below are grounded in the actual code.

---

## What is Already Well-Optimized

- `contiguous` on all `calc_flux_base.f90` wrapper arguments ✓
- `_internal` kernel variants (skip boundary branch divergence) already integrated for KEEP and SLAU in `calc_flux_base.f90:49–101` ✓
- `__shfl_down_sync` warp shuffle in `calc_steps.f90` for RK update ✓
- Jacobians/grid metrics pre-computed and stored on device ✓
- No heap allocations inside the RK inner loop ✓
- Thread block sizes already manually tuned per case ✓

---

## Tier 1: Low Effort, High Confidence (~1–2 hours each)

### 1A. Add `contiguous` to CUDA kernel dummy arguments

**Finding:** `calc_flux_base.f90` (the wrapper) already has `contiguous` on its arguments, but the individual CUDA kernel subroutines inside `calc_keep_kernel.f90`, `calc_slau_kernel.f90`, `calc_roe_kernel.f90`, `calc_hybrid_kernel.f90` do **not**. The compiler cannot assume contiguity across the call boundary without this, so it generates conservative gather/scatter loads.

**Change:** In each `attributes(global)` subroutine, add `contiguous` to all `device` array arguments:

```fortran
! Before (e.g. calc_keep_kernel.f90:39-41):
real(8), intent(in), device   :: Q(5,nx,ny,nz)
real(8), intent(in), device   :: T(nx,ny,nz)
real(8), intent(out), device  :: E(5,nx-1,ny-2,nz-2)

! After:
real(8), intent(in), device, contiguous   :: Q(5,nx,ny,nz)
real(8), intent(in), device, contiguous   :: T(nx,ny,nz)
real(8), intent(out), device, contiguous  :: E(5,nx-1,ny-2,nz-2)
```

**Files to modify:**
- `3D_solver/src/calc_keep_kernel.f90` — all 9 subroutines (`calc_keep_x6/4/2`, `y6/4/2`, `z6/4/2`)
- `3D_solver/src/calc_keep_kernel_internal.f90` — `calc_keep_x/y/z_in`
- `3D_solver/src/calc_slau_kernel.f90` — all 9 subroutines + `sensor(nx,ny,nz)` arg
- `3D_solver/src/calc_slau_kernel_internal.f90` — `calc_slau_x/y/z_in`
- `3D_solver/src/calc_roe_kernel.f90` — all 9 subroutines + `sensor(nx,ny,nz)` arg
- `3D_solver/src/calc_hybrid_kernel.f90` — all 9 subroutines + `sensor(nx,ny,nz)` arg

**Verify:** Build with `-Minfo=accel` and confirm the compiler report no longer shows "copy" operations for device arrays; load patterns should show vectorized or non-strided access.

**Estimated impact:** 3–8% per convective kernel. Largest benefit in 6th-order KEEP/SLAU where 6-point stencil reads 6 components from Q.

---

### 1B. Persist the `sensor` array (avoid repeated device allocation)

**Finding:** `calc_conv_slau`, `calc_conv_roe`, and `calc_conv_hybrid` in `calc_flux_base.f90` each declare:
```fortran
real(8), device :: sensor(nx,ny,nz)   ! line 85, 124, 151
```
This is a local variable that gets allocated and deallocated from device memory on every call — 3 times per RK stage, hundreds of times per time step.

**Change:** Move `sensor` to a module-level `save` variable in `calc_flux_base.f90`, allocated once during initialization and reused:

```fortran
! At module level in calc_flux_base:
real(8), allocatable, device, save :: sensor(:,:,:)

! In a new init_flux_base(nx,ny,nz) routine called from preprocess.f90:
allocate(sensor(nx,ny,nz))
```

Remove the local `real(8), device :: sensor(nx,ny,nz)` declarations from the three `calc_conv_*` subroutines. Call `init_flux_base` from `preprocess.f90` before the time loop.

**Files to modify:**
- `3D_solver/src/calc_flux_base.f90` — add module-level sensor, remove 3 local declarations, add init routine
- `3D_solver/src/preprocess.f90` — call `init_flux_base(nx, ny, nz)` during setup

**Estimated impact:** Eliminates ~9 device heap allocations per time step (3 RK stages × 3 schemes). On large grids (SBLI: 513×161×33 = 2.7M points × 8 bytes = 21 MB sensor), repeated allocation has measurable overhead in device memory allocator.

---

### 1C. Enable GPU-direct MPI (SBLI/TBL cases with rescaling)

**Finding:** `id_gpumpi` in `mod_globals.f90` uses Fortran kind-dispatch:
- `integer(kind=2)` → CPU path: allocates temp host buffer, D→H cudaMemcpy, MPI_SEND, free (repeated every RK stage during SBLI rescaling)
- `integer(kind=4)` → GPU-direct: passes device pointer directly to MPI

**Change:** In `3D_solver/SBLI/mod_globals.f90` and `3D_solver/TBL/mod_globals.f90`:
```fortran
! Before:
integer(kind=2), parameter :: id_gpumpi = 0
! After:
integer(kind=4), parameter :: id_gpumpi = 0
```

**Prerequisite check:** Verify MPI is CUDA-aware — `ompi_info | grep cuda` should show `true`. If using NVIDIA HPC SDK's bundled MPI (HPCX), this is typically available.

**Estimated impact:** 5–15% wall-clock reduction for SBLI with `id_rescale = kind4`. Eliminates CPU staging buffer and two PCIe round-trips per rescaling MPI transfer.

---

## Tier 2: Medium Effort, Significant Impact (~1–3 days each)

### 2A. CUDA streams to overlap E, F, G flux kernels

**Finding:** In `calc_conv_keep/slau/roe/hybrid`, the E, F, G kernels are launched sequentially to the default stream. These three kernels are **fully independent**: they all read `Q` (read-only) and write to separate `E`, `F`, `G` arrays with no overlap. They are safe to run concurrently.

**Design:**

Add stream handles to `calc_flux_base.f90` module state:
```fortran
integer, save :: stream_E, stream_F, stream_G
```

Initialize in `init_flux_base`:
```fortran
istat = cudaStreamCreate(stream_E)
istat = cudaStreamCreate(stream_F)
istat = cudaStreamCreate(stream_G)
```

In each `calc_conv_*` subroutine, pass stream as 4th `<<<>>>` argument:
```fortran
call calc_keep_x_in<<<blocksE,threadsE,1,stream_E>>>(...)
call calc_keep_y_in<<<blocksF,threadsF,2,stream_F>>>(...)
call calc_keep_z_in<<<blocksG,threadsG,3,stream_G>>>(...)
istat = cudaDeviceSynchronize()  ! after all 3, before caller uses E/F/G
```

For Roe/Hybrid with Ducros sensor: `calc_Ducros` must complete before E/F/G kernels read `sensor`. Use a CUDA event to express this dependency without a global sync:
```fortran
call calc_Ducros<<<blocks,threads,0,stream_E>>>(...)  ! use stream_E for Ducros
istat = cudaEventRecord(ducros_done, stream_E)
istat = cudaStreamWaitEvent(stream_F, ducros_done, 0)
istat = cudaStreamWaitEvent(stream_G, ducros_done, 0)
! Now launch all three:
call calc_roe_x<<<blocksE,threadsE,1,stream_E>>>(...)
call calc_roe_y<<<blocksF,threadsF,2,stream_F>>>(...)
call calc_roe_z<<<blocksG,threadsG,3,stream_G>>>(...)
istat = cudaDeviceSynchronize()
```

The viscous kernels (`calc_Ev4`, `calc_Fv4`, `calc_Gv4`) can also use the same three streams after the `cudaDeviceSynchronize` — the existing `stat = cudaDeviceSynchronize()` at line 224 of `calc_EFG_visc` already provides the necessary fence.

**Files to modify:**
- `3D_solver/src/calc_flux_base.f90` — add stream/event handles, modify all 4 `calc_conv_*` subroutines and the viscous section of `calc_EFG_visc` and `calc_EFG_LES`
- `3D_solver/src/preprocess.f90` — call stream/event creation

**Estimated impact:** Up to 33% reduction in flux computation time if all three direction kernels fully overlap. In practice, 15–25% on Ada (RTX 4090) with adequate SM occupancy. Largest for symmetric grids (NSTGV 513³ where all three kernels are equal size). Profile with nsys to confirm actual overlap.

---

### 2B. Extend `_internal` kernel pattern to Roe and Hybrid

**Finding:** `calc_conv_roe` (lines 126–128) and `calc_conv_hybrid` (lines 153–155) call kernels unconditionally — no `if (id_bc_x)` guard, no `_internal` variant. Every thread pays the cost of the 3-branch if-elseif-else stencil fallback (4th-order interior → 2nd-order near boundary → 1st-order at boundary) even for interior threads.

**Design:** Create `3D_solver/src/calc_roe_kernel_internal.f90` following `calc_slau_kernel_internal.f90` as template:
- Module name: `calc_roe_kernel_internal`, public procedures: `calc_roe_x_in`, `calc_roe_y_in`, `calc_roe_z_in`
- Use `integer, parameter :: io = kind(id_accuracy) / 3` for compile-time stencil half-width (2→0, 4→1, 8→2)
- Only compute output where `io+1 <= i <= nx-(io+1)` — eliminates all boundary branching
- Preserve the Roe kernel's 2D shaped shared memory: `dimension(-1:threadsE%x+3, threadsE%y, threadsE%z)` — do NOT switch to 1D linearized layout like KEEP-internal uses
- Keep the two-phase `syncthreads()` pattern (after fill, after left-state write) from the regular Roe kernel

Update `calc_flux_base.f90` for Roe (mirror the KEEP/SLAU pattern at lines 49–63):
```fortran
! In calc_conv_roe, replace lines 126-128 with:
if (id_bc_x) then
  call calc_roe_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, Q, sensor, E)
else
  call calc_roe_x_in<<<blocksE,threadsE,1>>>(nx, ny, nz, Q, sensor, E)
endif
! repeat for y and z
```

Add to `USE` list: `use calc_roe_kernel_internal` in `calc_flux_base.f90`.

Add to Makefile `OBJ` in each case that uses Roe:
```
calc_roe_kernel_internal.o \
```
With dependency:
```
calc_roe_kernel_internal.o: mod_globals.mod mod_constant.mod calc_muscl.mod calc_hybrid.mod calc_roe_3d.f90
```

Do the same for `calc_hybrid_kernel_internal.f90`. Hybrid is more complex since it branches on `fdx <= threshold` (physics-driven by Ducros sensor), but the stencil-order boundary fallback can still be eliminated for the interior.

**Estimated impact:** 5–15% reduction in convective kernel time for 4th and 6th-order schemes with many interior points. Most noticeable for NSTGV (513³, all periodic — currently all directions use `_internal` for KEEP/SLAU but fall back to boundary kernel for Roe/Hybrid).

---

### 2C. Asynchronous I/O (overlap output with compute)

**Finding:** Every `np` outer steps, `send_recv_for_print_even3` in `calc_time_dev.f90` performs:
1. `Q = QJ` — blocking D→H copy of full `5×nx×ny×nz` doubles (e.g. SBLI: ~4 GB)
2. `make_1d_for_print` — real(8)→real(4) conversion on CPU
3. `MPI_WAITALL` — blocks until I/O rank has received all data

GPU cannot proceed during step 1. This is on the critical path.

**Step 1 (pin host Q array):** Change the host Q allocation in `calc_time_dev.f90` to pinned memory to enable DMA transfers:
```fortran
! Before (somewhere in preprocess.f90/ calc_time_dev.f90):
allocate(Q(5,nx,ny,nz))
! After:
allocate(Q(5,nx,ny,nz), pinned=.true.)
! (CUDA Fortran syntax: allocate(arr, pinned=.true.) or use cudaHostAlloc)
```
Pinned memory enables 2–3× faster PCIe transfers at no cost to correctness.

**Step 2 (async double-buffer):** Instead of blocking `Q = QJ`, use:
```fortran
istat = cudaMemcpyAsync(Q, QJ, 5*nx*ny*nz, cudaMemcpyDeviceToHost, io_stream)
```
Then immediately start the next `nt` inner steps on the default stream while the copy drains in `io_stream`. A sync before the next I/O cycle ensures the copy has landed. This hides D→H transfer latency behind compute.

**Files to modify:**
- `3D_solver/src/calc_time_dev.f90` — add `io_stream`, restructure outer `t2` loop to double-buffer
- `src/print.f90` — ensure host Q array is pinned-compatible

---

## Tier 3: Architectural (High Effort, ~1–2 weeks)

### 3A. Array layout transposition: `Q(5,nx,ny,nz)` → `Q(nx,ny,nz,5)`

**Finding:** Current AoS layout means consecutive threads loading component 1 of adjacent i-points must jump 40-byte strides (5 doubles apart). SoA layout `Q(nx,ny,nz,5)` makes all 5 component accesses stride-1 for x-direction kernels, maximizing cache line utilization during the shared-memory load phase.

**Scope:** Every file that accesses Q: all 6 kernel files + `calc_visc2.f90`, `calc_visc4.f90`, `calc_les.f90`, `calc_steps.f90`, `calc_physical_quantities.f90`, `calc_rescale.f90`, `set_bc_common.f90`, `preprocess.f90`, `print.f90`, `calc_time_dev.f90`. All `mod_globals.f90` array declarations.

**Strategy:** Add a parallel `QT(nx,ny,nz,5)` device array alongside `QJ(5,nx,ny,nz)`. Write a transpose kernel. Port kernels one at a time, validating against baseline on NSTGV TGV decay. Remove `QJ` only after all kernels are ported.

**Correctness validation:** Run NSTGV to t=1 (normalized) and compare kinetic energy spectrum against the published reference curves in `img/Ek.png`.

**Estimated impact:** 15–30% end-to-end speedup for compute-bound cases with large grids. Largest for NSTGV 513³ (all three direction kernels benefit equally). Smaller for SBLI where z-direction is thin (nz=33/257).

---

## Implementation Order

1. **1A** (contiguous in kernel args) — build and verify bitwise identical output
2. **1B** (persist sensor) — simple refactor, verify first 10-step output matches
3. **2A** (CUDA streams) — profile with `nsys` before/after; check for NaN fields as stream ordering bug indicator
4. **1C** (GPU-direct MPI) — only for SBLI/TBL; validate rescaling data with known reference state
5. **2B** (Roe/Hybrid internal) — compare flux arrays against full-kernel run for first 10 steps
6. **2C** (async I/O) — pin first, then add async memcpy
7. **3A** (layout transposition) — last; full regression suite required

## Critical Files

| File | Tier | Change |
|------|------|--------|
| `3D_solver/src/calc_keep_kernel.f90` | 1A | Add `contiguous` to 9 kernel subroutines |
| `3D_solver/src/calc_slau_kernel.f90` | 1A | Same + sensor arg |
| `3D_solver/src/calc_roe_kernel.f90` | 1A, 2B | Add `contiguous`; source for `_internal` variant |
| `3D_solver/src/calc_hybrid_kernel.f90` | 1A, 2B | Add `contiguous`; source for `_internal` variant |
| `3D_solver/src/calc_keep_kernel_internal.f90` | 1A | Add `contiguous` |
| `3D_solver/src/calc_slau_kernel_internal.f90` | 1A | Add `contiguous` |
| `3D_solver/src/calc_flux_base.f90` | 1B, 2A, 2B | Persist sensor; add streams/events; Roe/Hybrid `_internal` dispatch |
| `3D_solver/src/preprocess.f90` | 1B, 2A | Call `init_flux_base` for stream/sensor init |
| `3D_solver/src/calc_time_dev.f90` | 2C | Async I/O double-buffer |
| `3D_solver/SBLI/mod_globals.f90` | 1C | `id_gpumpi` kind=4 |
| `3D_solver/TBL/mod_globals.f90` | 1C | `id_gpumpi` kind=4 |
| `3D_solver/src/calc_roe_kernel_internal.f90` | 2B | **New file** |
| `3D_solver/src/calc_hybrid_kernel_internal.f90` | 2B | **New file** |
| `3D_solver/SBLI/Makefile` (and others) | 2B | Add `_internal.o` to OBJ |
