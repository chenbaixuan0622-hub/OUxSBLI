# Plan: Fix regression + H100 optimizations for `load_smem_visc4`

## Context

The `pipeline` branch introduced 3-stage `pipelineMemcpyAsync`-based async loading in `load_smem_visc4.f90` to overlap velocity data loading with gradient computation. NCU profiling of the NSTGV case shows the `calc_visc4_internal_calc_gv4_in_` kernel regressed +37.6% (9.77 ms → 13.47 ms).

**Root cause**: The GPU (SM frequency 1.52–1.53 GHz, DRAM 2.62 GHz — consistent with A100 at reduced clock, or V100) does not benefit from `pipelineMemcpyAsync` / `cp.async` for this kernel. The staged pipeline adds net overhead instead of benefit:

| Metric | Stable | Pipeline | Δ |
|--------|--------|----------|---|
| Duration (calc_gv4_in) | 9.77 ms | 13.47 ms | +37.6% |
| Registers/thread | 60 | 72 | +20% |
| Achieved occupancy | 47.98% | 42.73% | −10.9% |
| DRAM throughput | 57.4% | 41.7% | −27.4% |
| Memory throughput | 71.3% | 58.2% | −18.4% |

**Why pipeline loses**: `cp.async` (`pipelineMemcpyAsync`) was supposed to bypass the register file and free up registers. Instead, the 3-stage structure forces the compiler to keep more live registers (boundary-path intermediates `mx, muy, mvy, muz, mwz`, and stage-persistent `mu3(3)`) alive across stage boundaries. This increases register count 60→72, reduces occupancy, fewer active warps → less latency hiding → lower DRAM throughput. The kernel is memory-latency-bound, so occupancy is the primary driver of performance.

## Files to Modify

- `3D_solver/src/load_smem_visc4.f90` — the only file that needs changes

## Implementation Plan

### Change: Replace `pipelineMemcpyAsync` with direct writes in all non-Hopper paths

The x-direction subroutines already have a `#if _CUDA_ARCH_ >= 900` guard that uses TMA bulk loads on Hopper. The `#else` block (and all y/z direction subroutines which have no guard) use `pipelineMemcpyAsync`. Replace the element-by-element async copy pattern with a simple direct write, keeping a single `syncthreads()` at the end.

**Pattern to replace** (blocking variants `_s1`, `_s2`, `_s3`):
```fortran
! BEFORE
do ii = ..
  call pipelineMemcpyAsync(vel(idx), Q(i,comp,j,k))
enddo
call pipelineCommit()
do ii = ..
  g1(idx) = [compute from Q]
enddo
call pipelineWaitPrior(0)
call syncthreads()

! AFTER
do ii = ..
  vel(idx) = Q(i,comp,j,k)   ! direct write
enddo
do ii = ..
  g1(idx) = [compute from Q]
enddo
call syncthreads()
```

**Pattern to replace** (non-blocking `_issue` variants `_s2_issue`, `_s3_issue`):
```fortran
! BEFORE
do ii = ..
  call pipelineMemcpyAsync(vel(idx), Q(i,comp,j,k))
enddo
call pipelineCommit()
do ii = ..
  g1(idx) = [compute from Q]
enddo
! (no pipelineWaitPrior/syncthreads — caller provides it)

! AFTER
do ii = ..
  vel(idx) = Q(i,comp,j,k)   ! direct write
enddo
do ii = ..
  g1(idx) = [compute from Q]
enddo
! (no syncthreads — caller provides it, same as before)
```

### Affected subroutines (all in the `#else` block or unconditional):

**X-direction `#else` block** (lines ~234–407):
- `load_smem_visc4_x_s1` — blocking: remove pipelineMemcpyAsync/Commit; change pipelineWaitPrior+sync → single sync
- `load_smem_visc4_x_s2` — blocking: same
- `load_smem_visc4_x_s3` — blocking: same
- `load_smem_visc4_x_s2_issue` — non-blocking: remove pipelineMemcpyAsync/Commit; no sync at end
- `load_smem_visc4_x_s3_issue` — non-blocking: same

**Y-direction** (lines ~413–582, unconditional):
- `load_smem_visc4_y_s1` — blocking: same pattern
- `load_smem_visc4_y_s2` — blocking: same
- `load_smem_visc4_y_s3` — blocking: same
- `load_smem_visc4_y_s2_issue` — non-blocking: same
- `load_smem_visc4_y_s3_issue` — non-blocking: same

**Z-direction** (lines ~589–757, unconditional):
- `load_smem_visc4_z_s1` — blocking: same pattern
- `load_smem_visc4_z_s2` — blocking: same
- `load_smem_visc4_z_s3` — blocking: same
- `load_smem_visc4_z_s2_issue` — non-blocking: same
- `load_smem_visc4_z_s3_issue` — non-blocking: same

### Note on `calc_visc4_internal.f90`

The `pipelineWaitPrior(0)` calls in `calc_visc4_internal.f90` (used after `_issue` calls in non-Hopper path) become no-ops since there are no more pending async copies. Leave them in place for now — they are harmless and the subsequent `syncthreads()` still correctly synchronizes the direct writes. (Removing them is a minor cleanup for a later commit.)

## Expected Outcome

- Register count: 72 → ~60 (restored, since no pipeline state tracking)
- Occupancy: 42.73% → ~48% (restored)
- DRAM throughput: 41.7% → ~57% (restored via occupancy recovery)
- Duration: 13.47 ms → ~9.8 ms (target: back to stable level)
- `calc_visc4.f90` boundary kernel (`calc_Ev4`, `calc_Fv4`, `calc_Gv4`) also benefits from same load_smem fix

## Current state (after user's replacement)

The user has already replaced the files with a cleaner design:
- **New API**: single `load_smem_visc4_x/y/z` call per direction, loading all 7 smem arrays at once (u, v, w + 4 gradients)
- `calc_visc4.f90` and `calc_visc4_internal.f90` both updated to use the new 7-array interface
- **H100 path** (`_CUDA_ARCH_ >= 900`): TMA bulk loads for u, v, w (x-direction only); gradients computed in parallel loops
- **Non-H100 path**: `pipelineMemcpyAsync` for u, v, w (x-direction); overlap with gradient computation. y/z directions use `pipelineMemcpyAsync` unconditionally

---

## H100 Optimization — Immediate Change

### A. Remove redundant `syncthreads()` before `barrier_arrive` in `load_smem_visc4_x`

**File**: `3D_solver/src/load_smem_visc4.f90`, inside the `#if _CUDA_ARCH_ >= 900` block of `load_smem_visc4_x`, around line 75.

**What to change**: Delete the `call syncthreads()` that appears immediately before `token = barrier_arrive(barrier)`.

**Why it is safe to remove**: H100's `mbarrier` provides `st.release` semantics at `barrier_arrive` and `ld.acquire` semantics at `barrier_try_wait_sleep`. Each thread reaches `barrier_arrive` only after completing its gradient loop body — sequential per-thread execution guarantees all shared memory writes precede the arrive. When `barrier_try_wait_sleep` returns, it guarantees both: all threads have arrived (i.e., completed their writes) AND the TMA transfer has finished. The explicit `syncthreads()` duplicates this guarantee and wastes a warp-level synchronization barrier.

Current code (lines ~74–80):
```fortran
    call syncthreads()           ! ← REMOVE this line
    token = barrier_arrive(barrier)
    do
      if (barrier_try_wait_sleep(barrier, token, 1000000) .ne. 0) exit
    enddo
```

After change:
```fortran
    token = barrier_arrive(barrier)
    do
      if (barrier_try_wait_sleep(barrier, token, 1000000) .ne. 0) exit
    enddo
```

The `syncthreads()` at line 48 (after `barrier_init`) must stay — it ensures the initialized barrier is visible to all threads before any TMA call proceeds.

---

### Future H100 Optimizations (not implemented now)

### Priority 2 (Medium effort — significant speedup on H100)

#### C. Add TMA bulk loads for y-direction velocity in H100 path

**File**: `load_smem_visc4.f90` — add `#if _CUDA_ARCH_ >= 900` block to `load_smem_visc4_y`

The x-direction TMA works because Q(i, comp, j, k) is contiguous along i (fastest Fortran index). For y-direction, Q(i, comp, j, k) with varying j is NOT contiguous (stride = nx*5 doubles). A 1D TMA bulk load cannot be used directly.

**Solution**: Use the H100 TMA 2D descriptor (`CUtensorMap`) to describe the strided tensor layout. This requires:
1. Adding a `CUtensorMap` parameter to the kernel (created at host level in `cpu_gpu_mpi.f90` or similar)
2. In `load_smem_visc4_y` H100 path: replacing `pipelineMemcpyAsync` loop with `tma_bulk_load` using the descriptor for the y-direction slice

This eliminates per-element cp.async instructions (replaced by a single bulk DMA), freeing registers and reducing instruction count.

**Trade-off**: Host-side setup complexity. The tensor descriptor must be created once and passed to each kernel invocation.

#### D. Add TMA bulk loads for z-direction velocity in H100 path

Same approach as C, but for z-direction. Stride = nx*5*ny doubles between z-elements.

---

### Priority 3 (Architecture-level)

#### E. Increase thread block size for viscous kernels on H100

Current: `threadsEv = dim3(32,1,1)`, `threadsFv = dim3(32,4,1)`, `threadsGv = dim3(32,1,4)` → 32, 128, 128 threads/block.

For H100 with 228 KB shared memory per SM and 7 arrays × 9 elements × 8 bytes ≈ 16 KB/block:
- Current: ~14 blocks × 128 threads = 1792 threads/SM (87.5% of peak 2048)
- With `dim3(32,2,4)` = 256 threads/block and ~32 KB smem: 7 blocks × 256 = 1792 threads (same occupancy but fewer, larger blocks → better ILP)
- With `dim3(32,1,8)` = 256 threads/block and sz = 8+5 = 13 → 32×1×13×8×7 = 23 KB/block: ~9 blocks × 256 = 2304 → capped at 2048 (100% occupancy)

**Recommendation**: Try `threadsGv = dim3(32,1,8)` and `threadsFv = dim3(32,8,1)` for production grid sizes. Benchmark with ncu before committing. Larger blocks improve TMA amortization (fewer barrier_init calls, more data per TMA transfer).

---

## Verification

1. `cd 3D_solver/NSTGV && make clean && make` — confirm compilation success
2. `bash test_cicd.sh` from repo root — verify all 12 build targets still pass
3. Re-run `bash profile.sh` in NSTGV and compare with `ncu/my_report_ncu_20260519_pipeline.json`
4. Check that `calc_visc4_internal_calc_gv4_in_` duration is back to ~9.8 ms (regression fixed) and ideally below (H100 gains)
5. `pytest ouxsbli/tests/` — confirm physical correctness is unchanged
