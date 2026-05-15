# Plan: Double Buffering in calc_visc4_internal

## Context

The 4th-order viscous flux kernels (`calc_Ev4_in`, `calc_Fv4_in`, `calc_Gv4_in`) currently use a single set of three shared-memory buffers (`smem_a`, `smem_b`, `smem_c`) reused across 3 sequential stages. Each stage: loads data → waits (pipelineWaitPrior + syncthreads) → computes stress. There is no overlap between loading and computing.

The goal is to introduce **double buffering**: start loading the data for `calc_tau_cross` (stage 2) into a second smem buffer set *while* `calc_tau_straight` (stage 1) is being computed, hiding the load latency. Then overlap stage 3 load with stage 2 computation similarly.

## Approach

Two smem buffer pairs in each kernel:
- **Buffer 0** (`smem_a0`, `smem_b0`, `smem_c0`): stage 1 data (primary velocity + cross-derivatives for tau_straight), then reused for stage 3
- **Buffer 1** (`smem_a1`, `smem_b1`): stage 2 data (secondary velocity + cross-derivative for first tau_cross)

New execution order per kernel:
```
1. load_s1      → buffer0          (full: pipelineWaitPrior + syncthreads inside)
2. load_s2_issue → buffer1         (issue-only: async load + compute gradient, NO wait, NO sync)
3. calc_tau_straight  ← buffer0   (overlaps with s2 load in flight)
4. pipelineWaitPrior(0) + syncthreads   (wait for buffer1)
5. load_s3_issue → buffer0         (reuse buf0: async load + compute gradient, NO wait, NO sync)
6. calc_tau_cross     ← buffer1   (overlaps with s3 load in flight)
7. pipelineWaitPrior(0) + syncthreads   (wait for buffer0)
8. calc_tau_cross     ← buffer0   (s3 data now ready)
9. write flux output
```

## Files to Modify

### 1. `3D_solver/src/load_smem_visc4.f90`

Add 6 new `_issue` subroutines (no `pipelineWaitPrior`, no `syncthreads` at end).  
For **each direction** (x, y, z), add `s2_issue` and `s3_issue` variants identical to the existing `s2`/`s3` except the final two lines (`call pipelineWaitPrior(0)` + `call syncthreads()`) are removed.

Update the `public` list to export the 6 new names.

**x-direction non-TMA (`#else` branch):**
- `load_smem_visc4_x_s2_issue(it,jt,kt, j,k, nx,ny,nz, inv_dy, Q, vel, g1)` — async load v→vel + sync compute uy→g1, no wait/sync
- `load_smem_visc4_x_s3_issue(it,jt,kt, j,k, nx,ny,nz, inv_dz, Q, vel, g1)` — async load w→vel + sync compute uz→g1, no wait/sync

**x-direction TMA (`#if _CUDA_ARCH_ >= 900` branch):**
Same approach: `s2_issue` and `s3_issue` that omit the `barrier_arrive` / `barrier_try_wait_sleep` loop.  
The two barrier waits move into the kernel body after the corresponding compute step (see kernel changes below).

**y-direction** (no TMA branch):
- `load_smem_visc4_y_s2_issue` — async load u→vel + compute vx→g1, no wait/sync  
- `load_smem_visc4_y_s3_issue` — async load w→vel + compute vz→g1, no wait/sync

**z-direction** (no TMA branch):
- `load_smem_visc4_z_s2_issue` — async load u→vel + compute wx→g1, no wait/sync  
- `load_smem_visc4_z_s3_issue` — async load v→vel + compute wy→g1, no wait/sync

---

### 2. `3D_solver/src/calc_visc4_internal.f90`

Apply the double-buffer restructuring to all three kernels: `calc_Ev4_in`, `calc_Fv4_in`, `calc_Gv4_in`.

**Shared memory declarations** (example for `calc_Ev4_in`):
```fortran
! Buffer 0: stage 1 (tau_straight) and reused for stage 3
real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: smem_a0  !< vel: u, later w
real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: smem_b0  !< g1:  vy, later uz
real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: smem_c0  !< g2:  wz (stage 1 only)
! Buffer 1: stage 2 (first tau_cross)
real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: smem_a1  !< vel: v
real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: smem_b1  !< g1:  uy
```

Shared memory cost doubles for `smem_a`/`smem_b` (smem_c stays single):
- x-direction: 912 B → ~1520 B per block (within typical 96–192 KB SM limit)
- y/z-direction: ~11.6 KB per block (two 289-element double-precision arrays added)

**For TMA path (`#if _CUDA_ARCH_ >= 900`)**:  
Declare two shared barriers in the kernel: `integer(8), shared :: barrier_s2, barrier_s3`.  
After `load_smem_visc4_x_s2_issue` (which omits the barrier wait), the kernel does:
```fortran
token = barrier_arrive(barrier_s2)
do; if (barrier_try_wait_sleep(barrier_s2, token, 1000000) .ne. 0) exit; enddo
call syncthreads()
```
And similarly after `load_smem_visc4_x_s3_issue` with `barrier_s3`.

For the non-TMA path: simply `call pipelineWaitPrior(0)` + `call syncthreads()` in the kernel body after the compute step.

---

### 3. `3D_solver/src/calc_visc4_les_internal.f90`

Same double-buffer restructuring for `calc_Ev_LES4_in`, `calc_Fv_LES4_in`, `calc_Gv_LES4_in`.  
The only difference from the non-LES kernels: `calc_tau_straight_LES` / `calc_tau_cross_LES` instead of `calc_tau_straight` / `calc_tau_cross`, and additional `mut3(3)` and `Hsgs` variables.  
Buffer layout and pipeline structure are identical.

---

## Verification

```bash
cd 3D_solver/NSTGV
make clean && make        # must compile without errors
bash calc.sh              # run; check VTK output is produced and identical to pre-change reference
```

Check shared memory usage doesn't exceed GPU limit by inspecting ptxinfo output from `make` (the `-gpu=ptxinfo` flag is already in the Makefile).
