# Code Quality & Performance Audit Plan

## Context

A comprehensive audit of the OUxSBLI GPU CFD solver to identify bugs, bottlenecks, and
redundancies. Focus areas: the declared bottleneck files (`calc_flux_base` call tree and
`calc_steps`), the viscous kernel chain (`calc_visc4`, shared-memory loaders), the
MPI halo exchange staging layer, and the SLAU convective flux kernels.
`cpu_gpu_mpi.f90` is excluded (under separate revision).

---

## Findings

### BUGS

---

#### BUG-1 · `calc_para.f90:178,214` — Large automatic (stack) arrays in halo exchange routines
**Severity: High**

```fortran
! exchange_cyclic (line 178) and exchange_rescale (line 214):
real(8), dimension(overlap*(ny-2)*(nz-6)*5) :: Qs_left, Qs_right, Qr_left, Qr_right
```

`Qs_left/Qs_right/Qr_left/Qr_right` are automatic host arrays whose size depends on
runtime parameters. With `overlap=3, ny=66, nz=66` this is `3×64×60×5 = 57600` doubles
= 460 KB per array, 1.8 MB total, allocated on the MPI host thread's stack on every call.
This risks stack overflow, especially on machines where thread stack size is limited.

**Fix:** Declare these as `allocatable`, `allocate` before use, and `deallocate` after.

---

#### BUG-2 · `calc_hybrid.f90:23,78` — `eps` declared as variable, not `parameter`
**Severity: Low–Medium**

```fortran
! calc_Ducros line 23, Albada line 78:
real(8) :: eps = 1.d-12
```

In CUDA `attributes(global/device)` code, a local variable with an initializer is **not**
a compile-time constant — the initializer is treated as assignment-on-entry, preventing
constant-folding and register re-use. The compiler cannot propagate `1.d-12` at the use
site (line 50).

**Fix:**
```fortran
real(8), parameter :: eps = 1.d-12
```

---

### BOTTLENECKS

---

#### BOTTLE-1 · `calc_hybrid_kernel.f90:47,156,265` — `rhor/ur/vr/wr/pr` in shared memory are per-thread private
**Severity: High**

```fortran
! calc_hybrid_x6 line 47 (and y6 line 156, z6 line 265):
real(8), dimension(threadsE%x, threadsE%y, threadsE%z), shared :: rhor, ur, vr, wr, pr
```

These 5 shared arrays are written and read exclusively at each thread's own index
`(it,jt,kt)`. No thread reads another thread's slot:
- `delta6(..., rhor(it,jt,kt))` / `rhor(it,jt,kt) = rho(it+1,jt,kt)` — written by thread `(it,jt,kt)` only
- `SLAU(..., rhor(it,jt,kt), ...)` — read by the same thread

Cross-thread shared memory access does not exist for these arrays. With typical
`threadsE = dim3(16,4,1)`, this wastes `5 × 16 × 4 × 1 × 8 = 2560 bytes` of shared
memory per block, reducing SM occupancy and increasing shared memory bank traffic.

**Fix:** Replace the 5 shared arrays with 5 local scalar variables in the kernel body:
```fortran
real(8) rho_r, u_r, v_r, w_r, p_r
```
Then replace `rhor(it,jt,kt)` → `rho_r`, `ur(it,jt,kt)` → `u_r`, etc. throughout the
`else` branch (delta6/delta4/2-pt) and the SLAU call. Apply identically to `y6` and `z6`.

---

#### BOTTLE-2 · `load_smem_visc4.f90:76-78` — Busy-wait loop on TMA barrier
**Severity: Medium–High (Hopper only)**

```fortran
! TMA path only (_CUDA_ARCH_ >= 900):
token = barrier_arrive(barrier)
do
  if (barrier_try_wait_sleep(barrier, token, 1000000) .ne. 0) exit
enddo
```

`barrier_try_wait_sleep` with timeout `1000000` ns (1 ms) retries in a loop until the TMA
barrier is satisfied. On Hopper, the mbarrier `.try_wait` instruction either completes
(barrier ready) or suspends the warp briefly. However, if the barrier is not signalled
within the timeout, the loop spins — silently hanging on TMA load failures (e.g., boundary
threads that skipped the `tma_bulk_load` call due to the `count > 0 .and. j <= ny ...`
guard at line 54). Threads that arrive at `barrier_arrive` but whose TMA was never
submitted (because `tid_linear != 1` or because the boundary guard fired) wait forever
for TMA completion that was never queued.

Additionally, the barrier arrival count at init (`blockDim%x * blockDim%y * blockDim%z`)
counts all thread arrivals, but `tma_bulk_load` internally registers expected bytes in the
same mbarrier — the arrival count must account for both thread arrivals and TMA byte
completions or the barrier may never reach its target.

**Fix:** Verify the mbarrier is initialized with
`count = nthreads + n_tma_transactions * bytes_per_transaction` as required by the
H100 PTX spec, or use the alternative `mbarrier.arrive.expect_tx` pattern where
`barrier_init` is called with count = nthreads and each `tma_bulk_load` registers its
byte count separately.

---

#### BOTTLE-3 · `calc_hybrid.f90:182` — Always-true `kind()` condition in `y6`/`z6` kernels
**Severity: Low**

```fortran
! calc_hybrid_y6 line 182 (and z6):
if (3 <= j .and. j <= ny-3 .and. 8 <= kind(id_accuracy)) then
```

`id_accuracy` is declared `integer(kind=8)` in `calc_hybrid_y6` — the overload is
dispatched at compile time for kind=8 callers. Therefore `kind(id_accuracy)` is always 8,
and `8 <= 8` is always `.true.`. The equivalent condition in `calc_hybrid_x6` (line 73)
simply uses `3 <= i .and. i <= nx-3` without the kind check. The dead condition generates
a redundant runtime comparison on every thread.

**Fix:** Remove `8 <= kind(id_accuracy)` from the y6 and z6 conditions (match x6 style).

---

#### BOTTLE-4 · `calc_visc4_les_internal.f90:51` / `calc_visc4.f90:48` — Redundant `syncthreads()` after `load_smem_visc4_x` on TMA path
**Severity: Low (Hopper only)**

```fortran
call load_smem_visc4_x(...)
call syncthreads()   ! ← redundant on TMA path
```

For the `#else` (non-TMA) path, `load_smem_visc4_x` contains no syncthreads, so the
caller's `syncthreads()` is required and correct.

For the TMA (`#if _CUDA_ARCH_ >= 900`) path, `load_smem_visc4_x` already executes:
- `syncthreads()` at line 47 (post-barrier-init)
- `syncthreads()` at line 74 (before barrier_arrive)
- full TMA barrier wait

The caller's subsequent `syncthreads()` is a third synchronization that is guaranteed to
be a no-op (all threads already synchronized). On Hopper, each `syncthreads()` has non-zero
overhead (bar.sync instruction). This affects all 6 callers: `calc_Ev4`, `calc_Ev_LES4`,
`calc_Fv4`, `calc_Fv_LES4`, `calc_Gv4`, `calc_Gv_LES4`, and all `_les_internal` variants.

**Fix:** Wrap the caller's `syncthreads()` in `#if _CUDA_ARCH_ < 900` to emit it only
for non-TMA builds, or move the syncthreads inside the `#else` branch of
`load_smem_visc4_x`.

---

#### BOTTLE-5 · `calc_steps.f90:30-32, 67-69, 109-111, 149-151` — Jacobian volume recomputed in all 4 RK kernels
**Severity: Low**

```fortran
! Verbatim in calc_step1, calc_step, calc_step2_3, calc_step4:
dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1))
dtdxdy = dt * 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1))
```

3 multiplications × 3 pairs × 4 kernels = 36 multiply-adds per cell per time step, all
computing the same values. For non-uniform grids, these values genuinely vary per cell and
cannot be precomputed globally, but the duplication increases code divergence risk.
Extracting to an `attributes(device)` helper ensures consistency across RK stages.

---

### REDUNDANCIES

---

#### REDUN-1 · `calc_physical_quantities.f90:7-116` — Four nearly identical routines
**Severity: Medium**

`calc_quantities_2D`, `calc_quantities_T_2D`, `calc_quantities_3D`, `calc_quantities_T_3D`
share identical kernels for `rho`, `u`, `v`, `p`. The `_T` variants only add:
```fortran
temp(i,j) = ...
mu(i,j)   = mu0_T0_S_over_T0_2_3 / (temp(i,j) + 111.d0)
```
110 lines of copy-pasted code. Extract viscosity computation into a separate
`calc_viscosity_3D` kernel called only when needed.

---

#### REDUN-2 · `calc_para.f90:14-168` — Seven nearly identical flatten/reconstruct kernels
**Severity: Low–Medium**

`flatten`, `flatten_left`, `flatten_right`, `flatten_rescale`, `reconstruct`,
`reconstruct_left`, `reconstruct_right`, `reconstruct_sbli_inlet` — all share the same
4-level loop nest and index formula `nj*ni*5*(k-1)+ni*5*(j-1)+ni*(l-1)+i`, differing only
in the source/destination `Q` index offset. A single parameterized kernel with an `x_offset`
argument would cover all cases with ~50 lines instead of ~170.

---

#### REDUN-3 · `calc_para.f90:193-205` — Unnecessary D2H/H2D round-trips in `exchange_cyclic`
**Severity: Medium**

```fortran
Qs_left  = Qs1d_left   ! device → host (blocking)
call MPI_SENDRECV(...)  ! host MPI
Qr1d_right = Qr_right  ! host → device (blocking)
```

The Fortran assignment between device and host arrays calls blocking `cudaMemcpy`.
Two copies per direction × 2 directions = 4 blocking D2H/H2D copies per halo exchange,
each separated by a blocking MPI call. This serializes all data movement and prevents
any overlap between transfers. The sends could be pipelined: start D2H for left, issue
async send, then D2H for right in parallel.

---

## SLAU Kernel Bugs and Optimizations

Audit of `3D_solver/src/calc_slau_kernel.f90` (603 lines, 9 kernels: x2/x4/x6, y2/y4/y6, z2/z4/z6)
and `3D_solver/src/calc_slau_3d.f90` (103 lines, `SLAU_common`, `SLAU1`, `HRSLAU2`).

---

### BUG-SLAU-A · `calc_slau_3d.f90:8` — `cr` never assigned in `SLAU_common`
**Severity: Critical**

```fortran
block
  real(8) cl, cr
  cl = sqrt(gamma * p1 * over_rho1)   ! line 7: left speed of sound ✓
  cl = sqrt(gamma * p2 * over_rho2)   ! line 8: BUG — should be `cr =`, overwrites cl
  c  = 0.5d0 * (cl + cr)              ! line 9: cr is undefined (→ 0 on GPU registers)
end block
```

`cr` is never set. GPU registers default to 0, so `c = 0.5 * sqrt(gamma*p2*over_rho2)`
instead of the correct arithmetic average `0.5*(cl+cr)`. This makes the interface sound
speed wrong by roughly 2×, corrupting Mach numbers `Mp`, `Mm` and all downstream
quantities: pressure splitting functions `bp`/`bm`, mass flux, and all 5 flux components.
This bug affects **both** `SLAU1` and `HRSLAU2`, and therefore all 9 SLAU kernels and
all SLAU calls within the Hybrid kernels.

**Fix:** Change line 8 from `cl = sqrt(gamma * p2 * over_rho2)` to `cr = sqrt(gamma * p2 * over_rho2)`.

---

### BUG-SLAU-B · `calc_slau_kernel.f90:416` — `jt` never assigned in `calc_slau_z4`
**Severity: Critical**

```fortran
it = threadIdx%x
it = threadIdx%y   ! line 416: BUG — should be `jt = threadIdx%y`
kt = threadIdx%z
```

`jt` is 0 (uninitialized register), `it` gets `threadIdx%y`. Then:
- `i = (blockIdx%x-1)*blockDim%x + it + 1` uses y-thread index instead of x → all x-index
  computations are wrong
- `j = (blockIdx%y-1)*blockDim%y + jt + 1` uses 0 → all threads compute `j = blockIdx%y*blockDim%y + 1`
- `offset_xy = (jt-1)*sz + (it-1)*sz*sy` → both dimensions wrong

Every cell written by `calc_slau_z4` is placed at the wrong (i, j) location.

**Fix:** Change line 416 from `it = threadIdx%y` to `jt = threadIdx%y`.

---

### BUG-SLAU-C · `calc_slau_kernel.f90:528` — `jt` never assigned in `calc_slau_y2`
**Severity: Critical**

```fortran
it = threadIdx%x
it = threadIdx%y   ! line 528: BUG — should be `jt = threadIdx%y`
kt = threadIdx%z
```

Same pattern as BUG-SLAU-B: `jt = 0`, `it = threadIdx%y`. All global j-indices and the
`offset_xz` computation are wrong. Every cell written by `calc_slau_y2` is misplaced.

**Fix:** Change line 528 from `it = threadIdx%y` to `jt = threadIdx%y`.

---

### BUG-SLAU-D · `calc_slau_kernel.f90:55-56` — `offset_yzr` never assigned in `calc_slau_x6`
**Severity: High**

```fortran
offset_yz  = (jt-1) * sx  + (kt-1) * sx  * sy   ! line 55
offset_yz  = (jt-1) * sx  + (kt-1) * sx  * sy   ! line 56: BUG — duplicate, offset_yzr is 0
```

`offset_yzr` stays 0, so `idx_r = it + 0 = it`. All threads with the same `it` but
different `jt` or `kt` read and write the same shared memory slot in `rhor/ur/vr/wr/pr`.
This is a race condition: two threads can write different values to `rhor(it)` and then
read back a corrupted result.

**Fix (preferred):** Applying PERF-SLAU-A (scalars) eliminates `idx_r` and `offset_yzr`
entirely. Alternatively change line 56 to
`offset_yzr = (jt-1) * sxr + (kt-1) * sxr * sy`.

---

### PERF-SLAU-A · `calc_slau_kernel.f90:46,122,200,274,343,412` — `rhor/ur/vr/wr/pr` shared arrays → scalars
**Severity: High**

```fortran
! In each of x4, x6, y4, y6, z4, z6 (e.g. x6 line 46):
real(8), dimension(sxr*sy*sz), shared :: rhor, ur, vr, wr, pr
```

These 5 arrays are written at `idx_r` and read back at `idx_r` by the **same thread** —
no cross-thread access occurs. Identical pattern to BOTTLE-1 already fixed in
`calc_hybrid_kernel.f90`. With `threadsE = dim3(16,4,1)`, this wastes 5×16×4×8 = 2560
bytes of shared memory per block per kernel.

Removing these arrays also resolves BUG-SLAU-D as a side effect (no more `idx_r`/`offset_yzr`).

**Fix:** In each affected kernel, replace:
```fortran
! REMOVE these 5 shared declarations and idx_r / offset_yzr variables
real(8), dimension(sxr*sy*sz), shared :: rhor, ur, vr, wr, pr
integer idx_r, offset_yzr            ! (naming varies per direction)

! ADD 5 scalar locals
real(8) rho_r, u_r, v_r, w_r, p_r
```
Replace every `rhor(idx_r)` → `rho_r`, `ur(idx_r)` → `u_r`, etc. throughout delta6/delta4/2-pt
branches and the `SLAU(...)` call. Apply to all 6 kernels: x4, x6, y4, y6, z4, z6.
The 2nd-order kernels (x2, y2, z2) read directly from the main shared tile and have no
`_r` arrays — they need no change.

---

### CLEAN-SLAU-A · `calc_slau_kernel.f90:196-197` — Duplicate parameter declarations in `calc_slau_z6`
**Severity: Low**

```fortran
integer, parameter :: sy  = threadsG%y       ! line 194
integer, parameter :: sz  = threadsG%z + 5   ! line 195
integer, parameter :: sy  = threadsG%y       ! line 196: DUPLICATE
integer, parameter :: sz  = threadsG%z + 5   ! line 197: DUPLICATE
```

Lines 196–197 are exact copies of lines 194–195. Remove them.

---

## Status of All Items

| Item | File | Status |
|------|------|--------|
| BUG-1 (allocatable halo arrays) | `calc_para.f90:178` | ✅ Done |
| BUG-2 (eps as parameter) | `calc_hybrid.f90:23,79` | ✅ Done |
| BOTTLE-1 (shared arrays → scalars, hybrid) | `calc_hybrid_kernel.f90:47,154,263` | ✅ Done |
| BOTTLE-2 (TMA mbarrier count) | `load_smem_visc4.f90:76-78` | ⏸ Deferred |
| BOTTLE-3 (dead kind() condition) | `calc_hybrid_kernel.f90:180` | ✅ Done |
| BOTTLE-4 (redundant syncthreads) | `calc_visc4.f90`, `calc_visc4_les_internal.f90` | ✅ N/A — no syncthreads present |
| BOTTLE-5 (Jacobian helper) | `calc_steps.f90` | ✅ Done |
| REDUN-1 (calc_physical_quantities) | `src/calc_physical_quantities.f90` | ⏸ Skipped |
| REDUN-2 (flatten/reconstruct) | `calc_para.f90:14-168` | ⏸ Skipped |
| REDUN-3 (exchange_cyclic pipeline) | `calc_para.f90:193-205` | ✅ Done |
| BUG-SLAU-A (cr never assigned) | `calc_slau_3d.f90:8` | ✅ Done |
| BUG-SLAU-B (jt uninitialized, z4) | `calc_slau_kernel.f90:416` | ✅ Done |
| BUG-SLAU-C (jt uninitialized, y2) | `calc_slau_kernel.f90:528` | ✅ Done |
| BUG-SLAU-D (offset_yzr uninitialized, x6) | `calc_slau_kernel.f90:55-56` | ✅ Done |
| PERF-SLAU-A (shared arrays → scalars, SLAU) | `calc_slau_kernel.f90` | ❌ Not done |
| CLEAN-SLAU-A (duplicate parameters, z6) | `calc_slau_kernel.f90:196-197` | ✅ Done |

Note: `calc_slau_kernel_internal.f90` (14 KB) likely has the same patterns and should
be audited separately.

---

## Recommended Fix Order (SLAU)

1. **BUG-SLAU-A** (`calc_slau_3d.f90:8`) — one-character fix; restores physically correct
   sound speed in ALL SLAU and Hybrid-SLAU fluxes. Fix first before any profiling.
2. **BUG-SLAU-B** (`calc_slau_kernel.f90:416`) — one-character fix; restores correct
   z-flux computation in 4th-order SLAU.
3. **BUG-SLAU-C** (`calc_slau_kernel.f90:528`) — one-character fix; restores correct
   y-flux computation in 2nd-order SLAU.
4. **PERF-SLAU-A** (`calc_slau_kernel.f90`) — 5 shared arrays → scalars in 6 kernels;
   also resolves BUG-SLAU-D as a side effect. Apply to x4/x6, y4/y6, z4/z6.
5. **CLEAN-SLAU-A** (`calc_slau_kernel.f90:196-197`) — remove 2 duplicate lines.

---

## Verification

After each fix, verify by:

```bash
cd 3D_solver/NSTGV
make clean && make          # must compile without new errors/warnings
bash calc.sh                # run to completion; confirm VTK output produced
# Numeric regression: compare Q norms against a known-good baseline
mpirun -n 2 ./a.out > out_new.txt
diff <(grep "step" out_baseline.txt) <(grep "step" out_new.txt)
# For BOTTLE-1 (shared memory): profile with nsys to verify shared memory reduction
nsys profile --stats=true ./a.out
# For BOTTLE-2 (TMA barrier): run on Hopper and check for hangs or cudaErrorUnknown
```
