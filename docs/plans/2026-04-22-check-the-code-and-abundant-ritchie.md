# Code Audit: Bugs & Latent Bottlenecks — Final Status

All Round 1 and Round 2 issues have been fixed by the user.
**Pending:** Round 3 src-file optimizations below.

---

## All Fixed Issues

### Round 1 — Initial Bugs

| Bug | File | Fix |
|---|---|---|
| `calc_Fv2`/`calc_Gv2`: early return before `syncthreads()` → deadlock | `calc_visc2.f90` | Out-of-bounds threads write 0 to shared memory, then call `syncthreads()`, then return |
| `calc_Ducros`: Q array `dimension(5,nx,ny,nz)` vs caller's `(nx,5,ny,nz)` → wrong sensor | `calc_hybrid.f90` | Changed to `dimension(nx,5,ny,nz)`, all 9 accesses updated to `Q(i,comp,j,k)` |
| 2nd-order LES Hsgs missing from energy flux (all 3 directions) | `calc_visc2.f90` | `+ Hsgs` added to E, F, G energy writes |
| `coun` vs `count` typo in TMA path | `load_smem_visc4.f90:51` | `coun` → `count` |
| `CPU_MPI_ISEND`: local buf_cpu freed on subroutine return | `cpu_gpu_mpi.f90` | `MPI_WAIT` added before `deallocate` |
| `CPU_MPI_IRECV`: buf copied and freed before receive completes | `cpu_gpu_mpi.f90` | Changed to automatic array + `MPI_WAIT` before `buf = buf_cpu` |

### Round 2 — Additional Bugs

| Bug | File | Fix |
|---|---|---|
| `calc_hybrid_x2/y2/z2`: `fdx/fdy/fdz` never computed before use → arbitrary KEEP/SLAU | `calc_hybrid_kernel.f90:661,710,759` | Added sensor interpolation after bounds check in each kernel |
| `calc_keep_x2/y2/z2`: `device` attribute on local kernel arrays → non-standard semantics | `calc_keep_kernel.f90:401,428,455` | Removed `device` attribute; arrays are now plain locals (registers) |

### Round 2 — Code Clarity & Performance

| Issue | File | Fix |
|---|---|---|
| Hsgs formula: `0.125 * one_third` should be `one_24` for readability (same math) | `calc_visc4.f90:162,380,597` + `calc_visc4_les_internal.f90:71,136,201` | Rewritten as `1.125d0 * ... - one_24 * ...` matching the adjacent kTx formula |
| `cudaDeviceSynchronize()` after async copies on separate streams → stalls whole device | `calc_para.f90:233,247` | Changed to `cudaStreamSynchronize(1)` / `cudaStreamSynchronize(2)` |

### Configuration Changes (NSTGV test setup)

| Parameter | Old | New | Note |
|---|---|---|---|
| `id_scheme` | `real(2)` → SLAU | `integer(2)` → KEEP | Correct scheme for smooth TGV flow |
| `nx/ny/nz` | 513 | 66 | Reduced for fast testing |
| `id_bc_x/y/z` | `.false.` | `.true.` | Interior-only → full boundary kernels |
| `np/nt` | 100 / computed | 1 / 1 | Single-step smoke test |

---

## Round 3 — Src-File Optimizations

### OPT-1: `calc_step2_3` — hoist FP64 division out of inner loop ✓ FIXED

**File:** `3D_solver/src/calc_steps.f90:121`

Replaced `/ coef4` inside the `do l = 1, 5` loop with `* coef4_inv` where `coef4_inv = 1.d0/coef4` is computed once before the loop. Saves 4 FP64 divides per thread per stage-3 invocation.

---

### OPT-2: `calc_quantities_T_3D/2D` — replace `temp**1.5d0` with `temp * sqrt(temp)` ✓ FIXED

**Files:** `src/calc_physical_quantities.f90:55` (2D) and `:114` (3D)

Replaced `temp**1.5d0` (general power via `exp(log)`, ~100+ GPU cycles) with `temp * sqrt(temp)` (hardware sqrt, ~20 cycles). Runs in `!$cuf kernel do(3)` over all cells, 3× per time step.

---

### Analysis of remaining src files (no further optimizations found in active NS path)

Reviewed: `calc_keep_kernel.f90` + `calc_keep_3d.f90` (KEEP6 flux), `calc_visc4.f90` + `calc_visc_me4_base.f90` (tau subroutines), `calc_time_dev.f90` (RK loop), `calc_slau_kernel.f90` (SLAU), `calc_keep_kernel_internal.f90`. All divisions in hot paths are either precomputed reciprocal constants (one_24, one_twelfth, etc.), compile-time Fortran parameters (`Prt`, `dt`) that constant-fold under `-fast`, or already hoisted. No further issues in the active NS paths.

---

### Round 3 — LES Mode Issues (id_visc kind=8)

#### BUG-C: `calc_mut` — unconditional out-of-bounds array read before bounds check

**File:** `3D_solver/src/calc_les.f90:106-108`

**Problem:**
```fortran
i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
u3 = Q(i-1:i+1,2,j-1:j+1,k-1:k+1)   ! UNCONDITIONAL: i+1 can be > nx
v3 = Q(i-1:i+1,3,j-1:j+1,k-1:k+1)
w3 = Q(i-1:i+1,4,j-1:j+1,k-1:k+1)
if (3<=i .and. i<=nx-2 ...) then       ! bounds check comes AFTER reads
```
When `ceil((nx-2)/32)*32 > nx-2` (any nx not exactly divisible by 32 after trimming), some threads have `i > nx-1`, causing `Q(i+1,...)` = `Q(nx+1,...)` out-of-bounds global memory read. For the current nx=66 → max i=65 = nx-1, so it's safe today, but any other nx can trigger it.

**Fix:** Add `if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return` before line 106.

---

#### PERF-C: `SMS` function — multiple expensive power operations and runtime `tan()` call

**File:** `3D_solver/src/calc_les.f90:35,61,68,69`

**Problem:**
```fortran
real(8) :: pi = acos(-1.d0), Cm = 0.06d0, alpha = 0.5d0  ! local vars, not parameters
...
rtheta = tan2 / (tan(10.d0 * pi / 180.d0))**2  ! tan() called per thread
...
delta = cosh(...) * (dx * dy * dz)**(-1.d0/3.d0)  ! general power: cube root
nut   = ftheta * Cm * (S2**(0.5d0 * alpha))       ! S2**0.25 → sqrt(sqrt)
             * (qc2**(0.5d0 * (1.d0 - alpha)))     ! qc2**0.25 → sqrt(sqrt)
             * (delta**(1.d0 + alpha))              ! delta**1.5 → delta*sqrt(delta)
```
`pi`, `Cm`, `alpha` are local non-parameter variables → compiler can't constant-fold `tan(10*pi/180)`. Every thread computes an expensive `tan()` call. Four `x**real` power operations use the general ~100-cycle `exp(log)` path; 3 can be replaced with `sqrt`.

**Fix:**
1. Change `pi`, `Cm`, `alpha` to `parameter` constants (or define named constants at module level)
2. Precompute `tan_10deg_sq = (tan(10.d0 * acos(-1.d0) / 180.d0))**2` as a local parameter
3. Replace `S2**0.25d0` → `sqrt(sqrt(S2))`
4. Replace `qc2**0.25d0` → `sqrt(sqrt(qc2))`
5. Replace `delta**1.5d0` → `delta * sqrt(delta)` (same pattern as Sutherland fix)
6. Replace `(dx*dy*dz)**(-1.d0/3.d0)` → `exp(-log(dx*dy*dz) * one_third)` to avoid general power with non-simple exponent

Note: `SMS` is called once per grid cell in LES mode. With these changes, 4 expensive FP64 pow calls (~100 cycles each) and 1 trig call (~50 cycles) are replaced with hardware sqrt operations (~20 cycles each).

---

### OPT-3: `calc_flux_base.f90` — viscous dispatch redundancy (code clarity)

**Files:** `3D_solver/src/calc_flux_base.f90:249-271` (NS) and `:305-324` (LES)

**Assessment:** Not extractable into a shared helper without procedure pointers or optional arguments — both add complexity. All branches are compile-time constants. **No change recommended.**

---

### OPT-4: `CODE_ANALYSIS.md` — confirmed deleted

`CLAUDE.md` covers all necessary architectural context. `CODE_ANALYSIS.md` has been deleted and should stay deleted.

---

## Notes for Future Work

- `id_bc_x/y/z = .true.` activates full boundary kernels instead of the faster interior-only (`_in`) variants. For NSTGV with periodic BCs, `.false.` is correct and more efficient.
- `CPU_MPI_ISEND`/`CPU_MPI_IRECV` are now effectively blocking (MPI_WAIT inside), defeating their non-blocking intent. These routines have no call sites in the current solver (only blocking `CPUGPU_MPI_SEND/RECV` is used), so this is academic.
- The deleted `CODE_ANALYSIS.md` was a previous analysis document — its removal is tracked in git.
