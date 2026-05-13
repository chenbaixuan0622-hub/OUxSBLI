# Fix: NaN after refactoring shared-memory loading to device subroutines

## Context

The y/z shared-memory loading loops were extracted into `load_smem_visc4_y` and `load_smem_visc4_z` (mirroring the already-working `load_smem_visc4_x`). The code compiled and ran, but produced NaN results. The x-direction worked fine; y and z did not. All implementation files (`load_smem_visc4.f90`, `calc_visc4.f90`, `calc_visc4_internal.f90`, `calc_visc4_les_internal.f90`) and the Makefile dependency lines have already been updated correctly.

---

## Root cause of NaN

The subroutine refactoring is logically identical to the inline code — same conditions, same formulas, same array indices, same `syncthreads()` call placement. The NaN is not a logic error.

**Key discriminating factor — warp count per block:**

| Direction | Thread block | Warps | x worked? |
|-----------|-------------|-------|-----------|
| x (Ev) | `threadsEv = dim3(32,1,1)` | **1 warp** | Yes |
| y (Fv) | `threadsFv = dim3(32,4,1)` | **4 warps** | No (NaN) |
| z (Gv) | `threadsGv = dim3(32,1,4)` | **4 warps** | No (NaN) |

With a single warp, all 32 threads always execute in lock-step — `syncthreads()` is trivially satisfied even without a proper barrier. With 4 warps (128 threads), the barrier is essential: warp A may start reading shared memory while warp B hasn't finished writing yet, producing uninitialized values → NaN.

**The actual compiler issue:** When `syncthreads()` appears inside a **non-inlined** `attributes(device)` subroutine (called from a `attributes(global)` kernel), the NVIDIA HPC SDK compiler does not reliably emit a true `__syncthreads()` PTX barrier across all warps. The x-direction happened to work because 1 warp masked the problem entirely.

---

## Fix — single line change in `3D_solver/NSTGV/Makefile`

Force-inline the three loading subroutines by adding them to the `-Minline` list.

**Current (line 3):**
```makefile
subroutines = name:calc_tau_straight,calc_tau_cross,flux4
```

**Change to:**
```makefile
subroutines = name:calc_tau_straight,calc_tau_cross,flux4,load_smem_visc4_x,load_smem_visc4_y,load_smem_visc4_z
```

This forces the compiler to inline the subroutine call sites into the calling `attributes(global)` kernels, so `syncthreads()` is generated inline as a proper block-level barrier — identical to what the original explicit loops produced.

No other files need to change.

---

## Verification

```bash
cd 3D_solver/NSTGV
make clean && make 2>&1 | grep -iE "inline|error"
```
Look for `load_smem_visc4_x`, `load_smem_visc4_y`, `load_smem_visc4_z` in the `-Minfo=inline` output — each should show "inlined into" the corresponding kernel. Then run the simulation and confirm results are no longer NaN.
