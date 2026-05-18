# Plan: Eliminate stack frame in calc_visc2.f90 by removing intent(out) scalar outputs from load_smem_visc2

## Root Cause

All 6 kernels in `3D_solver/src/calc_visc2.f90` call a device subroutine in
`3D_solver/src/load_smem_visc2.f90`. Each call has **4 `real(8), intent(out)` scalar arguments**:

| Load subroutine | intent(out) outputs |
|---|---|
| `load_smem_visc2_x` | `kTx, mux, mvx, mwx` |
| `load_smem_visc2_y` | `kTy, muy, mvy, mwy` |
| `load_smem_visc2_z` | `kTz, muz, mvz, mwz` |

In CUDA Fortran, Fortran pass-by-reference semantics require these scalars to live in **addressable
memory** (local memory = per-thread stack) in the caller. This creates an identical **32-byte stack
frame** (4 × 8 bytes) in all 6 kernels regardless of whether they are LES or not. LTO cannot inline
away the stack frame because `load_smem_visc2_*` uses `pipelineMemcpyAsync/pipelineCommit/
pipelineWaitPrior` from `use wmma` — CUDA built-ins the compiler cannot analyze for side effects,
blocking inlining. The subroutines are also absent from the Makefile's `-Minline` list.

**Why LES kernels use more registers but have the same stack frame**: The LES extra locals (sgs
variables) fit in registers; the 32-byte stack frame comes purely from the 4 by-reference scalars
shared by both LES and non-LES variants.

---

## Fix

Remove all 4 `intent(out)` scalars from each `load_smem_visc2_*` signature. Those subroutines
become pure shared-memory loaders: they fill `u/v/w` and call `syncthreads()`, nothing more.
Move the primary-direction stress/flux computation into the calling kernels after the call.

---

## Changes to `3D_solver/src/load_smem_visc2.f90`

Remove `use mod_constant, only : Cp_over_Pr` (no longer needed).

For each of the three subroutines, **remove** from the argument list: `dx`/`dy`/`dz`, `T`, `mu`,
and the 4 `intent(out)` scalars. Remove the corresponding `intent(in)` / `intent(out)` declarations
and the entire computation block after `pipelineCommit()` (the `mx`, `kTx/y/z`, `mu*/mv*/mw*`
calculations). Keep the pipeline load loop, `pipelineCommit`, `pipelineWaitPrior(0)`, and
`syncthreads()` unchanged.

**New signatures:**
```fortran
attributes(device) subroutine load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, Q, u, v, w)
attributes(device) subroutine load_smem_visc2_y(it, jt, kt, i, k, idx, nx, ny, nz, Q, u, v, w)
attributes(device) subroutine load_smem_visc2_z(it, jt, kt, i, j, idx, nx, ny, nz, Q, u, v, w)
```

**New body structure** (same for all three, shown for `_x`):
```fortran
! ... parameter sx/sy/sz, u/v/w intent(inout) declarations unchanged ...
integer i_base, ii, i, idx_l, offset_yz   ! (no mx, kTx, mux, mvx, mwx locals)
i_base    = (blockIdx%x-1)*blockDim%x
offset_yz = (jt-1)*sx + (kt-1)*sx*sy
do ii = it, threadsEv%x+1, blockDim%x
  i = i_base + ii
  idx_l = (ii-1) + offset_yz
  if (1 <= i .and. i <= nx .and. j <= ny .and. k <= nz) then
    call pipelineMemcpyAsync(u(idx_l), Q(i,2,j,k))
    call pipelineMemcpyAsync(v(idx_l), Q(i,3,j,k))
    call pipelineMemcpyAsync(w(idx_l), Q(i,4,j,k))
  else
    u(idx_l) = 0.d0; v(idx_l) = 0.d0; w(idx_l) = 0.d0
  endif
enddo
call pipelineCommit()
call pipelineWaitPrior(0)
call syncthreads()
```

---

## Changes to `3D_solver/src/calc_visc2.f90`

For each of the 6 subroutines, update the call site and add primary-direction computation inline.
The pattern is identical across all 6; the direction-specific variable names change.

### Pattern (shown for `calc_Ev2`, primary = x)

**Old call:**
```fortran
call load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, dx(i), Q, T, mu, u, v, w, viscous_work, mux, mvx, mwx)
if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
```

**New call + inlined computation:**
```fortran
call load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, Q, u, v, w)
if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
block
  real(8) mx
  mx           = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
  viscous_work = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
  mux          = mx * (-u(idx) + u(idx+1)) * dx(i)
  mvx          = mx * (-v(idx) + v(idx+1)) * dx(i)
  mwx          = mx * (-w(idx) + w(idx+1)) * dx(i)
end block
```

### Mapping for all 6 kernels

| Kernel | Load call | Primary metric | Outputs moved to caller |
|---|---|---|---|
| `calc_Ev2` / `calc_Ev_LES2` | `load_smem_visc2_x` | `dx(i)`, `mu(i,j,k)+mu(i+1,j,k)` | `viscous_work`, `mux`, `mvx`, `mwx` |
| `calc_Fv2` / `calc_Fv_LES2` | `load_smem_visc2_y` | `dy(j)`, `mu(i,j,k)+mu(i,j+1,k)` | `viscous_work`, `muy`, `mvy`, `mwy` |
| `calc_Gv2` / `calc_Gv_LES2` | `load_smem_visc2_z` | `dz(k)`, `mu(i,j,k)+mu(i,j,k+1)` | `viscous_work`, `muz`, `mvz`, `mwz` |

`Cp_over_Pr` is already imported in `calc_visc2.f90` via
`use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third` — no change needed.

---

## What does NOT change

- The async pipeline body of `load_smem_visc2_*` (loop, `pipelineMemcpyAsync`, `pipelineCommit`,
  `pipelineWaitPrior`, `syncthreads`) — unchanged
- All `u/v/w(idx)` shared memory accesses in `calc_visc2.f90`
- All cross-direction `mu/mut/Q` stencils (unchanged)
- The bounds check `if (nx-1 < i .or. ...) return` position (unchanged)

---

## Verification

```bash
cd 3D_solver/NSTGV
module load nvhpc-openmpi3
make clean && make 2>&1 | grep -E "Stack frame|calc_Ev2|calc_Fv2|calc_Gv2|calc_Ev_LES|calc_Fv_LES|calc_Gv_LES"
# Expect: stack frame = 0 bytes for all 6 kernels
bash calc.sh   # confirm no NaN and simulation completes
```
