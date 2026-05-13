# Optimization Plan: OUxSBLI GPU-Accelerated CFD Solver (Phase 3)

## Context

All convective and time-integration optimizations from Phase 1/2 are complete:

| Commit | Items |
|--------|-------|
| 4c65c86 | 1A: `contiguous` in all convective kernel args; CLAUDE.md added |
| 004a263 | Smem arg cleanup; DeviceSync removal; 2A (streams) denied |
| 39bc8af | 1D: `contiguous` in `calc_steps.f90`; 1E: `E_next` shuffle removed |
| 95adb01 | 1B: persistent `sensor` device array; 2B: Roe `_internal`; 2B-hybrid: Hybrid `_internal` |

**What remains:** The viscous flux kernels (`calc_visc2.f90`, `calc_visc4.f90`, `calc_les.f90`) have NOT been optimized. They are in the hot path via `calc_EFG_visc` / `calc_EFG_LES`, called at every RK stage. Two distinct bugs/inefficiencies exist.

---

## Status Summary

| Item | Status |
|------|--------|
| 1A–1E, 1B, 2B, 2B-hybrid | ✅ Done |
| 2A — CUDA streams | ❌ Denied (occupancy) |
| **3A — `contiguous` in viscous kernels** | 🆕 New |
| **3B — Remove `device` from local vars in viscous kernels** | 🆕 New (bug fix) |
| **3C — `_internal` variants for visc4** | 🆕 New |
| 1C — GPU-direct MPI | ⏳ Low priority |
| 2C — Async I/O | ⏳ Low priority |
| 4A — Array layout transposition | ⏳ Architectural |

---

## Implementation Plan (in order)

### Step 1: 3A — Add `contiguous` to viscous kernel dummy arguments

**Files:** `3D_solver/src/calc_visc2.f90`, `3D_solver/src/calc_visc4.f90`, `3D_solver/src/calc_les.f90`

All `attributes(global)` kernels in these files declare device array dummy arguments **without** `contiguous`. The compiler cannot assume contiguity and generates conservative gather/scatter loads. This is the same issue as 1A for convective kernels.

**`calc_visc2.f90` — 6 kernels:** `calc_Ev2`, `calc_Fv2`, `calc_Gv2`, `calc_Ev_LES2`, `calc_Fv_LES2`, `calc_Gv_LES2`

For each, change:
```fortran
real(8), intent(in), device    :: dx(nx-1)
real(8), intent(in), device    :: Q(5,nx,ny,nz)
real(8), intent(in), device    :: T(nx,ny,nz)
real(8), intent(in), device    :: mu(nx,ny,nz)
real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
```
to:
```fortran
real(8), intent(in), device, contiguous    :: dx(nx-1)
real(8), intent(in), device, contiguous    :: Q(5,nx,ny,nz)
real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
```
Apply to all array arguments (`dx`, `dy`, `dz`, `Q`, `T`, `mu`, `mut`, `qc2`, `E`/`F`/`G`).

**`calc_visc4.f90` — 6 kernels:** `calc_Ev4`, `calc_Fv4`, `calc_Gv4`, `calc_Ev_LES4`, `calc_Fv_LES4`, `calc_Gv_LES4`. Same change pattern.

**`calc_les.f90` — 1 kernel:** `calc_mut`. Confirm and add `contiguous` to its device array args.

---

### Step 2: 3B — Remove `device` attribute from local variables in viscous kernels

**This is a bug.** In CUDA Fortran, declaring a local variable with `device` attribute inside an `attributes(global)` subroutine allocates it in device global memory per thread — not in registers. This is far slower than register-resident scalars.

**`calc_visc2.f90` — LES kernels (e.g., `calc_Ev_LES2`, line 118):**
```fortran
! Before (device global memory — slow):
real(8), dimension(2), device :: my, mysgs, mz, mzsgs

! After (register/local memory — fast):
real(8), dimension(2) :: my, mysgs, mz, mzsgs
```
Apply the same fix in `calc_Fv_LES2` and `calc_Gv_LES2`.

**`calc_visc4.f90` — all 6 kernels.** Two patterns appear:
- Interior block (e.g., line 174): `real(8), device :: mu3(3), mut3(3)` → remove `device`
- Boundary fallback block (e.g., lines 193, 203): `real(8), device :: my(2)`, `real(8), device :: mz(2)` → remove `device`

Example fix for `calc_Ev4` interior block:
```fortran
! Before:
real(8), device :: mu3(3)
real(8), device :: kTx3(3)

! After:
real(8) :: mu3(3)
real(8) :: kTx3(3)
```
Example fix for `calc_Ev4` boundary block:
```fortran
! Before:
real(8), device :: my(2)
real(8), device :: mz(2)

! After:
real(8) :: my(2)
real(8) :: mz(2)
```
Systematically apply to ALL occurrences in all 6 visc4 kernels (x, y, z directions × non-LES + LES = 6).

---

### Step 3: 3C — Create `_internal` variants for visc4

**Finding:** For periodic domains (ETGV, NSTGV, KHI: `id_bc_x/y/z = .false.`), all threads in `calc_Ev4` etc. only enter the interior 4th-order branch (`if (3 <= i .and. i <= nx-3 ...)`) — the `if (id_bc_x)` boundary block never executes. Creating `_internal` variants eliminates the dead branch and its register pressure from the inactive code path.

**New file:** `3D_solver/src/calc_visc4_internal.f90`

Structure:
```fortran
module calc_visc4_internal
  use mod_globals, only : ...
  use mod_constant, only : ...
  implicit none
  private
  public calc_Ev4_in, calc_Fv4_in, calc_Gv4_in
  ! LES variants if needed: calc_Ev_LES4_in, calc_Fv_LES4_in, calc_Gv_LES4_in
contains
  attributes(global) subroutine calc_Ev4_in(nx, ny, nz, dx, dy, dz, Q, T, mu, E)
    ! Same args as calc_Ev4, all with contiguous
    ! Body: shared memory load + interior stencil only; NO if(id_bc_x) block
  end subroutine
  ! ... Fv4_in, Gv4_in similarly
end module calc_visc4_internal
```

**`3D_solver/src/calc_flux_base.f90`** — add `use calc_visc4_internal` and dispatch in `calc_EFG_visc`:
```fortran
! Before:
call calc_Ev4<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
call calc_Fv4<<<blocksFv,threadsFv>>>(...)
call calc_Gv4<<<blocksGv,threadsGv>>>(...)

! After:
if (id_bc_x) then
  call calc_Ev4<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
else
  call calc_Ev4_in<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
endif
! (repeat for Fv4, Gv4)
```
Apply same dispatch for LES variants in `calc_EFG_LES` if LES with 4th-order viscous is used.

**Makefiles** — add `calc_visc4_internal.o` to OBJ in all 6 cases, with dependency:
```makefile
calc_visc4_internal.o: mod_globals.mod mod_constant.mod
```

---

## Critical Files

| File | Change |
|------|--------|
| `3D_solver/src/calc_visc2.f90` | Step 1: `contiguous` on all dummy args; Step 2: remove `device` from LES local vars |
| `3D_solver/src/calc_visc4.f90` | Step 1: `contiguous` on all dummy args; Step 2: remove `device` from local block vars |
| `3D_solver/src/calc_les.f90` | Step 1: `contiguous` on `calc_mut` dummy args |
| `3D_solver/src/calc_visc4_internal.f90` | Step 3: new file — interior-only visc4 kernels |
| `3D_solver/src/calc_flux_base.f90` | Step 3: `use calc_visc4_internal`; id_bc dispatch for visc4 |
| `3D_solver/*/Makefile` (all 6) | Step 3: add `calc_visc4_internal.o` |

---

## Verification

Build and run NSTGV (periodic, NS, likely visc4) and SBLI (wall-bounded, NS):

```bash
cd 3D_solver/NSTGV
make clean && make
bash calc.sh

cd 3D_solver/SBLI
make clean && make
bash calc.sh
```

- Output values (kinetic energy, entropy) should match baseline within round-off.
- Profile with `bash profile.sh` (nsys) before/after to quantify improvement.
- No NaN/Inf in output VTK files.

---

## Deferred Items

- **1C** (GPU-direct MPI) — SBLI/TBL only; requires CUDA-aware MPI verification
- **2C** (Async I/O) — pin host Q, use `cudaMemcpyAsync`
- **4A** (Array layout transposition) — high effort, full regression suite required
