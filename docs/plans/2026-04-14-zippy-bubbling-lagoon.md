# Optimization Plan: OUxSBLI GPU-Accelerated CFD Solver (Phase 4)

## Context

Continuing GPU optimization of OUxSBLI. All prior convective and viscous kernel items from Phases 1–3 have now been closed or are pending below. The uncommitted change in `3D_solver/ETGV/mod_globals.f90` switches the ETGV case to Roe scheme (`real(4)`) at 6th-order accuracy (`integer(kind=8)`) with no TVD (`integer(kind=2)`) — this is a configuration experiment, not an optimization task.

---

## Completed Items (do not redo)

| Commit | Item | Description |
|--------|------|-------------|
| 4c65c86 | 1A | `contiguous` in all convective kernel args |
| 004a263 | — | Unnecessary smem arg removed; DeviceSync cleanup |
| 39bc8af | 1D/1E | `contiguous` in `calc_steps.f90`; warp shuffle removed |
| 95adb01 | 1B / 2B / 2B-hybrid | Persistent `sensor`; Roe `_internal`; Hybrid `_internal` |
| 1ea54b8 | **3A** | `contiguous` on all dummy args in `calc_visc2.f90`, `calc_visc4.f90`, `calc_les.f90` |
| 1ea54b8 | **3B** | `device` removed from local block variables in visc2/visc4 LES kernels (bug fix) |

2A (CUDA streams) was explicitly denied due to occupancy concerns.

---

## Current Status

| Item | Status |
|------|--------|
| 3A — `contiguous` in viscous kernels | ✅ Done |
| 3B — Remove `device` from local vars | ✅ Done |
| **3C — `_internal` variants for visc4** | 🔲 Next |
| 1C — GPU-direct MPI (SBLI/TBL) | ⏳ Low priority |
| 2C — Async I/O | ⏳ Low priority |
| 4A — Array layout transposition | ⏳ Architectural |

---

## Next Step: 3C — `_internal` Variants for Visc4

### Motivation

For periodic cases (ETGV, NSTGV, KHI: `id_bc_x/y/z = .false.`), all threads in `calc_Ev4` etc. only enter the interior 4th-order branch — the `if (id_bc_x)` boundary block never executes. Creating `_internal` variants eliminates the dead boundary-branching code and its register pressure, following the established `_internal` pattern already used for all convective kernels.

Current `calc_flux_base.f90` always calls full kernels (lines 249–251) with no `id_bc` dispatch for viscous.

### New File: `3D_solver/src/calc_visc4_internal.f90`

Mirror the structure of `calc_slau_kernel_internal.f90`. Expose three public procedures:

```fortran
module calc_visc4_internal
  use mod_globals, only : ...
  use mod_constant, only : ...
  implicit none
  private
  public :: calc_Ev4_in, calc_Fv4_in, calc_Gv4_in
contains
  attributes(global) subroutine calc_Ev4_in(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
    ! Same args as calc_Ev4 in calc_visc4.f90, all with contiguous
    ! Body: interior 4th-order stencil only — NO if(id_bc_x) boundary block
  end subroutine
  ! ... Fv4_in, Gv4_in similarly
end module
```

LES variants (`calc_Ev_LES4_in` etc.) can be added if LES+4th-order viscous cases are used; defer if not immediately needed.

### Changes to `3D_solver/src/calc_flux_base.f90`

Add `use calc_visc4_internal` to the module USE list.

In `calc_EFG_visc` (around line 249), replace the unconditional calls with `id_bc` dispatch:

```fortran
! Before:
call calc_Ev4<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
call calc_Fv4<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
call calc_Gv4<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)

! After:
if (id_bc_x) then
  call calc_Ev4<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
else
  call calc_Ev4_in<<<blocksEv,threadsEv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, E)
endif
if (id_bc_y) then
  call calc_Fv4<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
else
  call calc_Fv4_in<<<blocksFv,threadsFv>>>(nx, ny, nz, inv_dy, inv_dx, inv_dz, Q, T, mu, F)
endif
if (id_bc_z) then
  call calc_Gv4<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
else
  call calc_Gv4_in<<<blocksGv,threadsGv>>>(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, mu, G)
endif
```

### Makefile Changes

Add `calc_visc4_internal.o` to OBJ in all 6 case Makefiles (`ETGV`, `NSTGV`, `KHI`, `SBLI`, `TBL`, `IVST`):

```makefile
calc_visc4_internal.o: mod_globals.mod mod_constant.mod
```

---

## Critical Files

| File | Change |
|------|--------|
| `3D_solver/src/calc_visc4_internal.f90` | **New file** — interior-only visc4 kernels |
| `3D_solver/src/calc_flux_base.f90` | Add `use calc_visc4_internal`; `id_bc` dispatch at lines ~249–251 |
| `3D_solver/ETGV/Makefile` | Add `calc_visc4_internal.o` |
| `3D_solver/NSTGV/Makefile` | Add `calc_visc4_internal.o` |
| `3D_solver/KHI/Makefile` | Add `calc_visc4_internal.o` |
| `3D_solver/SBLI/Makefile` | Add `calc_visc4_internal.o` |
| `3D_solver/TBL/Makefile` | Add `calc_visc4_internal.o` |
| `3D_solver/IVST/Makefile` | Add `calc_visc4_internal.o` |

---

## Uncommitted Config Change in ETGV

`3D_solver/ETGV/mod_globals.f90` has an uncommitted change (not yet committed):
- `id_scheme`: `real(2)` (SLAU) → `real(4)` (Roe)
- `id_accuracy`: `integer(kind=4)` (4th-order) → `integer(kind=8)` (6th-order)
- `id_tvd`: `integer(kind=8)` (Hybrid TVD) → `integer(kind=2)` (no TVD)

This is a configuration experiment. Commit or revert it separately before starting 3C.

---

## Verification

```bash
cd 3D_solver/NSTGV
make clean && make
bash calc.sh
# Check output VTK: no NaN/Inf; kinetic energy matches baseline

cd 3D_solver/ETGV
make clean && make
bash calc.sh
# Periodic case — exercises the _internal visc4 path

cd 3D_solver/SBLI
make clean && make
bash calc.sh
# Wall-bounded case — exercises the full visc4 path (id_bc = .true.)
```

Profile with `bash profile.sh` (nsys) before and after on NSTGV to quantify improvement in `calc_Ev4` / `calc_Fv4` / `calc_Gv4` kernel time.

---

## Deferred Items

- **1C** — GPU-direct MPI for SBLI/TBL (`id_gpumpi = integer(kind=4)`)
- **2C** — Async I/O: pin host Q, use `cudaMemcpyAsync`
- **4A** — Array layout transposition `Q(5,nx,ny,nz)` → `Q(nx,ny,nz,5)`
