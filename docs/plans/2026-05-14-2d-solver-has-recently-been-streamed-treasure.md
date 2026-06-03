# Bug Fixes for 2D_solver

## Context

The 2D_solver was adapted from the 3D solver by reducing conservative variables from 5 components (ρ, ρu, ρv, ρw, ρE) to 4 (ρ, ρu, ρv, ρE) and removing the G/z-direction flux. All compilation bugs (Bugs 1–8) have been fixed. The remaining issue is a runtime segfault: the simulation crashes when rank 1 (the CPU print rank) tries to write Q00001.vtr. Q00000.vtr is written by `pre_calc` from rank 0's memory and exists; Q00001.vtr is never successfully written.

---

## Status of All Bugs

| # | File | Description | Status |
|---|------|-------------|--------|
| 1 | `2D_solver/src/calc_flux_base.f90:5-7` | Imported undefined 3D symbols: `id_bc_z, blocksG, threadsG, blocksGv, threadsGv` | **FIXED** |
| 2 | `2D_solver/src/calc_keep_kernel.f90:2` | Imported undefined `threadsG` | **FIXED** |
| 3 | `2D_solver/IVST/set.f90:2` and `IVST/data/set.f90:2` | Imported undefined `nz` | **FIXED** |
| 4 | `src/print.f90:229` | `u = QJ(i,3,j)` (copy-paste: u used v's index) | **FIXED** |
| 5 | `2D_solver/IVST/mod_globals.f90:30` and `IVST/data/mod_globals.f90:30` | `dimension = 3` instead of 2, causing `Q(nx,5,ny)` allocation | **FIXED** |
| 6 | `2D_solver/IVST/set.f90:50` and `IVST/data/set.f90:50` | `Jacobian(ny)` — missing x-dimension, should be `Jacobian(nx,ny)` | **FIXED** |
| 7 | `2D_solver/IVST/set.f90:83,89` and `IVST/data/set.f90:83,89` | `do l = 1, 5` out-of-bounds for `Q(nx,4,ny)`; should be `1, 4` | **FIXED** |
| 8 | `2D_solver/IVST/mod_globals.f90` and `IVST/data/mod_globals.f90` | `id_rescale` missing → compile error | **FIXED** |
| **9** | **`2D_solver/src/main.f90:28-34`** | **`allocate(Q, x, dx, y, dy, Jacobian)` + `set_grid` + `set_Jacobian_xy2` inside even-rank-only `if` block → rank 1 has null x, y → segfault in `print_vtk_2D`** | **TO FIX** |

---

## Bug 9 — Segfault: rank 1 never allocates `x`, `y`

### Root cause

In `2D_solver/src/main.f90`, the array allocation and grid setup are mistakenly guarded inside `if (mod(myrank,2) == 0)`:

```fortran
! WRONG (current 2D code)
if (mod(myrank,2) == 0) then
  call set_block2(...)
  allocate(Q(nx,dimension+2,ny), x(nx), dx(nx-1), y(ny), dy(ny-1), Jacobian(nx,ny))
  call set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
  call set_Jacobian_xy2(nx, ny, dx, dy, Jacobian)
endif
```

Rank 1 (print rank) skips this block entirely, so `x`, `y`, `Q`, `Jacobian` remain unallocated. Later, `RungeKutta` is called for all ranks (line 69), and inside `send_recv_for_print_odd2`, it calls:

```fortran
call print_vtk(step, nx, ny, x, y, rho1d, p1d, v1d)
```

`print_vtk_2D` calls `print_xml(... real(x), real(y) ...)`, which dereferences the null base address of the unallocated `x` and `y` → segfault at `(nil)`.

### How the 3D solver avoids this

In `3D_solver/src/main.f90` (line 29–34), only `set_block3` is inside the even-rank guard; allocation and `set_grid` are unconditional:

```fortran
! CORRECT (3D reference)
if (mod(myrank,2) == 0) then
  call set_block3(...)
endif
allocate(Q(nx,dimension+2,ny,nz), x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny))
call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
call set_Jacobian_xy3(nx, ny, nz, dx, dy, dz, Jacobian)
```

### Fix

**File:** `2D_solver/src/main.f90`

Move `allocate`, `set_grid`, and `set_Jacobian_xy2` out of the even-rank guard, matching the 3D pattern. Only `set_block2` stays inside the `if`.

**Before (lines 28–34):**
```fortran
if (mod(myrank,2) == 0) then
  call set_block2(nx, ny, threads, threadsE, threadsEv, threadsF, threadsFv, &
                  blocks, blocksE, blocksEv, blocksF, blocksFv)
  allocate(Q(nx,dimension+2,ny), x(nx), dx(nx-1), y(ny), dy(ny-1), Jacobian(nx,ny))
  call set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
  call set_Jacobian_xy2(nx, ny, dx, dy, Jacobian)
endif
```

**After:**
```fortran
if (mod(myrank,2) == 0) then
  call set_block2(nx, ny, threads, threadsE, threadsEv, threadsF, threadsFv, &
                  blocks, blocksE, blocksEv, blocksF, blocksFv)
endif
allocate(Q(nx,dimension+2,ny), x(nx), dx(nx-1), y(ny), dy(ny-1), Jacobian(nx,ny))
call set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
call set_Jacobian_xy2(nx, ny, dx, dy, Jacobian)
```

No other changes needed: `set_init` (initial conditions) correctly stays inside the even-rank guard, and `deallocate(Q, x, dx, y, dy, Jacobian)` at the end already applies to all ranks and will work correctly once they're all allocated.

---

## Files to Modify

| File | Change |
|------|--------|
| `2D_solver/src/main.f90` | Move `allocate` + `set_grid` + `set_Jacobian_xy2` outside the `if (mod(myrank,2)==0)` block |

## Verification

```bash
cd 2D_solver/IVST
make clean && make        # should compile without errors
bash calc.sh              # should run past the first output step (Q00001.vtr written)
ls data/Q000*.vtr         # confirm Q00001.vtr and Q00002.vtr exist
```
