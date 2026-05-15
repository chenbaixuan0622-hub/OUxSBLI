# Plan: Propagate `real(sp)` precision for the shock sensor

## Context

The Ducros shock sensor is computed at `real(sp)` (`real(4)`, single precision) and stored
in a `real(sp)` device array in `calc_flux_base.f90`. However, several kernel subroutines
that receive the sensor array still declare it as `real(8)`, causing a type mismatch on the
GPU. The user wants all sensor parameters and local sensor-holding variables to use `real(sp)`
consistently.

**Already correct (`real(sp)`):**
- `sensor` array in `calc_flux_base.f90:30`
- `fd` output of `calc_Ducros` in `calc_hybrid.f90:18`
- `sensor` param in `calc_slau_kernel_internal.f90:83,145,207`

**Needs fixing (currently `real(8)`):**
- `calc_hybrid_kernel_internal.f90` — sensor param, `fdx/fdy/fdz` locals, `delta_r*` sensor param
- `calc_hybrid_kernel.f90` — sensor param (merged with Q,T), `fdx/fdy/fdz` locals
- `calc_slau_kernel.f90` — sensor param, `fdx/fdy/fdz` locals
- `calc_slau_kernel_internal.f90` — `fd` in `interp2/4/6`, `fdx/fdy/fdz` locals
- `src/calc_muscl.f90` — `sensor` in `delta4`, `delta6`, and all `MUSCL3rd*/4th*`
- Five `mod_globals.f90` (NSTGV, SBLI, KHI, TBL, IVST) — missing `sp` definition

---

## Changes

### 1. Add `sp` to missing `mod_globals.f90` files
Add `integer, parameter :: sp = kind(1.e0)` to each:
- `3D_solver/NSTGV/mod_globals.f90`
- `3D_solver/SBLI/mod_globals.f90`
- `3D_solver/KHI/mod_globals.f90`
- `3D_solver/TBL/mod_globals.f90`
- `3D_solver/IVST/mod_globals.f90`

### 2. `3D_solver/src/calc_hybrid_kernel_internal.f90`
`sp` is already in scope via `use calc_hybrid` (line 12). No import change needed.
- Lines 80, 150, 215: `real(8), intent(in), device, contiguous :: sensor` → `real(sp), ...`
- Lines 85, 155, 220: Split `real(8) fdx/fdy/fdz` out from `rhol, ul, vl, wl, pl`:
  - `real(sp) fdx` (and similarly fdy, fdz)
  - `real(8) rhol, ul, vl, wl, pl`
- Lines 103, 173, 238: `0.5d0 * (sensor(...) + sensor(...))` → `0.5_sp * (...)`
- Lines 45, 54, 63: `real(8), intent(in), value :: sensor` in `delta_r2/4/6` → `real(sp), ...`

### 3. `3D_solver/src/calc_hybrid_kernel.f90`
- Line 2: Add `sp` to `use mod_globals` import
- Lines 43, 151, 259, 367, 455, 543, 631, 681, 731: Split the combined declaration
  `real(8), intent(in), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)`
  into two lines:
  ```fortran
  real(8),  intent(in), device, contiguous :: Q(nx,5,ny,nz), T(nx,ny,nz)
  real(sp), intent(in), device, contiguous :: sensor(nx,ny,nz)
  ```
- All `real(8) fdx/fdy/fdz` locals → `real(sp) fdx/fdy/fdz`
- All `fdx/fdy/fdz = 0.5d0 * (sensor(...) + sensor(...))` → `0.5_sp * (...)`

### 4. `3D_solver/src/calc_slau_kernel.f90`
- Line 2: Add `sp` to `use mod_globals` import
- All `real(8), intent(in), device, contiguous :: sensor(nx,ny,nz)` → `real(sp), ...`
- All `real(8) fdx/fdy/fdz` locals → `real(sp) fdx/fdy/fdz`
- All `fdx/fdy/fdz = 0.5d0 * (sensor(...) + sensor(...))` → `0.5_sp * (...)` (some already use `0.50_sp` constant but are stored in `real(8)` vars — those are fixed by changing the var type)

### 5. `3D_solver/src/calc_slau_kernel_internal.f90`
`sp` is already in scope via `use calc_hybrid` (line 5).
- Lines 30, 48, 66: `real(8), intent(inout) :: fd` → `real(sp), intent(inout) :: fd`
- Lines 95, 157, 219: `real(8) fdx/fdy/fdz` → `real(sp) fdx/fdy/fdz`

### 6. `src/calc_muscl.f90`
- Module level: Add `use mod_globals, only : sp`
- `delta4` (line 136): `use mod_globals, only : id_tvd` → `use mod_globals, only : id_tvd, sp`
  - Line 137: `real(8), intent(in), value :: sensor` → `real(sp), ...`
- `delta6` (line 149): same pattern for use + sensor param
- All `MUSCL3rdnonTVD`, `MUSCL3rdMinmod`, `MUSCL3rdThreshold`, `MUSCL4thnonTVD`,
  `MUSCL4thTVD`, `MUSCL4thThreshold`: change `real(8), intent(in) :: sensor` → `real(sp), ...`

---

## Critical files
- `3D_solver/src/calc_hybrid_kernel_internal.f90`
- `3D_solver/src/calc_hybrid_kernel.f90`
- `3D_solver/src/calc_slau_kernel.f90`
- `3D_solver/src/calc_slau_kernel_internal.f90`
- `src/calc_muscl.f90`
- `3D_solver/NSTGV/mod_globals.f90`, `3D_solver/SBLI/mod_globals.f90`,
  `3D_solver/KHI/mod_globals.f90`, `3D_solver/TBL/mod_globals.f90`,
  `3D_solver/IVST/mod_globals.f90`

## Verification
```bash
cd 3D_solver/ETGV && make clean && make   # ETGV case (uses Hybrid scheme)
cd 3D_solver/NSTGV && make clean && make  # NSTGV case
```
Both should compile without type-mismatch errors. Check that no `real(8)` sensor param
remains in the relevant subroutines:
```bash
grep -n "real(8).*sensor" 3D_solver/src/calc_hybrid_kernel*.f90 \
     3D_solver/src/calc_slau_kernel*.f90 src/calc_muscl.f90
```
Expected output: empty.
