# Plan: 4th-order G-flux in z for cyclic BCs (NS and LES koff paths)

## Context

`calc_Gv4_koff` falls back to 2nd-order for G-flux faces where k < 3 or k > nz-3, because the 4th-order z-stencil would index outside the array. For non-periodic z BCs those faces are real walls and 2nd-order is necessary. But for cyclic (periodic) z BCs, those faces are between ghost cells and do not need to be computed at all — the correct action is to skip them, not downgrade to 2nd-order. The non-koff path already does this correctly via `calc_Gv4_in` / `calc_Gv_LES4_in`, but those kernels are gated behind `kind(id_accuracy) == 8` even though `kind(id_accuracy) == 4` provides 2 ghost cells — exactly what the 4th-order viscous stencil requires. The koff path has no `_in_koff` G-flux variant at all.

## Changes

### 1. `3D_solver/src/calc_visc4_internal.f90`

Add `calc_Gv4_in_koff` to the `public` list and implement it. This is `calc_Gv4_in` (lines 132–184) with the koff machinery added:
- Add `k_lo, k_hi` value parameters
- Replace `k = (blockIdx%z-1)*blockDim%z + kt` with `k = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt`
- Replace `call load_smem_visc4_z(...)` with `call load_smem_visc4_z_koff(..., k_lo)` (already public in `load_smem_visc4`)
- Add early exit `if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return`
- Keep the 4th-order `if` block unchanged — **no `else` branch** (unlike `calc_Gv4_koff`)

### 2. `3D_solver/src/calc_visc4_les_internal.f90`

Add `calc_Gv_LES4_in_koff` to the `public` list and implement it. Same pattern as above but modeled on `calc_Gv_LES4_in` (lines 146–206), with `mut, qc2` parameters and `calc_tau_straight_LES` / `calc_tau_cross_LES` / LES enthalpy term.

### 3. `3D_solver/src/calc_flux_base.f90`

Six edits — all within the 4th-order accuracy branch (`kind(id_accuracy) >= 4`):

| Location | Current | Change |
|---|---|---|
| Line 399 (NS non-koff, G) | `kind(id_accuracy) == 8` | `kind(id_accuracy) >= 4` |
| Line 477 (LES non-koff, G) | `kind(id_accuracy) == 8` | `kind(id_accuracy) >= 4` |
| Lines 432–433 (NS koff interior G) | unconditional `calc_Gv4_koff` | `if (.not. id_bc_z)` → `calc_Gv4_in_koff`, `else` → `calc_Gv4_koff` |
| Lines 593–596 (NS koff halo G ×2) | unconditional `calc_Gv4_koff` | same if/else |
| Lines 512–513 (LES koff interior G) | unconditional `calc_Gv_LES4_koff` | `if (.not. id_bc_z)` → `calc_Gv_LES4_in_koff`, `else` → `calc_Gv_LES4_koff` |
| Lines 659–662 (LES koff halo G ×2) | unconditional `calc_Gv_LES4_koff` | same if/else |

`id_bc_z` is already imported from `mod_globals` at line 5.

## Why this is correct

- `kind(id_accuracy) >= 4` means `overlap_fb >= 2`: at least 2 ghost cells in z, exactly what the 4th-order viscous z-stencil needs. (`kind == 2` gives only 1 ghost cell — insufficient.)
- For periodic z, the halo koff calls (`k_lo=1, k_hi=overlap_fb`) land on ghost-zone G-flux faces (k<3). `calc_Gv4_in_koff` simply skips them (no write), leaving those G array entries at zero, which is correct since the flux divergence only loops over real cells.
- For non-periodic z, `calc_Gv4_koff` (with its 2nd-order fallback) is still used unchanged.
- No changes to `id_visc` value or `kind(id_accuracy)` in `mod_globals.f90` are required.

## Verification

After building in a case with cyclic z BCs and `kind(id_accuracy) == 4` (e.g., NSTGV or KHI):
1. Check that the build succeeds with no new warnings.
2. Run a short simulation and confirm the solution is smooth (no checkerboard or z-boundary artifacts).
3. Compare energy/enstrophy time traces against the `kind(id_accuracy) == 8` reference to confirm the G-flux now uses 4th-order instead of 2nd-order at the z ghost-zone boundaries.
