# Plan: Bug fixes and optimizations for z-split MPI overlap implementation

## Context

The z-decomposition with overlapped MPI communication was implemented and tested with STZ (Sod shock tube, 2nd-order Euler, 4 ranks). The test passes. However several bugs lurk in code paths not exercised by STZ, and two performance problems exist even in the Euler path.

---

## Confirmed Bugs

### Bug 1 — `intent(out)` on `G` in all three `_koff` convective kernels
**Files:** [calc_keep_kernel_internal.f90](3D_solver/src/calc_keep_kernel_internal.f90#L170), [calc_slau_kernel_internal.f90](3D_solver/src/calc_slau_kernel_internal.f90#L268), [calc_hybrid_kernel_internal.f90](3D_solver/src/calc_hybrid_kernel_internal.f90#L272)

All three `calc_*_z_in_koff` kernels declare `G` as `intent(out)`. Each only writes to `G(:,:,:,k)` for `k ∈ [k_lo, k_hi]`, leaving the rest untouched. With `intent(out)` the CUDA Fortran compiler may assume the kernel is the sole writer of the entire `G` array and generate incorrect alias-based optimizations. The correct attribute is `intent(inout)`.

**Fix:** Change `intent(out)` → `intent(inout)` for `G` in all three kernels (3 lines total).

---

### Bug 2 — LES interior phase writes viscous fluxes to wrong k-positions
**File:** [calc_flux_base.f90](3D_solver/src/calc_flux_base.f90#L480)  
**Subroutine:** `calc_EFG_z_interior_LES`

The non-koff `calc_Ev_LES2`/`calc_Fv_LES2` map `k = (blockIdx%z-1)*blockDim%z + kt + 1` (k starts at 2). The non-koff `calc_Gv_LES2` maps `k = (blockIdx%z-1)*blockDim%z + kt` (k starts at 1). With restricted interior block counts these produce:

| Kernel | Non-koff computes | Correct interior range | Error |
|---|---|---|---|
| Ev/Fv (overlap=1) | k = 2..nz-3 | k = 3..nz-2 | 1 off each end |
| Gv (overlap=1) | k = 1..nz-3 | k = 2..nz-2 | 1 off each end |
| Ev/Fv (overlap=2) | k = 2..nz-5 | k = 4..nz-3 | 2 off each end |

The non-koff kernels write viscous contributions to E/F/G at the wrong k-positions, corrupting flux values near both boundaries.

---

### Bug 3 — LES halo phase misses the top (hi) halo entirely
**File:** [calc_flux_base.f90](3D_solver/src/calc_flux_base.f90#L533)  
**Subroutine:** `calc_EFG_z_halo_LES`

The halo phase launches each non-koff LES kernel once with `n_halo` blocks. Since those kernels start at k=1 or k=2, they cover the bottom halo region by coincidence (exactly for 2nd order, approximately wrong for 4th+ order), but the **top halo (k = nz-overlap_fb..nz-1) is never computed**. The comment in the code acknowledges this.

Additionally, for 4th-order and above, even the bottom Gv halo is wrong: non-koff `calc_Gv_LES2` with `n_halo` blocks gives k=1..overlap_fb, but the correct Gv lo halo is k=overlap_fb..2*overlap_fb-1.

---

## Fix for Bugs 2 & 3: Add LES `_koff` variants

The only correct fix is to add k-restricted variants of all four LES viscous kernels:
- `calc_Ev_LES2_koff` in [calc_visc2.f90](3D_solver/src/calc_visc2.f90) (and export from module public list)
- `calc_Fv_LES2_koff` in [calc_visc2.f90](3D_solver/src/calc_visc2.f90)
- `calc_Gv_LES2_koff` in [calc_visc2.f90](3D_solver/src/calc_visc2.f90)
- `calc_Gv_LES4_koff` in [calc_visc4.f90](3D_solver/src/calc_visc4.f90) (for 4th-order LES)

Pattern to follow: copy the corresponding LES kernel, change `k = ...+ kt + 1` to `k = ...+ k_lo - 1 + kt` (for Ev/Fv-style) or `k = ...+ k_lo - 1 + kt` with `+0` offset (for Gv-style), add `k_lo`/`k_hi` parameters, and change the active guard from `nz-N < k` to `k_hi < k`.

Then fix `calc_EFG_z_interior_LES`:
```fortran
! Replace:
call calc_Ev_LES2<<<blocksEv_int,threadsEv>>>(...)
call calc_Fv_LES2<<<blocksFv_int,threadsFv>>>(...)
call calc_Gv_LES2<<<blocksGv_int,threadsGv>>>(...)
! With:
call calc_Ev_LES2_koff<<<blocksEv_int,threadsEv>>>(..., overlap_fb+2, nz-overlap_fb-1)
call calc_Fv_LES2_koff<<<blocksFv_int,threadsFv>>>(..., overlap_fb+2, nz-overlap_fb-1)
call calc_Gv_LES2_koff<<<blocksGv_int,threadsGv>>>(..., overlap_fb+1, nz-overlap_fb-1)
```

And fix `calc_EFG_z_halo_LES` — add **two** launches per kernel (lo + hi halo), mirroring the NS halo pattern:
```fortran
! lo halo Ev, Fv, Gv
call calc_Ev_LES2_koff<<<blocksEv_halo,threadsEv>>>(..., 2,         overlap_fb+1)
call calc_Ev_LES2_koff<<<blocksEv_halo,threadsEv>>>(..., nz-overlap_fb, nz-1)
call calc_Fv_LES2_koff<<<blocksFv_halo,threadsFv>>>(..., 2,         overlap_fb+1)
call calc_Fv_LES2_koff<<<blocksFv_halo,threadsFv>>>(..., nz-overlap_fb, nz-1)
call calc_Gv_LES2_koff<<<blocksGv_halo,threadsGv>>>(..., overlap_fb,   2*overlap_fb-1)
call calc_Gv_LES2_koff<<<blocksGv_halo,threadsGv>>>(..., nz-overlap_fb, nz-1)
```
(4th-order LES uses the same pattern with `calc_Gv_LES4_koff`.)

---

## Optimizations

### Opt 1 — Double `calc_quantities_*` per RK stage (most impactful)
**Files:** [calc_flux_base.f90](3D_solver/src/calc_flux_base.f90#L370), [calc_physical_quantities.f90](src/calc_physical_quantities.f90)

`calc_EFG_z_interior_*` calls `calc_quantities_3D`/`calc_quantities_T_3D` for ALL `nz` cells. Then `calc_EFG_z_halo_*` calls it again for all `nz` cells, even though only the ghost cells (`k=1..overlap_fb` and `k=nz-overlap_fb+1..nz`) changed between the two calls.

**Fix:** Add `calc_quantities_3D_koff(nx, ny, nz, Jacobian, QJ, Q, T, k_lo, k_hi)` and `calc_quantities_T_3D_koff(...)` in [calc_physical_quantities.f90](src/calc_physical_quantities.f90). The `_koff` variants launch with a z-restricted block grid (same pattern as other koff kernels). In the halo phase, call them twice (lo ghost + hi ghost) instead of once for the full domain:
```fortran
! Replace full-domain call in halo phase:
call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T)
! With two ghost-cell-only calls:
call calc_quantities_3D_koff(nx, ny, nz, Jacobian, QJ, Q, T, 1,              overlap_fb)
call calc_quantities_3D_koff(nx, ny, nz, Jacobian, QJ, Q, T, nz-overlap_fb+1, nz)
```
This eliminates ~(nz-2*overlap_fb)/nz ≈ 99% of the second `calc_quantities` call per RK stage.

---

### Opt 2 — Interior phase computes halo-k convective G faces (wasted work)
**File:** [calc_flux_base.f90](3D_solver/src/calc_flux_base.f90#L368)

`calc_conv(id_scheme, ...)` computes G for the full valid range `k = overlap_fb..nz-overlap_fb`. The halo-k faces (`k = overlap_fb..2*overlap_fb-1` and `k = nz-2*overlap_fb+1..nz-overlap_fb`) are immediately overwritten by the halo phase. These faces can be skipped in the interior phase.

**Fix:** Replace the `calc_conv` call in the interior subroutines with separate scheme-dispatched calls for E, F, and interior-only G. The interior G should use a `_koff` call with `k_lo = 2*overlap_fb, k_hi = nz-2*overlap_fb`. This requires adding an interior-G-only dispatch interface (similar to the existing `calc_conv_G_halo_koff`) or a simple integer(2)/real(2)/real(8) dispatch inline in each `_z_interior_*` subroutine. The savings scale as `2*overlap_fb / nz_local` per RK stage.

---

## Critical Files

| File | Change |
|------|--------|
| [3D_solver/src/calc_keep_kernel_internal.f90](3D_solver/src/calc_keep_kernel_internal.f90) | Bug 1: `intent(out)` → `intent(inout)` on G |
| [3D_solver/src/calc_slau_kernel_internal.f90](3D_solver/src/calc_slau_kernel_internal.f90) | Bug 1: same |
| [3D_solver/src/calc_hybrid_kernel_internal.f90](3D_solver/src/calc_hybrid_kernel_internal.f90) | Bug 1: same |
| [3D_solver/src/calc_visc2.f90](3D_solver/src/calc_visc2.f90) | Bugs 2&3: add `calc_Ev_LES2_koff`, `calc_Fv_LES2_koff`, `calc_Gv_LES2_koff` |
| [3D_solver/src/calc_visc4.f90](3D_solver/src/calc_visc4.f90) | Bugs 2&3: add `calc_Gv_LES4_koff` |
| [3D_solver/src/calc_flux_base.f90](3D_solver/src/calc_flux_base.f90) | Bugs 2&3: fix `calc_EFG_z_interior_LES` and `calc_EFG_z_halo_LES` |
| [src/calc_physical_quantities.f90](src/calc_physical_quantities.f90) | Opt 1: add `calc_quantities_3D_koff`, `calc_quantities_T_3D_koff` |
| [3D_solver/src/calc_flux_base.f90](3D_solver/src/calc_flux_base.f90) | Opt 1: use koff quantities in all halo subroutines |

---

## Recommended Order

1. **Bug 1** (3 one-line fixes, zero risk, do first)
2. **Opt 1** (`calc_quantities_*_koff`) — highest performance impact, needed before shipping NS/LES cases
3. **Bugs 2 & 3** (LES koff variants) — correctness for LES mode
4. **Opt 2** (skip interior halo-k convective G) — lower priority, small savings for large grids

---

## Verification

- **Bug 1:** Recompile STZ case, confirm same results. No behavioral change expected (CUDA Fortran doesn't initialize `intent(out)` device arrays), but removes incorrect aliasing assumption.
- **Opt 1:** Add a timing print around `calc_EFG_z_halo` before/after, run STZ; halo phase time should drop noticeably per step.
- **Bugs 2 & 3:** Switch STZ to NS (`id_visc = integer(4)`) and verify identical output vs non-zdec RungeKutta path. Then switch to LES and verify against non-zdec reference.
- **Opt 2:** Profile with `nsys` — confirm G-kernel time in interior phase is reduced proportional to `2*overlap_fb/nz`.
