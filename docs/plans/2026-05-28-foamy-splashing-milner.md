# Bug Survey — Remaining Files

## Context

Continuation of the prior bug survey. The previous plan left these areas unsurveyed:
- `3D_solver/src/calc_visc4.f90.fypp`, `calc_visc4_internal.f90.fypp`
- `3D_solver/src/calc_*_kernel*.f90.fypp` (KEEP, SLAU, Roe, Hybrid — all 8 files)
- `3D_solver_curv/src/calc_visc2_curv.f90` lines 300+ (LES/SGS section)
- `set.f90` for 3D cases: ETGV, DHIT, SBLI, TBL, NSTGV
- `src/set_compressible_bl.f90`

All files have now been read. Findings below.

---

## Confirmed Bug

### BUG 5 — `3D_solver/src/calc_roe_kernel_internal.f90.fypp`: flux array dimensions transposed

**Lines:** 64, 123, 177

```fortran
! Wrong — declared as (nx-1, 5, ny-2, nz-2):
real(8), intent(out), device, contiguous :: E(nx-1,5,ny-2,nz-2)
real(8), intent(out), device, contiguous :: F(nx-2,5,ny-1,nz-2)
real(8), intent(out), device, contiguous :: G(nx-2,5,ny-2,nz-1)
```

Every other flux kernel (KEEP, SLAU, Hybrid — both regular and internal variants) puts the 5-component variable index **first**:

```fortran
! Correct — from calc_keep_kernel_internal.f90.fypp line 34:
real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2)
```

The access pattern inside the Roe subroutines confirms the first-index-is-5 convention is expected:
```fortran
! Line 109 — accesses component 1..5 as first index:
call Roe(..., E(1,i,j-1,k-1), E(2,i,j-1,k-1), E(3,i,j-1,k-1), E(4,i,j-1,k-1), E(5,i,j-1,k-1))
```
With the wrong declaration `E(nx-1,5,...)`, this write scatters flux components into entirely wrong memory locations.

**Fix:** Change the three declarations:
```fortran
real(8), intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2)
real(8), intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2)
real(8), intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1)
```

**Note:** The Roe scheme is not currently used in production (CLAUDE.md), but the interior-only kernel is structurally broken and would produce silently wrong fluxes if enabled.

---

## Clean Files (No Bugs Found)

| File | Status |
|------|--------|
| `calc_keep_kernel.f90.fypp` | Clean |
| `calc_keep_kernel_internal.f90.fypp` | Clean |
| `calc_slau_kernel.f90.fypp` | Clean |
| `calc_slau_kernel_internal.f90.fypp` | Clean |
| `calc_hybrid_kernel.f90.fypp` | Clean |
| `calc_hybrid_kernel_internal.f90.fypp` | Clean |
| `calc_roe_kernel.f90.fypp` | Clean |
| `calc_visc2_curv.f90` (LES section, lines 300+) | Clean |
| `src/set_compressible_bl.f90` | Clean |
| `3D_solver/ETGV/set.f90` | Clean |
| `3D_solver/DHIT/set.f90` | Clean |
| `3D_solver/SBLI/set.f90` | Clean |
| `3D_solver/TBL/set.f90` | Clean |
| `3D_solver/NSTGV/set.f90` | Clean |

---

## False Positives Investigated and Ruled Out

| Reported | Verdict |
|----------|---------|
| `calc_visc2_curv.f90` line 711: extra `/dz` in SGS heat flux (`calc_Gv_LES2_curv`) | Not a bug — G-flux uses physical z spacing `dz`; E/F fluxes use unit computational spacing Δξ=Δη=1, so the asymmetry is correct |
| `SBLI/set.f90` line 112: `Jacobian(nx,ny,nz)` instead of `(nx,ny)` | Not a bug — `set_bc` correctly declares `Jacobian(nx,ny)` |
| `calc_roe_kernel_internal.f90.fypp`: `id_accuracy` not passed as subroutine argument | Not a bug — imported as a module-level compile-time constant from `mod_constant`; accessible in device code |

---

## Code Quality Notes (No Runtime Impact)

- **`calc_visc4.f90.fypp` lines 484–485**: `mz` and `kTz` are computed, then a cross-stencil block runs (lines 486–515), then `mz` and `kTz` are computed again identically (lines 516–517). The first pair (484–485) is dead code — the values are overwritten before use. No correctness impact; the second pair at 516–517 is the one that feeds `muz`, `mvz`, `mwz`.
- **`calc_visc4_internal.f90.fypp` lines 145–146**: Comments on `dx` and `dy` parameters are swapped (`dx` is labelled "inverse grid spacing in y" and vice versa). No runtime impact.

---

## Verification

After applying the fix:
```bash
cd 3D_solver/ETGV   # any case — need to verify Roe kernels compile cleanly
# Switch config.fypp to SCHEME = 'Roe', then:
cmake --build build -j
# Confirm no compile errors; check that ETGV energy decay is physically correct
```
Since the Roe scheme is not used in current test suite (`pytest ouxsbli/tests/`), compilation is the primary check. The existing test suite (test_etgv, test_evc, test_os, test_corn) uses KEEP/SLAU/Hybrid and will continue to pass unaffected.
