# Plan 4A: Array Layout Transposition

## Status
- LES Bug 6 (Hsgs using molecular instead of turbulent viscosity): confirmed fixed ✅
- `kind(id_accuracy) == 8` guard in viscous `_in` dispatch: intentionally correct ✅

## Problem
`Q(5,nx,ny,nz)`, `E(5,nx-1,ny-2,nz-2)`, `F(5,nx-2,ny-1,nz-2)`, `G(5,nx-2,ny-2,nz-1)` — components are the fastest-varying dimension. Warp threads with consecutive i access Q(l,i1..i32,j,k) at stride-5 (40-byte gaps) → ~10% cache line utilization.

## Fix
Transpose to `Q(nx,ny,nz,5)`, `E(nx-1,ny-2,nz-2,5)`, etc. Consecutive-i warp threads access Q(i1..i32,j,k,l) at stride-1 → fully coalesced.

## Files
All kernel files that access Q/E/F/G:
- `calc_visc4.f90`, `calc_visc4_internal.f90`, `calc_visc4_les_internal.f90`, `calc_visc2.f90`
- `calc_keep_kernel.f90` + `_internal.f90`, `calc_slau_kernel.f90` + `_internal.f90`
- `calc_hybrid_kernel.f90` + `_internal.f90`, `calc_roe_kernel.f90` + `_internal.f90`
- `calc_steps.f90`, `calc_flux_base.f90`
- `calc_physical_quantities.f90`, `calc_muscl.f90`, `print.f90`
- Per-case `set.f90` (initial conditions), `main.f90` (allocation)

## Verification
```bash
cd 3D_solver/NSTGV && make clean && make
bash calc.sh        # correctness: no NaN/Inf, KE matches baseline
bash profile.sh     # performance: kernel times for calc_Ev4_in, calc_keep_x_in
```
