# Layout Transposition: Plan 4A (complete)

## Context

Q conservative variable array transposed from `Q(5,nx,ny,nz)` to `Q(nx,5,ny,nz)` throughout the solver for GPU memory coalescing. Flux arrays similarly reordered to `E(nx-1,5,ny-2,nz-2)`, `F(nx-2,5,ny-1,nz-2)`, `G(nx-2,5,ny-2,nz-1)`.

All bugs have been fixed. Changes are staged or unstaged but present in the working tree.

---

## All Fixed Files

| File | Notes |
|------|-------|
| `3D_solver/src/calc_para.f90` | Q decls/accesses; 1D formula; loop reorder; QJ decls |
| `3D_solver/src/calc_rescale.f90` | QJ decls; `calc_mean`/`copy` accesses |
| `3D_solver/src/preprocess.f90` | E/F/G alloc; QJ decl; Q access |
| `3D_solver/src/set_bc_common.f90` | All `_init` variants; numbered variants with `device` + `!$cuf kernel` restored; `set_bc_cyclic_z` QJ decl |
| `3D_solver/src/set_bc_tbl_sbli.f90` | QJ decl + accesses |
| `3D_solver/src/set_init_common.f90` | Q decl + accesses |
| `3D_solver/src/calc_visc2.f90` | E/F/G writes; Bug F2: SGS stress fixed in `calc_Gv_LES2` (lines 480–482) |
| `3D_solver/src/calc_visc4.f90` | E/F/G writes; Bug F1: swapped Q indices fixed in `calc_Fv4` (lines 286–287) |
| `3D_solver/src/calc_flux_base.f90` | Bug F3: `id_bc_y`/`id_bc_z` restored in `calc_conv_roe` (lines 144, 149) |
| `3D_solver/src/calc_roe_kernel_internal.f90` | G write |
| `3D_solver/src/calc_slau_3d.f90` | SLAU1/HRSLAU2 now use 5 scalar `intent(out)` args |
| `3D_solver/src/calc_slau_kernel.f90` | All 9 SLAU call sites use individual element refs |
| `3D_solver/src/calc_hybrid_kernel.f90` | KEEP uses array assignment (valid); SLAU call sites use individual element refs |
| `3D_solver/src/calc_keep_kernel.f90` | Array assignment `E(i,:,j,k) = KEEP(...)` — valid, no change needed |
| `3D_solver/src/calc_keep_3d.f90` | KEEP2/4/6 return `F(5)` via array assignment — valid, no change needed |
| `3D_solver/src/calc_time_dev.f90` | QJ2/QJs allocations |
| `3D_solver/NSTGV/set.f90` | `set_init` + `set_bc` Q layout |
| `3D_solver/ETGV/set.f90` | `set_init` + `set_bc` Q layout |
| `3D_solver/KHI/set.f90` | `set_init` + `set_bc` Q layout |
| `3D_solver/TBL/set.f90` | `set_init` Q decl + `set_bc` full QJ sweep |
| `3D_solver/SBLI/set.f90` | `set_init` Q decl + all subroutines with device QJ |
| `3D_solver/IVST/set.f90` | `set_init` Q decl + array-slice accesses + `Qc(nx,5)` + all element accesses |
| `src/print.f90` | QJ/Q decls and accesses |
| `src/sbli.f90` | Q allocation + checkpoint write |

---

## Verification

```bash
cd 3D_solver/NSTGV   # periodic NS — tests cyclic BCs, all convective/viscous kernels
make clean && make
bash calc.sh          # kinetic energy at first output must match pre-change baseline

cd 3D_solver/TBL     # turbulent boundary layer — tests set_bc QJ kernel, Qre unpack
make clean && make
bash calc.sh

cd 3D_solver/SBLI    # wall-bounded with shock — tests set_bc_tbl_sbli, calc_rescale, calc_para
make clean && make
bash calc.sh
```

- No NaN/Inf in VTK output
- Kinetic energy matches pre-change baseline within round-off
- `bash profile.sh` (nsys) to confirm coalescing improvement in Q/QJ loads
