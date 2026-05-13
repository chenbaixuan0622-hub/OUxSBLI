# Layout Transposition: `Q(5,nx,ny,nz)` → `Q(nx,5,ny,nz)` (Plan 4A)

## Context

All Phase 3 items (3A/3B/3C) are committed. In-progress work (uncommitted) implements:

- `Q(5,nx,ny,nz)` → `Q(nx,5,ny,nz)` — component dim moves from 1st to 2nd
- `E(5,nx-1,ny-2,nz-2)` → `E(nx-1,5,ny-2,nz-2)`, same for `F`, `G`, `Rs`, `QJ`, `QJ2`, `QJs`

**Why:** With old layout, threads loading component `l` across x-positions `i,i+1,...` access stride-5 memory. New layout `Q(i,l,j,k)` makes that load stride-1 — fully coalesced.

---

## Fully Updated (Clean)

| File | Status |
|------|--------|
| `calc_flux_base.f90` | ✅ Wrapper arg declarations |
| `calc_keep_kernel.f90`, `_internal.f90` | ✅ |
| `calc_slau_kernel.f90`, `_internal.f90` | ✅ |
| `calc_roe_kernel.f90` | ✅ `E(i,1,...)` element writes correct |
| `calc_hybrid_kernel.f90`, `_internal.f90` | ✅ |
| `calc_steps.f90` | ✅ E/F/G/Q/Q2/Rs decls + accesses |
| `calc_les.f90` | ✅ Q section loads (`Q(i-1:i+1,2,...)`) |
| `calc_visc4_internal.f90` | ✅ All Q loads + E/F/G writes |
| `calc_visc4_les_internal.f90` | ✅ All Q loads + E/F/G writes |
| `calc_roe_kernel_internal.f90` — x/y kernels | ✅ E (line 107) + F (line 158) |
| `calc_physical_quantities.f90` | ✅ 3D kernels all updated |
| `src/print.f90` | ✅ 3D QJ reads updated |
| `src/sbli.f90` — allocation | ✅ `allocate(Q(nx,5,ny,nz), ...)` |
| `set_bc_common.f90` — non-init + device | ✅ Explicit `l`-loop variants correct |
| `preprocess.f90` | ✅ |

---

## Remaining Bugs

### Bug 1: `calc_visc2.f90` — 3 kernels with component-first E/F/G writes

Kernels `calc_Ev2`, `calc_Fv_LES2`, `calc_Gv_LES2` still write component-first:

```fortran
! WRONG — E declared as E(nx-1,5,ny-2,nz-2):
E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx    ! lines 94–97
! CORRECT:
E(i,2,j-1,k-1) = E(i,2,j-1,k-1) - txx
```
Same pattern: `F(2,i-1,j,k-1)` → `F(i-1,2,j,k-1)` (lines 334–337) and `G(2,i-1,j-1,k)` → `G(i-1,2,j-1,k)` (lines 488–491).

### Bug 2: `calc_visc4.f90` — 3 LES kernels with component-first E/F/G writes

Same component-first write pattern in `calc_Ev_LES4` (lines 244–247), `calc_Fv_LES4` (482–485), `calc_Gv_LES4` (719–722). Same fix as Bug 1.

### Bug 3: `calc_roe_kernel_internal.f90` — G write in `calc_roe_z_in` (line 213)

```fortran
! WRONG — G declared as G(nx-2,5,ny-2,nz-1):
call Roe(..., G(1,i-1,j-1,k), G(4,i-1,j-1,k), G(2,i-1,j-1,k), G(3,i-1,j-1,k), G(5,i-1,j-1,k))
! CORRECT (note: z-kernel arg order is rho,w,u,v,E so indices 1,4,2,3,5):
call Roe(..., G(i-1,1,j-1,k), G(i-1,4,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,5,j-1,k))
```

### Bug 4: `calc_time_dev.f90` — QJ2 and QJs allocated with old layout (lines 61, 149, 247, 333)

```fortran
! WRONG:
allocate(QJ(nx,5,ny,nz), QJ2(5,nx,ny,nz), ...)  ! lines 61, 149
allocate(QJ(nx,5,ny,nz), QJs(5,nx,ny,nz), ...)  ! lines 247, 333
! CORRECT:
allocate(QJ(nx,5,ny,nz), QJ2(nx,5,ny,nz), ...)
allocate(QJ(nx,5,ny,nz), QJs(nx,5,ny,nz), ...)
```
QJ2/QJs are passed to GPU kernels (calc_step1, calc_EFG, calc_step2_3) that expect new layout — mismatched allocation causes wrong data interpretation.

### Bug 5: `set_bc_common.f90` — `_init` subroutines have `:` in wrong position (lines 9–93)

The `_init` variants use Fortran array section notation but with `:` in dim 1 (the x-position dim in new layout) instead of dim 2 (the component dim):

```fortran
! WRONG — Q declared as Q(nx,5,ny,nz):
! Q(:,1,j,k) means all-x, comp=1, j, k — copies rho across all x (semantically wrong)
! Q(:,nx-1,j,k) means all-x, comp=nx-1 — OUT OF BOUNDS (comp dim is only 5)
Q(:,1,j,k)  = Q(:,nx-1,j,k)    ! set_bc_cyclic2_init line 16
Q(:,nx,j,k) = Q(:,2,j,k)

! CORRECT — swap : from dim 1 to dim 2, spatial index from dim 2 to dim 1:
Q(1,:,j,k)  = Q(nx-1,:,j,k)    ! i=1, all-comps, j, k
Q(nx,:,j,k) = Q(2,:,j,k)
```

Applies to all three `_init` subroutines (2nd, 4th, 6th order variants), in all three directions (x, y, z) and corners. Pattern: swap `:` and the spatial range index — positions 1 and 2 in the subscript list.

4th-order example:
```fortran
! WRONG:  Q(:,1:2,j,k)   = Q(:,nx-3:nx-2,j,k)
! CORRECT: Q(1:2,:,j,k)   = Q(nx-3:nx-2,:,j,k)
```

### Bug 6: `set_bc_tbl_sbli.f90` — old-layout QJ accesses + spurious `nx,` prefix

Two sub-issues:

**(a)** QJ element accesses still component-first (old layout):
```fortran
! WRONG:
over_QJ1 = 1.d0 / QJ(1,i,ny-1,k)   ! line 17
rhob = QJ(1,i,ny-1,k) * Jacobian(i,ny-1)
ub   = QJ(2,i,ny-1,k) * over_QJ1   ! etc.

! CORRECT:
over_QJ1 = 1.d0 / QJ(i,1,ny-1,k)
rhob = QJ(i,1,ny-1,k) * Jacobian(i,ny-1)
ub   = QJ(i,2,ny-1,k) * over_QJ1
```
Affects lines 17–26, 29–32, 56, 58, 65–66, and similar patterns throughout.

**(b)** Spurious `nx,` prefix in energy component accesses (lines 33, 54):
```fortran
! WRONG (5 indices for a 4D array — nonsensical):
p_wall = gamma_1 * (QJ(nx,5,i,2,k) - ...)   ! line 33
pin    = gamma_1 * (QJ(nx,5,i,ny-1,k) - ...) ! line 54

! CORRECT (new layout: QJ(i,comp,j,k)):
p_wall = gamma_1 * (QJ(i,5,2,k) - ...)
pin    = gamma_1 * (QJ(i,5,ny-1,k) - ...)
```

### Bug 7: `calc_rescale.f90` — QJ element accesses still old layout (lines 25–31, 47–53)

```fortran
! WRONG:
rhoinv = 1.d0 / QJ(1,i,j,k)
Q1 = Q1 + QJ(1,i,j,k) * Jacobian_tmp
Q2 = Q2 + QJ(2,i,j,k) * rhoinv
! CORRECT:
rhoinv = 1.d0 / QJ(i,1,j,k)
Q1 = Q1 + QJ(i,1,j,k) * Jacobian_tmp
Q2 = Q2 + QJ(i,2,j,k) * rhoinv
```
Same for components 3, 4 and the energy expression with squares.

### Bug 8: `set_init_common.f90` — host-side Q initialization (lines 279–292)

```fortran
! WRONG:
Q(1,i,j,k) = p0 / (R * (T(j) + Tstd(i,j,k)))   ! line 279
Q(2,i,j,k) = Q(1,i,j,k) * (u(j) + ustd(i,j,k)) ! line 280
! ...
Q(1,:,1,:) = Q(1,:,2,:)   ! line 287 — all i, comp=1, j=1 (WRONG direction)
! CORRECT:
Q(i,1,j,k) = p0 / (R * (T(j) + Tstd(i,j,k)))
Q(i,2,j,k) = Q(i,1,j,k) * (u(j) + ustd(i,j,k))
! ...
Q(:,1,1,:) = Q(:,1,2,:)   ! all i, comp=1, j=1 (but check semantics — this may need a loop)
```

### Bug 9: `calc_para.f90` — MPI halo exchange with old-layout Q accesses + stale declaration

**(a)** Q accesses in packing/unpacking still component-first (lines 31–32, 54, 76, 94, 113–114, 132, 150, 168):
```fortran
! WRONG:
Q1d_left(...+l) = Q(l, overlap+i, j+1, k+3)
! CORRECT:
Q1d_left(...+l) = Q(overlap+i, l, j+1, k+3)
```

**(b)** Local declaration in `reconstruct_sbli_inlet` (line 158):
```fortran
! WRONG:
real(8), intent(inout), device :: Q(5,nx,ny2,nz)
! CORRECT:
real(8), intent(inout), device :: Q(nx,5,ny2,nz)
```

### Bug 10: `src/sbli.f90` — checkpoint write loop (line 104)

```fortran
! WRONG — outer loop is m=1..5 (component), inner are i,j,l (spatial):
do m = 1, 5
  Q(m,i,j,l) = Jacobian(i,j) * Q(m,i,j,l)
! CORRECT:
  Q(i,m,j,l) = Jacobian(i,j) * Q(i,m,j,l)
```

---

## Implementation Order (remaining work)

1. Fix Bug 1 — `calc_visc2.f90`: 3 × 4 lines (Ev2, Fv_LES2, Gv_LES2 writes)
2. Fix Bug 2 — `calc_visc4.f90`: 3 × 4 lines (Ev_LES4, Fv_LES4, Gv_LES4 writes)
3. Fix Bug 3 — `calc_roe_kernel_internal.f90` line 213: 5 G elements
4. Fix Bug 4 — `calc_time_dev.f90`: change `QJ2(5,...)` → `(nx,5,...)` and `QJs(5,...)` → `(nx,5,...)` at lines 61, 149, 247, 333
5. Fix Bug 5 — `set_bc_common.f90` `_init` variants: swap `:` to position 2 across all 3 subroutines (x/y/z + corners)
6. Fix Bug 6 — `set_bc_tbl_sbli.f90`: old-layout QJ accesses + `QJ(nx,5,i,j,k)` 5-index correction
7. Fix Bug 7 — `calc_rescale.f90`: QJ element accesses
8. Fix Bug 8 — `set_init_common.f90`: Q init lines 279–292
9. Fix Bug 9 — `calc_para.f90`: Q accesses + local declaration
10. Fix Bug 10 — `src/sbli.f90` line 104

---

## Verification

```bash
cd 3D_solver/NSTGV   # periodic NS — all BC, convective, viscous kernels
make clean && make
bash calc.sh          # check first 10-step output against pre-change baseline

cd 3D_solver/SBLI    # wall-bounded — tests set_bc_tbl_sbli, calc_rescale, calc_para
make clean && make
bash calc.sh
```

- Kinetic energy at first output must match pre-change baseline within round-off
- No NaN/Inf in VTK output
- `bash profile.sh` (nsys) to confirm memory coalescing improvement in Q loads
