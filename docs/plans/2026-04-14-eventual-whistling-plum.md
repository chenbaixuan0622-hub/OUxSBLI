# Plan 4A-1: Array Layout — Fix MIXED Viscous Kernel Files

## Current Status

The layout transposition `Q(5,nx,ny,nz)` → `Q(nx,5,ny,nz)` is in progress.

| File | Status |
|------|--------|
| `calc_steps.f90` | ✅ Done |
| `calc_flux_base.f90` | ✅ Done (E/F/G local decls updated) |
| `calc_keep_kernel.f90` + `_internal.f90` | ✅ Done |
| `calc_slau_kernel.f90` + `_internal.f90` | ✅ Done |
| `calc_hybrid_kernel.f90` + `_internal.f90` | ✅ Done |
| `calc_roe_kernel.f90` + `_internal.f90` | ✅ Done |
| `calc_les.f90` | ✅ Done |
| **`calc_visc4.f90`** | ⚠️ MIXED — body has remaining OLD accesses |
| **`calc_visc4_internal.f90`** | ⚠️ MIXED — body has remaining OLD accesses |
| **`calc_visc4_les_internal.f90`** | ⚠️ MIXED — body has remaining OLD accesses |
| **`calc_visc2.f90`** | ⚠️ MIXED — body has remaining OLD accesses |
| `set_bc_common.f90` | 🔲 Phase 2 |
| `calc_rescale.f90` | 🔲 Phase 2 |
| `preprocess.f90` | 🔲 Phase 2 |
| `calc_time_dev.f90` | 🔲 Phase 2 |
| `src/print.f90`, `src/main.f90` | 🔲 Phase 2 |
| Per-case `set.f90` (×6) | 🔲 Phase 2 |
| Per-case `mod_globals.f90` (×6) | 🔲 Phase 2 (declaration flip — last) |

---

## Phase 1 Task: Fix the Four MIXED Files

All four files have the dummy argument declaration already updated to `Q(nx,5,ny,nz)` but the
cross-derivative and boundary stencil accesses in the body still use the old `Q(l,i,j,k)` order.

### Pattern to fix throughout

```fortran
! OLD (still present in body):
Q(2,i,j-1,k)   Q(3,i,j-1,k)   Q(4,i,j-1,k)
Q(2,i,j+1,k)   Q(3,i,j+1,k)   Q(4,i,j+1,k)
Q(2,i,j,k-1)   Q(4,i,j,k-1)   ...
Q(2,i,j,k+1)   Q(4,i,j,k+1)   ...
Q(2,i-2,j,k)   Q(3,i-2,j,k)   ...   (6th-order x ±2 terms)
Q(2,i+2,j,k)   etc.
! Array slices (visc2 only):
Q(2,i,j-1:j+1,k-1:k+1)    Q(3,i,j-1:j+1,k)   etc.

! NEW (correct form):
Q(i,2,j-1,k)   Q(i,3,j-1,k)   Q(i,4,j-1,k)
Q(i,2,j+1,k)   Q(i,3,j+1,k)   Q(i,4,j+1,k)
Q(i,2,j,k-1)   Q(i,4,j,k-1)   ...
Q(i,2,j-2,j,k) → Q(i,2,j,k-2) ...
Q(i-2,2,j,k)   Q(i-2,3,j,k)   ...   (6th-order x ±2 terms)
Q(i+2,2,j,k)   etc.
! Array slices (visc2):
Q(i,2,j-1:j+1,k-1:k+1)    Q(i,3,j-1:j+1,k)   etc.
```

---

### `3D_solver/src/calc_visc4.f90`

**x-kernel** (`calc_Ev4`, lines ~50–55): y and z cross-derivative loads:
```fortran
! OLD:
uy(ii,jt,kt) = (two_third*(-Q(2,i,j-1,k)+Q(2,i,j+1,k)) - one_twelfth*(-Q(2,i,j-2,k)+Q(2,i,j+2,k)))*dy(j)
vy(ii,jt,kt) = (two_third*(-Q(3,i,j-1,k)+Q(3,i,j+1,k)) - one_twelfth*(-Q(3,i,j-2,k)+Q(3,i,j+2,k)))*dy(j)
uz(ii,jt,kt) = (two_third*(-Q(2,i,j,k-1)+Q(2,i,j,k+1)) - one_twelfth*(-Q(2,i,j,k-2)+Q(2,i,j,k+2)))*dz(k)
wz(ii,jt,kt) = (two_third*(-Q(4,i,j,k-1)+Q(4,i,j,k+1)) - one_twelfth*(-Q(4,i,j,k-2)+Q(4,i,j,k+2)))*dz(k)
! NEW:
uy(ii,jt,kt) = (two_third*(-Q(i,2,j-1,k)+Q(i,2,j+1,k)) - one_twelfth*(-Q(i,2,j-2,k)+Q(i,2,j+2,k)))*dy(j)
vy(ii,jt,kt) = (two_third*(-Q(i,3,j-1,k)+Q(i,3,j+1,k)) - one_twelfth*(-Q(i,3,j-2,k)+Q(i,3,j+2,k)))*dy(j)
uz(ii,jt,kt) = (two_third*(-Q(i,2,j,k-1)+Q(i,2,j,k+1)) - one_twelfth*(-Q(i,2,j,k-2)+Q(i,2,j,k+2)))*dz(k)
wz(ii,jt,kt) = (two_third*(-Q(i,4,j,k-1)+Q(i,4,j,k+1)) - one_twelfth*(-Q(i,4,j,k-2)+Q(i,4,j,k+2)))*dz(k)
```

**Boundary region** (lines ~84–96, `my`/`mz` block): MIXED — finish updating all
`Q(2,i,j-1,k)` → `Q(i,2,j-1,k)` etc. Pattern: any remaining `Q(l,i,j,k)` with numeric first arg.

Apply the same fixes to `calc_Fv4` and `calc_Gv4` kernels. For `calc_Fv4` the perpendicular directions are x and z; for `calc_Gv4` they are x and y.

Also apply to `calc_Ev_LES4`, `calc_Fv_LES4`, `calc_Gv_LES4` (LES variants in the same file).

---

### `3D_solver/src/calc_visc4_internal.f90`

Same cross-derivative fixes as `calc_visc4.f90`, plus:

**6th-order x ±2 stencil terms** (e.g. `calc_Fv4_in` y-kernel, lines ~119–120):
```fortran
! OLD (one_twelfth terms):
ux(jj,it,kt) = (two_third*(-Q(i-1,2,j,k)+Q(i+1,2,j,k)) - one_twelfth*(-Q(2,i-2,j,k)+Q(2,i+2,j,k)))*dx(i)
vx(jj,it,kt) = (two_third*(-Q(i-1,3,j,k)+Q(i+1,3,j,k)) - one_twelfth*(-Q(3,i-2,j,k)+Q(3,i+2,j,k)))*dx(i)
! NEW:
ux(jj,it,kt) = (two_third*(-Q(i-1,2,j,k)+Q(i+1,2,j,k)) - one_twelfth*(-Q(i-2,2,j,k)+Q(i+2,2,j,k)))*dx(i)
vx(jj,it,kt) = (two_third*(-Q(i-1,3,j,k)+Q(i+1,3,j,k)) - one_twelfth*(-Q(i-2,3,j,k)+Q(i+2,3,j,k)))*dx(i)
```

Apply to all three direction kernels (`calc_Ev4_in`, `calc_Fv4_in`, `calc_Gv4_in`).

---

### `3D_solver/src/calc_visc4_les_internal.f90`

Identical fixes to `calc_visc4_internal.f90` — same OLD access patterns in y/z cross-derivatives and
x ±2 stencil terms. Apply the same changes to all three direction kernels (and their LES variants if present).

---

### `3D_solver/src/calc_visc2.f90`

**Array slice syntax** (x-kernel, lines ~36–42):
```fortran
! OLD (array slices with component first):
u(it,jt-1:jt+1,kt-1:kt+1) = Q(2,i,j-1:j+1,k-1:k+1)
v(it,jt-1:jt+1,kt)        = Q(3,i,j-1:j+1,k)
w(it,kt-1:kt+1,jt)        = Q(4,i,j,k-1:k+1)
! Also at i+1:
u(it+1,...) = Q(2,i+1,j-1:j+1,k-1:k+1)
! NEW:
u(it,jt-1:jt+1,kt-1:kt+1) = Q(i,2,j-1:j+1,k-1:k+1)
v(it,jt-1:jt+1,kt)        = Q(i,3,j-1:j+1,k)
w(it,kt-1:kt+1,jt)        = Q(i,4,j,k-1:k+1)
u(it+1,...) = Q(i+1,2,j-1:j+1,k-1:k+1)
```

**Boundary region** (lines ~136–158, `my`/`mz` block): finish updating the remaining
`Q(2,i,j-1,k)` etc. accesses. Apply to `calc_Fv2` and `calc_Gv2` (and their LES variants) as well.

---

## Verification (after Phase 2 completes)

```bash
cd 3D_solver/NSTGV && make clean && make && bash calc.sh
# KE decay must match pre-transposition baseline; no NaN/Inf
cd 3D_solver/SBLI  && make clean && make && bash calc.sh
```

---

## Phase 2 (plan-4A-2-host-code.md)

Remaining files: `set_bc_common.f90`, `calc_rescale.f90`, `preprocess.f90`,
`calc_time_dev.f90`, `src/print.f90`, `src/main.f90`, per-case `set.f90` (×6),
and the final `mod_globals.f90` (×6) declaration flip.
