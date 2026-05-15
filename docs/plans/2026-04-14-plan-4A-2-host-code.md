# Plan 4A-2: Array Layout — Host/Utility Code + Declaration Flip

## Context

Phase 2 of 2 for `Q(5,nx,ny,nz)` → `Q(nx,5,ny,nz)` transposition.
Phase 1 (`eventual-whistling-plum.md`) covers the GPU kernel files.
This phase covers host-side code and the final declaration flip.
The binary is not buildable until this phase completes.

Index change throughout: `Q(l,i,j,k)` → `Q(i,l,j,k)`, and same for `QJ`, `QJ2`, `QJs`, `ruvwp`, `E`, `F`, `G`, `Rs`.

---

## Files to Modify

### `3D_solver/src/calc_rescale.f90`
Uses `QJ(5,nx,ny,nz)` with accesses `QJ(1,i,j,k)` through `QJ(5,i,j,k)`.
- Dummy arg: `QJ(5,nx,ny,nz)` → `QJ(nx,5,ny,nz)`
- All `QJ(l,i,j,k)` → `QJ(i,l,j,k)`

### `3D_solver/src/set_bc_common.f90`
Declares `Q(5,nx,ny,nz)` and uses whole-component-slice syntax `Q(:,i,j,k)`.
- Dummy arg: `Q(5,nx,ny,nz)` → `Q(nx,5,ny,nz)`
- Change `Q(:,i,j,k)` → `Q(i,:,j,k)` throughout (copy between boundary points)
  ```fortran
  ! OLD:
  Q(:,1,j,k)  = Q(:,nx-1,j,k)
  Q(:,nx,j,k) = Q(:,2,j,k)
  ! NEW:
  Q(1,:,j,k)  = Q(nx-1,:,j,k)
  Q(nx,:,j,k) = Q(2,:,j,k)
  ```

### `3D_solver/src/preprocess.f90`
Two subroutines to update:

**`allocate_device_mem`** (line ~48): allocation of E/F/G and ruvwp:
```fortran
! OLD:
allocate(ruvwp(5,nx,ny,nz), E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1), ...)
! NEW:
allocate(ruvwp(nx,5,ny,nz), E(nx-1,5,ny-2,nz-2), F(nx-2,5,ny-1,nz-2), G(nx-2,5,ny-2,nz-1), ...)
```

**Host Q Jacobian scaling** (line ~82, 91): dummy arg and loop body:
```fortran
! OLD decls:
real(8), intent(inout) :: Q(5,nx,ny,nz)
real(8), intent(out), device, contiguous :: QJ(5,nx,ny,nz)
! NEW:
real(8), intent(inout) :: Q(nx,5,ny,nz)
real(8), intent(out), device, contiguous :: QJ(nx,5,ny,nz)
! OLD body (line ~102):
Q(l,i,j,k) = Q(l,i,j,k) / Jacobian_cpu(i,j)
! NEW:
Q(i,l,j,k) = Q(i,l,j,k) / Jacobian_cpu(i,j)
```

### `3D_solver/src/calc_time_dev.f90`
Three to four subroutines each with the same pattern.

**Dummy arg declarations** (lines 44, 128, 228, 310):
```fortran
! OLD:  real(8), intent(inout) :: Q(5,nx,ny,nz)
! NEW:  real(8), intent(inout) :: Q(nx,5,ny,nz)
```

**Device array allocations** (lines ~61, 149, 247, 333):
```fortran
! OLD:
allocate(QJ(5,nx,ny,nz), QJ2(5,nx,ny,nz), ...)
allocate(QJ(5,nx,ny,nz), QJs(5,nx,ny,nz), Rs(5,nx-2,ny-2,nz-2))
! NEW:
allocate(QJ(nx,5,ny,nz), QJ2(nx,5,ny,nz), ...)
allocate(QJ(nx,5,ny,nz), QJs(nx,5,ny,nz), Rs(nx-2,5,ny-2,nz-2))
```

Any direct QJ/Q index access (e.g. D→H copy, print path):
- `QJ(l,i,j,k)` → `QJ(i,l,j,k)` etc.

### `src/print.f90`
Host-side VTK output. Reads Q (or a host copy of QJ).
- Dummy arg declaration: `Q(5,nx,ny,nz)` → `Q(nx,5,ny,nz)`
- All `Q(l,i,j,k)` → `Q(i,l,j,k)` in the output loops

### `src/main.f90`
Host Q allocation (line ~36):
```fortran
! OLD:
allocate(Q(dimension+2,nx,ny,nz), ...)
! NEW:  (dimension+2 = 5, becomes second index)
allocate(Q(nx,dimension+2,ny,nz), ...)
```
Any direct Q index access in main → update.

### Per-case `set.f90` (×6: ETGV, NSTGV, KHI, SBLI, TBL, IVST)
Each sets ICs and calls BCs using `Q(5,nx,ny,nz)` with `Q(1,i,j,k)` style.
- Dummy arg: `Q(5,nx,ny,nz)` → `Q(nx,5,ny,nz)`
- All `Q(l,i,j,k)` → `Q(i,l,j,k)`

### Per-case `mod_globals.f90` (×6) — **do last, after all callers updated**
This is the actual declaration flip that makes the new layout take effect.
Each `mod_globals.f90` declares the main device arrays. Find and change:
```fortran
! OLD:
real(8), allocatable, device :: QJ(5,nx,ny,nz)   ! or similar
! NEW:
real(8), allocatable, device :: QJ(nx,5,ny,nz)
```
(Exact names may differ per case — search for `(5,nx` in each file.)

---

## Order Within Phase 2

1. `calc_rescale.f90`
2. `set_bc_common.f90`
3. `preprocess.f90`
4. `calc_time_dev.f90`
5. `src/print.f90`
6. `src/main.f90`
7. All 6 `set.f90` files (can do in parallel)
8. **All 6 `mod_globals.f90`** ← compile after this

---

## Verification

```bash
cd 3D_solver/NSTGV && make clean && make
bash calc.sh        # no NaN/Inf; KE decay matches pre-transposition baseline
bash profile.sh     # nsys: compare calc_keep_x_in, calc_Ev4_in kernel times

cd 3D_solver/SBLI  && make clean && make && bash calc.sh
cd 3D_solver/ETGV  && make clean && make && bash calc.sh
```

Expected: 15–30% end-to-end speedup on large periodic cases (NSTGV 513³).
