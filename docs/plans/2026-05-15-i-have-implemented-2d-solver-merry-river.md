# Bug Report: 2D_solver

## Context

The user implemented a 2D CFD solver under `2D_solver/` with four test cases (EVC, DSL, ST, BL). The code was adapted from the existing 3D solver. Several compilation-blocking and runtime bugs were introduced during the adaptation.

---

## Confirmed Bugs

### Bug 1 — CRITICAL (all cases): `calc_time_dev.f90:7–9` imports non-existent symbols

**File:** [2D_solver/src/calc_time_dev.f90](2D_solver/src/calc_time_dev.f90#L7-L9)

```fortran
use mod_globals, only : id_visc, nt, np, nre2, rerank, &
& blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
& blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
```

The 2D solver has no z-direction flux, so `blocksG`, `threadsG`, `blocksGv`, `threadsGv` are absent from every `mod_globals.f90`. The variables `nre2` and `rerank` do not exist in any 2D case's `mod_globals` either (3D-only leftovers). This prevents **all four cases** from compiling.

**Fix:** Remove `nre2, rerank, blocksG, threadsG, blocksGv, threadsGv` from the `use` list.

---

### Bug 2 — CRITICAL (EVC): `EVC/set.f90:2` imports undefined `dtn`

**File:** [2D_solver/EVC/set.f90](2D_solver/EVC/set.f90#L2)

```fortran
use mod_globals, only : id_accuracy, nx, ny, Lx, Ly, gamma, R, dtn, theta
```

`EVC/mod_globals.f90` defines `dt`, not `dtn`. The variable `dtn` does not exist in the module; Fortran will error on the `use` statement even though `dtn` is never referenced in the subroutine bodies.

**Fix:** Remove `dtn` from the `only:` list on line 2 of `EVC/set.f90`.

---

### Bug 3 — CRITICAL (BL): `BL/set.f90:4–5` has three import errors

**File:** [2D_solver/BL/set.f90](2D_solver/BL/set.f90#L4-L5)

```fortran
use mod_globals, only : id_rescale, ny1, nre2, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy, rf, Taw
use mod_constant, only : Cp, gamma_1, over_gamma, over_gamma_1
```

Three problems:
- `ny1` — not defined in `BL/mod_globals.f90` (3D leftover)
- `nre2` — not defined in `BL/mod_globals.f90` (3D leftover)
- `Cp` — appears in both `use mod_globals` (line 4) AND `use mod_constant` (line 5). `Cp` is not in `BL/mod_globals`, so the mod_globals import fails; and even if it were, importing the same name from two modules is ambiguous and will error.

**Fix:** Remove `ny1, nre2, Cp` from the `mod_globals` use statement on line 4. `Cp` is correctly available from `mod_constant` on line 5.

---

### Bug 4 — RUNTIME CRASH (EVC, DSL, ST, BL): Missing `$` in CUF kernel directives for corner cells

**File:** [2D_solver/src/set_bc_common.f90](2D_solver/src/set_bc_common.f90)

Three occurrences — lines 91, 123, and 171:

```fortran
!cuf kernel do(1)<<<*,*>>>   ← missing $; this is a plain comment, not a CUF directive
do l = 1, 4
  Q(1,l,1) = Q(nx-1,l,ny-1)
  ...
enddo
```

Because `!cuf` is not `!$cuf`, the loop runs on the **CPU** while `Q` is a `device` array. Accessing device memory from the CPU causes an **illegal memory access / runtime crash**. The surrounding face loops (lines 79–90, 107–122, 151–170) are correct (`!$cuf`), so only the corner-cell assignments are affected.

**Fix:** Change `!cuf` → `!$cuf` at lines 91, 123, and 171.

---

### Bug 5 — WRONG RESULTS (ST): Off-by-one in shock tube initial condition

**File:** [2D_solver/ST/set.f90](2D_solver/ST/set.f90#L41)

```fortran
do i = 1, nx / 2        ! sets left state for i = 1 … 2048
  Q(i,*,j) = left_state
enddo
do i = nx / 2, nx       ! BUG: starts at 2048, overwrites the left-state cell with right state
  Q(i,*,j) = right_state
enddo
```

For `nx = 4097`, `nx/2 = 2048`. Cell 2048 should be the last left-state cell, but the second loop overwrites it with the right state. The diaphragm location shifts by one cell.

**Fix:** Change `do i = nx / 2, nx` → `do i = nx / 2 + 1, nx` on line 41.

---

## Summary Table

| # | File | Lines | Severity | Impact |
|---|------|-------|----------|--------|
| 1 | `src/calc_time_dev.f90` | 7–9 | CRITICAL | All 4 cases fail to compile |
| 2 | `EVC/set.f90` | 2 | CRITICAL | EVC fails to compile |
| 3 | `BL/set.f90` | 4–5 | CRITICAL | BL fails to compile |
| 4 | `src/set_bc_common.f90` | 91, 123, 171 | RUNTIME CRASH | Corner ghost cells access device memory from CPU |
| 5 | `ST/set.f90` | 41 | WRONG RESULTS | Shock tube diaphragm off by one cell |

## Verification

After fixes, build each case and check for clean compilation:
```bash
cd 2D_solver/EVC && make clean && make
cd 2D_solver/DSL && make clean && make
cd 2D_solver/ST  && make clean && make
cd 2D_solver/BL  && make clean && make
```
Then run ST and confirm the shock front position at `t=0+ε` matches the analytic Sod solution contact surface at `x = nx/2 + 1` (not `nx/2`).
