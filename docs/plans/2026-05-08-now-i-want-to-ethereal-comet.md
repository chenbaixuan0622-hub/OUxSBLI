# Plan: Comprehensive CICD Directories for 3D_solver and 3D_solver_curv

## Context

The CI/CD suite should compile-test every supported (id_scheme × id_accuracy) combination for both solvers. Per CLAUDE.md, Roe is unused; the three active schemes are KEEP, SLAU, Hybrid. Accuracy has three kinds: 2nd (kind=2), 4th (kind=4), 6th (kind=8). That gives **9 combinations per solver = 18 total CICD build targets**.

`3D_solver/CICD/` already has 9 `CICD_*` subdirectories (untracked), but they contain two bugs:
1. **Wrong vpath**: `../src:../../src` was copied from the parent `CICD/` Makefile but does not reach `OUxSBLI/src/` (main.f90, mod_constant.f90, etc.) from one level deeper. Correct path: `../../src:../../../src`.
2. **Wrong kind parameters**: Only CICD_KEEP2 (id_scheme), CICD_SLAU2/4/6 (id_scheme), and CICD_SLAU6/CICD_KEEP6/CICD_Hybrid6 (id_accuracy) are partially correct. Summary of what each file currently has vs. needs:

| Directory  | id_scheme now   | id_scheme needed | id_accuracy now | id_accuracy needed |
|------------|-----------------|------------------|-----------------|-------------------|
| CICD_KEEP2 | `integer(2)` ✓  | `integer(2)`     | `kind=8` ✗      | `kind=2`          |
| CICD_KEEP4 | `real(2)` ✗     | `integer(2)`     | `kind=8` ✗      | `kind=4`          |
| CICD_KEEP6 | `real(2)` ✗     | `integer(2)`     | `kind=8` ✓      | `kind=8`          |
| CICD_SLAU2 | `real(2)` ✓     | `real(2)`        | `kind=8` ✗      | `kind=2`          |
| CICD_SLAU4 | `real(2)` ✓     | `real(2)`        | `kind=8` ✗      | `kind=4`          |
| CICD_SLAU6 | `real(2)` ✓     | `real(2)`        | `kind=8` ✓      | `kind=8`          |
| CICD_Hybrid2 | `real(2)` ✗   | `real(8)`        | `kind=8` ✗      | `kind=2`          |
| CICD_Hybrid4 | `real(2)` ✗   | `real(8)`        | `kind=8` ✗      | `kind=4`          |
| CICD_Hybrid6 | `real(2)` ✗   | `real(8)`        | `kind=8` ✓      | `kind=8`          |

`3D_solver_curv/CICD/` has only one flat case (Hybrid + NS + 2nd order). All 9 `CICD_*` subdirectories must be created from scratch, adapting the curvilinear `mod_globals.f90` for each combination.

The GitHub Actions workflow currently builds only `3D_solver/NSTGV` and `3D_solver_curv/NACA`. It must be updated to build all 18 CICD_* targets using a matrix strategy.

---

## Critical Files

- `3D_solver/CICD/CICD_*/Makefile` — 9 files, all need vpath fix
- `3D_solver/CICD/CICD_*/mod_globals.f90` — 9 files, 8 need id_scheme/id_accuracy fix
- `3D_solver/CICD/CICD_*/set.f90` — verify all 9 have it; copy from CICD_KEEP2 if any are missing
- `3D_solver_curv/CICD/CICD_*/Makefile` — 9 new files
- `3D_solver_curv/CICD/CICD_*/mod_globals.f90` — 9 new files
- `3D_solver_curv/CICD/CICD_*/set.f90` — 9 new files (copy of `3D_solver_curv/CICD/set.f90`)
- `.github/workflows/makefile.yml` — update to matrix build

---

## Implementation Steps

### Step 1 — Fix Cartesian CICD_* directories (3D_solver/CICD/)

For each of the 9 directories, update two files:

**Makefile**: change the vpath line from  
`vpath %f90 ../src:../../src`  
to  
`vpath %f90 ../../src:../../../src`

The rest of the Makefile (OBJ list, compile flags, dependency rules, `a.out` link rule) stays identical across all 9.

**mod_globals.f90**: update only the two kind-tagged parameter lines. All other parameters (physics, grid, BCs, threads, dt formula) stay fixed. The table of changes:

| Directory    | Line to change                       | Replace with                          |
|--------------|--------------------------------------|---------------------------------------|
| CICD_KEEP2   | `integer(kind=8), parameter :: id_accuracy = 0` | `integer(kind=2), parameter :: id_accuracy = 0` |
| CICD_KEEP4   | `real(2), parameter :: id_scheme = 0` | `integer(2), parameter :: id_scheme = 0` |
| CICD_KEEP4   | `integer(kind=8), parameter :: id_accuracy = 0` | `integer(kind=4), parameter :: id_accuracy = 0` |
| CICD_KEEP6   | `real(2), parameter :: id_scheme = 0` | `integer(2), parameter :: id_scheme = 0` |
| CICD_SLAU2   | `integer(kind=8), parameter :: id_accuracy = 0` | `integer(kind=2), parameter :: id_accuracy = 0` |
| CICD_SLAU4   | `integer(kind=8), parameter :: id_accuracy = 0` | `integer(kind=4), parameter :: id_accuracy = 0` |
| CICD_Hybrid2 | `real(2), parameter :: id_scheme = 0` | `real(8), parameter :: id_scheme = 0` |
| CICD_Hybrid2 | `integer(kind=8), parameter :: id_accuracy = 0` | `integer(kind=2), parameter :: id_accuracy = 0` |
| CICD_Hybrid4 | `real(2), parameter :: id_scheme = 0` | `real(8), parameter :: id_scheme = 0` |
| CICD_Hybrid4 | `integer(kind=8), parameter :: id_accuracy = 0` | `integer(kind=4), parameter :: id_accuracy = 0` |
| CICD_Hybrid6 | `real(2), parameter :: id_scheme = 0` | `real(8), parameter :: id_scheme = 0` |

**set.f90**: CICD_KEEP2 already has one. Verify the other 8 dirs have set.f90; if any are missing, copy CICD_KEEP2/set.f90 verbatim (it is scheme-agnostic: periodic TGV init + all-periodic BCs).

### Step 2 — Create curvilinear CICD_* directories (3D_solver_curv/CICD/)

Create 9 subdirectories. Each contains three files:

**Makefile** (identical for all 9):  
- Copy `3D_solver_curv/CICD/Makefile` and change:  
  `vpath %f90 ../src:../../src:../../3D_solver/src`  
  → `vpath %f90 ../../src:../../../src:../../../3D_solver/src`  
- Everything else (OBJ list, compile flags, dependency rules, `a.out` rule) is unchanged.

**set.f90** (identical for all 9):  
- Copy `3D_solver_curv/CICD/set.f90` verbatim (NACA O-grid IC + BCs, scheme-agnostic).

**mod_globals.f90** (9 variants):  
- Base: `3D_solver_curv/CICD/mod_globals.f90` (Hybrid, NS, kind=2 accuracy)  
- Change only `id_scheme` and `id_accuracy` per directory:

| Directory       | id_scheme                     | id_accuracy             |
|-----------------|-------------------------------|-------------------------|
| CICD_KEEP2      | `integer(2), parameter`       | `integer(kind=2), parameter` |
| CICD_KEEP4      | `integer(2), parameter`       | `integer(kind=4), parameter` |
| CICD_KEEP6      | `integer(2), parameter`       | `integer(kind=8), parameter` |
| CICD_SLAU2      | `real(2), parameter`          | `integer(kind=2), parameter` |
| CICD_SLAU4      | `real(2), parameter`          | `integer(kind=4), parameter` |
| CICD_SLAU6      | `real(2), parameter`          | `integer(kind=8), parameter` |
| CICD_Hybrid2    | `real(8), parameter`          | `integer(kind=2), parameter` |
| CICD_Hybrid4    | `real(8), parameter`          | `integer(kind=4), parameter` |
| CICD_Hybrid6    | `real(8), parameter`          | `integer(kind=8), parameter` |

All other parameters stay as in the curvilinear CICD base (NACA 0012 geometry, NS viscosity `integer(4)`, non-TVD `integer(2)`, TVD-RK3 `integer(2)`, Ma_inf=0.8, dt=2e-9, nt=1000, np=100).

### Step 3 — Update GitHub Actions workflow

File: `.github/workflows/makefile.yml`

Replace the two sequential build steps with two **matrix jobs** (one for Cartesian, one for curvilinear):

```yaml
jobs:
  build-3d-solver:
    runs-on: ubuntu-latest
    container:
      image: nvcr.io/nvidia/nvhpc:25.11-devel-cuda13.0-ubuntu24.04
    strategy:
      matrix:
        case: [CICD_KEEP2, CICD_KEEP4, CICD_KEEP6,
               CICD_SLAU2, CICD_SLAU4, CICD_SLAU6,
               CICD_Hybrid2, CICD_Hybrid4, CICD_Hybrid6]
    steps:
      - uses: actions/checkout@v4
      - name: Setup MPI path
        run: |
          MPI_DIR=$(find /opt/nvidia/hpc_sdk -type d -path "*/comm_libs/mpi/bin" | head -n 1)
          echo "$MPI_DIR" >> $GITHUB_PATH
      - name: Build ${{ matrix.case }}
        run: |
          cd 3D_solver/CICD/${{ matrix.case }}
          make

  build-3d-solver-curv:
    runs-on: ubuntu-latest
    container:
      image: nvcr.io/nvidia/nvhpc:25.11-devel-cuda13.0-ubuntu24.04
    strategy:
      matrix:
        case: [CICD_KEEP2, CICD_KEEP4, CICD_KEEP6,
               CICD_SLAU2, CICD_SLAU4, CICD_SLAU6,
               CICD_Hybrid2, CICD_Hybrid4, CICD_Hybrid6]
    steps:
      - uses: actions/checkout@v4
      - name: Setup MPI path
        run: |
          MPI_DIR=$(find /opt/nvidia/hpc_sdk -type d -path "*/comm_libs/mpi/bin" | head -n 1)
          echo "$MPI_DIR" >> $GITHUB_PATH
      - name: Build ${{ matrix.case }}
        run: |
          cd 3D_solver_curv/CICD/${{ matrix.case }}
          make
```

---

## Verification

After implementation, verify locally (requires HPC SDK):
```bash
# Cartesian
for d in CICD_KEEP2 CICD_KEEP4 CICD_KEEP6 CICD_SLAU2 CICD_SLAU4 CICD_SLAU6 CICD_Hybrid2 CICD_Hybrid4 CICD_Hybrid6; do
  (cd 3D_solver/CICD/$d && make clean && make) || echo "FAILED: $d"
done

# Curvilinear
for d in CICD_KEEP2 CICD_KEEP4 CICD_KEEP6 CICD_SLAU2 CICD_SLAU4 CICD_SLAU6 CICD_Hybrid2 CICD_Hybrid4 CICD_Hybrid6; do
  (cd 3D_solver_curv/CICD/$d && make clean && make) || echo "FAILED: $d"
done
```

CI verification: push to `stable` branch → GitHub Actions should show 18 green matrix jobs (9 + 9).
