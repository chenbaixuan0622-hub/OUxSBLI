# Plan: Update CLAUDE.md and README.md

## Context

New untracked files have been added to the repository:
- `2D_solver/OS/` — new 2D oblique shock case (Mach 2, 8° deflection, SLAU, Euler)
- `ouxsbli/tests/test_os.py` — integration test: builds OS, verifies Rankine-Hugoniot
- `ouxsbli/tests/test_evc.py` — grid-convergence test for Euler vortex convection (KEEP/SLAU)
- `ouxsbli/tests/test_corn.py` — integration test for 3D_solver_curv/CORN case
- `ouxsbli/tests/utils/oblique_shock.py` — analytical θ-β-M relations utility
- `tutorials/ouxsbli_bl/` — standalone supersonic boundary layer tutorial case

Neither `CLAUDE.md` nor `README.md` currently documents the 2D solver at all (the 2D solver has cases BL, DSL, EVC, OS, ST under `2D_solver/`). Both files also omit the CORN curvilinear case and the Python test suite details.

## Files to Modify

- `/home/jhatayama/ouxsbli/2d/OUxSBLI/CLAUDE.md`
- `/home/jhatayama/ouxsbli/2d/OUxSBLI/README.md`

---

## CLAUDE.md Changes

### 1. Add "2D Solver" section (after Curvilinear Solver section)

Content to add:

```
### 2D Solver (`2D_solver/`)

A standalone 2D solver sharing the same convective/viscous kernels as the 3D solver.
Source layout:
- `2D_solver/src/` — shared 2D utilities (main, grid, BCs)
- `2D_solver/<CASE>/` — per-case config: `mod_globals.f90`, `set.f90`, `Makefile`, `calc.sh`

Available cases:

| Case | Description |
|------|-------------|
| BL   | Supersonic laminar boundary layer |
| DSL  | Double shear layer |
| EVC  | Euler vortex convection (convergence study) |
| OS   | 2D oblique shock (M=2, θ=8°, SLAU, Euler) |
| ST   | Sod shock tube |

Build and run (same pattern as 3D):
```bash
cd 2D_solver/OS
make clean && make
bash calc.sh
```
```

### 2. Update Curvilinear Solver section

Add CORN to the list of per-case directories:
- Currently only mentions `3D_solver_curv/NACA` by name; add CORN (compression corner, M=2, θ=8°).

### 3. Add "Python Test Suite" section (after Compile Testing section)

```
## Python Test Suite

Integration and convergence tests live in `ouxsbli/tests/`. Run with pytest from the repo root:

```bash
pytest ouxsbli/tests/
```

| Test file | What it checks |
|-----------|----------------|
| `test_etgv.py` | Supersonic Taylor-Green vortex (energy decay) |
| `test_st.py`   | Sod shock tube (exact Riemann solution) |
| `test_evc.py`  | Euler vortex convergence (KEEP 2nd/4th/6th, SLAU 2nd; expected order ≥1.5/3.5) |
| `test_os.py`   | 2D oblique shock — pre/post state vs. Rankine-Hugoniot (tol 2%/5%) |
| `test_corn.py` | 3D_solver_curv/CORN — pressure and density ratios vs. θ-β-M theory (tol 5%) |

Analytical helpers are in `ouxsbli/tests/utils/`:
- `oblique_shock.py` — `beta_from_theta()`, `post_shock_state()` via bisection on θ-β-M
- `sod_exact.py` — exact Riemann solver for Sod tube
- `vtk_reader.py` — VTK output reader
```

### 4. Add "Tutorials" subsection (inside or after Test Cases section)

```
## Tutorials

`tutorials/ouxsbli_bl/` — Supersonic flat-plate boundary layer (M=2, dimensional parameters,
wall-normal grid stretching, Riemann-invariant top BC). Run like any case:

```bash
cd tutorials/ouxsbli_bl
make clean && make
bash calc.sh
```
```

---

## README.md Changes

### 1. Update first paragraph / feature list

Add "2D solver" to the one-liner description of solvers available.

### 2. Add 2D Oblique Shock to Validations section

Add a brief entry (similar to the SBLI entry) noting:
- 2D oblique shock: M=2, θ=8°, Rankine-Hugoniot verified

### 3. (Minor) Mention CORN in curvilinear section if a visualization is available; otherwise skip.

---

## Verification

After editing:
1. Read both files to confirm formatting is correct Markdown.
2. Check that no existing content was accidentally removed.
3. The CICD script (`test_cicd.sh`) does not need changes — it only covers compilation targets, not test content.
