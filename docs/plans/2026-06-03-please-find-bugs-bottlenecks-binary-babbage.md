# Refactoring Plan — OUxSBLI (Whole Codebase)

## Context

Full-codebase refactor driven by four goals: readability/naming, reduce duplication, modularity,
and performance. All changes preserve existing behaviour and GPU execution paths. Items are
ordered from safest (pure Python) to highest risk (GPU kernel code).

Round 1 and Round 2 bug/dead-code fixes are already applied. This plan covers structural
improvement only.

---

## Python Test Suite

### P1 — Create `ouxsbli/tests/conftest.py`

Every test file repeats the same three patterns. Extract them once into pytest fixtures:

```python
# Build + run a solver case
@pytest.fixture
def run_case(tmp_path):
    def _runner(source, workdir, **params):
        case = Case(source=source, workdir=workdir, **params)
        case.build()
        case.run(nranks=params.pop("nranks", 2))
        return pathlib.Path(workdir) / "data"
    return _runner

# Relative-error assertion (used in test_os.py and test_corn.py)
def assert_close_relative(computed, reference, rtol, label=""):
    rel_err = abs(computed - reference) / reference
    assert rel_err < rtol, \
        f"{label}: rel_err={rel_err:.3%}, computed={computed:.4g}, ref={reference:.4g}"
```

Update each test to import and use these instead of inline boilerplate.

### P2 — Refactor `ouxsbli/tests/utils/vtk_reader.py`

Two independent cleanups:

1. **Merge `extract_number` / `extract_number_vts`** — both parse `_(\d+)\.(vtr|vts)` from a
   filename; collapse to one function that accepts the extension.

2. **Unify `getQ` vts/vtr reshape** — lines 74–88 (vts path) and 89–107 (vtr path) both call
   the same reshape sequence. Extract:
   ```python
   def _reshape_q(rho, vel, p, Nz, Ny, Nx):
       rho = rho.reshape((Nz, Ny, Nx))
       u, v, w = [vel.reshape((Nz,Ny,Nx,3))[:,:,:,i] for i in range(3)]
       return rho, u, v, w, p.reshape((Nz,Ny,Nx))
   ```
   Call it from both branches.

### P3 — `test_evc.py`: convert inner grid-size loop to `pytest.mark.parametrize`

Currently loops `for nx in GRIDS` inside a single test. Convert so each (scheme, order, nx)
triple is a separate pytest item with its own pass/fail status.

---

## Core Fortran Utilities (host-side code — no GPU execution path changes)

### F1 — `3D_solver/src/calc_para.f90`: extract pack/unpack slab helpers

The three `flatten` variants and three `reconstruct` variants share identical loop bodies
differing only in the start x-index. The core index expression:

```fortran
nj*ni*5*(k-1) + ni*5*(j-1) + ni*(l-1) + i_offset
```

appears ~10 times verbatim. Extract two private subroutines:

```fortran
subroutine pack_slab(Q, Q1d, i_start, overlap, ny, nz)
subroutine unpack_slab(Q1d, Q, i_start, overlap, ny, nz)
```

Replace `flatten`, `flatten_left`, `flatten_right` with `pack_slab` calls with appropriate
`i_start`. Do the same for the `reconstruct` family.
Keep `reconstruct_sbli_inlet` as-is (distinct signature with `ny1` / `ny2` rescale logic).

### F2 — `src/set_coordinate.f90`: extract ghost-fill helper and fix naming

The six `set_grid_cyclic{2,4,6}_{2D,3D}` routines each repeat the same boundary ghost-cell
fill pattern for x (and y/z). Extract as a private `fill_ghost_x(arr, nx, ny, nz)` helper
(and `_y`, `_z` variants) and call from all six routines. Do NOT merge the six routines
themselves — different stencil widths and 2D/3D loop structures are worth keeping separate.

Also rename: `set_block2` → `set_block_2D`, `set_block3` → `set_block_3D`.

### F3 — `src/calc_physical_quantities.f90`: extract Sutherland formula

The expression

```fortran
mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
```

appears identically in both `calc_quantities_T_2D` and `calc_quantities_T_3D`. Extract as a
module-private pure function:

```fortran
pure function sutherland(temp) result(mu)
  real(8), intent(in), value :: temp
  real(8) :: mu
  mu = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
end function
```

### F4 — `src/calc_muscl.f90`: rename stencil-difference locals

The local variables `d1, d2, d3` (in `delta4`) and `d1..d5` (in `delta6`) encode stencil
positions but the names are opaque. Rename:

| Old | New | Meaning |
|-----|-----|---------|
| `d1` | `dm1` | left-of-left difference |
| `d2` | `d0`  | left difference |
| `d3` | `dp1` | centred difference |
| `d4` | `dp2` | right difference |
| `d5` | `dp3` | right-of-right difference |

Also rename the private MUSCL subroutines to consistent snake_case:
`MUSCL3rdnonTVD` → `muscl_3rd_none`, `MUSCL3rdMinmod` → `muscl_3rd_minmod`,
`MUSCL3rdThreshold` → `muscl_3rd_threshold` (and 4th-order equivalents). These are module-private;
no callers outside the file.

### F5 — `src/print.f90`: naming and constant cleanup

- Rename `rho1d`, `p1d`, `v1d` → `rho_flat`, `p_flat`, `vel_flat` (the `1d` suffix misleads;
  these are flattened 3D arrays used for I/O).
- Add `integer, parameter :: IO_UNIT_VTK = 10` and replace all hardcoded `10` unit literals.

---

## GPU Kernel Files (naming and comment changes only — no logic changes)

### G1 — `3D_solver/src/calc_les.f90`: rename cryptic locals

| Old | New | Rationale |
|-----|-----|-----------|
| `uh, vh, wh` | `u_test, v_test, w_test` | test-filtered (coarse) velocity |
| `vord` | `vorticity_diff` | resolved minus filtered vorticity |
| `qc2` | `subgrid_ke` | subgrid kinetic energy correction |

`qc2` is also used in `calc_flux_base.f90.fypp` and `preprocess.f90.fypp`, so rename
consistently across all three files.

### G2 — `3D_solver/src/calc_slau_kernel.f90.fypp`: document sensor precision choice

`fdx`, `fdy`, `fdz` are declared `real(sp)`. This is intentional — single precision reduces
register pressure and bandwidth on the sensor interpolation. Add a one-line comment at each
`real(sp) fdx` declaration in `calc_slau_x`, `calc_slau_y`, `calc_slau_z` (and the `_internal`
variant) to record the deliberate choice.

### G3 — `3D_solver/src/calc_visc4_internal.f90.fypp`: remove empty else branches

The interior-only kernels have:

```fortran
if (3 <= i .and. i <= nx-3 .and. ...) then
  ! high-order computation
else
  ! (nothing — interior variants have no low-order fallback)
endif
```

Remove the empty `else` clause in all three direction variants. Reduces dead code and removes
the theoretical register pressure from a dead branch.

---

## Files to Modify

| File | Change |
|------|--------|
| `ouxsbli/tests/conftest.py` | **Create** — run_case fixture, assert_close_relative |
| `ouxsbli/tests/utils/vtk_reader.py` | Merge extract_number; unify getQ reshape |
| `ouxsbli/tests/test_evc.py` | Parametrize grid-convergence loop |
| `ouxsbli/tests/test_os.py` | Use assert_close_relative from conftest |
| `ouxsbli/tests/test_corn.py` | Use assert_close_relative from conftest |
| `3D_solver/src/calc_para.f90` | Extract pack_slab / unpack_slab |
| `src/set_coordinate.f90` | Extract fill_ghost helpers; rename set_block_2D/3D |
| `src/calc_physical_quantities.f90` | Extract sutherland() |
| `src/calc_muscl.f90` | Rename stencil diffs and MUSCL subroutines |
| `src/print.f90` | Rename flat arrays; IO_UNIT_VTK constant |
| `3D_solver/src/calc_les.f90` | Rename uh/vh/wh, vord, qc2 |
| `3D_solver/src/calc_flux_base.f90.fypp` | Rename qc2 → subgrid_ke |
| `3D_solver/src/preprocess.f90.fypp` | Rename qc2 → subgrid_ke |
| `3D_solver/src/calc_slau_kernel.f90.fypp` | Add sp-precision comment |
| `3D_solver/src/calc_slau_kernel_internal.f90.fypp` | Add sp-precision comment |
| `3D_solver/src/calc_visc4_internal.f90.fypp` | Remove empty else branches |

---

## Out of Scope (explicitly excluded)

- **Generic kernel wrapper** — restructuring all `calc_*_kernel*.f90.fypp` files simultaneously
  carries high risk of GPU register-allocation regressions and subtle correctness bugs.
- **Q array layout** (`Q(nx,5,ny,nz)` → `Q(5,nx,ny,nz)`) — affects every file; requires a
  full benchmark pass to verify no performance regression.
- **`calc_time_dev.f90.fypp` RK variant unification** — four 100-line FYPP branches; the COMMZ
  overlapped-comms logic in particular must stay separate.
- **`preprocess.f90.fypp` allocation splitting** — changing device allocation sequencing risks
  breaking CUDA context setup order.

---

## Verification

| Change group | Command |
|---|---|
| Python (P1–P3) | `pytest ouxsbli/tests/ -v` |
| Fortran utilities (F1–F5) | `cd 3D_solver/NSTGV && cmake -B build && cmake --build build -j` |
| 2D solver sanity | `cd 2D_solver/ST && cmake -B build && cmake --build build -j` |
| calc_les / subgrid_ke rename (G1) | `cd 3D_solver/TBL && cmake -B build && cmake --build build -j` |
| SLAU sensor comment (G2) | `cd 3D_solver/SBLI && cmake -B build && cmake --build build -j` |
| calc_visc4_internal (G3) | `cd 3D_solver/TBL && cmake -B build && cmake --build build -j` |
