# Plan: Add pytest test suite for OUxSBLI validation cases

## Context
The user wants integration tests in `ouxsbli/tests/` covering five physical validation cases. Tests build and run real simulations via the `Case` class, then compare output against analytical solutions. A new 2D oblique shock case (`2D_solver/OS`) must also be created.

---

## Files to create

### Python utilities
| File | Purpose |
|------|---------|
| `ouxsbli/tests/utils/__init__.py` | empty namespace |
| `ouxsbli/tests/utils/vtk_reader.py` | binary VTK RectilinearGrid reader |
| `ouxsbli/tests/utils/sod_exact.py` | exact Sod shock-tube solver |
| `ouxsbli/tests/utils/oblique_shock.py` | θ-β-M analytical relations |
| `ouxsbli/tests/conftest.py` | `tmp_path`-based fixtures, `integration` mark |

### Test files
| File | Case | Assertion |
|------|------|-----------|
| `ouxsbli/tests/test_evc.py` | EVC | convergence order ≥ expected − 0.5 |
| `ouxsbli/tests/test_ivst.py` | IVST | L1(ρ,p) < 5 % vs Sod exact |
| `ouxsbli/tests/test_os.py` | OS (new) | post-shock ρ₂,p₂ within 5 % of analytical |
| `ouxsbli/tests/test_etgv.py` | ETGV | |Δs/s₀| < 1e-3, |ΔKE/KE₀| < 1e-3 |
| `ouxsbli/tests/test_corn.py` | CORN | post-shock p₂/p₁ within 5 % of analytical |

### New Fortran case
| File | Notes |
|------|-------|
| `2D_solver/OS/mod_globals.f90` | M=2, θ=8°, β≈37.2°, Euler, SLAU 2nd |
| `2D_solver/OS/set.f90` | diagonal-shock IC; left inflow / right extrapolation / periodic y |
| `2D_solver/OS/Makefile` | copy from `2D_solver/EVC`, set case dir |
| `2D_solver/OS/calc.sh` | `mpirun -n 2 ./a.out` |

---

## VTK reader (`vtk_reader.py`)

`print.f90` opens files with `form="unformatted", access="stream", convert="Little_ENDIAN"`. The binary layout after the XML header is:

```
<AppendedData encoding="raw">  _[int32 n_bytes][float32*ni][int32][float32*nj]
  [int32][float32*nk][int32][float32*ni*nj*nk rho][int32][float32*ni*nj*nk p]
  [int32][float32*3*ni*nj*nk velocity]
```

Grid dimensions come from `WholeExtent="0 NI 0 NJ 0 NK"` in the XML header (ni=NI+1, etc.). Data layout: i (x) fastest, then j, then k — reshape as `(nk, nj, ni)`.

For 2D cases: single file `data/Q%05d.vtr` written by rank 1 with rank 0's data (x ∈ [0, Lx/2]).  
For 3D with nranks < 4: same convention (`data/Q%05d.vtr`, rank 0's x-half).

---

## Analytical helpers

**`sod_exact.py`** — Newton-Raphson on the pressure function f(p*)=f_L+f_R+ΔU=0, then five-region reconstruction. Returns ρ(x,t), u(x,t), p(x,t) for any array of positions.

**`oblique_shock.py`** — Given M₁, θ (deflection), γ, solve θ-β-M relation with bisection to get β, then apply Rankine-Hugoniot on the normal Mach M_n1=M₁ sin β:
- ρ₂/ρ₁ = (γ+1)M_n1² / ((γ-1)M_n1²+2)
- p₂/p₁ = 1 + 2γ(M_n1²−1)/(γ+1)
- M_n2 = √(((γ-1)M_n1²+2)/(2γM_n1²−(γ-1)))
- M₂ = M_n2 / sin(β−θ)

Post-shock velocity: u_x2 = M₂·a₂·cos(θ), u_y2 = −M₂·a₂·sin(θ).

---

## Test details

### 1 — EVC grid convergence (`test_evc.py`)
Source: `2D_solver/EVC`. Build 4 schemes × 3 grid levels = 12 runs.

| Scheme | `id_scheme` | `id_accuracy` | Expected order |
|--------|-------------|---------------|----------------|
| KEEP2 | keep | 2 | 2 |
| KEEP4 | keep | 4 | 4 |
| KEEP6 | keep | 6 | 4 (capped by RK4) |
| SLAU2 | slau | 2 | 2 |

Grid levels: nx=ny ∈ [18, 34, 66] (small to keep CI fast). At t=T=Lx/u0 (one vortex period), compare ρ with analytical IC via L2 error. Compute log-ratio of errors → check order ≥ expected − 0.5.

Patch per build: `nx`, `ny`, `id_scheme`, `id_accuracy`. Keep `nt` small (patch to 1) to run only to t=T.

### 2 — Inviscid shock tube (`test_ivst.py`)
Source: `3D_solver/IVST`. One build, `nranks=2`. Reads `data/Q00001.vtr` (rank 0's x∈[0,0.5] portion, Fortran default step naming = last output index).

Sod initial state: ρ_L=1, p_L=1; ρ_R=0.125, p_R=0.1; t_final=0.2. Compare slice at j=4 (middle y) against `sod_exact`. Assert L1(ρ) < 5 % over the left half.

### 3 — 2D oblique shock (`test_os.py`)
Source: `2D_solver/OS` (new case, created as part of this plan).

**OS case physics** (M=2, θ=8°, γ=1.4):
- Solve β≈37.2° from θ-β-M; M_n1≈1.208
- ρ₂/ρ₁≈1.355, p₂/p₁≈1.536, M₂≈1.715

**IC** (set.f90): shock line `x·sin β − y·cos β = c` (c chosen to center shock). Cells satisfying the condition → post-shock; otherwise → pre-shock. BCs: left fixed inflow (pre-shock), right extrapolation, top/bottom periodic.

**mod_globals.f90** key parameters:
```fortran
integer(2), parameter :: id_visc    = 1    ! Euler
real(2),    parameter :: id_scheme  = 0    ! SLAU
integer(2), parameter :: id_accuracy= 0   ! 2nd order
logical, parameter :: id_bc_x = .true.
logical, parameter :: id_bc_y = .false.  ! periodic y
real(8), parameter :: Lx = 4.0d0, Ly = 2.0d0
integer, parameter :: nx = 257, ny = 129
```

Test reads final VTK, splits cells by shock-line condition, asserts that the post-shock region mean matches analytical ρ₂, p₂ within 5 %; pre-shock region mean matches ρ₁, p₁ within 2 %.

### 4 — ETGV entropy & KE preservation (`test_etgv.py`)
Source: `3D_solver/ETGV`. Patch `nt=5, np=20` (100 total steps, ~1 s wall-clock on GPU). Reads text files:
- `data/entropy.d`: columns `[t, (s0−s)/s0]`
- `data/kinetic_energy.d`: columns `[t, ke, ke/ke0]`

Assertions on every row:
- `|(s0−s)/s0| < 1e-3`  (entropy preserved)
- `|ke/ke0 − 1| < 1e-3`  (KE preserved)

### 5 — CORN compression corner (`test_corn.py`)
Source: `3D_solver_curv/CORN`. One build, `nranks=2`. VTK x-coords are computational ξ=1..nx; the wall corner is at ξ≈i_corner where the physical x=x_corner=0.5, corresponding to i_corner = round(0.5·(nx−2)/Lx)+2 ≈ round(0.5·190/2)+2 = 49.

Extract post-shock region: i ∈ [i_corner+20, nx-10], j ∈ [1, 5] (near lower wall, downstream of corner). Assert mean p/(p·γ⁻¹) / (p₁/(p₁·γ⁻¹)) ≈ p₂/p₁=1.536 within 5 % (same analytical from oblique_shock.py with M=2, θ=8°).

---

## conftest.py design

```python
import pytest

def pytest_configure(config):
    config.addinivalue_line("markers", "integration: requires GPU + HPC SDK")

@pytest.fixture(scope="session")
def repo_root():
    return pathlib.Path(__file__).resolve().parents[2]

# Per-test workdir under /tmp so builds don't pollute the repo
@pytest.fixture
def workdir(tmp_path):
    return tmp_path / "sim"
```

All integration tests are decorated `@pytest.mark.integration`.

---

## OS Makefile
Copy `2D_solver/EVC/Makefile` verbatim; only the `vpath` line changes (points to `2D_solver/src`). The Case class rewrites vpaths to absolute paths automatically, so no manual edits are needed beyond copying.

---

## Verification
Run after implementation:
```bash
cd /home/jhatayama/ouxsbli/2d/OUxSBLI
pytest ouxsbli/tests/ -v -m integration
```
Expected: all 5 test groups pass. For a quick smoke test without GPU:
```bash
pytest ouxsbli/tests/test_patcher.py -v   # existing unit tests, no GPU needed
```
