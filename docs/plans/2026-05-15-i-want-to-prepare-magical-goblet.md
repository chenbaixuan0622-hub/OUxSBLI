# Python API for OUxSBLI CFD Code — Approach B (Revised)

## Context

OUxSBLI is a CUDA Fortran + MPI CFD solver. **All parameters are Fortran compile-time constants** in `mod_globals.f90`. There is no runtime config file.

**Key constraint from user:** Original case directories (and their `mod_globals.f90`) must **not** be modified. Instead the Python API copies source cases into new work directories and edits only the copies.

Any Python API must:
1. Copy a source case directory to a new work directory
2. Modify `mod_globals.f90` **in the copy** via targeted text substitution
3. Recompile (`make clean && make`) inside the work directory
4. Launch (`mpirun -n N a.out`) inside the work directory
5. Parse outputs (VTK `.vtr`, text monitoring files)

Key dispatch mechanism: scheme/physics selection uses Fortran **type kind** (not value), e.g. `integer(2)` = KEEP, `real(2)` = SLAU.

---

## Directory Layout (new files only)

```
ouxsbli/                         # new Python package (at repo root)
├── __init__.py                  # public API: Case, Sweep
├── case.py                      # Case class — configure, build, run
├── patcher.py                   # mod_globals.f90 text patcher
├── io/
│   ├── __init__.py
│   ├── vtk.py                   # .vtr reader → numpy arrays
│   └── monitor.py               # kinetic_energy.d / entropy.d → pandas
└── tests/
    ├── test_patcher.py          # unit tests for text patcher
    └── test_io.py               # unit tests for VTK reader
```

No existing files are modified.

---

## How `patcher.py` Works

The original `mod_globals.f90` is the source of truth. The patcher does targeted line-by-line substitution of named parameter declarations. It does **not** use regex on the full file — it locates the line containing the parameter by name, then rewrites only that line.

```python
# ouxsbli/patcher.py

KIND_MAP = {
    "visc":     {"euler": "integer(2)", "ns": "integer(4)", "les": "integer(8)"},
    "scheme":   {"keep": "integer(2)", "slau": "real(2)", "roe": "real(4)", "hybrid": "real(8)"},
    "accuracy": {2: "integer(2)", 4: "integer(4)", 6: "integer(8)"},
    "tvd":      {"none": "integer(2)", "minmod": "integer(4)", "muscl4": "integer(8)"},
    "rk":       {"tvd_rk3": "integer(2)", "rk4": "integer(4)"},
    "rescale":  {False: "integer(2)", True: "integer(4)"},
}

def patch(src_text: str, changes: dict) -> str:
    """Apply parameter changes to mod_globals.f90 text.

    changes = {
        "Re": 1600.0,          # scalar value replacement
        "nx": 65,              # integer replacement
        "scheme": "slau",      # kind-based dispatch
        "visc": "ns",
        "accuracy": 4,
    }
    """
    lines = src_text.splitlines()
    for i, line in enumerate(lines):
        for param, value in changes.items():
            if _is_param_line(line, param):
                lines[i] = _rewrite_line(line, param, value)
    return "\n".join(lines) + "\n"

def _is_param_line(line: str, param: str) -> bool:
    """True if this line declares the named parameter."""
    return (f":: {param}" in line or f":: {param} " in line) and "parameter" in line

def _rewrite_line(line: str, param: str, value) -> str:
    """Rewrite one parameter declaration line."""
    if param in ("scheme", "visc", "accuracy", "tvd", "rk", "rescale"):
        kind = KIND_MAP[param][value]
        # Replace type-kind prefix: e.g. integer(4) → real(2)
        import re
        return re.sub(r"\b(integer|real)\(\d\)", kind, line)
    else:
        # Scalar replacement: find = <value> and replace
        import re
        if isinstance(value, int):
            return re.sub(r"(::.*=\s*)[\d]+", f"\\g<1>{value}", line)
        elif isinstance(value, float):
            return re.sub(r"(::.*=\s*)[\d.d+\-eE]+", f"\\g<1>{value:.6g}d0", line)
    return line
```

---

## `case.py` — Case Class

```python
# ouxsbli/case.py
import shutil, subprocess, pathlib
from .patcher import patch

class Case:
    def __init__(self, source: str, workdir: str, **params):
        """
        source  : relative path to original case (e.g. "3D_solver/NSTGV")
        workdir : new working directory to create  (e.g. "/tmp/run_001")
        params  : keyword args forwarded to patcher
                  e.g. Re=1600, nx=65, scheme="slau", visc="ns", accuracy=4
        """
        self.source  = pathlib.Path(source).resolve()
        self.workdir = pathlib.Path(workdir).resolve()
        self.params  = params

    def _setup(self):
        if self.workdir.exists():
            shutil.rmtree(self.workdir)
        shutil.copytree(self.source, self.workdir)
        # Apply parameter patches to the copy
        globals_path = self.workdir / "mod_globals.f90"
        original = globals_path.read_text()
        patched  = patch(original, self.params)
        globals_path.write_text(patched)

    def build(self):
        self._setup()
        result = subprocess.run(
            ["make", "clean"],  cwd=self.workdir, check=True, capture_output=True)
        result = subprocess.run(
            ["make"],           cwd=self.workdir, check=True, capture_output=True)

    def run(self, nranks: int = 2):
        data_dir = self.workdir / "data"
        data_dir.mkdir(exist_ok=True)
        subprocess.run(
            ["mpirun", "-n", str(nranks), "./a.out"],
            cwd=self.workdir, check=True)

    def read_vtk(self, step: int):
        from .io.vtk import read_vtr
        return read_vtr(self.workdir / "data" / f"Q{step:05d}.vtr")

    def read_kinetic_energy(self):
        from .io.monitor import read_monitor
        return read_monitor(self.workdir / "data" / "kinetic_energy.d")
```

---

## `io/vtk.py` — VTK Reader

Uses the `vtk` Python library (`pip install vtk`).

```python
# ouxsbli/io/vtk.py
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def read_vtr(path) -> dict:
    """Read a .vtr (RectilinearGrid) file. Returns dict of numpy arrays
    shaped (nz, ny, nx) matching mod_globals dimensions."""
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()

    dims = grid.GetDimensions()   # (nx, ny, nz)
    nx, ny, nz = dims

    point_data = grid.GetPointData()
    arrays = {}
    for i in range(point_data.GetNumberOfArrays()):
        name = point_data.GetArrayName(i)
        arr  = vtk_to_numpy(point_data.GetArray(i))
        # VTK stores in Fortran column-major (x varies fastest)
        arrays[name] = arr.reshape(nz, ny, nx, order="C")
    return arrays
```

---

## `io/monitor.py` — Monitoring File Reader

```python
# ouxsbli/io/monitor.py
import pandas as pd

def read_monitor(path, columns=("t", "value", "normalized")) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None, names=columns)
```

---

## Usage Example

```python
from ouxsbli import Case

# Create a new run from the NSTGV template — original is untouched
case = Case(
    source  = "3D_solver/NSTGV",
    workdir = "/tmp/run_Re1600_slau",
    Re      = 1600,
    M0      = 0.1,
    nx      = 65,  ny = 65,  nz = 65,
    nt      = 500, np = 50,
    scheme  = "slau",
    visc    = "ns",
    accuracy = 2,
)

case.build()          # copies + patches + make
case.run(nranks=2)    # mpirun -n 2 ./a.out

ke     = case.read_kinetic_energy()   # pandas DataFrame: t, value, normalized
fields = case.read_vtk(step=50)       # dict: {"rho": ..., "u": ..., "p": ...}
```

---

## Parameter Sweep (Approach C extension)

```python
from ouxsbli import Case, Sweep

sweep = Sweep(
    base   = {"source": "3D_solver/NSTGV", "nx": 65, "ny": 65, "nz": 65, "nt": 200},
    vary   = {"Re": [400, 800, 1600], "scheme": ["keep", "slau"]},
    outdir = "/tmp/sweep",
)
sweep.build_all(parallel=True)
sweep.run_all(nranks=2)
results = sweep.collect()   # pandas DataFrame with one row per case
```

---

## Critical Files to Create

| File | Purpose |
|------|---------|
| `ouxsbli/__init__.py` | Export `Case`, `Sweep` |
| `ouxsbli/case.py` | `Case` class (setup, build, run, read) |
| `ouxsbli/patcher.py` | `patch()` — text substitution on `mod_globals.f90` copy |
| `ouxsbli/io/vtk.py` | `.vtr` XML reader → numpy |
| `ouxsbli/io/monitor.py` | `kinetic_energy.d` / `entropy.d` → pandas |
| `ouxsbli/tests/test_patcher.py` | Unit tests for patcher (no GPU needed) |
| `ouxsbli/tests/test_io.py` | Unit tests for VTK reader (fixture file) |

**No existing files are modified.**

---

## Verification

1. `test_patcher.py`: take a real `mod_globals.f90`, apply changes, assert new lines match expected patterns
2. `test_io.py`: read an existing `.vtr` file from `2D_solver/BL/data/`, check array shape and value range
3. Full integration: `case.build()` → `a.out` exists in workdir; `case.run()` → VTK files appear; `case.read_vtk()` → shapes are `(nz, ny, nx)`
