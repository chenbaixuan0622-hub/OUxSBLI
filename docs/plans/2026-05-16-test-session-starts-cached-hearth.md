# Plan: Fit test files to updated vtk_reader.py

## Context

The old `vtk_reader.py` used a hand-rolled binary VTK parser (`read_vtr` / `_read`). It had a bug: when the z-coordinate block size was not a multiple of 4 bytes, `np.frombuffer(..., dtype="<f4")` raised `ValueError: buffer size must be a multiple of element size`. The user rewrote `vtk_reader.py` to use the proper VTK Python library instead.

The new API in `ouxsbli/tests/utils/vtk_reader.py`:
- `extract_number(filename)` — extracts the trailing number from a Q*.vtr filename
- `getGrid(file_path)` → `(nx, ny, nz, x, y, z)` — reads coordinate arrays
- `getVector(file_path, Nx, Ny, Nz, name)` → `(u, v, w)` shape `(nz,ny,nx)`
- `getScalar(file_path, Nx, Ny, Nz, name)` → array shape `(nz,ny,nx)`
- `getQ(file_path, Nx, Ny, Nz, reader=None)` → `(rho, u, v, w, p)` each shape `(nz,ny,nx)`

The four test files that use VTK data still import the old API (`latest_vtr`, `read_vtr`) and access a dict `d` with keys `"rho"`, `"p"`, `"x"`, `"y"`, `"ni"`, `"nj"`, `"nk"`. These need to be updated.

## Changes

### 1. `ouxsbli/tests/utils/vtk_reader.py` — add `latest_vtr`

`latest_vtr` is a utility that logically belongs alongside `extract_number` (which it uses). Add it at the end of the file:

```python
def latest_vtr(data_dir):
    import glob
    files = glob.glob(os.path.join(str(data_dir), "Q*.vtr"))
    if not files:
        raise FileNotFoundError(f"No Q*.vtr files found in {data_dir}")
    return max(files, key=lambda f: extract_number(os.path.basename(f)))
```

### 2. `ouxsbli/tests/test_os.py`

Replace import + usage of `read_vtr` with `getGrid` / `getQ`:

```python
from .utils.vtk_reader import latest_vtr, getGrid, getQ
...
vtr_path = latest_vtr(data_dir)
nx, ny, nz, x, y, z = getGrid(vtr_path)
rho, u, v, w, p = getQ(vtr_path, nx, ny, nz)
rho_2d = rho[0].astype(float)   # (ny, nx) — nz=1 for 2-D
p_2d   = p[0].astype(float)
x2d, y2d = np.meshgrid(x, y)
```

### 3. `ouxsbli/tests/test_ivst.py`

```python
from .utils.vtk_reader import latest_vtr, getGrid, getQ
...
vtr_path = latest_vtr(data_dir)
ni, nj, nk, x, y, z = getGrid(vtr_path)   # ni=nx, nj=ny, nk=nz
rho, u, v, w, p = getQ(vtr_path, ni, nj, nk)
# rho/p shape: (nk, nj, ni)
jmid = nj // 2
kmid = nk // 2
rho_num = rho[kmid, jmid, :].astype(float)
p_num   = p[kmid, jmid, :].astype(float)
x_arr   = x.astype(float)
```

### 4. `ouxsbli/tests/test_evc.py`

Inside `_run_evc`:

```python
from .utils.vtk_reader import latest_vtr, getGrid, getQ
...
vtr_path = latest_vtr(data_dir)
nx, ny, nz, x, y, z = getGrid(vtr_path)
rho, u, v, w, p = getQ(vtr_path, nx, ny, nz)
rho_num = rho[0].astype(float)   # (ny, nx)
x1d = x.astype(float)
y1d = y.astype(float)
```

### 5. `ouxsbli/tests/test_corn.py`

```python
from .utils.vtk_reader import latest_vtr, getGrid, getQ
...
vtr_path = latest_vtr(data_dir)
ni, nj, nk, x, y, z = getGrid(vtr_path)
rho, u, v, w, p = getQ(vtr_path, ni, nj, nk)
# rest of slicing uses ni, nj, nk directly (was d["ni"], d["nj"], d["nk"])
```

## Key dimension mapping

| old dict key | new variable | getGrid return position |
|---|---|---|
| `d["ni"]` | `nx` (or `ni`) | position 0 |
| `d["nj"]` | `ny` (or `nj`) | position 1 |
| `d["nk"]` | `nz` (or `nk`) | position 2 |
| `d["x"]`  | `x` | position 3 |
| `d["y"]`  | `y` | position 4 |
| `d["rho"]` | `rho` — shape `(nz,ny,nx)` | from `getQ` |
| `d["p"]`   | `p`   — shape `(nz,ny,nx)` | from `getQ` |

## Files to modify

- `ouxsbli/tests/utils/vtk_reader.py` — add `latest_vtr`
- `ouxsbli/tests/test_os.py`
- `ouxsbli/tests/test_ivst.py`
- `ouxsbli/tests/test_evc.py`
- `ouxsbli/tests/test_corn.py`

`ouxsbli/tests/test_etgv.py` — no VTK reads, no changes needed.

## Verification

After editing, confirm there are no import errors:
```bash
python -m pytest ouxsbli/tests/test_os.py --collect-only
python -m pytest ouxsbli/tests/test_ivst.py --collect-only
python -m pytest ouxsbli/tests/test_evc.py --collect-only
python -m pytest ouxsbli/tests/test_corn.py --collect-only
```
