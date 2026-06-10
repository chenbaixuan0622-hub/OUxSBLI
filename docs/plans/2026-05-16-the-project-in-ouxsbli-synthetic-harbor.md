# Plan: Remove IO and KE/Entropy Monitoring Utilities from ouxsbli

## Context
The `ouxsbli/` Python package has grown too large. The IO subpackage (`io/vtk.py`, `io/monitor.py`) and the methods that consume them in `case.py` are to be removed to slim the package down.

## Files to Delete
- `ouxsbli/io/__init__.py`
- `ouxsbli/io/monitor.py`
- `ouxsbli/io/vtk.py`
- `ouxsbli/tests/test_io.py`
- `ouxsbli/io/` directory itself (will be empty after above)

## Files to Modify

### `ouxsbli/__init__.py`
- Remove `from . import io` (line 26)
- Remove `"io"` from `__all__` (line 28)
- Remove `fields = case.read_vtk(step=0)` and `ke = case.read_kinetic_energy()` from the module docstring

### `ouxsbli/case.py`
- Remove module docstring examples referencing `read_vtk` / `read_kinetic_energy` (lines 19–20)
- Remove methods from `Case`:
  - `read_vtk()` (lines 132–149)
  - `read_kinetic_energy()` (lines 151–160)
  - `read_entropy()` (lines 162–171)
  - `vtk_files()` (lines 173–175)
- Remove `Sweep.collect()` method (lines 285–303), which aggregates KE data — it has no purpose without `read_kinetic_energy()`

## Verification
After the changes:
```bash
python -c "from ouxsbli import Case, Sweep, patch; print('OK')"
python -m pytest ouxsbli/tests/test_patcher.py -q
```
Both should pass. The `io` subpackage and its tests will no longer exist.
