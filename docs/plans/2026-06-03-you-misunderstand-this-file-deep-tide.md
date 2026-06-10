# Plan: Rename TVD option strings `minmod`→`tvd` and `muscl4`→`hybrid`

## Context

The kind-dispatch table for `id_tvd` uses the wrong canonical names. The correct mapping is:

| kind | canonical name | meaning |
|------|---------------|---------|
| `integer(2)` | `none` | no TVD (correct, unchanged) |
| `integer(4)` | `tvd` | TVD reconstruction (was `minmod`) |
| `integer(8)` | `hybrid` | hybrid high-order reconstruction (was `muscl4`) |

This affects the fypp template, all case `config.fypp` files, the Python normalization table, its tests, and documentation.

---

## Changes

### 1. `src/mod_constant.f90.fypp` — fypp condition strings
```fortran
#:elif TVD == 'tvd'      ← was 'minmod'
#:elif TVD == 'hybrid'   ← was 'muscl4'
```

### 2. `config.fypp` files — canonical value in each case
- `3D_solver/STZ/config.fypp` line 12: `TVD = 'minmod'` → `TVD = 'tvd'`
- `3D_solver/ETGV/config.fypp` line 12: `TVD = 'muscl4'` → `TVD = 'hybrid'`
- `3D_solver/DHIT/config.fypp` line 12: `TVD = 'muscl4'` → `TVD = 'hybrid'`
- `3D_solver/IVST/config.fypp` line 12: `TVD = 'muscl4'` → `TVD = 'hybrid'`
- `3D_solver/SBLI/config.fypp` line 12: `TVD = 'muscl4'` → `TVD = 'hybrid'`
- `3D_solver/NSTGV/config.fypp` line 12: `TVD = 'muscl4'` → `TVD = 'hybrid'`
- `ouxsbli/tests/tmp/etgv/config.fypp` line 12: `TVD = 'muscl4'` → `TVD = 'hybrid'`

### 3. `ouxsbli/case.py` — `_VALUE_NORMALIZE["TVD"]`
```python
"TVD": {
    "none":   "none",
    "tvd":    "tvd",
    "hybrid": "hybrid",
    "minmod": "tvd",    # backward-compat alias
    "muscl4": "hybrid", # backward-compat alias
},
```
Remove the old comment `# "hybrid" → muscl4 limiter`.

### 4. `ouxsbli/tests/test_patcher.py` — update TVD assertion
```python
assert _VALUE_NORMALIZE["TVD"]["hybrid"] == "hybrid"   # canonical name
```
Also add/update assertions for `"tvd"` → `"tvd"` and backward-compat aliases if tested.

### 5. Documentation
- `CLAUDE.md` lines 60, 213, 238: replace `'minmod'`/`'muscl4'` with `'tvd'`/`'hybrid'`; update kind-dispatch table row for `id_tvd`
- `docs/configuration.md` lines 31, 57, 117: same substitutions
- `docs/installation.md` line 179: update `id_tvd` table row
- `docs/index.md` line 55: update feature description

---

## Verification

```bash
# Confirm no stale 'minmod'/'muscl4' remain as TVD= values
grep -r "TVD\s*=\s*'minmod'\|TVD\s*=\s*'muscl4'" .

# Confirm fypp template uses new names
grep -n "minmod\|muscl4" src/mod_constant.f90.fypp

# Run Python tests
pytest ouxsbli/tests/test_patcher.py -v

# Rebuild one case to verify fypp expansion is clean
cd 3D_solver/NSTGV && cmake -B build && cmake --build build -j
```
