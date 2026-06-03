# Plan: Fix patcher.py and case.py for correct param routing

## Context

`pytest` fails for two reasons:

1. **ImportError**: `test_patcher.py` imports `_is_param_line` from `ouxsbli.patcher`, but that function does not exist, breaking every test collection.
2. **Wrong routing**: `case.py` routes *all* `Case()` kwargs to `config.fypp` via `patch()`. Variables like `nx`, `ny`, `Re`, `dt`, `nt`, `np`, `Lx`, `Ly`, `gamma`, etc. are **not** in `config.fypp` — they live in `mod_globals.f90` — so they are silently appended as useless `#:set NX = 65` macros instead of patching the real Fortran parameters.

Additionally `test_patcher.py` has two stale tests (`TestKindPatch3D`, `TestAlias`) that reference `id_*` variables and aliases which are being removed from the codebase. These tests must be updated to reflect the current design.

**Constraints:** No changes to any Fortran source files (`.f90` / `.fypp`).

---

## Part 1 — `ouxsbli/patcher.py`

### Add `_is_param_line(line, param)`

Detects a Fortran `parameter` declaration for the named variable:

```python
def _is_param_line(line: str, param: str) -> bool:
    """True if *line* is a Fortran parameter declaration for *param*.

    Matches:  [type[(kind)][, attribs]], parameter :: param  ...
    Stops at the '!' comment character; case-insensitive.
    """
    code = line.split("!")[0]
    pattern = r"parameter\s*::\s*" + re.escape(param) + r"\b"
    return bool(re.search(pattern, code, re.IGNORECASE))
```

### Add `_rewrite_param(line, param, value)`

Replaces the RHS of a Fortran parameter declaration with a new literal value.

```python
def _rewrite_param(line: str, param: str, value) -> str:
    """Replace the value of a Fortran parameter declaration.

    e.g.  'integer, parameter :: nx = 513'  →  'integer, parameter :: nx = 65'
          'real(8), parameter :: Re = 1600.d0' → 'real(8), parameter :: Re = 800.0d0'
          'integer, parameter :: nt = int(T/dt)' → 'integer, parameter :: nt = 500'
    """
    if isinstance(value, bool):
        new_val = ".true." if value else ".false."
    elif isinstance(value, int):
        new_val = str(value)
    elif isinstance(value, float):
        new_val = _to_fortran_double(value)
    else:
        new_val = str(value)

    # Split off trailing Fortran comment
    if "!" in line:
        excl = line.index("!")
        code_part, comment_part = line[:excl], line[excl:]
    else:
        code_part, comment_part = line, ""

    # Replace everything after  'parameter :: param ='  (handles expressions too)
    new_code = re.sub(
        r"(parameter\s*::\s*" + re.escape(param) + r"\s*=\s*).*$",
        lambda m: m.group(1) + new_val,
        code_part,
        flags=re.IGNORECASE,
    )
    return new_code + comment_part
```

### Modify `patch()`

Try fypp-macro match first, then Fortran-parameter match. Append new lines **only** for fypp files (those containing `#:set` directives).

```python
def patch(src_text: str, changes: dict) -> str:
    """Patch fypp macro definitions (#:set) or Fortran parameter declarations.

    Handles both config.fypp and mod_globals.f90 transparently.
    """
    lines = src_text.splitlines()
    for i, line in enumerate(lines):
        for param, value in changes.items():
            if _is_macro_line(line, param):
                lines[i] = _rewrite_macro(line, param, value)
                break
            elif _is_param_line(line, param):
                lines[i] = _rewrite_param(line, param, value)
                break

    # Append not-found params only for fypp files
    is_fypp = any("#:set" in l for l in lines)
    for param, value in changes.items():
        found = any(_is_macro_line(l, param) or _is_param_line(l, param) for l in lines)
        if not found and is_fypp:
            lines.append(_rewrite_macro(f"#:set {param} = ", param, value))

    return "\n".join(lines) + "\n"
```

---

## Part 2 — `ouxsbli/case.py`

### New module-level helpers

```python
# String → integer normalisation for the RK stage count expected by config.fypp
_RK_NORMALIZE: dict[str, int] = {
    "tvd_rk3": 3, "rk3": 3,
    "rk4": 4, "classical_rk4": 4,
}

# Correct config.fypp casing for string-valued parameters.
# fypp comparisons are case-sensitive; 'EULER' ≠ 'Euler'.
_VALUE_NORMALIZE: dict[str, dict[str, str]] = {
    "SCHEME":       {"keep": "KEEP", "slau": "SLAU", "roe": "Roe", "hybrid": "Hybrid"},
    "VISC":         {"euler": "Euler", "ns": "NS", "les": "LES"},
    "TVD":          {"none": "none", "minmod": "minmod",
                     "muscl4": "muscl4", "hybrid": "muscl4"},  # hybrid → muscl4
    "SLAU_VARIANT": {"slau": "SLAU", "hrslau2": "HRSLAU2"},
}
```

### Updated `_ALIAS`

Aliases stay pointed at config.fypp macro names (no `id_*`):

```python
_ALIAS = {
    "scheme":     "SCHEME",
    "visc":       "VISC",
    "accuracy":   "ORDER",
    "visc_order": "VISC_ORDER",
    "tvd":        "TVD",
    "slau":       "SLAU_VARIANT",
    "rescale":    "RESCALE",
    "recal":      "RESTART",
    "rk":         "RK",
    "gpumpi":     "GPUMPI",
    "bc_x":       "BC_X",
    "bc_y":       "BC_Y",
    "commz":      "COMMZ",
}
```

### Updated `__init__`

Remove the blanket `v.upper()` that was corrupting fypp string values like `'Euler'` and `'none'`:

```python
def __init__(self, source: str, workdir: str, **params: Any) -> None:
    ...
    self.params: dict[str, Any] = {
        _ALIAS.get(k.lower(), k): v   # alias lookup only; no blanket uppercase
        for k, v in params.items()
    }
    self._built = False
```

### Rewrite `_setup()`

```python
def _setup(self) -> None:
    if self.workdir.exists():
        shutil.rmtree(self.workdir)
    shutil.copytree(
        self.source, self.workdir,
        ignore=shutil.ignore_patterns("build", "CMakeCache.txt", "CMakeFiles", "*.cmake"),
    )

    config_path  = self.workdir / "config.fypp"
    globals_path = self.workdir / "mod_globals.f90"
    if not config_path.exists():
        raise FileNotFoundError(f"config.fypp not found in {self.source}")
    if not globals_path.exists():
        raise FileNotFoundError(f"mod_globals.f90 not found in {self.source}")

    config_text  = config_path.read_text()
    globals_text = globals_path.read_text()
    config_lines = config_text.splitlines()

    config_params:  dict[str, Any] = {}
    globals_params: dict[str, Any] = {}

    for k, v in self.params.items():
        # Normalise RK string → int
        if k == "RK" and isinstance(v, str):
            v = _RK_NORMALIZE.get(v.lower(), v)

        if any(_is_macro_line(line, k) for line in config_lines):
            # Param is a fypp macro in config.fypp → apply casing normalisation
            if isinstance(v, str) and k in _VALUE_NORMALIZE:
                v = _VALUE_NORMALIZE[k].get(v.lower(), v)
            config_params[k] = v
        else:
            # Param is not in config.fypp → route to mod_globals.f90
            globals_params[k] = v

    config_path.write_text(patch(config_text, config_params))
    globals_path.write_text(patch(globals_text, globals_params))
```

**Also import** `_is_macro_line` at the top of `case.py`:
```python
from .patcher import patch, _is_macro_line
```

---

## Part 3 — `ouxsbli/tests/test_patcher.py`

Remove the stale `id_*`-based tests (those parameters are being removed from the codebase) and fix `TestAlias` to reflect the current aliases.

### Remove `TestKindPatch3D` entirely

This class tested patching of `id_visc`, `id_scheme`, etc. in `mod_globals.f90`. Since those parameters do not exist there and are being removed, delete the whole class.

### Update `TestPatch2D`

Remove only the `test_scheme_keep_to_slau` method (references `id_scheme`). Keep `test_nx` and `test_M0`.

### Replace `TestAlias`

```python
class TestAlias:
    def test_alias_scheme(self):
        """Case maps 'scheme' → 'SCHEME' (config.fypp macro)."""
        from ouxsbli.case import _ALIAS
        assert _ALIAS["scheme"] == "SCHEME"
        assert _ALIAS["visc"]   == "VISC"

    def test_alias_accuracy(self):
        from ouxsbli.case import _ALIAS
        assert _ALIAS["accuracy"] == "ORDER"

    def test_value_normalize_visc(self):
        """'euler' string is normalised to the fypp-correct 'Euler'."""
        from ouxsbli.case import _VALUE_NORMALIZE
        assert _VALUE_NORMALIZE["VISC"]["euler"]  == "Euler"
        assert _VALUE_NORMALIZE["VISC"]["ns"]     == "NS"

    def test_value_normalize_tvd(self):
        from ouxsbli.case import _VALUE_NORMALIZE
        assert _VALUE_NORMALIZE["TVD"]["none"]    == "none"
        assert _VALUE_NORMALIZE["TVD"]["hybrid"]  == "muscl4"  # safe alias
```

---

## Routing summary

| `Case()` kwarg | Alias | Goes to | Value applied |
|---|---|---|---|
| `scheme="slau"` | `SCHEME` | config.fypp | `'SLAU'` |
| `visc="euler"` | `VISC` | config.fypp | `'Euler'` |
| `tvd="none"` | `TVD` | config.fypp | `'none'` |
| `tvd="hybrid"` | `TVD` | config.fypp | `'muscl4'` (via `_VALUE_NORMALIZE`) |
| `accuracy=4` | `ORDER` | config.fypp | `4` |
| `rk="tvd_rk3"` | `RK` → normalised | config.fypp | `3` |
| `nx=65` | none | mod_globals.f90 | `65` |
| `Re=800.0` | none | mod_globals.f90 | `800.0d0` |
| `nt=int(endT/dt)` | none | mod_globals.f90 | literal integer |
| `p_tot=100e3` | none | mod_globals.f90 | `100000.0d0` |

---

## Critical Files

- [ouxsbli/patcher.py](ouxsbli/patcher.py) — add `_is_param_line`, `_rewrite_param`; modify `patch()`
- [ouxsbli/case.py](ouxsbli/case.py) — add `_RK_NORMALIZE`, `_VALUE_NORMALIZE`; update `_ALIAS`; update `__init__`; rewrite `_setup()`
- [ouxsbli/tests/test_patcher.py](ouxsbli/tests/test_patcher.py) — remove `TestKindPatch3D`; update `TestPatch2D`, `TestAlias`

No Fortran source files are modified.

---

## Verification

```bash
cd /home/jhatayama/ouxsbli/fypp/OUxSBLI

# 1. Import smoke test
python -c "from ouxsbli.patcher import patch, _is_param_line, _to_fortran_double; print('import OK')"

# 2. Unit tests
pytest ouxsbli/tests/test_patcher.py -v
# Expected: ALL tests pass (TestToFortranDouble, TestScalarPatch3D, TestPatch2D, TestAlias)

# 3. Routing smoke test (no GPU needed)
python - <<'EOF'
from ouxsbli.patcher import patch, _is_param_line

# config.fypp patching
config = "#:set VISC = 'NS'\n#:set SCHEME = 'KEEP'\n#:set ORDER = 6\n#:set TVD = 'none'\n"
out = patch(config, {"VISC": "Euler", "ORDER": 4, "TVD": "none"})
assert "VISC = 'Euler'" in out, out
assert "ORDER = 4" in out, out
print("config.fypp patch: OK")

# mod_globals.f90 patching
mg = "  integer, parameter :: nx = 513\n  real(8), parameter :: Re = 1600.d0\n"
out2 = patch(mg, {"nx": 65, "Re": 800.0})
assert "nx = 65" in out2, out2
assert "800" in out2, out2
print("mod_globals.f90 patch: OK")
EOF
```
