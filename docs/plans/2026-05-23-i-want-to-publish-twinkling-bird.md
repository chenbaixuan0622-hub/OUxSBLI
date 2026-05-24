# Plan: Fix ouxinfo conda-forge recipe (PR #33447)

## Context

PR #33447 (`Add ouxinfo`) is open on conda-forge/staged-recipes. Linting passes, but all three platform builds fail with `OverLinkingError` because the package links against OpenMP runtime libraries that are not declared as run dependencies.

## Root Cause

The built extension links against platform-specific OpenMP runtime DSOs. `conda-build` requires these to be explicitly declared as dependencies.

| Platform | Missing DSO         | Package in conda-forge | Note |
|----------|---------------------|------------------------|------|
| Linux    | `libgomp.so.1`      | `_openmp_mutex`        | setup.py uses `-fopenmp` on Linux |
| macOS    | `libomp.dylib`      | `llvm-openmp`          | setup.py uses `-fopenmp` on macOS |
| Windows  | `VCOMP140.DLL`      | `vcomp14`              | setup.py does NOT pass `/openmp` on Windows, but VCOMP140.DLL is still linked (likely via MSVC `/GL` LTCG or boost math headers) |

macOS needs `llvm-openmp` in both **host** (for the linker at build time) and **run** (at runtime).
Linux and Windows only need the runtime package in **run**.

## File to Modify

`recipes/ouxinfo/meta.yaml`

## Changes Required

Add to `requirements.host`:
```yaml
    - llvm-openmp  # [osx]
```

Add to `requirements.run`:
```yaml
    - llvm-openmp  # [osx]
    - _openmp_mutex >=4.5  # [linux]
    - vcomp14  # [win]
```

## Final meta.yaml (requirements section)

```yaml
requirements:
  build:
    - {{ compiler('cxx') }}
    - {{ stdlib('c') }}
  host:
    - python
    - pip
    - setuptools
    - pybind11
    - llvm-openmp  # [osx]
  run:
    - python
    - numpy
    - matplotlib-base
    - scipy
    - tqdm
    - joblib
    - llvm-openmp  # [osx]
    - _openmp_mutex >=4.5  # [linux]
    - vcomp14  # [win]
```

## Verification

After pushing, the PR CI should re-run. Watch for:
1. Linter: must still pass (no new lint issues — the `# [osx]`/`# [linux]`/`# [win]` selector form is correct)
2. Build linux_64: `libgomp.so.1` overlinking error should be gone
3. Build osx_64: `libomp.dylib` missing DSO error should be gone
4. Build win_64: `VCOMP140.DLL` overlinking error should be gone
