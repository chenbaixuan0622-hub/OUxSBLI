# Plan: Fix `case.build()` failing when NVIDIA HPC SDK is not in PATH

## Context

`Case.build()` in `ouxsbli/case.py` calls `subprocess.run(["make"], ...)` without passing a custom environment. When invoked from a Jupyter notebook, the kernel's PATH resolves `mpif90` to the system gfortran (`/usr/bin/mpif90`) instead of the NVIDIA HPC SDK wrapper. The Makefiles use NVIDIA-only flags (`-cuda`, `-acc`, `-Mpreprocess`, etc.), so gfortran immediately errors out.

The NVIDIA HPC SDK 24.7 is installed and `nvfortran` exists at:
`/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/bin/nvfortran`

The correct NVIDIA `mpif90` is at:
`/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/hpcx/bin/mpif90`

## Fix

Modify `_run_cmd` in `ouxsbli/case.py` to pass a PATH-augmented environment to every `subprocess.run` call, prepending the NVIDIA HPC SDK compiler and MPI directories when they are found.

### File to change

**`ouxsbli/case.py`**

1. Add `import os` at the top (alongside existing stdlib imports).

2. Add a `_build_env()` helper method to `Case`:

```python
def _build_env(self) -> dict:
    """Return os.environ copy with NVIDIA HPC SDK paths prepended if found."""
    env = os.environ.copy()
    nvfortran = shutil.which("nvfortran")
    if nvfortran is None:
        hpcsdk_root = pathlib.Path("/opt/nvidia/hpc_sdk")
        candidates = sorted(
            hpcsdk_root.glob("Linux_x86_64/*/compilers/bin/nvfortran"),
            reverse=True,
        )
        if candidates:
            nvfortran = str(candidates[0])
    if nvfortran is None:
        raise RuntimeError(
            "NVIDIA HPC SDK not found. Install it from https://developer.nvidia.com/hpc-sdk "
            "or add nvfortran to PATH."
        )
    compiler_bin = pathlib.Path(nvfortran).parent
    sdk_root = compiler_bin.parent.parent          # .../Linux_x86_64/<ver>
    extra = [str(compiler_bin)]
    for mpi_subpath in ("comm_libs/hpcx/bin", "comm_libs/openmpi4/bin"):
        mpi_bin = sdk_root / mpi_subpath
        if mpi_bin.exists():
            extra.append(str(mpi_bin))
            break
    env["PATH"] = ":".join(extra) + ":" + env.get("PATH", "")
    return env
```

3. Change `_run_cmd` to pass `env=self._build_env()`:

```python
def _run_cmd(self, cmd: list[str]) -> None:
    result = subprocess.run(
        cmd,
        cwd=self.workdir,
        capture_output=True,
        text=True,
        env=self._build_env(),
    )
    if result.returncode != 0:
        raise RuntimeError(...)
```

## How it works

- If NVIDIA's `nvfortran` is already on PATH (e.g. in a properly configured shell), `shutil.which` finds it and the correct paths are prepended.
- If not (Jupyter kernel, CI), it falls back to globbing the well-known `/opt/nvidia/hpc_sdk` install prefix, picking the lexicographically latest version.
- If NVIDIA HPC SDK is absent entirely, `_build_env` raises a clear `RuntimeError` before any subprocess is launched:
  ```
  RuntimeError: NVIDIA HPC SDK not found. Install it from https://developer.nvidia.com/hpc-sdk
  or set nvfortran on PATH.
  ```

## Verification

After the fix, re-run cell 4 of `tutorials/02_run_simulation.ipynb`:
```python
case.build()
```
Expected: no `RuntimeError`; `a.out` is created in the workdir.
