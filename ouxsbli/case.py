"""
Case — lifecycle manager for a single OUxSBLI simulation run.

Usage::

    from ouxsbli import Case

    case = Case(
        source  = "3D_solver/NSTGV",
        workdir = "/tmp/run_01",
        Re      = 800.0,
        nx      = 65,
        scheme  = "slau",
        visc    = "ns",
        accuracy = 2,
    )
    case.build()
    case.run(nranks=2)
"""
import os
import shutil
import subprocess
import pathlib
from typing import Any
from .patcher import patch


# Mapping from friendly parameter alias → Fortran parameter name used in patcher
_ALIAS = {
    "scheme":   "id_scheme",
    "visc":     "id_visc",
    "accuracy": "id_accuracy",
    "tvd":      "id_tvd",
    "slau":     "id_slau",
    "rescale":  "id_rescale",
    "recal":    "id_recal",
    "rk":       "id_RungeKutta",
    "gpumpi":   "id_gpumpi",
}


class Case:
    """Manage one simulation run (copy source → patch → build → run → read)."""

    def __init__(self, source: str, workdir: str, **params: Any) -> None:
        """
        Parameters
        ----------
        source:
            Path to the original case directory (e.g. ``"3D_solver/NSTGV"``).
            This directory is never modified.
        workdir:
            Path for the new working directory that will be created.
        **params:
            Parameter overrides forwarded to :func:`~ouxsbli.patcher.patch`.
            Scalar parameters (``Re``, ``nx``, ``CFL``, …) accept int/float.
            Kind-dispatch parameters accept the friendly strings defined in
            ``patcher.KIND_MAP`` (e.g. ``scheme="slau"``, ``visc="ns"``).
        """
        repo_root = pathlib.Path(__file__).resolve().parent.parent
        self.source  = (repo_root / source).resolve()
        self.workdir = pathlib.Path(workdir).resolve()
        # Expand friendly aliases to the exact Fortran parameter names
        self.params: dict[str, Any] = {
            _ALIAS.get(k, k): v for k, v in params.items()
        }
        self._built = False

    # ------------------------------------------------------------------
    # Setup (copy + patch)
    # ------------------------------------------------------------------

    def _setup(self) -> None:
        """Copy source directory to workdir and patch mod_globals.f90."""
        if self.workdir.exists():
            shutil.rmtree(self.workdir)
        shutil.copytree(self.source, self.workdir)

        globals_path = self.workdir / "mod_globals.f90"
        if not globals_path.exists():
            raise FileNotFoundError(f"mod_globals.f90 not found in {self.source}")

        original = globals_path.read_text()
        patched  = patch(original, self.params)
        globals_path.write_text(patched)

        # The Makefile's vpath uses paths relative to the original source dir.
        # make runs from workdir, so rewrite them as absolute paths.
        makefile_path = self.workdir / "Makefile"
        if makefile_path.exists():
            lines = makefile_path.read_text().splitlines(keepends=True)
            for i, line in enumerate(lines):
                if line.startswith("vpath %f90 "):
                    rel_paths = line[len("vpath %f90 "):].strip().split(":")
                    abs_paths = [str((self.source / p).resolve()) for p in rel_paths]
                    lines[i] = "vpath %f90 " + ":".join(abs_paths) + "\n"
            makefile_path.write_text("".join(lines))

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def build(self) -> None:
        """Copy source, apply patches, and compile with make."""
        self._setup()
        self._run_cmd(["make", "clean"])
        self._run_cmd(["make"])
        self._built = True

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, nranks: int = 2) -> None:
        """Launch the simulation with mpirun.

        Parameters
        ----------
        nranks:
            Number of MPI ranks (default 2).
        """
        (self.workdir / "data").mkdir(exist_ok=True)
        self._run_cmd(["mpirun", "-n", str(nranks), "./a.out"])

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _build_env(self) -> dict:
        """Return os.environ copy with NVIDIA HPC SDK paths prepended."""
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
                "NVIDIA HPC SDK not found. Install it from "
                "https://developer.nvidia.com/hpc-sdk or add nvfortran to PATH."
            )
        compiler_bin = pathlib.Path(nvfortran).parent
        sdk_root = compiler_bin.parent.parent
        extra = [str(compiler_bin)]
        for mpi_subpath in ("comm_libs/hpcx/bin", "comm_libs/openmpi4/bin"):
            mpi_bin = sdk_root / mpi_subpath
            if mpi_bin.exists():
                extra.append(str(mpi_bin))
                break
        env["PATH"] = ":".join(extra) + ":" + env.get("PATH", "")
        return env

    def _run_cmd(self, cmd: list[str]) -> None:
        result = subprocess.run(
            cmd,
            cwd=self.workdir,
            capture_output=True,
            text=True,
            env=self._build_env(),
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Command {cmd} failed (exit {result.returncode}):\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )

    def __repr__(self) -> str:
        return (
            f"Case(source={self.source.name!r}, "
            f"workdir={str(self.workdir)!r}, "
            f"params={self.params})"
        )

