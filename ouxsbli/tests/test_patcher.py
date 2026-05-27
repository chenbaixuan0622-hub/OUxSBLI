"""
Unit tests for ouxsbli.patcher.

These tests use the real mod_globals.f90 from the repository as input, so no
GPU or build environment is needed.
"""

import pathlib
import pytest

from ouxsbli.patcher import patch, _is_param_line, _to_fortran_double

REPO_ROOT  = pathlib.Path(__file__).resolve().parents[2]
NSTGV_GLOB = REPO_ROOT / "3D_solver" / "NSTGV" / "mod_globals.f90"
BL_GLOB    = REPO_ROOT / "2D_solver" / "BL"   / "mod_globals.f90"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load(path: pathlib.Path) -> str:
    return path.read_text()


def _find_line(text: str, param: str) -> str:
    """Return the declaration line for *param* (raises if not found)."""
    for line in text.splitlines():
        if _is_param_line(line, param):
            return line
    raise AssertionError(f"Parameter '{param}' not found in patched text")


# ---------------------------------------------------------------------------
# _to_fortran_double
# ---------------------------------------------------------------------------

class TestToFortranDouble:
    def test_integer_valued(self):
        assert _to_fortran_double(1600.0) == "1600.0d0"

    def test_fractional(self):
        result = _to_fortran_double(0.03)
        assert result.endswith("d0")
        assert float(result.replace("d0", "e0")) == pytest.approx(0.03, rel=1e-6)

    def test_large(self):
        result = _to_fortran_double(1.25)
        assert "1.25" in result and result.endswith("d0")


# ---------------------------------------------------------------------------
# Scalar parameters (3D NSTGV)
# ---------------------------------------------------------------------------

class TestScalarPatch3D:
    @pytest.fixture
    def src(self):
        return _load(NSTGV_GLOB)

    def test_Re(self, src):
        out  = patch(src, {"Re": 800.0})
        line = _find_line(out, "Re")
        assert "800" in line

    def test_nx(self, src):
        out  = patch(src, {"nx": 65})
        line = _find_line(out, "nx")
        assert "65" in line
        assert "513" not in line

    def test_np(self, src):
        out  = patch(src, {"np": 50})
        line = _find_line(out, "np")
        assert "50" in line

    def test_nt_replaces_expression(self, src):
        # nt is computed from a formula in NSTGV; after patching it should be a literal
        out  = patch(src, {"nt": 500})
        line = _find_line(out, "nt")
        assert "500" in line
        assert "int(" not in line

    def test_CFL(self, src):
        out  = patch(src, {"CFL": 0.1})
        line = _find_line(out, "CFL")
        assert "0.1" in line

    def test_multiple_scalars(self, src):
        out = patch(src, {"Re": 400.0, "nx": 33, "ny": 33, "nz": 33})
        assert "400" in _find_line(out, "Re")
        assert "33"  in _find_line(out, "nx")
        assert "33"  in _find_line(out, "ny")
        assert "33"  in _find_line(out, "nz")

    def test_original_unchanged(self, src):
        out = patch(src, {"Re": 800.0})
        # patch() must not mutate the original string
        assert "1600" in src


# ---------------------------------------------------------------------------
# 2D BL case
# ---------------------------------------------------------------------------

class TestPatch2D:
    @pytest.fixture
    def src(self):
        return _load(BL_GLOB)

    def test_nx(self, src):
        out  = patch(src, {"nx": 129})
        line = _find_line(out, "nx")
        assert "129" in line

    def test_M0(self, src):
        out  = patch(src, {"M0": 1.5})
        line = _find_line(out, "M0")
        assert "1.5" in line


# ---------------------------------------------------------------------------
# Alias resolution and value normalisation (smoke tests, no build)
# ---------------------------------------------------------------------------

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
        """'euler' is normalised to the fypp-correct 'Euler' (case-sensitive)."""
        from ouxsbli.case import _VALUE_NORMALIZE
        assert _VALUE_NORMALIZE["VISC"]["euler"] == "Euler"
        assert _VALUE_NORMALIZE["VISC"]["ns"]    == "NS"

    def test_value_normalize_tvd(self):
        from ouxsbli.case import _VALUE_NORMALIZE
        assert _VALUE_NORMALIZE["TVD"]["none"]   == "none"
        assert _VALUE_NORMALIZE["TVD"]["hybrid"] == "muscl4"   # safe alias
