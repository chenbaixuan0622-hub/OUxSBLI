"""
Patch mod_globals.f90 parameter declarations in-place.

OUxSBLI uses Fortran type *kind* (not value) as a compile-time dispatch flag.
For example:
  integer(2) → KEEP scheme     real(2) → SLAU scheme
  integer(4) → Navier-Stokes   integer(8) → LES

KIND_PARAMS lists the parameters that use this kind-dispatch mechanism.
All other parameters are treated as plain scalar replacements.
"""
import re

# Maps each kind-dispatch parameter to a dict of {user-string → Fortran type+kind}.
KIND_MAP: dict[str, dict] = {
    "id_visc": {
        "euler": "integer(2)",
        "ns":    "integer(4)",
        "les":   "integer(8)",
    },
    "id_scheme": {
        "keep":   "integer(2)",
        "slau":   "real(2)",
        "roe":    "real(4)",
        "hybrid": "real(8)",
    },
    "id_accuracy": {
        2: "integer(2)",
        4: "integer(4)",
        6: "integer(8)",
    },
    "id_tvd": {
        "none":    "integer(2)",
        "minmod":  "integer(4)",
        "hybrid":  "integer(8)",
    },
    "id_slau": {
        "slau":    "integer(2)",
        "hr_slau2": "integer(4)",
    },
    "id_rescale": {
        False: "integer(2)",
        True:  "integer(4)",
    },
    "id_recal": {
        "init":    "integer(2)",
        "restart": "integer(4)",
    },
    "id_RungeKutta": {
        "tvd_rk3": "integer(2)",
        "rk4":     "integer(4)",
    },
    "id_gpumpi": {
        False: "integer(2)",
        True:  "integer(4)",
    },
}

KIND_PARAMS = set(KIND_MAP)

# Matches the type+kind prefix in both styles: integer(4)  and  integer(kind=4)
_KIND_PREFIX_RE = re.compile(r'\b(integer|real)\((kind=)?\d+\)')


def patch(src_text: str, changes: dict) -> str:
    """Return src_text with each parameter in *changes* rewritten.

    Scalar parameters (Re, nx, CFL, …) accept int or float values.
    Kind-dispatch parameters (id_visc, id_scheme, …) accept the string/int
    keys defined in KIND_MAP above.

    Example::
        new_src = patch(src, {"Re": 800.0, "nx": 65, "scheme": "slau"})
    """
    lines = src_text.splitlines()
    for i, line in enumerate(lines):
        for param, value in changes.items():
            if not _is_param_line(line, param):
                continue
            if param in KIND_PARAMS:
                lines[i] = _rewrite_kind(line, param, value)
            else:
                lines[i] = _rewrite_scalar(line, param, value)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _is_param_line(line: str, param: str) -> bool:
    """True if *line* is the Fortran parameter declaration for *param*."""
    # Strip inline comment so we don't match inside a comment block
    code = line.split("!")[0]
    if "parameter" not in code.lower():
        return False
    # Look for :: param_name followed by word boundary (not part of a longer name)
    return bool(re.search(r"::\s*" + re.escape(param) + r"\b", code, re.IGNORECASE))


def _rewrite_kind(line: str, param: str, value) -> str:
    """Replace the type+kind prefix to change the dispatch kind."""
    try:
        new_kind = KIND_MAP[param][value]
    except KeyError:
        valid = list(KIND_MAP[param].keys())
        raise ValueError(f"Invalid value {value!r} for {param}. Valid: {valid}") from None
    return _KIND_PREFIX_RE.sub(new_kind, line, count=1)


def _rewrite_scalar(line: str, param: str, value) -> str:
    """Replace the RHS of a scalar parameter declaration."""
    if isinstance(value, bool):
        new_val = ".true." if value else ".false."
    elif isinstance(value, int):
        new_val = str(value)
    elif isinstance(value, float):
        # Fortran double-precision literal: 1600.0 → 1600.0d0
        new_val = _to_fortran_double(value)
    else:
        new_val = str(value)

    # Separate code from trailing inline comment
    if "!" in line:
        excl = line.index("!")
        code_part    = line[:excl]
        comment_part = line[excl:]
    else:
        code_part    = line
        comment_part = ""

    # Replace everything after the = sign in the code part
    new_code = re.sub(
        r"(::\s*" + re.escape(param) + r"\s*=\s*)(.+?)\s*$",
        lambda m: m.group(1) + new_val,
        code_part,
        flags=re.IGNORECASE,
    )
    return new_code + comment_part


def _to_fortran_double(v: float) -> str:
    """Format a Python float as a Fortran double-precision literal (e.g. 1600.0d0)."""
    s = f"{v:.10g}"
    if "e" in s or "E" in s:
        s = s.replace("e", "d").replace("E", "d")
    else:
        if "." not in s:
            s += ".0"
        s += "d0"
    return s
