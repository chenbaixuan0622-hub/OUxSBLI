"""
Patch config.fypp macro definitions or mod_globals.f90 Fortran parameter declarations.
"""
import re


def patch(src_text: str, changes: dict) -> str:
    """Patch fypp macro definitions (#:set) or Fortran parameter declarations in *src_text*.

    Works transparently on both config.fypp (fypp macros) and mod_globals.f90
    (Fortran parameter declarations).  For params not found in the text, a new
    #:set line is appended **only** when the file is a fypp source (i.e. it
    already contains at least one '#:set' directive).
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

    # Append not-found params only for fypp files (those containing '#:set' directives)
    is_fypp = any("#:set" in l for l in lines)
    for param, value in changes.items():
        found = any(_is_macro_line(l, param) or _is_param_line(l, param) for l in lines)
        if not found and is_fypp:
            lines.append(_rewrite_macro(f"#:set {param} = ", param, value))

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Fortran parameter declaration helpers
# ---------------------------------------------------------------------------

def _is_param_line(line: str, param: str) -> bool:
    """True if *line* is a Fortran parameter declaration for *param*.

    Matches lines of the form::

        [type[(kind)][, attribs]], parameter :: param ...

    Comparison is case-insensitive; only the part before the first ``!``
    (Fortran comment character) is examined.
    """
    code = line.split("!")[0]
    pattern = r"parameter\s*::\s*" + re.escape(param) + r"\b"
    return bool(re.search(pattern, code, re.IGNORECASE))


def _rewrite_param(line: str, param: str, value) -> str:
    """Replace the value of a Fortran parameter declaration.

    Handles scalar types (int, float, bool) and replaces any existing RHS
    expression (including complex ones like ``int(T/dt)``).

    Examples::

        'integer, parameter :: nx = 513'         → 'integer, parameter :: nx = 65'
        'real(8), parameter :: Re = 1600.d0'     → 'real(8), parameter :: Re = 800.0d0'
        'integer, parameter :: nt = int(T/dt)'   → 'integer, parameter :: nt = 500'
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


# ---------------------------------------------------------------------------
# fypp macro helpers  (unchanged from original)
# ---------------------------------------------------------------------------

def _is_macro_line(line: str, param: str) -> bool:
    """True if *line* is the fypp macro definition for *param* (#:set param = ...)."""
    code = line.split("!")[0].split("#")[0] if "!" in line else line
    pattern = r"^\s*#:\s*set\s+" + re.escape(param) + r"\b"
    return bool(re.search(pattern, line, re.IGNORECASE))


def _rewrite_macro(line: str, param: str, value) -> str:
    """Replace the RHS of a fypp macro definition."""
    if isinstance(value, bool):
        new_val = "True" if value else "False"
    elif isinstance(value, int):
        new_val = str(value)
    elif isinstance(value, float):
        new_val = _to_fortran_double(value)
    elif isinstance(value, str):
        val_str = value.strip()
        if (val_str.startswith("'") and val_str.endswith("'")) or (val_str.startswith('"') and val_str.endswith('"')):
            new_val = val_str
        else:
            new_val = f"'{val_str}'"
    else:
        new_val = str(value)

    if "!" in line:
        excl = line.index("!")
        code_part    = line[:excl]
        comment_part = line[excl:]
    else:
        code_part    = line
        comment_part = ""

    new_code = re.sub(
        r"(#:\s*set\s+" + re.escape(param) + r"\s*=\s*)(.+?)\s*$",
        lambda m: m.group(1) + new_val,
        code_part,
        flags=re.IGNORECASE,
    )

    if new_code == code_part and "=" in code_part:
        prefix = code_part.split("=")[0] + "="
        new_code = f"{prefix} {new_val}"

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
