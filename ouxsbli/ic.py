"""
Initial condition and boundary condition specifications for OUxSBLI.

These classes translate Python-level IC/BC specs into the ``mod_globals.f90``
parameter values read by a case's :file:`set.f90` at runtime.  Pass an
instance as the ``ic=`` or ``bc=`` argument to :class:`Case`.
"""
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# IC type integer codes (must match select case(ic_type) in set.f90)
# ---------------------------------------------------------------------------
_IC_UNIFORM     = 1
_IC_TAYLOR_GREEN = 2
_IC_RIEMANN     = 3


@dataclass
class UniformIC:
    """Spatially uniform initial condition.

    Parameters
    ----------
    rho:
        Initial density.
    u, v, w:
        Initial velocity components.
    p:
        Initial pressure.  For a quiescent ideal gas at T=T_inf with
        ``R=287.03`` and ``T_inf=300 K``, a convenient default is
        ``rho * R * T_inf ≈ 86100 Pa`` (dimensional) or normalised ``p0 =
        rho/(gamma*M²)`` for a given Mach number.
    """
    rho: float = 1.0
    u:   float = 0.0
    v:   float = 0.0
    w:   float = 0.0
    p:   float = 0.71429  # ≈ 1/(gamma*M²) for M=1/√gamma, gamma=1.4

    def to_params(self) -> dict:
        return {
            "ic_type": _IC_UNIFORM,
            "rho0":    self.rho,
            "u0":      self.u,
            "v0":      self.v,
            "w0":      self.w,
            "p0":      self.p,
        }


@dataclass
class TaylorGreenIC:
    """3D Taylor-Green vortex initial condition.

    Initialises the incompressible-limit Taylor-Green vortex:

    .. math::

        u &=  \\text{Ma}\\sin x \\cos y \\cos z \\\\
        v &= -\\text{Ma}\\cos x \\sin y \\cos z \\\\
        w &= 0 \\\\
        p &= p_0 + \\frac{\\rho_0 \\text{Ma}^2}{16}
               (\\cos 2x + \\cos 2y)(\\cos 2z + 2)

    The domain should be ``[0, 2π]³`` (set ``Lx=Ly=Lz=2*math.pi``).

    Parameters
    ----------
    Ma:
        Peak velocity Mach number (incompressible limit: Ma → 0).
    rho0:
        Reference density.
    p0:
        Reference pressure.
    """
    Ma:   float = 0.1
    rho0: float = 1.0
    p0:   float = 0.71429

    def to_params(self) -> dict:
        return {
            "ic_type": _IC_TAYLOR_GREEN,
            "Ma":      self.Ma,
            "rho0":    self.rho0,
            "p0":      self.p0,
        }


@dataclass
class RiemannIC:
    """1D Riemann (shock-tube) initial condition.

    The domain is divided at ``x = x0 * Lx`` into a high-pressure left state
    and a low-pressure right state.  This is the standard Sod shock-tube
    setup in x with y and z treated as homogeneous.

    Parameters
    ----------
    rho_l, u_l, p_l:
        Left (high-pressure) state.
    rho_r, u_r, p_r:
        Right (low-pressure) state.
    x0:
        Membrane position as a fraction of the domain length Lx (0–1).
    """
    rho_l: float = 1.0
    u_l:   float = 0.0
    p_l:   float = 1.0
    rho_r: float = 0.125
    u_r:   float = 0.0
    p_r:   float = 0.1
    x0:    float = 0.5  # fraction of Lx

    def to_params(self) -> dict:
        return {
            "ic_type":  _IC_RIEMANN,
            "rho_l":    self.rho_l,
            "u_l":      self.u_l,
            "p_l":      self.p_l,
            "rho_r":    self.rho_r,
            "u_r":      self.u_r,
            "p_r":      self.p_r,
            "x0_frac":  self.x0,
        }


# ---------------------------------------------------------------------------
# Boundary condition
# ---------------------------------------------------------------------------

@dataclass
class PeriodicBC:
    """Fully periodic boundary conditions.

    All three directions use periodic (cyclic) ghost-cell exchange and the
    periodic flux stencil (BC_X=BC_Y=BC_Z=False).
    """
    def to_params(self) -> dict:
        return {"BC_X": False, "BC_Y": False, "BC_Z": False}


# ---------------------------------------------------------------------------
# Helper: convert an IC/BC spec (dataclass or dict) to a flat param dict
# ---------------------------------------------------------------------------

def ic_to_params(ic) -> dict:
    """Convert an IC specification to a flat parameter dict.

    Parameters
    ----------
    ic:
        An :class:`UniformIC`, :class:`TaylorGreenIC`, :class:`RiemannIC`
        instance, or a dict with a ``"type"`` key (``"uniform"``,
        ``"tgv"`` / ``"taylor_green"``, ``"riemann"``) plus any
        field overrides.

    Returns
    -------
    dict
        Flat mapping of Fortran parameter names → Python values, suitable
        for passing to :class:`Case`.
    """
    if isinstance(ic, dict):
        ic = _ic_from_dict(ic)
    if hasattr(ic, "to_params"):
        return ic.to_params()
    raise TypeError(f"Unsupported IC specification: {type(ic)!r}")


def bc_to_params(bc) -> dict:
    """Convert a BC specification to a flat parameter dict.

    Parameters
    ----------
    bc:
        A :class:`PeriodicBC` instance, the string ``"periodic"``, or a
        dict with a ``"type"`` key.

    Returns
    -------
    dict
        Flat mapping of fypp macro names → Python values.
    """
    if isinstance(bc, str):
        bc = bc.lower()
        if bc == "periodic":
            return PeriodicBC().to_params()
        raise ValueError(f"Unknown BC type string: {bc!r}.  Supported: 'periodic'.")
    if isinstance(bc, dict):
        bc_type = bc.get("type", "periodic").lower()
        if bc_type == "periodic":
            return PeriodicBC().to_params()
        raise ValueError(f"Unknown BC type: {bc_type!r}.")
    if hasattr(bc, "to_params"):
        return bc.to_params()
    raise TypeError(f"Unsupported BC specification: {type(bc)!r}")


def _ic_from_dict(d: dict):
    """Build an IC dataclass from a dict with a ``"type"`` key."""
    d = dict(d)  # copy
    ic_type = d.pop("type", "uniform").lower()
    if ic_type == "uniform":
        kwargs = {k: v for k, v in d.items() if k in UniformIC.__dataclass_fields__}
        return UniformIC(**kwargs)
    if ic_type in ("tgv", "taylor_green"):
        kwargs = {k: v for k, v in d.items() if k in TaylorGreenIC.__dataclass_fields__}
        return TaylorGreenIC(**kwargs)
    if ic_type == "riemann":
        kwargs = {k: v for k, v in d.items() if k in RiemannIC.__dataclass_fields__}
        return RiemannIC(**kwargs)
    raise ValueError(
        f"Unknown IC type: {ic_type!r}.  Supported: 'uniform', 'tgv', 'riemann'."
    )
