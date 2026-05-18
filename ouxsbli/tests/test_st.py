"""
Inviscid shock tube (Sod problem): compare the 2D_solver/ST output
with the exact Sod solution.

The L1 error of density is required to be < 1 %.
"""
import pathlib
import pytest
import numpy as np
from ouxsbli import Case
from .utils.vtk_reader import latest_vtr, getGrid, getQ
from .utils.sod_exact import solve as sod_solve


L1_TOL = 0.01   # 1 % relative L1 error


@pytest.mark.integration
def test_ivst_vs_sod_exact(tmp_path):
    workdir = "./tmp/st"
    case = Case(
        source   = "2D_solver/ST",
        workdir  = workdir,
        visc     = "euler",
        tvd      = "hybrid",
        scheme   = "slau",
        accuracy = 6,
        rk       = "tvd_rk3",
        nx = 513,
        ny = 7,
        np = 1,
    )
    case.build()
    case.run(nranks=2)
    data_dir = pathlib.Path(workdir) / "data"
    vtr_path = latest_vtr(data_dir)
    ni, nj, nk, x, y, z = getGrid(vtr_path)
    rho, u, v, w, p = getQ(vtr_path, ni, nj, nk)
    # ST is essentially 1-D in x; take a middle slice in j and k
    # rho/p shape: (nj, ni)
    rho_num = rho[0, -1, :].astype(float)  # (ni,)
    p_num   = p[0, -1, :].astype(float)
    x = x.astype(float)  # rank-0 x-coords, shape (ni,)
    # IVST parameters
    Rgas  = 287.03e0
    T     = 300.e0
    C     = 1.461e-6
    S     = 110.3e0
    mu0   = C * T**1.5 / (T + S)
    Re    = 25000.e0
    rho_l = 1.293e0
    p_l   = rho_l * Rgas * T
    rho_r = 0.125e0 * rho_l
    p_r   = 0.1e0 * p_l
    Lx    = Re * mu0 / np.sqrt(rho_l * p_l)
    endT  = 0.2136e0 * Lx / np.sqrt(p_l / rho_l)
    # analytical solution
    rho_ex, _, p_ex = sod_solve(x, endT,
                                rho_l=rho_l, u_l=0.0, p_l=p_l,
                                rho_r=rho_r, u_r=0.0, p_r=p_r,
                                gamma=1.4, x0=0.5*Lx)
    # L1 error normalised by mean exact density
    l1_rho = np.mean(np.abs(rho_num - rho_ex)) / np.mean(rho_ex)
    l1_p   = np.mean(np.abs(p_num   - p_ex))   / np.mean(p_ex)
    assert l1_rho < L1_TOL, (
        f"Density L1 error = {l1_rho:.3f} >= {L1_TOL}"
    )
    assert l1_p < L1_TOL, (
        f"Pressure L1 error = {l1_p:.3f} >= {L1_TOL}"
    )

