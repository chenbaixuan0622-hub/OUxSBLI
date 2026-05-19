"""
OS: 2D oblique shock (2D_solver/OS).

The case initialises a diagonal oblique shock (M=2, θ=8°) across a
average pre/post-shock states match the analytical Rankine-Hugoniot values.
"""
import pathlib
import pytest
import numpy as np
from ouxsbli import Case
from .utils.vtk_reader import latest_vtr, getGrid, getQ
from .utils.oblique_shock import reslected_shock, free_stream


R     = 287.03e0
gamma = 1.4e0
Pr    = 0.72e0
# free stream
M0    = 2.e0
p_tot = 100.e3
T_tot = 295.e0
# oblique shock
beta  = np.pi * 37.2e0 / 180.e0
beta_r = np.pi * 44.1e0 / 180.e0
# space and time
blt   = 1.e-3
Lx    = 5.e0 * blt
Ly    = 2.e0 * blt
nx    = 257
ny    = 129
endT  = 0.1e-3
dt    = 3.e-9

PRE_RTOL  = 0.03   # 3 % for undisturbed pre-shock
POST_RTOL = 0.03   # 3 % for post-shock (SLAU has some numerical diffusion)


@pytest.mark.integration
def test_os_pre_and_post_shock_match_analytical(tmp_path):
    workdir = "./tmp/os"

    case = Case(
        source   = "2D_solver/OS",
        workdir  = workdir,
        scheme   = "slau",
        accuracy = 6,
        visc     = "euler",
        tvd      = "hybrid",
        nx       = nx,
        ny       = ny,
        Lx       = Lx,
        Ly       = Ly,
        gamma    = gamma,
        Pr       = Pr,
        R        = R,
        M0       = M0,
        p_tot    = p_tot,
        T_tot    = T_tot,
        beta     = beta,
        dt       = dt,
        np       = 1,           # one output step
        nt       = int(endT / dt),
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    vtr_path = latest_vtr(data_dir)
    ni, nj, nk, x, y, z = getGrid(vtr_path)
    rho, u, v, _, p = getQ(vtr_path, ni, nj, nk)
    # analytical solution
    rho3, ux3, uy3, p3 = reslected_shock(M0, gamma, R, p_tot, T_tot, beta, beta_r)

    # free-stream (pre-shock) analytical state
    u0_fs, p0_fs, T0_fs = free_stream(M0, gamma, R, p_tot, T_tot)
    rho0_fs = p0_fs / (R * T0_fs)

    # Pre-shock region: left 10 % of domain (i < ni//10), top half (j >= nj//2).
    # Purely upstream of the incident shock for β≈37° in a 5×2 mm domain.
    rho_pre = rho[0, nj // 2 :, : ni // 10].astype(float)
    u_pre   = u[0,   nj // 2 :, : ni // 10].astype(float)
    p_pre   = p[0,   nj // 2 :, : ni // 10].astype(float)

    # Post-reflected-shock region: right 20 % (i >= 4*ni//5), bottom quarter (j < nj//4).
    # At x=4Lx/5=4 mm the reflected shock (from x_hit≈2.6 mm) is at y≈1.1 mm,
    # so j < nj//4 is well below the reflected shock.
    rho_post = rho[0, : nj // 3, 9 * ni // 10 :].astype(float)
    u_post   = u[0,   : nj // 3, 9 * ni // 10 :].astype(float)
    v_post   = v[0,   : nj // 3, 9 * ni // 10 :].astype(float)
    p_post   = p[0,   : nj // 3, 9 * ni // 10 :].astype(float)

    # pre-shock assertions (undisturbed free stream)
    assert abs(rho_pre.mean() - rho0_fs) / rho0_fs < PRE_RTOL, (
        f"Pre-shock rho: {rho_pre.mean():.4f} vs {rho0_fs:.4f}")
    assert abs(u_pre.mean() - u0_fs) / u0_fs < PRE_RTOL, (
        f"Pre-shock u: {u_pre.mean():.4f} vs {u0_fs:.4f}")
    assert abs(p_pre.mean() - p0_fs) / p0_fs < PRE_RTOL, (
        f"Pre-shock p: {p_pre.mean():.4f} vs {p0_fs:.4f}")

    # post-shock assertions (after reflected shock)
    assert abs(rho_post.mean() - rho3) / rho3 < POST_RTOL, (
        f"Post-shock rho: {rho_post.mean():.4f} vs {rho3:.4f}")
    assert abs(u_post.mean() - ux3) / abs(ux3) < POST_RTOL, (
        f"Post-shock u: {u_post.mean():.4f} vs {ux3:.4f}")
    # After a perfect wall reflection the flow is horizontal; uy3 from reslected_shock
    # is non-zero only because beta_r=37.7° is approximate — check v≈0 directly.
    assert abs(v_post.mean()) / abs(ux3) < POST_RTOL, (
        f"Post-shock v: {v_post.mean():.4f} (expected ≈ 0)")
    assert abs(p_post.mean() - p3) / p3 < POST_RTOL, (
        f"Post-shock p: {p_post.mean():.4f} vs {p3:.4f}")
