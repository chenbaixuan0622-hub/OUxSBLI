"""
ETGV: verify that the KEEP scheme on the inviscid Taylor-Green vortex

The test builds 3D_solver/ETGV, runs 100 timesteps (nt=5 intervals of np=20
steps), and reads the text diagnostics that the solver appends each interval:
    data/entropy.d          columns: [t, (s0-s)/s0]
    data/kinetic_energy.d   columns: [t, ke, ke/ke0]
"""
import pathlib
import pytest
import numpy as np
from ouxsbli import Case


ENTROPY_TOL = 3.e-3   # |(s0-s)/s0| threshold
KE_RES      = 0.75e0  # |ke/ke0|    residual


@pytest.mark.integration
def test_etgv_ke_and_entropy_preserved(tmp_path):
    workdir = "./tmp/etgv"
    case = Case(
        source   = "3D_solver/ETGV",
        workdir  = workdir,
        visc     = "euler",
        scheme   = "keep",
        accuracy = 2,
        rk       = "rk4",
        nx       = 66,
        ny       = 66,
        nz       = 66,
        nt       = 20000,
        np       = 1,
        dt       = 0.01e0,
    )
    case.build()
    case.run(nranks=2)
    data_dir = pathlib.Path(workdir) / "data"
    
    # --- entropy ---
    entropy_file = data_dir / "entropy.d"
    assert entropy_file.exists(), "entropy.d was not produced"
    ent_data  = np.loadtxt(entropy_file)
    ent_final = ent_data[-1, -1]
    assert np.abs(ent_final) < ENTROPY_TOL, (
        f"Entropy not preserved: max |(s0-s)/s0| = {ent_final:.2e} >= {ENTROPY_TOL}"
    )

    # --- kinetic energy ---
    ke_file = data_dir / "kinetic_energy.d"
    assert ke_file.exists(), "kinetic_energy.d was not produced"
    ke_data  = np.loadtxt(ke_file)
    ke_final = ke_data[-1, -1]
    assert ke_final > KE_RES, (
        f"Kinetic energy not preserved: ke/ke0 = {ke_final:.2e} < {KE_RES}"
    )

