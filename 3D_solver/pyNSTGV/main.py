import numpy as np
import time
import cufd
import os
import sys
sys.path.append(os.pardir)
from mod import set_grid, mod_time
from rk import RK44
from tgv import set_init, set_Jacobian


def main():
  Nx, Ny, Nz, x, y, z, dx, dy, dz, Jacobian = set_grid()
  Nt, Np, dt = mod_time()
  # physical properties
  L0    = 1.524e-3
  gamma = 1.4e0
  Rgas  = 287.03e0
  M0    = 0.1e0
  Re    = 1600.e0
  T     = 530.e0 * 5.e0 / 9.e0 
  S     = 111.e0
  mu0   = 1.716e-5 * (273.2e0 + S) / (T + S) * (T / 273.2e0)**1.5e0
  v0    = M0 * np.sqrt(gamma * Rgas * T)
  rho0  = mu0 * Re / (v0 * L0)
  p0    = rho0 * Rgas * T

  dir = "./data"
  Q = set_init(Nx, Ny, Nz, x, y, z, gamma, Rgas, rho0, v0, p0, T, L0)
  t1 = time.time()
  RK44(Q, x, y, z, gamma, dt, dx, dy, dz, Jacobian, Nt, Np, dir)
  t2 = time.time()
  print("Elapsed time:", t2-t1)


main()

