import numpy as np
from rk import RK44
from tgv import set_grid, set_init, set_Jacobian


def main():
  # grid info
  Nx    = 130
  Ny    = 130
  Nz    = 130
  L0    = 1.524e-3
  Lx    = 2.e0 * np.pi * L0
  Ly    = 2.e0 * np.pi * L0
  Lz    = 2.e0 * np.pi * L0

  # physical properties
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

  # time
  CFL   = 0.03e0
  dt    = CFL * (Lx / (Nx-1)) / v0
  dtn   = v0 * dt / L0
  Np    = 100
  Nt    = int(20.e0 / (Np * dtn))

  dir = "./data"
  x, y, z, dx, dy, dz = set_grid(Nx, Ny, Nz, Lx, Ly, Lz)
  Jacobian = set_Jacobian(Nx, Ny, Nz, dx, dy, dz)
  Q = set_init(Nx, Ny, Nz, x, y, z, gamma, Rgas, rho0, v0, p0, T, L0)
  RK44(Q, x, y, z, gamma, dt, dx, dy, dz, Jacobian, Nt, Np, dir)


main()

