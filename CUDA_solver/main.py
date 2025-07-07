import numpy as np
from rk import RK33, RK44
from tgv import set_grid, set_init


def main():
  # grid info
  Nx    = 65
  Ny    = 65
  Nz    = 65
  Lx    = 2.e0 * np.pi
  Ly    = 2.e0 * np.pi
  Lz    = 2.e0 * np.pi

  # physical properties
  gamma = 1.4e0
  Rgas  = 287.03e0
  M0    = 0.4e0
  rho0  = 1.e0 

  # time
  Nt    = 200
  Np    = 10#100
  dt    = 0.01e0

  dir = "./data"
  x, y, z, dx, dy, dz = set_grid(Nx, Ny, Nz, Lx, Ly, Lz)
  Q = set_init(Nx, Ny, Nz, x, y, z, gamma, Rgas, M0, rho0)

  #RK33(Q, x, y, z, gamma, dt, dx, dy, dz, Nt, Np, dir)
  RK44(Q, x, y, z, gamma, dt, dx, dy, dz, Nt, Np, dir)


main()

