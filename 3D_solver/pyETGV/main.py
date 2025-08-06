import numpy as np
from rk import RK44
from tgv import set_grid, set_init, set_Jacobian


def main():
  # grid info
  Nx    = 66
  Ny    = 66
  Nz    = 66
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
  Np    = 100
  dt    = 0.01e0

  dir = "./data"
  x, y, z, dx, dy, dz = set_grid(Nx, Ny, Nz, Lx, Ly, Lz)
  Jacobian = set_Jacobian(Nx, Ny, Nz, dx, dy, dz)
  Q = set_init(Nx, Ny, Nz, x, y, z, gamma, Rgas, M0, rho0)
  RK44(Q, x, y, z, gamma, dt, dx, dy, dz, Jacobian, Nt, Np, dir)


main()

