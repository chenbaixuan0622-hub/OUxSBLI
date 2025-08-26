import numpy as np
import cufd
import tgv


def mod_time():
  Nt = cufd.get_Nt()
  Np = cufd.get_Np()
  dt = cufd.get_dt()
  return Nt, Np, dt


def set_grid():
  Nx = cufd.get_Nx()
  Ny = cufd.get_Ny()
  Nz = cufd.get_Nz()
  Lx = cufd.get_Lx()
  Ly = cufd.get_Ly()
  Lz = cufd.get_Lz()
  x, y, z, dx, dy, dz = tgv.set_grid(Nx, Ny, Nz, Lx, Ly, Lz)
  Jacobian = tgv.set_Jacobian(Nx, Ny, Nz, dx, dy, dz)
  return Nx, Ny, Nz, x, y, z, dx, dy, dz, Jacobian

