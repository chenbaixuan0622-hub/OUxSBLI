import numpy as np
import cupy as cp
import cufd
from tqdm import tqdm
from print import print_vtk
from tgv import set_bc


def RK33(Q, x, y, z, gamma, dt, dx, dy, dz, Nt, Np, dir):
  print_vtk(x, y, z, gamma, Q, 0, dir) # initial condition
  Nz, Ny, Nx, _ = Q.shape
  q  = cp.zeros(Nz*Ny*Nx*5, dtype=cp.float64)
  q1 = cp.zeros(Nz*Ny*Nx*5, dtype=cp.float64)
  e  = cp.zeros((Nz-2)*(Ny-2)*(Nx-1)*5, dtype=cp.float64)
  f  = cp.zeros((Nz-2)*(Ny-1)*(Nx-2)*5, dtype=cp.float64)
  g  = cp.zeros((Nz-1)*(Ny-2)*(Nx-2)*5, dtype=cp.float64)
  q  = cp.asarray(np.ravel(Q))
  for itr in tqdm(range(1, Np+1)):
    for n in range(Nt):
      cufd.calc_flux(q, e, f, g, gamma, Nx, Ny, Nz)
      cufd.RK33_1(q, e, f, g, q1, dt, dx, dy, dz, Nx, Ny, Nz)
      set_bc(q1, Nx, Ny, Nz)
  
      cufd.calc_flux(q1, e, f, g, gamma, Nx, Ny, Nz)
      cufd.RK33_2(q, e, f, g, q1, dt, dx, dy, dz, Nx, Ny, Nz)
      set_bc(q1, Nx, Ny, Nz)
  
      cufd.calc_flux(q1, e, f, g, gamma, Nx, Ny, Nz)
      cufd.RK33_3(q1, e, f, g, q, dt, dx, dy, dz, Nx, Ny, Nz)
      set_bc(q, Nx, Ny, Nz)
    Q = q.reshape(Nz,Ny,Nx,5).get()
    print_vtk(x, y, z, gamma, Q, itr, dir)


def RK44(Q, x, y, z, gamma, dt, dx, dy, dz, Nt, Np, dir):
  print_vtk(x, y, z, gamma, Q, 0, dir) # initial condition
  Nz, Ny, Nx, _ = Q.shape
  q  = cp.zeros(Nz*Ny*Nx*5, dtype=cp.float64)
  q1 = cp.zeros(Nz*Ny*Nx*5, dtype=cp.float64)
  r1 = cp.zeros((Nz-2)*(Ny-2)*(Nx-2)*5, dtype=cp.float64)
  e  = cp.zeros((Nz-2)*(Ny-2)*(Nx-1)*5, dtype=cp.float64)
  f  = cp.zeros((Nz-2)*(Ny-1)*(Nx-2)*5, dtype=cp.float64)
  g  = cp.zeros((Nz-1)*(Ny-2)*(Nx-2)*5, dtype=cp.float64)
  q  = cp.asarray(np.ravel(Q))
  for itr in tqdm(range(1, Np+1)):
    for n in range(Nt):
      cufd.calc_flux(q, e, f, g, gamma, Nx, Ny, Nz)
      cufd.RK44_1(q, e, f, g, q1, r1, dt, dx, dy, dz, Nx, Ny, Nz)
      set_bc(q1, Nx, Ny, Nz)

      cufd.calc_flux(q1, e, f, g, gamma, Nx, Ny, Nz)
      cufd.RK44_2(q, e, f, g, q1, r1, dt, dx, dy, dz, Nx, Ny, Nz)
      set_bc(q1, Nx, Ny, Nz)

      cufd.calc_flux(q1, e, f, g, gamma, Nx, Ny, Nz)
      cufd.RK44_3(q, e, f, g, q1, r1, dt, dx, dy, dz, Nx, Ny, Nz)
      set_bc(q1, Nx, Ny, Nz)

      cufd.calc_flux(q1, e, f, g, gamma, Nx, Ny, Nz)
      cufd.RK44_4(q1, e, f, g, q, r1, dt, dx, dy, dz, Nx, Ny, Nz)
      set_bc(q, Nx, Ny, Nz)
    Q = q.reshape(Nz,Ny,Nx,5).get()
    print_vtk(x, y, z, gamma, Q, itr, dir)

