import numpy as np
import cupy as cp
import cufd
from tqdm import tqdm
from print import print_vtk
from tgv import set_bc


def RK44(Q, x, y, z, gamma, dt, dx, dy, dz, Jacobian, Nt, Np, dir):
  print_vtk(x, y, z, gamma, Q, 0, dir) # initial condition
  _, Nz, Ny, Nx = Q.shape
  QJ = cp.zeros(Nz*Ny*Nx*5, dtype=cp.float64)
  Q1 = cp.zeros(Nz*Ny*Nx*5, dtype=cp.float64)
  Rs = cp.zeros((Nz-2)*(Ny-2)*(Nx-2)*5, dtype=cp.float64)
  E  = cp.zeros((Nz-2)*(Ny-2)*(Nx-1)*5, dtype=cp.float64)
  F  = cp.zeros((Nz-2)*(Ny-1)*(Nx-2)*5, dtype=cp.float64)
  G  = cp.zeros((Nz-1)*(Ny-2)*(Nx-2)*5, dtype=cp.float64)
  QJ = cp.asarray(np.ravel(Q))
  xix   = cp.asarray(1.e0 / dx)
  etay  = cp.asarray(1.e0 / dy)
  zetaz = cp.asarray(1.e0 / dz)
  print(cp.cuda.runtime.runtimeGetVersion())
  print(cp.cuda.compiler.get_rocm_path())
  for itr in tqdm(range(1, Np+1)):
    for n in range(Nt):
      cufd.calc_EFG_Euler(Nx, Ny, Nz, xix, etay, zetaz, Jacobian, QJ, E, F, G)
      cufd.calc_step(Nx, Ny, Nz, 0.5e0, 1.e0, xix, etay, zetaz, E, F, G, QJ, Q1, Rs)
      set_bc(Q1, Nx, Ny, Nz)

      cufd.calc_EFG_Euler(Nx, Ny, Nz, xix, etay, zetaz, Jacobian, Q1, E, F, G)
      cufd.calc_step(Nx, Ny, Nz, 0.5e0, 2.e0, xix, etay, zetaz, E, F, G, QJ, Q1, Rs)
      set_bc(Q1, Nx, Ny, Nz)

      cufd.calc_EFG_Euler(Nx, Ny, Nz, xix, etay, zetaz, Jacobian, Q1, E, F, G)
      cufd.calc_step(Nx, Ny, Nz, 1.e0, 2.e0, xix, etay, zetaz, E, F, G, QJ, Q1, Rs)
      set_bc(Q1, Nx, Ny, Nz)

      cufd.calc_EFG_Euler(Nx, Ny, Nz, xix, etay, zetaz, Jacobian, Q1, E, F, G)
      cufd.calc_step_4(Nx, Ny, Nz, xix, etay, zetaz, E, F, G, Rs, QJ)
      set_bc(QJ, Nx, Ny, Nz)
    Q = QJ.reshape(5,Nz,Ny,Nx).get()
    print_vtk(x, y, z, gamma, Q, itr, dir)

