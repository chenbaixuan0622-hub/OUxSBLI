import numpy as np
import os
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_matrix import Poisson_cyclic
from mod.mod_nabla import Vector
from mod.mod_Helmholtz import save_Helmholtz_decomposition


set_Params()

# parameter
Q_dir   = "../3D_solver/TBL/data"
Lx1     = 24.e-3#4.e-3
Lx2     = 40.e-3#20.e-3
Ly1     = 0.e-3
Ly2     = 10.e-3
stridex = 4
stridez = 4
endT    = 0.1e-3


def main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, endT):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  
  Nx, Ny, Nz, X, Y, Z = getGrid(os.path.join(Q_dir, Q_files[0]))
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  indicesx = np.arange(nx1,  nx2, stridex)
  indicesy = np.arange(0, len(Y), 1)
  indicesz = np.arange(0, len(Z), stridez)
  x = X[indicesx]
  z = Z[indicesz]

  indicesx, indicesy, indicesz = np.meshgrid(indicesx, indicesy, indicesz)
  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = U[indicesz,indicesy,indicesx].transpose(2,0,1)
    v = V[indicesz,indicesy,indicesx].transpose(2,0,1)
    w = W[indicesz,indicesy,indicesx].transpose(2,0,1)

    y  = np.linspace(Ly1, Ly2, 100, dtype=np.float32)
    Ui = np.zeros((3,len(z),len(y),len(x)), dtype=np.float32)

    for k in range(len(z)):
      for i in range(len(x)):
        Ui[0,k,:,i] = np.interp(y, Y, u[k,:,i])
        Ui[1,k,:,i] = np.interp(y, Y, v[k,:,i])
        Ui[2,k,:,i] = np.interp(y, Y, w[k,:,i])
 
    div = Vector(Ui, x, y, z).divergence()
    phi = Poisson_cyclic(x, y, z, div)
    name = "HD" + str(itr).zfill(5)
    HD_dir = os.path.join(Q_dir, "HD")
    os.makedirs(HD_dir, exist_ok=True)
    save_Helmholtz_decomposition(Ui, phi, x, y, z, HD_dir, name)
    itr += 1


main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, endT)

