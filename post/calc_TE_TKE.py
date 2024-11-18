import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_TKE import BudgetTerms
from mod.mod_info import TE, NTE, cross_corr

Q_directory = "../../a100"
Q_files     = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
x_target    = 20e-3
yp_target   = 10

def TKE(Nx, Ny, Nz, nx0, nx, ny0, ny, nz, x, y, z, Rho, U, V, W, P, UF, VF, WF):
  Pro = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  T   = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  PI  = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  D   = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  eps = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)

  num_files = len(Q_files)
  value = np.zeros((5,num_files), dtype=np.float32)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr"):
      continue
    rhos       = getScalar(file_path, Nx, Ny, Nz, 'rho')
    us, vs, ws = getVector(file_path, Nx, Ny, Nz, 'velocity')
    ps         = getScalar(file_path, Nx, Ny, Nz, 'p')
    rho = rhos[:,ny0:ny0+ny,nx0:nx0+nx]
    u   =   us[:,ny0:ny0+ny,nx0:nx0+nx]
    v   =   vs[:,ny0:ny0+ny,nx0:nx0+nx]
    w   =   ws[:,ny0:ny0+ny,nx0:nx0+nx]
    p   =   ps[:,ny0:ny0+ny,nx0:nx0+nx]
    # TKE Budget
    BudgetTerms(x, y, z, rho, u, v, w, p, Rho, U, V, W, P, UF, VF, WF, Pro, T, PI, D, eps)
    value[0,itr] = np.mean(Pro)
    value[1,itr] = np.mean(T)
    value[2,itr] = np.mean(PI)
    value[3,itr] = np.mean(D)
    value[4,itr] = np.mean(eps)
    itr += 1
  del Pro, T, PI, D, eps
 
  # calc corr
  corr_map = np.zeros((5,5), dtype=np.float32)
  for j in range(5):
    for i in range(5):
      if i > j:
        corr_map[i,j] = cross_corr(value[i,:], value[j,:])
        corr_map[j,i] = corr_map[i,j]

  # calc TE
  causal_map = np.zeros((5,5), dtype=np.float32)
  for j in range(5):
    for i in range(5):
      if i > j:
        causal_map[i,j], causal_map[j,i] = TE(value[i,:], value[j,:], 10, 5)
  print(corr_map)
  print(causal_map)


def main(x_target, yp_target):
  first_path          = os.path.join(Q_directory, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  
  rho_path = os.path.join(Q_directory, "rho.npy")
  u_path   = os.path.join(Q_directory, "u.npy")
  v_path   = os.path.join(Q_directory, "v.npy")
  w_path   = os.path.join(Q_directory, "w.npy")
  p_path   = os.path.join(Q_directory, "p.npy")
  uF_path  = os.path.join(Q_directory, "uF.npy")
  vF_path  = os.path.join(Q_directory, "vF.npy")
  wF_path  = os.path.join(Q_directory, "wF.npy")
  rho = np.load(rho_path)
  u   = np.load(u_path)
  v   = np.load(v_path)
  w   = np.load(w_path)
  p   = np.load(p_path)
  uF  = np.load(uF_path)
  vF  = np.load(vF_path)
  wF  = np.load(wF_path)

  yp_path = os.path.join(Q_directory, "yp.npy")
  yp      = np.load(yp_path)

  for i in range(Nx):
    if X[i] > x_target:
      nx0 = i-1
      break

  for j in range(Ny):
    if yp[j] > yp_target:
      ny0 = j-1
      break

  nx  = 5
  ny  = 5
  nz  = Nz

  Rho = rho[:,ny0:ny0+ny,nx0:nx0+nx]
  U   =   u[:,ny0:ny0+ny,nx0:nx0+nx]
  V   =   v[:,ny0:ny0+ny,nx0:nx0+nx]
  W   =   w[:,ny0:ny0+ny,nx0:nx0+nx]
  P   =   p[:,ny0:ny0+ny,nx0:nx0+nx]
  UF  =  uF[:,ny0:ny0+ny,nx0:nx0+nx]
  VF  =  vF[:,ny0:ny0+ny,nx0:nx0+nx]
  WF  =  wF[:,ny0:ny0+ny,nx0:nx0+nx]

  x = X[nx0:nx0+nx]
  y = Y[ny0:ny0+ny]
  z = Z
  del rho, u, v, w, p, uF, vF, wF, X, Y, Z

  TKE(Nx, Ny, Nz, nx0, nx, ny0, ny, nz, x, y, z, Rho, U, V, W, P, UF, VF, WF)


main(x_target, yp_target)

