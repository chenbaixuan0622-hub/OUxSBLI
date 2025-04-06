import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland, tau
from mod.mod_TKE import BudgetTerms, print_TKE


Q_directory = "../../SBLI/SBLI_16delta/2.2ms/5"
Q_files = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

Rgas  = 287.03

# wall unit
rhow = 0.1886e0
ut   = 23.33e0
mu   = 1.74e-5


def TKE(Nx, Ny, Nz, x, y, z, Rho, U, V, W, P, UF, VF, WF, Norm):
  nz, ny, nx = Rho.shape
  Pro = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  T   = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  PI  = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  D   = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  eps = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)

  itr = 0.e0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") \
    or file_path == os.path.join(Q_directory, "ReynoldsStress.vtr") \
    or file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    p       = getScalar(file_path, Nx, Ny, Nz, 'p')
    # TKE Budget
    BudgetTerms(x, y, z, rho, u, v, w, p, Rho, U, V, W, P, UF, VF, WF, Pro, T, PI, D, eps)
    itr += 1.e0
  # normalize
  Pro = Pro / (Norm * itr)
  T   = T   / (Norm * itr)
  PI  = PI  / (Norm * itr)
  D   = D   / (Norm * itr)
  eps = eps / (Norm * itr)
  print_TKE(x[1:nx-1], y[1:ny-1], z[1:nz-1], Pro, T, PI, D, eps, Q_directory, "TKE")


def main(rhow=None, ut=None, mu=None, tau=None):
  first_path          = os.path.join(Q_directory, Q_files[0])
  Nx, Ny, Nz, x, y, z = getGrid(first_path)
  
  rho_path = os.path.join(Q_directory, "rho.npy")
  u_path   = os.path.join(Q_directory, "u.npy")
  v_path   = os.path.join(Q_directory, "v.npy")
  w_path   = os.path.join(Q_directory, "w.npy")
  p_path   = os.path.join(Q_directory, "p.npy")
  rho = np.load(rho_path)
  u   = np.load(u_path)
  v   = np.load(v_path)
  w   = np.load(w_path)
  p   = np.load(p_path)

  uF_path = os.path.join(Q_directory, "uF.npy")
  vF_path = os.path.join(Q_directory, "vF.npy")
  wF_path = os.path.join(Q_directory, "wF.npy")
  uF = np.load(uF_path)
  vF = np.load(vF_path)
  wF = np.load(wF_path)

  nx1 = int(0.05 * Nx)
  nx2 = int(0.25 * Nx)

  if rhow is not None and ut is not None and mu is not None:
    Norm = rhow**2 * ut**4 / mu
  else:
    Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
    Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rho, u, v, w, p
    yp, _, _, _, ut, _ = non_dim_tbl(Q[:,:,:,nx1:nx2], X[nx1:nx2], y, z)
    nu   = np.mean(Sutherland(p[:,0,nx1:nx2] / (Rgas * rho[:,0,nx1:nx2])) / (rho[:,0,nx1:nx2]))
    Norm = np.mean(rho[:,0,nx1:nx2]) * ut**4 / nu
    del Q
  TKE(Nx, Ny, Nz, x, y, z, rho, u, v, w, p, uF, vF, wF, Norm)


main(rhow, ut, mu)

