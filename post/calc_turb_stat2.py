import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland, tau
from mod.mod_TKE import BudgetTerms, print_TKE

Q_directory = "../../../../../../media/user/HD-EDS-E/hatayama/TBL/SBLI_05delta/stat03ms_09ms_SLAU"
Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

Rgas  = 287.03

def TKE(Nx, Ny, Nz, nx0, nx, ny, nz, x, y, z, Rho, U, V, W, P, UF, VF, WF, yp, ut, Norm):
  Pro = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  T   = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  PI  = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  D   = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  eps = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)

  itr = 0.e0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") \
    or file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    rhos       = getScalar(file_path, Nx, Ny, Nz, 'rho')
    us, vs, ws = getVector(file_path, Nx, Ny, Nz, 'velocity')
    ps         = getScalar(file_path, Nx, Ny, Nz, 'p')
    rho = rhos[:,:,nx0:nx0+nx]
    u   =   us[:,:,nx0:nx0+nx]
    v   =   vs[:,:,nx0:nx0+nx]
    w   =   ws[:,:,nx0:nx0+nx]
    p   =   ps[:,:,nx0:nx0+nx]
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


def main():
  first_path          = os.path.join(Q_directory, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  
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

  Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rho, u, v, w, p
  yp, _, _, _, ut, _ = non_dim_tbl(Q[:,:,:,nx1:nx2], X[nx1:nx2], Y, Z)
  save_path = os.path.join(Q_directory, "yp.npy")
  np.save(save_path, yp)
  del Q
  nu   = np.mean(Sutherland(p[:,0,nx1:nx2] / (Rgas * rho[:,0,nx1:nx2])) / (rho[:,0,nx1:nx2]))
  Norm = np.mean(rho[:,0,nx1:nx2]) * ut**4 / nu

  nx0 = 0#int(0.3 * Nx)
  nx  = Nx#int(0.6 * Nx)
  ny  = Ny
  nz  = Nz

  Rho = rho[:,:,nx0:nx0+nx]
  U   =   u[:,:,nx0:nx0+nx]
  V   =   v[:,:,nx0:nx0+nx]
  W   =   w[:,:,nx0:nx0+nx]
  P   =   p[:,:,nx0:nx0+nx]
  UF  =  uF[:,:,nx0:nx0+nx]
  VF  =  vF[:,:,nx0:nx0+nx]
  WF  =  wF[:,:,nx0:nx0+nx]

  x = X[nx0:nx0+nx]
  y = Y
  z = Z
  del rho, u, v, w, p, uF, vF, wF, X, Y, Z

  print(Norm)
  TKE(Nx, Ny, Nz, nx0, nx, ny, nz, x, y, z, Rho, U, V, W, P, UF, VF, WF, yp, ut, Norm)


main()

