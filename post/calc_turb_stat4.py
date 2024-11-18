import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl
from mod.mod_ReynoldsStress import ReynoldsStress, print_ReynoldsStress


Q_directory = "../../../../../../media/user/HD-EDS-E/hatayama/TBL/SBLI_05delta/stat03ms_09ms_SLAU"
#Q_directory = "../../SBLI_1delta/stat03ms_05ms_SLAU"
Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]


def RS(Nx, Ny, Nz, x, y, z, Rho, U, V, W, rhow, ut):
  ruu = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  rvv = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  rww = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  ruv = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  itr = 0.e0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") \
    or file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    p       = getScalar(file_path, Nx, Ny, Nz, 'p')
    uf = u - U
    vf = v - V
    wf = w - W
    # calc Reynolds Stress
    ReynoldsStress(Nx, Ny, Nz, Rho, uf, vf, wf, rhow, ut, ruu, rvv, rww, ruv)
    itr += 1.e0
  ruu /= itr
  rvv /= itr
  rww /= itr
  ruv /= itr
  print_ReynoldsStress(x, y, z, ruu, rvv, rww, ruv, Q_directory, "ReynoldsStress")

def main():
  first_path          = os.path.join(Q_directory, Q_files[0])
  Nx, Ny, Nz, x, y, z = getGrid(first_path)
  
  rho_path = os.path.join(Q_directory, "rho.npy")
  u_path   = os.path.join(Q_directory, "u.npy")
  v_path   = os.path.join(Q_directory, "v.npy")
  w_path   = os.path.join(Q_directory, "w.npy")
  p_path   = os.path.join(Q_directory, "p.npy")
  Rho  = np.load(rho_path)
  U    = np.load(u_path)
  V    = np.load(v_path)
  W    = np.load(w_path)
  P    = np.load(p_path)

  Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = Rho, U, V, W, P
  yp, _, _, _, ut, _ = non_dim_tbl(Q, x, y, z)
  save_path = os.path.join(Q_directory, "yp.npy")
  np.save(save_path, yp)
  del Q

  RS(Nx, Ny, Nz, x, y, z, Rho, U, V, W, np.mean(Rho[:,0,:]), ut)

main()

