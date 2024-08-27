import numpy as np
import os
import re
from numba import jit
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl
import mod.mod_plot as myplt

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL20240813"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

Qpm_path  = os.path.join(Q_directory, "Qpm.npy")
Qp2m_path = os.path.join(Q_directory, "Qp2m.npy")
yp_path   = os.path.join(Q_directory, "yp.npy")
if os.path.isfile(Qpm_path) & os.path.isfile(Qp2m_path) & os.path.isfile(yp_path):
  Qpm  = np.load(Qpm_path)
  Qp2m = np.load(Qp2m_path)
  yps  = np.load(yp_path)
else:
  Qpm  = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  Qp2m = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  yps  = np.zeros(Ny, dtype=np.float32)
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    Q[0,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'rho')
    Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:]  = getVector(file_path, Nx, Ny, Nz, 'velocity')
    Q[4,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'p')
    yp, tw, Qp, Qvd = non_dim_tbl(Q, x, y, z)
    Qpm  += Qp
    Qp2m += Qp**2
    yps  += yp
  Qpm  = Qpm  / float(num_files)
  Qp2m = Qp2m / float(num_files)
  yps  = yps  / float(num_files)
  np.save(Qpm_path, Qpm)
  np.save(Qp2m_path, Qp2m)
  np.save(yp_path, yps)

Qrms = np.sqrt(Qp2m - Qpm**2)

save_path = os.path.join(Q_directory, "u_v_rms.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# y       urms      vrms", file=f)
  for j in range(Ny):
    print(f'{y[j]/delta:.3e}', f'{yps[j]:.3e}', f'{np.mean(Qrms[1,:,j,:]):.3e}', f'{np.mean(Qrms[2,:,j,:]):.3e}', file=f)

