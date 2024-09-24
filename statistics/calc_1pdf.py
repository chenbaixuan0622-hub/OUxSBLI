import numpy as np
import os
import re
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar

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

u    = np.zeros(num_files, dtype=np.float32)
dudx = np.zeros(num_files, dtype=np.float32)

ix = int(0.5 * Nx)
iy = 50
iz = int(0.5 * Nz)

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  Q[0,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'rho')
  Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:]  = getVector(file_path, Nx, Ny, Nz, 'velocity')
  Q[4,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'p')
  u[itr]    = Q[1,iz,iy,ix]
  dudx[itr] = (-Q[1,iz,iy,ix-1] + Q[1,iz,iy,ix+1]) / (-x[ix-1] + x[ix+1])
  itr += 1

num_bins = 20

u_hist,  u_bins  = np.histogram(u,    bins=num_bins, density=True)
du_hist, du_bins = np.histogram(dudx, bins=num_bins, density=True)

save_path = os.path.join(Q_directory, "u_dudx_pdf.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# u_bins  u_hist  du_bins du_hist", file=f)
  for i in range(num_bins):
    print(f'{0.5e0 * (u_bins[i]  +  u_bins[i+1]):.3e}', f'{u_hist[i]:.3e}', \
          f'{0.5e0 * (du_bins[i] + du_bins[i+1]):.3e}', f'{du_hist[i]:.3e}', file=f)

