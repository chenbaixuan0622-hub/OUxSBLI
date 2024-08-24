import numpy as np
import os
import re
from numba import jit
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL20240813"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

Nx1 = 
Nx2 = 

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

pm_path  = os.path.join(Q_directory, "pm.npy")
p2m_path = os.path.join(Q_directory, "p2m.npy")
if os.path.isfile(pm_path) & os.path.isfile(p2m_path):
  pm  = np.load(pm_path)
  p2m = np.load(p2m_path)
else:
  pm  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  p2m = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    p = getScalar(file_path, Nx, Ny, Nz, 'p')
    pm  += p
    p2m += p**2
  pm  = pm  / float(num_files)
  p2m = p2m / float(num_files)
  np.save(pm_path,  pm)
  np.save(p2m_path, p2m)

prms = np.sqrt(p2m - pm**2)

save_path = os.path.join(Q_directory, "p_rms.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# (x-x0)/d prms", file=f)
  for i in range(Nx1, Nx2):
    print(f'{(x[i]-x[Nx1])/delta:.3e}', f'{np.mean(prms[:,0,i]):.3e}', file=f)

