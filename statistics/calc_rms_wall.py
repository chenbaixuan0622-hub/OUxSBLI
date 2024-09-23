import numpy as np
import os
import re
from numba import jit
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TOS/TOS20240830"

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

um_path   = os.path.join(Q_directory, "um.npy")
vm_path   = os.path.join(Q_directory, "vm.npy")
wm_path   = os.path.join(Q_directory, "wm.npy")
pm_path   = os.path.join(Q_directory, "pm.npy")
rhom_path = os.path.join(Q_directory, "rhom.npy")
if os.path.isfile(pm_path) & os.path.isfile(rhom_path) & os.path.isfile(um_path) & os.path.isfile(vm_path) & os.path.isfile(wm_path):
  pm   = np.load(pm_path)
  rhom = np.load(rhom_path)
  um   = np.load(um_path)
  vm   = np.load(vm_path)
  wm   = np.load(wm_path)
else:
  rhom = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  pm   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  p2m  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  um   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vm   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  wm   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    rho     = getScalar(file_path, Nx, Ny, Nz, "rho")
    p       = getScalar(file_path, Nx, Ny, Nz, 'p')
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    rhom += rho
    pm   += p
    p2m  += p**2
    um   += u
    vm   += v
    wm   += w
  rhom = rhom / float(num_files)
  pm   = pm   / float(num_files)
  p2m  = p2m  / float(num_files)
  um   = um   / float(num_files)
  vm   = vm   / float(num_files)
  wm   = wm   / float(num_files)
  np.save(rhom_path, rhom)
  np.save(pm_path,     pm)
  np.save(um_path,     um)
  np.save(vm_path,     vm)
  np.save(wm_path,     wm)

prms = np.sqrt(p2m - pm**2)
prms_path = os.path.join(Q_directory, "prms.npy")
np.save(prms_path, prms)

'''
save_path = os.path.join(Q_directory, "p_rms.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# (x-x0)/d prms", file=f)
  for i in range(Nx1, Nx2):
    print(f'{(x[i]-x[Nx1])/delta:.3e}', f'{np.mean(prms[:,0,i]):.3e}', file=f)
'''

