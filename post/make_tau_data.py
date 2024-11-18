import os
import re
import numpy as np
import matplotlib.pyplot as plt
import mod.mod_plot as myplt
from tqdm import tqdm
from mod.mod_read import getGrid, getScalar, getVector
from mod.mod_turb_stat import tau_2d

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240819_KEEP6thVisc4th"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

Q_path              = os.path.join(Q_directory, Q_files[50])
Nx, Ny, Nz, x, y, z = getGrid(Q_path)


data_directory = os.path.join(Q_directory, "learning_data")
os.makedirs(data_directory, exist_ok = True)

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)


um_path   = os.path.join(Q_directory, "um.npy"  )
vm_path   = os.path.join(Q_directory, "vm.npy"  )
wm_path   = os.path.join(Q_directory, "wm.npy"  )

um = np.load(um_path)
vm = np.load(vm_path)
wm = np.load(wm_path)

for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  Q[0,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'rho'     )
  Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(file_path, Nx, Ny, Nz, 'velocity')
  Q[4,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'p'       )
  # tw[nz,nx], ut[nz,nx]
  tw, ut = tau_2d(Q, x, y, z)
  # TKE
  uf = Q[1,:,:,:] - um
  vf = Q[2,:,:,:] - vm
  wf = Q[3,:,:,:] - wm

