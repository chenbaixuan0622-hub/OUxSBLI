import numpy as np
import os
import re
from numba import jit
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl
import mod.mod_plot as myplt

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240826_KEEP4thVisc2nd"

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

um_path   = os.path.join(Q_directory, "um.npy"  )
vm_path   = os.path.join(Q_directory, "vm.npy"  )
wm_path   = os.path.join(Q_directory, "wm.npy"  )
rhom_path = os.path.join(Q_directory, "rhom.npy")
pm_path   = os.path.join(Q_directory, "pm.npy"  )

if os.path.isfile(Qpm_path) & os.path.isfile(Qp2m_path) & os.path.isfile(yp_path):
  Qpm  = np.load(Qpm_path )
  Qp2m = np.load(Qp2m_path)
  yps  = np.load(yp_path  )
  um   = np.load(um_path  )
  vm   = np.load(vm_path  )
  wm   = np.load(wm_path  )
  rhom = np.load(rhom_path)
  pm   = np.load(pm_path  )
else:
  Qpm  = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  Qp2m = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  um   = np.zeros((Nz,Ny,Nx),   dtype=np.float32)
  vm   = np.zeros((Nz,Ny,Nx),   dtype=np.float32)
  wm   = np.zeros((Nz,Ny,Nx),   dtype=np.float32)
  pm   = np.zeros((Nz,Ny,Nx),   dtype=np.float32)
  rhom = np.zeros((Nz,Ny,Nx),   dtype=np.float32)
  yps  = np.zeros(Ny, dtype=np.float32)
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    Q[0,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'rho')
    Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:]  = getVector(file_path, Nx, Ny, Nz, 'velocity')
    Q[4,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'p')
    yp, _, _, Qp, Qvd = non_dim_tbl(Q, x, y, z)
    Qpm  += Qp
    Qp2m += Qp**2
    yps  += yp
    rhom += Q[0,:,:,:]
    um   += Q[1,:,:,:]
    vm   += Q[2,:,:,:] 
    wm   += Q[3,:,:,:]
    pm   += Q[4,:,:,:]
  Qpm  = Qpm  / float(num_files)
  Qp2m = Qp2m / float(num_files)
  yps  = yps  / float(num_files)
  um   = um   / float(num_files)
  vm   = vm   / float(num_files)
  wm   = wm   / float(num_files)
  rhom = rhom / float(num_files)
  pm   = pm   / float(num_files)
  np.save(Qpm_path,  Qpm )
  np.save(Qp2m_path, Qp2m)
  np.save(yp_path,   yps )
  np.save(um_path,   um  )
  np.save(vm_path,   vm  )
  np.save(wm_path,   wm  )
  np.save(rhom_path, rhom)
  np.save(pm_path,   pm  )

Qrms = np.sqrt(Qp2m - Qpm**2)

save_path = os.path.join(Q_directory, "u_v_rms.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# y       urms      vrms", file=f)
  for j in range(Ny):
    print(f'{y[j]/delta:.3e}', f'{yps[j]:.3e}', f'{np.mean(Qrms[1,:,j,:]):.3e}', f'{np.mean(Qrms[2,:,j,:]):.3e}', file=f)

