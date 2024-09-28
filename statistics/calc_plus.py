import numpy as np
import os
import re
import matplotlib.pyplot as plt
from numba import jit
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl

Q_directory = "../3D_solver/TBL/data"#"../../../../../media/user/HD-EDS-E/hatayama/TBL/yp1/coarse/HRSLAU2MUSCL4th"
target_path = os.path.join(Q_directory, "Q00271.vtr")
save_path   = os.path.join(Q_directory, "yplus.d")
hist_path   = os.path.join(Q_directory, "hist.png")

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

for Q_file in Q_files:
  file_path = os.path.join(Q_directory, Q_file)
  if file_path == target_path:
    Q[0,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'rho')
    Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:]  = getVector(file_path, Nx, Ny, Nz, 'velocity')
    Q[4,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'p')
    yp, tw, ut, twm, utm, up = non_dim_tbl(Q, x, y, z)
    print("tw is", twm)
    print("ut is", utm)
    with open(save_path, "w", encoding="UTF-8") as f:
      print("# yplus  uplus", file=f)
      for j in range(Ny):
        print(f'{yp[j]:.3e}', f'{up[j]:.3e}', file=f)
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.size'] = 16
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    fig, axs = plt.subplots(1, 2, figsize=(8,4))
    axs[0].hist(tw.ravel(), bins=30, density=True, color='blue')
    axs[0].set_xlabel(r'$\it{\tau_{w}}$', fontsize=18)
    axs[0].set_ylabel(r'$\it{Frequency}$', fontsize=18)

    axs[1].hist(ut.ravel(), bins=30, density=True, color='blue')
    axs[1].set_xlabel(r'$\it{u_{\tau}}$', fontsize=18)
    axs[1].set_ylabel(r'$\it{Frequency}$', fontsize=18)

    plt.tight_layout()
    plt.savefig(hist_path)

