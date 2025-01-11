import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number, Data
from mod.mod_POD import make_data3D, make_grid3D, snapshot_pod, calc_time_coef
from mod.mod_plot import print_VTK


set_Params()

# parameter
#Q_dir   = "../3D_solver/TBL/data"#"../../z6mm"
Q_dir   = "../../../../../../media/user/HD-EDS-E/TBL/SBLI_1delta/HD"
Lx1     = 24.e-3
Lx2     = 40.e-3
Ly1     = 0.e-3
Ly2     = 10.e-3
Lz1     = 0.e-3 
Lz2     = 2.e-3
stridex = 1
stridey = 1
stridez = 1
endT    = 0.8e-3


def plot_pod_results(POD_dir, x, y, z, t, data, eigenvalues, modes, time_coefficients, num_modes):
  # energy
  plt.figure(figsize=(4, 4))
  plt.plot(range(1,11), eigenvalues[:10] / np.sum(eigenvalues) * 100, 'o-')
  plt.title('Energy Contribution of Modes')
  plt.xlabel('Mode Index')
  plt.ylabel('Energy (%)')
  save_name = "POD3D_energy_contribution.png"
  save_path = os.path.join(POD_dir, save_name)
  plt.savefig(save_path)
  plt.close()

  # space mode
  nx = len(x)
  ny = len(y)
  nz = len(z)
  for i in range(num_modes):
    u = np.reshape(modes[:,i], [nz, ny, nx])
    # write VTK file
    name = "Mode" + str(i+1)
    print_VTK(x, y, z, u, POD_dir, name)

  # time coefficient
  row = 2
  col = num_modes//2
  fig, ax = plt.subplots(col, row, figsize=(12, 8))
  for i in range(num_modes):
    ax[i%col,i//col].plot(t*1e3, time_coefficients[:,i])
    ax[i%col,i//col].set_xlabel("t [ms]")
    ax[i%col,i//col].set_xlim(t[0]*1e3, t[-1]*1e3)
    ax[i%col,i//col].set_ylabel(f'Mode {i+1}')

  save_name = "POD_time_coef"
  save_path = os.path.join(POD_dir, save_name)
  np.save(save_path, time_coefficients)

  fig.tight_layout()
  save_name = "POD_time_coef.png"
  save_path = os.path.join(POD_dir, save_name)
  plt.savefig(save_path)
  plt.close()


def plot_reconstruction(POD_dir, x, y, z, modes, time_coef, num_modes, mean=None):
  Nt = len(time_coef[:,0])
  nx = len(x)
  ny = len(y)
  nz = len(z)
  POD = np.dot(time_coef[:,:num_modes], modes[:,:num_modes].T)
  if mean is not None:
    POD = POD + mean[None,:]
  POD = np.reshape(POD, [Nt,nz,ny,nx])
  for i in range(Nt):
    name = "POD" + str(i+1).zfill(5)
    print_VTK(x, y, z, POD[i,:,:,:], POD_dir, name)


def main(Q_dir, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez, endT):
  POD_dir = os.path.join(Q_dir, "vector")
  os.makedirs(POD_dir, exist_ok=True)
  D_path = os.path.join(POD_dir, "D.npy")

  if os.path.isfile(D_path):
    data = Data(Q_dir, endT, Lx1=None, Lx2=None, Ly1=None, Ly2=None, Lz1=None, Lz2=None)
    _, _, _, _, _, _, x, y, z, t = data.make_grid()
    D = np.load(D_path)
  else:
    # D[space, time]
    x, y, z, t, D = make_data3D(Q_dir, endT, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez)
    np.save(D_path, D)
    
  mean = np.mean(D[:,120:], axis=-1)
  D[:,120:] = D[:,120:] - mean[:,None]

  eigenvalues, eigenvectors, modes = snapshot_pod(D[:,120:])
  time_coef = calc_time_coef(D[:,120:], modes)

  plot_pod_results(POD_dir, x, y, z, t[:-120], D[:,120:], eigenvalues, modes, time_coef, num_modes=10)
  #plot_reconstruction(POD_dir, x, y, z, modes, time_coef, num_modes=10, mean=mean)


main(Q_dir, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez, endT)

