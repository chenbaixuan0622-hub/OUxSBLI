import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_POD import make_data3D, snapshot_pod, calc_time_coef
from mod.mod_plot import print_VTK


set_Params()

# parameter
#Q_dir   = "../3D_solver/TBL/data"#"../../z6mm"
Q_dir   = "../../../../../../media/user/HD-EDS-E/TBL/SBLI_1delta"
Lx1     = 24.e-3
Lx2     = 40.e-3
Ly1     = 0.e-3
Ly2     = 3.e-3
Lz1     = 0.e-3 
Lz2     = 2.e-3
stridex = 4
stridey = 8
stridez = stridex * 2
endT    = 0.1e-3


def plot_pod_results(POD_dir, x, y, z, data, eigenvalues, modes, time_coefficients, num_modes=6):
  # energy
  plt.figure(figsize=(6, 4))
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
  fig, ax = plt.subplots(num_modes, 1, figsize=(8, 8))
  for i in range(num_modes):
    ax[i].plot(time_coefficients[:,i])
    
  save_name = "POD_time_coef"
  save_path = os.path.join(POD_dir, save_name)
  np.save(save_path, time_coefficients)

  save_name = "POD_time_coef.png"
  save_path = os.path.join(POD_dir, save_name)
  plt.savefig(save_path)
  plt.close()


def plot_reconstruction(POD_dir, x, y, z, mean, data, modes, time_coef, num_modes):
  Nt = len(time_coef[:,0])
  nx = len(x)
  ny = len(y)
  nz = len(z)
  POD = np.dot(time_coef[:,:num_modes], modes[:,:num_modes].T)
  POD = np.reshape(POD, [Nt,nz,ny,nx])
  for i in range(Nt):
    name = "POD" + str(i+1).zfill(5)
    print_VTK(x, y, z, POD[i,:,:,:], POD_dir, name)


POD_dir = os.path.join(Q_dir, "POD3D")
os.makedirs(POD_dir, exist_ok=True)

x, y, z, t, D = make_data3D(Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez, endT, Q_dir)
mean = np.mean(D)
D    = D - mean
eigenvalues, eigenvectors, modes = snapshot_pod(D)
time_coef = calc_time_coef(D, modes)

plot_pod_results(POD_dir, x, y, z, D, eigenvalues, modes, time_coef, num_modes=6)
plot_reconstruction(POD_dir, x, y, z, mean, D, modes, time_coef, num_modes=6)

