import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number, Data
from mod.mod_POD import make_data3D, make_grid3D, snapshot_pod, calc_time_coef, plot_pod_results, plot_reconstruction
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

  plot_pod_results(POD_dir, x, y, z, t[:-120], eigenvalues, modes, time_coef, num_modes=10)
  plot_reconstruction(POD_dir, x, y, z, modes, time_coef, num_modes=10, mean=mean)


main(Q_dir, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez, endT)

