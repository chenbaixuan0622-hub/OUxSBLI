import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_POD import make_grid, make_data, snapshot_pod, calc_time_coef


set_Params()


# parameter
Q_dir   = "../3D_solver/TBL/data"
Lx1     = 8.e-3
Lx2     = 20.e-3
Ly1     = 0.e-3
Ly2     = 4.e-3
stridex = 8
stridey = 4
endT    = 0.1e-3
num_modes = 6

# D[space=nx*ny, time]
save_path = os.path.join(Q_dir, "D.npy")
if os.path.isfile(save_path):
  D = np.load(save_path)
  _, _, _, _, _, x, y, t = make_grid(Lx1, Lx2, Ly1, Ly2, stridex, stridey, endT, Q_dir)
else:
  x, y, t, D = make_data(Lx1, Lx2, Ly1, Ly2, stridex, stridey, endT, Q_dir)
  np.save(save_path, D)

print(len(x), len(y), np.shape(D))

mean = np.mean(D, axis=-1)
for i in range(len(D[0,:])):
  D[:,i] = D[:,i] - mean

eigenvalues, eigenvectors, modes = snapshot_pod(D)
time_coef = calc_time_coef(D, modes)

Nt = len(D[:,0])

POD  = np.dot(time_coef[:,:num_modes], modes[:,:num_modes].T)
POD  = np.reshape(POD, [Nt,nz,nx])

plt.imshow(POD[0,:,:])

