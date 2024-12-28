import numpy as np
from networkx import from_numpy_array, draw
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_POD import make_data, make_grid
from mod_AI.gnn import plot_graph_from_array
from mod.mod_info import TE


Q_dir = "../3D_solver/TBL/data"
Lx1   = 28.e-3
Lx2   = 40.e-3
Ly2   = 8.e-3
endT  = 0.1e-3
stridex = 8
stridey = 16

# D[space=nx*ny, time]
save_path = os.path.join(Q_dir, "D.npy")
if os.path.isfile(save_path):
  D = np.load(save_path)
  stridex = 8
  stridey = 16
  _, _, _, _, _, x, y, t = make_grid(Lx1, Lx2, Ly2, stridex, stridey, endT, Q_dir)
else:
  x, y, t, D = make_data(Lx1, Lx2, Ly2, stridex, stridey, endT, Q_dir)
  np.save(save_path, D)

nx = len(x)
ny = len(y)
Nt = len(t)
print("nx ", nx, " ny ", ny, " Nt ", Nt)

# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12


x, y = np.meshgrid(x*1e3, y*1e3)

save_path = os.path.join(Q_dir, "TE_matrix.npy")
A = np.load(save_path)

threshold = 0.12
A = np.where(A < threshold, 0, A)

plt.figure(figsize=(8,4))
plot_graph_from_array(A, x, y, num_edges=100)
save_path = os.path.join(Q_dir, "TE_causal_graph.png")
plt.savefig(save_path)
plt.close()

