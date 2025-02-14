import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getScalar, getVector, extract_number
from mod.mod_POD import snapshot_pod, calc_time_coef, plot_pod_results, plot_reconstruction
from mod.mod_shannon_nn import embedding_entropy_surrogate as embedding_entropy
from mod.mod_shannon_nn import transfer_entropy_surrogate as transfer_entropy
from mod.mod_plot import print_scalar


set_Params()

# parameter
dir   = "../../SBLI"
Q_dir = [os.path.join(dir, "9"), os.path.join(dir, "11")]
Lx1   = 29.e-3
Lx2   = 36.e-3
Ly1   = 0.e-3
Ly2   = 7.e-3
Lz1   = 0.e-3 
Lz2   = 32.e-3
endT  = 0.5e-3


def make_div_separation(Q_dir, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2):
  files1 = [f for f in os.listdir(Q_dir[0]) if f.endswith(".vtr")]
  files2 = [f for f in os.listdir(Q_dir[1]) if f.endswith(".vtr")]
  files1.sort(key=extract_number)
  files2.sort(key=extract_number)
  Nx1, Ny, Nz, X1, Y, Z = getGrid(os.path.join(Q_dir[0], files1[0]))
  Nx2, Ny, Nz, X2, Y, Z = getGrid(os.path.join(Q_dir[1], files2[0]))
  for i in range(Nx1):
    if Lx1 < X1[i]:
      nx1 = i
      break
  for i in range(Nx2):
    if Lx2 < X2[i]:
      nx2 = i
      break
  for j in range(Ny):
    if Ly2 < Y[j]:
      ny = j
      break
  y = Y[:ny]
  itr = 0
  for file in tqdm(files1):
    file_path1 = os.path.join(Q_dir[0], file)
    file_path2 = os.path.join(Q_dir[1], file)
    U1, V1, W1 = getVector(file_path1, Nx1, Ny, Nz, 'velocity')
    U2, V2, W2 = getVector(file_path2, Nx2, Ny, Nz, 'velocity')
    U = np.concatenate([U1[:,:ny,nx1:-3], U2[:,:ny,3:nx2]], axis=-1)
    V = np.concatenate([V1[:,:ny,nx1:-3], V2[:,:ny,3:nx2]], axis=-1)
    W = np.concatenate([W1[:,:ny,nx1:-3], W2[:,:ny,3:nx2]], axis=-1)
    x = np.linspace(X1[nx1], X2[nx2], U.shape[-1])
    # shock
    div = (-U[1:-1,1:-1,:-2] + U[1:-1,1:-1,2:]) / (-x[None,None,:-2] + x[None,None,2:]) \
        + (-V[1:-1,:-2,1:-1] + V[1:-1,2:,1:-1]) / (-y[None,:-2,None] + y[None,2:,None]) \
        + (-W[:-2,1:-1,1:-1] + W[2:,1:-1,1:-1]) / (-Z[:-2,None,None] + Z[2:,None,None])
    z = Z[1:-1]
    # masked shock
    div = np.where(div>-1e5, 0, div)
    # masked velocity
    u  = np.where(U > 0, 0, U)
    for j in range(len(y)):
      if 2.e-3 < y[j]:
        ny1 = j
        break
    for j in range(len(y)):
      if 4.e-3 < y[j]:
        ny2 = j
        break
    filename = "div" + str(itr+1)
    print_scalar(x[1:-1], y[ny1+1:-1], z[::2], div[::2,ny1:,:], div_dir, "div", filename)
    filename = "separation" + str(itr+1)
    print_scalar(x, y[:ny2], Z[::2], u[::2,:ny2,:], sepa_dir, "u", filename)
    itr += 1
 

def calc_POD(dir, name, endT, num_modes):
  files = [f for f in os.listdir(dir) if f.endswith(".vtr")]
  files.sort(key=extract_number)
  Nx, Ny, Nz, x, y, z = getGrid(os.path.join(dir, files[0]))
  Nt = len(files)
  t  = np.linspace(0.e0, endT, Nt)
  D  = np.zeros((Nx*Ny*Nz, Nt), dtype=np.float32)
  itr = 0
  for file in tqdm(files):
    file_path = os.path.join(dir, file)
    p = getScalar(file_path, Nx, Ny, Nz, name)
    D[:,itr] = p.flatten()
    itr += 1
  # D[space, time]
  mean = np.mean(D, axis=-1)
  D = D - mean[:,None]
  # POD main 
  eigenvalues, eigenvectors, modes = snapshot_pod(D)
  time_coef = calc_time_coef(D, modes)
  POD_dir  = os.path.join(dir, "POD")
  os.makedirs(POD_dir, exist_ok=True)
  plot_pod_results(POD_dir, x, y, z, t, eigenvalues, modes, time_coef, num_modes=num_modes)
  plot_reconstruction(POD_dir, x, y, z, modes, time_coef, num_modes, mean=mean)


def time_series(div_dir, sepa_dir, endT, lx1, lx2, ly, lz1, lz2):
  files1 = [f for f in os.listdir(div_dir)  if f.endswith(".vtr")]
  files2 = [f for f in os.listdir(sepa_dir) if f.endswith(".vtr")]
  files1.sort(key=extract_number)
  files2.sort(key=extract_number)
  Nx1, Ny1, Nz1, x1, y1, z1 = getGrid(os.path.join(div_dir,  files1[0]))
  Nx2, Ny2, Nz2, x2, y2, z2 = getGrid(os.path.join(sepa_dir, files2[0]))
  Nt = len(files1)
  t  = np.linspace(0.e0, endT, Nt)
  nx1 = int(Nx1 * lx1 / x1[-1])
  nx2 = int(Nx1 * lx2 / x1[-1])
  nz1 = int(Nz1 * lz1 / z1[-1])
  nz2 = int(Nz1 * lz2 / z1[-1])
  for j in range(len(y1)):
    if ly < y1[j]:
      ny = j
      break
  itr = 0
  x_div  = np.zeros(Nt)
  x_sepa = np.zeros(Nt)
  for file in tqdm(files1):
    file_path1 = os.path.join(div_dir,  file)
    file_path2 = os.path.join(sepa_dir, file)
    Div = getScalar(file_path1, Nx1, Ny1, Nz1, "div")
    U   = getScalar(file_path2, Nx2, Ny2, Nz2, "u")
    div = np.mean(Div[nz1:nz2,ny,nx1:nx2], axis=0)
    u   = np.mean(  U[nz1:nz2,1,:],        axis=0)
    #x_div[itr]  = x1[]
    #x_sepa[itr] = x2[]
  np.save(os.path.join(div_dir,  "div_time_series"),  x_div)
  np.save(os.path.join(sepa_dir, "sepa_time_series"), x_sepa)
  plt.plot(t, x_div)
  plt.plot(t, x_sepa)
  plt.show()


def calc_causality_POD(div_dir, sepa_dir, endT, num_modes):
  path_div  = os.path.join(div_dir,  "POD", "POD_time_coef.npy")
  path_sepa = os.path.join(sepa_dir, "POD", "POD_time_coef.npy")
  time_coef_div  = np.load(path_div)
  time_coef_sepa = np.load(path_sepa)
  data = np.concatenate((time_coef_div[:,:num_modes], time_coef_sepa[:,:num_modes]), axis=-1)
  for i in range(data.shape[0]):
    plt.plot(data[:,i])
    plt.show()
  # causal map
  n   = data.shape[1]
  Map = np.zeros((n,n))
  for j in range(n):
    for i in range(n):
      if i == j:
        Map[j,i] = np.inf
      else:
        #Map[j,i] = embedding_entropy(data[:,i], data[:,j], p=2*num_modes, k=5, Thei=1)
        Map[j,i] = transfer_entropy(data[:,i], data[:,j], k=5, Thei=1)
  plt.imshow(Map, cmap='jet', extent=None, origin='lower')
  plt.xlabel("Effect", fontsize=24)
  plt.ylabel("Cause", fontsize=24)
  plt.colorbar()
  plt.show()


div_dir  = os.path.join(dir, "div")
os.makedirs(div_dir, exist_ok=True)
sepa_dir = os.path.join(dir, "separation")
os.makedirs(sepa_dir, exist_ok=True)


num_modes = 9
make_div_separation(Q_dir, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2)
calc_POD(div_dir, "div", endT, num_modes)
calc_POD(sepa_dir, "u", endT, num_modes)
'''
lx1 =
lx2 =
ly  =
lz1 =
lz2 =
time_series(div_dir, sepa_dir, endT, lx1, lx2, ly, lz1, lz2)
'''
calc_causality_POD(div_dir, sepa_dir, endT, num_modes)

