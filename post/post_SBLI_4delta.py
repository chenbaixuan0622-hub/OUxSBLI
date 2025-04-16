import numpy as np
import os
from scipy.signal import welch
import pandas as pd
from statsmodels.tsa.stattools import grangercausalitytests
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getScalar, getVector, extract_number
from mod.mod_POD import snapshot_pod, calc_time_coef, plot_pod_results, reconstruct_POD
from mod.mod_ds import search_tau, fnn
from mod.mod_shannon import embedding_entropy_surrogate as embedding_entropy
from mod.mod_shannon import transfer_entropy_surrogate as transfer_entropy
from mod.mod_shannon import mutual_info
from mod.mod_plot import print_scalar


#set_Params()

# parameter
Q_dir = "../../SBLI/SBLI_4delta/2.2ms"
Lx1   = 25.e-3
Lx2   = 37.e-3
Ly1   = 0.e-3
Ly2   = 8.e-3
Lz1   = 0.e-3 
Lz2   = 8.e-3
endT  = 1.e-3
# wall unit
rhow  = 0.1886e0
ut    = 23.33e0
mu    = 1.74e-5
nu    = mu / rhow


def set_range(dir, files, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, ut=None, nu=None):
  Nx, Ny, Nz, X, Y, Z = getGrid(os.path.join(dir, files[0]))
  def search(X, N, L1, L2, ut=None, nu=None):
    n1 = 0
    n2 = N
    if ut is None or nu is None:
      for i in range(N):
        if X[i] >= L1:
          n1 = i
          break
      for i in range(N):
        if X[i] >= L2:
          n2 = i
          break
    else:
      for i in range(N):
        if ut * X[i] / nu >= L1:
          n1 = i
          break
      for i in range(N):
        if ut * X[i] / nu >= L2:
          n2 = i
          break
    return n1, n2
  nx1, nx2 = search(X, Nx, Lx1, Lx2)
  ny1, ny2 = search(Y, Ny, Ly1, Ly2, ut=ut, nu=nu)
  nz1, nz2 = search(Z, Nz, Lz1, Lz2)
  x = X[nx1:nx2]
  y = Y[ny1:ny2]
  z = Z[nz1:nz2]
  print("x from", nx1, " to ", nx2, " y from ", ny1, " to ", ny2, " z from ", nz1, " to ", nz2)
  return Nx, Ny, Nz, x, y, z, nx1, nx2, ny1, ny2, nz1, nz2


def make_data(Q_dir, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2):
  files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  files.sort(key=extract_number)
  Nx, Ny, Nz, X, Y, Z = getGrid(os.path.join(Q_dir, files[0]))
  for i in range(Nx):
    if Lx1 < X[i]:
      nx1 = i
      break
  for i in range(Nx):
    if Lx2 < X[i]:
      nx2 = i
      break
  for j in range(Ny):
    if Ly2 < Y[j]:
      ny = j
      break
  y = Y[:ny]
  z = Z[::2]
  itr = 0
  for file in tqdm(files):
    file_path = os.path.join(Q_dir, file)
    if file_path == os.path.join(Q_dir, "TKE.vtr") \
    or file_path == os.path.join(Q_dir, "ReynoldsStress.vtr") \
    or file_path == os.path.join(Q_dir, "Qmean.vtr"):
      continue
    U, _, _   = getVector(file_path, Nx, Ny, Nz, 'velocity')
    P         = getScalar(file_path, Nx, Ny, Nz, 'p')
    u = U[::2,:ny,nx1:nx2]
    p = P[::2,:ny,nx1:nx2]
    x = np.linspace(X[nx1], X[nx2], u.shape[-1])
    filename = "p" + str(itr+1).zfill(5)
    print_scalar(x, y, z, p, p_dir, "p", filename)
    filename = "u" + str(itr+1).zfill(5)
    print_scalar(x, y, z, u, u_dir, "u", filename)
    itr += 1
 

def calc_POD(dir, name, endT, num_modes, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, ut=None, nu=None):
  files = [f for f in os.listdir(dir) if f.endswith(".vtr")]
  files.sort(key=extract_number)
  Nx, Ny, Nz, x, y, z, nx1, nx2, ny1, ny2, nz1, nz2 = set_range(dir, files, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, ut=ut, nu=nu)
  Nt = len(files)
  t  = np.linspace(0.e0, endT, Nt)
  D  = np.zeros(((nx2-nx1)*(ny2-ny1)*(nz2-nz1), Nt), dtype=np.float32)
  for itr, file in tqdm(enumerate(files)):
    file_path = os.path.join(dir, file)
    p = getScalar(file_path, Nx, Ny, Nz, name)
    D[:,itr] = p[nz1:nz2,ny1:ny2,nx1:nx2].flatten()
  # separation
  #D = np.minimum(D, 0.e0)
  mean = np.mean(D, axis=-1)
  D -= mean[:,None]
  # POD main 
  Sigma, U, phi = snapshot_pod(D)
  time_coef = calc_time_coef(D, phi)
  del D
  POD_dir  = os.path.join(dir, "POD")
  os.makedirs(POD_dir, exist_ok=True)
  reconstruct_POD(POD_dir, x, y, z, phi, time_coef, num_modes, mean=mean)
  plot_pod_results(POD_dir, x, y, z, t, Sigma, phi, time_coef, num_modes=num_modes)


def calc_causality_POD(p_dir, u_dir, endT, num_modes):
  path_p_down = os.path.join(p_dir, "POD_down",       "POD_time_coef.npy")
  path_p_up   = os.path.join(p_dir, "POD_up",         "POD_time_coef.npy")
  path_sepa   = os.path.join(u_dir, "POD_separation", "POD_time_coef.npy")
  path_u_visc = os.path.join(u_dir, "POD_visc",       "POD_time_coef.npy")
  path_u_buf  = os.path.join(u_dir, "POD_buf",        "POD_time_coef.npy")
  path_u_log  = os.path.join(u_dir, "POD_log",        "POD_time_coef.npy")
  time_coef1  = np.load(path_p_down)
  time_coef2  = np.load(path_p_up)
  time_coef3  = np.load(path_sepa)
  time_coef4  = np.load(path_u_visc)
  time_coef5  = np.load(path_u_buf)
  time_coef6  = np.load(path_u_log)
  data = np.concatenate((time_coef1[:,:num_modes], time_coef2[:,:num_modes], \
                         time_coef3[:,:num_modes], time_coef4[:,:num_modes], \
                         time_coef5[:,:num_modes], time_coef6[:,:num_modes]), axis=-1)
  # causal map
  n = data.shape[1]
  Map_te = np.zeros((n,n))
  Map_gc = np.zeros((2,n,n)) # granger causality
  #Map_mi     = np.zeros((n,n))
  #Map_corr   = np.zeros((n,n))
  dim = np.ones(n, dtype=np.int32)
  tau = np.ones(n, dtype=np.int32)
  for i in range(n):
    tau[i] = int(search_tau(data[:,i]))
    dim[i] = int(fnn(data[:,i], tau=tau[i]))
  save_path = os.path.join(p_dir, "tau.npy")
  np.save(save_path, tau)
  save_path = os.path.join(p_dir, "dim.npy")
  np.save(save_path, dim)
  for j in tqdm(range(n)):
    for i in range(n):
      if i == j:
        #Map_te[j,i] = np.inf
        Map_gc[:,j,i] = np.inf
      else:
        #Map_te[j,i] = embedding_entropy(data[:,i], data[:,j], p=2*num_modes, k=5, Thei=1)
        #Map_te[j,i] = transfer_entropy(data[:,i], data[:,j], \
        #                               p=min(dim[i], dim[j]), tau=min(tau[i], tau[j]), k=5, Thei=1)
        df = pd.DataFrame({'x': data[:,i], 'y': data[:,j]}) # causality y -> x
        lag = min(tau[i], tau[j])
        gc = grangercausalitytests(df, [lag], verbose=False)
        Map_gc[0,j,i] = gc[lag][0]['ssr_ftest'][0] # f values
        Map_gc[1,j,i] = gc[lag][0]['ssr_ftest'][1] # p values
      #Map_mi[j,i]   = mutual_info(data[:,i].reshape(-1,1), data[:,j].reshape(-1,1), k=5, Thei=10)
      #Map_corr[j,i] = np.corrcoef(data[:,i], data[:,j])[0,1]
  #save_path = os.path.join(p_dir, "Map_te.npy")
  #np.save(save_path, Map_te)
  save_path = os.path.join(p_dir, "Map_gc.npy")
  np.save(save_path, Map_gc)


p_dir = os.path.join(Q_dir, "p")
os.makedirs(p_dir, exist_ok=True)
u_dir = os.path.join(Q_dir, "u")
os.makedirs(u_dir, exist_ok=True)


num_modes = 100
#make_data(Q_dir, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2)

# separation
#calc_POD(u_dir, "u",  endT, num_modes, Lx1=28.e-3, Lx2=37.e-3, Ly1=0.e0, Ly2=2e-3, Lz1=0.e0, Lz2=8.e-3)

# structures coming from upstream boundary layer (viscous sub-layer)
#calc_POD(u_dir, "u",  endT, num_modes, Lx1=20.e-3, Lx2=28.e-3, Ly1=0.e0, Ly2=5.e0, Lz1=0.e0, Lz2=8.e-3, ut=ut, nu=nu)
# structures coming from upstream boundary layer (buffer layer)
#calc_POD(u_dir, "u",  endT, num_modes, Lx1=20.e-3, Lx2=28.e-3, Ly1=5.e0, Ly2=30.e0, Lz1=0.e0, Lz2=8.e-3, ut=ut, nu=nu)
# structures coming from upstream boundary layer (log layer)
#calc_POD(u_dir, "u",  endT, num_modes, Lx1=20.e-3, Lx2=28.e-3, Ly1=30.e0, Ly2=1000.e0, Lz1=0.e0, Lz2=8.e-3, ut=ut, nu=nu)

# upstream pressure field
#calc_POD(p_dir, "p",  endT, num_modes, Lx1=28.e-3, Lx2=32.e-3, Ly1=1.e-3, Ly2=8.e-3, Lz1=0.e0, Lz2=8.e-3)
# compression wave
#calc_POD(p_dir, "p",  endT, num_modes, Lx1=32.e-3, Lx2=37.e-3, Ly1=1.e-3, Ly2=8.e-3, Lz1=0.e0, Lz2=8.e-3)

calc_causality_POD(p_dir, u_dir, endT, num_modes)

