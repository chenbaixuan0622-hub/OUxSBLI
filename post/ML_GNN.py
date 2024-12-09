import numpy as np
import os
import torch
from torch_geometric.data import Data
from tqdm import tqdm
import cv2
import matplotlib.pyplot as plt
from mod_AI.gnn import GATModel, trainGNN, plot_graph
from mod.mod_read import extract_number, getGrid, getScalar, getVector


def make_graph_data(Lx1, Lx2, device, Q_dir):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Nt = len(Q_files)

  first_path = os.path.join(Q_dir, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  nx  = (nx2 - nx1) // 8
  nz  = Nz // 16
  x   = np.linspace(X[nx1], X[nx2], nx)
  z   = np.linspace(Z[0],   Z[-1],  nz)
  
  print("nx = ", nx)
  print("nz = ", nz)
  # [xi, xj, rho, u, v, w, p]
  num_nodes    = nx * nz
  num_features = 7
  features     = np.zeros((num_nodes, Nt, num_features), dtype=np.float32)
  edge_index   = np.empty([2,0], dtype=np.int64)

  x_mean = np.mean(x)
  z_mean = np.mean(z)
  x_std  = np.std(x)
  z_std  = np.std(z)
  x      = (x - x_mean) / x_std
  z      = (z - z_mean) / z_std

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    Rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    P       = getScalar(file_path, Nx, Ny, Nz, 'p')
    Rho = np.reshape(Rho[:,0,nx1:nx2], [Nz, nx2-nx1])
    U   = np.reshape(  U[:,0,nx1:nx2], [Nz, nx2-nx1])
    V   = np.reshape(  V[:,0,nx1:nx2], [Nz, nx2-nx1])
    W   = np.reshape(  W[:,0,nx1:nx2], [Nz, nx2-nx1])
    P   = np.reshape(  P[:,0,nx1:nx2], [Nz, nx2-nx1])
    rho = cv2.resize(Rho, (nz, nx))
    u   = cv2.resize(  U, (nz, nx))
    v   = cv2.resize(  V, (nz, nx))
    w   = cv2.resize(  W, (nz, nx))
    p   = cv2.resize(  P, (nz, nx))
    if itr == 0:
      u0 = u
      v0 = v
    rho_mean = np.mean(rho)
    u_mean   = np.mean(u)
    v_mean   = np.mean(v)
    w_mean   = np.mean(w)
    p_mean   = np.mean(p)
    rho_std  = np.std(rho)
    u_std    = np.std(u)
    v_std    = np.std(v)
    w_std    = np.std(w)
    p_std    = np.std(p)
    rho = (rho - rho_mean) / rho_std
    u   = (  u -   u_mean) /   u_std
    v   = (  v -   v_mean) /   v_std
    w   = (  w -   w_mean) /   w_std
    p   = (  p -   p_mean) /   p_std
    for k in range(nz):
      for i in range(nx):
        idx = k * nx + i
        features[idx,itr,0] = x[i]
        features[idx,itr,1] = z[k]
        features[idx,itr,2] = rho[k,i]
        features[idx,itr,3] = u[k,i]
        features[idx,itr,4] = v[k,i]
        features[idx,itr,5] = w[k,i]
        features[idx,itr,6] = p[k,i]
    if itr == 0:
      for k in range(nz):
        for i in range(nx):
          idx = k * nx + i
          if i < nx-1:
            if 0.e0 < u0[k,i] + u0[k,i+1]:
              index = np.array([idx,idx+1]).reshape(2,1)
            else:
              index = np.array([idx+1,idx]).reshape(2,1)
            edge_index = np.append(edge_index, index, axis=-1)
          if k < nz-1:
            if 0.e0 < v0[k,i] + v0[k+1,i]:
              index = np.array([idx,idx+nx]).reshape(2,1)
            else:
              index = np.array([idx+nx,idx]).reshape(2,1)
            edge_index = np.append(edge_index, index, axis=-1)
          else:
            if 0.e0 < v0[k,i] + v0[0,i]:
              index = np.array([idx,i]).reshape(2,1)
            else:
              index = np.array([i,idx]).reshape(2,1)
            edge_index = np.append(edge_index, index, axis=-1)
    itr += 1

  num_edge     = len(edge_index[0,:])
  features     = torch.tensor(features,   device=device)
  edge_index   = torch.tensor(edge_index, device=device)
  edge_weights = torch.randn((num_edge,), requires_grad=True, device=device)
  data = Data(x=features[:,0,:], edge_index=edge_index, edge_attr=edge_weights)
  return data, features, Nt


def main():
  # parameter
  Q_dir = "../../yp269"
  Lx1   = 25.e-3
  Lx2   = 35.e-3
  dt    = 1

  # parameters for GNN
  in_channels      = 7
  out_channels     = 5
  hidden_channels  = 128
  threshold_remove = 0.1 
  threshold_add    = 0.9

  device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
  model = GATModel(in_channels, hidden_channels, out_channels).to(device)
  data, features, Nt = make_graph_data(Lx1, Lx2, device, Q_dir)

  plt.figure(figsize=(4,8))
  plot_graph(data)
  plt.show()
  plt.close()
  
  data = data.to(device)
  optimizer = torch.optim.Adam(list(model.parameters()) + [data.edge_attr], \
                               lr=0.005, weight_decay=1e-4)

  for i in range(10):
    data, loss_list = trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add)

  plt.figure(figsize=(4,8))
  plot_graph(data)
  plt.show()
  plt.close()

main()

