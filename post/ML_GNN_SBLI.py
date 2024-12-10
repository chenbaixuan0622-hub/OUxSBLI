import numpy as np
import os
import torch
from torch_geometric.data import Data
from scipy.stats import zscore
from tqdm import tqdm
import cv2
import matplotlib.pyplot as plt
from mod_AI.gnn import GATModel, trainGNN, plot_graph
from mod.mod_read import extract_number, getGrid, getScalar, getVector


def make_graph_data_xy(Lx1, Lx2, Ly2, device, Q_dir):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Nt = len(Q_files)

  stridex = 4
  stridey = 8

  first_path = os.path.join(Q_dir, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  for j in range(Ny):
    if Ly2 < Y[j]:
      ny2 = j
      break
  indicesx = np.arange(nx1, nx2, stridex)
  indicesy = np.arange(0,   ny2, stridey)
  nx  = len(indicesx)
  ny  = len(indicesy)
  x0  = X[indicesx]
  y0  = Y[indicesy]

  indicesx, indicesy = np.meshgrid(indicesx, indicesy)

  print("nx = ", nx, " ny = ", ny)
  # [xi, xj, rho, u, v, w, p]
  num_nodes    = nx * ny
  num_features = 7
  features     = np.zeros((num_nodes, Nt, num_features), dtype=np.float32)
  edge_index   = np.empty([2,0], dtype=np.int64)

  x = zscore(x0)
  y = zscore(y0)

  mean = np.zeros(5)
  std  = np.zeros(5)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    Rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    P       = getScalar(file_path, Nx, Ny, Nz, 'p')
    rho = Rho[0,indicesy,indicesx]
    u   =   U[0,indicesy,indicesx]
    v   =   V[0,indicesy,indicesx]
    w   =   W[0,indicesy,indicesx]
    p   =   P[0,indicesy,indicesx]

    if itr == 0:
      mean[0] = np.mean(rho)
      mean[1] = np.mean(u)
      mean[2] = np.mean(v)
      mean[3] = np.mean(w)
      mean[4] = np.mean(p)
      std[0]  = np.std(rho)
      std[1]  = np.std(u)
      std[2]  = np.std(v)
      std[3]  = np.std(w)
      std[4]  = np.std(p)
      u0 = u
      v0 = v
      X, Y = np.meshgrid(x0*1e3, y0*1e3)
      plt.contourf(X, Y, u0, levels=50, cmap='jet')
      save_path = os.path.join(Q_dir, "init.png")
      plt.savefig(save_path)
      plt.close()
    rho = zscore(rho, axis=None)
    u   = zscore(u,   axis=None)
    v   = zscore(v,   axis=None)
    w   = zscore(w,   axis=None)
    p   = zscore(p,   axis=None)
    for j in range(ny):
      for i in range(nx):
        idx = j * nx + i
        features[idx,itr,0] = x[i]
        features[idx,itr,1] = y[j]
        features[idx,itr,2] = rho[j,i]
        features[idx,itr,3] = u[j,i]
        features[idx,itr,4] = v[j,i]
        features[idx,itr,5] = w[j,i]
        features[idx,itr,6] = p[j,i]
    if itr == 0:
      for j in range(ny):
        for i in range(nx):
          idx = j * nx + i
          if i < nx-1:
            if 0.e0 < u0[j,i] + u0[j,i+1]:
              index = np.array([idx,idx+1]).reshape(2,1)
            else:
              index = np.array([idx+1,idx]).reshape(2,1)
            edge_index = np.append(edge_index, index, axis=-1)
          if j < ny-1:
            if 0.e0 < v0[j,i] + v0[j+1,i]:
              index = np.array([idx,idx+nx]).reshape(2,1)
            else:
              index = np.array([idx+nx,idx]).reshape(2,1)
            edge_index = np.append(edge_index, index, axis=-1)
    itr += 1

  num_edge     = len(edge_index[0,:])
  features     = torch.tensor(features,   device='cpu')
  edge_index   = torch.tensor(edge_index, device=device)
  edge_weights = torch.randn((num_edge,), requires_grad=True, device=device)
  features0    = features[:,0,:].to(device)
  data = Data(x=features0, edge_index=edge_index, edge_attr=edge_weights)
  return data, features, Nt, X, Y, mean, std


def main():
  # parameter
  Q_dir = "../../z6mm"
  Lx1   = 28.e-3
  Lx2   = 40.e-3
  Ly2   = 8.e-3
  dt    = 1

  # parameters for GNN
  in_channels       = 7
  out_channels      = 5
  hidden_channels   = 128
  threshold_remove  = 0.2 
  threshold_add     = 0.8
  max_edge_per_node = 6

  device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
  model = GATModel(in_channels, hidden_channels, out_channels).to(device)
  data, features, Nt, X, Y, mean, std = make_graph_data_xy(Lx1, Lx2, Ly2, device, Q_dir)

  plt.figure(figsize=(8,4))
  plot_graph(data)
  save_path = os.path.join(Q_dir, "graph_before.png")
  plt.savefig(save_path)
  plt.close()
  
  data = data.to(device)
  optimizer = torch.optim.Adam(list(model.parameters()) + [data.edge_attr], \
                               lr=0.005, weight_decay=1e-4)

  data, pred, loss_list = trainGNN(device, model, optimizer, Nt, dt, data, features, \
                                   threshold_remove, threshold_add, max_edge_per_node)

  # plot loss
  '''
  plt.plot(loss_list.cpu())
  plt.xscale('log')
  plt.yscale('log')
  plt.show()
  plt.close()
  '''

  # plot predicted contour
  pred_cpu = pred.detach().cpu()
  u = pred_cpu[:,1].numpy()
  u = u * std[1] + mean[1]
  u = np.reshape(u, X.shape)
  plt.contourf(X, Y, u, levels=50, cmap='jet')
  save_path = os.path.join(Q_dir, "pred.png")
  plt.savefig(save_path)
  plt.close()

  # plot graph
  plt.figure(figsize=(8,4))
  plot_graph(data)
  save_path = os.path.join(Q_dir, "graph_after.png")
  plt.show()
  plt.savefig(save_path)
  plt.close()

main()

