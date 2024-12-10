import numpy as np
import os
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from scipy.stats import zscore
from tqdm import tqdm
import matplotlib.pyplot as plt
from mod_AI.gnn import GATModel, plot_graph, update_edge_index
import matplotlib.pyplot as plt


def conv(nx, Lx):
  x    = np.linspace(0.e0, Lx, nx)
  dx   = -x[0] + x[1]
  c    = 1.e0
  CFL  = 0.1e0
  dt   = CFL * dx / c
  endT = 2.e0 * Lx / c
  nt   = int(endT / dt)

  U   = np.zeros((nt,nx), dtype=np.float32)
  u   = np.zeros(nx,      dtype=np.float32)
  u[int(0.1 / dx):int(0.2 / dx + 1)] = 1.e0

  for i in range(nt):
    u_old   = u.copy()
    u[1:-1] = u_old[1:-1] - c * dt / dx * (-u_old[:-2] + u_old[1:-1])
    u[0]    = u[-2]
    u[-1]   = u[1]
    U[i,:]  = u
  return x, U


def make_graph_data(nx, Lx, device):
  x, U  = conv(nx, Lx)

  Nt = len(U[:,0])

  stride = 4

  num_nodes    = nx
  num_features = 2
  features     = np.zeros((num_nodes, Nt, num_features), dtype=np.float32)
  edge_index   = np.empty([2,0], dtype=np.int64)

  u_mean = np.mean(U[0,:])
  u_std  = np.std(U[0,:])
  x_mean = np.mean(x)
  x_std  = np.std(x)
  for itr in range(Nt):
    for i in range(nx):
      features[i,itr,0] = (x[i] - x_mean) / x_std
      features[i,itr,1] = (U[itr,i] - u_mean) / u_std
      if itr == 0:
        if i < nx-1:
          index = np.array([i,i+1]).reshape(2,1)
          edge_index = np.append(edge_index, index, axis=-1)

  num_edge     = nx-1
  features     = torch.tensor(features,   device='cpu')
  edge_index   = torch.tensor(edge_index, device=device)
  edge_weights = torch.randn((num_edge,), requires_grad=True, device=device)
  features0    = features[:,0,:].to(device)
  data = Data(x=features0, edge_index=edge_index, edge_attr=edge_weights)
  return data, features, Nt, x, u_mean, u_std


def trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add):
  t = 0
  loss_list = []
  for epoch in range(Nt - dt):
    model.train()
    optimizer.zero_grad()
    data.x = features[:,t,:].to(device)

    node_embeddings = model(data.x, data.edge_index, data.edge_attr)
    data.x = features[:,t + dt,:].to(device)
  
    # [xi, u]
    loss = F.mse_loss(node_embeddings, data.x[:,1])
    loss.backward()
    optimizer.step()
    data.edge_index, data.edge_attr = update_edge_index(device, \
                                                        data.edge_index, \
                                                        data.edge_attr, \
                                                        node_embeddings, \
                                                        threshold_remove, \
                                                        threshold_add)
    t = (t + 1) % Nt
    print(f'Epoch {epoch+1}, Loss: {loss.item():.4f}')
    loss_list.append(loss)
  return data, node_embeddings, loss_list


def main():
  # parameters for conv eq
  nx = 33
  Lx = 1.e0

  # parameters for GNN
  dt               = 2
  in_channels      = 2
  out_channels     = 1
  hidden_channels  = 128
  threshold_remove = 0.1 
  threshold_add    = 0.9

  device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
  model = GATModel(in_channels, hidden_channels, out_channels).to(device)
  data, features, Nt, x, mean, std = make_graph_data(nx, Lx, device)

  print(data.x.size())
  print(data.edge_index.size())
  print(data.edge_attr.size())

  plt.figure(figsize=(8,4))
  plot_graph(data)
  plt.show()
  plt.close()
  
  data = data.to(device)
  optimizer = torch.optim.Adam(list(model.parameters()) + [data.edge_attr], \
                               lr=0.005, weight_decay=1e-4)

  for i in range(3):
    data, pred, loss_list = trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add)

  # plot predicted contour
  u = pred.detach().cpu().numpy()
  u = u * std + mean
  plt.plot(x, u)
  plt.show()
  plt.close()

  # plot graph
  plt.figure(figsize=(8,4))
  plot_graph(data)
  plt.show()
  plt.close()

main()

