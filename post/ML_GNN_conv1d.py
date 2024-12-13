import numpy as np
import os
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from scipy.stats import zscore
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mod_AI.gnn import GATModel, plot_graph, update_edge_index
import matplotlib.pyplot as plt


def conv(nx, Lx):
  x    = np.linspace(0.e0, Lx, nx)
  dx   = -x[0] + x[1]
  c    = 1.e0
  CFL  = 0.1e0
  dt   = CFL * dx / c
  endT = 100.e0 * Lx / c
  nt   = int(endT / dt)
  no   = nt // 100

  U = np.zeros((no,nx), dtype=np.float32)
  u = np.random.randn(nx)

  for i in range(nt):
    u_old   = u.copy()
    # 2nd order
    u[2:-2] = u_old[2:-2] - c * dt / dx * 0.5e0 * \
              (u_old[:-4] - 4.e0 * u_old[1:-3] + 3.e0 * u_old[2:-2])
    # 1st order
    u[1]  = u_old[1]  - c * dt / dx * (-u_old[0]  + u_old[1])
    u[-1] = u_old[-1] - c * dt / dx * (-u_old[-2] + u_old[-1])
    u[1:-1] = u_old[1:-1] - c * dt / dx * (-u_old[:-2] + u_old[1:-1])
    # boundary condition
    u[0]    = u_old[0] + 0.5e0 * np.random.randn(1)
    u[-1]   = u[-2]
    # random disturbance
    u += 0.1e0 * np.random.randn(nx)
    if i % 100 == 0:
      U[i//100,:]  = u
  U = zscore(U, axis=None)
  return x, U


def make_graph_data(nx, Lx, device):
  x, U  = conv(nx, Lx)

  Nt = len(U[:,0])

  stride = 4

  num_nodes    = nx
  num_features = 2
  features     = np.zeros((num_nodes, Nt, num_features), dtype=np.float32)
  edge_index   = np.empty([2,0], dtype=np.int64)

  x_mean = np.mean(x)
  x_std  = np.std(x)
  for itr in range(Nt):
    for i in range(nx):
      features[i,itr,0] = (x[i] - x_mean) / x_std
      features[i,itr,1] = U[itr,i]
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
  return data, features, Nt, x


def trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add, max_edges_per_node):
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
                                                        threshold_add, \
                                                        max_edges_per_node)
    t = (t + 1) % Nt
    print(f'Epoch {epoch+1}, Loss: {loss.item():.4f}')
    loss_list.append(loss)
  return data, node_embeddings, loss_list


def main():
  # parameters for conv eq
  nx = 257
  Lx = 2.e0

  x, U = conv(nx, Lx)

  # plot animation
  ims = []
  fig = plt.figure()

  for i in range(len(U[:,0])):
    line = plt.plot(x, U[i,:], color='black')
    ims.append(line)

  ani = animation.ArtistAnimation(fig, ims, interval=200)
  plt.show()

  # parameters for GNN
  dt                 = 2
  in_channels        = 2
  out_channels       = 1
  hidden_channels    = 128
  threshold_remove   = 0.1 
  threshold_add      = 0.9
  max_edges_per_node = 6

  device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
  model = GATModel(in_channels, hidden_channels, out_channels).to(device)
  data, features, Nt, x = make_graph_data(nx, Lx, device)

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
    data, pred, loss_list = trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add, max_edges_per_node)

  # plot predicted contour
  u = pred.detach().cpu().numpy()
  plt.figure(tight_layout=True)
  plt.subplot(121)
  plt.plot(x, u)

  # plot graph
  plt.subplot(122)
  #plt.figure(figsize=(8,4))
  plot_graph(data)
  plt.show()

main()

