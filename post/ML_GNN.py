import numpy as np
import torch
from torch_geometric.data import Data
import matplotlib.pyplot as plt
from mod_AI.gnn import GATModel, trainGNN, plot_graph


# grid info
Nx = 16
Nz = 8
Lx = 2.e0
Lz = 1.e0
x = np.linspace(0.e0, Lx, Nx)
z = np.linspace(0.e0, Lz, Nz)

# parameters for time-series graph data
num_nodes = Nx * Nz
# [xi, xj, rho, u, v, w, p, rho', u', v', w', p']
num_features = 12
num_initial_edges = 2*20
Nt = 100
dt = 10

# parameters for GNN
in_channels      = num_features
out_channels     = 5
hidden_channels  = 128
threshold_remove = 0.1 
threshold_add    = 0.9

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

features = torch.zeros((num_nodes, Nt, num_features), device=device)
for j in range(Nz):
  for i in range(Nx):
    for itr in range(Nt):
      idx = j * Nx + i 
      features[idx,itr,0]    = x[i]
      features[idx,itr,1]    = z[j]
      features[idx,itr,2:12] = torch.randn(10)

labels = torch.randint(0, 3, (num_nodes,), device=device)

edge_index   = torch.randint(0, num_nodes, (2, num_initial_edges), device=device)
edge_weights = torch.randn((num_initial_edges,), \
                           requires_grad=True, device=device)

t = 0
data = Data(x=features[:,t,:], y=labels, \
            edge_index=edge_index, edge_attr=edge_weights)

model = GATModel(in_channels, hidden_channels, out_channels).to(device)
data = data.to(device)
optimizer = torch.optim.Adam(list(model.parameters()) + [edge_weights], \
                             lr=0.005, weight_decay=1e-4)

data, loss_list = trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add)
data1 = data.clone().cpu()

plt.figure(figsize=(4,8))
plot_graph(data1)
plt.show()
plt.close()

