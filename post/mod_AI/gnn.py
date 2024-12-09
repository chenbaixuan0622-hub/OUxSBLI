import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.nn import GATConv
from torch_geometric.utils import to_networkx
import networkx as nx


class GATModel(torch.nn.Module):
  def __init__(self, in_channels, hidden_channels, out_channels):
    super(GATModel, self).__init__()
    self.gat1 = GATConv(in_channels, hidden_channels, heads=4)
    self.gat2 = GATConv(hidden_channels * 4, hidden_channels, heads=4)
    self.gat3 = GATConv(hidden_channels * 4, out_channels, heads=1, concat=False)
  
  def forward(self, x, edge_index, edge_attr):
    x1 = F.elu(self.gat1(x,       edge_index, edge_attr))
    x2 = F.elu(self.gat2(x1,      edge_index, edge_attr)) + x1
    x  = F.elu(self.gat3(x1 + x2, edge_index, edge_attr))
    return x


def update_edge_index(device, edge_index, edge_weights, node_embeddings, threshold_remove, threshold_add):
  # cut edge
  mask = edge_weights > threshold_remove
  filtered_edge_index   = edge_index[:, mask]
  filtered_edge_weights = edge_weights[mask]
  num_removed_edges     = edge_index.size(1) - filtered_edge_index.size(1)

  # calc similarity
  num_nodes = node_embeddings.size(0)
  row, col  = torch.triu_indices(num_nodes, num_nodes, offset=1, device=device)
  similarity_matrix = torch.matmul(node_embeddings, node_embeddings.T)
  similarity_scores = similarity_matrix[row, col]

  new_edge_mask   = similarity_scores > threshold_add
  new_edges       = torch.stack((row[new_edge_mask], col[new_edge_mask]), dim=0)
  new_edge_scores = similarity_scores[new_edge_mask]

  if new_edges.size(1) > num_removed_edges:
    top_indices     = torch.topk(new_edge_scores, num_removed_edges).indices
    new_edges       = new_edges[:,top_indices]
    new_edge_scores = new_edge_scores[top_indices]

  non_self_loop_mask = new_edges[0] != new_edges[1]
  new_edges = new_edges[:, non_self_loop_mask]
  new_edge_scores = new_edge_scores[non_self_loop_mask]

  updated_edge_index   = torch.cat([filtered_edge_index, new_edges], dim=1)
  updated_edge_weights = torch.cat([filtered_edge_weights, new_edge_scores])
  return updated_edge_index, updated_edge_weights


def trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add):
  t = 0
  loss_list = []
  for epoch in range(Nt - dt):
    model.train()
    optimizer.zero_grad()
    data.x = features[:,t,:].to(device)

    node_embeddings = model(data.x, data.edge_index, data.edge_attr)
    data.x = features[:,t + dt,:].to(device)
  
    # [xi, xj, rho, u, v, w, p]
    loss = F.mse_loss(node_embeddings, data.x[:,2:7])
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
  return data, loss_list


def plot_graph(data):
  data_cpu = data.clone().cpu()
  G = to_networkx(data_cpu)
  x = data_cpu.x[:,0].numpy()
  y = data_cpu.x[:,1].numpy()
  w = data_cpu.edge_attr.detach().cpu().numpy()

  num_edges = len(w)
  node_positions = {i: (x[i], y[i]) for i in G.nodes}
  edge_widths    = {i: (2.e0 * np.arctan(w[i])) for i in range(num_edges)}
  arrow_size     = 5.e0#[(5.e0 * np.arctan(w[i])) for i in range(num_edges)]
  nx.draw_networkx_nodes(G, pos=node_positions, node_size=10, node_color='blue', alpha=1.0)
  nx.draw_networkx_edges(G, pos=node_positions, edge_color='black', \
                         width=edge_widths, arrowstyle='->', arrowsize=arrow_size)
  nx.draw_networkx_labels(G, pos=node_positions, font_size=0, alpha=0.0)

