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
    self.gat2 = GATConv(hidden_channels * 4, out_channels, heads=1, concat=False)
  
  def forward(self, x, edge_index, edge_attr):
    x1 = F.elu(self.gat1(x,  edge_index, edge_attr))
    x  = F.elu(self.gat2(x1, edge_index, edge_attr))
    return x


def edge_sort_key(edge):
  return edge[1]


@torch.jit.script
def update_edge_index(
  device:             torch.device,
  edge_index:         torch.Tensor,
  edge_weights:       torch.Tensor,
  node_embeddings:    torch.Tensor,
  threshold_remove:   float,
  threshold_add:      float,
  max_edges_per_node: int
  ):

  # cut edges with weights below threshold
  mask = edge_weights > threshold_remove
  filtered_edge_index   = edge_index[:, mask]
  filtered_edge_weights = edge_weights[mask]
  num_removed_edges     = edge_index.size(1) - filtered_edge_index.size(1)

  # calculate similarity
  num_nodes    = node_embeddings.size(0)
  triu_indices = torch.triu_indices(num_nodes, num_nodes, offset=1, device=device)
  row = triu_indices[0]
  col = triu_indices[1]
  similarity_matrix = torch.matmul(node_embeddings, node_embeddings.T)
  similarity_scores = similarity_matrix[row, col]

  # add new edges based on similarity
  new_edge_mask   = similarity_scores > threshold_add
  new_edges       = torch.stack((row[new_edge_mask], col[new_edge_mask]), dim=0)
  new_edge_scores = similarity_scores[new_edge_mask]

  if new_edges.size(1) > num_removed_edges:
    top_indices     = torch.topk(new_edge_scores, num_removed_edges).indices
    new_edges       = new_edges[:,top_indices]
    new_edge_scores = new_edge_scores[top_indices]

  # prevent self-loops
  non_self_loop_mask = new_edges[0] != new_edges[1]
  new_edges = new_edges[:, non_self_loop_mask]
  new_edge_scores = new_edge_scores[non_self_loop_mask]

  # limit the number of edges per node
  final_edges   = torch.empty(0, 2, dtype=torch.long, device=device)
  final_weights = torch.empty(0,   dtype=torch.float, device=device)

  src_edges_dict   = {i: torch.empty(0, dtype=torch.long,  device=device) for i in range(num_nodes)}
  src_weights_dict = {i: torch.empty(0, dtype=torch.float, device=device) for i in range(num_nodes)}

  for i in range(new_edges.size(1)):
    src, tgt = new_edges[0, i].item(), new_edges[1, i].item()
    weight   = new_edge_scores[i].item()

    # collect edges and weights in separate dictionaries
    src_edges_dict[src] = torch.cat([src_edges_dict[src], torch.tensor([tgt], dtype=torch.long, device=device)])
    src_weights_dict[src] = torch.cat([src_weights_dict[src], torch.tensor([weight], dtype=torch.float, device=device)])

  for src in range(num_nodes):
    if src_edges_dict[src].size(0) > 0:
      sorted_indices = torch.argsort(
        src_weights_dict[src], descending=True,
      )[:max_edges_per_node]
      for idx in sorted_indices:
        tgt = src_edges_dict[src][idx]
        weight = src_weights_dict[src][idx]
        final_edges = torch.cat([final_edges, torch.tensor([[src, tgt]], dtype=torch.long, device=device)], dim=0)
        final_weights = torch.cat([final_weights, torch.tensor([weight], dtype=torch.float, device=device)], dim=0)

  # construct the final edge index and weights
  if final_edges.size(0) > 0:
    updated_edge_index = torch.cat([filtered_edge_index, final_edges.T], dim=1)
    updated_edge_weights = torch.cat([filtered_edge_weights, final_weights])
  else:
    updated_edge_index = filtered_edge_index
    updated_edge_weights = filtered_edge_weights

  return updated_edge_index, updated_edge_weights


def trainGNN(device, model, optimizer, Nt, dt, data, features, threshold_remove, threshold_add, max_edges_per_node):
  loss_list = []
  for epoch in range(10):
    for t in range(Nt - dt):
      model.train()
      optimizer.zero_grad()
      data.x = features[:,t,:].to(device)

      node_embeddings = model(data.x, data.edge_index, data.edge_attr)
      data.x = features[:,t + dt,:].to(device)
  
      # [xi, xj, rho, u, v, w, p]
      loss = F.mse_loss(node_embeddings, data.x[:,2:7])
      loss.backward()
      optimizer.step()
      data.edge_index, data.edge_attr = update_edge_index(
        device,
        data.edge_index,
        data.edge_attr,
        node_embeddings,
        threshold_remove,
        threshold_add,
        max_edges_per_node
      )

      print(f'Epoch {epoch+1}, t {t}, Loss: {loss.item():.4f}')
      loss_list.append(loss)
  return data, node_embeddings, loss_list


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

