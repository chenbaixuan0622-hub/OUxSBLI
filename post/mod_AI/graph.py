import numpy as np
import networkx
from sklearn.preprocessing import KBinsDiscretizer
import matplotlib.pyplot as plt


def plot_graph_from_array(A, x, y, num_edges):
  est = KBinsDiscretizer(n_bins=100, encode='ordinal', strategy='uniform')
  est.fit(A)
  A = est.transform(A)
  A[range(len(A[0,:])), range(len(A[0,:]))] = 0
  G = networkx.from_numpy_array(A, create_using=networkx.DiGraph)
  pos = {i: (x.flatten()[i], y.flatten()[i]) for i in G.nodes}

  a = np.sort(A.flatten())
  threshold = a[-num_edges]

  edges_to_remove = [(u, v) for u, v, data in G.edges(data=True) \
                    if data.get('weight', 0) < threshold]
  G.remove_edges_from(edges_to_remove)

  nodes_to_remove = [node for node in G.nodes if G.degree(node) == 0]
  G.remove_nodes_from(nodes_to_remove)

  networkx.draw_networkx_nodes(G, pos, node_color='blue')
  networkx.draw_networkx_edges(G, pos, edge_color='black', arrowstyle='->')
  plt.show()


Nx = 10
Ny = 10
n = Nx * Ny
A = np.random.randn(n,n)
x = np.linspace(0,1,Nx)
y = np.linspace(0,1,Ny)
x, y = np.meshgrid(x, y)
plot_graph_from_array(A, x, y, 100)

