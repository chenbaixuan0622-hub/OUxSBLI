import numpy as np
from networkx import from_numpy_array, draw
import os
import matplotlib.pyplot as plt
from mod.mod_POD import make_data, make_grid


Q_dir = "../3D_solver/TBL/data"
Lx1   = 28.e-3
Lx2   = 40.e-3
Ly2   = 8.e-3
endT  = 0.1e-3

# D[space=nx*ny, time]
save_path = os.path.join(Q_dir, "D.npy")
if os.path.isfile(save_path):
  D = np.load(save_path)
  _, _, _, _, _, x, y, t = make_grid(Lx1, Lx2, Ly2, endT, Q_dir)
else:
  x, y, t, D = make_data(Lx1, Lx2, Ly2, endT, Q_dir)
  np.save(save_path, D)

nx = len(x)
ny = len(y)
Nt = len(t)
print("nx ", nx, " ny ", ny, " Nt ", Nt)

# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12


# causal matrix
save_path = os.path.join(Q_dir, "GWN_causal_matrix.npy")
A = np.load(save_path)
# A[nx*ny, nx*ny]

X, Y = np.meshgrid(x, y)
X    = X.flatten()
Y    = Y.flatten()
A    = np.abs(A)
for j in range(nx*ny):
  for i in range(nx*ny):
    d = np.sqrt((X[i] - X[j])**2 + (Y[i] - Y[j])**2)
    if d == 0.e0:
      A[j,i] = 0.e0
    if A[j,i] > A[i,j]:
      A[j,i] = A[j,i] - A[i,j]
      A[i,j] = 0.e0
    elif A[j,i] < A[i,j]:
      A[i,j] = A[i,j] - A[j,i]
      A[j,i] = 0.e0

threshold = np.sort(A)[-50]
A = np.where(A <= threshold, 0, A)
Amax = np.max(A)
A = np.array(5.e0 * A / Amax, dtype=np.int32)

G   = from_numpy_array(A)
pos = {i: (X[i], Y[i]) for i in range(nx*ny)}

print("num of nodes ", G.number_of_nodes())
print("num of edges ", G.number_of_edges())

edges_to_remove = [(u, v) for u, v, data in G.edges(data=True) if data.get('weight', 0) == 0]
G.remove_edges_from(edges_to_remove)

nodes_to_remove = [node for node in G.nodes if G.degree(node) == 0]
G.remove_nodes_from(nodes_to_remove)

print("num of nodes ", G.number_of_nodes())
print("num of edges ", G.number_of_edges())

edges = G.edges(data=True)
edge_widths = [d['weight'] for _, _, d in edges]

'''
plt.figure(figsize=(6, 6))
draw(G, pos, with_labels=False, node_size=10, edge_color='blue', \
     width=edge_widths)
plt.show()
'''
