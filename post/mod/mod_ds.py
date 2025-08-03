import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import squareform, pdist
from sklearn.feature_selection import mutual_info_regression
from nolitsa import dimension
from nolitsa import delay
import matplotlib.pyplot as plt


def search_tau(time_series, plot=False):
  n = len(time_series)
  max_lag = n//2
  mutual_info_values = []
  if plot:
    for lag in range(1, max_lag + 1):
      x  = time_series[:-lag].reshape(-1,1)
      y  = time_series[lag:]
      mi = mutual_info_regression(x, y, n_neighbors=5)
      mutual_info_values.append(mi[0])
    mutual_info_values = np.array(mutual_info_values)
    dmi = -mutual_info_values[:-1] + mutual_info_values[1:]
    for lag in range(len(dmi)):
      if dmi[lag] >= 0.e0:
        tau = lag
        break
    plt.plot(range(1, max_lag+1), mutual_info_values)
    plt.xlabel("Lag")
    plt.ylabel("MI")
    plt.show()
  else:
    mi_old = 1e10
    for lag in range(1, max_lag + 1):
      x  = time_series[:-lag].reshape(-1,1)
      y  = time_series[lag:]
      mi = mutual_info_regression(x, y, n_neighbors=5)
      if mi - mi_old >= 0.e0:
        tau = lag
        break
      mi_old = mi
  return tau


def fnn(x, max_dim=10, tau=None):
  if tau is None:
    tau = search_tau(x)
  x = np.array(x, dtype=np.float64)
  dim = max_dim
  for i in range(max_dim):
    fnn_ratios = dimension.fnn(x, dim=[i+1], tau=tau)
    if fnn_ratios[0] == 0:
      dim = i
      break
  return dim


def embedding(x, dim=None, tau=None, Takens=False):
  if tau is None:
    tau = search_tau(x)
  if dim is None:
    if Takens is True:
      dim = 2 * fnn(x, tau=tau)
    else:
      dim = fnn(x, tau=tau)
  X = np.zeros((len(x)-tau*(dim-1),dim), dtype=np.float64)
  X[:,0] = x[:-tau*(dim-1)]
  for i in range(1,dim):
    X[:,i] = np.roll(x, -tau*i)[:-tau*(dim-1)]
  return X


def local_constant_pred(t, x, m, dim, k):
  n = len(x)
  tau, _ = search_tau(x, 5, n)
  X = Takens_embedding(x, tau, dim)
  N = len(X[:,0])
  distance = [euclidean(X[i,:], X[j,:]) if i != j else np.inf for j in range(i)]
  indices  = np.argsort(distance)[:k]
  print("t[i] ", t[i], " indices ", t[indices])
  xp = np.zeros(m, dtype=np.float32)
  xa = np.zeros(m, dtype=np.float32)
  #xp = x[indices[0]:indices[0]+m]
  #xa = x[i:i+m]
  for j in range(m):
    # prediction
    xsum = 0.e0
    for l in range(k):
      xsum += x[indices[l]+j]
    xp[j] = xsum / k
    # answer
    xa[j] = x[i+j]
  return xp, xa


def recurrence_plot(x, tau=None, dim=None):
  if tau is not None and dim is not None:
    X = embedding(x, dim=dim, tau=tau)
    n = X.shape[0]
  else:
    n = len(x)
  dist_matrix = squareform(pdist(X))
  eps = np.median(dist_matrix)
  adjacency = (dist_matrix < eps).astype(np.uint8)
  np.fill_diagonal(adjacency, 0)
  src, dst   = np.nonzero(adjacency)
  edge_index = np.stack([src, dst], axis=0)
  return dist_matrix, adjacency, edge_index

