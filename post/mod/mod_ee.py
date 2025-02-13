import numpy as np
from scipy.spatial import KDTree
from scipy.special import psi
from mod.mod_shannon_nn import mutual_info


def EE(x, y, p=1, k=2, Thei=None):
  '''
  x: ndarray, shape (dim_x, T) time series
  y: ndarray, shape (dim_y, T) time series
  p: int, order of the modelt to estimate causality
  k: int, k-th nearest neighbor number
  Thei: int, half-length of Theiler correction window
  out: float, embedding entropy x->y
  '''

  if Thei is None or Thei < p:
    Thei = p

  dy, T = y.shape
  dx    = x.shape[0]
  
  # X [T-p, dx*p], Y[T-p, dy*(p+1)]
  X = np.hstack([x[:,(p - i):(T - i)].T for i in range(1, p+1)])
  Y = np.hstack([y[:,(p - i):(T - i)].T for i in range(p+1)])
 
  N = T - p

  YpNN = np.zeros((N, Y.shape[1] * (dy * (p + 1) + 1)))
  tree_Y = KDTree(Y)

  for i in range(N):
    idx = np.ones(N, dtype=bool)
    idx[max(0, i - Thei):min(i + Thei + 1, N)] = False

    temp_Y = Y[idx]

    _, pnn_idx = tree_Y.query(Y[i], k=min(dy * (p + 1) + 1, temp_Y.shape[0]))
    YpNN[i] = temp_Y[pnn_idx].flatten()

  out = mutual_info(X, YpNN, k=k, Thei=Thei)
  return out

