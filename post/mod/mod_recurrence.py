import numpy as np
from numba import njit
from sklearn.metrics.pairwise import pairwise_distances
from scipy.stats import norm


#@njit(cache=True, fastmath=True, nogil=True)
def recurrence_plot(x):
  n = len(x)
  D = np.zeros((n,n), dtype=np.float32)
  d = np.zeros((n,n), dtype=np.int32)
  # calc distance
  for j in range(n):
    for i in range(n):
      if i != j:
        D[i,j] = -x[i] + x[j]
        D[j,i] = D[i,j]
      else:
        D[i,j] = 0.e0
  eps = np.median(D.flatten())
  # 0 or 1
  for j in range(n):
    for i in range(n):
      if D[i,j] < eps:
        d[i,j] = 1.e0
        d[j,i] = 1.e0
      else:
        d[i,j] = 0.e0
        d[j,i] = 0.e0
  
  return d


def excessive_recurrence_test(Rx, Ry):
  px  = Rx.mean()
  py  = Ry.mean()
  Rxy = Rx * Ry
  u_obs = Rxy.mean()

  u_expected = px * py
  variance   = px * py * (1.e0 - px * py)

  z = (u_obs - u_expected) / np.sqrt(variance / Rx.size)
  p = 1.e0 - norm.cdf(z)
  return Rxy, z, p

