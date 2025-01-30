import numpy as np
from scipy.spatial.distance import euclidean
from sklearn.feature_selection import mutual_info_regression
from pyinform.mutualinfo import mutual_info
import matplotlib.pyplot as plt


def Takens_embedding(x, tau, dim):
  xe = np.zeros((len(x)-tau*(dim-1),dim), dtype=np.float64)
  xe[:,0] = x[:-tau*(dim-1)]
  for i in range(1,dim):
    xe[:,i] = np.roll(x, -tau*i)[:-tau*(dim-1)]
  return xe


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


def kNN(k, X, Y=None):
  if X.shape[0] < k:
    raise Exception("time series data is too short")
  nt  = X.shape[0] - k + 1
  dim = X.shape[1]
  XNN = np.zeros((nt, k, dim), dtype=np.float64)
  for i in range(k, nt + k):
    d = [euclidean(X[i-1,:], X[j,:]) for j in range(i)]
    indices = np.argsort(d)[:k]
    XNN[i-k,:,:] = X[indices,:]
  if Y is not None:
    return XNN, X[k-1:,:], Y[k-1:,:]
  else:
    return XNN, X[k-1:,:]


def local_constant_pred(t, x, m, dim, k):
  n = len(x)
  tau, _ = search_tau(x, 5, n)
  X = Takens_embedding(x, tau, dim)
  N = len(X[:,0])

  i = n // 2
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

