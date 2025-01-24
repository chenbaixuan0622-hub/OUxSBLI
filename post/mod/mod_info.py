import numpy as np
from scipy.spatial import cKDTree
from scipy.special import digamma
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.feature_selection import mutual_info_regression
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy
from pyinform.mutualinfo import mutual_info
from mod.mod_ds import Takens_embedding, search_tau, kNN


def mutual_information(X, Y, k=3, theiler_window=10, eps=1e-10):
  '''
  X: [n_samples, dim_X]
  Y: [n_smaples, dim_Y]
  '''
  N, dx = X.shape
  _, dy = Y.shape

  XY = np.hstack((X, Y))
  tree_XY = cKDTree(XY)
  tree_X  = cKDTree(X)
  tree_Y  = cKDTree(Y)

  distances_XY = []
  distances_X  = []
  distances_Y  = []

  for i in range(N):
    # define the range of points to exclude
    exclude_indices = list(range(max(0, i-theiler_window), min(N, i+theiler_window+1)))
    # query neighbors for joint space
    dist_XY, indices_XY = tree_XY.query(XY[i], k=k+1, p=2)
    valid_distances_XY  = [dist_XY[j] for j in range(len(indices_XY)) if indices_XY[j] not in exclude_indices]
    distances_XY.append(valid_distances_XY[k-1] if len(valid_distances_XY) >= k else np.inf)

    # query neighbors for X space
    dist_X, indices_X = tree_X.query(X[i], k=k+1, p=2)
    valid_distances_X = [dist_X[j] for j in range(len(indices_X)) if indices_X[j] not in exclude_indices]
    distances_X.append(valid_distances_X[k-1] if len(valid_distances_X) >= k else np.inf)

    # query neighbors for Y space
    dist_Y, indices_Y = tree_Y.query(Y[i], k=k+1, p=2)
    valid_distances_Y = [dist_Y[j] for j in range(len(indices_Y)) if indices_Y[j] not in exclude_indices]
    distances_Y.append(valid_distances_Y[k-1] if len(valid_distances_Y) >= k else np.inf)

  distances_XY = np.array(distances_XY)
  distances_X  = np.array(distances_X)
  distances_Y  = np.array(distances_Y)

  distances_XY[distances_XY < eps] = eps
  distances_X[distances_X < eps] = eps
  distances_Y[distances_Y < eps] = eps

  print(np.min(distances_XY), np.max(distances_XY))
  print(np.min(distances_X),  np.max(distances_X))
  print(np.min(distances_Y),  np.max(distances_Y))

  c = digamma(k) - digamma(N)
  H_XY = c + (dx + dy) * np.mean(np.log(2.e0 * distances_XY[distances_XY < np.inf]))
  H_X  = c +       dx  * np.mean(np.log(2.e0 * distances_X[distances_X < np.inf]))
  H_Y  = c +       dy  * np.mean(np.log(2.e0 * distances_Y[distances_Y < np.inf]))
  MI = H_X + H_Y - H_XY
  return MI


def cross_corr(x, y):
  X = np.mean(x)
  Y = np.mean(y)
  corr = np.mean((x - X) * (y - Y)) \
         / (np.sqrt(np.mean((x - X)**2)) * np.sqrt(np.mean((y - Y)**2)))
  return corr


def X_embedding(x, p, tau=None):
  if tau is None:
    tau = search_tau(x)
  X = np.zeros((len(x)-tau*p,p+1), dtype=np.float64)
  for i in range(p+1):
    X[:,i] = np.roll(x, tau*i)[tau*p:]
  return X


def Y_embedding(y, p, tau=None):
  if tau is None:
    tau = search_tau(y)
    Y = X_embedding(y, p-1, tau)
  else:
    Y = X_embedding(y, p-1, tau)
  return np.delete(Y, 0, axis=0)


def EE(x, y, p):
  # x[time], y[time]
  # p: dimension
  X = X_embedding(x, p)
  Y = Y_embedding(y, p)

  # X, Y: embedded time series data
  # X[time, dim], Y[time, dim]

  # length of X must be longer than or equal to p+2
  XNN, X, Y = kNN(p+2, X, Y)
  # XNN[time, k, dim], Y[time, dim]
  #XNN = np.reshape(XNN, [XNN.shape[0],-1])
  nt  = min(XNN.shape[0], Y.shape[0])
  '''
  MI  = 0.e0
  for i in range(p):
    MI += np.mean(mutual_info_regression(XNN[:nt,:], Y[:nt,i], n_neighbors=5))
  MI /= p
  '''
  MI = 0.e0
  for k in range(p+2):
    MI += mutual_information(XNN[:nt,k,:], Y[:nt,:], k=3)
  MI /= (p+2)
  return MI


def KMeans(x, bin):
  est = KBinsDiscretizer(n_bins=bin, encode='onehot', strategy='kmeans')
  x = x.reshape(-1, 1)
  x_binned = est.fit_transform(x)
  return x_binned.indices


def Hxy(x, y, bin):
  # x and y are time-series data
  # H(y|x), H(x|y)
  if np.min(x) == np.max(x) or np.min(y) == np.max(y):
    xs = np.zeros_like(x, dtype=np.int32)
    ys = np.zeros_like(y, dtype=np.int32)
  else:
    xs = KMeans(x, bin)
    ys = KMeans(y, bin)
  Hyx = conditional_entropy(xs, ys)
  Hxy = conditional_entropy(ys, xs)
  return Hyx, Hxy


def TE(x, y, history_len, bin):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
  if np.min(x) == np.max(x) or np.min(y) == np.max(y):
    xs = np.zeros_like(x, dtype=np.int32)
    ys = np.zeros_like(y, dtype=np.int32)
  else:
    xs = KMeans(x, bin)
    ys = KMeans(y, bin)
  TEx_y = transfer_entropy(xs, ys, k=history_len)
  TEy_x = transfer_entropy(ys, xs, k=history_len)
  return TEx_y, TEy_x


def NTE(x, y, history_len, bin):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
  xs = KMeans(x, bin)
  ys = KMeans(y, bin)
  TEx_y = transfer_entropy(xs, ys, k=history_len)
  TEy_x = transfer_entropy(ys, xs, k=history_len)
  # H(Yt|Yt-1)
  Yt   = KMeans(y[1:],  bin)
  Yt_1 = KMeans(y[:-1], bin)
  Xt   = KMeans(x[1:],  bin)
  Xt_1 = KMeans(x[:-1], bin)
  Hy = conditional_entropy(Yt_1, Yt) # H(Yt|Yt-1)
  Hx = conditional_entropy(Xt_1, Xt) # H(Xt|Xt-1)
  # E[TExd_y] = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xdt-1:t-L)
  Ex_y = np.zeros(100, dtype=np.float32)
  Ey_x = np.zeros(100, dtype=np.float32)
  for i in range(100):
    x_shuffled = np.random.permutation(x)
    y_shuffled = np.random.permutation(y)
    xd = KMeans(x_shuffled, bin)
    yd = KMeans(y_shuffled, bin)
    Ex_y[i] = transfer_entropy(xd, ys, k=history_len)
    Ey_x[i] = transfer_entropy(yd, xs, k=history_len)
  # calc normalized TE
  NTEx_y = (TEx_y - np.mean(Ex_y)) / Hy
  NTEy_x = (TEy_x - np.mean(Ey_x)) / Hx
  return NTEx_y, NTEy_x


def causal_map(V, history_len, bin):
  # V[kind, time_series]
  N = len(V[:,0])
  map = np.zeros((N,N), dtype=np.float32)
  for j in range(N):
    for i in range(N):
      if i != j:
        # N[j,i]: causality from j to i
        map[j,i], map[i,j] = TE(V[j,:], V[i,:], history_len, bin)
      else:
        map[j,i] = np.nan
  return map

