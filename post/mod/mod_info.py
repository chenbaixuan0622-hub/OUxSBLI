import numpy as np
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.feature_selection import mutual_info_regression
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy
from pyinform.mutualinfo import mutual_info
from mod.mod_ds import Takens_embedding, search_tau, kNN


def cross_corr(x, y):
  X = np.mean(x)
  Y = np.mean(y)
  corr = np.mean((x - X) * (y - Y)) \
         / (np.sqrt(np.mean((x - X)**2)) * np.sqrt(np.mean((y - Y)**2)))
  return corr


def X_embedding(x, p, tau=None, bin=10):
  if tau is None:
    tau, _ = search_tau(x, bin, len(x))

  X = np.zeros((len(x)-tau*p,p+1), dtype=np.float32)
  for i in range(p+1):
    X[:,i] = np.roll(x, tau*i)[tau*p:]
  return X


def Y_embedding(y, p, tau=None, bin=10):
  if tau is None:
    tau, _ = search_tau(y, bin, len(y))
    Y = X_embedding(y, p-1, tau, bin=10)
  else:
    Y = X_embedding(y, p-1, tau)
  return np.delete(Y, 0, axis=0)


def EE(x, y, p, bin=10):
  # x[time], y[time]
  # p: dimension
  X = X_embedding(x, p, bin)
  Y = Y_embedding(y, p, bin)

  # X, Y: embedded time series data
  # X[time, dim], Y[time, dim]

  # length of X must be longer than or equal to p+2
  XNN, X, Y = kNN(p+2, X, Y)
  # XNN[time, k, dim], Y[time, dim]
  XNN = np.reshape(XNN, [XNN.shape[0],-1])
  nt  = min(XNN.shape[0], Y.shape[0])
  MI  = 0.e0
  for i in range(p):
    MI += np.mean(mutual_info_regression(XNN[:nt,:], Y[:nt,i]))
  MI /= p
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

