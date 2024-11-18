import numpy as np
from sklearn.preprocessing import KBinsDiscretizer
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy


def cross_corr(x, y):
  X = np.mean(x)
  Y = np.mean(y)
  corr = np.mean((x - X) * (y - Y)) \
         / (np.sqrt(np.mean((x - X)**2)) * np.sqrt(np.mean((y - Y)**2)))
  return corr


def KMeans(x, bin):
  est = KBinsDiscretizer(n_bins=bin, encode='onehot', strategy='kmeans')
  x = x.reshape(-1, 1)
  x_binned = est.fit_transform(x)
  return x_binned.indices


def TE(x, y, history_len, bin):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
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

