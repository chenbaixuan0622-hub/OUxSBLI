import numpy as np
from sklearn.preprocessing import KBinsDiscretizer
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy
from pyinform.mutualinfo import mutual_info


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

