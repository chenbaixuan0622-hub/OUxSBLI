import numpy as np
import lingam
import pandas as pd
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy

n = 1000

e1 = 2.e0 * (np.random.rand(n) - 0.5e0)
e2 = 2.e0 * (np.random.rand(n) - 0.5e0)
e3 = 2.e0 * (np.random.rand(n) - 0.5e0)

x = np.zeros((3, n), dtype=np.float32)

x[1,:] = e2
x[0,:] = 3.e0 * x[1,:] + e1
x[2,:] = 2.e0 * x[0,:] + 4.e0 * x[1,:] + e3

df = pd.DataFrame({"x1":x[0,:], "x2":x[1,:], "x3":x[2,:]})
df.head()

model = lingam.DirectLiNGAM()
model.fit(df)

mat_LiNGAM = np.array(model._adjacency_matrix)
print(mat_LiNGAM)


def NTE(x, y, history_len):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
  xs = bin_series(x, b=10)
  ys = bin_series(y, b=10)
  TEx_y = transfer_entropy(xs[0], ys[0], k=history_len)
  TEy_x = transfer_entropy(ys[0], xs[0], k=history_len)
  # H(Yt|Yt-1)
  Yt   = bin_series(y[1:],  b=10)
  Yt_1 = bin_series(y[:-1], b=10)
  Xt   = bin_series(x[1:],  b=10)
  Xt_1 = bin_series(x[:-1], b=10)
  Hy = conditional_entropy(Yt_1[0], Yt[0]) # H(Yt|Yt-1)
  Hx = conditional_entropy(Xt_1[0], Xt[0]) # H(Xt|Xt-1)
  # E[TExd_y] = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xdt-1:t-L)
  x_shuffled = np.random.permutation(x)
  y_shuffled = np.random.permutation(y)
  xd = bin_series(x_shuffled, b=10)
  yd = bin_series(y_shuffled, b=10)
  Ex_y = transfer_entropy(xd[0], ys[0], k=history_len)
  Ey_x = transfer_entropy(yd[0], xs[0], k=history_len)
  # calc normalized TE
  NTEx_y = (TEx_y - Ex_y) / Hy
  NTEy_x = (TEy_x - Ey_x) / Hx
  #return NTEx_y, NTEy_x
  return TEx_y, TEy_x


TE = np.zeros((3,3), dtype=np.float32)

for j in range(3):
  for i in range(3):
    TE[j,i], TE[i,j] = NTE(x[i,:], x[j,:], 5)

print(TE)

