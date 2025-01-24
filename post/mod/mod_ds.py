import numpy as np
from scipy.spatial.distance import euclidean
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.feature_selection import mutual_info_regression
from pyinform.mutualinfo import mutual_info
import matplotlib.pyplot as plt
from tqdm import tqdm


def lorenz(x, y, z, s=10.e0, r=28.e0, b=8.e0/3.e0):
  dot_x = s * (y - x)
  dot_y = r * x - y - x * z
  dot_z = x * y - b * z
  return dot_x, dot_y, dot_z


def case1(nt):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  x[0] = 0.1e0#np.random.rand(1)
  y[0] = 0.2e0#np.random.rand(1)
  for i in range(nt):
    x[i+1] = 3.81e0 * x[i] * (1.e0 - x[i])
    y[i+1] = 3.82e0 * y[i] * (1.e0 - y[i])
  return x, y


def case2(nt, w):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  x[0] = 0.1e0#np.random.rand(1)
  y[0] = 0.2e0#np.random.rand(1)
  for i in range(nt):
    x[i+1] = 3.81e0 * x[i] * (1.e0 - x[i])
    y[i+1] = (1.e0 - w) * 3.82e0 * y[i] * (1.e0 - y[i]) \
            + w * 3.81e0 * x[i] * (1.e0 - x[i])
  return x, y


def case3(nt, w):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  x[0] = 0.1e0#np.random.rand(1)
  y[0] = 0.2e0#np.random.rand(1)
  for i in range(nt):
    x[i+1] = (1.e0 - w) * 3.81e0 * x[i] * (1.e0 - x[i]) \
            + w * 3.82e0 * y[i] * (1.e0 - y[i])
    y[i+1] = (1.e0 - w) * 3.82e0 * y[i] * (1.e0 - y[i]) \
            + w * 3.81e0 * x[i] * (1.e0 - x[i])
  return x, y


def case4(nt, w):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  z = np.zeros(nt+1)
  x[0] = 0.1e0#np.random.rand(1)
  y[0] = 0.2e0#np.random.rand(1)
  z[0] = 0.2e0#np.random.rand(1)
  for i in range(nt):
    z[i+1] = 3.8e0 * z[i] * (1.e0 - z[i])
    x[i+1] = (1.e0 - w) * 3.81e0 * x[i] * (1.e0 - x[i]) \
            + w * 3.82e0 * z[i] * (1.e0 - z[i])
    y[i+1] = (1.e0 - w) * 3.82e0 * y[i] * (1.e0 - y[i]) \
            + w * 3.81e0 * z[i] * (1.e0 - z[i])
  return x, y, z


def coupling_system(nt, x0, y0, bxy, byx, gx=3.7e0, gy=3.72e0):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  x[0] = x0
  y[0] = y0
  for i in range(nt):
    epsx = np.random.normal(0.e0, 0.01e0)
    epsy = np.random.normal(0.e0, 0.01e0)
    x[i+1] = x[i] * (gx - (gx - byx) * x[i] - byx * y[i]) + epsx
    y[i+1] = y[i] * (gy - (gy - bxy) * y[i] - bxy * x[i]) + epsy
  return x, y


def discrete_logistic_map(nt, x0, y0, z0, bxy, byz, bxz, gx=3.7e0, gy=3.72e0, gz=3.78e0):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  z = np.zeros(nt+1)
  x[0] = x0
  y[0] = y0
  z[0] = z0
  for i in range(nt):
    epsx = np.random.normal(0.e0, 0.005e0)
    epsy = np.random.normal(0.e0, 0.005e0)
    epsz = np.random.normal(0.e0, 0.005e0)
    x[i+1] = gx * x[i] * (1.e0 - x[i]) + epsx
    y[i+1] = gy * y[i] * (1.e0 - (1.e0 - bxy / gy) * y[i] - bxy * x[i] / gy) + epsy
    z[i+1] = gz * z[i] * (1.e0 - (1.e0 - (bxz + byz) / gz) * z[i] - bxz * x[i] / gz - byz * y[i] / gz) + epsz
  return x, y, z
  

def Takens_embedding(x, tau, dim):
  xe = np.zeros((len(x)-(dim-1),dim), dtype=np.float64)
  xe[:,0] = x[:-(dim-1)]
  for i in range(1,dim):
    xe[:,i] = np.roll(x, -tau*i)[:-(dim-1)]
  return xe


def search_tau(time_series, plot=False):
  n = len(time_series)
  max_lag = n//2
  mutual_info_values = []

  for lag in range(1, max_lag + 1):
    x  = time_series[:-lag].reshape(-1,1)
    y  = time_series[lag:]
    mi = mutual_info_regression(x, y, n_neighbors=5)
    mutual_info_values.append(mi[0])
  
  mutual_info_values = np.array(mutual_info_values)
  dmi = -mutual_info_values[:-1] + mutual_info_values[1:]
  for lag in range(len(dmi)):
    if dmi[lag] >= 0.e0:
      tau = lag + 1
      break

  if plot:
    plt.plot(range(1, max_lag+1), mutual_info_values)
    plt.xlabel("Lag")
    plt.ylabel("MI")
    plt.show()
  return tau


def KMeans(x, bin):
  est = KBinsDiscretizer(n_bins=bin, encode='onehot', strategy='kmeans')
  x = x.reshape(-1, 1)
  x_binned = est.fit_transform(x)
  return x_binned.indices


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

