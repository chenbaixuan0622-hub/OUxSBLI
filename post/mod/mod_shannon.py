import numpy as np
from scipy.special import psi, gamma
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
from mod.mod_ds import search_tau


def count_points1(tree_X, tree_XY, X, XY, k, Thei=1):
  N   = X.shape[0]
  eps = np.zeros(N)
  nx  = np.zeros(N)
  for i in range(N):
    distance, _ = tree_XY.query(XY[i], k=k+1, p=np.inf)
    eps[i] = 2.e0 * distance[-1]
    idx = tree_X.query_ball_point(X[i], 0.5e0 * eps[i], p=np.inf)
    idx = [j for j in idx if abs(j - i) >= Thei]
    nx[i] = max(len(idx), 1.e0) # nx must be greater than 0 because psi(0) = -inf
  return nx, eps


def reshape_matrix(X):
  if X.ndim == 1:
    dx = 1.e0
    X  = X.reshape((-1,1))
  else:
    dx = X.shape[1]
  return X, dx


def shannon_entropy(X, k=3, Thei=1, Z=None):
  # Kozachenko Leonenko
  N = X.shape[0]
  X, dx = reshape_matrix(X)
  eps = np.zeros(N)
  tree_X = KDTree(X)
  # volume of unit ball in d*n
  Cdx = np.pi**(0.5e0*dx) / gamma(1.e0 + 0.5e0 * dx) / 2.e0**dx
  if Z is not None:
    Z, _ = reshape_matrix(Z)
    nx, eps = count_points1(tree_X, KDTree(Z), X, Z, k=k, Thei=Thei)
    H = np.mean(-psi(nx + 1.e0) + dx * np.log(eps)) + psi(N) + np.log(Cdx)
  else:
    for i in range(N):
      distance, _ = tree_X.query(X[i], k=k+1, p=np.inf)
      eps[i] = 2.e0 * distance[-1]
    H = -psi(k) + psi(N) + np.log(Cdx) + dx * np.mean(np.log(eps))
  return H


def mutual_info(X, Y, k=5, Thei=10):
  N = X.shape[0]
  nX, nY, nN = np.zeros(N), np.zeros(N), np.zeros(N)
  
  for i in range(N):
    idx = np.ones(N, dtype=bool)
    idx[max(0, i - Thei):min(i + Thei + 1, N)] = False
    
    tree_XY = KDTree(np.hstack((X[idx], Y[idx])))
    dist, _ = tree_XY.query(np.hstack((X[i], Y[i])), k=k, p=np.inf, workers=-1)
    half_epsilon_XYkNN = dist[-1]
    
    nX[i] = np.sum(cdist(X[idx], X[i].reshape(1,-1), metric='chebyshev') < half_epsilon_XYkNN)
    nY[i] = np.sum(cdist(Y[idx], Y[i].reshape(1,-1), metric='chebyshev') < half_epsilon_XYkNN)
    nN[i] = np.sum(idx)
  
  valid_idx = (nX > 0) & (nY > 0)
  I = psi(k) - np.mean(psi(nX[valid_idx] + 1)) - np.mean(psi(nY[valid_idx] + 1)) + np.mean(psi(nN[valid_idx] + 1))
  return I 


def joint_entropy(X, Y, k=3, Thei=1):
  X, _  = reshape_matrix(X)
  Y, _  = reshape_matrix(Y)
  H_X = shannon_entropy(X, k=k, Thei=Thei)
  H_Y = shannon_entropy(Y, k=k, Thei=Thei)
  I   = mutual_info(X, Y, k=k, Thei=Thei)
  H_XY = H_X + H_Y - I
  return H_XY


def conditional_entropy(X, Y, k=3, Thei=1):
  X, _ = reshape_matrix(X)
  Y, _ = reshape_matrix(Y)
  XY   = np.concatenate((X, Y), axis=-1)
  tree_XY = KDTree(XY)
  H_XY = joint_entropy(X, Y, k=k, Thei=Thei)
  H_Y  = shannon_entropy(Y,  k=k, Thei=Thei, Z=XY)
  return H_XY - H_Y


def joint_entropy3(X, Y, Z, k=3, Thei=1):
  X, _  = reshape_matrix(X)
  Y, _  = reshape_matrix(Y)
  H_X   = shannon_entropy(X, k=k, Thei=Thei)
  H_YX  = conditional_entropy(Y, X, k=k, Thei=Thei)
  XY    = np.concatenate((X, Y), axis=-1)
  H_ZXY = conditional_entropy(Z, XY, k=k, Thei=Thei)
  H_XYZ = H_X + H_YX + H_ZXY
  return H_XYZ

'''
def transfer_entropy(X, Y, Z, k=3, Thei=1):
  # informtion flow Z -> X
  X, _ = reshape_matrix(X)
  Y, _ = reshape_matrix(Y)
  Z, _ = reshape_matrix(Z)
  H_XY  = conditional_entropy(X, Y, k=k, Thei=Thei)
  H_XYZ = joint_entropy3(X, Y, Z,   k=k, Thei=Thei)
  H_YZ  = joint_entropy(Y, Z,       k=k, Thei=Thei)
  H_Z   = shannon_entropy(Z,        k=k, Thei=Thei)
  TE = H_XY - H_XYZ + H_YZ + H_Z
  return TE


def transfer_entropy_biased(X, Y, Z, k=3, Thei=1):
  # informtion flow Z -> X
  Z0  = np.random.permutation(Z)
  TE  = transfer_entropy(X, Y, Z,  k=k, Thei=Thei)
  TE0 = transfer_entropy(X, Y, Z0, k=k, Thei=Thei)
  return TE - TE0


def transfer_entropy_timeseries(X, Y, k=3):
  # informtion flow Y -> X
  tx = search_tau(X)
  ty = search_tau(Y)
  t  = min(tx, ty)
  TE = transfer_entropy_biased(X=X[1:], Y=X[:-1], Z=Y[:-1], k=k, Thei=t)
  return TE
'''

def transfer_entropy(x, y, p=1, tau=1, k=2, Thei=None):
  '''
  x   : ndarray, shape (T) time series
  y   : ndarray, shape (T) time series
  p   : int, order of the model to estimate causality
  tau : int, time delay for embedding
  k   : int, k-th nearest neighbor number
  Thei: int, half-length of Theiler correction window
  out : float, embedding entropy x->y
  '''
 
  if Thei is None or Thei < p * tau:
    Thei = p * tau
  
  if y.ndim == 1:
    dy = 1
    y  = y.reshape(1,-1)
  if x.ndim == 1:
    dx = 1
    x  = x.reshape(1,-1)
  T = y.shape[1]

  Y = y[:,p*tau:T].T
  X = np.hstack([x[:,(p-i)*tau:T-i*tau].T for i in range(1, p + 1)])
  Z = np.hstack([y[:,(p-i)*tau:T-i*tau].T for i in range(1, p + 1)])

  N = T - p * tau

  nZ, nYZ, nXZ = np.zeros(N), np.zeros(N), np.zeros(N)

  tree_XYZ = KDTree(np.hstack((X, Y, Z)))
  tree_Z   = KDTree(Z)
  tree_YZ  = KDTree(np.hstack((Y, Z)))
  tree_XZ  = KDTree(np.hstack((X, Z)))

  for i in range(N):
    idx = np.ones(N, dtype=bool)
    idx[max(0, i - Thei):min(i + Thei + 1, N)] = False

    temp_XYZ = np.hstack((X[idx], Y[idx], Z[idx]))
    temp_Z   = Z[idx]
    temp_YZ  = np.hstack((Y[idx], Z[idx]))
    temp_XZ  = np.hstack((X[idx], Z[idx]))

    dist, _ = tree_XYZ.query(np.hstack((X[i], Y[i], Z[i])), k=k, p=np.inf)
    half_epsilon_XYZ = dist[-1]

    nZ[i]  = np.sum(cdist(temp_Z, Z[i].reshape(1,-1), metric='chebyshev') < half_epsilon_XYZ)
    nYZ[i] = np.sum(cdist(temp_YZ, np.hstack((Y[i], Z[i])).reshape(1,-1), metric='chebyshev') < half_epsilon_XYZ)
    nXZ[i] = np.sum(cdist(temp_XZ, np.hstack((X[i], Z[i])).reshape(1,-1), metric='chebyshev') < half_epsilon_XYZ)

  valid_idx = (nZ > 0) & (nYZ > 0) & (nXZ > 0)
  TE = psi(k) + np.mean(psi(nZ[valid_idx] + 1)) - np.mean(psi(nXZ[valid_idx] + 1)) - np.mean(psi(nYZ[valid_idx] + 1))
  return TE


def transfer_entropy_surrogate(x, y, p=1, tau=1, k=5, Thei=10):
  # informtion flow y -> x
  x0  = np.random.permutation(x)
  TE  = transfer_entropy(y, x,  p=p, tau=tau, k=k, Thei=Thei)
  TE0 = transfer_entropy(y, x0, p=p, tau=tau, k=k, Thei=Thei)
  return TE - TE0


def embedding_entropy(x, y, p=1, tau=1, k=5, Thei=None):
  '''
  x   : ndarray, shape (T) time series
  y   : ndarray, shape (T) time series
  p   : int, order of the model to estimate causality
  tau : int, time delay for embedding
  k   : int, k-th nearest neighbor number
  Thei: int, half-length of Theiler correction window
  out : float, embedding entropy x->y
  '''
  
  if Thei is None or Thei < p * tau:
    Thei = p * tau
  
  if x.ndim == 1:
    dx = 1
    x  = x.reshape(1,-1)
  else:
    dx = x.shape[0]
  if y.ndim == 1:
    dy = 1
    y  = y.reshape(1,-1)
  else:
    dy = y.shape[0]
 
  T = x.shape[1]

  X = np.hstack([x[:,(p-i)*tau:T-i*tau].T for i in range(p + 1)])
  Y = np.hstack([y[:,(p-i)*tau:T-i*tau].T for i in range(1, p + 1)])
  N = T - p * tau

  XNN = np.zeros((N, X.shape[1] * (dx * (p + 1) + 1)))
  
  for i in range(N):
    idx = np.ones(N, dtype=bool)
    idx[max(0, i - Thei):min(i + Thei + 1, N)] = False
    idx[i] = True
    
    temp_X = X[idx]
    tree_X = KDTree(temp_X)
    
    _, pnn_idx = tree_X.query(X[i], k=p + 2, workers=-1)
    XNN[i] = temp_X[pnn_idx].flatten()
  
  EE = mutual_info(Y[:N], XNN[:N], k=k, Thei=Thei)
  return EE


def embedding_entropy_surrogate(x, y, p=1, tau=1, k=5, Thei=10):
  # informtion flow y -> x
  y0  = np.random.permutation(y)
  EE  = embedding_entropy(x, y,  p=p, tau=tau, k=k, Thei=Thei)
  EE0 = embedding_entropy(x, y0, p=p, tau=tau, k=k, Thei=Thei)
  return EE - EE0

