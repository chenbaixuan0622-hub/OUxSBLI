import numpy as np
from scipy.special import psi, gamma
from scipy.spatial import KDTree
from mod.mod_ds import search_tau, Takens_embedding


def count_points1(tree_X, tree_XY, X, XY, k, theiler_window=1):
  N   = X.shape[0]
  eps = np.zeros(N)
  nx  = np.zeros(N)
  for i in range(N):
    distance, _ = tree_XY.query(XY[i], k=k+1, p=np.inf)
    eps[i] = 2.e0 * distance[-1]
    idx = tree_X.query_ball_point(X[i], 0.5e0 * eps[i], p=np.inf)
    idx = [j for j in idx if abs(j - i) >= theiler_window]
    nx[i] = max(len(idx), 1.e0) # nx must be greater than 0 because psi(0) = -inf
  return nx, eps


# this function doesn't work
def count_points2(tree_X, tree_Y, tree_XY, X, Y, XY, k):
  N, dx = X.shape
  eps = np.zeros(N)
  nx  = np.zeros(N)
  ny  = np.zeros(N)
  for i in range(N):
    #distance, indexes = tree_XY.query(XY[i], k=k+1)
    #Z = XY[indexes[-1]]
    #epsx = 2.e0 * np.linalg.norm(Z[:dx] - X[i])
    #epsy = 2.e0 * np.linalg.norm(Z[dx:] - Y[i])
    distancex, _ = tree_X.query(X[i], k=k+1, p=np.inf)
    distancey, _ = tree_Y.query(Y[i], k=k+1, p=np.inf)
    epsx = 2.e0 * distancex[-1]
    epsy = 2.e0 * distancey[-1]
    nx[i] = len(tree_X.query_ball_point(X[i], 0.5e0 * epsx, p=np.inf)) - 1.e0
    ny[i] = len(tree_Y.query_ball_point(Y[i], 0.5e0 * epsy, p=np.inf)) - 1.e0
  return nx, ny


def reshape_matrix(X):
  if X.ndim == 1:
    dx = 1.e0
    X  = X.reshape((-1,1))
  else:
    dx = X.shape[1]
  return X, dx


def shannon_entropy(X, k=3, theiler_window=1, Z=None):
  # Kozachenko Leonenko
  N = X.shape[0]
  X, dx = reshape_matrix(X)
  eps = np.zeros(N)
  tree_X = KDTree(X)
  # volume of unit ball in d*n
  Cdx = np.pi**(0.5e0*dx) / gamma(1.e0 + 0.5e0 * dx) / 2.e0**dx
  if Z is not None:
    Z, _ = reshape_matrix(Z)
    nx, eps = count_points1(tree_X, KDTree(Z), X, Z, k=k, theiler_window=theiler_window)
    H = np.mean(-psi(nx + 1.e0) + dx * np.log(eps)) + psi(N) + np.log(Cdx)
  else:
    for i in range(N):
      distance, _ = tree_X.query(X[i], k=k+1, p=np.inf)
      eps[i] = 2.e0 * distance[-1]
    H = -psi(k) + psi(N) + np.log(Cdx) + dx * np.mean(np.log(eps))
  return H


def joint_entropy(X, Y, k=3, theiler_window=1):
  # Kozachenko Leonenko
  N = X.shape[0]
  X, dx = reshape_matrix(X)
  Y, dy = reshape_matrix(Y)
  eps   = np.zeros(N)
  XY    = np.concatenate((X, Y), axis=-1)
  tree_XY = KDTree(XY)
  # volume of unit ball in d*n
  Cdx = np.pi**(0.5e0*dx) / gamma(1.e0 + 0.5e0 * dx) / 2.e0**dx
  Cdy = np.pi**(0.5e0*dy) / gamma(1.e0 + 0.5e0 * dy) / 2.e0**dy
  '''
  if Z is not None:
    Z, _ = reshape_matrix(Z)
    n, eps = count_points1(tree_XY, KDTree(Z), XY, Z, k=k, theiler_window=theiler_window)
    H_XY = np.mean(-psi(n + 1.e0) + (dx + dy) * np.log(eps)) + psi(N) + np.log(Cdx*Cdy)
  '''
  for i in range(N):
    distance, _ = tree_XY.query(XY[i], k=k+1, p=np.inf)
    eps[i] = 2.e0 * distance[-1]
  H_XY = -psi(k) + psi(N) + np.log(Cdx*Cdy) + (dx + dy) * np.mean(np.log(eps))
  return H_XY


def mutual_information(X, Y, k=3, type=1, theiler_window=1):
  # Kozachenko Leonenko
  N = X.shape[0]
  X, dx = reshape_matrix(X)
  Y, dy = reshape_matrix(Y)
  XY = np.concatenate((X, Y), axis=-1)
  tree_XY = KDTree(XY)
  tree_X  = KDTree(X)
  tree_Y  = KDTree(Y)
  if type == 1:
    nx, _ = count_points1(tree_X, tree_XY, X, XY, k=k, theiler_window=theiler_window)
    ny, _ = count_points1(tree_Y, tree_XY, Y, XY, k=k, theiler_window=theiler_window)
    I = psi(k) - np.mean(psi(nx) + psi(ny)) + psi(N)
  elif type == 2:
    nx, ny = count_points2(tree_X, tree_Y, tree_XY, X, Y, XY, k)
    I = psi(k) - 1.e0 / k - np.mean(psi(nx) + psi(ny)) + psi(N)
  else:
    H_X  = shannon_entropy(X, k=k, theiler_window=theiler_window, Z=XY)
    H_Y  = shannon_entropy(Y, k=k, theiler_window=theiler_window, Z=XY)
    H_XY = joint_entropy(X, Y, k=k, theiler_window=theiler_window)
    I    = H_X + H_Y - H_XY
  return I


def conditional_entropy(X, Y, k=3, theiler_window=1):
  X, _ = reshape_matrix(X)
  Y, _ = reshape_matrix(Y)
  '''
  if Z is not None:
    Z, _ = reshape_matrix(Z)
    XYZ  = np.concatenate((X, Y, Z), axis=-1)
    tree_XYZ = KDTree(XYZ)
    H_XY = joint_entropy(X, Y, k=k, theiler_window=theiler_window, Z=XYZ, tree_Z=tree_XYZ)
    H_Y  = shannon_entropy(Y,  k=k, theiler_window=theiler_window, Z=XYZ)
  else:
  '''
  XY   = np.concatenate((X, Y), axis=-1)
  tree_XY = KDTree(XY)
  H_XY = joint_entropy(X, Y, k=k, theiler_window=theiler_window)
  H_Y  = shannon_entropy(Y,  k=k, theiler_window=theiler_window, Z=XY)
  return H_XY - H_Y


def joint_entropy3(X, Y, Z, k=3, theiler_window=1):
  X, _  = reshape_matrix(X)
  Y, _  = reshape_matrix(Y)
  H_X   = shannon_entropy(X, k=k, theiler_window=theiler_window)
  H_YX  = conditional_entropy(Y, X, k=k, theiler_window=theiler_window)
  XY    = np.concatenate((X, Y), axis=-1)
  H_ZXY = conditional_entropy(Z, XY, k=k, theiler_window=theiler_window)
  H_XYZ = H_X + H_YX + H_ZXY
  return H_XYZ


def transfer_entropy(X, Y, Z, k=3, theiler_window=1):
  # informtion flow Z -> X
  X, _ = reshape_matrix(X)
  Y, _ = reshape_matrix(Y)
  Z, _ = reshape_matrix(Z)
  H_XY  = conditional_entropy(X, Y, k=k, theiler_window=theiler_window)
  H_XYZ = joint_entropy3(X, Y, Z,   k=k, theiler_window=theiler_window)
  H_YZ  = joint_entropy(Y, Z,       k=k, theiler_window=theiler_window)
  H_Z   = shannon_entropy(Z,        k=k, theiler_window=theiler_window)
  TE = H_XY - H_XYZ + H_YZ + H_Z
  return TE


def transfer_entropy_biased(X, Y, Z, k=3, theiler_window=1):
  # informtion flow Z -> X
  Z0  = np.random.permutation(Z)
  TE  = transfer_entropy(X, Y, Z,  k=k, theiler_window=theiler_window)
  TE0 = transfer_entropy(X, Y, Z0, k=k, theiler_window=theiler_window)
  return TE - TE0


def transfer_entropy_timeseries(X, Y, k=3):
  # informtion flow Y -> X
  tx = search_tau(X)
  ty = search_tau(Y)
  t  = min(tx, ty)
  TE = transfer_entropy_biased(X=X[1:], Y=X[:-1], Z=Y[:-1], k=k, theiler_window=t)
  return TE


# this functio doesn't work
def embedding_transfer_entropy(X, Y, d=2, k=3):
  # informtion flow Y -> X
  # x t
  # X t-1:t-tau*d
  # Y t-1:t-tau*d
  tx = search_tau(X)
  ty = search_tau(Y)
  t  = min(tx, ty)
  x  = np.roll(X, -1)[:-1]
  Xe = Takens_embedding(X, t, d)
  Ye = Takens_embedding(Y, t, d)
  nt = min(Xe.shape[0], Ye.shape[0])
  TE = transfer_entropy_biased(X=x[:nt], Y=Xe[:nt], Z=Ye[:nt], k=k, theiler_window=theiler_window)
  return TE

