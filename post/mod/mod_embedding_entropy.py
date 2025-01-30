import numpy as np
from scipy.spatial import KDTree
from mod.mod_ds import search_tau
from mod.mod_shannon_nn import mutual_information


def X_embedding(x, p, tau=None):
  if tau is None:
    tau = search_tau(x)
  X = np.zeros((len(x)-tau*p,p+1))
  for i in range(p+1):
    X[:,i] = np.roll(x, tau*i)[tau*p:]
  return X, tau


def Y_embedding(y, p, tau=None):
  if tau is None:
    tau = search_tau(y)
    Y, _ = X_embedding(y, p-1, tau)
  else:
    Y, _ = X_embedding(y, p-1, tau)
  return np.delete(Y, 0, axis=0)


def EE(x, y, p, k=5, theiler_window=None):
  # causality from y to x
  # x[time], y[time], p: dimension
  X, tau = X_embedding(x, p)
  Y      = Y_embedding(y, p)

  if theiler_window is not None:
    tau = theiler_window

  nt, d = X.shape
  XNN   = np.zeros((nt, p+2, d))
  tree  = KDTree(X)
  for i in range(nt):
    _, indexes = tree.query(X[i], k=p+2)
    XNN[i] = X[indexes]
  nt = min(XNN.shape[0], Y.shape[0])
  
  MI = 0.e0
  for knn in range(p+2):
    MI += mutual_information(XNN[:nt,knn,:], Y[:nt,:], k=k, type=1, theiler_window=tau)
  MI /= (p+2)
  return MI

