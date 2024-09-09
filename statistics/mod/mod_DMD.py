import numpy as np
import scipy as sp

def rank_truncation(r, t, X, Y, U2, Sig2, Vh2):
  dt  = -t[0] + t[1]
  U   = U2[:,:r]
  Sig = np.diag(Sig2)[:r,:r]
  V   = Vh2.conj().T[:,:r]
  # build A tilde
  Anti  = np.dot(np.dot(np.dot(U.conj().T, Y), V), np.linalg.inv(Sig))
  mu, W = np.linalg.eig(Anti)
  # build DMD modes
  Phi   = np.dot(np.dot(np.dot(Y, V), np.linalg.inv(Sig)), W)
  # compute time evolution
  b   = np.dot(np.linalg.pinv(Phi), X[:,0])
  Psi = np.zeros([r, len(t)], dtype='complex')
  for i, _t in enumerate(t):
    Psi[:,i] = np.multiply(np.power(mu, _t / dt), b)

  # compute DMD reconstruction
  D2 = np.dot(Phi, Psi)
  return Psi, D2

def DMD(t, D, rank):
  # t[time]
  # D[spcae, time]
  X = D[:,:-1]
  Y = D[:,1:]
  # SVD
  U2, Sig2, Vh2 = sp.linalg.svd(X, False)
  Psi, D2       = rank_truncation(rank, t, X, Y, U2, Sig2, Vh2)
  return Psi, D2

