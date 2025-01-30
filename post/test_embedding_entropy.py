import numpy as np
import matplotlib.pyplot as plt
from mod.mod_plot import set_Params
from mod.mod_embedding_entropy import EE
from mod.mod_ds_test import coupling_system, discrete_logistic_map, lorenz


set_Params()


def causal_map(data, p, k, theiler_window):
  n   = len(data)
  Map = np.zeros((n,n))
  for j in range(n):
    for i in range(n):
      if i == j:
        Map[j,i] = np.inf
      else:
        Map[j,i] = EE(x=data[i], y=data[j], p=p, k=k, theiler_window=theiler_window)

  plt.imshow(Map, cmap='jet', extent=None, origin='lower')
  plt.xlabel("Effect", fontsize=24)
  plt.ylabel("Cause", fontsize=24)
  #plt.xticks([0,1,2,3], ['$\it{x}$', '$\it{y}$', '$\it{z}$', '$\it{\epsilon}$'], fontsize=20)
  #plt.yticks([0,1,2,3], ['$\it{x}$', '$\it{y}$', '$\it{z}$', '$\it{\epsilon}$'], fontsize=20)
  plt.colorbar()
  plt.show()


def test_coupling_system(nt, n, byx, x0, y0, trial):
  # coupling system
  bxy = np.linspace(0.e0, 0.3e0, n)

  EExy = np.zeros(n)
  EEyx = np.zeros(n)

  for j in range(n):
    corr = 0.e0
    for i in range(trial):
      x, y = coupling_system(nt, x0, y0, bxy[j], byx)
      EEyx[j] += EE(x=x, y=y, p=1, k=5, theiler_window=10)
      EExy[j] += EE(x=y, y=x, p=1, k=5, theiler_window=10)
      corr    += np.corrcoef(x, y)[0,1]
    EExy[j] /= trial
    EEyx[j] /= trial
    corr    /= trial
    print('y->x', EEyx[j], 'x->y', EExy[j], 'corr', corr)

  plt.plot(bxy, np.zeros_like(bxy), color='black')
  plt.plot(bxy, EExy, color='blue')
  plt.plot(bxy, EEyx, color='red')
  plt.plot(bxy, EExy - EEyx, color='green')
  plt.savefig("byx00.png")
  plt.show()


def test_discrete_logistic_map(nt, n, bxy, x0, y0, z0, trial):
  byz = bxy
  '''
  bxz = np.linspace(0.e0, 0.3e0, n)
  EExz = np.zeros(n)
  EEzx = np.zeros(n)

  for j in range(n):
    corr = 0.e0
    for i in range(trial):
      x, y, z = discrete_logistic_map(nt, x0, y0, z0, bxy, byz, bxz[j])
      EEzx[j] += EE(x=x, y=z, p=6, k=5, theiler_window=10)
      EExz[j] += EE(x=z, y=x, p=6, k=5, theiler_window=10)
      corr    += np.corrcoef(z, x)[0,1]
    EExz[j] /= trial
    EEzx[j] /= trial
    corr    /= trial
    print('x->z', EExz[j], 'z->x', EEzx[j], 'corr', corr)

  plt.plot(bxz, EExz, color='blue')
  plt.plot(bxz, EEzx, color='red')
  plt.savefig("bxy02.png")
  plt.show()
  ''' 
  x, y, z = discrete_logistic_map(nt, x0, y0, z0, bxy, byz, bxz=0.2e0)
  data = [x, y, z, np.random.permutation(x)]
  causal_map(data, p=6, k=5, theiler_window=10)


def test_lorenz_system():
  x, y, z = lorenz(nt=1000, dt=0.01e0, x0=0.e0, y0=1.e0, z0=1.05e0, s=10.e0, r=28.e0, b=8.e0/3.e0)
  data = [x[::2], y[::2], z[::2], np.random.permutation(x[::2])]
  causal_map(data, p=6, k=5, theiler_window=10)


#test_coupling_system(nt=500, n=31, byx=0.e0, x0=0.4e0, y0=0.6e0, trial=10)
test_discrete_logistic_map(nt=500, n=31, bxy=0.45e0, x0=0.4e0, y0=0.5e0, z0=0.6e0, trial=1)
#test_lorenz_system()

