import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import loadmat
from tqdm import tqdm
from mod.mod_shannon_nn import shannon_entropy, mutual_info, joint_entropy
from mod.mod_shannon_nn import transfer_entropy_surrogate as transfer_entropy
from mod.mod_shannon_nn import embedding_entropy_surrogate as embedding_entropy
from mod.mod_ds_test import coupling_system, discrete_logistic_map, lorenz, coupled_lorenz
from mod.mod_ds import search_tau
from mod.mod_plot import set_Params


set_Params()
path = "./data"


def causal_map(data, k=5, trial=1):
  n   = len(data)
  Map = np.zeros((n,n))
  for itr in tqdm(range(trial)):
    for j in range(n):
      for i in range(n):
        if i == j:
          Map[j,i] = np.inf
        else:
          #Map[j,i] += transfer_entropy_timeseries(data[i], data[j], k=k)
          Map[j,i] += embedding_entropy(data[i], data[j], p=6, k=5, Thei=10)
  Map /= trial
  print("x1 -> x2 ", Map[0,1], " x2 -> x1 ", Map[1,0])
  plt.imshow(Map, cmap='jet', extent=None, origin='lower')
  plt.xlabel("Effect", fontsize=24)
  plt.ylabel("Cause", fontsize=24)
  plt.colorbar()
  plt.savefig(os.path.join(path, "causal_map.png"))
  plt.show()


def test_shannon_mutual():
  sigma = np.linspace(1.e0, 5.e0, 50)
  H_t = np.zeros(len(sigma))
  H_n = np.zeros(len(sigma))
  for i, s in enumerate(sigma):
    x = np.random.normal(0.e0, s, 500)
    H_t[i] = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * s**2))
    H_n[i] = shannon_entropy(x, k=5)

  plt.figure(figsize=(6,6))
  plt.plot(sigma, H_t, color='black')
  plt.plot(sigma, H_n, "o", color='blue')
  plt.xlabel(r'$\sigma$', fontsize=18, style='italic')
  plt.ylabel("H(X)", fontsize=18)
  plt.savefig(os.path.join(path, "H_sigma.png"))
  plt.show()


  K = np.arange(1,50)
  s = 1.e0
  x = np.random.normal(0.e0, s, 500)
  H_t = np.zeros(len(K))
  H_n = np.zeros(len(K))
  for i, k in enumerate(K):
    H_t[i] = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * s**2))
    H_n[i] = shannon_entropy(x, k=k)

  plt.figure(figsize=(6,6))
  plt.plot(K, H_t, color='black')
  plt.plot(K, H_n, "o", color='blue')
  plt.xlabel("k", fontsize=18, style='italic')
  plt.ylabel("H(X)", fontsize=18)
  plt.savefig(os.path.join(path, "H_k.png"))
  plt.show()


  # variance
  var_x = 9.e0
  var_y = 25.e0
  covariance = np.linspace(0.e0, 10.e0, 50)
  mean = (0.e0, 0.e0)
  H_xy = np.zeros(len(covariance))
  H_XY = np.zeros(len(covariance))
  I    = np.zeros(len(covariance))
  for i, cov in enumerate(covariance):
    rho  = cov / (np.sqrt(var_x * var_y)) # corr coef
    Cov  = [[var_x, cov], \
            [cov, var_y]]
    x = np.random.multivariate_normal(mean, Cov, 10000)
    mi  = -0.5e0 * np.log(1.e0 - rho**2)
    H_x = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * np.sqrt(var_x)**2))
    H_y = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * np.sqrt(var_y)**2))
    H_xy[i] = H_x + H_y - mi
    H_XY[i] = joint_entropy(x[:,0], x[:,1], k=5)
    I[i] = mutual_info(x[:,0].reshape(-1,1), x[:,1].reshape(-1,1), k=10, Thei=1)
    print(mi, I[i])

  plt.figure(figsize=(6,6))
  plt.plot(covariance, H_xy, color='black')
  plt.plot(covariance, H_XY, "o", color='blue')
  plt.xlabel(r'$\sigma_{xy}$', fontsize=18, style='italic')
  plt.ylabel("H(X,Y)", fontsize=18)
  plt.savefig(os.path.join(path, "H_XY.png"))
  plt.show()


def test_coupling_system(nt, n, byx, x0, y0, trial):
  bxy = np.linspace(0.e0, 0.3e0, n)
  TExy = np.zeros(n)
  TEyx = np.zeros(n)
  EExy = np.zeros(n)
  EEyx = np.zeros(n)
  for j in tqdm(range(n)):
    for i in range(trial):
      x, y = coupling_system(nt, x0, y0, bxy[j], byx)
      TEyx[j] += transfer_entropy(x, y, k=5, Thei=10)
      TExy[j] += transfer_entropy(y, x, k=5, Thei=10)
      EExy[j] += embedding_entropy(y, x, p=1, k=5, Thei=10)
      EEyx[j] += embedding_entropy(x, y, p=1, k=5, Thei=10)
    TExy[j] /= trial
    TEyx[j] /= trial
    EExy[j] /= trial
    EEyx[j] /= trial
  plt.figure(figsize=(6,6))
  plt.plot(bxy, TExy, color='blue',  linestyle='solid')
  plt.plot(bxy, TEyx, color='blue',  linestyle='dashed')
  plt.plot(bxy, EExy, color='red', linestyle='solid')
  plt.plot(bxy, EEyx, color='red', linestyle='dashed')
  plt.xlabel(r'$\beta_{xy}$', fontsize=18, style='italic')
  plt.ylabel("Causality", fontsize=18)
  plt.ylim([-0.1, 1])
  plt.savefig(os.path.join(path, "coupling_system.png"))
  #plt.show()

def test_discrete_logistic_map(nt, n, bxy, x0, y0, z0, trial):
  byz = bxy
  bxz = np.linspace(0.e0, 0.3e0, n)
  TExz = np.zeros(n)
  TEzx = np.zeros(n)
  EExz = np.zeros(n)
  EEzx = np.zeros(n)
  for j in tqdm(range(n)):
    for i in range(trial):
      x, y, z = discrete_logistic_map(nt, x0, y0, z0, bxy, byz, bxz[j])
      TEzx[j] += transfer_entropy(x, z, k=5, Thei=10)
      TExz[j] += transfer_entropy(z, x, k=5, Thei=10)
      EEzx[j] += embedding_entropy(x, z, p=6, k=5, Thei=10)
      EExz[j] += embedding_entropy(z, x, p=6, k=5, Thei=10)
    TExz[j] /= trial
    TEzx[j] /= trial
    EExz[j] /= trial
    EEzx[j] /= trial
  plt.figure(figsize=(6,6))
  plt.plot(bxz, TExz, color='blue',  linestyle='solid')
  plt.plot(bxz, TEzx, color='blue',  linestyle='dashed')
  plt.plot(bxz, EExz, color='red', linestyle='solid')
  plt.plot(bxz, EEzx, color='red', linestyle='dashed')
  plt.xlabel(r'$\beta_{xz}$', fontsize=18, style='italic')
  plt.ylabel("Causality", fontsize=18)
  plt.ylim([-0.1, 1])
  plt.savefig(os.path.join(path, "discrete_logstic_map.png"))
  plt.show()


def test_coupled_lorenz():
  x1, x2, x3, y1, y2, y3, z1, z2, z3 = coupled_lorenz(nt=5000, dt=0.01e0, a=10.e0, b=28.e0, c=8.e0/3.e0, bxy=0.3e0, bxz=0.2e0)
  tau  = 1
  data = [x1[::tau], x2[::tau], x3[::tau], y1[::tau], y2[::tau], y3[::tau], z1[::tau], z2[::tau], z3[::tau]]
  causal_map(data, k=5, trial=10)


def test_tbl(dir):
  data1 = loadmat(os.path.join(dir, 'Inner_outer_u_z32_c1.mat'))['data']
  data2 = loadmat(os.path.join(dir, 'Inner_outer_u_z32_c2.mat'))['data']
  data3 = loadmat(os.path.join(dir, 'Inner_outer_u_z32_c3.mat'))['data']
 
  #tau = min(search_tau(data1[:1000,0]), search_tau(data1[:1000,1]))
  #print("delay: ", tau)

  #nt  = 50000
  '''
  # check data
  plt.plot(data1[:nt:tau,0], color='red')
  plt.plot(data2[:nt:tau,1], color='blue')
  plt.show()
  '''

  TEoi = transfer_entropy(data1[:,1], data1[:,0], k=5, Thei=10)
  TEio = transfer_entropy(data1[:,0], data1[:,1], k=5, Thei=10)
  print("outer -> inner:", TEoi, " inner -> outer:", TEio)


#test_shannon_mutual()
test_coupling_system(nt=500, n=31, byx=0.1e0, x0=0.4e0, y0=0.6e0, trial=100)
test_discrete_logistic_map(nt=1000, n=31, bxy=0.4e0, x0=0.4e0, y0=0.5e0, z0=0.6e0, trial=100)
#test_coupled_lorenz()
#test_tbl('../../SURD_data')

