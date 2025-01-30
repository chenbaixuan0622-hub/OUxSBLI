import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_shannon_nn import shannon_entropy, mutual_information, joint_entropy, joint_entropy3, conditional_entropy, transfer_entropy_timeseries, embedding_transfer_entropy
from mod.mod_ds_test import coupling_system, discrete_logistic_map, lorenz, coupled_lorenz


def causal_map(data, k=5, trial=1):
  n   = len(data)
  Map = np.zeros((n,n))
  for itr in tqdm(range(trial)):
    for j in range(n):
      for i in range(n):
        if i == j:
          Map[j,i] = np.inf
        else:
          Map[j,i] += transfer_entropy_timeseries(data[i], data[j], k=k)
  Map /= trial
  plt.imshow(Map, cmap='jet', extent=None, origin='lower')
  plt.xlabel("Effect", fontsize=24)
  plt.ylabel("Cause", fontsize=24)
  plt.colorbar()
  plt.show()


'''
# standard deviation
sigma_x = 3.e0
sigma_y = 5.e0

x = np.random.normal(0.e0, sigma_x, 500)
y = np.random.normal(0.e0, sigma_y, 500)
'''
'''
plt.hist(x, bins=50)
plt.show()
'''
'''
# theoretical
H_x = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * sigma_x**2))
H_y = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * sigma_y**2))
# numerical
H_X = shannon_entropy(x, k=5)
H_Y = shannon_entropy(y, k=5)
print("theoretical ", H_x, "numerical ", H_X)
print("theoretical ", H_y, "numerical ", H_Y)
# numerical
H_X = shannon_entropy(x, Y=y, k=5)
H_Y = shannon_entropy(y, Y=x, k=5)
print("theoretical ", H_x, "numerical ", H_X)
print("theoretical ", H_y, "numerical ", H_Y)
'''

# variance
var_x = 9.e0
var_y = 25.e0
# covariance
cov   = 4.e0
# corr coef
rho   = cov / (np.sqrt(var_x * var_y))

mean  = (0.e0, 0.e0)
cov   = [[var_x, cov], \
         [cov, var_y]]

x = np.random.multivariate_normal(mean, cov, 10000)

'''
plt.plot(x[:,0], x[:,1], '.', alpha=0.5)
plt.axis('equal')
plt.show()
'''

# theoretical
mi = -0.5e0 * np.log(1.e0 - rho**2)
H_x = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * np.sqrt(var_x)**2))
H_y = 0.5e0 * (1.e0 + np.log(2.e0 * np.pi * np.sqrt(var_y)**2))
# numerical
H_X = shannon_entropy(x[:,0], k=5)
H_Y = shannon_entropy(x[:,1], k=5)
MI1 = mutual_information(x[:,0], x[:,1], k=5, type=1)
MI3 = mutual_information(x[:,0], x[:,1], k=5, type=3)
Hj_XY  = joint_entropy(x[:,0], x[:,1], k=5)
Hc_XY  = conditional_entropy(x[:,0], x[:,1], k=5)
Hj_XYZ = joint_entropy3(x[:,0], x[:,1], x[:,0], k=5)
print("mutual information theoretical ", mi, "numerical1", MI1, "numerical2", MI3)
print("shannon entropy theoretical ", H_x, "numerical", H_X)
print("shannon entropy theoretical ", H_y, "numerical", H_Y)
print("joint entropy theoretical", H_x + H_y - mi, "numerical", Hj_XY)
print("conditional entropy theoretical", H_x - mi, "numerical", Hc_XY)
print("joint entropy 3 theoretical", Hj_XYZ)


def test_coupling_system(nt, n, byx, x0, y0, trial):
  bxy = np.linspace(0.e0, 0.3e0, n)
  TExy = np.zeros(n)
  TEyx = np.zeros(n)
  for j in range(n):
    for i in range(trial):
      x, y = coupling_system(nt, x0, y0, bxy[j], byx)
      TEyx[j] += transfer_entropy_timeseries(x, y, k=5)
      TExy[j] += transfer_entropy_timeseries(y, x, k=5)
    TExy[j] /= trial
    TEyx[j] /= trial
    print('x->y', TExy[j], 'y->x', TEyx[j])


def test_discrete_logistic_map(nt, n, bxy, x0, y0, z0, trial):
  byz = bxy
  bxz = np.linspace(0.e0, 0.3e0, n)
  TExz = np.zeros(n)
  TEzx = np.zeros(n)
  for j in range(n):
    for i in range(trial):
      x, y, z = discrete_logistic_map(nt, x0, y0, z0, bxy, byz, bxz[j])
      TEzx[j] += transfer_entropy_timeseries(x, z, k=5)
      TExz[j] += transfer_entropy_timeseries(z, x, k=5)
    TExz[j] /= trial
    TEzx[j] /= trial
    print('x->z', TExz[j], 'z->x', TEzx[j])


def test_coupled_lorenz():
  x1, x2, x3, y1, y2, y3, z1, z2, z3 = coupled_lorenz(nt=1000, dt=0.01e0, a=10.e0, b=28.e0, c=8.e0/3.e0, bxy=0.3e0, bxz=0.2e0)
  data = [x1, x2, x3, y1, y2, y3, z1, z2, z3]
  causal_map(data, k=5, trial=10)


#print("coupling system")
#test_coupling_system(nt=500, n=31, byx=0.e0, x0=0.4e0, y0=0.6e0, trial=10)
#print("discrete logistic map")
#test_discrete_logistic_map(nt=500, n=31, bxy=0.45e0, x0=0.4e0, y0=0.5e0, z0=0.6e0, trial=10)
test_coupled_lorenz()

