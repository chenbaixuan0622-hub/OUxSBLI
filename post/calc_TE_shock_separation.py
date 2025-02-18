import numpy as np
import os
import matplotlib.pyplot as plt
from mod.mod_shannon_nn import transfer_entropy_surrogate as transfer_entropy
from mod.mod_plot import set_Params


set_Params()
dir = "../../SBLI_data"

path   = os.path.join(dir, "39.17", "x_non_dimensional_y=5.dat")
shock5 = np.loadtxt(path)
path   = os.path.join(dir, "39.17", "x_non_dimensional_y=6.dat")
shock6 = np.loadtxt(path)
path   = os.path.join(dir, "39.17", "x_non_dimensional_y=7.dat")
shock7 = np.loadtxt(path)
path   = os.path.join(dir, "separate_position", "separate_position_39.17.txt")
sepa   = np.loadtxt(path)

shock5[:,1] *= 1e-3
shock6[:,1] *= 1e-3
shock7[:,1] *= 1e-3

nt = 2000
t  = shock5[nt:,0]

data = np.zeros((4,len(t)))
data[0] = shock5[nt:,1] - np.mean(shock5[nt:,1])
data[1] = shock6[nt:,1] - np.mean(shock6[nt:,1])
data[2] = shock7[nt:,1] - np.mean(shock7[nt:,1])
data[3] = sepa[nt:,1]   - np.mean(sepa[nt:,1])

plt.plot(t, data[0], color='blue')
plt.plot(t, data[1], color='green')
plt.plot(t, data[2], color='red')
plt.plot(t, data[3], color='black')
plt.show()


def causal_map(data, k=5, trial=1):
  n   = len(data)
  Map = np.zeros((n,n))
  for j in range(n):
    for i in range(n):
      if i == j:
        Map[j,i] = np.inf
      else:
        Map[j,i] += transfer_entropy(data[i], data[j], k=k)
  Map /= trial
  print("3 to 0:", Map[3,0], " 0 to 3:", Map[0,3])
  plt.imshow(Map, cmap='jet', extent=None, origin='lower')
  plt.xlabel("Effect", fontsize=24)
  plt.ylabel("Cause", fontsize=24)
  plt.colorbar()
  plt.show()


causal_map(data, k=5, trial=1)

