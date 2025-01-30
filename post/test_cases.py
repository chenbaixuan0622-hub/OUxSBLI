import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_embedding_entropy import EE
from mod.mod_ds_test import coupling_system, discrete_logistic_map

set_Params()


w   = 0.1e0 

nt  = 1000

# coupling system
n   = 41
byx = 0.e0
bxy = np.linspace(0.e0, 0.40e0, n)

x0  = 0.5e0
y0  = 0.5e0
z0  = 0.5e0

EExy = np.zeros(n)
EEyx = np.zeros(n)

trial = 1
 
for j in range(n):
  for i in range(trial):
    x, y = coupling_system(nt, x0, y0, bxy[j], byx)
    EEyx[j] += EE(x=x, y=y, p=5, theiler_window=5)
    EExy[j] += EE(x=y, y=x, p=5, theiler_window=5)
  EExy[j] /= trial
  EEyx[j] /= trial
  print('y->x', EEyx[j], 'x->y', EExy[j])

plt.plot(bxy, np.zeros_like(bxy), color='black')
plt.plot(bxy, EExy, color='blue')
plt.plot(bxy, EEyx, color='red')
plt.plot(bxy, EExy - EEyx, color='green')
plt.savefig("byx00.png")
plt.show()

'''
byx = 0.1e0

EExy = np.zeros(n)
EEyx = np.zeros(n)

for j in range(n):
  for i in range(trial):
    x, y = coupling_system(nt, x0, y0, bxy[j], byx)
    EExy[j] += EE(x, y, p=5)
    EEyx[j] += EE(y, x, p=5)
  EExy[j] /= trial
  EEyx[j] /= trial
  print(EExy[j], EEyx[j])

plt.plot(bxy, np.zeros_like(bxy), color='black')
plt.plot(bxy, EExy, color='blue')
plt.plot(bxy, EEyx, color='red')
plt.plot(bxy, EExy - EEyx, color='green')
plt.savefig("byx01.png")
plt.show()
'''

'''
nt  = 500

bxy = 0.2e0
byz = bxy
bxz = np.linspace(0.e0, 0.4e0, n)
EExz = np.zeros(n)
EEzx = np.zeros(n)

for j in range(n):
  for i in range(trial):
    x, y, z = discrete_logistic_map(nt, x0, y0, z0, bxy, byz, bxz[j])
    EExz[j] += EE(x, z, p=5)
    EEzx[j] += EE(z, x, p=5)
  EExz[j] /= trial
  EEzx[j] /= trial
  print(EExz[j], EEzx[j])

plt.plot(bxz, EExz, color='blue')
plt.plot(bxz, EEzx, color='red')
plt.savefig("bxy02.png")
plt.show()
'''
