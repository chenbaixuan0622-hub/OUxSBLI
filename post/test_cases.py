import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_info import EE
from mod.mod_ds import case1, case2, case3, case4, coupling_system, discrete_logistic_map

set_Params()


w   = 0.1e0 

'''
# case 1
x, y = case1(nt)
EE1xy = EE(x, y, p=5)
EE1yx = EE(y, x, p=5)

print("case 1 x->y:", EE1xy, "y->x", EE1yx)

# case 2
x, y = case2(nt, w)
EE2xy = EE(x, y, p=5)
EE2yx = EE(y, x, p=5)

print("case 2 x->y:", EE2xy, "y->x", EE2yx)

# case 3
x, y = case3(nt, w)
EE3xy = EE(x, y, p=5)
EE3yx = EE(y, x, p=5)

print("case 3 x->y:", EE3xy, "y->x", EE3yx)

# case 4
x, y, z = case4(nt, w)
EE4xy = EE(x, y, p=7)
EE4yx = EE(y, x, p=7)
EE4yz = EE(y, z, p=7)
EE4zy = EE(z, y, p=7)
EE4zx = EE(z, x, p=7)
EE4xz = EE(x, z, p=7)

print("case 4 x->y:", EE4xy, "y->x", EE4yx)
print("case 4 y->z:", EE4yz, "z->y", EE4zy)
print("case 4 z->x:", EE4zx, "x->z", EE4xz)
'''

nt  = 250

# coupling system
n   = 41
byx = 0.e0
bxy = np.linspace(0.e0, 0.4e0, n)

x0  = 0.5e0
y0  = 0.5e0
z0  = 0.5e0

EExy = np.zeros(n)
EEyx = np.zeros(n)

trial = 100

for j in range(n):
  for i in range(trial):
    x, y = coupling_system(nt, x0, y0, bxy[j], byx)
    EExy[j] += EE(x, y, p=5)
    EEyx[j] += EE(y, x, p=5)
  EExy[j] /= trial
  EEyx[j] /= trial
  print(EExy[j], EEyx[j])

plt.plot(bxy, EExy, color='blue')
plt.plot(bxy, EEyx, color='red')
plt.savefig("byx00.png")
plt.show()


byx = 0.1e0
EExy = 0.e0
EEyx = 0.e0

for j in range(n):
  for i in range(trial):
    x, y = coupling_system(nt, x0, y0, bxy[j], byx)
    EExy[j] += EE(x, y, p=5)
    EEyx[j] += EE(y, x, p=5)
  EExy[j] /= trial
  EEyx[j] /= trial
  print(EExy[j], EEyx[j])

plt.plot(bxy, EExy, color='blue')
plt.plot(bxy, EEyx, color='red')
plt.savefig("byx01.png")
plt.show()


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

