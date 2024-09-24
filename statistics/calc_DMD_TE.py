import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy
# read module
from mod.mod_read import getGrid
# DMD module
from mod.mod_DMD import SVD, DMD

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240819_KEEP6thVisc4th"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

info_directory = os.path.join(Q_directory, "info")

t0 = 0.385384e-3
t1 = 2.e0 * t0
t  = np.linspace(t0, t1, num_files)

tke_path = os.path.join(Q_directory, "tke_log.npy")
# TKE[space, time]
TKE = np.load(tke_path)

nx = len(TKE[:,0])
TE12     = np.zeros((nx), dtype=np.float32)
TE21     = np.zeros((nx), dtype=np.float32)
TE12_DMD = np.zeros((nx), dtype=np.float32)
TE21_DMD = np.zeros((nx), dtype=np.float32)

x = np.linspace(0, 2.e-3, nx)

def NTE(x, y, history_len):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
  xs = bin_series(x, b=10)
  ys = bin_series(y, b=10)
  TEx_y = transfer_entropy(xs[0], ys[0], k=history_len)
  TEy_x = transfer_entropy(ys[0], xs[0], k=history_len)
  # H(Yt|Yt-1)
  Yt   = bin_series(y[1:],  b=10)
  Yt_1 = bin_series(y[:-1], b=10)
  Xt   = bin_series(x[1:],  b=10)
  Xt_1 = bin_series(x[:-1], b=10)
  Hy = conditional_entropy(Yt_1[0], Yt[0]) # H(Yt|Yt-1)
  Hx = conditional_entropy(Xt_1[0], Xt[0]) # H(Xt|Xt-1)
  # E[TExd_y] = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xdt-1:t-L)
  x_shuffled = np.random.permutation(x)
  y_shuffled = np.random.permutation(y)
  xd = bin_series(x_shuffled, b=10)
  yd = bin_series(y_shuffled, b=10)
  Ex_y = transfer_entropy(xd[0], ys[0], k=history_len)
  Ey_x = transfer_entropy(yd[0], xs[0], k=history_len)
  # calc normalized TE
  NTEx_y = (TEx_y - Ex_y) / Hy
  NTEy_x = (TEy_x - Ey_x) / Hx
  return NTEx_y, NTEy_x

for i in range(len(TKE[:,0])):
  TE12[i], TE21[i] = NTE(TKE[0,:], TKE[i,:], 5)

U2, Sig2, Vh2 = SVD(t, TKE)
plt.scatter(range(0, len(Sig2)), Sig2)
plt.show()

rank = int(input("Enter the rank: "))
mu, Phi, Psi, TKE2 = DMD(rank, t, TKE, U2, Sig2, Vh2)

Re = np.linspace(-1.e0, 1.e0, 100)
plt.axis("equal")
plt.xlim([-1.2, 1.2])
plt.ylim([-1.2, 1.2])
plt.plot(Re,  np.sqrt(1.e0 - Re**2), '-', color="black")
plt.plot(Re, -np.sqrt(1.e0 - Re**2), '-', color="black")
for i in range(len(mu)):
  plt.scatter(np.real(mu[i]), np.imag(mu[i]))
plt.show()

fig, ax = plt.subplots()
for i in range(len(mu)):
  ax.plot(x, np.real(Phi[:,i]), label="mode" + str(i+1))
ax.set_xlabel("space")
ax.legend()
plt.show()

#print(np.shape(TKE2))
plt.plot(t, TKE[0,:])
plt.plot(t, TKE2[0,:].real)
plt.show()

'''
for i in range(nx):
  TE12_DMD[i], TE21_DMD[i] = NTE(TKE2[0,:], TKE2[i,:], 5)
'''
'''
file_path = os.path.join(info_directory, "TE_tke_DMD.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  print('# x     TE       TE       TE1      TE1       TE2       TE2      TE3      TE3      TE4       TE4       TE5       TE5', file=fo)
  for i in range(nx):
    print(f'{x[i]-x[0]:.3e}', f'{TE12[i]:.3e}', f'{TE12[i]:.3e}', \
          f'{TE12_DMD[i]:.3e}', f'{TE21_DMD[i]:.3e}', file=fo)
'''
'''
ax = []
plt.clf()
ax_count = 1
max_plot_num = 10
indx = 10
fig = plt.figure(figsize=(15, 40))
fig.subplots_adjust(hspace=1)

ax.append(fig.add_subplot((max_plot_num), 1, ax_count))
ax[-1].set_xlabel('t')
ax[-1].set_ylabel('y')
ax[-1].set_title('fusion')
ax[-1].plot(t, D.real[indx,:], color='blue', label='Real', marker='x')
ax[-1].plot(t, D.imag[indx,:], color='green', label='Complex', marker='x')
ax_count += 1

ax.append(fig.add_subplot((max_plot_num), 1, ax_count))
ax[-1].set_xlabel('t')
ax[-1].set_ylabel('y')
ax[-1].set_title('DMD (m1)')
ax[-1].plot(t, Psi.real[0,:], color='blue', label='Original', marker='x')
ax[-1].plot(t, Psi.imag[0,:], color='green', label='Complex', marker='x')
ax_count += 1

ax.append(fig.add_subplot((max_plot_num), 1, ax_count))
ax[-1].set_xlabel('t')
ax[-1].set_ylabel('y')
ax[-1].set_title('DMD (m2)')
ax[-1].plot(t, Psi.real[1,:], color='blue', label='Original', marker='x')
ax[-1].plot(t, Psi.imag[1,:], color='green', label='Complex', marker='x')
ax_count += 1

ax.append(fig.add_subplot((max_plot_num), 1, ax_count))
ax[-1].set_xlabel('t')
ax[-1].set_ylabel('y')
ax[-1].set_title('DMD (m3)')
ax[-1].plot(t, Psi.real[2,:], color='blue', label='Original', marker='x')
ax[-1].plot(t, Psi.imag[2,:], color='green', label='Complex', marker='x')
ax_count += 1

ax.append(fig.add_subplot((max_plot_num), 1, ax_count))
ax[-1].set_xlabel('t')
ax[-1].set_ylabel('y')
ax[-1].set_title('fusion')
ax[-1].plot(t, D2.real[indx,:], color='blue', label='Real', marker='x')
ax[-1].plot(t, D2.imag[indx,:], color='green', label='Complex', marker='x')
ax_count += 1

plt.savefig( 'output.png' )
'''

