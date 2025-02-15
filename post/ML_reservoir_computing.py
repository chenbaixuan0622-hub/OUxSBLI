import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mod_AI.esn import ESN
from mod_AI.utils import train_ESN
from mod.mod_ds_test import gen_data_norm
from mod.mod_plot import set_Params


set_Params()
np.random.seed(99)

# Generate data

N_washout = 100
N_train   = 3000
N_tstart  = 3000
N_test    = 500

N_x = 2
N_dim   = 3 * N_x # dimension of inputs (and outputs)

A = np.array([[0,1],[0,0]])

# Lorenz
X = gen_data_norm(N_x, A, N_tstart+N_test, sigma_dyn=1e-4, sigma_obs=1e-4)

gtcolors   = ['#C7A085', '#5D74A2']
plotcolors = ['#0E419C', '#CC1A1A']


fig = plt.figure(constrained_layout=True, figsize=(6, 2*N_dim))
axs = fig.subplots(N_dim, 1)
for i in range(N_dim):
  axs[i].plot(X.T[i], '-', label="close-loop", linewidth=2.0, color='grey', alpha=1, ms=10)
  axs[i].axis('off')
plt.show()


train_ESN(N_washout, N_train, N_tstart, N_test, ESN, X)

