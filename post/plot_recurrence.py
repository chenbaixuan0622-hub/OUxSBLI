import numpy as np
import os
from scipy.signal import welch
from nolitsa import lyapunov
import matplotlib.pyplot as plt
from mod.mod_recurrence import recurrence_plot
from mod.mod_ds import search_tau
from mod.mod_POD import make_data, make_grid
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

nt = 2000
t  = shock5[nt:,0] * 1e-3
shock5 = np.float64(shock5[nt:,1] - np.mean(shock5[nt:,1])) * 1e-3
shock6 = np.float64(shock6[nt:,1] - np.mean(shock6[nt:,1])) * 1e-3
shock7 = np.float64(shock7[nt:,1] - np.mean(shock7[nt:,1])) * 1e-3
sepa   = np.float64(sepa[nt:,1]   - np.mean(sepa[nt:,1]))


fs = 1.e0 / (-t[0] + t[1])
f, psd = welch(shock6, fs, nperseg=100)
plt.loglog(f, psd)
plt.xlabel('frequency')
plt.ylabel('PSD')
plt.show()


plt.imshow(recurrence_plot(sepa), cmap='grey', origin='lower')
plt.show()


tau = search_tau(sepa)
lyap_exp = lyapunov.mle_embed(sepa[::tau//2], dim=[10], tau=tau, window=tau)[0]
plt.plot(lyap_exp)
plt.show()

