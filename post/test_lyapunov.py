import numpy as np
from scipy.signal import welch
from nolitsa import lyapunov
import matplotlib.pyplot as plt
from mod.mod_ds import Takens_embedding, search_tau
from mod.mod_ds_test import lorenz
from mod.mod_FFT import PSD


nt = 1000
dt = 0.01e0
t  = np.linspace(0.e0, nt * dt, nt+1)
x, y, z = lorenz(nt=nt, dt=dt)
tau = search_tau(x)

lyap_exp = lyapunov.mle_embed(x, dim=[6], tau=tau, maxt=500, window=tau)[0]
plt.plot(lyap_exp)
plt.show()


#freq, psd = PSD(x, t)
freq, psd = welch(x, fs=1.e0 / dt, nperseg=512)

plt.plot(freq, psd)
plt.xscale('log')
plt.yscale('log')
plt.show()

