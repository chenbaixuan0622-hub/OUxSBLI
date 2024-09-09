import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
# DMD module
from mod.mod_DMD import DMD

x  = np.linspace(-10.e0, 10.e0, 100)
t  = np.linspace(0.e0, 6.e0 * np.pi, 80)
Xm, Tm = np.meshgrid(x, t)

f1 = np.multiply(20.e0 - 0.2e0 * np.power(Xm, 2), np.exp((2.3j)*Tm))
f2 = np.multiply(Xm, np.exp(0.6j * Tm))
f3 = np.multiply(5.e0 * np.multiply(1.e0 / np.cosh(Xm/2), np.tanh(Xm/2)), 2.e0 * np.exp((0.1+2.8j)*Tm))


# transpose for np.meshgrid
D = (f1 + f2 + f3).T

Psi, D2 = DMD(t, D, 3)

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
ax[-1].set_title('f1_raw')
ax[-1].plot(t, f1.T.real[indx,:], color='blue', label='Real', marker='x')
ax[-1].plot(t, f1.T.imag[indx,:], color='green', label='Complex', marker='x')
ax_count += 1

ax.append(fig.add_subplot((max_plot_num), 1, ax_count))
ax[-1].set_xlabel('t')
ax[-1].set_ylabel('y')
ax[-1].set_title('f2_raw')
ax[-1].plot(t, f2.T.real[indx,:], color='blue', label='Real', marker='x')
ax[-1].plot(t, f2.T.imag[indx,:], color='green', label='Complex', marker='x')
ax_count += 1

ax.append(fig.add_subplot((max_plot_num), 1, ax_count))
ax[-1].set_xlabel('t')
ax[-1].set_ylabel('y')
ax[-1].set_title('f3_raw')
ax[-1].plot(t, f3.T.real[indx,:], color='blue', label='Real', marker='x')
ax[-1].plot(t, f3.T.imag[indx,:], color='green', label='Complex', marker='x')
ax_count += 1

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

