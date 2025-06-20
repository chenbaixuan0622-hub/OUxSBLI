import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mod.mod_plot import set_Params
#from mod.mod_info import EE
#from mod.mod_recurrence import recurrence_plot, excessive_recurrence_test
#from mod.mod_ds import lorenz, search_tau, Takens_embedding, local_constant_pred
from mod.mod_ds_test import lorenz


set_Params()


dynamics = lorenz(s=10.e0, r=28.e0, b=8.e0/3.e0)
xf, yf, zf = dynamics.GaussRK2(nt=100000, dt= 0.001e0)
xb, yb, zb = dynamics.GaussRK2(nt=2500,   dt=-0.001e0, x0=xf[-1], y0=yf[-1], z0=zf[-1])

ax = plt.axes(projection='3d')
ax.view_init(elev=30, azim=-60)
ax.plot3D(xf, yf, zf, 'blue', lw=0.5)
ax.plot3D(xb, yb, zb, 'red', lw=0.5)
ax.set_xlabel('x', fontsize=16)
ax.set_ylabel('y', fontsize=16)
ax.set_zlabel('z', fontsize=16)
plt.show()
plt.close()


# causal patterns
# data1 data2 dim1 dim2 
#     x     y    0    1
#     x     z    0    2
#     y     z    0    2

'''
fig, ax = plt.subplots(1,3, figsize=(12,4), tight_layout=True)
ax[0].plot(t, xs)
ax[0].set_xlabel('$\it{t}$', fontsize=16)
ax[0].set_ylabel('$\it{x}$', fontsize=16)
ax[1].plot(t, ys)
ax[1].set_xlabel('$\it{t}$', fontsize=16)
ax[1].set_ylabel('$\it{y}$', fontsize=16)
ax[2].plot(t, zs)
ax[2].set_xlabel('$\it{t}$', fontsize=16)
ax[2].set_ylabel('$\it{z}$', fontsize=16)
plt.show()
plt.close()


taux, mix = search_tau(xs, bin=5, tau_max=nt-1)
tauy, miy = search_tau(ys, bin=5, tau_max=nt-1)
tauz, miz = search_tau(zs, bin=5, tau_max=nt-1)

print("taux ", taux)
print("tauy ", tauy)
print("tauz ", tauz)

fig, ax = plt.subplots(1,3, figsize=(12,4), tight_layout=True)
ax[0].plot(t[:len(t)//2], mix[:len(t)//2])
ax[0].set_xlabel(r'$\it{\tau}$', fontsize=16)
ax[0].set_ylabel('$\it{MI}$', fontsize=16)
ax[1].plot(t[:len(t)//2], miy[:len(t)//2])
ax[1].set_xlabel(r'$\it{\tau}$', fontsize=16)
ax[1].set_ylabel('$\it{MI}$', fontsize=16)
ax[2].plot(t[:len(t)//2], miz[:len(t)//2])
ax[2].set_xlabel(r'$\it{\tau}$', fontsize=16)
ax[2].set_ylabel('$\it{MI}$', fontsize=16)
plt.show()
plt.close()


xe = Takens_embedding(xs, taux, dim=3)
ye = Takens_embedding(ys, tauy, dim=3)
ze = Takens_embedding(zs, tauz, dim=3)


ax = plt.axes(projection='3d')
ax.view_init(elev=30, azim=-60)
ax.plot3D(xe[:,0], xe[:,1], xe[:,2], 'blue', lw=0.5)
ax.set_xlabel('$\it{x_{1}}$', fontsize=16)
ax.set_ylabel('$\it{x_{2}}$', fontsize=16)
ax.set_zlabel('$\it{x_{3}}$', fontsize=16)
fig.tight_layout()
plt.show()
plt.close()


ax = plt.axes(projection='3d')
ax.view_init(elev=30, azim=-60)
ax.plot3D(ye[:,0], ye[:,1], ye[:,2], 'blue', lw=0.5)
ax.set_xlabel('$\it{y_{1}}$', fontsize=16)
ax.set_ylabel('$\it{y_{2}}$', fontsize=16)
ax.set_zlabel('$\it{y_{3}}$', fontsize=16)
plt.show()
plt.close()


ax = plt.axes(projection='3d')
ax.view_init(elev=30, azim=-60)
ax.plot3D(ze[:,0], ze[:,1], ze[:,2], 'blue', lw=0.5)
ax.set_xlabel('$\it{z_{1}}$', fontsize=16)
ax.set_ylabel('$\it{z_{2}}$', fontsize=16)
ax.set_zlabel('$\it{z_{3}}$', fontsize=16)
plt.show()
plt.close()


fig, ax = plt.subplots(1,3, figsize=(12,4), tight_layout=True)
extent  = (0, t[-1], 0, t[-1])
ax[0].imshow(recurrence_plot(xs), extent=extent, cmap='Greys', origin='lower')
ax[0].set_xlabel('$\it{t}$')
ax[0].set_ylabel('$\it{t}$')
ax[1].imshow(recurrence_plot(ys), extent=extent, cmap='Greys', origin='lower')
ax[1].set_xlabel('$\it{t}$')
ax[1].set_ylabel('$\it{t}$')
ax[2].imshow(recurrence_plot(zs), extent=extent, cmap='Greys', origin='lower')
ax[2].set_xlabel('$\it{t}$')
ax[2].set_ylabel('$\it{t}$')
plt.show()
plt.close()
'''
'''
# local constatnt prediction
m = 20
taux, _ = search_tau(xs, bin=5, tau_max=len(xs)-1)
indices = np.arange(0, len(xs), taux)
T = t[indices]
x = xs[indices]

xp, xa = local_constant_pred(T, x, m, dim=7, k=3)

plt.plot(xa, color='blue')
plt.plot(xp, color='red')
plt.xlabel(r'$\it{\tau}$', fontsize=16)
plt.ylabel('$\it{x}$', fontsize=16)
plt.xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
plt.show()
'''
'''
data = [xs, ys, zs, noise]

Map = np.zeros((4,4), dtype=np.float32)

for j in range(4):
  for i in range(4):
    if i == j:
      Map[j,i] = np.inf
    else:
      Map[j,i] = EE(data[j], data[i], p=3)

plt.imshow(Map, cmap='jet', extent=None, origin='lower')
plt.xlabel("Effect", fontsize=24)
plt.ylabel("Cause", fontsize=24)
plt.xticks([0,1,2,3], ['$\it{x}$', '$\it{y}$', '$\it{z}$', '$\it{\epsilon}$'], fontsize=20)
plt.yticks([0,1,2,3], ['$\it{x}$', '$\it{y}$', '$\it{z}$', '$\it{\epsilon}$'], fontsize=20)
plt.show()
'''
