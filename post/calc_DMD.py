import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_POD import make_data


set_Params()

# parameter
#Q_dir = "../../z6mm"
Q_dir = "../3D_solver/TBL/data"#"../../z6mm"
Lx1   = 28.e-3
Lx2   = 40.e-3
Ly2   = 8.e-3
endT  = 0.1e-3
rank  = 20


def DMD(D, rank, t):
  X = D[:,:-1]
  Y = D[:,1:]
  U, S, Vh = np.linalg.svd(X, full_matrices=False)
  Ur  = U[:,:rank]
  Sr  = np.diag(S[:rank])
  Vhr = Vh[:rank,:]

  # obtain Atilda by computing the pseudo-inverse of X
  Atilda = np.linalg.solve(Sr, (Ur.T @ Y @ Vhr.T))

  # the spectral decomposition of Atilda
  mu, W = np.linalg.eig(Atilda)
  mu    = np.diag(mu)

  # reconstructed the high-dimensional DMD modes
  Phi = Y @ np.linalg.solve(Sr.T, Vhr).T @ W
  x1tilda = Sr @ Vhr[:,0]
  b     = np.linalg.solve(W @ mu, x1tilda)
  dt    = -t[0] + t[1]
  omega = np.log(np.diag(mu)) / dt

  # reconstruction Xr for all time steps
  Dr = np.zeros_like(D, dtype=np.float32)
  for i in range(len(t)):
    Dr[:,i] = np.real(Phi @ (b * np.exp(omega * t[i])))
  return S, Atilda, Phi, mu, b, omega, Dr


def plot_DMD_mode(x, y, rank, Phi, Q_dir):
  nx = len(x)
  ny = len(y)
  x  = x * 1e3
  y  = y * 1e3
  xmin = x[0]
  xmax = x[-1]
  ymin = y[0]
  ymax = y[-1]
  x, y = np.meshgrid(x, y)
  for i in range(rank):
    fig, ax = plt.subplots(1, figsize=(6,6))
    DMD_mode = np.real(Phi[:,i]).reshape([ny, nx]) 
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f'DMD MODE {i+1}', y=-0.2)
    im = ax.contourf(x, y, DMD_mode, levels=50, cmap='jet', extend='both')
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', '5%', pad='3%')
    #fig.colorbar(im, cax=cax, extendrect=True)
    #cbar = fig.colorbar(im)
    fig.tight_layout()

    save_name = "DMD_Mode" + str(i+1) + ".png" 
    save_path = os.path.join(Q_dir, "DMD_mode", save_name)
    plt.savefig(save_path)
    plt.close()
  

def plot_DMD_circle(mu, rank, Q_dir):
  def circle():
    x, y = [], []
    for _x in np.linspace(-180, 180, 360):
      x.append(np.sin(np.radians(_x)))
      y.append(np.cos(np.radians(_x)))
    return x, y

  c_x, c_y = circle()

  plt.plot(c_x, c_y, c='k', linestyle='dashed')
  for i in range(rank):
    if i == 0:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label="1st")
    elif i == 1:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label="2nd")
    elif i == 2:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label="3rd")
    elif i >= 3:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label=f'{i+1}th')
  plt.xlabel(r"$\it{Re}\,\mu$")
  plt.ylabel(r"$\it{Im}\,\mu$")
  plt.gca().set_aspect('equal', adjustable='box')
  save_name = "DMD_circle.png"
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()


def plot_DMD_reconstruct(x, y, D, D2, Q_dir):
  nx = len(x)
  ny = len(y)
  Nt = len(D[0,:])
  x, y = np.meshgrid(x*1e3, y*1e3)
  fig, ax = plt.subplots(1, 2, figsize=(12,8))
  crange  = np.linspace(0.e0, 500.e0, 50)
  im1 = ax[0].contourf(x, y,  D[:,0].reshape(ny, nx), crange, cmap='jet', extend='both')
  im2 = ax[1].contourf(x, y, D2[:,0].reshape(ny, nx), crange, cmap='jet', extend='both')
  cb  = fig.colorbar(im1, ax=ax, extendrect=True, orientation='horizontal', fraction=0.046, location='top')

  def update(frame):
    for a in ax:
      a.clear()
      a.set_aspect('equal', adjustable='box')
    im1 = ax[0].contourf(x, y,  D[:,frame].reshape(ny, nx), crange, cmap='jet', extend='both')
    im2 = ax[1].contourf(x, y, D2[:,frame].reshape(ny, nx), crange, cmap='jet', extend='both')
    return im1.collections + im2.collections

  ani = FuncAnimation(fig, update, frames=range(Nt), blit=True, interval=200)
  save_path = os.path.join(Q_dir, "DMD.gif")
  ani.save(save_path, writer='Pillow')

# D[space, time]
save_path = os.path.join(Q_dir, "D.npy")
if os.path.isfile(save_path):
  D = np.load(save_path)
else:
  x, y, t, D = make_data(Lx1, Lx2, Ly2, endT, Q_dir)
  np.save(save_path, D)

S, Atilda, Phi, mu, b, omega, Dr = DMD(D, rank, t)

'''
nx = len(x)
ny = len(y)
x, y = np.meshgrid(x*1e3, y*1e3)
fig, ax = plt.subplots(figsize=(8,8))
crange  = np.linspace(0.e0, 500.e0, 50)
im = ax.contourf(x, y, Xr.real.reshape(ny, nx), crange, cmap='jet', extend='both')
cb = fig.colorbar(im, ax=ax, extendrect=True, orientation='horizontal', fraction=0.046, location='top')
plt.show()
'''
plt.plot(range(1,rank+1), S[:rank].real / np.sum(S.real) * 100, 'o-')
plt.xlabel('Mode Index')
plt.ylabel('Energy (%)')
save_name = "DMD_energy_contribution.png"
save_path = os.path.join(Q_dir, save_name)
plt.savefig(save_path)
plt.close()

plot_DMD_mode(x, y, rank, Phi, Q_dir)
plot_DMD_circle(mu, rank, Q_dir)
plot_DMD_reconstruct(x, y, D, Dr, Q_dir)

