import numpy as np
from numpy.linalg import inv, eig, pinv
import os
from scipy.linalg import svd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_POD import plot_reconstruction_DMD


# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12

# parameter
Q_dir = "../../z6mm"
Lx1   = 28.e-3
Lx2   = 40.e-3
Ly2   = 8.e-3
endT  = 0.1e-3


def DMD(x, y, D, r, t, Q_dir):
  # D[space, time]
  X = D[:,:-1]
  Y = D[:,1:]
  U2, Sig2, Vh2 = svd(X, False)
  U   = U2[:,:r]
  Sig = np.diag(Sig2)[:r,:r]
  V   = Vh2.conj().T[:,:r]
 
  plt.plot(range(1,11), Sig2[:10].real / np.sum(Sig2[:].real) * 100, 'o-')
  plt.xlabel('Mode Index')
  plt.ylabel('Energy (%)')
  save_name = "DMD_energy_contribution.png"
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()

  # build A tilda
  Atil  = np.dot(np.dot(np.dot(U.conj().T, Y), V), inv(Sig))
  mu, W = eig(Atil)
  
  def circle():
    x, y = [], []
    for _x in np.linspace(-180, 180, 360):
      x.append(np.sin(np.radians(_x)))
      y.append(np.cos(np.radians(_x)))
    return x, y

  c_x, c_y = circle()

  plt.plot(c_x, c_y, c='k', linestyle='dashed')
  for i in range(r):
    if i == 0:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label="1st")
    elif i == 1:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label="2nd")
    elif i == 2:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label="3rd")
    elif i >= 3:
      plt.scatter(np.real(mu[i]), np.imag(mu[i]), label=f'{r+1}th')
  plt.xlabel(r"$\it{Re}\,\mu$")
  plt.ylabel(r"$\it{Im}\,\mu$")
  plt.gca().set_aspect('equal', adjustable='box')
  save_name = "DMD_circle.png"
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()

  # build DMD mode
  Phi = np.dot(np.dot(np.dot(Y, V), inv(Sig)), W)
    
  nx = len(x)
  ny = len(y)
  x  = x * 1e3
  y  = y * 1e3
  xmin = x[0]
  xmax = x[-1]
  ymin = y[0]
  ymax = y[-1]
  x, y = np.meshgrid(x, y)
  fig, ax = plt.subplots(r//3, 3, figsize=(10,6))
  for i in range(r):
    DMD_mode = np.real(Phi[:,i]).reshape([ny, nx]) 
    ax[i//3,i%3].set_xlim(xmin, xmax)
    ax[i//3,i%3].set_ylim(ymin, ymax)
    ax[i//3,i%3].set_aspect('equal', adjustable='box')
    ax[i//3,i%3].set_title(f'DMD MODE {i+1}', y=-0.5)
    im = ax[i//3,i%3].contourf(x, y, DMD_mode, levels=50, cmap='jet', extend='both')
    divider = make_axes_locatable(ax[i//3,i%3])
    cax = divider.append_axes('right', '5%', pad='3%')
    fig.colorbar(im, cax=cax, extendrect=True)
    fig.tight_layout()

  save_name = "DMD_Mode.png" 
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()
  
  # compute time evolution
  b   = np.dot(pinv(Phi), X[:,0])
  Psi = np.zeros([r, len(t)], dtype=np.complex64)
  dt  = -t[0] + t[1]
  for i, _t in enumerate(t):
    Psi[:,i] = np.multiply(np.power(mu, _t / dt), b)
  return Atil, Phi, Psi


def make_data(Lx1, Lx2, Ly2, endT, Q_dir):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Nt = len(Q_files)
  t  = np.linspace(0.e0, endT, Nt)

  stridex = 4
  stridey = 8

  first_path = os.path.join(Q_dir, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  for j in range(Ny):
    if Ly2 < Y[j]:
      ny2 = j
      break
  indicesx = np.arange(nx1, nx2, stridex)
  indicesy = np.arange(0,   ny2, stridey)
  nx = len(indicesx)
  ny = len(indicesy)
  x  = X[indicesx]
  y  = Y[indicesy]

  indicesx, indicesy = np.meshgrid(indicesx, indicesy)

  print("nx = ", nx, " ny = ", ny)
  
  D = np.zeros((nx*ny, Nt), dtype=np.float32)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = U[0,indicesy,indicesx]
    D[:,itr] = u.flatten()
    itr += 1

  rank = 6
  Atil, Phi, Psi = DMD(x, y, D, rank, t, Q_dir)
  D2 = np.dot(Phi, Psi)
  name = "DMD.gif"
  
  x, y = np.meshgrid(x,y)
  for i in range(Nt):
    u = np.real(D[:,i]).reshape([ny,nx])
    plt.contourf(x,y,u,levels=50,cmap='jet')
    plt.show()
  #plot_reconstruction_DMD(Q_dir, x, y, D, np.real(D2), name, interval=100)

make_data(Lx1, Lx2, Ly2, endT, Q_dir)

