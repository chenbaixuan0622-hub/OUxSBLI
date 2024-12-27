import numpy as np
import cv2
import os
from scipy.fft import ifft
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, extract_number


# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12

# parameters
#Q_dir = "../../z6mm"
Q_dir = "../3D_solver/TBL/data"#"../../z6mm"
Lx1   = 28.e-3
Lx2   = 38.e-3
Ly2   = 8.e-3
endT  = 0.1e-3

# parameters for SPOD
num_mode = 4
num_freq = 3

def plot_timewise_spod_result(Q_dir, x, y, modes, eigenvalues, freq, t, num_freq, num_mode):
  #        modes[Nf, Nx, Nm=Nb]
  #  eigenvalues[Nf, Nb]
  # eigenvectors[Nf, Nb, Nb]
  nf = 12

  print(eigenvalues[nf,:])

  save_name = "SPOD_energy_contribution_timewise.png"
  fig, ax = plt.subplots(num_freq, figsize=(8, 8))
  ax[0].plot(abs(eigenvalues[nf,:10]) / np.sum(abs(eigenvalues[nf,:])) * 100, 'o-')
  ax[0].set_title(f'Frequency {freq[nf]:.1f}')
  ax[0].set_xlabel('Mode Index')
  ax[0].set_ylabel('Energy (%)')

  fig.tight_layout()
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()

  nx   = len(x)
  ny   = len(y)
  x    = x * 1e3
  y    = y * 1e3
  x, y = np.meshgrid(x, y)
  fig, ax = plt.subplots(num_mode, figsize=(18, 18))
  for i in range(num_mode):
    mode = ifft(modes[nf,:,i])
    mode = mode.reshape([ny,nx])
    ax[i].set_xlabel('x [mm]')
    ax[i].set_ylabel('y [mm]')
    ax[i].set_title(f'Frequency {freq[nf]:.1f} [Hz] Mode {i+1}')
    ax[i].set_aspect('equal', adjustable='box')
    im = ax[i].contourf(x, y, mode.real, levels=50, cmap='jet')
    divider = make_axes_locatable(ax[i])
    cax = divider.append_axes('right', '5%', pad='3%')
    fig.colorbar(im, cax=cax)
    fig.tight_layout()

  save_name = "SPOD_Mode_timewise.png"
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()


Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
Q_files.sort(key=extract_number)
Nt = len(Q_files)

# grid info
first_path = os.path.join(Q_dir, Q_files[0])
Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
nx1 = int(Lx1 / X[-1] * Nx)
nx2 = int(Lx2 / X[-1] * Nx)
for j in range(Ny):
  if Ly2 < Y[j]:
    ny2 = j
    break
stridex  = 4
stridey  = 8
indicesx = np.arange(nx1, nx2, stridex)
indicesy = np.arange(0,   ny2, stridey)
nx = len(indicesx)
ny = len(indicesy)
x  = X[indicesx]
y  = Y[indicesy]

indicesx, indicesy = np.meshgrid(indicesx, indicesy)
t  = np.linspace(0.e0, endT, Nt)

save_path    = os.path.join(Q_dir, "SPOD_eigenvalues.npy")
eigenvalues  = np.load(save_path)
save_path    = os.path.join(Q_dir, "SPOD_eigenvectors.npy")
eigenvectors = np.load(save_path)
save_path    = os.path.join(Q_dir, "SPOD_modes.npy")
modes        = np.load(save_path)
save_path    = os.path.join(Q_dir, "SPOD_freq.npy")
freq         = np.load(save_path)

plot_timewise_spod_result(Q_dir, x, y, modes, eigenvalues, freq, t, num_freq, num_mode)

