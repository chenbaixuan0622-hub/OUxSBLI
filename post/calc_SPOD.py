import numpy as np
import cv2
import os
from scipy.fft import ifft
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_POD import spanwise_spod, timewise_spod, reconstruct_spanwise_spod, reconstruct_timewise_spod 
from mod.mod_POD import plot_energy_contribution, plot_reconstruction


# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12

# parameter
Q_dir = "../../z6mm"
Lx1   = 28.e-3
Lx2   = 38.e-3
Ly2   = 8.e-3
endT  = 0.1e-3


def plot_spanwise_spod_result(Q_dir, U, Sigma, freq, t, num_freq, num_modes):
  #     U[Nf, Nt, Nm=Nt]
  # Sigma[Nf, Nt]
  name = "energy_contribution_spanwise_spod.png"
  plot_energy_contribution(Q_dir, Sigma, freq, num_freq, name)

  fig, ax = plt.subplots(num_modes, num_freq, figsize=(12, 12))
  for j in range(num_freq):
    for i in range(num_modes):
      wl = 1.e3 / freq[j+1]
      ax[i,j].plot(t*1e3, abs(U[j+1,:,i]))
      ax[i,j].set_title(f'Wavelength {wl:.1f} [mm] Mode {i+1}')
      ax[i,j].set_xlabel('Time [ms]')
      ax[i,j].set_ylabel('Amplitude')

  fig.tight_layout()
  save_name = "time_evolution_mode_spanwise_spod.png"
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()


def plot_timewise_spod_result(Q_dir, x, y, modes, eigenvalues, freq, t, num_freq, num_modes):
  #        modes[Nf, Nx, Nm=Nb]
  #  eigenvalues[Nf, Nb]
  # eigenvectors[Nf, Nb, Nb]
  name = "SPOD_energy_contribution_timewise.png"
  plot_energy_contribution(Q_dir, eigenvalues, freq, num_freq, name)

  nx   = len(x)
  nz   = len(y)
  x    = x * 1e3
  y    = y * 1e3
  x, y = np.meshgrid(x, y)
  fig, ax = plt.subplots(num_modes, num_freq, figsize=(18, 18))
  for j in range(num_freq):
    for i in range(num_modes):
      mode = ifft(modes[j+1,:,i])
      mode = mode.reshape([nz,nx])
      ax[i,j].set_xlabel('x [mm]')
      ax[i,j].set_ylabel('z [mm]')
      ax[i,j].set_title(f'Frequency {freq[j+1]:.1f} [Hz] Mode {i+1}')
      ax[i,j].set_aspect('equal', adjustable='box')
      im = ax[i,j].contourf(x, y, mode.real, levels=50, cmap='jet')
      divider = make_axes_locatable(ax[i,j])
      cax = divider.append_axes('right', '5%', pad='3%')
      fig.colorbar(im, cax=cax)
      fig.tight_layout()

  save_name = "SPOD_Mode_timewise.png"
  save_path = os.path.join(Q_dir, save_name)
  plt.savefig(save_path)
  plt.close()


def calc_spanwise_spod(nx1, nx2, nx, nz, Nx, Ny, Nz, Nt, x, z, t, Q_dir, Q_files):
  data = np.zeros((Nt,nx,nz), dtype=np.float32)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    U = np.reshape(U[:,0,nx1:nx2], [Nz, nx2-nx1])
    u = cv2.resize(U, (nx, nz))
    data[itr,:,:] = u.T
    itr += 1

  # spanwise SPOD: data[Nt, nx, nz]
  num_divide = 2
  num_mode   = 6
  num_freq   = 3
  name       = "SPOD_spanwise.gif"
  U, Sigma, Vt, freq = spanwise_spod(data, num_divide, dz=-z[0]+z[1])
  SPOD = reconstruct_spanwise_spod(U, Sigma, Vt, num_divide, num_mode)
  plot_spanwise_spod_result(Q_dir, U, Sigma, freq, t, num_freq, num_mode)
  plot_reconstruction(Q_dir, x, z, data, SPOD, name, interval=200)


def calc_timewise_spod(indicesx, indicesy, Nx, Ny, Nz, x, y, t, Q_dir, Q_files):
  nx = len(x)
  ny = len(y)
  Nt = len(t)
  data = np.zeros((nx*ny,Nt), dtype=np.float32)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = U[0,indicesy,indicesx]
    data[:,itr] = u.flatten()
    itr += 1

  # timewise SPOD: data[nx*ny, Nt]
  num_divide = 2
  num_mode   = 6
  num_freq   = 3
  name       = "SPOD_timewise.gif"
  eigenvalues, eigenvectors, modes, freq = timewise_spod(data, num_divide, dt=-t[0]+t[1])
  SPOD = reconstruct_timewise_spod(eigenvalues, eigenvectors, modes, num_divide, num_mode)
  # data[Nt, nx, ny]
  Nt   = (Nt // num_divide) * num_divide
  data = data[:,:Nt].T
  SPOD = SPOD[:,:Nt].T
  data = data.reshape([Nt,ny,nx])
  SPOD = SPOD.reshape([Nt,ny,nx])
  data = data.transpose((0,2,1))
  SPOD = SPOD.transpose((0,2,1))
  plot_timewise_spod_result(Q_dir, x, y, modes, eigenvalues, freq, t, num_freq=3, num_modes=2)
  plot_reconstruction(Q_dir, x, y, data, SPOD, name, interval=200)


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

#calc_spanwise_spod(nx1, nx2, nx, nz, Nx, Ny, Nz, Nt, x, z, t, Q_dir, Q_files)
calc_timewise_spod(indicesx, indicesy, Nx, Ny, Nz, x, y, t, Q_dir, Q_files)

