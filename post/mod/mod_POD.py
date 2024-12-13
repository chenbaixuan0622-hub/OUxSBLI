import numpy as np
import os
from scipy.linalg import eigh, svd
from scipy.fft import rfft, irfft, rfftfreq
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def snapshot_pod(data):
  # data[space, time]
  C_s  = np.dot(data.T, data)
  eigenvalues, eigenvectors = eigh(C_s)
  del C_s
  idx = eigenvalues.argsort()[::-1]
  eigenvalues  = eigenvalues[idx]
  eigenvectors = eigenvectors[:,idx]

  modes  = np.dot(data, eigenvectors)
  modes /= np.linalg.norm(modes, axis=0)
    
  return eigenvalues, eigenvectors, modes


def calc_time_coef(data, modes):
  time_coef = np.dot(data.T, modes)
  return time_coef


def spanwise_spod(data, num_divide, dz):
  # data[Nt, Nx, Nz]
  Nt = len(data[:,0,0])
  Nx = len(data[0,:,0])
  Nz = len(data[0,0,:])
  nz = Nz // num_divide
  Nb = Nx * num_divide
  data = data[:,:,:nz*num_divide]
  data = data.reshape([Nt,Nb,nz])
  data = rfft(data, axis=-1)
  freq = rfftfreq(nz, dz)
  # data[Nt, Nb, Nf]
  Nf = len(data[0,0,:])

  U     = np.zeros((Nf,Nt,Nt), dtype=np.complex64)
  Sigma = np.zeros((Nf,Nt),    dtype=np.complex64)
  Vt    = np.zeros((Nf,Nb,Nb), dtype=np.complex64)
  for i in range(Nf):
    U[i,:,:], Sigma[i,:], Vt[i,:,:] = svd(data[:,:,i])

  return U, Sigma, Vt, freq


def timewise_spod(data, num_divide, dt):
  # data[Nx, Nt]
  Nt = len(data[0,:])
  Nx = len(data[:,0])
  nt = Nt // num_divide
  Nb = num_divide 
  data = data[:,:nt*num_divide]
  data = data.reshape([Nx,Nb,nt])
  data = rfft(data, axis=-1)
  freq = rfftfreq(nt, dt)
  # data[Nx, Nb, Nf]
  Nf = len(data[0,0,:])

  eigenvalues  = np.zeros((Nf,Nb),    dtype=np.complex64)
  eigenvectors = np.zeros((Nf,Nb,Nb), dtype=np.complex64)
  modes        = np.zeros((Nf,Nx,Nb), dtype=np.complex64)
  time_coef    = np.zeros((Nf,), dtype=np.complex64)
  for i in range(Nf):
    # snapshot POD
    C_s = np.dot(data[:,:,i].T, data[:,:,i])
    eigenvalue, eigenvector = eigh(C_s)
    del C_s
    idx = eigenvalue.argsort()[::-1]
    eigenvalue  = eigenvalue[idx]
    eigenvector = eigenvector[:,idx]

    mode  = np.dot(data[:,:,i], eigenvector)
    mode /= np.linalg.norm(mode, axis=0)
   
    eigenvalues[i,:]    = eigenvalue
    eigenvectors[i,:,:] = eigenvector
    modes[i,:,:]        = mode
  return eigenvalues, eigenvectors, modes, freq


def reconstruct_spanwise_spod(U, Sigma, Vt, num_divide, num_modes):
  #  data[Nt, Nb, Nf]
  #     U[Nf, Nt, Nt]
  # Sigma[Nf, Nt]
  #    Vt[Nf, Nb, Nb]
  Nt = len(Sigma[0,:])
  Nb = len(Vt[0,0,:])
  Nf = len(Sigma[:,0])
  data = np.zeros((Nt,Nb,Nf), dtype=np.complex64)
  for i in range(Nf):
    data[:,:,i] = U[i,:,:num_modes] @ np.diag(Sigma[i,:num_modes]) @ Vt[i,:num_modes,:]
  
  data = irfft(data, axis=-1)
  # data[Nt, Nb, nz]
  nz = len(data[0,0,:])
  Nz = nz * num_divide
  Nx = Nb // num_divide
  data = data.reshape([Nt,Nx,Nz])
  return data


def reconstruct_timewise_spod(eigenvalues, eigenvectors, modes, num_divide, num_modes):
  #         data[Nx, Nt]
  #  eigenvalues[Nf, Nb]
  # eigenvectors[Nf, Nb, Nb]
  #        modes[Nf, Nx, Nb]
  Nx = len(modes[0,:,0])
  Nb = len(eigenvalues[0,:])
  Nf = len(eigenvalues[:,0])
  data = np.zeros((Nx,Nb,Nf), dtype=np.complex64)
  for i in range(Nf):
    data[:,:,i] = np.dot(modes[i,:,:num_modes] * np.sqrt(eigenvalues[i,:num_modes]), eigenvectors[i,:,:num_modes].T)

  data = irfft(data, axis=-1)
  # data[Nx, Nb, nt]
  nt = len(data[0,0,:])
  Nt = nt * num_divide
  data = data.reshape([Nx,Nt])
  return data


def plot_reconstruction(Q_dir, x, z, data, SPOD, name, interval=200):
  # data[Nt, nx, nz]
  Nt   = len(data[:,0,0])
  nx   = len(x)
  nz   = len(z)
  x    = x * 1e3
  z    = z * 1e3
  x, z = np.meshgrid(x, z)

  fig, ax  = plt.subplots(1, 3, figsize=(18, 6))
  crange   = np.linspace(0, 500, 50)
  contour1 = ax[0].contourf(x, z, data[0,:,:].T, crange, cmap='jet', extend='both')
  contour2 = ax[1].contourf(x, z, SPOD[0,:,:].T, crange, cmap='jet', extend='both')
  contour3 = ax[2].contourf(x, z, SPOD[0,:,:].T - data[0,:,:].T, levels=50, cmap='jet')
  cbar1    = fig.colorbar(contour1, ax=ax[0:1], extendrect=True, \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar1.set_label("u [m/s]")
  cbar1.ax.xaxis.set_ticks_position('top')
  cbar1.ax.xaxis.set_label_position('top')
  cbar2    = fig.colorbar(contour3, ax=ax[2], extendrect=True, \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar2.set_label("SPOD - u [m/s]")
  cbar2.ax.xaxis.set_ticks_position('top')
  cbar2.ax.xaxis.set_label_position('top')

  def update(frame):
    for a in ax:
      a.clear()
      a.set_aspect('equal', adjustable='box')
    contour1 = ax[0].contourf(x, z, data[frame,:,:].T, crange, cmap="jet", extend='both')
    contour2 = ax[1].contourf(x, z, SPOD[frame,:,:].T, crange, cmap='jet', extend='both')
    contour3 = ax[2].contourf(x, z, SPOD[frame,:,:].T - data[frame,:,:].T, levels=50, cmap='jet')
    return contour1.collections + contour2.collections + contour3.collections

  ani = FuncAnimation(fig, update, frames=Nt, interval=interval, blit=False)

  save_path = os.path.join(Q_dir, name)
  ani.save(save_path, writer='Pillow')


def plot_reconstruction_DMD(Q_dir, x, y, D, D2, name, interval=200):
  # D, D2 [space, time]
  nx   = len(x)
  ny   = len(y)
  Nt   = len(D) // (nx * ny)
  x    = x * 1e3
  y    = y * 1e3
  x, y = np.meshgrid(x, y)

  fig, ax  = plt.subplots(1, 3, figsize=(18, 6))
  crange   = np.linspace(0, 500, 50)
  contour1 = ax[0].contourf(x, y,  D[:,0].reshape([ny,nx]), crange, cmap='jet', extend='both')
  contour2 = ax[1].contourf(x, y, D2[:,0].reshape([ny,nx]), crange, cmap='jet', extend='both')
  contour3 = ax[2].contourf(x, y, (D2[:,0] - D[:,0]).reshape([ny,nx]), levels=50, cmap='jet')
  cbar1    = fig.colorbar(contour1, ax=ax[0:1], extendrect=True, \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar1.set_label("u [m/s]")
  cbar1.ax.xaxis.set_ticks_position('top')
  cbar1.ax.xaxis.set_label_position('top')
  cbar2    = fig.colorbar(contour3, ax=ax[2], extendrect=True, \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar2.set_label("DMD - u [m/s]")
  cbar2.ax.xaxis.set_ticks_position('top')
  cbar2.ax.xaxis.set_label_position('top')

  def update(frame):
    for a in ax:
      a.clear()
      a.set_aspect('equal', adjustable='box')
    contour1 = ax[0].contourf(x, y,  D[:,frame].reshape([ny,nx]), crange, cmap="jet", extend='both')
    contour2 = ax[1].contourf(x, y, D2[:,frame].reshape([ny,nx]), crange, cmap='jet', extend='both')
    contour3 = ax[2].contourf(x, y, (D2[:,frame] - D[:,frame]).reshape([ny,nx]), levels=50, cmap='jet')
    return contour1.collections + contour2.collections + contour3.collections

  ani = FuncAnimation(fig, update, frames=Nt, interval=interval, blit=False)

  save_path = os.path.join(Q_dir, name)
  ani.save(save_path, writer='Pillow')


def plot_energy_contribution(Q_dir, Sigma, freq, num_freq, name):
  # Sigma[Nf, Nt]
  fig, ax = plt.subplots(num_freq, figsize=(8, 8))
  for i in range(num_freq):
    #wl = 1.e3 / freq[i+1]
    ax[i].plot(range(1,11), Sigma[i+3,:10].real / np.sum(Sigma[i+3,:].real) * 100, 'o-')
    ax[i].set_title(f'Frequency {freq[i+3]:.1f}')
    ax[i].set_xlabel('Mode Index')
    ax[i].set_ylabel('Energy (%)')

  fig.tight_layout()
  save_path = os.path.join(Q_dir, name)
  plt.savefig(save_path)
  plt.close()

