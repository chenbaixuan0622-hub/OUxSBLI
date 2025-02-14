import numpy as np
import os
from scipy.linalg import eigh, svd
from scipy.fft import rfft, irfft, rfftfreq
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, extract_number, Data
from mod.mod_plot import print_scalar


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


def plot_reconstruction(POD_dir, x, y, z, modes, time_coef, num_modes, mean=None):
  Nt = len(time_coef[:,0])
  nx = len(x)
  ny = len(y)
  nz = len(z)
  POD = np.dot(time_coef[:,:num_modes], modes[:,:num_modes].T)
  if mean is not None:
    POD = POD + mean[None,:]
  POD = np.reshape(POD, [Nt,nz,ny,nx])
  for i in range(Nt):
    name = "POD" + str(i+1).zfill(5)
    print_scalar(x, y, z, POD[i,:,:,:], POD_dir, "Mode", name)


def plot_pod_results(POD_dir, x, y, z, t, eigenvalues, modes, time_coefficients, num_modes):
  # energy
  plt.figure(figsize=(4, 4))
  plt.plot(range(1,11), eigenvalues[:num_modes*2] / np.sum(eigenvalues) * 100, 'o-')
  plt.title('Energy Contribution of Modes')
  plt.xlabel('Mode Index')
  plt.ylabel('Energy (%)')
  save_name = "energy_contribution.png"
  save_path = os.path.join(POD_dir, save_name)
  plt.savefig(save_path)
  plt.close()

  # space mode
  nx = len(x)
  ny = len(y)
  nz = len(z)
  for i in range(num_modes):
    u = np.reshape(modes[:,i], [nz, ny, nx])
    # write VTK file
    name = "Mode" + str(i+1)
    print_scalar(x, y, z, u, POD_dir, "Mode", name)

  save_name = "POD_time_coef"
  save_path = os.path.join(POD_dir, save_name)
  np.save(save_path, time_coefficients)


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

'''
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
'''

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


def make_grid(Lx1, Lx2, Ly1, Ly2, stridex, stridey, endT, Q_dir):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Nt = len(Q_files)
  t  = np.linspace(0.e0, endT, Nt)

  first_path = os.path.join(Q_dir, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  for j in range(Ny):
    if Ly1 < Y[j]:
      ny1 = j
      break
  ny2 = Ny
  for j in range(Ny):
    if Ly2 <= Y[j]:
      ny2 = j
      break
  indicesx = np.arange(nx1, min(nx2, Nx), stridex)
  indicesy = np.arange(ny1, min(ny2, Ny), stridey)
  nx = len(indicesx)
  ny = len(indicesy)
  x  = X[indicesx]
  y  = Y[indicesy]
  indicesx, indicesy = np.meshgrid(indicesx, indicesy)
  return Nx, Ny, Nz, indicesx, indicesy, x, y, t


def make_grid3D(Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez, endT, Q_dir):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Nt = len(Q_files)
  t  = np.linspace(0.e0, endT, Nt)

  first_path = os.path.join(Q_dir, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  for j in range(Ny):
    if Ly1 < Y[j]:
      ny1 = j
      break
  ny2 = Ny
  for j in range(Ny):
    if Ly2 <= Y[j]:
      ny2 = j
      break
  nz1 = int(Lz1 / Z[-1] * Nz)
  nz2 = int(Lz2 / Z[-1] * Nz)
  indicesx = np.arange(nx1, min(nx2, Nx), stridex)
  indicesy = np.arange(ny1, min(ny2, Ny), stridey)
  indicesz = np.arange(nz1, min(nz2, Nz), stridez)
  nx = len(indicesx)
  ny = len(indicesy)
  nz = len(indicesz)
  x  = X[indicesx]
  y  = Y[indicesy]
  z  = Z[indicesz]
  indicesx, indicesy, indicesz = np.meshgrid(indicesx, indicesy, indicesz)
  return Nx, Ny, Nz, indicesx, indicesy, indicesz, x, y, z, t


def make_data(Lx1, Lx2, Ly1, Ly2, stridex, stridey, endT, Q_dir):
  Nx, Ny, Nz, indicesx, indicesy, x, y, t = make_grid(Lx1, Lx2, Ly1, Ly2, stridex, stridey, endT, Q_dir)
  
  print("nx = ", len(x), " ny = ", len(y))
  D = np.zeros((len(x)*len(y), len(t)), dtype=np.float32)

  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = U[0,indicesy,indicesx]
    D[:,itr] = u.flatten()
    itr += 1
  return x, y, t, D


def make_data3D(Q_dir, endT, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez):
  data = Data(Q_dir, endT, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez)
  Nx, Ny, Nz, indicesx, indicesy, indicesz, x, y, z, t = data.make_grid()
  
  print("nx = ", len(x), " ny = ", len(y), " nz = ", len(z))
  D = np.zeros((len(x)*len(y)*len(z), len(t)), dtype=np.float32)

  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, _, _ = getVector(file_path, Nx, Ny, Nz, 'rotA')
    u = U[indicesz,indicesy,indicesx].transpose(2,0,1)
    D[:,itr] = u.flatten()
    itr += 1
  return x, y, z, t, D


def make_data_combine(Q_dir, endT, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, stridex, stridey, stridez):
  data1 = Data(Q_dir[0], endT, Lx1=Lx1, Lx2=None, Ly1=Ly1, Ly2=Ly2, Lz1=Lz1, Lz2=Lz2, \
               stridex=1, stridey=stridey, stridez=stridez)
  data2 = Data(Q_dir[1], endT, Lx1=None, Lx2=Lx2, Ly1=Ly1, Ly2=Ly2, Lz1=Lz1, Lz2=Lz2, \
               stridex=1, stridey=stridey, stridez=stridez)
  Nx1, Ny, Nz, indicesx1, indicesy1, indicesz1, x1, y, z, t = data1.make_grid()
  Nx2, _,  _,  indicesx2, indicesy2, indicesz2, x2, _, _, _ = data2.make_grid()
 
  def read_data(Q_dir, t, Nx, Ny, Nz, indicesx, indicesy, indicesz):
    Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
    Q_files.sort(key=extract_number)
    ny, nx, nz = indicesx.shape

    U = np.zeros((len(t), nz, ny, nx), dtype=np.float32)
    V = np.zeros((len(t), nz, ny, nx), dtype=np.float32)
    W = np.zeros((len(t), nz, ny, nx), dtype=np.float32)
    
    itr = 0
    for Q_file in tqdm(Q_files):
      file_path = os.path.join(Q_dir, Q_file)
      u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
      U[itr,:,:,:] = u[indicesz,indicesy,indicesx].transpose(2,0,1)
      V[itr,:,:,:] = v[indicesz,indicesy,indicesx].transpose(2,0,1)
      W[itr,:,:,:] = w[indicesz,indicesy,indicesx].transpose(2,0,1)
      itr += 1
    return U, V, W
  
  U1, V1, W1 = read_data(Q_dir[0], t, Nx1, Ny, Nz, indicesx1, indicesy1, indicesz1)
  U2, V2, W2 = read_data(Q_dir[1], t, Nx2, Ny, Nz, indicesx2, indicesy2, indicesz2)
  U = np.concatenate([U1, U2[:,:,:,3:]], axis=3)
  V = np.concatenate([V1, V2[:,:,:,3:]], axis=3)
  W = np.concatenate([W1, W2[:,:,:,3:]], axis=3)
  U = U[:,:,:,::stridex]
  V = V[:,:,:,::stridex]
  W = W[:,:,:,::stridex]
  x = np.linspace(x1[0], x2[-1], U.shape[-1])
  return x, y, z, t, U, V, W

