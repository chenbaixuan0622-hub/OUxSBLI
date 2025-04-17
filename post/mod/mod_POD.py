import numpy as np
import os
from scipy.linalg import eigh, svd
from scipy.fft import rfft, irfft, rfftfreq
from scipy.signal import welch
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, extract_number, Data
from mod.mod_plot import print_scalar


def snapshot_pod(D):
  # D[space, time]
  C_s  = np.dot(D.T, D)
  Sigma, U = eigh(C_s)
  del C_s
  idx   = Sigma.argsort()[::-1]
  Sigma = Sigma[idx]
  U     = U[:,idx]
  phi   = np.dot(D, U)
  phi  /= np.linalg.norm(phi, axis=0)
  return Sigma, U, phi


def spectral_pod(D, dir, dt, window_size, overlap=0.5e0):
  # D[space, time]
  _, Nt = D.shape
  step  = int((1.e0 - overlap) * window_size)
  freq  = rfftfreq(window_size, d=dt)
  Nb    = (Nt - window_size) // step + 1
  w     = np.array(np.hamming(window_size), dtype=np.float32)
  for j in tqdm(range(len(freq))):
    # window[Nx,Nb,window_size]
    window = np.array(
      [D[:,i*step:i*step+window_size] * w[None,:] for i in range(Nb)], dtype=np.float32
    ).transpose(1,0,2)
    Df = rfft(window, axis=-1)[:,:,j]
    del window
    Sigma, U = eigh(Df.conj().T @ Df / Nb)
    idx   = Sigma.argsort()[::-1]
    Sigma = Sigma[idx]
    U     = U[:,idx]
    Sigma_diag = np.sqrt(Nb) * np.sqrt(Sigma)
    phi = Df @ (U * (1.e0 / Sigma_diag))
    A   = Sigma_diag @ U.conj().T
    file_path = os.path.join(dir, "SPOD_" + str(j).zfill(4))
    np.savez(file_path, Sigma=Sigma, U=U, phi=phi, A=A)


def calc_time_coef(D, modes):
  time_coef = np.dot(D.T, modes)
  return time_coef


def calc_time_coef_conv(dir, D, num_modes, window_size, dt):
  files = [f for f in os.listdir(dir) if f.endswith(".npz")]
  files.sort(key=extract_number)
  Df     = rfft(D[:,:window_size], axis=-1)
  Nx, Nf = Df.shape
  phi = np.zeros((Nx,num_modes,Nf), dtype=np.complex64)
  for f, file in tqdm(enumerate(files)):
    SPOD = np.load(os.path.join(dir, file))
    phi[:,:,f] = SPOD['phi'][:,:num_modes]
  a = np.einsum('xrf,xf->rf', phi.conj(). Df)
  a = irfft(a, axis=-1)
  file_path = os.path.join(dir, "time_coef_conv")
  np.save(file_path, a)


def reconstruct_POD(POD_dir, x, y, z, modes, time_coef, num_modes, mean=None):
  Nt  = time_coef.shape[0]
  POD = np.dot(time_coef[:,:num_modes], modes[:,:num_modes].T)
  if mean is not None:
    POD = POD + mean[None,:]
  POD = np.reshape(POD, [Nt,len(z),len(y),len(x)])
  for i in tqdm(range(Nt)):
    name = "POD" + str(i+1).zfill(5)
    print_scalar(x, y, z, POD[i,:,:,:], POD_dir, "Mode", name)


# reconstruct in the frequency domain
def reconstruct_spod(dir, x, y, z, num_modes, mean, window_size):
  files = [f for f in os.listdir(dir) if f.endswith(".npz")]
  files.sort(key=extract_number)
  D = np.zeros((len(x)*len(y)*len(z),len(files)), dtype=np.complex64)
  for i, f in tqdm(enumerate(files)):
    SPOD = np.load(os.path.join(dir, f))
    A    = SPOD['A']
    phi  = SPOD['phi']
    for mode in range(num_modes):
      D[:,i] += A[mode] * phi[:,mode]
  D  = irfft(D, axis=-1) / np.hamming(window_size)[None,:]
  Nt = D.shape[-1]
  D += mean[:,None]
  D  = np.reshape(D, [len(z),len(y),len(x),Nt])
  for i in tqdm(range(Nt)):
    name = "SPOD" + str(i+1).zfill(5)
    print_scalar(x, y, z, D[:,:,:,i], dir, "Mode", name)


def plot_pod_results(POD_dir, x, y, z, t, Sigma, phi, time_coef, num_modes):
  filename  = "energy_contribution.d"
  save_path = os.path.join(POD_dir, filename)
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# Mode       energy", file=f)
    for i in range(len(Sigma)):
      print(f'{i+1:.0e}', f'{Sigma[i]/np.sum(Sigma)*100:.3e}', file=f)
  # space mode
  nx = len(x)
  ny = len(y)
  nz = len(z)
  for i in range(num_modes):
    u = np.reshape(phi[:,i], [nz, ny, nx])
    name = "Mode" + str(i+1)
    print_scalar(x, y, z, u, POD_dir, "Mode", name)
  save_name = "POD_time_coef"
  save_path = os.path.join(POD_dir, save_name)
  np.save(save_path, time_coef)
  dt = -t[0] + t[1]
  for mode in range(num_modes):
    freq, psd = welch(time_coef[:,mode], fs=1.e0 / dt, nperseg=len(t)//4)
    filename = "time_coef" + str(mode).zfill(4) + ".d"
    save_path = os.path.join(POD_dir, filename)
    with open(save_path, "w", encoding="UTF-8") as f:
      print('# time[ms]   time_coef', file=f)
      for i in range(len(t)):
        print(f'{t[i]*1e3:.3e}', f'{time_coef[i,mode]:.3e}', file=f)
    filename = "PSD" + str(mode).zfill(4) + ".d"
    save_path = os.path.join(POD_dir, filename)
    with open(save_path, "w", encoding="UTF-8") as f:
      print("# freq       psd", file=f)
      for i in range(len(freq)):
        print(f'{freq[i]:.3e}', f'{psd[i]:.3e}', file=f)


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

