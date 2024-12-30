import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def PSD(data, t):
  n    = len(t)
  dt   = -t[0] + t[1]
  fs   = 1.e0 / dt
  df   = fs / float(n)
  fft  = np.fft.fft(data)
  freq = np.fft.fftfreq(n, d=dt)
  amp  = abs(fft / (0.5e0 * float(n)))
  ps   = amp**2
  psd  = ps / df
  return freq, psd


def low_pass(U, dz, lz):
  fft  = np.fft.fft(U)
  Nz   = len(U)
  freq = np.linspace(0.e0, 1.e0/dz, Nz)
  fc   = 1.e0 / lz
  fft[(freq > fc)] = 0.e0
  u = np.fft.ifft(fft)
  return u.real


def high_pass(U, dz, lz):
  fft  = np.fft.fft(U)
  Nz   = len(U)
  freq = np.linspace(0.e0, 1.e0/dz, Nz)
  fc   = 1.e0 / lz
  fft[(freq < fc)] = 0.e0
  u = np.fft.ifft(fft)
  return u.real


def band_pass(U, dz, lz1, lz2):
  fft  = np.fft.fft(U)
  Nz   = len(U)
  freq = np.linspace(0.e0, 1.e0/dz, Nz)
  fc1  = 1.e0 / lz1
  fc2  = 1.e0 / lz2
  fft[(fc1 < freq) and (freq < fc2)] = 0.e0
  u = np.fft.ifft(fft)
  return u.real


def getFourierMode(u, v, w):
  # n: streamwise wavenumber
  # m: spanwise wavenumber
  # fft_img[streamwise wavenumber, spanwise wavenumber]
  fft_u = np.fft.fft2(u, s=[2,2], norm="forward")
  fft_u = np.fft.fftshift(fft_u)
  fft_v = np.fft.fft2(v, s=[2,2], norm="forward")
  fft_v = np.fft.fftshift(fft_v)
  fft_w = np.fft.fft2(w, s=[2,2], norm="forward")
  fft_w = np.fft.fftshift(fft_w)
  # fft_img[0]: mode0
  # fft_img[1]: mode1

  # v[u[0,1], v[0,1], w[1,0], u[1,1], v[1,1], w[1,1]]
  fft_mode = np.zeros((3,2,2), dtype=complex)

  for j in range(2):
    for i in range(2):
      fft_mode[0,j,i] = fft_u[j,i]
      fft_mode[1,j,i] = fft_u[j,i]
      fft_mode[2,j,i] = fft_u[j,i]
  return fft_mode


def IFFT(fft_U, index, U):
  fft_u = np.zeros_like(fft_U, dtype=complex)
  fft_u[index] = fft_U[index]
  fft_u = np.fft.ifftshift(fft_u)
  u = np.fft.ifft2(fft_u, s=U.shape, norm="forward")
  return u.real


def IFFT11(fft_U, U):
  fft_u = np.zeros_like(fft_U, dtype=complex)
  fft_u[1,0] = fft_U[1,0]
  fft_u[0,1] = fft_U[0,1]
  fft_u = np.fft.ifftshift(fft_u)
  u = np.fft.ifft2(fft_u, s=U.shape, norm="forward")
  return u.real


def getV(u, v, w):
  # V[u[0,1], v[0,1], w[1,0], u[1,1], v[1,1], w[1,1]]
  V = np.zeros(6, dtype=np.float32)

  fft_mode = getFourierComponent(u, v, w)

  V[0] = abs(fft_mode[0,0,1])**2
  V[1] = abs(fft_mode[1,0,1])**2
  V[2] = abs(fft_mode[2,1,0])**2
  V[3] = abs(fft_mode[0,1,1])**2
  V[4] = abs(fft_mode[1,1,1])**2
  V[5] = abs(fft_mode[2,1,1])**2
  return V


def reconstruct_u(x, z, U, V, W):
  fft_mode = getFourierMode(U, V, W)
 
  u00 = IFFT(fft_mode[0,:,:], (0,0), U)
  u01 = IFFT(fft_mode[0,:,:], (0,1), U)
  u10 = IFFT(fft_mode[0,:,:], (1,0), U)
  u11 = IFFT(fft_mode[0,:,:], (1,1), U)
  U11 = IFFT11(fft_mode[0,:,:], U)

  x, z = np.meshgrid(x, z)

  plt.figure(figsize=(12, 12))
  plt.subplot(321)
  plt.title("u")
  plt.axis('equal')
  plt.contourf(x, z, U, cmap='jet')
  plt.colorbar()

  plt.subplot(322)
  plt.title("u00")
  plt.axis('equal')
  plt.contourf(x, z, u00, cmap='jet')
  plt.colorbar()

  plt.subplot(323)
  plt.title("u01")
  plt.axis('equal')
  plt.contourf(x, z, u01, cmap='jet')
  plt.colorbar()

  plt.subplot(324)
  plt.title("u10")
  plt.axis('equal')
  plt.contourf(x, z, u10, cmap='jet')
  plt.colorbar()

  plt.subplot(325)
  plt.title("u11")
  plt.axis('equal')
  plt.contourf(x, z, u11, cmap='jet')
  plt.colorbar()

  plt.subplot(326)
  plt.title("u11")
  plt.axis('equal')
  plt.contourf(x, z, U11, cmap='jet')
  plt.colorbar()

  plt.show()
  plt.close()


def calc_FFT_save_animation(Nk, Nt, Nz, z, data, dir, name):
  # store data for animation and calc FFT
  dz   = -z[0] + z[1]
  freq = np.fft.fftfreq(Nz, d=dz)
  amp  = np.zeros((Nk,Nz))

  for k in range(Nk):
    for i in range(Nt):
      amp[k,:] += abs(np.fft.fft(data[k,i,:]) / (Nz / 2.e0))

  amp /= Nt

  # for plot
  plt.rcParams['font.family'] = 'Times New Roman'
  plt.rcParams['mathtext.fontset'] = 'stix'
  plt.rcParams['xtick.direction'] = 'in'
  plt.rcParams['ytick.direction'] = 'in'
  plt.rcParams['font.size'] = 12

  # plot FFT
  plt.xlabel("Wavelength [mm]", fontsize='14', style='italic', color='black')
  plt.ylabel("Amp [m/s]",  fontsize='14', style='italic', color='black')
  colors = ['black', 'blue', 'red', 'green', 'purple', 'orange', 'pink']
  for k in range(Nk):
    plt.plot(1.e3 / freq[1:Nz//2], amp[k,1:Nz//2], color=colors[k])
  save_path = os.path.join(dir, name+".png")
  plt.savefig(save_path)
  plt.show()
  plt.close()

  # plot animation
  ims = []

  fig = plt.figure()

  plt.xlabel("z [mm]", fontsize='14', style='italic', color='black')
  plt.ylabel("U [m/s]", fontsize='14', style='italic', color='black')

  for i in range(Nt):
    lines = []
    for k in range(Nk):
      line, = plt.plot(z*1e3, data[k,i,:], color=colors[k])
      lines.append(line)
    ims.append(lines)

  ani = animation.ArtistAnimation(fig, ims, interval=200)
  plt.show()
  save_path = os.path.join(dir, name+".gif")
  ani.save(save_path, writer='Pillow')
  plt.close()

