import numpy as np
import matplotlib.pyplot as plt

def FourierComponent(u, v, w, n, m):
  # n: streamwise wavenumber
  # m: spanwise wavenumber
  # fft_img[streamwise wavenumber, spanwise wavenumber]
  fft_u = np.fft.fft2(u, s=[2,2])
  fft_u = np.fft.fftshift(fft_u)
  fft_v = np.fft.fft2(v, s=[2,2])
  fft_v = np.fft.fftshift(fft_v)
  fft_w = np.fft.fft2(w, s=[2,2])
  fft_w = np.fft.fftshift(fft_w)
  # fft_img[0]: mean
  # fft_img[1]: fluctuating

  # V[u[0,1], v[0,1], w[1,0], u[1,1], v[1,1], w[1,1]]
  V = np.zeros(6, dtype=np.float32)
  V[0] = abs(fft_u[0,1])**2
  V[1] = abs(fft_v[0,1])**2
  V[2] = abs(fft_w[1,0])**2
  V[3] = abs(fft_u[1,1])**2
  V[4] = abs(fft_v[1,1])**2
  V[5] = abs(fft_w[1,1])**2
  return V

def reconstruct_energy():
  return K



  '''
  spanwise_amp = amp.sum(axis=1)
  print(spanwise_amp)
  spanwise_indices = np.argsort(spanwise_amp)[::-1]
  nth_strongest_index = spanwise_indices[n-1]

  streamwise_amp = amp[nth_strongest_index,:]
  print(streamwise_amp)
  streamwise_indices = np.argsort(streamwise_amp)[::-1]
  mth_strongest_index = streamwise_indices[m-1]
  '''

  #component_index = (nth_strongest_index, mth_strongest_index)
  isolated_fft_img = np.zeros_like(fft_img, dtype=complex)
  isolated_fft_img[0,1] = fft_img[0,1]
  isolated_fft_img[1,1] = fft_img[1,1]
  return isolated_fft_img

def display_selected_frequency_component(x, z, img, n, m):
  fft_img = FourierComponent2d(x, z, img, n, m)

  print(fft_img)

  fft_img = np.fft.ifftshift(fft_img)
  inverse_img = np.fft.ifft2(fft_img, s=x.shape)
  inverse_img = np.abs(inverse_img)

  plt.figure(figsize=(12,6))
  plt.subplot(1, 2, 1)
  plt.title("Original Image")
  plt.contourf(x, z, img, cmap="jet")
  plt.axis('off')

  plt.subplot(1, 2, 2)
  plt.title("Reconstruct Image")
  plt.contourf(x, z, inverse_img, cmap="jet")
  plt.axis('off')

  plt.show()

n, m = 1, 1

Nx = 32
Nz = 16
x = np.linspace(0.e0, 2.e0 * np.pi, Nx)
z = np.linspace(0.e0, 2.e0 * np.pi, Nz)
img = np.zeros((Nz,Nx))
for k in range(Nz):
  for i in range(Nx):
    #img[k,i] = np.abs(np.sin(x[i]) + np.cos(z[k]))
    img[k,i] = np.sin(z[k])
#img = img - np.mean(img)

X, Z = np.meshgrid(x, z)
plt.title("Original")
plt.contourf(X, Z, img, cmap="jet")
plt.show()
plt.close()

display_selected_frequency_component(X, Z, img, n, m)

