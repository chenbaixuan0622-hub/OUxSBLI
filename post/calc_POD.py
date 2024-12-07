import numpy as np
import cv2
import os
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_POD import snapshot_pod, calc_time_coef


# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12


def plot_pod_results(Q_dir, x, z, data, eigenvalues, modes, time_coefficients, num_modes=3):
    # energy
    plt.figure(figsize=(6, 4))
    plt.plot(eigenvalues[:10] / np.sum(eigenvalues) * 100, 'o-')
    plt.title('Energy Contribution of Modes')
    plt.xlabel('Mode Index')
    plt.ylabel('Energy (%)')
    save_name = "energy_contribution.png"
    save_path = os.path.join(Q_dir, save_name)
    plt.savefig(save_path)
    plt.close()

    # space mode
    nx = len(x)
    nz = len(z)
    x  = x * 1e3
    z  = z * 1e3
    xmin = x[0]
    xmax = x[-1]
    zmin = z[0]
    zmax = z[-1]
    x, z = np.meshgrid(x, z)
    fig, ax = plt.subplots(num_modes//2, 2, figsize=(8, 8))
    for i in range(num_modes):
      u = np.reshape(modes[:,i], [nz, nx])
      ax[i//2,i%2].set_xlim(xmin, xmax)
      ax[i//2,i%2].set_ylim(zmin, zmax)
      ax[i//2,i%2].set_aspect('equal', adjustable='box')
      ax[i//2,i%2].set_title(f'POD MODE {i+1}')
      im = ax[i//2,i%2].contourf(x, z, u, levels=50, cmap='jet', extend='both')
      divider = make_axes_locatable(ax[i//2,i%2])
      cax = divider.append_axes('right', '5%', pad='3%')
      fig.colorbar(im, cax=cax, extendrect=True)
      fig.tight_layout()

    save_name = "Mode.png" 
    save_path = os.path.join(Q_dir, save_name)
    plt.savefig(save_path)
    plt.close()

    # time coefficient
    fig, ax = plt.subplots(num_modes, 1, figsize=(8, 8))
    for i in range(num_modes):
      ax[i].plot(time_coefficients[:,i])
    
    save_name = "time_coef.png"
    save_path = os.path.join(Q_dir, save_name)
    plt.savefig(save_path)
    plt.close()


def plot_reconstruction(Q_dir, x, z, data, modes, time_coef, num_modes=3, interval=200):
  Nt   = len(time_coef[:,0])
  nx   = len(x)
  nz   = len(z)
  x    = x * 1e3
  z    = z * 1e3
  x, z = np.meshgrid(x, z)
  u    = np.reshape(data, [nz,nx,Nt])

  POD  = np.dot(time_coef[:,:num_modes], modes[:,:num_modes].T)
  POD  = np.reshape(POD, [Nt,nz,nx])

  fig, ax  = plt.subplots(1, 3, figsize=(18, 6))
  crange   = np.linspace(0, 500, 50)
  contour1 = ax[0].contourf(x, z,   u[:,:,0], crange, cmap='jet', extend='both')
  contour2 = ax[1].contourf(x, z, POD[0,:,:], crange, cmap='jet', extend='both')
  contour3 = ax[2].contourf(x, z, POD[0,:,:] - u[:,:,0], cmap='jet')
  cbar1    = fig.colorbar(contour1, ax=ax[0:1], extendrect=True, \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar1.set_label("u [m/s]")
  cbar1.ax.xaxis.set_ticks_position('top')
  cbar1.ax.xaxis.set_label_position('top')
  cbar2    = fig.colorbar(contour3, ax=ax[2], extendrect=True, \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar2.set_label("POD - u [m/s]")
  cbar2.ax.xaxis.set_ticks_position('top')
  cbar2.ax.xaxis.set_label_position('top')

  def update(frame):
    for a in ax:
      a.clear()
      a.set_aspect('equal', adjustable='box')
    contour1 = ax[0].contourf(x, z,   u[:,:,frame], crange, cmap="jet", extend='both')
    contour2 = ax[1].contourf(x, z, POD[frame,:,:], crange, cmap='jet', extend='both')
    contour3 = ax[2].contourf(x, z, POD[frame,:,:] - u[:,:,frame], cmap='jet')
    return contour1.collections + contour2.collections + contour3.collections

  ani = FuncAnimation(fig, update, frames=Nt, interval=interval, blit=False)

  save_name = "POD.gif"
  save_path = os.path.join(Q_dir, save_name)
  ani.save(save_path, writer='Pillow')


def main(Lx1, Lx2):
  Q_dir   = "../../yp269"
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Nt      = len(Q_files)

  first_path = os.path.join(Q_dir, Q_files[0])
  Nx, Ny, Nz, X, Y, Z = getGrid(first_path)
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  nx  = (nx2 - nx1) // 4
  nz  = Nz // 4
  
  x   = np.linspace(X[nx1], X[nx2], nx)
  z   = np.linspace(Z[0],   Z[-1],  nz)

  data = np.zeros((nx*nz,Nt), dtype=np.float32)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    U = np.reshape(U[:,0,nx1:nx2], [Nz, nx2-nx1])
    u = cv2.resize(U, (nx, nz))
    data[:,itr] = u.flatten()
    itr += 1

  eigenvalues, eigenvectors, modes = snapshot_pod(data)
  time_coef = calc_time_coef(data, modes)

  plot_pod_results(Q_dir, x, z, data, eigenvalues, modes, time_coef, num_modes=6)
  plot_reconstruction(Q_dir, x, z, data, modes, time_coef, num_modes=6, interval=200)

Lx1 = 25.e-3
Lx2 = 35.e-3
main(Lx1, Lx2)

