import numpy as np
import os
import matplotlib.pyplot as plt
from mod.mod_read import getGrid, getVector
from mod.mod_plot import print_scalar, print_vector
from mod.mod_Poisson import Poisson_CND, Poisson_Spectral


def simple_test():
  nx = 129
  ny = 129
  nz = 129
  Lx = 2.e0
  Ly = 2.e0
  Lz = 2.e0
  dx = Lx / (nx-1)
  dy = Ly / (ny-1)
  dz = Lz / (nz-1)
  x  = np.linspace(-Lx/2 + 0.5e0 * dx, Lx/2 - 0.5e0 * dx, nx, dtype=np.float64)
  y  = np.linspace(-Ly/2 + 0.5e0 * dy, Ly/2 - 0.5e0 * dy, ny, dtype=np.float64)
  z  = np.linspace(-Lz/2 + 0.5e0 * dz, Lz/2 - 0.5e0 * dz, nz, dtype=np.float64)
  xc = 0.5e0 * (x[:-1] + x[1:])
  yc = 0.5e0 * (y[:-1] + y[1:])
  zc = 0.5e0 * (z[:-1] + z[1:])
  f3D = np.zeros((nz-1,ny-1,nx-1), dtype=np.float64)
  for k in range(nz-1):
    for j in range(ny-1):
      for i in range(nx-1):
        f3D[k,j,i] = np.exp(-10.e0*(xc[i]**2 + yc[j]**2 + zc[k]**2))
  #p3D = Poisson_CND(f3D, dx=dx, loop_x='C', dy=dy, loop_y='C', dz=dz, loop_z='C', check=True)
  p3D = Poisson_Spectral(f3D, dx, dy, dz)
  print_scalar(xc, yc, zc, p3D, '.', 'T', 'T')


def Helmholtz_Decomposition(x, y, z, u, v, w, loop_x, loop_y, loop_z):
  # 2nd order accuracy solver
  # x[nx+2], y[ny+2], z[nz+2], u[nz+2,ny+2,nx+2], v[nz+2,ny+2,nx+2], w[nz+2,ny+2,nx+2]
  # div[nz,ny,nx], phi[nz,ny,nx]
  dx  = -x[0] + x[1]
  dy  = -y[0] + y[1]
  dz  = -z[0] + z[1]
  div = 0.5e0 * (-u[1:-1,1:-1,:-2] + u[1:-1,1:-1,2:]) / dx \
      + 0.5e0 * (-v[1:-1,:-2,1:-1] + v[1:-1,2:,1:-1]) / dy \
      + 0.5e0 * (-w[:-2,1:-1,1:-1] + w[2:,1:-1,1:-1]) / dz
  U   = np.stack((u[1:-1,1:-1,1:-1], v[1:-1,1:-1,1:-1], w[1:-1,1:-1,1:-1]), axis=-1)
  phi = Poisson_CND(div, dx=dx, loop_x=loop_x, dy=dy, loop_y=loop_y, dz=dz, loop_z=loop_z, check=True)
  Phi = np.zeros_like(u)
  
  def set_bc_x(Phi, bc_type):
    if bc_type == 'C' or bc_type == 'c':
      Phi[:,:,0]  = Phi[:,:,-2]
      Phi[:,:,-1] = Phi[:,:,1]
    elif bc_type == 'N' or bc_type == 'n':
      Phi[:,:,0]  = Phi[:,:,1]
      Phi[:,:,-1] = Phi[:,:,-2]
    elif bc_type == 'D' or bc_type == 'd':
      Phi[:,:,0]  = 0.e0
      Phi[:,:,-1] = 0.e0
    else:
      print("Invalid arg")
    return Phi

  def set_bc_y(Phi, bc_type):
    if bc_type == 'C' or bc_type == 'c':
      Phi[:,0,:]  = Phi[:,-2,:]
      Phi[:,-1,:] = Phi[:,1,:]
    elif bc_type == 'N' or bc_type == 'n':
      Phi[:,0,:]  = Phi[:,1,:]
      Phi[:,-1,:] = Phi[:,-2,:]
    elif bc_type == 'D' or bc_type == 'd':
      Phi[:,0,:]  = 0.e0
      Phi[:,-1,:] = 0.e0
    else:
      print("Invalid arg")
    return Phi
  
  def set_bc_z(Phi, bc_type):
    if bc_type == 'C' or bc_type == 'c':
      Phi[0,:,:]  = Phi[-2,:,:]
      Phi[-1,:,:] = Phi[1,:,:]
    elif bc_type == 'N' or bc_type == 'n':
      Phi[0,:,:]  = Phi[1,:,:]
      Phi[-1,:,:] = Phi[-2,:,:]
    elif bc_type == 'D' or bc_type == 'd':
      Phi[0,:,:]  = 0.e0
      Phi[-1,:,:] = 0.e0
    else:
      print("Invalid arg")
    return Phi
  
  Phi[1:-1,1:-1,1:-1] = phi
  Phi = set_bc_x(Phi, bc_type=loop_x)
  Phi = set_bc_y(Phi, bc_type=loop_y)
  Phi = set_bc_z(Phi, bc_type=loop_z)
  gradPhix = 0.5e0 * (-phi[1:-1,1:-1,:-2] + phi[1:-1,1:-1,2:]) / dx
  gradPhiy = 0.5e0 * (-phi[1:-1,:-2,1:-1] + phi[1:-1,2:,1:-1]) / dy
  gradPhiz = 0.5e0 * (-phi[:-2,1:-1,1:-1] + phi[2:,1:-1,1:-1]) / dz
  gradPhi  = np.stack((gradPhix, gradPhiy, gradPhiz), axis=-1)
  rotA     = U[1:-1,1:-1,1:-1,:] + gradPhi
  return gradPhi, rotA


simple_test()

'''
#Q_dir = "../../SBLI"
Q_dir = "../3D_solver/ETGV/data"
Nx, Ny, Nz, X, Y, Z = getGrid(os.path.join(Q_dir, "Q00100.vtr"))
u, v, w    = getVector(os.path.join(Q_dir, "Q00100.vtr"), Nx, Ny, Nz, 'velocity')
#um, vm, wm = getVector(os.path.join(Q_dir, "Qmean.vtr"), Nx, Ny, Nz, 'velocity')
#Uf = np.stack((u - um, v - vm, w - wm), axis=-1)


gradPhi, rotA = Helmholtz_Decomposition(X, Y, Z, u, v, w, loop_x='C', loop_y='C', loop_z='C')
print_vector(X[2:-2], Y[2:-2], Z[2:-2], gradPhi, '.', 'gradPhi', 'gradPhi')
print_vector(X[2:-2], Y[2:-2], Z[2:-2], rotA, '.', 'rotA', 'rotA')
'''
