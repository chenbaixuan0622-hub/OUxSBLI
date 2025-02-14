import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number
#from mod.mod_matrix import Poisson_x_Dirichlet_y_Neumann_z_cyclic
from mod.mod_Poisson import Poisson
from mod.mod_nabla import Scalar, Vector, bc
from mod.mod_Helmholtz import save_Helmholtz_decomposition


set_Params()

# parameter
# shock
# Lx1 = 24.e-3
# Lx2 = 40.e-3
# boundary layer
# Lx1 =  4.e-3
# Lx2 = 20.e-3
Q_dir   = "../../../../../../media/user/HD-EDS-E/TBL/SBLI_1delta_test"
Lx1     = 24.e-3
Lx2     = 40.e-3
Ly1     = 0.e-3
Ly2     = 10.e-3
stridex = 4
stridez = 4
endT    = 0.1e-3


def main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, endT):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  
  Nx, Ny, Nz, X, Y, Z = getGrid(os.path.join(Q_dir, Q_files[0]))
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  indicesx = np.arange(nx1,  nx2, stridex)
  indicesy = np.arange(0, len(Y), 1)
  indicesz = np.arange(0, len(Z), stridez)
  x = X[indicesx]
  z = Z[indicesz]
  nx = len(x)
  ny = len(Y)
  nz = len(z)

  indicesx, indicesy, indicesz = np.meshgrid(indicesx, indicesy, indicesz)

  itr = 0
  # calc mean value
  umean = np.zeros((nz,ny,nx), dtype=np.float32)
  vmean = np.zeros((nz,ny,nx), dtype=np.float32)
  wmean = np.zeros((nz,ny,nx), dtype=np.float32)
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = U[indicesz,indicesy,indicesx].transpose(2,0,1)
    v = V[indicesz,indicesy,indicesx].transpose(2,0,1)
    w = W[indicesz,indicesy,indicesx].transpose(2,0,1)
    umean += u
    vmean += v
    wmean += w
    itr += 1
  umean /= itr 
  vmean /= itr
  wmean /= itr

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = U[indicesz,indicesy,indicesx].transpose(2,0,1)
    v = V[indicesz,indicesy,indicesx].transpose(2,0,1)
    w = W[indicesz,indicesy,indicesx].transpose(2,0,1)
    u -= umean
    v -= vmean
    w -= wmean

    y  = np.linspace(Ly1, Ly2, 100, dtype=np.float32)
    Ui = np.zeros((3,len(z),len(y),len(x)), dtype=np.float32)

    for k in range(len(z)):
      for i in range(len(x)):
        Ui[0,k,:,i] = np.interp(y, Y, u[k,:,i])
        Ui[1,k,:,i] = np.interp(y, Y, v[k,:,i])
        Ui[2,k,:,i] = np.interp(y, Y, w[k,:,i])
 
    '''
    Div  = Vector(Ui, x, y, z).divergence(periodic_z=True)
    div  = Div[:,:-1,1:-1]
    x    = x[1:-1]
    y    = y[:-1]
    PDE  = Poisson(div, x, y, z)

    dx   = -x[0] + x[1]
    dpx1 = np.zeros_like(div[:,:,0])#Ui[0,:,:-1,0] * dx# * Fin[None,:,0]
    dpx2 = np.zeros_like(div[:,:,0])#Ui[0,:,:-1,-1] * dx# * Fout[None,:,0]
    dpy1 = np.zeros_like(div[:,0,:])
    dpy2 = np.zeros_like(div[:,0,:])
    phi  = PDE.x_Neumann_y_Neumann_z_periodic(dpx1, dpx2, dpy1, dpy2)
    name = "HD" + str(itr).zfill(5)
    HD_dir = os.path.join(Q_dir, "HD")
    os.makedirs(HD_dir, exist_ok=True)
    save_Helmholtz_decomposition(Ui[:,:,:-1,1:-1], phi, x, y, z, HD_dir, name)
    itr += 1
    '''

    U0  = Ui
    Phi = np.zeros_like(Ui[0,:,:,:])
    A   = np.zeros_like(Ui)
    for i_try in range(10):
      # gradPhi
      div  = Vector(Ui, x, y, z).divergence(periodic_z=True)
      div  = bc(div).Neumann(x1=True, x2=True, y2=True)
      div  = bc(div).Dirichlet(y1=0.e0)
      PDE1 = Poisson(div, x, y, z)

      dx   = -x[0] + x[1]
      dpx1 = np.zeros_like(div[:,:,0])
      dpx2 = np.zeros_like(div[:,:,0])
      dpy1 = np.zeros_like(div[:,0,:])
      dpy2 = np.zeros_like(div[:,0,:])
      phi  = PDE1.x_Dirichlet_y_Neumann_z_periodic(dpx1, dpx2, dpy1, dpy2)
      gradPhi = Scalar(phi, x, y, z).gradient()
      for i in range(3):
        gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).periodic(z=True)
        gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).Neumann(x1=True, x2=True, y2=True)
        gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).Dirichlet(y1=0.e0)

      '''
      # rotA
      rot  = Vector(Ui, x, y, z).rotation(periodic_z=True)
      for i in range(3):
        rot[i,:,:,:]  = bc(rot[i,:,:,:]).Neumann(x1=True, x2=True, y2=True)
        rot[i,:,:,:]  = bc(rot[i,:,:,:]).Dirichlet(y1=0.e0)
      PDEx = Poisson(rot[0,:,:,:], x, y, z)
      PDEy = Poisson(rot[1,:,:,:], x, y, z)
      PDEz = Poisson(rot[2,:,:,:], x, y, z)

      Ax1  = np.zeros_like(rot[0,:,:,0])
      Ax2  = np.zeros_like(rot[0,:,:,0])
      dAy1 = np.zeros_like(rot[0,:,0,:])
      dAy2 = np.zeros_like(rot[0,:,0,:])
      Ay1  = np.zeros_like(rot[0,:,0,:])
      Ay2  = np.zeros_like(rot[0,:,0,:])
      Ax   = PDEx.x_Dirichlet_y_Dirichlet_z_periodic(Ax1, Ax2, Ay1,   Ay2)
      Ay   = PDEy.x_Dirichlet_y_Neumann_z_periodic(Ax1,   Ax2, dAy1, dAy2)
      Az   = PDEz.x_Dirichlet_y_Dirichlet_z_periodic(Ax1, Ax2, Ay1,   Ay2)
      rotA = Vector(np.array([Ax, Ay, Az]), x, y, z).rotation(periodic_z=True) 
      '''
      # update based on gradPhi
      ud   = gradPhi
      ur   = U0 - ud
      Phi += phi
      Ui  -= gradPhi
      '''
      # update based on rotA
      ur   = rotA
      ud   = U0 - ur
      A   += np.array([Ax, Ay, Az])
      Ui  -= rotA
      '''
      rotud = Vector(ud, x, y, z).rotation()
      divur = Vector(ur, x, y, z).divergence()
      print('divrotA = ', np.mean(np.abs(divur)), ' rotgradPhi = ', np.mean(np.abs(rotud)))

    name = "HD" + str(itr).zfill(5)
    HD_dir = os.path.join(Q_dir, "HD")
    os.makedirs(HD_dir, exist_ok=True)
    save_Helmholtz_decomposition(U0, Phi, x, y, z, HD_dir, name)
    itr += 1


main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, endT)

