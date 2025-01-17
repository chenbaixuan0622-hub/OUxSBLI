import numpy as np
from scipy.sparse import diags
import jax
import jax.numpy as jnp
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number, Vector_Data
from mod.mod_POD import make_data, make_grid
from mod.mod_nabla import Vector, Scalar, bc
from mod.mod_Poisson import Poisson
from mod.mod_Helmholtz import HD, HD_jax, HD3D, HD3D_jax, plot_Helmholtz_decomposition, save_Helmholtz_decomposition
from mod.mod_BiCGStab import BiCGStab2D, BiCGStab3D


set_Params()

# parameter
Q_dir   = "../../SBLI_1delta"
Lx1     =  0.e-3#4.e-3
Lx2     = 50.e-3#20.e-3
Ly1     = 0.e-3
Ly2     = 12.e-3
stridex = 4
stridez = 4
ny      = 100
endT    = 0.1e-3


def main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny, endT):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Data = Vector_Data(Q_dir)
  
  mean_path = os.path.join(Q_dir, 'Mean', 'Qmean.vtr')
  umean, vmean, wmean, _, _, _ = Data.getMeanVector_interp('velocity', Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny, mean_path)

  itr  = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    u, v, w, x, y, z = Data.getVector_interp(file_path, 'velocity', Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny)
    U = np.array([u - umean, v - vmean, w - wmean])
  
    nx, ny, nz = len(x), len(y), len(z)
    '''
    x0 = np.ones(ny, dtype=np.float32)
    x1 = np.ones(ny, dtype=np.float32)
    y1 = np.ones(nx, dtype=np.float32)
    '''
  
    '''
    # 2D
    U   = U[0:2,5,:,:]

    div = Vector(U, x, y).divergence()

    #phi = HD(x, y, div, x0, x1, y1)
    div = jnp.array(div)
    phi = HD_jax(len(x), len(y), x, y, div, x0, x1, y1)
    phi = np.array(phi)
    plot_Helmholtz_decomposition(U, phi, x, y)
    '''

    U0  = U.copy()
    phi = np.zeros((nz, ny, nx), dtype=np.float32)
    for itr_try in range(2):
      div = Vector(U, x, y, z).divergence(TVD=True)
      #phi = HD3D(x, y, z, div, x0, x1, y1)
    
      '''
      X, Y = np.meshgrid(x, y)
      fig, ax = plt.subplots(1, 2, figsize=(12, 6))
      c1 = ax[0].contourf(X, Y, div[1,:,:],     levels=100, cmap='jet')
      c2 = ax[1].contourf(X, Y, div_TVD[1,:,:], levels=100, cmap='jet')
      plt.show()
      '''
      
      '''
      # LU decomposition
      Div = jnp.array(div)
      jax.block_until_ready(Div)
      Phi = HD3D_jax(len(x), len(y), len(z), x, y, z, Div, x0, x1, y1)
      jax.block_until_ready(Phi)
      phi = np.array(Phi)
      '''

      # Discrete cosine transform
      PDE  = Poisson(div, x, y, z)
      px1  = np.zeros_like(div[:,:,0])
      px2  = np.zeros_like(div[:,:,0])
      dpy1 = np.zeros_like(div[:,0,:])
      dpy2 = np.zeros_like(div[:,0,:])
      phi += PDE.x_Dirichlet_y_Neumann_z_periodic(px1, px2, dpy1, dpy2)

      '''
      @jax.jit
      def set_bc(p):
        p = p.at[:,:,0].set(p[:,:,100])
        p = p.at[:,:,-1].set(0.5e0)
        p = p.at[:,0,:].set(p[:,1,:])
        p = p.at[:,-1,:].set(0.5e0)
        p = p.at[0,:,:].set(p[-2,:,:])
        p = p.at[-1,:,:].set(p[1,:,:])
        return p

      Div = jnp.array(div)
      Phi = HD3D_jax(len(x), len(y), len(z), x, y, z, Div, x0, x1, y1)
      dx  = jnp.array(-x[:-1] + x[1:])
      dy  = jnp.array(-y[:-1] + y[1:])
      dz  = jnp.array(-z[:-1] + z[1:])
      err_tol = 1.e-20
      Phi, err = BiCGStab3D(Phi, Div, dx, dy, dz, err_tol, set_bc)
      phi = np.array(Phi)
      '''

      gradPhi = Scalar(phi, x, y, z).gradient(TVD=True)
      for i in range(3):
        gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).periodic(z=True)
        gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).Neumann(x1=True, x2=True, y2=True)
        gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).Dirichlet(y1=0.e0)
      U -= gradPhi


    name = "U" + str(itr).zfill(5)
    HD_dir = os.path.join(Q_dir, "HD")
    os.makedirs(HD_dir, exist_ok=True)
    save_Helmholtz_decomposition(U0, phi, x, y, z, HD_dir, name)
    itr += 1

main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny, endT)

