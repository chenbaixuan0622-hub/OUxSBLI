import numpy as np
import os
import matplotlib.pyplot as plt

def plot_causality(uu,title):
  x = np.linspace(1, 3, len(uu[:,0]))
  y = np.linspace(1, 3, len(uu[0,:]))
  x, y = np.meshgrid(x,y)
  plt.contourf(x, y, np.transpose(uu), cmap=plt.cm.jet, levels=100)
  plt.colorbar()
  plt.savefig(title)
  plt.close()

def plot_hist(x,px,filename):
  X = []
  for i in range(1, len(x)):
    X.append(0.5e0 * (x[i-1] + x[i]))
  plt.bar(X, px)
  plt.savefig(filename)
  plt.close()

def plot_corr(x,Nx1,length,R11,R22):
  plt.plot(figsize=(8,6))
  plt.plot(x[Nx1:Nx1+length-1],R11)
  plt.xlabel("r")
  plt.ylabel("f(r)")
  plt.savefig("data/longitudinal_corr.png")
  plt.close()
  plt.plot(figsize=(8,6))
  plt.plot(x[Nx1:Nx1+length-1],R22)
  plt.xlabel("r")
  plt.ylabel("g(r)")
  plt.savefig("data/lateral_corr.png")
  plt.close()

def plot_scalar(x,y,u,title):
  plt.axis("equal")
  plt.axis("off")
  plt.contourf(x, y, u, cmap=plt.cm.jet, levels=100)
  plt.colorbar()
  plt.savefig(title)
  plt.close()

'''
#def plot_fluxtuating_velocity(directory_path,vtk_files,Nx,Ny,Nz,x,y,z,umean,vmean,wmean): 
  i = 0
  for vtk_file in vtk_files:
    file_path = os.path.join(directory_path, vtk_file)
    u, v, w = readVTK.getVelocity(file_path,Nx,Ny,Nz)
    plt.axis("equal")
    plt.axis("off")
    plt.contourf(x, y, (u[0,:,:] - umean[0,:,:]), cmap=plt.cm.jet, levels=100)
    plt.colorbar()
    file_path = os.path.join("./img",u,str(i))
    plt.savefig(file_path)
    plt.contourf(x, y, (v[0,:,:] - vmean[0,:,:]), cmap=plt.cm.jet, levels=100)
    plt.colorbar()
    file_path = os.path.join("./img",v,str(i))
    plt.savefig(file_path)
    i += 1
'''

