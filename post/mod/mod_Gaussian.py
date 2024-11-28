import numpy as np
from scipy.ndimage import gaussian_filter
import vtk
import os
from mod.mod_read import getGrid, getVector


def Gaussian(nx1,nx2,ny1,ny2,nz1,nz2,Q_path,sigma):
  Nx, Ny, Nz, X, Y, Z = getGrid(Q_path)
  U, V, W = getVector(Q_path, Nx, Ny, Nz, "velocity")
  
  u = U[nz1:nz2,ny1:ny2,nx1:nx2]
  v = V[nz1:nz2,ny1:ny2,nx1:nx2]
  w = W[nz1:nz2,ny1:ny2,nx1:nx2]
  x = X[nx1:nx2]
  y = Y[ny1:ny2]
  z = Z[nz1:nz2]

  ug = gaussian_filter(u, sigma=[sigma,sigma,sigma], mode='nearest')
  vg = gaussian_filter(v, sigma=[sigma,sigma,sigma], mode='nearest')
  wg = gaussian_filter(w, sigma=[sigma,sigma,sigma], mode='nearest')
  
  basename = os.path.splitext(os.path.basename(Q_path))[0]
  name = basename + '_sigma' + str(int(sigma))
  dir = os.path.dirname(Q_path)
  print_Gaussian(x, y, z, ug, vg, wg, dir, name)


def print_Gaussian(x, y, z, u, v, w, directory, name):
  u1d = np.float32(u.flatten())
  v1d = np.float32(v.flatten())
  w1d = np.float32(w.flatten())

  os.makedirs(directory, exist_ok=True)
  filename  = name + ".vtr" 
  filepath  = os.path.join(directory, filename)

  x_coords = vtk.vtkFloatArray()
  y_coords = vtk.vtkFloatArray()
  z_coords = vtk.vtkFloatArray()
  x_coords.SetName("X-Axis")
  y_coords.SetName("Y-Axis")
  z_coords.SetName("Z-Axis")

  nx = len(x)
  ny = len(y)
  nz = len(z)
  
  for i in range(nx):
    x_coords.InsertNextValue(x[i])
  for j in range(ny):
    y_coords.InsertNextValue(y[j])
  for k in range(nz):
    z_coords.InsertNextValue(z[k])
  
  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions(nx, ny, nz)
  grid.SetXCoordinates(x_coords)
  grid.SetYCoordinates(y_coords)
  grid.SetZCoordinates(z_coords)

  velocity = vtk.vtkFloatArray()
  velocity.SetName("velocity")
  velocity.SetNumberOfComponents(3)
  for i in range(nx * ny * nz):
    velocity.InsertNextTuple3(u1d[i], v1d[i], w1d[i])
  grid.GetPointData().SetVectors(velocity)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()

