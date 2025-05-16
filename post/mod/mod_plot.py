import numpy as np
import os
import vtk
import matplotlib.pyplot as plt


def set_Params():
  plt.rcParams['font.family'] = 'Times New Roman'
  plt.rcParams['mathtext.fontset'] = 'stix'
  plt.rcParams['xtick.direction'] = 'in'
  plt.rcParams['ytick.direction'] = 'in'
  plt.rcParams['font.size'] = 12


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


def print_slice(x, y, z, rho, u, v, w, p, directory, name):
  rho1d = np.float32(rho.flatten())
  u1d   = np.float32(u.flatten())
  v1d   = np.float32(v.flatten())
  w1d   = np.float32(w.flatten())
  p1d   = np.float32(p.flatten())

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

  rho = vtk.vtkFloatArray()
  rho.SetName("rho")
  for i in range(nx * ny * nz):
    rho.InsertNextValue(rho1d[i])
  grid.GetPointData().AddArray(rho)

  velocity = vtk.vtkFloatArray()
  velocity.SetName("velocity")
  velocity.SetNumberOfComponents(3)
  for i in range(nx * ny * nz):
    velocity.InsertNextTuple3(u1d[i], v1d[i], w1d[i])
  grid.GetPointData().SetVectors(velocity)

  p = vtk.vtkFloatArray()
  p.SetName("p")
  for i in range(nx * ny * nz):
    p.InsertNextValue(p1d[i])
  grid.GetPointData().AddArray(p)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()


def print_scalar(x, y, z, p, dir, name, filename):
  p1d = np.float32(p.flatten())

  os.makedirs(dir, exist_ok=True)
  filename  = filename + ".vtr" 
  filepath  = os.path.join(dir, filename)

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

  p = vtk.vtkFloatArray()
  p.SetName(name)
  for i in range(nx * ny * nz):
    p.InsertNextValue(p1d[i])
  grid.GetPointData().AddArray(p)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()


def print_vector(x, y, z, V, dir, name, filename):
  u1d = np.float32(V[:,:,:,0].flatten())
  v1d = np.float32(V[:,:,:,1].flatten())
  w1d = np.float32(V[:,:,:,2].flatten())

  os.makedirs(dir, exist_ok=True)
  filename  = filename + ".vtr" 
  filepath  = os.path.join(dir, filename)

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

  vector = vtk.vtkFloatArray()
  vector.SetName(name)
  vector.SetNumberOfComponents(3)
  for i in range(nx * ny * nz):
    vector.InsertNextTuple3(u1d[i], v1d[i], w1d[i])
  grid.GetPointData().SetVectors(vector)
  
  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()

