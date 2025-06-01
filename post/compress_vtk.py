import numpy as np
import os
import vtk
from tqdm import tqdm
from mod.mod_read import getGrid, getScalar, getVector, extract_number


dir = "../3D_solver/SBLI/data/3/init"


def print_data(x, y, z, rho, u, v, w, p, directory, filename):
  rho1d = np.float32(rho.flatten())
  u1d   = np.float32(u.flatten())
  v1d   = np.float32(v.flatten())
  w1d   = np.float32(w.flatten())
  p1d   = np.float32(p.flatten())

  os.makedirs(directory, exist_ok=True)
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
  writer.SetDataModeToAppended()
  writer.EncodeAppendedDataOff()
  writer.Write()


def compress_data(dir):
  files = [f for f in os.listdir(dir) if f.endswith(".vtr")]
  files.sort(key=extract_number)
  Nx, Ny, Nz, x, y, z = getGrid(os.path.join(dir, files[0]))
  for file in tqdm(files):
    file_path = os.path.join(dir, file)
    rho       = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, v, w   = getVector(file_path, Nx, Ny, Nz, 'velocity')
    p         = getScalar(file_path, Nx, Ny, Nz, 'p')
    base_name = os.path.basename(file_path)
    temp_name = base_name + "_tmp" + ".vtr"
    temp_path = os.path.join(dir, temp_name)
    print_data(x, y, z, rho, u, v, w, p, dir, temp_name)
    os.remove(file_path)
    os.rename(temp_path, file_path)


compress_data(dir)

