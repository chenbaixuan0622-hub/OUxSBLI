import numpy as np
import os
import vtk
from mod.mod_read import getGrid

def print_vtk(x, y, z, rho, u, v, w, p, directory, name):
  rho1d  = np.float32(rho.flatten())
  u1d    = np.float32(u.flatten())
  v1d    = np.float32(v.flatten())
  w1d    = np.float32(w.flatten())
  p1d    = np.float32(p.flatten())

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


#Q_directory = "../3D_solver/TBL/data"
Q_directory = "../../SBLI/SBLI_4delta/stat03ms_04ms_SLAU"

Q_files = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

first_path       = os.path.join(Q_directory, Q_files[0])
_, _, _, x, y, z = getGrid(os.path.join(first_path))

rho_path  = os.path.join(Q_directory, "rho.npy")
u_path    = os.path.join(Q_directory, "u.npy"  )
v_path    = os.path.join(Q_directory, "v.npy"  )
w_path    = os.path.join(Q_directory, "w.npy"  )
p_path    = os.path.join(Q_directory, "p.npy"  )

rho  = np.load(rho_path)
u    = np.load(u_path  )
v    = np.load(v_path  )
w    = np.load(w_path  )
p    = np.load(p_path  )

print_vtk(np.float32(x), np.float32(y), np.float32(z), rho, u, v, w, p, Q_directory, "Qmean")

