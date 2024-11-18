import numpy as np
import os
import vtk
import re
from mod.mod_read import getGrid, getScalar, getVector
from mod.mod_turb_stat import tau_2d

def print_vtk(x, y, z, kp, directory, name):
  kp1d  = np.float32(kp.flatten())

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

  kp = vtk.vtkFloatArray()
  kp.SetName("kp")
  for i in range(nx * ny * nz):
    kp.InsertNextValue(kp1d[i])
  grid.GetPointData().AddArray(kp)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()


Q_directory = "../../../../../mnt/data1/TBL20240819_KEEP6thVisc4th"

Q_files = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

first_path       = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(os.path.join(first_path))

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

Q  = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
kp = np.zeros((Nz,Ny,Nx),   dtype=np.float32)

um_path = os.path.join(Q_directory, "um.npy")
vm_path = os.path.join(Q_directory, "vm.npy")
wm_path = os.path.join(Q_directory, "wm.npy")

um = np.load(um_path)
vm = np.load(vm_path)
wm = np.load(wm_path)

target_path = os.path.join(Q_directory, "Q01500.vtr")

for Q_file in Q_files:
  file_path = os.path.join(Q_directory, Q_file)
  if file_path == target_path:
    file_path = os.path.join(Q_directory, Q_file)
    Q[0,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'rho')
    Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(file_path, Nx, Ny, Nz, 'velocity')
    Q[4,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'p')
    _, ut = tau_2d(Q, x, y, z)
    print(np.mean(ut))
    kp = 0.5e0 * ((Q[1,:,:,:]-um)**2 + (Q[2,:,:,:]-vm)**2 + (Q[3,:,:,:]-wm)**2) / np.mean(ut**2)

print_vtk(np.float32(x), np.float32(y), np.float32(z), kp, Q_directory, "kp")

