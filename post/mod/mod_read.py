import numpy as np
import os
import re
import vtk
from vtk.util import numpy_support
from tqdm import tqdm


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


def gridInfo(data_directory):
  file_path = os.path.join(data_directory, "x.npy")
  x         = np.load(file_path)
  file_path = os.path.join(data_directory, "y.npy")
  y         = np.load(file_path)
  file_path = os.path.join(data_directory, "z.npy")
  z         = np.load(file_path)
  Nx        = len(x)
  Ny        = len(y)
  Nz        = len(z)
  return np.float32(x), np.float32(y), np.float32(z), Nx, Ny, Nz


def getGrid(file_path):
  # make VTK Structured Grid Reader
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()
  
  # get grid
  grid = reader.GetOutput()
  x = numpy_support.vtk_to_numpy(grid.GetXCoordinates())
  y = numpy_support.vtk_to_numpy(grid.GetYCoordinates())
  z = numpy_support.vtk_to_numpy(grid.GetZCoordinates())
  return len(x), len(y), len(z), np.float32(x), np.float32(y), np.float32(z)


def getVector(file_path,Nx,Ny,Nz,name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()

  # get dataset
  Q = reader.GetOutput()
  V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))

  u = np.reshape(V[:,0], [Nz,Ny,Nx])
  v = np.reshape(V[:,1], [Nz,Ny,Nx])
  w = np.reshape(V[:,2], [Nz,Ny,Nx])
  return np.float32(u), np.float32(v), np.float32(w)


def getScalar(file_path,Nx,Ny,Nz,name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()

  # get dataset
  Q = reader.GetOutput()
  a = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  a = np.reshape(a, [Nz,Ny,Nx])
  return np.float32(a)


def getMeanVector(directory_path,vtk_files,Nx,Ny,Nz,name):
  um = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vm = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  wm = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  for vtk_file in tqdm(vtk_files):
    file_path = os.path.join(directory_path, vtk_file)
    u, v, w   = getVector(file_path,Nx,Ny,Nz,name)
    um += u
    vm += v
    wm += w
  n  = float(len(vtk_files))
  um = um / n
  vm = vm / n
  wm = wm / n
  return um, vm, wm


def getMeanScalar(directory_path,vtk_files,Nx,Ny,Nz,name):
  am = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  for vtk_file in tqdm(vtk_files):
    file_path = os.path.join(directory_path, vtk_file)
    a  = getScalar(file_path,Nx,Ny,Nz,name)
    am += a
  n  = float(len(vtk_files))
  am = am / n
  return am

