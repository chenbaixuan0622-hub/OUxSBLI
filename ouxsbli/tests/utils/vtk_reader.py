import numpy as np
import os
import re
import vtk
from vtk.util import numpy_support


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vt[rs]$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


def get_ext(file_path):
  dirname, basename = os.path.split(file_path)
  basename_without_ext, ext = basename.split('.', 1)
  return ext


def getGrid_Rect(file_path):
  # make VTK Rectilinear Grid Reader
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()
  # get grid
  grid = reader.GetOutput()
  x = numpy_support.vtk_to_numpy(grid.GetXCoordinates())
  y = numpy_support.vtk_to_numpy(grid.GetYCoordinates())
  z = numpy_support.vtk_to_numpy(grid.GetZCoordinates())
  return len(x), len(y), len(z), x, y, z


def getGrid_Str(file_path):
  # make VTK Structured Grid Reader
  reader = vtk.vtkXMLStructuredGridReader()
  reader.SetFileName(file_path)
  reader.Update()
  # get grid
  grid = reader.GetOutput()
  dims = [0, 0, 0]
  grid.GetDimensions(dims)
  Nx, Ny, Nz = dims
  points = numpy_support.vtk_to_numpy(grid.GetPoints().GetData())
  points = points.reshape((Nz, Ny, Nx, 3))
  return Nx, Ny, Nz, points[:,:,:,0], points[:,:,:,1], points[:,:,:,2]


def getGrid(file_path):
  ext = get_ext(file_path)
  if ext == 'vtr':
    Nx, Ny, Nz, x, y, z = getGrid_Rect(file_path)
  elif ext == 'vts':
    Nx, Ny, Nz, x, y, z = getGrid_Str(file_path)
  else:
    raise ValueError("Invalid file type:", file_path)
  return Nx, Ny, Nz, x, y, z


def getVector(file_path, Nx, Ny, Nz, name):
  ext = get_ext(file_path)
  if ext == 'vtr':
    reader = vtk.vtkXMLRectilinearGridReader()
  elif ext == 'vts':
    reader = vtk.vtkXMLStructuredGridReader()
  else:
    raise ValueError("Invalid file type:", file_path)
  reader.SetFileName(file_path)
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  # get dataset
  Q = reader.GetOutput()
  V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  V = V.reshape((Nz,Ny,Nx,3))
  u = V[:,:,:,0]
  v = V[:,:,:,1]
  w = V[:,:,:,2]
  return u, v, w


def getScalar(file_path, Nx, Ny, Nz, name):
  ext = get_ext(file_path)
  if ext == 'vtr':
    reader = vtk.vtkXMLRectilinearGridReader()
  elif ext == 'vts':
    reader = vtk.vtkXMLStructuredGridReader()
  else:
    raise ValueError("Invalid file type:", file_path)
  reader.SetFileName(file_path)
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  # get dataset
  Q = reader.GetOutput()
  a = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  a = a.reshape((Nz,Ny,Nx))
  return a


def getQ(file_path, Nx, Ny, Nz, reader=None):
  if reader is None:
    ext = get_ext(file_path)
    if ext == 'vtr':
      reader = vtk.vtkXMLRectilinearGridReader()
    elif ext == 'vts':
      reader = vtk.vtkXMLStructuredGridReader()
    else:
      raise ValueError("Invalid file type:", file_path)
    reader.GetPointDataArraySelection().DisableAllArrays()
    reader.GetPointDataArraySelection().EnableArray("rho")
    reader.GetPointDataArraySelection().EnableArray("velocity")
    reader.GetPointDataArraySelection().EnableArray("p")
  reader.SetFileName(file_path)
  reader.Update()
  # get dataset
  Q   = reader.GetOutput()
  rho = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("rho"))
  rho = rho.reshape((Nz,Ny,Nx))
  V   = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("velocity"))
  V   = V.reshape((Nz,Ny,Nx,3))
  u   = V[:,:,:,0]
  v   = V[:,:,:,1]
  w   = V[:,:,:,2]
  p   = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("p"))
  p   = p.reshape((Nz,Ny,Nx))
  return rho, u, v, w, p


def initial_vtr(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vt[rs]"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr or Q*.vts files found in {data_dir}")
  return min(files, key=lambda f: extract_number(os.path.basename(f)))


def latest_vtr(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vt[rs]"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr or Q*.vts files found in {data_dir}")
  return max(files, key=lambda f: extract_number(os.path.basename(f)))

