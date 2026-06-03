import numpy as np
import os
import re
import vtk
from vtk.util import numpy_support


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


def extract_number_vts(filename):
  match = re.search(r'Q(\d+)\.vts$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


def getGrid(file_path):
  if str(file_path).endswith('.vts'):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(str(file_path))
    reader.Update()
    grid = reader.GetOutput()
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    ni, nj, nk = dims
    x = np.arange(ni, dtype=np.float64)
    y = np.arange(nj, dtype=np.float64)
    z = np.arange(nk, dtype=np.float64)
    return ni, nj, nk, x, y, z
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(str(file_path))
  reader.Update()
  grid = reader.GetOutput()
  x = numpy_support.vtk_to_numpy(grid.GetXCoordinates())
  y = numpy_support.vtk_to_numpy(grid.GetYCoordinates())
  z = numpy_support.vtk_to_numpy(grid.GetZCoordinates())
  return len(x), len(y), len(z), x, y, z


def getVector(file_path, Nx, Ny, Nz, name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(str(file_path))
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  Q = reader.GetOutput()
  V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  V = V.reshape((Nz,Ny,Nx,3))
  u = V[:,:,:,0]
  v = V[:,:,:,1]
  w = V[:,:,:,2]
  return u, v, w


def getScalar(file_path, Nx, Ny, Nz, name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(str(file_path))
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  Q = reader.GetOutput()
  a = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  a = a.reshape((Nz,Ny,Nx))
  return a


def getQ(file_path, Nx, Ny, Nz, reader=None):
  if str(file_path).endswith('.vts'):
    r = vtk.vtkXMLStructuredGridReader()
    r.SetFileName(str(file_path))
    r.Update()
    Q = r.GetOutput()
    rho = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("rho"))
    rho = rho.reshape((Nz, Ny, Nx))
    V   = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("velocity"))
    V   = V.reshape((Nz, Ny, Nx, 3))
    u   = V[:,:,:,0]
    v   = V[:,:,:,1]
    w   = V[:,:,:,2]
    p   = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("p"))
    p   = p.reshape((Nz, Ny, Nx))
    return rho, u, v, w, p
  if reader is None:
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.GetPointDataArraySelection().DisableAllArrays()
    reader.GetPointDataArraySelection().EnableArray("rho")
    reader.GetPointDataArraySelection().EnableArray("velocity")
    reader.GetPointDataArraySelection().EnableArray("p")
  reader.SetFileName(str(file_path))
  reader.Update()
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
  files = glob.glob(os.path.join(str(data_dir), "Q*.vtr"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr files found in {data_dir}")
  return min(files, key=lambda f: extract_number(os.path.basename(f)))


def latest_vtr(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vtr"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr files found in {data_dir}")
  return max(files, key=lambda f: extract_number(os.path.basename(f)))


def latest_vts(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vts"))
  if not files:
    raise FileNotFoundError(f"No Q*.vts files found in {data_dir}")
  return max(files, key=lambda f: extract_number_vts(os.path.basename(f)))
