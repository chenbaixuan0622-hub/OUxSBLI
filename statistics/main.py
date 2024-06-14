import os
import time
import pickle
import numpy as np
from scipy.stats import gaussian_kde
from scipy.fft import fft
import matplotlib.pyplot as plt
import readVTK
import turbulent_statistics as ts
import myplot as myplt

start_time = time.time()

# search vtk files
directory_path = "./data/"
vtk_files = [f for f in os.listdir(directory_path) if f.endswith(".vtk")]
num_files = len(vtk_files)

# get grid information
first_path = os.path.join(directory_path, vtk_files[0])
Nx, Ny, Nz, x, y, z = readVTK.getGrid(first_path)

# check save data
means = ["./data/umean.npy", "./data/vmean.npy", "./data/wmean.npy"]
if all(os.path.exists(mean) for mean in means):
  # read pickle
  umean = np.load("./data/umean.npy")
  vmean = np.load("./data/vmean.npy")
  wmean = np.load("./data/wmean.npy")
else:
  # get mean velocity
  umean, vmean, wmean = readVTK.getMean(directory_path,vtk_files,Nx,Ny,Nz)
  # save mean velocity
  np.save("./data/umean.npy", umean)
  np.save("./data/vmean.npy", vmean)
  np.save("./data/wmean.npy", wmean)
'''
myplt.plot_velocity(x,y,umean[0,:,:])
myplt.plot_velocity(x,y,vmean[0,:,:])
myplt.plot_velocity(x,y,wmean[0,:,:])
'''
'''
# calc 2 points longitudial / lateral velocity correlation
u1 = np.zeros(num_files)
v1 = np.zeros(num_files)
u2 = np.zeros(num_files)
v2 = np.zeros(num_files)

# choose grid points
Nx1 = 50
Ny1 = 20
Nx2 = 150
Ny2 = 20

t = 0
for vtk_file in vtk_files:
  # get path
  file_path = os.path.join(directory_path, vtk_file)
  u, v, w = readVTK.getVelocity(file_path,Nx,Ny,Nz)
  u1[t] = u[0,Ny1,Nx1] - umean[0,Ny1,Nx1]
  v1[t] = v[0,Ny1,Nx1] - vmean[0,Ny1,Nx1]
  u2[t] = u[0,Ny2,Nx2] - umean[0,Ny2,Nx2]
  v2[t] = v[0,Ny2,Nx2] - vmean[0,Ny2,Nx2]
  t += 1
  print("read No.", vtk_file, "file")

# velocity correlation
R11 = ts.longitudinal_corr(x[Nx1],y[Ny1],u1,v1,x[Nx2],y[Ny2],u2,v2)
R22 = ts.lateral_corr(x[Nx1],y[Ny1],u1,v1,x[Nx2],y[Ny2],u2,v2)

print("longitudinal velocity correlation",R11)
print("lateral velocity correlation",R22)
'''
'''
del u1
del v1
del u2
del v2
del R11
del R22
'''

# calc integral length scale
# reference point
Nx1 = 50
Ny1 = 20
Nz1 = 0

length = 100
span = Nz-1
us = np.zeros((num_files,span,length))
vs = np.zeros((num_files,span,length))

t = 0
for vtk_file in vtk_files:
  # get path
  file_path = os.path.join(directory_path, vtk_file)
  u, v, w = readVTK.getVelocity(file_path,Nx,Ny,Nz)
  for j in range(span):
    for i in range(length):
      I = Nx1 + i
      J = Nz1 + j
      us[t,j,i] = u[J,Ny1,I] - umean[J,Ny1,I]
      vs[t,j,i] = v[J,Ny1,I] - vmean[J,Ny1,I]
  t += 1
  print("read No.", vtk_file, "file")

R11, R22 = ts.longitudinal_and_lateral_corr(span,length,Nx1,Ny1,x,y,us,vs)
L11, L22 = ts.integral_scale(x[Nx1:Nx1+length-1],R11,R22)
print("longitudinal integral scale is",L11)
print("lateral integral scale is",L22)
myplt.plot_corr(x,Nx1,length,R11,R22)

# calculate probability density function
x_range = np.linspace(min(us[:,0,0]), max(us[:,0,0]), 1000)
u_kde = gaussian_kde(us[:,0,0])
u_pdf = u_kde(x_range)
myplt.plot_pdf(x_range,u_pdf)

end_time = time.time()

# write output file
with open("./data/output.txt", "w", encoding="UTF-8") as f:
  print("grid information",file=f)
  print("Nx =",Nx,file=f)
  print("Ny =",Ny,file=f)
  print("Nz =",Nz,file=f)
  print("longitudinal integral scale is",L11,file=f)
  print("lateral integral scale is",L22,file=f)
  print("elapsed time:", end_time - start_time,file=f)

