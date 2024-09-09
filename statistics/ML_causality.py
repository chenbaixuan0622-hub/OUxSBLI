import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import optim
from pytorch_memlab import profile
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import os
import sys
import cv2
import matplotlib.pyplot as plt
# AI modules
from mod_AI.plot  import plot_Dataset, plot_loss, plot_result
from mod_AI.utils import Dataset, divide_into_batch, trainNN 
from mod_AI.block import noBlock, BasicBlock, DenseBlock, SpectralTransform, FourierBlock
from mod_AI.cnn   import ResNet, DenseNet, UNet
# vtk modules
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import tau_2d


Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240819_KEEP6thVisc4th"
Q_files     = [f for f in os.listdir(Q_directory) if f.endswith("vtr")]
num_files   = len(Q_files)

# read grid information
first_path = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
yp_path = os.path.join(Q_directory, "yp.npy")
yp      = np.load(yp_path)

yp1 = 5
yp2 = 100

# calc yp
for j in range(Ny):
  if yp[j] > yp1:
    Ny_target1 = j
    break

for j in range(Ny):
  if yp[j] > yp2:
    Ny_target2 = j
    break

print("yp =", yp1, " when j = ", Ny_target1)
print("yp =", yp2, " when j = ", Ny_target2)

nx = 32
nz = 32

# set GPU and NN property
device       = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
net_instance = UNet(FourierBlock, [64,32,16,8]).to(device)
net          = torch.jit.script(net_instance).to(device)
criterion    = nn.MSELoss()
optimizer    = optim.Adam(net.parameters())
epoch        = 50
batchsize    = 100
print(net)
print("Device: {}".format(device))

num_split = 500

um_path = os.path.join(Q_directory, "um.npy"  )
vm_path = os.path.join(Q_directory, "vm.npy"  )
wm_path = os.path.join(Q_directory, "wm.npy"  )
um = np.load(um_path)
vm = np.load(vm_path)
wm = np.load(wm_path)
kp = np.zeros((num_files // num_split,Nz,2,Nx), dtype=np.float32)
teaching_data = np.zeros((num_files // num_split * Nz // nz * Nx // nx, 1, nz, nx), dtype=np.float32)
test_data     = np.zeros((num_files // num_split * Nz // nz * Nx // nx, 1, nz, nx), dtype=np.float32)

for i in range(num_split):
  '''
  # load NN parameters
  pdir_path = os.path.join("./data/", str(i-1))
  if os.path.isfile(os.path.join(pdir_path, "model.pth")):
    net.load_state_dict(torch.load(os.path.join(pdir_path, "model.pth")))
    print("load previous parameters")
  '''

  start = int(num_files / num_split) * i
  end   = int(num_files / num_split) * (i + 1)
  print("read from ", start, " to ", end)

  itr = 0
  print("reading VTK files")
  Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  for Q_file in tqdm(Q_files[start:end]):
    file_path = os.path.join(Q_directory, Q_file)
    Q[0,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, "rho")
    Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(file_path, Nx, Ny, Nz, "velocity")
    Q[4,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz,"p")
    tau, _ = tau_2d(Q, x, y, z)

    kp[itr,:,0,:] = 0.5e0 * ((Q[1,:,Ny_target1,:] - um[:,Ny_target1,:])**2 \
                           + (Q[2,:,Ny_target1,:] - vm[:,Ny_target1,:])**2 \
                           + (Q[3,:,Ny_target1,:] - wm[:,Ny_target1,:])**2) / tau**2
    kp[itr,:,1,:] = 0.5e0 * ((Q[1,:,Ny_target2,:] - um[:,Ny_target2,:])**2 \
                           + (Q[2,:,Ny_target2,:] - vm[:,Ny_target2,:])**2 \
                           + (Q[3,:,Ny_target2,:] - wm[:,Ny_target2,:])**2) / tau**2

    itr += 1

  itr = 0
  for l in range(num_files // num_split):
    for k in range(0, Nz // nz, nz):
      for i in range(0, Nx // nx, nx):
        teaching_data[itr,0,:,:] = kp[l,k:k+nz,0,i:i+nx]
        test_data[itr,0,:,:]     = kp[l,k:k+nz,1,i:i+nx]
        itr += 1

  del Q
  
  '''
  # dataset is torch tensor
  dataset                     = Dataset(teaching_data, test_data)
  train_dataset, test_dataset = train_test_split(dataset, test_size=0.3, shuffle=False)

  # make dataloader
  plot_Dataset(x,y,train_dataset[50],dir_path)
  train_batch, test_batch = divide_into_batch(train_dataset,test_dataset,batchsize)

  net, train_loss_list, test_loss_list = trainNN(net,device,optimizer,criterion,train_batch,test_batch,epoch)
  plot_loss(epoch,train_loss_list,test_loss_list,dir_path)
  torch.save(net.state_dict(), os.path.join(dir_path, 'model.pth'))

  plot_result(x,y,net,device,test_batch,dir_path,'p','Q')
  '''

