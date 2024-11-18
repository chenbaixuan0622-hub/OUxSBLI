import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import os
import re
from tqdm import tqdm
# AI modules
from mod_AI.utils import Dataset_classify, divide_into_batch, trainNN_classify
from mod_AI.block import noBlock
from mod_AI.cnn import Encoder
# info modules
from mod.mod_info import KMeans
# vtk modules
from mod.mod_read import getGrid, getVector, getScalar
# turb modules
from mod.mod_turb_stat import tau_2d



nx = 28
nz = 28
yp_target = 10



# make dataset
Q_dir     = "../../DNS"
Q_files   = [f for f in os.listdir(Q_dir) if f.endswith("vtr")]
num_files = len(Q_files)

Nx, Ny, Nz, x, y, z = getGrid(os.path.join(Q_dir, Q_files[0]))

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

'''
yp_path = os.path.join(Q_dir, "yp.npy")
yp      = np.load(yp_path)

for j in range(Ny):
  if yp[j] > yp_target:
    Ny_target = j
    break
'''
Ny_target = 15

n1 = nx * np.floor(Nx / nx).astype(int)
n3 = nz * np.floor(Nz / nz).astype(int)

Q   = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
u   = np.zeros((num_files,n3,n1), dtype=np.float32)
tau = np.zeros((num_files,n3,n1), dtype=np.float32)
itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_dir, Q_file)
  Q[0,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, "rho")
  Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(file_path, Nx, Ny, Nz, "velocity")
  Q[4,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, "p")
  u[itr,:,:]      = Q[1,:n3,Ny_target,:n1]
  tau[itr,:,:], _ = tau_2d(Q[:,:n3,:,:n1], x[:n1], y, z[:n3])
  itr += 1

u   = u.reshape([-1,1,nz,nx])
tau = tau.reshape([-1,nz,nx])

# make 1d data
tau1d = np.zeros(len(tau[:,0,0]), dtype=np.float32)
for i in range(len(tau[:,0,0])):
  tau1d[i] = np.mean(tau[i,int(nz/2)-1:int(nz/2)+1,int(nx/2)-1:int(nx/2)+1])

label  = KMeans(tau1d, bin=5)
images = (u - np.mean(u)) / np.std(u)

print(max(label))
print(min(label))

# 
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# plot histogram
fig, axes = plt.subplots(1,2, tight_layout=True)
axes[0].hist(tau1d, bins=100, density=True)
axes[0].set_xlabel("tau")
axes[0].set_ylabel("density")
axes[1].hist(label, bins=5, density=True)
axes[1].set_xlabel("tau")
axes[1].set_ylabel("density")
fig.savefig(os.path.join(Q_dir, "hist_tau.png"))
plt.close(fig)

'''
itr = 0
X, Z = np.meshgrid(x[:nx] * 1e3, z[:nz] * 1e3)
fig, axes = plt.subplots(2,5)
for i in range(len(label)):
  if label[i] == itr:
    cb = axes[itr//5,itr%5].contourf(X, Z, images[i,0,:,:], 100, cmap='jet')
    axes[itr//5,itr%5].set_title(str(label[i]))
    axes[itr//5,itr%5].axis('equal')
    axes[itr//5,itr%5].set_xlim(min(x[:nx] * 1e3), max(x[:nx]) * 1e3)
    axes[itr//5,itr%5].set_ylim(min(z[:nz] * 1e3), max(z[:nz]) * 1e3)
    if itr//5 == 1:
      axes[itr//5,itr%5].set_xlabel("x [mm]")
    else:
      axes[itr//5,itr%5].xaxis.set_visible('False')
      axes[itr//5,itr%5].xaxis.set_ticks([])
    if itr%5 == 0:
      axes[itr//5,itr%5].set_ylabel("z [mm]")
    else:
      axes[itr//5,itr%5].yaxis.set_visible('False')
      axes[itr//5,itr%5].yaxis.set_ticks([])
    itr += 1
  if itr == 10:
    break

cbar = fig.colorbar(cb, ax=axes.ravel().tolist(), pad=0.025)
cbar.ax.set_title("u [m/s]")
fig.savefig(os.path.join(Q_dir, "u.png"))
plt.close(fig)
'''


dataset                     = Dataset_classify(images, label)
train_dataset, test_dataset = train_test_split(dataset, test_size=0.3, shuffle=False)


# train NN
batch_size = 256

train_loader, test_loader = divide_into_batch(train_dataset, test_dataset, batch_size)

size_flatten = 3200

device = 'cuda' if torch.cuda.is_available() else 'cpu'
model = Encoder(noBlock, [64,32,16,8], size_flatten, 5).to(device)
print(model)

criterion = nn.CrossEntropyLoss()
optimizer = optim.SGD(model.parameters(), lr=0.001)

num_epochs = 500

train_loss_list, test_loss_list = trainNN_classify(model, device, optimizer, criterion, train_loader, test_loader, num_epochs)

plt.plot(range(len(train_loss_list)), train_loss_list, c='b', label='train loss')
plt.plot(range(len(test_loss_list)), test_loss_list, c='r', label='test loss')
plt.xscale("log")
plt.xlabel("epoch")
plt.ylabel("loss")
plt.legend()
plt.savefig(os.path.join(Q_dir, "loss.png"))
plt.show()

