import numpy as np
import torch
import torch.nn as nn
from torch import optim
from scipy.stats import zscore
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
from mod_AI.plot  import plot_Dataset, plot_loss, plot_result
from mod_AI.utils import Dataset, trainGWN 
from mod_AI.gwn import gwnet
from mod.mod_POD import make_data, make_grid


Q_dir   = "../3D_solver/TBL/data"
Lx1     = 28.e-3
Lx2     = 40.e-3
Ly2     = 8.e-3
stridex = 32
stridey = 16
endT    = 0.1e-3
nt      = 10
epoch   = 500
batch_size = 5

# D[space=nx*ny, time]
save_path = os.path.join(Q_dir, "D.npy")
if os.path.isfile(save_path):
  D = np.load(save_path)
  _, _, _, _, _, x, y, t = make_grid(Lx1, Lx2, Ly2, stridex, stridey, endT, Q_dir)
else:
  x, y, t, D = make_data(Lx1, Lx2, Ly2, stridex, stridey, endT, Q_dir)
  np.save(save_path, D)

nx = len(x)
ny = len(y)
Nt = len(t)
print("nx ", nx, " ny ", ny, " Nt ", Nt)

mean = np.mean(D)
std  = np.std(D)
D    = zscore(D, axis=None)

# calc initial matrix
sigma2 = 0.1e0 * np.sqrt((-x[0] + x[1])**2 + (-y[0] + y[1])**2)
eps    = 0.5e0
Ainit  = np.zeros((nx*ny, nx*ny), dtype=np.float32)
X, Y   = np.meshgrid(x, y)
X      = X.flatten()
Y      = Y.flatten()
for j in range(nx*ny):
  for i in range(nx*ny):
    d = np.sqrt((X[i] - X[j])**2 + (Y[i] - Y[j])**2)
    if d > 0.e0 and np.exp(-d**2 / sigma2) > eps:
      Ainit[j,i] = np.exp(-d**2 / sigma2)
    else:
      Ainit[j,i] = 0.e0

save_path = os.path.join(Q_dir, "GWN_initial_matrix")
np.save(save_path, Ainit)

# set GPU and NN property
device    = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
criterion = nn.L1Loss()
Ainit     = torch.tensor(Ainit, device=device)
net       = gwnet(device, in_dim=1, out_dim=1, residual_channels=32, \
                  dilation_channels=32, skip_channels=32, \
                  end_channels=32, dropout=0.1, nx=nx*ny, nt=nt, Ainit=Ainit)
net       = net.to(device)
optimizer = optim.RAdam(net.parameters(), lr=0.001)

# teaching_data[data_len, feature=1 (u), space=nx*ny, time]
#     test_data[data_len, feature=1 (u), space=nx*ny, time]
data_len = Nt - nt
test_data     = np.zeros((data_len,1,nx*ny,nt), dtype=np.float32)
teaching_data = np.zeros((data_len,1,nx*ny,1),  dtype=np.float32)
for i in range(data_len):
  start = i
  test_data[i,0,:,:]     = D[:,start:start+nt]
  teaching_data[i,0,:,0] = D[:,start+nt]

print("teaching data size: ", np.shape(teaching_data), " test data size: ", np.shape(test_data))

dataset                     = Dataset(teaching_data, test_data)
train_dataset, test_dataset = train_test_split(dataset, test_size=0.2, shuffle=False)

train_batch = torch.utils.data.DataLoader(dataset=train_dataset,
                                          batch_size=batch_size,
                                          shuffle=False)
test_batch  = torch.utils.data.DataLoader(dataset=test_dataset,
                                          batch_size=batch_size,
                                          shuffle=False)

print("train batch size: ", len(train_batch), " test batch size: ", len(test_batch))

train_loss_list, test_loss_list = trainGWN(net,device,optimizer,criterion,train_batch,test_batch,epoch)
plot_loss(epoch-1,train_loss_list[1:],test_loss_list[1:],Q_dir,0)

# plot result
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12

net.eval()
loss0 = 1e12
with torch.no_grad():
  for itr, tuple in enumerate(test_batch):
    teaching_data, test_data = tuple
    teaching_data = teaching_data.to(device)
    test_data     = test_data.to(device)
    y_pred, var   = net(test_data)
    loss          = criterion(y_pred, teaching_data)

    if loss < loss0:
      loss0 = loss
      best  = itr

  teaching_data, test_data = list(test_batch)[best]
  test_data                = test_data.to(device)
  pred, var                = net(test_data)
  
  teaching = np.array(teaching_data)
  pred     = np.array(pred.to('cpu'))
  # unnormalize
  teaching = teaching * std + mean
  pred     =     pred * std + mean
  contours = [teaching[-1,0,:,:], pred[-1,0,:,:]]
  x, y     = np.meshgrid(x*1e3, y*1e3)
  crange   = np.linspace(0, 500, 50)
  fig, ax  = plt.subplots(1, 2, figsize=(12, 6))
  im1 = ax[0].contourf(x, y, contours[0].reshape(ny,nx), crange, cmap='jet', extend='both')
  im2 = ax[1].contourf(x, y, contours[1].reshape(ny,nx), crange, cmap='jet', extend='both')
  for a in ax:
    a.set_xlabel("x")
    a.set_ylabel("y")
    a.set_aspect('equal', adjustable='box')

  cbar = plt.colorbar(im1, ax=ax, extendrect=True, orientation='horizontal', \
                      pad=0.1, fraction=0.046, location='top')
  save_path = os.path.join(Q_dir, "GWN.png")
  plt.savefig(save_path)
  plt.close()

  # check causal matrix
  A = net.A.detach().cpu().numpy()
  save_path = os.path.join(Q_dir, "GWN_causal_matrix")
  np.save(save_path, A)

