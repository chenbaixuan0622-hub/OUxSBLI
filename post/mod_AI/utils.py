import numpy as np
import torch
from torch import nn
import cv2
import matplotlib.pyplot as plt


class Dataset(torch.utils.data.Dataset):
  def __init__(self, teaching_data, test_data):
    self.data    = torch.tensor(teaching_data, dtype=torch.float32, requires_grad=False)
    self.targets = torch.tensor(test_data,     dtype=torch.float32, requires_grad=False)
  
  def __len__(self):
    return len(self.data)

  def __getitem__(self, index):
    x = self.data[index]
    y = self.targets[index]
    return x, y


class Dataset_classify(torch.utils.data.Dataset):
  def __init__(self, image, label):
    self.image = torch.tensor(image, dtype=torch.float32, requires_grad=False)
    self.label = torch.tensor(label, dtype=torch.long, requires_grad=False)
  
  def __len__(self):
    return len(self.image)

  def __getitem__(self, index):
    x = self.image[index]
    y = self.label[index]
    return x, y


def divide_into_batch(train_dataset,test_dataset,batchsize):
  train_batch = torch.utils.data.DataLoader(dataset=train_dataset,
                                          batch_size=batchsize,
                                          shuffle=True)
  test_batch  = torch.utils.data.DataLoader(dataset=test_dataset,
                                          batch_size=batchsize,
                                          shuffle=True)
  return train_batch, test_batch


def trainNN(net,device,optimizer,criterion,train_batch,test_batch,epoch):
  # make lists to store loss
  train_loss_list = []
  test_loss_list  = []

  # do machine learning
  for i in torch.arange(epoch):
    # progress var
    print('---------------------------------------------')
    print("Epoch: {}/{}".format(i+1, epoch))

    # initialize loss
    train_loss = 0.e0
    test_loss  = 0.e0

    # train NN
    net.train()
    # load mini batch
    for teaching_data, test_data in train_batch:
      teaching_data = teaching_data.to(device)
      test_data     = test_data.to(device)
      optimizer.zero_grad()
      y_pred = net(test_data)
      loss   = criterion(y_pred, teaching_data)
      loss.backward()
      optimizer.step()
      train_loss += loss.item()

    # calc mean loss
    batch_train_loss = train_loss / len(train_batch)

    # evaluate NN
    net.eval()
    with torch.no_grad():
      for teaching_data, test_data in test_batch:
        teaching_data = teaching_data.to(device)
        test_data     = test_data.to(device)
        y_pred  = net(test_data)
        loss = criterion(y_pred, teaching_data)
        test_loss += loss.item()

    # calc mean loss
    batch_test_loss = test_loss / len(test_batch)

    print("Train_Loss: {:E}".format(batch_train_loss))
    print("Test_Loss : {:E}".format(batch_test_loss))

    train_loss_list.append(batch_train_loss)
    test_loss_list.append(batch_test_loss)
  return train_loss_list, test_loss_list


def trainNN_classify(net,device,optimizer,criterion,train_batch,test_batch,epoch):
  # make lists to store loss
  train_loss_list = []
  test_loss_list  = []
  train_acc_list  = []
  test_acc_list   = []

  # do machine learning
  for i in torch.arange(epoch):
    # progress var
    print('---------------------------------------------')
    print("Epoch: {}/{}".format(i+1, epoch))

    # initialize loss
    train_loss    = 0.e0
    test_loss     = 0.e0
    correct_train = 0
    correct_test  = 0

    # train NN
    net.train()
    # load mini batch
    for images, labels in train_batch:
      images = images.to(device)
      labels = labels.to(device)
      optimizer.zero_grad()
      y_pred = net(images)
      loss   = criterion(y_pred, labels)
      loss.backward()
      optimizer.step()
      train_loss += loss.item()

      _, predicted = torch.max(y_pred, 1)
      correct_train += (predicted == labels).sum().item()

    # calc mean loss
    batch_train_loss = train_loss    / len(train_batch)
    batch_train_acc  = correct_train / len(train_batch.dataset)

    # evaluate NN
    net.eval()
    with torch.no_grad():
      for images, labels in test_batch:
        images = images.to(device)
        labels = labels.to(device)
        y_pred = net(images)
        loss = criterion(y_pred, labels)
        test_loss += loss.item()

        _, predicted = torch.max(y_pred, 1)
        correct_test += (predicted == labels).sum().item()

    # calc mean loss
    batch_test_loss = test_loss    / len(test_batch)
    batch_test_acc  = correct_test / len(test_batch.dataset)

    print("Train_Loss: {:E} | Train_Acc: {:.2%}".format(batch_train_loss, batch_train_acc))
    print("Test_Loss : {:E} | Test_Acc : {:.2%}".format(batch_test_loss,  batch_test_acc))

    train_loss_list.append(batch_train_loss)
    test_loss_list.append(batch_test_loss)
    train_acc_list.append(batch_train_acc)
    test_acc_list.append(batch_test_acc)
  return train_loss_list, test_loss_list


def trainGWN(net,device,optimizer,criterion,train_batch,test_batch,epoch):
  # make lists to store loss
  train_loss_list = []
  test_loss_list  = []

  # do machine learning
  for i in torch.arange(epoch):
    # progress var
    print('---------------------------------------------')
    print("Epoch: {}/{}".format(i+1, epoch))

    # initialize loss
    train_loss = 0.e0
    test_loss  = 0.e0

    # train NN
    net.train()
    # load mini batch
    for teaching_data, test_data in train_batch:
      teaching_data = teaching_data.to(device)
      test_data     = test_data.to(device)
      optimizer.zero_grad()
      y_pred, var   = net(test_data)
      loss          = criterion(y_pred, teaching_data) / var
      loss.backward()
      optimizer.step()
      train_loss += loss.item()

    # calc mean loss
    batch_train_loss = train_loss / len(train_batch)

    # evaluate NN
    net.eval()
    with torch.no_grad():
      for teaching_data, test_data in test_batch:
        teaching_data = teaching_data.to(device)
        test_data     = test_data.to(device)
        y_pred, var   = net(test_data)
        loss = criterion(y_pred, teaching_data) / var
        test_loss += loss.item()

    # calc mean loss
    batch_test_loss = test_loss / len(test_batch)

    print("Train_Loss: {:E}".format(batch_train_loss))
    print("Test_Loss : {:E}".format(batch_test_loss))

    train_loss_list.append(batch_train_loss)
    test_loss_list.append(batch_test_loss)
  return train_loss_list, test_loss_list


def train_ESN(N_washout, N_train, N_tstart, N_test, net, X):
  '''
    arg: X[time, dim]
  '''
  N_dim   = X.shape[1] # dimension of inputs (and outputs)
  N_units = 100 * N_dim #units in the reservoir 
  
  X_washout = X[:N_washout]
  X_t = X[N_washout:N_washout+N_train-1]
  Y_t = X[N_washout+1:N_washout+N_train]
  X_test = X[N_tstart:]

  esn = net(N_units, N_dim, N_dim)
  esn.train(X_washout, X_t, Y_t)

  r0s = []
  for pic in range(10):
    fig = plt.figure(constrained_layout=True, figsize=(10, 4))

    axs = fig.subplots(1, 2)

    N_evo = 20
    offset = np.random.randint(N_test - N_evo)
    X_test_washout = X[N_tstart - N_washout + offset:N_tstart + offset]
    Y = X[N_tstart + offset : N_tstart + offset + N_evo + 1]
    X_delay = Y[:0]

    # idle iteration
    r0 = esn.open_loop(X_test_washout, np.zeros(N_units))[-1]
    r0s.append(r0)

    # the unintervened sequence x_CL
    Yh = esn.evolve(r0, N_evo)

    pairs = [[1,0],[0,1]]
    locs = ['upper left', 'lower right']

    for p in range(len(pairs)):
      i = pairs[p][0]
      j = pairs[p][1]

      targetI = 3*i
      targetJ = 3*j+1

      # the intervened sequence x_j->i 
      Yhji = esn.evolve_edge_removal(r0, targetJ, targetI, N_evo)

      ax = axs[p]

      ax.plot(Yh.T[targetI], '-', label="closed-loop", linewidth=4.0, color='blue', alpha=1, ms=10)
      ax.plot(Yhji.T[targetI], '-', label="intervened-loop", linewidth=4.0, color='red', ms=10)
      ax.plot(Y.T[targetI], label="ground truth", linewidth=2.0, color='#989A9E', linestyle='--', alpha=1)

      ax.set_xlabel('Time step')
      ax.set_ylabel(f'$x_{i+1}$')
      ax.set_xticks(ticks=[0, (Yh.shape[0]-1) / 2, Yh.shape[0]-1])
      if p == 1:
        ax.yaxis.set_label_position('right')
        ax.yaxis.set_ticks_position('right')
      if p == 0:
        ax.legend(loc=locs[p])

    plt.show()

