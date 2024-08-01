import numpy as np
import torch
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


def divide_into_batch(train_dataset,test_dataset,batchsize):
  train_batch = torch.utils.data.DataLoader(dataset=train_dataset,
                                          batch_size=batchsize,
                                          shuffle=True,
                                          num_workers=2)
  test_batch  = torch.utils.data.DataLoader(dataset=test_dataset,
                                          batch_size=batchsize,
                                          shuffle=True,
                                          num_workers=2)
  return train_batch, test_batch


def trainNN(net,device,optimizer,criterion,train_batch,test_batch,epoch):
  # make lists to store MSE
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
      # transfer Tensor to GPU
      teaching_data = teaching_data.to(device)
      test_data     = test_data.to(device)
      # initialize grad
      optimizer.zero_grad()
      # calc pred
      y_pred = net(test_data)
      # calc loss
      loss   = criterion(y_pred, teaching_data)
      # calc grad
      loss.backward()
      # update parameters
      optimizer.step()
      # stock train loss
      train_loss += loss.item()

    # calc mean loss
    batch_train_loss = train_loss / len(train_batch)

    # evaluate NN
    net.eval()
    with torch.no_grad():
      for teaching_data, test_data in test_batch:
        # transfer Tensor to GPU
        teaching_data = teaching_data.to(device)
        test_data     = test_data.to(device)
        # calc pred
        y_pred = net(test_data)
        # calc loss
        loss   = criterion(y_pred, teaching_data)
        # stock test loss
        test_loss += loss.item()

    # calc mean loss
    batch_test_loss = test_loss / len(test_batch)

    print("Train_Loss: {:E}".format(batch_train_loss))
    print("Test_Loss : {:E}".format(batch_test_loss))

    train_loss_list.append(batch_train_loss)
    test_loss_list.append(batch_test_loss)
  return net, train_loss_list, test_loss_list

