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

