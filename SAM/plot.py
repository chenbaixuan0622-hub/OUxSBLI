import numpy as np
import torch
import matplotlib.pyplot as plt
import os

def plot_loss(epoch,train_loss_list,test_loss_list,dir_path):
  # save fig
  plt.figure()
  plt.title('Train and Test Loss')
  plt.xlabel('Epoch')
  plt.ylabel('Loss')
  plt.plot(range(1, epoch+1), train_loss_list, color='blue',
           linestyle='-', label='Train_Loss')
  plt.plot(range(1, epoch+1), test_loss_list,  color='red',
           linestyle='-', label='Test_Loss')
  plt.legend()
  plt.xscale('log')
  plt.yscale('log')
  plt.savefig(os.path.join(dir_path, 'loss.png'))
  # save text
  f = open(os.path.join(dir_path, 'train_loss.txt'), 'w')
  f.write(train_loss_list)
  f.close()
  f = open(os.path.join(dir_path, 'test_loss.txt'), 'w')
  f.write(test_loss_list)
  f.close()


def plot_Dataset(x,y,dataset,dir_path):
  Qplot = np.array(dataset)
  x, y  = np.meshgrid(x,y)
  fig, axes = plt.subplots(2, 1, figsize=(16,8))
  contours  = [Qplot[0,0,:,:],Qplot[1,0,:,:]]
  titles    = ['teaching','test']
  for i, ax in enumerate(axes.flat):
    CS = ax.contourf(x, y, contours[i], levels=100, cmap='turbo')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(titles[i])
  plt.tight_layout()
  plt.savefig(os.path.join(dir_path, 'dataset.png'))


def plot_result(x,y,net,device,test_batch,dir_path,input,output):
  net.eval()
  with torch.no_grad():
    for teaching_data, test_data in test_batch:
      teaching_data = teaching_data.to(device)
      test_data     = test_data.to(device)
      y_pred        = net(test_data)
    
      # plot Dataset
      test      = np.array(test_data.to('cpu'))
      teaching  = np.array(teaching_data.to('cpu'))
      pred      = np.array(y_pred.to('cpu'))
      fig, axes = plt.subplots(3, 1, figsize=(16,8))
      contours  = [test[-1,0,:,:], teaching[-1,0,:,:], pred[-1,0,:,:]]
      minval    = np.min(np.minimum(teaching[-1,0,:,:], pred[-1,0,:,:]))
      maxval    = np.max(np.maximum(teaching[-1,0,:,:], pred[-1,0,:,:]))
      mins      = [np.min(test[-1,0,:,:]), minval, minval]
      maxs      = [np.max(test[-1,0,:,:]), maxval, maxval]
      titles    = [input, output, 'predicted']
      for i, ax in enumerate(axes.flat):
        CS = ax.contourf(x, y, contours[i], levels=100, cmap='turbo', vmin=mins[i], vmax=maxs[i])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(titles[i])
        plt.tight_layout()
        plt.savefig(os.path.join(dir_path, 'trainedNN.png'))
      break

