import numpy as np
import torch
import matplotlib.pyplot as plt
import os

def plot_loss(epoch,train_loss_list,test_loss_list,dir_path,itr):
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
  filename1 = 'loss' + str(itr) + '.png'
  plt.savefig(os.path.join(dir_path, filename1))
  # save text
  filename2 = 'train_loss' + str(itr) + '.txt'
  file_path = os.path.join(dir_path, filename2)
  with open(file_path, "w", encoding="UTF-8") as fo:
    for i in range(len(train_loss_list)):
      print(f'{train_loss_list[i]}:.3e', file=fo)

  filename3 = 'test_loss' + str(itr) + '.txt'
  file_path = os.path.join(dir_path, filename3)
  with open(file_path, "w", encoding="UTF-8") as fo:
    for i in range(len(train_loss_list)):
      print(f'{test_loss_list[i]}:.3e', file=fo)


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
    cbar = plt.colorbar(CS, ax=ax)
    cbar.set_label('K')
  plt.tight_layout()
  plt.savefig(os.path.join(dir_path, 'dataset.png'))


def plot_result(x,y,net,device,test_batch,dir_path,input,output,itr,teaching_mean,test_mean,teaching_std,test_std):
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
      # unnormalize data
      test      = test * test_std         + test_mean
      teaching  = teaching * teaching_std + teaching_mean
      pred      = pred * teaching_std     + teaching_mean
      fig, axes = plt.subplots(1, 3, figsize=(24,8))
      contours  = [test[-1,0,:,:], teaching[-1,0,:,:], pred[-1,0,:,:]]
      minval    = np.min(np.minimum(teaching[-1,0,:,:], pred[-1,0,:,:]))
      maxval    = np.max(np.maximum(teaching[-1,0,:,:], pred[-1,0,:,:]))
      mins      = [np.min(test[-1,0,:,:]), minval, minval]
      maxs      = [np.max(test[-1,0,:,:]), maxval, maxval]
      titles    = [input, output, 'predicted']
      for i, ax in enumerate(axes.flat):
        CS = ax.contourf(x, y, contours[i], levels=10, cmap='turbo', vmin=mins[i], vmax=maxs[i])
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_title(titles[i])

        cbar = plt.colorbar(CS, ax=ax)
        cbar.set_label('K')

        plt.tight_layout()
        filename = 'trainedNN' + str(itr) + '.png'
        plt.savefig(os.path.join(dir_path, filename))
      break

