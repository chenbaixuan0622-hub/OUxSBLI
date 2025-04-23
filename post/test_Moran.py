import numpy as np
import matplotlib.pyplot as plt
from libpysal.weights import lat2W
from esda.moran import Moran, Moran_Local


rows, cols = 20, 20
data1 = np.zeros((rows, cols), dtype=int)
data1[:,cols//2:] = 1
data2 = np.fromfunction(lambda i, j: (i + j) % 2, (rows, cols), dtype=int)

plt.figure(figsize=(5, 4))
plt.imshow(data1, cmap="grey")
plt.show()

w = lat2W(rows, cols)
global_moran1 = Moran(data1.flatten(), w)
local_moran1  = Moran_Local(data1.flatten(), w)
local_moran1  = local_moran1.Is.reshape(rows, cols)

print(global_moran1.I)
plt.figure(figsize=(5, 4))
plt.imshow(local_moran1, cmap="RdBu_r")
plt.show()

plt.figure(figsize=(5, 4))
plt.imshow(data2, cmap="grey")
plt.show()

w = lat2W(rows, cols)
global_moran2 = Moran(data2.flatten(), w)
local_moran2  = Moran_Local(data2.flatten(), w)
local_moran2  = local_moran2.Is.reshape(rows, cols)

print(global_moran2.I)
plt.figure(figsize=(5, 4))
plt.imshow(local_moran2, cmap="RdBu_r")
plt.show()

