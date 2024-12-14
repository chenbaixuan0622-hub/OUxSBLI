import numpy as np
import os
import matplotlib.pyplot as plt
from mod.mod_info import KMeans, causal_map

Q_dir = ""
y_dir = ""
lx    = []
lz    = []

history_len = 10
bin         = 5

Nk = len(Q_path)
Nt = 
V  = np.zeros((Nk,Nt), dtype=np.float32)
for i in range(Nk):
  name = "amp_" + "x=" + str(lx[i]) + "_k=" + str(lz[i])
  path   = os.path.join(Q_dir, y_dir, name)
  amp    = np.loadtext(Q_path[i])
  binned = KMeans(amp[], bin)
  V[i,:] = binned

map = causal_map(V, history_len, bin)
save_path = os.path.join(Q_dir, "map.npy")
np.save(save_path, map)
print("causal map")
print(map)

