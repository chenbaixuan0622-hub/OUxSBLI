import numpy as np
import os
import matplotlib.pyplot as plt
from pyinform import transfer_entropy, conditional_entropy
from mod.mod_info import KMeans

dir = "../../time_series"
path_data1 = os.path.join(dir, "yp269", "amp_x=0.03_k=0.004.d")
path_data2 = os.path.join(dir, "yp590", "shock_k=0.004.d")

history_len = 10
bin         = 3

data1   = np.loadtxt(path_data1)
data2   = np.loadtxt(path_data2)
binned1 = KMeans(data1[:,1], bin)
binned2 = KMeans(data2[:,1], bin)

print(binned1)
print(binned2)

TE1_2 = transfer_entropy(binned1, binned2, k=history_len, local=True)
TE2_1 = transfer_entropy(binned2, binned1, k=history_len, local=True)
H12   = conditional_entropy(binned1, binned2)
H21   = conditional_entropy(binned2, binned1)
print("TE 1_2 ", TE1_2)
print("TE 2_1 ", TE2_1)
print("TE 1_2 ", np.mean(TE1_2))
print("TE 2_1 ", np.mean(TE2_1))
print("H(1|2) ", H12)
print("H(2|1) ", H21)

