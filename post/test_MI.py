import numpy as np
from sklearn.feature_selection import mutual_info_regression
import matplotlib.pyplot as plt
from mod.mod_plot import set_Params 

set_Params()

n = 1000

theta = np.random.randn(n)
x = np.cos(2.e0 * np.pi * theta)
y = np.sin(2.e0 * np.pi * theta)

plt.scatter(x,y)
plt.xlabel('$\it{x}$')
plt.ylabel('$\it{y}$')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
plt.close()

corr = np.corrcoef(x,y)
mi   = np.mean(mutual_info_regression(np.reshape(x, [-1,1]),y))

print("correlation: ", corr)
print("mutual information: ", mi)

noise = np.random.randn(n)
corr  = np.corrcoef(x,noise)
mi    = np.mean(mutual_info_regression(np.reshape(x, [-1,1]),noise))

plt.scatter(x,noise)
plt.xlabel('$\it{x}$')
plt.ylabel('$\it{y}$')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
plt.close()

print("correlation: ", corr)
print("mutual information: ", mi)

