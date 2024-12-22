import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

x = np.linspace(0, 10, 20).reshape(-1,1)
y = np.sin(x).ravel() + 0.1 * np.random.randn(x.shape[0])

x_pred = np.linspace(0, 10, 100).reshape(-1, 1)

kernel = C(1.0, (1e-2, 1e2)) * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2))

gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=0.01)

gpr.fit(x, y)

y_pred, sigma = gpr.predict(x_pred, return_std=True)

plt.scatter(x, y, c='red')
plt.plot(x_pred, np.sin(x_pred))
plt.plot(x_pred, y_pred)
plt.fill_between(x_pred.ravel(), y_pred - 2*sigma, y_pred + 2*sigma, alpha=0.2)
plt.show()

