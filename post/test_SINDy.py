import numpy as np
import pysindy as ps
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def lorenz(t, state, sigma=10, rho=28, beta=8/3):
  x, y, z = state
  dxdt = sigma * (y - x)
  dydt = x * (rho - z) -y
  dzdt = x * y - beta * z
  return [dxdt, dydt, dzdt]


t_span = (0, 10)
dt = 0.01e0
t_eval = np.arange(t_span[0], t_span[1], dt)

x0 = [1.e0, 1.e0, 1.e0]

sol = solve_ivp(lorenz, t_span, x0, t_eval=t_eval, method='RK45')
X = sol.y.T

sindy = ps.SINDy(
  optimizer=ps.STLSQ(threshold=0.1),
  feature_library=ps.PolynomialLibrary(degree=3),
  differentiation_method=ps.FiniteDifference()
)

sindy.fit(X, t=dt)
sindy.print()

X_sim = sindy.simulate(x0, t_eval)

fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

for i, label in enumerate(["x", "y", "z"]):
  ax[i].plot(t_eval, X[:,i], 'k', label=f"True {label}(t)")
  ax[i].plot(t_eval, X_sim[:,i], '--r', label=f"SINDy {label}(t)")
  ax[i].legend()
  ax[i].set_ylabel(label)

plt.xlabel("Time")
plt.suptitle("Lorrenz System: True vs SINDy")
plt.show()


