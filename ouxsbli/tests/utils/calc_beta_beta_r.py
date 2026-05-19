import numpy as np
from scipy.optimize import bisect


def beta(M, theta_deg, gamma=1.4):
    theta = np.radians(theta_deg)

    f = lambda b: (
        np.tan(theta)
        - 2/np.tan(b)
        * (M**2*np.sin(b)**2 - 1)
        / (M**2*(gamma + np.cos(2*b)) + 2)
    )
    mu = np.arcsin(1/M)
    # weak solution only
    beta_max = np.radians(45)
    return np.degrees(bisect(f, mu + 1e-6, beta_max))


def downstream_mach(M1, beta_deg, theta_deg, gamma=1.4):
    beta = np.radians(beta_deg)
    theta = np.radians(theta_deg)
    Mn1 = M1 * np.sin(beta)
    Mn2 = np.sqrt(
        (1 + 0.5*(gamma-1)*Mn1**2)
        / (gamma*Mn1**2 - 0.5*(gamma-1))
    )
    return Mn2 / np.sin(beta - theta)


M1 = 2.0
theta = 8.0

beta_i = beta(M1, theta)
M2 = downstream_mach(M1, beta_i, theta)
beta_r = beta(M2, theta)

print(f"beta_i = {beta_i:.3f} deg")
print(f"beta_r = {beta_r:.3f} deg")
