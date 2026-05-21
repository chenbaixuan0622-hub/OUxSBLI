import numpy as np
from scipy.optimize import bisect


def free_stream(M0, gamma, R, p_tot, T_tot):
  p0 = p_tot / ((1.e0 + 0.5e0 * (gamma - 1.e0) * M0**2)**(gamma/(gamma-1.e0)))
  T0 = T_tot /  (1.e0 + 0.5e0 * (gamma - 1.e0) * M0**2)
  u0 = M0 * np.sqrt(gamma * R * T0)
  return u0, p0, T0


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


def oblique_shock(M0, p0, T0, beta_rad, gamma, R):
  Ms   = M0 * np.sin(beta_rad)
  Ms2  = Ms**2
  T2   = T0 * (1.e0 + 2.e0 * (gamma - 1.e0) * (Ms2 - 1.e0) * (1.e0 + gamma * Ms2) / (Ms2 * (gamma + 1.e0)**2))
  p2   = p0 * (1.e0 + 2.e0 * gamma * (Ms2 - 1.e0) / (gamma + 1.e0))
  rho2 = p2 / (R * T2)
  return rho2, p2, T2


def reslected_shock(M0, gamma, R, p_tot, T_tot, beta, beta_r):
  u0, p0, T0   = free_stream(M0, gamma, R, p_tot, T_tot)
  rho2, p2, T2 = oblique_shock(M0, p0, T0, beta, gamma, R)
  # incident shock
  Ms    = M0 * np.sin(beta)
  Ms2   = Ms**2
  theta = np.atan(2.e0 * (1.e0 / np.tan(beta)) * (Ms2 - 1.e0) / (M0**2 * (gamma + np.cos(2.e0 * beta)) + 2.e0))
  u1    = u0 * np.sin(beta)
  v1    = u0 * np.cos(beta)
  a1    = u0 / M0
  u2    = u1 - 2.e0 * a1 * (Ms - 1.e0 / Ms) / (gamma + 1.e0)
  v2    = u0 * np.cos(beta)
  ux    =   np.sqrt(u2**2 + v2**2) * np.cos(theta)
  uy    = - np.sqrt(u2**2 + v2**2) * np.sin(theta)
  # reflected shock
  a2    = np.sqrt(gamma * R * T2)
  M2    = np.sqrt(u2**2 + v2**2) / a2
  Mr    = M2 * np.sin(beta_r)
  Mr2   = Mr**2
  rho3, p3, T3 = oblique_shock(M2, p2, T2, beta_r, gamma, R)
  un2   = ux * np.sin(beta_r) - uy * np.cos(beta_r)
  ut2   = ux * np.cos(beta_r) + uy * np.sin(beta_r)
  un3   = un2 - 2.e0 * a2 * (Mr - 1.e0 / Mr) / (gamma + 1.e0)
  ut3   = ut2
  ux3   = un3 * np.sin(beta_r) + ut3 * np.cos(beta_r)
  uy3   =-un3 * np.cos(beta_r) + ut3 * np.sin(beta_r)
  return rho3, ux3, uy3, p3

