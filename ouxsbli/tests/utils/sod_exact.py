"""Exact solution for the Sod shock-tube problem (Toro, Chapter 4)."""

import numpy as np


def solve(x, t, rho_l=1.0, u_l=0.0, p_l=1.0,
          rho_r=0.125, u_r=0.0, p_r=0.1,
          gamma=1.4, x0=0.5):
    """Return (rho, u, p) arrays at positions x and time t.

    The diaphragm is initially at x0.  Solves for p_star with Newton–Raphson,
    then reconstructs all five regions.
    """
    g  = gamma
    g1 = (g - 1.0) / (2.0 * g)
    g2 = (g + 1.0) / (2.0 * g)
    g3 = 2.0 * g / (g - 1.0)
    g4 = 2.0 / (g - 1.0)
    g5 = 2.0 / (g + 1.0)
    g6 = (g - 1.0) / (g + 1.0)
    g7 = (g - 1.0) / 2.0

    a_l = np.sqrt(g * p_l / rho_l)
    a_r = np.sqrt(g * p_r / rho_r)

    # Shock / rarefaction functions
    def f_k(p, p_k, rho_k, a_k):
        if p > p_k:
            A = g5 / rho_k
            B = g6 * p_k
            return (p - p_k) * np.sqrt(A / (p + B))
        else:
            return g4 * a_k * ((p / p_k) ** g1 - 1.0)

    def df_k(p, p_k, rho_k, a_k):
        if p > p_k:
            A = g5 / rho_k
            B = g6 * p_k
            return np.sqrt(A / (p + B)) * (1.0 - (p - p_k) / (2.0 * (p + B)))
        else:
            return 1.0 / (rho_k * a_k) * (p / p_k) ** (-(g + 1.0) / (2.0 * g))

    # Newton–Raphson for p_star
    p_star = max(1e-6, 0.5 * (p_l + p_r))
    for _ in range(100):
        f  = f_k(p_star, p_l, rho_l, a_l) + f_k(p_star, p_r, rho_r, a_r) + (u_r - u_l)
        df = df_k(p_star, p_l, rho_l, a_l) + df_k(p_star, p_r, rho_r, a_r)
        dp = -f / df
        p_star += dp
        if abs(dp) / (p_star + 1e-12) < 1e-10:
            break

    u_star = 0.5 * (u_l + u_r) + 0.5 * (f_k(p_star, p_r, rho_r, a_r)
                                          - f_k(p_star, p_l, rho_l, a_l))

    # Densities in star regions
    if p_star > p_l:
        rho_sl = rho_l * ((p_star / p_l + g6) / (g6 * p_star / p_l + 1.0))
    else:
        rho_sl = rho_l * (p_star / p_l) ** (1.0 / g)

    if p_star > p_r:
        rho_sr = rho_r * ((p_star / p_r + g6) / (g6 * p_star / p_r + 1.0))
    else:
        rho_sr = rho_r * (p_star / p_r) ** (1.0 / g)

    a_sl = np.sqrt(g * p_star / rho_sl)
    a_sr = np.sqrt(g * p_star / rho_sr)

    # Wave speeds
    if p_star > p_l:
        s_l = u_l - a_l * np.sqrt(g2 * p_star / p_l + g1)  # left shock
        s_hl = s_tl = s_l
    else:
        s_hl = u_l - a_l                  # left rarefaction head
        s_tl = u_star - a_sl              # left rarefaction tail

    if p_star > p_r:
        s_r = u_r + a_r * np.sqrt(g2 * p_star / p_r + g1)  # right shock
        s_hr = s_tr = s_r
    else:
        s_hr = u_r + a_r                  # right rarefaction head
        s_tr = u_star + a_sr              # right rarefaction tail

    s = (x - x0) / (t + 1e-300)           # similarity variable

    rho_out = np.empty_like(x, dtype=float)
    u_out   = np.empty_like(x, dtype=float)
    p_out   = np.empty_like(x, dtype=float)

    for idx in range(len(x)):
        si = s[idx]
        if si <= s_hl:                     # left undisturbed
            rho_out[idx], u_out[idx], p_out[idx] = rho_l, u_l, p_l
        elif si <= s_tl:                   # inside left rarefaction
            u_i   = g5 * (a_l + g7 * u_l + si)
            a_i   = g5 * (a_l + g7 * (u_l - si))
            rho_i = rho_l * (a_i / a_l) ** g4
            p_i   = p_l   * (a_i / a_l) ** g3
            rho_out[idx], u_out[idx], p_out[idx] = rho_i, u_i, p_i
        elif si <= u_star:                 # left star region
            rho_out[idx], u_out[idx], p_out[idx] = rho_sl, u_star, p_star
        elif si <= s_tr:                   # right star region
            rho_out[idx], u_out[idx], p_out[idx] = rho_sr, u_star, p_star
        elif si <= s_hr:                   # inside right rarefaction (if applicable)
            u_i   = g5 * (-a_r + g7 * u_r + si)
            a_i   = g5 * (a_r - g7 * (u_r - si))
            rho_i = rho_r * (a_i / a_r) ** g4
            p_i   = p_r   * (a_i / a_r) ** g3
            rho_out[idx], u_out[idx], p_out[idx] = rho_i, u_i, p_i
        else:                              # right undisturbed
            rho_out[idx], u_out[idx], p_out[idx] = rho_r, u_r, p_r
    return rho_out, u_out, p_out
