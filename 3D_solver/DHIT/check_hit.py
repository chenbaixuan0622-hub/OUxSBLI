import numpy as np
import matplotlib.pyplot as plt
from myvtk import getGrid, getVector


filepath = "./data/Q00016.vtr"
Nx, Ny, Nz, x, y, z = getGrid(filepath)
u, v, w = getVector(filepath, Nx, Ny, Nz, 'velocity')


def compute_3d_energy_spectrum(
    u,
    v=None,
    w=None,
    lx=2.0*np.pi,
    ly=2.0*np.pi,
    lz=2.0*np.pi,
    detrend=True,
):
    """
    Compute isotropic 3D energy spectrum E(k).

    Parameters
    ----------
    u, v, w : ndarray
        Velocity components with shape (nx, ny, nz).

        If only u is given:
            spectrum of u is computed.

        If u, v, w are all given:
            kinetic energy spectrum is computed.

    lx, ly, lz : float
        Domain lengths.

    detrend : bool
        Remove spatial mean before FFT.

    Returns
    -------
    k_bin : ndarray
        Wavenumber bins.

    E_bin : ndarray
        Energy spectrum.
    """

    nx, ny, nz = u.shape

    # ----------------------------------------
    # Remove mean
    # ----------------------------------------
    if detrend:
        u = u - np.mean(u)

    if v is not None:
        if detrend:
            v = v - np.mean(v)

    if w is not None:
        if detrend:
            w = w - np.mean(w)

    # ----------------------------------------
    # FFT
    # ----------------------------------------
    uhat = np.fft.fftn(u)

    if (v is not None) and (w is not None):

        vhat = np.fft.fftn(v)
        what = np.fft.fftn(w)

        energy_density = (
            np.abs(uhat)**2
            + np.abs(vhat)**2
            + np.abs(what)**2
        ) / 2.0

    else:

        energy_density = 0.5 * np.abs(uhat)**2

    # normalization
    energy_density /= (nx * ny * nz)**2

    # ----------------------------------------
    # Wavenumbers
    # ----------------------------------------
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=lx/nx)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=ly/ny)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=lz/nz)

    KX, KY, KZ = np.meshgrid(
        kx,
        ky,
        kz,
        indexing="ij",
    )

    kmag = np.sqrt(
        KX**2
        + KY**2
        + KZ**2
    )

    # ----------------------------------------
    # Radial binning
    # ----------------------------------------
    dk = min(
        2.0*np.pi/lx,
        2.0*np.pi/ly,
        2.0*np.pi/lz,
    )

    kmax = np.max(kmag)

    nbins = int(kmax / dk)

    k_bin = np.zeros(nbins)
    E_bin = np.zeros(nbins)

    kmag_flat = kmag.ravel()
    energy_flat = energy_density.ravel()

    inds = np.floor(kmag_flat / dk).astype(int)

    for i in range(nbins):

        mask = inds == i

        if np.any(mask):

            E_bin[i] = np.sum(
                energy_flat[mask]
            )

            k_bin[i] = (i + 0.5) * dk

    # remove empty bins
    mask = E_bin > 0.0

    k_bin = k_bin[mask]
    E_bin = E_bin[mask]

    return k_bin, E_bin


def plot_energy_spectrum(
    k,
    E1,
    E2,
    E3,
    ref_slope=-5.0/3.0,
    plot_reference=True,
):
    """
    Plot energy spectrum.

    Parameters
    ----------
    k : ndarray
        Wavenumbers.

    E : ndarray
        Energy spectrum.

    ref_slope : float
        Reference slope.

    label : str
        Legend label.

    plot_reference : bool
        Plot reference slope line.
    """

    plt.figure(figsize=(7, 5))

    # ----------------------------------------
    # Main spectrum
    # ----------------------------------------
    ylabel = r"$E(k)$"
    plt.loglog(k, E3, linewidth=2, label="w")
    plt.loglog(k, E1, linewidth=2, label="u")
    plt.loglog(k, E2, linewidth=2, label="v")
    plt.legend()
    # ----------------------------------------
    # Reference slope
    # ----------------------------------------
    if plot_reference:

        i0 = len(k) // 3

        kref = k[i0:]

        if len(kref) > 0:

            cref = E1[i0] / (k[i0]**ref_slope)

            yref = cref * kref**ref_slope

            plt.loglog(
                kref,
                yref,
                "--",
                linewidth=1.5,
                label=rf"$k^{{{ref_slope:.2f}}}$",
            )

    plt.xlabel(r"$k$")
    plt.ylabel(ylabel)

    plt.grid(True, which="both", alpha=0.3)

    plt.legend()

    plt.tight_layout()

    plt.show()




k, Eu = compute_3d_energy_spectrum(u, lx=2.0*np.pi, ly=2.0*np.pi, lz=2.0*np.pi)
k, Ev = compute_3d_energy_spectrum(u, lx=2.0*np.pi, ly=2.0*np.pi, lz=2.0*np.pi)
k, Ew = compute_3d_energy_spectrum(u, lx=2.0*np.pi, ly=2.0*np.pi, lz=2.0*np.pi)

plot_energy_spectrum(k, Eu, Ev, Ew, ref_slope=-5.0/3.0, plot_reference=True)

