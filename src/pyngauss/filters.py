import numpy as np
import healpy as hp

def apply_bandpass_filter(hp_map, ell_min, ell_max, delta_ell=10.0, lmax=None, iter=0):
    """
    Apply a smooth tanh bandpass filter to a HEALPix map.

    Parameters
    ----------
    hp_map : ndarray
        Input HEALPix map.
    ell_min : int or float
        Lower multipole cutoff.
    ell_max : int or float
        Upper multipole cutoff.
    delta_ell : float, optional
        Width of the tanh transition.
    lmax : int or None, optional
        Maximum multipole. Defaults to 3*nside-1.
    iter : int, optional
        Iteration count for healpy.map2alm.

    Returns
    -------
    ndarray
        Filtered HEALPix map.
    """
    nside = hp.get_nside(hp_map)
    if lmax is None:
        lmax = 3 * nside - 1

    alm = hp.map2alm(hp_map, lmax=lmax, iter=iter)

    bp_filter = np.zeros(lmax + 1)
    ells = np.arange(lmax + 1)

    x_low = (ells - ell_min) / delta_ell
    x_high = (ells - ell_max) / delta_ell

    if ell_min == 0:
        bp_filter = 0.5 * (1.0 - np.tanh(x_high))
    else:
        bp_filter = 0.5 * (1.0 + np.tanh(x_low)) * 0.5 * (1.0 - np.tanh(x_high))

    hp.almxfl(alm, bp_filter, inplace=True)
    return hp.alm2map(alm, nside=nside, lmax=lmax, verbose=False)