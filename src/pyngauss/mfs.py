import numpy as np
import healpy as hp

def compute_mfs(hp_map, mask=None, mask_threshold=0.9, n_nu=51, nu_min=-5.0, nu_max=5.0, dnu=0.2, iter=0):
    """
    Compute scalar Minkowski Functionals V0, V1, V2 for a HEALPix map.

    Parameters
    ----------
    hp_map : ndarray
        Input HEALPix map.
    mask : ndarray
        Analysis mask.
    mask_threshold : float, optional
        Pixels with mask > mask_threshold are used.
    n_nu : int, optional
        Number of threshold values.
    nu_min, nu_max : float, optional
        Range of normalized thresholds.
    dnu : float, optional
        Width used to approximate delta-function boundary terms.
    iter : int, optional
        Iteration count for healpy.map2alm.

    Returns
    -------
    nus : ndarray
        Threshold values.
    mfs : ndarray
        Array of shape (3, n_nu) containing V0, V1, V2.
    """

    if mask is None:
        mask = np.ones_like(hp_map)

    nside = hp.get_nside(hp_map)
    lmax = 3 * nside - 1
    npix = hp.nside2npix(nside)
    theta = hp.pix2ang(nside, np.arange(npix))[0]

    valid = mask > mask_threshold
    mean = np.mean(hp_map[valid])
    sigma0 = np.std(hp_map[valid])

    norm_map = (hp_map - mean) / sigma0
    nus = np.linspace(nu_min, nu_max, n_nu)

    alm = hp.map2alm(norm_map, lmax=lmax, iter=iter)
    map_der = hp.alm2map_der1(alm, nside, lmax=lmax)

    alm_ut = hp.map2alm(map_der[1], lmax=lmax, iter=iter)
    alm_up = hp.map2alm(map_der[2], lmax=lmax, iter=iter)

    map_der1 = hp.alm2map_der1(alm_ut, nside, lmax=lmax)
    map_der2 = hp.alm2map_der1(alm_up, nside, lmax=lmax)

    ut = map_der[1]
    up = map_der[2]
    cot = 1.0 / np.tan(theta)

    utt = map_der1[1]
    utp = map_der1[2] - cot * up
    upp = map_der2[2] + cot * ut

    grad = ut**2 + up**2
    B1 = np.sqrt(grad)
    B2 = (2 * ut * up * utp - ut**2 * upp - up**2 * utt) / grad

    mfs = np.zeros((3, n_nu))
    mask_sum = np.sum(valid)

    for i, nu in enumerate(nus):
        pix0 = (norm_map > nu) & valid
        pix1 = (norm_map >= (nu - 0.5 * dnu)) & (norm_map <= (nu + 0.5 * dnu)) & valid

        mfs[0, i] = np.sum(pix0) / mask_sum
        mfs[1, i] = 0.25 * np.sum(B1[pix1]) / (dnu * mask_sum)
        mfs[2, i] = (1.0 / 6.28) * np.sum(B2[pix1]) / (dnu * mask_sum)

    return nus, mfs