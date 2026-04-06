import numpy as np
import healpy as hp

def compute_sk_parameters(hp_map, mask, mask_threshold=0.9, iter=0):
    """
    Compute generalized skewness and kurtosis parameters for a HEALPix map.

    Returns
    -------
    tuple
        (sigma0, sigma1, S0, S1, S2, K0, K1, K2, K3)
    """
    nside = hp.get_nside(hp_map)
    lmax = 3 * nside - 1
    npix = len(hp_map)

    theta = hp.pix2ang(nside, np.arange(npix))[0]
    valid = mask > mask_threshold

    mean = np.mean(hp_map[valid])
    sigma0 = np.std(hp_map[valid])

    map1 = hp_map - mean

    alm = hp.map2alm(map1, lmax=lmax, iter=iter)
    map_1der = hp.alm2map_der1(alm, nside, lmax=lmax)

    alm_ut = hp.map2alm(map_1der[1], lmax=lmax, iter=iter)
    alm_up = hp.map2alm(map_1der[2], lmax=lmax, iter=iter)

    map_2der1 = hp.alm2map_der1(alm_ut, nside, lmax=lmax)
    map_2der2 = hp.alm2map_der1(alm_up, nside, lmax=lmax)

    ut = map_1der[1]
    up = map_1der[2]
    cot = 1.0 / np.tan(theta)

    utt = map_2der1[1]
    upp = map_2der2[2] + cot * ut

    grad = ut**2 + up**2
    laplacian = utt + upp + cot * ut

    sigma1 = np.sqrt(np.mean(grad[valid]))

    S0 = np.mean(map1[valid]**3) / sigma0**4
    S1 = np.mean((map1[valid]**2) * laplacian[valid]) / (sigma0**2 * sigma1**2)
    S2 = 2.0 * np.mean(grad[valid] * laplacian[valid]) / sigma1**4

    K0 = (np.mean(map1[valid]**4) - 3.0 * sigma0**4) / sigma0**6
    K1 = -3.0 * (np.mean(map1[valid]**2 * grad[valid]) - sigma0**2 * sigma1**2) / (sigma0**4 * sigma1**2)

    K21 = 2.0 * (np.mean(map1[valid] * grad[valid] * laplacian[valid]) + sigma1**4) / (sigma0**2 * sigma1**4)
    K22 = (np.mean(grad[valid]**2) - 2.0 * sigma1**4) / (sigma0**2 * sigma1**4)

    K2 = K21 + K22
    K3 = 0.5 * K22

    return sigma0, sigma1, S0, S1, S2, K0, K1, K2, K3