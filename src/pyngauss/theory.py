import numpy as np
from pyngauss.sk import compute_sk_parameters
from scipy.special import erf

def get_theory_mfs(map, mask=None, n_nu=51, nu_min=-5.0, nu_max=5.0, dnu=0.2):
    r"""
    Computes the theoretical Minkowski Functionals for a 2D Gaussian random field.
    
    Parameters:
    -----------
    nu : array-like
        Normalized thresholds (field value / s0).
    s0 : float
        Variance of the field (\sigma_0).
    s1 : float
        Variance of the first derivative of the field (\sigma_1).
        
    Returns:
    --------
    tuple : (nus, V0, V1, V2) theoretical curves.
    """

    if mask is None:
        mask = np.ones_like(map)

    nus = np.linspace(nu_min, nu_max, n_nu)

    s0, s1 = compute_sk_parameters(map, mask)[:2]

    # V0: Area Fraction
    v0 = 0.5 * (1.0 - erf(nus / np.sqrt(2.0)))
    
    # V1: Contour Length
    v1 = (1.0 / (8.0 * np.sqrt(2.0))) * (s1 / s0) * np.exp(-nus**2 / 2.0)
    
    # V2: Genus (Euler Characteristic)
    v2 = (1.0 / (2.0 * (2*np.pi)**1.5)) * (s1 / s0)**2 * nus * np.exp(-nus**2 / 2.0)
    
    return nus, v0, v1, v2