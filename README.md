# PyNGauss
PyNGauss is a Python package for computing non-Gaussian summary statistics of 2D random fields and sky maps. Its long-term goal is to provide a unified framework for morphology-based and harmonic-space estimators, including Minkowski Functionals, skewness–kurtosis statistics, bispectra, and related tools. The current release includes implementations of Minkowski Functionals and generalized skewness–kurtosis statistics.


## Installation

You can install `PyNGauss` directly from the source code. First, clone the repository to your local machine:

```bash
git clone [https://github.com/yourusername/pyngauss.git](https://github.com/yourusername/pyngauss.git)
cd pyngauss


# Standard installation
pip install .

# Editable/development installation (Recommended)
pip install -e .

References
If you use PyNGauss in your research, please consider citing the theoretical framework and methodological papers that form the basis of the algorithms used in this package:

Pixel Window Corrections (Lim-Simon): Lim, E. A., & Simon, D. (2012). Minkowski functionals of cosmic microwave background temperature and polarization fields. Journal of Cosmology and Astroparticle Physics, 2012(01), 048. arXiv:1103.4300

Continuous Sphere Theory: Schmalzing, J., & Górski, K. M. (1998). Minkowski functionals used in the morphological analysis of cosmic microwave background anisotropy maps. Monthly Notices of the Royal Astronomical Society, 297(2), 355-365. arXiv:astro-ph/9710185
