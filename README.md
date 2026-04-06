# PyNGauss

PyNGauss is a Python package for computing non-Gaussian summary statistics of 2D random fields and sky maps. Its long-term goal is to provide a unified framework for morphology-based and harmonic-space estimators, including Minkowski Functionals, skewness–kurtosis statistics, bispectra, and related tools. The current release includes implementations of Minkowski Functionals and generalized skewness–kurtosis statistics.

## Installation

You can install `PyNGauss` directly from the source code. First, clone the repository to your local machine:

```bash
git clone [https://github.com/fazlurahman/pyngauss.git](https://github.com/fazlurahman/pyngauss.git)
cd pyngauss
```

### Standard Installation

To install the package normally, run:

```bash
pip install .
```

### Editable/Development Installation (Recommended)

If you plan to modify the code or run notebooks directly from the repository, install it in editable mode so your changes take effect immediately:

```bash
pip install -e .
```

## Citation & References

If you use `PyNGauss` in your research, please cite the primary package paper:

* **Rahman, et al. (2026).** *Understanding the non-Gaussian nature of Galactic foreground emissions towards small scales* 

The theoretical framework and numerical algorithms implemented in this code are based on the following foundational papers:

* **Continuous Sphere Theory:** Schmalzing, J., & Górski, K. M. (1998). Minkowski Functionals used in the Morphological Analysis of Cosmic Microwave Background Anisotropy Maps *Monthly Notices of the Royal Astronomical Society*, 297(2), 355-365. [arXiv:astro-ph/9710185](https://arxiv.org/abs/astro-ph/9710185)
* **Pixel Window Corrections (Lim-Simon):** Lim, E. A., & Simon, D. (2012). Can we detect Hot or Cold spots in the CMB with Minkowski Functionals?. *Journal of Cosmology and Astroparticle Physics*, 2012(01), 048. [arXiv:1103.4300](https://arxiv.org/abs/1103.4300)
