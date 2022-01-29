# Contributing to pyEQUIB

The following guidelines are designed for contributors to the pyEQUIB Python package, which is 
hosted in the [pyEQUIB repository](https://github.com/equib/pyEQUIB) on GitHub. 

## Reporting Issues

The [issue tracker](https://github.com/equib/pyEQUIB/issues) is used to report bugs, request new functionality, and discuss improvements. 
For reporting a bug or a failed function or requesting a new feature, you can simply open an issue 
in the [issue tracker](https://github.com/equib/pyEQUIB/issues) of the 
[pyEQUIB repository](https://github.com/equib/pyEQUIB). If you are reporting a bug, please also include a minimal code
example that reproduces the problem, and Python version you are using.

## Contributing Code

Fo contributing code to pyEQUIB, you need to set up your [GitHub](https://github.com) 
account if you do not have and sign in, and request your change(s) or contribution via 
opening a pull request against the ``master``
branch in your fork of the [pyEQUIB repository](https://github.com/equib/pyEQUIB). 

To contribute to this package, you need to follow these steps:

- Open a new issue for new feature or failed function in the [Issue tracker](https://github.com/equib/pyEQUIB/issues).
- Fork the [pyEQUIB repository](https://github.com/equib/pyEQUIB) to your GitHub account.
- Clone your fork of the [pyEQUIB repository](https://github.com/equib/pyEQUIB):

      $ git clone git@github.com:your-username/pyEQUIB.git
      
- Make your change(s) in the `master` branch of your cloned fork.
- Make sure that it passes all tests and there is no error.
- Push yout change(s) to your fork in your GitHub account.
- [Submit a pull request][pr], mentioning what issue has been addressed.

[pr]: https://github.com/equib/pyEQUIB/compare/

Then, you are waiting, until your contribution is checked and merged into the original repository. 
We will contact you if there is a problem in your code.

While you are opening a pull request for your contribution, be sure that you have included:

* **Code** which you are contributing to this package.

* **Documentation** of this code if it provides new functionality. This should be a
  description of new functionality added to the API documentation (in ``docs/``). 

- **Tests** of this code to make sure that the previously failed function or the new functionality now works properly.

- **Revision history** if you fixed a bug in the previously failed function or add a code for new functionality, you should
well document your change(s) or addition in the *Revision History* entry of the changed or added function in your code.
---
title: "pyEQUIB Python Package, an addendum to proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis"
tags:
  - python
  - astrophysics
  - gaseous nebulae
  - plasma diagnostics
  - abundance analysis
authors:
  - name: Ashkbiz Danehkar
    orcid: 0000-0003-4552-5997
    affiliation: "1, 2, 3"
affiliations:
 - name: Department of Physics and Astronomy, Macquarie University, Sydney, NSW 2109, Australia
   index: 1
 - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA 
   index: 2
 - name: Department of Astronomy, University of Michigan, 1085 S. University Avenue, Ann Arbor, MI 48109, USA 
   index: 3
date: 20 October 2020
bibliography: paper.bib
---

# Addendum

`pyEQUIB` is a pure Python open-source package containing several application programming interface (API) functions that can be employed for plasma diagnostics and abundance analysis of nebular emission lines. This package is a Python implementation of the IDL library `proEQUIB` [@Danehkar:2018b] that is coupled to the IDL library `AtomNeb` [@Danehkar:2019]. The collisional excitation and recombination units of this package need to have the energy levels, collision strengths, transition probabilities, and recombination coefficients, which can be retrieved from the Python package `AtomNeb` for _Atomic Data of Ionized Nebulae_ [@Danehkar:2020]. The API functions of this package can be used to deduce the electron temperature, electron concentration, chemical elements from CELs and Rls, and the interstellar extinction from the Balmer decrements emitted from ionized gaseous nebulae. This package can simply be used by astronomers, who are familiar with the high-level, general-purpose programming language Python.

This package requires the Python packages `NumPy` [@Walt:2011; @Harris:2020], `SciPy` [@Virtanen:2020], and `AtomNeb` [@Danehkar:2020]. This package is released under the GNU General Public License. The source code is publicly available on the GitHub platform. The latest version of this package can be installed directly from its repository on the GitHub, and its stable version from the Python Package Index (PyPi) via ``pip install pyequib`` or alternatively from the Anaconda Python package distributor via ``conda install -c conda-forge pyequib``. The online documentation, tutorials and examples are available on its GitHub page (https://equib.github.io/pyEQUIB/) and its Read the Docs documentation page (https://pyequib.readthedocs.io/).

# Acknowledgements

The author acknowledges the support of Research Excellence Scholarship from Macquarie University.

# References
## Atomic Data

The atomic data are from the [AtomNeb database](https://github.com/atomneb/AtomNeb-py).

## pyEQUIB
[![PyPI version](https://badge.fury.io/py/pyequib.svg)](https://badge.fury.io/py/pyequib)
[![Build Status](https://travis-ci.org/equib/pyEQUIB.svg?branch=master)](https://travis-ci.org/equib/pyEQUIB)
[![Coverage Status](https://coveralls.io/repos/github/equib/pyEQUIB/badge.svg?)](https://coveralls.io/github/equib/pyEQUIB?branch=master)
[![GitHub license](https://img.shields.io/aur/license/yaourt.svg)](https://github.com/equib/pyEQUIB/blob/master/LICENSE)

Python package for atomic level populations and line emissivities in statistical equilibrium

### Description
**pyEQUIB** is a collection of [Python](https://www.python.org/) programs developed to perform plasma diagnostics and abundance analysis using emission line fluxes measured in ionzed nebulae. It uses the [AtomNeb Python Package](https://github.com/atomneb/AtomNeb-py) to read collision strengths and transition probabilities for collisionally excited lines (CEL), and recombination coefficients for recombination lines (RL). This Python package can be used to determine interstellar extinctions, electron temperatures, electron densities, and ionic abundances from the measured fluxes of emission lines.

### Installation
To install the last version, all you should need to do is

    python setup.py install

To install the stable version, you can use the preferred installer program (pip):

    pip install pyequib

### References

* Danehkar, A. (2018). proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis. *J. Open Source Softw.*, **3**, 899. doi:[10.21105/joss.00899](https://doi.org/10.21105/joss.00899) ads:[2018JOSS....3..899D](https://ui.adsabs.harvard.edu/abs/2018JOSS....3..899D).

This is a recipe for building the current development package into a conda binary.
## pyEQUIB
**pyEQUIB** - Python package for Plasma Diagnostics and Abundance Analysis

**pyEQUIB** is a collection of [Python](https://www.python.org/) programs developed to perform plasma diagnostics and abundance analysis using emission line fluxes measured in ionzed nebulae. It uses the [AtomNeb Python Package](https://github.com/atomneb/AtomNeb-py) to read collision strengths and transition probabilities for collisionally excited lines (CEL), and recombination coefficients for recombination lines (RL). This Python package can be used to determine interstellar extinctions, electron temperatures, electron densities, and ionic abundances from the measured fluxes of emission lines.

## Documentation

The [API Documentation](https://equib.github.io/pyEQUIB/doc/) is available on [equib.github.io/pyEQUIB](https://equib.github.io/pyEQUIB/). The documentation for AtomNeb is produced by [Sphinx](https://www.sphinx-doc.org) Python Documentation Generator using the [Sphinx RTD theme](https://pypi.org/project/sphinx-rtd-theme/). 





This sub-directory contains Python external packages which are required for use by the pyEQUIB Python Package.

The respective copyrights, restrictions and disclaimers of the original libraries apply for these procedures.

* ./atomneb/: required atomic data for ionized nebulae from the [AtomNeb Python Package](https://github.com/atomneb/AtomNeb-py)

