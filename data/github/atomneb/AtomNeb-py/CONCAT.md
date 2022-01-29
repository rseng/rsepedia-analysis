# Contributing to AtomNeb-py

The following guidelines are designed for contributors to the AtomNeb Python package, which is 
hosted in the [AtomNeb-py repository](https://github.com/atomneb/AtomNeb-py) on GitHub. 

## Reporting Issues

The [issue tracker](https://github.com/atomneb/AtomNeb-py/issues) is used to report any bugs, request new functionality, and discuss improvements. 
For reporting a bug or a failed function or requesting a new feature, you can simply open an issue 
in the [issue tracker](https://github.com/atomneb/AtomNeb-py/issues) of the 
[AtomNeb-py](https://github.com/atomneb/AtomNeb-py) repository. If you are reporting a bug, please also include a minimal code
example that makes the issue, and Python version used by you.

## Contributing Code

To make contributions to AtomNeb-py, you need to set up your [GitHub](https://github.com) 
account if you do not have and sign in, and request your change(s) or contribution via 
opening a pull request against the ``master``
branch in your fork of the [AtomNeb-py repository](https://github.com/atomneb/AtomNeb-py). 

Please use the following steps:

- Open a new issue for new feature or failed function in the [Issue tracker](https://github.com/atomneb/AtomNeb-py/issues).
- Fork the [AtomNeb-py repository](https://github.com/atomneb/AtomNeb-py) to your GitHub account.
- Clone your fork of the [AtomNeb-py repository](https://github.com/atomneb/AtomNeb-py):

      $ git clone git@github.com:your-username/AtomNeb-py.git
      
- Make your change(s) in the `master` branch of your cloned fork.
- Make sure that all tests are passed without any errors.
- Push yout change(s) to your fork in your GitHub account.
- [Submit a pull request][pr], mentioning what problem has been solved.

[pr]: https://github.com/atomneb/AtomNeb-py/compare/

Your contribution will be checked and merged into the original repository. You will be contacted if there is any problem in your contribution.

While you are opening a pull request for your contribution, be sure that you have included:

* **Code** which you are contributing to this package.

* **Documentation** of this code if it provides new functionality. This should be a
  description of new functionality added to the API documentation (in ``docs/``). 

- **Tests** of this code to make sure that the previously failed function or the new functionality now works properly.

- **Revision history** if you fixed a bug in the previously failed function or add a code for new functionality, you should
well document your change(s) or addition in the *Revision History* entry of the changed or added function in your code.
## AtomNeb-py
**AtomNeb-py** - Python package for Atomic Data of Ionized Nebulae

**AtomNeb-py** is a Python package for reading atomic data from **AtomNeb**, which is a database containing atomic data stored in the Flexible Image Transport System (FITS) file format for *collisionally excited lines* and *recombination lines* typically observed in spectra of ionized gaseous nebulae. The AtomNeb database were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB), [proEQUIB](https://github.com/equib/proEQUIB), and other nebular spectral analysis tools. 


---
title: "AtomNeb Python Package, an addendum to AtomNeb: IDL Library for Atomic Data of Ionized Nebulae"
tags:
  - python
  - astrophysics
  - atomic data
  - gaseous nebulae
  - spectral analysis
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

`AtomNeb` is a Python open-source package containing atomic data for gaseous nebulae stored in the Flexible Image Transport System (FITS) file format [@Wells:1981; @Hanisch:2001; @Pence:2010]. These FITS files offer easy access to the atomic data required for emissivity calculations in the collisional excitation and recombination processes usually occurred in ionized gases of planetary nebulae and H II regions. This package has several application programming interface (API) functions developed in Python for retrieving the energy levels, collision strengths, transition probabilities, and recombination coefficients from its FITS files. The previous library `AtomNeb` [@Danehkar:2019] coupled to the library `proEQUIB` [@Danehkar:2018b] needs the Interactive Data Language (IDL) compiler, so this package offers an identical package for the high-level programming language Python that can be used by those astrophysicists, who intend to analyze nebular emission lines by developing codes in Python. The `AtomNeb` Python functions can be used by the Python package `pyEQUIB` [@Danehkar:2020] to analyze emission-line spectra. 

`AtomNeb` uses the FITS handling routines of the Python package `Astropy` [@Astropy:2013; @Astropy:2018] to retrieve the atomic data from its FITS files. It also requires the Python packages `NumPy` [@Walt:2011; @Harris:2020] and `pandas` [@McKinney:2010; @McKinney:2011; @McKinney:2017]. This package is released under the GNU General Public License, and its source code is publicly available on its GitHub repository. Its latest version can be installed directly from its repository on the GitHub, and the stable version from the Python Package Index (PyPi) via ``pip install atomneb`` or alternatively from the Conda Python package manager via ``conda install -c conda-forge atomneb``. The online documentation, tutorials and examples are provided on the GitHub platform (https://github.com/atomneb/AtomNeb-py) and the Read the Docs documentation host (https://atomneb-py.readthedocs.io/).


# Acknowledgements

AD acknowledges the support of Research Excellence Scholarship from Macquarie University.

# References
## AtomNeb

**Collisionally Excitation Atomic Data for Ionized Nebulae**

The collisionally excitation atomic data were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB) and [proEQUIB](https://github.com/equib/proEQUIB):

* **Collection** from the [National Institute of Standards and Technology (NIST) Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database), the [CHIANTI atomic database](http://www.chiantidatabase.org/), and some improved atomic data from Cloudy v13.04 and pyNeb v1.0.

* **Chianti52** from the [CHIANTI atomic database](http://www.chiantidatabase.org/) version 5.2.

* **Chianti60** from the [CHIANTI atomic database](http://www.chiantidatabase.org/) version 6.0.

* **Chianti70** from the [CHIANTI atomic database](http://www.chiantidatabase.org/) version 7.0.
## AtomNeb

**Atomic Data for Ionized Nebulae**

The atomic data were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB) and [proEQUIB](https://github.com/equib/proEQUIB) obtained from the [CHIANTI atomic database](http://www.chiantidatabase.org/) version 7.0:

* **AtomElj.fits** contains *Energy Levels* (Ej).

* **AtomOmij.fits** contains *Collision Strengths* (Ωij).

* **AtomAij.fits** contains *Transition Probabilities* (Aij). 
## AtomNeb

**Atomic Data for Ionized Nebulae**

The atomic data were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB) and [proEQUIB](https://github.com/equib/proEQUIB) obtained from the [CHIANTI atomic database](http://www.chiantidatabase.org/) version 5.2:

* **AtomElj.fits** contains *Energy Levels* (Ej).

* **AtomOmij.fits** contains *Collision Strengths* (Ωij).

* **AtomAij.fits** contains *Transition Probabilities* (Aij). 
## AtomNeb

**Atomic Data for Ionized Nebulae**

The atomic data were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB) and [proEQUIB](https://github.com/equib/proEQUIB) obtained from the [CHIANTI atomic database](http://www.chiantidatabase.org/) version 9.0:

* **AtomElj.fits** contains *Energy Levels* (Ej).

* **AtomOmij.fits** contains *Collision Strengths* (Ωij).

* **AtomAij.fits** contains *Transition Probabilities* (Aij). 
## AtomNeb

**Atomic Data for Ionized Nebulae**

The collection of atomic data was generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB) and [proEQUIB](https://github.com/equib/proEQUIB):

* **AtomElj.fits** contains *Energy Levels* (Ej) obtained from the [National Institute of Standards and Technology (NIST) Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database), and some improvements from Cloudy v13.04 and pyNeb v1.0.

* **AtomOmij.fits** contains *Collision Strengths* (Ωij) obtained from the [National Institute of Standards and Technology (NIST) Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database), the [CHIANTI atomic database](http://www.chiantidatabase.org/), and some improved atomic data from Cloudy v13.04 and pyNeb v1.0. 

* **AtomAij.fits** contains *Transition Probabilities* (Aij) obtained from the [National Institute of Standards and Technology (NIST) Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database), the [CHIANTI atomic database](http://www.chiantidatabase.org/), and some improved atomic data from Cloudy v13.04 and pyNeb v1.0. 
## AtomNeb

**Atomic Data for Ionized Nebulae**

The atomic data were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB) and [proEQUIB](https://github.com/equib/proEQUIB) obtained from the [CHIANTI atomic database](http://www.chiantidatabase.org/) version 6.0:

* **AtomElj.fits** contains *Energy Levels* (Ej).

* **AtomOmij.fits** contains *Collision Strengths* (Ωij).

* **AtomAij.fits** contains *Transition Probabilities* (Aij). 
## AtomNeb

**Recombination Atomic Data for Ionized Nebulae**

The recombination coefficient atomic data were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB) and [proEQUIB](https://github.com/equib/proEQUIB):

* **rc_collection.fits**, effective recombination coefficients for C II ([Davey et al. 2000](http://adsabs.harvard.edu/abs/2000A%26AS..142...85D)), N II ([Escalante and Victor 1990](http://adsabs.harvard.edu/abs/1990ApJS...73..513E)), O II ([Storey 1994](http://adsabs.harvard.edu/abs/1994A%26A...282..999S); [Liu et al. 1995](http://adsabs.harvard.edu/abs/1995MNRAS.272..369L)), and Ne II ions ([Kisielius et al. 1998](http://adsabs.harvard.edu/abs/1998A%26AS..133..257K)), including Branching ratios (Br) for O II and N II ions.

* **rc_SH95.fits**, hydrogenic ions for Z=1 to 8, namely H I, He II, Li III, Be IV, B V, C VI, N VII, and O VIII ions from [Storey and Hummer (1995)](http://adsabs.harvard.edu/abs/1995MNRAS.272...41S).

* **rc_PPB91.fits**, effective recombination coefficients for H, He, C, N, O, Ne ions from [Pequignot, Petitjean and Boisson (1991)](http://adsabs.harvard.edu/abs/1991A%26A...251..680P).

* **rc_he_ii_PFSD12.fits**, effective He I recombination coefficients from [Porter et al (2012)](
http://adsabs.harvard.edu/abs/2012MNRAS.425L..28P) and [(2013)](http://adsabs.harvard.edu/abs/2013MNRAS.433L..89P).

* **rc_n_iii_FSL13.fits**, effective N II recombination coefficients from [Fang, Storey and Liu (2011)](
http://adsabs.harvard.edu/abs/2011A%26A...530A..18F) and [(2013)](http://adsabs.harvard.edu/abs/2013A%26A...550C...2F).

* **rc_o_iii_SSB17.fits**, effective O II recombination coefficients of 8889 recombination lines for Cases A, B, and C from [Storey, Sochi and Bastin (2017)](
http://adsabs.harvard.edu/abs/2017MNRAS.470..379S). (Use these commends to unpack: tar -xvf rc_o_iii_SSB17.fits.tar.gz)

* **rc_o_iii_SSB17_orl_case_b.fits**, effective O II recombination coefficients of 2433 optical recombination lines for Case B (Wavelength: 3500-9000Å) from [Storey, Sochi and Bastin (2017)](
http://adsabs.harvard.edu/abs/2017MNRAS.470..379S). (Use this commend to unpack: tar -xvf rc_o_iii_SSB17_orl_case_b.fits.tar.gz)
## atomneb
[![PyPI version](https://badge.fury.io/py/atomneb.svg)](https://badge.fury.io/py/pyemcee)
[![Build Status](https://travis-ci.org/atomneb/AtomNeb-py.svg?branch=master)](https://travis-ci.org/atomneb/AtomNeb-py)
[![Coverage Status](https://coveralls.io/repos/github/atomneb/AtomNeb-py/badge.svg?branch=master)](https://coveralls.io/github/atomneb/AtomNeb-py?branch=master)
[![GitHub license](https://img.shields.io/badge/license-GPL-blue.svg)](https://github.com/mcfit/pyemcee/blob/master/LICENSE)

Python Package for Atomic Data of Ionized Nebulae

### Description

**AtomNeb-py** is a Python package for reading atomic data from *AtomNeb*, which is a database containing atomic data stored in the Flexible Image Transport System (FITS) file format for *collisionally excited lines* and *recombination lines* typically observed in spectra of ionized gaseous nebulae. The AtomNeb database were generated for use in [pyEQUIB](https://github.com/equib/pyEQUIB), [proEQUIB](https://github.com/equib/proEQUIB), and other nebular spectral analysis tools. 

### Installation
To install the last version, all you should need to do is

    python setup.py install

To install the stable version, you can use the preferred installer program (pip):

    pip install atomneb


This is a recipe for building the current development package into a conda binary.
## Documentation

The [API Documentation](https://atomneb.github.io/AtomNeb-py/doc/) is available on [atomneb.github.io/AtomNeb-py](https://atomneb.github.io/AtomNeb-py/). The documentation for AtomNeb is produced by [Sphinx](https://www.sphinx-doc.org) Python Documentation Generator using the [Sphinx RTD theme](https://pypi.org/project/sphinx-rtd-theme/). 





