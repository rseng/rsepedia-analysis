# disksurf

<p align='center'>
  <img src="HD163296_zeroth.png" width="793" height="549">
  <br>
  <a href='https://disksurf.readthedocs.io/en/latest/?badge=latest'>
    <img src='https://readthedocs.org/projects/disksurf/badge/?version=latest' alt='Documentation Status' />
  </a>
  <a href='https://doi.org/10.21105/joss.03827'>
    <img src='https://joss.theoj.org/papers/10.21105/joss.03827/status.svg' alt='DOI'>
  </a>
  <a href="https://zenodo.org/badge/latestdoi/184391824">
    <img src="https://zenodo.org/badge/184391824.svg" alt="DOI">
  </a>
</p>

## What is it?

`disksurf` is a package which contains the functions to measure the height of optically thick emission, or photosphere, using the method presented in [Pinte et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...609A..47P/abstract) (with associated [example script](https://github.com/cpinte/CO_layers)).

## How do I install it?

Grab the latest version from PyPI:

```
$ pip install disksurf
```

This has a couple of dependencies, namely [astropy](https://github.com/astropy/astropy) and [GoFish](https://github.com/richteague/gofish), which should be installed automatically if you don't have them. To verify that everything was installed as it should, running through the [tutorials](https://disksurf.readthedocs.io/en/latest/tutorials/tutorial_1.html) should work without issue.

## How do I use it?

At its most basic, `disksurf` is as easy as:

```python
from disksurf import observation                        # import the module
cube = observations('path/to/cube.fits')                # load up the data
surface = cube.get_emission_surface(inc=30.0, PA=45.0)  # extract the surface
surface.plot_surface()                                  # plot the surface
```

Follow our [tutorials](https://disksurf.readthedocs.io/en/latest/tutorials/tutorial_1.html) for a quick guide on how to use `disksurf` with DSHARP data and some of the additional functions that will help you extract the best surface possible.

## Citation

If you use this software, please remember to cite both [Pinte et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...609A..47P/abstract) for the method, and [Teague et al. (2021)](https://joss.theoj.org/papers/10.21105/joss.03827#) for the software.

```
@article{2018A&A...609A..47P,
  doi = {10.1051/0004-6361/201731377},
  year = {2018},
  volume = {609},
  eid = {A47},
  pages = {A47},
  author = {{Pinte}, C. and {M{\'e}nard}, F. and {Duch{\^e}ne}, G. and {Hill}, T. and {Dent}, W.~R.~F. and {Woitke}, P. and {Maret}, S. and {van der Plas}, G. and {Hales}, A. and {Kamp}, I. and {Thi}, W.~F. and {de Gregorio-Monsalvo}, I. and {Rab}, C. and {Quanz}, S.~P. and {Avenhaus}, H. and {Carmona}, A. and {Casassus}, S.},
  title = "{Direct mapping of the temperature and velocity gradients in discs. Imaging the vertical CO snow line around IM Lupi}",
  journal = {\aap}
}

@article{disksurf,
  doi = {10.21105/joss.03827},
  url = {https://doi.org/10.21105/joss.03827},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {67},
  pages = {3827},
  author = {Richard Teague and Charles J. Law and Jane Huang and Feilong Meng},
  title = {disksurf: Extracting the 3D Structure of Protoplanetary Disks},
  journal = {Journal of Open Source Software}
}
```
# Contributing to `disksurf`

If you are interested in contributing, please submit a pull request. New functions should be documented in a manner consistent with existing code. Software bugs or mistakes in the documentation should be reported by [opening up an issue](https://github.com/richteague/disksurf/issues). If reporting a software bug, please provide a minimal reproducible example.
---
title: 'disksurf: Extracting the 3D Structure of Protoplanetary Disks'
tags:
  - Python
  - astronomy
  - protoplanetary disks
authors:
  - name: Richard Teague
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Charles J. Law
    orcid: 0000-0003-1413-1776
    affiliation: 1
  - name: Jane Huang
    orcid: 0000-0001-6947-6072
    affiliation: "2, 3"
  - name: Feilong Meng
    orcid: 0000-0003-0079-6830
    affiliation: 2
affiliations:
 - name: Center for Astrophysics | Harvard & Smithsonian, 60 Garden St., Cambridge, MA 02138, USA
   index: 1
 - name: Department of Astronomy, University of Michigan, 323 West Hall, 1085 South University Avenue, Ann Arbor, MI 48109, USA
   index: 2
 - name: NASA Hubble Fellowship Program Sagan Fellow
   index: 3
date: 29 September 2021
bibliography: paper.bib

---

# Summary

 `disksurf` implements the method presented in @Pinte:2018 to extract the molecular emission surface (i.e., the height above the midplane from which molecular emission arises) in moderately inclined protoplanetary disks. The Python-based code leverages the open-source `GoFish` [@GoFish] package to read in and interact with FITS data cubes used for essentially all astronomical observations at submillimeter wavelengths. The code also uses the open-source `detect_peaks.py` routine from @Duarte:2021 for peak detection. For a given set of geometric parameters specified by the user, `disksurf` will return a surface object containing both the disk-centric coordinates of the surface as well as the gas temperature and rotation velocity at those locations. The user is able to 'filter' the returned surface using a variety of clipping and smoothing functions. Several simple analytical forms commonly adopted in the protoplanetary disk literature can then be fit to this surface using either a chi-squared minimization with `scipy` [@Virtanen:2020] or through an Monte-Carlo Markov-Chain approach with `emcee` [@Foreman-Mackey:2013]. To verify the 3D geometry of the system is well constrained, `disksurf` also provides diagnostic functions to plot the emission surface over channel maps of line emission (i.e., the emission morphology for a specific frequency).

# Statement of need

The Atacama Millimeter/submillimeter Array (ALMA) has brought our view of protoplanetary disks, the formation environment of planets, into sharp focus. The unparalleled angular resolution now achievable with ALMA allows us to routinely resolve the 3D structure of these disks; detailing the vertical structure of the gas and dust from which planets are formed. Extracting the precise height from where emission arises is a key step towards understanding the conditions in which a planet is born, and, in turn, how the planet can affect the parental disk.

A method for extracting a 'scattering surface', the emission surface equivalent for small, submicron grains was described in @Stolker:2016 who provided the `diskmap` package. However, this approach is not suitable for molecular emission, which traces the gas component of the disk and has a strong frequency dependence due to Doppler shifting from the disk rotation. @Pinte:2018 presented an alternative method that could account for this frequency dependence and demonstrated that this could be used to trace key physical properties of the protoplanetary disk, namely the gas temperature and rotation velocity, along the emission surface. An example script was provided with this publication demonstrating the algorithm: https://github.com/cpinte/CO_layers.

While the measurement of the emission surface only requires simple geometrical transformations, the largest source of uncertainty arises through the handling of noisy data. As more works perform such analyses, for example @Teague:2019, @Rich:2021, and @Law:2021, the need for an easy-to-use package that implements this method was clear. Such a package would facilitate the rapid reproduction of published results, enable direct comparisons between numerical simulations and observations [@Calahan:2021; @Schwarz:2021], and ease benchmarking between different publications. `disksurf` provides this functionality, along with a tutorial to guide users through the process of extracting an emission surface. The code is developed in such a way that as the quality of observations improve, the extraction methods can be easily refined to maintain precise measurements of the emission surface.

# Acknowledgements

We acknowledge help from Christophe Pinte in benchmarking early versions of the code with those presented in the original paper detailing the method, @Pinte:2018. R.T. acknowledges support from the Smithsonian Institution as a Submillimeter Array (SMA) Fellow. C.J.L. acknowledges funding from the National Science Foundation Graduate Research Fellowship under Grant DGE1745303. Support for J.H. was provided by NASA through the NASA Hubble Fellowship grant #HST-HF2-51460.001- A awarded by the Space Telescope Science Institute, which is operated by the Association of Universities for Research in Astronomy, Inc., for NASA, under contract NAS5-26555.

# References
