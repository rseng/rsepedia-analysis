---
title: 'The Dusty Evolved Star Kit (DESK): A Python package for fitting the Spectral Energy Distribution of Evolved Stars'
tags:
  - Python
  - astronomy
  - asymptotic giant branch stars
  - radiative Transfer
  - stellar Mass Loss
  - spectral energy distribution fitting
authors:
  - name: Steven R. Goldman
    orcid: 0000-0002-8937-3844
    affiliation: 1
affiliations:
 - name: Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
   index: 1
date: 15 July 2020
bibliography: paper.bib
---

# Summary

One of the few ways that we can understand the environment around dusty stars and how much material they contribute back to the Universe, is by fitting their brightness at different wavelengths with models that account for how the energy transfers through the dust. Codes for creating models have been developed and refined [@Elitzur:2001; @Ueta:2003], but a code for easily fitting data to grids of realistic models has been up-to-this-point unavailable.


The ``DESK`` is a python package designed to compare the best fits of different stellar samples and model grids for a better understanding of the results and their uncertainties. The package fits the Spectral Energy Distribution (SED) of evolved stars, using photometry or spectra, to grids of radiative transfer models using a least-squares method. The package includes newly created grids using a variety of different dust species, and state-of-the-art dust growth grids [@Nanni:2019]. Early versions of the code have been used in [@Orosz:2017; @Goldman:2017; @Goldman:2018; @Goldman:2019b]

# Statement of need

To understand the ranges and estimated errors of fitted results, they must be compared to results from different model grids. Results from these grids (e.g. luminosity, mass-loss rate) can vary dramatically as a result of the unknown dust properties and geometry of evolved stars [@Sargent:2010; @Srinivasan:2011; @Wiegert:2020]. This is especially true of the oxygen-rich Asymptotic Giant Branch (AGB) stars. Adding to this challenge is the fact that models are calculated based on measured values of the dust (optical constants) which can not be interpolated over. A robust method for testing different model grids will be particularly important given the wealth of infrared data to come from the James Webb Space Telescope (JWST).


# User interface

The package can be installed using `pip` and imported within python. Using "entrypoints", the package can also be accessed from any terminal prompt once installed. The fitting method uses a brute-force technique to ensure a true best fit. New grids of multi-dimensional radiative transfer models will be added to the model grid library as they are developed. The available model grids for this version are listed in Table 1.

# Figures

![An example of three massive oxygen-rich AGB stars in the Large Magellanic Cloud (LMC) galaxy fit with the default oxygen-rich model grid (Oss-Orich-bb). These three example sources can be fit, and this figure can be created, using the command `desk fit` and then the command `desk sed`.  ](docs/example.png)

![](docs/paper/joss_table.png)

---
nocite: |
  @Aringer:2016, @Begemann:1997, @Henning:1995, @Jaeger:1998, @Ossenkopf:1992, @Zubko:1996
---

# References
Dusty-Evolved-Star-Kit<img align="left" width="100" height="100" src="docs/the_desk.png">
=========================================================================================
[![Build](https://github.com/s-goldman/Dusty-Evolved-Star-Kit/workflows/Python%20package/badge.svg?branch=master)](https://github.com/s-goldman/Dusty-Evolved-Star-Kit/actions)
[![Documentation Status](https://readthedocs.org/projects/dusty-evolved-star-kit/badge/?version=latest)](https://dusty-evolved-star-kit.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b6bd41e6d7db48e7b811a106015f2d82)](https://www.codacy.com/manual/s-goldman/Dusty-Evolved-Star-Kit?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=s-goldman/Dusty-Evolved-Star-Kit&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/s-goldman/Dusty-Evolved-Star-Kit/branch/master/graph/badge.svg)](https://codecov.io/gh/s-goldman/Dusty-Evolved-Star-Kit)
[![pypi](https://img.shields.io/badge/pypi-DESK-blue.svg)](https://pypi.org/project/desk/)
[![status](https://joss.theoj.org/papers/b78c206113fdb59a7a8839649786e9d8/status.svg)](https://joss.theoj.org/papers/b78c206113fdb59a7a8839649786e9d8)
[![Steve Goldman](https://img.shields.io/badge/STScI-Steve%20Goldman-blue.svg)](http://www.stsci.edu/~sgoldman/)

The DESK is an SED-fitting python package for fitting data from evolved stars (photometry or spectra) with radiative transfer model grids. The package is currently in development and all contributions are welcomed. For current progress, see the Issues tab at the top of the page. The package is ideal for fitting small samples of dusty evolved stars. It will soon utilize a bayesian-fitting strategy with mass-loss rate and luminosity distributions as inputs (priors), and will provide a better fit  to these broader sample properties.

**Input**: A csv file with the first column as wavelength in um and second column as flux in Jy. To fit multiple csv files, put them in a directory, and use the directory name as the input.

**Output**: A csv files with the best fit model and corresponding stellar parameters, as well as an optional figure of the fit SED.

**Available model grids**:
Several grids are **already available** upon installation. A range of other model grids, including 2D [GRAMS](https://ui.adsabs.harvard.edu/abs/2011ApJ...728...93S/abstract) model grids based on the [2DUST](https://ui.adsabs.harvard.edu/abs/2003ApJ...586.1338U/abstract) code, and state-of-the-art dust-growth models by [Nanni et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..502N/abstract), are downloaded automatically and used when selected. Descriptions of the [model grids](https://dusty-evolved-star-kit.readthedocs.io/en/latest/grids.html) can be found in the documentation.


A module for creating your own [DUSTY](https://github.com/ivezic/dusty) grid is under development, but for now, please email me ([Dr. Steven Goldman](http://www.stsci.edu/~sgoldman/)) directly for potential grid requests or for help with the package.

Documentation
-------------

The documentation can be found on [readthedocs](http://dusty-evolved-star-kit.readthedocs.io/en/latest/).

Install Using Python
--------------------

1). Install the package from [source](https://dusty-evolved-star-kit.readthedocs.io/en/latest/installation.html) or with [pip](https://pypi.org/project/pip/) using the command `pip install desk`.

![](docs/pip_install2.gif)

Using the DESK
--------------

2). Go to the directory where your target csv file (or target directory of files) is.  

3). Use the command (without starting python)

  `desk fit --source='target_name.csv'`

or if you have a folder of csv files

  `desk fit --source='folder_of_csvs'`

To fit the example sources use the command

  `desk fit`

additional options are:

`desk fit --source='target_name.csv' --distance=50 --grid='Oss-Orich-bb'`

The other important options are the distance (in kpc) and the grid of models you would like to use (options listed below). For other options see the [Usage](https://dusty-evolved-star-kit.readthedocs.io/en/latest/usage.html) page. For the model grids, you can select 'oxygen' or 'carbon' to use the default models. To see other available grids use:

`desk grids`

To create a figure showing all of the fits of the SED, use the following command in the same directory.

`desk sed`


This is an example of the output_sed.png file fitting three massive oxygen-rich AGB stars from the LMC.

<img src="docs/example.png"  width="400" height="500">

To produce individual figures for each SED instead use the command:

`desk sed_indiv`

The package can also be used within python (see the [docs](https://dusty-evolved-star-kit.readthedocs.io/en/latest/usage.html#use-in-python-environment)).

Retrieve Photometry
-------------------

Don't have the photometry? You can retrieve them from *Vizier* using the [vizier_sed](https://dusty-evolved-star-kit.readthedocs.io/en/latest/usage.html#use-in-python-environment) command if you have the source name or position in degrees:

`desk vizier_sed 'MSX LMC 807'`

or

`desk vizier_sed '(83.15482600, -67.11567600)'`


Afterwords, you can fit that data with the command:

`desk fit --source='MSX_LMC_807_sed.csv'`



Citation
-----------

Goldman, S. R. 2020, Journal of Open Source Software, 5, 2554, doi: 10.21105/joss.02554

or with bibtex:

<pre>@article{Goldman2020,
  doi = {10.21105/joss.02554},
  url = {https://doi.org/10.21105/joss.02554},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2554},
  author = {Steven R. Goldman},
  title = {The Dusty Evolved Star Kit (DESK): A Python package for fitting the Spectral Energy Distribution of Evolved Stars},
  journal = {Journal of Open Source Software}
}</pre>

Please also specify the options selected and make the data publicly available for reproducibility.

License
-------

This project is Copyright (c) [Dr. Steven Goldman](http://www.stsci.edu/~sgoldman/) and licensed under
the terms of the BSD 3-Clause license.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
 - OS: [e.g. iOS8.1]
 - Version [e.g. 22]
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
