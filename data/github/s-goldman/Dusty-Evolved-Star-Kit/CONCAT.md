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
  url = {https://joss.theoj.org/papers/10.21105/joss.02554},
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
.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/s-goldman/Dusty-Evolved-Star-Kit/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Dusty-Evolved-Star-Kit could always use more documentation, whether as part of the
official Dusty-Evolved-Star-Kit docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/s-goldman/Dusty-Evolved-Star-Kit/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `dusty_evolved_star_kit` for local development.

1. Fork the `dusty_evolved_star_kit` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/dusty_evolved_star_kit.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv dusty_evolved_star_kit
    $ cd dusty_evolved_star_kit/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass the
   tests, including testing other Python versions and operating systems::

    $ pytest


6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 2.7, 3.4, 3.5 and 3.6, and for PyPy. Check
   https://travis-ci.org/github/s-goldman/Dusty-Evolved-Star-Kit/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run the tests use the command::

    $ pytest
=======
History
=======

1.7.2 (2020-08-05)

Updates plots with now 17 source maximum for single SED figure. The sed and
sed_indivscripts now use a common set of functions to minimize duplicated code,
and ensure consistency. 


1.7.1 (2020-08-05)

Minor bug fixes to the data retrieval script.


1.7.0 (2020-07-24)

First major release of least-squares fitting DESK. The package is stable with
pytest testing through Github-actions (Ubuntu and OSX: Python 3.6, 3.7, 3.8),
documentation with Sphinx on readthedocs, coverage with codecov,
code quality checks with codacy, installation with pip, and hosted on Github.
Model grids are downloaded using direct-download links on Box.

1.3.1 (2019-04-29)
------------------

* First release on PyPI.
======================
Dusty-Evolved-Star-Kit
======================




The DESK is an SED-fitting python package for fitting data from evolved stars
(photometry or spectra) with radiative transfer model grids. The package
is currently in development and all contributions are welcomed. For current
progress, see the 'Issues' tab on the Github_ page. The package is ideal for
fitting small samples of dusty evolved stars. It will soon utilize a
bayesian-fitting strategy with mass-loss rate and luminosity distributions as
inputs (priors), and will provide a better fit to these broader sample
properties.

**Input**: A csv file with the first column as wavelength in um and second column
as flux in Jy. To fit multiple csv files, put them in a directory, and use the
directory name as the input.

**Output**: A csv file including the best fit model and corresponding
stellar parameters, as well as an optional figure of the fit SED.

Available model grids: Several grids are already available upon installation. A range of
other model grids, including state-of-the-art dust-growth models by `Nanni et al. (2019)`_
, are downloaded automatically and used when selected. Descriptions of the model
grids can be found in the Documentation_. A module for creating your own DUSTY_ grid
is under development, but for now, please email me (`Dr. Steven Goldman`_) directly
for potential grid requests or for help with the package.

* Free software: BSD license
* Documentation: https://dusty-evolved-star-kit.readthedocs.io.


* TODO

**Credits**: This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Github: https://github.com/s-goldman/Dusty-Evolved-Star-Kit
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _DUSTY : https://github.com/ivezic/dusty
.. _Documentation : https://dusty-evolved-star-kit.readthedocs.io/en/latest/grids.html
.. _Nanni et al. (2019) : https://ui.adsabs.harvard.edu/abs/ 2019MNRAS.487..502N/abstract
.. _GRAMS : https://2dust.stsci.edu/grams_models.cgi
.. _2DUST : https://2dust.stsci.edu/index.cgi
.. _`Dr. Steven Goldman` : http://www.stsci.edu/~sgoldman/
.. include:: ../CONTRIBUTING.rst
.. include:: ../HISTORY.rst
=====
Usage
=====

SED-Fitting
-----------

The Dusty-Evolved-Star-Kit can fit the spectra and photometry of evolved stars
with both carbon- and oxygen-rich grids of radiative transfer models.
This enables easy comparison of grids of models and different samples.

The input for this package is a csv file with wavelength in microns in the first
column, and flux in Jy in the second column.

After installation with pip, in any command line prompt, use the command:

.. code-block:: console

	> desk grids

This will display the available grids for fitting. Next you need to point the
package to your csv file, and specify the distance and grid of choice:

.. code-block:: console

	> desk fit --source='target_name.csv' --distance=50 --grid='H11-LMC'

or specify a directory with multiple csv files:

.. code-block:: console

	> desk fit --source='folder_of_csvs' --distance=30 --grid='silicates'


Additional fitting options
--------------------------

Users can also specify any-and-all of the following additional options. Requests
for more additional features can be submitted through the `Issues`_ tab on the
`Github`_ page.

Grid density
============

.. code-block:: console

	> desk fit --source='target_name.csv' --n=200

The strength of the 1D
DUSTY models is that they can be scaled to create more luminous models. The DESK
takes the initial grid and scales it n times (default: 50) to create a larger
denser grid of sources within the luinosity limits (default: 1,000 - 150,000 Msun).
As the shape of the model is the same for each scaling of the grid, the distance to
the source(s) is important for an accurate results. A more realistic Bayesian method
of fitting is under development.

Wavelength range
================
.. code-block:: console

	> desk fit --source='target_name.csv' --min_wavelength=0.1 --max_wavelength=30

The user may also specify a wavelength minimum and maximum. This will still show
the full photometry in the final SED figure, but fit only the wavelength range
specified.


Multiprocessing
===============
The user can specify whether to fit using multiprocessing
(using all but 1 computer cores), single core fitting (multiprocessing=False), or
specify the number of cores to use (multiprocessing=6).
Multiprocessing uses a core per source, and will have little affect on small samples
or individual sources:

.. code-block:: console

	> desk fit --source='target_name.csv' --multiprocessing=True

Output Figures and Model Spectra
--------------------------------

For the DUSTY grids, the code automatically outputs the best-fit model spectrum as
a csv file named as the target name, grid name, luminosity, effective temperature,
inner dust temperature, optical depth, and distance, all separated by underscores.
You also have the option of creating a figure with the data and best-fit SED using
the following commands.

.. image:: ./example.png
	:width: 400
	:alt: SED example

This is an example of the output_sed.png file fitting three massive oxygen-rich
AGB stars from the LMC created using

.. code-block:: console

	> desk sed

To produce individual figures subsequently run the command:

.. code-block:: console

	> desk sed_indiv

Additionally you can specify whether you want the output flux in the figure to
be in W/m2 or Jy (W/m2 is the default).

.. code-block:: console

	> desk sed --flux='Jy'


Retrieve photometry
-------------------
Users can retrieve all of the photometry hosted on `Vizier`_ for a given source name
or coordinates. Retrieving photometry using a source name is as simple as:

.. code-block:: console

	> desk vizier_sed 'MSX LMC 807'

In order to return photometry using a source position (RA and Decl. in degrees), use
the command:

.. code-block:: console

	> desk vizier_sed '(83.15482600, -67.11567600)'

Additionally, you can specify the radius (in arcseconds) you would like to
search for photometry using. To specify a 5 arcsecond radius use:

.. code-block:: console

	> desk vizier_sed 'MSX LMC 807' --r=5

or

.. code-block:: console

	> desk vizier_sed '(083.15482600, -67.11567600)' --r=5

Use in Python Environment
-------------------------

SED-fitting can be done with the DESK within the python environment. To do this
simply import the package and use the 'fit' function in a similar manner as the
console commands.


.. code-block:: console

	>>> from desk import *
	>>> fit(source="target.csv", distance=3, grid="oxygen")

One can also use the sed, vizier_sed, and grids in a similar fashion.

.. code-block:: console

	>>> sed()
	>>> sed(flux='Jy')
	>>> grids()
	>>> vizier_sed('MSX LMC 807', 5)

How reliable in SED-fitting?
----------------------------
The DESK is a tool designed to allow for the easy comparison of samples and model
grids. Taken at face value, the results for a given sample or model grid may give
incorrect results. For example, recent work by `Wiegert et al. 2019`_ has shown
that the assumed geometry can affect measured mass loss rates by several orders
of magnitude. It is up to the user to interpret the results, and I would urge those
interested in using the DESK to also take a look at the excellent `recent review`_ by Leen Decin.


Using Multi-epoch data
-----------------------
The continuum shape of an SED is very useful in constraining values like luminosity and mass-loss rate. For variable evolved stars, however, fluxes can change by orders of magnitude on scales of 200-2000 days. Data taken at different times can alter the observed shape dramtically, and thus using multi-epoch data is discouraged. If a user has a large sample with data in many overlapping bands, the DESK can attempt to fit the median SEDs giving an idea of the properties of the sample as a whole.


Package Testing
---------------
The desk uses continuous integration testing through Github actions. This
automatically runs the package tests for several commonly used operating systems
and python versions, before every change that is made to the code.
The current status of the `tests`_ and `coverage`_.
are available online. To run the tests locally, download/clone the package and
use the command 'pytest' within the pacakge directory.

.. _Vizier: http://vizier.cfa.harvard.edu/
.. _github: https://github.com/s-goldman/Dusty-Evolved-Star-Kit/
.. _Issues: https://github.com/s-goldman/Dusty-Evolved-Star-Kit/issues
.. _tests: https://github.com/s-goldman/Dusty-Evolved-Star-Kit/actions?query=workflow%3A%22Python+package%22
.. _coverage: https://codecov.io/gh/s-goldman/Dusty-Evolved-Star-Kit
.. _recent review: https://ui.adsabs.harvard.edu/abs/2020arXiv201113472D/abstract
.. _Wiegert et al. 2019: https://ui.adsabs.harvard.edu/abs/2020A%26A...642A.142W/abstract
desk.fitting package
====================

Submodules
----------

desk.fitting.dusty\_fit module
------------------------------

.. automodule:: desk.fitting.dusty_fit
   :members:
   :undoc-members:
   :show-inheritance:

desk.fitting.fitting\_tools module
----------------------------------

.. automodule:: desk.fitting.fitting_tools
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: desk.fitting
   :members:
   :undoc-members:
   :show-inheritance:
.. include:: ../README.rst
.. highlight:: shell

============
Installation
============


Stable release
--------------

To install Dusty-Evolved-Star-Kit, run this command in your terminal:

.. code-block:: console

    $ pip install desk

This is the preferred method to install Dusty-Evolved-Star-Kit, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for Dusty-Evolved-Star-Kit can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/s-goldman/Dusty-Evolved-Star-Kit

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/s-goldman/Dusty-Evolved-Star-Kit


Dependecies
-----------
astropy

numpy

ipdb

tqdm

matplotlib

h5py

wheel

twine

sphinx

sphinx_automodapi

pytest-cov

seaborn
desk
====

.. toctree::
   :maxdepth: 4

   desk
desk.set\_up package
====================

Submodules
----------

desk.set\_up.config module
--------------------------

.. automodule:: desk.set_up.config
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.create\_output\_files module
-----------------------------------------

.. automodule:: desk.set_up.create_output_files
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.error\_messages module
-----------------------------------

.. automodule:: desk.set_up.error_messages
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.full\_grid module
------------------------------

.. automodule:: desk.set_up.full_grid
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.get\_data module
-----------------------------

.. automodule:: desk.set_up.get_data
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.get\_inputs module
-------------------------------

.. automodule:: desk.set_up.get_inputs
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.get\_models module
-------------------------------

.. automodule:: desk.set_up.get_models
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.scale\_dusty module
--------------------------------

.. automodule:: desk.set_up.scale_dusty
   :members:
   :undoc-members:
   :show-inheritance:

desk.set\_up.scale\_external module
-----------------------------------

.. automodule:: desk.set_up.scale_external
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: desk.set_up
   :members:
   :undoc-members:
   :show-inheritance:
The Dusty-Evolved-Star-Kit
======================================

.. toctree::
   :maxdepth: 2
   :caption: Installation:

   readme
   installation


.. toctree::
   :maxdepth: 2
   :caption: User Documentation:

   usage
   grids
   modules
   contributing
   history

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
desk package
============

Subpackages
-----------

.. toctree::

   desk.fitting
   desk.outputs
   desk.probabilities
   desk.set_up

Submodules
----------

desk.console\_commands module
-----------------------------

.. automodule:: desk.console_commands
   :members:
   :undoc-members:
   :show-inheritance:

desk.main module
----------------

.. automodule:: desk.main
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: desk
   :members:
   :undoc-members:
   :show-inheritance:
desk.outputs package
====================

Submodules
----------

desk.outputs.interpolate\_dusty module
--------------------------------------

.. automodule:: desk.outputs.interpolate_dusty
   :members:
   :undoc-members:
   :show-inheritance:

desk.outputs.parameter\_ranges module
-------------------------------------

.. automodule:: desk.outputs.parameter_ranges
   :members:
   :undoc-members:
   :show-inheritance:

desk.outputs.plot\_pdf module
-----------------------------

.. automodule:: desk.outputs.plot_pdf
   :members:
   :undoc-members:
   :show-inheritance:

desk.outputs.plotting\_seds module
----------------------------------

.. automodule:: desk.outputs.plotting_seds
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: desk.outputs
   :members:
   :undoc-members:
   :show-inheritance:
desk.probabilities package
==========================

Submodules
----------

desk.probabilities.compute\_grid\_weights module
------------------------------------------------

.. automodule:: desk.probabilities.compute_grid_weights
   :members:
   :undoc-members:
   :show-inheritance:

desk.probabilities.create\_pdf module
-------------------------------------

.. automodule:: desk.probabilities.create_pdf
   :members:
   :undoc-members:
   :show-inheritance:

desk.probabilities.create\_prior module
---------------------------------------

.. automodule:: desk.probabilities.create_prior
   :members:
   :undoc-members:
   :show-inheritance:

desk.probabilities.resample\_prior\_to\_model\_grid module
----------------------------------------------------------

.. automodule:: desk.probabilities.resample_prior_to_model_grid
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: desk.probabilities
   :members:
   :undoc-members:
   :show-inheritance:
=====
Grids
=====

The DESK_ downloads and fits both new and commonly-used model grids, automatically.
The grids are also available for download on Zenodo_. Columns are in wavelength (um) and flux density (W*m-2).

DUSTY model grids
-----------------

The DUSTY_ models are 1-D radiative transfer models that exploit
scaling relations to generate models. The model shape can be scaled
to any luminosity, expansion velocity, or mass-loss rate so long as the
combination of these three parameters remains the same. The DESK model grids
contain 2,000--32,000 unique model shapes. Theses grids are scaled
to the desired number of luminosities between 1,000 and 150,000 using the
command option `--n`.

.. code-block:: console

	> desk fit --source='example.csv' --grid='corundum-mix' --n=10

This command will scale the 'corundum-mix' model grid of 3,000 unique model
shapes to luminosities at 1000, 17555, 34111, 50666, 67222, 83777, 100333, 116888,
133444, and 150000 solar luminosities (30,000 models).

The output model parameters of expansion velocity and gas mass-loss rate are
scaled using the relations from `Elitzur & Ivezić 2001`_. The model flux is scaled
using the brightness and distance to the sun, and scaling for the luminosity
which is already in solar luminosities.

Only the default carbon (amorphous-carbon) and oxygen-rich (silicates) model grids
are downloaded during installation. Additional grids are downloaded when selected.
The number of unique models and the file size of each grid are specified next to
the grid names below.

The models assume a standard `MRN`_ grain size distribution from
0.005 - 0.25 microns. The DUSTY code uses an analytical approximation of the wind density
distribution. This calculation extends to a distance of 10,000 times the inner
radius (density type = 4, within the DUSTY code).


Oxygen-rich DUSTY model grids
=============================

silicates (*N*\ =2,000; 3.4 MB): Uses oxygen-rich warm silicates from
`Ossenkopf et al. (1992)`_ and a black body. This grid provides ranges in
effective temperature (2600-3400 K: 200 K interval) inner dust
temperature (600-1200K: 200 K interval) and optical depth (0.1 - 30: 100
spaced logarithmically).

o-def-silicates (*N*\ =2,000; 3.4 MB): Same as the silicates grid but using
oxygen-deficient warm silicates from `Ossenkopf et al. (1992)`_.


corundum-mix (*N*\ =3,000; 5 MB): Same temperature range as the silicates grid
but with with more limited density of optical depths (0.1 - 30: 50 spaced
logarithmically) created for three scenarios of different fractions of warm
silicates and corundum grains from `Begemann et al. (1997)`_. The fractions
included are 70:30, 80:20, 90:10 for the warm silicates:corundum grain ratio.


crystalline-mix (*N*\ =3,000; 5 MB): Same as corundum-mix but with crystalline silicate
grains from from `Jaeger et al. (1994)`_.


draine-mix (*N*\ =3,000; 5 MB): Same as corundum-mix but with astronomical silicate
grains from from `Draine & Lee (1984)`_.


iron-mix (*N*\ =2,960; 5 MB): Same as corundum-mix but with iron grains from
`Henning et al. (1995)`_,.


carbon-mix (*N*\ =3,000; 5 MB): Same as corundum-mix but with amorphous carbon
grains from `Zubko et al. (1996)`_.


silicate-mix (*N*\ =71,435; 127 MB): This grid provides ranges in
effective temperature (2600-3400 K: 200 K interval) inner dust
temperature (600-1200K: 100 K interval) and optical depth (0.1 - 30: 100
spaced linearally). We set 4% iron grains from `Henning et al. (1995)`_, and different
fractions of the oxygen-rich and
oxygen-deficient grains from `Ossenkopf et al. (1992)`_ and crystalline silicates from
`Jaeger et al. (1994)`_. The fractions for each grain type are 0, 16, 24, 38, 42,
56, 64, 80% for the oxygen-rich and oxygen-deficient silicates with the remainder filled with
the crystalline silicates and 4% iron grains.

Explore this model grid using the new `interactive notebook`_


Carbon-rich DUSTY model grids
=============================

amorphous-carbon (*N*\ =1,340; 145 MB): Same as silicates but with
amorphous carbon grains from `Zubko et al. (1996)`_.


.. _the-dust-growth-model-grids-from-nanni-et-al-2019:

Dust growth grids
=================

The dust growth model grids from `Nanni et al. (2019)`_.

H11-LMC (*N*\ =90,899; 770 MB): A carbon-rich grid for the LMC metallicity (1/2
solar) using optical constants from `Hanner et al. (1988)`_.

H11-SMC (*N*\ =91,058; 772 MB): A carbon-rich grid for the SMC metallicity (1/5
solar) using optical constants from `Hanner et al. (1988)`_.

J1000-LMC (*N*\ =85,392; 723 MB): A carbon-rich grid for the LMC metallicity
(1/2 solar) using optical constants from `Jaeger et al. (1998)`_

J1000-SMC (*N*\ =85,546; 724 MB): A carbon-rich grid for the SMC metallicity
(1/5 solar) using optical constants from `Jaeger et al. (1998)`_


2-D model grids
-------------------------


The GRAMS model grids
=====================

The GRAMS model grids from `Sargent et al. (2011)`_ and `Srinivasan et al. (2011)`_.
These models assume a density distribution of 1/r\ :sup:`2`, a modified KMH dust grain
distribution, and an assumed expansion velocity of 10 km/s. Compared to the 1D DUSTY models,
the optical depths for the GRAMS models are more limited, with  optical depths at 10 microns
ranging from 0.001-0.4 for the carbon-rich grid, to 0.0001-26 for the oxygen-rich grid.


grams-carbon (*N*\ =12,244; 41.3 MB): A 2D carbon-rich grid using the `2DUST`_
code for the LMC metallicity (1/2 solar) using optical constants from
`Zubko et al. (1996)`_.

grams-oxygen (*N*\ =68,601; 200 MB): A 2D oxygen-rich grid using the `2DUST`_
code for the LMC metallicity (1/2 solar) using optical constants from
`Ossenkopf et al. (1992)`_.

grams: Both GRAMS datasets combined

.. code:: diff

   - Warning: results uncertain outside of a distance 20-150 kpc.

.. _DESK: https://github.com/s-goldman/Dusty-Evolved-Star-Kit
.. _Zenodo: https://zenodo.org/record/5574616
.. _DUSTY: https://github.com/ivezic/dusty
.. _Elitzur & Ivezić 2001: https://ui.adsabs.harvard.edu/abs/2001MNRAS.327..403E/abstract
.. _Sargent et al. (2011): https://ui.adsabs.harvard.edu/abs/2011ApJ...728...93S/abstract
.. _Srinivasan et al. (2011): https://ui.adsabs.harvard.edu/abs/2011A%26A...532A..54S/abstract
.. _2DUST: https://ui.adsabs.harvard.edu/abs/2003ApJ...586.1338U/abstract
.. _Zubko et al. (1996): https://ui.adsabs.harvard.edu/abs/1996MNRAS.282.1321Z/abstract
.. _Ossenkopf et al. (1992): https://ui.adsabs.harvard.edu/abs/1992A%26A...261..567O/abstract
.. _Aringer et al. (2016): https://ui.adsabs.harvard.edu/abs/2016MNRAS.457.3611A/abstract
.. _MRN: https://ui.adsabs.harvard.edu/abs/1977ApJ...217..425M/abstract
.. _Jaeger et al. (1994): https://ui.adsabs.harvard.edu/abs/1994A%26A...292..641J/abstract
.. _Jaeger et al. (1998): https://ui.adsabs.harvard.edu/abs/1998A%26A...339..904J/abstract
.. _Begemann et al. (1997): https://ui.adsabs.harvard.edu/abs/1997ApJ...476..199B/abstract
.. _Henning et al. (1995): https://ui.adsabs.harvard.edu/abs/1995A%26AS..112..143H/abstract
.. _Zubko et al. (1996): https://ui.adsabs.harvard.edu/abs/1996MNRAS.282.1321Z/abstract
.. _Nanni et al. (2019): https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..502N/abstract
.. _Hanner et al. (1988): https://ui.adsabs.harvard.edu/abs/1988ioch.rept.....H/abstract
.. _Groenewegen 2012: https://ui.adsabs.harvard.edu/abs/2012A&A...543A..36G/abstract
.. _Dorschner et al. (1995): https://ui.adsabs.harvard.edu/abs/1995A&A...300..503D/abstract
.. _Gustafsson et al. (2008): https://ui.adsabs.harvard.edu/abs/2008A%26A...486..951G/abstract
.. _Draine & Lee (1984): https://ui.adsabs.harvard.edu/abs/1984ApJ...285...89D/abstract
.. _interactive notebook: https://mybinder.org/v2/gh/s-goldman/Dusty-Evolved-Star-Kit_notebooks/main?labpath=silicate-mix_interactive.ipynb
