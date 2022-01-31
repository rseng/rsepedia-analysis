# Tools21cm

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02363/status.svg)](https://doi.org/10.21105/joss.02363)

A python package for analysing 21-cm signals from the Epoch of Reionization (EoR) and Cosmic Dawn (CD). Full documentation (with examples, installation instructions and complete module description) can be found at [readthedocs](https://tools21cm.readthedocs.io/).

## Package details

The package provides tools to analyse cosmological simulations of EoR and CD. It contains modules to create mock 21-cm observations for current and upcoming radio telescopes, such as LOFAR, MWA and SKA, and to construct statistical measures.

### Input

Currently, `Tools21cm` supports the following simulation codes:

* [CUBEP3M](https://github.com/jharno/cubep3m)
* [C2RAY](https://github.com/garrelt/C2-Ray3Dm)
* [GRIZZLY](https://arxiv.org/abs/1710.09397)
* [21cmFAST](https://21cmfast.readthedocs.io/en/latest/)
* [Simfast21](https://github.com/mariogrs/Simfast21)
* [sem_num](https://arxiv.org/abs/1403.0941)

### Outputs

There are various manipulation and analysis moduled in `Tools21cm`. 

* Angular coordinates: methods to convert data between physical (cMpc) coordinates and observational (angular-frequency) coordinates
* Bubble Statistics: methods to calcluate the sizes of the regions of interest and estimate the size distributions
* Cosmological calculators: various functions for calculating some cosmological stuff
* Identifying regions: methods to identify regions of interest in images
* Lightcone: methods to construct lightcones
* Power spectrum: contains functions to estimate various two point statistics
* Reading simuation outputs: methods to read simulations outputs
* Smoothing: methods to smooth or reduce resolution of the data to reduce noise
* Point statistics: contains various useful statistical methods
* Temperature: methods to estimate the brightness temperature

For detail documentation and how to use them, see [here](https://tools21cm.readthedocs.io/contents.html).

### Under Developement

* Foreground model: methods to simulate and analyse the foreground signal for 21 cm signal
* Primary beam: methods to simulate the primary beam of radio telescope
* Telescope noise: 
	* simulate the radio telescope observation strategy
	* simulate telescope noise
* Topology: methods to estimate the topology of the region of interest



## INSTALLATION

To install the package from source, one should clone this package running the following::

    git clone https://github.com/sambit-giri/tools21cm.git

To install the package in the standard location, run the following in the root directory::

    python setup.py install

In order to install it in a separate directory::

    python setup.py install --home=directory

One can also install it using pip by running the following command::

    pip install git+https://github.com/sambit-giri/tools21cm.git

The dependencies should be installed automatically during the installation process. If they fail for some reason, you can install them manually before installing tools21cm. The list of required packages can be found in the requirements.txt file present in the root directory.

### Tests

For testing, one can use [pytest](https://docs.pytest.org/en/stable/) or [nosetests](https://nose.readthedocs.io/en/latest/). Both packages can be installed using pip. To run all the test script, run the either of the following::

    python -m pytest tests
    
	nosetests -v

## CONTRIBUTING

If you find any bugs or unexpected behavior in the code, please feel free to open a [Github issue](https://github.com/sambit-giri/tools21cm/issues). The issue page is also good if you seek help or have suggestions for us. For more details, please see [here](https://tools21cm.readthedocs.io/contributing.html).
---
title: 'Tools21cm: A python package to analyse the large-scale 21-cm signal from the Epoch of Reionization and Cosmic Dawn'
tags:
  - Python
  - astronomy
  - early universe
  - 21-cm signal
  - reionization
authors:
  - name: Sambit K. Giri
    orcid: 0000-0002-2560-536X
    affiliation: "1, 2" 
  - name: Garrelt Mellema
    orcid: 0000-0002-2512-6748
    affiliation: 1
  - name: Hannes Jensen
    affiliation: 3
affiliations:
 - name: Department of Astronomy and Oskar Klein Centre, Stockholm University, AlbaNova, SE-106 91 Stockholm, Sweden
   index: 1
 - name: Institute for Computational Science, University of Zurich, Winterthurerstrasse 190, CH-8057 Zurich, Switzerland
   index: 2
 - name: Self
   index: 3
date: 8 June 2020
bibliography: paper.bib

---

# Summary

The Cosmic Dawn (CD) and Epoch of Reionization (EoR) are among the least understood epochs from the history of our Universe. Large-scale cosmological simulations are useful to understand this era and help prepare for future observations. `Tools21cm` is a data manipulation and analysis package for studying such cosmological simulations. It is written in the python programming language and designed to be very user-friendly.

The 21-cm signal, produced by the spin-flip transition of neutral hydrogen, is a unique tracer of matter in large-scale structures present during the EoR and CD [e.g. @Furlanetto2006CosmologyUniverse; @Pritchard201221cmCentury]. This signal will be cosmologically redshifted and thus found at low radio frequencies.
Much of the functionality of `Tools21cm` is focused on the construction and analysis of mock 21-cm observations in the context of current and upcoming radio telescopes, such as the Low Frequency Array [LOFAR; @vanHaarlem2013LOFAR:ARray], the Murchison Widefield Array [MWA; @Tingay2013MWA; @Wayth2018MWAdesign], Hydrogen Epoch of Reionization Array [HERA; @deboer2017hydrogen] and the Square Kilometre Array [SKA; @Mellema2013ReionizationArray]. `Tools21cm` post-processes cosmological simulation data to create mock 21-cm observations.

Radio telescopes typically observe the redshifted 21-cm signal over a range of frequencies and therefore produce 21-cm images at different redshifts. A sequence of such images from different redshifts is known as a tomographic data set and is three dimensional. `Tools21cm` can construct such tomographic data sets from simulation snapshots [@Datta2012Light-coneSpectrum; @Giri2018BubbleTomography]. When constructing these data sets, it can also add the impact of peculiar velocities, leading to an effect known as redshift space distortions [e.g. @Jensen2013ProbingDistortions; @Jensen2016TheMeasurements; @Giri2018BubbleTomography]. See @giri2019tomographic for a detailed description of tomographic 21-cm data sets.

`Tools21cm` also includes tools to calculate a wide range of statistical quantities from simulated 21-cm data sets. These include one-point statistics, such as the global or sky-averaged signal as a function of frequency, as well the variance, skewness and kurtosis [e.g. @Ross2017SimulatingDawn; @Ross2019EvaluatingDawn]. It can also characterise the spatial fluctuations in the signal through spherically and cylindrically averaged power spectra [@Jensen2013ProbingDistortions; @Ross2017SimulatingDawn; @Giri2019NeutralTomography] and position dependent power spectra [@Giri2019Position-dependentReionization].
It also has the capability to find interesting features, such as ionized regions, in (tomographic) image data [@Giri2018OptimalObservations; @Giri2019NeutralTomography] and from these to derive statistical quantities, such as size distributions [@Giri2018BubbleTomography] and topological quantities such as the Euler characteristic [@Giri2019NeutralTomography]. Such statistical characterisations of the data are required when comparing observations with simulations using a Bayesian inference framework in order to derive constraints on model parameters [e.g. @Greig201521CMMC:Signal].


# Acknowledgements

We acknowledge the contributions and feedback from Keri Dixon, Hannah Ross, Raghunath Ghara, Ilian Iliev, Catherine Watkinson and Michele Bianco. We are also thankful for support by Swedish Research Council grant 2016-03581.

# References
# Paper submission

Tools21cm paper for submission to the [Journal of Open Source Science](https://joss.theoj.org/).

[Paper submission instructions](https://joss.readthedocs.io/en/latest/submitting.html) -
this link also includes formatting instructions.

[Review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html)============
Contributing
============

Thank you so much for using ``tools21cm``. Contributions are welcome and are greatly appreciated. 

Bug reports/Feedback/Questions
===============================================
It is very helpful to us when users report bugs or unexpected behaviour in the code. Please feel free to open an issue at the repository's `issue <https://github.com/sambit-giri/tools21cm/issues>`_ page. It is also the best place to seek any help or would like to give us any feedback. 

If you report an issue, please give the necessary details (python environment, code version, the unexpected output or error message) and a "Working Example" to reproduce the problem.

Adding new features
===================
If you feel that there are some features missing, we welcome you to open an issue and discuss it. We can plan how this missing feature can be included in the code. 

Setting up your development environment
---------------------------------------
In case you want to make substantial change to the code, you can fork the repository and clone it to your local machine::

       git clone https://github.com/<YOURUSERNAME>/tools21cm.git
       cd tools21cm
       git checkout -b <BRANCHNAME


where ``<YOURUSERNAME>`` and ``<BRANCHNAME>`` are your Github username and the name of new branch for local development respectively.
.. include:: ../CONTRIBUTING.rst=============
Other modules
=============

Pre-defined constants
---------------------
.. automodule:: t2c.const
    :members:

For CUBEP\ :sup:`3`\M and C\ :sup:`2`\RAY simulations
-----------------------------------------------------

Conversion factors
~~~~~~~~~~~~~~~~~~
.. automodule:: t2c.conv
    :members:

Density file
~~~~~~~~~~~~
.. automodule:: t2c.density_file
    :members:

Ionisation rate file
~~~~~~~~~~~~~~~~~~~~
.. automodule:: t2c.irate_file
    :members:


Ionisation fraction file
~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: t2c.xfrac_file
    :members:

=======
Authors
=======

* `Sambit Giri <https://sambit-giri.github.io/>`_ - `GitHub <https://github.com/sambit-giri>`_
* `Garrelt Mellema <https://www.su.se/english/profiles/gmell-1.184545>`_ - `GitHub <https://github.com/garrelt>`_
* Hannes Jensen - `GitHub <https://github.com/hjens>`_

Contributors
============

* `Hannah Ross <https://crd.lbl.gov/departments/computational-science/c3/c3-people/hannah-ross/>`_
* `Raghunath Ghara <https://scholar.google.com/citations?user=WmNdlCkAAAAJ&hl=en>`_
* Ilian Iliev
* Michele Bianco
=========
Tutorials
=========


The following examples will help you get started with ``tools21cm``:

.. toctree::
   :maxdepth: 2
   :caption: Examples:

   examples/tutorials
   examples/lightcone
=========
Tools21cm
=========

.. image:: https://joss.theoj.org/papers/10.21105/joss.02363/status.svg
   :target: https://doi.org/10.21105/joss.02363

A python package for analysing 21-cm signals from the Epoch of Reionization (EoR) and Cosmic Dawn (CD). The source files can be found `here <https://github.com/sambit-giri/tools21cm>`_.

Note: There are some modules in the package that are still under active development. Therefore please contact the authors if you get erronous results.


Package details
===============

The package provides tools to analyse cosmological simulations of EoR and CD. It contains modules to create mock 21-cm observations for current and upcoming radio telescopes, such as LOFAR, MWA and SKA, and to construct statistical measures.

Input
-----

Currently, `Tools21cm` supports the following simulation codes:

* |cubep3m|_: a high performance cosmological N-body code
* |c2ray|_: a numerical method for calculating 3D radiative transfer and for simulating the EoR and CD
* `GRIZZLY <https://arxiv.org/abs/1710.09397>`_: an EoR and CD simulation code based on 1D radiative transfer 
* `21cmFAST <https://21cmfast.readthedocs.io/en/latest/>`_: a semi-numerical cosmological simulation code for the radio 21cm signal
* `Simfast21 <https://github.com/mariogrs/Simfast21>`_: a semi-numerical cosmological simulation code for the radio 21cm signal
* `sem_num <https://arxiv.org/abs/1403.0941>`_: a simple set of codes to semi-numerically simulate HI maps during reionization


.. |c2ray| replace:: C\ :sup:`2`\RAY
.. _c2ray: https://github.com/garrelt/C2-Ray3Dm

.. |cubep3m| replace:: CUBEP\ :sup:`3`\M
.. _cubep3m: https://github.com/jharno/cubep3m

Outputs
-------

There are various manipulation and analysis moduled in `Tools21cm`. 

* Angular coordinates: methods to convert data between physical (cMpc) coordinates and observational (angular-frequency) coordinates

* Bubble Statistics: methods to calcluate the sizes of the regions of interest and estimate the size distributions

* Cosmological calculators: various functions for calculating some cosmological stuff

* Identifying regions: methods to identify regions of interest in images

* Lightcone: methods to construct lightcones

* Power spectrum: contains functions to estimate various two point statistics

* Reading simuation outputs: methods to read simulations outputs

* Smoothing: methods to smooth or reduce resolution of the data to reduce noise

* Point statistics: contains various useful statistical methods

* Temperature: methods to estimate the brightness temperature

For detail documentation and how to use them, see `here <https://tools21cm.readthedocs.io/contents.html>`_.

Under Developement
------------------

* Foreground model: methods to simulate and analyse the foreground signal for 21 cm signal
* Primary beam: methods to simulate the primary beam of radio telescope
* Telescope noise: 
	* simulate the radio telescope observation strategy
	* simulate telescope noise
* Topology: methods to estimate the topology of the region of interest

============
Installation
============

We highly recommend the use of a virtual environement. It helps to keep dependencies required by different projects separate. Few example tools to create virtual environments are `anaconda <https://www.anaconda.com/distribution/>`_, `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ and `venv <https://docs.python.org/3/library/venv.html>`_.

The dependencies should be installed automatically during the installation process. If they fail for some reason, you can install them manually before installing `tools21cm <https://github.com/sambit-giri/tools21cm>`_. The list of required packages can be found in the *requirements.txt* file present in the root directory.

For a standard non-editable installation use::

    pip install git+git://github.com/sambit-giri/tools21cm.git [--user]

The --user is optional and only required if you don't have write permission to your main python installation.
If you wants to work on the code, you can download it directly from the `GitHub <https://github.com/sambit-giri/tools21cm>`_ page or clone the project using::

    git clone git://github.com/sambit-giri/tools21cm.git

Then, you can just install in place without copying anything using::

    pip install -e /path/to/tools21cm [--user]

The package can also be installed using the *setup.py* script. Find the file *setup.py* in the root directory. To install in the standard directory, run::

    python setup.py install

If you do not have write permissions, or you want to install somewhere else, you can specify some other installation directory, for example::

    python setup.py install --home=~/mydir

To see more options, run::

    python setup.py --help-commands

Or look `here <http://docs.python.org/2/install/>`_ for more details.

Tests
-----
For testing, one can use `pytest <https://docs.pytest.org/en/stable/>`_ or `nosetests <https://nose.readthedocs.io/en/latest/>`_. Both packages can be installed using pip. To run all the test script, run the either of the following::

    python -m pytest tests 
    nosetests -v
============
Main modules
============

Angular coordinates
-------------------
.. automodule:: t2c.angular_coordinates
    :members:

Bubble Statistics
-----------------
.. automodule:: t2c.bubble_stats
    :members:

Cosmological calculators
------------------------
.. automodule:: t2c.cosmology
    :members:

Foreground model
----------------
.. automodule:: t2c.foreground_model
    :members:

Identifying regions
-------------------
.. automodule:: t2c.identify_regions
    :members:

Lightcone
---------
.. automodule:: t2c.lightcone
    :members:

Telescope noise
---------------
.. automodule:: t2c.noise_model
    :members:

.. automodule:: t2c.telescope_functions
    :members:

Power spectrum
--------------
.. automodule:: t2c.power_spectrum
    :members:

Primary beam
------------
.. automodule:: t2c.primary_beam
    :members:

Reading simuation outputs
-------------------------
.. automodule:: t2c.read_files
    :members:

Smoothing
---------
.. automodule:: t2c.smoothing
    :members:

Point statistics
----------------
.. automodule:: t2c.statistics
    :members:

Temperature
-----------
.. automodule:: t2c.temperature
    :members:

Topology
--------
.. automodule:: t2c.topology
    :members:


.. include:: readme.rst

Contents
========
.. toctree::
   :maxdepth: 2

   installation
   tutorials
   contents
   contents_other
   contributing
   authors
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
=========
Changelog
=========

v2.1
----
* Modules to analyse 21 cm images added.
* Compatible with python 3 only

v1.1
----
* Sambit Giri and Garrelt Mellema made the package more general purpose in 2016 and renamed it to `tools21cm <https://tools21cm.readthedocs.io/>`_.
* Compatible with python 2 only.
* Post this version, we stop the development of the python 2 version.

v0.1
----
* Initial version of the package was called `c2raytools <https://ttt.astro.su.se/~gmell/c2raytools/build/>`_. 
* It was created by Hannes Jensen and Garrelt Mellema in 2014 to analyse data output from |c2ray|_ and |cubep3m|_ simulations.
* Compatible with python 2 only.

.. |c2ray| replace:: C\ :sup:`2`\RAY
.. _c2ray: https://github.com/garrelt/C2-Ray3Dm

.. |cubep3m| replace:: CUBEP\ :sup:`3`\M
.. _cubep3m: https://wiki.cita.utoronto.ca/index.php/CubePM
{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}

   {% block functions %}
   {% if functions %}
   .. rubric:: Functions

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: Classes

   .. autosummary::
      :toctree: {{ objname }}
      :template: class.rst
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: Exceptions

   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
