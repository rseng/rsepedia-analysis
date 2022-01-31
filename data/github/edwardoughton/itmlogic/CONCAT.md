# itmlogic – Longley-Rice Irregular Terrain Model

[![Build Status](https://travis-ci.org/edwardoughton/itmlogic.svg?branch=master)](https://travis-ci.org/edwardoughton/itmlogic)
[![Documentation Status](https://readthedocs.org/projects/itmlogic/badge/?version=latest)](https://itmlogic.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/edwardoughton/itmlogic/badge.svg?branch=master)](https://coveralls.io/github/edwardoughton/itmlogic?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02266/status.svg)](https://doi.org/10.21105/joss.02266)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3931350.svg)](https://doi.org/10.5281/zenodo.3931350)

**itmlogic** is a Python implementation of the classic Longley-Rice propagation model (v1.2.2)
and capable of estimating the signal propagation effects resulting from irregular terrain.

Software citation
-----------------

- Oughton, E.J., Russell, T., Johnson, J., Yardim, C., Kusuma, J., 2020. itmlogic:
  The Irregular Terrain Model by Longley and Rice. Journal of Open Source Software 5, 2266.
  https://doi.org/10.21105/joss.02266

Software purpose
----------------

This Python repo implements the model properties and algorithm defined in:

* Hufford, G. A., A. G. Longley, and W. A. Kissick (1982), A guide    to the use of the ITS
  Irregular Terrain Model in the area prediction mode, NTIA Report 82-100. (NTIS Order No.
  PB82-217977)
* Hufford, G. A. (1995) The ITS Irregular Terrain Model, version 1.2.2, the Algorithm.

**itmlogic** enables you to account for the radio propagation impacts occuring from irregular
terrain (hills, mountains etc.). For example, the image below shows the terrain undulation 
between the Crystal Palace (South London) transmitter and Mursley, Buckinghamshire, England.
Such estimates enable the engineering design of many types of wireless radio systems, including 
4G and 5G Radio Access Networks and wireless backhaul connections. 

Terrain profile slice: Crystal Palace (South London) to Mursley
---------------------------------------------------------------
![Example](/docs/_static/terrain_profile.png)


## Setup and configuration

All code for ``itmlogic`` is written in Python (Python>=3.7).

See requirements.txt for a full list of dependencies.


## Conda

The recommended installation method is to use conda, which handles packages and virtual
environments, along with the conda-forge channel which has a host of pre-built libraries
and packages.

Create a conda environment called ``itmlogic``:

    conda create --name itmlogic python=3.7 gdal

Activate it (run this each time you switch projects):

    conda activate itmlogic

First, install optional packages:

    conda install numpy fiona shapely rtree rasterio pyproj tqdm pytest rasterstats pandas matplotlib

Once in the new environment, to install ``itmlogic`` clone this repository and either run:

    python setup.py install

Or:

    python setup.py develop

You can first run the tests to make sure everything is working correctly:

    python -m pytest


Quick start
-----------

If you want to quickly generate results run using point-to-point mode run:

    python scripts/p2p.py

Or using area prediction mode run:

    python scripts/area.py

Results can then be visualized using:

    python vis/vis.py


Example results - Point-to-point mode
-------------------------------------
![Example](/docs/_static/p2p_results.png)


Example results - Area mode
---------------------------
![Example](/docs/_static/area_results.png)


Documentation
-------------

For more information, see the ``itmlogic`` [readthedocs documentation](https://itmlogic.readthedocs.io/en/latest/?badge=latest).


## Background

The model was developed by the Institute for Telecommunication Sciences (ITS) for frequencies
between 20 MHz and 20 GHz (named for Anita Longley & Phil Rice, 1968), and as a general
purpose model can be applied to a large variety of engineering problems. Based on
both electromagnetic theory and empirical statistical analyses of both terrain features and
radio measurements, the Longley-Rice Irregular Terrain Model predicts the median attenuation
of a radio signal as a function of distance and the variability of signal in time and in space.

The original NTIA disclaimer states:

> The ITM software was developed by NTIA. NTIA does not make any warranty of any kind, express,
implied or statutory, including, without limitation, the implied warranty of merchantability,
fitness for a particular purpose, non-infringement and data accuracy. NTIA does not warrant or
make any representations regarding the use of the software or the results thereof, including
but not limited to the correctness, accuracy, reliability or usefulness of the software or the
results. You can use, copy, modify, and redistribute the NTIA-developed software upon your
acceptance of these terms and conditions and upon your express agreement to provide appropriate
acknowledgments of NTIA's ownership of and development of the software by keeping this exact
text present in any copied or derivative works.


## Thanks for the support

The software repository **itmlogic** was written and developed at the [Environmental Change Institute, University of
Oxford](http://www.eci.ox.ac.uk) within the EPSRC-sponsored MISTRAL programme (EP/N017064/1),
as part of the [Infrastructure Transition Research Consortium](http://www.itrc.org.uk/)

## Contributors
- Edward J. Oughton (University of Oxford) (Software Engineering Lead)
- Tom Russell (University of Oxford) (Software Engineering)
- Joel Johnson (The Ohio State University) (ITM Modeling Lead)
- Caglar Yardim (The Ohio State University) (ITM Modeling)
- Julius Kusuma (Facebook Research) (ITM Modeling)

If you find an error or have a question, please submit an issue.

## Folder structure

The folder structure for the ``itmlogic`` package is summarized as follows, and matches the
box diagram highlighted in both the JOSS paper and the documentation:

    +---src
    |   +---itmlogic
    |   |   |   lrprop.py
    |   |   |   __init__.py
    |   |   |
    |   |   +---diffraction_attenuation
    |   |   |       adiff.py
    |   |   |       aknfe.py
    |   |   |       fht.py
    |   |   |
    |   |   +---los_attenuation
    |   |   |       alos.py
    |   |   |
    |   |   +---misc
    |   |   |       qerf.py
    |   |   |       qerfi.py
    |   |   |       qtile.py
    |   |   |
    |   |   +---preparatory_subroutines
    |   |   |       dlthx.py
    |   |   |       hzns.py
    |   |   |       qlra.py
    |   |   |       qlrpfl.py
    |   |   |       qlrps.py
    |   |   |       zlsq1.py
    |   |   |
    |   |   +---scatter_attenuation
    |   |   |       ahd.py
    |   |   |       ascat.py
    |   |   |       h0f.py
    |   |   |
    |   |   +---statistics
    |   |   |       avar.py
    |   |   |       curv.py
Contributing guidelines
=======================

We welcome contributions to ``itmlogic``.

When submitting a change to the repository, please first create an issue that covers the item
that you'd like to change, update or enhance. Once a discussion has yielded a vote of support
for that addition to the package, you are ready to submit a pull request.

If you are proposing a feature:

- Explain in detail how it would work.
- Keep the scope as narrow as possible, to make it easier to implement.

Create a branch for local development.
--------------------------------------

- Use the ``git checkout`` command to create your own branch, and pick a name that describes
the changes that you are making.

    $ git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.

Test the package
----------------

Ensure that the tests pass, and the documentation builds successfully::

    $ pytest
    $ make docs

Commit and push your changes
----------------------------

Once you are sure that all tests are passing, you can commit your changes and push to GitHub::

$ git add .
$ git commit -m "Your detailed description of your changes."
$ git push origin name-of-your-bugfix-or-feature
Submit a pull request on GitHub
When submitting a pull request:

All existing tests should pass. Please make sure that the test suite passes, both locally
and on Travis CI <https://travis-ci.org/nismod/itmlogic>_ Status on Travis will be visible on a
pull request. If you want to enable Travis CI on your own fork, please read the getting
started docs <https://docs.travis-ci.com/user/getting-started/>_.

New functionality should include tests. Please write reasonable tests for your code and make
sure that they pass on your pull request.

Classes, methods, functions, etc. should have docstrings. The first line of a docstring
should be a standalone summary. Parameters and return values should be documented explicitly.

The API documentation is automatically generated from docstrings, which should conform to
NumpPy styling. For examples, see the Napoleon docs <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>_.

Please note that tests are also run via Travis-CI on our documentation. So be sure that any
.rst file submissions are properly formatted and tests are passing.

Documentation Updates
=====================

Improving the documentation and testing for code already in itmlogic is a great way to get
started if you'd like to make a contribution. Please note that our documentation files are
in ReStructuredText (.rst) <http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>_
format and format your pull request accordingly.

To build the documentation, use the command::

    $ make docs

By default make docs will only rebuild the documentation if source files (e.g., .py or .rst
files) have changed. To force a rebuild, use make -B docs. You can preview the generated
documentation by opening docs/_build/html/index.html in a web browser.
---
title: 'itmlogic: The Irregular Terrain Model by Longley and Rice'
tags:
  - python
  - mobile telecommunications
  - propagation
  - longley-rice
authors:
  - name: Edward J Oughton
    orcid: 0000-0002-2766-008X
    affiliation:  "1, 2"
  - name: Tom Russell
    orcid: 0000-0002-0081-400X
    affiliation: 1
  - name: Joel Johnson
    affiliation: 3
  - name: Caglar Yardim
    affiliation: 3
  - name: Julius Kusuma
    affiliation: 4
affiliations:
  - name: Environmental Change Institute, University of Oxford
    index: 1
  - name: Computer Laboratory, University of Cambridge
    index: 2
  - name: ElectroScience Laboratory, The Ohio State University
    index: 3
  - name: Facebook Connectivity Lab, Facebook Research
    index: 4
date: 27 January 2020
bibliography: paper.bib
---

# Summary

Billions of people still do not have access to a reliable internet connection. One of the most effective ways to provide wide area access to a dispersed user base is via wireless radio technologies, such as cellular 4G or 5G [@Oughton:2018a]. The costs of wireless deployment are considerably lower than fixed alternatives, which is beneficial in areas with low per-capita income or adoption.

Data science methods can help us to more accurately identify unconnected groups and help to design least-cost internet access strategies. However, many of the statistical tools in the field are written in Python and therefore there is a language conflict with classic propagation models which have not yet been made available in this programming language.

The Longley-Rice Irregular Terrain Model is a classic propagation model developed by the Central Radio Propagation Laboratory during the 1960s in Colorado, USA, by A.G. Longley and P.L. Rice [@Longley:1968]. The model is still widely used throughout the cellular industry by Mobile Network Operators (MNOs) as it can predict long-term median transmission loss over irregular terrain. The original open-source model is available in Fortran or C++ [@ITS:2007].

This paper describes the ``itmlogic`` package, which provides a Python implementation of the Longley-Rice Irregular Terrain Model. It implements the classic model, enabling the quantification of propagation loss over irregular terrain. ``itmlogic`` is capable of predicting the the statistics of propagation loss given input parameters such as transmitter and receiver heights, frequency, surface permittivity, climate zone, and terrain information.

![Longley-Rice Irregular Terrain Model Scripts, Routines and Functions](lritm_box_diagram.png)

## Statement of Need

Smaller Mobile Network Operators may not have their own in-house engineering models. While other software packages are available, they need to be commercially licensed. Hence, this open-source package can help keep costs low for MNOs working to connect communities in rural and remote regions where costs are high and returns low.

## Uniqueness

The Longley-Rice model provides desirable features for a propagation model by predicting propagation loss given user parameters including information on the cumulative distribution function of the predicted loss so that the user can predict system performance as a function of reliability.

## Spatial Units

Inputs to the model are in the MKS system of units, so that the transmitter and receiver heights above the local terrain are specified in meters while ranges are specified in kilometers. Information on the terrain profile between transmitter and receiver is specified using a terrain profile with heights above sea level in meters. Information on expected atmospheric refractivity properties can also be input so that refractive effects are taken into account through the use of a modified Earth radius.

## The Longley-Rice Model

The model's purpose is to predict properties of the propagation loss in a communications link between a transmitter and receiver. The predicted propagation loss is described using cumulative distributions given the stochastic nature of radio wave propagation. The model originally was created in the 1960s when television broadcasting and terrestrial radio were important systems that required better engineering [@Hufford:1982]. The model is based on empirical curve fits to an extensive set of propagation measurements performed by the Institution for Telecommunication Sciences and other organizations.

Two modes of prediction are available: "area prediction" and "point-to-point". Area prediction mode uses a terrain irregularity parameter based on the inter-decile range of terrain elevations (the range after removing the top 10% and bottom 10% of elevations). Point-to-point mode uses a sample of up to 600 points from the terrain profile of the straight line between transmitter and receiver [@Hufford:1995].

## Applications

The median propagation loss estimates produced by ``itmlogic`` can be used with other link budget estimation models to assess the capacity, coverage and cost of 5G infrastructure [@Oughton:2019a]. For example, this could include application via the path loss module of the Python Simulator for Integrated Modelling of 5G, ``pysim5G`` [@Oughton:2019b]. The use of ``itmlogic`` is an improvement over previous analyses which have used propagation models which do not directly model the impacts of irregular
terrain for static [@Oughton:2018b] or moving users [@Oughton:2020].

## Acknowledgements

We would like to acknowledge two different funding sources which have enabled development of ``itmlogic``. Firstly, Oughton and Russell were funded by the UK EPSRC via the Infrastructure Transitions Research Consortium Mistral project (EP/N017064/1). Secondly, Johnson and Yardim were funded by Facebook Connectivity Lab's Rural Connectivity program.

# References
===============
Getting Started
===============

This tutorial provides explanation for the various modes available for ``itmlogic``.

Firstly a summary of the main primary parameters, secondary parameters and output values is
given.

Primary input parameters
------------------------

============= ============================
Parameters    Description
============= ============================
mdp           Controlling
dist          Distance
hg            Antenna structural heights
wn            Wave number
dh            Terrain irregularity
ens           Surface refractivity
gme           Earth's effective curvature
zgnd          Surface transfer impedance
ze            Effective antenna heights
dl            Horizon distances
the           Horizon elevation angles
============= ============================

Secondary parameters (computed in lrprop)
-----------------------------------------

=============== ============================
Parameters      Description
=============== ============================
dlsa            Line-of-sight distance
dx              Scatter distance
ael, ak1, ak2   Line-of-sight coefficients
aed, emd        Diffraction coefficients
aes, ems        Scatter coefficients
dls             Smooth earth horizon distances
dla             Total horizon distance
tha             Total bending angle
=============== ============================

Output values
-------------

================ ============================
Output values    Description
================ ============================
 kwx             Error indicator
 aref            Reference attenuation
================ ============================


Figure 1 provides an overview of the program flow, subroutines and statistics.


Longley-Rice Irregular Terrain Model Scripts, Routines and Functions
--------------------------------------------------------------------

.. image:: _static/lritm_box_diagram.png
    :target: _static/lritm_box_diagram.png


The model can run in one of two modes: 'area prediction mode' or 'point-to-point' prediction
mode.


Area Prediction Mode
--------------------

A reproducible example for the Crystal Palace radio transmitter (South London) is provided
using a single Digital Elevation Model (DEM) tile. Use the following to run the code:

.. code-block:: python

    python scripts/area.py

The repo already includes a DEM tile for London (see the .tif in the data folder).

For simplicity, this example specifies the coordinates of the ``transmitter`` as a point
feature. This is a standard `GeoJSON <https://geojson.org/>`_-like Python ``dict``, as you
would get from using `shapely
<https://shapely.readthedocs.io/en/latest/manual.html#geometric-objects>`_ to read point
features from a file:

.. code-block:: python

    transmitter = {
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': (-0.07491679518573545, 51.42413477117786)
            },
        'properties': {
            'id': 'Crystal Palace radio transmitter'
        }
    }

An estimated range (``cell_range``) is also provided as a maximum cell radius (in meters).

To assess landscape elevation the ``terrain_area`` function is imported from the
``terrain_module``. The function enables the estimation of the Terrain Irregularity Parameter
(``tip``), for a cell radius of 20,000 meters (20 km):

.. code-block:: python

    tip = terrain_area(dem_path, tx_coordinate_0, tx_coordinate_1, cell_range)

The ``tip`` is the inter-decile range for all elevation values (the range between the top
10% and bottom 10% of values). This parameter can then be passed to the ``itmlogic_area``
function:

.. code-block:: python

    output = itmlogic_area(tip)

As the ``itmlogic_area`` is used here to merely demonstrate the code functionality, a user will
need to adapt parameters to their specific scenario. For example, the user will want to
specify the specific antenna heights, frequency to be modelled and local atmospheric conditions.
The main user defined parameters can be set via the ``main_user_defined_parameters`` dict,
but environmental and statistical paramters will need to be adjusted by the user in the
``itmlogic_area`` function.

In the given scenario, the propagation loss across this terrain is estimated for a certain
distance, at a specific confidence level, and returned as a list of dicts named ``output``:

.. code-block:: python

    output = [
        {
            'distance_km': 10,
            'confidence_level_%': 50,
            'propagation_loss_dB': 111.6920084
        },
        {
            'distance_km': 10,
            'confidence_level_%': 90,
            'propagation_loss_dB': 121.5943795
        },
        ...
    ]

The results are then written to a csv file in the processed data folder ('uarea_output.csv).

We also provide an example which spans more than one coverage tile, as defined in:

.. code-block:: python

    python scripts/area_2tiles.py


Point-to-Point Mode
-------------------

In contrast to the area prediction mode, the point-to-point mode focuses on a single path
across an area of irregular terrain between a transmitter and receiver. To use the
reproducible example for p2p, run:

.. code-block:: bash

    python scripts/p2p.py

The example given is based on the original radio propagation scenario used which is between
the Crystal Palace radio transmitter in South London and a receiver in the small village of
Mursley in Buckinghamshire, England. For consistency, ``itmlogic`` also uses this example,
particularly for providing tests for the codebase, to guarantee reliability.

Like the area prediction function, the ``itmlogic_p2p`` is used here to merely demonstrate the
code functionality, so a user will need to adapt parameters to their specific scenario. For
example, the user will want to specify the specific antenna heights, frequency to be modelled
and local atmospheric conditions. The main user defined parameters can be set via the
``main_user_defined_parameters`` dict, but environmental and statistical paramters will need
to be adjusted by the user in the ``itmlogic_p2p`` function.

To begin, the transmitter is specified as a point feature:

.. code-block:: python

    transmitter = {
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': (-0.07491679518573545, 51.42413477117786)
            },
        'properties': {
            'id': 'Crystal Palace radio transmitter'
        }
    }

Along with the receiver:

.. code-block:: python

    receiver = {
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': (-0.8119433954872186, 51.94972494521946)
            },
        'properties': {
            'id': 'Mursley'
        }
    }

The terrain path is then specified as a line feature:

.. code-block:: python

    line = {
        'type': 'Feature',
        'geometry': {
            'type': 'LineString',
            'coordinates': [
                    (
                        transmitter['geometry']['coordinates'][0],
                        transmitter['geometry']['coordinates'][1]
                    ),
                    (
                        receiver['geometry']['coordinates'][0],
                        receiver['geometry']['coordinates'][1]
                    ),
                ]
            },
        'properties': {
            'id': 'terrain path'
        }
    }

Using the ``terrain_p2p`` function from the ``terrain_module`` we can get the terrain
profile, over a set distance, with each point across the terrain profile being returned as a
GeoJSON object.

.. code-block:: python

    measured_terrain_profile, distance_km, points = terrain_p2p(
        dem_folder, line
    )

A list of terrain elevation values (``measured_terrain_profile``) (in meters) is returned:

.. code-block:: python

    measured_terrain_profile = [
        109, 66, 28, 48, 29, 32, 29, 20, 13, 9...
    ]

These data can then be passed to the ``itmlogic_p2p`` function along with the distance (km)
of the link:

.. code-block:: python

    output = itmlogic_p2p(original_surface_profile_m, distance_km)

The results are returned in a list of dicts called ``output`` containing the path loss over
the link distance given certain reliability and confidence levels.

.. code-block:: python

    output = [
        {
            'distance_km': 77.8,
            'reliability_level_%': 1,
            'confidence_level_%': 50,
            'propagation_loss_dB': 128.5969039310673
        },
        {
            'distance_km': 77.8,
            'reliability_level_%': 1,
            'confidence_level_%': 90,
            'propagation_loss_dB': 137.64279211442656
        },
        ...
    ]

We also provide an example which spans more than one coverage tile, as defined in:

.. code-block:: python

    python scripts/p2p_2tiles.py
==========
Developers
==========

* Edward J. Oughton <e.oughton@ouce.ox.ac.uk>
* Tom Russell <tom.russell@ouce.ox.ac.uk>
* Joel Johnson <johnson.1374@osu.edu>
* Caglar Yardim <yardim.7@osu.edu>
* Julius Kusuma <jkusuma@fb.com>
.. itmlogic documentation master file, created by
   sphinx-quickstart on Mon Feb 17 18:32:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

itmlogic - Longley-Rice Irregular Terrain Model
=======================================================================

.. image:: https://readthedocs.org/projects/itmlogic/badge/?version=latest
    :target: https://itmlogic.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.com/edwardoughton/itmlogic.svg?branch=master
    :target: https://travis-ci.org/edwardoughton/itmlogic
    :alt: Build Status

.. image:: https://coveralls.io/repos/github/edwardoughton/itmlogic/badge.svg?branch=master
    :target: https://coveralls.io/github/edwardoughton/itmlogic?branch=master
    :alt: Coverage Status

.. image:: https://img.shields.io/badge/github-itmlogic-brightgreen
    :target: https://github.com/edwardoughton/itmlogic
    :alt: Source Code


Welcome to the documentation for itmlogic!

These docs provide an overview of the python implementation of the classic Longley-Rice
propagation model. For in-depth detail on the engineering model we advise you to refer to the
original Longley-Rice documents provided in the repo folder named `References`.

Description
===========

The Longley-Rice Irregular Terrain Model is one of the classic radio propagation models. The
model is still widely used, particularly in industry. In comparison with other radio models,
Longley-Rice accounts for various physical effects that result from irregular terrain.


Statement of Need
=================

Many engineering models for wireless networks completely ignore terrain effects. While
software packages are available to address this, they usually need to be commerically licensed.
Those that are open-source are usually in less commonly used languages such as Fortran.
This open-source python package overcomes these limitations by providing an easy-to-use Python
version. Given the popularity of Python, there is a need for an easy-to-use Python version of
this model.


Setup and configuration
=======================

All code for ``itmlogic`` is written in Python (Python>=3.7).

See requirements.txt for a full list of dependencies.


Conda
=====

The recommended installation method is to use conda, which handles packages and virtual
environments, along with the conda-forge channel which has a host of pre-built libraries
and packages.

Create a conda environment called ``itmlogic`` type::

    conda create --name itmlogic python=3.7 gdal

Activate it (run this each time you switch projects)::

    conda activate itmlogic

Install any optional packages::

    conda install numpy fiona shapely rtree rasterio pyproj tqdm pytest

For development purposes, clone this repository and run::

    python setup.py develop

Run the tests::

    python -m pytest


Quick start
===========

If you want to quickly generate results run using area prediction mode type::

    python scripts/area.py

Or using point-to-point mode type::

    python scripts/p2p.py


Contents
========

.. toctree::
   :maxdepth: 3

   Getting Started <getting-started>
   Reference <api/modules>

.. toctree::
   :maxdepth: 1

   License <license>
   Authors <authors>


Make Contact
============

- Report bugs, suggest features or view the source code `on GitHub`_.
    .. _on GitHub: https://github.com/edwardoughton/itmlogic
