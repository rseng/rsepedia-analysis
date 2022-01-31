---
title: 'HyRiver: Hydroclimate Data Retriever'
tags:
  - Python
  - hydrology
  - climate
  - web services
authors:
  - name: Taher Chegini
    orcid: 0000-0002-5430-6000
    affiliation: 1
  - name: Hong-Yi Li
    orcid: 0000-0002-9807-3851
    affiliation: 1
  - name: L. Ruby Leung
    orcid: 0000-0002-3221-9467
    affiliation: 2
affiliations:
 - name: University of Houston
   index: 1
 - name: Pacific Northwest National Laboratory
   index: 2
date: 26 February 2021
bibliography: paper.bib
---

# Summary

Over the last decade, increasing availability of web services for hydrology and
climatology data has facilitated publication of reproducible scientific research in hydrological
and climate studies. Such web services allow researchers to subset big databases and perform some
common data processing operations on the server-side. However, implementing such services increases
the technical complexity of code development as it requires sufficient understanding of their
underlying protocols to generate valid queries and filters. `HyRiver` bridges this gap
by providing a unified and simple Application Programming Interface (API) to web services that are
based on three of the most commonly used protocols for geo-spatial/temporal data publication:
REpresentational State Transfer (RESTful), Web Feature Services (WFS), and Web Map Services (WMS).
`HyRiver` is a software stack consisting of the following seven Python packages:

* [PyGeoHydro](https://github.com/cheginit/pygeohydro): Provides access to NWIS (National Water
  Information System), NID (National Inventory of Dams), HCDN-2009 (Hydro-Climatic Data Network),
  NLCD (National Land Cover Database), and SSEBop (operational Simplified Surface Energy Balance)
  databases. Moreover, it can generate an interactive map for exploring NWIS stations within a
  bounding box, compute categorical statistics of land use/land cover data, and plot five
  hydrologic signature graphs. There is also a helper function which returns a roughness
  coefficients lookup table for NLCD's land cover types. These coefficients can be
  useful for overland flow routing among other applications.
* [PyNHD](https://github.com/cheginit/pynhd): Provides the ability to navigate and subset
  National Hydrography Database [@Buto_2020], at medium- and high-resolution, using NLDI (Hydro
  Network-Linked Data Index), WaterData, and TNM (The National Map) web services. Additionally,
  it can retrieve over 30 catchment-scale attributes from
  [ScienceBase](https://www.sciencebase.gov/catalog/item/5669a79ee4b08895842a1d47)
  and [NHDPlus Value Added Attributes](https://www.hydroshare.org/resource/6092c8a62fac45be97a09bfd0b0bf726)
  that are linked to the NHDPlus database via Common Identifiers (ComIDs). `PyNHD` has some
  additional river network tools that use NHDPlus data for routing through a river network.
  This flow routing module is general and accepts any user-defined transport equation for
  computing flow accumulation through a given river network. It sorts the river network
  topologically from upstream to downstream, then accumulates a given attribute based on the
  user-defined transport equation.
* [Py3DEP](https://github.com/cheginit/py3dep): Gives access to topographic data through the
  3D Elevation Program (3DEP) service [@Thatcher_2020]. This package can pull 12 types of
  topographic data from the 3DEP service, such as Digital Elevation Model, slope, aspect, and
  hillshade.
* [PyDaymet](https://github.com/cheginit/pydaymet): Retrieves daily climate data as well as
  their monthly and annual summaries from the Daymet dataset [@Thornton_2020]. It is possible to
  request data for a single location as well as a grid (any valid geometrical shape) at 1-km
  spatial resolution.
* [PyGeoOGC](https://github.com/cheginit/pygeoogc): Generates valid queries for retrieving data
  from supported RESTful-, WMS-, and WFS-based services. Although these web services limit
  the number of features in a single query, under-the-hood, `PyGeoOGC` takes care of breaking down
  a large query into smaller queries according to specifications of the services. Additionally,
  this package offers several notable utilities, such as data re-projection and asynchronous data
  retrieval for speeding up sending/receiving queries.
* [PyGeoUtils](https://github.com/cheginit/pygeoutils): Converts responses from PyGeoOGC's
  supported web services to geo-dataframes (vector data type) or datasets (raster data type).
  Moreover, for gridded data, it can mask the output dataset based on any given geometry.
* [AsyncRetriever](https://github.com/cheginit/async_retriever): `AsyncRetriever` has only one
  purpose; asynchronously sending requests and retrieving responses as `text`, `binary`,
  or `json` objects. It uses persistent caching to speedup the retrieval even further.

These packages are standalone and users can install and work with them independently.
Furthermore, `PyGeoOGC`, `PyGeoUtils`, and `AsyncRetriever` are low-level engines of this software
stack that the other four packages utilize for providing access to some of the most popular
databases in the hydrology community. These two low-level packages are generic and developers can
use them for connecting and sending queries to any other web services that are based on the
protocols that `HyRiver` supports. Currently, `HyRiver` only supports hydrology and climatology
datasets within the US.

# Statement of need

Preparing input data for conducting studies, is often one of the most time-consuming steps. The
difficulties of processing such input data stem from the diverse data sources and types as well as
sizes. For example, hydrological modeling of watersheds might require climate data such as
precipitation and temperature, topology data such as Digital Elevation Model, and a river network.
Climate and topology data are available in raster format, and river network could be from a vector
data type. Additionally, these datasets often have large sizes and subsetting operations can be
computationally demanding. Geospatial web services can carry out subsetting and some common
geographic information system (GIS) operations on server-side. However, these services usually have
different specifications, thus implementing them can be technically challanging. Moreover,
since the underlying protocols of these web services are under active development by organizations
such as [Open Geospatial Consortium](https://www.ogc.org), keeping track of the latest developments
adds another level of complexity. `HyRiver` removes these barriers by providing
consistent and simple, yet configurable, interfaces to these web services. Since these
interfaces are web protocol-specific, not web service-specific, researchers can utilize `HyRiver`
to access a plethora of databases that are offered through RESTFul-, WFS-, and WMS-based services.

There are several open-source packages that offer similar functionalities. For example,
[hydrofunctions](https://github.com/mroberge/hydrofunctions) is a Python package that retrieves
streamflow data from NWIS and [ulmo](https://github.com/ulmo-dev/ulmo) is another Python package
that provides access to several public hydrology and climatology data.
[Sentinelhub-py](https://github.com/sentinel-hub/sentinelhub-py) can download and process
satellite images from six data sources through Python. `Dataretrieval` gives access to some of
the USGS (United States Geological Survey) databases and has two versions in
[R](https://github.com/USGS-R/dataRetrieval) and [Python](https://github.com/USGS-python/dataretrieval).
Another R Package called [HydroData](https://github.com/mikejohnson51/HydroData), provides access
to 15 earth system datasets. Although there are overlaps between `HyRiver` and these packages,
to the best of our knowledge, none of them offer access to the diverse data sources and types
that this software stack provides.

# Acknowledgements

T. Chegini was supported by the University of Houston’s internal funds and the US National Science
Foundation (EAR #1804560). H.-Y. Li and L. R. Leung were supported by the U.S. Department of
Energy Office of Science Biological and Environmental Research as part of the Earth System Model
Development and Regional and Global Model Analysis program areas through the collaborative,
multi-program Integrated Coastal Modeling (ICoM) project. PNNL is operated for the Department of
Energy by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References
.. highlight:: bash

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways to any of the packages that are included in HyRiver
project. The workflow is the same for all packages. In this page, a contribution workflow
for `PyGeoHydro <https://github.com/cheginit/pygeohydro>`__ is explained.

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/cheginit/pygeohydro/issues.

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

Other than new features that you might have in mind, you can look through
the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

PyGeoHydro could always use more documentation, whether as part of the
official PyGeoHydro docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/cheginit/pygeohydro/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up PyGeoHydro for local development.

1. Fork the PyGeoHydro repo through the GitHub website.
2. Clone your fork locally and add the main PyGeoHydro as the upstream remote:

.. code-block:: console

    $ git clone git@github.com:your_name_here/pygeohydro.git
    $ git remote add upstream git@github.com:cheginit/pygeohydro.git

3. Install your local copy into a virtualenv. Assuming you have Conda installed, this is how you
   can set up your fork for local development:

.. code-block:: console

    $ cd pygeohydro/
    $ conda env create -f ci/requirements/environment.yml
    $ conda activate pygeohydro-dev
    $ python -m pip install . --no-deps

4. Create a branch for local development:

.. code-block:: console

    $ git checkout -b bugfix-or-feature/name-of-your-bugfix-or-feature
    $ git push

5. Before you first commit, pre-commit hooks needs to be setup:

.. code-block:: console

    $ pre-commit install
    $ pre-commit run --all-files

6. Now you can make your changes locally, make sure to add a description of
   the changes to ``HISTORY.rst`` file and add extra tests, if applicable,
   to ``tests`` folder. Also, make sure to give yourself credit by adding
   your name at the end of the item(s) that you add in the history like this
   ``By `Taher Chegini <https://github.com/cheginit>`_``. Then,
   fetch the latest updates from the remote and resolve any merge conflicts:

.. code-block:: console

    $ git fetch upstream
    $ git merge upstream/name-of-your-branch

7. Then lint and test the code:

.. code-block:: console

    $ make lint

8. If you are making breaking changes make sure to reflect them in
   the documentation, ``README.rst``, and tests if necessary.

9. Commit your changes and push your branch to GitHub:

.. code-block:: console

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

10. Submit a pull request through the GitHub website.

Tips
----

To run a subset of tests:

.. code-block:: console

    $ pytest -k "test_name1 or test_name2"

Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run:

.. code-block:: console

    $ git tag -a vX.X.X -m "vX.X.X"
    $ git push --follow-tags

where ``X.X.X`` is the version number following the
`semantic versioning spec <https://semver.org>`__ i.e., MAJOR.MINOR.PATCH.
Then release the tag from Github and Github Actions will deploy it to PyPi.
Contributor Covenant Code of Conduct
====================================

Our Pledge
----------

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our
project and our community a harassment-free experience for everyone,
regardless of age, body size, disability, ethnicity, sex
characteristics, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance,
race, religion, or sexual identity and orientation.

Our Standards
-------------

Examples of behavior that contributes to creating a positive environment
include:

-  Using welcoming and inclusive language
-  Being respectful of differing viewpoints and experiences
-  Gracefully accepting constructive criticism
-  Focusing on what is best for the community
-  Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

-  The use of sexualized language or imagery and unwelcome sexual
   attention or advances
-  Trolling, insulting/derogatory comments, and personal or political
   attacks
-  Public or private harassment
-  Publishing others’ private information, such as a physical or
   electronic address, without explicit permission
-  Other conduct which could reasonably be considered inappropriate in a
   professional setting

Our Responsibilities
--------------------

Project maintainers are responsible for clarifying the standards of
acceptable behavior and are expected to take appropriate and fair
corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that they
deem inappropriate, threatening, offensive, or harmful.

Scope
-----

This Code of Conduct applies both within project spaces and in public
spaces when an individual is representing the project or its community.
Examples of representing a project or community include using an
official project e-mail address, posting via an official social media
account, or acting as an appointed representative at an online or
offline event. Representation of a project may be further defined and
clarified by project maintainers.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may
be reported by contacting the project team at tchegini@uh.edu. All
complaints will be reviewed and investigated and will result in a
response that is deemed necessary and appropriate to the circumstances.
The project team is obligated to maintain confidentiality with regard to
the reporter of an incident. Further details of specific enforcement
policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in
good faith may face temporary or permanent repercussions as determined
by other members of the project’s leadership.

Attribution
-----------

This Code of Conduct is adapted from the `Contributor
Covenant <https://www.contributor-covenant.org>`__, version 1.4,
available at
https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
=======
Credits
=======

Development Lead
----------------

* Taher Chegini <cheginit@gmail.com>

Contributors
------------

None yet. Why not be the first?
.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/hyriver_logo_text.png
    :target: https://github.com/cheginit/HyRiver-examples

|

.. |pygeohydro| image:: https://github.com/cheginit/pygeohydro/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pygeohydro/actions/workflows/test.yml
    :alt: Github Actions

.. |pygeoogc| image:: https://github.com/cheginit/pygeoogc/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pygeoogc/actions/workflows/test.yml
    :alt: Github Actions

.. |pygeoutils| image:: https://github.com/cheginit/pygeoutils/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pygeoutils/actions/workflows/test.yml
    :alt: Github Actions

.. |pynhd| image:: https://github.com/cheginit/pynhd/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pynhd/actions/workflows/test.yml
    :alt: Github Actions

.. |py3dep| image:: https://github.com/cheginit/py3dep/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/py3dep/actions/workflows/test.yml
    :alt: Github Actions

.. |pydaymet| image:: https://github.com/cheginit/pydaymet/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pydaymet/actions/workflows/test.yml
    :alt: Github Actions

.. |async| image:: https://github.com/cheginit/async_retriever/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/async_retriever/actions/workflows/test.yml
    :alt: Github Actions

.. |geoh_stat| image:: https://static.pepy.tech/personalized-badge/hydrodata?period=total&left_color=blue&right_color=yellowgreen&left_text=PyGeoHydro
    :target: https://github.com/cheginit/pygeohydro
    :alt: Download Stat

.. |ogc_stat| image:: https://static.pepy.tech/personalized-badge/pygeoogc?period=total&left_color=blue&right_color=yellowgreen&left_text=PyGeoOGC
    :target: https://github.com/cheginit/pygeoogc
    :alt: Download Stat

.. |utils_stat| image:: https://static.pepy.tech/personalized-badge/pygeoutils?period=total&left_color=blue&right_color=yellowgreen&left_text=PyGeoUtils
    :target: https://github.com/cheginit/pygeoutils
    :alt: Download Stat

.. |nhd_stat| image:: https://static.pepy.tech/personalized-badge/pynhd?period=total&left_color=blue&right_color=yellowgreen&left_text=PyNHD
    :target: https://github.com/cheginit/pynhd
    :alt: Download Stat

.. |3dep_stat| image:: https://static.pepy.tech/personalized-badge/py3dep?period=total&left_color=blue&right_color=yellowgreen&left_text=Py3DEP
    :target: https://github.com/cheginit/py3dep
    :alt: Download Stat

.. |day_stat| image:: https://static.pepy.tech/personalized-badge/pydaymet?period=total&left_color=blue&right_color=yellowgreen&left_text=PyDaymet
    :target: https://github.com/cheginit/pydaymet
    :alt: Download Stat

.. |async_stat| image:: https://static.pepy.tech/personalized-badge/async_retriever?period=total&left_color=blue&right_color=yellowgreen&left_text=AsyncRetriever
    :target: https://github.com/cheginit/async_retriever
    :alt: Download Stat

.. _PyGeoHydro: https://github.com/cheginit/pygeohydro
.. _PyGeoOGC: https://github.com/cheginit/pygeoogc
.. _PyGeoUtils: https://github.com/cheginit/pygeoutils
.. _PyNHD: https://github.com/cheginit/pynhd
.. _Py3DEP: https://github.com/cheginit/py3dep
.. _PyDaymet: https://github.com/cheginit/pydaymet

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

.. image:: https://readthedocs.org/projects/hyriver/badge/?version=latest
    :target: https://hyriver.readthedocs.io/en/latest/?badge=latest
    :alt: ReadTheDocs

.. image:: https://joss.theoj.org/papers/b0df2f6192f0a18b9e622a3edff52e77/status.svg
    :target: https://joss.theoj.org/papers/b0df2f6192f0a18b9e622a3edff52e77
    :alt: JOSS

=============== ==================================================================== ============
Package         Description                                                          CI
=============== ==================================================================== ============
|nhd_stat|      Navigate and subset NHDPlus (MR and HR) using web services           |pynhd|
|3dep_stat|     Access topographic data through National Map's 3DEP web service      |py3dep|
|geoh_stat|     Access NWIS, NID, WQP, HCDN 2009, NLCD, and SSEBop databases         |pygeohydro|
|day_stat|      Access Daymet for daily climate data both single pixel and gridded   |pydaymet|
|async_stat|    High-level API for asynchronous requests with persistent caching     |async|
|ogc_stat|      Send queries to any ArcGIS RESTful-, WMS-, and WFS-based services    |pygeoogc|
|utils_stat|    Convert responses from PyGeoOGC's supported web services to datasets |pygeoutils|
=============== ==================================================================== ============


HyRiver: Hydroclimate Data Retriever
=====================================

Features
--------

`HyRiver <https://hyriver.readthedocs.io>`__ is a software stack consisting of seven
Python libraries that are designed to aid in watershed analysis through web services.
Currently, this project only includes hydrology and climatology data
within the US. Some of the major capabilities of HyRiver are as follows:

* Easy access to many web services for subsetting data on server-side and returning the requests
  as masked Datasets or GeoDataFrames.
* Splitting large requests into smaller chunks, under-the-hood, since web services often limit
  the number of features per request. So the only bottleneck for subsetting the data
  is your local machine memory.
* Navigating and subsetting NHDPlus database (both medium- and high-resolution) using web services.
* Cleaning up the vector NHDPlus data, fixing some common issues, and computing vector-based
  accumulation through a river network.
* A URL inventory for some of the popular (and tested) web services.
* Some utilities for manipulating the obtained data and their visualization.

Please visit `examples <https://hyriver.readthedocs.io/en/latest/examples.html>`__
webpage to see some example notebooks. You can also try this project without installing
it on you system by clicking on the binder badge. A Jupyter Lab
instance with the HyRiver software stack pre-installed will be launched in your web browser
and you can start coding!

Please note that this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

Moreover, requests for additional databases and functionalities can be submitted via issue trackers
of packages.

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/flow_accumulation.png
    :target: https://github.com/cheginit/HyRiver-examples

Installation
------------

You can install all the packages using ``pip``:

.. code-block:: console

    $ pip install py3dep pynhd pygeohydro pydaymet pygeoogc pygeoutils async_retriever

or ``conda``:

.. code-block:: console

    $ conda install -c conda-forge py3dep pynhd pygeohydro pydaymet pygeoogc pygeoutils async_retriever

or ``mamba`` (recommended):

.. code-block:: console

    $ mamba install -c conda-forge --strict-channel-priority py3dep pynhd pygeohydro pydaymet pygeoogc pygeoutils async_retriever
Software Stack
==============

A detailed description of each component of the HyRiver software stack.

.. panels::
    :container: container-lg pb-3
    :body: text-center
    :column: col-lg-3 col-md-3 col-sm-4 col-xs-6 p-2

    .. link-button:: readme/pynhd
        :type: ref
        :text: PyNHD
        :classes: stretched-link

    ---

    .. link-button:: readme/pygeohydro
        :type: ref
        :text: PyGeoHydro
        :classes: stretched-link

    ---

    .. link-button:: readme/py3dep
        :type: ref
        :text: Py3DEP
        :classes: stretched-link

    ---

    .. link-button:: readme/pydaymet
        :type: ref
        :text: PyDaymet
        :classes: stretched-link

    ---

    .. link-button:: readme/async_retriever
        :type: ref
        :text: AsyncRetriever
        :classes: stretched-link

    ---

    .. link-button:: readme/pygeoogc
        :type: ref
        :text: PyGeoOGC
        :classes: stretched-link

    ---

    .. link-button:: readme/pygeoutils
        :type: ref
        :text: PyGeoUtils
        :classes: stretched-link

.. toctree::
    :maxdepth: 1
    :hidden:

    PyNHD <readme/pynhd>
    PyGeoHydro <readme/pygeohydro>
    Py3DEP <readme/py3dep>
    PyDaymet <readme/pydaymet>
    AsyncRetriever<readme/async_retriever>
    PyGeoOGC <readme/pygeoogc>
    PyGeoUtils <readme/pygeoutils>
.. include:: ../../CONTRIBUTING.rst
Example Gallery
===============

The following notebooks demonstrate capabilities of the HyRiver software stack.

.. nbgallery::
    :glob:

    notebooks/*
Changelogs
==========

.. panels::
    :container: container-lg pb-3
    :body: text-center
    :column: col-lg-3 col-md-3 col-sm-4 col-xs-6 p-2

    .. link-button:: changelogs/pynhd
        :type: ref
        :text: PyNHD
        :classes: stretched-link

    ---

    .. link-button:: changelogs/pygeohydro
        :type: ref
        :text: PyGeoHydro
        :classes: stretched-link

    ---

    .. link-button:: changelogs/py3dep
        :type: ref
        :text: Py3DEP
        :classes: stretched-link

    ---

    .. link-button:: changelogs/pydaymet
        :type: ref
        :text: PyDaymet
        :classes: stretched-link

    ---

    .. link-button:: changelogs/async_retriever
        :type: ref
        :text: AsyncRetriever
        :classes: stretched-link

    ---

    .. link-button:: changelogs/pygeoogc
        :type: ref
        :text: PyGeoOGC
        :classes: stretched-link

    ---

    .. link-button:: changelogs/pygeoutils
        :type: ref
        :text: PyGeoUtils
        :classes: stretched-link

.. toctree::
    :maxdepth: 1
    :hidden:

    PyNHD <changelogs/pynhd>
    PyGeoHydro <changelogs/pygeohydro>
    Py3DEP <changelogs/py3dep>
    PyDaymet <changelogs/pydaymet>
    AsyncRetriever<changelogs/async_retriever>
    PyGeoOGC <changelogs/pygeoogc>
    PyGeoUtils <changelogs/pygeoutils>
.. include:: ../../AUTHORS.rst
.. highlight:: bash

===============
Getting Started
===============

Why HyRiver?
------------

Some of the major capabilities of HyRiver are as follows:

* Easy access to many web services for subsetting data on server-side and returning the requests
  as masked Datasets or GeoDataFrames.
* Splitting large requests into smaller chunks, under-the-hood, since web services often limit
  the number of features per request. So the only bottleneck for subsetting the data
  is your local machine memory.
* Navigating and subsetting NHDPlus database (both medium- and high-resolution) using web services.
* Cleaning up the vector NHDPlus data, fixing some common issues, and computing vector-based
  accumulation through a river network.
* A URL inventory for some of the popular (and tested) web services.
* Some utilities for manipulating the obtained data and their visualization.

Installation
------------

You can install all the packages using ``pip``:

.. code-block:: console

    $ pip install py3dep pynhd pygeohydro pydaymet pygeoogc pygeoutils async_retriever

Please note that installation with ``pip`` fails if ``libgdal`` is not installed on your system.
You should install this package manually beforehand. For example, on Ubuntu-based distros
the required package is ``libgdal-dev``. If this package is installed on your system
you should be able to run ``gdal-config --version`` successfully.

Alternatively, you can use ``conda`` or ``mamba`` (recommended) to install these packages from
the ``conda-forge`` repository:

.. code-block:: console

    $ conda install -c conda-forge py3dep pynhd pygeohydro pydaymet pygeoogc pygeoutils async_retriever

or:

.. code-block:: console

    $ mamba install -c conda-forge --strict-channel-priority py3dep pynhd pygeohydro pydaymet pygeoogc pygeoutils async_retriever

Dependencies
------------

.. panels::
    :column: col-lg-12 col-md-12 col-sm-12 col-xs-12 p-2

    .. tabbed:: PyNHD

        - ``async_retriever``
        - ``cytoolz``
        - ``geopandas``
        - ``networkx``
        - ``numpy``
        - ``pandas``
        - ``pyarrow``
        - ``pygeoogc``
        - ``pygeoutils``
        - ``requests``
        - ``shapely``
        - ``simplejson``

    .. tabbed:: PyGeoHydro

        - ``async_retriever``
        - ``defusedxml``
        - ``folium``
        - ``geopandas``
        - ``lxml``
        - ``matplotlib``
        - ``numpy``
        - ``openpyxl``
        - ``pandas``
        - ``pygeoogc``
        - ``pygeoutils``
        - ``pynhd``
        - ``rasterio``
        - ``shapely``

    .. tabbed:: Py3DEP

        - ``async_retriever``
        - ``click``
        - ``cytoolz``
        - ``numpy``
        - ``pydantic``
        - ``pygeoogc``
        - ``pygeoutils``
        - ``rasterio``
        - ``scipy``
        - ``shapely``
        - ``xarray``

    .. tabbed:: PyDaymet

        - ``async_retriever``
        - ``click``
        - ``dask``
        - ``lxml``
        - ``numpy``
        - ``pandas``
        - ``py3dep``
        - ``pygeoogc``
        - ``pygeoutils``
        - ``rasterio``
        - ``scipy``
        - ``shapely``
        - ``xarray``

.. panels::
    :column: col-lg-12 col-md-12 col-sm-12 col-xs-12 p-2

    .. tabbed:: PyGeoOGC

        - ``async_retriever``
        - ``cytoolz``
        - ``defusedxml``
        - ``owslib``
        - ``pydantic``
        - ``pyproj``
        - ``pyyaml``
        - ``requests``
        - ``shapely``
        - ``simplejson``
        - ``urllib3``

    .. tabbed:: PyGeoUtils

        - ``affine``
        - ``dask``
        - ``geopandas``
        - ``netcdf4``
        - ``numpy``
        - ``orjson``
        - ``pygeoogc``
        - ``pyproj``
        - ``rasterio``
        - ``shapely``
        - ``xarray``

    .. tabbed:: AsyncRetriever

        - ``aiohttp-client-cache``
        - ``aiohttp[speedups]``
        - ``aiosqlite``
        - ``cytoolz``
        - ``nest-asyncio``
        - ``orjson``

Additionally, you can also install ``bottleneck``, ``pygeos``, and ``rtree`` to improve
performance of ``xarray`` and ``geopandas``. For handling vector and
raster data projections, ``cartopy`` and ``rioxarray`` are useful.
.. |async| image:: https://img.shields.io/pypi/v/async_retriever?label=AsyncRetriever&color=green
    :target: https://github.com/cheginit/async_retriever
    :alt: PyPi Version

.. |pygeohydro| image:: https://img.shields.io/pypi/v/pygeohydro?label=PyGeoHydro&color=green
    :target: https://github.com/cheginit/pygeohydro
    :alt: PyPi Version

.. |pygeoogc| image:: https://img.shields.io/pypi/v/pygeoogc?label=PyGeoOGC&color=green
    :target: https://github.com/cheginit/pygeoogc
    :alt: PyPi Version

.. |pygeoutils| image:: https://img.shields.io/pypi/v/pygeoutils?label=PyGeoUtils&color=green
    :target: https://github.com/cheginit/pygeoutils
    :alt: PyPi Version

.. |pynhd| image:: https://img.shields.io/pypi/v/pynhd?label=PyNHD&color=green
    :target: https://github.com/cheginit/pynhd
    :alt: PyPi Version

.. |py3dep| image:: https://img.shields.io/pypi/v/py3dep?label=Py3DEP&color=green
    :target: https://github.com/cheginit/py3dep
    :alt: PyPi Version

.. |pydaymet| image:: https://img.shields.io/pypi/v/pydaymet?label=PyDaymet&color=green
    :target: https://github.com/cheginit/pydaymet
    :alt: PyPi Version

.. |async_stat| image:: https://pepy.tech/badge/async_retriever
    :target: https://pepy.tech/project/async_retriever

.. |pygeohydro_stat| image:: https://pepy.tech/badge/hydrodata
    :target: https://pepy.tech/project/hydrodata
    :alt: Download Stat

.. |pygeoogc_stat| image:: https://pepy.tech/badge/pygeoogc
    :target: https://pepy.tech/project/pygeoogc
    :alt: Download Stat

.. |pygeoutils_stat| image:: https://pepy.tech/badge/pygeoutils
    :target: https://pepy.tech/project/pygeoutils
    :alt: Download Stat

.. |pynhd_stat| image:: https://pepy.tech/badge/pynhd
    :target: https://pepy.tech/project/pynhd
    :alt: Download Stat

.. |py3dep_stat| image:: https://pepy.tech/badge/py3dep
    :target: https://pepy.tech/project/py3dep
    :alt: Download Stat

.. |pydaymet_stat| image:: https://pepy.tech/badge/pydaymet
    :target: https://pepy.tech/project/pydaymet
    :alt: Download Stat

Hydroclimate Data Retriever
===========================

.. image:: https://img.shields.io/pypi/pyversions/pygeoogc.svg
    :target: https://pypi.python.org/pypi/pygeoogc
    :alt: Python Versions

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

.. image:: https://readthedocs.org/projects/hyriver/badge/?version=latest
    :target: https://hyriver.readthedocs.io/en/latest/?badge=latest
    :alt: ReadTheDocs

.. image:: https://joss.theoj.org/papers/b0df2f6192f0a18b9e622a3edff52e77/status.svg
    :target: https://joss.theoj.org/papers/b0df2f6192f0a18b9e622a3edff52e77
    :alt: JOSS

|

`HyRiver <https://hyriver.readthedocs.io>`__ (formerly named
`hydrodata <https://pypi.org/project/hydrodata>`__) is a suite of Python packages
that provides a unified API for retrieving geospatial/temporal data from various
web services. HyRiver includes two categories of packages:

- Low-level APIs for accessing any of the supported web services, i.e., ArcGIS RESTful,
  WMS, and WFS.
- High-level APIs for accessing some of the most commonly used datasets in hyrdology and
  climatology studies. Currently, this project only includes hydrology and climatology data
  within the US.

You can watch these videos for a quick overview of ``HyRiver``:

* `Pangeo Showcase <https://discourse.pangeo.io/t/may-26-2021-accessing-hydrology-and-climatology-database-using-web-services-through-python/1521>`__
* `ESIP IT&I <https://youtu.be/Wz8Y5G9oy-M?t=1838>`__

.. toctree::
    :maxdepth: 1
    :hidden:

    getting_started
    readmes
    examples

.. toctree::
    :maxdepth: 1
    :hidden:

    autoapi/index
    history
    contributing
    authors
    license

.. panels::
    :column: col-lg-12 p-2

    **High-level APIs for accessing some pre-configured web services**
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. tabbed:: PyNHD

        Navigate and subset mid- and high-res NHD, NHDPlus, and NHDPlus VAA
        using WaterData, NLDI, ScienceBase, and The National Map web services.

        |pynhd| |pynhd_stat|

    .. tabbed:: PyGeoHydro

        Access NWIS, NID, HCDN 2009, NLCD, and SSEBop databases.

        |pygeohydro| |pygeohydro_stat|

    .. tabbed:: Py3DEP

        Access topographic data through The National Map's 3DEP web service.

        |py3dep| |py3dep_stat|

    .. tabbed:: PyDaymet

        Access Daymet for daily, monthly and annual summaries of climate data
        at 1-km scale for both single pixels and gridded.

        |pydaymet| |pydaymet_stat|

.. panels::
    :column: col-lg-12 p-2

    **Low-level APIs for connecting to supported web service protocols**
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. tabbed:: PyGeoOGC

        Send queries to and receive responses from any ArcGIS RESTful-, WMS-,
        and WFS-based services.

        |pygeoogc| |pygeoogc_stat|

    .. tabbed:: PyGeoUtils

        Convert responses from PyGeoOGC's supported web services protocols into
        geospatial and raster datasets.

        |pygeoutils| |pygeoutils_stat|

    .. tabbed:: AsyncRetriever

        Asynchronous send/receive requests with persistent caching.

        |async| |async_stat|
API References
==============

.. panels::
    :container: container-lg pb-3
    :body: text-center
    :column: col-lg-3 col-md-3 col-sm-4 col-xs-6 p-2

    .. link-button:: /autoapi/pynhd/index
        :type: ref
        :text: PyNHD
        :classes: stretched-link

    |pynhd|

    ---

    .. link-button:: /autoapi/pygeohydro/index
        :type: ref
        :text: PyGeoHydro
        :classes: stretched-link

    |pygeohydro|

    ---

    .. link-button:: /autoapi/py3dep/index
        :type: ref
        :text: Py3DEP
        :classes: stretched-link

    |py3dep|

    ---

    .. link-button:: /autoapi/pydaymet/index
        :type: ref
        :text: PyDaymet
        :classes: stretched-link

    |pydaymet|

    ---

    .. link-button:: /autoapi/pygeoogc/index
        :type: ref
        :text: PyGeoOGC
        :classes: stretched-link

    |pygeoogc|

    ---

    .. link-button:: /autoapi/pygeoutils/index
        :type: ref
        :text: PyGeoUtils
        :classes: stretched-link

    |pygeoutils|

    ---

    .. link-button:: /autoapi/async_retriever/index
        :type: ref
        :text: AsyncRetriever
        :classes: stretched-link

    |async|

.. toctree::
   :titlesonly:
   :hidden:

   PyNHD </autoapi/pynhd/index>
   PyGeoHydro </autoapi/pygeohydro/index>
   Py3DEP </autoapi/py3dep/index>
   PyDaymet </autoapi/pydaymet/index>
   AsyncRetriever </autoapi/async_retriever/index>
   PyGeoOGC </autoapi/pygeoogc/index>
   PyGeoUtils </autoapi/pygeoutils/index>

.. |pygeohydro| image:: https://img.shields.io/pypi/v/pygeohydro?label=%20&color=blue
    :target: https://github.com/cheginit/pygeohydro
    :alt: PyPi Version

.. |pygeoogc| image:: https://img.shields.io/pypi/v/pygeoogc?label=%20&color=blue
    :target: https://github.com/cheginit/pygeoogc
    :alt: PyPi Version

.. |pygeoutils| image:: https://img.shields.io/pypi/v/pygeoutils?label=%20&color=blue
    :target: https://github.com/cheginit/pygeoutils
    :alt: PyPi Version

.. |pynhd| image:: https://img.shields.io/pypi/v/pynhd?label=%20&color=blue
    :target: https://github.com/cheginit/pynhd
    :alt: PyPi Version

.. |py3dep| image:: https://img.shields.io/pypi/v/py3dep?label=%20&color=blue
    :target: https://github.com/cheginit/py3dep
    :alt: PyPi Version

.. |pydaymet| image:: https://img.shields.io/pypi/v/pydaymet?label=%20&color=blue
    :target: https://github.com/cheginit/pydaymet
    :alt: PyPi Version

.. |async| image:: https://img.shields.io/pypi/v/async_retriever?label=%20&color=blue
    :target: https://github.com/cheginit/async_retriever
    :alt: PyPi Version
:py:mod:`pygeoutils`
====================

.. py:module:: pygeoutils

.. autoapi-nested-parse::

   Top-level package for PyGeoUtils.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   pygeoutils/index.rst
   utils/index.rst


Package Contents
----------------












:py:mod:`pygeoutils.pygeoutils`
===============================

.. py:module:: pygeoutils.pygeoutils

.. autoapi-nested-parse::

   Some utilities for manipulating GeoSpatial data.



Module Contents
---------------

.. py:class:: Coordinates

   Generate validated and normalized coordinates in WGS84.

   :Parameters: * **lon** (:class:`float` or :class:`list` of :class:`floats`) -- Longitude(s) in decimal degrees.
                * **lat** (:class:`float` or :class:`list` of :class:`floats`) -- Latitude(s) in decimal degrees.

   .. rubric:: Examples

   >>> from pygeoutils import Coordinates
   >>> c = Coordinates([460, 20, -30], [80, 200, 10])
   >>> c.points.x.tolist()
   [100.0, -30.0]

   .. py:method:: points(self)
      :property:

      Get validate coordinate as a ``geopandas.GeoSeries``.



.. py:class:: GeoBSpline(points, npts_sp, degree = 3)

   Create B-spline from a geo-dataframe of points.

   :Parameters: * **points** (:class:`geopandas.GeoDataFrame` or :class:`geopandas.GeoSeries`) -- Input points as a ``GeoDataFrame`` or ``GeoSeries`` in a projected CRS.
                * **npts_sp** (:class:`int`) -- Number of points in the output spline curve.
                * **degree** (:class:`int`, *optional*) -- Degree of the spline. Should be less than the number of points and
                  greater than 1. Default is 3.

   .. rubric:: Examples

   >>> from pygeoutils import GeoBSpline
   >>> import geopandas as gpd
   >>> xl, yl = zip(
   ...     *[
   ...         (-97.06138, 32.837),
   ...         (-97.06133, 32.836),
   ...         (-97.06124, 32.834),
   ...         (-97.06127, 32.832),
   ...     ]
   ... )
   >>> pts = gpd.GeoSeries(gpd.points_from_xy(xl, yl, crs="epsg:4326"))
   >>> sp = GeoBSpline(pts.to_crs("epsg:3857"), 5).spline
   >>> pts_sp = gpd.GeoSeries(gpd.points_from_xy(sp.x, sp.y, crs="epsg:3857"))
   >>> pts_sp = pts_sp.to_crs("epsg:4326")
   >>> list(zip(pts_sp.x, pts_sp.y))
   [(-97.06138, 32.837),
   (-97.06135, 32.83629),
   (-97.06131, 32.83538),
   (-97.06128, 32.83434),
   (-97.06127, 32.83319)]

   .. py:method:: spline(self)
      :property:

      Get the spline as a ``Spline`` object.



.. py:function:: arcgis2geojson(arcgis, id_attr = None)

   Convert ESRIGeoJSON format to GeoJSON.

   .. rubric:: Notes

   Based on `arcgis2geojson <https://github.com/chris48s/arcgis2geojson>`__.

   :Parameters: * **arcgis** (:class:`str` or :class:`binary`) -- The ESRIGeoJSON format str (or binary)
                * **id_attr** (:class:`str`, *optional*) -- ID of the attribute of interest, defaults to ``None``.

   :returns: :class:`dict` -- A GeoJSON file readable by GeoPandas.


.. py:function:: geo2polygon(geometry, geo_crs, crs)

   Convert a geometry to a Shapely's Polygon and transform to any CRS.

   :Parameters: * **geometry** (:class:`Polygon` or :class:`tuple` of :class:`length 4`) -- Polygon or bounding box (west, south, east, north).
                * **geo_crs** (:class:`str` or :class:`pyproj.CRS`) -- Spatial reference of the input geometry
                * **crs** (:class:`str` or :class:`pyproj.CRS`) -- Target spatial reference.

   :returns: :class:`Polygon` -- A Polygon in the target CRS.


.. py:function:: get_transform(ds, ds_dims = ('y', 'x'))

   Get transform of a ``xarray.Dataset`` or ``xarray.DataArray``.

   :Parameters: * **ds** (:class:`xarray.Dataset` or :class:`xarray.DataArray`) -- The dataset(array) to be masked
                * **ds_dims** (:class:`tuple`, *optional*) -- Names of the coordinames in the dataset, defaults to ``("y", "x")``.
                  The order of the dimension names must be (vertical, horizontal).

   :returns: :class:`rasterio.Affine`, :class:`int`, :class:`int` -- The affine transform, width, and height


.. py:function:: gtiff2xarray(r_dict, geometry = None, geo_crs = None, ds_dims = None, driver = None, all_touched = False, nodata = None, drop = True)

   Convert (Geo)Tiff byte responses to ``xarray.Dataset``.

   :Parameters: * **r_dict** (:class:`dict`) -- Dictionary of (Geo)Tiff byte responses where keys are some names that are used
                  for naming each responses, and values are bytes.
                * **geometry** (:class:`Polygon`, :class:`MultiPolygon`, or :class:`tuple`, *optional*) -- The geometry to mask the data that should be in the same CRS as the r_dict.
                  Defaults to ``None``.
                * **geo_crs** (:class:`str` or :class:`pyproj.CRS`, *optional*) -- The spatial reference of the input geometry, defaults to ``None``. This
                  argument should be given when ``geometry`` is given.
                * **ds_dims** (:class:`tuple` of :class:`str`, *optional*) -- The names of the vertical and horizontal dimensions (in that order)
                  of the target dataset, default to None. If None, dimension names are determined
                  from a list of common names.
                * **driver** (:class:`str`, *optional*) -- A GDAL driver for reading the content, defaults to automatic detection. A list of
                  the drivers can be found here: https://gdal.org/drivers/raster/index.html
                * **all_touched** (:class:`bool`, *optional*) -- Include a pixel in the mask if it touches any of the shapes.
                  If False (default), include a pixel only if its center is within one
                  of the shapes, or if it is selected by Bresenham's line algorithm.
                * **nodata** (:class:`float` or :class:`int`, *optional*) -- The nodata value of the raster, defaults to None, i.e., is determined from the raster.
                * **drop** (:class:`bool`, *optional*) -- If True, drop the data outside of the extent of the mask geometries.
                  Otherwise, it will return the same raster with the data masked.
                  Default is True.

   :returns: :class:`xarray.Dataset` or :class:`xarray.DataAraay` -- Parallel (with dask) dataset or dataarray.


.. py:function:: json2geodf(content, in_crs = DEF_CRS, crs = DEF_CRS)

   Create GeoDataFrame from (Geo)JSON.

   :Parameters: * **content** (:class:`dict` or :class:`list` of :class:`dict`) -- A (Geo)JSON dictionary e.g., response.json() or a list of them.
                * **in_crs** (:class:`str` or :class:`pyproj.CRS`) -- CRS of the content, defaults to ``epsg:4326``.
                * **crs** (:class:`str` or :class:`pyproj.CRS`, *optional*) -- The target CRS of the output GeoDataFrame, defaults to ``epsg:4326``.

   :returns: :class:`geopandas.GeoDataFrame` -- Generated geo-data frame from a GeoJSON


.. py:function:: xarray2geodf(da, dtype, mask_da = None, connectivity = 8)

   Vectorize a ``xarray.DataArray`` to a ``geopandas.GeoDataFrame``.

   :Parameters: * **da** (:class:`xarray.DataArray`) -- The dataarray to vectorize.
                * **dtype** (:class:`type`) -- The data type of the dataarray. Valid types are ``int16``, ``int32``,
                  ``uint8``, ``uint16``, and ``float32``.
                * **mask_da** (:class:`xarray.DataArray`, *optional*) -- The dataarray to use as a mask, defaults to ``None``.
                * **connectivity** (:class:`int`, *optional*) -- Use 4 or 8 pixel connectivity for grouping pixels into features,
                  defaults to 8.

   :returns: :class:`geopandas.GeoDataFrame` -- The vectorized dataarray.


.. py:function:: xarray_geomask(ds, geometry, crs, all_touched = False, drop = True, from_disk = False)

   Mask a ``xarray.Dataset`` based on a geometry.

   :Parameters: * **ds** (:class:`xarray.Dataset` or :class:`xarray.DataArray`) -- The dataset(array) to be masked
                * **geometry** (:class:`Polygon`, :class:`MultiPolygon`) -- The geometry to mask the data
                * **crs** (:class:`str` or :class:`pyproj.CRS`) -- The spatial reference of the input geometry
                * **all_touched** (:class:`bool`, *optional*) -- Include a pixel in the mask if it touches any of the shapes.
                  If False (default), include a pixel only if its center is within one
                  of the shapes, or if it is selected by Bresenham's line algorithm.
                * **drop** (:class:`bool`, *optional*) -- If True, drop the data outside of the extent of the mask geometries.
                  Otherwise, it will return the same raster with the data masked.
                  Default is True.
                * **from_disk** (:class:`bool`, *optional*) -- If True, it will clip from disk using rasterio.mask.mask if possible.
                  This is beneficial when the size of the data is larger than memory.
                  Default is False.

   :returns: :class:`xarray.Dataset` or :class:`xarray.DataArray` -- The input dataset with a mask applied (np.nan)


:py:mod:`pygeoutils.utils`
==========================

.. py:module:: pygeoutils.utils

.. autoapi-nested-parse::

   Some utilities for manipulating GeoSpatial data.



Module Contents
---------------

.. py:class:: Attrs



   Attributes of a GTiff byte response.


.. py:class:: Convert

   Functions to Convert an ArcGIS JSON object to a GeoJSON object.

   .. py:method:: attributes(self, arcgis, geojson)

      Get the attributes from an ArcGIS JSON object.


   .. py:method:: coords(self, arcgis, geojson)

      Get the bounds from an ArcGIS JSON object.


   .. py:method:: features(self, arcgis, geojson)

      Get the features from an ArcGIS JSON object.


   .. py:method:: geometry(arcgis, geojson)
      :staticmethod:

      Get the geometry from an ArcGIS JSON object.


   .. py:method:: get_outer_rings(rings)
      :staticmethod:

      Get outer rings and holes in a list of rings.


   .. py:method:: get_uncontained_holes(outer_rings, holes)
      :staticmethod:

      Get all the uncontstrained holes.


   .. py:method:: isnumber(nums)
      :staticmethod:

      Check if all items in a list are numbers.


   .. py:method:: paths(arcgis, geojson)
      :staticmethod:

      Get the paths from an ArcGIS JSON object.


   .. py:method:: points(arcgis, geojson)
      :staticmethod:

      Get the points from an ArcGIS JSON object.


   .. py:method:: rings(self, arcgis, _)

      Get the rings from an ArcGIS JSON object.


   .. py:method:: xy(self, arcgis, geojson)

      Get the xy coordinates from an ArcGIS JSON object.



.. py:function:: convert(arcgis, id_attr = None)

   Convert an ArcGIS JSON object to a GeoJSON object.


.. py:function:: get_bounds(ds, ds_dims = ('y', 'x'))

   Get bounds of a ``xarray.Dataset`` or ``xarray.DataArray``.

   :Parameters: * **ds** (:class:`xarray.Dataset` or :class:`xarray.DataArray`) -- The dataset(array) to be masked
                * **ds_dims** (:class:`tuple`, *optional*) -- Names of the coordinames in the dataset, defaults to ``("y", "x")``.
                  The order of the dimension names must be (vertical, horizontal).

   :returns: :class:`tuple` -- The bounds in the order of (left, bottom, right, top)


.. py:function:: get_dim_names(ds)

   Get vertical and horizontal dimension names.


.. py:function:: get_nodata(src)

   Get the nodata value of a GTiff byte response.


.. py:function:: transform2tuple(transform)

   Convert an affine transform to a tuple.

   :Parameters: **transform** (:class:`rio.Affine`) -- The affine transform to be converted

   :returns: :class:`tuple` -- The affine transform as a tuple (a, b, c, d, e, f)


:py:mod:`pynhd`
===============

.. py:module:: pynhd

.. autoapi-nested-parse::

   Top-level package for PyNHD.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   core/index.rst
   network_tools/index.rst
   nhdplus_derived/index.rst
   pynhd/index.rst


Package Contents
----------------

















:py:mod:`pynhd.pynhd`
=====================

.. py:module:: pynhd.pynhd

.. autoapi-nested-parse::

   Access NLDI and WaterData databases.



Module Contents
---------------

.. py:class:: NHD(layer, outfields = '*', crs = DEF_CRS, expire_after = EXPIRE, disable_caching = False)



   Access National Hydrography Dataset (NHD), both meduim and high resolution.

   .. rubric:: Notes

   For more info visit: https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer

   :Parameters: * **layer** (:class:`str`, *optional*) -- A valid service layer. Layer names with ``_hr`` are high resolution and
                  ``_mr`` are medium resolution. Also, layer names with ``_nonconus`` are for
                  non-conus areas, i.e., Alaska, Hawaii, Puerto Rico, the Virgin Islands , and
                  the Pacific Islands. Valid layers are:

                  - ``point``
                  - ``point_event``
                  - ``line_hr``
                  - ``flow_direction``
                  - ``flowline_mr``
                  - ``flowline_hr_nonconus``
                  - ``flowline_hr``
                  - ``area_mr``
                  - ``area_hr_nonconus``
                  - ``area_hr``
                  - ``waterbody_mr``
                  - ``waterbody_hr_nonconus``
                  - ``waterbody_hr``
                * **outfields** (:class:`str` or :class:`list`, *optional*) -- Target field name(s), default to "*" i.e., all the fields.
                * **crs** (:class:`str`, *optional*) -- Target spatial reference, default to ``EPSG:4326``
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.


.. py:class:: NHDPlusHR(layer, outfields = '*', crs = DEF_CRS, service = 'hydro')



   Access NHDPlus HR database through the National Map ArcGISRESTful.

   :Parameters: * **layer** (:class:`str`) -- A valid service layer. To see a list of available layers instantiate the class
                  with passing an empty string like so ``NHDPlusHR("")``.
                * **outfields** (:class:`str` or :class:`list`, *optional*) -- Target field name(s), default to "*" i.e., all the fields.
                * **crs** (:class:`str`, *optional*) -- Target spatial reference, default to EPSG:4326
                * **service** (:class:`str`, *optional*) -- Name of the web service to use, defaults to hydro. Supported web services are:

                  * hydro: https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer
                  * edits: https://edits.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/NHDPlus_HR/MapServer


.. py:class:: NLDI(expire_after = EXPIRE, disable_caching = False)

   Access the Hydro Network-Linked Data Index (NLDI) service.

   :Parameters: * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: comid_byloc(self, coords, loc_crs = DEF_CRS)

      Get the closest ComID(s) based on coordinates.

      :Parameters: * **coords** (:class:`tuple` or :class:`list`) -- A tuple of length two (x, y) or a list of them.
                   * **loc_crs** (:class:`str`, *optional*) -- The spatial reference of the input coordinate, defaults to EPSG:4326.

      :returns: :class:`geopandas.GeoDataFrame` or :class:`(geopandas.GeoDataFrame`, :class:`list)` -- NLDI indexed ComID(s) in EPSG:4326. If some coords don't return any ComID
                a list of missing coords are returned as well.


   .. py:method:: get_basins(self, station_ids, split_catchment = False)

      Get basins for a list of station IDs.

      :Parameters: * **station_ids** (:class:`str` or :class:`list`) -- USGS station ID(s).
                   * **split_catchment** (:class:`bool`, *optional*) -- If True, split the basin at the watershed outlet location. Default to False.

      :returns: :class:`geopandas.GeoDataFrame` or :class:`(geopandas.GeoDataFrame`, :class:`list)` -- NLDI indexed basins in EPSG:4326. If some IDs don't return any features
                a list of missing ID(s) are returned as well.


   .. py:method:: get_validchars(self, char_type)

      Get all the available characteristics IDs for a given characteristics type.


   .. py:method:: getcharacteristic_byid(self, comids, char_type, char_ids = 'all', values_only = True)

      Get characteristics using a list ComIDs.

      :Parameters: * **comids** (:class:`str` or :class:`list`) -- The ID of the feature.
                   * **char_type** (:class:`str`) -- Type of the characteristic. Valid values are ``local`` for
                     individual reach catchments, ``tot`` for network-accumulated values
                     using total cumulative drainage area and ``div`` for network-accumulated values
                     using divergence-routed.
                   * **char_ids** (:class:`str` or :class:`list`, *optional*) -- Name(s) of the target characteristics, default to all.
                   * **values_only** (:class:`bool`, *optional*) -- Whether to return only ``characteristic_value`` as a series, default to True.
                     If is set to False, ``percent_nodata`` is returned as well.

      :returns: :class:`pandas.DataFrame` or :class:`tuple` of :class:`pandas.DataFrame` -- Either only ``characteristic_value`` as a dataframe or
                or if ``values_only`` is Fale return ``percent_nodata`` as well.


   .. py:method:: getfeature_byid(self, fsource, fid)

      Get feature(s) based ID(s).

      :Parameters: * **fsource** (:class:`str`) -- The name of feature(s) source. The valid sources are:
                     comid, huc12pp, nwissite, wade, wqp
                   * **fid** (:class:`str` or :class:`list` of :class:`str`) -- Feature ID(s).

      :returns: :class:`geopandas.GeoDataFrame` or :class:`(geopandas.GeoDataFrame`, :class:`list)` -- NLDI indexed features in EPSG:4326. If some IDs don't return any features
                a list of missing ID(s) are returned as well.


   .. py:method:: navigate_byid(self, fsource, fid, navigation, source, distance = 500)

      Navigate the NHDPlus database from a single feature id up to a distance.

      :Parameters: * **fsource** (:class:`str`) -- The name of feature source. The valid sources are:
                     ``comid``, ``huc12pp``, ``nwissite``, ``wade``, ``WQP``.
                   * **fid** (:class:`str`) -- The ID of the feature.
                   * **navigation** (:class:`str`) -- The navigation method.
                   * **source** (:class:`str`, *optional*) -- Return the data from another source after navigating
                     the features using fsource, defaults to None.
                   * **distance** (:class:`int`, *optional*) -- Limit the search for navigation up to a distance in km,
                     defaults is 500 km. Note that this is an expensive request so you
                     have be mindful of the value that you provide. The value must be
                     between 1 to 9999 km.

      :returns: :class:`geopandas.GeoDataFrame` -- NLDI indexed features in EPSG:4326.


   .. py:method:: navigate_byloc(self, coords, navigation = None, source = None, loc_crs = DEF_CRS, distance = 500)

      Navigate the NHDPlus database from a coordinate.

      :Parameters: * **coords** (:class:`tuple`) -- A tuple of length two (x, y).
                   * **navigation** (:class:`str`, *optional*) -- The navigation method, defaults to None which throws an exception
                     if comid_only is False.
                   * **source** (:class:`str`, *optional*) -- Return the data from another source after navigating
                     the features using fsource, defaults to None which throws an exception
                     if comid_only is False.
                   * **loc_crs** (:class:`str`, *optional*) -- The spatial reference of the input coordinate, defaults to EPSG:4326.
                   * **distance** (:class:`int`, *optional*) -- Limit the search for navigation up to a distance in km,
                     defaults to 500 km. Note that this is an expensive request so you
                     have be mindful of the value that you provide. If you want to get
                     all the available features you can pass a large distance like 9999999.

      :returns: :class:`geopandas.GeoDataFrame` -- NLDI indexed features in EPSG:4326.



.. py:class:: PyGeoAPI

   Access `PyGeoAPI <https://labs.waterdata.usgs.gov/api/nldi/pygeoapi>`__ service.

   :Parameters: * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: cross_section(self, coord, width, numpts, crs = DEF_CRS)

      Return a GeoDataFrame from the xsatpoint service.

      :Parameters: * **coord** (:class:`tuple`) -- The coordinate of the point to extract the cross-section as a tuple,e.g., (lon, lat).
                   * **width** (:class:`float`) -- The width of the cross-section in meters.
                   * **numpts** (:class:`int`) -- The number of points to extract the cross-section from the DEM.
                   * **crs** (:class:`str`, *optional*) -- The coordinate reference system of the coordinates, defaults to EPSG:4326.

      :returns: :class:`geopandas.GeoDataFrame` -- A GeoDataFrame containing the cross-section at the requested point.

      .. rubric:: Examples

      >>> from pynhd import PyGeoAPI
      >>> pygeoapi = PyGeoAPI()
      >>> gdf = pygeoapi.cross_section((-103.80119, 40.2684), width=1000.0, numpts=101, crs=DEF_CRS)  # doctest: +SKIP
      >>> print(gdf.iloc[-1, 1])  # doctest: +SKIP
      1000.0


   .. py:method:: elevation_profile(self, coords, numpts, dem_res, crs = DEF_CRS)

      Return a GeoDataFrame from the xsatendpts service.

      :Parameters: * **coords** (:class:`list`) -- A list of two coordinates to trace as a list of tuples,e.g., [(lon, lat), (lon, lat)].
                   * **numpts** (:class:`int`) -- The number of points to extract the elevation profile from the DEM.
                   * **dem_res** (:class:`int`) -- The target resolution for requesting the DEM from 3DEP service.
                   * **crs** (:class:`str`, *optional*) -- The coordinate reference system of the coordinates, defaults to EPSG:4326.

      :returns: :class:`geopandas.GeoDataFrame` -- A GeoDataFrame containing the elevation profile along the requested endpoints.

      .. rubric:: Examples

      >>> from pynhd import PyGeoAPI
      >>> pygeoapi = PyGeoAPI()
      >>> gdf = pygeoapi.elevation_profile(
      ...     [(-103.801086, 40.26772), (-103.80097, 40.270568)], numpts=101, dem_res=1, crs=DEF_CRS
      ... )  # doctest: +SKIP
      >>> print(gdf.iloc[-1, 1])  # doctest: +SKIP
      411.5906


   .. py:method:: flow_trace(self, coord, crs = DEF_CRS, raindrop = False, direction = 'down')

      Return a GeoDataFrame from the flowtrace service.

      :Parameters: * **coord** (:class:`tuple`) -- The coordinate of the point to trace as a tuple,e.g., (lon, lat).
                   * **crs** (:class:`str`) -- The coordinate reference system of the coordinates, defaults to EPSG:4326.
                   * **raindrop** (:class:`bool`, *optional*) -- If True, use raindrop-based flowpaths, i.e. use raindrop trace web service
                     with direction set to "none", defaults to False.
                   * **direction** (:class:`str`, *optional*) -- The direction of flowpaths, either "down", "up", or "none". Defaults to "down".

      :returns: :class:`geopandas.GeoDataFrame` -- A GeoDataFrame containing the traced flowline.

      .. rubric:: Examples

      >>> from pynhd import PyGeoAPI
      >>> pygeoapi = PyGeoAPI()
      >>> gdf = pygeoapi.flow_trace(
      ...     (1774209.63, 856381.68), crs="ESRI:102003", raindrop=False, direction="none"
      ... )  # doctest: +SKIP
      >>> print(gdf.comid.iloc[0])  # doctest: +SKIP
      22294818


   .. py:method:: split_catchment(self, coord, crs = DEF_CRS, upstream = False)

      Return a GeoDataFrame from the splitcatchment service.

      :Parameters: * **coord** (:class:`tuple`) -- The coordinate of the point to trace as a tuple,e.g., (lon, lat).
                   * **crs** (:class:`str`, *optional*) -- The coordinate reference system of the coordinates, defaults to EPSG:4326.
                   * **upstream** (:class:`bool`, *optional*) -- If True, return all upstream catchments rather than just the local catchment,
                     defaults to False.

      :returns: :class:`geopandas.GeoDataFrame` -- A GeoDataFrame containing the local catchment or the entire upstream catchments.

      .. rubric:: Examples

      >>> from pynhd import PyGeoAPI
      >>> pygeoapi = PyGeoAPI()
      >>> gdf = pygeoapi.split_catchment((-73.82705, 43.29139), crs=DEF_CRS, upstream=False)  # doctest: +SKIP
      >>> print(gdf.catchmentID.iloc[0])  # doctest: +SKIP
      22294818



.. py:class:: WaterData(layer, crs = DEF_CRS)

   Access to `Water Data <https://labs.waterdata.usgs.gov/geoserver>`__ service.

   :Parameters: * **layer** (:class:`str`) -- A valid layer from the WaterData service. Valid layers are:
                  ``nhdarea``, ``nhdwaterbody``, ``catchmentsp``, ``nhdflowline_network``
                  ``gagesii``, ``huc08``, ``huc12``, ``huc12agg``, and ``huc12all``. Note that
                  the layers' worksapce for the Water Data service is ``wmadata`` which will
                  be added to the given ``layer`` argument if it is not provided.
                * **crs** (:class:`str`, *optional*) -- The target spatial reference system, defaults to ``epsg:4326``.

   .. py:method:: bybox(self, bbox, box_crs = DEF_CRS)

      Get features within a bounding box.


   .. py:method:: bydistance(self, coords, distance, loc_crs = DEF_CRS)

      Get features within a radius (in meters) of a point.


   .. py:method:: byfilter(self, cql_filter, method = 'GET')

      Get features based on a CQL filter.


   .. py:method:: bygeom(self, geometry, geo_crs = DEF_CRS, xy = True, predicate = 'INTERSECTS')

      Get features within a geometry.

      :Parameters: * **geometry** (:class:`shapely.geometry`) -- The input geometry
                   * **geo_crs** (:class:`str`, *optional*) -- The CRS of the input geometry, default to epsg:4326.
                   * **xy** (:class:`bool`, *optional*) -- Whether axis order of the input geometry is xy or yx.
                   * **predicate** (:class:`str`, *optional*) -- The geometric prediacte to use for requesting the data, defaults to
                     INTERSECTS. Valid predicates are:
                     ``EQUALS``, ``DISJOINT``, ``INTERSECTS``, ``TOUCHES``, ``CROSSES``, ``WITHIN``
                     ``CONTAINS``, ``OVERLAPS``, ``RELATE``, ``BEYOND``

      :returns: :class:`geopandas.GeoDataFrame` -- The requested features in the given geometry.


   .. py:method:: byid(self, featurename, featureids)

      Get features based on IDs.



:py:mod:`pynhd.network_tools`
=============================

.. py:module:: pynhd.network_tools

.. autoapi-nested-parse::

   Access NLDI and WaterData databases.



Module Contents
---------------

.. py:function:: flowline_xsection(flw, distance, width)

   Get cross-section of a river network at a given spacing.

   :Parameters: * **flw** (:class:`geopandas.GeoDataFrame`) -- A dataframe with ``geometry`` and ``comid`` columns and CRS attribute.
                * **distance** (:class:`float`) -- The distance between two consecutive cross-sections.
                * **width** (:class:`float`) -- The width of the cross-section.

   :returns: :class:`geopandas.GeoDataFrame` -- A dataframe with two columns: ``geometry`` and ``comid``. The ``geometry``
             column contains the cross-section of the river network and the ``comid``
             column contains the corresponding ``comid`` from the input dataframe.
             Note that each ``comid`` can have multiple cross-sections depending on
             the given spacing distance.


.. py:function:: network_xsection(flw, distance, width)

   Get cross-section of a river network at a given spacing.

   :Parameters: * **flw** (:class:`geopandas.GeoDataFrame`) -- A dataframe with ``geometry`` and ``comid`` columns and CRS attribute.
                * **distance** (:class:`float`) -- The distance between two consecutive cross-sections.
                * **width** (:class:`float`) -- The width of the cross-section.

   :returns: :class:`geopandas.GeoDataFrame` -- A dataframe with two columns: ``geometry`` and ``comid``. The ``geometry``
             column contains the cross-section of the river network and the ``comid``
             column contains the corresponding ``comid`` from the input dataframe.
             Note that each ``comid`` can have multiple cross-sections depending on
             the given spacing distance.


.. py:function:: nhdflw2nx(flowlines, id_col = 'comid', toid_col = 'tocomid', edge_attr = None)

   Convert NHDPlus flowline database to networkx graph.

   :Parameters: * **flowlines** (:class:`geopandas.GeoDataFrame`) -- NHDPlus flowlines.
                * **id_col** (:class:`str`, *optional*) -- Name of the column containing the node ID, defaults to "comid".
                * **toid_col** (:class:`str`, *optional*) -- Name of the column containing the downstream node ID, defaults to "tocomid".
                * **edge_attr** (:class:`str`, *optional*) -- Name of the column containing the edge attributes, defaults to ``None``.

   :returns: :class:`nx.DiGraph` -- Networkx directed graph of the NHDPlus flowlines.


.. py:function:: prepare_nhdplus(flowlines, min_network_size, min_path_length, min_path_size = 0, purge_non_dendritic = False, use_enhd_attrs = False, terminal2nan = True)

   Clean up and fix common issues of NHDPlus flowline database.

   Ported from `nhdplusTools <https://github.com/USGS-R/nhdplusTools>`__.

   :Parameters: * **flowlines** (:class:`geopandas.GeoDataFrame`) -- NHDPlus flowlines with at least the following columns:
                  ``comid``, ``lengthkm``, ``ftype``, ``terminalfl``, ``fromnode``, ``tonode``,
                  ``totdasqkm``, ``startflag``, ``streamorde``, ``streamcalc``, ``terminalpa``,
                  ``pathlength``, ``divergence``, ``hydroseq``, ``levelpathi``.
                * **min_network_size** (:class:`float`) -- Minimum size of drainage network in sqkm
                * **min_path_length** (:class:`float`) -- Minimum length of terminal level path of a network in km.
                * **min_path_size** (:class:`float`, *optional*) -- Minimum size of outlet level path of a drainage basin in km.
                  Drainage basins with an outlet drainage area smaller than
                  this value will be removed. Defaults to 0.
                * **purge_non_dendritic** (:class:`bool`, *optional*) -- Whether to remove non dendritic paths, defaults to False
                * **use_enhd_attrs** (:class:`bool`, *optional*) -- Whether to replace the attributes with the ENHD attributes, defaults to False.
                  For more information, see
                  `this <https://www.sciencebase.gov/catalog/item/60c92503d34e86b9389df1c9>`__.
                * **terminal2nan** (:class:`bool`, *optional*) -- Whether to replace the COMID of the terminal flowline of the network with NaN,
                  defaults to True. If False, the terminal COMID will be set from the
                  ENHD attributes i.e. use_enhd_attrs will be set to True.

   :returns: :class:`geopandas.GeoDataFrame` -- Cleaned up flowlines. Note that all column names are converted to lower case.


.. py:function:: topoogical_sort(flowlines, edge_attr = None)

   Topological sorting of a river network.

   :Parameters: * **flowlines** (:class:`pandas.DataFrame`) -- A dataframe with columns ID and toID
                * **edge_attr** (:class:`str` or :class:`list`, *optional*) -- Names of the columns in the dataframe to be used as edge attributes, defaults to None.

   :returns: :class:`(list`, dict , :class:`networkx.DiGraph)` -- A list of topologically sorted IDs, a dictionary
             with keys as IDs and values as its upstream nodes,
             and the generated networkx object. Note that the
             terminal node ID is set to pd.NA.


.. py:function:: vector_accumulation(flowlines, func, attr_col, arg_cols, id_col = 'comid', toid_col = 'tocomid')

   Flow accumulation using vector river network data.

   :Parameters: * **flowlines** (:class:`pandas.DataFrame`) -- A dataframe containing comid, tocomid, attr_col and all the columns
                  that ara required for passing to ``func``.
                * **func** (:class:`function`) -- The function that routes the flow in a single river segment.
                  Positions of the arguments in the function should be as follows:
                  ``func(qin, *arg_cols)``
                  ``qin`` is computed in this function and the rest are in the order
                  of the ``arg_cols``. For example, if ``arg_cols = ["slope", "roughness"]``
                  then the functions is called this way:
                  ``func(qin, slope, roughness)``
                  where slope and roughness are elemental values read from the flowlines.
                * **attr_col** (:class:`str`) -- The column name of the attribute being accumulated in the network.
                  The column should contain the initial condition for the attribute for
                  each river segment. It can be a scalar or an array (e.g., time series).
                * **arg_cols** (:class:`list` of :class:`strs`) -- List of the flowlines columns that contain all the required
                  data for a routing a single river segment such as slope, length,
                  lateral flow, etc.
                * **id_col** (:class:`str`, *optional*) -- Name of the flowlines column containing IDs, defaults to ``comid``
                * **toid_col** (:class:`str`, *optional*) -- Name of the flowlines column containing ``toIDs``, defaults to ``tocomid``

   :returns: :class:`pandas.Series` -- Accumulated flow for all the nodes. The dataframe is sorted from upstream
             to downstream (topological sorting). Depending on the given initial
             condition in the ``attr_col``, the outflow for each river segment can be
             a scalar or an array.


:py:mod:`pynhd.nhdplus_derived`
===============================

.. py:module:: pynhd.nhdplus_derived

.. autoapi-nested-parse::

   Access NLDI and WaterData databases.



Module Contents
---------------

.. py:function:: enhd_attrs(parquet_path = None, expire_after = EXPIRE, disable_caching = False)

   Get updated NHDPlus attributes from ENHD.

   .. rubric:: Notes

   This downloads a 140 MB ``parquet`` file from
   `here <https://www.sciencebase.gov/catalog/item/60c92503d34e86b9389df1c9>`__ .
   Although this dataframe does not include geometry, it can be linked to other geospatial
   NHDPlus dataframes through ComIDs.

   :Parameters: * **parquet_path** (:class:`str` or :class:`~~pathlib.Path`, *optional*) -- Path to a file with ``.parquet`` extension for storing the file, defaults to
                  ``./cache/enhd_attrs.parquet``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`pandas.DataFrame` -- A dataframe that includes ComID-level attributes for 2.7 million NHDPlus flowlines.


.. py:function:: nhd_fcode()

   Get all the NHDPlus FCodes.


.. py:function:: nhdplus_attrs(name = None, parquet_path = None, expire_after = EXPIRE, disable_caching = False)

   Access NHDPlus V2.1 Attributes from ScienceBase over CONUS.

   More info can be found `here <https://www.sciencebase.gov/catalog/item/5669a79ee4b08895842a1d47>`_.

   :Parameters: * **name** (:class:`str`, *optional*) -- Name of the NHDPlus attribute, defaults to None which returns a dataframe containing
                  metadata of all the available attributes in the database.
                * **parquet_path** (:class:`str` or :class:`~~pathlib.Path`, *optional*) -- Path to a file with ``.parquet`` extension for saving the processed to disk for
                  later use. Defaults to ``./cache/nhdplus_attrs.parquet``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`pandas.DataFrame` -- Either a dataframe containing the database metadata or the requested attribute over CONUS.


.. py:function:: nhdplus_vaa(parquet_path = None, expire_after = EXPIRE, disable_caching = False)

   Get NHDPlus Value Added Attributes with ComID-level roughness and slope values.

   .. rubric:: Notes

   This downloads a 200 MB ``parquet`` file from
   `here <https://www.hydroshare.org/resource/6092c8a62fac45be97a09bfd0b0bf726>`__ .
   Although this dataframe does not include geometry, it can be linked to other geospatial
   NHDPlus dataframes through ComIDs.

   :Parameters: * **parquet_path** (:class:`str` or :class:`~~pathlib.Path`, *optional*) -- Path to a file with ``.parquet`` extension for storing the file, defaults to
                  ``./cache/nldplus_vaa.parquet``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`pandas.DataFrame` -- A dataframe that includes ComID-level attributes for 2.7 million NHDPlus flowlines.

   .. rubric:: Examples

   >>> vaa = nhdplus_vaa() # doctest: +SKIP
   >>> print(vaa.slope.max()) # doctest: +SKIP
   4.6


:py:mod:`pynhd.core`
====================

.. py:module:: pynhd.core

.. autoapi-nested-parse::

   Base classes for PyNHD functions.



Module Contents
---------------

.. py:class:: AGRBase(base_url, layer = None, outfields = '*', crs = DEF_CRS, outformat = 'json', expire_after = EXPIRE, disable_caching = False)

   Base class for getting geospatial data from a ArcGISRESTful service.

   :Parameters: * **base_url** (:class:`str`, *optional*) -- The ArcGIS RESTful service url. The URL must either include a layer number
                  after the last ``/`` in the url or the target layer must be passed as an argument.
                * **layer** (:class:`str`, *optional*) -- A valid service layer. To see a list of available layers instantiate the class
                  without passing any argument.
                * **outfields** (:class:`str` or :class:`list`, *optional*) -- Target field name(s), default to "*" i.e., all the fields.
                * **crs** (:class:`str`, *optional*) -- Target spatial reference, default to EPSG:4326
                * **outformat** (:class:`str`, *optional*) -- One of the output formats offered by the selected layer. If not correct
                  a list of available formats is shown, defaults to ``json``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: bygeom(self, geom, geo_crs = DEF_CRS, sql_clause = '', distance = None, return_m = False, return_geom = True)

      Get feature within a geometry that can be combined with a SQL where clause.

      :Parameters: * **geom** (:class:`Polygon` or :class:`tuple`) -- A geometry (Polygon) or bounding box (tuple of length 4).
                   * **geo_crs** (:class:`str`) -- The spatial reference of the input geometry.
                   * **sql_clause** (:class:`str`, *optional*) -- A valid SQL 92 WHERE clause, defaults to an empty string.
                   * **distance** (:class:`int`, *optional*) -- The buffer distance for the input geometries in meters, default to None.
                   * **return_m** (:class:`bool`, *optional*) -- Whether to activate the Return M (measure) in the request, defaults to False.
                   * **return_geom** (:class:`bool`, *optional*) -- Whether to return the geometry of the feature, defaults to ``True``.

      :returns: :class:`geopandas.GeoDataFrame` -- The requested features as a GeoDataFrame.


   .. py:method:: byids(self, field, fids, return_m = False, return_geom = True)

      Get features based on a list of field IDs.

      :Parameters: * **field** (:class:`str`) -- Name of the target field that IDs belong to.
                   * **fids** (:class:`str` or :class:`list`) -- A list of target field ID(s).
                   * **return_m** (:class:`bool`) -- Whether to activate the Return M (measure) in the request, defaults to False.
                   * **return_geom** (:class:`bool`, *optional*) -- Whether to return the geometry of the feature, defaults to ``True``.

      :returns: :class:`geopandas.GeoDataFrame` -- The requested features as a GeoDataFrame.


   .. py:method:: bysql(self, sql_clause, return_m = False, return_geom = True)

      Get feature IDs using a valid SQL 92 WHERE clause.

      .. rubric:: Notes

      Not all web services support this type of query. For more details look
      `here <https://developers.arcgis.com/rest/services-reference/query-feature-service-.htm#ESRI_SECTION2_07DD2C5127674F6A814CE6C07D39AD46>`__

      :Parameters: * **sql_clause** (:class:`str`) -- A valid SQL 92 WHERE clause.
                   * **return_m** (:class:`bool`) -- Whether to activate the measure in the request, defaults to False.
                   * **return_geom** (:class:`bool`, *optional*) -- Whether to return the geometry of the feature, defaults to ``True``.

      :returns: :class:`geopandas.GeoDataFrame` -- The requested features as a GeoDataFrame.


   .. py:method:: get_validlayers(self, url)

      Get a list of valid layers.

      :Parameters: **url** (:class:`str`) -- The URL of the ArcGIS REST service.

      :returns: :class:`dict` -- A dictionary of valid layers.



.. py:class:: ScienceBase

   Access and explore files on ScienceBase.

   :Parameters: * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: get_children(self, item)

      Get children items of an item.


   .. py:method:: get_file_urls(self, item)

      Get download and meta URLs of all the available files for an item.



.. py:function:: stage_nhdplus_attrs(parquet_path = None, expire_after = EXPIRE, disable_caching = False)

   Stage the NHDPlus Attributes database and save to nhdplus_attrs.parquet.

   More info can be found `here <https://www.sciencebase.gov/catalog/item/5669a79ee4b08895842a1d47>`_.

   :Parameters: * **parquet_path** (:class:`str` or :class:`~~pathlib.Path`) -- Path to a file with ``.parquet`` extension for saving the processed to disk for
                  later use.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`pandas.DataFrame` -- The staged data as a DataFrame.


:py:mod:`pygeoogc`
==================

.. py:module:: pygeoogc

.. autoapi-nested-parse::

   Top-level package for PyGeoOGC.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   core/index.rst
   pygeoogc/index.rst
   utils/index.rst


Package Contents
----------------






:py:mod:`pygeoogc.pygeoogc`
===========================

.. py:module:: pygeoogc.pygeoogc

.. autoapi-nested-parse::

   Base classes and function for REST, WMS, and WMF services.



Module Contents
---------------

.. py:class:: ArcGISRESTful(base_url, layer = None, outformat = 'geojson', outfields = '*', crs = DEF_CRS, max_workers = 1, verbose = False, disable_retry = False, expire_after = EXPIRE, disable_caching = False)

   Access to an ArcGIS REST service.

   .. rubric:: Notes

   By default, all retrieval methods retry to get the missing feature IDs,
   if there are any. You can disable this behavior by setting ``disable_retry``
   to ``True``. If there are any missing feature IDs after the retry,
   they are saved to a text file, path of which can be accessed by
   ``self.client.failed_path``.

   :Parameters: * **base_url** (:class:`str`, *optional*) -- The ArcGIS RESTful service url. The URL must either include a layer number
                  after the last ``/`` in the url or the target layer must be passed as an argument.
                * **layer** (:class:`int`, *optional*) -- Target layer number, defaults to None. If None layer number must be included as after
                  the last ``/`` in ``base_url``.
                * **outformat** (:class:`str`, *optional*) -- One of the output formats offered by the selected layer. If not correct
                  a list of available formats is shown, defaults to ``geojson``.
                * **outfields** (:class:`str` or :class:`list`) -- The output fields to be requested. Setting ``*`` as outfields requests
                  all the available fields which is the default behaviour.
                * **crs** (:class:`str`, *optional*) -- The spatial reference of the output data, defaults to ``epsg:4326``.
                * **max_workers** (:class:`int`, *optional*) -- Number of simultaneous download, default to 1, i.e., no threading. Note
                  that some services might face issues when several requests are sent
                  simultaneously and will return the requests partially. It's recommended
                  to avoid using too many workers unless you are certain the web service
                  can handle it.
                * **verbose** (:class:`bool`, *optional*) -- If True, prints information about the requests and responses,
                  defaults to False.
                * **disable_retry** (:class:`bool`, *optional*) -- If ``True`` in case there are any failed queries, no retrying attempts
                  is done and object IDs of the failed requests is saved to a text file
                  which its path can be accessed via ``self.client.failed_path``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: get_features(self, featureids, return_m = False, return_geom = True)

      Get features based on the feature IDs.

      :Parameters: * **featureids** (:class:`list`) -- List of feature IDs.
                   * **return_m** (:class:`bool`, *optional*) -- Whether to activate the Return M (measure) in the request,
                     defaults to ``False``.
                   * **return_geom** (:class:`bool`, *optional*) -- Whether to return the geometry of the feature, defaults to ``True``.

      :returns: :class:`dict` -- (Geo)json response from the web service.


   .. py:method:: oids_byfield(self, field, ids)

      Get Object IDs based on a list of field IDs.

      :Parameters: * **field** (:class:`str`) -- Name of the target field that IDs belong to.
                   * **ids** (:class:`str` or :class:`list`) -- A list of target ID(s).

      :returns: :class:`list` of :class:`tuples` -- A list of feature IDs partitioned by ``self.max_nrecords``.


   .. py:method:: oids_bygeom(self, geom, geo_crs = DEF_CRS, spatial_relation = 'esriSpatialRelIntersects', sql_clause = None, distance = None)

      Get feature IDs within a geometry that can be combined with a SQL where clause.

      :Parameters: * **geom** (:class:`LineString`, :class:`Polygon`, :class:`Point`, :class:`MultiPoint`, :class:`tuple`, or :class:`list` of :class:`tuples`) -- A geometry (LineString, Polygon, Point, MultiPoint), tuple of length two
                     (``(x, y)``), a list of tuples of length 2 (``[(x, y), ...]``), or bounding box
                     (tuple of length 4 (``(xmin, ymin, xmax, ymax)``)).
                   * **geo_crs** (:class:`str` or :class:`pyproj.CRS`) -- The spatial reference of the input geometry.
                   * **spatial_relation** (:class:`str`, *optional*) -- The spatial relationship to be applied on the input geometry
                     while performing the query. If not correct a list of available options is shown.
                     It defaults to ``esriSpatialRelIntersects``. Valid predicates are:

                     * ``esriSpatialRelIntersects``
                     * ``esriSpatialRelContains``
                     * ``esriSpatialRelCrosses``
                     * ``esriSpatialRelEnvelopeIntersects``
                     * ``esriSpatialRelIndexIntersects``
                     * ``esriSpatialRelOverlaps``
                     * ``esriSpatialRelTouches``
                     * ``esriSpatialRelWithin``
                     * ``esriSpatialRelRelation``
                   * **sql_clause** (:class:`str`, *optional*) -- Valid SQL 92 WHERE clause, default to None.
                   * **distance** (:class:`int`, *optional*) -- Buffer distance in meters for the input geometries, default to None.

      :returns: :class:`list` of :class:`tuples` -- A list of feature IDs partitioned by ``self.max_nrecords``.


   .. py:method:: oids_bysql(self, sql_clause)

      Get feature IDs using a valid SQL 92 WHERE clause.

      .. rubric:: Notes

      Not all web services support this type of query. For more details look
      `here <https://developers.arcgis.com/rest/services-reference/query-feature-service-.htm#ESRI_SECTION2_07DD2C5127674F6A814CE6C07D39AD46>`__.

      :Parameters: **sql_clause** (:class:`str`) -- A valid SQL 92 WHERE clause.

      :returns: :class:`list` of :class:`tuples` -- A list of feature IDs partitioned by ``self.max_nrecords``.


   .. py:method:: partition_oids(self, oids)

      Partition feature IDs based on ``self.max_nrecords``.

      :Parameters: **oids** (:class:`list` of :class:`int` or :class:`int`) -- A list of feature ID(s).

      :returns: :class:`list` of :class:`tuples` -- A list of feature IDs partitioned by ``self.max_nrecords``.



.. py:class:: ServiceURL

   Base URLs of the supported services.

   .. py:method:: http(self)
      :property:

      Read HTTP URLs from the source yml file.


   .. py:method:: restful(self)
      :property:

      Read RESTful URLs from the source yml file.


   .. py:method:: wfs(self)
      :property:

      Read WFS URLs from the source yml file.


   .. py:method:: wms(self)
      :property:

      Read WMS URLs from the source yml file.



.. py:class:: WFS(url, layer = None, outformat = None, version = '2.0.0', crs = DEF_CRS, read_method = 'json', max_nrecords = 1000, validation = True, expire_after = EXPIRE, disable_caching = False)



   Data from any WFS service within a geometry or by featureid.

   :Parameters: * **url** (:class:`str`) -- The base url for the WFS service, for examples:
                  https://hazards.fema.gov/nfhl/services/public/NFHL/MapServer/WFSServer
                * **layer** (:class:`str`) -- The layer from the service to be downloaded, defaults to None which throws
                  an error and includes all the available layers offered by the service.
                * **outformat** (:class:`str`) --

                  The data format to request for data from the service, defaults to None which
                   throws an error and includes all the available format offered by the service.
                * **version** (:class:`str`, *optional*) -- The WFS service version which should be either 1.0.0, 1.1.0, or 2.0.0.
                  Defaults to 2.0.0.
                * **crs** (:class:`str`, *optional*) -- The spatial reference system to be used for requesting the data, defaults to
                  ``epsg:4326``.
                * **read_method** (:class:`str`, *optional*) -- Method for reading the retrieved data, defaults to ``json``. Valid options are
                  ``json``, ``binary``, and ``text``.
                * **max_nrecords** (:class:`int`, *optional*) -- The maximum number of records in a single request to be retrieved from the service,
                  defaults to 1000. If the number of records requested is greater than this value,
                  it will be split into multiple requests.
                * **validation** (:class:`bool`, *optional*) -- Validate the input arguments from the WFS service, defaults to True. Set this
                  to False if you are sure all the WFS settings such as layer and crs are correct
                  to avoid sending extra requests.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: getfeature_bybox(self, bbox, box_crs = DEF_CRS, always_xy = False)

      Get data from a WFS service within a bounding box.

      :Parameters: * **bbox** (:class:`tuple`) -- A bounding box for getting the data: [west, south, east, north]
                   * **box_crs** (:class:`str`, or :class:`pyproj.CRS`, *optional*) -- The spatial reference system of the input bbox, defaults to
                     ``epsg:4326``.
                   * **always_xy** (:class:`bool`, *optional*) -- Whether to always use xy axis order, defaults to False. Some services change the axis
                     order from xy to yx, following the latest WFS version specifications but some don't.
                     If the returned value does not have any geometry, it indicates that most probably the
                     axis order does not match. You can set this to True in that case.

      :returns: :class:`str` or :class:`bytes` or :class:`dict` -- WFS query response within a bounding box.


   .. py:method:: getfeature_byfilter(self, cql_filter, method = 'GET')

      Get features based on a valid CQL filter.

      .. rubric:: Notes

      The validity of the input CQL expression is user's responsibility since
      the function does not perform any checks and just sends a request using
      the input filter.

      :Parameters: * **cql_filter** (:class:`str`) -- A valid CQL filter expression.
                   * **method** (:class:`str`) -- The request method, could be GET or POST (for long filters).

      :returns: :class:`str` or :class:`bytes` or :class:`dict` -- WFS query response


   .. py:method:: getfeature_bygeom(self, geometry, geo_crs = DEF_CRS, always_xy = False, predicate = 'INTERSECTS')

      Get features based on a geometry.

      :Parameters: * **geometry** (:class:`shapely.geometry`) -- The input geometry
                   * **geo_crs** (:class:`str`, or :class:`pyproj.CRS`, *optional*) -- The CRS of the input geometry, default to ``epsg:4326``.
                   * **always_xy** (:class:`bool`, *optional*) -- Whether to always use xy axis order, defaults to False. Some services change the axis
                     order from xy to yx, following the latest WFS version specifications but some don't.
                     If the returned value does not have any geometry, it indicates that most probably the
                     axis order does not match. You can set this to True in that case.
                   * **predicate** (:class:`str`, *optional*) -- The geometric predicate to use for requesting the data, defaults to ``INTERSECTS``.
                     Valid predicates are:

                     * ``EQUALS``
                     * ``DISJOINT``
                     * ``INTERSECTS``
                     * ``TOUCHES``
                     * ``CROSSES``
                     * ``WITHIN``
                     * ``CONTAINS``
                     * ``OVERLAPS``
                     * ``RELATE``
                     * ``BEYOND``

      :returns: :class:`str` or :class:`bytes` or :class:`dict` -- WFS query response based on the given geometry.


   .. py:method:: getfeature_byid(self, featurename, featureids)

      Get features based on feature IDs.

      :Parameters: * **featurename** (:class:`str`) -- The name of the column for searching for feature IDs.
                   * **featureids** (:class:`str` or :class:`list`) -- The feature ID(s).

      :returns: :class:`str` or :class:`bytes` or :class:`dict` -- WMS query response.



.. py:class:: WMS(url, layers, outformat, version = '1.3.0', crs = DEF_CRS, validation = True, expire_after = EXPIRE, disable_caching = False)



   Get data from a WMS service within a geometry or bounding box.

   :Parameters: * **url** (:class:`str`) -- The base url for the WMS service e.g., https://www.mrlc.gov/geoserver/mrlc_download/wms
                * **layers** (:class:`str` or :class:`list`) -- A layer or a list of layers from the service to be downloaded. You can pass an empty
                  string to get a list of available layers.
                * **outformat** (:class:`str`) -- The data format to request for data from the service. You can pass an empty
                  string to get a list of available output formats.
                * **crs** (:class:`str`, *optional*) -- The spatial reference system to be used for requesting the data, defaults to
                  ``epsg:4326``.
                * **version** (:class:`str`, *optional*) -- The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.
                * **validation** (:class:`bool`, *optional*) -- Validate the input arguments from the WMS service, defaults to True. Set this
                  to False if you are sure all the WMS settings such as layer and crs are correct
                  to avoid sending extra requests.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: getmap_bybox(self, bbox, resolution, box_crs = DEF_CRS, always_xy = False, max_px = 8000000, kwargs = None)

      Get data from a WMS service within a geometry or bounding box.

      :Parameters: * **bbox** (:class:`tuple`) -- A bounding box for getting the data.
                   * **resolution** (:class:`float`) -- The output resolution in meters. The width and height of output are computed in pixel
                     based on the geometry bounds and the given resolution.
                   * **box_crs** (:class:`str`, or :class:`pyproj.CRS`, *optional*) -- The spatial reference system of the input bbox, defaults to
                     ``epsg:4326``.
                   * **always_xy** (:class:`bool`, *optional*) -- Whether to always use xy axis order, defaults to False. Some services change the axis
                     order from xy to yx, following the latest WFS version specifications but some don't.
                     If the returned value does not have any geometry, it indicates that most probably the
                     axis order does not match. You can set this to True in that case.
                   * **max_px** (:class:`int`, :class:`opitonal`) -- The maximum allowable number of pixels (width x height) for a WMS requests,
                     defaults to 8 million based on some trial-and-error.
                   * **kwargs** (:class:`dict`, *optional*) -- Optional additional keywords passed as payload, defaults to None.
                     For example, ``{"styles": "default"}``.

      :returns: :class:`dict` -- A dict where the keys are the layer name and values are the returned response
                from the WMS service as bytes.



:py:mod:`pygeoogc.utils`
========================

.. py:module:: pygeoogc.utils

.. autoapi-nested-parse::

   Some utilities for PyGeoOGC.



Module Contents
---------------

.. py:class:: ESRIGeomQuery

   Generate input geometry query for ArcGIS RESTful services.

   :Parameters: * **geometry** (:class:`tuple` or :class:`sgeom.Polygon` or :class:`sgeom.Point` or :class:`sgeom.LineString`) -- The input geometry which can be a point (x, y), a list of points [(x, y), ...],
                  bbox (xmin, ymin, xmax, ymax), or a Shapely's sgeom.Polygon.
                * **wkid** (:class:`int`) -- The Well-known ID (WKID) of the geometry's spatial reference e.g., for EPSG:4326,
                  4326 should be passed. Check
                  `ArcGIS <https://developers.arcgis.com/rest/services-reference/geographic-coordinate-systems.htm>`__
                  for reference.

   .. py:method:: bbox(self)

      Query for a bbox.


   .. py:method:: multipoint(self)

      Query for a multi-point.


   .. py:method:: point(self)

      Query for a point.


   .. py:method:: polygon(self)

      Query for a polygon.


   .. py:method:: polyline(self)

      Query for a polyline.



.. py:class:: RetrySession(retries = 3, backoff_factor = 0.3, status_to_retry = (500, 502, 504), prefixes = ('https://', ), cache_name = None)

   Configures the passed-in session to retry on failed requests.

   The fails can be due to connection errors, specific HTTP response
   codes and 30X redirections. The code is was originally based on:
   https://github.com/bustawin/retry-requests

   :Parameters: * **retries** (:class:`int`, *optional*) -- The number of maximum retries before raising an exception, defaults to 5.
                * **backoff_factor** (:class:`float`, *optional*) -- A factor used to compute the waiting time between retries, defaults to 0.5.
                * **status_to_retry** (:class:`tuple`, *optional*) -- A tuple of status codes that trigger the reply behaviour, defaults to (500, 502, 504).
                * **prefixes** (:class:`tuple`, *optional*) -- The prefixes to consider, defaults to ("http://", "https://")
                * **cache_name** (:class:`str`, *optional*) -- Path to a folder for caching the session, default to None which uses
                  system's temp directory.

   .. py:method:: get(self, url, payload = None, headers = None)

      Retrieve data from a url by GET and return the Response.


   .. py:method:: post(self, url, payload = None, headers = None)

      Retrieve data from a url by POST and return the Response.



.. py:function:: bbox_decompose(bbox, resolution, box_crs = DEF_CRS, max_px = 8000000)

   Split the bounding box vertically for WMS requests.

   :Parameters: * **bbox** (:class:`tuple`) -- A bounding box; (west, south, east, north)
                * **resolution** (:class:`float`) -- The target resolution for a WMS request in meters.
                * **box_crs** (:class:`str`, *optional*) -- The spatial reference of the input bbox, default to EPSG:4326.
                * **max_px** (:class:`int`, :class:`opitonal`) -- The maximum allowable number of pixels (width x height) for a WMS requests,
                  defaults to 8 million based on some trial-and-error.

   :returns: :class:`list` of :class:`tuples` -- Each tuple includes the following elements:

             * Tuple of length 4 that represents a bounding box (west, south, east, north) of a cell,
             * A label that represents cell ID starting from bottom-left to top-right, for example a
               2x2 decomposition has the following labels::

               |---------|---------|
               |         |         |
               |   0_1   |   1_1   |
               |         |         |
               |---------|---------|
               |         |         |
               |   0_0   |   1_0   |
               |         |         |
               |---------|---------|

             * Raster width of a cell,
             * Raster height of a cell.


.. py:function:: bbox_resolution(bbox, resolution, bbox_crs = DEF_CRS)

   Image size of a bounding box WGS84 for a given resolution in meters.

   :Parameters: * **bbox** (:class:`tuple`) -- A bounding box in WGS84 (west, south, east, north)
                * **resolution** (:class:`float`) -- The resolution in meters
                * **bbox_crs** (:class:`str`, *optional*) -- The spatial reference of the input bbox, default to EPSG:4326.

   :returns: :class:`tuple` -- The width and height of the image


.. py:function:: check_bbox(bbox)

   Check if an input inbox is a tuple of length 4.


.. py:function:: check_response(resp)

   Extract error message from a response, if any.


.. py:function:: match_crs(geom, in_crs, out_crs)

   Reproject a geometry to another CRS.

   :Parameters: * **geom** (:class:`list` or :class:`tuple` or :class:`geometry`) -- Input geometry which could be a list of coordinates such as ``[(x1, y1), ...]``,
                  a bounding box like so ``(xmin, ymin, xmax, ymax)``, or any valid ``shapely``'s
                  geometry such as ``Polygon``, ``MultiPolygon``, etc..
                * **in_crs** (:class:`str`) -- Spatial reference of the input geometry
                * **out_crs** (:class:`str`) -- Target spatial reference

   :returns: :class:`same type as the input geometry` -- Transformed geometry in the target CRS.

   .. rubric:: Examples

   >>> from pygeoogc.utils import match_crs
   >>> from shapely.geometry import Point
   >>> point = Point(-7766049.665, 5691929.739)
   >>> match_crs(point, "epsg:3857", "epsg:4326").xy
   (array('d', [-69.7636111130079]), array('d', [45.44549114818127]))
   >>> bbox = (-7766049.665, 5691929.739, -7763049.665, 5696929.739)
   >>> match_crs(bbox, "epsg:3857", "epsg:4326")
   (-69.7636111130079, 45.44549114818127, -69.73666165448431, 45.47699468552394)
   >>> coords = [(-7766049.665, 5691929.739)]
   >>> match_crs(coords, "epsg:3857", "epsg:4326")
   [(-69.7636111130079, 45.44549114818127)]


.. py:function:: traverse_json(obj, path)

   Extract an element from a JSON file along a specified path.

   This function is based on `bcmullins <https://bcmullins.github.io/parsing-json-python/>`__.

   :Parameters: * **obj** (:class:`dict`) -- The input json dictionary
                * **path** (:class:`list`) -- The path to the requested element

   :returns: :class:`list` -- The items founds in the JSON

   .. rubric:: Examples

   >>> from pygeoogc.utils import traverse_json
   >>> data = [{
   ...     "employees": [
   ...         {"name": "Alice", "role": "dev", "nbr": 1},
   ...         {"name": "Bob", "role": "dev", "nbr": 2}],
   ...     "firm": {"name": "Charlie's Waffle Emporium", "location": "CA"},
   ... },]
   >>> traverse_json(data, ["employees", "name"])
   [['Alice', 'Bob']]


.. py:function:: valid_wms_crs(url)

   Get valid CRSs from a WMS service version 1.3.0.


.. py:function:: validate_crs(val)

   Validate a CRS.

   :Parameters: **val** (:class:`str` or :class:`int`) -- Input CRS.

   :returns: :class:`str` -- Validated CRS as a string.


:py:mod:`pygeoogc.core`
=======================

.. py:module:: pygeoogc.core

.. autoapi-nested-parse::

   Base classes and function for REST, WMS, and WMF services.



Module Contents
---------------

.. py:class:: ArcGISRESTfulBase(base_url, layer = None, outformat = 'geojson', outfields = '*', crs = DEF_CRS, max_workers = 1, verbose = False, disable_retry = False, expire_after = EXPIRE, disable_caching = False)

   Access to an ArcGIS REST service.

   :Parameters: * **base_url** (:class:`str`, *optional*) -- The ArcGIS RESTful service url. The URL must either include a layer number
                  after the last ``/`` in the url or the target layer must be passed as an
                  argument.
                * **layer** (:class:`int`, *optional*) -- Target layer number, defaults to None. If None layer number must be
                  included as after the last ``/`` in ``base_url``.
                * **outformat** (:class:`str`, *optional*) -- One of the output formats offered by the selected layer. If not correct
                  a list of available formats is shown, defaults to ``geojson``.
                  It defaults to ``esriSpatialRelIntersects``.
                * **outfields** (:class:`str` or :class:`list`) -- The output fields to be requested. Setting ``*`` as outfields requests
                  all the available fields which is the default setting.
                * **crs** (:class:`str`, *optional*) -- The spatial reference of the output data, defaults to EPSG:4326
                * **max_workers** (:class:`int`, *optional*) -- Max number of simultaneous requests, default to 2. Note
                  that some services might face issues when several requests are sent
                  simultaneously and will return the requests partially. It's recommended
                  to avoid using too many workers unless you are certain the web service
                  can handle it.
                * **verbose** (:class:`bool`, *optional*) -- If True, prints information about the requests and responses,
                  defaults to False.
                * **disable_retry** (:class:`bool`, *optional*) -- If ``True`` in case there are any failed queries, no retrying attempts
                  is done and object IDs of the failed requests is saved to a text file
                  which its path can be accessed via ``self.failed_path``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: esri_query(self, geom, geo_crs = DEF_CRS)

      Generate geometry queries based on ESRI template.


   .. py:method:: get_features(self, featureids, return_m = False, return_geom = True)

      Get features based on the feature IDs.

      :Parameters: * **featureids** (:class:`list`) -- List of feature IDs.
                   * **return_m** (:class:`bool`, *optional*) -- Whether to activate the Return M (measure) in the request,
                     defaults to ``False``.
                   * **return_geom** (:class:`bool`, *optional*) -- Whether to return the geometry of the feature, defaults to ``True``.

      :returns: :class:`dict` -- (Geo)json response from the web service.


   .. py:method:: get_response(self, url, payloads, method = 'GET')

      Send payload and get the response.


   .. py:method:: initialize_service(self)

      Initialize the RESTFul service.


   .. py:method:: partition_oids(self, oids)

      Partition feature IDs based on ``self.max_nrecords``.


   .. py:method:: retry_failed_requests(self)

      Retry failed requests.



.. py:class:: RESTValidator(__pydantic_self__, **data)



   Validate ArcGISRESTful inputs.

   :Parameters: * **base_url** (:class:`str`, *optional*) -- The ArcGIS RESTful service url. The URL must either include a layer number
                  after the last ``/`` in the url or the target layer must be passed as an argument.
                * **layer** (:class:`int`, *optional*) -- Target layer number, defaults to None. If None layer number must be included as after
                  the last ``/`` in ``base_url``.
                * **outformat** (:class:`str`, *optional*) -- One of the output formats offered by the selected layer. If not correct
                  a list of available formats is shown, defaults to ``geojson``.
                * **outfields** (:class:`str` or :class:`list`) -- The output fields to be requested. Setting ``*`` as outfields requests
                  all the available fields which is the default setting.
                * **crs** (:class:`str`, *optional*) -- The spatial reference of the output data, defaults to EPSG:4326
                * **max_workers** (:class:`int`, *optional*) -- Max number of simultaneous requests, default to 2. Note
                  that some services might face issues when several requests are sent
                  simultaneously and will return the requests partially. It's recommended
                  to avoid using too many workers unless you are certain the web service
                  can handle it.
                * **verbose** (:class:`bool`, *optional*) -- If True, prints information about the requests and responses,
                  defaults to False.
                * **disable_retry** (:class:`bool`, *optional*) -- If ``True`` in case there are any failed queries, no retrying attempts
                  is done and object IDs of the failed requests is saved to a text file
                  which its path can be accessed via ``self.failed_path``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.


.. py:class:: WFSBase(__pydantic_self__, **data)



   Base class for WFS service.

   :Parameters: * **url** (:class:`str`) -- The base url for the WFS service, for examples:
                  https://hazards.fema.gov/nfhl/services/public/NFHL/MapServer/WFSServer
                * **layer** (:class:`str`) -- The layer from the service to be downloaded, defaults to None which throws
                  an error and includes all the available layers offered by the service.
                * **outformat** (:class:`str`) --

                  The data format to request for data from the service, defaults to None which
                   throws an error and includes all the available format offered by the service.
                * **version** (:class:`str`, *optional*) -- The WFS service version which should be either ``1.0.0``, ``1.1.0``, or
                  ``2.0.0``. Defaults to ``2.0.0``.
                * **crs** (:class:`str`, *optional*) -- The spatial reference system to be used for requesting the data, defaults to
                  ``epsg:4326``.
                * **read_method** (:class:`str`, *optional*) -- Method for reading the retrieved data, defaults to ``json``. Valid options are
                  ``json``, ``binary``, and ``text``.
                * **max_nrecords** (:class:`int`, *optional*) -- The maximum number of records in a single request to be retrieved from the service,
                  defaults to 1000. If the number of requested records is greater than this value,
                  the query will be split into multiple requests.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to ``False``.

   .. py:method:: get_validnames(self)

      Get valid column names for a layer.


   .. py:method:: validate_wfs(self)

      Validate input arguments with the WFS service.



.. py:class:: WMSBase(__pydantic_self__, **data)



   Base class for accessing a WMS service.

   :Parameters: * **url** (:class:`str`) -- The base url for the WMS service e.g., https://www.mrlc.gov/geoserver/mrlc_download/wms
                * **layers** (:class:`str` or :class:`list`) -- A layer or a list of layers from the service to be downloaded. You can pass an empty
                  string to get a list of available layers.
                * **outformat** (:class:`str`) -- The data format to request for data from the service. You can pass an empty
                  string to get a list of available output formats.
                * **version** (:class:`str`, *optional*) -- The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.
                * **crs** (:class:`str`, *optional*) -- The spatial reference system to be used for requesting the data, defaults to
                  ``epsg:4326``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: get_validlayers(self)

      Get the layers supported by the WMS service.


   .. py:method:: validate_wms(self)

      Validate input arguments with the WMS service.



.. py:function:: validate_version(val, valid_versions)

   Validate version from a list of valid versions.

   :Parameters: * **val** (:class:`str`) -- Input version value.
                * **valid_versions** (:class:`list` of :class:`str`) -- List of valid versions.

   :returns: :class:`str` -- Validated version value.


:orphan:

:py:mod:`pygeohydro`
====================

.. py:module:: pygeohydro


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   helpers/index.rst
   plot/index.rst
   pygeohydro/index.rst
   waterdata/index.rst


Package Contents
----------------













:py:mod:`pygeohydro.pygeohydro`
===============================

.. py:module:: pygeohydro.pygeohydro

.. autoapi-nested-parse::

   Accessing data from the supported databases through their APIs.



Module Contents
---------------

.. py:class:: NID(expire_after = EXPIRE, disable_caching = False)

   Retrieve data from the National Inventory of Dams web service.

   :Parameters: * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: get_byfilter(self, query_list)

      Query dams by filters from the National Inventory of Dams web service.

      :Parameters: **query_list** (:class:`list` of :class:`dict`) -- List of dictionary of query parameters. For an exhaustive list of the parameters,
                   use the advanced fields dataframe that can be accessed via ``NID().fields_meta``.
                   Some filter require min/max values such as ``damHeight`` and ``drainageArea``.
                   For such filters, the min/max values should be passed like so:
                   ``{filter_key: ["[min1 max1]", "[min2 max2]"]}``.

      :returns: :class:`geopandas.GeoDataFrame` -- Query results.

      .. rubric:: Examples

      >>> from pygeohydro import NID
      >>> nid = NID()
      >>> query_list = [
      ...    {"huc6": ["160502", "100500"], "drainageArea": ["[200 500]"]},
      ...    {"nidId": ["CA01222"]},
      ... ]
      >>> dam_dfs = nid.get_byfilter(query_list)
      >>> print(dam_dfs[0].name[0])
      Stillwater Point Dam


   .. py:method:: get_bygeom(self, geometry, geo_crs)

      Retrieve NID data within a geometry.

      :Parameters: * **geometry** (:class:`Polygon`, :class:`MultiPolygon`, or :class:`tuple` of :class:`length 4`) -- Geometry or bounding box (west, south, east, north) for extracting the data.
                   * **geo_crs** (:class:`list` of :class:`str`) -- The CRS of the input geometry, defaults to epsg:4326.

      :returns: :class:`geopandas.GeoDataFrame` -- GeoDataFrame of NID data

      .. rubric:: Examples

      >>> from pygeohydro import NID
      >>> nid = NID()
      >>> dams = nid.get_bygeom((-69.77, 45.07, -69.31, 45.45), "epsg:4326")
      >>> print(dams.name.iloc[0])
      Little Moose


   .. py:method:: get_suggestions(self, text, context_key = '')

      Get suggestions from the National Inventory of Dams web service.

      .. rubric:: Notes

      This function is useful for exploring and/or narrowing down the filter fields
      that are needed to query the dams using ``get_byfilter``.

      :Parameters: * **text** (:class:`str`) -- Text to query for suggestions.
                   * **context_key** (:class:`str`, *optional*) -- Suggestion context, defaults to empty string, i.e., all context keys.
                     For a list of valid context keys, see ``NID().fields_meta``.

      :returns: :class:`tuple` of :class:`pandas.DataFrame` -- The suggestions for the requested text as two DataFrames:
                First, is suggestions found in the dams properties and
                second, those found in the query fields such as states, huc6, etc.

      .. rubric:: Examples

      >>> from pygeohydro import NID
      >>> nid = NID()
      >>> dams, contexts = nid.get_suggestions("texas", "huc2")
      >>> print(contexts.loc["HUC2", "value"])
      12


   .. py:method:: inventory_byid(self, dam_ids)

      Get extra attributes for dams based on their dam ID.

      .. rubric:: Notes

      This function is meant to be used for getting extra attributes for dams.
      For example, first you need to use either ``get_bygeom`` or ``get_byfilter``
      to get basic attributes of the target dams. Then you can use this function
      to get extra attributes using the ``id`` column of the ``GeoDataFrame``
      that ``get_bygeom`` or ``get_byfilter`` returns.

      :Parameters: **dam_ids** (:class:`list` of :class:`int` or :class:`str`) -- List of the target dam IDs (digists only). Note that the dam IDs are not the
                   same as the NID IDs.

      :returns: :class:`pandas.DataFrame` -- Dams with extra attributes in addition to the standard NID fields
                that other ``NID`` methods return.

      .. rubric:: Examples

      >>> from pygeohydro import NID
      >>> nid = NID()
      >>> dams = nid.inventory_byid([514871, 459170, 514868, 463501, 463498])
      >>> print(dams.damHeight.max())
      120.0



.. py:function:: cover_statistics(ds)

   Percentages of the categorical NLCD cover data.

   :Parameters: **ds** (:class:`xarray.DataArray`) -- Cover DataArray from a LULC Dataset from the ``nlcd`` function.

   :returns: :class:`dict` -- Statistics of NLCD cover data


.. py:function:: nlcd(geometry, resolution, years = None, region = 'L48', geo_crs = DEF_CRS, crs = DEF_CRS)

   Get data from NLCD database (2019).

   .. deprecated:: 0.11.5
       Use :func:`nlcd_bygeom` or :func:`nlcd_bycoords`  instead.

   :Parameters: * **geometry** (:class:`Polygon`, :class:`MultiPolygon`, or :class:`tuple` of :class:`length 4`) -- The geometry or bounding box (west, south, east, north) for extracting the data.
                * **resolution** (:class:`float`) -- The data resolution in meters. The width and height of the output are computed in pixel
                  based on the geometry bounds and the given resolution.
                * **years** (:class:`dict`, *optional*) -- The years for NLCD layers as a dictionary, defaults to
                  ``{'impervious': [2019], 'cover': [2019], 'canopy': [2019], "descriptor": [2019]}``.
                  Layers that are not in years are ignored, e.g., ``{'cover': [2016, 2019]}`` returns
                  land cover data for 2016 and 2019.
                * **region** (:class:`str`, *optional*) -- Region in the US, defaults to ``L48``. Valid values are ``L48`` (for CONUS),
                  ``HI`` (for Hawaii), ``AK`` (for Alaska), and ``PR`` (for Puerto Rico).
                  Both lower and upper cases are acceptable.
                * **geo_crs** (:class:`str`, *optional*) -- The CRS of the input geometry, defaults to epsg:4326.
                * **crs** (:class:`str`, *optional*) -- The spatial reference system to be used for requesting the data, defaults to
                  epsg:4326.

   :returns: :class:`xarray.Dataset` -- NLCD within a geometry


.. py:function:: nlcd_bycoords(coords, years = None, region = 'L48', expire_after = EXPIRE, disable_caching = False)

   Get data from NLCD database (2019).

   :Parameters: * **coords** (:class:`list` of :class:`tuple`) -- List of coordinates in the form of (longitude, latitude).
                * **years** (:class:`dict`, *optional*) -- The years for NLCD layers as a dictionary, defaults to
                  ``{'impervious': [2019], 'cover': [2019], 'canopy': [2019], "descriptor": [2019]}``.
                  Layers that are not in years are ignored, e.g., ``{'cover': [2016, 2019]}`` returns
                  land cover data for 2016 and 2019.
                * **region** (:class:`str`, *optional*) -- Region in the US, defaults to ``L48``. Valid values are ``L48`` (for CONUS),
                  ``HI`` (for Hawaii), ``AK`` (for Alaska), and ``PR`` (for Puerto Rico).
                  Both lower and upper cases are acceptable.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`geopandas.GeoDataFrame` -- A GeoDataFrame with the NLCD data and the coordinates.


.. py:function:: nlcd_bygeom(geometry, resolution, years = None, region = 'L48', crs = DEF_CRS, expire_after = EXPIRE, disable_caching = False)

   Get data from NLCD database (2019).

   :Parameters: * **geometry** (:class:`geopandas.GeoDataFrame` or :class:`geopandas.GeoSeries`) -- A GeoDataFrame or GeoSeries with the geometry to query. The indices are used
                  as keys in the output dictionary.
                * **resolution** (:class:`float`) -- The data resolution in meters. The width and height of the output are computed in pixel
                  based on the geometry bounds and the given resolution.
                * **years** (:class:`dict`, *optional*) -- The years for NLCD layers as a dictionary, defaults to
                  ``{'impervious': [2019], 'cover': [2019], 'canopy': [2019], "descriptor": [2019]}``.
                  Layers that are not in years are ignored, e.g., ``{'cover': [2016, 2019]}`` returns
                  land cover data for 2016 and 2019.
                * **region** (:class:`str`, *optional*) -- Region in the US, defaults to ``L48``. Valid values are ``L48`` (for CONUS),
                  ``HI`` (for Hawaii), ``AK`` (for Alaska), and ``PR`` (for Puerto Rico).
                  Both lower and upper cases are acceptable.
                * **crs** (:class:`str`, *optional*) -- The spatial reference system to be used for requesting the data, defaults to
                  epsg:4326.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`dict` of :class:`xarray.Dataset` or :class:`xarray.Dataset` -- A single or a ``dict`` of NLCD datasets. If dict, the keys are indices
             of the input ``GeoDataFrame``.


.. py:function:: ssebopeta_bycoords(coords, dates, crs = DEF_CRS)

   Daily actual ET for a dataframe of coords from SSEBop database in mm/day.

   :Parameters: * **coords** (:class:`pandas.DataFrame`) -- A dataframe with ``id``, ``x``, ``y`` columns.
                * **dates** (:class:`tuple` or :class:`list`, *optional*) -- Start and end dates as a tuple (start, end) or a list of years [2001, 2010, ...].
                * **crs** (:class:`str`, *optional*) -- The CRS of the input coordinates, defaults to epsg:4326.

   :returns: :class:`xarray.Dataset` -- Daily actual ET in mm/day as a dataset with ``time`` and ``location_id`` dimensions.
             The ``location_id`` dimension is the same as the ``id`` column in the input dataframe.


.. py:function:: ssebopeta_bygeom(geometry, dates, geo_crs = DEF_CRS)

   Get daily actual ET for a region from SSEBop database.

   .. rubric:: Notes

   Since there's still no web service available for subsetting SSEBop, the data first
   needs to be downloaded for the requested period then it is masked by the
   region of interest locally. Therefore, it's not as fast as other functions and
   the bottleneck could be the download speed.

   :Parameters: * **geometry** (:class:`shapely.geometry.Polygon` or :class:`tuple`) -- The geometry for downloading clipping the data. For a tuple bbox,
                  the order should be (west, south, east, north).
                * **dates** (:class:`tuple` or :class:`list`, *optional*) -- Start and end dates as a tuple (start, end) or a list of years [2001, 2010, ...].
                * **geo_crs** (:class:`str`, *optional*) -- The CRS of the input geometry, defaults to epsg:4326.

   :returns: :class:`xarray.DataArray` -- Daily actual ET within a geometry in mm/day at 1 km resolution


.. py:function:: ssebopeta_byloc(coords, dates)

   Daily actual ET for a location from SSEBop database in mm/day.

   .. deprecated:: 0.11.5
       Use :func:`ssebopeta_bycoords` instead. For now, this function calls
       :func:`ssebopeta_bycoords` but retains the same functionality, i.e.,
       returns a dataframe and accepts only a single coordinate. Whereas the
       new function returns a ``xarray.Dataset`` and accepts a dataframe
       containing coordinates.

   :Parameters: * **coords** (:class:`tuple`) -- Longitude and latitude of a single location as a tuple (lon, lat)
                * **dates** (:class:`tuple` or :class:`list`, *optional*) -- Start and end dates as a tuple (start, end) or a list of years [2001, 2010, ...].

   :returns: :class:`pandas.Series` -- Daily actual ET for a location


:py:mod:`pygeohydro.plot`
=========================

.. py:module:: pygeohydro.plot

.. autoapi-nested-parse::

   Plot hydrological signatures.

   Plots includes  daily, monthly and annual hydrograph as well as regime
   curve (monthly mean) and flow duration curve.



Module Contents
---------------

.. py:class:: PlotDataType



   Data structure for plotting hydrologic signatures.


.. py:function:: cover_legends()

   Colormap (cmap) and their respective values (norm) for land cover data legends.


.. py:function:: descriptor_legends()

   Colormap (cmap) and their respective values (norm) for land cover data legends.


.. py:function:: exceedance(daily)

   Compute Flow duration (rank, sorted obs).


.. py:function:: prepare_plot_data(daily)

   Generae a structured data for plotting hydrologic signatures.

   :Parameters: **daily** (:class:`pandas.Series` or :class:`pandas.DataFrame`) -- The data to be processed

   :returns: :class:`PlotDataType` -- Containing ``daily, ``monthly``, ``annual``, ``mean_monthly``, ``ranked`` fields.


.. py:function:: signatures(daily, precipitation = None, title = None, title_ypos = 1.02, figsize = (14, 13), threshold = 0.001, output = None)

   Plot hydrological signatures with w/ or w/o precipitation.

   Plots includes daily, monthly and annual hydrograph as well as
   regime curve (mean monthly) and flow duration curve. The input
   discharges are converted from cms to mm/day based on the watershed
   area, if provided.

   :Parameters: * **daily** (:class:`pd.DataFrame` or :class:`pd.Series`) -- The streamflows in mm/day. The column names are used as labels
                  on the plot and the column values should be daily streamflow.
                * **precipitation** (:class:`pd.Series`, *optional*) -- Daily precipitation time series in mm/day. If given, the data is
                  plotted on the second x-axis at the top.
                * **title** (:class:`str`, *optional*) -- The plot supertitle.
                * **title_ypos** (:class:`float`) -- The vertical position of the plot title, default to 1.02
                * **figsize** (:class:`tuple`, *optional*) -- Width and height of the plot in inches, defaults to (14, 13) inches.
                * **threshold** (:class:`float`, *optional*) -- The threshold for cutting off the discharge for the flow duration
                  curve to deal with log 0 issue, defaults to :math:`1^{-3}` mm/day.
                * **output** (:class:`str`, *optional*) -- Path to save the plot as png, defaults to ``None`` which means
                  the plot is not saved to a file.


:py:mod:`pygeohydro.helpers`
============================

.. py:module:: pygeohydro.helpers

.. autoapi-nested-parse::

   Some helper function for PyGeoHydro.



Module Contents
---------------

.. py:function:: nlcd_helper()

   Get legends and properties of the NLCD cover dataset.

   .. rubric:: Notes

   The following references have been used:
       - https://github.com/jzmiller1/nlcd
       - https://www.mrlc.gov/data-services-page
       - https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend

   :returns: :class:`dict` -- Years where data is available and cover classes and categories, and roughness estimations.


.. py:function:: nwis_errors()

   Get error code lookup table for USGS sites that have daily values.


:py:mod:`pygeohydro.waterdata`
==============================

.. py:module:: pygeohydro.waterdata

.. autoapi-nested-parse::

   Accessing data from the supported databases through their APIs.



Module Contents
---------------

.. py:class:: NWIS(expire_after = EXPIRE, disable_caching = False)

   Access NWIS web service.

   :Parameters: * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: get_info(self, queries, expanded = False)

      Send multiple queries to USGS Site Web Service.

      :Parameters: * **queries** (:class:`dict` or :class:`list` of :class:`dict`) -- A single or a list of valid queries.
                   * **expanded** (:class:`bool`, *optional*) -- Whether to get expanded sit information for example drainage area,
                     default to False.

      :returns: :class:`geopandas.GeoDataFrame` -- A correctly typed ``GeoDataFrame`` containing site(s) information.


   .. py:method:: get_parameter_codes(self, keyword)

      Search for parameter codes by name or number.

      .. rubric:: Notes

      NWIS guideline for keywords is as follows:

          By default an exact search is made. To make a partial search the term
          should be prefixed and suffixed with a % sign. The % sign matches zero
          or more characters at the location. For example, to find all with "discharge"
          enter %discharge% in the field. % will match any number of characters
          (including zero characters) at the location.

      :Parameters: **keyword** (:class:`str`) -- Keyword to search for parameters by name of number.

      :returns: :class:`pandas.DataFrame` -- Matched parameter codes as a dataframe with their description.

      .. rubric:: Examples

      >>> from pygeohydro import NWIS
      >>> nwis = NWIS()
      >>> codes = nwis.get_parameter_codes("%discharge%")
      >>> codes.loc[codes.parameter_cd == "00060", "parm_nm"][0]
      'Discharge, cubic feet per second'


   .. py:method:: get_streamflow(self, station_ids, dates, freq = 'dv', mmd = False, to_xarray = False)

      Get mean daily streamflow observations from USGS.

      :Parameters: * **station_ids** (:class:`str`, :class:`list`) -- The gage ID(s)  of the USGS station.
                   * **dates** (:class:`tuple`) -- Start and end dates as a tuple (start, end).
                   * **freq** (:class:`str`, *optional*) -- The frequency of the streamflow data, defaults to ``dv`` (daily values).
                     Valid frequencies are ``dv`` (daily values), ``iv`` (instantaneous values).
                     Note that for ``iv`` the time zone for the input dates is assumed to be UTC.
                   * **mmd** (:class:`bool`, *optional*) -- Convert cms to mm/day based on the contributing drainage area of the stations.
                     Defaults to False.
                   * **to_xarray** (:class:`bool`, *optional*) -- Whether to return a xarray.Dataset. Defaults to False.

      :returns: :class:`pandas.DataFrame` or :class:`xarray.Dataset` -- Streamflow data observations in cubic meter per second (cms). The stations that
                don't provide the requested discharge data in the target period will be dropped.
                Note that when frequency is set to ``iv`` the time zone is converted to UTC.


   .. py:method:: retrieve_rdb(self, url, payloads)

      Retrieve and process requests with RDB format.

      :Parameters: * **url** (:class:`str`) -- Name of USGS REST service, valid values are ``site``, ``dv``, ``iv``,
                     ``gwlevels``, and ``stat``. Please consult USGS documentation
                     `here <https://waterservices.usgs.gov/rest>`__ for more information.
                   * **payloads** (:class:`list` of :class:`dict`) -- List of target payloads.

      :returns: :class:`pandas.DataFrame` -- Requested features as a pandas's DataFrame.



.. py:class:: WaterQuality(expire_after = EXPIRE, disable_caching = False)

   Water Quality Web Service https://www.waterqualitydata.us.

   .. rubric:: Notes

   This class has a number of convenience methods to retrieve data from the
   Water Quality Data. Since there are many parameter combinations that can be
   used to retrieve data, a general method is also provided to retrieve data from
   any of the valid endpoints. You can use ``get_json`` to retrieve stations info
   as a ``geopandas.GeoDataFrame`` or ``get_csv`` to retrieve stations data as a
   ``pandas.DataFrame``. You can construct a dictionary of the parameters and pass
   it to one of these functions. For more information on the parameters, please
   consult the
   `Water Quality Data documentation <https://www.waterqualitydata.us/webservices_documentation>`__.

   :Parameters: * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   .. py:method:: data_bystation(self, station_ids, wq_kwds)

      Retrieve data for a single station.

      :Parameters: * **station_ids** (:class:`str` or :class:`list` of :class:`str`) -- Station ID(s). The IDs should have the format "Agency code-Station ID".
                   * **wq_kwds** (:class:`dict`, *optional*) -- Water Quality Web Service keyword arguments. Default to None.

      :returns: :class:`pandas.DataFrame` -- DataFrame of data for the stations.


   .. py:method:: get_csv(self, endpoint, kwds, request_method = 'GET')

      Get the CSV response from the Water Quality Web Service.

      :Parameters: * **endpoint** (:class:`str`) -- Endpoint of the Water Quality Web Service.
                   * **kwds** (:class:`dict`) -- Water Quality Web Service keyword arguments.
                   * **request_method** (:class:`str`, *optional*) -- HTTP request method. Default to GET.

      :returns: :class:`pandas.DataFrame` -- The web service response as a DataFrame.


   .. py:method:: get_json(self, endpoint, kwds, request_method = 'GET')

      Get the JSON response from the Water Quality Web Service.

      :Parameters: * **endpoint** (:class:`str`) -- Endpoint of the Water Quality Web Service.
                   * **kwds** (:class:`dict`) -- Water Quality Web Service keyword arguments.
                   * **request_method** (:class:`str`, *optional*) -- HTTP request method. Default to GET.

      :returns: :class:`geopandas.GeoDataFrame` -- The web service response as a GeoDataFrame.


   .. py:method:: get_param_table(self)

      Get the parameter table from the USGS Water Quality Web Service.


   .. py:method:: lookup_domain_values(self, endpoint)

      Get the domain values for the target endpoint.


   .. py:method:: station_bybbox(self, bbox, wq_kwds)

      Retrieve station info within bounding box.

      :Parameters: * **bbox** (:class:`tuple` of :class:`float`) -- Bounding box coordinates (west, south, east, north) in epsg:4326.
                   * **wq_kwds** (:class:`dict`, *optional*) -- Water Quality Web Service keyword arguments. Default to None.

      :returns: :class:`geopandas.GeoDataFrame` -- GeoDataFrame of station info within the bounding box.


   .. py:method:: station_bydistance(self, lon, lat, radius, wq_kwds)

      Retrieve station within a radius (decimal miles) of a point.

      :Parameters: * **lon** (:class:`float`) -- Longitude of point.
                   * **lat** (:class:`float`) -- Latitude of point.
                   * **radius** (:class:`float`) -- Radius (decimal miles) of search.
                   * **wq_kwds** (:class:`dict`, *optional*) -- Water Quality Web Service keyword arguments. Default to None.

      :returns: :class:`geopandas.GeoDataFrame` -- GeoDataFrame of station info within the radius of the point.



.. py:function:: interactive_map(bbox, crs = DEF_CRS, nwis_kwds = None, expire_after = EXPIRE, disable_caching = False)

   Generate an interactive map including all USGS stations within a bounding box.

   :Parameters: * **bbox** (:class:`tuple`) -- List of corners in this order (west, south, east, north)
                * **crs** (:class:`str`, *optional*) -- CRS of the input bounding box, defaults to EPSG:4326.
                * **nwis_kwds** (:class:`dict`, *optional*) -- Optional keywords to include in the NWIS request as a dictionary like so:
                  ``{"hasDataTypeCd": "dv,iv", "outputDataTypeCd": "dv,iv", "parameterCd": "06000"}``.
                  Default to None.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`folium.Map` -- Interactive map within a bounding box.

   .. rubric:: Examples

   >>> import pygeohydro as gh
   >>> nwis_kwds = {"hasDataTypeCd": "dv,iv", "outputDataTypeCd": "dv,iv"}
   >>> m = gh.interactive_map((-69.77, 45.07, -69.31, 45.45), nwis_kwds=nwis_kwds)
   >>> n_stations = len(m.to_dict()["children"]) - 1
   >>> n_stations
   10


:py:mod:`py3dep`
================

.. py:module:: py3dep

.. autoapi-nested-parse::

   Top-level package for Py3DEP.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   py3dep/index.rst
   utils/index.rst


Package Contents
----------------








:py:mod:`py3dep.utils`
======================

.. py:module:: py3dep.utils

.. autoapi-nested-parse::

   Utilities for Py3DEP.



Module Contents
---------------

.. py:function:: deg2mpm(slope)

   Convert slope from degrees to meter/meter.

   :Parameters: **slope** (:class:`xarray.DataArray`) -- Slope in degrees.

   :returns: :class:`xarray.DataArray` -- Slope in meter/meter. The name is set to ``slope`` and the ``units`` attribute
             is set to ``m/m``.


.. py:function:: fill_depressions(dem_da)

   Fill depressions and adjust flat areas in a DEM using `RichDEM <https://richdem.readthedocs.io>`__.

   :Parameters: **dem** (:class:`xarray.DataArray` or :class:`numpy.ndarray`) -- Digital Elevation Model.

   :returns: :class:`xarray.DataArray` -- Conditioned DEM after applying
             `depression filling <https://richdem.readthedocs.io/en/latest/depression_filling.html>`__
             and
             `flat area resolution <https://richdem.readthedocs.io/en/latest/flat_resolution.html>`__
             operations.


:py:mod:`py3dep.py3dep`
=======================

.. py:module:: py3dep.py3dep

.. autoapi-nested-parse::

   Get data from 3DEP database.



Module Contents
---------------

.. py:function:: check_3dep_availability(bbox, crs = DEF_CRS)

   Query 3DEP's resolution availability within a bounding box.

   This function checks availability of 3DEP's at the following resolutions:
   1 m, 3 m, 5 m, 10 m, 30 m, and 60 m.

   :Parameters: * **bbox** (:class:`tuple`) -- Bounding box as tuple of ``(min_x, min_y, max_x, max_y)``.
                * **crs** (:class:`str` or :class:`pyproj.CRS`, *optional*) -- Spatial reference (CRS) of bbox, defaults to ``EPSG:4326``.

   :returns: :class:`dict` -- True if bbox intersects 3DEP elevation for each available resolution.
             Keys are the supported resolutions and values are their availability.

   .. rubric:: Examples

   >>> import py3dep
   >>> bbox = (-69.77, 45.07, -69.31, 45.45)
   >>> py3dep.check_3dep_availability(bbox)
   {'1m': True, '3m': False, '5m': False, '10m': True, '30m': True, '60m': False}


.. py:function:: elevation_bycoords(coords, crs = DEF_CRS, source = 'tep', expire_after = EXPIRE, disable_caching = False)

   Get elevation for a list of coordinates.

   :Parameters: * **coords** (:class:`list` of :class:`tuple`) -- Coordinates of target location as list of tuples ``[(x, y), ...]``.
                * **crs** (:class:`str` or :class:`pyproj.CRS`, *optional*) -- Spatial reference (CRS) of coords, defaults to ``EPSG:4326``.
                * **source** (:class:`str`, *optional*) -- Data source to be used, default to ``airmap``. Supported sources are
                  ``airmap`` (30 m resolution), ``tnm`` (using The National Map's Bulk Point
                  Query Service with 10 m resolution) and ``tep`` (using 3DEP's WMS service
                  at 10 m resolution). The ``tnm`` and ``tep`` sources are more accurate since they
                  use the 1/3 arc-second DEM layer from 3DEP service but it is limited to the US.
                  They both tend to be slower than the Airmap service. Note that ``tnm`` is bit unstable.
                  It's recommended to use ``tep`` unless 10-m resolution accuracy is not necessary which
                  in that case ``airmap`` is more appropriate.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`list` of :class:`float` -- Elevation in meter.


.. py:function:: elevation_bygrid(xcoords, ycoords, crs, resolution, depression_filling = False, expire_after = EXPIRE, disable_caching = False)

   Get elevation from DEM data for a grid.

   This function is intended for getting elevations for a gridded dataset.

   :Parameters: * **xcoords** (:class:`list`) -- List of x-coordinates of a grid.
                * **ycoords** (:class:`list`) -- List of y-coordinates of a grid.
                * **crs** (:class:`str` or :class:`pyproj.CRS`) -- The spatial reference system of the input grid, defaults to ``EPSG:4326``.
                * **resolution** (:class:`float`) -- The accuracy of the output, defaults to 10 m which is the highest
                  available resolution that covers CONUS. Note that higher resolution
                  increases computation time so chose this value with caution.
                * **depression_filling** (:class:`bool`, *optional*) -- Fill depressions before sampling using
                  `RichDEM <https://richdem.readthedocs.io/en/latest/>`__ package, defaults to False.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`xarray.DataArray` -- Elevations of the input coordinates as a ``xarray.DataArray``.


.. py:function:: get_map(layers, geometry, resolution, geo_crs = DEF_CRS, crs = DEF_CRS, expire_after = EXPIRE, disable_caching = False)

   Access to `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__ service.

   The 3DEP service has multi-resolution sources, so depending on the user
   provided resolution the data is resampled on server-side based
   on all the available data sources. The following layers are available:

   - ``DEM``
   - ``Hillshade Gray``
   - ``Aspect Degrees``
   - ``Aspect Map``
   - ``GreyHillshade_elevationFill``
   - ``Hillshade Multidirectional``
   - ``Slope Map``
   - ``Slope Degrees``
   - ``Hillshade Elevation Tinted``
   - ``Height Ellipsoidal``
   - ``Contour 25``
   - ``Contour Smoothed 25``

   :Parameters: * **layers** (:class:`str` or :class:`list` of :class:`str`) -- A valid 3DEP layer or a list of them.
                * **geometry** (:class:`Polygon`, :class:`MultiPolygon`, or :class:`tuple`) -- A shapely Polygon or a bounding box of the form ``(west, south, east, north)``.
                * **resolution** (:class:`float`) -- The target resolution in meters. The width and height of the output are computed in
                  pixels based on the geometry bounds and the given resolution.
                * **geo_crs** (:class:`str`, *optional*) -- The spatial reference system of the input geometry, defaults to ``EPSG:4326``.
                * **crs** (:class:`str`, *optional*) -- The spatial reference system to be used for requesting the data, defaults to
                  ``EPSG:4326``. Valid values are ``EPSG:4326``, ``EPSG:3576``, ``EPSG:3571``,
                  ``EPSG:3575``, ``EPSG:3857``, ``EPSG:3572``, ``CRS:84``, ``EPSG:3573``,
                  and ``EPSG:3574``.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`xarray.DataArray` or :class:`xarray.Dataset` -- The requested topographic data as an ``xarray.DataArray`` or ``xarray.Dataset``.


:py:mod:`async_retriever`
=========================

.. py:module:: async_retriever

.. autoapi-nested-parse::

   Top-level package.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   async_retriever/index.rst
   utils/index.rst


Package Contents
----------------






:py:mod:`async_retriever.utils`
===============================

.. py:module:: async_retriever.utils

.. autoapi-nested-parse::

   Core async functions.



Module Contents
---------------

.. py:class:: BaseRetriever(urls, read, request_kwds = None, request_method = 'GET', cache_name = None, family = 'both')

   Base class for async retriever.

   .. py:method:: generate_requests(urls, request_kwds)
      :staticmethod:

      Generate urls and keywords.



.. py:function:: create_cachefile(db_name = None)

   Create a cache folder in the current working directory.


.. py:function:: delete_url(url, method = 'GET', cache_name = None, **kwargs)
   :async:

   Delete cached response associated with ``url``.


.. py:function:: get_event_loop()

   Create an event loop.


.. py:function:: retriever(uid, url, s_kwds, session, read_type, r_kwds)
   :async:

   Create an async request and return the response as binary.

   :Parameters: * **uid** (:class:`int`) -- ID of the URL for sorting after returning the results
                * **url** (:class:`str`) -- URL to be retrieved
                * **s_kwds** (:class:`dict`) -- Arguments to be passed to requests
                * **session** (:class:`ClientSession`) -- A ClientSession for sending the request
                * **read_type** (:class:`str`) -- Return response as text, bytes, or json.
                * **r_kwds** (:class:`dict`) -- Keywords to pass to the response read function.
                  It is ``{"content_type": None}`` if ``read`` is ``json``
                  else an empty ``dict``.

   :returns: :class:`bytes` -- The retrieved response as binary.


:py:mod:`async_retriever.async_retriever`
=========================================

.. py:module:: async_retriever.async_retriever

.. autoapi-nested-parse::

   Core async functions.



Module Contents
---------------

.. py:function:: delete_url_cache(url, request_method = 'GET', cache_name = None, **kwargs)

   Delete cached response associated with ``url``, along with its history (if applicable).

   :Parameters: * **url** (:class:`str`) -- URL to be deleted from the cache
                * **request_method** (:class:`str`, *optional*) -- HTTP request method to be deleted from the cache, defaults to ``GET``.
                * **cache_name** (:class:`str`, *optional*) -- Path to a file for caching the session, defaults to
                  ``./cache/aiohttp_cache.sqlite``.
                * **kwargs** (:class:`dict`, *optional*) -- Keywords to pass to the ``cache.delete_url()``.


.. py:function:: retrieve(urls, read, request_kwds = None, request_method = 'GET', max_workers = 8, cache_name = None, family = 'both', timeout = 5.0, expire_after = EXPIRE, ssl = None, disable = False)

   Send async requests.

   :Parameters: * **urls** (:class:`list` of :class:`str`) -- List of URLs.
                * **read** (:class:`str`) -- Method for returning the request; ``binary``, ``json``, and ``text``.
                * **request_kwds** (:class:`list` of :class:`dict`, *optional*) -- List of requests keywords corresponding to input URLs (1 on 1 mapping),
                  defaults to ``None``. For example, ``[{"params": {...}, "headers": {...}}, ...]``.
                * **request_method** (:class:`str`, *optional*) -- Request type; ``GET`` (``get``) or ``POST`` (``post``). Defaults to ``GET``.
                * **max_workers** (:class:`int`, *optional*) -- Maximum number of async processes, defaults to 8.
                * **cache_name** (:class:`str`, *optional*) -- Path to a file for caching the session, defaults to ``./cache/aiohttp_cache.sqlite``.
                * **family** (:class:`str`, *optional*) -- TCP socket family, defaults to both, i.e., IPv4 and IPv6. For IPv4
                  or IPv6 only pass ``ipv4`` or ``ipv6``, respectively.
                * **timeout** (:class:`float`, *optional*) -- Timeout for the request, defaults to 5.0.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **ssl** (:class:`bool` or :class:`SSLContext`, *optional*) -- SSLContext to use for the connection, defaults to None. Set to False to disable
                  SSL cetification verification.
                * **disable** (:class:`bool`, *optional*) -- If ``True`` temporarily disable caching requests and get new responses
                  from the server, defaults to False.

   :returns: :class:`list` -- List of responses in the order of input URLs.

   .. rubric:: Examples

   >>> import async_retriever as ar
   >>> stations = ["01646500", "08072300", "11073495"]
   >>> url = "https://waterservices.usgs.gov/nwis/site"
   >>> urls, kwds = zip(
   ...     *[
   ...         (url, {"params": {"format": "rdb", "sites": s, "siteStatus": "all"}})
   ...         for s in stations
   ...     ]
   ... )
   >>> resp = ar.retrieve(urls, "text", request_kwds=kwds)
   >>> resp[0].split('\n')[-2].split('\t')[1]
   '01646500'


.. py:function:: retrieve_binary(urls, request_kwds = None, request_method = 'GET', max_workers = 8, cache_name = None, family = 'both', timeout = 5.0, expire_after = EXPIRE, ssl = None, disable = False)

   Send async requests and get the response as ``bytes``.

   :Parameters: * **urls** (:class:`list` of :class:`str`) -- List of URLs.
                * **request_kwds** (:class:`list` of :class:`dict`, *optional*) -- List of requests keywords corresponding to input URLs (1 on 1 mapping),
                  defaults to ``None``. For example, ``[{"params": {...}, "headers": {...}}, ...]``.
                * **request_method** (:class:`str`, *optional*) -- Request type; ``GET`` (``get``) or ``POST`` (``post``). Defaults to ``GET``.
                * **max_workers** (:class:`int`, *optional*) -- Maximum number of async processes, defaults to 8.
                * **cache_name** (:class:`str`, *optional*) -- Path to a file for caching the session, defaults to ``./cache/aiohttp_cache.sqlite``.
                * **family** (:class:`str`, *optional*) -- TCP socket family, defaults to both, i.e., IPv4 and IPv6. For IPv4
                  or IPv6 only pass ``ipv4`` or ``ipv6``, respectively.
                * **timeout** (:class:`float`, *optional*) -- Timeout for the request, defaults to 5.0.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **ssl** (:class:`bool` or :class:`SSLContext`, *optional*) -- SSLContext to use for the connection, defaults to None. Set to False to disable
                  SSL cetification verification.
                * **disable** (:class:`bool`, *optional*) -- If ``True`` temporarily disable caching requests and get new responses
                  from the server, defaults to False.

   :returns: :class:`bytes` -- List of responses in the order of input URLs.


.. py:function:: retrieve_json(urls, request_kwds = None, request_method = 'GET', max_workers = 8, cache_name = None, family = 'both', timeout = 5.0, expire_after = EXPIRE, ssl = None, disable = False)

   Send async requests and get the response as ``json``.

   :Parameters: * **urls** (:class:`list` of :class:`str`) -- List of URLs.
                * **request_kwds** (:class:`list` of :class:`dict`, *optional*) -- List of requests keywords corresponding to input URLs (1 on 1 mapping),
                  defaults to ``None``. For example, ``[{"params": {...}, "headers": {...}}, ...]``.
                * **request_method** (:class:`str`, *optional*) -- Request type; ``GET`` (``get``) or ``POST`` (``post``). Defaults to ``GET``.
                * **max_workers** (:class:`int`, *optional*) -- Maximum number of async processes, defaults to 8.
                * **cache_name** (:class:`str`, *optional*) -- Path to a file for caching the session, defaults to ``./cache/aiohttp_cache.sqlite``.
                * **family** (:class:`str`, *optional*) -- TCP socket family, defaults to both, i.e., IPv4 and IPv6. For IPv4
                  or IPv6 only pass ``ipv4`` or ``ipv6``, respectively.
                * **timeout** (:class:`float`, *optional*) -- Timeout for the request, defaults to 5.0.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **ssl** (:class:`bool` or :class:`SSLContext`, *optional*) -- SSLContext to use for the connection, defaults to None. Set to False to disable
                  SSL cetification verification.
                * **disable** (:class:`bool`, *optional*) -- If ``True`` temporarily disable caching requests and get new responses
                  from the server, defaults to False.

   :returns: :class:`dict` -- List of responses in the order of input URLs.

   .. rubric:: Examples

   >>> import async_retriever as ar
   >>> urls = ["https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/position"]
   >>> kwds = [
   ...     {
   ...         "params": {
   ...             "f": "json",
   ...             "coords": "POINT(-68.325 45.0369)",
   ...         },
   ...     },
   ... ]
   >>> r = ar.retrieve_json(urls, kwds)
   >>> print(r[0]["features"][0]["properties"]["identifier"])
   2675320


.. py:function:: retrieve_text(urls, request_kwds = None, request_method = 'GET', max_workers = 8, cache_name = None, family = 'both', timeout = 5.0, expire_after = EXPIRE, ssl = None, disable = False)

   Send async requests and get the response as ``text``.

   :Parameters: * **urls** (:class:`list` of :class:`str`) -- List of URLs.
                * **request_kwds** (:class:`list` of :class:`dict`, *optional*) -- List of requests keywords corresponding to input URLs (1 on 1 mapping),
                  defaults to ``None``. For example, ``[{"params": {...}, "headers": {...}}, ...]``.
                * **request_method** (:class:`str`, *optional*) -- Request type; ``GET`` (``get``) or ``POST`` (``post``). Defaults to ``GET``.
                * **max_workers** (:class:`int`, *optional*) -- Maximum number of async processes, defaults to 8.
                * **cache_name** (:class:`str`, *optional*) -- Path to a file for caching the session, defaults to ``./cache/aiohttp_cache.sqlite``.
                * **family** (:class:`str`, *optional*) -- TCP socket family, defaults to both, i.e., IPv4 and IPv6. For IPv4
                  or IPv6 only pass ``ipv4`` or ``ipv6``, respectively.
                * **timeout** (:class:`float`, *optional*) -- Timeout for the request in seconds, defaults to 5.0.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **ssl** (:class:`bool` or :class:`SSLContext`, *optional*) -- SSLContext to use for the connection, defaults to None. Set to False to disable
                  SSL cetification verification.
                * **disable** (:class:`bool`, *optional*) -- If ``True`` temporarily disable caching requests and get new responses
                  from the server, defaults to False.

   :returns: :class:`list` -- List of responses in the order of input URLs.

   .. rubric:: Examples

   >>> import async_retriever as ar
   >>> stations = ["01646500", "08072300", "11073495"]
   >>> url = "https://waterservices.usgs.gov/nwis/site"
   >>> urls, kwds = zip(
   ...     *[
   ...         (url, {"params": {"format": "rdb", "sites": s, "siteStatus": "all"}})
   ...         for s in stations
   ...     ]
   ... )
   >>> resp = ar.retrieve_text(urls, kwds)
   >>> resp[0].split('\n')[-2].split('\t')[1]
   '01646500'


:orphan:

:mod:`noxfile`
==============

.. py:module:: noxfile


:py:mod:`pydaymet`
==================

.. py:module:: pydaymet

.. autoapi-nested-parse::

   Top-level package for PyDaymet.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   core/index.rst
   pet/index.rst
   pydaymet/index.rst


Package Contents
----------------





:py:mod:`pydaymet.pet`
======================

.. py:module:: pydaymet.pet

.. autoapi-nested-parse::

   Core class for the Daymet functions.



Module Contents
---------------

.. py:function:: potential_et(clm, coords = None, crs = DEF_CRS, method = 'hargreaves_samani', params = None)

   Compute Potential EvapoTranspiration for both gridded and a single location.

   :Parameters: * **clm** (:class:`pandas.DataFrame` or :class:`xarray.Dataset`) -- The dataset must include at least the following variables:

                  * Minimum temperature in degree celsius
                  * Maximum temperature in degree celsius
                  * Solar radiation in in W/m2
                  * Daylight duration in seconds

                  Optionally, relative humidity and wind speed at 2-m level will be used if available.

                  Table below shows the variable names that the function looks for in the input data.

                  ==================== ========
                  DataFrame            Dataset
                  ==================== ========
                  ``tmin (degrees C)`` ``tmin``
                  ``tmax (degrees C)`` ``tmax``
                  ``srad (W/m2)``      ``srad``
                  ``dayl (s)``         ``dayl``
                  ``rh (-)``           ``rh``
                  ``u2 (m/s)``         ``u2``
                  ==================== ========

                  If relative humidity and wind speed at 2-m level are not available,
                  actual vapour pressure is assumed to be saturation vapour pressure at daily minimum
                  temperature and 2-m wind speed is considered to be 2 m/s.
                * **coords** (:class:`tuple` of :class:`floats`, *optional*) -- Coordinates of the daymet data location as a tuple, (x, y). This is required when ``clm``
                  is a ``DataFrame``.
                * **crs** (:class:`str`, *optional*) -- The spatial reference of the input coordinate, defaults to ``epsg:4326``. This is only used
                  when ``clm`` is a ``DataFrame``.
                * **method** (:class:`str`, *optional*) -- Method for computing PET. Supported methods are
                  ``penman_monteith``, ``priestley_taylor``, ``hargreaves_samani``, and
                  None (don't compute PET). The ``penman_monteith`` method is based on
                  :footcite:t:`Allen_1998` assuming that soil heat flux density is zero.
                  The ``priestley_taylor`` method is based on
                  :footcite:t:`Priestley_1972` assuming that soil heat flux density is zero.
                  The ``hargreaves_samani`` method is based on :footcite:t:`Hargreaves_1982`.
                  Defaults to ``hargreaves_samani``.
                * **params** (:class:`dict`, *optional*) -- Model-specific parameters as a dictionary, defaults to ``None``.

   :returns: :class:`pandas.DataFrame` or :class:`xarray.Dataset` -- The input DataFrame/Dataset with an additional variable named ``pet (mm/day)`` for
             DataFrame and ``pet`` for Dataset.

   .. rubric:: References

   .. footbibliography::


:py:mod:`pydaymet.core`
=======================

.. py:module:: pydaymet.core

.. autoapi-nested-parse::

   Core class for the Daymet functions.



Module Contents
---------------

.. py:class:: Daymet(variables = None, pet = None, time_scale = 'daily', region = 'na')

   Base class for Daymet requests.

   :Parameters: * **variables** (:class:`str` or :class:`list` or :class:`tuple`, *optional*) -- List of variables to be downloaded. The acceptable variables are:
                  ``tmin``, ``tmax``, ``prcp``, ``srad``, ``vp``, ``swe``, ``dayl``
                  Descriptions can be found `here <https://daymet.ornl.gov/overview>`__.
                  Defaults to None i.e., all the variables are downloaded.
                * **pet** (:class:`str`, *optional*) -- Method for computing PET. Supported methods are
                  ``penman_monteith``, ``priestley_taylor``, ``hargreaves_samani``, and
                  None (don't compute PET). The ``penman_monteith`` method is based on
                  :footcite:t:`Allen_1998` assuming that soil heat flux density is zero.
                  The ``priestley_taylor`` method is based on
                  :footcite:t:`Priestley_1972` assuming that soil heat flux density is zero.
                  The ``hargreaves_samani`` method is based on :footcite:t:`Hargreaves_1982`.
                  Defaults to ``None``.
                * **time_scale** (:class:`str`, *optional*) -- Data time scale which can be daily, monthly (monthly summaries),
                  or annual (annual summaries). Defaults to daily.
                * **region** (:class:`str`, *optional*) -- Region in the US, defaults to na. Acceptable values are:

                  * na: Continental North America
                  * hi: Hawaii
                  * pr: Puerto Rico

   .. rubric:: References

   .. footbibliography::

   .. py:method:: check_dates(dates)
      :staticmethod:

      Check if input dates are in correct format and valid.


   .. py:method:: dates_todict(self, dates)

      Set dates by start and end dates as a tuple, (start, end).


   .. py:method:: dates_tolist(self, dates)

      Correct dates for Daymet accounting for leap years.

      Daymet doesn't account for leap years and removes Dec 31 when
      it's leap year.

      :Parameters: **dates** (:class:`tuple`) -- Target start and end dates.

      :returns: :class:`list` -- All the dates in the Daymet database within the provided date range.


   .. py:method:: years_todict(self, years)

      Set date by list of year(s).


   .. py:method:: years_tolist(self, years)

      Correct dates for Daymet accounting for leap years.

      Daymet doesn't account for leap years and removes Dec 31 when
      it's leap year.

      :Parameters: **years** (:class:`list`) -- A list of target years.

      :returns: :class:`list` -- All the dates in the Daymet database within the provided date range.



:py:mod:`pydaymet.pydaymet`
===========================

.. py:module:: pydaymet.pydaymet

.. autoapi-nested-parse::

   Access the Daymet database for both single single pixel and gridded queries.



Module Contents
---------------

.. py:function:: get_bycoords(coords, dates, crs = DEF_CRS, variables = None, region = 'na', time_scale = 'daily', pet = None, pet_params = None, ssl = None, expire_after = EXPIRE, disable_caching = False)

   Get point-data from the Daymet database at 1-km resolution.

   This function uses THREDDS data service to get the coordinates
   and supports getting monthly and annual summaries of the climate
   data directly from the server.

   :Parameters: * **coords** (:class:`tuple`) -- Coordinates of the location of interest as a tuple (lon, lat)
                * **dates** (:class:`tuple` or :class:`list`, *optional*) -- Start and end dates as a tuple (start, end) or a list of years ``[2001, 2010, ...]``.
                * **crs** (:class:`str`, *optional*) -- The CRS of the input geometry, defaults to ``"epsg:4326"``.
                * **variables** (:class:`str` or :class:`list`) -- List of variables to be downloaded. The acceptable variables are:
                  ``tmin``, ``tmax``, ``prcp``, ``srad``, ``vp``, ``swe``, ``dayl``
                  Descriptions can be found `here <https://daymet.ornl.gov/overview>`__.
                * **region** (:class:`str`, *optional*) -- Target region in the US, defaults to ``na``. Acceptable values are:

                  * ``na``: Continental North America
                  * ``hi``: Hawaii
                  * ``pr``: Puerto Rico
                * **time_scale** (:class:`str`, *optional*) -- Data time scale which can be ``daily``, ``monthly`` (monthly summaries),
                  or ``annual`` (annual summaries). Defaults to ``daily``.
                * **pet** (:class:`str`, *optional*) -- Method for computing PET. Supported methods are
                  ``penman_monteith``, ``priestley_taylor``, ``hargreaves_samani``, and
                  None (don't compute PET). The ``penman_monteith`` method is based on
                  :footcite:t:`Allen_1998` assuming that soil heat flux density is zero.
                  The ``priestley_taylor`` method is based on
                  :footcite:t:`Priestley_1972` assuming that soil heat flux density is zero.
                  The ``hargreaves_samani`` method is based on :footcite:t:`Hargreaves_1982`.
                  Defaults to ``None``.
                * **pet_params** (:class:`dict`, *optional*) -- Model-specific parameters as a dictionary that is passed to the PET function.
                  Defaults to ``None``.
                * **ssl** (:class:`bool` or :class:`SSLContext`, *optional*) -- SSLContext to use for the connection, defaults to None. Set to False to disable
                  SSL cetification verification.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`pandas.DataFrame` -- Daily climate data for a location.

   .. rubric:: Examples

   >>> import pydaymet as daymet
   >>> coords = (-1431147.7928, 318483.4618)
   >>> dates = ("2000-01-01", "2000-12-31")
   >>> clm = daymet.get_bycoords(
   ...     coords,
   ...     dates,
   ...     crs="epsg:3542",
   ...     pet="hargreaves_samani",
   ...     ssl=False
   ... )
   >>> clm["pet (mm/day)"].mean()
   3.713

   .. rubric:: References

   .. footbibliography::


.. py:function:: get_bygeom(geometry, dates, crs = DEF_CRS, variables = None, region = 'na', time_scale = 'daily', pet = None, pet_params = None, ssl = None, expire_after = EXPIRE, disable_caching = False)

   Get gridded data from the Daymet database at 1-km resolution.

   :Parameters: * **geometry** (:class:`Polygon`, :class:`MultiPolygon`, or :class:`bbox`) -- The geometry of the region of interest.
                * **dates** (:class:`tuple` or :class:`list`, *optional*) -- Start and end dates as a tuple (start, end) or a list of years [2001, 2010, ...].
                * **crs** (:class:`str`, *optional*) -- The CRS of the input geometry, defaults to epsg:4326.
                * **variables** (:class:`str` or :class:`list`) -- List of variables to be downloaded. The acceptable variables are:
                  ``tmin``, ``tmax``, ``prcp``, ``srad``, ``vp``, ``swe``, ``dayl``
                  Descriptions can be found `here <https://daymet.ornl.gov/overview>`__.
                * **region** (:class:`str`, *optional*) -- Region in the US, defaults to na. Acceptable values are:

                  * na: Continental North America
                  * hi: Hawaii
                  * pr: Puerto Rico
                * **time_scale** (:class:`str`, *optional*) -- Data time scale which can be daily, monthly (monthly average),
                  or annual (annual average). Defaults to daily.
                * **pet** (:class:`str`, *optional*) -- Method for computing PET. Supported methods are
                  ``penman_monteith``, ``priestley_taylor``, ``hargreaves_samani``, and
                  None (don't compute PET). The ``penman_monteith`` method is based on
                  :footcite:t:`Allen_1998` assuming that soil heat flux density is zero.
                  The ``priestley_taylor`` method is based on
                  :footcite:t:`Priestley_1972` assuming that soil heat flux density is zero.
                  The ``hargreaves_samani`` method is based on :footcite:t:`Hargreaves_1982`.
                  Defaults to ``None``.
                * **pet_params** (:class:`dict`, *optional*) -- Model-specific parameters as a dictionary that is passed to the PET function.
                  Defaults to ``None``.
                * **ssl** (:class:`bool` or :class:`SSLContext`, *optional*) -- SSLContext to use for the connection, defaults to None. Set to False to disable
                  SSL cetification verification.
                * **expire_after** (:class:`int`, *optional*) -- Expiration time for response caching in seconds, defaults to -1 (never expire).
                * **disable_caching** (:class:`bool`, *optional*) -- If ``True``, disable caching requests, defaults to False.

   :returns: :class:`xarray.Dataset` -- Daily climate data within the target geometry.

   .. rubric:: Examples

   >>> from shapely.geometry import Polygon
   >>> import pydaymet as daymet
   >>> geometry = Polygon(
   ...     [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
   ... )
   >>> clm = daymet.get_bygeom(geometry, 2010, variables="tmin", time_scale="annual")
   >>> clm["tmin"].mean().compute().item()
   1.361

   .. rubric:: References

   .. footbibliography::


.. include:: ../../../pynhd/HISTORY.rst
.. include:: ../../../pygeoutils/HISTORY.rst
.. include:: ../../../pygeoogc/HISTORY.rst
.. include:: ../../../pydaymet/HISTORY.rst
.. include:: ../../../async_retriever/HISTORY.rst
.. include:: ../../../py3dep/HISTORY.rst
.. include:: ../../../pygeohydro/HISTORY.rst
PyNHD: Navigate and subset NHDPlus database
-------------------------------------------

.. image:: https://img.shields.io/pypi/v/pynhd.svg
    :target: https://pypi.python.org/pypi/pynhd
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/pynhd.svg
    :target: https://anaconda.org/conda-forge/pynhd
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/pynhd/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pynhd
    :alt: CodeCov

.. image:: https://img.shields.io/pypi/pyversions/pynhd.svg
    :target: https://pypi.python.org/pypi/pynhd
    :alt: Python Versions

.. image:: https://pepy.tech/badge/pynhd
    :target: https://pepy.tech/project/pynhd
    :alt: Downloads

|

.. image:: https://www.codefactor.io/repository/github/cheginit/pynhd/badge
   :target: https://www.codefactor.io/repository/github/cheginit/pynhd
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

|

Features
--------

PyNHD is a part of `HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services.

This package provides access to
`WaterData <https://labs.waterdata.usgs.gov/geoserver/web/wicket/bookmarkable/org.geoserver.web.demo.MapPreviewPage?1>`__,
the National Map's `NHDPlus HR <https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer>`__,
`NLDI <https://labs.waterdata.usgs.gov/about-nldi/>`__,
and `PyGeoAPI <https://labs.waterdata.usgs.gov/api/nldi/pygeoapi>`__ web services. These web services
can be used to navigate and extract vector data from NHDPlus V2 (both medium- and
high-resolution) database such as catchments, HUC8, HUC12, GagesII, flowlines, and water bodies.
Moreover, PyNHD gives access to an item on `ScienceBase <https://sciencebase.usgs.gov>`_ called
`Select Attributes for NHDPlus Version 2.1 Reach Catchments and Modified Network Routed Upstream Watersheds for the Conterminous United States <https://www.sciencebase.gov/catalog/item/5669a79ee4b08895842a1d47>`_.
This item provides over 30 attributes at catchment-scale based on NHDPlus ComIDs.
These attributes are available in three categories:

1. Local (`local`): For individual reach catchments,
2. Total (`upstream_acc`): For network-accumulated values using total cumulative drainage area,
3. Divergence (`div_routing`): For network-accumulated values using divergence-routed.

Moreover, the PyGeoAPI service provides four functionalities:

1. ``flow_trace``: Trace flow from a starting point to up/downstream direction.
2. ``split_catchment``: Split the local catchment of a point of interest at the point's location.
3. ``elevation_profile``: Extract elevation profile along a flow path between two points.
4. ``cross_section``: Extract cross-section at a point of interest along a flow line.

A list of these attributes for each characteristic type can be accessed using ``nhdplus_attrs``
function.

Similarly, PyNHD uses `this <https://www.hydroshare.org/resource/6092c8a62fac45be97a09bfd0b0bf726/>`__
item on Hydroshare to get ComID-linked NHDPlus Value Added Attributes. This dataset includes
slope and roughness, among other attributes, for all the flowlines. You can use ``nhdplus_vaa``
function to get this dataset.

Additionally, PyNHD offers some extra utilities for processing the flowlines:

- ``prepare_nhdplus``: For cleaning up the dataframe by, for example, removing tiny networks,
  adding a ``to_comid`` column, and finding a terminal flowlines if it doesn't exist.
- ``topoogical_sort``: For sorting the river network topologically which is useful for routing
  and flow accumulation.
- ``vector_accumulation``: For computing flow accumulation in a river network. This function
  is generic, and any routing method can be plugged in.

These utilities are developed based on an ``R`` package called
`nhdplusTools <https://github.com/USGS-R/nhdplusTools>`__.

All functions and classes that request data from web services use ``async_retriever``
that offers response caching. By default, the expiration time is set to never expire.
All these functions and classes have two optional parameters for controlling the cache:
``expire_after`` and ``disable_caching``. You can use ``expire_after`` to set the expiration
time in seconds. If ``expire_after`` is set to ``-1``, the cache will never expire (default).
You can use ``disable_caching`` if you don't want to use the cached responses. The cached
responses are stored in the ``./cache/aiohttp_cache.sqlite`` file.

You can find some example notebooks `here <https://github.com/cheginit/HyRiver-examples>`__.

Furthermore, you can try using PyNHD without even installing it on your system by
clicking on the binder badge below the PyNHD banner. A JupyterLab instance
with the software stack pre-installed and all example notebooks will be launched
in your web browser, and you can start coding!

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pynhd/issues>`__.

Installation
------------

You can install PyNHD using ``pip`` after installing ``libgdal`` on your system
(for example, in Ubuntu run ``sudo apt install libgdal-dev``):

.. code-block:: console

    $ pip install pynhd

Alternatively, PyNHD can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__
or `Mamba <https://github.com/conda-forge/miniforge>`__:

.. code-block:: console

    $ conda install -c conda-forge pynhd

Quick start
-----------

Let's explore the capabilities of ``NLDI``. We need to instantiate the class first:

.. code:: python

    from pynhd import NLDI, WaterData, NHDPlusHR
    import pynhd as nhd

First, let’s get the watershed geometry of the contributing basin of a
USGS station using ``NLDI``:

.. code:: python

    nldi = NLDI()
    station_id = "01031500"

    basin = nldi.get_basins(station_id)

The ``navigate_byid`` class method can be used to navigate NHDPlus in
both upstream and downstream of any point in the database. Let’s get ComIDs and flowlines
of the tributaries and the main river channel in the upstream of the station.

.. code:: python

    flw_main = nldi.navigate_byid(
        fsource="nwissite",
        fid=f"USGS-{station_id}",
        navigation="upstreamMain",
        source="flowlines",
        distance=1000,
    )

    flw_trib = nldi.navigate_byid(
        fsource="nwissite",
        fid=f"USGS-{station_id}",
        navigation="upstreamTributaries",
        source="flowlines",
        distance=1000,
    )

We can get other USGS stations upstream (or downstream) of the station
and even set a distance limit (in km):

.. code:: python

    st_all = nldi.navigate_byid(
        fsource="nwissite",
        fid=f"USGS-{station_id}",
        navigation="upstreamTributaries",
        source="nwissite",
        distance=1000,
    )

    st_d20 = nldi.navigate_byid(
        fsource="nwissite",
        fid=f"USGS-{station_id}",
        navigation="upstreamTributaries",
        source="nwissite",
        distance=20,
    )

Now, let’s get the
`HUC12 pour points <https://www.sciencebase.gov/catalog/item/5762b664e4b07657d19a71ea>`__:

.. code:: python

    pp = nldi.navigate_byid(
        fsource="nwissite",
        fid=f"USGS-{station_id}",
        navigation="upstreamTributaries",
        source="huc12pp",
        distance=1000,
    )

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/nhdplus_navigation.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/nhdplus.ipynb
    :align: center

Also, we can get the slope data for each river segment from NHDPlus VAA database:

.. code:: python

    vaa = nhd.nhdplus_vaa("input_data/nhdplus_vaa.parquet")

    flw_trib["comid"] = pd.to_numeric(flw_trib.nhdplus_comid)
    slope = gpd.GeoDataFrame(
        pd.merge(flw_trib, vaa[["comid", "slope"]], left_on="comid", right_on="comid"),
        crs=flw_trib.crs,
    )
    slope[slope.slope < 0] = np.nan

Now, let's explore the PyGeoAPI capabilities:

.. code:: python

    pygeoapi = PyGeoAPI()

    trace = pygeoapi.flow_trace(
        (1774209.63, 856381.68), crs="ESRI:102003", raindrop=False, direction="none"
    )

    split = pygeoapi.split_catchment((-73.82705, 43.29139), crs="epsg:4326", upstream=False)

    profile = pygeoapi.elevation_profile(
        [(-103.801086, 40.26772), (-103.80097, 40.270568)], numpts=101, dem_res=1, crs="epsg:4326"
    )

    section = pygeoapi.cross_section((-103.80119, 40.2684), width=1000.0, numpts=101, crs="epsg:4326")

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/split_catchment.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/pygeoapi.ipynb
    :align: center

Next, we retrieve the medium- and high-resolution flowlines within the bounding box of our
watershed and compare them. Moreover, Since several web services offer access to NHDPlus database,
``NHDPlusHR`` has an argument for selecting a service and also an argument for automatically
switching between services.

.. code:: python

    mr = WaterData("nhdflowline_network")
    nhdp_mr = mr.bybox(basin.geometry[0].bounds)

    hr = NHDPlusHR("networknhdflowline", service="hydro", auto_switch=True)
    nhdp_hr = hr.bygeom(basin.geometry[0].bounds)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/hr_mr.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/nhdplus.ipynb
    :align: center

Moreover, ``WaterData`` can find features within a given radius (in meters) of a point:

.. code:: python

    eck4 = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    coords = (-5727797.427596455, 5584066.49330473)
    rad = 5e3
    flw_rad = mr.bydistance(coords, rad, loc_crs=eck4)
    flw_rad = flw_rad.to_crs(eck4)

Instead of getting all features within a radius of the coordinate, we can snap to the closest
flowline using NLDI:

.. code:: python

    comid_closest = nldi.comid_byloc((x, y), eck4)
    flw_closest = nhdp_mr.byid("comid", comid_closest.comid.values[0])


.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/nhdplus_radius.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/nhdplus.ipynb
    :align: center

Since NHDPlus HR is still at the pre-release stage let's use the MR flowlines to
demonstrate the vector-based accumulation.
Based on a topological sorted river network
``pynhd.vector_accumulation`` computes flow accumulation in the network.
It returns a dataframe which is sorted from upstream to downstream that
shows the accumulated flow in each node.

PyNHD has a utility called ``prepare_nhdplus`` that identifies such
relationship among other things such as fixing some common issues with
NHDPlus flowlines. But first we need to get all the NHDPlus attributes
for each ComID since ``NLDI`` only provides the flowlines’ geometries
and ComIDs which is useful for navigating the vector river network data.
For getting the NHDPlus database we use ``WaterData``. Let’s use the
``nhdflowline_network`` layer to get required info.

.. code:: python

    wd = WaterData("nhdflowline_network")

    comids = flw_trib.nhdplus_comid.to_list()
    nhdp_trib = wd.byid("comid", comids)
    flw = nhd.prepare_nhdplus(nhdp_trib, 0, 0, purge_non_dendritic=False)

To demonstrate the use of routing, let's use ``nhdplus_attrs`` function to get list of available
NHDPlus attributes

.. code:: python

    char = "CAT_RECHG"
    area = "areasqkm"

    local = nldi.getcharacteristic_byid(comids, "local", char_ids=char)
    flw = flw.merge(local[char], left_on="comid", right_index=True)


    def runoff_acc(qin, q, a):
        return qin + q * a


    flw_r = flw[["comid", "tocomid", char, area]]
    runoff = nhd.vector_accumulation(flw_r, runoff_acc, char, [char, area])


    def area_acc(ain, a):
        return ain + a


    flw_a = flw[["comid", "tocomid", area]]
    areasqkm = nhd.vector_accumulation(flw_a, area_acc, area, [area])

    runoff /= areasqkm

Since these are catchment-scale characteristic, let’s get the catchments
then add the accumulated characteristic as a new column and plot the
results.

.. code:: python

    wd = WaterData("catchmentsp")
    catchments = wd.byid("featureid", comids)

    c_local = catchments.merge(local, left_on="featureid", right_index=True)
    c_acc = catchments.merge(runoff, left_on="featureid", right_index=True)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/flow_accumulation.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/nhdplus.ipynb
    :align: center

More examples can be found `here <https://pygeohydro.readthedocs.io/en/latest/examples.html>`__.
PyGeoUtils: Utilities for (Geo)JSON and (Geo)TIFF Conversion
------------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/pygeoutils.svg
    :target: https://pypi.python.org/pypi/pygeoutils
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/pygeoutils.svg
    :target: https://anaconda.org/conda-forge/pygeoutils
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/pygeoutils/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pygeoutils
    :alt: CodeCov

.. image:: https://img.shields.io/pypi/pyversions/pygeoutils.svg
    :target: https://pypi.python.org/pypi/pygeoutils
    :alt: Python Versions

.. image:: https://pepy.tech/badge/pygeoutils
    :target: https://pepy.tech/project/pygeoutils
    :alt: Downloads

|

.. image:: https://www.codefactor.io/repository/github/cheginit/pygeoutils/badge
   :target: https://www.codefactor.io/repository/github/cheginit/pygeoutils
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

|

Features
--------

PyGeoUtils is a part of `HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services. This package provides
utilities for manipulating (Geo)JSON and (Geo)TIFF responses from web services.
These utilities are:

- ``json2geodf``: For converting (Geo)JSON objects to GeoPandas dataframe.
- ``arcgis2geojson``: For converting ESRIGeoJSON to the standard GeoJSON format.
- ``gtiff2xarray``: For converting (Geo)TIFF objects to `xarray <https://xarray.pydata.org/>`__.
  datasets.
- ``xarray2geodf``: For converting ``xarray.DataArray`` to a ``geopandas.GeoDataFrame``, i.e.,
  vectorization.
- ``xarray_geomask``: For masking a ``xarray.Dataset`` or ``xarray.DataArray`` using a polygon.

All these functions handle all necessary CRS transformations.

You can find some example notebooks `here <https://github.com/cheginit/HyRiver-examples>`__.

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pygeoutils/issues>`__.

Installation
------------

You can install PyGeoUtils using ``pip`` after installing ``libgdal`` on your system
(for example, in Ubuntu run ``sudo apt install libgdal-dev``). Moreover, PyGeoUtils has an optional
dependency for using persistent caching, ``requests-cache``. We highly recommend to install
this package as it can significantly speedup send/receive queries. You don't have to change
anything in your code, since PyGeoUtils under-the-hood looks for ``requests-cache`` and
if available, it will automatically use persistent caching:

.. code-block:: console

    $ pip install pygeoutils

Alternatively, PyGeoUtils can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__:

.. code-block:: console

    $ conda install -c conda-forge pygeoutils

Quick start
-----------

To demonstrate capabilities of PyGeoUtils let's use
`PyGeoOGC <https://github.com/cheginit/pygeoogc>`__ to access
`National Wetlands Inventory <https://www.fws.gov/wetlands/>`__ from WMS, and
`FEMA National Flood Hazard <https://www.fema.gov/national-flood-hazard-layer-nfhl>`__
via WFS, then convert the output to ``xarray.Dataset`` and ``GeoDataFrame``, respectively.

.. code-block:: python

    import pygeoutils as geoutils
    from pygeoogc import WFS, WMS, ServiceURL
    from shapely.geometry import Polygon


    geometry = Polygon(
        [
            [-118.72, 34.118],
            [-118.31, 34.118],
            [-118.31, 34.518],
            [-118.72, 34.518],
            [-118.72, 34.118],
        ]
    )
    crs = "epsg:4326"

    wms = WMS(
        ServiceURL().wms.mrlc,
        layers="NLCD_2011_Tree_Canopy_L48",
        outformat="image/geotiff",
        crs=crs,
    )
    r_dict = wms.getmap_bybox(
        geometry.bounds,
        1e3,
        box_crs=crs,
    )
    canopy = geoutils.gtiff2xarray(r_dict, geometry, crs)

    mask = canopy > 60
    canopy_gdf = geoutils.xarray2geodf(canopy, "float32", mask)

    url_wfs = "https://hazards.fema.gov/gis/nfhl/services/public/NFHL/MapServer/WFSServer"
    wfs = WFS(
        url_wfs,
        layer="public_NFHL:Base_Flood_Elevations",
        outformat="esrigeojson",
        crs="epsg:4269",
    )
    r = wfs.getfeature_bybox(geometry.bounds, box_crs=crs)
    flood = geoutils.json2geodf(r.json(), "epsg:4269", crs)
PyGeoOGC: Retrieve Data from RESTful, WMS, and WFS Services
-----------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/pygeoogc.svg
    :target: https://pypi.python.org/pypi/pygeoogc
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/pygeoogc.svg
    :target: https://anaconda.org/conda-forge/pygeoogc
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/pygeoogc/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pygeoogc
    :alt: CodeCov

.. image:: https://img.shields.io/pypi/pyversions/pygeoogc.svg
    :target: https://pypi.python.org/pypi/pygeoogc
    :alt: Python Versions

.. image:: https://pepy.tech/badge/pygeoogc
    :target: https://pepy.tech/project/pygeoogc
    :alt: Downloads

|

.. image:: https://img.shields.io/badge/security-bandit-green.svg
    :target: https://github.com/PyCQA/bandit
    :alt: Security Status

.. image:: https://www.codefactor.io/repository/github/cheginit/pygeoogc/badge
   :target: https://www.codefactor.io/repository/github/cheginit/pygeoogc
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

|

Features
--------

PyGeoOGC is a part of `HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services. This package provides
general interfaces to web services that are based on
`ArcGIS RESTful <https://en.wikipedia.org/wiki/Representational_state_transfer>`__,
`WMS <https://en.wikipedia.org/wiki/Web_Map_Service>`__, and
`WFS <https://en.wikipedia.org/wiki/Web_Feature_Service>`__. Although
all these web service have limits on the number of features per requests (e.g., 1000
object IDs for a RESTful request or 8 million pixels for a WMS request), PyGeoOGC divides
requests into smaller chunks, under-the-hood, and then merges the results.

All functions and classes that request data from web services use ``async_retriever``
that offers response caching. By default, the expiration time is set to never expire.
All these functions and classes have two optional parameters for controlling the cache:
``expire_after`` and ``disable_caching``. You can use ``expire_after`` to set the expiration
time in seconds. If ``expire_after`` is set to ``-1``, the cache will never expire (default).
You can use ``disable_caching`` if you don't want to use the cached responses. The cached
responses are stored in the ``./cache/aiohttp_cache.sqlite`` file.

There is also an inventory of URLs for some of these web services in form of a class called
``ServiceURL``. These URLs are in four categories: ``ServiceURL().restful``,
``ServiceURL().wms``, ``ServiceURL().wfs``, and ``ServiceURL().http``. These URLs provide you
with some examples of the services that PyGeoOGC supports. All the URLs are read from a YAML
file located `here <pygeoogc/static/urls.yml>`_. If you have success using PyGeoOGC with a web
service please consider submitting a request to be added to this URL inventory, located at
``pygeoogc/static/urls.yml``.

PyGeoOGC has three main classes:

* ``ArcGISRESTful``: This class can be instantiated by providing the target layer URL.
  For example, for getting Watershed Boundary Data we can use ``ServiceURL().restful.wbd``.
  By looking at the web service's
  `website <https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer>`_
  we see that there are nine layers. For example, 1 for 2-digit HU (Region), 6 for 12-digit HU
  (Subregion), and so on. We can pass the URL to the target layer directly, like this
  ``f"{ServiceURL().restful.wbd}/6"`` or as a separate argument via ``layer``.

  Afterward, we request for the data in two steps. First, we need to get
  the target object IDs using ``oids_bygeom`` (within a geometry), ``oids_byfield`` (specific
  field IDs), or ``oids_bysql`` (any valid SQL 92 WHERE clause) class methods. Then, we can get
  the target features using ``get_features`` class method. The returned response can be converted
  into a GeoDataFrame using ``json2geodf`` function from
  `PyGeoUtils <https://github.com/cheginit/pygeoutils>`__.

* ``WMS``: Instantiation of this class requires at least 3 arguments: service URL, layer
  name(s), and output format. Additionally, target CRS and the web service version can be provided.
  Upon instantiation, we can use ``getmap_bybox`` method class to get the target raster data
  within a bounding box. The box can be in any valid CRS and if it is different from the default
  CRS, ``EPSG:4326``, it should be passed using ``box_crs`` argument. The service response can be
  converted into a ``xarray.Dataset`` using ``gtiff2xarray`` function from PyGeoUtils.

* ``WFS``: Instantiation of this class is similar to ``WMS``. The only difference is that
  only one layer name can be passed. Upon instantiation there are three ways to get the data:

  - ``getfeature_bybox``: Get all the target features within a bounding box in any valid CRS.
  - ``getfeature_byid``: Get all the target features based on the IDs. Note that two arguments
    should be provided: ``featurename``, and ``featureids``. You can get a list of valid feature
    names using ``get_validnames`` class method.
  - ``getfeature_byfilter``: Get the data based on any valid
    `CQL <https://docs.geoserver.org/latest/en/user/tutorials/cql/cql_tutorial.html>`__ filter.

  You can convert the returned response of this function to a ``GeoDataFrame`` using ``json2geodf``
  function from PyGeoUtils package.

You can find some example notebooks `here <https://github.com/cheginit/HyRiver-examples>`__.

Furthermore, you can try using PyGeoOGC without even installing it on your system by
clicking on the binder badge below the PyGeoOGC banner. A JupyterLab instance
with the software stack pre-installed and all example notebooks will be launched
in your web browser, and you can start coding!

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback.
Contributions are most welcome.

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pygeoogc/issues>`__.

Installation
------------

You can install PyGeoOGC using ``pip``:

.. code-block:: console

    $ pip install pygeoogc

Alternatively, PyGeoOGC can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__
or `Mamba <https://github.com/conda-forge/miniforge>`__:

.. code-block:: console

    $ conda install -c conda-forge pygeoogc

Quick start
-----------

We can access
`NHDPlus HR <https://edits.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/NHDPlus_HR/MapServer>`__
via RESTful service,
`National Wetlands Inventory <https://www.fws.gov/wetlands/>`__ from WMS, and
`FEMA National Flood Hazard <https://www.fema.gov/national-flood-hazard-layer-nfhl>`__
via WFS. The output for these functions are of type ``requests.Response`` that
can be converted to ``GeoDataFrame`` or ``xarray.Dataset`` using
`PyGeoUtils <https://github.com/cheginit/pygeoutils>`__.

Let's start the National Map's NHDPlus HR web service. We can query the flowlines that are
within a geometry as follows:

.. code-block:: python

    from pygeoogc import ArcGISRESTful, WFS, WMS, ServiceURL
    import pygeoutils as geoutils
    from pynhd import NLDI

    basin_geom = NLDI().get_basins("01031500").geometry[0]

    hr = ArcGISRESTful(ServiceURL().restful.nhdplushr, 2, outformat="json")

    resp = hr.get_features(hr.oids_bygeom(basin_geom, "epsg:4326"))
    flowlines = geoutils.json2geodf(resp)

Note ``oids_bygeom`` has three additional arguments: ``sql_clause``, ``spatial_relation``,
and ``distance``. We can use ``sql_clause`` for passing any valid SQL WHERE clauses and
``spatial_relation`` for specifying the target predicate such as
intersect, contain, cross, etc. The default predicate is intersect
(``esriSpatialRelIntersects``). Additionally, we can use ``distance`` for specifying the buffer
distance from the input geometry for getting features.

We can also submit a query based on IDs of any valid field in the database. If the measure
property is desired you can pass ``return_m`` as ``True`` to the ``get_features`` class method:

.. code-block:: python

    oids = hr.oids_byfield("PERMANENT_IDENTIFIER", ["103455178", "103454362", "103453218"])
    resp = hr.get_features(oids, return_m=True)
    flowlines = geoutils.json2geodf(resp)

Additionally, any valid SQL 92 WHERE clause can be used. For more details look
`here <https://developers.arcgis.com/rest/services-reference/query-feature-service-.htm#ESRI_SECTION2_07DD2C5127674F6A814CE6C07D39AD46>`__.
For example, let's limit our first request to only include catchments with
areas larger than 0.5 sqkm.

.. code-block:: python

    oids = hr.oids_bygeom(basin_geom, geo_crs="epsg:4326", sql_clause="AREASQKM > 0.5")
    resp = hr.get_features(oids)
    catchments = geoutils.json2geodf(resp)

A WMS-based example is shown below:

.. code-block:: python

    wms = WMS(
        ServiceURL().wms.fws,
        layers="0",
        outformat="image/tiff",
        crs="epsg:3857",
    )
    r_dict = wms.getmap_bybox(
        basin_geom.bounds,
        1e3,
        box_crs="epsg:4326",
    )
    wetlands = geoutils.gtiff2xarray(r_dict, basin_geom, "epsg:4326")

Query from a WFS-based web service can be done either within a bounding box or using
any valid `CQL filter <https://docs.geoserver.org/stable/en/user/tutorials/cql/cql_tutorial.html>`__.

.. code-block:: python

    wfs = WFS(
        ServiceURL().wfs.fema,
        layer="public_NFHL:Base_Flood_Elevations",
        outformat="esrigeojson",
        crs="epsg:4269",
    )
    r = wfs.getfeature_bybox(basin_geom.bounds, box_crs="epsg:4326")
    flood = geoutils.json2geodf(r.json(), "epsg:4269", "epsg:4326")

    layer = "wmadata:huc08"
    wfs = WFS(
        ServiceURL().wfs.waterdata,
        layer=layer,
        outformat="application/json",
        version="2.0.0",
        crs="epsg:4269",
    )
    r = wfs.getfeature_byfilter(f"huc8 LIKE '13030%'")
    huc8 = geoutils.json2geodf(r.json(), "epsg:4269", "epsg:4326")

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/sql_clause.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/webservices.ipynb

PyDaymet: Daily climate data through Daymet
-------------------------------------------

.. image:: https://img.shields.io/pypi/v/pydaymet.svg
    :target: https://pypi.python.org/pypi/pydaymet
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/pydaymet.svg
    :target: https://anaconda.org/conda-forge/pydaymet
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/pydaymet/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pydaymet
    :alt: CodeCov

.. image:: https://img.shields.io/pypi/pyversions/pydaymet.svg
    :target: https://pypi.python.org/pypi/pydaymet
    :alt: Python Versions

.. image:: https://pepy.tech/badge/pydaymet
    :target: https://pepy.tech/project/pydaymet
    :alt: Downloads

|

.. image:: https://www.codefactor.io/repository/github/cheginit/pydaymet/badge
   :target: https://www.codefactor.io/repository/github/cheginit/pydaymet
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

|

Features
--------

PyDaymet is a part of `HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services. This package provides
access to climate data from
`Daymet V4 <https://daac.ornl.gov/DAYMET/guides/Daymet_Daily_V4.html>`__ database using NetCDF
Subset Service (NCSS). Both single pixel (using ``get_bycoords`` function) and gridded data (using
``get_bygeom``) are supported which are returned as
``pandas.DataFrame`` and ``xarray.Dataset``, respectively. Climate data is available for North
America, Hawaii from 1980, and Puerto Rico from 1950 at three time scales: daily, monthly,
and annual. Additionally, PyDaymet can compute Potential EvapoTranspiration (PET)
using three methods: ``penman_monteith``, ``priestley_taylor``, and ``hargreaves_samani`` for
both single pixel and gridded data.

To fully utilize the capabilities of the NCSS, under-the-hood, PyDaymet uses
`AsyncRetriever <https://github.com/cheginit/async_retriever>`__
for retrieving Daymet data asynchronously with persistent caching. This improves the reliability
and speed of data retrieval significantly.

You can try using PyDaymet without installing it on you system by clicking on the binder badge
below the PyDaymet banner. A Jupyter notebook instance with the stack
pre-installed will be launched in your web browser and you can start coding!

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pydaymet/issues>`__.

Installation
------------

You can install PyDaymet using ``pip`` after installing ``libgdal`` on your system
(for example, in Ubuntu run ``sudo apt install libgdal-dev``):

.. code-block:: console

    $ pip install pydaymet

Alternatively, PyDaymet can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__:

.. code-block:: console

    $ conda install -c conda-forge pydaymet

Quick start
-----------

You can use PyDaymet using command-line or as a Python library. The commanda-line
provides access to two functionality:

- Getting gridded climate data: You must create a ``geopandas.GeoDataFrame`` that contains
  the geometries of the target locations. This dataframe must have four columns:
  ``id``, ``start``, ``end``, ``geometry``. The ``id`` column is used as
  filenames for saving the obtained climate data to a NetCDF (``.nc``) file. The ``start``
  and ``end`` columns are starting and ending dates of the target period. Then,
  you must save the dataframe as a shapefile (``.shp``) or geopackage (``.gpkg``) with
  CRS attribute.
- Getting single pixel climate data: You must create a CSV file that
  contains coordinates of the target locations. This file must have at four columns:
  ``id``, ``start``, ``end``, ``lon``, and ``lat``. The ``id`` column is used as filenames
  for saving the obtained climate data to a CSV (``.csv``) file. The ``start`` and ``end``
  columns are the same as the ``geometry`` command. The ``lon`` and ``lat`` columns are
  the longitude and latitude coordinates of the target locations.

.. code-block:: console

    $ pydaymet -h
    Usage: pydaymet [OPTIONS] COMMAND [ARGS]...

    Command-line interface for PyDaymet.

    Options:
    -h, --help  Show this message and exit.

    Commands:
    coords    Retrieve climate data for a list of coordinates.
    geometry  Retrieve climate data for a dataframe of geometries.

The ``coords`` sub-command is as follows:

.. code-block:: console

    $ pydaymet coords -h
    Usage: pydaymet coords [OPTIONS] FPATH

    Retrieve climate data for a list of coordinates.

    FPATH: Path to a csv file with four columns:
        - ``id``: Feature identifiers that daymet uses as the output netcdf filenames.
        - ``start``: Start time.
        - ``end``: End time.
        - ``lon``: Longitude of the points of interest.
        - ``lat``: Latitude of the points of interest.
        - ``time_scale``: (optional) Time scale, either ``daily`` (default), ``monthly`` or ``annual``.
        - ``pet``: (optional) Method to compute PET. Suppoerted methods are:
                   ``penman_monteith``, ``hargreaves_samani``, ``priestley_taylor``, and ``none`` (default).
        - ``alpha``: (optional) Alpha parameter for Priestley-Taylor method for computing PET. Defaults to 1.26.

    Examples:
        $ cat coords.csv
        id,lon,lat,start,end,pet
        california,-122.2493328,37.8122894,2012-01-01,2014-12-31,hargreaves_samani
        $ pydaymet coords coords.csv -v prcp -v tmin

    Options:
    -v, --variables TEXT  Target variables. You can pass this flag multiple
                            times for multiple variables.

    -s, --save_dir PATH   Path to a directory to save the requested files.
                            Extension for the outputs is .nc for geometry and .csv
                            for coords.

    -h, --help            Show this message and exit.

And, the ``geometry`` sub-command is as follows:

.. code-block:: console

    $ pydaymet geometry -h
    Usage: pydaymet geometry [OPTIONS] FPATH

    Retrieve climate data for a dataframe of geometries.

    FPATH: Path to a shapefile (.shp) or geopackage (.gpkg) file.
    This file must have four columns and contain a ``crs`` attribute:
        - ``id``: Feature identifiers that daymet uses as the output netcdf filenames.
        - ``start``: Start time.
        - ``end``: End time.
        - ``geometry``: Target geometries.
        - ``time_scale``: (optional) Time scale, either ``daily`` (default), ``monthly`` or ``annual``.
        - ``pet``: (optional) Method to compute PET. Suppoerted methods are:
                   ``penman_monteith``, ``hargreaves_samani``, ``priestley_taylor``, and ``none`` (default).
        - ``alpha``: (optional) Alpha parameter for Priestley-Taylor method for computing PET. Defaults to 1.26.

    Examples:
        $ pydaymet geometry geo.gpkg -v prcp -v tmin

    Options:
    -v, --variables TEXT  Target variables. You can pass this flag multiple
                            times for multiple variables.

    -s, --save_dir PATH   Path to a directory to save the requested files.
                            Extension for the outputs is .nc for geometry and .csv
                            for coords.

    -h, --help            Show this message and exit.

Now, let's see how we can use PyDaymet as a library.

PyDaymet offers two functions for getting climate data; ``get_bycoords`` and ``get_bygeom``.
The arguments of these functions are identical except the first argument where the latter
should be polygon and the former should be a coordinate (a tuple of length two as in (x, y)).
The input geometry or coordinate can be in any valid CRS (defaults to EPSG:4326). The ``dates``
argument can be either a tuple of length two like ``(start_str, end_str)`` or a list of years
like ``[2000, 2005]``. It is noted that both functions have a ``pet`` flag for computing PET.
Additionally, we can pass ``time_scale`` to get daily, monthly or annual summaries. This flag
by default is set to daily.

.. code-block:: python

    from pynhd import NLDI
    import pydaymet as daymet

    geometry = NLDI().get_basins("01031500").geometry[0]

    var = ["prcp", "tmin"]
    dates = ("2000-01-01", "2000-06-30")

    daily = daymet.get_bygeom(geometry, dates, variables=var, pet="priestley_taylor")
    monthly = daymet.get_bygeom(geometry, dates, variables=var, time_scale="monthly")

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/daymet_grid.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/daymet.ipynb

If the input geometry (or coordinate) is in a CRS other than EPSG:4326, we should pass
it to the functions.

.. code-block:: python

    coords = (-1431147.7928, 318483.4618)
    crs = "epsg:3542"
    dates = ("2000-01-01", "2006-12-31")
    annual = daymet.get_bycoords(coords, dates, variables=var, loc_crs=crs, time_scale="annual")

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/daymet_loc.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/daymet.ipynb

Also, we can use the ``potential_et`` function to compute PET by passing the daily climate data.
We can either pass a ``pandas.DataFrame`` or a ``xarray.Dataset``. Note that, ``penman_monteith``
and ``priestley_taylor`` methods have parameters that can be passed via the ``params`` argument,
if any value other than the default values are needed. For example, default value of ``alpha``
for ``priestley_taylor`` method is 1.26 (humid regions), we can set it to 1.74 (arid regions)
as follows:

.. code-block:: python

    pet_hs = daymet.potential_et(daily, methods="priestley_taylor", params={"alpha": 1.74})

Next, let's get annual total precipitation for Hawaii and Puerto Rico for 2010.

.. code-block:: python

    hi_ext = (-160.3055, 17.9539, -154.7715, 23.5186)
    pr_ext = (-67.9927, 16.8443, -64.1195, 19.9381)
    hi = daymet.get_bygeom(hi_ext, 2010, variables="prcp", region="hi", time_scale="annual")
    pr = daymet.get_bygeom(pr_ext, 2010, variables="prcp", region="pr", time_scale="annual")

Some example plots are shown below:

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/hi.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/daymet.ipynb

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/pr.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/daymet.ipynb
AsyncRetriever: Asynchronous requests with persistent caching
-------------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/async_retriever.svg
    :target: https://pypi.python.org/pypi/async_retriever
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/async_retriever.svg
    :target: https://anaconda.org/conda-forge/async_retriever
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/async_retriever/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/async_retriever
    :alt: CodeCov

.. image:: https://img.shields.io/pypi/pyversions/async_retriever.svg
    :target: https://pypi.python.org/pypi/async_retriever
    :alt: Python Versions

.. image:: https://github.com/cheginit/async_retriever/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/async_retriever/actions/workflows/test.yml
    :alt: Github Actions

|

.. image:: https://img.shields.io/badge/security-bandit-green.svg
    :target: https://github.com/PyCQA/bandit
    :alt: Security Status

.. image:: https://www.codefactor.io/repository/github/cheginit/async_retriever/badge
   :target: https://www.codefactor.io/repository/github/cheginit/async_retriever
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

|

Features
--------

AsyncRetriever is a part of `HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services. This package has only one purpose;
asynchronously sending requests and retrieving responses as ``text``, ``binary``, or ``json``
objects. It uses persistent caching to speed up the retrieval even further. Moreover, thanks
to `nest_asyncio <https://github.com/erdewit/nest_asyncio>`__ you can use this package in
Jupyter notebooks. Although this package is in the HyRiver software stack, it's
applicable to any HTTP requests.

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/async_retriever/issues>`__.

Installation
------------

You can install ``async_retriever`` using ``pip``:

.. code-block:: console

    $ pip install async_retriever

Alternatively, ``async_retriever`` can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__:

.. code-block:: console

    $ conda install -c conda-forge async_retriever

Quick start
-----------

AsyncRetriever has two public function: ``retrieve`` for sending requests and ``delete_url_cache``
for removing all requests from the cache file that contain a given URL. By default, ``retrieve``
creates and/or uses ``./cache/aiohttp_cache.sqlite`` as the cache that you can customize it
by the ``cache_name`` argument. Also, by default, the cache doesn't have any expiration date and
the ``delete_url_cache`` function should be used if you know that a database on a server was
updated, and you want to retrieve the latest data. Alternatively, you can use the ``expire_after``
argument to set the expiration date for the cache.

As an example for retrieving a ``binary`` response, let's use the DAAC server to get
`NDVI <https://daac.ornl.gov/VEGETATION/guides/US_MODIS_NDVI.html>`_.
The responses can be directly passed to ``xarray.open_mfdataset`` to get the data as
a ``xarray`` Dataset. We can also disable SSL certificate verification by setting
``ssl=False``.

.. code-block:: python

    import io
    import xarray as xr
    import async_retriever as ar
    from datetime import datetime

    west, south, east, north = (-69.77, 45.07, -69.31, 45.45)
    base_url = "https://thredds.daac.ornl.gov/thredds/ncss/ornldaac/1299"
    dates_itr = ((datetime(y, 1, 1), datetime(y, 1, 31)) for y in range(2000, 2005))
    urls, kwds = zip(
        *[
            (
                f"{base_url}/MCD13.A{s.year}.unaccum.nc4",
                {
                    "params": {
                        "var": "NDVI",
                        "north": f"{north}",
                        "west": f"{west}",
                        "east": f"{east}",
                        "south": f"{south}",
                        "disableProjSubset": "on",
                        "horizStride": "1",
                        "time_start": s.strftime("%Y-%m-%dT%H:%M:%SZ"),
                        "time_end": e.strftime("%Y-%m-%dT%H:%M:%SZ"),
                        "timeStride": "1",
                        "addLatLon": "true",
                        "accept": "netcdf",
                    }
                },
            )
            for s, e in dates_itr
        ]
    )
    resp = ar.retrieve(urls, "binary", request_kwds=kwds, max_workers=8, ssl=False)
    data = xr.open_mfdataset(io.BytesIO(r) for r in resp)

We can remove these requests and their responses from the cache like so:

.. code-block:: python

    ar.delete_url_cache(base_url)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/ndvi.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/async.ipunb

For a ``json`` response example, let's get water level recordings of a NOAA's water level station,
8534720 (Atlantic City, NJ), during 2012, using CO-OPS API. Note that this CO-OPS product has a 31-day
limit for a single request, so we have to break the request down accordingly.

.. code-block:: python

    import pandas as pd

    station_id = "8534720"
    start = pd.to_datetime("2012-01-01")
    end = pd.to_datetime("2012-12-31")

    s = start
    dates = []
    for e in pd.date_range(start, end, freq="m"):
        dates.append((s.date(), e.date()))
        s = e + pd.offsets.MonthBegin()

    url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

    urls, kwds = zip(
        *[
            (
                url,
                {
                    "params": {
                        "product": "water_level",
                        "application": "web_services",
                        "begin_date": f'{s.strftime("%Y%m%d")}',
                        "end_date": f'{e.strftime("%Y%m%d")}',
                        "datum": "MSL",
                        "station": f"{station_id}",
                        "time_zone": "GMT",
                        "units": "metric",
                        "format": "json",
                    }
                },
            )
            for s, e in dates
        ]
    )

    resp = ar.retrieve(urls, read="json", request_kwds=kwds, cache_name="~/.cache/async.sqlite")
    wl_list = []
    for rjson in resp:
        wl = pd.DataFrame.from_dict(rjson["data"])
        wl["t"] = pd.to_datetime(wl.t)
        wl = wl.set_index(wl.t).drop(columns="t")
        wl["v"] = pd.to_numeric(wl.v, errors="coerce")
        wl_list.append(wl)
    water_level = pd.concat(wl_list).sort_index()
    water_level.attrs = rjson["metadata"]

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/water_level.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/async.ipunb

Now, let's see an example without any payload or headers. Here's how we can retrieve
harmonic constituents of several NOAA stations from CO-OPS:

.. code-block:: python

    stations = [
        "8410140",
        "8411060",
        "8413320",
        "8418150",
        "8419317",
        "8419870",
        "8443970",
        "8447386",
    ]

    base_url = "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations"
    urls = [f"{base_url}/{i}/harcon.json?units=metric" for i in stations]
    resp = ar.retrieve(urls, "json")

    amp_list = []
    phs_list = []
    for rjson in resp:
        sid = rjson["self"].rsplit("/", 2)[1]
        const = pd.DataFrame.from_dict(rjson["HarmonicConstituents"]).set_index("name")
        amp = const.rename(columns={"amplitude": sid})[sid]
        phase = const.rename(columns={"phase_GMT": sid})[sid]
        amp_list.append(amp)
        phs_list.append(phase)

    amp = pd.concat(amp_list, axis=1)
    phs = pd.concat(phs_list, axis=1)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/tides.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/async.ipunb
Py3DEP: Topographic data through 3DEP
-------------------------------------

.. image:: https://img.shields.io/pypi/v/py3dep.svg
    :target: https://pypi.python.org/pypi/py3dep
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/py3dep.svg
    :target: https://anaconda.org/conda-forge/py3dep
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/py3dep/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/py3dep
    :alt: CodeCov

.. image:: https://img.shields.io/pypi/pyversions/py3dep.svg
    :target: https://pypi.python.org/pypi/py3dep
    :alt: Python Versions

.. image:: https://pepy.tech/badge/py3dep
    :target: https://pepy.tech/project/py3dep
    :alt: Downloads

|

.. image:: https://www.codefactor.io/repository/github/cheginit/py3dep/badge
   :target: https://www.codefactor.io/repository/github/cheginit/py3dep
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

|

Features
--------

Py3DEP is a part of `HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services. This package provides
access to the `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__
database which is a part of the
`National Map services <https://viewer.nationalmap.gov/services/>`__.
The 3DEP service has multi-resolution sources and depending on the user provided resolution,
the data is resampled on the server-side based on all the available data sources. Py3DEP returns
the requests as `xarray <https://xarray.pydata.org/en/stable>`__ dataset. Moreover,
under-the-hood, this package uses ``requests-cache`` for persistent caching that can improve
the performance significantly. The 3DEP web service includes the following layers:

- DEM
- Hillshade Gray
- Aspect Degrees
- Aspect Map
- GreyHillshade Elevation Fill
- Hillshade Multidirectional
- Slope Map
- Hillshade Elevation Tinted
- Height Ellipsoidal
- Contour 25
- Contour Smoothed 25

Moreover, Py3DEP offers some additional utilities:

- ``elevation_bygrid``: For getting elevations of all the grid points in a 2D grid.
- ``elevation_bycoords``: For getting elevation of a list of ``x`` and ``y`` coordinates.
- ``deg2mpm``: For converting slope dataset from degree to meter per meter.

You can also try using PyGeoHydro without installing
it on your system by clicking on the binder badge. A Jupyter Lab
instance with the HyRiver stack pre-installed will be launched in your web browser, and you
can start coding!

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/py3dep/issues>`__.


Installation
------------

You can install Py3DEP using ``pip`` after installing ``libgdal`` on your system
(for example, in Ubuntu run ``sudo apt install libgdal-dev``). Moreover, Py3DEP has an optional
dependency for using persistent caching, ``requests-cache``. We highly recommend installing
this package as it can significantly speed up send/receive queries. You don't have to change
anything in your code, since Py3DEP under-the-hood looks for ``requests-cache`` and if available,
it will automatically use persistent caching:

.. code-block:: console

    $ pip install py3dep

Alternatively, Py3DEP can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__:

.. code-block:: console

    $ conda install -c conda-forge py3dep

Quick start
-----------

You can use Py3DEP using command-line or as a Python library. The command line interface
provides access to two functionality:

- Getting topographic data: You must create a ``geopandas.GeoDataFrame`` that contains
  the geometries of the target locations. This dataframe must have at least three columns:
  ``id``, ``res``, and ``geometry``. The ``id`` column is used as filenames for saving
  the obtained topographic data to a NetCDF (``.nc``) file. The ``res`` column must be
  the target resolution in meter. Then, you must save the dataframe to a file with extensions
  such as ``.shp`` or ``.gpkg`` (whatever that ``geopandas.read_file`` can read).
- Getting elevation: You must create a ``pandas.DataFrame`` that contains coordinates of the
  target locations. This dataframe must have at least two columns: ``x`` and ``y``. The elevations
  are obtained using ``airmap`` service in meters. The data are saved as a ``csv`` file with the
  same filename as the input file with an ``_elevation`` appended, e.g., ``coords_elevation.csv``.

.. code-block:: console

    $ py3dep --help
    Usage: py3dep [OPTIONS] COMMAND [ARGS]...

    Command-line interface for Py3DEP.

    Options:
    -h, --help  Show this message and exit.

    Commands:
    coords    Retrieve topographic data for a list of coordinates.
    geometry  Retrieve topographic data within geometries.

The ``coords`` sub-command is as follows:

.. code-block:: console

    $ py3dep coords -h
    Usage: py3dep coords [OPTIONS] FPATH

    Retrieve topographic data for a list of coordinates.

    FPATH: Path to a csv file with two columns named ``lon`` and ``lat``.

    Examples:
        $ cat coords.csv
        lon,lat
        -122.2493328,37.8122894
        $ py3dep coords coords.csv -q airmap -s topo_dir

    Options:
    -q, --query_source [airmap|tnm|tep]
                                    Source of the elevation data.
    -s, --save_dir PATH             Path to a directory to save the requested
                                    files. Extension for the outputs is either
                                    `.nc` for geometry or `.csv` for coords.

    -h, --help                      Show this message and exit.

And, the ``geometry`` sub-command is as follows:

.. code-block:: console

    $ py3dep geometry -h
    Usage: py3dep geometry [OPTIONS] FPATH

    Retrieve topographic data within geometries.

    FPATH: Path to a shapefile (.shp) or geopackage (.gpkg) file.
    This file must have three columns and contain a ``crs`` attribute:
        - ``id``: Feature identifiers that py3dep uses as the output netcdf/csv filenames.
        - ``res``: Target resolution in meters.
        - ``geometry``: A Polygon or MultiPloygon.

    Examples:
        $ py3dep geometry ny_geom.gpkg -l "Slope Map" -l DEM -s topo_dir

    Options:
    -l, --layers [DEM|Hillshade Gray|Aspect Degrees|Aspect Map|GreyHillshade_elevationFill|Hillshade Multidirectional|Slope Map|Slope Degrees|Hillshade Elevation Tinted|Height Ellipsoidal|Contour 25|Contour Smoothed 25]
                                    Target topographic data layers
    -s, --save_dir PATH             Path to a directory to save the requested
                                    files.Extension for the outputs is either
                                    `.nc` for geometry or `.csv` for coords.

    -h, --help                      Show this message and exit.


Now, let's see how we can use Py3DEP as a library.

Py3DEP accepts `Shapely <https://shapely.readthedocs.io/en/latest/manual.html>`__'s
Polygon or a bounding box (a tuple of length four) as an input geometry.
We can use PyNHD to get a watershed's geometry, then use it to get the DEM and slope
in meters/meters from Py3DEP using ``get_map`` function.

The ``get_map`` has a ``resolution`` argument that sets the target resolution
in meters. Note that the highest available resolution throughout the CONUS is about 10 m,
though higher resolutions are available in limited parts of the US. Note that the input
geometry can be in any valid spatial reference (``geo_crs`` argument). The ``crs`` argument,
however, is limited to ``CRS:84``, ``EPSG:4326``, and ``EPSG:3857`` since 3DEP only supports
these spatial references.

.. code-block:: python

    import py3dep
    from pynhd import NLDI

    geom = NLDI().get_basins("01031500").geometry[0]
    dem = py3dep.get_map("DEM", geom, resolution=30, geo_crs="epsg:4326", crs="epsg:3857")
    slope = py3dep.get_map("Slope Degrees", geom, resolution=30)
    slope = py3dep.deg2mpm(slope)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/dem_slope.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/3dep.ipynb
    :align: center

We can use `rioxarray <https://github.com/corteva/rioxarray>`__ package to save the obtained
dataset as a raster file:

.. code-block:: python

    import rioxarray

    dem.rio.to_raster("dem_01031500.tif")

Moreover, we can get the elevations of set of x- and y- coordinates on a grid. For example,
let's get the minimum temperature data within this watershed from Daymet using PyDaymet then
add the elevation as a new variable to the dataset:

.. code-block:: python

    import pydaymet as daymet
    import xarray as xr
    import numpy as np

    clm = daymet.get_bygeom(geometry, ("2005-01-01", "2005-01-31"), variables="tmin")
    elev = py3dep.elevation_bygrid(clm.x.values, clm.y.values, clm.crs, clm.res[0] * 1000)
    attrs = clm.attrs
    clm = xr.merge([clm, elev])
    clm["elevation"] = clm.elevation.where(~np.isnan(clm.isel(time=0).tmin), drop=True)
    clm.attrs.update(attrs)

Now, let's get street network data using `osmnx <https://github.com/gboeing/osmnx>`__ package
and add elevation data for its nodes using ``elevation_bycoords`` function.

.. code-block:: python

    import osmnx as ox

    G = ox.graph_from_place("Piedmont, California, USA", network_type="drive")
    x, y = nx.get_node_attributes(G, "x").values(), nx.get_node_attributes(G, "y").values()
    elevation = py3dep.elevation_bycoords(zip(x, y), crs="epsg:4326")
    nx.set_node_attributes(G, dict(zip(G.nodes(), elevation)), "elevation")

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/street_elev.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/3dep.ipynb
    :align: center
PyGeoHydro: Retrieve Geospatial Hydrology Data
----------------------------------------------

.. image:: https://img.shields.io/pypi/v/pygeohydro.svg
    :target: https://pypi.python.org/pypi/pygeohydro
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/pygeohydro.svg
    :target: https://anaconda.org/conda-forge/pygeohydro
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/pygeohydro/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pygeohydro
    :alt: CodeCov

.. image:: https://img.shields.io/pypi/pyversions/pygeohydro.svg
    :target: https://pypi.python.org/pypi/pygeohydro
    :alt: Python Versions

.. image:: https://pepy.tech/badge/hydrodata
    :target: https://pepy.tech/project/hydrodata
    :alt: Downloads

|

.. image:: https://www.codefactor.io/repository/github/cheginit/pygeohydro/badge/main
    :target: https://www.codefactor.io/repository/github/cheginit/pygeohydro/overview/main
    :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/HyRiver-examples/main?urlpath=lab/tree/notebooks
    :alt: Binder

|

Features
--------

PyGeoHydro (formerly named `hydrodata <https://pypi.org/project/hydrodata>`__) is a part of
`HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services. This package provides
access to some public web services that offer geospatial hydrology data. It has three
main modules: ``pygeohydro``, ``plot``, and ``helpers``.

The ``pygeohydro`` module can pull data from the following web services:

* `NWIS <https://nwis.waterdata.usgs.gov/nwis>`__ for daily mean streamflow observations
  (returned as a ``pandas.DataFrame`` or ``xarray.Dataset`` with station attributes),
* `Water Quality Portal <https://www.waterqualitydata.us/>`__ for accessing current and
  historical water quality data from more than 1.5 million sites across the US,
* `NID <https://nid.sec.usace.army.mil>`__ for accessing the National Inventory of Dams
  web service,
* `HCDN 2009 <https://www2.usgs.gov/science/cite-view.php?cite=2932>`__ for identifying sites
  where human activity affects the natural flow of the watercourse,
* `NLCD 2019 <https://www.mrlc.gov/>`__ for land cover/land use, imperviousness, imperviousness
  descriptor, and canopy data. You can get data using both geometries and coordinates.
* `SSEBop <https://earlywarning.usgs.gov/ssebop/modis/daily>`__ for daily actual
  evapotranspiration, for both single pixel and gridded data.

Also, it has two other functions:

* ``interactive_map``: Interactive map for exploring NWIS stations within a bounding box.
* ``cover_statistics``: Categorical statistics of land use/land cover data.

The ``plot`` module includes two main functions:

* ``signatures``: Hydrologic signature graphs.
* ``cover_legends``: Official NLCD land cover legends for plotting a land cover dataset.
* ``descriptor_legends``: Color map and legends for plotting an imperviousness descriptor dataset.

The ``helpers`` module includes:

* ``nlcd_helper``: A roughness coefficients lookup table for each land cover and imperviousness
  descriptor type which is useful for overland flow routing among other applications.
* ``nwis_error``: A dataframe for finding information about NWIS requests' errors.

Moreover, requests for additional databases and functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pygeohydro/issues>`__.

You can find some example notebooks `here <https://github.com/cheginit/HyRiver-examples>`__.

You can also try using PyGeoHydro without installing
it on your system by clicking on the binder badge. A Jupyter Lab
instance with the HyRiver stack pre-installed will be launched in your web browser, and you
can start coding!

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pygeohydro/issues>`__.

Installation
------------

You can install PyGeoHydro using ``pip`` after installing ``libgdal`` on your system
(for example, in Ubuntu run ``sudo apt install libgdal-dev``). Moreover, PyGeoHydro has an optional
dependency for using persistent caching, ``requests-cache``. We highly recommend installing
this package as it can significantly speed up send/receive queries. You don't have to change
anything in your code, since PyGeoHydro under-the-hood looks for ``requests-cache`` and
if available, it will automatically use persistent caching:

.. code-block:: console

    $ pip install pygeohydro

Alternatively, PyGeoHydro can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__:

.. code-block:: console

    $ conda install -c conda-forge pygeohydro

Quick start
-----------

We can explore the available NWIS stations within a bounding box using ``interactive_map``
function. It returns an interactive map and by clicking on a station some of the most
important properties of stations are shown.

.. code-block:: python

    import pygeohydro as gh

    bbox = (-69.5, 45, -69, 45.5)
    gh.interactive_map(bbox)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/interactive_map.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/nwis.ipynb
    :alt: Interactive Map

We can select all the stations within this boundary box that have daily mean streamflow data from
``2000-01-01`` to ``2010-12-31``:

.. code-block:: python

    from pygeohydro import NWIS

    nwis = NWIS()
    query = {
        **nwis.query_bybox(bbox),
        "hasDataTypeCd": "dv",
        "outputDataTypeCd": "dv",
    }
    info_box = nwis.get_info(query)
    dates = ("2000-01-01", "2010-12-31")
    stations = info_box[
        (info_box.begin_date <= dates[0]) & (info_box.end_date >= dates[1])
    ].site_no.tolist()

Then, we can get the daily streamflow data in mm/day (by default the values are in cms)
and plot them:

.. code-block:: python

    from pygeohydro import plot

    qobs = nwis.get_streamflow(stations, dates, mmd=True)
    plot.signatures(qobs)

By default, ``get_streamflow`` returns a ``pandas.DataFrame`` that has a ``attrs`` method
containing metadata for all the stations. You can access it like so ``qobs.attrs``.
Moreover, we can get the same data as ``xarray.Dataset`` as follows:

.. code-block:: python

    qobs_ds = nwis.get_streamflow(stations, dates, to_xarray=True)

This ``xarray.Dataset`` has two dimensions: ``time`` and ``station_id``. It has
10 variables including ``discharge`` with two dimensions while other variables
that are station attitudes are one dimensional.

We can also get instantaneous streamflow data using ``get_streamflow``. This method assumes
that the input dates are in UTC time zone and returns the data in UTC time zone as well.

.. code-block:: python

    date = ("2005-01-01 12:00", "2005-01-12 15:00")
    qobs = nwis.get_streamflow("01646500", date, freq="iv")

The ``WaterQuality`` has a number of convenience methods to retrieve data from the
web service. Since there are many parameter combinations that can be
used to retrieve data, a general method is also provided to retrieve data from
any of the valid endpoints. You can use ``get_json`` to retrieve stations info
as a ``geopandas.GeoDataFrame`` or ``get_csv`` to retrieve stations data as a
``pandas.DataFrame``. You can construct a dictionary of the parameters and pass
it to one of these functions. For more information on the parameters, please
consult the `Water Quality Data documentation <https://www.waterqualitydata.us/webservices_documentation>`__.
For example, let's find all the stations within a bounding box that have Caffeine data:

.. code-block:: python

    from pynhd import WaterQuality

    bbox = (-92.8, 44.2, -88.9, 46.0)
    kwds = {"characteristicName": "Caffeine"}
    wq = WaterQuality()
    stations = wq.station_bybbox(bbox, kwds)

Or the same criterion but within a 30-mile radius of a point:

.. code-block:: python

    stations = wq.station_bydistance(-92.8, 44.2, 30, kwds)

Then we can get the data for all these stations the data like this:

.. code-block:: python

    sids = stations.MonitoringLocationIdentifier.tolist()
    caff = wq.data_bystation(sids, kwds)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/water_quality.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/water_quality.ipynb
    :alt: Water Quality

Moreover, we can get land use/land cove data using ``nlcd_bygeom`` or ``nlcd_bycoods`` functions
and percentages of land cover types using ``cover_statistics``.
The ``nlcd_bycoords`` function returns a ``geopandas.GeoDataFrame`` with the NLCD
layers as columns and input coordinates as the ``geometry`` column. Moreover, The ``nlcd_bygeom``
function accepts both a single geometry or a ``geopandas.GeoDataFrame`` as the input.

.. code-block:: python

    from pynhd import NLDI

    basins = NLDI().get_basins(["01031450", "01031500", "01031510"])
    lulc = gh.nlcd_bygeom(basins, 100, years={"cover": [2016, 2019]})
    stats = gh.cover_statistics(lulc["01031450"].cover_2016)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/lulc.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/nlcd.ipynb
    :alt: Land Use/Land Cover

Next, let's use ``ssebopeta_bygeom`` to get actual ET data for a basin. Note that there's a
``ssebopeta_bycoords`` function that returns an ETA time series for a single coordinate.

.. code-block:: python

    geometry = NLDI().get_basins("01315500").geometry[0]
    eta = gh.ssebopeta_bygeom(geometry, dates=("2005-10-01", "2005-10-05"))

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/eta.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/ssebop.ipynb
    :alt: Actual ET

Additionally, we can pull all the US dams data using ``NID``. Let's get dams that are within this
bounding box and have a maximum storage larger than 200 acre-feet.

.. code-block:: python

    nid = NID()
    dams = nid.get_bygeom((-65.77, 43.07, -69.31, 45.45), "epsg:4326")
    dams = nid.inventory_byid(dams.id.to_list())
    dams = dams[dams.maxStorage > 200]

We can get also all dams within CONUS in NID with maximum storage larger than 200 acre-feet:

.. code-block:: python

    import geopandas as gpd

    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    conus = world[world.name == "United States of America"].geometry.iloc[0].geoms[0]

    dam_list = nid.get_byfilter([{"maxStorage": ["[200 5000]"]}])
    dams = dam_list[0][dam_list[0].is_valid]
    dams = dams[dams.within(conus)]

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/dams.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/nid.ipynb
    :alt: Dams
