<!---
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
--->

Atlite's contributor guidelines can be found in the official [documentation](https://atlite.readthedocs.io/en/master/contributing.html).<!--
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
-->

Closes # (if applicable).

## Change proposed in this Pull Request

<!--- Provide a general, short summary of your changes in the title above -->

## Description
<!--- Describe your changes in detail -->

## Motivation and Context
<!--- Why is this change required? What problem does it solve? -->
<!--- If it fixes an open issue, please link to the issue here. -->

## How Has This Been Tested?
<!--- Please describe in detail how you tested your changes. -->
<!--- Include details of your testing environment, and the tests you ran to -->
<!--- see how your change affects other areas of the code, etc. -->

## Type of change
<!--- What types of changes does your code introduce? Put an `x` in all the boxes that apply: -->
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to change)

## Checklist
<!--- Go over all the following points, and put an `x` in all the boxes that apply. -->
<!--- If you're unsure about any of these, don't hesitate to ask. We're here to help! -->
- [ ] I tested my contribution locally and it seems to work fine.
- [ ] I locally ran `pytest` inside the repository and no unexpected problems came up.
- [ ] I have adjusted the docstrings in the code appropriately.
- [ ] I have documented the effects of my code changes in the documentation `doc/`.
- [ ] I have added newly introduced dependencies to `environment.yaml` file.
- [ ] I have added a note to release notes `doc/release_notes.rst`.
- [ ] I have used `pre-commit run --all` to lint/format/check my contribution
---
name: Bug report
about: Create a report if something doesn't work quite right.
title: ''
labels: bug
assignees: ''
---
<!---
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
--->


<!-- Provide a general summary of the issue -->

## Description
<!-- Provide a more detailed introduction to the issue itself, and why you consider it to be a bug -->
<!-- If you can, add a minimal example which reproduces the bug -->

## Expected Behavior
<!-- Tell us what should happen -->

## Actual Behavior
<!-- Tell us what goes wrong and happens instead -->

## Error Message
<!-- Paste any terminal output and error message you encounter here to help illustrate the problem -->

## Your Environment
<!-- Include relevant details about the environment you experienced the bug in -->
* The `atlite` version used:
* How you installed `atlite` (`conda`, `pip` or `github`):
* Operating System:
* My environment:
    <details>
      <summary>(output of `conda list`)</summary>
      ```
        <!-- output of `conda list` -->
      ```
    </details>


---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---
<!---
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
--->

<!-- Provide a general summary of the feature you would like to see -->

## Detailed Description
<!-- Provide a detailed description of the change or addition you are proposing -->

## Context
<!-- Why is this change important to you? -->
<!-- How would you use it? -->

## Possible Implementation
<!-- Not obligatory, but suggest an idea for implementing addition or change -->
..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

#############
Release Notes
#############


.. Upcoming Release
.. =================

* Atlite now supports calculating dynamic line ratings based on the IEEE-738 standard (https://github.com/PyPSA/atlite/pull/189).
* The wind feature provided by ERA5 now also calculates the wind angle `wnd_azimuth` in range [0 - 2π) spanning the cirlce from north in clock-wise direction (0 is north, π/2 is east, -π is south, 3π/2 is west).
* A new intersection matrix function was added, which works similarly to incidence matrix but has boolean values.
* Atlite now supports two CSP (concentrated solar power) technologies, solar tower and parabolic trough. See (https://atlite.readthedocs.io/en/latest/examples/working-with-csp.html) for details.
* The solar position (azimuth and altitude) are now part of the cutout feature `influx`. Cutouts created with earlier versions will become incompatible with the next major.
* Automated upload of code coverage reports via Codecov.
* DataArrays returned by `.pv(...)` and `.wind(...)` now have a clearer name and 'units' attribute.
* If the `matrix` argument in conversion functions (`.pv(...)`, `.wind(...)` etc.) is a `DataArray`, the alignment of the coordinate axis with the cutout grid is double-checked. 
* Due to ambiguity, conversion functions (`.pv(...)`, `.wind(...)` etc.) now raise an `ValueError` if shapes and matrix are given. 
* Atlite now supports calculating of heat pump coefficients of performance (https://github.com/PyPSA/atlite/pull/145).
* Enabled the GitHub feature "Cite this repository" to generate a BibTeX file (Added a `CITATION.cff` file to the repository).

**Bug fixes**
* The solar position for ERA5 cutouts is now calculated for half a time step earlier (time-shift by `cutout.dt/2`) to account for the aggregated nature of
  ERA5 variables (see https://github.com/PyPSA/atlite/issues/158). The fix is only applied to newly created cutouts. Previously created cutouts do not profit
  from this fix and need to be recreated `cutout.prepare(overwrite=True)`.


Version 0.2.5 
==============

* Clarification for ``ExclusionContainer.add_raster(..)`` that ``codes=..`` does not accept ``lambda``-functions in combination with ``multiprocessing``.
* Internal change: We are moving to `black` for internal code formatting.
* Fix ignored keywords in convert_and_aggregate(...) for capacity_layout=True.

Version 0.2.4 
==============

* Fix cutout merge and update for xarray ``>=v0.18.0`` (https://github.com/PyPSA/atlite/issues/147)
* Set multiprocessing context to ``spawn`` for ensuring equal computation across all platforms. 

Version 0.2.3 
==============

* The progressbar used in ``atlite.gis.availability_matrix`` is now a `tqdm` progressbar which displays better in parallel executions.
* The function ``layout_from_capacity_list`` was added to the cutout class. It is a convenience function that calculates the aggregated capacities per cutout grid cells (layout) based on a list of capacities with coordinates, e.g. list of wind turbines.    
* The dask version was fixed to a xarray-compatible versions (see https://github.com/dask/dask/issues/7583)

Version 0.2.2 
==============

This update is mainly due to fixes in the data handling of the SARAH module. If you work with the SARAH data, we encourage you to update. 

* Fixed compatibility with xarray v0.17.
* Fixed sarah data for ``dx = dy = 0.05``. Due to the float32 dtype of the sarah coordinates, the cutout coordinates were corrupted when merging. This was fixed in the sarah module by converting the coordinates to float64. This also speeds up the cutout creation for more coarse grained cutouts.  
* Fixed sarah data for a time frequency of 30 minutes. This was raising an assertion error as the (new) pandas frequency string for 30 minutes is '30T' not '30min'.
* Fix the ``regrid`` function in ``atlite.gis`` for target coords which are not having the same bounds as the original ``xarray.Dataset``. The previous implementation was leading to a small shift of coordinates in the preparation of SARAH data.



Version 0.2.1
==============
* The `regrid` function in `atlite.gis` was fixed. The previous implementation set an affine transform starting at the center of a cell at the origin. The corrected transform starts at the real origin (origin of the origin cell). Further a padding of the extent ensures that all values are taken into account in the target projection.  
* Exclusion Calculation is now possible with `atlite` (find an usage example at Examples -> Calculate Landuse Availability), Therefore 

  - a new class  `atlite.gis.ExclusionContainer`  was added. It serves as a container of rasters and geometries which should be excluded from the landuse availability.  
  - `Cutout` has a new `availabilitymatrix` function which calculates the overlap of weather cells with shapes while excluding areas based on an `ExclusionContainer`.  
  - `Cutout` has now a affine transform property (`rasterio.Affine`). 
* Fix resolution for dx and dy unequal to 0.25: Due to floating point precision errors, loading data with ERA5 corrupted the cutout coordinates. This was fixed by converting the dtype of era5 coordinates to float64 and rounding. Corresponding tests were added.
* Round cutout.dx and cutout.dy in order to prevent precision errors.    
* Allow passing keyword arguments to `dask.compute` in `convert_and_aggregate` functions. 
* The Cutout class has a new property `bounds` (same as extent but in different order).

**Breaking Change**
* `Cutout.extent` was adjusted to cover the whole cutout area. The extent is now a numpy array. Before, it indicated the coordinates of the centers of the corner cells. 

Version 0.2
===============

**Major changes**


* Atlite now **requires Python 3.6 or higher**.
* We changed the Atlite backend for storing cutout data.
  Existing cutouts either need to be migrated with the
  appropriate functions or (what we recommended) recreated.
* The backend change also includes some changes to the API.
  Most notably:
  
  - The `xarray` for cutouts is now exposed as `Cutout.data`
  - The `Cutout.meta` attribute was deprecated in favour of
    `Cutout.data.attrs`
  - `xarray` and `dask` can now handle some data caching
    automatically.
    If you wish to preload some data before your calculation,
    you can now use `Cutout.data.load()` to load all of the
    cutouts data into memory.  
    *(Warning: Requires a large enough memory.)*
  - The `Cutout` class has a new property `grid`, a GeoPandas DataFrame 
    which combines and deprecates `grid_cells()` and `grid_coordinates()`
* The order of coordinates (indices) for `Cutouts` changed: `x` and `y` (e.g. longitude and latitude) are now both ascending (before: `x` ascending and `y` descending).
* Following the lead of geopandas, pyproj, cartopy and rasterio, atlite now uses Coordinate Reference System (`CRS`) instead of the old   fashioned projection strings. 

**New features**


* You can now use wind turbine configurations as stored in the
  `Open Energy Database <https://openenergy-platform.org/dataedit/view/supply/turbine_library>`_
  using the string prefix `"oedb:"` when specifying a turbine,
  e.g. `"oedb:Enercon_E-141/4200"`.
* Atlite now has and uses a new configuration system.
  See the new section on `configuration <https://atlite.readthedocs.io/en/latest/configuration.html>`_
  for details.
* It is possible to merge two cutouts together, using `Cutout.merge`


**Breaking changes**

* The argument `show_progress` of function `atlite.convert.convert_and_aggregate` does not take strings anymore. 
* The argument `layout` of function `atlite.convert.convert_and_aggregate` must be a `xarray.DataArray`.
* Due to the change of the order of coordinates in cutouts the order of coordinates in `matrix` passed to `convert_*` functions
    changed likewise: `x` and `y` are both ascending now.
* Due to the change of the order of coordinates in cutouts the order of elements returned by `grid_coordinates()` has changed.
* Due to the change of the order of coordinates in cutouts the order of elements in the attribute `grid_cells` has changed.


Version 0.0.4
===============

* support negative latitudes to PV panel orientation
* add support for ERA5 back extension to 1950
* add PROJ>=7 valid 'aea' projection string 



Version 0.0.3
==============

Brings a minor bug fix and prepares for the next version jump to version 0.2.

* Fix heat demand hourshift for xarray 0.15.1
* Add Travis CI and simplified release management..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0


=========================
Authors and Contributions
=========================

..
  Use this marker to reference the table of authors from other files.

.. headline-marker

Major contributions until now from:

+--------------------+----------------------+----------------------+
| Contributions      | Contributor          | Affiliation          |
+====================+======================+======================+
| 2016-2020          | Jonas Hörsch         | * FIAS Frankfurt     |
|                    |                      | * KIT Karlsruhe      |
|                    |                      | * RLI Berlin         |
+--------------------+----------------------+----------------------+
| 2019-2021          | Fabian Hofmann       | FIAS Frankfurt       |
+--------------------+----------------------+----------------------+
| 2019-2020          | Johannes Hampp       | University Giessen   |
+--------------------+----------------------+----------------------+
| 2016-2019          | Tom Brown            | * FIAS Frankfurt     |
|                    |                      | * KIT Karlsruhe      |
+--------------------+----------------------+----------------------+
| 2016-2017          | Gorm Andresen        | Aarhus University    |
+--------------------+----------------------+----------------------+
| 2016-2017          | David Schlachtberger | FIAS Frankfurt       |
+--------------------+----------------------+----------------------+
| 2016-2017          | Markus Schlott       | FIAS Frankfurt       |
+--------------------+----------------------+----------------------+
  .. SPDX-FileCopyrightText: 2016-2021 The Atlite Authors

  .. SPDX-License-Identifier: CC-BY-4.0

======
Atlite
======

|PyPI version| |Conda version| |Documentation Status| |ci| |codecov| |standard-readme compliant| |GPL-3-or-later-image| |reuse| |black| |pre-commit.ci| |joss|

Atlite is a `free software`_, `xarray`_-based Python library for
converting weather data (like wind speeds, solar influx) into energy systems data.
It has a  lightweight design and works with big weather datasets
while keeping the resource requirements especially on CPU and RAM
resources low.


.. Atlite is designed to be modular, so that it can work with any weather
.. datasets. It currently has modules for the following datasets: 

.. * `NCEP Climate Forecast System <http://rda.ucar.edu/datasets/ds094.1/>`_ hourly
..   historical reanalysis weather data available on a 0.2 x 0.2 degree global grid
.. * `ECMWF ERA5
..   <https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation>`_ hourly
..   historical reanalysis weather data on an approximately 0.25 x 0.25 deg global
..   grid
.. * `EURO-CORDEX Climate Change Projection <http://www.euro-cordex.net/>`_
..   three-hourly up until 2100, available on a 0.11 x 0.11 degree grid for Europe
.. * `CMSAF SARAH-2
..   <https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>`_
..   half-hourly historical surface radiation on a 0.05 x 0.05 deg grid available
..   for Europe and Africa (automatically interpolated to a 0.2 deg grid and
..   combined with ERA5 temperature).


Atlite can process the following weather data fields and can convert them into following power-system relevant time series for any subsets of a full weather database.

.. image:: doc/workflow_chart.png

.. * Temperature
.. * Downward short-wave radiation
.. * Upward short-wave radiation
.. * Wind 
.. * Runoff
.. * Surface roughness
.. * Height maps
.. * Soil temperature


.. * Wind power generation for a given turbine type
.. * Solar PV power generation for a given panel type
.. * Solar thermal collector heat output
.. * Hydroelectric inflow (simplified)
.. * Heating demand (based on the degree-day approximation)


Atlite was initially developed by the `Renewable Energy Group
<https://fias.uni-frankfurt.de/physics/schramm/renewable-energy-system-and-network-analysis/>`_
at `FIAS <https://fias.uni-frankfurt.de/>`_ to carry out simulations
for the `CoNDyNet project <http://condynet.de/>`_, financed by the
`German Federal Ministry for Education and Research (BMBF)
<https://www.bmbf.de/en/index.html>`_ as part of the `Stromnetze
Research Initiative
<http://forschung-stromnetze.info/projekte/grundlagen-und-konzepte-fuer-effiziente-dezentrale-stromnetze/>`_.


Installation
============

To install you need a working installation running Python 3.6 or above
and we strongly recommend using either miniconda or anaconda for package
management.

To install the current stable version:

with ``conda`` from `conda-forge`_

.. code:: shell

       conda install -c conda-forge atlite

with ``pip`` from `pypi`_

.. code:: shell

       pip install atlite

to install the most recent upstream version from `GitHub`_

.. code:: shell

       pip install git+https://github.com/pypsa/atlite.git


Documentation
===============
.. * Install atlite from conda-forge or pypi.
.. * Download one of the weather datasets listed above (ERA5 is downloaded
..   automatically on-demand after the ECMWF
..   `cdsapi<https://cds.climate.copernicus.eu/api-how-to>` client is 
..   properly installed)
.. * Create a cutout, i.e. a geographical rectangle and a selection of
..   times, e.g. all hours in 2011 and 2012, to narrow down the scope -
..   see `examples/create_cutout.py <examples/create_cutout.py>`_
.. * Select a sparse matrix of the geographical points inside the cutout
..   you want to aggregate for your time series, and pass it to the
..   appropriate converter function - see `examples/ <examples/>`_


Please check the `documentation <https://atlite.readthedocs.io/en/latest>`_.

Contributing
============

If you have any ideas, suggestions or encounter problems, feel invited
to file issues or make pull requests.

Authors and Copyright
---------------------

Copyright (C) 2016-2021 The Atlite Authors.

See the `AUTHORS`_ for details.

Licence
=======

|GPL-3-or-later-image|

This work is licensed under multiple licences:

-  All original source code is licensed under `GPL-3.0-or-later`_.
-  Auxiliary code from SPHINX is licensed under `BSD-2-Clause`_.
-  The documentation is licensed under `CC-BY-4.0`_.
-  Configuration and data files are mostly licensed under `CC0-1.0`_.

See the individual files for license details.

.. _free software: http://www.gnu.org/philosophy/free-sw.en.html
.. _xarray: http://xarray.pydata.org/en/stable/

.. _conda-forge: https://anaconda.org/conda-forge/atlite
.. _pypi: https://pypi.org/project/atlite/%3E
.. _GitHub: https://github.com/pypsa/atlite

.. _documentation on getting started: https://atlite.readthedocs.io/en/latest/getting-started.html

.. _AUTHORS: AUTHORS.rst

.. _GPL-3.0-or-later: LICENSES/GPL-3.0-or-later.txt
.. _BSD-2-Clause: LICENSES/BSD-2-Clause.txt
.. _CC-BY-4.0: LICENSES/CC-BY-4.0.txt
.. _CC0-1.0: LICENSES/CC0-1.0.txt

.. |PyPI version| image:: https://img.shields.io/pypi/v/atlite.svg
   :target: https://pypi.python.org/pypi/atlite
.. |Conda version| image:: https://img.shields.io/conda/vn/conda-forge/atlite.svg
   :target: https://anaconda.org/conda-forge/atlite
.. |Documentation Status| image:: https://readthedocs.org/projects/atlite/badge/?version=master
   :target: https://atlite.readthedocs.io/en/master/?badge=master
.. |standard-readme compliant| image:: https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat
   :target: https://github.com/RichardLitt/standard-readme
.. |GPL-3-or-later-image| image:: https://img.shields.io/pypi/l/atlite.svg
   :target: LICENSES/GPL-3.0-or-later.txt
.. |codecov| image:: https://codecov.io/gh/PyPSA/atlite/branch/master/graph/badge.svg?token=TEJ16CMIHJ
   :target: https://codecov.io/gh/PyPSA/atlite
.. |ci| image:: https://github.com/PyPSA/atlite/actions/workflows/CI.yaml/badge.svg
   :target: https://github.com/PyPSA/atlite/actions/workflows/CI.yaml
.. |reuse| image:: https://api.reuse.software/badge/github.com/pypsa/atlite
   :target: https://api.reuse.software/info/github.com/pypsa/atlite
.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Code style: black
.. |pre-commit.ci| image:: https://results.pre-commit.ci/badge/github/PyPSA/atlite/master.svg
   :target: https://results.pre-commit.ci/latest/github/PyPSA/atlite/master
   :alt: pre-commit.ci status
.. |joss| image:: https://joss.theoj.org/papers/10.21105/joss.03294/status.svg
   :target: https://doi.org/10.21105/joss.03294..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0


============
Contributing
============

We welcome anyone interested in contributing to this project,
be it with new ideas, suggestions, by filing bug reports or
contributing code.

You are invited to submit pull requests / issues to our
`Github repository <https://github.com/pypsa/atlite>`_.

For linting, formatting and checking your code contributions
against our guidelines (e.g. we use `Black <https://github.com/psf/black>`_ as code style
and aim for `REUSE compliance <https://reuse.software/>`_,
use `pre-commit <https://pre-commit.com/index.html>`_:

1. Installation ``conda install -c conda-forge pre-commit`` or ``pip install pre-commit``
2. Usage:
    * To automatically activate ``pre-commit`` on every ``git commit``: Run ``pre-commit install``
    * To manually run it: ``pre-commit run --all``

Contributing examples
=====================

Nice examples are always welcome.

You can even submit your `Jupyter notebook`_ (``.ipynb``) directly
as an example.
For contributing notebooks (and working with notebooks in `git`
in general) we have compiled a workflow for you which we suggest
you follow:

* Locally install `this precommit hook for git`_

This obviously has to be done only once.
The hook checks if any of the notebooks you are including in a commit
contain a non-empty output cells.

Then for every notebook:

1. Write the notebook (let's call it ``foo.ipynb``) and place it
   in ``examples/foo.ipynb``.
2. Ask yourself: Is the output in each of the notebook's cells
   relevant for to example?

    * Yes: Leave it there.
      Just make sure to keep the amount of pictures/... to a minimum.
    * No: Clear the output of all cells,
      e.g. `Edit -> Clear all output` in JupyterLab.

3. Provide a link to the documentation:
   Include a file ``foo.nblink`` located in ``doc/examples/foo.nblink``
   with this content

   .. code-block:
        {
            'path' : '../../examples/foo.ipynb'
        }

    Adjust the path for your file's name.
    This ``nblink`` allows us to link your notebook into the documentation.
4. Link your file in the documentation:

   Either

    * Include your ``examples/foo.nblink`` directly into one of
      the documentations toctrees; or
    * Tell us where in the documentation you want your example to show up

5. Commit your changes.
   If the precommit hook you installed above kicks in, confirm
   your decision ('y') or go back ('n') and delete the output
   of the notebook's cells.
6. Create a pull request for us to accept your example.

For the reasoning and details behind this workflow, see `atlite/issues#38`_.

The support for the the ``.ipynb`` notebook format in our documentation
is realised via the extensions `nbsphinx`_ and `nbsphinx_link`_.

.. _Jupyter notebook: https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/what_is_jupyter.html
.. _this precommit hook for git: https://jamesfolberth.org/articles/2017/08/07/git-commit-hook-for-jupyter-notebooks/
.. _atlite/issues#38: https://github.com/PyPSA/atlite/issues/38
.. _nbsphinx: https://nbsphinx.readthedocs.io/en/0.4.2/installation.html
.. _nbsphinx_link: https://nbsphinx.readthedocs.io/en/latest/

Authors, Contributions and Copyright
------------------------------------

.. include:: ../AUTHORS.rst
    :start-after: headline-marker
..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

.. include:: ../RELEASE_NOTES.rst..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

#############
API Reference
#############

Cutout
------

.. autoclass:: atlite.Cutout
    :members:


Data
------

.. automodule:: atlite.data
    :members:


Convert
-------

.. automodule:: atlite.convert
    :members:


Resource
--------

.. automodule:: atlite.resource
    :members:

Wind
----

.. automodule:: atlite.wind
    :members:

GIS
-------

.. automodule:: atlite.gis
    :members:


Utils
-----

.. automodule:: atlite.utils
    :members:..
  SPDX-FileCopyrightText: 2016-2021 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

############
Conventions
############

Atlite uses the following conventions which are applied for processing geo-spatial and temporal data. 

Grid coordinates 
================

According to the ``xarray`` conventions, grid coordinates are ordered such that both ``x`` and ``y`` are ascending.  

The coordinates represent points in the **center** of the corresponding grid cells. Given a cutout, the geographical data of the grid is given by :py:attr:`atlite.Cutout.grid`, which returns a `GeoPandas` dataframe with coordinates and geometries. The coordinates are the geographical centroids of the geometries, the grid cells. When initializing a cutout with e.g. 

>>> cutout = atlite.Cutout('example', module='era5', x=slice(5,10), y=slice(30,35), time='2013')

the cutout is built with centroids starting at :math:`5^\circ` for ``x`` (longitude) and :math:`30^\circ` for ``y`` (latitude). That means the effective covered area by the cutout spans from longitude :math:`4.875^\circ` to :math:`10.125^\circ` and from latitude :math:`29.875^\circ` to :math:`35.125^\circ`, given the default resolution of :math:`0.25^\circ\times0.25^\circ` per grid cell. This information is accessible via the the extent of the cutout:


>>> cutout.extent 
array([ 4.875, 10.125, 29.875, 35.125])


Time Points 
===========

Following the ERA5 convention, the time-index of time-dependent data refers to the end of the time-span over which was averaged. So, given a time resolution of 1 hour (=averaging window), the time-index 12:00 refers to the time-averaged values from 11:00 to 12:00. For datasets other than ERA5 this convention is not necessarily fulfilled. For example, the SARAH data refers to instantaneous data-points, i.e. data for the time-index 12:00 refers to the momentaneous value of the variable at 12:00. In the implementation, we try to consider this circumstance in order to appropriately align the datasets in order to merge them.  
..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

############
Installation
############

There are three possibilities to install atlite:

* conda::

    conda install -c conda-forge atlite


* pypi::

    pip install atlite

* or directly from GitHub for the most recent version::

    pip install git+https://github.com/pypsa/atlite/

Requirements
============

Conda environment
-----------------

We provide a `conda environment file <https://github.com/PyPSA/atlite/blob/documentation/environment.yaml>`_
for conveniently setting up all required and optional packages.
For information on setting up conda environments from file,
`click here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_.

Python version
--------------

With the new version 0.2, Atlite requires an installed version of
**Python 3.6 or above**.

Required packages
-----------------

* bottleneck
* cdsapi
* dask (0.18.0+)
* geopandas
* netcdf4
* numexpr
* numpy
* pandas
* progressbar2
* rasterio
* rtree
* scipy
* shapely
* xarray (0.11.2+)


Computational resources
-----------------------

As for requirements on your computing equipment, we tried to keep
the resource requirements low.
The requirements obviously depend on the size of the cutouts and
datasets your parse and use.

We run our conversions on our laptops and this usually works fine
and can run in the background without clogging our computers.

With regards to

* CPU: While Atlite does some number crunching, it does not require
  special or large multicore CPUs
* Memory: For the ERA5 dataset you should be fine running Atlite with
  even 2-4 GiB.
  Other datasets can require more memory, as they sometimes need to be
  loaded fully or partially into memory for creating a cutout.
* Disk space: Is really all about the cutout and dataset sizes
  (in time and space) you use.
  We can only provide two examples for reference:
  
    - Small cutout (Republic of Ireland + UK + some atlantic ocean),
      1 month with hourly resolution using ERA5: 60 MiB
    - Large cutout (Western Asia), 
      1 year with hourly resolution using ERA5: 6 GiB

We guess you do not need to worry about your computer being able to handle
common small or medium scenarios.

Rule of thumb:
    The requirements on your machine are increasing with the
    size of the cutout (in time and space).
..
  SPDX-FileCopyrightText: 2016-2021 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

Atlite: Convert weather data to energy systems data
===================================================

.. image:: https://img.shields.io/pypi/v/atlite.svg
    :target: https://pypi.python.org/pypi/atlite
    :alt: PyPI version

.. image:: https://img.shields.io/conda/vn/conda-forge/atlite.svg
    :target: https://anaconda.org/conda-forge/atlite
    :alt: Conda version

.. image:: https://github.com/PyPSA/atlite/actions/workflows/CI.yaml/badge.svg
   :target: https://github.com/PyPSA/atlite/actions/workflows/CI.yaml

.. image:: https://codecov.io/gh/PyPSA/atlite/branch/master/graph/badge.svg?token=TEJ16CMIHJ
   :target: https://codecov.io/gh/PyPSA/atlite

.. image:: https://readthedocs.org/projects/atlite/badge/?version=master
    :target: https://atlite.readthedocs.io/en/latest/?badge=master
    :alt: Documentation Status

.. image:: https://img.shields.io/pypi/l/atlite.svg
    :target: License

.. image:: https://api.reuse.software/badge/github.com/pypsa/atlite
    :target: https://api.reuse.software/info/github.com/pypsa/atlite

.. image:: https://joss.theoj.org/papers/10.21105/joss.03294/status.svg
    :target: https://doi.org/10.21105/joss.03294

Atlite is a `free software
<http://www.gnu.org/philosophy/free-sw.en.html>`_, `xarray
<http://xarray.pydata.org/en/stable/>`_-based Python library for
converting weather data (such as wind speeds, solar radiation,
temperature and runoff) into power systems data (such as wind
power, solar power, hydro power and heating demand time series).
It is designed to work with big datasets as is common with weather
reanalysis, while maintaining low computational requirements.

The spatial and time resolution of the obtainable power time series
depends on the resolutions of the original weather reanalysis dataset.
E.g. using our recommended dataset `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_
we can obtain time-series with hourly resolution on a 30 km x 30 km
grid.

The time series derived with Atlite are used in energy system models
like e.g. `PyPSA-EUR <https://github.com/PyPSA/pypsa-eur>`_
or projects like `model.energy <https://model.energy/>`_.

Maintainers
===========

Atlite is currently maintained by volunteers from different institutions
with no dedicated funding for developing this package.

Atlite was initially developed by the `Renewable Energy Group
<https://fias.uni-frankfurt.de/physics/schramm/renewable-energy-system-and-network-analysis/>`_
at `FIAS <https://fias.uni-frankfurt.de/>`_ to carry out simulations
for the `CoNDyNet project <http://condynet.de/>`_, financed by the
`German Federal Ministry for Education and Research (BMBF)
<https://www.bmbf.de/en/index.html>`_ as part of the `Stromnetze
Research Initiative
<http://forschung-stromnetze.info/projekte/grundlagen-und-konzepte-fuer-effiziente-dezentrale-stromnetze/>`_.

Origin
======

Atlite was originally conceived as a light-weight version of the Aarhus
University RE Atlas (`original publication doi:10.1016/j.energy.2015.09.071 <http://dx.doi.org/10.1016/j.energy.2015.09.071>`_).
It has since been extended to use weather datasets simulated with projected
climate change and to compute other time series, such as hydro power,
solar thermal collectors and heating demand.


Citing Atlite
=============

If you would like to cite the Atlite software, please refer to `this paper <https://doi.org/10.21105/joss.03294>`_ published in `JOSS <https://joss.theoj.org/>`_.



.. Documentation
.. =============

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Getting Started

   introduction
   installation
   conventions

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Examples

   examples/create_cutout.ipynb
   examples/create_cutout_SARAH.ipynb
   examples/historic-comparison-germany.ipynb
   examples/landuse-availability.ipynb
   examples/using_gebco_heightmap.ipynb
   examples/plotting_with_atlite.ipynb
   examples/logfiles_and_messages.ipynb
   examples/more_examples.rst

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: References

   ref_api
   release_notes
   contributing
   license

License
=======

Atlite is released and licensed under the 
`GPLv3 <http://www.gnu.org/licenses/gpl-3.0.en.html>`_.
..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

############
Introduction
############


Atlite processes weather data and converts it into energy system 
relevant quantities, mainly you can create data for 


* **Wind** power generation: Using predefined or custom turbine properties
  and smoothing options for modelling more realistic results.
  *New:* Turbines can also be imported from the 
  `Open Energy Database <https://openenergy-platform.org/dataedit/view/supply/wind_turbine_library>`_.
* **Solar (PV)** power generation: Using predefined or custom panel properties.
* **Solar (thermal)** heat generation from solar collectors.
* **Hydro (run-off)** power generation.
* **Heating demand** (based on degree-day approx.).

How it works
========================

Starting from a global weather dataset, *e.g.* our standard data source ECMWF's ERA5 dataset,  

.. image:: img/worldmap.png
    :align: center
    :alt: alternate text


Atlite enables you to create a **cutout**, a spatial and 
temporal subset of the original data which includes all relevant data 
such as wind velocity, influx, temperature etc.

.. image:: img/cutout.png
    :align: center
    :alt: alternate text

The **cutout** consists of grid cells of a certain resolution (depending on your input data, here 30 km x 30 km) 
each having its own weather timeseries. 
The **cutout** builds the starting point to your calculations. 

From there, you can extract various quantities and general properties of the data, *e.g.* wind capacity factors per grid cell 

.. image:: img/capfactors.png
    :width: 400pt
    :align: center
    :alt: alternate text

for a specific turbine type (this gives you information on the share of capacity which is in average running and producing power). 

Further, you can set power plants to specific spots and let Atlite calculate their actual power production. The **capacity layout**  specifies which grid cell contains what amount of capacity of *one* production type. Atlite comes along with a small `library of wind turbine configurations and PV panel configurations <https://github.com/PyPSA/atlite/tree/master/atlite/resources>`_  which you can directly use, *e.g.* 'Vestas_V112_3MW' wind turbines.

.. image:: img/layout.png
    :width: 400pt
    :align: center
    :alt: alternate text

Atlite then calculates the power generation data for each cell and either aggregates them to buses

.. image:: img/produced_power.png
    :width: 400pt
    :align: center
    :alt: alternate text

or to geometrical shapes


.. image:: img/production_per_country.png
    :align: center


Whereas for the first case, grid cells must directly be assigned to buses by passing a matrix of size :math:`N_{cell} \times N_{bus}`, for the second case, the aggregation to shapes takes place in Atlite itself: It creates the mentioned matrix, the so-called **indicator matrix**, which contains the spatial overlap of each grid cell (weighted by the capacity layout if present) with each shape. This is why the shapes can the very refined and even smaller than the grid cells. In our example the **indicator matrix** for the shape of United Kingdom without being weighted by the **capacity layout** looks like this


.. image:: img/indicator_matrix.png
    :align: center
    :width: 400pt



Datasets
==================

The standard data source we currently employ is ECMWF's ERA5 dataset
(reanalysis weather data in a ca. 30 km x 30 km and hourly resolution).
This dataset is easily available at no additional costs and requires only
minimal setup from the user in comparison to other datasets.
It is downloaded automatically on-demand after the 
`ECMWF ADS API <https://cds.climate.copernicus.eu/api-how-to>`_
(European Centre for Medium-Range Weather Forecasts Climate Data Store
Application Program Interface) client is properly installed. See separate,
linked installation guide for details, especially for correctly setting up
your CDS API key.

Previously and in the future other datasets where and (hopefully) will 
again be usable, including

* the *NCEP Climate Forecast System* dataset 
* the *EURO-CORDEX Climate Change Projection* dataset
* the *CMSAF SARAH-2* dataset
* Satellite based radiation observations, e.g. SARAH-2.
* Weather data forecasts from climate models.

Their support however is currently on hold (time limitation on developer
side).

If you need to process these (or other) data sources, feel free to
file an issue on our `GitHub <https://github.com/PyPSA/atlite>`_ or (even better) create a pull request!




What Atlite does not cover (yet)
=================================

* Atlite does not provide and **graphical user interface** (GUI) and relies on prior knowledge on working with Python commands.

* Atlite does not provide **exact prediction** of the time-series generation at high resolution in a **future point** in time. The spatial resolution of the  results is limited by the input data used. The accuracy of the results is in parts limited by the methodologies used for translating weather data into generation and the underlying assumptions. With the current assumptions Atlite is not suited for predicting the output of single wind turbines or solar panels.

* As the results of Atlite are theoretical and are not validated per se, and while usually a good approximation, can **deviate significantly from reality**. While in the past and also at the moment datasets generate by packages similar to Atlite where commonly used without a comparison and validation with reality, there is currently a trend to validate the datasets before using them to make sure that results are at least plausible. The Atlite team is planning to include auxiliary functions which help to validate generated datasets.
..
  SPDX-FileCopyrightText: 2016-2019 The Atlite Authors

  SPDX-License-Identifier: CC-BY-4.0

#############
More Examples
#############

All examples are available as Jupyter notebooks in our 
`Atlite GitHub repository <https://github.com/pypsa/atlite/docs/examples>`_.


If you are missing an example for a certain
aspect of atlite or if you have created some
examples of your own we would be happy to hear
from you.
Feel free to open a 
`Pull Request <https://github.com/PyPSA/atlite/pulls>`_
or an
`Issue <https://github.com/PyPSA/atlite/issues>`_
for your suggestion.

