# In general
Contributions are very welcome. Please make sure there is a github issue
associated with with every pull request. Creating an issue is also a good
way to propose new features.

# Readthedocs
For technical assistance with your contribution, please check the [contributing
guidelines on
readthedocs](https://era5cli.readthedocs.io/en/latest/contribute.html).

# Making a release

## Author information
Ensure all authors are present in:

- `.zenodo.json`
- `CITATION.cff`
- `era5cli/__version__.py`

## Confirm release info
Ensure the right date and upcoming version number (including release candidate
tag, if applicable) is set in:

- `CITATION.cff`
- `era5cli/__version__.py`

## Update the changelog
Update `CHANGELOG.rst` with new features and fixes in the upcoming version.
Confirm that `README.rst` is up to date with new features as well.

## Release on GitHub
Open [releases](https://github.com/eWaterCycle/era5cli/releases) and draft a new
release. Copy the changelog for this version into the description (though note
that the description is in Markdown, so reformat from Rst if necessary).

Tag the release according to [semantic versioning
guidelines](https://semver.org/), preceded with a `v` (e.g.: v1.0.0). The
release title is the tag and the release date together (e.g.: v1.0.0
(2019-07-25)).

### Release candidate
When releasing a release candidate on Github, tick the pre-release box, and
amend the version tag with `-rc` and the candidate number. Ensure the release
candidate version is accurate in `CITATION.cff` and `era5cli/__version__.py`.
If the version number in these files is not updated, Zenodo release
workflows will fail.

Releasing a release candidate is not required, but can help detect bugs early.

## PyPI release workflow
Publishing a new release in github triggers the github Action workflow that
builds and publishes the package to test.PyPI or PyPI. Versions with "rc"
(release candidate) in their version tag will only be published to test.PyPI.
Other version tags will trigger a PyPI release. Inspect
`.github/workflows/publish-to-pypi.yml` for more information.

Confirm a release candidate on test.PyPI with:
```
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ era5cli
```

## Release on Zenodo
Confirm the new release on [Zenodo](https://doi.org/10.5281/zenodo.3252665).

## Release on the Research Software Directory
Wait a few hours, then confirm the addition of a new release on the
[RSD](https://www.research-software.nl/software/era5cli).
Release Notes
*************

1.3.1 (2021-12-01)
~~~~~~~~~~~~~~~~~~
* Automatic Zenodo/RSD release failed; updated contribution guidelines `#106 <https://github.com/eWaterCycle/era5cli/pull/106>`_

1.3.0 (2021-11-30)
~~~~~~~~~~~~~~~~~~
* Fix compatibility with changed CDS variables geopotential/orography `#98 <https://github.com/eWaterCycle/era5cli/pull/98>`_
* Add integration testing `#102 <https://github.com/eWaterCycle/era5cli/pull/102>`_

1.2.1 (2021-04-21)
~~~~~~~~~~~~~~~~~~
* Automatic PyPI release for 1.2.0 failed; updated github action workflow `#91 <https://github.com/eWaterCycle/era5cli/pull/91>`_

1.2.0 (2021-04-21)
~~~~~~~~~~~~~~~~~~
* Add support for ERA5-Land data `#67 <https://github.com/eWaterCycle/era5cli/pull/67>`_
* Add functionality to download subregions `#70 <https://github.com/eWaterCycle/era5cli/pull/70>`_
* Update variables available for ERA5 datasets `#84 <https://github.com/eWaterCycle/era5cli/pull/84>`_

1.1.1 (2020-12-15)
~~~~~~~~~~~~~~~~~~
* Patch to fix the github actions publish automation `#64 <https://github.com/eWaterCycle/era5cli/pull/64>`_

1.1.0 (2020-12-14)
~~~~~~~~~~~~~~~~~~
* The stable 1.1.0 era5cli minor release.
* Add support for ERA5 preliminary back extension `#58 <https://github.com/eWaterCycle/era5cli/pull/58>`_
* Update tests `#57 <https://github.com/eWaterCycle/era5cli/pull/57>`_
* Add automated PyPI package building and publishing with github Actions `#62 <https://github.com/eWaterCycle/era5cli/pull/62>`_

1.0.0 (2019-07-25)
~~~~~~~~~~~~~~~~~~
* The stable 1.0.0 era5cli release.
* Adding more useful information to netCDF history `#48 <https://github.com/eWaterCycle/era5cli/pull/48>`_

1.0.0rc3 (2019-07-16)
~~~~~~~~~~~~~~~~~~~~~
* Third Release Candidate for the stable 1.0.0 era5cli release.
* Improve documentation `#21 <https://github.com/eWaterCycle/era5cli/issues/21>`_ `#29 <https://github.com/eWaterCycle/era5cli/issues/29>`_.
* Cleanup command line options `#19 <https://github.com/eWaterCycle/era5cli/issues/19>`_ `#20 <https://github.com/eWaterCycle/era5cli/issues/20>`_.
* Append era5cli version to history of downloaded netCDF file `#17 <https://github.com/eWaterCycle/era5cli/issues/17>`_.

1.0.0rc2 (2019-07-01)
~~~~~~~~~~~~~~~~~~~~~
* Second Release Candidate for the stable 1.0.0 era5cli release.
* Fix downloading all variables when requesting multiple variables and using --split `#23 <https://github.com/eWaterCycle/era5cli/issues/23>`_.
* Fix link to PyPI package in documentation `#22 <https://github.com/eWaterCycle/era5cli/issues/22>`_.

1.0.0rc1 (2019-06-22)
~~~~~~~~~~~~~~~~~~~~~
* First Release Candidate for the stable 1.0.0 era5cli release: A commandline utility to download ERA-5 data using cdsapi.
era5cli
=======
.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://opensource.org/licenses/Apache-2.0
    :alt: License

.. image:: https://img.shields.io/badge/docs-stable-brightgreen.svg
   :target: http://era5cli.readthedocs.io/en/stable/?badge=stable
   :alt: Stable Documentation

.. image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :target: http://era5cli.readthedocs.io/en/latest/?badge=latest
   :alt: Latest Documentation

.. image:: https://github.com/eWaterCycle/era5cli/actions/workflows/test_codecov.yml/badge.svg
   :target: https://github.com/eWaterCycle/era5cli/actions/workflows/test_codecov.yml
   :alt: Github Actions

.. image:: https://codecov.io/gh/eWaterCycle/era5cli/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/eWaterCycle/era5cli
   :alt: Test coverage

.. image:: https://badge.fury.io/py/era5cli.svg
    :target: https://badge.fury.io/py/era5cli

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3252665.svg
   :target: https://doi.org/10.5281/zenodo.3252665

.. image:: https://img.shields.io/badge/rsd-era5cli-00a3e3.svg
   :target: https://www.research-software.nl/software/era5cli
   :alt: Research Software Directory Badge

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow
   :target: https://fair-software.eu
   
.. inclusion-marker-start-do-not-remove

A command line interface to download ERA5 data from the `Copernicus Climate Data Store <https://climate.copernicus.eu/>`_.

With era5cli you can:

- download meteorological data in GRIB/NetCDF, including ERA5 data from the preliminary back extension, and ERA5-Land data.
- list and retrieve information on available variables and pressure levels
- select multiple variables for several months and years
- split outputs by years, producing a separate file for every year instead of merging them in one file
- download multiple files at once
- extract data for a sub-region of the globe

.. inclusion-marker-end-do-not-remove

| Free software: Apache Software License 2.0
| Documentation: https://era5cli.readthedocs.io
.. include:: ../../CHANGELOG.rst
Contribute
**********

Bug reports
===========

File bug reports or feature requests, and make contributions (e.g. code
patches), by `opening a new issue on GitHub <https://github.com/ewatercycle/era5cli/issues>`_.

Please give as much information as you can in the issue. It is very useful if
you can supply the command or a small self-contained code snippet that
reproduces the problem. Please also specify the era5cli version that you are
using and the version of the cdsapi library installed.

Contribute to the tool
======================

Create and activate the development environment:
::

    python3 -m venv env
    . env/bin/activate


Populate the development environment with the required dependencies:
::

    pip install -U pip
    pip install -r requirements-dev.txt
    pip install -e .

Before pushing a new addition, run flake8 and pytest to confirm that the code
is up to standard.

Use flake8 to check for code style issues:
::

   flake8 era5cli/

Use pytest to run the test suite:
::

   pytest era5cli/

Deactivate the environment with:
::

   deactivate


Contribute to the documentation
===============================

When updating the documentation, use the environment created above.

Build the documentation with:
::

   sphinx-build docs/source docs/buildInstructions
------------

Register at Copernicus Climate Data Service
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  You have to register at Copernicus Climate Data Service:
   `copernicus <https://cds.climate.copernicus.eu/user/register?destination=%2F%23!%2Fhome>`__.
   After activating your account use your new account to log in. In you
   profile page you can find your user ID and your API key.

-  Copy your user ID and API key.

Create a key ascii file
~~~~~~~~~~~~~~~~~~~~~~~

Linux
#####
In Linux create a new file called .cdsapirc in the home directory of your user and add the following two lines:

::

   url: https://cds.climate.copernicus.eu/api/v2

   key: UID:KEY

Replace UID with your user ID and KEY with your API key

Windows
#######
In Windows create a new file called .cdsapirc (e.g. with Notepad) where in your windows environment, %USERPROFILE% is usually located at C:\Users\Username folder). And add the following two lines

::

   url: https://cds.climate.copernicus.eu/api/v2

   key: UID:KEY

Replace UID with your user ID and KEY with your API key

MacOS
#####
In MacOS create a new file called .cdsapirc in the home directory of your user and add the following two lines:


::

   url: https://cds.climate.copernicus.eu/api/v2

   key: UID:KEY

Replace UID with your user ID and KEY with your API key

Info on available variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :module: era5cli.cli
   :func: _build_parser
   :prog: era5cli
   :path: info

Tip: search in variables
########################

To quickly search the list of variables for a specific word, you can use the
built-in ``grep`` command. For example:

::

   era5cli info "2Dvars" | grep temperature


should list all single level variables that contain the word 'temperature', so
they can be easily identified for an era5cli request.

Running era5cli from the command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
era5cli can be used to fetch both hourly data and monthly averaged data.

Fetching hourly data
####################

Fetch hourly data through an cdsapi call via command line. More information on the available data and options can be found on:

| `Era5 hourly single levels download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels>`_.
| `Era5 hourly pressure levels download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels>`_.
| `Era5 hourly single levels preliminary back extension download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-preliminary-back-extension>`_.
| `Era5 hourly pressure levels preliminary back extension download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-preliminary-back-extension>`_.
| `Era5-Land hourly download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land>`_.


.. argparse::
   :module: era5cli.cli
   :func: _build_parser
   :prog: era5cli
   :path: hourly


Fetching monthly data
#####################

Fetch monthly data through an cdsapi call via command line. More information on the available data and options can be found on:

| `Era5 monthly single levels download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means>`_.
| `Era5 monthly pressure levels download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means>`_.
| `Era5 monthly single levels preliminary back extension download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means-preliminary-back-extension>`_.
| `Era5 monthly pressure levels preliminary back extension download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension>`_.
| `Era5-Land monthly download page <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means>`_.


For the monthly data, some of the variables are not available. Exceptions on the single level data can be found in table 8 of
`ERA5 parameter listings <https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Parameterlistings>`_

.. argparse::
   :module: era5cli.cli
   :func: _build_parser
   :prog: era5cli
   :path: monthly


Removing or canceling requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERA-5 download requests will be saved in the `Your requests <https://cds.climate.copernicus.eu/cdsapp#!/yourrequests>`_ section in your profile on the Copernicus Climate Data Store. Here you can re-download the requested data, cancel active requests, or remove old requests.

Note that it is currently not possible to cancel active requests from the command line: Killing the process will not download the data to your local machine but still add it to your Copernicus account.
Installation
------------

Install from PyPI
~~~~~~~~~~~~~~~~~
era5cli is available as a package in `PyPI <https://pypi.org/project/era5cli/>`_.

To install from PyPI, the following command can be used:
::

   pip install era5cli

Install from GitHub
~~~~~~~~~~~~~~~~~~~
To install directly from GitHub, the following command can be used:
::

   pip install -U  git+https://github.com/eWaterCycle/era5cli.git
Welcome to era5cli's documentation!
===================================
.. include:: ../../README.rst
  :start-after: inclusion-marker-start-do-not-remove
  :end-before: inclusion-marker-end-do-not-remove

Content
~~~~~~~

.. toctree::
   :maxdepth: 2

   installation
   instructions
   api


Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Meta Information
~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   contribute
   changes
   license
.. _api:

API
============

.. contents::
    :local:
    :depth: 1

era5cli.info
~~~~~~~~~~~~
.. automodule:: info
    :members:
    :undoc-members:
    :show-inheritance:

era5cli.fetch
~~~~~~~~~~~~~
.. automodule:: fetch
    :members:
    :undoc-members:
    :show-inheritance:
