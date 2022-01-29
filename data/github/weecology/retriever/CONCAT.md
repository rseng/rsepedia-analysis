![Retriever logo](http://i.imgur.com/se7TtrK.png)


[![Python package](https://github.com/weecology/retriever/actions/workflows/python-package.yml/badge.svg)](https://github.com/weecology/retriever/actions/workflows/python-package.yml)
[![Build Status (windows)](https://ci.appveyor.com/api/projects/status/qetgo4jxa5769qtb/branch/main?svg=true)](https://ci.appveyor.com/project/ethanwhite/retriever/branch/main)
[![Research software impact](http://depsy.org/api/package/pypi/retriever/badge.svg)](http://depsy.org/package/python/retriever)
[![codecov.io](https://codecov.io/github/weecology/retriever/coverage.svg?branch=main)](https://codecov.io/github/weecology/retriever?branch=main)
[![Documentation Status](https://readthedocs.org/projects/retriever/badge/?version=latest)](http://retriever.readthedocs.io/en/latest/?badge=latest)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/weecology/retriever/main/LICENSE)
[![Join the chat at https://gitter.im/weecology/retriever](https://badges.gitter.im/weecology/retriever.svg)](https://gitter.im/weecology/retriever?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1038272.svg)](https://doi.org/10.5281/zenodo.1038272)
[![JOSS Publication](http://joss.theoj.org/papers/10.21105/joss.00451/status.svg)](https://doi.org/10.21105/joss.00451)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/retriever/badges/downloads.svg)](https://anaconda.org/conda-forge/retriever)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/retriever/badges/version.svg)](https://anaconda.org/conda-forge/retriever)
[![Version](https://img.shields.io/pypi/v/retriever.svg)](https://pypi.python.org/pypi/retriever)
<a href="https://numfocus.org/sponsored-projects">
<img alt="NumFOCUS"
   src="https://i0.wp.com/numfocus.org/wp-content/uploads/2019/06/AffiliatedProject.png" width="100" height="18">
</a>

Finding data is one thing. Getting it ready for analysis is another. Acquiring,
cleaning, standardizing and importing publicly available data is time consuming
because many datasets lack machine readable metadata and do not conform to
established data structures and formats. The Data Retriever automates the first
steps in the data analysis pipeline by downloading, cleaning, and standardizing
datasets, and importing them into relational databases, flat files, or
programming languages. The automation of this process reduces the time for a
user to get most large datasets up and running by hours, and in some cases days.

## Installing the Current Release

If you have Python installed you can install the current release using either `pip`:

```bash
pip install retriever
```

or `conda` after adding the `conda-forge` channel (`conda config --add channels conda-forge`):

```bash
conda install retriever
```

Depending on your system configuration this may require `sudo` for `pip`:

```bash
sudo pip install retriever
```

Precompiled binary installers are also available for Windows, OS X, and
Ubuntu/Debian on
the [releases page](https://github.com/weecology/retriever/releases). These do
not require a Python installation.

[List of Available Datasets](https://retriever.readthedocs.io/en/latest/datasets_list.html)
----------------------------

Installing From Source
----------------------

To install the Data Retriever from source, you'll need Python 3.6.8+ with the following packages installed:

* xlrd

The following packages are optionally needed to interact with associated
database management systems:

* PyMySQL (for MySQL)
* sqlite3 (for SQLite)
* psycopg2-binary (for PostgreSQL), previously psycopg2.
* pyodbc (for MS Access - this option is only available on Windows)
* Microsoft Access Driver (ODBC for windows)

### To install from source

Either use `pip` to install directly from GitHub:

```shell
pip install git+https://git@github.com/weecology/retriever.git
```

or:

1. Clone the repository
2. From the directory containing setup.py, run the following command: `pip
   install .`. You may need to include `sudo` at the beginning of the
   command depending on your system (i.e., `sudo pip install .`).

More extensive documentation for those that are interested in developing can be found [here](http://retriever.readthedocs.io/en/latest/?badge=latest)

Using the Command Line
----------------------
After installing, run `retriever update` to download all of the available dataset scripts.
To see the full list of command line options and datasets run `retriever --help`.
The output will look like this:

```shell
usage: retriever [-h] [-v] [-q]
                 {download,install,defaults,update,new,new_json,edit_json,delete_json,ls,citation,reset,help}
                 ...

positional arguments:
  {download,install,defaults,update,new,new_json,edit_json,delete_json,ls,citation,reset,help}
                        sub-command help
    download            download raw data files for a dataset
    install             download and install dataset
    defaults            displays default options
    update              download updated versions of scripts
    new                 create a new sample retriever script
    new_json            CLI to create retriever datapackage.json script
    edit_json           CLI to edit retriever datapackage.json script
    delete_json         CLI to remove retriever datapackage.json script
    ls                  display a list all available dataset scripts
    citation            view citation
    reset               reset retriever: removes configuration settings,
                        scripts, and cached data
    help

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -q, --quiet           suppress command-line output
```

To install datasets, use `retriever install`:

```shell
usage: retriever install [-h] [--compile] [--debug]
                         {mysql,postgres,sqlite,msaccess,csv,json,xml} ...

positional arguments:
  {mysql,postgres,sqlite,msaccess,csv,json,xml}
                        engine-specific help
    mysql               MySQL
    postgres            PostgreSQL
    sqlite              SQLite
    msaccess            Microsoft Access
    csv                 CSV
    json                JSON
    xml                 XML

optional arguments:
  -h, --help            show this help message and exit
  --compile             force re-compile of script before downloading
  --debug               run in debug mode
```


### Examples

These examples are using the [*Iris* flower dataset](https://en.wikipedia.org/wiki/Iris_flower_data_set).
More examples can be found in the Data Retriever documentation.

Using Install

```shell
retriever install -h   (gives install options)
```

Using specific database engine, retriever install {Engine}

```shell
retriever install mysql -h     (gives install mysql options)
retriever install mysql --user myuser --password ******** --host localhost --port 8888 --database_name testdbase iris
```
install data into an sqlite database named iris.db you would use:

```shell
retriever install sqlite iris -f iris.db
```

Using download

```shell
retriever download -h    (gives you help options)
retriever download iris
retriever download iris --path C:\Users\Documents
```

Using citation

```shell
retriever citation   (citation of the retriever engine)
retriever citation iris  (citation for the iris data)
```

Spatial Dataset Installation
----------------------------

**Set up Spatial support**

To set up spatial support for Postgres using Postgis please
refer to the [spatial set-up docs](https://retriever.readthedocs.io/en/latest/spatial_dbms.html).

```shell
retriever install postgres harvard-forest # Vector data
retriever install postgres bioclim # Raster data
# Install only the data of USGS elevation in the given extent
retriever install postgres usgs-elevation -b -94.98704597353938 39.027001800158615 -94.3599408119917 40.69577051867074

```

Website
-------

For more information see the
[Data Retriever website](http://www.data-retriever.org/).

Acknowledgments
---------------

Development of this software was funded by the [Gordon and Betty Moore
Foundation's Data-Driven Discovery
Initiative](https://www.moore.org/initiative-strategy-detail?initiativeId=data-driven-discovery) through
[Grant GBMF4563](http://www.moore.org/grants/list/GBMF4563) to Ethan White and
the [National Science Foundation](http://nsf.gov/) as part of a [CAREER award to
Ethan White](http://nsf.gov/awardsearch/showAward.do?AwardNumber=0953694).
# Contributing to the Data Retriever

We welcome contributions of all kinds including improvements to the core code,
addition of new dataset scripts to add new datasets to the Retriever,
improvements to the documentation, bug reports, or anything else you can think
of. We strive to be supportive of anyone who wants to contribute, so don't be
shy, give it a go, and we'll do our best to help.  One way to ease into
contributing is to
[add datasets](https://retriever.readthedocs.io/en/latest/recipes.html) to the
Retriever.

## Process for contributing changes

We use a standard
[GitHub flow](https://guides.github.com/introduction/flow/index.html) for
development and reviewing contributions. Fork the repository. Make changes to a
branch of your fork and then submit a pull request.


## Running the tests

We use [pytest](https://docs.pytest.org/en/latest/) for testing. To run the
tests first install nose using pip:

`pip install pytest`

Then from the root of the repository install the Retriever:

`python setup.py install`

and run the tests:

`pytest`

You should see a bunch of output from the Retriever showing the test results.

Tests for MySQL and PostgreSQL require properly configured database management
systems for testing.

### Postgres setup

Requires that the `postgres` user has permissions on a database named `testdb_retriever`
from `localhost`. This login information should be stored in the [postgreSQL
password file](http://www.postgresql.org/docs/9.1/static/libpq-pgpass.html).


### MySQL setup

Requires that the `travis` user has permissions on a database named `testdb_retriever`
from `localhost`.


## Continuous integration

We use GitHub actions and appveyor for continuous integration
testing. All pull requests will automatically report whether the tests are
passing.

# v3.0.0

### Major changes

- Add provenance support to the Data Retriever
- Use utf-8 as default
- Move scripts from Retriever to retriever-recipes repository
- Adapt google code style and add linters, use yapf. Test linters
- Extend CSV field size limit
- Improve output when connection is not made
- Add __version__ to the interface
- Prompt user if a newer version of script is available
- Add all the recipes datasets
- Add test for installation of committed dataset
- Add function to commit dataset

### Minor changes

- Improve "argcomplete-command"
- Add NUMFOCUS logo in README

# v2.4.0

### Minor changes

- Update long description
- Remove Python 2 utilities

### New datasets

- Catalogos-dados-brasil
- Transparencia-dados-abertos-brasil
- biotimesql

### Dataset changes

- Elton-traits and mt-st-helens-veg change encoding

### v2.3.1

### Minor changes

- Update PyPi description

# v2.3.0

### Major changes

- Change Psycopg2 to psycopg2-binary
- Add Spatial data testing on Docker
- Add option for pretty json
- keep order of fetched tables and order of processing resources
- Add reset to specific dataset and script function
- Use tqdm 4.30.0
- Install data into custom director using data_dir option
- Download data into custom directory using sub_dir

### Minor changes

- Add tests for reset script
- Add smaller samples of GIS data for testing
- Reactivate MySQL tests on Travis
- Allow custom arguments for psql
- Add docs and examples for Postgis support
- Change testdb name to testdb_retriever
- Improve Pypi retriever description
- Update documentation for passwordless setup of Postgres on Windows
- Setting up infrastructure for automating script creation

### New datasets

- USA eco legions, `ecoregions-us`
- LTREB `Prairie-forest` ecotone of eastern Kansas/Foster Lab dataset
- Sonoran Desert, `sonoran-desert`
- Adding Acton Lake dataset `acton-lake`

### Dataset changes

- MammalSuperTree.py to `mammal_super_tree.py`
- lakecats_finaltables.json to `lakecats_final_tables`
- harvard_forests.json to `harvard_forest.json`
- macroalgal_communities to `macroalgal-communities`

# v2.2.0

## Major changes

- Use requests package to for handling downloads
- Add support for spatial datasets using postGIS and PostgreSQL
- Update ls to include more details about datasets
- Update license lookup for datasets
- Update keyword lookup for datasets
- Use tqdm for cleaner progress tracking
- Add `fetch` function that installs a dataset and returns it as a dictionary of dataframes

## Minor changes

- Documentation refinement
- Connect to MySQL using preferred encoding.
- License search and keyword search added.
- Conda_Forge docs
- Add Zenodo badge to link to archive
- Add test for extracting data
- Changed all "-" in JSON files to "_"

## New datasets

- Add Noaa Fisheries trade, noaa-fisheries-trade.
- Add Fishery Statistical Collections data, fao-global-capture-product.
- Add bupa liver disorders dataset, bupa-liver-disorders.
- Add GLOBI interactions data. globi-interaction.
- Addition of the National Aquatic Resource Surveys (NARS), nla.
- Addition of partners in flight dataset, partners-in-flight.
- Add the ND-GAIN Country Index. nd-gain.
- Add world GDP in current US Dollars. dgp.
- Add airports dataset, airports.
- Repair aquatic animal excretion.
- Add Biotime dataset. 
- Add lakecats final tables dataset, lakecats-final-tables.
- Add harvard forests data, harvard forests.
- Add USGS elevation data, usgs-elevation.

# v2.1.0

## Major changes
- Add Python interface
- Add Retriever to conda
- Auto complete of Retriever commands on Unix systems

## Minor changes

- Add license to datasets
- Change the structure of raw data from string to list
- Add testing on any modified dataset
- Improve memory usage in cross-tab processing
- Add capabilitiy for datasets to use custom Encoding
- Use new Python interface for regression testing
- Use Frictionless Data specification terminology for internals

## New datasets
- Add ant dataset and weather data to the portal dataset
- NYC TreesCount
- PREDICTS
- aquatic_animal_excretion
- biodiversity_response
- bird_migration_data
- chytr_disease_distr
- croche_vegetation_data
- dicerandra_frutescens
- flensburg_food_web
- great_basin_mammal_abundance
- macroalgal_communities
- macrocystis_variation
- marine_recruitment_data
- mediter_basin_plant_traits
- nematode_traits
- ngreatplains-flowering-dates
- portal-dev
- portal
- predator_prey_body_ratio
- predicts
- socean_diet_data
- species_exctinction_rates
- streamflow_conditions
- tree_canopy_geometries
- turtle_offspring_nesting
- Add vertnet individual datasets
  vertnet_amphibians
  vertnet_birds
  vertnet_fishes
  vertnet_mammals
  vertnet_reptiles

# v2.0.0

## Major changes

* Add Python 3 support, python 2/3 compatibility
* Add json and xml as output formats
* Switch to using the frictionless data datapackage json standard. This a
  **backwards incompatible change** as the form of dataset description files the
  retriever uses to describe the location and processing of simple datasets has
  changed.
* Add CLI for creating, editing, deleting datapackage.json scripts
* Broaden scope to include non-ecological data and rename to Data Retriever
* Major expansion of documentation and move documentation to Read the Docs
* Add developer documentation
* Remove the GUI
* Use csv module for reading of raw data to improve handling of newlines in fields
* Major expansion of integration testing
* Refactor regression testing to produce a single hash for a dataset regardless
  of output format
* Add continuous integration testing for Windows


## Minor changes

* Use pyinstaller for creating exe for windows and app for mac and remove py2app
* Use 3 level semantic versioning for both scripts and core code
* Rename datasets with more descriptive names
* Add a retriever minimum version for each dataset
* Rename dataset description files to follow python modules conventions
* Switch to py.test from nose
* Expand unit testing
* Add version requirements for sqlite and postgresql
* Default to latin encoding
* Improve UI for updating user on downloading and processing progress


## New datasets

* Added machine Learning datasets from UC Irvine's machine learning data sets

# v1.8.3

* Fixed regression in GUI

# v1.8.2

* Improved cleaning of column names
* Fixed thread bug causing Gentry dataset to hang when installed via GUI
* Removed support for 32-bit only Macs in binaries
* Removed unused code

# v1.8.0

* Added scripts for 21 new datasets: leaf herbivory, biomass allocation,
  community dynamics of shortgrass steppe plants, mammal and bird foraging
  attributes, tree demography in Indian, small mammal community dynamics in
  Chile, community dynamics of Sonoran Desert perennials, biovolumes of
  freshwater phytoplankton, plant dynamics in Montana, Antarctic Site Inventory
  breeding bird survey, community abundance data compiled from the literature,
  spatio-temporal population data for butterflies, fish parasite host ecological
  characteristics, eBird, Global Wood Density Database, multiscale community
  data on vascular plants in a North Carolina, vertebrate home range sizes,
  PRISM climate data, Amniote life history database, woody plan Biomass And
  Allometry Database, Vertnet data on amphibians, birds, fishes, mammals,
  reptiles
* Added `reset` command to allow resetting database configuration settings,
  scripts, and cached raw data
* Added Dockerfile for building docker containers of each version of the
  software for reproducibility
* Added support for wxPython 3.0
* Added support for `tar` and `gz` archives
* Added support for archive files whose contents don't fit in memory
* Added checks for and use of system proxies
* Added ability to download archives from web services
* Added tests for regressions in download engine
* Added `citation` command to provide information on citing datasets
* Improved column name cleanup
* Improved whitespace consistency
* Improved handling of Excel files
* Improved function documentation
* Improved unit testing and added coverage analysis
* Improved the sample script by adding a url field
* Improved script loading behavior by only loading a script the first time it is
  discovered
* Improved operating system identification
* Improved download engine by allowing ability to maintain archive and
  subdirectory structure (particular relevant for spatial data)
* Improved cross-platform directory and line ending handling
* Improved testing across platforms
* Improved checking for updated scripts so that scripts are only downloaded if
  the current version isn't available
* Improved metadata in setup.py
* Fixed type issues in Portal dataset
* Fixed GUI always downloading scripts instead of checking if it needed to
* Fixed bug that sometimes resulted in `.retriever` directories not belonging to
  the user who did the installation
* Fixed issues with downloading files to specific paths
* Fixed BBS50 script to match newer structure of the data
* Fixed bug where csv files were not being closed after installation
* Fixed errors when closing the GUI
* Fixed issue where enclosing quotes in csv files were not being respected
  during cross-tab restructuring
* Fixed bug causing v1.6 to break when newer scripts were added to `version.txt`
* Fixed Bioclim script to include `hdr` files
* Fixed missing icon images on Windows
* Removed unused code

# v1.7.0

* Added ability to download files directly for non-tabular data
* Added scripts to download Bioclim and Mammal Supertree data
* Added a script for the MammalDIET database
* Fixed bug where some nationally standardized FIA surveys where not included
* Added check for wxpython on installation to allow non-gui installs
* Fixed several minor issues with Gentry script including a missing site and a column in one file that was misnamed
* Windows install now adds the retriever to the path to facilitate command line use
* Fixed a bug preventing installation from PyPI
* Added icons to installers
* Fixed the retriever failing when given a script it couldn't handle

# v1.6.0

* Added full OS X support to the Retriever
* Added a proper Windows installer
* Fixed a number of bugs
---
title: 'Retriever: Data Retrieval Tool'
tags:
  - data retrieval, data processing, python, data, data science, datasets
authors:
 - name: Henry Senyondo
   orcid: 0000-0001-7105-5808
   affiliation: 1
 - name: Benjamin D. Morris
   orcid: 0000-0003-4418-1360
 - name: Akash Goel
   orcid: 0000-0001-9878-0401
   affiliation: 3
 - name: Andrew Zhang
   orcid: 0000-0003-1148-7734
   affiliation: 4
 - name: Akshay Narasimha
   orcid: 0000-0002-3901-2610
   affiliation: 5
 - name: Shivam Negi
   orcid: 0000-0002-5637-0479
   affiliation: 6
 - name: David J. Harris
   orcid: 0000-0003-3332-9307
   affiliation: 4
 - name: Deborah Gertrude Digges
   orcid: 0000-0002-7840-5054
   affiliation: 10
 - name: Kapil Kumar
   orcid: 0000-0002-2292-1033
   affiliation: 7
 - name: Amritanshu Jain
   orcid: 0000-0003-1187-7900
   affiliation: 5
 - name: Kunal Pal
   orcid: 0000-0002-9657-0053
   affiliation: 8
 - name: Kevinkumar Amipara
   orcid: 0000-0001-5021-2018
   affiliation: 9
 - name: Ethan P. White
   orcid: 0000-0001-6728-7745
   affiliation: 1, 2
affiliations:
 - name: Department of Wildlife Ecology and Conservation, University of Florida
   index: 1
 - name: Informatics Institute, University of Florida
   index: 2
 - name: Delhi Technological University, Delhi
   index: 3
 - name: The University of Florida
   index: 4
 - name: Birla Institute of Technology and Science, Pilani
   index: 5
 - name: Manipal Institute of Technology, Manipal
   index: 6
 - name: National Institute of Technology, Delhi
   index: 7
 - name: RWTH Aachen University, Aachen, Germany
   index: 8
 - name: Sardar Vallabhbhai National Institute of Technology, Surat
   index: 9
 - name: PES Institute of Technology, Bengaluru
   index: 10


date: 16 September 2017 
bibliography: paper.bib
---

# Summary

The Data Retriever automates the first steps in the data analysis workflow by downloading, cleaning, and standardizing tabular datasets, and importing them into relational databases, flat files, or programming languages [@Morris2013]. The automation of this process reduces the time for a user to get most large datasets up and running by hours to days. The retriever uses a plugin infrastructure for both datasets and storage backends. New datasets that are relatively well structured can be added adding a JSON file following the Frictionless Data tabular data metadata package standard [@frictionlessdata_specs]. More complex datasets can be added using a Python script to handle complex data cleaning and merging tasks and then defining the metadata associated with the cleaned tables. New storage backends can be added by overloading a general class with details for storing the data in new file formats or database management systems. The retriever has both a Python API and a command line interface. An R package that wraps the command line interface and a Julia package that wraps the Python API are also available.

The 2.0 and 2.1 releases add extensive new functionality. This includes the Python API, the use of the Frictionless Data metadata standard, Python 3 support, JSON and XML backends, and autocompletion for the command line interface.

# References
