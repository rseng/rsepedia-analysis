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
======================
Spatial database setup
======================

Supporting installation of spatial data into `Postgres DBMS`.

1. **Install Postgres**

  For Mac the easiest way to get started with PostgreSQL is with `Postgres.app`_.

  For Debain and Ubuntu, install `PostgresSQL and PostGis` please ref to `Postgres installation`_.

  Otherwise you can try package installers for WINDOWS, MAC, Linux and MacOS from the main `PostgreSQL download`_ page

  For simplicity, use `.pgpass` file(`pgpass.conf` file for Microsoft Windows) to avoid supplying the password every time
  as decribed in `Passwordless configuration`_.

  After installation, Make sure you have the paths to these tools added to your system's `PATHS`.

  Note: Older version of this `raster2pgsql` was a python script that you had to download and manually include in Postgres's directory.
  Please consult an operating system expert for help on how to change or add the `PATH` variables.

  **For example, this could be a sample of paths exported on Mac:**

  #~/.bash_profile file, Postgres PATHS and tools .

.. code-block::

  export PATH="/Applications/Postgres.app/Contents/MacOS/bin:${PATH}"
  export PATH="$PATH:/Applications/Postgres.app/Contents/Versions/10/bin"


2. **Enable PostGIS extensions**

  If you have Postgres set up, enable `PostGIS extensions`.
  This is done by using either `Postgres CLI` or `GUI(PgAdmin)` and run the commands below.

**For psql CLI**

.. code-block::

  psql -d yourdatabase -c "CREATE EXTENSION postgis;"
  psql -d yourdatabase -c "CREATE EXTENSION postgis_topology;"

**For GUI(PgAdmin)**

.. code-block::

  `CREATE EXTENSION postgis;`
  `CREATE EXTENSION postgis_topology`

For more details refer to the `Postgis docs`_.

.. note::


    PostGIS excluded the raster types and functions from the main extension as of version 3.x;
    A separate CREATE EXTENSION postgis_raster; is then needed to get raster support.

    Versions 2.x have full raster support as part of the main extension environment,
    so CREATE EXTENSION postgis; is all that you need`


.. _PostgreSQL download: https://www.postgresql.org/download/
.. _Postgres.app: https://postgresapp.com/
.. _Postgres installation: https://trac.osgeo.org/postgis/wiki/UsersWikiPostGIS21UbuntuPGSQL93Apt
.. _Postgis docs: https://postgis.net/docs/postgis_installation.html#install_short_version
.. _Passwordless configuration: developer.html#passwordless-configuration
===========================
Using the retriever-recipes
===========================

retriever-recipes
-----------------

The `Data Retriever`_ earlier used a simple CLI for developing new dataset scripts. This allowed users with no programming experience to quickly add most standard datasets to the Retriever by specifying the names and locations of the tables along with additional information about the configuration of the data. The script is saved as a JSON file, that follows the DataPackage_ standards.

.. _Data Retriever: http://data-retriever.org
.. _DataPackage: https://specs.frictionlessdata.io/data-package/

This functionality has been moved to the ``retriever-recipes`` repository to separate the scripts from the core ``retriever`` functionalities to help with organization, maintenance, and testing. The `retriever recipes`_ repository thus holds all the scripts which were earlier shipped with ``retriever`` and also all the script adding/editing functionalities.

.. _retriever recipes: https://github.com/weecology/retriever-recipes

Installation
------------

The ``retriever-recipes`` project can be installed from Github using the following steps:

::

  git clone https://www.github.com/weecology/retriever-recipes.git
  cd retriever-recipes
  python setup.py install

Script Creation
---------------

To create a new script, there are 2 methods :- 

1. Use retriever autocreate to automatically create a script template. Specify the type of data using -dt, the default data type is tabular. Download the files to a folder. In case of tabular data, the files should be CSV files. Autocreate can create a script template for each file using -f or use -d to create a single script template for all files in the directory.

::


  usage: retriever autocreate [-h] [-dt [{raster,vector,tabular}]] [-f] [-d]
                              [-o [O]] [--skip-lines SKIP_LINES]
                              path

  positional arguments:
    path                  path to the data file(s)

  optional arguments:
    -h, --help            show this help message and exit
    -dt [{raster,vector,tabular}]
                          datatype for files
    -f                    turn files into scripts
    -d                    turn a directory and subdirectories into scripts
    -o [O]                write scripts out to a designated directory
    --skip-lines SKIP_LINES
                          skip a set number of lines before processing data



2. Manual script creation using ``retriever-recipes new_json``, which starts the CLI tool for new script creation.

``Required``

#. **name:** A one word name for the dataset

``Strongly recommended``

#. **title:** Give the name of the dataset
#. **description:** A brief description of the dataset of ~25 words.
#. **citation:** Give a citation if available
#. **homepage:** A reference to the data or the home page
#. **keywords:** Helps in classifying the type of data (i.e using Taxon, Data Type, Spatial Scale, etc.)


``Mandatory for any table added; Add Table? (y/N)``

#. **table-name:** Name of the table, URL to the table
#. **table-url:** Name of the table, URL to the table

Basic Scripts
-------------

The most basic scripts structure requires only some general metadata about the
dataset, i.e., the shortname of the database and table, and the location of the
table.

**Example of a basic script, example.script**

``Creating script from the CLI``

::

  name (a short unique identifier; only lowercase letters and - allowed): example-mammal
  title: Mammal Life History Database - Ernest, et al., 2003
  description:
  citation: S. K. Morgan Ernest. 2003. Life history characteristics of placental non-volant mammals. Ecology 84:3402.
  homepage (for the entire dataset):
  keywords (separated by ';'): mammals ; compilation

  Add Table? (y/N): y
  table-name: species
  table-url: http://esapubs.org/archive/ecol/E084/093/Mammal_lifehistories_v2.txt
  missing values (separated by ';'):
  replace_columns (separated by ';'):
  delimiter:
  do_not_bulk_insert (bool = True/False):
  contains_pk (bool = True/False):
  escape_single_quotes (bool = True/False):
  escape_double_quotes (bool = True/False):
  fixed_width (bool = True/False):
  header_rows (int):
  Enter columns [format = name, type, (optional) size]:


  Add crosstab columns? (y,N): n

  Add Table? (y/N): n

``Created script``

::

  {
      "citation": "S. K. Morgan Ernest. 2003. Life history characteristics of placental non-volant mammals. Ecology 84:3402.",
      "description": "",
      "homepage": "",
      "keywords": [
          "Mammals",
          "Compilation"
      ],
      "name": "example-mammal",
      "resources": [
          {
              "dialect": {},
              "name": "species",
              "schema": {
                  "fields": []
              },
              "url": "http://esapubs.org/archive/ecol/E084/093/Mammal_lifehistories_v2.txt"
          }
      ],
      "retriever": "True",
      "retriever_minimum_version": "2.0.dev",
      "title": "Mammal Life History Database - Ernest, et al., 2003"
      "version": "1.0.0"
  }

Explanation for the keys:

- ``citation``: Citation for the dataset
- ``description``: Description for the dataset
- ``homepage``: Homepage or website where the data is hosted
- ``keywords``: Keywords/tags for the dataset (for searching and classification)
- ``name``: Shortname for the dataset. Unique, URL-identifiable
- ``resources``: List of tables within the dataset

  - ``dialect``: Metadata for retriever to process the table

    - ``missingValues``: (Optional) List of strings which represents missing values in tables
    - ``delimiter``: (Optional) a character which represent boundary between two separate value(ex. ',' in csv files)
    - ``header_rows``: (Optional) number of header rows in table.
  - ``name``: Name of the table
  - ``schema``: List of the columns in the table

    - ``fields``: (Optional-Recommended) List of columns and their types and (optional) size values
    - ``ct_column``: (Optional) Cross-tab column with column names from dataset

  - ``url``: URL of the table

- ``retriever``: Auto generated tag for script identification
- ``retriever_minimum_version``: Minimum version that supports this script
- ``title``: Title/Name of the dataset
- ``urls``: dictionary of table names and the respective urls
- ``version``: "1.0.0"

Multiple Tables
---------------

A good example of data with multiple tables is Ecological Archives E091-124-D1, `McGlinn et al. 2010`_. ``plant-comp-ok`` Vascular plant composition data.
Since there are several csv files, we create a table for each of the files.

Assuming we want to call our dataset McGlinn2010, below is an example of the script that will handle this data

.. _`McGlinn et al. 2010`: http://esapubs.org/archive/ecol/E091/124/

::

  ...
    "name": "McGlinn2010",
    "resources": [
        {
            "dialect": {},
            "name": "pres",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E091/124/TGPP_pres.csv"
        },
        {
            "dialect": {},
            "name": "cover",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E091/124/TGPP_cover.csv"
        },
        {
            "dialect": {},
            "name": "richness",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E091/124/TGPP_rich.csv"
        },
        {
            "dialect": {},
            "name": "species",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E091/124/TGPP_specodes.csv"
        },
        {
            "dialect": {},
            "name": "environment",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E091/124/TGPP_env.csv"
        },
        {
            "dialect": {},
            "name": "climate",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E091/124/TGPP_clim.csv"
        }
    ],
    "retriever": "True",
    "retriever_minimum_version": "2.0.dev",
    "title": "Vascular plant composition - McGlinn, et al., 2010",
    ...

Null Values
-----------

The Retriever can replace non-standard null values by providing a semi-colon separated list of those null values
after the table in which the null values occur.

::

  ...
  Table name: species
  Table URL: http://esapubs.org/archive/ecol/E084/093/Mammal_lifehistories_v2.txt
  nulls (separated by ';'): -999 ; 'NA'
  ...

For example, the `Adler et al. 2010`_. ``mapped-plant-quads-ks`` script uses -9999 to indicate null values.

.. _`Adler et al. 2010`: http://esapubs.org/archive/ecol/E088/161/

::

  ...
        {
            "dialect": {},
            "name": "quadrat_info",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E088/161/quadrat_info.csv"
        },
        {
            "dialect": {
                "missingValues": [
                    "NA"
                ]
            },
  ...


Headers
-------

If the first row of a table is the headers then naming the columns will, be default, be handled automatically.
If you want to rename an existing header row for some reason, e.g.,
it includes reserved keywords for a database management system,
you can do so by adding a list of semi-colon separated column names,
with the new columns provided after a comma for each such column.

::

  ...
  Add Table? (y/N): y
  Table name: species
  Table URL: http://esapubs.org/archive/ecol/E091/124/TGPP_specodes.csv
  replace_columns (separated by ';', with comma-separated values): jan, january ; feb, february ; mar, march
  ...


The ``mapped-plant-quads-ks`` script for the `Adler et al. 2007`_. dataset from Ecological Archives
includes this functionality:


.. _`Adler et al. 2007`: http://esapubs.org/archive/ecol/E088/161/

::

  ...
   "name": "mapped-plant-quads-ks",
    "resources": [
        {
            "dialect": {},
            "name": "main",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E088/161/allrecords.csv"
        },
        {
            "dialect": {},
            "name": "quadrat_info",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E088/161/quadrat_info.csv"
        },
        {
            "dialect": {
                "missingValues": [
                    "NA"
                ]
            },
            "name": "quadrat_inventory",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E088/161/quadrat_inventory.csv"
        },
        {
            "dialect": {},
            "name": "species",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E088/161/species_list.csv"
        },
        {
            "dialect": {
                "missingValues": [
                    "NA"
                ],
                "replace_columns": [
                    [
                        "jan",
                        "january"
                    ],
                    [
                        "feb",
                        "february"
                    ],
                    [
                        "mar",
                        "march"
                    ],
                    [
                        "apr",
                        "april"
                    ],
                    [
                        "jun",
                        "june"
                    ],
                    [
                        "jul",
                        "july"
                    ],
                    [
                        "aug",
                        "august"
                    ],
                    [
                        "sep",
                        "september"
                    ],
                    [
                        "oct",
                        "october"
                    ],
                    [
                        "nov",
                        "november"
                    ],
                    [
                        "dec",
                        "december"
                    ]
                ]
            },
            "name": "monthly_temp",
            "schema": {},
            "url": "http://esapubs.org/archive/ecol/E088/161/monthly_temp.csv"
        },
    ...


Data Format
-----------

Data packages for different data formats can been added to Retriever now.
To add data format add keys in the script for Data sources except in the case of csv.

Data formats which can be added are :-

1. JSON Data :- For JSON raw data, add the key word ``json_data`` to the
resource. To add data formats for a given data package(`nuclear-power-plants`),
add keys to the resource part as described below.

::

  ...
   "name": "nuclear-power-plants",
    "resources": [
        {
            "dialect": {
                "delimiter": ","
            },
            "name": "nuclear_power_plants",
            "path": "nuclear_power_plants.csv",
            "json_data": "nuclear_power_plants.json",
            "schema": {
                "fields": [
                    {
                        "name": "id",
                        "type": "int"
                    },
                    {
                        "name": "name",
                        "size": "40",
                        "type": "char"
                    },
    ...

2. XML Data :- For XML raw data, add the key words ``xml_data`` and ``empty_rows``
to the resource. To add data formats for a given data package(`county-emergency-management-offices`), add keys to
the resource part as described below.

::

  ...
  "name": "county-emergency-management-offices",
    "resources": [
        {
            "dialect": {
                "delimiter": ","
            },
            "name": "county_emergency_management_offices",
            "path": "County_Emergency_Management_Offices.csv",
            "xml_data": "emergency_offices.xml",
            "empty_rows": 1,
            "schema": {
                "fields": [
                    {
                        "name": "county",
                        "size": "11",
                        "type": "char"
                    },
                    {
                        "name": "emergency_manager",
                        "size": "20",
                        "type": "char"
    ...

3. SQlite Data :- For SQlite raw data, add the key word ``sqlite_data`` to the
resource. To add data formats for a given data package(`portal-project-test`),
add keys to the resource part as described below.

::

  ...
   "name": "portal-project-test",
    "resources": [
        {
            "dialect": {
                "delimiter": ","
            },
            "name": "species",
            "path": "species.csv",
            "sqlite_data": ["species","portal_project.sqlite"],
            "schema": {
                "fields": [
                    {
                        "name": "species_id",
                        "size": "2",
                        "type": "char"
                    },
                    {
                        "name": "genus",
                        "size": "16",
                        "type": "char"
                    },
    ...


4. GeoJSON Data :- For GeoJSON raw data, add the key word ``geojson_data`` to the
resource. To add data formats for a given data package(`lake-county-illinois-cancer-rates`),
add keys to the resource part as described below.

::

  ...
   "name": "lake-county-illinois-cancer-rates",
    "resources": [
        {
            "dialect": {
                "delimiter": ","
            },
            "name": "lakecounty_health",
            "path": "LakeCounty_Health.csv",
            "format": "tabular",
            "geojson_data": "mytest.geojson",
            "schema": {
                "fields": [
                    {
                        "name": "fid",
                        "type": "int"
                    },
                    {
                        "name": "zip",
                        "type": "int"
                    },

    ...

5. HDF5 Data :- For HDF5 raw data, add the key word ``hdf5_data`` to the
resource. To add data formats for a given data package(`sample-hdf`),
add keys to the resource part as described below.

::

  ...
   "name": "sample-hdf",
  "title": "Test data raster bio1",
  "description": "Test data sampled from bioclim bio1",
  "citation": "N/A",
  "keywords": [
    "tests"
  ],
  "encoding": "latin-1",
  "url": "https://github.com/ashishpriyadarshiCIC/My_data/raw/main/Test_table_image_data.h5",
  "ref": "N/A",
  "version": "1.0.0",
  "retriever_minimum_version": "2.1.dev",
  "driver": "GTiff",
  "colums": 284,
  "rows": 249,
  "band_count": 1,
  "datum": "N/A - Coordinate Reference System",
  "projection": "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]]",
  "file_size": "N/A",
  "group_count": "N/A",
  "dataset_count": "N/A",
  "transform": {
    "xOrigin": -121.6250000000029,
    "pixelWidth": 0.041666666666667,
    "rotation_2": 0.0,
    "yOrigin": 42.79166666666632,
    "rotation_4": 0.0,
    "pixelHeight": -0.041666666666667
  },
  "resources": [
    {
      "dialect": {
        "delimiter": ","
      },
      "name": "table_data",
      "path": "table_data.csv",
      "hdf5_data": [
        "Test_table_image_data.h5",
        "csv",
        "G1/table_data"
      ],
      "schema": {
        "fields": [
          {
            "name": "region",
            "size": "33",
            "type": "char"
          },
          {
            "name": "country",
            "size": "32",
            "type": "char"
          },
          {
            "name": "item_type",
            "size": "15",
            "type": "char"
          },
          {
            "name": "sales_channel",
            "size": "7",
            "type": "char"
          },
          {
            "name": "order_id",
            "type": "int"
          },
          {
            "name": "total_profit",
            "type": "double"
          }
        ]
      },
      "url": "https://github.com/ashishpriyadarshiCIC/My_data/raw/main/Test_table_image_data.h5"
    },
    {
      "name": "test_image",
      "format": "raster",
      "hdf5_data": [
        "Test_table_image_data.h5",
        "image",
        "G1/G2/im"
      ],
      "path": "test_raster_bio1.tif",
      "extensions": [
        "tif"
      ],
      "no_data_value": -9999.0,
      "min": null,
      "max": null,
      "scale": 1.0,
      "color_table": null,
      "statistics": {
        "minimum": 0.0,
        "maximum": 0.0,
        "mean": 0.0,
        "stddev": -1.0
      },
    ...


Full control over column names and data types
---------------------------------------------

By default the Retriever automatically detects both column names and data types, but you can also exercise complete
control over the structure of the resulting database by adding column names and types.

It is recommended to describe the schema of the table while creating the JSON file. This enables processing of the data faster since column detection increases the processing time.

These values are stored in the ``fields`` array of the ``schema`` dict of the JSON script.

The ``fields`` value enables full control of the columns, which includes, renaming columns, skipping unwanted columns, mentioning primary key and combining columns.

The basic format for ``fields`` is as shown below:

::

  ...
  Enter columns [format = name, type, (optional) size]:

  count, int
  name, char, 40
  year, int
  ...


where ``name`` represents name of the column and ``type`` represents the type of data present in the column. The following can be used to describe the data type:

::

  pk-auto: Auto generated primary key starting from 1
  pk-[char,int,double]: primary key with data type
  char: strings
  int: integers
  double:floats/decimals
  ct-[int,double,char]:Cross tab data
  skip: used to skip the column in database


``pk-auto`` is used to create an additional column of type int which acts as a primary key with values starting from 1. While ``pk-[char,int,double]`` is used to make a primary key from existing columns of the table having data type of char/int/double.

The Smith et al. Masses of Mammals ``mammal-masses`` dataset script includes this type of functionality.

::

  ...
     "name": "mammal-masses",
    "resources": [
        {
            "dialect": {
                "missingValues": [
                    -999
                ],
                "header_rows": 0
            },
            "name": "MammalMasses",
            "schema": {
                "fields": [
                    {
                        "name": "record_id",
                        "type": "pk-auto"
                    },
                    {
                        "name": "continent",
                        "size": "20",
                        "type": "char"
                    },
                    {
                        "name": "status",
                        "size": "20",
                        "type": "char"
                    },
                    {
                        "name": "sporder",
                        "size": "20",
                        "type": "char"
                    },
                    {
                        "name": "family",
                        "size": "20",
                        "type": "char"
                    },
                    {
                        "name": "genus",
                        "size": "20",
                        "type": "char"
                    },
                    {
                        "name": "species",
                        "size": "20",
                        "type": "char"
                    },
                    {
                        "name": "log_mass_g",
                        "type": "double"
                    },
                    {
                        "name": "comb_mass_g",
                        "type": "double"
                    },
                    {
                        "name": "reference",
                        "type": "char"
                    }
                ]
            },
            "url": "http://www.esapubs.org/Archive/ecol/E084/094/MOMv3.3.txt"
        }
    ],
    "retriever": "True",
    "retriever_minimum_version": "2.0.dev",
    "title": "Masses of Mammals (Smith et al. 2003)",
  ...

Restructuring cross-tab data
----------------------------

It is common in ecology to see data where the rows indicate one level of grouping (e.g., by site),
the columns indicate another level of grouping (e.g., by species), and the values in each cell indicate
the value for the group indicated by the row and column (e.g., the abundance of species x at site y).
This is referred as cross-tab data and cannot be easily handled by database management systems,
which are based on a one record per line structure. The Retriever can restructure this type of
data into the appropriate form.
In scripts this involves telling the retriever the name of the column to store the data in
and the names of the columns to be restructured.

::

  ...
  Add crosstab columns? (y,N): y
  Crosstab column name: <name of column to store cross-tab data>
  Enter names of crosstab column values (Press return after each name):

  ct column 1
  ct column 2
  ct column 3
  ...

The `Moral et al 2010 script`_. ``mt-st-helens-veg`` takes advantage of this functionality.

.. _`Moral et al 2010 script`: http://esapubs.org/archive/ecol/E091/152/

::

  ...
  "name": "mt-st-helens-veg",
    "resources": [
        {
            "dialect": {
                "delimiter": ","
            },
            "name": "species_plot_year",
            "schema": {
                "ct_column": "species",
                "ct_names": [
                    "Abilas",
                    "Abipro",
                    "Achmil",
                    "Achocc",
                    "Agoaur",
                    "Agrexa",
                    "Agrpal",
                    "Agrsca",
                    "Alnvir",
                    "Anamar",
                    "Antmic",
                    "Antros",
                    "Aqifor",
                    "Arcnev",
                    "Arnlat",
                    "Astled",
                    "Athdis",
                    "Blespi",
                    "Brocar",
                    "Brosit",
                    "Carmer",
                    "Carmic",
                    "Carpac",
                    "Carpay",
                    "Carpha",
                    "Carros",
                    "Carspe",
                    "Casmin",
                    "Chaang",
                    "Cirarv",
                    "Cisumb",
                    "Crycas",
                    "Danint",
                    "Descae",
                    "Elyely",
                    "Epiana",
                    "Eriova",
                    "Eripyr",
                    "Fesocc",
                    "Fravir",
                    "Gencal",
                    "Hiealb",
                    "Hiegra",
                    "Hyprad",
                    "Junmer",
                    "Junpar",
                    "Juncom",
                    "Leppun",
                    "Lommar",
                    "Luepec",
                    "Luihyp",
                    "Luplat",
                    "Luplep",
                    "Luzpar",
                    "Maiste",
                    "Pencar",
                    "Pencon",
                    "Penser",
                    "Phahas",
                    "Phlalp",
                    "Phldif",
                    "Phyemp",
                    "Pincon",
                    "Poasec",
                    "Poldav",
                    "Polmin",
                    "Pollon",
                    "Poljun",
                    "Popbal",
                    "Potarg",
                    "Psemen",
                    "Raccan",
                    "Rumace",
                    "Salsit",
                    "Saxfer",
                    "Senspp",
                    "Sibpro",
                    "Sorsit",
                    "Spiden",
                    "Trispi",
                    "Tsumer",
                    "Vacmem",
                    "Vervir",
                    "Vioadu",
                    "Xerten"
                ],
                "fields": [
                    {
                        "name": "record_id",
                        "type": "pk-auto"
                    },
                    {
                        "name": "plot_id_year",
                        "size": "20",
                        "type": "char"
                    },
                    {
                        "name": "plot_name",
                        "size": "4",
                        "type": "char"
                    },
                    {
                        "name": "plot_number",
                        "type": "int"
                    },
                    {
                        "name": "year",
                        "type": "int"
                    },
                    {
                        "name": "count",
                        "type": "ct-double"
                    }
                ]
            },
            "url": "http://esapubs.org/archive/ecol/E091/152/MSH_SPECIES_PLOT_YEAR.csv"
        },
  ...

Script Editing
--------------
**Note:** Any time a script gets updated, the minor version number must be incremented from within the script.

The JSON scripts created using the retriever-recipes CLI can also be edited using the CLI.

To edit a script, use the ``retriever-recipes edit_json`` command, followed by the script's shortname;

For example, editing the ``mammal-life-hist`` (Mammal Life History Database - Ernest, et al., 2003)
dataset, the editing tool will ask a series a questions for each of the keys and values of the script,
and act according to the input.


The tool describes the values you want to edit.
In the script below the first keyword is citation, ``citation ( <class 'str'> )``
and it is of class string or expects a string.

::

  dev@retriever:~$ retriever-recipes edit_json mammal-life-hist

    ->citation ( <class 'str'> ) :

    S. K. Morgan Ernest. 2003. Life history characteristics of placental non-volant mammals. Ecology 84:3402

    Select one of the following for the key 'citation'

    1. Modify value
    2. Remove from script
    3. Continue (no changes)


    Your choice: 3

      ->homepage ( <class 'str'> ) :

      http://esapubs.org/archive/ecol/E084/093/


    Select one of the following for the key 'homepage':

    1. Modify value
    2. Remove from script
    3. Continue (no changes)


    Your choice: 3

      ->description ( <class 'str'> ) :

      The purpose of this data set was to compile general life history characteristics for a variety of mammalian
      species to perform comparative life history analyses among different taxa and different body size groups.


    Select one of the following for the key 'description':

    1. Modify value
    2. Remove from script
    3. Continue (no changes)


    Your choice: 3

      ->retriever_minimum_version ( <class 'str'> ) :

      2.0.dev


    Select one of the following for the key 'retriever_minimum_version':

    1. Modify value
    2. Remove from script
    3. Continue (no changes)


    Your choice: 3

      ->version ( <class 'str'> ) :

      1.1.0


    Select one of the following for the key 'version':

    1. Modify value
    2. Remove from script
    3. Continue (no changes)


    Your choice: 3

      ->resources ( <class 'list'> ) :

      {'dialect': {}, 'schema': {}, 'name': 'species', 'url': 'http://esapubs.org/archive/ecol/E084/093/Mammal_lifehistories_v2.txt'}


    1 .  {'dialect': {}, 'schema': {}, 'name': 'species', 'url': 'http://esapubs.org/archive/ecol/E084/093/Mammal_lifehistories_v2.txt'}

    Edit this dict in 'resources'? (y/N): n
    Select one of the following for the key 'resources':

    1. Add an item
    2. Delete an item
    3. Remove from script
    4. Continue (no changes)
    ...
===========================
Data Retriever using Python
===========================


`Data Retriever <http://data-retriever.org>`_ is written purely in `python <http://www.python.org/>`_.
The Python interface provides the core functionality supported by the CLI (Command Line Interface).



Installation
============

The installation instructions for the CLI and module are the same. Links have been provided below for convenience.

- Instructions for installing from binaries `project website <http://data-retriever.org>`_.
- Instructions for installing from source  `install from Source <https://github.com/weecology/retriever>`_.

Note: The python interface requires version 2.1 and above.

Tutorial
========

Importing retriever

.. code-block:: python

  >>> import retriever

In this tutorial, the module will be referred to as ``rt``.

.. code-block:: python

  >>> import retriever as rt

List Datasets
=============

Listing available datasets using ``dataset_names`` function.
The function returns a list of all the currently available scripts.

.. code-block:: python

  >>> rt.dataset_names()

  ['abalone-age',
   'antarctic-breed-bird',
   .
   .
   'wine-composition',
   'wine-quality']


For a more detailed description of the scripts installed in retriever, the ``datasets`` function can be used. This function returns a list of ``Scripts`` objects.
From these objects, we can access the available Script's attributes as follows.

.. code-block:: python

  >>> for dataset in rt.datasets():
        print(dataset.name)

  abalone-age
  airports
  amniote-life-hist
  antarctic-breed-bird
  aquatic-animal-excretion
  .
  .

There are a lot of different attributes provided in the Scripts class. Some notably useful ones are:

.. code-block:: python

  - name
  - citation
  - description
  - keywords
  - title
  - urls
  - version

You can add more datasets locally by yourself.
`Adding dataset <http://retriever.readthedocs.io/en/latest/scripts.html>`_ documentation.

Update Datasets
===============

If there are no scripts available, or you want to update scripts to the latest version,
``check_for_updates`` will download the most recent version of all scripts.


.. code-block:: python

  >>> rt.check_for_updates()

  Downloading scripts...
  Download Progress: [####################] 100.00%
  The retriever is up-to-date


Downloading recipes for all datasets can take a while depending on the internet connection.

Download Datasets
=================

To directly download datasets without cleaning them use the ``download`` function

.. code-block:: python

  def download(dataset, path='./', quiet=False, subdir=False, debug=False):

A simple download for the ``iris`` dataset can be done using the following.

.. code-block:: python

  >>> rt.download("iris")

Output:

.. code-block:: python

  => Downloading iris

  Downloading bezdekIris.data...
  100%  0 seconds Copying bezdekIris.data

The files will be downloaded into your current working directory by default.
You can change the default download location by using the ``path`` parameter.
Here, we are downloading the ``NPN`` dataset to our ``Desktop`` directory

.. code-block:: python

  >>> rt.download("NPN","/Users/username/Desktop")

Output:

.. code-block:: python

  => Downloading NPN

  Downloading 2009-01-01.xml...
  11  MBB
  Downloading 2009-04-02.xml...
  42  MBB
  .
  .


.. code-block:: python

  path (String): Specify dataset download path.

  quiet  (Bool): Setting True minimizes the console output.

  subdir (Bool): Setting True keeps the subdirectories for archived files.

  debug  (Bool): Setting True helps in debugging in case of errors.

Install Datasets
================

Retriever supports installation of datasets into 7 major databases and file formats.

.. code-block:: python

  - csv
  - json
  - msaccess
  - mysql
  - postgres
  - sqlite
  - xml


There are separate functions for installing into each of the 7 backends:

.. code-block:: python

    def install_csv(dataset, table_name=None, compile=False, debug=False,
                quiet=False, use_cache=True):

    def install_json(dataset, table_name=None, compile=False,
                 debug=False, quiet=False, use_cache=True, pretty=False):

    def install_msaccess(dataset, file=None, table_name=None,
                     compile=False, debug=False, quiet=False, use_cache=True):

    def install_mysql(dataset, user='root', password='', host='localhost',
                  port=3306, database_name=None, table_name=None,
                  compile=False, debug=False, quiet=False, use_cache=True):

    def install_postgres(dataset, user='postgres', password='',
                     host='localhost', port=5432, database='postgres',
                     database_name=None, table_name=None,
                     compile=False, debug=False, quiet=False, use_cache=True):

    def install_sqlite(dataset, file=None, table_name=None,
                   compile=False, debug=False, quiet=False, use_cache=True):

    def install_xml(dataset, table_name=None, compile=False, debug=False,
                quiet=False, use_cache=True):

A description of default parameters mentioned above:

.. code-block:: python

  compile         (Bool): Setting True recompiles scripts upon installation.

  database_name (String): Specify database name. For postgres, mysql users.

  debug           (Bool): Setting True helps in debugging in case of errors.

  file          (String): Enter file_name for database. For msaccess, sqlite users.

  host          (String): Specify host name for database. For postgres, mysql users.

  password      (String): Specify password for database. For postgres, mysql users.

  port             (Int): Specify the port number for installation. For postgres, mysql users.

  pretty          (Bool): Setting True adds indentation in JSON files.

  quiet           (Bool): Setting True minimizes the console output.

  table_name    (String): Specify the table name to install.

  use_cache       (Bool): Setting False reinstalls scripts even if they are already installed.

  user          (String): Specify the username. For postgres, mysql users.

Examples to Installing Datasets:

Here, we are installing the dataset wine-composition as a CSV file in our current working directory.

.. code-block:: python

  rt.install_csv("wine-composition")

  => Installing wine-composition

  Downloading wine.data...
  100%  0 seconds Progress: 178/178 rows inserted into ./wine_composition_WineComposition.csv totaling 178

The installed file is called ``wine_composition_WineComposition.csv``

Similarly, we can download any available dataset as a JSON file:

.. code-block:: python

  rt.install_json("wine-composition")

  => Installing wine-composition

  Progress: 178/178 rows inserted into ./wine_composition_WineComposition.json totaling 17

The wine-composition dataset is now installed as a JSON file called wine_composition_WineComposition.json in our current working directory.
=================
Developer's guide
=================

1. Quickstart by forking the main repository https://github.com/weecology/retriever
2. Clone your copy of the repository

    - Using https ``git clone https://github.com/henrykironde/retriever.git``
    - Using ssh ``git clone git@github.com:henrykironde/retriever.git``

3. Link or point your cloned copy to the main repository. (I always name it upstream)

    - ``git remote add upstream https://github.com/weecology/retriever.git``

5. Check/confirm your settings using ``git remote -v``

::

    origin	git@github.com:henrykironde/retriever.git (fetch)
    origin	git@github.com:henrykironde/retriever.git (push)
    upstream	https://github.com/weecology/retriever.git (fetch)
    upstream	https://github.com/weecology/retriever.git (push)

6. Install the package from the main directory.
use `-U or --upgrade` to upgrade or overwrite any previously installed versions.

::

    pip install . -U

7. Check if the package was installed

::

    retriever ls
    retriever -v

8. Run sample test on  CSV engine only, with the option `-k`

::

   pip install pytest
   pytest -k "CSV" -v


Required Modules
================

You will need Python 3.6.8+
Make sure the required modules are installed: ``Pip install -r requirements.txt``

Developers need to install these extra packages.

::

   pip install codecov
   pip install pytest-cov
   pip install pytest-xdist
   pip install pytest
   pip install yapf
   pip install pylint
   pip install flake8
   Pip install pypyodbc # For only Windows(MS Access)

Setting up servers
==================

You need to install all the database infrastructures to enable local testing.


`PostgresSQL`_
`MySQL`_
`SQLite`_
`MSAccess`_ (For only Windows, MS Access)

After installation, configure passwordless access to MySQL and PostgresSQL Servers

Passwordless configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^

To avoid supplying the passwords when using the tool, use the config files
`.pgpass`(`pgpass.conf` for Microsoft Windows) for Postgres and `.my.cnf`
for MySQL.

Create if not exists, and add/append the configuration details as below.
PostgresSQL conf file `~/.pgpass` file.

For more information regarding Passwordless configuration you can visit `PostgreSQL Password File`_ 
and `MySQL Password File`_ 

::

  localhost:*:*:postgres:Password12!

**Postgress:**

(Linux / Macos):- A `.pgpass` file in your HOME directory(~)

(WINDOWS 10-) - A `pgpass.conf` in your HOME directory(~)

(WINDOWS 10+):- Entering `%APPDATA%` will take you to `C:/\Users/\username/\AppData/\Roaming`.

In this directory create a new subdirectory named `postgresql`. Then create the `pgpass.conf` file inside it. On Microsoft Windows, it is assumed that the file is stored in a secure directory, hence no special permissions setting is needed.

Make sure you set the file permissions to 600

::

  # Linux / Macos
  chmod 600 ~/.pgpass
  chmod 600 ~/.my.cnf

For most of the recent versions of **Postgress server 10+**, you need to find `pg_hba.conf`. This file is located in the installed Postgres directory.
One way to find the location of the file `pg_hba.conf` is using ``psql -t -P format=unaligned -c 'show hba_file';``
To allow passwordless login to Postgres, change peer to `trust` in `pg_hba.conf` file.

::

  # Database administrative login by Unix domain socket
  local   all             postgres                                trust

Run commands in terminal to create user
::

  PostgreSQL
  ----------
  psql -c "CREATE USER postgres WITH PASSWORD 'Password12!'"
  psql -c 'CREATE DATABASE testdb_retriever'
  psql -c 'GRANT ALL PRIVILEGES ON DATABASE testdb_retriever to postgres'

Restart the server and test Postgress passwordless setup using retriever without providing the password

``retriever install postgres iris``

**MySQL:** Create if not exists `.my.cnf` in your HOME directory(~).
Add the configuration info to the MySQL conf file `~.my.cnf` file.

::

  [client]
  user="travis"
  password="Password12!"
  host="mysqldb"
  port="3306"

Run commands in terminal to create user
::

  MySQL
  -----
  mysql -e "CREATE USER 'travis'@'localhost';" -uroot
  mysql -e "GRANT ALL PRIVILEGES ON *.* TO 'travis'@'localhost';" -uroot
  mysql -e "GRANT FILE ON *.* TO 'travis'@'localhost';" -uroot

 Restart the server and test MySQL passwordless setup using retriever without providing the password

``retriever install mysql iris``

Testing
=======

Before running the tests make sure Postgis is set up `Spatial database setup`_.

Follow these instructions to run a complete set of tests for any branch
Clone the branch you want to test.

Two ways of installing the program using the `setup tools`_.

we can either install from source as

.. code-block:: bash

  pip install . --upgrade or python setup.py install

or install in development mode.

.. code-block:: bash

  python setup.py develop

For more about `installing`_ refer to the python setuptools `documentation`_.

you can also install from Git.

.. code-block:: bash

  # Local repository
  pip install git+file:///path/to/your/git/repo #  test a PIP package located in a local git repository
  pip install git+file:///path/to/your/git/repo@branch  # checkout a specific branch by adding @branch_name at the end

  # Remote GitHub repository
  pip install git+git://github.com/myuser/myproject  #  package from a GitHub repository
  pip install git+git://github.com/myuser/myproject@my_branch # github repository Specific branch


Running tests locally
^^^^^^^^^^^^^^^^^^^^^

Services Used
-------------

`Read The Docs`_,
`codecov`_,
`AppVeyor`_

From the source top-level directory, Use Pytest as examples below

.. code-block:: sh

  $   py.test -v # All tests
  $   py.test -v -k"csv" # Specific test with expression csv
  $   py.test ./test/test_retriever.py # Specific file

In case ``py.test`` requests for Password (even after Passwordless configuration), change the owner and group
permissions for the config files ``~/.pgpass, ~/.my.cnf``

Style Guide for Python Code
---------------------------

Use ``yapf -d --recursive retriever/ --style=.style.yapf`` to check style.

Use ``yapf -i --recursive retriever/ --style=.style.yapf`` refactor style

Continuous Integration
^^^^^^^^^^^^^^^^^^^^^^

The main GitHub repository runs the test on both the GitHub Actions (Linux) and AppVeyor
(Windows) continuous-integration platforms.

Pull requests submitted to the repository will automatically be tested using
these systems and results reported in the ``checks`` section of the pull request
page.


Create Release
==============

Start
^^^^^

1. **Run the tests**. Seriously, do it now.
2. Update ``CHANGES.md`` with major updates since the last release
3. Run ``python version.py`` (this will update ``version.txt``)
4. In the `main` branch update the version number and create a tag, run `bumpversion release`
5. Push the release commit and the tag
6. After the release, update the version to dev, run `bumpversion patch`

   ::

       git push upstream main
       git push upstream --tags

Pypi
^^^^

You will need to create an API key on PyPI and store it in ~/.pypirc to upload to PyPI.

1. `sudo python setup.py sdist bdist_wheel`
2. `sudo python -m twine upload -r pypi dist/*`

Cleanup
^^^^^^^

1. Bump the version numbers as needed. The version number is located in the ``setup.py``,
   ``retriever_installer.iss``, ``version.txt`` and ``retriever/_version.py``

Mac OSX Build
=============

Building the Retriever on OSX.

Python binaries
^^^^^^^^^^^^^^^

This build will allow you to successfully build the Mac App for
distribution to other systems.

1. Install the Python 3 Installer (or Python 2 if you have a specific reason for doing so)
   from the `Python download site`_.
2. Use pip to install any desired optional dependencies ``pip install pymysql psycopg2-binary pyinstaller pytest``
   You will need all of these dependencies, for example pyinstaller, if you want to build the Mac App for distribution

Homebrew
^^^^^^^^

Homebrew works great if you just want to install the Retriever from
source on your machine, but at least based on this recipe it does
not support the distribution of the Mac App to other versions of OS X (i.e.,
if you build the App on OS X 10.9 it will only run on 10.9)

1.  Install Homebrew
    ``ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go)"``
2.  Install Xcode
3.  Install Python ``brew install python``
4.  Install the Xcode command-line tools ``xcode-select --install``
5.  Make brew’s Python the default
    ``echo export PATH='usr/local/bin:$PATH' >> ~/.bash_profile``
6.  Install xlrd via pip ``pip install xlrd``. No ``sudo`` is necessary since we’re using brew.
7.  Clone the Retriever
    ``git clone git@github.com:weecology/retriever.git``
8. Switch directories ``cd retriever``
9. Standard install ``pip install . --upgrade``

If you also want to install the dependencies for MySQL and PostgreSQL
this can be done using a combination of homebrew and pip.

1. ``brew install mysql``
2. Follow the instructions from ``brew`` for starting MySQL
3. ``brew install postgresql``
4. Follow the instructions from ``brew`` for starting Postgres
5. ``sudo pip install pymysql MySQL-python psycopg2-binary``

``MySQL-python`` should be installed in addition to ``pymysql`` for
building the ``.app`` file since pymysql is not currently working
properly in the ``.app``.

Conda
^^^^^

-  This hasn’t been tested yet

.. _Python download site: http://www.python.org/download/



Creating or Updating a Conda Release
====================================

To create or update a Conda Release, first fork the conda-forge `retriever-feedstock repository <https://github.com/conda-forge/retriever-feedstock>`_.

Once forked, open a pull request to the retriever-feedstock repository. Your package will be tested on Windows, Mac, and Linux.

When your pull request is merged, the package will be rebuilt and become automatically available on conda-forge.

All branches in the conda-forge/retriever-feedstock are created and uploaded immediately, so PRs should be based on branches in forks. Branches in the main repository shall be used to build distinct package versions only.

For producing a uniquely identifiable distribution:

 - If the version of a package is not being incremented, then the build/number can be added or increased.
 - If the version of a package is being incremented, then remember to change the build/number back to 0.

Documentation
=============

We are using `Sphinx`_ and `Read the Docs`_. for the documentation.
Sphinx uses reStructuredText as its markup language.
Source Code documentation is automatically included after committing to the main.
Other documentation (not source code) files are added as new reStructuredText in the docs folder

In case you want to change the organization of the Documentation, please refer to `Sphinx`_

**Update Documentation**

The documentation is automatically updated for changes within modules.
However, the documentation should be updated after the addition of new modules in the engines or lib directory.
Change to the docs directory and create a temporary directory, i.e. ``source``.
Run

.. code-block:: bash

  cd  docs
  mkdir source
  sphinx-apidoc -f  -o ./source /Users/../retriever/

The ``source`` is the destination folder for the source rst files. ``/Users/../retriever/`` is the path to where
the retriever source code is located.
Copy the ``.rst`` files that you want to update to the docs directory, overwriting the old files.
Make sure you check the changes and edit if necessary to ensure that only what is required is updated.
Commit and push the new changes.
Do not commit the temporary source directory.

**Test Documentation locally**

.. code-block:: bash

  cd  docs  # go the docs directory
  make html && python3 -m http.server --directory _build/html
  # Makes the html files and hosts a HTTP server on localhost:8000 to view the documentation pages locally

.. note::
  Do not commit the _build directory after making HTML.

**Read The Docs configuration**

Configure read the docs (advanced settings) so that the source is first installed then docs are built.
This is already set up but could be changed if need be.

Collaborative Workflows with GitHub
===================================

First fork the `Data Retriever repository`_.
Then Clone your forked version with either HTTPS or SSH

   ::

      # Clone with HTTPS
      git clone https://github.com/[myusername]/retriever.git
      # Clone with SSH
      git clone git@github.com:[myusername]/retriever.git

This will update your `.git/config` to point to your repository copy of the Data Retriever as `remote "origin"`

   ::

       [remote "origin"]
       url = git@github.com:[myusername]/retriever.git
       fetch = +refs/heads/*:refs/remotes/origin/*

Point to Weecology `Data Retriever repository`_ repo.
This will enable you to update your main(origin) and you can then push to your origin main.
In our case, we can call this upstream().

   ::

      git remote add upstream https://github.com/weecology/retriever.git

This will update your `.git/config` to point to the Weecology `Data Retriever repository`_.

.. code-block:: bash

  [remote "upstream"]
  url = https://github.com/weecology/retriever.git
  fetch = +refs/heads/*:refs/remotes/upstream/*
  # To fetch pull requests add
  fetch = +refs/pull/*/head:refs/remotes/origin/pr/*

Fetch upstream main and create a branch to add the contributions to.

.. code-block:: bash

  git fetch upstream
  git checkout main
  git reset --hard upstream main
  git checkout -b [new-branch-to-fix-issue]

**Submitting issues**

Categorize the issues based on labels. For example (Bug, Dataset Bug, Important, Feature Request, etc..)
Explain the issue explicitly with all details, giving examples and logs where applicable.

**Commits**

From your local branch of retriever, commit to your origin.
Once tests have passed you can then make a pull request to the retriever main (upstream)
For each commit, add the issue number at the end of the description with the tag ``fixes #[issue_number]``.

Example
::

  Add version number to postgres.py to enable tracking

  Skip a line and add more explanation if needed
  fixes #3

**Clean history**

Make one commit for each issue.
As you work on a particular issue, try adding all the commits into one general commit rather than several commits.

Use ``git commit --amend`` to add new changes to a branch.

Use ``-f`` flag to force pushing changes to the branch. ``git push -f origin [branch_name]``


.. _codecov: https://codecov.io/
.. _project website: http://data-retriever.org
.. _Sphinx: http://www.sphinx-doc.org/en/stable/
.. _Read The Docs: https://readthedocs.org/
.. _AppVeyor: https://www.appveyor.com/
.. _documentation: https://pythonhosted.org/an_example_pypi_project/setuptools.html
.. _installing: https://docs.python.org/3.6/install/
.. _installing the wheel: http://www.lfd.uci.edu/~gohlke/pythonlibs/
.. _setup tools: https://pythonhosted.org/an_example_pypi_project/setuptools.html
.. _Data Retriever repository: https://github.com/weecology/retriever
.. _Spatial database setup: developer.html#Spatial-database-setup
.. _PostgresSQL: https://www.postgresql.org/download/
.. _SQlite: https://sqlite.org/download.html
.. _MySQL: https://www.mysql.com/downloads/
.. _MSAccess: https://www.microsoft.com/en-ww/microsoft-365/access
.. _PostgreSQL Password File : https://www.postgresql.org/docs/current/libpq-pgpass.html
.. _MySQL Password File : https://dev.mysql.com/doc/refman/8.0/en/option-files.html=======================
Using the Rdatasets API
=======================

This tutorial explains the usage of the Rdatasets API in Data Retriever. It includes both the
CLI (Command Line Interface) commands as well as the Python interface for the same.


Command Line Interface
======================

Listing the Rdatasets
---------------------

The ``retriever ls rdataset`` command displays the Rdatasets. 

$ ``retriever ls rdataset -h`` (gives listing options)

::

    usage: retriever ls rdataset [-h] [-p P [P ...]] all

    positional arguments:
        all           display all the packages present in rdatasets

    optional arguments:
        -h, --help    show this help message and exit
        -p P [P ...]  display a list of all rdatasets present in the package(s)

**Examples**

This example will display all the Rdatasets present with their package name, dataset name and script name

$ ``retriever ls rdataset``

::

    List of all available Rdatasets

    Package: aer              Dataset: affairs                   Script Name: rdataset-aer-affairs
    Package: aer              Dataset: argentinacpi              Script Name: rdataset-aer-argentinacpi
    Package: aer              Dataset: bankwages                 Script Name: rdataset-aer-bankwages
    ...
    Package: vcd              Dataset: vonbort                   Script Name: rdataset-vcd-vonbort
    Package: vcd              Dataset: weldondice                Script Name: rdataset-vcd-weldondice
    Package: vcd              Dataset: womenqueue                Script Name: rdataset-vcd-womenqueue


This example will display all the Rdatasets present in the packages ``vcd`` and ``aer``

$ ``retriever ls rdataset -p vcd aer``

::

    List of all available Rdatasets in packages: ['vcd', 'aer']
    Package: vcd              Dataset: arthritis                 Script Name: rdataset-vcd-arthritis
    Package: vcd              Dataset: baseball                  Script Name: rdataset-vcd-baseball
    Package: vcd              Dataset: brokenmarriage            Script Name: rdataset-vcd-brokenmarriage
    ...
    Package: aer              Dataset: affairs                   Script Name: rdataset-aer-affairs
    Package: aer              Dataset: argentinacpi              Script Name: rdataset-aer-argentinacpi
    Package: aer              Dataset: bankwages                 Script Name: rdataset-aer-bankwages
    ...


This example will display all the Rdatasets present in the package ``vcd``

$ ``retriever ls rdataset -p vcd``

::

    List of all available Rdatasets in packages: ['vcd', 'aer']
    Package: vcd              Dataset: arthritis                 Script Name: rdataset-vcd-arthritis
    Package: vcd              Dataset: baseball                  Script Name: rdataset-vcd-baseball
    Package: vcd              Dataset: brokenmarriage            Script Name: rdataset-vcd-brokenmarriage
    ...


This example will display all the packages present in rdatasets

$ ``retriever ls rdataset all``

::

    List of all the packages present in Rdatasets

    aer         cluster   dragracer  fpp2           gt        islr     mass        multgee         plyr      robustbase  stevedata  
    asaur       count     drc        gap            histdata  kmsurv   mediation   nycflights13    pscl      rpart       survival   
    boot        daag      ecdat      geepack        hlmdiag   lattice  mi          openintro       psych     sandwich    texmex     
    cardata     datasets  evir       ggplot2        hsaur     lme4     mosaicdata  palmerpenguins  quantreg  sem         tidyr      
    causaldata  dplyr     forecast   ggplot2movies  hwde      lmec     mstate      plm             reshape2  stat2data   vcd 


Downloading the Rdatasets
-------------------------

The ``retriever download rdataset-<package>-<dataset>`` command downloads the Rdataset ``dataset`` which exists in the package ``package``.
You can also copy the script name from the output of ``retriever ls rdataset``.

**Example**

This example downloads the ``rdataset-vcd-bundesliga`` dataset.

$ ``retriever download rdataset-vcd-bundesliga``

::

    => Installing rdataset-vcd-bundesliga
    Downloading Bundesliga.csv: 60.0B [00:00, 117B/s]                                                                                                 
    Done!

The downloaded raw data files are stored in the ``raw_data`` directory in the ``~/.retriever`` directory.


Installing the Rdatasets
------------------------

The ``retriever install <engine> rdataset-<package>-<dataset>`` command downloads the raw data, creates the script for it and then installs
the Rdataset ``dataset`` present in the package ``package`` into the provided ``engine``.

**Example**

This example install the ``rdataset-aer-usmoney`` dataset into the ``postgres`` engine.

$ ``retriever install postgres rdataset-aer-usmoney``

::

    => Installing rdataset-aer-usmoney
    Downloading USMoney.csv: 1.00B [00:00, 2.52B/s]
    Processing... USMoney.csv
    Successfully wrote scripts to /home/user/.retriever/rdataset-scripts/usmoney.csv.json
    Updating script name to rdataset-aer-usmoney.json
    Updating the contents of script rdataset-aer-usmoney
    Successfully updated rdataset_aer_usmoney.json
    Updated the script rdataset-aer-usmoney
    Creating database rdataset_aer_usmoney...

    Installing rdataset_aer_usmoney.usmoney
    Progress: 100%|█████████████████████████████████████████████████████████████████████████████████████████████| 136/136 [00:00<00:00, 2225.09rows/s]
    Done!

The script created for the Rdataset is stored in the ``rdataset-scripts`` directory in the ``~/.retriever`` directory.


Python Interface in Data Retriever
==================================

Updating Rdatasets Catalog
--------------------------

The function ``update_rdataset_catalog`` creates/updates the ``datasets_url.json`` in the ``~/.retriever/rdataset-scripts`` directory,
which contains the information about all the Rdatasets.

.. code-block:: python

  >>> import retriever as rt
  >>> rt.update_rdataset_catalog()

.. note::

  The ``update_rdataset_catalog`` function has a default argument ``test`` which is set to ``False``.
  If ``test`` is set to ``True``, then the contents of the ``datasets_url.json`` file would be returned as
  a dict. 

Listing Rdatasets
-----------------

The function ``display_all_rdataset_names`` prints the package, dataset name and the script name for the Rdatasets present in the package(s) requested.
If no package is specified, it prints all the rdatasets, and if ``all`` is passed as the function argument then all the package names are displayed.

.. note::

  The function argument ``package_name`` takes a list as an input when you want to display rdatasets based on the packages.
  If you want to display all packages names, set ``package_name`` argument to ``all`` (refer to the example below).


.. code-block:: python

  >>> import retriever as rt
  >>>
  >>> # Display all Rdatasets
  >>> rt.display_all_rdataset_names()
  List of all available Rdatasets

  Package: aer              Dataset: affairs                   Script Name: rdataset-aer-affairs
  Package: aer              Dataset: argentinacpi              Script Name: rdataset-aer-argentinacpi
  Package: aer              Dataset: bankwages                 Script Name: rdataset-aer-bankwages
  ...
  Package: vcd              Dataset: vonbort                   Script Name: rdataset-vcd-vonbort
  Package: vcd              Dataset: weldondice                Script Name: rdataset-vcd-weldondice
  Package: vcd              Dataset: womenqueue                Script Name: rdataset-vcd-womenqueue
  >>>
  >>> # Display all the Rdatasets present in packages 'aer' and 'drc'
  >>> rt.display_all_rdataset_names(['aer', 'drc'])
  List of all available Rdatasets in packages: ['aer', 'drc']
  Package: aer              Dataset: affairs                   Script Name: rdataset-aer-affairs
  Package: aer              Dataset: argentinacpi              Script Name: rdataset-aer-argentinacpi
  Package: aer              Dataset: bankwages                 Script Name: rdataset-aer-bankwages
  ...
  Package: drc              Dataset: spinach                   Script Name: rdataset-drc-spinach
  Package: drc              Dataset: terbuthylazin             Script Name: rdataset-drc-terbuthylazin
  Package: drc              Dataset: vinclozolin               Script Name: rdataset-drc-vinclozolin
  >>>
  >>> # Display all the packages in Rdatasets
  >>> rt.display_all_rdataset_names('all')
  List of all the packages present in Rdatasets

  aer         cluster   dragracer  fpp2           gt        islr     mass        multgee         plyr      robustbase  stevedata  
  asaur       count     drc        gap            histdata  kmsurv   mediation   nycflights13    pscl      rpart       survival   
  boot        daag      ecdat      geepack        hlmdiag   lattice  mi          openintro       psych     sandwich    texmex     
  cardata     datasets  evir       ggplot2        hsaur     lme4     mosaicdata  palmerpenguins  quantreg  sem         tidyr      
  causaldata  dplyr     forecast   ggplot2movies  hwde      lmec     mstate      plm             reshape2  stat2data   vcd 


Downloading a Rdataset
----------------------

.. code-block:: python

  >>> import retriever as rt
  >>> rt.download('rdataset-drc-earthworms')


Installing a Rdataset
---------------------

.. code-block:: python

  >>> import retriever as rt
  >>> rt.install_postgres('rdataset-mass-galaxies')

.. note::

  For downloading or installing the Rdatasets, the script name should follow the syntax given below.
  The script name should be ``rdataset-<package name>-<dataset name>``. The ``package name`` and ``dataset name``
  should be valid.

  Example:
    - Correct: ``rdataset-drc-earthworms``

    - Incorrect:  ``rdataset-drcearthworms``, ``rdatasetdrcearthworms``===============================
Using the Data Retriever from R
===============================

rdataretriever_
===============

The rdataretriever provides an R interface to the Data Retriever so
that the ``retriever``'s data handling can easily be integrated into R workflows.

Installation
============

To use the R package rdataretriever_, you first need to install the `Data Retriever`_.

The rdataretriever_ can then be installed using
``install.packages("rdataretriever")``

To install the development version, use ``devtools``

::

  # install.packages("devtools")
  library(devtools)
  install_github("ropensci/rdataretriever")

Note: The R package takes advantage of the Data Retriever's command line
interface, which must be available in the path. This path is given to the
rdataretriever_ using the function ``use_RetrieverPath()``. The location of
``retriever`` is dependent on the Python installation (Python.exe, Anaconda, Miniconda),
the operating system and the presence of virtual environments in the system. The following instances
exemplify this reliance and how to find retriever's path.

Ubuntu OS with default Python:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If ``retriever`` is installed in default Python, it can be found out in the system with the help
of ``which`` command in the terminal. For example:

::

  $ which retriever
  /home/<system_name>/.local/bin/retriever

The path to be given as input to ``use_RetrieverPath()`` function is */home/<system_name>/.local/bin/*
as shown below:

::

  library(rdataretriever)
  use_RetrieverPath("/home/<system_name>/.local/bin/")

The ``which`` command in the terminal finds the location of ``retriever`` including the name
of the program, but the path required by the function is the directory that contains ``retriever``.
Therefore, the `retriever` needs to be removed from the path before using it.

Ubuntu OS with Anaconda environment:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``retriever`` is installed in an virtual environment, the user can track its location only
when that particular environment is activated. To illustrate, assume the virtual environment is *py27*:

::

  $ conda activate py27
  (py27) $ which retriever
  /home/<system_name>/anaconda2/envs/py27/bin/retriever

This path can be used for ``rdataretriever`` after removing `retriever` as follows:

::

  library(rdataretriever)
  use_RetrieverPath("/home/<system_name>/anaconda2/envs/py27/bin/")

Note: ``rdataretriever`` will be able to locate ``retriever`` even if the virtual environment is
deactivated.

rdataretriever functions:
=========================

datasets()
^^^^^^^^^^
**Description** : The function returns a list of available datasets.

**Arguments** : No arguments needed.

**Example** :

::

  rdataretriever::datasets()

fetch()
^^^^^^^
**Description** : Each datafile in a given dataset is downloaded to a temporary directory and then imported as a
data.frame as a member of a named list.

**Arguments** :

- ``dataset`` (String): Name of dataset to be downloaded
- ``quiet`` (Bool): The argument decides if warnings need to be displayed (TRUE/FALSE)
- ``data_name`` (String): Name assigned to dataset once it is downloaded

**Example** :

::

  rdataretriever :: fetch(dataset = 'portal')

download()
^^^^^^^^^^
**Description** : Used to download datasets directly without cleaning them and when user does not
have a specific preference for the format of the data and the kind of database.


**Arguments** :

- ``dataset`` (String): Name of the dataset to be downloaded.

- ``path`` (String): Specify dataset download path.

- ``quiet``  (Bool): Setting TRUE minimizes the console output.

- ``sub_dir`` (String): sub_dir downloaded dataset is stored into a custom subdirectory.

- ``debug``  (Bool): Setting TRUE helps in debugging in case of errors.

- ``use_cache``  (Bool): Setting FALSE reinstalls scripts even if they are already installed.

**Example** :

::

  rdataretriever :: download("iris","/Users/username/Desktop")

Installation functions
^^^^^^^^^^^^^^^^^^^^^^
Format specific installation
----------------------------
**Description** : ``rdataretriever`` supports installation of datasets in three file formats through different functions:

- csv (``install_csv``)
- json (``install_json``)
- xml (``install_xml``)

**Arguments** : These functions require same arguments.

- ``dataset`` (String): Name of the dataset to install.

- ``table_name`` (String): Specify the table name to install.

- ``data_dir`` (String): Specify the dir path to store data, defaults to working dir

- ``debug`` (Bool): Setting TRUE helps in debugging in case of errors.

- ``use_cache`` (Bool): Setting FALSE reinstalls scripts even if they are already installed.

**Example** :

::

  rdataretriever :: install_csv("bird-size",table_name = "Bird_Size",debug = TRUE)

Database specific installation
------------------------------
**Description** : ``rdataretriever`` supports installation of datasets in four different databses through different functions:

- MySQL (``install_mysql``)
- PostgreSQL (``install_postgres``)
- SQLite (``install_sqlite``)
- MSAccess (``install_msaccess``)

**Arguments for PostgreSQL and MySQL** :

- ``database_name`` (String): Specify database name.

- ``debug``           (Bool): Setting True helps in debugging in case of errors.

- ``host``          (String): Specify host name for database.

- ``password``      (String): Specify password for database.

- ``port``             (Int): Specify the port number for installation.

- ``quiet``           (Bool): Setting True minimizes the console output.

- ``table_name``    (String): Specify the table name to install.

- ``use_cache``       (Bool): Setting False reinstalls scripts even if they are already installed.

- ``user``          (String): Specify the username.

**Example** :

::

  rdataretriever :: install_postgres(dataset = 'portal', user='postgres', password='abcdef')

**Arguments for MSAccess and SQLite** :

- ``file`` (String): Enter file_name for database.

- ``table_name`` (String): Specify the table name to install.

- ``debug`` (Bool): Setting True helps in debugging in case of errors.

- ``use_cache`` (Bool): Setting False reinstalls scripts even if they are already installed.

**Example** :

::

  rdataretriever :: install_sqlite(dataset = 'iris', file = 'sqlite.db',debug=FALSE, use_cache=TRUE)

get_updates()
^^^^^^^^^^^^^
**Description** : This function will check if the version of the retriever’s scripts in your local directory ‘
~/.retriever/scripts/' is up-to-date with the most recent official retriever release.

**Example** :

::

  rdataretriever :: get_updates()

reset()
^^^^^^^
**Description** : The function will Reset the components of rdataretriever using scope [ all, scripts, data, connection]

**Arguments** :

- ``scope`` : Specifies what components to reset.  Options include:  ’scripts’, ’data’, ’connection’ and ’all’, where ’all’ is the default setting that resets all components.

**Example** :

::

  rdataretriever :: reset(scope = 'data')


Examples
========

::

 library(rdataretriever)

 # List the datasets available via the retriever
 rdataretriever::datasets()

 # Install the Gentry forest transects dataset into csv files in your working directory
 rdataretriever::install('gentry-forest-transects', 'csv')

 # Download the raw Gentry dataset files without any processing to the
 # subdirectory named data
 rdataretriever::download('gentry-forest-transects', './data/')

 # Install and load a dataset as a list
 Gentry = rdataretriever::fetch('gentry-forest-transects')
 names(gentry-forest-transects)
 head(gentry-forest-transects$counts)


To get citation information for the ``rdataretriever`` in R use ``citation(package = 'rdataretriever')``:


.. _Data Retriever: https://github.com/weecology/retriever
.. _rdataretriever: https://github.com/ropensci/rdataretriever
====================
Retriever Provenance
====================

Retriever allows committing of datasets and installation of the committed dataset into the database of your choice at a
later date.
This ensures that the previous outputs/results can be produced easily.

Provenance Directory
====================
The directory to save your committed dataset can be defined by setting the environment variable ``PROVENANCE_DIR``.
However, you can still save the committed dataset in a directory of your choice by defining the ``path`` while committing
the dataset.

Commit Datasets
===============

Retriever supports committing of a dataset into a compressed archive.


.. code-block:: python

  def commit(dataset, commit_message='', path=None, quiet=False):

A description of the default parameters mentioned above:

.. code-block:: python

  dataset               (String): Name of the dataset.

  commit_message        (String): Specify commit message for a commit.

  path                  (String): Specify the directory path to store the compressed archive file.

  quiet                   (Bool): Setting True minimizes the console output.


Example to commit dataset:

.. code-block:: bash

   retriever commit abalone-age -m "Example commit" --path .
   Committing dataset abalone-age
   Successfully committed.


.. code-block:: python

  >>> from retriever import commit
  >>> commit('abalone-age', commit_message='Example commit', path='/home/')

If the path is not provided the committed dataset is saved in the ``provenance directory``.

Log Of Committed Datasets
=========================
You can view the log of commits of the datasets stored in the provenance directory.

.. code-block:: python

  def commit_log(dataset):

A description of the parameter mentioned above:

.. code-block:: python

  dataset       (String): Name of the dataset.

Example:

.. code-block:: bash

  retriever log abalone-age

  Commit message: Example commit
  Hash: 02ee77
  Date: 08/16/2019, 16:12:28


.. code-block:: python

  >>> from retriever import commit_log
  >>> commit_log('abalone-age')


Installing Committed Dataset
============================
You can install committed datasets by using the hash-value or by providing the path of the compressed archive.
Installation using hash-value is supported only for datasets stored in the provenance directory.

For installing dataset from a committed archive you can provide the path to the archive in place of dataset name:

.. code-block:: bash

  retriever install sqlite abalone-age-02ee77.zip

.. code-block:: python

  >>> from retriever import install_sqlite
  >>> install_sqlite('abalone-age-02ee77.zip')

Also, you can install using the hash-value of the datasets stored in provenance directory. You can always look up the
hash-value of your previous commits using the command ``retriever log dataset_name``.

For installing dataset from provenance directory provide the ``hash-value`` of the commit.

.. code-block:: bash

  retriever install sqlite abalone-age --hash-value 02ee77

.. code-block:: python

  >>> from retriever import install_sqlite
  >>> install_sqlite('abalone-age', hash_value='02ee77')
retriever package
=================

Subpackages
-----------

.. toctree::

    retriever.engines
    retriever.lib
===========
Quick Start
===========

The Data Retriever is written in Python and has a Python interface, a command
line interface or an associated R package. It installs publicly available data
into a variety of databases (MySQL, PostgreSQL, SQLite) and flat file formats
(csv, json, xml).

Installation
~~~~~~~~~~~~

Using conda:

$ ``conda install retriever -c conda-forge``

or pip:

$ ``pip install retriever``

To install the associated R package:

$ ``install.packages('rdataretriever')``

Python interface
~~~~~~~~~~~~~~~~

Import:

$ ``import retriever as rt``

List available datasets:

$ ``rt.dataset_names()``

Load data on GDP from the World bank:

$ ``rt.fetch('gdp')``

Install the World Bank data on GDP into an SQLite databased named
"gdp.sqlite":

$ ``rt.install_sqlite('gdp', file='gdp.sqlite)``

Command line interface
~~~~~~~~~~~~~~~~~~~~~~

List available datasets:

$ ``retriever ls``

Install the `Portal dataset <https://github.com/weecology/portaldata>`_ into a
set of json files:

$ ``retriever install json portal``

Install the Portal dataset into an SQLite database named "portal.sqlite":

$ ``retriever install sqlite portal -f portal.sqlite``

R interface
~~~~~~~~~~~

List available datasets:

$ ``rdataretriever::datasets()``

Load data on GDP from the World bank:

$ ``rdataretriever::fetch(dataset = 'gdp')``

Install the GDP dataset into SQLite:

$ ``rdataretriever::install('gdp', 'sqlite')``

Learn more
~~~~~~~~~~

Check out the rest of the documentation for more commands, details, and
datasets.

Available install formats for all interfaces are: mysql, postgres, sqlite,
csv, json, and xml.
=====================
Using the Socrata API
=====================

This tutorial explains the usage of the Socrata API in Data Retriever. It includes both the
CLI (Command Line Interface) commands as well as the Python interface for the same.

.. note::
  Currently Data Retriever only supports tabular Socrata datasets (tabular Socrata datasets which are of type map are not supported).


Command Line Interface
======================

Listing the Socrata Datasets
----------------------------

The ``retriever ls -s`` command displays the Socrata datasets which contain the provided keywords in their title.

$ ``retriever ls -h`` (gives listing options)

::

    usage: retriever ls [-h] [-l L [L ...]] [-k K [K ...]] [-v V [V ...]]
                    [-s S [S ...]]

    optional arguments:
      -h, --help    show this help message and exit
      -l L [L ...]  search datasets with specific license(s)
      -k K [K ...]  search datasets with keyword(s)
      -v V [V ...]  verbose list of specified dataset(s)
      -s S [S ...]  search socrata datasets with name(s)

**Example**

This example will list the names of the socrata datasets which contain the word ``fishing``.

$ ``retriever ls -s fishing``

::

    Autocomplete suggestions : Total 34 results

    [?] Select the dataset name: Recommended Fishing Rivers And Streams
     > Recommended Fishing Rivers And Streams
       Recommended Fishing Rivers And Streams API
       Iowa Fishing Report
       Recommended Fishing Rivers, Streams, Lakes and Ponds Map
       Public Fishing Rights Parking Areas Map
       Fishing Atlas
       Cook County - Fishing Lakes
       [ARCHIVED] Fishing License Sellers
       Public Fishing Rights Parking Areas
       Recommended Fishing Lakes and Ponds Map
       Recommended Fishing Lakes and Ponds
       Delaware Fishing Licenses and Trout Stamps
       Cook County - Fishing Lakes - KML

Here the user is prompted to select a dataset name. After selecting a dataset, the command returns some information related to the dataset selected.

Let's select the ``Public Fishing Rights Parking Areas`` dataset, after pressing Enter, the command returns
some information regarding the dataset selected.

::

    Autocomplete suggestions : Total 34 results

    [?] Select the dataset name: Public Fishing Rights Parking Areas
       Iowa Fishing Report
       Recommended Fishing Rivers, Streams, Lakes and Ponds Map
       Fishing Atlas
       Public Fishing Rights Parking Areas Map
       [ARCHIVED] Fishing License Sellers
       Cook County - Fishing Lakes
     > Public Fishing Rights Parking Areas
       Recommended Fishing Lakes and Ponds Map
       Recommended Fishing Lakes and Ponds
       Delaware Fishing Licenses and Trout Stamps
       Cook County - Fishing Lakes - KML
       General Fishing and Salmon Licence Sales
       Hunting and Fishing License Sellers

    Dataset Information of Public Fishing Rights Parking Areas: Total 1 results

    1. Public Fishing Rights Parking Areas
      ID : 9vef-6whi
      Type : {'dataset': 'tabular'}
      Description : The New York State Department of Environmental Con...
      Domain : data.ny.gov
      Link : https://data.ny.gov/Recreation/Public-Fishing-Rights-Parking-Areas/9vef-6whi


Downloading the Socrata Datasets
--------------------------------

The ``retriever download socrata-<socrata id>`` command downloads the Socrata dataset which matches the provided ``socrata id``.

**Example**

From the example in ``Listing the Socrata Datasets`` section, we selected the *Public Fishing Rights Parking Areas* dataset.
Since the dataset is of type ``tabular``, we can download it. The information received in the previous example contains the ``socrata id``.
We use this ``socrata id`` to download the dataset.

$ ``retriever download socrata-9vef-6whi``

::

    => Installing socrata-9vef-6whi
    Downloading 9vef-6whi.csv: 10.0B [00:03, 2.90B/s]
    Done!

The downloaded raw data files are stored in the ``raw_data`` directory in the ``~/.retriever`` directory.


Installing the Socrata Datasets
-------------------------------

The ``retriever install <engine> socrata-<socrata id>`` command downloads the raw data, creates the script for it and then installs
the Socrata dataset which matches the provided ``socrata id`` into the provided ``engine``.

**Example**

From the example in ``Listing the Socrata Datasets`` section, we selected the *Public Fishing Rights Parking Areas* dataset.
Since the dataset is of type ``tabular``, we can install it. The information received in that section contains the ``socrata id``.
We use this ``socrata id`` to install the dataset.

$ ``retriever install postgres socrata-9vef-6whi``

::

    => Installing socrata-9vef-6whi
    Downloading 9vef-6whi.csv: 10.0B [00:03, 2.69B/s]
    Processing... 9vef-6whi.csv
    Successfully wrote scripts to /home/user/.retriever/socrata-scripts/9vef_6whi.csv.json
    Updating script name to socrata-9vef-6whi.json
    Updating the contents of script socrata-9vef-6whi
    Successfully updated socrata_9vef_6whi.json
    Creating database socrata_9vef_6whi...

    Bulk insert on ..  socrata_9vef_6whi.socrata_9vef_6whi
    Done!

The script created for the Socrata dataset is stored in the ``socrata-scripts`` directory in the ``~/.retriever`` directory.


Python Interface in Data Retriever
==================================

Searching Socrata Datasets
--------------------------

The function ``socrata_autocomplete_search`` takes a list of strings as input and returns a list of strings which are the autocompleted names.

.. code-block:: python

  >>> import retriever as rt
  >>> names = rt.socrata_autocomplete_search(['clinic', '2015', '2016'])
  >>> for name in names:
  ...     print(name)
  ...
  2016 & 2015 Clinic Quality Comparisons for Clinics with Five or More Service Providers
  2015 - 2016 Clinical Quality Comparison (>=5 Providers) by Geography
  2016 & 2015 Clinic Quality Comparisons for Clinics with Fewer than Five Service Providers


Socrata Dataset Info by Dataset Name
------------------------------------

The input argument for the function ``socrata_dataset_info`` should be a string (valid dataset name returned by ``socrata_autocomplete_search``).
It returns a list of dicts, because there are multiple datasets on socrata with same name (e.g. ``Building Permits``).

.. code-block:: python

  >>> import retriever as rt
  >>> resource = rt.socrata_dataset_info('2016 & 2015 Clinic Quality Comparisons for Clinics with Five or More Service Providers')
  >>> from pprint import pprint
  >>> pprint(resource)
  [{'description': 'This data set includes comparative information for clinics '
                   'with five or more physicians for medical claims in 2015 - '
                   '2016. \r\n'
                   '\r\n'
                   'This data set was calculated by the Utah Department of '
                   'Health, Office of Healthcare Statistics (OHCS) using Utah’s '
                   'All Payer Claims Database (APCD).',
    'domain': 'opendata.utah.gov',
    'id': '35s3-nmpm',
    'link': 'https://opendata.utah.gov/Health/2016-2015-Clinic-Quality-Comparisons-for-Clinics-w/35s3-nmpm',
    'name': '2016 & 2015 Clinic Quality Comparisons for Clinics with Five or '
            'More Service Providers',
    'type': {'dataset': 'tabular'}}]


Finding Socrata Dataset by Socrata ID
-------------------------------------

The input argument of the function ``find_socrata_dataset_by_id`` should be the four-by-four socrata dataset identifier (e.g. ``35s3-nmpm``).
The function returns a dict which contains metadata about the dataset.

.. code-block:: python

  >>> import retriever as rt
  >>> from pprint import pprint
  >>> resource = rt.find_socrata_dataset_by_id('35s3-nmpm')
  >>> pprint(resource)
  {'datatype': 'tabular',
   'description': 'This data set includes comparative information for clinics '
                  'with five or more physicians for medical claims in 2015 - '
                  '2016. \r\n'
                  '\r\n'
                  'This data set was calculated by the Utah Department of '
                  'Health, Office of Healthcare Statistics (OHCS) using Utah’s '
                  'All Payer Claims Database (APCD).',
   'domain': 'opendata.utah.gov',
   'homepage': 'https://opendata.utah.gov/Health/2016-2015-Clinic-Quality-Comparisons-for-Clinics-w/35s3-nmpm',
   'id': '35s3-nmpm',
   'keywords': ['socrata'],
   'name': '2016 & 2015 Clinic Quality Comparisons for Clinics with Five or More '
           'Service Providers'}


Downloading a Socrata Dataset
-----------------------------

.. code-block:: python

  import retriever as rt
  rt.download('socrata-35s3-nmpm')


Installing a Socrata Dataset
----------------------------

.. code-block:: python

  import retriever as rt
  rt.install_postgres('socrata-35s3-nmpm')


.. note::

  For downloading or installing the Socrata Datasets, the dataset should follow the syntax given.
  The dataset name should be ``socrata-<socrata id>``. The ``socrata id`` should be the four-by-four
  socrata dataset identifier (e.g. ``35s3-nmpm``).

  Example:
    - Correct: ``socrata-35s3-nmpm``

    - Incorrect:  ``socrata35s3-nmpm``, ``socrata35s3nmpm``
==============
APIs Available
==============

**Socrata API**
---------------

**Total number of datasets supported : 85,244 out of 213,965**

**Rdatasets**
-------------


The list of rdatasets is generated using conf.py.
The file can't be edited on GitHub because it is created in runtime.

Look at the python `conf`_. module

.. _conf: https://github.com/weecology/retriever/blob/main/docs/conf.pyretriever.lib package
=====================

Submodules
----------

retriever.lib.cleanup module
----------------------------

.. automodule:: retriever.lib.cleanup
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.create\_scripts module
------------------------------------

.. automodule:: retriever.lib.create_scripts
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.datapackage module
--------------------------------

.. automodule:: retriever.lib.datapackage
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.datasets module
-----------------------------

.. automodule:: retriever.lib.datasets
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.defaults module
-----------------------------

.. automodule:: retriever.lib.defaults
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.download module
-----------------------------

.. automodule:: retriever.lib.download
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.dummy module
--------------------------

.. automodule:: retriever.lib.dummy
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.engine module
---------------------------

.. automodule:: retriever.lib.engine
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.engine\_tools module
----------------------------------

.. automodule:: retriever.lib.engine_tools
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.excel module
--------------------------

.. automodule:: retriever.lib.excel
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.fetch module
--------------------------

.. automodule:: retriever.lib.fetch
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.get\_opts module
------------------------------

.. automodule:: retriever.lib.get_opts
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.install module
----------------------------

.. automodule:: retriever.lib.install
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.load\_json module
-------------------------------

.. automodule:: retriever.lib.load_json
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.models module
---------------------------

.. automodule:: retriever.lib.models
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.provenance module
-------------------------------

.. automodule:: retriever.lib.provenance
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.provenance\_tools module
--------------------------------------

.. automodule:: retriever.lib.provenance_tools
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.rdatasets module
------------------------------

.. automodule:: retriever.lib.rdatasets
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.repository module
-------------------------------

.. automodule:: retriever.lib.repository
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.scripts module
----------------------------

.. automodule:: retriever.lib.scripts
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.socrata module
----------------------------

.. automodule:: retriever.lib.socrata
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.table module
--------------------------

.. automodule:: retriever.lib.table
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.templates module
------------------------------

.. automodule:: retriever.lib.templates
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.tools module
--------------------------

.. automodule:: retriever.lib.tools
    :members:
    :undoc-members:
    :show-inheritance:

retriever.lib.warning module
----------------------------

.. automodule:: retriever.lib.warning
    :members:
    :undoc-members:
    :show-inheritance:
retriever.engines package
=========================

Submodules
----------

retriever.engines.csvengine module
----------------------------------

.. automodule:: retriever.engines.csvengine
    :members:
    :undoc-members:
    :show-inheritance:

retriever.engines.download\_only module
---------------------------------------

.. automodule:: retriever.engines.download_only
    :members:
    :undoc-members:
    :show-inheritance:

retriever.engines.jsonengine module
-----------------------------------

.. automodule:: retriever.engines.jsonengine
    :members:
    :undoc-members:
    :show-inheritance:

retriever.engines.msaccess module
---------------------------------

.. automodule:: retriever.engines.msaccess
    :members:
    :undoc-members:
    :show-inheritance:

retriever.engines.mysql module
------------------------------

.. automodule:: retriever.engines.mysql
    :members:
    :undoc-members:
    :show-inheritance:

retriever.engines.postgres module
---------------------------------

.. automodule:: retriever.engines.postgres
    :members:
    :undoc-members:
    :show-inheritance:

retriever.engines.sqlite module
-------------------------------

.. automodule:: retriever.engines.sqlite
    :members:
    :undoc-members:
    :show-inheritance:

retriever.engines.xmlengine module
----------------------------------

.. automodule:: retriever.engines.xmlengine
    :members:
    :undoc-members:
    :show-inheritance:

Retriever API
=============

.. toctree::
   :maxdepth: 4

   retriever
.. retriever documentation master file, created by
   sphinx-quickstart on Tue Jan 12 14:27:41 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Data retriever's documentation!
=============================================

Contents:

.. toctree::
   :maxdepth: 2

   introduction
   quickstart
   rdataretriever.rst
   recipes.rst
   pyretriever.rst
   provenance.rst
   datasets_list.rst
   available_apis.rst
   socrata_api.rst
   rdatasets_api.rst
   developer.rst
   spatial_dbms.rst
   code_of_conduct.rst
   modules.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

============
User Guide
============


We handle the data so you can focus on the science
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finding data is one thing. Getting it ready for analysis is another. Acquiring,
cleaning, standardizing and importing publicly available data is time consuming
because many datasets lack machine readable metadata and do not conform to
established data structures and formats.

The Data Retriever automates the first
steps in the data analysis pipeline by downloading, cleaning, and standardizing
datasets, and importing them into relational databases, flat files, or
programming languages. The automation of this process reduces the time for a
user to get most large datasets up and running by hours, and in some cases days.


What data tasks does the Retriever handle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Data Retriever handles a number of common tasks including:
 #. Creating the underlying database structures, including automatically determining the data types
 #. Downloading the data
 #. Transforming data into appropriately normalized forms for database management systems (e.g., "wide" data into "long" data and splitting tables into proper sub-tables to reduce duplication)
 #. Converting heterogeneous null values (e.g., 999.0, -999, NaN) into standard null values
 #. Combining multiple data files into single tables; and 6) placing all related tables in a single database or schema.

A couple of examples on the more complicated end include the Breeding Bird
Survey of North America (breed-bird-survey) and the Alwyn Gentry Tree Transect
data(gentry-forest-transects):

- Breeding bird survey data consists of multiple tables. The main table is divided
  into one file per region in 70 individual compressed files. Supplemental tables
  required to work with the data are posted in a variety of locations and formats.
  The Data Retriever automates: downloading all data files, extracting data from
  region-specific raw data files into single tables, correcting typographic
  errors, replacing non-standard null values, and adding a Species table that
  links numeric identifiers to actual species names.
- The Gentry forest transects data is stored in over 200 Excel spreadsheets, each
  representing an individual study site, and compressed in a zip archive.
  Each spreadsheet contains counts of individuals found at a given site and all stems
  measured from that individual; each stem measurement is placed in a separate column,
  resulting in variable numbers of columns across rows, a format that is
  difficult to work with in both database and analysis software. There is no
  information on the site in the data files themselves, it is only present in
  the names of the files. The Retriever downloads the archive, extracts the
  files, and splits the data they contain into four tables: Sites, Species,
  Stems, and Counts, keeping track of which file each row of count data
  originated from in the Counts table and placing a single stem on each row in
  the Stems table.

*Adapted from* `Morris & White 2013`_.

Installing (binaries)
~~~~~~~~~~~~~~~~~~~~~

Precompiled binaries of the most recent release are available for Windows,
OS X, and Ubuntu/Debian at the `project website`_.

Installing From Source
~~~~~~~~~~~~~~~~~~~~~~

**Required packages**

To install the Data Retriever from source, you’ll need Python 3.6.8+
with the following packages installed:

-  xlrd

**The following packages are optional**

-  PyMySQL (for MySQL)
-  sqlite3 (for SQLite, v3.8 or higher required)
-  psycopg2-binary (for PostgreSQL)
-  pypyodbc (for MS Access)

**Steps to install from source**

1. Clone the repository
2. From the directory containing setup.py, run the following command:
   ``python setup.py install`` or use pip ``pip install . --upgrade`` to install and
   ``pip uninstall retriever`` to uninstall the retriever

3. After installing, type ``retriever`` from a command prompt to see the available options of
   the Data Retriever. Use ``retriever --version`` to confirm the version installed on your system.

Using the Data Retriever Commands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After installing, run ``retriever update`` to download all of the
available dataset scripts. Run ``retriever ls`` to see the available datasets

To see the full list of command line options
and datasets run ``retriever --help``. The output will look like this:

::

    usage: retriever [-h] [-v] [-q]
                     {download,install,defaults,update,new,new_json,edit_json,delete_json,ls,citation,reset,help,commit}
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
        commit              commit dataset to a zipped file
        log                 see log of a committed dataset

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      -q, --quiet           suppress command-line output

To install datasets, use the ``install`` command.

Examples
~~~~~~~~

**Using install**

The install command downloads the datasets and installs them in the desired engine.

$ ``retriever install -h`` (gives install options)

::

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


**Examples using install**


These examples use Breeding Bird Survey data (breed-bird-survey).
The retriever has support for various databases and flat file
formats (mysql, postgres, sqlite, msaccess, csv, json, xml).
All the engines have a variety of options or flags. Run ```retriever defaults`` to see the defaults.
For example, the default options for mysql and postgres engines are given below.

::

    retriever defaults

    Default options for engine  MySQL
    user   root
    password
    host   localhost
    port   3306
    database_name   {db}
    table_name   {db}.{table}

    Default options for engine  PostgreSQL
    user   postgres
    password
    host   localhost
    port   5432
    database   postgres
    database_name   {db}
    table_name   {db}.{table}

Help information for a particular engine can be obtained by running
retriever install [engine name] [-h] [--help], for example, ``retriever install mysql -h``.
Both mysql and postgres require the database user name ``--user [USER], -u [USER]``
and password ``--password [PASSWORD], -p [PASSWORD]``.
MySQL and PostgreSQL database management systems support the use of configuration files.
The configuration files provide a mechanism to support using the engines without providing authentication directly.
To set up the configuration files please refer to the respective database management systems documentation.

Install data into Mysql::

   retriever install mysql –-user myusername –-password ***** –-host localhost –-port 8888 –-database_name testdbase breed-bird-survey
   retriever install mysql –-user myusername breed-bird-survey (using attributes in the client authentication configuration file)

Install data into postgres::

   retriever install postgres –-user myusername –-password ***** –-host localhost –-port 5432 –-database_name testdbase breed-bird-survey
   retriever install postgres breed-bird-survey (using attributes in the client authentication configuration file)

Install data into sqlite::

   retriever install sqlite breed-bird-survey -f mydatabase.db (will use mydatabase.db)
   retriever install sqlite breed-bird-survey (will use or create default sqlite.db in working directory)

Install data into csv::

   retriever install csv breed-bird-survey --table_name  "BBS_{table}.csv"
   retriever install csv breed-bird-survey

**Using download**

The ``download`` command downloads the raw data files exactly as they occur at the
source without any clean up or modification. By default the files will be stored in the working directory.

``--path`` can be used to specify a location other than the working directory to download the files to. E.g., ``--path ./data``

``--subdir`` can be used to maintain any subdirectory structure that is present in the files being downloaded.

::

   retriever download -h (gives you help options)
   retriever download breed-bird-survey (download raw data files to the working directory)
   retriever download breed-bird-survey –path  C:\Users\Documents (download raw data files to path)


**Using citation**

The ``citation`` command show the citation for the retriever and for the scripts.

::

   retriever citation (citation of the Data retriever)
   retriever citation breed-bird-survey (citation of Breed bird survey data)

**To create new, edit, delete scripts please read the documentation on scripts**


Storing database connection details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The retriever reads from the standard configuration files for the database
management systems. If you want to store connection details they should be
stored in those files. Make sure to secure these files appropriately.

For postgreSQL, create or modify `~/.pgpass`. This is a file named `.pgpass`
located in the users home directory. On Microsoft Windows, the file is named
`%APPDATA%\postgresql\pgpass.conf` (where `%APPDATA%` refers to the Application
Data subdirectory in the user's profile). It should take the general form:

``hostname:port:database:username:password``

where each word is replaced with the correct information for your database
connection or replaced with an ``*`` to apply to all values for that section.

For MySQL, create or modify `~/.my.cnf`. This is a file named `.my.cnf` located
in the users home directory. The relevant portion of this file for the retriever
is the `client` section which should take the general form:

::

   [client]
   host=hostname
   port=port
   user=username
   password=password

where each word to the right of the `=` is replaced with the correct information
for your database connection. Remove or comment out the lines for any values you
don't want to set.


Acknowledgments
~~~~~~~~~~~~~~~

Development of this software was funded by `the Gordon and Betty Moore
Foundation’s Data-Driven Discovery Initiative`_ through `Grant
GBMF4563`_ to Ethan White and the `National Science Foundation`_ as part
of a `CAREER award to Ethan White`_.


.. _the Gordon and Betty Moore Foundation’s Data-Driven Discovery Initiative: http://www.moore.org/programs/science/data-driven-discovery
.. _Grant GBMF4563: http://www.moore.org/grants/list/GBMF4563
.. _National Science Foundation: http://nsf.gov/
.. _CAREER award to Ethan White: http://nsf.gov/awardsearch/showAward.do?AwardNumber=0953694
.. _project website: http://data-retriever.org
.. _Morris & White 2013: https://dx.doi.org/10.1371/journal.pone.0065848
==================
Datasets Available
==================


The list of datasets is generated using conf.py.
The file can't be edited on GitHub because it is created in runtime.

Look at the python `conf`_. module

.. _conf: https://github.com/weecology/retriever/blob/main/docs/conf.py
Contributor Code of Conduct
===========================


Our Pledge
^^^^^^^^^^

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

Our Standards
^^^^^^^^^^^^^

Examples of behavior that contributes to creating a positive environment include:

Using welcoming and inclusive language
Being respectful of differing viewpoints and experiences
Gracefully accepting constructive criticism
Focusing on what is best for the community
Showing empathy towards other community members
Examples of unacceptable behavior by participants include:

The use of sexualized language or imagery and unwelcome sexual attention or advances
Trolling, insulting/derogatory comments, and personal or political attacks
Public or private harassment
Publishing others' private information, such as a physical or electronic address, without explicit permission
Other conduct which could reasonably be considered inappropriate in a professional setting

Our Responsibilities
^^^^^^^^^^^^^^^^^^^^

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

Scope
^^^^^

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

Enforcement
^^^^^^^^^^^

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at `ethan@weecology.org
<mailto:ethan@weecology.org>`__. All complaints will be reviewed and investigated and will result in a response that is deemed necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

Attribution
^^^^^^^^^^^

This Code of Conduct is adapted from the  `Contributor Covenant, version 1.4`_.


.. _Contributor Covenant, version 1.4: http://contributor-covenant.org/version/1/4
