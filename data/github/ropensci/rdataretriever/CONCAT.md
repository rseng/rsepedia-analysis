# rdataretriever

<!-- badges: start -->
[![Build Status](https://travis-ci.org/ropensci/rdataretriever.svg?branch=main)](https://travis-ci.org/ropensci/rdataretriever)
[![Build status](https://ci.appveyor.com/api/projects/status/de1badmnrt6goamh?svg=true)](https://ci.appveyor.com/project/ethanwhite/rdataretriever)[![cran version](https://www.r-pkg.org/badges/version/rdataretriever)](https://CRAN.R-project.org/package=rdataretriever)
[![Documentation Status](https://readthedocs.org/projects/retriever/badge/?version=latest)](https://retriever.readthedocs.io/en/latest/rdataretriever.html#)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rdataretriever)](https://CRAN.R-project.org/package=rdataretriever) +
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ecoretriever)](https://CRAN.R-project.org/package=ecoretriever)
(old package name)
[![DOI](https://zenodo.org/badge/18155090.svg)](https://zenodo.org/badge/latestdoi/18155090)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02800/status.svg)](https://doi.org/10.21105/joss.02800)

<!-- badges: end -->

R interface to the [Data Retriever](https://retriever.readthedocs.io/en/latest/).

The `rdataretriever` provides access to cleaned versions of hundreds of commonly used public datasets with a single line of code. 

These datasets come from many different sources and most of them require some cleaning and restructuring prior to analysis.
The `rdataretriever` uses a set of actively maintained recipes for downloading, cleaning, and restructing these datasets using a combination of the [Frictionless Data Specification](https://specs.frictionlessdata.io/) and custom data cleaning scripts. 

The `rdataretriever` also facilitates the automatic storage of these datasets in a choice of database management systems (PostgreSQL, SQLite, MySQL, MariaDB) or flat file formats (CSV, XML, JSON) for later use and integration with large data analysis pipelines.

The `rdatretriever` also facilitates reproducibile science by providing tools to archive and rerun the precise version of a dataset and associated cleaning steps that was used for a specific analysis.

The `rdataretriever` handles the work of cleaning, storing, and archiving data so that you can focus on analysis, inference and visualization.


## Table of Contents

  - [Installation](#installation)
      - [Basic Installation (no Python experience needed)](#basic-installation)
      - [Advanced Installation for Python Users](#advanced-installation-for-python-users)
  - [Installing Tabular Datasets](#installing-tabular-datasets)
  - [Installing Spatial Datasets](#installing-spatial-datasets)
  - [Using Docker Containers](#using-docker-containers)
  - [Provenance](#provenance)
  - [Acknowledgements](#acknowledgements)

## Installation

The `rdataretriever` is an R wrapper for the Python package, [Data Retriever](https://retriever.readthedocs.io/en/latest/). This means
that *Python* and the `retriever` Python package need to be installed first.

### Basic Installation

If you just want to use the Data Retriever from within R follow these
instuctions run the following commands in R. This will create a local Python
installation that will only be used by R and install the needed Python package
for you.

```coffee
install.packages('reticulate') # Install R package for interacting with Python
reticulate::install_miniconda() # Install Python
reticulate::py_install('retriever') # Install the Python retriever package
install.packages('rdataretriever') # Install the R package for running the retriever
rdataretriever::get_updates() # Update the available datasets
```

**After running these commands restart R.**

### Advanced Installation for Python Users

If you are using Python for other tasks you can use `rdataretriever` with your
existing Python installation (though the [basic installation](#basic-installation)
above will also work in this case by creating a separate miniconda install and
Python environment).

#### Install the `retriever` Python package

Install the `retriever` Python package into your prefered Python environment
using either `conda` (64-bit conda is required):

  ```bash
  conda install -c conda-forge retriever
  ```

  or `pip`:

  ```bash
  pip install retriever
  ```

#### Select the Python environment to use in R

`rdataretriever` will try to find Python environments with `retriever` (see the
`reticulate` documentation on
[order of discovery](https://rstudio.github.io/reticulate/articles/versions.html#order-of-discovery-1)
for more details) installed. Alternatively you can select a Python environment
to use when working with `rdataretriever` (and other packages using
`reticulate`).

The most robust way to do this is to set the `RETICULATE_PYTHON` environment
variable to point to the preferred Python executable:

```coffee
Sys.setenv(RETICULATE_PYTHON = "/path/to/python")
```

This command can be run interactively or placed in `.Renviron` in your home
directory.

Alternatively you can do select the Python environment through the `reticulate`
package for either `conda`:

```coffee
library(reticulate)
use_conda('name_of_conda_environment')
```

or `virtualenv`:

```coffee
library(reticulate)
use_virtualenv("path_to_virtualenv_environment")
```

You can check to see which Python environment is being used with:

```coffee
py_config()
```

#### Install the `rdataretriever` R package

```coffee
install.packages("rdataretriever") # latest release from CRAN
```

```coffee
remotes::install_github("ropensci/rdataretriever") # development version from GitHub
```

## Installing Tabular Datasets

```coffee
library(rdataretriever)

# List the datasets available via the Retriever
rdataretriever::datasets()

# Install the portal into csv files in your working directory
rdataretriever::install_csv('portal')

# Download the raw portal dataset files without any processing to the
# subdirectory named data
rdataretriever::download('portal', './data/')

# Install and load a dataset as a list
portal = rdataretriever::fetch('portal')
names(portal)
head(portal$species)

```

## Installing Spatial Datasets

**Set-up and Requirements**

**Tools**

-  PostgreSQL with PostGis, psql(client), raster2pgsql, shp2pgsql, gdal,

The `rdataretriever` supports installation of spatial data into `Postgres DBMS`.

1. **Install PostgreSQL and PostGis**

	To install `PostgreSQL` with `PostGis` for use with spatial data please refer to the
	[OSGeo Postgres installation instructions](https://trac.osgeo.org/postgis/wiki/UsersWikiPostGIS21UbuntuPGSQL93Apt).

	We recommend storing your PostgreSQL login information in a `.pgpass` file to
	avoid supplying the password every time.
	See the [`.pgpass` documentation](https://wiki.postgresql.org/wiki/Pgpass) for more details.

	After installation, Make sure you have the paths to these tools added to your system's `PATHS`.
	Please consult an operating system expert for help on how to change or add the `PATH` variables.

	**For example, this could be a sample of paths exported on Mac:**

	```shell
	#~/.bash_profile file, Postgres PATHS and tools.
	export PATH="/Applications/Postgres.app/Contents/MacOS/bin:${PATH}"
	export PATH="$PATH:/Applications/Postgres.app/Contents/Versions/10/bin"

	```

2. **Enable PostGIS extensions**

	If you have `Postgres` set up, enable `PostGIS` extensions.
	This is done by using either `Postgres CLI` or `GUI(PgAdmin)` and run

	**For psql CLI**
	```shell
	psql -d yourdatabase -c "CREATE EXTENSION postgis;"
	psql -d yourdatabase -c "CREATE EXTENSION postgis_topology;"
	```

	**For GUI(PgAdmin)**

	```sql
	CREATE EXTENSION postgis;
	CREATE EXTENSION postgis_topology
	```
	For more details refer to the
	[PostGIS docs](https://postgis.net/docs/postgis_installation.html#install_short_version).

**Sample commands**

```R
rdataretriever::install_postgres('harvard-forest') # Vector data
rdataretriever::install_postgres('bioclim') # Raster data

# Install only the data of USGS elevation in the given extent
rdataretriever::install_postgres('usgs-elevation', list(-94.98704597353938, 39.027001800158615, -94.3599408119917, 40.69577051867074))

```


## Provenance

To ensure reproducibility the `rdataretriever` supports creating snapshots of the data and the script in time.

Use the commit function to create and store the snapshot image of the data in time. Provide a descriptive message for the created commit. This is comparable to a git commit, however the function bundles the data and scripts used as a backup.

With provenace, you will be able to reproduce the same analysis in the future.

**Commit a dataset**

By default commits will be stored in the provenance directory `.retriever_provenance`, but this directory can be changed by setting the environment variable `PROVENANCE_DIR`.

```coffee
rdataretriever::commit('abalone-age',
                       commit_message='A snapshot of Abalone Dataset as of 2020-02-26')
```

You can also set the path for an individual commit:

```coffee
rdataretriever::commit('abalone-age',
                       commit_message='Data and recipe archive for Abalone Data on 2020-02-26',
					   path='.')
```

**View a log of committed datasets in the provenance directory**

```coffee
rdataretriever::commit_log('abalone-age')
```

**Install a committed dataset**

To reanalyze a committed dataset, rdataretriever will obtain the data and script from the history and rdataretriever will install this particular data into the given back-end. For example, SQLite:

```coffee
rdataretriever::install_sqlite('abalone-age-a76e77.zip') 
```
Datasets stored in provenance directory can be installed directly using hash value
```coffee
rdataretriever::install_sqlite('abalone-age', hash_value='a76e77')
``` 

## Using Docker Containers

To run the image interactively

`docker-compose run --service-ports rdata /bin/bash`

To run tests

`docker-compose  run rdata Rscript load_and_test.R`

Release
-------

Make sure you have tests passing on R-oldrelease, current R-release and R-devel

To check the package

```Shell
R CMD Build #build the package
R CMD check  --as-cran --no-manual rdataretriever_[version]tar.gz
```

To Test

```R
setwd("./rdataretriever") # Set working directory
# install all deps
# install.packages("reticulate")
library(DBI)
library(RPostgreSQL)
library(RSQLite)
library(reticulate)
library(RMariaDB)
install.packages(".", repos = NULL, type="source")
roxygen2::roxygenise()
devtools::test()
```

To get citation information for the `rdataretriever` in R use `citation(package = 'rdataretriever')`

Acknowledgements
----------------
A big thanks to Ben Morris for helping to develop the Data Retriever.
Thanks to the rOpenSci team with special thanks to Gavin Simpson,
Scott Chamberlain, and Karthik Ram who gave helpful advice and fostered
the development of this R package.
Development of this software was funded by the [National Science Foundation](https://nsf.gov/)
as part of a [CAREER award to Ethan White](https://nsf.gov/awardsearch/showAward.do?AwardNumber=0953694).

---
[![ropensci footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

rdataretriever 3.0.0
====================

### New Coauthors

We welcome Apoorva Pandey and Hao Ye as coauthors on the package.
Thanks for the great contributions!

### New Features

* Add provenance support
* Add support for new online/offline dataset feature in retriever
* Automatically update dataset status
* Check that the retriever Python package is sufficiently up-to-date
* Update and simplify installation instructions
* Add get citation function
* Update and simplify internals to benefit from updates to reticulate

## Minor Improvements

* Improvements to code style and documentation
* Use Cloud CRAN mirror in Docker-based testing for stability
* Add get_version function, returns the version of the Data Retriever

rdataretriever 2.0.0
====================

### New Coauthors

We welcome Harshit Bansal as a coauthor on the package. 
Thanks for the great help!

### New Features

* Add a customize installation to directory using data dir
* Use Dockers for testing
* Add spatial support using Postgis
* Using reticulate
* Add get_citation function
* Update Reset retriever to include reseting specific scripts or data

### Minor Improvements

* Improve the test platform and use reticulate in the tests
* Test using custom service names specific to project
* Add potential path for retriever in Windows

### Bug Fixes
* No scripts when using reticulate based install on retriever


rdataretriever 1.1.0
====================

NEW COAUTHORS

* We welcome Pranita Sharma and David J. Harris as coauthors on the package. Thanks
for the great help Pranita and David!

NEW FEATURES

* Use reticulate to integrate with retriever Python
* Spatial dataset processing Beta version using PostGis

MINOR IMPROVEMENTS

* Use Docker compose for testing
* Switch to Python 3 for testing
* Enhance path search on Windows
* Add a debug script for reporting the environment variables
* Update documentation to match latest version
* Add license file

rdataretriever 1.0.0
====================

NEW PACKAGE NAME

* The `EcoData Retriever` has been renamed to `Data Retriever` to reflect its
utility outside of ecological data and consequently we have renamed the R package
from `ecoretriever` to `rdataretriever`

NEW COAUTHORS

* We welcome Henry Senyondo and Shawn Taylor as coauthors on the package. Thanks
for the great help Henry and Shawn!

NEW FEATURES

* Add `reset` which allows a user to delete all the `Data Retriever` downloaded
files
* Add `json` and `xml` as output options

MINOR IMPROVEMENTS

* Accommodate new retriever naming conventions in fetch
* Don't change the class or return the update log
* Specify in documentation which functions are for internal use.
* Change dataset names in source and README.md

BUG FIXES

* Search for Anaconda installs of the `Data Retriever`
* Obtain correct home path in RStudio on Windows

ecoretriever 0.3.0
==================

MINOR IMPROVEMENTS

* Improve documentation for using the connection file

BUG FIXES

* Fix issues with running on some Windows machines by using `shell()` instead of
  `system()` on Windows
* Fix new `--subdir` functionality (released in 0.2.2)


ecoretriever 0.2
================

NEW FEATURES
* We added a new function `get_updates` which can be used to update the `retriever` scripts. This is a big improvement for users because it avoids automatically updating the scripts every time the package is imported. The log of the scripts update can be printed in a cleaner format as well.
* Added support for maintaining subdirectory structure when using the function `download`.

MINOR IMPROVEMENTS

* default data_dir argument is now set to working directory rather than NULL for the function `install()`

BUG FIXES

* On windows machine if the data directory was not specified for a dataset install an error would occur. Now the dataset directory is always specified in external calls to `retriever install ...`

# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at rdataretriever@weecology.org. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [https://contributor-covenant.org/version/1/4][version]

[homepage]: https://contributor-covenant.org
[version]: https://contributor-covenant.org/version/1/4/
---
title: 'Rdataretriever: R Interface to the Data Retriever'
tags:
  - data retrieval, data processing, R, data, data science, datasets
authors:
 - name: Henry Senyondo
   orcid: 0000-0001-7105-5808
   affiliation: 1
 - name: Daniel J. McGlinn
   orcid: 0000-0003-2359-3526
   affiliation: 2
 - name: Pranita Sharma
   affiliation: 3
   orcid: 0000-0002-5871-380X
 - name: David J. Harris
   affiliation: 1
   orcid: 0000-0003-3332-9307
 - name: Hao Ye
   affiliation: 4
   orcid: 0000-0002-8630-1458
 - name: Shawn D. Taylor
   affiliation: "1, 5"
   orcid: 0000-0002-6178-6903
 - name: Jeroen Ooms
   affiliation: 6
   orcid: 0000-0002-4035-0289
 - name: Francisco Rodríguez-Sánchez
   affiliation: 7
   orcid: 0000-0002-7981-1599
 - name: Karthik Ram
   affiliation: 6
   orcid: 0000-0002-0233-1757
 - name: Apoorva Pandey
   affiliation: 8
   orcid: 0000-0001-9834-4415
 - name: Harshit Bansal
   affiliation: 9
   orcid: 0000-0002-3285-812X
 - name: Max Pohlman
   affiliation: 10
 - name: Ethan P. White
   affiliation: "1, 11, 12"
   orcid: 0000-0001-6728-7745


affiliations:
 - name: Department of Wildlife Ecology and Conservation, University of Florida
   index: 1
 - name: Department of Biology, College of Charleston
   index: 2
 - name: North Carolina State University, Department of Computer Science
   index: 3
 - name: Health Science Center Libraries, University of Florida
   index: 4
 - name: USDA-ARS Jornada Experimental Range
   index: 5
 - name: Berkeley Institute for Data Science, University of California, Berkeley
   index: 6
 - name: Department of Agricultural Economics, Sociology, and Education, Penn State University
   index: 7
 - name: Department of Electronics and Communication, Indian Institute of Technology, Roorkee
   index: 8
 - name: Ajay Kumar Garg Engineering College, Ghaziabad
   index: 9
 - name: Departamento de Biología Vegetal y Ecología, Universidad de Sevilla. 
   index: 10
 - name: Informatics Institute, University of Florida
   index: 11
 - name: Biodiversity Institute, University of Florida
   index: 12

date: 16 September 2020 
bibliography: paper.bib
---

# rdataretriever: An R package for downloading, cleaning, and installing publicly available datasets

## Summary

The rdataretriever provides an R interface to the Python-based Data Retriever software. The Data Retriever automates the multiple steps of data analysis including downloading, cleaning, standardizing, and importing datasets into a variety of relational databases and flat file formats. It also supports provenance tracking for these steps of the analysis workflow by allowing datasets to be committed at the time of installation and allowing them to be reinstalled with the same data and processing steps in the future. Finally, it supports the installation of spatial datasets into relational databases with spatial support. The rdataretriever provides an R interface to this functionality and also supports importing of datasets directly into R for immediate analysis. The system also supports the use of custom data processing routines to support complex datasets that require custom data manipulation steps. The Data Retriever and rdataretriever are focused on scientific data applications including a number of widely used, but difficult to work with, datasets in ecology and the environmental sciences.

## Statement of Need

Finding, cleaning, standardizing, and importing data into efficient data structures for modeling and visualization represents a major component of most research workflows. This is a time-consuming process for researchers even when working with relatively simple datasets. For more complex datasets, these steps can be so complex as to prevent domain experts from engaging with the dataset at all. Systems that operate like package managers for scientific data can overcome these barriers, allowing researchers to move quickly to the final steps in the data analysis workflow (visualization and modeling) and allowing domain experts to leverage the most complex data appropriate to their research questions. The rdataretriever allows R users to automatically conduct these early steps of the analysis workflow for over 200 datasets including a number of the most widely used and difficult to work with datasets in the environmental sciences including the North American Breeding Bird Survey and the Forest Inventory and Analysis datasets. This actively facilitates research on important ecological and environmental questions that would otherwise be limited.

## Implementation

The main Data Retriever software is written in Python [@Morris2013], [@Senyondo2017]. The rdataretriever allows R users to access this data processing workflow through a combination of the reticulate package [@reticulate] and custom features developed for working in R [@R]. Because many R users, including the domain researchers most strongly supported by this package, are not familiar with Python and its package management systems, a strong emphasis has been placed on simplifying the installation process for this package so that it can be done entirely from R. Installation requires no direct use of Python or the command line. Detailed documentation has been developed to support users in both installation and use of the software. A Docker-based testing system and associated test suite has also been implemented to ensure that the interoperability of the R package and Python package are maintained, which is challenging due to frequent changes in reticulate and complexities in supporting cross-language functionality across multiple operating systems (Windows, Mac OS, Linux) and R programming environments (terminal-based R and RStudio).

For tabular datasets requiring relatively simple workflows the software uses the JSON based Frictionless Data tabular data metadata package standard [@frictionlessdata_specs]. For more complex data processing workflows, custom Python code is used to process the data into cleaned and standardized formats. Spatial data support is available for PostgreSQL using PostGIS. The information required for handling these datasets is based on a customized version of the Frictionless Data Geo Data schema [@frictionlessdata_specs] that also supports raster datasets.

## Acknowledgements

Development of this software was funded by the Gordon and Betty Moore Foundation's Data-Driven Discovery Initiative through Grant GBMF4563 to Ethan White and the National Science Foundation CAREER Award 0953694 to Ethan White.

## References

