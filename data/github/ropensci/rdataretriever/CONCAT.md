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

---
title: "Provenance & Reproducibility Using the rdataretriever"
output: rmarkdown::html_vignette
bibliography: refs.bibtex
vignette: >
  %\VignetteIndexEntry{Provenance & Reproducibility Using the rdataretriever}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

Datasets that are regularly updated are increasingly common [@yenni2019].
This presents two challenges for reprodubility.
First, if the underlying structure of the dataset changes then previously written code for processing the data will often cease to run properly.
Second, if the version of the data used in a particularly analysis isn't archived then if the data changes it will be difficult to reproduce the original analysis. 

The `retriever` and `rdataretriever` address both of these limitations.
The centrally maintained scripts for processing datasets are updated when datasets change structure and so as long as `rdataretriever::get_updates()` is run before installing the dataset all data code for downloading, cleaning, and installing the data will continue to work.
While the regularly updated data processing recipes ensure that code analyzing the datasets will always continue to run, it is important for reproducibility that we be able to rerun the exact data processing steps on the exact data that was used for the original analysis.
The `rdataretriever` has built in provenance functionality to support this.

To store the data and processing script in their current state we use the `commit()` function to store both components of the data processing in a zip file for future reuse.
This is logically similar to a git commit in that we store the state of the data and the process at a moment in time using a hash.

For example, the `portal-dev` dataset is updated weekly.
If we want to be able rerun our original analysis after the reviews for a paper come back we'll need to store that version of the data.

```{r, eval=FALSE}
rdataretriever::commit('portal-dev', commit_message='Archive Portal data processing for initial submission on 2020-02-26', path = '.')
```

When we want to reanalyze this exact state of the dataset we can load it back into SQLite (or any of the other backends). Use the hash number related to the commit.

```{r eval = FALSE}
rdataretriever::install_sqlite("portal-dev-326d87.zip")
```

## References
---
title: "Using the rdataretriever to quickly analyze Breeding Bird Survey data"
output: rmarkdown::html_vignette
bibliography: refs.bibtex
vignette: >
  %\VignetteIndexEntry{Using the rdataretriever to quickly analyze Breeding Bird Survey data}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

## Introduction

The Breeding Bird Survey of North America (BBS) is a widely used dataset for understanding geospatial variation and dynamic changes in bird communities.
The dataset is a continental scale community science project where thousands of birders count birds at locations across North America.
It has been used in hundreds of research projects including research on biodiversity gradients [@hurlbert2003], ecological forecasting [@harris2018], and bird declines [@rosenberg2019].

However working with the Breeding Bird Survey data can be challenging because it is composed of roughly 100 different files that need to be cleaned and then combined in multiple ways.
These initial phases of the data analysis pipeline can require hours of work for even experienced users to understand the detailed layout of the data and either manually assemble it or write code to combine the data.
The data structure and location also changes regularly meaning that code for this work that is not regularly tested quickly stops working.

This vignette demonstrates how using the `rdataretriever` can eliminate hours of work on data cleaning and restructing allowing researchers to quickly begin addressing interesting scientific questions.
To do this is demonstrates analyzing the BBS data to evaluate correlates of biodiversity (in the form of species richness).

## Load R Packages

We start by loading the necessary packages.
In addition to the `rdataretriever` in this demo we'll also use `DBI`, `RSQLite` and `dplyr` to work with the data, `raster` for working with environmental data, and `ggplot2` for visualization.

```{r, eval=FALSE}
library(rdataretriever)
library(DBI)
library(dplyr)
library(dbplyr)
library(raster)
library(ggplot2)
```

## Install And Connect To The Breeding Bird Survey Data

First we'll update the `rdataretriever` to make sure we have the newest data processing recipes in case something about the structure or location of the dataset has changed.
Centrally updated data recipes reproducible research with this dataset because something changes every year when the newest data is released meaning that any custom code for processing the data stops working.
When this happens the `retriever` recipe is updated and after running `get_updates()` data analysis code continues to run as it always has.

```{r, echo = FALSE, results = "hide", eval=FALSE}
rdataretriever::get_updates()
```

Next install the BBS data into an SQLite database named `bbs.sqlite`.

```{r, echo = FALSE, results = "hide", eval=FALSE}
rdataretriever::install_sqlite('breed-bird-survey', 'bbs.sqlite')
```

We could also load the data straight into R (using `rdataretriever::fetch('breed-bird-survey')`), store it as flat files (CSV, JSON, or XML), or load it into other database management systems (PostgreSQL, MariaDB, MySQL).
The data is moderately large (~1GB) so SQLite represents a nice compromise between efficiently conducting the first steps in the data manipulation pipeline while requiring no additional setup or expertise.
The large number of storage backends makes the `rdataretriever` easy to integrate into existing data processing workflows and to implement designs appropriate to the scale of the data with no additional work. 

Having installed the data into SQLite we can then connect to the database to start analyzing the data.

```{r, eval=FALSE}
bbs_db <- dbConnect(RSQLite::SQLite(), 'bbs.sqlite')
```

The two key tables for this analysis are the `surveys` and `sites` tables, so let's create connections to those tables.

```{r, eval=FALSE}
surveys <- tbl(bbs_db, "breed_bird_survey_counts")
sites <- tbl(bbs_db, "breed_bird_survey_routes")
```

The `surveys` table holds the data on how many individuals of each species are sampled at each site.
The `sites` table holds information on where each site is which we'll use to link the data to environmental variables.

## Analyze The Data

To calculate the measure of biodiversity, which is species richness or the number of species, we'll use `dplyr` to determine the number of species observed at each site in a recent year.

```{r, eval=FALSE}
rich_data <- surveys %>%
  filter(year == 2016) %>%
  group_by(statenum, route) %>%
  summarize(richness = n()) %>%
  collect()
rich_data
```

The data is now smaller than the original ~1 GB, so we used `collect` to load the summarized data directly into R.

Next we need to get environmental data for each site, which we'll get from the `worldclim` dataset.

```{r, eval=FALSE}
bioclim <- getData('worldclim', var = 'bio', res = 10)
```

To extract the environmental data we first make our sites data spatial and add them to our map.

```{r, eval=FALSE}
sites <- as.data.frame(sites)
sites_spatial <- SpatialPointsDataFrame(sites[c('longitude', 'latitude')], sites)
```

We can then extract the environmental data for each site from the bioclim raster and add it to the data on biodiversity.

```{r, eval=FALSE}
bioclim_bbs <- extract(bioclim, sites_spatial) %>%
  cbind(sites)
richness_w_env <- inner_join(rich_data, bioclim_bbs)
richness_w_env
```

Now let's see how richness relates to the precipitation.
Annual precipition is stored in `bio12`.

```{r, eval=FALSE}
ggplot(richness_w_env, aes(x = bio12, y = richness)) +
  geom_point(alpha = 0.5) +
  labs(x = "Annual Precipitation", y = "Number of Species")
```

It looks like there's a pattern here, so let's fit a smoother through it.

```{r, eval=FALSE}
ggplot(richness_w_env, aes(x = bio12, y = richness)) +
  geom_point(alpha = 0.5) +
  geom_smooth()
```

This shows that there is low bird biodiversity in really dry areas, biodiversity peaks at intermediate precipitations, and then drops off at the highest precipitation values.

If we wanted to use this kind of information to inform conservation decisions at the state level, could look at the patterns within each state after filtering to ensure enough data points.

```{r, fig.height = 8, fig.width = 8, warning = FALSE, eval=FALSE}
richness_w_env_high_n <- richness_w_env %>%
  group_by(statenum) %>%
  filter(n() >= 50)

ggplot(richness_w_env_high_n, aes(x = bio12, y = richness)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~statenum, scales = 'free') +
  labs(x = "Annual Precipitation", y = "Number of Species")
```

Looking back at this demo there is only one line of code directly involving the `rdataretriever` (other than installation).
This demonstrates the strength that the `rdataretriever` brings to the the early phases of the data acquistion and processing pipeline by distilling those steps to a single line that provides the data in a ready-to-analyze form so that researchers can focus on the analysis of the data itself.

## Conclusion

Thanks to the `rdataretriever` we can generate meaningful information about large scale bird biodiversity patterns in about 15 minutes.
If we'd been working with the raw BBS data we would likely have spent hours manipulating and cleaning data before we could start this analysis.

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{fetch}
\alias{fetch}
\title{Fetch a dataset via the Data Retriever}
\usage{
fetch(dataset, quiet = TRUE, data_names = NULL)
}
\arguments{
\item{dataset}{the names of the dataset that you wish to download}

\item{quiet}{logical, if true retriever runs in quiet mode}

\item{data_names}{the names you wish to assign to cells of the list which
stores the fetched dataframes. This is only relevant if you are
downloading more than one dataset.}
}
\description{
Each datafile in a given dataset is downloaded to a temporary directory and
then imported as a data.frame as a member of a named list.
}
\examples{
\dontrun{
## fetch the portal Database
portal <- rdataretriever::fetch("portal")
class(portal)
names(portal)
## preview the data in the portal species datafile
head(portal$species)
vegdata <- rdataretriever::fetch(c("plant-comp-ok", "plant-occur-oosting"))
names(vegdata)
names(vegdata$plant_comp_ok)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{get_script_citation}
\alias{get_script_citation}
\title{Get citation}
\usage{
get_script_citation(dataset = NULL)
}
\arguments{
\item{dataset}{dataset to obtain citation}
}
\description{
Get citation
}
\examples{
\dontrun{
rdataretriever::get_script_citation(dataset = "")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{use_RetrieverPath}
\alias{use_RetrieverPath}
\title{Setting path of retriever}
\usage{
use_RetrieverPath(path)
}
\arguments{
\item{path}{location of retriever in the system}
}
\description{
Setting path of retriever
}
\examples{
\dontrun{
rdataretriever::use_RetrieverPath("/home/<system_name>/anaconda2/envs/py27/bin/")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{reload_scripts}
\alias{reload_scripts}
\title{Update the retriever's global_script_list with the scripts present
in the ~/.retriever directory.}
\usage{
reload_scripts()
}
\description{
Update the retriever's global_script_list with the scripts present
in the ~/.retriever directory.
}
\examples{
\dontrun{
rdataretriever::reload_scripts()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_msaccess}
\alias{install_msaccess}
\title{Install datasets via the Data Retriever.}
\usage{
install_msaccess(
  dataset,
  file = "access.mdb",
  table_name = "[{db} {table}]",
  debug = FALSE,
  use_cache = TRUE,
  force = FALSE,
  hash_value = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to install or path to a committed dataset zip file}

\item{file}{file name for database}

\item{table_name}{table name for installing of dataset}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{Setting FALSE reinstalls scripts even if they are already installed}

\item{force}{setting TRUE doesn't prompt for confirmation while installing committed datasets when changes are discovered in environment}

\item{hash_value}{the hash value of committed dataset when installing from provenance directory}
}
\description{
Data is stored in MSAccess database
}
\examples{
\dontrun{
rdataretriever::install_msaccess(dataset = "iris", file = "sqlite.db")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{reset}
\alias{reset}
\title{Reset the scripts or data(raw_data) directory or both}
\usage{
reset(scope = "all")
}
\arguments{
\item{scope}{All resets both scripst and data directory}
}
\description{
Reset the scripts or data(raw_data) directory or both
}
\examples{
\dontrun{
rdataretriever::reset("iris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install}
\alias{install}
\title{Install datasets via the Data Retriever (deprecated).}
\usage{
install(
  dataset,
  connection,
  db_file = NULL,
  conn_file = NULL,
  data_dir = ".",
  log_dir = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to download}

\item{connection}{what type of database connection should be used.
The options include: mysql, postgres, sqlite, msaccess, or csv'}

\item{db_file}{the name of the datbase file the dataset should be loaded
into}

\item{conn_file}{the path to the .conn file that contains the connection
configuration options for mysql and postgres databases. This defaults to
mysql.conn or postgres.conn respectively. The connection file is a file that
is formated in the following way:
\tabular{ll}{
  host     \tab my_server@my_host.com\cr
  port     \tab my_port_number       \cr
  user     \tab my_user_name         \cr
  password \tab my_password
}}

\item{data_dir}{the location where the dataset should be installed.
Only relevant for csv connection types. Defaults to current working directory}

\item{log_dir}{the location where the retriever log should be stored if
the progress is not printed to the console}
}
\description{
Data is stored in either CSV files or one of the following database management
systems: MySQL, PostgreSQL, SQLite, or Microsoft Access.
}
\examples{
\dontrun{
rdataretriever::install("iris", "csv")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{socrata_dataset_info}
\alias{socrata_dataset_info}
\title{Get socrata dataset info}
\usage{
socrata_dataset_info(dataset_name)
}
\arguments{
\item{dataset_name}{dataset name to obtain info}
}
\description{
Get socrata dataset info
}
\examples{
\dontrun{
rdataretriever::socrata_dataset_info()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{dataset_names}
\alias{dataset_names}
\title{Name all available dataset scripts.}
\usage{
dataset_names()
}
\value{
returns a character vector with the available datasets for download
}
\description{
Additional information on the available datasets can be found at url https://retriever.readthedocs.io/en/latest/datasets.html
}
\examples{
\dontrun{
rdataretriever::dataset_names()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_xml}
\alias{install_xml}
\title{Install datasets via the Data Retriever.}
\usage{
install_xml(
  dataset,
  table_name = "{db}_{table}.xml",
  data_dir = getwd(),
  debug = FALSE,
  use_cache = TRUE,
  force = FALSE,
  hash_value = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to install or path to a committed dataset zip file}

\item{table_name}{the name of the database file to store data}

\item{data_dir}{the dir path to store data, defaults to working dir}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{Setting FALSE reinstalls scripts even if they are already installed}

\item{force}{setting TRUE doesn't prompt for confirmation while installing committed datasets when changes are discovered in environment}

\item{hash_value}{the hash value of committed dataset when installing from provenance directory}
}
\description{
Data is stored in XML files
}
\examples{
\dontrun{
rdataretriever::install_xml("iris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{update_rdataset_catalog}
\alias{update_rdataset_catalog}
\title{Updates the datasets_url.json from the github repo}
\usage{
update_rdataset_catalog(test = FALSE)
}
\arguments{
\item{test}{flag set when testing}
}
\description{
Updates the datasets_url.json from the github repo
}
\examples{
\dontrun{
rdataretriever::update_rdataset_catalog()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{get_retriever_citation}
\alias{get_retriever_citation}
\title{Get retriever citation}
\usage{
get_retriever_citation()
}
\description{
Get retriever citation
}
\examples{
\dontrun{
rdataretriever::get_retriever_citation()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{get_dataset_names_upstream}
\alias{get_dataset_names_upstream}
\title{Get dataset names from upstream}
\usage{
get_dataset_names_upstream(keywords = "", licenses = "", repo = "")
}
\arguments{
\item{keywords}{filter datasets based on keywords}

\item{licenses}{filter datasets based on license}

\item{repo}{path to the repository}
}
\description{
Get dataset names from upstream
}
\examples{
\dontrun{
rdataretriever::get_dataset_names_upstream(keywords = "", licenses = "", repo = "")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{get_rdataset_names}
\alias{get_rdataset_names}
\title{Returns a list of all the available RDataset names present}
\usage{
get_rdataset_names()
}
\description{
Returns a list of all the available RDataset names present
}
\examples{
\dontrun{
rdataretriever::get_rdataset_names()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{download}
\alias{download}
\title{Download datasets via the Data Retriever.}
\usage{
download(
  dataset,
  path = "./",
  quiet = FALSE,
  sub_dir = "",
  debug = FALSE,
  use_cache = TRUE
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to download}

\item{path}{the path where the data should be downloaded to}

\item{quiet}{logical, if true retriever runs in quiet mode}

\item{sub_dir}{downloaded dataset is stored into a custom subdirectory.}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{Setting FALSE reinstalls scripts even if they are already installed}
}
\description{
Directly downloads data files with no processing, allowing downloading of
non-tabular data.
}
\examples{
\dontrun{
rdataretriever::download("plant-comp-ok")
# downloaded files will be copied to your working directory
# when no path is specified
dir()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_retriever}
\alias{install_retriever}
\title{install the python module `retriever`}
\usage{
install_retriever(method = "auto", conda = "auto")
}
\arguments{
\item{method}{Installation method. By default, "auto" automatically finds a
method that will work in the local environment. Change the default to force
a specific installation method. Note that the "virtualenv" method is not
available on Windows.}

\item{conda}{The path to a \code{conda} executable. Use \code{"auto"} to allow
\code{reticulate} to automatically find an appropriate \code{conda} binary. See
\strong{Finding Conda} for more details.}
}
\description{
install the python module `retriever`
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{display_all_rdataset_names}
\alias{display_all_rdataset_names}
\title{Displays the list of rdataset names present in the list of packages provided}
\usage{
display_all_rdataset_names(package_name = NULL)
}
\arguments{
\item{package_name}{print datasets in the package, default to print rdataset and all to print all}
}
\description{
Can take a list of packages, or NULL or a string 'all' for all rdataset packages and datasets
}
\examples{
\dontrun{
rdataretriever::display_all_rdataset_names()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{socrata_autocomplete_search}
\alias{socrata_autocomplete_search}
\title{Returns the list of dataset names after autocompletion}
\usage{
socrata_autocomplete_search(dataset)
}
\arguments{
\item{dataset}{the name of the dataset}
}
\description{
Returns the list of dataset names after autocompletion
}
\examples{
\dontrun{
rdataretriever::socrata_autocomplete_search()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_csv}
\alias{install_csv}
\title{Install datasets via the Data Retriever.}
\usage{
install_csv(
  dataset,
  table_name = "{db}_{table}.csv",
  data_dir = getwd(),
  debug = FALSE,
  use_cache = TRUE,
  force = FALSE,
  hash_value = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to install or path to a committed dataset zip file}

\item{table_name}{the name of the database file to store data}

\item{data_dir}{the dir path to store data, defaults to working dir}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{Setting FALSE reinstalls scripts even if they are already installed}

\item{force}{setting TRUE doesn't prompt for confirmation while installing committed datasets when changes are discovered in environment}

\item{hash_value}{the hash value of committed dataset when installing from provenance directory}
}
\description{
Data is stored in CSV files
}
\examples{
\dontrun{
rdataretriever::install_csv("iris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{commit_log}
\alias{commit_log}
\title{See the log of committed dataset stored in provenance directory}
\usage{
commit_log(dataset)
}
\arguments{
\item{dataset}{name of the dataset stored in provenance directory}
}
\description{
See the log of committed dataset stored in provenance directory
}
\examples{
\dontrun{
rdataretriever::commit_log("iris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{find_socrata_dataset_by_id}
\alias{find_socrata_dataset_by_id}
\title{Returns metadata for the following dataset id}
\usage{
find_socrata_dataset_by_id(dataset_id)
}
\arguments{
\item{dataset_id}{id of the dataset}
}
\description{
Returns metadata for the following dataset id
}
\examples{
\dontrun{
rdataretriever::socrata_dataset_info()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_postgres}
\alias{install_postgres}
\title{Install datasets via the Data Retriever.}
\usage{
install_postgres(
  dataset,
  user = "postgres",
  password = "",
  host = "localhost",
  port = 5432,
  database = "postgres",
  database_name = "{db}",
  table_name = "{db}.{table}",
  bbox = list(),
  debug = FALSE,
  use_cache = TRUE,
  force = FALSE,
  hash_value = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to install or path to a committed dataset zip file}

\item{user}{username for database connection}

\item{password}{password for database connection}

\item{host}{hostname for connection}

\item{port}{port number for connection}

\item{database}{the database name default is postres}

\item{database_name}{database schema name in which dataset will be installed}

\item{table_name}{table name specified especially for datasets
containing one file}

\item{bbox}{optional extent values used to fetch data from the spatial dataset}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{setting FALSE reinstalls scripts even if they are already installed}

\item{force}{setting TRUE doesn't prompt for confirmation while installing committed datasets when changes are discovered in environment}

\item{hash_value}{the hash value of committed dataset when installing from provenance directory}
}
\description{
Data is stored in PostgreSQL database
}
\examples{
\dontrun{
rdataretriever::install_postgres(dataset = "portal", user = "postgres", password = "abcdef")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{datasets}
\alias{datasets}
\title{Name all available dataset scripts.}
\usage{
datasets(keywords = "", licenses = "")
}
\arguments{
\item{keywords}{search all datasets by keywords}

\item{licenses}{search all datasets by licenses}
}
\value{
returns a character vector with the available datasets for download
}
\description{
Additional information on the available datasets can be found at url https://retriever.readthedocs.io/en/latest/datasets.html
}
\examples{
\dontrun{
rdataretriever::datasets()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{data_retriever_version}
\alias{data_retriever_version}
\title{Get Data Retriever version}
\usage{
data_retriever_version(clean = TRUE)
}
\arguments{
\item{clean}{boolean return cleaned version appropriate for semver}
}
\value{
returns a string with the version information
}
\description{
Get Data Retriever version
}
\examples{
\dontrun{
rdataretriever::data_retriever_version()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{commit}
\alias{commit}
\title{Commit a dataset}
\usage{
commit(dataset, commit_message = "", path = NULL, quiet = FALSE)
}
\arguments{
\item{dataset}{name of the dataset}

\item{commit_message}{commit message for the commit}

\item{path}{path to save the committed dataset, if no path given save in provenance directory}

\item{quiet}{logical, if true retriever runs in quiet mode}
}
\description{
Commit a dataset
}
\examples{
\dontrun{
rdataretriever::commit("iris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{get_updates}
\alias{get_updates}
\title{Update the retriever's dataset scripts to the most recent versions.}
\usage{
get_updates()
}
\description{
This function will check if the version of the retriever's scripts in your local
directory \file{~/.retriever/scripts/} is up-to-date with the most recent official
retriever release. Note it is possible that even more updated scripts exist
at the retriever repository \url{https://github.com/weecology/retriever/tree/main/scripts}
that have not yet been incorperated into an official release, and you should
consider checking that page if you have any concerns.
}
\examples{
\dontrun{
rdataretriever::get_updates()
}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_mysql}
\alias{install_mysql}
\title{Install datasets via the Data Retriever.}
\usage{
install_mysql(
  dataset,
  user = "root",
  password = "",
  host = "localhost",
  port = 3306,
  database_name = "{db}",
  table_name = "{db}.{table}",
  debug = FALSE,
  use_cache = TRUE,
  force = FALSE,
  hash_value = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to install or path to a committed dataset zip file}

\item{user}{username for database connection}

\item{password}{password for database connection}

\item{host}{hostname for connection}

\item{port}{port number for connection}

\item{database_name}{database name in which dataset will be installed}

\item{table_name}{table name specified especially for datasets
containing one file}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{setting FALSE reinstalls scripts even if they are already installed}

\item{force}{setting TRUE doesn't prompt for confirmation while installing committed datasets when changes are discovered in environment}

\item{hash_value}{the hash value of committed dataset when installing from provenance directory}
}
\description{
Data is stored in MySQL database
}
\examples{
\dontrun{
rdataretriever::install_mysql(dataset = "portal", user = "postgres", password = "abcdef")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{check_for_updates}
\alias{check_for_updates}
\title{Check for updates}
\usage{
check_for_updates(repo = "")
}
\arguments{
\item{repo}{path to the repository}
}
\description{
Check for updates
}
\examples{
\dontrun{
rdataretriever::check_for_updates()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_json}
\alias{install_json}
\title{Install datasets via the Data Retriever.}
\usage{
install_json(
  dataset,
  table_name = "{db}_{table}.json",
  data_dir = getwd(),
  debug = FALSE,
  use_cache = TRUE,
  force = FALSE,
  hash_value = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to install or path to a committed dataset zip file}

\item{table_name}{the name of the database file to store data}

\item{data_dir}{the dir path to store data, defaults to working dir}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{setting FALSE reinstalls scripts even if they are already installed}

\item{force}{setting TRUE doesn't prompt for confirmation while installing committed datasets when changes are discovered in environment}

\item{hash_value}{the hash value of committed dataset when installing from provenance directory}
}
\description{
Data is stored in JSON files
}
\examples{
\dontrun{
rdataretriever::install_json("iris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{install_sqlite}
\alias{install_sqlite}
\title{Install datasets via the Data Retriever.}
\usage{
install_sqlite(
  dataset,
  file = "sqlite.db",
  table_name = "{db}_{table}",
  data_dir = getwd(),
  debug = FALSE,
  use_cache = TRUE,
  force = FALSE,
  hash_value = NULL
)
}
\arguments{
\item{dataset}{the name of the dataset that you wish to install or path to a committed dataset zip file}

\item{file}{Sqlite database file name or path}

\item{table_name}{table name for installing of dataset}

\item{data_dir}{the dir path to store the db, defaults to working dir}

\item{debug}{setting TRUE helps in debugging in case of errors}

\item{use_cache}{setting FALSE reinstalls scripts even if they are already installed}

\item{force}{setting TRUE doesn't prompt for confirmation while installing committed datasets when changes are discovered in environment}

\item{hash_value}{the hash value of committed dataset when installing from provenance directory}
}
\description{
Data is stored in SQLite database
}
\examples{
\dontrun{
rdataretriever::install_sqlite(dataset = "iris", file = "sqlite.db")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{check_retriever_availability}
\alias{check_retriever_availability}
\title{Check to see if minimum version of retriever Python package is installed}
\usage{
check_retriever_availability()
}
\value{
boolean
}
\description{
Check to see if minimum version of retriever Python package is installed
}
\examples{
\dontrun{
rdataretriever::check_retriever_availability()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdataretriever.R
\name{get_script_upstream}
\alias{get_script_upstream}
\title{Get script upstream}
\usage{
get_script_upstream(dataset, repo = "")
}
\arguments{
\item{dataset}{name of the dataset}

\item{repo}{path to the repository}
}
\description{
Get script upstream
}
\examples{
\dontrun{
rdataretriever::get_script_upstream("iris")
}
}
