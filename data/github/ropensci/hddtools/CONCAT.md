# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# hddtools: Hydrological Data Discovery Tools

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/22423032.svg)](https://zenodo.org/badge/latestdoi/22423032)
[![status](https://joss.theoj.org/papers/10.21105/joss.00056/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00056)

[![CRAN Status
Badge](http://www.r-pkg.org/badges/version/hddtools)](https://cran.r-project.org/package=hddtools)
[![CRAN Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/hddtools)](https://cran.r-project.org/package=hddtools)
[![CRAN Monthly
Downloads](http://cranlogs.r-pkg.org/badges/hddtools)](https://cran.r-project.org/package=hddtools)

[![R-CMD-check](https://github.com/ropensci/hddtools/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/hddtools/actions)
[![codecov.io](https://codecov.io/github/ropensci/hddtools/coverage.svg?branch=master)](https://codecov.io/github/ropensci/hddtools?branch=master)
[![](https://badges.ropensci.org/73_status.svg)](https://github.com/ropensci/software-review/issues/73)
<!-- badges: end -->

`hddtools` stands for Hydrological Data Discovery Tools. This R package
is an open source project designed to facilitate access to a variety of
online open data sources relevant for hydrologists and, in general,
environmental scientists and practitioners.

This typically implies the download of a metadata catalogue, selection
of information needed, a formal request for dataset(s), de-compression,
conversion, manual filtering and parsing. All those operations are made
more efficient by re-usable functions.

Depending on the data license, functions can provide offline and/or
online modes. When redistribution is allowed, for instance, a copy of
the dataset is cached within the package and updated twice a year. This
is the fastest option and also allows offline use of package’s
functions. When re-distribution is not allowed, only online mode is
provided.

## Installation

Get the stable version from CRAN:

``` r
install.packages("hddtools")
```

Or the development version from GitHub using the package `remotes`:

``` r
install.packages("remotes")
remotes::install_github("ropensci/hddtools")
```

Load the `hddtools` package:

``` r
library("hddtools")
```

## Data sources and Functions

The package contains functions to interact with the data providers
listed below. For examples of the various functionalities see the
[vignette](https://github.com/ropensci/hddtools/blob/master/vignettes/hddtools_vignette.Rmd).

  - [KGClimateClass](http://koeppen-geiger.vu-wien.ac.at/): The Koppen
    Climate Classification map is used for classifying the world’s
    climates based on the annual and monthly averages of temperature and
    precipitation.

  - [GRDC](http://www.bafg.de/GRDC/): The Global Runoff Data Centre
    (GRDC) provides datasets for all the major rivers in the world.

  - [Data60UK](http://tdwg.catchment.org/datasets.html): The Data60UK
    initiative collated datasets of areal precipitation and streamflow
    discharge across 61 gauging sites in England and Wales (UK).

  - [MOPEX](http://tdwg.catchment.org/datasets.html): This dataset
    contains historical hydrometeorological data and river basin
    characteristics for hundreds of river basins in the US.

  - [SEPA](https://www2.sepa.org.uk/WaterLevels/): The Scottish
    Environment Protection Agency (SEPA) provides river level data for
    hundreds of gauging stations in the UK.

## Meta

  - This package and functions herein are part of an experimental
    open-source project. They are provided as is, without any guarantee.
  - Please note that this project is released with a [Contributor Code
    of Conduct](https://github.com/ropensci/hddtools/blob/master/CONDUCT.md).
    By participating in this project you agree to abide by its terms.
  - Please [report any issues or bugs](https://github.com/ropensci/hddtools/issues).
  - License: [GPL-3](https://opensource.org/licenses/GPL-3.0)
  - This package was reviewed by [Erin Le
    Dell](https://github.com/ledell) and [Michael
    Sumner](https://github.com/mdsumner) for submission to ROpenSci (see
    review [here](https://github.com/ropensci/software-review/issues/73)) and
    the Journal of Open Source Software (see review status
    [here](https://github.com/openjournals/joss-reviews/issues/56)).
  - Cite `hddtools`: `citation(package = "rdefra")`

<br/>

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# hddtools 0.9.4

* Removed dependency from the rnrfa package

# hddtools 0.9.3

* Updated URLs
* SEPA tests now skipped if source is temporarily unavailable

# hddtools 0.9.0

* Updated URLs of SEPA and MOPEX data sources
* Removed obsolete dependencies
* Amended URLs
* Simplified API
* Removed evaluation of GRDC catalogue in vignette
* Using getBinaryURL() to retrieve zip file from ftp

# hddtools 0.8.3

* Updated URLs of SEPA and MOPEX data sources
* Removed obsolete dependencies
* Amended URLs
* Simplified API
* Switched from download.file to downloader::download() to better work cross platform

# hddtools 0.8.2

* Connection to TRMM database removed as the web service is no longer available.
* Fixed problem with Rd file which prevented the manual to be created.
* Fixed problem with URL in DESCRIPTION file

# hddtools 0.8.1

* Any reference to TRMM database was removed as the web service is no longer available.

# hddtools 0.8

* tsGRDC function: this function now returns 6 tables (see documentation) according to latest updates on the GRDC side.
* GRDC data catalogue was updated in October 2017.

# hddtools 0.7

* TRMM function: this function has been temporarily removed from hddtools v0.7 as the ftp at NASA containing the data has been migrated. A new function is under development.

# hddtools 0.6

* TRMM function: set method for download.file to "auto" and added to arguments
* TRMM function: downloaded files are in temporary folder
* TRMM function: removed inputfileLocation
* TRMM function: added support for 3B42
* HadDAILY function removed as the service is no longer available
* The output of catalogueGRDC() is now a tibble table 

# hddtools 0.5

* Updated all the links to show the new repository (rOpenSci)

# hddtools 0.4

* Added a vignette
* Added paper for submission to JOSS
* Fixed bug with the TRMM function (for Windows users)

# hddtools 0.3

* Added unist tests (using testthat framework), 
* Added Travis CI
* Added Appveyor

# hddtools 0.2

* Added functions to get data from the following source:

  - TRMM
  - DATA60UK
  - MOPEX
  - GRDC
  - SEPA
  - Met Office Hadley Centre

# hddtools 0.1

* Initial release
This is a resubmission after the package was archived.
New release (hddtools v0.9.4).

---------------------------------

## Release Summary

* Removed dependency from the rnrfa package

## Test environment
* Ubuntu 18.04, R 3.6.3

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.
---
title: 'hddtools: Hydrological Data Discovery Tools'
authors:
- affiliation: 1
  name: Claudia Vitolo
  orcid: 0000-0002-4252-1176
date: "27 December 2016"
output: pdf_document
bibliography: paper.bib
csl: hydrology-and-earth-system-sciences.csl
tags:
- open data
- hydrology
- R
affiliations:
- index: 1
  name: European Centre for Medium-range Weather Forecasts
---

# Summary

The hddtools [@hddtoolsCRAN] (**h**ydrological **d**ata **d**iscovery **tools**) is an R package [@R-base] designed to facilitate access to a variety of online open data sources relevant for hydrologists and, in general, environmental scientists and practitioners. This typically implies the download of a metadata catalogue, selection of information needed, formal request for dataset(s), de-compression, conversion, manual filtering and parsing. All those operation are made more efficient by re-usable functions. 

Depending on the data license, functions can provide offline and/or online modes. When redistribution is allowed, for instance, a copy of the dataset is cached within the package and updated twice a year. This is the fastest option and also allows offline use of package's functions. When re-distribution is not allowed, only online mode is provided.

Datasets for which functions are provided include: the Global Runoff Data Center (GRDC), the Scottish Environment Protection Agency (SEPA), the Top-Down modelling Working Group (Data60UK and MOPEX), Met Office Hadley Centre Observation Data (HadUKP Data) and NASA's Tropical Rainfall Measuring Mission (TRMM). 

This package follows a logic similar to other packages such as rdefra [@rdefraJOSS] and rnrfa [@rnrfa]: sites are first identified through a catalogue (if available), data are imported via the station identification number, then data are visualised and/or used in analyses. The metadata related to the monitoring stations are accessible through the functions: `catalogueGRDC()`, `catalogueSEPA()`, `catalogueData60UK()` and  `catalogueMOPEX()`. Time series data can be obtained using the functions: `tsGRDC()`, `tsSEPA()`, `tsData60UK()`, `tsMOPEX()` and `HadDAILY()`. Geospatial information can be retrieved using the functions: `KGClimateClass()` returning the Koppen-Greiger climate zone and `TRMM()` which retrieves global historical rainfall estimations.

The retrieved hydrological time series (e.g. using `tsData60UK()`) can be used to feed hydrological models such as fuse [@fuseGitHub; @fuseJOSS], topmodel [@topmodel] and hydromad [@Andrews20111171; @hydromad].

For more details and examples, please refer to the help pages and vignette.

# References
