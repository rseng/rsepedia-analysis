---
title: 'OpenTripPlanner for R'
authors:
- affiliation: 1
  name: Malcolm Morgan
  orcid: 0000-0002-9488-9183
- affiliation: 2
  name: Marcus Young
  orcid: 0000-0003-4627-1116
- affiliation: 1
  name: Robin Lovelace
  orcid: 0000-0001-5679-6536
- affiliation: 1
  name: Layik Hama
  orcid: 0000-0003-1912-4890
date: "27 Nov 2019"
output:
  html_document: default
  pdf_document: default
tags:
- R
- transport
- geospatial
- transit
- OpenTripPlanner
- OpenStreetMap
affiliations:
- index: 1
  name: Institute for Transport Studies, University of Leeds, UK
- index: 2
  name: Transportation Research Group, University of Southampton, UK
bibliography: paper.bib
---

<!--
generate citations (in R)
refs = RefManageR::ReadZotero(group = "418217", .params = list(collection = "JFR868KJ", limit = 100))
RefManageR::WriteBib(refs, "paper.bib")
citr::tidy_bib_file(rmd_file = "paper.md", messy_bibliography = "paper.bib")
-->

# Summary

**opentripplanner** provides functions that enable multi-modal routing in R.
It provides an interface to the OpenTripPlanner (OTP) Java routing service, which allows calculation of routes on large transport networks, locally or via calls to a remote server.
The package contains three groups of functions for: (1) setting up and managing a local instance of OTP; (2) connecting to OTP locally or remotely; and (3) sending requests to the OTP API and importing the results and converting them into appropriate classes for further analysis.

# Motivation

Routing, the process of calculating paths between points in geographic space that follow a transport network, is a fundamental part of transport planning and vital to solving many real-world transport problems.
The outputs of a routing service, which can calculate many routes over a large area, comprises coordinates and other data representing a movement that is in some way 'optimal', based on a range of criteria (that the user should understand and be able to change).
Routing services, such as the one provided by Google Maps, are a well-known and increasingly vital component of personal travel planning for many people [@bast_fast_2010].
Less well-known, but perhaps equally important, is that routing services are also key to understanding aggregate travel patterns and guiding policy and commercial decisions [@giusti_new_2018].
To meet this need for route planning capabilities, a wide range of both proprietary and open source tools have been created.

# Functionality

[OpenTripPlanner](https://www.opentripplanner.org/) (OTP) is written in Java and designed to work with [Open Street Map](https://www.openstreetmap.org) (OSM) data for road-based modes (Car, Bike, Walking) and [General Transit Feed Specification]( https://developers.google.com/transit/gtfs/) (GTFS) data for public transit (Bus, Train, Metro).
OTP is unusual among open source routing tools in its ability to account for a wide range of modes of travel, e.g. bicycle hire, and support for complex multi-stage multi-modal journeys such as park and ride. 
However, OTP’s primary purpose is to support public-facing websites such as TriMet; thus its analytical capabilities are limited.
Conversely, the R language is well suited to statistical and spatial analysis but has no route planning capabilities.

The OpenTripPlanner for R package aims to bridge the gap between OTP and R by supplying simple ways for R to connect to OTP either on a local machine or on a remote server, via OTP’s API.
The package has been designed to ease bulk routing by allowing the input of multiple origins and destinations as two column matrices of longitude-latitude pairs.
The package also supports multi-core operation to take advantage of OTP’s multicore functionality.
Results are returned in the widely used [`sf` data frame](https://CRAN.R-project.org/package=sf) format.
Although performance is dependant on the size of the map being routed over, it typically can achieve more than 10 routes per second per core.

The package has been developed from a set of R functions that formed part of an intermediate-level [OTP tutorial](https://github.com/marcusyoung/otp-tutorial/raw/master/intro-otp.pdf) as part of research at [Centre for Research into Energy Demand Solutions]( https://www.creds.ac.uk/) and the [Institute of Transport Studies](https://environment.leeds.ac.uk/transport).

# Reproducible demonstration

Example data for the Isle of Wight, UK is provided with the package. The example below uses this data to demonstrate the basic functionality of the package. A full explanation is provided in the [package vignettes](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html)

First, download the demo data and OTP.

``` r
library(opentripplanner)                 # Load Package
path_data <- file.path(tempdir(), "OTP") # Create folder for data
dir.create(path_data)
path_otp <- otp_dl_jar(path_data)        # Download OTP
otp_dl_demo(path_data)                   # Download demo data
```

Second, build the OTP graph, start up OTP server and connect to the server

``` r
log <- otp_build_graph(otp = path_otp,     # Build Graph
                       dir = path_data)
otp_setup(otp = path_otp, dir = path_data) # Start OTP
otpcon <- otp_connect()                    # Connect R to OTP
```

Finally, find routes

``` r
route <- otp_plan(otpcon, 
                 fromPlace = c(-1.17502, 50.64590), 
                 toPlace = c(-1.15339, 50.72266))
```

# References


<!-- README.md is generated from README.Rmd. Please edit that file -->

# OpenTripPlanner for R <a href='https://itsleeds.github.io/'><img src='man/figures/logo.png' align="right" height=180/></a>

[![R build
status](https://github.com/ropensci/opentripplanner/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/opentripplanner/actions)
[![codecov](https://codecov.io/gh/ropensci/opentripplanner/branch/master/graph/badge.svg?token=iLEB77PnMk)](https://app.codecov.io/gh/ropensci/opentripplanner)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://badges.ropensci.org/295_status.svg)](https://github.com/ropensci/software-review/issues/295)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3558311.svg)](https://doi.org/10.5281/zenodo.3558311)
[![status](https://joss.theoj.org/papers/10.21105/joss.01926/status.svg)](https://joss.theoj.org/papers/10.21105/joss.01926)
[![](https://cranlogs.r-pkg.org/badges/grand-total/opentripplanner)](https://cran.r-project.org/package=opentripplanner)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/opentripplanner)](https://cran.r-project.org/package=opentripplanner)

**opentripplanner** is an R package that provides a simple yet flexible
interface to [OpenTripPlanner (OTP)](https://www.opentripplanner.org/).
OTP is a multimodal trip planning service written in Java. For more
information on what OTP is, see the [prerequisites
vignette](https://docs.ropensci.org/opentripplanner/articles/prerequisites.html).

**opentripplanner** can be used to interface with a remote instance of
OTP (e.g. a website) or help you set up and manage a local version of
OTP for private use. Basic setup and routing functions are outlined in
the [getting started
vignette](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html),
while advanced functionality such as batch routing, isochrones, and
customised setup is described in the [advanced features
vignette](https://docs.ropensci.org/opentripplanner/articles/advanced_features.html).

## Installation

### OpenTripPlanner

To use OpenTripPlanner on your local computer you will need to install
Java 8 and download the latest version of OTP. Instructions on
installing Java and setting up OTP can be found in the [prerequisites
vignette](https://docs.ropensci.org/opentripplanner/articles/prerequisites.html).

### R Package

To install the stable CRAN version:

``` r
install.packages("opentripplanner") # Install Package
library(opentripplanner)            # Load Package
```

Install the development version using **remotes**:

``` r
# If you do not already have the remotes package
install.packages("remotes")
# Install the package from GitHub
remotes::install_github("ropensci/opentripplanner")
# Load the package
library(opentripplanner)
```

#### RcppSimdJson

From version 0.3.0 of `opentripplanner` the package `RcppSimdJson` is
used for JSON parsing. This package is not supported on some older
versions of R (\<= 3.6) and some older Operating Systems. To meet CRAN
requirements version 0.3.1 added a legacy mode for older versions of R.
This legacy mode has reduced functionality and users on old systems may
get better results using version 0.2.3 of the package. You can install
older versions using **remotes**.

``` r
remotes::install_version("opentripplanner", "0.2.3")
```

## Usage

The package contains three groups of functions:

Functions for setting up a local instance of OTP:

1.  `otp_dl_jar()` To download the OTP Jar file;
2.  `otp_dl_demo()` To download the demo data for the Isle of Wight;
3.  `otp_check_java()` To check you have the correct version of Java;
4.  `otp_build_graph()` To make a OTP graph from raw data;
5.  `otp_setup()` To start up a local instance of OTP;
6.  `otp_make_config()` To make a config object;
7.  `otp_validate_config()` To validate a config object;
8.  `otp_write_config()` To save a config object as a json file.

Functions for connecting to a local or remote instance of OTP:

1.  `otp_connect()` To connect to OTP.

Functions for retrieving data from OTP:

1.  `otp_plan()` To get routes from A to B;
2.  `otp_geocode()` To get the locations of named places e.g. road
    names.
3.  `otp_isochrone()` To get isochrone maps (OTP 1.x only);
4.  `otp_make_surface()` To make an analyst surface (OTP 1.x only);
5.  `otp_surface()` To evaluate a analyst surface (OTP 1.x only);
6.  `otp_traveltime()` To make a travel time matrix (OTP 1.x only);
7.  `otp_surface_isochrone()` To make a raster isochrone map (OTP 1.x
    only);

Results are returned as [sf
objects](https://CRAN.R-project.org/package=sf).

## Acknowledgement

This package was built off the [tutorial by Marcus
Young](https://github.com/marcusyoung/otp-tutorial).

## Contribution

Please note that the `opentripplanner` project is released with a
[Contributor Code of
Conduct](https://github.com/ropensci/opentripplanner/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms. Bug
reports and comments are welcome as Github
[Issues](https://github.com/ropensci/opentripplanner/issues) and code
submissions as [Pull
Requests](https://github.com/ropensci/opentripplanner/pulls).

## Citation

Please cite the JOSS paper in publications:

Morgan et al., (2019). OpenTripPlanner for R. Journal of Open Source
Software, 4(44), 1926, <https://doi.org/10.21105/joss.01926>

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# opentripplanner 0.4.0

* Fix broken or moved URLs
* `tibble` moved from imports to suggests
* Support for OTP 2.0
* Support for analyst features in OTP 1.x with new functions
* Support for `s2` features in the `sf` package

# opentripplanner 0.3.2

* Fix a bug with parsing fare data
* Fix bug with latest version of Java defaulting to 32Bit even when 64bit is available.
* Added more `try()` functions to reduce risk of crashes in large scale batch routing
* Added `flag64bit` argument to `otp_build_graph()` and `otp_setup()`
* Added `quiet` argument to `otp_build_graph()`
* Updated the Known Issues Vignette
* reduce initial wait time for `otp_setup()` from 60 seconds to 30 seconds
* switch to OTP 1.5.0 as default version in `otp_dl_jar()`

# opentripplanner 0.3.1

Limited support for version of R than can't install `RcppSimdJson`

# opentripplanner 0.3.0

Note that this version makes minor changes to how results are returned, for example column order. These changes are due to the new json parser and should not affect the overall results but may affect any dependent code.

* Significant refactor of code giving up to 50% reduction in routing time
* replaced `dplyr` with `data.table`
* replaced `httr` with `curl`
* replaced many `rjson` functions with `RcppSimdJson` equivalents
* replaced many `sf` functions with `sfheaders` equivalents
* fixed bug in `otp_plan` when `distance_balancing = TRUE`
* fixed bug #69
* In `otp_plan` will now return `fromPlace` and `toPlace` as the first two columns
* In `otp_plan` set `get_elevation = FALSE` as default this boosts performance
* Fixed bug in `distance_balancing` that gave sub-optimal balancing
* When `distance_balancing = TRUE` zero distance routes will not be found, as OTP will reject these in any case, this saves time with no impact on results.


# opentripplanner 0.2.3.0

* Added `distance_balancing` argument to `otp_plan` gives a small performance boost to multicore routing
* Added `get_elevation` argument to `otp_plan` default TRUE, when FALSE returns XY coordinates rather than XYZ coordinates and gives a 4% performance boost.
* Removed helper code for `dplyr::bind_rows` as no longer required for `dplyr 1.0.0`

# opentripplanner 0.2.2.0

* Changes to support `dplyr 1.0.0`, package now needs `vctrs 0.3.1`
* Added timezone support to `otp_connect`, `otp_plan`, and `otp_isochrone` fixing issue #54, see docs for details.
* Added `quiet` argument to `otp_dl_jar` and `otp_dl_demo`
* Fixed error in advanced features vignette, issue #57
* Switch from dplyr to data.table, issue #60

# opentripplanner 0.2.1.0

* Batch isochrones support added
* Fix bug in `correct_distances()` when input is of length <= 2 or the distances never decrease
* Fix bug in `polyline2linestring()` when elevation is length <= 2
* New input argument to `otp_plan()` and otp_isochrone routingOptions this allows support
    for many more routing options to be set. Arguments walkReluctance, transferPenalty, and
    minTransferTime have been removed and replaced with routingOptions.
* New functions `otp_routing_options()`, `otp_validate_routing_options()`, `otp_check_java()`
* New data for internal package checking
* Added support of the UK way property set (for OTP 1.5)

# opentripplanner 0.2.0.8

* Disabled CRAN tests that fail on solaris OS, due to different wording of error messages

# opentripplanner 0.2.0.6

* Fixed bug where routing fails due to missing fare data

# opentripplanner 0.2.0.4

* Added a `NEWS.md` file to track changes to the package.
* Fixed CRAN test error on certain OS
* Added CRAN badges 

# opentripplanner 0.2.0.3

* First CRAN release

# opentripplanner 0.2.0.0

* Release peer-reviewed by ROpenSci
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
---
name: Bug report
about: To report that the R package is not working correctly
title: "[BUG]"
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**System**
 Please paste `sessionInfo()` here:


**Additional context**
Add any other context about the problem here.
---
name: Need Help
about: For people who need help getting OTP working on their computer
title: "[Need Help] "
labels: help wanted
assignees: ''

---

Thank you for using OpenTripPlanner for R. This issue template is to help you provide the necessary information for us to help you.

## Pre Issue Check List

Before opening a new issue please confirm you have done the following:

- [ ] I have searched through [old issues](https://github.com/ropensci/opentripplanner/issues?q=is%3Aissue) for people with similar problems

- [ ] I have read the [Prerequisites Vignette](https://docs.ropensci.org/opentripplanner/articles/prerequisites.html) and installed the required software

- [ ] I have worked through the [Getting Started Vignette](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html) with the demo data.

- [ ] I have read the [Known Issues Vignette](https://docs.ropensci.org/opentripplanner/articles/known_issues.html) and checked if the solution to my problem already exists

- [ ] I have the latest version of the `opentripplanner` R package installed

- [ ] I have read any error messages including in the logs

## Problem Type

Tell us where you are having problems: Please tick all that apply:

- [ ]  Installing the R Package

- [ ] Installing Java

- [ ] Building the Graph: i.e. `otp_build_graph()`

- [ ] Starting OTP: i.e. `otp_setup()`

- [ ] Routing: i.e. `otp_plan()`

- [ ] Other

## Demo Data vs Own Data

Are you having a problem using the demo data or your own data?

- [ ]  Demo data: The problem might be with the `opentripplanner` package

- [ ] Own data: The problem might be with your data, so we may not be able to help you (you can still post and issue)

**If you are having a problem with your own data:**

1.  Can you reproduce the problem using the demo data?

1. Have you checked that your data is not missing / corrupted / faulty?

1. If you have multiple input datasets have you tried using a minimal subset of the data to identify which dataset is causing the problem?

## Session Info

Please run `sessionInfo()` and paste the results here:

## Your Issue

Please describe your problem here:
---
name: Feature request
about: Suggest an idea for this project
title: "[Feature Request]"
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# OpenTripPlanner for R <a href='https://itsleeds.github.io/'><img src='man/figures/logo.png' align="right" height=180/></a>

[![R build status](https://github.com/ropensci/opentripplanner/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/opentripplanner/actions)
[![codecov](https://codecov.io/gh/ropensci/opentripplanner/branch/master/graph/badge.svg?token=iLEB77PnMk)](https://app.codecov.io/gh/ropensci/opentripplanner) 
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://badges.ropensci.org/295_status.svg)](https://github.com/ropensci/software-review/issues/295)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3558311.svg)](https://doi.org/10.5281/zenodo.3558311)
[![status](https://joss.theoj.org/papers/10.21105/joss.01926/status.svg)](https://joss.theoj.org/papers/10.21105/joss.01926)
[![](https://cranlogs.r-pkg.org/badges/grand-total/opentripplanner)](https://cran.r-project.org/package=opentripplanner)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/opentripplanner)](https://cran.r-project.org/package=opentripplanner)


**opentripplanner** is an R package that provides a simple yet flexible interface to [OpenTripPlanner (OTP)](https://www.opentripplanner.org/). OTP is a multimodal trip planning service written in Java. For more information on what OTP is, see the [prerequisites vignette](https://docs.ropensci.org/opentripplanner/articles/prerequisites.html).

**opentripplanner** can be used to interface with a remote instance of OTP (e.g. a website) or help you set up and manage a local version of OTP for private use.  Basic setup and routing functions are outlined in the [getting started vignette](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html), while advanced functionality such as batch routing, isochrones, and customised setup is described in the [advanced features vignette](https://docs.ropensci.org/opentripplanner/articles/advanced_features.html).

## Installation

### OpenTripPlanner

To use OpenTripPlanner on your local computer you will need to install Java 8 and download the latest version of OTP. Instructions on installing Java and setting up OTP can be found in the [prerequisites vignette](https://docs.ropensci.org/opentripplanner/articles/prerequisites.html).

### R Package

To install the stable CRAN version:

```{r installCRAN, eval=FALSE}
install.packages("opentripplanner") # Install Package
library(opentripplanner)            # Load Package
```

Install the development version using **remotes**:

```{r install, eval=FALSE}
# If you do not already have the remotes package
install.packages("remotes")
# Install the package from GitHub
remotes::install_github("ropensci/opentripplanner")
# Load the package
library(opentripplanner)
```

#### RcppSimdJson

From version 0.3.0 of `opentripplanner` the package `RcppSimdJson` is used for JSON parsing. This package is not supported on some older versions of R (<= 3.6) and some older Operating Systems. To meet CRAN requirements version 0.3.1 added a legacy mode for older versions of R. This legacy mode has reduced functionality and users on old systems may get better results using version 0.2.3 of the package. You can install older versions using **remotes**.

```{r, eval=FALSE}
remotes::install_version("opentripplanner", "0.2.3")
```

## Usage

The package contains three groups of functions:

Functions for setting up a local instance of OTP:

1. `otp_dl_jar()`          To download the OTP Jar file;
1. `otp_dl_demo()`         To download the demo data for the Isle of Wight;
1. `otp_check_java()`      To check you have the correct version of Java;
1. `otp_build_graph()`     To make a OTP graph from raw data;
1. `otp_setup()`           To start up a local instance of OTP;
1. `otp_make_config()`     To make a config object;
1. `otp_validate_config()` To validate a config object;
1. `otp_write_config()`    To save a config object as a json file.

Functions for connecting to a local or remote instance of OTP:

1. `otp_connect()`      To connect to OTP.

Functions for retrieving data from OTP:

1. `otp_plan()`              To get routes from A to B;
1. `otp_geocode()`           To get the locations of named places e.g. road names.
1. `otp_isochrone()`         To get isochrone maps (OTP 1.x only);
1. `otp_make_surface()`      To make an analyst surface (OTP 1.x only);
1. `otp_surface()`           To evaluate a analyst surface (OTP 1.x only);
1. `otp_traveltime()`        To make a travel time matrix (OTP 1.x only);
1. `otp_surface_isochrone()` To make a raster isochrone map (OTP 1.x only);

Results are returned as [sf objects](https://CRAN.R-project.org/package=sf).


## Acknowledgement

This package was built off the [tutorial by Marcus Young](https://github.com/marcusyoung/otp-tutorial).

## Contribution

Please note that the `opentripplanner` project is released with a [Contributor Code of Conduct](https://github.com/ropensci/opentripplanner/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms. Bug reports and comments are welcome as Github [Issues](https://github.com/ropensci/opentripplanner/issues) and code submissions as [Pull Requests](https://github.com/ropensci/opentripplanner/pulls).


## Citation

Please cite the JOSS paper in publications:

Morgan et al., (2019). OpenTripPlanner for R. Journal of Open Source Software,
  4(44), 1926, https://doi.org/10.21105/joss.01926



[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Advanced Features"
author: "Malcolm Morgan"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{opentripplanner-advanced}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

The vignette introduces some of the more advanced features of OTP and gives some examples of the types of analysis that are possible when using OTP and R together.

### Recap

For this vignette, we will use the same data as the [Getting Started vignette](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html) vignette. If you have not yet created the example graph you can set it up with the following commands. If you are using non-default settings see the Getting Started vignette for full details.

```{r eval =FALSE}
library(opentripplanner)
# Path to a folder containing the OTP.jar file, change to where you saved the file.
path_data <- file.path(tempdir(), "OTP")
dir.create(path_data)
path_otp <- otp_dl_jar()
otp_dl_demo(path_data)
# Build Graph and start OTP
log1 <- otp_build_graph(otp = path_otp, dir = path_data)
log2 <- otp_setup(otp = path_otp, dir = path_data)
otpcon <- otp_connect(timezone = "Europe/London")
```

## Batch Routing

The `otp_plan()` function can produce multiple routes at once. In this example, we will gather data on travel times between each of the [LSOAs](https://www.ons.gov.uk/methodology/geography/ukgeographies/censusgeography) on the Isle of White and the [Ryde Ferry](https://www.wightlink.co.uk/plan-your-journey/routes).

`otp_plan()` accepts three types of input for the `fromPlace` and `toPlace`: a numeric longitude/latitude pair; a 2 x m matrix where each row is a longitude/latitude pair; or an SF data.frame of only POINTS. The number of `fromPlace` and `toPlace` must be the same or equal one (in which case `otp_plan()` will repeat the single location to match the length of the longer locations.

We'll start by importing the locations of the LSOA points.

```{r, eval=FALSE}
lsoa <- sf::st_read("https://github.com/ropensci/opentripplanner/releases/download/0.1/centroids.gpkg",
                    stringsAsFactors = FALSE)
head(lsoa)
```

Then we will define our destination as the Ryde Ferry:

```{r, eval=FALSE}
toPlace <- c(-1.159494,50.732429)
```

Now we can use the `otp_plan()` to find the routes

```{r, eval=FALSE}
routes <- otp_plan(otpcon = otpcon,
                    fromPlace = lsoa,
                    toPlace = toPlace)
```

You may get some warning messages returned as OTP is unable to find some of the routes. The `otp_plan()` will skip over errors and return all the routes it can get. It will then print any messages to the console. You will have also noticed the handy progress bar.

You can plot the routes using the `tmap` package.

If you do plot all the routes it should look something like this:

```{r, echo=FALSE, fig.align='center', fig.cap="\\label{fig:route2airport}Driving Routes to Ryde Ferry"}
knitr::include_graphics("images/routes_to_ferry.jpg")
```

### All to All routing

It is sometimes useful to find the route between every possible origin and destination for example when producing an Origin-Destination (OD) matrix. If you wished to route from every LSOA to every other LSOA point this can easily be done by repeating the points. 

```{r, eval=FALSE}
toPlace   = lsoa[rep(seq(1, nrow(lsoa)), times = nrow(lsoa)),]
fromPlace = lsoa[rep(seq(1, nrow(lsoa)), each  = nrow(lsoa)),]
```

**Warning** routing from all points to all other point increases the total number of routes to calculate exponentially. In this case, 89 points results in 89 x 89 = 7921 routes, this will take a while.

For an OD matrix, you may only be interested in the total travel time and not require the route geometry. By setting `get_geometry = FALSE` in `otp_plan()` R will just return the meta-data and discard the geometry. This is slightly faster than when using `get_geometry = TRUE` and uses less memory.

For example to make a travel time matrix:

```{r, eval=FALSE}
routes <- otp_plan(otpcon = otpcon,
                   fromPlace = fromPlace,
                   toPlace = toPlace,
                   fromID = fromPlace$geo_code,
                   toID = toPlace$geo_code,
                   get_geometry = FALSE)
routes <- routes[,c("fromPlace","toPlace","duration")]
# Use the tidyr package to go from long to wide format
routes_matrix <- tidyr::pivot_wider(routes, 
                               names_from = "toPlace", 
                               values_from = "duration") 
```

Notice the use of `fromID` and `toID` this allows `otp_plan` to return the LSOA `geo_code` with the routes. This can be useful when producing many routes. If no IDs are provided `otp_plan` will return the latitude/longitude of the fromPlace and toPlace.

### Multicore Support

OTP supports multicore routing out of the box. This is based on one core per route, so is only suited to finding a large number of routes. The `otp_plan()` function has the argument `ncores = 1` this can be changed to any positive integer to enable multicore processing e.g. `ncores = 4`. It is recommended that the maximum value for `ncores` is one less than the number of cores on your system. This allows one core to be left for the operating system and any other tasks. 

This graph demonstrates the reduction in time taken to route between all LSOA pairs on the Isle of Wight demo, using one to six cores.

```{r, echo=FALSE, fig.align='center', fig.cap="\\label{fig:multicore} Multicore performance improvements"}
knitr::include_graphics("images/multicore.jpeg")
```

### Distance Balancing

When using multicore routing in `otp_plan` you can optionally set `distance_balance = TRUE`. Distance Balancing sorts the routes by decreasing euclidean distance before sending them to OTP to route. This results in more efficient [load balancing](https://en.wikipedia.org/wiki/Load_balancing_(computing)) between the cores and thus a small reduction in routing time (around five percent). As the original order of the inputs is lost `fromID` and `toID` must be specified to use distance balancing.




## Elevation Profiles

For walking and cycling routes the hilliness of the route matters. If elevation data is available OTP will return the elevation profile of the route. By default, OTP returns the elevation separately from the XY coordinates, but for convenience `otp_plan()` has the argument `get_elevation` which matches the Z coordinates to the XY coordinates. This may result in some minor misalignments. To demonstrate this, let's get a walking route.

```{r, eval=FALSE}
route <- otp_plan(otpcon = otpcon,
                    fromPlace = c(-1.18968, 50.60096),
                    toPlace = c(-1.19105, 50.60439),
                    mode = "WALK",
                    get_elevation = TRUE
                    full_elevation = TRUE)
```

Notice the use of `full_elevation = TRUE` this will return the raw elevation profile from OTP.

We can view the raw profile. It is a data.frame of 3 columns, `first` is the distance along a leg of the route, `second` is the elevation, and `distance` is calculated by `otp_plan()` as the cumulative distance along the whole route. 

As of version 0.3.0.0 the `get_elevation` argument in `otp_plan` is set to FALSE by default, this speeds up routing by only returning XY coordinates rather than XYZ coordinates. 

```{r, eval=FALSE}
profile_raw <- route$elevation[[1]]
plot(profile_raw$distance, profile_raw$second, type = "p",
     xlab = "distance along route", ylab = "elevation")
```

```{r, echo=FALSE, fig.align='center', fig.cap="\\label{fig:ele1}Elevation profile from raw data"}
knitr::include_graphics("images/elevation1.png")
```

To get an elevation profile from the XYZ coordinates is a little more complicated. The `sf::st_coordinates` function returns a matrix of the XYZ coordinates that make up the line. The `geodist` package provides a quick way to calculate the lengths in metres between lng/lat points. 

```{r, eval=FALSE}
profile_xyz <- sf::st_coordinates(route)
dists <- geodist::geodist(profile_xyz[,c("X","Y")], sequential = TRUE)
dists <- cumsum(dists)
plot(dists, profile_xyz[2:nrow(profile_xyz),"Z"], type = "p",
     xlab = "distance along route", ylab = "elevation")
```

```{r, echo=FALSE, fig.align='center', fig.cap="\\label{fig:ele2}Elevation profile from XZY coordinates"}
knitr::include_graphics("images/elevation2.png")
```

Notice that there is less detail in the XYZ graph as the Z coordinates are only matched to a change in XY coordinates, i.e. you only check the elevation when there is a turn in the road.

## Isochrones

Isochrones are lines of equal time. Suppose we are interested in visualising how long it takes to access Ryde ferry using public transport from different parts of the island. We will do this by requesting isochrones from OTP for 15, 30, 45, 60, 75 and 90 minutes. This can be achieved with a single function `otp_isochrone()`.

**Note** isochrones are currently only reliable for "Walk" and "Transit" modes.

```{r eval=FALSE}
ferry_current  <- otp_isochrone(otpcon = otpcon,
            fromPlace = c(-1.159494, 50.732429), # lng/lat of Ryde ferry
            mode = c("WALK","TRANSIT"),
            maxWalkDistance = 2000,
            date_time = as.POSIXct(strptime("2018-06-03 13:30", "%Y-%m-%d %H:%M")),
            cutoffSec = c(15, 30, 45, 60, 75, 90) * 60 ) # Cut offs in seconds
ferry_current$minutes = ferry_current$time / 60 # Convert back to minutes

```

We can visualise the isochrones on a map using the `tmap` package.

```{r, eval=FALSE}
library(tmap)                       # Load the tmap package
tmap_mode("view")                   # Set tmap to interative viewing
map <- tm_shape(ferry_current) +  # Build the map
  tm_fill("minutes",
          breaks = c(0, 15.01, 30.01, 45.01, 60.01, 75.01, 90.01),
          style = "fixed",
          palette = "-RdYlBu") +
  tm_borders()
map                                 # Plot the map
```

You should see a map like this.

```{r, echo=FALSE, fig.align='center', fig.cap="\\label{fig:otpgui}Isochrones from Ryde ferry"}
knitr::include_graphics("images/isochrone.jpg")
```


## Geo-coding

OTP has a built in geo-coder to allow you to search for places by names.
```{r, eval=FALSE}
stations <- otp_geocode(otpcon = otpcon, query = "station")
head(stations)
```

## Debug Layers

For troubleshooting routing issues, you can visualise the traversal permissions of street edges, the bike safety of edges, and how transit stops are linked to streets. For these additional debug layers to be available, add `?debug_layers=true` to the URL, like this:  `http://localhost:8080?debug_layers=true`. The extra layers will be listed in the layer stack menu. 

You can read more about the different debug layers in the official OTP documentation: [http://docs.opentripplanner.org/en/latest/Troubleshooting-Routing/#debug-layers](http://docs.opentripplanner.org/en/latest/Troubleshooting-Routing/#debug-layers).

## Analyst

Older versions of OTP has some limited analytical features built-in which needed to be enabled during graph build and startup. These features are accessible via the `analyst = TRUE` arguments of `otp_build_graph()` and `otp_setup()`. For more information see the [OTP documentation](https://docs.ropensci.org/opentripplanner/articles/Analyst.html)

## Configuring OpenTripPlanner

How OTP works can be configured using JSON files. `build-config.json` is used during graph building (i.e. `otp_build_graph()`). While `router-config.json` is used during setup (i.e. `otp_setup()`). These files must be saved with the rest of your data and each router can have a unique configuration.

To help configure OTP there are several useful functions. `otp_make_config()` makes a default config object and fills it with default values. It is simply a named list, so you can easily modify the values. `otp_validate_config()` does basic checks on a config object to make sure it is valid. Finally `otp_write_config()` exports the config object as a properly formatted JSON file.

A simple example of changing the default walking speed.

```{r, eval=FALSE}
router_config <- otp_make_config("router")     # Make a config object
router_config$routingDefaults$walkSpeed        # Currently 1.34 m/s
router_config$routingDefaults$walkSpeed <- 1.5 # Increase the walking speed
otp_validate_config(router_config)             # Check the new config is valid
otp_write_config(router_config,                # Save the config file
                 dir = path_data,
                 router = "default")

```

There is much more information about configuring OpenTripPlanner at https://opentripplanner.readthedocs.io/en/latest/Configuration/ 

## Running an OTP instance in Docker
We have been able to run OTP version 1.5.0 from `https://repo1.maven.org/maven2/org/opentripplanner/otp/` in a Dockerfile and query it via the package. The basic Dockerfile looks like:

```Docker

FROM java:8-alpine

ENV OTP_VERSION=1.5.0
ENV JAVA_OPTIONS=-Xmx1G

ADD https://repo1.maven.org/maven2/org/opentripplanner/otp/$OTP_VERSION/otp-$OTP_VERSION-shaded.jar /usr/local/share/java/otp.jar

COPY otp /usr/local/bin/

EXPOSE 8080

ENTRYPOINT ["otp"]
CMD ["--help"]
```
And then you can build the image using a default Docker build command like `docker build -t <name> .` where of course "." is your working directory with your Dockerfile in. Then just running the instance like

```sh
docker run \
  -p 8080:8080 \
  -v $PWD/graphs:/var/otp/graphs \
  -e JAVA_OPTIONS=-Xmx4G \
  urbica/otp --server --autoScan --verbose
```
That of course let us place our graphs in the docker volume `$PWD/graphs`. This is slightly edited version of the work described [here](https://github.com/urbica/docker-otp).
---
title: "Analyst"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Analyst}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The Analyst is a set of optional features in OTP 1.x used for transport analysis. These features were removed from OTP 2.0.  This vignette will explain how to enable and use these features.

## Analyst Features

The Analyst adds the following features to OTP:

1. **PointSets** - A more efficient way to specify a batch of fromPlace / toPlaces for routing.
1. **Surfaces** - A efficient way to get travel time from one to many places.
1. **Isochrones** - An alternative implementation of the Isochrone features.

## Loading OTP with the Analyst

The analyst is an optional feature in OTP 1.x and can be loaded in `otp_setup`.

```{r eval=FALSE}
# Get OTP and the Demo Data in the normal way
library(opentripplanner)
path_data <- file.path(tempdir(), "OTP")
dir.create(path_data)
path_otp <- otp_dl_jar(version = "1.5.0") #Must use OTP 1.x
otp_dl_demo(path_data)
log1 <- otp_build_graph(otp = path_otp, dir = path_data)

# Setup OTP with the analyst and pointsets enabled
log2 <- otp_setup(otp = path_otp, dir = path_data, analyst = TRUE, pointsets = TRUE)
otpcon <- otp_connect(timezone = "Europe/London")
```
You can see the analyst is working in the web UI.

```{r, echo = FALSE, fig.align='center', fig.cap="\\label{fig:analystui}OTP Web Interface in Analyst Mode"}
knitr::include_graphics("images/analyst.jpg")
```
## Creating a Pointset

You will have noticed that when we started OTP we also enabled the Pointsets feature. Pointsets are `.csv` files of latitude/longitude coordinates that OTP can route to. PointSets can be created from SF objects using the `otp_poinset()` function.

We will use the LSOA points from the [Advanced Features Vignette](https://docs.ropensci.org/opentripplanner/) to create a PointSet. We will also add a column of data called jobs to analyse. For any numerical values in a PointSet OTP can provide the count and sum of that value based on travel criteria. In this case we might want to know how many jobs can be accessed within a 30 minute drive of a given location.

```{r, eval=FALSE}
lsoa <- sf::st_read("https://github.com/ropensci/opentripplanner/releases/download/0.1/centroids.gpkg", stringsAsFactors = FALSE)
lsoa$jobs <- sample(100:500, nrow(lsoa))
otp_pointset(lsoa, "lsoa", path_data)
```

Here the `otp_pointset()` function takes the `lsoa` object and creates a PointSet that can be accessed by OTP using the name "lsoa".

## Creating a Surface

A surface evaluates the travel times from a given point to all locations within 120 minutes in one minute increments. So it quite similar to an isochrone.

The fist step is to make a surface from a given location using `otp_make_surface`

```{r eval=FALSE}
surfaceid <- otp_make_surface(otpcon, c(-1.17502, 50.64590))
```

The `otp_make_surface` function returns a list containing the surface ID and the parameters used to create the surface. You can use this object in other functions that support surfaces.

Once the surface has been made we can evaluate the travel times to our pointset.

```{r eval=FALSE}
ttimes <- otp_surface(otpcon, surfaceid, "lsoa")
```

The `otp_surface` function returns a list of two objects. The first object `data` is a data frame summarising the sum and count of any values in the pointset in one minute increments from the `fromPlace` specified by the `surfaceid` object. The second object `times` is a vector of travel times in seconds to each of the locations in the pointset. NA values are returned if OTP was unable to find a route to the point or if the point is more than 120 minutes from the `fromPlace`.

If you are not using the sum and count features you can specify `get_data = FALSE` to get just the travel times.

```{r eval=FALSE}
ttimes <- otp_surface(otpcon, surfaceid, "lsoa", get_data = FALSE)
```

We can visualise the travel times using the `tmap` package.

```{r eval=FALSE}
lsoa$time <- ttimes$times / 60
library(tmap)
tmap_mode("view")
tm_shape(lsoa) +
  tm_dots(col = "time")
```

```{r, echo = FALSE, fig.align='center', fig.cap="\\label{fig:ttimes}Travel times to LSOA points"}
knitr::include_graphics("images/ttimes.jpg")
```

Notice that there are a few locations where OTP has failed to find a travel time. It is likely that point is too far from the road network to route, or the `mode` does not allow access (e.g. driving up a pedestrianised road). You could try modifying your point locations slightly to allow OTP to find a valid route.

## Producing a travel time matrix

We can use surfaces as a more efficient way to produce a travel time matrix than using `otp_plan`. A special function `otp_traveltime` creates batches of surfaces and pointsets and summarises the results as a travel time matrix.

```{r eval=FALSE}
ttmatrix <- otp_traveltime(otpcon, 
                           path_data, 
                           fromPlace = lsoa,
                           toPlace = lsoa,
                           fromID = lsoa$geo_code,
                           toID = lsoa$geo_code)
```
The function will return the travel times in seconds between the fromPlaces and toPlaces as a data frame. With columns names from `fromID` and row names from `toID`. Any failed routes will return and NA, as will any greater than 120 minutes.

## Producing an isochrone from a surface

The surfaces can also be used to produce isochrones from a surface.

```{r eval=FALSE}
isochone <- otp_surface_isochrone(otpcon, surfaceid)
```

The function `otp_surface_isochrone` returns a raster image with the value of the travel time in minutes. You can view the raster using the `raster` package.

```{r eval=FALSE}
library(raster)
plot(isochone)
```
```{r, echo = FALSE, fig.align='center', fig.cap="\\label{fig:ttimes}Isochone Raster"}
knitr::include_graphics("images/raster.jpg")
```
---
title: "OTP OpenTripPlanner Version 2.0"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{OTPv2}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

The `opentripplanner` package was originally developed to support OpenTripPlanner v1.3, and has been updated to support each subsequent release. OpenTripPlanner v2.0 is the next major release, and this vignette will document the changes within the package to support OpenTripPlanner v2.0. While most changes for v2.0 are improvements some, like the removal of the isochrone feature are not. So the package will continue to support both version 1.5 (the last v1.x) and subsequent v2.x releases. This vignette will be updated to help users select the best version for their needs.

**Major Changes in 2.0**

* Switch from Java 8 to Java 11
* Support for Netex transit data
* Support for SIRI transit data
* Switch from A* routing algorithm to Range Raptor when searching transit routes
* Removal of Isochrone, Geocode, and other analysis features

## Getting Java 11
At the time of writing Java 8 is still the default on may systems. To get Java 11 you must download the Java Development Kit https://www.oracle.com/java/technologies/javase-jdk11-downloads.html this requires a free account with the Oracle website.

On windows you will also need to change the PATH variable to point to the new version of Java

## Checking your Java version

You can check the version of Java accessible to R by using `otp_check_java` and specifying the version of OTP you want to use.

For example if you wanted to use OTP 2.0 but have Java 8 installed then you would see:

```{r, eval = FALSE}
otp_check_java(2)
[1] FALSE
Warning message:
In otp_check_java(2) : You have OTP 2.x but the version of Java for OTP 1.x

```

## Building a graph and starting up OTP 2

You can select OTP 2.0 using the `otp_dl_jar` function.


```{r setup, eval=FALSE}
library(opentripplanner)
# Path to a folder containing the OTP.jar file, change to where you saved the file.
path_data <- file.path(tempdir(), "OTP")
dir.create(path_data)
path_otp <- otp_dl_jar(version = "2.0.0")
otp_dl_demo(path_data)
# Build Graph and start OTP
log1 <- otp_build_graph(otp = path_otp, dir = path_data)
log2 <- otp_setup(otp = path_otp, dir = path_data)
otpcon <- otp_connect(timezone = "Europe/London")
```
We can download the demo data and build a graph in the usual way.



---
title: "Prerequisites"
author: "Malcolm Morgan"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{opentripplanner-prerequisites}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


## Introduction

This vignette provides guidance for people who new to OpenTripPlanner (OTP) and R.
After you have R and OTP installed and have become familiar with them (e.g. by reading this vignette), we recommend working through [Getting Started vignette](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html).

### Why would I want to use OpenTripPlanner for R?

[OpenTripPlanner](http://www.opentripplanner.org/) (OTP) is a free, open-source, and cross-platform multi-modal route planner written in Java.  It is like having your own private version of Google Maps. OTP can be used to run public route-finding services such as https://ride.trimet.org but it can also be run on your own computer or server. If you want to analyse a transport networks OTP is a very useful tool. However, OTP is, as the name suggests, a Trip Planner, not analytical software.  So, while OTP can find the fastest route from A to B or find all the places that are within a 20-minute walk of C, it cannot answer a question like “how many people live within 10-minutes of a park?”  This is where R can help.  R can process multiple spatial datasets such as population densities and park locations but does not have a built-in journey planner. 

This package allows you to use the trip planning power of OTP with the analytical might of R.

## What is this package for, and what is a package anyway?

The `opentripplanner` R package makes it easier for R and OpenTripPlanner to communicate.
Specifically, it allows you to do use R to control OTP and use it as a local routing service.
For more on local versus remote routing services, see the Transportation chapter in [Geocomputation with R](https://geocompr.robinlovelace.net/transport.html).

### What are R and RStudio?

[R](https://www.r-project.org/) is an open-source programming language and free software environment for statistical computing and graphics.  R has many capabilities for analysing data and writing software, but in this context, its ability to produce and analyse spatial data and maps is most relevant. [RStudio](https://www.rstudio.com/) is an Integrated Development Environment (IDE) for R which is free for personal use.

### What is an R package?

An R package is a small piece of software that extends the basic capabilities of R. It is a bit like how a new phone can do some things out of the box (make phone calls, send email) but you have to install apps to add extra abilities.

### Help with R

To get started with R, see [An Introduction to R](https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf) or from within R by running `help.start()`, introductory tutorials such as DataCamp's free Introduction to R course, or the [R tutor website](http://www.r-tutor.com/r-introduction). The [Geocomputation with R](https://geocompr.robinlovelace.net) book covers the packages and skills required to analyse spatial datasets such as those produced by OpenTripPlanner.

A video tutorial for installing R and RStudio 

* On [Mac](https://www.youtube.com/watch?v=cX532N_XLIs) 
* On [Windows](https://www.youtube.com/watch?v=9-RrkJQQYqY)


## What do I need to get started?

You will need a modern computer and to install some software to use OTP and R on your local computer.

### Hardware

Running your own trip planner is computationally expensive an so the best results will be had on a modern desktop computer. This package comes with some demonstration data for the Isle of Wight, which as a very small area will run on most modern laptops and desktops.

**Recommended Hardware**

* Modern CPU (2010 or later)
* 64 Bit Operating System (Windows/Mac/Linux)
* 8+ GB RAM
* 64 GB+ Hard Disk

**Minimum Hardware**

* 32 Bit Operating System (Windows/Mac/Linux)
* 2 GB RAM

OTP requires at least 1 GB of space to run the demonstration dataset, but by default will request 2 GB. On low-end machines, it is necessary to change the default memory allocations and minimise memory use by other programmes. 

For larger areas OTP will need more memory. An approximate guide to memory use is:

* 2 GB - Small Town / Region
* 4 GB - Large Town / City
* 8 GB - Region / Very Large City
* 20 GB - Country
* 50 GB - Continent

OTP is optimised for city-scale routing and performance will degrade with larger areas, although OTP is used successfully with several small European Countries e.g. the Netherlands.

**Note** If you use a 32 Bit version of Java the maximum amount of memory that can be used by OTP is 4 GB.

### Software

The OpenTripPlanner for R package requires:

* R - [download R](https://cran.r-project.org/mirrors.html) selecting your country (or nearest available country)
* RStudio - [download free version](https://www.rstudio.com/products/rstudio/download/) for personal use. RStudio is not essential but is strongly recommended.
* Java 8 **Note** OTP does not work with later versions of Java such as Java 11.

Download and install Java 8. You will need the correct version for your operating system. If possible the 64 Bit version of Java is preferable, especially if you want to use OTP over large areas.

The package includes a simple function for checking if you have the correct version of Java `otp_check_java()`.

To get Java:

* **Windows/Mac** https://www.java.com/en/download/ This link defaults to the 32 Bit version, so also check https://java.com/en/download/manual.jsp.

* **Linux**  we recommend instructions at [StackOverflow](https://askubuntu.com/questions/740757/switch-between-multiple-java-versions).

For Debian based Linux including Ubuntu and Linux Mint, the following commands in the terminal will install the correct version.

```{eval = FALSE}
sudo apt install openjdk-8-jdk
```

## Next Steps

Now that you have installed Java and R go to the [Getting Started vignette](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html), to find how to install the package, create a graph, and use it to plan trips.
---
title: "Known Issues"
author: "Malcolm Morgan"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{opentripplanner-known-issues}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This vignette highlights some known limitations or bugs with OpenTripPlanner as well as suggested solutions.

## Getting help and reporting bugs

### OpenTripPlanner

* [OTP Users Group](https://groups.google.com/forum/#!forum/opentripplanner-users)
* [OTP GitHub Issues](https://github.com/opentripplanner/OpenTripPlanner/issues)

### OpenTripPlanner for R

* [OTP GitHub Issues](https://github.com/ropensci/opentripplanner/issues)

## Reasons for Graph building to Fail

The `otp_build_graph()` function returns the log from OTP, so the reasons for failures should be documented in the log. Common problems include.

### Can't allocate memory

If you get the error:

`Invalid maximum heap size: The specified size exceeds the maximum representable size.`

Then the `memory` argument in `otp_build_graph()` is too large. This may be because you have set it to more than the amount of RAM on your computer or because you have the 32 bit version of Java rather than the 64 Bit version installed.

### Ran out of memory

A graph build can fail due to lack of memory. You can increase the memory allocated using the `memory` argument in `otp_build_graph()`.

### Bad input data

If you have multiple input files (e.g. GTFS timetables, elevation, config files) and your graph build is failing, try removing all but the PBF file and seeing if the graph can now successfully build. If it does, add files one at a time until you find the file that is causing the build to fail.

## Reasons for Routing to Fail

If you find OTP can not find a route here are some common reasons to check:

### Start or End is too far from the road network

OTP will snap fromPlace and toPlace coordinates to the road network, but only for a limited distance. If your points are far from the road network (e.g. in a lake or middle of a park) then OTP will fail to find a route.

### Mode Specific Limits

OTP does not support all mode combinations (e.g. walk + drive) so some places may only be accessible by certain modes. For example, you can't walk on a motorway or drive on a path. Use the debug layers to check for mode restrictions.

### Maximum walk to transit

By default, OTP caps the maximum walking distance to a transit stop at a low level, so some areas are unreachable by transit. Increase the maximum walking distance to get better results.

### Driving on Roads with cycleway tag

A known bug that stops driving on any road with cycling infrastructure. https://groups.google.com/forum/#!topic/opentripplanner-users/BOv1J32k6Sc 

### Multicore instability

If you get the error:

`Error in unserialize(socklist[[n]]) : error reading from connection`

It means that one of the parallel workers has crashed. Try running your routing in smaller batches like this:

This code will break the routing into 5 batches and put the results into a list called `routes_list` before binding them together into the finished `routes` data frame.

```
routes_list <- list()
n <- split(1:nrow(fromPlace), cut(1:nrow(fromPlace), 5))

for(i in 1:5){
  vals <- n[[i]]
  routes_sub <- otp_plan(otpcon,
                        fromPlace = fromPlace[vals,],
                        toPlace = toPlace[vals,])
  routes_list[[i]] <- routes_sub
}

routes <- dplyr::bind_rows(routes_list)

```

## Speeding up routing

If you are doing a large amount of routing consider the following options.

1. Use multicore routing using `ncores` argument of `otp_plan()`
2. Reduce the area, for example, if you are only routing in one city you don't need a whole country graph.
3. Remove unneeded roads using [osmfilter](https://wiki.openstreetmap.org/wiki/Osmfilter) or [OSMtools](https://github.com/ITSLeeds/OSMtools) for example if you are only interested in driving you can remove footpaths and cycleways. This will simplify your graph and provide a small performance boost.
4. Upgrading your computer, a new fast CPU with more cores will route faster.


## GTFS Data for the UK

See https://github.com/ITSLeeds/UK2GTFS

## Elevation Data

It is common for GeoTIFF to have a no data value often the maximum possible value. OTP can misinterpret this as an elevation value. So set your no data values in your elevation data to something more plausible like 0.

Note that OTP does not support all types of GeoTIFF compression so you may have to change the compression type of the image if you are experiencing problems. 
---
title: "opentripplanner: getting started"
author: "Marcus Young, modified by Malcolm Morgan and Robin Lovelace"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{opentripplanner-get-started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

This vignette is an introduction to OpenTripPlanner (OTP) - an open-source and cross-platform multi-modal route planner written in Java. It uses imported OpenStreetMap (OSM) data for routing on the street and path network and supports
multi-agency public transport routing through imported  [General Transit Feed Specification](https://developers.google.com/transit/gtfs/) (GTFS) feeds. It can also apply a digital elevation model to the OSM street network, allowing, for example, cycle-friendly routes to be requested. OTP has a web front-end that can be used by end-users and a sophisticated routing API.

OTP works worldwide as long as there is OSM map coverage in the area. Support for transit timetables and elevation calculations are dependant on the available data. For more information on getting data see the [getting data for OTP](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html#getting-data-for-otp) section.

A significant advantage of running your own multi-modal route planner is the ability to carry out analysis using amended transport data. Services such as Google Maps or TransportAPI are based on current public transport schedules and the existing road network. OTP enables you to modify transit schedules and/or make changes to the underlying street network. By editing a local copy of OSM, you can model the effects of opening new roads, closing roads or imposing other restrictions. You can also look back in time. For example, you might want to examine the effect of reductions in rural bus services on the accessibility of health facilities. To do this, you would need a network with bus schedules as they were in previous years.

## Prerequisites

You will need to have installed R, RStudio, and Java 8. For more details on the required software, see the [prerequisites vignette](https://docs.ropensci.org/opentripplanner/articles/prerequisites.html) included with this package.

## Installation

Install the stable version of the package from CRAN as follows:

```{r installCRAN, eval=FALSE}
install.packages("opentripplanner") # Install Package
library(opentripplanner)            # Load Package
```

Or you can install the development version from GitHub using the **remotes** package:

```{r installGitHub, eval=FALSE}
remotes::install_github("ropensci/opentripplanner") # Install Package
library(opentripplanner)                            # Load Package
```

## Downloading OTP and the demonstration data

We will build an example graph for the Isle of Wight using some example data provided for the package. A graph is what OTP uses to find routes, and must be built out of the raw data provided. Please note the demo data has been modified for teaching purposes and should not be used for analysis.

### Creating a folder to store OTP and its data

OTP expects its data to be stored in a specific structure, see [building your own otp graph](https://docs.ropensci.org/opentripplanner/articles/opentripplanner.html#building-your-own-otp-graph) for more details. We will make a folder called OTP in the temporary directory. 

```{r, eval=FALSE}
path_data <- file.path(tempdir(), "OTP")
dir.create(path_data) 
```

**You may wish to change this** to keep your files after closing R. Otherwise you will need to download your files and build your graph every time. For example:

```{r, eval=FALSE}
path_data <- file.path("C:/Users/Public", "OTP")
dir.create(path_data) 
```

### Downloading OTP

As of version 0.3 the `otp_dl_jar` function will download and cache the OTP jar file. So you can simply get the path to the OTP JAR file like this:

```{r, eval=FALSE}
path_otp <- otp_dl_jar()
```
`otp_dl_jar` will first check its internal cache for the JAR file and only download it if it is unavailable.

If you wish to specify the location the JAR file is saved (as in version 0.2)  use:

```{r, eval=FALSE}
path_otp <- otp_dl_jar(path_data, cache = FALSE)
```

### Downloading Example Data

The `otp_dl_demo` function downloads the demonstration data and puts it in the correct structure.

```{r eval=FALSE}
otp_dl_demo(path_data)
```

## Building an OTP Graph

Now we can build the graph. This code will create a new file `Graph.obj` that will be saved in the location defined by `path_data`. 

```{r eval=FALSE}
log1 <- otp_build_graph(otp = path_otp, dir = path_data) 
```

By default, R will assign OTP 2GB of memory to build the graph. For larger areas, you may need more memory. You can use the `memory` argument to set the memory allocation in MB. For example, to allocate 10GB, you would use:

```{r eval=FALSE}
log1 <- otp_build_graph(otp = path_otp, dir = path_data, memory = 10240) 
```

Note that you cannot allocate more memory than you have RAM, and if you use 32 Bit Java you cannot allocate more than 3GB. It is possible to run OTP in just 1GB of memory for very small areas (including the demo dataset). 

## Launch OTP and load the graph

The next step is to start up your OTP server, running the router called 'default'. OTP will load the graph you created into memory, and you will then be able to plan multi-modal routes using the web interface. Run the following command:

```{r, eval = FALSE}
log2 <- otp_setup(otp = path_otp, dir = path_data)
```

OTP has a built-in web server called Grizzly which runs on port 8080 (http) and 8081 (https). If you have another application running on your computer that uses these ports then you will need to specify alternative ports using the `port` and `securePort` options, for example:

```{r, eval = FALSE}
log2 <- otp_setup(otp = path_otp, dir = path_data, port = 8801, securePort = 8802)
```

It should only take a minute or two for OTP to load the graph and start the Grizzly server. If it has worked, you should see the message: `OTP is ready to use`, and R will open your web browser at the OTP.

You can also access the web interface using the URL: `http://localhost:8080`. You can now zoom into the Isle of Wight area and request a route by setting an origin and a destination directly on the map (by right-clicking your mouse). You can specify travel dates, times and modes using the 'Trip Options' window (see Figure \ref{fig:otpgui}). You can change the background map from the layer stack icon at the top right.

```{r, echo = FALSE, fig.align='center', fig.cap="\\label{fig:otpgui}OTP Web Interface"}
knitr::include_graphics("images/otpwebgui.jpg")
```

**Note:** The web interface does not work correctly in Internet Explorer - use Firefox or Chrome.

## Connecting to the OTP from R

Now you have the OTP running on your computer you can let R connect to the OTP. `otp_connect()` creates an OTP connection object which will allow R to connect to the OTP. By default, your local timezone is used, but if you are not in the UK, you need to specify the UK timezone to get the correct results with the demo data.

```{r, eval = FALSE}
otpcon <- otp_connect(timezone = "Europe/London")
```

The connection is created and tested, and a message will be returned, saying if the connection exists or not. 
If you have not used the default settings, such as a different `port` you can specify those settings in `otp_connect()`.

```{r, eval = FALSE}
otpcon <- otp_connect(hostname =  "localhost",
                      router = "default",
                      port = 8801)
```

You can also use the `otp_connect()` function to make a connection to OTP on a remote server, for example [digitransit](https://digitransit.fi/en/) run an OTP server for Finland. See their [documentation](https://digitransit.fi/en/developers/) for more details. Digitransit use a non-standard URL structure for their API, so we provide the full URL directly to `otp_connect()`. See the `otp_connect()` help for more details.

```{r, eval = FALSE}
otpcon <- otp_connect(url = "https://api.digitransit.fi/routing/v1/routers/hsl")
```

## Getting a route from the OTP

Now we can use R to get a route from the OTP. OTP accepts pairs of longitude/latitude coordinates for a `fromPlace` (start of the journey) and `toPlace` (end of the journey).

```{r, eval = FALSE}
route <- otp_plan(otpcon, 
                  fromPlace = c(-1.17502, 50.64590), 
                  toPlace = c(-1.15339, 50.72266))
```

If you have the `tmap` package installed, you can view the route using.

```{r, eval = FALSE}
# install.packages("tmap") # Only needed if you don't have tmap
library(tmap)              # Load the tmap package
tmap_mode("view")          # Set tmap to interactive viewing
qtm(route)      # Plot the route on a map
```


## Stopping the OTP

As the OTP is running in Java, it will continue to run after you close R.

You can stop the OTP running using the command. **NOTE: This will stop all running JAVA applications!**

```{r, eval = FALSE}
otp_stop()
```

Congratulations, you now have your own multi-modal router planner!

## Building your own OTP Graph

If you want to build your own graph for a different location, follow these steps, and change your `path_data` variable to the folder with your data.

An OTP graph specifies every location in the region covered and how to travel between them, and is compiled by OTP using OSM data for the street and path network (used for walking, bicycle, and drive modes) and GTFS data for transit scheduling.

Our first task is to create the folder and file structure expected by OTP. This is a base directory called `otp` which contains a sub-directory called `graphs`. Directories created under `graphs` are known as OTP routers and contain all the files required to build a graph. A single OTP instance can host several routers, for example, covering different regions. 

The basic file structure is shown below.

```{r, engine='bash', eval=FALSE}
/ otp                         # Your top folder for storing all OTP data
  /graphs                     
     /default                 # Subfolder with the name of the router
         osm.pbf              # Required OSM road map
         router-config.json   # Optional config file
         build-config.json    # Optional config file
         gtfs.zip             # Optional GTFS data
         dem.tif              # Optional Elevation data
         
```

`router-config.json` is read when the OTP server is started, while `build-config.json` is read during the graph building stage. For more information on the config files, see the [advanced features vignette](https://docs.ropensci.org/opentripplanner/articles/advanced_features.html).

### Multiple Routers 

OTP supports multiple routes by having them in separate folders. For example:

```{r, engine='bash', eval=FALSE}
/ otp                        # Your top folder for storing all OTP data
  /graphs                     
     /london                 # router called london
         osm.pbf             # map of London
     /manchester             # router called manchester
         osm.pbf             # map of Manchester
```

Having multiple routers may be useful to support different regions, or different scenarios, e.g. past, present, future.

Each router requires its own graph via `otp_build_graph()` but only one version of the OTP jar files is needed (i.e. `otp_dl_jar()`) Each route must be started via `otp_setup()` before use. While it is technically possible to run multiple routers simultaneously, it is not recommended for performance reasons.

### Getting Data for OTP

To use OTP on your local computer, you will need data about your location of interest, such as where the roads are, what is the public transport timetable etc.

#### Road Map Data

OTP uses road maps from the [Open Street Map (OSM)](https://www.openstreetmap.org) in the [.osm or .pbf format](https://wiki.openstreetmap.org/wiki/OSM_file_formats). OSM is a free map that anybody can edit. You can download all or part of the OSM from [Geofabrik](https://download.geofabrik.de/).

Note that OTP only really needs the road network from the OSM, but the Geofabrik downloads include non-road data such as buildings, parks, water etc. So for large areas, you may wish to remove unneeded data from the file before building the graph. You can do this using [osmconvert](https://wiki.openstreetmap.org/wiki/Osmconvert) and [osmfilter](https://wiki.openstreetmap.org/wiki/Osmfilter). These are command-line utilities for editing large OSM files and are documented on the OSM Wiki. For small areas (e.g. city scale) this is not necessary. 

#### Public Transport Timetables

OTP can use public transport timetable data in the [GTFS format](https://developers.google.com/transit/gtfs/). You can find GTFS data for many regions at [transitland](https://www.transit.land/). For users in the UK see [UK2GTFS](https://itsleeds.github.io/UK2GTFS/)

#### Elevation Data

You can add terrain information to your routes, especially useful for walking and cycling, using [GeoTIFF images](https://trac.osgeo.org/geotiff). You can find worldwide elevation data from [NASA](https://www2.jpl.nasa.gov/srtm/).

**Warning** It is common for GeoTIFF to have a no data value often the maximum possible value. OTP can misinterpret this as an elevation value. So set your no data values in your elevation data to something more plausible like 0.

**Note** that OTP does not support all types of GeoTIFF compression so you may have to change the compression type of the image if you are experiencing problems. 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-routing-options.R
\name{otp_routing_options}
\alias{otp_routing_options}
\title{Make routingOptions object}
\usage{
otp_routing_options()
}
\description{
OTP supports a wide selection of routing options `otp_plan()` accepts a
named list of these options. This function produces an empty named list
of valid options supported by both this package and OTP.
}
\details{
Supports almost all of the possible options in OTP 1.4. Note that some
of the most popular option (mode, date, time, etc.) are set directly
in `otp_plan()`. If you want to permenaty set an option many are supported
in the config files, see help on `otp_make_config()`.

http://dev.opentripplanner.org/apidoc/1.4.0/resource_PlannerResource.html
}
\examples{
\dontrun{
routingOptions <- otp_routing_options()
routingOptions$walkSpeed <- 1.5
routingOptions <- otp_validate_routing_options(routingOptions)
}
}
\seealso{
Other routing: 
\code{\link{otp_geocode}()},
\code{\link{otp_isochrone}()},
\code{\link{otp_plan}()},
\code{\link{otp_pointset}()},
\code{\link{otp_validate_routing_options}()}
}
\concept{routing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-connect.R
\name{otp_connect}
\alias{otp_connect}
\title{Set up and confirm a connection to an OTP instance.}
\usage{
otp_connect(
  hostname = "localhost",
  router = "default",
  url = NULL,
  port = 8080,
  ssl = FALSE,
  check = TRUE,
  timezone = Sys.timezone(),
  otp_version = 1.5
)
}
\arguments{
\item{hostname}{A string, e.g.
"ec2-34-217-73-26.us-west-2.compute.amazonaws.com". Optional, default is
"localhost".}

\item{router}{A string, e.g. "UK2018". Optional, default is "default". OTP
can support multiple routers see advanced vignette for details.}

\item{url}{If a non-standard URL structure is used provide a full url,
default is NULL}

\item{port}{A positive integer. Optional, default is 8080.}

\item{ssl}{Logical, indicates whether to use https. Optional, default is
FALSE.}

\item{check}{Logical. If TRUE connection object is only returned if OTP
instance and router are confirmed reachable. Optional, default is TRUE.}

\item{timezone}{Character, timezone, defaults to local timezone}

\item{otp_version}{Numeric, OTP Version in use default 1.5, if check is TRUE
then `otp_check_version()` is called and otp_version is updated}
}
\value{
Returns an S3 object of class otpconnect. If \code{check} is TRUE and
  the router is not reachable the object is not returned.
}
\description{
Defines the parameters required to connect to a router on an OTP instance
and, if required, confirms that the instance and router are query-able.
}
\details{
The default URL structure for the OTP API is:
http://<hostname>:<port>/otp/routers/<router> For example:
http://localhost:8080/otp/routers/default

Functions construct the URL from the parameters provided in otpconnect
objects. However some websites hosting OTP have modified the default URL
structure. If this is the case you can use the \code{url} parameter to bypass
the URL construction and provide a fully formed URL. In this case the
\code{hostname}, \code{router}, \code{port}, and \code{ssl} are ignored.
}
\examples{
\dontrun{
otpcon <- otp_connect()
otpcon <- otp_connect(
  router = "UK2018",
  ssl = TRUE
)
otpcon <- otp_connect(
  hostname = "ec2.us-west-2.compute.amazonaws.com",
  router = "UK2018",
  port = 8888,
  ssl = TRUE
)
otpcon <- otp_connect(
  url = "https://api.digitransit.fi:443/routing/v1/routers/hsl"
)
}
}
\concept{connect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-traveltime.R
\name{otp_traveltime}
\alias{otp_traveltime}
\title{Get travel times between points}
\usage{
otp_traveltime(
  otpcon = NA,
  path_data = NULL,
  fromPlace = NA,
  toPlace = NA,
  fromID = NULL,
  toID = NULL,
  mode = "CAR",
  date_time = Sys.time(),
  arriveBy = FALSE,
  maxWalkDistance = 1000,
  numItineraries = 3,
  routeOptions = NULL,
  ncores = 1,
  timezone = otpcon$timezone
)
}
\arguments{
\item{otpcon}{OTP connection object produced by otp_connect()}

\item{path_data}{Path to data used in otp_build_graph()}

\item{fromPlace}{Numeric vector, Longitude/Latitude pair, e.g.
`c(-0.134649,51.529258)`, or 2 column matrix of Longitude/Latitude pairs,
or sf data frame of POINTS with CRS 4326}

\item{toPlace}{Numeric vector, Longitude/Latitude pair, e.g.
`c(-0.088780,51.506383)`, or 2 column matrix of Longitude/Latitude pairs,
or sf data frame of POINTS with CRS 4326}

\item{fromID}{character vector same length as fromPlace}

\item{toID}{character vector same length as toPlace}

\item{mode}{character vector of one or more modes of travel valid values
TRANSIT, WALK, BICYCLE, CAR, BUS, RAIL, default CAR. Not all combinations
are valid e.g. c("WALK","BUS") is valid but c("WALK","CAR") is not.}

\item{date_time}{POSIXct, a date and time, defaults to current date and time}

\item{arriveBy}{Logical, Whether the trip should depart or arrive at the
specified date and time, default FALSE}

\item{maxWalkDistance}{Numeric passed to OTP in metres}

\item{numItineraries}{The maximum number of possible itineraries to return}

\item{routeOptions}{Named list of values passed to OTP use
`otp_route_options()` to make template object.}

\item{ncores}{Numeric, number of cores to use when batch processing, default
1, see details}

\item{timezone}{Character, what timezone to use, see as.POSIXct, default is
local timezone}
}
\value{
Returns an  data frame
}
\description{
This function requires OTP 1.x and the analyst
}
\details{
Make a travel time matrix using the analyst features in OPT 1.x
}
\seealso{
Other analyst: 
\code{\link{otp_make_surface}()}
}
\concept{analyst}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-surface.R
\name{otp_make_surface}
\alias{otp_make_surface}
\title{Make a Surface}
\usage{
otp_make_surface(
  otpcon = NULL,
  fromPlace = c(-1.17502, 50.6459),
  mode = "CAR",
  date_time = Sys.time(),
  maxWalkDistance = 1000,
  arriveBy = FALSE,
  routeOptions = NULL,
  timezone = otpcon$timezone
)
}
\arguments{
\item{otpcon}{OTP connection object produced by otp_connect()}

\item{fromPlace}{Numeric vector, Longitude/Latitude pair, e.g.
`c(-0.134649,51.529258)`, or 2 column matrix of Longitude/Latitude pairs,
or sf data frame of POINTS with CRS 4326}

\item{mode}{character vector of one or more modes of travel valid values
TRANSIT, WALK, BICYCLE, CAR, BUS, RAIL, SUBWAY, TRAM, FERRY, default CAR.
Not all combinations are valid e.g. c("WALK","BUS") is valid but
c("WALK","CAR") is not.}

\item{date_time}{POSIXct, a date and time, defaults to current date and time}

\item{maxWalkDistance}{Numeric passed to OTP in metres}

\item{arriveBy}{Logical, Whether the trip should depart or arrive at the
specified date and time, default FALSE}

\item{routeOptions}{Named list of values passed to OTP use}

\item{timezone}{Character, what timezone to use, see as.POSIXct, default is
local timezone}
}
\value{
Returns a list with information about the surface created
}
\description{
Requires OTP 1.x and the analyst
}
\details{
THis function requires the analysis and pointset features to be
  enabled during `otp_setup()`. Thus it will only work with OTP 1.x. For more
  detail see the analyst vignette.
}
\examples{
\dontrun{
surface <- otp_make_surface(otpcon, c(-1.17502, 50.64590))
}
}
\seealso{
Other analyst: 
\code{\link{otp_traveltime}()}
}
\concept{analyst}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-setup.R
\name{otp_setup}
\alias{otp_setup}
\title{Set up an OTP instance.}
\usage{
otp_setup(
  otp = NULL,
  dir = NULL,
  memory = 2048,
  router = "default",
  port = 8080,
  securePort = 8081,
  analyst = FALSE,
  pointsets = FALSE,
  wait = TRUE,
  flag64bit = TRUE,
  quiet = TRUE,
  otp_version = NULL,
  open_browser = TRUE
)
}
\arguments{
\item{otp}{A character string, path to the OTP .jar file}

\item{dir}{A character string, path to a directory containing the
necessary files, see details}

\item{memory}{A positive integer. Amount of memory to assign to
the OTP in MB, default is 2048}

\item{router}{A character for the name of the router to use, must
be subfolder of dir/graphs, default "default". See
vignettes for details.}

\item{port}{A positive integer. Optional, default is 8080.}

\item{securePort}{A positive integer. Optional, default is 8081.}

\item{analyst}{Logical. Should the analyst features be loaded?
Default FALSE}

\item{pointsets}{Logical. Should the pointsets be loaded?
Default FALSE}

\item{wait}{Logical, Should R wait until OTP has loaded before
running next line of code, default TRUE}

\item{flag64bit}{Logical, if true the -d64 flag is added to Java instructions,
ignored if otp_version >= 2}

\item{quiet}{Logical, if FALSE the Java commands will be printed to console}

\item{otp_version}{Numeric, version of OTP to build, default NULL when version
is auto-detected}

\item{open_browser}{Logical, if TRUE web browser is loaded when OTP is ready}
}
\value{
This function does not return a value to R. If wait is TRUE R
will wait until OTP is running (maximum of 5 minutes).
After 5 minutes (or if wait is FALSE) the function will return
R to your control, but the OTP will keep loading.
}
\description{
OTP is run in Java and requires Java commands to be typed into the
command line. The function allows the parameters to be defined in
R and automatically passed to Java. This function sets up a local
instance of OTP, for remote versions see documentation.

The function assumes you have run otp_build_graph()
}
\details{
To run an OTP graph must have been created using otp_build_graph
and the following files to be in the directory specified by the
dir variable.

/graphs - A sub-directory

  /default - A sub-directory with the name of the OTP router used in 'router' variable

    graph.obj  OTP graph
}
\examples{
\dontrun{
otp_setup(
  otp = "C:/otp/otp.jar",
  dir = "C:/data"
)
otp_setup(
  otp = "C:/otp/otp.jar",
  dir = "C:/data",
  memory = 5000,
  analyst = TRUE
)
}
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-routing-options.R
\name{otp_validate_routing_options}
\alias{otp_validate_routing_options}
\title{Validate routingOptions object}
\usage{
otp_validate_routing_options(opts)
}
\arguments{
\item{opts}{a named list of options possibly from `otp_routing_options()`}
}
\description{
OTP supports a wide selection of routing options `otp_plan()` accepts a
named list of these options. This function validates a named list of inputs
and removes any empty inputs.
}
\details{
Supports almost all of the possible options in OTP 1.4. Note that some
of the most popular option (mode, date, time, etc.) are set directly
in `otp_plan()`. If you want to permenaty set an option many are supported
in the config files, see help on `otp_make_config()`.
http://dev.opentripplanner.org/apidoc/1.4.0/resource_PlannerResource.html
}
\examples{
\dontrun{
routingOptions <- otp_routing_options()
routingOptions$walkSpeed <- 1.5
routingOptions <- otp_validate_routing_options(routingOptions)
}
}
\seealso{
Other routing: 
\code{\link{otp_geocode}()},
\code{\link{otp_isochrone}()},
\code{\link{otp_plan}()},
\code{\link{otp_pointset}()},
\code{\link{otp_routing_options}()}
}
\concept{routing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{json_example_long_drive}
\alias{json_example_long_drive}
\title{Example JSON for driving long distance}
\format{
json
}
\usage{
json_example_long_drive
}
\description{
Example JSON response from OTP
This is used for internal testing and has no use
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-geocode.R
\name{otp_geocode}
\alias{otp_geocode}
\title{Use OTP Geo-coder to find a location}
\usage{
otp_geocode(
  otpcon = NULL,
  query = NULL,
  autocomplete = FALSE,
  stops = TRUE,
  clusters = FALSE,
  corners = TRUE,
  type = "SF"
)
}
\arguments{
\item{otpcon}{OTP connection object produced by otp_connect()}

\item{query}{Character, The query string we want to geocode}

\item{autocomplete}{logical Whether we should use the query
string to do a prefix match, default FALSE}

\item{stops}{Logical, Search for stops, either by name or
stop code, default TRUE}

\item{clusters}{Logical, Search for clusters by their name,
default FALSE}

\item{corners}{Logical, Search for street corners using at
least one of the street names, default TRUE}

\item{type}{Character, How should results be returned can
be "SF" or "Coordinates" or "Both", Default "SF"}
}
\value{
Returns a data.frame of SF POINTS or Coordinates of all
    the locations that match `query`
}
\description{
Geo-coding converts a named place, such as a street name into a
     lng/lat pair.
}
\details{
OTP will return a maximum of 10 results
}
\examples{
\dontrun{
locations <- otp_geocode(otpcon, "High Street")
}
}
\seealso{
Other routing: 
\code{\link{otp_isochrone}()},
\code{\link{otp_plan}()},
\code{\link{otp_pointset}()},
\code{\link{otp_routing_options}()},
\code{\link{otp_validate_routing_options}()}
}
\concept{routing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-surface.R
\name{otp_surface}
\alias{otp_surface}
\title{Evaluate a surface against a pointset}
\usage{
otp_surface(otpcon = NULL, surface = NULL, pointsset = NULL, get_data = TRUE)
}
\arguments{
\item{otpcon}{OTP connection object produced by otp_connect()}

\item{surface}{A suface list from otp_make_surface()}

\item{pointsset}{character, name of pointset}

\item{get_data}{Logical, should data be returned or just travel times.}
}
\value{
Returns a data.frame of travel times
}
\description{
Geo-coding converts a named place, such as a street name into a lng/lat pair.
}
\details{
THis function requires the analysis and pointset features to be
enabled during `otp_setup()`. Thus it will only work with OTP 1.x. For more
detail see the analyst vignette.
}
\examples{
\dontrun{
times <- otp_surface(otpcon, c(-1.17502, 50.64590), "lsoa", path_data)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-plan.R
\name{otp_plan}
\alias{otp_plan}
\title{Get get a route or routes from the OTP}
\usage{
otp_plan(
  otpcon = NA,
  fromPlace = NA,
  toPlace = NA,
  fromID = NULL,
  toID = NULL,
  mode = "CAR",
  date_time = Sys.time(),
  arriveBy = FALSE,
  maxWalkDistance = 1000,
  numItineraries = 3,
  routeOptions = NULL,
  full_elevation = FALSE,
  get_geometry = TRUE,
  ncores = 1,
  timezone = otpcon$timezone,
  distance_balance = FALSE,
  get_elevation = FALSE
)
}
\arguments{
\item{otpcon}{OTP connection object produced by otp_connect()}

\item{fromPlace}{Numeric vector, Longitude/Latitude pair, e.g.
`c(-0.134649,51.529258)`, or 2 column matrix of Longitude/Latitude pairs,
or sf data frame of POINTS with CRS 4326}

\item{toPlace}{Numeric vector, Longitude/Latitude pair, e.g.
`c(-0.088780,51.506383)`, or 2 column matrix of Longitude/Latitude pairs,
or sf data frame of POINTS with CRS 4326}

\item{fromID}{character vector same length as fromPlace}

\item{toID}{character vector same length as toPlace}

\item{mode}{character vector of one or more modes of travel valid values
TRANSIT, WALK, BICYCLE, CAR, BUS, RAIL, default CAR. Not all combinations
are valid e.g. c("WALK","BUS") is valid but c("WALK","CAR") is not.}

\item{date_time}{POSIXct, a date and time, defaults to current date and time}

\item{arriveBy}{Logical, Whether the trip should depart or arrive at the
specified date and time, default FALSE}

\item{maxWalkDistance}{Numeric passed to OTP in metres}

\item{numItineraries}{The maximum number of possible itineraries to return}

\item{routeOptions}{Named list of values passed to OTP use
`otp_route_options()` to make template object.}

\item{full_elevation}{Logical, should the full elevation profile be returned,
default FALSE}

\item{get_geometry}{Logical, should the route geometry be returned, default
TRUE, see details}

\item{ncores}{Numeric, number of cores to use when batch processing, default
1, see details}

\item{timezone}{Character, what timezone to use, see as.POSIXct, default is
local timezone}

\item{distance_balance}{Logical, use distance balancing, default false, see
details}

\item{get_elevation}{Logical, default FALSE, if true XYZ coordinates returned
else XY coordinates returned.}
}
\value{
Returns an SF data frame of LINESTRINGs
}
\description{
This is the main routing function for OTP and can find single or
  multiple routes between `fromPlace` and `toPlace`.
}
\details{
This function returns a SF data.frame with one row for each leg of
  the journey (a leg is defined by a change in mode). For transit, more than
  one route option may be returned and is indicated by the `route_option`
  column. The number of different itineraries can be set with the
  `numItineraries` variable.

  ## Batch Routing

  When passing a matrix or SF data frame object to fromPlace and toPlace
  `otp_plan` will route in batch mode. In this case the `ncores` variable
  will be used. Increasing `ncores` will enable multicore routing, the max
  `ncores` should be the number of cores on your system - 1.

  ## Distance Balancing

  When using multicore routing each task does not take the same amount of
  time. This can result in wasted time between batches. Distance Balancing
  sorts the routing by the euclidean distance between fromPlace and toPlace
  before routing. This offers a small performance improvement of around five
  percent. As the original order of the inputs is lost fromID and toID must
  be provided.

  ## Elevation

  OTP supports elevation data and can return the elevation profile of the
  route if available. OTP returns the elevation profile separately from the
  XY coordinates, this means there is not direct match between the number of
  XY points and the number of Z points.  OTP also only returns the elevation
  profile for the first leg of the route (this appears to be a bug). If
  `get_elevation` is TRUE the otp_plan function matches the elevation profile
  to the XY coordinates to return an SF linestring with XYZ coordinates. If
  you require a more detailed elevation profile, the full_elevation parameter
  will return a nested data.frame with three columns. first and second are
  returned from OTP, while distance is the cumulative distance along the
  route and is derived from First.

  ## Route Geometry

  The `get_geometry` provides the option to not return the route geometry,
  and just return the meta-data (e.g. journey time). This may be useful when
  creating an Origin:Destination matrix and also provides a small performance
  boost by reduced processing of geometries.
}
\examples{
\dontrun{
otpcon <- otp_connect()
otp_plan(otpcon, c(0.1, 55.3), c(0.6, 52.1))
otp_plan(otpcon, c(0.1, 55.3), c(0.6, 52.1),
  mode = c("WALK", "TRANSIT")
)
otp_plan(otpcon, c(0.1, 55.3), c(0.6, 52.1),
  mode = "BICYCLE", arriveBy = TRUE,
  date_time = as.POSIXct(strptime("2018-06-03 13:30", "\%Y-\%m-\%d \%H:\%M"))
)
}

}
\seealso{
Other routing: 
\code{\link{otp_geocode}()},
\code{\link{otp_isochrone}()},
\code{\link{otp_pointset}()},
\code{\link{otp_routing_options}()},
\code{\link{otp_validate_routing_options}()}
}
\concept{routing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opentripplanner-package.R
\docType{package}
\name{opentripplanner-package}
\alias{opentripplanner}
\alias{opentripplanner-package}
\title{OpenTripPlanner of R}
\description{
The goal of OpenTripPlanner for R is to provide a simple R interface
to OpenTripPlanner (OTP).
The OTP is a multimodal trip planning service.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/opentripplanner}
  \item \url{https://docs.ropensci.org/opentripplanner/}
  \item Report bugs at \url{https://github.com/ropensci/opentripplanner/issues}
}

}
\author{
\strong{Maintainer}: Malcolm Morgan \email{m.morgan1@leeds.ac.uk} (\href{https://orcid.org/0000-0002-9488-9183}{ORCID})

Authors:
\itemize{
  \item Marcus Young \email{M.A.Young@soton.ac.uk} (\href{https://orcid.org/0000-0003-4627-1116}{ORCID})
  \item Robin Lovelace \email{rob00x@gmail.com} (\href{https://orcid.org/0000-0001-5679-6536}{ORCID})
}

Other contributors:
\itemize{
  \item Layik Hama \email{layik.hama@gmail.com} (\href{https://orcid.org/0000-0003-1912-4890}{ORCID}) [contributor]
}

}
\keyword{mulitmodal}
\keyword{opentripplanner}
\keyword{routing}
\keyword{transport}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-setup.R
\name{otp_build_graph}
\alias{otp_build_graph}
\title{Build an OTP Graph}
\usage{
otp_build_graph(
  otp = NULL,
  dir = NULL,
  memory = 2048,
  router = "default",
  flag64bit = TRUE,
  quiet = TRUE,
  otp_version = NULL
)
}
\arguments{
\item{otp}{A character string, path to the OTP .jar file}

\item{dir}{A character string, path to a directory containing the necessary
files, see details}

\item{memory}{A positive integer. Amount of memory to assign to the OTP in
MB, default is 2048}

\item{router}{A character string for the name of the router, must subfolder
of  dir/graphs, default "default". See vignettes for details.}

\item{flag64bit}{Logical, if true the -d64 flag is added to Java instructions,
ignored if otp_version >= 2}

\item{quiet}{Logical, if FALSE the Java commands will be printed to console}

\item{otp_version}{Numeric, version of OTP to build, default NULL when version
is auto-detected}
}
\value{
Character vector of messages produced by OTP, and will return the
  message "Graph built" if successful
}
\description{
OTP is run in Java and requires Java commands to be typed into
  the command line. The function allows the parameters to be defined in R and
  automatically passed to Java. This function builds a OTP graph from the
  Open Street Map and other files.
}
\details{
The OTP .jar file can be downloaded from
  https://repo1.maven.org/maven2/org/opentripplanner/otp/

  To build an OTP graph requires the following files to be in the directory
  specified by the dir variable.

  /graphs - A sub-directory

  /default - A sub-directory with the name of the OTP router used in router'
  variable

  osm.pbf - Required, pbf file containing the Open Street Map

  router-config.json - Required, json file containing configurations settings
  for the OTP

  gtfs.zip - Optional, and number of GTFS files with transit timetables

  terrain.tif - Optional, GeoTiff image of terrain map

  The function will accept any file name for the .jar file, but it must be
  the only .jar file in that directory OTP can support multiple routers (e.g.
  different regions), each router must have its own sub-directory in the
  graphs directory
}
\examples{
\dontrun{
log <- otp_build_graph(otp = "C:/otp/otp.jar", dir = "C:/data")
}
}
\seealso{
Other setup: 
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-surface.R
\name{otp_pointset}
\alias{otp_pointset}
\title{Create a pointset}
\usage{
otp_pointset(points = NULL, name = NULL, dir = NULL)
}
\arguments{
\item{points}{sf data frame of POINTS with CRS 4326}

\item{name}{Character, name for pointset}

\item{dir}{A character string, path to a directory containing the necessary
files, see details}
}
\value{
Returns a data.frame of SF POINTS or Coordinates of all
    the locations that match `query`
}
\description{
Create a pointset
}
\details{
OTP will return a maximum of 10 results
}
\examples{
\dontrun{
locations <- otp_geocode(otpcon, "High Street")
}
}
\seealso{
Other routing: 
\code{\link{otp_geocode}()},
\code{\link{otp_isochrone}()},
\code{\link{otp_plan}()},
\code{\link{otp_routing_options}()},
\code{\link{otp_validate_routing_options}()}
}
\concept{routing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-config.R
\name{otp_make_config}
\alias{otp_make_config}
\title{Make Config Object}
\usage{
otp_make_config(type)
}
\arguments{
\item{type}{Which type of config file to create, "otp", "build", "router"}
}
\description{
OTP can be configured using three json files `otp-config.json`,
`build-config.json`, and `router-config.json`. This function
creates a named list for each config file and
populates the defaults values.
}
\details{
For more details see:
http://docs.opentripplanner.org/en/latest/Configuration
}
\examples{
{
  conf <- otp_make_config("build")
  conf <- otp_make_config("router")
}
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{json_example_drive}
\alias{json_example_drive}
\title{Example JSON for driving}
\format{
json
}
\usage{
json_example_drive
}
\description{
Example JSON response from OTP
This is used for internal testing and has no use
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-download.R
\name{otp_dl_demo}
\alias{otp_dl_demo}
\title{Download Demo Data}
\usage{
otp_dl_demo(
  path_data = NULL,
  url = paste0("https://github.com/ropensci/opentripplanner/",
    "releases/download/0.1/isle-of-wight-demo.zip"),
  quiet = FALSE
)
}
\arguments{
\item{path_data}{path to folder where data for OTP is to be stored}

\item{url}{URL to data}

\item{quiet}{logical, passed to download.file, default FALSE}
}
\description{
Download the demonstration data for the Isle of Wight
}
\examples{
\dontrun{
otp_dl_demo(tempdir())
}
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-download.R
\name{otp_dl_jar}
\alias{otp_dl_jar}
\title{Download OTP Jar File}
\usage{
otp_dl_jar(
  path = NULL,
  version = "1.5.0",
  file_name = paste0("otp-", version, "-shaded.jar"),
  url = "https://repo1.maven.org/maven2/org/opentripplanner/otp",
  quiet = FALSE,
  cache = TRUE
)
}
\arguments{
\item{path}{path to folder where OTP is to be stored}

\item{version}{a character string of the version number default is "1.5.0"}

\item{file_name}{file name to give the otp default "otp.jar"}

\item{url}{URL to the download server}

\item{quiet}{logical, passed to download.file, default FALSE}

\item{cache}{logical, default TRUE, see details}
}
\value{
The path to the OTP file
}
\description{
Download the OTP jar file from maven.org
}
\details{
As of version 0.3.0.0 `otp_dl_jar` will cache the JAR file within
the package and ignore the `path` argument. You can force a new download to
be saved in the `path` location by setting `cache = FALSE`.
}
\examples{
\dontrun{
otp_dl_jar(tempdir())
}
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{json_example_transit}
\alias{json_example_transit}
\title{Example JSON for transit}
\format{
json
}
\usage{
json_example_transit
}
\description{
Example JSON response from OTP
This is used for internal testing and has no use
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-config.R
\name{otp_write_config}
\alias{otp_write_config}
\title{Write config object as json file}
\usage{
otp_write_config(config, dir = NULL, router = "default")
}
\arguments{
\item{config}{A named list made/modified from `otp_make_config()`}

\item{dir}{Path to folder where data for OTP is to be stored}

\item{router}{name of the router, default is "default", must be a subfolder
of dir/graphs}
}
\description{
Takes a config list produced by `otp_make_config()` and saves it
as json file for OTP
}
\examples{
\dontrun{
conf <- otp_make_config("build")
otp_write_config(conf, dir = tempdir())
}
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-isochrone-batch.R
\name{otp_isochrone}
\alias{otp_isochrone}
\title{Get the Isochrones from a location}
\usage{
otp_isochrone(
  otpcon = NA,
  fromPlace = NA,
  fromID = NULL,
  mode = "TRANSIT",
  date_time = Sys.time(),
  arriveBy = FALSE,
  maxWalkDistance = 1000,
  routingOptions = NULL,
  cutoffSec = c(600, 1200, 1800, 2400, 3000, 3600),
  ncores = 1,
  timezone = otpcon$timezone
)
}
\arguments{
\item{otpcon}{OTP connection object produced by otp_connect()}

\item{fromPlace}{Numeric vector, Longitude/Latitude pair,
    e.g. `c(-0.134649,51.529258)`,
or 2 column matrix of Longitude/Latitude pairs, or sf
    data frame of POINTS}

\item{fromID}{character vector same length as fromPlace}

\item{mode}{character vector of one or more modes of travel valid values
TRANSIT, WALK, BICYCLE, CAR, BUS, RAIL, default CAR. Not all
combinations are valid e.g. c("WALK","BUS") is valid but
c("WALK","CAR") is not.}

\item{date_time}{POSIXct, a date and time, defaults to current
date and time}

\item{arriveBy}{Logical, Whether the trip should depart or
arrive at the specified date and time, default FALSE}

\item{maxWalkDistance}{maximum distance to walk in metres}

\item{routingOptions}{named list passed to OTP see `otp_routing_options()`}

\item{cutoffSec}{Numeric vector, number of seconds to define
the break points of each Isochrone}

\item{ncores}{number of cores to use in parallel processing}

\item{timezone}{character, timezone to use, default from otpcon}
}
\value{
Returns a SF data.frame of POLYGONs
}
\description{
Get the Isochrones from a location
}
\details{
Isochrones are maps of equal travel time,
for a given location a map is produced showing how long it takes to reach
each location.

Isochrones are only available from OTP v1.x and will not work with v2.0
}
\examples{
\dontrun{
isochrone1 <- otp_isochrone(otpcon, fromPlace = c(-0.1346, 51.5292))
isochrone2 <- otp_isochrone(otpcon,
  fromPlace = c(-0.1346, 51.5292),
  mode = c("WALK", "TRANSIT"), cutoffSec = c(600, 1200, 1800)
)
}
}
\seealso{
Other routing: 
\code{\link{otp_geocode}()},
\code{\link{otp_plan}()},
\code{\link{otp_pointset}()},
\code{\link{otp_routing_options}()},
\code{\link{otp_validate_routing_options}()}
}
\concept{routing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-config.R
\name{otp_validate_config}
\alias{otp_validate_config}
\title{Validate Config Object}
\usage{
otp_validate_config(config, type = attributes(config)$config_type)
}
\arguments{
\item{config}{A named list made/modified from `otp_make_config()`}

\item{type}{type of config file}
}
\description{
Checks if the list of OTP configuration options is valid
}
\details{
Performs basic validity checks on class, max/min values etc as appropriate,
some of more complex parameters are not checked. For more details see:

http://docs.opentripplanner.org/en/latest/Configuration
http://dev.opentripplanner.org/javadoc/1.3.0/org/opentripplanner/routing/core/RoutingRequest.html
}
\examples{
\dontrun{
conf <- otp_make_config("build")
otp_validate_config(conf)
}
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-setup.R
\name{otp_stop}
\alias{otp_stop}
\title{Stop and OTP Instance}
\usage{
otp_stop(warn = TRUE, kill_all = TRUE)
}
\arguments{
\item{warn}{Logical, should you get a warning message}

\item{kill_all}{Logical, should all Java instances be killed?}
}
\value{
This function return a message but no object
}
\description{
OTP is run in Java and requires Java commands to be typed
into the command line. The function allows the parameters
to be defined in R and automatically passed to Java.
This function stops an already running OTP instance
}
\details{
The function assumes you have run otp_setup()
}
\examples{
\dontrun{
otp_stop(kill_all = FALSE)
}
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-surface.R
\name{otp_surface_isochrone}
\alias{otp_surface_isochrone}
\title{Make an isochrone from a surface}
\usage{
otp_surface_isochrone(otpcon = NULL, surface = NULL)
}
\arguments{
\item{otpcon}{OTP connection object produced by otp_connect()}

\item{surface}{A suface list from otp_make_surface()}
}
\value{
Returns a data.frame of travel times
}
\description{
Geo-coding converts a named place, such as a street name into a lng/lat pair.
}
\details{
THis function requires the analysis and pointset features to be
enabled during `otp_setup()`. Thus it will only work with OTP 1.x. For more
detail see the analyst vignette.
}
\examples{
\dontrun{
times <- otp_surface(otpcon, c(-1.17502, 50.64590), "lsoa", path_data)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-connect.R
\name{otp_check_version}
\alias{otp_check_version}
\title{Check the what version of OTP the server is running}
\usage{
otp_check_version(otpcon, warn = TRUE)
}
\arguments{
\item{otpcon}{otpcon object from otp_connect()}

\item{warn}{logical, if TRUE will check that OTP version matches contents of otpcon}
}
\description{
Check the what version of OTP the server is running
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_java}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otp-setup.R
\name{otp_check_java}
\alias{otp_check_java}
\title{Check Java version}
\usage{
otp_check_java(otp_version = 1.5)
}
\arguments{
\item{otp_version}{numeric, OTP version number default 1.5}
}
\description{
Check if you have the correct version of Java for running OTP locally
}
\seealso{
Other setup: 
\code{\link{otp_build_graph}()},
\code{\link{otp_check_version}()},
\code{\link{otp_dl_demo}()},
\code{\link{otp_dl_jar}()},
\code{\link{otp_make_config}()},
\code{\link{otp_setup}()},
\code{\link{otp_stop}()},
\code{\link{otp_validate_config}()},
\code{\link{otp_write_config}()}
}
\concept{setup}
