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
