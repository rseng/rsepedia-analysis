# rdefra: Interact with the UK AIR Pollution Database from DEFRA

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.593187.svg)](https://doi.org/10.5281/zenodo.593187)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.00051/status.svg)](https://doi.org/10.21105/joss.00051)

[![R-CMD-check](https://github.com/ropensci/rdefra/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rdefra/actions)
[![codecov.io](https://codecov.io/gh/ropensci/rdefra/coverage.svg?branch=master)](https://codecov.io/gh/ropensci/rdefra?branch=master)

[![CRAN Status
Badge](http://www.r-pkg.org/badges/version/rdefra)](https://cran.r-project.org/package=rdefra)
[![CRAN Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/rdefra)](https://cran.r-project.org/package=rdefra)
[![CRAN Monthly
Downloads](http://cranlogs.r-pkg.org/badges/rdefra)](https://cran.r-project.org/package=rdefra)
[![](https://badges.ropensci.org/68_status.svg)](https://github.com/ropensci/software-review/issues/68)

The package [rdefra](https://cran.r-project.org/package=rdefra) allows
to retrieve air pollution data from the Air Information Resource
[UK-AIR](https://uk-air.defra.gov.uk/) of the Department for
Environment, Food and Rural Affairs in the United Kingdom. UK-AIR does
not provide a public API for programmatic access to data, therefore this
package scrapes the HTML pages to get relevant information.

This package follows a logic similar to other packages such as
[waterData](https://cran.r-project.org/package=waterData) and
[rnrfa](https://cran.r-project.org/package=rnrfa): sites are first
identified through a catalogue, data are imported via the station
identification number, then data are visualised and/or used in analyses.
The metadata related to the monitoring stations are accessible through
the function `ukair_catalogue()`, missing stations’ coordinates can be
obtained using the function `ukair_get_coordinates()`, and time series
data related to different pollutants can be obtained using the function
`ukair_get_hourly_data()`.

DEFRA’s servers can handle multiple data requests, therefore concurrent
calls can be sent simultaneously using the
[parallel](https://www.R-project.org/) package. Although the limit rate
depends on the maximum number of concurrent calls, traffic and available
infrustracture, data retrieval is very efficient. Multiple years of data
for hundreds of sites can be downloaded in only few minutes.

For similar functionalities see also the
[openair](https://cran.r-project.org/package=openair) package, which
relies on a local copy of the data on servers at King’s College (UK),
and the [ropenaq](https://CRAN.R-project.org/package=ropenaq) which
provides UK-AIR latest measured levels (see
<https://uk-air.defra.gov.uk/latest/currentlevels>) as well as data from
other countries.

## Installation

Get the released version from CRAN:

``` r
install.packages("rdefra")
```

Or the development version from GitHub using the package `remotes`:

``` r
install.packages("remotes")
remotes::install_github("ropensci/rdefra")
```

Load the rdefra package:

``` r
library("rdefra")
```

## Functions

The package logic assumes that users access the UK-AIR database in the
fllowing steps:

1.  Browse the catalogue of available stations and selects some stations
    of interest (see function `ukair_catalogue()`).
2.  Get missing coordinates (see function `ukair_get_coordinates()`).
3.  Retrieves data for the selected stations (see functions
    `ukair_get_site_id()` and `ukair_get_hourly_data()`).

For an in-depth description of the various functionalities andexample
applications, please refer to the package
[vignette](https://github.com/ropensci/rdefra/blob/master/vignettes/rdefra_vignette.Rmd).

## Meta

  - This package and functions herein are part of an experimental open-source project. They are provided as is, without any guarantee.
  - Please [report any issues or
    bugs](https://github.com/ropensci/rdefra/issues).
  - License: [GPL-3](https://opensource.org/licenses/GPL-3.0)
  - This package was reviewed by [Maëlle
    Salmon](https://github.com/maelle) and [Hao
    Zhu](https://github.com/haozhu233) for submission to ROpenSci (see
    review [here](https://github.com/ropensci/software-review/issues/68)) and
    the Journal of Open Source Software (see review
    [here](https://github.com/openjournals/joss-reviews/issues/51)).
  - Cite `rdefra`: `citation(package = "rdefra")`

<br/>

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
rdefra 0.3.8
==============

This release corresponds to the latest CRAN submission.

This is a resubmission due to bug fixing.

## BUG FIXES
* Updated tests to be compatible with PROJ6 [#9](https://github.com/ropensci/rdefra/issues/9)

## MINOR CHANGES
* Removed obsolete packages in 'Suggests'
* Fixed invalid URIs
* The following directory looks like a leftover from knitr

rdefra 0.3.6
==============

## BUG FIXES
* SiteID = NA causes hanging errors [#6](https://github.com/ropensci/rdefra/issues/6)
* Different variables for ukair_get_coordinates() when inputs are fed in differently [#7](https://github.com/ropensci/rdefra/issues/7)

## MINOR CHANGES
* function ukair_get_coords back to original name ukair_get_coordinates

rdefra 0.3.5
==============

## MINOR FIXES
Functions are updated due to a recent modification of the catalogue API.

## MINOR IMPROVEMENTS
Changes made after scanning the package using goodpractice:
* lines no longer than 80 characters
* 84% code coverage
* function names shorter than 30 characters
  - function ukair_get_coordinates now renamed ukair_get_coords

rdefra 0.3.4
==============

## MINOR IMPROVEMENTS
Changes to DESCRIPTION file:

* Added rmarkdown and knitr in Suggests
* Added entry VignetteBuilder: knitr
* Changed all links to the new ropenscilabs github account (ropenscilabs instead of kehraProject)

This repo is now transferred to the ropenscilabs github account.

rdefra 0.3.3
==============

In this release the package was moved to the root directory (needed based on the ropensci review) and the related adjustments made.

## MINOR FIXES

* Corrected units in README (#2)
* Merged README_files folder to assets? (#3)

rdefra 0.3.2
==============

Accepted for pubblication on JOSS

rdefra 0.3.1
==============

Minor changes

rdefra 0.3.0
==============

Minor fixes

rdefra 0.2.0
==============

* Added unist tests (using testthat framework), 
* Added Travis CI integration
* Added a vignette
* Added paper for submission to JOSS
* Added a document to review the package under the ropensci project

rdefra 0.1.0
==============

* Initial release
This is a resubmission due to bug fixing.

---------------------------------

## Release Summary

## BUG FIXES
* Code is now updated to ensure that in case of failure, the code fails gracefully and returns an informative message.

## NOTES
Please note that:
* The following words are not misspelled: DEFRA, Vitolo, al, et, rdefra

## Test environment
* Ubuntu 18.04, R 4.0.4

## R CMD check results

There were no other ERRORs, WARNINGs or NOTEs.
---
title: 'rdefra: Interact with the UK AIR Pollution Database from DEFRA'
bibliography: paper.bib
date: "3 August 2016"
output: pdf_document
tags:
- open data
- air pollution
- R
authors:
- affiliation: Brunel University London
  name: Claudia Vitolo
  orcid: 0000-0002-4252-1176
- affiliation: Brunel University London
  name: Andrew Russell
  orcid: 0000-0001-7120-8499
- affiliation: Brunel University London
  name: Allan Tucker
  orcid: 0000-0001-5105-3506
---

# Summary

Rdefra [@rdefra-archive] is an R package [@R-base] to retrieve air pollution data from the Air Information Resource (UK-AIR) of the Department for Environment, Food and Rural Affairs in the United Kingdom. UK-AIR does not provide a public API for programmatic access to data, therefore this package scrapes the HTML pages to get relevant information.

This package follows a logic similar to other packages such as waterData[@waterdata] and rnrfa[@rnrfa]: sites are first identified through a catalogue, data are imported via the station identification number, then data are visualised and/or used in analyses. The metadata related to the monitoring stations are accessible through the function `ukair_catalogue()`, missing stations' coordinates can be obtained using the function `ukair_get_coordinates()`, and time series data related to different pollutants can be obtained using the function `ukair_get_hourly_data()`.

The package is designed to collect data efficiently. It allows to download multiple years of data for a single station with one line of code and, if used with the parallel package [@R-base], allows the acquisition of data from hundreds of sites in only few minutes.

The figure below shows the 6566 stations with valid coordinates within the UK-AIR (blue circles) database, for 225 of them hourly data is available and their location is shown as red circles.

![UK-AIR monitoring stations (August 2016)](MonitoringStations.png)

# References
