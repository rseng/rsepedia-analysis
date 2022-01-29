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

# rnpn

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran
checks](https://cranchecks.info/badges/worst/rnpn)](https://cranchecks.info/pkgs/rnpn)
[![codecov.io](https://codecov.io/github/ropensci/rnpn/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rnpn?branch=master)
[![R build
status](https://github.com/usa-npn/rnpn//workflows/R-CMD-check/badge.svg)](https://github.com/usa-npn/rnpn//actions)

`rnpn` is an R client for interacting with the USA National Phenology
Network data web services. These services include access to a rich set
of observer-contributed, point-based phenology records as well as
geospatial data products including gridded phenological model and
climatological data.

Documentation is available for the National Phenology Network [API
documentation](https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit?hl=en_US),
which describes the full set of REST services this package wraps.

There is no need for an API key to grab data from the National Phenology
Network but users are required to self identify, on an honor system,
against requests that may draw upon larger datasets. Simply populate the
request\_source parameter, as necessary, with your name or the name of
your institution.

Currently there are services for writing to the database but those
endpoints do require user authentication and are not accessible through
this R wrapper. Please contact the package authors for more information
if that’s what you’re trying to do.

## Installation

This package has evolved slowly and is currently managed in a few
locations, with varying degrees of available functionality.

The original, v. 0.1, iteration of the package is available through CRAN
but has limited functionality and a number of endpoints that have been
deprecated.

CRAN version

``` r
install.packages("rnpn")
```

There’s a newer iteration of the package that includes a lot more
functionality, including the ability to access geospatial data, and is
up-to-date with the backing data services. This version of the package
is actively maintained, but is not managed through CRAN yet and as such
is more bug-prone. This version of the package must be installed through
devtools.

Development version

``` r
install.packages("devtools")
library('devtools')
devtools::install_github("usa-npn/rnpn")
```

``` r
library('rnpn')
```

## The Basics

Many of the functions to search for data require knowing the internal
unique identifiers of some of the database entities to filter the data
down efficiently. For example, if you want to search by species, then
you must know the internal identifier of the species. To get a list of
all available species use the following:

``` r
species_list <- npn_species()
```

Similarly, for phenophases:

``` r
phenophases <- npn_phenophases()
```

### Getting Observational Data

There are four main functions for accessing observational data, at
various levels of aggregation. At the most basic level you can download
the raw status and intensity data.

``` r
some_data <- npn_download_status_data(request_source='Your Name or Org Here',years=c(2015),species_id=c(35),states=c('AZ','IL'))
```

Note that through this API, data can only be filtered chronologically by
full calendar years. You can specify any number of years in each API
call. Also note that request\_source is a required parameter and should
be populated with your name or the name of the organization you
represent. All other parameters are optional but it is highly
recommended that you filter your data search further.

### Getting Geospatial Data

This package wraps around standard WCS endpoints to facilitate the
transfer of raster data. Generally, this package does not focus on
interacting with WMS services, although they are available. To get a
list of all available data layers, use the following:

``` r
layers <- npn_get_layer_details()
```

You can then use the name of the layers to select and download
geospatial data as a raster.

``` r
npn_download_geospatial(coverage_id = 'si-x:lilac_leaf_ncep_historic',date='2016-12-31',format='geotiff',output_path='./six-test-raster.tiff')
```

If you’re looking for a grid value at a specific latitude/longitude,
that is also possible.

``` r
point_value <- npn_get_point_data('si-x:lilac_leaf_ncep_historic',date='2016-12-31',lat=38.5,long=-110.7)
```

## What’s Next

Please read and review the vignettes for this package to get further
information about the full scope of functionality available.

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/rnpn/issues).
  - License: MIT
  - Get citation information for `rnpn` in R doing `citation(package =
    'rnpn')`
  - Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.

[![image](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
rnpn 1.1.1 (2020-10-27)
=======================

### NEW FEATURES

* Total overhaul of the rNPN package
* Added functions for directly downloading different observation record data types
* Added additional utility and lookup type functions
* Added functions for downloading USA-NPN raster data and geospatial values by latitutde/longitude
* Deprecated the following functions: lookup_names, npn_allobssp, npn_indsatstations, npn_indspatstations, npn_species_comm, npn_species_itis, npn_species_sci, npn_stationsbystate, npn_stationswithspp



rnpn 0.1.0
==========

### NEW FEATURES

* released to CRAN
I have read and agree to the the CRAN policies at
http://cran.r-project.org/web/packages/policies.html


## Test environments

R CMD CHECK passed on local Windows 10 and Ubuntu 18 using R 4.0.2 
Also passed checks on Github Actions on macOS, Windows with R 3.6
and ubuntu 16 using R 3.5.

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

Changing maintainers. This project is transitioning to a new maintainer.

## Downstream dependencies

There were no downstream dependencies.


## Other Notes

This version is a total overhaul of the previous version of the package currently available on CRAN.

Thanks! 
