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

<!-- README.md is generated from README.Rmd. Please edit that file -->

![](https://github.com/ropensci/dbhydroR/raw/master/inst/images/profile.png)

# Programmatic access to the South Florida Water Management District’s [DBHYDRO database](https://www.sfwmd.gov/science-data/dbhydro)

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/dbhydroR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/dbhydroR/actions)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/dbhydroR)](https://cran.r-project.org/package=dbhydroR)
[![](https://badges.ropensci.org/61_status.svg)](https://github.com/ropensci/software-review/issues/61)
[![DOI](https://zenodo.org/badge/64503356.svg)](https://zenodo.org/badge/latestdoi/64503356)

`dbhydroR` provides scripted access to the South Florida Water
Management District’s DBHYDRO database which holds over 35 million
hydrologic and water quality records from the Florida Everglades and
surrounding areas.

## Installation

### Stable version from CRAN

`install.packages("dbhydroR")`

### or development version from Github

`install.packages("devtools") # Requires RTools if using Windows`

`devtools::install_github("ropensci/dbhydroR")`

## Usage

### Load dbhydroR

`library("dbhydroR")`

### Water Quality Data

Station IDs and date ranges can be viewed in the [Environmental
Monitoring Location
Maps](https://www.sfwmd.gov/documents-by-tag/emmaps). Test names can be
viewed in the Data Types Metadata Table on the DBHYDRO website.

#### One variable at one station

    get_wq(station_id = "FLAB08", date_min = "2011-03-01", 
          date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")

#### One variable at multiple stations

    get_wq(station_id = c("FLAB08","FLAB09"), date_min = "2011-03-01",
          date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")

#### One variable at a wildcard station

    get_wq(station_id = c("FLAB0%"), date_min = "2011-03-01", 
          date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")

#### Multiple variables at multiple stations

    get_wq(station_id = c("FLAB08","FLAB09"), date_min = "2011-03-01",
          date_max = "2012-05-01", test_name = c("CHLOROPHYLL-A, SALINE",
          "SALINITY"))

#### Operate on raw data

    raw_data <- get_wq(station_id = "FLAB08", date_min = "2011-03-01", 
          date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE", raw = TRUE)

    clean_wq(raw_data)

### Hydrologic data

Station IDs and date ranges can be viewed in the [Environmental
Monitoring Location
Maps](https://www.sfwmd.gov/documents-by-tag/emmaps).

#### Identify unique time series (dbkeys) before-hand

    get_dbkey(stationid = "C111%", stat = 'MEAN', category = "WQ", detail.level = "full")
    get_hydro(dbkey = 38104, date_min = "2009-01-01", date_max = "2009-01-12")

#### Pass station info on-the-fly

    get_hydro(date_min = "2013-01-01", date_max = "2013-02-02",
             stationid = "JBTS", category = "WEATHER", param = "WNDS",
             freq = "DA", stat = "MEAN", recorder = "CR10", agency = "WMD")

#### Operate on raw data

    raw_data <- get_hydro(date_min = "2013-01-01", date_max = "2013-02-02",
             stationid = "JBTS", category = "WEATHER", param = "WNDS",
             freq = "DA", stat = "MEAN", recorder = "CR10", agency = "WMD", raw = TRUE)
             
    clean_hydro(raw_data)

## References

`vignette("dbhydroR", package = "dbhydroR")`

[DBHYDRO User’s
Guide](https://www.sfwmd.gov/sites/default/files/documents/dbhydrobrowseruserdocumentation.pdf)

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/dbhydroR/issues).

-   Get citation information for `dbhydroR` in R by running
    `citation(package = 'dbhydroR')`

-   Please note that this project is released with a [Contributor Code
    of
    Conduct](https://github.com/ropensci/dbhydroR/blob/master/CONDUCT.md).
    By participating in this project you agree to abide by its terms

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# dbhydroR 0.2-9

## Minor changes

* Maintenance release to fix testing on CRAN

# dbhydroR 0.2-8 

## Minor changes

* Switched all URLs from http to https
* Now caching test web requests

# dbhydroR 0.2-7 (2019-02-15)

## Bug fixes

* Fixed critical bug in `get_hydro` causing data parsing failure in all cases (#16)

# dbhydroR 0.2-6 (2018-07-19)

## Bug fixes

* Fixed critical bug in `get_hydro` causing data parsing failure in all cases

# dbhydroR 0.2-5 (2018-05-21)

## Bug fixes

* `get_dbkey` was incorrectly processing data headers

## Minor changes

* Rebranded from ropenscilabs to ropensci
* Converted vignette to rmarkdown

# dbhydroR 0.2-4 (2017-10-30)

## Bug fixes

* The ArcGIS online station map no longer resolves. Links have been updated.
* Sweave sty files are excluded in CRAN build.

# dbhydroR 0.2-3 (2017-08-02)

## Bug fixes

* `get_hydro()` now resolves multiple matching of on-the-fly dbkeys to the one with the longest period of record.

## Minor changes

* Fixed broken links
* Add rOpenSci badge

# dbhydroR 0.2-2 (2017-02-03)

## Bug fixes

`get_hydro()` now works if a `dbkey` contains leading zeros

# dbhydroR 0.2-1 (2016-11-23)

## Minor changes

* Improved installation instructions in vignette.
* Added package level documentation.
* Added rOpenSci branding.
* Use https. #6

# dbhydroR 0.2

## Major changes

* The package API has been changed to underscored function names. `getwq()`, `gethydro()`, and `getdbkey()` are now deprecated in favor of `get_wq()`, `get_hydro()`, `get_dbkey()`.

## Bug fixes

* `getdbkey()` is no longer limited to < 100 results
* MDL (Minumum Detection Limit) handling now occurs in `getwq()` regardless of how the `raw` parameter is set
* `getwq()` returns a no data warning even if the `raw` parameter is set to `TRUE`
* `gethydro()` and `getwq()` date/time stamps are now forced to the `EST` timezone independently of the user environment
* The character encoding of function results is forced to `UTF-8` regardless of the user environment

## Minor changes

* Documentation formatting is now consistent with CRAN policies
* Added links to the ArcGIS Online Station Map in the README and vignette
* `getdbkey()` coordinates are now in decimal degree format

# dbhydroR 0.1-6

## Minor changes

* Added argument to handle MDLs (Minimum Detection Limits) in `getwq()`


# dbhydroR 0.1-5

## Major changes

* Added ability to pass a vector of values to `getdbkey()` arguments
* Added ability to fully define a unique dbkey in `getdbkey()`

## Minor changes

* Document MDL handling in `cleanwq()`
* Added unit tests
* Remove standalone plotting functions
* Cleanup source code formatting

# dbhydroR 0.1-4

## Bug fixes

* Improvements to `gethydro()` to guess missing column names of instantaneous data
## Test environments

* ubuntu 18 (on ghactions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note
