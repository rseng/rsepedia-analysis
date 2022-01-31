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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

![](https://github.com/ropensci/dbhydroR/raw/master/inst/images/profile.png)

# Programmatic access to the South Florida Water Management District's [DBHYDRO database](https://www.sfwmd.gov/science-data/dbhydro)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/dbhydroR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/dbhydroR/actions)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/dbhydroR)](https://cran.r-project.org/package=dbhydroR) 
[![](https://badges.ropensci.org/61_status.svg)](https://github.com/ropensci/software-review/issues/61)
[![DOI](https://zenodo.org/badge/64503356.svg)](https://zenodo.org/badge/latestdoi/64503356)

`dbhydroR` provides scripted access to the South Florida Water Management District's DBHYDRO database which holds over 35 million hydrologic and water quality records from the Florida Everglades and surrounding areas. 

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

Station IDs and date ranges can be viewed in the [Environmental Monitoring Location Maps](https://www.sfwmd.gov/documents-by-tag/emmaps). Test names can be viewed in the Data Types Metadata Table on the DBHYDRO website.

#### One variable at one station
```
get_wq(station_id = "FLAB08", date_min = "2011-03-01", 
      date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")
```

#### One variable at multiple stations 
```
get_wq(station_id = c("FLAB08","FLAB09"), date_min = "2011-03-01",
      date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")
```

#### One variable at a wildcard station
```
get_wq(station_id = c("FLAB0%"), date_min = "2011-03-01", 
      date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")
```

#### Multiple variables at multiple stations
```
get_wq(station_id = c("FLAB08","FLAB09"), date_min = "2011-03-01",
      date_max = "2012-05-01", test_name = c("CHLOROPHYLL-A, SALINE",
      "SALINITY"))
```

#### Operate on raw data
```
raw_data <- get_wq(station_id = "FLAB08", date_min = "2011-03-01", 
      date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE", raw = TRUE)

clean_wq(raw_data)
```

### Hydrologic data
Station IDs and date ranges can be viewed in the [Environmental Monitoring Location Maps](https://www.sfwmd.gov/documents-by-tag/emmaps). 

#### Identify unique time series (dbkeys) before-hand
```
get_dbkey(stationid = "C111%", stat = 'MEAN', category = "WQ", detail.level = "full")
get_hydro(dbkey = 38104, date_min = "2009-01-01", date_max = "2009-01-12")
```

#### Pass station info on-the-fly
```
get_hydro(date_min = "2013-01-01", date_max = "2013-02-02",
         stationid = "JBTS", category = "WEATHER", param = "WNDS",
         freq = "DA", stat = "MEAN", recorder = "CR10", agency = "WMD")
```

#### Operate on raw data
```
raw_data <- get_hydro(date_min = "2013-01-01", date_max = "2013-02-02",
         stationid = "JBTS", category = "WEATHER", param = "WNDS",
         freq = "DA", stat = "MEAN", recorder = "CR10", agency = "WMD", raw = TRUE)
         
clean_hydro(raw_data)
```

## References

`vignette("dbhydroR", package = "dbhydroR")`

[DBHYDRO User's Guide](https://www.sfwmd.gov/sites/default/files/documents/dbhydrobrowseruserdocumentation.pdf)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/dbhydroR/issues).

* Get citation information for `dbhydroR` in R by running `citation(package = 'dbhydroR')`

* Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/dbhydroR/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "dbhydroR: An R package to access the DBHYDRO Environmental Database"
author: "Jemma Stachelek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: bib.bib
vignette: >
  %\VignetteIndexEntry{An R package to access the DBHYDRO Environmental Database}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

This document introduces the `dbhydroR` package and its associated functions. These functions are aimed at improving programmatic workflows that query the DBHYDRO Environmental Database which holds over 35 million hydrologic and water quality records from the Florida Everglades and surrounding areas.  

## Package installation

Computers running the Windows operating system can only install binary package archive files unless they have additional [compiler software](https://cran.r-project.org/bin/windows/Rtools/) installed. Without this software, `dbhydroR` can be installed from CRAN by running the following command in the `R` console:

```{r echo = FALSE}
options(width = 50, useFancyQuotes = FALSE, prompt = " ", continue = " ")
```

### Stable version from CRAN

```{r eval = FALSE}
install.packages("dbhydroR")
```

Otherwise, the `dbhydroR` can be installed by running the following command in the `R` console:

### or development version from Github

```{r eval = FALSE}
devtools::install_github("ropensci/dbhydroR")
```

\vspace{8pt}
\noindent Once installed, the package can be loaded using the following command:


```{r message=FALSE}
library(dbhydroR)
```


## Composing database queries
### Water quality data

Water quality data can be retrieved using the `get_wq` function which takes four required arguments. The user must specify a station ID, a test name, and a date range. Station IDs can be located on the [SFWMD Station Maps](https://www.sfwmd.gov/documents-by-tag/emmaps). An abbreviated list of available test names can be found in the [appendix](#appendix) to this document while a full listing can be found on the DBHYDRO metadata page. Dates must be specified in YYYY-MM-DD format (e.g. 2015-02-26).   The following set of examples retrieve measurements between March 2011 and May 2012. They can be run from the R console by issuing the command:

```{r eval = FALSE}
example(get_wq)
```

 * One variable at one station

```{r eval = FALSE}
get_wq(station_id = "FLAB08", date_min = "2011-03-01", 
      date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")
```

 * One variable at multiple stations

```{r eval = FALSE}
get_wq(station_id = c("FLAB08","FLAB09"), date_min = "2011-03-01",
      date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")
```

 * One variable at a wildcard station

```{r eval = FALSE}
get_wq(station_id = c("FLAB0%"), date_min = "2011-03-01", 
      date_max = "2012-05-01", test_name = "CHLOROPHYLL-A, SALINE")
```

 * Multiple variables at multiple stations

```{r eval = FALSE}
get_wq(station_id = c("FLAB08","FLAB09"), date_min = "2011-03-01",
      date_max = "2012-05-01", test_name = c("CHLOROPHYLL-A, SALINE",
      "SALINITY"))
```

\noindent By default, `get_wq` returns a *cleaned output*. First, the cleaning function \verb|clean_wq| converts the raw output from native DBHYDRO *long* format (each piece of data on its own row) to *wide* format (each site x variable combination in its own column) using the reshape2 package [@reshape2]. Next, the extra columns associated with QA flags, LIMS, and District receiving are removed. Finally, row entries associated with QA field *blanks*, which are used to check on potential sources of contamination, are removed.  Setting the `raw` flag to `TRUE` will force \verb|get_wq| to retain information on QA field blanks as well as the other QA fields. An example query that retains this information and the original *long* formatting is shown below.

```{r eval = FALSE}
raw_wq <- get_wq(station_id = "FLAB08", date_min = "2011-03-01", 
      date_max = "2011-05-01", test_name = "CHLOROPHYLL-A, SALINE",
      raw = TRUE)
```

This raw data can then be cleaned using the \verb|clean_wq| function:

```{r eval = FALSE}
clean_wq(raw_wq)
```

### Hydrologic data

Hydrologic time series data can be retrieved using the `get_hydro` function. The first task to accomplish prior to running `get_hydro` is to identify one or more dbkeys which correspond to unique site x variable time-series. This can be done before-hand using the `get_dbkey` function, the [SFWMD Station Maps](https://www.sfwmd.gov/documents-by-tag/emmaps) or the DBHYDRO Browser. One useful strategy for finding desired dbkeys is to run the `get_dbkey`   function interactively using progressively narrower search terms. For example, suppose we are interested in daily average wind data at Joe Bay but we have no alphanumeric `dbkey`. Initially we could run `get_dbkey` with the `detail.level` set to "summary".

```{r eval = FALSE}
get_dbkey(stationid = "JBTS", category = "WEATHER", param = "WNDS",
         detail.level = "summary")
```

\noindent Our search returns two results but only one of them has a daily average (DA) measurement frequency. We can verify the remaining attributes of our likely dbkey by setting the `freq` parameter to "DA" and the `detail.level` parameter to "full".

```{r eval = FALSE}
get_dbkey(stationid = "JBTS", category = "WEATHER", param = "WNDS",
         freq = "DA", detail.level = "full")
```

\noindent This exact dbkey can only be returned reliably by specifying all of the `get_dbkey` parameters applicable to the "WEATHER" category.

```{r eval = FALSE}
get_dbkey(stationid = "JBTS", category = "WEATHER", param = "WNDS",
         freq = "DA", stat = "MEAN", recorder = "CR10", agency = "WMD",
         detail.level = "dbkey")
```

\noindent Now that we have our dbkey in hand, we can use is as input to `get_hydro`. In addition to a dbkey, we must specify a date range. Dates must be entered in YYYY-MM-DD format (e.g. 2015-02-26).

```{r eval = FALSE}
get_hydro(dbkey = "15081",
         date_min = "2013-01-01", date_max = "2013-02-02")
```

\noindent Alternatively, we can specify a set of arguments in our call to `get_hydro` that will be passed to `get_dbkey` on-the-fly. Use caution when using this strategy as complex stationid/category/parameter combinations can easily cause errors or return unexpected results. It is good practice to pre-screen your parameter values using `get_dbkey`.

```{r eval = FALSE}
get_hydro(date_min = "2013-01-01", date_max = "2013-02-02",
         stationid = "JBTS", category = "WEATHER", param = "WNDS",
         freq = "DA", stat = "MEAN", recorder = "CR10", agency = "WMD")
```

\noindent The contents of multiple data streams can be returned by specifying multiple dbkeys or entering on-the-fly `get_dbkey` queries that return multiple dbkeys.

```{r eval = FALSE}
get_hydro(dbkey = c("15081", "15069"), date_min = "2013-01-01",
         date_max = "2013-02-02")
```

```{r eval = FALSE}
get_hydro(date_min = "2013-01-01", date_max = "2013-02-02",
         category = "WEATHER", stationid = c("JBTS", "MBTS"),
         param = "WNDS", freq = "DA", stat = "MEAN")
```

\noindent More `get_hydro` examples including queries of other `category` values ("SW", "GW", and "WQ") can be viewed by issuing the following commands from the `R` console:

```{r eval = FALSE}
example(get_dbkey)
example(get_hydro)
```

\noindent By default, `get_hydro` returns a *cleaned output*. First, the cleaning function `clean_hydro` converts the raw output from native DBHYDRO *long* format (each piece of data on its own row) to *wide* format (each site x variable combination in its own column) using the reshape2 package [@reshape2]. Next, some extra columns are removed that are associated with measurement location (longitude/latitude), frequency, and QA flags are removed. Setting the `raw` flag to `TRUE` will force `get_hydro` to retain the original formatting and metadata fields. An example query that retains this information and the original *long* formatting is shown below.

```{r eval = FALSE}
raw_data <- get_hydro(date_min = "2013-01-01", date_max = "2013-02-02",
         stationid = "JBTS", category = "WEATHER", param = "WNDS",
         freq = "DA", stat = "MEAN", recorder = "CR10", agency = "WMD", raw = TRUE)
         
clean_hydro(raw_data)
```

## Appendix \label{sec:appendix}  
### Test names
There are many test names available in DBHYDRO. A subset of these are detailed in the following table.

| Code |
| --- |
| AMMONIA-N |
| CARBON, TOTAL ORGANIC |
| CHLOROPHYLL-A(LC) |
| CHLOROPHYLL-B(LC) |
| CHLOROPHYLL-A, SALINE |
| DISSOLVED OXYGEN |
| KJELDAHL NITROGEN,TOTAL |
| NITRATE+NITRITE-N |
| NITRITE-N |
| PHEOPHYTIN-A(LC) |
| PHOSPHATE,ORTHO AS P |
| PHOSPHATE,TOTAL AS P |
| SALINITY |
| SILICA |
| SP CONDUCTIVITY, FIELD |
| TEMP |
| TOTAL NITROGEN |
| TURBIDITY |

### Further reading
See section on URL-based data access in the [DBHYDRO Browser User's Guide](https://www.sfwmd.gov/sites/default/files/documents/dbhydrobrowseruserdocumentation.pdf)

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbhydro_get.R
\name{get_wq}
\alias{get_wq}
\alias{getwq}
\title{Retrieve water quality data from the DBHYDRO Environmental Database}
\usage{
get_wq(
  station_id = NA,
  date_min = NA,
  date_max = NA,
  test_name = NA,
  mdl_handling = "raw",
  raw = FALSE,
  qc_strip = "N",
  qc_field = "N",
  test_number = NA,
  v_target_code = "file_csv",
  sample_id = NA,
  project_code = NA
)
}
\arguments{
\item{station_id}{character string of station id(s). See the SFWMD station
search utility for specific options}

\item{date_min}{character date must be in POSIXct YYYY-MM-DD format}

\item{date_max}{character date must be in POSIXct YYYY-MM-DD format}

\item{test_name}{character string of test name(s). See the SFWMD
Station Maps at \url{https://www.sfwmd.gov/documents-by-tag/emmaps}
for specific options}

\item{mdl_handling}{character string specifying the handling of measurement
values below the minimum detection limit (MDL). Example choices for this
argument include:
\itemize{
\item \code{raw}: Returns values exactly as they are stored in the database.
Current practice is to return values below the MDL as 0 minus the
uncertainty estimate.
\item \code{half}: Returns values below the MDL as half the MDL
\item \code{full}: Returns values below the MDL as the MDL
}}

\item{raw}{logical default is FALSE, set to TRUE to return data in "long"
format with all comments, qa information, and database codes included}

\item{qc_strip}{logical set TRUE to avoid returning QAQC flagged data entries}

\item{qc_field}{logical set TRUE to avoid returning field QC results}

\item{test_number}{numeric test name alternative (not implemented)}

\item{v_target_code}{string print to file? (not implemented)}

\item{sample_id}{numeric (not implemented)}

\item{project_code}{numeric (not implemented)}
}
\description{
Retrieve water quality data from the
DBHYDRO Environmental Database
}
\details{
By default, \code{get_wq} returns a cleaned output. First, the
cleaning function \code{\link{clean_wq}} converts the raw output from native
DBHYDRO long format (each piece of data on its own row) to wide format (each
site x variable combination in its own column). Next, the extra columns
associated with QA flags, LIMS, and District receiving are removed. Finally,
row entries associated with QA field blanks, which are used to check on
potential sources of contamination, are removed. Setting the raw flag to TRUE
will force getwq to retain information on QA field blanks as well as the
other QA fields.
}
\examples{

\dontrun{
#one variable and one station
get_wq(station_id = "FLAB08",
date_min = "2011-03-01", date_max = "2012-05-01",
test_name = "CHLOROPHYLLA-SALINE")

#one variable at multiple stations
get_wq(station_id = c("FLAB08", "FLAB09"),
date_min = "2011-03-01", date_max = "2012-05-01",
test_name = "CHLOROPHYLLA-SALINE")

#One variable at a wildcard station
get_wq(station_id = c("FLAB0\%"),
date_min = "2011-03-01",
date_max = "2012-05-01",
test_name = "CHLOROPHYLLA-SALINE")

#multiple variables at multiple stations
get_wq(station_id = c("FLAB08", "FLAB09"),
date_min = "2011-03-01", date_max = "2012-05-01",
test_name = c("CHLOROPHYLLA-SALINE", "SALINITY"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbhydro_clean.R
\name{clean_wq}
\alias{clean_wq}
\alias{cleanwq}
\title{Clean raw water quality DBHYDRO data retrievals}
\usage{
clean_wq(dt, raw = FALSE, mdl_handling = "raw")
}
\arguments{
\item{dt}{data.frame output of \code{\link{getwq}}}

\item{raw}{logical default is FALSE, set to TRUE to return data in "long"
format with all comments, qa information, and database codes included}

\item{mdl_handling}{character string specifying the handling of measurement
values below the minimum detection limit (MDL). Example choices for this
argument include:
\itemize{
\item \code{raw}: Returns values exactly as they are stored in the database.
Current practice is to return values below the MDL as 0 minus the uncertainty
estimate.
\item \code{half}: Returns values below the MDL as half the MDL
\item \code{full}: Returns values below the MDL as the MDL
}}
}
\description{
Removes extra columns associated with QA flags and QA blanks
which are used to check on potential sources of contamination. If raw is set
to TRUE, \code{\link{get_wq}} results are converted from long (each piece of
data on its own row) to \code{wide} format (each site x variable combination
in its own column).
}
\examples{
\dontrun{
#check handling of values below MDL
dt <- getwq("FLAB01", "2014-09-14", "2014-09-18", "NITRATE+NITRITE-N",
 raw = TRUE)
clean_wq(dt, mdl_handling = "raw")
clean_wq(dt, mdl_handling = "half")
}

dt <- read.csv(system.file("extdata", "testwq.csv", package = "dbhydroR"))
clean_wq(dt)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbhydoR-package.R
\docType{package}
\name{dbhydroR-package}
\alias{dbhydroR-package}
\alias{dbhydroR}
\title{dbhydroR}
\description{
dbhydroR is an R interface to the South Florida Water Management District's DBHYDRO database which holds over 35 million hydrologic and water quality records from the Florida Everglades and surrounding areas.
}
\section{Hydrologic Data}{

The \code{\link[dbhydroR]{get_hydro}} function provides the capability to return:
\itemize{
  \item weather data
  \item surfacewater data
  \item groundwater data
  \item water-quality sonde data
}
}

\section{Water Quality Data}{

The \code{\link[dbhydroR]{get_wq}} function provides the capability to return:
\itemize{
  \item water quality data
}
}

\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbhydro_clean.R
\name{clean_hydro}
\alias{clean_hydro}
\alias{cleanhydro}
\title{Clean raw hydrologic DBHYDRO data retrievals}
\usage{
clean_hydro(dt)
}
\arguments{
\item{dt}{data.frame output of \code{\link[dbhydroR]{gethydro}}}
}
\description{
Converts output of \code{\link{get_hydro}} from long (each piece
of data on its own row) to wide format (each site x variable combination in
its own column). Metadata (station-name, variable, measurement units) is
parsed so that it is wholly contained in column names.
}
\examples{
\dontrun{
clean_hydro(gethydro(dbkey = "15081", date_min = "2013-01-01", date_max = "2013-02-02", raw = TRUE))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbhydro_get.R
\name{get_dbkey}
\alias{get_dbkey}
\alias{getdbkey}
\title{Query dbkey information}
\usage{
get_dbkey(
  category,
  stationid = NA,
  param = NA,
  freq = NA,
  longest = FALSE,
  stat = NA,
  recorder = NA,
  agency = NA,
  strata = NA,
  detail.level = "summary",
  ...
)
}
\arguments{
\item{category}{character string, choice of "WEATHER", "SW", "GW", or "WQ"}

\item{stationid}{character string specifying station name}

\item{param}{character string specifying desired parameter name}

\item{freq}{character string specifying collection frequency (daily = "DA")}

\item{longest}{logical limit results to the longest period-of-record?}

\item{stat}{character string specifying statistic type}

\item{recorder}{character string specifying recorder information}

\item{agency}{character string specifying collector agency}

\item{strata}{numeric vector of length 2 specifying a range of z-coordinates
relative to local ground elevation. Only applicable for queries in the
"WEATHER" and "GW" categories.}

\item{detail.level}{character string specifying the level of detail to return.
Choices are "full", "summary", and "dbkey".}

\item{...}{Options passed as named parameters}
}
\description{
Retrieve a data.frame summary including dbkeys or a vector of
dbkeys corresponding to specified parameters
}
\details{
A \code{dbkey} represents a unique station x variable time-series. A
value in the "Recorder" field of "PREF" should be used whenever possible.
This indicates that the dataset has been checked by the SFWMD modelling
group.
}
\examples{
\dontrun{
# Weather
get_dbkey(stationid = "JBTS", category = "WEATHER", param = "WNDS",
detail.level = "summary")
get_dbkey(stationid = "JBTS", category = "WEATHER", param = "WNDS",
detail.level = "dbkey")

# query on multiple values
get_dbkey(stationid = c("MBTS", "JBTS"), category = "WEATHER",
param = "WNDS", freq = "DA", detail.level = "dbkey")


# Surfacewater
get_dbkey(stationid = "C111\%", category = "SW")
get_dbkey(category = "SW", stationid = "LAKE\%", detail.level = "full")

# Groundwater
get_dbkey(stationid = "C111\%", category = "GW")
get_dbkey(stationid = "C111AE", category = "GW", param = "WELL",
freq = "DA", stat = "MEAN", strata = c(9, 22), recorder = "TROL",
 agency = "WMD", detail.level = "full")

# Water Quality
get_dbkey(stationid = "C111\%", category = "WQ")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbhydro_get.R
\name{get_hydro}
\alias{get_hydro}
\alias{gethydro}
\title{Retrieve hydrologic data from the DBHYDRO Environmental Database}
\usage{
get_hydro(dbkey = NA, date_min = NA, date_max = NA, raw = FALSE, ...)
}
\arguments{
\item{dbkey}{character string specifying a unique data series.
See \code{\link[dbhydroR]{get_dbkey}}}

\item{date_min}{character date must be in YYYY-MM-DD format}

\item{date_max}{character date must be in YYYY-MM-DD format}

\item{raw}{logical default is FALSE, set to TRUE to return data in "long"
format with all comments, qa information, and database codes included.}

\item{...}{Options passed on to \code{\link[dbhydroR]{get_dbkey}}}
}
\description{
Retrieve hydrologic data from the DBHYDRO Environmental Database
}
\details{
\code{get_hydro} can be run in one of two ways.

\itemize{

\item The first, is to identify one or more \code{dbkeys} before-hand that
correspond to unique data series and are passed to the \code{dbkey}
argument. \code{dbkeys} can be found by:
\itemize{ \item iterative calls to \code{\link{get_dbkey}} (see example)
\item using the Environmental Monitoring Location Maps
(\url{https://www.sfwmd.gov/documents-by-tag/emmaps})
\item using the DBHYDRO Browser.
}

\item The second way to run \code{get_hydro} is to specify additional
arguments to \code{...} which are passed to \code{\link{get_dbkey}}
on-the-fly.

}
By default, \code{get_hydro} returns a cleaned output where metadata
(station-name, variable, measurement units) is wholly contained in the column
name. This is accomplished internally by the \code{\link{clean_hydro}}
function. If additional metadata such as latitude and longitude are desired
set the \code{raw} argument to \code{TRUE}.
}
\examples{
\dontrun{
#One variable/station time series
get_hydro(dbkey = "15081", date_min = "2013-01-01", date_max = "2013-02-02")

#Multiple variable/station time series
get_hydro(dbkey = c("15081", "15069"),
date_min = "2013-01-01", date_max = "2013-02-02")

#Instantaneous hydro retrieval
get_hydro(dbkey = "IY639", date_min = "2015-11-01", date_max = "2015-11-04")

#Looking up unknown dbkeys on the fly
get_hydro(stationid = "JBTS", category = "WEATHER",
param = "WNDS", freq = "DA", date_min = "2013-01-01",
date_max = "2013-02-02")
}
}
