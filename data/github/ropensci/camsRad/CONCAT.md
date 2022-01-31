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
camsRad
=======

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Travis-CI Build Status](https://travis-ci.org/ropensci/camsRad.svg?branch=master)](https://travis-ci.org/ropensci/camsRad) [![codecov.io](https://codecov.io/gh/ropenscilabs/camsRad/coverage.svg?branch=master)](https://codecov.io/gh/ropenscilabs/camsRad) [![](https://badges.ropensci.org/72_status.svg)](https://github.com/ropensci/onboarding/issues/72)

`camsRad` is a R client for [CAMS Radiation Service](http://www.soda-pro.com/web-services/radiation/cams-radiation-service). CAMS Radiation Service provides time series of global, direct, and diffuse irradiations on horizontal surface, and direct irradiation on normal plane for the actual weather conditions as well as for clear-sky conditions. The geographical coverage is the field-of-view of the Meteosat satellite, roughly speaking Europe, Africa, Atlantic Ocean, Middle East (-66° to 66° in both latitudes and longitudes). The time coverage of data is from 2004-02-01 up to 2 days ago. Data are available with a time step ranging from 15 min to 1 month. Target audience are researchers, developers and consultants in need of high resolution solar radiations time series.

Quick start
-----------

### Install

Dev version from GitHub.

``` r
# CRAN version
install.packages("camsRad")

# Or Github version
if (!require('devtools')) install.packages('devtools')
devtools::install_github("ropensci/camsRad")
```

``` r
library("camsRad")
```

### Authentication

To access the CAMS Radiation Service you need to register at <http://www.soda-pro.com/web-services/radiation/cams-radiation-service>. The email you use at the registration step will be used for authentication, and need to be set with `cams_set_user()`.

``` r
# Authentication
cams_set_user("your@email.com") # An email registered at soda-pro.com
```

### Example 1

Get hourly CAMS solar data into a R data frame. For the location 60° latitude and 15° longitude, and for period 2016-01-01 to 2016-01-15.

``` r

df <- cams_get_radiation(
  lat=60, lng=15, 
  date_begin="2016-07-01", 
  date_end="2016-07-01")
print(df)
```

### Example 2

Retrieve daily CAMS solar data in netCDF format. You need to have the `ncdf4` package installed.

``` r
library(ncdf4)

filename <- paste0(tempfile(), ".nc")

r <- cams_api(
  60, 15, "2016-06-01", "2016-06-10", 
  format="application/x-netcdf",
  time_step = "P01D",
  filename=filename)

# Access the on disk stored ncdf4 file 
nc <- nc_open(filename)

# list names of available variables
names(nc$var)

# create data.frame with timestamp and global horizontal irradiation and plot it
df <- data.frame(
  timestamp = as.POSIXct(nc$dim$time$vals, "UTC", origin="1970-01-01"),
  GHI = ncvar_get(nc, "GHI"))

plot(df, type="l")

nc_close(nc)
```

Meta
----

-   This package and functions herein are provided as is, without any guarantee.
-   Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
-   Please [report any issues or bugs](https://github.com/ropensci/camsRad/issues).

<!--[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org) 
doesn´t knit. add following to the .md file
[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
-->
[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# camsRad 0.3.0 (2016-11-07)
=========================
* Compliance with rOpenSci review
### MINOR IMPROVEMENTS
* covr::codecov() badge included
* more testing
* removed dependency to readr package
* return data.frame instead of tibble
* new authentication method
* more examples and improved documentation
* cams_api stops if error in httr calls or errir in returned content

# camsRad 0.2.0 (2016-08-22)
=========================
* Realese for onboarding to rOpenSci
### MINOR IMPROVEMENTS
* Exported functions are documented with examples
* README extended with examples
* vignette added
* Travis-CI added
* Test with testhat added


# camsRad 0.1.0 (2016-08-18)
=========================
* Initial public release
## Resubmission IV

This is a resubmission. In this version I have:

* Added brackets to the url in the Description field of the DESCRIPTION file.


This is the first submission of this packe.
[Review at rOpenSci accepted](https://github.com/ropensci/onboarding/issues/72)

## Test environments
* local Windows 7 x64 install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "inst/img/README-",
  warning = FALSE
)
if(Sys.getenv("CAMS_USERNAME")=="") {
  knitr::opts_chunk$set(eval = FALSE) # doesn't build README if CAMS_USERNAME is not set
}
```

camsRad
=======
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/camsRad.svg?branch=master)](https://travis-ci.org/ropensci/camsRad)
[![codecov.io](https://codecov.io/gh/ropenscilabs/camsRad/coverage.svg?branch=master)](https://codecov.io/gh/ropenscilabs/camsRad)
[![](https://badges.ropensci.org/72_status.svg)](https://github.com/ropensci/onboarding/issues/72)


`camsRad` is a R client for [CAMS Radiation Service](http://www.soda-pro.com/web-services/radiation/cams-radiation-service). CAMS Radiation Service provides time series of global, direct, and diffuse irradiations on horizontal surface, and direct irradiation on normal plane for the actual weather conditions as well as for clear-sky conditions. The geographical coverage is the field-of-view of the Meteosat satellite, roughly speaking Europe, Africa, Atlantic Ocean, Middle East (-66° to 66° in both latitudes and longitudes). The time coverage of data is from 2004-02-01 up to 2 days ago. Data are available with a time step ranging from 15 min to 1 month. Target audience are researchers, developers and consultants in need of high resolution solar radiations time series. 

## Quick start

### Install

Dev version from GitHub.

```{r eval=FALSE}
# CRAN version
install.packages("camsRad")

# Or Github version
if (!require('devtools')) install.packages('devtools')
devtools::install_github("ropensci/camsRad")
```

```{r}
library("camsRad")
```

### Authentication

To access the CAMS Radiation Service you need to register at [http://www.soda-pro.com/web-services/radiation/cams-radiation-service](http://www.soda-pro.com/web-services/radiation/cams-radiation-service). The email you use at the registration step will be used for authentication, and need to be set with  `cams_set_user()`.

```{r eval=FALSE}
# Authentication
cams_set_user("your@email.com") # An email registered at soda-pro.com
```


### Example 1

Get hourly CAMS solar data into a R data frame. For the location 60° latitude and 15° longitude, and for period 2016-01-01 to 2016-01-15.
```{r, eval=FALSE}

df <- cams_get_radiation(
  lat=60, lng=15, 
  date_begin="2016-07-01", 
  date_end="2016-07-01")
print(df)
```

### Example 2
Retrieve daily CAMS solar data in netCDF format. You need to have the `ncdf4` package installed.
```{r}
library(ncdf4)

filename <- paste0(tempfile(), ".nc")

r <- cams_api(
  60, 15, "2016-06-01", "2016-06-10", 
  format="application/x-netcdf",
  time_step = "P01D",
  filename=filename)

# Access the on disk stored ncdf4 file 
nc <- nc_open(filename)

# list names of available variables
names(nc$var)

# create data.frame with timestamp and global horizontal irradiation and plot it
df <- data.frame(
  timestamp = as.POSIXct(nc$dim$time$vals, "UTC", origin="1970-01-01"),
  GHI = ncvar_get(nc, "GHI"))

plot(df, type="l")

nc_close(nc)
```

## Meta

* This package and functions herein are provided as is, without any guarantee.
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
* Please [report any issues or bugs](https://github.com/ropensci/camsRad/issues).

<!--[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org) 
doesn´t knit. add following to the .md file
[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
-->
---
title: "Working with CAMS solar data"
author: "Lukas Lundström"
date: "2016-11-01"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Working with CAMS solar data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>", 
  collapse = TRUE,
  warning = FALSE
)
if(Sys.getenv("CAMS_USERNAME")=="") {
  knitr::opts_chunk$set(eval = FALSE) # doesn't build vignette if CAMS_USERNAME is not set
}
```

## CAMS Radiation Service
The CAMS Radiation Service provides time series of global, direct, and diffuse irradiations on horizontal surface, and direct Irradiation on normal plane for the actual weather conditions as well as for clear-sky conditions. The data can be accessed manually on the CAMS Radiation Service [site](http://www.soda-pro.com/web-services/radiation/cams-radiation-service). The service is part of the [Copernicus Atmosphere Monitoring Service](http://atmosphere.copernicus.eu/) (CAMS).

* The geographical coverage is the field-of-view of the Meteosat satellite, roughly speaking Europe, Africa, Atlantic Ocean, Middle East (-66 to 66 degrees in both latitudes and longitudes).
* The time coverage of data is from 2004-02-01 up to 2 days ago. Data are available with a time step ranging from 15 min to 1 month. The service allows time steps of 1 min as well, these will be interpolated values from the 15 min time steps (accounting for changes in clear-sky conditions). 
* The horizontal resolution is the original pixel of the Meteosat images, 3-5 km.  The CAMS Radiation Service has currently (as of 2016-11-01) limited the amount of requests to 15 per day. This limit may evolve. 

## Authentication

To access the CAMS Radiation Service you need to register at [http://www.soda-pro.com/web-services/radiation/cams-radiation-service](http://www.soda-pro.com/web-services/radiation/cams-radiation-service). The email you use at the registration step will be used for authentication, and need to be set with  `cams_set_user()`.

```{r eval=FALSE}
# Authentication
cams_set_user("your@email.com") # An email registered at soda-pro.com
```


## Retrieve hourly solar data

`cams_get_radiation()` and `cams_get_radiation()` are convenience wrappers that retrieves CAMS solar data straight into a R data frame. The example bellow retrieves hourly radiation data for the location 60° latitude and 15° longitude for the period 2016-01-01 to 2016-01-15.
```{r}  
library(camsRad)

df <- cams_get_radiation(
  lat=60, lng=15, # Latitude & longitude coordinates 
  date_begin="2016-07-01", date_end="2016-07-01", # Date range
  time_step="PT01H") # Use hourly time step

```
As seen above the `cams_get_radiation()` prints additional information about the data, these can be suppressed by wrapping the call with `suppressMessages()`. Next the date frame is printed.

```{r}
print(df)
```
The first column holds the timestamp information. It follows the convention of giving solar radiation as the sum during the previous hour. E.g. the timestamp of 14:00 shows the solar radiation during 13:00-14:00. 

## Advanced retrievals
To use other data formats and to save data to the disk we need to use the `cams_api()`. The example bellow writes daily solar radiation in netCDF format to the disk. You need to have the `ncdf4` package installed.

```{r, warning = FALSE, fig.width=7}
library(ncdf4)

filename <- paste0(tempfile(), ".nc")

r <- cams_api(
  60, 15, "2016-06-01", "2016-07-1", # Latitude/longitude and date range 
  format="application/x-netcdf", # specifies output format as netCDF
  time_step = "P01D", # daily sum is specified
  filename=filename)

# Access the on disk stored netCDF file
nc <- nc_open(filename)  

# List names of available variables
names(nc$var)

# Create data.frame with datetime and global horizontal irradiation
df <- data.frame(
  timestamp = as.POSIXct(nc$dim$time$vals, "UTC", origin="1970-01-01"),
  GHI = ncvar_get(nc, "GHI"))
df$timestamp <- df$timestamp-3600*24 # shift timestamp 24 hours backwards

nc_close(nc) # close connection

# And plot it
par(mar=c(4.5,4,0.8,1))
plot(df, type="b", ylab="GHI, Wh/m2,day", xlab="2016")


```
Note that the *timestamp* follows the convention of giving solar radiation as the sum during the previous time step. This is often correct when working with hourly data. But when working with daily (or monthly) data it is more common to have the *timestamp* at the starting point of summation. The `df$timestamp-3600*24`part achieves this for daily data.

To get the data in *csv* or *json* format instead of *netCDF*, just change the *format* parameter to "application/csv" or "application/json" (and the filename extension to *.csv* or *.json* respectively).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cams_api.R
\name{cams_api}
\alias{cams_api}
\title{API client for
\href{http://www.soda-pro.com/web-services/radiation/cams-radiation-service}{CAMS
radiation service}}
\usage{
cams_api(lat, lng, date_begin, date_end, alt = -999,
  time_step = "PT01H", time_ref = "UT", verbose = FALSE,
  service = "get_cams_radiation", format = "application/csv",
  filename = "")
}
\arguments{
\item{lat}{Latitude, in decimal degrees. Required}

\item{lng}{Longitude, in decimal degrees. Required}

\item{date_begin}{Start date as 'yyyy-mm-dd' string. Required}

\item{date_end}{End date as 'yyyy-mm-dd' string. Required}

\item{alt}{Altitude in meters, use -999 to let CAMS decide. Default -99}

\item{time_step}{Aggregation: 'PT01M' for minutes, 'PT15M' for 15 minutes,
'PT01H' for hourly, 'P01D' for daily, 'P01M' for monthly. Deafult 'PT01H'}

\item{time_ref}{Time reference:'UT' for universal time, 'TST' for true solar
time. Default 'UT'}

\item{verbose}{TRUE for verbose output. Default "FALSE"}

\item{service}{'get_mcclear' for CAMS McClear data, 'get_cams_radiation' for
CAMS radiation data. Default 'get_cams_radiation'}

\item{format}{'application/csv', 'application/json', 'application/x-netcdf'
or 'text/csv'. Default 'application/csv'}

\item{filename}{path to file on disk to write to. If empty, data is kept in
memory. Default empty}
}
\value{
list(ok=TRUE/FALSE, response=response). If ok=TRUE, response is the
  response from httr::GET. If ok=FALSE, response holds exception text
}
\description{
API client for
\href{http://www.soda-pro.com/web-services/radiation/cams-radiation-service}{CAMS
radiation service}
}
\examples{
\dontrun{
library(ncdf4)

filename <- paste0(tempfile(), ".nc")

# API call to CAMS
r <- cams_api(
  60, 15,                       # latitude=60, longitude=15
  "2016-06-01", "2016-06-10",   # for 2016-06-01 to 2016-06-10
  time_step="PT01H",            # hourly data
  service="get_cams_radiation", # CAMS radiation
  format="application/x-netcdf",# netCDF format
  filename=filename)            # file to save to

# Access the on disk stored ncdf4 file
nc <- nc_open(filename)
# list names of available variables
names(nc$var)

# create data.frame with timestamp and global horizontal irradiation
df <- data.frame(datetime=as.POSIXct(nc$dim$time$vals, "UTC",
                                     origin="1970-01-01"),
                 GHI = ncvar_get(nc, "GHI"))

plot(df, type="l")

nc_close(nc)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cams_get.R
\name{cams_get_radiation}
\alias{cams_get_radiation}
\title{Retrieve CAMS solar radiation data}
\usage{
cams_get_radiation(lat, lng, date_begin, date_end, time_step = "PT01H",
  alt = -999, verbose = FALSE)
}
\arguments{
\item{lat}{Latitude, in decimal degrees. Required}

\item{lng}{Longitude, in decimal degrees. Required}

\item{date_begin}{Start date as 'yyyy-mm-dd' string. Required}

\item{date_end}{End date as 'yyyy-mm-dd' string. Required}

\item{time_step}{Aggregation: 'PT01M' for minutes, 'PT15M' for 15 minutes,
'PT01H' for hourly, 'P01D' for daily, 'P01M' for monthly. Deafult 'PT01H'}

\item{alt}{Altitude in meters, use -999 to let CAMS decide. Default -99}

\item{verbose}{TRUE for verbose output. Default "FALSE"}
}
\value{
A data frame with requested solar data
}
\description{
Retrieve CAMS solar radiation data
}
\examples{
\dontrun{
# Get hourly solar radiation data
df <- cams_get_radiation(
  lat=60, lng=15,
  date_begin="2016-06-01", date_end="2016-06-15")
head(df)

# Get daily solar radiation data
df <- cams_get_radiation(
  lat=60, lng=15,
  date_begin="2016-06-01", date_end="2016-06-15",
  time_step="P01D")
head(df)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cams_auth.R
\name{cams_set_user}
\alias{cams_set_user}
\title{Set username used for authentication by CAMS radiation service}
\usage{
cams_set_user(username)
}
\arguments{
\item{username}{Email registered at soda-pro.com. Required}
}
\description{
Set username used for authentication by CAMS radiation service
}
\examples{
\dontrun{
# cams_set_user("your@email.com") # An email registered at soda-pro.com
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cams_get.R
\name{cams_get_mcclear}
\alias{cams_get_mcclear}
\title{Retrieve McClear clear sky solar radiation data}
\usage{
cams_get_mcclear(lat, lng, date_begin, date_end, time_step = "PT01H",
  alt = -999, verbose = FALSE)
}
\arguments{
\item{lat}{Latitude, in decimal degrees. Required}

\item{lng}{Longitude, in decimal degrees. Required}

\item{date_begin}{Start date as 'yyyy-mm-dd' string. Required}

\item{date_end}{End date as 'yyyy-mm-dd' string. Required}

\item{time_step}{Aggregation: 'PT01M' for minutes, 'PT15M' for 15 minutes,
'PT01H' for hourly, 'P01D' for daily, 'P01M' for monthly. Deafult 'PT01H'}

\item{alt}{Altitude in meters, use -999 to let CAMS decide. Default -99}

\item{verbose}{TRUE for verbose output. Default "FALSE"}
}
\value{
A data frame with requested solar data
}
\description{
Retrieve McClear clear sky solar radiation data
}
\examples{
\dontrun{
df <- cams_get_mcclear(
  lat=60, lng=15, date_begin="2016-01-01", date_end="2016-01-15")
print(head(df))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/camsRad.R
\docType{package}
\name{camsRad}
\alias{camsRad}
\alias{camsRad-package}
\title{R client for CAMS radiation service}
\description{
CAMS radiation service provides time series of global, direct, and diffuse
irradiations on horizontal surface, and direct irradiation on normal plane
for the actual weather conditions as well as for clear-sky conditions. The
geographical coverage is the field-of-view of the Meteosat satellite, roughly
speaking Europe, Africa, Atlantic Ocean, Middle East (-66 to 66 degrees in
both latitudes and longitudes). The time coverage of data is from 2004-02-01
up to 2 days ago. Data are available with a time step ranging from 15 min to
1 month.
}
