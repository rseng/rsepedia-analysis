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
---
title: "rdefra: Interact with the UK AIR Pollution Database from DEFRA"
author: "Claudia Vitolo"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{rdefra: Interact with the UK AIR Pollution Database from DEFRA}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
  
```{r setup, echo = FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  eval = FALSE
)
```

# Introduction
The package rdefra allows to retrieve air pollution data from the Air Information Resource (UK-AIR, https://uk-air.defra.gov.uk/) of the Department for Environment, Food and Rural Affairs (DEFRA) in the United Kingdom. UK-AIR does not provide a public API for programmatic access to data, therefore this package scrapes the HTML pages to get relevant information.

This package follows a logic similar to other packages such as waterData and rnrfa: sites are first identified through a catalogue, data are imported via the station identification number, then visualised and/or used in analyses. The information related to the monitoring stations is accessible through the function `ukair_catalogue()`. Some station may have missing coordinates, which can be recovered using the function `ukair_get_coordinates()`. Lastly, time series data related to different pollutants can be obtained using the function `ukair_get_hourly_data()`.

The package is designed to collect data efficiently. It allows to download multiple years of data for a single station with one line of code and, if used with the parallel package, allows the acquisition of data from hundreds of sites in only few minutes.

## Installation

Get the released version from CRAN:
  
```{r installation_cran}
install.packages("rdefra")
```

Or the development version from GitHub using the package `remotes`:

```{r installation_github}
install.packages("remotes")
remotes::install_github("ropensci/rdefra")
```

Load the rdefra package:
  
```{r load_library, eval = TRUE}
library("rdefra")
```

## Functions

The package logic assumes that users access the UK-AIR database in two steps:
  
  1. Browse the catalogue of available stations and selects some stations of interest.
  2. Retrieves data for the selected stations.

### Get stations catalogue

The list of monitoring stations can be downloaded using the function `ukair_catalogue()` with no input parameters, as in the example below. 

```{r catalogue_full, eval = TRUE}
# Get full catalogue
stations <- ukair_catalogue()
head(stations)
```

There are currently `r dim(stations)[1]` stations in UK-AIR. The same function, can be used to filter the catalogue using the following input parameters:

  * `site_name` IDs of specific site (UK.AIR.ID). By default this is left blank to get info on all the available sites.
  * `pollutant` This is an integer between 1 and 10. Default is 9999, which means all the pollutants.
  * `group_id` This is the identification number of a group of stations. Default is 9999 which means all available networks.
  * `closed` This is set to TRUE to include closed stations, FALSE otherwise.
  * `country_id` This is the identification number of the country, it can be an integer between 1 and 6. Default is 9999, which means all the countries.
  * `region_id` This is the identification number of the region. 1 = Aberdeen City, etc. (for the full list see https://uk-air.defra.gov.uk/). Default is 9999, which means all the local authorities.

```{r catalogue_filter, eval = TRUE}
stations_EnglandOzone <- ukair_catalogue(pollutant = 1, country_id = 1)
head(stations_EnglandOzone)
```

The example above shows how to retrieve the `r dim(stations_EnglandOzone)[1]` stations in England in which ozone is measured.

### Get missing coordinates

Locating a station is extremely important to be able to carry out any spatial analysis. If coordinates are missing, for some stations in the catalogue, it might be possible to retrieve Easting and Northing coordinates (British National Grid) from DEFRA web pages, transform them to latitude and longitude and populate the missing coordinates as shown below.

```{r get_coords1, eval = TRUE}
# How many stations have missing coordinates?
length(which(is.na(stations$Latitude) | is.na(stations$Longitude)))
```

```{r get_coords2, eval = FALSE}
# Scrape DEFRA website to get Easting/Northing (if available)
stations <- ukair_get_coordinates(stations)

# How many stations still have missing coordinates?
length(which(is.na(stations$Latitude) | is.na(stations$Longitude)))
#> [1] 2
```

### Check hourly data availability

Pollution data started to be collected in 1972 and consists of hourly concentration of various species (in micrograms/m<sup>3</sup>), such as ozone (O<sub>3</sub>), particulate matters (PM<sub>2.5</sub> and PM<sub>10</sub>), nitrogen dioxide (NO<sub>2</sub>), sulphur dioxide (SO<sub>2</sub>), and so on.

The ID under which these data are available differs from the UK.AIR.ID. The catalogue does not contain this additional station ID (called SiteID hereafter) but DEFRA's web pages contain references to both the UK.AIR.ID and the SiteID. The function below uses as input the UK.AIR.ID and outputs the SiteID, if available. 

```{r get_site_id, eval = FALSE}
stations$SiteID <- ukair_get_site_id(stations$UK.AIR.ID)
```

Please note this function takes several minutes to run.

### Cached catalogue

For convenience, a cached version of the catalogue (last updated in April 2021) is included in the package and can be loaded using the following command:

```{r load_dataset_stations}
data("stations")
```

The cached catalogue contains all the available siteIDs and coordinates and can be used offline as lookup table to find out the correspondence between the UK.AIR.ID and SiteID, as well as to investigate station characteristics.

### Get hourly data

Once the SiteID is known, time series for a given station can be retrieved in one line of code:
  
```{r get_hourly_data, eval = FALSE, fig.width = 7, fig.height = 5, fig.cap = "\\label{fig:hdata}Hourly ozone data from London Marylebone Road monitoring station in 2015"}
# Get 1 year of hourly ozone data from London Marylebone Road monitoring station
df <- ukair_get_hourly_data("MY1", years = 2015)

# Aggregate to daily means and plot
# please note we use the zoo package here because time series could be irregular
library("zoo")
my1 <- zoo(x = df$Ozone, order.by = as.POSIXlt(df$datetime))

daily_means <- aggregate(my1, as.Date(as.POSIXlt(df$datetime)), mean)

plot(daily_means, main = "", xlab = "",
     ylab = expression(paste("Ozone concentration [", mu, "g/", m^3, "]")))
```
![get_hourly_data.png](get_hourly_data.png)

The above figure \ref{fig:hdata} shows the highest concentrations happen in late spring and at the beginning of summer. In order to check whether this happens every year, we can download multiple years of data and then compare them.

The code below explores the distribution of ozone by month. The resulting box plots show that the highest concentrations usually occurr during April/May and that these vary year-by-year.

```{r ozone_data, eval = FALSE}
# Get 15 years of hourly ozone data from the same monitoring station
library("ggplot2")
library("dplyr")
library("lubridate")

df <- ukair_get_hourly_data("MY1", years = 2000:2015)

df %>%
  mutate(year = year(datetime),
         month = month(datetime),
         year_month = strftime(datetime, "%Y-%m")) %>%
  group_by(month, year_month) %>%
  summarize(ozone = mean(Ozone, na.rm=TRUE)) %>%
  na.omit %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month), y = ozone, group = month),
               outlier.shape = NA) +
  xlab("Month of the year") +
  ylab(expression(paste("Ozone concentration (", mu, "g/",m^3,")"))) +
  ggtitle("15 years of hourly ozone data from London Marylebone Road monitoring station")
```
![ozone_data.png](ozone_data.png)

## Applications

### Plotting stations' locations 

After scraping DEFRA's web pages, almost all the stations have valid coordinates. You can create an interactive map using leaflet. The code below generates a map where blue circles are all the stations with valid coordinates, while red circles show locations with available hourly data.

```{r map_data, eval = FALSE}
# Keep only station with coordinates
stations_with_coords <- stations[complete.cases(stations[, c("Longitude",
                                                             "Latitude")]), ]
# Keep only station with known SiteID
stations_with_SiteID <- which(!is.na(stations_with_coords$SiteID))

# An interactive map
library("leaflet")
leaflet(data = stations_with_coords) %>% addTiles() %>% 
  addCircleMarkers(lng = ~Longitude, 
                   lat = ~Latitude,  
                   popup = ~SiteID,
                   radius = 1, color="blue", fill = FALSE) %>%
  addCircleMarkers(lng = ~Longitude[stations_with_SiteID], 
                   lat = ~Latitude[stations_with_SiteID], 
                   radius = 0.5, color="red", 
                   popup = ~SiteID[stations_with_SiteID])
```
![map_data.png](map_data.png)

### Analyse the spatial distribution of the monitoring stations

Below are two plots showing the spatial distribution of the monitoring stations. These are concentrated largely in urban areas and mostly estimate the background level of concentration of pollutants.

```{r dotchart1, eval = FALSE, fig.width = 7, fig.height = 10, fig.cap = "\\label{fig:dotchart1}Spatial distribution of the monitoring stations across zones."}
# Zone
dotchart(as.matrix(table(stations$Zone))[,1])
```
![dotchart_zone.png](dotchart_zone.png)

```{r dotchart2, eval = FALSE, fig.width = 7, fig.height = 5, fig.cap = "\\label{fig:dotchart2}Spatial distribution of the monitoring stations across environment types."}
# Environment.Type
dotchart(as.matrix(table(stations$Environment.Type[stations$Environment.Type != "Unknown Unknown"]))[,1])
```
![dotchart_envtype.png](dotchart_envtype.png)

### Use multiple cores to speed up data retrieval from numerous sites

The acquisition of data from hundreds of sites takes only few minutes:

```{r parallel_example, eval = FALSE}
library("parallel")

# Use detectCores() to find out many cores are available on your machine
cl <- makeCluster(getOption("cl.cores", detectCores()))

system.time(myList <- parLapply(cl, stations$SiteID[stations_with_SiteID], 
                                ukair_get_hourly_data, years=1999:2016))

stopCluster(cl)

df <- bind_rows(myList)
```
\name{regions}
\alias{regions}
\docType{data}
\title{
English regions
}
\description{
This is a SpatialPolygonDataFrame containing the English regions.
}
\usage{data("regions")}
\format{
  A data frame with 9 observations and 1 variable: NAME.
}
\details{
The spatial object is in latitude and longitude coordinates.
}
\examples{
data(regions)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdefra-package.R
\docType{package}
\name{rdefra}
\alias{rdefra}
\title{rdefra: Interact with the UK AIR Pollution Database from DEFRA}
\description{
The R package rdefra allows to retrieve air pollution data from the Air
Information Resource (UK-AIR) of the Department for Environment, Food and
Rural Affairs in the United Kingdom (see \url{https://uk-air.defra.gov.uk/}).
UK-AIR does not provide public APIs for programmatic access to data,
therefore this package scrapes the HTML pages to get relevant information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ukair_get_coordinates.R
\name{ukair_get_coordinates}
\alias{ukair_get_coordinates}
\title{Get Easting and Northing coordinates from DEFRA}
\usage{
ukair_get_coordinates(ids)
}
\arguments{
\item{ids}{contains the station identification code defined by DEFRA. It can
be: a) an alphanumeric string, b) a vector of strings or c) a data frame. In
the latter case, the column containing the codes should be named "UK.AIR.ID",
all the other columns will be ignored.}
}
\value{
A data.frame containing at least five columns named "UK.AIR.ID",
"Easting", "Northing", "Latitude" and "Longitude".
}
\description{
This function takes as input the UK AIR ID and returns Easting
and Northing coordinates (British National Grid, EPSG:27700).
}
\details{
If the input is a data frame with some of the columns named
"UK.AIR.ID", "Northing" and "Easting", the function only infills missing
Northing/Easting values (if available on the relevant webpage).
}
\examples{
 \dontrun{
 # Case a: alphanumeric string
 ukair_get_coordinates("UKA12536")

 # Case b: vector of strings
 ukair_get_coordinates(c("UKA15910", "UKA15956", "UKA16663", "UKA16097"))

 # Case c: data frame
 ukair_get_coordinates(ukair_catalogue()[1:10,])
 }

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ukair_get_hourly_data.R
\name{ukair_get_hourly_data}
\alias{ukair_get_hourly_data}
\title{Get hourly data for DEFRA stations}
\usage{
ukair_get_hourly_data(site_id = NULL, years = NULL)
}
\arguments{
\item{site_id}{This is the ID of a specific site.}

\item{years}{Years for which data should be downloaded.}
}
\value{
A data.frame containing hourly pollution data.
}
\description{
This function fetches hourly data from DEFRA's air pollution
monitoring stations.
}
\details{
The measurements are generally in \eqn{\mu g/m^3} (micrograms per
cubic metre). To check the units, refer to the table of attributes (see
example below). Please double check the units on the DEFRA website, as they
might change over time.
}
\examples{
 \dontrun{
 # Get data for 1 year
 output <- ukair_get_hourly_data("ABD", 2014)

 # Get data for multiple years
 output <- ukair_get_hourly_data("ABD", 2014:2016)

 # Get units
 attributes(output)$units

 }

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ukair_get_site_id.R
\name{ukair_get_site_id}
\alias{ukair_get_site_id}
\title{Get site identification numbers for DEFRA stations}
\usage{
ukair_get_site_id(id_s)
}
\arguments{
\item{id_s}{An alphanumeric string (or vector of strings) containing the UK
AIR ID defined by DEFRA.}
}
\value{
A named vector containing the site id_s.
}
\description{
Given the UK AIR ID (from the \code{ukair_catalogue()}), this
function fetches the catalogue of monitoring stations from DEFRA's website.
}
\examples{
 \dontrun{
 ukair_get_site_id("UKA00399")
 }

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdefra-package.R
\docType{data}
\name{stations}
\alias{stations}
\title{List of all the DEFRA air quality monitoring stations with complete
coordinates}
\format{
A data frame with 6561 observations on the following 14 variables.
\describe{
  \item{\code{UK.AIR.ID}}{ID reference for monitoring stations}
  \item{\code{EU.Site.ID}}{EU.Site.ID}
  \item{\code{EMEP.Site.ID}}{EMEP.Site.ID}
  \item{\code{Site.Name}}{Site name}
  \item{\code{Environment.Type}}{a factor with levels \code{Background Rural}
  \code{Background Suburban} \code{Background Urban}
  \code{Industrial Suburban} \code{Industrial Unknown}
  \code{Industrial Urban} \code{Traffic Urban} \code{Unknown Unknown}}
  \item{\code{Zone}}{Zone}
  \item{\code{Start.Date}}{Start date}
  \item{\code{End.Date}}{End date}
  \item{\code{Latitude}}{Latitude (WGS 84)}
  \item{\code{Longitude}}{Longitude (WGS 84)}
  \item{\code{Northing}}{Northing coordinate (British National Grid)}
  \item{\code{Easting}}{Easting coordinate (British National Grid)}
  \item{\code{Altitude..m.}}{Altitude in metres above sea level}
  \item{\code{Networks}}{Monitoring Networks}
  \item{\code{AURN.Pollutants.Measured}}{Pollutant measured}
  \item{\code{Site.Description}}{Description of the site.}
  \item{\code{SiteID}}{Site ID, used to retrieve time series data.}
}
}
\source{
\url{https://uk-air.defra.gov.uk/}
}
\usage{
data("stations")
}
\description{
This is the list of all the air quality monitoring stations ever
installed in the UK and operated by DEFRA networks (as per February 2016).
As the network expands, metadata for new stations will be added.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ukair_catalogue.R
\name{ukair_catalogue}
\alias{ukair_catalogue}
\title{Get DEFRA UK-AIR stations metadata}
\usage{
ukair_catalogue(
  site_name = "",
  pollutant = 9999,
  group_id = 9999,
  closed = "true",
  country_id = 9999,
  region_id = 9999,
  location_type = 9999
)
}
\arguments{
\item{site_name}{This is the name of a specific site. By default this is left
blank to get info on all the available sites.}

\item{pollutant}{This is a number from 1 to 10. Default is 9999, which means
all the pollutants.}

\item{group_id}{This is the identification number of a group of stations.
Default is 9999 which means all available networks.}

\item{closed}{This is "true" to include closed stations, "false" otherwise.}

\item{country_id}{This is the identification number of the country, it can be
a number from 1 to 6. Default is 9999, which means all the countries.}

\item{region_id}{This is the identification number of the region. 1 =
Aberdeen City, etc. (for the full list see
\url{https://uk-air.defra.gov.uk/}). Default is 9999, which means all the
local authorities.}

\item{location_type}{This is the identification number of the location type.
Default is 9999, which means all the location types.}
}
\value{
A dataframe listing stations and related information.
}
\description{
This function fetches the catalogue of monitoring stations from
DEFRA's website.
}
\details{
The argument \code{Pollutant} is defined based on the following convention:
\itemize{
 \item{1 = Ozone (O3)}
 \item{2 = Nitrogen oxides (NOx)}
 \item{3 = Carbon monoxide (CO)}
 \item{4 = Sulphur dioxide (SO2)}
 \item{5 = Particulate Matter (PM10)}
 \item{6 = Particulate Matter (PM2.5)}
 \item{7 = PAHs}
 \item{8 = Metals in PM10}
 \item{9 = Benzene}
 \item{10 = Black Carbon}
}

The argument \code{group_id} is defined based on the following convention:
\itemize{
 \item{1 = UKEAP: Precip-Net}
 \item{2 = Air Quality Strategy Pollutants}
 \item{3 = Ammonia and Nitric Acid}
 \item{4 = Automatic Urban and Rural Monitoring Network (AURN)}
 \item{5 = Dioxins and Furans}
 \item{6 = Black Smoke & SO2}
 \item{7 = Automatic Hydrocarbon Network}
 \item{8 = Heavy Metals}
 \item{9 = Nitrogen Dioxide Diffusion Tube}
 \item{10 = PAH Andersen}
 \item{11 = Particle Size Composition}
 \item{12 = PCBs}
 \item{13 = TOMPs}
 \item{14 = Non-Automatic Hydrocarbon Network}
 \item{15 = 1,3-Butadiene Diffusion Tube}
 \item{16 = Black Carbon}
 \item{17 = Automatic Urban and Rural Monitoring Network (AURN)}
 \item{18 = Defra NO2 Diffusion Tube}
 \item{19 = PAH Digitel (solid phase)}
 \item{20 = PAH Digitel (solid+vapour)}
 \item{21 = PAH Deposition}
 \item{22 = Particle size and number}
 \item{23 = Rural Automatic Mercury network}
 \item{24 = Urban Sulphate}
 \item{25 = UKEAP: Rural NO2}
 \item{26 = Automatic Urban and Rural Monitoring Network (AURN)}
 \item{27 = UKEAP: National Ammonia Monitoring Network}
 \item{28 = UKEAP: Acid Gases & Aerosol Network}
 \item{29 = Particle Speciation (MARGA)}
 \item{30 = UKEAP: Historic Aerosol measurements}
}

The argument \code{country_id} is defined based on the following convention:
\itemize{
 \item{1 = England}
 \item{2 = Wales}
 \item{3 = Scotland}
 \item{4 = Northern Ireland}
 \item{5 = Republic of Ireland}
 \item{6 = Channel Islands}
 }
}
\examples{
 \dontrun{
 stations <- ukair_catalogue()
 }

}
