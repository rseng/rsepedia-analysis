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
---
title: "hddtools: Hydrological Data Discovery Tools"
author: "Claudia Vitolo"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{hddtools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r init, echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  cache = FALSE,
  eval = FALSE
)
```

# Introduction

**hddtools** stands for Hydrological Data Discovery Tools. This R package is an open source project designed to facilitate access to a variety of online open data sources relevant for hydrologists and, in general, environmental scientists and practitioners. 

This typically implies the download of a metadata catalogue, selection of information needed, formal request for dataset(s), de-compression, conversion, manual filtering and parsing. All those operation are made more efficient by re-usable functions. 

Depending on the data license, functions can provide offline and/or online modes. When redistribution is allowed, for instance, a copy of the dataset is cached within the package and updated twice a year. This is the fastest option and also allows offline use of package's functions. When re-distribution is not allowed, only online mode is provided.

## Installation
Get the released version from CRAN:
  
```{r installation1}
install.packages("hddtools")
```

Or the development version from github using [devtools](https://github.com/r-lib/devtools):
  
```{r installation2}
devtools::install_github("ropensci/hddtools")
```

Load the hddtools package:
  
```{r loading, eval = TRUE}
library("hddtools")
```

## Data sources and Functions

The functions provided can retrieve hydrological information from a variety of data providers.
To filter the data, it is advisable to use the package `dplyr`.

```{r dplyr, eval = TRUE}
library("dplyr")
```

### The Koppen Climate Classification map
The Koppen Climate Classification is the most widely used system for classifying the world's climates. Its categories are based on the annual and monthly averages of temperature and precipitation. It was first updated by Rudolf Geiger in 1961, then by Kottek et al. (2006), Peel et al. (2007) and then by Rubel et al. (2010). 

The package hddtools contains a function to identify the updated Koppen-Greiger climate zone, given a bounding box.

```{r KGClimateClass1, eval = TRUE}
# Define a bounding box
areaBox <- raster::extent(-10, 5, 48, 62)

# Extract climate zones from Peel's map:
KGClimateClass(areaBox = areaBox, updatedBy = "Peel")
```

```{r KGClimateClass2, eval = TRUE}
# Extract climate zones from Kottek's map:
KGClimateClass(areaBox = areaBox, updatedBy = "Kottek")
```

### The Global Runoff Data Centre
The Global Runoff Data Centre (GRDC) is an international archive hosted by the Federal Institute of Hydrology in Koblenz, Germany. The Centre operates under the auspices of the World Meteorological Organisation and retains services and datasets for all the major rivers in the world. Catalogue, kml files and the product Long-Term Mean Monthly Discharges are open data and accessible via the hddtools.

Information on all the GRDC stations can be retrieved using the function `catalogueGRDC` with no input arguments, as in the examle below:  
```{r catalogueGRDC1, eval = FALSE}
# GRDC full catalogue
GRDC_catalogue <- catalogueGRDC()
```

It is advisable to use the package `dplyr` for convenient filtering, some examples are provided below.

```{r catalogueGRDC2, eval = FALSE}
# Filter GRDC catalogue based on a country code
GRDC_catalogue %>%
  filter(country == "IT")

# Filter GRDC catalogue based on rivername
GRDC_catalogue %>%
  filter(river == "PO, FIUME")

# Filter GRDC catalogue based on which daily data is available since 2000
GRDC_catalogue %>%
  filter(d_start >= 2000)

# Filter the catalogue based on a geographical bounding box
GRDC_catalogue %>%
  filter(between(x = long, left = -10, right = 5),
         between(x = lat, left = 48, right = 62))

# Combine filtering criteria
GRDC_catalogue %>%
  filter(between(x = long, left = -10, right = 5),
         between(x = lat, left = 48, right = 62),
         d_start >= 2000,
         area > 1000)
```

The GRDC catalogue (or a subset) can be used to create a map.
```{r catalogueGRDC7, eval = FALSE}
# Visualise outlets on an interactive map
library(leaflet)
leaflet(data = GRDC_catalogue %>% filter(river == "PO, FIUME")) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(~long, ~lat, popup = ~station)
```
![](leaflet.png)

### Top-Down modelling Working Group (Data60UK and MOPEX)
The Top-Down modelling Working Group (TDWG) for the Prediction in Ungauged Basins (PUB) Decade (2003-2012) is an initiative of the International Association of Hydrological Sciences (IAHS) which collected datasets for hydrological modelling free-of-charge, available [here](http://tdwg.catchment.org/datasets.html). This package provides a common interface to retrieve, browse and filter information.

#### The Data60UK dataset
The Data60UK initiative collated datasets of areal precipitation and streamflow discharge across 61 gauging sites in England and Wales (UK). The database was prepared from source databases for research purposes, with the intention to make it re-usable. This is now available in the public domain free of charge. 

The hddtools contain two functions to interact with this database: one to retreive the catalogue and another to retreive time series of areal precipitation and streamflow discharge.

```{r catalogueData60UK1, eval = TRUE}
# Data60UK full catalogue
Data60UK_catalogue_all <- catalogueData60UK()

# Filter Data60UK catalogue based on bounding box
areaBox <- raster::extent(-4, -3, 51, 53)
Data60UK_catalogue_bbox <- catalogueData60UK(areaBox = areaBox)
```

```{r catalogueData60UK2}
# Visualise outlets on an interactive map
library(leaflet)
leaflet(data = Data60UK_catalogue_bbox) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(~Longitude, ~Latitude, popup = ~Location)
```

![](leaflet2.png)

```{r catalogueData60UK3, eval = TRUE, message = FALSE, fig.width = 7, fig.height = 7}
# Extract time series 
id <- catalogueData60UK()$id[1]

# Extract only the time series
MorwickTS <- tsData60UK(id)
```

#### MOPEX
The MOPEX dataset contains river basin characteristics and data for hundreds of river basins from a range of climates in the US. The catalogue can be retrieved as follows:
```{r MOPEX_meta, eval = FALSE, message = FALSE, fig.width = 7, fig.height = 7}
# MOPEX full catalogue
MOPEX_catalogue <- catalogueMOPEX()

# Extract data within a geographic bounding box
MOPEX_catalogue %>%
  filter(dplyr::between(x = Longitude, left = -95, right = -92),
         dplyr::between(x = Latitude, left = 37, right = 41))
```

```{r MOPEX_meta2, eval = FALSE, message = FALSE, fig.width = 7, fig.height = 7}
# Get stations with recondings in the period 1st Jan to 31st Dec 1995
MOPEX_catalogue %>%
  filter(Date_start <= as.Date("1995-01-01"),
         Date_end >= as.Date("1995-12-31"))

# Get only catchments within NC
MOPEX_catalogue %>%
  filter(State == "NC")
```

For each station, historical hydrometeorological data  can also be retrieved.

```{r MOPEX_data, eval = FALSE, message = FALSE, fig.width = 7, fig.height = 7}
# Take the first record in the catalogue
river_metadata <- MOPEX_catalogue[1,]

# Get corresponding time series
river_ts <- tsMOPEX(id = river_metadata$USGS_ID)

# Extract data between 1st Jan and 31st December 1948
river_ts_shorter <- window(river_ts,
                           start = as.Date("1948-01-01"),
                           end = as.Date("1948-12-31"))

# Plot
plot(river_ts_shorter,
     main = river_metadata$Name,
     xlab = "",
     ylab = c("P [mm/day]","E [mm/day]", "Q [mm/day]", "Tmax [C]","Tmin [C]"))
```

![](mopex.png)

### SEPA river level data
The Scottish Environment Protection Agency (SEPA) manages river level data for hundreds of gauging stations in the UK. The catalogue of stations is derived from the list here:
https://www2.sepa.org.uk/waterlevels/CSVs/SEPA_River_Levels_Web.csv.

```{r SEPA1, eval = FALSE}
# SEPA catalogue
SEPA_catalogue <- catalogueSEPA()
```

The time series of the last few days is available from SEPA website and can be downloaded using the following function:
```{r SEPA2, eval = FALSE, message = FALSE, fig.width = 7}
# Take the first record in the catalogue
Perth_metadata <- SEPA_catalogue[1,]

# Single time series extraction
Perth_ts <- tsSEPA(id = Perth_metadata$LOCATION_CODE)

# Plot
plot(Perth_ts,
     main = Perth_metadata$STATION_NAME,
     xlab = "",
     ylab = "Water level [m]")

# Get only catchments with area above 4000 Km2
SEPA_catalogue %>%
  filter(CATCHMENT_AREA >= 4000)

# Get only catchments within river Ayr
SEPA_catalogue %>%
  filter(RIVER_NAME == "Ayr")
```

Plese note that these data are updated every 15 minutes and the code will always generate different plots. 
```{r SEPA3, eval=FALSE, message = FALSE, fig.width = 7}
# Multiple time series extraction
y <- tsSEPA(id = c("234253", "234174", "234305"))

plot(y[[1]], ylim = c(0, max(y[[1]], y[[2]], y[[3]])), 
     xlab = "", ylab = "Water level [m]")
lines(y[[2]], col = "red")
lines(y[[3]], col = "blue")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRDC.R
\name{catalogueGRDC}
\alias{catalogueGRDC}
\title{Data source: Global Runoff Data Centre catalogue}
\usage{
catalogueGRDC()
}
\value{
This function returns a data frame made with the following columns:
\itemize{
  \item{\code{grdc_no}}{: GRDC station number}
  \item{\code{wmo_reg}}{: WMO region}
  \item{\code{sub_reg}}{: WMO subregion}
  \item{\code{river}}{: river name}
  \item{\code{station}}{: station name}
  \item{\code{country}}{: 2-letter country code (ISO 3166)}
  \item{\code{lat}}{: latitude, decimal degree}
  \item{\code{long}}{: longitude, decimal degree}
  \item{\code{area}}{: catchment size, km2}
  \item{\code{altitude}}{: height of gauge zero, m above sea level}
  \item{\code{d_start}}{: daily data available from year}
  \item{\code{d_end}}{: daily data available until year}
  \item{\code{d_yrs}}{: length of time series, daily data}
  \item{\code{d_miss}}{: percentage of missing values (daily data)}
  \item{\code{m_start}}{: monthly data available from}
  \item{\code{m_end}}{: monthly data available until}
  \item{\code{m_yrs}}{: length of time series, monthly data}
  \item{\code{m_miss}}{: percentage of missing values (monthly data)}
  \item{\code{t_start}}{: earliest data available}
  \item{\code{t_end}}{: latest data available}
  \item{\code{t_yrs}}{: maximum length of time series, daily and monthly
  data}
  \item{\code{lta_discharge}}{: mean annual streamflow, m3/s}
  \item{\code{r_volume_yr}}{: mean annual volume, km3}
  \item{\code{r_height_yr}}{: mean annual runoff depth, mm}
}
}
\description{
This function interfaces the Global Runoff Data Centre database
which provides river discharge data for almost 1000 sites over 157 countries.
}
\examples{
\dontrun{
  # Retrieve the catalogue
  GRDC_catalogue_all <- catalogueGRDC()
}

}
\author{
Claudia Vitolo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hddtools-package.R
\docType{data}
\name{grdcLTMMD}
\alias{grdcLTMMD}
\title{Data set: The grdcLTMMD look-up table}
\format{
A data frame with 6 rows and 4 columns.
\describe{
  \item{\code{WMO_Region}}{an integer between 1 and 6}
  \item{\code{Coverage}}{}
  \item{\code{Number_of_stations}}{}
  \item{\code{Archive}}{url to spreadsheet}
}
}
\source{
\url{http://www.bafg.de/GRDC}
}
\usage{
data("grdcLTMMD")
}
\description{
The grdcLTMMD look-up table
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hddtools-package.R
\docType{package}
\name{hddtools}
\alias{hddtools}
\title{hddtools: Hydrological Data Discovery Tools}
\description{
Many governmental bodies and institutions are currently committed to publish open data as the result of a trend of increasing transparency, based on which a wide variety of information produced at public expense is now becoming open and freely available to improve public involvement in the process of decision and policy making. Discovery, access and retrieval of information is, however, not always a simple task. Especially when access to data APIs is not allowed, downloading a metadata catalogue, selecting the information needed, requesting datasets, de-compression, conversion, manual filtering and parsing can become rather tedious. The R package hddtools is an open source project, designed to make all the above operations more efficient by means of reusable functions.

The package facilitate access to various online data sources such as:
\itemize{
 \item{\strong{KGClimateClass} (\url{http://koeppen-geiger.vu-wien.ac.at/}): The Koppen Climate Classification map is used for classifying the world's climates based on the annual and monthly averages of temperature and precipitation}
 \item{\strong{GRDC} (\url{http://www.bafg.de/GRDC/EN/Home/homepage_node.html}): The Global Runoff Data Centre (GRDC) provides datasets for all the major rivers in the world}
 \item{\strong{Data60UK} (\url{http://tdwg.catchment.org/datasets.html}): The Data60UK initiative collated datasets of areal precipitation and streamflow discharge across 61 gauging sites in England and Wales (UK).}
 \item{\strong{MOPEX} (\url{https://www.nws.noaa.gov/ohd/mopex/mo_datasets.htm}): This dataset contains historical hydrometeorological data and river basin characteristics for hundreds of river basins in the US.}
 \item{\strong{SEPA} (\url{https://www2.sepa.org.uk/WaterLevels/}): The Scottish Environment Protection Agency (SEPA) provides river level data for hundreds of gauging stations in the UK.}}
This package complements R's growing functionality in environmental web technologies by bridging the gap between data providers and data consumers. It is designed to be an initial building block of scientific workflows for linking data and models in a seamless fashion.
}
\references{
Vitolo C, Buytaert W, 2014, HDDTOOLS: an R package serving Hydrological Data
Discovery Tools, AGU Fall Meeting, 15-19 December 2014, San Francisco, USA.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data60UK.R
\name{catalogueData60UK}
\alias{catalogueData60UK}
\title{Data source: Data60UK catalogue}
\source{
\url{http://nrfaapps.ceh.ac.uk/datauk60/data.html}
}
\usage{
catalogueData60UK(areaBox = NULL)
}
\arguments{
\item{areaBox}{bounding box, a list made of 4 elements: minimum longitude
(lonMin), minimum latitude (latMin), maximum longitude (lonMax), maximum
latitude (latMax)}
}
\value{
This function returns a data frame containing the following columns:
\describe{
  \item{\code{id}}{Station id number.}
  \item{\code{River}}{String describing the river's name.}
  \item{\code{Location}}{String describing the location.}
  \item{\code{gridReference}}{British National Grid Reference.}
  \item{\code{Latitude}}{}
  \item{\code{Longitude}}{}
}
}
\description{
This function interfaces the Data60UK database catalogue
listing 61 gauging stations.
}
\examples{
\dontrun{
  # Retrieve the whole catalogue
  Data60UK_catalogue_all <- catalogueData60UK()

  # Filter the catalogue based on a bounding box
  areaBox <- raster::extent(-4, -2, +52, +53)
  Data60UK_catalogue_bbox <- catalogueData60UK(areaBox)
}

}
\author{
Claudia Vitolo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOPEX.R
\name{tsMOPEX}
\alias{tsMOPEX}
\title{Interface for the MOPEX database of Daily Time Series}
\usage{
tsMOPEX(id, MAP = TRUE)
}
\arguments{
\item{id}{String for the station ID number (USGS_ID)}

\item{MAP}{Boolean, TRUE by default. If FALSE it looks for data through all
the 1861 potential MOPEX basins.
If TRUE, it looks for data through the 438 MOPEX basins with MAP estimates.}
}
\value{
If MAP = FALSE, this function returns a time series of daily
streamflow discharge (Q, in mm). If MAP = TRUE, this function returns a data
frame containing the following columns (as zoo object):
\describe{
  \item{\code{Date}}{Format is "yyyymmdd"}
  \item{\code{P}}{Mean areal precipitation (mm)}
  \item{\code{E}}{Climatic potential evaporation (mm, based NOAA Freewater Evaporation Atlas)}
  \item{\code{Q}}{Daily streamflow discharge (mm)}
  \item{\code{T_max}}{Daily maximum air temperature (Celsius)}
  \item{\code{T_min}}{Daily minimum air temperature (Celsius)}
}
}
\description{
This function extract the dataset containing daily rainfall and
streamflow discharge at one of the MOPEX locations.
}
\examples{
\dontrun{
  BroadRiver <- tsMOPEX(id = "01048000")
}

}
\author{
Claudia Vitolo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bboxSpatialPolygon.R
\name{bboxSpatialPolygon}
\alias{bboxSpatialPolygon}
\title{Convert a bounding box to a SpatialPolygons object
Bounding box is first created (in lat/lon) then projected if specified}
\usage{
bboxSpatialPolygon(boundingbox, proj4stringFrom = NULL, proj4stringTo = NULL)
}
\arguments{
\item{boundingbox}{Bounding box: a 2x2 numerical matrix of lat/lon
coordinates}

\item{proj4stringFrom}{Projection string for the current boundingbox
coordinates (defaults to lat/lon, WGS84)}

\item{proj4stringTo}{Projection string, or NULL to not project}
}
\value{
A SpatialPolygons object of the bounding box
}
\description{
Convert a bounding box to a SpatialPolygons object
Bounding box is first created (in lat/lon) then projected if specified
}
\examples{
\dontrun{
  boundingbox <- raster::extent(-180, +180, -50, +50)
  bbSP <- bboxSpatialPolygon(boundingbox = boundingbox)
}

}
\references{
\url{https://gis.stackexchange.com/questions/46954/clip-spatial-object-to-bounding-box-in-r}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data60UK.R
\name{tsData60UK}
\alias{tsData60UK}
\title{Interface for the Data60UK database of Daily Time Series}
\usage{
tsData60UK(id)
}
\arguments{
\item{id}{String which identifies the station ID number}
}
\value{
The function returns a data frame containing 2 time series (as zoo objects): "P" (precipitation) and "Q" (discharge).
}
\description{
This function extract the dataset containing daily rainfall and streamflow discharge at one of the Data60UK locations.
}
\examples{
\dontrun{
  Morwick <- tsData60UK(id = "22001")
}

}
\author{
Claudia Vitolo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEPA.R
\name{tsSEPA}
\alias{tsSEPA}
\title{Interface for the MOPEX database of Daily Time Series}
\usage{
tsSEPA(id)
}
\arguments{
\item{id}{hydrometric reference number (string)}
}
\value{
The function returns river level data in metres, as a zoo object.
}
\description{
This function extract the dataset containing daily rainfall and
streamflow discharge at one of the MOPEX locations.
}
\examples{
\dontrun{
  sampleTS <- tsSEPA(id = "10048")
}

}
\author{
Claudia Vitolo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KGClimateClass.R
\name{KGClimateClass}
\alias{KGClimateClass}
\title{Function to identify the updated Koppen-Greiger climate zone (on a 0.1 x 0.1 degrees resolution map).}
\usage{
KGClimateClass(areaBox = NULL, updatedBy = "Peel", verbose = FALSE)
}
\arguments{
\item{areaBox}{bounding box, a list made of 4 elements: minimum longitude (lonMin), minimum latitude (latMin), maximum longitude (lonMax), maximum latitude (latMax)}

\item{updatedBy}{this can either be "Kottek" or "Peel"}

\item{verbose}{if TRUE more info are printed on the screen}
}
\value{
List of overlapping climate zones.
}
\description{
Given a bounding box, the function identifies the overlapping climate zones.
}
\examples{
\dontrun{
  # Define a bounding box
  areaBox <- raster::extent(-3.82, -3.63, 52.41, 52.52)
  # Get climate classes
  KGClimateClass(areaBox = areaBox)
}

}
\references{
Kottek et al. (2006): \url{http://koeppen-geiger.vu-wien.ac.at/}. Peel et al. (2007): \url{https://people.eng.unimelb.edu.au/mpeel/koppen.html}.
}
\author{
Claudia Vitolo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOPEX.R
\name{catalogueMOPEX}
\alias{catalogueMOPEX}
\title{Data source: MOPEX catalogue}
\source{
https://hydrology.nws.noaa.gov/pub/gcip/mopex/US_Data/Documentation/
}
\usage{
catalogueMOPEX(MAP = TRUE)
}
\arguments{
\item{MAP}{Boolean, TRUE by default. If FALSE it returns a list of the USGS
station ID’s and the gage locations of all 1861 potential MOPEX basins.
If TRUE, it return a list of the USGS station ID’s and the gage locations of
the 438 MOPEX basins with MAP estimates.}
}
\value{
This function returns a data frame containing the following columns:
\describe{
  \item{\code{USGS_ID}}{Station id number}
  \item{\code{Longitude}}{Decimal degrees East}
  \item{\code{Latitude}}{Decimal degrees North}
  \item{\code{Drainage_Area}}{Square Miles}
  \item{\code{R_gauges}}{Required number of precipitation gages to meet MAP accuracy criteria}
  \item{\code{N_gauges}}{Number of gages in total gage window used to estimate MAP}
  \item{\code{A_gauges}}{Avaliable number of gages in the basin}
  \item{\code{Ratio_AR}}{Ratio of Available to Required number of gages in the basin}
  \item{\code{Date_start}}{Date when recordings start}
  \item{\code{Date_end}}{Date when recordings end}
  \item{\code{State}}{State of the basin}
  \item{\code{Name}}{Name of the basin}
}
Columns Date_start, Date_end, State, Name are taken from:
https://hydrology.nws.noaa.gov/pub/gcip/mopex/US_Data/Basin_Characteristics/usgs431.txt
Date_start and Date_end are conventionally set to the first of the month
here, however actual recordings my differ. Always refer to the recording date
obtained as output of \code{tsMOPEX()}.
}
\description{
This function retrieves the list of the MOPEX basins.
}
\examples{
\dontrun{
  # Retrieve the MOPEX catalogue
  catalogue <- catalogueMOPEX()
}

}
\author{
Claudia Vitolo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEPA.R
\name{catalogueSEPA}
\alias{catalogueSEPA}
\title{Data source: SEPA catalogue}
\usage{
catalogueSEPA()
}
\value{
This function returns a data frame containing the following columns:
\describe{
  \item{\code{SEPA_HYDROLOGY_OFFICE}}{}
  \item{\code{STATION_NAME}}{}
  \item{\code{LOCATION_CODE}}{Station id number.}
  \item{\code{NATIONAL_GRID_REFERENCE}}{}
  \item{\code{CATCHMENT_NAME}}{}
  \item{\code{RIVER_NAME}}{}
  \item{\code{GAUGE_DATUM}}{}
  \item{\code{CATCHMENT_AREA}}{in Km2}
  \item{\code{START_DATE}}{}
  \item{\code{END_DATE}}{}
  \item{\code{SYSTEM_ID}}{}
  \item{\code{LOWEST_VALUE}}{}
  \item{\code{LOW}}{}
  \item{\code{MAX_VALUE}}{}
  \item{\code{HIGH}}{}
  \item{\code{MAX_DISPLAY}}{}
  \item{\code{MEAN}}{}
  \item{\code{UNITS}}{}
  \item{\code{WEB_MESSAGE}}{}
  \item{\code{NRFA_LINK}}{}
}
}
\description{
This function provides the official SEPA database catalogue of
river level data
(from https://www2.sepa.org.uk/waterlevels/CSVs/SEPA_River_Levels_Web.csv)
containing info for hundreds of stations. Some are NRFA stations.
The function has no input arguments.
}
\examples{
\dontrun{
  # Retrieve the whole catalogue
  SEPA_catalogue_all <- catalogueSEPA()
}

}
\author{
Claudia Vitolo
}
