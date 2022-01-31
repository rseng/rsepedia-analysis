
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN
version](https://www.r-pkg.org/badges/version/FedData)](https://cran.r-project.org/package=FedData)
[![CRAN downloads per
month](https://cranlogs.r-pkg.org/badges/FedData)](https://github.com/r-hub/cranlogs.app)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/FedData)](https://github.com/r-hub/cranlogs.app)
[![Build
Status](https://api.travis-ci.org/ropensci/FedData.png)](https://travis-ci.org/ropensci/FedData)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/FedData/master.svg)](https://codecov.io/github/ropensci/FedData?branch=master)
[![Zenodo
DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596344.svg)](https://doi.org/10.5281/zenodo.596344)
[![ROpenSci
Status](https://badges.ropensci.org/13_status.svg)](https://github.com/ropensci/software-review/issues/13)

**FedData version 3.0 is about to be released to CRAN. There are several
breaking changes in the FedData API from version 2.x. Please see
\[NEWS.md\] for a list of changes.**

`FedData` is an *R* package implementing functions to automate
downloading geospatial data available from several federated data
sources.

Currently, the package enables extraction from seven datasets:

-   The [National Elevation Dataset (NED)](https://ned.usgs.gov) digital
    elevation models (1 and 1/3 arc-second; USGS)
-   The [National Hydrography Dataset (NHD)](https://nhd.usgs.gov)
    (USGS)
-   The [Soil Survey Geographic (SSURGO)
    database](https://websoilsurvey.sc.egov.usda.gov/) from the National
    Cooperative Soil Survey (NCSS), which is led by the Natural
    Resources Conservation Service (NRCS) under the USDA
-   The [Global Historical Climatology Network
    (GHCN)](https://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn),
    coordinated by National Climatic Data Center at NOAA
-   The [Daymet](https://daymet.ornl.gov/) gridded estimates of daily
    weather parameters for North America, version 4, available from the
    Oak Ridge National Laboratory’s Distributed Active Archive Center
    (DAAC)
-   The [International Tree Ring Data Bank
    (ITRDB)](https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring),
    coordinated by National Climatic Data Center at NOAA
-   The [National Land Cover Database (NLCD)](https://www.mrlc.gov/)
-   The [NASS Cropland Data
    Layer](https://www.nass.usda.gov/Research_and_Science/Cropland/SARS1a.php)
    from the National Agricultural Statistics Service

This package is designed with the large-scale geographic information
system (GIS) use-case in mind: cases where the use of dynamic
web-services is impractical due to the scale (spatial and/or temporal)
of analysis. It functions primarily as a means of downloading tiled or
otherwise spatially-defined datasets; additionally, it can preprocess
those datasets by extracting data within an area of interest (AoI),
defined spatially. It relies heavily on the
[**sf**](https://cran.r-project.org/package=sf) and
[**raster**](https://cran.r-project.org/package=raster) packages.

This package has been built and tested on a binary install of *R* on
macOS 11.5 (Big Sur), and has been successfully run on Ubuntu via
[rocker/geospatial](https://hub.docker.com/r/rocker/geospatial) and on
Windows 10.

### Development

-   [Kyle Bocinsky](https://www.bocinsky.io) - Montana Climate Office,
    Missoula, MT

### Contributors

-   Dylan Beaudette - USDA-NRCS Soil Survey Office, Sonora, CA
-   Jeffrey Hollister - US EPA Atlantic Ecology Division, Narragansett,
    RI
-   Scott Chamberlain - ROpenSci and Museum of Paleontology at UC
    Berkeley

### Install `FedData`

-   From CRAN:

``` r
install.packages("FedData")
```

-   Development version from GitHub:

``` r
install.packages("devtools")
devtools::install_github("ropensci/FedData")
```

-   Linux: Follow instructions for installing `sf` available at
    <https://r-spatial.github.io/sf/>.

### Demonstration

This demonstration script is available as an R Markdown document in the
GitHub repository: <https://github.com/ropensci/FedData>.

#### Load `FedData` and define a study area

``` r
# FedData Tester
library(FedData)
library(magrittr)

# FedData comes loaded with the boundary of Mesa Verde National Park, for testing
FedData::meve
```

#### Get and plot the National Elevation Dataset for the study area

``` r
# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(
  template = FedData::meve,
  label = "meve"
)
# Plot with raster::plot
raster::plot(NED)
```

<img src="man/figures/README-NED-1.png" width="100%" />

#### Get and plot the Daymet dataset for the study area

``` r
# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp", "tmax"),
  years = 1980:1985
)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)
```

<img src="man/figures/README-DAYMET-1.png" width="100%" />

#### Get and plot the daily GHCN precipitation data for the study area

``` r
# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp")
)
#> Warning in if (!is.null(template) & !(class(template) %in%
#> c("SpatialPolygonsDataFrame", : the condition has length > 1 and only the first
#> element will be used
#> Warning: `select_()` was deprecated in dplyr 0.7.0.
#> Please use `select()` instead.
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.prcp$spatial,
  pch = 1,
  add = TRUE
)
legend("bottomleft",
  pch = 1,
  legend = "GHCN Precipitation Records"
)
```

<img src="man/figures/README-GHCN-precipitation-1.png" width="100%" />

#### Get and plot the daily GHCN temperature data for the study area

``` r
# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(
  template = FedData::meve,
  label = "meve",
  elements = c("tmin", "tmax"),
  years = 1980:1985,
  standardize = TRUE
)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.temp$spatial,
  add = TRUE,
  pch = 1
)
legend("bottomleft",
  pch = 1,
  legend = "GHCN Temperature Records"
)
```

<img src="man/figures/README-GHCN-temperature-1.png" width="100%" />

#### Get and plot the National Hydrography Dataset for the study area

``` r
# Get the NHD (USA ONLY)
get_nhd(
  template = FedData::meve,
  label = "meve"
) %>%
  plot_nhd(template = FedData::meve)
```

<img src="man/figures/README-NHD-1.png" width="100%" />

#### Get and plot the NRCS SSURGO data for the study area

``` r
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.MEVE <- get_ssurgo(
  template = FedData::meve,
  label = "meve"
)
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.MEVE$spatial$geom,
  lwd = 0.1,
  add = TRUE
)
```

<img src="man/figures/README-SSURGO-1.png" width="100%" />

#### Get and plot the NRCS SSURGO data for particular soil survey areas

``` r
# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(
  template = c("CO670", "CO075"),
  label = "CO_TEST"
)
#> Warning: One or more parsing issues, see `problems()` for details

#> Warning: One or more parsing issues, see `problems()` for details

#> Warning: One or more parsing issues, see `problems()` for details

#> Warning: One or more parsing issues, see `problems()` for details

#> Warning: One or more parsing issues, see `problems()` for details

#> Warning: One or more parsing issues, see `problems()` for details

#> Warning: One or more parsing issues, see `problems()` for details

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <-
  SSURGO.areas$spatial %>%
  dplyr::filter(AREASYMBOL == "CO075")

# And get the NED data under them for pretty plotting
NED.CO675 <- get_ned(
  template = SSURGO.areas.CO675,
  label = "SSURGO_CO675"
)

# Plot the SSURGO mapunit polygons, but only for CO675
raster::plot(NED.CO675)
plot(SSURGO.areas.CO675$geom,
  lwd = 0.1,
  add = TRUE
)
```

<img src="man/figures/README-SSURGO-area-1.png" width="100%" />

#### Get and plot the ITRDB chronology locations in the study area

``` r
# Get the ITRDB records
# Buffer MEVE, because there aren't any chronologies in the Park
ITRDB <- get_itrdb(
  template = FedData::meve %>%
    sf::st_buffer(50000),
  label = "meve",
  measurement.type = "Ring Width",
  chronology.type = "Standard"
)
#> Warning in eval(jsub, SDenv, parent.frame()): NAs introduced by coercion
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries

# Plot the MEVE buffer
plot(
  FedData::meve %>%
    sf::st_buffer(50000) %>%
    sf::st_transform(4326)
)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata$geometry,
  pch = 1,
  add = TRUE
)
legend("bottomleft",
  pch = 1,
  legend = "ITRDB chronologies"
)
```

<img src="man/figures/README-ITRDB-1.png" width="100%" />

#### Get and plot the National Land Cover Dataset for the study area

``` r
# Get the NLCD (USA ONLY)
# Returns a raster
NLCD <- get_nlcd(
  template = FedData::meve,
  year = 2011,
  label = "meve"
)

# Plot with raster::plot
raster::plot(NLCD)
```

<img src="man/figures/README-NLCD-1.png" width="100%" /><img src="man/figures/README-NLCD-2.png" width="100%" />

#### Get and plot the NASS Cropland Data Layer for the study area

``` r
# Get the NASS (USA ONLY)
# Returns a raster
NASS_CDL <- get_nass_cdl(
  template = FedData::meve,
  year = 2016,
  label = "meve"
)
# Plot with raster::plot
raster::plot(NASS_CDL)
```

<img src="man/figures/README-NASS-CDL-1.png" width="100%" /><img src="man/figures/README-NASS-CDL-2.png" width="100%" />

``` r
# Get the NASS CDL classification table
raster::levels(NASS_CDL)[[1]]

# Also, a convenience function loading the NASS CDL categories and hex colors
cdl_colors()
```

------------------------------------------------------------------------

### Acknowledgements

This package is a product of SKOPE ([Synthesizing Knowledge of Past
Environments](https://www.openskope.org/)) and the [Village Ecodynamics
Project](https://crowcanyon.github.io/veparchaeology/) through grants
awarded to the [Crow Canyon Archaeological
Center](https://www.crowcanyon.org) and Washington State University by
the National Science Foundation. This software is licensed under the
[MIT license](https://opensource.org/licenses/MIT). Continuing
development is supported by the [Montana Climate
Office](https://climate.umt.edu).

FedData was reviewed for [rOpenSci](https://ropensci.org) by
[@jooolia](https://github.com/jooolia), and was greatly improved as a
result. [rOpenSci](https://ropensci.org) on-boarding was coordinated by
[@sckott](https://github.com/sckott).

<!-- [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org) -->
# FedData 3.0.0.9000
* All `get_*()` functions now return `sf` or `raster` objects.
* Changed `SDA_query()` to `soils_query()` to avoid namespace masking with SoilDB.
* Added `get_nass_cdl()` to retrieve the NASS Cropland Data Layer
* Updated `get_daymet()` to pull from ORNL WCS, and fixed bug in [Issue #49](https://github.com/ropensci/FedData/issues/49)
* Fixed issue where `soils_query()` was only returning first SSURGO study area ([Issue #51](https://github.com/ropensci/FedData/issues/51))
* Fixed issue where date parsing was USA-specific by leveraging `lubridate::parse_date_time` ([Issue #61](https://github.com/ropensci/FedData/issues/61))
* `get_ssurgo()` now saves in the [GeoPackage file format](http://www.geopackage.org)
* `get_nhd()` now access ESRI web services. Users can optionally get data from NHDPlus.
* `get_nlcd()` now provides data in native CRS (CONUS Albers), rather than web-mercator ([Issue #77](https://github.com/ropensci/FedData/issues/77)) by pulling from self-hosted Cloud-Optimized GeoTiffs.
* Added Mesa Verde National Park as exemplar region, and removed use of `paleocar` package.
* `get_ned()` now pulls from USGS NED Cloud-Optimized GeoTiffs available at [https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/](https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/).

# FedData 2.5.7
* Removing many internet resource tests from CRAN, to satisfy: 'Packages which use Internet resources should fail gracefully with an informative message if the resource is not available (and not give a check warning nor error).'

# FedData 2.5.6
* Built-in access to the Soils Data Analysis query service to remove dependency on
soilDB package.

# FedData 2.5.5
* Fixed issue (#41) that occurs when mosaicking NLCD tiles that are not cropped. 
When they aren't cropped, the NLCD data is never read into memory, and the temporary 
file that the raster was created from gets destroyed.
Solution: Force NLCD data into memory prior to mosaicking.
* Added (non-CRAN) test for issue #41

# FedData 2.5.4
* Fixed issue in downloading NED tiles.

# FedData 2.5.3
* Added httr to package imports.

# FedData 2.5.2
* Updated NHD HUC4 to copy stored on Github.
* Fixed bug in ITRDB that caused some chronologies not to be read.

# FedData 2.5.1
* Switch to laze-loading data.
* Updated NHD paths to new National Map directory structure.

# FedData 2.5.0
* Added functions for the National Land Cover Database.

# FedData 2.4.7
* SSURGO fixed test where supplying an unavailable survey area now returns NULL instead of an error.
* SSURGO zip directory encoding changes as of late October 2017 forced changes in the FedData:::get_ssurgo_study_area function.
* Fixed issue where NHD template wouldn't load because they added a jpeg preview to the directory.

# FedData 2.4.6
* DAYMET functions now do *not* operate in parallel. This was breaking the download functions.
* Final update for version 2 of FedData.
* Accepted to ROpenSci! Migrating to the ROpenSci organization on GitHub.

# FedData 2.4.3
* writeOGR for SSURGO and NHD were failing on Windows when the `extraction.dir` included a trailing slash. Paths are now normalized to remove the trailing slash.

# FedData 2.4.2
* Updated the `get_ned` function to provide more useful errors and warnings when downloads are unsuccessful.

# FedData 2.4.1
* Added pkgdown site.
* SSURGO functions (e.g., `get_ssurgo`) now doesn't bomb on large (> 1 billion sq meter) requests. Now, the area of interest is broken into smaller chunks to build the download list.

# FedData 2.4.0
* Added a `NEWS.md` file to track changes to the package.
* Updated DAYMET functions to fix a bug that downloaded only one tile at a time.
* Linted all code.



## Note to CRAN maintainer

## Test environments
* local macOS install, R 4.1.0  (2021-05-18)
* ubuntu 14.04.5 (on Docker), R 3.5.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.---
name: FedData-release
layout: post
title: FedData - Getting assorted geospatial data into R
date: 2017-11-06
authors:
  - name: Kyle Bocinsky
categories:
  - technotes
tags:
  - R
  - package
  - data-access
  - packages
  - spatial
  - geospatial
  - review
  - onboarding
  - FedData
---

The package [FedData](https://github.com/ropensci/FedData) is now part of [rOpenSci](https://ropensci.org/). FedData includes functions to automate downloading geospatial data available from several federated data sources (mainly sources maintained by the US Federal government).

Currently, the package enables extraction from six datasets:

-   The [National Elevation Dataset (NED)](http://ned.usgs.gov) digital elevation models (1 and 1/3 arc-second; USGS)
-   The [National Hydrography Dataset (NHD)](http://nhd.usgs.gov) (USGS)
-   The [Soil Survey Geographic (SSURGO) database](http://websoilsurvey.sc.egov.usda.gov/) from the National Cooperative Soil Survey (NCSS), which is led by the Natural Resources Conservation Service (NRCS) under the USDA,
-   The [Global Historical Climatology Network (GHCN)](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn), coordinated by National Climatic Data Center at NOAA,
-   The [Daymet](https://daymet.ornl.gov/) gridded estimates of daily weather parameters for North America, version 3, available from the Oak Ridge National Laboratory's Distributed Active Archive Center (DAAC), and
-   The [International Tree Ring Data Bank (ITRDB)](http://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring), coordinated by National Climatic Data Center at NOAA.

FedData is designed with the large-scale geographic information system (GIS) use-case in mind: cases where the use of dynamic web-services is impractical due to the scale (spatial and/or temporal) of analysis. It functions primarily as a means of downloading tiled or otherwise spatially-defined datasets; additionally, it can preprocess those datasets by extracting data within an area of interest (AoI), defined spatially. It relies heavily on the [**sp**](https://cran.r-project.org/package=sp), [**raster**](https://cran.r-project.org/package=raster), and [**rgdal**](https://cran.r-project.org/package=rgdal) packages.

Acknowledgements
----------------

FedData is a product of SKOPE ([Synthesizing Knowledge of Past Environments](http://www.openskope.org)) and the [Village Ecodynamics Project](http://veparchaeology.org/).

FedData was reviewed for [rOpenSci](https://ropensci.org) by [@jooolia](https://github.com/jooolia), and was greatly improved as a result. [rOpenSci](https://ropensci.org) onboarding was coordinated by [@sckott](https://github.com/sckott).

TODO
----

**The current CRAN version of FedData, v2.4, will be the final minor CRAN release of FedData 2. FedData 3 will be released in the coming months, but some code built on FedData 2 will not be compatible with FedData 3.**

FedData was initially developed prior to widespread use of modern web mapping services and RESTful APIs by many Federal data-holders. Future releases of FedData will limit data transfer by utilizing server-side geospatial and data queries. We will also implement of data grammars from [dplyr](https://github.com/hadley/dplyr), tidy data structures, piping throughout ([magrittr](https://github.com/tidyverse/magrittr)), functional programming using [purrr](https://github.com/hadley/purrr), simple features for spatial data from [sf](https://github.com/edzer/sfr), and local data storage in OGC-compliant data formats (probably geojson and netCDF). I am also aiming for 100% testing coverage!

All that being said, much of the functionality of the FedData package could be spun off into more domain-specific packages. For example, ITRDB download functions could be part of the [dplR](https://r-forge.r-project.org/projects/dplr/) dendrochronology package; concepts/functions having to do with the GHCN data integrated into [rnoaa](https://github.com/ropensci/rnoaa); and Daymet concepts integrated into [daymetr](https://github.com/khufkens/daymetr). I welcome any and all suggestions about how to improve the utility of FedData; please [submit an issue](https://github.com/ropensci/FedData/issues).

Examples
--------

### Load `FedData` and define a study area

``` r
# FedData Tester
library(FedData)
library(magrittr)

# Extract data for the Village Ecodynamics Project "VEPIIN" study area:
# http://veparchaeology.org
vepPolygon <- polygon_from_extent(raster::extent(672800, 740000, 4102000, 4170000),
                                  proj4string = "+proj=utm +datum=NAD83 +zone=12")
```

### Get and plot the National Elevation Dataset for the study area

``` r
# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(template = vepPolygon,
               label = "VEPIIN")
# Plot with raster::plot
raster::plot(NED)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-6-1.png)

### Get and plot the Daymet dataset for the study area

``` r
# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(template = vepPolygon,
               label = "VEPIIN",
               elements = c("prcp","tmax"),
               years = 1980:1985)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-7-1.png)

### Get and plot the daily GHCN precipitation data for the study area

``` r
# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(template = vepPolygon, 
                            label = "VEPIIN", 
                            elements = c('prcp'))
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.prcp$spatial,
         pch = 1,
         add = TRUE)
legend('bottomleft',
       pch = 1,
       legend="GHCN Precipitation Records")
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-8-1.png)

### Get and plot the daily GHCN temperature data for the study area

``` r
# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(template = vepPolygon, 
                            label = "VEPIIN", 
                            elements = c('tmin','tmax'), 
                            years = 1980:1985,
                            standardize = TRUE)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.temp$spatial,
         add = TRUE,
         pch = 1)
legend('bottomleft',
       pch = 1,
       legend = "GHCN Temperature Records")
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-9-1.png)

### Get and plot the National Hydrography Dataset for the study area

``` r
# Get the NHD (USA ONLY)
NHD <- get_nhd(template = vepPolygon, 
               label = "VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the NHD data
NHD %>%
  lapply(sp::plot,
         col = 'black',
         add = TRUE)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-10-1.png)

### Get and plot the NRCS SSURGO data for the study area

``` r
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.VEPIIN <- get_ssurgo(template = vepPolygon, 
                     label = "VEPIIN")
#> Warning: 1 parsing failure.
#> row # A tibble: 1 x 5 col     row     col               expected actual expected   <int>   <chr>                  <chr>  <chr> actual 1  1276 slope.r no trailing characters     .5 file # ... with 1 more variables: file <chr>
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.VEPIIN$spatial,
     lwd = 0.1,
     add = TRUE)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-11-1.png)

### Get and plot the NRCS SSURGO data for particular soil survey areas

``` r
# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(template = c("CO670","CO075"), 
                           label = "CO_TEST")

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <- SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL=="CO075",]

# And get the NED data under them for pretty plotting
NED.CO675 <- get_ned(template = SSURGO.areas.CO675,
                            label = "SSURGO_CO675")
               
# Plot the SSURGO mapunit polygons, but only for CO675
plot(NED.CO675)
plot(SSURGO.areas.CO675,
     lwd = 0.1,
     add = TRUE)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-12-1.png)

### Get and plot the ITRDB chronology locations in the study area

``` r
# Get the ITRDB records
ITRDB <- get_itrdb(template = vepPolygon,
                        label = "VEPIIN",
                        makeSpatial = TRUE)
# Plot the NED again
raster::plot(NED)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata,
     pch = 1,
     add = TRUE)
legend('bottomleft',
       pch = 1,
       legend = "ITRDB chronologies")
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-13-1.png)
---
output: github_document
editor_options: 
  chunk_output_type: console
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%",
  dpi = 72
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN version](https://www.r-pkg.org/badges/version/FedData)](https://cran.r-project.org/package=FedData)
[![CRAN downloads per month](https://cranlogs.r-pkg.org/badges/FedData)](https://github.com/r-hub/cranlogs.app)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/FedData)](https://github.com/r-hub/cranlogs.app)
[![Build Status](https://api.travis-ci.org/ropensci/FedData.png)](https://travis-ci.org/ropensci/FedData)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/FedData/master.svg)](https://codecov.io/github/ropensci/FedData?branch=master)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596344.svg)](https://doi.org/10.5281/zenodo.596344)
[![ROpenSci Status](https://badges.ropensci.org/13_status.svg)](https://github.com/ropensci/software-review/issues/13)

**FedData version 3.0 is about to be released to CRAN. There are several breaking changes in the FedData API from version 2.x. Please see [NEWS.md] for a list of changes.**

`FedData` is an *R* package implementing functions to automate downloading geospatial data available from several federated data sources.

Currently, the package enables extraction from seven datasets: 

* The [National Elevation Dataset (NED)](https://ned.usgs.gov) digital elevation models (1 and 1/3 arc-second; USGS)
* The [National Hydrography Dataset (NHD)](https://nhd.usgs.gov) (USGS)
* The [Soil Survey Geographic (SSURGO) database](https://websoilsurvey.sc.egov.usda.gov/) from the National Cooperative Soil Survey (NCSS), which is led by the Natural Resources Conservation Service (NRCS) under the USDA
* The [Global Historical Climatology Network (GHCN)](https://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn), coordinated by National Climatic Data Center at NOAA
* The [Daymet](https://daymet.ornl.gov/) gridded estimates of daily weather parameters for North America, version 4, available from the Oak Ridge National Laboratory's Distributed Active Archive Center (DAAC)
* The [International Tree Ring Data Bank (ITRDB)](https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring), coordinated by National Climatic Data Center at NOAA
* The [National Land Cover Database (NLCD)](https://www.mrlc.gov/)
* The [NASS Cropland Data Layer](https://www.nass.usda.gov/Research_and_Science/Cropland/SARS1a.php) from the National Agricultural Statistics Service

This package is designed with the large-scale geographic information system (GIS) use-case in mind: cases where the use of dynamic web-services is impractical due to the scale (spatial and/or temporal) of analysis. It functions primarily as a means of downloading tiled or otherwise spatially-defined datasets; additionally, it can preprocess those datasets by extracting data within an area of interest (AoI), defined spatially. It relies heavily on the [**sf**](https://cran.r-project.org/package=sf) and [**raster**](https://cran.r-project.org/package=raster) packages.

This package has been built and tested on a binary install of *R* on macOS 11.5 (Big Sur), and has been successfully run on Ubuntu via [rocker/geospatial](https://hub.docker.com/r/rocker/geospatial) and on Windows 10.

### Development
+ [Kyle Bocinsky](https://www.bocinsky.io) - Montana Climate Office, Missoula, MT

### Contributors
+ Dylan Beaudette - USDA-NRCS Soil Survey Office, Sonora, CA
+ Jeffrey Hollister - US EPA Atlantic Ecology Division, Narragansett, RI
+ Scott Chamberlain - ROpenSci and Museum of Paleontology at UC Berkeley

### Install `FedData`
+ From CRAN:
```{r install, eval = FALSE}
install.packages("FedData")
```

+ Development version from GitHub:
```{r install dev, eval = FALSE}
install.packages("devtools")
devtools::install_github("ropensci/FedData")
```
+ Linux:
Follow instructions for installing `sf` available at https://r-spatial.github.io/sf/.

### Demonstration
This demonstration script is available as an R Markdown document in the GitHub repository: [https://github.com/ropensci/FedData](https://github.com/ropensci/FedData).

#### Load `FedData` and define a study area
```{r load, results="hide", message=FALSE}
# FedData Tester
library(FedData)
library(magrittr)

# FedData comes loaded with the boundary of Mesa Verde National Park, for testing
FedData::meve
```

#### Get and plot the National Elevation Dataset for the study area
```{r NED, results="hide", message=FALSE}
# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(
  template = FedData::meve,
  label = "meve"
)
# Plot with raster::plot
raster::plot(NED)
```

#### Get and plot the Daymet dataset for the study area
```{r DAYMET, results="hide", message=FALSE, warning=FALSE}
# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp", "tmax"),
  years = 1980:1985
)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)
```

#### Get and plot the daily GHCN precipitation data for the study area
```{r GHCN-precipitation, results="hide", message=FALSE}
# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp")
)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.prcp$spatial,
  pch = 1,
  add = TRUE
)
legend("bottomleft",
  pch = 1,
  legend = "GHCN Precipitation Records"
)
```

#### Get and plot the daily GHCN temperature data for the study area
```{r GHCN-temperature, results="hide", message=FALSE}
# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(
  template = FedData::meve,
  label = "meve",
  elements = c("tmin", "tmax"),
  years = 1980:1985,
  standardize = TRUE
)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.temp$spatial,
  add = TRUE,
  pch = 1
)
legend("bottomleft",
  pch = 1,
  legend = "GHCN Temperature Records"
)
```

#### Get and plot the National Hydrography Dataset for the study area
```{r NHD, results="hide", message=FALSE}
# Get the NHD (USA ONLY)
get_nhd(
  template = FedData::meve,
  label = "meve"
) %>%
  plot_nhd(template = FedData::meve)
```

#### Get and plot the NRCS SSURGO data for the study area
```{r SSURGO, results="hide", message=FALSE, warning=FALSE}
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.MEVE <- get_ssurgo(
  template = FedData::meve,
  label = "meve"
)
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.MEVE$spatial$geom,
  lwd = 0.1,
  add = TRUE
)
```

#### Get and plot the NRCS SSURGO data for particular soil survey areas
```{r SSURGO-area, results="hide", message=FALSE}
# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(
  template = c("CO670", "CO075"),
  label = "CO_TEST"
)

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <-
  SSURGO.areas$spatial %>%
  dplyr::filter(AREASYMBOL == "CO075")

# And get the NED data under them for pretty plotting
NED.CO675 <- get_ned(
  template = SSURGO.areas.CO675,
  label = "SSURGO_CO675"
)

# Plot the SSURGO mapunit polygons, but only for CO675
raster::plot(NED.CO675)
plot(SSURGO.areas.CO675$geom,
  lwd = 0.1,
  add = TRUE
)
```

#### Get and plot the ITRDB chronology locations in the study area
```{r ITRDB, results="hide", message=FALSE}
# Get the ITRDB records
# Buffer MEVE, because there aren't any chronologies in the Park
ITRDB <- get_itrdb(
  template = FedData::meve %>%
    sf::st_buffer(50000),
  label = "meve",
  measurement.type = "Ring Width",
  chronology.type = "Standard"
)

# Plot the MEVE buffer
plot(
  FedData::meve %>%
    sf::st_buffer(50000) %>%
    sf::st_transform(4326)
)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata$geometry,
  pch = 1,
  add = TRUE
)
legend("bottomleft",
  pch = 1,
  legend = "ITRDB chronologies"
)
```

#### Get and plot the National Land Cover Dataset for the study area
```{r NLCD, results="hide", message=FALSE}
# Get the NLCD (USA ONLY)
# Returns a raster
NLCD <- get_nlcd(
  template = FedData::meve,
  year = 2011,
  label = "meve"
)

# Plot with raster::plot
raster::plot(NLCD)
```

#### Get and plot the NASS Cropland Data Layer for the study area
```{r NASS-CDL, results="hide", message=FALSE}
# Get the NASS (USA ONLY)
# Returns a raster
NASS_CDL <- get_nass_cdl(
  template = FedData::meve,
  year = 2016,
  label = "meve"
)
# Plot with raster::plot
raster::plot(NASS_CDL)

# Get the NASS CDL classification table
raster::levels(NASS_CDL)[[1]]

# Also, a convenience function loading the NASS CDL categories and hex colors
cdl_colors()
```

-----------

### Acknowledgements
This package is a product of SKOPE ([Synthesizing Knowledge of Past Environments](https://www.openskope.org/)) and the [Village Ecodynamics Project](https://crowcanyon.github.io/veparchaeology/) through grants awarded to the [Crow Canyon Archaeological Center](https://www.crowcanyon.org) and Washington State University by the National Science Foundation. This software is licensed under the [MIT license](https://opensource.org/licenses/MIT). Continuing development is supported by the [Montana Climate Office](https://climate.umt.edu).

FedData was reviewed for [rOpenSci](https://ropensci.org) by [\@jooolia](https://github.com/jooolia), and was greatly improved as a result. [rOpenSci](https://ropensci.org) on-boarding was coordinated by [\@sckott](https://github.com/sckott).

<!-- [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org) -->
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GHCN_FUNCTIONS.R
\name{get_ghcn_daily_station}
\alias{get_ghcn_daily_station}
\title{Download and extract the daily data for a GHCN weather station.}
\usage{
get_ghcn_daily_station(
  ID,
  elements = NULL,
  years = NULL,
  raw.dir,
  standardize = F,
  force.redo = F
)
}
\arguments{
\item{ID}{A character string giving the station ID.}

\item{elements}{A character vector of elements to extract.\cr
The five core elements are:\cr
PRCP = Precipitation (tenths of mm)\cr
SNOW = Snowfall (mm)\cr
SNWD = Snow depth (mm)\cr
TMAX = Maximum temperature (tenths of degrees C)\cr
TMIN = Minimum temperature (tenths of degrees C)\cr
\cr
The other elements are:\cr

ACMC = Average cloudiness midnight to midnight from 30-second
ceilometer data (percent)\cr
ACMH = Average cloudiness midnight to midnight from
manual observations (percent)\cr
ACSC = Average cloudiness sunrise to sunset from 30-second
ceilometer data (percent)\cr
ACSH = Average cloudiness sunrise to sunset from manual
observations (percent)\cr
AWDR = Average daily wind direction (degrees)\cr
AWND = Average daily wind speed (tenths of meters per second)\cr
DAEV = Number of days included in the multiday evaporation
total (MDEV)\cr
DAPR = Number of days included in the multiday precipitation
total (MDPR)\cr
DASF = Number of days included in the multiday snowfall
total (MDSF)\cr
DATN = Number of days included in the multiday minimum temperature
(MDTN)\cr
DATX = Number of days included in the multiday maximum temperature
(MDTX)\cr
DAWM = Number of days included in the multiday wind movement
(MDWM)\cr
DWPR = Number of days with non-zero precipitation included in
multiday precipitation total (MDPR)\cr
EVAP = Evaporation of water from evaporation pan (tenths of mm)\cr
FMTM = Time of fastest mile or fastest 1-minute wind
(hours and minutes, i.e., HHMM)\cr
FRGB = Base of frozen ground layer (cm)\cr
FRGT = Top of frozen ground layer (cm)\cr
FRTH = Thickness of frozen ground layer (cm)\cr
GAHT = Difference between river and gauge height (cm)\cr
MDEV = Multiday evaporation total (tenths of mm; use with DAEV)\cr
MDPR = Multiday precipitation total (tenths of mm; use with DAPR and
                                     DWPR, if available)\cr
MDSF = Multiday snowfall total \cr
MDTN = Multiday minimum temperature (tenths of degrees C; use with DATN)\cr
MDTX = Multiday maximum temperature (tenths of degrees C; use with DATX)\cr
MDWM = Multiday wind movement (km)\cr
MNPN = Daily minimum temperature of water in an evaporation pan
(tenths of degrees C)\cr
MXPN = Daily maximum temperature of water in an evaporation pan
(tenths of degrees C)\cr
PGTM = Peak gust time (hours and minutes, i.e., HHMM)\cr
PSUN = Daily percent of possible sunshine (percent)\cr
SN*# = Minimum soil temperature (tenths of degrees C)
  where * corresponds to a code
for ground cover and # corresponds to a code for soil
depth.\cr
\cr
Ground cover codes include the following:\cr
0 = unknown\cr
1 = grass\cr
2 = fallow\cr
3 = bare ground\cr
4 = brome grass\cr
5 = sod\cr
6 = straw multch\cr
7 = grass muck\cr
8 = bare muck\cr
\cr
Depth codes include the following:\cr
1 = 5 cm\cr
2 = 10 cm\cr
3 = 20 cm\cr
4 = 50 cm\cr
5 = 100 cm\cr
6 = 150 cm\cr
7 = 180 cm\cr
\cr
SX*# = Maximum soil temperature (tenths of degrees C)
  where * corresponds to a code for ground cover
and # corresponds to a code for soil depth.\cr
See SN*# for ground cover and depth codes. \cr
TAVG = Average temperature (tenths of degrees C)
[Note that TAVG from source 'S' corresponds
to an average for the period ending at
2400 UTC rather than local midnight]\cr
THIC = Thickness of ice on water (tenths of mm)\cr
TOBS = Temperature at the time of observation (tenths of degrees C)\cr
TSUN = Daily total sunshine (minutes)\cr
WDF1 = Direction of fastest 1-minute wind (degrees)\cr
WDF2 = Direction of fastest 2-minute wind (degrees)\cr
WDF5 = Direction of fastest 5-second wind (degrees)\cr
WDFG = Direction of peak wind gust (degrees)\cr
WDFI = Direction of highest instantaneous wind (degrees)\cr
WDFM = Fastest mile wind direction (degrees)\cr
WDMV = 24-hour wind movement (km)\cr
WESD = Water equivalent of snow on the ground (tenths of mm)\cr
WESF = Water equivalent of snowfall (tenths of mm)\cr
WSF1 = Fastest 1-minute wind speed (tenths of meters per second)\cr
WSF2 = Fastest 2-minute wind speed (tenths of meters per second)\cr
WSF5 = Fastest 5-second wind speed (tenths of meters per second)\cr
WSFG = Peak gust wind speed (tenths of meters per second)\cr
WSFI = Highest instantaneous wind speed (tenths of meters per second)\cr
WSFM = Fastest mile wind speed (tenths of meters per second)\cr
WT** = Weather Type where ** has one of the following values:\cr
  \cr
01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
02 = Heavy fog or heaving freezing fog (not always
                                        distinguished from fog)\cr
03 = Thunder\cr
04 = Ice pellets, sleet, snow pellets, or small hail \cr
05 = Hail (may include small hail)\cr
06 = Glaze or rime \cr
07 = Dust, volcanic ash, blowing dust, blowing sand, or
blowing obstruction\cr
08 = Smoke or haze \cr
09 = Blowing or drifting snow\cr
10 = Tornado, waterspout, or funnel cloud \cr
11 = High or damaging winds\cr
12 = Blowing spray\cr
13 = Mist\cr
14 = Drizzle\cr
15 = Freezing drizzle \cr
16 = Rain (may include freezing rain, drizzle, and freezing drizzle) \cr
17 = Freezing rain \cr
18 = Snow, snow pellets, snow grains, or ice crystals\cr
19 = Unknown source of precipitation \cr
21 = Ground fog \cr
22 = Ice fog or freezing fog\cr
\cr
WV** = Weather in the Vicinity where ** has one of the following
values:\cr
01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
03 = Thunder\cr
07 = Ash, dust, sand, or other blowing obstruction\cr
18 = Snow or ice crystals\cr
20 = Rain or snow shower}

\item{years}{A numeric vector indicating which years to get.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.}

\item{standardize}{Select only common year/month/day? Defaults to FALSE.}

\item{force.redo}{If this weather station has been downloaded before, should it be updated? Defaults to FALSE.}
}
\value{
A named list of \code{\link{data.frame}s}, one for each \code{elements}.
}
\description{
\code{get_ghcn_daily_station} returns a named list of \code{\link{data.frame}s}, one for
each \code{elements}. If \code{elements} is undefined, it returns all available weather
tables for the station
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GHCN_FUNCTIONS.R
\name{download_ghcn_daily_station}
\alias{download_ghcn_daily_station}
\title{Download the daily data for a GHCN weather station.}
\usage{
download_ghcn_daily_station(ID, raw.dir, force.redo = F)
}
\arguments{
\item{ID}{A character string giving the station ID.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.}

\item{force.redo}{If this weather station has been downloaded before, should it be updated? Defaults to FALSE.}
}
\value{
A character string representing the full local path of the GHCN station data.
}
\description{
Download the daily data for a GHCN weather station.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GHCN_FUNCTIONS.R
\name{station_to_data_frame}
\alias{station_to_data_frame}
\title{Convert a list of station data to a single data frame.}
\usage{
station_to_data_frame(station.data)
}
\arguments{
\item{station.data}{A named list containing station data}
}
\value{
A \code{data.frame} of the containing the unwrapped station data
}
\description{
\code{station_to_data_frame} returns a \code{data.frame} of the GHCN station data list.
}
\details{
This function unwraps the station data and merges all data into a single data frame,
with the first column being in the \code{Date} class.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITRDB_FUNCTIONS.R
\name{get_itrdb}
\alias{get_itrdb}
\title{Download the latest version of the ITRDB, and extract given parameters.}
\usage{
get_itrdb(
  template = NULL,
  label = NULL,
  recon.years = NULL,
  calib.years = NULL,
  species = NULL,
  measurement.type = NULL,
  chronology.type = NULL,
  raw.dir = paste0(tempdir(), "/FedData/raw/itrdb"),
  extraction.dir = ifelse(!is.null(label), paste0(tempdir(),
    "/FedData/extractions/itrdb/", label, "/"), paste0(tempdir(),
    "/FedData/extractions/itrdb")),
  force.redo = FALSE
)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for selecting chronologies. If missing,
all available global chronologies are returned.}

\item{label}{A character string naming the study area.}

\item{recon.years}{A numeric vector of years over which reconstructions are needed;
if missing, the union of all years in the available chronologies are given.}

\item{calib.years}{A numeric vector of all required years---chronologies without these years will be discarded;
if missing, all available chronologies are given.}

\item{species}{A character vector of 4-letter tree species identifiers;
if missing, all available chronologies are given.}

\item{measurement.type}{A character vector of measurement type identifiers. Options include:
\itemize{
\item 'Total Ring Density'
\item 'Earlywood Width'
\item 'Earlywood Density'
\item 'Latewood Width'
\item 'Minimum Density'
\item 'Ring Width'
\item 'Latewood Density'
\item 'Maximum Density'
\item 'Latewood Percent'
}
if missing, all available chronologies are given.}

\item{chronology.type}{A character vector of chronology type identifiers. Options include:
\itemize{
\item 'ARSTND'
\item 'Low Pass Filter'
\item 'Residual'
\item 'Standard'
\item 'Re-Whitened Residual'
\item 'Measurements Only'
}
if missing, all available chronologies are given.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing.}

\item{extraction.dir}{A character string indicating where the extracted and cropped ITRDB dataset should be put.
The directory will be created if missing.}

\item{force.redo}{If an extraction already exists, should a new one be created? Defaults to FALSE.}
}
\value{
A named list containing the 'metadata', 'widths', and 'depths' data.
}
\description{
\code{get_itrdb} returns a named list of length 3:
\enumerate{
\item 'metadata': A data.table or \code{SpatialPointsDataFrame} (if \code{makeSpatial==TRUE}) of the locations
and names of extracted ITRDB chronologies,
\item 'widths': A matrix of tree-ring widths/densities given user selection, and
\item 'depths': A matrix of tree-ring sample depths.
}
}
\examples{
\dontrun{
# Get the ITRDB records
ITRDB <- get_itrdb(template = FedData::meve, label = "meve", makeSpatial = T)

# Plot the VEP polygon
plot(meve$geometry)

# Map the locations of the tree ring chronologies
plot(ITRDB$metadata, pch = 1, add = T)
legend("bottomleft", pch = 1, legend = "ITRDB chronologies")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITRDB_FUNCTIONS.R
\name{read_crn}
\alias{read_crn}
\title{Read a Tucson-format chronology file.}
\usage{
read_crn(file)
}
\arguments{
\item{file}{A character string path pointing to a \code{*.crn} file to be read.}
}
\value{
A list containing the metadata and chronology.
}
\description{
This function includes improvements to the \code{read.crn} function from the
\pkg{dplR} library. The principle changes are better parsing of metadata, and support
for the Schweingruber-type Tucson format. Chronologies that are unable to be read
are reported to the user. This function automatically recognizes Schweingruber-type files.
}
\details{
This wraps two other functions: \code{\link{read_crn_metadata}} \code{\link{read_crn_data}}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITRDB_FUNCTIONS.R
\name{download_itrdb}
\alias{download_itrdb}
\title{Download the latest version of the ITRDB.}
\usage{
download_itrdb(
  raw.dir = paste0(tempdir(), "/FedData/raw/itrdb"),
  force.redo = FALSE
)
}
\arguments{
\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing. Defaults to './RAW/ITRDB/'.}

\item{force.redo}{If a download already exists, should a new one be created? Defaults to FALSE.}
}
\value{
A data.table containing all of the ITRDB data.
}
\description{
Downloads and parses the latest zipped (numbered) version of the ITRDB.
This function includes improvements to the \code{\link{read_crn}} function from the
\pkg{dplR} library. The principle changes are better parsing of metadata, and support
for the Schweingruber-type Tucson format. Chronologies that are unable to be read
are reported to the user.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{polygon_from_extent}
\alias{polygon_from_extent}
\title{Turn an extent object into a polygon}
\usage{
polygon_from_extent(x, proj4string = NULL)
}
\arguments{
\item{x}{An \code{\link{extent}} object, or an object from which an extent object can be retrieved.}

\item{proj4string}{A PROJ.4 formatted string defining the required projection. If NULL,
the function will attempt to get the projection from x using \code{\link{projection}}}
}
\value{
A SpatialPolygons object.
}
\description{
Turn an extent object into a polygon
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GHCN_FUNCTIONS.R
\name{get_ghcn_daily}
\alias{get_ghcn_daily}
\title{Download and crop the Global Historical Climate Network-Daily data.}
\usage{
get_ghcn_daily(
  template = NULL,
  label = NULL,
  elements = NULL,
  years = NULL,
  raw.dir = paste0(tempdir(), "/FedData/raw/ghcn"),
  extraction.dir = paste0(tempdir(), "/FedData/extractions/ghcn/", label, "/"),
  standardize = F,
  force.redo = F
)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping. Alternatively, a character vector providing GHCN station IDs. If missing, all stations
will be downloaded!}

\item{label}{A character string naming the study area.}

\item{elements}{A character vector of elements to extract.\cr
The five core elements are:\cr
PRCP = Precipitation (tenths of mm)\cr
SNOW = Snowfall (mm)\cr
SNWD = Snow depth (mm)\cr
TMAX = Maximum temperature (tenths of degrees C)\cr
TMIN = Minimum temperature (tenths of degrees C)\cr
\cr
The other elements are:\cr

ACMC = Average cloudiness midnight to midnight from 30-second
ceilometer data (percent)\cr
ACMH = Average cloudiness midnight to midnight from
manual observations (percent)\cr
ACSC = Average cloudiness sunrise to sunset from 30-second
ceilometer data (percent)\cr
ACSH = Average cloudiness sunrise to sunset from manual
observations (percent)\cr
AWDR = Average daily wind direction (degrees)\cr
AWND = Average daily wind speed (tenths of meters per second)\cr
DAEV = Number of days included in the multiday evaporation
total (MDEV)\cr
DAPR = Number of days included in the multiday precipitation
total (MDPR)\cr
DASF = Number of days included in the multiday snowfall
total (MDSF)\cr
DATN = Number of days included in the multiday minimum temperature
(MDTN)\cr
DATX = Number of days included in the multiday maximum temperature
(MDTX)\cr
DAWM = Number of days included in the multiday wind movement
(MDWM)\cr
DWPR = Number of days with non-zero precipitation included in
multiday precipitation total (MDPR)\cr
EVAP = Evaporation of water from evaporation pan (tenths of mm)\cr
FMTM = Time of fastest mile or fastest 1-minute wind
(hours and minutes, i.e., HHMM)\cr
FRGB = Base of frozen ground layer (cm)\cr
FRGT = Top of frozen ground layer (cm)\cr
FRTH = Thickness of frozen ground layer (cm)\cr
GAHT = Difference between river and gauge height (cm)\cr
MDEV = Multiday evaporation total (tenths of mm; use with DAEV)\cr
MDPR = Multiday precipitation total (tenths of mm; use with DAPR and
                                     DWPR, if available)\cr
MDSF = Multiday snowfall total \cr
MDTN = Multiday minimum temperature (tenths of degrees C; use with DATN)\cr
MDTX = Multiday maximum temperature (tenths of degrees C; use with DATX)\cr
MDWM = Multiday wind movement (km)\cr
MNPN = Daily minimum temperature of water in an evaporation pan
(tenths of degrees C)\cr
MXPN = Daily maximum temperature of water in an evaporation pan
(tenths of degrees C)\cr
PGTM = Peak gust time (hours and minutes, i.e., HHMM)\cr
PSUN = Daily percent of possible sunshine (percent)\cr
SN*# = Minimum soil temperature (tenths of degrees C)
  where * corresponds to a code
for ground cover and # corresponds to a code for soil
depth.\cr
\cr
Ground cover codes include the following:\cr
0 = unknown\cr
1 = grass\cr
2 = fallow\cr
3 = bare ground\cr
4 = brome grass\cr
5 = sod\cr
6 = straw multch\cr
7 = grass muck\cr
8 = bare muck\cr
\cr
Depth codes include the following:\cr
1 = 5 cm\cr
2 = 10 cm\cr
3 = 20 cm\cr
4 = 50 cm\cr
5 = 100 cm\cr
6 = 150 cm\cr
7 = 180 cm\cr
\cr
SX*# = Maximum soil temperature (tenths of degrees C)
  where * corresponds to a code for ground cover
and # corresponds to a code for soil depth.\cr
See SN*# for ground cover and depth codes. \cr
TAVG = Average temperature (tenths of degrees C)
[Note that TAVG from source 'S' corresponds
to an average for the period ending at
2400 UTC rather than local midnight]\cr
THIC = Thickness of ice on water (tenths of mm)\cr
TOBS = Temperature at the time of observation (tenths of degrees C)\cr
TSUN = Daily total sunshine (minutes)\cr
WDF1 = Direction of fastest 1-minute wind (degrees)\cr
WDF2 = Direction of fastest 2-minute wind (degrees)\cr
WDF5 = Direction of fastest 5-second wind (degrees)\cr
WDFG = Direction of peak wind gust (degrees)\cr
WDFI = Direction of highest instantaneous wind (degrees)\cr
WDFM = Fastest mile wind direction (degrees)\cr
WDMV = 24-hour wind movement (km)\cr
WESD = Water equivalent of snow on the ground (tenths of mm)\cr
WESF = Water equivalent of snowfall (tenths of mm)\cr
WSF1 = Fastest 1-minute wind speed (tenths of meters per second)\cr
WSF2 = Fastest 2-minute wind speed (tenths of meters per second)\cr
WSF5 = Fastest 5-second wind speed (tenths of meters per second)\cr
WSFG = Peak gust wind speed (tenths of meters per second)\cr
WSFI = Highest instantaneous wind speed (tenths of meters per second)\cr
WSFM = Fastest mile wind speed (tenths of meters per second)\cr
WT** = Weather Type where ** has one of the following values:\cr
  \cr
01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
02 = Heavy fog or heaving freezing fog (not always
                                        distinguished from fog)\cr
03 = Thunder\cr
04 = Ice pellets, sleet, snow pellets, or small hail \cr
05 = Hail (may include small hail)\cr
06 = Glaze or rime \cr
07 = Dust, volcanic ash, blowing dust, blowing sand, or
blowing obstruction\cr
08 = Smoke or haze \cr
09 = Blowing or drifting snow\cr
10 = Tornado, waterspout, or funnel cloud \cr
11 = High or damaging winds\cr
12 = Blowing spray\cr
13 = Mist\cr
14 = Drizzle\cr
15 = Freezing drizzle \cr
16 = Rain (may include freezing rain, drizzle, and freezing drizzle) \cr
17 = Freezing rain \cr
18 = Snow, snow pellets, snow grains, or ice crystals\cr
19 = Unknown source of precipitation \cr
21 = Ground fog \cr
22 = Ice fog or freezing fog\cr
\cr
WV** = Weather in the Vicinity where ** has one of the following
values:\cr
01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
03 = Thunder\cr
07 = Ash, dust, sand, or other blowing obstruction\cr
18 = Snow or ice crystals\cr
20 = Rain or snow shower}

\item{years}{A numeric vector indicating which years to get.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing. Defaults to './RAW/GHCN/'.}

\item{extraction.dir}{A character string indicating where the extracted and cropped GHCN shapefiles should be put.
The directory will be created if missing. Defaults to './EXTRACTIONS/GHCN/'.}

\item{standardize}{Select only common year/month/day? Defaults to FALSE.}

\item{force.redo}{If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.}
}
\value{
A named list containing the 'spatial' and 'tabular' data.
}
\description{
\code{get_ghcn_daily} returns a named list of length 2:
\enumerate{
\item 'spatial': A \code{SpatialPointsDataFrame} of the locations of GHCN weather stations
in the template, and
\item 'tabular': A named list of \code{\link{data.frame}s} with the daily weather data for each station.
The name of each list item is the station ID.
}
}
\examples{
\dontrun{
# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <-
  get_ghcn_daily(
    template = FedData::meve,
    label = "meve",
    elements = c("prcp")
  )

# Plot the VEP polygon
plot(meve$geometry)

# Plot the spatial locations
plot(GHCN.prcp$spatial, pch = 1, add = T)
legend("bottomleft", pch = 1, legend = "GHCN Precipitation Records")

# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(
  template = FedData::meve,
  label = "meve",
  elements = c("tmin", "tmax"),
  standardize = T
)

# Plot the VEP polygon
plot(meve$geometry)

# Plot the spatial locations
plot(GHCN.temp$spatial, pch = 1, add = T)
legend("bottomleft", pch = 1, legend = "GHCN Temperature Records")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{split_bbox}
\alias{split_bbox}
\title{Splits a bbox into a list of bboxes less than a certain size}
\usage{
split_bbox(bbox, x, y = x)
}
\arguments{
\item{x}{The maximum x size of the resulting bounding boxes}

\item{y}{The maximum y size of the resulting bounding boxes; defaults to x}
}
\value{
A list of bbox objects
}
\description{
Splits a bbox into a list of bboxes less than a certain size
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{check_service}
\alias{check_service}
\title{Check whether a web service is unavailable, and stop function if necessary.}
\usage{
check_service(x)
}
\arguments{
\item{x}{The path to the web service.}
}
\value{
Error if service unavailable.
}
\description{
Check whether a web service is unavailable, and stop function if necessary.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NED_FUNCTIONS.R
\name{get_ned}
\alias{get_ned}
\title{Download and crop the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.}
\usage{
get_ned(
  template,
  label,
  res = "1",
  raw.dir = paste0(tempdir(), "/FedData/raw/ned"),
  extraction.dir = paste0(tempdir(), "/FedData/extractions/ned/", label, "/"),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9"),
  force.redo = F
)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping.}

\item{label}{A character string naming the study area.}

\item{res}{A character string representing the desired resolution of the NED. '1'
indicates the 1 arc-second NED (the default), while '13' indicates the 1/3 arc-second dataset.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing. Defaults to './RAW/NED/'.}

\item{extraction.dir}{A character string indicating where the extracted and cropped DEM should be put.
The directory will be created if missing. Defaults to './EXTRACTIONS/NED/'.}

\item{raster.options}{a vector of options for raster::writeRaster.}

\item{force.redo}{If an extraction for this template and label already exists, should a new one be created?}
}
\value{
A \code{RasterLayer} DEM cropped to the extent of the template.
}
\description{
\code{get_ned} returns a \code{RasterLayer} of elevation data cropped to a given
template study area.
}
\examples{
\dontrun{
# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(template = FedData::meve, label = "meve")

# Plot with raster::plot
plot(NED)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{get_ssurgo_inventory}
\alias{get_ssurgo_inventory}
\title{Download and crop a shapefile of the SSURGO study areas.}
\usage{
get_ssurgo_inventory(template = NULL, raw.dir)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing.}
}
\value{
A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
the specified \code{template}.
}
\description{
\code{get_ssurgo_inventory} returns a \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
the specified \code{template}. If template is not provided, returns the entire SSURGO inventory of study areas.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NHD_FUNCTIONS.R
\name{get_wbd}
\alias{get_wbd}
\title{Download and crop the Watershed Boundary Dataset.}
\usage{
get_wbd(
  template,
  label,
  extraction.dir = paste0(tempdir(), "/FedData/extractions/nhd/", label, "/"),
  force.redo = FALSE
)
}
\arguments{
\item{template}{A `Raster*`, `Spatial*`, or `sf` object to serve
as a template for cropping.}

\item{label}{A character string naming the study area.}

\item{extraction.dir}{A character string indicating where the extracted and cropped NHD data should be put.}

\item{force.redo}{If an extraction for this template and label already exists, should a new one be created?}
}
\value{
A  `sf` collection of the HUC 12 regions within
the specified \code{template}.
}
\description{
\code{get_wbd} returns an `sf` collection of the HUC 12 regions within
the specified \code{template}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{get_ssurgo_study_area}
\alias{get_ssurgo_study_area}
\title{Download and crop the spatial and tabular data for a SSURGO study area.}
\usage{
get_ssurgo_study_area(template = NULL, area, date, raw.dir)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping. If missing, whose study area is returned}

\item{area}{A character string indicating the SSURGO study area to be downloaded.}

\item{date}{A character string indicating the date of the most recent update to the SSURGO
area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{get_ssurgo_inventory}}).}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing.}
}
\value{
A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
the specified \code{template}.
}
\description{
\code{get_ssurgo_study_area} returns a named list of length 2:
\enumerate{
\item 'spatial': A \code{SpatialPolygonsDataFrame} of soil mapunits
in the template, and
\item 'tabular': A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{url_base}
\alias{url_base}
\title{Strip query parameters from a URL}
\usage{
url_base(x)
}
\arguments{
\item{url}{The URL to be modified}
}
\value{
The URL without parameters
}
\description{
Strip query parameters from a URL
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GHCN_FUNCTIONS.R
\name{get_ghcn_inventory}
\alias{get_ghcn_inventory}
\title{Download and crop the inventory of GHCN stations.}
\usage{
get_ghcn_inventory(template = NULL, elements = NULL, raw.dir)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping.}

\item{elements}{A character vector of elements to extract.
Common elements include 'tmin', 'tmax', and 'prcp'.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing.}
}
\value{
A \code{SpatialPolygonsDataFrame} of the GHCN stations within
the specified \code{template}
}
\description{
\code{get_ghcn_inventory} returns a \code{SpatialPolygonsDataFrame} of the GHCN stations within
the specified \code{template}. If template is not provided, returns the entire GHCN inventory.
}
\details{
Stations with multiple elements will have multiple points. This allows for easy mapping of stations
by element availability.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NASS_FUNCTIONS.R
\name{get_nass_cdl}
\alias{get_nass_cdl}
\alias{get_nass}
\alias{get_cdl}
\alias{cdl_colors}
\title{Download and crop the NASS Cropland Data Layer.}
\usage{
get_nass_cdl(
  template,
  label,
  year = 2019,
  extraction.dir = paste0(tempdir(), "/FedData/"),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
  force.redo = FALSE,
  progress = TRUE
)

get_nass(template, label, ...)

get_cdl(template, label, ...)

cdl_colors()
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping.}

\item{label}{A character string naming the study area.}

\item{year}{An integer representing the year of desired NASS Cropland Data Layer product.
Acceptable values are 2007--the last year.}

\item{extraction.dir}{A character string indicating where the extracted and cropped NASS data should be put.
The directory will be created if missing.}

\item{raster.options}{a vector of options for raster::writeRaster.}

\item{force.redo}{If an extraction for this template and label already exists, should a new one be created?}

\item{progress}{Draw a progress bar when downloading?}

\item{...}{Other parameters passed on to [get_nass_cdl].}
}
\value{
A \code{RasterLayer} cropped to the bounding box of the template.
}
\description{
\code{get_nass_cdl} returns a \code{RasterLayer} of NASS Cropland Data Layer cropped to a given
template study area.
}
\examples{
\dontrun{
# Extract data for the Mesa Verde National Park:

# Get the NASS CDL (USA ONLY)
# Returns a raster
NASS <-
  get_nass_cdl(
    template = FedData::meve,
    label = "meve",
    year = 2011
  )

# Plot with raster::plot
plot(NASS)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{substr_right}
\alias{substr_right}
\title{Get the rightmost 'n' characters of a character string.}
\usage{
substr_right(x, n)
}
\arguments{
\item{x}{A character string.}

\item{n}{The number of characters to retrieve.}
}
\value{
A character string.
}
\description{
Get the rightmost 'n' characters of a character string.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITRDB_FUNCTIONS.R
\name{read_crn_metadata}
\alias{read_crn_metadata}
\title{Read metadata from a Tucson-format chronology file.}
\usage{
read_crn_metadata(file, SCHWEINGRUBER)
}
\arguments{
\item{file}{A character string path pointing to a \code{*.crn} file to be read.}

\item{SCHWEINGRUBER}{Is the file in the Schweingruber-type Tucson format?}
}
\value{
A data.frame containing the metadata.
}
\description{
This function includes improvements to the \code{\link{read_crn}} function from the
\pkg{dplR} library. The principle changes are better parsing of metadata, and support
for the Schweingruber-type Tucson format. Chronologies that are unable to be read
are reported to the user. The user (or \code{\link{read_crn}}) must tell the function whether
the file is a Schweingruber-type chronology.
}
\details{
Location information is converted to decimal degrees.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{spdf_from_polygon}
\alias{spdf_from_polygon}
\title{Turn a SpatialPolygons object into a SpatialPolygonsDataFrame.}
\usage{
spdf_from_polygon(x)
}
\arguments{
\item{x}{An SpatialPolygons object.}
}
\value{
A SpatialPolygonsDataFrame object.
}
\description{
Turn a SpatialPolygons object into a SpatialPolygonsDataFrame.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{download_ssurgo_inventory}
\alias{download_ssurgo_inventory}
\title{Download a zipped directory containing a shapefile of the SSURGO study areas.}
\usage{
download_ssurgo_inventory(raw.dir, ...)
}
\arguments{
\item{raw.dir}{A character string indicating where raw downloaded files should be put.}
}
\value{
A character string representing the full local path of the SSURGO study areas zipped directory.
}
\description{
Download a zipped directory containing a shapefile of the SSURGO study areas.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAYMET_FUNCTIONS.R
\name{get_daymet}
\alias{get_daymet}
\title{Download and crop the 1-km DAYMET v4 daily weather dataset.}
\usage{
get_daymet(
  template,
  label,
  elements = c("dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp"),
  years = 1980:(lubridate::year(Sys.time()) - 1),
  region = "na",
  tempo = "day",
  extraction.dir = paste0(tempdir(), "/FedData/extractions/daymet/", label, "/"),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
  force.redo = F,
  progress = TRUE
)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping.}

\item{label}{A character string naming the study area.}

\item{elements}{A character vector of elements to extract.\cr
The available elements are:\cr
dayl = Duration of the daylight period in seconds per day. This calculation is based on the period of the day during which the sun is above a hypothetical flat horizon.\cr
prcp = Daily total precipitation in millimeters per day, sum of all forms converted to water-equivalent. Precipitation occurrence on any given day may be ascertained.\cr
srad = Incident shortwave radiation flux density in watts per square meter, taken as an average over the daylight period of the day. NOTE: Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)\cr
swe = Snow water equivalent in kilograms per square meter. The amount of water contained within the snowpack.\cr
tmax = Daily maximum 2-meter air temperature in degrees Celsius.\cr
tmin = Daily minimum 2-meter air temperature in degrees Celsius.\cr
vp = Water vapor pressure in pascals. Daily average partial pressure of water vapor.\cr}

\item{years}{A numeric vector of years to extract.}

\item{region}{The name of a region. The available regions are:\cr
na = North America\cr
hi = Hawaii\cr
pr = Puerto Rico\cr}

\item{tempo}{The frequency of the data. The available tempos are:\cr
day = Daily data\cr
mon = Monthly summary data\cr
ann = Annual summary data\cr}

\item{extraction.dir}{A character string indicating where the extracted and cropped DEM should be put.
Defaults to a temporary directory.}

\item{raster.options}{a vector of options for raster::writeRaster.}

\item{force.redo}{If an extraction for this template and label already exists in extraction.dir,
should a new one be created?}

\item{progress}{Draw a progress bar when downloading?}
}
\value{
A named list of \code{RasterBrick}s of weather data cropped to the extent of the template.
}
\description{
\code{get_daymet} returns a \code{RasterBrick} of weather data cropped to a given
template study area.
}
\examples{
\dontrun{

# Get the DAYMET (North America only)
# Returns a list of raster bricks
DAYMET <- get_daymet(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp", "tmin", "tmax"),
  years = 1980:1985
)

# Plot with raster::plot
plot(DAYMET$tmin$X1985.10.23)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{download_ssurgo_study_area}
\alias{download_ssurgo_study_area}
\title{Download a zipped directory containing the spatial and tabular data for a SSURGO study area.}
\usage{
download_ssurgo_study_area(area, date, raw.dir)
}
\arguments{
\item{area}{A character string indicating the SSURGO study area to be downloaded.}

\item{date}{A character string indicating the date of the most recent update to the SSURGO
area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{get_ssurgo_inventory}}).}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.}
}
\value{
A character string representing the full local path of the SSURGO study areas zipped directory.
}
\description{
\code{download_ssurgo_study_area} first tries to download data including a state-specific Access
template, then the general US template.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{download_data}
\alias{download_data}
\title{Use curl to download a file.}
\usage{
download_data(
  url,
  destdir = getwd(),
  timestamping = T,
  nc = F,
  verbose = F,
  progress = F
)
}
\arguments{
\item{url}{The location of a file.}

\item{destdir}{Where the file should be downloaded to.}

\item{timestamping}{Should only newer files be downloaded?}

\item{nc}{Should files of the same type not be clobbered?}

\item{verbose}{Should cURL output be shown?}

\item{progress}{Should a progress bar be shown with cURL output?}
}
\value{
A character string of the file path to the downloaded file.
}
\description{
This function makes it easy to implement timestamping and no-clobber of files.
}
\details{
If both \code{timestamping} and \code{nc} are TRUE, nc behavior trumps timestamping.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{meve}
\alias{meve}
\title{The boundary of Mesa Verde National Park}
\format{
Simple feature collection with 1 feature and a geometry field.
}
\usage{
meve
}
\description{
A dataset containing the spatial polygon defining the boundary
of Mesa Verde National Park in Montana.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{get_ssurgo}
\alias{get_ssurgo}
\title{Download and crop data from the NRCS SSURGO soils database.}
\usage{
get_ssurgo(
  template,
  label,
  raw.dir = paste0(tempdir(), "/FedData/raw/ssurgo"),
  extraction.dir = paste0(tempdir(), "/FedData/"),
  force.redo = FALSE
)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping; optionally, a vector of area names [e.g., c('IN087','IN088')] may be provided.}

\item{label}{A character string naming the study area.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing. Defaults to './RAW/SSURGO/'.}

\item{extraction.dir}{A character string indicating where the extracted and cropped SSURGO shapefiles should be put.
The directory will be created if missing. Defaults to './EXTRACTIONS/SSURGO/'.}

\item{force.redo}{If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.}
}
\value{
A named list containing the 'spatial' and 'tabular' data.
}
\description{
This is an efficient method for spatially merging several different soil survey areas
as well as merging their tabular data.
}
\details{
\code{get_ssurgo} returns a named list of length 2:
\enumerate{
\item 'spatial': A \code{SpatialPolygonsDataFrame} of soil mapunits
in the template, and
\item 'tabular': A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
}
}
\examples{
\dontrun{
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.MEVE <- get_ssurgo(template = FedData::meve, label = "meve")

# Plot the VEP polygon
plot(meve$geometry)

# Plot the SSURGO mapunit polygons
plot(SSURGO.MEVE$spatial, lwd = 0.1, add = T)

# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(template = c("CO670", "CO075"), label = "CO_TEST")

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <- SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL == "CO075", ]

# And get the NED data under them for pretty plotting
NED.CO675 <- get_ned(template = SSURGO.areas.CO675, label = "SSURGO_CO675")

# Plot the SSURGO mapunit polygons, but only for CO675
plot(NED.CO675)
plot(SSURGO.areas.CO675, lwd = 0.1, add = T)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{soils_query}
\alias{soils_query}
\title{Submit a Soil Data Access (SDA) Query}
\usage{
soils_query(q)
}
\arguments{
\item{q}{A character string representing a SQL query to the SDA service}
}
\value{
A tibble returned from the SDA service
}
\description{
\code{soils_query} submit an SQL query to retrieve data from the Soil Data Mart.
Please see https://sdmdataaccess.sc.egov.usda.gov/Query.aspx for guidelines
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{replace_null}
\alias{replace_null}
\title{Replace NULLs}
\usage{
replace_null(x)
}
\arguments{
\item{x}{A list}
}
\description{
Replace all the empty values in a list
}
\examples{
list(a = NULL, b = 1, c = list(foo = NULL, bar = NULL)) \%>\% replace_null()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NHD_FUNCTIONS.R
\name{plot_nhd}
\alias{plot_nhd}
\title{A basic plotting function for NHD data.}
\usage{
plot_nhd(x, template = NULL)
}
\arguments{
\item{x}{The result of [get_nhd].}

\item{template}{template A `Raster*`, `Spatial*`, or `sf` object to serve
as a template for cropping.}
}
\value{
A [ggplot2] panel of plots
}
\description{
This is more of an example than anything
}
\examples{
\dontrun{
# Get the NHD (USA ONLY)
NHD <- get_nhd(
  template = FedData::meve,
  label = "meve"
)
NHD
NHD \%>\%
  plot_nhd(template = FedData::meve)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAYMET_FUNCTIONS.R
\name{download_daymet_thredds}
\alias{download_daymet_thredds}
\title{Download the 1-km DAYMET daily weather dataset for a region as a netcdf.}
\usage{
download_daymet_thredds(bbox, element, year, region, tempo)
}
\arguments{
\item{bbox}{the bounding box in WGS84 coordinates as a comma-separated character vector
"xmin,ymin,xmax,ymax"}

\item{element}{An element to extract.\cr
The available elements are:\cr
dayl = Duration of the daylight period in seconds per day. This calculation is based on the period of the day during which the sun is above a hypothetical flat horizon.\cr
prcp = Daily total precipitation in millimeters per day, sum of all forms converted to water-equivalent. Precipitation occurrence on any given day may be ascertained.\cr
srad = Incident shortwave radiation flux density in watts per square meter, taken as an average over the daylight period of the day. NOTE: Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)\cr
swe = Snow water equivalent in kilograms per square meter. The amount of water contained within the snowpack.\cr
tmax = Daily maximum 2-meter air temperature in degrees Celsius.\cr
tmin = Daily minimum 2-meter air temperature in degrees Celsius.\cr
vp = Water vapor pressure in pascals. Daily average partial pressure of water vapor.\cr}

\item{year}{An integer year to extract.}

\item{region}{The name of a region. The available regions are:\cr
na = North America\cr
hi = Hawaii\cr
pr = Puerto Rico\cr}

\item{tempo}{The frequency of the data. The available tempos are:\cr
day = Daily data\cr
mon = Monthly summary data\cr
ann = Annual summary data\cr}
}
\value{
A named list of character vectors, each representing the full local paths of the tile downloads.
}
\description{
Data are downloaded in the NetCDF format. \code{download_daymet_thredds} returns the path to the downloaded NetCDF file.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NED_FUNCTIONS.R
\name{get_ned_tile}
\alias{get_ned_tile}
\title{Load and crop tile from the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.}
\usage{
get_ned_tile(template = NULL, res = "1", tileNorthing, tileWesting, raw.dir)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping. If missing, entire tile is returned.}

\item{res}{A character string representing the desired resolution of the NED. '1'
indicates the 1 arc-second NED (the default), while '13' indicates the 1/3 arc-second dataset.}

\item{tileNorthing}{An integer representing the northing (latitude, in degrees north of the equator) of the northwest corner of the tile to
be downloaded.}

\item{tileWesting}{An integer representing the westing (longitude, in degrees west of the prime meridian) of the northwest corner of the tile to
be downloaded.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing. Defaults to './RAW/NED/'.}
}
\value{
A \code{RasterLayer} cropped within the specified \code{template}.
}
\description{
\code{get_ned_tile} returns a \code{RasterLayer} cropped within the specified \code{template}.
If template is not provided, returns the entire NED tile.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NLCD_FUNCTIONS.R
\name{get_nlcd}
\alias{get_nlcd}
\alias{nlcd_colors}
\alias{pal_nlcd}
\title{Download and crop the National Land Cover Database.}
\usage{
get_nlcd(
  template,
  label,
  year = 2019,
  dataset = c("landcover", "impervious", "canopy"),
  landmass = "L48",
  extraction.dir = paste0(tempdir(), "/FedData/"),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9"),
  force.redo = F
)

nlcd_colors()

pal_nlcd()
}
\arguments{
\item{template}{A sf, Raster* or Spatial* object to serve
as a template for cropping.}

\item{label}{A character string naming the study area.}

\item{year}{An integer representing the year of desired NLCD product.
Acceptable values are 2019 (default), 2016, 2011, 2008, 2006, 2004, and 2001.}

\item{dataset}{A character string representing type of the NLCD product.
Acceptable values are 'landcover' (default), 'impervious', and
'canopy' (2016 and 2011, L48 only).}

\item{landmass}{A character string representing the landmass to be extracted
Acceptable values are 'L48' (lower 48 US states, the default),
'AK' (Alaska, 2011 and 2016 only), 'HI' (Hawaii, 2001 only), and
'PR' (Puerto Rico, 2001 only).}

\item{extraction.dir}{A character string indicating where the extracted
and cropped NLCD data should be put. The directory will be created if missing.}

\item{raster.options}{a vector of options for raster::writeRaster.}

\item{force.redo}{If an extraction for this template and label already exists,
should a new one be created?}
}
\value{
A \code{RasterLayer} cropped to the bounding box of the template.
}
\description{
\code{get_nlcd} returns a \code{RasterLayer} of NLCD data cropped to a given
template study area. \code{nlcd_colors} and \code{pal_nlcd} return the NLCD
legend and color palette, as available through the
[MLRC website](https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend).
}
\details{
NOTE: Prior to FedData version 3.0.0.9000, the `get_nlcd` function returned
data in the web Mercator coordinate reference system available through
the [MRLC web mapping services](https://www.mrlc.gov/geoserver/web/), rather
than data in the NLCD's native projection (a flavor of North American Albers).
Until the MRLC web services return data in the original projection, these
data are being served from a Google Cloud bucket of pre-processed cloud-optimized
GeoTIFFs. The script used to prepare the GeoTIFFs is available at
[https://github.com/bocinsky/feddata-nlcd](https://github.com/bocinsky/feddata-nlcd).
}
\examples{
\dontrun{
# Extract data for the Mesa Verde National Park:

# Get the NLCD (USA ONLY)
# Returns a raster
NLCD <-
  get_nlcd(
    template = FedData::meve,
    label = "meve",
    year = 2016
  )

# Plot with raster::plot
plot(NLCD)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITRDB_FUNCTIONS.R
\name{read_crn_data}
\alias{read_crn_data}
\title{Read chronology data from a Tucson-format chronology file.}
\usage{
read_crn_data(file, SCHWEINGRUBER)
}
\arguments{
\item{file}{A character string path pointing to a \code{*.crn} file to be read.}

\item{SCHWEINGRUBER}{Is the file in the Schweingruber-type Tucson format?}
}
\value{
A data.frame containing the data, or if \code{SCHWEINGRUBER==T}, a list containing four types of data.
}
\description{
This function includes improvements to the \code{\link{read_crn}} function from the
\pkg{dplR} library. The principle changes are better parsing of metadata, and support
for the Schweingruber-type Tucson format. Chronologies that are unable to be read
are reported to the user. The user (or \code{\link{read_crn}}) must tell the function whether
the file is a Schweingruber-type chronology.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{sequential_duplicated}
\alias{sequential_duplicated}
\title{Get a logical vector of which elements in a vector are sequentially duplicated.}
\usage{
sequential_duplicated(x, rows = F)
}
\arguments{
\item{x}{An vector of any type, or, if \code{rows}, a matrix.}

\item{rows}{Is x a matrix?}
}
\value{
A logical vector of the same length as x.
}
\description{
Get a logical vector of which elements in a vector are sequentially duplicated.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NHD_FUNCTIONS.R
\name{get_nhd}
\alias{get_nhd}
\title{Download and crop the National Hydrography Dataset.}
\usage{
get_nhd(
  template,
  label,
  nhdplus = FALSE,
  extraction.dir = paste0(tempdir(), "/FedData/"),
  force.redo = FALSE
)
}
\arguments{
\item{template}{A `Raster*`, `Spatial*`, or `sf` object to serve
as a template for cropping.}

\item{label}{A character string naming the study area.}

\item{nhdplus}{Extract data from the USGS NHDPlus High Resolution service (experimental)}

\item{extraction.dir}{A character string indicating where the extracted and cropped NHD data should be put.}

\item{force.redo}{If an extraction for this template and label already exists, should a new one be created?}
}
\value{
A list of `sf` collections extracted from the National Hydrography Dataset.
}
\description{
\code{get_nhd} returns a list of `sf` objects extracted
from the National Hydrography Dataset.
}
\examples{
\dontrun{
# Get the NHD (USA ONLY)
NHD <- get_nhd(
  template = FedData::meve,
  label = "meve"
)
NHD
NHD \%>\%
  plot_nhd(template = FedData::meve)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{extract_ssurgo_data}
\alias{extract_ssurgo_data}
\title{Extract data from a SSURGO database pertaining to a set of mapunits.}
\usage{
extract_ssurgo_data(tables, mapunits)
}
\arguments{
\item{tables}{A list of SSURGO tabular data.}

\item{mapunits}{A character vector of mapunits (likely dropped from SSURGO spatial data)
defining which mapunits to retain.}
}
\value{
A list of extracted SSURGO tabular data.
}
\description{
\code{extract_ssurgo_data} creates a directed graph of the joins in a SSURGO tabular dataset,
and then iterates through the tables, only retaining data pertinent to a set of mapunits.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UTILITY_FUNCTIONS.R
\name{unwrap_rows}
\alias{unwrap_rows}
\title{Unwraps a matrix and only keep the first n elements.}
\usage{
unwrap_rows(mat, n)
}
\arguments{
\item{mat}{A matrix}

\item{n}{A numeric vector}
}
\value{
A logical vector of the same length as x
}
\description{
A function that unwraps a matrix and only keeps the first n elements
n can be either a constant (in which case it will be repeated), or a vector
}
\keyword{internal}
