
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
