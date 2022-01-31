<!-- badges: start -->
[![peer-review](https://badges.ropensci.org/357_status.svg)](https://github.com/ropensci/software-review/issues/357)
[![status](https://joss.theoj.org/papers/3367fdbff2db55a60c1ab7d611017940/status.svg)](https://joss.theoj.org/papers/3367fdbff2db55a60c1ab7d611017940)
[![CRAN status](https://www.r-pkg.org/badges/version/chirps)](https://cran.r-project.org/package=chirps)
[![codecov](https://codecov.io/gh/ropensci/chirps/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/chirps)
[![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/225693680.svg)](https://zenodo.org/badge/latestdoi/225693680)
[![tic](https://github.com/ropensci/chirps/workflows/tic/badge.svg?branch=master)](https://github.com/ropensci/chirps/actions)
<!-- badges: end -->

# *chirps*: API Client for CHIRPS and CHIRTS <img align="right" src="man/figures/logo.png">

## Overview

**chirps** provides the API Client for the Climate Hazards Center 'CHIRPS' and 'CHIRTS'. The 'CHIRPS' data is a quasi-global (50°S – 50°N) high-resolution (0.05 arc-degrees) rainfall data set, which incorporates satellite imagery 
  and in-situ station data to create gridded rainfall time series for trend analysis and seasonal drought monitoring. 'CHIRTS' is a quasi-global (60°S – 70°N), high-resolution data set of daily maximum and minimum temperatures. For more details on 'CHIRPS' and 'CHIRTS' data please visit its official home page <https://www.chc.ucsb.edu/data>.

## Quick start

### From CRAN

The stable version is available through CRAN.

```r
install.packages("chirps")
```

### From GitHub

A development version that may have new features or bug fixes is available through GitHub.

``` r
library("remotes")

install_github("ropensci/chirps", build_vignettes = TRUE)
```

## Example

Fetch CHIRPS data from three points across the *Tapajós* National Forest (Brazil) from Jan-2017 to Dec-2017. The default procedure will download the COG files from the CHIRPS server and handle it internally using the package `terra`. This is more interesting when dealing with hundreds of points and days.

```r
library("chirps")

lonlat <- data.frame(lon = c(-55.0281,-54.9857, -55.0714),
                     lat = c(-2.8094, -2.8756, -3.5279))

dates <- c("2017-01-01", "2017-12-31")

dat <- get_chirps(lonlat, dates, server = "CHC")

```

For a faster download of few datapoints (~ 10 datapoints), the argument `server = "ClimateSERV"` can be used  

```r
library("chirps")

lonlat <- data.frame(lon = c(-55.0281,-54.9857, -55.0714),
                     lat = c(-2.8094, -2.8756, -3.5279))

dates <- c("2017-01-01", "2017-12-31")

dat <- get_chirps(lonlat, dates, server = "ClimateSERV")

```

## Going further

The full functionality of **chirps** is illustrated in the package vignette. The vignette can be found on the [package website](https://docs.ropensci.org/chirps/) or from within `R` once the package has been installed, e.g. via: 

``` r
vignette("Overview", package = "chirps")
```

## Use of CHIRPS data

While *chirps* does not redistribute the data or provide it in any way, we encourage users to cite Funk et al. (2015) when using CHIRPS and Funk et al. (2019) when using CHIRTS

> Funk C., Peterson P., Landsfeld M., … Michaelsen J. (2015). The climate hazards infrared precipitation with stations—a new environmental record for monitoring extremes. *Scientific Data*, 2, 150066. <https://doi.org/10.1038/sdata.2015.66>

> Funk, C., Peterson, P., Peterson, S., … Mata, N. (2019). A high-resolution 1983–2016 TMAX climate data record based on infrared temperatures and stations by the climate hazard center. *Journal of Climate*, 32(17), 5639–5658. <https://doi.org/10.1175/JCLI-D-18-0698.1>

## Meta

  - Please [report any issues or bugs](https://github.com/ropensci/chirps/issues).

  - License: MIT.

  - Get citation information for *chirps* in R by typing `citation(package = "chirps")`.

  - You are welcome to contribute to the *chirps* project. Please read our [contribution guidelines](CONTRIBUTING.md).

  - Please note that the *chirps* project is released with a a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
chirps 0.1.4.001 (2022-01-15) 
=========================

### ENHACEMENTS 

* Methods for objects of class "SpatExtent" in `get_chirps()` and `get_chirts()` to return a raster within a given area 


chirps 0.1.4 (2022-01-13) 
=========================

### ENHANCEMENTS 

* Add new function `get_chirts()` to fetch temperature data from CHC server (https://data.chc.ucsb.edu/products/CHIRTSdaily/v1.0/global_cogs_p05/)
* Implement data fetching from CHC server in `get_chirps()` which offers a better alternative for requests with multiple data points using GoC files from CHC server (https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/cogs/) and the `terra` package  
* New S3 methods in `get_chirps()` for objects of class 'SpatVector' and 'SpatRaster' from the `terra` package
* Data can be returned as an object of class 'matrix' when using the argument `as.matrix = TRUE` in the S3 methods for objects of class 'default', 'SpatVector' and 'SpatRaster'
* Updates the URL to request data from ClimateSERV
 

### CHANGES IN BEHAVIOUR

* New argument `server = ` is added to indicate from which server the function should send the request, either 'CHC' or 'ClimateSERV'. Please use the argument `server = "ClimateSERV"` for backward compatibility with previous versions of the package. 
* API requests to ClimateSERV use package httr instead of curl
* Argument `operation = ` in `get_chirps()` is only required when `server = "ClimateSERV"`
* Updates function `as.geojson()` to matches with the new requirements for ClimateSERV

chirps 0.1.3 (2021-07-10)
=========================

* GitHub version with ongoing updates and changes in behaviour. 


chirps 0.1.2 (2020-07-12)
=========================

* Add citation info for JOSS paper
* Fix vignette build
* A S3 method `as.geojson()` is added to replace the functions `data.frame_to_geojson()` and `sf_to_geojson()`

chirps 0.1.0 (2020-07-01)
=========================

* rOpenSci release version

### ENHANCEMENTS 

* Add `get_imerg()` to fetch IMERG data https://disasters.nasa.gov/instruments/imerg


chirps 0.0.8 (2020-05-22)
=========================

### CHANGES IN BEHAVIOUR

* Remove Imports of pkg 'tibble' which was basically to provide a "cool" print method. 
* A new `print()` method is added for objects that inherits the class 'chirps_df'
* Pkg 'methods' was moved from Imports to Depends
* Comments/suggestions given by Jake Zwart in rOpenSci pkg review are added


chirps 0.0.6 (2020-01-29)
=========================

### ENHANCEMENTS 

* Comments/suggestions given by Claudia Vitolo in rOpenSci pkg review are added
* `dataframe_to_geojson()`, `sf_to_geojson()` are added as exported functions avoiding `chirps:::`
* documentation for `tapajos` is given avoiding `chirps:::`


chirps 0.0.5 (2020-01-09)
=========================

### ENHANCEMENTS

* Fix comments given by good practice `goodpractice::gp()`. Avoid long code lines, it is bad for readability. Avoid 1:length(...), 1:nrow(...), 1:ncol(...), 1:NROW(...) and 1:NCOL(...) expressions. Not import packages as a whole, as this can cause name clashes between the imported packages. 


chirps 0.0.4 (2020-01-03)
=========================

### NEW FEATURES

* S3 methods for objects of class "geojson" in `get_chirps()` and `get_esi()`

* Package vignette

* Prepare for submission to ropensci

### ENHANCEMENTS

* Validations in internal functions to transform 'sf' into geojson

* Add properties features to geojson output in `get_chirps()` and `get_esi()` via `.add_geojson_properties()`


chirps 0.0.3 (2019-12-31)
=========================

### NEW FEATURES

* `get_esi` is added to retrieve Evaporative Stress Index with S3 methods for "data.frame" and "sf"

* S3 methods for objects of class "data.frame" and "sf" in `get_chirps`

### CHANGES IN BEHAVIOUR

* `.get_request_progress` and a `while` condition are added to check the progress of large requests and prevent the function to fail.

* `.GET` is added as a general function to retrieve other datasets from ClimateSERV

* improvements in internal functions documentation 


chirps 0.0.2 (2019-12-05)
=========================

### NEW FEATURES

* Calculate precipitation indices with `precip_indices` over a time span


chirps 0.0.1 (2019-12-03)
=========================

* GitHub-only release of prototype package.
# CONTRIBUTING #

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html) for further details.

### Code of Conduct

Please note that the *chirps* project is released with a a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 

---
title: 'chirps: API Client for the CHIRPS Precipitation Data in R'
tags:
- CHIRPS
- climate data
- climatology
- earth science
- evapotranspiration
- precipitation data
- R
- reproducibility
- weather data
authors:
  - name: Kauê de Sousa
    orcid: 0000-0002-7571-7845
    affiliation: "1, 2"
  - name: Adam H. Sparks
    orcid: 0000-0002-0061-8359
    affiliation: 3
  - name: William Ashmall
    affiliation: 4
  - name: Jacob van Etten
    orcid: 0000-0001-7554-2558
    affiliation: 2
  - name: Svein Ø. Solberg
    orcid: 0000-0002-4491-4483
    affiliation: 1
affiliations:
  - name: Department of Agricultural Sciences, Inland Norway University of Applied Sciences, Hamar, Norway
    index: 1
  - name: Bioversity International, Rome, Italy
    index: 2
  - name: Centre for Crop Health, University of Southern Queensland, Toowoomba, Australia
    index: 3
  - name: Universities Space Research Association, National Aeronautics and Space Administration (NASA), Huntsville, USA
    index: 4
citation_author: de Sousa et. al.
date: "22 May 2020"
year: 2020
bibliography: paper.bib
output: rticles::joss_article
journal: JOSS
---

# Summary

The *chirps* package provides functionalities for reproducible analysis in R [@RCoreTeam] using the CHIRPS [@Funk2015] data. CHIRPS is daily precipitation data set developed by the Climate Hazards Group [@Funk2015] for high resolution precipitation gridded data. Spanning 50$^{\circ}$ S to 50$^{\circ}$ N (and all longitudes) and ranging from 1981 to near-present (normally with a 45 day lag), CHIRPS incorporates 0.05 arc-degree resolution satellite imagery, and in-situ station data to create gridded precipitation time series for trend analysis and seasonal drought monitoring [@Funk2015]. Additionally, the package provides the API client for the IMERG [@Huffman2014] and ESI [@esi] data.
The Integrated Multi-satelliE Retrievals for GPM (IMERG) data provides near-real time global observations of rainfall at 0.5 arc-degree resolution, which can be used to estimate total rainfall accumulation from storm systems and quantify the intensity of rainfall and flood impacts from tropical cyclones and other storm systems. IMERG is a daily precipitation dataset available from 2015 to near-present. The evaporative stress index (ESI) data describes temporal anomalies in evapotranspiration produced weekly at 0.25 arc-degree resolution for the entire globe [@Anderson2011]. The ESI data is based on satellite observations of land surface temperature, which are used to estimate water loss due to evapotranspiration (the sum of evaporation and plant transpiration from the Earth's land and ocean surface to the atmosphere). The ESI data is available from 2001 to near-present. When using these data sets in publications please cite @Funk2015 for CHIRPS, @Huffman2014 for IMERG and @esi for ESI.

# Implementation

Four main functions are provided, `get_chirps()`, `get_imerg()`, `get_esi()` and `precip_indices()`. The `get_chirps()` function provides access to CHIRPS data via the ClimateSERV API Client [@ClimateSERV] with methods to handle objects of class 'data.frame', 'geojson' and 'sf' via the package *methods* [@RCoreTeam]. To accept the query, ClimateSERV requires a geojson object of type 'Polygon' (one single polygon per request). Using the package *sf* [@sf] internally, the input provided in `get_chirps()` is transformed into a list of polygons with a small buffer area (0.0001 arc-sec by default) around the point and transformed into a list of geojson strings. *chirps* uses *crul* [@crul] to interface with ClimateSERV API. The query returns a JSON object parsed to *jsonlite* [@jsonlite] to obtain the data frame for the time series required. `get_chirps()` returns a data.frame, which also inherits the classes 'chirps' and 'chirps_df', where each id represents the index for the rows in the in-putted 'object'. The function `get_imerg()` returns the precipitation data from the IMERG data set. The function works with the same parameters described for `get_chirps()` and also inherits the classes 'chirps' and 'chirps_df'. The function `get_esi()` returns the evaporative stress index (ESI) data [@Anderson2011], and works similarly to `get_chirps()` returning a data.frame which inherit the class 'chirps_df'. Users providing objects of class 'sf' and 'geojson' in `get_chirps()`, `get_imerg()` and `get_esi()` can also choose to return an object with the same class as the object provided using the arguments 'as.sf = TRUE' or 'as.geojson = TRUE'. With the function `precip_indices()` users can assess how the precipitation changes across the requested time series using precipitation variability indices [@Aguilar2005], computed using *stats* [@RCoreTeam], the main input is an object of class 'chirps'. Extended documentation is provided with examples on how to increase the buffer area and draw quadrants for the geojson polygon using *sf* [@sf].

# Application: a case study in the Tapajós National Forest

The *Tapajós* National Forest is a protected area in the Brazilian Amazon. Located within the coordinates -55.4$^{\circ}$ and -54.8$^{\circ}$ E and -4.1$^{\circ}$ and -2.7$^{\circ}$ S with ~ 527,400 ha of multiple Amazonian ecosystems. We take twenty random points across its area to get the precipitation from Jan-2008 to Dec-2018 using `get_chirps()`. We use an object of class 'sf' which is passed to the method `get_chirps.sf()`. Then, we compute the precipitation indices for the time series with intervals of 30 days using `precip_indices()`.

```r
library("chirps")
library("sf")

data("tapajos", package = "chirps")
set.seed(1234)
tp <- st_sample(tapajos, 20)
tp <- st_as_sf(tp)

dt <- get_chirps(tp, dates = c("2008-01-01","2018-01-31"))

p_ind <- precip_indices(dt, timeseries = TRUE, intervals = 30)

```

We selected four indices for the visualization using *tidyverse* [@tidyverse]. Plots were ensembled together using *gridExtra* [@gridExtra]. Here we see how these indices are changing across the time series (Figure 1). In this quick assessment, we note an increasing extent of consecutive dry days (MLDS) across the time series, with also a decrease in the number of consecutive rainy days (MLWS), which stays above the historical average for MLDS and bellow the historical average for MLWS. The trends also show a decrease in the total rainfall in the 30-days intervals, staying below the average after 2014. Finally, we note a decrease in maximum consecutive 5-days precipitation, which also stays bellow the historical average. 

\begin{figure}
\includegraphics[width=0.9\linewidth]{Fig1} \caption{Trends in precipitation variability across the Tapajós National Forest, Brazil, for the period of 01-Jan-2010 to 31-Dec-2018 with four precipitation indices. MLDS, maximum length of consecutive dry days (days), MLWS, maximum length of consecutive wet days (days), Rtotal, total precipitation (mm), Rx5day, maximum consecutive 5-days precipitation (mm). Red lines indicates the historical mean of each index in the time series. Blue line indicates the smoothed trends in each index using the 'loess' method.}\label{fig:fig1}
\end{figure}

# Other applications and conclusion

Deriving precipitation indices that can be obtained from CHIRPS proved to be an excellent approach to evaluate the climate variability using precipitation data [@DeSousa2018] and the effects of climate change on a continental analysis [@Aguilar2005]. Additionally, these indices can be used to register specific effects of climate variability on crop varietal performance. In crop modelling, @Kehel2016 applied this to assess the interactions of wheat varieties with the environment, showing how severe drought, assessed with the maximum length of dry spell (MLDS), can affect the plant development and the yield. These indices can also be useful to improve variety recommendation for climate adaptation in marginal production environments [@vanEtten2019]. 

Overall, CHIRPS data can be used in many applications and currently has over 800 citations from studies using this tool. Many applications are the field of agriculture, hydrologic modelling and drought monitoring, but also some studies using this in disease control programs (e.g. @Thomson2017, @Horn2018). The *chirps* package aims to make it possible for researchers in these fields to implement this tool into a replicable and reproducible workflow in R. 

# Acknowledgements

This work was supported by The Nordic Joint Committee for Agricultural and Food Research (grant num. 202100-2817). The idea for this package was conceived during the course "Analysing Spatial Data" at the Norwegian School of Economics (NHH), we thank Professor Roger Bivand for his insights.

# References
---
title: "Introduction to chirps"
package: chirps
author:
- name: Kauê de Sousa
  affiliation: Department of Agricultural Sciences, Inland Norway University, Hamar, Norway </br> The Alliance of Bioversity International and CIAT, Montpellier, France
- name: Adam H. Sparks 
  affiliation: Centre for Crop Health, University of Southern Queensland, Toowoomba, Australia
- name: Aniruddha Ghosh
  affiliation: The Alliance of Bioversity International and CIAT, Nairobi, Kenya
output: html_document
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Introduction to chirps}
  %\usepackage[UTF-8]{inputenc}
  %\VignetteEncoding{UTF-8}
bibliography: ["chirps.bib"]
csl: citation_style.csl
---



# Summary

The **chirps** package [@chirps] provides functionalities for reproducible analysis using the CHIRPS [@Funk2015] and CHIRTS data [@Funk2019] . CHIRPS is a daily precipitation data set developed by the [Climate Hazards Group](https://www.chc.ucsb.edu/) for high resolution precipitation gridded data. Spanning 50°S - 50°N (and all longitudes) and ranging from 1981 to near-present (normally with a 45 day lag), CHIRPS incorporates 0.05 arc-degree resolution satellite imagery, and in-situ station data to create gridded precipitation time series for trend analysis and seasonal drought monitoring. CHIRTS is a quasi-global (60°S – 70°N), high-resolution data set of daily maximum and minimum temperatures. 

Other functionalities of **chirps** are the computation of precipitation indices, the retrieval of the evaporative stress index (ESI) which describes temporal anomalies in evapotranspiration produced weekly at 0.25 arc-degree resolution for the entire globe, and the retrieval of IMERG data which provides near-real time global observations of rainfall at 0.5 arc-degree resolution.

# CHIRPS (precipitation data)

The *Tapajós* National Forest is a protected area in the Brazilian Amazon. Located within the coordinates -55.4° and -54.8°E and -4.1° and -2.7°S with ~527,400 ha of multiple Amazonian ecosystems. We take three points within its area to get the precipitation from Jan-2013 to Dec-2018 using `get_chirps()`.

<img src="map.png" title="plot of chunk map" alt="plot of chunk map" style="display: block; margin: auto;" />

For this example we fetch the data from the server "ClimateSERV" using the argument `server = "ClimateSERV"`. This option is recommended when working with few data points as the request could be faster. The default `server = "CHC"` is used for multiple data points and years. 


```r

library("chirps")
library("sf")

data("tapajos", package = "chirps")

# sample three points within the Tapajos area
set.seed(1234)
tp_point <- st_sample(tapajos, 3)

# coerce as sf points
tp_point <- st_as_sf(tp_point)

dat <- get_chirps(tp_point,
                  dates = c("2013-01-01","2018-12-31"), 
                  server = "ClimateSERV")
#> Fetching data from ClimateSERV
#> Getting your request...
```

## Precipitation indices

By default, the function `get_chirps()` returns a data.frame which inherits the classes 'chirps' and 'chirps_df', where each id represents the index for the rows in the in-putted 'object'. It is possible to return the data as a matrix using the argument `as.matrix = TRUE`.


```r
dat
#>          id    lon   lat       date chirps
#>       <int>  <dbl> <dbl>     <date>  <dbl>
#> 1:        1 -55.03 -3.80 2013-01-01   0.00
#> 2:        1 -55.03 -3.80 2013-01-02  12.36
#> 3:        1 -55.03 -3.80 2013-01-03  24.72
#> 4:        1 -55.03 -3.80 2013-01-04   0.00
#> 5:        1 -55.03 -3.80 2013-01-05   0.00
#> ---                                       
#> 6569:     3 -55.03 -3.41 2018-12-27   0.00
#> 6570:     3 -55.03 -3.41 2018-12-28   0.00
#> 6571:     3 -55.03 -3.41 2018-12-29   0.00
#> 6572:     3 -55.03 -3.41 2018-12-30   0.00
#> 6573:     3 -55.03 -3.41 2018-12-31   0.00
```

With `precip_indices()` is possible to assess how the precipitation changes across a time series using precipitation variability indices [@Aguilar2005]. Here, we take the indices for intervals of 15 days and compute the indices for the time series (from Jan-2013 to Dec-2018).


```r
p_ind <- precip_indices(dat, timeseries = TRUE, intervals = 15)

p_ind
#>          id       date    lon   lat  index  value
#>       <int>     <date>  <dbl> <dbl>  <chr>  <dbl>
#> 1:        1 2013-01-01 -55.03 -3.80   MLDS   7.00
#> 2:        1 2013-01-01 -55.03 -3.80   MLWS   2.00
#> 3:        1 2013-01-01 -55.03 -3.80  R10mm   1.00
#> 4:        1 2013-01-01 -55.03 -3.80  R20mm   3.00
#> 5:        1 2013-01-01 -55.03 -3.80 Rx1day  45.70
#> ---                                              
#> 3446:     3 2018-12-16 -55.03 -3.41 Rx5day  53.90
#> 3447:     3 2018-12-16 -55.03 -3.41   R95p  34.53
#> 3448:     3 2018-12-16 -55.03 -3.41   R99p  34.53
#> 3449:     3 2018-12-16 -55.03 -3.41 Rtotal  80.49
#> 3450:     3 2018-12-16 -55.03 -3.41   SDII  13.42
```

The function `precip_indices()` returns a data.frame with the precipitation indices. Each date corresponds to the first day in the time series intervals as defined by the argument 'intervals'. When `timeseries = FALSE` the function returns a single precipitation index for the entire time series.

# CHIRTS (temperature data)

Maximum and minimum temperature and relative humidity data are available with the function `get_chirts()`. Data is requested to the server CHC as default and is currently available from 1983 to 2016. We use the same random points from the Tapajós National Forest but for few days to speed up the call. 



```r

dates <- c("2010-12-15","2010-12-31")

temp1 <- get_chirts(tp_point, dates, var = "Tmax", as.matrix = TRUE)

temp2 <- get_chirts(tp_point, dates, var = "Tmin", as.matrix = TRUE)

rhu <- get_chirts(tp_point, dates, var = "RHum", as.matrix = TRUE)

```


# Going further

## Evapotranspiration 

The **chirps** package also retrieves the Evaporative Stress Index (ESI) using the function `get_esi()` which behaves similarly as `get_chirps()`. 


```r

dt <- get_esi(tp_point, c("2016-05-01","2016-12-31"))

```

The function `get_esi()` may return `NA`s due to cloudiness in the dataset. Which will return an error message:


```r
set.seed(123)
lonlat <- data.frame(lon = runif(1, -55, -54),
                     lat = runif(1, -3, -2.7))

get_esi(lonlat, c("2017-12-01","2018-01-20"))

```

One way to deal with this is increase the buffer area around the in-putted object with the argument `dist` passed to `st_buffer()` from *sf*[@sf] through the `...` functionality in `get_esi()`. The argument `nQuadSegs` defines the number of segments per quadrant in the buffer.  


```r

get_esi(lonlat, c("2017-12-01","2018-01-20"), dist = 0.1, nQuadSegs = 6)

```

## Objects of class sf

To return an object with the same class (`sf`), the argument `as.sf = TRUE` is used.



```r

get_chirps(tapajos, dates = c("2017-12-15","2017-12-31"), as.sf = TRUE)

```
## Objects of class geojson

`get_chirps()` and `get_esi()` also contains a method for objects of class geojson with geometries 'Point' and 'Polygon'. To return an object with the same class (`geojson`), the argument `as.geojson = TRUE` is used.



```r

tp_gjson <- sf_to_geojson(tp_point)

dt <- get_esi(tp_gjson, dates = c("2017-12-15","2017-12-31"), dist = 0.1)

```

# References




---
title: "Introduction to chirps"
package: chirps
author:
- name: Kauê de Sousa
  affiliation: Department of Agricultural Sciences, Inland Norway University, Hamar, Norway </br> The Alliance of Bioversity International and CIAT, Montpellier, France
- name: Adam H. Sparks 
  affiliation: Centre for Crop Health, University of Southern Queensland, Toowoomba, Australia
- name: Aniruddha Ghosh
  affiliation: The Alliance of Bioversity International and CIAT, Nairobi, Kenya
output: html_document
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Introduction to chirps}
  %\usepackage[UTF-8]{inputenc}
  %\VignetteEncoding{UTF-8}
bibliography: ["chirps.bib"]
csl: citation_style.csl
---



# Summary

The **chirps** package [@chirps] provides functionalities for reproducible analysis using the CHIRPS [@Funk2015] and CHIRTS data [@Funk2019] . CHIRPS is a daily precipitation data set developed by the [Climate Hazards Group](https://www.chc.ucsb.edu/) for high resolution precipitation gridded data. Spanning 50°S - 50°N (and all longitudes) and ranging from 1981 to near-present (normally with a 45 day lag), CHIRPS incorporates 0.05 arc-degree resolution satellite imagery, and in-situ station data to create gridded precipitation time series for trend analysis and seasonal drought monitoring. CHIRTS is a quasi-global (60°S – 70°N), high-resolution data set of daily maximum and minimum temperatures. 

Other functionalities of **chirps** are the computation of precipitation indices, the retrieval of the evaporative stress index (ESI) which describes temporal anomalies in evapotranspiration produced weekly at 0.25 arc-degree resolution for the entire globe, and the retrieval of IMERG data which provides near-real time global observations of rainfall at 0.5 arc-degree resolution.

# CHIRPS (precipitation data)

The *Tapajós* National Forest is a protected area in the Brazilian Amazon. Located within the coordinates -55.4° and -54.8°E and -4.1° and -2.7°S with ~527,400 ha of multiple Amazonian ecosystems. We take three points within its area to get the precipitation from Jan-2013 to Dec-2018 using `get_chirps()`.

<img src="map.png" title="plot of chunk map" alt="plot of chunk map" style="display: block; margin: auto;" />

For this example we fetch the data from the server "ClimateSERV" using the argument `server = "ClimateSERV"`. This option is recommended when working with few data points as the request could be faster. The default `server = "CHC"` is used for multiple data points and years. 


```r

library("chirps")
library("sf")

data("tapajos", package = "chirps")

# sample three points within the Tapajos area
set.seed(1234)
tp_point <- st_sample(tapajos, 3)

# coerce as sf points
tp_point <- st_as_sf(tp_point)

dat <- get_chirps(tp_point,
                  dates = c("2013-01-01","2018-12-31"), 
                  server = "ClimateSERV")
#> Fetching data from ClimateSERV
#> Getting your request...
```

## Precipitation indices

By default, the function `get_chirps()` returns a data.frame which inherits the classes 'chirps' and 'chirps_df', where each id represents the index for the rows in the in-putted 'object'. It is possible to return the data as a matrix using the argument `as.matrix = TRUE`.


```r
dat
#>          id    lon   lat       date chirps
#>       <int>  <dbl> <dbl>     <date>  <dbl>
#> 1:        1 -55.03 -3.80 2013-01-01   0.00
#> 2:        1 -55.03 -3.80 2013-01-02  12.36
#> 3:        1 -55.03 -3.80 2013-01-03  24.72
#> 4:        1 -55.03 -3.80 2013-01-04   0.00
#> 5:        1 -55.03 -3.80 2013-01-05   0.00
#> ---                                       
#> 6569:     3 -55.03 -3.41 2018-12-27   0.00
#> 6570:     3 -55.03 -3.41 2018-12-28   0.00
#> 6571:     3 -55.03 -3.41 2018-12-29   0.00
#> 6572:     3 -55.03 -3.41 2018-12-30   0.00
#> 6573:     3 -55.03 -3.41 2018-12-31   0.00
```

With `precip_indices()` is possible to assess how the precipitation changes across a time series using precipitation variability indices [@Aguilar2005]. Here, we take the indices for intervals of 15 days and compute the indices for the time series (from Jan-2013 to Dec-2018).


```r
p_ind <- precip_indices(dat, timeseries = TRUE, intervals = 15)

p_ind
#>          id       date    lon   lat  index  value
#>       <int>     <date>  <dbl> <dbl>  <chr>  <dbl>
#> 1:        1 2013-01-01 -55.03 -3.80   MLDS   7.00
#> 2:        1 2013-01-01 -55.03 -3.80   MLWS   2.00
#> 3:        1 2013-01-01 -55.03 -3.80  R10mm   1.00
#> 4:        1 2013-01-01 -55.03 -3.80  R20mm   3.00
#> 5:        1 2013-01-01 -55.03 -3.80 Rx1day  45.70
#> ---                                              
#> 3446:     3 2018-12-16 -55.03 -3.41 Rx5day  53.90
#> 3447:     3 2018-12-16 -55.03 -3.41   R95p  34.53
#> 3448:     3 2018-12-16 -55.03 -3.41   R99p  34.53
#> 3449:     3 2018-12-16 -55.03 -3.41 Rtotal  80.49
#> 3450:     3 2018-12-16 -55.03 -3.41   SDII  13.42
```

The function `precip_indices()` returns a data.frame with the precipitation indices. Each date corresponds to the first day in the time series intervals as defined by the argument 'intervals'. When `timeseries = FALSE` the function returns a single precipitation index for the entire time series.

# CHIRTS (temperature data)

Maximum and minimum temperature and relative humidity data are available with the function `get_chirts()`. Data is requested to the server CHC as default and is currently available from 1983 to 2016. We use the same random points from the Tapajós National Forest but for few days to speed up the call. 



```r

dates <- c("2010-12-15","2010-12-31")

temp1 <- get_chirts(tp_point, dates, var = "Tmax", as.matrix = TRUE)

temp2 <- get_chirts(tp_point, dates, var = "Tmin", as.matrix = TRUE)

rhu <- get_chirts(tp_point, dates, var = "RHum", as.matrix = TRUE)

```


# Going further

## Evapotranspiration 

The **chirps** package also retrieves the Evaporative Stress Index (ESI) using the function `get_esi()` which behaves similarly as `get_chirps()`. 


```r

dt <- get_esi(tp_point, c("2016-05-01","2016-12-31"))

```

The function `get_esi()` may return `NA`s due to cloudiness in the dataset. Which will return an error message:


```r
set.seed(123)
lonlat <- data.frame(lon = runif(1, -55, -54),
                     lat = runif(1, -3, -2.7))

get_esi(lonlat, c("2017-12-01","2018-01-20"))

```

One way to deal with this is increase the buffer area around the in-putted object with the argument `dist` passed to `st_buffer()` from *sf*[@sf] through the `...` functionality in `get_esi()`. The argument `nQuadSegs` defines the number of segments per quadrant in the buffer.  


```r

get_esi(lonlat, c("2017-12-01","2018-01-20"), dist = 0.1, nQuadSegs = 6)

```

## Objects of class sf

To return an object with the same class (`sf`), the argument `as.sf = TRUE` is used.



```r

get_chirps(tapajos, dates = c("2017-12-15","2017-12-31"), as.sf = TRUE)

```
## Objects of class geojson

`get_chirps()` and `get_esi()` also contains a method for objects of class geojson with geometries 'Point' and 'Polygon'. To return an object with the same class (`geojson`), the argument `as.geojson = TRUE` is used.



```r

tp_gjson <- sf_to_geojson(tp_point)

dt <- get_esi(tp_gjson, dates = c("2017-12-15","2017-12-31"), dist = 0.1)

```

# References




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chirps.R
\name{get_chirps}
\alias{get_chirps}
\alias{get_chirps.default}
\alias{get_chirps.SpatVector}
\alias{get_chirps.SpatRaster}
\alias{get_chirps.sf}
\alias{get_chirps.geojson}
\alias{get_chirps.SpatExtent}
\title{Get CHIRPS precipitation data}
\usage{
get_chirps(object, dates, server, ...)

\method{get_chirps}{default}(object, dates, server, as.matrix = FALSE, ...)

\method{get_chirps}{SpatVector}(object, dates, server = "CHC", as.raster = TRUE, ...)

\method{get_chirps}{SpatRaster}(
  object,
  dates,
  server = "CHC",
  as.matrix = TRUE,
  as.raster = FALSE,
  ...
)

\method{get_chirps}{sf}(object, dates, server, as.sf = FALSE, ...)

\method{get_chirps}{geojson}(object, dates, server, as.geojson = FALSE, ...)

\method{get_chirps}{SpatExtent}(object, dates, server = "CHC", as.raster = TRUE, ...)
}
\arguments{
\item{object}{input, an object of class \code{\link[base]{data.frame}} (or
any other object that can be coerced to data.frame), \code{\link[terra]{SpatVector}}, 
\code{\link[terra]{SpatRaster}}, \code{\link[terra]{SpatExtent}},
\code{\link[sf]{sf}} or \code{geojson}}

\item{dates}{a character of start and end dates in that order in the format
"YYYY-MM-DD"}

\item{server}{a character that represent the server source "CHC" or
"ClimateSERV"}

\item{...}{additional arguments passed to \code{\link[terra]{terra}} 
or \code{\link[sf]{sf}} methods
See details}

\item{as.matrix}{logical, returns an object of class \code{matrix}}

\item{as.raster}{logical, returns an object of class \code{\link[terra]{SpatRaster}}}

\item{as.sf}{logical, returns an object of class \code{\link[sf]{sf}}}

\item{as.geojson}{logical, returns an object of class \code{geojson}}
}
\value{
A matrix, raster or a data frame of \acronym{CHIRPS} data:
\describe{
  \item{id}{the index for the rows in \code{object}}
  \item{dates}{the dates from which \acronym{CHIRPS} was requested}
  \item{lon}{the longitude as provided in \code{object}}
  \item{lat}{the latitude as provided in \code{object}}
  \item{chirps}{the \acronym{CHIRPS} value in mm}
}
}
\description{
Get daily precipitation data from the "Climate Hazards Group". Two server 
 sources are available. The first, "CHC" (default) is recommended for 
 multiple data-points, while "ClimateSERV" is recommended when 
 few data-points are required (~ 50).
}
\details{
Data description at 
\url{https://data.chc.ucsb.edu/products/CHIRPS-2.0/README-CHIRPS.txt}

\strong{Additional arguments when using server = "CHC"}

\bold{resolution}: numeric, resolution of CHIRPS tiles either 
 0.05 (default) or 0.25 degrees

\strong{Additional arguments when using server = "ClimateSERV"}

\bold{dist}: numeric, buffer distance for each \code{object} coordinate

\bold{nQuadSegs}: integer, number of segments per buffer quadrant

\bold{operation}: supported operations for ClimateSERV are:
 \tabular{rll}{
 \bold{operation}      \tab    \tab \bold{value}\cr
 max                   \tab =  \tab 0\cr
 min                   \tab =  \tab 1\cr
 median                \tab =  \tab 2\cr
 sum                   \tab =  \tab 4\cr
 average               \tab =  \tab 5 (\emph{default value})\cr
 }
}
\note{
get_chirps may return some warning messages given by 
\code{\link[sf]{sf}}, please look sf documentation for 
possible issues.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library("chirps")
library("terra")

# Case 1: return as a data.frame
dates <- c("2017-12-15","2017-12-31")
lonlat <- data.frame(lon = c(-55.0281,-54.9857), lat = c(-2.8094, -2.8756))

r1 <- get_chirps(lonlat, dates, server = "CHC")

# Case 2: return a matrix
r2 <- get_chirps(lonlat, dates, server = "CHC", as.matrix = TRUE)

# Case 3: input SpatVector and return raster
f <- system.file("ex/lux.shp", package = "terra")
v <- vect(f)
r3 <- get_chirps(v, dates, server = "CHC", as.raster = TRUE)

# Case 4: input SpatExtent and return a raster within the extent
area <- ext(c(-66, -64, -6, -4))

dates <- c("2017-12-15", "2017-12-31")

r4 <- get_chirps(area, dates, server = "CHC")

# Case 5: using the server "ClimateSERV"
r5 <- get_chirps(lonlat, dates, server = "ClimateSERV")

# Case 6: from "ClimateSERV" and return as a matrix
r6 <- get_chirps(lonlat, dates, server = "ClimateSERV", as.matrix = TRUE)

\dontshow{\}) # examplesIf}
}
\references{
Funk C. et al. (2015). Scientific Data, 2, 150066.
 \cr\doi{10.1038/sdata.2015.66}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precip_indices.R
\name{precip_indices}
\alias{precip_indices}
\title{Compute precipitation indices over a time series.}
\usage{
precip_indices(object, timeseries = FALSE, intervals = NULL)
}
\arguments{
\item{object}{an object of class \code{chirps} as provided by
\code{\link{get_chirps}}}

\item{timeseries}{logical, \code{FALSE} for a single point time series
observation or \code{TRUE} for a time series based on \var{intervals}}

\item{intervals}{integer no lower than 5, for the days intervals when
\var{timeseries} = \code{TRUE}}
}
\value{
A dataframe with precipitation indices:
\item{MLDS}{maximum length of consecutive dry day, rain < 1 mm (days)}
\item{MLWS}{maximum length of consecutive wet days, rain >= 1 mm (days)}
\item{R10mm}{number of heavy precipitation days 10 >= rain < 20 mm (days)}
\item{R20mm}{number of very heavy precipitation days rain >= 20 (days)}
\item{Rx1day}{maximum 1-day precipitation (mm)}
\item{Rx5day}{maximum 5-day precipitation (mm)}
\item{R95p}{total precipitation when rain > 95th percentile (mm)}
\item{R99p}{total precipitation when rain > 99th percentile (mm)}
\item{Rtotal}{total precipitation (mm) in wet days, rain >= 1 (mm)}
\item{SDII}{simple daily intensity index, total precipitation divided by the
 number of wet days (mm/days)}
}
\description{
Compute precipitation indices over a time series.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
lonlat <- data.frame(lon = c(-55.0281,-54.9857),
                     lat = c(-2.8094, -2.8756))

dates <- c("2017-12-15", "2017-12-31")

dt <- get_chirps(lonlat, dates, server = "ClimateSERV")

# take the indices for the entire period
precip_indices(dt, timeseries = FALSE)

# take the indices for periods of 7 days
precip_indices(dt, timeseries = TRUE, intervals = 7)
\dontshow{\}) # examplesIf}
}
\references{
Aguilar E., et al. (2005). Journal of Geophysical Research, 110(D23), D23107.

Kehel Z., et al. (2016). In: Applied Mathematics and Omics to Assess Crop
 Genetic Resources for Climate Change Adaptive Traits (eds Bari A., Damania
 A. B., Mackay M., Dayanandan S.), pp. 151–174. CRC Press.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chirts.R
\name{get_chirts}
\alias{get_chirts}
\alias{get_chirts.default}
\alias{get_chirts.SpatVector}
\alias{get_chirts.SpatRaster}
\alias{get_chirts.SpatExtent}
\title{Get CHIRTS temperature data data}
\usage{
get_chirts(object, dates, var, ...)

\method{get_chirts}{default}(object, dates, var, as.matrix = FALSE, ...)

\method{get_chirts}{SpatVector}(object, dates, var, as.raster = TRUE, ...)

\method{get_chirts}{SpatRaster}(object, dates, var, as.raster = TRUE, ...)

\method{get_chirts}{SpatExtent}(object, dates, var, as.raster = TRUE, ...)
}
\arguments{
\item{object}{an object of class \code{\link[base]{data.frame}} (or any other 
object that can be coerced to a \code{data.frame}), \code{\link[terra]{SpatExtent}},
\code{\link[terra]{SpatVector}}, or \code{\link[terra]{SpatRaster}}}

\item{dates}{a character of start and end dates in that order in the format
"YYYY-MM-DD"}

\item{var}{character, A valid variable from the options: \dQuote{Tmax},
\dQuote{Tmin}, \dQuote{RHum} and \dQuote{HeatIndex}}

\item{...}{additional arguments passed to \code{\link[terra]{terra}}}

\item{as.matrix}{logical, returns an object of class \code{matrix}}

\item{as.raster}{logical, returns an object of class \code{\link[terra]{SpatRaster}}}
}
\value{
A SpatRaster object if \code{as.raster=TRUE}, else \code{matrix}, 
\code{list}, or \code{data.frame}
}
\description{
Get daily maximum and minimum temperature data from the "Climate Hazards
 Group". CHIRTS-daily is a global 2-m temperature product that combines the
 monthly CHIRTSmax data set with the ERA5 reanalysis to produce routinely
 updated data to support the monitoring of temperature extreme. Data is
 currently available from 1983 to 2016. Soon available to near-present.
}
\details{
Variable description from 
\url{https://data.chc.ucsb.edu/products/CHIRTSdaily/aaa.Readme.txt}
\describe{
  \item{Tmax}{Daily average maximum air temperature at 2 m above ground}
  \item{Tmin}{Daily average minimum air temperature at 2 m above ground}
  \item{RHum}{Daily average relative humidity}
  \item{HeatIndex}{Daily average heat index}
  }
}
\section{Additional arguments}{
 
\bold{interval}: supported intervals are \dQuote{daily}, \dQuote{pentad},
 \dQuote{dekad}, \dQuote{monthly}, \dQuote{2-monthly}, \dQuote{3-monthly},
 and \dQuote{annual}. Currently hard coded to \dQuote{daily}.
}

\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library("chirps")
library("terra")

# Case 1: input a data frame return a data frame in the long format
dates <- c("2010-12-15","2010-12-31")
lonlat <- data.frame(lon = c(-55.0281,-54.9857),
                     lat = c(-2.8094, -2.8756))

temp1 <- get_chirts(lonlat, dates, var = "Tmax")

# Case 2: input a data frame return a matrix
temp2 <- get_chirts(lonlat, dates, "Tmax", as.matrix = TRUE)

# Case 3: input a raster and return raster
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
temp3 <- get_chirts(v, dates, var = "Tmax", as.raster = TRUE)

# Case 4: input a raster and return raster
temp4 <- get_chirts(v, dates, var = "Tmax", as.matrix = TRUE)
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tapajos.R
\docType{data}
\name{tapajos}
\alias{tapajos}
\title{Tapajos National Forest}
\format{
An object of class 'sfc_POLYGON' within the bounding box 
xmin: -55.41127 ymin: -4.114584 
xmax: -54.7973 ymax: -2.751706
}
\source{
The data was provided by the Chico Mendes Institute via
\url{https://www.protectedplanet.net/en}
}
\usage{
tapajos
}
\description{
Geometries for the Tapajos National Forest, a protected 
area in the Brazilian Amazon
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_imerg.R
\name{get_imerg}
\alias{get_imerg}
\alias{get_imerg.default}
\alias{get_imerg.sf}
\alias{get_imerg.geojson}
\title{Get Integrated Multisatellite Retrievals for GPM (IMERG) data}
\usage{
get_imerg(object, dates, operation = 5, ...)

\method{get_imerg}{default}(object, dates, operation = 5, ...)

\method{get_imerg}{sf}(object, dates, operation = 5, as.sf = FALSE, ...)

\method{get_imerg}{geojson}(object, dates, operation = 5, as.geojson = FALSE, ...)
}
\arguments{
\item{object}{input, an object of class \code{\link[base]{data.frame}} (or
any other object that can be coerced to data.frame), \code{\link[terra]{SpatVector}}, 
\code{\link[terra]{SpatRaster}}, \code{\link[terra]{SpatExtent}},
\code{\link[sf]{sf}} or \code{geojson}}

\item{dates}{a character of start and end dates in that order in the format
"YYYY-MM-DD"}

\item{operation}{optional, an integer that represents which type of
statistical operation to perform on the dataset}

\item{...}{additional arguments passed to \code{\link[terra]{terra}} 
or \code{\link[sf]{sf}} methods
See details}

\item{as.sf}{logical, returns an object of class \code{\link[sf]{sf}}}

\item{as.geojson}{logical, returns an object of class \code{geojson}}
}
\value{
A data frame of \acronym{imerg} data:
\item{id}{the index for the rows in \code{object}}
\item{dates}{the dates from which imerg was requested}
\item{lon}{the longitude as provided in \code{object}}
\item{lat}{the latitude as provided in \code{object}}
\item{imerg}{the IMERG value}
}
\description{
The IMERG dataset provides near-real time global observations of 
 rainfall at 10km resolution, which can be used to estimate total
 rainfall accumulation from storm systems and quantify the intensity 
 of rainfall and flood impacts from tropical cyclones and other storm 
 systems. IMERG is a daily precipitation dataset available from 2015 
 to present within the latitudes 70 and -70.
}
\details{
\bold{operation}: supported operations are:  
 \tabular{rll}{
 \bold{operation}      \tab    \tab \bold{value}\cr
 max                   \tab =  \tab 0\cr
 min                   \tab =  \tab 1\cr
 median                \tab =  \tab 2\cr
 sum                   \tab =  \tab 4\cr
 average               \tab =  \tab 5 (\emph{default value})\cr
 }

\bold{dist}: numeric, buffer distance for each \code{object} coordinate

\bold{nQuadSegs}: integer, number of segments per buffer quadrant
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
lonlat <- data.frame(lon = c(-55.0281,-54.9857),
                     lat = c(-2.8094, -2.8756))

dates <- c("2017-12-15", "2017-12-31")

dt <- get_imerg(lonlat, dates)

dt
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chirps.R
\docType{package}
\name{chirps}
\alias{chirps}
\alias{chirps-package}
\title{API Client for CHIRPS and CHIRTS}
\description{
API Client for the Climate Hazards Center 'CHIRPS' and 'CHIRTS'. The 'CHIRPS' data is a quasi-global (50°S – 50°N) high-resolution (0.05 arc-degrees) rainfall data set, which incorporates satellite imagery and in-situ station data to create gridded rainfall time series for trend analysis and seasonal drought monitoring. 'CHIRTS' is a quasi-global (60°S – 70°N), high-resolution data set of daily maximum and minimum temperatures. For more details on 'CHIRPS' and 'CHIRTS' data please visit its official home page <https://www.chc.ucsb.edu/data>.
}
\note{
While \CRANpkg{chirps} does not redistribute the data or provide it in any
 way, we encourage users to cite Funk et al. (2015) when using
 \acronym{CHIRPS} and Funk et al. (2019) when using \acronym{CHIRTS}.

Funk et al. (2015). Scientific Data, 2, 150066. 
\doi{10.1038/sdata.2015.66}

Funk et al. (2019). Journal of Climate, 32(17), 5639–5658. 
\doi{10.1175/JCLI-D-18-0698.1}
}
\seealso{
\strong{Useful links:}
\itemize{
\item{JOSS paper: 
 \doi{10.21105/joss.02419}}
\item{Development repository: 
 \url{https://github.com/ropensci/chirps}}
\item{Static documentation: 
 \url{https://docs.ropensci.org/chirps/}}
\item{Report bugs: 
 \url{https://github.com/ropensci/chirps/issues}}
\item{CHC website: 
 \url{https://www.chc.ucsb.edu}}
}
}
\author{
Kauê de Sousa and Adam H. Sparks and Aniruddha Ghosh
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.geojson.R
\name{as.geojson}
\alias{as.geojson}
\alias{as.geojson.default}
\alias{as.geojson.sf}
\title{Methods to coerce geographical coordinates into a geojson polygon}
\usage{
as.geojson(lonlat, dist = 1e-05, nQuadSegs = 2L, ...)

\method{as.geojson}{default}(lonlat, dist = 1e-05, nQuadSegs = 2L, ...)

\method{as.geojson}{sf}(lonlat, dist = 1e-05, nQuadSegs = 2L, ...)
}
\arguments{
\item{lonlat}{a data.frame or matrix with geographical coordinates lonlat, in 
that order, or an object of class 'sf' with geometry type 'POINT' or 'POLYGON'}

\item{dist}{numeric, buffer distance for all \code{lonlat}}

\item{nQuadSegs}{integer, number of segments per quadrant}

\item{...}{further arguments passed to \code{\link[sf]{sf}} methods}
}
\value{
An object of class 'geosjon' for each row in \code{lonlat}
}
\description{
Take single points from geographical coordinates and coerce into a
geojson of geometry 'Polygon'
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Default S3 Method
# random geographic points within bbox(10, 12, 45, 47)
library("sf")

set.seed(123)
lonlat <- data.frame(lon = runif(1, 10, 12),
                     lat = runif(1, 45, 47))

gjson <- as.geojson(lonlat)

#################

# S3 Method for objects of class 'sf'
# random geographic points within bbox(10, 12, 45, 47)
library("sf")

set.seed(123)
lonlat <- data.frame(lon = runif(5, 10, 12),
                     lat = runif(5, 45, 47))

lonlat <- st_as_sf(lonlat, coords = c("lon","lat"))

gjson <- as.geojson(lonlat)
\dontshow{\}) # examplesIf}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_esi.R
\name{get_esi}
\alias{get_esi}
\alias{get_esi.default}
\alias{get_esi.sf}
\alias{get_esi.geojson}
\title{Get evaporative stress index (ESI) data}
\usage{
get_esi(object, dates, operation = 5, period = 1, ...)

\method{get_esi}{default}(object, dates, operation = 5, period = 1, ...)

\method{get_esi}{sf}(object, dates, operation = 5, period = 1, as.sf = FALSE, ...)

\method{get_esi}{geojson}(object, dates, operation = 5, period = 1, as.geojson = FALSE, ...)
}
\arguments{
\item{object}{input, an object of class \code{\link[base]{data.frame}} (or
any other object that can be coerced to data.frame), \code{\link[terra]{SpatVector}}, 
\code{\link[terra]{SpatRaster}}, \code{\link[terra]{SpatExtent}},
\code{\link[sf]{sf}} or \code{geojson}}

\item{dates}{a character of start and end dates in that order in the format
"YYYY-MM-DD"}

\item{operation}{optional, an integer that represents which type of
statistical operation to perform on the dataset}

\item{period}{an integer value for the period of ESI data, 
four weeks period = 1, twelve weeks = 2}

\item{...}{additional arguments passed to \code{\link[terra]{terra}} 
or \code{\link[sf]{sf}} methods
See details}

\item{as.sf}{logical, returns an object of class \code{\link[sf]{sf}}}

\item{as.geojson}{logical, returns an object of class \code{geojson}}
}
\value{
A data frame of \acronym{ESI} data:
\item{id}{the index for the rows in \code{object}}
\item{dates}{the dates from which ESI was requested}
\item{lon}{the longitude as provided in \code{object}}
\item{lat}{the latitude as provided in \code{object}}
\item{esi}{the ESI value}
}
\description{
Get evaporative stress index (\acronym{ESI}) from \acronym{SERVIR} Global
via ClimateSERV \acronym{API} Client. \acronym{ESI} is available every four
 (or twelve) weeks from 2001 to present.
The dataset may contain cloudy data which is returned as \code{NA}s.
ClimateSERV works with geojson of type 'Polygon'. The input \code{object} is
 then transformed into polygons with a small buffer area around the point.
}
\details{
\bold{operation}: supported operations are:  
 \tabular{rll}{
 \bold{operation}      \tab    \tab \bold{value}\cr
 max                   \tab =  \tab 0\cr
 min                   \tab =  \tab 1\cr
 median                \tab =  \tab 2\cr
 sum                   \tab =  \tab 4\cr
 average               \tab =  \tab 5 (\emph{default value})\cr
 }

\bold{dist}: numeric, buffer distance for each \code{object} coordinate

\bold{nQuadSegs}: integer, number of segments per buffer quadrant
}
\note{
get_esi may return some warning messages given by 
\code{\link[sf]{sf}}, please look sf documentation for 
possible issues.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

lonlat <- data.frame(lon = c(-55.0281,-54.9857),
                     lat = c(-2.8094, -2.8756))

dates <- c("2017-12-15","2018-06-20")

# by default the function set a very small buffer around the points
# which can return NAs due to cloudiness in ESI data

dt <- get_esi(lonlat, dates = dates)

# the argument dist passed through sf increase the buffer area

dt <- get_esi(lonlat, dates = dates, dist = 0.1)
\dontshow{\}) # examplesIf}
}
