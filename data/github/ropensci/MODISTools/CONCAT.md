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

# MODISTools <a href='https://github.com/ropensci/MODISTools'><img src='https://raw.githubusercontent.com/ropensci/MODISTools/master/MODISTools-logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/MODISTools/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/MODISTools/actions)
[![codecov](https://codecov.io/gh/ropensci/MODISTools/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/MODISTools)
![Status](https://www.r-pkg.org/badges/version/MODISTools) [![rOpenSci
Peer
Review](https://badges.ropensci.org/246_status.svg)](https://github.com/ropensci/software-review/issues/246)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/MODISTools)
<!-- badges: end -->

Programmatic interface to the [‘MODIS Land Products Subsets’ web
services](https://modis.ornl.gov/data/modis_webservice.html). Allows for
easy downloads of [‘MODIS’](http://modis.gsfc.nasa.gov/) time series
directly to your R workspace or your computer. When using the package
please cite the manuscript as referenced below. Keep in mind that the
original manuscript describes versions prior to release 1.0 of the
package. Functions described in this manuscript do not exist in the
current package, please consult [the
documentation](https://docs.ropensci.org/MODISTools/reference/index.html)
to find matching functionality.

## Installation

### stable release

To install the current stable release use a CRAN repository:

``` r
install.packages("MODISTools")
library("MODISTools")
```

### development release

To install the development releases of the package run the following
commands:

``` r
if(!require(devtools)){install.package("devtools")}
devtools::install_github("khufkens/MODISTools")
library("MODISTools")
```

Vignettes are not rendered by default, if you want to include additional
documentation please use:

``` r
if(!require(devtools)){install.package("devtools")}
devtools::install_github("khufkens/MODISTools", build_vignettes = TRUE)
library("MODISTools")
```

## Use

### Downloading MODIS time series

To extract a time series of modis data for a given location and its
direct environment use the mt_subset() function.

<details>
<summary>
detailed parameter description (click to expand)
</summary>
<p>

| Parameter | Description                                                                                                                     |
|-----------|---------------------------------------------------------------------------------------------------------------------------------|
| product   | a MODIS product                                                                                                                 |
| band      | a MODIS product band (if NULL all bands are downloaded)                                                                         |
| lat       | latitude of the site                                                                                                            |
| lon       | longitude of the site                                                                                                           |
| start     | start year of the time series (data start in 1980)                                                                              |
| end       | end year of the time series (current year - 2 years, use force = TRUE to override)                                              |
| internal  | logical, TRUE or FALSE, if true data is imported into R workspace otherwise it is downloaded into the current working directory |
| out_dir   | path where to store the data when not used internally, defaults to tempdir()                                                    |
| km_lr     | force “out of temporal range” downloads (integer)                                                                               |
| km_ab     | suppress the verbose output (integer)                                                                                           |
| site_name | a site identifier                                                                                                               |
| site_id   | a site_id for predefined locations (not required)                                                                               |
| progress  | logical, TRUE or FALSE (show download progress)                                                                                 |

</p>
</details>

``` r
# load the library
library(MODISTools)

# download data
subset <- mt_subset(product = "MOD11A2",
                    lat = 40,
                    lon = -110,
                    band = "LST_Day_1km",
                    start = "2004-01-01",
                    end = "2004-02-01",
                    km_lr = 1,
                    km_ab = 1,
                    site_name = "testsite",
                    internal = TRUE,
                    progress = FALSE)
print(str(subset))
#> 'data.frame':    36 obs. of  21 variables:
#>  $ xllcorner    : chr  "-9370962.97" "-9370962.97" "-9370962.97" "-9370962.97" ...
#>  $ yllcorner    : chr  "4446875.49" "4446875.49" "4446875.49" "4446875.49" ...
#>  $ cellsize     : chr  "926.625433055834" "926.625433055834" "926.625433055834" "926.625433055834" ...
#>  $ nrows        : int  3 3 3 3 3 3 3 3 3 3 ...
#>  $ ncols        : int  3 3 3 3 3 3 3 3 3 3 ...
#>  $ band         : chr  "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" ...
#>  $ units        : chr  "Kelvin" "Kelvin" "Kelvin" "Kelvin" ...
#>  $ scale        : chr  "0.02" "0.02" "0.02" "0.02" ...
#>  $ latitude     : num  40 40 40 40 40 40 40 40 40 40 ...
#>  $ longitude    : num  -110 -110 -110 -110 -110 -110 -110 -110 -110 -110 ...
#>  $ site         : chr  "testsite" "testsite" "testsite" "testsite" ...
#>  $ product      : chr  "MOD11A2" "MOD11A2" "MOD11A2" "MOD11A2" ...
#>  $ start        : chr  "2004-01-01" "2004-01-01" "2004-01-01" "2004-01-01" ...
#>  $ end          : chr  "2004-02-01" "2004-02-01" "2004-02-01" "2004-02-01" ...
#>  $ complete     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ modis_date   : chr  "A2004001" "A2004009" "A2004017" "A2004025" ...
#>  $ calendar_date: chr  "2004-01-01" "2004-01-09" "2004-01-17" "2004-01-25" ...
#>  $ tile         : chr  "h09v05" "h09v05" "h09v05" "h09v05" ...
#>  $ proc_date    : chr  "2015212185706" "2015212201022" "2015212213103" "2015213005429" ...
#>  $ pixel        : int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ value        : int  13135 13120 13350 13354 13123 13100 13324 13331 13098 13069 ...
#> NULL
```

The output format is a *tidy* data frame, as shown above. When witten to
a csv with the parameter `internal = FALSE` this will result in a flat
file on disk.

Note that when a a region is defined using km_lr and km_ab multiple
pixels might be returned. These are indexed using the `pixel` column in
the data frame containing the time series data. The remote sensing
values are listed in the `value` column. When no band is specified all
bands of a given product are returned, be mindful of the fact that
different bands might require different multipliers to represent their
true values. To list all available products, bands for particular
products and temporal coverage see function descriptions below.

### Batch downloading MODIS time series

When a large selection of locations is needed you might benefit from
using the batch download function `mt_batch_subset()`, which provides a
wrapper around the `mt_subset()` function in order to speed up large
download batches. This function has a similar syntax to `mt_subset()`
but requires a data frame defining site names (site_name) and locations
(lat / lon) (or a comma delimited file with the same structure) to
specify a list of download locations.

Below an example is provided on how to batch download data for a data
frame of given site names and locations (lat / lon).

``` r
# create data frame with a site_name, lat and lon column
# holding the respective names of sites and their location
df <- data.frame("site_name" = paste("test",1:2))
df$lat <- 40
df$lon <- -110
  
# test batch download
subsets <- mt_batch_subset(df = df,
                     product = "MOD11A2",
                     band = "LST_Day_1km",
                     internal = TRUE,
                     start = "2004-01-01",
                     end = "2004-02-01")

print(str(subsets))
#> 'data.frame':    8 obs. of  21 variables:
#>  $ xllcorner    : chr  "-9370036.39" "-9370036.39" "-9370036.39" "-9370036.39" ...
#>  $ yllcorner    : chr  "4447802.08" "4447802.08" "4447802.08" "4447802.08" ...
#>  $ cellsize     : chr  "926.625433055834" "926.625433055834" "926.625433055834" "926.625433055834" ...
#>  $ nrows        : int  1 1 1 1 1 1 1 1
#>  $ ncols        : int  1 1 1 1 1 1 1 1
#>  $ band         : chr  "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" ...
#>  $ units        : chr  "Kelvin" "Kelvin" "Kelvin" "Kelvin" ...
#>  $ scale        : chr  "0.02" "0.02" "0.02" "0.02" ...
#>  $ latitude     : num  40 40 40 40 40 40 40 40
#>  $ longitude    : num  -110 -110 -110 -110 -110 -110 -110 -110
#>  $ site         : chr  "test 1" "test 1" "test 1" "test 1" ...
#>  $ product      : chr  "MOD11A2" "MOD11A2" "MOD11A2" "MOD11A2" ...
#>  $ start        : chr  "2004-01-01" "2004-01-01" "2004-01-01" "2004-01-01" ...
#>  $ end          : chr  "2004-02-01" "2004-02-01" "2004-02-01" "2004-02-01" ...
#>  $ complete     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ modis_date   : chr  "A2004001" "A2004009" "A2004017" "A2004025" ...
#>  $ calendar_date: chr  "2004-01-01" "2004-01-09" "2004-01-17" "2004-01-25" ...
#>  $ tile         : chr  "h09v05" "h09v05" "h09v05" "h09v05" ...
#>  $ proc_date    : chr  "2015212185706" "2015212201022" "2015212213103" "2015213005429" ...
#>  $ pixel        : int  1 1 1 1 1 1 1 1
#>  $ value        : int  13098 13062 13297 13323 13098 13062 13297 13323
#> NULL
```

### Listing products

To list all available products use the mt_products() function.

``` r
products <- mt_products()
head(products)
#>        product
#> 1       Daymet
#> 2 ECO4ESIPTJPL
#> 3      ECO4WUE
#> 4       GEDI03
#> 5      MCD12Q1
#> 6      MCD12Q2
#>                                                                       description
#> 1 Daily Surface Weather Data (Daymet) on a 1-km Grid for North America, Version 4
#> 2            ECOSTRESS Evaporative Stress Index PT-JPL (ESI) Daily L4 Global 70 m
#> 3                       ECOSTRESS Water Use Efficiency (WUE) Daily L4 Global 70 m
#> 4             GEDI Gridded Land Surface Metrics (LSM) L3 1km EASE-Grid, Version 2
#> 5           MODIS/Terra+Aqua Land Cover Type (LC) Yearly L3 Global 500 m SIN Grid
#> 6      MODIS/Terra+Aqua Land Cover Dynamics (LCD) Yearly L3 Global 500 m SIN Grid
#>   frequency resolution_meters
#> 1     1 day              1000
#> 2    Varies                70
#> 3    Varies                70
#> 4  One time              1000
#> 5    1 year               500
#> 6    1 year               500
```

### Listing bands

To list all available bands for a given product use the mt_bands()
function.

``` r
bands <- mt_bands(product = "MOD11A2")
head(bands)
#>               band                          description valid_range fill_value
#> 1   Clear_sky_days               Day clear-sky coverage    1 to 255          0
#> 2 Clear_sky_nights             Night clear-sky coverage    1 to 255          0
#> 3    Day_view_angl View zenith angle of day observation    0 to 130        255
#> 4    Day_view_time        Local time of day observation    0 to 240        255
#> 5          Emis_31                   Band 31 emissivity    1 to 255          0
#> 6          Emis_32                   Band 32 emissivity    1 to 255          0
#>    units scale_factor add_offset
#> 1   <NA>         <NA>       <NA>
#> 2   <NA>         <NA>       <NA>
#> 3 degree            1        -65
#> 4    hrs          0.1          0
#> 5   <NA>        0.002       0.49
#> 6   <NA>        0.002       0.49
```

### listing dates

To list all available dates (temporal coverage) for a given product and
location use the mt_dates() function.

``` r
dates <- mt_dates(product = "MOD11A2", lat = 42, lon = -110)
head(dates)
#>   modis_date calendar_date
#> 1   A2000049    2000-02-18
#> 2   A2000057    2000-02-26
#> 3   A2000065    2000-03-05
#> 4   A2000073    2000-03-13
#> 5   A2000081    2000-03-21
#> 6   A2000089    2000-03-29
```

## References

Tuck et al. (2014). [MODISTools - downloading and processing MODIS
remotely sensed data in R Ecology &
Evolution](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.1273),
4(24), 4658 - 4668.

## Acknowledgements

Original development was supported by the UK Natural Environment
Research Council (NERC; grants NE/K500811/1 and NE/J011193/1), and the
Hans Rausing Scholarship. Refactoring was supported through the Belgian
Science Policy office COBECORE project (BELSPO; grant
BR/175/A3/COBECORE). Logo design elements are taken from the FontAwesome
library according to [these terms](https://fontawesome.com/license),
where the globe element was inverted and intersected.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# MODISTools 1.1.2

* removing parallel processing in batch tool due to timeout issues

# MODISTools 1.1.1

* added raster conversion function `mt_to_raster()`
* force integer values on buffer values
* removed keywords

# MODISTools 1.1.0

* ropensci code review
* replaced NULL parameters with missing()
* output a tidy data frame
* include coordinate transform functions to sf objects
* pkgdown documentation generation

# MODISTools 1.0.0

* ropensci refactoring
* version increment follows original package
* push to CRAN

# MODISTools 0.1.0

* first draft
I have read and agree to the the CRAN policies at
http://cran.r-project.org/web/packages/policies.html

## test environments, local, CI and r-hub

- local OSX / Ubuntu 16.04 install on R 3.6.3
- Ubuntu 16.04 on Travis-CI (devel / release)
- rhub check on windows OK
- codecove.io code coverage at ~90%

## local / Travis CI R CMD check results

0 errors | 0 warnings | 0 notes
---
output: github_document
---
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# MODISTools <a href='https://github.com/ropensci/MODISTools'><img src='https://raw.githubusercontent.com/ropensci/MODISTools/master/MODISTools-logo.png' align="right" height="139" /></a>

<!-- badges: start -->
[![R build status](https://github.com/ropensci/MODISTools/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/MODISTools/actions) [![codecov](https://codecov.io/gh/ropensci/MODISTools/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/MODISTools) ![Status](https://www.r-pkg.org/badges/version/MODISTools)
[![rOpenSci Peer Review](https://badges.ropensci.org/246_status.svg)](https://github.com/ropensci/software-review/issues/246)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/MODISTools)
<!-- badges: end -->

Programmatic interface to the ['MODIS Land Products Subsets' web services](https://modis.ornl.gov/data/modis_webservice.html). Allows for easy downloads of ['MODIS'](http://modis.gsfc.nasa.gov/) time series directly to your R workspace or your computer. When using the package please cite the manuscript as referenced below. Keep in mind that the original manuscript describes versions prior to release 1.0 of the package. Functions described in this manuscript do not exist in the current package, please consult [the documentation](https://docs.ropensci.org/MODISTools/reference/index.html) to find matching functionality.

## Installation

### stable release

To install the current stable release use a CRAN repository:

```{r eval = FALSE}
install.packages("MODISTools")
library("MODISTools")
```

### development release

To install the development releases of the package run the following commands:

```{r eval = FALSE}
if(!require(devtools)){install.package("devtools")}
devtools::install_github("khufkens/MODISTools")
library("MODISTools")
```

Vignettes are not rendered by default, if you want to include additional documentation please use:

```{r eval = FALSE}
if(!require(devtools)){install.package("devtools")}
devtools::install_github("khufkens/MODISTools", build_vignettes = TRUE)
library("MODISTools")
```

## Use

### Downloading MODIS time series

To extract a time series of modis data for a given location and its direct environment use the mt_subset() function.

<details><summary>detailed parameter description (click to expand)</summary>
<p>

Parameter     | Description                      
------------- | ------------------------------ 	
product	      | a MODIS product
band	      | a MODIS product band (if NULL all bands are downloaded)
lat           | latitude of the site
lon           | longitude of the site
start      | start year of the time series (data start in 1980)
end        | end year of the time series (current year - 2 years, use force = TRUE to override)
internal      | logical, TRUE or FALSE, if true data is imported into R workspace otherwise it is downloaded into the current working directory
out_dir | path where to store the data when not used internally, defaults to tempdir()
km_lr | force "out of temporal range" downloads (integer)
km_ab | suppress the verbose output (integer)
site_name | a site identifier
site_id | a site_id for predefined locations (not required)
progress | logical, TRUE or FALSE (show download progress)

</p>
</details>

```{r eval = TRUE}
# load the library
library(MODISTools)

# download data
subset <- mt_subset(product = "MOD11A2",
                    lat = 40,
                    lon = -110,
                    band = "LST_Day_1km",
                    start = "2004-01-01",
                    end = "2004-02-01",
                    km_lr = 1,
                    km_ab = 1,
                    site_name = "testsite",
                    internal = TRUE,
                    progress = FALSE)
print(str(subset))
```

The output format is a *tidy* data frame, as shown above. When witten to a csv with the parameter `internal = FALSE` this will result in a flat file on disk.

Note that when a a region is defined using km_lr and km_ab multiple pixels might be returned. These are indexed using the `pixel` column in the data frame containing the time series data. The remote sensing values are listed in the `value` column. When no band is specified all bands of a given product are returned, be mindful of the fact that different bands might require different multipliers to represent their true values. To list all available products, bands for particular products and temporal coverage see function descriptions below.

### Batch downloading MODIS time series

When a large selection of locations is needed you might benefit from using the batch download function `mt_batch_subset()`, which provides a wrapper around the `mt_subset()` function in order to speed up large download batches. This function has a similar syntax to `mt_subset()` but requires a data frame defining site names (site_name) and locations (lat / lon) (or a comma delimited file with the same structure) to specify a list of download locations.

Below an example is provided on how to batch download data for a data frame of given site names and locations (lat / lon).

```{r eval = TRUE}
# create data frame with a site_name, lat and lon column
# holding the respective names of sites and their location
df <- data.frame("site_name" = paste("test",1:2))
df$lat <- 40
df$lon <- -110
  
# test batch download
subsets <- mt_batch_subset(df = df,
                     product = "MOD11A2",
                     band = "LST_Day_1km",
                     internal = TRUE,
                     start = "2004-01-01",
                     end = "2004-02-01")

print(str(subsets))
```

### Listing products

To list all available products use the mt_products() function.

```{r eval = TRUE}
products <- mt_products()
head(products)
```

### Listing bands

To list all available bands for a given product use the mt_bands() function.

```{r eval = TRUE}
bands <- mt_bands(product = "MOD11A2")
head(bands)
```

### listing dates

To list all available dates (temporal coverage) for a given product and location use the mt_dates() function.

```{r eval = TRUE}
dates <- mt_dates(product = "MOD11A2", lat = 42, lon = -110)
head(dates)
```

## References

Tuck et al. (2014). [MODISTools - downloading and processing MODIS remotely sensed data in R Ecology & Evolution](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.1273), 4(24), 4658 - 4668.

## Acknowledgements

Original development was supported by the UK Natural Environment Research Council (NERC; grants NE/K500811/1 and NE/J011193/1), and the Hans Rausing Scholarship. Refactoring was supported through the Belgian Science Policy office COBECORE project (BELSPO; grant BR/175/A3/COBECORE). Logo design elements are taken from the FontAwesome library according to [these terms](https://fontawesome.com/license), where the globe element was inverted and intersected.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "MODISTools"
author: "Koen Hufkens"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{MODISTools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# load the library
library(MODISTools)
library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)
library(sf)

# pre-load data
data("arcachon_lai")
data("arcachon_lc")

```

The MODISTools package has as goal to facilitate the interface between R and the MODIS Land Product Subset API at the Oak Ridge National Laboratory Distributed Active Archive Center (DAAC). This programmatic interface to the ['MODIS Land Products Subsets' web services](https://modis.ornl.gov/data/modis_webservice.html) allows for easy downloads of 'MODIS' time series (of single pixels or small regions of interest) directly to your R workspace or your computer. Below an example is provided on how to download a MODIS time series as well as list ancillary data.

# Listing products / bands / dates

In order to assess which products are available, which product bands are provided and which temporal range is covered one has to list these ancillary data. All these options can be queried using the `mt_*()` functions.

To list all available products use the `mt_products()` function.

```{r eval = TRUE}
products <- mt_products()
head(products)
```

Once you have settled on the product you want to use in your analysis you will have to select a band, or multiple bands you want to download for a given location. To list all available bands for a given product use the mt_bands() function. You can use the `mt_bands()` function to list all available bands for a given product. Below I list all bands for the MOD13Q1 vegetation index product.

```{r eval = TRUE}
bands <- mt_bands(product = "MOD13Q1")
head(bands)
```

> Note the band names (listed in band column) you want to download, this variable will need to be passed in the final download statement.

Similarly you can list all available dates (temporal coverage) for a given product and location defined using a latitude and longitude with the `mt_dates()` function.

```{r eval = TRUE}
dates <- mt_dates(product = "MOD13Q1", lat = 42, lon = -110)
head(dates)
```

# Downloading MODIS time series

Once you decide on which data to download using the above functions you can use these parameters to download a time series using the `mt_subset()` function. The below query downloads MOD15A2H based leaf area index (LAI) data for the year 2014 for an area around the Arcachon basin in the south west of France. We will also download land cover data (MCD12Q1, IGBP) at a similar scale. The location is named 'arcachon'. The output is saved to a variables called subset and LC in the R workspace (as defined by the parameter internal = TRUE, when set to FALSE the data is written to file). Keep in mind that this operation might take a while.

```{r eval = FALSE}
# download the MODIS land cover (IGBP) and NDVI data
# for a region around the French city and basin of Arcachon
arcachon_lai <- mt_subset(product = "MOD15A2H",
                    lat = 44.656286,
                    lon =  -1.174748,
                    band = "Lai_500m",
                    start = "2004-01-01",
                    end = "2004-12-30",
                    km_lr = 20,
                    km_ab = 20,
                    site_name = "arcachon",
                    internal = TRUE,
                    progress = FALSE)

arcachon_lc <- mt_subset(product = "MCD12Q1",
  lat = 44.656286,
  lon =  -1.174748,
  band = "LC_Type1",
  start = "2004-01-01",
  end = "2004-3-20",
  km_lr = 20,
  km_ab = 20,
  site_name = "arcachon",
  internal = TRUE,
  progress = FALSE)
```

The output format is a *tidy* data frame, as shown above. When witten to a csv with the parameter `internal = FALSE` this will result in a flat file on disk.

```{r}
head(arcachon_lai)
head(arcachon_lc)
```

Note that when a a region is defined using km_lr and km_ab multiple pixels might be returned. These are indexed using the `pixel` column in the data frame containing the time series data. The remote sensing values are listed in the `value` column. When no band is specified all bands of a given product are returned, be mindful of the fact that different bands might require different multipliers to represent their true values. 

When a large selection of locations is needed you might benefit from using the batch download function `mt_batch_subset()`, which provides a wrapper around the `mt_subset()` function in order to speed up large download batches. This function has a similar syntax to `mt_subset()` but requires a data frame defining site names (site_name) and locations (lat / lon) (or a comma delimited file with the same structure) to specify a list of download locations.

```{r eval = TRUE}
# create data frame with a site_name, lat and lon column
# holding the respective names of sites and their location
df <- data.frame("site_name" = paste("test",1:2), stringsAsFactors = FALSE)
df$lat <- 40
df$lon <- -110

# an example batch download data frame
head(df)
```

```{r eval = FALSE}
# test batch download
subsets <- mt_batch_subset(df = df,
                     product = "MOD13Q1",
                     band = "250m_16_days_NDVI",
                     km_lr = 1,
                     km_ab = 1,
                     start = "2004-01-01",
                     end = "2004-12-30",
                     internal = TRUE)
```

# Worked example using LAI values around the bay of Arcachon

The below example processes the data downloaded above to look at differences in the seasonal changes in leaf area index (LAI, or the amount of leaves per unit ground area) for the Arcachon bay in south-west France. To do this we merge the land cover and LAI data on a pixel by pixel basis.

```{r}
# merge land cover and lai data
arcachon <- arcachon_lc %>%
  rename("lc" = "value") %>%
  select("lc","pixel") %>%
  right_join(arcachon_lai, by = "pixel")
```

Then, filter out all non valid values (> 100), only select evergreen and deciduous land cover classes (1 and 5, or, ENF and DBF respectivelly), convert them to more readable labels, and across these land cover classes take the median per acquisition date.

```{r}
# create a plot of the data - accounting for the multiplier (scale) component
arcachon <- arcachon %>%
  filter(value <= 100,
         lc %in% c("1","5")) %>% # retain everything but fill values
  mutate(lc = ifelse(lc == 1, "ENF","DBF")) %>%
  group_by(lc, calendar_date) %>% # group by lc and date
  summarize(doy = as.numeric(format(as.Date(calendar_date)[1],"%j")),
            lai_mean = median(value * as.double(scale)))
```

Finally, the plot will show you the seasonal time series of LAI for both land cover classes (ENF and DBF). Note the difference in timing and amplitude between both these forest types, where the evergreen (ENF) pixels show lower LAI values and a more gradual seasonal pattern compared to the deciduous trees.

```{r fig.width = 7, fig.height=3}
# plot LAI by date and per land cover class
ggplot(arcachon, aes(x = doy, y = lai_mean)) +
  geom_point() +
  geom_smooth(span = 0.3, method = "loess") +
  labs(x = "day of year (DOY)",
       y = "leaf area index (LAI)") +
  theme_minimal() +
  facet_wrap(~ lc)
```

# Conversion of corner coordinates

Corner coordinates of the pixel area extracted are provided, these can be used to calculate the coverage of the extracted area. Coordinates are provided in the original sinusoidal grid coordinates and first have to be transformed into latitude longitude (for convenience).

```{r }
# convert the coordinates
lat_lon <- sin_to_ll(arcachon_lc$xllcorner, arcachon_lc$yllcorner)

# bind with the original dataframe
subset <- cbind(arcachon_lc, lat_lon)

head(subset)
```

Together with meta-data regarding cell size, number of columns and rows the bounding box of the extracted data can be calculated.

```{r fig.width = 5, fig.height=5}
# convert to bounding box
bb <- apply(arcachon_lc, 1, function(x){
  mt_bbox(xllcorner = x['xllcorner'],
          yllcorner = x['yllcorner'],
           cellsize = x['cellsize'],
           nrows = x['nrows'],
           ncols = x['ncols'])
})

# plot one bounding box
plot(bb[[1]])

# add the location of the queried coordinate within the polygon
points(arcachon_lc$longitude[1],
       arcachon_lc$latitude[1],
       pch = 20,
       col = "red")
```
 
# Conversion to (gridded) raster data

Although the package is often used to deal with single pixel locations the provisions to download a small region of interest defined by kilometers left-right (west-east) or top-bottom (north-south) allows you to grab small geographic regions for further analysis. The default tidy dataframe format isn't ideal for visualizing this inherently spatial data. Therefore a helper function `mt_to_raster()` is available to convert the tidy dataframe to a gridded (georeferenced) raster format.

Below a small region (20x20km) is downloaded around the seaside town of Arcachon, France for the [MODIS land cover product (MCD12Q1)](https://lpdaac.usgs.gov/products/mcd12q1v006/). The data is converted using `mt_to_raster()` with a reproject parameter set to true to plot latitude and longitude coordinates (instead of the default sinusoidal ones).


```{r fig.width = 5, fig.height=5}
# convert to raster, when reproject is TRUE
# the data is reprojected to lat / lon if FALSE
# the data is shown in its original sinuidal projection
LC_r <- mt_to_raster(df = arcachon_lc, reproject = TRUE)

# plot the raster data as a map
plot(LC_r)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{arcachon_lai}
\alias{arcachon_lai}
\title{arcachon_lai}
\format{
A MODISTools tidy data frame
}
\usage{
arcachon_lai
}
\description{
MODIS leaf area index (LAI) around the French town of Arcachon
derived from the MODIS MOD15A2H product (band Lai_500m).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{arcachon_lc}
\alias{arcachon_lc}
\title{arcachon_lc}
\format{
A MODISTools tidy data frame
}
\usage{
arcachon_lc
}
\description{
MODIS land cover (IGBP) around the French town of Arcachon
derived from the MODIS MCD12Q2 product (band LC_Type1).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_to_raster.R
\name{mt_to_raster}
\alias{mt_to_raster}
\title{Convert tidy MODISTools data to raster (stack)}
\usage{
mt_to_raster(df, reproject = FALSE)
}
\arguments{
\item{df}{a valid MODISTools data frame with a single band (filter for a
particular band using the dplyr \code{filter()} function or base \code{subset()}}

\item{reproject}{reproject output to lat / long (default = \code{FALSE})}
}
\value{
A raster stack populated with the tidy dataframe values
}
\description{
Convert tidy MODISTools data to a raster (stack)
}
\examples{

\donttest{
# list all available MODIS Land Products Subsets products
# download data
LC <- mt_subset(product = "MCD12Q1",
 lat = 48.383662,
 lon = 2.610250,
 band = "LC_Type1",
 start = "2005-01-01",
 end = "2005-12-30",
 km_lr = 2,
 km_ab = 2,
 site_name = "testsite",
 internal = TRUE,
 progress = FALSE)

head(LC)

# convert to raster
LC_r <- mt_to_raster(df = LC)
}

}
\seealso{
\code{\link[MODISTools]{mt_subset}}
\code{\link[MODISTools]{mt_batch_subset}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_bands.R
\name{mt_bands}
\alias{mt_bands}
\title{Download all available bands}
\usage{
mt_bands(product)
}
\arguments{
\item{product}{a valid MODIS product name}
}
\value{
A data frame of all available bands for a MODIS Land
Products Subsets products
}
\description{
Lists all available bands for a MODIS Land Products Subset product.
}
\examples{

\donttest{
# list all available MODIS Land Products Subsets products
bands <- mt_bands(product = "MCD12Q2")
head(bands)

}

}
\seealso{
\code{\link[MODISTools]{mt_products}}
\code{\link[MODISTools]{mt_sites}} \code{\link[MODISTools]{mt_dates}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_batch_subset.R
\name{mt_batch_subset}
\alias{mt_batch_subset}
\title{Batch download MODIS Land Products subsets}
\usage{
mt_batch_subset(
  df,
  product,
  band,
  start = "2000-01-01",
  end = format(Sys.time(), "\%Y-\%m-\%d"),
  km_lr = 0,
  km_ab = 0,
  out_dir = tempdir(),
  internal = TRUE
)
}
\arguments{
\item{df}{a CSV file or data frame holding locations and their sitenames to
batch process with column names site_name, lat, lon holding the respective
sitenames, latitude and longitude. When providing a CSV make sure that the
data are comma separated.}

\item{product}{a valid MODIS product name}

\item{band}{band to download}

\item{start}{start date}

\item{end}{end date}

\item{km_lr}{km left-right to sample}

\item{km_ab}{km above-below to sample}

\item{out_dir}{location where to store all data}

\item{internal}{should the data be returned as an internal data structure
\code{TRUE} or \code{FALSE} (default = \code{TRUE})}
}
\value{
A data frame combining meta-data and actual data values, data from
different sites is concatenated into one large dataframe. Subsets can be
created by searching on sitename.
}
\description{
Lists all available dates for a MODIS Land Products Subset product
at a particular location.
}
\examples{

\dontrun{
# create data frame with a site_name, lat and lon column
# holding the respective names of sites and their location
df <- data.frame("site_name" = paste("test",1:2))
df$lat <- 40
df$lon <- -110

print(df)

# test batch download
subsets <- mt_batch_subset(df = df,
                        product = "MOD11A2",
                        band = "LST_Day_1km",
                        internal = TRUE,
                        start = "2004-01-01",
                        end = "2004-03-31")

# the same can be done using a CSV file with
# a data structure similar to the dataframe above

write.table(df, file.path(tempdir(),"my_sites.csv"),
 quote = FALSE,
 row.names = FALSE,
 col.names = TRUE,
 sep = ",")

# test batch download form CSV
subsets <- mt_batch_subset(df = file.path(tempdir(),"my_sites.csv"),
                        product = "MOD11A2",
                        band = "LST_Day_1km",
                        internal = TRUE,
                        start = "2004-01-01",
                        end = "2004-03-31"
                        )

head(subsets)
}
}
\seealso{
\code{\link[MODISTools]{mt_sites}}
\code{\link[MODISTools]{mt_dates}} \code{\link[MODISTools]{mt_bands}}
\code{\link[MODISTools]{mt_products}}
\code{\link[MODISTools]{mt_subset}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate_conversion.R
\name{sin_to_ll}
\alias{sin_to_ll}
\title{Convert sinusoidal coordinates to lat / lon}
\usage{
sin_to_ll(x, y)
}
\arguments{
\item{x}{sinusoidal x coordinate (vector)}

\item{y}{sinusoidal y coordinate (vector)}
}
\description{
A full description of the sinusoidal projection is provided on the
lpdaac page:
https://lpdaac.usgs.gov/dataset_discovery/modis
and wikipedia:
https://en.wikipedia.org/wiki/Sinusoidal_projection
}
\examples{

\donttest{
# Download some test data
subset <- mt_subset(product = "MOD11A2",
                        lat = 40,
                        lon = -110,
                        band = "LST_Day_1km",
                        start = "2004-01-01",
                        end = "2004-03-31",
                        progress = FALSE)

# convert sinusoidal to lat / lon
lat_lon <- sin_to_ll(subset$xllcorner, subset$yllcorner)

# bind with the original dataframe
subset <- cbind(subset, lat_lon)
head(subset)
}
}
\seealso{
\code{\link[MODISTools]{mt_bbox}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_subset.R
\name{mt_subset}
\alias{mt_subset}
\title{Download MODIS Land Products subsets}
\usage{
mt_subset(
  product,
  band,
  lat,
  lon,
  start = "2000-01-01",
  end = format(Sys.time(), "\%Y-\%m-\%d"),
  km_lr = 0,
  km_ab = 0,
  site_id,
  network,
  site_name = "sitename",
  out_dir = tempdir(),
  internal = TRUE,
  progress = TRUE
)
}
\arguments{
\item{product}{a valid MODIS product name}

\item{band}{band or bands (as a character vector) to download}

\item{lat}{latitude in decimal degrees}

\item{lon}{longitude in decimal degrees}

\item{start}{start date}

\item{end}{end date}

\item{km_lr}{km left-right to sample (rounded to the nearest integer)}

\item{km_ab}{km above-below to sample (rounded to the nearest integer)}

\item{site_id}{site id (overides lat / lon)}

\item{network}{the network for which to generate the site list,
when not provided the complete list is provided}

\item{site_name}{arbitrary site name used in writing data to file
(default = sitename)}

\item{out_dir}{path where to store the data if writing to disk
(default = tempdir())}

\item{internal}{should the data be returned as an internal data structure
\code{TRUE} or \code{FALSE} (default = \code{TRUE})}

\item{progress}{show download progress}
}
\value{
A data frame combining meta-data and actual data values.
}
\description{
Download a MODIS Land Products Subset product
for a given point location buffered with a given amount of kilometers
left-right, top-bottom for a given location (provided as latitude and
longitude values).
}
\examples{

\donttest{
# list all available MODIS Land Products Subsets products
# download data
subset <- mt_subset(product = "MOD11A2",
                        lat = 40,
                        lon = -110,
                        band = "LST_Day_1km",
                        start = "2004-01-01",
                        end = "2004-03-31",
                        progress = FALSE)
 head(subset)
}
}
\seealso{
\code{\link[MODISTools]{mt_sites}}
\code{\link[MODISTools]{mt_dates}} \code{\link[MODISTools]{mt_bands}}
\code{\link[MODISTools]{mt_products}}
\code{\link[MODISTools]{mt_batch_subset}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_products.R
\name{mt_products}
\alias{mt_products}
\title{Download all available products}
\usage{
mt_products()
}
\value{
A data frame of all available MODIS Land Products Subsets products
}
\description{
Lists all available MODIS Land Products Subset products.
}
\examples{

\donttest{
# list all available MODIS Land Products Subsets products
products <- mt_products()
head(products)
}

}
\seealso{
\code{\link[MODISTools]{mt_bands}}
\code{\link[MODISTools]{mt_sites}} \code{\link[MODISTools]{mt_dates}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_sites.R
\name{mt_sites}
\alias{mt_sites}
\title{Download all available fixed sites}
\usage{
mt_sites(network)
}
\arguments{
\item{network}{the network for which to generate the site list,
when not provided the complete list is provided}
}
\value{
A data frame of all available MODIS Land Products Subsets
pre-processed sites
}
\description{
Lists all available MODIS Land Products Subset pre-processed sites
}
\examples{

\donttest{
# list all available MODIS Land Products Subsets products
sites <- mt_sites()
print(head(sites))
}

}
\seealso{
\code{\link[MODISTools]{mt_products}}
\code{\link[MODISTools]{mt_bands}} \code{\link[MODISTools]{mt_dates}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_dates.R
\name{mt_dates}
\alias{mt_dates}
\title{Download all available dates}
\usage{
mt_dates(product, lat, lon, site_id, network)
}
\arguments{
\item{product}{a valid MODIS product name}

\item{lat}{latitude in decimal degrees}

\item{lon}{longitude in decimal degrees}

\item{site_id}{site id (overides lat / lon)}

\item{network}{the network for which to generate the site list,
when not provided the complete list is provided}
}
\value{
A data frame of all available dates for a MODIS Land
Products Subsets products at the given location.
}
\description{
Lists all available dates for a MODIS Land Products Subset product
at a particular location.
}
\examples{

\donttest{
# list all available MODIS Land Products Subsets products
bands <- mt_dates(product = "MOD11A2", lat = 40, lon = -110)
head(bands)
}
}
\seealso{
\code{\link[MODISTools]{mt_products}}
\code{\link[MODISTools]{mt_sites}} \code{\link[MODISTools]{mt_bands}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate_conversion.R
\name{mt_bbox}
\alias{mt_bbox}
\title{Converts lower-left sinusoidal coordinates to lat-lon sf bounding box}
\usage{
mt_bbox(xllcorner, yllcorner, cellsize, nrows, ncols, transform = TRUE)
}
\arguments{
\item{xllcorner}{lower left x coordinate as provided by
\code{\link[MODISTools]{mt_subset}}}

\item{yllcorner}{lower left y coordinate as provided by
\code{\link[MODISTools]{mt_subset}}}

\item{cellsize}{cell size provided by \code{\link[MODISTools]{mt_subset}}}

\item{nrows}{cell size provided by \code{\link[MODISTools]{mt_subset}}}

\item{ncols}{cell size provided by \code{\link[MODISTools]{mt_subset}}}

\item{transform}{transform the bounding box from sin to lat long coordinates,
\code{TRUE} or \code{FALSE} (default = \code{TRUE})}
}
\description{
Converts lower-left sinusoidal coordinates to lat-lon sf bounding box
}
\examples{

\donttest{
# Download some test data
subset <- mt_subset(product = "MOD11A2",
                        lat = 40,
                        lon = -110,
                        band = "LST_Day_1km",
                        start = "2004-01-01",
                        end = "2004-03-31",
                        progress = FALSE)

# convert sinusoidal to lat / lon
lat_lon <- sin_to_ll(subset$xllcorner, subset$yllcorner)

# bind with the original dataframe
subset <- cbind(subset, lat_lon)

# convert to bounding box
bb <- apply(subset, 1, function(x){
  mt_bbox(xllcorner = x['xllcorner'],
          yllcorner = x['yllcorner'],
          cellsize = x['cellsize'],
          nrows = x['nrows'],
          ncols = x['ncols'])
})

head(bb)
}
}
\seealso{
\code{\link[MODISTools]{sin_to_ll}},
\code{\link[MODISTools]{mt_subset}}
}
