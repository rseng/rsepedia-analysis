# rsat

<!-- badges: start -->

[![CRAN version](https://www.r-pkg.org/badges/version/rsat)](https://CRAN.R-project.org/package=rsat) 
[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/437_status.svg)](https://github.com/ropensci/software-review/issues/437) 
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) 
[![Codecov test coverage](https://codecov.io/gh/ropensci/rsat/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/rsat?branch=master)
[![R build status](https://github.com/ropensci/rsat/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rsat/actions)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/rsat)](https://CRAN.R-project.org/package=rsat)
<!-- badges: end -->

The goal of `rsat` is to help you handling time-series of satellite images from multiple platforms in a local, efficient, and standardized way. The package provides tools to;

1.  Search (run `vignette("rsat1_search", package = "rsat")` command)
2.  Download (run `vignette("rsat2_download", package = "rsat")` command)
3.  Customize, and (run `vignette("rsat3_customize", package = "rsat")` command)
4.  Process (run `vignette("rsat4_process", package = "rsat")` command)

satellite images from Landsat, MODIS, and Sentinel for a region and time of interest.

## Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# check and install devtools
if(!require("devtools")){
   install.packages("devtools")
}
# check and install rmarkdown
if(!require("rmarkdown")){
  install.packages("rmarkdown")
}

devtools::install_github("spatialstatisticsupna/rsat", build_vignettes=TRUE)
```

### Linux

In Linux, you need to install additional libraries before starting with `rsat`. Use the following commands for:

-   **Debian/Ubuntu**

<!-- -->

    sudo apt update 
    sudo apt install r-cran-rcpp gdal-bin libgdal-dev libproj-dev openssl libssl-dev xml2 libxml2-dev libmagick++-dev

-   **RedHat/Fedora**

<!-- -->

    sudo dnf install gdal gdal-devel proj-devel xml2 libxml2-devel libcurl-devel openssl-devel ImageMagick-c++-devel R-devel udunits2-devel sqlite-devel geos-devel pandoc

## Log-in profiles

The registration in the following online portals is required to get a full access to satellite images with `rsat`;
-   [USGS](https://ers.cr.usgs.gov/register/) USGS is the sole science agency for the Department of the Interior of United States. Provide access to Modis Images. More information about USGS can be found [Here](https://www.usgs.gov/).
-   [EarthData](https://urs.earthdata.nasa.gov): A repository of NASA's earth observation data-sets. More information about EarthData can be found [here](https://earthdata.nasa.gov/earth-observation-data).
-   [SciHub](https://scihub.copernicus.eu/dhus/#/self-registration), a web service giving access to Copernicus' scientific data hub. Please go [here](https://scihub.copernicus.eu/) to find more details about the data hub.

For convenience, try to use the same username and password for all of them. To satisfy the criteria of all web services make sure that the username is $4$ characters long and includes a period, number or underscore. The password must be $12$ character long and should include characters with at least one capital letter, and numbers.

## Example

This is a basic example which shows you how to compute the Normalized Difference Vegetation Index from a MODIS image captured on January 11th, 2020 in northern Spain (Navarre):

``` r
library(rsat)

# replace with your own "username" and "password"
set_credentials("username", "password")

# region and time of interest: rtoi
roi <- ex.navarre
toi <- as.Date("2020-01-11")
rtp <- tempdir()

set_database(file.path(tempdir(), "DATABASE"))

navarre <- new_rtoi("Navarre", roi, rtp)

# search, acquire, customize, and process
rsat_search(region = navarre, product = "mod09ga", dates = toi)
rsat_download(navarre)
rsat_mosaic(navarre, overwrite = TRUE)

rsat_derive(navarre, 
            product = "mod09ga", 
            variable = "NDVI")

# plot the results
plot(navarre, "view" , 
      product = "mod09ga", 
      variable = "NDVI", 
      breaks = seq(0, 1, 0.1))
      
plot(navarre,"dates")
```

See the vignettes for more examples:

    browseVignettes("rsat")

## Related similar packages

R has become an outstanding tool for remote sensing image analysis. There are several tools for the search and acquisition of satellite images, however, rsat is the first package that standardizes all the procedures in data acquisition to provide an unique workflow for any multispectral satellite.

Currently there are several packages dedicated to remote sensing topic, but they are usually ad-hoc packages for each satellite. Here is a list of some of the most popular R packages dedicated to satellite imagery:

### Multi satellite packages

-   [RGISTools](https://github.com/spatialstatisticsupna/RGISTools)

-   [getSpatialData](https://github.com/16EAGLE/getSpatialData)

-   [luna](https://github.com/rspatial/luna)

The closest package to `rsat` is RGISTools. `rsat` is the redefinition of the RGISTools package reprogrammed from scrach in the object-oriented programming paradigm. Many of the RGISTools code lines have been used to develop `rsat`, but these have been optimized and redundancies in the code have been removed in order to facilitate its maintenance. In addition, `rsat` contains new features and R classes to make it more user-friendly.

`getSpatialData` is another package very similar to `rsat`. The package has the same philosophy of having a single package for searching and downloading satellite images. However, the development of `rsat` goes a bit further and in addition to search and download, the package helps you to organize all the downloaded information in a structured database. `rsat` allows you to use the metadata of the images to see the direct relation with your region of interest before downloading it. Also all image processing standardization is not developed in `getSpatialData`.

The last package dedicated to image downloading is `luna`.Searching and downloading images compared to rsat is a bit more complicated. It is only able to search and download Modis and Landsat images, and does not help you in organizing the image products.

### Single satellite packages

-   [rLandsat](https://github.com/atlanhq/rLandsat)

-   [getLandsat](https://github.com/ropensci/getlandsat)

-   [sen2r](https://github.com/ranghetti/sen2r)

`rLandsat` makes it easy to search for Landsat8 product IDs, place an order on USGS-ESPA and download the data. `rsat` on the other hand is able to do the image search without knowing all the ids, just using a polygon of the region of interest, making the search process much easier.

`getlandsat` provides access to Landsat 8 metadata and images hosted on AWS S3 at. The package only data for the users, and does not help in further use, as rsat does.

`sen2r` is an R library which helps to download and preprocess Sentinel-2 optical images. This is done through a GUI, something that can be very interesting for users but limits the analysis of the information prior to downloading, which can be done with `rsat`.

### Raster processing packages

-   [landsat](https://cran.r-project.org/package=landsat)

-   [satellite](https://github.com/environmentalinformatics-marburg/satellite)

-   [OpenImageR](https://github.com/mlampros/OpenImageR)

-   [RSToolbox](https://github.com/bleutner/RStoolbox)

-   [sits](https://github.com/e-sensing/sits)

`rsat` helps you to search, download and pre-process the images, but once these procedures are done it allows you to extract all the processed information into the most used raster classes in R (`raster`, `stars` or `spatRaster`). The image processing packages can be used for further analysis in these R classes.

## Contributing

We accept contributions to improve the package. Before contributing, please follow these steps:

-    Contributions should be thoroughly tested with testthat.
-    Code style should attempt to follow the tidyverse style guide.
-    Please attempt to describe what you want to do prior to contributing by submitting an issue.
-    Please follow the typical github fork - pull-request workflow.
-    Make sure you use roxygen and run Check before contributing. More on this front as the package matures.


## Citation

``` r
citation("rsat")[1]
```

To cite the package:

U. Pérez-Goya, M. Montesino-SanMartin, A F Militino, M D Ugarte (2021). rsat: Dealing with Multiplatform Satellite Images from Landsat, MODIS, and Sentinel. R package version 0.1.16. <https://github.com/ropensci/rsat>.

## Acknowledgements

This work has been financed by projects MTM2017-82553-R (AEI/FEDER, UE) and PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033.
rsat 0.1.17 (2021-12-13)
=========================
### NEW FEATURES
  * Fixed minor bugs in the code
  * Corrected some examples from the manual
  * Fixed vignettes
  
rsat 0.1.16 (2021-10-06)
=========================
### NEW FEATURES
  * Added support to new landsat products: landsat_ot_c2_l2, 
  landsat_ot_c2_l1, lsr_landsat_8_c1, landsat_8_c1,landsat_etm_c2_l2, 
  landsat_etm_c2_l1, lsr_landsat_etm_c1, landsat_etm_c1, 
  landsat_tm_c2_l2, landsat_tm_c2_l1, landsat_mss_c2_l1, 
  lsr_landsat_tm_tm_c1, landsat_tm_c1, landsat_mss_c1
  
  
rsat 0.1.15 (2021-09-09)
=========================
### NEW FEATURES
  
  * Added database path to global environment with the function 
  `set_database("character")` and `get_database()`
  * Added `get_SpatRaster` function to get `SpatRaster` 
  class from `terra` package 

### MINOR IMPROVEMENTS

  * Improved the examples in the documentation and added new datasets
  * Added vignettes enumeration
  * Improved the performance of search function
  * Fixed the warnings produced by mosaic and derive functions
  * Updated `sf` objects to the last version
  * Use `terra` package instead of `raster` internally
  * Improved stored `rtoi` files
  * Added `rsat` section in the manual

### BUG FIXES

  * Fixed error searching modis images
  * Fixed `unique` function in records
  * Fixed `rsat_list_data`
  * Fixed multiple errors around the code
  
### DEPRECATED AND DEFUNCT

  * Function name change `sat_search()` -> `rsat_search()`
  * Function name change `download()` -> `rsat_download()`
  * Function name change `mosaic()` -> `rsat_mosaic()` 
  * Function name change `derive()` -> `rsat_derive()`
  * Function name change `preview()` -> `rsat_preview()`
  * Function name change `cloud_mask()` -> `rsat_cloudMask()`
  * Function name change `smoothing_images()` -> `rsat_smoothing_images()`
  * Function name change `list_data()` -> `rsat_list_data()`
  * Function name change `get_raster()` -> `rsat_get_raster()`
  * Function name change `get_stars()` -> `rsat_get_stars()`

### DOCUMENTATION

  * Added rsat package section
  * Added examples in several functions with new datasets
  * Improved documentation in general
  * Grouped some functions in one manual entry


rsat 0.1.14 (2021-01-01)
=========================

### NEW FEATURES

  * ROpenSci review version
---
output: github_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# rsat

<!-- badges: start -->

[![CRAN version](https://www.r-pkg.org/badges/version/rsat)](https://CRAN.R-project.org/package=rsat)  
[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/437_status.svg)](https://github.com/ropensci/software-review/issues/437)
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Codecov test coverage](https://codecov.io/gh/ropensci/rsat/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/rsat?branch=master)
[![R build status](https://github.com/ropensci/rsat/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rsat/actions)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/rsat)](https://CRAN.R-project.org/package=rsat)
<!-- badges: end -->

The goal of `rsat` is to help you handling time-series of satellite images from multiple platforms in a local, efficient, and standardized way. The package provides tools to;

1.  Search (run `vignette("rsat1_search", package = "rsat")`)
2.  Download (run `vignette("rsat2_download", package = "rsat")`)
3.  Customize, and (run `vignette("rsat3_customize", package = "rsat")`)
4.  Process (run `vignette("rsat4_process", package = "rsat")`)

satellite images from Landsat, MODIS, and Sentinel for a region and time of interest.

## Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# check and install devtools
if(!require("devtools")){
   install.packages("devtools")
}
# check and install rmarkdown
if(!require("rmarkdown")){
  install.packages("rmarkdown")
}

devtools::install_github("spatialstatisticsupna/rsat", build_vignettes=TRUE)
```

### Linux

In Linux, you need to install additional libraries before starting with `rsat`. Use the following commands for:

-   **Debian/Ubuntu**

<!-- -->

    sudo apt update sudo apt install r-cran-rcpp gdal-bin libgdal-dev libproj-dev openssl openssl-dev xml2 libxml2-dev libmagick++-dev

-   **RedHat/Fedora**

<!-- -->

    sudo dnf install gdal gdal-devel proj-devel xml2 libxml2-devel libcurl-devel openssl-devel ImageMagick-c++-devel R-devel udunits2-devel sqlite-devel geos-devel pandoc

## Log-in profiles

The registration in the following online portals is required to get a full access to satellite images with `rsat`;
-   [USGS](https://ers.cr.usgs.gov/register/) USGS is the sole science agency for the Department of the Interior of United States. Provide access to Modis Images. More information about USGS can be found [Here](https://www.usgs.gov/).
-   [EarthData](https://urs.earthdata.nasa.gov): A repository of NASA's earth observation data-sets. More information about EarthData can be found [here](https://earthdata.nasa.gov/earth-observation-data).
-   [SciHub](https://scihub.copernicus.eu/dhus/#/self-registration), a web service giving access to Copernicus' scientific data hub. Please go [here](https://scihub.copernicus.eu/) to find more details about the data hub.

For convenience, try to use the same username and password for all of them. To satisfy the criteria of all web services make sure that the username is $4$ characters long and includes a period, number or underscore. The password must be $12$ character long and should include characters with at least one capital letter, and numbers.

## Example

This is a basic example which shows you how to compute the Normalized Difference Vegetation Index from a MODIS image captured on January 11th, 2020 in northern Spain (Navarre):

``` r
library(rsat)

# replace with your own "username" and "password"
set_credentials("username", "password")

# region and time of interest: rtoi
roi <- ex.navarre
toi <- as.Date("2020-01-11")
rtp <- tempdir()

set_database(file.path(tempdir(), "DATABASE"))

navarre <- new_rtoi("Navarre", roi, rtp)

# search, acquire, customize, and process
rsat_search(region = navarre, product = "mod09ga", dates = toi)
rsat_download(navarre)
rsat_mosaic(navarre, overwrite = TRUE)

rsat_derive(navarre, 
            product = "mod09ga", 
            variable = "NDVI")

# plot the results
plot(navarre, "view" , 
      product = "mod09ga", 
      variable = "NDVI", 
      breaks = seq(0, 1, 0.1))
      
plot(navarre,"dates")
```
See the vignettes for more examples:
```
browseVignettes("rsat")
```
## Related similar packages

R has become an outstanding tool for remote sensing image analysis.
There are several tools for the search and acquisition of satellite
images, however, rsat is the first package that standardizes all the
procedures in data acquisition to provide an unique workflow for any
multispectral satellite.

Currently there are several packages dedicated to remote sensing topic,
but they are usually ad-hoc packages for each satellite. Here is a list
of some of the most popular R packages dedicated to satellite imagery:

### Multi satellite packages

-   [RGISTools](https://github.com/spatialstatisticsupna/RGISTools)

-   [getSpatialData](https://github.com/16EAGLE/getSpatialData)

-   [luna](https://github.com/rspatial/luna)

The closest package to `rsat` is RGISTools. `rsat` is the redefinition
of the RGISTools package reprogrammed from scrach in the object-oriented
programming paradigm. Many of the RGISTools code lines have been used to
develop `rsat`, but these have been optimized and redundancies in the
code have been removed in order to facilitate its maintenance. In
addition, `rsat` contains new features and R classes to make it more
user-friendly.

`getSpatialData` is another package very similar to `rsat`. The package
has the same philosophy of having a single package for searching and
downloading satellite images. However, the development of `rsat` goes a
bit further and in addition to search and download, the package helps
you to organize all the downloaded information in a structured database.
`rsat` allows you to use the metadata of the images to see the direct
relation with your region of interest before downloading it. Also all
image processing standardization is not developed in `getSpatialData`.

The last package dedicated to image downloading is `luna`.Searching and
downloading images compared to rsat is a bit more complicated. It is
only able to search and download Modis and Landsat images, and does not
help you in organizing the image products.

### Single satellite packages

-   [rLandsat](https://github.com/atlanhq/rLandsat)

-   [getLandsat](https://github.com/ropensci/getlandsat)

-   [sen2r](https://github.com/ranghetti/sen2r)

`rLandsat` makes it easy to search for Landsat8 product IDs, place an
order on USGS-ESPA and download the data. `rsat` on the other hand is
able to do the image search without knowing all the ids, just using a
polygon of the region of interest, making the search process much
easier.

`getlandsat` provides access to Landsat 8 metadata and images hosted on
AWS S3 at. The package only data for the users, and does not help in
further use, as rsat does.

`sen2r` is an R library which helps to download and preprocess
Sentinel-2 optical images. This is done through a GUI, something that
can be very interesting for users but limits the analysis of the
information prior to downloading, which can be done with `rsat`.

### Raster processing packages

-   [landsat](https://cran.r-project.org/web/packages/landsat/index.html)

-   [satellite](https://github.com/environmentalinformatics-marburg/satellite)

-   [OpenImageR](https://github.com/mlampros/OpenImageR)

-   [RSToolbox](https://github.com/bleutner/RStoolbox)

-   [sits](https://github.com/e-sensing/sits)

`rsat` helps you to search, download and pre-process the images, but
once these procedures are done it allows you to extract all the
processed information into the most used raster classes in R (`raster`,
`stars` or `spatRaster`). The image processing packages can be used for
further analysis in these R classes.

## Contributing
We accept contributions to improve the package. Before contributing, please follow these steps:

-    Contributions should be thoroughly tested with testthat.
-    Code style should attempt to follow the tidyverse style guide.
-    Please attempt to describe what you want to do prior to contributing by submitting an issue.
-    Please follow the typical github fork - pull-request workflow.
-    Make sure you use roxygen and run Check before contributing. More on this front as the package matures.
## Code of conduct
Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By
contributing to this project, you agree to abide by its terms.
## Citation

``` r
citation("rsat")[1]
```

To cite the package:

U. Pérez-Goya, M. Montesino-SanMartin, A F Militino, M D Ugarte (2021). rsat: Dealing with Multiplatform Satellite Images from Landsat, MODIS, and Sentinel. R package version 0.1.16. <https://github.com/ropensci/rsat>.

## Acknowledgements

This work has been financed by projects MTM2017-82553-R (AEI/FEDER, UE) and PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033.

---
title: "2. Download"
output: rmarkdown::html_vignette
bibliography: '`r system.file("REFERENCES.bib", package="rsat")`'
vignette: >
  %\VignetteIndexEntry{rsat2_download}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

# Download

**Downloading** implies acquiring and saving the list of satellite images in a `records` on your machine. This demo builds on the showcase from the search vignette and so, the first section reviews the most important code from the previous vignette. The second section explains how to obtain satellite images with `rsat`. The last section mentions how `rtoi`s are designed to favor collaborative efforts to save downloading time.

------------------------------------------------------------------------

## Review

As a first step of `rsat`'s workflow is specifying the credentials for the the web services:

```{r}
library(rsat)
set_credentials("rsat.package","UpnaSSG.2021")
```

The showcase aims at assessing the effect of the [Snowstorm Filomena](https://en.wikipedia.org/wiki/2020%E2%80%9321_European_windstorm_season#Storm_Filomena) on the Iberian peninsula during January $10^{th}$ and $15^{th}$, $2021$. Hence, the *roi* and *toi* correspond to an `sf` polygon around the peninsula (`ip`) and a vector of dates (`toi`) covering the time-span:

```{r search_review}
ip <- st_sf(st_as_sfc(st_bbox(c(
  xmin = -9.755859,
  xmax =  4.746094,
  ymin = 35.91557,
  ymax = 44.02201 
), crs = 4326)))
toi <- seq(as.Date("2021-01-10"),as.Date("2021-01-15"),1)
```

The folders for the database and dataset can be created programmatically as follows:

```{r}
db.path <- "C:/database"
ds.path <- "C:/datasets"
dir.create(db.path)
dir.create(ds.path)
```

The minimum information to generate a new `rtoi` is the `name`, a polygon of the *roi*, and the paths to database and dataset:

```{r}
filomena <- new_rtoi(name = "filomena",
                     region = ip,
                     db_path = db.path,
                     rtoi_path = ds.path)
```

To limit the amount of data and processing times, the assessment is conducted over MODIS imagery. A total number of $24$ images are found for the region over the $6$-day period:

```{r}
rsat_search(region = filomena, product = c("mod09ga"), dates = toi)
```

------------------------------------------------------------------------

## Image acquisition

Downloading is straightforward with the function `rsat_download()`. The simplest way to use this function is passing the `rtoi` as an input. Depending on the speed of the internet connection, the following instruction may take from few to several minutes to run:

```{r download_rtoi}
rsat_download(filomena)
```

The function saves the satellite images automatically in the database. The path to the database is provided by the `rtoi`:

```{r download_database}
list.files(get_database(filomena), recursive = TRUE)
```

Another way to download images is using a `records`. This variant requires defining a path for saving the resulting files (`out.dir`). The next line is equivalent to the `rtoi` version but using its `records` class object:

```{r download_records}
rsat_download(records(filomena), out.dir = get_database(filomena))
```

This second time, the message reveals that the function reads the database first and checks which images in the `rtoi` are already available in the destination path. If it is available, the function skips its download. This feature becomes handy when teams share a common database.

------------------------------------------------------------------------

## Collaborative `rtoi`s

The `rtoi` leverages the collective use of the package by different working groups within the same or different institutions. Teams working on separate studies or `rtoi`s can refer to the same database increasing its size over time. The database can progressively turn into a local repository. Eventually, new `rtoi`s may find the requested information in the local database, skipping their download, and saving processing time. We encourage `rsat` users to develop common databases when possible on shared machines.

![The `rtoi` architecture and its role in collaborative work](images/rtoi_collaborative.PNG "The rtoi architecture and its role in collaborative work"){width="348"}

The following vignette explains how to customize the satellite images to turn raw data into valuable information for the aim of the analysis.
---
title: "1. Search"
output: rmarkdown::html_vignette
bibliography: '`r system.file("REFERENCES.bib", package="rsat")`'
vignette: >
  %\VignetteIndexEntry{rsat1_search}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

# Search

In `rsat`, **searching** means locating the satellite images of a desired data product for a given region and time of interest (referred to as *roi* and *toi* henceforth). Searching requires **getting access** to satellite imagery and **doing a search** (Section 2), i.e. finding, previewing, and filtering the available information for a *roi* and *toi*.

------------------------------------------------------------------------

## Getting access {#Getting-access}

### Profiles

`rsat` communicates with two public data archives through their *Application Programming Interfaces* (*APIs*). The registration in the following online portals is required to get a full access to satellite images with `rsat`;

-   [USGS](https://ers.cr.usgs.gov/register/) USGS is the sole science agency for the Department of the Interior of United States. Provide access to Modis Images. More information about USGS can be found [Here](https://www.usgs.gov/).
-   [EarthData](https://urs.earthdata.nasa.gov): A repository of NASA's earth observation data-sets. More information about EarthData can be found [here](https://earthdata.nasa.gov/earth-observation-data).
-   [SciHub](https://scihub.copernicus.eu/dhus/#/self-registration), a web service giving access to Copernicus' scientific data hub. Please go [here](https://scihub.copernicus.eu/) to find more details about the data hub.

For convenience, try to use the same *username* and *password* for all of them. To satisfy the criteria in both web services make sure that the [*username* is $4$ characters long and includes a period, number or underscore]{.ul}. The [*password* must be $12$ character long and should include characters with at least one capital letter, and numbers]{.ul}.

### Credentials management

Two functions deal with the *usernames* and *passwords* for the various *APIs* considered in `rsat`; `set_credentials()` and `print_credentials()`*.* The former defines which `username` (first argument) and `password` (second argument) should be used for which `portal` (third argument):

```{r credentials_individual}
library(rsat)
set_credentials("rsat.package","UpnaSSG.2021", "scihub")
set_credentials("rsat.package","UpnaSSG.2021", "earthdata")
```

The latter prints the *usernames* and *passwords* saved during the current `R` session:

```{r credentials_print}
print_credentials()
```

If *usernames* and *passwords* are the same for all the portals, they can be set with a single instruction:

```{r credentials_simultaneous}
set_credentials("rsat.package","UpnaSSG.2021")
```

------------------------------------------------------------------------

## Doing a search {#Doing-a-search}

### Simple search

#### Finding data records

The function `rsat_search()` performs the search of satellite imagery. One way to use this function is providing; (1) a geo-referenced polygon (`sf` class object), (2) the initial and final dates of the relevant time period (`date` vector), and (3) the name(s) of the satellite data product(s) (`character` vector). `rsat` gives access to several satellite data products. Among them, the *surface reflectance* products are frequently used in applied research studies. Their names are:

-   for the Landsat program:

    -   Landsat 4-5 mission: `"LANDSAT_TM_C1"`
    -   Landsat-7 mission: `"LANDSAT_ETM_C1"`
    -   Landsat-8 mission: `"LANDSAT_8_C1"`

-   for the MODIS program:

    -   Terra satellite: `"mod09ga"`
    -   Aqua satellite: `"myd09ga"`

-   for the Sentinel program:

    -   Sentinel-2 mission: `"S2MSI2A"`

    -   Sentinel-3 mission: `"SY_2_SYN___"`

More details about other products and their names can be found [here](https://www.usgs.gov/landsat-missions/product-information), [here](https://modis.gsfc.nasa.gov/data/dataprod/), and [here](https://sentinel.esa.int/web/sentinel/missions) for Landsat, MODIS, and Sentinel respectively. The package can search and download other data than multispectral images (e.g., radar) although this goes beyond the scope of the current status of the package.

The following code looks for satellite images (*surface reflectance*) of a region in norther Spain captured during the first half of January $2021$, captured by MODIS:

```{r records_search}
data("ex.navarre")
toi <- as.Date("2021-01-01") + 0:15
rcd <- rsat_search(region = ex.navarre,
                   product = c("mod09ga"),
                   dates = toi)
```

A search can involve several products simultaneously. In addition to MODIS, the instruction below also considers images from Sentinel-3:

```{r records_search_multiple}
rcd <- rsat_search(region = ex.navarre,
                   product = c("mod09ga", "SY_2_SYN___"),
                   dates = toi)
class(rcd)
```

The output from `rsat_search()` is a `records` object. This object stores the search results metadata from the various *APIs* in standardized manner. A `records` works as a vector and every element of the vector refers to an image. Information about an image is saved in the following *slots*:

-   sat (`character`): name of the satellite (e.g., Landsat, MODIS, Sentinel-2, Sentinel-3, etc.).
-   name (`character`): name of the file.
-   date (`Date`): capturing date of the image.
-   product (`character`): name of the data product.
-   path (`numeric`): horizontal location of the tile (MODIS/Landsat only).
-   row (`numeric`): vertical location of the tile (MODIS/Landsat only).
-   tileid (`character`): tile identification number (Sentinel only).
-   download (`character`): download URL.
-   file_path (`character`): the relative path for local store of the image.
-   preview (`character`): preview URL.
-   api_name (`character`): name of the API internally used by `rsat`.
-   order (`boolean`): if `TRUE`, the image needs to be ordered before the download.
-   extent_crs (`extent_crs`): coordinate reference system of the preview.

You can extract the most relevant *slots* from a `records` using specific functions such as `sat_name()`, `names()`, or `dates()`:

```{r records_slots}
unique(sat_name(rcd))
names(rcd)[1]
unique(dates(rcd))
```

### Filtering search results

As mentioned earlier, a `records` object behaves like a vector. It works with common `R` methods such as `c()`, `[]`, `length()`, `subset()`, or `unique()`. These methods allow to filter and tinker with the search results. Selecting and filtering satellite images is specially powerful when combined with visualization (see below). Discarding useless images at this stage of the process can save memory space and processing time.

For instance, the code below counts the results found from each mission:

```{r records_filter_basic}
mod.rcd <- subset(rcd, "sat", "Modis")
sn3.rcd <- subset(rcd, "sat", "Sentinel-3")
length(mod.rcd)
length(sn3.rcd)
```

The first, second, and third argument in `subset()` correspond to, (1) the `records` object, (2) the name of the slot used for subsetting, and (3) the value of the slot we are interested in. The next lines of code show a more advance filtering example where a new `records` is built from images captured by both missions on the same dates:

```{r records_filter_advanced}
mod.mtch <- mod.rcd[dates(mod.rcd) %in% dates(sn3.rcd)]
sn3.mtch <- sn3.rcd[dates(sn3.rcd) %in% dates(mod.rcd)]
rcd <- c(mod.mtch, sn3.mtch)
```

A `records` can be coerced to a `data.frame` with `as.data.frame()` and transformed back with `as.records()`. The rows and columns of the `data.frame` correspond to the satellite images and their slots respectively:

```{r records_dataframe}
rcd.df <- as.data.frame(rcd)
dim(rcd.df)
names(rcd.df)
# rcd <- as.records(rcd.df)
# class(rcd)
```

### Previewing search results

The `rsat` package provides a `plot()` method for satellite `records`. It displays a preview of the satellite images, i.e a low-resolution versions of the satellite images. Previewing can help to decide whether the actual image, much heavier than the preview, is worth downloading.

Cloudy images are frequently useless and they can be removed from the `records` object using vector-like methods (see previous section). For instance, the code below displays the first $10$ records and, based on visual evaluation of cloud coverage, it removes the $9^{th}$ image from the `records`:

```{r}
prview <- plot(rcd[1:12])
prview
rcd <- rcd[-9]
```

Sometimes, it might be difficult to spot the location of the region of interest on the previews. The *roi* polygon can be passed to the `plot()` function using the `region` argument. Additionally, the `plot()` method in `rsat` is a wrapper function of `tmap` and it accepts other inputs from `tm_raster()` and `tm_polygon()`. These arguments should be preceded by the `tm.raster` and `tm.polygon` identifiers.

The instruction below leverages the preview to its fullest. It shows the *RGB* satellite image preview with the *roi* superposed (), with a red border color (`tm.polygon.region.border.col = "red"`) and no fill (`tm.polygon.region.alpha = 0`). The compass (`compass.rm = TRUE`) and scale bar (`scale.bar.rm = TRUE`) are removed:

```{r}
plot(rcd[1:6],
     region = ex.navarre,
     tm.polygon.region.border.col = "red",
     tm.polygon.region.alpha = 0,
     compass.rm = T,
     scale.bar.rm = T)
```

------------------------------------------------------------------------

## Advanced search

### The `rtoi`

Sometimes, regions are located at the intersection of images or they are too large to fit in a single scene. In these situations, working with separate files might be cumbersome. As a solution, `rsat` provides the `rtoi` class object. An `rtoi` is a project or study-case which is associated with a *region and time of interest* (hence its name). An `rtoi` consists of the following elements:

1.  `name`: a project identifier.
2.  `region`: a georeferenced polygon (`sf` class object).
3.  `db_path`: the path to a database, i.e. a folder for raw satellite images which have generic purposes.
4.  `rtoi_path`: the path to a dataset, i.e. a folder for customized/processed.
5.  `records`: a vector of satellite records relevant for the study.

As a showcase, we'll assess the effects of the [Storm Filomena](https://en.wikipedia.org/wiki/2020%E2%80%9321_European_windstorm_season#Storm_Filomena) over the Iberian peninsula in terms of snow coverage. The storm was an extreme meteorological event (largest since 1971, according to [*AEMET*](http://www.aemet.es/en/portada)). The storm swept the peninsula from January $6^{th}$ and $11^{th}$, $2021$. The code below generates a bounding box around the peninsula (`ip`) and limits the study period (`toi`) to the immediate dates after the storm:

```{r roi_toi}
ip <- st_sf(st_as_sfc(st_bbox(c(
  xmin = -9.755859,
  xmax =  4.746094,
  ymin = 35.91557,
  ymax = 44.02201 
), crs = 4326)))

toi <- seq(as.Date("2021-01-10"),as.Date("2021-01-15"),1)
```

A new `rtoi` requires at least the elements from $1$ to $4$. The following lines generate programmatically the database and dataset folders:

```{r database_dataset}
db.path <- "C:/database"
ds.path <- "C:/datasets"
dir.create(db.path)
dir.create(ds.path)
```

The function `new_rtoi()` builds a new `rtoi` from the elements mentioned above:

```{r rtoi_search}
filomena <- new_rtoi(name = "filomena",
                     region = ip,
                     db_path = db.path,
                     rtoi_path = ds.path)
```

The `new_rtoi()` function writes a set of files representing a copy of the `R` object (filomena.rtoi) and a shapefile of the *roi* (region):

```{r rtoi_file}
rtoi.files <- list.files(file.path(ds.path, "filomena"), full.name = TRUE)
rtoi.files
```

This file is updated whenever the `rtoi` is modified. Thus, if the `R` session accidentally breaks or closes, the last version of the `rtoi` is saved in you hard drive and can be loaded with `read_rtoi()`:

```{r rtoi_read}
filomena <- read_rtoi(file.path(ds.path, "filomena"))
```

Printing the `rtoi` provides a summary about the region and time of interest, a summary of the relevant `records` for the analysis, and the paths to the database and the dataset:

```{r rtoi_print}
print(filomena)
```

### Search with an `rtoi`

You can use an `rtoi` to search satellite images with `rsat_search()`:

```{r rtoi_search_2}
toi <- as.Date("2021-01-10") + 0:5
rsat_search(region = filomena,
            product = c("mod09ga", "SY_2_SYN___"),
            dates = toi)
```

An `rtoi` is an S6 class object, so it behaves like any object in other programming languages. That is, if the `rtoi` is passed as an argument and modifies the `rtoi`, the object in the main environment is also updated. The function `rsat_search()` places the resulting `records` into the `rtoi`:

```{r rtoi_r_update}
print(filomena)
```

Since the `rtoi` has been modified, the `rtoi.file` has also been updated:

```{r rtoi_file_update}
file.info(rtoi.files[1])$ctime # creation time
file.info(rtoi.files[1])$mtime # modification time
```

To extract the `records` from the `rtoi` use `records()`:

```{r rtoi_records}
rcds <- records(filomena)
class(rcds)
```

If you want to apply further filters on the results, as shown in the previous section, you need to extract the `records` from the `rtoi` first and then apply `c()`, `[]`, `subset()`, and `length()` at your convenience.

### Displaying results

There is a `plot()` method too for `rtoi`s but it is more sophisticated. There is a `"preview"` mode that shows the low-resolution *RGB* version of the satellite images. The plot binds together the satellite images available for a given date and `roi`, and the function allows other arguments, such as `product` or `dates`, to better select the information being shown:

```{r rtoi_preview}
plot(filomena,
     "preview",
     product = "mod09ga",
     dates = "2021-01-11")
```

`rtoi`s are designed save information from several missions. The `mode = "dates"` prints a calendar indicating the availability of satellite images from one or multiple missions during the time of interest:

```{r rtoi_calendar}
plot(filomena, "dates")
```

The following vignette explains how to download the satellite images and the role of the *database* when managing multiple `rtoi`s. To reduce the processing times and memory space, the showcase is restricted to MODIS images from here onward:

```{r rtoi_lightening}
rcd <- records(filomena)
records(filomena) <- subset(rcd, "product", "mod09ga")
```
---
title: "4. Process"
output: rmarkdown::html_vignette
bibliography: '`r system.file("REFERENCES.bib", package="rsat")`'
vignette: >
  %\VignetteIndexEntry{rsat4_process}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

# Processing

Processing encapsulates the operations applied to satellite images to solve particular issues derived from errors, anomalies or discrepancies from sensors. For instance, cloud removal or sensor failures can lead to data gaps in time-series of satellite images. Additionally, noise from aerosols, dust, and sensor measurement errors can reduce the quality of the remotely sensed data. The ***Interpolation of Mean Anomalies (IMA***) is a gap-filling and smoothing approach that mitigates these issues [@militino2019interpolation].

### Filling missing data {#filling}

#### Theoretical background

`rsat` implements a generic version of the algorithm with the function `smoothing_image()`.

![Neighborhood definition in the Interpolation of Mean Anomalies, assuming `nDays = 1` and `nYears = 1`](images/neighborhood_ima.png "Neighborhood defintion in IMA"){width="456"}

For the theoretical explanation, let's assume that $15$ images are available in the dataset (squares in the image above). The imagery corresponds to $5$ consecutive days over $3$ years. Consider that we are willing to fill/smooth the target image (red square). *IMA* fills the gaps borrowing information from an *adaptable* temporal neighborhood (yellow squares). Two parameters determine the size of the neighborhood; the number of days before and after the target image (`nDays`) and the number of previous and subsequent years (`nYears`). Both parameters should be adjusted based on the temporal resolution of the of the time-series of images. We recommend that the neighborhood extends over days rather than years, when there is little resemblance between seasons. Also, cloudy series may require larger neighborhoods.

![Summarized algorithm of the Interpolation of Mean Anomalies](images/ima_technique.png "IMA summarized algorithm")

*IMA* gives the following steps; (1) creates a representative image of the neighborhood ignoring missing values e.g., doing the mean, median, etc. for each pixel's time-series (`fun`), (2) the target and representative images are subtracted giving an image of *anomalies*, (3) the *anomalies* falling outside the quantile limits (`aFilter`) are considered outliers and therefore removed, (4) it aggregates the anomaly image into a coarser resolution (`fact`) to reveal potential spatial dependencies, (5) the procedure fits a spatial model (thin plate splines or *TPS*) to the anomalies which is then used to interpolate the values at the original resolution, and (6) the output is the sum of the interpolated anomalies and the average image. The process is encapsulated in `smoothing_image()` and the section below shows its usage.

------------------------------------------------------------------------

#### Hands-on demonstration

##### Data

The showcase uses a time-series of $6$ *NDVI* images. *NDVI* stands for *Normalized Difference Vegetation Index* and it is an indicator of vegetation's vigor. Values close to $1$ indicate the presence of green and dense vegetation. Images correspond to a region in norther Spain between August $2^{nd}$ and $4^{th}$, in $2017$ and $2018$. The images are part of the package and can be loaded into the environment as follows:

```{r advanced_ima_data}
library(rsat)
library(terra)

data("ex.ndvi.navarre")
ex.ndvi.navarre <- rast(ex.ndvi.navarre)
```

Let's display the series of images to spot the gaps of data:

```{r advanced_ima_show}
library(tmap)
tm_shape(ex.ndvi.navarre) + tm_raster(title = "NDVI", style = "cont")
```

The maps show a hole on August $3^{rd}$, $2017$ in the upper part of the image. The aim of *IMA* is to predict the values of the missing parts to provide a complete and continuous dataset.

##### Interpolation of Mean Anomalies

The instruction below applies the *IMA* to the entire series. It uses a temporal neighborhood of $2$ days and $1$ year and the representative image is the `mean`. Outliers are the anomalies outside $1-99$ quantile range. Anomalies are aggregated at a factor of $10$ (every $10$ pixels are aggregated into $1$) to fit an approximate spatial model. We recommend to tune the parameters accordingly to obtain the best performance possible in each situation:

```{r advanced_ima}
library(rsat)
ndvi.fill <- rsat_smoothing_images(method = "IMA",
                                   ex.ndvi.navarre,
                                   nDays = 2,
                                   nYears = 1,
                                   fun = mean,
                                   aFilter = c(0.01,0.99),
                                   fact = 10,
                                   only.na = TRUE)
```

The algorithm predicts the value for the missing and the existing pixels. The last argument, `only.na = TRUE`, indicates that the results should preserve the original values wherever available. *IMA* is able to run when neighborhoods are incomplete. In those situations, the algorithm simply uses the imagery available within the temporal neighborhood. The function prints a message specifying the actual number of images being used in each case.

Let's make a comparison between *before* and *after* the application of *IMA.*

```{r advanced_ima_result}
before <- ex.ndvi.navarre[[1:3]]
after <- ndvi.fill[[1:3]]
tm_shape(c(before,after)) + tm_raster(title = "NDVI", style = "cont")
```

The intention is to keep growing the number of algorithms devoted to processing in future versions of the package.
---
title: "3. Customize"
output: rmarkdown::html_vignette
bibliography: '`r system.file("REFERENCES.bib", package="rsat")`'
vignette: >
  %\VignetteIndexEntry{rsat3_customize}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

# Customize

**Customizing** means transforming raw images into useful information for particular needs. Customizing includes; (1) mosaicking, i.e. joining scenes of a region and date in a single file , (2) calculating and index to highlight the presence of process/material in the satellite image, and (3) mask the cloudy pixels to avoid misinterpreting the surface reflectance values. This demo builds on the showcase from the search vignette and so, the first section reviews the most important code from the previous vignette.

------------------------------------------------------------------------

## Review

As a first step of `rsat`'s workflow is specifying the credentials for the the web services:

```{r}
library(rsat)
set_credentials("rsat.package","UpnaSSG.2021")
```

The showcase aims at assessing the effect of the [Snowstorm Filomena](https://en.wikipedia.org/wiki/2020%E2%80%9321_European_windstorm_season#Storm_Filomena) on the Iberian peninsula during January $10^{th}$ and $15^{th}$, $2021$. Hence, the *roi* and *toi* correspond to an `sf` polygon around the peninsula (`ip`) and a vector of dates (`toi`) covering the time-span:

```{r search_review}
ip <- st_sf(st_as_sfc(st_bbox(c(
  xmin = -9.755859,
  xmax =  4.746094,
  ymin = 35.91557,
  ymax = 44.02201 
), crs = 4326)))
toi <- seq(as.Date("2021-01-10"),as.Date("2021-01-15"),1)
```

The folders for the database and dataset can be created programmatically as follows:

```{r}
db.path <- "C:/database"
ds.path <- "C:/datasets"
dir.create(db.path)
dir.create(ds.path)
```

The minimum information to generate a new `rtoi` is a `name` for the object, a polygon of the region of interest (*roi*), and the paths to the database and dataset:

```{r}
filomena <- new_rtoi(name = "filomena",
                     region = ip,
                     db_path = db.path,
                     rtoi_path = ds.path)
```

To limit the amount of data and processing times, the assessment is conducted over MODIS imagery. A total number of $24$ images are found for the region over the $6$-day period:

```{r}
rsat_search(region = filomena, product = c("mod09ga"), dates = toi)
```

The way to download the search results is as follows:

```{r}
rsat_download(filomena)
```

------------------------------------------------------------------------

## Mosaic

***Mosaicking*** involves binding together several images of a region from the same date. The function `mosaic()` finds automatically the relevant images in the database and joins them together in a single file. Additionally, by default, the function crops around the *roi* of the `rtoi` to remove unnecessary information and save space on your hard disk drive:

```{r, eval=FALSE}
rsat_mosaic(filomena)
```

The cropping option can be disabled with the argument `warp = NULL`. Here, cropping is appropriate since images extend far beyond our region of interest.

The results are saved under `rtoi_path` inside Modis/mod09ga/mosaic (i.e. *constellation_name/data_product_name/mosaic*). This is the first time in the workflow that the `rtoi_path` is being used. The reason is that mosaicking is the first transformation applied to the raw images to better fit the particular needs of the analysis. The outcomes from the mosaic are compressed (*zip)* to minimize their size:

```{r, eval=FALSE}
list.files(file.path(ds.path, "filomena", "Modis/mod09ga/mosaic"), full.name = TRUE)
```

At this point of the workflow, *RGB* representations of the satellite images can be displayed on the fly, without loading the entire image in `R`. By default, the `plot()` method for an `rtoi` displays the mosaicked images when the object is not followed by the words `"preview"` or `"date"` (see other types of plots for an `rtoi` in the search vignette):

```{r}
plot(filomena, as.Date("2021-01-11"))
```

The map shows a sample of pixels from the original image. This is a strategy to save space in *RAM* memory. The downside is the reduction in the level of detail of the image. By default, the number of pixels on the horizontal and vertical axis is $250$. The arguments `xsize` and `ysize` change the size of the sample to vary the crispness of the image. The code below doubles the number of pixels on each axis of the image;

```{r}
plot(filomena, as.Date("2021-01-11"),xsize = 500, ysize = 500)
```

Clouds and snow are visually similar. False color images are frequently used in the field in order to ease the detection of features. False color images switch the natural color in the *RGB* image to other bands. These kind of representations can be built using the `band_name` argument followed by the name of the bands replacing the *red*, *green*, and *blue* colors. For instance, the rendering of *swir1*, *nir*, and *blue* highlights the snow in blue color:

```{r}
plot(filomena,
     as.Date("2021-01-11"),
     xsize = 500,
     ysize = 500,
     band_name = c("swir1", "nir", "blue"))
```

## Index calculation {#index-calculation}

### Definition

A ***remote sensing index*** is an indicator that reveals the presence of a material in a satellite image. Indexes are the result of simple math applied to the bands of an image. The computation involves the bands with a distinctively high or low reflectance for the feature at hand. Over the years, researchers have developed a wide variety of indexes for different materials or processes which can be consulted [here](https://www.indexdatabase.de/db/i.php).

For instance, the *Normalized Difference Snow Index* (NDSI) (see e.g., [@salomonson2004estimating]) highlights the snow using the *green* and *shortwave-infrared* bands (around $1.5 \mu m$). The subtraction of this two bands gives a large number for those pixels depicting snow. The denominator ensures that values oscillate between $-1$ and $1$.

$$ NDSI = \frac{Green - SWIR1}{Green + SWIR1}$$

### Calculation

In `R` we can create a function replicating the calculation of the *NDSI* index*:*

```{r basic_ndsi, eval = FALSE}
NDSI = function(green, swir1){
  ndsi <- (green - swir1)/(green + swir1)
  return(ndsi)
}
```

`rsat` demands that these formulas use the band names, e.g. *red, green, blue, etc.* rather than band number. Band names and numbers differ among mission/satellites. For instance, the *green* corresponds to the band number $4$ in MODIS and Landsat-7, number $3$ in Landsat-8 and Sentinel-2, and number $6$ in Sentinel-3 (see [here](https://drive.google.com/file/d/1cSw4LaTLPlGBHmG8v7uwH54f-m9jZz1N/view?usp=sharing)). Using their names enables the use of a unique custom function across satellites/missions. `rsat` functions are responsible for linking the name to the corresponding band number according to the mission. Some widespread variables are built-in the package and the list of variables can be printed using;

```{r basic_variables}
show_variables()
```

To use the `NDSI()` function over the series of satellite images of the Iberian peninsula, use the function `derive()` as follows;

```{r basic_derive, eval = FALSE}
rsat_derive(filomena, product = "mod09ga", variable = "ndsi", fun = NDSI)
```

Again, you can plot the results without loading the scenes in `R`:

```{r}
plot(filomena,
     as.Date("2021-01-11"),
     variable = "ndsi",
     xsize = 500,
     ysize = 500,
     zlim = c(-1,1))
```

The *NDSI* index improves the separability between clouds and snow. However, there might be some difficulties distinguishing between them in certain parts of the image. As a solution, the next step removes cloud-covered pixels.

## Cloud removal {#cloud-removal}

Some data providers apply algorithms over their data-sets to detect the presence of clouds (Level 1/2 products). The analysis is part of the quality assessment done during pre-processing and the results are included in the ***Quality Assurance*** (*QA*) band of the image. In addition to cloud coverage, the band provides information about over-saturated or filled pixels. The information is packed in this band using the bit format.

The function `cloud_mask()` interprets the *QA* band to obtain images showing the presence/absence of clouds. Its application is straightforward;

```{r basic_cloud, eval=FALSE}
rsat_cloudMask(filomena)
```

To apply the cloud mask, we need to import the *NDSI* images into `R` using the `get_raster()` function. The values of the index must be truncated between $-1$ and $1$ to avoid values outside the feasible range (sun reflections on mirror-like surfaces, such as water, can lead to misleading results):

```{r basic_ndsi_import, eval = FALSE}
ndsi.img <- rsat_get_raster(filomena, "mod09ga", "ndsi")
ndsi.img <- clamp(ndsi.img, -1, 1)
```

For every image in the `rtoi`, the `cloud_mask()` function generates a new image, called ***mask***, in which $1$s and $NA$s indicate clear and covered pixels. The function identifies the mission/program and applies the appropriate interpretation of bits to create the cloud mask. To import the result run;

```{r basic_mask, eval = FALSE}
clds.msk <- rsat_get_raster(filomena, "mod09ga", "CloudMask")
```

In MODIS, cloud-masks have a different resolution than the multispectral image. To adjust the resolution, *resample* the cloud mask to match the resolution of the *NDSI* images (`resample()`) using the nearest neighbor method (`"ngb"`):

```{r basic_mask_resample}
clds.msk <- resample(clds.msk, ndsi.img, method = "ngb")
```

To apply the cloud mask, we just multiply both series of pixels. Dot multiplications are performed pixel-wise. *NDSI* values multiplied by $1$ remain unaltered but those multiplied by $NA$ become missing:

```{r basic_mask_apply}
ndsi.filt <- ndsi.img * clds.msk
names(ndsi.filt) <- names(clds.msk) # keep the names
```

As an attempt to obtain a ***composite image***, we extract maximum value of the *NDSI* for each pixel in the time series. Maximum value compositions are frequent in this field [@holben1986characteristics]. Compositing is as a way to summarize the information in a time-lapse and ignore the presence of clouds:

```{r basic_composite}
snow.spain <- calc(ndsi.filt, max, na.rm = TRUE)
```

Represent the results:

```{r basic_ndsi_map}
library(tmap)
tm_shape(snow.spain) + tm_raster(style = "cont")
```

The results are not completely satisfactory. The image shows unfilled gaps and this is because for those pixels there is no valid information along the time-series. The following vignette explains how to process images to fill gaps and smooth outliers.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\docType{class}
\name{rtoi-class}
\alias{rtoi-class}
\title{Region and Time Of Interest (\code{rtoi})}
\description{
It is a proxy object to store metadata about satellite imagery
covering a spatial region over a time period. Images can come from
multiple missions/programs and its purpose is to help managing
heterogeneous datasets.
}
\details{
An \code{rtoi} object manages two main folders called database and rtoi.
The database is meant to work as a local, generic, and organized archive
of raw satellite data retrieved with the \code{download()} function.
The rtoi folder contains processed information for
a particular region and time of interest. When \code{mosaic()}
is called, the function crops and mosaics the relevant raw images from
the database and saves the results in the rtoi folder. This folder also
contains a \code{region.rtoi} file which saves metadata about the
region/time of interest and satellite imagery available.
}
\section{Fields}{

\describe{
\item{\code{name}}{a character with the name of the region of interest.}

\item{\code{rtoi_path}}{a character with the path to the rtoi folder.}

\item{\code{region}}{an sf with the region of interest.}

\item{\code{records}}{the satellite records available for
your region and time of interest.}

\item{\code{db_path}}{a character with the path to the database.}
}}


\examples{
data(ex.navarre)
## Create an rtoi with database
# path where the region is stored
rtoi.path <- tempdir()

# path where downloads are stored
db.path <- file.path(tempdir(), "DATABASE")
navarre <- new_rtoi(
  name = "Navarre_rtoi",
  region = ex.navarre,
  rtoi_path = rtoi.path,
  db_path = db.path
)

print(navarre)

## Create an rtoi without database
navarre2 <- new_rtoi(
  name = "Navarre_rtoi2",
  region = ex.navarre,
  rtoi_path = rtoi.path
)

print(navarre2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_list_data.R
\name{rsat_list_data}
\alias{rsat_list_data}
\alias{rsat_list_data,rtoi-method}
\alias{rsat_list_data,rtoi}
\title{List the information available for an \code{rtoi}}
\usage{
rsat_list_data(x, ...)

\S4method{rsat_list_data}{rtoi}(x, ...)
}
\arguments{
\item{x}{an \code{rtoi} object.}

\item{...}{additional arguments.}
}
\value{
a \code{data.frame} of the available information.
}
\description{
Displays the existing products, bands,
and processing levels for a given \code{rtoi}
}
\examples{
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

print(navarre)

# print empty rtoi
rsat_list_data(navarre)

file.copy(from=system.file("ex/Pamplona",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
pamplona <- read_rtoi(file.path(tempdir(),"Pamplona"))

print(pamplona)

rtoi.data <- rsat_list_data(pamplona)
# print mosaicked bands
print(rtoi.data)

# print mosaicked bands + derived NDVI
file.copy(from=system.file("ex/PamplonaDerived",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
pamplona.derived <- read_rtoi(file.path(tempdir(),"PamplonaDerived"))
rsat_list_data(pamplona.derived)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api.R, R/extent_crs.R, R/records.R, R/rtoi.R,
%   R/variables.R
\name{print,api-method}
\alias{print,api-method}
\alias{print,api}
\alias{print,extent_crs-method}
\alias{print,extent_crs}
\alias{print,records-method}
\alias{print,records}
\alias{print,rtoi-method}
\alias{print,variables-method}
\title{Prints the values}
\usage{
\S4method{print}{api}(x)

\S4method{print}{extent_crs}(x)

\S4method{print}{records}(x)

\S4method{print}{rtoi}(x)

\S4method{print}{variables}(x, ...)
}
\arguments{
\item{x}{an object to be printed..}

\item{...}{additional arguments.}
}
\value{
prints rtoi metadata
}
\description{
prints an object and returns it invisibly (via invisible(x)).
}
\examples{
\dontrun{
library(rsat)

# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

print(navarre)

# get records
rcrds <- records(navarre)

print(rcrds)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_mosaic.R
\name{rsat_mosaic}
\alias{rsat_mosaic}
\alias{rsat_mosaic,rtoi-method}
\alias{rsat_mosaic,sf,character}
\alias{rsat_mosaic,records-method}
\alias{rsat_mosaic,records}
\title{Mosaic the tiles intersecting the region of interest}
\usage{
rsat_mosaic(x, ...)

\S4method{rsat_mosaic}{rtoi}(x, ...)

\S4method{rsat_mosaic}{records}(
  x,
  out_path,
  db_path,
  bfilter,
  warp = "extent",
  region,
  overwrite = FALSE,
  ...
)
}
\arguments{
\item{x}{an \code{rtoi} object.}

\item{...}{additional arguments.}

\item{out_path}{path to save the mosaicked images. By default, the path
is defined by \code{x}.}

\item{db_path}{path to the database. By default, the path
is defined by \code{x}.}

\item{bfilter}{a vector of bands to. If not supplied, all are used.}

\item{warp}{character. If equal to "extent", it also crops the images
around the \code{rtoi}. Use "" otherwise.}

\item{region}{an sf object. Region for cropping the images around.
By default, the path is defined by \code{x}.}

\item{overwrite}{logical argument. If \code{TRUE}, overwrites the existing
images with the same name.}
}
\value{
nothing. Mosaics the downloaded images and stored them on the hard disk
}
\description{
Satellite measurements are divided into indivisible units called tiles. The
mosaic function binds and crops the tiles to generate a single image
of the region of interest for each date.
}
\examples{
\dontrun{
library(rsat)

# load navarre sf from the package
data(ex.navarre)

# set the credentials
set_credentials("username", "password")

# path where the region is stored
rtoi.path <- tempdir()
# path where downloads are stored
db.path <- file.path(tempdir(), "DATABASE")
navarre <- new_rtoi(
  "Navarre",
  ex.navarre,
  rtoi.path,
  db.path
) #'
# Landsat-5
rsat_search(
  region = navarre,
  product = "LANDSAT_TM_C1",
  dates = as.Date("1988-08-01") + seq(1, 35)
)
rsat_download(navarre)

rsat_mosaic(navarre, overwrite = T)

rsat_list_data(navarre)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ex.ndvi.navarre}
\alias{ex.ndvi.navarre}
\title{A time series of NDVI in Navarre (Spain)}
\format{
The \code{RasterBrick} contains 6 images, from the 2nd to the 4th of
August in 2017 and 2018. The \code{RasterBrick} coordinates are in the
Sinusoidal projection:

\describe{
  \item{name}{layer names contain the date of the
  image in the format "\code{YYYYJJJ}"}.
  \item{size}{each layer contains 113 rows and 105 columns}.
}
}
\description{
Geographically projected \code{RasterBrick} object of the normalized
difference vegetation index (NDVI) in Navarre.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_preview.R
\name{rsat_preview}
\alias{rsat_preview}
\alias{rsat_preview,rtoi,Date-method}
\alias{rsat_preview,rtoi,date}
\alias{rsat_preview,rtoi,missing-method}
\alias{rsat_preview,rtoi,missing}
\alias{rsat_preview,records,Date-method}
\alias{rsat_preview,records,date}
\alias{rsat_preview,records,numeric-method}
\alias{rsat_preview,rtoi,numeric}
\title{Preview a \code{records} or an \code{rtoi} object}
\usage{
rsat_preview(x, n, ...)

\S4method{rsat_preview}{rtoi,Date}(x, n, lpos = c(3, 2, 1), add.layer = FALSE, verbose = FALSE, ...)

\S4method{rsat_preview}{rtoi,missing}(x, n, lpos = c(3, 2, 1), add.layer = FALSE, verbose = FALSE, ...)

\S4method{rsat_preview}{records,Date}(
  x,
  n,
  lpos = c(3, 2, 1),
  tmp_dir = file.path(tempdir()),
  add.layer = FALSE,
  verbose = FALSE,
  get.map = TRUE,
  ...
)

\S4method{rsat_preview}{records,numeric}(
  x,
  n,
  lpos = c(3, 2, 1),
  tmp_dir = file.path(tempdir()),
  add.layer = FALSE,
  verbose = FALSE,
  get.map = TRUE,
  ...
)
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}

\item{n}{the date expressed as the temporal index in the time series.}

\item{...}{additional arguments}

\item{lpos}{vector argument. Defines the position of the red-green-blue
layers to enable a false color visualization.}

\item{add.layer}{logical argument. If \code{TRUE}, the function plots the
image on an existing map.}

\item{verbose}{logical argument. If \code{TRUE}, the function prints the
running steps and warnings.}

\item{tmp_dir}{character argument. The directory where preview images
are located.}

\item{get.map}{logical argument. If \code{TRUE}, the function
return the leaflet map.}
}
\value{
nothing. Previews the region in the viewer.
}
\description{
Preview a \code{records} or an \code{rtoi} object
}
\examples{
\dontrun{
library(rsat)

# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

set_credentials("username", "password")
set_database(file.path(tempdir(), "DATABASE"))

# by default the first date in rtoi is previewed
rsat_preview(navarre)


preview.dates <- dates(navarre)
# use add.layer to preview images of several days
rsat_preview(navarre,preview.dates[2],add.layer = TRUE)

# you can also preview records
rcrds <- records(navarre)
rsat_preview(rcrds, n = 1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_smoothing_images.R
\name{rsat_smoothing_images}
\alias{rsat_smoothing_images}
\alias{rsat_smoothing_images,rtoi,character-method}
\alias{rsat_smoothing_images,rtoi,character}
\alias{rsat_smoothing_images,SpatRaster,character-method}
\title{Fill data gaps and smooth outliers in a time series of satellite images}
\usage{
rsat_smoothing_images(x, method, ...)

\S4method{rsat_smoothing_images}{rtoi,character}(
  x,
  method,
  product = "ALL",
  satellite = "ALL",
  stage = "ALL",
  variable = "ALL",
  test.mode = FALSE,
  ...
)

\S4method{rsat_smoothing_images}{SpatRaster,character}(x, method, ...)
}
\arguments{
\item{x}{\code{rtoi} or \code{RastespatRaster} containing
a time series of satellite images.}

\item{method}{character argument. Defines the method used
for processing the images, e.a. "IMA".}

\item{...}{arguments for nested functions:
\itemize{
  \item \code{Img2Fill}  a \code{vector} defining the
  images to be filled/smoothed.
  \item \code{r.dates} a \code{vector} of dates for the layers in \code{x}.
Mandatory when layer names of \code{x} do not contain their
 capturing dates
"\code{YYYYJJJ}" format.
  \item \code{nDays} a \code{numeric} argument with the number
  of previous and subsequent days of the temporal neighborhood.
  \item \code{nYears} a \code{numeric} argument with the number
  of previous and subsequent years of the temporal neighborhood.
  \item \code{aFilter} a \code{vector} of lower and upper
  quantiles defining the outliers in the anomalies. Ex. c(0.05,0.95).
  \item \code{fact} a \code{numeric} argument specifying the aggregation
  factor of the anomalies.
  \item \code{fun} a \code{function} used to aggregate the
  image of anomalies. Both \code{mean} (default) or \code{median} are
  accepted.
  \item \code{snow.mode} logical argument. If \code{TRUE},
  the process is parallelized using the functionalities from the `
  \code{raster}' package.
  \item \code{predictSE} calculate the standard error instead
  the prediction.
  \item \code{factSE} the \code{fact} used in the standard error
  prediction.
  \item \code{out.name} the name of the folder containing the
  smoothed/filled images when saved in the Hard Disk Device (HDD).
  \item \code{only.na} logical argument. If \code{TRUE}
  only fills the \code{NA} values.
\code{FALSE}  by default.
}}

\item{product}{character argument. The name of the product to
be processed. Check the name of the parameter with \code{\link{rsat_list_data}}
function. Check the name of the parameter with
\code{\link{rsat_list_data}} function. By default, "ALL".}

\item{satellite}{character argument. The name of the satellite to
be processed. Check the name of the parameter with
\code{\link{rsat_list_data}} function. By default, "ALL".}

\item{stage}{character argument. The name of the processed stage
of the data. Check the name of the parameter with
\code{\link{rsat_list_data}} function. By default, "ALL".}

\item{variable}{character argument.The name of the variable to
be processed. Check the name of the parameter with
\code{\link{rsat_list_data}} function. By default, "ALL".}

\item{test.mode}{logical argument. If \code{TRUE}, the function runs some
lines to test \code{rsat_smoothing_images} with rtoi object.}
}
\value{
a \code{RastespatRaster} with the filled/smoothed images.
}
\description{
\code{apply_ima} is the implementation of a spatio-temporal method
called Interpolation of Mean Anomalies(IMA) for gap filling and smoothing
satellite data \insertCite{militino2019interpolation}{rsat}.
\code{smoothing_images} is the implementation of a spatio temporal method
called image mean anomaly (IMA) for gap filling and smoothing satellite
data \insertCite{militino2019interpolation}{rsat}.
}
\details{
This filling/smoothing method was developed by
\insertCite{militino2019interpolation;textual}{rsat}. IMA fills the gaps
borrowing information from an adaptable temporal neighborhood. Two
parameters determine the size of the neighborhood; the number of days
 before and after the target image (\code{nDays}) and the number of previous
and subsequent years (\code{nYears}). Both parameters should be adjusted
based on the temporal resolution of the of the time-series of images. We
recommend that the neighborhood extends over days rather than years, when
there is little resemblance between seasons. Also, cloudy series may require
larger neighborhoods.

IMA gives the following steps; (1) creates a representative image from the
temporal neighborhood of the target image (image to be filled/smoothed) e.g.,
doing the mean, median, etc. for each pixel's time-series (\code{fun}), (2)
the target and representative images are subtracted giving an image of
anomalies, (3) the anomalies falling outside the quantile limits
(\code{aFilter}) are considered outliers and therefore removed, (4) it
aggregates the anomaly image into a coarser resolution (\code{fact}) to
reveal potential spatial dependencies, (5) the procedure fits a spatial
model (thin plate splines or TPS) to the anomalies which is then used to
interpolate the values at the original resolution, and (6) the output
is the sum of the interpolated anomalies and the average image.
}
\examples{
## Smooth data in rtoi
library(rsat)
require(terra)

# create a copy of pamplona in temp file
file.copy(from=system.file("ex/PamplonaDerived",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
pamplona <- read_rtoi(file.path(tempdir(),"PamplonaDerived"))
rsat_smoothing_images(pamplona,
                      method = "IMA",
                      variable="NDVI"
)
rsat_list_data(pamplona)
# get smoothed
smoothed <- rsat_get_SpatRaster(pamplona,p="mod09ga",v="NDVI",s="IMA")
plot(smoothed)

# get original
original <- rsat_get_SpatRaster(pamplona,p="mod09ga",v="NDVI",s="variables")
plot(original)
plot(smoothed[[1]]-original[[1]])

## smooth user defined SpatRaster dataset
require(terra)
data(ex.ndvi.navarre)

# load an example of NDVI time series in Navarre
ex.ndvi.navarre <- rast(ex.ndvi.navarre)
# the raster stack with the date in julian format as name
plot(ex.ndvi.navarre)

# smoothin and fill all the time series
tiles.mod.ndvi.filled <- rsat_smoothing_images(ex.ndvi.navarre,
  method = "IMA"
)
# show the filled images
plot(tiles.mod.ndvi.filled)

# plot comparison of the cloud and the filled images
tiles.mod.ndvi.comp <- c(
  ex.ndvi.navarre[[1]], tiles.mod.ndvi.filled[[1]],
  ex.ndvi.navarre[[2]], tiles.mod.ndvi.filled[[2]]
)
plot(tiles.mod.ndvi.comp)
}
\references{
\insertRef{militino2019interpolation}{rsat}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extent_crs.R, R/records.R
\name{as.data.frame,extent_crs-method}
\alias{as.data.frame,extent_crs-method}
\alias{as.data.frame,extent_crs}
\alias{as.data.frame,records-method}
\title{Coerce to a Data Frame}
\usage{
\S4method{as.data.frame}{extent_crs}(x)

\S4method{as.data.frame}{records}(x)
}
\arguments{
\item{x}{Any R object.}
}
\value{
returns a data frame, normally with all row names
}
\description{
Functions to check if an object is a data frame, or coerce it if possible.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)
# coerce the records to rtoi
df <- as.data.frame(rcds)
# print the dataframe
print(df)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{as.records}
\alias{as.records}
\alias{as.records,data.frame-method}
\alias{as.records,data.frame}
\title{Create records object from data frame}
\usage{
as.records(x)

\S4method{as.records}{data.frame}(x)
}
\arguments{
\item{x}{a \code{data.frame} with columns representing the slots of
records.}
}
\value{
returns a records objects with the columns values in \code{x}
}
\description{
Create records object from data frame
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)
# coerce the records to dataframr
df <- as.data.frame(rcds)
# print the dataframe
print(df)

# coerce the dataframe to records
rcds2 <- as.records(df)
# check the conversion
identical(rcds,rcds2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connections.R
\name{print_credentials}
\alias{print_credentials}
\alias{print_credentials,ANY-method}
\title{Prints the credentials for the web services}
\usage{
print_credentials(...)

\S4method{print_credentials}{ANY}()
}
\arguments{
\item{...}{additional arguments.}
}
\value{
print the credentials asigned in the package environment variable
}
\description{
Prints the credentials for the web services
}
\examples{
print_credentials()
set_credentials("example", "example", "earthdata")
print_credentials()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_derive.R
\name{rsat_derive}
\alias{rsat_derive}
\alias{rsat_derive,rtoi,character-method}
\alias{rsat_derive,rtoi,character}
\title{Computes a remote sensing index from an \code{rtoi}}
\usage{
rsat_derive(x, variable, ...)

\S4method{rsat_derive}{rtoi,character}(
  x,
  variable,
  product,
  dates,
  fun,
  overwrite = FALSE,
  verbose = FALSE,
  suppressWarnings = TRUE,
  ...
)
}
\arguments{
\item{x}{an \code{rtoi} as the source of images.}

\item{variable}{the name of the variable.}

\item{...}{additional argument for variable deriving}

\item{product}{the name of the product from which the index is computed.}

\item{dates}{a vector of dates being considered (optional).}

\item{fun}{a \code{function} that computes the remote sensing index.}

\item{overwrite}{logical argument. If \code{TRUE}, overwrites the existing
images with the same name.}

\item{verbose}{logical argument. If \code{TRUE}, the function prints the
running steps and warnings.}

\item{suppressWarnings}{evaluates its expression in a context that ignores all warnings.}
}
\value{
nothing. The derived variables will be save in the hard drive.
Use get_stars to get the variables.
}
\description{
Combines the bands from multispectral satellite products through simple
math to highlight a process or material in the image.
}
\details{
The package contemplates some pre-defined indexes, which can be displayed
using the \code{show_variables()} function. To compute one of those, write
its name in the \code{variable} argument. Custom indexes can be
supplied through the \code{fun} argument. The function should use the
name of the bands as inputs (red, green, blue, nir, swir1, or swir2) and
return a single element. For instance, the Normalized Difference Snow
Index would be;

NDSI = function(green, swir1){
ndsi <- (green - swir1)/(green + swir1)
return(ndsi)
}
}
\examples{
library(rsat)

# create a copy of pamplona in temp file
file.copy(from=system.file("ex/Pamplona",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
pamplona <- read_rtoi(file.path(tempdir(),"Pamplona"))

rsat_list_data(pamplona)
# show prefedined varibles
show_variables()
rsat_derive(pamplona, "NDVI", product = "mod09ga")
# now NDVI is processed
rsat_list_data(pamplona)

# ad-hoc variable
NDSI = function(green, swir1){
ndsi <- (green - swir1)/(green + swir1)
return(ndsi)
}
rsat_derive(pamplona, "NDSI", product = "mod09ga",fun=NDSI)
# now NDVI is processed
rsat_list_data(pamplona)
plot(pamplona, product="mod09ga",variable="NDSI")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\name{rsat_get_raster}
\alias{rsat_get_raster}
\alias{rsat_get_raster,rtoi-method}
\alias{rsat_get_raster,rtoi}
\alias{rsat_get_SpatRaster}
\alias{rsat_get_SpatRaster,rtoi-method}
\alias{rsat_get_SpatRaster,rtoi}
\alias{rsat_get_stars}
\alias{rsat_get_stars,rtoi-method}
\alias{rsat_get_stars,rtoi}
\title{Loads into R a time series of images regarding an rtoi, satellite product,
and remote sensing index.}
\usage{
rsat_get_raster(x, p, v, s, ...)

\S4method{rsat_get_raster}{rtoi}(x, p, v, s, ...)

rsat_get_SpatRaster(x, p, v, s, ...)

\S4method{rsat_get_SpatRaster}{rtoi}(x, p, v, s, ...)

rsat_get_stars(x, p, v, s, ...)

\S4method{rsat_get_stars}{rtoi}(x, p, v, s, ...)
}
\arguments{
\item{x}{an rtoi.}

\item{p}{a character with the name of the satellite data product.}

\item{v}{a character with the name of the index.}

\item{s}{a character with the name of the stage wanted.}

\item{...}{additional arguments.}
}
\value{
a raster stack.
}
\description{
Loads into R a time series of images regarding an rtoi, satellite product,
and remote sensing index.
}
\examples{
\dontrun{
library(rsat)
# load example rtoi
file.copy(from=system.file("ex/PamplonaDerived",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
pamplona.derived <- read_rtoi(file.path(tempdir(),"PamplonaDerived"))

# print available variables
rsat_list_data(pamplona.derived)

# get RasterStack from raster package
suppressWarnings(mod.ndvi.raster <-
           rsat_get_raster(pamplona.derived, "mod09ga", "NDVI"))
plot(mod.ndvi.raster)

# get spatraster from terra package
mod.ndvi.rast <- rsat_get_SpatRaster(pamplona.derived, "mod09ga", "NDVI")
plot(mod.ndvi.rast)

# get stars from stars package
suppressWarnings(mod.ndvi.stars <-
rsat_get_stars(pamplona.derived, "mod09ga", "NDVI"))
plot(mod.ndvi.stars)


## get any band in rtoi
# list available data
rsat_list_data(pamplona.derived)
# select band 1: MODIS_Grid_500m_2D_sur_refl_b01_1
mod.ndvi.rast <- rsat_get_SpatRaster(pamplona.derived,
                                     "mod09ga",
                                     "MODIS_Grid_500m_2D_sur_refl_b01_1")
plot(mod.ndvi.rast)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extent_crs.R, R/records.R
\name{[,extent_crs,ANY,ANY,ANY-method}
\alias{[,extent_crs,ANY,ANY,ANY-method}
\alias{sub,extent_crs}
\alias{[<-,extent_crs,ANY,ANY,ANY-method}
\alias{sub,extent_crs,extent_crs}
\alias{[,records,ANY,ANY,ANY-method}
\alias{[<-,records,ANY,ANY,ANY-method}
\alias{'[<-',records,records}
\title{Extract or replace parts of an object}
\usage{
\S4method{[}{extent_crs,ANY,ANY,ANY}(x, i)

\S4method{[}{extent_crs,ANY,ANY,ANY}(x, i) <- value

\S4method{[}{records,ANY,ANY,ANY}(x, i)

\S4method{[}{records,ANY,ANY,ANY}(x, i) <- value
}
\arguments{
\item{x}{object from which to extract element(s) or in which to
replace element(s).}

\item{i}{numeric argument. The the position of the element to
select/modify.}

\item{value}{a \code{records} argument. The slot of the records
to be changed.}
}
\value{
returns a selected value
}
\description{
Operators acting on vectors, matrices, arrays and lists to
extract or replace parts.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_search.R
\name{rsat_search}
\alias{rsat_search}
\alias{rsat_search,rtoi,character-method}
\alias{rsat_search,rtoi,character}
\alias{rsat_search,sf,character-method}
\alias{rsat_search,sf,character}
\title{Search satellite images}
\usage{
rsat_search(region, product, ...)

\S4method{rsat_search}{rtoi,character}(region, product, verbose = FALSE, test.mode = FALSE, ...)

\S4method{rsat_search}{sf,character}(region, product, verbose = FALSE, test.mode = FALSE, ...)
}
\arguments{
\item{region}{a \code{Spatial*}, \code{Raster*}, \code{sf} or \code{rtoi}
class objects defining the region of interest.}

\item{product}{a character vector of product names.}

\item{...}{additional arguments for searching}

\item{verbose}{logical argument. If \code{TRUE}, the function prints the
running steps and warnings.}

\item{test.mode}{logical argument. If \code{TRUE}, the function gets test
data from github.}
}
\value{
nothing if x is an rtoi, records class if you search a region.
}
\description{
Search satellite images concerning a particular location, data product, and
date interval. The function returns a \code{records} object if the
\code{region} is a \code{sf}. If an \code{rtoi} is used, the
function returns nothing and the records are added to the \code{rtoi}.
}
\details{
MODIS images are found through the
\href{https://lpdaacsvc.cr.usgs.gov/services/inventory}{NASA Common
 Metadata Repository}
(CMR). The inventory of MODIS products can be found
\href{https://modis.gsfc.nasa.gov/data/dataprod/}{here}.
The catalog shows the product short names and detailed information.
MODIS surface reflectance products are named `mod09ga' and `myd09ga' for
Terra and Aqua satellites. By the time \code{rsat} is
released, NASA carries out the maintenance of its website on Wednesdays,
which may cause an error when connecting to their server.

We use \href{https://scihub.copernicus.eu/}{ESA's powered API} (`SciHub') to
find Sentinel images. The catalog of Sentinel-2 and -3 products can be found
\href{https://sentinel.esa.int/web/sentinel/missions/sentinel-2/data-products}{here}
and
\href{https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-3/data-products}{here},
respectively. Sentinel-2 and -3 surface reflectance product names are
referred to as `S2MSI2A' and `SY_2_SYN___'.

Landsat images are accessed via the
\href{https://m2m.cr.usgs.gov/}{Machine-to-Machine API}.
Details about the Landsat products can be found
\href{https://www.usgs.gov/landsat-missions/product-information}{here}.
The names of Landsat products are `LANDSAT_TM_C1', `LANDSAT_ETM_C1', and
`LANDSAT_8_C1' for missions 4-5, 7, and 8.
}
\examples{
\dontrun{
library(rsat)
set_credentials("username", "password")

# search navarre images using sf
record.list <- rsat_search(
  region = ex.navarre,
  product = "mod09ga",
  dates = as.Date("2011-01-01") + seq(1, 10, 1)
)

# creating a new rtoi
rtoi.path <- tempdir()
navarre <- new_rtoi(
  "Navarre", # name of the region
  ex.navarre, # sf of the region
  rtoi.path
) # path for the rtoi

# see the number of records in navarre
print(navarre)

# search modis images using rtoi
rsat_search(
  region = navarre,
  product = "mod09ga",
  dates = as.Date("2011-01-01") + seq(1, 10, 1)
)

# see the number of records in navarre
print(navarre)

# search landsat images using rtoi
rsat_search(
  region = navarre,
  product = "LANDSAT_8_C1",
  dates = as.Date("2016-01-01") + seq(1, 30, 1)
)

# see the number of records in navarre
print(navarre)

# search sentinel-2 (level 1 and level 2) images using rtoi
rsat_search(
  region = navarre,
  product = c("S2MSI1C", "S2MSI2A"),
  dates = as.Date("2016-01-01") + seq(1, 30, 1)
)

# see the number of records in navarre
print(navarre)

# search sentinel-3 level-2 images using rtoi
rsat_search(
  region = navarre,
  product = "OL_2_LFR___",
  dates = as.Date("2019-01-01") + seq(1, 2, 1)
)

# search sentinel-1 level-2 images using rtoi
rsat_search(
  region = navarre,
  product = "GRD",
  dates = as.Date("2019-01-01") + seq(1, 2, 1)
)

# search Landsat-5 images using rtoi
rsat_search(
  region = navarre,
  product = "LANDSAT_TM_C1",
  dates = as.Date("1988-08-01") + seq(1, 35)
)

print(navarre)

# get all records from rtoi
navarre.records <- records(navarre)

print(navarre.records)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_download.R
\name{rsat_download}
\alias{rsat_download}
\alias{rsat_download,rtoi-method}
\alias{rsat_download,rtoi}
\alias{rsat_download,records-method}
\alias{rsat_download,records}
\title{Download the images from a \code{records} or an \code{rtoi} object}
\usage{
rsat_download(x, ...)

\S4method{rsat_download}{rtoi}(x, db_path, verbose = FALSE, test.mode = FALSE, ...)

\S4method{rsat_download}{records}(x, out.dir, verbose = FALSE, test.mode = FALSE, ...)
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}

\item{...}{additional arguments.}

\item{db_path}{path to the database. By default, the path
is defined by the \code{rtoi}.}

\item{verbose}{logical argument. If \code{TRUE}, the function prints the
running steps and warnings.}

\item{test.mode}{logical argument. If \code{TRUE}, the function gets test
data from github.}

\item{out.dir}{path where the outputs are stored when using a \code{records}.}
}
\value{
nothing. Downloads the images into your database
}
\description{
The function saves the raw images in the database or the specified directory.
It skips the images that already exist in the database or directory.
}
\examples{
\dontrun{
library(rsat)

# create a copy of navarre in temp file
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# assign the path of the database
set_database(file.path(tempdir(),"DATABASE"))
rsat_download(navarre)

rcrds <-  records(navarre)

rsat_download(rcrds)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat_cloud_mask.R
\name{rsat_cloudMask}
\alias{rsat_cloudMask}
\alias{rsat_cloudMask,rtoi-method}
\alias{cloud_mask,rtoi}
\title{Create cloud mask from an rtoi}
\usage{
rsat_cloudMask(x, ...)

\S4method{rsat_cloudMask}{rtoi}(x, products = "ALL", verbose = FALSE, overwrite = FALSE, ...)
}
\arguments{
\item{x}{rtoi object from which cloud masks are computed.}

\item{...}{additional arguments}

\item{products}{the name of the dataset from which cloud masks are computed.}

\item{verbose}{logical argument. If \code{TRUE}, the function prints the
running steps and warnings.}

\item{overwrite}{logical argument. If \code{TRUE}, overwrites the existing
images with the same name.}
}
\value{
nothing. The cloud masks will be save in the hard drive.
Use get_stars to get the variables.
}
\description{
Create cloud mask from an rtoi
}
\examples{
## Smooth data in rtoi
library(rsat)

# create a copy of pamplona in temp file
file.copy(from=system.file("ex/Pamplona",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
pamplona <- read_rtoi(file.path(tempdir(),"Pamplona"))

rsat_cloudMask(pamplona)

rsat_list_data(pamplona)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{new_record}
\alias{new_record}
\alias{new_record,character,character,Date,character,character,character,numeric,numeric,character,character,character,logical,extent_crs-method}
\alias{new_record,}
\alias{character,}
\alias{Date,}
\alias{numeric,}
\alias{logical,}
\alias{extent_crs}
\alias{new_record,character,character,Date,character,character,character,numeric,numeric,character,character,character,logical,missing-method}
\alias{missing}
\title{Create a new \code{records} object}
\usage{
new_record(
  sat,
  name,
  date,
  product,
  download,
  file_path,
  path,
  row,
  tileid,
  preview,
  api_name,
  order,
  extent_crs
)

\S4method{new_record}{character,character,Date,character,character,character,numeric,numeric,character,character,character,logical,extent_crs}(
  sat,
  name,
  date,
  product,
  download,
  file_path,
  path,
  row,
  tileid,
  preview,
  api_name,
  order,
  extent_crs
)

\S4method{new_record}{character,character,Date,character,character,character,numeric,numeric,character,character,character,logical,missing}(
  sat,
  name,
  date,
  product,
  download,
  file_path,
  path,
  row,
  tileid,
  preview,
  api_name,
  order
)
}
\arguments{
\item{sat}{the name of the satellite to which the record belongs.}

\item{name}{the name of the record.}

\item{date}{the date of the record.}

\item{product}{the product.}

\item{download}{the url to download the satellite record.}

\item{file_path}{the saving directory for the satellite record.}

\item{path}{the path of the tiling system.}

\item{row}{the row of the tiling system.}

\item{tileid}{the tile id.}

\item{preview}{the url of the preview of the satellite record.}

\item{api_name}{the api name.}

\item{order}{boolean, defines if the image must be requested or not.}

\item{extent_crs}{extent (used to project the preview).}
}
\value{
records object
}
\description{
Create a new \code{records} object from scratch
}
\examples{
# create a new record from scrach
rcds <- new_record(
  sat = "modis",
  name = "mod09a",
  date = as.Date("2011087", "\%Y\%j"),
  product = "product",
  download = "url/aaa/download",
  file_path = "file_path",
  path = 1,
  row = 1,
  tileid = "exampleid",
  preview = "url",
  api_name = "nasa_inventory",
  order = FALSE
)
rcds

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\name{records}
\alias{records}
\alias{records,rtoi-method}
\alias{records,rtoi}
\alias{records<-}
\alias{records<-,rtoi,records-method}
\alias{records<-,rtoi,records}
\title{Extracts the satellite records}
\usage{
records(x)

\S4method{records}{rtoi}(x)

records(x) <- value

\S4method{records}{rtoi,records}(x) <- value
}
\arguments{
\item{x}{an rtoi object}

\item{value}{a records object to be set to x.}
}
\value{
a set of records in x rtoi
}
\description{
returns the object records from an rtoi.
}
\examples{
#' library(rsat)
# create a copy of navarre
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)
# load example rtoi
navarre <- read_rtoi(file.path(tempdir(),"Navarre"))
print(navarre)

rcrds <- records(navarre)

records(navarre)<-rcrds[1]
print(navarre)

records(navarre) <- rcrds
print(navarre)
unlink(file.path(tempdir(),"Navarre"),recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ex.dem.navarre}
\alias{ex.dem.navarre}
\title{A Digital Elevation Model (DEM) of the region of Navarre (Spain)}
\format{
The \code{RasterStack} contains 6 layers with the same DEM, one for
every image in \code{\link{ex.ndvi.navarre}}.

The \code{RasterStack} coordinates are in the Sinusoidal projection.

\describe{
  \item{name}{layer names contain the capturing date of the
  corresponding image in the format "\code{YYYYJJJ}"}.
  \item{size}{113 rows by 105 columns and 6 layers}.
}
}
\description{
Geographically projected \code{RasterStack} with the digital elevation model
(DEM) of the region of Navarre (Spain). The DEM was obtained from the
\href{http://centrodedescargas.cnig.es/CentroDescargas/locale?request_locale=en}{National Center for Geographic Information}
of Spain. The DEM is used as a covariate in the Image Mean Anomaly (IMA)
algorithm (\code{\link{rsat_smoothing_images}}).
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ex.manhattan}
\alias{ex.manhattan}
\title{A polygon with the border of Manhattan (USA)}
\description{
Spatial feature (\code{sf}) representing the border of Manhattan with
coordinates in the NAD83 format.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R, R/rtoi.R
\name{sat_name}
\alias{sat_name}
\alias{sat_name,records-method}
\alias{sat_name,records}
\alias{sat_name,rtoi-method}
\alias{sat_name,rtoi}
\title{Get the name of the satellite(s) from a \code{records} or an \code{rtoi}}
\usage{
sat_name(x)

\S4method{sat_name}{records}(x)

\S4method{sat_name}{rtoi}(x)
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}
}
\value{
the name of the satellite
}
\description{
Get the name of the satellite(s) from a \code{records} or an \code{rtoi}
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))
# get the records
rcds <- records(navarre)
# coerce the records to dataframe
sat_name(rcds)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R, R/rtoi.R
\name{dates}
\alias{dates}
\alias{dates,records-method}
\alias{dates<-}
\alias{dates<-,records-method}
\alias{dates,rtoi-method}
\title{Get/set the dates from a \code{records} or an \code{rtoi}}
\usage{
dates(x)

\S4method{dates}{records}(x)

dates(x) <- value

\S4method{dates}{records}(x) <- value

\S4method{dates}{rtoi}(x)
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}

\item{value}{the new date to asign}
}
\value{
returns a vector of \code{Date} class
}
\description{
Get/set the dates from a \code{records} or an \code{rtoi}
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get a vector of dates includes in rtoi
dates(navarre)

# get the records
rcds <- records(navarre)

# coerce the records to dataframr
dates(rcds)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{get_download}
\alias{get_download}
\title{Extract the url to download a data record}
\usage{
get_download(x)
}
\arguments{
\item{x}{a \code{records} object.}
}
\value{
download url of a records
}
\description{
It returns a character with the url to download the image.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)
# coerce the records to rtoi
get_download(rcds)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{unique,records,ANY-method}
\alias{unique,records,ANY-method}
\alias{unique}
\title{Extract unique elements}
\usage{
\S4method{unique}{records,ANY}(x)
}
\arguments{
\item{x}{a \code{records} object.}
}
\value{
unique elements in records class
}
\description{
It returns a \code{records} like \code{x} but with duplicate
elements/rows removed.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)

duplicate.records <- c(rcds[1],rcds[1])
length(duplicate.records)
print(duplicate.records)
single.record <- unique(duplicate.records)
length(single.record)
print(single.record)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{get_order}
\alias{get_order}
\alias{get_order<-}
\alias{get_order<-,records-method}
\alias{get_order<-,records}
\title{Get the slot called order from a \code{records} or an \code{rtoi}}
\usage{
get_order(x)

get_order(x) <- value

\S4method{get_order}{records}(x) <- value
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}

\item{value}{logical argument. The new value for \code{x}.}
}
\value{
the value of called order
}
\description{
Get the slot called order from a \code{records} or an \code{rtoi}
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)

# gets a boolean
get_order(rcds)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ex.navarre}
\alias{ex.navarre}
\title{A polygon with the border of Navarre (Spain)}
\description{
Spatial feature (\code{sf}) representing the border of Navarre with
coordinates in the longitude/latitude format.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extent_crs.R, R/records.R
\name{length,extent_crs-method}
\alias{length,extent_crs-method}
\alias{length,records-method}
\title{Length of an object}
\usage{
\S4method{length}{extent_crs}(x)

\S4method{length}{records}(x)
}
\arguments{
\item{x}{a \code{records} object to compute its length.}
}
\value{
Length currently returns a non-negative integer of length 1
}
\description{
Get or set the length of vectors (including lists) and factors,
and of any other R object for which a method has been defined.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)

length(rcds)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\name{rename}
\alias{rename}
\alias{rename,rtoi,character-method}
\alias{rename,rtoi,character}
\title{Renames an \code{rtoi}}
\usage{
rename(x, newname)

\S4method{rename}{rtoi,character}(x, newname)
}
\arguments{
\item{x}{an rtoi object}

\item{newname}{a character class to rename the \code{rtoi}.}
}
\value{
nothing. the  changes the internal name of the rtoi
}
\description{
Renames all parameters and folder name of an \code{rtoi}.
}
\examples{
\dontrun{
myrtoi <- read_rtoi("file_path/rtoir_name")
rename(myrtoi, "Navarre_BACK")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\name{new_rtoi}
\alias{new_rtoi}
\alias{new_rtoi,character,sf,character,missing,missing,missing-method}
\alias{new_rtoi,character,sf,character,missing,missing}
\alias{new_rtoi,character,sf,character,character,missing,missing-method}
\alias{character,sf,character,character}
\alias{new_rtoi,character,sf,character,character,records,missing-method}
\alias{character,sf,character,character,records}
\alias{new_rtoi,character,sf,character,character,records,numeric-method}
\alias{character,sf,character,character,records,size}
\title{Creates a new \code{rtoi} object}
\usage{
new_rtoi(name, region, rtoi_path, db_path, records, size)

\S4method{new_rtoi}{character,sf,character,missing,missing,missing}(name, region, rtoi_path)

\S4method{new_rtoi}{character,sf,character,character,missing,missing}(name, region, rtoi_path, db_path)

\S4method{new_rtoi}{character,sf,character,character,records,missing}(name, region, rtoi_path, db_path, records)

\S4method{new_rtoi}{character,sf,character,character,records,numeric}(name, region, rtoi_path, db_path, records, size)
}
\arguments{
\item{name}{the name of the region of interest.}

\item{region}{an sf object.}

\item{rtoi_path}{the path to the \code{rtoi} folder.}

\item{db_path}{the path to the database.}

\item{records}{a records object.}

\item{size}{the size of \code{rtoi} folder. By default,
the size is computed from \code{rtoi_path}.}
}
\value{
the reference of the \code{rtoi} object
}
\description{
Creates a new \code{rtoi} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R, R/rtoi.R
\name{names,records-method}
\alias{names,records-method}
\alias{names,rtoi-method}
\alias{names,rtoi}
\alias{names<-,rtoi,character-method}
\alias{names<-,rtoi,character}
\title{Get the name of the object}
\usage{
\S4method{names}{records}(x)

\S4method{names}{rtoi}(x)

\S4method{names}{rtoi,character}(x) <- value
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}

\item{value}{character argument. The new value for \code{x}.}
}
\value{
a character vector containing the name of
all the names in \code{x}.
}
\description{
A function to get or set the names of an object.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

names(navarre)
names(navarre) <- "New name"
names(navarre)

rcrds <- records(navarre)

names(rcrds)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testing_function.R
\name{test_function}
\alias{test_function}
\title{Testing function}
\usage{
test_function()
}
\description{
Function used for testing some internal functions in continuous integration.
}
\examples{
test_function()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extent_crs.R, R/records.R
\name{c,extent_crs-method}
\alias{c,extent_crs-method}
\alias{c}
\alias{c,records-method}
\title{Combine values into a vector or a list}
\usage{
\S4method{c}{extent_crs}(x, ...)

\S4method{c}{records}(x, ...)
}
\arguments{
\item{x}{a \code{records} object.}

\item{...}{additional arguments.}
}
\value{
a combination of 'x' class elements
}
\description{
This is a generic function which combines its arguments.
}
\details{
The default method combines its arguments to form a vector.
All arguments are coerced to a common type which is the type
 of the returned value. All attributes except names are removed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connections.R
\name{set_credentials}
\alias{set_credentials}
\alias{set_credentials,character,character,missing-method}
\alias{set_credentials,character,character,missing}
\alias{set_credentials,character,character,character-method}
\alias{set_credentials,character,character,character}
\title{Saves the credentials for the web services}
\usage{
set_credentials(user, pass, credential)

\S4method{set_credentials}{character,character,missing}(user, pass)

\S4method{set_credentials}{character,character,character}(user, pass, credential)
}
\arguments{
\item{user}{character argument. Defines the username of an api platform
to search or download images}

\item{pass}{character argument. Defines the password of an api platform
to search and download images}

\item{credential}{optional argument to specify the name of the platform.
Valid names are earthdata, scihub, scihubs5p, or ALL}
}
\value{
nothing. set the credentials in the package environment variable
}
\description{
Saves the credentials for the web services
}
\examples{
print_credentials()
set_credentials("example", "example")
print_credentials()
set_credentials("example", "example", "earthdata")
print_credentials()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{subset,records-method}
\alias{subset,records-method}
\title{Filter the satellite records of a \code{records} or an \code{rtoi}}
\usage{
\S4method{subset}{records}(x, subset, select)
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}

\item{subset}{character argument indicating the name of the slot.}

\item{select}{character with the value for subsetting.}
}
\value{
filtered records class
}
\description{
Filter the satellite records of a \code{records} or an \code{rtoi}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsat.R
\docType{package}
\name{rsat}
\alias{rsat}
\title{`rsat'}
\description{
The goal of `rsat` is to help you handling time-series of satellite
images from multiple platforms in a local, efficient, and standardized
way. The package provides tools to;
\enumerate{
  \item Search (run \code{vignette("rsat1_search", package = "rsat")})
  \item Download (run \code{vignette("rsat2_download", package = "rsat")})
  \item Customize, and (run \code{vignette("rsat3_customize", package = "rsat")})
  \item Process (run \code{vignette("rsat4_process", package = "rsat")})
}
satellite images from Landsat, MODIS, and Sentinel for a region and time
of interest.
}
\section{Registration}{

The registration in the following online portals is required to get a
full access to satellite images with `rsat`;
\itemize{
\item \href{https://ers.cr.usgs.gov/register}{USGS} USGS is the sole
science agency for the Department of the Interior of United States.
Provide access to Modis Images. More information about EarthData
can be found \href{https://www.usgs.gov}{Here}.
\item \href{https://urs.earthdata.nasa.gov}{EarthData}: A repository
of NASA's earth observation data-sets. More information about
EarthData can be found
\href{https://earthdata.nasa.gov/earth-observation-data}{here}.
\item \href{https://scihub.copernicus.eu/dhus/#/self-registration}{SciHub},
a web service giving access to Copernicus' scientific data hub.
Please go \href{https://scihub.copernicus.eu}{here} to find more
details about the data hub.
}
For convenience, try to use the same username and password for all of them.
To satisfy the criteria of all web services make sure that the username is 4
 characters long and includes a period, number or underscore. The password must
 be 12 character long and should include characters with at least one capital
 letter, and numbers.
}

\section{Contributing}{

If you want to contribute by adding new features or fixing bugs in the package you can do it from our github address:
\url{https://github.com/spatialstatisticsupna/rsat}
Bug report: \url{https://github.com/spatialstatisticsupna/rsat/issues}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\name{region}
\alias{region}
\alias{region,rtoi-method}
\alias{region,rtoi}
\alias{region<-}
\alias{region<-,rtoi}
\alias{region<-,rtoi,sf-method}
\alias{region<-,rtoi,sf}
\alias{region<-,rtoi,NULL-method}
\alias{region<-,rtoi,NULL}
\title{Extracts region from an rtoi}
\usage{
region(x)

\S4method{region}{rtoi}(x)

region(x) <- value

\S4method{region}{rtoi,sf}(x) <- value

\S4method{region}{rtoi,`NULL`}(x) <- value
}
\arguments{
\item{x}{an rtoi object.}

\item{value}{an sf object to define the region in x.}
}
\value{
the sf class with the region of an rtoi
}
\description{
gets the sf that specifies the region of an rtoi.
}
\examples{
library(rsat)
# create a copy of navarre
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the region from rtoi
sf.obj <-  region(navarre)
plot(sf.obj)

# asign new region value
region(navarre)<-NULL

region(navarre)<-sf.obj
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variables.R
\name{show_variables}
\alias{show_variables}
\alias{show_variables,ANY-method}
\alias{show_variables-method}
\title{List the variables and satellites supported by \code{rsat}}
\usage{
show_variables(...)

\S4method{show_variables}{ANY}()
}
\arguments{
\item{...}{arguments for nestering functions}
}
\value{
prints supported satellites and derived variables information.
}
\description{
Displays the satellites and variable method
}
\examples{
show_variables()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\name{get_database}
\alias{get_database}
\alias{get_database,rtoi-method}
\alias{get_database,rtoi}
\alias{get_database,missing-method}
\alias{set_database}
\alias{set_database,rtoi-method}
\alias{set_database,character-method}
\title{Extracts or assign the path of the database}
\usage{
get_database(x)

\S4method{get_database}{rtoi}(x)

\S4method{get_database}{missing}()

set_database(x, ...)

\S4method{set_database}{rtoi}(x, value)

\S4method{set_database}{character}(x)
}
\arguments{
\item{x}{an rtoi object.}

\item{...}{additional arguments.}

\item{value}{character argument. The value for
change the database directory of x.}
}
\value{
the database path of an rtoi
}
\description{
Extracts the path to the database from an rtoi/package environment.
If both, environment and rtoi database are defined the rtoi
database is used.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the databse used by navarre
get_database(navarre)

# set the a new database path
set_database(navarre,"new_path")

# get the database used by rsat by default
get_database()

# set the a new database path for the entire environment
set_database("new_path")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtoi.R
\name{read_rtoi}
\alias{read_rtoi}
\alias{read_rtoi,character-method}
\alias{read_rtoi,character}
\title{Reads an rtoi from the hard drive}
\usage{
read_rtoi(path, ...)

\S4method{read_rtoi}{character}(path, ...)
}
\arguments{
\item{path}{an rtoi object.}

\item{...}{additional arguments.}
}
\value{
rtoi object readed from disk.
}
\description{
Reads an rtoi from the hard drive
}
\examples{
library(rsat)

# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))
print(navarre)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R, R/rtoi.R
\name{get_dir}
\alias{get_dir}
\alias{get_dir,records-method}
\alias{get_dir,records}
\alias{get_order,records-method}
\alias{get_dir,rtoi-method}
\alias{get_dir,rtoi}
\title{Get the file path of a \code{records} or an \code{rtoi}}
\usage{
get_dir(x)

\S4method{get_dir}{records}(x)

\S4method{get_order}{records}(x)

\S4method{get_dir}{rtoi}(x)
}
\arguments{
\item{x}{.}
}
\value{
the file path in the records
}
\description{
Get the file path of a \code{records} or an \code{rtoi}
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the path of the
get_dir(navarre)

# get the records
rcds <- records(navarre)

# gets the relative path to store records data
get_dir(rcds)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot,rtoi,Date-method}
\alias{plot,rtoi,Date-method}
\alias{plot,rtoi,character-method}
\alias{plot,character}
\alias{plot,records,ANY-method}
\alias{plot,records}
\alias{plot,rtoi,missing-method}
\alias{plot,rtoi,missing}
\title{Plot an \code{rtoi} object}
\usage{
\S4method{plot}{rtoi,Date}(
  x,
  y,
  ...,
  variable = "rgb",
  band_name = c("red", "green", "blue"),
  verbose = FALSE,
  xsize = 250,
  ysize = 250
)

\S4method{plot}{rtoi,character}(
  x,
  y,
  ...,
  variable = "rgb",
  product = "ALL",
  band_name = c("red", "green", "blue"),
  dates = NULL,
  verbose = FALSE,
  xsize = 250,
  ysize = 250
)

\S4method{plot}{records,ANY}(x, y, verbose = FALSE, ...)

\S4method{plot}{rtoi,missing}(x, y, verbose = FALSE, ...)
}
\arguments{
\item{x}{an \code{rtoi} or \code{records}.}

\item{y}{character argument. The valid values are "dates", "preview", or
"view".}

\item{...}{additional arguments.}

\item{variable}{character argument. The variable to be plotted. By default,
a color (RGB) variable is selected .}

\item{band_name}{character vector argument. Enables false color plots. By
default, usual bands are selected \code{c("red","green","blue")}.}

\item{verbose}{logical argument. If \code{TRUE}, the function prints the
running steps and warnings.}

\item{xsize}{the number of samples on the horizontal axis.}

\item{ysize}{the number of samples on the vertical axis.}

\item{product}{character argument. The product name to be plotted.}

\item{dates}{date vector argument. The dates to be plotted.}
}
\value{
\code{tmap} plot.
}
\description{
Plot (a map of) the values of an \code{rtoi} or \code{records} object.
}
\examples{
library(rsat)
 \dontrun{

# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

print(navarre)

# plot the calendar
plot(navarre, "dates")



# replace with your own "username" and "password"
set_credentials("username", "password")

# plot the quicklook images before the download
# needs credentials to download preview images
plot(navarre, y = "preview")

# select partially cloud free
rcds <- records(navarre)
rcds <- rcds[dates(rcds) \%in\% as.Date(c("20210310", "20210313"), "\%Y\%m\%d")]
records(navarre) <- rcds

plot(navarre, "preview")

file.copy(from=system.file("ex/Pamplona",package="rsat"),
         to=tempdir(),
         recursive = TRUE)
# plot already mosaicked rtoi ("view" mode)
pamplona <- read_rtoi(file.path(tempdir(),"Pamplona"))

rsat_list_data(pamplona)

# plot can compute the rgb image on the fly from mosaicek bands
plot(pamplona, "view", product="mod09ga")

# plot on the fly with false color
plot(pamplona, "view",
     product = "mod09ga",
     band_name = c("nir", "red", "green"))

file.copy(from=system.file("ex/PamplonaDerived",package="rsat"),
         to=tempdir(),
         recursive = TRUE)
# plot already mosaicked rtoi ("view" mode)
pamplona.derived <- read_rtoi(file.path(tempdir(),"PamplonaDerived"))

rsat_list_data(pamplona.derived)

# plot derived variables
plot(pamplona.derived, "view",
     product = "mod09ga",
     variable = "NDVI")

# Set the max and min value in plot
plot(pamplona.derived,"view",
     variable="NDVI",
     product="mod09ga",
     zlim=c(0,1))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\docType{class}
\name{records-class}
\alias{records-class}
\title{A class object for satellite image metadata}
\description{
This class object organizes the attributes of satellite images' metadata
from several missions/programs uniformly. Structuring the information
facilitates managing, previewing, and downloading data records.
}
\details{
\code{records} works as vector. It accepts usual R methods such as
\code{c}, \code{[]}, \code{length()}, \code{subset()} or \code{unique()}.
Each record (vector element) contains several parameters or slots.

The object can be coerced into a \code{data.frame} by
using the function \code{as.data.frame()}. The \code{data.frame} can
be transformed back into a \code{records} with the function
\code{as.records()}.
}
\section{Slots}{

\describe{
\item{\code{sat}}{the name of the satellite.}

\item{\code{name}}{the name of the file.}

\item{\code{date}}{capturing date of the image.}

\item{\code{product}}{name of the data product.}

\item{\code{path}}{the path of the tiling system.}

\item{\code{row}}{the row of the tiling system.}

\item{\code{tileid}}{the tile identification number.}

\item{\code{download}}{the download url.}

\item{\code{file_path}}{the saving directory for the satellite record.}

\item{\code{preview}}{the preview url.}

\item{\code{api_name}}{the name of the API.}

\item{\code{order}}{boolean, whether the image needs to be ordered.}

\item{\code{extent_crs}}{coordinate reference system of the preview.}
}}

\examples{
\dontrun{
library(rsat)
# create a copy of navarre
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

# load example rtoi
navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

rcrds <- records(navarre)

modis.rcrds <- rcrds[sat_name(rcrds)\%in\%"Modis"]
ls8.rcrds <- rcrds[sat_name(rcrds)\%in\%"Landsat-8"]

class(modis.rcrds)
dates(ls8.rcrds)
modis.ls8.records <- c(ls8.rcrds, modis.rcrds)

print(modis.ls8.records)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ex.madrid}
\alias{ex.madrid}
\title{A polygon with the border of Madrid (Spain)}
\description{
Spatial feature (\code{sf}) representing the border of Madrid with
coordinates in the longitude/latitude format.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R, R/rtoi.R
\name{product}
\alias{product}
\alias{product,records-method}
\alias{product,records}
\alias{product,rtoi-method}
\alias{product,rtoi}
\title{Get the name of the product from a \code{records} or an \code{rtoi}}
\usage{
product(x)

\S4method{product}{records}(x)

\S4method{product}{rtoi}(x)
}
\arguments{
\item{x}{a \code{records} or an \code{rtoi} object.}
}
\value{
the name of the product in the records
}
\description{
Get the name of the product from a \code{records} or an \code{rtoi}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{get_api_name}
\alias{get_api_name}
\alias{get_api_name,records-method}
\alias{get_api_name,records}
\title{Get the API name of a \code{records}}
\usage{
get_api_name(x)

\S4method{get_api_name}{records}(x)
}
\arguments{
\item{x}{a \code{records} object.}
}
\value{
a character vector containing the API names of the
elements in \code{x}.
}
\description{
A function to get or set the api names of an object.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)

# get a vector with the api name of each records
get_api_name(rcds)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extent_crs.R, R/records.R, R/rtoi.R,
%   R/variables.R
\name{show,extent_crs-method}
\alias{show,extent_crs-method}
\alias{show}
\alias{show,records-method}
\alias{show,rtoi-method}
\alias{show,rtoi}
\alias{show,variables-method}
\alias{show,records}
\title{Show an Object}
\usage{
\S4method{show}{extent_crs}(object)

\S4method{show}{records}(object)

\S4method{show}{rtoi}(object)

\S4method{show}{variables}(object)
}
\arguments{
\item{object}{Any R object}
}
\value{
show returns an invisible NULL.
}
\description{
Display the object, by printing, plotting or whatever suits its class. This
function exists to be specialized by methods.
The default method calls showDefault.
}
\examples{
## load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

## The method will now be used for automatic printing of navarre
navarre

## get records
rcds <- records(navarre)

rcds
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records.R
\name{get_preview}
\alias{get_preview}
\alias{get_preview,records-method}
\alias{get_preview,records}
\alias{get_download,records-method}
\title{Extract the url of the preview}
\usage{
get_preview(x)

\S4method{get_preview}{records}(x)

\S4method{get_download}{records}(x)
}
\arguments{
\item{x}{a \code{records} object.}
}
\value{
preview url of a records
}
\description{
It returns a character vector of urls to preview the data records.
}
\examples{
# load example rtoi
file.copy(from=system.file("ex/Navarre",package="rsat"),
         to=tempdir(),
         recursive = TRUE)

navarre <- read_rtoi(file.path(tempdir(),"Navarre"))

# get the records
rcds <- records(navarre)

# get a vector with the preview url of each record
get_api_name(rcds)
}
