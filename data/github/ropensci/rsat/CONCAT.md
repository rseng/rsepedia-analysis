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
