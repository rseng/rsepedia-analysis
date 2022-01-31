
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tiler <img src="man/figures/logo.png" style="margin-left:10px;margin-bottom:5px;" width="120" align="right">

**Author:** [Matthew Leonawicz](https://github.com/leonawicz)
<a href="https://orcid.org/0000-0001-9452-2771" target="orcid.widget">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
<br/> **License:** [MIT](https://opensource.org/licenses/MIT)<br/>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Travis build
status](https://travis-ci.org/ropensci/tiler.svg?branch=master)](https://travis-ci.org/ropensci/tiler)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/tiler?branch=master&svg=true)](https://ci.appveyor.com/project/leonawicz/tiler)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/tiler/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/tiler?branch=master)

[![](https://badges.ropensci.org/226_status.svg)](https://github.com/ropensci/software-review/issues/226)
[![CRAN
status](http://www.r-pkg.org/badges/version/tiler)](https://cran.r-project.org/package=tiler)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/tiler)](https://cran.r-project.org/package=tiler)
[![Github
Stars](https://img.shields.io/github/stars/ropensci/tiler.svg?style=social&label=Github)](https://github.com/ropensci/tiler)

![](https://github.com/ropensci/tiler/blob/master/data-raw/ne.jpg?raw=true)

## Create geographic and non-geographic map tiles

The `tiler` package provides a tile generator function for creating map
tile sets for use with packages such as `leaflet`. In addition to
generating map tiles based on a common raster layer source, it also
handles the non-geographic edge case, producing map tiles from arbitrary
images. These map tiles, which have a a non-geographic simple coordinate
reference system, can also be used with `leaflet` when applying the
simple CRS option.

Map tiles can be created from an input file with any of the following
extensions: `tif`, `grd` and `nc` for spatial maps and `png`, `jpg` and
`bmp` for basic images.

## Motivation

This package helps R users who wish to create geographic and
non-geographic map tiles easily and seamlessly with only a single line
of R code. The intent is to do this with a package that has

  - minimal heavy package dependencies.
  - minimal extraneous general features and functions that do not have
    to do with tile generation.
  - to create tiles without having to code directly in other software,
    interact directly with Python, or make calls at the command line;
    allowing the R user to remain comfortably within the familiar R
    environment.
  - to support the creation on map tiles from raw images for users who
    wish to create non-standard maps, which may also be followed by
    georeferencing locations of interest in the simplified coordinate
    reference system of the map image.

## Installation

Install `tiler` from CRAN with

``` r
install.packages("tiler")
```

Install the development version from GitHub with

``` r
# install.packages("remotes")
remotes::install_github("ropensci/tiler")
```

For non-geographic tiles, using a `png` file is recommended for quality
and file size. `jpg` may yield a lower quality result, while a large,
high resolution `bmp` file may have an enormous file size compared to
`png`.

`jpg` and `bmp` are *optionally* supported by `tiler`. This means they
are not installed and imported with `tiler`. It is assumed the user will
provide `png` images. If using `jpg` or `bmp` and the packages `jpeg` or
`bmp` are not installed, respectively, `tile` will print a message to
the console notifying of the required package installations.

### System requirements

This package requires Python and the `gdal` library for Python. Windows
users are recommended to install
[OSGeo4W](https://trac.osgeo.org/osgeo4w/) as an easy way to obtain the
required `gdal` support for Python in Windows. See `?tiler_options` or
the package vignette for more information.

-----

Please note that the `tiler` project is released with a [Contributor
Code of
Conduct](https://github.com/ropensci/tiler/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# tiler 0.2.5

* Documentation updates.

# tiler 0.2.4

* Precompiling of vignette depending on external data.

# tiler 0.2.3

* Minor fixes for CRAN release.

# tiler 0.2.2

* Breaking change: no longer using a `format` argument. All tiles are TMS.
* Updated `gdal2tiles` to version 2.4 release.
* Bug fix.
* Updated documentation.

# tiler 0.2.1

* Improved and simplified instructions and expectations for Windows use. Windows users must add `OSGeo4W.bat` path to system path.
* Added `leaflet` examples of remotely hosted tiles generated by `tiler` to vignette.
* Bug fix related to system path to `OSGeo4W.bat` being ignored on Windows.

# tiler 0.2.0

* All three `gdal2tiles*` scripts have been updated to accept a command line argument when called by R that provides a path for any temporary files, i.e., `tmp.*.vrt` files created by the `gdal2tiles*` scripts. These were previously accumulating in the system temp folder. The new temporary directory is a sub-directory inside `tempdir()`. Therefore, it is cleaned up when exiting R. Nevertheless, `tile` also force deletes this subdirectory immediately after its internal system call to one of the `gdal2tiles*` scripts returns, so the temporary sub-directory does not even exist for the full duration of the `tile` call.
* Added functions `tile_viewer` and `view_tiles` and other supporting functions for generating HTML Leaflet tile preview web page.
* Added arguments to `tile`. `tile` now generates previewer by default.
* Added unit tests.
* Updated vignette.

# tiler 0.1.6

* Made minor formatting changes per CRAN request for resubmission.

# tiler 0.1.5

* Refactored `tile`, added arguments including `resume` and `format`, changed some argument names.
* Added default support for XYZ format tiles in addition to TMS. This brings in another version of `gdal2tiles`.
* Updated documentation.
* Added unit tests.

# tiler 0.1.0

* Created `tile` function for generating map tiles from geographic or non-geographic maps/images.
* Created readme with basic example.
* Added initial introduction vignette.
* Added robust unit tests and other external examples and spot testing of edge cases.

# tiler 0.0.0.9000

* Added initial package scaffolding.
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
## Test environments

* local Windows 10 install, R 4.0.3
* Windows 10 (AppVeyor), R 4.0.3
* Ubuntu 16.04 (Travis CI), R-devel, R-release, R-oldrel
* Mac OSX (Travis CI) R-release (4.0.3)
* win-builder

## Update release

* This update includes a maintainer email address update.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

Checks pass on all packages. (https://github.com/ropensci/tiler/tree/master/revdep)

## System Requirements details

This package requies Python and python-gdal. These are specified in the DESCRIPTION System Requirements field as well as the Description text.
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 3.6.1 (2019-07-05) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |English_United States.1252   |
|ctype    |English_United States.1252   |
|tz       |America/Denver               |
|date     |2019-11-25                   |

# Dependencies

|package |old   |new   |<U+0394>  |
|:-------|:-----|:-----|:--|
|tiler   |0.2.3 |0.2.4 |*  |
|png     |0.1-7 |0.1-7 |   |
|raster  |3.0-7 |3.0-7 |   |
|Rcpp    |1.0.3 |1.0.3 |   |
|rgdal   |1.4-7 |1.4-7 |   |
|sp      |1.3-2 |1.3-2 |   |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.path = "man/figures/README-",
  message = FALSE, warning = FALSE, error = FALSE
)
```

# tiler <img src="man/figures/logo.png" style="margin-left:10px;margin-bottom:5px;" width="120" align="right">
**Author:** [Matthew Leonawicz](https://github.com/leonawicz) <a href="https://orcid.org/0000-0001-9452-2771" target="orcid.widget">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Travis build status](https://travis-ci.org/ropensci/tiler.svg?branch=master)](https://travis-ci.org/ropensci/tiler)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/tiler?branch=master&svg=true)](https://ci.appveyor.com/project/leonawicz/tiler)
[![Codecov test coverage](https://codecov.io/gh/ropensci/tiler/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/tiler?branch=master)

[![](https://badges.ropensci.org/226_status.svg)](https://github.com/ropensci/software-review/issues/226)
[![CRAN status](http://www.r-pkg.org/badges/version/tiler)](https://cran.r-project.org/package=tiler)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/tiler)](https://cran.r-project.org/package=tiler)
[![Github Stars](https://img.shields.io/github/stars/ropensci/tiler.svg?style=social&label=Github)](https://github.com/ropensci/tiler)

![](https://github.com/ropensci/tiler/blob/master/data-raw/ne.jpg?raw=true)

## Create geographic and non-geographic map tiles

The `tiler` package provides a tile generator function for creating map tile sets for use with packages such as `leaflet`. In addition to generating map tiles based on a common raster layer source, it also handles the non-geographic edge case, producing map tiles from arbitrary images. These map tiles, which have a a non-geographic simple coordinate reference system, can also be used with `leaflet` when applying the simple CRS option. 

Map tiles can be created from an input file with any of the following extensions: `tif`, `grd` and `nc` for spatial maps and `png`, `jpg` and `bmp` for basic images.

## Motivation

This package helps R users who wish to create geographic and non-geographic map tiles easily and seamlessly with only a single line of R code. The intent is to do this with a package that has

*    minimal heavy package dependencies.
*    minimal extraneous general features and functions that do not have to do with tile generation.
*    to create tiles without having to code directly in other software, interact directly with Python, or make calls at the command line; allowing the R user to remain comfortably within the familiar R environment.
*    to support the creation on map tiles from raw images for users who wish to create non-standard maps, which may also be followed by georeferencing locations of interest in the simplified coordinate reference system of the map image.

## Installation

Install `tiler` from CRAN with

``` r
install.packages("tiler")
```

Install the development version from GitHub with

``` r
# install.packages("remotes")
remotes::install_github("ropensci/tiler")
```

For non-geographic tiles, using a `png` file is recommended for quality and file size. `jpg` may yield a lower quality result, while a large, high resolution `bmp` file may have an enormous file size compared to `png`. 

`jpg` and `bmp` are *optionally* supported by `tiler`. This means they are not installed and imported with `tiler`. It is assumed the user will provide `png` images. If using `jpg` or `bmp` and the packages `jpeg` or `bmp` are not installed, respectively, `tile` will print a message to the console notifying of the required package installations.

### System requirements

This package requires Python and the `gdal` library for Python.
Windows users are recommended to install [OSGeo4W](https://trac.osgeo.org/osgeo4w/) as an easy way to obtain the required `gdal` support for Python in Windows. See `?tiler_options` or the package vignette for more information.

---

Please note that the `tiler` project is released with a [Contributor Code of Conduct](https://github.com/ropensci/tiler/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction to tiler"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to tiler}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



The `tiler` package provides a map tile-generator function for creating map tile sets for use with packages such as `leaflet`.

In addition to generating map tiles based on a common raster layer source, it also handles the non-geographic edge case, producing map tiles from arbitrary images. These map tiles, which have a "simple CRS", a non-geographic simple Cartesian coordinate reference system, can also be used with `leaflet` when applying the simple CRS option. Map tiles can be created from an input file with any of the following extensions: `tif`, `grd` and `nc` for spatial maps and `png`, `jpg` and `bmp` for basic images.

## Setup (Windows users)

Python is required as well as the `gdal` library for Python. For Windows users, this is not the same as simply having `rgdal` installed through R. You need to install `gdal` so that it is accessible by Python. One of the easiest and typically recommended ways to do this in Windows is to install [OSGeo4W](https://trac.osgeo.org/osgeo4w/). It will bring all the required Python `gdal` library functionality along with it. OSGeo4W is also commonly installed along with [QGIS](https://qgis.org/en/site/forusers/download.html).

On Windows, `tiler_options` is set on package load with `osgeo4w = "OSgeo4W.bat"`. Make sure to add the path to this file to your system path variable. Otherwise it will not be found when calling `tile`. You can set the full path to `OSGeo4W.bat` directly in your R session with `tiler_options`. However, it is recommended to add the directory (e.g. `C:/OSGeo64` or wherever `OSGeo4W.bat` is located) to your system path so that you never had to deal with it in R.

Linux and Mac users should not have to do any additional setup as long as Python and `gdal` for Python are installed and available.

## Geographic map tiles

### Context

For the sake of simple, expedient examples, the map tiles generated below all use small zoom level ranges. There is also no reason to attempt displaying the tiles here. To make these examples more informative, each raster is loaded and plotted for context, though this is not necessary to the tiling process. Loading the `raster` package is only needed here for the print and plot calls.

The example maps packaged with `tiler` are not representative of large, high resolution imagery that benefits from tiling. These maps are very small in order to minimize package size and ensure examples run quickly. But the tiling procedures demonstrated are the same as would be applied to larger images.

Lastly, consider the power of your system before attempting to make a ton of tiles for large images at very high resolutions. You could find that the system could hang at any one of a number of choke points. If you are attempting to make thousands of tiles for a large, high resolution image and your system is struggling, it is recommended to (1) try making tiles for only one zoom level at a time, starting from zero and then increasing while monitoring your system resources. (2) If this is not enough, find a better system.

### Basic example

Map tiles are generated with `tile`. `tile` takes an input file name `file` for the source map and a `tiles` destination path for where to save map tiles. The only other required argument is `zoom`, the range of zoom levels for the tiles. This is a string of the form, e.g. `3:7`. In this and the subsequent examples zoom levels 0-3 are used.


```r
library(tiler)
library(raster)
tile_dir <- file.path(tempdir(), "tiles")
map <- system.file("maps/map_wgs84.tif", package = "tiler")
(r <- raster(map))
#> class      : RasterLayer 
#> dimensions : 32, 71, 2272  (nrow, ncol, ncell)
#> resolution : 0.8333333, 0.8333333  (x, y)
#> extent     : -125.0208, -65.85417, 23.27083, 49.9375  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
#> source     : C:/github/tiler/inst/maps/map_wgs84.tif 
#> names      : map_wgs84 
#> values     : -0.7205295, 5.545086  (min, max)
plot(r)
```

![plot of chunk ex1](vignette-ex1-1.png)

```r

tile(map, tile_dir, "0-3")
#> Coloring raster...
#> Preparing for tiling...
#> [1] "C:\\Users\\Matt\\AppData\\Local\\Temp\\RtmpQLDRP8"
#> Creating tiles. Please wait...
#> Creating tile viewer...
#> Complete.

list.files(tile_dir)
#> [1] "0"                   "1"                   "2"                   "3"                   "doc.kml"             "tilemapresource.xml"
```



Listing the files in `tile_dir` shows the top level map tiles directories, 0-3 as expected. This is not printed in subsequent examples since it is not going to change.

Note that these examples rendered to HTML here do not capture the parts of the log output that result from the internal system call made by `tile`. When you run this example yourself you will see a bit more information at the console.

### Projected maps

The previous example used a map with a geospatial coordinate reference system (CRS) but it was not projected. That map would be ready to view with the `leaflet` package for example, as would the tiles generated based on it. The next example uses a similar map that is projected. In order to generate map tiles that can be used with `leaflet` in the standard CRS, the map must be reprojected first. Then the same map tiles are generated as before.


```r
map <- system.file("maps/map_albers.tif", package = "tiler")
(r <- raster(map))
#> class      : RasterLayer 
#> dimensions : 32, 71, 2272  (nrow, ncol, ncell)
#> resolution : 85011.74, 103363.8  (x, y)
#> extent     : -2976297, 3059536, -1577300, 1730342  (xmin, xmax, ymin, ymax)
#> crs        : +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
#> source     : C:/github/tiler/inst/maps/map_albers.tif 
#> names      : map_albers 
#> values     : -0.593978, 5.544724  (min, max)
plot(r)
```

![plot of chunk ex2](vignette-ex2-1.png)

```r

tile(map, tile_dir, "0-3")
#> Reprojecting raster...
#> Coloring raster...
#> Preparing for tiling...
#> [1] "C:\\Users\\Matt\\AppData\\Local\\Temp\\RtmpQLDRP8"
#> Creating tiles. Please wait...
#> Creating tile viewer...
#> Complete.
```



The tiles generated this time are the same as before, that is, ready for `leaflet`. `tile` reprojects the raster layer internally. This can be seen in the log output printed to the console.

### Missing CRS

If the CRS of the raster is `NA`, there are two options. By default, `tile` will fall back on processing the raster layer as if it were a simple image file with no geospatial projection information (see the next section on simple CRS/non-geographic map tiles). These tiles are not the same as the previous sets.


```r
map <- system.file("maps/map_albers_NA.tif", package = "tiler")
(r <- raster(map))
#> class      : RasterLayer 
#> dimensions : 32, 71, 2272  (nrow, ncol, ncell)
#> resolution : 85011.74, 103363.8  (x, y)
#> extent     : -2976297, 3059536, -1577300, 1730342  (xmin, xmax, ymin, ymax)
#> crs        : NA 
#> source     : C:/github/tiler/inst/maps/map_albers_NA.tif 
#> names      : map_albers_NA 
#> values     : -0.593978, 5.544724  (min, max)
plot(r)
```

![plot of chunk ex3](vignette-ex3-1.png)

```r

tile(map, tile_dir, "0-3")
#> Warning in .proj_check(file, crs, ...): Projection expected but is missing. Continuing as non-geographic image.
#> Coloring raster...
#> Preparing for tiling...
#> [1] "C:\\Users\\Matt\\AppData\\Local\\Temp\\RtmpQLDRP8"
#> Creating tiles. Please wait...
#> Creating tile viewer...
#> Complete.
```



This is not likely what is wanted. The other option is to force set a known CRS if it is missing from the file or was dropped for whatever reason. Now reprojection can proceed and the expected tiles are generated.


```r
crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"
tile(map, tile_dir, "0-3", crs)
#> Reprojecting raster...
#> Coloring raster...
#> Preparing for tiling...
#> [1] "C:\\Users\\Matt\\AppData\\Local\\Temp\\RtmpQLDRP8"
#> Creating tiles. Please wait...
#> Creating tile viewer...
#> Complete.
```



A note on reprojection: Depending on the nature of the data in a raster, the `...` argument to `tile` allows you to pass through the `method` argument to `raster::projectRaster`. This is `bilinear` by default for bilinear interpolation, appropriate for continuous data. It can be set to `ngb` for nearest neighbor, appropriate for discrete or categorical data. If more control is needed over the reprojection, you should just prepare your raster file first before using `tile`. `tiler` is not intended to substitute for or wrap general geospatial processing tasks that can easily be done with other packages.

### Coloring tiles

Being able to change the default color palette or specify color breaks is important. All other optional `...` arguments to `tile` are passed to `raster::RGB` to provide control over the creation of an intermediary RGB raster. Most arguments to `RGB` can usually be ignored. The most useful ones are `col` and `colNA` for the data values color palette and the `noData` color, respectively. Coloring tiles differently for the original example is as simple as the following.


```r
map <- system.file("maps/map_wgs84.tif", package = "tiler")
pal <- colorRampPalette(c("darkblue", "lightblue"))(20)
nodata <- "tomato"
tile(map, tile_dir, "0-3", col = pal, colNA = nodata)
#> Coloring raster...
#> Preparing for tiling...
#> [1] "C:\\Users\\Matt\\AppData\\Local\\Temp\\RtmpQLDRP8"
#> Creating tiles. Please wait...
#> Creating tile viewer...
#> Complete.
```



### RGB and RGBA rasters

Multi-band rasters are supported as long as they have three or four layers, in which case `tile` assumes these represent red, green, blue and alpha channel, respectively. Internally, single-layer raster files are colored and converted to a three- or four-band RGB/A raster object prior to tile generation. If `file` is already such a raster, this step is simply skipped. Optional arguments like data and `noData` color, break points, etc., are ignored since this type of raster contains its own color information.


```r
map <- system.file("maps/map_albers_rgb.tif", package = "tiler")
(r <- brick(map))
#> class      : RasterBrick 
#> dimensions : 32, 71, 2272, 3  (nrow, ncol, ncell, nlayers)
#> resolution : 85011.74, 103363.8  (x, y)
#> extent     : -2976297, 3059536, -1577300, 1730342  (xmin, xmax, ymin, ymax)
#> crs        : +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
#> source     : C:/github/tiler/inst/maps/map_albers_rgb.tif 
#> names      : map_albers_rgb.1, map_albers_rgb.2, map_albers_rgb.3 
#> min values :                0,                0,                0 
#> max values :              253,              255,              244
plot(r)
```

![plot of chunk ex6](vignette-ex6-1.png)

```r

tile(map, tile_dir, "0-3")
#> Reprojecting raster...
#> Preparing for tiling...
#> [1] "C:\\Users\\Matt\\AppData\\Local\\Temp\\RtmpQLDRP8"
#> Creating tiles. Please wait...
#> Creating tile viewer...
#> Complete.
```



## Non-geographic map tiles

Almost all map tiles you encounter are for geographic maps. They have a geographic coordinate reference system (CRS). Software used to display these map tiles, such as Leaflet, is similarly focused on these kinds of map tiles. Everything is geared towards the dominant use case involving geospatial coordinate systems.

However, there are edge cases where non-geographic maps are required. These can be maps of outer space, game board maps, etc. The base map used to generate map tiles is usually a simple image like a `png`, `jpg` or `bmp` file. The coordinate reference system is a simple Cartesian coordinate system based on the matrix of pixels or grid cells that represent the image.

There is no longitude or latitude or more complex geospatial projection associated with these maps, which is why they are said to have a "simple CRS". Simple does not necessarily mean easier to work with, however, because geospatial tools, like Leaflet for example, do not cater naturally to non-geographic coordinate systems. Using these map tiles in Leaflet is possible, but takes a bit of non-standard effort.

### Basic example

One example was shown previously where a spatial map lacking critical spatial reference information was processed as a simple image. In the example below, this is the intent. Here, the map is a `png` file. It is a previously saved plot of the Albers-projected US map used in the earlier projected geotiff example. You can see it has a color key legend. As a simple image, all of this will be tiled.


```r
map <- system.file("maps/map.png", package = "tiler")
plotRGB(brick(map))
```

![plot of chunk ex7](vignette-ex7-1.png)

```r

tile(map, tile_dir, "0-3")
#> [1] "C:\\Users\\Matt\\AppData\\Local\\Temp\\RtmpQLDRP8"
#> Creating tiles. Please wait...
#> Creating tile viewer...
#> Complete.
```



The `tile` function will automatically process simple image files differently. There is no concept of projection, and coloring tiles is irrelevant because the image file has its own coloring already. Map tiles generated from regular image files can be used with `leaflet` if done properly. The generated tiles have a simple CRS that is based on the pixel dimensions of the image file. If you were to use these tiles in `leaflet` for example and you wanted to overlay map markers, you would have to first georeference your locations of interest based on the matrix rows and columns of the image. This is outside the scope of `tiler`. See the Leaflet JS and `leaflet` package documentation for details on using simple CRS/non-geographic tiles.

Using a `png` file is recommended for quality and file size. `jpg` may yield a lower quality result, while a large, high resolution `bmp` file may have an enormous file size compared to `png`.

`jpg` and `bmp` are *optionally* supported by `tiler`. This means they are not installed and imported with `tiler`. It is assumed the user will provide `png` images. If using `jpg` or `bmp` and the packages `jpeg` or `bmp` are not installed, respectively, `tile` will print a message to the console notifying of the required package installations.

## Additional arguments

Other arguments to `tile` include `format`, `resume`, `viewer` and `georef`.

`format` is either `xyz` (default) or `tms`. `gdal2tiles` generates TMS tiles, but XYZ may be more familiar. Tile format only applies to geographic maps. All simple image-based tiles are XYZ format. If setting `format = "tms"` you may need to do something similar in your Leaflet JavaScript or `leaflet` package R code for tiles to display with the proper orientation.

`resume = TRUE` simply avoids overwriting tiles by picking up where you left off without changing your set zoom levels and output path.

`viewer = TRUE` creates `preview.html` adjacent to the `tiles` directory for previewing tiles in the browser using Leaflet. `georef = TRUE` adds mouse click feedback to the Leaflet widget. Map markers with matrix index coordinate labels appear on mouse click to assist with georeferencing. `georef` only applies to non-geographic tiles. This allows for interactive georeferencing of pixels in the image.

Finally, `...` can pass along some additional arguments. See help documentation for details.

## Serving map tiles

Map tiles must be served online to be of much use. Serving map tiles is not the purpose of `tiler` but using your GitHub account is an easy way to do this. Create a GitHub repository, enable GitHub pages for the repository in the repository settings. If the repository is exclusively for serving your map tiles, just set the master branch as the source for your GitHub pages. After committing and pushing your tiles to GitHub, you can access them using a URL of the form

`https://<your account name>.github.io/maptiles/tiles/{z}/{x}/{y}.png`

if you store your tiles in a folder named `tiles` in a repository named `maptiles` for example.

Here are some [examples of non-geographic tile sets](https://github.com/leonawicz/tiles/) hosted on GitHub using Star Trek galaxy maps generated with `tiler` and here they are used in [Leaflet maps](https://leonawicz.github.io/rtrek/articles/sc.html).

### Leaflet examples using remotely hosted tiles

The following two examples use Leaflet to display interactive maps. The maps use remotely hosted map tiles that were created with `tiler`.

Geographic provider tiles based on the low resolution example map included in `tiler` of the 48 contiguous US states. These map tiles were originally generated with the following code


```r
file <- system.file("maps/map_wgs84.tif", package = "tiler")
tile(file, "tiles", "0-7")
```

and are hosted [here](https://github.com/leonawicz/tiles/tree/master/us48lr).


```r
library(leaflet)
tiles <- "https://leonawicz.github.io/tiles/us48lr/tiles/{z}/{x}/{y}.png"
leaflet(
  options = leafletOptions(minZoom = 0, maxZoom = 7), width = "100%") %>%
  addProviderTiles("Stamen.Toner") %>%
  addTiles(tiles, options = tileOptions(opacity = 0.8)) %>% setView(-100, 40, 3)
```

Non-geographic provider tiles based on a Star Trek fictional galaxy map. These map tiles were generated with


```r
tile("st2.png", "tiles", "0-7")
```

and are hosted [here](https://github.com/leonawicz/tiles/tree/master/st2), including the source file `st2.png`.


```r
tiles <- "https://leonawicz.github.io/tiles/st2/tiles/{z}/{x}/{y}.png"
leaflet(
  options = leafletOptions(
    crs = leafletCRS("L.CRS.Simple"), minZoom = 0, maxZoom = 7, attributionControl = FALSE), width = "100%") %>%
  addTiles(tiles) %>% setView(71, -60, 3)
```

*Note: These hosted tiles were made when `tiler` made XYZ-format tiles. It now strictly makes TMS tiles. For new tiles in `leaflet` switch the end of the url from `{y}` to `{-y}`.

### Local preview

By default `tile` also creates `preview.html` as noted above. This can also be created or re-created directly using `tile_viewer`. In either case, as long as `preview.html` exists alongside the tiles directory, it can easily be loaded in the browser with `view_tiles`. See help documentation for details.

If you have tiles in a directory `"project/tiles"`, creating `preview.html` directly can be done as follows. The arguments shown are just for illustration.


```r
tile_viewer("project/tiles", "3-7") # geographic tiles
tile_viewer("project/tiles", "3-7", width = 1000, height = 1000) # non-geographic tiles
```

However the tile preview document is created, it can be viewed by passing the same `tiles` directory to `view_tiles`.


```r
view_tiles("project/tiles")
```

### Details

The `leaflet` R code needed in order to use custom non-geographic tiles like these requires setting `leafletOptions(crs = leafletCRS("L.CRS.Simple"))` as well as calling `addTiles(urlTemplate = url)` where `url` is like the example URL shown above. Setting the focus of the map can be a bit tricky for non-geographic map tiles based on an arbitrary image file. It may take some trial and error to get a sense for the custom coordinate system.


---
title: "Introduction to tiler"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to tiler}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", message = FALSE, error = FALSE
)
set.seed(1)
tmpfiles <- list.files(tempdir(), full.names = TRUE) # any pre-existing temp files
```

The `tiler` package provides a map tile-generator function for creating map tile sets for use with packages such as `leaflet`.

In addition to generating map tiles based on a common raster layer source, it also handles the non-geographic edge case, producing map tiles from arbitrary images. These map tiles, which have a "simple CRS", a non-geographic simple Cartesian coordinate reference system, can also be used with `leaflet` when applying the simple CRS option. Map tiles can be created from an input file with any of the following extensions: `tif`, `grd` and `nc` for spatial maps and `png`, `jpg` and `bmp` for basic images.

## Setup (Windows users)

Python is required as well as the `gdal` library for Python. For Windows users, this is not the same as simply having `rgdal` installed through R. You need to install `gdal` so that it is accessible by Python. One of the easiest and typically recommended ways to do this in Windows is to install [OSGeo4W](https://trac.osgeo.org/osgeo4w/). It will bring all the required Python `gdal` library functionality along with it. OSGeo4W is also commonly installed along with [QGIS](https://qgis.org/en/site/forusers/download.html).

On Windows, `tiler_options` is set on package load with `osgeo4w = "OSgeo4W.bat"`. Make sure to add the path to this file to your system path variable. Otherwise it will not be found when calling `tile`. You can set the full path to `OSGeo4W.bat` directly in your R session with `tiler_options`. However, it is recommended to add the directory (e.g. `C:/OSGeo64` or wherever `OSGeo4W.bat` is located) to your system path so that you never had to deal with it in R.

Linux and Mac users should not have to do any additional setup as long as Python and `gdal` for Python are installed and available.

## Geographic map tiles

### Context

For the sake of simple, expedient examples, the map tiles generated below all use small zoom level ranges. There is also no reason to attempt displaying the tiles here. To make these examples more informative, each raster is loaded and plotted for context, though this is not necessary to the tiling process. Loading the `raster` package is only needed here for the print and plot calls.

The example maps packaged with `tiler` are not representative of large, high resolution imagery that benefits from tiling. These maps are very small in order to minimize package size and ensure examples run quickly. But the tiling procedures demonstrated are the same as would be applied to larger images.

Lastly, consider the power of your system before attempting to make a ton of tiles for large images at very high resolutions. You could find that the system could hang at any one of a number of choke points. If you are attempting to make thousands of tiles for a large, high resolution image and your system is struggling, it is recommended to (1) try making tiles for only one zoom level at a time, starting from zero and then increasing while monitoring your system resources. (2) If this is not enough, find a better system.

### Basic example

Map tiles are generated with `tile`. `tile` takes an input file name `file` for the source map and a `tiles` destination path for where to save map tiles. The only other required argument is `zoom`, the range of zoom levels for the tiles. This is a string of the form, e.g. `3:7`. In this and the subsequent examples zoom levels 0-3 are used.

```{r ex1}
library(tiler)
library(raster)
tile_dir <- file.path(tempdir(), "tiles")
map <- system.file("maps/map_wgs84.tif", package = "tiler")
(r <- raster(map))
plot(r)

tile(map, tile_dir, "0-3")

list.files(tile_dir)
```

```{r unlink1, echo=FALSE}
unlink(tile_dir, recursive = TRUE, force = TRUE)
```

Listing the files in `tile_dir` shows the top level map tiles directories, 0-3 as expected. This is not printed in subsequent examples since it is not going to change.

Note that these examples rendered to HTML here do not capture the parts of the log output that result from the internal system call made by `tile`. When you run this example yourself you will see a bit more information at the console.

### Projected maps

The previous example used a map with a geospatial coordinate reference system (CRS) but it was not projected. That map would be ready to view with the `leaflet` package for example, as would the tiles generated based on it. The next example uses a similar map that is projected. In order to generate map tiles that can be used with `leaflet` in the standard CRS, the map must be reprojected first. Then the same map tiles are generated as before.

```{r ex2}
map <- system.file("maps/map_albers.tif", package = "tiler")
(r <- raster(map))
plot(r)

tile(map, tile_dir, "0-3")
```

```{r unlink2, echo=FALSE}
unlink(tile_dir, recursive = TRUE, force = TRUE)
```

The tiles generated this time are the same as before, that is, ready for `leaflet`. `tile` reprojects the raster layer internally. This can be seen in the log output printed to the console.

### Missing CRS

If the CRS of the raster is `NA`, there are two options. By default, `tile` will fall back on processing the raster layer as if it were a simple image file with no geospatial projection information (see the next section on simple CRS/non-geographic map tiles). These tiles are not the same as the previous sets.

```{r ex3}
map <- system.file("maps/map_albers_NA.tif", package = "tiler")
(r <- raster(map))
plot(r)

tile(map, tile_dir, "0-3")
```

```{r unlink3, echo=FALSE}
unlink(tile_dir, recursive = TRUE, force = TRUE)
```

This is not likely what is wanted. The other option is to force set a known CRS if it is missing from the file or was dropped for whatever reason. Now reprojection can proceed and the expected tiles are generated.

```{r ex4}
crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"
tile(map, tile_dir, "0-3", crs)
```

```{r unlink4, echo=FALSE}
unlink(tile_dir, recursive = TRUE, force = TRUE)
```

A note on reprojection: Depending on the nature of the data in a raster, the `...` argument to `tile` allows you to pass through the `method` argument to `raster::projectRaster`. This is `bilinear` by default for bilinear interpolation, appropriate for continuous data. It can be set to `ngb` for nearest neighbor, appropriate for discrete or categorical data. If more control is needed over the reprojection, you should just prepare your raster file first before using `tile`. `tiler` is not intended to substitute for or wrap general geospatial processing tasks that can easily be done with other packages.

### Coloring tiles

Being able to change the default color palette or specify color breaks is important. All other optional `...` arguments to `tile` are passed to `raster::RGB` to provide control over the creation of an intermediary RGB raster. Most arguments to `RGB` can usually be ignored. The most useful ones are `col` and `colNA` for the data values color palette and the `noData` color, respectively. Coloring tiles differently for the original example is as simple as the following.

```{r ex5}
map <- system.file("maps/map_wgs84.tif", package = "tiler")
pal <- colorRampPalette(c("darkblue", "lightblue"))(20)
nodata <- "tomato"
tile(map, tile_dir, "0-3", col = pal, colNA = nodata)
```

```{r unlink5, echo=FALSE}
unlink(tile_dir, recursive = TRUE, force = TRUE)
```

### RGB and RGBA rasters

Multi-band rasters are supported as long as they have three or four layers, in which case `tile` assumes these represent red, green, blue and alpha channel, respectively. Internally, single-layer raster files are colored and converted to a three- or four-band RGB/A raster object prior to tile generation. If `file` is already such a raster, this step is simply skipped. Optional arguments like data and `noData` color, break points, etc., are ignored since this type of raster contains its own color information.

```{r ex6}
map <- system.file("maps/map_albers_rgb.tif", package = "tiler")
(r <- brick(map))
plot(r)

tile(map, tile_dir, "0-3")
```

```{r unlink6, echo=FALSE}
unlink(tile_dir, recursive = TRUE, force = TRUE)
```

## Non-geographic map tiles

Almost all map tiles you encounter are for geographic maps. They have a geographic coordinate reference system (CRS). Software used to display these map tiles, such as Leaflet, is similarly focused on these kinds of map tiles. Everything is geared towards the dominant use case involving geospatial coordinate systems. 

However, there are edge cases where non-geographic maps are required. These can be maps of outer space, game board maps, etc. The base map used to generate map tiles is usually a simple image like a `png`, `jpg` or `bmp` file. The coordinate reference system is a simple Cartesian coordinate system based on the matrix of pixels or grid cells that represent the image.

There is no longitude or latitude or more complex geospatial projection associated with these maps, which is why they are said to have a "simple CRS". Simple does not necessarily mean easier to work with, however, because geospatial tools, like Leaflet for example, do not cater naturally to non-geographic coordinate systems. Using these map tiles in Leaflet is possible, but takes a bit of non-standard effort.

### Basic example

One example was shown previously where a spatial map lacking critical spatial reference information was processed as a simple image. In the example below, this is the intent. Here, the map is a `png` file. It is a previously saved plot of the Albers-projected US map used in the earlier projected geotiff example. You can see it has a color key legend. As a simple image, all of this will be tiled.

```{r ex7}
map <- system.file("maps/map.png", package = "tiler")
plotRGB(brick(map))

tile(map, tile_dir, "0-3")
```

```{r unlink7, echo=FALSE}
unlink(tile_dir, recursive = TRUE, force = TRUE)
```

The `tile` function will automatically process simple image files differently. There is no concept of projection, and coloring tiles is irrelevant because the image file has its own coloring already. Map tiles generated from regular image files can be used with `leaflet` if done properly. The generated tiles have a simple CRS that is based on the pixel dimensions of the image file. If you were to use these tiles in `leaflet` for example and you wanted to overlay map markers, you would have to first georeference your locations of interest based on the matrix rows and columns of the image. This is outside the scope of `tiler`. See the Leaflet JS and `leaflet` package documentation for details on using simple CRS/non-geographic tiles.

Using a `png` file is recommended for quality and file size. `jpg` may yield a lower quality result, while a large, high resolution `bmp` file may have an enormous file size compared to `png`.

`jpg` and `bmp` are *optionally* supported by `tiler`. This means they are not installed and imported with `tiler`. It is assumed the user will provide `png` images. If using `jpg` or `bmp` and the packages `jpeg` or `bmp` are not installed, respectively, `tile` will print a message to the console notifying of the required package installations.

## Additional arguments

Other arguments to `tile` include `format`, `resume`, `viewer` and `georef`.

`format` is either `xyz` (default) or `tms`. `gdal2tiles` generates TMS tiles, but XYZ may be more familiar. Tile format only applies to geographic maps. All simple image-based tiles are XYZ format. If setting `format = "tms"` you may need to do something similar in your Leaflet JavaScript or `leaflet` package R code for tiles to display with the proper orientation.

`resume = TRUE` simply avoids overwriting tiles by picking up where you left off without changing your set zoom levels and output path.

`viewer = TRUE` creates `preview.html` adjacent to the `tiles` directory for previewing tiles in the browser using Leaflet. `georef = TRUE` adds mouse click feedback to the Leaflet widget. Map markers with matrix index coordinate labels appear on mouse click to assist with georeferencing. `georef` only applies to non-geographic tiles. This allows for interactive georeferencing of pixels in the image.

Finally, `...` can pass along some additional arguments. See help documentation for details.

## Serving map tiles

Map tiles must be served online to be of much use. Serving map tiles is not the purpose of `tiler` but using your GitHub account is an easy way to do this. Create a GitHub repository, enable GitHub pages for the repository in the repository settings. If the repository is exclusively for serving your map tiles, just set the master branch as the source for your GitHub pages. After committing and pushing your tiles to GitHub, you can access them using a URL of the form

`https://<your account name>.github.io/maptiles/tiles/{z}/{x}/{y}.png`

if you store your tiles in a folder named `tiles` in a repository named `maptiles` for example.

Here are some [examples of non-geographic tile sets](https://github.com/leonawicz/tiles/) hosted on GitHub using Star Trek galaxy maps generated with `tiler` and here they are used in [Leaflet maps](https://leonawicz.github.io/rtrek/articles/sc.html).

### Leaflet examples using remotely hosted tiles

The following two examples use Leaflet to display interactive maps. The maps use remotely hosted map tiles that were created with `tiler`.

Geographic provider tiles based on the low resolution example map included in `tiler` of the 48 contiguous US states. These map tiles were originally generated with the following code

```{r map1_tiles, eval=FALSE}
file <- system.file("maps/map_wgs84.tif", package = "tiler")
tile(file, "tiles", "0-7")
```

and are hosted [here](https://github.com/leonawicz/tiles/tree/master/us48lr).

```{r map1}
library(leaflet)
tiles <- "https://leonawicz.github.io/tiles/us48lr/tiles/{z}/{x}/{y}.png"
leaflet(
  options = leafletOptions(minZoom = 0, maxZoom = 7), width = "100%") %>% 
  addProviderTiles("Stamen.Toner") %>% 
  addTiles(tiles, options = tileOptions(opacity = 0.8)) %>% setView(-100, 40, 3)
```

Non-geographic provider tiles based on a Star Trek fictional galaxy map. These map tiles were generated with

```{r map2_tiles, eval=FALSE}
tile("st2.png", "tiles", "0-7")
``` 

and are hosted [here](https://github.com/leonawicz/tiles/tree/master/st2), including the source file `st2.png`.

```{r map2}
tiles <- "https://leonawicz.github.io/tiles/st2/tiles/{z}/{x}/{y}.png"
leaflet(
  options = leafletOptions(
    crs = leafletCRS("L.CRS.Simple"), minZoom = 0, maxZoom = 7, attributionControl = FALSE), width = "100%") %>% 
  addTiles(tiles) %>% setView(71, -60, 3)
```

*Note: These hosted tiles were made when `tiler` made XYZ-format tiles. It now strictly makes TMS tiles. For new tiles in `leaflet` switch the end of the url from `{y}` to `{-y}`.

### Local preview

By default `tile` also creates `preview.html` as noted above. This can also be created or re-created directly using `tile_viewer`. In either case, as long as `preview.html` exists alongside the tiles directory, it can easily be loaded in the browser with `view_tiles`. See help documentation for details.

If you have tiles in a directory `"project/tiles"`, creating `preview.html` directly can be done as follows. The arguments shown are just for illustration.

```{r ex8, eval=FALSE}
tile_viewer("project/tiles", "3-7") # geographic tiles
tile_viewer("project/tiles", "3-7", width = 1000, height = 1000) # non-geographic tiles
```

However the tile preview document is created, it can be viewed by passing the same `tiles` directory to `view_tiles`.

```{r ex9, eval=FALSE}
view_tiles("project/tiles")
```

### Details

The `leaflet` R code needed in order to use custom non-geographic tiles like these requires setting `leafletOptions(crs = leafletCRS("L.CRS.Simple"))` as well as calling `addTiles(urlTemplate = url)` where `url` is like the example URL shown above. Setting the focus of the map can be a bit tricky for non-geographic map tiles based on an arbitrary image file. It may take some trial and error to get a sense for the custom coordinate system.

```{r cleanup, echo=FALSE}
# supplemental check for excess temp files
extrafiles <- setdiff(list.files(tempdir(), full.names = TRUE),
                      c(tmpfiles, list.files(tempdir(), "\\.js$", full.names = TRUE)))
if(length(extrafiles)) unlink(extrafiles, recursive = TRUE, force = TRUE)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tiler.R
\docType{package}
\name{tiler}
\alias{tiler}
\title{tiler: Create map tiles from R}
\description{
The tiler package creates geographic map tiles from geospatial map files or
non-geographic map tiles from simple image files.
}
\details{
This package provides a tile generator function for creating map tile sets
for use with packages such as \code{leaflet}.
In addition to generating map tiles based on a common raster layer source,
it also handles the non-geographic edge case, producing map tiles from
arbitrary images.
These map tiles, which have a non-geographic simple coordinate reference
system (CRS), can also be used with \code{leaflet} when applying the simple
CRS option.
\cr\cr
Map tiles can be created from an input file with any of the following
extensions: tif, grd and nc for spatial maps and png, jpg and bmp for basic
images.
\cr\cr
This package requires Python and the \code{gdal} library for Python.
Windows users are recommended to install
\code{OSGeo4W}: \code{https://trac.osgeo.org/osgeo4w/} as an easy way to
obtain the required \code{gdal} support for Python in Windows.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/options.R
\name{tiler_options}
\alias{tiler_options}
\title{Options}
\usage{
tiler_options(...)
}
\arguments{
\item{...}{a list of options.}
}
\value{
The function prints all set options if called with no arguments.
When setting options, nothing is returned.
}
\description{
Options for tiler package.
}
\details{
On Windows systems, if the system paths for \code{python.exe} and
\code{OSGeo4W.bat} are not added to the system PATH variable, they must be
provided by the user after loading the package.
It is recommended to add these to the system path so they do not need to be
specified for every R session.

As long as you are using OSGeo4W, you can ignore the Python path
specification and do not even need to install it on your system separately;
OSGeo4W will use its own built-in version.

The recommended way to have GDAL available to Python in Windows is to
install \href{https://trac.osgeo.org/osgeo4w/}{OSGeo4W}. This is commonly
installed along with
\href{https://qgis.org/en/site/forusers/download.html}{QGIS}.

By default, \code{tiler_options} is set on package load with
\code{osgeo4w = "OSGeo4W.bat"}. It is expected that the user has added the
path to this file to the system PATH variable in Windows.
For example, if it is installed to \code{C:/OSGeo4W64/OSGeo4W.bat}, add
\code{C:/OSGeo4W64} to your PATH.
If you do want to specify the path in the R session using
\code{tiler_options}, provide the full path including the filename.
See the example.

None of this applies to other systems. As long as the system requirements,
Python and GDAL, are installed, then \code{tile} should generate tiles
without getting or setting any \code{tiler_options}.
}
\examples{
tiler_options()
tiler_options(osgeo4w = "C:/OSGeo4W64/OSGeo4W.bat")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viewer.R
\name{view_tiles}
\alias{view_tiles}
\title{View map tiles with Leaflet}
\usage{
view_tiles(tiles)
}
\arguments{
\item{tiles}{character, directory where tiles are stored.}
}
\value{
nothing is returned, but the default browser is launched.
}
\description{
View map tiles in the browser using leaflet.
}
\details{
This function opens \code{preview.html} in a web browser. This file displays
map tiles in a Leaflet widget.
The file is created when \code{tile} is called to generate the map tiles,
unless \code{viewer = FALSE}.
Alternatively, it is created (or re-created) subsequent to tile creation
using \code{tile_viewer}.
}
\examples{
# launches browser; requires an existing tile set
\dontrun{view_tiles(file.path(tempdir(), "tiles"))}
}
\seealso{
\code{\link{tile_viewer}}, \code{\link{tile}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tile.R
\name{tile}
\alias{tile}
\title{Create map tiles}
\usage{
tile(
  file,
  tiles,
  zoom,
  crs = NULL,
  resume = FALSE,
  viewer = TRUE,
  georef = TRUE,
  ...
)
}
\arguments{
\item{file}{character, input file.}

\item{tiles}{character, output directory for generated tiles.}

\item{zoom}{character, zoom levels. Example format: \code{"3-7"}.
See details.}

\item{crs}{character, Proj4 string. Use this to force set the CRS of a
loaded raster object from \code{file} in cases where the CRS is missing but
known, to avoid defaulting to non-geographic tiling.}

\item{resume}{logical, only generate missing tiles.}

\item{viewer}{logical, also create \code{preview.html} adjacent to
\code{tiles} directory for previewing tiles in the browser using Leaflet.}

\item{georef}{logical, for non-geographic tiles only. If
\code{viewer = TRUE}, then the Leaflet widget in \code{preview.html} will
add map markers with coordinate labels on mouse click to assist with
georeferencing of non-geographic tiles.}

\item{...}{additional arguments for projected maps: reprojection method or
any arguments to \code{raster::RGB}, e.g. \code{col} and \code{colNA}. See
details. Other additional arguments \code{lng} and \code{lat} can also be
passed to the tile previewer. See \code{\link{tile_viewer}} for details.}
}
\value{
nothing is returned but tiles are written to disk.
}
\description{
Create geographic and non-geographic map tiles from a file.
}
\details{
This function supports both geographic and non-geographic tile generation.
When \code{file} is a simple image file such as \code{png}, \code{tile}
generates non-geographic, simple CRS tiles.
Files that can be loaded by the \code{raster} package yield geographic
tiles as long as \code{file} has projection information.
If the raster object's proj4 string is \code{NA}, it falls back on
non-geographic tile generation and a warning is thrown.

Choice of appropriate zoom levels for non-geographic image files may depend
on the size of the image. A \code{zoom} value may be partially ignored for
image files under certain conditions.
For instance using the example \code{map.png} below, when passing strictly
\code{zoom = n} where \code{n} is less than 3, this still generates tiles
for zoom \code{n} up through 3.

\subsection{Supported file types}{
Supported simple CRS/non-geographic image file types include \code{png},
\code{jpg} and \code{bmp}.
For projected map data, supported file types include three types readable by
the \code{raster} package: \code{grd}, \code{tif}, and \code{nc} (requires
\code{ncdf4}).
Other currently unsupported file types passed to \code{file} throw an error.
}
\subsection{Raster file inputs}{
If a map file loadable by \code{raster} is a single-layer raster object,
tile coloring is applied.
To override default coloring of data and \code{noData} pixels, pass the
additional arguments \code{col} and \code{colNA} to \code{...}.
Multi-layer raster objects are rejected with an error message. The only
exception is a three- or four-band raster, which is assumed to represent
red, green, blue and alpha channels, respectively.
In this case, processing will continue but coloring arguments are ignored as
unnecessary.
\cr\cr
Prior to tiling, a geographically-projected raster layer is reprojected to
EPSG:4326 only if it has some other projection.
The only reprojection argument available through \code{...} is
\code{method}, which can be \code{"bilinear"} (default) or\code{"ngb"}.
If complete control over reprojection is required, this should be done prior
to passing the rasterized file to the \code{tile} function. Then no
reprojection is performed by \code{tile}.
When \code{file} consists of RGB or RGBA bands, \code{method} is ignored if
provided and reprojection uses nearest neighbor.
\cr\cr
It is recommended to avoid using a projected 4-band RGBA raster file.
However, the alpha channel appears to be ignored anyway.
\code{gdal2tiles} gives an internal warning.
Instead, create your RGBA raster file in unprojected form and it should
pass through to \code{gdal2tiles} without any issues.
Three-band RGB raster files are unaffected by reprojection.
}
\subsection{Tiles and Leaflet}{
\code{gdal2tiles} generates TMS tiles. If expecting XYZ, for example when
using with Leaflet, you can change the end of the URL to your hosted tiles
from \code{{z}/{x}/{y}.png} to \code{{z}/{x}/{-y}.png}.
\cr\cr
This function is supported by two different versions of \code{gdal2tiles}.
There is the standard version, which generates geospatial tiles in TMS
format. The other version of \code{gdal2tiles} handles basic image files
like a matrix of rows and columns, using a simple Cartesian coordinate
system based on pixel dimensions of the image file.
See the Leaflet JS library and \code{leaflet} package documentation for
working with custom tiles in Leaflet.
}
}
\examples{
\dontshow{tmpfiles <- list.files(tempdir(), full.names = TRUE)}
# non-geographic/simple CRS
x <- system.file("maps/map.png", package = "tiler")
tiles <- file.path(tempdir(), "tiles")
tile(x, tiles, "2-3")

# projected map
x <- system.file("maps/map_wgs84.tif", package = "tiler")
tile(x, tiles, 0)
\dontshow{
unlink(c(tiles, file.path(tempdir(), "preview.html")), recursive = TRUE,
       force = TRUE)
extrafiles <- setdiff(list.files(tempdir(), full.names = TRUE), tmpfiles)
if(length(extrafiles)) unlink(extrafiles, recursive = TRUE, force = TRUE)
}
}
\seealso{
\code{\link{view_tiles}}, \code{\link{tile_viewer}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viewer.R
\name{tile_viewer}
\alias{tile_viewer}
\title{Create an HTML tile preview}
\usage{
tile_viewer(tiles, zoom, width = NULL, height = NULL, georef = TRUE, ...)
}
\arguments{
\item{tiles}{character, directory where tiles are stored.}

\item{zoom}{character, zoom levels full range. Example format: \code{"3-7"}.}

\item{width}{\code{NULL} (default) for geospatial map tiles. The original
image width in pixels for non-geographic, simple CRS tiles.}

\item{height}{\code{NULL} (default) for geospatial map tiles. The original
image height in pixels for non-geographic, simple CRS tiles.}

\item{georef}{logical, for non-geographic tiles only.
If \code{viewer = TRUE}, then the Leaflet widget in \code{preview.html} will
add map markers with coordinate labels on mouse click to assist with
georeferencing of non-geographic tiles.}

\item{...}{additional optional arguments include \code{lng} and \code{lat}
for setting the view longitude and latitude. These three arguments only
apply to geographic tiles. Viewer centering is \code{0, 0} by default.}
}
\value{
nothing is returned, but a file is written to disk.
}
\description{
Create an HTML file that displays a tile preview using Leaflet.
}
\details{
This function creates a file \code{preview.html} adjacent to the
\code{tiles} base directory.
When loaded in the browser, this file displays map tiles from the adjacent
folder.
For example, if tiles are stored in \code{project/tiles}, this function
creates \code{project/preview.html}.

By default, \code{tile} creates this file. The only reasons to call
\code{tile_viewer} directly after producing map tiles are:
(1) if \code{viewer = FALSE} was set in the call to \code{tile},
(2) if \code{tile} was called multiple times, e.g., for different batches of
zoom levels, and thus the most recent call did not use the full zoom range,
or (3) \code{preview.html} was deleted for some other reason.

If calling this function directly, ensure that the min and max zoom, and
original image pixel dimensions if applicable, match the generated tiles.
These arguments are passed to \code{tile_viewer} automatically when called
within \code{tile}, based on the source file provided to \code{tile}.
}
\examples{
tile_viewer(file.path(tempdir(), "tiles"), "3-7") # requires existing tiles
}
\seealso{
\code{\link{view_tiles}}, \code{\link{tile}}
}
