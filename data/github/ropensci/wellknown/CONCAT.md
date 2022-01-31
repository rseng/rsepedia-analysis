wellknown
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/summary/wellknown)](https://cranchecks.info/pkgs/wellknown)
[![R-check](https://github.com/ropensci/wellknown/workflows/R-check/badge.svg)](https://github.com/ropensci/wellknown/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/wellknown/coverage.svg?branch=master)](https://codecov.io/github/ropensci/wellknown?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/wellknown)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/wellknown)](https://cran.r-project.org/package=wellknown)

`wellknown` - convert WKT to GeoJSON and vice versa, and other WKT utilities.

Inspiration partly comes from Python's geomet (https://github.com/geomet/geomet) - and the name from Javascript's wellknown (https://github.com/mapbox/wellknown) (it's a good name).

Docs: https://docs.ropensci.org/wellknown/


## Install

Stable version


```r
install.packages("wellknown")
```

Dev version


```r
pak::pkg_install("ropensci/wellknown")
# OR
install.packages("wellknown", repos="https://dev.ropensci.org")
```


```r
library("wellknown")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/wellknown/issues).
* License: MIT
* Get citation information for `wellknown` in R doing `citation(package = 'wellknown')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
wellknown 0.7.4
===============

### MINOR IMPROVEMENTS

* fix a few examples for `geojson2wkt()` that had calls to `matrix()` that were leading to warnings because the input vectors were not a multiple of the nrow/ncol (#35)


wellknown 0.7.2
===============

### MINOR IMPROVEMENTS

* fixes to make this package compatible with an upcoming version of BH (v1.75.0) (#33)


wellknown 0.7.0
===============

### NEW FEATURES

* Dropped import package V8 used in `wkt_wkb()` and `wkb_wkt()` functions; now using package `wk` for those functions. Javascript no longer used in the package; should make installation of this package easier on some platforms that had trouble installing V8  (#24) (#31)
* Gains new functions for working with WKT: `bounding_wkt`, `wkt_bounding`, `sf_convert`, `validate_wkt`, `wkt_centroid`, `wkt_coords`, `wkt_reverse`. As part of this, package now uses Rcpp and BH (boost headers), so installation from source requires compilation (#32)

### MINOR IMPROVEMENTS

* vignette available on docs site only now


wellknown 0.6.0
===============

### BUG FIXES

* fix to `wkt_wkb` method; support new version of V8 that converts JS buffers to raw vectors (#29)


wellknown 0.5.0
===============

### NEW FEATURES

* New functions `wkt_wkb()` and `wkb_wkt()` for converting WKT to WKB, and WKB to WKT. Depends on `V8` for doing the conversion. (#5)
* New function `get_centroid()` to get a centroid (lon, lat) for a WKT character string or a GeoJSON list object (#14) (#15)
* `wkt2geojson()` gains a new parameter `numeric`. It is `TRUE` by default, meaning that we convert values to numeric unless you set `numeric=FALSE` in which case we maintain numbers as strings, useful when you want to retain zero digits after the decimal (#14)
* `wkt2geojson()` gains a new parameter `simplify`, which if `TRUE` attempts to simplify from a multi- geometry type to a single type (e.g., mulitpolygon to polygon) when there's really only a single object. Applies to multi features only. (#20)
* Throughout package now we account for 3D and 4D WKT. For `wkt2geojson()` GeoJSON doesn't support a 4th dimension so we drop the 4th value, but for `geojson2wkt()` you can have GeoJSON with a 4th value so that you can convert it and any 3D data to WKT. We've added checks to make sure not more than 4D is used, and we follow `sf` by filling in zeros for any objects that are shorter in number of dimensions than the object with the largest number of dimensions (#18) (#23)
* `geojson2wkt()` inputs it accepts have changed. The function now accepts two different formats of GeoJSON like data. 1) The old format of full GeoJSON as a list like `list('type' = 'Point', 'coordinates' = c(116.4, 45.2))`, and 2) a simplified format `list(Point = c(116.4, 52.2))` (#17) (#19)

### MINOR IMPROVEMENTS

* Removed `magrittr` package. Simply load the package to have access to pipes (#25)
* Fixes to `lint()` function for validating WKT to make it work in more cases (#9)

### BUG FIXES

* Fixed bug in `wkt2geojson()` to not be case-sensitive to object names (e.g. , now `point`, `Point`, and `POINT` are all fine) (#16)


wellknown 0.1.0
===============

### NEW FEATURES

* Releasd to CRAN.
## Test environments

* local macOS, R 4.1.0
* ubuntu 16.04 (on GitHub Actions), R 4.1.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 5 reverse dependencies. No problems were
found

--------

This version makes the package compatible with the upcoming BH version and fixes a matrix warning from an example.

Thanks!
Scott Chamberlain
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->

# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/wellknown/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/wellknown.git`
* Make sure to track progress upstream (i.e., on our version of `wellknown` at `ropensci/wellknown`) by doing `git remote add upstream https://github.com/ropensci/wellknown.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/wellknown`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.3 (2020-10-10) |
|os       |macOS Catalina 10.15.7       |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |US/Pacific                   |
|date     |2020-10-20                   |

# Dependencies

|package   |old   |new   |Δ  |
|:---------|:-----|:-----|:--|
|wellknown |0.6.0 |0.7.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*wellknown
=========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  eval = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/summary/wellknown)](https://cranchecks.info/pkgs/wellknown)
[![R-check](https://github.com/ropensci/wellknown/workflows/R-check/badge.svg)](https://github.com/ropensci/wellknown/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/wellknown/coverage.svg?branch=master)](https://codecov.io/github/ropensci/wellknown?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/wellknown)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/wellknown)](https://cran.r-project.org/package=wellknown)

`wellknown` - convert WKT to GeoJSON and vice versa, and other WKT utilities.

Inspiration partly comes from Python's geomet (https://github.com/geomet/geomet) - and the name from Javascript's wellknown (https://github.com/mapbox/wellknown) (it's a good name).

Docs: https://docs.ropensci.org/wellknown/


## Install

Stable version

```{r}
install.packages("wellknown")
```

Dev version

```{r}
pak::pkg_install("ropensci/wellknown")
# OR
install.packages("wellknown", repos="https://dev.ropensci.org")
```

```{r}
library("wellknown")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/wellknown/issues).
* License: MIT
* Get citation information for `wellknown` in R doing `citation(package = 'wellknown')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Introduction to the wellknown package"
author: "Scott Chamberlain"
date: "2020-10-20"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{wellknown introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`wellknown` - convert WKT to GeoJSON and vice versa.

Inspiration partly comes from Python's [geomet/geomet](https://github.com/geomet/geomet) - and the name from Javascript's [wellknown](https://github.com/mapbox/wellknown) (it's a good name).

## Different interfaces

### WKT from R stuctures

There's a family of functions that make it easy to go from familiar R objects like lists and data.frames to WKT, including:

* `point()` - make a point, e.g., `POINT (-116 45)`
* `multipoint()` - make a multipoint, e.g., `MULTIPOINT ((100 3), (101 2))`
* `linestring()` - make a linestring, e.g., `LINESTRING (100 0, 101 1)`
* `polygon()` - make a polygon, e.g., `POLYGON ((100 0), (101 0), (101 1), (100 0))`
* `multipolygon()` - make a multipolygon, e.g., `MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))`

The above currently accept (depending on the fxn) `numeric`, `list`, and `data.frame` (and `character` for special case of `EMPTY` WKT objects).

### Geojson to WKT and vice versa

`geojson2wkt()` and `wkt2geojson()` cover a subset of the various formats available:

* `Point`
* `MultiPoint`
* `Polygon`
* `MultiPolygon`
* `LineString`
* `MultilineString`
* `Geometrycollection`

#### Geojson to WKT

`geojson2wkt()` converts any geojson as a list to a WKT string (the same format )

#### WKT to Geojson

`wkt2geojson()` converts any WKT string into geojson as a list. This list format for geojson can be used downstream e.g., in the `leaflet` package.

#### WKT to WKB, and vice versa

`wkt_wkb()` converts WKT to WKB, while `wkb_wkt()` converts WKB to WKT

## Install

Stable version


```r
install.packages("wellknown")
```

Dev version


```r
remotes::install_github("ropensci/wellknown")
# OR
install.packages("wellknown", repos="https://dev.ropensci.org")
```


```r
library("wellknown")
```

## GeoJSON to WKT

### Point


```r
point <- list(Point = c(116.4, 45.2, 11.1))
geojson2wkt(point)
#> [1] "POINT Z(116.4000000000000057  45.2000000000000028  11.0999999999999996)"
```

### Multipoint


```r
mp <- list(
  MultiPoint = matrix(c(100, 101, 3.14, 3.101, 2.1, 2.18), 
    ncol = 2)
)
geojson2wkt(mp)
#> [1] "MULTIPOINT ((100.0000000000000000 3.1010000000000000), (101.0000000000000000 2.1000000000000001), (3.1400000000000001 2.1800000000000002))"
```

### LineString


```r
st <- list(
  LineString = matrix(c(0.0, 2.0, 4.0, 5.0,
                         0.0, 1.0, 2.0, 4.0), ncol = 2)
)
geojson2wkt(st, fmt=0)
#> [1] "LINESTRING (0 0, 2 1, 4 2, 5 4)"
```

### Multilinestring


```r
multist <- list(
  MultiLineString = list(
   matrix(c(0, -2, -4, -1, -3, -5), ncol = 2),
   matrix(c(1.66, 10.9999, 10.9, 0, -31.5, 3.0, 1.1, 0), ncol = 2)
 )
)
geojson2wkt(multist)
#> [1] "MULTILINESTRING ((0.0000000000000000 -1.0000000000000000, -2.0000000000000000 -3.0000000000000000, -4.0000000000000000 -5.0000000000000000), (1.6599999999999999 -31.5000000000000000, 10.9999000000000002 3.0000000000000000, 10.9000000000000004 1.1000000000000001, 0.0000000000000000 0.0000000000000000))"
```

### Polygon


```r
poly <- list(
  Polygon = list(
    matrix(c(100.001, 101.1, 101.001, 100.001, 0.001, 0.001, 1.001, 0.001),
      ncol = 2),
    matrix(c(100.201, 100.801, 100.801, 100.201, 0.201, 0.201, 0.801, 0.201),
      ncol = 2)
  )
)
geojson2wkt(poly)
#> [1] "POLYGON ((100.0010000000000048 0.0010000000000000, 101.0999999999999943 0.0010000000000000, 101.0010000000000048 1.0009999999999999, 100.0010000000000048 0.0010000000000000), (100.2009999999999934 0.2010000000000000, 100.8010000000000019 0.2010000000000000, 100.8010000000000019 0.8010000000000000, 100.2009999999999934 0.2010000000000000))"
```

### Multipolygon


```r
mpoly <- list(
  MultiPolygon = list(
    list(
      matrix(c(100, 101, 101, 100, 0.001, 0.001, 1.001, 0.001), ncol = 2),
      matrix(c(100.2, 100.8, 100.8, 100.2, 0.2, 0.2, 0.8, 0.2), ncol = 2)
    ),
    list(
      matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 1.0), ncol = 3),
      matrix(c(9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 9.0), ncol = 3)
    )
  )
)
geojson2wkt(mpoly, fmt=1)
#> [1] "MULTIPOLYGON Z(((100.000 0.001 0.000, 101.000 0.001 0.000, 101.000 1.001 0.000, 100.000 0.001 0.000), (100.2 0.2 0.0, 100.8 0.2 0.0, 100.8 0.8 0.0, 100.2 0.2 0.0)), ((1.0 4.0 7.0, 2.0 5.0 8.0, 3.0 6.0 1.0), (9.0 12.0 3.0, 10.0 1.0 4.0, 11.0 2.0 9.0)))"
```

### GeometryCollection


```r
gmcoll <- list(
 GeometryCollection = list(
   list(type = 'Point', coordinates = c(0.0, 1.0)),
   list(type = 'LineString', coordinates = matrix(c(0.0, 2.0, 4.0, 5.0,
                           0.0, 1.0, 2.0, 4.0),
                           ncol = 2)),
   list(type = 'Polygon', coordinates = list(
     matrix(c(100.001, 101.1, 101.001, 100.001, 0.001, 0.001, 1.001, 0.001),
       ncol = 2),
     matrix(c(100.201, 100.801, 100.801, 100.201, 0.201, 0.201, 0.801, 0.201),
       ncol = 2)
  ))
 )
)
geojson2wkt(gmcoll, fmt=0)
#> [1] "GEOMETRYCOLLECTION (POINT (0 1), LINESTRING (0 0, 2 1, 4 2, 5 4), POLYGON ((100.001 0.001, 101.100 0.001, 101.001 1.001, 100.001 0.001), (100.201 0.201, 100.801 0.201, 100.801 0.801, 100.201 0.201)))"
```

### Convert json or character objects

You can convert directly from an object of class `json`, which is output from `jsonlite::toJSON()`.


```r
library("jsonlite")
(json <- toJSON(list(Point = c(-105, 39)), auto_unbox = TRUE))
#> {"Point":[-105,39]}
```


```r
geojson2wkt(json)
#> [1] "POINT (-105   39)"
```

And you can convert from a geojson character string:


```r
str <- '{"type":"LineString","coordinates":[[0,0,10],[2,1,20],[4,2,30],[5,4,40]]}'
geojson2wkt(str)
#> [1] "LINESTRING Z(0 0 10, 2 1 20, 4 2 30, 5 4 40)"
```

## WKT to GeoJSON

### Point

As a `Feature`


```r
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
wkt2geojson(str)
#> $type
#> [1] "Feature"
#> 
#> $geometry
#> $geometry$type
#> [1] "Point"
#> 
#> $geometry$coordinates
#> [1] -116.4   45.2
#> 
...
```

Not `Feature`


```r
wkt2geojson(str, feature=FALSE)
#> $type
#> [1] "Point"
#> 
#> $coordinates
#> [1] -116.4   45.2
#> 
#> attr(,"class")
#> [1] "geojson"
```

### Multipoint


```r
str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180))'
wkt2geojson(str, feature=FALSE)
#> $type
#> [1] "MultiPoint"
#> 
#> $coordinates
#>        [,1]  [,2]
#> [1,] 100.00 3.101
#> [2,] 101.00 2.100
#> [3,]   3.14 2.180
#> 
#> attr(,"class")
...
```

### Polygon


```r
str <- "POLYGON ((100 0.1, 101.1 0.3, 101 0.5, 100 0.1), (103.2 0.2, 104.8 0.2, 100.8 0.8, 103.2 0.2))"
wkt2geojson(str, feature=FALSE)
#> $type
#> [1] "Polygon"
#> 
#> $coordinates
#> $coordinates[[1]]
#>       [,1] [,2]
#> [1,] 100.0  0.1
#> [2,] 101.1  0.3
#> [3,] 101.0  0.5
#> [4,] 100.0  0.1
...
```

### MultiPolygon


```r
str <- "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)),
    ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35), (30 20, 20 25, 20 15, 30 20)))"
wkt2geojson(str, feature=FALSE)
#> $type
#> [1] "MultiPolygon"
#> 
#> $coordinates
#> $coordinates[[1]]
#> $coordinates[[1]][[1]]
#>      [,1] [,2]
#> [1,]   40   40
#> [2,]   20   45
#> [3,]   45   30
...
```

### Linestring


```r
wkt2geojson("LINESTRING (0 -1, -2 -3, -4 5)", feature=FALSE)
#> $type
#> [1] "LineString"
#> 
#> $coordinates
#>      [,1] [,2]
#> [1,]    0   -1
#> [2,]   -2   -3
#> [3,]   -4    5
#> 
#> attr(,"class")
...
```

## lint WKT


```r
lint("POINT (1 2)")
#> [1] TRUE
lint("LINESTRING EMPTY")
#> [1] TRUE
lint("MULTIPOINT ((1 2), (3 4), (-10 100))")
#> [1] TRUE
lint("POLYGON((20.3 28.6, 20.3 19.6, 8.5 19.6, 8.5 28.6, 20.3 28.6))")
#> [1] TRUE
lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))")
#> [1] TRUE
lint("POINT (1 2 3 4 5)")
#> [1] FALSE
lint("LINESTRING (100)")
#> [1] FALSE
lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, a b, 10 20, 5 10, 15 5)))")
#> [1] FALSE
```

## WKT <--> WKB

WKT to WKB


```r
## point
wkt_wkb("POINT (-116.4 45.2)")
#>  [1] 01 01 00 00 00 9a 99 99 99 99 19 5d c0 9a 99 99 99 99 99 46 40

## polygon
wkt_wkb("POLYGON ((100.0 0.0, 101.1 0.0, 101.0 1.0, 100.0 0.0))")
#>  [1] 01 03 00 00 00 01 00 00 00 04 00 00 00 00 00 00 00 00 00 59 40 00 00 00 00
#> [26] 00 00 00 00 66 66 66 66 66 46 59 40 00 00 00 00 00 00 00 00 00 00 00 00 00
#> [51] 40 59 40 00 00 00 00 00 00 f0 3f 00 00 00 00 00 00 59 40 00 00 00 00 00 00
#> [76] 00 00
```

WKB to WKT


```r
## point
(x <- wkt_wkb("POINT (-116.4 45.2)"))
#>  [1] 01 01 00 00 00 9a 99 99 99 99 19 5d c0 9a 99 99 99 99 99 46 40
wkb_wkt(x)
#> [1] "POINT (-116.4 45.2)"

## polygon
(x <- wkt_wkb("POLYGON ((100.0 0.0, 101.1 0.0, 101.0 1.0, 100.0 0.0))"))
#>  [1] 01 03 00 00 00 01 00 00 00 04 00 00 00 00 00 00 00 00 00 59 40 00 00 00 00
#> [26] 00 00 00 00 66 66 66 66 66 46 59 40 00 00 00 00 00 00 00 00 00 00 00 00 00
#> [51] 40 59 40 00 00 00 00 00 00 f0 3f 00 00 00 00 00 00 59 40 00 00 00 00 00 00
#> [76] 00 00
wkb_wkt(x)
#> [1] "POLYGON ((100 0, 101.1 0, 101 1, 100 0))"
```


## Bounding boxes

A bounding box is a very simple concept: a representation of the smallest area in which all the points in a dataset lie. In WKT, bounding boxes look like:

```
POLYGON((10 14,10 16,12 16,12 14,10 14))
```

Sometimes you've got WKT data like this - a Polygon, a LineString, whatever - and you want a bounding box in a format R can understand. The answer is `wkt_bounding`, which takes a vector of valid WKT objects and produces a data.frame or matrix of R representations, whichever you'd prefer:


```r
wkt <- c("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))",
         "LINESTRING (30 10, 10 90, 40 40)")
wkt_bounding(wkt)
#>   min_x min_y max_x max_y
#> 1    10    10    40    40
#> 2    10    10    40    90
```

Alternately you might want to go in the other direction and turn R bounding boxes into WKT objects. You can do that with, appropriately, `bounding_wkt`:


```r
bounding_wkt(min_x = 10, min_y = 10, max_x = 40, max_y = 40)
#> [1] "POLYGON((10 10,10 40,40 40,40 10,10 10))"
```

This accepts either a series of vectors, one for each min or max value, or a list of length-4 vectors. Either way, it produces a nice WKT representation of the R data you give it.

## WKT validation

`validate_wkt` takes a vector of WKT objects and spits out a data.frame containing whether each object is valid, and any comments the parser has in the case that it isn't:


```r
wkt <- c("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))",
         "ARGHLEFLARFDFG",
         "LINESTRING (30 10, 10 90, 40 some string)")
validate_wkt(wkt)
#>   is_valid
#> 1    FALSE
#> 2    FALSE
#> 3    FALSE
#>                                                                                                                          comments
#> 1                                           The WKT object has a different orientation from the default. Use ?wkt_correct to fix.
#> 2                                                                          Object could not be recognised as a supported WKT type
#> 3 bad lexical cast: source type value could not be interpreted as target at 'some' in 'linestring (30 10, 10 90, 40 some string)'
```

With this you can check and clean your data before you rely on it and watch all your code fall down in a heap.

## Coordinate and centroid extraction

WKT POLYGONs are often used to store latitude and longitude coordinates - and you can use `wkt_coords` to get them:


```r
wkt_coords(("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))"))
#>   object  ring lng lat
#> 1      1 outer  30  10
#> 2      1 outer  40  40
#> 3      1 outer  20  40
#> 4      1 outer  10  20
#> 5      1 outer  30  10
```

The result of a `wkt_coords` call is a data.frame of four columns - `object`, identifying which of the input WKT objects the row refers to, `ring` referring to the layer in that object, and then `lat` and `lng`.

Extracting centroids is also useful, and can be performed with `wkt_centroid`. Again, it's entirely vectorised and produces a data.frame:


```r
wkt_centroid(("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))"))
#>        lng     lat
#> 1 25.45455 26.9697
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{wkt_centroid}
\alias{wkt_centroid}
\title{Extract Centroid}
\usage{
wkt_centroid(wkt)
}
\arguments{
\item{wkt}{a character vector of WKT objects, represented as strings}
}
\value{
a data.frame of two columns, \code{lat} and \code{lng},
with each row containing the centroid from the corresponding wkt
object. In the case that the object is NA (or cannot be decoded)
the resulting values will also be NA
}
\description{
\code{get_centroid} identifies the 2D centroid
in a WKT object (or vector of WKT objects). Note that it assumes
cartesian values.
}
\examples{
wkt_centroid("POLYGON((2 1.3,2.4 1.7))")
}
\seealso{
\code{\link[=wkt_coords]{wkt_coords()}} to extract all coordinates, and
\code{\link[=wkt_bounding]{wkt_bounding()}} to extract a bounding box.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygon.R
\name{polygon}
\alias{polygon}
\title{Make WKT polygon objects}
\usage{
polygon(..., fmt = 16, third = "z")
}
\arguments{
\item{...}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{third}{(character) Only applicable when there are three dimensions.
If \code{m}, assign a \code{M} value for a measurement, and if \code{z} assign a
\code{Z} value for three-dimenionsal system. Case is ignored. An \code{M} value
represents  a measurement, while a \code{Z} value usually represents altitude
(but can be something like depth in a water based location).}
}
\description{
Make WKT polygon objects
}
\details{
You can create nested polygons with \code{list} and
\code{data.frame} inputs, but not from \code{numeric} inputs. See examples.
}
\examples{
## empty polygon
polygon("empty")
# polygon("stuff")

# numeric
polygon(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001),
  c(100.001, 0.001), fmt=2)

# data.frame
## single polygon
df <- us_cities[2:5,c('long','lat')]
df <- rbind(df, df[1,])
wktview(polygon(df, fmt=2), zoom=4)
## nested polygons
df2 <- data.frame(long = c(-85.9, -85.9, -93, -93, -85.9),
                  lat = c(37.5, 35.3, 35.3, 37.5, 37.5))
wktview(polygon(df, df2, fmt=2), zoom=4)

# matrix
mat <- matrix(c(df$long, df$lat), ncol = 2)
polygon(mat)

# list
# single list - creates single polygon
ply <- list(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001),
  c(100.001, 0.001))
wktview(polygon(ply, fmt=2), zoom=7)
# nested list - creates nested polygon
vv <- polygon(list(c(35, 10), c(45, 45), c(15, 40), c(10, 20), c(35, 10)),
   list(c(20, 30), c(35, 35), c(30, 20), c(20, 30)), fmt=2)
wktview(vv, zoom=3)
# multiple lists nested within a list
zz <- polygon(list(list(c(35, 10), c(45, 45), c(15, 40), c(10, 20), c(35, 10)),
   list(c(20, 30), c(35, 35), c(30, 20), c(20, 30))), fmt=2)
wktview(zz, zoom=3)


## a 3rd point is included
### numeric
polygon(c(100, 0, 30), c(101, 0, 30), c(101, 1, 30),
  c(100, 0, 30), fmt = 2)
polygon(c(100, 0, 30), c(101, 0, 30), c(101, 1, 30),
  c(100, 0, 30), fmt = 2, third = "m")

### data.frame
df <- us_cities[2:5, c('long','lat')]
df <- rbind(df, df[1,])
df$altitude <- round(runif(NROW(df), 10, 50))
polygon(df, fmt=2)
polygon(df, fmt=2, third = "m")

### matrix
mat <- matrix(c(df$long, df$lat, df$altitude), ncol = 3)
polygon(mat, fmt=2)
polygon(mat, fmt=2, third = "m")

### list
ply <- list(c(100, 0, 1), c(101, 0, 1), c(101, 1, 1),
  c(100, 0, 1))
polygon(ply, fmt=2)
polygon(ply, fmt=2, third = "m")


## a 4th point is included
### numeric
polygon(c(100, 0, 30, 3.5), c(101, 0, 30, 3.5), c(101, 1, 30, 3.5),
  c(100, 0, 30, 3.5), fmt = 2)

### data.frame
df <- us_cities[2:5, c('long','lat')]
df <- rbind(df, df[1,])
df$altitude <- round(runif(NROW(df), 10, 50))
df$weight <- round(runif(NROW(df), 0, 1), 1)
polygon(df, fmt=2)

### matrix
mat <- matrix(unname(unlist(df)), ncol = 4)
polygon(mat, fmt=2)

### list
ply <- list(c(100, 0, 1, 40), c(101, 0, 1, 44), c(101, 1, 1, 45),
  c(100, 0, 1, 49))
polygon(ply, fmt=2)
}
\seealso{
Other R-objects: 
\code{\link{circularstring}()},
\code{\link{geometrycollection}()},
\code{\link{linestring}()},
\code{\link{multilinestring}()},
\code{\link{multipoint}()},
\code{\link{multipolygon}()},
\code{\link{point}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wellknown-package.R
\docType{data}
\name{us_cities}
\alias{us_cities}
\title{This is the same data set from the maps library, named differently}
\format{
A list with 6 components, namely "name", "country.etc", "pop",
"lat", "long", and "capital", containing the city name, the state
abbreviation, approximate population (as at January 2006), latitude,
longitude and capital status indication (0 for non-capital, 1 for
capital, 2 for state capital.
}
\description{
This database is of us cities of population greater than about 40,000.
Also included are state capitals of any population size.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_featurecollection.R
\name{as_featurecollection}
\alias{as_featurecollection}
\title{As featurecollection}
\usage{
as_featurecollection(x)
}
\arguments{
\item{x}{(list) GeoJSON as a list}
}
\description{
Helper function to make a FeatureCollection list object for use
in vizualizing, e.g., with \code{leaflet}
}
\examples{
str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180),
(31.140 6.180), (31.140 78.180))'
x <- wkt2geojson(str, fmt = 2)
as_featurecollection(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_convert.R
\name{sf_convert}
\alias{sf_convert}
\title{Convert spatial objects to WKT}
\usage{
sf_convert(x)
}
\arguments{
\item{x}{an sf or sfc object. one or more can be submitted}
}
\value{
a character vector of WKT objects
}
\description{
\code{sp_convert} turns objects from the \code{sp} package
(SpatialPolygons, SpatialPolygonDataFrames) or the \code{sf} package
(sf, sfc, POLYGON, MULTIPOLYGON) - into WKT POLYGONs or MULTIPOLYGONs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wkb.R
\name{wkb}
\alias{wkb}
\alias{wkt_wkb}
\alias{wkb_wkt}
\title{Convert WKT to WKB}
\usage{
wkt_wkb(x, ...)

wkb_wkt(x, ...)
}
\arguments{
\item{x}{For \code{wkt_wkb()}, a \code{character} string representing a WKT object;
for \code{wkb_wkt()}, an of class \code{raw} representing a WKB object}

\item{...}{arguments passed on to \code{\link[wk:deprecated]{wk::wkt_translate_wkb()}} or
\code{\link[wk:deprecated]{wk::wkb_translate_wkt()}}}
}
\value{
\code{wkt_wkb} returns an object of class \code{raw}, a WKB
reprsentation. \code{wkb_wkt} returns an object of class \code{character},
a WKT representation
}
\description{
Convert WKT to WKB
}
\examples{
# WKT to WKB
## point
wkt_wkb("POINT (-116.4 45.2)")

## linestring
wkt_wkb("LINESTRING (-116.4 45.2, -118.0 47.0)")

## multipoint
### only accepts the below format, not e.g., ((1 2), (3 4))
wkt_wkb("MULTIPOINT (100.000 3.101, 101.00 2.10, 3.14 2.18)")

## polygon
wkt_wkb("POLYGON ((100.0 0.0, 101.1 0.0, 101.0 1.0, 100.0 0.0))")

# WKB to WKT
## point
(x <- wkt_wkb("POINT (-116.4 45.2)"))
wkb_wkt(x)

## linestring
(x <- wkt_wkb("LINESTRING (-116.4 45.2, -118.0 47.0)"))
wkb_wkt(x)

## multipoint
(x <- wkt_wkb("MULTIPOINT (100.000 3.101, 101.00 2.10, 3.14 2.18)"))
wkb_wkt(x)

## polygon
(x <- wkt_wkb("POLYGON ((100.0 0.0, 101.1 0.0, 101.0 1.0, 100.0 0.0))"))
wkb_wkt(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wellknown-package.R
\docType{package}
\name{wellknown-package}
\alias{wellknown-package}
\alias{wellknown}
\title{wellknown}
\description{
WKT to GeoJSON and vice versa
}
\examples{
# GeoJSON to WKT
point <- list(Point = c(116.4, 45.2, 11.1))
geojson2wkt(point)

# WKT to GeoJSON
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
wkt2geojson(str)

## lint WKT
lint("POINT (1 2)")
lint("POINT (1 2 3 4 5)")

# WKT <--> WKB
wkt_wkb("POINT (-116.4 45.2)")
wkb_wkt(wkt_wkb("POINT (-116.4 45.2)"))
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multipoint.R
\name{multipoint}
\alias{multipoint}
\title{Make WKT multipoint objects}
\usage{
multipoint(..., fmt = 16, third = "z")
}
\arguments{
\item{...}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{third}{(character) Only applicable when there are three dimensions.
If \code{m}, assign a \code{M} value for a measurement, and if \code{z} assign a
\code{Z} value for three-dimenionsal system. Case is ignored. An \code{M} value
represents  a measurement, while a \code{Z} value usually represents altitude
(but can be something like depth in a water based location).}
}
\description{
Make WKT multipoint objects
}
\examples{
## empty multipoint
multipoint("empty")
# multipoint("stuff")

# numeric
multipoint(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180))

# data.frame
df <- us_cities[1:25, c('long', 'lat')]
multipoint(df)

# matrix
mat <- matrix(c(df$long, df$lat), ncol = 2)
multipoint(mat)

# list
multipoint(list(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180)))


## a 3rd point is included
multipoint(c(100, 3, 0), c(101, 2, 0), c(3, 2, 0), 
  third = "z", fmt = 1)
multipoint(c(100, 3, 0), c(101, 2, 0), c(3, 2, 0), 
  third = "m", fmt = 1)

df <- us_cities[1:25, c('long', 'lat')]
df$altitude <- round(runif(25, 100, 500))
multipoint(df, fmt = 2)
multipoint(df, fmt = 2, third = "m")

mat <- matrix(1:9, 3)
multipoint(mat)
multipoint(mat, third = "m")

x <- list(c(100, 3, 0), c(101, 2, 1), c(3, 2, 5))
multipoint(x)


## a 4th point is included
multipoint(
  c(100, 3, 0, 500), c(101, 2, 0, 505), c(3, 2, 0, 600), 
  fmt = 1)

df <- us_cities[1:25, c('long', 'lat')]
df$altitude <- round(runif(25, 100, 500))
df$weight <- round(runif(25, 1, 100))
multipoint(df, fmt = 2)

mat <- matrix(1:12, 3)
multipoint(mat)

x <- list(c(100, 3, 0, 300), c(101, 2, 1, 200), c(3, 2, 5, 100))
multipoint(x, fmt = 3)
}
\seealso{
Other R-objects: 
\code{\link{circularstring}()},
\code{\link{geometrycollection}()},
\code{\link{linestring}()},
\code{\link{multilinestring}()},
\code{\link{multipolygon}()},
\code{\link{point}()},
\code{\link{polygon}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{validate_wkt}
\alias{validate_wkt}
\title{Validate WKT objects}
\usage{
validate_wkt(x)
}
\arguments{
\item{x}{a character vector of WKT objects.}
}
\value{
a data.frame of two columns, \code{is_valid} (containing
\code{TRUE} or \code{FALSE} values for whether the WKT object is parseable and
valid) and \code{comments} (containing any error messages
in the case that the WKT object is not). If the objects are simply NA,
both fields will contain NA.
}
\description{
\code{validate_wkt} takes a vector of WKT objects and validates
them, returning a data.frame containing the status of each entry and
(in the case it cannot be parsed) any comments as to what, in particular,
may be wrong with it. It does not, unfortunately, check whether the
object meets the WKT spec - merely that it is formatted correctly.
}
\examples{
wkt <- c("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))",
 "ARGHLEFLARFDFG",
 "LINESTRING (30 10, 10 90, 40 some string)")
validate_wkt(wkt)
}
\seealso{
\code{\link[=wkt_correct]{wkt_correct()}} for correcting WKT objects
that fail validity checks due to having a non-default orientation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{wkt_coords}
\alias{wkt_coords}
\title{Extract Latitude and Longitude from WKT polygons}
\usage{
wkt_coords(wkt)
}
\arguments{
\item{wkt}{a character vector of WKT objects}
}
\value{
a data.frame of four columns; \code{object} (containing which object
the row refers to), \code{ring} containing which layer of the object the row
refers to, \code{lng} and \code{lat}.
}
\description{
\code{wkt_coords} extracts lat/long values from WKT polygons,
specifically the outer shell of those polygons (working on the assumption
that said outer edge is what you want).

Because it assumes \strong{coordinates}, it also assumes a sphere - say, the
earth - and uses spherical coordinate values.
}
\examples{
wkt_coords("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))")
}
\seealso{
\code{\link[=wkt_bounding]{wkt_bounding()}} to extract a bounding box, and \code{\link[=wkt_centroid]{wkt_centroid()}}
to extract the centroid.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multilinestring.R
\name{multilinestring}
\alias{multilinestring}
\title{Make WKT multilinestring objects}
\usage{
multilinestring(..., fmt = 16, third = "z")
}
\arguments{
\item{...}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{third}{(character) Only applicable when there are three dimensions.
If \code{m}, assign a \code{M} value for a measurement, and if \code{z} assign a
\code{Z} value for three-dimenionsal system. Case is ignored. An \code{M} value
represents  a measurement, while a \code{Z} value usually represents altitude
(but can be something like depth in a water based location).}
}
\description{
Make WKT multilinestring objects
}
\details{
There is no \code{numeric} input option for multilinestring.
There is no way as of yet to make a nested multilinestring with
\code{data.frame} input, but you can do so with list input. See examples.
}
\examples{
## empty multilinestring
multilinestring("empty")
# multilinestring("stuff")

# character string
x <- "MULTILINESTRING ((30 20, 45 40, 10 40), (15 5, 40 10, 10 20))"
multilinestring(x)

# data.frame
df <- data.frame(long = c(30, 45, 10), lat = c(20, 40, 40))
df2 <- data.frame(long = c(15, 40, 10), lat = c(5, 10, 20))
multilinestring(df, df2, fmt=0)
lint(multilinestring(df, df2, fmt=0))
wktview(multilinestring(df, df2), zoom=3)

# matrix
mat <- matrix(c(df$long, df$lat), ncol = 2)
mat2 <- matrix(c(df2$long, df2$lat), ncol = 2)
multilinestring(mat, mat2, fmt=0)

# list
x1 <- list(c(30, 20), c(45, 40), c(10, 40))
x2 <- list(c(15, 5), c(40, 10), c(10, 20))
multilinestring(x1, x2, fmt=2)

polys <- list(
  list(c(30, 20), c(45, 40), c(10, 40)),
  list(c(15, 5), c(40, 10), c(10, 20))
)
wktview(multilinestring(polys, fmt=2), zoom=3)

# 3D
## data.frame
df <- data.frame(long = c(30, 45, 10), lat = c(20, 40, 40), altitude = 1:3)
df2 <- data.frame(long = c(15, 40, 10), lat = c(5, 10, 20), altitude = 1:3)
multilinestring(df, df2, fmt=0)
multilinestring(df, df2, fmt=0, third = "m")
## matrix
mat <- matrix(unname(unlist(df)), ncol = 3)
mat2 <- matrix(unname(unlist(df2)), ncol = 3)
multilinestring(mat, mat2, fmt=0)
multilinestring(mat, mat2, fmt=0, third = "m")
## list
x1 <- list(c(30, 20, 1), c(45, 40, 1), c(10, 40, 1))
x2 <- list(c(15, 5, 0), c(40, 10, 3), c(10, 20, 4))
multilinestring(x1, x2, fmt=2)
multilinestring(x1, x2, fmt=2, third = "m")


# 4D
## data.frame
df <- data.frame(long = c(30, 45, 10), lat = c(20, 40, 40), 
  altitude = 1:3, weight = 4:6)
df2 <- data.frame(long = c(15, 40, 10), lat = c(5, 10, 20), 
  altitude = 1:3, weight = 4:6)
multilinestring(df, df2, fmt=0)
## matrix
mat <- matrix(unname(unlist(df)), ncol = 4)
mat2 <- matrix(unname(unlist(df2)), ncol = 4)
multilinestring(mat, mat2, fmt=0)
## list
x1 <- list(c(30, 20, 1, 40), c(45, 40, 1, 40), c(10, 40, 1, 40))
x2 <- list(c(15, 5, 0, 40), c(40, 10, 3, 40), c(10, 20, 4, 40))
multilinestring(x1, x2, fmt=2)
}
\seealso{
Other R-objects: 
\code{\link{circularstring}()},
\code{\link{geometrycollection}()},
\code{\link{linestring}()},
\code{\link{multipoint}()},
\code{\link{multipolygon}()},
\code{\link{point}()},
\code{\link{polygon}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_centroid.R
\name{get_centroid}
\alias{get_centroid}
\title{Get a centroid from WKT or geojson}
\usage{
get_centroid(x)
}
\arguments{
\item{x}{Input, a wkt character string or geojson class object}
}
\value{
A length 2 numeric vector, with longitude first, latitude second
}
\description{
Get a centroid from WKT or geojson
}
\examples{
# WKT
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
get_centroid(str)
str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180))'
get_centroid(str)
str <- "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)),
 ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35), (30 20, 20 25, 20 15, 30 20)))"
get_centroid(str)

# Geojson as geojson class
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
get_centroid(wkt2geojson(str))
str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180))'
get_centroid(wkt2geojson(str))
str <- "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)),
 ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35), (30 20, 20 25, 20 15, 30 20)))"
get_centroid(wkt2geojson(str))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bounding_wkt.R
\name{bounding_wkt}
\alias{bounding_wkt}
\title{Generate Bounding Boxes}
\usage{
bounding_wkt(min_x, min_y, max_x, max_y, values = NULL)
}
\arguments{
\item{min_x}{a numeric vector of the minimum value for \code{x} coordinates.}

\item{min_y}{a numeric vector of the minimum value for \code{y} coordinates.}

\item{max_x}{a numeric vector of the maximum value for \code{x} coordinates.}

\item{max_y}{a numeric vector of the maximum value for \code{y} coordinates.}

\item{values}{as an alternative to specifying the various values as vectors,
a list of length-4 numeric vectors containing min and max x and y values, or
just a single vector fitting that spec. NULL (meaning that the other
parameters will be expected) by default.}
}
\value{
a character vector of WKT POLYGON objects
}
\description{
\code{bounding_wkt} takes bounding boxes, in various formats,
and turns them into WKT POLYGONs.
}
\examples{
# With individual columns
bounding_wkt(10, 12, 14, 16)

# With a list
bounding_wkt(values = list(c(10, 12, 14, 16)))
}
\seealso{
\code{\link[=wkt_bounding]{wkt_bounding()}}, to turn WKT objects of various types into
a matrix or data.frame of bounding boxes.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.json.R
\name{as_json}
\alias{as_json}
\title{Convert geojson R list to JSON}
\usage{
as_json(x, pretty = TRUE, auto_unbox = TRUE, ...)
}
\arguments{
\item{x}{Output from \code{\link[=wkt2geojson]{wkt2geojson()}}}

\item{pretty}{(logical) Adds indentation whitespace to JSON output. Can
be \code{TRUE}/\code{FALSE} or a number specifying the number of spaces to indent.
See \code{\link[jsonlite:prettify]{jsonlite::prettify()}}. Default: \code{TRUE}. Having \code{TRUE} as default
makes it easy to copy paste to a text editor, etc.}

\item{auto_unbox}{(logical) Automatically unbox all atomic vectors of
length 1. Default: \code{TRUE}}

\item{...}{Further args passed on to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}
}
\description{
Convert geojson R list to JSON
}
\examples{
str <- "POLYGON ((100 0.1, 101.1 0.3, 101 0.5, 100 0.1),
   (103.2 0.2, 104.8 0.2, 100.8 0.8, 103.2 0.2))"
as_json(wkt2geojson(str))
as_json(wkt2geojson(str), FALSE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{wkt_reverse}
\alias{wkt_reverse}
\title{Reverses the points within a geometry.}
\usage{
wkt_reverse(x)
}
\arguments{
\item{x}{a character vector of WKT objects, represented as strings}
}
\value{
a string, same length as given
}
\description{
\code{wkt_reverse} reverses the points in any of
point, multipoint, linestring, multilinestring, polygon, or
multipolygon
}
\details{
segment, box, and ring types not supported
}
\examples{
wkt_reverse("POLYGON((42 -26,42 -13,52 -13,52 -26,42 -26))")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lint.R
\name{lint}
\alias{lint}
\title{Validate WKT strings}
\usage{
lint(str)
}
\arguments{
\item{str}{A WKT string}
}
\value{
A logical (\code{TRUE} or \code{FALSE})
}
\description{
Validate WKT strings
}
\details{
This function uses R regex - there's no error messages about
what is wrong in the WKT.
}
\examples{
lint("POINT (1 2)")
lint("POINT (1 2 3)")
lint("LINESTRING EMPTY")
lint("LINESTRING (100 0, 101 1)")
lint("MULTIPOINT EMPTY")
lint("MULTIPOINT ((1 2), (3 4))")
lint("MULTIPOINT ((1 2), (3 4), (-10 100))")
lint("POLYGON ((1 2, 3 4, 0 5, 1 2))")
lint("POLYGON((20.3 28.6, 20.3 19.6, 8.5 19.6, 8.5 28.6, 20.3 28.6))")
lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))")
lint("TRIANGLE ((0 0, 0 1, 1 1, 0 0))")
lint("TRIANGLE ((0.1 0.1, 0.1 1.1, 1.1 1.1, 0.1 0.1))")
lint("CIRCULARSTRING (1 5, 6 2, 7 3)")
lint("CIRCULARSTRING (1 5, 6 2, 7 3, 5 6, 4 3)")
lint('COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 1, -1 0), (-1 0, 2 0))')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/point.R
\name{point}
\alias{point}
\title{Make WKT point objects}
\usage{
point(..., fmt = 16, third = "z")
}
\arguments{
\item{...}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{third}{(character) Only applicable when there are three dimensions.
If \code{m}, assign a \code{M} value for a measurement, and if \code{z} assign a
\code{Z} value for three-dimenionsal system. Case is ignored. An \code{M} value
represents  a measurement, while a \code{Z} value usually represents altitude
(but can be something like depth in a water based location).}
}
\description{
Make WKT point objects
}
\details{
The \code{third} parameter is used only when there are sets of
three points, and you can toggle whether the object gets a \code{Z} or \code{M}.

When four points are included, the object automatically gets
assigned \code{ZM}
}
\examples{
## empty point
point("empty")
# point("stuff")

## single point
point(-116.4, 45.2)
point(0, 1)

## single point, from data.frame
df <- data.frame(lon=-116.4, lat=45.2)
point(df)

## many points, from a data.frame
ussmall <- us_cities[1:5, ]
df <- data.frame(long = ussmall$long, lat = ussmall$lat)
point(df)

## many points, from a matrix
mat <- matrix(c(df$long, df$lat), ncol = 2)
point(mat)

## single point, from a list
point(list(c(100.0, 3.101)))

## many points, from a list
point(list(c(100.0, 3.101), c(101.0, 2.1), c(3.14, 2.18)))

## when a 3rd point is included
point(1:3, third = "m")
point(1:3, third = "z")
point(list(1:3, 4:6), third = "m")
point(list(1:3, 4:6), third = "z")
point(matrix(1:9, ncol = 3), third = "m")
point(matrix(1:9, ncol = 3), third = "z")
point(data.frame(1, 2, 3), third = "m")
point(data.frame(1, 2, 3), third = "z")
point(data.frame(1:3, 4:6, 7:9), third = "m")

## when a 4th point is included
point(1:4)
point(list(1:4, 5:8))
point(matrix(1:12, ncol = 4))
point(data.frame(1, 2, 3, 4))
point(data.frame(1:3, 4:6, 7:9, 10:12))
}
\seealso{
Other R-objects: 
\code{\link{circularstring}()},
\code{\link{geometrycollection}()},
\code{\link{linestring}()},
\code{\link{multilinestring}()},
\code{\link{multipoint}()},
\code{\link{multipolygon}()},
\code{\link{polygon}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wktview.R
\name{wktview}
\alias{wktview}
\title{Visualize geojson from a character string or list}
\usage{
wktview(x, center = NULL, zoom = 5, fmt = 16)
}
\arguments{
\item{x}{Input, a geojson character string or list}

\item{center}{(numeric) A length two vector of the form:
\verb{longitude, latitude}}

\item{zoom}{(integer) A number between 1 and 18 (1 zoomed out, 18 zoomed in)}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}
}
\value{
Opens a map with the geojson object(s) using \code{leaflet}
}
\description{
Visualize geojson from a character string or list
}
\examples{
\dontrun{
# point
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
wktview(str)

# multipoint
df <- us_cities[1:5,c('long','lat')]
str <- multipoint(df)
wktview(str, center = c(-100,40))
wktview(str, center = c(-100,40), zoom = 3)

# linestring
wktview(linestring(c(100.000, 0.000), c(101.000, 1.000), fmt=2),
  center = c(100, 0))

# polygon
a <- polygon(list(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001),
  c(100.001, 0.001)))
wktview(a, center = c(100, 0))
wktview(a, center = c(100.5, 0), zoom=9)
}
}
\seealso{
\code{\link[=as_featurecollection]{as_featurecollection()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/properties.R
\name{properties}
\alias{properties}
\title{Add properties to a GeoJSON object}
\usage{
properties(x, style = NULL, popup = NULL)
}
\arguments{
\item{x}{(list) GeoJSON as a list}

\item{style}{(list) named list of color, fillColor, etc. attributes.
Default: \code{NULL}}

\item{popup}{(list) named list of popup values. Default: \code{NULL}}
}
\description{
Add properties to a GeoJSON object
}
\examples{
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
x <- wkt2geojson(str)
properties(x, style = list(color = "red"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linestring.R
\name{linestring}
\alias{linestring}
\title{Make WKT linestring objects}
\usage{
linestring(..., fmt = 16, third = "z")
}
\arguments{
\item{...}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{third}{(character) Only applicable when there are three dimensions.
If \code{m}, assign a \code{M} value for a measurement, and if \code{z} assign a
\code{Z} value for three-dimenionsal system. Case is ignored. An \code{M} value
represents  a measurement, while a \code{Z} value usually represents altitude
(but can be something like depth in a water based location).}
}
\description{
Make WKT linestring objects
}
\examples{
## empty linestring
linestring("empty")
# linestring("stuff")

## character string
linestring("LINESTRING (-116.4 45.2, -118.0 47.0)")

# numeric
## 2D
linestring(c(100.000, 0.000), c(101.000, 1.000), fmt=2)
linestring(c(100.0, 0.0), c(101.0, 1.0), c(120.0, 5.00), fmt=2)
## 3D
linestring(c(0.0, 0.0, 10.0), c(2.0, 1.0, 20.0),
           c(4.0, 2.0, 30.0), c(5.0, 4.0, 40.0), fmt=2)
## 4D
linestring(c(0.0, 0.0, 10.0, 5.0), c(2.0, 1.0, 20.0, 5.0),
           c(4.0, 2.0, 30.0, 5.0), c(5.0, 4.0, 40.0, 5.0), fmt=2)

# data.frame
df <- data.frame(lon=c(-116.4,-118), lat=c(45.2,47))
linestring(df, fmt=1)
df <- data.frame(lon=c(-116.4,-118,-120), lat=c(45.2,47,49))
linestring(df, fmt=1)
## 3D
df$altitude <- round(runif(NROW(df), 10, 50))
linestring(df, fmt=1)
linestring(df, fmt=1, third = "m")
## 4D
df$weight <- round(runif(NROW(df), 0, 1), 1)
linestring(df, fmt=1)


# matrix
mat <- matrix(c(-116.4,-118, 45.2, 47), ncol = 2)
linestring(mat, fmt=1)
mat2 <- matrix(c(-116.4, -118, -120, 45.2, 47, 49), ncol = 2)
linestring(mat2, fmt=1)
## 3D
mat <- matrix(c(df$long, df$lat, df$altitude), ncol = 3)
polygon(mat, fmt=2)
polygon(mat, fmt=2, third = "m")
## 4D
mat <- matrix(unname(unlist(df)), ncol = 4)
polygon(mat, fmt=2)

# list
linestring(list(c(100.000, 0.000), c(101.000, 1.000)), fmt=2)
## 3D
line <- list(c(100, 0, 1), c(101, 0, 1), c(101, 1, 1),
  c(100, 0, 1))
linestring(line, fmt=2)
linestring(line, fmt=2, third = "m")
## 4D
line <- list(c(100, 0, 1, 40), c(101, 0, 1, 44), c(101, 1, 1, 45),
  c(100, 0, 1, 49))
linestring(line, fmt=2)
}
\seealso{
Other R-objects: 
\code{\link{circularstring}()},
\code{\link{geometrycollection}()},
\code{\link{multilinestring}()},
\code{\link{multipoint}()},
\code{\link{multipolygon}()},
\code{\link{point}()},
\code{\link{polygon}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wkt2geojson.R
\name{wkt2geojson}
\alias{wkt2geojson}
\title{Convert WKT to GeoJSON-like objects.}
\usage{
wkt2geojson(str, fmt = 16, feature = TRUE, numeric = TRUE, simplify = FALSE)
}
\arguments{
\item{str}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{feature}{(logical) Make a feature geojson object. Default: \code{TRUE}}

\item{numeric}{(logical) Give back values as numeric. Default: \code{TRUE}}

\item{simplify}{(logical) Attempt to simplify from a multi- geometry type
to a single type. Applies to multi features only. Default: \code{FALSE}}
}
\description{
Convert WKT to GeoJSON-like objects.
}
\details{
Should be robust against a variety of typing errors, including
extra spaces between coordinates, no space between WKT type and coordinates.
However, some things won't pass, including lowercase WKT types, no
spaces between coordinates.

WKT with a 3rd value and when Z is found will be left as is and assumed to
be a altitude or similar value. WKT with a 3rd value and when M is found
will be discarded as the GeoJSON spec says to do so. WKT with a 4th value
as (presumably as a measurement) will also be discarded.
}
\examples{
# point
str <- "POINT (-116.4000000000000057 45.2000000000000028)"
wkt2geojson(str)
wkt2geojson(str, feature=FALSE)
wkt2geojson(str, numeric=FALSE)
wkt2geojson("POINT (-116 45)")
wkt2geojson("POINT (-116 45 0)")
## 3D
wkt2geojson("POINT Z(100 3 35)")
wkt2geojson("POINT M(100 3 35)") # dropped if M
## 4D
wkt2geojson("POINT ZM(100 3 35 1.5)") # Z retained

# multipoint
str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180))'
wkt2geojson(str, fmt = 2)
wkt2geojson(str, fmt = 2, feature=FALSE)
wkt2geojson(str, numeric=FALSE)
wkt2geojson("MULTIPOINT ((100 3), (101 2), (3 2))")
wkt2geojson("MULTIPOINT ((100 3 0), (101 2 0), (3 2 0))")
wkt2geojson("MULTIPOINT ((100 3 0 1), (101 2 0 1), (3 2 0 1))") 
## 3D
wkt2geojson("MULTIPOINT Z((100 3 35), (101 2 45), (3 2 89))")
wkt2geojson("MULTIPOINT M((100 3 1.3), (101 2 1.4), (3 2 1.9))")
## 4D
wkt2geojson("MULTIPOINT ZM((100 3 35 0), (101 2 45 0), (3 2 89 0))")

## simplify
wkt2geojson("MULTIPOINT ((100 3))", simplify = FALSE)
wkt2geojson("MULTIPOINT ((100 3))", simplify = TRUE)


# polygon
str <- "POLYGON ((100 0.1, 101.1 0.3, 101 0.5, 100 0.1),
   (103.2 0.2, 104.8 0.2, 100.8 0.8, 103.2 0.2))"
wkt2geojson(str)
wkt2geojson(str, feature=FALSE)
wkt2geojson(str, numeric=FALSE)
## 3D
str <- "POLYGON Z((100 0.1 3, 101.1 0.3 1, 101 0.5 5, 100 0.1 8),
   (103.2 0.2 3, 104.8 0.2 4, 100.8 0.8 5, 103.2 0.2 9))"
wkt2geojson(str)
## 4D
str <- "POLYGON ZM((100 0.1 3 0, 101.1 0.3 1 0, 101 0.5 5 0, 100 0.1 8 0),
   (103.2 0.2 3 0, 104.8 0.2 4 0, 100.8 0.8 5 0, 103.2 0.2 9 0))"
wkt2geojson(str)


# multipolygon
str <- "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)),
 ((20 35, 45 20, 30 5, 10 10, 10 30, 20 35), (30 20, 20 25, 20 15, 30 20)))"
wkt2geojson(str)
wkt2geojson(str, feature=FALSE)
wkt2geojson(str, numeric=FALSE)
## 3D
str <- "MULTIPOLYGON Z(((40 40 1, 20 45 3, 45 30 10, 40 40 0)),
 ((20 35 5, 45 20 67, 30 5 890, 10 10 2, 10 30 0, 20 35 4), 
 (30 20 4, 20 25 54, 20 15 56, 30 20 89)))"
wkt2geojson(str)
## 4D
str <- "MULTIPOLYGON ZM(((40 40 1 0, 20 45 3 4, 45 30 10 45, 40 40 0 1)),
 ((20 35 5 8, 45 20 67 9, 30 5 890 89, 10 10 2 234, 10 30 0 5, 20 35 4 1), 
 (30 20 4 0, 20 25 54 5, 20 15 56 55, 30 20 89 78)))"
wkt2geojson(str)

# simplify multipolygon to polygon if possible
str <- "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)))"
wkt2geojson(str, simplify = FALSE)
wkt2geojson(str, simplify = TRUE)


# linestring
str <- "LINESTRING (100.000 0.000, 101.000 1.000)"
wkt2geojson(str)
wkt2geojson(str, feature = FALSE)
wkt2geojson("LINESTRING (0 -1, -2 -3, -4 5)")
wkt2geojson("LINESTRING (0 1 2, 4 5 6)")
wkt2geojson(str, numeric = FALSE)
## 3D
wkt2geojson("LINESTRING Z(100.000 0.000 3, 101.000 1.000 5)")
wkt2geojson("LINESTRING M(100.000 0.000 10, 101.000 1.000 67)")
## 4D
wkt2geojson("LINESTRING ZM(100 0 1 4, 101 1 5 78)")


# multilinestring
str <- "MULTILINESTRING ((30 1, 40 30, 50 20)(10 0, 20 1))"
wkt2geojson(str)
wkt2geojson(str, numeric=FALSE)

str <- "MULTILINESTRING (
   (-105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5,
     -105.0 39.5, -105.0 39.5),
   (-105.0 39.5, -105.0 39.5, -105.0 39.5),
   (-105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5),
   (-105.0 39.5, -105.0 39.5, -105.0 39.5, -105.0 39.5))"
wkt2geojson(str)
wkt2geojson(str, numeric=FALSE)

## 3D
str <- "MULTILINESTRING Z((30 1 0, 40 30 0, 50 20 0)(10 0 1, 20 1 1))"
wkt2geojson(str)
str <- "MULTILINESTRING M((30 1 0, 40 30 0, 50 20 0)(10 0 1, 20 1 1))"
wkt2geojson(str)
## 4D
str <- "MULTILINESTRING ZM((30 1 0 5, 40 30 0 7, 50 20 0 1)(10 0 1 1, 20 1 1 1))"
wkt2geojson(str)

# simplify multilinestring to linestring if possible
str <- "MULTILINESTRING ((30 1, 40 30, 50 20))"
wkt2geojson(str, simplify = FALSE)
wkt2geojson(str, simplify = TRUE)


# Geometrycollection
str <- "GEOMETRYCOLLECTION (
 POINT Z(0 1 4),
 LINESTRING (-100 0, -101 -1),
 POLYGON ((100.001 0.001, 101.1235 0.0010, 101.001 1.001, 100.001 0.001),
           (100.201 0.201, 100.801 0.201, 100.801 0.801, 100.201 0.201)),
 MULTIPOINT ((100.000 3.101), (101.0 2.1), (3.14 2.18)),
 MULTILINESTRING ((0 -1, -2 -3, -4 -5),
       (1.66 -31.50, 10.0 3.0, 10.9 1.1, 0.0 4.4)),
 MULTIPOLYGON (((100.001 0.001, 101.001 0.001, 101.001 1.001, 100.001 0.001),
               (100.201 0.201, 100.801 0.201, 100.801 0.801, 100.201 0.201)),
                 ((1 2 3, 5 6 7, 9 10 11, 1 2 3))))"
wkt2geojson(str)
wkt2geojson(str, numeric=FALSE)

# case doesn't matter
str <- "point (-116.4000000000000057 45.2000000000000028)"
wkt2geojson(str)
}
\references{
\url{https://tools.ietf.org/html/rfc7946}
}
\seealso{
\code{\link[=geojson2wkt]{geojson2wkt()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circularstring.R
\name{circularstring}
\alias{circularstring}
\title{Make WKT circularstring objects}
\usage{
circularstring(..., fmt = 16)
}
\arguments{
\item{...}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}
}
\description{
Make WKT circularstring objects
}
\examples{
## empty circularstring
circularstring("empty")
# circularstring("stuff")

# Character string
circularstring("CIRCULARSTRING(1 5, 6 2, 7 3)")

# data.frame
df <- data.frame(lon = c(-116.4, -118), lat = c(45.2, 47))
circularstring(df, fmt=1)
df <- data.frame(lon=c(-116.4, -118, -120), lat=c(45.2, 47, 49))
circularstring(df, fmt=1)

# matrix
mat <- matrix(c(-116.4,-118, 45.2, 47), ncol = 2)
circularstring(mat, fmt=1)
mat2 <- matrix(c(-116.4, -118, -120, 45.2, 47, 49), ncol = 2)
circularstring(mat2, fmt=1)

# list
x <- list(c(1, 5), c(6, 2), c(7, 3))
circularstring(x, fmt=2)
}
\seealso{
Other R-objects: 
\code{\link{geometrycollection}()},
\code{\link{linestring}()},
\code{\link{multilinestring}()},
\code{\link{multipoint}()},
\code{\link{multipolygon}()},
\code{\link{point}()},
\code{\link{polygon}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{wkt_bounding}
\alias{wkt_bounding}
\title{Convert WKT Objects into Bounding Boxes}
\usage{
wkt_bounding(wkt, as_matrix = FALSE)
}
\arguments{
\item{wkt}{a character vector of WKT objects.}

\item{as_matrix}{whether to return the results as a matrix (\code{TRUE})
or data.frame (\code{FALSE}). Set to \code{FALSE} by default.}
}
\value{
either a data.frame or matrix, depending on the value of
\code{as_matrix}, containing four columns - \code{min_x}, \code{min_y}, \code{max_x} and
\code{max_y} - representing the various points of the bounding box. In the
event that a valid bounding box cannot be generated
(due to the invalidity or incompatibility of the WKT object), NAs will
be returned.
}
\description{
\code{wkt_bounding} turns WKT objects (specifically points,
linestrings, polygons, and multi-points/linestrings/polygons) into
bounding boxes.
}
\examples{
wkt_bounding("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))")
}
\seealso{
\code{\link[=bounding_wkt]{bounding_wkt()}}, to turn R-size bounding boxes into WKT objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrycollection.R
\name{geometrycollection}
\alias{geometrycollection}
\title{Make WKT geometrycollection objects}
\usage{
geometrycollection(...)
}
\arguments{
\item{...}{Character string WKT objects representing a Point, LineString,
Polygon, etc.}
}
\description{
Make WKT geometrycollection objects
}
\details{
This is different from the other functions that create WKT from R
objects, in that we can't do the same thing for GeometryCollection's since
many different WkT object could be created from the same input. So,
this function accepts WKT strings already formed and attempts to creat a
GeommetryCollection from them.
}
\examples{
## empty geometrycollection
geometrycollection("empty")
# geometrycollection("stuff")

# Character string, returns itself
geometrycollection("GEOMETRYCOLLECTION(POINT(4 6), LINESTRING(4 6, 7 10))")

# From a point
geometrycollection(point(-116.4, 45.2))

# From two points
geometrycollection(point(-116.4, 45.2), point(-118.4, 49.2))

# From various object types
geometrycollection(point(-116.4, 45.2),
 linestring("LINESTRING (-116.4 45.2, -118.0 47.0)"),
 circularstring(list(c(1, 5), c(6, 2), c(7, 3)), fmt = 2)
)
}
\seealso{
Other R-objects: 
\code{\link{circularstring}()},
\code{\link{linestring}()},
\code{\link{multilinestring}()},
\code{\link{multipoint}()},
\code{\link{multipolygon}()},
\code{\link{point}()},
\code{\link{polygon}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{wkt_correct}
\alias{wkt_correct}
\title{Correct Incorrectly Oriented WKT Objects}
\usage{
wkt_correct(x)
}
\arguments{
\item{x}{a character vector of WKT objects to correct}
}
\value{
a character vector, the same length as \code{x}, containing
either the original value (if there was no correction to make, or if
the object was invalid for other reasons) or the corrected WKT
value.
}
\description{
\code{wkt_correct} does precisely what it says on the tin,
correcting the orientation of WKT objects that are improperly oriented
(say, back to front). It can be applied to WKT objects that,
when validated with \code{\link[=validate_wkt]{validate_wkt()}}, fail for that reason.
}
\examples{
# A WKT object
wkt <- "POLYGON((30 20, 10 40, 45 40, 30 20), (15 5, 5 10, 10 20, 40 10, 15 5))"

# That's invalid due to a non-default orientation
validate_wkt(wkt)

# And suddenly isn't!
wkt_correct(wkt)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multipolygon.R
\name{multipolygon}
\alias{multipolygon}
\title{Make WKT multipolygon objects}
\usage{
multipolygon(..., fmt = 16, third = "z")
}
\arguments{
\item{...}{A GeoJSON-like object representing a Point, LineString, Polygon,
MultiPolygon, etc.}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{third}{(character) Only applicable when there are three dimensions.
If \code{m}, assign a \code{M} value for a measurement, and if \code{z} assign a
\code{Z} value for three-dimenionsal system. Case is ignored. An \code{M} value
represents  a measurement, while a \code{Z} value usually represents altitude
(but can be something like depth in a water based location).}
}
\description{
Make WKT multipolygon objects
}
\details{
There is no \code{numeric} input option for multipolygon. There
is no way as of yet to make a nested multipolygon with \code{data.frame}
input, but you can do so with list input. See examples.
}
\examples{
## empty multipolygon
multipolygon("empty")
# multipolygon("stuff")

# data.frame
df <- data.frame(long = c(30, 45, 10, 30), lat = c(20, 40, 40, 20))
df2 <- data.frame(long = c(15, 40, 10, 5, 15), lat = c(5, 10, 20, 10, 5))
multipolygon(df, df2, fmt=0)
lint(multipolygon(df, df2, fmt=0))
wktview(multipolygon(df, df2), zoom=3)

# matrix
mat <- matrix(c(df$long, df$lat), ncol = 2)
mat2 <- matrix(c(df2$long, df2$lat), ncol = 2)
multipolygon(mat, mat2, fmt=0)

# list
multipolygon(list(c(30, 20), c(45, 40), c(10, 40), c(30, 20)),
  list(c(15, 5), c(40, 10), c(10, 20), c(5, 10), c(15, 5)), fmt=2)

polys <- list(
  list(c(30, 20), c(45, 40), c(10, 40), c(30, 20)),
  list(c(15, 5), c(40, 10), c(10, 20), c(5, 10), c(15, 5))
)
wktview(multipolygon(polys, fmt=2), zoom=3)

## nested polygons
polys <- list(
  list(c(40, 40), c(20, 45), c(45, 30), c(40, 40)),
  list(
    list(c(20, 35), c(10, 30), c(10, 10), c(30, 5), c(45, 20), c(20, 35)),
    list(c(30, 20), c(20, 15), c(20, 25), c(30, 20))
  )
)
multipolygon(polys, fmt=0)
lint(multipolygon(polys, fmt=0))



# 3D
## data.frame
df <- data.frame(long = c(30, 45, 10, 30), lat = c(20, 40, 40, 20), 
  altitude = 1:4)
df2 <- data.frame(long = c(15, 40, 10, 5, 15), lat = c(5, 10, 20, 10, 5), 
  altitude = 1:5)
multipolygon(df, df2, fmt=0)
multipolygon(df, df2, fmt=0, third = "m")
## matrix
mat <- matrix(unname(unlist(df)), ncol = 3)
mat2 <- matrix(unname(unlist(df2)), ncol = 3)
multipolygon(mat, mat2, fmt=0)
multipolygon(mat, mat2, fmt=0, third = "m")
## list
l1 <- list(c(30, 20, 2), c(45, 40, 2), c(10, 40, 2), c(30, 20, 2))
l2 <- list(c(15, 5, 5), c(40, 10, 5), c(10, 20, 5), c(5, 10, 5), 
  c(15, 5, 5))
multipolygon(l1, l2, fmt=2)
multipolygon(l1, l2, fmt=2, third = "m")


# 4D
## data.frame
df <- data.frame(long = c(30, 45, 10, 30), lat = c(20, 40, 40, 20), 
  altitude = 1:4, weigjht = 20:23)
df2 <- data.frame(long = c(15, 40, 10, 5, 15), lat = c(5, 10, 20, 10, 5), 
  altitude = 1:5, weigjht = 200:204)
multipolygon(df, df2, fmt=0)
## matrix
mat <- matrix(unname(unlist(df)), ncol = 4)
mat2 <- matrix(unname(unlist(df2)), ncol = 4)
multipolygon(mat, mat2, fmt=0)
## list
l1 <- list(c(30, 20, 2, 1), c(45, 40, 2, 1), c(10, 40, 2, 1), c(30, 20, 2, 1))
l2 <- list(c(15, 5, 5, 1), c(40, 10, 5, 1), c(10, 20, 5, 1), c(5, 10, 5, 1), 
  c(15, 5, 5, 1))
multipolygon(l1, l2, fmt=2)
}
\seealso{
Other R-objects: 
\code{\link{circularstring}()},
\code{\link{geometrycollection}()},
\code{\link{linestring}()},
\code{\link{multilinestring}()},
\code{\link{multipoint}()},
\code{\link{point}()},
\code{\link{polygon}()}
}
\concept{R-objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson2wkt.R
\name{geojson2wkt}
\alias{geojson2wkt}
\title{Convert GeoJSON-like objects to WKT}
\usage{
geojson2wkt(obj, fmt = 16, third = "z", ...)
}
\arguments{
\item{obj}{(list/json/character) A GeoJSON-like object representing a
Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon,
or GeometryCollection}

\item{fmt}{Format string which indicates the number of digits to display
after the decimal point when formatting coordinates. Max: 20}

\item{third}{(character) Only applicable when there are three dimensions.
If \code{m}, assign a \code{M} value for a measurement, and if \code{z} assign a
\code{Z} value for three-dimenionsal system. Case is ignored. An \code{M} value
represents  a measurement, while a \code{Z} value usually represents altitude
(but can be something like depth in a water based location).}

\item{...}{Further args passed on to \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} only
in the event of json passed as a character string (can also be json
of class \code{json} as returned from \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}} or simply coerced
to \code{json} by adding the class manually)}
}
\description{
Convert GeoJSON-like objects to WKT
}
\section{Inputs}{

Input to \code{obj} parameter can take two forms:
\itemize{
\item A list with named elements \code{type} and \code{coordinates} OR
\code{type} and \code{geometries} (only in the case of GeometryCollection).
e.g., \code{list(type = "Point", coordinates = c(1, 0))}
\item A list with single named element in the set Point, Multipoint,
Polygon, Multipolygon, Linestring, Multilinestring,or Geometrycollection,
e.g., \code{list(Point = c(1, 0))} - Note that this format is not proper
GeoJSON, but is less verbose than the previous format, so should save
the user time and make it easier to use.
}
}

\section{Each point}{

For any one point, 2 to 4 values can be used:
\itemize{
\item 2 values: longitude, latitude
\item 3 values: longitude, latitude, altitude
\item 4 values: longitude, latitude, altitude, measure
}

The 3rd value is typically altitude though can be depth in an
aquatic context.

The 4th value is a measurement of some kind.

The GeoJSON spec \url{https://tools.ietf.org/html/rfc7946} actually
doesn't allow a 4th value for a point, but we allow it here since
we're converting to WKT which does allow a 4th value for a point.
}

\section{Coordinates data formats}{

Coordinates data should follow the following formats:
\itemize{
\item Point: a vector or list, with a single point (2-4 values)
\item MultiPoint: a matrix, with N points
\item Linestring: a matrix, with N points
\item MultiLinestring: the top most level is a list, containing N matrices
\item Polygon: the top most level is a list, containing N matrices
\item MultiPolygon: the top most level is a list, the next level is N
lists, each of them containing N matrices
\item Geometrycollection: a list containing any combination and number
of the above types
}

Matrices by definition can not have unequal lengths in their columns,
so we don't have to check for that user error.

Each matrix can have any number of rows, and from 2 to 4 columns.
If > 5 columns we stop with an error message.
}

\examples{
# point
## new format
point <- list(Point = c(116.4, 45.2))
geojson2wkt(point)
## old format, warns
point <- list(type = 'Point', coordinates = c(116.4, 45.2))
geojson2wkt(point)

# multipoint
## new format
mp <- list(MultiPoint = matrix(c(100, 101, 3.14, 3.101, 2.1, 2.18),
   ncol = 2))
geojson2wkt(mp)
## 3D
mp <- list(MultiPoint = matrix(c(100, 101, 3, 3, 2, 2, 4, 5, 6),
   ncol = 3))
geojson2wkt(mp)
## old format, warns
mp <- list(
  type = 'MultiPoint',
  coordinates = matrix(c(100, 101, 3.14, 3.101, 2.1, 2.18), ncol = 2)
)
geojson2wkt(mp)

# linestring
## new format
st <- list(LineString = matrix(c(0.0, 2.0, 4.0, 5.0,
                               0.0, 1.0, 2.0, 4.0),
                               ncol = 2))
geojson2wkt(st)
## 3D
st <- list(LineString = matrix(
  c(0.0, 0, 0, 
   2, 1, 5,
   100, 300, 800), nrow = 3))
geojson2wkt(st, fmt = 2)
geojson2wkt(st, fmt = 2, third = "m")
## old format, warns
st <- list(
  type = 'LineString',
  coordinates = matrix(c(0.0, 2.0, 4.0, 5.0,
                         0.0, 1.0, 2.0, 4.0), ncol = 2)
)
geojson2wkt(st)
## 3D
st <- list(LineString = matrix(c(0.0, 2.0, 4.0, 5.0,
                               0.0, 1.0, 2.0, 4.0,
                               10, 20, 30, 40),
                               ncol = 3))
geojson2wkt(st, fmt = 2)

## 4D
st <- list(LineString = matrix(c(0.0, 2.0, 4.0, 5.0,
                               0.0, 1.0, 2.0, 4.0,
                               10, 20, 30, 40,
                               1, 2, 3, 4),
                               ncol = 4))
geojson2wkt(st, fmt = 2)


# multilinestring
## new format
multist <- list(MultiLineString = list(
   matrix(c(0, -2, -4, -1, -3, -5), ncol = 2),
   matrix(c(1.66, 10.9999, 10.9, 0, -31.5, 3.0, 1.1, 0), ncol = 2)
 )
)
geojson2wkt(multist)
## 3D
multist <- list(MultiLineString = list(
   matrix(c(0, -2, -4, -1, -3, -5, 100, 200, 300), ncol = 3),
   matrix(c(1, 10, 10.9, 0, -31.5, 3.0, 1.1, 0, 3, 3, 3, 3), ncol = 3)
 )
)
geojson2wkt(multist, fmt = 2)
geojson2wkt(multist, fmt = 2, third = "m")
## old format, warns
multist <- list(
  type = 'MultiLineString',
  coordinates = list(
   matrix(c(0, -2, -4, -1, -3, -5), ncol = 2),
   matrix(c(1.66, 10.9999, 10.9, 0, -31.5, 3.0, 1.1, 0), ncol = 2)
 )
)
geojson2wkt(multist)

## points within MultiLineString that differ 
## -> use length of longest 
## -> fill with zeros
# 3D and 2D
multist <- list(MultiLineString = list(
    matrix(1:6, ncol = 3), matrix(1:8, ncol = 2)))
geojson2wkt(multist, fmt = 0)
# 4D and 2D
multist <- list(MultiLineString = list(
    matrix(1:8, ncol = 4), matrix(1:8, ncol = 2)))
geojson2wkt(multist, fmt = 0)
# 2D and 2D
multist <- list(MultiLineString = list(
    matrix(1:4, ncol = 2), matrix(1:8, ncol = 2)))
geojson2wkt(multist, fmt = 0)
# 5D and 2D - FAILS
# multist <- list(MultiLineString = list(
#    matrix(1:10, ncol = 5), matrix(1:8, ncol = 2)))
# geojson2wkt(multist, fmt = 0)


# polygon
## new format
poly <- list(Polygon = list(
   matrix(c(100.001, 101.1, 101.001, 100.001, 0.001, 0.001, 1.001, 0.001), ncol = 2),
   matrix(c(100.201, 100.801, 100.801, 100.201, 0.201, 0.201, 0.801, 0.201), ncol = 2)
))
geojson2wkt(poly)
geojson2wkt(poly, fmt=6)
## 3D
poly <- list(Polygon = list(
   matrix(c(100.1, 101.1, 101.1, 100.1, 0.1, 0.1, 1.1, 0.1, 1, 1, 1, 1), ncol = 3),
   matrix(c(100.2, 100.8, 100.8, 100.2, 0.2, 0.2, 0.80, 0.2, 3, 3, 3, 3), ncol = 3)
))
geojson2wkt(poly, fmt = 2)
geojson2wkt(poly, fmt = 2, third = "m")
## old format, warns
poly <- list(
  type = 'Polygon',
  coordinates = list(
    matrix(c(100.001, 101.1, 101.001, 100.001, 0.001, 0.001, 1.001, 0.001),
      ncol = 2),
    matrix(c(100.201, 100.801, 100.801, 100.201, 0.201, 0.201, 0.801, 0.201),
      ncol = 2)
  )
)
geojson2wkt(poly)
geojson2wkt(poly, fmt=6)

## points within Polygon that differ 
## -> use length of longest 
## -> fill with zeros
# 3D and 2D
poly <- list(Polygon = list(
   matrix(c(100, 101, 101, 100, 0.1, 0.2, 0.3, 0.1, 5, 6, 7, 8), ncol = 3),
   matrix(c(40, 41, 61, 40, 0.1, 0.2, 0.3, 0.1), ncol = 2)
))
geojson2wkt(poly, fmt = 0)
# 4D and 2D
poly <- list(Polygon = list(
   matrix(c(100, 101, 101, 100, 0.1, 0.2, 0.3, 0.1, 5, 6, 7, 8, 1, 1, 1, 1), 
     ncol = 4),
   matrix(c(40, 41, 61, 40, 0.1, 0.2, 0.3, 0.1), ncol = 2)
))
geojson2wkt(poly, fmt = 0)
# 5D and 2D - FAILS
# multist <- list(Polygon = list(
#    matrix(1:10, ncol = 5), matrix(1:8, ncol = 2)))
# geojson2wkt(poly, fmt = 0)



# multipolygon
## new format
mpoly <- list(MultiPolygon = list(
  list(
    matrix(c(100, 101, 101, 100, 0.001, 0.001, 1.001, 0.001), ncol = 2),
    matrix(c(100.2, 100.8, 100.8, 100.2, 0.2, 0.2, 0.8, 0.2), ncol = 2)
  ),
  list(
    matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0), ncol = 2),
    matrix(c(9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0), ncol = 2)
  )
))
geojson2wkt(mpoly, fmt=2)
## 3D
mpoly <- list(MultiPolygon = list(
  list(
    matrix(c(100, 101, 101, 100, 0.001, 0.001, 1.001, 0.001, 1, 1, 1, 1), 
      ncol = 3),
    matrix(c(100.2, 100.8, 100.8, 100.2, 0.2, 0.2, 0.8, 0.2, 3, 4, 5, 6), 
      ncol = 3)
  ),
  list(
    matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 1, 1, 1, 1),
      ncol = 3),
    matrix(c(9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 9, 9, 9, 9),
      ncol = 3)
  )
))
geojson2wkt(mpoly, fmt=2)
geojson2wkt(mpoly, fmt=2, third = "m")
## old format, warns
mpoly <- list(
  type = 'MultiPolygon',
  coordinates = list(
    list(
      matrix(c(100, 101, 101, 100, 0.001, 0.001, 1.001, 0.001), ncol = 2),
      matrix(c(100.2, 100.8, 100.8, 100.2, 0.2, 0.2, 0.8, 0.2), ncol = 2)
    ),
    list(
      matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 1.0), ncol = 3),
      matrix(c(9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 9.0), ncol = 3)
    )
  )
)
geojson2wkt(mpoly, fmt=2)

mpoly2 <- list(
  type = "MultiPolygon",
  coordinates = list(
    list(list(c(30, 20), c(45, 40), c(10, 40), c(30, 20))),
    list(list(c(15, 5), c(40, 10), c(10, 20), c(5 ,10), c(15, 5)))
  )
)

mpoly2 <- list(
  type = "MultiPolygon",
  coordinates = list(
    list(
      matrix(c(30, 45, 10, 30, 20, 40, 40, 20), ncol = 2)
    ),
    list(
      matrix(c(15, 40, 10, 5, 15, 5, 10, 20, 10, 5), ncol = 2)
    )
  )
)
geojson2wkt(mpoly2, fmt=1)

## points within MultiPolygon that differ 
## -> use length of longest 
## -> fill with zeros
# 3D and 2D
mpoly <- list(MultiPolygon = list(
  list(
    matrix(c(40, 130, 155, 40, 20, 34, 34, 20), ncol = 2),
    matrix(c(30, 40, 54, 30, 0.1, 42, 62, 0.1, 1, 1, 1, 1), ncol = 3)
  ),
  list(
    matrix(c(9, 49, 79, 9, 11, 35, 15, 11), ncol = 2),
    matrix(c(1, 33, 59, 1, 5, 16, 36, 5), ncol = 2)
  )
))
geojson2wkt(mpoly, fmt = 0)
# 4D and 2D
mpoly <- list(MultiPolygon = list(
  list(
    matrix(c(40, 130, 155, 40, 20, 34, 34, 20), ncol = 2),
    matrix(c(30, 40, 54, 30, 0.1, 42, 62, 0.1, 1, 1, 1, 1, 0, 0, 0, 0), ncol = 4)
  ),
  list(
    matrix(c(9, 49, 79, 9, 11, 35, 15, 11), ncol = 2),
    matrix(c(1, 33, 59, 1, 5, 16, 36, 5), ncol = 2)
  )
))
geojson2wkt(mpoly, fmt = 0)
# 5D and 2D - FAILS
mpoly <- list(MultiPolygon = list(
  list(
    matrix(c(40, 130, 155, 40, 20, 34, 34, 20), ncol = 2),
    matrix(c(30, 40, 54, 30, 
             0.1, 42, 62, 0.1, 
             1, 1, 1, 1, 
             0, 0, 0, 0,
             0, 0, 0, 0), ncol = 5)
  ),
  list(
    matrix(c(9, 49, 79, 9, 11, 35, 15, 11), ncol = 2),
    matrix(c(1, 33, 59, 1, 5, 16, 36, 5), ncol = 2)
  )
))
# geojson2wkt(mpoly, fmt = 0)



# geometrycollection
## new format
gmcoll <- list(GeometryCollection = list(
 list(Point = c(0.0, 1.0)),
 list(LineString = matrix(c(0.0, 2.0, 4.0, 5.0,
                           0.0, 1.0, 2.0, 4.0),
                           ncol = 2)),
 list(Polygon = list(
   matrix(c(100.001, 101.1, 101.001, 100.001, 0.001, 0.001, 1.001, 0.001),
     ncol = 2),
   matrix(c(100.201, 100.801, 100.801, 100.201, 0.201, 0.201, 0.801, 0.201),
     ncol = 2)
  ))
))
geojson2wkt(gmcoll, fmt=0)
## old format, warns
gmcoll <- list(
 type = 'GeometryCollection',
 geometries = list(
   list(type = 'Point', coordinates = c(0.0, 1.0)),
   list(type = 'LineString', coordinates = matrix(c(0.0, 2.0, 4.0, 5.0,
                           0.0, 1.0, 2.0, 4.0),
                           ncol = 2)),
   list(type = 'Polygon', coordinates = list(
     matrix(c(100.001, 101.1, 101.001, 100.001, 0.001, 0.001, 1.001, 0.001),
       ncol = 2),
     matrix(c(100.201, 100.801, 100.801, 100.201, 0.201, 0.201, 0.801, 0.201),
       ncol = 2)
  ))
 )
)
geojson2wkt(gmcoll, fmt=0)

# Convert geojson as character string to WKT
# new format
str <- '{
   "Point": [
       -105.01621,
       39.57422
   ]
}'
geojson2wkt(str)

## old format, warns
str <- '{
   "type": "Point",
   "coordinates": [
       -105.01621,
       39.57422
   ]
}'
geojson2wkt(str)

## new format
str <- '{"LineString":[[0,0,10],[2,1,20],[4,2,30],[5,4,40]]}'
geojson2wkt(str)
## old format, warns
str <-
'{"type":"LineString","coordinates":[[0,0,10],[2,1,20],[4,2,30],[5,4,40]]}'
geojson2wkt(str)

# From a jsonlite json object
library("jsonlite")
json <- toJSON(list(Point=c(-105,39)), auto_unbox=TRUE)
geojson2wkt(json)
## old format, warns
json <- toJSON(list(type="Point", coordinates=c(-105,39)), auto_unbox=TRUE)
geojson2wkt(json)
}
\references{
\url{https://tools.ietf.org/html/rfc7946},
\url{https://en.wikipedia.org/wiki/Well-known_text}
}
\seealso{
\code{\link[=wkt2geojson]{wkt2geojson()}}
}
