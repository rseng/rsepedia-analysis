parzer
================

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran
checks](https://cranchecks.info/badges/worst/parzer)](https://cranchecks.info/pkgs/parzer)
[![R-CMD-check](https://github.com/ropensci/parzer/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/parzer/actions/)
[![codecov.io](https://codecov.io/github/ropensci/parzer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/parzer?branch=master)
[![](https://badges.ropensci.org/341_status.svg)](https://github.com/ropensci/software-review/issues/341)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/parzer?color=C9A115)](https://github.com/r-hub/cranlogs.app)
[![cran
version](https://www.r-pkg.org/badges/version/parzer)](https://cran.r-project.org/package=parzer)

`parzer` parses messy geographic coordinates

Docs: <https://docs.ropensci.org/parzer/>

You may get data from a published study or a colleague, and the
coordinates may be in some messy character format that you’d like to
clean up to have all decimal degree numeric data.

`parzer` API:

-   `parse_hemisphere`
-   `parse_lat`
-   `parse_llstr`
-   `parse_lon`
-   `parse_lon_lat`
-   `parse_parts_lat`
-   `parse_parts_lon`
-   `pz_d`
-   `pz_degree`
-   `pz_m`
-   `pz_minute`
-   `pz_s`
-   `pz_second`

## Usage

For example, parse latitude and longitude from messy character vectors.

``` r
parse_lat(c("45N54.2356", "-45.98739874", "40.123°"))
#> [1]  45.90393 -45.98740  40.12300
```

``` r
parse_lon(c("45W54.2356", "-45.98739874", "40.123°"))
#> [1] -45.90393 -45.98740  40.12300
```

See more in the [Introduction to the `parzer` package
vignette](https://docs.ropensci.org/parzer/articles/parzer.html).

## Installation

Stable version

``` r
install.packages("parzer")
```

Development version

``` r
remotes::install_github("ropensci/parzer")
```

``` r
library("parzer")
```

## Similar art

-   `sp::char2dms`: is most similar to `parzer::parse_lat` and
    `parzer::parse_lon`. However, with `sp::char2dms` you have to
    specify the termination character for each of degree, minutes and
    seconds. `parzer` does this for the user.
-   `biogeo::dms2dd`: very unlike functions in this package. You must
    pass separate degrees, minutes, seconds and direction to `dms2dd`.
    No exact analog is found in `parzer`, whose main focus is parsing
    messy geographic coordinates in strings to a more machine readable
    version

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/parzer/issues).
-   License: MIT
-   Get citation information for `parzer` in R doing
    `citation(package = 'parzer')`
-   Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
parzer 0.4.1
============

### MINOR IMPROVEMENTS

* documentation and package description describe more clearly `parzer` core objective of parsing messy coordinates in 
character strings to convert them to decimal numeric values. Suggestion and work by @robitalec

### ACKNOWLEDGEMENTS CHANGES
* new contributors to the package: @robitalec, @maelle and @yutannihilation
* new maintainer: @AlbanSagouis

parzer 0.4.0
============

### MINOR IMPROVEMENTS

* performance improvement for internal function `scrub()`, used in most exported functions in parzer (#30) work by @AlbanSagouis
* work around for non-UTF8 MBCS locales: now all exported functions go through a modified `.Call()` in which we use `withr::with_locale()` if the user is on a Windows operating system (#31) (#32) work by @yutannihilation

parzer 0.3.0
============

### BUG FIXES

* fix problem in `parse_llstr()`: on older R versions where `stringsAsFactors=TRUE` by default this function was returning strings as factors from an internal function that caused a problem in a subsequent step in the function (#29)

parzer 0.2.0
============

### NEW FEATURES

* new contributor to the package @AlbanSagouis
* gains new function `parse_llstr()` to parse a string that contains both latitude and longitude (#3) (#24) (#26) (#28) work by @AlbanSagouis

### MINOR IMPROVEMENTS

* updated `scrub()` internal function that strips certain characters to include more things to scrub (#25) work by @AlbanSagouis

parzer 0.1.4
============

### MINOR IMPROVEMENTS

* add support to internal function for additional degree like symbols (#21)
* fix issue with `parse_parts_lat()`/`parse_parts_lon()` functions where an NA was causing warnings on the cpp side; on cpp side, now check for NA and return list of NAs instead of NAs passing through other code (#23)

parzer 0.1.0
============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local win install, R 4.1.1
* on github actions
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (devel)
  - ubuntu-latest (release)
  - ubuntu-latest (oldrel-1)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

---

This version 0.4.1.

Former maintainer Scott Chamberlain (sckott@protonmail.com)
New maintainer: Alban Sagouis (sagouis@pm.me)  
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

* Submit an issue on the [Issues page](https://github.com/ropensci/parzer/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/parzer.git`
* Make sure to track progress upstream (i.e., on our version of `parzer` at `ropensci/parzer`) by doing `git remote add upstream https://github.com/ropensci/parzer.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/parzer`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-09-01 r79114) |
|os       |macOS Catalina 10.15.6                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-10-07                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|parzer  |0.1.4 |0.2.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
title: "parzer"
output: github_document
---


```{r echo=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/parzer)](https://cranchecks.info/pkgs/parzer)
[![R-CMD-check](https://github.com/ropensci/parzer/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/parzer/actions/)
[![codecov.io](https://codecov.io/github/ropensci/parzer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/parzer?branch=master)
[![](https://badges.ropensci.org/341_status.svg)](https://github.com/ropensci/software-review/issues/341)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/parzer?color=C9A115)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/parzer)](https://cran.r-project.org/package=parzer)

`parzer` parses messy geographic coordinates

Docs: https://docs.ropensci.org/parzer/

You may get data from a published study or a colleague, and the coordinates
may be in some messy character format that you'd like to clean up to have 
all decimal degree numeric data.

`parzer` API:

```{r echo=FALSE, comment=NA, results='asis'}
cat(paste(" -", paste(sprintf("`%s`", sort(getNamespaceExports("parzer"))), collapse = "\n - ")))
```

## Usage
For example, parse latitude and longitude from messy character vectors. 

```{r eval=TRUE, echo=FALSE}
library("parzer")
```

```{r}
parse_lat(c("45N54.2356", "-45.98739874", "40.123°"))
```

```{r}
parse_lon(c("45W54.2356", "-45.98739874", "40.123°"))
```

See more in the [Introduction to the `parzer` package vignette](https://docs.ropensci.org/parzer/articles/parzer.html). 

## Installation

Stable version

```{r eval=FALSE}
install.packages("parzer")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/parzer")
```

```{r}
library("parzer")
```


## Similar art

- `sp::char2dms`: is most similar to `parzer::parse_lat` and `parzer::parse_lon`. However,
with `sp::char2dms` you have to specify the termination character for each of degree,
minutes and seconds. `parzer` does this for the user.
- `biogeo::dms2dd`: very unlike functions in this package. You must pass separate degrees,
minutes, seconds and direction to `dms2dd`. No exact analog is found in `parzer`, whose
main focus is parsing messy geographic coordinates in strings to a more machine readable
version

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/parzer/issues).
* License: MIT
* Get citation information for `parzer` in R doing `citation(package = 'parzer')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction to the parzer package"
author: "Scott Chamberlain"
date: "2020-03-26"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to the parzer package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




`parzer` parses messy coordinates

You may get data from a published study or a colleague, and the coordinates
may be in some messy character format that you'd like to clean up to have 
all decimal degree numeric data.

`parzer` API:

 - `parse_hemisphere`
 - `parse_lat`
 - `parse_lon`
 - `parse_lon_lat`
 - `parse_parts_lat`
 - `parse_parts_lon`
 - `pz_d`
 - `pz_degree`
 - `pz_m`
 - `pz_minute`
 - `pz_s`
 - `pz_second`


## Install

Stable version

```r
install.packages("parzer")
```

Development version


```r
remotes::install_github("ropensci/parzer")
```


## parse

```r
library("parzer")
```

### latitudes


```r
parse_lat("45N54.2356")
#> [1] 45.90393
parse_lat("-45.98739874")
#> [1] -45.9874
parse_lat("40.123°")
#> [1] 40.123
parse_lat("40.123N")
#> [1] 40.123
parse_lat("N45 04.25764")
#> [1] 45.07096

# bad values -> NaN
parse_lat("191.89")
#> Warning in pz_parse_lat(lat): not within -90/90 range, got: 191.89
#>   check that you did not invert lon and lat
#> [1] NaN

# many inputs
x <- c("40.123°", "40.123N74.123W", "191.89", 12, "N45 04.25764")
parse_lat(x)
#> Warning in pz_parse_lat(lat): invalid characters, got: 40.123n74.123w

#> Warning in pz_parse_lat(lat): not within -90/90 range, got: 191.89
#>   check that you did not invert lon and lat
#> [1] 40.12300      NaN      NaN 12.00000 45.07096

# parse_lat("N455698735", "HDDMMmmmmm") # custom formats not ready yet
```

### longitudes


```r
parse_lon("45W54.2356")
#> [1] -45.90393
parse_lon("-45.98739874")
#> [1] -45.9874
parse_lon("40.123°")
#> [1] 40.123
parse_lon("74.123W")
#> [1] -74.123
parse_lon("W45 04.25764")
#> [1] -45.07096

# bad values
parse_lon("361")
#> Warning in pz_parse_lon(lon): not within -180/360 range, got: 361
#> [1] NaN

# many inputs
x <- c("45W54.2356", "181", 45, 45.234234, "-45.98739874")
parse_lon(x)
#> [1] -45.90393 181.00000  45.00000  45.23423 -45.98740

# parse_lon("N455698735", "HDDMMmmmmm") # custom formats not ready yet
```

### both lon and lat together


```r
lons <- c("45W54.2356", "181", 45, 45.234234, "-45.98739874")
lats <- c("40.123°", "40.123N", "191.89", 12, "N45 04.25764")
parse_lon_lat(lons, lats)
#>         lon      lat
#> 1 -45.90393 40.12300
#> 2 181.00000 40.12300
#> 3  45.00000      NaN
#> 4  45.23423 12.00000
#> 5 -45.98740 45.07096
```

### parse into degree, min, sec parts


```r
parse_parts_lat("45N54.2356")
#>   deg min      sec
#> 1  45  54 14.13674
parse_parts_lon("-74.6411133")
#>   deg min      sec
#> 1 -74  38 28.00784
# many inputs
x <- c("40.123°", "40.123N74.123W", "191.89", 12, "N45 04.25764")
parse_parts_lon(x)
#>   deg min      sec
#> 1  40   7 22.80395
#> 2  NA  NA      NaN
#> 3 191  53 23.99783
#> 4  12   0  0.00000
#> 5  NA  NA      NaN
```

### get hemisphere from lat/lon coords


```r
parse_hemisphere("74.123E", "45N54.2356")
#> [1] "NE"
parse_hemisphere("-120", "40.4183318")
#> [1] "NW"
parse_hemisphere("-120", "-40.4183318")
#> [1] "SW"
parse_hemisphere("120", "-40.4183318")
#> [1] "SE"
```

### get degree, minutes, or seconds separately


```r
coords <- c(45.23323, "40:25:6N", "40° 25´ 5.994\" N")
pz_degree(lat = coords)
#> [1] 45 40 40
pz_minute(lat = coords)
#> [1] 13 25 25
pz_second(lat = coords)
#> [1] 59.630119  6.005895  5.992162

coords <- c(15.23323, "40:25:6E", "192° 25´ 5.994\" E")
pz_degree(lon = coords)
#> [1]  15  40 192
pz_minute(lon = coords)
#> [1] 13 25 25
pz_second(lon = coords)
#> [1] 59.626686  6.005895  6.005895
```

### add or subtract degrees, minutes, or seconds


```r
pz_d(31)
#> 31
pz_d(31) + pz_m(44)
#> 31.73333
pz_d(31) - pz_m(44)
#> 30.26667
pz_d(31) + pz_m(44) + pz_s(59)
#> 31.74972
pz_d(-121) + pz_m(1) + pz_s(33)
#> -120.9742
```
---
title: "parzer use cases"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{parzer use cases}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```

## Use case: working with spatial R packages

```{r}
library("parzer")
if (!requireNamespace("sf")) install.packages("sf")
library("sf")
```

One may find themselves having to clean up messy coordinates as part of
their project/work/etc. How does this look when fit into a workflow going 
all the way to visualization?

Let's say you have the following messy coordinates that you've compiled
from different places, leading to a variety of messy formats:

```{r}
lats <- c(
  "46.4183",
  "46.4383° N",
  "46.5683° N",
  "46° 27´ 5.4\" N",
  "46° 25.56’",
  "N46°24’4.333"
)
lons <- c(
  "25.7391",
  "E25°34’6.4533",
  "25.3391° E",
  "25.8391° E",
  "25° 35.56’",
  "E25°34’4.333"
)
```

Parse messy coordinates


```{r}
dat <- tibble::tibble(
  longitude = parse_lon(lons),
  latitude = parse_lat(lats)
)
dat
```

Combine coordinates with other data


```{r}
dat$shape <- c("round", "square", "triangle", "round", "square", "square")
dat$color <- c("blue", "yellow", "green", "red", "green", "yellow")
dat
```

Coerce to an sf object

```{r}
datsf <- sf::st_as_sf(dat, coords = c("longitude", "latitude"))
datsf
```

Calculate the center of the plot view

```{r}
center_lon <- mean(dat$longitude)
center_lat <- mean(dat$latitude)
```

Plot data using the `leaflet` package

```{r out.with="100%"}
if (!requireNamespace("leaflet")) install.packages("leaflet")
library("leaflet")
leaflet() %>%
  addTiles() %>%
  addMarkers(data = datsf) %>%
  setView(center_lon, center_lat, zoom = 10)
```

We'd like to have data only for a certain area, e.g., a political boundary or
a park boundary. We can clip the data to a bounding box using `sf::st_crop()`.

First, define the bounding box, and visualize


```{r}
bbox <- c(
  xmin = 25.42813, ymin = 46.39455,
  xmax = 25.68769, ymax = 46.60346
)
leaflet() %>%
  addTiles() %>%
  addRectangles(bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]]) %>%
  setView(center_lon, center_lat, zoom = 10)
```

Crop the data to the bounding box

```{r}
datsf_c <- st_crop(datsf, bbox)
```

Plot the new data


```{r out.with="100%"}
leaflet() %>%
  addTiles() %>%
  addMarkers(data = datsf_c) %>%
  setView(center_lon, center_lat, zoom = 10)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dms-fxns.R
\name{dms}
\alias{dms}
\alias{pz_degree}
\alias{pz_minute}
\alias{pz_second}
\alias{print.pz}
\alias{pz_d}
\alias{pz_m}
\alias{pz_s}
\alias{+.pz}
\alias{-.pz}
\alias{/.pz}
\alias{*.pz}
\title{extract degree, minutes, and seconds}
\usage{
pz_degree(lon = NULL, lat = NULL)

pz_minute(lon = NULL, lat = NULL)

pz_second(lon = NULL, lat = NULL)

\method{print}{pz}(x, ...)

pz_d(x)

pz_m(x)

pz_s(x)

\method{+}{pz}(e1, e2)

\method{-}{pz}(e1, e2)

\method{/}{pz}(e1, e2)

\method{*}{pz}(e1, e2)
}
\arguments{
\item{lon, lat}{(numeric/integer/character) one or more longitude or
latitude values. values are internally validated. only one of
lon or lat accepted}

\item{x}{(integer) an integer representing a degree, minute or second}

\item{...}{print dots}

\item{e1, e2}{objects of class pz, from using \code{pz_d()}, \code{pz_m()}, or \code{pz_s()}}
}
\value{
\code{pz_degree}: integer, \code{pz_minute}: integer, \code{pz_second}: numeric,
\code{pz_d}: numeric, \code{pz_m}: numeric, \code{pz_s}: numeric (adding/subtracting
these also gives numeric)
}
\description{
extract degree, minutes, and seconds
}
\details{
Mathematics operators are exported for \code{+}, \code{-}, \code{/}, and \code{*},
but \code{/} and \code{*} are only exported with a stop message to say it's not
supported; otherwise you'd be allow to divide degrees by minutes, leading
to nonsense.
}
\examples{
# extract parts of a coordinate value
pz_degree(-45.23323)
pz_minute(-45.23323)
pz_second(-45.23323)

pz_degree(lon = 178.23423)
pz_minute(lon = 178.23423)
pz_second(lon = 178.23423)
\dontrun{
pz_degree(lat = c(45.23323, "40:25:6N", "40° 25´ 5.994 S"))
pz_minute(lat = c(45.23323, "40:25:6N", "40° 25´ 5.994 S"))
pz_second(lat = c(45.23323, "40:25:6N", "40° 25´ 5.994 S"))

# invalid
pz_degree(445.23323)

# add or subtract
pz_d(31)
pz_m(44)
pz_s(3)
pz_d(31) + pz_m(44)
pz_d(-31) - pz_m(44)
pz_d(-31) + pz_m(44) + pz_s(59)
pz_d(31) - pz_m(44) + pz_s(59)
pz_d(-121) + pz_m(1) + pz_s(33)
unclass(pz_d(31) + pz_m(44) + pz_s(59))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parzer-package.R
\docType{package}
\name{parzer-package}
\alias{parzer-package}
\alias{parzer}
\title{parzer}
\description{
parse geographic coordinates
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_lat.R
\name{parse_lat}
\alias{parse_lat}
\title{Parse latitude values}
\usage{
parse_lat(lat, format = NULL)
}
\arguments{
\item{lat}{(numeric/integer/character) one or more latitude values}

\item{format}{(character) format, default often works}
}
\value{
numeric vector
}
\description{
Parse latitude values
}
\section{Errors}{

Throws warnings on parsing errors, and returns \code{NaN} in each case

Types of errors:
\itemize{
\item invalid argument: e.g., letters passed instead of numbers,
see \url{https://en.cppreference.com/w/cpp/error/invalid_argument}
\item out of range: numbers of out acceptable range, see
\url{https://en.cppreference.com/w/cpp/error/out_of_range}
\item out of latitude range: not within -90/90 range
}
}

\examples{
parse_lat("")
\dontrun{
parse_lat("-91")
parse_lat("95")
parse_lat("asdfaf")

parse_lat("45")
parse_lat("-45")
parse_lat("-45.2323")

# out of range with std::stod?
parse_lat("-45.23232e24")
parse_lat("-45.23232e2")

# numeric input
parse_lat(1:10)
parse_lat(85:94)

# different formats
parse_lat("40.4183318 N")
parse_lat("40.4183318 S")
parse_lat("40 25 5.994") # => 40.41833

parse_lat("40.4183318N")
parse_lat("N40.4183318")
parse_lat("40.4183318S")
parse_lat("S40.4183318")

parse_lat("N 39 21.440") # => 39.35733
parse_lat("S 56 1.389") # => -56.02315

parse_lat("N40°25’5.994") # => 40.41833
parse_lat("40° 25´ 5.994\" N") # => 40.41833
parse_lat("40:25:6N")
parse_lat("40:25:5.994N")
parse_lat("40d 25’ 6\" N")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_hemisphere.R
\name{parse_hemisphere}
\alias{parse_hemisphere}
\title{get hemisphere from long/lat coordinates}
\usage{
parse_hemisphere(lon, lat)
}
\arguments{
\item{lon}{(character/numeric/integer) one or more longitude values}

\item{lat}{(character/numeric/integer) one or more latitude values}
}
\value{
character vector of quadrants, one of: NE, NW, SE, SW.
if one of the coordinate values is invalid, and one is valid, you get
a length 1 string. if both coordinate values are bad, you get
a zero length string.

Warnings are thrown on invalid values
}
\description{
BEWARE: EXPERIMENTAL
}
\details{
length(lon) == length(lat)
}
\examples{
# NE
parse_hemisphere("74.123E", "45N54.2356")
\dontrun{
# NW
parse_hemisphere(-120, 40.4183318)
# SW
parse_hemisphere(-120, -40.4183318)
# SE
parse_hemisphere(120, -40.4183318)

# bad inputs, get one of the two strings
parse_hemisphere(-181, -40.4183318)
parse_hemisphere(-120, -192.4183318)

# many inputs
library(randgeo)
pts <- rg_position(count = 1000)
lons <- as.character(vapply(pts, "[[", 1, 1))
lats <- as.character(vapply(pts, "[[", 1, 2))
parse_hemisphere(lons, lats)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_llstr.R
\name{parse_llstr}
\alias{parse_llstr}
\title{parse string with lat and lon together}
\usage{
parse_llstr(str)
}
\arguments{
\item{str}{(character) string with latitude and longitude, one or more in a vector.}
}
\value{
A data.frame with parsed latitude and longitude in decimal degrees.
}
\description{
parse string with lat and lon together
}
\examples{
parse_llstr("N 04.1683, E 101.5823")
parse_llstr("N04.82344, E101.61320")
parse_llstr("N 04.25164, E 101.70695")
parse_llstr("N05.03062, E101.75172")
parse_llstr("N05.03062,E101.75172")
parse_llstr("N4.9196, E101.345")
parse_llstr("N4.9196, E101.346")
parse_llstr("N4.9196, E101.347")
# no comma
parse_llstr("N4.9196 E101.347")
# no space
parse_llstr("N4.9196E101.347")

# DMS
parse_llstr("N4 51'36\", E101 34'7\"")
parse_llstr(c("4 51'36\"S, 101 34'7\"W", "N4 51'36\", E101 34'7\""))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_lat_lon.R
\name{parse_lon_lat}
\alias{parse_lon_lat}
\title{parse longitude and latitude}
\usage{
parse_lon_lat(lon, lat)
}
\arguments{
\item{lon}{(character/numeric/integer) one or more longitude values}

\item{lat}{(character/numeric/integer) one or more latitude values}
}
\value{
data.frame, with columns lon, lat. on an invalid values, an \code{NA}
is returned. In addition, warnings are thrown on invalid values
}
\description{
parse longitude and latitude
}
\details{
length(lon) == length(lat)
}
\examples{
parse_lon_lat(-120.43, 49.12)
\dontrun{
parse_lon_lat(-120.43, 93)
parse_lon_lat(-190, 49.12)
parse_lon_lat(240, 49.12)
parse_lon_lat(-190, 92)
# many
lons <- c("45W54.2356", "181", 45, 45.234234, "-45.98739874")
lats <- c("40.123°", "40.123N74.123W", "191.89", 12, "N45 04.25764")
parse_lon_lat(lons, lats)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_lon.R
\name{parse_lon}
\alias{parse_lon}
\title{Parse longitude values}
\usage{
parse_lon(lon, format = NULL)
}
\arguments{
\item{lon}{(numeric/integer/character) one or more longitude values}

\item{format}{(character) format, default often works}
}
\value{
numeric vector
}
\description{
Parse longitude values
}
\section{Errors}{

Throws warnings on parsing errors, and returns \code{NaN} in each case

Types of errors:
\itemize{
\item invalid argument: e.g., letters passed instead of numbers,
see \url{https://en.cppreference.com/w/cpp/error/invalid_argument}
\item out of range: numbers of out acceptable range, see
\url{https://en.cppreference.com/w/cpp/error/out_of_range}
\item out of longitude range: not within -180/360 range
}
}

\examples{
parse_lon("")
\dontrun{
parse_lon("-181")
parse_lon("-361")
parse_lon("95")
parse_lon("asdfaf")

parse_lon("45")
parse_lon("-45")
parse_lon("-45.2323")
parse_lon("334")

# out of range with std::stod?
parse_lon("-45.23232e24")
parse_lon("-45.23232e2")
parse_lon("-45.23232")

# numeric input
parse_lon(1:10)
parse_lon(85:94)

# different formats
parse_lon("40.4183318 E")
parse_lon("40.4183318 W")
parse_lon("40 25 5.994") # => 40.41833

parse_lon("40.4183318W")
parse_lon("W40.4183318")
parse_lon("E40.4183318")
parse_lon("40.4183318E")

parse_lon("E 39 21.440") # => 39.35733
parse_lon("W 56 1.389") # => -56.02315

parse_lon("E40°25’5.994") # => 40.41833
parse_lon("40° 25´ 5.994\" E") # => 40.41833
parse_lon("40:25:6E")
parse_lon("40:25:5.994E")
parse_lon("40d 25’ 6\" E")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_parts.R
\name{parse_parts}
\alias{parse_parts}
\alias{parse_parts_lon}
\alias{parse_parts_lat}
\title{parse coordinates into degrees, minutes and seconds}
\usage{
parse_parts_lon(str)

parse_parts_lat(str)
}
\arguments{
\item{str}{(character) string including longitude or latitude}
}
\value{
data.frame with columns for:
\itemize{
\item deg (integer)
\item min (integer)
\item sec (numeric)
}

NA/NaN given upon error
}
\description{
parse coordinates into degrees, minutes and seconds
}
\examples{
parse_parts_lon("140.4183318")
\dontrun{
parse_parts_lon("174.6411133")
parse_parts_lon("-45.98739874")
parse_parts_lon("40.123W")

parse_parts_lat("45N54.2356")
parse_parts_lat("40.4183318")
parse_parts_lat("-74.6411133")
parse_parts_lat("-45.98739874")
parse_parts_lat("40.123N")
parse_parts_lat("N40°25’5.994")

# not working, needs format input
parse_parts_lat("N455698735")

# multiple
x <- c("40.123°", "40.123N74.123W", "191.89", 12, "N45 04.25764")
parse_parts_lat(x)
system.time(parse_parts_lat(rep(x, 10^2)))
}

}
