getlandsat
==========


[![cran checks](https://cranchecks.info/badges/worst/getlandsat)](https://cranchecks.info/pkgs/getlandsat)
[![Build Status](https://travis-ci.org/ropensci/getlandsat.svg?branch=master)](https://travis-ci.org/ropensci/getlandsat)
[![codecov](https://codecov.io/gh/ropensci/getlandsat/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/getlandsat)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/getlandsat)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/getlandsat)](https://cran.r-project.org/package=getlandsat)
[![](https://badges.ropensci.org/58_status.svg)](https://github.com/ropensci/onboarding/issues/58)

`getlandsat`: Get Landsat 8 data from AWS public data sets

`getlandsat` provides access to Landsat <https://landsat.usgs.gov> 8 metadata and images hosted on AWS S3 at <https://registry.opendata.aws/landsat-8/>. The package only fetches data. It does not attempt to aid users in downstream usage, but some additional functionality may be added.

A new function `lsat_search()` lets you search for Landsat images by using the API from Development Seed documented at <https://github.com/sat-utils/sat-api>

Potential users are probably anyone from scientists asking questions about biodiversity or land use change, to software developers creating tools for users to vizualize their data.

## Install

Stable version


```r
install.packages("getlandsat")
```

Dev version


```r
devtools::install_github("ropensci/getlandsat")
```


```r
library("getlandsat")
```

## Search for images

```{r eval=FALSE}
x <- lsat_search(collection = "landsat-8", cloud_cover = c(0, 20), limit = 3)$features
names(x)
#> [1] "type"       "properties" "bbox"       "geometry"   "assets"     "links"
x$properties
#>                                         id         c:id                 datetime eo:cloud_cover eo:sun_azimuth
#> 1 LC08_L1TP_183023_20160625_20170323_01_T1 landsat-8-l1 2016-06-25T09:00:16.825Z              0      150.61964
#> 2 LC08_L1TP_183037_20160625_20170323_01_T1 landsat-8-l1 2016-06-25T09:05:51.253Z             19      110.95730
#> 3 LC08_L1TP_183041_20160625_20170323_01_T1 landsat-8-l1 2016-06-25T09:07:26.830Z              0       95.69133
#>   eo:sun_elevation landsat:path landsat:row
#> 1         57.71269          183          23
#> 2         68.16356          183          37
#> 3         68.55517          183          41
```

## List scenes


```r
(res <- lsat_scenes(n_max = 10))
#> # A tibble: 10 x 11
#>    entityId acquisitionDate     cloudCover processingLevel  path   row
#>    <chr>    <dttm>                   <dbl> <chr>           <int> <int>
#>  1 LC80101… 2015-01-02 15:49:05      80.8  L1GT               10   117
#>  2 LC80260… 2015-01-02 16:56:51      90.8  L1GT               26    39
#>  3 LC82270… 2015-01-02 13:53:02      83.4  L1GT              227    74
#>  4 LC82270… 2015-01-02 13:52:38      52.3  L1T               227    73
#>  5 LC82270… 2015-01-02 13:48:14      38.8  L1T               227    62
#>  6 LC82111… 2015-01-02 12:30:31      22.9  L1GT              211   115
#>  7 LC81791… 2015-01-02 09:14:45       7.67 L1GT              179   120
#>  8 LC82111… 2015-01-02 12:28:55      43.4  L1GT              211   111
#>  9 LC81950… 2015-01-02 10:17:20      21.0  L1T               195    29
#> 10 LC81790… 2015-01-02 08:44:49       1.92 L1T               179    45
#> # ... with 5 more variables: min_lat <dbl>, min_lon <dbl>, max_lat <dbl>,
#> #   max_lon <dbl>, download_url <chr>
```

## List scene files


```r
lsat_scene_files(x = res$download_url[1])
#>                                 file    size
#> 2   LC80101172015002LGN00_B4.TIF.ovr   7.7MB
#> 26 LC80101172015002LGN00_B11.TIF.ovr  17.0KB
#> 3       LC80101172015002LGN00_B5.TIF  56.8MB
#> 4      LC80101172015002LGN00_BQA.TIF   2.7MB
#> 5      LC80101172015002LGN00_MTL.txt   7.5KB
#> 6   LC80101172015002LGN00_B5.TIF.ovr   7.8MB
#> 7   LC80101172015002LGN00_B2.TIF.ovr   7.5MB
#> 8   LC80101172015002LGN00_B1.TIF.ovr   7.5MB
#> 9   LC80101172015002LGN00_B7.TIF.ovr   7.9MB
#> 10      LC80101172015002LGN00_B4.TIF  55.4MB
#> 11      LC80101172015002LGN00_B8.TIF 212.3MB
#> 12  LC80101172015002LGN00_B3.TIF.ovr   7.6MB
#> 13      LC80101172015002LGN00_B3.TIF  54.4MB
#> 14      LC80101172015002LGN00_B2.TIF  54.0MB
#> 15 LC80101172015002LGN00_B10.TIF.ovr  17.0KB
#> 16  LC80101172015002LGN00_B6.TIF.ovr   7.9MB
#> 17  LC80101172015002LGN00_B9.TIF.ovr   7.0MB
#> 18     LC80101172015002LGN00_B11.TIF   0.1MB
#> 19  LC80101172015002LGN00_B8.TIF.ovr  29.0MB
#> 20      LC80101172015002LGN00_B1.TIF  54.2MB
#> 21     LC80101172015002LGN00_B10.TIF   0.1MB
#> 22      LC80101172015002LGN00_B6.TIF  58.0MB
#> 23 LC80101172015002LGN00_BQA.TIF.ovr   0.6MB
#> 24      LC80101172015002LGN00_B7.TIF  58.0MB
#> 25      LC80101172015002LGN00_B9.TIF  49.6MB
```

## Get an image

Returns path to the image


```r
lsat_image(x = "LC80101172015002LGN00_B5.TIF")
#> [1] "/Users/sckott/Library/Caches/landsat-pds/L8/010/117/LC80101172015002LGN00/LC80101172015002LGN00_B5.TIF"
```

### Caching

When requesting an image, we first check if you already have that image. If you do,
we return the path to the file. If not, we get the image, and return the file path.


```r
lsat_image(x = "LC80101172015002LGN00_B5.TIF")
#> File in cache
#> [1] "/Users/sckott/Library/Caches/landsat-pds/L8/010/117/LC80101172015002LGN00/LC80101172015002LGN00_B5.TIF"
```

Note the message given.

See `?lsat_cache` for cache management functions.

## Visualize


```r
library("raster")
x <- lsat_cache_details()[[1]]
img <- raster(x$file)
plot(img)
```

![plot of chunk unnamed-chunk-10](inst/img/unnamed-chunk-10-1.png)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/getlandsat/issues).
* License: MIT
* Get citation information for `getlandsat` in R doing `citation(package = 'getlandsat')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
getlandsat 0.2.0
================

### MINOR IMPROVEMENTS

* moved to using markdown docs (#21)
* changed to using `crul` for HTTP requests (#20)


getlandsat 0.1.0
================

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 3.5.0
* ubuntu 12.04 (on travis-ci), R 3.5.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

License components with restrictions and base license permitting such:
   MIT + file LICENSE
 File 'LICENSE':
   YEAR: 2018
   COPYRIGHT HOLDER: Scott Chamberlain

## Reverse dependencies

There are no reverse dependencies.

---

This version updates to using markdown in documentation and swaps out a dependency.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/getlandsat/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/getlandsat.git`
* Make sure to track progress upstream (i.e., on our version of `getlandsat` at `ropensci/getlandsat`) by doing `git remote add upstream https://github.com/ropensci/getlandsat.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/getlandsat`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
getlandsat
==========

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.path = "inst/img/"
)
```

[![cran checks](https://cranchecks.info/badges/worst/getlandsat)](https://cranchecks.info/pkgs/getlandsat)
[![Build Status](https://travis-ci.org/ropensci/getlandsat.svg?branch=master)](https://travis-ci.org/ropensci/getlandsat)
[![codecov](https://codecov.io/gh/ropensci/getlandsat/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/getlandsat)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/getlandsat)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/getlandsat)](https://cran.r-project.org/package=getlandsat)
[![](https://badges.ropensci.org/58_status.svg)](https://github.com/ropensci/onboarding/issues/58)

`getlandsat`: Get Landsat 8 data from AWS public data sets

`getlandsat` provides access to Landsat <https://landsat.usgs.gov> 8 metadata and images hosted on AWS S3 at <https://registry.opendata.aws/landsat-8/>. The package only fetches data. It does not attempt to aid users in downstream usage, but some additional functionality may be added.

A new function `lsat_search()` lets you search for Landsat images by using the API from Development Seed documented at <https://github.com/sat-utils/sat-api>

Potential users are probably anyone from scientists asking questions about biodiversity or land use change, to software developers creating tools for users to vizualize their data.

## Install

Stable version

```{r eval=FALSE}
install.packages("getlandsat")
```

Dev version

```{r eval=FALSE}
devtools::install_github("ropensci/getlandsat")
```

```{r}
library("getlandsat")
```

## Search for images

```{r eval=FALSE}
x <- lsat_search(collection = "landsat-8", cloud_cover = c(0, 20), limit = 3)$features
names(x)
#> [1] "type"       "properties" "bbox"       "geometry"   "assets"     "links"
x$properties
#>                                         id         c:id                 datetime eo:cloud_cover eo:sun_azimuth
#> 1 LC08_L1TP_183023_20160625_20170323_01_T1 landsat-8-l1 2016-06-25T09:00:16.825Z              0      150.61964
#> 2 LC08_L1TP_183037_20160625_20170323_01_T1 landsat-8-l1 2016-06-25T09:05:51.253Z             19      110.95730
#> 3 LC08_L1TP_183041_20160625_20170323_01_T1 landsat-8-l1 2016-06-25T09:07:26.830Z              0       95.69133
#>   eo:sun_elevation landsat:path landsat:row
#> 1         57.71269          183          23
#> 2         68.16356          183          37
#> 3         68.55517          183          41
```

## List scenes

```{r}
(res <- lsat_scenes(n_max = 10))
```

## List scene files

```{r}
lsat_scene_files(x = res$download_url[1])
```

## Get an image

Returns path to the image

```{r}
lsat_image(x = "LC80101172015002LGN00_B5.TIF")
```

### Caching

When requesting an image, we first check if you already have that image. If you do,
we return the path to the file. If not, we get the image, and return the file path.

```{r message=TRUE}
lsat_image(x = "LC80101172015002LGN00_B5.TIF")
```

Note the message given.

See `?lsat_cache` for cache management functions.

## Visualize

```{r}
library("raster")
x <- lsat_cache_details()[[1]]
img <- raster(x$file)
plot(img)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/getlandsat/issues).
* License: MIT
* Get citation information for `getlandsat` in R doing `citation(package = 'getlandsat')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Introduction to getlandsat"
author: "Scott Chamberlain"
date: "2020-12-30"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{getlandsat introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



getlandsat introduction
=======================

`getlandsat` provides access to Landsat <https://landsat.usgs.gov> 8 metadata and images hosted on AWS S3 at <https://aws.amazon.com/public-data-sets/landsat>. The package only fetches data. It does not attempt to aid users in downstream usage, but some additional functionality may be added.

Potential users are probably anyone from scientists asking questions about biodiversity or land use change, to software developers creating tools for users to vizualize their data.

## Install

Dev version


```r
remotes::install_github("ropensci/getlandsat")
```


```r
library("getlandsat")
```


## List scenes


```r
(res <- lsat_scenes(n_max = 10))
#> # A tibble: 10 x 11
#>                 entityId     acquisitionDate cloudCover processingLevel
#>                    <chr>              <time>      <dbl>           <chr>
#> 1  LC80101172015002LGN00 2015-01-02 15:49:05      80.81            L1GT
#> 2  LC80260392015002LGN00 2015-01-02 16:56:51      90.84            L1GT
#> 3  LC82270742015002LGN00 2015-01-02 13:53:02      83.44            L1GT
#> 4  LC82270732015002LGN00 2015-01-02 13:52:38      52.29             L1T
#> 5  LC82270622015002LGN00 2015-01-02 13:48:14      38.85             L1T
#> 6  LC82111152015002LGN00 2015-01-02 12:30:31      22.93            L1GT
#> 7  LC81791202015002LGN00 2015-01-02 09:14:45       7.67            L1GT
#> 8  LC82111112015002LGN00 2015-01-02 12:28:55      43.43            L1GT
#> 9  LC81950292015002LGN00 2015-01-02 10:17:20      21.02             L1T
#> 10 LC81790452015002LGN00 2015-01-02 08:44:49       1.92             L1T
#> # ... with 7 more variables: path <int>, row <int>, min_lat <dbl>,
#> #   min_lon <dbl>, max_lat <dbl>, max_lon <dbl>, download_url <chr>
```

## List scene files


```r
lsat_scene_files(x = res$download_url[1])
#>                                 file    size
#> 1   LC80101172015002LGN00_B4.TIF.ovr   7.7MB
#> 2  LC80101172015002LGN00_B11.TIF.ovr  17.0KB
#> 3       LC80101172015002LGN00_B5.TIF  56.8MB
#> 4      LC80101172015002LGN00_BQA.TIF   2.7MB
#> 5      LC80101172015002LGN00_MTL.txt   7.5KB
#> 6   LC80101172015002LGN00_B5.TIF.ovr   7.8MB
#> 7   LC80101172015002LGN00_B2.TIF.ovr   7.5MB
#> 8   LC80101172015002LGN00_B1.TIF.ovr   7.5MB
#> 9   LC80101172015002LGN00_B7.TIF.ovr   7.9MB
#> 10      LC80101172015002LGN00_B4.TIF  55.4MB
#> 11      LC80101172015002LGN00_B8.TIF 212.3MB
#> 12  LC80101172015002LGN00_B3.TIF.ovr   7.6MB
#> 13      LC80101172015002LGN00_B3.TIF  54.4MB
#> 14      LC80101172015002LGN00_B2.TIF  54.0MB
#> 15 LC80101172015002LGN00_B10.TIF.ovr  17.0KB
#> 16  LC80101172015002LGN00_B6.TIF.ovr   7.9MB
#> 17  LC80101172015002LGN00_B9.TIF.ovr   7.0MB
#> 18     LC80101172015002LGN00_B11.TIF   0.1MB
#> 19  LC80101172015002LGN00_B8.TIF.ovr  29.0MB
#> 20      LC80101172015002LGN00_B1.TIF  54.2MB
#> 21     LC80101172015002LGN00_B10.TIF   0.1MB
#> 22      LC80101172015002LGN00_B6.TIF  58.0MB
#> 23 LC80101172015002LGN00_BQA.TIF.ovr   0.6MB
#> 24      LC80101172015002LGN00_B7.TIF  58.0MB
#> 25      LC80101172015002LGN00_B9.TIF  49.6MB
```

## Get an image

Returns path to the image


```r
lsat_image(x = "LC80101172015002LGN00_B5.TIF")
#> [1] "/Users/sacmac/Library/Caches/landsat-pds/L8/010/117/LC80101172015002LGN00/LC80101172015002LGN00_B5.TIF"
```

Another one


```r
lsat_image("LC80010032014272LGN00_B10.TIF")
#> [1] "/Users/sacmac/Library/Caches/landsat-pds/L8/001/003/LC80010032014272LGN00/LC80010032014272LGN00_B10.TIF"
```

### Caching

When requesting an image, we first check if you already have that image. If you do, 
we return the path to the file. If not, we get the image, and return the file path.


```r
lsat_image(x = "LC80101172015002LGN00_B5.TIF")
#> File in cache
#> [1] "/Users/sacmac/Library/Caches/landsat-pds/L8/010/117/LC80101172015002LGN00/LC80101172015002LGN00_B5.TIF"
```

Note the message given.

See `?lsat_cache` for cache management functions.

## Visualize


```r
library("raster")
x <- lsat_cache_details()
img <- raster(x[[1]]$file)
plot(img)
```

![plot of chunk unnamed-chunk-10](img/unnamed-chunk-10-1.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getlandsat-package.R
\docType{package}
\name{getlandsat-package}
\alias{getlandsat-package}
\alias{getlandsat}
\title{getlandsat - get Landsat 8 data from AWS public data sets}
\description{
\pkg{getlandsat} provides access to Landsat \url{https://landsat.usgs.gov} 8
metadata and images hosted on AWS S3 at
\url{https://registry.opendata.aws/landsat-8/}. The package only
fetches data. It does not attempt to aid users in downstream usage.
}
\examples{
\dontrun{
## List scenes
(res <- lsat_scenes(n_max = 10))

## List scene files
lsat_scene_files(x = res$download_url[1])

## Get an image
### Returns path to the image
lsat_image(x = "LC80101172015002LGN00_B5.TIF")

## Visualize
if (requireNamespace("raster")) {
  library("raster")
  x <- lsat_cache_details()[[1]]
  img <- raster(x$file)
  plot(img)
}
}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsat_scenes.R
\name{lsat_scenes}
\alias{lsat_scenes}
\title{List Landsat scenes}
\usage{
lsat_scenes(...)
}
\arguments{
\item{...}{Further args passed on to \code{\link[readr:read_delim]{readr::read_csv()}}}
}
\description{
List Landsat scenes
}
\details{
We use \code{\link[readr:read_delim]{readr::read_csv()}} to read the scene list file
from \url{http://landsat-pds.s3.amazonaws.com/scene_list.gz}. See the help
file for \code{\link[readr:read_delim]{readr::read_csv()}} for what parameter you can pass to modify it's
behavior.

This is an alternative to using \code{\link[=lsat_list]{lsat_list()}}. This function
downloads the up to date compressed csv file, while \code{\link[=lsat_list]{lsat_list()}}
uses the AWS S3 API.
}
\examples{
\dontrun{
res <- lsat_scenes()
res

# read only N rows
lsat_scenes(n_max = 10)
}
}
\seealso{
\code{\link[=lsat_scene_files]{lsat_scene_files()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{lsat_cache}
\alias{lsat_cache}
\alias{lsat_cache_list}
\alias{lsat_cache_delete}
\alias{lsat_cache_delete_all}
\alias{lsat_cache_details}
\title{Manage cached files}
\usage{
lsat_cache_list()

lsat_cache_delete(files, force = TRUE)

lsat_cache_delete_all(force = TRUE)

lsat_cache_details(files = NULL)
}
\arguments{
\item{files}{(character) one or more complete file names}

\item{force}{(logical) Should files be force deleted? Default: \code{TRUE}}
}
\description{
Manage cached files
}
\details{
\code{cache_delete} only accepts 1 file name, while \code{cache_delete_all}
doesn't accept any names, but deletes all files. For deleting many
specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}} type call

We cache using \code{\link[rappdirs:user_cache_dir]{rappdirs::user_cache_dir()}}, find your cache
folder by executing \code{rappdirs::user_cache_dir("landsat-pds")}
}
\section{Functions}{

\itemize{
\item \code{lsat_cache_list()} returns a character vector of full path
file names
\item \code{lsat_cache_delete()} deletes one or more files, returns nothing
\item \code{lsat_cache_delete_all()} delete all files, returns nothing
\item \code{lsat_cache_details()} prints file name and file size for each file,
supply with one or more files, or no files (and get details for
all available)
}
}

\examples{
\dontrun{
# list files in cache
lsat_cache_list()

# List info for single files
lsat_cache_details(files = lsat_cache_list()[1])
lsat_cache_details(files = lsat_cache_list()[2])

# List info for all files
lsat_cache_details()

# delete files by name in cache
# lsat_cache_delete(files = lsat_cache_list()[1])

# delete all files in cache
# lsat_cache_delete_all()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsat_image.R
\name{lsat_image}
\alias{lsat_image}
\title{Get Landsat image(s)}
\usage{
lsat_image(x, overwrite = FALSE, ...)
}
\arguments{
\item{x}{(character) A file name for a geotif file, will be more general soon.}

\item{overwrite}{(logical) Will only overwrite existing path if \code{TRUE}.
Deprecated, will be removed in the next version. If file exists we return
that path so there's no chance of overwriting}

\item{...}{Curl args passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
Path to the file, whether found in cache or new file
requested.
}
\description{
Get Landsat image(s)
}
\examples{
\dontrun{
# pass an image name
(res <- lsat_list(max = 40))
tifs <- grep("\\\\.TIF$", res$Key, value = TRUE)
lsat_image(tifs[5])
lsat_image(tifs[6])
lsat_image(tifs[9])

# caching
## requesting an image you already have will return path if found
lsat_image(tifs[5])
}
}
\seealso{
\code{\link[=lsat_cache]{lsat_cache()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsat_list.R
\name{lsat_list}
\alias{lsat_list}
\title{List Landsat images}
\usage{
lsat_list(max = NULL, marker = NULL, prefix = NULL, delimiter = NULL, ...)
}
\arguments{
\item{max}{(integer) number indicating the maximum number of keys to
return (max 1000, default 1000).}

\item{marker}{(character) string that pecifies the key to start with
when listing objects in a AWS bucket. Amazon S3 returns object keys in
alphabetical order, starting with key after the marker in order}

\item{prefix}{(character) string that limits the response to keys
that begin with the specified prefix}

\item{delimiter}{(character) string used to group keys. Read the AWS
doc for more detail.}

\item{...}{curl args passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
List Landsat images
}
\details{
This is an alternative to using \code{\link[=lsat_scenes]{lsat_scenes()}}. This
function uses the AWS S3 API, while \code{\link[=lsat_scenes]{lsat_scenes()}} downloads
the up to date compressed csv file.
}
\examples{
\dontrun{
lsat_list(max = 10)

# paging, start a specific key string
res <- lsat_list(max = 10)
lsat_list(marker = res$Key[10], max = 10)

# curl options
lsat_list(max = 3, verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsat_scenes_files.R
\name{lsat_scene_files}
\alias{lsat_scene_files}
\title{List files for a Landsat scene}
\usage{
lsat_scene_files(x, ...)
}
\arguments{
\item{x}{(character) A URL to a scene html file}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A data.frame with two columns:
\itemize{
\item file - file name
\item size - file size
}
}
\description{
List files for a Landsat scene
}
\details{
This function fetches files available in a scene, while
\code{\link[=lsat_scenes]{lsat_scenes()}} lists the scenes, but not their files
}
\examples{
\dontrun{
res <- lsat_scenes(n_max = 10)
lsat_scene_files(x = res$download_url[1])
lsat_scene_files(x = res$download_url[2])
}
}
\seealso{
\code{\link[=lsat_scenes]{lsat_scenes()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsat_search.R
\name{lsat_search}
\alias{lsat_search}
\title{Search Landsat images via the sat-api}
\usage{
lsat_search(
  collection = NULL,
  datetime = NULL,
  cloud_cover = NULL,
  sun_azimuth = NULL,
  sun_elevation = NULL,
  row = NULL,
  path = NULL,
  limit = 10,
  page = 1,
  as = "df",
  ...
)
}
\arguments{
\item{collection}{(character) collection name}

\item{datetime}{(character) one or two dates. if one date, it's treated
as a single point date time. if two dates, its treated as a date
period}

\item{cloud_cover}{(integer) cloud cover number}

\item{sun_azimuth}{(numeric) sun azimuth number}

\item{sun_elevation}{(numeric) sun elevation number}

\item{row}{(integer) row}

\item{path}{(integer) path}

\item{limit}{(integer) limit number of results. default: 10}

\item{page}{(integer) page number to return. default: 1}

\item{as}{(character) if \code{df} we parse from JSON
to an R list with data.frame's where possible. if 'list' then an R list;
if \code{json} we give you back the JSON as character}

\item{...}{curl args passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
an R list with data.fraem's if \code{as=df}; an R list w/o attempting to
parse to data.frame's where possible; otherwise JSON as character
}
\description{
Search Landsat images via the sat-api
}
\examples{
\dontrun{
# collection input
lsat_search(collection = "sentinel-2")
lsat_search(collection = "landsat-8")

# dates
lsat_search(datetime = "2017-08")
z <- lsat_search(datetime = c("2016-12-31", "2017-01-01"))
z

# cloud cover
lsat_search(collection = "landsat-8", cloud_cover = 0)
lsat_search(collection = "landsat-8", cloud_cover = c(0, 20))

# sun azimuth
lsat_search(sun_azimuth = 162)

# row/path
lsat_search(row = 31, path = 128)

# pagination
lsat_search(limit = 0)
lsat_search(limit = 1)
lsat_search(limit = 3, page = 2)

# parsed to list, or not gives json
lsat_search(datetime = "2017-08", limit = 3, as = 'df')
lsat_search(datetime = "2017-08", limit = 3, as = 'list')
lsat_search(datetime = "2017-08", limit = 3, as = 'json')

# curl options
lsat_search(verbose = TRUE)
}
}
\references{
\url{https://github.com/sat-utils/sat-api}
}
