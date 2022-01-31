rerddap
=====



[![cran checks](https://cranchecks.info/badges/worst/rerddap)](https://cranchecks.info/pkgs/rerddap)
[![R-check](https://github.com/ropensci/rerddap/workflows/R-check/badge.svg)](https://github.com/ropensci/rerddap/actions)
[![codecov.io](https://codecov.io/github/ropensci/rerddap/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rerddap?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rerddap)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rerddap)](https://cran.r-project.org/package=rerddap)

`rerddap` is a general purpose R client for working with ERDDAP servers.

Package Docs: <https://docs.ropensci.org/rerddap/>

## Installation

From CRAN


```r
install.packages("rerddap")
```

Or development version from GitHub


```r
remotes::install_github("ropensci/rerddap")
```


```r
library("rerddap")
```

Some users may experience an installation error, stating to install 1 or more 
packages, e.g., you may need `DBI`, in which case do, for example, 
`install.packages("DBI")` before installing `rerddap`.

## Background

ERDDAP is a server built on top of OPenDAP, which serves some NOAA data. You can get gridded data (griddap (<https://upwell.pfeg.noaa.gov/erddap/griddap/documentation.html>)), which lets you query from gridded datasets, or table data (tabledap (<https://upwell.pfeg.noaa.gov/erddap/tabledap/documentation.html>)) which lets you query from tabular datasets. In terms of how we interface with them, there are similarities, but some differences too. We try to make a similar interface to both data types in `rerddap`.

## NetCDF

`rerddap` supports NetCDF format, and is the default when using the `griddap()` function. NetCDF is a binary file format, and will have a much smaller footprint on your disk than csv. The binary file format means it's harder to inspect, but the `ncdf4` package makes it easy to pull data out and write data back into a NetCDF file. Note the the file extension for NetCDF files is `.nc`. Whether you choose NetCDF or csv for small files won't make much of a difference, but will with large files.

## Caching

Data files downloaded are cached in a single directory on your machine determined by the `hoardr` package. When you use `griddap()` or `tabledap()` functions, we construct a MD5 hash from the base URL, and any query parameters - this way each query is separately cached. Once we have the hash, we look in the cache directory for a matching hash. If there's a match we use that file on disk - if no match, we make a http request for the data to the ERDDAP server you specify.

## ERDDAP servers

You can get a data.frame of ERDDAP servers using the function `servers()`. Most I think serve some kind of NOAA data, but there are a few that aren't NOAA data.  If you know of more ERDDAP servers, send a pull request, or let us know.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rerddap/issues).
* License: MIT
* Get citation information for `rerddap` in R doing `citation(package = 'rerddap')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rerddap 0.8.0
=============

* Added global search function
* fixed bug when dataset has a decreasing coordinate that
is not latitude or longitude


rerddap 0.7.6
=============

### MINOR IMPROVEMENTS

* fixed a bug in dealing with trailing slashes in URLs

rerddap 0.7.4
=============

### MINOR IMPROVEMENTS

* fix a broken test

rerddap 0.7.0
=============

### MINOR IMPROVEMENTS

* vignettes only on package documentation site now  (#87)
* `server()` (to fetch known ERDDAP server URLs) now uses the list maintained by `irishmarineinstitute/awesome-erddap` on GitHub (#86)
* better error handling for `griddap()`: if no dimension arguments passed, we error saying so (and no http requests made); in addition, if a dataset is passed to `griddap()`, to which the output of `info()` was also passed, then we can check if the dataset has griddap data or not, and fail saying so if not (#91)
* `griddap()` and `tabledap()`: if `info()` output passed to these two funcitons, we will now use the url within that info output, and use a message telling the user we are doing so; now you don't have to set the url if you pass info output  (#92)


rerddap 0.6.5
=============

### BUG FIXES

* fix a `convert_units` test that was failing because remote service had changed the response


rerddap 0.6.4
=============

### BUG FIXES

* fix to internal fxn `err_handle()` for handling http errors - ERDDAP servers changed to some weird JSON-ish type format (#85)


rerddap 0.6.0
=============

### MINOR IMPROVEMENTS

* change all `tibble::as_data_frame`/`tibble::data_frame` to `tibble::as_tibble` (#79)
* `info()` gains new element in its output list, `base_url`, the base url for the ERDDAP server under consideration (#80)
* improved docs for `griddap()` with respect to what's returned from the function  (#81)
* fix some test fixtures to use preserve exact bytes so that cran checks on debian clang devel don't fail (#83)
* add .github files: contributing, issue template, pull request template

### BUG FIXES

* fix for lat/lon parsing within `griddap()` to account for cases when min and max are reversed from the order they should be in (#78)
* fix to `griddap()` to parse additioanl dimensions returned; previously we were only returning time, lat, and lon, plus one more (#82) thanks @afredstonhermann


rerddap 0.5.0
=============

### MINOR IMPROVEMENTS

* added new `Caching` section to package level manual file (`?rerddap`) about caching  (#52)
* use markdown docs in package (#75)
* replace `httr` with `crul` (#54)
* cache most tests with HTTP requests using `vcr` (#76)
* add test for `read` parameter in `griddap()` (#47)
* use default url via `eurl()`; used as default in main functions; set default url with env vars, see `?eurl`  (#41)
* improve handling and reporting back to user of ERDDAP server errors (#70) (#73)
* change to `griddap()`: when nc format gridded datasets have latitude and longitude we "melt" them into a data.frame for easy downstream consumption. When nc format gridded datasets do not have latitude and longitude components, we do not read in the data, throw a warning saying so. You can readin the nc file yourself with the file path (#74)
* for for `griddap()` to support cases in wihch lat/lon runs north to south and south to north (#68)

### BUG FIXES

* `memory()` usage in `griddap()` wasn't working. fixed now (#77)


rerddap 0.4.2
=============

### NEW FEATURES

* Now using `hoardr` to manage caching paths and such (#60). Also
now asking users where they want to cache files, either in a 
`rappdirs` user cache dir or a temp directory. Now on tests and examples
we use temp dirs.
* Related to above, new functions `cache_info()` to get cache path and 
number of cached files, and `cache_setup()` to set cache path.
* Related to above, `cache_details()`, `cache_list()`, and `cache_delete()`
lose their `cache_path` parameter - now cache path is set package wide and 
we use the same cache path, so no need to set in the fxn call.

### MINOR IMPROVEMENTS

* Fixes to a number of `griddap()` and `tabledap()` examples to use 
datasets that still exist (previous examples used datasets that are no
gone)


rerddap 0.4.0
=============

### NEW FEATURES

* New vignette added that goes in to much more depth than 
the original vignette (#51) thx to @rmendels
* `info()` function gains new attribute `url` with the 
base url for the ERDDAP server used (#42)
* Replaced usage of internal compact data.frame code to 
use `tibble` package (#45)

### MINOR IMPROVEMENTS

* Added another ERDDAP server to `servers()` function (#49)
* Changed base URLs for default ERDDAP server from `http` 
to `https`  (#50)
* Added note to docs for `griddap()` and `tabledap()` for how
to best deal with 500 server errors (#48)
* Replaced all `dplyr::rbind_all` uses with `dplyr::bind_rows` (#46)


rerddap 0.3.4
=============

### MINOR IMPROVEMENTS

* Removed use of `ncdf` package, which has been taken off CRAN.
Using `ncdf4` now for all NetCDF file manipulation. (#35)
* Failing better now with custom error catching (#31)
* Added many internal checks for parameter inputs, warning or
stopping as necessary - ERDDAP servers silently drop with no
informative messages (#32)

### BUG FIXES

* Using now `file.info()$size` instead of `file.size()` to be
backwards compatible with R versions < 3.2


rerddap 0.3.0
=============

### NEW FEATURES

* Cache functions accept the outputs of `griddap()` and `tabledap()`
so that the user can easily see cache details or delete the file from
the cache without having to manually get the file name. (#30)

### MINOR IMPROVEMENTS

* All package dependencies now use `importFrom` so we only import
functions we need instead of their global namespaces.

### BUG FIXES

* Fixed bug in parsing data from netcdf files, affected the
`griddap()` function (#28)


rerddap 0.2.0
=============

### NEW FEATURES

* Added a suite of functions to manage local cached files (#17)

### MINOR IMPROVEMENTS

* Added new ERDDAP server to list of servers in the `servers()` function (#21)

### BUG FIXES

* Fixed a few cases across a number of functions in which an empty list
passed to `query` parmaeter in `httr::GET()` caused an error (#23)
* Fixed retrieval of path to file written to disk by `httr::write_disk()` (#24)
* `last` is a value accepted by ERDDAP servers, but internal functions
weren't checking correctly, fixed now. (#25)
* `as.info()` wasn't passing on the `url` parameter to the `info()` function.
fixed now. (#26)


rerddap 0.1.0
=============

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.1.1
* rhub ubuntu, debian, solaris, Apple Silicon (M1), macOS 11.6 Big Sur
* win-builder (devel and release)

## R CMD check results

OK from all checks

## Reverse dependencies

* No probelms with reverse dependencies.

---

This version fixes a bug dealing with coordinates that
are not lat-lon and are in decreasing order.  Adds a new
search function.

Thanks! 
Roy Mendelssohn
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

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

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  CORRECT: you edit a roxygen comment in a `.R` file below `R/`.
*  INCORRECT: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes. The `README` should contain badges for any continuous integration services used by the package.
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the `rerddadp` project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html) for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) DO NOT SHARE SCREENSHOTS OF CODE; SHARE ACTUAL CODE IN TEXT FORMAT -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.4 Patched (2021-02-17 r80031) |
|os       |macOS Big Sur 10.16                         |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-03-05                                  |

# Dependencies

|package   |old   |new      |Δ  |
|:---------|:-----|:--------|:--|
|rerddap   |0.7.0 |0.7.1.92 |*  |
|cli       |NA    |2.3.1    |*  |
|crul      |NA    |1.1.0    |*  |
|dplyr     |NA    |1.0.4    |*  |
|lifecycle |NA    |1.0.0    |*  |
|pillar    |NA    |1.5.1    |*  |
|tibble    |NA    |3.1.0    |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*rerddap
=====

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![cran checks](https://cranchecks.info/badges/worst/rerddap)](https://cranchecks.info/pkgs/rerddap)
[![R-check](https://github.com/ropensci/rerddap/workflows/R-check/badge.svg)](https://github.com/ropensci/rerddap/actions)
[![codecov.io](https://codecov.io/github/ropensci/rerddap/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rerddap?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rerddap)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rerddap)](https://cran.r-project.org/package=rerddap)

`rerddap` is a general purpose R client for working with ERDDAP servers.

Package Docs: <https://docs.ropensci.org/rerddap/>

## Installation

From CRAN

```{r eval=FALSE}
install.packages("rerddap")
```

Or development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/rerddap")
```

```{r}
library("rerddap")
```

Some users may experience an installation error, stating to install 1 or more 
packages, e.g., you may need `DBI`, in which case do, for example, 
`install.packages("DBI")` before installing `rerddap`.

## Background

ERDDAP is a server built on top of OPenDAP, which serves some NOAA data. You can get gridded data (griddap (<https://upwell.pfeg.noaa.gov/erddap/griddap/documentation.html>)), which lets you query from gridded datasets, or table data (tabledap (<https://upwell.pfeg.noaa.gov/erddap/tabledap/documentation.html>)) which lets you query from tabular datasets. In terms of how we interface with them, there are similarities, but some differences too. We try to make a similar interface to both data types in `rerddap`.

## NetCDF

`rerddap` supports NetCDF format, and is the default when using the `griddap()` function. NetCDF is a binary file format, and will have a much smaller footprint on your disk than csv. The binary file format means it's harder to inspect, but the `ncdf4` package makes it easy to pull data out and write data back into a NetCDF file. Note the the file extension for NetCDF files is `.nc`. Whether you choose NetCDF or csv for small files won't make much of a difference, but will with large files.

## Caching

Data files downloaded are cached in a single directory on your machine determined by the `hoardr` package. When you use `griddap()` or `tabledap()` functions, we construct a MD5 hash from the base URL, and any query parameters - this way each query is separately cached. Once we have the hash, we look in the cache directory for a matching hash. If there's a match we use that file on disk - if no match, we make a http request for the data to the ERDDAP server you specify.

## ERDDAP servers

You can get a data.frame of ERDDAP servers using the function `servers()`. Most I think serve some kind of NOAA data, but there are a few that aren't NOAA data.  If you know of more ERDDAP servers, send a pull request, or let us know.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rerddap/issues).
* License: MIT
* Get citation information for `rerddap` in R doing `citation(package = 'rerddap')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Using rerddap to Access Data from ERDDAP Servers"
author: "Roy Mendelssohn and Scott Chamberlain"
date: "2021-11-18"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{usimg rerddap}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---




## Introduction

`rerddap` is a general purpose <span style="color:red">R</span> client for working with <span style="color:red">ERDDAP</span> servers.  <span style="color:red">ERDDAP</span> is a web service developed by Bob Simons of NOAA.  At the time of this writing, there are over sixty <span style="color:red">ERDDAP</span> servers  (though not all are public facing) providing access to literally petabytes of data and model output relevant to oceanography, meteorology, fisheries and marine mammals, among other areas. <span style="color:red">ERDDAP</span> is a simple to use, RESTful web service, that allows data to be subsetted and returned in a variety of formats.

This vignette goes over some of the nuts and bolts of using the `rerddap` package, and shows the power of the combination of the `rerddap` package with <span style="color:red">ERDDAP</span> servers.  Some of the examples are taken from the `xtractomatic` package (available from CRAN - https://cran.r-project.org/package=xtractomatic), and some from the `rerddapXtracto` package available on Github (https://github.com/rmendels/rerddapXtracto),  but reworked to use `rerddap` directly.   Other examples are new to this vignette, and include both gridded and non-gridded datasets from several <span style="color:red">ERDDAPs</span>.

## Installation

The first step is to install the `rerddap` package, the stable version is available from CRAN:



```r
install.packages("rerddap")
```

or the development version can be installed from GitHub:


```r
remotes::install_github("ropensci/rerddap")
```

and to load the library:


```r
library("rerddap")
```

Besides `rerddap` the following libraries are used in this vignette:





```r
library("akima")
library("dplyr")
library("ggplot2")
library("mapdata")
library("ncdf4")
library("plot3D")
```

Code chunks are always given with the required libraries so that the chunks are more standalone in nature. Many of the plots use an early version of the `cmocean` colormaps designed by Kristen Thyng (see https://matplotlib.org/cmocean/ and https://github.com/matplotlib/cmocean) translatd from the original Python implementation. However, there is now a `cmocean` package for <span style="color:red">R</span>.  As past scripts may use the built-in color palette it is maintained here,  but use of the `cmocean` package is advisable as it is kept up-to-date.  In the present version of the `cmocean` package the names have changed, what is "temperature" here is now "thermal", "chlorophyll" is now "algae", and "salinity" is now "haline".

## The main `rerddap` functions

The complete list of rerddap functions can be seen by looking at he `rerddap` package help:


```r
?rerddap
```

and selecting the index of the package.  The main functions used here are:

* the list of servers `rerddap` knows about - `servers()`
* search an <span style="color:red">ERDDAP</span> server for terms - `ed_search(query, page = NULL, page_size = NULL, which = "griddap", url = eurl(), ...)`
* search a list of  <span style="color:red">ERDDAP</span> servers for terms - `global_search(query, server_list, which_service)`
* get a list of datasets on an <span style="color:red">ERDDAP</span> server - `ed_datasets(which = "tabledap", url = eurl())`
* obtain information about a dataset - `info(datasetid, url = eurl(), ...)`
* extract data from a griddap dataset - `griddap(x, ..., fields = "all", stride = 1, fmt = "nc", url = eurl(), store = disk(), read = TRUE, callopts = list())`
* extract data from a tabledap dataset - `tabledap(x, ..., fields = NULL, distinct = FALSE, orderby = NULL, orderbymax = NULL, orderbymin = NULL, orderbyminmax = NULL, units = NULL, url = eurl(), store = disk(), callopts = list())`

Be careful when using the functions `ed_search()`, `ed_datasets()` and `global_search()`.  The default <span style="color:red">ERDDAP</span> has over 9,000 datasets,  most of which are grids, so that a list of all the gridded datasets can be quite long.  A seemly reasonable search:


```r
whichSST <- ed_search(query = "SST")
```

returns about 1000 responses.  The more focused query:


```r
whichSST <- ed_search(query = "SST MODIS")
```

still returns 172 responses.  If the simple search doesn't narrow things enough,  look at the advanced search function `ed_search_adv()`.

## Finding the Data You Want

The first way to find a dataset is to browse the builtin web page for a particular <span style="color:red">ERDDAP</span> server.
A list of some of the public available <span style="color:red">ERDDAP</span> servers can be obtained from the `rerddap` command:


```r
servers()
#> # A tibble: 58 × 4
#>    name                                 short_name  url                   public
#>    <chr>                                <chr>       <chr>                 <lgl> 
#>  1 CoastWatch West Coast Node           CSWC        https://coastwatch.p… TRUE  
#>  2 ERDDAP at the Asia-Pacific Data-Res… APDRC       https://apdrc.soest.… TRUE  
#>  3 NOAA's National Centers for Environ… NCEI        https://www.ncei.noa… TRUE  
#>  4 Biological and Chemical Oceanograph… BCODMO      https://erddap.bco-d… TRUE  
#>  5 European Marine Observation and Dat… EMODnet     https://erddap.emodn… TRUE  
#>  6 Marine Institute - Ireland           MII         https://erddap.marin… TRUE  
#>  7 CoastWatch Caribbean/Gulf of Mexico… CSCGOM      https://cwcgom.aoml.… TRUE  
#>  8 NOAA IOOS Sensors ERDDAP             IOOS-Senso… https://erddap.senso… TRUE  
#>  9 NOAA IOOS CeNCOOS (Central and Nort… CeNCOOS     https://erddap.axiom… TRUE  
#> 10 CeNCOOS (Central and Northern Calif… CeNCOOS ER… http://erddap.cencoo… TRUE  
#> # … with 48 more rows
```

The list of <span style="color:red">ERDDAP</span> servers is based on the list maintained at the [Awesome ERDDAP site](https://irishmarineinstitute.github.io/awesome-erddap).

The second way to find and obtain the desired data is to use functions in `rerddap`.  The basic steps are:

1. Find the dataset on an <span style="color:red">ERDDAP</span> server (`rerddap::servers()`, `rerddap::ed_search()`, `rerddap::ed_datasets()` ).
2. Get the needed information about the dataset (`rerddap::info()` )
3. Think about what you are going to do.
4. Make the request for the data  (`rerddap::griddap()` or `rerddap::tabledap()` ).

We discuss each of these steps in more detail, and then look at some realistic uses of the package.


### Think about what you are going to do.

This may seem to be a strange step in the process, but it is important because many of the datasets are high-resolution, and data requests can get very large, very quickly.  As an example, based on a real use case.  The MUR SST (	Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1, see https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html )  is a daily, high-quality, high-resolution sea surface temperature product.  The user wanted the MUR data for a 2x2-degree box, daily for a year.  That seems innocuous enough.  Except that MURsst is at a one-hundreth of degree resolution.  If we assume just a binary representation of the data, assuming 8-bytes per value, and do the math:


```r
100*100*4*8*365
#> [1] 116800000
```


Yes, 116,800,000 bytes or roughly 115MB for that request.  Morever the user wanted the data as a .csv file, which usually makes the resulting file 8-10 times larger,  so now we are over a 1GB for the request. Even more so, there are four parameters in that dataset, and in `rerddap::griddap()` if "fields" is not specified,  all fields are downloaded, therefore the resulting files will be four times as large as given above.

So the gist of this is to think about your request before you make it.  Do a little mental math to get a rough estimate of the size of the download.  There are times receiving the data as a .csv file is convenient,  but make certain the request will not be too large.  For larger requests, obtain the data as netCDF files.  By default, `rerddap::griddap()` "melts"" the data into a dataframe,  so a .csv only provides a small convenience.  But for really large downloads,  you should select the option in `rerddap::griddap()` to not read in the data, and use instead the `netcdf4` package to read in the data, as this allows for only reading in parts of the data at a time.  [Below](#ncdf4) we provide a brief tutorial on reading in data using the `ncdf4` package.

## Some ERDDAP Basics

One of the main advantages of a service such as <span style="color:red">ERDDAP</span> is that you only need to download the subset of the data you desire,  rather than the entire dataset,  which is both convenient and essential for large datasets. The underlying data model in <span style="color:red">ERDDAP</span> is quite simple - everything is either a (multi-dimensional) grid  (think <span style="color:red">R</span> array) or a table (think a simple spreadsheet or data table). Grids are subsetted using the function `griddap()` and tables are subset using the function `tabledap()`.

If you know the datasetID of the data you are after, and are unsure if it is a grid or a table,  there are several ways to find out. We will look at two datasets, 'jplMURSST41' and 'siocalcofiHydroCasts'.  The first method is to use the `rerddap` function `browse()`


```r
browse('jplMURSST41')
browse('siocalcofiHydroCasts')
```

which brings up information on the datasets in the browser, in the first case the "data" link is under "griddap", the second is under "tabledap".

The other method is to use the `rerddap` function `info`:


```r
info('jplMURSST41')
#> <ERDDAP info> jplMURSST41 
#>  Base URL: https://upwell.pfeg.noaa.gov/erddap 
#>  Dataset Type: griddap 
#>  Dimensions (range):  
#>      time: (2002-06-01T09:00:00Z, 2021-11-17T09:00:00Z) 
#>      latitude: (-89.99, 89.99) 
#>      longitude: (-179.99, 180.0) 
#>  Variables:  
#>      analysed_sst: 
#>          Units: degree_C 
#>      analysis_error: 
#>          Units: degree_C 
#>      mask: 
#>      sea_ice_fraction: 
#>          Units: 1
info('siocalcofiHydroCasts')
#> <ERDDAP info> siocalcofiHydroCasts 
#>  Base URL: https://upwell.pfeg.noaa.gov/erddap 
#>  Dataset Type: tabledap 
#>  Variables:  
#>      ac_line: 
#>      ac_sta: 
#>      barometer: 
#>          Units: millibars 
#>      bottom_d: 
#>          Units: meters 
#>      cast_id: 
#>      civil_t: 
#>          Units: seconds since 1970-01-01T00:00:00Z 
#>      cloud_amt: 
#>          Units: oktas 
#>      cloud_typ: 
#>      cruise_id: 
#>      cruz_leg: 
#>      cruz_num: 
#>      cruz_sta: 
#>      cst_cnt: 
#>      data_or: 
#>      data_type: 
#>      date: 
#>      dbsta_id: 
#>      distance: 
#>      dry_t: 
#>          Units: degree C 
#>      event_num: 
#>      forelu: 
#>          Units: Forel-Ule scale 
#>      inc_end: 
#>          Units: seconds since 1970-01-01T00:00:00Z 
#>      inc_str: 
#>          Units: seconds since 1970-01-01T00:00:00Z 
#>      intc14: 
#>          Units: milligrams Carbon per square meter 
#>      intchl: 
#>      julian_date: 
#>      julian_day: 
#>      latitude: 
#>          Range: 18.417, 47.917 
#>          Units: degrees_north 
#>      latitude_degrees: 
#>      latitude_hemisphere: 
#>      latitude_minutes: 
#>      longitude: 
#>          Range: -164.083, -105.967 
#>          Units: degrees_east 
#>      longitude_degrees: 
#>      longitude_hemisphere: 
#>      longitude_minutes: 
#>      month: 
#>          Range: 1, 12 
#>      order_occ: 
#>      orig_sta_id: 
#>      pst_lan: 
#>          Units: seconds since 1970-01-01T00:00:00Z 
#>      quarter: 
#>      rpt_line: 
#>      rpt_sta: 
#>      secchi: 
#>          Units: meters 
#>      ship_code: 
#>      ship_name: 
#>      st_line: 
#>      st_station: 
#>      sta_code: 
#>      sta_id: 
#>      time: 
#>          Range: -6.5759508E8, 1.423140778E9 
#>          Units: seconds since 1970-01-01T00:00:00Z 
#>      time_ascii: 
#>      timezone: 
#>      visibility: 
#>      wave_dir: 
#>          Units: degrees 
#>      wave_ht: 
#>          Units: feet 
#>      wave_prd: 
#>          Units: seconds 
#>      wea: 
#>      wet_t: 
#>          Units: degree C 
#>      wind_dir: 
#>          Units: degrees 
#>      wind_spd: 
#>          Units: knots 
#>      year: 
#>          Range: 1949, 2015
```

Notice the information on 'jplMURSST41' lists the dimensions (that is a grid) while that of 'siocalcofiHydroCasts' does not (that is a table).

### Subsetting griddap()

Like an <span style="color:red">R</span> array, <span style="color:red">ERDDAP</span> grids are subsetted by setting limits on the dimension variables, the difference being that a subset is defined in coordinate space  (latitude values, longitude values, time values) rather than array space as is done with <span style="color:red">R</span> arrays.  Thus for 'jplMURSST41' if the desired area of the data extract is latitude limits of (22N, 51N), longitude limits of (140W, 105W), and with time limits of (2017-01-01, 2017-01-02) the following would be passed to the function `griddap()`:


```r
latitude = c(22., 51.)
longitude = c(-140., -105)
time = c("2017-01-01", "2017-01-02")
```

A full `griddap()` request to retrieve "analysed_sst" with these constraints would be:


```r
sstInfo <- info('jplMURSST41')
murSST <- griddap(sstInfo, latitude = c(22., 51.), longitude = c(-140., -105), time = c("2017-01-01", "2017-01-02"), fields = 'analysed_sst')

```

#### Strides

Strides allow to retrieve data within a coordinate bound only at every "n" values,  where "n" is an integer - think of the "by" part of the <span style="color:red">R</span> function `seq()`.  This is useful say for a monthly dataset where only the December values are desired, or you want to subsample a very high resolution dataset.  The default is a stride of 1 for all dimensions.  If you want to change the stride value in any dimension, then the value must be given for all dimensions.  So in the previous example, if only every fifth longitude is desired, the call would be:


```r
murSST <- griddap(sstInfo, latitude = c(22., 51.), longitude = c(-140., -105), time = c("2017-01-01", "2017-01-02"), stride = c(1,1,5), fields = 'analysed_sst')

```

Strides are done in array space, not in coordinate space - so it is not skipping say a number of degrees of longitude, but is skipping a number of values of the array index - if longitude is thought of as an array, then every fifth value is used.  There are many cases where having strides work in coordinate space would be handy,  but it can cause a lot of problems.  Consider the case where neither the starting longitude, nor the ending longitude in the request lie on the actual data grid,  and the stride is in coordinate units in such a way that no value requested actually lies on an actual grid value. This would be equivalent to the more complicated problem of data regridding and data interpolation.

When <span style="color:red">ERDDAP</span> receives a request where the bounds are not on the actual dataset grid,  <span style="color:red">ERDDAP</span> finds the closest values on the grid to the requested bounds, and returns those points and all grid points in between.  If a stride is added of a value greater than one but resricted to array space, this guarantees that every value returned lies on the dataset grid.




### Subsetting tabledap()

Tables in <span style="color:red">ERDDAP</span> are subset by using "constraint expressions" on any variable in the table.  The valid operators are =, != (not equals), =~ (a regular expression test), <, <=, >, and >=.  The constraint is constructed as the parameter on the left, value on the right, and the operator in the middle, all within a set of quotes. For example, if in the SWFSC/FRD trawl catch data (datasetID 'FRDCPSTrawlLHHaulCatch'), only sardines for 2010  were desired,  the following constraints would be set in the `tabledap()` call:


```r
'time>=2010-01-01'
'time<=2010-12-31'
'scientific_name="Sardinops sagax"'

```

Note that in `tabledap()` character strings usually must be passed as "double-quoted", as seen in the example with the scientific name. A full `tabledap()` request to retrieve 'latitude',  'longitude', 'time', 'scientific_name', and 'subsample_count' with these constraints would be:


```r
CPSinfo <- info('FRDCPSTrawlLHHaulCatch')
sardines <- tabledap(CPSinfo, fields = c('latitude',  'longitude', 'time', 'scientific_name', 'subsample_count'), 'time>=2010-01-01', 'time<=2010-12-31', 'scientific_name="Sardinops sagax"' )
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
```

##  Searching

The functions in `rerddap` work with any <span style="color:red">ERDDAP</span> server as long as the base URL is provided.  A list of advertised <span style="color:red">ERDDAP</span> servers is provided at https://irishmarineinstitute.github.io/awesome-erddap and the `rerddap` function `servers()` will download this list.

Given a base <span style="color:red">ERDDAP</span> URL, the function `ed_search()` will search the server for the given search terms.  The default <span style="color:red">ERDDAP</span> server is https://upwell.pfeg.noaa.gov/erddap/.  Alternately, the function `global_search()` will search a list of <span style="color:red">ERDDAP</span> servers for the given search terms.  The function `ed_datasets()` lists the datasets available from the given <span style="color:red">ERDDAP</span> server.


## griddap



### MUR SST

MUR (Multi-scale Ultra-high Resolution) is an analyzed SST product at 0.01-degree resolution going back to 2002, providing one of the longest satellite based time series at such high resolution (see https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1). The latest data available for a region off the west coast can be extracted and plotted by:


```r
require("ggplot2")
require("mapdata")
require("rerddap")
sstInfo <- info('jplMURSST41')
# get latest daily sst
murSST <- griddap(sstInfo, latitude = c(22., 51.), longitude = c(-140., -105), time = c('last','last'), fields = 'analysed_sst')
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
mycolor <- colors$temperature
w <- map_data("worldHires", ylim = c(22., 51.), xlim = c(-140, -105))
ggplot(data = murSST$data, aes(x = lon, y = lat, fill = analysed_sst)) +
    geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
    geom_raster(interpolate = FALSE) +
    scale_fill_gradientn(colours = mycolor, na.value = NA) +
    theme_bw() + ylab("latitude") + xlab("longitude") +
    coord_fixed(1.3, xlim = c(-140, -105),  ylim = c(22., 51.)) + ggtitle("Latest MUR SST")
```

<img src="man/figures/MUR-1.png" title="plot of chunk MUR" alt="plot of chunk MUR" style="display: block; margin: auto;" />



### VIIRS SST and Chlorophyll

VIIRS (Visible Infrared Imaging Radiometer Suite)  is a scanning radiometer, that collects visible and infrared imagery and radiometric measurements of the land, atmosphere, cryosphere, and oceans. VIIRS data is used to measure cloud and aerosol properties, ocean color, sea and land surface temperature, ice motion and temperature, fires, and Earth's albedo.   Both NASA and NOAA provide VIIRS-based high resolution SST and chlorophyll products.

ERD provides a 3-day composite SST product at 750 meter resolution developed from a real-time NOAA product. The most recent values can be obtained by setting "time" to be "last".  (Note that <span style="color:red">R</span> sees the latitude-longitude grid as slightly uneven (even though it is in fact even), and that produces artificial lines in `ggplot2::geom_raster()`.  In order to remove those lines, the latitude-longitude grid is remapped to an evenly-space grid.)


```r
require("ggplot2")
require("mapdata")
require("rerddap")
sstInfo <- info('erdVHsstaWS3day')
# get latest 3-day composite sst
viirsSST <- griddap(sstInfo, latitude = c(41., 31.), longitude = c(-128., -115), time = c('last','last'), fields = 'sst')
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
# remap latitiudes and longitudes to even grid
myLats <- unique(viirsSST$data$lat)
myLons <- unique(viirsSST$data$lon)
myLats <- seq(range(myLats)[1], range(myLats)[2], length.out = length(myLats))
myLons <- seq(range(myLons)[1], range(myLons)[2], length.out = length(myLons))
# melt these out to full grid
mapFrame <- expand.grid(x = myLons, y = myLats)
mapFrame$y <- rev(mapFrame$y)
# form a frame with the new values and the data
tempFrame <- data.frame(sst = viirsSST$data$sst, lat = mapFrame$y, lon = mapFrame$x)
mycolor <- colors$temperature
w <- map_data("worldHires", ylim = c(30., 42.), xlim = c(-128, -114))
ggplot(data = tempFrame, aes(x = lon, y = lat, fill = sst)) +
    geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
    geom_raster(interpolate = FALSE) +
    scale_fill_gradientn(colours = mycolor, na.value = NA) +
    theme_bw() + ylab("latitude") + xlab("longitude") +
    coord_fixed(1.3, xlim = c(-128, -114),  ylim = c(30., 42.)) + ggtitle("Latest VIIRS 3-day SST")
```

<img src="man/figures/VIIRS-1.png" title="plot of chunk VIIRS" alt="plot of chunk VIIRS" style="display: block; margin: auto;" />

A time series from the same dataset at a given location,  here (36., -126.):


```r
require("ggplot2")
require("rerddap")
viirsSST1 <- griddap(sstInfo, latitude = c(36., 36.), 
                     longitude = c(-126., -126.), 
                     time = c('2015-01-01','2015-12-31'), fields = 'sst')
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
tempTime <- as.Date(viirsSST1$data$time, origin = '1970-01-01', tz = "GMT")
tempFrame <- data.frame(time = tempTime, sst = viirsSST1$data$sst)
```


```r
ggplot(tempFrame, aes(time, sst)) + 
  geom_line() + 
  theme_bw() + 
  ylab("sst") +
  ggtitle("VIIRS SST at (36N, 126W)")
```

![plot of chunk unnamed-chunk-18](man/figures/unnamed-chunk-18-1.png)



A similar 3-day composite for chloropyll for the same region from a scientific quality product developed by NOAA:


```r
require("ggplot2")
require("mapdata")
require("rerddap")
chlaInfo <- info('erdVHNchla3day')
viirsCHLA <- griddap(chlaInfo, latitude = c(41., 31.), 
                     longitude = c(-128., -115), time = c('last','last'), 
                     fields = 'chla')
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
```


```r
# get latest 3-day composite sst
mycolor <- colors$chlorophyll
w <- map_data("worldHires", ylim = c(30., 42.), xlim = c(-128, -114))
ggplot(data = viirsCHLA$data, aes(x = lon, y = lat, fill = log(chla))) +
  geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_raster(interpolate = FALSE) +
  scale_fill_gradientn(colours = mycolor, na.value = NA) +
  theme_bw() + ylab("latitude") + xlab("longitude") +
  coord_fixed(1.3, xlim = c(-128, -114),  ylim = c(30., 42.)) + 
  ggtitle("Latest VIIRS 3-day Chla")
#> Warning: Removed 2015680 rows containing missing values (geom_raster).
```

![plot of chunk unnamed-chunk-19](man/figures/unnamed-chunk-19-1.png)

### Temperature at 70m in the north Pacific from the SODA model output

This is an example of an extract from a 4-D dataset (results from the "Simple Ocean Data Assimilation (SODA)" model - see http://www.atmos.umd.edu/~ocean/), and illustrate the case where the z-coordinate does not have the default name "altitude".  Water temperature at 70m depth is extracted for the North Pacific Ocean.



```r
require("rerddap")
dataInfo <- rerddap::info('hawaii_d90f_20ee_c4cb')
xpos <- c(135.25, 240.25)
ypos <- c(20.25, 60.25)
zpos <- c(70.02, 70.02)
tpos <- c('2010-12-15', '2010-12-15')
soda70 <- griddap(dataInfo,  longitude = xpos, latitude = ypos, 
                  time = tpos, depth = zpos, fields = 'temp' )
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
str(soda70$data)
#> 'data.frame':	17091 obs. of  5 variables:
#>  $ time : chr  "2010-12-15T00:00:00Z" "2010-12-15T00:00:00Z" "2010-12-15T00:00:00Z" "2010-12-15T00:00:00Z" ...
#>  $ lat  : num  20.2 20.2 20.2 20.2 20.2 ...
#>  $ lon  : num  135 136 136 137 137 ...
#>  $ depth: num  70 70 70 70 70 ...
#>  $ temp : num  27.1 27.1 27.5 27.8 28.2 ...
```

Since the data cross the dateline, it is necessary to use the new "world2Hires" continental outlines in the package `mapdata` which is Pacific Ocean centered.  Unfortunatley there is a small problem where the outlines from certain countries wrap and mistakenly appear in plots, and those countries must be removed,  see code below.



```r
require("ggplot2")
require("mapdata")
xlim <- c(135, 240)
ylim <- c(20, 60)
my.col <- colors$temperature
## Must do a kludge to remove countries that wrap and mess up the plot
w1 <- map("world2Hires", xlim = c(135, 240), ylim = c(20, 60), fill = TRUE, plot = FALSE)
remove <- c("UK:Great Britain", "France", "Spain", "Algeria", "Mali", "Burkina Faso", "Ghana", "Togo")
w <- map_data("world2Hires", regions = w1$names[!(w1$names %in% remove)], ylim = ylim, xlim = xlim)
myplot <- ggplot() +
    geom_raster(data = soda70$data, aes(x = lon, y = lat, fill = temp), interpolate = FALSE) +
    geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
    theme_bw() + scale_fill_gradientn(colours = my.col, na.value = NA, limits = c(-3,30), name = "temperature") +
    ylab("latitude") + xlab("longitude") +
    coord_fixed(1.3, xlim = xlim, ylim = ylim) +
    ggtitle(paste("70m temperature ", soda70$data$time[1]))
myplot
```

<img src="man/figures/soda70Plot-1.png" title="plot of chunk soda70Plot" alt="plot of chunk soda70Plot" style="display: block; margin: auto;" />



### Irish Marine Institute {#hourly}

The Irish Marine Institute has an <span style="color:red">ERDDAP</span> server at http://erddap.marine.ie/erddap.  Among other datasets, they have hourly output from a model of the North Altantic ocean, with a variety of ocean related parameters, see http://erddap.marine.ie/erddap/griddap/IMI_NEATL.html.  To obtain the latest sea surface salinity for the domain of the model:


```r
require("rerddap")
urlBase <- "http://erddap.marine.ie/erddap/"
parameter <- "sea_surface_salinity"
sssTimes <- c("last", "last")
sssLats <- c(48.00625, 57.50625)
sssLons <- c(-17.99375, -1.00625)
dataInfo <- rerddap::info("IMI_NEATL", url = urlBase)
NAtlSSS <- griddap(dataInfo, longitude = sssLons, latitude = sssLats, time = sssTimes, fields = parameter, url = urlBase)
#> info() output passed to x; setting base url to: http://erddap.marine.ie/erddap
str(NAtlSSS$data)
#> 'data.frame':	1034960 obs. of  4 variables:
#>  $ time                : chr  "2021-08-14T00:00:00Z" "2021-08-14T00:00:00Z" "2021-08-14T00:00:00Z" "2021-08-14T00:00:00Z" ...
#>  $ lat                 : num  48 48 48 48 48 ...
#>  $ lon                 : num  -18 -18 -18 -18 -17.9 ...
#>  $ sea_surface_salinity: num  35.5 35.5 35.5 35.5 35.5 ...
```


and the extracted data plotted:


```r
require("ggplot2")
require("mapdata")
xlim <- c(-17.99375, -1.00625)
ylim <- c(48.00625, 57.50625)
my.col <- colors$salinity
w <- map_data("worldHires", ylim = ylim, xlim = xlim)
myplot <- ggplot() +
    geom_raster(data = NAtlSSS$data, aes(x = lon, y = lat, fill = sea_surface_salinity), interpolate = FALSE) +
    geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
    theme_bw() + scale_fill_gradientn(colours = my.col, na.value = NA, limits = c(34, 36), name = "salinity") +
    ylab("latitude") + xlab("longitude") +
    coord_fixed(1.3, xlim = xlim, ylim = ylim) +
    ggtitle(paste("salinity", NAtlSSS$data$time[1]))
myplot
```

<img src="man/figures/NAtlSSSplot-1.png" title="plot of chunk NAtlSSSplot" alt="plot of chunk NAtlSSSplot" style="display: block; margin: auto;" />



### IFREMER

The French agency IFREMER also has an <span style="color:red">ERDDAP</span> server. Here salinity data at 75 meters from "Global Ocean, Coriolis Observation Re-Analysis CORA4.1" model off the west coast of the United States is extracted and plotted.


```r
require("rerddap")
urlBase <- "http://www.ifremer.fr/erddap/"
parameter <- "PSAL"
ifrTimes <- c("2013-05-15", "2013-05-15")
ifrLats <- c(30., 50.)
ifrLons <- c(-140., -110.)
ifrDepth <- c(75., 75.)
dataInfo <- rerddap::info("ifremer_tds0_6080_109e_ed80", url = urlBase)
#> Error: 'Error {
#>     code=404;
#>     message="Not Found: Currently unknown datasetID=ifremer_tds0_6080_109e_ed80";
#> }
#> ' does not exist in current working directory ('/Users/rmendels/WorkFiles/rerddap/vignettes').
ifrPSAL <- griddap(dataInfo, longitude = ifrLons, latitude = ifrLats, time = ifrTimes, depth = ifrDepth,  fields = parameter, url = urlBase)
#> info() output passed to x; setting base url to: http://erddap.marine.ie/erddap
#> Error: Some input dimensions (longitude, latitude, time, depth) don't match those in dataset (time, latitude, longitude)
str(ifrPSAL$data)
#> Error in str(ifrPSAL$data): object 'ifrPSAL' not found
```

The `ggplot2` function `geom_raster()` is not designed for unevenly spaced coordinates, as are the latitudes from this model.  The function `interp()` from the package `akima` is used to interpolate the data which are then plotted.



```r
## ggplot2 has trouble with unequal y's
 require("akima")
 require("dplyr")
 require("ggplot2")
 require("mapdata")
  xlim <- c(-140, -110)
  ylim <- c(30, 51)
## ggplot2 has trouble with unequal y's
  my.col <- colors$salinity
  tempData1 <- ifrPSAL$data$PSAL
#> Error in eval(expr, envir, enclos): object 'ifrPSAL' not found
  tempData <- array(tempData1 , 61 * 54)
#> Error in array(tempData1, 61 * 54): object 'tempData1' not found
  tempFrame <- data.frame(x = ifrPSAL$data$lon, y = ifrPSAL$data$lat)
#> Error in data.frame(x = ifrPSAL$data$lon, y = ifrPSAL$data$lat): object 'ifrPSAL' not found
  tempFrame$temp <- tempData
#> Error in eval(expr, envir, enclos): object 'tempData' not found
  tempFrame1 <- dplyr::filter(tempFrame, !is.nan(temp))
#> Error: Problem with `filter()` input `..1`.
#> ℹ Input `..1` is `!is.nan(temp)`.
#> ✖ object 'temp' not found
  myinterp <- akima::interp(tempFrame1$x, tempFrame1$y, tempFrame1$temp, xo = seq(min(tempFrame1$x), max(tempFrame1$x), length = 61), yo = seq(min(tempFrame1$y), max(tempFrame1$y), length = 54))
#> Error in akima::interp(tempFrame1$x, tempFrame1$y, tempFrame1$temp, xo = seq(min(tempFrame1$x), : object 'tempFrame1' not found
  myinterp1 <- expand.grid(x = myinterp$x, y = myinterp$y)
#> Error in expand.grid(x = myinterp$x, y = myinterp$y): object 'myinterp' not found
  myinterp1$temp <- array(myinterp$z, 61 * 54)
#> Error in array(myinterp$z, 61 * 54): object 'myinterp' not found
  w <- map_data("worldHires", ylim = ylim, xlim = xlim)
 myplot <- ggplot() +
    geom_raster(data = myinterp1, aes(x = x, y = y, fill = temp), interpolate = FALSE) +
    geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
    theme_bw() + scale_fill_gradientn(colours = my.col, na.value = NA, limits = c(32, 35), name = "salinity") +
    ylab("latitude") + xlab("longitude") +
    coord_fixed(1.3, xlim = xlim, ylim = ylim) + ggtitle(paste("salinity at 75 meters",ifrPSAL$data$time[1] ))
#> Error in fortify(data): object 'myinterp1' not found
 myplot
```

<img src="man/figures/ifrPSALplot-1.png" title="plot of chunk ifrPSALplot" alt="plot of chunk ifrPSALplot" style="display: block; margin: auto;" />


## tabledap


### CalCOFI data

CalCOFI (California Cooperative Oceanic Fisheries Investigations - http://www.calcofi.org) is a multi-agency partnership formed in 1949 to investigate the collapse of the sardine population off California. The organization's members are from NOAA Fisheries Service, Scripps Institution of Oceanography, and California Department of Fish and Wildlife. The scope of this research has evolved into the study of marine ecosystems off California and the management of its fisheries resources.  The nearly complete CalCOFI data, both physical and biological, are available through <span style="color:red">ERDDAP</span>.

The following example is a modification of a script developed by Dr. Andrew Leising of the Southwest Fisheries Science Center.  The original script has been used to automate the generation of several yearly reports about the California Current Ecosystem.   The script gets chlorophyll data and a measure of primary productivity from the hydrocasts,and then calculates a seasoanlly adjusted chlorophyll anomaly as well as a seasonally adjusted primary productivity anomaly.  The first step is to get the information about the particular dataset:


```r
require("rerddap")
hydroInfo <- info('siocalcofiHydroCasts')
```

and then get the desired data from 1984 through 2014:


```r
require("rerddap")
calcofi.df <- tabledap(hydroInfo, fields = c('cst_cnt',  'date', 'year', 'month', 'julian_date', 'julian_day', 'rpt_line', 'rpt_sta', 'cruz_num', 'intchl', 'intc14', 'time'), 'time>=1984-01-01T00:00:00Z', 'time<=2014-04-17T05:35:00Z')
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
str(calcofi.df)
#> Classes 'tabledap' and 'data.frame':	11072 obs. of  12 variables:
#>  $ cst_cnt    : int  22522 22523 22524 22525 22526 22527 22528 22529 22530 22531 ...
#>  $ date       : chr  "01/04/1984" "01/05/1984" "01/05/1984" "01/05/1984" ...
#>  $ year       : int  1984 1984 1984 1984 1984 1984 1984 1984 1984 1984 ...
#>  $ month      : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ julian_date: int  30685 30686 30686 30686 30686 30686 30687 30687 30687 30687 ...
#>  $ julian_day : int  4 5 5 5 5 5 6 6 6 6 ...
#>  $ rpt_line   : num  93.3 90 90 90 90 90 90 90 90 90 ...
#>  $ rpt_sta    : num  29 35 30 28 32 37 53 53 53 53 ...
#>  $ cruz_num   : chr  "8401" "8401" "8401" "8401" ...
#>  $ intchl     : num  NaN 30.3 27.9 29.9 39.3 36.7 30.3 30.6 21.3 31.1 ...
#>  $ intc14     : chr  "NaN" "NaN" "NaN" "NaN" ...
#>  $ time       : chr  "1984-01-04T22:08:00Z" "1984-01-05T03:20:00Z" "1984-01-05T09:03:00Z" "1984-01-05T12:26:00Z" ...
#>  - attr(*, "datasetid")= chr "siocalcofiHydroCasts"
#>  - attr(*, "path")= chr "~/Library/Caches/R/rerddap/a7a35af414deaac89811bf13fc934978.csv"
#>  - attr(*, "url")= chr "https://upwell.pfeg.noaa.gov/erddap/tabledap/siocalcofiHydroCasts.csv?cst_cnt,date,year,month,julian_date,julia"| __truncated__
```

Both "intchl" and "intC14" are returned as character strings, and are easier to work with as numbers:


```r
calcofi.df$cruz_num <- as.numeric(calcofi.df$cruz_num)
#> Warning: NAs introduced by coercion
calcofi.df$intc14 <- as.numeric(calcofi.df$intc14)
calcofi.df$time <- as.Date(calcofi.df$time, origin = '1970-01-01', tz = "GMT")
```

At this point the requested data are in the <span style="color:red">R</span> workspace - the rest of the code performs the calculations to derive the seasonally adjusted values and plot them.


```r
require("dplyr")

# calculate cruise means
by_cruznum <- group_by(calcofi.df, cruz_num)
tempData <- select(by_cruznum, year, month, cruz_num, intchl, intc14)
CruiseMeans <- summarize(by_cruznum, cruisechl = mean(intchl, na.rm = TRUE), cruisepp = mean(intc14, na.rm = TRUE), year = median(year, na.rm = TRUE), month = median(month, na.rm = TRUE))
tempTimes <- paste0(CruiseMeans$year,'-',CruiseMeans$month,'-1')
cruisetimes <- as.Date(tempTimes, origin = '1970-01-01', tz = "GMT")
CruiseMeans$cruisetimes <- cruisetimes
# calculate monthly "climatologies"
byMonth <- group_by(CruiseMeans, month)
climate <- summarize(byMonth, ppClimate = mean(cruisepp, na.rm = TRUE), chlaClimate = mean(cruisechl, na.rm = TRUE))
# calculate anomalies
CruiseMeans$chlanom <- CruiseMeans$cruisechl - climate$chlaClimate[CruiseMeans$month]
CruiseMeans$ppanom <- CruiseMeans$cruisepp - climate$ppClimate[CruiseMeans$month]
# calculate mean yearly anomaly
byYear <- select(CruiseMeans, year)
tempData <- select(CruiseMeans, year, chlanom, ppanom )
byYear <- group_by(tempData, year)
yearlyAnom <- summarize(byYear, ppYrAnom = mean(ppanom, na.rm = TRUE), chlYrAnom = mean(chlanom, na.rm = TRUE))
yearlyAnom$year <- ISOdate(yearlyAnom$year, 01, 01, hour = 0)
ggplot(yearlyAnom, aes(year, chlYrAnom)) + geom_line() +
  theme_bw() + ggtitle('yearly chla anom')
ggplot(yearlyAnom, aes(year, ppYrAnom)) + geom_line() +
  theme_bw() + ggtitle('yearly pp anom')
```

<img src="man/figures/calCOFIPlot-1.png" title="plot of chunk calCOFIPlot" alt="plot of chunk calCOFIPlot" style="display: block; margin: auto;" /><img src="man/figures/calCOFIPlot-2.png" title="plot of chunk calCOFIPlot" alt="plot of chunk calCOFIPlot" style="display: block; margin: auto;" />



### CPS Trawl Surveys


The CPS (Coastal Pelagic Species) Trawl Life History Length Frequency Data contains the length distribution of a subset of individuals from a species (mainly non-target) caught during SWFSC-FRD fishery independent trawl surveys of coastal pelagic species. Measured lengths for indicated length type (fork, standard, total, or mantle) were grouped in 10 mm bins (identified by the midpoint of the length class) and counts are recorded by sex.

The number and location of sardines (Sardinops sagax) in the tows in March 2010 and 2011 are extracted, and compared with monthly SST from satellites.  First, query the <span style="color:red">ERDDAP</span> server to see if CPS Trawl data are available through the <span style="color:red">ERDDAP</span> server, and if so, obtain the datasetID for the dataset.


```r
require("rerddap")
CPSquery <- ed_search(query = 'CPS Trawl')
CPSquery$alldata[[1]]$summary
#> [1] "Weight in kilograms for all species (identified to lowest taxonomic criteria) caught during SWFSC-FRD fishery independent surveys (including DEPM, ATM, SaKe) of coastal pelagic species using mid-water trawls (with most tows performed near the surface) at position and times listed. Additional information for a subset of individuals from some species can be found in either CPS Trawl Life History Length Frequency or the CPS Trawl Life History Specimen datasets.\n\ncdm_data_type = Point\nVARIABLES:\ncruise\nship\nhaul (Haul Number)\ncollection\nlatitude (Start Latitude, degrees_north)\nlongitude (Start Longitude, degrees_east)\nstop_latitude\nstop_longitude\ntime (Equilibrium Time, seconds since 1970-01-01T00:00:00Z)\nhaulback_time (Haul Back Time, seconds since 1970-01-01T00:00:00Z)\nsurface_temp (Surface Temperature, degree C)\nsurface_temp_method (Surface Temperature Method)\nship_spd_through_water (Ship Speed Through Water, knot)\nitis_tsn (ItisTSN)\nscientific_name\nsubsample_count\nsubsample_weight (kg)\nremaining_weight (kg)\n"
CPSquery$alldata[[1]]$tabledap
#> [1] "https://upwell.pfeg.noaa.gov/erddap/tabledap/FRDCPSTrawlLHHaulCatch"
CPSquery$alldata[[1]]$dataset_id
#> [1] "FRDCPSTrawlLHHaulCatch"
```

Then get the information for the CPS dataset:


```r
require("rerddap")
(CPSinfo <- info('FRDCPSTrawlLHHaulCatch'))
#> <ERDDAP info> FRDCPSTrawlLHHaulCatch 
#>  Base URL: https://upwell.pfeg.noaa.gov/erddap 
#>  Dataset Type: tabledap 
#>  Variables:  
#>      collection: 
#>          Range: 2003, 4319 
#>      cruise: 
#>          Range: 200307, 202103 
#>      haul: 
#>          Range: 1, 175 
#>      haulback_time: 
#>          Range: 1.05772524E9, 1.61794296E9 
#>          Units: seconds since 1970-01-01T00:00:00Z 
#>      itis_tsn: 
#>      latitude: 
#>          Range: 30.7001, 54.3997 
#>          Units: degrees_north 
#>      longitude: 
#>          Range: -134.0793, -117.2888 
#>          Units: degrees_east 
#>      remaining_weight: 
#>          Units: kg 
#>      scientific_name: 
#>      ship: 
#>      ship_spd_through_water: 
#>          Units: knot 
#>      stop_latitude: 
#>          Range: 30.6663, 54.4157 
#>      stop_longitude: 
#>          Range: -134.0793, -117.318 
#>      subsample_count: 
#>      subsample_weight: 
#>          Units: kg 
#>      surface_temp: 
#>          Units: degree C 
#>      surface_temp_method: 
#>      time: 
#>          Range: 1.05772338E9, 1.61796546E9 
#>          Units: seconds since 1970-01-01T00:00:00Z
```

extract the desired CPS data:


```r
require("dplyr")
require("rerddap")
sardines <- tabledap(CPSinfo, fields = c('latitude',  'longitude', 'time', 'scientific_name', 'subsample_count'), 'time>=2010-01-01', 'time<=2012-01-01', 'scientific_name="Sardinops sagax"' )
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
sardines$time <- as.Date(sardines$time, origin = '1970-01-01', tz = "GMT")
sardines$latitude <- as.numeric(sardines$latitude)
sardines$longitude <- as.numeric(sardines$longitude)
sardine2010 <- filter(sardines, time < as.Date('2010-12-01'))
```

and plot the data versus monthly SST values:


```r
# get the dataset info
sstInfo <- info('erdMWsstdmday')
# get 201004 monthly sst
sst201004 <- griddap('erdMWsstdmday', latitude = c(22., 51.), longitude = c(220., 255), time = c('2010-04-16','2010-04-16'), fields = 'sst')
# get 201104 monthly sst
sst201104 <- griddap('erdMWsstdmday', latitude = c(22., 51.), longitude = c(220., 255), time = c('2011-04-16','2011-04-16'), fields = 'sst')
```


```r
# get polygons of coast for this area
w <- map_data("worldHires", ylim = c(22., 51.), xlim = c(220 - 360, 250 - 360))
# plot 201004 sst on the map
sardine2010 <- filter(sardines, time < as.Date('2010-12-01', origin = '1970-01-01', tz = "GMT"))
sardine2011 <- filter(sardines, time > as.Date('2010-12-01', origin = '1970-01-01', tz = "GMT"))
mycolor <- colors$temperature
p1 <- ggplot() +
  geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_raster(data = sst201004$data, aes(x = lon - 360, y = lat, fill = sst), interpolate = FALSE) +
  scale_fill_gradientn(colours = mycolor, na.value = NA, limits = c(5,30)) +
  theme_bw() + ylab("latitude") + xlab("longitude") +
  coord_fixed(1.3, xlim = c(220 - 360, 250 - 360),  ylim = c(22., 51.))

# plot 201104 sst on the map
p2 <- ggplot() +
  geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_raster(data = sst201104$data, aes(x = lon - 360, y = lat, fill = sst), interpolate = FALSE) +
  geom_point(data = sardine2011, aes(x = longitude, y = latitude, colour = subsample_count)) +
  scale_fill_gradientn(colours = mycolor, na.value = NA, limits = c(5,30)) +
  theme_bw() + ylab("latitude") + xlab("longitude") +
  coord_fixed(1.3, xlim = c(220 - 360, 250 - 360),  ylim = c(22., 51.))
p1 + geom_point(data = sardine2010, aes(x = longitude, y = latitude, colour = subsample_count)) + scale_colour_gradient(space = "Lab", na.value = NA, limits = c(0,80))

p2 +   geom_point(data = sardine2011, aes(x = longitude, y = latitude, colour = subsample_count)) + scale_colour_gradient(space = "Lab", na.value = NA, limits = c(0,80))
```

<img src="man/figures/CPSPlot-1.png" title="plot of chunk CPSPlot" alt="plot of chunk CPSPlot" style="display: block; margin: auto;" /><img src="man/figures/CPSPlot-2.png" title="plot of chunk CPSPlot" alt="plot of chunk CPSPlot" style="display: block; margin: auto;" />

Also of interest is the distribution of sardines through the years:


```r
sardinops <- tabledap(CPSinfo, fields = c('longitude', 'latitude', 'time'),  'scientific_name="Sardinops sagax"')
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
sardinops$time <- as.Date(sardinops$time, origin = '1970-01-01', tz = "GMT")
sardinops$year <- as.factor(format(sardinops$time, '%Y'))
sardinops$latitude <- as.numeric(sardinops$latitude)
sardinops$longitude <- as.numeric(sardinops$longitude)
```


```r
xlim <- c(-135, -110)
ylim <- c(30, 51)
coast <- map_data("worldHires", ylim = ylim, xlim = xlim)
ggplot() +
    geom_point(data = sardinops, aes(x = longitude, y = latitude, colour = year)) +
    geom_polygon(data = coast, aes(x = long, y = lat, group = group), fill = "grey80") +
    theme_bw() + ylab("latitude") + xlab("longitude") +
    coord_fixed(1.3, xlim = xlim, ylim = ylim) +
    ggtitle("Location of sardines by year in EPM Trawls")
```

<img src="man/figures/sardinesPlot-1.png" title="plot of chunk sardinesPlot" alt="plot of chunk sardinesPlot" style="display: block; margin: auto;" />

### NDBC Buoys

NOAA's National Data Buoy Center (NDBC) collects world-wide data from buoys in the ocean. <span style="color:red">ERDDAP</span> can be searched for the location of all buoys in a bounding box with latitudes(37N, 47N) and longitudes (124W, 121W) and the results plotted:


```r
# get location and station ID of NDBC buoys in a region
BuoysInfo <- info('cwwcNDBCMet')
locationBuoys <- tabledap(BuoysInfo, distinct = TRUE, fields = c("station", "longitude", "latitude"), "longitude>=-124", "longitude<=-121", "latitude>=37", "latitude<=47")
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
locationBuoys$latitude <- as.numeric(locationBuoys$latitude)
locationBuoys$longitude <- as.numeric(locationBuoys$longitude)
```


```r
xlim <- c(-130, -110)
ylim <- c(35, 50)
coast <- map_data("worldHires", ylim = ylim, xlim = xlim)
ggplot() +
   geom_point(data = locationBuoys, aes(x = longitude , y = latitude, colour = factor(station) )) +
   geom_polygon(data = coast, aes(x = long, y = lat, group = group), fill = "grey80") +
   theme_bw() + ylab("latitude") + xlab("longitude") +
   coord_fixed(1.3, xlim = xlim, ylim = ylim) +
   ggtitle("Location of buoys in given region")
```

<img src="man/figures/NDBC-1.png" title="plot of chunk NDBC" alt="plot of chunk NDBC" style="display: block; margin: auto;" />

Looking at wind speed for 2012 for station "46012"


```r
buoyData <- tabledap(BuoysInfo, fields = c("time", "wspd"), 'station="46012"', 'time>=2012-01-01', 'time<=2013-01-01')
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
buoyData$wspd <- as.numeric(buoyData$wspd)
buoyData$time <- as.Date(buoyData$time, origin = '1970-01-01', tz = "GMT")
```


```r
ggplot(buoyData, aes(time, wspd)) + 
  geom_line() + 
  theme_bw() + 
  ylab("wind speed") +
  ggtitle("Wind Speed in 2012 from buoy 46012 ")
```

<img src="man/figures/NDBCTS-1.png" title="plot of chunk NDBCTS" alt="plot of chunk NDBCTS" style="display: block; margin: auto;" />

###  IOOS Glider Data

The mission of the IOOS Glider DAC is to provide glider operators with a simple process for submitting glider data sets to a centralized location, enabling the data to be visualized, analyzed, widely distributed via existing web services and the Global Telecommunications System (GTS) and archived at the National Centers for Environmental Information (NCEI).
The IOOS Glider Dac is accessible through `rerddap` (http://data.ioos.us/gliders/erddap/).  Extracting and plotting salinity from part of the path of one glider deployed by the Scripps Institution of Oceanography:


```r
urlBase <- "https://data.ioos.us/gliders/erddap/"
gliderInfo <- info("sp064-20161214T1913",  url = urlBase)
glider <- tabledap(gliderInfo, fields = c("longitude", "latitude", "depth", "salinity"), 'time>=2016-12-14', 'time<=2016-12-23', url = urlBase)
#> info() output passed to x; setting base url to: https://data.ioos.us/gliders/erddap
glider$longitude <- as.numeric(glider$longitude)
glider$latitude <- as.numeric(glider$latitude)
glider$depth <- as.numeric(glider$depth)
```


```r
require("plot3D")
scatter3D(x = glider$longitude , y = glider$latitude , z = -glider$depth, colvar = glider$salinity,              col = colors$salinity, phi = 40, theta = 25, bty = "g", type = "p",
           ticktype = "detailed", pch = 10, clim = c(33.2,34.31), clab = 'Salinity',
           xlab = "longitude", ylab = "latitude", zlab = "depth",
           cex = c(0.5, 1, 1.5))
```

<img src="man/figures/glider-1.png" title="plot of chunk glider" alt="plot of chunk glider" style="display: block; margin: auto;" />

## Animal Telemetry Network (ATN)

The Integrated Ocean Observing System Animal Telemetry Network (IOOS ATN) is designed to serve as an access point to search, discover and access animal telemetry data, and associated oceanographic datasets, from a wide variety of species and platforms.  The track of one of the tagged animals is extracted, and plotted with SST value from satellite along the track.



```r
atnURL <- 'http://oceanview.pfeg.noaa.gov/erddap/'
atnInfo <- info('gtoppAT', url = atnURL)
atnData <- tabledap(atnInfo, fields = c("time", "longitude", "latitude"), 'toppID="1807001"', url = atnURL)
#> info() output passed to x; setting base url to: http://oceanview.pfeg.noaa.gov/erddap
atnData$latitude <- as.numeric(atnData$latitude)
atnData$longitude <- as.numeric(atnData$longitude)
ncdcSST = array(NA_real_, dim = length(atnData$time))
ncdcSSTInfo = info('ncdcOisst2Agg')
for (i in 1:length(atnData$time)) {
  extract <- griddap(ncdcSSTInfo, fields = 'sst', latitude = c(atnData$latitude[i], atnData$latitude[i]), longitude = c(atnData$longitude[i], atnData$longitude[i]), time = c(atnData$time[i], atnData$time[i]))
ncdcSST[i] <- extract$data$sst
}
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
#> info() output passed to x; setting base url to: https://upwell.pfeg.noaa.gov/erddap
```


```r
ylim <- c(32.5, 34)
xlim <- c(-119, -116.5)
mycolor <- colors$temperature
w <- map_data("worldHires", ylim = ylim, xlim = xlim)
alldata <- data.frame(sst = ncdcSST, longitude = atnData$longitude - 360, latitude = atnData$latitude)
z <- ggplot(alldata, aes(x = longitude, y = latitude)) +
   geom_point(aes(colour = sst), size = .5)
z + geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
  theme_bw() +
  scale_colour_gradientn(colours = mycolor, limits = c(16.9, 17.3), "SST") +
  coord_fixed(1.3, xlim = xlim, ylim = ylim) + ggtitle("SST Along Track")
```

<img src="man/figures/unnamed-chunk-28-1.png" title="plot of chunk unnamed-chunk-28" alt="plot of chunk unnamed-chunk-28" style="display: block; margin: auto;" />



## California Current System Integrated Ecosystem Assessment (CCSIEA)

The primary goals of the CCIEA are to better understand the web of interactions that drive patterns and trends of components within the California Current ecosystem, and forecast how changing environmental conditions and management actions affect the status of these components. The conceptual model of the social-ecological system of the California Current illustrates how humans and their social systems are inextricably linked to these marine, coastal, and upland environments (see https://www.integratedecosystemassessment.noaa.gov/regions/california-current-region/index.html).

The over 300 indices developed for the CCSIEA are available through `rerddap`.  Here an index of coho abundance in California is compared with the February value of an index of the strength and location of the North Pacific High, developed in "The North Pacific High and wintertime pre-conditioning of California current productivity" by Schroeder et al.  (Geophys. Res. Lett., 40, 541–546, https://doi.org/10.1002/grl.50100 ).


```r
urlBase <- 'http://coastwatch.pfeg.noaa.gov/erddap/'
nphInfo <- info('erdNph', url = urlBase)
nphData <- tabledap(nphInfo, fields = c("year", "maxSLP" ), 'month=2', 'year>=1987', url = urlBase)
#> info() output passed to x; setting base url to: http://coastwatch.pfeg.noaa.gov/erddap
nphData$maxSLP <- as.numeric(nphData$maxSLP)
urlBase <- 'http://oceanview.pfeg.noaa.gov/erddap/'
cohoInfo <- info('cciea_SM_CA_CO_ABND', url = urlBase)
cohoData <- tabledap(cohoInfo, fields = c("abundance", "time"),  url = urlBase)
#> info() output passed to x; setting base url to: http://oceanview.pfeg.noaa.gov/erddap
#> Error: Query error: Unrecognized variable="abundance"type Status reportmessage Query error: Unrecognized variable="abundance"description The server encountered an internal error that prevented it from fulfilling this request.Apache Tomcat/8.5.11
cohoData$abundance <- as.numeric(cohoData$abundance)
#> Error in eval(expr, envir, enclos): object 'cohoData' not found
alldata <- data.frame(coho = log(cohoData$abundance[1:27]), maxSLP = nphData$maxSLP, year = nphData$year)
#> Error in data.frame(coho = log(cohoData$abundance[1:27]), maxSLP = nphData$maxSLP, : object 'cohoData' not found
```


```r
ggplot(alldata) + geom_line(aes(x = year, y = coho), colour = 'blue') + theme_bw() + ggtitle("Log(coho)")
#> Error in FUN(X[[i]], ...): object 'year' not found
ggplot(alldata) + geom_line(aes(x = year, y = maxSLP), colour = 'red') + theme_bw() + ggtitle("MaxSLP")
#> Error in FUN(X[[i]], ...): object 'year' not found
```

<img src="man/figures/unnamed-chunk-30-1.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" style="display: block; margin: auto;" />



## Cacheing, "last", "now", idempotency, and a gotcha

`rerddap` by default caches the requests you make, so that if you happen to make the same request again, the data is restored from the cache, rather than having to go out and retrieve it remotely.  For most applications, this a boon (such as when "knitting" and "reknitting" this document), as it can speed things up when doing a lot of request in a script, and works because in most cases an <span style="color:red">ERDDAP</span> request is "idempotent".  This means that the the request will always return the same thing no matter what requests came before - it doesn't depend on state. However this is not true if the script uses either "last" in `griddap()` or "now" in `tabledap()` as these will return different values as time elapses and data are added to the datasets.  While it is desirable to have <span style="color:red">ERDDAP</span> purely idempotent,  the "last" and "now" constructs are very helpful for people using <span style="color:red">ERDDAP</span> in dashboards, webpages, regular input to models and the like, and the benefits far outweigh the problems.  However, if you are using either "last" or "now" in an `rerddap` based script, you want to be very careful to clear the `rerddap` cache, otherwise the request will be viewed as the same,  and the data from the last request, rather than the latest data, will be returned.  Note that several examples in this vignette use "last", and therefore the graphics may look different depending on when you "knitted" the vignette.

For help in dealing with the cache, see:


```r
?cache_delete
?cache_delete_all
?cache_details
?cache_list
```



## Reading data from a netCDF file. {#ncdf4}

Here we give a brief summary of how to read in part of the data from a netCDF file.  The basic steps are:

* Open the netCDF file
* Map coordinate values to array indices
* Extract the data

A sample netCDF file, "MWsstd1day.nc" is included.  This is a small file and is a toy example,  but the basic principles remain the same for a larger file.

Open the netCDF file:


```r
require("ncdf4")
exampleFile <- system.file("extdata", "MWsstd1day.nc", package = "rerddap")
sstFile <- nc_open(exampleFile)
```

While it is possible to obtain the coordinate values by extracting them from the file,  `ncdf4` by default does so automatically in `nc_open`. The names of the coordinate dimensions can be found from:


```r
names(sstFile$dim)
#> [1] "time"      "altitude"  "latitude"  "longitude"
```

(note that the coordinate names here are given in 'C' order,  while for an extract the coordinates will be in the opposite, "Fortran" order) and the values of any coordinate, say longitude, can be found from:


```r
sstFile$dim$longitude$vals
#>  [1] 210.0000 210.0125 210.0250 210.0375 210.0500 210.0625 210.0750 210.0875
#>  [9] 210.1000 210.1125 210.1250 210.1375 210.1500 210.1625 210.1750 210.1875
#> [17] 210.2000 210.2125 210.2250 210.2375 210.2500 210.2625 210.2750 210.2875
#> [25] 210.3000 210.3125 210.3250 210.3375 210.3500 210.3625 210.3750 210.3875
#> [33] 210.4000 210.4125 210.4250 210.4375 210.4500 210.4625 210.4750 210.4875
#> [41] 210.5000 210.5125 210.5250 210.5375 210.5500 210.5625 210.5750 210.5875
#> [49] 210.6000 210.6125 210.6250 210.6375 210.6500 210.6625 210.6750 210.6875
#> [57] 210.7000 210.7125 210.7250 210.7375 210.7500 210.7625 210.7750 210.7875
#> [65] 210.8000 210.8125 210.8250 210.8375 210.8500 210.8625 210.8750 210.8875
#> [73] 210.9000 210.9125 210.9250 210.9375 210.9500 210.9625 210.9750 210.9875
#> [81] 211.0000
```

The names of the variables in the netCDF file can be found from:


```r
names(sstFile$var)
#> [1] "sst"
```


An extract is done by giving the pointer to the netCDF file ("sstFile" in this instance), the name of the variable ot be extracted ("sst" in this instance), the  starting index value (in array coordinates) for each dimension, and the the count (the number of other index values) to include.  If all values of a particular dimension are wanted, then "-1" can be used for the count.  So for example to extract all of the values for the first day:


```r
require("ncdf4")
day1SST <- ncvar_get(sstFile, "sst", start = c(1, 1, 1, 1), count = c(1, 1, -1, -1))
```

Suppose we only want the latitudes from (30.0, 30.5) and longitudes (210.0, 210.5) for the first day.  We need to find the array indices that match those coordinate values:


```r
latMin <- which(sstFile$dim$latitude$vals == 30.0)
latMax <- which(sstFile$dim$latitude$vals == 30.5)
lonMin <- which(sstFile$dim$longitude$vals == 210.0)
lonMax <- which(sstFile$dim$longitude$vals == 210.5)
```

and then extract the data:


```r
require("ncdf4")
day1SST <- ncvar_get(sstFile, "sst", start = c(lonMin, latMin, 1, 1), count = c( (lonMax - lonMin + 1), (latMax - latMin + 1), 1, 1 ))
```


If we wanted a time series at (30N, 210E) it would be:


```r
require("ncdf4")
day1SST <- ncvar_get(sstFile, "sst", start = c(lonMax, latMin, 1, 1), count = c(1, 1, 1, -1 ))
```

For this example, it is easy to visually peruse the dimension values but for a large extract this might not be the possible.  Suppose we wanted all values in a latitude range or (30.1, 30.3) and since the range of values we are interested in might not be on a grid boundary we want the smallest range that would include these values  (that is if not on the grid then the smallest value may be less than 30.1)  and similarly for longitude, say (210.1, 210.3):


```r
latMin <- max(which(sstFile$dim$latitude$vals <= 30.1))
latMax <- min(which(sstFile$dim$latitude$vals >= 30.3))
lonMin <- max(which(sstFile$dim$longitude$vals <= 210.1))
lonMax <- min(which(sstFile$dim$longitude$vals >= 210.3))
```

and then perform the extract as before:


```r
require("ncdf4")
day1SST <- ncvar_get(sstFile, "sst", start = c(lonMin, latMin, 1, 1), count = c( (lonMax - lonMin + 1), (latMax - latMin + 1), 1, 1 ))
```
---
title: rerddap introduction
author: Scott Chamberlain
date: "2021-11-18"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{rerddap introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---




`rerddap` is a general purpose R client for working with ERDDAP servers. ERDDAP is a server built on top of OPenDAP, which serves some NOAA data. You can get gridded data ([griddap](https://upwell.pfeg.noaa.gov/erddap/griddap/documentation.html)), which lets you query from gridded datasets, or table data ([tabledap](https://upwell.pfeg.noaa.gov/erddap/tabledap/documentation.html)) which lets you query from tabular datasets. In terms of how we interface with them, there are similarties, but some differences too. We try to make a similar interface to both data types in `rerddap`.

## NetCDF

`rerddap` supports NetCDF format, and is the default when using the `griddap()` function. NetCDF is a binary file format, and will have a much smaller footprint on your disk than csv. The binary file format means it's harder to inspect, but the `ncdf4` package makes it easy to pull data out and write data back into a NetCDF file. Note the the file extension for NetCDF files is `.nc`. Whether you choose NetCDF or csv for small files won't make much of a difference, but will with large files.

## Caching

Data files downloaded are cached in a single hidden directory `~/.rerddap` on your machine. It's hidden so that you don't accidentally delete the data, but you can still easily delete the data if you like.

When you use `griddap()` or `tabledap()` functions, we construct a MD5 hash from the base URL, and any query parameters - this way each query is separately cached. Once we have the hash, we look in `~/.rerddap` for a matching hash. If there's a match we use that file on disk - if no match, we make a http request for the data to the ERDDAP server you specify.

## ERDDAP servers

You can get a data.frame of ERDDAP servers using the function `servers()`. The list of ERDDAP servers is drawn from the [Awesome ERDDAP site](https://irishmarineinstitute.github.io/awesome-erddap).  If you know of more ERDDAP servers, follow the instructions on that page to add the server.

## Install

Stable version from CRAN


```r
install.packages("rerddap")
```

Or, the development version from GitHub


```r
remotes::install_github("ropensci/rerddap")
```


```r
library("rerddap")
```

## Search

First, you likely want to search for data, specify either `griddadp` or `tabledap`


```r
ed_search(query = 'size', which = "table")
#> # A tibble: 9 × 2
#>   title                                                         dataset_id      
#>   <chr>                                                         <chr>           
#> 1 Channel Islands, Kelp Forest Monitoring, Size and Frequency,… erdCinpKfmSFNH  
#> 2 CalCOFI Larvae Sizes                                          erdCalCOFIlrvsiz
#> 3 CalCOFI Larvae Counts Positive Tows                           erdCalCOFIlrvcn…
#> 4 CalCOFI Tows                                                  erdCalCOFItows  
#> 5 GLOBEC NEP MOCNESS Plankton (MOC1) Data, 2000-2002            erdGlobecMoc1   
#> 6 GLOBEC NEP Vertical Plankton Tow (VPT) Data, 1997-2001        erdGlobecVpt    
#> 7 File Names from the AWS S3 noaa-goes16 Bucket                 awsS3NoaaGoes16 
#> 8 AN EXPERIMENTAL DATASET: Underway Sea Surface Temperature an… nodcPJJU        
#> 9 PacIOOS Beach Camera 001: Waikiki, Oahu, Hawaii               BEACHCAM-001
```


```r
ed_search(query = 'size', which = "grid")
#> # A tibble: 5 × 2
#>   title                                                      dataset_id         
#>   <chr>                                                      <chr>              
#> 1 Extended AVHRR Polar Pathfinder Fundamental Climate Data … noaa_ngdc_da08_dcd…
#> 2 Extended AVHRR Polar Pathfinder Fundamental Climate Data … noaa_ngdc_0fe5_a4b…
#> 3 Extended AVHRR Polar Pathfinder Fundamental Climate Data … noaa_ngdc_5253_bf9…
#> 4 Extended AVHRR Polar Pathfinder Fundamental Climate Data … noaa_ngdc_0f24_2f8…
#> 5 SST and SST Anomaly, NOAA Global Coral Bleaching Monitori… NOAA_DHW_monthly
```

There is now a convenience function to search over a list of ERDDAP servers,  designed to work with the function `servers()`


```r
global_search(query, server_list, which_service)
#> Error in global_search(query, server_list, which_service): could not find function "global_search"
```

## Information

Then you can get information on a single dataset


```r
info('erdCalCOFIlrvsiz')
#> <ERDDAP info> erdCalCOFIlrvsiz 
#>  Base URL: https://upwell.pfeg.noaa.gov/erddap 
#>  Dataset Type: tabledap 
#>  Variables:  
#>      calcofi_species_code: 
#>          Range: 19, 946 
#>      common_name: 
#>      cruise: 
#>      itis_tsn: 
#>      larvae_10m2: 
...
```

## griddap (gridded) data

First, get information on a dataset to see time range, lat/long range, and variables.


```r
(out <- info('erdMBchla1day'))
#> <ERDDAP info> erdMBchla1day 
#>  Base URL: https://upwell.pfeg.noaa.gov/erddap 
#>  Dataset Type: griddap 
#>  Dimensions (range):  
#>      time: (2006-01-01T12:00:00Z, 2021-11-16T12:00:00Z) 
#>      altitude: (0.0, 0.0) 
#>      latitude: (-45.0, 65.0) 
#>      longitude: (120.0, 320.0) 
#>  Variables:  
#>      chlorophyll: 
#>          Units: mg m-3
```

Then query for gridded data using the `griddap()` function


```r
(res <- griddap(out,
  time = c('2015-01-01','2015-01-03'),
  latitude = c(14, 15),
  longitude = c(125, 126)
))
#> <ERDDAP griddap> erdMBchla1day
#>    Path: [/Users/rmendels/Library/Caches/R/rerddap/4d844aa48552049c3717ac94ced5f9b8.nc]
#>    Last updated: [2021-11-18 15:58:06]
#>    File size:    [0.03 mb]
#>    Dimensions (dims/vars):   [4 X 1]
#>    Dim names: time, altitude, latitude, longitude
#>    Variable names: Chlorophyll Concentration in Sea Water
#>    data.frame (rows/columns):   [5043 X 5]
#> # A tibble: 5,043 × 5
#>    time                   lat   lon altitude chlorophyll
#>    <chr>                <dbl> <dbl>    <dbl>       <dbl>
#>  1 2015-01-01T12:00:00Z    14  125         0          NA
#>  2 2015-01-01T12:00:00Z    14  125.        0          NA
#>  3 2015-01-01T12:00:00Z    14  125.        0          NA
#>  4 2015-01-01T12:00:00Z    14  125.        0          NA
#>  5 2015-01-01T12:00:00Z    14  125.        0          NA
#>  6 2015-01-01T12:00:00Z    14  125.        0          NA
#>  7 2015-01-01T12:00:00Z    14  125.        0          NA
#>  8 2015-01-01T12:00:00Z    14  125.        0          NA
#>  9 2015-01-01T12:00:00Z    14  125.        0          NA
#> 10 2015-01-01T12:00:00Z    14  125.        0          NA
#> # … with 5,033 more rows
```

The output of `griddap()` is a list that you can explore further. Get the summary


```r
res$summary
#> $filename
#> [1] "/Users/rmendels/Library/Caches/R/rerddap/4d844aa48552049c3717ac94ced5f9b8.nc"
#> 
#> $writable
#> [1] FALSE
#> 
#> $id
#> [1] 65536
#> 
#> $safemode
#> [1] FALSE
#> 
#> $format
#> [1] "NC_FORMAT_CLASSIC"
#> 
...
```

Get the dimension variables


```r
names(res$summary$dim)
#> [1] "time"      "altitude"  "latitude"  "longitude"
```

Get the data.frame (beware: you may want to just look at the `head` of the data.frame if large)


```r
head(res$data)
#>                   time lat     lon altitude chlorophyll
#> 1 2015-01-01T12:00:00Z  14 125.000        0          NA
#> 2 2015-01-01T12:00:00Z  14 125.025        0          NA
#> 3 2015-01-01T12:00:00Z  14 125.050        0          NA
#> 4 2015-01-01T12:00:00Z  14 125.075        0          NA
#> 5 2015-01-01T12:00:00Z  14 125.100        0          NA
#> 6 2015-01-01T12:00:00Z  14 125.125        0          NA
```

## tabledap (tabular) data


```r
(out <- info('erdCalCOFIlrvsiz'))
#> <ERDDAP info> erdCalCOFIlrvsiz 
#>  Base URL: https://upwell.pfeg.noaa.gov/erddap 
#>  Dataset Type: tabledap 
#>  Variables:  
#>      calcofi_species_code: 
#>          Range: 19, 946 
#>      common_name: 
#>      cruise: 
#>      itis_tsn: 
#>      larvae_10m2: 
...
```


```r
(dat <- tabledap('erdCalCOFIlrvsiz', fields=c('latitude','longitude','larvae_size',
  'scientific_name'), 'time>=2011-01-01', 'time<=2011-12-31'))
#> <ERDDAP tabledap> erdCalCOFIlrvsiz
#>    Path: [/Users/rmendels/Library/Caches/R/rerddap/db7389db5b5b0ed9c426d5c13bc43d18.csv]
#>    Last updated: [2021-11-18 15:58:09]
#>    File size:    [0.05 mb]
#> # A tibble: 1,304 × 4
#>    latitude  longitude  larvae_size scientific_name     
#>    <chr>     <chr>      <chr>       <chr>               
#>  1 32.956665 -117.305   4.5         Engraulis mordax    
#>  2 32.91     -117.4     5.0         Merluccius productus
#>  3 32.511665 -118.21167 2.0         Merluccius productus
#>  4 32.511665 -118.21167 3.0         Merluccius productus
#>  5 32.511665 -118.21167 5.5         Merluccius productus
#>  6 32.511665 -118.21167 6.0         Merluccius productus
#>  7 32.511665 -118.21167 2.8         Merluccius productus
#>  8 32.511665 -118.21167 3.0         Sardinops sagax     
#>  9 32.511665 -118.21167 5.0         Sardinops sagax     
#> 10 32.511665 -118.21167 2.5         Engraulis mordax    
#> # … with 1,294 more rows
```

Since both `griddap()` and `tabledap()` give back data.frame's, it's easy to do downstream manipulation. For example, we can use `dplyr` to filter, summarize, group, and sort:


```r
library("dplyr")
dat$larvae_size <- as.numeric(dat$larvae_size)
dat %>%
  group_by(scientific_name) %>%
  summarise(mean_size = mean(larvae_size)) %>%
  arrange(desc(mean_size))
#> # A tibble: 7 × 2
#>   scientific_name       mean_size
#>   <chr>                     <dbl>
#> 1 Anoplopoma fimbria        23.3 
#> 2 Engraulis mordax           9.26
#> 3 Sardinops sagax            7.28
#> 4 Merluccius productus       5.48
#> 5 Tactostoma macropus        5   
#> 6 Scomber japonicus          3.4 
#> 7 Trachurus symmetricus      3.29
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_setup.R
\name{cache_setup}
\alias{cache_setup}
\alias{cache_info}
\title{Setup cache path}
\usage{
cache_setup(path_suffix = NULL, temp_dir = FALSE)

cache_info()
}
\arguments{
\item{path_suffix}{(character) the path suffix to use for storing cached
files, appended to user cache dir.}

\item{temp_dir}{(logical) if \code{TRUE} use a randomly assigned
\code{tempdir} (and \code{path_suffix} is ignored), if \code{FALSE}, you
can use \code{path_suffix}.}
}
\value{
the full cache path, a directory (character)
}
\description{
Setup cache path
}
\details{
Looks first if the user has set a cache path suffix in an
env var or R option. If not found, proceeds to use a temp directory
if not in interactive mode, but if interactive, asks user to setup a
default cache location that will work across sessions (but user can say
no, in which case a temp directory will be used, and each package
start will require cache setup again)
}
\examples{
\dontrun{
# default path
cache_setup()

# you can define your own path
cache_setup(path = "foobar")

# set a tempdir - better for programming with to avoid prompt
cache_setup(temp_dir = TRUE)

# cache info
cache_info()
}
}
\seealso{
Other cache: 
\code{\link{cache_delete}()},
\code{\link{cache_details}()},
\code{\link{cache_list}()}
}
\concept{cache}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_units.R
\name{convert_units}
\alias{convert_units}
\title{Convert a CF Standard Name to/from a GCMD Science Keyword}
\usage{
convert_units(udunits = NULL, ucum = NULL, url = eurl(), ...)
}
\arguments{
\item{udunits}{character; A UDUNITS character string
https://www.unidata.ucar.edu/software/udunits/}

\item{ucum}{character; A UCUM character string
https://ucum.org/ucum.html}

\item{url}{Base URL of the ERDDAP server. See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Convert a CF Standard Name to/from a GCMD Science Keyword
}
\examples{
 \dontrun{
convert_units(udunits = "degree_C meter-1")
convert_units(ucum = "Cel.m-1")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse}
\alias{browse}
\title{Browse a dataset webpage.}
\usage{
browse(x, url = eurl(), ...)
}
\arguments{
\item{x}{datasetid or an object associated with a datasetid such
\code{\link[=info]{info()}}, \code{\link[=griddap]{griddap()}} or \code{\link[=tabledap]{tabledap()}}}

\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Further args passed on to \code{utils::browseURL}
(must be a named parameter)}
}
\value{
if in interactive mode, opens a URL in your default browser;
if not, then prints the URL in the console
}
\description{
Note that it is an error to call this when \code{base::interactive()}
returns \code{FALSE}
}
\examples{
\dontrun{
if (interactive()) {
# browse by dataset_id
browse('erdATastnhday')

# browse info class
my_info <- info('erdATastnhday')
browse(my_info)

# browse tabledap class
my_tabledap <- tabledap('erdCalCOFIlrvsiz', fields=c('latitude','longitude','larvae_size',
   'itis_tsn'), 'time>=2011-10-25', 'time<=2011-10-31')
browse(my_tabledap)
}}
}
\author{
Ben Tupper \email{btupper@bigelow.org}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerddap-package.r
\docType{data}
\name{variablenames}
\alias{variablenames}
\title{variablenames}
\format{
A character vector
}
\description{
variablenames
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_details.R
\name{cache_details}
\alias{cache_details}
\title{Get details of cached files}
\usage{
cache_details(x)
}
\arguments{
\item{x}{File names}
}
\description{
Get details of cached files
}
\details{
Can be used to list details for all files, both .nc and .csv
types, or details for just individual files of class \code{tabledap},
\code{griddap_nc}, and \code{griddap_csv}
}
\examples{
\dontrun{
# List details for all cached files
cache_details()
}
}
\seealso{
Other cache: 
\code{\link{cache_delete}()},
\code{\link{cache_list}()},
\code{\link{cache_setup}()}
}
\concept{cache}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/servers.R
\name{servers}
\alias{servers}
\title{ERDDAP server URLS and other info}
\usage{
servers(...)
}
\arguments{
\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
data.frame with 3 columns:
\itemize{
\item name (character): ERDDAP name
\item url (character): ERDDAP url
\item public (logical): whether it's public or not
}
}
\description{
ERDDAP server URLS and other info
}
\examples{
\dontrun{
servers()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerddap-package.r
\docType{data}
\name{standardnames}
\alias{standardnames}
\title{standardnames}
\format{
A character vector
}
\description{
standardnames
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/version.R
\name{version}
\alias{version}
\title{Get ERDDAP version}
\usage{
version(url = eurl(), ...)
}
\arguments{
\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Get ERDDAP version
}
\examples{
\dontrun{
version()
ss <- servers()
version(ss$url[2])
version(ss$url[3])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global_search.R
\name{global_search}
\alias{global_search}
\title{global_search}
\usage{
global_search(query, server_list, which_service)
}
\arguments{
\item{query}{(character) Search terms}

\item{server_list}{(list of character) List of ERDDAP servers to search}

\item{which_service}{character) One of tabledep or griddap.}
}
\value{
If successful a dataframe wih columns:
\itemize{
\item title - the dataset title
\item dataset_id - the datasetid on that ERDDAP server
\item url - base url of dataset ERDDAP server
}
if urls are valid,  no match is found,  will return no match found
else returns error message
}
\description{
Search for ERDDAP tabledap or griddap datasets from a list
of ERDDAP servers based on search terms.
}
\details{
Uses the 'reddap' function ed_search() to search over
the list of servers
}
\examples{
# get list of servers know by
# https://irishmarineinstitute.github.io/awesome-erddap
# e_servers <- servers()$url
# select a couple to search
# e_servers <- e_servers[c(1, 40)]
# to meet CRAN time limits will only search 1 place
e_servers <- "https://coastwatch.pfeg.noaa.gov/erddap/"
test_query <- 'NOAA/NCDC Blended Monthly'
query_results <- global_search(test_query, e_servers, "griddap")
}
\seealso{
\code{\link[crul]{HttpClient}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerddap-package.r
\docType{package}
\name{rerddap-package}
\alias{rerddap-package}
\alias{rerddap}
\title{rerddap}
\description{
General purpose R client for ERDDAP servers
}
\section{ERDDAP info}{

NOAA's ERDDAP service holds many datasets of interest. It's built on top of
OPenDAP. You can search for datasets via
\code{\link[=ed_search]{ed_search()}}, list datasets via \code{\link[=ed_datasets]{ed_datasets()}},
get information on a single dataset via \code{\link[=info]{info()}}, then get
data you want for either tabledap type via \code{\link[=tabledap]{tabledap()}}, or
for griddap type via \code{\link[=griddap]{griddap()}}
}

\section{tabledap/griddap}{

tabledap and griddap have different interfaces to query for data, so
\code{\link[=tabledap]{tabledap()}} and \code{\link[=griddap]{griddap()}} are separated out as
separate functions even though some of the internals are the same. In
particular, with tabledap you can query on/subset all variables, whereas
with gridddap, you can only query on/subset the dimension varibles (e.g.,
latitude, longitude, altitude).
}

\section{Data size}{

With griddap data via \code{\link[=griddap]{griddap()}} you can get a lot of
data quickly. Try small searches of a dataset to start to get a sense for
the data, then you can increase the amount of data you get. See
\code{\link[=griddap]{griddap()}} for more details.
}

\section{Caching}{

\pkg{rerddap} by default caches the requests you make, so that if you happen to
make the same request again, the data is restored from the cache, rather than
having to go out and retrieve it remotely.  For most applications, this is good,
as it can speed things up when doing a lot of request in a script, and works
because in most cases an ERDDAP request is "idempotent".  This means that the
the request will always return the same thing no matter what requests came
before - it doesn't depend on state. However this is not true if the
script uses either "last" in \code{\link[=griddap]{griddap()}} or "now" in \code{\link[=tabledap]{tabledap()}} as these
will return different values as time elapses and data are added to the
datasets.  While it is desirable to have ERDDAP purely idempotent,  the
"last" and "now" constructs are very helpful for people using ERDDAP
in dashboards, webpages, regular input to models and the like, and the
benefits far outweigh the problems.  However, if you are using either "last"
or "now" in an \pkg{rerddap} based script, you want to be very careful to
clear the \pkg{rerddap} cache, otherwise the request will be viewed as the
same,  and the data from the last request, rather than the latest data,
will be returned.
}

\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{ed_search}
\alias{ed_search}
\alias{ed_datasets}
\title{Search for ERDDAP tabledep or griddap datasets}
\usage{
ed_search(
  query,
  page = NULL,
  page_size = NULL,
  which = "griddap",
  url = eurl(),
  ...
)

ed_datasets(which = "tabledap", url = eurl())
}
\arguments{
\item{query}{(character) Search terms}

\item{page}{(integer) Page number}

\item{page_size}{(integer) Results per page}

\item{which}{(character) One of tabledep or griddap.}

\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET} (must be
named parameters)}
}
\description{
Search for ERDDAP tabledep or griddap datasets
}
\examples{
\dontrun{
(out <- ed_search(query='temperature'))
out$alldata[[1]]
(out <- ed_search(query='size'))
out$info

# List datasets
ed_datasets('table')
ed_datasets('grid')

# use a different ERDDAP server
## Marine Institute (Ireland)
ed_search("temperature", url = "http://erddap.marine.ie/erddap/")
}
}
\references{
https://upwell.pfeg.noaa.gov/erddap/index.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table.R
\name{tabledap}
\alias{tabledap}
\title{Get ERDDAP tabledap data.}
\usage{
tabledap(
  x,
  ...,
  fields = NULL,
  distinct = FALSE,
  orderby = NULL,
  orderbymax = NULL,
  orderbymin = NULL,
  orderbyminmax = NULL,
  units = NULL,
  url = eurl(),
  store = disk(),
  callopts = list()
)
}
\arguments{
\item{x}{Anything coercable to an object of class info. So the output of
a call to \code{\link[=info]{info()}}, or a datasetid, which will internally be passed
through \code{\link[=info]{info()}}}

\item{...}{Any number of key-value pairs in quotes as query constraints.
See Details & examples}

\item{fields}{Columns to return, as a character vector}

\item{distinct}{If \code{TRUE} ERDDAP will sort all of the rows in the results
table (starting with the first requested variable, then using the second
requested variable if the first variable has a tie, ...), then remove all
non-unique rows of data. In many situations, ERDDAP can return distinct
values quickly and efficiently. But in some cases, ERDDAP must look through
all rows of the source dataset.}

\item{orderby}{If used, ERDDAP will sort all of the rows in the results
table (starting with the first variable, then using the second variable
if the first variable has a tie, ...). Normally, the rows of data in the
response table are in the order they arrived from the data source. orderBy
allows you to request that the results table be sorted in a specific way.
For example, use \code{orderby=c("stationID,time")} to get the results
sorted by stationID, then time. The orderby variables MUST be included in
the list of requested variables in the fields parameter.}

\item{orderbymax}{Give a vector of one or more fields, that must be included
in the fields parameter as well. Gives back data given constraints. ERDDAP
will sort all of the rows in the results table (starting with the first
variable, then using the second variable if the first variable has a
tie, ...) and then just keeps the rows where the value of the last sort
variable is highest (for each combination of other values).}

\item{orderbymin}{Same as \code{orderbymax} parameter, except returns
minimum value.}

\item{orderbyminmax}{Same as \code{orderbymax} parameter, except returns
two rows for every combination of the n-1 variables: one row with the
minimum value, and one row with the maximum value.}

\item{units}{One of 'udunits' (units will be described via the UDUNITS
standard (e.g.,degrees_C)) or 'ucum' (units will be described via the
UCUM standard (e.g., Cel)).}

\item{url}{A URL for an ERDDAP server.
Default: https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{store}{One of \code{disk} (default) or \code{memory}. You can pass
options to \code{disk}}

\item{callopts}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET} (must be
named parameters)}
}
\value{
An object of class \code{tabledap}. This class is a thin wrapper
around a data.frame, so the data you get back is a data.frame with metadata
attached as attributes (datasetid, path (path where the csv is stored on
your machine), url (url for the request))
}
\description{
Get ERDDAP tabledap data.
}
\details{
For key-value pair query constraints, the valid operators are =,
!= (not equals), =~ (a regular expression test), <, <=, >, and >= . For
regular expressions you need to add a regular expression. For others, nothing
more is needed. Construct the entry like \code{'time>=2001-07-07'} with the
parameter on the left, value on the right, and the operator in the middle,
all within a set of quotes. Since ERDDAP accepts values other than \code{=},
we can't simply do \code{time = '2001-07-07'} as we normally would.

Server-side functionality: Some tasks are done server side. You don't have
to worry about what that means. They are provided via parameters in this
function. See \code{distinct}, \code{orderby}, \code{orderbymax},
\code{orderbymin}, \code{orderbyminmax}, and \code{units}.

Data is cached based on all parameters you use to get a dataset, including
base url, query parameters. If you make the same exact call in the same or
a different R session, as long you don't clear the cache, the function only
reads data from disk, and does not have to request the data from the web
again.

If you run into an error like "HTTP Status 500 - There was a (temporary?)
problem. Wait a minute, then try again.". it's likely they are hitting
up against a size limit, and they should reduce the amount of data they
are requesting either via space, time, or variables. Pass in
\code{config = verbose()} to the request, and paste the URL into your
browser to see if the output is garbled to examine if there's a problem
with servers or this package
}
\examples{
\dontrun{
# Just passing the datasetid without fields gives all columns back
tabledap('erdCinpKfmBT')

# Pass time constraints
tabledap('erdCinpKfmBT', 'time>=2006-08-24')

# Pass in fields (i.e., columns to retrieve) & time constraints
tabledap('erdCinpKfmBT',
  fields = c('longitude', 'latitude', 'Aplysia_californica_Mean_Density'),
  'time>=2006-08-24'
)

# Get info on a datasetid, then get data given information learned
info('erdCalCOFIlrvsiz')$variables
tabledap('erdCalCOFIlrvsiz', fields=c('latitude','longitude','larvae_size',
   'itis_tsn'), 'time>=2011-10-25', 'time<=2011-10-31')

# An example workflow
## Search for data
(out <- ed_search(query='fish', which = 'table'))
## Using a datasetid, search for information on a datasetid
id <- out$alldata[[1]]$dataset_id
vars <- info(id)$variables
## Get data from the dataset
vars$variable_name[1:3]
tabledap(id, fields = vars$variable_name[1:3])

# Time constraint
## Limit by time with date only
(info <- info('erdCinpKfmBT'))
tabledap(info, fields = c(
  'latitude','longitude','Haliotis_fulgens_Mean_Density'),
  'time>=2001-07-14')

# Use distinct parameter - compare to distinct = FALSE
tabledap('sg114_3',
   fields=c('longitude','latitude','trajectory'),
   'time>=2008-12-05', distinct = TRUE)

# Use units parameter
## In this example, values are the same, but sometimes they can be different
## given the units value passed
tabledap('erdCinpKfmT', fields=c('longitude','latitude','time','temperature'),
   'time>=2007-09-19', 'time<=2007-09-21', units='udunits')
tabledap('erdCinpKfmT', fields=c('longitude','latitude','time','temperature'),
   'time>=2007-09-19', 'time<=2007-09-21', units='ucum')

# Use orderby parameter
tabledap('erdCinpKfmT', fields=c('longitude','latitude','time','temperature'),
   'time>=2007-09-19', 'time<=2007-09-21', orderby='temperature')
# Use orderbymax parameter
tabledap('erdCinpKfmT', fields=c('longitude','latitude','time','temperature'),
   'time>=2007-09-19', 'time<=2007-09-21', orderbymax='temperature')
# Use orderbymin parameter
tabledap('erdCinpKfmT', fields=c('longitude','latitude','time','temperature'),
   'time>=2007-09-19', 'time<=2007-09-21', orderbymin='temperature')
# Use orderbyminmax parameter
tabledap('erdCinpKfmT', fields=c('longitude','latitude','time','temperature'),
   'time>=2007-09-19', 'time<=2007-09-21', orderbyminmax='temperature')
# Use orderbymin parameter with multiple values
tabledap('erdCinpKfmT',
   fields=c('longitude','latitude','time','depth','temperature'),
   'time>=2007-06-10', 'time<=2007-09-21',
   orderbymax=c('depth','temperature')
)

# Integrate with taxize
out <- tabledap('erdCalCOFIlrvcntHBtoHI',
   fields = c('latitude','longitude','scientific_name','itis_tsn'),
   'time>=2007-06-10', 'time<=2007-09-21'
)
tsns <- unique(out$itis_tsn[1:100])
library("taxize")
classif <- classification(tsns, db = "itis")
head(rbind(classif)); tail(rbind(classif))

# Write to memory (within R), or to disk
(out <- info('erdCinpKfmBT'))
## disk, by default (to prevent bogging down system w/ large datasets)
## the 2nd call is much faster as it's mostly just the time of reading
## in the table from disk
system.time( tabledap('erdCinpKfmBT', store = disk()) )
system.time( tabledap('erdCinpKfmBT', store = disk()) )
## memory
tabledap('erdCinpKfmBT', store = memory())

# use a different ERDDAP server
## NOAA IOOS NERACOOS
url <- "http://www.neracoos.org/erddap/"
tabledap("E01_optics_hist", url = url)
}
}
\references{
https://upwell.pfeg.noaa.gov/erddap/index.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keywords.R
\name{key_words}
\alias{key_words}
\title{Convert a CF Standard Name to/from a GCMD Science Keyword}
\usage{
key_words(cf = NULL, gcmd = NULL, url = eurl(), ...)
}
\arguments{
\item{cf}{character; A cf standard name
http://cfconventions.org/Data/cf-standard-names/27/build/cf-standard-name-table.html}

\item{gcmd}{character; A GCMD science keyword
http://gcmd.gsfc.nasa.gov/learn/keyword_list.html}

\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/. See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Convert a CF Standard Name to/from a GCMD Science Keyword
}
\examples{
 \dontrun{
key_words(cf = "air_pressure")
cat(key_words(cf = "air_pressure"))

# a different ERDDAP server
# key_words(cf = "air_pressure", url = servers()$url[6])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colors.R
\docType{data}
\name{colors}
\alias{colors}
\title{cmocean colors
The cmocean color palette by Kristen Thyng as implemented in the R package "oce"}
\format{
An object of class \code{list} of length 13.
}
\usage{
colors
}
\description{
str(colors)
List of 13
$ viridis
$ cdom
$ chlorophyll
$ density
$ freesurface
$ oxygen
$ par
$ phase
$ salinity
$ temperature
$ turbidity
$ velocity
$ vorticity
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_list.R
\name{cache_list}
\alias{cache_list}
\title{List cached files}
\usage{
cache_list()
}
\description{
List cached files
}
\examples{
\dontrun{
# list files in cache
cache_list()

# List info for files
## download some data first
tabledap('erdCinpKfmBT')
griddap('erdVHNchlamday',
 time = c('2015-04-01','2015-04-10'),
 latitude = c(18, 21),
 longitude = c(-120, -119)
)

(x <- cache_list())
cache_details(x$nc[1])
cache_details(x$csv[1])
cache_details()

# delete files by name in cache
# cache_delete(x$nc[1])
# cache_delete(x$nc[2:3])
}
}
\seealso{
Other cache: 
\code{\link{cache_delete}()},
\code{\link{cache_details}()},
\code{\link{cache_setup}()}
}
\concept{cache}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fipscounty.R
\name{fipscounty}
\alias{fipscounty}
\title{Convert a FIPS County Code to/from a County Name}
\usage{
fipscounty(county = NULL, code = NULL, url = eurl(), ...)
}
\arguments{
\item{county}{character; A county name.}

\item{code}{numeric; A FIPS code.}

\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Convert a FIPS County Code to/from a County Name
}
\examples{
 \dontrun{
fipscounty(code = "06053")
fipscounty(county = "CA, Monterey")
fipscounty(county = "OR, Multnomah")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default-url.R
\name{eurl}
\alias{eurl}
\title{Default ERDDAP server URL}
\usage{
eurl()
}
\description{
Default ERDDAP server URL
}
\details{
default url is https://upwell.pfeg.noaa.gov/erddap/

You can set a default using an environment variable so you
don't have to pass anything to the URL parameter in your
function calls.

In your .Renviron file or similar set a URL for the environment
variable \code{RERDDAP_DEFAULT_URL}, like
\verb{RERDDAP_DEFAULT_URL=https://upwell.pfeg.noaa.gov/erddap/}

It's important that you include a trailing slash in your URL
}
\examples{
eurl()
Sys.setenv(RERDDAP_DEFAULT_URL = "https://google.com")
Sys.getenv("RERDDAP_DEFAULT_URL")
eurl()
Sys.unsetenv("RERDDAP_DEFAULT_URL")
eurl()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerddap-package.r
\docType{data}
\name{keywords}
\alias{keywords}
\title{keywords}
\format{
A character vector
}
\description{
keywords
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_time.R
\name{convert_time}
\alias{convert_time}
\title{Convert a UDUNITS compatible time to ISO time}
\usage{
convert_time(
  n = NULL,
  isoTime = NULL,
  units = "seconds since 1970-01-01T00:00:00Z",
  url = eurl(),
  method = "local",
  ...
)
}
\arguments{
\item{n}{numeric; A unix time number.}

\item{isoTime}{character; A string time representation.}

\item{units}{character; Units to return. Default:
"seconds since 1970-01-01T00:00:00Z"}

\item{url}{Base URL of the ERDDAP server. See \code{\link[=eurl]{eurl()}} for
more information}

\item{method}{(character) One of local or web. Local simply uses
\code{\link[=as.POSIXct]{as.POSIXct()}}, while web method uses the ERDDAP time conversion service
\verb{/erddap/convert/time.txt}}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Convert a UDUNITS compatible time to ISO time
}
\details{
When \code{method = "web"} time zone is GMT/UTC
}
\examples{
 \dontrun{
# local conversions
convert_time(n = 473472000)
convert_time(isoTime = "1985-01-02T00:00:00Z")

# using an erddap web service
convert_time(n = 473472000, method = "web")
convert_time(isoTime = "1985-01-02T00:00:00Z", method = "web")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{disk}
\alias{disk}
\alias{memory}
\title{Options for saving ERDDAP datasets.}
\usage{
disk(path = NULL, overwrite = TRUE)

memory()
}
\arguments{
\item{path}{Path to store files in. A directory, not a file.
Default: the root cache path, see \code{\link{cache_setup}}}

\item{overwrite}{(logical) Overwrite an existing file of the same name?
Default: \code{TRUE}}
}
\description{
Options for saving ERDDAP datasets.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerddap-package.r
\docType{data}
\name{longnames}
\alias{longnames}
\title{longnames}
\format{
A character vector
}
\description{
longnames
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerddap-package.r
\docType{data}
\name{ioos_categories}
\alias{ioos_categories}
\title{ioos_categories}
\format{
A character vector
}
\description{
ioos_categories
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_delete.R
\name{cache_delete}
\alias{cache_delete}
\alias{cache_delete_all}
\title{Delete cached files}
\usage{
cache_delete(x, force = FALSE)

cache_delete_all(force = FALSE)
}
\arguments{
\item{x}{File names}

\item{force}{(logical) Should files be force deleted? Default: \code{FALSE}}
}
\description{
Delete cached files
}
\examples{
\dontrun{
# delete files by name in cache
# cache_delete('9911750294a039b8b517c8bf288978ea.csv')
# cache_delete(c('9911750294a039b8b517c8bf288978ea.csv',
#                'b26825b6737da13d6a52c28c8dfe690f.csv'))

# You can delete from the output of griddap or tabledap fxns
## tabledap
(table_res <- tabledap('erdCinpKfmBT'))
cache_delete(table_res)

## griddap
(out <- info('erdQMekm14day'))
(grid_res <- griddap(out,
 time = c('2015-12-28','2016-01-01'),
 latitude = c(24, 23),
 longitude = c(88, 90)
))
cache_delete(grid_res)
}
}
\seealso{
Other cache: 
\code{\link{cache_details}()},
\code{\link{cache_list}()},
\code{\link{cache_setup}()}
}
\concept{cache}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid.R
\name{griddap}
\alias{griddap}
\title{Get ERDDAP gridded data}
\usage{
griddap(
  x,
  ...,
  fields = "all",
  stride = 1,
  fmt = "nc",
  url = eurl(),
  store = disk(),
  read = TRUE,
  callopts = list()
)
}
\arguments{
\item{x}{Anything coercable to an object of class info. So the output of a
call to \code{\link{info}}, or a datasetid, which will internally be passed
through \code{\link{info}}}

\item{...}{Dimension arguments. See examples. Can be any 1 or more of the
dimensions for the particular dataset - and the dimensions vary by dataset.
For each dimension, pass in a vector of length two, with min and max value
desired. at least 1 required.}

\item{fields}{(character) Fields to return, in a character vector.}

\item{stride}{(integer) How many values to get. 1 = get every value, 2 = get
every other value, etc. Default: 1 (i.e., get every value)}

\item{fmt}{(character) One of csv or nc (for netcdf). Default: nc}

\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{store}{One of \code{\link{disk}} (default) or \code{\link{memory}}. You
can pass options to \code{\link{disk}}. Beware: if you choose \code{fmt="nc"},
we force \code{store=disk()} because nc files have to be written to disk.}

\item{read}{(logical) Read data into memory or not. Does not apply when
\code{store} parameter is set to memory (which reads data into memory).
For large csv, or especially netcdf files, you may want to set this to
\code{FALSE}, which simply returns a summary of the dataset - and you can
read in data piecemeal later. Default: \code{TRUE}}

\item{callopts}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
An object of class \code{griddap_csv} if csv chosen or
\code{griddap_nc} if nc file format chosen.

\itemize{
\item \code{griddap_csv}: a data.frame created from the downloaded csv
data
\item \code{griddap_nc}: a list, with slots "summary" and "data". "summary"
is the unclassed output from \code{ncdf4::nc_open}, from which you can
do any netcdf operations you like. "data" is a data.frame created
from the netcdf data. the data.frame may be empty if there were problems
parsing the netcdf data
}

Both have the attributes: datasetid (the dataset id), path (the path on file
for the csv or nc file), url (the url requested to the ERDDAP server)

If \code{read=FALSE}, the data.frame for \code{griddap_csv}
and the data.frame in the "data" slot is empty for \code{griddap_nc}
}
\description{
Get ERDDAP gridded data
}
\details{
Details:

If you run into an error like "HTTP Status 500 - There was a (temporary?)
problem. Wait a minute, then try again.". it's likely they are hitting
up against a size limit, and they should reduce the amount of data they
are requesting either via space, time, or variables. Pass in
\code{config = verbose()} to the request, and paste the URL into your
browser to see if the output is garbled to examine if there's a problem
with servers or this package
}
\section{Dimensions and Variables}{

ERDDAP grid dap data has this concept of dimenions vs. variables. Dimensions
are things like time, latitude, longitude, altitude, and depth. Whereas
variables are the measured variables, e.g., temperature, salinity, air.

You can't separately adjust values for dimensions for different variables.
So, here's how it's gonna work:

Pass in lower and upper limits you want for each dimension as a vector
(e.g., \code{c(1,2)}), or leave to defaults (i.e., don't pass anything to
a dimension). Then pick which variables you want returned via the
\code{fields} parameter. If you don't pass in options to the \code{fields}
parameter, you get all variables back.

To get the dimensions and variables, along with other metadata for a
dataset, run \code{\link{info}}, and each will be shown, with their min
and max values, and some other metadata.
}

\section{Where does the data go?}{

You can choose where data is stored. Be careful though. You can easily get a
single file of hundreds of MB's (upper limit: 2 GB) in size with a single
request. To the \code{store} parameter, pass \code{\link{memory}} if you
want to store the data in memory (saved as a data.frame), or pass
\code{\link{disk}} if you want to store on disk in a file. Note that
\code{\link{memory}} and \code{\link{disk}} are not character strings, but
function calls. \code{\link{memory}} does not accept any inputs, while
\code{\link{disk}} does. Possibly will add other options, like
\dQuote{sql} for storing in a SQL database.
}

\section{Non-lat/lon grid data}{

Some gridded datasets have latitude/longitude components, but some do not.
When nc format gridded datasets have latitude and longitude we "melt" them into a
data.frame for easy downstream consumption. When nc format gridded datasets do
not have latitude and longitude components, we do not read in the data, throw
a warning saying so. You can readin the nc file yourself with the file path.
CSV format is not affected by this issue as CSV data is easily turned into
a data.frame regardless of whether latitude/longitude data are present.
}

\examples{
\dontrun{
# single variable dataset
## You can pass in the outpu of a call to info
(out <- info('erdVHNchlamday'))
## Or, pass in a dataset id
(res <- griddap('erdVHNchlamday',
 time = c('2015-04-01','2015-04-10'),
 latitude = c(18, 21),
 longitude = c(-120, -119)
))

# multi-variable dataset
(out <- info('erdQMekm14day'))
(res <- griddap(out,
 time = c('2015-12-28','2016-01-01'),
 latitude = c(24, 23),
 longitude = c(88, 90)
))
(res <- griddap(out, time = c('2015-12-28','2016-01-01'),
   latitude = c(24, 23), longitude = c(88, 90), fields = 'mod_current'))
(res <- griddap(out, time = c('2015-12-28','2016-01-01'),
   latitude = c(24, 23), longitude = c(88, 90), fields = 'mod_current',
   stride = c(1,2,1,2)))
(res <- griddap(out, time = c('2015-12-28','2016-01-01'),
   latitude = c(24, 23), longitude = c(88, 90),
   fields = c('mod_current','u_current')))


# Write to memory (within R), or to disk
(out <- info('erdQSwindmday'))
## disk, by default (to prevent bogging down system w/ large datasets)
## you can also pass in path and overwrite options to disk()
(res <- griddap(out,
 time = c('2006-07-11','2006-07-20'),
 longitude = c(166, 170),
 store = disk()
))
## the 2nd call is much faster as it's mostly just the time of reading in
## the table from disk
system.time( griddap(out,
 time = c('2006-07-11','2006-07-15'),
 longitude = c(10, 15),
 store = disk()
) )
system.time( griddap(out,
 time = c('2006-07-11','2006-07-15'),
 longitude = c(10, 15),
 store = disk()
) )

## memory - you have to choose fmt="csv" if you use memory
(res <- griddap("erdMBchla1day",
 time = c('2015-01-01','2015-01-03'),
 latitude = c(14, 15),
 longitude = c(125, 126),
 fmt = "csv", store = memory()
))

## Use ncdf4 package to parse data
info("erdMBchla1day")
(res <- griddap("erdMBchla1day",
 time = c('2015-01-01','2015-01-03'),
 latitude = c(14, 15),
 longitude = c(125, 126)
))

# Get data in csv format
## by default, we get netcdf format data
(res <- griddap('erdMBchla1day',
 time = c('2015-01-01','2015-01-03'),
 latitude = c(14, 15),
 longitude = c(125, 126),
 fmt = "csv"
))

# Use a different ERDDAP server url
## NOAA IOOS PacIOOS
url = "https://cwcgom.aoml.noaa.gov/erddap/"
out <- info("miamiacidification", url = url)
(res <- griddap(out,
 time = c('2019-11-01','2019-11-03'),
 latitude = c(15, 16),
 longitude = c(-90, -88)
))
## pass directly into griddap() - if you pass a datasetid string directly
## you must pass in the url or you'll be querying the default ERDDAP url,
## which isn't the one you want if you're not using the default ERDDAP url
griddap("miamiacidification", url = url,
 time = c('2019-11-01','2019-11-03'),
 latitude = c(15, 16),
 longitude = c(-90, -88)
)

# Using 'last'
## with time
griddap('erdVHNchlamday',
 time = c('last-5','last'),
 latitude = c(18, 21),
 longitude = c(-120, -119)
)
## with latitude
griddap('erdVHNchlamday',
  time = c('2015-04-01','2015-04-10'),
  latitude = c('last', 'last'),
  longitude = c(-120, -119)
)
## with longitude
griddap('erdVHNchlamday',
  time = c('2015-04-01','2015-04-10'),
  latitude = c(18, 21),
  longitude = c('last', 'last')
)

# datasets without lat/lon grid and with fmt=nc
# FIXME: this dataset is gone
# (x <- info('glos_tds_5912_ca66_3f41'))
# res <- griddap(x,
#   time = c('2018-04-01','2018-04-10'),
#   ny = c(1, 2),
#   nx = c(3, 5)
# )
## data.frame is empty
# res$data
## read in from the nc file path
# ncdf4::nc_open(res$summary$filename)
}
}
\references{
https://upwell.pfeg.noaa.gov/erddap/rest.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_adv.R
\name{ed_search_adv}
\alias{ed_search_adv}
\title{Advanced search for ERDDAP tabledep or griddap datasets}
\usage{
ed_search_adv(
  query = NULL,
  page = 1,
  page_size = 1000,
  protocol = NULL,
  cdm_data_type = NULL,
  institution = NULL,
  ioos_category = NULL,
  keywords = NULL,
  long_name = NULL,
  standard_name = NULL,
  variableName = NULL,
  maxLat = NULL,
  minLon = NULL,
  maxLon = NULL,
  minLat = NULL,
  minTime = NULL,
  maxTime = NULL,
  url = eurl(),
  ...
)
}
\arguments{
\item{query}{(character) Search terms}

\item{page}{(integer) Page number. Default: 1}

\item{page_size}{(integer) Results per page: Default: 1000}

\item{protocol}{(character) One of any (default), tabledep or griddap}

\item{cdm_data_type}{(character) One of grid, other, point, profile,
timeseries, timeseriesprofile, trajectory, trajectoryprofile}

\item{institution}{(character) An institution. See the dataset
\code{institutions}}

\item{ioos_category}{(character) An ioos category See the dataset
\code{ioos_categories}}

\item{keywords}{(character) A keywords. See the dataset \code{keywords}}

\item{long_name}{(character) A long name. See the dataset \code{longnames}}

\item{standard_name}{(character) A standar dname. See the dataset
\code{standardnames}}

\item{variableName}{(character) A variable name. See the dataset
\code{variablenames}}

\item{minLon, maxLon}{(numeric) Minimum and maximum longitude. Some datasets
have longitude values within -180 to 180, others use 0 to 360. If you
specify min and max Longitude within -180 to 180 (or 0 to 360), ERDDAP will
only find datasets that match the values you specify. Consider doing one
search: longitude -180 to 360, or two searches: longitude -180 to 180,
and 0 to 360.}

\item{minLat, maxLat}{(numeric) Minimum and maximum latitude, between -90
and 90}

\item{minTime, maxTime}{(numeric/character) Minimum and maximum time. Time
string with the format "yyyy-MM-ddTHH:mm:ssZ, (e.g., 2009-01-21T23:00:00Z).
If you specify something, you must include at least yyyy-MM-dd; you can
omit Z, :ss, :mm, :HH, and T. Always use UTC (GMT/Zulu) time. Or specify
the number of seconds since 1970-01-01T00:00:00Z.}

\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET} (must be
named parameters)}
}
\description{
Advanced search for ERDDAP tabledep or griddap datasets
}
\examples{
\dontrun{
ed_search_adv(query = 'temperature')
ed_search_adv(query = 'temperature', protocol = "griddap")
ed_search_adv(query = 'temperature', protocol = "tabledap")
ed_search_adv(maxLat = 63, minLon = -107, maxLon = -87, minLat = 50,
  protocol = "griddap")
ed_search_adv(maxLat = 63, minLon = -107, maxLon = -87, minLat = 50,
  protocol = "tabledap")
ed_search_adv(minTime = "2010-01-01T00:00:00Z",
  maxTime="2010-02-01T00:00:00Z")
(out <- ed_search_adv(maxLat = 63, minLon = -107, maxLon = -87, minLat = 50,
             minTime = "2010-01-01T00:00:00Z",
             maxTime="2010-02-01T00:00:00Z"))
out$alldata[[1]]
ed_search_adv(variableName = 'upwelling')
ed_search_adv(query = 'upwelling', protocol = "tabledap")

# use a different URL
ed_search_adv(query = 'temperature', url = servers()$url[6])
}
}
\references{
https://upwell.pfeg.noaa.gov/erddap/index.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerddap-package.r
\docType{data}
\name{institutions}
\alias{institutions}
\title{institutions}
\format{
A character vector
}
\description{
institutions
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/info.R
\name{info}
\alias{info}
\alias{as.info}
\title{Get information on an ERDDAP dataset.}
\usage{
info(datasetid, url = eurl(), ...)

as.info(x, url)
}
\arguments{
\item{datasetid}{Dataset id}

\item{url}{A URL for an ERDDAP server. Default:
https://upwell.pfeg.noaa.gov/erddap/ - See \code{\link[=eurl]{eurl()}} for
more information}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET} (must be a
named parameter)}

\item{x}{A datasetid or the output of \code{info}}
}
\value{
Prints a summary of the data on return, but you can index to
various information.

The data is a list of length two with:
\itemize{
\item variables - Data.frame of variables and their types
\item alldata - List of data variables and their full attributes
}

Where \code{alldata} element has many data.frame's, one for each variable,
with metadata for that variable. E.g., for griddap dataset
\code{noaa_pfeg_696e_ec99_6fa6}, \code{alldata} has:
\itemize{
\item NC_GLOBAL
\item time
\item latitude
\item longitude
\item sss
}
}
\description{
Get information on an ERDDAP dataset.
}
\examples{
\dontrun{
# grid dap datasets
info('erdATastnhday')

(out <- ed_search(query='temperature'))
info(out$info$dataset_id[5])
info(out$info$dataset_id[15])
info(out$info$dataset_id[25])
info(out$info$dataset_id[150])
info(out$info$dataset_id[400])
info(out$info$dataset_id[678])

out <- info(datasetid='erdMBchla1day')
## See brief overview of the variables and range of possible values, if given
out$variables
## all information on longitude
out$alldata$longitude
## all information on chlorophyll
out$alldata$chlorophyll

# table dap datasets
(out <- ed_search(query='temperature', which = "table"))
info(out$info$dataset_id[1])
info(out$info$dataset_id[2])
info(out$info$dataset_id[3])
info(out$info$dataset_id[4])

info('erdCinpKfmBT')
out <- info('erdCinpKfmBT')
## See brief overview of the variables and range of possible values, if given
out$variables
## all information on longitude
out$alldata$longitude
## all information on Haliotis_corrugata_Mean_Density
out$alldata$Haliotis_corrugata_Mean_Density

# use a different ERDDAP server
## Marine Institute (Ireland)
info("IMI_CONN_2D", url = "http://erddap.marine.ie/erddap/")
}
}
\references{
https://upwell.pfeg.noaa.gov/erddap/index.html
}
