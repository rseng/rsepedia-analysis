scrubr
======



[![R-check](https://github.com/ropensci/scrubr/workflows/R-check/badge.svg)](https://github.com/ropensci/scrubr/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/scrubr)](https://cranchecks.info/pkgs/scrubr)
[![codecov.io](https://codecov.io/github/ropensci/scrubr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/scrubr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/scrubr?color=ff69b4)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/scrubr)](https://cran.r-project.org/package=scrubr)

__Clean Biological Occurrence Records__

Clean using the following use cases (checkmarks indicate fxns exist - not necessarily complete):

- [x] Impossible lat/long values: e.g., latitude 75
- [x] Incomplete cases: one or the other of lat/long missing
- [x] Unlikely lat/long values: e.g., points at 0,0
- [x] Deduplication: try to identify duplicates, esp. when pulling data from multiple sources, e.g., can try to use occurrence IDs, if provided
- [x] Date based cleaning
* [x] Outside political boundary: User input to check for points in the wrong country, or points outside of a known country
* [x] Taxonomic name based cleaning: via `taxize` (one method so far)
* Political centroids: unlikely that occurrences fall exactly on these points, more likely a
default position (Draft function started, but not exported, and commented out). see [issue #6](https://github.com/ropensci/scrubr/issues/6)
* Herbaria/Museums: many specimens may have location of the collection they are housed in, see [issue #20](https://github.com/ropensci/scrubr/issues/20)
* Habitat type filtering: e.g., fish should not be on land; marine fish should not be in fresh water
* Check for contextually wrong values: That is, if 99 out of 100 lat/long coordinates are within the continental US, but 1 is in China, then perhaps something is wrong with that one point
* Collector/recorder names: see [issue #19](https://github.com/ropensci/scrubr/issues/19)
* ...

A note about examples: We think that using a piping workflow with `%>%` makes code easier to
build up, and easier to understand. However, in some examples we provide examples without the pipe
to demonstrate traditional usage.

## Install

Stable CRAN version


```r
install.packages("scrubr")
```

Development version


```r
remotes::install_github("ropensci/scrubr")
```


```r
library("scrubr")
```

## Coordinate based cleaning


```r
data("sampledata1")
```

Remove impossible coordinates (using sample data included in the pkg)


```r
# coord_impossible(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_impossible()
#> # A tibble: 1,500 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,490 more rows
```

Remove incomplete coordinates


```r
# coord_incomplete(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_incomplete()
#> # A tibble: 1,306 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,296 more rows
```

Remove unlikely coordinates (e.g., those at 0,0)


```r
# coord_unlikely(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_unlikely()
#> # A tibble: 1,488 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,478 more rows
```

Do all three


```r
dframe(sample_data_1) %>%
  coord_impossible() %>%
  coord_incomplete() %>%
  coord_unlikely()
#> # A tibble: 1,294 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,284 more rows
```

Don't drop bad data


```r
dframe(sample_data_1) %>% coord_incomplete(drop = TRUE) %>% NROW
#> [1] 1306
dframe(sample_data_1) %>% coord_incomplete(drop = FALSE) %>% NROW
#> [1] 1500
```


## Deduplicate


```r
smalldf <- sample_data_1[1:20, ]
# create a duplicate record
smalldf <- rbind(smalldf, smalldf[10,])
row.names(smalldf) <- NULL
# make it slightly different
smalldf[21, "key"] <- 1088954555
NROW(smalldf)
#> [1] 21
dp <- dframe(smalldf) %>% dedup()
NROW(dp)
#> [1] 20
attr(dp, "dups")
#> # A tibble: 1 x 5
#>   name             longitude latitude date                       key
#>   <chr>                <dbl>    <dbl> <dttm>                   <dbl>
#> 1 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954555
```

## Dates

Standardize/convert dates


```r
df <- sample_data_1
# date_standardize(dframe(df), "%d%b%Y") # w/o pipe
dframe(df) %>% date_standardize("%d%b%Y")
#> # A tibble: 1,500 x 5
#>    name             longitude latitude date             key
#>    <chr>                <dbl>    <dbl> <chr>          <int>
#>  1 Ursus americanus     -79.7     38.4 14Jan2015 1065590124
#>  2 Ursus americanus     -82.4     35.7 13Jan2015 1065588899
#>  3 Ursus americanus     -99.1     23.7 20Feb2015 1098894889
#>  4 Ursus americanus     -72.8     43.9 13Feb2015 1065611122
#>  5 Ursus americanus     -72.3     43.9 01Mar2015 1088908315
#>  6 Ursus americanus    -109.      32.7 29Mar2015 1088932238
#>  7 Ursus americanus    -109.      32.7 29Mar2015 1088932273
#>  8 Ursus americanus    -124.      40.1 28Mar2015 1132403409
#>  9 Ursus americanus     -78.3     36.9 20Mar2015 1088923534
#> 10 Ursus americanus     -76.8     35.5 05Apr2015 1088954559
#> # … with 1,490 more rows
```

Drop records without dates


```r
NROW(df)
#> [1] 1500
NROW(dframe(df) %>% date_missing())
#> [1] 1498
```

Create date field from other fields


```r
dframe(sample_data_2) %>% date_create(year, month, day)
#> # A tibble: 1,500 x 8
#>    name             longitude latitude        key year  month day   date      
#>    <chr>                <dbl>    <dbl>      <int> <chr> <chr> <chr> <chr>     
#>  1 Ursus americanus     -79.7     38.4 1065590124 2015  01    14    2015-01-14
#>  2 Ursus americanus     -82.4     35.7 1065588899 2015  01    13    2015-01-13
#>  3 Ursus americanus     -99.1     23.7 1098894889 2015  02    20    2015-02-20
#>  4 Ursus americanus     -72.8     43.9 1065611122 2015  02    13    2015-02-13
#>  5 Ursus americanus     -72.3     43.9 1088908315 2015  03    01    2015-03-01
#>  6 Ursus americanus    -109.      32.7 1088932238 2015  03    29    2015-03-29
#>  7 Ursus americanus    -109.      32.7 1088932273 2015  03    29    2015-03-29
#>  8 Ursus americanus    -124.      40.1 1132403409 2015  03    28    2015-03-28
#>  9 Ursus americanus     -78.3     36.9 1088923534 2015  03    20    2015-03-20
#> 10 Ursus americanus     -76.8     35.5 1088954559 2015  04    05    2015-04-05
#> # … with 1,490 more rows
```

## Ecoregion

Filter by FAO areas


```r
wkt <- 'POLYGON((72.2 38.5,-173.6 38.5,-173.6 -41.5,72.2 -41.5,72.2 38.5))'
manta_ray <- rgbif::name_backbone("Mobula alfredi")$usageKey
res <- rgbif::occ_data(manta_ray, geometry = wkt, limit=300, hasCoordinate = TRUE)
dat <- sf::st_as_sf(res$data, coords = c("decimalLongitude", "decimalLatitude"))
dat <- sf::st_set_crs(dat, 4326)
mapview::mapview(dat)
tmp <- eco_region(dframe(res$data), dataset = "fao", region = "OCEAN:Indian")
tmp <- tmp[!is.na(tmp$decimalLongitude), ]
tmp2 <- sf::st_as_sf(tmp, coords = c("decimalLongitude", "decimalLatitude"))
tmp2 <- sf::st_set_crs(tmp2, 4326)
mapview::mapview(tmp2)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/scrubr/issues).
* License: MIT
* Get citation information for `scrubr` in R doing `citation(package = 'scrubr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
scrubr 0.4.0
============

### MINOR IMPROVEMENTS

* only run `coord_within()` examples in docs if `sf` and `rworldmap` are present (and in interactive mode) (#35)
* allow `eco_region()` to filter more than one region (#36)


scrubr 0.3.2
============

### BUG FIXES

* fix some failing cran checks - changed CRS specification when using sf package (#34)


scrubr 0.3.0
============

### NEW FEATURES

* gains function `fix_names()`, ported over from the `spocc` package; a helper function to change taxonomic names in a data.frame to make plotting simpler (#29)
* gains new function `eco_region()` to filter data by ecoregions; also exported are `regions_fao()` and `regions_meow()` that fetch the data used in `eco_region()` so the user can better figure out what variables to use (#30)
* gains new function `coord_imprecise()` to clean imprecise coordinates (#18)
* gains new function `coord_uncertain()` to clean uncertain occurrences, as determined through the `coordinateUncertaintyInMeters` variable reported in Darwin Core records
* now importing `data.table`, `fastmatch`, `crul`, `jsonlite`, `tibble`, `curl` and `hoardr` (`sf` and `mapview` in Suggests)

### MINOR IMPROVEMENTS

* `coord_within()` now uses `sf` instead of `sp` (#31)
* using tibble now for compact, easier to handle, data.frame's (#21)

### BUG FIXES

* fix to `dedup()` to remove duplicate entries (#27)


scrubr 0.1.1
============

### MINOR IMPROVEMENTS

* Fixed examples to be conditional on presence of `rgbif` (#17)
* Fix `as.matrix()` to use `Matrix::as.matrix()`


scrubr 0.1.0
============

### NEW FEATURES

* Releasd to CRAN.
## Test environments

* local macOS install, R 4.1.0
* ubuntu 16.04 (on github actions), R 4.1.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 downstream dependency; there were no errors related to scrubr.

--------

This version fixes a failing cran check, and another minor issue.

Thanks!
Scott Chamberlain
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

### Bugs?

* Submit an issue on the Issues page [here](https://github.com/ropensci/scrubr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/scrubr.git`
* Make sure to track progress upstream (i.e., on our version of `scrubr` at `ropensci/scrubr`) by doing `git remote add upstream https://github.com/ropensci/scrubr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new branch)
* If you alter package functionality at all (e.g., the code itself, not just documentation)
please do write some tests to cover the new functionality.
* Push up to your account
* Submit a pull request to home base at `ropensci/scrubr`
<!-- YOU WILL BE ASKED FOR YOUR SESSION INFO IF YOU DO NOT PROVIDE IT -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
location data
=============

* herbaria.csv - modified from https://github.com/ejedwards/reanalysis_zanne2014
scrubr
======

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![R-check](https://github.com/ropensci/scrubr/workflows/R-check/badge.svg)](https://github.com/ropensci/scrubr/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/scrubr)](https://cranchecks.info/pkgs/scrubr)
[![codecov.io](https://codecov.io/github/ropensci/scrubr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/scrubr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/scrubr?color=ff69b4)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/scrubr)](https://cran.r-project.org/package=scrubr)

__Clean Biological Occurrence Records__

Clean using the following use cases (checkmarks indicate fxns exist - not necessarily complete):

- [x] Impossible lat/long values: e.g., latitude 75
- [x] Incomplete cases: one or the other of lat/long missing
- [x] Unlikely lat/long values: e.g., points at 0,0
- [x] Deduplication: try to identify duplicates, esp. when pulling data from multiple sources, e.g., can try to use occurrence IDs, if provided
- [x] Date based cleaning
* [x] Outside political boundary: User input to check for points in the wrong country, or points outside of a known country
* [x] Taxonomic name based cleaning: via `taxize` (one method so far)
* Political centroids: unlikely that occurrences fall exactly on these points, more likely a
default position (Draft function started, but not exported, and commented out). see [issue #6](https://github.com/ropensci/scrubr/issues/6)
* Herbaria/Museums: many specimens may have location of the collection they are housed in, see [issue #20](https://github.com/ropensci/scrubr/issues/20)
* Habitat type filtering: e.g., fish should not be on land; marine fish should not be in fresh water
* Check for contextually wrong values: That is, if 99 out of 100 lat/long coordinates are within the continental US, but 1 is in China, then perhaps something is wrong with that one point
* Collector/recorder names: see [issue #19](https://github.com/ropensci/scrubr/issues/19)
* ...

A note about examples: We think that using a piping workflow with `%>%` makes code easier to
build up, and easier to understand. However, in some examples we provide examples without the pipe
to demonstrate traditional usage.

## Install

Stable CRAN version

```{r eval=FALSE}
install.packages("scrubr")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/scrubr")
```

```{r}
library("scrubr")
```

## Coordinate based cleaning

```{r}
data("sampledata1")
```

Remove impossible coordinates (using sample data included in the pkg)

```{r}
# coord_impossible(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_impossible()
```

Remove incomplete coordinates

```{r}
# coord_incomplete(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_incomplete()
```

Remove unlikely coordinates (e.g., those at 0,0)

```{r}
# coord_unlikely(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_unlikely()
```

Do all three

```{r}
dframe(sample_data_1) %>%
  coord_impossible() %>%
  coord_incomplete() %>%
  coord_unlikely()
```

Don't drop bad data

```{r}
dframe(sample_data_1) %>% coord_incomplete(drop = TRUE) %>% NROW
dframe(sample_data_1) %>% coord_incomplete(drop = FALSE) %>% NROW
```


## Deduplicate

```{r}
smalldf <- sample_data_1[1:20, ]
# create a duplicate record
smalldf <- rbind(smalldf, smalldf[10,])
row.names(smalldf) <- NULL
# make it slightly different
smalldf[21, "key"] <- 1088954555
NROW(smalldf)
dp <- dframe(smalldf) %>% dedup()
NROW(dp)
attr(dp, "dups")
```

## Dates

Standardize/convert dates

```{r}
df <- sample_data_1
# date_standardize(dframe(df), "%d%b%Y") # w/o pipe
dframe(df) %>% date_standardize("%d%b%Y")
```

Drop records without dates

```{r}
NROW(df)
NROW(dframe(df) %>% date_missing())
```

Create date field from other fields

```{r}
dframe(sample_data_2) %>% date_create(year, month, day)
```

## Ecoregion

Filter by FAO areas

```{r eval=FALSE}
wkt <- 'POLYGON((72.2 38.5,-173.6 38.5,-173.6 -41.5,72.2 -41.5,72.2 38.5))'
manta_ray <- rgbif::name_backbone("Mobula alfredi")$usageKey
res <- rgbif::occ_data(manta_ray, geometry = wkt, limit=300, hasCoordinate = TRUE)
dat <- sf::st_as_sf(res$data, coords = c("decimalLongitude", "decimalLatitude"))
dat <- sf::st_set_crs(dat, 4326)
mapview::mapview(dat)
tmp <- eco_region(dframe(res$data), dataset = "fao", region = "OCEAN:Indian")
tmp <- tmp[!is.na(tmp$decimalLongitude), ]
tmp2 <- sf::st_as_sf(tmp, coords = c("decimalLongitude", "decimalLatitude"))
tmp2 <- sf::st_set_crs(tmp2, 4326)
mapview::mapview(tmp2)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/scrubr/issues).
* License: MIT
* Get citation information for `scrubr` in R doing `citation(package = 'scrubr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "scrubr introduction"
date: "2020-04-06"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{scrubr introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




`scrubr` is a general purpose toolbox for cleaning biological occurrence records. Think
of it like `dplyr` but specifically for occurrence data. It includes functionality for
cleaning based on various aspects of spatial coordinates, unlikely values due to political
centroids, taxonomic names, and more.

## Installation

Install from CRAN


```r
install.packages("scrubr")
```

Or install the development version from GitHub


```r
remotes::install_github("ropensci/scrubr")
```

Load scrubr


```r
library("scrubr")
```

We'll use sample datasets included with the package, they are lazy loaded,
and available via `sample_data_1` and `sample_data_2`

## data.frame's

All functions expect data.frame's as input, and output data.frame's

## Pipe vs. no pipe

We think that using a piping workflow with `%>%` makes code easier to
build up, and easier to understand. However, in some examples below we provide
commented out examples without the pipe to demonstrate traditional usage - which
you can use if you remove the comment `#` at beginning of the line.

## dframe

`dframe()` is a utility function to create a compact data.frame representation. You 
don't have to use it. If you do, you can work with `scrubr` functions with a compact
data.frame, making it easier to see the data quickly. If you don't use `dframe()`
we just use your regular data.frame. Problem is with large data.frame's you deal with 
lots of stuff printed to the screen, making it hard to quickly wrangle data.

## Coordinate based cleaning

Remove impossible coordinates (using sample data included in the pkg)


```r
# coord_impossible(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_impossible()
#> # A tibble: 1,500 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,490 more rows
```

Remove incomplete coordinates


```r
# coord_incomplete(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_incomplete()
#> # A tibble: 1,306 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,296 more rows
```

Remove unlikely coordinates (e.g., those at 0,0)


```r
# coord_unlikely(dframe(sample_data_1)) # w/o pipe
dframe(sample_data_1) %>% coord_unlikely()
#> # A tibble: 1,488 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,478 more rows
```

Do all three


```r
dframe(sample_data_1) %>%
  coord_impossible() %>%
  coord_incomplete() %>%
  coord_unlikely()
#> # A tibble: 1,294 x 5
#>    name             longitude latitude date                       key
#>  * <chr>                <dbl>    <dbl> <dttm>                   <int>
#>  1 Ursus americanus     -79.7     38.4 2015-01-14 16:36:45 1065590124
#>  2 Ursus americanus     -82.4     35.7 2015-01-13 00:25:39 1065588899
#>  3 Ursus americanus     -99.1     23.7 2015-02-20 23:00:00 1098894889
#>  4 Ursus americanus     -72.8     43.9 2015-02-13 16:16:41 1065611122
#>  5 Ursus americanus     -72.3     43.9 2015-03-01 20:20:45 1088908315
#>  6 Ursus americanus    -109.      32.7 2015-03-29 17:06:54 1088932238
#>  7 Ursus americanus    -109.      32.7 2015-03-29 17:12:50 1088932273
#>  8 Ursus americanus    -124.      40.1 2015-03-28 23:00:00 1132403409
#>  9 Ursus americanus     -78.3     36.9 2015-03-20 21:11:24 1088923534
#> 10 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954559
#> # … with 1,284 more rows
```

Don't drop bad data


```r
dframe(sample_data_1) %>% coord_incomplete(drop = TRUE) %>% NROW
#> [1] 1306
dframe(sample_data_1) %>% coord_incomplete(drop = FALSE) %>% NROW
#> [1] 1500
```

## Deduplicate


```r
smalldf <- sample_data_1[1:20, ]
# create a duplicate record
smalldf <- rbind(smalldf, smalldf[10,])
row.names(smalldf) <- NULL
# make it slightly different
smalldf[21, "key"] <- 1088954555
NROW(smalldf)
#> [1] 21
dp <- dframe(smalldf) %>% dedup()
NROW(dp)
#> [1] 20
attr(dp, "dups")
#> # A tibble: 1 x 5
#>   name             longitude latitude date                       key
#>   <chr>                <dbl>    <dbl> <dttm>                   <dbl>
#> 1 Ursus americanus     -76.8     35.5 2015-04-05 23:00:00 1088954555
```

## Dates

Standardize/convert dates


```r
# date_standardize(dframe(df), "%d%b%Y") # w/o pipe
dframe(sample_data_1) %>% date_standardize("%d%b%Y")
#> # A tibble: 1,500 x 5
#>    name             longitude latitude date             key
#>    <chr>                <dbl>    <dbl> <chr>          <int>
#>  1 Ursus americanus     -79.7     38.4 14Jan2015 1065590124
#>  2 Ursus americanus     -82.4     35.7 13Jan2015 1065588899
#>  3 Ursus americanus     -99.1     23.7 20Feb2015 1098894889
#>  4 Ursus americanus     -72.8     43.9 13Feb2015 1065611122
#>  5 Ursus americanus     -72.3     43.9 01Mar2015 1088908315
#>  6 Ursus americanus    -109.      32.7 29Mar2015 1088932238
#>  7 Ursus americanus    -109.      32.7 29Mar2015 1088932273
#>  8 Ursus americanus    -124.      40.1 28Mar2015 1132403409
#>  9 Ursus americanus     -78.3     36.9 20Mar2015 1088923534
#> 10 Ursus americanus     -76.8     35.5 05Apr2015 1088954559
#> # … with 1,490 more rows
```

Drop records without dates


```r
NROW(sample_data_1)
#> [1] 1500
NROW(dframe(sample_data_1) %>% date_missing())
#> [1] 1498
```

Create date field from other fields


```r
dframe(sample_data_2) %>% date_create(year, month, day)
#> # A tibble: 1,500 x 8
#>    name             longitude latitude        key year  month day   date      
#>    <chr>                <dbl>    <dbl>      <int> <chr> <chr> <chr> <chr>     
#>  1 Ursus americanus     -79.7     38.4 1065590124 2015  01    14    2015-01-14
#>  2 Ursus americanus     -82.4     35.7 1065588899 2015  01    13    2015-01-13
#>  3 Ursus americanus     -99.1     23.7 1098894889 2015  02    20    2015-02-20
#>  4 Ursus americanus     -72.8     43.9 1065611122 2015  02    13    2015-02-13
#>  5 Ursus americanus     -72.3     43.9 1088908315 2015  03    01    2015-03-01
#>  6 Ursus americanus    -109.      32.7 1088932238 2015  03    29    2015-03-29
#>  7 Ursus americanus    -109.      32.7 1088932273 2015  03    29    2015-03-29
#>  8 Ursus americanus    -124.      40.1 1132403409 2015  03    28    2015-03-28
#>  9 Ursus americanus     -78.3     36.9 1088923534 2015  03    20    2015-03-20
#> 10 Ursus americanus     -76.8     35.5 1088954559 2015  04    05    2015-04-05
#> # … with 1,490 more rows
```

## Taxonomy

Only one function exists for taxonomy cleaning, it removes rows where taxonomic names are 
either missing an epithet, or are missing altogether  (`NA` or `NULL`).

Get some data from GBIF, via `rgbif`


```r
if (requireNamespace("rgbif", quietly = TRUE)) {
  library("rgbif")
  res <- occ_data(limit = 500)$data
} else {
  res <- sample_data_3
}
```

Clean names


```r
NROW(res)
#> [1] 500
df <- dframe(res) %>% tax_no_epithet(name = "name")
NROW(df)
#> [1] 481
attr(df, "name_var")
#> [1] "name"
attr(df, "tax_no_epithet")
#> # A tibble: 19 x 107
#>    key   scientificName decimalLatitude decimalLongitude issues datasetKey
#>    <chr> <chr>                    <dbl>            <dbl> <chr>  <chr>     
#>  1 1637… Aves                     -34.5           136.   ""     40c0f670-…
#>  2 1989… Psychidae                 25.0           122.   "txma… e0b8cb67-…
#>  3 2542… Agaricales                55.9            12.3  "cdro… 84d26682-…
#>  4 2542… Corticiaceae              56.5             9.84 ""     84d26682-…
#>  5 2542… Corticiaceae              55.1            10.5  "cdro… 84d26682-…
#>  6 2542… Fungi                     56.5             9.84 "cdro… 84d26682-…
#>  7 2542… Agaricales                55.9            12.5  "cdro… 84d26682-…
#>  8 2542… Polyporales               55.9            12.3  "cdro… 84d26682-…
#>  9 2542… Trichiaceae               56.7             9.87 "cdro… 84d26682-…
#> 10 2542… Xylariales                56.2            10.6  "cdro… 84d26682-…
#> 11 2542… Corticiaceae              56.5             9.84 "cdro… 84d26682-…
#> 12 2542… Hymenochaetac…            56.5             9.84 "cdro… 84d26682-…
#> 13 2542… Polyporales               56.0            12.3  "cdro… 84d26682-…
#> 14 2542… Fungi                     55.8            12.5  "cdro… 84d26682-…
#> 15 2542… Fungi                     55.9            12.4  "cdro… 84d26682-…
#> 16 2542… Fungi                     55.9            12.3  "cdro… 84d26682-…
#> 17 2542… Fungi                     55.9            12.3  "cdro… 84d26682-…
#> 18 2542… Hyaloscyphace…            54.9            11.5  "cdro… 84d26682-…
#> 19 2542… Physarales                54.9            11.5  "cdro… 84d26682-…
#> # … with 101 more variables: publishingOrgKey <chr>, installationKey <chr>,
#> #   publishingCountry <chr>, protocol <chr>, lastCrawled <chr>,
#> #   lastParsed <chr>, crawlId <int>, basisOfRecord <chr>,
#> #   individualCount <int>, taxonKey <int>, kingdomKey <int>, phylumKey <int>,
#> #   classKey <int>, orderKey <int>, familyKey <int>, genusKey <int>,
#> #   speciesKey <int>, acceptedTaxonKey <int>, acceptedScientificName <chr>,
#> #   kingdom <chr>, phylum <chr>, order <chr>, family <chr>, genus <chr>,
#> #   species <chr>, genericName <chr>, specificEpithet <chr>, taxonRank <chr>,
#> #   taxonomicStatus <chr>, year <int>, month <int>, eventDate <chr>,
#> #   modified <chr>, lastInterpreted <chr>, references <chr>, license <chr>,
#> #   class <chr>, countryCode <chr>, recordedByIDs <list>,
#> #   identifiedByIDs <list>, rightsHolder <chr>, identifier <chr>,
#> #   nomenclaturalCode <chr>, dynamicProperties <chr>, language <chr>,
#> #   collectionCode <chr>, gbifID <chr>, occurrenceID <chr>, type <chr>,
#> #   taxonRemarks <chr>, preparations <chr>, recordedBy <chr>,
#> #   catalogNumber <chr>, vernacularName <chr>, institutionCode <chr>,
#> #   previousIdentifications <chr>, ownerInstitutionCode <chr>,
#> #   occurrenceRemarks <chr>, bibliographicCitation <chr>, accessRights <chr>,
#> #   higherClassification <chr>, dateIdentified <chr>, elevation <dbl>,
#> #   elevationAccuracy <dbl>, stateProvince <chr>, day <int>,
#> #   geodeticDatum <chr>, country <chr>, recordNumber <chr>, municipality <chr>,
#> #   locality <chr>, datasetName <chr>, identifiedBy <chr>, eventID <chr>,
#> #   occurrenceStatus <chr>, locationRemarks <chr>, dataGeneralizations <chr>,
#> #   taxonConceptID <chr>, coordinateUncertaintyInMeters <dbl>, lifeStage <chr>,
#> #   infraspecificEpithet <chr>, associatedReferences <chr>, county <chr>,
#> #   verbatimElevation <chr>, fieldNumber <chr>, continent <chr>,
#> #   identificationVerificationStatus <chr>, taxonID <chr>, eventTime <chr>,
#> #   behavior <chr>, informationWithheld <chr>, endDayOfYear <chr>,
#> #   originalNameUsage <chr>, startDayOfYear <chr>, datasetID <chr>,
#> #   habitat <chr>, associatedTaxa <chr>, locationAccordingTo <chr>,
#> #   locationID <chr>, verbatimLocality <chr>, …
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecoregion.R
\name{eco_region}
\alias{eco_region}
\alias{regions_meow}
\alias{regions_fao}
\title{Filter points within ecoregions}
\usage{
eco_region(x, dataset = "meow", region, lat = NULL, lon = NULL, drop = TRUE)

regions_meow()

regions_fao()
}
\arguments{
\item{x}{(data.frame) A data.frame}

\item{dataset}{(character) the dataset to use. one of: "meow" (Marine
Ecoregions of the World), "fao" (). See Details.}

\item{region}{(character) one or more region names. has the form \code{a:b} where
\code{a} is a variable name (column in the sf object) and \code{b} is the value you want
to filter to within that variable. See Details.}

\item{lat, lon}{(character) name of the latitude and longitude column to use}

\item{drop}{(logical) Drop bad data points or not. Either way, we parse out
bad data points as an attribute you can access. Default: \code{TRUE}
#param ignore.na (logical) To consider NA values as a bad point or not.
Default: \code{FALSE}}
}
\value{
Returns a data.frame, with attributes
}
\description{
Filter points within ecoregions
}
\details{
see \code{scrubr_cache} for managing the cache of data
}
\section{dataset options}{

\itemize{
\item Marine Ecoregions of the World (meow):
\itemize{
\item data from: https://opendata.arcgis.com/datasets/ed2be4cf8b7a451f84fd093c2e7660e3_0.geojson
}
\item Food and Agriculture Organization (fao):
\itemize{
\item data from: http://www.fao.org/geonetwork/srv/en/main.home?uuid=ac02a460-da52-11dc-9d70-0017f293bd28
}
}
}

\section{region options}{

\itemize{
\item within meow:
\itemize{
\item ECOREGION: many options, see \code{regions_meow()}
\item ECO_CODE: many options, see \code{regions_meow()}
\item and you can use others as well; run \code{regions_meow()} to get the data used
within \code{eco_region()} and see what variables/columns can be used
}
\item within fao:
\itemize{
\item OCEAN: Atlantic, Pacific, Indian, Arctic
\item SUBOCEAN: 1 through 11 (inclusive)
\item F_AREA (fishing area): 18, 21, 27, 31, 34, 37, 41, 47, 48, 51, 57, 58,
61, 67, 71, 77, 81, 87, 88
\item and you can use others as well; run \code{regions_fao()} to get the data used
within \code{eco_region()} and see what variables/columns can be used
}
}
}

\examples{
\dontrun{
if (requireNamespace("mapview") && requireNamespace("sf") && interactive()) {
## Marine Ecoregions of the World
wkt <- 'POLYGON((-119.8 12.2, -105.1 11.5, -106.1 21.6, -119.8 20.9, -119.8 12.2))'
res <- rgbif::occ_data(geometry = wkt, limit=300)$data
res2 <- sf::st_as_sf(res, coords = c("decimalLongitude", "decimalLatitude"))
res2 <- sf::st_set_crs(res2, 4326)
mapview::mapview(res2)
tmp <- eco_region(dframe(res), dataset = "meow",
   region = "ECOREGION:Mexican Tropical Pacific")
tmp2 <- sf::st_as_sf(tmp, coords = c("decimalLongitude", "decimalLatitude"))
tmp2 <- sf::st_set_crs(tmp2, 4326)
mapview::mapview(tmp2)
## filter many regions at once
out <- eco_region(dframe(res), dataset = "meow",
   region = c("ECOREGION:Mexican Tropical Pacific", "ECOREGION:Seychelles"))
out
out2 <- sf::st_as_sf(out, coords = c("decimalLongitude", "decimalLatitude"))
out2 <- sf::st_set_crs(out2, 4326)
mapview::mapview(out2)

## FAO
## FIXME - this needs fixing, broken
wkt <- 'POLYGON((72.2 38.5,-173.6 38.5,-173.6 -41.5,72.2 -41.5,72.2 38.5))'
manta_ray <- rgbif::name_backbone("Mobula alfredi")$usageKey
res <- rgbif::occ_data(manta_ray, geometry = wkt, limit=300, hasCoordinate = TRUE)
dat <- sf::st_as_sf(res$data, coords = c("decimalLongitude", "decimalLatitude"))
dat <- sf::st_set_crs(dat, 4326)
mapview::mapview(dat)
tmp <- eco_region(dframe(res$data), dataset = "fao", region = "OCEAN:Indian")
tmp <- tmp[!is.na(tmp$decimalLongitude), ]
tmp2 <- sf::st_as_sf(tmp, coords = c("decimalLongitude", "decimalLatitude"))
tmp2 <- sf::st_set_crs(tmp2, 4326)
mapview::mapview(tmp2)
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\name{scrubr_datasets}
\alias{scrubr_datasets}
\title{scrubr datasets}
\description{
\itemize{
\item \link{sample_data_1}
\item \link{sample_data_2}
\item \link{sample_data_3}
\item \link{sample_data_4}
\item \link{sample_data_6}
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\docType{data}
\name{sample_data_3}
\alias{sample_data_3}
\title{Sample data.frame number 3}
\format{
A data frame with 200 rows and 5 variables:
\describe{
\item{name}{taxonomic name}
\item{decimalLatitude}{latitude, decimal degree}
\item{decimalLongitude}{longitude, decimal degree}
\item{key}{occurrence key}
\item{eventDate}{date, date the occurrence was recorded}
}

Data originally collected from GBIF
}
\description{
Sample data.frame number 3
}
\keyword{datasets}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coord-funs.R
\name{coords}
\alias{coords}
\alias{coord_incomplete}
\alias{coord_imprecise}
\alias{coord_impossible}
\alias{coord_unlikely}
\alias{coord_within}
\alias{coord_pol_centroids}
\alias{coord_uncertain}
\title{Coordinate based cleaning}
\usage{
coord_incomplete(x, lat = NULL, lon = NULL, drop = TRUE)

coord_imprecise(x, which = "both", lat = NULL, lon = NULL, drop = TRUE)

coord_impossible(x, lat = NULL, lon = NULL, drop = TRUE)

coord_unlikely(x, lat = NULL, lon = NULL, drop = TRUE)

coord_within(
  x,
  field = NULL,
  country = NULL,
  lat = NULL,
  lon = NULL,
  drop = TRUE
)

coord_pol_centroids(x, lat = NULL, lon = NULL, drop = TRUE)

coord_uncertain(
  x,
  coorduncertainityLimit = 30000,
  drop = TRUE,
  ignore.na = FALSE
)
}
\arguments{
\item{x}{(data.frame) A data.frame}

\item{lat, lon}{(character) Latitude and longitude column to use. See Details.}

\item{drop}{(logical) Drop bad data points or not. Either way, we parse out
bad data points as an attribute you can access. Default: \code{TRUE}}

\item{which}{(character) one of "has_dec", "no_zeros", or "both" (default)}

\item{field}{(character) Name of field in input data.frame x with country
names}

\item{country}{(character) A single country name}

\item{coorduncertainityLimit}{(numeric) numeric threshold for the
coordinateUncertainityInMeters variable. Default: 30000}

\item{ignore.na}{(logical) To consider NA values as a bad point or not.
Default: \code{FALSE}}
}
\value{
Returns a data.frame, with attributes
}
\description{
Coordinate based cleaning
}
\details{
Explanation of the functions:
\itemize{
\item coord_impossible - Impossible coordinates
\item coord_incomplete - Incomplete coordinates
\item coord_imprecise - Imprecise coordinates
\item coord_pol_centroids - Points at political centroids
\item coord_unlikely - Unlikely coordinates
\item coord_within - Filter points within user input
political boundaries
\item coord_uncertain - Uncertain occurrances of measured through
coordinateUncertaintyInMeters default limit= 30000
}

If either lat or lon (or both) given, we assign the given column name
to be standardized names of "latitude", and "longitude". If not given, we
attempt to guess what the lat and lon column names are and assign the
same standardized names. Assigning the same standardized names makes
downstream processing easier so that we're dealing with consistent column
names. On returning the data, we return the original names.

For \code{coord_within}, we use \code{countriesLow} dataset from the
\pkg{rworldmap} package to get country borders.
}
\section{coord_pol_centroids}{

Right now, this function only deals with city centroids, using the
\link[maps:world.cities]{maps::world.cities} dataset of more than 40,000 cities.
We'll work on adding country centroids, and perhaps others (e.g.,
counties, states, provinces, parks, etc.).
}

\examples{
df <- sample_data_1

# Remove impossible coordinates
NROW(df)
df[1, "latitude"] <- 170
df <- dframe(df) \%>\% coord_impossible()
NROW(df)
attr(df, "coord_impossible")

# Remove incomplete cases
NROW(df)
df_inc <- dframe(df) \%>\% coord_incomplete()
NROW(df_inc)
attr(df_inc, "coord_incomplete")


# Remove imprecise cases
df <- sample_data_5
NROW(df)
## remove records that don't have decimals at all
df_imp <- dframe(df) \%>\% coord_imprecise(which = "has_dec")
NROW(df_imp)
attr(df_imp, "coord_imprecise")
## remove records that have all zeros
df_imp <- dframe(df) \%>\% coord_imprecise(which = "no_zeros")
NROW(df_imp)
attr(df_imp, "coord_imprecise")
## remove both records that don't have decimals at all and those that
## have all zeros
df_imp <- dframe(df) \%>\% coord_imprecise(which = "both")
NROW(df_imp)
attr(df_imp, "coord_imprecise")


# Remove unlikely points
NROW(df)
df_unlikely <- dframe(df) \%>\% coord_unlikely()
NROW(df_unlikely)
attr(df_unlikely, "coord_unlikely")

# Remove points not within correct political borders
if (requireNamespace("rgbif", quietly = TRUE) && interactive()) {
   library("rgbif")
   wkt <- 'POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))'
   res <- rgbif::occ_data(geometry = wkt, limit=300)$data
} else {
   res <- sample_data_4
}

## By specific country name
if (
  interactive() &&
  requireNamespace("sf", quietly=TRUE) && 
  requireNamespace("s2", quietly=TRUE) && 
  requireNamespace("rworldmap", quietly=TRUE)
) {
NROW(res)
df_within <- dframe(res) \%>\% coord_within(country = "Israel")
NROW(df_within)
attr(df_within, "coord_within")

## By a field in your data - makes sure your points occur in one
## of those countries
NROW(res)
df_within <- dframe(res) \%>\% coord_within(field = "country")
NROW(df_within)
head(df_within)
attr(df_within, "coord_within")
}

# Remove those very near political centroids
## not ready yet
# NROW(df)
# df_polcent <- dframe(df) \%>\% coord_pol_centroids()
# NROW(df_polcent)
# attr(df_polcent, "coord_polcent")

## lat/long column names can vary
df <- sample_data_1
head(df)
names(df)[2:3] <- c('mylon', 'mylat')
head(df)
df[1, "mylat"] <- 170
dframe(df) \%>\% coord_impossible(lat = "mylat", lon = "mylon")

df <- sample_data_6

# Remove uncertain occurances

NROW(df)
df1<-df \%>\% coord_uncertain()
NROW(df1)
attr(df, "coord_uncertain")

NROW(df)
df2<-df \%>\% coord_uncertain(coorduncertainityLimit = 20000)
NROW(df2)

NROW(df)
df3<-df \%>\% coord_uncertain(coorduncertainityLimit = 20000,ignore.na=TRUE)
NROW(df3)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_df.R
\name{dframe}
\alias{dframe}
\title{Compact data.frame}
\usage{
dframe(x)
}
\arguments{
\item{x}{Input data.frame}
}
\description{
Compact data.frame
}
\examples{
dframe(sample_data_1)
dframe(mtcars)
dframe(iris)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\docType{data}
\name{sample_data_1}
\alias{sample_data_1}
\title{Sample data.frame 1}
\format{
A data frame with 1500 rows and 5 variables:
\describe{
\item{name}{taxonomic name}
\item{longitude}{longitude, decimal degree}
\item{latitude}{latitude, decimal degree}
\item{date}{date, date the occurrence was recorded}
\item{key}{occurrence key}
}

Data originally collected from GBIF
}
\description{
Sample data.frame 1
}
\keyword{datasets}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\docType{data}
\name{sample_data_4}
\alias{sample_data_4}
\title{Sample data.frame number 4}
\format{
A data frame with 100 rows and 6 variables:
\describe{
\item{name}{taxonomic name}
\item{decimalLatitude}{latitude, decimal degree}
\item{decimalLongitude}{longitude, decimal degree}
\item{key}{occurrence key}
\item{eventDate}{date, date the occurrence was recorded}
\item{country}{country of occurrence record}
}

Data originally collected from GBIF
}
\description{
Sample data.frame number 4
}
\keyword{datasets}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collector-funs.R
\name{collectors}
\alias{collectors}
\alias{coll_clean}
\title{Collector based cleaning}
\usage{
coll_clean(x, collector = NULL)
}
\arguments{
\item{x}{(data.frame) A data.frame}

\item{collector}{(character) Collector field to use. See Details.}

\item{drop}{(logical) Drop bad data points or not. Either way, we parse
out bade data points as an attribute you can access. Default: \code{TRUE}}
}
\value{
Returns a data.frame, with attributes
}
\description{
Collector based cleaning
}
\details{
Explanation of the functions:
\itemize{
\item coll_clean - Standardize collector names
}
}
\examples{
# df <- data.frame(
#   coll = c('K.F.P. Martius', 'C. F. P. Martius', 'C. F. P. von Martius'),
#   species = 'Poa annua',
#   lat = 1:3,
#   lon = 4:6,
#  stringsAsFactors = FALSE
# )

# Standardize names
# NROW(df)
# df <- dframe(df) \%>\% coll_clean()
# NROW(df)
# attr(df, "coll_clean")
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\docType{data}
\name{sample_data_7}
\alias{sample_data_7}
\title{Sample data.frame number 7}
\format{
A data frame with 50 rows and 91 variables

Data originally collected from GBIF
}
\description{
Sample data.frame number 7
}
\keyword{datasets}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harvard-botanists.R
\name{bot_search}
\alias{bot_search}
\title{Harvard botanist index functions}
\usage{
bot_search(
  name = NULL,
  individual = FALSE,
  start = NULL,
  fuzzy = FALSE,
  remarks = NULL,
  speciality = NULL,
  country = NULL,
  is_collector = FALSE,
  is_author = FALSE,
  team = FALSE,
  error = stop,
  ...
)
}
\description{
Harvard botanist index functions
}
\examples{
\dontrun{
# bot_search(name = "Asa Gray")
# bot_search(name = "A. Gray")
# bot_search(remarks = "harvard")
# bot_search(name = "Gray", fuzzy = TRUE)

## FIXME - this leads to a JSON parsing error because they give
##   bad JSON in some results, including this example
# bot_search(country = "China")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy-funs.R
\name{taxonomy}
\alias{taxonomy}
\alias{tax_no_epithet}
\title{Taxonomy based cleaning}
\usage{
tax_no_epithet(x, name = NULL, drop = TRUE)
}
\arguments{
\item{x}{(data.frame) A data.frame1}

\item{name}{(character) Taxonomic name field Optional. See Details.}

\item{drop}{(logical) Drop bad data points or not. Either way, we parse
out bade data points as an attribute you can access. Default: \code{TRUE}}
}
\value{
Returns a data.frame, with attributes
}
\description{
Taxonomy based cleaning
}
\examples{
if (requireNamespace("rgbif", quietly = TRUE) && interactive()) {
   library("rgbif")
   res <- rgbif::occ_data(limit = 200)$data
} else {
   res <- sample_data_3
}

# Remove records where names don't have genus + epithet
## so removes those with only genus and those with no name (NA or NULL)
NROW(res)
df <- dframe(res) \%>\% tax_no_epithet(name = "name")
NROW(df)
attr(df, "name_var")
attr(df, "tax_no_epithet")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
Pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr-package.R
\docType{package}
\name{scrubr-package}
\alias{scrubr-package}
\alias{scrubr}
\title{scrubr}
\description{
Clean biological occurrence data
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\docType{data}
\name{sample_data_6}
\alias{sample_data_6}
\title{Sample data.frame number 6}
\format{
A data frame with 50 rows and 5 variables:
\describe{
\item{name}{taxonomic name}
\item{key}{GBIF occurrence key}
\item{decimalLatitude}{latitude, decimal degree}
\item{decimalLongitude}{longitude, decimal degree}
\item{coordinateUncertaintyInMeters}{Uncertainity, the point-radius representation of the location}
}

Data originally collected from GBIF
}
\description{
Sample data.frame number 6
}
\keyword{datasets}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dedup.R
\name{dedup}
\alias{dedup}
\title{Deduplicate records}
\usage{
dedup(x, how = "one", tolerance = 0.9)
}
\arguments{
\item{x}{(data.frame) A data.frame, tibble, or data.table}

\item{how}{(character) How to deal with duplicates. The default of
"one" keeps one record of each group of duplicates, and drops the
others, putting them into the \code{dups} attribute. "all" drops all
duplicates, in case e.g., you don't want to deal with any records that are
duplicated, as e.g., it may be hard to tell which one to remove.}

\item{tolerance}{(numeric) Score (0 to 1) at which to determine a match.
You'll want to inspect outputs closely to tweak this value based on your
data, as results can vary.}
}
\value{
Returns a data.frame, optionally with attributes
}
\description{
Deduplicate records
}
\examples{
df <- sample_data_1
smalldf <- df[1:20, ]
smalldf <- rbind(smalldf, smalldf[10,])
smalldf[21, "key"] <- 1088954555
NROW(smalldf)
dp <- dframe(smalldf) \%>\% dedup()
NROW(dp)
attr(dp, "dups")

# Another example - more than one set of duplicates
df <- sample_data_1
twodups <- df[1:10, ]
twodups <- rbind(twodups, twodups[c(9, 10), ])
rownames(twodups) <- NULL
NROW(twodups)
dp <- dframe(twodups) \%>\% dedup()
NROW(dp)
attr(dp, "dups")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fix_names.R
\name{fix_names}
\alias{fix_names}
\title{Change taxonomic names to be the same for each taxon}
\usage{
fix_names(x, how = "shortest", replace = NULL)
}
\arguments{
\item{x}{(data.frame) A data.frame. the target taxonomic name
column should be 'name'}

\item{how}{One of a few different methods:
\itemize{
\item shortest - Takes the shortest name string that is likely to be the
prettiest to display name, and replaces alll names with that one, better
for maps, etc.
\item supplied - If this method, supply a vector of names to replace the
names with.
}}

\item{replace}{A data.frame of names to replace names in the occurrence
data.frames with. Only used if how="supplied". The data.frame should have
two columns: the first is the names to match in the input \code{x} data.frame,
and the second column is the name to replace with. The column names don't
matter.}
}
\value{
a data.frame
}
\description{
That is, this function attempts to take all the names that are synonyms,
for whatever reason (e.g., some names have authorities on them), and
collapses them to the same string - making data easier to deal with for
making maps, etc. OR - you can think of this as a tool for
}
\examples{
\dontrun{
df <- sample_data_7

# method: shortest
fix_names(df, how="shortest")$name

# method: supplied
(replace_df <- data.frame(
 one = unique(df$name), 
 two = c('P. contorta', 'P.c. var. contorta',
         'P.c. subsp bolanderi', 'P.c. var. murrayana'),
 stringsAsFactors = FALSE))
fix_names(df, how="supplied", replace = replace_df)$name
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\docType{data}
\name{sample_data_2}
\alias{sample_data_2}
\title{Sample data.frame number 2}
\format{
A data frame with 1500 rows and 7 variables:
\describe{
\item{name}{taxonomic name}
\item{longitude}{longitude, decimal degree}
\item{latitude}{latitude, decimal degree}
\item{key}{occurrence key}
\item{year}{year of occurrence}
\item{month}{month of occurrence}
\item{day}{day of occurrence}
}

Data originally collected from GBIF
}
\description{
Sample data.frame number 2
}
\keyword{datasets}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/date-funs.R
\name{date}
\alias{date}
\alias{date_standardize}
\alias{date_missing}
\alias{date_create}
\alias{date_create_}
\title{Date based cleaning}
\usage{
date_standardize(x, format = "\%Y-\%m-\%d", date_column = "date", ...)

date_missing(x, date_column = "date", drop = TRUE, ...)

date_create(x, ...)

date_create_(x, ..., .dots, format = "\%Y-\%m-\%d", date_column = "date")
}
\arguments{
\item{x}{(data.frame) A data.frame}

\item{format}{(character) Date format. See \code{\link[=as.Date]{as.Date()}}}

\item{date_column}{(character) Name of the date column}

\item{...}{Comma separated list of unquoted variable names}

\item{drop}{(logical) Drop bad data points or not. Either way, we parse
out bade data points as an attribute you can access. Default: \code{TRUE}}

\item{.dots}{Used to work around non-standard evaluation}
}
\value{
Returns a data.frame, with attributes
}
\description{
Date based cleaning
}
\details{
\itemize{
\item date_standardize - Converts dates to a specific format
\item date_missing - Drops records that do not have dates, either via being
NA or being a zero length character string
\item date_create - Create a date field from
}
}
\examples{
df <- sample_data_1
# Standardize dates
dframe(df) \%>\% date_standardize()
dframe(df) \%>\% date_standardize("\%Y/\%m/\%d")
dframe(df) \%>\% date_standardize("\%d\%b\%Y")
dframe(df) \%>\% date_standardize("\%Y")
dframe(df) \%>\% date_standardize("\%y")

# drop records without dates
NROW(df)
NROW(dframe(df) \%>\% date_missing())

# Create date field from other fields
df <- sample_data_2
## NSE
dframe(df) \%>\% date_create(year, month, day)
## SE
date_create_(dframe(df), "year", "month", "day")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrubr_datasets.R
\docType{data}
\name{sample_data_5}
\alias{sample_data_5}
\title{Sample data.frame number 5}
\format{
A data frame with 39 rows and 5 variables:
\describe{
\item{name}{taxonomic name}
\item{longitude}{longitude, decimal degree}
\item{latitude}{latitude, decimal degree}
\item{date}{date, date the occurrence was recorded}
\item{key}{GBIF occurrence key}
}

Data originally collected from GBIF
}
\description{
Sample data.frame number 5
}
\keyword{datasets}
\keyword{internal}
