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
