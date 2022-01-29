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

*Wow, no problems at all. :)**Wow, no problems at all. :)*