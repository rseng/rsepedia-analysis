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

*Wow, no problems at all. :)**Wow, no problems at all. :)*