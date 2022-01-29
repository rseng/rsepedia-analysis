mapr
====



[![R-check](https://github.com/ropensci/mapr/workflows/R-check/badge.svg)](https://github.com/ropensci/mapr/actions/)
[![cran checks](https://cranchecks.info/badges/worst/mapr)](https://cranchecks.info/pkgs/mapr)
[![codecov](https://codecov.io/gh/ropensci/mapr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/mapr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/mapr?color=FAB657)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/mapr)](https://cran.r-project.org/package=mapr)


Helper for making maps of species occurrence data, including for:

* spocc (https://github.com/ropensci/spocc)
* rgbif (https://github.com/ropensci/rgbif)
* and some `sp` classes

This package has utilities for making maps with:

* base R
* ggplot2
* Leaflet - via `leaflet` pkg
* GitHub Gists - via `gistr` package

Get started with the docs: https://docs.ropensci.org/mapr/

## Installation

Install `mapr`


```r
install.packages("mapr")
```

Or the development version from GitHub


```r
remotes::install_github("ropensci/mapr")
```


```r
library("mapr")
library("spocc")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/mapr/issues).
* License: MIT
* Get citation information for `mapr` in R doing `citation(package = 'mapr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
mapr 0.5.2
==========

### MINOR IMPROVEMENTS

* fix vignette issue (#40)


mapr 0.5.0
==========

### DEFUNCT

* `map_ggmap()` is defunct. authentication has changed, more trouble than its worth (#34)

### NEW FEATURES

* all mapping functions gain the `name` parameter to specify the column that holds the taxon name - if not given, we look for a column "name" (#32)

### MINOR IMPROVEMENTS

* fix vignette missing title (#37)

### BUG FIXES

* fix non-ascii strings in the two package datasets - and script added to make those datasets reproducible, including fixing non-ascii strings (#39)
* remove linked references to pkgs in docs that are not imported/suggested (#38)
* `map_plot()` speed up: using `maps::map()` instead of `rworldmap::getMap()`, faster (#35)
* improve internal handling of `name` parameter users can pass down through mapping functions (#36)
* `rgbif` added to Suggests - was used in examples but wasn't in Suggests - used conditionally in examples now


mapr 0.4.0
==========

### MINOR IMPROVEMENTS

* All `map_*()` functions now support `gbif_data` class from the `rgbif` package, which is the output of `rgbif::occ_data()` (#29)
* Added `.github` files for contributing, issue and PR templates (#30)


mapr 0.3.4
==========

### MINOR IMPROVEMENTS

* Now using markdown for docs (#27)
* Replaced `httr` with `crul` as http client (#26)

### Problem with ggmap

* Note that there is a problem with `map_ggmap` due to a bug in 
`ggmap`. It is fixed in the `ggmap` dev version, so should be fixed
in the CRAN version soon, hopefully.


mapr 0.3.0
==========

### NEW FEATURES

* Now in all functions, when there's more than 1 taxon, we'll do a separate
color for each taxon and draw a legend if applicable (#21) (#22)
* Added support for adding convex hulls to some of the plot types (#23)
thanks to @rossmounce for the feature request
* `map_leaflet()` now adds metadata as a popup to each marker (#18) (#25)


mapr 0.2.0
==========

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.3 RC
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes
     
## Reverse dependencies

There are no reverse dependencies.

---

This version includes a fix to the vignette that was causing a warning on two CRAN platforms.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/mapr/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/mapr.git`
* Make sure to track progress upstream (i.e., on our version of `mapr` at `ropensci/mapr`) by doing `git remote add upstream https://github.com/ropensci/mapr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/mapr`

### Vignette changes

If you want to contribute to the vignette, or add a new one, those are kept in the `inst/vign/` directory. The main vignette is in `inst/vign/mapr_vignette.Rmd`.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```

</details>
