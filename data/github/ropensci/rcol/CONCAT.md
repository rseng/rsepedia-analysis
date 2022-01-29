rcol
====



<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rcol)](https://cranchecks.info/pkgs/rcol)
[![R-CMD-check](https://github.com/ropensci/rcol/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rcol/actions)
[![codecov](https://codecov.io/gh/ropensci/rcol/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rcol)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcol)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcol)](https://cran.r-project.org/package=rcol)

`rcol` is a R client for the Catalogue of Life

Package documentation: https://docs.ropensci.org/rcol/

* Catalogue of Life: http://www.catalogueoflife.org/
* Catalogue of Life GitHub repository: https://github.com/CatalogueOfLife/general
* COL API docs: https://api.catalogueoflife.org/
* Web portal for COL: https://data.catalogueoflife.org/

## Installation


```r
install.packages("rcol")
```

Dev version


```r
pak::pkg_install("ropensci/rcol")
# OR
install.packages("rcol", repos="https://dev.ropensci.org")
```


```r
library("rcol")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rcol/issues).
* License: MIT
* Get citation information for `rcol` in R doing `citation(package = 'rcol')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rcol 0.2.0
==========

### MINOR IMPROVEMENTS

* fix broken example: required a fix in an internal funciton to sort data.frame columns that wasn't robust to missing columns (#8)
* fixed tests

rcol 0.1.0
==========

### NEW FEATURES

* First submission to CRAN.
## Test environments
* local R installation on macOS, R 4.1.0
* ubuntu 16.04 (on github actions), R 4.1.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

-----

This submission fixes broken examples.

Thanks!
Scott Chamberlain
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

* Submit an issue on the [Issues page](https://github.com/ropensci/rcol/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rcol.git`
* Make sure to track progress upstream (i.e., on our version of `rcol` at `ropensci/rcol`) by doing `git remote add upstream https://github.com/ropensci/rcol.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/rcol`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
