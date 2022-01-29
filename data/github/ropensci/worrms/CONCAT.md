worrms
======



<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/worrms)](https://cranchecks.info/pkgs/worrms)
[![R-check](https://github.com/ropensci/worrms/workflows/R-check/badge.svg)](https://github.com/ropensci/worrms/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/worrms/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/worrms)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/worrms)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/worrms)](https://cran.r-project.org/package=worrms)

`worrms` is a R client for the World Register of Marine Species

* World Register of Marine Species (WoRMS) http://www.marinespecies.org/
* WoRMS REST API docs: http://www.marinespecies.org/rest/

See the taxize book (https://taxize.dev) for taxonomically focused work
in this and similar packages.

## Installation

More stable CRAN version


```r
install.packages("worrms")
```

Development version


```r
remotes::install_github("ropensci/worrms")
```


```r
library("worrms")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/worrms/issues).
* License: MIT
* Get citation information for `worrms` in R doing `citation(package = 'worrms')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
worrms 0.4.2
============

### MINOR IMPROVEMENTS

* fix a few failing tests on cran (#22)


worrms 0.4.0
============

### NEW FEATURES

* new functions `wm_ranks_id()` and `wm_ranks_name()` for getting taxonomic ranks by rank identifier or rank name (#20)
* new function `wm_records_rank()` for getting AphiaRecords for a given rank id (#20)

### MINOR IMPROVEMENTS

* `wm_synonyms()` gains `offset` parameter to allow pagination (#20)
* `tibble::as_data_frame()` replaced with `tibble::as_tibble()`

### DEPRECATED AND DEFUNCT

* `wm_record_()` is deprecated; `wm_record()` now handles 1 or more AphiaID's

### BUG FIXES

* fix `wm_children` test that was failing on cran checks (#21)


worrms 0.3.2
============

### MINOR IMPROVEMENTS

* add link to taxize book in vignette and README (#12)

### BUG FIXES

* fix bug in test regarding date (#19)

worrms 0.3.0
============

### MINOR IMPROVEMENTS

* fix to most functions throughout the package (those that have two versions, with and without an underscore): underscore versions of functions now do not error when an input is not found, but instead warn the user and move on - to facilitate working with many inputs. the non-underscore version of each function still only accepts 1 input and errors if you give more than 1 (#14) (#18)

### BUG FIXES

* make sure that functions that accept only 1 input for the first parameter error well with an informative message (#15)


worrms 0.2.8
============

### NEW FEATURES

* Integration with `vcr` and `webmockr` packages for unit test stubbing
* gains new functions for getting WORMS traits data (they call them "attributes"): `wm_attr_aphia`, `wm_attr_aphia_`, `wm_attr_category`, `wm_attr_category_`, `wm_attr_data`, `wm_attr_data_`, `wm_attr_def`, `wm_attr_def_`  (#3)


worrms 0.2.0
============

### NEW FEATURES

* Added additional sister functions to most exported functions in the 
package, all with trailing underscore. For example, `wm_children` and 
`wm_children_`. These underscore methods take in many inputs, typically
of a AphiaID or a taxonomic or vernacular name. We decided to make 
separate functions so that we minimize any disturbance to the existing 
package API. (#4) (#6)

### MINOR IMPROVEMENTS

* Moved to using markdown docs (#5)
* All functions now state what they return (#9)


worrms 0.1.0
============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 reverse dependency.
  (Summary at <https://github.com/ropensci/worrms/blob/master/revdep/README.md>). No problems were found.

---

This version fixes some failing cran tests. 

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

* Submit an issue on the [Issues page](https://github.com/ropensci/worrms/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/worrms.git`
* Make sure to track progress upstream (i.e., on our version of `worrms` at `ropensci/worrms`) by doing `git remote add upstream https://github.com/ropensci/worrms.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/worrms`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.2 (2020-06-22) |
|os       |macOS Catalina 10.15.5       |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |US/Pacific                   |
|date     |2020-07-07                   |

# Dependencies

|package |old   |new        |Δ  |
|:-------|:-----|:----------|:--|
|worrms  |0.4.0 |0.4.2      |*  |
|glue    |NA    |1.4.1.9000 |*  |
|Rcpp    |NA    |1.0.4.6    |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*