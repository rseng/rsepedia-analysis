ritis
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/ritis)](https://cranchecks.info/pkgs/ritis)
[![R-check](https://github.com/ropensci/ritis/workflows/R-check/badge.svg)](https://github.com/ropensci/ritis/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/ritis/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/ritis)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ritis)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/ritis)](https://cran.r-project.org/package=ritis)

An interface to the Integrated Taxonomic Information System (ITIS)

* ITIS API Docs: <https://www.itis.gov/ws_description.html>
* Solr service: <https://www.itis.gov/solr_documentation.html>
* taxize book: <https://taxize.dev/>
* ritis docs: <https://docs.ropensci.org/ritis/>

How to cite ITIS. From <https://itis.gov/citation.html>

To cite data obtained from ITIS, the following citation format is offered as a suggestion:

    Retrieved [month, day, year], from the Integrated Taxonomic Information System on-line database, http://www.itis.gov.


ITIS is one of many different taxonomic data sources. Other include: Catalogue of Life (and COL+), NCBI taxonomy, International Plant Names Index, Index Fungorum, and more. The Wikipedia entry (<https://en.wikipedia.org/wiki/Integrated_Taxonomic_Information_System>) states that ITIS has a North American focus, but includes many taxa not in North America.

## Terminology

* "mononomial": a taxonomic name with one part, e.g, _Poa_
* "binomial": a taxonomic name with two parts, e.g, _Poa annua_
* "trinomial": a taxonomic name with three parts, e.g, _Poa annua annua_

## Installation

Stable, CRAN version


```r
install.packages("ritis")
```

Dev version


```r
remotes::install_github("ropensci/ritis")
```


```r
library("ritis")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/ritis/issues).
* License: MIT
* Get citation information for `ritis` in R doing `citation(package = 'ritis')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
ritis 1.0
=========

### MINOR IMPROVEMENTS

* reduce code duplication in `terms()` (#18)
* add more unit tests (#21)

### BUG FIXES

* fix for `publications()` function: parsing error fixed (#19) (#20)


ritis 0.9.0
===========

### MINOR IMPROVEMENTS

* an example added to `itis_search()` for how to do a search for Class Aves, and how to drill down from Class Aves to genera within Aves; added text to readme and vignette about how to cite ITIS and a brief comparison of ITIS to other taxonomic data sources; added brief terminology section to readme and vignette with 3 terms thus far (mononomial, binomial, trinomial) (#16) thanks to @TrashBirdEcology for the prompt
* change `tibble::data_frame` useage to `tibble::tibble` (#17)


ritis 0.8.0
===========

### MINOR IMPROVEMENTS

* updated docs and examples for `itis_search()` to demonstate how to search appropriately with spaces and other characters  (#14)

### BUG FIXES

* `itis_group()` was failing on a parsing error (after retrieving the payload), via an error in parsing in `solrium` package; fixed now (#15)


ritis 0.7.6
===========

### MINOR IMPROVEMENTS

* improve docs for solr functions, pointing to appropriate docs in `solrium` package (#12)
* give link to taxize book in readme, vignette, and pkg level manual file (#13)

### BUG FIXES

* fixed bug in `search_anymatch()`: we weren't correctly handling cases where no results were returned (#11)

Full diff: https://github.com/ropensci/ritis/compare/v0.7.2...v0.7.6


ritis 0.7.2
===========

### NEW FEATURES

* Integration with `vcr` and `webmockr` packages for unit test stubbing


ritis 0.7.0
===========

### NEW FEATURES

* Now using new version of `solrium` package - users shouldn't
see any differences (#9)


ritis 0.6.0
===========

### NEW FEATURES

* Now using `crul` as underlying HTTP client (#5)

### BUG FIXES

* Base URL change for Solr service from `http` to `https` (#8)
* Fixed JSON parsing problem (#6)


ritis 0.5.4
===========

### BUG FIXES

* Base URL changed from `http` to `https`, was causing problems in some
requests, but not others. Changed to `https` (#4)


ritis 0.5.0
===========

### NEW FEATURES

* Re-released to CRAN
* Complete overhaul of the package API, simplifying all function
interfaces, using JSON by default, shorter names, reduce code reuse.
* Added functions for interacting with ITIS's new Solr
interface via use of `solrium`


ritis 0.0.3
===========

### BUG FIXES

* Removed dependency on plyr - moved from laply to lapply across functions.


ritis 0.0.2
===========

### BUG FIXES

* Temporarily removed all tests until they can be fixed and updated, and so that package passes checks.


ritis 0.0.1
===========

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.3 patched
* ubuntu 16.04 (on github actions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

I have run R CMD check on the 2 downstream dependencies.
There were no problems related to this package. Summary at
<https://github.com/ropensci/ritis/blob/master/revdep/README.md>

---

This version xxxx.

Thanks!
Scott Chamberlain
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

* Submit an issue on the [Issues page](https://github.com/ropensci/ritis/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/ritis.git`
* Make sure to track progress upstream (i.e., on our version of `ritis` at `ropensci/ritis`) by doing `git remote add upstream https://github.com/ropensci/ritis.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/ritis`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.3 Patched (2021-01-08 r79819) |
|os       |macOS Catalina 10.15.7                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-02-01                                  |

# Dependencies

|package |old   |new        |Δ  |
|:-------|:-----|:----------|:--|
|ritis   |0.9.0 |1.0.0      |*  |
|crayon  |NA    |1.4.0.9000 |*  |

# Revdeps

*Wow, no problems at all. :)*# Check times

|package  |version | check_time|
|:--------|:-------|----------:|
|camtrapR |0.99.9  |      124.6|
|taxize   |0.9.0   |       63.3|


*Wow, no problems at all. :)*