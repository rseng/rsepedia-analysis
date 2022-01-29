rredlist
========



[![cran checks](https://cranchecks.info/badges/worst/rredlist)](https://cranchecks.info/pkgs/rredlist)
[![R-check](https://github.com/ropensci/rredlist/workflows/R-check/badge.svg)](https://github.com/ropensci/rredlist/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/rredlist/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rredlist?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rredlist)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rredlist)](https://cran.r-project.org/package=rredlist)

`rredlist` is an R client for the IUCN Red List (https://apiv3.iucnredlist.org/api/v3/docs). The IUCN Red List is a global list of threatened and endangered species. IUCN Red List docs: https://apiv3.iucnredlist.org/api/v3/docs

## Installation

CRAN


```r
install.packages("rredlist")
```

Development version


```r
remotes::install_github("ropensci/rredlist")
# OR
install.packages("rredlist", repos="https://dev.ropensci.org")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rredlist/issues).
* License: MIT
* Get citation information for `rredlist` in R doing `citation(package = 'rredlist')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[token]: https://apiv3.iucnredlist.org/api/v3/token
[redlistr]: https://github.com/red-list-ecosystem/redlistr
rredlist 0.7.0
===================

### MINOR IMPROVEMENTS

* vignette added, but only available on the docs site (#24)
* when testing, if a iucm redlist key not found, set a dummy key (#41)
* readme improvements (#42)
* change base url for Red List API to https from http

rredlist 0.6.0
===================

### MINOR IMPROVEMENTS

* note in docs about how result may differ in website vs. in this package through the API  (#35)
* fail with useful message when NA's passed to parameters in package functions (#38)


rredlist 0.5.0
===================

### NEW FEATURES 

* gains new function `rl_use_iucn` to help with API key setup (#31) by @maelle
* gains new functions `rl_comp_groups` and `rl_comp_groups_` to interface with the comprehensive groups API route (#26)
* `rl_sp` gains two new parameters: `all` (logical) to toggle getting all results or not, if selected we do paging internally; `quiet` parameter (logical) suppresses progress (#29)

### MINOR IMPROVEMENTS

* mention `redlistr` package in README to help users decide which package to use for which use cases (#30)
* now using `webmockr` and `vcr` to do unit test caching (#33) (#34)



rredlist 0.4.0
==============

### NEW FEATURES

* Gains new functions `rl_growth_forms()` and `rl_growth_forms_()`. added 
tests for them as well (#20) thanks @stevenpbachman

### MINOR IMPROVEMENTS

* Now using markdown documentation (#22)
* Fixed many man files which for `region` parameter described 
requiring a taxonomic name - fixed to describe accurately. Also 
improved docs in general (#21)
* Added the options for `category` parameter in `rl_sp_category()` function 
* Added in docs for `rl_sp_country` how to get acceptable country codes to 
pass to `country` parameter
* Added to package level manual file `?rredlist-package` a note from the 
IUCN Redlist API documentation about that they suggest using taxonomic 
names instead of IDs because IDs can change through time



rredlist 0.3.0
==============

### NEW FEATURES

* New functions `rl_occ_country` and `rl_occ_country_` for 
getting country occurrences by species name or ID (#13)
* Replaced `httr` with `crul`. Please note this only affects use 
of curl options. See `crul` docs for how to use curl options (#14)

### MINOR IMPROVEMENTS

* User agent string like `r-curl/2.3 crul/0.2.0 rOpenSci(rredlist/0.3.0)` 
sent in all requests now to help IUCN API maintainers know 
how often requests come from R and this package (#19)
* Taxon names are now given back in `rl_threats` - we didn't do 
anything in the package - the API now gives the names back and 
adds them in a column (#10)
* Type checking all parameter inputs now both in terms of class
and length - with helpful error messages on fail (#17)
* Simplify package codebase by having single internal function for a 
suite of half a dozen or so functions that have similar pattern (#18)
* Removed `key` parameter from `rl_version()` and `rl_citation()` as
API key not required for those methods
* More thorough test suite


rredlist 0.2.0
==============

### NEW FEATURES

* New methods added to get historical assessments: `rl_history()`
and `rl_history_()` (#8)

### MINOR IMPROVEMENTS

* Fixed description of what `rl_common_names` does. In addition, 
clarified descriptino of what other functions do as well, whenever
it was unclear (#12)

### BUG FIXES

* Some API tokens were being blocked, fixed now (#7)
* On some operating systems (at least some versions of Windows), queries 
that included taxonomic names weren't being processed correctly. It 
is fixed now (#11)


rredlist 0.1.0
==============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (evel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 downstream dependeny. There is one error in the reverse dependency taxize, but it's a simple check of a URL returned that has changed because the base url for the API wrapped in this package has changed. A fix is ready in the reverse dependency taxize and will be submitted soon. Summary at: 
<https://github.com/ropensci/rredlist/blob/master/revdep/README.md>

---

This version makes minor improvements to documentation.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rredlist/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rredlist.git`
* Make sure to track progress upstream (i.e., on our version of `rredlist` at `ropensci/rredlist`) by doing `git remote add upstream https://github.com/ropensci/rredlist.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Run tests!
* Push up to your account
* Submit a pull request to home base at `ropensci/rredlist`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- Do not share your Redlist API key in this issue - the maintainer has their own key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

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

|package  |old   |new   |Δ  |
|:--------|:-----|:-----|:--|
|rredlist |0.6.0 |0.7.0 |*  |

# Revdeps

## New problems (1)

|package                      |version |error  |warning |note |
|:----------------------------|:-------|:------|:-------|:----|
|[taxize](problems.md#taxize) |0.9.98  |__+1__ |        |     |

# taxize

<details>

* Version: 0.9.98
* Source code: https://github.com/cran/taxize
* URL: https://docs.ropensci.org/taxize/ (website), https://github.com/ropensci/taxize (devel), https://taxize.dev (user manual)
* BugReports: https://github.com/ropensci/taxize/issues
* Date/Publication: 2020-09-18 17:40:02 UTC
* Number of recursive dependencies: 100

Run `revdep_details(,"taxize")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
      > test_check("taxize")
      taxize options
        taxon_state_messages: TRUE
      ── 1. Failure: use_iucn produces expected URL and message (@test-key_helpers.R#4
      use_iucn() not equal to "http://apiv3.iucnredlist.org/api/v3/token".
      1/1 mismatches
      x[1]: "https://apiv3.iucnredlist.org/api/v3/token"
      y[1]: "http://apiv3.iucnredlist.org/api/v3/token"
      
      ══ testthat results  ═══════════════════════════════════════════════════════════
      [ OK: 74 | SKIPPED: 216 | WARNINGS: 0 | FAILED: 1 ]
      1. Failure: use_iucn produces expected URL and message (@test-key_helpers.R#4) 
      
      Error: testthat unit tests failed
      Execution halted
    ```

# Check times

|package |version | check_time|
|:-------|:-------|----------:|
|taxize  |0.8.9   |       67.8|


*Wow, no problems at all. :)*