rcoreoa
=======



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/rcoreoa/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rcoreoa/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/rcoreoa/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rcoreoa?branch=master)
[![cran checks](https://cranchecks.info/badges/worst/rcoreoa)](https://cranchecks.info/pkgs/rcoreoa)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcoreoa)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcoreoa)](https://cran.r-project.org/package=rcoreoa)

CORE API R client

CORE API docs: https://core.ac.uk/docs/

rcoreoa docs: https://docs.ropensci.org/rcoreoa/

Get an API key at https://core.ac.uk/api-keys/register. You'll need one, 
so do this now if you haven't yet. Once you have the key, you can pass it 
into the `key` parameter, or as a much better option store your key as an 
environment variable with the name `CORE_KEY` or an R option as `core_key`. 
See `?Startup` for how to work with env vars and R options

## About CORE

CORE's tagline is: "Aggregating the world's open access research papers"

CORE offers seamless access to millions of open access research papers, enrich
the collected data for text-mining and provide unique services to the research
community.

For more infos on CORE, see https://core.ac.uk/about

## Install


```r
install.packages("rcoreoa")
```

Development version


```r
remotes::install_github("ropensci/rcoreoa")
```


```r
library("rcoreoa")
```

## Get started

Get started with an introduction to the package: https://docs.ropensci.org/rcoreoa/articles/rcoreoa

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Aristotelis Charalampous](https://github.com/aristotelisxs)
* [Simon Goring](https://github.com/SimonGoring)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rcoreoa/issues).
* License: MIT
* Get citation information for `rcoreoa` in R doing `citation(package = 'rcoreoa')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rcoreoa 0.4.0
=============

### NEW FEATURES

* `core_articles_dedup()` function added for article deduplication (#23)
* `core_advanced_search()` overhauled because the CORE advanced search fields and syntax have completely changed. The interface changed from passing in a data.frame because we needed more flexibility to let users define how to combine multiple queries  (#20)

### MINOR IMPROVEMENTS

* using vcr for test suite http request caching (#21)
* add max value of `limit` (number of records per request) parameter to documentation (#18)

rcoreoa 0.3.0
=============

### NEW FEATURES

* gains new methods `core_articles_search()` and `core_articles_search_()` for searching for articles (#15)
* gains new manual file `?core_cache` for an overview of caching
* package `hoardr` now used for managing file caching (#13)
* function `core_advanced_search()` gains support for the `language` filter (#16) thanks @chreman

### MINOR IMPROVEMENTS

* support for passing in more than one journal/article/etc. identifier for `core_articles_pdf()`/`core_articles_history()` - already supported in other functions (#10)
* `overwrite` parameter was being ignored in `core_articles_pdf()`, now is passed on internally (#12)
* `parse` parameter dropped in `core_articles_pdf()`, only used internally (#11)
* Improved docs on how to get and use API keys (#14)
* Improve documentation about what a 404 error response means - that the thing reqeusted does not exist (#9)


rcoreoa 0.1.0
=============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.2 Patched
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version adds a function, changes a function interface due to changes in the remove service, and improves documentation.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rcoreoa/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rcoreoa.git`
* Make sure to track progress upstream (i.e., on our version of `rcoreoa` at `ropensci/rcoreoa`) by doing `git remote add upstream https://github.com/ropensci/rcoreoa.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rcoreoa`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) 

Make sure not to share any secrets. -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
