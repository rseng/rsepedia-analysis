handlr
======



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/handlr)](https://cranchecks.info/pkgs/handlr)
[![R-check](https://github.com/ropensci/handlr/workflows/R-check/badge.svg)](https://github.com/ropensci/handlr/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/handlr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/handlr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/handlr)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/handlr)](https://cran.r-project.org/package=handlr)


a tool for converting among citation formats.

heavily influenced by, and code ported from the Ruby gem `bolognese`

supported readers:

- citeproc
- ris
- bibtex
- codemeta
- cff

supported writers:

- citeproc
- ris
- bibtex
- schemaorg
- rdfxml
- codemeta
- cff

not supported yet, but plan to:

- crosscite

## Installation

stable version


```r
install.packages("handlr")
```

dev version


```r
remotes::install_github("ropensci/handlr")
```


```r
library("handlr")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/handlr/issues).
* License: MIT
* Get citation information for `handlr` in R doing `citation(package = 'handlr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
handlr 0.3.0
============

### DEPENDENCIES

* drop `RefManageR` from package Imports as it will likely be archived soon - add package `bibtex` to Suggests for reading/writing bibtex (can't be in Imports because it's Orphaned on CRAN) (#22)

### NEW FEATURES

* handlr gains support for Citation File Format (CFF), "plain text files with human- and machine-readable citation information for software". See https://citation-file-format.github.io/ for more info - new functions: `cff_reader()` and `cff_writer()` and associated changes in `HandlrClient`. Associated with CFF support, handlr gains new Import package `yaml`  (#16)

### MINOR IMPROVEMENTS

* improvements to Citeproc parsing: previously dropped many fields that we didn't support; now including all Citeproc fields that we don't specifically parse into extra fields prefixed with `csl_`  (#20)
* nothing changed, but see discussion of bibtex errors in case you run into them (#9)


handlr 0.2.0
============

### NEW FEATURES

* gains function `handl_to_df()`; converts any `handl` object (output from `HandlClient` or any `*_reader()` functions) to a data.frame for easier downstream data munging; `HandlClient` gains `$as_df()` method which runs `handl_to_df()`; to support this, now importing data.table package (#15) (#19) feature request by @GeraldCNelson

### MINOR IMPROVEMENTS

* now exporting the `print.handl` method. it only affects how a `handl` class object prints in the console, but is useful for making output more brief/concise (#14)
* filled out a lot more details of what a `handl` object contains. see `?handl` for the documentation (#17)


handlr 0.1.0
============

### NEW FEATURES

* Released to CRAN
## Test environments

* local OS X install, R 4.0.3 patched
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 reverse dependency. No problems were found. Summary at <https://github.com/ropensci/handlr/blob/master/revdep/README.md> 

---

This version includes dropping RefManageR package, scheduled to be archived soon, along with fixes and some new functions.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/handlr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/handlr.git`
* Make sure to track progress upstream (i.e., on our version of `handlr` at `ropensci/handlr`) by doing `git remote add upstream https://github.com/ropensci/handlr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/handlr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

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
|date     |2020-10-14                   |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|handlr  |0.2.0 |0.3.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*