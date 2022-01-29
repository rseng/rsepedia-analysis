ghql
====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/ghql)](https://cranchecks.info/pkgs/ghql)
[![R-check](https://github.com/ropensci/ghql/workflows/R-check/badge.svg)](https://github.com/ropensci/ghql/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/ghql/coverage.svg?branch=master)](https://codecov.io/github/ropensci/ghql?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ghql)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/ghql)](https://cran.r-project.org/package=ghql)

`ghql` - a GraphQL client for R

GraphQL - <https://graphql.org>

Examples of GraphQL APIs:

* GitHub: <https://docs.github.com/en/graphql/guides/introduction-to-graphql>
* Opentargets: <https://genetics-docs.opentargets.org/technical-pipeline/graphql-api>
* Countries GraphQL API: <https://github.com/trevorblades/countries>

Other GraphQL R packages:

* [graphql][] - GraphQL query parser
* [gqlr][] - GraphQL server and query methods

## Install

CRAN version


```r
install.packages("ghql")
```

Development version


```r
remotes::install_github("ropensci/ghql")
```


```r
library("ghql")
library("jsonlite")
library("dplyr")
```

## Package Documentation

<https://docs.ropensci.org/ghql/>

## Meta

* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[gqlr]: https://github.com/schloerke/gqlr
[graphql]: https://github.com/ropensci/graphql
[libgraphqlparser]: https://github.com/graphql/libgraphqlparser
ghql 0.1.0
==========

### NEW FEATURES

* Released to CRAN
## Test environments

* local OS X install, R 3.6.2 patched
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

New submission

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.

---

This is a new release. I have read and agree to the the 
CRAN policies at https://cran.r-project.org/web/packages/policies.html

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

* Submit an issue on the [Issues page](https://github.com/ropensci/ghql/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/ghql.git`
* Make sure to track progress upstream (i.e., on our version of `ghql` at `ropensci/ghql`) by doing `git remote add upstream https://github.com/ropensci/ghql.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/ghql`

### Check out our [discussion forum](https://discuss.ropensci.org)
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
