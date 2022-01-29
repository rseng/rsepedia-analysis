rdatacite
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rdatacite)](https://cranchecks.info/pkgs/rdatacite)
[![R-CMD-check](https://github.com/ropensci/rdatacite/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rdatacite/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/rdatacite/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rdatacite?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rdatacite)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rdatacite)](https://cran.r-project.org/package=rdatacite)
[![DOI](https://zenodo.org/badge/2521192.svg)](https://zenodo.org/badge/latestdoi/2521192)

`rdatacite` provides programmatic accesses to DataCite (https://datacite.org/) metadata

* REST API. Docs: https://support.datacite.org/docs/api and https://support.datacite.org/reference

`rdatacite` docs: https://docs.ropensci.org/rdatacite

Package API:

 - `dc_providers`
 - `dc_reports`
 - `dc_check`
 - `dc_events`
 - `dc_dois`
 - `dc_clients`
 - `dc_client_prefixes`
 - `dc_provider_prefixes`
 - `dc_status`
 - `dc_prefixes`
 - `dc_activities`

## Installation

Stable CRAN version


```r
install.packages("rdatacite")
```

Development version from github


```r
pak::pkg_install("ropensci/rdatacite")
```


```r
library('rdatacite')
```

## Result objects

Outputs from nearly all `rdatacite` functions will be of class `dc`, an S3 class that's 
simply a named list of results. You can easily remove the class via `unclass()`.
The `print.dc` method prints the data.frame for the `data`, `included`, and `reports`
slots if they exist, but hides the `meta` named list. You can get to the metadata by
indexing to it like `$meta`.

## Searching

You may want to start with `dc_dois()`.


```r
dc_dois(query = "climate change")
#> datacite: dois
#> found: 85075, pages: 400, page: 1
#> slots: data, meta, links
#> $data
#> # A tibble: 25 x 4
#>    id    type  attributes$doi $identifiers $creators $titles $publisher
#>    <chr> <chr> <chr>          <list>       <list>    <list>  <chr>     
#>  1 10.1… dois  10.15786/20.5… <list [0]>   <df[,6] … <df[,1… Mountain …
#>  2 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  3 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  4 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  5 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  6 10.2… dois  10.25676/1112… <list [0]>   <df[,6] … <df[,1… Mountain …
#>  7 10.2… dois  10.25676/1112… <list [0]>   <df[,6] … <df[,1… Mountain …
#>  8 10.2… dois  10.25675/1021… <list [0]>   <df[,6] … <df[,2… Mountain …
#>  9 10.2… dois  10.25675/1021… <list [0]>   <df[,6] … <df[,1… Mountain …
#> 10 10.2… dois  10.25675/1021… <list [0]>   <df[,6] … <df[,1… Mountain …
#> # … with 15 more rows, and 42 more variables: $container <df[,0]>,
#> #   $publicationYear <int>, $subjects <list>, $contributors <list>,
#> #   $dates <list>, $language <chr>, $types$ris <chr>, $$bibtex <chr>,
#> #   $$citeproc <chr>, $$schemaOrg <chr>, $$resourceType <chr>,
#> #   $$resourceTypeGeneral <chr>, $relatedIdentifiers <list>, $sizes <list>,
#> #   $formats <list>, $version <lgl>, $rightsList <list>, $descriptions <list>,
#> #   $geoLocations <list>, $fundingReferences <list>, $url <chr>,
#> #   $contentUrl <lgl>, $metadataVersion <int>, $schemaVersion <chr>,
#> #   $source <chr>, $isActive <lgl>, $state <chr>, $reason <lgl>,
#> #   $viewCount <int>, $downloadCount <int>, $referenceCount <int>,
#> #   $citationCount <int>, $partCount <int>, $partOfCount <int>,
#> #   $versionCount <int>, $versionOfCount <int>, $created <chr>,
#> #   $registered <chr>, $published <lgl>, $updated <chr>,
#> #   relationships$client$data$id <chr>, $$$type <chr>
#> 
#> $included
#> NULL
```

The `query` parameter supports Elasticearch query string queries. Some examples:


```r
# search within a field
dc_dois(query = "publicationYear:2016")
# fuzzy search (via *) on a nested field
dc_dois(query = "creators.familyName:mil*")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rdatacite/issues).
* License: MIT
* Get citation information for `rdatacite` in R doing `citation(package = 'rdatacite')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rdatacite 0.5.2
===============

### MINOR IMPROVEMENTS

* fix breaking test on one of the cran checks (#30)


rdatacite 0.5.0
===============

### NEW FEATURES

* Major refactor to work with the new DataCite API: all functions from the previous version are defunct; all OAI-PMH functions are gone; new functions all start with `dc_` (#24) (#29)

### MINOR IMPROVEMENTS

* all examples check if DataCite API is up before running (#28)


rdatacite 0.4.2
===============

### MINOR IMPROVEMENTS

* fix to two fixtures that had non-ascii text in them, that were causing tests to fail (#25)


rdatacite 0.4.0
===============

### MINOR IMPROVEMENTS

* pagination fixes (#18)
* fix unused httr package warning, flagged by cran team (#21)
* add .github PR and issue templates


rdatacite 0.3.0
===============

### NEW FEATURES

* Gains new functions for working with the DataCite REST API:
`dc_data_center`, `dc_data_centers`, `dc_member`, `dc_members`,
`dc_work`, `dc_works` (#13)
* Now using new version of solrium package - users shouldn't see any differences (#16)

### BUG FIXES

* Fix scientific notation (#15)
* Fix `vapply` error (#14)



rdatacite 0.1.0
===============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 3.6.3 RC
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version fixes a problem in the test suite causing a R CMD CHECK failure on R devel linux.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rdatacite/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rdatacite.git`
* Make sure to track progress upstream (i.e., on our version of `rdatacite` at `ropensci/rdatacite`) by doing `git remote add upstream https://github.com/ropensci/rdatacite.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/rdatacite`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
