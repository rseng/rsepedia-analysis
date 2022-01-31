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


*Wow, no problems at all. :)*rredlist
========

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![cran checks](https://cranchecks.info/badges/worst/rredlist)](https://cranchecks.info/pkgs/rredlist)
[![R-check](https://github.com/ropensci/rredlist/workflows/R-check/badge.svg)](https://github.com/ropensci/rredlist/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/rredlist/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rredlist?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rredlist)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rredlist)](https://cran.r-project.org/package=rredlist)

`rredlist` is an R client for the IUCN Red List (https://apiv3.iucnredlist.org/api/v3/docs). The IUCN Red List is a global list of threatened and endangered species. IUCN Red List docs: https://apiv3.iucnredlist.org/api/v3/docs

## Installation

CRAN

```{r eval=FALSE}
install.packages("rredlist")
```

Development version

```{r eval=FALSE}
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
---
title: "rredlist"
author: "Scott Chamberlain"
date: "2020-10-14"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to rredlist}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`rredlist` is an R client for the IUCN Red List (https://apiv3.iucnredlist.org/api/v3/docs). 
The IUCN Red List is a global list of threatened and endangered species.

IUCN Red List docs: http://apiv3.iucnredlist.org/api/v3/docs The web API [needs authentication](#authentication).

> What rredlist is not: [redlistr][] is a different package - not working with the IUCN Red List API; Furthermore, rredlist does not include support for the spatial API, described at
https://apiv3.iucnredlist.org/spatial.


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

## Authentication

IUCN requires you to get your own API key, an alphanumeric string that you need to send in every request. There's an helper function in the package helping you getting it at https://apiv3.iucnredlist.org/api/v3/token and storing it.


```r
rredlist::rl_use_iucn()
```

Keep this key private. You can pass the key in to each function via the `key` parameter, but it's better to store the key either as a environment variable (`IUCN_REDLIST_KEY`) or an R option (`iucn_redlist_key`) - we recommend using the former option.


## High level interface


```r
library("rredlist")
```

High level functions do the HTTP request and parse data to a data.frame for ease
of downstream use. The high level functions have no underscore on the end
of the function name, e.g., `rl_search()`


```r
rl_search('Fratercula arctica')
#> $name
#> [1] "Fratercula arctica"
#> 
#> $result
#>    taxonid    scientific_name  kingdom   phylum class           order  family
#> 1 22694927 Fratercula arctica ANIMALIA CHORDATA  AVES CHARADRIIFORMES ALCIDAE
#>        genus main_common_name        authority published_year assessment_date
#> 1 Fratercula  Atlantic Puffin (Linnaeus, 1758)           2018      2018-08-07
#>   category criteria population_trend marine_system freshwater_system
#> 1       VU  A4abcde       Decreasing          TRUE             FALSE
#>   terrestrial_system               assessor        reviewer aoo_km2  eoo_km2
#> 1               TRUE BirdLife International Westrip, J.R.S.      NA 20800000
#>   elevation_upper elevation_lower depth_upper depth_lower errata_flag
#> 1              NA              NA          NA          NA          NA
#>   errata_reason amended_flag amended_reason
#> 1            NA           NA             NA
```

> Note: there can sometimes be a discrepancy between what you get on the IUCN website and what
> you get with this package; we don't know why, the IUCN API is not an open book.

Likely a bit faster is to parse to a list only, and not take the extra data.frame parsing time


```r
rl_search('Fratercula arctica', parse = FALSE)
#> $name
#> [1] "Fratercula arctica"
#> 
#> $result
#> $result[[1]]
#> $result[[1]]$taxonid
#> [1] 22694927
#> 
#> $result[[1]]$scientific_name
#> [1] "Fratercula arctica"
...
```

For even more speed, use the low level package interface.

## Low level interface

The parsing to data.frame in the high level functions does take extra time. The low level functions
only do the HTTP request, and give back JSON without doing any more parsing. The low level functions DO have an underscore on the end
of the function name, e.g., `rl_search_()`


```r
rl_search_('Fratercula arctica')
#> [1] "{\"name\":\"Fratercula arctica\",\"result\":[{\"taxonid\":22694927,\"scientific_name\":\"Fratercula arctica\",\"kingdom\":\"ANIMALIA\",\"phylum\":\"CHORDATA\",\"class\":\"AVES\",\"order\":\"CHARADRIIFORMES\",\"family\":\"ALCIDAE\",\"genus\":\"Fratercula\",\"main_common_name\":\"Atlantic Puffin\",\"authority\":\"(Linnaeus, 1758)\",\"published_year\":2018,\"assessment_date\":\"2018-08-07\",\"category\":\"VU\",\"criteria\":\"A4abcde\",\"population_trend\":\"Decreasing\",\"marine_system\":true,\"freshwater_system\":false,\"terrestrial_system\":true,\"assessor\":\"BirdLife International\",\"reviewer\":\"Westrip, J.R.S.\",\"aoo_km2\":null,\"eoo_km2\":\"20800000\",\"elevation_upper\":null,\"elevation_lower\":null,\"depth_upper\":null,\"depth_lower\":null,\"errata_flag\":null,\"errata_reason\":null,\"amended_flag\":null,\"amended_reason\":null}]}"
```

To consume this JSON, you can use `jsonlite`


```r
library("jsonlite")
jsonlite::fromJSON(rl_search_('Fratercula arctica'))
#> $name
#> [1] "Fratercula arctica"
#> 
#> $result
#>    taxonid    scientific_name  kingdom   phylum class           order  family
#> 1 22694927 Fratercula arctica ANIMALIA CHORDATA  AVES CHARADRIIFORMES ALCIDAE
#>        genus main_common_name        authority published_year assessment_date
#> 1 Fratercula  Atlantic Puffin (Linnaeus, 1758)           2018      2018-08-07
#>   category criteria population_trend marine_system freshwater_system
#> 1       VU  A4abcde       Decreasing          TRUE             FALSE
#>   terrestrial_system               assessor        reviewer aoo_km2  eoo_km2
#> 1               TRUE BirdLife International Westrip, J.R.S.      NA 20800000
#>   elevation_upper elevation_lower depth_upper depth_lower errata_flag
#> 1              NA              NA          NA          NA          NA
#>   errata_reason amended_flag amended_reason
#> 1            NA           NA             NA
```

Or other tools, e.g., `jq` via the `jqr` R client


```r
# remotes::install_github("ropensci/jqr")
library("jqr")
rl_search_('Fratercula arctica') %>% dot()
#> $data
#> [1] "{\"name\":\"Fratercula arctica\",\"result\":[{\"taxonid\":22694927,\"scientific_name\":\"Fratercula arctica\",\"kingdom\":\"ANIMALIA\",\"phylum\":\"CHORDATA\",\"class\":\"AVES\",\"order\":\"CHARADRIIFORMES\",\"family\":\"ALCIDAE\",\"genus\":\"Fratercula\",\"main_common_name\":\"Atlantic Puffin\",\"authority\":\"(Linnaeus, 1758)\",\"published_year\":2018,\"assessment_date\":\"2018-08-07\",\"category\":\"VU\",\"criteria\":\"A4abcde\",\"population_trend\":\"Decreasing\",\"marine_system\":true,\"freshwater_system\":false,\"terrestrial_system\":true,\"assessor\":\"BirdLife International\",\"reviewer\":\"Westrip, J.R.S.\",\"aoo_km2\":null,\"eoo_km2\":\"20800000\",\"elevation_upper\":null,\"elevation_lower\":null,\"depth_upper\":null,\"depth_lower\":null,\"errata_flag\":null,\"errata_reason\":null,\"amended_flag\":null,\"amended_reason\":null}]}"
#> 
#> $args
#> $args[[1]]
#> [1] "."
#> attr(,"type")
#> [1] "dot"
#> 
#> 
#> attr(,"class")
#> [1] "jqr"
```

## Usage best practice

### Citing the IUCN Red List API

Use the function `rl_citation()`:

```r
rl_citation()
#> [1] "IUCN 2015. IUCN Red List of Threatened Species. Version 2020-2 <www.iucnredlist.org>"
```

### Rate Limiting

From the IUCN folks: "Too many frequent calls, or too many calls per day
might get your access blocked temporarily. If you're a heavy API user, the
Red List Unit asked that you contact them, as there might be better options.
They suggest a 2-second delay between your calls if you plan to make a
lot of calls."
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_history.R
\name{rl_history}
\alias{rl_history}
\title{Get historical assessments by taxon name, IUCN id, and region}
\usage{
rl_history(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get historical assessments by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_history('Loxodonta africana')
rl_history('Ursus maritimus', region = 'europe')
rl_history(id = 12392)
rl_history(id = 22823, region = 'europe')

rl_history_('Loxodonta africana')
rl_history_(id = 12392)
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_synonyms.R
\name{rl_synonyms}
\alias{rl_synonyms}
\alias{rl_synonyms_}
\title{Get species synonym information by taxonomic name}
\usage{
rl_synonyms(name = NULL, key = NULL, parse = TRUE, ...)

rl_synonyms_(name = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) Binomial taxonomic name}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get species synonym information by taxonomic name
}
\examples{
\dontrun{
rl_synonyms('Loxodonta africana')
rl_synonyms('Loxodonta africana', parse = FALSE)
rl_synonyms_('Loxodonta africana')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_narrative.R
\name{rl_narrative}
\alias{rl_narrative}
\alias{rl_narrative_}
\title{Get species narrative information by taxon name, IUCN id, and region}
\usage{
rl_narrative(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)

rl_narrative_(name = NULL, id = NULL, region = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get species narrative information by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_narrative('Fratercula arctica')
rl_narrative('Fratercula arctica', region = 'europe')
rl_narrative(id = 12392)
rl_narrative(id = 22694927, region = 'europe')

rl_narrative_('Fratercula arctica')
rl_narrative_('Fratercula arctica', region = 'europe')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_measures.R
\name{rl_measures}
\alias{rl_measures}
\alias{rl_measures_}
\title{Get species conservation measures by taxon name, IUCN id, and region}
\usage{
rl_measures(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)

rl_measures_(name = NULL, id = NULL, region = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get species conservation measures by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_measures('Fratercula arctica')
rl_measures('Fratercula arctica', region = 'europe')
rl_measures(id = 12392)
rl_measures(id = 22694927, region = 'europe')

rl_measures_('Fratercula arctica')
rl_measures_(id = 22694927, region = 'europe')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_sp_category.R
\name{rl_sp_category}
\alias{rl_sp_category}
\alias{rl_sp_category_}
\title{Get species by category}
\usage{
rl_sp_category(category, key = NULL, parse = TRUE, ...)

rl_sp_category_(category, key = NULL, parse = TRUE, ...)
}
\arguments{
\item{category}{(character) A two-letter category code. One of
"DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX", "LRlc", "LRnt", "LRcd"}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Get species by category
}
\examples{
\dontrun{
rl_sp_category('VU')
rl_sp_category('LRlc')
rl_sp_category('EN')
rl_sp_category('EX')
rl_sp_category('EX', parse = FALSE)
rl_sp_category_('EX')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_version.R
\name{rl_version}
\alias{rl_version}
\title{Get the Red List API version}
\usage{
rl_version(...)
}
\arguments{
\item{...}{Curl options passed to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
API version as character string
}
\description{
Get the Red List API version
}
\examples{
\dontrun{
rl_version()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_occ_country.R, R/rl_search.R
\name{rl_occ_country_}
\alias{rl_occ_country_}
\alias{rl_search}
\alias{rl_search_}
\title{Search by taxon name, IUCN id, and region}
\usage{
rl_occ_country_(name = NULL, id = NULL, region = NULL, key = NULL, ...)

rl_search(name = NULL, id = NULL, region = NULL, key = NULL, parse = TRUE, ...)

rl_search_(name = NULL, id = NULL, region = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Search by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_search('Fratercula arctica')
rl_search('Fratercula arctica', region = 'europe')
rl_search(id = 12392)
rl_search(id = 22694927, region = 'europe')

rl_search('Fratercula arctica', parse = FALSE)
rl_search_('Fratercula arctica')
rl_search_('Fratercula arctica', region = 'europe')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_occ_country.R
\name{rl_occ_country}
\alias{rl_occ_country}
\title{Get country occurrence by species name or ID}
\usage{
rl_occ_country(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get country occurrence by species name or ID
}
\examples{
\dontrun{
rl_occ_country('Loxodonta africana')
rl_occ_country('Fratercula arctica', region = 'europe')
rl_occ_country(id = 12392)
rl_occ_country(id = 22694927, region = 'europe')

rl_occ_country('Fratercula arctica', parse = FALSE)
rl_occ_country_('Fratercula arctica')
rl_occ_country_('Fratercula arctica', region = 'europe')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_sp_count.R
\name{rl_sp_count}
\alias{rl_sp_count}
\alias{rl_sp_count_}
\title{Get total species count of taxa in the Red List}
\usage{
rl_sp_count(key = NULL, parse = TRUE, ...)

rl_sp_count_(key = NULL, ...)
}
\arguments{
\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get total species count of taxa in the Red List
}
\examples{
\dontrun{
rl_sp_count()
rl_sp_count_()
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_growth_forms.R
\name{rl_growth_forms}
\alias{rl_growth_forms}
\alias{rl_growth_forms_}
\title{Get plant species growth forms by taxon name, IUCN id, and region}
\usage{
rl_growth_forms(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)

rl_growth_forms_(name = NULL, id = NULL, region = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get plant species growth forms by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_growth_forms('Quercus robur')
rl_growth_forms('Quercus robur', region = 'europe')
rl_growth_forms(id = 63532)
rl_growth_forms(id = 63532, region = 'europe')

rl_growth_forms('Mucuna bracteata')
rl_growth_forms('Abarema villifera')
rl_growth_forms('Adansonia perrieri')
rl_growth_forms('Adenostemma harlingii')

rl_growth_forms_('Quercus robur')
rl_growth_forms_(id = 63532, region = 'europe')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_key.R
\name{rl_use_iucn}
\alias{rl_use_iucn}
\title{Helper to get and save IUCN API key}
\usage{
rl_use_iucn()
}
\description{
Browse IUCN Red List API key request URL and
provides instruction on how to store the key.
}
\details{
Note that after filling the online form, you should
receive an API key shortly but not immediately.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_countries.R
\name{rl_countries}
\alias{rl_countries}
\alias{rl_countries_}
\title{Get countries}
\usage{
rl_countries(key = NULL, parse = TRUE, ...)

rl_countries_(key = NULL, ...)
}
\arguments{
\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get countries
}
\examples{
\dontrun{
rl_countries()
rl_countries_()
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_citation.R
\name{rl_citation}
\alias{rl_citation}
\title{Get the citation Red List API version}
\usage{
rl_citation(...)
}
\arguments{
\item{...}{Curl options passed to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
API citation as character string
}
\description{
Get the citation Red List API version
}
\examples{
\dontrun{
rl_citation()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_sp_country.R
\name{rl_sp_country}
\alias{rl_sp_country}
\alias{rl_sp_country_}
\title{Get species by country}
\usage{
rl_sp_country(country, key = NULL, parse = TRUE, ...)

rl_sp_country_(country, key = NULL, ...)
}
\arguments{
\item{country}{(character) A two-letter country code. See \code{isocode} column
in result of \code{\link[=rl_countries]{rl_countries()}} request for country codes.}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get species by country
}
\examples{
\dontrun{
rl_sp_country('AZ')
rl_sp_country('NZ')

# don't parse to data.frame, gives list
rl_sp_country('NZ', parse = FALSE)
# don't parse at all, get json back
rl_sp_country_('NZ')

# curl options
res <- rl_sp_country('NZ', verbose = TRUE)
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_comp_groups.R
\name{rl_comp_groups}
\alias{rl_comp_groups}
\alias{rl_comp_groups_}
\title{Information about comprehensive groups}
\usage{
rl_comp_groups(group = NULL, key = NULL, parse = TRUE, ...)

rl_comp_groups_(group = NULL, key = NULL, ...)
}
\arguments{
\item{group}{(character) A comprehensive group name.
Call \code{rl_comp_groups()} without passing this parameter
to get the list of comprehensive groups}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Information about comprehensive groups
}
\examples{
\dontrun{
rl_comp_groups()
rl_comp_groups('mammals')
rl_comp_groups('groupers')

rl_comp_groups_()
rl_comp_groups_('groupers')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_common_names.R
\name{rl_common_names}
\alias{rl_common_names}
\alias{rl_common_names_}
\title{Get common names for a given taxonomic name}
\usage{
rl_common_names(name = NULL, key = NULL, parse = TRUE, ...)

rl_common_names_(name = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) Binomial taxonomic name}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get common names for a given taxonomic name
}
\examples{
\dontrun{
rl_common_names('Loxodonta africana')
rl_common_names_('Loxodonta africana')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rredlist-package.R
\docType{package}
\name{rredlist-package}
\alias{rredlist-package}
\alias{rredlist}
\title{rredlist}
\description{
IUCN Red List R Client
}
\section{Taxonomic Names vs. IUCN IDs}{

From the documentation (quoting): "It is advisable wherever possible to use
the taxon name (species name) to make your API calls, rather than using IDs.
IDs are not immovable are expected to be used mainly by organisations
that work closely with the IUCN Red List."
}

\section{Authentication}{

IUCN requires you to get your own API key, an alphanumeric string that you
need to send in every request. See key A IUCN API token. See \code{\link[=rl_use_iucn]{rl_use_iucn()}}
for help getting and storing it. Get it at
https://apiv3.iucnredlist.org/api/v3/token
Keep this key private. You can pass the key in to each function via the
\code{key} parameter, but it's better to store the key either as a
environment variable (\code{IUCN_REDLIST_KEY}) or an R option
(\code{iucn_redlist_key}) - we recommend using the former option.
}

\section{High vs. Low level package APIs}{

\strong{High level API}
High level functions do the HTTP request and parse data to a data.frame for
ease of downstream use. The high level functions have no underscore on
the end of the function name, e.g., \code{\link[=rl_search]{rl_search()}}

\strong{Low level API}
The parsing to data.frame in the high level API does take extra time.
The low level API only does the HTTP request, and gives back JSON without
doing any more parsing. The low level functions DO have an underscore on
the end of the function name, e.g., \code{\link[=rl_search_]{rl_search_()}}
}

\section{No Spatial}{

This package does not include support for the spatial API, described at
https://apiv3.iucnredlist.org/spatial
}

\section{Citing the Red List API}{

Get the proper citation for the version of the Red List you are using
by programatically running \code{\link[=rl_citation]{rl_citation()}}
}

\section{Red List API Terms of Use}{

See https://www.iucnredlist.org/terms/terms-of-use
}

\section{Rate limiting}{

From the IUCN folks: Too many frequent calls, or too many calls per day
might get your access blocked temporarily. If you're a heavy API user, the
Red List Unit asked that you contact them, as there might be better options.
They suggest a 2-second delay between your calls if you plan to make a
lot of calls.
}

\section{Citing the IUCN Red List API}{

See https://apiv3.iucnredlist.org/about
}

\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_habitats.R, R/rl_history.R
\name{rl_habitats}
\alias{rl_habitats}
\alias{rl_habitats_}
\alias{rl_history_}
\title{Get species habitats by taxon name, IUCN id, and region}
\usage{
rl_habitats(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)

rl_habitats_(name = NULL, id = NULL, region = NULL, key = NULL, ...)

rl_history_(name = NULL, id = NULL, region = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get species habitats by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_habitats('Fratercula arctica')
rl_habitats('Fratercula arctica', region = 'europe')
rl_habitats(id = 12392)
rl_habitats(id = 22694927, region = 'europe')

rl_habitats_('Fratercula arctica')
rl_habitats_(id = 12392)
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_regions.R
\name{rl_regions}
\alias{rl_regions}
\alias{rl_regions_}
\title{Get regions}
\usage{
rl_regions(key = NULL, parse = TRUE, ...)

rl_regions_(key = NULL, ...)
}
\arguments{
\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get regions
}
\examples{
\dontrun{
rl_regions()
rl_regions(parse = FALSE)
rl_regions_()
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_sp_citation.R
\name{rl_sp_citation}
\alias{rl_sp_citation}
\alias{rl_sp_citation_}
\title{Get citations by taxon name, IUCN id, and region}
\usage{
rl_sp_citation(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)

rl_sp_citation_(name = NULL, id = NULL, region = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get citations by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_sp_citation('Balaena mysticetus')
rl_sp_citation('Balaena mysticetus', region = 'europe')
rl_sp_citation(id = 12392)

rl_sp_citation(id = 2467, region = 'europe')
rl_sp_citation(id = 2467, region = 'europe', parse = FALSE)
rl_sp_citation_(id = 2467, region = 'europe')
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_sp.R
\name{rl_sp}
\alias{rl_sp}
\alias{rl_sp_}
\title{Get species}
\usage{
rl_sp(page = 0, key = NULL, parse = TRUE, all = FALSE, quiet = FALSE, ...)

rl_sp_(page, key = NULL, all = FALSE, quiet = FALSE, ...)
}
\arguments{
\item{page}{(integer/numeric) Page to get. Default: 0. you can
get up to 10,000 records per page. Paging is required because
it's too much burden on a server to just "get all the data"
in one request}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{all}{(logical) to get all results or not. Default: \code{FALSE}.
this means we do the paging internally for you. result is a list
of results, so you have to bind them together yourself into
a data.frame, see example}

\item{quiet}{(logical) give progress for download or not.
Default: \code{FALSE} (that is, give progress). ignored if
\code{all = FALSE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Get species
}
\examples{
\dontrun{
rl_sp(page = 3)

# get all results
out <- rl_sp(all = TRUE)
length(out)
vapply(out, "[[", 1, "count")
all_df <- do.call(rbind, lapply(out, "[[", "result"))
head(all_df)
NROW(all_df)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rl_threats.R
\name{rl_threats}
\alias{rl_threats}
\alias{rl_threats_}
\title{Get species threats by taxon name, IUCN id, and region}
\usage{
rl_threats(
  name = NULL,
  id = NULL,
  region = NULL,
  key = NULL,
  parse = TRUE,
  ...
)

rl_threats_(name = NULL, id = NULL, region = NULL, key = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(character) An IUCN identifier}

\item{region}{(character) A region name, see \code{\link{rl_regions}} for
acceptable region identifiers (use the entries in the \code{identifier}
column)}

\item{key}{A IUCN API token. See \code{\link{rl_use_iucn}}.}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
A list, with the data in the \code{result} slot, unless using
a function with a trailing underscore, in which case json as character
string is returned.
}
\description{
Get species threats by taxon name, IUCN id, and region
}
\examples{
\dontrun{
rl_threats('Fratercula arctica')
rl_threats('Fratercula arctica', region = 'europe')
rl_threats(id = 12392)
rl_threats(id = 22694927, region = 'europe')
rl_threats(name = 'Abies numidica')
rl_threats_('Fratercula arctica')

rl_threats(id = 62290750)
}
}
\references{
API docs at https://apiv3.iucnredlist.org/api/v3/docs
}
