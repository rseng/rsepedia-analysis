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


*Wow, no problems at all. :)*ritis
=====

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("ritis")
```

Dev version

```{r eval=FALSE}
remotes::install_github("ropensci/ritis")
```

```{r}
library("ritis")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/ritis/issues).
* License: MIT
* Get citation information for `ritis` in R doing `citation(package = 'ritis')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "ritis introduction"
author: "Scott Chamberlain"
date: "2021-02-01"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{ritis introduction}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



An interface to the Integrated Taxonomic Information System (ITIS)

See also the taxize book (https://taxize.dev/) for 
a manual on working with taxonomic data in R, including with ITIS data.

How to cite ITIS. From https://itis.gov/citation.html 

To cite data obtained from ITIS, the following citation format is offered as a suggestion:

    Retrieved [month, day, year], from the Integrated Taxonomic Information System on-line database, http://www.itis.gov.


ITIS is one of many different taxonomic data sources. Other include: Catalogue of Life (and COL+), NCBI taxonomy, International Plant Names Index, Index Fungorum, and more. The Wikipedia entry (https://en.wikipedia.org/wiki/Integrated_Taxonomic_Information_System) states that ITIS has a North American focus, but includes many taxa not in North America.

## Terminology

* "mononomial": a taxonomic name with one part, e.g, _Poa_
* "binomial": a taxonomic name with two parts, e.g, _Poa annua_
* "trinomial": a taxonomic name with three parts, e.g, _Poa annua annua_

## Installation

Install from CRAN


```r
install.packages("ritis")
```

Or install the development version from GitHub


```r
remotes::install_github("ropensci/ritis")
```

Load `ritis`


```r
library("ritis")
```

## ITIS Solr interface

There are four methods.

* `itis_search()` - Search
* `itis_group()` - Group
* `itis_highlight()` - Hightlight
* `itis_facet()` - Facet

These four methods use the equivalent functions in the package `solrium`, e.g.,
`ritis::itis_search()` uses `solrium::solr_search()`, etc. The `itis_*()` functions
simply use `...` to allow users to pass on parameters to the wrapped `solrium`
functions. So do read the `solrium` docs.

ITIS Solr API docs: https://www.itis.gov/solr_documentation.html

Some examples:

matches only monomials


```r
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{0,0}*/")
#> # A tibble: 10 x 25
#>    tsn   nameWInd nameWOInd unit1 usage credibilityRati… completenessRat…
#>    <chr> <chr>    <chr>     <chr> <chr> <chr>            <chr>           
#>  1 1526… Ornitho… Ornithoi… Orni… valid No review; untr… unknown         
#>  2 1526… Ornitho… Ornithoc… Orni… valid No review; untr… unknown         
#>  3 1526… Ornitho… Ornithom… Orni… valid No review; untr… unknown         
#>  4 1526… Pseudol… Pseudoly… Pseu… valid No review; untr… unknown         
#>  5 1527… Stilbom… Stilbome… Stil… valid No review; untr… unknown         
#>  6 1527… Nycteri… Nycterib… Nyct… valid No review; untr… unknown         
#>  7 1527… Nycteri… Nycterib… Nyct… inva… No review; untr… unknown         
#>  8 1527… Nycteri… Nycterib… Nyct… valid No review; untr… unknown         
#>  9 1527… Basilia  Basilia   Basi… valid No review; untr… unknown         
#> 10 1527… Strebli… Streblin… Stre… valid No review; untr… unknown         
#> # … with 18 more variables: currencyRating <chr>, kingdom <chr>,
#> #   parentTSN <chr>, rankID <chr>, rank <chr>, hierarchySoFar <chr>,
#> #   hierarchySoFarWRanks <chr>, hierarchyTSN <chr>, otherSource <chr>,
#> #   createDate <chr>, updateDate <chr>, hierarchicalSort <chr>,
#> #   `_version_` <dbl>, synonyms <chr>, synonymTSNs <chr>, vernacular <chr>,
#> #   unacceptReason <chr>, acceptedTSN <chr>
```

matches only binomials


```r
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{1,1}[A-Za-z0-9]*/")
#> # A tibble: 10 x 24
#>    tsn   nameWInd nameWOInd unit1 unit2 usage unacceptReason credibilityRati…
#>    <chr> <chr>    <chr>     <chr> <chr> <chr> <chr>          <chr>           
#>  1 1526… Feronia… Feronia … Fero… spin… inva… junior synonym No review; untr…
#>  2 1526… Ornitho… Ornithoi… Orni… conf… valid <NA>           No review; untr…
#>  3 1526… Ornitho… Ornithom… Orni… conf… inva… junior synonym No review; untr…
#>  4 1526… Ornitho… Ornithoi… Orni… vici… valid <NA>           No review; untr…
#>  5 1526… Ornitho… Ornithoi… Orni… prom… inva… junior synonym No review; untr…
#>  6 1526… Ornitho… Ornithom… Orni… vici… inva… junior synonym No review; untr…
#>  7 1526… Ornitho… Ornithoc… Orni… eryt… valid <NA>           No review; untr…
#>  8 1526… Ornitho… Ornithom… Orni… bute… inva… junior synonym No review; untr…
#>  9 1526… Ornitho… Ornithom… Orni… eryt… inva… junior synonym No review; untr…
#> 10 1526… Ornitho… Ornithom… Orni… nebu… inva… junior synonym No review; untr…
#> # … with 16 more variables: taxonAuthor <chr>, kingdom <chr>, rankID <chr>,
#> #   rank <chr>, hierarchySoFar <chr>, hierarchySoFarWRanks <chr>,
#> #   hierarchyTSN <chr>, synonyms <chr>, synonymTSNs <chr>, otherSource <chr>,
#> #   acceptedTSN <chr>, createDate <chr>, updateDate <chr>, `_version_` <dbl>,
#> #   parentTSN <chr>, hierarchicalSort <chr>
```

The syntax for `itis_search()` can be a bit hard to grasp. See this ITIS page https://itis.gov/solr_examples.html for help on generating the syntax they want for specific searches.

## ITIS REST API interface

ITIS REST API docs: http://www.itis.gov/ws_description.html

The following are some example uses. There are many more methods not shown below

-------

Get accepted names for a TSN


```r
accepted_names(tsn = 504239)
#> # A tibble: 1 x 3
#>   acceptedName        acceptedTsn author    
#>   <chr>               <chr>       <chr>     
#> 1 Dasiphora fruticosa 836659      (L.) Rydb.
```

Get common names for a TSN


```r
common_names(tsn = 183833)
#> # A tibble: 3 x 3
#>   commonName          language tsn   
#>   <chr>               <chr>    <chr> 
#> 1 African hunting dog English  183833
#> 2 African Wild Dog    English  183833
#> 3 Painted Hunting Dog English  183833
```

Full hierarchy for a TSN


```r
hierarchy_full(tsn = 37906)
#> # A tibble: 60 x 5
#>    parentname        parenttsn rankname      taxonname       tsn   
#>    <chr>             <chr>     <chr>         <chr>           <chr> 
#>  1 ""                ""        Kingdom       Plantae         202422
#>  2 "Plantae"         "202422"  Subkingdom    Viridiplantae   954898
#>  3 "Viridiplantae"   "954898"  Infrakingdom  Streptophyta    846494
#>  4 "Streptophyta"    "846494"  Superdivision Embryophyta     954900
#>  5 "Embryophyta"     "954900"  Division      Tracheophyta    846496
#>  6 "Tracheophyta"    "846496"  Subdivision   Spermatophytina 846504
#>  7 "Spermatophytina" "846504"  Class         Magnoliopsida   18063 
#>  8 "Magnoliopsida"   "18063"   Superorder    Asteranae       846535
#>  9 "Asteranae"       "846535"  Order         Asterales       35419 
#> 10 "Asterales"       "35419"   Family        Asteraceae      35420 
#> # … with 50 more rows
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_any_match_paged.R
\name{search_any_match_paged}
\alias{search_any_match_paged}
\title{Search for any matched page}
\usage{
search_any_match_paged(
  x,
  pagesize = NULL,
  pagenum = NULL,
  ascend = NULL,
  wt = "json",
  raw = FALSE,
  ...
)
}
\arguments{
\item{x}{text or taxonomic serial number (TSN) (character or numeric)}

\item{pagesize}{An integer containing the page size (numeric)}

\item{pagenum}{An integer containing the page number (numeric)}

\item{ascend}{A boolean containing true for ascending sort order or false
for descending (logical)}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame

a data.frame
}
\description{
Search for any matched page
}
\examples{
\dontrun{
search_any_match_paged(x=202385, pagesize=100, pagenum=1, ascend=FALSE)
search_any_match_paged(x="Zy", pagesize=100, pagenum=1, ascend=FALSE)
}
}
\seealso{
\code{\link{search_anymatch}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itis_highlight.R
\name{itis_highlight}
\alias{itis_highlight}
\title{ITIS Solr highlight}
\usage{
itis_highlight(..., proxy = NULL, callopts = list())
}
\arguments{
\item{...}{Arguments passed on to the \code{params} parameter of
the \code{\link[solrium:solr_highlight]{solrium::solr_highlight()}} function. See \link{solr_fields} for possible
parameters, and examples below}

\item{proxy}{List of arguments for a proxy connection,
including one or more of: url, port, username, password,
and auth. See \code{\link[crul:proxies]{crul::proxy()}} for  help, which is used to
construct the proxy connection.}

\item{callopts}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
ITIS Solr highlight
}
\examples{
\dontrun{
itis_highlight(q = "rank:Species", hl.fl = 'rank', rows=10)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accepted_names.R
\name{accepted_names}
\alias{accepted_names}
\title{Get accepted names from tsn}
\usage{
accepted_names(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
Zero row data.frame if the name is accepted, otherwise a data.frame
with information on the currently accepted name
}
\description{
Get accepted names from tsn
}
\examples{
\dontrun{
# TSN accepted - good name, empty data.frame returned
accepted_names(tsn = 208527)

# TSN not accepted - input TSN is old name, non-empty data.frame returned
accepted_names(tsn = 504239)

# raw json
accepted_names(tsn = 208527, raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsn_by_vernacular_language.R
\name{tsn_by_vernacular_language}
\alias{tsn_by_vernacular_language}
\title{Get tsn by vernacular language}
\usage{
tsn_by_vernacular_language(language, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{language}{A string containing the language. This is a language string,
not the international language code (character)}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Get tsn by vernacular language
}
\examples{
\dontrun{
tsn_by_vernacular_language(language = "french")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_sources.R
\name{other_sources}
\alias{other_sources}
\title{Returns a list of the other sources used for the TSN.}
\usage{
other_sources(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Returns a list of the other sources used for the TSN.
}
\examples{
\dontrun{
# results
other_sources(tsn=182662)
# no results
other_sources(tsn=2085272) 
# get xml
other_sources(tsn=182662, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_authorship.R
\name{taxon_authorship}
\alias{taxon_authorship}
\title{Returns the author information for the TSN.}
\usage{
taxon_authorship(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Returns the author information for the TSN.
}
\examples{
\dontrun{
taxon_authorship(tsn = 183671)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/currency.R
\name{currency}
\alias{currency}
\title{Get currency from tsn}
\usage{
currency(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Get currency from tsn
}
\examples{
\dontrun{
# currency data
currency(tsn=28727)
currency(tsn=28727, wt = "xml")
# no currency dat
currency(526852)
currency(526852, raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global_species_completeness.R
\name{global_species_completeness}
\alias{global_species_completeness}
\title{Get global species completeness from tsn}
\usage{
global_species_completeness(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get global species completeness from tsn
}
\examples{
\dontrun{
global_species_completeness(tsn = 180541)
global_species_completeness(180541, wt = "xml")
global_species_completeness(180541, wt = "json", raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terms.R
\name{terms}
\alias{terms}
\title{Get ITIS terms, i.e., tsn's, authors, common names, and scientific names}
\usage{
terms(query, what = "both", wt = "json", raw = FALSE, ...)
}
\arguments{
\item{query}{One or more common or scientific names, or partial names}

\item{what}{One of both (search common and scientific names), common
(search just common names), or scientific (search just scientific names)}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get ITIS terms, i.e., tsn's, authors, common names, and scientific names
}
\examples{
\dontrun{
# Get terms searching both common and scientific names
terms(query='bear')

# Get terms searching just common names
terms(query='tarweed', "common")

# Get terms searching just scientific names
terms(query='Poa annua', "scientific")

# many at once
terms(query=c('Poa annua', 'Pinus contorta'), "scientific")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ritis-package.R
\docType{data}
\name{solr_fields}
\alias{solr_fields}
\title{List of fields that can be used in \link{solr} functions}
\format{
A list of length 36
}
\source{
\url{https://www.itis.gov/solr_documentation.html}
}
\description{
Each element in the list has a list of length tree, with:
}
\details{
\itemize{
\item field: the field name, this is the name you can use in your
queries
\item definition: the definition of the field
\item example: an example value
}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank_names.R
\name{rank_names}
\alias{rank_names}
\title{Provides a list of all the unique rank names contained in the database and
their kingdom and rank ID values.}
\usage{
rank_names(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame, with columns:
\itemize{
\item kingdomname
\item rankid
\item rankname
}
}
\description{
Provides a list of all the unique rank names contained in the database and
their kingdom and rank ID values.
}
\examples{
\dontrun{
rank_names()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scientific_name.R
\name{scientific_name}
\alias{scientific_name}
\title{Returns the scientific name for the TSN. Also returns the component parts
(names and indicators) of the scientific name.}
\usage{
scientific_name(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Returns the scientific name for the TSN. Also returns the component parts
(names and indicators) of the scientific name.
}
\examples{
\dontrun{
scientific_name(tsn = 531894)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geographic_divisions.R
\name{geographic_divisions}
\alias{geographic_divisions}
\title{Get geographic divisions from tsn}
\usage{
geographic_divisions(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get geographic divisions from tsn
}
\examples{
\dontrun{
geographic_divisions(tsn = 180543)

geographic_divisions(tsn = 180543, wt = "xml")

geographic_divisions(tsn = 180543, wt = "json", raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsn2lsid.R
\name{tsn2lsid}
\alias{tsn2lsid}
\title{Gets the unique LSID for the TSN, or an empty result if there is no match.}
\usage{
tsn2lsid(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a character string, an LSID, or \code{NULL} if nothing found
}
\description{
Gets the unique LSID for the TSN, or an empty result if there is no match.
}
\examples{
\dontrun{
tsn2lsid(tsn = 155166)
tsn2lsid(tsn = 333333333)
tsn2lsid(155166, raw = TRUE)
tsn2lsid(155166, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itis_group.R
\name{itis_group}
\alias{itis_group}
\title{ITIS Solr group search}
\usage{
itis_group(..., proxy = NULL, callopts = list())
}
\arguments{
\item{...}{Arguments passed on to the \code{params} parameter of
the \code{\link[solrium:solr_group]{solrium::solr_group()}} function. See \link{solr_fields} for possible
parameters, and examples below}

\item{proxy}{List of arguments for a proxy connection,
including one or more of: url, port, username, password,
and auth. See \code{\link[crul:proxies]{crul::proxy()}} for  help, which is used to
construct the proxy connection.}

\item{callopts}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
ITIS Solr group search
}
\examples{
\dontrun{
x <- itis_group(q = "nameWOInd:/[A-Za-z0-9]*[\%20]{0,0}*/",
   group.field = 'rank', group.limit = 3)
head(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/publications.R
\name{publications}
\alias{publications}
\title{Returns a list of the pulications used for the TSN.}
\usage{
publications(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Returns a list of the pulications used for the TSN.
}
\examples{
\dontrun{
publications(tsn = 70340)
publications(tsn = 70340, wt = "xml")

publications(tsn = 70340, verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_record.R
\name{full_record}
\alias{full_record}
\title{Get full record from TSN or lsid}
\usage{
full_record(tsn = NULL, lsid = NULL, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{lsid}{lsid for a taxonomic group (character)}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get full record from TSN or lsid
}
\examples{
\dontrun{
# from tsn
full_record(tsn = 50423)
full_record(tsn = 202385)
full_record(tsn = 183833)

full_record(tsn = 183833, wt = "xml")
full_record(tsn = 183833, raw = TRUE)

# from lsid
full_record(lsid = "urn:lsid:itis.gov:itis_tsn:180543")
full_record(lsid = "urn:lsid:itis.gov:itis_tsn:180543")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experts.R
\name{experts}
\alias{experts}
\title{Get expert information for the TSN.}
\usage{
experts(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get expert information for the TSN.
}
\examples{
\dontrun{
experts(tsn = 180544)
experts(180544, wt = "xml")
experts(180544, raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage.R
\name{coverage}
\alias{coverage}
\title{Get coverge from tsn}
\usage{
coverage(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get coverge from tsn
}
\examples{
\dontrun{
# coverage data
coverage(tsn=28727)
# no coverage data
coverage(526852)
coverage(526852, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ritis-package.R
\docType{package}
\name{ritis-package}
\alias{ritis-package}
\alias{ritis}
\title{ritis}
\description{
Interface to Integrated Taxonomic Information (ITIS)
}
\section{ritis package API}{

All functions that start with \code{itis_} work with the ITIS Solr
API described at \url{https://www.itis.gov/solr_documentation.html},
which uses the package \pkg{solrium}, and these functions have you
use the \pkg{solrium} function interfaces, so you can pass on parameters
to the \pkg{solrium} functions - so the \pkg{solrium} docs are important
here.

All other functions work with the ITIS REST API described at
\url{https://www.itis.gov/ws_description.html}. For these methods,
they can grab data in either JSON or XML format. JSON is the default.
We parse the JSON to R native format, either data.frame, character
string, or list. You can get raw JSON as a character string back,
or raw XML as a character string, and then parse yourself with
\pkg{jsonlite} or \pkg{xml2}

You'll also be interested in the taxize book
\url{https://taxize.dev/}
}

\section{Terminology}{

\itemize{
\item "mononomial": a taxonomic name with one part, e.g, \emph{Poa}
\item "binomial": a taxonomic name with two parts, e.g, \emph{Poa annua}
\item "trinomial": a taxonomic name with three parts, e.g, \emph{Poa annua annua}
}
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/usage.R
\name{usage}
\alias{usage}
\title{Returns the usage information for the TSN.}
\usage{
usage(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Returns the usage information for the TSN.
}
\examples{
\dontrun{
usage(tsn = 526852)
usage(tsn = 526852, raw = TRUE)
usage(tsn = 526852, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_metadata.R
\name{core_metadata}
\alias{core_metadata}
\title{Get core metadata from tsn}
\usage{
core_metadata(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get core metadata from tsn
}
\examples{
\dontrun{
# coverage and currrency data
core_metadata(tsn=28727)
core_metadata(tsn=28727, wt = "xml")
# no coverage or currrency data
core_metadata(183671)
core_metadata(183671, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credibility.R
\name{credibility}
\alias{credibility}
\alias{credibility_rating}
\alias{credibility_ratings}
\title{Get credibility rating from tsn}
\usage{
credibility_rating(tsn, wt = "json", raw = FALSE, ...)

credibility_ratings(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Get credibility rating from tsn
}
\details{
methods:
\itemize{
\item credibility_rating: Get credibility rating for a tsn
\item credibility_ratings: Get possible credibility ratings
}
}
\examples{
\dontrun{
credibility_rating(tsn = 526852)
credibility_rating(526852, wt = "xml")
credibility_rating(526852, raw = TRUE)

credibility_ratings()
credibility_ratings(wt = "xml")
credibility_ratings(raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itis_search.R
\name{itis_search}
\alias{itis_search}
\title{ITIS Solr search}
\usage{
itis_search(..., proxy = NULL, callopts = list())
}
\arguments{
\item{...}{Arguments passed on to the \code{params} parameter of
the \code{\link[solrium:solr_search]{solrium::solr_search()}} function. See \link{solr_fields} for possible
parameters, and examples below}

\item{proxy}{List of arguments for a proxy connection,
including one or more of: url, port, username, password,
and auth. See \code{\link[crul:proxies]{crul::proxy()}} for  help, which is used to
construct the proxy connection.}

\item{callopts}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
ITIS Solr search
}
\details{
The syntax for this function can be a bit hard to grasp. See
https://itis.gov/solr_examples.html for help on generating the
syntax ITIS wants for specific searches.
}
\examples{
\dontrun{
itis_search(q = "tsn:182662")

# get all orders within class Aves (birds)
z <- itis_search(q = "rank:Class AND nameWOInd:Aves")
hierarchy_down(z$tsn)

# get taxa "downstream" from a target taxon
## taxize and taxizedb packages have downstream() fxns, but
## you can do a similar thing here by iteratively drilling down
## the taxonomic hierarchy
## here, we get families within Aves
library(data.table)
aves <- itis_search(q = "rank:Class AND nameWOInd:Aves")
aves_orders <- hierarchy_down(aves$tsn)
aves_families <- lapply(aves_orders$tsn, hierarchy_down)
rbindlist(aves_families)

# the tila operator
itis_search(q = "nameWOInd:Liquidamber\\\\ styraciflua~0.4")

# matches only monomials
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{0,0}*/")

# matches only binomials
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{1,1}[A-Za-z0-9]*/")

# matches only trinomials
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{1,1}[A-Za-z0-9]*[ ]{1,1}[A-Za-z0-9]*/")

# matches binomials or trinomials
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{1,1}[A-Za-z0-9]*[ ]{0,1}[A-Za-z0-9]*/")

itis_search(q = "nameWOInd:Poa\\\\ annua")

# pagination
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{0,0}*/", rows = 2)
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{0,0}*/", rows = 200)

# select fields to return
itis_search(q = "nameWOInd:/[A-Za-z0-9]*[ ]{0,0}*/",
   fl = c('nameWInd', 'tsn'))
}
}
\references{
\url{https://www.itis.gov/solr_documentation.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parent_tsn.R
\name{parent_tsn}
\alias{parent_tsn}
\title{Returns the parent TSN for the entered TSN.}
\usage{
parent_tsn(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Returns the parent TSN for the entered TSN.
}
\examples{
\dontrun{
parent_tsn(tsn = 202385)
parent_tsn(tsn = 202385, raw = TRUE)
parent_tsn(tsn = 202385, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_common.R
\name{search_common}
\alias{search_common}
\title{Search for tsn by common name}
\usage{
search_common(x, from = "all", wt = "json", raw = FALSE, ...)
}
\arguments{
\item{x}{text or taxonomic serial number (TSN) (character or numeric)}

\item{from}{(character) One of "all", "begin", or "end". See Details.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Search for tsn by common name
}
\details{
The \code{from} parameter:
\itemize{
\item all - Search against the \code{searchByCommonName} API route, which
searches entire name string
\item begin - Search against the \code{searchByCommonNameBeginsWith} API
route, which searches for a match at the beginning of a name string
\item end - Search against the \code{searchByCommonNameEndsWith} API route,
which searches for a match at the end of a name string
}
}
\examples{
\dontrun{
search_common("american bullfrog")
search_common("ferret-badger")
search_common("polar bear")

# comparison: all, begin, end
search_common("inch")
search_common("inch", from = "begin")
search_common("inch", from = "end")

# end
search_common("snake", from = "end")
}
}
\seealso{
\code{\link[=search_scientific]{search_scientific()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/last_change_date.R
\name{last_change_date}
\alias{last_change_date}
\title{Provides the date the ITIS database was last updated}
\usage{
last_change_date(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
character value with a date
}
\description{
Provides the date the ITIS database was last updated
}
\examples{
\dontrun{
last_change_date()
last_change_date(wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsid2tsn.R
\name{lsid2tsn}
\alias{lsid2tsn}
\title{Gets the TSN corresponding to the LSID, or an empty result if there is no match.}
\usage{
lsid2tsn(lsid, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{lsid}{(character) lsid for a taxonomic group. Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Gets the TSN corresponding to the LSID, or an empty result if there is no match.
}
\examples{
\dontrun{
lsid2tsn(lsid="urn:lsid:itis.gov:itis_tsn:28726")
lsid2tsn(lsid="urn:lsid:itis.gov:itis_tsn:28726", wt = "xml")
lsid2tsn("urn:lsid:itis.gov:itis_tsn:0")
lsid2tsn("urn:lsid:itis.gov:itis_tsn:0", wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common_names.R
\name{common_names}
\alias{common_names}
\title{Get common names from tsn}
\usage{
common_names(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Get common names from tsn
}
\examples{
\dontrun{
common_names(tsn=183833)
common_names(tsn=183833, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solr.R
\name{solr}
\alias{solr}
\title{ITIS Solr Methods}
\description{
ITIS provides access to their data via their Solr service described at
\url{https://www.itis.gov/solr_documentation.html}.  This is a powerful
interace to ITIS data as you have access to a very flexible query
interface.
}
\details{
See \link{solr_fields} and \url{https://www.itis.gov/solr_documentation.html}
for guidance on available fields.
}
\section{Functions}{

\itemize{
\item \code{\link[=itis_search]{itis_search()}} - Search
\item \code{\link[=itis_group]{itis_group()}} - Group
\item \code{\link[=itis_highlight]{itis_highlight()}} - Highlight
\item \code{\link[=itis_facet]{itis_facet()}} - Facet
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank_name.R
\name{rank_name}
\alias{rank_name}
\title{Returns the kingdom and rank information for the TSN.}
\usage{
rank_name(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame, with rank name and other info
}
\description{
Returns the kingdom and rank information for the TSN.
}
\examples{
\dontrun{
rank_name(tsn = 202385)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synonym_names.R
\name{synonym_names}
\alias{synonym_names}
\title{Returns a list of the synonyms (if any) for the TSN.}
\usage{
synonym_names(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Returns a list of the synonyms (if any) for the TSN.
}
\examples{
\dontrun{
synonym_names(tsn=183671) # tsn not accepted
synonym_names(tsn=526852) # tsn accepted
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comment_detail.R
\name{comment_detail}
\alias{comment_detail}
\title{Get comment detail from TSN}
\usage{
comment_detail(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A data.frame with results.
}
\description{
Get comment detail from TSN
}
\examples{
\dontrun{
comment_detail(tsn=180543)
comment_detail(tsn=180543, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geographic_values.R
\name{geographic_values}
\alias{geographic_values}
\title{Get all possible geographic values}
\usage{
geographic_values(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
character vector of geographic names
}
\description{
Get all possible geographic values
}
\examples{
\dontrun{
geographic_values()
geographic_values(wt = "xml")
geographic_values(wt = "json", raw = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/record.R
\name{record}
\alias{record}
\title{Gets a record from an LSID}
\usage{
record(lsid, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{lsid}{lsid for a taxonomic group (character). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Gets a record from an LSID
}
\details{
Gets the partial ITIS record for the TSN in the LSID, found by comparing the
TSN in the search key to the TSN field. Returns an empty result set if
there is no match or the TSN is invalid.
}
\examples{
\dontrun{
record(lsid = "urn:lsid:itis.gov:itis_tsn:180543")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itis_facet.R
\name{itis_facet}
\alias{itis_facet}
\title{ITIS Solr facet}
\usage{
itis_facet(..., proxy = NULL, callopts = list())
}
\arguments{
\item{...}{Arguments passed on to the \code{params} parameter of
the \code{\link[solrium:solr_facet]{solrium::solr_facet()}} function. See \link{solr_fields} for possible
parameters, and examples below}

\item{proxy}{List of arguments for a proxy connection,
including one or more of: url, port, username, password,
and auth. See \code{\link[crul:proxies]{crul::proxy()}} for  help, which is used to
construct the proxy connection.}

\item{callopts}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
ITIS Solr facet
}
\examples{
\dontrun{
itis_facet(q = "rank:Species", rows = 0, facet.field = "kingdom")$facet_fields

x <- itis_facet(q = "hierarchySoFar:*$Aves$* AND rank:Species AND usage:valid",
   facet.pivot = "nameWInd,vernacular", facet.limit = -1, facet.mincount = 1,
   rows = 0)
head(x$facet_pivot$`nameWInd,vernacular`)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/review_year.R
\name{review_year}
\alias{review_year}
\title{Returns the review year for the TSN.}
\usage{
review_year(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Returns the review year for the TSN.
}
\examples{
\dontrun{
review_year(tsn = 180541)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hierarchy.R
\name{hierarchy}
\alias{hierarchy}
\alias{hierarchy_down}
\alias{hierarchy_up}
\alias{hierarchy_full}
\title{Get hierarchy down from tsn}
\usage{
hierarchy_down(tsn, wt = "json", raw = FALSE, ...)

hierarchy_up(tsn, wt = "json", raw = FALSE, ...)

hierarchy_full(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get hierarchy down from tsn
}
\details{
Hierarchy methods:
\itemize{
\item hierarchy_down: Get hierarchy down from tsn
\item hierarchy_up: Get hierarchy up from tsn
\item hierarchy_full: Get full hierarchy from tsn
}
}
\examples{
\dontrun{
## Full down (class Mammalia)
hierarchy_down(tsn=179913)

## Full up (genus Agoseris)
hierarchy_up(tsn=36485)

## Full hierarchy
### genus Liatris
hierarchy_full(tsn=37906)
### get raw data back
hierarchy_full(tsn=37906, raw = TRUE)
### genus Baetis, get xml back
hierarchy_full(100800, wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unacceptability_reason.R
\name{unacceptability_reason}
\alias{unacceptability_reason}
\title{Returns the unacceptability reason, if any, for the TSN.}
\usage{
unacceptability_reason(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Returns the unacceptability reason, if any, for the TSN.
}
\examples{
\dontrun{
unacceptability_reason(tsn = 183671)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/any_match_count.R
\name{any_match_count}
\alias{any_match_count}
\title{Get any match count.}
\usage{
any_match_count(x, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{x}{text or taxonomic serial number (TSN) (character or numeric)}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
An integer containing the number of matches the search will return.
}
\description{
Get any match count.
}
\examples{
\dontrun{
any_match_count(x = 202385)
any_match_count(x = "dolphin")
any_match_count(x = "dolphin", wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/description.R
\name{description}
\alias{description}
\title{Get description of the ITIS service}
\usage{
description(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a string, the ITIS web service description
}
\description{
Get description of the ITIS service
}
\examples{
\dontrun{
description()
description(wt = "xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_scientific.R
\name{search_scientific}
\alias{search_scientific}
\title{Search by scientific name}
\usage{
search_scientific(x, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{x}{text or taxonomic serial number (TSN) (character or numeric)}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Search by scientific name
}
\examples{
\dontrun{
search_scientific("Tardigrada")
search_scientific("Quercus douglasii")
}
}
\seealso{
\code{\link{search_common}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jurisdiction.R
\name{jurisdiction}
\alias{jurisdiction}
\alias{jurisdictional_origin}
\alias{jurisdiction_origin_values}
\alias{jurisdiction_values}
\title{Get jurisdictional origin from tsn}
\usage{
jurisdictional_origin(tsn, wt = "json", raw = FALSE, ...)

jurisdiction_origin_values(wt = "json", raw = FALSE, ...)

jurisdiction_values(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
\itemize{
\item jurisdictional_origin: data.frame
\item jurisdiction_origin_values: data.frame
\item jurisdiction_values: character vector
}
}
\description{
Get jurisdictional origin from tsn
}
\details{
Jurisdiction methods:
\itemize{
\item jurisdictional_origin: Get jurisdictional origin from tsn
\item jurisdiction_origin_values: Get jurisdiction origin values
\item jurisdiction_values: Get all possible jurisdiction values
}
}
\examples{
\dontrun{
jurisdictional_origin(tsn=180543)
jurisdictional_origin(tsn=180543, wt = "xml")

jurisdiction_origin_values()

jurisdiction_values()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vernacular_languages.R
\name{vernacular_languages}
\alias{vernacular_languages}
\title{Provides a list of the unique languages used in the vernacular table.}
\usage{
vernacular_languages(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a character vector of verncular names
}
\description{
Provides a list of the unique languages used in the vernacular table.
}
\examples{
\dontrun{
vernacular_languages()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_anymatch.R
\name{search_anymatch}
\alias{search_anymatch}
\title{Search for any match}
\usage{
search_anymatch(x, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{x}{text or taxonomic serial number (TSN) (character or numeric)}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a data.frame
}
\description{
Search for any match
}
\examples{
\dontrun{
search_anymatch(x = 202385)
search_anymatch(x = "dolphin")
# no results
search_anymatch(x = "Pisces")
}
}
\seealso{
\code{\link{search_any_match_paged}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kingdoms.R
\name{kingdoms}
\alias{kingdoms}
\alias{kingdom_name}
\alias{kingdom_names}
\title{Get kingdom names from tsn}
\usage{
kingdom_name(tsn, wt = "json", raw = FALSE, ...)

kingdom_names(wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get kingdom names from tsn
}
\details{
\itemize{
\item kingdom_name: Get kingdom name for a TSN
\item kingdom_names: Get all possible kingdom names
}
}
\examples{
\dontrun{
kingdom_name(202385)
kingdom_name(202385, wt = "xml")
kingdom_names()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/date_data.R
\name{date_data}
\alias{date_data}
\title{Get date data from tsn}
\usage{
date_data(tsn, wt = "json", raw = FALSE, ...)
}
\arguments{
\item{tsn}{TSN for a taxonomic group (numeric). Required.}

\item{wt}{(character) One of "json" or "xml". Required.}

\item{raw}{(logical) Return raw JSON or XML as character string. Required.
Default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get date data from tsn
}
\examples{
\dontrun{
date_data(tsn = 180543)
date_data(180543, wt = "xml")
date_data(180543, wt = "json", raw = TRUE)
}
}
