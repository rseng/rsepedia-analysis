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
rcol
====

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("rcol")
```

Dev version

```{r eval=FALSE}
pak::pkg_install("ropensci/rcol")
# OR
install.packages("rcol", repos="https://dev.ropensci.org")
```

```{r}
library("rcol")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rcol/issues).
* License: MIT
* Get citation information for `rcol` in R doing `citation(package = 'rcol')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "rcol"
author: "Scott Chamberlain"
date: "2021-01-08"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to rcol}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`rcol` is a R client for the new Catalogue of Life (COL)

The API is completely different from the  old one, so don't expect
the same interface. The `taxize` package used to have some COL
functionality but was 
deprecated when the old COL API was just too unreliable. COL
has a lot of different routes, many of which aren't relavant to
us here (e.g., creating/editing data in COL, user profile actions,
etc.). 

Package documentation: https://docs.ropensci.org/rcol/

Check out COL API docs https://api.catalogue.life/ where you
can get more info on each API route and try the routes.

The following are a few examples.

## Installation


```r
remotes::install_github("ropensci/rcol")
# OR
install.packages("rcol", repos="https://dev.ropensci.org")
```


```r
library("rcol")
```

## Search

Two function `cp_nu_suggest()` and `cp_nu_search()` can be used for searching
for taxa. The former returns very minimal information and is meant to return
a response very quickly, while the latter is a bit slower but more comprehensive.


```r
cp_nu_suggest(q="Apis", dataset_key = 3)
#> # A tibble: 10 x 8
#>    match   parentOrAccepted… usageId    rank  status score suggestion    nomCode
#>    <chr>   <chr>             <chr>      <chr> <chr>  <dbl> <chr>         <chr>  
#>  1 Apisti… Scorpaeniformes   382e8442-… fami… accep…     0 Apistidae (f… <NA>   
#>  2 Apisto… Spionida          9d89b038-… fami… accep…     0 Apistobranch… <NA>   
#>  3 Apis    Apini             5eecbda1-… genus accep…     0 Apis (genus … zoolog…
#>  4 Apisa   Arctiidae         da6bee8f-… genus accep…     0 Apisa (genus… <NA>   
#>  5 Apispi… Mangeliidae       4a2734aa-… genus accep…     0 Apispiralia … <NA>   
#>  6 Apista… Notodontidae      2deb9907-… genus accep…     0 Apistaeschra… <NA>   
#>  7 Apisto… Apistobranchidae  9f5ef30a-… genus accep…     0 Apistobranch… <NA>   
#>  8 Apisto… Buthidae          e9ffbc3d-… genus accep…     0 Apistobuthus… zoolog…
#>  9 Apisto… Cichlidae         a3d969ca-… genus accep…     0 Apistogramma… <NA>   
#> 10 Apisto… Cichlidae         d6ded63a-… genus accep…     0 Apistogrammo… <NA>
cp_nu_search(q="Apis", rank = "genus")
#> $result
#> # A tibble: 10 x 5
#>    id    classification usage$created $createdBy $modified $modifiedBy
#>    <chr> <list>         <chr>              <int> <chr>           <int>
#>  1 1271… <df[,3] [8 × … 2020-07-02T1…         10 2020-07-…          10
#>  2 1647… <df[,3] [10 ×… 2020-07-02T0…         10 2020-07-…          10
#>  3 x6KS  <df[,3] [6 × … 2020-07-02T1…         10 2020-07-…          10
#>  4 x32YJ <df[,3] [6 × … 2020-07-03T0…         10 2020-07-…          10
#>  5 4e10… <df[,3] [8 × … 2020-03-18T2…        103 2020-03-…         103
#>  6 x3KBG <df[,3] [7 × … 2020-07-02T0…         10 2020-07-…          10
#>  7 1024… <df[,3] [9 × … 2020-07-01T2…         10 2020-07-…          10
#>  8 1005… <df[,3] [6 × … 2020-07-03T0…         10 2020-07-…          10
#>  9 553e… <df[,3] [10 ×… 2020-08-20T1…        103 2020-08-…         103
#> 10 4WXV8 <df[,3] [10 ×… 2020-10-28T1…        102 2020-10-…         102
#> # … with 38 more variables: $datasetKey <int>, $id <chr>, $verbatimKey <int>,
#> #   $name$created <chr>, $$createdBy <int>, $$modified <chr>,
#> #   $$modifiedBy <int>, $$datasetKey <int>, $$id <chr>, $$verbatimKey <int>,
#> #   $$homotypicNameId <chr>, $$scientificName <chr>, $$authorship <chr>,
#> #   $$rank <chr>, $$uninomial <chr>, $$combinationAuthorship$authors <list>,
#> #   $$$year <chr>, $$code <chr>, $$publishedInId <chr>, $$origin <chr>,
#> #   $$type <chr>, $$link <chr>, $$parsed <lgl>, $$sectorKey <int>,
#> #   $$nomStatus <chr>, $status <chr>, $origin <chr>, $parentId <chr>,
#> #   $link <chr>, $referenceIds <list>, $label <chr>, $labelHtml <chr>,
#> #   $sectorKey <int>, $scrutinizerDate <chr>, $extinct <lgl>,
#> #   $scrutinizer <chr>, issues <list>, sectorDatasetKey <int>
#> 
#> $meta
#> # A tibble: 1 x 5
#>   offset limit total last  empty
#>    <int> <int> <int> <lgl> <lgl>
#> 1      0    10  4343 FALSE FALSE
```

## Parsing names

`cp_parser()` parses scientitic names into their components. See also the 
package `rgnparser` (docs: https://ropensci.github.io/rgnparser/) for much
more powerul scientitic name parsing.


```r
cp_parser(names = c("Apis mellifera", "Homo sapiens var. sapiens"))
#> # A tibble: 2 x 7
#>   scientificName    rank   genus specificEpithet type   parsed infraspecificEpi…
#>   <chr>             <chr>  <chr> <chr>           <chr>  <lgl>  <chr>            
#> 1 Apis mellifera    speci… Apis  mellifera       scien… TRUE   <NA>             
#> 2 Homo sapiens var… varie… Homo  sapiens         scien… TRUE   sapiens
```

## Vocab

If you're curious about COL vocabularies see the function `cp_vocab()`


```r
cp_vocab("rank")
#>  [1] "domain"                "realm"                 "subrealm"             
#>  [4] "superkingdom"          "kingdom"               "subkingdom"           
#>  [7] "infrakingdom"          "superphylum"           "phylum"               
#> [10] "subphylum"             "infraphylum"           "superclass"           
#> [13] "class"                 "subclass"              "infraclass"           
#> [16] "parvclass"             "superdivision"         "division"             
#> [19] "subdivision"           "infradivision"         "superlegion"          
#> [22] "legion"                "sublegion"             "infralegion"          
#> [25] "supercohort"           "cohort"                "subcohort"            
#> [28] "infracohort"           "gigaorder"             "magnorder"            
#> [31] "grandorder"            "mirorder"              "superorder"           
#> [34] "order"                 "nanorder"              "hypoorder"            
#> [37] "minorder"              "suborder"              "infraorder"           
#> [40] "parvorder"             "megafamily"            "grandfamily"          
#> [43] "superfamily"           "epifamily"             "family"               
#> [46] "subfamily"             "infrafamily"           "supertribe"           
#> [49] "tribe"                 "subtribe"              "infratribe"           
#> [52] "suprageneric name"     "genus"                 "subgenus"             
#> [55] "infragenus"            "supersection"          "section"              
#> [58] "subsection"            "superseries"           "series"               
#> [61] "subseries"             "infrageneric name"     "species aggregate"    
#> [64] "species"               "infraspecific name"    "grex"                 
#> [67] "subspecies"            "cultivar group"        "convariety"           
#> [70] "infrasubspecific name" "proles"                "natio"                
#> [73] "aberration"            "morph"                 "variety"              
#> [76] "subvariety"            "form"                  "subform"              
#> [79] "pathovar"              "biovar"                "chemovar"             
#> [82] "morphovar"             "phagovar"              "serovar"              
#> [85] "chemoform"             "forma specialis"       "cultivar"             
#> [88] "strain"                "other"                 "unranked"
```

## Datasets

To search datasets, see the function `cp_datasets()`


```r
cp_datasets(q = "life")
#> $result
#> # A tibble: 10 x 28
#>    title   key alias group license   size created createdBy modified modifiedBy
#>    <chr> <int> <chr> <chr> <chr>    <int> <chr>       <int> <chr>         <int>
#>  1 Cata…  2242 COL2… Biota other   3.98e6 2020-1…       101 2020-12…        101
#>  2 Cata…     3 col-… Biota other   3.99e6 2019-1…       102 2020-12…        102
#>  3 Cata…  2230 COL2… Biota cc by   3.98e6 2020-1…       101 2020-11…        101
#>  4 Cata…  1000 ColH  <NA>  cc by   2.84e3 2019-1…       101 2020-10…        101
#>  5 Cata…  2123 CoL2… Biota cc by   4.52e6 2020-0…       101 2020-10…        101
#>  6 Cata…  2079 CoL2… Biota cc by   4.01e6 2020-0…       101 2020-10…        103
#>  7 Cata…  2165 CoL2… Biota cc by   3.95e6 2020-0…       103 2020-10…        103
#>  8 Cata…  2206 CoL2… Biota cc by   4.00e6 2020-0…       103 2020-10…        103
#>  9 Cata…  2083 CoL2… Biota cc by   4.22e6 2020-0…       101 2020-10…        101
#> 10 Cata…  2081 CoL2… Biota cc by   3.91e6 2020-0…       101 2020-10…        101
#> # … with 22 more variables: sourceKey <int>, importAttempt <int>, type <chr>,
#> #   origin <chr>, description <chr>, organisations <list>,
#> #   contact$givenName <chr>, $familyName <chr>, $email <chr>, $orcid <chr>,
#> #   $name <chr>, editors <list>, version <chr>, released <chr>, citation <chr>,
#> #   geographicScope <chr>, website <chr>, confidence <int>, completeness <int>,
#> #   imported <chr>, private <lgl>, authors <list>
#> 
#> $meta
#> # A tibble: 1 x 5
#>   offset limit total last  empty
#>    <int> <int> <int> <lgl> <lgl>
#> 1      0    10    93 FALSE FALSE
```

There are A LOT of datasets API routes. Instead of making
an R function for each route, we have R functions for some of the
"more important" routes, then `cp_ds()` will allow you to make
requests to the remainder of the datasets API routes. The `route`
parameter accepts the route as the examples below, with the 
preceding `/dataset` part of the route. Then pass in named parameters
to the function to fill in the templated route.


```r
cp_ds(route = "{key}/tree", key = "1000")
#> $offset
#> [1] 0
#> 
#> $limit
#> [1] 100
#> 
#> $total
#> [1] 2
#> 
#> $result
#>   datasetKey  id         rank   status childCount  labelHtml       name
#> 1       1000 343 superkingdom accepted          5  Eukaryota  Eukaryota
#> 2       1000   1 superkingdom accepted          2 Prokaryota Prokaryota
#> 
#> $last
#> [1] TRUE
#> 
#> $empty
#> [1] FALSE
cp_ds(route = "{key}/name/{id}", key = 1005, id = 100003)
#> $created
#> [1] "2020-07-02T23:51:27.62508"
#> 
#> $createdBy
#> [1] 10
#> 
#> $modified
#> [1] "2020-09-30T02:54:17.968034"
#> 
#> $modifiedBy
#> [1] 10
#> 
#> $datasetKey
#> [1] 1005
#> 
#> $id
#> [1] "100003"
#> 
#> $verbatimKey
#> [1] 1810
#> 
#> $homotypicNameId
#> [1] "100003"
#> 
#> $scientificName
#> [1] "Cylindrotoma aurantia"
#> 
#> $authorship
#> [1] "Alexander, 1935"
#> 
#> $rank
#> [1] "species"
#> 
#> $genus
#> [1] "Cylindrotoma"
#> 
#> $specificEpithet
#> [1] "aurantia"
#> 
#> $combinationAuthorship
#> $combinationAuthorship$authors
#> [1] "Alexander"
#> 
#> $combinationAuthorship$year
#> [1] "1935"
#> 
#> 
#> $code
#> [1] "zoological"
#> 
#> $origin
#> [1] "source"
#> 
#> $type
#> [1] "scientific"
#> 
#> $parsed
#> [1] TRUE
```

Alternatively, pass a named list to the `.list` parameter.


```r
args <- list(key = 1005, id = 100003)
cp_ds("{key}/name/{id}", .list = args)
#> $created
#> [1] "2020-07-02T23:51:27.62508"
#> 
#> $createdBy
#> [1] 10
#> 
#> $modified
#> [1] "2020-09-30T02:54:17.968034"
#> 
#> $modifiedBy
#> [1] 10
#> 
#> $datasetKey
#> [1] 1005
#> 
#> $id
#> [1] "100003"
#> 
#> $verbatimKey
#> [1] 1810
#> 
#> $homotypicNameId
#> [1] "100003"
#> 
#> $scientificName
#> [1] "Cylindrotoma aurantia"
#> 
#> $authorship
#> [1] "Alexander, 1935"
#> 
#> $rank
#> [1] "species"
#> 
#> $genus
#> [1] "Cylindrotoma"
#> 
#> $specificEpithet
#> [1] "aurantia"
#> 
#> $combinationAuthorship
#> $combinationAuthorship$authors
#> [1] "Alexander"
#> 
#> $combinationAuthorship$year
#> [1] "1935"
#> 
#> 
#> $code
#> [1] "zoological"
#> 
#> $origin
#> [1] "source"
#> 
#> $type
#> [1] "scientific"
#> 
#> $parsed
#> [1] TRUE
```

## Classification


```r
cp_classification(dataset_key=1000, taxon_id=10)
#> # A tibble: 4 x 16
#>   scientificName rank  id    status created createdBy modified modifiedBy
#>   <chr>          <chr> <chr> <chr>  <chr>       <int> <chr>         <int>
#> 1 Thermoprotei   class 6     accep… 2020-0…        10 2020-07…         10
#> 2 Crenarchaeota  phyl… 3     accep… 2020-0…        10 2020-07…         10
#> 3 Archaea        king… 2     accep… 2020-0…        10 2020-07…         10
#> 4 Prokaryota     supe… 1     accep… 2020-0…        10 2020-07…         10
#> # … with 8 more variables: datasetKey <int>, verbatimKey <int>,
#> #   homotypicNameId <chr>, uninomial <chr>, origin <chr>, type <chr>,
#> #   parsed <lgl>, parentId <chr>
```

## Children


```r
cp_children(dataset_key=1000, taxon_id='1')
#> $result
#> # A tibble: 2 x 16
#>   scientificName rank  id    status created createdBy modified modifiedBy
#>   <chr>          <chr> <chr> <chr>  <chr>       <int> <chr>         <int>
#> 1 Archaea        king… 2     accep… 2020-0…        10 2020-07…         10
#> 2 Bacteria       king… 41    accep… 2020-0…        10 2020-07…         10
#> # … with 8 more variables: datasetKey <int>, verbatimKey <int>,
#> #   homotypicNameId <chr>, uninomial <chr>, origin <chr>, type <chr>,
#> #   parsed <lgl>, parentId <chr>
#> 
#> $meta
#> # A tibble: 1 x 5
#>   offset limit total last  empty
#>    <int> <int> <int> <lgl> <lgl>
#> 1      0    10     2 TRUE  FALSE
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_ds.R
\name{cp_ds}
\alias{cp_ds}
\title{Datasets API route catch all method}
\usage{
cp_ds(route, ..., .list = list())
}
\arguments{
\item{route}{(character) an API route. the \verb{/dataset} route
part is added internally; so just include the route following that.
required.}

\item{...}{named parameters, passed on to \code{\link[glue:glue]{glue::glue()}}. required.
param names must match must match names given in the \code{route}. For
example, if you have \verb{route = \\\\\{key\\\\\}/name\\\\\{id\\\\\}}, then you need
to pass in a \code{key} and an \code{id} parameter. The names in the route
(here, key and id) don't have to match the names in the API route you
are trying to use - they just need to match the named parameters you
pass in. Having said that, it may be easier to remember what you're
doing if you match the names to the route parts.}

\item{.list}{a named list. instead of passing in named parameters
through \code{...}, you can pre-prepare a named list and give to this
parameter}
}
\value{
output varies depending on the route requested, but output will
always be a named list. when no results found, an error message
will be returned
}
\description{
Datasets API route catch all method
}
\details{
There are A LOT of datasets API routes. Instead of making
an R function for each route, we have R functions for some of the
"more important" routes, then \code{cp_ds()} will allow you to make
requests to the remainder of the datasets API routes.
}
\section{Not supported dataset routes}{

Some dataset routes do not return JSON so we don't support those.
Thus far, the only route we don't support is \verb{/dataset/\\\\\{key\\\\\}/logo}
}

\examples{
\dontrun{
cp_ds(route = "{key}/tree", key = "1000")
cp_ds(route = "{key}/tree", key = "1014")
cp_ds(route = "{key}/name/{id}", key = 1005, id = 100003)

# pass a named list to the .list parameter
args <- list(key = 1005, id = 100003)
cp_ds("{key}/name/{id}", .list = args)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_datasets.R
\name{cp_datasets}
\alias{cp_datasets}
\alias{cp_dataset}
\title{Datasets}
\usage{
cp_datasets(q = NULL, start = 0, limit = 10, ...)

cp_dataset(dataset_keys, ...)
}
\arguments{
\item{q}{(character) main query string. optional}

\item{start}{(integer) requested number of offset records. Default: 0}

\item{limit}{(integer) requested number of maximum records to be returned.
Default: 10; max: 1000}

\item{...}{curl options passed on to \code{\link[crul]{verb-GET}}}

\item{dataset_keys}{(character) one or more dataset keys. required}
}
\value{
list with two slots
\itemize{
\item \code{result} (data.frame/tibble): results, a zero row data.frame
if no results found
\item \code{meta} (data.frame/tibble): number of results found
}
}
\description{
Datasets
}
\details{
for \code{cp_dataset()}, separate http requests are made for each
dataset key. unfortunately, the output of \code{cp_dataset()} is a list for
each dataset key because the nested structure of the data is hard to
rectangularize
}
\examples{
if (cp_up("/dataset")) {
cp_datasets(limit = 1)
}
\dontrun{
cp_datasets(q = "life")
cp_dataset(dataset_keys = 1000)
cp_dataset(dataset_keys = c(3, 1000, 1014))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_importer.R
\name{cp_importer}
\alias{cp_importer}
\title{Importer metrics}
\usage{
cp_importer(
  dataset_key = NULL,
  state = NULL,
  running = FALSE,
  start = 0,
  limit = 10,
  ...
)
}
\arguments{
\item{dataset_key}{(character) a dataset key to filter by. optional}

\item{state}{(character) filter listed import metrics by their state,
e.g. the last failed import. one of: downloading, processing, inserting,
unchanged, finished, canceled, failed. optional}

\item{running}{(logical) if only a list of running imports should
be returned. default: \code{FALSE}. optional}

\item{start}{(integer) requested number of offset records. Default: 0}

\item{limit}{(integer) requested number of maximum records to be returned.
Default: 10; max: 1000}

\item{...}{curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
a named list, with slots \code{offset} (integer), \code{limit} (integer),
\code{total} (integer), \code{result} (list), \code{empty} (boolean),
and \code{last} (boolean). The \code{result} slot is a list itself, with any
number of results as named lists
}
\description{
Importer metrics
}
\examples{
if (cp_up("/importer?limit=1")) {
cp_importer(limit = 1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_nu_search.R
\name{cp_nu_search}
\alias{cp_nu_search}
\title{Name Usage: Search}
\usage{
cp_nu_search(
  q = NULL,
  dataset_key = NULL,
  min_rank = NULL,
  max_rank = NULL,
  content = NULL,
  highlight = NULL,
  reverse = NULL,
  fuzzy = NULL,
  type = NULL,
  nomstatus = NULL,
  status = NULL,
  issue = NULL,
  published_in = NULL,
  facet = NULL,
  sortBy = NULL,
  start = 0,
  limit = 10,
  ...
)
}
\arguments{
\item{q}{(character) vector of one or more scientific names}

\item{dataset_key}{(character) dataset key}

\item{min_rank, max_rank}{(character) filter by rank. one of: domain, superkingdom,
kingdom, subkingdom, infrakingdom, superphylum, phylum, subphylum,
infraphylum, superclass, class, subclass, infraclass, parvclass,
superlegion, legion, sublegion, infralegion, supercohort, cohort,
subcohort, infracohort, magnorder, superorder, grandorder, order,
suborder, infraorder, parvorder, superfamily, family, subfamily,
infrafamily, supertribe, tribe, subtribe, infratribe, suprageneric name,
genus, subgenus, infragenus, supersection, section, subsection, superseries,
series, subseries, infrageneric name, species aggregate, species,
infraspecific name, grex, subspecies, cultivar group, convariety,
infrasubspecific name, proles, natio, aberration, morph, variety,
subvariety, form, subform, pathovar, biovar, chemovar, morphovar, phagovar,
serovar, chemoform, forma specialis, cultivar, strain, other, unranked}

\item{content}{(character) one of: 'scientific_name' or 'authorship'}

\item{highlight}{(logical) \code{TRUE} or \code{FALSE}. default: \code{NULL}}

\item{reverse}{(logical) \code{TRUE} or \code{FALSE}. default: \code{NULL}}

\item{fuzzy}{(logical) \code{TRUE} or \code{FALSE}. default: \code{NULL}}

\item{type}{(character) one of: 'prefix', 'whole_words', 'exact'}

\item{nomstatus}{(character) filter by nomenclatural status. one of: ok,
unavailable, illegitimate, variant, conserved, rejected, doubtful, unevaluated}

\item{status}{(character) filter by taxonomic status. one of: accepted,
doubtful, ambiguous synonym}

\item{issue}{(character) filter by issue found}

\item{published_in}{(character) reference id to filter names by}

\item{facet}{(character) request a facet to be returned. one of:
dataset_key, rank, nom_status, status, issue, type, field. facet
limit default: 50}

\item{sortBy}{(character) one of: "relevance", "name", "taxonomic",
"index_name_id", or "native"}

\item{start}{(integer) requested number of offset records. Default: 0}

\item{limit}{(integer) requested number of maximum records to be returned.
Default: 10; max: 1000}

\item{...}{curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
list with two slots
\itemize{
\item \code{result} (data.frame/tibble): results, a zero row data.frame
if no results found
\item \code{meta} (data.frame/tibble): number of results found
}
}
\description{
Name Usage: Search
}
\examples{
if (cp_up("/nameusage/search?q=Apis")) {
cp_nu_search(q="Apis", limit = 1)
}
\dontrun{
cp_nu_search(q="Agapostemon")
cp_nu_search(q="Agapostemon", dataset_key = 3)
cp_nu_search(q="Agapostemon", min_rank = "genus")
cp_nu_search(q="Agapostemon", nomstatus = "doubtful")
cp_nu_search(q="Agapostemon", status = "accepted")
cp_nu_search(q="Bombus", facet = "rank")
cp_nu_search(q="Agapostemon", dataset_key = 3, hasField="uninomial")

x <- cp_nu_search(q="Poa")
x
x$result
x$result$usage
x$result$usage$name
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_children.R
\name{cp_children}
\alias{cp_children}
\title{Children}
\usage{
cp_children(dataset_key, taxon_id, ...)
}
\arguments{
\item{dataset_key}{(character/integer/numeric) dataset identifier}

\item{taxon_id}{(character/integer/numeric) taxon identifier}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
list with two slots
\itemize{
\item \code{result} (data.frame/tibble): results, a zero row data.frame
if no results found
\item \code{meta} (data.frame/tibble): number of results found
}
}
\description{
Children
}
\examples{
chk <- function(x) {
  z <- tryCatch(crul::ok(x), error = function(e) e)
  if (inherits(z, "error")) FALSE else z
}
if (chk("https://api.catalogueoflife.org/version")) {
z <- cp_children(dataset_key=1000, taxon_id='1')
z
z$result
if (NROW(z$result) > 0) {
z$result$scientificName
z$result$created
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcol-package.R
\docType{package}
\name{rcol-package}
\alias{rcol-package}
\alias{rcol}
\title{rcol}
\description{
Catalogue of Life (CoL) Client
}
\note{
CoL API docs: https://api.catalogueoflife.org/
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_parser.R
\name{cp_parser}
\alias{cp_parser}
\title{Name Parser}
\usage{
cp_parser(names, ...)
}
\arguments{
\item{names}{(character) one or more scientific names to parse}

\item{...}{curl options passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\value{
tibble, with one row for each parsed name
}
\description{
Name Parser
}
\examples{
\dontrun{
cp_parser(names = "Apis mellifera")
cp_parser(names = c("Apis mellifera", "Homo sapiens var. sapiens"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_name_match.R
\name{cp_name_match}
\alias{cp_name_match}
\title{Name Matching}
\usage{
cp_name_match(
  q = NULL,
  rank = NULL,
  code = NULL,
  trusted = NULL,
  ver_bose = NULL,
  start = 0,
  limit = 10,
  ...
)
}
\arguments{
\item{q}{(character) scientific name to match}

\item{rank}{(character) rank to restrict matches to. one of: domain,
realm, subrealm, superkingdom, kingdom, subkingdom, infrakingdom,
superphylum, phylum, subphylum, infraphylum, superclass, class,
subclass, infraclass, parvclass, superdivision, division, subdivision,
infradivision, superlegion, legion, sublegion, infralegion, supercohort,
cohort, subcohort, infracohort, gigaorder, magnorder, grandorder, mirorder,
superorder, order, nanorder, hypoorder, minorder, suborder, infraorder,
parvorder, megafamily, grandfamily, superfamily, epifamily, family,
subfamily, infrafamily, supertribe, tribe, subtribe, infratribe,
suprageneric_name, genus, subgenus, infragenus, supersection, section,
subsection, superseries, series, subseries, infrageneric_name,
species_aggregate, species, infraspecific_name, grex, subspecies,
cultivar_group, convariety, infrasubspecific_name, proles, natio,
aberration, morph, variety, subvariety, form, subform, pathovar, biovar,
chemovar, morphovar, phagovar, serovar, chemoform, forma_specialis,
cultivar, strain, other, unranked}

\item{code}{(character) nomenclatural code to restrict matches to. one of:
bacterial, botanical, cultivars, viral, zoological, phytosociological}

\item{trusted}{(logical) if \code{TRUE}, unmatched name will be inserted into
the names index. default: \code{FALSE}}

\item{ver_bose}{(logical) if \code{TRUE}, list alternatively considered name
matches. default: \code{FALSE}}

\item{start}{(integer) requested number of offset records. Default: 0}

\item{limit}{(integer) requested number of maximum records to be returned.
Default: 10; max: 1000}

\item{...}{curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
a named list, with slots \code{name} (list), \code{type} (character),
\code{alternatives} (data.frame), and \code{nameKey} (integer)
}
\description{
Match name against the name index
}
\details{
Matches by the canonical name, it's authorship and rank.
Authorship matching is somewhat loose, but name matching is quite strict
and only allows for a few common misspellings frequently found in
epithets (silent h, gender suffix, double letters, i/y), but not in
uninomials. Suprageneric ranks are all considered to be the same,
otherwise a different rank results in a different match.
}
\examples{
if (cp_up("/name/matching?q=Apis")) { 
cp_name_match(q="Apis")
}
\dontrun{
cp_name_match(q="Agapostemon")
cp_name_match(q="Apis mellifera")
cp_name_match(q="Apis mellifer") # no fuzzy match apparently
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_vocab.R
\name{cp_vocab}
\alias{cp_vocab}
\title{CoL Vocabularies}
\usage{
cp_vocab(vocab, ...)
}
\arguments{
\item{vocab}{(character) a vocabulary name}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
character vector of words
}
\description{
CoL Vocabularies
}
\examples{
\dontrun{
cp_vocab("rank")
cp_vocab("datasetorigin")
cp_vocab("datasettype")
cp_vocab("matchtype")
cp_vocab("taxonomicstatus")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{cp_up}
\alias{cp_up}
\title{check if api up in examples}
\usage{
cp_up(x)
}
\value{
boolean. \code{TRUE} if http request succeeds, \code{FALSE} if not
}
\description{
check if api up in examples
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_classification.R
\name{cp_classification}
\alias{cp_classification}
\title{Classification}
\usage{
cp_classification(dataset_key, taxon_id, ...)
}
\arguments{
\item{dataset_key}{(character/integer/numeric) dataset identifier}

\item{taxon_id}{(character/integer/numeric) taxon identifier}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a data.frame/tibble with results, a zero row data.frame
if no results found
}
\description{
Classification
}
\examples{
if (cp_up("/dataset/1000/taxon/10/classification")) {
cp_classification(dataset_key=1000, taxon_id=10)
}
\dontrun{
cp_classification(dataset_key=1000, taxon_id=20)
cp_classification(dataset_key=3,
 taxon_id="6565450e-1cf2-4dc2-acbb-db728e42e635")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_nu_suggest.R
\name{cp_nu_suggest}
\alias{cp_nu_suggest}
\title{Name Usage: Suggest}
\usage{
cp_nu_suggest(
  q,
  dataset_key,
  fuzzy = FALSE,
  min_rank = NULL,
  max_rank = NULL,
  sort = NULL,
  reverse = FALSE,
  accepted = FALSE,
  limit = 10,
  ...
)
}
\arguments{
\item{q}{(character) main query string. required}

\item{dataset_key}{(character) a dataset key. required}

\item{fuzzy}{(logical) Whether or not to do fuzzy search (default: \code{FALSE})}

\item{min_rank, max_rank}{(character) See rank options in \code{\link[=cp_name_match]{cp_name_match()}}}

\item{sort}{(character) one of name, taxonomic, index_name_id, native,
relevance}

\item{reverse}{(logical) reverse order i assume (default: \code{FALSE})}

\item{accepted}{(logical) limit to accepted names (default: \code{FALSE})}

\item{limit}{(integer) requested number of maximum records to be returned.
Default: 10; max: 1000}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a data.frame/tibble of results. a zero row data.frame if no results
}
\description{
Name Usage: Suggest
}
\examples{
if (cp_up("/dataset/3/nameusage/suggest?q=Apis")) {
cp_nu_suggest(q="Apis", 3)
}
}
