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

*Wow, no problems at all. :)**Wow, no problems at all. :)*worrms
======

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

```{r eval=FALSE}
install.packages("worrms")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/worrms")
```

```{r}
library("worrms")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/worrms/issues).
* License: MIT
* Get citation information for `worrms` in R doing `citation(package = 'worrms')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: Introduction to worrms
author: Scott Chamberlain
date: "2020-07-07"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Introduction to worrms}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



`worrms` is an R client for the World Register of Marine Species (http://www.marinespecies.org/)

See the taxize book (https://taxize.dev) for taxonomically focused work in this and similar packages.

## Install

Stable version from CRAN


```r
install.packages("worrms")
```

Development version from GitHub


```r
install.packages("remotes")
remotes::install_github("ropensci/worrms")
```


```r
library("worrms")
```

## Get records

WoRMS 'records' are taxa, not specimen occurrences or something else.

by date


```r
wm_records_date('2016-12-23T05:59:45+00:00')
#> # A tibble: 50 x 27
#>    AphiaID url   scientificname authority status unacceptreason taxonRankID
#>      <int> <chr> <chr>          <chr>     <chr>  <lgl>                <int>
#>  1  894302 http… Paleopolymorp… Vasilenk… accep… NA                     220
#>  2  894296 http… Parapachyphlo… Miklukho… accep… NA                     220
#>  3  894298 http… Parapachyphlo… Miklukho… accep… NA                     220
#>  4  894301 http… Ovulina radia… Seguenza… accep… NA                     220
#>  5  894299 http… Parafissurina… Petri, 1… accep… NA                     220
#>  6  894297 http… Parapachyphlo… Miklukho… accep… NA                     220
#>  7  919684 http… Flintina serr… Cuvillie… accep… NA                     220
#>  8  906465 http… Verneuilinoid… Bartenst… accep… NA                     220
#>  9  903640 http… Quinqueloculi… Hussey, … accep… NA                     220
#> 10  917580 http… Orbitolina ra… Sahni & … accep… NA                     220
#> # … with 40 more rows, and 20 more variables: rank <chr>, valid_AphiaID <int>,
#> #   valid_name <chr>, valid_authority <chr>, parentNameUsageID <int>,
#> #   kingdom <chr>, phylum <chr>, class <chr>, order <chr>, family <chr>,
#> #   genus <chr>, citation <chr>, lsid <chr>, isMarine <int>, isBrackish <lgl>,
#> #   isFreshwater <lgl>, isTerrestrial <lgl>, isExtinct <int>, match_type <chr>,
#> #   modified <chr>
```

by a taxonomic name


```r
wm_records_name(name = 'Leucophaeus scoresbii')
#> # A tibble: 1 x 27
#>   AphiaID url   scientificname authority status unacceptreason taxonRankID rank 
#>     <int> <chr> <chr>          <chr>     <chr>  <lgl>                <int> <chr>
#> 1  344089 http… Leucophaeus s… Traill, … accep… NA                     220 Spec…
#> # … with 19 more variables: valid_AphiaID <int>, valid_name <chr>,
#> #   valid_authority <chr>, parentNameUsageID <int>, kingdom <chr>,
#> #   phylum <chr>, class <chr>, order <chr>, family <chr>, genus <chr>,
#> #   citation <chr>, lsid <chr>, isMarine <int>, isBrackish <lgl>,
#> #   isFreshwater <lgl>, isTerrestrial <lgl>, isExtinct <lgl>, match_type <chr>,
#> #   modified <chr>
```

by many names


```r
wm_records_names(name = c('Leucophaeus scoresbii', 'Coryphaena'))
#> [[1]]
#> # A tibble: 1 x 27
#>   AphiaID url   scientificname authority status unacceptreason taxonRankID rank 
#>     <int> <chr> <chr>          <chr>     <chr>  <lgl>                <int> <chr>
#> 1  344089 http… Leucophaeus s… Traill, … accep… NA                     220 Spec…
#> # … with 19 more variables: valid_AphiaID <int>, valid_name <chr>,
#> #   valid_authority <chr>, parentNameUsageID <int>, kingdom <chr>,
#> #   phylum <chr>, class <chr>, order <chr>, family <chr>, genus <chr>,
#> #   citation <chr>, lsid <chr>, isMarine <int>, isBrackish <lgl>,
#> #   isFreshwater <lgl>, isTerrestrial <lgl>, isExtinct <lgl>, match_type <chr>,
#> #   modified <chr>
#> 
#> [[2]]
#> # A tibble: 2 x 27
#>   AphiaID url   scientificname authority status unacceptreason taxonRankID rank 
#>     <int> <chr> <chr>          <chr>     <chr>  <chr>                <int> <chr>
#> 1  125960 http… Coryphaena     Linnaeus… accep… <NA>                   180 Genus
#> 2  843430 <NA>  <NA>           <NA>      quara… synonym                 NA <NA> 
#> # … with 19 more variables: valid_AphiaID <int>, valid_name <chr>,
#> #   valid_authority <chr>, parentNameUsageID <int>, kingdom <chr>,
#> #   phylum <chr>, class <chr>, order <chr>, family <chr>, genus <chr>,
#> #   citation <chr>, lsid <chr>, isMarine <int>, isBrackish <int>,
#> #   isFreshwater <int>, isTerrestrial <int>, isExtinct <lgl>, match_type <chr>,
#> #   modified <chr>
```

by common name


```r
wm_records_common(name = 'clam')
#> # A tibble: 4 x 27
#>   AphiaID url   scientificname authority status unacceptreason taxonRankID rank 
#>     <int> <chr> <chr>          <chr>     <chr>  <lgl>                <int> <chr>
#> 1  141919 http… Mercenaria me… (Linnaeu… accep… NA                     220 Spec…
#> 2  140431 http… Mya truncata   Linnaeus… accep… NA                     220 Spec…
#> 3  141936 http… Venus verruco… Linnaeus… accep… NA                     220 Spec…
#> 4  575771 http… Verpa penis    (Linnaeu… accep… NA                     220 Spec…
#> # … with 19 more variables: valid_AphiaID <int>, valid_name <chr>,
#> #   valid_authority <chr>, parentNameUsageID <int>, kingdom <chr>,
#> #   phylum <chr>, class <chr>, order <chr>, family <chr>, genus <chr>,
#> #   citation <chr>, lsid <chr>, isMarine <int>, isBrackish <lgl>,
#> #   isFreshwater <lgl>, isTerrestrial <lgl>, isExtinct <lgl>, match_type <chr>,
#> #   modified <chr>
```

using the TAXMATCH algorithm


```r
wm_records_taxamatch(name = 'Leucophaeus scoresbii')
#> [[1]]
#> # A tibble: 1 x 27
#>   AphiaID url   scientificname authority status unacceptreason taxonRankID rank 
#>     <int> <chr> <chr>          <chr>     <chr>  <lgl>                <int> <chr>
#> 1  344089 http… Leucophaeus s… Traill, … accep… NA                     220 Spec…
#> # … with 19 more variables: valid_AphiaID <int>, valid_name <chr>,
#> #   valid_authority <chr>, parentNameUsageID <int>, kingdom <chr>,
#> #   phylum <chr>, class <chr>, order <chr>, family <chr>, genus <chr>,
#> #   citation <chr>, lsid <chr>, isMarine <int>, isBrackish <lgl>,
#> #   isFreshwater <lgl>, isTerrestrial <lgl>, isExtinct <lgl>, match_type <chr>,
#> #   modified <chr>
```

## APHIA ID <--> name


```r
wm_name2id(name = "Rhincodon")
#> [1] 105749
```


```r
wm_id2name(id = 105706)
#> [1] "Rhincodontidae"
```

## Get AphiaID via an external ID


```r
wm_external(id = 1080)
#> [1] 85257
wm_external(id = 105706)
#> [1] 159854
```

## Get vernacular names from an AphiaID


```r
wm_common_id(id = 156806)
#> # A tibble: 2 x 3
#>   vernacular          language_code language
#>   <chr>               <chr>         <chr>   
#> 1 gilded wedgeclam    eng           English 
#> 2 Turton's wedge clam eng           English
```

## Children

Get direct taxonomic children for an AphiaID


```r
wm_classification(id = 105706)
#> # A tibble: 11 x 3
#>    AphiaID rank       scientificname  
#>      <int> <chr>      <chr>           
#>  1       2 Kingdom    Animalia        
#>  2    1821 Phylum     Chordata        
#>  3  146419 Subphylum  Vertebrata      
#>  4    1828 Superclass Gnathostomata   
#>  5   11676 Superclass Pisces          
#>  6   10193 Class      Elasmobranchii  
#>  7  368407 Subclass   Neoselachii     
#>  8  368408 Infraclass Selachii        
#>  9  368410 Superorder Galeomorphi     
#> 10   10208 Order      Orectolobiformes
#> 11  105706 Family     Rhincodontidae
```

## Classification

Get classification for an AphiaID


```r
wm_classification(id = 105706)
#> # A tibble: 11 x 3
#>    AphiaID rank       scientificname  
#>      <int> <chr>      <chr>           
#>  1       2 Kingdom    Animalia        
#>  2    1821 Phylum     Chordata        
#>  3  146419 Subphylum  Vertebrata      
#>  4    1828 Superclass Gnathostomata   
#>  5   11676 Superclass Pisces          
#>  6   10193 Class      Elasmobranchii  
#>  7  368407 Subclass   Neoselachii     
#>  8  368408 Infraclass Selachii        
#>  9  368410 Superorder Galeomorphi     
#> 10   10208 Order      Orectolobiformes
#> 11  105706 Family     Rhincodontidae
```

## Synonyms

Get synonyms for an AphiaID


```r
wm_synonyms(id = 105706)
#> # A tibble: 1 x 27
#>   AphiaID url   scientificname authority status unacceptreason taxonRankID rank 
#>     <int> <chr> <chr>          <chr>     <chr>  <chr>                <int> <chr>
#> 1  148832 http… Rhiniodontidae Müller &… unacc… synonym                140 Fami…
#> # … with 19 more variables: valid_AphiaID <int>, valid_name <chr>,
#> #   valid_authority <chr>, parentNameUsageID <int>, kingdom <chr>,
#> #   phylum <chr>, class <chr>, order <chr>, family <chr>, genus <lgl>,
#> #   citation <chr>, lsid <chr>, isMarine <lgl>, isBrackish <lgl>,
#> #   isFreshwater <lgl>, isTerrestrial <lgl>, isExtinct <lgl>, match_type <chr>,
#> #   modified <chr>
```

## attributes (i.e., traits)

attribute definition by ID


```r
wm_attr_def(id = 1)
#> # A tibble: 1 x 4
#>   measurementTypeID measurementType        CategoryID children        
#>               <int> <chr>                       <int> <list>          
#> 1                 1 IUCN Red List Category          1 <df[,4] [2 × 4]>
```

attribute data by AphiaID


```r
wm_attr_data(id = 127160)
#> # A tibble: 24 x 10
#>    AphiaID measurementType… measurementType measurementValue source_id reference
#>    <chr>              <int> <chr>           <chr>                <int> <chr>    
#>  1 127160                23 Species import… FAO-ASFIS: Spec…    197354 "FAO Fis…
#>  2 127160                23 Species import… MSFD indicators     197546 "Daniel …
#>  3 127160                23 Species import… MSFD indicators     197549 "ICES. 2…
#>  4 127160                23 Species import… MSFD indicators     197615 "List of…
#>  5 127160                23 Species import… MSFD indicators     197615 "List of…
#>  6 127160                23 Species import… MSFD indicators     197615 "List of…
#>  7 127160                23 Species import… MSFD indicators     197615 "List of…
#>  8 127160                23 Species import… MSFD indicators     197616 "List of…
#>  9 127160                23 Species import… MSFD indicators     197616 "List of…
#> 10 127160                23 Species import… MSFD indicators     197549 "ICES. 2…
#> # … with 14 more rows, and 4 more variables: qualitystatus <chr>,
#> #   AphiaID_Inherited <int>, CategoryID <int>, children <list>
```

attributes grouped by a CategoryID


```r
wm_attr_category(id = 7)
#> # A tibble: 6 x 4
#>   measurementValueID measurementValue measurementValueCode children        
#>                <int> <chr>            <chr>                <list>          
#> 1                183 benthos          <NA>                 <df[,4] [8 × 4]>
#> 2                184 plankton         <NA>                 <df[,4] [7 × 4]>
#> 3                194 nekton           <NA>                 <df[,0] [0 × 0]>
#> 4                323 neuston          <NA>                 <df[,0] [0 × 0]>
#> 5                378 edaphofauna      <NA>                 <df[,4] [2 × 4]>
#> 6                331 not applicable   N/A                  <df[,0] [0 × 0]>
```

AphiaIDs by attribute definition ID


```r
wm_attr_aphia(id = 4)
#> # A tibble: 50 x 2
#>    AphiaID Attributes        
#>      <int> <list>            
#>  1      11 <df[,10] [1 × 10]>
#>  2      55 <df[,10] [2 × 10]>
#>  3      57 <df[,10] [2 × 10]>
#>  4      58 <df[,10] [2 × 10]>
#>  5      59 <df[,10] [2 × 10]>
#>  6      63 <df[,10] [2 × 10]>
#>  7      64 <df[,10] [2 × 10]>
#>  8      69 <df[,10] [2 × 10]>
#>  9      90 <df[,10] [2 × 10]>
#> 10      91 <df[,10] [2 × 10]>
#> # … with 40 more rows
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_ranks.R
\name{wm_ranks}
\alias{wm_ranks}
\alias{wm_ranks_id}
\alias{wm_ranks_name}
\title{Get taxonomic ranks by their identifier}
\usage{
wm_ranks_id(rank_id, id = NULL, offset = 1, ...)

wm_ranks_name(rank_name, id = NULL, offset = 1, ...)
}
\arguments{
\item{rank_id}{(numeric/integer) a rank identifier. length==1}

\item{id}{an AphiaID. length==1}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{rank_name}{(character) a rank name. length==1}
}
\value{
A tibble/data.frame
}
\description{
Get taxonomic ranks by their identifier
}
\examples{
\dontrun{
wm_ranks_id(220)
wm_ranks_id(180)
wm_ranks_id(180, id = 4)

wm_ranks_name("genus")
wm_ranks_name("genus", id = 4)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_record.R
\name{wm_record}
\alias{wm_record}
\alias{wm_record_}
\title{Get complete AphiaRecord for an AphiaID}
\usage{
wm_record(id, ...)

wm_record_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_record} it's
required and must be \code{length(id) == 1}, for \code{wm_record_} it's
optional and can be \code{length(id) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A named list. When using underscore method, each output is named
by the input ID, and can be separated by the list names
}
\description{
Get complete AphiaRecord for an AphiaID
}
\note{
\code{wm_record_} is defunct, \code{wm_record} can do plural requests now
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_record(id = 105706)
wm_record(id = c(105706, 126436))
wm_record_(id = c(105706, 126436))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_record_by_external.R
\name{wm_record_by_external}
\alias{wm_record_by_external}
\alias{wm_record_by_external_}
\title{Get record by external ID}
\usage{
wm_record_by_external(id, type = "tsn", ...)

wm_record_by_external_(id = NULL, name = NULL, type = "tsn", ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_record_by_external}
it's required and must be \code{length(id) == 1}, for
\code{wm_record_by_external_} it's optional and can be \code{length(id) >= 1}}

\item{type}{(character) the type of external id. one of: tsn, bold,
dyntaxa, eol, fishbase, iucn, lsid, ncbi, gisd. default: tsn}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A named list. When using underscore method, each output is named
by the input ID, and can be separated by the list names
}
\description{
Get record by external ID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_record_by_external(id = 85257)
wm_record_by_external(id = 159854)

wm_record_by_external_(id = c(85257, 159854))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_records_name.R
\name{wm_records_name}
\alias{wm_records_name}
\title{Get records by single name, optional fuzzy matching}
\usage{
wm_records_name(name, fuzzy = TRUE, marine_only = TRUE, offset = 1, ...)
}
\arguments{
\item{name}{(character) a taxonomic name, required.}

\item{fuzzy}{(logical) fuzzy search. default: \code{TRUE}}

\item{marine_only}{(logical) marine only or not. default: \code{TRUE}}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
A tibble/data.frame
}
\description{
Get records by single name, optional fuzzy matching
}
\note{
there is no underscore method like other functions in this package
as there is already a plural version: \code{\link[=wm_records_names]{wm_records_names()}}
}
\examples{
\dontrun{
wm_records_name(name = 'Leucophaeus')
wm_records_name(name = 'Leucophaeus', fuzzy = FALSE)
wm_records_name(name = 'Leucophaeus', marine_only = FALSE)
wm_records_name(name = 'Platanista', marine_only = FALSE)
wm_records_name(name = 'Platanista', marine_only = FALSE, offset = 5)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_attr_aphia.R
\name{wm_attr_aphia}
\alias{wm_attr_aphia}
\alias{wm_attr_aphia_}
\title{Get AphiaIDs by attribute definition ID}
\usage{
wm_attr_aphia(id, offset = 1, ...)

wm_attr_aphia_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) a attribute ID. For \code{wm_attr_aphia} it's
required and must be \code{length(id) == 1}, for \code{wm_attr_aphia_} it's
optional and can be \code{length(id) >= 1}}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get AphiaIDs by attribute definition ID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_attr_aphia(id = 7)
wm_attr_aphia(id = 4)
wm_attr_aphia(id = 4, offset = 50)

wm_attr_aphia_(id = c(7, 2))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_records_taxamatch.R
\name{wm_records_taxamatch}
\alias{wm_records_taxamatch}
\title{Get records for onen or more taxonomic name(s) using
the TAXAMATCH fuzzy matching algorithm}
\usage{
wm_records_taxamatch(name, marine_only = TRUE, ...)
}
\arguments{
\item{name}{(character) taxon name. required.}

\item{marine_only}{(logical) marine only or not. default: \code{TRUE}}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
A list of tibble's/data.frame's, one for each of the input names
}
\description{
Get records for onen or more taxonomic name(s) using
the TAXAMATCH fuzzy matching algorithm
}
\note{
there is no underscore method like other functions in this package
as this function already accepts many names
}
\examples{
\dontrun{
wm_records_taxamatch(name = 'Leucophaeus')
wm_records_taxamatch(name = c('Leucophaeus', 'Coryphaena'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_records_date.R
\name{wm_records_date}
\alias{wm_records_date}
\title{Get records by date}
\usage{
wm_records_date(
  start_date,
  end_date = NULL,
  marine_only = TRUE,
  offset = 1,
  ...
)
}
\arguments{
\item{start_date}{(character) start date. required.}

\item{end_date}{(character) end date. optional}

\item{marine_only}{(logical) marine only or not. default: \code{TRUE}}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
A tibble/data.frame
}
\description{
Get records by date
}
\examples{
\dontrun{
a_date <- format(Sys.Date() - 1, "\%Y-\%m-\%dT\%H:\%M:\%S+00:00")
wm_records_date(a_date)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_records_common.R
\name{wm_records_common}
\alias{wm_records_common}
\alias{wm_records_common_}
\title{Get records by vernacular name, optional fuzzy matching}
\usage{
wm_records_common(name, fuzzy = FALSE, offset = 1, ...)

wm_records_common_(name, fuzzy = FALSE, offset = 1, ...)
}
\arguments{
\item{name}{(character) a species common name. required. For
\code{wm_records_common} must be \code{length(name) == 1}; for \code{wm_records_common_}
can be \code{length(name) >= 1}}

\item{fuzzy}{(logical) fuzzy search. default: \code{FALSE}}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get records by vernacular name, optional fuzzy matching
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_records_common(name = 'dolphin')
wm_records_common(name = 'clam')

wm_records_common_(name = c('dolphin', 'clam'))

wm_records_common(name = 'dolphin', fuzzy = TRUE)
wm_records_common(name = 'clam', fuzzy = TRUE, offset = 5)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_distribution.R
\name{wm_distribution}
\alias{wm_distribution}
\alias{wm_distribution_}
\title{Get distribution data by AphiaID}
\usage{
wm_distribution(id, ...)

wm_distribution_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_distribution} it's
required and must be \code{length(id) == 1}, for \code{wm_distribution_} it's
optional and can be \code{length(id) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get distribution data by AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_distribution(id = 156806)
wm_distribution(id = 126436)

wm_distribution_(id = c(156806, 126436))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_attr_category.R
\name{wm_attr_category}
\alias{wm_attr_category}
\alias{wm_attr_category_}
\title{Get attributes grouped by a CategoryID}
\usage{
wm_attr_category(id, ...)

wm_attr_category_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) a CategoryID. For \code{wm_attr_category} it's
required and must be \code{length(id) == 1}, for \code{wm_attr_category_} it's
optional and can be \code{length(id) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get attributes grouped by a CategoryID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_attr_category(id = 7)
wm_attr_category(id = 2)

wm_attr_category_(id = c(7, 2))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_externalid.R
\name{wm_external}
\alias{wm_external}
\alias{wm_external_}
\title{Get an external ID via an AphiaID}
\usage{
wm_external(id, type = "tsn", ...)

wm_external_(id = NULL, name = NULL, type = "tsn", ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_external} it's
required and must be \code{length(id) == 1}, for \code{wm_external_} it's
optional and can be \code{length(id) >= 1}}

\item{type}{(character) the type of external id. one of: tsn, bold,
dyntaxa, eol, fishbase, iucn, lsid, ncbi, gisd. default: tsn}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
An integer that is the ID. When using underscore method,
a list, named by the input IDs
}
\description{
Get an external ID via an AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
# by default, get a TSN (an ITIS code)
wm_external(id = 1080)

## get many
wm_external_(id = c(1080, 126436))

# BOLD code
wm_external(id = 278468, type = "bold")

# NCBI code
wm_external(id = 278468, type = "ncbi")

# fishbase code
wm_external(id = 278468, type = "fishbase")

# curl options
library(crul)
wm_external(id = 105706, verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_attr_data.R
\name{wm_attr_data}
\alias{wm_attr_data}
\alias{wm_attr_data_}
\title{Get attribute data by AphiaID}
\usage{
wm_attr_data(id, include_inherited = FALSE, ...)

wm_attr_data_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_attr_data} it's
required and must be \code{length(id) == 1}, for \code{wm_attr_data_} it's
optional and can be \code{length(id) >= 1}}

\item{include_inherited}{(logical) Include attributes inherited from
its parent taxon. Default: \code{FALSE}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get attribute data by AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_attr_data(id = 127160)
wm_attr_data(id = 126436)

wm_attr_data_(id = c(127160, 126436))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_name2id.R
\name{wm_name2id}
\alias{wm_name2id}
\alias{wm_name2id_}
\title{Get AphiaID from a taxonomic name}
\usage{
wm_name2id(name, ...)

wm_name2id_(name, ...)
}
\arguments{
\item{name}{(character) a taxonomic name, required. For
\code{wm_name2id} must be \code{length(name) == 1}, but for \code{wm_name2id_}
can be \code{length(name) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
An integer that is the AphiaID. When using underscore method,
a list, named by the input names
}
\description{
Get AphiaID from a taxonomic name
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_name2id(name = "Rhincodon")
wm_name2id_(name = c("Rhincodon", "Gadus morhua"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_records_names.R
\name{wm_records_names}
\alias{wm_records_names}
\title{Get records for onen or more taxonomic name(s)}
\usage{
wm_records_names(name, fuzzy = FALSE, marine_only = TRUE, ...)
}
\arguments{
\item{name}{(character) start date. required.}

\item{fuzzy}{(logical) fuzzy search. default: \code{FALSE}}

\item{marine_only}{(logical) marine only or not. default: \code{TRUE}}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
A list of tibble's/data.frame's, one for each of the input names
}
\description{
Get records for onen or more taxonomic name(s)
}
\note{
there is no underscore method like other functions in this package
as this is the plural version for \code{\link[=wm_records_name]{wm_records_name()}}
}
\examples{
\dontrun{
wm_records_names(name = 'Leucophaeus scoresbii')
wm_records_names(name = 'Leucophaeus scoresbii', fuzzy = TRUE)
wm_records_names(name = c('Leucophaeus scoresbii', 'Coryphaena'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_records_rank.R
\name{wm_records_rank}
\alias{wm_records_rank}
\title{Get AphiaRecords for a given taxonRankID}
\usage{
wm_records_rank(rank_id, id = NULL, offset = 1, ...)
}
\arguments{
\item{rank_id}{(numeric/integer) a rank id}

\item{id}{(character) a single AphiaID}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get AphiaRecords for a given taxonRankID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_records_rank(rank_id = 180, id = 106776)
wm_records_rank(rank_id = 180, id = 106776, offset = 50)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_classification.R
\name{wm_classification}
\alias{wm_classification}
\alias{wm_classification_}
\title{Get classification for an AphiaID}
\usage{
wm_classification(id, ...)

wm_classification_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_children} it's
required and must be \code{length(id) == 1}, for \code{wm_children_} it's
optional and can be \code{length(id) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get classification for an AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_classification(id = 105706)
wm_classification(id = 126436)

wm_classification(254967)
wm_classification(344089)

# plural version, via id or name
wm_classification_(id = c(254967, 344089))
wm_classification_(name = c('Platanista gangetica', 'Leucophaeus scoresbii'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_common_id.R
\name{wm_common_id}
\alias{wm_common_id}
\alias{wm_common_id_}
\title{Get vernacular names from an AphiaID}
\usage{
wm_common_id(id, ...)

wm_common_id_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_common_id} it's
required and must be \code{length(id) == 1}, for \code{wm_common_id_} it's
optional and can be \code{length(id) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get vernacular names from an AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_common_id(id = 105706)
wm_common_id(id = 156806)
wm_common_id(id = 397065)

wm_common_id_(id = c(105706, 156806, 397065))
nms <- c("Rhincodontidae", "Mesodesma deauratum", "Cryptomya californica")
wm_common_id_(name = nms)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_id2name.R
\name{wm_id2name}
\alias{wm_id2name}
\alias{wm_id2name_}
\title{Get taxonomic name for an AphiaID}
\usage{
wm_id2name(id, ...)

wm_id2name_(id, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID, required. For \code{wm_id2name}
must be \code{length(id) == 1}, but for \code{wm_id2name_} can be
\code{length(id) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}
}
\value{
An character string that is the taxnomic name. When using underscore
method, a list, named by the input IDs
}
\description{
Get taxonomic name for an AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_id2name(id = 105706)
wm_id2name_(id = c(105706, 126436))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_synonyms.R
\name{wm_synonyms}
\alias{wm_synonyms}
\alias{wm_synonyms_}
\title{Get synonyms for an AphiaID}
\usage{
wm_synonyms(id, offset = 1, ...)

wm_synonyms_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_synonyms} it's required
and must be \code{length(id) == 1}, for \code{wm_synonyms_} it's optional and
can be \code{length(id) >= 1}}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get synonyms for an AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_synonyms(id = 105706)
wm_synonyms_(id = 105706)
wm_synonyms(id = 126436)
wm_synonyms(id = 126436, offset = 10)
wm_synonyms_(id = c(105706, 126436))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/worrms-package.R
\docType{package}
\name{worrms-package}
\alias{worrms-package}
\alias{worrms}
\title{worrms}
\description{
World Register of Marine Species Client
}
\section{Fail behavior}{

The WoRMS REST API doesn't have sophisticated error messaging, so
most errors will result in a \verb{(204) - No Content} or
in \verb{(400) - Bad Request}

Because WoRMS doesn't do comprehensive error reporting, we do a fair
amount of checking user inputs to help prevent errors that will be
meaningless to the user. Let us know if we can improve on this.
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_sources.R
\name{wm_sources}
\alias{wm_sources}
\alias{wm_sources_}
\title{Get sources for an AphiaID}
\usage{
wm_sources(id, ...)

wm_sources_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_sources} it's required
and must be \code{length(id) == 1}, for \code{wm_sources_} it's optional and
can be \code{length(id) >= 1}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get sources for an AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_sources(id = 105706)
wm_sources_(id = 105706)
wm_sources_(id = c(105706, 126436))
wm_sources_(name = c("Rhincodontidae", "Gadus morhua"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_attr_def.R
\name{wm_attr_def}
\alias{wm_attr_def}
\alias{wm_attr_def_}
\title{Get attribute definition by ID}
\usage{
wm_attr_def(id, include_inherited = FALSE, ...)

wm_attr_def_(id = NULL, name = NULL, ...)
}
\arguments{
\item{id}{(numeric/integer) an attribute ID. For \code{wm_attr_def} it's
required and must be \code{length(id) == 1}, for \code{wm_attr_def_} it's
optional and can be \code{length(id) >= 1}}

\item{include_inherited}{(logical) Include attributes inherited from
its parent taxon. Default: \code{FALSE}}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get attribute definition by ID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_attr_def(id = 1)
wm_attr_def(id = 4)
wm_attr_def(id = 4, include_inherited = TRUE)

wm_attr_def_(id = c(4, 1))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wm_children.R
\name{wm_children}
\alias{wm_children}
\alias{wm_children_}
\title{Get children for an AphiaID}
\usage{
wm_children(id, marine_only = TRUE, offset = 1, ...)

wm_children_(id = NULL, name = NULL, marine_only = TRUE, offset = 1, ...)
}
\arguments{
\item{id}{(numeric/integer) an AphiaID. For \code{wm_children} it's
required and must be \code{length(id) == 1}, for \code{wm_children_} it's
optional and can be \code{length(id) >= 1}}

\item{marine_only}{(logical) marine only or not. default: \code{TRUE}}

\item{offset}{(integer) record to start at. default: 1}

\item{...}{named curl options. see \code{curl::curl_options}}

\item{name}{(character) one or more taxonomic names. optional}
}
\value{
A tibble/data.frame. when using underscore method, outputs from
each input are binded together, but can be split by \code{id} column
}
\description{
Get children for an AphiaID
}
\section{Singular vs. plural}{

Of the two sister functions, the one without the underscore is the original
function that wraps the relavant WoRMS API method - and only accepts
one thing (i.e., name or AphiaID) per request.

The sister function with the underscore at the end is the plural version,
accepting more than one input. Internally this function loops over
the non-underscore method, and labels output (whether it's a list or
data.frame rows) with the input names or IDs so that you can easily
parse output by your inputs.
}

\examples{
\dontrun{
wm_children(343613)
wm_children(id = 105706)
wm_children(id = 105706, FALSE)
wm_children(id = 105706, offset = 5)

# plural version, via id or name
wm_children_(id = c(105706, 343613))
wm_children_(name = c('Mesodesma', 'Leucophaeus'))
}
}
