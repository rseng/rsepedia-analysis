natserv
=======



[![cran checks](https://cranchecks.info/badges/worst/natserv)](https://cranchecks.info/pkgs/natserv)
[![R-check](https://github.com/ropensci/natserv/workflows/R-check/badge.svg)](https://github.com/ropensci/natserv/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/natserv/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/natserv)
[![cran version](https://www.r-pkg.org/badges/version/natserv)](https://cran.r-project.org/package=natserv)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/natserv)](https://github.com/metacran/cranlogs.app)


`natserv` NatureServe R client

NatureServe is a non-profit organization that provides wildlife conservation related data to various groups including the public.

* NatureServe site: https://services.natureserve.org
* NatureServe API docs: https://explorer.natureserve.org/api-docs/

All functions in this package are prefixed with `ns_` to prevent
collision with other pkgs.

You no longer need an API key.

See also the taxize book (https://taxize.dev/) for 
a manual on working with taxonomic data in R, including with NatureServe data.

## Installation

Stable version from CRAN


```r
install.packages("natserv")
```

Development version


```r
remotes::install_github("ropensci/natserv")
```


```r
library('natserv')
```

## Search


```r
ns_search_comb(text = "robin", page = 0, per_page = 5)
#> $results
#> # A tibble: 5 x 15
#>   recordType elementGlobalId uniqueId nsxUrl elcode scientificName
#>   <chr>                <int> <chr>    <chr>  <chr>  <chr>         
#> 1 SPECIES             100637 ELEMENT… /Taxo… ABPBJ… Copsychus sau…
#> 2 SPECIES             102323 ELEMENT… /Taxo… ABPBJ… Turdus grayi  
#> 3 SPECIES             102179 ELEMENT… /Taxo… ABPBJ… Turdus migrat…
#> 4 SPECIES             105536 ELEMENT… /Taxo… ABPBJ… Turdus migrat…
#> 5 SPECIES             105850 ELEMENT… /Taxo… ABPBJ… Turdus rufopa…
#> # … with 23 more variables: formattedScientificName <chr>,
#> #   primaryCommonName <chr>, primaryCommonNameLanguage <chr>,
#> #   roundedGRank <chr>, nations <list>, lastModified <chr>,
#> #   classificationStatus <chr>, speciesGlobal$usesaCode <lgl>,
#> #   $cosewicCode <lgl>, $saraCode <lgl>, $synonyms <list>,
#> #   $otherCommonNames <list>, $kingdom <chr>, $phylum <chr>, $taxclass <chr>,
#> #   $taxorder <chr>, $family <chr>, $genus <chr>, $taxonomicComments <chr>,
#> #   $informalTaxonomy <chr>, $infraspecies <lgl>, $completeDistribution <lgl>,
#> #   gRank <chr>
#> 
#> $resultsSummary
#>                            name value
#> 1                          page     0
#> 2                recordsPerPage     5
#> 3                    totalPages    26
#> 4                  totalResults   126
#> 5                 species_total   103
#> 6                         total    23
#> 7                       classes     0
#> 8                    subclasses     0
#> 9                    formations     0
#> 10                    divisions     0
#> 11                  macrogroups     1
#> 12                       groups     1
#> 13                    alliances     3
#> 14                 associations    18
#> 15 terrestrialEcologicalSystems     0
```

See the vignette (https://docs.ropensci.org/natserv/articles/natserv.html) for more examples.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/natserv/issues).
* License: MIT
* Get citation information for `natserv` in R doing `citation(package = 'natserv')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
natserv 1.0.0
=============

### BREAKING CHANGES

* All old functions are gone and have been replaced by new functions. NatureServe has released a new API and is ending support for the old one on 15 June 2020, so there's no reason to support the old API. See the docs (https://docs.ropensci.org/natserv/) for example usage (#21)


natserv 0.4.0
=============

### MINOR IMPROVEMENTS

* update vignette to use rmarkdown style, make sure vignette has a title (#20)


natserv 0.3.0
=============

### NEW FEATURES

* gains a vignette: "Introduction to `natserve` - An R package to wrap NatureServe's database API from rOpenSci" (#6) (#16) all work done by @mairindeith - thanks!

### MINOR IMPROVEMENTS

* add `vcr` for tests so that http requests are cached (#18)
* link to taxize book in readme and vignette (#15)
* improve failure behavior. the NatureServe API returns 200 on no result found, so we have to do gymnastics to catch that scenario (#9)

### BUG FIXES

* bug in `ns_data()` fixed related to subnations (#10)
* a bunch of checks added within `ns_data()`, `ns_images()`, and `ns_search()` to catch common bad inputs to the functions (#11)
* fixed bug in internal function `check_uid()` that was throwing warnings in `ns_data()` when passing more than one input (#14)


natserv 0.1.4
=============

### MINOR IMPROVEMENTS

* `natserv` now requires `crul` `>= 0.2.0`, which fixed URL encoding
to make our work in `natserv` easier.


natserv 0.1.0
=============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.0
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 reverse dependency. (Summary at <https://github.com/ropensci/natserv/blob/master/revdep/README.md>). There is a problem 
caused with the 1 dependency (taxize) because of function name changes in this version; a new release of taxize will follow shortly after this one.

---

This version improves changes all exported functions following release of a new API by NatureServe.

Thanks!
Scott Chamberlain
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

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

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the `natserv` project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.0 (2020-04-24) |
|os       |macOS Catalina 10.15.4       |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |US/Pacific                   |
|date     |2020-05-15                   |

# Dependencies

|package |old   |new      |Δ  |
|:-------|:-----|:--------|:--|
|natserv |0.4.0 |0.4.9.99 |*  |

# Revdeps

## New problems (1)

|package                      |version |error |warning |note   |
|:----------------------------|:-------|:-----|:-------|:------|
|[taxize](problems.md#taxize) |0.9.95  |      |        |__+1__ |

# taxize

<details>

* Version: 0.9.95
* Source code: https://github.com/cran/taxize
* URL: https://docs.ropensci.org/taxize (website), https://github.com/ropensci/taxize (devel), https://taxize.dev (user manual)
* BugReports: https://github.com/ropensci/taxize/issues
* Date/Publication: 2020-04-27 20:20:02 UTC
* Number of recursive dependencies: 99

Run `revdep_details(,"taxize")` for more info

</details>

## Newly broken

*   checking dependencies in R code ... NOTE
    ```
    Missing or unexported objects:
      ‘natserv::ns_data’ ‘natserv::ns_search’
    ```

*Wow, no problems at all. :)*natserv
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![cran checks](https://cranchecks.info/badges/worst/natserv)](https://cranchecks.info/pkgs/natserv)
[![R-check](https://github.com/ropensci/natserv/workflows/R-check/badge.svg)](https://github.com/ropensci/natserv/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/natserv/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/natserv)
[![cran version](https://www.r-pkg.org/badges/version/natserv)](https://cran.r-project.org/package=natserv)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/natserv)](https://github.com/metacran/cranlogs.app)


`natserv` NatureServe R client

NatureServe is a non-profit organization that provides wildlife conservation related data to various groups including the public.

* NatureServe site: https://services.natureserve.org
* NatureServe API docs: https://explorer.natureserve.org/api-docs/

All functions in this package are prefixed with `ns_` to prevent
collision with other pkgs.

You no longer need an API key.

See also the taxize book (https://taxize.dev/) for 
a manual on working with taxonomic data in R, including with NatureServe data.

## Installation

Stable version from CRAN

```{r eval=FALSE}
install.packages("natserv")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/natserv")
```

```{r}
library('natserv')
```

## Search

```{r}
ns_search_comb(text = "robin", page = 0, per_page = 5)
```

See the vignette (https://docs.ropensci.org/natserv/articles/natserv.html) for more examples.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/natserv/issues).
* License: MIT
* Get citation information for `natserv` in R doing `citation(package = 'natserv')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "natserv introduction"
author: "Scott Chamberlain"
date: "2020-05-15"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{natserv introduction}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



`natserv` is an R package that interacts with the API services of the non-profit organization NatureServe (https://services.natureserve.org/). If you want to read the full API documentation, you can find it at https://explorer.natureserve.org/api-docs/

See also the taxize book (https://taxize.dev/) for 
a manual on working with taxonomic data in R, including with NatureServe data.

This tutorial will walk you through installing `natserv` and using its functions. 

## A quick introduction to NatureServe

NatureServe is a non-profit organization that provides biodiversity data freely online. 
They maintain a database comprised of data from natural heritage programs and conservation data centers - this database includes information about the conservation status, taxonomy, geographic distribution, and life history information for over 70,000 species of plants, animals, and fungi in Canada and the United States.
You can find information about their data coverage (https://explorer.natureserve.org/AboutTheData/DataCoverage), and data sources (https://explorer.natureserve.org/AboutTheData/Sources) on their website. NatureServe also hosts data on ecological communities/systems and their conservation status.  

While small amounts of data can be easily collected using their online NatureServe explorer site (http://explorer.natureserve.org/), downloading species data this way would be incredibly slow. 
Thus `natserv` was born. 
This R package can access NatureServe's online API for rapid downloading of conservation data, allows for easy access to multiple species' data sets, and  loads the data directly into your R session. 

## Installing `natserv` from CRAN or GitHub

Stable version:


```r
install.packages("natserv")
```

Development version:


```r
remotes::install_github("ropensci/natserv")
```

After successful installation, load the package into the environment:


```r
library(natserv)
```

## `natserv` functions

All of `natserv`'s functions are prefixed with `ns_` to avoid confusion with other packages - here are the functions provided by `natserv`:


```r
cat(paste(" -", paste(sprintf("`%s`", sort(getNamespaceExports("natserv"))), collapse = "\n - ")))
#>  - `ns_altid`
#>  - `ns_ecohier`
#>  - `ns_export`
#>  - `ns_export_status`
#>  - `ns_id`
#>  - `ns_search_comb`
#>  - `ns_search_eco`
#>  - `ns_search_spp`
```

## Search

There's three functions for search:

- combined: supports searching for both species and ecosystems using search criteria which are applicable to both types of records
- species: supports searching for only species; extends the search criteria available through the Combined Search to include support for additional criteria which are only applicable to species records
- ecosystem: supports searching for only ecosystems; extends the search criteria available through the Combined Search to include support for additional criteria which are only applicable to ecosystems records

Combined search


```r
ns_search_comb(text = "robin", page = 0, per_page = 5)
#> $results
#> # A tibble: 5 x 14
#>   recordType elementGlobalId uniqueId nsxUrl elcode scientificName
#>   <chr>                <int> <chr>    <chr>  <chr>  <chr>         
#> 1 SPECIES             100637 ELEMENT… /Taxo… ABPBJ… Copsychus sau…
#> 2 SPECIES             102323 ELEMENT… /Taxo… ABPBJ… Turdus grayi  
#> 3 SPECIES             102179 ELEMENT… /Taxo… ABPBJ… Turdus migrat…
#> 4 SPECIES             105536 ELEMENT… /Taxo… ABPBJ… Turdus migrat…
#> 5 SPECIES             105850 ELEMENT… /Taxo… ABPBJ… Turdus rufopa…
#> # … with 22 more variables: formattedScientificName <chr>,
#> #   primaryCommonName <chr>, primaryCommonNameLanguage <chr>,
#> #   roundedGRank <chr>, nations <list>, lastModified <chr>,
#> #   speciesGlobal$usesaCode <lgl>, $cosewicCode <lgl>, $saraCode <lgl>,
#> #   $synonyms <list>, $otherCommonNames <list>, $kingdom <chr>, $phylum <chr>,
#> #   $taxclass <chr>, $taxorder <chr>, $family <chr>, $genus <chr>,
#> #   $taxonomicComments <chr>, $informalTaxonomy <chr>, $infraspecies <lgl>,
#> #   $completeDistribution <lgl>, gRank <chr>
#> 
#> $resultsSummary
#>                            name value
#> 1                          page     0
#> 2                recordsPerPage     5
#> 3                    totalPages    26
#> 4                  totalResults   126
#> 5                 species_total   103
#> 6                         total    23
#> 7                       classes     0
#> 8                    subclasses     0
#> 9                    formations     0
#> 10                    divisions     0
#> 11                  macrogroups     1
#> 12                       groups     1
#> 13                    alliances     3
#> 14                 associations    18
#> 15 terrestrialEcologicalSystems     0
```

Species search


```r
ns_search_spp(species_taxonomy = list(scientificTaxonomy = "Animalia", level = "kingdom"))
#> $results
#> # A tibble: 20 x 14
#>    recordType elementGlobalId uniqueId nsxUrl elcode scientificName
#>    <chr>                <int> <chr>    <chr>  <chr>  <chr>         
#>  1 SPECIES             828458 ELEMENT… /Taxo… AAABC… Acris blancha…
#>  2 SPECIES             828419 ELEMENT… /Taxo… AAABC… Acris crepita…
#>  3 SPECIES             103292 ELEMENT… /Taxo… AAABC… Acris gryllus 
#>  4 SPECIES             104324 ELEMENT… /Taxo… AAAAA… Ambystoma ann…
#>  5 SPECIES             100100 ELEMENT… /Taxo… AAAAA… Ambystoma bar…
#>  6 SPECIES             802300 ELEMENT… /Taxo… AAAAA… Ambystoma bis…
#>  7 SPECIES             104488 ELEMENT… /Taxo… AAAAA… Ambystoma cal…
#>  8 SPECIES             802301 ELEMENT… /Taxo… AAAAA… Ambystoma cin…
#>  9 SPECIES             103251 ELEMENT… /Taxo… AAAAA… Ambystoma gra…
#> 10 SPECIES             100401 ELEMENT… /Taxo… AAAAA… Ambystoma jef…
#> 11 SPECIES             102149 ELEMENT… /Taxo… AAAAA… Ambystoma lat…
#> 12 SPECIES            1143335 ELEMENT… /Taxo… AAAAA… Ambystoma lat…
#> 13 SPECIES            1143337 ELEMENT… /Taxo… AAAAA… Ambystoma lat…
#> 14 SPECIES             101261 ELEMENT… /Taxo… AAAAA… Ambystoma mab…
#> 15 SPECIES             106403 ELEMENT… /Taxo… AAAAA… Ambystoma mac…
#> 16 SPECIES             101632 ELEMENT… /Taxo… AAAAA… Ambystoma mac…
#> 17 SPECIES             103392 ELEMENT… /Taxo… AAAAA… Ambystoma mac…
#> 18 SPECIES             104239 ELEMENT… /Taxo… AAAAA… Ambystoma mac…
#> 19 SPECIES             100115 ELEMENT… /Taxo… AAAAA… Ambystoma mac…
#> 20 SPECIES             103039 ELEMENT… /Taxo… AAAAA… Ambystoma mac…
#> # … with 22 more variables: formattedScientificName <chr>,
#> #   primaryCommonName <chr>, primaryCommonNameLanguage <chr>,
#> #   roundedGRank <chr>, nations <list>, lastModified <chr>,
#> #   speciesGlobal$usesaCode <chr>, $cosewicCode <chr>, $saraCode <chr>,
#> #   $synonyms <list>, $otherCommonNames <list>, $kingdom <chr>, $phylum <chr>,
#> #   $taxclass <chr>, $taxorder <chr>, $family <chr>, $genus <chr>,
#> #   $taxonomicComments <chr>, $informalTaxonomy <chr>, $infraspecies <lgl>,
#> #   $completeDistribution <lgl>, gRank <chr>
#> 
#> $resultsSummary
#>             name value
#> 1           page     0
#> 2 recordsPerPage    20
#> 3     totalPages  2586
#> 4   totalResults 51705
#> 5  species_total 51705
```

Ecosystem search


```r
ns_search_eco(ecosystem_taxonomy = "M067")
#> $results
#> # A tibble: 20 x 14
#>    recordType elementGlobalId uniqueId nsxUrl elcode scientificName
#>    <chr>                <int> <chr>    <chr>  <chr>  <chr>         
#>  1 ECOSYSTEM           899525 ELEMENT… /Taxo… A3401  Eleocharis el…
#>  2 ECOSYSTEM           899121 ELEMENT… /Taxo… A1372  Fimbristylis …
#>  3 ECOSYSTEM           899126 ELEMENT… /Taxo… A1389  Spartina bake…
#>  4 ECOSYSTEM           899744 ELEMENT… /Taxo… A3692  Spartina pate…
#>  5 ECOSYSTEM           899523 ELEMENT… /Taxo… A3399  Typha dominge…
#>  6 ECOSYSTEM           899526 ELEMENT… /Taxo… A3402  Andropogon ca…
#>  7 ECOSYSTEM           899517 ELEMENT… /Taxo… A3393  Aristida palu…
#>  8 ECOSYSTEM           899132 ELEMENT… /Taxo… A1429  Eleocharis sp…
#>  9 ECOSYSTEM           899512 ELEMENT… /Taxo… A3388  Hypericum cha…
#> 10 ECOSYSTEM           899133 ELEMENT… /Taxo… A1430  Juncus milita…
#> 11 ECOSYSTEM           899122 ELEMENT… /Taxo… A1379  Panicum hemit…
#> 12 ECOSYSTEM           899124 ELEMENT… /Taxo… A1383  Rhynchospora …
#> 13 ECOSYSTEM           899518 ELEMENT… /Taxo… A3394  Rhynchospora …
#> 14 ECOSYSTEM          1146914 ELEMENT… /Taxo… A4404  Rhynchospora …
#> 15 ECOSYSTEM           899521 ELEMENT… /Taxo… A3397  Rhynchospora …
#> 16 ECOSYSTEM           899125 ELEMENT… /Taxo… A1384  Rhynchospora …
#> 17 ECOSYSTEM           899530 ELEMENT… /Taxo… A3406  Cladium maris…
#> 18 ECOSYSTEM           899120 ELEMENT… /Taxo… A1369  Cladium maris…
#> 19 ECOSYSTEM           899720 ELEMENT… /Taxo… A3668  Hottonia infl…
#> 20 ECOSYSTEM           900115 ELEMENT… /Taxo… A4065  Orontium aqua…
#> # … with 18 more variables: formattedScientificName <chr>,
#> #   primaryCommonName <chr>, primaryCommonNameLanguage <chr>,
#> #   roundedGRank <chr>, nations <list>, lastModified <chr>,
#> #   ecosystemGlobal$translatedScientificName <chr>, $taxclassCode <chr>,
#> #   $taxsubclassCode <chr>, $formationCode <chr>, $divisionCode <chr>,
#> #   $macrogroupKey <chr>, $taxgroupKey <chr>, $allianceKey <lgl>,
#> #   $ecosystemType <chr>, $classificationCode <chr>, $parentName <chr>,
#> #   gRank <chr>
#> 
#> $resultsSummary
#>                            name value
#> 1                          page     0
#> 2                recordsPerPage    20
#> 3                    totalPages    11
#> 4                  totalResults   210
#> 5                         total   210
#> 6                       classes     0
#> 7                    subclasses     0
#> 8                    formations     0
#> 9                     divisions     0
#> 10                  macrogroups     1
#> 11                       groups     6
#> 12                    alliances    33
#> 13                 associations   154
#> 14 terrestrialEcologicalSystems    16
```

## Get taxon by id

By UID


```r
w <- ns_id("ELEMENT_GLOBAL.2.154701")
str(w, max.level = 1)
#> List of 40
#>  $ elementGlobalId                        : int 154701
#>  $ circumscripConfidence                  : NULL
#>  $ classificationLevel                    :List of 4
#>  $ classificationStatus                   :List of 4
#>  $ iucn                                   :List of 5
#>  $ nameCategory                           :List of 6
#>  $ rankMethodUsed                         :List of 4
#>  $ formattedScientificName                : chr "<i>Hydrastis canadensis</i>"
#>  $ scientificName                         : chr "Hydrastis canadensis"
#>  $ scientificNameAuthor                   : chr "L."
#>  $ primaryCommonName                      : chr "Goldenseal"
#>  $ relatedItisNames                       : chr "<i>Hydrastis canadensis</i> L. (TSN 18781)"
#>  $ uniqueId                               : chr "ELEMENT_GLOBAL.2.154701"
#>  $ elcode                                 : chr "PDRAN0F010"
#>  $ conceptRefFullCitation                 : chr "Kartesz, J.T. 1994. A synonymized checklist of the vascular flora of the United States, Canada, and Greenland. "| __truncated__
#>  $ conceptName                            : chr "<i>Hydrastis canadensis</i>"
#>  $ taxonomicComments                      : chr "<i>Hydrastis canadensis</i> occurs in eastern North America and is a monotypic genus. In the most current taxon"| __truncated__
#>  $ roundedGRank                           : chr "G3"
#>  $ conservationStatusFactorsEditionDate   : chr "2013-04-29"
#>  $ conservationStatusFactorsEditionAuthors: chr "Oliver, L."
#>  $ primaryCommonNameLanguage              : chr "EN"
#>  $ recordType                             : chr "SPECIES"
#>  $ elementNationals                       :'data.frame':	2 obs. of  7 variables:
#>  $ lastModified                           : chr "2020-05-14T04:31:58.480462Z"
#>  $ lastPublished                          : chr "2020-05-14T02:14:14.569584Z"
#>  $ nsxUrl                                 : chr "/Taxon/ELEMENT_GLOBAL.2.154701/Hydrastis_canadensis"
#>  $ grank                                  : chr "G3G4"
#>  $ grankReviewDate                        : chr "2012-11-30"
#>  $ grankChangeDate                        : chr "2012-11-30"
#>  $ grankReasons                           : chr "Goldenseal, <i>Hydrastis canadensis, </i>an herbaceous understory species of the eastern deciduous forest, with"| __truncated__
#>  $ rankInfo                               :List of 27
#>  $ animalCharacteristics                  : NULL
#>  $ occurrenceDelineations                 :'data.frame':	1 obs. of  17 variables:
#>  $ plantCharacteristics                   :List of 7
#>  $ elementManagement                      :List of 14
#>  $ occurrenceViabilities                  :'data.frame':	1 obs. of  12 variables:
#>  $ references                             :'data.frame':	97 obs. of  6 variables:
#>  $ otherCommonNames                       :'data.frame':	6 obs. of  3 variables:
#>  $ speciesGlobal                          :List of 29
#>  $ speciesCharacteristics                 :List of 13
```

By alternate id


```r
x <- ns_altid(uid = "ELEMENT_GLOBAL.2.154701")
str(x, max.level = 1)
#> List of 40
#>  $ elementGlobalId                        : int 154701
#>  $ circumscripConfidence                  : NULL
#>  $ classificationLevel                    :List of 4
#>  $ classificationStatus                   :List of 4
#>  $ iucn                                   :List of 5
#>  $ nameCategory                           :List of 6
#>  $ rankMethodUsed                         :List of 4
#>  $ formattedScientificName                : chr "<i>Hydrastis canadensis</i>"
#>  $ scientificName                         : chr "Hydrastis canadensis"
#>  $ scientificNameAuthor                   : chr "L."
#>  $ primaryCommonName                      : chr "Goldenseal"
#>  $ relatedItisNames                       : chr "<i>Hydrastis canadensis</i> L. (TSN 18781)"
#>  $ uniqueId                               : chr "ELEMENT_GLOBAL.2.154701"
#>  $ elcode                                 : chr "PDRAN0F010"
#>  $ conceptRefFullCitation                 : chr "Kartesz, J.T. 1994. A synonymized checklist of the vascular flora of the United States, Canada, and Greenland. "| __truncated__
#>  $ conceptName                            : chr "<i>Hydrastis canadensis</i>"
#>  $ taxonomicComments                      : chr "<i>Hydrastis canadensis</i> occurs in eastern North America and is a monotypic genus. In the most current taxon"| __truncated__
#>  $ roundedGRank                           : chr "G3"
#>  $ conservationStatusFactorsEditionDate   : chr "2013-04-29"
#>  $ conservationStatusFactorsEditionAuthors: chr "Oliver, L."
#>  $ primaryCommonNameLanguage              : chr "EN"
#>  $ recordType                             : chr "SPECIES"
#>  $ elementNationals                       :'data.frame':	2 obs. of  7 variables:
#>  $ lastModified                           : chr "2020-05-14T04:31:58.480462Z"
#>  $ lastPublished                          : chr "2020-05-14T02:14:14.569584Z"
#>  $ nsxUrl                                 : chr "/Taxon/ELEMENT_GLOBAL.2.154701/Hydrastis_canadensis"
#>  $ grank                                  : chr "G3G4"
#>  $ grankReviewDate                        : chr "2012-11-30"
#>  $ grankChangeDate                        : chr "2012-11-30"
#>  $ grankReasons                           : chr "Goldenseal, <i>Hydrastis canadensis, </i>an herbaceous understory species of the eastern deciduous forest, with"| __truncated__
#>  $ rankInfo                               :List of 27
#>  $ animalCharacteristics                  : NULL
#>  $ occurrenceDelineations                 :'data.frame':	1 obs. of  17 variables:
#>  $ plantCharacteristics                   :List of 7
#>  $ elementManagement                      :List of 14
#>  $ occurrenceViabilities                  :'data.frame':	1 obs. of  12 variables:
#>  $ references                             :'data.frame':	97 obs. of  6 variables:
#>  $ otherCommonNames                       :'data.frame':	6 obs. of  3 variables:
#>  $ speciesGlobal                          :List of 29
#>  $ speciesCharacteristics                 :List of 13
```

## Get a summary of the upper level hierarchy for an Ecosystem record


```r
ns_ecohier("ELEMENT_GLOBAL.2.683060")
#>                  uniqueId
#> 1 ELEMENT_GLOBAL.2.860217
#> 2 ELEMENT_GLOBAL.2.860227
#> 3 ELEMENT_GLOBAL.2.860241
#> 4 ELEMENT_GLOBAL.2.860284
#> 5 ELEMENT_GLOBAL.2.838501
#> 6 ELEMENT_GLOBAL.2.833279
#> 7 ELEMENT_GLOBAL.2.899395
#>                                                            name
#> 1                                             Forest & Woodland
#> 2                          Temperate & Boreal Forest & Woodland
#> 3                              Cool Temperate Forest & Woodland
#> 4                      Eastern North American Forest & Woodland
#> 5         Southern & South-Central Oak - Pine Forest & Woodland
#> 6 South-Central Interior Shortleaf Pine - Oak Forest & Woodland
#> 7                  Ozark-Ouachita Shortleaf Pine - Oak Woodland
#>                                                                                                                 nsxUrl
#> 1                                                     /Taxon/ELEMENT_GLOBAL.2.860217/Mesomorphic_Tree_Vegetation_Class
#> 2                                             /Taxon/ELEMENT_GLOBAL.2.860227/Temperate_Boreal_Forest_Woodland_Subclass
#> 3                                              /Taxon/ELEMENT_GLOBAL.2.860241/Cool_Temperate_Forest_Woodland_Formation
#> 4           /Taxon/ELEMENT_GLOBAL.2.860284/Acer_saccharum_-_Fagus_grandifolia_-_Quercus_rubra_Forest_Woodland_Division
#> 5            /Taxon/ELEMENT_GLOBAL.2.838501/Quercus_alba_-_Quercus_falcata_-_Pinus_echinata_Forest_Woodland_Macrogroup
#> 6             /Taxon/ELEMENT_GLOBAL.2.833279/Pinus_echinata_-_Quercus_falcata_-_Quercus_stellata_Forest_Woodland_Group
#> 7 /Taxon/ELEMENT_GLOBAL.2.899395/Pinus_echinata_-_Quercus_stellata_-_Quercus_velutina_Ozark-Ouachita_Woodland_Alliance
#>   ecosystemType classificationCode
#> 1         CLASS                  1
#> 2      SUBCLASS                1.B
#> 3     FORMATION              1.B.2
#> 4      DIVISION           1.B.2.Na
#> 5    MACROGROUP               M016
#> 6         GROUP               G012
#> 7      ALLIANCE              A3271
```

## Search exports

`ns_export()` uses the same search interface as the `ns_search*` functions, but instead of downloading data immediately, `ns_export()` creates a "download job", which eventually provides a compressed JSON file that you can download.


```r
x <- ns_export(text = "robin")
x
#> [1] "7a107bea-b98d-4b5a-87b3-456ea2194f07"
```

You can pass the output of `ns_export()` to `ns_export_status()` to get the status of the job


```r
res <- ns_export_status(x)
#> $state
#> [1] "Finished"
#> 
#> $data
#> $data$success
#> [1] TRUE
#> $percentComplete
#> [1] 100
#> $successful
#> [1] TRUE
#> $error
#> ...
```

When state equals "Finished", you can read the data into R, e.g, with `jsonlite`:


```r
res$data$url
#> [1] "https://explorer-downloads.natureserve.org/shortTerm/explorer/taxaSearchExports/2020-05-15/7a107bea-b98d-4b5a-87b3-456ea2194f07.json"
tibble::as_tibble(jsonlite::fromJSON(res$data$url))
#> # A tibble: 126 x 14
#>    elementGlobalId uniqueId nsxUrl elcode scientificName formattedScient… primaryCommonNa… primaryCommonNa… roundedGRank
#>              <int> <chr>    <chr>  <chr>  <chr>          <chr>            <chr>            <chr>            <chr>
#>  1          100637 ELEMENT… /Taxo… ABPBJ… Copsychus sau… <i>Copsychus sa… Oriental Magpie… EN               G5
#>  2          102323 ELEMENT… /Taxo… ABPBJ… Turdus grayi   <i>Turdus grayi… Clay-colored Th… EN               G5
#>  3          102179 ELEMENT… /Taxo… ABPBJ… Turdus migrat… <i>Turdus migra… American Robin   EN               G5
#>  4          105536 ELEMENT… /Taxo… ABPBJ… Turdus migrat… <i>Turdus migra… Western America… EN               TU
#>  5          105850 ELEMENT… /Taxo… ABPBJ… Turdus rufopa… <i>Turdus rufop… Rufous-backed R… EN               G5
#>  6          100589 ELEMENT… /Taxo… AFC4B… Peristedion g… <i>Peristedion … Slender Searobin EN               GNR
#>  7          105826 ELEMENT… /Taxo… AFC4B… Prionotus ala… <i>Prionotus al… Spiny Searobin   EN               GNR
#>  8          101394 ELEMENT… /Taxo… AFC4B… Prionotus car… <i>Prionotus ca… Northern Searob… EN               G5
#>  9          100276 ELEMENT… /Taxo… AFC4B… Prionotus evo… <i>Prionotus ev… Striped Searobin EN               G5
#> 10          103595 ELEMENT… /Taxo… AFC4B… Prionotus lon… <i>Prionotus lo… Bigeye Searobin  EN               G5
#> # … with 116 more rows, and 29 more variables: nations <list>, lastModified <chr>, speciesGlobal$usesaCode <chr>,
#> #   $cosewicCode <chr>, $saraCode <chr>, $synonyms <list>, $otherCommonNames <list>, $kingdom <chr>, $phylum <chr>,
#> #   $taxclass <chr>, $taxorder <chr>, $family <chr>, $genus <chr>, $taxonomicComments <chr>, $informalTaxonomy <chr>,
#> #   $infraspecies <lgl>, $completeDistribution <lgl>, gRank <chr>, ecosystemGlobal$translatedScientificName <chr>,
#> #   $taxclassCode <chr>, $taxsubclassCode <chr>, $formationCode <chr>, $divisionCode <chr>, $macrogroupKey <chr>,
#> #   $taxgroupKey <chr>, $allianceKey <chr>, $ecosystemType <chr>, $classificationCode <chr>, $parentName <chr>
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ns_search_spp.R
\name{ns_search_spp}
\alias{ns_search_spp}
\title{Species search}
\usage{
ns_search_spp(
  text = NULL,
  text_adv = NULL,
  status = NULL,
  location = NULL,
  species_taxonomy = NULL,
  record_subtype = NULL,
  modified_since = NULL,
  page = NULL,
  per_page = NULL,
  ...
)
}
\arguments{
\item{text}{(character) basic text search, equiavalent to \code{text_adv}
with \code{matchAgainst="allNames"} and \code{operator="similarTo"}}

\item{text_adv}{(list) advanced search, must specify the following three
elements: \code{searchToken}, \code{matchAgainst}, and \code{operator}. see
https://explorer.natureserve.org/api-docs/#_advanced_text_search_parameter}

\item{status}{(character) conservation status, one of G1, G2, G3, G4,
G5, GH, GX, GNR, GNA, GU. case insensitive}

\item{location}{(list) location, country and sub-country. specify either
\code{nation} OR \code{nation} and \code{subnation}. each expects a two-letter ISO code}

\item{species_taxonomy}{(list) species taxonomy. either a list with
\code{level} and \code{scientificTaxonomy} (a scientific name), or with just
\code{informalTaxonomy} (a vernacular name). possible \code{level} values:
"kingdom", "phylum", "class", "order", "family", "genus"}

\item{record_subtype}{(character) limit results by record sub-type, one of:
"class", "subclass", "formation", "division", "macrogroup", "group",
"alliance", "association", "terrestrial_ecological_system"}

\item{modified_since}{(character) search for records modified since a
given time. value must be a date and time with a UTC offset in ISO 8601
format. optional}

\item{page}{(integer) Zero-indexed page number; default: 0. optional}

\item{per_page}{(integer) Records per page; default: 20. optional}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\description{
Species search
}
\examples{
\dontrun{
ns_search_spp(text = "robin")
ns_search_spp(text_adv = list(searchToken = "bird", 
  matchAgainst = "allNames", operator="similarTo"))
ns_search_spp(status = "G1")
ns_search_spp(location = list(nation = "US"))
ns_search_spp(location = list(nation = "US", subnation = "VA"))
ns_search_spp(species_taxonomy = list(scientificTaxonomy = "Animalia", level = "kingdom"))
ns_search_spp(species_taxonomy = list(informalTaxonomy = "birds"))
ns_search_spp(record_subtype = "macrogroup")
ns_search_spp(modified_since = "2020-04-30T00:00:00+0000")
ns_search_spp(page = 0, per_page = 2)
}
}
\references{
https://explorer.natureserve.org/api-docs/
}
\seealso{
Other search: 
\code{\link{ns_search_comb}()},
\code{\link{ns_search_eco}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ns_id.R
\name{ns_id}
\alias{ns_id}
\title{Get taxon by uid}
\usage{
ns_id(uid, ...)
}
\arguments{
\item{uid}{(character) A NatureServe taxon id (The taxon’s Element Global
UID). required.}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
A list with lots of elements
}
\description{
Get taxon by uid
}
\details{
see https://explorer.natureserve.org/api-docs/#_taxon_data_model
for details on the response data
}
\examples{
\dontrun{
ns_id("ELEMENT_GLOBAL.2.154701")
}
}
\references{
https://explorer.natureserve.org/api-docs/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ns_search_eco.R
\name{ns_search_eco}
\alias{ns_search_eco}
\title{Ecosystem search}
\usage{
ns_search_eco(
  text = NULL,
  text_adv = NULL,
  status = NULL,
  location = NULL,
  ecosystem_taxonomy = NULL,
  record_subtype = NULL,
  modified_since = NULL,
  page = NULL,
  per_page = NULL,
  ...
)
}
\arguments{
\item{text}{(character) basic text search, equiavalent to \code{text_adv}
with \code{matchAgainst="allNames"} and \code{operator="similarTo"}}

\item{text_adv}{(list) advanced search, must specify the following three
elements: \code{searchToken}, \code{matchAgainst}, and \code{operator}. see
https://explorer.natureserve.org/api-docs/#_advanced_text_search_parameter}

\item{status}{(character) conservation status, one of G1, G2, G3, G4,
G5, GH, GX, GNR, GNA, GU. case insensitive}

\item{location}{(list) location, country and sub-country. specify either
\code{nation} OR \code{nation} and \code{subnation}. each expects a two-letter ISO code}

\item{ecosystem_taxonomy}{(character) the classification code of the
higher level (ancestor) ecosystem. E.g.'s: "1" (Class code),
"1.B" (Subclass code), "1.B.2" (Formation code), "1.B.2.Nd" (Division code),
"M886" (Macrogroup key), "G206" (Group key), "A3328" (Alliance Key)}

\item{record_subtype}{(character) limit results by record sub-type, one of:
"class", "subclass", "formation", "division", "macrogroup", "group",
"alliance", "association", "terrestrial_ecological_system"}

\item{modified_since}{(character) search for records modified since a
given time. value must be a date and time with a UTC offset in ISO 8601
format. optional}

\item{page}{(integer) Zero-indexed page number; default: 0. optional}

\item{per_page}{(integer) Records per page; default: 20. optional}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\description{
Ecosystem search
}
\examples{
\dontrun{
ns_search_eco(text = "robin")
ns_search_eco(text_adv = list(searchToken = "bird",
  matchAgainst = "allNames", operator="similarTo"))
ns_search_eco(status = "G1")
ns_search_eco(location = list(nation = "US"))
ns_search_eco(location = list(nation = "US", subnation = "VA"))
ns_search_eco(ecosystem_taxonomy = "M067")
ns_search_eco(record_subtype = "macrogroup")
ns_search_eco(modified_since = "2020-04-30T00:00:00+0000")
ns_search_eco(page = 0, per_page = 2)
}
}
\references{
https://explorer.natureserve.org/api-docs/
}
\seealso{
Other search: 
\code{\link{ns_search_comb}()},
\code{\link{ns_search_spp}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/natserv-package.R
\docType{package}
\name{natserv-package}
\alias{natserv-package}
\alias{natserv}
\title{natserv}
\description{
Interface to NatureServe https://www.natureserve.org/
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ns_altid.R
\name{ns_altid}
\alias{ns_altid}
\title{Get taxon by uid, id, or elCode}
\usage{
ns_altid(uid = NULL, id = NULL, el_code = NULL, ...)
}
\arguments{
\item{uid}{(character) A NatureServe taxon id (The taxon’s Element Global
UID)}

\item{id}{The primary key value (ELEMENT_GLOBAL_ID) of the record within
Central Biotics}

\item{el_code}{The Biotics Element Code (ELCODE_BCD) of the record}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
A list with lots of elements
}
\description{
Get taxon by uid, id, or elCode
}
\details{
see https://explorer.natureserve.org/api-docs/#_taxon_data_model
for details on the response data
}
\examples{
\dontrun{
ns_altid(uid = "ELEMENT_GLOBAL.2.154701")
ns_altid(id = "154701")
ns_altid(el_code = "PDRAN0F010")
}
}
\references{
https://explorer.natureserve.org/api-docs/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ns_export.R
\name{ns_export}
\alias{ns_export}
\alias{ns_export_status}
\title{Search exports}
\usage{
ns_export(
  text = NULL,
  text_adv = NULL,
  status = NULL,
  location = NULL,
  record_type = NULL,
  record_subtype = NULL,
  modified_since = NULL,
  format = "json",
  lang = "en",
  ...
)

ns_export_status(id, ...)
}
\arguments{
\item{text}{(character) basic text search, equiavalent to \code{text_adv}
with \code{matchAgainst="allNames"} and \code{operator="similarTo"}}

\item{text_adv}{(list) advanced search, must specify the following three
elements: \code{searchToken}, \code{matchAgainst}, and \code{operator}. see
https://explorer.natureserve.org/api-docs/#_advanced_text_search_parameter}

\item{status}{(character) conservation status, one of G1, G2, G3, G4,
G5, GH, GX, GNR, GNA, GU. case insensitive}

\item{location}{(list) location, country and sub-country. specify either
\code{nation} OR \code{nation} and \code{subnation}. each expects a two-letter ISO code}

\item{record_type}{(character) limit results by record type, one of
"species" or "ecosystem"}

\item{record_subtype}{(character) limit results by record sub-type, one of:
"class", "subclass", "formation", "division", "macrogroup", "group",
"alliance", "association", "terrestrial_ecological_system"}

\item{modified_since}{(character) search for records modified since a
given time. value must be a date and time with a UTC offset in ISO 8601
format. optional}

\item{format}{(character) output format, one of "json" or "xlsx"}

\item{lang}{(character) language, one of "en", "es", or "fr"}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}

\item{id}{(character) a job id, from output of \code{ns_export()}}
}
\value{
\code{ns_export()} returns a single character string (a job id)
\code{ns_export_status()} returns a list of metadata concerning the
status of the export
}
\description{
Search exports
}
\examples{
\dontrun{
x <- ns_export(text = "robin")
res <- ns_export_status(x)
str(res)
res$state
res$data$errorMessage
res$data$url

w <- ns_export(text_adv = list(searchToken = "western",
  matchAgainst="allScientificNames", operator="startsWith"))
m <- ns_export_status(w)
head(jsonlite::fromJSON(m$data$url))
}
}
\references{
https://explorer.natureserve.org/api-docs/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ns_search_comb.R
\name{ns_search_comb}
\alias{ns_search_comb}
\title{Combined search}
\usage{
ns_search_comb(
  text = NULL,
  text_adv = NULL,
  status = NULL,
  location = NULL,
  record_type = NULL,
  record_subtype = NULL,
  modified_since = NULL,
  page = NULL,
  per_page = NULL,
  ...
)
}
\arguments{
\item{text}{(character) basic text search, equiavalent to \code{text_adv}
with \code{matchAgainst="allNames"} and \code{operator="similarTo"}}

\item{text_adv}{(list) advanced search, must specify the following three
elements: \code{searchToken}, \code{matchAgainst}, and \code{operator}. see
https://explorer.natureserve.org/api-docs/#_advanced_text_search_parameter}

\item{status}{(character) conservation status, one of G1, G2, G3, G4,
G5, GH, GX, GNR, GNA, GU. case insensitive}

\item{location}{(list) location, country and sub-country. specify either
\code{nation} OR \code{nation} and \code{subnation}. each expects a two-letter ISO code}

\item{record_type}{(character) limit results by record type, one of
"species" or "ecosystem"}

\item{record_subtype}{(character) limit results by record sub-type, one of:
"class", "subclass", "formation", "division", "macrogroup", "group",
"alliance", "association", "terrestrial_ecological_system"}

\item{modified_since}{(character) search for records modified since a
given time. value must be a date and time with a UTC offset in ISO 8601
format. optional}

\item{page}{(integer) Zero-indexed page number; default: 0. optional}

\item{per_page}{(integer) Records per page; default: 20. optional}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\description{
Combined search
}
\examples{
\dontrun{
ns_search_comb(text = "robin")
ns_search_comb(text_adv = list(searchToken = "western",
  matchAgainst="allScientificNames", operator="startsWith"))
ns_search_comb(status = "G1")
ns_search_comb(location = list(nation = "US"))
ns_search_comb(location = list(nation = "US", subnation = "VA"))
ns_search_comb(record_type = "species")
ns_search_comb(record_subtype = "macrogroup")
ns_search_comb(modified_since = "2020-04-30T00:00:00+0000")
ns_search_comb(page = 0, per_page = 2)
}
}
\references{
https://explorer.natureserve.org/api-docs/
}
\seealso{
Other search: 
\code{\link{ns_search_eco}()},
\code{\link{ns_search_spp}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/natserv-package.R
\docType{data}
\name{nat_states}
\alias{nat_states}
\title{A data.frame with 49 rows and 2 columns}
\description{
\itemize{
\item state (character) state 2 letter abbreviation
\item state_name (character) state full name
}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ns_ecohier.R
\name{ns_ecohier}
\alias{ns_ecohier}
\title{Get a summary of the upper level hierarchy for an Ecosystem record}
\usage{
ns_ecohier(uid, ...)
}
\arguments{
\item{uid}{(character) A NatureServe taxon id (The taxon’s Element Global
UID)}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
A list with lots of elements
}
\description{
Get a summary of the upper level hierarchy for an Ecosystem record
}
\details{
see https://explorer.natureserve.org/api-docs/#_taxon_data_model
for details on the response data
}
\examples{
\dontrun{
ns_ecohier("ELEMENT_GLOBAL.2.683060")
}
}
\references{
https://explorer.natureserve.org/api-docs/
}
