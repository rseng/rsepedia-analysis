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

*Wow, no problems at all. :)*