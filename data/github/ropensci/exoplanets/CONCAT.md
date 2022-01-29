
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exoplanets <img src="man/figures/logo.png" align="right" height=150/>

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/exoplanets/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/exoplanets/actions/workflows/check-pak.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/exoplanets/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/exoplanets?branch=master)
[![Peer
review](https://badges.ropensci.org/309_status.svg)](https://github.com/ropensci/software-review/issues/309)
[![CRAN
status](https://www.r-pkg.org/badges/version/exoplanets)](https://CRAN.R-project.org/package=exoplanets)
[![CRAN\_Download\_Badge](https://cranlogs.r-pkg.org/badges/exoplanets)](https://cran.r-project.org/package=exoplanets)
<!-- badges: end -->

The goal of exoplanets is to provide access to [NASA’s Exoplanet Archive
TAP
Service](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html).
The functionality of this package is minimal and is simply an R
interface to access exoplanet data.

<img src="man/figures/README-unnamed-chunk-2-1.png" title="Exoplanets color coded by discovery method" alt="Exoplanets color coded by discovery method" width="100%" />

## Installation

Install the released version of `exoplanets` from CRAN:

``` r
install.packages("exoplanets")
```

Or you can install from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/exoplanets")
```

## Example

This is a basic example which shows you how to access data from the
[k2names](https://exoplanetarchive.ipac.caltech.edu/docs/API_k2names_columns.html)
table:

``` r
library(exoplanets)

options(
  exoplanets.progress = FALSE, # hide progress
  readr.show_types = FALSE     # hide col spec, requires readr 2.0.0 >=
)

exoplanets("k2names")
#> • https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+k2names&format=csv
#> # A tibble: 449 x 3
#>    epic_id        k2_name  pl_name     
#>    <chr>          <chr>    <chr>       
#>  1 EPIC 246199087 K2-112 f TRAPPIST-1 f
#>  2 EPIC 246199087 K2-112 h TRAPPIST-1 h
#>  3 EPIC 211331236 K2-117 c K2-117 c    
#>  4 EPIC 212398486 K2-125 b K2-125 b    
#>  5 EPIC 217941732 K2-130 b K2-130 b    
#>  6 EPIC 228754001 K2-132 b K2-132 b    
#>  7 EPIC 247887989 K2-133 d K2-133 d    
#>  8 EPIC 247589423 K2-136 b K2-136 b    
#>  9 EPIC 247589423 K2-136 d K2-136 d    
#> 10 EPIC 201912552 K2-18 c  K2-18 c     
#> # … with 439 more rows
```

If you wish, you can select only the columns you need:

``` r
exoplanets("ps", c("pl_name", "hostname"))
#> • https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,hostname+from+ps&format=csv
#> # A tibble: 29,683 x 2
#>    pl_name      hostname  
#>    <chr>        <chr>     
#>  1 Kepler-11 c  Kepler-11 
#>  2 Kepler-11 f  Kepler-11 
#>  3 HAT-P-1 b    HAT-P-1   
#>  4 OGLE-TR-10 b OGLE-TR-10
#>  5 TrES-2 b     TrES-2    
#>  6 WASP-3 b     WASP-3    
#>  7 HD 210702 b  HD 210702 
#>  8 BD-08 2823 b BD-08 2823
#>  9 BD-08 2823 c BD-08 2823
#> 10 HAT-P-30 b   HAT-P-30  
#> # … with 29,673 more rows
```

You can also specify the number of rows returned using `limit`:

``` r
exoplanets("keplernames", limit = 5)
#> • https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+keplernames+top+5&format=csv
#> # A tibble: 5 x 4
#>     kepid koi_name  kepler_name   pl_name      
#>     <dbl> <chr>     <chr>         <chr>        
#> 1 7515212 K00679.02 Kepler-212 b  Kepler-212 b 
#> 2 8210018 K02762.01 Kepler-1341 b Kepler-1341 b
#> 3 9008737 K02768.01 Kepler-404 b  Kepler-404 b 
#> 4 4833421 K00232.05 Kepler-122 f  Kepler-122 f 
#> 5 9963524 K00720.02 Kepler-221 d  Kepler-221 d
```

Information on the tables and columns available can be found with:

``` r
tableinfo
#> # A tibble: 546 x 13
#>    table database_column_… table_label description  in_ps_table in_ps_comp_pars…
#>    <chr> <chr>             <chr>       <chr>        <lgl>       <lgl>           
#>  1 ps    default_flag      Default Pa… Boolean fla… TRUE        FALSE           
#>  2 ps    soltype           Solution T… Disposition… TRUE        FALSE           
#>  3 ps    pl_controv_flag   Controvers… Flag indica… TRUE        TRUE            
#>  4 ps    pl_name           Planet Name Planet name… TRUE        TRUE            
#>  5 ps    hostname          Host Name   Stellar nam… TRUE        TRUE            
#>  6 ps    pl_letter         Planet Let… Letter assi… TRUE        TRUE            
#>  7 ps    hd_name           HD ID       Name of the… TRUE        TRUE            
#>  8 ps    hip_name          HIP ID      Name of the… TRUE        TRUE            
#>  9 ps    tic_id            TIC ID      Name of the… TRUE        TRUE            
#> 10 ps    gaia_id           GAIA ID     Name of the… TRUE        TRUE            
#> # … with 536 more rows, and 7 more variables:
#> #   uncertainties_column_positive_negative <chr>, limit_column <chr>,
#> #   default <lgl>, notes <chr>, displayed_string_name <chr>, flag_column <lgl>,
#> #   number_of_measurements <lgl>
```

## Capabilities

At one time, this package used the *Exoplanet Archive Application
Programming Interface (API)*. Since then, a handful of tables have been
transitioned to the *Table Access Protocol (TAP) service*. More tables
will be transitioned to TAP and as such, this package only supports
queries from TAP. For more information, you can read
<https://exoplanetarchive.ipac.caltech.edu/docs/exonews_archive.html#29April2021>.

## Contributing

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# exoplanets (development version)

# exoplanets 0.2.2

* `tableinfo` has been updated to include additional tables.

# exoplanets 0.2.1

## Breaking changes

* The `progress` parameter in `exoplanets` is no longer available. It has been replaced as an option that can be set with `options`.

## Minor improvements

* Added `forget_exoplanets` to clear the `exoplanets` cache.
* Added package level documentation with `usethis::use_package_doc`.
* Added `limit` parameter to `exoplanets`.
* Added `quiet` option to suppress progress and query message.
* Cleaned up `tableinfo` to remove trailing spaces and other small improvements.

# exoplanets 0.2.0

* Complete rewrite to work with the new TAP service: https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html
* No longer supporting the older API: https://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html
* Use of `httr` for progress indicators and better request handling.
* Added `tableinfo` dataset for table and column information.
* Removed all `exo_` functions, added `exoplanets` as a single function for extracting data.
* Added memoization to `exoplanets`.
* Using `httptest` to avoid potential CRAN errors when checking the package.
* Precompiled vignettes to avoid potential CRAN errors when checking the package.
* Transferred ownership to ropensci organization.

# exoplanets 0.1.0

* Added a `NEWS.md` file to track changes to the package.
## Test environments

* macOS-latest, 'release', github actions
* windows-latest, 'release', github actions
* windows-latest, '3.6', github actions
* ubuntu-18.04, 'devel', github actions
* ubuntu-18.04, 'release', github actions
* ubuntu-18.04, 'oldrel', github actions
* ubuntu-18.04, '3.5', github actions
* ubuntu-18.04, '3.4', github actions
* ubuntu-18.04, '3.3', github actions

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release. This release updates the `tableinfo` dataset.
# parse_url works

    argument "format" is missing, with no default

---

    argument "format" is missing, with no default

---

    argument "format" is missing, with no default

---

    Expected formats are:
    * csv
    * tsv
    * json

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
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the exoplanets project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
