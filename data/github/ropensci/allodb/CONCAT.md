
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://i.imgur.com/39pvr4n.png" align="left" height=44 /> allodb: An R package for biomass estimation at extratropical forest plots

<!-- badges: start -->

[![peer-review](https://badges.ropensci.org/436_status.svg)](https://github.com/ropensci/software-review/issues/436)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/allodb/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/allodb?branch=master)
[![R-CMD-check](https://github.com/ropensci/allodb/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/allodb/actions)
<!-- badges: end -->

## Introduction

Allometric equations for calculation of tree aboveground biomass (AGB)
form the basis for estimates of forest carbon storage and exchange with
the atmosphere. While standard models exist to calculate forest biomass
across the tropics, we lack a standardized tool for computing AGB across
the global extratropics.

*allodb* was conceived as a framework to standardize and simplify the
biomass estimation process across globally distributed extratropical
forests (mainly temperate and boreal forests). With *allodb* we aimed
to: a) compile relevant published and unpublished allometries, focusing
on AGB but structured to handle other variables (e.g., height); b)
objectively select and integrate appropriate available equations across
the full range of tree sizes; and c) serve as a platform for future
updates and expansion to other research sites.

The *allodb* package contains a dataset of systematically selected
published allometric equations. This dataset was built based on 701
woody species identified at 24 large [ForestGEO forest dynamic
plots](https://forestgeo.si.edu/) representing all major extratropical
forest types. A total of 570 parsed allometric equations to estimate
individual tree biomass were retrieved, checked, and combined using a
weighting function designed to ensure optimal equation selection over
the full tree size range with smooth transitions across equations. The
equation dataset used can be customized with built-in functions that
subset the original dataset and add new equations.

The package provides functions to estimate tree biomass based on
user-provided census data (tree diameter, taxonomic identification, and
plot coordinates). New allometric equations are calibrated for each
species and location by resampling the original equations; equations
with a larger sample size and/or higher taxonomic and climatic
similarity with the species and location in question are given a higher
weight in this process.

## Installation

Install the development version of *allodb* from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/allodb")
```

## Examples

Prior to calculating tree biomass using *allodb*, users need to provide
a table (i.e. dataframe) with DBH (cm), parsed species Latin names, and
site(s) coordinates. In the following examples we use data from the
Smithsonian Conservation Biology Institute, USA (SCBI) ForestGEO
dynamics plot (trees from 1 hectare surveyed in 2008). Full tree census
data can be requested through the [ForestGEO
portal](https://forestgeo.si.edu/explore-data).

``` r
library(allodb)
data(scbi_stem1)
```

The biomass of all trees in one (or several) censuses can be estimated
using the `get_biomass` function.

``` r
scbi_stem1$agb <-
  get_biomass(
    dbh = scbi_stem1$dbh,
    genus = scbi_stem1$genus,
    species = scbi_stem1$species,
    coords = c(-78.2, 38.9)
  )
```

Biomass for a single tree can be estimated given dbh and species
identification (results in kilograms).

``` r
get_biomass(
  dbh = 50,
  genus = "liriodendron",
  species = "tulipifera",
  coords = c(-78.2, 38.9)
)
#> [1] 1578.644
```

Users can modify the set of equations that will be used to estimate the
biomass using the `new_equations` function. The default option is the
entire *allodb* equation table. Users can also work on a subset of those
equations, or add new equations to the table (see
`?allodb::new_equations`). This new equation table should be provided as
an argument in the `get_biomass` function.

``` r
show_cols <- c("equation_id", "equation_taxa", "equation_allometry")
eq_tab_acer <- new_equations(subset_taxa = "Acer")
head(eq_tab_acer[, show_cols])
#> # A tibble: 6 × 3
#>   equation_id equation_taxa       equation_allometry                            
#>   <chr>       <chr>               <chr>                                         
#> 1 a4e4d1      Acer saccharum      exp(-2.192-0.011*dbh+2.67*(log(dbh)))         
#> 2 dfc2c7      Acer rubrum         2.02338*(dbh^2)^1.27612                       
#> 3 eac63e      Acer rubrum         5.2879*(dbh^2)^1.07581                        
#> 4 f49bcb      Acer pseudoplatanus exp(-5.644074+(2.5189*(log(pi*dbh))))         
#> 5 14bf3d      Acer mandshuricum   0.0335*(dbh)^1.606+0.0026*(dbh)^3.323+0.1222*…
#> 6 0c7cd6      Acer mono           0.0202*(dbh)^1.810+0.0111*(dbh)^2.740+0.1156*…
```

Within the `get_biomass` function, this equation table is used to
calibrate a new allometric equation for all species/site combinations in
the user-provided dataframe. This is done by attributing a weight to
each equation based on its sampling size, and taxonomic and climatic
similarity with the species/site combination considered.

``` r
allom_weights <-
  weight_allom(
    genus = "Acer",
    species = "rubrum",
    coords = c(-78, 38)
  )

## visualize weights
equ_tab_acer <- new_equations()
equ_tab_acer$weights <- allom_weights
keep_cols <-
  c(
    "equation_id",
    "equation_taxa",
    "sample_size",
    "weights"
  )
order_weights <- order(equ_tab_acer$weights, decreasing = TRUE)
equ_tab_acer <- equ_tab_acer[order_weights, keep_cols]
head(equ_tab_acer)
#> # A tibble: 6 × 4
#>   equation_id equation_taxa        sample_size weights
#>   <chr>       <chr>                      <dbl>   <dbl>
#> 1 138258      Acer rubrum                  150   0.415
#> 2 d6be5c      Sapindaceae                  243   0.383
#> 3 a2fbbb      Sapindaceae                  200   0.349
#> 4 2630d5      Trees (Angiosperms)          886   0.299
#> 5 d4c590      Trees (Angiosperms)          549   0.289
#> 6 ed748f      Broad-leaved species        2223   0.270
```

Equations are then resampled within their original DBH range: the number
of resampled values for each equation is proportional to its weight (as
attributed by the `weight_allom` function).

``` r
df_resample <-
  resample_agb(
    genus = "Acer",
    species = "rubrum",
    coords = c(-78, 38)
  )

plot(
  df_resample$dbh,
  df_resample$agb,
  xlab = "DBH (cm)",
  ylab = "Resampled AGB values (kg)"
)
```

![](man/figures/README-resample-acer-1.png)<!-- -->

The resampled values are then used to fit the following nonlinear model:
<img src="https://render.githubusercontent.com/render/math?math=AGB = a * dbh ^ b %2B e">,
with i.i.d.
<img src="https://render.githubusercontent.com/render/math?math=e ~N(0, sigma^2)">.
The parameters (*a*, *b*, and *sigma*) are returned by the
`est_params()` function.

The resampled values (dots) and new fitted equation (red dotted line)
can be visualized with the `illustrate_allodb()` function.

``` r
pars_acer <- est_params(
  genus = "Acer",
  species = "rubrum",
  coords = c(-78, 38)
)
illustrate_allodb(
  genus = "Acer",
  species = "rubrum",
  coords = c(-78, 38)
)
```

![](man/figures/README-est-params-acer-1.png)<!-- -->

The `est_params` function can be used for all species/site combinations
in the dataset at once.

``` r
params <- est_params(
  genus = scbi_stem1$genus,
  species = scbi_stem1$species,
  coords = c(-78.2, 38.9)
)
head(params)
#> # A tibble: 6 × 7
#>   genus       species      long   lat      a     b sigma
#>   <chr>       <chr>       <dbl> <dbl>  <dbl> <dbl> <dbl>
#> 1 Acer        negundo     -78.2  38.9 0.0762  2.55  433.
#> 2 Acer        rubrum      -78.2  38.9 0.0768  2.55  412.
#> 3 Ailanthus   altissima   -78.2  38.9 0.0995  2.48  377.
#> 4 Amelanchier arborea     -78.2  38.9 0.0690  2.56  359.
#> 5 Asimina     triloba     -78.2  38.9 0.0995  2.48  377.
#> 6 Carpinus    caroliniana -78.2  38.9 0.0984  2.48  317.
```

AGB is then recalculated as `agb = a * dbh^b` within the `get_biomass`
function.

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# allodb (development version)

# allodb 0.0.1

* Add support for R 3.4.

* Incorporate peer-review suggestions as per
https://github.com/ropensci/software-review/issues/436.
## Test environments

* ubuntu 18.04 (local), R 4.0.3
* ubuntu 18.04 (github actions), R 3.5, R 3.6, R-oldrel, R-release, R-devel
* macOS-latest (github actions), R-release
* windows-latest (github actions), R 3.6, R-release
* win-builder, R-release, R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
# w/ `new_taxa` and `NULL` `new_allometry` errors gracefully

    `new_allometry` must not be `NULL`
    i Did you forget to add the new allometry?

# with `new_allometry` and NULL `new_taxa` errors gracefully

    You must provide the taxa, coordinates, DBH range
             and sample size of you new allometries.

# with a `new_coords` 'matrix' 2x1 errors gracefully

    `coords` must be a numeric vector or matrix, with 2 values or columns.

# with arguments of different lenght errors gracefully

    All of these arguments must have the same length:
    * `new_taxa`
    * `new_allometry`
    * `new_min_dbh`
    * `new_max_dbh`
    * `new_sample_size`

---

    All of these arguments must have the same length:
    * `new_taxa`
    * `new_allometry`
    * `new_min_dbh`
    * `new_max_dbh`
    * `new_sample_size`

---

    All of these arguments must have the same length:
    * `new_taxa`
    * `new_allometry`
    * `new_min_dbh`
    * `new_max_dbh`
    * `new_sample_size`

---

    All of these arguments must have the same length:
    * `new_taxa`
    * `new_allometry`
    * `new_min_dbh`
    * `new_max_dbh`
    * `new_sample_size`

# if `new_allometry` isn't of type character errors gracefully

    The equation allometry should be a character vector.

# if `new_allometry` conains an assignment errors gracefully

    `new_allometry` must be a function of dbh (e.g. '0.5 * dbh^2').

---

    `new_allometry` must be a function of dbh (e.g. '0.5 * dbh^2').

# height must be in meters

    Height allometries outputs must be in 'm'.

---

    Height allometries outputs must be in 'm'.

# with bad coordinates errors gracefully

    Longitude must be between -180 and 180, and latitude between 90 and 0.

# with equation not a function of DBH errors gracefully

    Each new allometry must contain DBH as a dependent variable.

# with bad `new_unit_dbh` throws no error

    `new_unit_dbh` must be in 'cm', 'mm' or 'inch'.

# with bad `new_unit_output` throws no error

    `new_unit_output` must be 'g', 'kg', 'Mg' or 'lbs', or 'm'.

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

Please note that the {{{ package }}} project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
The file data-raw/available_random_ids.csv contains random ids that you can use to identify new equations. Every time you use an id from this list, you should remove it from this list.

This chunk shows how the random ids were created.

```
# WARNING
# This chunk should not be re-run. If you re-run it, you will change the random
# ids. 
available_ids <- tibble::tibble(random_id = ids::random_id(2000, bytes = 3))
write_csv(available_ids, here("data-raw/available_random_ids.csv"))
```
