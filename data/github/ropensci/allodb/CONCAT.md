
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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  fig.width = 6,
  fig.height = 5
)
```

# <img src="https://i.imgur.com/39pvr4n.png" align="left" height=44 /> allodb: An R package for biomass estimation at extratropical forest plots

<!-- badges: start -->
[![peer-review](https://badges.ropensci.org/436_status.svg)](https://github.com/ropensci/software-review/issues/436)
[![Codecov test coverage](https://codecov.io/gh/ropensci/allodb/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/allodb?branch=master)
[![R-CMD-check](https://github.com/ropensci/allodb/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/allodb/actions)
<!-- badges: end -->

## Introduction

Allometric equations for calculation of tree aboveground biomass (AGB) form the basis for estimates of forest carbon storage and exchange with the atmosphere. While standard models exist to calculate forest biomass across the tropics, we lack a standardized tool for computing AGB across the global extratropics.

_allodb_ was conceived as a framework to standardize and simplify the biomass estimation process across globally distributed extratropical forests (mainly temperate and boreal forests). 
With _allodb_ we aimed to: a) compile relevant published and unpublished allometries, focusing on AGB but structured to handle other variables (e.g., height); b) objectively select and integrate appropriate available equations across the full range  of tree sizes; and c) serve as a platform for future updates and expansion to other research sites.

The _allodb_ package contains a dataset of systematically selected published allometric equations. This dataset was built based on 701 woody species identified at 24 large [ForestGEO forest dynamic plots](https://forestgeo.si.edu/) representing all major extratropical forest types. A total of 570 parsed allometric equations to estimate individual tree biomass were retrieved, checked, and combined using a weighting function designed to ensure optimal equation selection over the full tree size range with smooth transitions across equations. The equation dataset used can be customized with built-in functions that subset the original dataset and add new equations.

The package provides functions to estimate tree biomass based on user-provided census data (tree diameter, taxonomic identification, and plot coordinates). New allometric equations are calibrated for each species and location by resampling the original equations; equations with a larger sample size and/or higher taxonomic and climatic similarity with the species and location in question are given a higher weight in this process. 

## Installation

Install the development version of _allodb_ from GitHub:

```R
# install.packages("remotes")
remotes::install_github("ropensci/allodb")
```

## Examples

Prior to calculating tree biomass using _allodb_, users need to provide a table (i.e. dataframe) with DBH (cm), parsed species Latin names, and site(s) coordinates. In the following examples we use data from the Smithsonian Conservation Biology Institute, USA (SCBI) ForestGEO dynamics plot (trees from 1 hectare surveyed in 2008). Full tree census data can be requested through the [ForestGEO portal](https://forestgeo.si.edu/explore-data).

```{r open-data}
library(allodb)
data(scbi_stem1)
``` 

The biomass of all trees in one (or several) censuses can be estimated using the `get_biomass` function. 

```{r calc-agb-all}
scbi_stem1$agb <-
  get_biomass(
    dbh = scbi_stem1$dbh,
    genus = scbi_stem1$genus,
    species = scbi_stem1$species,
    coords = c(-78.2, 38.9)
  )
```

Biomass for a single tree can be estimated given dbh and species identification (results in kilograms).

```{r calc-agb-poplar}
get_biomass(
  dbh = 50,
  genus = "liriodendron",
  species = "tulipifera",
  coords = c(-78.2, 38.9)
)
```

Users can modify the set of equations that will be used to estimate the biomass using the `new_equations` function. The default option is the entire _allodb_ equation table. Users can also work on a subset of those equations, or add new equations to the table (see `?allodb::new_equations`). This new equation table should be provided as an argument in the `get_biomass` function.  

```{r}
show_cols <- c("equation_id", "equation_taxa", "equation_allometry")
eq_tab_acer <- new_equations(subset_taxa = "Acer")
head(eq_tab_acer[, show_cols])
```

Within the `get_biomass` function, this equation table is used to calibrate a new allometric equation for all species/site combinations in the user-provided dataframe. This is done by attributing a weight to each equation based on its sampling size, and taxonomic and climatic similarity with the species/site combination considered. 

```{r weights}
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
```

Equations are then resampled within their original DBH range: the number of resampled values for each equation is proportional to its weight (as attributed by the `weight_allom` function). 

```{r resample-acer, eval = TRUE}
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


The resampled values are then used to fit the following nonlinear model: <img src="https://render.githubusercontent.com/render/math?math=AGB = a * dbh ^ b %2B e">, with i.i.d. <img src="https://render.githubusercontent.com/render/math?math=e ~N(0, sigma^2)">. The parameters (_a_, _b_, and _sigma_) are returned by the `est_params()` function.

The resampled values (dots) and new fitted equation (red dotted line) can be visualized with the `illustrate_allodb()` function. 

```{r est-params-acer, eval = TRUE, fig.height=4, fig.width=10}
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


The `est_params` function can be used for all species/site combinations in the dataset at once. 

```{r est-params-all, eval = TRUE}
params <- est_params(
  genus = scbi_stem1$genus,
  species = scbi_stem1$species,
  coords = c(-78.2, 38.9)
)
head(params)
```

AGB is then recalculated as `agb = a * dbh^b` within the `get_biomass` function.

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.
---
title: "Using allodb to estimate aboveground biomass"
output: rmarkdown::html_vignette
vignette: >
 %\VignetteIndexEntry{Using allodb to estimate aboveground biomass}
 %\VignetteEngine{knitr::rmarkdown}
 %\VignetteEncoding{UTF-8}
---
  
```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE
)
```

## Installation

Install the development version of _allodb_ from GitHub:

```{r install, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/allodb")
```

## Load the census data

Prior to calculating tree biomass using _allodb_, users need to provide a table (i.e. dataframe) with DBH (cm), parsed species Latin names, and site(s) coordinates. In the following examples we use data from the Smithsonian Conservation Biology Institute, USA (SCBI) ForestGEO dynamics plot (1st census in 2008, trees from 1 hectare). Data can be requested through the ForestGEO portal (https://forestgeo.si.edu/)

```{r load-census-data}
require(allodb)
data(scbi_stem1)
str(scbi_stem1)
``` 

## Load and modify the equation table 

`allodb` provides a dataframe containing `r nrow(allodb::equations)` parsed allometric equations. 

```{r load-equations}
data(equations)
``` 

Additional information about the equation table can be found in the equation metadata. 

```{r load-metadata}
data(equations_metadata)
``` 

This equation table is the default in all functions of the package. Users can modify the set of equations that will be used to estimate biomass using the `new_equations()` function: users can work on a subset of those equations, or add new equations to the table (see `?allodb::new_equations`). The customized equation table should be provided as an argument in the `get_biomass()` (or other) function (argument name: `new_eqtable`). 

### Subset the equation table 

```{r subset-table}
show_cols <- c("equation_id", "equation_taxa", "equation_allometry")
eq_tab_acer <- new_equations(subset_taxa = "Acer")
head(eq_tab_acer[, show_cols])
```

### Add new equations

```{r add-equations}
eq_tab_add <- new_equations(
  new_taxa = c("Quercus ilex", "Castanea sativa"),
  new_allometry = c("0.12*dbh^2.5", "0.15*dbh^2.7"),
  new_coords = c(4, 44),
  new_min_dbh = c(5, 10),
  new_max_dbh = c(35, 68),
  new_sample_size = c(143, 62)
)
## show added equations - they contain "new" in their equation_id
head(eq_tab_add[grepl("new", eq_tab_add$equation_id), ])
```


## Estimate the aboveground biomass

The aboveground biomass (AGB) can be estimated using the `get_biomass()` function: the required arguments are the diameter at breast height (DBH, in cm), the taxonomic identification (to the genus or species level), and the location (long-lat coordinates). The output is the aboveground biomass of the tree in kg. 

```{r calc-agb-poplar}
get_biomass(
  dbh = 50,
  genus = "liriodendron",
  species = "tulipifera",
  coords = c(-78.2, 38.9)
)
```

The `get_biomass()` function can also be used to estimate the AGB of all trees in one (or several) censuses. 

```{r calc-agb-all}
scbi_stem1$agb <-
  get_biomass(
    dbh = scbi_stem1$dbh,
    genus = scbi_stem1$genus,
    species = scbi_stem1$species,
    coords = c(-78.2, 38.9)
  )
```

```{r plot-agb-scbi, fig.width = 6, fig.height = 5}
plot(
  x = scbi_stem1$dbh,
  y = scbi_stem1$agb,
  col = factor(scbi_stem1$genus),
  xlab = "DBH (cm)",
  ylab = "AGB (kg)"
)
```

## How AGB is estimated 

### Attribute a weight to each equation in the equation table for each taxon/location combination

Within the `get_biomass()` function, new allometric equations are calibrated for each taxon/location combinations in the user-provided dataframe. This is done by attributing a weight to each equation in the equation table, based on its sampling size, and taxonomic and climatic similarity with the taxon/location combination considered. 

```{r weights}
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
```

### Resample equations 

Equations are then resampled within their original DBH range: the number of resampled values for each equation is proportional to its weight (as attributed by the `weight_allom()` function). 

```{r resample-acer, fig.width = 6, fig.height = 5}
df_resample <-
  resample_agb(
    genus = "Acer",
    species = "rubrum",
    coords = c(-78, 38),
    nres = 1e4
  )

plot(
  df_resample$dbh,
  df_resample$agb,
  xlab = "DBH (cm)",
  ylab = "Resampled AGB values (kg)"
)
```

### Calibrate a new equation for each taxon/location combination

The resampled values are then used to fit the following nonlinear model: $AGB = a \cdot dbh ^ b + e$, with i.i.d. $e \sim \mathcal{N}(0, sigma^2)$. The parameters (_a_, _b_, and _sigma_) are returned by the `est_params()` function. In other words, this function calibrates new allometric equations from sampling previous ones. New allometric equations are calibrated for each species and location by resampling the original compiled equations; equations with a larger sample size, and/or higher taxonomic rank, and climatic similarity with the species and location in question are given a higher weight in this process.

```{r est-params-acer, fig.width = 6, fig.height = 5}
pars_acer <- est_params(
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
curve(pars_acer$a * x^pars_acer$b,
  add = TRUE, col = 2, lwd = 2
)
```

The `est_params()` function can be used for all species/site combinations in the dataset at once. 

```{r est-params-all}
params <- est_params(
  genus = scbi_stem1$genus,
  species = scbi_stem1$species,
  coords = c(-78.2, 38.9)
)
head(params)
```

AGB is then recalculated as `AGB = a * dbh^b` within the `get_biomass()` function. 

## Visualize the recalibration of equations

The recalibrated equation for one taxon/location combination can be easily visualized with the `illustrate_allodb()` function, which returns a ggplot (see package `ggplot2`) with all resampled values, and the top equations used displayed in the legend. The user can control the number of equations and equation information shown in the legend. The red dotted line is the recalibrated equation used in the `get_biomass()` function. 

```{r illustrate-allodb, fig.height = 4, fig.width=8}
illustrate_allodb(
  genus = "Acer",
  species = "rubrum",
  coords = c(-78, 38)
)
```

```{r illustrate-allodb2, fig.height = 5, fig.width=10}
illustrate_allodb(
  genus = "Acer",
  species = "rubrum",
  coords = c(-78, 38),
  neq = 15,
  eqinfo = c("equation_taxa", "geographic_area")
)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{shrub_species}
\alias{shrub_species}
\title{Shrub species identified in selected ForestGEO sites}
\format{
An object of class \code{character} of length 1.
}
\usage{
shrub_species
}
\description{
Genus and species of shrubby plants identified in the 24 extratropical
ForestGEO sites used in allodb.
}
\examples{
# preview the dataset
print(head(shrub_species))
}
\seealso{
Other datasets: 
\code{\link{genus_family}},
\code{\link{gymno_genus}},
\code{\link{koppenMatrix}},
\code{\link{scbi_stem1}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_biomass.R
\name{get_biomass}
\alias{get_biomass}
\title{Compute tree aboveground biomass (AGB) based on allometric equations}
\usage{
get_biomass(
  dbh,
  genus,
  coords,
  species = NULL,
  new_eqtable = NULL,
  wna = 0.1,
  w95 = 500,
  nres = 10000
)
}
\arguments{
\item{dbh}{a numeric vector containing tree diameter at breast height (dbh)
measurements, in cm.}

\item{genus}{a character vector (same length as dbh), containing the genus
(e.g. "Quercus") of each tree.}

\item{coords}{a numeric vector of length 2 with longitude and latitude (if
all trees were measured in the same location) or a matrix with 2 numerical
columns giving the coordinates of each tree.}

\item{species}{a character vector (same length as dbh), containing the
species (e.g. "rubra")  of each tree. Default is \code{NULL}, when no species
identification is available.}

\item{new_eqtable}{Optional. An equation table created with the
\code{\link[=new_equations]{new_equations()}} function.}

\item{wna}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the dbh-related weight attributed to equations
without a specified dbh range. Default is 0.1.}

\item{w95}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the value at which the sample-size-related weight
reaches 95\% of its maximum value (max=1). Default is 500.}

\item{nres}{number of resampled values. Default is "1e4".}
}
\value{
A "numeric" vector of the same length as dbh, containing AGB value
(in kg) for every stem.
}
\description{
This function calculates the aboveground biomass (or other tree components)
of a given tree based on published allometric equations. Users need to
provide a table (i.e. dataframe) with DBH (cm), parsed species Latin names,
and site(s) coordinates. The biomass of all trees in one (or several)
censuses can be estimated using this function.
}
\details{
\code{allodb} estimates AGB by calibrating a new allometric equation for each
taxon (arguments \code{genus} and  \code{species}) and location (argument \code{coords}) in
the user-provided census data. The new allometric equation is based on a set
of allometric equations that can be customized using the \code{new_eqtable}
argument. Each equation is then given a weight with the \code{\link[=weight_allom]{weight_allom()}}
function, based on: 1) its original sample size (numbers of trees used to
develop a given allometry), 2) its climatic similarity with the target
location, and 3) its taxonomic similarity with the target taxon (see
documentation of the \code{\link[=weight_allom]{weight_allom()}} function). The final weight attributed
to each equation is the product of those three weights. Equations are then
resampled with the\code{\link[=resample_agb]{resample_agb()}} funtion: the number of samples per
equation is proportional to its weight, and the total number of samples is
provided by the argument \code{nres}. The resampling is done by drawing DBH values
from a uniform distribution on the DBH range of the equation, and estimating
the AGB with the equation. The couples of values (DBH, AGB) obtained are then
used in the function \code{\link[=est_params]{est_params()}} to calibrate a new allometric equation,
by applying a linear regression to the log-transformed data. The parameters
of the new allometric equations are then used in the \code{\link[=get_biomass]{get_biomass()}} function
by back-transforming the AGB predictions based on the user-provided DBHs.
}
\section{Warning}{

The function can run into some memory problems when used on large datasets
(usually several hundred thousand observations).
}

\examples{
# Estimate biomass of all individuals from the Lauraceae family at the SCBI
# plot
lau <- subset(scbi_stem1, Family == "Lauraceae")
lau$agb <- get_biomass(lau$dbh, lau$genus, lau$species,
  coords = c(-78.2, 38.9)
)
lau

# Estimate biomass from multiple sites (using scbi_stem1 as example with
# multiple coord)
dat <- scbi_stem1[1:100, ]
dat$long <- c(rep(-78, 50), rep(-80, 50))
dat$lat <- c(rep(40, 50), rep(41, 50))
dat$biomass <- get_biomass(
  dbh = dat$dbh,
  genus = dat$genus,
  species = dat$species,
  coords = dat[, c("long", "lat")]
)
dat
}
\seealso{
\code{\link[=weight_allom]{weight_allom()}}, \code{\link[=new_equations]{new_equations()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{references}
\alias{references}
\alias{references_metadata}
\title{Equation references and associated metadata}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 57 rows and 6 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 7 rows and 4 columns.
}
\usage{
references

references_metadata
}
\description{
\itemize{
\item \link{references}: A data frame listing all references used in \code{equation} table.
\item \link{references_metadata}: Metadata for \code{reference} table.
}
}
\details{
Bibliographical information for sourced equations. Links to the \link{equations}
table by \code{ref_id}.
}
\examples{
# preview the datasets
print(head(references))
print(head(references_metadata))
}
\seealso{
Other database datasets: 
\code{\link{equations}},
\code{\link{missing_values}},
\code{\link{sites_info}},
\code{\link{sitespecies}}
}
\concept{database datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/illustrate_allodb.R
\name{illustrate_allodb}
\alias{illustrate_allodb}
\title{Illustrate the resampling of AGB values used in \emph{allodb}}
\usage{
illustrate_allodb(
  genus,
  coords,
  species = NULL,
  new_eqtable = NULL,
  logxy = FALSE,
  neq = 10,
  eqinfo = "equation_taxa",
  wna = 0.1,
  w95 = 500,
  nres = 10000
)
}
\arguments{
\item{genus}{A character value, containing the genus (e.g. "Quercus") of the
tree.}

\item{coords}{A numeric vector of length 2 with longitude and latitude.}

\item{species}{A character value, containing the species (e.g. "rubra") of
the tree. Default is \code{NULL}, when no species identification is available.}

\item{new_eqtable}{Optional. An equation table created with the
\code{\link[=new_equations]{new_equations()}} function. Default is the base \emph{allodb} equation
table.}

\item{logxy}{Logical: should values be plotted on a log scale? Default is
\code{FALSE}.}

\item{neq}{Number of top equations in the legend. Default is 10, meaning that
the 10 equations with the highest weights are shown in the legend.}

\item{eqinfo}{Which column(s) of the equation table should be used in the
legend? Default is "equation_taxa".}

\item{wna}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the dbh-related and sample-size related weights
attributed to equations without a specified dbh range or sample size,
respectively. Default is 0.1.}

\item{w95}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the value at which the sample-size-related weight
reaches 95\% of its maximum value (max=1). Default is 500.}

\item{nres}{number of resampled values. Default is "1e4".}
}
\value{
An object of class "ggplot" showing all resampled dbh-agb values. The
top equations used are shown in the legend. The red curve on the graph
represents the final fitted equation.
}
\description{
This function illustrates the resampling of AGB values used in \emph{allodb}. It
creates objects of class "ggplot".
}
\examples{
illustrate_allodb(
  genus = "Quercus",
  species = "rubra",
  coords = c(-78.2, 38.9)
)
}
\seealso{
\code{\link[=weight_allom]{weight_allom()}}, \code{\link[=new_equations]{new_equations()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_params.R
\name{est_params}
\alias{est_params}
\title{Calibrate new allometric equations}
\usage{
est_params(
  genus,
  coords,
  species = NULL,
  new_eqtable = NULL,
  wna = 0.1,
  w95 = 500,
  nres = 10000
)
}
\arguments{
\item{genus}{a character vector, containing the genus (e.g. "Quercus") of
each tree.}

\item{coords}{a numeric vector of length 2 with longitude and latitude (if
all trees were measured in the same location) or a matrix with 2 numerical
columns giving the coordinates of each tree.}

\item{species}{a character vector (same length as genus), containing the
species (e.g. "rubra")  of each tree. Default is \code{NULL}, when no species
identification is available.}

\item{new_eqtable}{Optional. An equation table created with the
\code{\link[=new_equations]{new_equations()}} function. Default is the compiled \emph{allodb} equation
table.}

\item{wna}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the dbh-related and sample-size related weights
attributed to equations without a specified dbh range or sample size,
respectively. Default is 0.1.}

\item{w95}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the value at which the sample-size-related weight
reaches 95\% of its maximum value (max=1). Default is 500.}

\item{nres}{number of resampled values. Default is "1e4".}
}
\value{
An object of class "data.frame" of fitted coefficients (columns) of
the non-linear least-square regression:
\deqn{AGB = a * dbh ^ b + e, \space
  \mathit{with} \space e ~ N(0, sigma^2)}{AGB = a * dbh ^ b + e, with e ~
  N(0, sigma^2)}
}
\description{
This function calibrates new allometric equations from sampling previous
ones. New allometric equations are calibrated for each species and location
by resampling the original compiled equations; equations with a larger sample
size, and/or higher taxonomic rank, and climatic similarity with the species
and location in question are given a higher weight in this process.
}
\examples{
# calibrate new allometries for all Lauraceae species
lauraceae <- subset(scbi_stem1, Family == "Lauraceae")
est_params(
  genus = lauraceae$genus,
  species = lauraceae$species,
  coords = c(-78.2, 38.9)
)
}
\seealso{
\code{\link[=weight_allom]{weight_allom()}}, \code{\link[=new_equations]{new_equations()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{genus_family}
\alias{genus_family}
\title{Genus and family table for selected ForestGEO sites}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 248 rows and 2 columns.
}
\usage{
genus_family
}
\description{
Genus and their associated family identified in the extratropical ForestGEO
sites used in allodb.
}
\examples{
# preview the dataset
print(head(genus_family))
}
\seealso{
Other datasets: 
\code{\link{gymno_genus}},
\code{\link{koppenMatrix}},
\code{\link{scbi_stem1}},
\code{\link{shrub_species}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{missing_values}
\alias{missing_values}
\title{Explanations of missing values codes}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 4 rows and 3 columns.
}
\usage{
missing_values
}
\description{
Explanation of the codes used to indicate missing information in equation
table.
}
\examples{
# preview the dataset
print(head(missing_values))
}
\seealso{
Other database datasets: 
\code{\link{equations}},
\code{\link{references}},
\code{\link{sites_info}},
\code{\link{sitespecies}}
}
\concept{database datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sitespecies}
\alias{sitespecies}
\alias{sitespecies_metadata}
\title{Sites and tree species used in allodb and associated metadata}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 1113 rows and 10 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 10 rows and 4 columns.
}
\usage{
sitespecies

sitespecies_metadata
}
\description{
\itemize{
\item \link{sitespecies}: Table of extratropical ForestGEO sites in allodb (n=24) and
their tree species.
\item \link{sitespecies_metadata}: Metadata for \link{sitespecies} table.
}
}
\examples{
# preview the datasets
print(head(sitespecies))
print(head(sitespecies_metadata))
}
\seealso{
Other database datasets: 
\code{\link{equations}},
\code{\link{missing_values}},
\code{\link{references}},
\code{\link{sites_info}}
}
\concept{database datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weight_allom.R
\name{weight_allom}
\alias{weight_allom}
\title{Attribute weights to equations}
\usage{
weight_allom(
  genus,
  coords,
  species = NULL,
  new_eqtable = NULL,
  wna = 0.1,
  w95 = 500
)
}
\arguments{
\item{genus}{a character value, containing the genus (e.g. "Quercus") of the
tree.}

\item{coords}{a numeric vector of length 2 with longitude and latitude.}

\item{species}{a character vector (same length as genus), containing the
species (e.g. "rubra") of the tree. Default is \code{NULL}, when no species
identification is available.}

\item{new_eqtable}{Optional. An equation table created with the
\code{\link[=new_equations]{new_equations()}} function.}

\item{wna}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the sample-size related weights attributed to
equations without a specified sample size. Default is 0.1.}

\item{w95}{a numeric vector, this parameter is used to determine the value at
which the sample-size-related weight reaches 95\% of its maximum value
(max=1). Default is 500.}
}
\value{
A named "numeric" vector with one weight for each equation.
}
\description{
This function attributes a weight to each equation based on its sampling
size, and taxonomic and climatic similarity with the species/site combination
considered.
}
\details{
Each equation is given a weight by the function \code{\link[=weight_allom]{weight_allom()}}, calculated
as the product of the following components:

(1) sample-size weight, calculated as:

\deqn{1-exp(-n*(log(20)/w95))}

where n is the sample size of the equation; the weight given to equations
with no sample size information is determined by argument \code{wna} (0.1 by
default).

(2) climate weight, based on the similarity between the climatic conditions
of the equation site and the target location, using the three-letter system
of Koppen's climate scheme. Climate weights associated with each combination
of two Koppen climates are provided in \code{data("koppenMatrix")}. The resulting
weight has a value between 1e-6 (different climate groups) and 1 (exactly the
same climate classification). When an equation was calibrated with trees from
several locations with different Koppen climates, the maximum value out of
all pairwise equation-site climate weights is used.

(3) taxonomic weight: equal to 1 for same species equations, 0.8 for same
genus equations, 0.5 for same family equations and for equations calibrated
for the same broad functional or taxonomic group (e.g. shrubs, conifers,
angiosperms). All other equations are given a low taxonomic weight of 1e-6:
these equations will have a significant relative weight in the final
prediction only when no other more specific equation is available.
}
\examples{
x <- weight_allom(
  genus = "Acer",
  species = "negundo",
  coords = c(-78.2, 38.9)
)
str(x)
head(x)
}
\seealso{
\code{\link[=get_biomass]{get_biomass()}}, \code{\link[=new_equations]{new_equations()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_equations.R
\name{new_equations}
\alias{new_equations}
\title{Modify the original equation table}
\usage{
new_equations(
  subset_taxa = "all",
  subset_climate = "all",
  subset_region = "all",
  subset_ids = "all",
  subset_output = c("Total aboveground biomass", "Whole tree (above stump)"),
  new_taxa = NULL,
  new_allometry = NULL,
  new_coords = NULL,
  new_min_dbh = NULL,
  new_max_dbh = NULL,
  new_sample_size = NULL,
  new_unit_dbh = "cm",
  new_unit_output = "kg",
  new_input_var = "DBH",
  new_output_var = "Total aboveground biomass",
  use_height_allom = TRUE
)
}
\arguments{
\item{subset_taxa}{character vector with taxa to be kept. Default is "all",
in which case all taxa are kept.}

\item{subset_climate}{character vector with Koppen climate classification to
be kept. Default is "all", in which case all climates are kept.}

\item{subset_region}{character vector with name of location(s) or
country(ies) or broader region(s) (eg. "Europe", "North America") to be
kept. Default is "all", in which case all regions/countries are kept.}

\item{subset_ids}{character vector with equation IDs to be kept. Default is
"all", in which case all equations are kept.}

\item{subset_output}{What dependent variable(s) should be provided in the
output? Default is "Total aboveground biomass" and "Whole tree (above
stump)", other possible values are: "Bark biomass", "Branches (dead)",
"Branches (live)", "Branches total (live, dead)", "Foliage total",
"Height", "Leaves", "Stem (wood only)", "Stem biomass", "Stem biomass (with
bark)", "Stem biomass (without bark)", "Whole tree (above and
belowground)". Be aware that currently only a few equations represent those
other variables, so estimated values might not be very accurate.}

\item{new_taxa}{character string or vector specifying the taxon (or taxa) for
which the allometry has been calibrated.}

\item{new_allometry}{a character string with the allometric equation.}

\item{new_coords}{a vector or matrix of coordinates (longitude, latitude) of
the calibration data.}

\item{new_min_dbh}{numerical value, minimum DBH for which the equation is
valid (in cm). Default is \code{NULL} (nothing is added).}

\item{new_max_dbh}{numerical value, maximum DBH for which the equation is
valid (in cm). Default is \code{NULL} (nothing is added).}

\item{new_sample_size}{number of measurements with which the allometry was
calibrated. Default is \code{NULL} (nothing is added).}

\item{new_unit_dbh}{character string with unit of DBH in the equation (either
\code{cm}, \code{mm} or \code{inch}). Default is "cm".}

\item{new_unit_output}{character string with unit of equation output (either
"g", "kg", "Mg" or "lbs" if the output is a mass, or "m" if the output is a
height).}

\item{new_input_var}{independent variable(s) needed in the allometry. Default
is "DBH", other option is "DBH, H".}

\item{new_output_var}{dependent variable estimated by the allometry. Default
is "Total aboveground biomass".}

\item{use_height_allom}{a logical value. In \emph{allodb} we use Bohn et al.
(2014) for European sites. User need to provide height allometry when
needed to calculate AGB. Default is \code{TRUE}.}
}
\value{
An object of class "data.frame" of new equations.
}
\description{
This function modifies the original equation table to be used in other
functions of the package including: subset the original equation table, add
new equations, and choose whether to include equations with a height
allometry.
}
\examples{
new_equations(
  new_taxa = "Faga",
  new_allometry = "exp(-2+log(dbh)*2.5)",
  new_coords = c(-0.07, 46.11),
  new_min_dbh = 5,
  new_max_dbh = 50,
  new_sample_size = 50
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{equations}
\alias{equations}
\alias{equations_metadata}
\title{Tables of allometric equations and associated metadata}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 570 rows and 47 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 47 rows and 7 columns.
}
\source{
See \link{references} for equations original sources.
}
\usage{
equations

equations_metadata
}
\description{
\itemize{
\item \link{equations}: Table of allometric equations.
\item \link{equations_metadata}: Explanation of columns for \link{equations} table.
}
}
\details{
A compilation of best available allometry equations to calculate tree
above-ground biomass (AGB) per species based on extratropical ForestGEO
sites.
}
\examples{
# preview the datasets
print(head(equations))
print(head(equations_metadata))
}
\seealso{
Other database datasets: 
\code{\link{missing_values}},
\code{\link{references}},
\code{\link{sites_info}},
\code{\link{sitespecies}}
}
\concept{database datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resample_agb.R
\name{resample_agb}
\alias{resample_agb}
\title{Resample \emph{allodb} equations to calibrate new allometries}
\usage{
resample_agb(
  genus,
  coords,
  species = NULL,
  new_eqtable = NULL,
  wna = 0.1,
  w95 = 500,
  nres = 10000
)
}
\arguments{
\item{genus}{a character value, containing the genus (e.g. "Quercus") of the
tree.}

\item{coords}{a numeric vector of length 2 with longitude and latitude.}

\item{species}{a character value, containing the species (e.g. "rubra") of
the tree. Default is "NULL", when no species identification is available.}

\item{new_eqtable}{Optional. An equation table created with the
\code{\link[=new_equations]{new_equations()}} function. Default is the original \emph{allodb} equation
table.}

\item{wna}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the dbh-related and sample-size related weights
attributed to equations without a specified dbh range or sample size,
respectively. Default is 0.1.}

\item{w95}{a numeric vector, this parameter is used in the \code{\link[=weight_allom]{weight_allom()}}
function to determine the value at which the sample-size-related weight
reaches 95\% of its maximum value (max=1). Default is 500.}

\item{nres}{number of resampled values. Default is "1e4".}
}
\value{
An object of class "data.frame" of resampled DBHs and associated AGB
from the equation table; the number of  resampled DBHs is proportional to
the weight provided by the \code{\link[=weight_allom]{weight_allom()}} function.
}
\description{
After attributing a weight to each equation in \emph{allodb} using the
\code{\link[=weight_allom]{weight_allom()}} function, equations are then resampled within their original
DBH range using \code{\link[=resample_agb]{resample_agb()}}: the number of resampled values for each
equation is proportional to its weight. It creates S3 objects of class
"numeric".
}
\examples{
resample_agb(
  genus = "Quercus",
  species = "rubra",
  coords = c(-78.2, 38.9)
)
}
\seealso{
\code{\link[=weight_allom]{weight_allom()}}, \code{\link[=new_equations]{new_equations()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allodb-package.R
\docType{package}
\name{allodb-package}
\alias{allodb}
\alias{allodb-package}
\title{allodb: Tree Biomass Estimation at Extra-Tropical Forest Plots}
\description{
Standardize and simplify the tree biomass estimation process across globally distributed extratropical forests.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/allodb/}
  \item \url{https://github.com/ropensci/allodb}
  \item Report bugs at \url{https://github.com/ropensci/allodb/issues}
}

}
\author{
\strong{Maintainer}: Erika Gonzalez-Akre \email{GonzalezEB@si.edu} (\href{https://orcid.org/0000-0001-8305-6672}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Camille Piponiot \email{camille.piponiot@gmail.com} (\href{https://orcid.org/0000-0002-3473-1982}{ORCID})
  \item Mauro Lepore \email{maurolepore@gmail.com} (\href{https://orcid.org/0000-0002-1986-7988}{ORCID})
  \item Kristina Anderson-Teixeira \email{TeixeiraK@si.edu} (\href{https://orcid.org/0000-0001-8461-9713}{ORCID})
}

Other contributors:
\itemize{
  \item Jeffrey Hanson \email{jeffrey.hanson@uqconnect.edu.au} [reviewer]
  \item Jonas Stillhard \email{jonas.stillhard@wsl.ch} [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{scbi_stem1}
\alias{scbi_stem1}
\title{Tree census data from SCBI ForestGEO plot}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 2287 rows and 6 columns.
}
\source{
Full datasets for tree census data at SCBI can be requested through
the ForestGEO portal (\url{https://forestgeo.si.edu/}). Census 1, 2, and 3 can
also be accessed at the public GitHub repository for SCBI-ForestGEO Data
(\url{https://github.com/SCBI-ForestGEO}).
}
\usage{
scbi_stem1
}
\description{
A table with tree data from the Smithsonian Conservation Biology Institute,
USA (SCBI) ForestGEO dynamics plot. This dataset is an extract from the first
tree census in 2008, only covering 1 hectare (SCBI plot is 25.6 ha). DBH in
cm.
}
\examples{
# preview the datasets
print(head(scbi_stem1))
}
\seealso{
Other datasets: 
\code{\link{genus_family}},
\code{\link{gymno_genus}},
\code{\link{koppenMatrix}},
\code{\link{shrub_species}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{koppenMatrix}
\alias{koppenMatrix}
\title{Koppen climate classification matrix}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 900 rows and 3 columns.
}
\usage{
koppenMatrix
}
\description{
A table built to facilitate the comparison between the Koppen climate of a
site and the allometric equation in question.
}
\details{
The value of column \code{we} is the weight given to the combination of Koppen
climates in columns \code{zone1} and \code{zone2}; the table is symmetric: \code{zone1} and
\code{zone2} can be interchanged. This weight is calculated in 3 steps: (1) if the
main climate group (first letter) is the same, the climate weight starts at
0.4; if one of the groups is "C" (temperate climate) and the other is "D"
(continental climate), the climate weight starts at 0.2 because the 2 groups
are considered similar enough; otherwise, the weight is 0; (2) if the
equation and site belong to the same group, the weight is incremented by an
additional value between 0 and 0.3 based on precipitation pattern similarity
(second letter of the Koppen zone), and (3) by an additional value between 0
and 0.3 based on temperature pattern similarity (third letter of the Koppen
zone). The resulting weight has a value between 0 (different climate groups)
and 1 (exactly the same climate classification).
}
\examples{
# preview the dataset
print(head(koppenMatrix))
}
\seealso{
Other datasets: 
\code{\link{genus_family}},
\code{\link{gymno_genus}},
\code{\link{scbi_stem1}},
\code{\link{shrub_species}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gymno_genus}
\alias{gymno_genus}
\title{Gymnosperms identified in selected ForestGEO sites}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 95 rows and 3 columns.
}
\usage{
gymno_genus
}
\description{
Table with genus and their associated family for Gymnosperms identified in
the ForestGEO sites used in allodb.
}
\examples{
# preview the dataset
print(head(gymno_genus))
}
\seealso{
Other datasets: 
\code{\link{genus_family}},
\code{\link{koppenMatrix}},
\code{\link{scbi_stem1}},
\code{\link{shrub_species}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sites_info}
\alias{sites_info}
\title{ForestGEO sites used in allodb}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 24 rows and 6 columns.
}
\usage{
sites_info
}
\description{
Table with geographical information for extratropical ForestGEO sites used in
allodb (n=24).
}
\details{
More details on geographical aspects of these ForestGEO sites can be found in
the accompanying manuscript.
}
\examples{
# preview the datasets
print(head(sites_info))
}
\seealso{
Other database datasets: 
\code{\link{equations}},
\code{\link{missing_values}},
\code{\link{references}},
\code{\link{sitespecies}}
}
\concept{database datasets}
\keyword{datasets}
