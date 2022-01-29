
<!-- README.md is generated from README.Rmd. Please edit that file -->

# grainchanger <img src="man/figures/logo.png" align="right" width="150" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/ropensci/grainchanger/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/grainchanger/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/grainchanger/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/grainchanger?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/grainchanger)](https://cran.r-project.org/package=grainchanger)
[![CranLogs](https://cranlogs.r-pkg.org/badges/grainchanger)](https://cran.r-project.org/package=grainchanger)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

The `grainchanger` package provides functionality for data aggregation
to a coarser resolution via moving-window or direct methods.

As landscape ecologists and macroecologists, we often need to aggregate
data in order to harmonise datasets. In doing so, we often lose a lot of
information about the spatial structure and environmental heterogeneity
of data measured at finer resolution. An example of this is when the
response data (e.g. species’ atlas data) are available at a coarser
resolution to the predictor data (e.g. land-use data). We developed this
method and R package in order to overcome some of these issues.

For more information on the background to and motivation for the
development of this method, see [Graham *et al.* 2019 in *Methods in
Ecology and Evolution*](https://doi.org/10.1111/2041-210X.13177).

# Package overview

The primary functions of the `grainchanger` package are those which
facilitate moving-window (`winmove_agg`) and direct (`nomove_agg`) data
aggregation. These functions aggregate fine-grain data (`fine_dat`) to a
coarse-grain (`coarse_dat`) using a function specified by the user
(`agg_fun`). The moving-window method takes in an additional function
(`win_fun`) which smooths the fine-grain data prior to aggregation.

The moving-window smoothing function is also available in the package
(`winmove`), as well as several built-in functions, and an additional
utility function for use with simulated landscapes (`create_torus`).

The `winmove` function acts as a convenient wrapper to
`raster::focalWeight` and `raster::focal` which takes advantage of
optimised functions built into the `grainchanger` package.

# Installation

    # Install release version from CRAN
    install.packages("grainchanger")
    
    # Install development version from GitHub
    devtools::install_github("ropensci/grainchanger")

# Examples

## Moving-window data aggregation

The below example shows the moving-window data aggregation in action. It
aggregates a categorical raster (`fine_dat`) to a grid using Shannon
evenness (specified by `win_fun`) as the function calculated within a
square moving window of 5 units. The value returned is the mean
(specified by `agg_fun`) of the smoothed value for each cell of
`coarse_dat`. This value is included as a column on the grid `sf`
object.

``` r
library(grainchanger)
library(ggplot2)
library(landscapetools)

# categorical landscape
show_landscape(cat_ls, discrete = TRUE)

# moving-window aggregation using Shannon evenness
g_sf$mwda <- winmove_agg(coarse_dat = g_sf,
                         fine_dat = cat_ls, 
                         d = 5,
                         type = "rectangle",
                         win_fun = shei,
                         agg_fun = mean,
                         lc_class = 1:4,
                         quiet = TRUE)

ggplot(g_sf) + 
  geom_sf(aes(fill = mwda)) + 
  scale_fill_viridis_c() +
  theme_bw()
```

<img src="man/figures/README-mwda_example-1.png" width="100%" /><img src="man/figures/README-mwda_example-2.png" width="100%" />

## Direct data aggregation

The below example shows the direct data aggregation in action. It
aggregates a continuous raster to a raster with a coarser resolution
using the range as the function calculated for each cell of the larger
grid. The resulting output is a raster of the coarser resolution.
`var_range` is an inbuilt function in the `grainchanger` package.

``` r
library(raster)

# continuous landscape
show_landscape(cont_ls)

# load the coarse resolution raster
g_raster <- raster(system.file("raster/g_raster.tif", package = "grainchanger"))

# direct aggregation using range
dda <- nomove_agg(coarse_dat = g_raster,
                       fine_dat = cont_ls, 
                       agg_fun = var_range)

show_landscape(dda)
```

<img src="man/figures/README-dda_example-1.png" width="100%" /><img src="man/figures/README-dda_example-2.png" width="100%" />

# Functions

There are a number of inbuilt functions in the grainchanger package,
with their usage outlined below. While it is possible to use
user-defined functions within both `winmove_agg` and `nomove_agg`, we
welcome suggestions for additional functions. Please [add as an
issue](https://github.com/ropensci/grainchanger/issues) - doing it this
way means we can maximise the speed of the function.

| Function.Name | Description                               | Additional.arguments |
| :------------ | :---------------------------------------- | :------------------- |
| prop          | Calculate the proportion of a given class | lc\_class (numeric)  |
| shdi          | Calculate the Shannon diversity           | lc\_class (numeric)  |
| shei          | Calculate the Shannon evenness            | lc\_class (numeric)  |
| range         | Calculate the range of values             |                      |

# Additional utilities

## Create torus

The `create_torus` function takes as input a raster and pads it by a
specified radius, creating the effect of a torus. We developed this
function in order to avoid edge effects when testing methods on
simulated rasters (such as those from
[NLMR](https://ropensci.github.io/NLMR/)).

``` r
torus <- create_torus(cat_ls, 5)

show_landscape(torus, discrete = TRUE)
```

<img src="man/figures/README-torus-1.png" width="100%" />

# Contributing

We welcome contributions to this package. To contribute, submit a [pull
request](https://help.github.com/en/articles/about-pull-requests) making
sure `develop` is the destination branch on the [`grainchanger`
repository](https://github.com/ropensci/grainchanger).

# Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/grainchanger/issues/new/).
  - License: GPL3
  - Get citation information for `grainchanger` in R doing
    `citation(package = 'grainchanger')`
  - Please note that the `grainchanger` project is released with a
    [Contributor Code of
    Conduct](https://github.com/ropensci/grainchanger/blob/master/CODE_OF_CONDUCT.md).
    By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# grainchanger 0.3.2
* skip var_range calculation is correct test on CRAN

# grainchanger 0.3.1
* Add GitHub actions for testing

# grainchanger 0.3.0
* Fix failing test on linux

# grainchanger 0.2.0

* successful review through rOpenSci
* added in functionality to aggregate data to polygon (use parameter `is_grid = FALSE`) 
* added in ability for user to specify aggregation function (`agg_fun`) in `winmove_agg`
* note that many parameter names have changed - check docs for more info: 
    * `fn` is now `win_fun`
    * `g` is now `coarse_dat`
    * `dat` is now `fine_dat`
* `agg_fun` and `win_fun` should be entered as unquoted function names (previously character input)

# grainchanger 0.1.0

* First stable release of grainchanger 

# grainchanger 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
## Test environments
* l Windows Server 2012 R2 x64 (build 9600) (on AppVeyor), R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0
* mac os x 10.12.6 (on travis-ci), R 3.6.0

## R CMD check results

0 errors | 0 warnings | 0 notes