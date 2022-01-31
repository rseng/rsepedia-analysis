
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

0 errors | 0 warnings | 0 notes---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# grainchanger <img src="man/figures/logo.png" align="right" width="150" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/grainchanger/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/grainchanger/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/grainchanger/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/grainchanger?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/grainchanger)](https://cran.r-project.org/package=grainchanger)
[![CranLogs](https://cranlogs.r-pkg.org/badges/grainchanger)](https://cran.r-project.org/package=grainchanger)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

The `grainchanger` package provides functionality for data aggregation to a coarser resolution via moving-window or direct methods. 

As landscape ecologists and macroecologists, we often need to aggregate data in order to harmonise datasets. In doing so, we often lose a lot of information about the spatial structure and environmental heterogeneity of data measured at finer resolution. An example of this is when the response data (e.g. species' atlas data) are available at a coarser resolution to the predictor data (e.g. land-use data). We developed this method and R package in order to overcome some of these issues. 

For more information on the background to and motivation for the development of this method, see [Graham *et al.* 2019 in *Methods in Ecology and Evolution*](https://doi.org/10.1111/2041-210X.13177). 

# Package overview

The primary functions of the `grainchanger` package are those which facilitate moving-window (`winmove_agg`) and direct (`nomove_agg`) data aggregation. These functions aggregate fine-grain data (`fine_dat`) to a coarse-grain (`coarse_dat`) using a function specified by the user (`agg_fun`). The moving-window method takes in an additional function (`win_fun`) which smooths the fine-grain data prior to aggregation. 

The moving-window smoothing function is also available in the package (`winmove`), as well as several built-in functions, and an additional utility function for use with simulated landscapes (`create_torus`).

The `winmove` function acts as a convenient wrapper to `raster::focalWeight` and `raster::focal` which takes advantage of optimised functions built into the `grainchanger` package. 

# Installation

```
# Install release version from CRAN
install.packages("grainchanger")

# Install development version from GitHub
devtools::install_github("ropensci/grainchanger")
```

# Examples

## Moving-window data aggregation

The below example shows the moving-window data aggregation in action. It aggregates a categorical raster (`fine_dat`) to a grid using Shannon evenness (specified by `win_fun`) as the function calculated within a square moving window of 5 units. The value returned is the mean (specified by `agg_fun`) of the smoothed value for each cell of `coarse_dat`. This value is included as a column on the grid `sf` object. 

```{r mwda_example, fig.show = "hold"}
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

## Direct data aggregation

The below example shows the direct data aggregation in action. It aggregates a continuous raster to a raster with a coarser resolution using the range as the function calculated for each cell of the larger grid. The resulting output is a raster of the coarser resolution. `var_range` is an inbuilt function in the `grainchanger` package.

```{r dda_example, fig.show = "hold"}
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

# Functions

There are a number of inbuilt functions in the grainchanger package, with their usage outlined below. While it is possible to use user-defined functions within both `winmove_agg` and `nomove_agg`, we welcome suggestions for additional functions. Please [add as an issue](https://github.com/ropensci/grainchanger/issues) - doing it this way means we can maximise the speed of the function. 

```{r functions, echo = FALSE}
function_overview <- data.frame(
  `Function Name` = c("prop", "shdi", "shei", "range"),
  `Description` = c("Calculate the proportion of a given class", 
                    "Calculate the Shannon diversity", 
                    "Calculate the Shannon evenness", 
                    "Calculate the range of values"),
  `Additional arguments` = c("lc_class (numeric)", 
                             "lc_class (numeric)",
                             "lc_class (numeric)",
                             "")
)

knitr::kable(function_overview)
```

# Additional utilities

## Create torus

The `create_torus` function takes as input a raster and pads it by a specified radius, creating the effect of a torus. We developed this function in order to avoid edge effects when testing methods on simulated rasters (such as those from [NLMR](https://ropensci.github.io/NLMR/)). 

```{r torus}
torus <- create_torus(cat_ls, 5)

show_landscape(torus, discrete = TRUE)
```

# Contributing

We welcome contributions to this package. To contribute, submit a [pull request](https://help.github.com/en/articles/about-pull-requests) making sure `develop` is the destination branch on the [`grainchanger` repository](https://github.com/ropensci/grainchanger).

# Meta

* Please [report any issues or bugs](https://github.com/ropensci/grainchanger/issues/new/).
* License: GPL3
* Get citation information for `grainchanger` in R doing `citation(package = 'grainchanger')`
* Please note that the `grainchanger` project is released with a [Contributor Code of Conduct](https://github.com/ropensci/grainchanger/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Background & Motivation"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The `grainchanger` package provides functionality for data aggregation to a grid via moving-window or direct methods.

As landscape ecologists and macroecologists, we often need to aggregate data in order to harmonise datasets. In doing so, we often lose a lot of information about the spatial structure and environmental heterogeneity of data measured at finer resolution. An example of this is when the response data (e.g. species' atlas data) are available at a coarser resolution to the predictor data (e.g. land-use data). We developed this method and R package in order to overcome some of these issues. 

For more information on the background to and motivation for the development of this method, see [Graham *et al.* 2019 in *Methods in Ecology and Evolution*](https://doi.org/10.1111/2041-210X.13177). 

## Moving-window data aggregation

The moving-window data aggregation (MWDA) method smooths an input raster using a specified function within a moving window of a specified size and shape prior to aggregation. This acts as a convenient wrapper for the `focalWeight()` and `focal()` functions in the `raster` package. Additionally, we have aimed to write efficient functions for some oft-used metrics within landscape ecology for use within the moving window.

```{r}
knitr::include_graphics("../man/figures/mwda_schematic.png")
```


The above is a graphical representation of the MWDA method. In calculating the MWDA measure, three aspects of scale are considered. Predictor grain is the characteristic spatial scale of a predictor variable, that is, the resolution of the environmental data; scale‐of‐effect determines the appropriate scale of the relationship between predictor and response, for example, an ecological neighbourhood; response grain is the grain of the unit into which you are predicting, that is, the resolution of the response variable (represented by the black lines). Note that the colour scale is unitless. Yellow cells represent ‘high’ values and dark blue cells ‘low’ values. Panel 1 shows a close up of one of the response grain cells in panel 2, whereas panel 2 shows all response grain cells for the study region. Panel 3 shows the study region after aggregation. From [*Graham et al. 2019*](https://doi.org/10.1111/2041-210X.13177).

## Direct data aggregation 

The direct method simply aggregates to the coarse data using the specified function. For example, say we want to calculate proportion of forest at the municipality or county level. ---
title: "Built-in functions"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

There are a number of built-in functions in the `grainchanger` package, with their usage outlined below. While it is possible to use user-defined functions within `winmove_agg`, `nomove_agg`, and `winmove`, we welcome suggestions for additional functions. Please [add as an issue](https://github.com/laurajanegraham/grainchanger/issues) - doing it this way means we can maximise the speed of the function. 

All functions can also be used on their own, either on an object of class `winmove` or `numeric`. 

When functions are used within `winmove_agg`, `winmove`, or directly on an object of class `winmove`, they are calculated relative to within a moving window. 

When functions are used within `nomove_agg` all cells of `fine_dat` within a given cell of `coarse_dat` are aggregated using the function. 

# Current functions

```{r functions, echo = FALSE}
function_overview <- data.frame(
  `Function Name` = c("prop", "shdi", "shei", "var_range"),
  `Description` = c("Calculate the proportion of a given class", 
                    "Calculate the Shannon diversity", 
                    "Calculate the Shannon evenness", 
                    "Calculate the size of the range of values"),
  `Additional arguments` = c("lc_class (numeric)", 
                             "lc_class (numeric)",
                             "lc_class (numeric)",
                             "")
)

knitr::kable(function_overview)
```

# Shannon diversity and evenness

Shannon diversity is calculated as $$SHDI = -\sum_{i = 1}^m p_i lnp_i$$ where $p_i$ is the proportion of a given class $i$ of a total $m$ classes. 

Shannon evenness is calculated as $$SHEI = \frac{S}{ln(m)}$$

# Additional functions

We plan to add other useful functions to this small set of built-in functions, such as relevant metrics from [FRAGSTATS](https://www.umass.edu/landeco/research/fragstats/documents/fragstats.help.4.2.pdf). 

We also welcome suggestions for additional functions. Please [add as an issue](https://github.com/laurajanegraham/grainchanger/issues) - doing it this way means we can maximise the speed of the function. ---
title: "Using Grainchanger"
author: "Laura J. Graham"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using Grainchanger}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>"
)
```

# Package overview

The primary functions of the `grainchanger` package are those which facilitate moving-window (`winmove_agg`) and direct (`nomove_agg`) data aggregation. These functions aggregate fine-grain data (`fine_dat`) to a coarse-grain (`coarse_dat`) using a function specified by the user (`agg_fun`). The moving-window method takes in an additional function (`win_fun`) which smooths the fine-grain data prior to aggregation. 

The moving-window smoothing function is also available in the package (`winmove`), as well as several built-in functions, and an additional utility function for use with simulated landscapes (`create_torus`).

The `winmove` function acts as a convenient wrapper to `raster::focalWeight` and `raster::focal` which takes advantage of optimised functions built into the `grainchanger` package. 

# Data aggregation & `sf`

The functions in `grainchanger` have been written to be compatible with the [`sf`](https://r-spatial.github.io/sf/) package. The online textbook [Geocomputation with R](https://geocompr.robinlovelace.net/) provides an introduction to spatial data analysis with the `sf` package. 

## Moving-window data aggregation

In this example, we show how to create a coarse grain grid over the fine-grain data using `st_make_grid` from the `sf` package. We then aggregate the package data `cat_ls` to this grid using moving-window data aggregation and plot using `ggplot2`. 

```{r}
library(grainchanger)
library(sf)
library(ggplot2)

coarse_dat <- cat_ls %>% 
  # get the bounding box
  st_bbox() %>% 
  # turn into an sfc object
  st_as_sfc() %>% 
  # negative buffer 
  st_buffer(-4) %>% 
  # make a square grid
  st_make_grid(cellsize = 19) %>% 
  # turn into sf object
  st_sf()

# we can plot this grid on top of the fine data
landscapetools::show_landscape(cat_ls) + 
  geom_sf(data = coarse_dat, alpha = 0.5)


coarse_dat$shdi_3 <- winmove_agg(coarse_dat = coarse_dat, 
                                 fine_dat = cat_ls,
                                 d = 3,
                                 type = "rectangle", 
                                 win_fun = shdi, 
                                 agg_fun = mean,
                                 is_grid = FALSE,
                                 lc_class = 1:4)

ggplot(coarse_dat, aes(fill = shdi_3)) + 
  geom_sf() + 
  theme_bw()
```

Note that when creating the grid, we made it 4 units smaller (the size of the moving window) than the fine-resolution data. The reason for this is that we avoid edge effects created when the moving window goes beyond the extent of the data. A warning is thrown if `fine_dat` is the same size as (or smaller than) `coarse_dat`. The below is an example of this, using `g_sf` from the package. 

```{r, warning = TRUE}
g_sf$shei_4 <- winmove_agg(coarse_dat = g_sf, 
                                 fine_dat = cat_ls,
                                 d = 4,
                                 type = "rectangle", 
                                 win_fun = shei, 
                                 agg_fun = mean,
                                 is_grid = FALSE,
                                 lc_class = 1:4)
```


## Direct data aggregation

In this example we show how to read in a shapefile as an `sf` object, apply direct data aggregation, and plot using `ggplot2`. `cont_ls` is a fine-resolution RasterLayer provided with the package. 

```{r}
library(sf)
library(ggplot2)

# coarse_dat <- st_read("your_file.shp")
coarse_dat <- st_read(system.file("shape/poly_sf.shp", package="grainchanger"))

coarse_dat$var_range <- nomove_agg(coarse_dat = coarse_dat,
                                   fine_dat = cont_ls,
                                   agg_fun = var_range,
                                   is_grid = FALSE)

ggplot(coarse_dat, aes(fill = var_range)) + 
  geom_sf() + 
  theme_bw()
```

# Functions

There are a number of inbuilt functions in the grainchanger package, with their usage outlined below. While it is possible to use user-defined functions within both `winmove_agg` and `nomove_agg`, we welcome suggestions for additional functions. Please [add as an issue](https://github.com/ropensci/grainchanger/issues) - doing it this way means we can maximise the speed of the function. 

```{r functions, echo = FALSE}
function_overview <- data.frame(
  `Function Name` = c("prop", "shdi", "shei", "range"),
  `Description` = c("Calculate the proportion of a given class", 
                    "Calculate the Shannon diversity", 
                    "Calculate the Shannon evenness", 
                    "Calculate the range of values"),
  `Additional arguments` = c("lc_class (numeric)", 
                             "lc_class (numeric)",
                             "lc_class (numeric)",
                             "")
)

knitr::kable(function_overview)
```

# Additional utilities

## Create torus

The `create_torus` function takes as input a square or rectangular landscape and pads it by a specified radius, creating the effect of a torus. We developed this function in order to avoid edge effects when testing methods on simulated landscapes (such as those from [NLMR](https://ropensci.github.io/NLMR/)). 

```{r torus}
torus <- create_torus(cat_ls, 5)

landscapetools::show_landscape(torus)
```

## Running aggregation in parallel

The `winmove_agg` and `nomove_agg` functions have been designed to work with the [`future`](https://github.com/HenrikBengtsson/future) package. This means that with just one additional line of code, the functions can be run in parallel while taking advantage of the resources available to the user. 

```
library(future)
plan(multisession)
```
See the [`future`](https://github.com/HenrikBengtsson/future) documentation for more information. 
---
title: "User-defined functions"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{User-defined functions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Within the `winmove`, `winmove_agg` and `nomove_agg` functions, it is possible to use user-defined functions for both `win_fun` and `agg_fun` arguments. 

> **WARNING** User-defined functions can be slower within the `grainchanger` functions because they have not been optimised. This is likely to be of particular issue with large datasets.

# User-defined `win_fun` example

Any user-defined `win_fun` should follow the rules of the `fun` argument in `raster::focal`:

> The function fun should take multiple numbers, and return a single number. For example mean, modal, min or max. It should also accept a na.rm argument (or ignore it, e.g. as one of the 'dots' arguments. For example, length will fail, but function(x, ...){na.omit(length(x))} works.

In this example, we define a function which counts the number of cells of a given class within a moving window. 

```{r}
library(grainchanger)
library(landscapetools)

num_cells <- function(x, lc_class, ...) {
  return(sum(x == lc_class))
}
d <- winmove(cat_ls, 4, "rectangle", num_cells, lc_class = 2)
show_landscape(d) 
```

This can also be used within `winmove_agg`

```{r}
library(ggplot2)
g_sf$num_cells <- winmove_agg(g_sf, cat_ls, 4, "rectangle", num_cells, lc_class = 2)

ggplot(g_sf, aes(fill = num_cells)) + 
  scale_fill_viridis_c() + 
  geom_sf() + 
  theme_bw()
```

# User-defined `agg_fun`

In this example, we define a function which calculates the number of land cover classes within each coarse grain cell. 

```{r}
num_classes <- function(x, ...) {
  length(unique(x))
}

g_sf$num_classes <- nomove_agg(g_sf, cat_ls, num_classes)

ggplot(g_sf, aes(fill = as.factor(num_classes))) +
  scale_fill_viridis_d("num_classes") + 
  geom_sf() + 
  theme_bw()
```

We can also define functions which work on continuous landscapes. For example, below we calculate the coefficient of variation for each coarse cell. 

```{r}
cv <- function(x) {
  sd(x) / mean(x)
}

poly_sf$cv <- nomove_agg(poly_sf, cont_ls, cv)

ggplot(poly_sf, aes(fill = cv)) +
  scale_fill_viridis_c() + 
  geom_sf() + 
  theme_bw()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_functions.R
\name{var_range}
\alias{var_range}
\alias{var_range.winmove}
\alias{var_range.numeric}
\title{Size of range of values}
\usage{
var_range(x, ...)

\method{var_range}{winmove}(x, d, type, na.rm = TRUE, ...)

\method{var_range}{numeric}(x, na.rm = TRUE, ...)
}
\arguments{
\item{x}{RasterLayer. The data over which to calculate the range size}

\item{...}{further arguments passed to or from other methods}

\item{d}{numeric. If \code{type=circle}, the radius of the circle (in units of the
CRS). If \code{type=rectangle} the dimension of the rectangle (one or two numbers)}

\item{type}{character. The shape of the moving window}

\item{na.rm}{logical. indicates whether \code{NA} values should be stripped before the
computation proceeds. \code{na.rm = TRUE} is the default}
}
\value{
If \code{class(x) == "winmove"}, a smoothed raster with the size of the range of values calculated within the specified
  moving window
  
  If \code{class(x) == "numeric"}, a single value representing the size of the range of values in \code{x}
}
\description{
Calculates the difference between the maximum and minimum value
}
\examples{

# load required data
data(cat_ls)
data(cont_ls)

# convert data to object of class winmove
cat_ls <- new("winmove", cat_ls)

# aggregate using a rectangular window with dimensions c(2,3)
d <- range(cont_ls, d = c(2,3), type = "rectangle")

# convert data to object of class numeric
cont_ls <- raster::values(cont_ls)
d <- range(cont_ls)

}
\keyword{focal}
\keyword{range}
\keyword{spatial}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_functions.R
\name{mean}
\alias{mean}
\alias{mean.winmove}
\title{Arithmetic mean}
\usage{
mean(x, ...)

\method{mean}{winmove}(x, d, type, ...)
}
\arguments{
\item{x}{RasterLayer. The data over which to calculate the mean value within a moving
window}

\item{...}{further arguments passed to or from other methods}

\item{d}{numeric. If \code{type=circle}, the radius of the circle (in units of the
CRS). If \code{type=rectangle} the dimension of the rectangle (one or two numbers)}

\item{type}{character. The shape of the moving window}
}
\value{
RasterLayer. A smoothed raster with the mean calculated within the specified
  moving window
}
\description{
An extension to \code{mean} for objects of class \code{winmove}
}
\examples{

# load required data
data(cont_ls)

# convert data to object of class winmove
cont_ls <- new("winmove", cont_ls)

# aggregate using a circular window with radius 3
d <- mean(cont_ls, d = 3, type = "circle")
}
\keyword{focal}
\keyword{mean}
\keyword{spatial}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cat_ls}
\alias{cat_ls}
\title{Example categorical raster (fine_dat)}
\format{
A raster layer object.
}
\source{
Sciaini M, Fritsch M, Scherer C, Simpkins CE. NLMR and landscapetools: An integrated environment
    for simulating and modifying neutral landscape models in R. Methods in Ecology and Evolution. 2018;
    00:1-9. https://doi.org/10.1111/2041-210X.13076

Marco Sciaini and Matthias Fritsch (2018). landscapetools: Landscape Utility Toolbox. R package version 0.4.0.
    https://CRAN.R-project.org/package=landscapetools
}
\usage{
cat_ls
}
\description{
An example map to show functionality on categorical surfaces.
}
\details{
Generated with \code{nlm_mpd()} from \code{NLMR} and classified with \code{util_classify()} from \code{landscapetools}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/winmove.R
\docType{class}
\name{winmove-class}
\alias{winmove-class}
\title{An S4 class for use with winmove functions (extends RasterLayer)}
\description{
An S4 class for use with winmove functions (extends RasterLayer). Objects
  will need to be set to this class in order to be used with the inbuilt \code{winmove}
  functions (e.g. \code{mean}, \code{prop}, \code{var_range}, \code{shdi}, \code{shei})
}
\section{Slots}{

Slots for RasterLayer and RasterBrick objects
	\describe{
    \item{\code{title}:}{Character} 
    \item{\code{file}:}{Object of class \code{".RasterFile"} }
    \item{\code{data}:}{Object of class \code{".SingleLayerData"} or \code{".MultipleLayerData"}}
    \item{\code{history}:}{To record processing history, not yet in use }
    \item{\code{legend}:}{Object of class \code{.RasterLegend}, Default legend. Should store preferences for plotting. Not yet implemented except that it stores the color table of images, if available}
    \item{\code{extent}:}{Object of \code{\link[raster]{Extent-class}} }
    \item{\code{ncols}:}{Integer} 
    \item{\code{nrows}:}{Integer} 
    \item{\code{crs}:}{Object of class \code{"CRS"}, i.e. the coordinate reference system. In Spatial* objects this slot is called 'proj4string' }
  }

}

\examples{
# load required data
data(cat_ls)

# set \code{cat_ls} to object of class winmove
new("winmove", cat_ls)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/winmove.R
\name{winmove}
\alias{winmove}
\title{Create moving window surface}
\usage{
winmove(fine_dat, d, type = c("circle", "rectangle"), win_fun, ...)
}
\arguments{
\item{fine_dat}{The raster dataset on which to calculate the moving window function}

\item{d}{numeric. If \code{type=circle}, the radius of the circle (in units of the
CRS). If \code{type=rectangle} the dimension of the rectangle (one or two numbers).}

\item{type}{The shape of the moving window}

\item{win_fun}{function. The function to apply. If not choosing one of the inbuilt
grainchanger functions, the function should take multiple numbers, and return a
single number. For example \code{mean}, \code{modal}, \code{min} or \code{max}. It should also accept a \code{na.rm}
argument (or ignore it, e.g. as one of the 'dots' arguments. For example, length will
fail, but \code{function(x, ...){na.omit(length(x))}} works. See Details}

\item{...}{further arguments passed to or from other methods}
}
\value{
RasterLayer. A smoothed raster with the moving window values calculated
}
\description{
Smooth a raster surface using a moving window with a given function, radius and shape.
}
\details{
\code{grainchanger} has several built-in functions. Functions currently
  included are: \itemize{ \item \code{wm_shei} - Shannon evenness, requires the
  additional argument \code{lc_class} (vector or scalar) \item \code{wm_prop} -
  Proportion, requires the additional argument \code{lc_class} (scalar) \item
  \code{wm_classes} - Unique number of classes in a categorical landscape \item
  \code{var_range} - Range (max - min) }
}
\examples{
# load required data
data(cat_ls)
data(cont_ls)

# calculate the moving window mean
d <- winmove(cont_ls, 5, "rectangle", mean)

# calculate the moving window Shannon evenness
d <- winmove(cat_ls, 5, "rectangle", shei, lc_class = 1:4)
}
\keyword{focal}
\keyword{spatial}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_functions.R
\name{prop}
\alias{prop}
\alias{prop.winmove}
\alias{prop.numeric}
\title{Calculate proportion of a given value}
\usage{
prop(x, lc_class, ...)

\method{prop}{winmove}(x, lc_class, d, type, ...)

\method{prop}{numeric}(x, lc_class, ...)
}
\arguments{
\item{x}{numeric, winmove. The data over which to calculate the proportion}

\item{lc_class}{numeric. The class value to calculate the proportion of}

\item{...}{further arguments passed to or from other methods}

\item{d}{numeric. If \code{type=circle}, the radius of the circle (in units of the
CRS). If \code{type=rectangle} the dimension of the rectangle (one or two numbers)}

\item{type}{character. The shape of the moving window}
}
\value{
If \code{class(x) == "winmove"}, a smoothed raster with the proportion of
  cells of the given class calculated within the specified moving window

  If \code{class(x) == "numeric"}, a single value representing the proportion of values
  of a given class in \code{x}
}
\description{
Calculate the proportion of a given value present within a raster. Useful for
calculating land-cover or soil type proportions. Should be used with a categorical
raster
}
\examples{

# load required data
data(cat_ls)

# convert data to object of class winmove
cat_ls <- new("winmove", cat_ls)

# aggregate using a rectangular window with dimension 5 for class 3
d <- prop(cat_ls, d = 5, type = "rectangle", lc_class = 3)

# convert data to object of class numeric
cat_ls <- raster::values(cat_ls)
d <- prop(cat_ls, lc_class = 2)
}
\keyword{focal}
\keyword{mean}
\keyword{spatial}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{poly_sf}
\alias{poly_sf}
\title{Example polygon (coarse_dat)}
\format{
An sf object.
}
\usage{
poly_sf
}
\description{
An example non-gridded coarse data to show functionality when aggregating using an sf object.
}
\details{
Generated with \code{sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(cont_ls)), cellsize = 13, square = FALSE)}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nomove_agg.R
\name{nomove_agg}
\alias{nomove_agg}
\title{Direct data aggregation}
\usage{
nomove_agg(coarse_dat, fine_dat, agg_fun, is_grid = TRUE, quiet = FALSE, ...)
}
\arguments{
\item{coarse_dat}{sf, Raster* or Spatial* object. The coarse grain data
(response data) across which to calculate the aggregated function}

\item{fine_dat}{Raster* object. Raster* object. The fine grain data
(predictor / covariate data) to aggregate}

\item{agg_fun}{function The function to apply. The function fun should take
multiple numbers, and return a single number. For example mean, modal, min
or max. It should also accept a na.rm argument (or ignore it, e.g. as one
of the 'dots' arguments. For example, length will fail, but function(x,
...){na.omit(length(x))} works. See Details}

\item{is_grid}{logical. Use \code{TRUE} (default) if \code{g} contains only
rectangular cells (i.e. a grid). If \code{g} is any other polygon file,
this should be set to false}

\item{quiet}{logical. If \code{FALSE} (default) and \code{is_grid == TRUE}
the user gets a warning that the aggregation assumes all cells are
rectangular}

\item{...}{further arguments passed to or from other methods}
}
\value{
Raster (if input is Raster) or numeric vector (if input is sp or sf
  object) containing values calculated for each coarser cell
}
\description{
Calculate the value for a given function for each cell in a larger resolution
grid.
}
\details{
\code{grainchanger} has several built-in functions. Functions
  currently included are: 
  \itemize{ 
     \item \code{shdi} - Shannon diversity, requires the additional argument \code{lc_class} (vector or scalar) 
     \item \code{shei} - Shannon evenness, requires the additional argument \code{lc_class} (vector or scalar) 
     \item \code{prop} - Proportion, requires the additional argument \code{lc_class} (scalar)
     \item \code{var_range} - Range (max - min) 
     }
     
 Note that \code{nomove_agg} can be run in parallel using \code{plan(multiprocess)} from the \code{future} package.
}
\examples{
# load required data
data(g_sf)
data(cont_ls)
data(cat_ls)

# aggregate using mean
d <- nomove_agg(g_sf, cont_ls, mean)

# aggregate using Shannon evenness
d <- nomove_agg(g_sf, cont_ls, shei, lc_class = 1:4)
}
\keyword{aggregate}
\keyword{spatial}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/winmove_agg.R
\name{winmove_agg}
\alias{winmove_agg}
\title{Moving-window data aggregation}
\usage{
winmove_agg(
  coarse_dat,
  fine_dat,
  d,
  type = c("circle", "rectangle"),
  win_fun,
  agg_fun = mean,
  is_grid = TRUE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{coarse_dat}{sf, Raster* or Spatial* object. The coarse grain data
(response data) across which to calculate the aggregated moving window
function}

\item{fine_dat}{Raster* object. The fine grain data (predictor / covariate
data) to aggregate}

\item{d}{numeric. If \code{type=circle}, the radius of the circle (in units
of the CRS). If \code{type=rectangle} the dimension of the rectangle (one
or two numbers).}

\item{type}{character. The shape of the moving window}

\item{win_fun}{character. The function to apply to the moving window. The
function \code{win_fun} should take multiple numbers, and return a single number. For
example \code{mean}, \code{modal}, \code{min} or \code{max}. It should also accept a \code{na.rm} argument (or
ignore it, e.g. as one of the 'dots' arguments. For example, \code{length} will
fail, but \code{function(x, ...){na.omit(length(x))}} works. See Details}

\item{agg_fun}{character. The function by which to aggregate. By default this
is set to \code{mean}}

\item{is_grid}{logical. Use \code{TRUE} (default) if \code{g} contains only
rectangular cells (i.e. a grid). If \code{g} is any other polygon file,
this should be set to false}

\item{quiet}{logical. If \code{FALSE} (default) and \code{is_grid == TRUE}
the user gets a warning that the aggregation assumes all cells are
rectangular}

\item{...}{further arguments passed to or from other methods}
}
\value{
Numeric vector containing moving window values calculated for each
  grid cell
}
\description{
Calculate the mean moving window value for a given radius, shape and function
for each cell in a larger resolution grid.
}
\details{
\code{grainchanger} has several built-in functions. Functions
  currently included are: 
  \itemize{ 
     \item \code{shdi} - Shannon diversity, requires the additional argument \code{lc_class} (vector or scalar) 
     \item \code{shei} - Shannon evenness, requires the additional argument \code{lc_class} (vector or scalar) 
     \item \code{prop} - Proportion, requires the additional argument \code{lc_class} (scalar)
     \item \code{var_range} - Range (max - min) 
     }
     
 Note that \code{winmove_agg} can be run in parallel using \code{plan(multiprocess)} from the \code{future} package.
}
\examples{
\dontrun{
# load required data
data(g_sf)
data(cont_ls)
data(cat_ls)

# aggregate using mean
d <- winmove_agg(g_sf, cont_ls, 5, "rectangle", mean)

# aggregate using Shannon evenness
d <- winmove_agg(g_sf, cat_ls, 5, "rectangle", shei, lc_class = 1:4)
}

}
\keyword{aggregate}
\keyword{focal}
\keyword{spatial}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grainchanger-package.R
\docType{package}
\name{grainchanger-package}
\alias{grainchanger}
\alias{grainchanger-package}
\title{grainchanger: Moving-Window and Direct Data Aggregation}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Data aggregation via moving window or direct methods. Aggregate a 
    fine-resolution raster to a grid. The moving window method smooths the surface 
    using a specified function within a moving window of a specified size and shape 
    prior to aggregation. The direct method simply aggregates to the grid using the 
    specified function.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/grainchanger/}
  \item \url{https://github.com/ropensci/grainchanger}
  \item Report bugs at \url{https://github.com/ropensci/grainchanger/issues}
}

}
\author{
\strong{Maintainer}: Laura Graham \email{LauraJaneEGraham@gmail.com} (\href{https://orcid.org/0000-0002-3611-7281}{ORCID})

Other contributors:
\itemize{
  \item Felix Eigenbrod \email{f.eigenbrod@soton.ac.uk} (Input on initial conceptual development) [contributor]
  \item Marco Sciaini \email{sciaini.marco@gmail.com} (Input on package development and structure) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_functions.R
\name{diversity-metrics}
\alias{diversity-metrics}
\alias{shdi}
\alias{shdi.winmove}
\alias{shdi.numeric}
\alias{shei}
\alias{shei.winmove}
\alias{shei.numeric}
\title{Diversity metrics}
\usage{
\method{shdi}{winmove}(x, lc_class, d, type, ...)

\method{shdi}{numeric}(x, lc_class, ...)

\method{shei}{winmove}(x, lc_class, d, type, ...)

\method{shei}{numeric}(x, lc_class, ...)
}
\arguments{
\item{x}{numeric, winmove. The data over which to calculate the diversity metrics}

\item{lc_class}{numeric. The class values to include in the diversity metric
calculation}

\item{d}{numeric. If \code{type=circle}, the radius of the circle (in units of the
CRS). If \code{type=rectangle} the dimension of the rectangle (one or two numbers)}

\item{type}{character. The shape of the moving window}

\item{...}{further arguments passed to or from other methods}
}
\value{
If \code{class(x) == "winmove"}, a smoothed raster with the diversity
  metric calculated within the specified moving window

  If \code{class(x) == "numeric"}, a single value representing the diversity metric in
  \code{x}
}
\description{
A range of functions to calculate well known landcover diversity metrics
}
\details{
Currently provided diversity metrics are Shannon diversity and Shannon
  evenness. Open a new issue (https://github.com/laurajanegraham/grainchanger/issues)
  to request additional diversity metrics.
}
\examples{
# load required data
data(cat_ls)

# convert data to object of class winmove
cat_ls <- new("winmove", cat_ls)

# calculate Shannon diversity in a rectangular window of dimension 5
d <- shdi(cat_ls, d = 5, type = "rectangle", lc_class = 1:4)

# convert data to object of class numeric
cat_ls <- raster::values(cat_ls)

# calculate Shannon evenness
d <- shei(cat_ls, lc_class = 1:4)
}
\references{
McGarigal, K. and Marks, B.J., 1995. FRAGSTATS: spatial pattern analysis
  program for quantifying landscape structure. \emph{Gen. Tech. Rep. PNW-GTR-351. Portland,
  OR: US Department of Agriculture, Forest Service, Pacific Northwest Research Station.
  122 p, 351.}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_torus.R
\name{create_torus}
\alias{create_torus}
\title{Pad a raster by a specified radius}
\usage{
create_torus(dat, dpad)
}
\arguments{
\item{dat}{The raster dataset to pad}

\item{dpad}{The amount by which to pad the raster (in the same units as the
raster)}
}
\value{
raster. Original raster padded by r cells with torus effect (see
  Details)
}
\description{
This function pads a raster by a specified number of cells, creating the
effect of a torus. This function is intended for use on simulated landscapes,
in order to avoid edge effects
}
\details{
A torus is an infinite surface where the top joins the bottom, and
  the left side meets the right side. See https://en.wikipedia.org/wiki/Torus
  for a full mathematical description.

  In this function, the torus effect is achieved by adding the specified
  number of rows of the top of the raster to the bottom (and vice versa) and
  the specified number of rows of the right of the raster to the left (and
  vice versa)
}
\examples{
data(cat_ls)
d <- create_torus(dat = cat_ls, dpad = 5)
}
\keyword{raster}
\keyword{torus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cont_ls}
\alias{cont_ls}
\title{Example continuous raster (fine_dat)}
\format{
A raster layer object.
}
\source{
Sciaini M, Fritsch M, Scherer C, Simpkins CE. NLMR and landscapetools: An integrated environment
    for simulating and modifying neutral landscape models in R. Methods in Ecology and Evolution. 2018;
    00:1-9. https://doi.org/10.1111/2041-210X.13076
}
\usage{
cont_ls
}
\description{
An example map to show functionality on continuous surfaces.
}
\details{
Generated with \code{nlm_mpd()} from \code{NLMR}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{g_sf}
\alias{g_sf}
\title{Example grid (coarse_dat)}
\format{
An sf object.
}
\source{
Sciaini M, Fritsch M, Scherer C, Simpkins CE. NLMR and landscapetools: An integrated environment
    for simulating and modifying neutral landscape models in R. Methods in Ecology and Evolution. 2018;
    00:1-9. https://doi.org/10.1111/2041-210X.13076
}
\usage{
g_sf
}
\description{
An example grid to show functionality when aggregating using an sf object.
}
\details{
Generated with \code{nlm_mpd()} and converted to sf.
}
\keyword{datasets}
