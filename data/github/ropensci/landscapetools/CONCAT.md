
[![Travis build
status](https://travis-ci.org/ropensci/landscapetools.svg?branch=master)](https://travis-ci.org/ropensci/landscapetools)
[![Build
status](https://ci.appveyor.com/api/projects/status/aehfkxfb5r4vjlm9?svg=true)](https://ci.appveyor.com/project/ropensci/landscapetools)
[![codecov](https://codecov.io/gh/ropensci/landscapetools/branch/develop/graph/badge.svg)](https://codecov.io/gh/ropensci/landscapetools)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
status](https://www.r-pkg.org/badges/version/landscapetools)](https://cran.r-project.org/package=landscapetools)
[![](http://cranlogs.r-pkg.org/badges/grand-total/landscapetools)](http://cran.rstudio.com/web/packages/landscapetools/index.html)
[![](https://badges.ropensci.org/188_status.svg)](https://github.com/ropensci/onboarding/issues/188)
[![DOI:10.1111/2041-210X.13076](https://zenodo.org/badge/DOI/10.1111/2041-210X.13076.svg)](https://doi.org/10.1111/2041-210X.13076)

# landscapetools

`landscapetools` provides utility functions for some of the
less-glamorous tasks involved in landscape analysis:

#### Utilities:

  - `util_binarize`: Binarize continuous raster values, if \> 1 breaks
    are given, return a RasterBrick.
  - `util_classify`: Classify a raster into proportions based upon a
    vector of class weightings.
  - `util_merge`: Merge a primary raster with other rasters weighted by
    scaling factors.
  - `util_raster2tibble`, `util_tibble2raster`: Coerce raster\* objects
    to tibbles and vice versa.
  - `util_rescale`: Linearly rescale element values in a raster to a
    range between 0 and 1.
  - `util_writeESRI`: Export raster objects as ESRI asciis (with Windows
    linebreaks).

#### Visualization

  - `show_landscape`: Plot a Raster\* object with the landscapetools
    default theme (as ggplot) or multiple raster (RasterStack, -brick or
    list of raster) side by side as facets.
  - `show_shareplot`: Plot the landscape share in subsequential buffers
    around a/multiple point(s) of interest

#### Themes:

  - `theme_nlm`, `theme_nlm_grey`: Opinionated ggplot2 theme to
    visualize raster (continuous data).
  - `theme_nlm_discrete`, `theme_nlm_grey_discrete`: Opinionated ggplot2
    theme to visualize raster (discrete data).
  - `theme_faceplot`: Opinionated ggplot2 theme to visualize raster in a
    facet wrap.

## Installation

You can install the released version from CRAN with:

``` r
install.packages("landscapetools")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/landscapetools")
```

## Utilities

### Classify

``` r
# Classify the landscape into land uses
classified_landscape <- util_classify(fractal_landscape,
                                      n = 3,
                                      level_names = c("Land Use 1", 
                                                      "Land Use 2",
                                                      "Land Use 3"))

show_landscape(classified_landscape, discrete = TRUE)
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

### Merge

``` r
# Merge all landscapes into one
merged_landscape <- util_merge(fractal_landscape,
                               c(gradient_landscape, random_landscape),
                               scalingfactor = 1)

# Plot an overview
merge_vis <- list(
    "1) Primary" = fractal_landscape,
    "2) Secondary 1" = gradient_landscape,
    "3) Secondary 2" = random_landscape,
    "4) Result" = merged_landscape
)

show_landscape(merge_vis)
#> Warning: Removed 1196 rows containing missing values (geom_raster).
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## See also

In the examples above we make heavy use of the `NLMR` package. Both
packages were developed together until we split them into pure landscape
functionality and utility tools. If you are interested in generating
neutral landscapes via a multitude of available algorithms take a closer
look at the [NLMR](https://github.com/ropensci/NLMR/) package.

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/landscapetools/issues/new/).
  - License: GPL3
  - Get citation information for `landscapetools` in R doing
    `citation(package = 'landscapetools')`
  - We are very open to contributions - if you are interested check
    [Contributing](CONTRIBUTING.md).
      - Please note that this project is released with a [Contributor
        Code of Conduct](CODE_OF_CONDUCT.md). By participating in this
        project you agree to abide by its
terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# landscapetools 0.6.2
- Bugfix in `util_classify`

# landscapetools 0.6.0
- `util_raster2tibble` can now return a wide tibble
- New function `show_shareplot`
- `util_as_integer` now returns integer values from 1:n instead of rounding numeric values

# landscapetools 0.5.0
- new interface for `util_classify`
    - now takes argument n to specify number of classes
    - n argument implemented in C++
- Removed Roboto font and `util_import_roboto`
- Removed `util_plot_grey`
- Renamed:
    - `util_plot` to `show_landscape`
- new function `util_writeESRI` that produces a replica of esris ascii file format

# landscapetools 0.4.0

* minor bug fixes
* util_facetplot now better handles lists of raster
* improved theme_facetplot
* util_classify can now reclassify based on real landscapes, the classification then overwrites the weightings with the proportions from this landscape
* util_classify now has an mask argument, that allows for the classification only outside this mask
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# CONTRIBUTING #

### Please contribute!

We love collaboration.

### Bugs?

* Submit an issue on the Issues page [here](https://github.com/marcosci/landscapetools/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/landscapetools.git`
* Make sure to track progress upstream (i.e., on our version of `landscapetools` at `marcosci/landscapetools`) by doing `git remote add upstream https://github.com/marcosci/landscapetools.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new branch)
* If you alter package functionality at all (e.g., the code itself, not just documentation)
please do write some tests to cover the new functionality.
* Push up to your account
* Submit a pull request to home base at `marcosci/landscapetools`

### Questions? Get in touch: [sciaini.marco@gmail.com](mailto:sciaini.marco@gmail.com)

### Thanks for contributing!
## Submission

New version includes massive dependency trimming and new functionality.

## Test environments

* local Ubuntu Linux 16.04 LTS install, R 3.4.1
* Ubuntu 14.04 (on travis-ci), R 3.4.1
* Windows Server 2012 R2 x64 (build 9600) (on appveyor), R 3.4.2
* Rhub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran
* macOS 10.11 El Capitan, R-release
* macOS 10.9 Mavericks, R-oldrel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Reverse dependencies

There are currently no reverse dependencies.
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/{owner}/{repo}/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/{repo}.git`
* Make sure to track progress upstream (i.e., on our version of `{repo}` at `{owner}/{repo}`) by doing `git remote add upstream https://github.com/{owner}/{repo}.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `{owner}/{repo}`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email?

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
output:
  github_document:
    html_preview: false
---
[![Travis build status](https://travis-ci.org/ropensci/landscapetools.svg?branch=master)](https://travis-ci.org/ropensci/landscapetools)
[![Build status](https://ci.appveyor.com/api/projects/status/aehfkxfb5r4vjlm9?svg=true)](https://ci.appveyor.com/project/ropensci/landscapetools)
[![codecov](https://codecov.io/gh/ropensci/landscapetools/branch/develop/graph/badge.svg)](https://codecov.io/gh/ropensci/landscapetools)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN status](https://www.r-pkg.org/badges/version/landscapetools)](https://cran.r-project.org/package=landscapetools)
[![](http://cranlogs.r-pkg.org/badges/grand-total/landscapetools)](http://cran.rstudio.com/web/packages/landscapetools/index.html)
[![](https://badges.ropensci.org/188_status.svg)](https://github.com/ropensci/onboarding/issues/188)
[![DOI:10.1111/2041-210X.13076](https://zenodo.org/badge/DOI/10.1111/2041-210X.13076.svg)](https://doi.org/10.1111/2041-210X.13076)

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# landscapetools

`landscapetools` provides utility functions for some of the less-glamorous tasks involved
in landscape analysis:

#### Utilities:

- `util_binarize`: Binarize continuous raster values, if > 1 breaks are given, return a RasterBrick.
- `util_classify`: Classify a raster into proportions based upon a vector of class weightings.
- `util_merge`: Merge a primary raster with other rasters weighted by scaling factors.
- `util_raster2tibble`, `util_tibble2raster`: Coerce raster* objects to tibbles and vice versa.
- `util_rescale`: Linearly rescale element values in a raster to a range between 0 and 1.
- `util_writeESRI`: Export raster objects as ESRI asciis (with Windows linebreaks).

#### Visualization

- `show_landscape`: Plot a Raster* object with the landscapetools default theme (as ggplot) or multiple raster (RasterStack, -brick or list of raster) side by side as facets.
- `show_shareplot`: Plot the landscape share in subsequential buffers around a/multiple point(s) of interest

#### Themes:

- `theme_nlm`, `theme_nlm_grey`: Opinionated ggplot2 theme to visualize raster (continuous data).
- `theme_nlm_discrete`, `theme_nlm_grey_discrete`: Opinionated ggplot2 theme to visualize raster (discrete data).
- `theme_faceplot`: Opinionated ggplot2 theme to visualize raster in a facet wrap.

## Installation

You can install the released version from CRAN with:

```r
install.packages("landscapetools")
```

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/landscapetools")
```

## Utilities

```{r load_libraries_hidden, eval=TRUE, echo=FALSE, message=FALSE, results='hide'}
library(landscapetools)
```

### Classify

```{r fig.retina=2, message=FALSE}
# Classify the landscape into land uses
classified_landscape <- util_classify(fractal_landscape,
                                      n = 3,
                                      level_names = c("Land Use 1", 
                                                      "Land Use 2",
                                                      "Land Use 3"))

show_landscape(classified_landscape, discrete = TRUE)
```

### Merge

```{r fig.retina=2, message=FALSE}
# Merge all landscapes into one
merged_landscape <- util_merge(fractal_landscape,
                               c(gradient_landscape, random_landscape),
                               scalingfactor = 1)

# Plot an overview
merge_vis <- list(
    "1) Primary" = fractal_landscape,
    "2) Secondary 1" = gradient_landscape,
    "3) Secondary 2" = random_landscape,
    "4) Result" = merged_landscape
)

show_landscape(merge_vis)
```

## See also

In the examples above we make heavy use of the `NLMR` package.
Both packages were developed together until we split them into pure landscape functionality and utility tools.
If you are interested in generating neutral landscapes via a multitude of available algorithms take a closer look at the [NLMR](https://github.com/ropensci/NLMR/) package.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/landscapetools/issues/new/).
* License: GPL3
* Get citation information for `landscapetools` in R doing `citation(package = 'landscapetools')`
* We are very open to contributions - if you are interested check [Contributing](CONTRIBUTING.md).
    * Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Short walkthrough and overview of landscapetools"
author: "Marco Sciaini"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(landscapetools)
```

*landscapetools* is not a coherent package designed for a specific scientific purpose, it is rather a collection of functions to perform  some of the less-glamorous tasks involved in landscape analysis.

It is basically designed to accompany all the packages in [r-spatialecology](https://github.com/r-spatialecology) and keep them lightweight. Hence, the functionality has a broad spectrum and we try to cover here some things one might miss about *landscapetools*.

# Visualize 

There are a plethora of R packages to visualize spatial data, all of them covering unique aspects and ways to do that (find a short introduction [here](https://docs.ropensci.org/NLMR/articles/articles/visualize_nlms.html)). With [NLMR](https://docs.ropensci.org/NLMR), we needed a way to visualize landscapes without much fuss and also have a way to visualize many of them in a way we found sufficient.

## General raster plotting
```{r fig.retina=2, message=FALSE, warning=FALSE}
# Plot continous landscapes
show_landscape(gradient_landscape)

# Plot continous landscapes 
show_landscape(classified_landscape, discrete = TRUE)

# RasterStack/RasterBrick
show_landscape(raster::stack(gradient_landscape, random_landscape), discrete = TRUE)

# Plot a list of raster (list names become facet text)
show_landscape(list("Gradient landscape" = gradient_landscape,
                    "Random landscape" = random_landscape))

# Plot multiple raster with unique scales
show_landscape(raster::stack(gradient_landscape, random_landscape, classified_landscape), unique_scales = TRUE)
```

# Scaling 
## Binarize

In landscape ecology, many people often work with landscapes that reflect a matrix / habitat context.
If you work with simulated landscale, `util_binarize` is a convienent wrapper to achieve this.
You can define a value in the range of your landscape values and get a binary reflection of it:

```{r fig.retina=2}
# Binarize the landscape into habitat and matrix
binarized_raster <- util_binarize(fractal_landscape, breaks = 0.31415)
show_landscape(binarized_raster, discrete = TRUE)

# You can also provide a vector with thresholds and get a RasterStack with multiple binarized maps
binarized_raster <- util_binarize(fractal_landscape, breaks = c(0.25, 0.5, 0.7))
show_landscape(binarized_raster)
```

## Classify

Complementary to `util_binarize`, `util_classify` classifies a 
landscape with continuous values into *n* discrete classes.
The function is quite the workhorse, so I will spent some more details here
to explain everything:

```{r fig.retina=2}
# Mode 1: Classify landscape into 3 classes based on the Fisher-Jenks algorithm:
mode_1 <- util_classify(fractal_landscape, n = 3)

# Mode 2: Classify landscapes into landscape with exact proportions:
mode_2 <- util_classify(fractal_landscape, weighting = c(0.5, 0.25, 0.25))

# Mode 3: Classify landscapes based on a real dataset (which we first create here)
#         and the distribution of values in this real dataset
mode_3 <- util_classify(gradient_landscape, n = 3)

## Mode 3a: ... now we just have to provide the "real landscape" (mode_3)
mode_3a <- util_classify(fractal_landscape, real_land = mode_3)

## Mode 3b: ... and we can also say that certain values are not important for our classification:
mode_3b <- util_classify(fractal_landscape, real_land = mode_3, mask_val = 1)

landscapes <- list(
'Mode 1'  = mode_1,
'Mode 2'  = mode_2,
'Mode 3'  = mode_3,
'Mode 3a' = mode_3a,
'Mode 3b' = mode_3b
)

show_landscape(landscapes, unique_scales = TRUE, nrow = 1)

# ... you can also name the classes:
classified_raster <- util_classify(fractal_landscape,
                                   n = 3,
                                   level_names = c("Land Use 1",
                                                   "Land Use 2",
                                                   "Land Use 3"))
show_landscape(classified_raster, discrete = TRUE)
```

## Rescale
`util_rescale` l linearly rescale element values in a raster to a range between 0 and 1.

```{r fig.retina=2}
library(raster) 
landscape <- raster(matrix(1:100, 10, 10))
summary(landscape)

scaled_landscape <- util_rescale(landscape)
summary(scaled_landscape)
```

## Merge
`util_merge` most likely makes sense in the context of [NLMR](https://docs.ropensci.org/NLMR). If you merge multiple 
neutral landscapes models, you can create more feasible landscape
patterns for certain questions, or come up with ecotones if you
merge fractal patterns with gradients.

```{r fig.retina=2}
# Merge all maps into one
merg <- util_merge(fractal_landscape, c(gradient_landscape, random_landscape), scalingfactor = 1)

# Plot an overview
merge_vis <- list(
    "1) Primary" = fractal_landscape,
    "2) Secondary 1" = gradient_landscape,
    "3) Secondary 2" = random_landscape,
    "4) Result" = merg
)
show_landscape(merge_vis)
```

# Export
Some propriatery requires that .asc files have the same
line breaks as ESRI ArcMap produces. As we didn't find a 
correct parser in R, we wrote our on:

```{r eval=FALSE}
util_rescale(fractal_landscape, "fractal.asc")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gradient_landscape}
\alias{gradient_landscape}
\title{Example map (planar gradient).}
\format{A raster layer object.}
\source{
Simulated neutral landscape models with R. \url{https://github.com/ropensci/NLMR/}
}
\usage{
gradient_landscape
}
\description{
An example map to show landscapetools functionality
generated with the nlm_planargradient() algorithm.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theme_nlm.R
\name{theme_nlm}
\alias{theme_nlm}
\alias{theme_nlm_discrete}
\alias{theme_nlm_grey}
\alias{theme_nlm_grey_discrete}
\alias{theme_facetplot}
\alias{theme_facetplot_discrete}
\title{theme_nlm}
\usage{
theme_nlm(base_family = NA, base_size = 11.5,
  plot_title_family = base_family, plot_title_size = 18,
  plot_title_face = "bold", plot_title_margin = 10,
  subtitle_family = NA, subtitle_size = 13, subtitle_face = "plain",
  subtitle_margin = 15, strip_text_family = base_family,
  strip_text_size = 12, strip_text_face = "plain",
  strip.background = "grey80", caption_family = NA, caption_size = 9,
  caption_face = "plain", caption_margin = 10,
  axis_text_size = base_size, axis_title_family = base_family,
  axis_title_size = 9, axis_title_face = "plain",
  axis_title_just = "rt", plot_margin = ggplot2::unit(c(0, 0, 0, 0),
  "lines"), grid_col = "#cccccc", grid = TRUE, axis_col = "#cccccc",
  axis = FALSE, ticks = FALSE, legend_title = "Z",
  legend_labels = NULL, legend_text_size = 8, legend_title_size = 10,
  ratio = 1, viridis_scale = "D", ...)

theme_nlm_discrete(base_family = NA, base_size = 11.5,
  plot_title_family = base_family, plot_title_size = 18,
  plot_title_face = "bold", plot_title_margin = 10,
  subtitle_family = NA, subtitle_size = 13, subtitle_face = "plain",
  subtitle_margin = 15, strip_text_family = base_family,
  strip_text_size = 12, strip_text_face = "plain",
  strip.background = "grey80", caption_family = NA, caption_size = 9,
  caption_face = "plain", caption_margin = 10,
  axis_text_size = base_size, axis_title_family = base_family,
  axis_title_size = 9, axis_title_face = "plain",
  axis_title_just = "rt", plot_margin = ggplot2::unit(c(0, 0, 0, 0),
  "lines"), grid_col = "#cccccc", grid = TRUE, axis_col = "#cccccc",
  axis = FALSE, ticks = FALSE, legend_title = "Z",
  legend_labels = NULL, legend_text_size = 8, legend_title_size = 10,
  ratio = 1, viridis_scale = "D", ...)

theme_nlm_grey(base_family = NA, base_size = 11.5,
  plot_title_family = base_family, plot_title_size = 18,
  plot_title_face = "bold", plot_title_margin = 10,
  subtitle_family = NA, subtitle_size = 13, subtitle_face = "plain",
  subtitle_margin = 15, strip_text_family = base_family,
  strip_text_size = 12, strip_text_face = "plain",
  strip.background = "grey80", caption_family = NA, caption_size = 9,
  caption_face = "plain", caption_margin = 10,
  axis_text_size = base_size, axis_title_family = base_family,
  axis_title_size = 9, axis_title_face = "plain",
  axis_title_just = "rt", plot_margin = ggplot2::unit(c(0, 0, 0, 0),
  "lines"), grid_col = "#cccccc", grid = TRUE, axis_col = "#cccccc",
  axis = FALSE, ticks = FALSE, legend_title = "Z",
  legend_labels = NULL, legend_text_size = 8, legend_title_size = 10,
  ratio = 1, ...)

theme_nlm_grey_discrete(base_family = NA, base_size = 11.5,
  plot_title_family = base_family, plot_title_size = 18,
  plot_title_face = "bold", plot_title_margin = 10,
  subtitle_family = NA, subtitle_size = 13, subtitle_face = "plain",
  subtitle_margin = 15, strip_text_family = base_family,
  strip_text_size = 12, strip_text_face = "plain",
  strip.background = "grey80", caption_family = NA, caption_size = 9,
  caption_face = "plain", caption_margin = 10,
  axis_text_size = base_size, axis_title_family = base_family,
  axis_title_size = 9, axis_title_face = "plain",
  axis_title_just = "rt", plot_margin = ggplot2::unit(c(0, 0, 0, 0),
  "lines"), grid_col = "#cccccc", grid = TRUE, axis_col = "#cccccc",
  axis = FALSE, ticks = FALSE, legend_title = "Z",
  legend_labels = NULL, legend_text_size = 8, legend_title_size = 10,
  ratio = 1, ...)

theme_facetplot(base_family = NA, base_size = 11.5,
  plot_title_family = base_family, plot_title_size = 18,
  plot_title_face = "bold", plot_title_margin = 10,
  subtitle_family = NA, subtitle_size = 13, subtitle_face = "plain",
  subtitle_margin = 15, strip.background = "grey80",
  caption_family = NA, caption_size = 9, caption_face = "plain",
  caption_margin = 10, ratio = 1, viridis_scale = "D", ...)

theme_facetplot_discrete(base_family = NA, base_size = 11.5,
  plot_title_family = base_family, plot_title_size = 18,
  plot_title_face = "bold", plot_title_margin = 10,
  subtitle_family = NA, subtitle_size = 13, subtitle_face = "plain",
  subtitle_margin = 15, strip.background = "grey80",
  caption_family = NA, caption_size = 9, caption_face = "plain",
  caption_margin = 10, ratio = 1, viridis_scale = "D", ...)
}
\arguments{
\item{base_family}{base font family size}

\item{base_size}{base font size}

\item{plot_title_family}{plot title family}

\item{plot_title_size}{plot title size}

\item{plot_title_face}{plot title face}

\item{plot_title_margin}{plot title ggplot2::margin}

\item{subtitle_family}{plot subtitle family}

\item{subtitle_size}{plot subtitle size}

\item{subtitle_face}{plot subtitle face}

\item{subtitle_margin}{plot subtitle ggplot2::margin bottom (single numeric value)}

\item{strip_text_family}{facet facet label font family}

\item{strip_text_size}{facet label font family, face and size}

\item{strip_text_face}{facet facet label font face}

\item{strip.background}{strip background}

\item{caption_family}{plot caption family}

\item{caption_size}{plot caption size}

\item{caption_face}{plot caption face}

\item{caption_margin}{plot caption ggplot2::margin}

\item{axis_text_size}{axis text size}

\item{axis_title_family}{axis title family}

\item{axis_title_size}{axis title size}

\item{axis_title_face}{axis title face}

\item{axis_title_just}{axis title justification}

\item{plot_margin}{plot ggplot2::margin (specify with `ggplot2::margin``)}

\item{grid_col}{grid color}

\item{grid}{grid TRUE/FALSE}

\item{axis_col}{axis color}

\item{axis}{axis TRUE/FALSE}

\item{ticks}{ticks TRUE/FALSE}

\item{legend_title}{Title of the legend (default \code{"Z"})}

\item{legend_labels}{Labels for the legend ticks, if
used with \code{\link{show_landscape}} they are automatically derived.}

\item{legend_text_size}{legend text size, default 8}

\item{legend_title_size}{legend text size, default 10}

\item{ratio}{ratio for tiles (default 1, if your raster is not a square the ratio should
be \code{raster::nrow(x) / raster::ncol(x)})}

\item{viridis_scale}{Five options are available: "viridis - magma" (= "A"),
"viridis - inferno" (= "B"),
"viridis - plasma" (= "C"),
"viridis - viridis" (= "D",  the default option),
"viridis - cividis" (= "E")}

\item{...}{optional arguments to ggplot2::theme}
}
\description{
Opinionated ggplot2 theme to visualize NLM raster.
}
\details{
A focused theme to visualize raster data that sets a lot of defaults for the
\code{ggplot2::theme}.

The functions are setup in such a way that you can customize your own one by
just wrapping the call and changing the parameters.
The theme itself is heavily influenced by hrbrmstr and his package
hrbrthemes (\url{https://github.com/hrbrmstr/hrbrthemes/}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fractal_landscape}
\alias{fractal_landscape}
\title{Example map (fractional brownian motion).}
\format{A raster layer object.}
\source{
Simulated neutral landscape models with R. \url{https://github.com/ropensci/NLMR/}
}
\usage{
fractal_landscape
}
\description{
An example map to show landscapetools functionality
generated with the nlm_fbm() algorithm.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_raster2tibble.R
\name{util_raster2tibble}
\alias{util_raster2tibble}
\title{Converts raster data into tibble}
\usage{
util_raster2tibble(x, format = "long")

util_raster2tibble(x, format = "long")
}
\arguments{
\item{x}{Raster* object}

\item{format}{Either \emph{"long"} (default) or \emph{"wide"} output for the resulting tibble}
}
\value{
a tibble
}
\description{
Writes spatial raster values into tibble and adds coordinates.
}
\details{
You will loose any resolution, extent or reference system.
The output is raw tiles.
}
\examples{
maptib <- util_raster2tibble(fractal_landscape)
\dontrun{
library(ggplot2)
ggplot(maptib, aes(x,y)) +
    coord_fixed() +
    geom_raster(aes(fill = z))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_writeESRI.R
\name{util_writeESRI}
\alias{util_writeESRI}
\alias{util_writeESRI.RasterLayer}
\title{util_writeESRI}
\usage{
util_writeESRI(x, filepath)

\method{util_writeESRI}{RasterLayer}(x, filepath)
}
\arguments{
\item{x}{Raster* object}

\item{filepath}{path where to write the raster to file}
}
\description{
Export raster objects as ESRI ascii files.
}
\details{
\code{raster::writeRaster} or \code{SDMTools::write.asc} both
export files that are recognised by most GIS software, nevertheless
they both have UNIX linebreaks.
Some proprietary software (like SPIP for example) require an exact 1:1
replica of the output of ESRI's ArcMap, which as a Windows software
has no carriage returns at the end of each line.
\code{util_writeESRI} should therefore only be used if you need this,
otherwise \code{raster::writeRaster} is the better fit for exporting
raster data in R.
}
\examples{
\dontrun{
util_writeESRI(gradient_landscape, "gradient_landscape.asc")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_calcboundaries.R
\name{util_calc_boundaries}
\alias{util_calc_boundaries}
\title{util_calc_boundaries}
\usage{
util_calc_boundaries(x, cumulative_proportions)
}
\arguments{
\item{x}{vector of data values.}

\item{cumulative_proportions}{Vector of class cumulative proportions, as generated by \code{w2cp}.}
}
\value{
Numerical vector with boundaries for matrix classification
}
\description{
Determine upper class boundaries for classification of a vector with values ranging 0-1 based upon an
vector of cumulative proportions.
}
\examples{
x <- matrix(runif(100,0,1),10,10)
y <- util_w2cp(c(0.5, 0.25, 0.25)) #cumulative proportion
util_calc_boundaries(x,y)

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_merge.R
\name{util_merge}
\alias{util_merge}
\alias{util_merge.RasterLayer}
\title{util_merge}
\usage{
util_merge(primary_nlm, secondary_nlm, scalingfactor = 1, rescale)

\method{util_merge}{RasterLayer}(primary_nlm, secondary_nlm,
  scalingfactor = 1, rescale = TRUE)
}
\arguments{
\item{primary_nlm}{Primary \code{Raster* object}}

\item{secondary_nlm}{A list or stack of \code{Raster* object}s that are merged with the primary \code{Raster* object}}

\item{scalingfactor}{Weight for the secondary \code{Raster* objects}}

\item{rescale}{If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
Rectangular matrix with values ranging from 0-1
}
\description{
Merge a primary raster with other rasters weighted by scaling factors.
}
\examples{
x <- util_merge(gradient_landscape, random_landscape)
show_landscape(x)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_shareplot.R
\name{show_shareplot}
\alias{show_shareplot}
\title{show_shareplot}
\usage{
show_shareplot(landscape, points, buffer_width, max_width,
  return_df = FALSE)

show_shareplot(landscape, points, buffer_width, max_width,
  return_df = FALSE)
}
\arguments{
\item{landscape}{Raster* object}

\item{points}{Point(s) represented by a two-column matrix or data.frame; SpatialPoints*; SpatialPolygons*; SpatialLines; Extent; a numeric vector representing cell numbers; or sf* POINT object}

\item{buffer_width}{Buffer width in which landscape share is measured}

\item{max_width}{Max distance to which buffer_width is summed up; the x axis in the plot}

\item{return_df}{Logical value indicating if a tibble with the underlying data should be returned}
}
\value{
ggplot2 Object
}
\description{
Plot the landscape share in subsequential buffers around a/multiple point(s) of interest
}
\examples{
# create single point
new_point = matrix(c(75,75), ncol = 2)

# show landscape and point of interest
show_landscape(classified_landscape, discrete = TRUE) +
ggplot2::geom_point(data = data.frame(x = new_point[,1], y = new_point[,2]),
                    ggplot2::aes(x = x, y = y),
                    col = "grey", size = 3)

# show single point share
show_shareplot(classified_landscape, new_point, 10, 50)

# show multiple points share
new_points = matrix(c(75, 110, 75, 30), ncol = 2)
show_shareplot(classified_landscape, new_points, 10, 50)

# get data frame with results back
result <- show_shareplot(classified_landscape, new_points, 10, 50, return_df = TRUE)
result$share_df

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_as_integer.R
\name{util_as_integer}
\alias{util_as_integer}
\alias{util_as_integer.RasterLayer}
\title{util_as_integer}
\usage{
util_as_integer(x)

\method{util_as_integer}{RasterLayer}(x)
}
\arguments{
\item{x}{raster}
}
\value{
RasterLayer
}
\description{
Coerces raster values to integers
}
\details{
Coerces raster values to integers, which is sometimes needed if you want further
methods that rely on integer values.
}
\examples{
# Mode 1
util_as_integer(fractal_landscape)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_binarize.R
\name{util_binarize}
\alias{util_binarize}
\alias{util_binarize.RasterLayer}
\title{Binarize continuous raster values}
\usage{
util_binarize(x, breaks)

\method{util_binarize}{RasterLayer}(x, breaks)
}
\arguments{
\item{x}{Raster* object}

\item{breaks}{Vector with one or more break percentages}
}
\value{
RasterLayer / RasterBrick
}
\description{
Classify continuous raster values into binary map cells based upon given
break(\code{s}).
}
\details{
Breaks are considered to be habitat percentages (\code{p}). If more than
one percentage is given multiple layers are written in the same brick.
}
\examples{
breaks <- c(0.3, 0.5)
binary_maps <- util_binarize(gradient_landscape, breaks)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{classified_landscape}
\alias{classified_landscape}
\title{Example map (factor).}
\format{A raster layer object.}
\source{
Simulated neutral landscape models with R. \url{https://github.com/ropensci/NLMR/}
}
\usage{
classified_landscape
}
\description{
An example map to show landscapetools functionality
generated with the nlm_random() algorithm with factorial values.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_tibble2raster.R
\name{util_tibble2raster}
\alias{util_tibble2raster}
\title{Converts tibble data into a raster}
\usage{
util_tibble2raster(x)

util_tibble2raster(x)
}
\arguments{
\item{x}{a tibble}
}
\value{
Raster* object
}
\description{
Writes spatial tibble values into a raster.
}
\details{
Writes tiles with coordinates from a tibble into a raster.
Resolution is set to 1 and the extent will be c(0, max(x), 0, max(y)).

You can directly convert back the result from 'util_raster2tibble()' without
problems. If you have altered the coordinates or otherwise played with the
data, be careful while using this function.
}
\examples{
maptib <- util_raster2tibble(random_landscape)
mapras <- util_tibble2raster(maptib)
all.equal(random_landscape, mapras)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_landscape.R
\name{show_landscape}
\alias{show_landscape}
\alias{show_landscape.RasterLayer}
\alias{show_landscape.list}
\alias{show_landscape.RasterStack}
\alias{show_landscape.RasterBrick}
\title{show_landscape}
\usage{
show_landscape(x, xlab, ylab, discrete, unique_scales, n_col, n_row, ...)

\method{show_landscape}{RasterLayer}(x, xlab = "Easting",
  ylab = "Northing", discrete = FALSE, ...)

\method{show_landscape}{list}(x, xlab = "Easting", ylab = "Northing",
  discrete = FALSE, unique_scales = FALSE, n_col = NULL,
  n_row = NULL, ...)

\method{show_landscape}{RasterStack}(x, xlab = "Easting",
  ylab = "Northing", discrete = FALSE, unique_scales = FALSE,
  n_col = NULL, n_row = NULL, ...)

\method{show_landscape}{RasterBrick}(x, xlab = "Easting",
  ylab = "Northing", discrete = FALSE, unique_scales = FALSE,
  n_col = NULL, n_row = NULL, ...)
}
\arguments{
\item{x}{Raster* object}

\item{xlab}{x axis label, default "Easting"}

\item{ylab}{y axis label, default "Northing"}

\item{discrete}{If TRUE, the function plots a raster with
a discrete legend.}

\item{unique_scales}{If TRUE and multiple raster are to be visualized, each facet can have a unique color scale for its fill}

\item{n_col}{If multiple rasters are to be visualized, n_col controls the number of columns for the facet}

\item{n_row}{If multiple rasters are to be visualized, n_row controls the number of rows for the facet}

\item{...}{Arguments for  \code{\link{theme_nlm}}}
}
\value{
ggplot2 Object
}
\description{
Plot a Raster* object with the NLMR default theme (as ggplot).
}
\examples{
\dontrun{
x <- gradient_landscape

# classify
y <- util_classify(gradient_landscape,
                   n = 3,
                   level_names = c("Land Use 1", "Land Use 2", "Land Use 3"))

show_landscape(x)
show_landscape(y, discrete = TRUE)

show_landscape(list(gradient_landscape, random_landscape))
show_landscape(raster::stack(gradient_landscape, random_landscape))

show_landscape(list(gradient_landscape, y), unique_scales = TRUE)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_rescale.R
\name{util_rescale}
\alias{util_rescale}
\title{util_rescale}
\usage{
util_rescale(x)

util_rescale(x)
}
\arguments{
\item{x}{Raster* object}
}
\value{
Raster* object with values ranging from 0-1
}
\description{
Linearly rescale element values in a raster to a range between 0 and 1.
}
\details{
Rasters generated by \code{nlm_} functions are scaled between 0 and 1 as default, this option can be set to \code{FALSE} if needed.
}
\examples{
unscaled_landscape <- gradient_landscape + fractal_landscape
util_rescale(unscaled_landscape)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landscapetools.R
\docType{package}
\name{landscapetools-package}
\alias{landscapetools}
\alias{landscapetools-package}
\title{landscapetools}
\description{
\emph{landscapetools} provides utility functions to work with landscape data
(raster* Objects).
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/landscapetools/}
  \item Report bugs at \url{https://github.com/ropensci/landscapetools/issues}
}

}
\author{
\strong{Maintainer}: Marco Sciaini \email{sciaini.marco@gmail.com} (0000-0002-3042-5435)

Authors:
\itemize{
  \item Matthias Fritsch \email{matthias.fritsch@forst.uni-goettingen.de}
  \item Maximillian H.K. Hesselbarth \email{maximilian.hesselbarth@uni-goettingen.de} (0000-0003-1125-9918)
  \item Jakub Nowosad \email{nowosad.jakub@gmail.com} (0000-0002-1057-3721)
}

Other contributors:
\itemize{
  \item Laura Graham (Laura reviewed the package for rOpenSci, see 
                  https://github.com/ropensci/onboarding/issues/188) [reviewer]
  \item Jeffrey Hollister (Jeffrey reviewed the package for rOpenSci, see 
                  https://github.com/ropensci/onboarding/issues/188) [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_w2cp.R
\name{util_w2cp}
\alias{util_w2cp}
\title{util_w2cp}
\usage{
util_w2cp(weighting)
}
\arguments{
\item{weighting}{A list of numeric values}
}
\value{
Rectangular matrix with values ranging from 0-1
}
\description{
Convert a list of category  weighting  into a 1D array of cumulative proportions.
}
\examples{
util_w2cp(c(0.2, 0.4, 0.6, 0.9))


}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_classify.R
\name{util_classify}
\alias{util_classify}
\alias{util_classify.RasterLayer}
\title{util_classify}
\usage{
util_classify(x, n, weighting, level_names, real_land, mask_val)

\method{util_classify}{RasterLayer}(x, n = NULL, weighting = NULL,
  level_names = NULL, real_land = NULL, mask_val = NULL)
}
\arguments{
\item{x}{raster}

\item{n}{Number of classes}

\item{weighting}{Vector of numeric values that are considered to be habitat percentages (see details)}

\item{level_names}{Vector of names for the factor levels.}

\item{real_land}{Raster with real landscape (see details)}

\item{mask_val}{Value to mask (refers to real_land)}
}
\value{
RasterLayer
}
\description{
Classify continuous landscapes into landscapes with discrete classes
}
\details{
Mode 1: Calculate the optimum breakpoints using Jenks natural
breaks optimization, the number of classes is determined with \code{n}.
The Jenks optimization seeks to minimize the variance within categories,
while maximizing the variance between categories.

Mode 2: The number of elements in the weighting vector determines the number of classes
in the resulting matrix. The classes start with the value 1.
If non-numerical levels are required, the user can specify a vector to turn the
numerical factors into other data types, for example into character strings (i.e. class labels).
If the numerical vector of weightings does not sum up to 1, the sum of the
weightings is divided by the number of elements in the weightings vector and this is then used for the classificat#'     .

Mode 3: For a given 'real' landscape the number of classes and the weightings are
extracted and used to classify the given landscape (any given weighting parameter is
overwritten in this case!). If an optional mask value is given the corresponding
class from the 'real' landscape is cut from the landscape beforehand.
}
\examples{
\dontrun{
# Mode 1
util_classify(fractal_landscape,
              n = 3,
              level_names = c("Land Use 1", "Land Use 2", "Land Use 3"))

# Mode 2
util_classify(fractal_landscape,
              weighting = c(0.5, 0.25, 0.25),
              level_names = c("Land Use 1", "Land Use 2", "Land Use 3"))

# Mode 3
real_land <- util_classify(gradient_landscape,
              n = 3,
              level_names = c("Land Use 1", "Land Use 2", "Land Use 3"))

fractal_landscape_real <- util_classify(fractal_landscape, real_land = real_land)
fractal_landscape_mask <- util_classify(fractal_landscape, real_land = real_land, mask_val = 1)

landscapes <- list(
'1 nlm' = fractal_landscape,
'2 real' = real_land,
'3 result' = fractal_landscape_real,
'4 result with mask' = fractal_landscape_mask
)

show_landscape(landscapes, unique_scales = TRUE, nrow = 1)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{random_landscape}
\alias{random_landscape}
\title{Example map (random).}
\format{A raster layer object.}
\source{
Simulated neutral landscape models with R. \url{https://github.com/ropensci/NLMR/}
}
\usage{
random_landscape
}
\description{
An example map to show landscapetools functionality
generated with the nlm_random() algorithm.
}
\keyword{datasets}
