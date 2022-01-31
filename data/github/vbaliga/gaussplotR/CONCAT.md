
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gaussplotR <img src='man/figures/logo.png' align="right" height="138.5" />

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build
status](https://github.com/vbaliga/gaussplotR/workflows/R-CMD-check/badge.svg)](https://github.com/vbaliga/gaussplotR/actions)
[![Codecov test
coverage](https://codecov.io/gh/vbaliga/gaussplotR/graph/badge.svg)](https://codecov.io/gh/vbaliga/gaussplotR?branch=master)  
[![status](https://joss.theoj.org/papers/10.21105/joss.03074/status.svg)](https://joss.theoj.org/papers/10.21105/joss.03074)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4041073.svg)](https://doi.org/10.5281/zenodo.4041073)
[![CRAN
status](https://www.r-pkg.org/badges/version/gaussplotR)](https://CRAN.R-project.org/package=gaussplotR)
<!-- badges: end -->

`gaussplotR` provides functions to fit two-dimensional Gaussian
functions, predict values from such functions, and produce plots of
predicted data.

## Installation

You can install `gaussplotR` from CRAN via:

``` r
install.packages("gaussplotR")
```

Or to get the latest (developmental) version through GitHub, use:

``` r
devtools::install_github("vbaliga/gaussplotR")
```

## Example

The function `fit_gaussian_2D()` is the workhorse of `gaussplotR`. It
uses `stats::nls()` to find the best-fitting parameters of a 2D-Gaussian
fit to supplied data based on one of three formula choices. The function
`autofit_gaussian_2D()` can be used to automatically figure out the best
formula choice and arrive at the best-fitting parameters.

The `predict_gaussian_2D()` function can then be used to predict values
from the Gaussian over a supplied grid of X- and Y-values (generated
here via `expand.grid()`). This is useful if the original data is
relatively sparse and interpolation of values is desired.

Plotting can then be achieved via `ggplot_gaussian_2D()`, but note that
the `data.frame` created by `predict_gaussian_2D()` can be supplied to
other plotting frameworks such as `lattice::levelplot()`. A 3D plot can
also be produced via `rgl_gaussian_2D()` (not shown here).

``` r
library(gaussplotR)

## Load the sample data set
data(gaussplot_sample_data)

## The raw data we'd like to use are in columns 1:3
samp_dat <-
  gaussplot_sample_data[,1:3]


#### Example 1: Unconstrained elliptical ####
## This fits an unconstrained elliptical by default
gauss_fit_ue <-
  fit_gaussian_2D(samp_dat)

## Generate a grid of X- and Y- values on which to predict
grid <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -1, to = 4, by = 0.1))

## Predict the values using predict_gaussian_2D
gauss_data_ue <-
  predict_gaussian_2D(
    fit_object = gauss_fit_ue,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
library(ggplot2); library(metR)
#> Warning: package 'ggplot2' was built under R version 4.0.5
#> Warning: package 'metR' was built under R version 4.0.5
ggplot_gaussian_2D(gauss_data_ue)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
## And another example plot via lattice::levelplot()
library(lattice)
lattice::levelplot(
  predicted_values ~ X_values * Y_values,
  data = gauss_data_ue,
  col.regions = colorRampPalette(
    c("white", "blue")
    )(100),
  asp = 1
)
```

<img src="man/figures/README-example-2.png" width="100%" />

``` r
#### Example 2: Constrained elliptical_log ####
## This fits a constrained elliptical, as in Priebe et al. 2003
gauss_fit_cel <-
  fit_gaussian_2D(
    samp_dat,
    method = "elliptical_log",
    constrain_orientation = -1
  )

## Generate a grid of x- and y- values on which to predict
grid <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -1, to = 4, by = 0.1))

## Predict the values using predict_gaussian_2D
gauss_data_cel <-
  predict_gaussian_2D(
    fit_object = gauss_fit_cel,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
ggplot_gaussian_2D(gauss_data_cel)
```

<img src="man/figures/README-example-3.png" width="100%" />

Should you be interested in having `gaussplotR` try to automatically
determine the best choice of `method` for `fit_gaussian_2D()`, the
`autofit_gaussian_2D()` function can come in handy. The default is to
select the `method` that produces a fit with the lowest `rmse`, but
other choices include `rss` and `AIC`.

``` r
## Use autofit_gaussian_2D() to automatically decide the best 
## model to use
gauss_auto <-
  autofit_gaussian_2D(
    samp_dat,
    comparison_method = "rmse", 
    simplify = TRUE
    )

## The output has the same components as `fit_gaussian_2D()` 
## but for the automatically-selected best-fitting method only:
summary(gauss_auto)
#> Model coefficients
#>   A_o   Amp theta X_peak Y_peak    a    b
#>  0.83 32.25  3.58  -2.64   2.02 0.91 0.96
#> Model error stats
#>     rss rmse deviance AIC
#>  156.23 2.08   156.23 171
#> Fitting methods
#>          method       amplitude     orientation 
#>    "elliptical" "unconstrained" "unconstrained"
```

## Contributing and/or raising Issues

Feedback on bugs, improvements, and/or feature requests are all welcome.
Please see the Issues templates on GitHub to make a bug fix request or
feature request.

To contribute code via a pull request, please consult the Contributing
Guide first.

## Citation

Baliga, VB. 2021. gaussplotR: Fit, predict, and plot 2D-Gaussians in R.
Journal of Open Source Software, 6(60), 3074.
<https://doi.org/10.21105/joss.03074>

## License

GPL (&gt;= 3) + file LICENSE

🐢
# gaussplotR 0.2.6
* Added R2 and adjusted R2 computations. To get these, predictions from 
Gaussians are simply regressed against original response data. 

# gaussplotR 0.2.5
* Citation information updated  
* Ready for updating CRAN

# gaussplotR 0.2.4
* Post-peer review at the Journal of Open Source Software.
* `print()` and `summary()` on objects that are output from `fit_gaussian_2D()`
now round parameter estimates to two decimal places.
* R version requirement downgraded to R >= 3.3.0 


# gaussplotR 0.2.2
* Newest functions enhance the automation of comparing among various Gaussian 
fit models.  
* The `autofit_gaussian_2D()` function can be used to find the best-fitting 
model for a given data set.  
* The `compare_gaussian_fits()` function compares models via criteria such as 
rmse or rss.  
* A `characterize_gaussian_fits()` analyzes the orientation and partial 
correlations of Gaussian data. Features include computation of partial 
correlations between response variables and independent and diagonally-tuned 
predictions, along with Z-difference scoring.  

# gaussplotR 0.2.0

* Added a `fit_gaussian_2D()` function to apply any of several available methods
to fit a 2D-Gaussian to supplied data
* `predict_gaussian_2D()` has been enhanced to allow any method to be used in
`fit_gaussian_2D()` and then fed into the predict function

# gaussplotR 0.1.6

* Added a `compute_gaussian_volume()` function to compute the volume under a 
given 2D-Gaussian

# gaussplotR 0.1.4

* Miscellaneous formatting fixes for CRAN checks. No changes to functions.

# gaussplotR 0.1.3

* Added a `NEWS.md` file to track changes to the package.
* Package now seems to be ready for submission to CRAN
## Test environments
* local R installation, R 4.0.2
* windows-latest (release) on GitHub Actions
* macOS-latest (release) on GitHub Actions
* ubuntu-20.04 (release) on GitHub Actions
* ubuntu-20.04 (devel) on GitHub Actions

## R CMD check results

0 errors | 0 warnings | 0 notes

---
title: 'gaussplotR: Fit, Predict and Plot 2D-Gaussians in R'
authors:
- affiliation: 1
  name: Vikram B. Baliga
  orcid: 0000-0002-9367-8974
date: "03 February 2021"
bibliography: paper.bib
tags:
- R
- 2D-Gaussian
- Gaussian fit
- Gaussian orientation
- nonlinear least squares
affiliations:
- index: 1
  name: Department of Zoology, University of British Columbia, Vancouver, British
    Colombia, Canada V6T 1Z4
---

# Summary

Should the need to model the relationship between bivariate data and a response
variable arise, two-dimensional (2D) Gaussian models are often the most
appropriate choice. For example, @Priebe2003 characterized motion-sensitive
neurons in the brains of macaques by fitting 2D-Gaussian functions to neurons'
response rates as spatial and temporal frequencies of visual stimuli were varied.
The width and orientation of these fitted 2D-Gaussian surfaces
provides insight on whether a neuron is "tuned" to particular spatial or
temporal domains. Two-dimensional Gaussians are also used in other scientific
disciplines such as physics [@Wu1998; @Kravtsov2004], materials sciences
[@Riekel1999], and image processing [@Hanuman2013; @Ketenci2013], particularly
in medical imaging [@Wu2019;
@Qadir2020].

Fitting 2D-Gaussian models to data is not always a straightforward process, as
finding appropriate values for the model's parameters relies on complex
procedures such as non-linear least-squares. `gaussplotR` is an R package that
is designed to fit 2D-Gaussian surfaces to data. Should a user supply bivariate
data (i.e., x-values and y-values) along with a univariate response variable,
functions within `gaussplotR` will allow for the automatic fitting of a
2D-Gaussian model to the data. Fitting the model then enables the user to
characterize various properties of the Gaussian surface (e.g., computing the
total volume under the surface). Further, new data can be predicted from models
fit via `gaussplotR`, which in combination with the package's plotting
functions, can enable smoother-looking plots from relatively sparse input data.
In principle, tools within `gaussplotR` have broad applicability to a variety of
scientific disciplines.


# Statement of Need

At the time of writing, we know of no other packages in the `R` ecosystem that 
automatically handle the fitting of 2D-Gaussians to supplied data. The `R` 
package `imagefx` [@imagefx] does offer the capability to predict
data from a 2D-Gaussian model, but only if the parameters of the model are known 
*a priori*. Further, although base `R` functions such as `stats::nls()` provide 
the capability to determine the non-linear least-squares estimates of the 
parameters for a non-linear model, the burden of determining the formula for a 
2D-Gaussian falls upon the user.  

To counter these issues, `gaussplotR` provides users with the capability to fit
2D-Gaussian models using one of three possible formulas, along with the ability
to apply constraints to the amplitude and/or orientation of the fitted Gaussian,
if desired. Coupled with the ability to characterize various properties of the
fitted model, along with plotting functions (as the name of the package
implies), `gaussplotR` is intended to be a feature-rich package for users
interested in 2D-Gaussian modeling. These capabilities are briefly explained
in the next section; vignettes supplied in the package delve into even further
detail.


# Overview and getting started

A series of vignettes that provides detailed guidance are available on
[gaussplotR's GitHub page](https://vbaliga.github.io/gaussplotR/).

The function `fit_gaussian_2D()` is the workhorse of `gaussplotR`. It uses
`stats::nls()` to find the best-fitting parameters of a 2D-Gaussian fit to
supplied data based on one of three formula choices. Each of these formula
choices is designed for a specific use case. The most generic method (and the
default) is `method = "elliptical"`. This allows the fitted 2D-Gaussian to take
an ellipsoid shape, and this will likely be the best option for most use cases.
A slightly-altered method to fit an ellipsoid 2D-Gaussian is available in
`method = "elliptical_log"`. This method follows @Priebe2003
and is geared towards use with log2-transformed data. A third option is `method
= "circular"`. This produces a very simple 2D-Gaussian that is constrained to
have to have a roughly circular shape (i.e. spread in X- and Y- are roughly
equal). Rather than place the burden on the user to determine formula choice,
the function `autofit_gaussian_2D()` can be used to automatically figure out the
best formula choice and arrive at the best-fitting parameters.

In some cases, the researcher may be interested in characterizing the
orientation of the fitted 2D-Gaussian and comparing it to theoretical
predictions. For example, studies of visual neuroscience often describe the
properties of individual motion-sensitive neurons based on whether they are
"speed-tuned" or whether they show independence from the speed of visual
stimuli. Assessing such properties can be done via fitting a 2D-Gaussian to the
response rate of a neuron for a grid of investigated spatial (X-axis) and
temporal frequencies (Y-axis). Should the orientation of the fitted 2D-Gaussian
lie along the diagonal of the plot, the neuron can be classified as
"speed-tuned". The function `characterize_gaussian_fits()` allows for such
analysis within `gaussplotR`. Following methods used in studies of visual
neuroscience [@Levitt1994; @Priebe2003; @Winship2006], the orientation and
partial correlations of 2D-Gaussian data are analyzed. Features include
computation of partial correlations between response variables and independent
and diagonally-tuned predictions, along with Z-difference scoring.

The `predict_gaussian_2D()` function can be used to predict values from the
fitted 2D-Gaussian over a supplied grid of X- and Y-values (usually generated
via `expand.grid()`). This is useful if the original data are relatively sparse
and interpolation of values is desired, e.g. to attain smoother-looking contours
in plots.

Plotting can then be achieved via `ggplot_gaussian_2D()`, but note that the 
`data.frame` created by `predict_gaussian_2D()` can be supplied to other 
plotting frameworks such as `lattice::levelplot()`. A 3D plot can also be 
produced via `rgl_gaussian_2D()`.

`gaussplotR` was designed for broad applicability; there are many disciplines
in which a 2D-Gaussian surface would be a useful model for describing a response
to a bivariate set of inputs. Functions in `gaussplotR` are being used in an 
in-prep article to determine the extent of spatiotemporal tuning of 
motion-sensitive neurons in hummingbirds and other avian species.


# Acknowledgements

We thank Douglas R. Wylie, Douglas L. Altshuler, and Graham Smyth for help in 
working with 2D-Gaussian data.

# References
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
*  Look at the GitHub Actions build status before and after making changes.
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

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
and the ropensci contributing guide available at
https://raw.githubusercontent.com/ropensci/dotgithubfiles/blob/master/dotgithub/CONTRIBUTING.md
---
name: 'Bug fix request'
about: 'Notify us of a bug so we can fix it'
title: "[bug_fix_request]"
labels: bug
assignees: ''

---

**Which function(s) is failing?**
A clear and concise description of what the problem is.

**Can you provide a reproducible example of the failure?**
Data on public repositories are preferred with steps to reproduce the failure. If that is not an option, please contact the author directly.

**What error message(s) are you receiving, if any?**

**Expected behavior or describe the solution you'd like**
A clear and concise description of what you expected to happen or what you would like to happen.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: 'Feature Request: General'
about: 'Let us know if you have an idea to improve `gaussplotR`'
title: "[feature_request_general]"
labels: enhancement
assignees: ''

---

**Please describe the nature of your feature request.**
Is the request to improve the tool(s) currently within `gaussplotR` or to
provide a tool that adds to the package's capabilities?

**Please explain the feature you would like to use in `gaussplotR`.**
Briefly describe the goal of the feature and/or which current aspects of `gaussplotR` relate to your request.

**Can you provide an example file to help explain your request?**
Data on public repositories are preferred. If that is not an option, please contact the author directly.

**Are you interested in contributing a function and/or code yourself?**
If so, please feel free to make a pull request for us to review.

**Additional context** 
Add any other context or screenshots about the feature request here.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r opts, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%",
  dpi = 300
)
```

# gaussplotR <img src='man/figures/logo.png' align="right" height="138.5" />

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build status](https://github.com/vbaliga/gaussplotR/workflows/R-CMD-check/badge.svg)](https://github.com/vbaliga/gaussplotR/actions)
[![Codecov test coverage](https://codecov.io/gh/vbaliga/gaussplotR/graph/badge.svg)](https://codecov.io/gh/vbaliga/gaussplotR?branch=master)  
[![status](https://joss.theoj.org/papers/10.21105/joss.03074/status.svg)](https://joss.theoj.org/papers/10.21105/joss.03074)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4041073.svg)](https://doi.org/10.5281/zenodo.4041073)
[![CRAN status](https://www.r-pkg.org/badges/version/gaussplotR)](https://CRAN.R-project.org/package=gaussplotR)
<!-- badges: end -->

`gaussplotR` provides functions to fit two-dimensional Gaussian functions,
predict values from such functions, and produce plots of predicted data.

## Installation

You can install `gaussplotR` from CRAN via:

``` {r install_cran, eval = FALSE}
install.packages("gaussplotR")
```

Or to get the latest (developmental) version through GitHub, use:
  
``` {r install_github, eval = FALSE}
devtools::install_github("vbaliga/gaussplotR")
```


## Example

The function `fit_gaussian_2D()` is the workhorse of `gaussplotR`. It uses
`stats::nls()` to find the best-fitting parameters of a 2D-Gaussian fit to
supplied data based on one of three formula choices. The function
`autofit_gaussian_2D()` can be used to automatically figure out the best formula
choice and arrive at the best-fitting parameters.

The `predict_gaussian_2D()` function can then be used to predict values from
the Gaussian over a supplied grid of X- and Y-values (generated here via 
`expand.grid()`). This is useful if the original data is relatively sparse and
interpolation of values is desired.

Plotting can then be achieved via `ggplot_gaussian_2D()`, but note that the 
`data.frame` created by `predict_gaussian_2D()` can be supplied to other 
plotting frameworks such as `lattice::levelplot()`. A 3D plot can also be 
produced via `rgl_gaussian_2D()` (not shown here).

```{r example}
library(gaussplotR)

## Load the sample data set
data(gaussplot_sample_data)

## The raw data we'd like to use are in columns 1:3
samp_dat <-
  gaussplot_sample_data[,1:3]


#### Example 1: Unconstrained elliptical ####
## This fits an unconstrained elliptical by default
gauss_fit_ue <-
  fit_gaussian_2D(samp_dat)

## Generate a grid of X- and Y- values on which to predict
grid <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -1, to = 4, by = 0.1))

## Predict the values using predict_gaussian_2D
gauss_data_ue <-
  predict_gaussian_2D(
    fit_object = gauss_fit_ue,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
library(ggplot2); library(metR)
ggplot_gaussian_2D(gauss_data_ue)

## And another example plot via lattice::levelplot()
library(lattice)
lattice::levelplot(
  predicted_values ~ X_values * Y_values,
  data = gauss_data_ue,
  col.regions = colorRampPalette(
    c("white", "blue")
    )(100),
  asp = 1
)

#### Example 2: Constrained elliptical_log ####
## This fits a constrained elliptical, as in Priebe et al. 2003
gauss_fit_cel <-
  fit_gaussian_2D(
    samp_dat,
    method = "elliptical_log",
    constrain_orientation = -1
  )

## Generate a grid of x- and y- values on which to predict
grid <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -1, to = 4, by = 0.1))

## Predict the values using predict_gaussian_2D
gauss_data_cel <-
  predict_gaussian_2D(
    fit_object = gauss_fit_cel,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
ggplot_gaussian_2D(gauss_data_cel)

```

Should you be interested in having `gaussplotR` try to automatically determine
the best choice of `method` for `fit_gaussian_2D()`, the `autofit_gaussian_2D()`
function can come in handy. The default is to select the `method` that 
produces a fit with the lowest `rmse`, but other choices include `rss` and 
`AIC`.

```{r autofit}
## Use autofit_gaussian_2D() to automatically decide the best 
## model to use
gauss_auto <-
  autofit_gaussian_2D(
    samp_dat,
    comparison_method = "rmse", 
    simplify = TRUE
    )

## The output has the same components as `fit_gaussian_2D()` 
## but for the automatically-selected best-fitting method only:
summary(gauss_auto)

```

## Contributing and/or raising Issues

Feedback on bugs, improvements, and/or feature requests are all welcome. 
Please see the Issues templates on GitHub to make a bug fix request or feature 
request.

To contribute code via a pull request, please consult the Contributing Guide 
first.


## Citation

Baliga, VB. 2021. gaussplotR: Fit, predict, and plot 2D-Gaussians in R. Journal of Open Source Software, 6(60), 3074. https://doi.org/10.21105/joss.03074


## License

GPL (>= 3) + file LICENSE

🐢
---
title: "Fitting 2D-Gaussians to Data"
author: "Vikram B. Baliga"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Fitting 2D-Gaussians to Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteDepends{ggplot2}
  %\VignetteDepends{metR}
  %\VignetteDepends{lattice}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Overview

The function `fit_gaussian_2D()` can be used fit 2D-Gaussians to data, and has
several methods for how the fitting is implemented. This vignette will run you
through what these methods mean with worked examples.

We'll begin by loading `gaussplotR` and loading the sample data set provided
within. The raw data we'd like to use are in columns 1:3, so we'll shave the
data set down to those columns before running through the examples.

```{r setup}
library(gaussplotR)

## We'll also use lattice, ggplot2 and metR
library(lattice); library(ggplot2); library(metR)

## Load the sample data set
data(gaussplot_sample_data)

## The raw data we'd like to use are in columns 1:3
samp_dat <-
  gaussplot_sample_data[,1:3]
```

It generally helps to plot the data beforehand to get a sense of its overall 
shape. We'll simply produce a contour plot.

```{r raw_data_contour}
lattice::levelplot(
  response ~ X_values * Y_values,
  data = samp_dat,
  col.regions = colorRampPalette(
    c("white", "blue")
    )(100),
  xlim = c(-5, 0),
  ylim = c(-1, 4), 
  asp = 1
)

```


## The `method` and `constrain_orientation` arguments

### `method`

`gaussplotR::fit_gaussian_2D()` has three main options for its `method` 
argument: 1) `"elliptical"`, 2) `"elliptical_log"`, or 3) `"circular"`. 

The most generic method (and the default) is `method = "elliptical"`. This
allows the fitted 2D-Gaussian to take an ellipsoid shape. If you would like the
best-fitting 2D-Gaussian, this is most likely your best bet.

A slightly-altered method to fit an ellipsoid Gaussian is available in 
`method = "elliptical_log"`. This method follows Priebe et al. 2003^[Priebe NJ, 
Cassanello CR, Lisberger SG. The neural representation of speed in macaque area 
MT/V5. J Neurosci. 2003 Jul 2;23(13):5650-61. doi: 
10.1523/JNEUROSCI.23-13-05650.2003.] and is geared towards use with 
log2-transformed data.

A third option is `method = "circular"`. This produces a very simple 2D-Gaussian
that is constrained to have to have a roughly circular shape (i.e. spread in X-
and Y- are roughly equal).

An additional argument, `constrain_orientation` gives additional control over the
orientation of the fitted Gaussian.  By default, the `constrain_orientation` is
`"unconstrained"`, meaning that the best-fit orientation is returned. 

### `constrain_orientation`
Setting constrain_orientation to a numeric (e.g. `constrain_orientation = pi/2`)
will force the orientation of the Gaussian to the specified value, but this is
only available when using `method = "elliptical"` or `method = "elliptical_log"`

Note that supplying a numeric to `constrain_orientation ` is handled differently 
by `method = "elliptical"` vs `method = "elliptical_log"`. With 
`method = "elliptical"`, a `theta` parameter dictates the rotation, in 
radians, from the x-axis in the clockwise direction. Thus, using 
`method = "elliptical", constrain_orientation = pi/2` will return parameters for
an elliptical 2D-Gaussian that is constrained to a 90-degree (pi/2) orientation. 
In contrast, the `method = "elliptical_log"` procedure uses a `Q` parameter to 
determine the orientation of the 2D-Gaussian. 
Setting `method = "elliptical_log", constrain_orientation = 0` will result 
in a diagonally-oriented Gaussian, whereas setting `constrain_orientation = -1`
will result in horizontal orientation. Again, see Priebe et al. 2003 for more 
details.

## Example 1: Unconstrained elliptical

Unconstrained ellipticals are the default option and are generally recommended
for most purposes. Here's an example:

```{r u_e}
gauss_fit_ue <-
    fit_gaussian_2D(samp_dat)

gauss_fit_ue
attributes(gauss_fit_ue)
```
Fitting an unconstrained ellipse returns an object (here: `gauss_fit_ue`) that
is a `data.frame` with one column per fitted parameter. The fitted parameters
are: `A_o` (a constant term), `Amp` (amplitude), `theta` (rotation, in radians,
from the x-axis in the clockwise direction), `X_peak` (x-axis peak location),
`Y_peak` (y-axis peak location), `a` (width of Gaussian along x-axis), and `b`
(width of Gaussian along y-axis).

Note that the `data.frame` in `gauss_fit_ue$fit_method` indicates the 
fitting method and whether amplitude and/or orientation were constrained.
This `data.frame` is used by `predict_gaussian_2D()` to 
automatically determine what method (and therefore, identity of parameters) was 
used and then sample points from that fitted Gaussian. 

We can elect to sample more points from the fitted Gaussian by feeding in a
grid of x- and y- values on which to predict (via `expand.grid()`. 

Then, the fitted object `gauss_fit_ue` along with the grid of points can be 
supplied to `predict_gaussian_2D` to sample more points from the fit, which
can be useful for plotting.

```{r predict_and_plot_ue}
## Generate a grid of x- and y- values on which to predict
grid <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -1, to = 4, by = 0.1))

## Predict the values using predict_gaussian_2D
gauss_data_ue <-
  predict_gaussian_2D(
    fit_object = gauss_fit_ue,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
ggplot_gaussian_2D(gauss_data_ue)
```

## Example 2: Constrained elliptical

As noted above, the `constrain_orientation` can be used to dictate the 
orientation. Please note that this will very likely result in a poorer fit, but
may be useful for certain types of analyses.

Here we'll force the Gaussian to be horizontally-oriented. 

```{r c_e}
gauss_fit_ce <-
    fit_gaussian_2D(samp_dat, 
                    constrain_orientation = 0)

gauss_fit_ce
```

We'll use the same grid of x- and y- points as above

```{r predict_and_plot_ce}
## Predict the values using predict_gaussian_2D
gauss_data_ce <-
  predict_gaussian_2D(
    fit_object = gauss_fit_ce,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
ggplot_gaussian_2D(gauss_data_ce)
```

## Example 3: Unconstrained elliptical_log

This procedure follows the formula used in Priebe et al. 2003 and is geared 
towards log2-transformed data (which the example data are).

Parameters from this model include: `Amp` (amplitude), `Q` (orientation
parameter), `X_peak` (x-axis peak location), `Y_peak` (y-axis peak location),
`X_sig` (spread along x-axis), and `Y_sig` (spread along y-axis).

```{r uel}
gauss_fit_uel <-
    fit_gaussian_2D(samp_dat, 
                    method = "elliptical_log")

gauss_fit_uel

## Predict the values using predict_gaussian_2D
gauss_data_uel <-
  predict_gaussian_2D(
    fit_object = gauss_fit_uel,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
ggplot_gaussian_2D(gauss_data_uel)
```

## Example 4: Constrained elliptical_log

Similar to the above, but here the `constrain_orientation` can be used to dictate
the value of the `Q` parameter used in Priebe et al. 2003. Setting `Q` to 0 will
result in a diagonally-oriented Gaussian, whereas setting `Q` to -1 will result
in horizontal orientation. `Q` is a continuous parameter, so values in between
may be used as well, such as in this example:

```{r cel}
gauss_fit_cel <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_orientation = -0.66
    )

gauss_fit_cel

## Predict the values using predict_gaussian_2D
gauss_data_cel <-
  predict_gaussian_2D(
    fit_object = gauss_fit_cel,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
ggplot_gaussian_2D(gauss_data_cel)
```

Again, setting the value of `Q` via `constrain_orientation` will very likely 
result in poorer-fitting Gaussians. See the analyses in Priebe et al. 2003 to
get a sense of useful applications of this approach. Forcing Q = -0.66 in 
the above example isn't all that useful, but goes to show that it can be done.

## Example 5: Circular

Using `method = "circular"` constrains the Gaussian to have a roughly circular
shape (i.e. spread in X- and Y- are roughly equal).

If this method is used, the fitted parameters are: `Amp` (amplitude), `X_peak`
(x-axis peak location), `Y_peak` (y-axis peak location), `X_sig` (spread along
x-axis), and `Y_sig `(spread along y-axis).

```{r cir}
gauss_fit_cir <-
    fit_gaussian_2D(samp_dat, 
                    method = "circular")

gauss_fit_cir

## Predict the values using predict_gaussian_2D
gauss_data_cir <-
  predict_gaussian_2D(
    fit_object = gauss_fit_cir,
    X_values = grid$X_values,
    Y_values = grid$Y_values,
  )

## Plot via ggplot2 and metR
ggplot_gaussian_2D(gauss_data_cir)
```

That's all!  

🐢
---
title: "Troubleshooting model fits"
author: "Vikram B. Baliga"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Troubleshooting model fits}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Automatically finding good initial values for parameters in a nonlinear model
(i.e. `stats::nls()`) is an art. Given that each of the formulas represented by
the `model` argument of `fit_gaussian_2D()` contains 5 to 7 parameters,
`stats::nls()` will often encounter singular gradients or step size errors. 

Code within `fit_gaussian_2D()` will first scan the supplied dataset to
guesstimate sensible initial parameters, which hopefully sidesteps these issues.
But there is no guarantee this strategy will always work.

This vignette will offer some guidance on what to do when `stats::nls()` fails
to converge, including the use of optional parameters in `fit_gaussian_2D()`
that are meant to help you address these issues.

Let's start by loading `gaussplotR` and getting some sample data loaded up:
```{r setup}
library(gaussplotR)

## Load the sample data set
data(gaussplot_sample_data)

## The raw data we'd like to use are in columns 1:3
samp_dat <-
  gaussplot_sample_data[,1:3]
```

## Singular gradients
One common problem is that of singular gradients. I will intentionally
comment out the next block of code because running it will produce the singular
gradient error, and generating errors in an R Markdown file will prevent its
rendering. Please un-comment the example below to see.

```{r singular_gradient}
## Un-comment this example if you'd like to see a singular gradient error
# gauss_fit_cir <-
#   fit_gaussian_2D(samp_dat,
#                   constrain_amplitude = TRUE,
#                   method = "circular")
```
The output from the above example should be:  
```{r error_msg_1}
#> Error in stats::nls(response ~ Amp_init * exp(-((((X_values - X_peak)^2)/(2 *  : 
#>   singular gradient
#> Called from: stats::nls(response ~ Amp_init * exp(-((((X_values - X_peak)^2)/(2 * 
#>     X_sig^2) + ((Y_values - Y_peak)^2)/(2 * Y_sig^2)))), start = c(X_peak = #> _peak_init, >
#>     Y_peak = Y_peak_init, X_sig = X_sig_init, Y_sig = Y_sig_init), 
#>     data = data, trace = verbose, control = list(maxiter = maxiter, 
#>         ...))
#> Error during wrapup: unimplemented type (29) in 'eval'
#> 
#> Error: no more error handlers available (recursive errors?); invoking 'abort' restart
#> Error during wrapup: INTEGER() can only be applied to a 'integer', not a 'unknown type #> #29'
#> Error: no more error handlers available (recursive errors?); invoking 'abort' restart
```


There are a couple tools in `gaussplotR` that can help you address this problem.

A good first step is to enable the optional argument `print_initial_params` in
`fit_gaussian_2D()` by setting it to `TRUE`. Again, please un-comment this next
block, since it will still produce an error:

```{r singular_gradient_print_params}
## Un-comment this example if you'd like to see a singular gradient error
# gauss_fit_cir <-
#   fit_gaussian_2D(samp_dat,
#                   constrain_amplitude = TRUE,
#                   method = "circular",
#                   print_initial_params = TRUE)
```

Though this block of code will not work, you will at least see something helpful
at the beginning of the error message:
```{r error_msg_2}
#> Initial parameters:
#>       Amp    X_peak    Y_peak     X_sig     Y_sig 
#> 25.725293 -2.000000  3.000000  2.482892  2.500000 
#> Error in stats::nls(response ~ Amp_init * exp(-((((X_values - X_peak)^2)/(2 *  : 
#>   singular gradient
#> Called from: stats::nls(response ~ Amp_init * exp(-((((X_values - X_peak)^2)/(2 * 
#>     X_sig^2) + ((Y_values - Y_peak)^2)/(2 * Y_sig^2)))), start = c(X_peak = #> _peak_init, >
#>     Y_peak = Y_peak_init, X_sig = X_sig_init, Y_sig = Y_sig_init), 
#>     data = data, trace = verbose, control = list(maxiter = maxiter, 
#>         ...))
#> Error during wrapup: unimplemented type (29) in 'eval'
#> 
#> Error: no more error handlers available (recursive errors?); invoking 'abort' restart
#> Error during wrapup: INTEGER() can only be applied to a 'integer', not a 'unknown type #> #29'
#> Error: no more error handlers available (recursive errors?); invoking 'abort' restart
```

Those first three lines indicate the initial values that were used. Often,
singular gradients will arise when initial values for parameters were poorly
chosen (sorry!). 

What you can do is supply your own set of initial values. To do this, use the 
optional argument `user_init` within `fit_gaussian_2D()`. You will need to 
supply a numeric vector that is of the same length as the number of parameters
for your chosen model. The values you supply must be provided in the same order
they appear in the Initial parameters message. They do not need to be named; 
values alone will suffice. 

This next example will work. I'll keep `print_initial_params` on too, since
it can be nice to see:
```{r no_singular_user_init}
## This should run with no errors
gauss_fit_cir_user <-
  fit_gaussian_2D(samp_dat,
                  constrain_amplitude = TRUE,
                  method = "circular",
                  user_init = c(25.72529, -2.5, 1.7, 1.3, 1.6),
                  print_initial_params = TRUE)
gauss_fit_cir_user
```

Note that although we are constraining the amplitude, the value of `Amp` must
still be provided (here it is `25.72529`). 

It may take some trial and error to find a set of `user_init` values that gets
your model to converge. It often makes sense to think about what values are
feasible for each parameter. For example, it should be relatively straightforward
to think of ranges of possible values for `X_peak` and `Y_peak`. I often find
that finding good initial values for the "spread" parameters (such as `X_sig`
and `Y_sig`) is the tough nut to crack, so I recommend tweaking those parameters
first.

## Additional control arguments to `nls()`
The `fit_gaussian_2D()` function also allows you to pass additional control 
arguments to `stats::nls.control()` via the `...` argument. 

To put this in more technical terms, arguments supplied to `...` are handled as:  
`stats::nls(control = list(maxiter, ...))`

Therefore, if you are interested in changing e.g. `minFactor` to `1/2048`:  
`fit_gaussian_2D(data, model, minFactor = 1/2048)`

See the Help file for `stats::nls.control()` for further guidance on what these
control arguments are.

Please also note that that tweaking `maxiter` should not be handled via `...`
but rather by the `maxiter` argument to `fit_gaussian_2D()`.


## Use parameter constraints with caution

Our scapegoat here is setting `constrain_amplitude = TRUE`. Often, when
constraining parameters in a nonlinear model, you'll find yourself in a scenario
where the QR decomposition of the gradient matrix is not of full column rank.

Constraining parameters (amplitude or orientation) will lead to poorer-fitting
Gaussians anyway, so these features should only be used if you have an *a
priori* reason to do so (see examples in Priebe et al. 2003^[Priebe NJ,
Cassanello CR, Lisberger SG. The neural representation of speed in macaque area
MT/V5. J Neurosci. 2003 Jul 2;23(13):5650-61. doi:
10.1523/JNEUROSCI.23-13-05650.2003.])

Turning off the `constrain_amplitude` constraint alleviates the problem in this
particular case:

```{r no_singular}
## This should run with no errors
gauss_fit_cir <-
  fit_gaussian_2D(samp_dat,
                  method = "circular")
gauss_fit_cir
```

Hope this helps!

🐢
---
title: "Formulas used by fit_gaussian_2D"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Formulas used by fit_gaussian_2D}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Overview

This document will provide specific details of 2D-Gaussian equations used by
the different `method` options within `gaussplotR::fit_gaussian_2D()`.

## `method = "elliptical"`

Using `method = "elliptical"` fits a two-dimensional, elliptical Gaussian
equation to gridded data.

$$G(x,y) = A_o + A * e^{-U/2}$$

where G is the value of the 2D-Gaussian at each ${(x,y)}$ point, $A_o$ is a 
constant term, and $A$ is the amplitude (i.e. scale factor).

The elliptical function, $U$, is:

$$U = (x'/a)^{2} + (y'/b)^{2}$$

where $a$ is the spread of Gaussian along the x-axis and $b$ is the spread of
Gaussian along the y-axis. 

$x'$ and $y'$ are defined as:

$$x' = (x - x_0)cos(\theta) - (y - y_0)sin(\theta)$$
$$y' = (x - x_0)sin(\theta) + (y - y_0)cos(\theta)$$
where $x_0$ is the center (peak) of the Gaussian along the x-axis, $y_0$ is the
center (peak) of the Gaussian along the y-axis, and $\theta$ is the rotation of
the ellipse from the x-axis in radians, counter-clockwise.

Therefore, all together:

$$G(x,y) = A_o + A * e^{-((((x - x_0)cos(\theta) - (y - y_0)sin(\theta))/a)^{2}+ 
(((x - x_0)sin(\theta) + (y - y_0)cos(\theta))/b)^{2})/2}$$

Setting the `constrain_orientation` argument to a numeric will optionally
constrain the value of $\theta$ to a user-specified value. If a numeric is
supplied here, please note that the value will be interpreted as a value in
radians. Constraining $\theta$ to a user-supplied value can lead to considerably
poorer-fitting Gaussians and/or trouble with converging on a stable solution; in
most cases `constrain_orientation` should remain its default: `"unconstrained"`.


## `method = "elliptical_log"`

The formula used in `method = "elliptical_log"` uses the modification of a 2D
Gaussian fit used by Priebe et al. 2003^[Priebe NJ, Cassanello CR, Lisberger SG.
The neural representation of speed in macaque area MT/V5. J Neurosci. 2003 Jul
2;23(13):5650-61. doi: 10.1523/JNEUROSCI.23-13-05650.2003.].

$$G(x,y) = A * e^{(-(x - x_0)^2)/\sigma_x^2} * e^{(-(y - y'(x)))/\sigma_y^2}$$

and 

$$y'(x) = 2^{(Q+1) * (x - x_0) + y_0}$$
where $A$ is the amplitude (i.e. scale factor), $x_0$ is the center (peak) of
the Gaussian along the x-axis, $y_0$ is the center (peak) of the Gaussian along
the y-axis, $\sigma_x$ is the spread along the x-axis, $\sigma_y$ is the spread
along the y-axis and $Q$ is an orientation parameter.

Therefore, all together:

$$G(x,y) = A * e^{(-(x - x_0)^2)/\sigma_x^2} * e^{(-(y - (2^{(Q+1) * (x - x_0) +
y_0})))/\sigma_y^2}$$

This formula is intended for use with log2-transformed data. 

Setting the `constrain_orientation` argument to a numeric will optionally
constrain the value of $Q$ to a user-specified value, which can be useful for
certain kinds of analyses (see Priebe et al. 2003 for more). Keep in mind that 
constraining $Q$ to a user-supplied value can lead to considerably
poorer-fitting Gaussians and/or trouble with converging on a stable solution; in
most cases `constrain_orientation` should remain its default: `"unconstrained"`.


## `method = "circular"`

This method uses a relatively simple formula:

$$G(x,y) = A * e^{(-(
((x-x_0)^2/2\sigma_x^2) + ((y-y_0)^2/2\sigma_y^2))
)}$$

where $A$ is the amplitude (i.e. scale factor), $x_0$ is the center (peak) of
the Gaussian along the x-axis, $y_0$ is the center (peak) of the Gaussian along
the y-axis, $\sigma_x$ is the spread along the x-axis, and $\sigma_y$ is the
spread along the y-axis.

That's all!  

🐢
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot_gaussian_2D.R
\name{ggplot_gaussian_2D}
\alias{ggplot_gaussian_2D}
\title{Plot a 2D-Gaussian via ggplot}
\usage{
ggplot_gaussian_2D(
  gauss_data,
  normalize = TRUE,
  contour_thickness = 0.04,
  contour_color = "black",
  bins = 15,
  viridis_dir = 1,
  viridis_opt = "B",
  x_lab = "X values",
  y_lab = "Y values",
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 7),
  axis.ticks = element_line(size = 0.3),
  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
  ...
)
}
\arguments{
\item{gauss_data}{Data.frame with X_values, Y_values, and predicted_values,
e.g. exported from \code{predict_gaussian_2D()}}

\item{normalize}{Default TRUE, should predicted_values be normalized on a 0
to 1 scale?}

\item{contour_thickness}{Thickness of contour lines}

\item{contour_color}{Color of the contour lines}

\item{bins}{Number of bins for the contour plot}

\item{viridis_dir}{See "direction" in scale_fill_viridis_c()}

\item{viridis_opt}{See "option" in scale_fill_viridis_c()}

\item{x_lab}{Arguments passed to xlab()}

\item{y_lab}{Arguments passed to ylab()}

\item{axis.text}{Arguments passed to axis.text}

\item{axis.title}{Arguments passed to axis.title}

\item{axis.ticks}{Arguments passed to axis.ticks}

\item{plot.margin}{Arguments passed to plot.margin}

\item{...}{Other arguments supplied to \code{ggplot2::theme()}}
}
\value{
A ggplot object that uses metR::geom_contour_fill() to display the
2D-Gaussian
}
\description{
Plot a 2D-Gaussian via ggplot
}
\examples{
if (interactive()) {
  ## Load the sample data set
  data(gaussplot_sample_data)

  ## The raw data we'd like to use are in columns 1:3
  samp_dat <-
    gaussplot_sample_data[,1:3]


  #### Example 1: Unconstrained elliptical ####
  ## This fits an unconstrained elliptical by default
  gauss_fit <-
    fit_gaussian_2D(samp_dat)

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  library(ggplot2); library(metR)
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)


  #### Example 2: Constrained elliptical_log ####
  ## This fits a constrained elliptical, as in Priebe et al. 2003
  gauss_fit <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_orientation = -1
    )

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)
}
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_gaussian_2D.R
\name{fit_gaussian_2D}
\alias{fit_gaussian_2D}
\title{Determine the best-fit parameters for a specific 2D-Gaussian model}
\usage{
fit_gaussian_2D(
  data,
  method = "elliptical",
  constrain_amplitude = FALSE,
  constrain_orientation = "unconstrained",
  user_init = NULL,
  maxiter = 1000,
  verbose = FALSE,
  print_initial_params = FALSE,
  ...
)
}
\arguments{
\item{data}{A data.frame that contains the raw data (generally rectilinearly
gridded data, but this is not a strict requirement). Columns must be named
\code{"X_values"}, \code{"Y_values"} and \code{"response"}.}

\item{method}{Choice of \code{"elliptical"}, \code{"elliptical_log"}, or
\code{"circular"}. Determine which specific implementation of 2D-Gaussian
to use. See Details for more.}

\item{constrain_amplitude}{Default FALSE; should the amplitude of the
Gaussian be set to the maximum value of the \code{"response"} variable
(\code{TRUE}), or should the amplitude fitted by \code{stats::nls()}
(\code{FALSE})?}

\item{constrain_orientation}{If using \code{"elliptical"} or \code{method =
  "elliptical_log"}, should the orientation of the Gaussian be unconstrained
(i.e. the best-fit orientation is returned) or should it be pre-set by the
user? See Details for more. Defaults to \code{"unconstrained"}.}

\item{user_init}{Default NULL; if desired, the user can supply initial values
for the parameters of the chosen model. See Details for more.}

\item{maxiter}{Default 1000. A positive integer specifying the maximum number
of iterations allowed. See \code{stats::nls.control()} for more details.}

\item{verbose}{TRUE or FALSE; should the trace of the iteration be printed?
See the \code{trace} argument of \code{stats::nls()} for more detail.}

\item{print_initial_params}{TRUE or FALSE; should the set of initial
parameters supplied to \code{stats::nls()} be printed to the console? Set
to FALSE by default to avoid confusion with the fitted parameters attained
after using \code{stats::nls()}.}

\item{...}{Additional arguments passed to \code{stats::nls.control()}}
}
\value{
A list with the components:
\itemize{
\item{"coefs"} {A data.frame of fitted model parameters.}
\item{"model"} {The model object, fitted by \code{stats::nls()}.}
\item{"model_error_stats"} {A data.frame detailing the rss, rmse, deviance,
AIC, R2 and adjusted R2 of the fitted model.}
\item{"fit_method"} {A character vector that indicates which method and
orientation strategy was used by this function.}
}
}
\description{
Determine the best-fit parameters for a specific 2D-Gaussian model
}
\details{
\code{stats::nls()} is used to fit parameters for a 2D-Gaussian to
the supplied data. Each method uses (slightly) different sets of
parameters. Note that for a small (but non-trivial) proportion of data
sets, nonlinear least squares may fail due to singularities or other
issues. Most often, this occurs because of the starting parameters that are
fed in. By default, this function attempts to set default parameters by
making an educated guess about the major aspects of the supplied data.
Should this strategy fail, the user can make use of the \code{user_init}
argument to supply an alternate set of starting values.

The simplest method is \code{method = "circular"}. Here, the 2D-Gaussian is
constrained to have a roughly circular shape (i.e. spread in X- and Y- are
roughly equal). If this method is used, the fitted parameters are: Amp
(amplitude), X_peak (x-axis peak location), Y_peak (y-axis peak location),
X_sig (spread along x-axis), and Y_sig (spread along y-axis).

A more generic method (and the default) is \code{method = "elliptical"}.
This allows the fitted 2D-Gaussian to take a more ellipsoid shape (but note
that \code{method = "circular"} can be considered a special case of this).
If this method is used, the fitted parameters are: A_o (a constant term),
Amp (amplitude), theta (rotation, in radians, from the x-axis in the
clockwise direction), X_peak (x-axis peak location), Y_peak (y-axis peak
location), a (width of Gaussian along x-axis), and b (width of Gaussian
along y-axis).

A third method is \code{method = "elliptical_log"}. This is a further
special case in which log2-transformed data may be used. See Priebe et al.
2003 for more details. Parameters from this model include: Amp (amplitude),
Q (orientation parameter), X_peak (x-axis peak location), Y_peak (y-axis
peak location), X_sig (spread along x-axis), and Y_sig (spread along
y-axis).

If using either \code{method = "elliptical"} or \code{method =
  "elliptical_log"}, the \code{"constrain_orientation"} argument can be used
to specify how the orientation is set. In most cases, the user should use
the default "unconstrained" setting for this argument. Doing so will
provide the best-fit 2D-Gaussian (assuming that the solution yielded by
\code{stats::nls()} converges on the global optimum).

Setting \code{constrain_orientation} to a numeric (e.g.
\code{constrain_orientation = pi/2}) will force the orientation of the
Gaussian to the specified value. Note that this is handled differently by
\code{method = "elliptical"} vs \code{method = "elliptical_log"}. In
\code{method = "elliptical"}, the theta parameter dictates the rotation, in
radians, from the x-axis in the clockwise direction. In contrast, the
\code{method = "elliptical_log"} procedure uses a Q parameter to determine
the orientation of the 2D-Gaussian. Setting \code{constrain_orientation =
  0} will result in a diagonally-oriented Gaussian, whereas setting
\code{constrain_orientation = -1} will result in horizontal orientation.
See Priebe et al. 2003 for more details.

The \code{user_init} argument can also be used to supply a vector of
initial values for the A, Q, X_peak, Y_peak, X_var, and Y_var parameters.
If the user chooses to make use of \code{user_init}, then a vector
containing all parameters must be supplied in a particular order.

Additional arguments to the \code{control} argument in \code{stats::nls()}
can be supplied via \code{...}.
}
\examples{
if (interactive()) {
  ## Load the sample data set
  data(gaussplot_sample_data)

  ## The raw data we'd like to use are in columns 1:3
  samp_dat <-
    gaussplot_sample_data[,1:3]


  #### Example 1: Unconstrained elliptical ####
  ## This fits an unconstrained elliptical by default
  gauss_fit <-
    fit_gaussian_2D(samp_dat)

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  library(ggplot2); library(metR)
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)


  #### Example 2: Constrained elliptical_log ####
  ## This fits a constrained elliptical, as in Priebe et al. 2003
  gauss_fit <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_orientation = -1
    )

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)
}
}
\references{
Priebe NJ, Cassanello CR, Lisberger SG. The neural representation
of speed in macaque area MT/V5. J Neurosci. 2003 Jul 2;23(13):5650-61. doi:
10.1523/JNEUROSCI.23-13-05650.2003.
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autofit_gaussian_2D.R
\name{autofit_gaussian_2D}
\alias{autofit_gaussian_2D}
\title{Automatically determine the best-fitting 2D-Gaussian for a data set}
\usage{
autofit_gaussian_2D(
  data,
  comparison_method = "rmse",
  maxiter = 1000,
  simplify = TRUE
)
}
\arguments{
\item{data}{A data.frame that contains the raw data (generally rectilinearly
gridded data, but this is not a strict requirement). Columns must be named
\code{"X_values"}, \code{"Y_values"} and \code{"response"}.}

\item{comparison_method}{One of "rmse", "rss", or "AIC"; what metric should
be used to determine the "best-fitting" Gaussian?}

\item{maxiter}{Default 1000. A positive integer specifying the maximum number
of iterations allowed. See \code{stats::nls.control()} for more details.}

\item{simplify}{TRUE or FALSE. If TRUE, return only the coefficients, model,
model_error_stats, and fit_method for the best-fitting model. If FALSE, a
model comparison table is also included in the returned list as
\code{$model_comparison}. This table is obtained via
\code{compare_gaussian_fits()}.}
}
\value{
If \code{simplify = TRUE}, a list with the components:
\itemize{
\item{"coefs"} {A data.frame of fitted model parameters.}
\item{"model"} {The model object, fitted by \code{stats::nls()}.}
\item{"model_error_stats"} {A data.frame detailing the rss, rmse, deviance,
and AIC of the fitted model.}
\item{"fit_method"} {A character vector that indicates which method and
orientation strategy was used by this function.}
}

If \code{simplify = FALSE}, a model comparison table is also included
in the returned list as \code{$model_comparison}. This table is obtained
via \code{compare_gaussian_fits()}.
}
\description{
Automatically determine the best-fitting 2D-Gaussian for a data set
}
\details{
This function runs \code{fit_gaussian_2D()} three times: once for
each of the "main" types of models: 1) elliptical, unconstrained; 2)
elliptical, log; 3) circular. In all three cases, amplitudes and
orientations are unconstrained. The function \code{compare_gaussian_fits()}
is then used to determine which of these three models is the best-fitting,
using the \code{comparison_method} argument to make the decision.
}
\examples{
if (interactive()) {
}
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_gaussian_2D.R
\name{predict_gaussian_2D}
\alias{predict_gaussian_2D}
\title{Predict values from a fitted 2D-Gaussian}
\usage{
predict_gaussian_2D(fit_object, X_values, Y_values, ...)
}
\arguments{
\item{fit_object}{Either the output of \code{gaussplotR::fit_gaussian_2D()}
or a list that contains coefficients and fit methods (see Details).}

\item{X_values}{vector of numeric values for the x-axis}

\item{Y_values}{vector of numeric values for the y-axis}

\item{...}{Additional arguments}
}
\value{
A data.frame with the supplied \code{X_values} and \code{Y_values}
along with the predicted values of the 2D-Gaussian
(\code{predicted_values})
}
\description{
Predict values from a fitted 2D-Gaussian
}
\details{
This function assumes Gaussian parameters have been fitted beforehand. No
fitting of parameters is done within this function; these can be
supplied via the object created by \code{gaussplotR::fit_gaussian_2D()}.

If \code{fit_object} is not an object created by
\code{gaussplotR::fit_gaussian_2D()}, \code{predict_gaussian_2D()} attempts
to parse \code{fit_object} as a list of two items. The coefficients of the
fit must be supplied as a one-row, named data.frame within
\code{fit_object$coefs}, and details of the methods for fitting the Gaussian
must be contained as a character vector in \code{fit_object$fit_method}. This
character vector in \code{fit_object$fit_method} must be a named vector that
provides information about the method, amplitude constraint choice, and
orientation constraint choice, using the names \code{method},
\code{amplitude}, and \code{orientation}. \code{method} must be one of:
\code{"elliptical"}, \code{"elliptical_log"}, or \code{"circular"}.
\code{amplitude} and \code{orientation} must each be either
\code{"unconstrained"} or \code{"constrained"}. For example, \code{c(method =
"elliptical", amplitude = "unconstrained", orientation = "unconstrained")}.
One exception to this is when \code{method = "circular"}, in which case
\code{orientation} must be \code{NA}, e.g.: \code{c(method = "circular",
amplitude = "unconstrained", orientation = NA)}.
}
\examples{
if (interactive()) {
  ## Load the sample data set
  data(gaussplot_sample_data)

  ## The raw data we'd like to use are in columns 1:3
  samp_dat <-
    gaussplot_sample_data[,1:3]


  #### Example 1: Unconstrained elliptical ####
  ## This fits an unconstrained elliptical by default
  gauss_fit <-
    fit_gaussian_2D(samp_dat)

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  library(ggplot2); library(metR)
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)


  #### Example 2: Constrained elliptical_log ####
  ## This fits a constrained elliptical, as in Priebe et al. 2003
  gauss_fit <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_orientation = -1
    )

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)
}
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgl_gaussian_2D.R
\name{rgl_gaussian_2D}
\alias{rgl_gaussian_2D}
\title{Produce a 3D plot of the 2D-Gaussian via rgl}
\usage{
rgl_gaussian_2D(
  gauss_data,
  normalize = TRUE,
  viridis_dir = 1,
  viridis_opt = "B",
  x_lab = "X values",
  y_lab = "Y values",
  box = FALSE,
  aspect = TRUE,
  ...
)
}
\arguments{
\item{gauss_data}{Data.frame with X_values, Y_values, and predicted_values,
e.g. exported from \code{predict_gaussian_2D()}}

\item{normalize}{Default TRUE, should predicted_values be normalized on a 0
to 1 scale?}

\item{viridis_dir}{See "direction" in scale_fill_viridis_c()}

\item{viridis_opt}{See "option" in scale_fill_viridis_c()}

\item{x_lab}{Arguments passed to xlab()}

\item{y_lab}{Arguments passed to ylab()}

\item{box}{Whether to draw a box; see \code{rgl::plot3d()}}

\item{aspect}{Whether to adjust the aspect ratio; see \code{rgl::plot3d()}}

\item{...}{Other arguments supplied to \code{rgl::plot3d()}}
}
\value{
An rgl object (i.e. of the class 'rglHighlevel'). See
\code{rgl::plot3d()} for details.
}
\description{
Produce a 3D plot of the 2D-Gaussian via rgl
}
\examples{
if (interactive()) {
  ## Load the sample data set
  data(gaussplot_sample_data)

  ## The raw data we'd like to use are in columns 1:3
  samp_dat <-
    gaussplot_sample_data[,1:3]


  #### Example 1: Unconstrained elliptical ####
  ## This fits an unconstrained elliptical by default
  gauss_fit <-
    fit_gaussian_2D(samp_dat)

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  library(ggplot2); library(metR)
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)


  #### Example 2: Constrained elliptical_log ####
  ## This fits a constrained elliptical, as in Priebe et al. 2003
  gauss_fit <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_orientation = -1
    )

  ## Generate a grid of x- and y- values on which to predict
  grid <-
    expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
                Y_values = seq(from = -1, to = 4, by = 0.1))

  ## Predict the values using predict_gaussian_2D
  gauss_data <-
    predict_gaussian_2D(
      fit_object = gauss_fit,
      X_values = grid$X_values,
      Y_values = grid$Y_values,
    )

  ## Plot via ggplot2 and metR
  ggplot_gaussian_2D(gauss_data)

  ## Produce a 3D plot via rgl
  rgl_gaussian_2D(gauss_data)
}
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/characterize_gaussian_fits.R
\name{characterize_gaussian_fits}
\alias{characterize_gaussian_fits}
\title{Characterize the orientation of fitted 2D-Gaussians}
\usage{
characterize_gaussian_fits(
  fit_objects_list = NULL,
  data = NULL,
  constrain_amplitude = FALSE,
  ...
)
}
\arguments{
\item{fit_objects_list}{A list of outputs from \code{fit_gaussian_2D()}. See
Details for more. This is the preferred input object for this function.}

\item{data}{A data.frame that contains the raw data (generally rectilinearly
gridded data, but this is not a strict requirement). Columns must be named
\code{"X_values"}, \code{"Y_values"} and \code{"response"}. See
\code{fit_gaussian_2D()} for details.}

\item{constrain_amplitude}{Default FALSE; should the amplitude of the
Gaussian be set to the maximum value of the \code{"response"} variable
(\code{TRUE}), or should the amplitude fitted by \code{stats::nls()}
(\code{FALSE})? See \code{fit_gaussian_2D()} for details.}

\item{...}{Additional arguments that can be passed to
\code{fit_gaussian_2D()} if data are supplied.}
}
\value{
A list with the following:
\itemize{
\item{"model_comparison"} {A model comparison output (i.e. what is produced
by \code{compare_gaussian_fits()}), which indicates the relative preference
of each of the three models.}
\item{"Q_table"} {A data.frame that provides information on the value of Q
from the best-fitting model, along with the 5-95\% confidence intervals of
this estimate.}
\item{"r_i"} {A numeric, the correlation of the data with the independent
(Q = -1) prediction.}
\item{"r_s"} {A numeric, the correlation of the data with the diagonally-
oriented (Q = 0) prediction.}
\item{"r_is"} {A numeric, the correlation between the independent
(Q = -1) prediction and the the diagonally-oriented (Q = 0) prediction.}
\item{"R_indp"} {A numeric, partial correlation of the response variable
with the independent (Q = -1) prediction.}
\item{"R_diag"} {A numeric, partial correlation of the response variable
with the diagonally-oriented (Q = 0) prediction.}
\item{"ZF_indp"} {A numeric, the Fisher Z-transform of the R_indp
coefficient. See Winship et al. 2006 for details.}
\item{"ZF_diag"} {A numeric, the Fisher Z-transform of the R_diag
coefficient. See Winship et al. 2006 for details.}
\item{"Z_diff"} {A numeric, the Z-difference between ZF_indp and ZF_diag.
See Winship et al. 2006 for details.}

}
}
\description{
The orientation and partial correlations of Gaussian data are analyzed
according to Levitt et al. 1994 and Priebe et al. 2003. Features include
computation of partial correlations between response variables and
independent and diagonally-tuned predictions, along with Z-difference
scoring.
}
\details{
This function accepts either a list of objects output from
\code{fit_gaussian_2D()} (preferred) or a data.frame that contains the raw
data.

The supplied fit_objects_list must be a list that contains objects returned
by \code{fit_gaussian_2D()}. This list must contain exactly three models. All
three models must have been run using \code{method = "elliptical_log"}. The
models must be: 1) one in which orientation is unconstrained, 2) one in which
orientation is constrained to Q = 0 (i.e. a diagonally-oriented Gaussian),
and 3) one in which orientation is constrained to Q = -1 (i.e. a
horizontally-oriented Gaussian). See this function's Examples for guidance.

Should raw data be provided instead of the fit_objects_list, the
\code{characterize_gaussian_fits()} runs \code{fit_gaussian_2D()} internally.
This is generally not recommended, as difficulties in fitting models via
\code{stats::nls()} are more easily troubleshot by the optional arguments in
\code{fit_gaussian_2D()}. Nevertheless, supplying raw data instead of a list
of fitted models is feasible, though your mileage may vary.
}
\examples{
if (interactive()) {
  library(gaussplotR)

  ## Load the sample data set
  data(gaussplot_sample_data)

  ## The raw data we'd like to use are in columns 1:3
  samp_dat <-
    gaussplot_sample_data[,1:3]

  ## Fit the three required models
  gauss_fit_uncn <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_amplitude = FALSE,
      constrain_orientation = "unconstrained"
    )

  gauss_fit_diag <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_amplitude = FALSE,
      constrain_orientation = 0
    )

  gauss_fit_indp <-
    fit_gaussian_2D(
      samp_dat,
      method = "elliptical_log",
      constrain_amplitude = FALSE,
      constrain_orientation = -1
    )

  ## Combine the outputs into a list
  models_list <-
    list(
      gauss_fit_uncn,
      gauss_fit_diag,
      gauss_fit_indp
    )

  ## Now characterize
  out <-
    characterize_gaussian_fits(models_list)
  out

  ## Alternatively, the raw data itself can be supplied.
  ## This is less preferred, as fitting of models may fail
  ## internally.
  out2 <-
    characterize_gaussian_fits(data = samp_dat)

  ## This produces the same output, assuming models are fit without error
  identical(out, out2)
}
}
\references{
Levitt JB, Kiper DC, Movshon JA. Receptive fields and functional architecture
of macaque V2. J Neurophysiol. 1994 71:2517–2542.

Priebe NJ, Cassanello CR, Lisberger SG. The neural representation of speed in
macaque area MT/V5. J Neurosci. 2003 Jul 2;23(13):5650-61. doi:
10.1523/JNEUROSCI.23-13-05650.2003.

Winship IR, Crowder N, Wylie DRW. Quantitative reassessment of speed tuning
in the accessory optic system and pretectum of pigeons. J Neurophysiol. 2006
95(1):546-551. doi: 10.1152/jn.00921.2005
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_gaussian_fits.R
\name{compare_gaussian_fits}
\alias{compare_gaussian_fits}
\title{Compare fitted 2D-Gaussians and determine the best-fitting model}
\usage{
compare_gaussian_fits(fit_objects_list, comparison_method = "rmse")
}
\arguments{
\item{fit_objects_list}{A list of outputs from \code{fit_gaussian_2D()}. See
Details for more}

\item{comparison_method}{One of "rmse", "rss", or "AIC"; what metric should
be used to determine the "best-fitting" Gaussian?}
}
\value{
A list with the components:
\itemize{
\item{"preferred_model"} {A character indicating the name of the preferred
model (or if a named list was not provided, a model number is given in
the order of the original supplied list). }
\item{"comparison_table"} {A data.frame detailing the rss, rmse, deviance,
, AIC, R2, and adjusted R2 of the fitted models. The data.frame is sorted
by the comparison_method that was selected.}
}
}
\description{
Compare fitted 2D-Gaussians and determine the best-fitting model
}
\details{
For the argument \code{fit_objects_list}, a list of fitted model
objects (output from \code{fit_gaussian_2D()}) can simply be combined via
\code{list()}. Naming the list is optional; should you supply names, the
output of \code{compare_gaussian_fits()} will refer to specific models by
these names.
}
\examples{
if (interactive()) {
  library(gaussplotR)

  ## Load the sample data set
  data(gaussplot_sample_data)

  ## The raw data we'd like to use are in columns 1:3
  samp_dat <-
    gaussplot_sample_data[,1:3]

  ## Fit a variety of different models
  gauss_fit_ue <-
    fit_gaussian_2D(samp_dat)
  gauss_fit_uel <-
    fit_gaussian_2D(samp_dat, method = "elliptical_log")
  gauss_fit_cir <-
    fit_gaussian_2D(samp_dat, method = "circular")

  ## Combine the outputs into a list
  models_list <-
    list(
      unconstrained_elliptical = gauss_fit_ue,
      unconstrained_elliptical_log = gauss_fit_uel,
      circular = gauss_fit_cir
    )

  ## Compare via rmse
  models_compared <-
    compare_gaussian_fits(
      fit_objects_list = models_list,
      comparison_method = "rmse" ## the default
    )
}
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussplot_sample_data.R
\docType{data}
\name{gaussplot_sample_data}
\alias{gaussplot_sample_data}
\title{Sample data set}
\format{
A data frame with 36 rows and 11 variables:
\describe{
\item{X_values}{vector of numeric values for the x-axis}
\item{Y_values}{vector of numeric values for the y-axis}
\item{response}{vector of numeric values for the response variable}
\item{norm_g_resp}{normalized values from the 2D-Gaussian fit}
\item{g_resp}{values from the 2D-Gaussian fit}
\item{A}{amplitude of 2D-Gaussian (repeated)}
\item{X_peak}{location of peak x-axis value (repeated)}
\item{X_var}{variance in x (repeated)}
\item{Q}{orientation parameter of the gaussian (repeated)}
\item{Y_peak}{location of peak y-axis value (repeated)}
\item{Y_var}{variance in y (repeated)}
}
}
\usage{
gaussplot_sample_data
}
\description{
A \code{data.frame} of raw data and fitted 2D-Gaussian parameters; intended
for use with \code{predict_gaussian_2D()}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_volume_gaussian_2D.R
\name{get_volume_gaussian_2D}
\alias{get_volume_gaussian_2D}
\title{Compute volume under 2D-Gaussian}
\usage{
get_volume_gaussian_2D(X_sig, Y_sig)
}
\arguments{
\item{X_sig}{numeric value(s) of the x-axis spread (sigma)}

\item{Y_sig}{numeric value(s) of the y-axis spread (sigma)}
}
\value{
Numeric value(s) indicating the computed volume(s)
}
\description{
Compute volume under 2D-Gaussian
}
\details{
Volume under the 2D-Gaussian is computed as:
\code{2 * pi * sqrt(abs(X_sig)) * sqrt(abs(Y_sig))}

Numeric vectors can be supplied to \code{X_sig} and \code{Y_sig}. If vectors
of length greater than 1 are given, the function computes volume for each
sequential pair of \code{X_sig}, \code{Y_sig} values. The lengths of these
supplied vectors must be identical.
}
\examples{
library(gaussplotR)

get_volume_gaussian_2D(5, 3) #24.33467
}
\author{
Vikram B. Baliga
}
