
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
