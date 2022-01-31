
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ppdiag

<!-- badges: start -->

[![R build
status](https://github.com/OwenWard/ppdiag/workflows/R-CMD-check/badge.svg)](https://github.com/OwenWard/ppdiag/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/ppdiag)](https://CRAN.R-project.org/package=ppdiag)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Codecov test
coverage](https://codecov.io/gh/OwenWard/ppdiag/branch/main/graph/badge.svg)](https://codecov.io/gh/OwenWard/ppdiag?branch=main)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03133/status.svg)](https://doi.org/10.21105/joss.03133)
<!-- badges: end -->

`ppdiag` is an `R` package which provides a collection of tools which
can be used to assess the fit of temporal point processes to data.

These currently include:

-   Simulating data from a specified point process
-   Fitting a specified point process model to data
-   Evaluating the fit of a point process model to data using several
    diagnostic tools

# Installation

You can install the released version of ppdiag from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ppdiag")
```

The current development version of this package is available from
[GitHub](https://github.com/OwenWard/ppdiag) with:

``` r
# install.packages("remotes")
remotes::install_github("OwenWard/ppdiag")
```

# Example

To illustrate some of the basic functionality of this package, we can
simulate data from a specified Hawkes process and examine our diagnostic
results when we fit a homogeneous Poisson process to this data.

``` r
library(ppdiag)

hp_obj <- pp_hp(lambda0 = 0.2, alpha = 0.35, beta = 0.8)
sim_hp <- pp_simulate(hp_obj, end = 200)
sim_hp
#>  [1]   1.275239   4.783765   5.594645   8.598805  12.615358  13.236031
#>  [7]  16.646178  17.963423  18.111810  22.084071  26.666076  34.308807
#> [13]  34.356333  34.495016  34.951780  35.092074  36.397702  37.473565
#> [19]  37.846293  37.999420  54.822306  54.960122  55.721565  56.161485
#> [25]  56.825700  57.272058  59.441202  67.259184  67.951046  73.067622
#> [31]  73.376321  73.864017  74.351076  78.675743  84.520527  86.082185
#> [37]  86.547380  89.036582  99.411569 100.220569 101.941447 104.265342
#> [43] 106.553315 115.496473 126.077235 126.679327 126.705392 130.829836
#> [49] 134.226466 135.613464 135.633740 149.809218 156.366308 156.732731
#> [55] 157.273463 160.788531 161.764239 166.976330 187.590412 187.737997
#> [61] 187.994724 195.153693
```

We can readily evaluate the fit of a homogeneous Poisson process to this
data.

``` r
est_hpp <- fithpp(sim_hp)
est_hpp
#> Homogeneous Poisson Process 
#> lambda  
#> events 1.275239 4.783765 5.594645 8.598805 12.61536 13.23603 16.64618 17.96342 18.11181 22.08407 26.66608 34.30881 34.35633 34.49502 34.95178 35.09207 36.3977 37.47357 37.84629 37.99942 54.82231 54.96012 55.72156 56.16149 56.8257 57.27206 59.4412 67.25918 67.95105 73.06762 73.37632 73.86402 74.35108 78.67574 84.52053 86.08219 86.54738 89.03658 99.41157 100.2206 101.9414 104.2653 106.5533 115.4965 126.0772 126.6793 126.7054 130.8298 134.2265 135.6135 135.6337 149.8092 156.3663 156.7327 157.2735 160.7885 161.7642 166.9763 187.5904 187.738 187.9947 195.1537

pp_diag(est_hpp, events = sim_hp)
```

<img src="man/figures/README-fit_hpp-1.png" width="75%" />

    #> 
    #> Raw residual: -7.105427e-15
    #> Pearson residual: 1.421085e-14
    #> 
    #>  One-sample Kolmogorov-Smirnov test
    #> 
    #> data:  r
    #> D = 0.20838, p-value = 0.007712
    #> alternative hypothesis: two-sided

``` r
hp_est <- fithp(events = sim_hp)
pp_diag(hp_est, events = sim_hp)
```

<img src="man/figures/README-fit_hp-1.png" width="75%" />

    #> Raw residual: -0.002513104
    #> Pearson residual: 0.2040545
    #> 
    #>  One-sample Kolmogorov-Smirnov test
    #> 
    #> data:  r
    #> D = 0.075428, p-value = 0.846
    #> alternative hypothesis: two-sided

## Markov Modulated Hawkes Process Example

This is particularly useful for more complex point processes, such as
the Markov Modulated Hawkes Process (MMHP). We can simulate events from
this model and examine the fit of simpler point processes to this data.

``` r
Q <- matrix(c(-0.2, 0.2, 0.1, -0.1), ncol = 2, byrow = TRUE)

mmhp_obj <- pp_mmhp(Q, delta = c(1 / 3, 2 / 3), 
          lambda0 = 0.2,
          lambda1 = .75,
          alpha = 0.4,
          beta = 0.8)

mmhp_obj
#> Markov Modulated Hawkes Process 
#> lambda0  0.2 
#> lambda1  0.75 
#> alpha  0.4 
#> beta  0.8 
#> Q  -0.2 0.1 0.2 -0.1 
#> delta 0.3333333 0.6666667
mmhp_events <- pp_simulate(mmhp_obj, n = 50)
```

We can easily fit a homogeneous Poisson process and visualise the
goodness of fit.

``` r
est_hpp <- fithpp(events = mmhp_events$events)
pp_diag(est_hpp,mmhp_events$events)
```

<img src="man/figures/README-fit_hpp_to_mmhp-1.png" width="75%" />

    #> 
    #> Raw residual: -1
    #> Pearson residual: -1.270479
    #> 
    #>  One-sample Kolmogorov-Smirnov test
    #> 
    #> data:  r
    #> D = 0.30169, p-value = 0.000156
    #> alternative hypothesis: two-sided

Similarly for a Hawkes process.

``` r
est_hp <- fithp(events = mmhp_events$events)
pp_diag(est_hp,mmhp_events$events)
```

<img src="man/figures/README-fit_hp_to_mmhp-1.png" width="75%" />

    #> Raw residual: -0.3695538
    #> Pearson residual: -1.850818
    #> 
    #>  One-sample Kolmogorov-Smirnov test
    #> 
    #> data:  r
    #> D = 0.081193, p-value = 0.87
    #> alternative hypothesis: two-sided

We can then compare to the true point process model.

``` r
pp_diag(mmhp_obj, mmhp_events$events)
```

<img src="man/figures/README-fit_mmhp-1.png" width="75%" />

    #> Raw residual: 6.507162
    #> Pearson residual: 3.326864
    #> 
    #>  One-sample Kolmogorov-Smirnov test
    #> 
    #> data:  r
    #> D = 0.1025, p-value = 0.6324
    #> alternative hypothesis: two-sided

# Getting help and contributing

Please file any issues
[here](https://github.com/OwenWard/ppdiag/issues). Similarly, we would
be delighted if anyone would like to contribute to this package (such as
adding other point processes, kernel functions). Feel free to take a
look
[here](https://github.com/OwenWard/ppdiag/blob/main/CONTRIBUTING.md) and
reach out with any questions.

# References

-   Sun et al., (2021). ppdiag: Diagnostic Tools for Temporal Point
    Processes. Journal of Open Source Software, 6(61), 3133,
    <https://doi.org/10.21105/joss.03133>
-   Wu et al. (2021), Diagnostics and Visualization of Point Process
    Models for Event Times on a Social Network, In Applied Modeling
    Techniques and Data Analysis 1 (eds Y. Dimotikalis, A.
    Karagrigoriou, C. Parpoula and C.H. Skiadas).
    <https://doi.org/10.1002/9781119821588.ch7>
# ppdiag 0.1.1

Some small updates which arose from the JOSS review process which 

- added print methods for point process objects
- tidy up the returned information for diagnostic functions

Also updated a citation in the description to published version.
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or imagery, derogatory comments or personal attacks, trolling, public or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by opening an issue or contacting the package maintainer(ogw2103@columbia.edu).

This Code of Conduct is adapted from the [Contributor Covenant](http:contributor-covenant.org), version 1.0.0, available at https://www.contributor-covenant.org/version/1/0/0/code-of-conduct.html
# Contributing to ppdiag development

Bug reports and any potential problems with any existing code can be addressed
by filing an issue giving some brief details and if possible
a small reproducible example. If you would like to fix the bug 
yourself and it is a small fix, feel free to go ahead with a direct
pull request. Larger issues might be better proceeded by some discussion.

We would be delighted for anyone who is interested in contributing to
this package to get involved also. Potential additions would be to extend the
included diagnostic tests to other point processes, such as different
Hawkes kernel functions. 

To contribute, please look at the existing code and feel free to file
an issue to discuss further any potential additions.

Please note that ppdiag is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, 
you agree to abide by its terms.
## Test environments
* local Windows R 4.1.0
* ubuntu 20.04 (github actions) R 4.1.1
* macOS latest (github actions) R 4.1.0
* win-builder (devel and release)

## R CMD check results 

There were no ERRORs or WARNINGs

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "75%"
)
```

# ppdiag

<!-- badges: start -->
[![R build status](https://github.com/OwenWard/ppdiag/workflows/R-CMD-check/badge.svg)](https://github.com/OwenWard/ppdiag/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/ppdiag)](https://CRAN.R-project.org/package=ppdiag)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Codecov test coverage](https://codecov.io/gh/OwenWard/ppdiag/branch/main/graph/badge.svg)](https://codecov.io/gh/OwenWard/ppdiag?branch=main)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03133/status.svg)](https://doi.org/10.21105/joss.03133)
<!-- badges: end -->


``ppdiag`` is an `R` package which provides a collection of tools which
can be used to assess the fit of temporal point processes to data.

These currently include:

- Simulating data from a specified point process
- Fitting a specified point process model to data
- Evaluating the fit of a point process model to data using 
several diagnostic tools

# Installation

You can install the released version of ppdiag from [CRAN](https://CRAN.R-project.org) with:

```{r, eval=FALSE}
install.packages("ppdiag")
```

The current development version of this
package is available from [GitHub](https://github.com/OwenWard/ppdiag) with:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("OwenWard/ppdiag")
```
# Example

To illustrate some of the basic functionality of this package,
we can simulate data from a specified Hawkes process and examine
our diagnostic results when we fit a homogeneous Poisson process
to this data.

```{r example}
library(ppdiag)

hp_obj <- pp_hp(lambda0 = 0.2, alpha = 0.35, beta = 0.8)
sim_hp <- pp_simulate(hp_obj, end = 200)
sim_hp

```
We can readily evaluate the fit of a homogeneous Poisson process to this
data.

```{r fit_hpp}
est_hpp <- fithpp(sim_hp)
est_hpp

pp_diag(est_hpp, events = sim_hp)
```


```{r fit_hp}
hp_est <- fithp(events = sim_hp)
pp_diag(hp_est, events = sim_hp)

```


## Markov Modulated Hawkes Process Example

This is particularly useful for more complex point processes, such as the 
Markov Modulated Hawkes Process (MMHP).
We can simulate events from this model
and examine the fit of simpler point processes to this data.

```{r mmhp example}
Q <- matrix(c(-0.2, 0.2, 0.1, -0.1), ncol = 2, byrow = TRUE)

mmhp_obj <- pp_mmhp(Q, delta = c(1 / 3, 2 / 3), 
          lambda0 = 0.2,
          lambda1 = .75,
          alpha = 0.4,
          beta = 0.8)

mmhp_obj
mmhp_events <- pp_simulate(mmhp_obj, n = 50)
```

We can easily fit a homogeneous Poisson process and visualise the goodness of 
fit.

```{r fit_hpp_to_mmhp}
est_hpp <- fithpp(events = mmhp_events$events)
pp_diag(est_hpp,mmhp_events$events)
```

Similarly for a Hawkes process.

```{r fit_hp_to_mmhp}
est_hp <- fithp(events = mmhp_events$events)
pp_diag(est_hp,mmhp_events$events)
```

We can then compare to the true point process model.

```{r fit_mmhp}
pp_diag(mmhp_obj, mmhp_events$events)
```


# Getting help and contributing

Please file any issues [here](https://github.com/OwenWard/ppdiag/issues). 
Similarly, we would be delighted if anyone would like to contribute to
this package (such as adding other point processes, kernel functions). 
Feel free to take a look
[here](https://github.com/OwenWard/ppdiag/blob/main/CONTRIBUTING.md)
and reach out with any questions.

# References

- Sun et al., (2021). ppdiag: Diagnostic Tools for Temporal Point Processes. Journal of Open Source Software, 6(61), 3133, https://doi.org/10.21105/joss.03133
- Wu et al. (2021), Diagnostics and Visualization of Point Process Models for Event Times on a Social Network, 
In Applied Modeling Techniques and Data Analysis 1 (eds Y. Dimotikalis, A. Karagrigoriou, C. Parpoula and C.H. Skiadas).
https://doi.org/10.1002/9781119821588.ch7

---
title: "`ppdiag`, diagnostic tools for temporal Point Processes"
output: rmarkdown::html_vignette
author: Sally Sun, Owen G. Ward, Xiaoxi Zhao, Jing Wu, Tian Zheng.
vignette: >
  %\VignetteIndexEntry{`ppdiag`, diagnostic tools for temporal Point Processes}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  markdown: 
    wrap: 72
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(100) # to make it reproducible
```

```{r setup}
# remotes::install_github("OwenWard/ppdiag")
library(ppdiag)
```

This vignette provides an introduction to the functions available in
`ppdiag` to evaluate the fit of univariate temporal point processes.

To achieve this, we currently include a range of functions which allow a
user to:

-   Simulate data from a range of common univariate point processes.
-   Fit a range of univariate point processes to data.
-   After fitting a point process to some data, evaluate the ability of
    that point process to capture the temporal structure present in this
    data.

## Classes

We create classes for each of the point process models included in the 
package. Currently, these are:

-   Homogeneous Poisson Process

`pp_hpp(lambda)` creates a `hpp` object with
rate parameter `lambda`.

```{r}
hpp_obj <- pp_hpp(lambda = 1)
hpp_obj
```

-   Hawkes Process:

`pp_hp(lambda0, alpha, beta, events = NULL)` creates a `hp` object.

```{r}
hp_obj <- pp_hp(lambda0 = 0.5, alpha = 0.2, beta = 0.5)
hp_obj
```

-   Markov Modulated Poisson Process:
`pp_mmpp(lambda0, lambda1, alpha, beta, Q, delta)` creates an `mmpp` object.

```{r}
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)

mmpp_obj <- pp_mmpp(Q, delta = c(1 / 3, 2 / 3), 
          lambda0 = 0.8,
          c = 1.2)

mmpp_obj
```

-   Markov-Modulated Hawkes Process:

`pp_mmhp(lambda0, lambda1, alpha, beta, Q, delta)` creates an `mmhp` object.

```{r}
mmhp_obj <- pp_mmhp(Q, delta = c(1 / 3, 2 / 3), 
          lambda0 = 0.2,
          lambda1 = .75,
          alpha = 0.1,
          beta = 0.2)

mmhp_obj

```

## Simulating data

To simulate data from a given point process, we use the function
`pp_simulate(pp_obj, ...)`. Here the first argument specifies one of the
above point processes, while the remaining arguments specify either the
number of events simulated or the length of the observation period for
possible events.

For example, we can simulate events up to a specified end time.

```{r}
hpp_events <- pp_simulate(hpp_obj, end = 10)
hpp_events
```

Alternatively, we can specify the number of events we wish to simulate.

```{r}
hp_events <- pp_simulate(hp_obj, start = 0, n = 20)
hp_events

```

This returns the simulated events of the specified point process. For
Markov Modulated processes, the states (and the times of these states)
are also returned. In this scenario only a specified number of events
can be simulated (currently).

```{r}
mmhp_events <- pp_simulate(object = mmhp_obj, n = 20)
mmhp_events
```

## Fitting a point process

For completeness, we include functions for fitting both homogeneous
Poisson and Hawkes processes to data. Fitting a Markov modulated model
is more complex, although we describe this procedure in an included
vignette.

`fithpp(hpp_events)` returns an object of class `hpp`, estimating the
MLE of a homogenous Poisson process for `hpp_events`

```{r}
fit_hpp <- fithpp(hpp_events)
fit_hpp
```

Similarly, `fithp(hp_events)` returns an object of class `hp`,
estimating the three parameters of the Hawkes process from `hp_events`
using `constrOptim`. This ensures that the returned solution (if one can
be obtained), satisfies the stationary condition of a Hawkes process.

```{r}
hp_events <- pp_simulate(hp_obj, n = 500)
fit_hp <- fithp(hp_events)
fit_hp$lambda0
fit_hp$alpha
fit_hp$beta
```

## Diagnosing the fit of a point process to data

The main goal of this package is to provide users with tools
to examine the fit of a specified point process to some data.
There are several methods which can be used to assess the
goodness of fit of a point process to temporal data. In this package we
allow a user to:

-   Visually inspect the estimated intensity of the point process.
-   Examine the fitted intensity along with the distribution of rescaled 
    inter-event times to help identify causes for lack of fit.
-   Examine the distribution of the rescaled inter-event times, by
    utilising the time rescaling theorem.
-   Examine the residual process of an estimated point process, in
    particular computing the raw and Pearson residuals for a given point
    process fit to data.


### Visualize the intensity function

`drawHPPIntensity(hpp, events)`
plots the intensity of a homogeneous Poisson process.

```{r,fig.width=4,fig.height=4}
drawHPPIntensity(fit_hpp, events = hpp_events,
                 color = "red")
```

Similarly, `drawHPIntensity(hp, events)`
plots the intensity of a Hawkes process. 

```{r,fig.width=4,fig.height=5}
drawHPIntensity(fit_hp, events = hp_events)
```

To plot the fitted intensity on the input events, set `fit=TRUE`.

```{r,fig.width=4,fig.height=5}
drawHPIntensity(events = hp_events, fit = TRUE)
```


Similarly, 
`drawUniMMHPIntensity(mmhp, mmhp_events)`
plots the intensity of a Markov modulated 
Hawkes process, with a similar 
function for Markov modulated Poisson processes. This requires both the point
process object and the output from `pp_simulate` which describes
the latent process.

```{r,fig.width=4,fig.height=5}
drawUniMMHPIntensity(mmhp_obj, mmhp_events)
```

### Visualize intensity and goodness of fit jointly

<!-- The main goal of this package is to provide tools to  -->
<!-- diagnose the fit of a given point process to data.  -->
<!-- To do this, we include several functions: -->

- `intensityqqplot` displays the estimated intensity of a given
point process along with a QQ-plot of the rescaled inter-event times. 
These together can often be useful in identifying issues with model fit for a chosen point process.

```{r intensityqqplot, fig.width=6, fig.height=5}
intensityqqplot(object = fit_hp, events = hp_events )
```



```{r intqqpot mmhp, eval=FALSE}
# this gives an error currently
intensityqqplot(object = mmhp_obj, markov_states = mmhp_events)
```


### Residual Analysis

- `pp_residual` returns both raw and Pearson residuals from fitting
the specified point process to the given events.

```{r mmhp_residual}
pp_residual(object = mmhp_obj, events = mmhp_events$events)

pp_residual(object = fit_hp, events = hp_events)

```


### Overall summary of fit

- Finally, `pp_diag` summarises (both graphically and numerically) the fit
of a specified point process to the data. For a given point process
it computes the residuals (both raw and Pearson) obtained from fitting
that point process to the data, performs a goodness of fit test
based on the rescaled inter-event times, and displays graphical
summaries of this diagnostic.

```{r ppdiag hp, fig.width = 6, fig.height=4}
pp_diag(object = fit_hp, events = hp_events)
```


<!-- -   Homogeneous Poisson Process -->

<!-- `pp_diag(object, events)` gives diagnostics of the model, including a qq -->
<!-- plot, a ks plot and the corresponding ks test, -->
<!-- along with raw and Pearson residuals in one function. -->

<!-- ```{r,fig.width=7,fig.height=4} -->
<!-- pp_diag(hpp_obj,hpp_events) -->
<!-- ``` -->

<!-- `pp_residual(object, events)` -->
<!-- returns only the raw and Pearson residuals. -->

<!-- ```{r} -->
<!-- pp_residual(hpp_obj,hpp_events) -->
<!-- ``` -->

<!-- `intensityqqplot(object, events)` gives both qqplot and intensity plot. -->

<!-- ```{r,fig.width=7,fig.height=4} -->
<!-- intensityqqplot(hpp_obj, hpp_events) -->
<!-- ``` -->

<!-- -   Hawkes Process -->

<!-- `pp_diag(object, events)` gives diagnostics of the model, including a qq -->
<!-- plot, a ks plot, ks test, raw and pearson residuals in one function. -->

<!-- ```{r,fig.width=7,fig.height=4} -->
<!-- pp_diag(hp_obj,hp_events) -->
<!-- ``` -->

<!-- `pp_residual(object, events)` gives raw and pearson residuals. -->

<!-- ```{r} -->
<!-- pp_residual(hp_obj,hp_events) -->
<!-- ``` -->

<!-- `intensityqqplot(object, events)` gives both qqplot and intensity plots. -->

<!-- ```{r,fig.width=7,fig.height=4} -->
<!-- intensityqqplot(hp_obj, hp_events) -->
<!-- ``` -->

<!-- ### MMHP -->

<!-- ```{r,fig.width=7,fig.height=4} -->
<!-- pp_diag(mmhp_obj,mmhp_events$events) -->
<!-- ``` -->
---
title: "Fitting Markov Modulated Point Process Models"
output: rmarkdown::html_vignette
author: Owen G. Ward.
vignette: >
  %\VignetteIndexEntry{Fitting Markov Modulated Point Process Models}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
  
```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = FALSE,
  comment = "#>"
)
set.seed(100) # to make it reproducible
options(rmarkdown.html_vignette.check_title = FALSE)
```

```{r setup}
library(ppdiag)
library(rstan)
library(cmdstanr)
library(tidyverse)
library(bayesplot)
```


Functions for fitting standard Hawkes and Poisson
point processes to data are included in `ppdiag`. However, currently, we do not
include fitting algorithms for Markov modulated point processes,
as these rely on the use of `rstan` or (more recently), `cmdstanr` to
perform Bayesian inference for the model parameters.

Here we provide instructions on how to fit these models so 
that they can easily be used in conjunction with the diagnostic tools
of  `ppdiag`.

# Installing RStan or Cmdstanr

- Instructions for installing RStan can be found at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started 
- Instructions for installing CmdStanR can be found at https://mc-stan.org/cmdstanr/articles/cmdstanr.html. As mentioned in this
link there are some differences between RStan and CmdStanR. In our experience,
CmdStanR does seem to be maybe more reliable and less prone to crashing.
- If you are having trouble installing or using Stan there is lots of great help
available at https://discourse.mc-stan.org/, where someone else has probably
already answered your question!


# The Stan code for MMPP/MMHP

We first include the stan code to fit each of these Markov modulated
models to a temporal point process.

```{r mmpp-stan-code}
mmpp_stan_code <- "
data{
  int<lower=1> num_events; //maximum of number of events for each pair each window => max(unlist(lapply(return_df$event.times,length))))
  vector[num_events+1] time_matrix; // include termination time in the last entry
}
parameters{
  real<lower=0> lambda0; //baseline rate for each pair
  real<lower=0> c; //baseline rate for each pair
  real<lower=0,upper=1> w1; //CTMC transition rate
  real<lower=0,upper=1> w2; //CTMC transition rate
}
transformed parameters{
  real<lower=0> q1;
  real<lower=0> q2;
  q1 = (lambda0).*w1;
  q2 = (lambda0).*w2;
}
model{
  real integ; // Placeholder variable for calculating integrals
  row_vector[2] forward[num_events]; // Forward variables from forward-backward algorithm
  row_vector[2] forward_termination; 
  row_vector[2] probs_1[num_events+1]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[num_events+1]; // Probability vector for transition to state 2 (inactive state)
  vector[num_events+1] interevent;
  
  //priors
  c ~ lognormal(0,1);
  w1 ~ beta(0.5,0.5);
  w2 ~ beta(0.5,0.5);
  
  lambda0 ~ gamma(1, 1);
  interevent = time_matrix;
  // ---- prepare for forward algorithm
  // --- log probability of Markov transition logP_ij(t)
  for(n in 1:(num_events + 1)){
    probs_1[n][1] = log(q2/(q1+q2)+q1/(q1+q2)*exp(-(q1+q2)*interevent[n])); //1->1
    probs_2[n][2] = log(q1/(q1+q2)+q2/(q1+q2)*exp(-(q1+q2)*interevent[n])); //2->2
    probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
    probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
  }
  //consider n = 1
  integ = interevent[1]*lambda0;
  forward[1][1] = log_sum_exp(probs_1[1]) + log(lambda0*(1+c)) - integ*(1+c); 
  forward[1][2] = log_sum_exp(probs_2[1]) + log(lambda0) - integ; 
  
  if(num_events>1){
    for(n in 2:num_events){
      integ = interevent[n]*lambda0;
      forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log(lambda0*(1+c))- integ*(1+c);
      forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(lambda0) - integ;
    }
  }
  
  integ = interevent[num_events]*lambda0;
  forward_termination[1] = log_sum_exp(forward[num_events] + probs_1[num_events]) - integ*(1+c);
  forward_termination[2] = log_sum_exp(forward[num_events] + probs_2[num_events]) - integ;
  
  target += log_sum_exp(forward_termination);
}
"

```

```{r mmhp-stan-code}
mmhp_stan_code <- "
data{
  int<lower=1> num_events; //number of events
  vector[num_events+1] time_matrix; // include termination time as last entry
}
parameters{
  real<lower=0> lambda0; //baseline rate for each pair
  real<lower=0> w_lambda;
  real<lower=0, upper=1> w_q1; //CTMC transition rate
  real<lower=0, upper=1> w_q2; //
  real<lower=0> alpha;
  real<lower=0> beta_delta;
  real<lower=0,upper=1> delta_1; // P(initial state = 1)
}
transformed parameters{
  real<lower=0> lambda1;
  real<lower=0> q1;
  real<lower=0> q2;
  row_vector[2] log_delta;
  real<lower=0> beta;
  lambda1 = (lambda0).*(1+w_lambda); 
  q2 = (lambda0).*w_q2;
  q1 = (lambda0).*w_q1;
  log_delta[1] = log(delta_1);
  log_delta[2] = log(1-delta_1);
  beta = alpha*(1+beta_delta);
}
model{
  real integ; // Placeholder variable for calculating integrals
  row_vector[2] forward[num_events]; // Forward variables from forward-backward algorithm
  row_vector[2] forward_termination; // Forward variables at termination time
  row_vector[2] probs_1[num_events+1]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[num_events+1]; // Probability vector for transition to state 2 (inactive state)
  row_vector[2] int_1[num_events+1]; // Integration of lambda when state transit to 1 (active state)
  row_vector[2] int_2[num_events+1]; // Integration of lambda when state transit to 2 (inactive state)
  real R[num_events+1]; // record variable for Hawkes process
  vector[num_events+1] interevent;
  real K0;
  real K1;
  real K2;
  real K3;
  real K4;
  real K5;
  //priors
  w_lambda ~ gamma(1,1);
  alpha ~ gamma(1,1);//lognormal(0,1);
  beta_delta ~ lognormal(0,2);//normal(0,10);
  //delta_1 ~ beta(2,2);
  w_q1 ~ beta(2,2);
  w_q2 ~ beta(2,2);
  

  lambda0 ~ gamma(1,1);
  interevent = time_matrix;
  if(num_events==0){ // there is no event occured in this period
    //--- prepare for forward calculation
    probs_1[1][1] = log(q2/(q1+q2)+q1/(q1+q2)*exp(-(q1+q2)*interevent[1])); //1->1
    probs_2[1][2] = log(q1/(q1+q2)+q2/(q1+q2)*exp(-(q1+q2)*interevent[1])); //2->2
    probs_1[1][2] = log1m_exp(probs_2[1][2]); //2->1
    probs_2[1][1] = log1m_exp(probs_1[1][1]); //1->2
    R[1] = 0;
    K0 = exp(-(q1+q2)*interevent[1]);
    K1 = (1-exp(-(q1+q2)*interevent[1]))/(q1+q2);
    K2 = (1-exp(-(q1+q2)*interevent[1]))/(q1+q2);
    int_1[1][1] = ((q2^2*lambda1+q2*q1*lambda0)*interevent[1] +
                     (q1^2*lambda1+q2*q1*lambda0)*K0*interevent[1] +
                     (lambda1-lambda0)*q2*q1*K1 + (lambda1-lambda0)*q2*q1*K2)/(q1+q2)^2/exp(probs_1[1][1]); //1->1
    int_1[1][2] = ((q2^2*lambda1+lambda0*q1*q2)*interevent[1] -
                     (lambda1*q1*q2+lambda0*q2^2)*K0*interevent[1] +
                     (lambda0-lambda1)*q2^2*K1 + (lambda1-lambda0)*q1*q2*K2)/(q1+q2)^2/exp(probs_1[1][2]); //2->1
    int_2[1][1] = ((q1*q2*lambda1+q1^2*lambda0)*interevent[1] -
                     (q1^2*lambda1+q1*q2*lambda0)*K0*interevent[1] +
                     (lambda1-lambda0)*q1^2*K1 + q1*q2*(lambda0-lambda1)*K2)/(q1+q2)^2/exp(probs_2[1][1]); //1->2
    int_2[1][2] = ((q1*q2*lambda1+lambda0*q1^2)*interevent[1] +
                     (q1*q2*lambda1+lambda0*q2^2)*K0*interevent[1] +
                     (lambda0-lambda1)*q1*q2*K1 + (lambda0-lambda1)*q1*q2*K2)/(q1+q2)^2/exp(probs_2[1][2]); //2->2
    
    forward_termination[1] = log_sum_exp(log_delta + probs_1[1] - int_1[1]);
    forward_termination[2] = log_sum_exp(log_delta + probs_2[1] - int_2[1]);
    target += log_sum_exp(forward_termination);
    //target += -lambda0*interevent[1]*delta_1-lambda1*interevent[1]*(1-delta_1);
  }else{ // there is event occured
    // ---- prepare for forward algorithm
    // --- log probability of Markov transition logP_ij(t)
    for(n in 1:(num_events + 1)){ //changed this
      probs_1[n][1] = log(q2/(q1+q2)+q1/(q1+q2)*exp(-(q1+q2)*interevent[n])); //1->1
      probs_2[n][2] = log(q1/(q1+q2)+q2/(q1+q2)*exp(-(q1+q2)*interevent[n])); //2->2
      probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
      probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
    }
    
    // --- R for Hawkes
    R[1] = 0;
    for(n in 2:(num_events + 1)){ // and this
      R[n] = exp(-beta*interevent[n])*(R[n-1]+1);
    }
    
    // Integration of lambda
    for(n in 1:(num_events)){ //and this
      K0 = exp(-(q1+q2)*interevent[n]);
      K1 = (1-exp(-(q1+q2)*interevent[n]))/(q1+q2);
      K2 = (1-exp(-(q1+q2)*interevent[n]))/(q1+q2);
      K3 = R[n]*(exp(beta*interevent[n])-1)/beta;
      K4 = R[n]*(1-exp(-(beta+q1+q2)*interevent[n]))*exp(beta*interevent[n])/(beta+q1+q2);
      K5 = R[n]*(1-exp(-(q1+q2-beta)*interevent[n]))/(q1+q2-beta);
      int_1[n][1] = ((q2^2*lambda1+q2*q1*lambda0)*interevent[n] +
                       (q1^2*lambda1+q2*q1*lambda0)*K0*interevent[n] +
                       (lambda1-lambda0)*q2*q1*K1 + (lambda1-lambda0)*q2*q1*K2 +
                       alpha*K3*(q2^2+q1^2*K0) +
                       alpha*q1*q2*K4 + alpha*q1*q2*K5)/(q1+q2)^2/exp(probs_1[n][1]); //1->1
      int_1[n][2] = ((q2^2*lambda1+lambda0*q1*q2)*interevent[n] -
                       (lambda1*q1*q2+lambda0*q2^2)*K0*interevent[n] +
                       (lambda0-lambda1)*q2^2*K1 + (lambda1-lambda0)*q1*q2*K2 +
                       alpha*q2*K3*(q2-q1*K0) -
                       alpha*q2^2*K4 + alpha*q1*q2*K5)/(q1+q2)^2/exp(probs_1[n][2]); //2->1
      int_2[n][1] = ((q1*q2*lambda1+q1^2*lambda0)*interevent[n] -
                       (q1^2*lambda1+q1*q2*lambda0)*K0*interevent[n] +
                       (lambda1-lambda0)*q1^2*K1 + q1*q2*(lambda0-lambda1)*K2 +
                       alpha*q1*K3*(q2-q1*K0) +
                       alpha*q1^2*K4 - alpha*q2*q1*K5)/(q1+q2)^2/exp(probs_2[n][1]); //1->2
      int_2[n][2] = ((q1*q2*lambda1+lambda0*q1^2)*interevent[n] +
                       (q1*q2*lambda1+lambda0*q2^2)*K0*interevent[n] +
                       (lambda0-lambda1)*q1*q2*K1 + (lambda0-lambda1)*q1*q2*K2 +
                       alpha*q1*q2*K3*(1+K0) -
                       alpha*q1*q2*K4 - alpha*q1*q2*K5)/(q1+q2)^2/exp(probs_2[n][2]); //2->2
    }
    
    //consider n = 1
    forward[1][1] = log(lambda1) + log_sum_exp(probs_1[1]-int_1[1]+log_delta); 
    forward[1][2] = log(lambda0) + log_sum_exp(probs_2[1]-int_2[1]+log_delta); 
    
    if(num_events>1){
      for(n in 2:num_events){
        forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n] - int_1[n]) + log(lambda1+alpha*R[n]);
        forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n] - int_2[n]) + log(lambda0);
      }
    }
    
    forward_termination[1] = log_sum_exp(forward[num_events] + probs_1[num_events] - int_1[num_events]);
    forward_termination[2] = log_sum_exp(forward[num_events] + probs_2[num_events] - int_2[num_events]);
    // lots of places with max_Nm and Nm got rid of the +1
    target += log_sum_exp(forward_termination);
  }
}
"

```


# Fitting these models in RStan

To demonstrate how to use these models, we will simulate some data
from each of these models and fit the included `stan` models.


```{r sim-mmpp}
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
mmpp_obj <- pp_mmpp(Q = Q, lambda0 = 1, c = 1.5, delta = c(1/3, 2/3))

sim_mmpp <- pp_simulate(mmpp_obj, n = 50)
```




We also show how to run these models using `rstan` for the MMPP 
example.

```{r rstan-mmpp}
mmpp_data <- list(num_events = length(sim_mmpp$events) - 1, 
                  # as first event is start time (0)
                  time_matrix = diff(c(sim_mmpp$events, sim_mmpp$end))
                  #interevent arrival time
                  )

mmpp_rstan <- stan(model_code = mmpp_stan_code, 
                   data = mmpp_data,
                   chains = 2)
mmpp_sim <- rstan::extract(mmpp_rstan)
```

The fit of this model can be evaluated in several ways, including
some of the tools in the `bayesplot` package, which are important to
ensure the fit is reasonable. 

```{r mmpp-plot}
bayesplot::mcmc_hist(as.matrix(mmpp_rstan), pars = c("lambda0", "c"))
```

To use the results of this fit with 
the functions of `ppdiag`, we want to extract an estimate
of each of the parameters. This can be achieved using the
posterior mean.

The desired posterior quantities can then be extracted from this object.

```{r mmpp-post}
mmpp_post <- lapply(mmpp_sim, mean)
```

Similarly, we can fit the MMHP model to simulated data also.

```{r sim-mmhp}
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
mmhp_obj <- pp_mmhp(Q = Q, lambda0 = 0.5, lambda1 = 1.5,
                    alpha = 0.5, beta = 0.75, delta = c(1/3, 2/3))

sim_mmhp <- pp_simulate(mmhp_obj, n = 50)


mmhp_data <- list(num_events = length(sim_mmhp$events) - 1, 
                  # as first event is start time (0)
                  time_matrix = diff(c(sim_mmhp$events, sim_mmhp$end))
                  #interevent arrival time
                  )
```

```{r rstan-mmhp}
mmhp_rstan <- stan(model_code = mmhp_stan_code, 
                   data = mmhp_data,
                   chains = 2)
mmhp_sim <- rstan::extract(mmhp_rstan)
mmhp_post <- lapply(mmhp_sim, mean)
```

```{r plot-draws-mmhp}
bayesplot::mcmc_hist(as.matrix(mmhp_rstan),
                     pars = c("lambda0", "lambda1", "alpha", "beta"))
```



# Fitting these models in Cmdstanr

We can also fit these models using `cmdstanr`. We include the code
to fit these examples below (which is not run here).

```{r mmpp-fit, eval=FALSE}
mmpp_file <- write_stan_file(mmpp_stan_code)
mmpp_stan <- cmdstan_model(stan_file = mmpp_file)

fit_mmpp <- mmpp_stan$sample(data = mmpp_data,
                             seed = 123,
                             chains = 4,
                             parallel_chains = 4,
                             refresh = 500)
```

Similarly, we can fit the MMHP model in `cmdstanr`.

```{r compile-mmhp, eval=FALSE}
mmhp_file <- write_stan_file(mmhp_stan_code)
mmhp_stan <- cmdstan_model(stan_file = mmhp_file)

fit_mmhp <- mmhp_stan$sample(data = mmhp_data,
                             seed = 123,
                             chains = 4,
                             parallel_chains = 4,
                             refresh = 500)
```


# Using the Stan output with `ppdiag`

Having obtained estimates from fitting these stan models, we
can now pass these parameter estimates into `ppdiag` to 
evaluate model fit.

```{r mmpp-ppdiag, fig.width=8}
mmpp_fit_obj <- pp_mmpp(lambda0 = mmpp_post$lambda0,
                       c = mmpp_post$c,
                       Q = matrix( c(-mmpp_post$q1,
                                     mmpp_post$q1,
                                     mmpp_post$q2,
                                     -mmpp_post$q2),
                                   nrow = 2, ncol = 2, byrow = T) )


pp_diag(mmpp_fit_obj, events = sim_mmpp$events)

```

```{r mmhp-ppdiag, fig.width=8}
mmhp_fit_obj <- pp_mmhp(lambda0 = mmhp_post$lambda0,
                        lambda1 = mmhp_post$lambda1,
                        alpha = mmhp_post$alpha,
                        beta = mmhp_post$beta,
                        Q = matrix( c(-mmhp_post$q1,
                                      mmhp_post$q1,
                                      mmhp_post$q2,
                                      -mmhp_post$q2),
                                 nrow = 2, ncol = 2, byrow = T) )

pp_diag(mmhp_fit_obj, events = sim_mmhp$events)

```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_diag.R
\name{pp_diag}
\alias{pp_diag}
\alias{pp_diag.default}
\alias{pp_diag.hp}
\alias{pp_diag.mmhp}
\alias{pp_diag.mmpp}
\alias{pp_diag.hpp}
\title{Summarise diagnostics for point process models}
\usage{
pp_diag(object, events)

\method{pp_diag}{default}(object, events)

\method{pp_diag}{hp}(object, events)

\method{pp_diag}{mmhp}(object, events)

\method{pp_diag}{mmpp}(object, events)

\method{pp_diag}{hpp}(object, events)
}
\arguments{
\item{object}{a point process model}

\item{events}{event times}
}
\value{
Invisibly returns NULL. Outputs plots and summary of diagnostics to
console
}
\description{
Generate diagnostic tools for different point process models,
including quantile-quantile plot, ks plot,
raw residual and pearson residual.
}
\examples{
hpp_obj <- pp_hpp(lambda = 1)
events <- pp_simulate(hpp_obj, end = 50)
pp_diag(hpp_obj, events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drawHPIntensity.R
\name{drawHPIntensity}
\alias{drawHPIntensity}
\title{Draw the intensity of Hawkes Process}
\usage{
drawHPIntensity(
  hp = NULL,
  events,
  int_title = "Hawkes Intensity",
  start = 0,
  end = max(events),
  history = NULL,
  color = 1,
  i = 1,
  add = FALSE,
  fit = FALSE,
  plot_events = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{hp}{object parameters for Hawkes process.}

\item{events}{the event times happened in this state}

\item{int_title}{title of the intensity plot}

\item{start}{the start time of current state}

\item{end}{the end time of current state}

\item{history}{the past event times}

\item{color}{specify the default plotting color.}

\item{i}{state number, used only for drawUniMMHPIntensity}

\item{add}{whether to add the hawkes intensity to an existing plot, used
for drawUniMMHPIntensity}

\item{fit}{a boolean indicating whether to fit a new HP to events}

\item{plot_events}{indicate whether events will be plotted}

\item{verbose}{whether to output informative messages as running}
}
\value{
no return value, intensity plot of Hawkes process
}
\description{
Draw the intensity of a Hawkes Process
}
\examples{
set.seed(100)
hp_obj <- pp_hp(lambda0 = 0.5, alpha = 0.45, beta = 0.5)
events <- pp_simulate(hp_obj, start = 0, end = 20)
drawHPIntensity(hp_obj, events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_residual.R
\name{pp_residual}
\alias{pp_residual}
\title{Compute raw and pearson residuals for point process models}
\usage{
pp_residual(object, events, start = 0, end = max(events), steps = 1000)
}
\arguments{
\item{object}{point process model containing the parameters}

\item{events}{vector of event times}

\item{start}{start of observation period (default 0)}

\item{end}{end of observation period (default final event)}

\item{steps}{number of steps for numeric integration (if needed)}
}
\value{
the raw and pearson residuals
}
\description{
Compute raw and pearson residuals for point process models
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
x <- pp_mmhp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9,
  lambda1 = 1.1, alpha = 0.8, beta = 1.2
)
y <- pp_simulate(x, n = 10)
pp_residual(x, events = y$events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_compensator.R
\name{pp_compensator}
\alias{pp_compensator}
\alias{pp_compensator.default}
\alias{pp_compensator.mmpp}
\alias{pp_compensator.hp}
\alias{pp_compensator.mmhp}
\alias{pp_compensator.hpp}
\title{Compensators for point processes}
\usage{
pp_compensator(object, events)

\method{pp_compensator}{default}(object, events)

\method{pp_compensator}{mmpp}(object, events)

\method{pp_compensator}{hp}(object, events)

\method{pp_compensator}{mmhp}(object, events)

\method{pp_compensator}{hpp}(object, events)
}
\arguments{
\item{object}{a point process model}

\item{events}{event times, which can have first value as 0}
}
\value{
compensator vector of rescaled interevent times
}
\description{
Computes the compensator for included point processes
}
\examples{
hpp_obj <- pp_hpp(lambda = 1)
events <- pp_simulate(hpp_obj, end = 10)
comp <- pp_compensator(hpp_obj, events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_mmhp.R
\name{pp_mmhp}
\alias{pp_mmhp}
\title{Create a Markov-modulated Hawkes Process(MMHP) object}
\usage{
pp_mmhp(lambda0, lambda1, alpha, beta, Q = NULL, delta = NULL, events = NULL)
}
\arguments{
\item{lambda0}{intensity for homogeneous Poisson process.}

\item{lambda1}{base intensity for Hawkes process.}

\item{alpha}{jump size of the increase in intensity in the hawkes process}

\item{beta}{exponential decrease of intensity in the hawkes process}

\item{Q}{transition probability matrix.}

\item{delta}{initial state probability.}

\item{events}{vector containing the event times.
Note that the first event is at time zero.
Alternatively, events could be specified as NULL,
 meaning that the data will be added later (e.g. simulated).}
}
\value{
mmhp object
}
\description{
Create a Markov-modulated Hawkes Process(MMHP) model
according to the given parameters: lambda0, lambda1,
alpha, beta, event times and transition probability matrix.
If event time events is missing,
 then it means that data will be added later(e.g. simulated)
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
pp_mmhp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9, lambda1 = 1.1,
  alpha = 0.8, beta = 1.2
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drawUniMMHPIntensity.R
\name{drawUniMMHPIntensity}
\alias{drawUniMMHPIntensity}
\title{Draw the intensity of the Markov-modulated Hawkes Process(MMHP)}
\usage{
drawUniMMHPIntensity(
  mmhp,
  simulation,
  int_title = "Intensity of MMHP",
  leg_location = "topright",
  color = 1,
  add = FALSE
)
}
\arguments{
\item{mmhp}{a mmhp object including its state, state_time,
events, lambda0, lambda1, beta and alpha.}

\item{simulation}{the simulated Markov-modulated Hawkes Process(MMHP)}

\item{int_title}{title of the plot.}

\item{leg_location}{location of legend, if moving needed}

\item{color}{A specification for the default plotting color.}

\item{add}{logical; if TRUE add to an already existing plot;
if NA start a new plot taking the defaults
for the limits and log-scaling of the x-axis from the previous plot.
 Taken as FALSE (with a warning if a different value is supplied)
 if no graphics device is open.}
}
\value{
no return value, intensity plot of Markov-modulated Hawkes process
}
\description{
Take a mmhp object and draw its intensity accordingly
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
x <- pp_mmhp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9, lambda1 = 1.1,
  alpha = 0.8, beta = 1.2
)
y <- pp_simulate(x, n = 25)
drawUniMMHPIntensity(x, y)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_mmpp.R
\name{pp_mmpp}
\alias{pp_mmpp}
\title{Create a Markov-modulated Poisson Process(MMPP) object}
\usage{
pp_mmpp(lambda0, c, Q, events = NULL, delta = NULL)
}
\arguments{
\item{lambda0}{parameters for Poisson process.}

\item{c}{the proportion of intensity 1 over intensity 2}

\item{Q}{transition probability matrix}

\item{events}{vector containing the event times.
Note that the first event is often specified as zero.
Alternatively, events could be specified as NULL,
 meaning that the data will be added later (e.g. simulated).}

\item{delta}{initial state probability.}
}
\value{
mmpp object
}
\description{
Create a Markov-modulated Poisson Process(MMPP) model
 according to the given parameters: lambda0, c, q1, q2 and event times.
If event time tau is missing,
then it means that data will be added later(e.g. simulated)
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
pp_mmpp(Q = Q, lambda0 = 1, c = 1.5, delta = c(1 / 3, 2 / 3))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_qqexp.R
\name{pp_qqexp}
\alias{pp_qqexp}
\title{Plot QQ-plot for rescaled-inter-event-times of fitted point process}
\usage{
pp_qqexp(r, ...)
}
\arguments{
\item{r}{rescaled-inter-event-times}

\item{...}{other arguments for plots}
}
\value{
no return value, quantile-quantile plot for rescaled-inter-event-times
}
\description{
Generate Quantile-quantile plot for rescaled-inter-event-times,
which are independently
and identically distributed as exponential random variables with rate 1
under the true point process.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmhp_event_state.R
\name{mmhp_event_state}
\alias{mmhp_event_state}
\title{Estimate the latent state for events of an MMHP}
\usage{
mmhp_event_state(
  params = list(lambda0, lambda1, alpha, beta, Q),
  events,
  start = 0
)
}
\arguments{
\item{params}{the parameters of the chosen MMHP}

\item{events}{the event times}

\item{start}{the start of observation}
}
\value{
probability of being in state 1 (active state) at each event,
along with most likely state
}
\description{
Given a MMHP object and events from that MMHP, infer the
most likely state of the Markov Process and each event time, along
with the probability of being in the active state.
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
mmhp_obj <- pp_mmhp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9,
  lambda1 = 1.1,
  alpha = 0.8, beta = 1.2
)
## evaluate at some fake event times
ppdiag:::mmhp_event_state(params = mmhp_obj, events = c(1, 2, 3, 5))
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmpp_event_state.R
\name{mmpp_event_state}
\alias{mmpp_event_state}
\title{Estimate the latent state for events of an MMPP}
\usage{
mmpp_event_state(params = list(lambda0, c, Q), events, start = 0)
}
\arguments{
\item{params}{the parameters of the chosen MMPP}

\item{events}{the event times}

\item{start}{the start of observation}
}
\value{
probability of being in state 1 (active state) at each event,
along with most likely state
}
\description{
Given a MMPP object and events from that MMPP, infer the
most likely state of the Markov Process at each event time, along
with the probability of being in the active state.
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
mmpp_obj <- pp_mmpp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9,
  c = 1.1
)
## evaluate at some fake event times
ppdiag:::mmpp_event_state(params = mmpp_obj, events = c(1, 2, 3, 5))
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_hpp.R
\name{pp_hpp}
\alias{pp_hpp}
\title{Create a homogeneous Poisson process object}
\usage{
pp_hpp(lambda, events = NULL)
}
\arguments{
\item{lambda}{rate of the Poisson process}

\item{events}{event times, optional}
}
\value{
hpp object
}
\description{
Create a homogeneous Poisson object according to given parameters:
lambda, and events.
If events are missing, then it means that data will be
added later(simulated from this process).
}
\examples{
pp_hpp(lambda = 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_ksplot.R
\name{pp_ksplot}
\alias{pp_ksplot}
\title{KS plot of empirical and theoretical cdf curve of fitted point process}
\usage{
pp_ksplot(r, ...)
}
\arguments{
\item{r}{rescaled-inter-event-times}

\item{...}{other arguments for plots}
}
\value{
no return value, KS plot for rescaled-inter-event-times and exponential cdf curve
}
\description{
Plot empirical cdf plot for rescaled-inter-event-times and
 exponential cdf as a reference curve
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rawresidual.R
\name{rawresidual}
\alias{rawresidual}
\title{Compute raw residuals for point process models}
\usage{
rawresidual(object, events, start, end, steps = 1000)
}
\arguments{
\item{object}{point process model containing the parameters}

\item{events}{vector of event times}

\item{start}{start of observation period (default 0)}

\item{end}{end of observation period (default final event)}

\item{steps}{number of steps for numeric integration (if needed)}
}
\value{
the raw residual
}
\description{
Compute raw residuals for for point processes
with specified parameters and events.
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
x <- pp_mmhp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9,
  lambda1 = 1.1, alpha = 0.8, beta = 1.2
)
y <- pp_simulate(x, n = 10)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fithp.R
\name{fithp}
\alias{fithp}
\title{Determine the MLE of Hawkes process numerically}
\usage{
fithp(events, end = max(events), vec = c(0.1, 0.2, 0.3))
}
\arguments{
\item{events}{event times}

\item{end}{end of observation period starting from 0 (default last event)}

\item{vec}{vector of initial parameter values}
}
\value{
a hp object indicating the maximum
likelihood parameter values (lambda0,alpha,beta) for Hawkes process.
This is a non-convex problem and a (unique) solution is not guaranteed.
}
\description{
Determine the MLE of Hawkes process numerically
}
\examples{
hp_obj <- pp_hp(lambda0 = 0.1, alpha = 0.45, beta = 0.5)
sims <- pp_simulate(hp_obj, start = 0, n = 10)
fithp(sims)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_simulate.R
\name{pp_simulate}
\alias{pp_simulate}
\alias{pp_simulate.default}
\alias{pp_simulate.hpp}
\alias{pp_simulate.hp}
\alias{pp_simulate.mmpp}
\alias{pp_simulate.mmhp}
\title{Simulate events from a temporal point process}
\usage{
pp_simulate(object, start = 0, end = 1, n = NULL, verbose = FALSE)

\method{pp_simulate}{default}(object, start = 0, end = 1, n = NULL, verbose = FALSE)

\method{pp_simulate}{hpp}(object, start = 0, end = 1, n = NULL, verbose = FALSE)

\method{pp_simulate}{hp}(object, start = 0, end = 1, n = NULL, verbose = FALSE)

\method{pp_simulate}{mmpp}(object, start = 0, end = 1, n = NULL, verbose = FALSE)

\method{pp_simulate}{mmhp}(object, start = 0, end = 1, n = NULL, verbose = FALSE)
}
\arguments{
\item{object}{point process model object of type hpp, hp, mmhp, or mmpp}

\item{start}{start time of events simulated. Not used for Markov modulated
models}

\item{end}{end time of events simulated. Not used for Markov modulated models}

\item{n}{number of events simulated. Required for Markov modulated models,
optional otherwise}

\item{verbose}{whether to output informative messages as running}
}
\value{
a vector of event times for all models. For Markov modulated models,
also returns details on the underlying latent process
}
\description{
Currently available point processes are homogeneous Poisson,
Hawkes with exponential kernel, MMHP and MMPP
}
\examples{
hpp_obj <- pp_hpp(lambda = 1)
s <- pp_simulate(hpp_obj, n = 50)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fithpp.R
\name{fithpp}
\alias{fithpp}
\title{Fit a homogeneous poisson process to event data}
\usage{
fithpp(events, end = max(events))
}
\arguments{
\item{events}{vector containing the event times.}

\item{end}{end of observation period, starting from 0 (default is last event)}
}
\value{
a hpp object containing the events and the estimated parameter
}
\description{
Compute maximum likelihood estimator of the rate of a homogeneous Poisson
process for the given events.
}
\examples{
pois_y <- pp_hpp(lambda = 1)
events <- pp_simulate(pois_y, end = 10)
fithpp(events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pearsonresidual.R
\name{pearsonresidual}
\alias{pearsonresidual}
\title{Compute Pearson residuals for point process models}
\usage{
pearsonresidual(object, events, start, end, steps = 1000)
}
\arguments{
\item{object}{point process model}

\item{events}{vector of event times}

\item{start}{start of observation period (default 0)}

\item{end}{termination time (default final event)}

\item{steps}{number of steps for numeric integration (if needed)}
}
\value{
the Pearson residual
}
\description{
Compute Pearson residuals for point processes
with specified parameters and events.
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
x <- pp_mmhp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9,
  lambda1 = 1.1, alpha = 0.8, beta = 1.2
)
y <- pp_simulate(x, n = 10)
ppdiag:::pearsonresidual.mmhp(x, events = y$events[-1])
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/negloglik.R
\name{negloglik}
\alias{negloglik}
\title{Compute negative log likelihood for point process models}
\usage{
negloglik(object, events, end)
}
\arguments{
\item{object}{point process object containing the parameters}

\item{events}{vector containing the event times.}

\item{end}{the end time of event times}
}
\value{
a scalar indicating the negative log likelihood
}
\description{
Compute negative log likelihood for point process models
 with model specified time events or simulated time events
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_hp.R
\name{pp_hp}
\alias{pp_hp}
\title{Create a Hawkes process object}
\usage{
pp_hp(lambda0, alpha, beta, events = NULL)
}
\arguments{
\item{lambda0}{initial intensity at the start time}

\item{alpha}{jump size in increase of intensity}

\item{beta}{exponential decay of intensity}

\item{events}{vector containing the event times.
Note that the first event is at time zero.
Alternatively, events could be specified as NULL,
meaning that the data will be added later (e.g. simulated).}
}
\value{
hp object
}
\description{
Create a Hawkes Process with an exponential
kernel according to the given parameters:
lambda0, alpha, beta and events.
If events are missing, then it means that data will be
added later(simulated from this process)
}
\examples{
pp_hp(lambda0 = 0.1, alpha = 0.45, beta = 0.5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intensityqqplot.R
\name{intensityqqplot}
\alias{intensityqqplot}
\alias{intensityqqplot.default}
\alias{intensityqqplot.hp}
\alias{intensityqqplot.hpp}
\alias{intensityqqplot.mmpp}
\alias{intensityqqplot.mmhp}
\title{Draw intensity of fitted point process and QQ-Plot of rescaled events}
\usage{
intensityqqplot(object, events, markov_states)

\method{intensityqqplot}{default}(object, events, markov_states)

\method{intensityqqplot}{hp}(object, events, markov_states = NULL)

\method{intensityqqplot}{hpp}(object, events, markov_states = NULL)

\method{intensityqqplot}{mmpp}(object, events = markov_states$events, markov_states)

\method{intensityqqplot}{mmhp}(object, events = markov_states$events, markov_states)
}
\arguments{
\item{object}{parameters for the models: hp, hpp, and mmhp}

\item{events}{event times}

\item{markov_states}{only for mmhp and mmpp, markov states simulation output}
}
\value{
no return value, intensity and qq-plot in a single plot
}
\description{
Draw the intensity and q-q plot for models
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fithp.R
\name{negloglik_hp}
\alias{negloglik_hp}
\title{Fit a Hawkes process with exponential kernel to event data}
\usage{
negloglik_hp(vec, events, end = max(events))
}
\arguments{
\item{vec}{vector containing initial values for the
object parameters (lambda0,alpha,beta) to be optimized.}

\item{events}{vector containing event times.}

\item{end}{the end time of event times.}
}
\description{
Compute the negative log likelihood parameter values for hawkes process.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hawkes_max_intensity.R
\name{hawkes_max_intensity}
\alias{hawkes_max_intensity}
\title{max intensity of a hawkes process}
\usage{
hawkes_max_intensity(object, events)
}
\arguments{
\item{object}{hawkes process object containing parameters for Hawkes process.}

\item{events}{events for Hawkes process.}
}
\value{
max of intensity
}
\description{
max intensity of a hawkes process
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drawUniMMPPIntensity.R
\name{drawUniMMPPIntensity}
\alias{drawUniMMPPIntensity}
\title{Draw the intensity of the Markov-modulated Poisson Process(MMPP)}
\usage{
drawUniMMPPIntensity(
  mmpp,
  simulation,
  add = FALSE,
  color = 1,
  fit = FALSE,
  int_title = "Intensity Plot of MMPP"
)
}
\arguments{
\item{mmpp}{a mmpp object including its transition probability matrix, lambda0, delta, and c.}

\item{simulation}{the simulated Markov-modulated Poisson Process(MMPP)}

\item{add}{logical; if TRUE add to an already existing plot;
if NA start a new plot taking the defaults
for the limits and log-scaling of the x-axis from the previous plot.
 Taken as FALSE (with a warning if a different value is supplied)
 if no graphics device is open.}

\item{color}{A specification for the default plotting color.}

\item{fit}{a boolean indicating whether to fit the events provided}

\item{int_title}{title of the plot.}
}
\value{
no return value, intensity plot of Markov-modulated Poisson process
}
\description{
Take a mmpp object and draw its intensity accordingly
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
x <- pp_mmpp(Q, delta = c(1 / 3, 2 / 3), lambda0 = 0.9, c = 1.2)
y <- pp_simulate(x, n = 10)
drawUniMMPPIntensity(x, y)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolate_mmhp_latent.R
\name{interpolate_mmhp_latent}
\alias{interpolate_mmhp_latent}
\title{Interpolate latent process of MMHP}
\usage{
interpolate_mmhp_latent(params, events, zt)
}
\arguments{
\item{params}{parameters of the MMHP, an MMHP object}

\item{events}{events (not including 0, but assumes start at 0)}

\item{zt}{inferred latent state of events}
}
\value{
list of the states of the Markov process (z.hat), including
initial state
and the times of the transitions between these times (x.hat).
}
\description{
Interpolate latent process of MMHP
}
\examples{
Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2, byrow = TRUE)
mmhp_obj <- pp_mmhp(Q,
  delta = c(1 / 3, 2 / 3), lambda0 = 0.9, lambda1 = 1.1,
  alpha = 0.8, beta = 1.2
)
ppdiag:::interpolate_mmhp_latent(
  params = mmhp_obj,
  events = c(1, 2, 3, 5), zt = c(2, 1, 1, 2)
)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drawHPPIntensity.R
\name{drawHPPIntensity}
\alias{drawHPPIntensity}
\title{Draw intensity of homogeneous Poisson process}
\usage{
drawHPPIntensity(
  hpp = NULL,
  events,
  int_title = "Homogeneous Poisson Process",
  start = 0,
  end = max(events),
  color = "red",
  plot_events = TRUE,
  fit = FALSE,
  add = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{hpp}{object for homogeneous Poisson process}

\item{events}{event times input}

\item{int_title}{the plot title}

\item{start}{start of events}

\item{end}{end of events}

\item{color}{a specification for the default plotting color.}

\item{plot_events}{a boolean indicating whether input events will be plotted}

\item{fit}{a boolean indicating whether to fit a hpp or
use the passed object}

\item{add}{whether to add the hpp intensity to an existing plot}

\item{verbose}{whether to output informative messages as running}
}
\value{
no return value, intensity plot of homogeneous Poisson process
}
\description{
Draw the intensity for a homogeneous Poisson process
}
\examples{
pois_y <- pp_hpp(lambda = 1)
drawHPPIntensity(pois_y, events = pp_simulate(pois_y, end = 10))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmpp_latent.R
\name{mmpp_latent}
\alias{mmpp_latent}
\title{Estimate Latent process of MMPP}
\usage{
mmpp_latent(params = list(lambda0, c, Q), events, zt, start = 0)
}
\arguments{
\item{params}{parameters of the MMPP, an MMPP object}

\item{events}{events (not including 0, but assumes start at 0)}

\item{zt}{inferred latent state of events}
}
\value{
list of the states of the Markov process (z.hat)
and the times of the transitions between these times (x.hat).
}
\description{
Estimate Latent process of MMPP
}
\examples{
Q <- matrix(c(-0.04, 0.04, 0.02, -0.02), ncol = 2, byrow = TRUE)
mmpp_obj <- pp_mmpp(Q, delta = c(1 / 3, 2 / 3), lambda0 = 1, c = 1)
ppdiag:::mmpp_latent(params = mmpp_obj, events = c(1, 2, 3), zt = c(2, 1, 1))
}
\keyword{internal}
