
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`hal9001`

[![R-CMD-check](https://github.com/tlverse/hal9001/workflows/R-CMD-check/badge.svg)](https://github.com/tlverse/hal9001/actions)
[![Coverage
Status](https://codecov.io/gh/tlverse/hal9001/branch/master/graph/badge.svg)](https://app.codecov.io/gh/tlverse/hal9001)
[![CRAN](https://www.r-pkg.org/badges/version/hal9001)](https://www.r-pkg.org/pkg/hal9001)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/hal9001)](https://CRAN.R-project.org/package=hal9001)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/hal9001)](https://CRAN.R-project.org/package=hal9001)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3558313.svg)](https://doi.org/10.5281/zenodo.3558313)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02526/status.svg)](https://doi.org/10.21105/joss.02526)

> The *Scalable* Highly Adaptive Lasso

**Authors:** [Jeremy Coyle](https://github.com/tlverse), [Nima
Hejazi](https://nimahejazi.org), [Rachael
Phillips](https://github.com/rachaelvp), [Lars van der
Laan](https://github.com/Larsvanderlaan), and [Mark van der
Laan](https://vanderlaan-lab.org/)

-----

## What’s `hal9001`?

`hal9001` is an R package providing an implementation of the scalable
*highly adaptive lasso* (HAL), a nonparametric regression estimator that
applies L1-regularized lasso regression to a design matrix composed of
indicator functions corresponding to the support of the functional over
a set of covariates and interactions thereof. HAL regression allows for
arbitrarily complex functional forms to be estimated at fast
(near-parametric) convergence rates under only global smoothness
assumptions (van der Laan 2017a; Bibaut and van der Laan 2019). For
detailed theoretical discussions of the highly adaptive lasso estimator,
consider consulting, for example, van der Laan (2017a), van der Laan
(2017b), and van der Laan and Bibaut (2017). For a computational
demonstration of the versatility of HAL regression, see Benkeser and van
der Laan (2016). Recent theoretical works have demonstrated success in
building efficient estimators of complex parameters when particular
variations of HAL regression are used to estimate nuisance parameters
(e.g., van der Laan, Benkeser, and Cai 2019; Ertefaie, Hejazi, and van
der Laan 2020).

-----

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=hal9001) via

``` r
install.packages("hal9001")
```

To contribute, install the *development version* of `hal9001` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("tlverse/hal9001")
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tlverse/hal9001/issues).

-----

## Example

Consider the following minimal example in using `hal9001` to generate
predictions via Highly Adaptive Lasso regression:

``` r
# load the package and set a seed
library(hal9001)
#> Loading required package: Rcpp
#> hal9001 v0.4.2: The Scalable Highly Adaptive Lasso
#> note: fit_hal defaults have changed. See ?fit_hal for details
set.seed(385971)

# simulate data
n <- 100
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

# fit the HAL regression
hal_fit <- fit_hal(X = x, Y = y, yolo = TRUE)
#> [1] "I'm sorry, Dave. I'm afraid I can't do that."
hal_fit$times
#>                   user.self sys.self elapsed user.child sys.child
#> enumerate_basis       0.016    0.000   0.017          0         0
#> design_matrix         0.005    0.000   0.005          0         0
#> reduce_basis          0.000    0.000   0.000          0         0
#> remove_duplicates     0.000    0.000   0.000          0         0
#> lasso                 3.767    0.012   3.782          0         0
#> total                 3.789    0.012   3.805          0         0

# training sample prediction
preds <- predict(hal_fit, new_data = x)
mean(hal_mse <- (preds - y)^2)
#> [1] 0.03754093
```

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/tlverse/hal9001/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `hal9001` R package, please cite both of the following:

``` 
    @software{coyle2022hal9001-rpkg,
      author = {Coyle, Jeremy R and Hejazi, Nima S and Phillips, Rachael V
        and {van der Laan}, Lars and {van der Laan}, Mark J},
      title = {{hal9001}: The scalable highly adaptive lasso},
      year  = {2022},
      url = {https://doi.org/10.5281/zenodo.3558313},
      doi = {10.5281/zenodo.3558313}
      note = {{R} package version 0.4.2}
    }

    @article{hejazi2020hal9001-joss,
      author = {Hejazi, Nima S and Coyle, Jeremy R and {van der Laan}, Mark
        J},
      title = {{hal9001}: Scalable highly adaptive lasso regression in
        {R}},
      year  = {2020},
      url = {https://doi.org/10.21105/joss.02526},
      doi = {10.21105/joss.02526},
      journal = {Journal of Open Source Software},
      publisher = {The Open Journal}
    }
```

-----

## License

© 2017-2022 [Jeremy R. Coyle](https://github.com/tlverse) & [Nima S.
Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

-----

## References

<div id="refs" class="references">

<div id="ref-benkeser2016hal">

Benkeser, David, and Mark J van der Laan. 2016. “The Highly Adaptive
Lasso Estimator.” In *2016 IEEE International Conference on Data Science
and Advanced Analytics (DSAA)*. IEEE.
<https://doi.org/10.1109/dsaa.2016.93>.

</div>

<div id="ref-bibaut2019fast">

Bibaut, Aurélien F, and Mark J van der Laan. 2019. “Fast Rates for
Empirical Risk Minimization over Càdlàg Functions with Bounded Sectional
Variation Norm.” <https://arxiv.org/abs/1907.09244>.

</div>

<div id="ref-ertefaie2020nonparametric">

Ertefaie, Ashkan, Nima S Hejazi, and Mark J van der Laan. 2020.
“Nonparametric Inverse Probability Weighted Estimators Based on the
Highly Adaptive Lasso.” <https://arxiv.org/abs/2005.11303>.

</div>

<div id="ref-vdl2017generally">

van der Laan, Mark J. 2017a. “A Generally Efficient Targeted Minimum
Loss Based Estimator Based on the Highly Adaptive Lasso.” *The
International Journal of Biostatistics*.
<https://doi.org/10.1515/ijb-2015-0097>.

</div>

<div id="ref-vdl2017finite">

———. 2017b. “Finite Sample Inference for Targeted Learning.”
<https://arxiv.org/abs/1708.09502>.

</div>

<div id="ref-vdl2019efficient">

van der Laan, Mark J, David Benkeser, and Weixin Cai. 2019. “Efficient
Estimation of Pathwise Differentiable Target Parameters with the
Undersmoothed Highly Adaptive Lasso.”
<https://arxiv.org/abs/1908.05607>.

</div>

<div id="ref-vdl2017uniform">

van der Laan, Mark J, and Aurélien F Bibaut. 2017. “Uniform Consistency
of the Highly Adaptive Lasso Estimator of Infinite-Dimensional
Parameters.” <https://arxiv.org/abs/1709.06256>.

</div>

</div>
# hal9001 0.4.2
* Version bump for CRAN resubmission following archiving.

# hal9001 0.4.1
* Minor adjustments to speed up unit tests and examples.
* Version bump for CRAN resubmission.

# hal9001 0.4.0

As of September 2021:
* Minor change to how binning is performed when `num_knots = 1`, ensuring that
  the minimal number of knots is chosen when `num_knots = 1`. This results in
  HAL agreeing with (main terms) `glmnet` when `smoothness_orders = 1` and
  `num_knots = 1`.
* Revised formula interface with enhanced capabilities, allowing specifciation
  of penalization factors, smoothness_orders, and the number of knots for each
  variable, for every single term separately using the new `h` function. It is
  possible to specify, e.g., `h(X) + h(W)` which will generate and concatenate
  the two basis function terms.

As of April 2021:
* The default of `fit_hal` is now a first order smoothed HAL with binning.
* Updated documentation for `formula_hal`, `fit_hal` and `predict`; and
  added `fit_control` and `formula_control` lists for arguments. Moved much of
  the text to details sections, and shortened the argument descriptions.
* Updated `summary` to support higher-order HAL fit interpretations.
* Added checks to `fit_hal` for missingness and dimensionality correspondence
  between `X`, `Y`, and `X_unpenalized`. These checks lead to quickly-produced
  errors, opposed to enumerating the basis list and then letting `glmnet` error
  on something trivial like this.
* Modified formula interface in `fit_hal`, so `formula` is now provided
  directly to `fit_hal` and `formula_hal` is run within `fit_hal`. Due to these
  changes, it no longer made sense for `formula_hal` to accept `data`, so it
  now takes as input `X`. Also, the `formula_fit_hal` function was removed as
  it is no longer needed.
* Support for the custom lasso procedure implemented in `Rcpp` has been
  discontinued. Accordingly, the `"lassi"` option and argument `fit_type` have
  been removed from `fit_hal`.
* Re-added `lambda.min.ratio` as a `fit_control` argument to `fit_hal`. We've
  seen that not setting `lambda.min.ratio` in `glmnet` can lead to no `lambda`
  values that fit the data sufficiently well, so it seems appropriate to
  override the `glmnet` default.

# hal9001 0.3.0

As of February 2021:
* Support _higher order_ HAL via the new `smoothness_orders` argument
   * `smoothness_orders` is a vector of length 1 or length `ncol(X)`.
  * If `smoothness_orders` is of length 1 then its values are recycled to form
      a vector of length `ncol(X)`.
  * Given such a vector of length `ncol(X)`, the ith element gives the level of
    smoothness for the variable corresponding to the ith column in `X`.
* Degree-dependant binning. Higher order terms are binned more coarsely; the
  `num_knots` argument is a vector up to `max_degree` controlling the
  degree-specific binning.
* Adds `formula_hal` which allows a formula specification of a HAL model.

# hal9001 0.2.8

As of November 2020:
* Allow support for Poisson family to `glmnet()`.
* Begins consideration of supporting arbitrary `stats::family()` objects to be
  passed through to calls to `glmnet()`.
* Simplifies output of `fit_hal()` by unifying the redundant `hal_lasso` and
  `glmnet_lasso` slots into the new `lasso_fit` slot.
* Cleans up of methods throughout and improves documentation, reducing a few
  redundancies for cleaner/simpler code in `summary.hal9001`.
* Adds link to DOI of the published _Journal of Open Source Software_ paper in
  `DESCRIPTION`.

# hal9001 0.2.7

As of September 2020:
* Adds a `summary` method for interpreting HAL regressions
  (https://github.com/tlverse/hal9001/pull/64).
* Adds a software paper for publication in the _Journal of Open Source
  Software_ (https://github.com/tlverse/hal9001/pull/71).

# hal9001 0.2.6

As of June 2020:
* Address bugs/inconsistencies reported in the prediction method when trying to
  specify a value of lambda not included in initial fitting.
* Addresses a bug arising from a silent failure in `glmnet` in which it ignores
  the argument `lambda.min.ratio` when `family = "gaussian"` is not set.
* Adds a short software paper for submission to JOSS.
* Minor documentation updates.

# hal9001 0.2.5

As of March 2020
* First CRAN release.
# Contributing to `hal9001` development

We, the authors of the `hal9001` R package, use the same guide as is used for
contributing to the development of the popular `tidyverse` ecosystem of R
packages. This document is simply a formal re-statement of that fact.

The goal of this guide is to help you get up and contributing to `hal9001` as
quickly as possible. The guide is divided into two main pieces:

* Filing a bug report or feature request in an issue.
* Suggesting a change via a pull request.

## Issues

When filing an issue, the most important thing is to include a minimal
reproducible example so that we can quickly verify the problem, and then figure
out how to fix it. There are three things you need to include to make your
example reproducible: required packages, data, code.

1.  **Packages** should be loaded at the top of the script, so it's easy to
    see which ones the example needs.

2.  The easiest way to include **data** is to use `dput()` to generate the R
    code to recreate it.

3.  Spend a little bit of time ensuring that your **code** is easy for others to
    read:

    * make sure you've used spaces and your variable names are concise, but
      informative

    * use comments to indicate where your problem lies

    * do your best to remove everything that is not related to the problem.
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a
fresh R session and pasting your script in.

(Unless you've been specifically asked for it, please don't include the output
of `sessionInfo()`.)

## Pull requests

To contribute a change to `hal9001`, you follow these steps:

1. Create a branch in git and make your changes.
2. Push branch to GitHub and issue pull request (PR).
3. Discuss the pull request.
4. Iterate until either we accept the PR or decide that it's not a good fit for
   `hal9001`.

Each of these steps are described in more detail below. This might feel
overwhelming the first time you get set up, but it gets easier with practice.

If you're not familiar with git or GitHub, please start by reading
<http://r-pkgs.had.co.nz/git.html>

Pull requests will be evaluated against a checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivates the
   need for change. Please describe the problem your PR addresses and show
   how your pull request solves it as concisely as possible.

   Also include this motivation in `NEWS` so that when a new release of
   `hal9001` comes out it's easy for users to see what's changed. Add your
   item at the top of the file and use markdown for formatting. The
   news item should end with `(@yourGithubUsername, #the_issue_number)`.

2.  __Only related changes__. Before you submit your pull request, please
    check to make sure that you haven't accidentally included any unrelated
    changes. These make it harder to see exactly what's changed, and to
    evaluate any unexpected side effects.

    Each PR corresponds to a git branch, so if you expect to submit
    multiple changes make sure to create multiple branches. If you have
    multiple changes that depend on each other, start with the first one
    and don't submit any others until the first one has been processed.

3.  __Use `hal9001` coding style__. To do so, please follow the [official
    `tidyverse` style guide](http://style.tidyverse.org). Maintaining a
    consistent style across the whole code base makes it much easier to jump
    into the code. If you're modifying existing `hal9001` code that doesn't
    follow the style guide, a separate pull request to fix the style would be
    greatly appreciated.

4.  If you're adding new parameters or a new function, you'll also need
    to document them with [`roxygen2`](https://github.com/klutometis/roxygen).
    Make sure to re-run `devtools::document()` on the code before submitting.

This seems like a lot of work but don't worry if your pull request isn't
perfect. It's a learning process. A pull request is a process, and unless
you've submitted a few in the past it's unlikely that your pull request will be
accepted as is. Please don't submit pull requests that change existing
behaviour. Instead, think about how you can add a new feature in a minimally
invasive way.

## Test environments
* ubuntu 20.04 (local + GitHub Actions), R 4.1.1
* macOS 10.15 (local + GitHub Actions), R 4.1.1
* windows 2019 (on GitHub Actions), R 4.1.1

## R CMD check results
There were no ERRORs or WARNINGs.
* There was 1 NOTE:
    installed size is 8.2Mb
      sub-directories of 1Mb or more:
         libs   8.0Mb

## Downstream dependencies
* None at present. There were two (`haldensify`, `txshift`) that were also
  removed when this package was archived.

## Resubmission
* This is a resubmission of a package removed from CRAN due to a Solaris build
  failure, though that OS is no longer tested against. There were no issues on
  any other OS's.
---
title: "`hal9001`: Scalable highly adaptive lasso regression in `R`"
tags:
  - machine learning
  - targeted learning
  - causal inference
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1, 2, 4
  - name: Jeremy R. Coyle
    orcid: 0000-0002-9874-6649
    affiliation: 2
  - name: Mark J. van der Laan
    orcid: 0000-0003-1432-5511
    affiliation: 2, 3, 4
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Division of Biostatistics, School of Public Health, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Center for Computational Biology, University of California, Berkeley
    index: 4
date: 25 September 2020
bibliography: refs.bib
---

# Summary

The `hal9001` `R` package provides a computationally efficient implementation of
the _highly adaptive lasso_ (HAL), a flexible nonparametric regression and
machine learning algorithm endowed with several theoretically convenient
properties. `hal9001` pairs an implementation of this estimator with an array of
practical variable selection tools and sensible defaults in order to improve the
scalability of the algorithm. By building on existing `R` packages for lasso
regression and leveraging compiled code in key internal functions, the `hal9001`
`R` package provides a family of highly adaptive lasso estimators suitable for
use in both modern large-scale data analysis and cutting-edge research efforts
at the intersection of statistics and machine learning, including the emerging
subfield of computational causal inference [@wong2020computational].

# Background

The highly adaptive lasso (HAL) is a nonparametric regression function capable
of estimating complex (e.g., possibly infinite-dimensional) functional
parameters at a fast $n^{-1/3}$ rate under only relatively mild conditions
[@vdl2017generally; @vdl2017uniform; @bibaut2019fast]. HAL requires that the
space of the functional parameter be a subset of the set of càdlàg (right-hand
continuous with left-hand limits) functions with sectional variation norm
bounded by a constant. In contrast to the wealth of data adaptive regression
techniques that make strong local smoothness assumptions on the true form of the
target functional, HAL regression's assumption of a finite sectional variation
norm constitutes only a _global_ smoothness assumption, making it a powerful and
versatile approach. The `hal9001` package primarily implements a zeroth-order
HAL estimator, which constructs and selects by lasso penalization a linear
combination of indicator basis functions, minimizing the loss-specific empirical
risk under the constraint that the $L_1$-norm of the resultant vector of
coefficients be bounded by a finite constant. Importantly, the estimator is
formulated such that this finite constant is the sectional variation norm of the
target functional.

Intuitively, construction of a HAL estimator proceeds in two steps. First,
a design matrix composed of basis functions is generated based on the available
set of covariates. The zeroth-order HAL makes use of indicator basis functions,
resulting in a large, sparse matrix with binary entries; higher-order HAL
estimators, which replace the use of indicator basis functions with splines,
have been formulated, with implementation in a nascent stage. Representation of
the target functional $f$ in terms of indicator basis functions partitions the
support of $f$ into knot points, with such basis functions placed over subsets
of sections of $f$. Generally, numerous basis functions are created, with an
appropriate set of indicator bases then selected through lasso penalization.
Thus, the second step of fitting a HAL model is performing $L_1$-penalized
regression on the large, sparse design matrix of indicator bases. The selected
HAL regression model approximates the sectional variation norm of the target
functional as the absolute sum of the estimated coefficients of indicator basis
functions. The $L_1$ penalization parameter $\lambda$ can be data adaptively
chosen via a cross-validation selector [@vdl2003unified; @vdv2006oracle];
however, alternative selection criteria may be more appropriate when the
estimand functional is not the target parameter but instead a nuisance function
of a possibly complex parameter [e.g., @vdl2019efficient;
@ertefaie2020nonparametric]. An extensive set of simulation experiments were
used to assess the prediction performance of HAL regression [@benkeser2016hal];
these studies relied upon the subsequently deprecated [`halplus` `R`
package](https://github.com/benkeser/halplus).

# `hal9001`'s core functionality

The `hal9001` package, for the `R` language and environment for statistical
computing [@R], aims to provide a scalable implementation of the HAL
nonparametric regression function. To provide a single, unified interface, the
principal user-facing function is `fit_hal()`, which, at minimum, requires
a matrix of predictors `X` and an outcome `Y`. By default, invocation of
`fit_hal()` will build a HAL model using indicator basis functions for up to
a limited number of interactions of the variables in `X`, fitting the penalized
regression model via the lasso procedure available in the extremely popular
`glmnet` `R` package [@friedman2009glmnet]. As creation of the design matrix of
indicator basis functions can be computationally expensive, several utility
functions (e.g., `make_design_matrix()`, `make_basis_list()`, `make_copy_map()`)
have been written in C++ and integrated into the package via the `Rcpp`
framework [@eddelbuettel2011rcpp; @eddelbuettel2013seamless]. `hal9001`
additionally supports the fitting of standard (Gaussian), logistic, and Cox
proportional hazards models (via the `family` argument), including variations
that accommodate offsets (via the `offset` argument) and partially penalized
models (via the `X_unpenalized` argument).

Over several years of development and usage, it was found that the performance
of HAL regression can suffer in high-dimensional settings. To alleviate these
computational limitations, several screening and filtering approaches were
investigated and implemented. These include screening of variables prior to
creating the design matrix and filtering of indicator basis functions (via the
`reduce_basis` argument) as well as early stopping when fitting the sequence of
HAL models in the $L_1$-norm penalization parameter $\lambda$. Future software
development efforts will continue to improve upon the computational aspects and
performance of the HAL regression options supported by `hal9001`. Currently,
stable releases of the `hal9001` package are made available on the Comprehensive
`R` Archive Network at https://CRAN.R-project.org/package=hal9001, while both
stable (branch `master`) and development (branch `devel`) versions of the
package are hosted at https://github.com/tlverse/hal9001. Releases of the
package use both GitHub and Zenodo (https://doit.org/10.5281/zenodo.3558313).

# Applications

As `hal9001` is the canonical implementation of the highly adaptive lasso, the
package has been relied upon in a variety of statistical applications. Speaking
generally, HAL regression is often used in order to develop efficient estimation
strategies in challenging estimation and inference problems; thus, we interpret
_statistical applications_ of HAL regression chiefly as examples of novel
theoretical developments that have been thoroughly investigated in simulation
experiments and with illustrative data analysis examples. In the sequel, we
briefly point out a few recently successful examples:

* @ju2020robust formulate a procedure based on HAL regression that allows the
  construction of asymptotically normal and efficient estimators of causal
  effects that are robust to the presence of instrumental variables, which can
  often lead to severe issues for estimation and inference [@hernan2020causal].
  While a variety of procedures have been proposed to overcome the issues posed
  by instrumental variables, a particularly successful idea was given by
  @shortreed2017outcome, who proposed standard lasso regression to select
  covariates for the exposure model based on an estimated outcome model. The
  work of @ju2020robust replaces the standard lasso with HAL regression,
  effectively screening for _infinitesimal instrumental basis functions_ rather
  than instrumental variables, providing much enhanced flexibility. Here, the
  authors demonstrate how HAL regression provides exceptionally fine-grained
  control over screening problematic covariates while simultaneously
  facilitating the construction of causal effect estimators with desirable
  asymptotic properties.
* @diaz2020causal introduce novel mediation effects based on joint stochastic
  interventions on exposure and mediator variables. To complement the new causal
  effects outlined in their work, these authors introduce efficient estimators
  that rely upon a fast rate of convergence of nuisance parameter estimators to
  their true counterparts. As the authors note, HAL is currently the only
  machine learning algorithm for which such a fast rate of convergence can
  rigorously be proven under minimal global smoothness assumptions. By relying
  upon HAL regression for the construction of their proposed estimators,
  @diaz2020causal advance not only the state-of-the-art in causal mediation
  analysis but also provide evidence, in both simulation experiments and an
  illustrative data analysis, of how HAL regression can be brought to bear on
  challenging causal inference problems to develop flexible and robust
  estimation strategies.
* @hejazi2020efficient develop novel theoretical insights for building efficient
  estimators of causal effects under two-phase sampling designs, relying upon
  the flexibility and fast convergence of HAL regression at the core of their
  theoretical contributions. Corrections for two-phase sampling, a family of
  procedures for developing efficient estimators of full-sample effects in spite
  of censoring introduced by the second-phase subsample, have received much
  attention, though developments applicable to large, unrestricted statistical
  models have been limited. These authors provide a formulation and theory for
  utilizing causal effect estimators, based on data subject to two-phase
  sampling, that attain asymptotic efficiency by way of the fast convergence
  rate of HAL regression. In effect, this works demonstrates that HAL regression
  has properties suitable for both flexible estimation and efficient inference
  in settings with complex data structures. The authors make their methodology
  available in the `txshift` R package [@hejazi2020txshift-rpkg;
  @hejazi2020txshift-joss], which relies upon `hal9001`. These authors
  additionally provide examples in simulation experiments and a re-analysis of
  a recent HIV-1 vaccine efficacy trial using their proposed statistical
  approach.
* @ertefaie2020nonparametric provide a considered study of the construction of
  inverse probability weighted (IPW) estimators that rely upon HAL regression in
  the estimation of the exposure mechanism. While IPW estimators classically
  require the use of parametric models of the exposure mechanism, these authors
  propose and investigate novel variants of these estimators that instead rely
  upon the fast convergence rate of HAL regression for the required nuisance
  parameter functional. In particular, @ertefaie2020nonparametric show through
  theoretical advances, several simulation experiments, and an illustrative data
  analysis of data from the well-documented NHEFS study that IPW estimators
  based on HAL regression can be made asymptotically linear and even efficient
  under an undersmoothing-based debiasing procedure. In so doing, the authors
  simultaneously advance the literatures on HAL regression and on IPW
  estimation, establishing the interface between the two as an area of viable
  future research. Notably, in demonstrating their proposed IPW estimators with
  the NHEFS data, the authors show that IPW estimators based on HAL regression
  can yield meaningful substantive conclusions without the typically restrictive
  parametric assumptions required for IPW estimation.

As further theoretical advances continue to be made with HAL regression, and the
resultant statistical methodology explored, we expect both the number and
variety of such examples to steadily increase.

# References

