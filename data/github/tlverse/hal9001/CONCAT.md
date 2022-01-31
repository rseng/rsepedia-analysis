
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

---
output:
  rmarkdown::github_document
bibliography: "inst/REFERENCES.bib"
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# R/`hal9001`

[![R-CMD-check](https://github.com/tlverse/hal9001/workflows/R-CMD-check/badge.svg)](https://github.com/tlverse/hal9001/actions)
[![Coverage Status](https://codecov.io/gh/tlverse/hal9001/branch/master/graph/badge.svg)](https://app.codecov.io/gh/tlverse/hal9001)
[![CRAN](https://www.r-pkg.org/badges/version/hal9001)](https://www.r-pkg.org/pkg/hal9001)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/hal9001)](https://CRAN.R-project.org/package=hal9001)
[![CRAN total downloads](http://cranlogs.r-pkg.org/badges/grand-total/hal9001)](https://CRAN.R-project.org/package=hal9001)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3558313.svg)](https://doi.org/10.5281/zenodo.3558313)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02526/status.svg)](https://doi.org/10.21105/joss.02526)

> The _Scalable_ Highly Adaptive Lasso

__Authors:__ [Jeremy Coyle](https://github.com/tlverse), [Nima
Hejazi](https://nimahejazi.org), [Rachael
Phillips](https://github.com/rachaelvp), [Lars van der
Laan](https://github.com/Larsvanderlaan), and [Mark van der
Laan](https://vanderlaan-lab.org/)

---

## What's `hal9001`?

`hal9001` is an R package providing an implementation of the scalable _highly
adaptive lasso_ (HAL), a nonparametric regression estimator that applies
L1-regularized lasso regression to a design matrix composed of indicator
functions corresponding to the support of the functional over a set of
covariates and interactions thereof. HAL regression allows for arbitrarily
complex functional forms to be estimated at fast (near-parametric) convergence
rates under only global smoothness assumptions [@vdl2017generally;
@bibaut2019fast]. For detailed theoretical discussions of the highly adaptive
lasso estimator, consider consulting, for example, @vdl2017generally,
@vdl2017finite, and @vdl2017uniform. For a computational demonstration of the
versatility of HAL regression, see @benkeser2016hal. Recent theoretical works
have demonstrated success in building efficient estimators of complex
parameters when particular variations of HAL regression are used to estimate
nuisance parameters [e.g., @vdl2019efficient; @ertefaie2020nonparametric].

---

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=hal9001) via

```{r cran-installation, eval = FALSE}
install.packages("hal9001")
```

To contribute, install the _development version_ of `hal9001` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

```{r gh-master-installation, eval = FALSE}
remotes::install_github("tlverse/hal9001")
```

---

## Issues

If you encounter any bugs or have any specific feature requests, please [file an
issue](https://github.com/tlverse/hal9001/issues).

---

## Example

Consider the following minimal example in using `hal9001` to generate
predictions via Highly Adaptive Lasso regression:

```{r example}
# load the package and set a seed
library(hal9001)
set.seed(385971)

# simulate data
n <- 100
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

# fit the HAL regression
hal_fit <- fit_hal(X = x, Y = y, yolo = TRUE)
hal_fit$times

# training sample prediction
preds <- predict(hal_fit, new_data = x)
mean(hal_mse <- (preds - y)^2)
```

---

## Contributions

Contributions are very welcome. Interested contributors should consult our
[contribution
guidelines](https://github.com/tlverse/hal9001/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

---

## Citation

After using the `hal9001` R package, please cite both of the following:

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

---

## License

&copy; 2017-2022 [Jeremy R. Coyle](https://github.com/tlverse) & [Nima S.
Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the GPL-3 license. See
file `LICENSE` for details.

---

## References

---
title: "Benchmarks for the `hal9001` Package"
author: "[Nima Hejazi](https://nimahejazi.org) and [Jeremy
  Coyle](https://github.com/jeremyrcoyle)"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: vignette-refs.bib
vignette: >
  %\VignetteIndexEntry{Benchmarks for the `hal9001` Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE, results='hide'}
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
set.seed(7194568)

library(ggplot2)
library(stringr)
library(tidyverse)

library(future)
library(microbenchmark)

library(SuperLearner)
library(hal9001)
```

## Introduction

This document consists of some simple benchmarks for various choices of the
`hal9001` implementation. The purpose of this document is two-fold:

1. Compare the computational performance of these methods
2. Illustrate the use of these different methods

```{r mse_fun}
# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y)^2)
}
```

---

## The Continuous Outcome Case

### Standard use-cases with simple data structures

```{r sim}
# generate simple test data
n = 100
p = 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

test_n = 10000
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(test_n, mean = 0,
                                                      sd = 0.2)
```

```{r hal_glmnet}
# glmnet implementation
glmnet_hal_fit <- fit_hal(X = x, Y = y, fit_type = "glmnet")
glmnet_hal_fit$times

# training sample prediction
preds_glmnet <- predict(glmnet_hal_fit, new_data = x)
glmnet_hal_mse <- mse(preds_glmnet, y)

# out-of-bag prediction
oob_preds_glmnet <- predict(glmnet_hal_fit, new_data = test_x)
oob_glmnet_hal_mse <- mse(oob_preds_glmnet, y = test_y)
```

```{r hal_lassi}
# lassi implementation
lassi_hal_fit <- fit_hal(X = x, Y = y, fit_type = "lassi")
lassi_hal_fit$times

# training sample prediction
preds_lassi <- predict(lassi_hal_fit, new_data = x)
lassi_hal_mse <- mse(preds_lassi, y)

# out-of-bag prediction
oob_preds_lassi <- predict(lassi_hal_fit, new_data = test_x)
oob_lassi_hal_mse <- mse(oob_preds_lassi, y = test_y)
```

```{r hal_microbenchmark}
mb_hal_lassi <- microbenchmark(unit = "s", times = 20,
  hal_fit_lassi <- fit_hal(Y = y, X = x, fit_type = "lassi", yolo = FALSE)
)
summary(mb_hal_lassi)

mb_hal_glmnet <- microbenchmark(unit = "s", times = 20,
  hal_fit_glmnet <- fit_hal(Y = y, X = x, fit_type = "glmnet", yolo = FALSE)
)
summary(mb_hal_glmnet)

# visualize
p_fit_times <- rbind(mb_hal_lassi, mb_hal_glmnet) %>%
  data.frame %>%
  transmute(
    type = ifelse(str_sub(as.character(expr), 9, 14) == "glmnet", "HAL-glmnet",
                  "HAL-lassi"),
    time = time / 1e9
  ) %>%
  group_by(type) %>%
  ggplot(., aes(x = type, y = time, colour = type)) + geom_boxplot() +
    geom_point() + stat_boxplot(geom = 'errorbar') +
    scale_color_manual(values = wes_palette("GrandBudapest")) +
    xlab("") + ylab("time (sec.)") + ggtitle("") +
    theme_bw() + theme(legend.position = "none")
p_fit_times
```

### Advanced use-cases, with the `drtmle` R package

```{r sim_drtmle}
makeData <- function(n = n) {
    L0 <- data.frame(x.0 = runif(n, -1, 1))
    A0 <- rbinom(n, 1, plogis(L0$x.0 ^ 2))
    L1 <- data.frame(x.1 = L0$x.0 ^ 2 * A0 + runif(n))
    A1 <- rbinom(n, 1, plogis(L0$x.0 * L1$x.1))
    L2 <- rnorm(n, L0$x.0 ^ 2 * A0 * A1 + L1$x.1)
    return(list(L0 = L0, L1 = L1, L2 = L2, A0 = A0, A1 = A1))
}
dat <- makeData(n = 200)

# add drtmle defaults
abar <- c(1, 1)
stratify <- TRUE
```

```{r hal_benchmark_drtmle_estimateG}
mb_hal_lassi_estimateG <- microbenchmark(unit = "s", times = 20,
  hal_lassi_estimateG <- fit_hal(Y = as.numeric(dat$A0 == abar[1]),
                                 X = dat$L0,
                                 fit_type = "lassi",
                                 yolo = FALSE)
)
summary(mb_hal_lassi_estimateG)

mb_hal_glmnet_estimateG <- microbenchmark(unit = "s", times = 20,
  hal_glmnet_estimateG <- fit_hal(Y = as.numeric(dat$A0 == abar[1]),
                                  X = dat$L0,
                                  fit_type = "glmnet",
                                  yolo = FALSE)
)
summary(mb_hal_glmnet_estimateG)

# visualize
p_estimateG_times <- rbind(mb_hal_lassi_estimateG,
                           mb_hal_glmnet_estimateG) %>%
  data.frame %>%
  transmute(
    type = ifelse(str_detect(as.character(expr), "glmnet"), "HAL-glmnet",
                  "HAL-lassi"),
    time = time / 1e9
  ) %>%
  group_by(type) %>%
  ggplot(., aes(x = type, y = time, colour = type)) + geom_boxplot() +
    geom_point() + stat_boxplot(geom = 'errorbar') +
    scale_color_manual(values = wes_palette("GrandBudapest")) +
    xlab("") + ylab("time (sec.)") + ggtitle("") +
    theme_bw() + theme(legend.position = "none")
p_estimateG_times
```

```{r hal_benchmark_drtmle_estimateQ}
Y_estQ <- if (stratify) {
  dat$L2[dat$A0 == abar[1] & dat$A1 == abar[2]]
} else {
  dat$L2
}
X_estQ = as.numeric(unlist(ifelse(stratify,
                                  cbind(dat$L0, dat$L1)[dat$A0 == abar[1] &
                                                        dat$A1 == abar[2],],
                                  cbind(dat$L0, dat$L1, dat$A0, dat$A1))))

mb_hal_lassi_estimateQ <- microbenchmark(unit = "s", times = 20,
  hal_lassi_estimateQ <- fit_hal(Y = Y_estQ,
                                 X = X_estQ,
                                 fit_type = "lassi",
                                 yolo = FALSE)
)
summary(mb_hal_lassi_estimateQ)

mb_hal_glmnet_estimateQ <- microbenchmark(unit = "s", times = 20,
  hal_glmnet_estimateQ <- fit_hal(Y = Y_estQ,
                                  X = X_estQ,
                                  fit_type = "glmnet",
                                  yolo = FALSE)
)
summary(mb_hal_glmnet_estimateQ)

# visualize
p_estimateQ_times <- rbind(mb_hal_lassi_estimateQ,
                           mb_hal_glmnet_estimateQ) %>%
  data.frame %>%
  transmute(
    type = ifelse(str_detect(as.character(expr), "glmnet"), "HAL-glmnet",
                  "HAL-lassi"),
    time = time / 1e9
  ) %>%
  group_by(type) %>%
  ggplot(., aes(x = type, y = time, colour = type)) + geom_boxplot() +
    geom_point() + stat_boxplot(geom = 'errorbar') +
    scale_color_manual(values = wes_palette("GrandBudapest")) +
    xlab("") + ylab("time (sec.)") + ggtitle("") +
    theme_bw() + theme(legend.position = "none")
p_estimateQ_times
```

### Using `hal9001` with `SuperLearner`

The `hal9001` R package includes a learner wrapper for use with the
`SuperLearner` R package. In order to use this learner, the user need only add
"SL.hal9001" to their `SL.library` input. Here, we are interested in the
difference in the performance of `SuperLearner` when the `hal9001` wrapper uses
the `lassi` backend versus the `glmnet` backend for the LASSO computation. To
facilitate this, we create two new SL wrappers from the one included:

```{r hal_sl_wrappers}
SL.hal9001.lassi <- function (..., fit_type = "lassi") {
  SL.hal9001(..., fit_type = fit_type)
}

SL.hal9001.glmnet <- function (..., fit_type = "glmnet") {
  SL.hal9001(..., fit_type = fit_type)
}
```


```{r hal_benchmark_drtmle_estimateG_SL}
mb_hal_lassi_sl_estimateG <- microbenchmark(unit = "s", times = 20,
  hal_lassi_sl_estimateG <- SuperLearner(Y = as.numeric(dat$A0 == abar[1]),
                                         X = dat$L0,
                                         SL.library = c("SL.speedglm",
                                                        "SL.hal9001.lassi")
                                        )
)
summary(mb_hal_lassi_sl_estimateG)

mb_hal_glmnet_sl_estimateG <- microbenchmark(unit = "s", times = 20,
  hal_glmnet_sl_estimateG <- SuperLearner(Y = as.numeric(dat$A0 == abar[1]),
                                          X = dat$L0,
                                          SL.library = c("SL.speedglm",
                                                         "SL.hal9001.glmnet")
                                         )
)
summary(mb_hal_glmnet_sl_estimateG)

# visualize
p_sl_estimateG_times <- rbind(mb_hal_lassi_sl_estimateG,
                              mb_hal_glmnet_sl_estimateG) %>%
  data.frame %>%
  transmute(
    type = ifelse(str_detect(as.character(expr), "glmnet"), "HAL-glmnet",
                  "HAL-lassi"),
    time = time / 1e9
  ) %>%
  group_by(type) %>%
  ggplot(., aes(x = type, y = time, colour = type)) + geom_boxplot() +
    geom_point() + stat_boxplot(geom = 'errorbar') +
    scale_color_manual(values = wes_palette("GrandBudapest")) +
    xlab("") + ylab("time (sec.)") +
    ggtitle("SuperLearner with HAL glmnet v. lassi") +
    theme_bw() + theme(legend.position = "none")
p_sl_estimateG_times
```


## The Binary Outcome Case

```{r sim_binary}
# generate simple test data
n = 100
p = 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y_p <- plogis(sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2))
y <- rbinom(n, size = 1, prob = y_p)

glmnet_hal_bin <- fit_hal(Y = y, X = x, fit_type = "glmnet",
                          family = "binomial")
glmnet_hal_bin$times

mb_hal_bin_glmnet <- microbenchmark(unit = "s", times = 20,
  hal_bin_glmnet <- fit_hal(Y = y, X = x,
                            fit_type = "glmnet",
                            family = "binomial",
                            yolo = FALSE)
)
summary(mb_hal_bin_glmnet)
```

---
title: "Evaluating HAL: R vs. Rcpp"
author: "[Nima Hejazi](http://nimahejazi.org) & Jeremy Coyle"
date: "`r Sys.Date()`"
output: html_document
---

```{r, echo=FALSE}
set.seed(479253)
library(dplyr)
library(BH)
library(Rcpp)
library(RcppArmadillo)
library(Rcereal)
library(microbenchmark)
```

## Introduction

...

---

## Analysis

### `hal`: A pure R implementation

...

### `hal9000`: An Rcpp implementation

```{Rcpp, ref.label=knitr::all_rcpp_labels(), cache=TRUE, include=FALSE}
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mangolassi_types.h"
using namespace Rcpp;
```

```{Rcpp dedupe, eval=FALSE}
// returns an index vector indicating index of first copy of column for each X
// [[Rcpp::export]]
IntegerVector index_first_copy(const MSpMat& X){
  int p = X.cols();

  ColMap col_map;
  IntegerVector copy_index(p);

  for (int j = 0; j < p; j++){
    MSpMatCol current_col(X, j);

    //https://stackoverflow.com/questions/97050/stdmap-insert-or-stdmap-find
    ColMap::iterator match = col_map.lower_bound(current_col);
    if (match != col_map.end() && !(col_map.key_comp()(current_col, match -> first)))
    {
      // column already exists
      copy_index[j] = match -> second + 1;  // just using 1-indexing
    } else {
      // column not yet in map
      col_map.insert(match, ColMap::value_type(current_col, j));
      copy_index[j] = j + 1;  //just using 1-indexing
    }
  }
  return(copy_index);
}

// returns true iff col_1 is strictly less than col_2 in the ordering scheme
// [[Rcpp::export]]
bool column_compare(const MSpMat& X, int col_1, int col_2) {
  ColMap cmap;
  MSpMatCol X_1(X, col_1);
  MSpMatCol X_2(X, col_2);
  return(cmap.key_comp()(X_1, X_2));
}

// ORs the columns of X listed in cols and places the result in column col[1]
// [[Rcpp::export]]
void or_duplicate_columns(MSpMat& X, const IntegerVector& cols) {
  int first = cols[0] - 1;  //cols is 1-indexed
  int p_cols = cols.length();
  int n = X.rows();
  for (int i = 0; i < n; i++){
    if(X.coeffRef(i, first) == 1) {
      continue;  // this is already 1
    }

    //search remaining columns for 1, inserting into first if found
    for (int j = 1; j < p_cols; j++) {
      int j_col = cols[j] - 1;  //cols is 1-indexed
      if (X.coeffRef(i, j_col) == 1) {
        X.coeffRef(i, j_col) = 1;
        break;
      }
    }
  }
}
```

```{Rcpp lassi-1, eval=FALSE}
// [[Rcpp::export]]
NumericVector lassi_predict(MSpMat X, NumericVector beta) {
  int n = X.rows();
  NumericVector pred(n, beta[0]);  // initialize with intercept
  int k = 0;
  double current_beta;

  for (k = 0; k < X.outerSize(); ++k) {
    current_beta = beta[k + 1];

    for (MInIterMat it_(X, k); it_; ++it_) {
      pred[it_.row()] += current_beta;
    }
  }
  return(pred);
}

double soft_threshold(double beta, double lambda) {
  if (beta > lambda) {
    beta -= lambda;
  } else if (beta < -1*lambda) {
    beta += lambda;
  } else {
    beta = 0;
  }
  return(beta);
}
```

_Code for coordinate descent:_
based on http://www.stanford.edu/~hastie/Papers/glmnet.pdf
TODO: implement scaling and centering:
Coordinate descent is ideally set up to exploit such sparsity, in an obvious
way. The O(N) inner-product operations in either the naive or covariance
updates can exploit the sparsity, by summing over only the non-zero entries.
Note that in this case scaling of the variables will not alter the sparsity,
but centering will. So scaling is performed up front, but the centering is
incorporated in the algorithm in an efficient and obvious manner.

```{Rcpp lassi-2, eval=FALSE}
// [[Rcpp::export]]
void update_coord(MSpMat X, NumericVector resids, NumericVector beta, double lambda, int j) {

  int n = resids.length();
  double beta_j = beta[j + 1];  //+1 for intercept
  double new_beta = 0;
  double resid_sum = 0;

  for (MInIterMat i_(X, j); i_; ++i_) {
    resid_sum += resids[i_.index()];
  }

  new_beta = resid_sum / n + beta_j;
  new_beta = soft_threshold(new_beta, lambda);

  // if we changed this beta, we must update the residuals
  if (new_beta != beta_j) {
    // Rcout << "Changed beta " << j << std::endl;
    double beta_diff = new_beta-beta_j;
    for (MInIterMat i_(X, j); i_; ++i_) {
      resids[i_.index()] -= beta_diff;
    }
    beta[j + 1] = new_beta;
  }
}

void update_coords(MSpMat X, NumericVector resids, NumericVector beta, double lambda){
  //update coordinates one-by-one
  int k;
  for (k = 0; k < X.outerSize(); ++k) {
    update_coord(X, resids, beta, lambda, k);
  }

  //update intercept to center predictions
  double mean_resid = mean(resids);
  resids = resids-mean_resid;
  beta[0] += mean_resid;
}
```

```{Rcpp lassi-3, eval=FALSE}
// [[Rcpp::export]]
NumericVector lassi_fit_cd(MSpMat X, NumericVector y, double lambda, int nsteps){
  int p = X.cols();
  NumericVector beta(p + 1);
  NumericVector resids = y-lassi_predict(X, beta);

  int step_num = 0;

  double mse = mean(resids*resids);
  double last_mse = mse;
  double ratio = 0;

  for (step_num = 0; step_num < nsteps; step_num++) {
    last_mse = mse;

    update_coords(X, resids, beta, lambda);
    mse = mean(resids*resids);
    ratio = (last_mse - mse) / last_mse;

    Rcout << "Step " << step_num << ", mse " << mse << ", ratio " << ratio << std::endl;
    if (ratio < 0.001) {
      break;
    }
  }
  return(beta);
}
```

The following are functions to enumerate basis functions:

```{Rcpp hal-basis-1, eval=FALSE}
// populates a map with unique basis functions based on data in xsub
// values are thresholds, keys are column indicies
BasisMap enumerate_basis(const NumericMatrix& X_sub, const NumericVector& cols){
  BasisMap bmap;

  //find unique basis functions
  int n = X_sub.rows();
  for (int i = 0; i < n; i++) {
    NumericVector cutoffs = X_sub.row(i);
    bmap.insert(std::pair<NumericVector, NumericVector>(cutoffs, cols));
  }
  return(bmap);
}
```

returns a sorted list of unique basis functions based on columns in cols (so basis order=cols.length())
each basis function is a list(cols,cutoffs)
X_sub is a subset of the columns of X (the original design matrix)
cols is an index of the columns that were subsetted

```{Rcpp hal-basis-2, eval=FALSE}
// [[Rcpp::export]]
List make_basis_list(const NumericMatrix& X_sub, const NumericVector& cols) {

  BasisMap bmap = enumerate_basis(X_sub, cols);
  List basis_list(bmap.size());
  int index = 0;
  for (BasisMap::iterator it = bmap.begin(); it != bmap.end(); ++it) {
    List basis = List::create(
      Rcpp::Named("cols") = it -> second,
      Rcpp::Named("cutoffs") = it -> first
    );

    basis_list[index++] = basis;
  }
  return(basis_list);
}
```

Functions to make a design matrix based on a list of basis functions

```{Rcpp hal-design-1, eval=FALSE}
// returns the indicator value for the basis described by cols,cutoffs for X[row_num,]
// X is the original design matrix
// row_num is a row index to evaluate
// cols are the column incides of the basis function
// cutoffs are thresholds
// [[Rcpp::export]]
bool meets_basis(const NumericMatrix& X, const int row_num, const IntegerVector& cols, const NumericVector& cutoffs) {
  int p = cols.length();
  for (int i = 0; i < p; i++) {
    double obs = X(row_num, cols[i] - 1);  //we're using 1-indexing for basis columns
    if (!(obs >= cutoffs[i])) {
      return(false);
    }
  }
  return(true);
}
```

```{Rcpp hal-design-2, eval=FALSE}
// populates a column (indexed by basis_col) of x_basis with basis indicators
// basis is the basis function
// X is the original design matrix
// x_basis is the hal design matrix
// basis_col indicates which column to populate
// [[Rcpp::export]]
void evaluate_basis(const List& basis, const NumericMatrix& X, SpMat& x_basis, int basis_col){
  int n=X.rows();
  //split basis into x[1] x[-1]
  //find sub-basises
  //intersect

  IntegerVector cols=as<IntegerVector>(basis["cols"]);
  NumericVector cutoffs=as<NumericVector>(basis["cutoffs"]);
  for (int row_num = 0; row_num < n; row_num++) {

    if (meets_basis(X, row_num, cols, cutoffs)) {
      //we can add a positive indicator for this row, basis
      x_basis.insert(row_num, basis_col) = 1;
    }
  }
}
```

```{Rcpp hal-design-3, eval=FALSE}
// makes a hal design matrix based on original design matrix X and
// a list of basis functions in blist
// [[Rcpp::export]]
SpMat make_design_matrix(NumericMatrix X, List blist) {
  //now generate an indicator vector for each
  int n = X.rows();
  int basis_p = blist.size();

  SpMat x_basis(n,basis_p);
  x_basis.reserve(0.5*n*basis_p);

  List basis;
  NumericVector cutoffs, current_row;
  IntegerVector last_cols, cols;
  NumericMatrix X_sub;

  //for each basis function
  for (int basis_col = 0; basis_col < basis_p; basis_col++) {
    last_cols = cols;

    basis = (List) blist[basis_col];
    evaluate_basis(basis, X, x_basis, basis_col);
  }
  x_basis.makeCompressed();
  return(x_basis);
}
```

---

## Discussion

...

### Future Work

...

---

## References

---
title: "hal90001 Benchmarks"
author: "Jeremy Coyle"
date: "10/5/2017"
output: html_document
---

```{r setup, include=FALSE, results='hide'}
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
library(sl3)
library(delayed)
library(SuperLearner)
library(future)
library(ggplot2)
library(data.table)
library(stringr)
library(scales)
```

## Introduction

This document consists of some simple benchmarks for various choices of SuperLearner implementation, wrapper functions, and parallelization schemes. The purpose of this document is two-fold: 

1. Compare the computational performance of these methods
2. Illustrate the use of these different methods

## Test Setup


### Test System

```{r systemInfo, echo=FALSE, results="asis"}
uname <- system("uname -a", intern = TRUE)
os <- sub(" .*", "", uname)
if(os=="Darwin"){
  cpu_model <- system("sysctl -n machdep.cpu.brand_string", intern = TRUE)
  cpus_physical <- as.numeric(system("sysctl -n hw.physicalcpu", intern = TRUE))
  cpus_logical <- as.numeric(system("sysctl -n hw.logicalcpu", intern = TRUE))
  cpu_clock <- system("sysctl -n hw.cpufrequency_max", intern = TRUE)
  memory <- system("sysctl -n hw.memsize", intern = TRUE)
} else if(os=="Linux"){
  cpu_model <- system("lscpu | grep 'Model name'", intern = TRUE)
  cpu_model <- gsub("Model name:[[:blank:]]*","", cpu_model)
  cpus_logical <- system("lscpu | grep '^CPU(s)'", intern = TRUE)
  cpus_logical <- as.numeric(gsub("^.*:[[:blank:]]*","", cpus_logical))
  tpc <- system("lscpu | grep '^Thread(s) per core'", intern = TRUE)
  tpc <- as.numeric(gsub("^.*:[[:blank:]]*","", tpc))
  cpus_physical <- cpus_logical/tpc
  cpu_clock <- as.numeric(gsub("GHz","",gsub("^.*@","",cpu_model)))*10^9
  memory <- system("cat /proc/meminfo | grep '^MemTotal'", intern = TRUE)
  memory <- as.numeric(gsub("kB","",gsub("^.*:","",memory)))*2^10
} else {
  stop("unsupported OS")
}
```

* CPU model: `r cpu_model`
* Physical cores: `r as.numeric(cpus_physical)`
* Logical cores: `r as.numeric(cpus_logical)`
* Clock speed: `r as.numeric(cpu_clock)/10^9`GHz
* Memory: `r round(as.numeric(memory)/2^30, 1)`GB

### Test Data

### Tests

```{r lassi}
microbenchmark({
   glmnet::glmnet(x = x_basis, y = y, intercept = TRUE, nlambda = 100,
   lambda.min.ratio = 0.01, family = "gaussian", alpha = 1, standardize = TRUE)
 }, times = 10)
microbenchmark({
   lassi(x_basis, y, nlambda=100, lambda_min_ratio = 0.01, center = FALSE)
 }, times = 10)

microbenchmark({
   lassi(x_basis, y, nlambda=100, lambda_min_ratio = 0.01, center = TRUE)
 }, times = 10)

microbenchmark({
   glmnet::cv.glmnet(x = x_basis, y = y, intercept = TRUE, nlambda = 100,
   lambda.min.ratio = 0.01, family = "gaussian", alpha = 1, standardize = TRUE)
 }, times = 1)

microbenchmark({
   cv_lasso(x_basis, y, center = FALSE)
 }, times = 1)

microbenchmark({
   cv_lasso(x_basis, y, center = TRUE)
 }, times = 1)

set.seed(1234)
cv_l_full <- cv_lasso(x_basis, y, center = FALSE)
set.seed(1234)
cv_l_es <- cv_lasso_early_stopping(x_basis, y)
plot(cv_l_es)
plot(cv_l_full$lambdas_cvmse)
microbenchmark({
   cv_lasso_early_stopping(x_basis, y)
 }, times = 1)

```

## Session Information

```{r sessionInfo, echo=FALSE, results="asis"}
sessionInfo()
```---
title: "Fitting the Highly Adaptive Lasso with `hal9001`"
author: "[Nima Hejazi](https://nimahejazi.org), [Jeremy
  Coyle](https://github.com/jeremyrcoyle), Rachael Phillips, Lars van der Laan"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: ../inst/REFERENCES.bib
vignette: >
  %\VignetteIndexEntry{Fitting the Highly Adaptive Lasso with hal9001}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

The _highly adaptive Lasso_ (HAL) is a flexible machine learning algorithm that
nonparametrically estimates a function based on available data by embedding a
set of input observations and covariates in an extremely high-dimensional space
(i.e., generating basis functions from the available data). For an input data
matrix of $n$ observations and $d$ covariates, the maximum number of zero-order
basis functions generated is approximately $n \cdot 2^{d - 1}$. To select a set
of basis functions from among the (possibly reduced/screener) set that's
generated, the lasso is employed. The `hal9001` R package [@hejazi2020hal9001;
@coyle-gh-hal9001] provides an efficient implementation of this routine, relying
on the `glmnet` R package [@friedman2010glmnet] for compatibility with the
canonical Lasso implementation and using lasso regression with an input matrix
composed of basis functions. Consult @benkeser2016hal, @vdl2015generally,
@vdl2017finite for detailed theoretical descriptions of HAL and its various
optimality properties.

---

## Preliminaries

```{r sim-data}
library(data.table)
library(ggplot2)
# simulation constants
set.seed(467392)
n_obs <- 500
n_covars <- 3

# make some training data
x <- replicate(n_covars, rnorm(n_obs))
y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n_obs, mean = 0, sd = 0.2)

# make some testing data
test_x <- replicate(n_covars, rnorm(n_obs))
test_y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n_obs, mean = 0, sd = 0.2)
```

Let's look at simulated data:

```{r sim-view}
head(x)
head(y)
```

## Using the Highly Adaptive Lasso

```{r}
library(hal9001)
```

### Fitting the model

HAL uses the popular `glmnet` R package for the lasso step:

```{r fit-hal-glmnet}
hal_fit <- fit_hal(X = x, Y = y)
hal_fit$times
```

### Summarizing the model

While the raw output object may be examined, it has (usually large) slots that
make quick examination challenging. The `summary` method provides an
interpretable table of basis functions with non-zero coefficients. All terms
(i.e., including the terms with zero coefficient) can be included by setting
`only_nonzero_coefs` to `FALSE` when calling `summary` on a `hal9001` model
object.

```{r results-hal-glmnet}
print(summary(hal_fit))
```

Note the length and width of these tables! The R environment might not be the
optimal location to view the summary. Tip: Tables can be exported from R to
LaTeX with the `xtable` R package. Here's an example:
`print(xtable(summary(fit)$table, type = "latex"), file = "haltbl_meow.tex")`.

### Obtaining model predictions

```{r eval-mse}
# training sample prediction for HAL vs HAL9000
mse <- function(preds, y) {
    mean((preds - y)^2)
}

preds_hal <- predict(object = hal_fit, new_data = x)
mse_hal <- mse(preds = preds_hal, y = y)
mse_hal
```

```{r eval-oob}
oob_hal <- predict(object = hal_fit, new_data = test_x)
oob_hal_mse <- mse(preds = oob_hal, y = test_y)
oob_hal_mse
```

### Reducing basis functions

As described in @benkeser2016hal, the HAL algorithm operates by first
constructing a set of basis functions and subsequently fitting a Lasso model
with this set of basis functions as the design matrix. Several approaches are
considered for reducing this set of basis functions:
1. Removing duplicated basis functions (done by default in the `fit_hal`
   function),
2. Removing basis functions that correspond to only a small set of observations;
   a good rule of thumb is to scale with $\frac{1}{\sqrt{n}}$, and that is the
   default.

The second of these two options may be modified by specifying the `reduce_basis`
argument to the `fit_hal` function:

```{r fit-hal-reduced}
hal_fit_reduced <- fit_hal(X = x, Y = y, reduce_basis = 0.1)
hal_fit_reduced$times
```

In the above, all basis functions with fewer than 10% of observations meeting
the criterion imposed are automatically removed prior to the Lasso step of
fitting the HAL regression. The results appear below

```{r results-hal-reduced}
summary(hal_fit_reduced)$table
```

Other approaches exist for reducing the set of basis functions *before* they are
actually created, which is essential for most real-world applications with HAL.
Currently, we provide this "pre-screening" via `num_knots` argument in
`hal_fit`. The `num_knots` argument is akin to binning: it increases the
coarseness of the approximation. `num_knots` allows one to specify the number of
knot points used to generate the basis functions for each/all interaction
degree(s). This reduces the total number of basis functions generated, and thus
the size of the optimization problem, and it can dramatically decrease runtime.
One can pass in a vector of length `max_degree` to `num_knots`, specifying the
number of knot points to use by interaction degree for each basis function.
Thus, one can specify if interactions of higher degrees (e.g., two- or three-
way interactions) should be more coarse.  Increasing the coarseness of more
complex basis functions helps prevent a combinatorial explosion of basis
functions, which can easily occur when basis functions are generated for all
possible knot points. We will show an example with `num_knots` in the section
that follows.

### Specifying smoothness of the HAL model

One might wish to enforce smoothness on the functional form of the HAL fit.
This can be done using the `smoothness_orders` argument. Setting
`smoothness_orders = 0` gives a piece-wise constant fit (via zero-order basis
functions), allowing for discontinuous jumps in the function. This is useful if
one does not want to assume any smoothness or continuity of the "true" function.
Setting `smoothness_orders = 1` gives a piece-wise linear fit (via first-order
basis functions), which is continuous and mostly differentiable. In general,
`smoothness_orders = k` corresponds to a piece-wise polynomial fit of degree
$k$. Mathematically, `smoothness_orders = k` corresponds with finding the best
fit under the constraint that the total variation of the function's
$k^{\text{th}}$ derivative is bounded by some constant, which is selected with
cross-validation.

Let's see this in action.

```{r}
set.seed(98109)
num_knots <- 100 # Try changing this value to see what happens.
n_covars <- 1
n_obs <- 250
x <- replicate(n_covars, runif(n_obs, min = -4, max = 4))
y <- sin(x[, 1]) + rnorm(n_obs, mean = 0, sd = 0.2)
ytrue <- sin(x[, 1])

hal_fit_0 <- fit_hal(
  X = x, Y = y, smoothness_orders = 0, num_knots = num_knots
)

hal_fit_smooth_1 <- fit_hal(
  X = x, Y = y, smoothness_orders = 1, num_knots = num_knots
)

hal_fit_smooth_2_all <- fit_hal(
  X = x, Y = y, smoothness_orders = 2, num_knots = num_knots, 
  fit_control = list(cv_select = FALSE)
)

hal_fit_smooth_2 <- fit_hal(
  X = x, Y = y, smoothness_orders = 2, num_knots = num_knots
)

pred_0 <- predict(hal_fit_0, new_data = x)
pred_smooth_1 <- predict(hal_fit_smooth_1, new_data = x)
pred_smooth_2 <- predict(hal_fit_smooth_2, new_data = x)

pred_smooth_2_all <- predict(hal_fit_smooth_2_all, new_data = x)
dt <- data.table(x = as.vector(x))
dt <- cbind(dt, pred_smooth_2_all)
long <- melt(dt, id = "x")
ggplot(long, aes(x = x, y = value, group = variable)) + geom_line()
```

Comparing the mean squared error (MSE) between the predictions and the true
(denoised) outcome, the first- and second- order smoothed HAL is able to recover
from the coarseness of the basis functions caused by the small `num_knots`
argument. Also, the HAL with second-order smoothness is able to fit the true
function very well (as expected, since sin(x) is a very smooth function). The
main benefit of imposing higher-order smoothness is that fewer knot points are
required for a near-optimal fit. Therefore, one can safely pass a smaller value
to `num_knots` for a big decrease in runtime without sacrificing performance.

```{r}
mean((pred_0 - ytrue)^2)
mean((pred_smooth_1- ytrue)^2)
mean((pred_smooth_2 - ytrue)^2)

dt <- data.table(x = as.vector(x),
                 ytrue = ytrue,
                 y = y,
                 pred0 = pred_0,
                 pred1 = pred_smooth_1,
                 pred2 = pred_smooth_2)
long <- melt(dt, id = "x")
ggplot(long, aes(x = x, y = value, color = variable)) + geom_line()
plot(x, pred_0, main = "Zero order smoothness fit")
plot(x, pred_smooth_1, main = "First order smoothness fit")
plot(x, pred_smooth_2, main = "Second order smoothness fit")
```

In general, if the basis functions are not coarse, then the performance for
different smoothness orders is similar. Notice how the runtime is a fair bit
slower when more knot points are considered. In general, we recommend either
zero- or first- order smoothness. Second-order smoothness tends to be less
robust and suffers from extrapolation on new data. One can also use
cross-validation to data-adaptively choose the optimal smoothness (invoked in
`fit_hal` by setting `adaptive_smoothing = TRUE`). Comparing the following
simulation and the previous one, the HAL with second-order smoothness performed
better when there were fewer knot points.

```{r}
set.seed(98109)
n_covars <- 1
n_obs <- 250
x <- replicate(n_covars, runif(n_obs, min = -4, max = 4))
y <- sin(x[, 1]) + rnorm(n_obs, mean = 0, sd = 0.2)
ytrue <- sin(x[, 1])

hal_fit_0 <- fit_hal(
  X = x, Y = y, smoothness_orders = 0, num_knots = 100
)
hal_fit_smooth_1 <- fit_hal(
  X = x, Y = y, smoothness_orders = 1, num_knots = 100
)

hal_fit_smooth_2 <- fit_hal(
  X = x, Y = y, smoothness_orders = 2, num_knots = 100
)

pred_0 <- predict(hal_fit_0, new_data = x)
pred_smooth_1 <- predict(hal_fit_smooth_1, new_data = x)
pred_smooth_2 <- predict(hal_fit_smooth_2, new_data = x)
```

```{r}
mean((pred_0 - ytrue)^2)
mean((pred_smooth_1- ytrue)^2)
mean((pred_smooth_2 - ytrue)^2)
plot(x, pred_0, main = "Zero order smoothness fit")
plot(x, pred_smooth_1, main = "First order smoothness fit")
plot(x, pred_smooth_2, main = "Second order smoothness fit")
```

### Formula interface

One might wish to specify the functional form of the HAL fit further. This can
be done using the formula interface. Specifically, the formula interface allows
one to specify monotonicity constraints on components of the HAL fit. It also
allows one to specify exactly which basis functions (e.g., interactions) one
wishes to model. The `formula_hal` function generates a `formula` object from a
user-supplied character string, and this `formula` object contains the necessary
specification information for `fit_hal` and `glmnet`.  The `formula_hal`
function is intended for use within `fit_hal`, and the user-supplied character
string is inputted into `fit_hal`. Here, we call `formula_hal` directly for
illustrative purposes.

```{r}
set.seed(98109)
num_knots <- 100

n_obs <- 500
x1 <-  runif(n_obs, min = -4, max = 4)
x2 <- runif(n_obs, min = -4, max = 4)
A <- runif(n_obs, min = -4, max = 4)
X <- data.frame(x1 = x1, x2 = x2, A = A)
Y <- rowMeans(sin(X)) + rnorm(n_obs, mean = 0, sd = 0.2)
```

We can specify an additive model in a number of ways. 

The formula below includes the outcome, but `formula_hal` doesn't fit a HAL
model, and doesn't need the outcome (actually everything before "$\tilde$" is
ignored in `formula_hal`). This is why `formula_hal` takes the input `X` matrix
of covariates, and not `X` and `Y`. In what follows, we include formulas with
and without "y" in the character string.

```{r}
# The `h` function is used to specify the basis functions for a given term
# h(x1) generates one-way basis functions for the variable x1
# This is an additive model:
formula <- ~h(x1) + h(x2) + h(A)
#We can actually evaluate the h function as well. We need to specify some tuning parameters in the current environment:
smoothness_orders <- 0
num_knots <- 10
# It will look in the parent environment for `X` and the above tuning parameters
form_term <- h(x1) + h(x2) + h(A)
form_term$basis_list[[1]]
# We don't need the variables in the parent environment if we specify them directly:
rm(smoothness_orders)
rm(num_knots)
# `h` excepts the arguments `s` and `k`. `s` stands for smoothness and is equivalent to smoothness_orders in use. `k` specifies the number of knots. ` 
form_term_new <- h(x1, s = 0, k = 10) + h(x2, s = 0, k = 10) + h(A, s = 0, k = 10)
# They are the same!
length( form_term_new$basis_list) == length(form_term$basis_list)

#To evaluate a unevaluated formula object like:
formula <- ~h(x1) + h(x2) + h(A)
# we can use the formula_hal function:
formula <- formula_hal(
  ~ h(x1) + h(x2) + h(A), X = X, smoothness_orders = 1, num_knots = 10
)
# Note that the arguments smoothness_orders and/or num_knots will not be used if `s` and/or `k` are specified in `h`.
formula <- formula_hal(
  Y ~ h(x1, k=1) + h(x2,  k=1) + h(A, k=1), X = X, smoothness_orders = 1, num_knots = 10
)
 
 
 
 
```

The `.` argument. We can generate an additive model for all or a subset of variables using the `.` variable and `.` argument of `h`. By default, `.` in `h(.)` is treated as a wildcard and basis functions are generated by replacing the `.` with all variables in `X`. 
```{r}
smoothness_orders <- 1
num_knots <- 5
# A additive model
colnames(X)
# Shortcut:
formula1 <- h(.)
# Longcut:
formula2 <- h(x1) + h(x2) + h(A)
# Same number of basis functions
length(formula1$basis_list ) == length(formula2$basis_list)

# Maybe we only want an additive model for x1 and x2
# Use the `.` argument
formula1 <- h(., . = c("x1", "x2"))
formula2 <- h(x1) + h(x2) 
length(formula1$basis_list ) == length(formula2$basis_list)

```

We can specify interactions as follows.
```{r}
#  Two way interactions
formula1 <-  h(x1) + h(x2)  + h(A) + h(x1, x2)
formula2 <-  h(.)  + h(x1, x2)
length(formula1$basis_list ) == length(formula2$basis_list)
#
formula1 <-  h(.) + h(x1, x2)    + h(x1,A) + h(x2,A)
formula2 <-  h(.)  + h(., .)
length(formula1$basis_list ) == length(formula2$basis_list)

#  Three way interactions
formula1 <-  h(.) + h(.,.)    + h(x1,A,x2)  
formula2 <-  h(.)  + h(., .)+ h(.,.,.)  
length(formula1$basis_list ) == length(formula2$basis_list)

 
```

Sometimes, one might want to build an additive model, but include all two-way
interactions with one variable (e.g., treatment "A"). This can be done in a
variety of ways. The `.` argument allows you to specify a subset of variables.

```{r}
# Write it all out
formula <-  h(x1) + h(x2) + h(A) + h(A, x1) + h(A,x2)
 

# Use the "h(.)" which stands for add all additive terms and then manually add
# interactions
formula <- y ~ h(.) + h(A,x1) + h(A,x2)

 

# Use the "wildcard" feature for when "." is included in the "h()" term. This
# useful when you have many variables and do not want to write out every term.
formula <-  h(.) + h(A,.) 


formula1 <-   h(A,x1) 
formula2 <-   h(A,., . = c("x1")) 
 length(formula1$basis_list) == length(formula2$basis_list)
```

 


A key feature of the HAL formula is **monotonicity constraints**. Specifying
these constraints is achieved by specifying the `monotone` argument of `h`. Note if smoothness_orders = 0 then this is a monotonicity constrain on the function, but if if smoothness_orders = 1 then this is a monotonicity constraint on the function's derivative (e.g. a convexity constraint). We can also specify that certain terms are not penalized in the LASSO/glmnet using the `pf` argument of `h` (stands for penalty factor).

```{r}
# An additive monotone increasing model
formula <- formula_hal(
  y ~ h(., monotone = "i"), X, smoothness_orders = 0, num_knots = 100
)

# An additive unpenalized monotone increasing model (NPMLE isotonic regressio)
# Set the penalty factor argument `pf` to remove L1 penalization
formula <- formula_hal(
  y ~ h(., monotone = "i", pf = 0), X, smoothness_orders = 0, num_knots = 100
)

# An additive unpenalized convex model (NPMLE convex regressio)
# Set the penalty factor argument `pf` to remove L1 penalization
# Note the second term is equivalent to adding unpenalized and unconstrained main-terms (e.g. main-term glm)
formula <- formula_hal(
  ~ h(., monotone = "i", pf = 0, k=200, s=1) + h(., monotone = "none", pf = 0, k=1, s=1), X)
 
# A bi-additive monotone decreasing model
formula <- formula_hal(
  ~ h(., monotone = "d") + h(.,., monotone = "d"), X, smoothness_orders = 1, num_knots = 100
)
 

 
 
```


The penalization feature can be used to reproduce glm

```{r}
# Additive glm
# One knot (at the origin) and first order smoothness
formula <- h(., s = 1, k = 1, pf = 0)
# Running HAL with this formula will be equivalent to running glm with the formula Y ~ .

# intraction glm
formula <- h(., ., s = 1, k = 1, pf = 0) + h(., s = 1, k = 1, pf = 0)
# Running HAL with this formula will be equivalent to running glm with the formula Y ~ .^2

```


Now, that we've illustrated the options with `formula_hal`, let's show how to
fit a HAL model with the specified formula.
```{r}
# get formula object
fit <- fit_hal(
  X = X, Y = Y, formula = ~ h(.), smoothness_orders = 1, num_knots = 100
)
print(summary(fit), 10) # prints top 10 rows, i.e., highest absolute coefs
plot(predict(fit, new_data = X), Y)
```

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/squash_hal.R
\name{squash_hal_fit}
\alias{squash_hal_fit}
\title{Squash HAL objects}
\usage{
squash_hal_fit(object)
}
\arguments{
\item{object}{An object of class \code{hal9001}, containing the results of
fitting the Highly Adaptive LASSO, as produced by a call to \code{fit_hal}.}
}
\value{
Object of class \code{hal9001}, similar to the input object but
reduced such that coefficients belonging to bases with coefficients equal
to zero removed.
}
\description{
Reduce footprint by dropping basis functions with coefficients of zero
}
\examples{
\donttest{
# generate simple test data
n <- 100
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

# fit HAL model and squash resulting object to reduce footprint
hal_fit <- fit_hal(X = x, Y = y, yolo = FALSE)
squashed <- squash_hal_fit(hal_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{print.summary.hal9001}
\alias{print.summary.hal9001}
\title{Print Method for Summary Class of HAL fits}
\usage{
\method{print}{summary.hal9001}(x, length = NULL, ...)
}
\arguments{
\item{x}{An object of class \code{summary.hal9001}.}

\item{length}{The number of ranked coefficients to be summarized.}

\item{...}{Other arguments (ignored).}
}
\description{
Print Method for Summary Class of HAL fits
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula_hal9001.R
\name{h}
\alias{h}
\title{HAL Formula term: Generate a single term of the HAL basis}
\usage{
h(
  ...,
  k = NULL,
  s = NULL,
  pf = 1,
  monotone = c("none", "i", "d"),
  . = NULL,
  dot_args_as_string = FALSE,
  X = NULL
)
}
\arguments{
\item{...}{Variables for which to generate multivariate interaction basis
function where the variables can be found in a matrix \code{X} in a parent
environment/frame. Note, just like standard \code{formula} objects, the
variables should not be characters (e.g. do h(W1,W2) not h("W1", "W2"))
h(W1,W2,W3) will generate three-way HAL basis functions between W1, W2, and
W3. It will \code{not} generate the lower dimensional basis functions.}

\item{k}{The number of knots for each univariate basis function used to
generate the tensor product basis functions. If a single value then this
value is used for the univariate basis functions for each variable.
Otherwise, this should be a variable named list that specifies for each
variable how many knots points should be used.
\code{h(W1,W2,W3, k = list(W1 = 3, W2 = 2, W3=1))} is equivalent to first
binning the variables \code{W1}, \code{W2} and \code{W3} into \code{3}, \code{2} and \code{1} unique
values and then calling \code{h(W1,W2,W3)}. This coarsening of the data ensures
that fewer basis functions are generated, which can lead to substantial
computational speed-ups. If not provided and the variable \code{num_knots}
is in the parent environment, then \code{s} will be set to
\code{num_knots}`.}

\item{s}{The \code{smoothness_orders} for the basis functions. The possible
values are \code{0} for piece-wise constant zero-order splines or \code{1} for
piece-wise linear first-order splines. If not provided and the variable
\code{smoothness_orders} is in the parent environment, then \code{s} will
be set to \code{smoothness_orders}.}

\item{pf}{A \code{penalty.factor} value the generated basis functions that is
used by \code{glmnet} in the LASSO penalization procedure. \code{pf = 1}
(default) is the standard penalization factor used by \code{glmnet} and
\code{pf = 0} means the generated basis functions are unpenalized.}

\item{monotone}{Whether the basis functions should enforce monotonicity of
the interaction term. If \verb{\code{s} = 0}, this is monotonicity of the
function, and, if \verb{\code{s} = 1}, this is monotonicity of its derivative
(e.g., enforcing a convex fit). Set \code{"none"} for no constraints, \code{"i"} for
a monotone increasing constraint, and \code{"d"} for a monotone decreasing
constraint. Using \code{"i"} constrains the basis functions to have positive
coefficients in the fit, and \code{"d"} constrains the basis functions to have
negative coefficients.}

\item{.}{Just like with \code{formula}, \code{.} as in \code{h(.)} or \code{h(.,.)} is
treated as a wildcard variable that generates terms using all variables in
the data. The argument \code{.} should be a character vector of variable
names that \code{.} iterates over. Specifically,
\code{h(., k=1, . = c("W1", "W2", "W3"))} is equivalent to
\code{h(W1, k=1) + h(W2, k=1) + h(W3, k=1)}, and
\code{h(., .,  k=1, . = c("W1", "W2", "W3"))} is equivalent to
\code{h(W1,W2, k=1) + h(W2,W3, k=1) + h(W1, W3, k=1)}}

\item{dot_args_as_string}{Whether the arguments \code{...} are characters or
character vectors and should thus be evaluated directly. When \code{TRUE}, the
expression h("W1", "W2") can be used.}

\item{X}{An optional design matrix where the variables given in \code{...}
can be found. Otherwise, \code{X} is taken from the parent environment.}
}
\description{
HAL Formula term: Generate a single term of the HAL basis
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_pnz}
\alias{calc_pnz}
\title{Calculate Proportion of Nonzero Entries}
\usage{
calc_pnz(X)
}
\description{
Calculate Proportion of Nonzero Entries
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_basis.R
\name{enumerate_edge_basis}
\alias{enumerate_edge_basis}
\title{Enumerate Basis Functions at Generalized Edges}
\usage{
enumerate_edge_basis(
  x,
  max_degree = 3,
  smoothness_orders = rep(0, ncol(x)),
  include_zero_order = FALSE,
  include_lower_order = FALSE
)
}
\arguments{
\item{x}{An input \code{matrix} containing observations and covariates
following standard conventions in problems of statistical learning.}

\item{max_degree}{The highest order of interaction terms for which the basis
functions ought to be generated. The default (\code{NULL}) corresponds to
generating basis functions for the full dimensionality of the input matrix.}

\item{smoothness_orders}{An integer vector of length \code{ncol(x)}
specifying the desired smoothness of the function in each covariate. k = 0
is no smoothness (indicator basis), k = 1 is first order smoothness, and so
on. For an additive model, the component function for each covariate will
have the degree of smoothness as specified by smoothness_orders. For
non-additive components (tensor products of univariate basis functions),
the univariate basis functions in each tensor product have smoothness
degree as specified by smoothness_orders.}

\item{include_zero_order}{A \code{logical}, indicating whether the zeroth
order basis functions are included for each covariate (if \code{TRUE}), in
addition to the smooth basis functions given by \code{smoothness_orders}.
This allows the algorithm to data-adaptively choose the appropriate degree
of smoothness.}

\item{include_lower_order}{A \code{logical}, like \code{include_zero_order},
except including all basis functions of lower smoothness degrees than
specified via \code{smoothness_orders}.}
}
\description{
For degrees of smoothness greater than 1, we must generate the lower order
smoothness basis functions using the knot points at the "edge" of the
hypercube. For example, consider f(x) = x^2 + x, which is second-order
smooth, but will not be generated by purely quadratic basis functions. We
also need to include the y = x function (which corresponds to first-order
HAL basis functions at the left most value/edge of x).
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula_hal9001.R
\name{formula_helpers}
\alias{formula_helpers}
\alias{fill_dots_helper}
\alias{fill_dots}
\title{Formula Helpers}
\usage{
fill_dots_helper(var_names, .)

fill_dots(var_names, .)
}
\arguments{
\item{var_names}{A \code{character} vector of variable names.}

\item{.}{Specification of variables for use in the formula.}
}
\description{
Formula Helpers
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sl_hal9001.R
\name{predict.SL.hal9001}
\alias{predict.SL.hal9001}
\title{predict.SL.hal9001}
\usage{
\method{predict}{SL.hal9001}(object, newdata, ...)
}
\arguments{
\item{object}{A fitted object of class \code{hal9001}.}

\item{newdata}{A matrix of new observations on which to obtain predictions.}

\item{...}{Not used.}
}
\value{
A \code{numeric} vector of predictions from a \code{SL.hal9001}
object based on the provide \code{newdata}.
}
\description{
Predict method for objects of class \code{SL.hal9001}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hal9000.R
\name{hal9000}
\alias{hal9000}
\title{HAL 9000 Quotes}
\usage{
hal9000()
}
\description{
Prints a quote from the HAL 9000 robot from 2001: A Space Odyssey
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_lasso_early_stopping.R
\name{cv_lasso_early_stopping}
\alias{cv_lasso_early_stopping}
\title{Cross-validated LASSO on Indicator Bases}
\usage{
cv_lasso_early_stopping(x_basis, y, n_lambda = 100, n_folds = 10)
}
\arguments{
\item{x_basis}{A \code{dgCMatrix} object corresponding to a sparse matrix of
the basis functions generated for the HAL algorithm.}

\item{y}{A \code{numeric} vector of the observed outcome variable values.}

\item{n_lambda}{A \code{numeric} scalar indicating the number of values of
the L1 regularization parameter (lambda) to be obtained from fitting the
LASSO to the full data. Cross-validation is used to select an optimal
lambda (that minimizes the risk) from among these.}

\item{n_folds}{A \code{numeric} scalar for the number of folds to be used in
the cross-validation procedure to select an optimal value of lambda.}
}
\description{
Fits the LASSO regression using a customized procedure with cross-validation
based on \pkg{origami}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_basis.R
\name{basis_of_degree}
\alias{basis_of_degree}
\title{Compute Degree of Basis Functions}
\usage{
basis_of_degree(
  x,
  degree,
  smoothness_orders,
  include_zero_order,
  include_lower_order
)
}
\arguments{
\item{x}{An input \code{matrix} containing observations and covariates
following standard conventions in problems of statistical learning.}

\item{degree}{The highest order of interaction terms for which the basis
functions ought to be generated. The default (\code{NULL}) corresponds to
generating basis functions for the full dimensionality of the input matrix.}

\item{smoothness_orders}{An integer vector of length \code{ncol(x)}
specifying the desired smoothness of the function in each covariate. k = 0
is no smoothness (indicator basis), k = 1 is first order smoothness, and so
on. For an additive model, the component function for each covariate will
have the degree of smoothness as specified by smoothness_orders. For
non-additive components (tensor products of univariate basis functions),
the univariate basis functions in each tensor product have smoothness
degree as specified by smoothness_orders.}

\item{include_zero_order}{A \code{logical}, indicating whether the zeroth
order basis functions are included for each covariate (if \code{TRUE}), in
addition to the smooth basis functions given by \code{smoothness_orders}.
This allows the algorithm to data-adaptively choose the appropriate degree
of smoothness.}

\item{include_lower_order}{A \code{logical}, like \code{include_zero_order},
except including all basis functions of lower smoothness degrees than
specified via \code{smoothness_orders}.}
}
\value{
A \code{list} containing  basis functions and cutoffs generated from
a set of input columns up to a particular pre-specified degree.
}
\description{
Find the full list of basis functions up to a particular degree
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula_hal9001.R
\name{print.formula_hal9001}
\alias{print.formula_hal9001}
\title{Print formula_hal9001 object}
\usage{
\method{print}{formula_hal9001}(x, ...)
}
\arguments{
\item{x}{A formula_hal9001 object.}

\item{...}{Other arguments (ignored).}
}
\description{
Print formula_hal9001 object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{as_dgCMatrix}
\alias{as_dgCMatrix}
\title{Fast Coercion to Sparse Matrix}
\usage{
as_dgCMatrix(XX_)
}
\arguments{
\item{XX_}{An object of class \code{Matrix} that has a sparse structure
suitable for coercion to a sparse matrix format of \code{dgCMatrix}.}
}
\value{
An object of class \code{dgCMatrix}, coerced from input \code{XX_}.
}
\description{
Fast and efficient coercion of standard matrix objects to sparse matrices.
Borrowed from http://gallery.rcpp.org/articles/sparse-matrix-coercion/.
INTERNAL USE ONLY.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sl_hal9001.R
\name{SL.hal9001}
\alias{SL.hal9001}
\title{Wrapper for Classic SuperLearner}
\usage{
SL.hal9001(
  Y,
  X,
  newX = NULL,
  family = stats::gaussian(),
  obsWeights = rep(1, length(Y)),
  id = NULL,
  max_degree = ifelse(ncol(X) >= 20, 2, 3),
  smoothness_orders = 1,
  num_knots = ifelse(smoothness_orders >= 1, 25, 50),
  reduce_basis = 1/sqrt(length(Y)),
  lambda = NULL,
  ...
)
}
\arguments{
\item{Y}{A \code{numeric} vector of observations of the outcome variable.}

\item{X}{An input \code{matrix} with dimensions number of observations -by-
number of covariates that will be used to derive the design matrix of basis
functions.}

\item{newX}{A matrix of new observations on which to obtain predictions. The
default of \code{NULL} computes predictions on training inputs \code{X}.}

\item{family}{A \code{\link[stats]{family}} object (one that is supported
by \code{\link[glmnet]{glmnet}}) specifying the error/link family for a
generalized linear model.}

\item{obsWeights}{A \code{numeric} vector of observational-level weights.}

\item{id}{A \code{numeric} vector of IDs.}

\item{max_degree}{The highest order of interaction terms for which basis
functions ought to be generated.}

\item{smoothness_orders}{An \code{integer} vector of length 1 or greater,
specifying the smoothness of the basis functions. See the argument
\code{smoothness_orders} of \code{\link{fit_hal}} for more information.}

\item{num_knots}{An \code{integer} vector of length 1 or \code{max_degree},
specifying the maximum number of knot points (i.e., bins) for each
covariate for generating basis functions. See \code{num_knots} argument in
\code{\link{fit_hal}} for more information.}

\item{reduce_basis}{A \code{numeric} value bounded in the open unit interval
indicating the minimum proportion of 1's in a basis function column needed
for the basis function to be included in the procedure to fit the lasso.
Any basis functions with a lower proportion of 1's than the cutoff will be
removed.}

\item{lambda}{A user-specified sequence of values of the regularization
parameter for the lasso L1 regression. If \code{NULL}, the default sequence
in \code{\link[glmnet]{cv.glmnet}} will be used. The cross-validated
optimal value of this regularization parameter will be selected with
\code{\link[glmnet]{cv.glmnet}}.}

\item{...}{Not used.}
}
\value{
An object of class \code{SL.hal9001} with a fitted \code{hal9001}
object and corresponding predictions based on the input data.
}
\description{
Wrapper for \pkg{SuperLearner} for objects of class \code{hal9001}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reduce_basis_filter.R
\name{make_reduced_basis_map}
\alias{make_reduced_basis_map}
\title{Mass-based reduction of basis functions}
\usage{
make_reduced_basis_map(x_basis, reduce_basis_crit)
}
\arguments{
\item{x_basis}{A matrix of basis functions with all redundant basis
functions already removed.}

\item{reduce_basis_crit}{A scalar \code{numeric} value bounded in the open
interval (0,1) indicating the minimum proportion of 1's in a basis function
column needed for the basis function to be included in the procedure to fit
the Lasso. Any basis functions with a lower proportion of 1's than the
specified cutoff will be removed. This argument defaults to \code{NULL}, in
which case all basis functions are used in the lasso-fitting stage of the
HAL algorithm.}
}
\value{
A binary \code{numeric} vector indicating which columns of the
matrix of basis functions to keep (given a one) and which to discard (given
a zero).
}
\description{
A helper function that finds which basis functions to keep (and equivalently
which to discard) based on the proportion of 1's (observations, i.e.,
"mass") included in a given basis function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hal.R
\name{fit_hal}
\alias{fit_hal}
\title{HAL: The Highly Adaptive Lasso}
\usage{
fit_hal(
  X,
  Y,
  formula = NULL,
  X_unpenalized = NULL,
  max_degree = ifelse(ncol(X) >= 20, 2, 3),
  smoothness_orders = 1,
  num_knots = num_knots_generator(max_degree = max_degree, smoothness_orders =
    smoothness_orders, base_num_knots_0 = 200, base_num_knots_1 = 50),
  reduce_basis = 1/sqrt(length(Y)),
  family = c("gaussian", "binomial", "poisson", "cox"),
  lambda = NULL,
  id = NULL,
  offset = NULL,
  fit_control = list(cv_select = TRUE, n_folds = 10, foldid = NULL, use_min = TRUE,
    lambda.min.ratio = 1e-04, prediction_bounds = "default"),
  basis_list = NULL,
  return_lasso = TRUE,
  return_x_basis = FALSE,
  yolo = FALSE
)
}
\arguments{
\item{X}{An input \code{matrix} with dimensions number of observations -by-
number of covariates that will be used to derive the design matrix of basis
functions.}

\item{Y}{A \code{numeric} vector of observations of the outcome variable.}

\item{formula}{A character string formula to be used in
\code{\link{formula_hal}}. See its documentation for details.}

\item{X_unpenalized}{An input \code{matrix} with the same number of rows as
\code{X}, for which no L1 penalization will be performed. Note that
\code{X_unpenalized} is directly appended to the design matrix; no basis
expansion is performed on \code{X_unpenalized}.}

\item{max_degree}{The highest order of interaction terms for which basis
functions ought to be generated.}

\item{smoothness_orders}{An \code{integer}, specifying the smoothness of the
basis functions. See details for \code{smoothness_orders} for more
information.}

\item{num_knots}{An \code{integer} vector of length 1 or \code{max_degree},
specifying the maximum number of knot points (i.e., bins) for any covariate
for generating basis functions. If \code{num_knots} is a unit-length
vector, then the same \code{num_knots} are used for each degree (this is
not recommended). The default settings for \code{num_knots} are
recommended, and these defaults decrease \code{num_knots} with increasing
\code{max_degree} and \code{smoothness_orders}, which prevents (expensive)
combinatorial explosions in the number of higher-degree and higher-order
basis functions generated. This allows the complexity of the optimization
problem to grow scalably. See details of \code{num_knots} more information.}

\item{reduce_basis}{A \code{numeric} value bounded in the open unit interval
indicating the minimum proportion of 1's in a basis function column needed
for the basis function to be included in the procedure to fit the lasso.
Any basis functions with a lower proportion of 1's than the cutoff will be
removed. When \code{reduce_basis} is set to \code{NULL}, all basis
functions are used in the lasso-fitting stage of \code{fit_hal}.}

\item{family}{A \code{character} or a \code{\link[stats]{family}} object
(supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
family for a generalized linear model. \code{character} options are limited
to "gaussian" for fitting a standard penalized linear model, "binomial" for
penalized logistic regression, "poisson" for penalized Poisson regression,
and "cox" for a penalized proportional hazards model. Note that passing in
family objects leads to slower performance relative to passing in a
character family (if supported). For example, one should set
\code{family = "binomial"} instead of \code{family = binomial()} when
calling \code{fit_hal}.}

\item{lambda}{User-specified sequence of values of the regularization
parameter for the lasso L1 regression. If \code{NULL}, the default sequence
in \code{\link[glmnet]{cv.glmnet}} will be used. The cross-validated
optimal value of this regularization parameter will be selected with
\code{\link[glmnet]{cv.glmnet}}. If \code{fit_control}'s \code{cv_select}
argument is set to \code{FALSE}, then the lasso model will be fit via
\code{\link[glmnet]{glmnet}}, and regularized coefficient values for each
lambda in the input array will be returned.}

\item{id}{A vector of ID values that is used to generate cross-validation
folds for \code{\link[glmnet]{cv.glmnet}}. This argument is ignored when
\code{fit_control}'s \code{cv_select} argument is \code{FALSE}.}

\item{offset}{a vector of offset values, used in fitting.}

\item{fit_control}{List of arguments for fitting. Includes the following
arguments, and any others to be passed to \code{\link[glmnet]{cv.glmnet}}
or \code{\link[glmnet]{glmnet}}.
\itemize{
\item \code{cv_select}: A \code{logical} specifying if the sequence of
specified \code{lambda} values should be passed to
\code{\link[glmnet]{cv.glmnet}} in order for a single, optimal value of
\code{lambda} to be selected according to cross-validation. When
\code{cv_select = FALSE}, a \code{\link[glmnet]{glmnet}} model will be
used to fit the sequence of (or single) \code{lambda}.
\item \code{n_folds}: Integer for the number of folds to be used when splitting
the data for V-fold cross-validation. Only used when
\code{cv_select = TRUE}.
\item \code{foldid}: An optional \code{numeric} containing values between 1 and
\code{n_folds}, identifying the fold to which each observation is
assigned. If supplied, \code{n_folds} can be missing. In such a case,
this vector is passed directly to \code{\link[glmnet]{cv.glmnet}}. Only
used when \code{cv_select = TRUE}.
\item \code{use_min}: Specify the choice of lambda to be selected by
\code{\link[glmnet]{cv.glmnet}}. When \code{TRUE}, \code{"lambda.min"} is
used; otherwise, \code{"lambda.1se"}. Only used when
\code{cv_select = TRUE}.
\item \code{lambda.min.ratio}: A \code{\link[glmnet]{glmnet}} argument specifying
the smallest value for \code{lambda}, as a fraction of \code{lambda.max},
the (data derived) entry value (i.e. the smallest value for which all
coefficients are zero). We've seen that not setting \code{lambda.min.ratio}
can lead to no \code{lambda} values that fit the data sufficiently well.
\item \code{prediction_bounds}: A vector of size two that provides the lower and
upper bounds for predictions. When \code{prediction_bounds = "default"},
the predictions are bounded between \code{min(Y) - sd(Y)} and
\code{max(Y) + sd(Y)}. Bounding ensures that there is no extrapolation,
and it is necessary for cross-validation selection and/or Super Learning.
}}

\item{basis_list}{The full set of basis functions generated from \code{X}.}

\item{return_lasso}{A \code{logical} indicating whether or not to return
the \code{\link[glmnet]{glmnet}} fit object of the lasso model.}

\item{return_x_basis}{A \code{logical} indicating whether or not to return
the matrix of (possibly reduced) basis functions used in \code{fit_hal}.}

\item{yolo}{A \code{logical} indicating whether to print one of a curated
selection of quotes from the HAL9000 computer, from the critically
acclaimed epic science-fiction film "2001: A Space Odyssey" (1968).}
}
\value{
Object of class \code{hal9001}, containing a list of basis
functions, a copy map, coefficients estimated for basis functions, and
timing results (for assessing computational efficiency).
}
\description{
Estimation procedure for HAL, the Highly Adaptive Lasso
}
\details{
The procedure uses a custom C++ implementation to generate a design
matrix of spline basis functions of covariates and interactions of
covariates. The lasso regression is fit to this design matrix via
\code{\link[glmnet]{cv.glmnet}} or a custom implementation derived from
\pkg{origami}. The maximum dimension of the design matrix is \eqn{n} -by-
\eqn{(n * 2^(d-1))}, where where \eqn{n} is the number of observations and
\eqn{d} is the number of covariates.

For \code{smoothness_orders = 0}, only zero-order splines (piece-wise
constant) are generated, which assume the true regression function has no
smoothness or continuity. When \code{smoothness_orders = 1}, first-order
splines (piece-wise linear) are generated, which assume continuity of the
true regression function. When \code{smoothness_orders = 2}, second-order
splines (piece-wise quadratic and linear terms) are generated, which assume
a the true regression function has a single order of differentiability.

\code{num_knots} argument specifies the number of knot points for each
covariate and for each \code{max_degree}. Fewer knot points can
significantly decrease runtime, but might be overly simplistic. When
considering \code{smoothness_orders = 0}, too few knot points (e.g., < 50)
can significantly reduce performance. When \code{smoothness_orders = 1} or
higher, then fewer knot points (e.g., 10-30) is actually better for
performance. We recommend specifying \code{num_knots} with respect to
\code{smoothness_orders}, and as a vector of length \code{max_degree} with
values decreasing exponentially. This prevents combinatorial explosions in
the number of higher-degree basis functions generated. The default behavior
of \code{num_knots} follows this logic --- for \code{smoothness_orders = 0},
\code{num_knots} is set to \eqn{500 / 2^{j-1}}, and for
\code{smoothness_orders = 1} or higher, \code{num_knots} is set to
\eqn{200 / 2^{j-1}}, where \eqn{j} is the interaction degree. We also
include some other suitable settings for \code{num_knots} below, all of
which are less complex than default \code{num_knots} and will thus result
in a faster runtime:
\itemize{
\item Some good settings for little to no cost in performance:
\itemize{
\item If \code{smoothness_orders = 0} and \code{max_degree = 3},
\code{num_knots = c(400, 200, 100)}.
\item If \code{smoothness_orders = 1+} and \code{max_degree = 3},
\code{num_knots = c(100, 75, 50)}.
}
\item Recommended settings for fairly fast runtime:
\itemize{
\item If \code{smoothness_orders = 0} and \code{max_degree = 3},
\code{num_knots = c(200, 100, 50)}.
\item If \code{smoothness_orders = 1+} and \code{max_degree = 3},
\code{num_knots = c(50, 25, 15)}.
}
\item Recommended settings for fast runtime:
\itemize{
\item If \code{smoothness_orders = 0} and \code{max_degree = 3},
\code{num_knots = c(100, 50, 25)}.
\item If \code{smoothness_orders = 1+} and \code{max_degree = 3},
\code{num_knots = c(40, 15, 10)}.
}
\item Recommended settings for very fast runtime:
\itemize{
\item If \code{smoothness_orders = 0} and \code{max_degree = 3},
\code{num_knots = c(50, 25, 10)}.
\item If \code{smoothness_orders = 1+} and \code{max_degree = 3},
\code{num_knots = c(25, 10, 5)}.
}
}
}
\examples{
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
y <- rbinom(n = n, size = 1, prob = y_prob)
hal_fit <- fit_hal(X = x, Y = y, family = "binomial")
preds <- predict(hal_fit, new_data = x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lassi.R
\name{predict.lassi}
\alias{predict.lassi}
\title{Predict Method for Lasso on Indicator Bases}
\usage{
\method{predict}{lassi}(fit, new_x_basis, lambdas = NULL)
}
\arguments{
\item{fit}{...}

\item{new_x_basis}{...}

\item{lambdas}{...}
}
\description{
Predict Method for Lasso on Indicator Bases
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{index_first_copy}
\alias{index_first_copy}
\title{Find Copies of Columns}
\usage{
index_first_copy(X)
}
\arguments{
\item{X}{Sparse matrix containing columns of indicator functions.}
}
\description{
Index vector that, for each column in X, indicates the index of the first
copy of that column
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{evaluate_basis}
\alias{evaluate_basis}
\title{Generate Basis Functions}
\usage{
evaluate_basis(basis, X, x_basis, basis_col)
}
\arguments{
\item{basis}{The basis function.}

\item{X}{The design matrix, containing the original data.}

\item{x_basis}{The HAL design matrix, containing indicator functions.}

\item{basis_col}{Numeric indicating which column to populate.}
}
\description{
Populates a column (indexed by basis_col) of x_basis with basis indicators.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{meets_basis}
\alias{meets_basis}
\title{Compute Values of Basis Functions}
\usage{
meets_basis(X, row_num, cols, cutoffs, orders)
}
\arguments{
\item{X}{The design matrix, containing the original data.}

\item{row_num}{Numeri for  a row index over which to evaluate.}

\item{cols}{Numeric for the column indices of the basis function.}

\item{cutoffs}{Numeric providing thresholds.}

\item{orders}{Numeric providing smoothness orders}
}
\description{
Computes and returns the indicator value for the basis described by
cols and cutoffs for a given row of X
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{make_basis_list}
\alias{make_basis_list}
\title{Sort Basis Functions}
\usage{
make_basis_list(X_sub, cols, order_map)
}
\arguments{
\item{X_sub}{A subset of the columns of X, the original design matrix.}

\item{cols}{An index of the columns that were reduced to by sub-setting.}

\item{order_map}{A vector with length the original unsubsetted matrix X which specifies the smoothness of the function in each covariate.}
}
\description{
Build a sorted list of unique basis functions based on columns, where each
basis function is a list
}
\details{
Note that sorting of columns is performed such that the basis order
equals cols.length() and each basis function is a list(cols, cutoffs).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula_hal9001.R
\name{+.formula_hal9001}
\alias{+.formula_hal9001}
\title{HAL Formula addition: Adding formula term object together into a single
formula object term.}
\usage{
\method{+}{formula_hal9001}(x, y)
}
\arguments{
\item{x}{A \code{formula_hal9001} object as outputted by \code{h}.}

\item{y}{A \code{formula_hal9001} object as outputted by \code{h}.}
}
\description{
HAL Formula addition: Adding formula term object together into a single
formula object term.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hal_quotes.R
\docType{data}
\name{hal_quotes}
\alias{hal_quotes}
\title{HAL9000 Quotes from "2001: A Space Odyssey"}
\format{
A vector of quotes.
}
\usage{
hal_quotes
}
\description{
Curated selection of quotes from the HAL9000 computer, from the critically
acclaimed epic science-fiction film "2001: A Space Odyssey" (1968).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_dupes.R
\name{make_copy_map}
\alias{make_copy_map}
\title{Build Copy Maps}
\usage{
make_copy_map(x_basis)
}
\arguments{
\item{x_basis}{A design matrix consisting of basis (indicator) functions for
covariates (X) and terms for interactions thereof.}
}
\value{
A \code{list} of \code{numeric} vectors indicating indices of basis
functions that are identical in the training set.
}
\description{
Build Copy Maps
}
\examples{
\donttest{
gendata <- function(n) {
  W1 <- runif(n, -3, 3)
  W2 <- rnorm(n)
  W3 <- runif(n)
  W4 <- rnorm(n)
  g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
  A <- rbinom(n, 1, g0)
  Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
  Y <- rbinom(n, 1, Q0)
  data.frame(A, W1, W2, W3, W4, Y)
}
set.seed(1234)
data <- gendata(100)
covars <- setdiff(names(data), "Y")
X <- as.matrix(data[, covars, drop = FALSE])
basis_list <- enumerate_basis(X)
x_basis <- make_design_matrix(X, basis_list)
copy_map <- make_copy_map(x_basis)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula_hal9001.R
\name{formula_hal}
\alias{formula_hal}
\title{HAL Formula: Convert formula or string to \code{formula_HAL} object.}
\usage{
formula_hal(formula, smoothness_orders, num_knots, X = NULL)
}
\arguments{
\item{formula}{A \code{formula_hal9001} object as outputted by \code{h}.}

\item{smoothness_orders}{A default value for \code{s} if not provided
explicitly to the function \code{h}.}

\item{num_knots}{A default value for \code{k} if not provided explicitly to
the function \code{h}.}

\item{X}{Controls inheritance of the variable \code{X} from parent environment.
When \code{NULL} (the default), such a variable is inherited.}
}
\description{
HAL Formula: Convert formula or string to \code{formula_HAL} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hal.R
\name{num_knots_generator}
\alias{num_knots_generator}
\title{A default generator for the \code{num_knots} argument for each degree of
interactions and the smoothness orders.}
\usage{
num_knots_generator(
  max_degree,
  smoothness_orders,
  base_num_knots_0 = 500,
  base_num_knots_1 = 200
)
}
\arguments{
\item{smoothness_orders}{see \code{\link{fit_hal}}.}

\item{base_num_knots_0}{The base number of knots for zeroth-order smoothness
basis functions. The number of knots by degree interaction decays as
\code{base_num_knots_0/2^(d-1)} where \code{d} is the interaction degree of the basis
function.}

\item{base_num_knots_1}{The base number of knots for 1 or greater order
smoothness basis functions. The number of knots by degree interaction
decays as \code{base_num_knots_1/2^(d-1)} where \code{d} is the interaction degree of
the basis function.}

\item{d}{interaction degree.}
}
\description{
A default generator for the \code{num_knots} argument for each degree of
interactions and the smoothness orders.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lassi.R
\name{lassi}
\alias{lassi}
\title{Custom Lasso implementation for matrices of indicator functions}
\usage{
lassi(
  x,
  y,
  lambdas = NULL,
  nlambda = 100,
  lambda_min_ratio = 0.01,
  center = FALSE
)
}
\arguments{
\item{x}{The covariate matrix}

\item{y}{The outcome vector}

\item{lambdas}{A sequence of values for the L1 regularization parameter
(lambda) to be used in fitting the LASSO. Defaults to \code{NULL}.}

\item{nlambda}{number of lambdas to fit.}

\item{lambda_min_ratio}{ratio of largest to smallest lambda to fit.}

\item{center}{...}
}
\description{
Custom Lasso implementation for matrices of indicator functions
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.hal9001}
\alias{summary.hal9001}
\title{Summary Method for HAL fit objects}
\usage{
\method{summary}{hal9001}(
  object,
  lambda = NULL,
  only_nonzero_coefs = TRUE,
  include_redundant_terms = FALSE,
  round_cutoffs = 3,
  ...
)
}
\arguments{
\item{object}{An object of class \code{hal9001}, containing the results of
fitting the Highly Adaptive Lasso, as produced by \code{\link{fit_hal}}.}

\item{lambda}{Optional \code{numeric} value of the lambda tuning
parameter, for which corresponding coefficient values will be summarized.
Defaults to \code{\link{fit_hal}}'s optimal value, \code{lambda_star}, or
the minimum value of \code{lambda_star}.}

\item{only_nonzero_coefs}{A \code{logical} specifying whether the summary
should include only terms with non-zero coefficients.}

\item{include_redundant_terms}{A \code{logical} specifying whether the
summary should remove so-called "redundant terms". We define a redundant
term (say x1) as a term (1) with basis function corresponding to an
existing basis function, a duplicate; and (2) the duplicate contains the
x1 term as part of its term, so that x1 terms inclusion would be redundant.
For example, say the same coefficient corresponds to these three terms:
(1) "I(age >= 50)*I(bmi >= 18)", (2) "I(age >= 50)", and (3)
"I(education >= 16)". When \code{include_redundant_terms} is
\code{FALSE} (default), the second basis function is omitted.}

\item{round_cutoffs}{An \code{integer} indicating the number of decimal
places to be used for rounding cutoff values in the term. For example, if
"bmi" was numeric that was rounded to the third decimal, in the example
above we would have needed to specify \code{round_cutoffs = 0} in order to
yield a term like "I(bmi >= 18)" opposed to something like
"I(bmi >= 18.111)". This rounding is intended to simplify the term-wise
part of the output and only rounds the basis cutoffs, the \code{hal9001}
model's coefficients are not rounded.}

\item{...}{Additional arguments passed to \code{summary}, not supported.}
}
\value{
A list summarizing a \code{hal9001} object's coefficients.
}
\description{
Summary Method for HAL fit objects
}
\details{
Method for summarizing the coefficients of the Highly Adaptive
Lasso estimator in terms of the basis functions corresponding to covariates
and interactions of covariates, returned as a single S3 object of class
\code{hal9001}.

Due to the nature of the basis function terms, the summary tables can be
extremely wide. The R environment might not be the optimal location to view
the summary. Tables can be exported from R to LaTeX with \pkg{xtable}
package (or similar). Here's an example:
\code{print(xtable(summary(fit)$table, type = "latex"), file = "dt.tex")}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.hal9001}
\alias{predict.hal9001}
\title{Prediction from HAL fits}
\usage{
\method{predict}{hal9001}(
  object,
  new_data,
  new_X_unpenalized = NULL,
  offset = NULL,
  type = c("response", "link"),
  p_reserve = 0.75,
  ...
)
}
\arguments{
\item{object}{An object of class \code{hal9001}, containing the results of
fitting the Highly Adaptive Lasso, as produced by \code{\link{fit_hal}}.}

\item{new_data}{A \code{matrix} or \code{data.frame} containing new data
(i.e., observations not used for fitting the \code{hal9001} object that's
passed in via the \code{object} argument) for which the \code{hal9001}
object will compute predicted values.}

\item{new_X_unpenalized}{If the user supplied \code{X_unpenalized} during
training, then user should also supply this matrix with the same number of
observations as \code{new_data}.}

\item{offset}{A vector of offsets. Must be provided if provided at training.}

\item{type}{Either "response" for predictions of the response, or "link" for
un-transformed predictions (on the scale of the link function).}

\item{p_reserve}{Sparse matrix pre-allocation proportion, which is the
anticipated proportion of 1's in the design matrix. Default value is
recommended in most settings. If a dense design matrix is expected, it
would be useful to set \code{p_reserve} to a higher value.}

\item{...}{Additional arguments passed to \code{predict} as necessary.}
}
\value{
A \code{numeric} vector of predictions from a \code{hal9001} object.
}
\description{
Prediction from HAL fits
}
\details{
Method for computing and extracting predictions from fits of the
Highly Adaptive Lasso estimator, returned as a single S3 objects of class
\code{hal9001}.
}
\note{
This prediction method does not function similarly to the equivalent
method from \pkg{glmnet}. In particular, this procedure will not return a
subset of lambdas originally specified in calling \code{\link{fit_hal}}
nor result in re-fitting. Instead, it will return predictions for all of
the lambdas specified in the call to \code{\link{fit_hal}} that constructs
\code{object}, when \code{fit_control}'s \code{cv_select} is set to
\code{FALSE}. When \code{fit_control}'s \code{cv_select} is set to
\code{TRUE}, predictions will only be returned for the value of lambda
selected by cross-validation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_lasso.R
\name{cv_lasso}
\alias{cv_lasso}
\title{Cross-validated Lasso on Indicator Bases}
\usage{
cv_lasso(x_basis, y, n_lambda = 100, n_folds = 10, center = FALSE)
}
\arguments{
\item{x_basis}{A \code{dgCMatrix} object corresponding to a sparse matrix of
the basis functions generated for the HAL algorithm.}

\item{y}{A \code{numeric} vector of the observed outcome variable values.}

\item{n_lambda}{A \code{numeric} scalar indicating the number of values of
the L1 regularization parameter (lambda) to be obtained from fitting the
Lasso to the full data. Cross-validation is used to select an optimal
lambda (that minimizes the risk) from among these.}

\item{n_folds}{A \code{numeric} scalar for the number of folds to be used in
the cross-validation procedure to select an optimal value of lambda.}

\item{center}{binary. If \code{TRUE}, covariates are centered. This is much
slower, but matches the \code{glmnet} implementation. Default \code{FALSE}.}
}
\description{
Fits Lasso regression using a customized procedure, with cross-validation
based on \pkg{origami}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_basis.R
\name{quantizer}
\alias{quantizer}
\title{Discretize Variables into Number of Bins by Unique Values}
\usage{
quantizer(X, bins)
}
\arguments{
\item{X}{A \code{numeric} vector to be discretized.}

\item{bins}{A \code{numeric} scalar indicating the number of bins into which
\code{X} should be discretized..}
}
\description{
Discretize Variables into Number of Bins by Unique Values
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{lassi_predict}
\alias{lassi_predict}
\title{Prediction from a Lassi Model}
\usage{
lassi_predict(X, beta, intercept)
}
\arguments{
\item{X}{A sparse matrix of HAL basis functions.}

\item{beta}{A vector of coefficient values for the HAL basis functions.}

\item{intercept}{A numeric value giving the intercept of the HAL model.}
}
\description{
Prediction from a Lassi Model
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{make_design_matrix}
\alias{make_design_matrix}
\title{Build HAL Design Matrix}
\usage{
make_design_matrix(X, blist, p_reserve = 0.5)
}
\arguments{
\item{X}{Matrix of covariates containing observed data in the columns.}

\item{blist}{List of basis functions with which to build HAL design matrix.}

\item{p_reserve}{Sparse matrix pre-allocation proportion. Default value is 0.5.
If one expects a dense HAL design matrix, it is useful to set p_reserve to a higher value.}
}
\value{
A \code{dgCMatrix} sparse matrix of indicator basis functions
corresponding to the design matrix in a zero-order highly adaptive lasso.
}
\description{
Make a HAL design matrix based on original design matrix X and a list of
basis functions in argument blist
}
\examples{
\donttest{
gendata <- function(n) {
  W1 <- runif(n, -3, 3)
  W2 <- rnorm(n)
  W3 <- runif(n)
  W4 <- rnorm(n)
  g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
  A <- rbinom(n, 1, g0)
  Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
  Y <- rbinom(n, 1, Q0)
  data.frame(A, W1, W2, W3, W4, Y)
}
set.seed(1234)
data <- gendata(100)
covars <- setdiff(names(data), "Y")
X <- as.matrix(data[, covars, drop = FALSE])
basis_list <- enumerate_basis(X)
x_basis <- make_design_matrix(X, basis_list)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{apply_copy_map}
\alias{apply_copy_map}
\title{Apply copy map}
\usage{
apply_copy_map(X, copy_map)
}
\arguments{
\item{X}{Sparse matrix containing columns of indicator functions.}

\item{copy_map}{the copy map}
}
\value{
A \code{dgCMatrix} sparse matrix corresponding to the design matrix
for a zero-th order highly adaptive lasso, but with all duplicated columns
(basis functions) removed.
}
\description{
OR duplicate training set columns together
}
\examples{
\donttest{
gendata <- function(n) {
  W1 <- runif(n, -3, 3)
  W2 <- rnorm(n)
  W3 <- runif(n)
  W4 <- rnorm(n)
  g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
  A <- rbinom(n, 1, g0)
  Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
  Y <- rbinom(n, 1, Q0)
  data.frame(A, W1, W2, W3, W4, Y)
}
set.seed(1234)
data <- gendata(100)
covars <- setdiff(names(data), "Y")
X <- as.matrix(data[, covars, drop = FALSE])
basis_list <- enumerate_basis(X)
x_basis <- make_design_matrix(X, basis_list)
copy_map <- make_copy_map(x_basis)
x_basis_uniq <- apply_copy_map(x_basis, copy_map)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_basis.R
\name{basis_list_cols}
\alias{basis_list_cols}
\title{List Basis Functions}
\usage{
basis_list_cols(
  cols,
  x,
  smoothness_orders,
  include_zero_order,
  include_lower_order = FALSE
)
}
\arguments{
\item{cols}{Index or indices (as \code{numeric}) of covariates (columns) of
interest in the data matrix \code{x} for which basis functions ought to be
generated. Note that basis functions for interactions of these columns are
computed automatically.}

\item{x}{A \code{matrix} containing observations in the rows and covariates
in the columns. Basis functions are computed for these covariates.}

\item{smoothness_orders}{An integer vector of length \code{ncol(x)}
specifying the desired smoothness of the function in each covariate. k = 0
is no smoothness (indicator basis), k = 1 is first order smoothness, and so
on. For an additive model, the component function for each covariate will
have the degree of smoothness as specified by smoothness_orders. For
non-additive components (tensor products of univariate basis functions),
the univariate basis functions in each tensor product have smoothness
degree as specified by smoothness_orders.}

\item{include_zero_order}{A \code{logical}, indicating whether the zeroth
order basis functions are included for each covariate (if \code{TRUE}), in
addition to the smooth basis functions given by \code{smoothness_orders}.
This allows the algorithm to data-adaptively choose the appropriate degree
of smoothness.}

\item{include_lower_order}{A \code{logical}, like \code{include_zero_order},
except including all basis functions of lower smoothness degrees than
specified via \code{smoothness_orders}.}
}
\value{
A \code{list} containing the basis functions generated from a set of
input columns.
}
\description{
Build a list of basis functions from a set of columns
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lassi.R
\name{lassi_fit_module}
\alias{lassi_fit_module}
\title{Rcpp module: lassi_fit_module}
\description{
Rcpp module: lassi_fit_module
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_xscale}
\alias{calc_xscale}
\title{Calculating Centered and Scaled Matrices}
\usage{
calc_xscale(X, xcenter)
}
\arguments{
\item{X}{A sparse matrix, to be centered.}

\item{xcenter}{A vector of column means to be used for centering X.}
}
\description{
Calculating Centered and Scaled Matrices
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hal9001-package.R
\name{hal9001}
\alias{hal9001}
\title{hal9001}
\description{
Package for fitting the Highly Adaptive LASSO (HAL) estimator
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_lasso.R
\name{lassi_origami}
\alias{lassi_origami}
\title{Single Lasso estimation for cross-validation with Origami}
\usage{
lassi_origami(fold, data, lambdas, center = FALSE)
}
\arguments{
\item{fold}{A \code{fold} object produced by a call to \code{make_folds}
from the \pkg{origami}.}

\item{data}{A \code{dgCMatrix} object containing the outcome values (Y) in
its first column and vectors corresponding to the basis functions of HAL in
all other columns. Consult the description of HAL regression for details.}

\item{lambdas}{A \code{numeric} vector corresponding to a sequence of lambda
values obtained by fitting the Lasso on the full data.}

\item{center}{binary. If \code{TRUE}, covariates are centered. This is much
slower, but matches the \code{glmnet} implementation. Default \code{FALSE}.}
}
\description{
Fits Lasso regression over a single fold of a cross-validated data set. This
is meant to be called using \code{\link[origami]{cross_validate}}, which is
done through \code{\link{cv_lasso}}. Note that this procedure is NOT meant
to be invoked by itself. INTERNAL USE ONLY.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_basis.R
\name{enumerate_basis}
\alias{enumerate_basis}
\title{Enumerate Basis Functions}
\usage{
enumerate_basis(
  x,
  max_degree = NULL,
  smoothness_orders = rep(0, ncol(x)),
  include_zero_order = FALSE,
  include_lower_order = FALSE,
  num_knots = NULL
)
}
\arguments{
\item{x}{An input \code{matrix} containing observations and covariates
following standard conventions in problems of statistical learning.}

\item{max_degree}{The highest order of interaction terms for which the basis
functions ought to be generated. The default (\code{NULL}) corresponds to
generating basis functions for the full dimensionality of the input matrix.}

\item{smoothness_orders}{An integer vector of length \code{ncol(x)}
specifying the desired smoothness of the function in each covariate. k = 0
is no smoothness (indicator basis), k = 1 is first order smoothness, and so
on. For an additive model, the component function for each covariate will
have the degree of smoothness as specified by smoothness_orders. For
non-additive components (tensor products of univariate basis functions),
the univariate basis functions in each tensor product have smoothness
degree as specified by smoothness_orders.}

\item{include_zero_order}{A \code{logical}, indicating whether the zeroth
order basis functions are included for each covariate (if \code{TRUE}), in
addition to the smooth basis functions given by \code{smoothness_orders}.
This allows the algorithm to data-adaptively choose the appropriate degree
of smoothness.}

\item{include_lower_order}{A \code{logical}, like \code{include_zero_order},
except including all basis functions of lower smoothness degrees than
specified via \code{smoothness_orders}.}

\item{num_knots}{A vector of length \code{max_degree}, which determines how
granular the knot points to generate basis functions should be for each
degree of basis function. The first entry of \code{num_knots} determines
the number of knot points to be used for each univariate basis function.
More generally, The kth entry of \code{num_knots} determines the number of
knot points to be used for the kth degree basis functions. Specifically,
for a kth degree basis function, which is the tensor product of k
univariate basis functions, this determines the number of knot points to be
used for each univariate basis function in the tensor product.}
}
\value{
A \code{list} of basis functions generated for all covariates and
interaction thereof up to a pre-specified degree.
}
\description{
Generate basis functions for all covariates and interaction terms thereof up
to a specified order/degree.
}
\examples{
\donttest{
gendata <- function(n) {
  W1 <- runif(n, -3, 3)
  W2 <- rnorm(n)
  W3 <- runif(n)
  W4 <- rnorm(n)
  g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
  A <- rbinom(n, 1, g0)
  Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
  Y <- rbinom(n, 1, Q0)
  data.frame(A, W1, W2, W3, W4, Y)
}
set.seed(1234)
data <- gendata(100)
covars <- setdiff(names(data), "Y")
X <- as.matrix(data[, covars, drop = FALSE])
basis_list <- enumerate_basis(X)
}

}
