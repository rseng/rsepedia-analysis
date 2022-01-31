
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`txshift`

<!-- badges: start -->

[![R-CMD-check](https://github.com/nhejazi/txshift/workflows/R-CMD-check/badge.svg)](https://github.com/nhejazi/txshift/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/txshift/master.svg)](https://codecov.io/github/nhejazi/txshift?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/txshift)](https://www.r-pkg.org/pkg/txshift)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/txshift)](https://CRAN.R-project.org/package=txshift)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/txshift)](https://CRAN.R-project.org/package=txshift)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447)
<!-- badges: end -->

> Efficient Estimation of the Causal Effects of Stochastic Interventions

**Authors:** [Nima Hejazi](https://nimahejazi.org) and [David
Benkeser](https://www.sph.emory.edu/faculty/profile/#!dbenkes)

-----

## What’s `txshift`?

The `txshift` R package is designed to provide facilities for the
construction of efficient estimators of the counterfactual mean of an
outcome under stochastic interventions that depend on the natural value
of treatment (Dı́az and van der Laan 2012; Haneuse and Rotnitzky 2013).
`txshift`implements and builds upon a simplified algorithm for the
targeted maximum likelihood (TML) estimator of such a causal parameter,
originally proposed by Dı́az and van der Laan (2018), and makes use of
analogous machinery to compute an efficient one-step estimator (Pfanzagl
and Wefelmeyer 1985). `txshift` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) (Coyle, Hejazi, Malenica, et
al. 2022) to allow for ensemble machine learning to be leveraged in the
estimation procedure.

For many practical applications (e.g., vaccine efficacy trials),
observed data is often subject to a two-phase sampling mechanism (i.e.,
through the use of a two-stage design). In such cases, efficient
estimators (of both varieties) must be augmented to construct unbiased
estimates of the population-level causal parameter. Rose and van der
Laan (2011) first introduced an augmentation procedure that relies on
introducing inverse probability of censoring (IPC) weights directly to
an appropriate loss function or to the efficient influence function
estimating equation. `txshift` extends this approach to compute
IPC-weighted one-step and TML estimators of the counterfactual mean
outcome under a shift stochastic treatment regime. The package is
designed to implement the statistical methodology described in Hejazi et
al. (2020) and extensions thereof.

-----

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=txshift) via

``` r
install.packages("txshift")
```

*Note:* If `txshift` is installed from
[CRAN](https://CRAN.R-project.org/package=txshift), the `sl3`, an
enhancing dependency that allows ensemble machine learning to be used
for nuisance parameter estimation, won’t be included. We highly
recommend additionally installing `sl3` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("tlverse/sl3@master")
```

For the latest features, install the most recent *stable version* of
`txshift` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/txshift@master")
```

To contribute, install the *development version* of `txshift` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/txshift@devel")
```

-----

## Example

To illustrate how `txshift` may be used to ascertain the effect of a
treatment, consider the following example:

``` r
library(txshift)
#> txshift v0.3.7: Efficient Estimation of the Causal Effects of Stochastic
#> Interventions
library(sl3)
set.seed(429153)

# simulate simple data
n_obs <- 500
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))

# now, let's introduce a a two-stage sampling process
C_samp <- rbinom(n_obs, 1, plogis(W + Y))

# fit the full-data TMLE (ignoring two-phase sampling)
tmle <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "tmle",
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: tmle
#> Estimate: 0.7685
#> Std. Error: 0.019
#> 95% CI: [0.7292, 0.8037]

# fit a full-data one-step estimator for comparison (again, no sampling)
os <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "onestep",
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
os
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: onestep
#> Estimate: 0.7685
#> Std. Error: 0.019
#> 95% CI: [0.7292, 0.8037]

# fit an IPCW-TMLE to account for the two-phase sampling process
tmle_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = 0.5, C_samp = C_samp, V = c("W", "Y"),
  estimator = "tmle", max_iter = 5, eif_reg_type = "glm",
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle_ipcw
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: tmle
#> Estimate: 0.7603
#> Std. Error: 0.0204
#> 95% CI: [0.718, 0.798]

# compare with an IPCW-agumented one-step estimator under two-phase sampling
os_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = 0.5, C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep", eif_reg_type = "glm",
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
os_ipcw
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: onestep
#> Estimate: 0.7601
#> Std. Error: 0.0204
#> 95% CI: [0.7178, 0.7979]
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/txshift/issues). Further
details on filing issues are provided in our [contribution
guidelines](https://github.com/nhejazi/txshift/blob/master/CONTRIBUTING.md).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/txshift/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `txshift` R package, please cite the following:

``` 
    @article{hejazi2020efficient,
      author = {Hejazi, Nima S and {van der Laan}, Mark J and Janes, Holly
        E and Gilbert, Peter B and Benkeser, David C},
      title = {Efficient nonparametric inference on the effects of
        stochastic interventions under two-phase sampling, with
        applications to vaccine efficacy trials},
      year = {2020},
      doi = {10.1111/biom.13375},
      url = {https://doi.org/10.1111/biom.13375},
      journal = {Biometrics},
      publisher = {Wiley Online Library}
    }

    @article{hejazi2020txshift-joss,
      author = {Hejazi, Nima S and Benkeser, David C},
      title = {{txshift}: Efficient estimation of the causal effects of
        stochastic interventions in {R}},
      year  = {2020},
      doi = {10.21105/joss.02447},
      url = {https://doi.org/10.21105/joss.02447},
      journal = {Journal of Open Source Software},
      publisher = {The Open Journal}
    }

    @software{hejazi2022txshift-rpkg,
      author = {Hejazi, Nima S and Benkeser, David C},
      title = {{txshift}: Efficient Estimation of the Causal Effects of
        Stochastic Interventions},
      year  = {2022},
      doi = {10.5281/zenodo.4070042},
      url = {https://CRAN.R-project.org/package=txshift},
      note = {R package version 0.3.7}
    }
```

-----

## Related

  - [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R
    package providing an independent implementation of the same core
    routines for the TML estimation procedure and statistical
    methodology as is made available here, through reliance on a unified
    interface for Targeted Learning provided by the
    [`tmle3`](https://github.com/tlverse/tmle3) engine of the [`tlverse`
    ecosystem](https://github.com/tlverse).

  - [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
    providing facilities to estimate the causal effect of stochastic
    treatment regimes in the mediation setting, including classical
    (IPW) and augmented double robust (one-step) estimators. This is an
    implementation of the methodology explored by Dı́az and Hejazi
    (2020).

  - [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
    package for estimating the conditional density treatment mechanism
    component of this parameter based on using the [highly adaptive
    lasso](https://github.com/tlverse/hal9001) (Coyle, Hejazi, Phillips,
    et al. 2022; Hejazi, Coyle, and van der Laan 2020) in combination
    with a pooled hazard regression. This package implements a variant
    of the approach advocated by Dı́az and van der Laan (2011).

-----

## Funding

The development of this software was supported in part through grants
from the National Library of Medicine (award no. [T32
LM012417](https://reporter.nih.gov/project-details/9248418)) and the
National Institute of Allergy and Infectious Diseases (award no. [R01
AI074345](https://reporter.nih.gov/project-details/9926564)) of the
National Institutes of Health, as well as by the National Science
Foundation (award no.
[DMS 2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

-----

## License

© 2017-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2017-2022 Nima S. Hejazi
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

-----

## References

<div id="refs" class="references">

<div id="ref-coyle-sl3-rpkg">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2022. *sl3: Modern Machine Learning Pipelines for Super
Learning*. <https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-coyle-hal9001-rpkg">

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars W van der Laan,
and Mark J van der Laan. 2022. *hal9001: The Scalable Highly Adaptive
Lasso*. <https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-diaz2012population">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2): 541–49.

</div>

<div id="ref-diaz2018stochastic">

———. 2018. “Stochastic Treatment Regimes.” In *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*, 167–80.
Springer Science & Business Media.

</div>

<div id="ref-haneuse2013estimation">

Haneuse, Sebastian, and Andrea Rotnitzky. 2013. “Estimation of the
Effect of Interventions That Modify the Received Treatment.” *Statistics
in Medicine* 32 (30): 5260–77.

</div>

<div id="ref-hejazi2020hal9001-joss">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020. “hal9001:
Scalable Highly Adaptive Lasso Regression in R.” *Journal of Open Source
Software* 5 (53): 2526. <https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-hejazi2020efficient">

Hejazi, Nima S, Mark J van der Laan, Holly E Janes, Peter B Gilbert, and
David C Benkeser. 2020. “Efficient Nonparametric Inference on the
Effects of Stochastic Interventions Under Two-Phase Sampling, with
Applications to Vaccine Efficacy Trials.” *Biometrics* 77 (4): 1241–53.
<https://doi.org/10.1111/biom.13375>.

</div>

<div id="ref-pfanzagl1985contributions">

Pfanzagl, J, and W Wefelmeyer. 1985. “Contributions to a General
Asymptotic Statistical Theory.” *Statistics & Risk Modeling* 3 (3-4):
379–88.

</div>

<div id="ref-rose2011targeted2sd">

Rose, Sherri, and Mark J van der Laan. 2011. “A Targeted Maximum
Likelihood Estimator for Two-Stage Designs.” *International Journal of
Biostatistics* 7 (1): 1–21.

</div>

</div>
# txshift 0.3.6

As of October 2021:
* Minor updates to ensure compatibility with v0.4.1 of `hal9001` and v0.2.1 of
  `haldensify`, both recently updated on CRAN.
* Removal of the `LazyData` field from the `DESCRIPTION`, since no `data`
  directory is included with the package.
* Minor tweaks to existing unit tests to remove `rlang` from the `Suggests`
  field of the `DESCRIPTION`.
* Vignettes for the standard and IPCW-augmented estimation procedures have been
  combined to reduce redundancy and reduce build time per CRAN requests.

As of May 2021:
* The use of `hal9001::fit_hal()` internally for evaluation of a conditional
  mean of the full-data EIF has been revised for compatibility with v0.4.0+ of
  the `hal9001` package.
* Defaults passed in through the argument `g_exp_fit_args`, and to the function
  `est_g_exp()`, have been updated for compatibility with v0.1.5+ of the
  `haldensify` package.

As of April 2021:
* The `print()` methods have been updated to remove the use of [`cli`
  functions](https://github.com/r-lib/cli), which, for simplicity, has been
  replaced by the use of `message()`.
* Addition of a hidden slot `.eif_mat` to the `txshift_msm` class, supporting
  export of the matrix of EIF estimates for each shift in `delta_grid`.

# txshift 0.3.5

As of February 2021:
* Remove cross-linking to `sl3` functions as per request from CRAN. This can be
  reversed once `sl3` is available on CRAN.

As of January 2021:
* Simulation experiments testing the performance of the procedures in the
  presence of loss to follow-up censoring indicate that the TML estimator
  outperforms the one-step for the EIF-based two-phase sampling correction.
  Generally, we recommend use of the TML estimator (the default) across all
  settings, though performance of the one-step estimator is much worse.

As of December 2020:
* A `delta` slot has been added to the `txshift` class to record the shift.
* Hidden slots have been similarly added to the `txshift_msm` class.
* The `summary` method has been removed, with the functionality now supported
  by the `print` methods for the `txshift` and `txshift_msm` classes.
* The `plot` method has been amended to support simultaneous confidence bands.

As of October 2020:
* Changes all references to the argument `C` to `C_samp` for the indicator of
  inclusion in the second-stage sample.
* Adds the new argument `C_cens` to denote censoring due to loss to follow-up,
  i.e., prior to the occurrence of the outcome.
* Adds a nuisance regression for censoring `C_cens` and adjusts the estimation
  procedure so as to use inverse censoring weights in the full-data EIF
  procedure (NOTE: these are not updated in the two-phase sampling correction).
* Renaming of arguments to internal functions and functions themselves:
  * From `est_g` to `est_g_exp` for the exposure mechanism density estimation
  * From `est_ipcw` to `est_samp` for the two-phase sampling mechanism
  * Add `est_g_cens` for the loss to follow-up censoring mechanism

# txshift 0.3.4

As of September 2020:
* Moved `sl3` dependency to an `Enhances` designation for CRAN submission.
* As above, removed `sl3` from `Remotes` and added installation safety checks.

As of June 2020:
* Add single-knot spline to MSM summarization (`msm_vimshift`).
* Add class and `plot` method for MSM summarization (`msm_vimshift`).
* Fix bug in `msm_vimshift` for computing CIs for binary outcomes by switching
  from manually computing CIs to internally using custom `confint` method.
* Fix bug in `msm_vimshift` for building `lm` model objects through weighted
  regression; move models from `plot` method to `msm_vimshift`.
* Finish drafting brief paper for _Journal of Open Source Software_.

# txshift 0.3.3

As of April 2020:
* Change export status of internal functions (e.g., no longer exporting
  `onestep_txshift` and `tmle_txshift`).
* Finish adding Roxygen "details" and "return" slots throughout functions.
* Add examples to main estimation functions (`txshift`, `vimshift_msm`).
* Update argument names and add several `assert_that` checks.
* Change `fit_spec` terminology to `fit_ext` for external fits.
* Add unit tests for MSM functionality and nuisance parameter estimation.

As of March 2020:
* Extensive documentation, including fixing estimation terminology (e.g.,
  one-step instead of AIPW) and adding Roxygen "details" and "return" slots.
* Begin adding examples to exported functions.

# txshift 0.3.2

As of March 2020:
* Corrections to dependencies in preparation for eventual CRAN release.
* Change several previously exported functions to internal, including `eif`,
  `est_Hn`, `est_Q`, `est_g`, `est_ipcw`, `fit_fluctuation`, `ipcw_eif_update`).
* Remove/reduce GitHub-only dependencies (now only `sl3`).
* Change title partially (from "Targeted Learning" to "Efficient Estimation").
* Lock dependency versions (e.g., `sl3` >= v1.3.7)
* Extensive documentation updates.

# txshift 0.3.1

As of December 2019:
* Changes arguments of `hal9001::fit_hal` in pseudo-outcome regression for
    efficient estimation by explicitly including `max_degree = NULL`.
* Change to TMLE convergence criterion: use a less strict criterion such that
     | Pn D | \leq sigma / (sqrt(n) \cdot max(10, log(n))) instead of \leq 1/n.
    Empirical studies suggest this curbs issues addressed by over-agressive
    updates from the targeting step.
* Remove pinning of `sl3` dependency to a specific tag (formerly v1.2.0).
* Lock dependency version: `sl3` >= v1.3.6 and `hal9001` >= v0.2.5.

# txshift 0.3.0

As of October 2019:
* Change use of `as.data.table` to `data.table` in internal functions to catch
    up with changes in dependencies.

# txshift 0.2.9

As of September 2019:
* Remove errant intercept term and lower iterations for fluctuation models.
* Change weighting scheme in marginal structural model summarization to weight
    all estimates identically rather than by inverse variance as a default.
* Updates to documentation.

# txshift 0.2.8

As of September 2019:
* Add safety checks for convergence of fluctuation regressions based on those
    appearing in `drtmle` and/or `survtmle`.
* Change default confidence interval type to use marginal CIs across multiple
    parameters instead of a simultaneous confidence band.
* Switch internal parametric regressions to use `sl3::Lrnr_glm` instead of
    `sl3::Lrnr_glm_fast`.

# txshift 0.2.7

As of July 2019:
* Improve argument names for clarity and update documentation.
* Addition of tighter unit tests for both one-step and TML estimators.

As of June 2019:
* Pin `sl3` dependency to version 1.2.0 of that package for stability.

# txshift 0.2.6

As of June 2019:
* Changes to arguments of `hal9001::fit_hal` for pseudo-outcome EIF regression.
* Addition of clarifying notes to core internal functions.
* Removal of outdated (and commented out) code in core internal functions.
* Clarifying alterations to internal function and argument names.
* Renaming internal function `tx_shift` to `shift_additive`.

# txshift 0.2.5

As of June 2019:
* Remove inverse weights from estimated efficient influence function necessary
    for pseudo-outcome regression for efficient IPCW-augmented estimators.
* Reduce use of redundant variables across core functions, reorganize functions
    across files, clarify documentation.
* Tweak arguments for fitting pseudo-outcome regression with HAL in order to
    diagnose performance issues revealed by simulation.
* Fix how inverse weights are passed to full-data estimators.
* Pare down arguments for the one-step estimation routine.

# txshift 0.2.4

As of June 2019:
* Introduce `bound_propensity` function to bound the propensity score away
    from zero by a factor 1/n, rather than to numerical precision.

As of April 2019:
* Minor improvements to documentation and vignettes.
* Fix a bug in the output of the IPCW one-step estimator.
* Pare down packages listed in imports, moving several to suggests.
* Introduce option to compute simultaneous confidence band for working MSMs.
* Fix a bug introduced by newly added imputation functionality in `sl3`.

# txshift 0.2.3

As of March 2019:
* Introduce functionality for computing one-step estimators to complement the
    the available TMLEs.
* Add initial functionality for summarizing estimated effects across a grid of
    shifts via working marginal structural models (MSMs).

# txshift 0.2.2

As of February 2019:
* Added helper functions and caught edge cases in auxiliary covariate for TMLE
    fluctuation models.
* Fixed a bug in how the auxiliary covariate for TMLEs is computed by keeping
    track of an extra shift g(a+2*delta|w).
* Revised inference machinery to create confidence intervals on the logit scale
    in the case of binary outcomes.

# txshift 0.2.0

As of May 2018:
* An initial public release of this package, version 0.2.0.
* This version including complete functionality for both standard TML and
    IPCW-TML estimators.
# Contributing to `txshift` development

We, the authors of the `txshift` R package, use the same guide as is used for
contributing to the development of the popular `tidyverse` ecosystem of R
packages. This document is simply a formal re-statement of that fact.

The goal of this guide is to help you get up and contributing to `txshift` as
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

To contribute a change to `txshift`, you follow these steps:

1. Create a branch in git and make your changes.
2. Push branch to GitHub and issue pull request (PR).
3. Discuss the pull request.
4. Iterate until either we accept the PR or decide that it's not a good fit for
   `txshift`.

Each of these steps are described in more detail below. This might feel
overwhelming the first time you get set up, but it gets easier with practice.

If you're not familiar with git or GitHub, please start by reading
<http://r-pkgs.had.co.nz/git.html>

Pull requests will be evaluated against a checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivates the
   need for change. Please describe the problem your PR addresses and show
   how your pull request solves it as concisely as possible.

   Also include this motivation in `NEWS` so that when a new release of
   `txshift` comes out it's easy for users to see what's changed. Add your
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

3.  __Use `txshift` coding style__. To do so, please follow the [official
    `tidyverse` style guide](http://style.tidyverse.org). Maintaining a
    consistent style across the whole code base makes it much easier to jump
    into the code. If you're modifying existing `txshift` code that doesn't
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
* There were 0 ERRORs.
* There were 0 WARNINGs.
* There were 0 NOTEs.

## Downstream dependencies
* There are no downstream dependencies.

## Additional notes
* This is a resubmission of a package removed from CRAN due to a Solaris build
  failure (no longer part of the CRAN build/check suite) that affected one of
  its upstream dependencies (`hal9001`).
* The upstream dependencies `hal9001` and `haldensify` have both recently been
  re-uploaded to CRAN.
---
title: "`txshift`: Efficient estimation of the causal effects of stochastic interventions in `R`"
tags:
  - causal inference
  - machine learning
  - two-phase sampling
  - efficient estimation
  - targeted learning
  - stochastic intervention
  - modified treatment policy
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1, 2
  - name: David Benkeser
    orcid: 0000-0002-1019-8343
    affiliation: 3
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Center for Computational Biology, University of California, Berkeley
    index: 2
  - name: Department of Biostatistics and Bioinformatics, Rollins School of Public Health, Emory University
    index: 3
date: 06 October 2020
bibliography: refs.bib
---

# Summary

Statistical causal inference has traditionally focused on effects defined by
inflexible static interventions, applicable only to binary or categorical
exposures. The evaluation of such interventions is often plagued by many
problems, both theoretical (e.g., non-identification) and practical (e.g.,
positivity violations); however, stochastic interventions provide a promising
solution to these fundamental issues [@diaz2018stochastic]. The `txshift` `R`
package provides researchers in (bio)statistics, epidemiology, health policy,
economics, and related disciplines with access to state-of-the-art statistical
methodology for evaluating the causal effects of stochastic shift interventions
on _continuous-valued_ exposures. `txshift` estimates the causal effects of
modified treatment policies (or "feasible interventions"), which take into
account the natural value of an exposure in assigning an intervention level. To
accommodate use in study designs incorporating outcome-dependent two-phase
sampling (e.g., case-control), the package provides two types of modern
corrections, both rooted in semiparametric theory, for constructing unbiased and
efficient estimates, despite the significant limitations induced by such
designs. Thus, `txshift` makes possible the estimation of the causal effects of
stochastic interventions in experimental and observational study settings
subject to real-world design limitations that commonly arise in modern
scientific practice.

# Statement of Need

Researchers seeking to build upon or apply cutting-edge statistical approaches
for causal inference often face significant obstacles: such methods are usually
not accompanied by robust, well-tested, and well-documented software packages.
Yet coding such methods from scratch is often impractical for the applied
researcher, as understanding the theoretical underpinnings of these methods
requires advanced training, severely complicating the assessment and testing of
bespoke causal inference software. What's more, even when such software tools
exist, they are usually minimal implementations, providing support only for
deploying the statistical method in problem settings untouched by the
complexities of real-world data. The `txshift` `R` package solves this problem
by providing an open source tool for evaluating the causal effects of flexible,
stochastic interventions, applicable to categorical or continuous-valued
exposures, while providing corrections for appropriately handling data generated
by commonly used but complex two-phase sampling designs.

# Background

Causal inference has traditionally focused on the effects of static
interventions, under which the magnitude of the exposure is set to a fixed,
prespecified value for each unit. The evaluation of such interventions faces
a host of issues, among them non-identification, violations of the assumption of
positivity, and inefficiency. Stochastic interventions provide a promising
solution to these fundamental issues by allowing for the target parameter to be
defined as the mean counterfactual outcome under a hypothetically shifted
version of the observed exposure distribution [@diaz2012population].
Modified treatment policies, a particular class of such interventions, may be
interpreted as shifting the natural exposure level at the level of a given
observational unit [@haneuse2013estimation; @diaz2018stochastic].

Despite the promise of such advances in causal inference, real data analyses are
often further complicated by economic constraints, such as when the primary
variable of interest is far more expensive to collect than auxiliary covariates.
Two-phase sampling is often used to bypass these limitations -- unfortunately,
these sampling schemes produce side effects that require further adjustment when
formal statistical inference is the principal goal of a study. Among the rich
literature on two-phase designs, @rose2011targeted2sd stand out for providing
a study of nonparametric efficiency theory under a broad class of two-phase
designs. Their work provides guidance on constructing efficient estimators of
causal effects under general two-phase sampling designs.

# `txshift`'s Scope

Building on these prior works, @hejazi2020efficient outlined a novel approach
for use in such settings: augmented targeted minimum loss (TML) and one-step
estimators for the causal effects of stochastic interventions, with guarantees
of consistency, efficiency, and multiple robustness despite the presence of
two-phase sampling. These authors further outlined a technique that summarizes
the effect of shifting an exposure variable on the outcome of interest via
a nonparametric working marginal structural model, analogous to a dose-response
analysis. The `txshift` software package, for the `R` language and environment
for statistical computing [@R], implements this methodology.

`txshift` is designed to facilitate the construction of TML and one-step
estimators of the causal effects of modified treatment policies that shift the
observed exposure value up (or down) by an arbitrary scalar $\delta$, which may
possibly take into account the natural value of the exposure (and, in future
versions, the covariates). The `R` package includes tools for deploying these
efficient estimators under outcome-dependent two-phase sampling designs, with
two types of corrections: (1) a reweighting procedure that introduces inverse
probability of censoring weights directly into relevant loss functions, as
discussed in @rose2011targeted2sd; as well as (2) an augmented efficient
influence function estimating equation, studied more thoroughly by
@hejazi2020efficient. `txshift` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) [@coyle2020sl3] to allow for ensemble
machine learning to be leveraged in the estimation of nuisance parameters.
What's more, the `txshift` package draws on both the `hal9001`
[@coyle2020hal9001; @hejazi2020hal9001] and `haldensify` [@hejazi2020haldensify]
`R` packages to allow each of the efficient estimators to be constructed in
a manner consistent with the methodological and theoretical advances of
@hejazi2020efficient, which require fast convergence rates of nuisance
parameters to their true counterparts for efficiency of the resultant estimator.

# Availability

The `txshift` package has been made publicly available both [via
GitHub](https://github.com/nhejazi/txshift) and the [Comprehensive `R` Archive
Network](https://CRAN.R-project.org/package=txshift). Use of the `txshift`
package has been extensively documented in the package's `README`, two
vignettes, and its [`pkgdown` documentation
website](https://code.nimahejazi.org/txshift).

# Acknowledgments

Nima Hejazi's contributions to this work were supported in part by a grant from
the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).

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


# R/`txshift`

<!-- badges: start -->
[![R-CMD-check](https://github.com/nhejazi/txshift/workflows/R-CMD-check/badge.svg)](https://github.com/nhejazi/txshift/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/nhejazi/txshift/master.svg)](https://codecov.io/github/nhejazi/txshift?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/txshift)](https://www.r-pkg.org/pkg/txshift)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/txshift)](https://CRAN.R-project.org/package=txshift)
[![CRAN total downloads](http://cranlogs.r-pkg.org/badges/grand-total/txshift)](https://CRAN.R-project.org/package=txshift)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447)
<!-- badges: end -->

> Efficient Estimation of the Causal Effects of Stochastic Interventions

__Authors:__ [Nima Hejazi](https://nimahejazi.org) and [David
Benkeser](https://www.sph.emory.edu/faculty/profile/#!dbenkes)

---

## What's `txshift`?

The `txshift` R package is designed to provide facilities for the construction
of efficient estimators of the counterfactual mean of an outcome under
stochastic interventions that depend on the natural value of treatment
[@diaz2012population; @haneuse2013estimation]. `txshift `implements and builds
upon a simplified algorithm for the targeted maximum likelihood (TML) estimator
of such a causal parameter, originally proposed by @diaz2018stochastic, and
makes use of analogous machinery to compute an efficient one-step estimator
[@pfanzagl1985contributions]. `txshift` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) [@coyle-sl3-rpkg] to allow for ensemble
machine learning to be leveraged in the estimation procedure.

For many practical applications (e.g., vaccine efficacy trials), observed data
is often subject to a two-phase sampling mechanism (i.e., through the use of a
two-stage design). In such cases, efficient estimators (of both varieties) must
be augmented to construct unbiased estimates of the population-level causal
parameter. @rose2011targeted2sd first introduced an augmentation procedure that
relies on introducing inverse probability of censoring (IPC) weights directly to
an appropriate loss function or to the efficient influence function estimating
equation. `txshift` extends this approach to compute IPC-weighted one-step and
TML estimators of the counterfactual mean outcome under a shift stochastic
treatment regime. The package is designed to implement the statistical
methodology described in @hejazi2020efficient and extensions thereof.

---

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=txshift) via

```{r cran-installation, eval = FALSE}
install.packages("txshift")
```

_Note:_ If `txshift` is installed from
[CRAN](https://CRAN.R-project.org/package=txshift), the `sl3`, an enhancing
dependency that allows ensemble machine learning to be used for nuisance
parameter estimation, won't be included. We highly recommend additionally
installing `sl3` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

```{r sl3-gh-master-installation, eval = FALSE}
remotes::install_github("tlverse/sl3@master")
```

For the latest features, install the most recent _stable version_  of `txshift`
from GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

```{r gh-master-installation, eval = FALSE}
remotes::install_github("nhejazi/txshift@master")
```

To contribute, install the _development version_ of `txshift` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

```{r gh-devel-installation, eval = FALSE}
remotes::install_github("nhejazi/txshift@devel")
```

---

## Example

To illustrate how `txshift` may be used to ascertain the effect of a treatment,
consider the following example:

```{r example, warning=FALSE}
library(txshift)
library(sl3)
set.seed(429153)

# simulate simple data
n_obs <- 500
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))

# now, let's introduce a a two-stage sampling process
C_samp <- rbinom(n_obs, 1, plogis(W + Y))

# fit the full-data TMLE (ignoring two-phase sampling)
tmle <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "tmle",
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle

# fit a full-data one-step estimator for comparison (again, no sampling)
os <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "onestep",
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
os

# fit an IPCW-TMLE to account for the two-phase sampling process
tmle_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = 0.5, C_samp = C_samp, V = c("W", "Y"),
  estimator = "tmle", max_iter = 5, eif_reg_type = "glm",
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle_ipcw

# compare with an IPCW-agumented one-step estimator under two-phase sampling
os_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = 0.5, C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep", eif_reg_type = "glm",
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
os_ipcw
```

---

## Issues

If you encounter any bugs or have any specific feature requests, please [file an
issue](https://github.com/nhejazi/txshift/issues). Further details on filing
issues are provided in our [contribution
guidelines](https://github.com/nhejazi/txshift/blob/master/CONTRIBUTING.md).

---

## Contributions

Contributions are very welcome. Interested contributors should consult our
[contribution
guidelines](https://github.com/nhejazi/txshift/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

---

## Citation

After using the `txshift` R package, please cite the following:

        @article{hejazi2020efficient,
          author = {Hejazi, Nima S and {van der Laan}, Mark J and Janes, Holly
            E and Gilbert, Peter B and Benkeser, David C},
          title = {Efficient nonparametric inference on the effects of
            stochastic interventions under two-phase sampling, with
            applications to vaccine efficacy trials},
          year = {2020},
          doi = {10.1111/biom.13375},
          url = {https://doi.org/10.1111/biom.13375},
          journal = {Biometrics},
          publisher = {Wiley Online Library}
        }

        @article{hejazi2020txshift-joss,
          author = {Hejazi, Nima S and Benkeser, David C},
          title = {{txshift}: Efficient estimation of the causal effects of
            stochastic interventions in {R}},
          year  = {2020},
          doi = {10.21105/joss.02447},
          url = {https://doi.org/10.21105/joss.02447},
          journal = {Journal of Open Source Software},
          publisher = {The Open Journal}
        }

        @software{hejazi2022txshift-rpkg,
          author = {Hejazi, Nima S and Benkeser, David C},
          title = {{txshift}: Efficient Estimation of the Causal Effects of
            Stochastic Interventions},
          year  = {2022},
          doi = {10.5281/zenodo.4070042},
          url = {https://CRAN.R-project.org/package=txshift},
          note = {R package version 0.3.7}
        }

---

## Related

* [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R package
  providing an independent implementation of the same core routines for the TML
  estimation procedure and statistical methodology as is made available here,
  through reliance on a unified interface for Targeted Learning provided by the
  [`tmle3`](https://github.com/tlverse/tmle3) engine of the [`tlverse`
  ecosystem](https://github.com/tlverse).

* [R/`medshift`](https://github.com/nhejazi/medshift) - An R package providing
  facilities to estimate the causal effect of stochastic treatment regimes in
  the mediation setting, including classical (IPW) and augmented double robust
  (one-step) estimators. This is an implementation of the methodology explored
  by @diaz2020causal.

* [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal package
  for estimating the conditional density treatment mechanism component of this
  parameter based on using the [highly adaptive
  lasso](https://github.com/tlverse/hal9001) [@coyle-hal9001-rpkg;
  @hejazi2020hal9001-joss] in combination with a pooled hazard regression. This
  package implements a variant of the approach advocated by @diaz2011super.

---

## Funding

The development of this software was supported in part through grants from the
National Library of Medicine (award no. [T32
LM012417](https://reporter.nih.gov/project-details/9248418)) and the National
Institute of Allergy and Infectious Diseases (award no.  [R01
AI074345](https://reporter.nih.gov/project-details/9926564)) of the National
Institutes of Health, as well as by the National Science Foundation (award no.
[DMS 2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

---

## License

&copy; 2017-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license. See below
for details:
```
MIT License

Copyright (c) 2017-2022 Nima S. Hejazi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## References

---
title: "Evaluating Causal Effects of Modified Treatment Policies"
author: "[Nima Hejazi](https://nimahejazi.org) and [David
  Benkeser](https://www.sph.emory.edu/faculty/profile/#!dbenkes)"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: ../inst/REFERENCES.bib
vignette: >
  %\VignetteIndexEntry{Evaluating Causal Effects of Modified Treatment Policies}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

Stochastic treatment regimes constitute a flexible framework for evaluating the
effects of continuous-valued exposures/treatments. _Modified treatment
policies_, one such technique within this framework, examine the effects
attributable to shifting the observed ("natural") value of a treatment, usually
up or down by some scalar $\delta$. The `txshift` package implements algorithms
for computing one-step or targeted minimum loss-based (TML) estimates of the
counterfactual means induced by additive modified treatment policies (MTPs),
defined by a shifting function $\delta(A,W)$. For a technical presentation, the
interested reader is invited to consult @diaz2018stochastic or the earlier work
of @diaz2012population and @haneuse2013estimation. For background on Targeted
Learning, consider consulting @vdl2011targeted, @vdl2018targeted, and
@vdl2022targeted.

To start, let's load the packages we'll need and set a seed for simulation:

```{r setup}
library(data.table)
library(haldensify)
library(txshift)
set.seed(11249)
```

---

## Data and Notation

We'll consider $n$ observed units $O_1, \ldots, O_n$, where each random variable
$O = (W, A, Y)$ corresponds to the data available on a single unit. Within $O$,
$W$ denotes baseline covariates (e.g., age, biological sex, BMI), $A \in
\mathbb{R}$ a continuous-valued exposure (e.g., dosage of nutritional
supplements taken), and $Y$ an outcome of interest (e.g., disease status).  To
minimize unjustifiable assumptions, we let $O \sim \mathcal{P} \in \mathcal{M}$,
where $\mathcal{P}$ is simply any distribution within the nonparametric
statistical model $\mathcal{M}$. To formalize the definition of stochastic
interventions and their corresponding causal effects, we consider a
nonparametric structural equation model (NPSEM), introduced by
@pearl2000causality, to define how the system changes under interventions of
interest:
\begin{align*}\label{eqn:npsem}
  W &= f_W(U_W) \\ A &= f_A(W, U_A) \\ Y &= f_Y(A, W, U_Y),
\end{align*}
We denote the observed data structure $O = (W, A, Y)$

Assuming that the distribution of $A$ conditional on $W = w$ has support in the
interval $(l(w), u(w))$ -- for convenience, we assume that the minimum natural
value of treatment $A$ for an individual with covariates $W = w$ is $l(w)$,
while, similarly, the maximum is $u(w)$ -- a simple MTP based on a shift
$\delta$, is
\begin{equation}\label{eqn:shift}
  \delta(a, w) =
  \begin{cases}
    a - \delta & \text{if } a > l(w) + \delta \\
    a & \text{if } a \leq l(w) + \delta,
  \end{cases}
\end{equation}
where $0 \leq \delta$ is an arbitrary pre-specified value that defines the
degree to which the observed value $A$ is to be shifted, where possible.

In case-cohort studies, it is common practice to make use of outcome-dependent
two-phase sampling designs, which allow for expensive measurements made on the
exposure (e.g., genomic sequencing of immune markers) to be avoided. As a
complication, such sampling schemes alter the observed data structure from the
simpler $O = (W, A, Y)$ to $O = (W, \Delta A, Y, \Delta)$, where the sampling
indicator $\Delta$ may itself be a function of the variables $\{W, Y\}$. In this
revised data structure, the value of $A$ is only observed for units in the
two-phase sample, for whom $\Delta = 1$. @hejazi2020efficient provide a detailed
investigation of the methodological details of efficient estimation under such
designs in the context of vaccine efficacy trials; their work was used in the
analysis of immune correlates of protection for HIV-1 [@hejazi2020efficient] and
COVID-19 [@gilbert2021covpn]. Of course, one may also account for loss to
follow-up (i.e., censoring), which the `txshift` package supports (through the
`C_cens` argument of the eponymous `txshift` function), though we avoid this
complication in our subsequent examples in the interest of clarity of
exposition. Corrections for both censoring and two-phase sampling make use of
inverse probability of censoring weighting (IPCW), leading to IPCW-augmented
one-step and TML estimators.

### Simulate Data

```{r make_data}
# parameters for simulation example
n_obs <- 200  # number of observations

# baseline covariate -- simple, binary
W <- rbinom(n_obs, 1, 0.5)

# create treatment based on baseline W
A <- rnorm(n_obs, mean = 2 * W, sd = 0.5)

# create outcome as a linear function of A, W + white noise
Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)

# two-phase sampling based on covariates
Delta_samp <- rbinom(n_obs, 1, plogis(W))

# treatment shift parameter
delta <- 0.5
```

## Estimating the Effects of Additive MTPs

The simplest way to compute an efficient estimator for an additive MTP is to fit
each of the nuisance parameters internally. This procedure can be sped up by
using generalized linear models (GLMs) to fit the outcome regression $Q_n$. The
`txshift()` function provides a simple interface.

```{r est_simple}
est_shift <- txshift(
  W = W, A = A, Y = Y, delta = delta,
  g_exp_fit_args = list(
    fit_type = "hal", n_bins = 5,
    grid_type = "equal_mass",
    lambda_seq = exp(seq(-1, -10, length = 100))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ .^2"
  )
)
est_shift
```

### Interlude: Ensemble Machine Learning with `sl3`

To easily incorporate ensemble machine learning into the estimation procedure,
`txshift` integrates with the [`sl3` R
package](https://tlverse.org/sl3/) [@coyle-sl3-rpkg]. For a complete guide on
using the `sl3` R package, consider consulting [the chapter on Super
Learning](https://tlverse.org/tlverse-handbook/sl3.html) in @vdl2022targeted.

```{r setup_sl, eval = FALSE}
library(sl3)

# SL learners to be used for most fits (e.g., IPCW, outcome regression)
mean_learner <- Lrnr_mean$new()
glm_learner <- Lrnr_glm$new()
rf_learner <- Lrnr_ranger$new()
Q_lib <- Stack$new(mean_learner, glm_learner, rf_learner)
sl_learner <- Lrnr_sl$new(learners = Q_lib, metalearner = Lrnr_nnls$new())

# SL learners for fitting the generalized propensity score fit
hose_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = glm_learner
)
hese_learner <- make_learner(Lrnr_density_semiparametric,
  mean_learner = rf_learner,
  var_learner = glm_learner
)
g_lib <- Stack$new(hose_learner, hese_learner)
sl_learner_density <- Lrnr_sl$new(
  learners = g_lib,
  metalearner = Lrnr_solnp_density$new()
)
```

### Efficient Effect Estimates with Machine Learning

Using the framework provided by the [`sl3`
package](https://github.com/tlverse/sl3/), the nuisance functions required for
our efficient estimators may be fit with ensemble machine learning. The Super
Learner algorithm [@vdl2007super] implemented in `sl3` uses the asymptotic
optimality of V-fold cross-validation [@dudoit2005asymptotics;
@vdl2004asymptotic; @vdv2006oracle] to select an optimal prediction functions
from a library or to construct an optimal combination of prediction functions.

```{r est_with_sl, eval = FALSE}
est_shift_sl <- txshift(
  W = W, A = A, Y = Y, delta = delta,
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = sl_learner_density
  ),
  Q_fit_args = list(
    fit_type = "sl",
    sl_learners = sl_learner
  )
)
est_shift_sl
```

## Estimating the Effects of Additive MTPs Under Two-Phase Sampling

In case-cohort studies in which two-phase sampling is used, the data structure
takes the from $O = (W, \Delta A, Y, \Delta)$ as previously discussed. Under
such sampling, the `txshift()` function may still be used to estimate the causal
effect of an additive MTP -- only a few additional arguments need to be
specified:

```{r ipcw_est_shift}
est_shift_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = delta,
  C_samp = Delta_samp, V = c("W", "Y"),
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "hal", n_bins = 5,
    grid_type = "equal_mass",
    lambda_seq = exp(seq(-1, -10, length = 100))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ .^2"
  ),
  eif_reg_type = "glm"
)
est_shift_ipcw
```

Note that we specify a few additional arguments in the call to `txshift()`,
including `C_samp`, the indicator of inclusion in the two-phase sample; `V`, the
set of other variables that may affect the sampling decision (in this case, both
the baseline covariates and the outcome); `samp_fit_args`, which indicates how
the sampling mechanism ought to be estimated; and `eif_reg_type`, which
indicates how a particular reduced-dimension nuisance regression ought to be
estimated (see @rose2011targeted2sd and @hejazi2020efficient for details). This
last argument only has options for using a GLM or the highly adaptive lasso
(HAL), a nonparametric regression estimator, using the [`hal9001`
package](https://github.com/tlverse/hal9001/) [@coyle-hal9001-rpkg;
hejazi2020hal9001-joss]. In practice, we recommend leaving this argument to the
default and fitting this nuisance function with HAL; however, this is
significantly more costly computationally.

## Statistical Inference for One-step and TML Estimators

The efficient estimators implemented in `txshift` are asymptotically linear;
thus, the estimator $\psi_n$ converges to the true parameter value $\psi_0$:
$$\psi_n - \psi_0 = (P_n - P_0) \cdot D(\bar{Q}_n^{\star}, g_n) +
  R(\hat{P}^{\star}, P_0),$$
which yields
$$\psi_n - \psi_0 = (P_n - P_0) \cdot D(P_0) + o_P \left( \frac{1}{\sqrt{n}}
 \right),$$
provided the following conditions,

1. if $D(\bar{Q}_n^{\star}, g_n)$ converges to $D(P_0)$ in $L_2(P_0)$ norm;
2. the size of the class of functions considered for estimation of
   $\bar{Q}_n^{\star}$ and $g_n$ is bounded (technically, $\exists \mathcal{F}$
   such that $D(\bar{Q}_n^{\star}, g_n) \in \mathcal{F}$ with high probability,
   where $\mathcal{F}$ is a Donsker class); and
3. the remainder term $R(\hat{P}^{\star}, P_0)$ decays as $o_P
   \left( \frac{1}{\sqrt{n}} \right)$.

By the central limit theorem, the estimators then have a Gaussian limiting
distribution,
$$\sqrt{n}(\psi_n - \psi) \to N(0, V(D(P_0))),$$
where $V(D(P_0))$ is the variance of the efficient influence function (or
canonical gradient).

The above implies that $\psi_n$ is a $\sqrt{n}$-consistent estimator of $\psi$,
that it is asymptotically normal (as given above), and that it is locally
efficient. This allows us to build Wald-type confidence intervals in a
straightforward manner:

$$\psi_n \pm z_{\alpha} \cdot \frac{\sigma_n}{\sqrt{n}},$$
where $\sigma_n^2$ is an estimator of $V(D(P_0))$. The estimator $\sigma_n^2$
may be obtained using the bootstrap or computed directly via the following

$$\sigma_n^2 = \frac{1}{n} \sum_{i = 1}^{n} D^2(\bar{Q}_n^{\star}, g_n)(O_i).$$

Such confidence intervals may easily be created with the `confint` method:

```{r confint}
(ci_est_shift <- confint(est_shift))
```

## _Advanced Usage:_ User-Specified Nuisance Regressions

In some special cases it may be useful for the experienced user to compute the
treatment mechanism, censoring mechanism, outcome regression, and sampling
mechanism fits separately (i.e., outside of the `txshift` wrapper function),
instead applying the wrapper only to construct an efficient one-step or TML
estimator. In such cases, the optional arguments ending in `_ext`.

```{r fit_external, eval=FALSE}
# compute censoring mechanism and produce IPC weights externally
pi_mech <- plogis(W)
ipcw_out <- pi_mech

# compute treatment mechanism (propensity score) externally
## first, produce the down-shifted treatment data
gn_downshift <- dnorm(A - delta, mean = tx_mult * W, sd = 1)
## next, initialize and produce the up-shifted treatment data
gn_upshift <- dnorm(A + delta, mean = tx_mult * W, sd = 1)
## now, initialize and produce the up-up-shifted (2 * delta) treatment data
gn_upupshift <- dnorm(A + 2 * delta, mean = tx_mult * W, sd = 1)
## then, initialize and produce the un-shifted treatment data
gn_noshift <- dnorm(A, mean = tx_mult * W, sd = 1)
## finally, put it all together into an object like what's produced internally
gn_out <- as.data.table(cbind(gn_downshift, gn_noshift, gn_upshift,
                              gn_upupshift))[C == 1, ]
setnames(gn_out, c("downshift", "noshift", "upshift", "upupshift"))

# compute outcome regression externally
# NOTE: transform Y to lie in the unit interval and bound predictions such that
#       no values fall near the bounds of the interval
Qn_noshift <- (W + A - min(Y)) / diff(range(Y))
Qn_upshift <- (W + A + delta - min(Y)) / diff(range(Y))
Qn_noshift[Qn_noshift < 0] <- 0.025
Qn_noshift[Qn_noshift > 1] <- 0.975
Qn_upshift[Qn_upshift < 0] <- 0.025
Qn_upshift[Qn_upshift > 1] <- 0.975
Qn_out <- as.data.table(cbind(Qn_noshift, Qn_upshift))[C == 1, ]
setnames(Qn_out, c("noshift", "upshift"))

# construct efficient estimator by applying wrapper function 
est_shift_spec <- txshift(
  W = W, A = A, Y = Y, delta = delta,
 samp_fit_args = NULL,
 samp_fit_ext = ipcw_out,
 g_exp_fit_args = list(fit_type = "external"),
 Q_fit_args = list(fit_type = "external"),
 gn_exp_fit_ext = gn_out,
 Qn_fit_ext = Qn_out
)
```

## References

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{est_Q}
\alias{est_Q}
\title{Estimate the Outcome Mechanism}
\usage{
est_Q(
  Y,
  C_cens = rep(1, length(Y)),
  A,
  W,
  delta = 0,
  samp_weights = rep(1, length(Y)),
  fit_type = c("sl", "glm"),
  glm_formula = "Y ~ .",
  sl_learners = NULL
)
}
\arguments{
\item{Y}{A \code{numeric} vector of observed outcomes.}

\item{C_cens}{A \code{numeric} vector of loss to follow-up indicators.}

\item{A}{A \code{numeric} vector of observed exposure values.}

\item{W}{A \code{numeric} matrix of observed baseline covariate values.}

\item{delta}{A \code{numeric} indicating the magnitude of the shift to be
computed for the exposure \code{A}. This is passed to the internal
\code{\link{shift_additive}} and is currently limited to additive shifts.}

\item{samp_weights}{A \code{numeric} vector of observation-level sampling
weights, as produced by the internal procedure to estimate the two-phase
sampling mechanism \code{\link{est_samp}}.}

\item{fit_type}{A \code{character} indicating whether to use GLMs or Super
Learner to fit the outcome regression. If the option "glm" is selected, the
argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
model formula (as per \code{\link[stats]{glm}}) as a \code{character}. If
the option "sl" is selected, the argument \code{sl_learners} must NOT be
\code{NULL}; instead, an instantiated \pkg{sl3} \code{Lrnr_sl} object,
specifying learners and a metalearner for the Super Learner fit, must be
provided. Consult the documentation of \pkg{sl3} for details.}

\item{glm_formula}{A \code{character} giving a \code{\link[stats]{formula}}
for fitting a (generalized) linear model via \code{\link[stats]{glm}}.}

\item{sl_learners}{Object containing a set of instantiated learners from the
\pkg{sl3}, to be used in fitting an ensemble model.}
}
\value{
A \code{data.table} with two columns, containing estimates of the
 outcome mechanism at the natural value of the exposure Q(A, W) and an
 upshift of the exposure Q(A + delta, W).
}
\description{
Estimate the Outcome Mechanism
}
\details{
Compute the outcome regression for the observed data, including
 with the shift imposed by the intervention. This returns the outcome
 regression for the observed data (at A) and under the counterfactual shift
 shift (at A + delta).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bound.R
\name{scale_to_unit}
\alias{scale_to_unit}
\title{Transform values by scaling to the unit interval}
\usage{
scale_to_unit(vals)
}
\arguments{
\item{vals}{A \code{numeric} vector corresponding to the observed values of
the variable of interest, to be re-scaled to the unit interval [0,1].}
}
\value{
A \code{numeric} vector of the same length as \code{vals}, where the
 values are re-scaled to lie in unit interval [0, 1].
}
\description{
Transform values by scaling to the unit interval
}
\details{
A transformation that scales an arbitrary set of input values to
 the unit interval. See \code{\link{scale_to_original}} for a corresponding
 backtransformation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{print.txshift_msm}
\alias{print.txshift_msm}
\title{Print Method for Marginal Structural Models}
\usage{
\method{print}{txshift_msm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{txshift_msm}.}

\item{...}{Other options (not currently used).}
}
\value{
None. Called for the side effect of printing an informative summary
 of slots of objects of class \code{txshift_msm}.
}
\description{
Print Method for Marginal Structural Models
}
\details{
The \code{print} method for objects of class \code{txshift_msm}.
}
\examples{
if (require("sl3")) {
  set.seed(3287)
  n_obs <- 1000
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(2 * A - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25)
  )
  print(msm)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bound.R
\name{bound_precision}
\alias{bound_precision}
\title{Bound Precision}
\usage{
bound_precision(vals)
}
\arguments{
\item{vals}{\code{numeric} vector of values in the interval [0, 1] to be
bounded within arbitrary machine precision. The most common use of this
functionality is to avoid indeterminate or non-finite values after the
application \code{stats::qlogis}.}
}
\value{
A \code{numeric} vector of the same length as \code{vals}, where
 the returned values are bounded to machine precision. This is intended to
 avoid numerical instability issues.
}
\description{
Bound Precision
}
\details{
Bound values in the unit interval to machine precision in order to
 avoid numerical instability issues in downstream computation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.txshift_msm}
\alias{plot.txshift_msm}
\title{Plot working MSM for causal effects of an intervention grid}
\usage{
\method{plot}{txshift_msm}(x, ...)
}
\arguments{
\item{x}{Object of class \code{txshift_msm} as produced by a call to
\code{\link{msm_vimshift}}.}

\item{...}{Additional arguments passed to \code{plot} as necessary.}
}
\description{
Plot working MSM for causal effects of an intervention grid
}
\details{
Creates a visualization of the intervention-specific counterfactual
 means as well as the working marginal structural model summarizing the
 trend across posited values of the intervention.
}
\examples{
if (require("sl3")) {
  set.seed(3287)
  n_obs <- 1000
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(2 * A - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25)
  )
  plot(msm)

  # fit a linear spline with knot at 0
  set.seed(8293)
  n_obs <- 1000
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(0.1 * A * (A >= 0) - 3 * A * (A < 0) - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25),
    msm_form = list(type = "piecewise", knot = 0)
  )
  plot(msm)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msm.R
\name{msm_vimshift}
\alias{msm_vimshift}
\title{Working marginal structural model for causal effects of an intervention grid}
\usage{
msm_vimshift(
  W,
  A,
  C_cens = rep(1, length(Y)),
  Y,
  C_samp = rep(1, length(Y)),
  V = NULL,
  delta_grid = seq(-0.5, 0.5, 0.5),
  msm_form = list(type = "linear", knot = NA),
  estimator = c("tmle", "onestep"),
  weighting = c("identity", "variance"),
  ci_level = 0.95,
  ci_type = c("marginal", "simultaneous"),
  ...
)
}
\arguments{
\item{W}{A \code{matrix}, \code{data.frame}, or similar containing a set of
baseline covariates.}

\item{A}{A \code{numeric} vector corresponding to a treatment variable. The
parameter of interest is defined as a location shift of this quantity.}

\item{C_cens}{A \code{numeric} indicator for whether a given observation was
subject to censoring by way of loss to follow-up. The default assumes no
censoring due to loss to follow-up.}

\item{Y}{A \code{numeric} vector of the observed outcomes.}

\item{C_samp}{A \code{numeric} indicator for whether a given observation was
subject to censoring by being omitted from the second-stage sample, used to
compute an inverse probability of censoring weighted estimator in such
cases. The default assumes no censoring due to two-phase sampling.}

\item{V}{The covariates that are used in determining the sampling procedure
that gives rise to censoring. The default is \code{NULL} and corresponds to
scenarios in which there is no censoring (in which case all values in the
preceding argument \code{C} must be uniquely 1. To specify this, pass in a
NAMED \code{list} identifying variables amongst W, A, Y that are thought to
have played a role in defining the sampling/censoring mechanism (C).}

\item{delta_grid}{A \code{numeric} vector giving the individual values of
the shift parameter used in computing each of the estimates.}

\item{msm_form}{A \code{list} specifying the type of working MSM to fit to
summarize the counterfactual means. The \code{list} has two components:
(1) \code{"type"}, which may be either "linear" or "piecewise", and (2)
\code{"knot"}, which, if specified, must be a value in \code{delta_grid}.
See examples for its use.}

\item{estimator}{The type of estimator to be fit, either \code{"tmle"} for
targeted maximum likelihood estimation or \code{"onestep"} for a one-step
augmented inverse probability weighted (AIPW) estimator.}

\item{weighting}{Whether to weight each parameter estimate by the inverse of
its variance (in order to improve stability of the resultant MSM fit) or to
simply weight all parameter estimates equally. The default is the option
\code{"identity"}, weighting all estimates identically.}

\item{ci_level}{A \code{numeric} indicating the desired coverage level of
the confidence interval to be computed.}

\item{ci_type}{Whether to construct a simultaneous confidence band covering
all parameter estimates at once or marginal confidence intervals covering
each parameter estimate separately. The default is to construct marginal
confidence intervals for each parameter estimate rather than a simultaneous
confidence band.}

\item{...}{Additional arguments to be passed to \code{\link{txshift}}.}
}
\value{
A \code{list} containing estimates of the individual counterfactual
 means over a grid in the shift parameters (\code{delta_grid}), alongside
 the estimate of a marginal structural model that summarizes a trend through
 these counterfactual means.
}
\description{
Working marginal structural model for causal effects of an intervention grid
}
\details{
Computes estimates of the counterfactual mean over a grid of shift
 stochastic interventions and fits a working marginal structural model to
 summarize the trend through the counterfactual means as a function of the
 specified shift intervention. The working marginal structural model may be
 linear in the shift parameter or piecewise linear with a single knot point.
 Provides support for two weighting schemes, may be used with either of the
 one-step or TML estimators, and also allows the construction of marginal or
 simultaneous confidence intervals.
}
\examples{
if (require("sl3")) {
  n_obs <- 100
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(2 * A - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25)
  )

  # fit a linear spline with knot at 0
  n_obs <- 100
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(0.1 * A * (A >= 0) - 3 * A * (A < 0) - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25),
    msm_form = list(type = "piecewise", knot = 0)
  )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{print.txshift}
\alias{print.txshift}
\title{Print Method for Counterfactual Mean of Stochastic Shift Intervention}
\usage{
\method{print}{txshift}(x, ..., ci_level = 0.95)
}
\arguments{
\item{x}{An object of class \code{txshift}.}

\item{...}{Other options (not currently used).}

\item{ci_level}{A \code{numeric} indicating the level of the confidence
interval to be computed.}
}
\value{
None. Called for the side effect of printing an informative summary
 of slots of objects of class \code{txshift}.
}
\description{
Print Method for Counterfactual Mean of Stochastic Shift Intervention
}
\details{
The \code{print} method for objects of class \code{txshift}.
}
\examples{
if (require("sl3")) {
  set.seed(429153)
  n_obs <- 100
  W <- replicate(2, rbinom(n_obs, 1, 0.5))
  A <- rnorm(n_obs, mean = 2 * W, sd = 1)
  Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
  txout <- txshift(
    W = W, A = A, Y = Y, delta = 0.5,
    estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    )
  )
  print(txout)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{est_g_exp}
\alias{est_g_exp}
\title{Estimate the Exposure Mechanism via Generalized Propensity Score}
\usage{
est_g_exp(
  A,
  W,
  delta = 0,
  samp_weights = rep(1, length(A)),
  fit_type = c("hal", "sl"),
  sl_learners_density = NULL,
  haldensify_args = list(grid_type = "equal_range", lambda_seq = exp(seq(-1, -13,
    length = 300)))
)
}
\arguments{
\item{A}{A \code{numeric} vector of observed exposure values.}

\item{W}{A \code{numeric} matrix of observed baseline covariate values.}

\item{delta}{A \code{numeric} value identifying a shift in the observed
value of the exposure under which observations are to be evaluated.}

\item{samp_weights}{A \code{numeric} vector of observation-level sampling
weights, as produced by the internal procedure to estimate the two-phase
sampling mechanism \code{\link{est_samp}}.}

\item{fit_type}{A \code{character} specifying whether to use Super Learner
(from \pkg{sl3}) or the Highly Adaptive Lasso (from \pkg{hal9001}) to
estimate the conditional exposure density.}

\item{sl_learners_density}{Object containing a set of instantiated learners
from \pkg{sl3}, to be used in fitting an ensemble model.}

\item{haldensify_args}{A \code{list} of argument to be directly passed to
\code{\link[haldensify]{haldensify}} when \code{fit_type} is set to
\code{"hal"}. Note that this invokes the Highly Adaptive Lasso instead of
Super Learner and is thus only feasible for relatively small data sets.}
}
\value{
A \code{data.table} with four columns, containing estimates of the
 generalized propensity score at a downshift (g(A - delta | W)), no shift
 (g(A | W)), an upshift (g(A + delta) | W), and an upshift of magnitude two
 (g(A + 2 delta) | W).
}
\description{
Estimate the Exposure Mechanism via Generalized Propensity Score
}
\details{
Compute the propensity score (exposure mechanism) for the observed
 data, including the shift. This gives the propensity score for the observed
 data (at the observed A) the counterfactual shifted exposure levels (at
 {A - delta}, {A + delta}, and {A + 2 * delta}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eifs.R
\name{eif}
\alias{eif}
\title{Compute the Shift Parameter Estimate and the Efficient Influence Function}
\usage{
eif(
  Y,
  Qn,
  Hn,
  estimator = c("tmle", "onestep"),
  fluc_mod_out = NULL,
  C_samp = rep(1, length(Y)),
  ipc_weights = rep(1, length(Y))
)
}
\arguments{
\item{Y}{A \code{numeric} vector of the observed outcomes.}

\item{Qn}{An object providing the value of the outcome evaluated after
imposing a shift in the treatment. This object is passed in after being
constructed by a call to the internal function \code{est_Q}.}

\item{Hn}{An object providing values of the auxiliary ("clever") covariate,
constructed from the treatment mechanism and required for targeted minimum
loss-based estimation. This object object should be passed in after being
constructed by a call to the internal function \code{est_Hn}.}

\item{estimator}{The type of estimator to be fit, either \code{"tmle"} for
targeted maximum likelihood estimation or \code{"onestep"} for a one-step
estimator.}

\item{fluc_mod_out}{An object giving values of the logistic tilting model
for targeted minimum loss estimation. This type of object should be the
output of the internal routines to perform this step of the TML estimation
procedure, as given by \code{\link{fit_fluctuation}}.}

\item{C_samp}{Indicator for missingness due to exclusion from second-phase
sample. Used for compatibility with the IPCW-TML estimation routine.}

\item{ipc_weights}{A \code{numeric} vector that gives inverse probability of
censoring weights for each observation. These are generated by invoking the
routines for estimating the censoring mechanism.}
}
\value{
A \code{list} containing the parameter estimate, estimated variance
 based on the efficient influence function (EIF), the estimate of the EIF
 incorporating inverse probability of censoring weights, and the estimate of
 the EIF without the application of such weights.
}
\description{
Compute the Shift Parameter Estimate and the Efficient Influence Function
}
\details{
Estimate the value of the causal parameter alongside statistical
 inference for the parameter estimate based on the efficient influence
 function of the target parameter, which takes the following form:
 %D(P)(o) = H(a,w)(y - \bar{Q}(a,w)) + \bar{Q}(d(a,w)) - \psi(P)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{est_Hn}
\alias{est_Hn}
\title{Estimate Auxiliary Covariate of Full Data Efficient Influence Function}
\usage{
est_Hn(gn_exp)
}
\arguments{
\item{gn_exp}{An estimate of the exposure density (a generalized propensity
score) using the output provided by \code{\link{est_g_exp}}.}
}
\value{
A \code{data.table} with two columns, containing estimates of the
 auxiliary covariate at the natural value of the exposure H(A, W) and at the
 shifted value of the exposure H(A + delta, W).
}
\description{
Estimate Auxiliary Covariate of Full Data Efficient Influence Function
}
\details{
Compute an estimate of the auxiliary covariate of the efficient
 influence function required to update initial estimates through logistic
 tilting models for targeted minimum loss estimation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmle_txshift.R
\name{fit_fluctuation}
\alias{fit_fluctuation}
\title{Fit One-Dimensional Fluctuation Model for Updating Initial Estimates}
\usage{
fit_fluctuation(
  Y,
  Qn_scaled,
  Hn,
  ipc_weights = rep(1, length(Y)),
  method = c("standard", "weighted"),
  flucmod_tol = 50
)
}
\arguments{
\item{Y}{A \code{numeric} vector corresponding to an outcome variable.}

\item{Qn_scaled}{An object providing the value of the outcome evaluate
after inducing a shift in the exposure. This object should be passed in
after being constructed by a call to \code{\link{est_Q}}.}

\item{Hn}{An object providing values of the auxiliary ("clever") covariate,
constructed from the treatment mechanism and required for targeted minimum
loss estimation. This object object should be passed in after being
constructed by a call to \code{\link{est_Hn}}.}

\item{ipc_weights}{A \code{numeric} vector that gives inverse probability of
censoring weights for each observation. These are generated by invoking the
routines for estimating the censoring mechanism.}

\item{method}{A \code{character} giving the type of regression to be used in
traversing the fluctuation sub-model. The available choices are "weighted"
and "standard". Consult the literature for details on the differences.}

\item{flucmod_tol}{A \code{numeric} indicating the largest value to be
tolerated in the fluctuation model for the targeted minimum loss estimator.}
}
\value{
A \code{list} containing the fluctuation model (a \code{glm} object)
 produced by logistic regression, a \code{character} vector indicating the
 type of fluctuation (whether the auxiliary covariates was used as a weight
 or included directly in the model formula), the updated estimates of the
 outcome regression under the shifted value of the exposure, and the updated
 estimates of the outcome regression under the natural value of exposure.
}
\description{
Fit One-Dimensional Fluctuation Model for Updating Initial Estimates
}
\details{
Procedure for fitting a one-dimensional fluctuation model to update
 the initial estimates of the outcome regression based on the auxiliary
 covariate. These updated estimates are subsequently used to construct the
 TML estimator of the counterfactual mean under a modified treatment policy.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eifs.R
\name{ipcw_eif_update}
\alias{ipcw_eif_update}
\title{Iterative IPCW Update Procedure of Augmented Efficient Influence Function}
\usage{
ipcw_eif_update(
  data_internal,
  C_samp,
  V,
  ipc_mech,
  ipc_weights,
  Qn_estim,
  Hn_estim,
  estimator = c("tmle", "onestep"),
  fluctuation = NULL,
  flucmod_tol = 50,
  eif_reg_type = c("hal", "glm")
)
}
\arguments{
\item{data_internal}{A \code{data.table} containing of the observations
selected into the second-phase sample.}

\item{C_samp}{A \code{numeric} indicator for missingness due to exclusion
the from second-stage sample.}

\item{V}{A \code{data.table} giving the values across all observations of
all variables that play a role in the censoring mechanism.}

\item{ipc_mech}{A \code{numeric} vector of the censoring mechanism estimates
all of the observations, only for the two-phase sampling mechanism. Note
well that these values do NOT account for censoring from loss to follow-up.}

\item{ipc_weights}{A \code{numeric} vector of inverse probability of
censoring weights, including such weights for censoring due to loss to
follow-up. Without loss to follow-up, these are equivalent to \code{C_samp
/ ipc_mech} in an initial run of this procedure.}

\item{Qn_estim}{A \code{data.table} corresponding to the outcome regression.
This is produced by invoking the internal function \code{est_Q}.}

\item{Hn_estim}{A \code{data.table} corresponding to values produced in the
computation of the auxiliary ("clever") covariate. This is produced easily
by invoking the internal function \code{est_Hn}.}

\item{estimator}{The type of estimator to be fit, either \code{"tmle"} for
targeted maximum likelihood estimation or \code{"onestep"} for a one-step
estimator.}

\item{fluctuation}{A \code{character} giving the type of regression to be
used in traversing the fluctuation submodel. The choices are "weighted" and
"standard".}

\item{flucmod_tol}{A \code{numeric} indicating the largest value to be
tolerated in the fluctuation model for the targeted minimum loss estimator.}

\item{eif_reg_type}{Whether a flexible nonparametric function ought to be
used in the dimension-reduced nuisance regression of the targeting step for
the censored data case. By default, the method used is a nonparametric
regression based on the Highly Adaptive Lasso (from \pkg{hal9001}). Set
this to \code{"glm"} to instead use a simple linear regression model. In
this step, the efficient influence function (EIF) is regressed against
covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).}
}
\value{
A \code{list} containing the estimated outcome mechanism, the fitted
 fluctuation model for TML updates, the updated inverse probability of
 censoring weights (IPCW), the updated estimate of the efficient influence
 function, and the estimated IPCW component of the EIF.
}
\description{
Iterative IPCW Update Procedure of Augmented Efficient Influence Function
}
\details{
An adaptation of the IPCW-TMLE for iteratively constructing an
 efficient inverse probability of censoring weighted TML or one-step
 estimator. The efficient influence function of the parameter and updating
 the IPC weights in an iterative process, until a convergence criteria is
 satisfied.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bound.R
\name{scale_to_original}
\alias{scale_to_original}
\title{Transform values from the unit interval back to their original scale}
\usage{
scale_to_original(scaled_vals, max_orig, min_orig)
}
\arguments{
\item{scaled_vals}{A \code{numeric} vector corresponding to re-scaled values
in the unit interval, to be re-scaled to the original interval.}

\item{max_orig}{A \code{numeric} scalar value giving the maximum of the
values on the original scale.}

\item{min_orig}{A \code{numeric} scalar value giving the minimum of the
values on the original scale.}
}
\value{
A \code{numeric} vector of the same length as \code{scaled_vals},
 where the values are re-scaled to lie in their original/natural interval.
}
\description{
Transform values from the unit interval back to their original scale
}
\details{
A back-transformation that returns values computed in the unit
 interval to their original scale. This is used in re-scaling updated TML
 estimates back to their natural scale. Undoes \code{\link{scale_to_unit}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onestep_txshift.R
\name{onestep_txshift}
\alias{onestep_txshift}
\title{One-Step Estimate of Counterfactual Mean of Stochastic Shift Intervention}
\usage{
onestep_txshift(
  data_internal,
  C_samp = rep(1, nrow(data_internal)),
  V = NULL,
  delta,
  samp_estim,
  gn_cens_weights,
  Qn_estim,
  Hn_estim,
  eif_reg_type = c("hal", "glm"),
  samp_fit_args,
  ipcw_efficiency = TRUE
)
}
\arguments{
\item{data_internal}{A \code{data.table} constructed internally by a call to
\code{\link{txshift}}. This contains the data elements needed for computing
the one-step estimator.}

\item{C_samp}{A \code{numeric} indicator for whether a given observation was
included in the second-stage sample, used to compute an IPC-weighted
one-step estimator in cases where two-stage sampling is performed. Default
assumes no censoring due to sampling.}

\item{V}{The covariates that are used in determining the sampling procedure
that gives rise to censoring. The default is \code{NULL} and corresponds to
scenarios in which there is no censoring (in which case all values in the
preceding argument \code{C} must be uniquely 1. To specify this, pass in a
NAMED \code{list} identifying variables amongst W, A, Y that are thought to
have played a role in defining the sampling/censoring mechanism (C).}

\item{delta}{A \code{numeric} value indicating the shift in the treatment to
be used in defining the target parameter. This is defined with respect to
the scale of the treatment (A).}

\item{samp_estim}{An object providing the value of the censoring mechanism
evaluated across the full data. This object is passed in after being
constructed by a call to the internal function \code{\link{est_samp}}.}

\item{gn_cens_weights}{TODO: document}

\item{Qn_estim}{An object providing the value of the outcome evaluated after
imposing a shift in the treatment. This object is passed in after being
constructed by a call to the internal function \code{est_Q}.}

\item{Hn_estim}{An object providing values of the auxiliary ("clever")
covariate, constructed from the treatment mechanism and required for
targeted minimum loss estimation. This object object should be passed in
after being constructed by a call to the internal function \code{est_Hn}.}

\item{eif_reg_type}{Whether a flexible nonparametric function ought to be
used in the dimension-reduced nuisance regression of the targeting step for
the censored data case. By default, the method used is a nonparametric
regression based on the Highly Adaptive Lasso (from \pkg{hal9001}). Set
this to \code{"glm"} to instead use a simple linear regression model.
In this step, the efficient influence function (EIF) is regressed against
covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).}

\item{samp_fit_args}{A \code{list} of arguments, all but one of which are
passed to \code{\link{est_samp}}. For details, consult the documentation
for \code{\link{est_samp}}. The first element (i.e., \code{fit_type}) is
used to determine how this regression is fit: "glm" for generalized linear
model, "sl" for a Super Learner, and "external" for a user-specified input
of the form produced by \code{\link{est_samp}}.}

\item{ipcw_efficiency}{Whether to invoke an augmentation of the IPCW-TMLE
procedure that performs an iterative process to ensure efficiency of the
resulting estimate. The default is \code{TRUE}; set to \code{FALSE} to use
an IPC-weighted loss rather than the IPC-augmented influence function.}
}
\value{
S3 object of class \code{txshift} containing the results of the
 procedure to compute a one-step estimate of the treatment shift parameter.
}
\description{
One-Step Estimate of Counterfactual Mean of Stochastic Shift Intervention
}
\details{
Invokes the procedure to construct a one-step estimate of the
 counterfactual mean under a modified treatment policy.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bound.R
\name{bound_propensity}
\alias{bound_propensity}
\title{Bound Generalized Propensity Score}
\usage{
bound_propensity(vals)
}
\arguments{
\item{vals}{\code{numeric} vector of propensity score estimate values. Note
that, for this parameter, the propensity score is (conditional) density and
so it ought not be bounded from above.}
}
\value{
A \code{numeric} vector of the same length as \code{vals}, where the
 returned values are bounded such that the minimum is no lower than 1/n, for
 the sample size n.
}
\description{
Bound Generalized Propensity Score
}
\details{
Bound estimated values of the generalized propensity score (a
 conditional density) to avoid numerical instability issues arising from
 practical violations of the assumption of positivity.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{est_samp}
\alias{est_samp}
\title{Estimate Probability of Censoring by Two-Phase Sampling}
\usage{
est_samp(V, C_samp, fit_type = c("sl", "glm"), sl_learners = NULL)
}
\arguments{
\item{V}{A \code{numeric} vector, \code{matrix}, \code{data.frame} or
similar object giving the observed values of the covariates known to
potentially inform the sampling mechanism.}

\item{C_samp}{A \code{numeric} vector of observed values of the indicator
for inclusion in the second-phase sample.}

\item{fit_type}{A \code{character} indicating whether to perform the fit
using GLMs or a Super Learner ensemble model. If use of Super Learner is
desired, then the argument \code{sl_learners} must be provided.}

\item{sl_learners}{An \pkg{sl3} \code{Lrnr_sl} object, a Super Learner
ensemble or learner instantiated externally using \pkg{sl3}.}
}
\value{
A \code{numeric} vector of the estimated sampling mechanism.
}
\description{
Estimate Probability of Censoring by Two-Phase Sampling
}
\details{
Compute estimates of the sampling probability for inclusion in the
 the second-phase via the two-phase sampling mechanism. These estimates are
 used for the creation of inverse probability weights.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift_funs.R
\name{shift_additive}
\alias{shift_additive}
\title{Simple Additive Modified Treatment Policy}
\usage{
shift_additive(A, W = NULL, delta)
}
\arguments{
\item{A}{A \code{numeric} vector of observed treatment values.}

\item{W}{A \code{numeric} matrix of observed baseline covariate values.}

\item{delta}{A \code{numeric} indicating the magnitude of the shift to be
computed for the treatment \code{A}.}
}
\value{
A \code{numeric} vector containing the shifted exposure values.
}
\description{
Simple Additive Modified Treatment Policy
}
\details{
A simple modified treatment policy that modifes the observed value
 of the exposure by shifting it by a value \code{delta}. Note that this
 shifting function assumes support of A|W across all strata of W.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmle_txshift.R
\name{tmle_txshift}
\alias{tmle_txshift}
\title{Targeted Minimum Loss Estimate of Counterfactual Mean of Stochastic Shift
Intervention}
\usage{
tmle_txshift(
  data_internal,
  C_samp = rep(1, nrow(data_internal)),
  V = NULL,
  delta,
  samp_estim,
  gn_cens_weights,
  Qn_estim,
  Hn_estim,
  fluctuation = c("standard", "weighted"),
  max_iter = 10,
  eif_reg_type = c("hal", "glm"),
  samp_fit_args,
  ipcw_efficiency = TRUE
)
}
\arguments{
\item{data_internal}{A \code{data.table} constructed internally by a call to
\code{\link{txshift}}. This contains most of the data for computing the
targeted minimum loss (TML) estimator.}

\item{C_samp}{A \code{numeric} indicator for whether a given observation was
included in the second-stage sample, used to compute an IPC-weighted
one-step estimator in cases where two-stage sampling is performed. Default
assumes no censoring due to sampling.}

\item{V}{The covariates that are used in determining the sampling procedure
that gives rise to censoring. The default is \code{NULL} and corresponds to
scenarios in which there is no censoring (in which case all values in the
preceding argument \code{C_samp} must be 1. To specify this, pass in a
NAMED \code{list} identifying variables amongst W, A, Y that are thought to
have played a role in defining the sampling mechanism.}

\item{delta}{A \code{numeric} value indicating the shift in the treatment to
be used in defining the target parameter. This is defined with respect to
the scale of the treatment (A).}

\item{samp_estim}{An object providing the value of the sampling mechanism
evaluated across the full data. This object is passed in after being
constructed by a call to the internal function \code{\link{est_samp}}.}

\item{gn_cens_weights}{TODO: document}

\item{Qn_estim}{An object providing the value of the outcome evaluated after
imposing a shift in the treatment. This object is passed in after being
constructed by a call to the internal function \code{\link{est_Q}}.}

\item{Hn_estim}{An object providing values of the auxiliary ("clever")
covariate, constructed from the treatment mechanism and required for
targeted minimum loss-based estimation. This object object should be passed
in after being constructed by a call to \code{\link{est_Hn}}.}

\item{fluctuation}{The method to be used in the submodel fluctuation step
(targeting step) to compute the TML estimator. The choices are "standard"
and "weighted" for where to place the auxiliary covariate in the logistic
tilting regression.}

\item{max_iter}{A \code{numeric} integer giving the maximum number of steps
to be taken in iterating to a solution of the efficient influence function.}

\item{eif_reg_type}{Whether a flexible nonparametric function ought to be
used in the dimension-reduced nuisance regression of the targeting step for
the censored data case. By default, the method used is a nonparametric
regression based on the Highly Adaptive Lasso (from \pkg{hal9001}).
Set this to \code{"glm"} to instead use a simple linear regression model.
In this step, the efficient influence function (EIF) is regressed against
covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).}

\item{samp_fit_args}{A \code{list} of arguments, all but one of which are
passed to \code{\link{est_samp}}. For details, consult the documentation
for \code{\link{est_samp}}. The first element (i.e., \code{fit_type}) is
used to determine how this regression is fit: "glm" for generalized linear
model, "sl" for a Super Learner, and "external" for a user-specified input
of the form produced by \code{\link{est_samp}}.}

\item{ipcw_efficiency}{Whether to invoke an augmentation of the IPCW-TMLE
procedure that performs an iterative process to ensure efficiency of the
resulting estimate. The default is \code{TRUE}; set to \code{FALSE} to use
an IPC-weighted loss rather than the IPC-augmented influence function.}
}
\value{
S3 object of class \code{txshift} containing the results of the
 procedure to compute a TML estimate of the treatment shift parameter.
}
\description{
Targeted Minimum Loss Estimate of Counterfactual Mean of Stochastic Shift
Intervention
}
\details{
Invokes the procedure to construct a targeted minimum loss estimate
 (TMLE) of the counterfactual mean under a modified treatment policy.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{est_g_cens}
\alias{est_g_cens}
\title{Estimate the Censoring Mechanism}
\usage{
est_g_cens(
  C_cens,
  A,
  W,
  samp_weights = rep(1, length(C_cens)),
  fit_type = c("sl", "glm"),
  glm_formula = "C_cens ~ .",
  sl_learners = NULL
)
}
\arguments{
\item{C_cens}{A \code{numeric} vector of loss to follow-up indicators.}

\item{A}{A \code{numeric} vector of observed exposure values.}

\item{W}{A \code{numeric} matrix of observed baseline covariate values.}

\item{samp_weights}{A \code{numeric} vector of observation-level sampling
weights, as produced by the internal procedure to estimate the two-phase
sampling mechanism \code{\link{est_samp}}.}

\item{fit_type}{A \code{character} indicating whether to use GLMs or Super
Learner to fit the censoring mechanism. If option "glm" is selected, the
argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
model formula (as per \code{\link[stats]{glm}}) as a \code{character}. If
the option "sl" is selected, the argument \code{sl_learners} must NOT be
\code{NULL}; instead, an instantiated \pkg{sl3} \code{Lrnr_sl} object,
specifying learners and a metalearner for the Super Learner fit, must be
provided. Consult the documentation of \pkg{sl3} for details.}

\item{glm_formula}{A \code{character} giving a \code{\link[stats]{formula}}
for fitting a (generalized) linear model via \code{\link[stats]{glm}}.}

\item{sl_learners}{Object containing a set of instantiated learners from the
\pkg{sl3}, to be used in fitting an ensemble model.}
}
\value{
A \code{numeric} vector of the propensity score for censoring.
}
\description{
Estimate the Censoring Mechanism
}
\details{
Compute the censoring mechanism for the observed data, in order to
 apply a joint intervention for removing censoring by re-weighting.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint.R
\name{confint.txshift}
\alias{confint.txshift}
\title{Confidence Intervals for Counterfactual Mean Under Stochastic Intervention}
\usage{
\method{confint}{txshift}(object, parm = seq_len(object$psi), level = 0.95, ..., ci_mult = NULL)
}
\arguments{
\item{object}{An object of class \code{txshift}, as produced by invoking
the function \code{\link{txshift}}, for which a confidence interval is to
be computed.}

\item{parm}{A \code{numeric} vector indicating indices of \code{object$est}
for which to return confidence intervals.}

\item{level}{A \code{numeric} indicating the level of the confidence
interval to be computed.}

\item{...}{Other arguments. Not currently used.}

\item{ci_mult}{Pre-computed multipliers for generating confidence intervals.
The default of \code{NULL} should generally NOT be changed and is only used
by the internal machinery for creating simultaneous confidence bands.}
}
\value{
A named \code{numeric} vector containing the parameter estimate from
 a \code{txshift} object, alongside lower and upper Wald-style confidence
 intervals at a specified coverage level.
}
\description{
Confidence Intervals for Counterfactual Mean Under Stochastic Intervention
}
\details{
Compute confidence intervals for estimates produced by
 \code{\link{txshift}}.
}
\examples{
set.seed(429153)
n_obs <- 100
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
txout <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "tmle",
  g_exp_fit_args = list(
    fit_type = "hal", n_bins = 5,
    grid_type = "equal_mass",
    lambda_seq = exp(-1:-9)
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  )
)
confint(txout)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/txshift.R
\name{txshift}
\alias{txshift}
\title{Efficient Estimate of Counterfactual Mean of Stochastic Shift Intervention}
\usage{
txshift(
  W,
  A,
  C_cens = rep(1, length(A)),
  Y,
  C_samp = rep(1, length(Y)),
  V = NULL,
  delta = 0,
  estimator = c("tmle", "onestep"),
  fluctuation = c("standard", "weighted"),
  max_iter = 10,
  samp_fit_args = list(fit_type = c("glm", "sl", "external"), sl_learners = NULL),
  g_exp_fit_args = list(fit_type = c("hal", "sl", "external"), lambda_seq = exp(seq(-1,
    -13, length = 300)), sl_learners_density = NULL),
  g_cens_fit_args = list(fit_type = c("glm", "sl", "external"), glm_formula =
    "C_cens ~ .^2", sl_learners = NULL),
  Q_fit_args = list(fit_type = c("glm", "sl", "external"), glm_formula = "Y ~ .^2",
    sl_learners = NULL),
  eif_reg_type = c("hal", "glm"),
  ipcw_efficiency = TRUE,
  samp_fit_ext = NULL,
  gn_exp_fit_ext = NULL,
  gn_cens_fit_ext = NULL,
  Qn_fit_ext = NULL
)
}
\arguments{
\item{W}{A \code{matrix}, \code{data.frame}, or similar containing a set of
baseline covariates.}

\item{A}{A \code{numeric} vector corresponding to a treatment variable. The
parameter of interest is defined as a location shift of this quantity.}

\item{C_cens}{A \code{numeric} indicator for whether a given observation was
subject to censoring by way of loss to follow-up. The default assumes no
censoring due to loss to follow-up.}

\item{Y}{A \code{numeric} vector of the observed outcomes.}

\item{C_samp}{A \code{numeric} indicator for whether a given observation was
subject to censoring by being omitted from the second-stage sample, used to
compute an inverse probability of censoring weighted estimator in such
cases. The default assumes no censoring due to two-phase sampling.}

\item{V}{The covariates that are used in determining the sampling procedure
that gives rise to censoring. The default is \code{NULL} and corresponds to
scenarios in which there is no censoring (in which case all values in the
preceding argument \code{C_samp} must be uniquely 1). To specify this, pass
in a \code{character} vector identifying variables amongst W, A, Y thought
to have impacted the definition of the sampling mechanism (C_samp). This
argument also accepts a \code{data.table} (or similar) object composed of
combinations of variables W, A, Y; use of this option is NOT recommended.}

\item{delta}{A \code{numeric} value indicating the shift in the treatment to
be used in defining the target parameter. This is defined with respect to
the scale of the treatment (A).}

\item{estimator}{The type of estimator to be fit, either \code{"tmle"} for
targeted maximum likelihood or \code{"onestep"} for a one-step estimator.}

\item{fluctuation}{The method to be used in the submodel fluctuation step
(targeting step) to compute the TML estimator. The choices are "standard"
and "weighted" for where to place the auxiliary covariate in the logistic
tilting regression.}

\item{max_iter}{A \code{numeric} integer giving the maximum number of steps
to be taken in iterating to a solution of the efficient influence function.}

\item{samp_fit_args}{A \code{list} of arguments, all but one of which are
passed to \code{\link{est_samp}}. For details, consult the documentation of
\code{\link{est_samp}}. The first element (i.e., \code{fit_type}) is used
to determine how this regression is fit: generalized linear model ("glm")
or Super Learner ("sl"), and "external" a user-specified input of the form
produced by \code{\link{est_samp}}.}

\item{g_exp_fit_args}{A \code{list} of arguments, all but one of which are
passed to \code{\link{est_g_exp}}. For details, see the documentation of
\code{\link{est_g_exp}}. The 1st element (i.e., \code{fit_type}) specifies
how this regression is fit: \code{"hal"} to estimate conditional densities
via the highly adaptive lasso (via \pkg{haldensify}), \code{"sl"} for
\pkg{sl3} learners used to fit Super Learner ensembles to densities via
\pkg{sl3}'s \code{Lrnr_haldensify} or similar, and \code{"external"} for
user-specified input of the form produced by \code{\link{est_g_exp}}.}

\item{g_cens_fit_args}{A \code{list} of arguments, all but one of which are
passed to \code{\link{est_g_cens}}. For details, see the documentation of
\code{\link{est_g_cens}}. The 1st element (i.e., \code{fit_type}) specifies
how this regression is fit: \code{"glm"} for a generalized linear model
or \code{"sl"} for \pkg{sl3} learners used to fit a Super Learner ensemble
for the censoring mechanism, and \code{"external"} for user-specified input
of the form produced by \code{\link{est_g_cens}}.}

\item{Q_fit_args}{A \code{list} of arguments, all but one of which are
passed to \code{\link{est_Q}}. For details, consult the documentation for
\code{\link{est_Q}}. The first element (i.e., \code{fit_type}) is used to
determine how this regression is fit: \code{"glm"} for a generalized linear
model for the outcome mechanism, \code{"sl"} for \pkg{sl3} learners used
to fit a Super Learner for the outcome mechanism, and \code{"external"}
for user-specified input of the form produced by \code{\link{est_Q}}.}

\item{eif_reg_type}{Whether a flexible nonparametric function ought to be
used in the dimension-reduced nuisance regression of the targeting step for
the censored data case. By default, the method used is a nonparametric
regression based on the Highly Adaptive Lasso (from \pkg{hal9001}). Set
this to \code{"glm"} to instead use a simple linear regression model. In
this step, the efficient influence function (EIF) is regressed against
covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).}

\item{ipcw_efficiency}{Whether to use an augmented inverse probability of
censoring weighted EIF estimating equation to ensure efficiency of the
resultant estimate. The default is \code{TRUE}; the inefficient estimation
procedure specified by \code{FALSE} is only supported for completeness.}

\item{samp_fit_ext}{The results of an external fitting procedure used to
estimate the two-phase sampling mechanism, to be used in constructing the
inverse probability of censoring weighted TML or one-step estimator. The
input provided must match the output of \code{\link{est_samp}} exactly.}

\item{gn_exp_fit_ext}{The results of an external fitting procedure used to
estimate the exposure mechanism (generalized propensity score), to be used
in constructing the TML or one-step estimator. The input provided must
match the output of \code{\link{est_g_exp}} exactly.}

\item{gn_cens_fit_ext}{The results of an external fitting procedure used to
estimate the censoring mechanism (propensity score for missingness), to be
used in constructing the TML or one-step estimator. The input provided must
match the output of \code{\link{est_g_cens}} exactly.}

\item{Qn_fit_ext}{The results of an external fitting procedure used to
estimate the outcome mechanism, to be used in constructing the TML or
one-step estimator. The input provided must match the output of
\code{\link{est_Q}} exactly; use of this argument is only recommended for
power users.}
}
\value{
S3 object of class \code{txshift} containing the results of the
 procedure to compute a TML or one-step estimate of the counterfactual mean
 under a modified treatment policy that shifts a continuous-valued exposure
 by a scalar amount \code{delta}. These estimates can be augmented to be
 consistent and efficient when two-phase sampling is performed.
}
\description{
Efficient Estimate of Counterfactual Mean of Stochastic Shift Intervention
}
\details{
Construct a one-step estimate or targeted minimum loss estimate of
 the counterfactual mean under a modified treatment policy, automatically
 making adjustments for two-phase sampling when a censoring indicator is
 included. Ensemble machine learning may be used to construct the initial
 estimates of nuisance functions using \pkg{sl3}.
}
\examples{
set.seed(429153)
n_obs <- 100
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
C_samp <- rbinom(n_obs, 1, plogis(W + Y)) # two-phase sampling
C_cens <- rbinom(n_obs, 1, plogis(rowSums(W) + 0.5))

# construct a TML estimate, ignoring censoring
tmle <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "onestep",
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  )
)
\dontrun{
# construct a TML estimate, accounting for censoring
tmle <- txshift(
  W = W, A = A, C_cens = C_cens, Y = Y, delta = 0.5,
  estimator = "onestep",
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  g_cens_fit_args = list(
    fit_type = "glm",
    glm_formula = "C_cens ~ ."
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  )
)

# construct a TML estimate under two-phase sampling, ignoring censoring
ipcwtmle <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep", max_iter = 3,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  ),
  eif_reg_type = "glm"
)

# construct a TML estimate acconting for two-phase sampling and censoring
ipcwtmle <- txshift(
  W = W, A = A, C_cens = C_cens, Y = Y, delta = 0.5,
  C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep", max_iter = 3,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  g_cens_fit_args = list(
    fit_type = "glm",
    glm_formula = "C_cens ~ ."
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  ),
  eif_reg_type = "glm"
)
}
}
