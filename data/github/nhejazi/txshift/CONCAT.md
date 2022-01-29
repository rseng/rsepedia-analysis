
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

