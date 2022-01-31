
# effectsize <img src="man/figures/logo.png" align="right" width="120" />

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02815/status.svg/)](https://doi.org/10.21105/joss.02815)
[![downloads](https://cranlogs.r-pkg.org/badges/effectsize)](https://cran.r-project.org/package=effectsize/)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/effectsize)](https://cran.r-project.org/package=effectsize/)
[![status](https://tinyverse.netlify.com/badge/effectsize/)](https://CRAN.R-project.org/package=effectsize/)

***Significant is just not enough\!***

The goal of this package is to provide utilities to work with indices of
effect size and standardized parameters, allowing computation and
conversion of indices such as Cohen’s *d*, *r*, odds-ratios, etc.

## Installation

[![CRAN](https://www.r-pkg.org/badges/version/effectsize)](https://cran.r-project.org/package=effectsize/)
[![effectsize status
badge](https://easystats.r-universe.dev/badges/effectsize/)](https://easystats.r-universe.dev/)
[![R-check](https://github.com/easystats/effectsize/workflows/R-check/badge.svg/)](https://github.com/easystats/effectsize/actions/)
[![pkgdown](https://github.com/easystats/effectsize/workflows/pkgdown/badge.svg/)](https://github.com/easystats/effectsize/actions/)
[![Codecov test
coverage](https://codecov.io/gh/easystats/effectsize/branch/main/graph/badge.svg/)](https://app.codecov.io/gh/easystats/effectsize?branch=main/)

Run the following to install the stable release of **effectsize** from
CRAN:

``` r
install.packages("effectsize")
```

Or you can install the latest development version `0.6.0.2` from
[*R-universe*](https://easystats.r-universe.dev):

``` r
install.packages("effectsize", repos = "https://easystats.r-universe.dev/")
```

<!-- Or from *GitHub*: -->

<!-- ```{r, warning=FALSE, message=FALSE, eval=FALSE} -->

<!-- if (!require("remotes")) install.packages("remotes") -->

<!-- remotes::install_github("easystats/effectsize") -->

<!-- ``` -->

## Documentation

[![Documentation](https://img.shields.io/badge/documentation-effectsize-orange.svg?colorB=E91E63/)](https://easystats.github.io/effectsize/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800/)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-effectsize-orange.svg?colorB=2196F3/)](https://easystats.github.io/effectsize/reference/index.html)

Click on the buttons above to access the package
[**documentation**](https://easystats.github.io/effectsize/) and the
[**easystats blog**](https://easystats.github.io/blog/posts/), and
check-out these vignettes:

  - **Effect Sizes**
      - [**Parameter and Model
        Standardization**](https://easystats.github.io/effectsize/articles/standardize_parameters.html)
      - [**ANOVA Effect
        Sizes**](https://easystats.github.io/effectsize/articles/anovaES.html)
      - [**Effect Sizes in Bayesian
        Models**](https://easystats.github.io/effectsize/articles/bayesian_models.html)  
      - [**For Simple Hypothesis
        Tests**](https://easystats.github.io/effectsize/articles/simple_htests.html)  
  - **Effect Sizes Conversion**
      - [**Between Effect
        Sizes**](https://easystats.github.io/effectsize/articles/convert.html)
      - [**Effect Size from Test
        Statistics**](https://easystats.github.io/effectsize/articles/from_test_statistics.html)
  - [**Automated Interpretation of Indices of Effect
    Size**](https://easystats.github.io/effectsize/articles/interpret.html)

# Features

This package is focused on indices of effect size. Check out the package
website for [**a full list of features and functions** provided by
`effectsize`](https://easystats.github.io/effectsize/reference/index.html).

``` r
library(effectsize)
```

## Effect Size Computation

### Standardized Differences (Cohen’s *d*, Hedges’ *g*, Glass’ *delta*)

The package provides functions to compute indices of effect size.

``` r
cohens_d(mpg ~ am, data = mtcars)
## Cohen's d |         95% CI
## --------------------------
## -1.48     | [-2.27, -0.67]
## 
## - Estimated using pooled SD.

hedges_g(mpg ~ am, data = mtcars)
## Hedges' g |         95% CI
## --------------------------
## -1.44     | [-2.21, -0.65]
## 
## - Estimated using pooled SD.

glass_delta(mpg ~ am, data = mtcars)
## Glass' delta |         95% CI
## -----------------------------
## -1.17        | [-1.93, -0.39]
```

`effectsize` also provides effect sizes for *contingency tables*, *rank
tests*, and more…

### ANOVAs (Eta<sup>2</sup>, Omega<sup>2</sup>, …)

``` r
model <- aov(mpg ~ factor(gear), data = mtcars)

eta_squared(model)
## # Effect Size for ANOVA
## 
## Parameter    | Eta2 |       95% CI
## ----------------------------------
## factor(gear) | 0.43 | [0.18, 1.00]
## 
## - One-sided CIs: upper bound fixed at (1).

omega_squared(model)
## # Effect Size for ANOVA
## 
## Parameter    | Omega2 |       95% CI
## ------------------------------------
## factor(gear) |   0.38 | [0.14, 1.00]
## 
## - One-sided CIs: upper bound fixed at (1).

epsilon_squared(model)
## # Effect Size for ANOVA
## 
## Parameter    | Epsilon2 |       95% CI
## --------------------------------------
## factor(gear) |     0.39 | [0.14, 1.00]
## 
## - One-sided CIs: upper bound fixed at (1).
```

And more…

### Regression Models (Standardized Parameters)

Importantly, `effectsize` also provides [advanced
methods](https://easystats.github.io/effectsize/articles/standardize_parameters.html)
to compute standardized parameters for regression models.

``` r
m <- lm(rating ~ complaints + privileges + advance, data = attitude)

standardize_parameters(m)
## # Standardization method: refit
## 
## Parameter   | Coefficient (std.) |        95% CI
## ------------------------------------------------
## (Intercept) |          -9.57e-16 | [-0.22, 0.22]
## complaints  |               0.85 | [ 0.58, 1.13]
## privileges  |              -0.04 | [-0.33, 0.24]
## advance     |              -0.02 | [-0.26, 0.22]
```

Also, models can be re-fit with standardized data:

``` r
standardize(m)
## 
## Call:
## lm(formula = rating ~ complaints + privileges + advance, data = data_std)
## 
## Coefficients:
## (Intercept)   complaints   privileges      advance  
##   -9.57e-16     8.55e-01    -4.35e-02    -2.19e-02
```

## Effect Size Conversion

The package also provides ways of converting between different effect
sizes.

``` r
d_to_r(d = 0.2)
## [1] 0.0995

oddsratio_to_riskratio(2.6, p0 = 0.4)
## [1] 1.59
```

And for recovering effect sizes from test statistics.

``` r
F_to_d(15, df = 1, df_error = 60)
## d    |       95% CI
## -------------------
## 1.00 | [0.46, 1.53]

F_to_r(15, df = 1, df_error = 60)
## r    |       95% CI
## -------------------
## 0.45 | [0.22, 0.61]

F_to_eta2(15, df = 1, df_error = 60)
## Eta2 (partial) |       95% CI
## -----------------------------
## 0.20           | [0.07, 1.00]
## 
## - One-sided CIs: upper bound fixed at (1).
```

## Effect Size Interpretation

The package allows for an automated interpretation of different indices.

``` r
interpret_r(r = 0.3)
## [1] "large"
## (Rules: funder2019)
```

Different sets of “rules of thumb” are implemented ([**guidelines are
detailed
here**](https://easystats.github.io/effectsize/articles/interpret.html))
and can be easily changed.

``` r
interpret_cohens_d(d = 0.45, rules = "cohen1988")
## [1] "small"
## (Rules: cohen1988)

interpret_cohens_d(d = 0.45, rules = "gignac2016")
## [1] "moderate"
## (Rules: gignac2016)
```

### Citation

In order to cite this package, please use the following citation:

  - Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation
    of Effect Size Indices and Standardized Parameters. *Journal of Open
    Source Software*, *5*(56), 2815. doi: 10.21105/joss.02815

Corresponding BibTeX entry:

    @Article{,
      title = {{e}ffectsize: Estimation of Effect Size Indices and Standardized Parameters},
      author = {Mattan S. Ben-Shachar and Daniel Lüdecke and Dominique Makowski},
      year = {2020},
      journal = {Journal of Open Source Software},
      volume = {5},
      number = {56},
      pages = {2815},
      publisher = {The Open Journal},
      doi = {10.21105/joss.02815},
      url = {https://doi.org/10.21105/joss.02815}
    }

# Contributing and Support

If you have any questions regarding the the functionality of the
package, you may either contact us via email or also [file an
issue](https://github.com/easystats/effectsize/issues/). Anyone wishing
to contribute to the package by adding functions, features, or in
another way, please follow [this
guide](https://github.com/easystats/effectsize/blob/main/.github/CONTRIBUTING.md/)
and our [code of
conduct](https://github.com/easystats/effectsize/blob/main/.github/CODE_OF_CONDUCT.md/).
# effectsize 0.6.0.2



# effectsize 0.6.0.1

*This is a patch release.*

## Bug fixes

- `interpret.performance_lavaan()` now works without attaching `effectsize` ( #410 ).  
- `eta_squared()` now fully support multi-variate `car` ANOVAs (class `Anova.mlm`; #406 ).

# effectsize 0.6.0

## Breaking Changes

- `pearsons_c()` effect size column name changed to `Pearsons_c` for consistency. 

## New features

### New API

See [*Support functions for model extensions* vignette](https://easystats.github.io/effectsize/articles/effectsize_API.html).

### Other features

- `eta_squared()` family now supports `afex::mixed()` models.
- `cles()` for estimating common language effect sizes.
- `rb_to_cles()` for converting rank-biserial correlation to Probability of superiority.

## Changes

- `effectsize()` for `BayesFactor` objects returns the same standardized output as for `htest`.

## Bug fixes

- `eta_squared()` for MLM return effect sizes in the correct order of the responses.  
- `eta_squared()` family no longer fails when CIs fail due to non-finite *F*s / degrees of freedom.  
- `standardize()` for multivariate models standardizes the (multivariate) response.
- `standardize()` for models with offsets standardizes offset variables according to `include_response` and `two_sd` ( #396 ).
- `eta_squared()`: fixed a bug that caused `afex_aov` models with more than 2 within-subject factors to return incorrect effect sizes for the lower level factors ( #389 ).

# effectsize 0.5.0

## Breaking Changes

- `cramers_v()` correctly does not work with 1-dimentional tables (for goodness-of-fit tests).
- `interpret_d()`, `interpret_g()`, and `interpret_delta()` are now `interpret_cohens_d()`, `interpret_hedges_g()`, and `interpret_glass_delta()`.
- `interpret_parameters()` was removed. Use `interpret_r()` instead (with caution!).
- Phi, Cohen's *w*, Cramer's *V*, ANOVA effect sizes, rank Epsilon squared, Kendall's *W* - CIs default to 95% one-sided CIs (`alternative = "greater"`). (To restore previous behavior, set `ci = .9, alternative = "two.sided"`.)
- `adjust()`, `change_scale()`, `normalize()`, `ranktransform()`, `standardize()` (data), and `unstandardize()` have moved to the new [`{datawizard}`](https://easystats.github.io/datawizard/) package!

## New features

- `pearsons_c()` (and `chisq_to_pearsons_c()`) for estimating Pearson's contingency coefficient.
- `interpret_vif()` for interpretation of *variance inflation factors*.
- `oddsratio_to_riskratio()` can now convert OR coefficients to RR coefficients from a logistic GLM(M). 
- All effect-size functions gain an `alternative` argument which can be used to make one- or two-sided CIs.
- `interpret()` now accepts as input the results from `cohens_d()`, `eta_squared()`, `rank_biserial()`, etc.
- `interpret_pd()` for the interpretation of the [*Probability of Direction*](https://easystats.github.io/bayestestR/reference/p_direction.html).

## Bug fixes

- `kendalls_w()` CIs now correctly bootstrap samples from the raw data (previously the rank-transformed data was sampled from).
- `cohens_d()`, `sd_pooled()` and `rank_biserial()` now properly respect when `y` is a grouping character vector.
- `effectsize()` for Chi-squared test of goodness-of-fit now correctly respects non-uniform expected probabilities ( #352 ).

## Changes

- `interpret_bf()` now accepts *`log(BF)`* as input.


# effectsize 0.4.5

## New features

- `eta_squared()` family now indicate the type of sum-of-squares used.
- `rank_biserial()` estimates CIs using the normal approximation (previously used bootstrapping).  
- `hedges_g()` now used exact bias correction (thanks to @mdelacre for the suggestion!)  
- `glass_delta()` now estimates CIs using the NCP method based on Algina et al (2006).

## Bug fixes

- `eta_squared()` family returns correctly returns the type 2/3 effect sizes for mixed ANOVAs fit with `afex`.
- `cohens_d()` family now correctly deals with missing factor levels ( #318 )
- `cohens_d()` / `hedges_g()` minor fix for CI with unequal variances.  


## Changes

- `mad_pooled()` (the robust version of `sd_pooled()`) now correctly pools the the two samples.

# effectsize 0.4.4-1

## New features

- `standardize_parameters()` + `eta_sqaured()` support `tidymodels` (when that the underlying model is supported; #311 ).
- `cohens_d()` family now supports `Pairs()` objects as input.
- `standardize_parameters()` gains the `include_response` argument (default to `TRUE`) ( #309 ).

## Bug fixes

- `kendalls_w()` now actually returns correct effect size. Previous estimates were incorrect, and based on transposing the groups and blocks.

# effectsize 0.4.4

`effectsize` now supports `R >= 3.4`.

## New features

- `standardize_parameters()` now supports bootstrapped estimates (from `parameters::bootstrap_model()` and `parameters::bootstrap_parameters()`).
- `unstandardize()` which will reverse the effects of `standardize()`.
- `interpret_kendalls_w()` to interpret Kendall's coefficient of concordance.
- `eta_squared()` family of functions can now also return effect sizes for the intercept by setting `include_intercept = TRUE` ( #156 ).

## Bug fixes

- `standardize()` can now deal with dates ( #300 ).

# effectsize 0.4.3

## Breaking Changes

- `oddsratio()` and `riskratio()` - order of groups has been changed (the
  *first* groups is now the **treatment group**, and the *second* group is the
  **control group**), so that effect sizes are given as *treatment over control*
  (treatment / control) (previously was reversed). This is done to be consistent
  with other functions in R and in `effectsize`.

## New features

- `cohens_h()` effect size for comparing two independent proportions.

- `rank_biserial()`, `cliffs_delta()`, `rank_epsilon_squared()` and
  `kendalls_w()` functions for effect sizes for rank-based tests.

- `adjust()` gains `keep_intercept` argument to keep the intercept.

- `eta_squared()` family of functions supports `Anova.mlm` objects (from the
  `car` package).

- `effectsize()`:

  - supports Cohen's *g* for McNemar's test.

  - Extracts OR from Fisher's Exact Test in the 2x2 case.

- `eta2_to_f2()` / `f2_to_eta2()` to convert between two types of effect sizes
  for ANOVA ( #240 ).

- `cohens_d()` family of functions gain `mu` argument.

## Bug fixes

- `adjust()` properly works when `multilevel = TRUE`.

- `cohens_d()` family / `sd_pooled()` now properly fails when given a missing
  column name.

## Changes

- `effectsize()` for `htest` objects now tries first to extract the data used
  for testing, and computed the effect size directly on that data.

- `cohens_d()` family / `sd_pooled()` now respect any transformations (e.g.
  `I(log(x) - 3) ~ factor(y)`) in a passed formula.

- `eta_squared()` family of functions gains a `verbose` argument.

- `verbose` argument more strictly respected.

- `glass_delta()` returns CIs based on the bootstrap.

# effectsize 0.4.1

## Breaking Changes

- `cohens_d()` and `glass_delta()`: The `correction` argument has been
  deprecated, in favor of it being correctly implemented in `hedges_g()` ( #222
  ).

- `eta_squared_posterior()` no longer uses `car::Anova()` by default.

## New features

- `effectsize()` gains `type = ` argument for specifying which effect size to
  return.

- `eta_squared_posterior()` can return a generalized Eta squared.

- `oddsratio()` and `riskratio()` functions for 2-by-2 contingency tables.

- `standardize()` gains support for `mediation::mediate()` models.

- `eta_squared()` family available for `manova` objects.

## Changes

- `eta_squared()` family of functions returns non-partial effect size for
  one-way between subjects design (#180).

## Bug fixes

- `hedges_g()` correctly implements the available bias correction methods ( #222
  ).

- Fixed width of CI for Cohen's *d* and Hedges' *g* when using *non*-pooled SD.

# effectsize 0.4.0

## Breaking Changes

- `standardize_parameters()` for multi-component models (such as zero-inflated)
  now returns the unstandardized parameters in some cases where standardization
  is not possible (previously returned `NA`s).

- Column name changes:

  - `eta_squared()` / `F_to_eta2` families of function now has the `Eta2`
    format, where previously was `Eta_Sq`.

  - `cramers_v` is now `Cramers_v`

## New features

- `effectsize()` added support for `BayesFactor` objects (Cohen's *d*, Cramer's
  *v*, and *r*).

- `cohens_g()` effect size for paired contingency tables.

- Generalized Eta Squared now available via `eta_squared(generalized = ...)`.

- `eta_squared()`, `omega_squared()` and `epsilon_squared()` fully support
  `aovlist`, `afex_aov` and `mlm` (or `maov`) objects.

- `standardize_parameters()` can now return Odds ratios / IRRs (or any
  exponentiated parameter) by setting `exponentiate = TRUE`.

- Added `cohens_f_squared()` and `F_to_f2()` for Cohen's *f*-squared.

- `cohens_f()` / `cohens_f_squared()`can be used to estimate Cohen's *f* for the
  R-squared change between two models.

- `standardize()` and `standardize_info()` work with weighted models / data (
  #82 ).

- Added `hardlyworking` (simulated) dataset, for use in examples.

- `interpret_*` ( #131 ):

  - `interpret_omega_squared()` added `"cohen1992"` rule.

  - `interpret_p()` added *Redefine statistical significance* rules.

- `oddsratio_to_riskratio()` for converting OR to RR.

## Changes

- CIs for Omega-/Epsilon-squared and Adjusted Phi/Cramer's V return 0s instead
  of negative values.

- `standardize()` for data frames gains the `remove_na` argument for dealing
  with `NA`s ( #147 ).

- `standardize()` and `standardize_info()` now (and by extension,
  `standardize_parameters()`) respect the weights in weighted models when
  standardizing ( #82 ).

- Internal changes to `standardize_parameters()` (reducing co-dependency with
  `parameters`) - argument `parameters` has been dropped.

## Bug fixes

- `ranktransform(sign = TURE)` correctly (doesn't) deal with zeros.

- `effectsize()` for `htest` works with Spearman and Kendall correlations ( #165
  ).

- `cramers_v()` and `phi()` now work with goodness-of-fit data ( #158 )

- `standardize_parameters()` for post-hoc correctly standardizes transformed
  outcome.

- Setting `two_sd = TRUE` in `standardize()` and `standardize_parameters()`
  (correctly) on uses 2-SDs of the predictors (and not the response).

- `standardize_info()` / `standardize_parameters(method = "posthoc")` work for
  zero-inflated models ( #135 )

- `standardize_info(include_pseudo = TRUE)` / `standardize_parameters(method =
  "pseudo")` are less sensitive in detecting between-group variation of
  within-group variables.

- `interpret_oddsratio()` correctly treats extremely small odds the same as
  treats extremely large ones.

# effectsize 0.3.3

## New features

- `standardize_parameters(method = "pseudo")` returns pseudo-standardized
  coefficients for (G)LMM models.

- `d_to_common_language()` for common language measures of standardized
  differences (a-la Cohen's d).

## Changes

- `r_to_odds()` family is now deprecated in favor of `r_to_oddsratio()`.

- `interpret_odds()` is now deprecated in favor of `interpret_oddsratio()`

## Bug fixes

- `phi()` and `cramers_v()` did not respect the CI argument ( #111 ).

- `standardize()` / `standardize_parameters()` properly deal with transformed
  data in the model formula ( #113 ).

- `odds_to_probs()` was mis-treating impossible odds (NEVER TELL ME THE ODDS!
  #123 )

# effectsize 0.3.2

## New features

- `eta_squared_posterior()` for estimating Eta Squared for Bayesian models.

- `eta_squared()`, `omega_squared()` and `epsilon_squared()` now works with

  - `ols` / `rms` models.

- `effectsize()` for class `htest` supports `oneway.test(...)`.

## Bug fixes

- Fix minor miss-calculation of Chi-squared for 2*2 table with small samples (
  #102 ).

- Fixed miss-calculation of signed rank in `ranktransform()` ( #87 ).

- Fixed bug in `standardize()` for standard objects with non-standard
  class-attributes (like vectors of class `haven_labelled` or `vctrs_vctr`).

- Fix `effectsize()` for one sample `t.test(...)` ( #95 ; thanks to pull request
  by @mutlusun )

# effectsize 0.3.1

## New features

- `standardize_parameters()` now returns CIs ( #72 )

- `eta_squared()`, `omega_squared()` and `epsilon_squared()` now works with

  - `gam` models.

  - `afex` models.

  - `lme` and `anova.lme` objects.

- New function `equivalence_test()` for effect sizes.

- New plotting methods in the `see` package.

# effectsize 0.3.0

## New features

- New general purpose `effectsize()` function.

- Effectsize for differences have CI methods, and return a data frame.

- Effectsize for ANOVA all have CI methods, and none are based on
  bootstrapping.

- New effect sizes for contingency tables (`phi()` and `cramers_v()`).

- `chisq_to_phi()` / `cramers_v()` functions now support CIs (via the ncp
  method), and return a data frame.

- `F_to_eta2()` family of functions now support CIs (via the ncp method), and
  return a data frame.

- `t_to_d()` and `t_to_r()` now support CIs (via the ncp method), and return a
  data frame.

- `standardize()` for model-objects has a default-method, which usually accepts
  all models. Exception for model-objects that do not work will be added if
  missing.

- `standardize.data.frame()` gets `append` and `suffix` arguments, to add
  (instead of replace) standardized variables to the returned data frame.

- `eta_squared()`, `omega_squared()` and `epsilon_squared()` now works

  - output from `parameters::model_parameters()`.

  - `mlm` models.

## Bug fixes

- Fix `cohens_d()`'s dealing with formula input (#44).

- `sd_pooled()` now returns the... pooled sd (#44).

## Changes

- In `t_to_d()`, argument `pooled` is now `paired`.

# effectsize 0.2.0

## Bug fixes

- `standardize.data.frame()` did not work when variables had missing values.

- Fixed wrong computation in `standardize()` when `two_sd = TRUE`.

- Fixed bug with missing column names in `standardize_parameters()` for models
  with different components (like count and zero-inflation).

# effectsize 0.1.1

## Changes

- News are hidden in an air of mystery...

# effectsize 0.1.0

## New features

- `standardize_parameters()` and `standardize()` now support models from
  packages *brglm*, *brglm2*, *mixor*, *fixest*, *cgam*, *cplm*, *cglm*,
  *glmmadmb* and *complmrob*.

## Bug fixes

- Fix CRAN check issues.

All URL issues have been resolved.

## Test environments

* local installation: R 4.1.1 on Windows
* GitHub Actions
    - Windows:        devel, release, oldrel
    - macOS:          devel, release, oldrel
    - ubuntu-16.04:   devel, release, oldrel, 3.6, 3.5, 3.4
* win-builder:        release


## R CMD check results

0 errors | 0 warnings | 0 notes


### Known issues

- Failed handshake with *shinyapps.io* is a false positive.


## revdepcheck results

We checked 16 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
---
title: "effectsize: Estimation of Effect Size Indices and Standardized Parameters"
authors:
- affiliation: 1
  name: Mattan S. Ben-Shachar
  orcid: 0000-0002-4287-4801
- affiliation: 2
  name: Daniel Lüdecke
  orcid: 0000-0002-8895-3206
- affiliation: 3
  name: Dominique Makowski
  orcid: 0000-0001-5375-9967

date: "28 October, 2020"
output: 
  pdf_document:
    latex_engine: xelatex
bibliography: paper.bib
csl: apa.csl
tags:
- R
- easystats
- effect size
- regression
- linear models
- standardized coefficients
affiliations:
- index: 1
  name: Ben-Gurion University of the Negev, Israel
- index: 2
  name:  University Medical Center Hamburg-Eppendorf, Germany
- index: 3
  name: Nanyang Technological University, Singapore
---

# Aims of the Package

In both theoretical and applied research, it is often of interest to assess the strength of an observed association. This is typically done to allow the judgment of the magnitude of an effect [especially when units of measurement are not meaningful, e.g., in the use of estimated latent variables; @bollen1989structural], to facilitate comparing between predictors' importance within a given model, or both. Though some indices of effect size, such as the correlation coefficient (itself a standardized covariance coefficient) are readily available, other measures are often harder to obtain. **effectsize** is an R package [@rcore] that fills this important gap, providing utilities for easily estimating a wide variety of standardized effect sizes (i.e., effect sizes that are not tied to the units of measurement of the variables of interest) and their confidence intervals (CIs), from a variety of statistical models. **effectsize** provides easy-to-use functions, with full documentation and explanation of the various effect sizes offered, and is also used by developers of other R packages as the back-end for effect size computation, such as **parameters** [@ludecke2020extracting], **ggstatsplot** [@patil2020ggstatsplot], **gtsummary** [@sjoberg2020gtsummary] and more.

# Comparison to Other Packages

**effectsize**'s functionality is in part comparable to packages like **lm.beta** [@behrendt2014lmbeta], **MOTE** [@buchanan2019MOTE], and **MBESS** [@kelley2020MBESS]. Yet, there are some notable differences, e.g.:

- **lm.beta** provides standardized regression coefficients for linear models, based on post-hoc model matrix standardization. However, the functionality is available only for a limited number of models (models inheriting from the `lm` class), whereas **effectsize** provides support for many types of models, including (generalized) linear mixed models, Bayesian models, and more. Additionally, in additional to post-hoc model matrix standardization, **effectsize** offers other methods of standardization (see below).  
- Both **MOTE** and **MBESS** provide functions for computing effect sizes such as Cohen's *d* and effect sizes for ANOVAs [@cohen1988statistical], and their confidence intervals. However, both require manual input of *F*- or *t*-statistics, *degrees of freedom*, and *sums of squares* for the computation the effect sizes, whereas **effectsize** can automatically extract this information from the provided models, thus allowing for better ease-of-use as well as reducing any potential for error.  
- Finally, in **base R**, the function `scale()` can be used to standardize vectors, matrices and data frame, which can be used to standardize data prior to model fitting. The coefficients of a linear model fit on such data are in effect standardized regression coefficients. **effectsize** expands an this, allowing for robust standardization (using the median and the MAD, instead of the mean and SD), post-hoc parameter standardization, and more.

# Examples of Features

**effectsize** provides various functions for extracting and estimating effect sizes and their confidence intervals [estimated using the noncentrality parameter method; @steiger2004beyond]. In this article, we provide basic usage examples for estimating some of the most common effect size. A comprehensive overview, including in-depth examples and [a full list of features and functions](https://easystats.github.io/effectsize/reference/index.html), are accessible via a dedicated website (https://easystats.github.io/effectsize/).

## Indices of Effect Size

### Standardized Differences

**effectsize** provides functions for estimating the common indices of standardized differences such as Cohen's *d* (`cohens_d()`), Hedges' *g* (`hedges_g()`) for both paired and independent samples [@cohen1988statistical; @hedges1985statistical], and Glass' $\Delta$ (`glass_delta()`) for independent samples with different variances [@hedges1985statistical].

``` r
cohens_d(mpg ~ am, data = mtcars)
#> Cohen's d |         95% CI
#> --------------------------
#>     -1.48 | [-2.27, -0.67]
#>
#> - Estimated using pooled SD.

```

### Contingency Tables

Pearson's $\phi$ (`phi()`) and Cramér's *V* (`cramers_v()`) can be used to estimate the strength of association between two categorical variables [@cramer1946mathematical], while Cohen's *g* (`cohens_g()`) estimates the deviance between paired categorical variables [@cohen1988statistical].

``` r
M <- rbind(c(150, 130, 35, 55),
           c(100, 50,  10, 40),
           c(165, 65,  2,  25))

cramers_v(M)
#> Cramer's V |       95% CI
#> -------------------------
#>       0.18 | [0.12, 0.22]
```

## Parameter and Model Standardization

Standardizing parameters (i.e., coefficients) can allow for their comparison within and between models, variables and studies. To this end, two functions are available: `standardize()`, which returns an updated model, re-fit with standardized data, and `standardize_parameters()`, which returns a table of standardized coefficients from a provided model [for a list of supported models, see the *insight* package; @luedecke2019insight].

``` r
model <- lm(mpg ~ cyl * am, 
            data = mtcars)

standardize(model)
#> 
#> Call:
#> lm(formula = mpg ~ cyl * am, data = data_std)
#> 
#> Coefficients:
#> (Intercept)          cyl           am       cyl:am  
#>     -0.0977      -0.7426       0.1739      -0.1930 

standardize_parameters(model)
#> Parameter   | Coefficient (std.) |         95% CI
#> -------------------------------------------------
#> (Intercept) |              -0.10 | [-0.30,  0.11]
#> cyl         |              -0.74 | [-0.95, -0.53]
#> am          |               0.17 | [-0.04,  0.39]
#> cyl:am      |              -0.19 | [-0.41,  0.02]
#> 
#> # Standardization method: refit
```

Standardized parameters can also be produced for generalized linear models (GLMs; where only the predictors are standardized):

``` r
model <- glm(am ~ cyl + hp,
             family = "binomial",
             data = mtcars)

standardize_parameters(model, exponentiate = TRUE)
#> Parameter   | Odds Ratio (std.) |        95% CI
#> -----------------------------------------------
#> (Intercept) |              0.53 | [0.18,  1.32]
#> cyl         |              0.05 | [0.00,  0.29]
#> hp          |              6.70 | [1.32, 61.54]
#> 
#> # Standardization method: refit
```

`standardize_parameters()` provides several standardization methods, such as robust standardization, or *pseudo*-standardized coefficients for (generalized) linear mixed models [@hoffman2015longitudinal]. A full review of these methods can be found in the [*Parameter and Model Standardization* vignette](https://easystats.github.io/effectsize/articles/standardize_parameters.html).

## Effect Sizes for ANOVAs

Unlike standardized parameters, the effect sizes reported in the context of ANOVAs (analysis of variance) or ANOVA-like tables represent the amount of variance explained by each of the model's terms, where each term can be represented by one or more parameters. `eta_squared()` can produce such popular effect sizes as Eta-squared ($\eta^2$), its partial version ($\eta^2_p$), as well as the generalized $\eta^2_G$ [@cohen1988statistical; @olejnik2003generalized]:


``` r
options(contrasts = c('contr.sum', 'contr.poly'))

data("ChickWeight")
# keep only complete cases and convert `Time` to a factor
ChickWeight <- subset(ChickWeight, ave(weight, Chick, FUN = length) == 12)
ChickWeight$Time <- factor(ChickWeight$Time)

model <- aov(weight ~ Diet * Time + Error(Chick / Time),
             data = ChickWeight) 

eta_squared(model, partial = TRUE)
#> Group      | Parameter | Eta2 (partial) |       90% CI
#> ------------------------------------------------------
#> Chick      |      Diet |           0.27 | [0.06, 0.42]
#> Chick:Time |      Time |           0.87 | [0.85, 0.88]
#> Chick:Time | Diet:Time |           0.22 | [0.11, 0.23]

eta_squared(model, generalized = "Time")
#> Group      | Parameter | Eta2 (generalized) |       90% CI
#> ----------------------------------------------------------
#> Chick      |      Diet |               0.04 | [0.00, 0.09]
#> Chick:Time |      Time |               0.74 | [0.71, 0.77]
#> Chick:Time | Diet:Time |               0.03 | [0.00, 0.00]
```

**effectsize** also offers $\epsilon^2_p$ (`epsilon_squared()`) and $\omega^2_p$ (`omega_squared()`), which are less biased estimates of the variance explained in the population [@kelley1935unbiased; @olejnik2003generalized]. For more details about the various effect size measures and their applications, see the [*Effect sizes for ANOVAs* vignette](https://easystats.github.io/effectsize/articles/anovaES.html).

## Effect Size Conversion

### From Test Statistics

In many real world applications there are no straightforward ways of obtaining standardized effect sizes. However, it is possible to get approximations of most of the effect size indices (*d*, *r*, $\eta^2_p$...) with the use of test statistics [@friedman1982simplified]. These conversions are based on the idea that test statistics are a function of effect size and sample size (or more often of degrees of freedom). Thus it is possible to reverse-engineer indices of effect size from test statistics (*F*, *t*, $\chi^2$, and *z*). 

``` r
F_to_eta2(f = c(40.72, 33.77),
          df = c(2, 1), df_error = c(18, 9))
#> Eta2 (partial) |       90% CI
#> -----------------------------
#>           0.82 | [0.66, 0.89]
#>           0.79 | [0.49, 0.89]

t_to_d(t = -5.14, df_error = 22)
#>     d |         95% CI
#> ----------------------
#> -2.19 | [-3.23, -1.12]

t_to_r(t = -5.14, df_error = 22)
#>     r |         95% CI
#> ----------------------
#> -0.74 | [-0.85, -0.49]
```

These functions also power the `effectsize()` convenience function for estimating effect sizes from R's `htest`-type objects. For example:

``` r
aov1 <- oneway.test(salary ~ n_comps, 
                    data = hardlyworking, var.equal = TRUE)
effectsize(aov1)
#> Eta2 |       90% CI
#> -------------------
#> 0.20 | [0.14, 0.24]

xtab <- rbind(c(762, 327, 468), c(484, 239, 477), c(484, 239, 477))
Xsq <- chisq.test(xtab)
effectsize(Xsq)
#> Cramer's V |       95% CI
#> -------------------------
#>       0.07 | [0.05, 0.09]
```

These functions also power our *Effect Sizes From Test Statistics* shiny app (https://easystats4u.shinyapps.io/statistic2effectsize/).

### Between Effect Sizes

For comparisons between different types of designs and analyses, it is useful to be able to convert between different types of effect sizes [*d*, *r*, Odds ratios and Risk ratios; @borenstein2009converting; @grant2014converting].

``` r
r_to_d(0.7)
#> [1] 1.960392

d_to_oddsratio(1.96)
#> [1] 34.98946

oddsratio_to_riskratio(34.99, p0 = 0.4)
#> [1] 2.397232

oddsratio_to_r(34.99)
#> [1] 0.6999301
```

## Effect Size Interpretation

Finally, **effectsize** provides convenience functions to apply existing or custom interpretation rules of thumb, such as for instance Cohen's (1988). Although we strongly advocate for the cautious and parsimonious use of such judgment-replacing tools, we provide these functions to allow users and developers to explore and hopefully gain a deeper understanding of the relationship between data values and their interpretation. More information is available in the [*Automated Interpretation of Indices of Effect Size* vignette](https://easystats.github.io/effectsize/articles/interpret.html).

``` r
interpret_d(c(0.02, 0.52, 0.86), rules = "cohen1988")
#> [1] "very small" "medium"     "large"     
#> (Rules: cohen1988)
```


# Licensing and Availability

**effectsize** is licensed under the GNU General Public License (v3.0), with all source code stored at GitHub (https://github.com/easystats/effectsize), and with a corresponding issue tracker for bug reporting and feature enhancements. In the spirit of honest and open science, we encourage requests/tips for fixes, feature updates, as well as general questions and concerns via direct interaction with contributors and developers, by [filing an issue](https://github.com/easystats/effectsize/issues/). See the package's [*Contribution Guidelines*](https://github.com/easystats/effectsize/blob/main/.github/CONTRIBUTING.md/).

# Acknowledgments

**effectsize** is part of the [*easystats*](https://github.com/easystats/easystats/) ecosystem, a collaborative project created to facilitate the usage of R for statistical analyses. Thus, we would like to thank the [members of easystats](https://github.com/orgs/easystats/people/) as well as the users.

# References
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at matanshm@post.bgu.ac.il. 
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at <https://www.contributor-covenant.org/version/2/0/code_of_conduct.html>.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
<https://www.contributor-covenant.org/faq>. Translations are available at <https://www.contributor-covenant.org/translations>.
# Contribution Guidelines 

<sup>easystats guidelines 0.1.0</sup>

**All people are very much welcome to contribute to code, documentation, testing and suggestions.**

This package aims at being beginner-friendly. Even if you're new to this open-source way of life, new to coding and github stuff, we encourage you to try submitting pull requests (PRs). 

- **"I'd like to help, but I'm not good enough with programming yet"**

It's alright, don't worry! You can always dig in the code, in the documentation or tests. There are always some typos to fix, some docs to improve, some details to add, some code lines to document, some tests to add... **Even the smaller PRs are appreciated**.

- **"I'd like to help, but I don't know where to start"**

You can look around the **issue section** to find some features / ideas / bugs to start working on. You can also open a new issue **just to say that you're there, interested in helping out**. We might have some ideas adapted to your skills.

- **"I'm not sure if my suggestion or idea is worthwile"**

Enough with the impostor syndrom! All suggestions and opinions are good, and even if it's just a thought or so, it's always good to receive feedback.

- **"Why should I waste my time with this? Do I get any credit?"**

Software contributions are getting more and more valued in the academic world, so it is a good time to collaborate with us! Authors of substantial contributions will be added within the **authors** list. We're also very keen on including them to eventual academic publications.


**Anyway, starting is the most important! You will then enter a *whole new world, a new fantastic point of view*... So fork this repo, do some changes and submit them. We will then work together to make the best out of it :)**


## Code

- Please document and comment your code, so that the purpose of each step (or code line) is stated in a clear and understandable way.
- Before submitting a change, please read the [**R style guide**](https://style.tidyverse.org/) and in particular our [**easystats convention of code-style**](https://github.com/easystats/easystats#convention-of-code-style) to keep some consistency in code formatting.
- Regarding the style guide, note this exception: we put readability and clarity before everything. Thus, we like underscores and full names (prefer `model_performance` over `modelperf` and `interpret_odds_logistic` over `intoddslog`).
- Before you start to code, make sure you're on the `dev` branch (the most "advanced"). Then, you can create a new branch named by your feature (e.g., `feature_lightsaber`) and do your changes. Finally, submit your branch to be merged into the `dev` branch. Then, every now and then, the dev branch will merge into `main`, as a new package version.

## Checks to do before submission

- Make sure **documentation** (roxygen) is good
- Make sure to add **tests** for the new functions
- Run:

  - `styler::style_pkg()`: Automatic style formatting
  - `lintr::lint_package()`: Style checks
  - `devtools::check()`: General checks



## Useful Materials

- [Understanding the GitHub flow](https://guides.github.com/introduction/flow/)


# Description

This PR aims at adding this feature...

# Proposed Changes

I changed the `foo` function so that ...
---
name: Bug report
about: Create a report to help us improve

---

**Describe the bug**
A description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

[conciser using `reprex::reprex()`]

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Specifiations (please complete the following information):**
 - `R` Version [version]
 - `effectsize` Version [version]
  
---
name: Feature idea
about: Suggest an idea for this project

---

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**How could we do it?**
A description of actual ways of implementing a feature.
---
name: Question
about: You didn't understand something?

---

**Question and context**
---
output: 
  github_document:
    toc: false
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, standardized]
vignette: >
  %\VignetteIndexEntry{README}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
csl: vignettes/apa.csl
editor_options: 
  chunk_output_type: console
---

# effectsize <img src="man/figures/logo.png" align="right" width="120" />

```{r setup, echo = FALSE, warning=FALSE, message=FALSE}
library(effectsize)

options(digits=3)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi=450,
  fig.path = "man/figures/"
)

set.seed(111)
```

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02815/status.svg/)](https://doi.org/10.21105/joss.02815)
[![downloads](https://cranlogs.r-pkg.org/badges/effectsize)](https://cran.r-project.org/package=effectsize/)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/effectsize)](https://cran.r-project.org/package=effectsize/)
[![status](https://tinyverse.netlify.com/badge/effectsize/)](https://CRAN.R-project.org/package=effectsize/)


***Significant is just not enough!***

The goal of this package is to provide utilities to work with indices of effect size and standardized parameters, allowing computation and conversion of indices such as Cohen's *d*, *r*, odds-ratios, etc.


## Installation

[![CRAN](https://www.r-pkg.org/badges/version/effectsize)](https://cran.r-project.org/package=effectsize/)
[![effectsize status badge](https://easystats.r-universe.dev/badges/effectsize/)](https://easystats.r-universe.dev/)
[![R-check](https://github.com/easystats/effectsize/workflows/R-check/badge.svg/)](https://github.com/easystats/effectsize/actions/)
[![pkgdown](https://github.com/easystats/effectsize/workflows/pkgdown/badge.svg/)](https://github.com/easystats/effectsize/actions/)
[![Codecov test coverage](https://codecov.io/gh/easystats/effectsize/branch/main/graph/badge.svg/)](https://app.codecov.io/gh/easystats/effectsize?branch=main/)

Run the following to install the stable release of **effectsize** from CRAN:

```{r install-CRAN, warning=FALSE, message=FALSE, eval=FALSE}
install.packages("effectsize")
```

Or you can install the latest development version ``r available.packages(repos = "https://easystats.r-universe.dev/")["effectsize","Version"]`` from [*R-universe*](https://easystats.r-universe.dev):

```{r install-R-universe, warning=FALSE, message=FALSE, eval=FALSE}
install.packages("effectsize", repos = "https://easystats.r-universe.dev/")
```

<!-- Or from *GitHub*: -->

<!-- ```{r, warning=FALSE, message=FALSE, eval=FALSE} -->
<!-- if (!require("remotes")) install.packages("remotes") -->
<!-- remotes::install_github("easystats/effectsize") -->
<!-- ``` -->

## Documentation

[![Documentation](https://img.shields.io/badge/documentation-effectsize-orange.svg?colorB=E91E63/)](https://easystats.github.io/effectsize/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800/)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-effectsize-orange.svg?colorB=2196F3/)](https://easystats.github.io/effectsize/reference/index.html)

Click on the buttons above to access the package [**documentation**](https://easystats.github.io/effectsize/) and the [**easystats blog**](https://easystats.github.io/blog/posts/), and check-out these vignettes:

- **Effect Sizes**  
  - [**Parameter and Model Standardization**](https://easystats.github.io/effectsize/articles/standardize_parameters.html)
  - [**ANOVA Effect Sizes**](https://easystats.github.io/effectsize/articles/anovaES.html)
  - [**Effect Sizes in Bayesian Models**](https://easystats.github.io/effectsize/articles/bayesian_models.html)  
  - [**For Simple Hypothesis Tests**](https://easystats.github.io/effectsize/articles/simple_htests.html)  
- **Effect Sizes Conversion**    
  - [**Between Effect Sizes**](https://easystats.github.io/effectsize/articles/convert.html)
  - [**Effect Size from Test Statistics**](https://easystats.github.io/effectsize/articles/from_test_statistics.html)
- [**Automated Interpretation of Indices of Effect Size**](https://easystats.github.io/effectsize/articles/interpret.html)



# Features

This package is focused on indices of effect size. Check out the package website for [**a full list of features and functions** provided by `effectsize`](https://easystats.github.io/effectsize/reference/index.html).

```{r load, message=FALSE, warning=FALSE}
library(effectsize)
```

## Effect Size Computation

### Standardized Differences (Cohen's *d*, Hedges' *g*, Glass' *delta*)

The package provides functions to compute indices of effect size.

```{r d, warning=FALSE, message=FALSE}
cohens_d(mpg ~ am, data = mtcars)

hedges_g(mpg ~ am, data = mtcars)

glass_delta(mpg ~ am, data = mtcars)
```

`effectsize` also provides effect sizes for *contingency tables*, *rank tests*, and more...

### ANOVAs (Eta<sup>2</sup>, Omega<sup>2</sup>, ...)

```{r aov, warning=FALSE, message=FALSE}
model <- aov(mpg ~ factor(gear), data = mtcars)

eta_squared(model)

omega_squared(model)

epsilon_squared(model)
```

And more...


### Regression Models (Standardized Parameters)

Importantly, `effectsize` also provides [advanced methods](https://easystats.github.io/effectsize/articles/standardize_parameters.html) to compute standardized parameters for regression models.

```{r beta, warning=FALSE, message=FALSE}
m <- lm(rating ~ complaints + privileges + advance, data = attitude)

standardize_parameters(m)
```

Also, models can be re-fit with standardized data:

```{r std-model, warning=FALSE, message=FALSE}
standardize(m)
```

## Effect Size Conversion

The package also provides ways of converting between different effect sizes.

```{r convert-between, warning=FALSE, message=FALSE}
d_to_r(d = 0.2)

oddsratio_to_riskratio(2.6, p0 = 0.4)
```

And for recovering effect sizes from test statistics.

```{r convert-stat, warning=FALSE, message=FALSE}
F_to_d(15, df = 1, df_error = 60)

F_to_r(15, df = 1, df_error = 60)

F_to_eta2(15, df = 1, df_error = 60)
```

## Effect Size Interpretation

The package allows for an automated interpretation of different indices. 

```{r interp-r, warning=FALSE, message=FALSE}
interpret_r(r = 0.3)
```

Different sets of "rules of thumb" are implemented ([**guidelines are detailed here**](https://easystats.github.io/effectsize/articles/interpret.html)) and can be easily changed.


```{r interp-d, warning=FALSE, message=FALSE}
interpret_cohens_d(d = 0.45, rules = "cohen1988")

interpret_cohens_d(d = 0.45, rules = "gignac2016")
```


### Citation

In order to cite this package, please use the following citation:

  * Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. _Journal of Open Source Software_, *5*(56), 2815. doi: 10.21105/joss.02815
  
Corresponding BibTeX entry:

```
@Article{,
  title = {{e}ffectsize: Estimation of Effect Size Indices and Standardized Parameters},
  author = {Mattan S. Ben-Shachar and Daniel Lüdecke and Dominique Makowski},
  year = {2020},
  journal = {Journal of Open Source Software},
  volume = {5},
  number = {56},
  pages = {2815},
  publisher = {The Open Journal},
  doi = {10.21105/joss.02815},
  url = {https://doi.org/10.21105/joss.02815}
}
```

# Contributing and Support

If you have any questions regarding the the functionality of the package, you may either contact us via email or also [file an issue](https://github.com/easystats/effectsize/issues/). Anyone wishing to contribute to the package by adding functions, features, or in another way, please follow [this guide](https://github.com/easystats/effectsize/blob/main/.github/CONTRIBUTING.md/) and our [code of conduct](https://github.com/easystats/effectsize/blob/main/.github/CODE_OF_CONDUCT.md/).
---
title: "Effect sizes for ANOVAs"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, ANOVA]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{Effect sizes for ANOVAs}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---


```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)
knitr::opts_chunk$set(comment = ">", warning = FALSE)

set.seed(1)
pkgs <- c("effectsize", "parameters", "car", "afex")
if (!all(sapply(pkgs, require, quietly = TRUE, character.only = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}
```

## Eta<sup>2</sup>

In the context of ANOVA-like tests, it is common to report ANOVA-like effect sizes. Unlike [standardized parameters](https://easystats.github.io/effectsize/articles/standardize_parameters.html), these effect sizes represent the amount of variance explained by each of the model's terms, where each term can be represented by 1 *or more* parameters.

For example, in the following case, the parameters for the `treatment` term represent specific contrasts between the factor's levels (treatment groups) - the difference between each level and the reference level (`obk.long == 'control'`).

```{r}
data(obk.long, package = "afex")
# modify the data slightly for the demonstration:
obk.long <- obk.long[1:240 %% 3 == 0, ]
obk.long$id <- seq_len(nrow(obk.long))

m <- lm(value ~ treatment, data = obk.long)

parameters::model_parameters(m)
```

But we can also ask about the overall effect of `treatment` - how much of the
variation in our dependent variable `value` can be predicted by (or explained
by) the variation between the `treatment` groups. Such a question can be
answered with an ANOVA test:

```{r}
parameters::model_parameters(anova(m))
```

As we can see, the variance in `value` (the *sums-of-squares*, or *SS*) has been
split into pieces:

- The part associated with `treatment`.
- The unexplained part (The Residual-*SS*).

We can now ask what is the percent of the total variance in `value` that is
associated with `treatment`. This measure is called Eta-squared (written as
$\eta^2$):

$$
\eta^2 = \frac{SS_{effect}}{SS_{total}} = \frac{72.23}{72.23 + 250.96} = 0.22
$$

and can be accessed via the `eta_squared()` function:

```{r}
library(effectsize)

eta_squared(m, partial = FALSE)
```


### Adding More Terms

When we add more terms to our model, we can ask two different questions about
the percent of variance explained by a predictor - how much variance is
accounted by the predictor in *total*, and how much is accounted when
*controlling* for any other predictors. The latter questions is answered by the
*partial*-Eta squared ($\eta^2_p$), which is the percent of the **partial**
variance (after accounting for other predictors in the model) associated with a
term:

$$
\eta^2_p = \frac{SS_{effect}}{SS_{effect} + SS_{error}}
$$
which can also be accessed via the `eta_squared()` function:

```{r}
m <- lm(value ~ gender + phase + treatment, data = obk.long)

eta_squared(m, partial = FALSE)

eta_squared(m) # partial = TRUE by default
```

*(`phase` is a repeated-measures variable, but for simplicity it is not modeled as such.)*

In the calculation above, the *SS*s were computed sequentially - that is the
*SS* for `phase` is computed after controlling for `gender`, and the *SS* for
`treatment` is computed after controlling for both `gender` and `phase`. This
method of sequential *SS* is called also *type-I* test. If this is what you
want, that's great - however in many fields (and other statistical programs) it
is common to use "simultaneous" sums of squares (*type-II* or *type-III* tests),
where each *SS* is computed controlling for all other predictors, regardless of
order. This can be done with `car::Anova(type = ...)`:

```{r}
eta_squared(car::Anova(m, type = 2), partial = FALSE)

eta_squared(car::Anova(m, type = 3)) # partial = TRUE by default
```

$\eta^2_p$ will always be larger than $\eta^2$. The idea is to simulate the
effect size in a design where only the term of interest was manipulated. This
terminology assumes some causal relationship between the predictor and the
outcome, which reflects the experimental world from which these analyses and
measures hail; However, $\eta^2_p$ can also simply be seen as a
**signal-to-noise- ratio**, as it only uses the term's *SS* and the error-term's
*SS*.[^in repeated-measure designs the term-specific residual-*SS* is used for
the computation of the effect size].

(Note that in a one-way fixed-effect designs $\eta^2 = \eta^2_p$.)

### Adding Interactions

Type II and type III treat interaction differently.
Without going into the weeds here, keep in mind that **when using type III SS, it is important to center all of the predictors**;
for numeric variables this can be done by mean-centering the predictors;
for factors this can be done by using orthogonal coding (such as `contr.sum` for *effects-coding*) for the dummy variables (and *NOT* treatment coding, which is the default in R).
This unfortunately makes parameter interpretation harder, but *only* when this is does do the *SS*s associated with each lower-order term (or lower-order interaction) represent the ***SS*** of the **main effect** (with treatment coding they represent the *SS* of the simple effects).


```{r}
# compare
m_interaction1 <- lm(value ~ treatment * gender, data = obk.long)

# to:
m_interaction2 <- lm(
  value ~ treatment * gender,
  data = obk.long,
  contrasts = list(
    treatment = "contr.sum",
    gender = "contr.sum"
  )
)

eta_squared(car::Anova(m_interaction1, type = 3))
eta_squared(car::Anova(m_interaction2, type = 3))
```

If all of this type-III-effects-coding seems like a hassle, you can use the `afex` package, which takes care of all of this behind the scenes:

```{r}
library(afex)
m_afex <- aov_car(value ~ treatment * gender + Error(id), data = obk.long)

eta_squared(m_afex)
```


## Other Measures of Effect Size

### Unbiased Effect Sizes

These effect sizes are unbiased estimators of the population's $\eta^2$:

- **Omega Squared** ($\omega^2$)
- **Epsilon Squared** ($\epsilon^2$), also referred to as *Adjusted Eta Squared*.

```{r}
omega_squared(m_afex)

epsilon_squared(m_afex)
```


Both $\omega^2$ and $\epsilon^2$ (and their partial counterparts, $\omega^2_p$ &
$\epsilon^2_p$) are unbiased estimators of the population's $\eta^2$ (or
$\eta^2_p$, respectively), which is especially important is small samples.
Though $\omega^2$ is the more popular choice [@albers2018power], $\epsilon^2$ is
analogous to adjusted-$R^2$ [@allen2017statistics, p. 382], and has been found
to be less biased [@carroll1975sampling].

### Generalized Eta<sup>2</sup>

*Partial* Eta squared aims at estimating the effect size in a design where only
the term of interest was manipulated, assuming all other terms are have also
manipulated. However, not all predictors are always manipulated - some can only
be observed. For such cases, we can use *generalized* Eta squared ($\eta^2_G$),
which like $\eta^2_p$ estimating the effect size in a design where only the term
of interest was manipulated, accounting for the fact that some terms cannot be
manipulated (and so their variance would be present in such a design).

```{r}
eta_squared(m_afex, generalized = "gender")
```

$\eta^2_G$ is useful in repeated-measures designs, as it can estimate what a
*within-subject* effect size would have been had that predictor been manipulated
*between-subjects* [@olejnik2003generalized].


### Cohen's *f*

Finally, we have the forgotten child - Cohen's $f$. Cohen's $f$ is a
transformation of $\eta^2_p$, and is the ratio between the term-*SS* and the
error-*SS*.

$$\text{Cohen's} f_p = \sqrt{\frac{\eta^2_p}{1-\eta^2_p}} = \sqrt{\frac{SS_{effect}}{SS_{error}}}$$

It can take on values between zero, when the population means are all equal, and
an indefinitely large number as the means are further and further apart. It is
analogous to Cohen's $d$ when there are only two groups.

```{r}
cohens_f(m_afex)
```

## When Sum-of-Squares are Hard to Come By

Until now we've discusses effect sizes in fixed-effect linear model and
repeated-measures ANOVA's - cases where the *SS*s are readily available, and so
the various effect sized presented can easily be estimated. How ever this is not
always the case.

For example, in linear mixed models (LMM/HLM/MLM), the estimation of all required *SS*s is not straightforward. However, we can still *approximate* these effect sizes (only their partial versions) based on the **test-statistic approximation method** (learn more in the [*Effect Size from Test Statistics* vignette](https://easystats.github.io/effectsize/articles/from_test_statistics.html)).

```{r, eval=require(lmerTest)}
library(lmerTest)

fit_lmm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

anova(fit_lmm) # note the type-3 errors

F_to_eta2(45.8, df = 1, df_error = 17)
```

Or directly with `eta_squared() and co.:


```{r, eval=require(lmerTest)}
eta_squared(fit_lmm)
epsilon_squared(fit_lmm)
omega_squared(fit_lmm)
```

Another case where *SS*s are not available is when use Bayesian models. `effectsize` has Bayesian solutions for Bayesian models, about which you can read in the [*Effect Sizes for Bayesian Models* vignette](https://easystats.github.io/effectsize/articles/bayesian_models.html).


# References
---
title: "Effect Sizes: Getting Started"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, rules of thumb, guidelines, conversion]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{effectsize}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Aims of the Package

In both theoretical and applied research, it is often of interest to assess the strength of an observed association. This is typically done to allow the judgment of the magnitude of an effect [especially when units of measurement are not meaningful, e.g., in the use of estimated latent variables; @bollen1989structural], to facilitate comparing between predictors' importance within a given model, or both. Though some indices of effect size, such as the correlation coefficient (itself a standardized covariance coefficient) are readily available, other measures are often harder to obtain. **effectsize** is an R package [@rcore] that fills this important gap, providing utilities for easily estimating a wide variety of standardized effect sizes (i.e., effect sizes that are not tied to the units of measurement of the variables of interest) and their confidence intervals (CIs), from a variety of statistical models. **effectsize** provides easy-to-use functions, with full documentation and explanation of the various effect sizes offered, and is also used by developers of other R packages as the back-end for effect size computation, such as **parameters** [@ludecke2020extracting], **ggstatsplot** [@patil2020ggstatsplot], **gtsummary** [@sjoberg2020gtsummary] and more.

# Comparison to Other Packages

**effectsize**'s functionality is in part comparable to packages like **lm.beta** [@behrendt2014lmbeta], **MOTE** [@buchanan2019MOTE], and **MBESS** [@kelley2020MBESS]. Yet, there are some notable differences, e.g.:

- **lm.beta** provides standardized regression coefficients for linear models, based on post-hoc model matrix standardization. However, the functionality is available only for a limited number of models (models inheriting from the `lm` class), whereas **effectsize** provides support for many types of models, including (generalized) linear mixed models, Bayesian models, and more. Additionally, in additional to post-hoc model matrix standardization, **effectsize** offers other methods of standardization (see below).  
- Both **MOTE** and **MBESS** provide functions for computing effect sizes such as Cohen's *d* and effect sizes for ANOVAs [@cohen1988statistical], and their confidence intervals. However, both require manual input of *F*- or *t*-statistics, *degrees of freedom*, and *sums of squares* for the computation the effect sizes, whereas **effectsize** can automatically extract this information from the provided models, thus allowing for better ease-of-use as well as reducing any potential for error.  
<!-- - Finally, in **base R**, the function `scale()` can be used to standardize vectors, matrices and data frame, which can be used to standardize data prior to model fitting. The coefficients of a linear model fit on such data are in effect standardized regression coefficients. **effectsize** expands an this, allowing for robust standardization (using the median and the MAD, instead of the mean and SD), post-hoc parameter standardization, and more. -->

# Examples of Features

**effectsize** provides various functions for extracting and estimating effect sizes and their confidence intervals [estimated using the noncentrality parameter method; @steiger2004beyond]. In this article, we provide basic usage examples for estimating some of the most common effect size. A comprehensive overview, including in-depth examples and [a full list of features and functions](https://easystats.github.io/effectsize/reference/index.html), are accessible via a dedicated website (https://easystats.github.io/effectsize/).

## Indices of Effect Size

### Standardized Differences

**effectsize** provides functions for estimating the common indices of standardized differences such as Cohen's *d* (`cohens_d()`), Hedges' *g* (`hedges_g()`) for both paired and independent samples [@cohen1988statistical; @hedges1985statistical], and Glass' $\Delta$ (`glass_delta()`) for independent samples with different variances [@hedges1985statistical].

```{r}
library(effectsize)

cohens_d(mpg ~ am, data = mtcars)
```

### Contingency Tables

Pearson's $\phi$ (`phi()`) and Cramér's *V* (`cramers_v()`) can be used to estimate the strength of association between two categorical variables [@cramer1946mathematical], while Cohen's *g* (`cohens_g()`) estimates the deviance between paired categorical variables [@cohen1988statistical].

```{r}
M <- rbind(c(150, 130, 35, 55),
           c(100, 50,  10, 40),
           c(165, 65,  2,  25))

cramers_v(M)
```

## Parameter and Model Standardization

Standardizing parameters (i.e., coefficients) can allow for their comparison within and between models, variables and studies. To this end, two functions are available: `standardize()`, which returns an updated model, re-fit with standardized data, and `standardize_parameters()`, which returns a table of standardized coefficients from a provided model [for a list of supported models, see the *insight* package; @luedecke2019insight].

```{r}
model <- lm(mpg ~ cyl * am, 
            data = mtcars)

standardize(model)

standardize_parameters(model)
```

Standardized parameters can also be produced for generalized linear models (GLMs; where only the predictors are standardized):

```{r}
model <- glm(am ~ cyl + hp,
             family = "binomial",
             data = mtcars)

standardize_parameters(model, exponentiate = TRUE)
```

`standardize_parameters()` provides several standardization methods, such as robust standardization, or *pseudo*-standardized coefficients for (generalized) linear mixed models [@hoffman2015longitudinal]. A full review of these methods can be found in the [*Parameter and Model Standardization* vignette](https://easystats.github.io/effectsize/articles/standardize_parameters.html).

## Effect Sizes for ANOVAs

Unlike standardized parameters, the effect sizes reported in the context of ANOVAs (analysis of variance) or ANOVA-like tables represent the amount of variance explained by each of the model's terms, where each term can be represented by one or more parameters. `eta_squared()` can produce such popular effect sizes as Eta-squared ($\eta^2$), its partial version ($\eta^2_p$), as well as the generalized $\eta^2_G$ [@cohen1988statistical; @olejnik2003generalized]:

```{r}
options(contrasts = c('contr.sum', 'contr.poly'))

data("ChickWeight")
# keep only complete cases and convert `Time` to a factor
ChickWeight <- subset(ChickWeight, ave(weight, Chick, FUN = length) == 12)
ChickWeight$Time <- factor(ChickWeight$Time)

model <- aov(weight ~ Diet * Time + Error(Chick / Time),
             data = ChickWeight) 

eta_squared(model, partial = TRUE)

eta_squared(model, generalized = "Time")
```

**effectsize** also offers $\epsilon^2_p$ (`epsilon_squared()`) and $\omega^2_p$ (`omega_squared()`), which are less biased estimates of the variance explained in the population [@kelley1935unbiased; @olejnik2003generalized]. For more details about the various effect size measures and their applications, see the [*Effect sizes for ANOVAs* vignette](https://easystats.github.io/effectsize/articles/anovaES.html).

## Effect Size Conversion

### From Test Statistics

In many real world applications there are no straightforward ways of obtaining standardized effect sizes. However, it is possible to get approximations of most of the effect size indices (*d*, *r*, $\eta^2_p$...) with the use of test statistics [@friedman1982simplified]. These conversions are based on the idea that test statistics are a function of effect size and sample size (or more often of degrees of freedom). Thus it is possible to reverse-engineer indices of effect size from test statistics (*F*, *t*, $\chi^2$, and *z*). 

```{r}
F_to_eta2(f = c(40.72, 33.77),
          df = c(2, 1), df_error = c(18, 9))

t_to_d(t = -5.14, df_error = 22)

t_to_r(t = -5.14, df_error = 22)
```

These functions also power the `effectsize()` convenience function for estimating effect sizes from R's `htest`-type objects. For example:

```{r}
data(hardlyworking, package = "effectsize")

aov1 <- oneway.test(salary ~ n_comps, 
                    data = hardlyworking, var.equal = TRUE)
effectsize(aov1)

xtab <- rbind(c(762, 327, 468), c(484, 239, 477), c(484, 239, 477))
Xsq <- chisq.test(xtab)
effectsize(Xsq)
```

These functions also power our *Effect Sizes From Test Statistics* shiny app (https://easystats4u.shinyapps.io/statistic2effectsize/).

### Between Effect Sizes

For comparisons between different types of designs and analyses, it is useful to be able to convert between different types of effect sizes [*d*, *r*, Odds ratios and Risk ratios; @borenstein2009converting; @grant2014converting].

```{r}
r_to_d(0.7)

d_to_oddsratio(1.96)

oddsratio_to_riskratio(34.99, p0 = 0.4)

oddsratio_to_r(34.99)
```

## Effect Size Interpretation

Finally, **effectsize** provides convenience functions to apply existing or custom interpretation rules of thumb, such as for instance Cohen's (1988). Although we strongly advocate for the cautious and parsimonious use of such judgment-replacing tools, we provide these functions to allow users and developers to explore and hopefully gain a deeper understanding of the relationship between data values and their interpretation. More information is available in the [*Automated Interpretation of Indices of Effect Size* vignette](https://easystats.github.io/effectsize/articles/interpret.html).

```{r}
interpret_cohens_d(c(0.02, 0.52, 0.86), rules = "cohen1988")
```


# Licensing and Availability

**effectsize** is licensed under the GNU General Public License (v3.0), with all source code stored at GitHub (https://github.com/easystats/effectsize), and with a corresponding issue tracker for bug reporting and feature enhancements. In the spirit of honest and open science, we encourage requests/tips for fixes, feature updates, as well as general questions and concerns via direct interaction with contributors and developers, by [filing an issue](https://github.com/easystats/effectsize/issues/). See the package's [*Contribution Guidelines*](https://github.com/easystats/effectsize/blob/main/.github/CONTRIBUTING.md/).

# Acknowledgments

**effectsize** is part of the [*easystats*](https://github.com/easystats/easystats/) ecosystem, a collaborative project created to facilitate the usage of R for statistical analyses. Thus, we would like to thank the [members of easystats](https://github.com/orgs/easystats/people/) as well as the users.

# References
---
title: "Effect Sizes for Simple Hypothesis Tests"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, rules of thumb, guidelines, conversion]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{Effect Sizes for Simple Hypothesis Tests}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r setup, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
knitr::opts_chunk$set(comment = ">")
options(digits = 3)

pkgs <- c("effectsize", "BayesFactor")
if (!all(sapply(pkgs, require, quietly = TRUE, character.only = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}
set.seed(7)
```

This vignette provides a short review of effect sizes for common hypothesis
tests (in **`R`** these are usually achieved with various `*.test()`
functions).

```{r}
library(effectsize)
library(BayesFactor)
```

In most cases, the effect sizes can be automagically extracted from the `htest`
object via the `effectsize()` function.

## Standardized Differences

For *t*-tests, it is common to report an effect size representing a standardized
difference between the two compared samples' means. These measures range from
$-\infty$ to $+\infty$, with negative values indicating the second group's mean
is larger (and vice versa).

### Two Independent Samples

For two independent samples, the difference between the means is standardized
based on the pooled standard deviation of both samples (assumed to be equal in
the population):

```{r}
t.test(mpg ~ am, data = mtcars, var.equal = TRUE)

cohens_d(mpg ~ am, data = mtcars)
```

Hedges' *g* provides a bias correction for small sample sizes ($N < 20$).

```{r}
hedges_g(mpg ~ am, data = mtcars)
```

If variances cannot be assumed to be equal, it is possible to get estimates that
are not based on the pooled standard deviation:

```{r}
t.test(mpg ~ am, data = mtcars, var.equal = FALSE)

cohens_d(mpg ~ am, data = mtcars, pooled_sd = FALSE)

hedges_g(mpg ~ am, data = mtcars, pooled_sd = FALSE)
```

In cases where the differences between the variances are substantial, it is also
common to standardize the difference based only on the standard deviation of one
of the groups (usually the "control" group); this effect size is known as Glass'
$\Delta$ (delta) (Note that the standard deviation is taken from the *second* sample).

```{r}
glass_delta(mpg ~ am, data = mtcars)
```

For a one-sided hypothesis, it is also possible to construct one-sided confidence intervals:

```{r}
t.test(mpg ~ am, data = mtcars, var.equal = TRUE, alternative = "less")

cohens_d(mpg ~ am, data = mtcars, pooled_sd = TRUE, alternative = "less")
```

#### Common Language Effect Sizes

Related effect sizes are the *common language effect sizes* which present information about group differences in terms of probability.

```{r}
cles(mpg ~ am, data = mtcars)
```

### One Sample and Paired Samples

In the case of a one-sample test, the effect size represents the standardized
distance of the mean of the sample from the null value. For paired-samples, the
difference between the paired samples is used:

```{r}
t.test(extra ~ group, data = sleep, paired = TRUE)

cohens_d(extra ~ group, data = sleep, paired = TRUE)

hedges_g(extra ~ group, data = sleep, paired = TRUE)
```

### For a Bayesian *t*-test

```{r}
(BFt <- ttestBF(mtcars$mpg[mtcars$am == 0], mtcars$mpg[mtcars$am == 1]))

effectsize(BFt, test = NULL)
```

## One way ANOVA

For more details, see [ANOVA
vignette](https://easystats.github.io/effectsize/articles/anovaES.html).

```{r, message=FALSE}
onew <- oneway.test(mpg ~ gear, data = mtcars, var.equal = TRUE)

eta_squared(onew)
```

## Contingency Tables and Proportions

For contingency tables Cramér's *V*, $\phi$ (Phi, also known as Cohen's *w*) and Pearson's contingency coefficient indicate the strength of association with 0 indicating no association between the variables.
While Cramér's *V* and Pearson's *C* are capped at 1 (perfect association), $\phi$ can be larger than 1.

```{r}
(Music <- matrix(
  c(
    150, 130, 35, 55,
    100, 50, 10, 40,
    165, 65, 2, 25
  ),
  byrow = TRUE, nrow = 3,
  dimnames = list(
    Study = c("Psych", "Econ", "Law"),
    Music = c("Pop", "Rock", "Jazz", "Classic")
  )
))

chisq.test(Music)

cramers_v(Music)

phi(Music)

pearsons_c(Music)
```

Pearson's *C* and $\phi$ are also applicable to tests of goodness-of-fit, where small values indicate no deviation from the hypothetical probabilities and large values indicate... large deviation from the hypothetical probabilities.

```{r}
O <- c(89,  37,  130, 28,  2) # observed group sizes
E <- c(.40, .20, .20, .15, .05) # expected group freq

chisq.test(O, p = E, rescale.p = TRUE)

pearsons_c(O, p = E, rescale.p = TRUE)

phi(O, p = E, rescale.p = TRUE)
```

These can also be extracted from the equivalent Bayesian test:

```{r}
(BFX <- contingencyTableBF(Music, sampleType = "jointMulti"))

effectsize(BFX, type = "cramers_v", test = NULL)

effectsize(BFX, type = "phi", test = NULL)

effectsize(BFX, type = "pearsons_c", test = NULL)
```

### Comparing Two Proportions (2x2 tables)

For $2\times 2$ tables, in addition to Cramér's *V*, $\phi$ and Pearson's *C*, we can also
compute the Odds-ratio (OR), where each column represents a different group.
Values larger than 1 indicate that the odds are higher in the first group (and
vice versa).

```{r}
(RCT <- matrix(
  c(
    71, 30,
    50, 100
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    Diagnosis = c("Sick", "Recovered"),
    Group = c("Treatment", "Control")
  )
))

chisq.test(RCT) # or fisher.test(RCT)

oddsratio(RCT)
```

We can also compute the Risk-ratio (RR), which is the ratio between the
proportions of the two groups - a measure which some claim is more intuitive.

```{r}
riskratio(RCT)
```

Additionally, Cohen's *h* can also be computed, which uses the *arcsin*
transformation. Negative values indicate smaller proportion in the first group
(and vice versa).

```{r}
cohens_h(RCT)
```

### Paired Contingency Tables

For dependent (paired) contingency tables, Cohen's *g* represents the symmetry
of the table, ranging between 0 (perfect symmetry) and 0.5 (perfect asymmetry).

```{r}
(Performance <- matrix(
  c(
    794, 86,
    150, 570
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    "1st Survey" = c("Approve", "Disapprove"),
    "2nd Survey" = c("Approve", "Disapprove")
  )
))

mcnemar.test(Performance)

cohens_g(Performance)
```

## Rank Based tests

Rank based tests get rank based effect sizes!

### Difference in Ranks

For two independent samples, the rank-biserial correlation ($r_{rb}$) is a
measure of relative superiority - i.e., larger values indicate a higher
probability of a randomly selected observation from *X* being larger than
randomly selected observation from *Y*. A value of $(-1)$ indicates that all
observations in the second group are larger than the first, and a value of
$(+1)$ indicates that all observations in the first group are larger than the
second.

```{r, warning=FALSE}
A <- c(48, 48, 77, 86, 85, 85)
B <- c(14, 34, 34, 77)

wilcox.test(A, B) # aka Mann–Whitney U test

rank_biserial(A, B)
```

Here too we have a *common language effect size*:

```{r}
cles(A, B, rank = TRUE)
```

For one sample, $r_{rb}$ measures the symmetry around $\mu$ (mu; the null
value), with 0 indicating perfect symmetry, $(-1)$ indicates that all
observations fall below $\mu$, and $(+1)$ indicates that all observations fall
above $\mu$. For paired samples the difference between the paired samples is
used:

```{r}
x <- c(1.15, 0.88, 0.90, 0.74, 1.21, 1.36, 0.89)

wilcox.test(x, mu = 1) # aka Signed-Rank test

rank_biserial(x, mu = 1)


x <- c(1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

wilcox.test(x, y, paired = TRUE) # aka Signed-Rank test

rank_biserial(x, y, paired = TRUE)
```

### Rank One way ANOVA

The Rank-Epsilon-Squared ($\varepsilon^2$) is a measure of association
for the rank based one-way ANOVA. Values range between 0 (no relative
superiority between any of the groups) to 1 (complete separation - with no
overlap in ranks between the groups).

```{r}
group_data <- list(
  g1 = c(2.9, 3.0, 2.5, 2.6, 3.2), # normal subjects
  g2 = c(3.8, 2.7, 4.0, 2.4), # with obstructive airway disease
  g3 = c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
)

kruskal.test(group_data)

rank_epsilon_squared(group_data)
```

### Rank One way Repeated-Measures ANOVA

For a rank based repeated measures one-way ANOVA, Kendall's *W* is a measure of
agreement on the effect of condition between various "blocks" (the subjects), or
more often conceptualized as a measure of reliability of the rating / scores of
observations (or "groups") between "raters" ("blocks").

```{r}
# Subjects are COLUMNS
(ReactionTimes <- matrix(
  c(398, 338, 520,
    325, 388, 555,
    393, 363, 561,
    367, 433, 470,
    286, 492, 536,
    362, 475, 496,
    253, 334, 610),
  nrow = 7, byrow = TRUE,
  dimnames = list(
    paste0("Subject", 1:7),
    c("Congruent", "Neutral", "Incongruent")
  )
))

friedman.test(ReactionTimes)

kendalls_w(ReactionTimes)
```

# References
---
title: "Effect Sizes for Bayesian Models"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, bayesian, effect size]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{Effect Sizes for Bayesian Models}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)
knitr::opts_chunk$set(comment = ">")

set.seed(1)
pkgs <- c("effectsize", "parameters", "rstanarm", "bayestestR", "car")
if (!all(sapply(pkgs, require, quietly = TRUE, character.only = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}
```

## Standardized Parameters

### Introduction

Like in OLS / ML or other frequentists methods of model parameter estimation,
standardizing the parameters of Bayesian (generalized) linear regression models
can allow for the comparison of so-called "effects" within and between models,
variables and studies.

As with frequentists methods, standardizing parameters should not be the only
method of examining the role different predictors play in a particular Bayesian
model, and this vignette generally assumes that the issues of model convergence,
goodness of fit and model selection have already been taken care of. (Learn more
about how to become a Bayesian master with [the `bayestestR`
package](https://easystats.github.io/bayestestR/).)

### Setup

We will examine the predictive role of overtime (`xtra_hours`), number of
compliments given to the boss (`n_comps`) and seniority in predicting workers
salaries. Let's fit the model:

```{r, warning=FALSE}
library(rstanarm)

data("hardlyworking", package = "effectsize")

head(hardlyworking)


mod <- stan_glm(salary ~ xtra_hours + n_comps + seniority,
  data = hardlyworking,
  prior = normal(0, scale = c(1, 0.5, 0.5), autoscale = TRUE), # set some priors
  refresh = 0
)

parameters::model_parameters(mod, test = NULL)
```

Looking at the un-standardized ("raw") parameters, it looks like all predictors
positively predict workers' salaries, but which has the highest predictive
power? Unfortunately, the predictors are not on the same scale (hours,
compliments, years), so comparing them is hard when looking at the raw data.
This is where standardization comes in.

Like with [frequentists
models](https://easystats.github.io/effectsize/articles/standardize_parameters.html)
we can choose from the same standardization methods. Let's use the (slow)
`"refit"` method.

```{r}
library(effectsize)

standardize_parameters(mod, method = "refit", ci = 0.89)
```

Note that the central tendency of the posterior distribution is still the
*median* - the median of the standardized posterior distribution. We can easily
change this, of the type of credible interval used:

```{r}
library(effectsize)

standardize_parameters(mod,
  method = "basic", ci = 0.89,
  centrality = "MAP", ci_method = "eti"
)
```

As we can see, working harder (or at least for longer hours) has stronger
predictive power than complementing or seniority. (Do note, however, that this
does not mean that if you wish to have a higher salary you should work overtime
- the raw parameters seem to suggest that complementing your boss is the way to
go, with one compliment worth almost 3.5 times **more** than a full hours'
work!)

## Eta<sup>2</sup>

### Introduction

In classical frequentists models, the computation of $\eta^2$ or $\eta^2_p$ is
straightforward: based on the right combinations of sums-of-squares (*SS*s), we
get the correct proportion of variance accounted for by some predictor term.
However such a computation is not as straightforward for Bayesian models, for
various reasons (e.g., the model-*SS* and the residual-*SS* don't necessarily
sum to the total-*SS*). Although some have proposed Bayesian methods of
estimating explained variance in ANOVA designs [@marsman2019bayesian], these are
not yet easy to implement with `stan`-based models.

An alternative route to obtaining effect sizes of explained variance, is via the
use of the ***posterior predictive distribution*** (*PPD*). The PPD is the
Bayesian expected distribution of possible unobserved values. Thus, after
observing some data, we can estimate not just the expected mean values (the
conditional marginal means), but also the full *distribution* of data around
these values [@gelman2014bayesian, chapter 7].

By sampling from the PPD, we can decompose the sample to the various *SS*s
needed for the computation of explained variance measures. By repeatedly
sampling from the PPD, we can generate a posterior distribution of explained
variance estimates. But note that **these estimates are conditioned not only on
the location-parameters of the model, but also on the scale-parameters of the
model!** So it is vital to [validate the
PPD](https://mc-stan.org/docs/2_23/stan-users-guide/meta-models-part.html#meta-models.part/)
before using it to estimate explained variance measures.

### Setup

Let's factorize out data from above:

```{r}
hardlyworking$age_f <- cut(hardlyworking$age,
  breaks = c(25, 35, 45), right = FALSE,
  labels = c("Young", "Less_young")
)
hardlyworking$comps_f <- cut(hardlyworking$n_comps,
  breaks = c(0, 1, 2, 3),
  include.lowest = TRUE,
  right = FALSE
)

table(hardlyworking$age_f, hardlyworking$comps_f)
```

And fit our model:

```{r}
# use (special) effects coding
contrasts(hardlyworking$age_f) <- bayestestR::contr.bayes
contrasts(hardlyworking$comps_f) <- bayestestR::contr.bayes

modAOV <- stan_glm(salary ~ age_f * comps_f,
  data = hardlyworking, family = gaussian(),
  refresh = 0
)
```

We can use `eta_squared_posterior()` to get the posterior distribution of
$eta^2$ or $eta^2_p$ for each effect. Like an ANOVA table, we must make sure to
use the right effects-coding and *SS*-type:

```{r}
pes_posterior <- eta_squared_posterior(modAOV,
  draws = 500, # how many samples from the PPD?
  partial = TRUE, # partial eta squared
  # type 3 SS
  ss_function = car::Anova, type = 3
)

head(pes_posterior)

bayestestR::describe_posterior(pes_posterior, rope_range = c(0, 0.1), test = "rope")
```

Compare to:

```{r}
modAOV_f <- lm(salary ~ age_f * comps_f,
  data = hardlyworking
)

eta_squared(car::Anova(modAOV_f, type = 3))
```

# References
---
title: "Additional Vignettes"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    toc: false
vignette: >
  %\VignetteIndexEntry{additional vignettes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Additional Vignettes

Due to constraints on the size of the `R` package, all available vignettes are
only available on the website for this package: <br>
<https://easystats.github.io/effectsize/>---
title: "Support functions for model extensions"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, ANOVA, standardization, standardized coefficients]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{Support functions for model extensions}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
knitr::opts_chunk$set(comment = ">",
                      warning = FALSE,
                      message = FALSE)
options(digits = 2)
options(knitr.kable.NA = '')


pkgs <- c("effectsize")
if (!all(sapply(pkgs, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

set.seed(333)
```

```{r}
library(effectsize)
```

## Supporting ANOVA Effect Sizes

To add support for you model, create a new `.anova_es()` method function. This functions should generally do 3 things:

1. Build a data frame with all the required information.
2. Pass the data frame to one of the 3 functions.
3. Set some attributes to the output.

### Simple ANOVA tables

The input data frame must have these columns:
- `Parameter` (char) - The name of the parameter or, more often, the term.
- `Sum_Squares` (num) - The sum of squares.
- `df` (num) - The degrees of freedom associated with the `Sum_Squares`.
- `Mean_Square_residuals` (num; *optional*) - if *not* present, is calculated as `Sum_Squares / df`.
(Any other column is ignored.)

And exactly *1* row Where `Parameter` is `Residual`.

Optionally, one of the rows can have a `(Intercept)` value for `Parameter`.

An example of a minimally valid data frame:

```{r}
min_aov <- data.frame(
  Parameter = c("(Intercept)", "A", "B", "Residuals"),
  Sum_Squares = c(30, 40, 10, 100),
  df = c(1, 1, 2, 50)
)
```

Pass the data frame to `.es_aov_simple()`:

```{r}
.es_aov_simple(
  min_aov,
  type = "eta", partial = TRUE, generalized = FALSE,
  include_intercept = FALSE,
  ci = 0.95, alternative = "greater",
  verbose = TRUE
)
```

The output is a data frame with the columns: `Parameter`, the effect size, and (optionally) `CI` + `CI_low` + `CI_high`,

And with the following attributes: `partial`, `generalized`, `ci`, `alternative`, `anova_type` (`NA` or `NULL`), `approximate`.

You can then set the `anova_type` attribute to {1, 2, 3, or `NA`} and return the output.

### ANOVA Tables with Multiple Error Strata

(e.g., `aovlist` models.)

The input data frame must have these columns:

- `Group` (char) - The strata
- `Parameter` (char)
- `Sum_Squares` (num)
- `df` (num)
- `Mean_Square_residuals` (num; *optional*)

And exactly *1* row ***per `Group`*** Where `Parameter` is `Residual`.

Optionally, one of the rows can have a `(Intercept)` value for `Parameter`.

An example of a minimally valid data frame:

```{r}
min_aovlist <- data.frame(
  Group = c("S", "S", "S:A", "S:A"),
  Parameter = c("(Intercept)", "Residuals", "A", "Residuals"),
  Sum_Squares = c(34, 21, 34, 400),
  df = c(1, 12, 4, 30)
)
```

Pass the data frame to `.es_aov_strata()`, along with a list of predictors (including the stratifying variables) to the `DV_names` argument:

```{r}
.es_aov_strata(
  min_aovlist, DV_names = c("S", "A"),
  type = "omega", partial = TRUE, generalized = FALSE,
  ci = 0.95, alternative = "greater",
  verbose = TRUE,
  include_intercept = TRUE
)
```

The output is a data frame with the columns: `Group`, `Parameter`, the effect size, and (optionally) `CI` + `CI_low` + `CI_high`,

And with the following attributes: `partial`, `generalized`, `ci`, `alternative`, `approximate`.

You can then set the `anova_type` attribute to {1, 2, 3, or `NA`} and return the output.


### Approximate Effect sizes

When *sums of squares* cannot be extracted, we can still get *approximate* effect sizes based on the `F_to_eta2()` family of functions.

The input data frame must have these columns:

- `Parameter` (char)
- `F` (num) - The *F* test statistic.
- `df` (num) - effect degrees of freedom.
- (Can also have a `t` col instead, in which case `df` is set to 1, and `F` is `t^2`).
- `df_error` (num) - error degrees of freedom.

Optionally, one of the rows can have `(Intercept)` as the `Parameter`.

An example of a minimally valid data frame:

```{r}
min_anova <- data.frame(
  Parameter = c("(Intercept)", "A", "B"),
  F = c(4, 7, 0.7),
  df = c(1, 1, 2),
  df_error = 34
)
```

Pass the table to `.es_aov_table()`:

```{r}
.es_aov_table(
  min_anova,
  type = "eta", partial = TRUE, generalized = FALSE,
  include_intercept = FALSE,
  ci = 0.95, alternative = "greater",
  verbose = TRUE
)
```

The output is a data frame with the columns: `Parameter`, the effect size, and (optionally) `CI` + `CI_low` + `CI_high`,

And with the following attributes: `partial`, `generalized`, `ci`, `alternative`, `approximate`.

You can then set the `anova_type` attribute to {1, 2, 3, or `NA`} and return the output, and optionally the `approximate` attribute, and return the output.

### *Example*

Let's fit a simple linear model and change its class:

```{r}
mod <- lm(mpg ~ factor(cyl) + am, mtcars)

class(mod) <- "superMODEL"
```

We now need a new `.anova_es.superMODEL` function:

```{r}
.anova_es.superMODEL <- function(model, ...) {
  # Get ANOVA table
  anov <- suppressWarnings(stats:::anova.lm(model))
  anov <- as.data.frame(anov)
  
  # Clean up
  anov[["Parameter"]] <- rownames(anov)
  colnames(anov)[2:1] <- c("Sum_Squares", "df")
  
  # Pass
  out <- .es_aov_simple(anov, ...)
  
  # Set attribute
  attr(out, "anova_type") <- 1
  
  out
}
```

And... that's it! Our new `superMODEL` class of models is fully supported!

```{r}
eta_squared(mod)

eta_squared(mod, partial = FALSE)

omega_squared(mod)

# Etc...
```


## Supporting Model Re-Fitting with Standardized Data

`effectsize::standardize.default()` should support your model if you have methods for:

1. `{insight}` functions.
2. An `update()` method that can take the model and a data frame via the `data = ` argument.

Or you can make your own `standardize.my_class()` function, DIY-style (possibly using `datawizard::standardize.data.frame()` or `datawizard::standardize.numeric()`). This function should return a fiffed model of the same class as the input model.

## Supporting Standardized Parameters

`standardize_parameters.default()` offers a few methods of parameter standardization:

- For `method = "refit"` all you need is to have `effectsize::standardize()` support (see above) as well as `parameters::model_parameters()`.  
- ***API for post-hoc methods coming soon...***  

<!-- `standardize_parameters.default()` should support your model if it is already supported by `{parameters}` and `{insight}`. -->

<!-- - For `method = "refit"`, to have `effectsize::standardize()` support (see above). -->
<!-- - For the post-hoc methods, you will need to have a method for `standardize_info()` (or use the default method). See next section. -->

<!-- Or you can make your own `standardize_parameters.my_class()` and/or `standardize_info.my_class()` functions. -->

<!-- ## Extracting Post-Hoc Standardization Information (`standardize_info`) -->

<!-- The `standardize_info()` function computes the standardized units needed for standardization; In order to standardize some slope $b_{xi}$, we need to multiply it by a scaling factor: -->

<!-- $$b^*_{xi} = \frac{\text{Deviation}_{xi}}{\text{Deviation}_{y}}\times b_{xi}$$ -->

<!-- These "deviations" are univariate scaling factors of the response and the specific parameter (usin its corresponding feature in the design matrix). Most often these are a single standard deviation (*SD*), but depending on the `robust` and `two_sd` arguments, these can be also be two *MAD*s, etc. -->

<!-- Let's look at an example: -->

<!-- ```{r} -->
<!-- m <- lm(mpg ~ factor(cyl) * am, data = mtcars) -->

<!-- standardize_info(m) -->
<!-- ``` -->

<!-- - The first 4 columns (`Parameter`, `Type`, `Link`, `Secondary_Parameter`) are taken from `parameters::parameters_type()`.   -->
<!-- - The `EffectSize_Type` column is not used here, but is used in the the `{report}` package.   -->
<!-- - `Deviation_Response_Basic` and `Deviation_Response_Smart` correspond to the $\text{Deviation}_{y}$ scalar using two different methods of post-hoc standardization (see `standardize_parameters()` docs for more details).   -->
<!--     - Note then when the response is not standardized (either due to `standardize_parameters(include_response = FALSE)` or because the model uses a non-continuous response), both methods are fixed at **1** (i.e., no standardization with respect to the outcome).   -->
<!-- - `Deviation_Basic` and `Deviation_Smart` correspond to the $\text{Deviation}_{xi}$ scaler using two different methods of post-hoc standardization. -->

<!-- This information is then used by the `standardize_parameters()` to standardize the parameters. -->
    

# References

---
title: "Parameter and Model Standardization"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, standardization, effect size, cohen d, standardized coefficients]
vignette: >
  %\VignetteIndexEntry{Parameter and Model Standardization}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
knitr::opts_chunk$set(comment = ">",
                      warning = FALSE,
                      message = FALSE)
options(digits = 2)
options(knitr.kable.NA = '')

pkgs <- c("effectsize", "parameters", "correlation")
if (!all(sapply(pkgs, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

set.seed(333)
```

<!-- centering and interactions! -->

# Introduction

Standardizing parameters (*i.e.*, coefficients) can allow for their comparison
within and between models, variables and studies. Moreover, as it returns
coefficients expressed in terms of **change of variance** (for instance,
coefficients expressed in terms of SD of the response variable), it can allow
for the usage of [effect size interpretation
guidelines](https://easystats.github.io/effectsize/articles/interpret.html),
such as Cohen's (1988) famous rules of thumb.

However, standardizing a model's parameters should *not* be automatically and
mindlessly done: for some research fields, particular variables or types of
studies (*e.g.*, replications), it sometimes makes more sense to keep, use and
interpret the original parameters, especially if they are well known or easily
understood.

Critically, **parameters standardization is not a trivial process**. Different
techniques exist, that can lead to drastically different results. Thus, it is
critical that the standardization method is explicitly documented and detailed.

<!-- **`parameters` include different techniques of parameters
standardization**, described below
[@bring1994standardize;@menard2004six;@gelman2008scaling;@schielzeth2010simple;@menard2011standards].
-->

## Standardizing Parameters of Simple Models

### Standardized Associations

```{r}
library(effectsize)

m <- lm(rating ~ complaints, data = attitude)

standardize_parameters(m)
```

Standardizing the coefficient of this *simple* linear regression gives a value
of `0.87`, but did you know that for a simple regression this is actually the
**same as a correlation**? Thus, you can eventually apply some (*in*)famous
interpretation guidelines (e.g., Cohen's rules of thumb).

```{r}
correlation::correlation(attitude, select = c("rating", "complaints"))
```

### Standardized Differences

How does it work in the case of differences, when **factors** are entered and
differences between a given level and a reference level? You might have heard
that it is similar to a **Cohen's *d***. Well, let's see.

```{r include=FALSE}
mtcars <- datasets::mtcars
```

```{r}
# Select portion of data containing the two levels of interest
mtcars$am <- factor(mtcars$am, labels = c("Manual", "Automatic"))

m <- lm(mpg ~ am, data = mtcars)
standardize_parameters(m)
```

This linear model suggests that the *standardized* difference between *Manual*
(the reference level - the model's intercept) and *Automatic* is of 1.20
standard deviation of `mpg` (because the response variable was standardized,
right?). Let's compute the **Cohen's *d*** between these two levels:

```{r}
cohens_d(mpg ~ am, data = mtcars) 
```

***It is larger!*** Why? How? Both differences should be expressed in units of
SD! But which SDs? Different SDs!

When looking at the difference between groups as a **slope**, the standardized
parameter is the difference between the means in $SD_{mpg}$. That is, the
*slope* between `Manual` and `Automatic` is a change of 1.20 $SD_{mpg}$s.

However, when looking a the difference as a **distance between two populations**, Cohen's d is the distance between the means in units of [**pooled SDs**](https://easystats.github.io/effectsize/reference/sd_pooled.html). That
is, the *distance* between `Manual` and `Automatic` is of 1.48 SDs of *each of
the groups* (here assumed to be equal).

In this simple model, the pooled SD is the residual SD, so we can also estimate
Cohen's *d* as:

```{r}
coef(m)[2] / sigma(m)
```

And we can also get an approximation of Cohen's *d* by converting the
$t$-statistic from the regression model via `t_to_d()`:

```{r}
parameters::model_parameters(m)

t_to_d(4.11, df_error = 30)
```

It is also interesting to note that using the `smart` method (explained in
detail below) when standardizing parameters will give you indices equivalent to
**Glass' *delta***, which is a standardized difference expressed in terms of SD
of the reference group.

```{r}
m <- lm(mpg ~ am, data = mtcars)

standardize_parameters(m, method = "smart")

glass_delta(mpg ~ am, data = mtcars)
```

***... So note that some standardized differences are different than others!
:)***

## Standardizing Parameters of Linear Models

As mentioned above, standardization of parameters can also be used to compare
among parameters within the same model. Essentially, what prevents us from
normally being able to compare among different parameters is that their
underlying variables are on different scales.[^But also as noted above, this is
not always an issue. For example, when the variables scale is important for the
interpretation of results, standardization might in fact hinder interpretation!]

For example, in the following example, we use a liner regression model to
predict a worker's salary (in Shmekels) from their age (years), seniority
(years), overtime (`xtra_hours`) and how many compliments they give their boss
(`n_comps`).

Let us explore the different parameter standardization methods provided by
`effectsize`.

### Standardized Slopes are Not (Always) Correlations

We saw that in simple linear models, the standardized slope is equal to the
correlation between the outcome and predictor - does this hold for **multiple
regression** as well? As in each effect in a regression model is "adjusted" for
the other ones, we might expect coefficients to be somewhat alike to **partial
correlations**. Let's first start by computing the partial correlation between
numeric predictors and the outcome.

```{r}
data("hardlyworking", package = "effectsize")

head(hardlyworking)

correlation::correlation(
  hardlyworking,
  select = "salary",
  select2 = c("xtra_hours", "n_comps", "age", "seniority"),
  partial = TRUE # get partial correlations
) 
```

Let's compare these to the standardized slopes:

```{r}
mod <- lm(salary ~ xtra_hours + n_comps + age + seniority,
          data = hardlyworking)

standardize_parameters(mod)
```

They are quite different! It seems then that ***standardized slopes in multiple
linear regressions are not the same a correlations or partial correlations***
:(

However, not all hope is lost yet - we can still try and recover the partial
correlations from our model, in another way: by converting the *t*-statistics
(and their degrees of freedom, *df*) into a partial correlation coefficient
*r*.

```{r}
params <- parameters::model_parameters(mod)

t_to_r(params$t[-1], df_error = params$df_error[-1])
```

Wow, the retrieved correlations coefficients from the regression model are
**exactly** the same as the partial correlations we estimated above! So these
"*r*" effect sizes can also be used.

### Methods of Standardizing Parameters

Let's convert `age` into a 3-level factor:

```{r}
hardlyworking$age_g <- cut(hardlyworking$age,
                           breaks = c(25,30,35,45))

mod <- lm(salary ~ xtra_hours + n_comps + age_g + seniority,
          data = hardlyworking)

parameters::model_parameters(mod)
```

It seems like the best or most important predictor is `n_comps` as it has the
coefficient. However, it is hard to compare among predictors, as they are on
different scales. To address this issue, we must have all the predictors on the
same scale - usually in the arbitrary unit of *standard deviations*.

#### **`"refit"`**: Re-fitting the model with standardized data

**This method is based on a complete model re-fit with a standardized version of
data**. Hence, this method is equal to standardizing the variables *before*
fitting the model. It is the "purest" and the most accurate [@neter1989applied],
but it is also the most computationally costly and long (especially for heavy
models such as Bayesian models, or complex mixed models). This method is
particularly recommended for models that include interactions or transformations
(e.g., exponentiation, log, polynomial or spline terms).

```{r}
standardize_parameters(mod, method = "refit")
```

`standardize_parameters` also has a `robust` argument (default to `FALSE`),
which enables a **robust standardization of the data**, *i.e.*, based on the
**median** and **MAD** instead of the **mean** and **SD**:

```{r}
standardize_parameters(mod, method = "refit", robust = TRUE)
```

Note that since `age_g` is a factor, it is not numerically standardized, and so
it standardized parameter is still not directly comparable to those of numeric
variables. To address this, we can set `two_sd = TRUE`, thereby scaling
parameters on 2 SDs (or MADs) of the predictors [@gelman2008scaling].

```{r}
standardize_parameters(mod, method = "refit", two_sd = TRUE)
```

`effectsize` also comes with a helper function that returns the re-fit model,
without summarizing it, which can then be used as the original model would:

```{r}
mod_z <- standardize(mod, two_sd = FALSE, robust = FALSE)
mod_z

parameters::model_parameters(mod_z)
```

#### **`"posthoc"`**: Refit without refitting

Post-hoc standardization of the parameters aims at emulating the results
obtained by `"refit"` without refitting the model. The coefficients are divided
by the standard deviation (or MAD if `robust`) of the outcome (which becomes
their expression 'unit'). Then, the coefficients related to numeric variables
are additionally multiplied by the standard deviation (or MAD if `robust`) of
the related terms, so that they correspond to changes of 1 SD of the predictor
(e.g., "A change in 1 SD of *x* is related to a change of 0.24 of the SD of
*y*). This does not apply to binary variables or factors, so the coefficients
are still related to changes in levels. This method is not accurate and tend to
give aberrant results when interactions are specified.

```{r}
standardize_parameters(mod, method = "posthoc")
```

#### **`"smart"`**: Standardization of Model's parameters with Adjustment, Reconnaissance and Transformation

> Experimental

Similar to `method = "posthoc"` in that it does not involve model refitting. The
difference is that the SD of the response is computed on the relevant section of
the data. For instance, if a factor with 3 levels A (the intercept), B and C is
entered as a predictor, the effect corresponding to B vs. A will be scaled by
the variance of the response at the intercept only. As a results, the
coefficients for effects of factors are similar to a Glass' *delta*.

```{r}
standardize_parameters(mod, method = "smart")
```

#### **`"basic"`**: Raw scaling of the model frame

This method is similar to `method = "posthoc"`, but treats all variables as
continuous: it scales the coefficient by the standard deviation of model's
matrix' parameter of factors levels (transformed to integers) or binary
predictors. Although it can be argued that this might be inappropriate for these
cases, this method allows for easier importance judgment across all predictor
type (numeric, factor, interactions...). It is also the type of standardization
implemented by default in other software packages (also `lm.beta::lm.beta()`),
and, such as can be used for reproducibility and replication purposes.

```{r}
standardize_parameters(mod, method = "basic")
```

### Standardizing Parameters In Mixed Models

Linear mixed models (LMM/HLM/MLM) offer an additional conundrum to
standardization - how does one even calculate the SDs of the various predictors?
Or of the response - is it the deviations within each group? Or perhaps between
them?

The solution: standardize according to level of the predictor
[@hoffman2015longitudinal, page 342]! Level 1 parameters are standardized
according to variance *within* groups, while level 2 parameters are standardized
according to variance *between* groups. The resulting standardized coefficient
are also called *pseudo*-standardized coefficients.[^Note that like method
`"basic"`, these are based on the model matrix.]

```{r, eval=knitr::opts_chunk$get("eval") && require(lme4) && require(lmerTest), warning=FALSE}
m <- lme4::lmer(Reaction ~ Days + (Days|Subject), data = lme4::sleepstudy)

standardize_parameters(m, method = "pseudo", ci_method = "satterthwaite")

# compare to:
standardize_parameters(m, method = "basic", ci_method = "satterthwaite")
```

### Standardizing Parameters In Generalized Linear Models

Unlike linear (/mixed) models, in generalized linear (/mixed) models (GLMs)
there is *less* of a need for standardization. Why? Because in many GLMs the
estimated coefficients are themselves measures of effect size, such as
*odds-ratios* (OR) in logistic regression, or *incidence rate ratios* (IRR) in
Poisson regressions. This is because in such model the outcome is **not** on an
arbitrary scale - that is, the meaning of rates and probabilities are changed by
arbitrary linear transformations.

But still, some standardization is sometimes needed, for the predictors.
Luckily, `standardize_parameters()` (and `standardize()`) are smart enough to
know when GLMs are passed so as to only standardize according to the
predictors:

```{r}
mod_b <- glm(am ~ mpg + factor(cyl),
             data = mtcars,
             family = binomial())

standardize_parameters(mod_b, method = "refit", two_sd = TRUE)
# standardize_parameters(mod_b, method = "posthoc", two_sd = TRUE)
# standardize_parameters(mod_b, method = "basic")
```

These can then be converted to OR (with `exp()`) and discussed as the "*change
in Odds as a function of a change in one SD of x*".

```{r}
std <- standardize_parameters(mod_b, method = "refit", two_sd = TRUE)
exp(std$Std_Coefficient)
```

Or we can directly ask for the coefficients to be exponentiated:

```{r}
standardize_parameters(mod_b, method = "refit", two_sd = TRUE, exponentiate = TRUE)
```

## Cohen's *f*

Cohen's $f$ (of [ANOVA fame](https://easystats.github.io/effectsize/articles/anovaES.html)) can be used as a measure of effect size in the context of sequential multiple regression (i.e., [**nested models**](https://easystats.github.io/performance/reference/test_performance.html)).
That is, when comparing two models, we can examine the ratio between the
increase in $R^2$ and the unexplained variance:

$$
f^{2}={R_{AB}^{2}-R_{A}^{2} \over 1-R_{AB}^{2}}
$$

```{r}
m1 <- lm(salary ~ xtra_hours, data = hardlyworking)
m2 <- lm(salary ~ xtra_hours + n_comps + seniority, data = hardlyworking)

cohens_f_squared(m1, model2 = m2)
```

<!-- ## Methods Comparison -->

<!-- We will use the "refit" method as the baseline. We will then compute the
differences between these standardized parameters and the ones provided by the
other functions. The **bigger the (absolute) number, the worse it is**. -->

<!-- > **SPOILER ALERT: the standardization implemented in `effectsize` is the
most accurate and the most flexible.** -->

<!-- ### Convenience function -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- library(effectsize) --> <!--
library(lm.beta) --> <!-- library(MuMIn) -->

<!-- comparison <- function(model, robust=FALSE){ --> <!-- out <-
standardize_parameters(model, method="refit", robust=robust)[1:2] -->

<!-- out$posthoc <- tryCatch({ --> <!-- out[, 2] - standardize_parameters(model,
method="posthoc", robust=robust)[, 2] --> <!-- }, error =
function(error_condition) { --> <!-- "Error" --> <!-- }) --> <!-- out$basic <-
tryCatch({ --> <!-- out[, 2] - standardize_parameters(model, method="basic",
robust=robust)[, 2] --> <!-- }, error = function(error_condition) { --> <!--
"Error" --> <!-- }) -->

<!-- out$lm.beta <- tryCatch({ --> <!-- out[, 2] -
lm.beta::lm.beta(model)$standardized.coefficients --> <!-- }, error =
function(error_condition) { --> <!-- "Error" --> <!-- }, warning =
function(warning_condition) { --> <!-- "Error" --> <!-- }) -->

<!-- out$MuMIn <- tryCatch({ --> <!-- out[, 2] - MuMIn::std.coef(model,
partial.sd=FALSE)[, 1] --> <!-- }, error = function(error_condition) { --> <!--
"Error" --> <!-- }) -->

<!-- ### Data -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- data <- iris --> <!--
data$Group_Sepal.Width <- as.factor(ifelse(data$Sepal.Width > 3, "High", "Low"))
--> <!-- data$Binary_Sepal.Width <- as.factor(ifelse(data$Sepal.Width > 3, 1,
0)) -->

<!-- summary(data) --> <!-- ``` -->

<!-- ### Models with only numeric predictors --> <!-- HERE -->

<!-- #### Linear Model -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <- lm(Sepal.Length ~
Petal.Width + Sepal.Width, data=data) --> <!-- comparison(model) --> <!-- ```
-->

<!-- #### Linear Mixed Model -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- library(lme4) -->

<!-- model <- lme4::lmer(Sepal.Length ~ Petal.Width + Sepal.Width + (1|Species),
--> <!-- data=data) --> <!-- comparison(model) --> <!-- ``` -->

<!-- #### Bayesian Models -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- library(rstanarm) -->

<!-- model <- stan_glm(Sepal.Length ~ Petal.Width + Sepal.Width, data=data) -->
<!-- comparison(model) --> <!-- ``` -->

<!-- For these simple models, **all methods return results equal to the "refit"
method** (although the other packages fail). -->

<!-- #### Transformation -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <- lm(Sepal.Length ~
poly(Petal.Width, 2) + poly(Sepal.Width, 2), data=data) --> <!--
comparison(model) --> <!-- ``` -->

<!-- When transformation are involved (e.g., polynomial transformations), **the
basic method becomes very unreliable**. -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <- lm(Sepal.Length ~
Petal.Width + Group_Sepal.Width, data=data) --> <!-- comparison(model) --> <!--
``` -->

<!-- #### Logistic Model -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <-
glm(Binary_Sepal.Width ~ Petal.Width + Species, data=data, family="binomial")
--> <!-- comparison(model) --> <!-- ``` -->

<!-- #### Linear Mixed Model -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- library(lme4) -->

<!-- model <- lme4::lmer(Sepal.Length ~ Petal.Length + Group_Sepal.Width +
(1|Species), data=data) --> <!-- comparison(model) --> <!-- ``` -->

<!-- #### Bayesian Models -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- library(rstanarm) -->

<!-- model <- stan_lmer(Sepal.Length ~ Petal.Width + Group_Sepal.Width +
(1|Species), --> <!-- data=data) --> <!-- comparison(model) --> <!-- ``` -->

<!-- When factors are involved, the basic method (that standardizes the numeric
transformation of factors) give again different results. --> <!-- HERE -->

<!-- ### Models with interactions -->

<!-- Long story short, coeffcient obtained via **posthoc** standardization
(without refitting the model) go berserk when interactions are involved.
However, **this is "normal"**: a regression model estimates coefficient between
two variables when the other predictors are at 0 (are *fixed* at 0, that people
interpret as *"adjusted for"*). When a standardized data is passed (in the
*refit* method), the effects and interactions are estimated at the **means** of
the other predictors (because 0 is the mean for a standardized variable).
Whereas in posthoc standardization, this coefficient correspond to something
different (because the 0 corresponds to something different in standardzed and
non-standardized data). In other words, when it comes to interaction, passing
standardized data results in a different model, which coefficient have an
intrinsically different meaning from unstandardized data. And as [for
now](https://github.com/easystats/effectsize/issues/6/), we are unable to
retrieve one from another. -->

<!-- #### Between continuous -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <- lm(Sepal.Length ~
Petal.Width * Sepal.Width, data=data) --> <!-- comparison(model) --> <!-- ```
-->

<!-- #### Between factors -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <- lm(Sepal.Length ~
Species * Group_Sepal.Width, data=data) --> <!-- comparison(model) --> <!-- ```
-->

<!-- #### Between factors and continuous -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <- lm(Sepal.Length ~
Petal.Width * Group_Sepal.Width, data=data) --> <!-- comparison(model) --> <!--
``` -->

<!-- ```{r message=FALSE, warning=FALSE} --> <!-- model <- lm(Sepal.Length ~
Group_Sepal.Width * Petal.Width, data=data) --> <!-- comparison(model) --> <!--
``` -->

<!-- ## Conclusion -->

<!-- Use `refit` if possible, but if no interactions, can use `posthoc` or
`smart`. -->

# References

---
title: "Effect Size from Test Statistics"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, standardization, cohen d]
vignette: >
  %\VignetteIndexEntry{Effect Size from Test Statistics}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
library(effectsize)

knitr::opts_chunk$set(comment = ">")
options(digits = 2)
options(knitr.kable.NA = "")

pkgs <- c("effectsize", "afex", "lmerTest", "emmeans", "parameters")
if (!all(sapply(pkgs, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(afex)
  library(lmerTest)
  library(emmeans)
  library(parameters)
}

set.seed(747)
```

# Introduction

In many real world applications there are no straightforward ways of obtaining
standardized effect sizes. However, it is possible to get approximations of most
of the effect size indices ($d$, $r$, $\eta^2_p$...) with the use of test
statistics. These conversions are based on the idea that **test statistics are a
function of effect size and sample size**. Thus information about samples size
(or more often of degrees of freedom) is used to reverse-engineer indices of
effect size from test statistics. This idea and these functions also power our
[***Effect Sizes From Test Statistics*** *shiny app*](https://easystats4u.shinyapps.io/statistic2effectsize/).

The measures discussed here are, in one way or another, ***signal to noise
ratios***, with the "noise" representing the unaccounted variance in the outcome
variable^[Note that for generalized linear models (Poisson, Logistic...), where
the outcome is never on an arbitrary scale, estimates themselves **are** indices
of effect size! Thus this vignette is relevant only to general linear models.].

The indices are:

- Percent variance explained ($\eta^2_p$, $\omega^2_p$, $\epsilon^2_p$).

- Measure of association ($r$).

- Measure of difference ($d$).

## (Partial) Percent Variance Explained

These measures represent the ratio of $Signal^2 / (Signal^2 + Noise^2)$, with
the "noise" having all other "signals" partial-ed out (be they of other fixed or
random effects). The most popular of these indices is $\eta^2_p$ (Eta; which is
equivalent to $R^2$).

The conversion of the $F$- or $t$-statistic is based on
@friedman1982simplified.

Let's look at an example:

```{r}
library(afex)

data(md_12.1)

aov_fit <- aov_car(rt ~ angle * noise + Error(id / (angle * noise)),
  data = md_12.1,
  anova_table = list(correction = "none", es = "pes")
)
aov_fit
```

Let's compare the $\eta^2_p$ (the `pes` column) obtained here with ones
recovered from `F_to_eta2()`:

```{r}
library(effectsize)

F_to_eta2(
  f = c(40.72, 33.77, 45.31),
  df = c(2, 1, 2),
  df_error = c(18, 9, 18)
)
```

**They are identical!**^[Note that these are *partial* percent variance
explained, and so their sum can be larger than 1.] (except for the fact that
`F_to_eta2()` also provides confidence intervals^[Confidence intervals for all
indices are estimated using the non-centrality parameter method; These methods
search for a the best non-central parameter of the non-central $F$/$t$
distribution for the desired tail-probabilities, and then convert these ncps to
the corresponding effect sizes.] :)

In this case we were able to easily obtain the effect size (thanks to `afex`!),
but in other cases it might not be as easy, and using estimates based on test
statistic offers a good approximation.

For example:

### In Simple Effect and Contrast Analysis

```{r}
library(emmeans)

joint_tests(aov_fit, by = "noise")

F_to_eta2(
  f = c(5, 79),
  df = 2,
  df_error = 29
)
```

We can also use `t_to_eta2()` for contrast analysis:

```{r}
pairs(emmeans(aov_fit, ~angle))

t_to_eta2(
  t = c(-5.7, -8.9, -3.2),
  df_error = 18
)
```

### In Linear Mixed Models

```{r}
library(lmerTest)

fit_lmm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

anova(fit_lmm)

F_to_eta2(45.8, 1, 17)
```

We can also use `t_to_eta2()` for the slope of `Days` (which in this case gives
the same result).

```{r}
parameters::model_parameters(fit_lmm, effects = "fixed", ci_method = "satterthwaite")

t_to_eta2(6.77, df_error = 17)
```

### Bias-Corrected Indices

Alongside $\eta^2_p$ there are also the less biased $\omega_p^2$ (Omega) and
$\epsilon^2_p$ (Epsilon; sometimes called $\text{Adj. }\eta^2_p$, which is
equivalent to $R^2_{adj}$; @albers2018power, @mordkoff2019simple).

```{r}
F_to_eta2(45.8, 1, 17)
F_to_epsilon2(45.8, 1, 17)
F_to_omega2(45.8, 1, 17)
```

## Measure of Association

Similar to $\eta^2_p$, $r$ is a signal to noise ratio, and is in fact equal to
$\sqrt{\eta^2_p}$ (so it's really a *partial* $r$). It is often used instead of
$\eta^2_p$ when discussing the *strength* of association (but I suspect people
use it instead of $\eta^2_p$ because it gives a bigger number, which looks
better).

### For Slopes

```{r}
parameters::model_parameters(fit_lmm, effects = "fixed", ci_method = "satterthwaite")

t_to_r(6.77, df_error = 17)
```

In a fixed-effect linear model, this returns the **partial** correlation.
Compare:

```{r}

fit_lm <- lm(rating ~ complaints + critical, data = attitude)

parameters::model_parameters(fit_lm)

t_to_r(
  t = c(7.46, 0.01),
  df_error = 27
)
```

to:

```{r, eval=require(correlation, quietly = TRUE)}
correlation::correlation(attitude[, c(1, 2, 6)], partial = TRUE)[1:2, c(2, 3, 7, 8)]
```

### In Contrast Analysis

This measure is also sometimes used in contrast analysis, where it is called the
point bi-serial correlation - $r_{pb}$ [@cohen1965some; @rosnow2000contrasts]:

```{r}
pairs(emmeans(aov_fit, ~angle))

t_to_r(
  t = c(-5.7, -8.9, -3.2),
  df_error = 18
)
```

## Measures of Difference

These indices represent $Signal/Noise$ with the "signal" representing the
difference between two means. This is akin to Cohen's $d$, and is a close
approximation when comparing two groups of equal size [@wolf1986meta;
@rosnow2000contrasts].

These can be useful in contrast analyses.

### Between-Subject Contrasts

```{r}
m <- lm(breaks ~ tension, data = warpbreaks)

em_tension <- emmeans(m, ~tension)
pairs(em_tension)

t_to_d(
  t = c(2.53, 3.72, 1.20),
  df_error = 51
)
```

However, these are merely approximations of a *true* Cohen's *d*. It is advised
to directly estimate Cohen's *d*, whenever possible. For example, here with
`emmeans::eff_size()`:

```{r}
eff_size(em_tension, sigma = sigma(m), edf = df.residual(m))
```

### Within-Subject Contrasts

```{r}

pairs(emmeans(aov_fit, ~angle))

t_to_d(
  t = c(-5.7, -5.9, -3.2),
  df_error = 18,
  paired = TRUE
)
```

(Note set `paired = TRUE` to not over estimate the size of the effect;
@rosenthal1991meta; @rosnow2000contrasts)

# References
---
title: "Converting Between Indices of Effect Size"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, rules of thumb, guidelines, conversion]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{Converting Between Indices of Effect Size}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
knitr::opts_chunk$set(comment = ">")
options(digits = 3)

pkgs <- c("effectsize", "ggplot2", "correlation", "parameters", "bayestestR")
if (!all(sapply(pkgs, require, quietly = TRUE, character.only = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}
```

The `effectsize` package contains function to convert among indices of effect
size. This can be useful for meta-analyses, or any comparison between different
types of statistical analyses.

## Converting Between *d*, *r*, and *OR*

The most basic conversion is between *r* values, a measure of standardized
association between two continuous measures, and *d* values (such as Cohen's
*d*), a measure of standardized differences between two groups / conditions.

Let's simulate some data:

```{r}
set.seed(1)
data <- bayestestR::simulate_difference(
  n = 10,
  d = 0.2,
  names = c("Group", "Outcome")
)
```

```{r, echo=FALSE}
print(data, digits = 3)
```

We can compute Cohen's *d* between the two groups:

```{r}
cohens_d(Outcome ~ Group, data = data)
```

But we can also treat the 2-level `group` variable as a numeric variable, and
compute Pearson's *r*:

```{r, warning=FALSE}
correlation::correlation(data, include_factors = TRUE)[2, ]
```

But what if we only have summary statistics? Say, we only have $d=-0.37$ and we
want to know what the *r* would have been? We can approximate *r* using the
following formula [@borenstein2009converting]:

$$
r \approx \frac{d}{\sqrt{d^2 + 4}}
$$
And indeed, if we use `d_to_r()`, we get a pretty decent approximation:

```{r}
d_to_r(-0.31)
```

(Which also works in the other way, with `r_to_d(0.17)` gives `r
round(r_to_d(-0.17),3)`)

As we can see, these are rough approximations, but they can be useful when we
don't have the raw data on hand.

### In multiple regression

Although not exactly a classic Cohen's d, we can also approximate a partial-*d*
value (that is, the standardized difference between two groups / conditions,
with variance from other predictors partilled out). For example:

```{r}
fit <- lm(mpg ~ am + hp, data = mtcars)

parameters::model_parameters(fit)

# A couple of ways to get partial-d:
5.28 / sigma(fit)
t_to_d(4.89, df_error = 29)[[1]]
```

We can convert these semi-*d* values to *r* values, but in this case these
represent the *partial* correlation:

```{r}
t_to_r(4.89, df_error = 29)

correlation::correlation(mtcars[, c("mpg", "am", "hp")], partial = TRUE)[1, ]

# all close to:
d_to_r(1.81)
```

## From Odds ratios

In binomial regression (more specifically in logistic regression), Odds ratios
(OR) are themselves measures of effect size; they indicate the expected change
in the odds of a some event.

In some fields, it is common to dichotomize outcomes in order to be able to
analyze them with logistic models. For example, if the outcome is the count of
white blood cells, it can be more useful (medically) to predict the crossing of
the threshold rather than the raw count itself. And so, where some scientists
would maybe analyze the above data with a *t*-test and present Cohen's *d*,
others might analyze it with a logistic regression model on the dichotomized
outcome, and present OR. So the question can be asked: given such a OR, what
would Cohen's *d* have been?

Fortunately, there is a formula to approximate this [@sanchez2003effect]:

$$
d = log(OR) \times \frac{\sqrt{3}}{\pi}
$$

which is implemented in the `oddsratio_to_d()` function.

Let's give it a try:

```{r}
# 1. Set a threshold
thresh <- 0

# 2. dichotomize the outcome
data$Outcome_binom <- data$Outcome < thresh

# 3. Fit a logistic regression:
fit <- glm(Outcome_binom ~ Group,
  data = data,
  family = binomial()
)

parameters::model_parameters(fit)

# Convert log(OR) (the coefficient) to d
oddsratio_to_d(-0.81, log = TRUE)
```

### Odds ratios to Risk Ratios

Odds ratio, although popular, are not very intuitive in their interpretations.
We don't often think about the chances of catching a disease in terms of *odds*,
instead we instead tend to think in terms of *probability* or some event - or
the *risk*. Talking about *risks* we can also talk about the *change in risk*,
knows as the *risk ratio* (*RR*).

For example, if we find that for individual suffering from a migraine, for every
bowl of brussels sprouts they eat, they're odds of reducing the migraine
increase by an $OR = 3.5$ over a period of an hour. So, should people eat
brussels sprouts to effectively reduce pain? Well, hard to say... Maybe if we
look at *RR* we'll get a clue.

We can convert between *OR* and *RR* for the following formula
[@grant2014converting]:

$$
RR = \frac{OR}{(1 - p0 + (p0 \times OR))}  
$$

Where $p0$ is the base-rate risk - the probability of the event without the
intervention (e.g., what is the probability of the migraine subsiding within an
hour without eating any brussels sprouts). If it the base-rate risk is, say,
85%, we get a *RR* of:

```{r}
OR <- 3.5
baserate <- 0.85

oddsratio_to_riskratio(OR, baserate)
```

That is - for every bowl of brussels sprouts, we increase the chances of
reducing the migraine by a mere 12%! Is if worth it? Depends on you affinity to
brussels sprouts...

Note that the base-rate risk is crucial here. If instead of 85% it was only 4%,
then the *RR* would be:

```{r}
OR <- 3.5
baserate <- 0.04

oddsratio_to_riskratio(OR, baserate)
```

That is - for every bowl of brussels sprouts, we increase the chances of
reducing the migraine by a whopping 318%! Is if worth it? I guess that still
depends on your affinity to brussels sprouts...

# References
---
title: "Automated Interpretation of Indices of Effect Size"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, effect size, rules of thumb, guidelines, interpretation]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{Automated Interpretation of Indices of Effect Size}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
knitr::opts_chunk$set(comment = ">")
options(digits = 2)

pkgs <- c("effectsize")
if (!all(sapply(pkgs, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}
```

## Why?

The metrics used in statistics (indices of fit, model performance, or parameter
estimates) can be very abstract. A long experience is required to intuitively
***feel*** the meaning of their values. In order to facilitate the understanding
of the results they are facing, many scientists use (often implicitly) some set
of **rules of thumb**. Some of these rules of thumb have been standardize and
validated and subsequently published as guidelines. Understandably then, such
rules of thumb are just suggestions and there is nothing universal about them.
The interpretation of **any** effect size measures is always going to be
relative to the discipline, the specific data, and the aims of the analyst. This
is important because what might be considered a small effect in psychology might
be large for some other field like public health.

One of the most famous interpretation grids was proposed by **Cohen (1988)** for
a series of widely used indices, such as the correlation **r** (*r* = .20,
small; *r* = .40, moderate and *r* = .60, large) or the **standardized difference** (*Cohen's d*). However, there is now a clear evidence that Cohen's
guidelines (which he himself later disavowed; Funder, 2019) are much too
stringent and not particularly meaningful taken out of context
[@funder2019evaluating]. This led to the emergence of a literature discussing
and creating new sets of rules of thumb.

Although **everybody** agrees on the fact that effect size interpretation in a
study should be justified with a rationale (and depend on the context, the
field, the literature, the hypothesis, etc.), these pre-baked rules can
nevertheless be useful to give a rough idea or frame of reference to understand
scientific results.

The package **`effectsize`** catalogs such sets of rules of thumb for a
variety of indices in a flexible and explicit fashion, helping you understand
and report your results in a scientific yet meaningful way. Again, readers
should keep in mind that these thresholds, as ubiquitous as they may be,
**remain arbitrary**. Thus, their use should be discussed on a case-by-case
basis depending on the field, hypotheses, prior results, and so on, to avoid
their crystallization, as for the infamous $p < .05$ criterion of hypothesis
testing.

Moreover, some authors suggest the counter-intuitive idea that *very large effects*, especially in the context of psychological research, is likely to be a
"gross overestimate that will rarely be found in a large sample or in a
replication" [@funder2019evaluating]. They suggest that smaller effect size are
worth taking seriously (as they can be potentially consequential), as well as
more believable.

## Correlation *r*

#### @funder2019evaluating

```r
interpret_r(x, rules = "funder2019")
```

- **r < 0.05** - Tiny

- **0.05 <= r < 0.1** - Very small

- **0.1 <= r < 0.2** - Small

- **0.2 <= r < 0.3** - Medium

- **0.3 <= r < 0.4** - Large

- **r >= 0.4** - Very large

#### @gignac2016effect

Gignac's rules of thumb are actually one of few interpretation grid justified
and based on actual data, in this case on the distribution of effect magnitudes
in the literature.

```r
interpret_r(x, rules = "gignac2016")
```

- **r < 0.1** - Very small

- **0.1 <= r < 0.2** - Small

- **0.2 <= r < 0.3** - Moderate

- **r >= 0.3** - Large

#### @cohen1988statistical

```r
interpret_r(x, rules = "cohen1988")
```

- **r < 0.1** - Very small

- **0.1 <= r < 0.3** - Small

- **0.3 <= r < 0.5** - Moderate

- **r >= 0.5** - Large

#### @evans1996straightforward

```r
interpret_r(x, rules = "evans1996")
```

- **r < 0.2** - Very weak

- **0.2 <= r < 0.4** - Weak

- **0.4 <= r < 0.6** - Moderate

- **0.6 <= r < 0.8** - Strong

- **r >= 0.8** - Very strong

#### @lovakov2021empirically

```r
interpret_r(x, rules = "lovakov2021")
```

- **r < 0.12** - Very small

- **0.12 <= r < 0.24** - Small

- **0.24 <= r < 0.41** - Moderate

- **r >= 0.41** - Large


## Standardized Difference *d* (Cohen's *d*)

The standardized difference can be obtained through the standardization of
linear model's parameters or data, in which they can be used as indices of
effect size.

#### @cohen1988statistical

```r
interpret_cohens_d(x, rules = "cohen1988")
```

- **d < 0.2** - Very small

- **0.2 <= d < 0.5** - Small

- **0.5 <= d < 0.8** - Medium

- **d >= 0.8** - Large

#### @sawilowsky2009new

```r
interpret_cohens_d(x, rules = "sawilowsky2009")
```

- **d < 0.1** - Tiny

- **0.1 <= d < 0.2** - Very small

- **0.2 <= d < 0.5** - Small

- **0.5 <= d < 0.8** - Medium

- **0.8 <= d < 1.2** - Large

- **1.2 <= d < 2** - Very large

- **d >= 2** - Huge

#### @gignac2016effect

Gignac's rules of thumb are actually one of few interpretation grid justified
and based on actual data, in this case on the distribution of effect magnitudes
in the literature. These is in fact the same grid used for *r*, based on the
conversion of *r* to *d*:

```r
interpret_cohens_d(x, rules = "gignac2016")
```

- **d < 0.2** - Very small

- **0.2 <= d < 0.41** - Small

- **0.41 <= d < 0.63** - Moderate

- **d >= 0.63** - Large

#### @lovakov2021empirically

```r
interpret_cohens_d(x, rules = "lovakov2021")
```

- **r < 0.15** - Very small

- **0.15 <= r < 0.36** - Small

- **0.36 <= r < 0.65** - Moderate

- **r >= 0.65** - Large


## Odds Ratio (OR)

Odds ratio, and *log* odds ratio, are often found in epidemiological studies.
However, they are also the parameters of ***logistic*** regressions, where they
can be used as indices of effect size. Note that the (log) odds ratio from
logistic regression coefficients are *unstandardized*, as they depend on the
scale of the predictor. In order to apply the following guidelines, make sure
you
[*standardize*](https://easystats.github.io/effectsize/articles/standardize_parameters.html)
your predictors!

Keep in mind that these apply to Odds *ratios*, so Odds ratio of 10 is as
extreme as a Odds ratio of 0.1 (1/10).

#### @chen2010big

```r
interpret_oddsratio(x, rules = "chen2010")
```

- **OR < 1.68** - Very small

- **1.68 <= OR < 3.47** - Small

- **3.47 <= OR < 6.71** - Medium

- **OR >= 6.71 ** - Large

#### @cohen1988statistical

```r
interpret_oddsratio(x, rules = "cohen1988")
```

- **OR < 1.44** - Very small

- **1.44 <= OR < 2.48** - Small

- **2.48 <= OR < 4.27** - Medium

- **OR >= 4.27 ** - Large

This converts (log) odds ratio to standardized difference *d* using the
following formula [@cohen1988statistical;@sanchez2003effect]:

$$
d = log(OR) \times \frac{\sqrt{3}}{\pi}
$$

## Coefficient of determination  (R<sup>2</sup>)

### For Linear Regression

#### @cohen1988statistical

```r
interpret_r2(x, rules = "cohen1988")
```

- **R2 < 0.02** - Very weak

- **0.02 <= R2 < 0.13** - Weak

- **0.13 <= R2 < 0.26** - Moderate

- **R2 >= 0.26** - Substantial

#### @falk1992primer

```r
interpret_r2(x, rules = "falk1992")
```

- **R2 < 0.1** - Negligible

- **R2 >= 0.1** - Adequate

### For PLS / SEM R-Squared of *latent* variables

#### @chin1998partial

```r
interpret_r2(x, rules = "chin1998")
```

- **R2 < 0.19** - Very weak

- **0.19 <= R2 < 0.33** - Weak

- **0.33 <= R2 < 0.67** - Moderate

- **R2 >= 0.67** - Substantial

#### @hair2011pls

```r
interpret_r2(x, rules = "hair2011")
```

- **R2 < 0.25** - Very weak

- **0.25 <= R2 < 0.50** - Weak

- **0.50 <= R2 < 0.75** - Moderate

- **R2 >= 0.75** - Substantial

## Omega / Eta / Epsilon Squared

The Omega squared is a measure of effect size used in ANOVAs. It is an estimate
of how much variance in the response variables are accounted for by the
explanatory variables. Omega squared is widely viewed as a lesser biased
alternative to eta-squared, especially when sample sizes are small.

#### @field2013discovering

```r
interpret_omega_squared(x, rules = "field2013")
```

- **ES < 0.01** - Very small

- **0.01 <= ES < 0.06** - Small

- **0.16 <= ES < 0.14** - Medium

- **ES >= 0.14 ** - Large

#### @cohen1992power

These are applicable to one-way ANOVAs, or to *partial* Eta / Omega / Epsilon
Squared in a multi-way ANOVA.

```r
interpret_omega_squared(x, rules = "cohen1992")
```

- **ES < 0.02** - Very small

- **0.02 <= ES < 0.13** - Small

- **0.13 <= ES < 0.26** - Medium

- **ES >= 0.26** - Large


##  Kendall's coefficient of concordance

The interpretation of Kendall's coefficient of concordance (*w*) is a measure of
effect size used in non-parametric ANOVAs (the Friedman rank sum test). It is an
estimate of agreement among multiple raters.

#### @landis1977measurement

```r
interpret_omega_squared(w, rules = "landis1977")
```

- **0.00 <= w < 0.20** - Slight agreement
- **0.20 <= w < 0.40** - Fair agreement
- **0.40 <= w < 0.60** - Moderate agreement
- **0.60 <= w < 0.80** - Substantial agreement
- **w >= 0.80**        - Almost perfect agreement

## Cohen's *g*

Cohen's *g* is a measure of effect size used for McNemar's test of agreement in
selection - when repeating a multiple chose selection, is the percent of matches
(first response is equal to the second response) different than 50%?

#### @cohen1988statistical

```r
interpret_cohens_g(x, rules = "cohen1988")
```

- **d < 0.05** - Very small

- **0.05 <= d < 0.15** - Small

- **0.15 <= d < 0.25** - Medium

- **d >= 0.25** - Large

## Interpretation of other Indices

`effectsize` also offers functions for interpreting other statistical indices:

- `interpret_gfi()`, `interpret_agfi()`, `interpret_nfi()`, `interpret_nnfi()`,
  `interpret_cfi()`, `interpret_rmsea()`, `interpret_srmr()`, `interpret_rfi()`,
  `interpret_ifi()`, and `interpret_pnfi()` for interpretation CFA / SEM
  goodness of fit.

- `interpret_p()` for interpretation of *p*-values.

- `interpret_direction()` for interpretation of direction.

- `interpret_bf()` for interpretation of Bayes factors.

- `interpret_rope()` for interpretation of Bayesian ROPE tests.

- `interpret_ess()` and `interpret_rhat()` for interpretation of Bayesian
  diagnostic indices.

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common_language.R
\name{cles}
\alias{cles}
\alias{common_language}
\alias{cohens_u3}
\alias{p_superiority}
\alias{p_overlap}
\title{Estimate Common Language Effect Sizes (CLES)}
\usage{
cles(
  x,
  y = NULL,
  data = NULL,
  mu = 0,
  ci = 0.95,
  alternative = "two.sided",
  parametric = TRUE,
  verbose = TRUE,
  iterations = 200,
  ...
)

common_language(
  x,
  y = NULL,
  data = NULL,
  mu = 0,
  ci = 0.95,
  alternative = "two.sided",
  parametric = TRUE,
  verbose = TRUE,
  iterations = 200,
  ...
)

cohens_u3(...)

p_superiority(...)

p_overlap(...)
}
\arguments{
\item{x}{A formula, a numeric vector, or a character name of one in \code{data}.}

\item{y}{A numeric vector, a grouping (character / factor) vector, a or a
character  name of one in \code{data}. Ignored if \code{x} is a formula.}

\item{data}{An optional data frame containing the variables.}

\item{mu}{a number indicating the true value of the mean (or
    difference in means if you are performing a two sample test).}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"two.sided"} (default, two-sided CI),
\code{"greater"} or \code{"less"} (one-sided CI). Partial matching is allowed (e.g.,
\code{"g"}, \code{"l"}, \code{"two"}...). See \emph{One-Sided CIs} in \link{effectsize_CIs}.}

\item{parametric}{Use parametric estimation (see \code{\link[=cohens_d]{cohens_d()}}) or
non-parametric estimation (see \code{\link[=rank_biserial]{rank_biserial()}}).}

\item{verbose}{Toggle warnings and messages on or off.}

\item{iterations}{The number of bootstrap replicates for computing confidence
intervals. Only applies when \code{ci} is not \code{NULL} and \code{parametric = FALSE}.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame containing the common language effect sizes (and
optionally their CIs).
}
\description{
\code{cohens_u3()}, \code{p_superiority()}, and \code{p_overlap()} give only one of the
CLESs.
}
\details{
These measures of effect size present group differences in probabilistic
terms:
\itemize{
\item \strong{Probability of superiority} is the probability that, when sampling an
observation from each of the groups at random, that the observation from
the second group will be larger than the sample from the first group.
\item \strong{Cohen's U3} is the proportion of the second group that is smaller than
the median of the first group.
\item \strong{Overlap} (OVL) is the proportional overlap between the distributions.
(When \code{parametric = FALSE}, \code{\link[bayestestR:overlap]{bayestestR::overlap()}} is used.)
}

For unequal group sizes, it is recommended to use the non-parametric based
CLES (\code{parametric = FALSE}).
}
\section{Confidence Intervals (CIs)}{

For parametric CLES, the CIs are transformed CIs for Cohen's \emph{d}
(\code{\link[=d_to_cles]{d_to_cles()}}). For non-parametric (\code{parametric = FALSE}) CLES, the CI of
\emph{Pr(superiority)} is a transformed CI of the rank-biserial correlation
(\code{\link[=rb_to_cles]{rb_to_cles()}}), while for Cohen's \emph{U3} and the Overlap coefficient the
confidence intervals are bootstrapped (requires the \code{boot} package).
}

\examples{
cles(mpg ~ am, data = mtcars)

set.seed(4)
cles(mpg ~ am, data = mtcars, parametric = FALSE)

\dontrun{
## Individual CLES
p_superiority(extra ~ group, data = sleep)

cohens_u3(extra ~ group, data = sleep, parametric = FALSE)

p_overlap(extra ~ group, data = sleep)
}

}
\references{
\itemize{
\item Cohen, J. (1977). Statistical power analysis for the behavioral sciences.
New York: Routledge.
\item Reiser, B., & Faraggi, D. (1999). Confidence intervals for the overlapping
coefficient: the normal equal variance case. Journal of the Royal Statistical
Society, 48(3), 413-418.
\item Ruscio, J. (2008). A probability-based measure of effect size: robustness
to base rates and other factors. Psychological methods, 13(1), 19–30.
}
}
\seealso{
\code{\link[=d_to_cles]{d_to_cles()}} \code{\link[=sd_pooled]{sd_pooled()}}

Other effect size indices: 
\code{\link{cohens_d}()},
\code{\link{effectsize.BFBayesFactor}()},
\code{\link{eta_squared}()},
\code{\link{phi}()},
\code{\link{rank_biserial}()},
\code{\link{standardize_parameters}()}
}
\concept{effect size indices}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_between_anova.R
\name{eta2_to_f2}
\alias{eta2_to_f2}
\alias{eta2_to_f}
\alias{f2_to_eta2}
\alias{f_to_eta2}
\title{Convert between ANOVA effect sizes}
\usage{
eta2_to_f2(es)

eta2_to_f(es)

f2_to_eta2(f2)

f_to_eta2(f)
}
\arguments{
\item{es}{Any measure of variance explained such as Eta-, Epsilon-, Omega-,
or R-Squared, partial or otherwise. See details.}

\item{f, f2}{Cohen's \emph{f} or \emph{f}-squared.}
}
\description{
Convert between ANOVA effect sizes
}
\details{
Any measure of variance explained can be converted to a corresponding Cohen's
\emph{f} via:
\cr\cr
\deqn{f^2 = \frac{\eta^2}{1 - \eta^2}}
\cr\cr
\deqn{\eta^2 = \frac{f^2}{1 + f^2}}
\cr\cr
If a partial Eta-Squared is used, the resulting Cohen's \emph{f} is a
partial-Cohen's \emph{f}; If a less biased estimate of variance explained is used
(such as Epsilon- or Omega-Squared), the resulting Cohen's \emph{f} is likewise a
less biased estimate of Cohen's \emph{f}.
}
\references{
\itemize{
\item Cohen, J. (1988). Statistical power analysis for the behavioral sciences
(2nd Ed.). New York: Routledge.
\item Steiger, J. H. (2004). Beyond the F test: Effect size confidence intervals
and tests of close fit in the analysis of variance and contrast analysis.
Psychological Methods, 9, 164-182.
}
}
\seealso{
\code{\link[=eta_squared]{eta_squared()}} for more details.

Other convert between effect sizes: 
\code{\link{d_to_cles}()},
\code{\link{d_to_r}()},
\code{\link{odds_to_probs}()},
\code{\link{oddsratio_to_riskratio}()}
}
\concept{convert between effect sizes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize_info.R
\name{standardize_info}
\alias{standardize_info}
\title{Get Standardization Information}
\usage{
standardize_info(
  model,
  robust = FALSE,
  two_sd = FALSE,
  include_pseudo = FALSE,
  ...
)
}
\arguments{
\item{model}{A statistical model.}

\item{robust}{Logical, if \code{TRUE}, centering is done by subtracting the
median from the variables and dividing it by the median absolute deviation
(MAD). If \code{FALSE}, variables are standardized by subtracting the
mean and dividing it by the standard deviation (SD).}

\item{two_sd}{If \code{TRUE}, the variables are scaled by two times the deviation
(SD or MAD depending on \code{robust}). This method can be useful to obtain
model coefficients of continuous parameters comparable to coefficients
related to binary predictors, when applied to \strong{the predictors} (not the
outcome) (Gelman, 2008).}

\item{include_pseudo}{(For (G)LMMs) Should Pseudo-standardized information be
included?}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame with information on each parameter (see
\link[parameters:parameters_type]{parameters::parameters_type}), and various standardization coefficients
for the post-hoc methods (see \code{\link[=standardize_parameters]{standardize_parameters()}}) for the predictor
and the response.
}
\description{
This function extracts information, such as the deviations (SD or MAD) from
parent variables, that are necessary for post-hoc standardization of
parameters. This function gives a window on how standardized are obtained,
i.e., by what they are divided. The "basic" method of standardization uses.
}
\examples{
model <- lm(mpg ~ ., data = mtcars)
standardize_info(model)
standardize_info(model, robust = TRUE)
standardize_info(model, two_sd = TRUE)
}
\seealso{
Other standardize: 
\code{\link{standardize.default}()},
\code{\link{standardize_parameters}()}
}
\concept{standardize}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_r.R
\name{interpret_r}
\alias{interpret_r}
\alias{interpret_phi}
\alias{interpret_cramers_v}
\alias{interpret_rank_biserial}
\title{Interpret correlation coefficient}
\usage{
interpret_r(r, rules = "funder2019")

interpret_phi(r, rules = "funder2019")

interpret_cramers_v(r, rules = "funder2019")

interpret_rank_biserial(r, rules = "funder2019")
}
\arguments{
\item{r}{Value or vector of correlation coefficient.}

\item{rules}{Can be \code{"funder2019"} (default), \code{"gignac2016"}, \code{"cohen1988"},
\code{"evans1996"}, \code{"lovakov2021"} or a custom set of \code{\link[=rules]{rules()}}.}
}
\description{
Interpret correlation coefficient
}
\note{
As \eqn{\phi}{\phi} can be larger than 1 - it is recommended to compute
and interpret Cramer's \emph{V} instead.
}
\section{Rules}{


Rules apply positive and negative \emph{r} alike.
\itemize{
\item Funder & Ozer (2019) (\code{"funder2019"}; default)
\itemize{
\item \strong{r < 0.05} - Tiny
\item \strong{0.05 <= r < 0.1} - Very small
\item \strong{0.1 <= r < 0.2} - Small
\item \strong{0.2 <= r < 0.3} - Medium
\item \strong{0.3 <= r < 0.4} - Large
\item \strong{r >= 0.4} - Very large
}
\item Gignac & Szodorai (2016) (\code{"gignac2016"})
\itemize{
\item \strong{r < 0.1} - Very small
\item \strong{0.1 <= r < 0.2} - Small
\item \strong{0.2 <= r < 0.3} - Moderate
\item \strong{r >= 0.3} - Large
}
\item Cohen (1988) (\code{"cohen1988"})
\itemize{
\item \strong{r < 0.1} - Very small
\item \strong{0.1 <= r < 0.3} - Small
\item \strong{0.3 <= r < 0.5} - Moderate
\item \strong{r >= 0.5} - Large
}
\item Lovakov & Agadullina (2021) (\code{"lovakov2021"})
\itemize{
\item \strong{r < 0.12} - Very small
\item \strong{0.12 <= r < 0.24} - Small
\item \strong{0.24 <= r < 0.41} - Moderate
\item \strong{r >= 0.41} - Large
}
\item Evans (1996) (\code{"evans1996"})
\itemize{
\item \strong{r < 0.2} - Very weak
\item \strong{0.2 <= r < 0.4} - Weak
\item \strong{0.4 <= r < 0.6} - Moderate
\item \strong{0.6 <= r < 0.8} - Strong
\item \strong{r >= 0.8} - Very strong
}
}
}

\examples{
interpret_r(.015)
interpret_r(c(.5, -.02))
interpret_r(.3, rules = "lovakov2021")
}
\references{
\itemize{
\item Lovakov, A., & Agadullina, E. R. (2021). Empirically Derived Guidelines for
Effect Size Interpretation in Social Psychology. European Journal of Social
Psychology.
\item Funder, D. C., & Ozer, D. J. (2019). Evaluating effect size in
psychological research: sense and nonsense. Advances in Methods and Practices
in Psychological Science.
\item Gignac, G. E., & Szodorai, E. T. (2016). Effect size guidelines for
individual differences researchers. Personality and individual differences,
102, 74-78.
\item Cohen, J. (1988). Statistical power analysis for the behavioral sciences
(2nd Ed.). New York: Routledge.
\item Evans, J. D. (1996). Straightforward statistics for the behavioral
sciences. Thomson Brooks/Cole Publishing Co.
}
}
\seealso{
Page 88 of APA's 6th Edition.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_cfa_fit.R
\name{interpret_gfi}
\alias{interpret_gfi}
\alias{interpret_agfi}
\alias{interpret_nfi}
\alias{interpret_nnfi}
\alias{interpret_cfi}
\alias{interpret_rmsea}
\alias{interpret_srmr}
\alias{interpret_rfi}
\alias{interpret_ifi}
\alias{interpret_pnfi}
\alias{interpret.lavaan}
\alias{interpret.performance_lavaan}
\title{Interpret of indices of CFA / SEM goodness of fit}
\usage{
interpret_gfi(x, rules = "default")

interpret_agfi(x, rules = "default")

interpret_nfi(x, rules = "byrne1994")

interpret_nnfi(x, rules = "byrne1994")

interpret_cfi(x, rules = "default")

interpret_rmsea(x, rules = "default")

interpret_srmr(x, rules = "default")

interpret_rfi(x, rules = "default")

interpret_ifi(x, rules = "default")

interpret_pnfi(x, rules = "default")

\method{interpret}{lavaan}(x, ...)

\method{interpret}{performance_lavaan}(x, ...)
}
\arguments{
\item{x}{vector of values, or an object of class \code{lavaan}.}

\item{rules}{Can be \code{"default"} or custom set of \code{\link[=rules]{rules()}}.}

\item{...}{Currently not used.}
}
\description{
Interpretation of indices of fit found in confirmatory analysis or structural
equation modelling, such as RMSEA, CFI, NFI, IFI, etc.
}
\details{
\subsection{Indices of fit}{
\itemize{
\item \strong{Chisq}: The model Chi-squared assesses overall fit and the discrepancy
between the sample and fitted covariance matrices. Its p-value should be >
.05 (i.e., the hypothesis of a perfect fit cannot be rejected). However, it
is quite sensitive to sample size.
\item \strong{GFI/AGFI}: The (Adjusted) Goodness of Fit is the proportion of variance
accounted for by the estimated population covariance. Analogous to R2. The
GFI and the AGFI should be > .95 and > .90, respectively.
\item \strong{NFI/NNFI/TLI}: The (Non) Normed Fit Index. An NFI of 0.95, indicates the
model of interest improves the fit by 95\\% relative to the null model. The
NNFI (also called the Tucker Lewis index; TLI) is preferable for smaller
samples. They should be > .90 (Byrne, 1994) or > .95 (Schumacker & Lomax,
2004).
\item \strong{CFI}: The Comparative Fit Index is a revised form of NFI. Not very
sensitive to sample size (Fan, Thompson, & Wang, 1999). Compares the fit of a
target model to the fit of an independent, or null, model. It should be >
.90.
\item \strong{RMSEA}: The Root Mean Square Error of Approximation is a
parsimony-adjusted index. Values closer to 0 represent a good fit. It should
be < .08 or < .05. The p-value printed with it tests the hypothesis that
RMSEA is less than or equal to .05 (a cutoff sometimes used for good fit),
and thus should be not significant.
\item \strong{RMR/SRMR}: the (Standardized) Root Mean Square Residual represents the
square-root of the difference between the residuals of the sample covariance
matrix and the hypothesized model. As the RMR can be sometimes hard to
interpret, better to use SRMR. Should be < .08.
\item \strong{RFI}: the Relative Fit Index, also known as RHO1, is not guaranteed to
vary from 0 to 1. However, RFI close to 1 indicates a good fit.
\item \strong{IFI}: the Incremental Fit Index (IFI) adjusts the Normed Fit Index (NFI)
for sample size and degrees of freedom (Bollen's, 1989). Over 0.90 is a good
fit, but the index can exceed 1.
\item \strong{PNFI}: the Parsimony-Adjusted Measures Index. There is no commonly
agreed-upon cutoff value for an acceptable model for this index. Should be >
0.50.
}

See the documentation for \code{\link[lavaan:fitmeasures]{fitmeasures()}}.
}

\subsection{What to report}{

For structural equation models (SEM), Kline (2015) suggests that at a minimum
the following indices should be reported: The model \strong{chi-square}, the
\strong{RMSEA}, the \strong{CFI} and the \strong{SRMR}.
}
}
\note{
When possible, it is recommended to report dynamic cutoffs of fit
indices. See https://dynamicfit.app/cfa/.
}
\examples{
interpret_gfi(c(.5, .99))
interpret_agfi(c(.5, .99))
interpret_nfi(c(.5, .99))
interpret_nnfi(c(.5, .99))
interpret_cfi(c(.5, .99))
interpret_rmsea(c(.07, .04))
interpret_srmr(c(.5, .99))
interpret_rfi(c(.5, .99))
interpret_ifi(c(.5, .99))
interpret_pnfi(c(.5, .99))

# Structural Equation Models (SEM)
if (require("lavaan")) {
  structure <- " ind60 =~ x1 + x2 + x3
                 dem60 =~ y1 + y2 + y3
                 dem60 ~ ind60 "
  model <- lavaan::sem(structure, data = lavaan::PoliticalDemocracy)
  interpret(model)
}

}
\references{
\itemize{
\item Awang, Z. (2012). A handbook on SEM. Structural equation modeling.
\item Byrne, B. M. (1994). Structural equation modeling with EQS and EQS/Windows.
Thousand Oaks, CA: Sage Publications.
\item Tucker, L. R., \& Lewis, C. (1973). The reliability coefficient for maximum
likelihood factor analysis. Psychometrika, 38, 1-10.
\item Schumacker, R. E., \& Lomax, R. G. (2004). A beginner's guide to structural
equation modeling, Second edition. Mahwah, NJ: Lawrence Erlbaum Associates.
\item Fan, X., B. Thompson, \& L. Wang (1999). Effects of sample size, estimation
method, and model specification on structural equation modeling fit indexes.
Structural Equation Modeling, 6, 56-83.
\item Kline, R. B. (2015). Principles and practice of structural equation
modeling. Guilford publications.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xtab.R
\name{phi}
\alias{phi}
\alias{cohens_w}
\alias{cramers_v}
\alias{pearsons_c}
\alias{oddsratio}
\alias{riskratio}
\alias{cohens_h}
\alias{cohens_g}
\title{Effect size for contingency tables}
\usage{
phi(x, y = NULL, ci = 0.95, alternative = "greater", adjust = FALSE, ...)

cohens_w(x, y = NULL, ci = 0.95, alternative = "greater", adjust = FALSE, ...)

cramers_v(x, y = NULL, ci = 0.95, alternative = "greater", adjust = FALSE, ...)

pearsons_c(
  x,
  y = NULL,
  ci = 0.95,
  alternative = "greater",
  adjust = FALSE,
  ...
)

oddsratio(x, y = NULL, ci = 0.95, alternative = "two.sided", log = FALSE, ...)

riskratio(x, y = NULL, ci = 0.95, alternative = "two.sided", log = FALSE, ...)

cohens_h(x, y = NULL, ci = 0.95, alternative = "two.sided", ...)

cohens_g(x, y = NULL, ci = 0.95, alternative = "two.sided", ...)
}
\arguments{
\item{x}{a numeric vector or matrix. \code{x} and \code{y} can also
    both be factors.}

\item{y}{a numeric vector; ignored if \code{x} is a matrix.  If
    \code{x} is a factor, \code{y} should be a factor of the same length.}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"greater"} (two-sided CI; default for
Cramer's \emph{V}, phi (\eqn{\phi}), and Cohen's \emph{w}), \code{"two.sided"} (default
for OR, RR, Cohen's \emph{h} and Cohen's \emph{g}) or \code{"less"} (one-sided CI).
Partial matching is allowed (e.g., \code{"g"}, \code{"l"}, \code{"two"}...). See
\emph{One-Sided CIs} in \link{effectsize_CIs}.}

\item{adjust}{Should the effect size be bias-corrected? Defaults to \code{FALSE}.}

\item{...}{Arguments passed to \code{\link[stats:chisq.test]{stats::chisq.test()}}, such as \code{p}. Ignored
for \code{cohens_g()}.}

\item{log}{Take in or output the log of the ratio (such as in logistic models).}
}
\value{
A data frame with the effect size (\code{Cramers_v}, \code{phi} (possibly with
the suffix \verb{_adjusted}), \code{Odds_ratio}, \code{Risk_ratio} (possibly with the
prefix \code{log_}), \code{Cohens_h}, or \code{Cohens_g}) and its CIs (\code{CI_low} and
\code{CI_high}).
}
\description{
Compute Cramer's \emph{V}, phi (\eqn{\phi}), Cohen's \emph{w} (an alias of phi),
Pearson's contingency coefficient, Odds ratios, Risk ratios, Cohen's \emph{h} and
Cohen's \emph{g} for contingency tables or goodness-of-fit. See details.
}
\details{
Cramer's \emph{V}, phi (\eqn{\phi}) and Pearson's \emph{C} are effect sizes for tests
of independence in 2D contingency tables. For 2-by-k tables, Cramer's \emph{V} and
phi are identical, and are equal to the simple correlation between two
dichotomous variables, ranging between  0 (no dependence) and 1 (perfect
dependence). For larger tables, Cramer's \emph{V} or Pearson's \emph{C} should be used,
as they are bounded between 0-1, whereas phi can be larger than 1 (upper
bound is \verb{sqrt(min(nrow, ncol) - 1))}).
\cr\cr
For goodness-of-fit in 1D tables Pearson's \emph{C} or phi can be used. Phi has no
upper bound (can be arbitrarily large, depending on the expected
distribution), while Pearson's \emph{C} is bounded between 0-1.
\cr\cr
For 2-by-2 contingency tables, Odds ratios, Risk ratios and Cohen's \emph{h} can
also be estimated. Note that these are computed with each \strong{column}
representing the different groups, and the first column representing the
treatment group and the second column baseline (or control). Effects are
given as \code{treatment / control}. If you wish you use rows as groups you must
pass a transposed table, or switch the \code{x} and \code{y} arguments.
\cr\cr
Cohen's \emph{g} is an effect size for dependent (paired) contingency tables
ranging between 0 (perfect symmetry) and 0.5 (perfect asymmetry) (see
\code{\link[stats:mcnemar.test]{stats::mcnemar.test()}}).
}
\section{Confidence Intervals for Cohen's g, OR, RR and Cohen's h}{
For Cohen's \emph{g}, confidence intervals are based on the proportion (\eqn{P = g
+ 0.5}) confidence intervals returned by \code{\link[stats:prop.test]{stats::prop.test()}} (minus 0.5),
which give a good close approximation.
\cr\cr
For Odds ratios, Risk ratios and Cohen's \emph{h}, confidence intervals are
estimated using the standard normal parametric method (see Katz et al., 1978;
Szumilas, 2010).
\cr\cr
See \emph{Confidence (Compatibility) Intervals (CIs)}, \emph{CIs and Significance
Tests}, and \emph{One-Sided CIs} sections for \emph{phi}, Cohen's \emph{w}, Cramer's \emph{V} and
Pearson's \emph{C}.
}

\section{Confidence (Compatibility) Intervals (CIs)}{

Unless stated otherwise, confidence (compatibility) intervals (CIs) are
estimated using the noncentrality parameter method (also called the "pivot
method"). This method finds the noncentrality parameter ("\emph{ncp}") of a
noncentral \emph{t}, \emph{F}, or \eqn{\chi^2} distribution that places the observed
\emph{t}, \emph{F}, or \eqn{\chi^2} test statistic at the desired probability point of
the distribution. For example, if the observed \emph{t} statistic is 2.0, with 50
degrees of freedom, for which cumulative noncentral \emph{t} distribution is \emph{t} =
2.0 the .025 quantile (answer: the noncentral \emph{t} distribution with \emph{ncp} =
.04)? After estimating these confidence bounds on the \emph{ncp}, they are
converted into the effect size metric to obtain a confidence interval for the
effect size (Steiger, 2004).
\cr\cr
For additional details on estimation and troubleshooting, see \link{effectsize_CIs}.
}

\section{CIs and Significance Tests}{

"Confidence intervals on measures of effect size convey all the information
in a hypothesis test, and more." (Steiger, 2004). Confidence (compatibility)
intervals and p values are complementary summaries of parameter uncertainty
given the observed data. A dichotomous hypothesis test could be performed
with either a CI or a p value. The 100 (1 - \eqn{\alpha})\% confidence
interval contains all of the parameter values for which \emph{p} > \eqn{\alpha}
for the current data and model. For example, a 95\% confidence interval
contains all of the values for which p > .05.
\cr\cr
Note that a confidence interval including 0 \emph{does not} indicate that the null
(no effect) is true. Rather, it suggests that the observed data together with
the model and its assumptions combined do not provided clear evidence against
a parameter value of 0 (same as with any other value in the interval), with
the level of this evidence defined by the chosen \eqn{\alpha} level (Rafi &
Greenland, 2020; Schweder & Hjort, 2016; Xie & Singh, 2013). To infer no
effect, additional judgments about what parameter values are "close enough"
to 0 to be negligible are needed ("equivalence testing"; Bauer & Kiesser,
1996).
}

\examples{
M <-
  matrix(c(150, 100, 165,
           130, 50, 65,
           35, 10, 2,
           55, 40, 25), nrow = 4,
         dimnames = list(
           Music = c("Pop", "Rock", "Jazz", "Classic"),
           Study = c("Psych", "Econ", "Law")))
M

# Note that Phi is not bound to [0-1], but instead
# the upper bound for phi is sqrt(min(nrow, ncol) - 1)
phi(M)

cramers_v(M)

pearsons_c(M)


## 2-by-2 tables
## -------------
RCT <-
  matrix(c(71, 30,
           50, 100), nrow = 2, byrow = TRUE,
         dimnames = list(
           Diagnosis = c("Sick", "Recovered"),
           Group = c("Treatment", "Control")))
RCT # note groups are COLUMNS

oddsratio(RCT)
oddsratio(RCT, alternative = "greater")

riskratio(RCT)

cohens_h(RCT)



## Dependent (Paired) Contingency Tables
## -------------------------------------
Performance <-
  matrix(c(794, 150,
           86, 570), nrow = 2,
         dimnames = list(
           "1st Survey" = c("Approve", "Disapprove"),
           "2nd Survey" = c("Approve", "Disapprove")))
Performance

cohens_g(Performance)

}
\references{
\itemize{
\item Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd Ed.). New York: Routledge.
\item Katz, D. J. S. M., Baptista, J., Azen, S. P., & Pike, M. C. (1978). Obtaining confidence intervals for the risk ratio in cohort studies. Biometrics, 469-474.
\item Szumilas, M. (2010). Explaining odds ratios. Journal of the Canadian academy of child and adolescent psychiatry, 19(3), 227.
}
}
\seealso{
\code{\link[=chisq_to_phi]{chisq_to_phi()}} for details regarding estimation and CIs.

Other effect size indices: 
\code{\link{cles}()},
\code{\link{cohens_d}()},
\code{\link{effectsize.BFBayesFactor}()},
\code{\link{eta_squared}()},
\code{\link{rank_biserial}()},
\code{\link{standardize_parameters}()}
}
\concept{effect size indices}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohens_d.R
\name{cohens_d}
\alias{cohens_d}
\alias{hedges_g}
\alias{glass_delta}
\title{Effect size for differences}
\usage{
cohens_d(
  x,
  y = NULL,
  data = NULL,
  pooled_sd = TRUE,
  mu = 0,
  paired = FALSE,
  ci = 0.95,
  alternative = "two.sided",
  verbose = TRUE,
  ...
)

hedges_g(
  x,
  y = NULL,
  data = NULL,
  pooled_sd = TRUE,
  mu = 0,
  paired = FALSE,
  ci = 0.95,
  alternative = "two.sided",
  verbose = TRUE,
  ...,
  correction
)

glass_delta(
  x,
  y = NULL,
  data = NULL,
  mu = 0,
  ci = 0.95,
  alternative = "two.sided",
  verbose = TRUE,
  ...,
  iterations
)
}
\arguments{
\item{x}{A formula, a numeric vector, or a character name of one in \code{data}.}

\item{y}{A numeric vector, a grouping (character / factor) vector, a or a
character  name of one in \code{data}. Ignored if \code{x} is a formula.}

\item{data}{An optional data frame containing the variables.}

\item{pooled_sd}{If \code{TRUE} (default), a \code{\link[=sd_pooled]{sd_pooled()}} is used (assuming equal
variance). Else the mean SD from both groups is used instead.}

\item{mu}{a number indicating the true value of the mean (or
    difference in means if you are performing a two sample test).}

\item{paired}{If \code{TRUE}, the values of \code{x} and \code{y} are considered as paired.
This produces an effect size that is equivalent to the one-sample effect
size on \code{x - y}.}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"two.sided"} (default, two-sided CI),
\code{"greater"} or \code{"less"} (one-sided CI). Partial matching is allowed (e.g.,
\code{"g"}, \code{"l"}, \code{"two"}...). See \emph{One-Sided CIs} in \link{effectsize_CIs}.}

\item{verbose}{Toggle warnings and messages on or off.}

\item{...}{Arguments passed to or from other methods.}

\item{iterations, correction}{deprecated.}
}
\value{
A data frame with the effect size ( \code{Cohens_d}, \code{Hedges_g},
\code{Glass_delta}) and their CIs (\code{CI_low} and \code{CI_high}).
}
\description{
Compute effect size indices for standardized differences: Cohen's \emph{d},
Hedges' \emph{g} and Glass’s \emph{delta} (\eqn{\Delta}). (This function returns the
\strong{population} estimate.)
\cr\cr
Both Cohen's \emph{d} and Hedges' \emph{g} are the estimated the standardized
difference between the means of two populations. Hedges' \emph{g} provides a bias
correction (using the exact method) to Cohen's \emph{d} for small sample sizes.
For sample sizes > 20, the results for both statistics are roughly
equivalent. Glass’s \emph{delta} is appropriate when the standard deviations are
significantly different between the populations, as it uses only the \emph{second}
group's standard deviation.
}
\details{
Set \code{pooled_sd = FALSE} for effect sizes that are to accompany a Welch's
\emph{t}-test (Delacre et al, 2021).
}
\note{
The indices here give the population estimated standardized difference.
Some statistical packages give the sample estimate instead (without
applying Bessel's correction).
}
\section{Confidence (Compatibility) Intervals (CIs)}{

Unless stated otherwise, confidence (compatibility) intervals (CIs) are
estimated using the noncentrality parameter method (also called the "pivot
method"). This method finds the noncentrality parameter ("\emph{ncp}") of a
noncentral \emph{t}, \emph{F}, or \eqn{\chi^2} distribution that places the observed
\emph{t}, \emph{F}, or \eqn{\chi^2} test statistic at the desired probability point of
the distribution. For example, if the observed \emph{t} statistic is 2.0, with 50
degrees of freedom, for which cumulative noncentral \emph{t} distribution is \emph{t} =
2.0 the .025 quantile (answer: the noncentral \emph{t} distribution with \emph{ncp} =
.04)? After estimating these confidence bounds on the \emph{ncp}, they are
converted into the effect size metric to obtain a confidence interval for the
effect size (Steiger, 2004).
\cr\cr
For additional details on estimation and troubleshooting, see \link{effectsize_CIs}.
}

\section{CIs and Significance Tests}{

"Confidence intervals on measures of effect size convey all the information
in a hypothesis test, and more." (Steiger, 2004). Confidence (compatibility)
intervals and p values are complementary summaries of parameter uncertainty
given the observed data. A dichotomous hypothesis test could be performed
with either a CI or a p value. The 100 (1 - \eqn{\alpha})\% confidence
interval contains all of the parameter values for which \emph{p} > \eqn{\alpha}
for the current data and model. For example, a 95\% confidence interval
contains all of the values for which p > .05.
\cr\cr
Note that a confidence interval including 0 \emph{does not} indicate that the null
(no effect) is true. Rather, it suggests that the observed data together with
the model and its assumptions combined do not provided clear evidence against
a parameter value of 0 (same as with any other value in the interval), with
the level of this evidence defined by the chosen \eqn{\alpha} level (Rafi &
Greenland, 2020; Schweder & Hjort, 2016; Xie & Singh, 2013). To infer no
effect, additional judgments about what parameter values are "close enough"
to 0 to be negligible are needed ("equivalence testing"; Bauer & Kiesser,
1996).
}

\examples{
\donttest{
data(mtcars)
mtcars$am <- factor(mtcars$am)

# Two Independent Samples ----------

(d <- cohens_d(mpg ~ am, data = mtcars))
# Same as:
# cohens_d("mpg", "am", data = mtcars)
# cohens_d(mtcars$mpg[mtcars$am=="0"], mtcars$mpg[mtcars$am=="1"])

# More options:
cohens_d(mpg ~ am, data = mtcars, pooled_sd = FALSE)
cohens_d(mpg ~ am, data = mtcars, mu = -5)
cohens_d(mpg ~ am, data = mtcars, alternative = "less")
hedges_g(mpg ~ am, data = mtcars)
glass_delta(mpg ~ am, data = mtcars)


# One Sample ----------

cohens_d(wt ~ 1, data = mtcars)

# same as:
# cohens_d("wt", data = mtcars)
# cohens_d(mtcars$wt)

# More options:
cohens_d(wt ~ 1, data = mtcars, mu = 3)
hedges_g(wt ~ 1, data = mtcars, mu = 3)


# Paired Samples ----------

data(sleep)

cohens_d(Pair(extra[group == 1], extra[group == 2]) ~ 1, data = sleep)

# same as:
# cohens_d(sleep$extra[sleep$group == 1], sleep$extra[sleep$group == 2], paired = TRUE)

# More options:
cohens_d(Pair(extra[group == 1], extra[group == 2]) ~ 1, data = sleep, mu = -1)
hedges_g(Pair(extra[group == 1], extra[group == 2]) ~ 1, data = sleep)


# Interpretation -----------------------
interpret_cohens_d(-1.48, rules = "cohen1988")
interpret_hedges_g(-1.48, rules = "sawilowsky2009")
interpret_glass_delta(-1.48, rules = "gignac2016")
# Or:
interpret(d, rules = "sawilowsky2009")

# Common Language Effect Sizes
d_to_cles(1.48)
# Or:
print(d, append_CLES = TRUE)
}

}
\references{
\itemize{
\item Algina, J., Keselman, H. J., & Penfield, R. D. (2006). Confidence intervals
for an effect size when variances are not equal. Journal of Modern Applied
Statistical Methods, 5(1), 2.
\item Cohen, J. (1988). Statistical power analysis for the behavioral
sciences (2nd Ed.). New York: Routledge.
\item Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021, May 7). Why
Hedges’ g*s based on the non-pooled standard deviation should be reported
with Welch's t-test. https://doi.org/10.31234/osf.io/tu6mp
\item Hedges, L. V. & Olkin, I. (1985). Statistical methods for
meta-analysis. Orlando, FL: Academic Press.
\item Hunter, J. E., & Schmidt, F. L. (2004). Methods of meta-analysis:
Correcting error and bias in research findings. Sage.
}
}
\seealso{
\code{\link[=d_to_cles]{d_to_cles()}} \code{\link[=sd_pooled]{sd_pooled()}}

Other effect size indices: 
\code{\link{cles}()},
\code{\link{effectsize.BFBayesFactor}()},
\code{\link{eta_squared}()},
\code{\link{phi}()},
\code{\link{rank_biserial}()},
\code{\link{standardize_parameters}()}
}
\concept{effect size indices}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_r2.R
\name{interpret_r2}
\alias{interpret_r2}
\title{Interpret coefficient of determination (R2)}
\usage{
interpret_r2(r2, rules = "cohen1988")
}
\arguments{
\item{r2}{Value or vector of R2 values.}

\item{rules}{Can be \code{"cohen1988"} (default), \code{"falk1992"}, \code{"chin1998"},
\code{"hair2011"}, or custom set of \code{\link[=rules]{rules()}}].}
}
\description{
Interpret coefficient of determination (R2)
}
\section{Rules}{

\subsection{For Linear Regression}{
\itemize{
\item Cohen (1988) (\code{"cohen1988"}; default)
\itemize{
\item \strong{R2 < 0.02} - Very weak
\item \strong{0.02 <= R2 < 0.13} - Weak
\item \strong{0.13 <= R2 < 0.26} - Moderate
\item \strong{R2 >= 0.26} - Substantial
}
\item Falk & Miller (1992) (\code{"falk1992"})
\itemize{
\item \strong{R2 < 0.1} - Negligible
\item \strong{R2 >= 0.1} - Adequate
}
}
}

\subsection{For PLS / SEM R-Squared of \emph{latent} variables}{
\itemize{
\item Chin, W. W. (1998) (\code{"chin1998"})
\itemize{
\item \strong{R2 < 0.19} - Very weak
\item \strong{0.19 <= R2 < 0.33} - Weak
\item \strong{0.33 <= R2 < 0.67} - Moderate
\item \strong{R2 >= 0.67} - Substantial
}
\item Hair et al. (2011) (\code{"hair2011"})
\itemize{
\item \strong{R2 < 0.25} - Very weak
\item \strong{0.25 <= R2 < 0.50} - Weak
\item \strong{0.50 <= R2 < 0.75} - Moderate
\item \strong{R2 >= 0.75} - Substantial
}
}
}
}

\examples{
interpret_r2(.02)
interpret_r2(c(.5, .02))
}
\references{
\itemize{
\item Cohen, J. (1988). Statistical power analysis for the behavioral sciences
(2nd Ed.). New York: Routledge.
\item Falk, R. F., & Miller, N. B. (1992). A primer for soft modeling. University
of Akron Press.
\item Chin, W. W. (1998). The partial least squares approach to structural
equation modeling. Modern methods for business research, 295(2), 295-336.
\item Hair, J. F., Ringle, C. M., & Sarstedt, M. (2011). PLS-SEM: Indeed a silver
bullet. Journal of Marketing theory and Practice, 19(2), 139-152.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_between_OR_to_RR.R
\name{oddsratio_to_riskratio}
\alias{oddsratio_to_riskratio}
\alias{riskratio_to_oddsratio}
\title{Convert between Odds ratios and Risk ratios}
\usage{
oddsratio_to_riskratio(OR, p0, log = FALSE, ...)

riskratio_to_oddsratio(RR, p0, log = FALSE)
}
\arguments{
\item{OR, RR}{Risk ratio of \code{p1/p0} or Odds ratio of \code{odds(p1)/odds(p0)},
possibly log-ed. \code{OR} can also be a logistic regression model.}

\item{p0}{Baseline risk}

\item{log}{Take in or output the log of the ratio (such as in logistic models).}

\item{...}{Arguments passed to and from other methods.}
}
\value{
Converted index, or if \code{OR} is a logistic regression model, a
parameter table with the converted indices.
}
\description{
Convert between Odds ratios and Risk ratios
}
\examples{
p0 <- 0.4
p1 <- 0.7

(OR <- probs_to_odds(p1) / probs_to_odds(p0))
(RR <- p1 / p0)

riskratio_to_oddsratio(RR, p0 = p0)
oddsratio_to_riskratio(OR, p0 = p0)

m <- glm(am ~ factor(cyl), data = mtcars,
         family = binomial())
oddsratio_to_riskratio(m)
}
\references{
Grant, R. L. (2014). Converting an odds ratio to a range of plausible
relative risks for better communication of research findings. Bmj, 348,
f7450.
}
\seealso{
Other convert between effect sizes: 
\code{\link{d_to_cles}()},
\code{\link{d_to_r}()},
\code{\link{eta2_to_f2}()},
\code{\link{odds_to_probs}()}
}
\concept{convert between effect sizes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret.R
\name{rules}
\alias{rules}
\alias{is.rules}
\title{Interpretation Grid}
\usage{
rules(values, labels = NULL, name = NULL, right = TRUE)

is.rules(x)
}
\arguments{
\item{values}{Vector of reference values (edges defining categories or
critical values).}

\item{labels}{Labels associated with each category. If \code{NULL}, will try to
infer it from \code{values} (if it is a named vector or a list), otherwise, will
return the breakpoints.}

\item{name}{Name of the set of rules (will be printed).}

\item{right}{logical, for threshold-type rules, indicating if the thresholds
themselves should be included in the interval to the right (lower values)
or in the interval to the left (higher values).}

\item{x}{An arbitrary R object.}
}
\description{
Create a container for interpretation rules of thumb. Usually used in conjunction with \link{interpret}.
}
\examples{
rules(c(0.05), c("significant", "not significant"), right = FALSE)
rules(c(0.2, 0.5, 0.8), c("small", "medium", "large"))
rules(c("small" = 0.2, "medium" = 0.5), name = "Cohen's Rules")
}
\seealso{
interpret
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize_parameters.R
\name{standardize_parameters}
\alias{standardize_parameters}
\alias{standardize_posteriors}
\title{Parameters standardization}
\usage{
standardize_parameters(
  model,
  method = "refit",
  ci = 0.95,
  robust = FALSE,
  two_sd = FALSE,
  include_response = TRUE,
  verbose = TRUE,
  parameters,
  ...
)

standardize_posteriors(
  model,
  method = "refit",
  robust = FALSE,
  two_sd = FALSE,
  include_response = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{A statistical model.}

\item{method}{The method used for standardizing the parameters. Can be
\code{"refit"} (default), \code{"posthoc"}, \code{"smart"}, \code{"basic"} or \code{"pseudo"}. See
'Details'.}

\item{ci}{Confidence Interval (CI) level}

\item{robust}{Logical, if \code{TRUE}, centering is done by subtracting the
median from the variables and dividing it by the median absolute deviation
(MAD). If \code{FALSE}, variables are standardized by subtracting the
mean and dividing it by the standard deviation (SD).}

\item{two_sd}{If \code{TRUE}, the variables are scaled by two times the deviation
(SD or MAD depending on \code{robust}). This method can be useful to obtain
model coefficients of continuous parameters comparable to coefficients
related to binary predictors, when applied to \strong{the predictors} (not the
outcome) (Gelman, 2008).}

\item{include_response}{If \code{TRUE} (default), the response value will also be
standardized. If \code{FALSE}, only the predictors will be standardized. For
GLMs the response value will never be standardized (see \emph{Generalized Linear
Models} section).}

\item{verbose}{Toggle warnings and messages on or off.}

\item{parameters}{Deprecated.}

\item{...}{For \code{standardize_parameters()}, arguments passed to
\link[parameters:model_parameters]{parameters::model_parameters}, such as:
\itemize{
\item \code{ci_method}, \code{centrality} for Mixed models and Bayesian models...
\item \code{exponentiate}, ...
\item etc.
}}
}
\value{
A data frame with the standardized parameters (\verb{Std_*}, depending on
the model type) and their CIs (\code{CI_low} and \code{CI_high}). Where applicable,
standard errors (SEs) are returned as an attribute (\code{attr(x, "standard_error")}).
}
\description{
Compute standardized model parameters (coefficients).
}
\section{Standardization Methods:}{
\itemize{
\item \strong{refit}: This method is based on a complete model re-fit with a
standardized version of the data. Hence, this method is equal to
standardizing the variables before fitting the model. It is the "purest" and
the most accurate (Neter et al., 1989), but it is also the most
computationally costly and long (especially for heavy models such as Bayesian
models). This method is particularly recommended for complex models that
include interactions or transformations (e.g., polynomial or spline terms).
The \code{robust} (default to \code{FALSE}) argument enables a robust standardization
of data, i.e., based on the \code{median} and \code{MAD} instead of the \code{mean} and
\code{SD}. \strong{See \code{\link[=standardize]{standardize()}} for more details.}
\itemize{
\item \strong{Note} that \code{standardize_parameters(method = "refit")} may not return
the same results as fitting a model on data that has been standardized with
\code{standardize()}; \code{standardize_parameters()} used the data used by the model
fitting function, which might not be same data if there are missing values.
see the \code{remove_na} argument in \code{standardize()}.
}
\item \strong{posthoc}: Post-hoc standardization of the parameters, aiming at
emulating the results obtained by "refit" without refitting the model. The
coefficients are divided by the standard deviation (or MAD if \code{robust}) of
the outcome (which becomes their expression 'unit'). Then, the coefficients
related to numeric variables are additionally multiplied by the standard
deviation (or MAD if \code{robust}) of the related terms, so that they correspond
to changes of 1 SD of the predictor (e.g., "A change in 1 SD of \code{x} is
related to a change of 0.24 of the SD of \code{y}). This does not apply to binary
variables or factors, so the coefficients are still related to changes in
levels. This method is not accurate and tend to give aberrant results when
interactions are specified.
\item \strong{basic}: This method is similar to \code{method = "posthoc"}, but treats all
variables as continuous: it also scales the coefficient by the standard
deviation of model's matrix' parameter of factors levels (transformed to
integers) or binary predictors. Although being inappropriate for these cases,
this method is the one implemented by default in other software packages,
such as \code{\link[lm.beta:lm.beta]{lm.beta::lm.beta()}}.
\item \strong{smart} (Standardization of Model's parameters with Adjustment,
Reconnaissance and Transformation - \emph{experimental}): Similar to \code{method = "posthoc"} in that it does not involve model refitting. The difference is
that the SD (or MAD if \code{robust}) of the response is computed on the relevant
section of the data. For instance, if a factor with 3 levels A (the
intercept), B and C is entered as a predictor, the effect corresponding to B
vs. A will be scaled by the variance of the response at the intercept only.
As a results, the coefficients for effects of factors are similar to a Glass'
delta.
\item \strong{pseudo} (\emph{for 2-level (G)LMMs only}): In this (post-hoc) method, the
response and the predictor are standardized based on the level of prediction
(levels are detected with \code{\link[performance:check_heterogeneity_bias]{performance::check_heterogeneity_bias()}}): Predictors
are standardized based on their SD at level of prediction (see also
\code{\link[datawizard:demean]{datawizard::demean()}}); The outcome (in linear LMMs) is standardized based
on a fitted random-intercept-model, where \code{sqrt(random-intercept-variance)}
is used for level 2 predictors, and \code{sqrt(residual-variance)} is used for
level 1 predictors (Hoffman 2015, page 342). A warning is given when a
within-group varialbe is found to have access between-group variance.
}
}

\section{Transformed Variables}{
When the model's formula contains transformations (e.g. \code{y ~ exp(X)}) \code{method = "refit"} will give different results compared to \code{method = "basic"}
(\code{"posthoc"} and \code{"smart"} do not support such transformations): While
\code{"refit"} standardizes the data \emph{prior} to the transformation (e.g.
equivalent to \code{exp(scale(X))}), the \code{"basic"} method standardizes the
transformed data (e.g. equivalent to \code{scale(exp(X))}).
\cr\cr
See the \emph{Transformed Variables} section in \code{\link[=standardize.default]{standardize.default()}} for more
details on how different transformations are dealt with when \code{method = "refit"}.
}

\section{Confidence Intervals}{
The returned confidence intervals are re-scaled versions of the
unstandardized confidence intervals, and not "true" confidence intervals of
the standardized coefficients (cf. Jones & Waller, 2015).
}

\section{Generalized Linear Models}{
Standardization for generalized linear models (GLM, GLMM, etc) is done only
with respect to the predictors (while the outcome remains as-is,
unstandardized) - maintaining the interpretability of the coefficients (e.g.,
in a binomial model: the exponent of the standardized parameter is the OR of
a change of 1 SD in the predictor, etc.)
}

\section{Dealing with Factors}{
\code{standardize(model)} or \code{standardize_parameters(model, method = "refit")} do
\emph{not} standardized categorical predictors (i.e. factors) / their
dummy-variables, which may be a different behaviour compared to other R
packages (such as \pkg{lm.beta}) or other software packages (like SPSS). To
mimic such behaviours, either use \code{standardize_parameters(model, method = "basic")} to obtain post-hoc standardized parameters, or standardize the data
with \code{datawizard::standardize(data, force = TRUE)} \emph{before} fitting the
model.
}

\examples{
library(effectsize)

model <- lm(len ~ supp * dose, data = ToothGrowth)
standardize_parameters(model, method = "refit")
\donttest{
standardize_parameters(model, method = "posthoc")
standardize_parameters(model, method = "smart")
standardize_parameters(model, method = "basic")

# Robust and 2 SD
standardize_parameters(model, robust = TRUE)
standardize_parameters(model, two_sd = TRUE)


model <- glm(am ~ cyl * mpg, data = mtcars, family = "binomial")
standardize_parameters(model, method = "refit")
standardize_parameters(model, method = "posthoc")
standardize_parameters(model, method = "basic", exponentiate = TRUE)
}

\donttest{
if (require("lme4")) {
  m <- lmer(mpg ~ cyl + am + vs + (1 | cyl), mtcars)
  standardize_parameters(m, method = "pseudo", ci_method = "satterthwaite")
}


\dontrun{
if (require("rstanarm")) {
  model <- stan_glm(rating ~ critical + privileges, data = attitude, refresh = 0)
  standardize_posteriors(model, method = "refit")
  standardize_posteriors(model, method = "posthoc")
  standardize_posteriors(model, method = "smart")
  head(standardize_posteriors(model, method = "basic"))
}
}
}

}
\references{
\itemize{
\item Hoffman, L. (2015). Longitudinal analysis: Modeling within-person fluctuation and change. Routledge.
\item Jones, J. A., & Waller, N. G. (2015). The normal-theory and asymptotic distribution-free (ADF) covariance matrix of standardized regression coefficients: theoretical extensions and finite sample behavior. Psychometrika, 80(2), 365-378.
\item Neter, J., Wasserman, W., & Kutner, M. H. (1989). Applied linear regression models.
\item Gelman, A. (2008). Scaling regression inputs by dividing by two standard deviations. Statistics in medicine, 27(15), 2865-2873.
}
}
\seealso{
Other standardize: 
\code{\link{standardize.default}()},
\code{\link{standardize_info}()}

Other effect size indices: 
\code{\link{cles}()},
\code{\link{cohens_d}()},
\code{\link{effectsize.BFBayesFactor}()},
\code{\link{eta_squared}()},
\code{\link{phi}()},
\code{\link{rank_biserial}()}
}
\concept{effect size indices}
\concept{standardize}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_pd.R
\name{interpret_pd}
\alias{interpret_pd}
\title{Interpret Probability of Direction (pd)}
\usage{
interpret_pd(pd, rules = "default", ...)
}
\arguments{
\item{pd}{Value or vector of probabilities of direction.}

\item{rules}{Can be \code{"default"}, \code{"makowski2019"} or a custom set of
\code{\link[=rules]{rules()}}.}

\item{...}{Not directly used.}
}
\description{
Interpret Probability of Direction (pd)
}
\section{Rules}{

\itemize{
\item Default (i.e., equivalent to p-values)
\itemize{
\item \strong{pd <= 0.975} - not significant
\item \strong{pd > 0.975} - significant
}
\item Makowski et al. (2019) (\code{"makowski2019"})
\itemize{
\item \strong{pd <= 0.95} - uncertain
\item \strong{pd > 0.95} - possibly existing
\item \strong{pd > 0.97} - likely existing
\item \strong{pd > 0.99} - probably existing
\item \strong{pd > 0.999} - certainly existing
}
}
}

\examples{
interpret_pd(.98)
interpret_pd(c(.96, .99), rules = "makowski2019")
}
\references{
\itemize{
\item Makowski, D., Ben-Shachar, M. S., Chen, S. H., \& Lüdecke, D. (2019). Indices of effect existence and significance in the Bayesian framework. Frontiers in psychology, 10, 2767.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_direction.R
\name{interpret_direction}
\alias{interpret_direction}
\title{Interpret direction}
\usage{
interpret_direction(x)
}
\arguments{
\item{x}{Numeric value.}
}
\description{
Interpret direction
}
\examples{
interpret_direction(.02)
interpret_direction(c(.5, -.02))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_omega_squared.R
\name{interpret_omega_squared}
\alias{interpret_omega_squared}
\alias{interpret_eta_squared}
\alias{interpret_epsilon_squared}
\title{Interpret ANOVA effect size}
\usage{
interpret_omega_squared(es, rules = "field2013", ...)

interpret_eta_squared(es, rules = "field2013", ...)

interpret_epsilon_squared(es, rules = "field2013", ...)
}
\arguments{
\item{es}{Value or vector of eta / omega / epsilon squared values.}

\item{rules}{Can be \code{"field2013"} (default), \code{"cohen1992"} or custom set of \code{\link[=rules]{rules()}}.}

\item{...}{Not used for now.}
}
\description{
Interpret ANOVA effect size
}
\section{Rules}{

\itemize{
\item Field (2013) (\code{"field2013"}; default)
\itemize{
\item \strong{ES < 0.01} - Very small
\item \strong{0.01 <= ES < 0.06} - Small
\item \strong{0.16 <= ES < 0.14} - Medium
\item **ES >= 0.14 ** - Large
}
\item Cohen (1992) (\code{"cohen1992"}) applicable to one-way anova, or to \emph{partial}
eta / omega / epsilon squared in multi-way anova.
\itemize{
\item \strong{ES < 0.02} - Very small
\item \strong{0.02 <= ES < 0.13} - Small
\item \strong{0.13 <= ES < 0.26} - Medium
\item \strong{ES >= 0.26} - Large
}
}
}

\examples{
interpret_eta_squared(.02)
interpret_eta_squared(c(.5, .02), rules = "cohen1992")
}
\references{
\itemize{
\item Field, A (2013) Discovering statistics using IBM SPSS Statistics. Fourth
Edition. Sage:London.
\item Cohen, J. (1992). A power primer. Psychological bulletin, 112(1), 155.
}
}
\seealso{
https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/effectSize/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank_effectsizes.R
\name{rank_biserial}
\alias{rank_biserial}
\alias{cliffs_delta}
\alias{rank_epsilon_squared}
\alias{kendalls_w}
\title{Effect size for non-parametric (rank sum) tests}
\usage{
rank_biserial(
  x,
  y = NULL,
  data = NULL,
  mu = 0,
  ci = 0.95,
  alternative = "two.sided",
  paired = FALSE,
  verbose = TRUE,
  ...,
  iterations
)

cliffs_delta(
  x,
  y = NULL,
  data = NULL,
  mu = 0,
  ci = 0.95,
  alternative = "two.sided",
  verbose = TRUE,
  ...
)

rank_epsilon_squared(
  x,
  groups,
  data = NULL,
  ci = 0.95,
  alternative = "greater",
  iterations = 200,
  ...
)

kendalls_w(
  x,
  groups,
  blocks,
  data = NULL,
  ci = 0.95,
  alternative = "greater",
  iterations = 200,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{Can be one of:
\itemize{
\item A numeric vector, or a character name of one in \code{data}.
\item A formula in to form of \code{DV ~ groups} (for \code{rank_biserial()} and
\code{rank_epsilon_squared()}) or \code{DV ~ groups | blocks} (for \code{kendalls_w()};
See details for the \code{blocks} and \code{groups} terminology used here).
\item A list of vectors (for \code{rank_epsilon_squared()}).
\item A matrix of \verb{blocks x groups} (for \code{kendalls_w()}). See details for the
\code{blocks} and \code{groups} terminology used here.
}}

\item{y}{An optional numeric vector of data values to compare to \code{x}, or a
character name of one in \code{data}. Ignored if \code{x} is not a vector.}

\item{data}{An optional data frame containing the variables.}

\item{mu}{a number indicating the value around which (a-)symmetry (for
one-sample or paired samples) or shift (for independent samples) is to be
estimated. See \link[stats:wilcox.test]{stats::wilcox.test}.}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"two.sided"} (two-sided CI; default for
rank-biserial correlation and Cliff's \emph{delta}), \code{"greater"} (default for
rank epsilon squared and Kendall's \emph{W}) or \code{"less"} (one-sided CI). Partial
matching is allowed (e.g., \code{"g"}, \code{"l"}, \code{"two"}...). See \emph{One-Sided CIs}
in \link{effectsize_CIs}.}

\item{paired}{If \code{TRUE}, the values of \code{x} and \code{y} are considered as paired.
This produces an effect size that is equivalent to the one-sample effect
size on \code{x - y}.}

\item{verbose}{Toggle warnings and messages on or off.}

\item{...}{Arguments passed to or from other methods.}

\item{iterations}{The number of bootstrap replicates for computing confidence
intervals. Only applies when \code{ci} is not \code{NULL}. (Deprecated for
\code{rank_biserial()}).}

\item{groups, blocks}{A factor vector giving the group / block for the
corresponding elements of \code{x}, or a character name of one in \code{data}.
Ignored if \code{x} is not a vector.}
}
\value{
A data frame with the effect size (\code{r_rank_biserial},
\code{rank_epsilon_squared} or \code{Kendalls_W}) and its CI (\code{CI_low} and
\code{CI_high}).
}
\description{
Compute the rank-biserial correlation (\eqn{r_{rb}}{r_rb}), Cliff's \emph{delta}
(\eqn{\delta}), rank epsilon squared (\eqn{\varepsilon^2}{\epsilon^2}), and
Kendall's \emph{W} effect sizes for non-parametric (rank sum) tests.
}
\details{
The rank-biserial correlation is appropriate for non-parametric tests of
differences - both for the one sample or paired samples case, that would
normally be tested with Wilcoxon's Signed Rank Test (giving the
\strong{matched-pairs} rank-biserial correlation) and for two independent samples
case, that would normally be tested with Mann-Whitney's \emph{U} Test (giving
\strong{Glass'} rank-biserial correlation). See \link[stats:wilcox.test]{stats::wilcox.test}. In both
cases, the correlation represents the difference between the proportion of
favorable and unfavorable pairs / signed ranks (Kerby, 2014). Values range
from \code{-1} (\emph{all} values of the second sample are larger than \emph{all} the values
of the first sample) to \code{+1} (\emph{all} values of the second sample are smaller
than \emph{all} the values of the first sample). Cliff's \emph{delta} is an alias to
the rank-biserial correlation in the two sample case.
\cr\cr
The rank epsilon squared is appropriate for non-parametric tests of
differences between 2 or more samples (a rank based ANOVA). See
\link[stats:kruskal.test]{stats::kruskal.test}. Values range from 0 to 1, with larger values
indicating larger differences between groups.
\cr\cr
Kendall's \emph{W} is appropriate for non-parametric tests of differences between
2 or more dependent samples (a rank based rmANOVA), where each \code{group} (e.g.,
experimental condition) was measured for each \code{block} (e.g., subject). This
measure is also common as a measure of reliability of the rankings of the
\code{groups} between raters (\code{blocks}). See \link[stats:friedman.test]{stats::friedman.test}. Values range
from 0 to 1, with larger values indicating larger differences between groups
/ higher agreement between raters.
\subsection{Ties}{

When tied values occur, they are each given the average of the ranks that
would have been given had no ties occurred. No other corrections have been
implemented yet.
}
}
\section{Confidence Intervals}{
Confidence intervals for the rank-biserial correlation (and Cliff's \emph{delta})
are estimated using the normal approximation (via Fisher's transformation).
Confidence intervals for rank Epsilon squared, and Kendall's \emph{W} are
estimated using the bootstrap method (using the \code{{boot}} package).
}

\examples{
\donttest{
data(mtcars)
mtcars$am <- factor(mtcars$am)
mtcars$cyl <- factor(mtcars$cyl)

# Rank Biserial Correlation
# =========================

# Two Independent Samples ----------
(rb <- rank_biserial(mpg ~ am, data = mtcars))
# Same as:
# rank_biserial("mpg", "am", data = mtcars)
# rank_biserial(mtcars$mpg[mtcars$am=="0"], mtcars$mpg[mtcars$am=="1"])

# More options:
rank_biserial(mpg ~ am, data = mtcars, mu = -5)
print(rb, append_CLES = TRUE)


# One Sample ----------
rank_biserial(wt ~ 1, data = mtcars, mu = 3)
# same as:
# rank_biserial("wt", data = mtcars, mu = 3)
# rank_biserial(mtcars$wt, mu = 3)


# Paired Samples ----------
dat <- data.frame(Cond1 = c(1.83, 0.5, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.3),
                  Cond2 = c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29))
(rb <- rank_biserial(Pair(Cond1, Cond2) ~ 1, data = dat, paired = TRUE))

# same as:
# rank_biserial(dat$Cond1, dat$Cond2, paired = TRUE)

interpret_rank_biserial(0.78)
interpret(rb, rules = "funder2019")


# Rank Epsilon Squared
# ====================

rank_epsilon_squared(mpg ~ cyl, data = mtcars)



# Kendall's W
# ===========
dat <- data.frame(cond = c("A", "B", "A", "B", "A", "B"),
                  ID = c("L", "L", "M", "M", "H", "H"),
                  y = c(44.56, 28.22, 24, 28.78, 24.56, 18.78))
(W <- kendalls_w(y ~ cond | ID, data = dat, verbose = FALSE))

interpret_kendalls_w(0.11)
interpret(W, rules = "landis1977")
}

}
\references{
\itemize{
\item Cureton, E. E. (1956). Rank-biserial correlation. Psychometrika, 21(3),
287-290.
\item Glass, G. V. (1965). A ranking variable analogue of biserial correlation:
Implications for short-cut item analysis. Journal of Educational Measurement,
2(1), 91-95.
\item Kendall, M.G. (1948) Rank correlation methods. London: Griffin.
\item Kerby, D. S. (2014). The simple difference formula: An approach to teaching
nonparametric correlation. Comprehensive Psychology, 3, 11-IT.
\item King, B. M., & Minium, E. W. (2008). Statistical reasoning in the
behavioral sciences. John Wiley & Sons Inc.
\item Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal
questions. Psychological bulletin, 114(3), 494.
\item Tomczak, M., & Tomczak, E. (2014). The need to report effect size estimates
revisited. An overview of some recommended measures of effect size.
}
}
\seealso{
Other effect size indices: 
\code{\link{cles}()},
\code{\link{cohens_d}()},
\code{\link{effectsize.BFBayesFactor}()},
\code{\link{eta_squared}()},
\code{\link{phi}()},
\code{\link{standardize_parameters}()}
}
\concept{effect size indices}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{hardlyworking}
\alias{hardlyworking}
\title{Workers' salary and other information}
\format{
A data frame with 500 rows and 5 variables:
\describe{
\item{salary}{Salary, in Shmekels}
\item{xtra_hours}{Number of overtime hours (on average, per week)}
\item{n_comps}{Number of compliments given to the boss (observed over the last week)}
\item{age}{Age in years}
\item{seniority}{How many years with the company}
}
}
\description{
A sample (simulated) dataset, used in tests and some examples.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_p.R
\name{interpret_p}
\alias{interpret_p}
\title{Interpret p-values}
\usage{
interpret_p(p, rules = "default")
}
\arguments{
\item{p}{Value or vector of p-values.}

\item{rules}{Can be \code{"default"}, \code{"rss"} (for \emph{Redefine statistical
significance} rules) or custom set of \code{\link[=rules]{rules()}}.}
}
\description{
Interpret p-values
}
\section{Rules}{

\itemize{
\item Default
\itemize{
\item \strong{p >= 0.05} - Not significant
\item \strong{p < 0.05} - Significant
}
\item Benjamin et al. (2018) (\code{"rss"})
\itemize{
\item \strong{p >= 0.05} - Not significant
\item \strong{0.005 <= p < 0.05} - Suggestive
\item \strong{p < 0.005} - Significant
}
}
}

\examples{
interpret_p(c(.5, .02, 0.001))
interpret_p(c(.5, .02, 0.001), rules = "rss")
}
\references{
\itemize{
\item Benjamin, D. J., Berger, J. O., Johannesson, M., Nosek, B. A., Wagenmakers, E. J., Berk, R., ... & Cesarini, D. (2018). Redefine statistical significance. Nature Human Behaviour, 2(1), 6-10.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R, R/print.effectsize_table.R
\name{plot.effectsize_table}
\alias{plot.effectsize_table}
\alias{plot.equivalence_test_effectsize}
\alias{print.effectsize_table}
\alias{format.effectsize_table}
\alias{print.effectsize_difference}
\title{Methods for \code{effectsize} tables}
\usage{
\method{plot}{effectsize_table}(x, ...)

\method{plot}{equivalence_test_effectsize}(x, ...)

\method{print}{effectsize_table}(x, digits = 2, ...)

\method{format}{effectsize_table}(x, digits = 2, ...)

\method{print}{effectsize_difference}(x, digits = 2, append_CLES = FALSE, ...)
}
\arguments{
\item{x}{Object to print.}

\item{...}{Arguments passed to or from other functions.}

\item{digits}{Number of digits for rounding or significant figures. May also
be \code{"signif"} to return significant figures or \code{"scientific"}
to return scientific notation. Control the number of digits by adding the
value as suffix, e.g. \code{digits = "scientific4"} to have scientific
notation with 4 decimal places, or \code{digits = "signif5"} for 5
significant figures (see also \code{\link[=signif]{signif()}}).}

\item{append_CLES}{Should the Common Language Effect Sizes be printed as well?
Only applicable to Cohen's \emph{d}, Hedges' \emph{g} for independent samples of
equal variance (pooled sd) or for the rank-biserial correlation for
independent samples (See \code{\link[=d_to_cles]{d_to_cles()}})}
}
\description{
Printing, formatting and plotting methods for \code{effectsize} tables.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_kendalls_w.R
\name{interpret_kendalls_w}
\alias{interpret_kendalls_w}
\title{Interpret Kendall's coefficient of concordance}
\usage{
interpret_kendalls_w(w, rules = "landis1977")
}
\arguments{
\item{w}{Value or vector of Kendall's coefficient of concordance.}

\item{rules}{Can be \code{"landis1977"} (default) or a custom set of \code{\link[=rules]{rules()}}.}
}
\description{
Interpret Kendall's coefficient of concordance
}
\section{Rules}{

\itemize{
\item Landis & Koch (1977) (\code{"landis1977"}; default)
\itemize{
\item \strong{0.00 <= w < 0.20} - Slight agreement
\item \strong{0.20 <= w < 0.40} - Fair agreement
\item \strong{0.40 <= w < 0.60} - Moderate agreement
\item \strong{0.60 <= w < 0.80} - Substantial agreement
\item \strong{w >= 0.80}        - Almost perfect agreement
}
}
}

\references{
\itemize{
\item Landis, J. R., & Koch G. G. (1977). The measurement of observer agreement
for categorical data. Biometrics, 33:159-74.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_stat_to_anova.R
\name{F_to_eta2}
\alias{F_to_eta2}
\alias{t_to_eta2}
\alias{F_to_epsilon2}
\alias{t_to_epsilon2}
\alias{F_to_eta2_adj}
\alias{t_to_eta2_adj}
\alias{F_to_omega2}
\alias{t_to_omega2}
\alias{F_to_f}
\alias{t_to_f}
\alias{F_to_f2}
\alias{t_to_f2}
\title{Convert test statistics (F, t) to indices of \strong{partial} variance explained
(\strong{partial} Eta / Omega / Epsilon squared and Cohen's f)}
\usage{
F_to_eta2(f, df, df_error, ci = 0.95, alternative = "greater", ...)

t_to_eta2(t, df_error, ci = 0.95, alternative = "greater", ...)

F_to_epsilon2(f, df, df_error, ci = 0.95, alternative = "greater", ...)

t_to_epsilon2(t, df_error, ci = 0.95, alternative = "greater", ...)

F_to_eta2_adj(f, df, df_error, ci = 0.95, alternative = "greater", ...)

t_to_eta2_adj(t, df_error, ci = 0.95, alternative = "greater", ...)

F_to_omega2(f, df, df_error, ci = 0.95, alternative = "greater", ...)

t_to_omega2(t, df_error, ci = 0.95, alternative = "greater", ...)

F_to_f(
  f,
  df,
  df_error,
  ci = 0.95,
  alternative = "greater",
  squared = FALSE,
  ...
)

t_to_f(t, df_error, ci = 0.95, alternative = "greater", squared = FALSE, ...)

F_to_f2(
  f,
  df,
  df_error,
  ci = 0.95,
  alternative = "greater",
  squared = TRUE,
  ...
)

t_to_f2(t, df_error, ci = 0.95, alternative = "greater", squared = TRUE, ...)
}
\arguments{
\item{df, df_error}{Degrees of freedom of numerator or of the error estimate
(i.e., the residuals).}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"greater"} (default) or \code{"less"}
(one-sided CI), or \code{"two.sided"} (default, two-sided CI). Partial matching
is allowed (e.g., \code{"g"}, \code{"l"}, \code{"two"}...). See \emph{One-Sided CIs} in
\link{effectsize_CIs}.}

\item{...}{Arguments passed to or from other methods.}

\item{t, f}{The t or the F statistics.}

\item{squared}{Return Cohen's \emph{f} or Cohen's \emph{f}-squared?}
}
\value{
A data frame with the effect size(s) between 0-1 (\code{Eta2_partial},
\code{Epsilon2_partial}, \code{Omega2_partial}, \code{Cohens_f_partial} or
\code{Cohens_f2_partial}), and their CIs (\code{CI_low} and \code{CI_high}). (Note that
for \eqn{\omega_p^2}{Omega2p} and \eqn{\epsilon_p^2}{Epsilon2p} it is possible to compute a
negative number; even though this doesn't make any practical sense, it is
recommended to report the negative number and not a 0).
}
\description{
These functions are convenience functions to convert F and t test statistics
to \strong{partial} Eta- (\eqn{\eta}), Omega- (\eqn{\omega}) Epsilon-
(\eqn{\epsilon}) squared (an alias for the adjusted Eta squared) and Cohen's
f. These are useful in cases where the various Sum of Squares and Mean
Squares are not easily available or their computation is not straightforward
(e.g., in liner mixed models, contrasts, etc.). For test statistics derived
from \code{lm} and \code{aov} models, these functions give exact results. For all other
cases, they return close approximations.
\cr
See \href{https://easystats.github.io/effectsize/articles/from_test_statistics.html}{Effect Size from Test Statistics vignette.}
}
\details{
These functions use the following formulae:
\cr
\deqn{\eta_p^2 = \frac{F \times df_{num}}{F \times df_{num} + df_{den}}}{\eta^2_p = F * df1 / (F * df1 + df2)}
\cr
\deqn{\epsilon_p^2 = \frac{(F - 1) \times df_{num}}{F \times df_{num} + df_{den}}}{\epsilon^2_p = (F - 1) * df1 / (F * df1 + df2)}
\cr
\deqn{\omega_p^2 = \frac{(F - 1) \times df_{num}}{F \times df_{num} + df_{den} + 1}}{\omega^2_p=(F - 1) * df1 / (F * df1 + df2 + 1)}
\cr
\deqn{f_p = \sqrt{\frac{\eta_p^2}{1-\eta_p^2}}}{f = \eta^2 / (1 - \eta^2)}
\cr\cr
For \emph{t}, the conversion is based on the equality of \eqn{t^2 = F} when \eqn{df_{num}=1}{df1 = 1}.
\subsection{Choosing an Un-Biased Estimate}{

Both Omega and Epsilon are unbiased estimators of the population Eta. But
which to choose? Though Omega is the more popular choice, it should be noted
that:
\enumerate{
\item The formula given above for Omega is only an approximation for complex
designs.
\item Epsilon has been found to be less biased (Carroll & Nordholm, 1975).
}
}
}
\note{
Adjusted (partial) Eta-squared is an alias for (partial) Epsilon-squared.
}
\section{Confidence (Compatibility) Intervals (CIs)}{

Unless stated otherwise, confidence (compatibility) intervals (CIs) are
estimated using the noncentrality parameter method (also called the "pivot
method"). This method finds the noncentrality parameter ("\emph{ncp}") of a
noncentral \emph{t}, \emph{F}, or \eqn{\chi^2} distribution that places the observed
\emph{t}, \emph{F}, or \eqn{\chi^2} test statistic at the desired probability point of
the distribution. For example, if the observed \emph{t} statistic is 2.0, with 50
degrees of freedom, for which cumulative noncentral \emph{t} distribution is \emph{t} =
2.0 the .025 quantile (answer: the noncentral \emph{t} distribution with \emph{ncp} =
.04)? After estimating these confidence bounds on the \emph{ncp}, they are
converted into the effect size metric to obtain a confidence interval for the
effect size (Steiger, 2004).
\cr\cr
For additional details on estimation and troubleshooting, see \link{effectsize_CIs}.
}

\section{CIs and Significance Tests}{

"Confidence intervals on measures of effect size convey all the information
in a hypothesis test, and more." (Steiger, 2004). Confidence (compatibility)
intervals and p values are complementary summaries of parameter uncertainty
given the observed data. A dichotomous hypothesis test could be performed
with either a CI or a p value. The 100 (1 - \eqn{\alpha})\% confidence
interval contains all of the parameter values for which \emph{p} > \eqn{\alpha}
for the current data and model. For example, a 95\% confidence interval
contains all of the values for which p > .05.
\cr\cr
Note that a confidence interval including 0 \emph{does not} indicate that the null
(no effect) is true. Rather, it suggests that the observed data together with
the model and its assumptions combined do not provided clear evidence against
a parameter value of 0 (same as with any other value in the interval), with
the level of this evidence defined by the chosen \eqn{\alpha} level (Rafi &
Greenland, 2020; Schweder & Hjort, 2016; Xie & Singh, 2013). To infer no
effect, additional judgments about what parameter values are "close enough"
to 0 to be negligible are needed ("equivalence testing"; Bauer & Kiesser,
1996).
}

\examples{
\donttest{
if (require("afex")) {
  data(md_12.1)
  aov_ez("id", "rt", md_12.1,
    within = c("angle", "noise"),
    anova_table = list(correction = "none", es = "pes")
  )
}
# compare to:
(etas <- F_to_eta2(
  f = c(40.72, 33.77, 45.31),
  df = c(2, 1, 2),
  df_error = c(18, 9, 18)
))

if (require(see)) plot(etas)


if (require("lmerTest")) { # for the df_error
  fit <- lmer(extra ~ group + (1 | ID), sleep)
  # anova(fit)
  # #> Type III Analysis of Variance Table with Satterthwaite's method
  # #>       Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)
  # #> group 12.482  12.482     1     9  16.501 0.002833 **
  # #> ---
  # #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

  F_to_eta2(16.501, 1, 9)
  F_to_omega2(16.501, 1, 9)
  F_to_epsilon2(16.501, 1, 9)
  F_to_f(16.501, 1, 9)
}


## Use with emmeans based contrasts
## --------------------------------
if (require(emmeans)) {
  warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)

  jt <- joint_tests(warp.lm, by = "wool")
  F_to_eta2(jt$F.ratio, jt$df1, jt$df2)
}
}
}
\references{
\itemize{
\item Albers, C., & Lakens, D. (2018). When power analyses based on pilot data
are biased: Inaccurate effect size estimators and follow-up bias. Journal of
experimental social psychology, 74, 187-195. \doi{10.31234/osf.io/b7z4q}
\item Carroll, R. M., & Nordholm, L. A. (1975). Sampling Characteristics of
Kelley's epsilon and Hays' omega. Educational and Psychological Measurement,
35(3), 541-554.
\item Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and noncentral
distributions. Educational and Psychological Measurement, 61(4), 532-574.
\item Friedman, H. (1982). Simplified determinations of statistical power,
magnitude of effect and research sample sizes. Educational and Psychological
Measurement, 42(2), 521-526. \doi{10.1177/001316448204200214}
\item Mordkoff, J. T. (2019). A Simple Method for Removing Bias From a Popular
Measure of Standardized Effect Size: Adjusted Partial Eta Squared. Advances
in Methods and Practices in Psychological Science, 2(3), 228-232.
\doi{10.1177/2515245919855053}
\item Morey, R. D., Hoekstra, R., Rouder, J. N., Lee, M. D., & Wagenmakers, E. J.
(2016). The fallacy of placing confidence in confidence intervals.
Psychonomic bulletin & review, 23(1), 103-123.
\item Steiger, J. H. (2004). Beyond the F test: Effect size confidence intervals
and tests of close fit in the analysis of variance and contrast analysis.
Psychological Methods, 9, 164-182.
}
}
\seealso{
\code{\link[=eta_squared]{eta_squared()}} for more details.

Other effect size from test statistic: 
\code{\link{chisq_to_phi}()},
\code{\link{t_to_d}()}
}
\concept{effect size from test statistic}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effectsize.BFBayesFactor.R, R/effectsize.R,
%   R/effectsize.htest.R
\name{effectsize.BFBayesFactor}
\alias{effectsize.BFBayesFactor}
\alias{effectsize}
\alias{effectsize.aov}
\alias{effectsize.htest}
\title{Effect Size}
\usage{
\method{effectsize}{BFBayesFactor}(model, type = NULL, verbose = TRUE, test = NULL, ...)

effectsize(model, ...)

\method{effectsize}{aov}(model, type = NULL, ...)

\method{effectsize}{htest}(model, type = NULL, verbose = TRUE, ...)
}
\arguments{
\item{model}{An object of class \code{htest}, or a statistical model. See details.}

\item{type}{The effect size of interest. See details.}

\item{verbose}{Toggle off warnings.}

\item{test}{The indices of effect existence to compute. Character (vector) or
list with one or more of these options: \code{"p_direction"} (or \code{"pd"}),
\code{"rope"}, \code{"p_map"}, \code{"equivalence_test"} (or \code{"equitest"}),
\code{"bayesfactor"} (or \code{"bf"}) or \code{"all"} to compute all tests.
For each "test", the corresponding \pkg{bayestestR} function is called
(e.g. \code{\link[bayestestR:rope]{rope()}} or \code{\link[bayestestR:p_direction]{p_direction()}}) and its results
included in the summary output.}

\item{...}{Arguments passed to or from other methods. See details.}
}
\value{
A data frame with the effect size (depending on input) and and its
CIs (\code{CI_low} and \code{CI_high}).
}
\description{
This function tries to return the best effect-size measure for the provided
input model. See details.
}
\details{
\itemize{
\item For an object of class \code{htest}, data is extracted via \code{\link[insight:get_data]{insight::get_data()}}, and passed to the relevant function according to:
\itemize{
\item A \strong{t-test} depending on \code{type}: \code{"cohens_d"} (default), \code{"hedges_g"}, or \code{"cles"}.
\item A \strong{Chi-squared tests of independence or goodness-of-fit}, depending on \code{type}: \code{"cramers_v"} (default), \code{"phi"}, \code{"cohens_w"}, \code{"pearsons_c"}, \code{"cohens_h"}, \code{"oddsratio"}, or \code{"riskratio"}.
\item A \strong{One-way ANOVA test}, depending on \code{type}: \code{"eta"} (default), \code{"omega"} or \code{"epsilon"} -squared, \code{"f"}, or \code{"f2"}.
\item A \strong{McNemar test} returns \emph{Cohen's g}.
\item A \strong{Wilcoxon test} depending on \code{type}: returns "\code{rank_biserial}" correlation (default) or \code{"cles"}.
\item A \strong{Kruskal-Wallis test} returns \emph{rank Epsilon squared}.
\item A \strong{Friedman test} returns \emph{Kendall's W}.
(Where applicable, \code{ci} and \code{alternative} are taken from the \code{htest} if not otherwise provided.)
}
\item For an object of class \code{BFBayesFactor}, using \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}},
\itemize{
\item A \strong{t-test} depending on \code{type}: "cohens_d"\verb{(default) or}"cles"`.
\item A \strong{correlation test} returns \emph{r}.
\item A \strong{contingency table test}, depending on \code{type}: \code{"cramers_v"} (default), \code{"phi"}, \code{"cohens_w"}, \code{"pearsons_c"}, \code{"cohens_h"}, \code{"oddsratio"}, or \code{"riskratio"}.
\item A \strong{proportion test} returns \emph{p}.
}
\item Objects of class \code{anova}, \code{aov}, or \code{aovlist}, depending on \code{type}: \code{"eta"} (default), \code{"omega"} or \code{"epsilon"} -squared, \code{"f"}, or \code{"f2"}.
\item Other objects are passed to \code{\link[=standardize_parameters]{standardize_parameters()}}.
}

\strong{For statistical models it is recommended to directly use the listed
functions, for the full range of options they provide.}
}
\examples{

## Hypothesis Testing
## ------------------
contingency_table <- as.table(rbind(c(762, 327, 468), c(484, 239, 477), c(484, 239, 477)))
Xsq <- chisq.test(contingency_table)
effectsize(Xsq)
effectsize(Xsq, type = "phi")

Tt <- t.test(1:10, y = c(7:20), alternative = "less")
effectsize(Tt)

Aov <- oneway.test(extra ~ group, data = sleep, var.equal = TRUE)
effectsize(Aov)
effectsize(Aov, type = "omega")

Wt <- wilcox.test(1:10, 7:20, mu = -3, alternative = "less")
effectsize(Wt)
effectsize(Wt, type = "cles")

## Bayesian Hypothesis Testing
## ---------------------------
\donttest{
if (require(BayesFactor)) {
  bf_prop <- proportionBF(3, 7, p = 0.3)
  effectsize(bf_prop)

  bf_corr <- correlationBF(attitude$rating, attitude$complaints)
  effectsize(bf_corr)

  data(raceDolls)
  bf_xtab <- contingencyTableBF(raceDolls, sampleType = "poisson", fixedMargin = "cols")
  effectsize(bf_xtab)
  effectsize(bf_xtab, type = "oddsratio")

  bf_ttest <- ttestBF(sleep$extra[sleep$group==1],
                      sleep$extra[sleep$group==2],
                      paired = TRUE, mu = -1)
  effectsize(bf_ttest)
}
}

## Models and Anova Tables
## -----------------------
fit <- lm(mpg ~ factor(cyl) * wt + hp, data = mtcars)
effectsize(fit)

anova_table <- anova(fit)
effectsize(anova_table)
effectsize(anova_table, type = "epsilon")
}
\seealso{
Other effect size indices: 
\code{\link{cles}()},
\code{\link{cohens_d}()},
\code{\link{eta_squared}()},
\code{\link{phi}()},
\code{\link{rank_biserial}()},
\code{\link{standardize_parameters}()}
}
\concept{effect size indices}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sd_pooled.R
\name{sd_pooled}
\alias{sd_pooled}
\alias{mad_pooled}
\title{Pooled Standard Deviation}
\usage{
sd_pooled(x, y = NULL, data = NULL, verbose = TRUE)

mad_pooled(x, y = NULL, data = NULL, constant = 1.4826, verbose = TRUE)
}
\arguments{
\item{x}{A formula, a numeric vector, or a character name of one in \code{data}.}

\item{y}{A numeric vector, a grouping (character / factor) vector, a or a
character  name of one in \code{data}. Ignored if \code{x} is a formula.}

\item{data}{An optional data frame containing the variables.}

\item{verbose}{Toggle warnings and messages on or off.}

\item{constant}{scale factor.}
}
\value{
Numeric, the pooled standard deviation.
}
\description{
The Pooled Standard Deviation is a weighted average of standard deviations
for two or more groups, \emph{assumed to have equal variance}. It represents the
common deviation among the groups, around each of their respective means.
}
\details{
The standard version is calculated as:
\deqn{\sqrt{\frac{\sum (x_i - \bar{x})^2}{n_1 + n_2 - 2}}}{sqrt(sum(c(x - mean(x), y - mean(y))^2) / (n1 + n2 - 2))}
The robust version is calculated as:
\deqn{1.4826 \times Median(|\left\{x - Median_x,\,y - Median_y\right\}|)}{mad(c(x - median(x), y - median(y)), constant = 1.4826)}
}
\examples{
sd_pooled(mpg ~ am, data = mtcars)
mad_pooled(mtcars$mpg, factor(mtcars$am))
}
\seealso{
\code{\link[=cohens_d]{cohens_d()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_standardize.R
\name{format_standardize}
\alias{format_standardize}
\title{Transform a standardized vector into character}
\usage{
format_standardize(
  x,
  reference = x,
  robust = FALSE,
  digits = 1,
  protect_integers = TRUE,
  ...
)
}
\arguments{
\item{x}{A standardized numeric vector.}

\item{reference}{The reference vector from which to compute the mean and SD.}

\item{robust}{Logical, if \code{TRUE}, centering is done by subtracting the
median from the variables and dividing it by the median absolute deviation
(MAD). If \code{FALSE}, variables are standardized by subtracting the
mean and dividing it by the standard deviation (SD).}

\item{digits}{Number of digits for rounding or significant figures. May also
be \code{"signif"} to return significant figures or \code{"scientific"}
to return scientific notation. Control the number of digits by adding the
value as suffix, e.g. \code{digits = "scientific4"} to have scientific
notation with 4 decimal places, or \code{digits = "signif5"} for 5
significant figures (see also \code{\link[=signif]{signif()}}).}

\item{protect_integers}{Should integers be kept as integers (i.e., without
decimals)?}

\item{...}{Other arguments to pass to \code{\link[insight:format_value]{insight::format_value()}} such as \code{digits}, etc.}
}
\description{
Transform a standardized vector into character, e.g., \code{c("-1 SD", "Mean", "+1 SD")}.
}
\examples{
format_standardize(c(-1, 0, 1))
format_standardize(c(-1, 0, 1, 2), reference = rnorm(1000))
format_standardize(c(-1, 0, 1, 2), reference = rnorm(1000), robust = TRUE)

format_standardize(standardize(mtcars$wt), digits = 1)
format_standardize(standardize(mtcars$wt, robust = TRUE), digits = 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eta_squared.R, R/eta_squared_posterior.R
\name{eta_squared}
\alias{eta_squared}
\alias{omega_squared}
\alias{epsilon_squared}
\alias{cohens_f}
\alias{cohens_f_squared}
\alias{eta_squared_posterior}
\title{Effect size for ANOVA}
\usage{
eta_squared(
  model,
  partial = TRUE,
  generalized = FALSE,
  ci = 0.95,
  alternative = "greater",
  verbose = TRUE,
  ...
)

omega_squared(
  model,
  partial = TRUE,
  ci = 0.95,
  alternative = "greater",
  verbose = TRUE,
  ...
)

epsilon_squared(
  model,
  partial = TRUE,
  ci = 0.95,
  alternative = "greater",
  verbose = TRUE,
  ...
)

cohens_f(
  model,
  partial = TRUE,
  ci = 0.95,
  alternative = "greater",
  squared = FALSE,
  verbose = TRUE,
  model2 = NULL,
  ...
)

cohens_f_squared(
  model,
  partial = TRUE,
  ci = 0.95,
  alternative = "greater",
  squared = TRUE,
  verbose = TRUE,
  model2 = NULL,
  ...
)

eta_squared_posterior(
  model,
  partial = TRUE,
  generalized = FALSE,
  ss_function = stats::anova,
  draws = 500,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{A model, ANOVA object, or the result of \code{parameters::model_parameters}.}

\item{partial}{If \code{TRUE}, return partial indices.}

\item{generalized}{If TRUE, returns generalized Eta Squared, assuming all
variables are manipulated. Can also be a character vector of observed
(non-manipulated) variables, in which case generalized Eta Squared is
calculated taking these observed variables into account. For \code{afex_aov}
model, when \code{generalized = TRUE}, the observed variables are extracted
automatically from the fitted model, if they were provided then.}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"greater"} (default) or \code{"less"}
(one-sided CI), or \code{"two.sided"} (default, two-sided CI). Partial matching
is allowed (e.g., \code{"g"}, \code{"l"}, \code{"two"}...). See \emph{One-Sided CIs} in
\link{effectsize_CIs}.}

\item{verbose}{Toggle warnings and messages on or off.}

\item{...}{Arguments passed to or from other methods.
\itemize{
\item Can be \code{include_intercept = TRUE} to include the effect size for the intercept.
\item For Bayesian models, arguments passed to \code{ss_function}.
}}

\item{squared}{Return Cohen's \emph{f} or Cohen's \emph{f}-squared?}

\item{model2}{Optional second model for Cohen's f (/squared). If specified,
returns the effect size for R-squared-change between the two models.}

\item{ss_function}{For Bayesian models, the function used to extract
sum-of-squares. Uses \code{\link[=anova]{anova()}} by default, but can also be \code{car::Anova()}
for simple linear models.}

\item{draws}{For Bayesian models, an integer indicating the number of draws
from the posterior predictive distribution to return. Larger numbers take
longer to run, but provide estimates that are more stable.}
}
\value{
A data frame with the effect size(s) between 0-1 (\code{Eta2}, \code{Epsilon2},
\code{Omega2}, \code{Cohens_f} or \code{Cohens_f2}, possibly with the \code{partial} or
\code{generalized} suffix), and their CIs (\code{CI_low} and \code{CI_high}).
\cr\cr
For \code{eta_squared_posterior()}, a data frame containing the ppd of the Eta
squared for each fixed effect, which can then be passed to
\code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}} for summary stats.

A data frame containing the effect size values and their confidence
intervals.
}
\description{
Functions to compute effect size measures for ANOVAs, such as Eta-
(\eqn{\eta}), Omega- (\eqn{\omega}) and Epsilon- (\eqn{\epsilon}) squared,
and Cohen's f (or their partialled versions) for ANOVA tables. These indices
represent an estimate of how much variance in the response variables is
accounted for by the explanatory variable(s).
\cr\cr
When passing models, effect sizes are computed using the sums of squares
obtained from \code{anova(model)} which might not always be appropriate. See
details.
}
\details{
For \code{aov}, \code{aovlist} and \code{afex_aov} models, and for \code{anova} objects that
provide Sums-of-Squares, the effect sizes are computed directly using
Sums-of-Squares (for \code{mlm} / \code{maov} models, effect sizes are computed for
each response separately). For all other model, effect sizes are approximated
via test statistic conversion of the omnibus \emph{F} statistic provided by the
appropriate \code{anova()} method (see \code{\link[=F_to_eta2]{F_to_eta2()}} for more details.)
\subsection{Type of Sums of Squares}{

The sums of squares (or \emph{F} statistics) used for the computation of the
effect sizes is based on those returned by \code{anova(model)} (whatever those may
be - for \code{aov} and \code{aovlist} these are \emph{type-1} sums of squares; for
\code{lmerMod} (and \code{lmerModLmerTest}) these are \emph{type-3} sums of squares). Make
sure these are the sums of squares you are interested in; You might want to
pass the result of \code{car::Anova(mode, type = 2)} or \code{type = 3} instead of the
model itself, or use the \code{afex} package to fit ANOVA models.
\cr\cr
For type 3 sum of squares, it is generally recommended to fit models with
\emph{\code{contr.sum} factor weights} and \emph{centered covariates}, for sensible results.
See examples and the \code{afex} package.
}

\subsection{Un-Biased Estimate of Eta}{

Both \emph{\strong{Omega}} and \emph{\strong{Epsilon}} are unbiased estimators of the
population's \emph{\strong{Eta}}, which is especially important is small samples. But
which to choose?
\cr\cr
Though Omega is the more popular choice (Albers \& Lakens, 2018), Epsilon is
analogous to adjusted R2 (Allen, 2017, p. 382), and has been found to be less
biased (Carroll & Nordholm, 1975).
\cr\cr
(Note that for Omega- and Epsilon-squared it is possible to compute a
negative number; even though this doesn't make any practical sense, it is
recommended to report the negative number and not a 0.)
}

\subsection{Cohen's f}{

Cohen's f can take on values between zero, when the population means are all
equal, and an indefinitely large number as standard deviation of means
increases relative to the average standard deviation within each group.
\cr\cr
When comparing two models in a sequential regression analysis, Cohen's f for
R-square change is the ratio between the increase in R-square
and the percent of unexplained variance.
\cr\cr
Cohen has suggested that the values of 0.10, 0.25, and 0.40 represent small,
medium, and large effect sizes, respectively.
}

\subsection{Eta Squared from Posterior Predictive Distribution}{

For Bayesian models (fit with \code{brms} or \code{rstanarm}),
\code{eta_squared_posterior()} simulates data from the posterior predictive
distribution (ppd) and for each simulation the Eta Squared is computed for
the model's fixed effects. This means that the returned values are the
population level effect size as implied by the posterior model (and not the
effect size in the sample data). See \code{\link[rstantools:posterior_predict]{rstantools::posterior_predict()}} for
more info.
}
}
\section{Confidence (Compatibility) Intervals (CIs)}{

Unless stated otherwise, confidence (compatibility) intervals (CIs) are
estimated using the noncentrality parameter method (also called the "pivot
method"). This method finds the noncentrality parameter ("\emph{ncp}") of a
noncentral \emph{t}, \emph{F}, or \eqn{\chi^2} distribution that places the observed
\emph{t}, \emph{F}, or \eqn{\chi^2} test statistic at the desired probability point of
the distribution. For example, if the observed \emph{t} statistic is 2.0, with 50
degrees of freedom, for which cumulative noncentral \emph{t} distribution is \emph{t} =
2.0 the .025 quantile (answer: the noncentral \emph{t} distribution with \emph{ncp} =
.04)? After estimating these confidence bounds on the \emph{ncp}, they are
converted into the effect size metric to obtain a confidence interval for the
effect size (Steiger, 2004).
\cr\cr
For additional details on estimation and troubleshooting, see \link{effectsize_CIs}.
}

\section{CIs and Significance Tests}{

"Confidence intervals on measures of effect size convey all the information
in a hypothesis test, and more." (Steiger, 2004). Confidence (compatibility)
intervals and p values are complementary summaries of parameter uncertainty
given the observed data. A dichotomous hypothesis test could be performed
with either a CI or a p value. The 100 (1 - \eqn{\alpha})\% confidence
interval contains all of the parameter values for which \emph{p} > \eqn{\alpha}
for the current data and model. For example, a 95\% confidence interval
contains all of the values for which p > .05.
\cr\cr
Note that a confidence interval including 0 \emph{does not} indicate that the null
(no effect) is true. Rather, it suggests that the observed data together with
the model and its assumptions combined do not provided clear evidence against
a parameter value of 0 (same as with any other value in the interval), with
the level of this evidence defined by the chosen \eqn{\alpha} level (Rafi &
Greenland, 2020; Schweder & Hjort, 2016; Xie & Singh, 2013). To infer no
effect, additional judgments about what parameter values are "close enough"
to 0 to be negligible are needed ("equivalence testing"; Bauer & Kiesser,
1996).
}

\examples{
\donttest{
data(mtcars)
mtcars$am_f <- factor(mtcars$am)
mtcars$cyl_f <- factor(mtcars$cyl)

model <- aov(mpg ~ am_f * cyl_f, data = mtcars)

(eta2 <- eta_squared(model))

# More types:
eta_squared(model, partial = FALSE)
eta_squared(model, generalized = "cyl_f")
omega_squared(model)
epsilon_squared(model)
cohens_f(model)

if (require(see)) plot(eta2)

model0 <- aov(mpg ~ am_f + cyl_f, data = mtcars) # no interaction
cohens_f_squared(model0, model2 = model)

## Interpretation of effect sizes
## -------------------------------------

interpret_omega_squared(0.10, rules = "field2013")
interpret_eta_squared(0.10, rules = "cohen1992")
interpret_epsilon_squared(0.10, rules = "cohen1992")

interpret(eta2, rules = "cohen1992")

# Recommended: Type-3 effect sizes + effects coding
# -------------------------------------------------
if (require(car, quietly = TRUE)) {
  contrasts(mtcars$am_f) <- contr.sum
  contrasts(mtcars$cyl_f) <- contr.sum

  model <- aov(mpg ~ am_f * cyl_f, data = mtcars)
  model_anova <- car::Anova(model, type = 3)

  eta_squared(model_anova)
}

# afex takes care of both type-3 effects and effects coding:
if (require(afex)) {
  data(obk.long, package = "afex")
  model <- aov_car(value ~ treatment * gender + Error(id / (phase)),
    data = obk.long, observed = "gender"
  )
  eta_squared(model)
  epsilon_squared(model)
  omega_squared(model)
  eta_squared(model, partial = FALSE)
  epsilon_squared(model, partial = FALSE)
  omega_squared(model, partial = FALSE)
  eta_squared(model, generalized = TRUE) # observed vars are pulled from the afex model.
}



## Approx. effect sizes for mixed models
## -------------------------------------
if (require(lmerTest, quietly = TRUE)) {
  model <- lmer(mpg ~ am_f * cyl_f + (1 | vs), data = mtcars)
  omega_squared(model)
}




## Bayesian Models (PPD)
## ---------------------
\dontrun{
if (require(rstanarm) && require(bayestestR) && require(car)) {
  fit_bayes <- stan_glm(mpg ~ factor(cyl) * wt + qsec,
    data = mtcars,
    family = gaussian(),
    refresh = 0
  )

  es <- eta_squared_posterior(fit_bayes,
    ss_function = car::Anova, type = 3
  )
  bayestestR::describe_posterior(es)


  # compare to:
  fit_freq <- lm(mpg ~ factor(cyl) * wt + qsec,
    data = mtcars
  )
  aov_table <- car::Anova(fit_freq, type = 3)
  eta_squared(aov_table)
}
}
}

}
\references{
\itemize{
\item Albers, C., \& Lakens, D. (2018). When power analyses based on pilot data
are biased: Inaccurate effect size estimators and follow-up bias. Journal of
experimental social psychology, 74, 187-195.
\item Allen, R. (2017). Statistics and Experimental Design for Psychologists: A
Model Comparison Approach. World Scientific Publishing Company.
\item Carroll, R. M., & Nordholm, L. A. (1975). Sampling Characteristics of
Kelley's epsilon and Hays' omega. Educational and Psychological Measurement,
35(3), 541-554.
\item Kelley, T. (1935) An unbiased correlation ratio measure. Proceedings of the
National Academy of Sciences. 21(9). 554-559.
\item Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: measures of effect size for some common research designs.
Psychological methods, 8(4), 434.
\item Steiger, J. H. (2004). Beyond the F test: Effect size confidence intervals
and tests of close fit in the analysis of variance and contrast analysis.
Psychological Methods, 9, 164-182.
}
}
\seealso{
\code{\link[=F_to_eta2]{F_to_eta2()}}

Other effect size indices: 
\code{\link{cles}()},
\code{\link{cohens_d}()},
\code{\link{effectsize.BFBayesFactor}()},
\code{\link{phi}()},
\code{\link{rank_biserial}()},
\code{\link{standardize_parameters}()}
}
\concept{effect size indices}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz_deprecated.R
\name{effectsize_deprecated}
\alias{effectsize_deprecated}
\alias{interpret_d}
\alias{interpret_g}
\alias{interpret_delta}
\alias{interpret_parameters}
\title{Deprecated functions}
\usage{
interpret_d(...)

interpret_g(...)

interpret_delta(...)

interpret_parameters(...)
}
\arguments{
\item{...}{Arguments to the deprecated function.}
}
\description{
Deprecated functions
}
\details{
\itemize{
\item \code{interpret_d} is now \code{\link{interpret_cohens_d}}.
\item \code{interpret_g} is now \code{\link{interpret_hedges_g}}.
\item \code{interpret_delta} is now \code{\link{interpret_glass_delta}}.
\item \code{interpret_parameters} for \emph{standardized parameters} was incorrect. Use \code{\link{interpret_r}} instead.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret.R
\name{interpret}
\alias{interpret}
\alias{interpret.numeric}
\alias{interpret.effectsize_table}
\title{Generic function for interpretation}
\usage{
interpret(x, ...)

\method{interpret}{numeric}(x, rules, name = attr(rules, "rule_name"), ...)

\method{interpret}{effectsize_table}(x, rules, ...)
}
\arguments{
\item{x}{Vector of value break points (edges defining categories), or a data
frame of class \code{effectsize_table}.}

\item{...}{Currently not used.}

\item{rules}{Set of \code{\link[=rules]{rules()}}. When \code{x} is a data frame, can be a name of an
established set of rules.}

\item{name}{Name of the set of rules (will be printed).}
}
\value{
\itemize{
\item For numeric input: A character vector of interpertations.
\item For data frames: the \code{x} input with an additional \code{Interpretation} column.
}
}
\description{
Interpret a value based on a set of rules. See \code{\link[=rules]{rules()}}.
}
\examples{
rules_grid <- rules(c(0.01, 0.05), c("very significant", "significant", "not significant"))
interpret(0.001, rules_grid)
interpret(0.021, rules_grid)
interpret(0.08, rules_grid)
interpret(c(0.01, 0.005, 0.08), rules_grid)

interpret(c(0.35, 0.15), c("small" = 0.2, "large" = 0.4), name = "Cohen's Rules")
interpret(c(0.35, 0.15), rules(c(0.2, 0.4), c("small", "medium", "large")))

# ----------
d <- cohens_d(mpg ~ am, data = mtcars)
interpret(d, rules = "cohen1988")

d <- glass_delta(mpg ~ am, data = mtcars)
interpret(d, rules = "gignac2016")

interpret(d, rules = rules(1, c("tiny", "yeah okay")))

m <- lm(formula = wt ~ am * cyl, data = mtcars)
eta2 <- eta_squared(m)
interpret(eta2, rules = "field2013")

X <- chisq.test(mtcars$am, mtcars$cyl == 8)
interpret(oddsratio(X), rules = "chen2010")
interpret(cramers_v(X), "lovakov2021")
}
\seealso{
rules
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_ess_rhat.R
\name{interpret_ess}
\alias{interpret_ess}
\alias{interpret_rhat}
\title{Interpret Bayesian diagnostic indices}
\usage{
interpret_ess(ess, rules = "burkner2017")

interpret_rhat(rhat, rules = "vehtari2019")
}
\arguments{
\item{ess}{Value or vector of Effective Sample Size (ESS) values.}

\item{rules}{A character string (see \emph{Rules}) or a custom set of \code{\link[=rules]{rules()}}.}

\item{rhat}{Value or vector of Rhat values.}
}
\description{
Interpretation of Bayesian diagnostic indices, such as Effective Sample Size (ESS) and Rhat.
}
\section{Rules}{

\subsection{ESS}{
\itemize{
\item Bürkner, P. C. (2017) (\code{"burkner2017"}; default)
\itemize{
\item \strong{ESS < 1000} - Insufficient
\item \strong{ESS >= 1000} - Sufficient
}
}
}

\subsection{Rhat}{
\itemize{
\item Vehtari et al. (2019) (\code{"vehtari2019"}; default)
\itemize{
\item \strong{Rhat < 1.01} - Converged
\item \strong{Rhat >= 1.01} - Failed
}
\item Gelman & Rubin (1992) (\code{"gelman1992"})
\itemize{
\item \strong{Rhat < 1.1} - Converged
\item \strong{Rhat >= 1.1} - Failed
}
}
}
}

\examples{
interpret_ess(1001)
interpret_ess(c(852, 1200))

interpret_rhat(1.00)
interpret_rhat(c(1.5, 0.9))
}
\references{
\itemize{
\item Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models
using Stan. Journal of Statistical Software, 80(1), 1-28.
\item Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
using multiple sequences. Statistical science, 7(4), 457-472.
\item Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P. C.
(2019). Rank-normalization, folding, and localization: An improved Rhat for
assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_icc.R
\name{interpret_icc}
\alias{interpret_icc}
\title{Interpret Intraclass Correlation Coefficient (ICC)}
\usage{
interpret_icc(icc, rules = "koo2016", ...)
}
\arguments{
\item{icc}{Value or vector of Intraclass Correlation Coefficient (ICC) values.}

\item{rules}{Can be \code{"koo2016"} (default) or custom set of \code{\link[=rules]{rules()}}.}

\item{...}{Not used for now.}
}
\description{
The value of an ICC lies between 0 to 1, with 0 indicating no reliability among raters and 1 indicating perfect reliability.
}
\section{Rules}{

\itemize{
\item Koo (2016) (\code{"koo2016"}; default)
\itemize{
\item \strong{ICC < 0.50} - Poor reliability
\item \strong{0.5 <= ICC < 0.75} - Moderate reliability
\item \strong{0.75 <= ICC < 0.9} - Good reliability
\item **ICC >= 0.9 ** - Excellent reliability
}
}
}

\examples{
interpret_icc(0.6)
interpret_icc(c(0.4, 0.8))
}
\references{
\itemize{
\item Koo, T. K., \& Li, M. Y. (2016). A guideline of selecting and reporting intraclass correlation coefficients for reliability research. Journal of chiropractic medicine, 15(2), 155-163.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_oddsratio.R
\name{interpret_oddsratio}
\alias{interpret_oddsratio}
\title{Interpret Odds ratio}
\usage{
interpret_oddsratio(OR, rules = "chen2010", log = FALSE, ...)
}
\arguments{
\item{OR}{Value or vector of (log) odds ratio values.}

\item{rules}{Can be "\verb{chen2010"} (default), \code{"cohen1988"} (through
transformation to standardized difference, see \code{\link[=oddsratio_to_d]{oddsratio_to_d()}}) or custom set
of \code{\link[=rules]{rules()}}.}

\item{log}{Are the provided values log odds ratio.}

\item{...}{Currently not used.}
}
\description{
Interpret Odds ratio
}
\section{Rules}{


Rules apply to OR as ratios, so OR of 10 is as extreme as a OR of 0.1 (1/10).
\itemize{
\item Chen et al. (2010) (\code{"chen2010"}; default)
\itemize{
\item \strong{OR < 1.68} - Very small
\item \strong{1.68 <= OR < 3.47} - Small
\item \strong{3.47 <= OR < 6.71} - Medium
\item **OR >= 6.71 ** - Large
}
\item Cohen (1988) (\code{"cohen1988"}, based on the \code{\link[=oddsratio_to_d]{oddsratio_to_d()}} conversion, see \code{\link[=interpret_cohens_d]{interpret_cohens_d()}})
\itemize{
\item \strong{OR < 1.44} - Very small
\item \strong{1.44 <= OR < 2.48} - Small
\item \strong{2.48 <= OR < 4.27} - Medium
\item **OR >= 4.27 ** - Large
}
}
}

\examples{
interpret_oddsratio(1)
interpret_oddsratio(c(5, 2))

}
\references{
\itemize{
\item Cohen, J. (1988). Statistical power analysis for the behavioral sciences
(2nd Ed.). New York: Routledge.
\item Chen, H., Cohen, P., & Chen, S. (2010). How big is a big odds ratio?
Interpreting the magnitudes of odds ratios in epidemiological studies.
Communications in Statistics-Simulation and Computation, 39(4), 860-864.
\item Sánchez-Meca, J., Marín-Martínez, F., & Chacón-Moscoso, S. (2003).
Effect-size indices for dichotomized outcomes in meta-analysis. Psychological
methods, 8(4), 448.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_between_common_language.R
\name{d_to_cles}
\alias{d_to_cles}
\alias{convert_d_to_common_language}
\alias{d_to_common_language}
\alias{rb_to_cles}
\alias{rb_to_common_language}
\alias{convert_rb_to_common_language}
\title{Convert Standardized Mean Difference to Common Language Effect Sizes}
\usage{
d_to_cles(d)

rb_to_cles(rb)
}
\arguments{
\item{d, rb}{A numeric value of Cohen's d / rank-biserial correlation \emph{or}
the output from \code{\link[=cohens_d]{cohens_d()}} / \code{\link[=rank_biserial]{rank_biserial()}}.}
}
\value{
A list of \verb{Cohen's U3}, \code{Overlap}, \code{Pr(superiority)}, a
numeric vector of \code{Pr(superiority)}, or a data frame, depending
on the input.
}
\description{
Convert Standardized Mean Difference to Common Language Effect Sizes
}
\details{
This function use the following formulae for Cohen's \emph{d}:
\deqn{Pr(superiority) = \Phi(d/\sqrt{2})}{Pr(superiority) = pnorm(d / sqrt(2))}
\cr
\deqn{Cohen's U_3 = \Phi(d)}{U3 = pnorm(d)}
\cr
\deqn{Overlap = 2 \times \Phi(-|d|/2)}{Overlap = 2 * pnorm(-abs(d) / 2)}
\cr
And the following for the rank-biserial correlation:
\deqn{Pr(superiority) = (r_{rb} + 1)/2}{Pr(superiority) = (rb + 1)/2}
}
\note{
These calculations assume that the populations have equal variance and are
normally distributed.
}
\references{
\itemize{
\item Cohen, J. (1977). Statistical power analysis for the behavioral sciences.
New York: Routledge.
\item Reiser, B., & Faraggi, D. (1999). Confidence intervals for the overlapping
coefficient: the normal equal variance case. Journal of the Royal Statistical
Society, 48(3), 413-418.
\item Ruscio, J. (2008). A probability-based measure of effect size: robustness
to base rates and other factors. Psychological methods, 13(1), 19–30.
}
}
\seealso{
\code{\link[=cohens_d]{cohens_d()}}, \code{\link[=rank_biserial]{rank_biserial()}}

Other convert between effect sizes: 
\code{\link{d_to_r}()},
\code{\link{eta2_to_f2}()},
\code{\link{odds_to_probs}()},
\code{\link{oddsratio_to_riskratio}()}
}
\concept{convert between effect sizes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_effectsize_name.R
\name{is_effectsize_name}
\alias{is_effectsize_name}
\alias{get_effectsize_name}
\alias{get_effectsize_label}
\title{Checks if character is of a supported effect size}
\usage{
is_effectsize_name(x, ignore_case = TRUE)

get_effectsize_name(x, ignore_case = TRUE)

get_effectsize_label(x, ignore_case = TRUE)
}
\arguments{
\item{x}{A character, or a vector.}

\item{ignore_case}{Should case of input be ignored?}
}
\description{
For use by other functions and packages.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_cohens_g.R
\name{interpret_cohens_g}
\alias{interpret_cohens_g}
\title{Interpret Cohen's g}
\usage{
interpret_cohens_g(g, rules = "cohen1988", ...)
}
\arguments{
\item{g}{Value or vector of effect size values.}

\item{rules}{Can be \code{"cohen1988"} (default) or a custom set of \code{\link[=rules]{rules()}}.}

\item{...}{Not directly used.}
}
\description{
Interpret Cohen's g
}
\note{
"\emph{Since \strong{g} is so transparently clear a unit, it is expected that
workers in any given substantive area of the behavioral sciences will very
frequently be able to set relevant [effect size] values without the
proposed conventions, or set up conventions of their own which are suited
to their area of inquiry.}" - Cohen, 1988, page 147.
}
\section{Rules}{


Rules apply to equally to positive and negative \emph{g} (i.e., they are given as
absolute values).
\itemize{
\item Cohen (1988) (\code{"cohen1988"}; default)
\itemize{
\item \strong{d < 0.05} - Very small
\item \strong{0.05 <= d < 0.15} - Small
\item \strong{0.15 <= d < 0.25} - Medium
\item \strong{d >= 0.25} - Large
}
}
}

\examples{
interpret_cohens_g(.02)
interpret_cohens_g(c(.3, .15))

}
\references{
\itemize{
\item Cohen, J. (1988). Statistical power analysis for the behavioral sciences
(2nd Ed.). New York: Routledge.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexports.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{equivalence_test}
\alias{standardize}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{bayestestR}{\code{\link[bayestestR]{equivalence_test}}}

  \item{datawizard}{\code{\link[datawizard]{standardize}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_effectsize_name.R
\docType{data}
\name{es_info}
\alias{es_info}
\title{List of effect size names}
\format{
An object of class \code{data.frame} with 40 rows and 6 columns.
}
\usage{
es_info
}
\description{
Can always add more info here if need be...
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_extra.R
\name{effectsize_CIs}
\alias{effectsize_CIs}
\title{Confidence (Compatibility) Intervals}
\description{
More information regarding Confidence (Compatibiity) Intervals and how
they are computed in \emph{effectsize}.
}
\section{Confidence (Compatibility) Intervals (CIs)}{

Unless stated otherwise, confidence (compatibility) intervals (CIs) are
estimated using the noncentrality parameter method (also called the "pivot
method"). This method finds the noncentrality parameter ("\emph{ncp}") of a
noncentral \emph{t}, \emph{F}, or \eqn{\chi^2} distribution that places the observed
\emph{t}, \emph{F}, or \eqn{\chi^2} test statistic at the desired probability point of
the distribution. For example, if the observed \emph{t} statistic is 2.0, with 50
degrees of freedom, for which cumulative noncentral \emph{t} distribution is \emph{t} =
2.0 the .025 quantile (answer: the noncentral \emph{t} distribution with \emph{ncp} =
.04)? After estimating these confidence bounds on the \emph{ncp}, they are
converted into the effect size metric to obtain a confidence interval for the
effect size (Steiger, 2004).
\cr\cr
For additional details on estimation and troubleshooting, see \link{effectsize_CIs}.
}

\section{CIs and Significance Tests}{

"Confidence intervals on measures of effect size convey all the information
in a hypothesis test, and more." (Steiger, 2004). Confidence (compatibility)
intervals and p values are complementary summaries of parameter uncertainty
given the observed data. A dichotomous hypothesis test could be performed
with either a CI or a p value. The 100 (1 - \eqn{\alpha})\% confidence
interval contains all of the parameter values for which \emph{p} > \eqn{\alpha}
for the current data and model. For example, a 95\% confidence interval
contains all of the values for which p > .05.
\cr\cr
Note that a confidence interval including 0 \emph{does not} indicate that the null
(no effect) is true. Rather, it suggests that the observed data together with
the model and its assumptions combined do not provided clear evidence against
a parameter value of 0 (same as with any other value in the interval), with
the level of this evidence defined by the chosen \eqn{\alpha} level (Rafi &
Greenland, 2020; Schweder & Hjort, 2016; Xie & Singh, 2013). To infer no
effect, additional judgments about what parameter values are "close enough"
to 0 to be negligible are needed ("equivalence testing"; Bauer & Kiesser,
1996).
}

\section{One-Sided CIs}{

Typically, CIs are constructed as two-tailed intervals, with an equal
proportion of the cumulative probability distribution above and below the
interval. CIs can also be constructed as \emph{one-sided} intervals,
giving only a lower bound or upper bound. This is analogous to computing a
1-tailed \emph{p} value or conducting a 1-tailed hypothesis test.
\cr\cr
Significance tests conducted using CIs (whether a value is inside the interval)
and using \emph{p} values (whether p < alpha for that value) are only guaranteed
to agree when both are constructed using the same number of sides/tails.
\cr\cr
Most effect sizes are not bounded by zero (e.g., \emph{r}, \emph{d}, \emph{g}), and as such
are generally tested using 2-tailed tests and 2-sided CIs.
\cr\cr
Some effect sizes are strictly positive--they do have a minimum value, of 0.
For example, \eqn{R^2}, \eqn{\eta^2}, and other variance-accounted-for effect
sizes, as well as Cramer's \emph{V} and multiple \emph{R}, range from 0 to 1. These
typically involve \emph{F}- or \eqn{\chi^2}-statistics and are generally tested
using \emph{1-tailed} tests which test whether the estimated effect size is
\emph{larger} than the hypothesized null value (e.g., 0). In order for a CI to
yield the same significance decision it must then by a \emph{1-sided} CI,
estimating only a lower bound. This is the default CI computed by
\emph{effectsize} for these effect sizes, where \code{alternative = "greater"} is set.
\cr\cr
This lower bound interval indicates the smallest effect size that is not
significantly different from the observed effect size. That is, it is the
minimum effect size compatible with the observed data, background model
assumptions, and \eqn{\alpha} level. This type of interval does not indicate
a maximum effect size value; anything up to the maximum possible value of the
effect size (e.g., 1) is in the interval.
\cr\cr
One-sided CIs can also be used to test against a maximum effect size value
(e.g., is \eqn{R^2} significantly smaller than a perfect correlation of 1.0?)
can by setting \code{alternative = "less"}. This estimates a CI with only an
\emph{upper} bound; anything from the minimum possible value of the effect size
(e.g., 0) up to this upper bound is in the interval.
\cr\cr
We can also obtain a 2-sided interval by setting \code{alternative = "two-sided"}.
These intervals can be interpreted in the same way as other 2-sided
intervals, such as those for \emph{r}, \emph{d}, or \emph{g}.
\cr\cr
An alternative approach to aligning significance tests using CIs and 1-tailed
\emph{p} values that can often be found in the literature is to construct a
2-sided CI at a lower confidence level (e.g., 100(1-2\eqn{\alpha})\% = 100 -
2*5\% = 90\%. This estimates the lower bound and upper bound for the above
1-sided intervals simultaneously. These intervals are commonly reported when
conducting \strong{equivalence tests}. For example, a 90\% 2-sided interval gives
the bounds for an equivalence test with \eqn{\alpha} = .05. However, be aware
that this interval does not give 95\% coverage for the underlying effect size
parameter value. For that, construct a 95\% 2-sided CI.\if{html}{\out{<div class="sourceCode r">}}\preformatted{data("hardlyworking")
fit <- lm(salary ~ n_comps + age, data = hardlyworking)
eta_squared(fit) # default, ci = 0.95, alternative = "greater"
}\if{html}{\out{</div>}}\preformatted{## # Effect Size for ANOVA (Type I)
## 
## Parameter | Eta2 (partial) |       95\% CI
## -----------------------------------------
## n_comps   |           0.21 | [0.16, 1.00]
## age       |           0.10 | [0.06, 1.00]
## 
## - One-sided CIs: upper bound fixed at (1).
}\if{html}{\out{<div class="sourceCode r">}}\preformatted{eta_squared(fit, alternative = "less") # Test is eta is smaller than some value
}\if{html}{\out{</div>}}\preformatted{## # Effect Size for ANOVA (Type I)
## 
## Parameter | Eta2 (partial) |       95\% CI
## -----------------------------------------
## n_comps   |           0.21 | [0.00, 0.26]
## age       |           0.10 | [0.00, 0.14]
## 
## - One-sided CIs: lower bound fixed at (0).
}\if{html}{\out{<div class="sourceCode r">}}\preformatted{eta_squared(fit, alternative = "two.sided") # 2-sided bounds for alpha = .05
}\if{html}{\out{</div>}}\preformatted{## # Effect Size for ANOVA (Type I)
## 
## Parameter | Eta2 (partial) |       95\% CI
## -----------------------------------------
## n_comps   |           0.21 | [0.15, 0.27]
## age       |           0.10 | [0.06, 0.15]
}\if{html}{\out{<div class="sourceCode r">}}\preformatted{eta_squared(fit, ci = 0.9, alternative = "two.sided") # both 1-sided bounds for alpha = .05
}\if{html}{\out{</div>}}\preformatted{## # Effect Size for ANOVA (Type I)
## 
## Parameter | Eta2 (partial) |       90\% CI
## -----------------------------------------
## n_comps   |           0.21 | [0.16, 0.26]
## age       |           0.10 | [0.06, 0.14]
}
}

\section{CI Does Not Contain the Estimate}{

For very large sample sizes or effect sizes, the width of the CI can be
smaller than the tolerance of the optimizer, resulting in CIs of width 0.
This can also result in the estimated CIs excluding the point estimate.

For example:\if{html}{\out{<div class="sourceCode r">}}\preformatted{t_to_d(80, df_error = 4555555)
}\if{html}{\out{</div>}}\preformatted{## d    |       95\% CI
## -------------------
## 0.07 | [0.08, 0.08]
}

In these cases, consider an alternative optimizer, or an alternative method
for computing CIs, such as the bootstrap.
}

\references{
Bauer, P., & Kieser, M. (1996).
A unifying approach for confidence intervals and testing of equivalence and difference.
\emph{Biometrika, 83}(4), 934-–937.
\doi{10.1093/biomet/83.4.934}

Rafi, Z., & Greenland, S. (2020).
Semantic and cognitive tools to aid statistical science: Replace confidence and significance by compatibility and surprise.
\emph{BMC Medical Research Methodology, 20}(1), Article 244.
\doi{10.1186/s12874-020-01105-9}

Schweder, T., & Hjort, N. L. (2016).
\emph{Confidence, likelihood, probability: Statistical inference with confidence distributions.}
Cambridge University Press.
\doi{10.1017/CBO9781139046671}

Steiger, J. H. (2004).
Beyond the \emph{F} test: Effect size confidence intervals and tests of close fit in the analysis of variance and contrast analysis.
\emph{Psychological Methods, 9}(2), 164--182.
\doi{10.1037/1082-989x.9.2.164}

Xie, M., & Singh, K. (2013).
Confidence distribution, the frequentist distribution estimator of a parameter: A review.
\emph{International Statistical Review, 81}(1), 3–-39.
\doi{10.1111/insr.12000}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equivalence_test.R
\name{equivalence_test.effectsize_table}
\alias{equivalence_test.effectsize_table}
\title{Test for Practical Equivalence}
\usage{
\method{equivalence_test}{effectsize_table}(
  x,
  range = "default",
  rule = c("classic", "cet", "bayes"),
  ...
)
}
\arguments{
\item{x}{An effect size table, such as returned by \code{\link[=cohens_d]{cohens_d()}},
\code{\link[=eta_squared]{eta_squared()}}, \code{\link[=F_to_r]{F_to_r()}}, etc.}

\item{range}{The range of practical equivalence of an effect. For one-sides
CIs, a single value can be proved for the lower / upper bound to test
against (but see more details below). For two-sided CIs, a single value is
duplicated to \code{c(-range, range)}. If \code{"default"}, will be set to \verb{[-.1, .1]}.}

\item{rule}{How should acceptance and rejection be decided? See details.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame with the results of the equivalence test.
}
\description{
Perform a \strong{Test for Practical Equivalence} for indices of
effect size.
}
\details{
The CIs used in the equivalence test are the ones in the provided effect size
table. For results equivalent (ha!) to those that can be obtained using the
TOST approach (e.g., Lakens, 2017), appropriate CIs should be extracted using
the function used to make the effect size table (\code{cohens_d}, \code{eta_squared},
\code{F_to_r}, etc), with \code{alternative = "two.sided"}. See examples.
\subsection{The Different Rules}{
\itemize{
\item \code{"classic"} - \strong{the classic method}:
\itemize{
\item If the CI is completely within the ROPE - \emph{Accept H0}
\item Else, if the CI does not contain 0 - \emph{Reject H0}
\item Else - \emph{Undecided}
}
\item \code{"cet"} - \strong{conditional equivalence testing}:
\itemize{
\item If the CI does not contain 0 - \emph{Reject H0}
\item Else, If the CI is completely within the ROPE - \emph{Accept H0}
\item Else - \emph{Undecided}
}
\item \code{"bayes"} - \strong{The Bayesian approach}, as put forth by Kruschke:
\itemize{
\item If the CI does is completely outside the ROPE - \emph{Reject H0}
\item Else, If the CI is completely within the ROPE - \emph{Accept H0}
\item Else - \emph{Undecided}
}
}
}
}
\examples{
\donttest{

model <- aov(mpg ~ hp + am * factor(cyl), data = mtcars)
es <- eta_squared(model, ci = 0.9, alternative = "two.sided")
equivalence_test(es, range = 0.30) # TOST

RCT <- matrix(c(71, 101,
                50, 100), nrow = 2)
OR <- oddsratio(RCT, alternative = "greater")
equivalence_test(OR, range = 1)

ds <- t_to_d(
  t = c(0.45, -0.65, 7, -2.2, 2.25),
  df_error = c(675, 525, 2000, 900, 1875),
  ci = 0.9, alternative = "two.sided" # TOST
)
# Can also plot
if (require(see)) plot(equivalence_test(ds, range = 0.2))
if (require(see)) plot(equivalence_test(ds, range = 0.2, rule = "cet"))
if (require(see)) plot(equivalence_test(ds, range = 0.2, rule = "bayes"))
}

}
\references{
\itemize{
\item Campbell, H., & Gustafson, P. (2018). Conditional equivalence testing: An
alternative remedy for publication bias. PLOS ONE, 13(4), e0195145.
https://doi.org/10.1371/journal.pone.0195145
\item Kruschke, J. K. (2014). Doing Bayesian data analysis: A tutorial with R,
JAGS, and Stan. Academic Press
\item Kruschke, J. K. (2018). Rejecting or accepting parameter values in Bayesian
estimation. Advances in Methods and Practices in Psychological Science, 1(2),
270-280. doi: 10.1177/2515245918771304
\item Lakens, D. (2017). Equivalence Tests: A Practical Primer for t Tests,
Correlations, and Meta-Analyses. Social Psychological and Personality
Science, 8(4), 355–362. https://doi.org/10.1177/1948550617697177
}
}
\seealso{
For more details, see \code{\link[bayestestR:equivalence_test]{bayestestR::equivalence_test()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_stat_chisq.R
\name{chisq_to_phi}
\alias{chisq_to_phi}
\alias{chisq_to_cohens_w}
\alias{chisq_to_cramers_v}
\alias{chisq_to_pearsons_c}
\alias{phi_to_chisq}
\title{Conversion Chi-Squared to Phi or Cramer's V}
\usage{
chisq_to_phi(
  chisq,
  n,
  nrow,
  ncol,
  ci = 0.95,
  alternative = "greater",
  adjust = FALSE,
  ...
)

chisq_to_cohens_w(
  chisq,
  n,
  nrow,
  ncol,
  ci = 0.95,
  alternative = "greater",
  adjust = FALSE,
  ...
)

chisq_to_cramers_v(
  chisq,
  n,
  nrow,
  ncol,
  ci = 0.95,
  alternative = "greater",
  adjust = FALSE,
  ...
)

chisq_to_pearsons_c(
  chisq,
  n,
  nrow,
  ncol,
  ci = 0.95,
  alternative = "greater",
  ...
)

phi_to_chisq(phi, n, ...)
}
\arguments{
\item{chisq}{The Chi-squared statistic.}

\item{n}{Total sample size.}

\item{nrow, ncol}{The number of rows/columns in the contingency table (ignored
for Phi when \code{adjust=FALSE} and \code{CI=NULL}).}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"greater"} (default) or \code{"less"}
(one-sided CI), or \code{"two.sided"} (default, two-sided CI). Partial matching
is allowed (e.g., \code{"g"}, \code{"l"}, \code{"two"}...). See \emph{One-Sided CIs} in
\link{effectsize_CIs}.}

\item{adjust}{Should the effect size be bias-corrected? Defaults to \code{FALSE}.}

\item{...}{Arguments passed to or from other methods.}

\item{phi}{The Phi statistic.}
}
\value{
A data frame with the effect size(s), and confidence interval(s). See
\code{\link[=cramers_v]{cramers_v()}}.
}
\description{
Convert between Chi square (\eqn{\chi^2}), Cramer's V, phi (\eqn{\phi}),
Cohen's \emph{w} and Pearson's \emph{C} for contingency tables or goodness of fit.
}
\details{
These functions use the following formulae:
\cr
\deqn{\phi = \sqrt{\chi^2 / n}}{phi = sqrt(\chi^2 / n)}
\cr
\deqn{Cramer's V = \phi / \sqrt{min(nrow,ncol)-1}}{Cramer's V = \phi / sqrt(min(nrow,ncol)-1)}
\cr
\deqn{Pearson's C = \sqrt{\chi^2 / (\chi^2 + n)}}{Pearson's C = sqrt(\chi^2 / (\chi^2 + n))}
\cr\cr
For adjusted versions, see Bergsma, 2013.
}
\note{
Cohen's \emph{w} is equivalent to \emph{Phi}.
}
\section{Confidence (Compatibility) Intervals (CIs)}{

Unless stated otherwise, confidence (compatibility) intervals (CIs) are
estimated using the noncentrality parameter method (also called the "pivot
method"). This method finds the noncentrality parameter ("\emph{ncp}") of a
noncentral \emph{t}, \emph{F}, or \eqn{\chi^2} distribution that places the observed
\emph{t}, \emph{F}, or \eqn{\chi^2} test statistic at the desired probability point of
the distribution. For example, if the observed \emph{t} statistic is 2.0, with 50
degrees of freedom, for which cumulative noncentral \emph{t} distribution is \emph{t} =
2.0 the .025 quantile (answer: the noncentral \emph{t} distribution with \emph{ncp} =
.04)? After estimating these confidence bounds on the \emph{ncp}, they are
converted into the effect size metric to obtain a confidence interval for the
effect size (Steiger, 2004).
\cr\cr
For additional details on estimation and troubleshooting, see \link{effectsize_CIs}.
}

\section{CIs and Significance Tests}{

"Confidence intervals on measures of effect size convey all the information
in a hypothesis test, and more." (Steiger, 2004). Confidence (compatibility)
intervals and p values are complementary summaries of parameter uncertainty
given the observed data. A dichotomous hypothesis test could be performed
with either a CI or a p value. The 100 (1 - \eqn{\alpha})\% confidence
interval contains all of the parameter values for which \emph{p} > \eqn{\alpha}
for the current data and model. For example, a 95\% confidence interval
contains all of the values for which p > .05.
\cr\cr
Note that a confidence interval including 0 \emph{does not} indicate that the null
(no effect) is true. Rather, it suggests that the observed data together with
the model and its assumptions combined do not provided clear evidence against
a parameter value of 0 (same as with any other value in the interval), with
the level of this evidence defined by the chosen \eqn{\alpha} level (Rafi &
Greenland, 2020; Schweder & Hjort, 2016; Xie & Singh, 2013). To infer no
effect, additional judgments about what parameter values are "close enough"
to 0 to be negligible are needed ("equivalence testing"; Bauer & Kiesser,
1996).
}

\examples{
contingency_table <- as.table(rbind(c(762, 327, 468), c(484, 239, 477), c(484, 239, 477)))

# chisq.test(contingency_table)
#>
#>         Pearson's Chi-squared test
#>
#> data:  contingency_table
#> X-squared = 41.234, df = 4, p-value = 2.405e-08

chisq_to_phi(41.234,
  n = sum(contingency_table),
  nrow = nrow(contingency_table),
  ncol = ncol(contingency_table)
)
chisq_to_cramers_v(41.234,
  n = sum(contingency_table),
  nrow = nrow(contingency_table),
  ncol = ncol(contingency_table)
)
chisq_to_pearsons_c(41.234,
  n = sum(contingency_table),
  nrow = nrow(contingency_table),
  ncol = ncol(contingency_table)
)
}
\references{
\itemize{
\item Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and noncentral
distributions. Educational and Psychological Measurement, 61(4), 532-574.
\item Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
Journal of the Korean Statistical Society, 42(3), 323-328.
}
}
\seealso{
Other effect size from test statistic: 
\code{\link{F_to_eta2}()},
\code{\link{t_to_d}()}
}
\concept{effect size from test statistic}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_stat_to_d.R, R/convert_stat_to_r.R
\name{t_to_d}
\alias{t_to_d}
\alias{z_to_d}
\alias{F_to_d}
\alias{t_to_r}
\alias{z_to_r}
\alias{F_to_r}
\title{Convert test statistics (t, z, F) to effect sizes of differences (Cohen's d)
or association (\strong{partial} r)}
\usage{
t_to_d(
  t,
  df_error,
  paired = FALSE,
  ci = 0.95,
  alternative = "two.sided",
  pooled,
  ...
)

z_to_d(z, n, paired = FALSE, ci = 0.95, alternative = "two.sided", pooled, ...)

F_to_d(
  f,
  df,
  df_error,
  paired = FALSE,
  ci = 0.95,
  alternative = "two.sided",
  ...
)

t_to_r(t, df_error, ci = 0.95, alternative = "two.sided", ...)

z_to_r(z, n, ci = 0.95, alternative = "two.sided", ...)

F_to_r(f, df, df_error, ci = 0.95, alternative = "two.sided", ...)
}
\arguments{
\item{t, f, z}{The t, the F or the z statistics.}

\item{paired}{Should the estimate account for the t-value being testing the
difference between dependent means?}

\item{ci}{Confidence Interval (CI) level}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"two.sided"} (default, two-sided CI),
\code{"greater"} or \code{"less"} (one-sided CI). Partial matching is allowed (e.g.,
\code{"g"}, \code{"l"}, \code{"two"}...). See \emph{One-Sided CIs} in \link{effectsize_CIs}.}

\item{pooled}{Deprecated. Use \code{paired}.}

\item{...}{Arguments passed to or from other methods.}

\item{n}{The number of observations (the sample size).}

\item{df, df_error}{Degrees of freedom of numerator or of the error estimate
(i.e., the residuals).}
}
\value{
A data frame with the effect size(s)(\code{r} or \code{d}), and their CIs
(\code{CI_low} and \code{CI_high}).
}
\description{
These functions are convenience functions to convert t, z and F test
statistics to Cohen's d and \strong{partial} r. These are useful in cases where
the data required to compute these are not easily available or their
computation is not straightforward (e.g., in liner mixed models, contrasts,
etc.).
\cr
See \href{https://easystats.github.io/effectsize/articles/from_test_statistics.html}{Effect Size from Test Statistics vignette.}
}
\details{
These functions use the following formulae to approximate \emph{r} and \emph{d}:
\cr\cr
\deqn{r_{partial} = t / \sqrt{t^2 + df_{error}}}
\cr\cr
\deqn{r_{partial} = z / \sqrt{z^2 + N}}
\cr\cr
\deqn{d = 2 * t / \sqrt{df_{error}}}
\cr\cr
\deqn{d_z = t / \sqrt{df_{error}}}
\cr\cr
\deqn{d = 2 * z / \sqrt{N}}

The resulting \code{d} effect size is an \emph{approximation} to Cohen's \emph{d}, and
assumes two equal group sizes. When possible, it is advised to directly
estimate Cohen's \emph{d}, with \code{\link[=cohens_d]{cohens_d()}}, \code{emmeans::eff_size()}, or similar
functions.
}
\section{Confidence (Compatibility) Intervals (CIs)}{

Unless stated otherwise, confidence (compatibility) intervals (CIs) are
estimated using the noncentrality parameter method (also called the "pivot
method"). This method finds the noncentrality parameter ("\emph{ncp}") of a
noncentral \emph{t}, \emph{F}, or \eqn{\chi^2} distribution that places the observed
\emph{t}, \emph{F}, or \eqn{\chi^2} test statistic at the desired probability point of
the distribution. For example, if the observed \emph{t} statistic is 2.0, with 50
degrees of freedom, for which cumulative noncentral \emph{t} distribution is \emph{t} =
2.0 the .025 quantile (answer: the noncentral \emph{t} distribution with \emph{ncp} =
.04)? After estimating these confidence bounds on the \emph{ncp}, they are
converted into the effect size metric to obtain a confidence interval for the
effect size (Steiger, 2004).
\cr\cr
For additional details on estimation and troubleshooting, see \link{effectsize_CIs}.
}

\section{CIs and Significance Tests}{

"Confidence intervals on measures of effect size convey all the information
in a hypothesis test, and more." (Steiger, 2004). Confidence (compatibility)
intervals and p values are complementary summaries of parameter uncertainty
given the observed data. A dichotomous hypothesis test could be performed
with either a CI or a p value. The 100 (1 - \eqn{\alpha})\% confidence
interval contains all of the parameter values for which \emph{p} > \eqn{\alpha}
for the current data and model. For example, a 95\% confidence interval
contains all of the values for which p > .05.
\cr\cr
Note that a confidence interval including 0 \emph{does not} indicate that the null
(no effect) is true. Rather, it suggests that the observed data together with
the model and its assumptions combined do not provided clear evidence against
a parameter value of 0 (same as with any other value in the interval), with
the level of this evidence defined by the chosen \eqn{\alpha} level (Rafi &
Greenland, 2020; Schweder & Hjort, 2016; Xie & Singh, 2013). To infer no
effect, additional judgments about what parameter values are "close enough"
to 0 to be negligible are needed ("equivalence testing"; Bauer & Kiesser,
1996).
}

\examples{
## t Tests
res <- t.test(1:10, y = c(7:20), var.equal = TRUE)
t_to_d(t = res$statistic, res$parameter)
t_to_r(t = res$statistic, res$parameter)
t_to_r(t = res$statistic, res$parameter, alternative = "less")

res <- with(sleep, t.test(extra[group == 1], extra[group == 2], paired = TRUE))
t_to_d(t = res$statistic, res$parameter, paired = TRUE)
t_to_r(t = res$statistic, res$parameter)
t_to_r(t = res$statistic, res$parameter, alternative = "greater")

\donttest{
## Linear Regression
model <- lm(rating ~ complaints + critical, data = attitude)
(param_tab <- parameters::model_parameters(model))

(rs <- t_to_r(param_tab$t[2:3], param_tab$df_error[2:3]))

if (require(see)) plot(rs)

# How does this compare to actual partial correlations?
if (require("correlation")) {
  correlation::correlation(attitude[, c(1, 2, 6)], partial = TRUE)[1:2, c(2, 3, 7, 8)]
}
}

}
\references{
\itemize{
\item Friedman, H. (1982). Simplified determinations of statistical power,
magnitude of effect and research sample sizes. Educational and Psychological
Measurement, 42(2), 521-526. \doi{10.1177/001316448204200214}
\item Wolf, F. M. (1986). Meta-analysis: Quantitative methods for research
synthesis (Vol. 59). Sage.
\item Rosenthal, R. (1994) Parametric measures of effect size. In H. Cooper and
L.V. Hedges (Eds.). The handbook of research synthesis. New York: Russell
Sage Foundation.
\item Steiger, J. H. (2004). Beyond the F test: Effect size confidence intervals
and tests of close fit in the analysis of variance and contrast analysis.
Psychological Methods, 9, 164-182.
\item Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and noncentral
distributions. Educational and Psychological Measurement, 61(4), 532-574.
}
}
\seealso{
Other effect size from test statistic: 
\code{\link{F_to_eta2}()},
\code{\link{chisq_to_phi}()}
}
\concept{effect size from test statistic}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_extra.R, R/eta_squared.R
\name{effectsize_API}
\alias{effectsize_API}
\alias{.es_aov_simple}
\alias{.es_aov_strata}
\alias{.es_aov_table}
\title{\code{effectsize} API}
\usage{
.es_aov_simple(
  aov_table,
  type = c("eta", "omega", "epsilon"),
  partial = TRUE,
  generalized = FALSE,
  ci = 0.95,
  alternative = "greater",
  verbose = TRUE,
  include_intercept = FALSE
)

.es_aov_strata(
  aov_table,
  DV_names,
  type = c("eta", "omega", "epsilon"),
  partial = TRUE,
  generalized = FALSE,
  ci = 0.95,
  alternative = "greater",
  verbose = TRUE,
  include_intercept = FALSE
)

.es_aov_table(
  aov_table,
  type = c("eta", "omega", "epsilon"),
  partial = TRUE,
  generalized = FALSE,
  ci = 0.95,
  alternative = "greater",
  verbose = TRUE,
  include_intercept = FALSE
)
}
\arguments{
\item{aov_table}{Input data frame}

\item{type}{Which effect size to compute?}

\item{partial, generalized, ci, alternative, verbose}{See \code{\link[=eta_squared]{eta_squared()}}.}

\item{include_intercept}{Should the intercept (\code{(Intercept)}) be included?}

\item{DV_names}{A character vector with the names of all the predictors,
including the grouping variable (e.g., \code{"Subject"}).}
}
\description{
Read the \href{https://easystats.github.io/effectsize/articles/effectsize_API.html}{\emph{Support functions for model extensions}} vignette.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_cohens_d.R
\name{interpret_cohens_d}
\alias{interpret_cohens_d}
\alias{interpret_hedges_g}
\alias{interpret_glass_delta}
\title{Interpret standardized differences}
\usage{
interpret_cohens_d(d, rules = "cohen1988", ...)

interpret_hedges_g(g, rules = "cohen1988")

interpret_glass_delta(delta, rules = "cohen1988")
}
\arguments{
\item{d, g, delta}{Value or vector of effect size values.}

\item{rules}{Can be \code{"cohen1988"} (default), \code{"gignac2016"},
\code{"sawilowsky2009"}, \code{"lovakov2021"} or a custom set of \code{\link[=rules]{rules()}}.}

\item{...}{Not directly used.}
}
\description{
Interpretation of standardized differences using different sets of rules of
thumb.
}
\section{Rules}{


Rules apply to equally to positive and negative \emph{d} (i.e., they are given as
absolute values).
\itemize{
\item Cohen (1988) (\code{"cohen1988"}; default)
\itemize{
\item \strong{d < 0.2} - Very small
\item \strong{0.2 <= d < 0.5} - Small
\item \strong{0.5 <= d < 0.8} - Medium
\item \strong{d >= 0.8} - Large
}
\item Sawilowsky (2009) (\code{"sawilowsky2009"})
\itemize{
\item \strong{d < 0.1} - Tiny
\item \strong{0.1 <= d < 0.2} - Very small
\item \strong{0.2 <= d < 0.5} - Small
\item \strong{0.5 <= d < 0.8} - Medium
\item \strong{0.8 <= d < 1.2} - Large
\item \strong{1.2 <= d < 2} - Very large
\item \strong{d >= 2} - Huge
}
\item Lovakov & Agadullina (2021) (\code{"lovakov2021"})
\itemize{
\item \strong{d < 0.15} - Very small
\item \strong{0.15 <= d < 0.36} - Small
\item \strong{0.36 <= d < 0.65} - Medium
\item \strong{d >= 0.65} - Large
}
\item Gignac & Szodorai (2016) (\code{"gignac2016"}, based on the \code{\link[=d_to_r]{d_to_r()}} conversion, see \code{\link[=interpret_r]{interpret_r()}})
\itemize{
\item \strong{d < 0.2} - Very small
\item \strong{0.2 <= d < 0.41} - Small
\item \strong{0.41 <= d < 0.63} - Moderate
\item \strong{d >= 0.63} - Large
}
}
}

\examples{
interpret_cohens_d(.02)
interpret_cohens_d(c(.5, .02))
interpret_cohens_d(.3, rules = "lovakov2021")
}
\references{
\itemize{
\item Lovakov, A., & Agadullina, E. R. (2021). Empirically Derived Guidelines for
Effect Size Interpretation in Social Psychology. European Journal of Social
Psychology.
\item Gignac, G. E., & Szodorai, E. T. (2016). Effect size guidelines for
individual differences researchers. Personality and individual differences,
102, 74-78.
\item Cohen, J. (1988). Statistical power analysis for the behavioral sciences
(2nd Ed.). New York: Routledge.
\item Sawilowsky, S. S. (2009). New effect size rules of thumb.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_rope.R
\name{interpret_rope}
\alias{interpret_rope}
\title{Interpret Bayesian diagnostic indices}
\usage{
interpret_rope(rope, ci = 0.9, rules = "default")
}
\arguments{
\item{rope}{Value or vector of percentages in ROPE.}

\item{ci}{The Credible Interval (CI) probability, corresponding to the proportion of HDI, that was used. Can be \code{1} in the case of "full ROPE".}

\item{rules}{A character string (see details) or a custom set of \code{\link[=rules]{rules()}}.}
}
\description{
Interpretation of Bayesian indices of percentage in ROPE.
}
\section{Rules}{

\itemize{
\item Default
\itemize{
\item For CI < 1
\itemize{
\item \strong{Rope = 0} - Significant
\item \strong{0 < Rope < 1} - Undecided
\item \strong{Rope = 1} - Negligible
}
\item For CI = 1
\itemize{
\item \strong{Rope < 0.01} - Significant
\item \strong{0.01 < Rope < 0.025} - Probably significant
\item \strong{0.025 < Rope < 0.975} - Undecided
\item \strong{0.975 < Rope < 0.99} - Probably negligible
\item \strong{Rope > 0.99} - Negligible
}
}
}
}

\examples{
interpret_rope(0, ci = 0.9)
interpret_rope(c(0.005, 0.99), ci = 1)
}
\references{
\href{https://easystats.github.io/bayestestR/articles/guidelines.html}{BayestestR's reporting guidelines}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_bf.R
\name{interpret_bf}
\alias{interpret_bf}
\title{Interpret Bayes Factor (BF)}
\usage{
interpret_bf(
  bf,
  rules = "jeffreys1961",
  log = FALSE,
  include_value = FALSE,
  protect_ratio = TRUE,
  exact = TRUE
)
}
\arguments{
\item{bf}{Value or vector of Bayes factor (BF) values.}

\item{rules}{Can be \code{"jeffreys1961"} (default), \code{"raftery1995"} or custom set
of \code{\link[=rules]{rules()}} (for the \emph{absolute magnitude} of evidence).}

\item{log}{Is the \code{bf} value \code{log(bf)}?}

\item{include_value}{Include the value in the output.}

\item{protect_ratio}{Should values smaller than 1 be represented as ratios?}

\item{exact}{Should very large or very small values be reported with a
scientific format (e.g., 4.24e5), or as truncated values (as "> 1000" and
"< 1/1000").}
}
\description{
Interpret Bayes Factor (BF)
}
\details{
Argument names can be partially matched.
}
\section{Rules}{


Rules apply to BF as ratios, so BF of 10 is as extreme as a BF of 0.1 (1/10).
\itemize{
\item Jeffreys (1961) (\code{"jeffreys1961"}; default)
\itemize{
\item \strong{BF = 1} - No evidence
\item \strong{1 < BF <= 3} - Anecdotal
\item \strong{3 < BF <= 10} - Moderate
\item \strong{10 < BF <= 30} - Strong
\item \strong{30 < BF <= 100} - Very strong
\item \strong{BF > 100} - Extreme.
}
\item Raftery (1995) (\code{"raftery1995"})
\itemize{
\item \strong{BF = 1} - No evidence
\item \strong{1 < BF <= 3} - Weak
\item \strong{3 < BF <= 20} - Positive
\item \strong{20 < BF <= 150} - Strong
\item \strong{BF > 150} - Very strong
}
}
}

\examples{
interpret_bf(1)
interpret_bf(c(5, 2))
}
\references{
\itemize{
\item Jeffreys, H. (1961), Theory of Probability, 3rd ed., Oxford University
Press, Oxford.
\item Raftery, A. E. (1995). Bayesian model selection in social research.
Sociological methodology, 25, 111-164.
\item Jarosz, A. F., & Wiley, J. (2014). What are the odds? A practical guide to
computing and reporting Bayes factors. The Journal of Problem Solving, 7(1),
}
\enumerate{
\item 
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_between_d_to_r.R
\name{d_to_r}
\alias{d_to_r}
\alias{convert_d_to_r}
\alias{r_to_d}
\alias{convert_r_to_d}
\alias{oddsratio_to_d}
\alias{convert_oddsratio_to_d}
\alias{logoddsratio_to_d}
\alias{convert_logoddsratio_to_d}
\alias{d_to_oddsratio}
\alias{convert_d_to_oddsratio}
\alias{oddsratio_to_r}
\alias{convert_oddsratio_to_r}
\alias{logoddsratio_to_r}
\alias{convert_logoddsratio_to_r}
\alias{r_to_oddsratio}
\alias{convert_r_to_oddsratio}
\title{Convert between \emph{d}, \emph{r} and \emph{Odds ratio}}
\usage{
d_to_r(d, ...)

r_to_d(r, ...)

oddsratio_to_d(OR, log = FALSE, ...)

logoddsratio_to_d(OR, log = TRUE, ...)

d_to_oddsratio(d, log = FALSE, ...)

oddsratio_to_r(OR, log = FALSE, ...)

logoddsratio_to_r(OR, log = TRUE, ...)

r_to_oddsratio(r, log = FALSE, ...)
}
\arguments{
\item{d}{Standardized difference value (Cohen's d).}

\item{...}{Arguments passed to or from other methods.}

\item{r}{Correlation coefficient r.}

\item{OR}{\emph{Odds ratio} values in vector or data frame.}

\item{log}{Take in or output the log of the ratio (such as in logistic models).}
}
\value{
Converted index.
}
\description{
Enables a conversion between different indices of effect size, such as
standardized difference (Cohen's d), correlation r or (log) odds ratios.
}
\details{
Conversions between \emph{d} and \emph{OR} or \emph{r} is done through these formulae.
\itemize{
\item \eqn{d = \frac{2 * r}{\sqrt{1 - r^2}}}{d = 2 * r / sqrt(1 - r^2)}
\item \eqn{r = \frac{d}{\sqrt{d^2 + 4}}}{r = d / sqrt(d^2 + 4)}
\item \eqn{d = \frac{\log(OR)\times\sqrt{3}}{\pi}}{d = log(OR) * sqrt(3) / pi}
\item \eqn{log(OR) = d * \frac{\pi}{\sqrt(3)}}{log(OR) = d * pi / sqrt(3)}
}

The conversion from \emph{d} to \emph{r} assumes equally sized groups. The resulting
\emph{r} is also called the binomial effect size display (BESD; Rosenthal et al.,
1982).
}
\examples{
r_to_d(0.5)
d_to_oddsratio(1.154701)
oddsratio_to_r(8.120534)

d_to_r(1)
r_to_oddsratio(0.4472136, log = TRUE)
oddsratio_to_d(1.813799, log = TRUE)

}
\references{
\itemize{
\item Sánchez-Meca, J., Marín-Martínez, F., & Chacón-Moscoso, S. (2003).
Effect-size indices for dichotomized outcomes in meta-analysis. Psychological
methods, 8(4), 448.
\item Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
(2009). Converting among effect sizes. Introduction to meta-analysis, 45-49.
\item Rosenthal, R., & Rubin, D. B. (1982). A simple, general purpose display of
magnitude of experimental effect. Journal of educational psychology, 74(2), 166.
}
}
\seealso{
Other convert between effect sizes: 
\code{\link{d_to_cles}()},
\code{\link{eta2_to_f2}()},
\code{\link{odds_to_probs}()},
\code{\link{oddsratio_to_riskratio}()}
}
\concept{convert between effect sizes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_vif.R
\name{interpret_vif}
\alias{interpret_vif}
\title{Interpret the Variance Inflation Factor (VIF)}
\usage{
interpret_vif(vif, rules = "default")
}
\arguments{
\item{vif}{Value or vector of VIFs.}

\item{rules}{Can be \code{"default"} or a custom set of \code{\link[=rules]{rules()}}.}
}
\description{
Interpret VIF index of multicollinearity.
}
\section{Rules}{

\itemize{
\item Default
\itemize{
\item \strong{VIF < 5} - Low
\item \strong{5 <= VIF < 10} - Moderate
\item \strong{VIF >= 10} - High
}
}
}

\examples{

interpret_vif(c(1.4, 30.4))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_between_odds_to_probs.R
\name{odds_to_probs}
\alias{odds_to_probs}
\alias{convert_odds_to_probs}
\alias{odds_to_probs.data.frame}
\alias{probs_to_odds}
\alias{convert_probs_to_odds}
\alias{probs_to_odds.data.frame}
\title{Convert between Odds and Probabilities}
\usage{
odds_to_probs(odds, log = FALSE, ...)

\method{odds_to_probs}{data.frame}(odds, log = FALSE, select = NULL, exclude = NULL, ...)

probs_to_odds(probs, log = FALSE, ...)

\method{probs_to_odds}{data.frame}(probs, log = FALSE, select = NULL, exclude = NULL, ...)
}
\arguments{
\item{odds}{The \emph{Odds} (or \code{log(odds)} when \code{log = TRUE}) to convert.}

\item{log}{Take in or output log odds (such as in logistic models).}

\item{...}{Arguments passed to or from other methods.}

\item{select}{When a data frame is passed, character or list of of column
names to be transformed.}

\item{exclude}{When a data frame is passed, character or list of column names
to be excluded from transformation.}

\item{probs}{Probability values to convert.}
}
\value{
Converted index.
}
\description{
Convert between Odds and Probabilities
}
\examples{
odds_to_probs(3)
odds_to_probs(1.09, log = TRUE)

probs_to_odds(0.95)
probs_to_odds(0.95, log = TRUE)
}
\seealso{
\code{\link[stats:Logistic]{stats::plogis()}}

Other convert between effect sizes: 
\code{\link{d_to_cles}()},
\code{\link{d_to_r}()},
\code{\link{eta2_to_f2}()},
\code{\link{oddsratio_to_riskratio}()}
}
\concept{convert between effect sizes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.models.R
\name{standardize.default}
\alias{standardize.default}
\alias{standardize_models}
\alias{standardize.models}
\title{Re-fit a model with standardized data}
\usage{
\method{standardize}{default}(
  x,
  robust = FALSE,
  two_sd = FALSE,
  weights = TRUE,
  verbose = TRUE,
  include_response = TRUE,
  ...
)
}
\arguments{
\item{x}{A statistical model.}

\item{robust}{Logical, if \code{TRUE}, centering is done by subtracting the
median from the variables and dividing it by the median absolute deviation
(MAD). If \code{FALSE}, variables are standardized by subtracting the
mean and dividing it by the standard deviation (SD).}

\item{two_sd}{If \code{TRUE}, the variables are scaled by two times the deviation
(SD or MAD depending on \code{robust}). This method can be useful to obtain
model coefficients of continuous parameters comparable to coefficients
related to binary predictors, when applied to \strong{the predictors} (not the
outcome) (Gelman, 2008).}

\item{weights}{If \code{TRUE} (default), a weighted-standardization is carried out.}

\item{verbose}{Toggle warnings and messages on or off.}

\item{include_response}{If \code{TRUE} (default), the response value will also be
standardized. If \code{FALSE}, only the predictors will be standardized.
\itemize{
\item Note that for GLMs and models with non-linear link functions, the
response value will not be standardized, to make re-fitting the model work.
\item If the model contains an \code{\link[stats:offset]{stats::offset()}}, the offset variable(s) will
be standardized only if the response is standardized. If \code{two_sd = TRUE},
offsets are standardized by one-sd (similar to the response).
\item (For \code{mediate} models, the \code{include_response} refers to the outcome in
the y model; m model's response will always be standardized when possible).
}}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A statistical model fitted on standardized data
}
\description{
Performs a standardization of data (z-scoring) using
\code{\link[datawizard:standardize]{datawizard::standardize()}} and then re-fits the model to the standardized
data.
\cr\cr
Standardization is done by completely refitting the model on the standardized
data. Hence, this approach is equal to standardizing the variables \emph{before}
fitting the model and will return a new model object. This method is
particularly recommended for complex models that include interactions or
transformations (e.g., polynomial or spline terms). The \code{robust} (default to
\code{FALSE}) argument enables a robust standardization of data, based on the
\code{median} and the \code{MAD} instead of the \code{mean} and the \code{SD}.
}
\section{Generalized Linear Models}{
Standardization for generalized linear models (GLM, GLMM, etc) is done only
with respect to the predictors (while the outcome remains as-is,
unstandardized) - maintaining the interpretability of the coefficients (e.g.,
in a binomial model: the exponent of the standardized parameter is the OR of
a change of 1 SD in the predictor, etc.)
}

\section{Dealing with Factors}{
\code{standardize(model)} or \code{standardize_parameters(model, method = "refit")} do
\emph{not} standardized categorical predictors (i.e. factors) / their
dummy-variables, which may be a different behaviour compared to other R
packages (such as \pkg{lm.beta}) or other software packages (like SPSS). To
mimic such behaviours, either use \code{standardize_parameters(model, method = "basic")} to obtain post-hoc standardized parameters, or standardize the data
with \code{datawizard::standardize(data, force = TRUE)} \emph{before} fitting the
model.
}

\section{Transformed Variables}{
When the model's formula contains transformations (e.g. \code{y ~ exp(X)}) the
transformation effectively takes place after standardization (e.g.,
\code{exp(scale(X))}). Since some transformations are undefined for none positive
values, such as \code{log()} and \code{sqrt()}, the releven variables are shifted (post
standardization) by \code{Z - min(Z) + 1} or \code{Z - min(Z)} (respectively).
}

\examples{
model <- lm(Infant.Mortality ~ Education * Fertility, data = swiss)
coef(standardize(model))

}
\seealso{
Other standardize: 
\code{\link{standardize_info}()},
\code{\link{standardize_parameters}()}
}
\concept{standardize}
