
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
