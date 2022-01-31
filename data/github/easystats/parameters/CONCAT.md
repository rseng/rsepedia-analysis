
# parameters <img src="man/figures/logo.png" align="right" height="139" />

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02445/status.svg)](https://doi.org/10.21105/joss.02445)
[![downloads](http://cranlogs.r-pkg.org/badges/parameters)](https://cran.r-project.org/package=parameters)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/parameters)](https://cranlogs.r-pkg.org/)
[![status](https://tinyverse.netlify.com/badge/parameters)](https://CRAN.R-project.org/package=parameters)

------------------------------------------------------------------------

:warning: For Bayesian models, we changed the default the CI width!
Please make an [informed
decision](https://easystats.github.io/bayestestR/articles/credible_interval.html)
and set it explicitly (`ci = 0.89`, `ci = 0.95`, or anything else that
you decide) :warning:

------------------------------------------------------------------------

***Describe and understand your model’s parameters!***

**parameters**’ primary goal is to provide utilities for processing the
parameters of various statistical models (see
[here](https://easystats.github.io/insight/) for a list of supported
models). Beyond computing *p-values*, *CIs*, *Bayesian indices* and
other measures for a wide variety of models, this package implements
features like *bootstrapping* of parameters and models, *feature
reduction* (feature extraction and variable selection), or tools for
data reduction like functions to perform cluster, factor or principal
component analysis.

Another important goal of the **parameters** package is to facilitate
and streamline the process of reporting results of statistical models,
which includes the easy and intuitive calculation of standardized
estimates or robust standard errors and p-values. **parameters**
therefor offers a simple and unified syntax to process a large variety
of (model) objects from many different packages.

## Installation

[![CRAN](http://www.r-pkg.org/badges/version/parameters)](https://cran.r-project.org/package=parameters)
[![parameters status
badge](https://easystats.r-universe.dev/badges/parameters)](https://easystats.r-universe.dev)
[![R-check](https://github.com/easystats/parameters/workflows/R-check/badge.svg?branch=main)](https://github.com/easystats/parameters/actions)

Run the following to install the stable release of **parameters** from
CRAN:

``` r
install.packages("parameters")
```

Or this one to install the latest development version from R-universe…

``` r
install.packages("parameters", repos = "https://easystats.r-universe.dev")
```

…or from GitHub:

``` r
install.packages("remotes")
remotes::install_github("easystats/parameters")
```

## Documentation

[![Documentation](https://img.shields.io/badge/documentation-parameters-orange.svg?colorB=E91E63)](https://easystats.github.io/parameters/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-parameters-orange.svg?colorB=2196F3)](https://easystats.github.io/parameters/reference/index.html)

Click on the buttons above to access the package
[documentation](https://easystats.github.io/parameters/) and the
[easystats blog](https://easystats.github.io/blog/posts/), and check-out
these vignettes:

-   [Summary of Model
    Parameters](https://easystats.github.io/parameters/articles/model_parameters.html)
-   [Standardized Model
    Parameters](https://easystats.github.io/parameters/articles/model_parameters_standardized.html)
-   [Robust Estimation of Standard Errors, Confidence Intervals and
    p-values](https://easystats.github.io/parameters/articles/model_parameters_robust.html)
-   [Model Parameters and Missing
    Data](https://easystats.github.io/parameters/articles/model_parameters_mice.html)
-   [Feature reduction (PCA, cMDS,
    ICA…)](https://easystats.github.io/parameters/articles/parameters_reduction.html)
-   [Structural models (EFA, CFA,
    SEM…)](https://easystats.github.io/parameters/articles/efa_cfa.html)
-   [Parameters
    selection](https://easystats.github.io/parameters/articles/parameters_selection.html)
-   [A Practical Guide for Panel Data
    Analysis](https://easystats.github.io/datawizard/articles/demean.html)

## Contributing and Support

In case you want to file an issue or contribute in another way to the
package, please follow [this
guide](https://github.com/easystats/parameters/blob/master/.github/CONTRIBUTING.md).
For questions about the functionality, you may either contact us via
email or also file an issue.

# Features

## Model’s parameters description

<img src="man/figures/figure1.png" style="display: block; margin: auto;" />

The
[`model_parameters()`](https://easystats.github.io/parameters/articles/model_parameters.html)
function (that can be accessed via the `parameters()` shortcut) allows
you to extract the parameters and their characteristics from various
models in a consistent way. It can be considered as a lightweight
alternative to [`broom::tidy()`](https://github.com/tidymodels/broom),
with some notable differences:

-   The column names of the returned data frame are *specific* to their
    content. For instance, the column containing the statistic is named
    following the statistic name, i.e., *t*, *z*, etc., instead of a
    generic name such as *statistic* (however, you can get standardized
    (generic) column names using
    [`standardize_names()`](https://easystats.github.io/insight/reference/standardize_names.html)).
-   It is able to compute or extract indices not available by default,
    such as *p-values*, *CIs*, etc.
-   It includes *feature engineering* capabilities, including parameters
    [bootstrapping](https://easystats.github.io/parameters/reference/bootstrap_parameters.html).

### Classical Regression Models

``` r
model <- lm(Sepal.Width ~ Petal.Length * Species + Petal.Width, data = iris)

# regular model parameters
model_parameters(model)
#> Parameter                           | Coefficient |   SE |         95% CI | t(143) |      p
#> -------------------------------------------------------------------------------------------
#> (Intercept)                         |        2.89 | 0.36 | [ 2.18,  3.60] |   8.01 | < .001
#> Petal Length                        |        0.26 | 0.25 | [-0.22,  0.75] |   1.07 | 0.287 
#> Species [versicolor]                |       -1.66 | 0.53 | [-2.71, -0.62] |  -3.14 | 0.002 
#> Species [virginica]                 |       -1.92 | 0.59 | [-3.08, -0.76] |  -3.28 | 0.001 
#> Petal Width                         |        0.62 | 0.14 | [ 0.34,  0.89] |   4.41 | < .001
#> Petal Length * Species [versicolor] |       -0.09 | 0.26 | [-0.61,  0.42] |  -0.36 | 0.721 
#> Petal Length * Species [virginica]  |       -0.13 | 0.26 | [-0.64,  0.38] |  -0.50 | 0.618

# standardized parameters
model_parameters(model, standardize = "refit")
#> Parameter                           | Coefficient |   SE |         95% CI | t(143) |      p
#> -------------------------------------------------------------------------------------------
#> (Intercept)                         |        3.59 | 1.30 | [ 1.01,  6.17] |   2.75 | 0.007 
#> Petal Length                        |        1.07 | 1.00 | [-0.91,  3.04] |   1.07 | 0.287 
#> Species [versicolor]                |       -4.62 | 1.31 | [-7.21, -2.03] |  -3.53 | < .001
#> Species [virginica]                 |       -5.51 | 1.38 | [-8.23, -2.79] |  -4.00 | < .001
#> Petal Width                         |        1.08 | 0.24 | [ 0.59,  1.56] |   4.41 | < .001
#> Petal Length * Species [versicolor] |       -0.38 | 1.06 | [-2.48,  1.72] |  -0.36 | 0.721 
#> Petal Length * Species [virginica]  |       -0.52 | 1.04 | [-2.58,  1.54] |  -0.50 | 0.618
```

### Mixed Models

``` r
library(lme4)

model <- lmer(Sepal.Width ~ Petal.Length + (1|Species), data = iris)

# model parameters with CI, df and p-values based on Wald approximation
model_parameters(model, effects = "all")
#> # Fixed Effects
#> 
#> Parameter    | Coefficient |   SE |       95% CI | t(146) |      p
#> ------------------------------------------------------------------
#> (Intercept)  |        2.00 | 0.56 | [0.89, 3.11] |   3.56 | < .001
#> Petal Length |        0.28 | 0.06 | [0.16, 0.40] |   4.75 | < .001
#> 
#> # Random Effects
#> 
#> Parameter               | Coefficient
#> -------------------------------------
#> SD (Intercept: Species) |        0.89
#> SD (Residual)           |        0.32

# model parameters with CI, df and p-values based on Kenward-Roger approximation
model_parameters(model, df_method = "kenward")
#> # Fixed Effects
#> 
#> Parameter    | Coefficient |   SE |       95% CI |    t |     df |      p
#> -------------------------------------------------------------------------
#> (Intercept)  |        2.00 | 0.57 | [0.07, 3.93] | 3.53 |   2.67 | 0.046 
#> Petal Length |        0.28 | 0.06 | [0.16, 0.40] | 4.58 | 140.98 | < .001
#> 
#> # Random Effects
#> 
#> Parameter               | Coefficient
#> -------------------------------------
#> SD (Intercept: Species) |        0.89
#> SD (Residual)           |        0.32
```

### Structural Models

Besides many types of regression models and packages, it also works for
other types of models, such as [**structural
models**](https://easystats.github.io/parameters/articles/efa_cfa.html)
(EFA, CFA, SEM…).

``` r
library(psych)

model <- psych::fa(attitude, nfactors = 3)
model_parameters(model)
#> # Rotated loadings from Factor Analysis (oblimin-rotation)
#> 
#> Variable   |  MR1  |  MR2  |  MR3  | Complexity | Uniqueness
#> ------------------------------------------------------------
#> rating     | 0.90  | -0.07 | -0.05 |    1.02    |    0.23   
#> complaints | 0.97  | -0.06 | 0.04  |    1.01    |    0.10   
#> privileges | 0.44  | 0.25  | -0.05 |    1.64    |    0.65   
#> learning   | 0.47  | 0.54  | -0.28 |    2.51    |    0.24   
#> raises     | 0.55  | 0.43  | 0.25  |    2.35    |    0.23   
#> critical   | 0.16  | 0.17  | 0.48  |    1.46    |    0.67   
#> advance    | -0.11 | 0.91  | 0.07  |    1.04    |    0.22   
#> 
#> The 3 latent factors (oblimin rotation) accounted for 66.60% of the total variance of the original data (MR1 = 38.19%, MR2 = 22.69%, MR3 = 5.72%).
```

## Variable and parameters selection

<img src="man/figures/figure2.png" style="display: block; margin: auto;" />

[`select_parameters()`](https://easystats.github.io/parameters/articles/parameters_selection.html)
can help you quickly select and retain the most relevant predictors
using methods tailored for the model type.

``` r
library(poorman)

lm(disp ~ ., data = mtcars) %>% 
  select_parameters() %>% 
  model_parameters()
#> Parameter   | Coefficient |     SE |            95% CI | t(26) |      p
#> -----------------------------------------------------------------------
#> (Intercept) |      141.70 | 125.67 | [-116.62, 400.02] |  1.13 | 0.270 
#> cyl         |       13.14 |   7.90 | [  -3.10,  29.38] |  1.66 | 0.108 
#> hp          |        0.63 |   0.20 | [   0.22,   1.03] |  3.18 | 0.004 
#> wt          |       80.45 |  12.22 | [  55.33, 105.57] |  6.58 | < .001
#> qsec        |      -14.68 |   6.14 | [ -27.31,  -2.05] | -2.39 | 0.024 
#> carb        |      -28.75 |   5.60 | [ -40.28, -17.23] | -5.13 | < .001
```

## Miscellaneous

This packages also contains a lot of [other useful
functions](https://easystats.github.io/parameters/reference/index.html):

### Describe a Distribution

``` r
data(iris)
describe_distribution(iris)
#> Variable     | Mean |   SD |  IQR |        Range | Skewness | Kurtosis |   n | n_Missing
#> ----------------------------------------------------------------------------------------
#> Sepal.Length | 5.84 | 0.83 | 1.30 | [4.30, 7.90] |     0.31 |    -0.55 | 150 |         0
#> Sepal.Width  | 3.06 | 0.44 | 0.52 | [2.00, 4.40] |     0.32 |     0.23 | 150 |         0
#> Petal.Length | 3.76 | 1.77 | 3.52 | [1.00, 6.90] |    -0.27 |    -1.40 | 150 |         0
#> Petal.Width  | 1.20 | 0.76 | 1.50 | [0.10, 2.50] |    -0.10 |    -1.34 | 150 |         0
```

### Citation

In order to cite this package, please use the following command:

``` r
citation("parameters")

Lüdecke D, Ben-Shachar M, Patil I, Makowski D (2020). "Extracting, Computing and
Exploring the Parameters of Statistical Models using R." _Journal of Open Source
Software_, *5*(53), 2445. doi: 10.21105/joss.02445 (URL:
https://doi.org/10.21105/joss.02445).

A BibTeX entry for LaTeX users is

  @Article{,
    title = {Extracting, Computing and Exploring the Parameters of Statistical Models using {R}.},
    volume = {5},
    doi = {10.21105/joss.02445},
    number = {53},
    journal = {Journal of Open Source Software},
    author = {Daniel Lüdecke and Mattan S. Ben-Shachar and Indrajeet Patil and Dominique Makowski},
    year = {2020},
    pages = {2445},
  }
```
# parameters 0.15.1

## General

* Improved speed performance for `model_parameters()`, in particular for glm's
  and mixed models where random effect variances were calculated.

* Added more options for printing `model_parameters()`. See also revised vignette:
  https://easystats.github.io/parameters/articles/model_parameters_print.html

## Changes to functions

### `model_parameters()`

* `model_parameters()` for mixed models gains an `include_sigma` argument. If
  `TRUE`, adds the residual variance, computed from the random effects variances,
  as an attribute to the returned data frame. Including sigma was the default 
  behaviour, but now defaults to `FALSE` and is only included when 
  `include_sigma = TRUE`, because the calculation was very time consuming.

* `model_parameters()` for `merMod` models now also computes CIs for the random
  SD parameters when `ci_method="boot"` (previously, this was only possible when
  `ci_method` was `"profile"`).

* `model_parameters()` for `glmmTMB` models now computes CIs for the random SD
  parameters. Note that these are based on a Wald-z-distribution.

* Similar to `model_parameters.htest()`, the `model_parameters.BFBayesFactor()`
  method gains `cohens_d` and `cramers_v` arguments to control if you need to
  add frequentist effect size estimates to the returned summary data frame.
  Previously, this was done by default.

* Column name for coefficients from *emmeans* objects are now more specific.

* `model_prameters()` for `MixMod` objects (package *GLMMadaptive*) gains a
  `robust` argument, to compute robust standard errors.

## Bug fixes 

* Fixed bug with `ci()` for class `merMod` when `method="boot"`.

* Fixed issue with correct association of components for ordinal models of 
  classes `clm` and `clm2`.

* Fixed issues in `random_parameters()` and  `model_parameters()` for mixed 
  models without random intercept.

* Confidence intervals for random parameters in `model_parameters()` failed for
  (some?) `glmer` models.

* Fix issue with default `ci_type` in `compare_parameters()` for Bayesian models.

# parameters 0.15.0

## Breaking changes

* Following functions were moved to the new *datawizard* package and are now
  re-exported from *parameters* package:

  - `center()`

  - `convert_data_to_numeric()`

  - `data_partition()`

  - `demean()` (and its aliases `degroup()` and `detrend()`)

  - `kurtosis()`

  - `rescale_weights()`

  - `skewness()`

  - `smoothness()`

Note that these functions will be removed in the next release of *parameters*
package and they are currently being re-exported only as a convenience for the
package developers. This release should provide them with time to make the
necessary changes before this breaking change is implemented.

* Following functions were moved to the *performance* package:

  - `check_heterogeneity()`

  - `check_multimodal()`

## General

* The handling to approximate the degrees of freedom in `model_parameters()`,
  `ci()` and `p_value()` was revised and should now be more consistent. Some
  bugs related to the previous computation of confidence intervals and p-values
  have been fixed. Now it is possible to change the method to approximate
  degrees of freedom for CIs and p-values using the `ci_method`, resp. `method`
  argument. This change has been documented in detail in `?model_parameters`,
  and online here:
  https://easystats.github.io/parameters/reference/model_parameters.html

* Minor changes to `print()` for *glmmTMB* with dispersion parameter.

* Added vignette on printing options for model parameters.

## Changes to functions

### `model_parameters()`

* The `df_method` argument in `model_parameters()` is deprecated. Please use
  `ci_method` now.

* `model_parameters()` with `standardize = "refit"` now returns random effects
  from the standardized model.

* `model_parameters()` and `ci()` for `lmerMod` models gain a `"residuals"`
  option for the `ci_method` (resp. `method`) argument, to explicitly calculate
  confidence intervals based on the residual degrees of freedom, when present.

* `model_parameters()` supports following new objects: `trimcibt`, `wmcpAKP`,
  `dep.effect` (in *WRS2* package), `systemfit`

* `model_parameters()` gains a new argument `table_wide` for ANOVA tables. This
  can be helpful for users who may wish to report ANOVA table in wide format
  (i.e., with numerator and denominator degrees of freedom on the same row).

* `model_parameters()` gains two new arguments, `keep` and `drop`. `keep` is the
  new names for the former `parameters` argument and can be used to filter
  parameters. While `keep` selects those parameters whose names match the
  regular expression pattern defined in `keep`, `drop` is the counterpart and
  excludes matching parameter names.

* When `model_parameters()` is called with `verbose = TRUE`, and `ci_method` is
  not the default value, the printed output includes a message indicating which
  approximation-method for degrees of freedom was used.

* `model_parameters()` for mixed models with `ci_method = "profile` computes
  (profiled) confidence intervals for both fixed and random effects. Thus,
  `ci_method = "profile` allows to add confidence intervals to the random effect
  variances.

* `model_parameters()` should longer fail for supported model classes when
  robust standard errors are not available.

### Other functions

* `n_factors()` the methods based on fit indices have been fixed and can be
  included separately (`package = "fit"`). Also added a `n_max` argument to crop
  the output.

* `compare_parameters()` now also accepts a list of model objects.

* `describe_distribution()` gets `verbose` argument to toggle warnings and
  messages.

* `format_parameters()` removes dots and underscores from parameter names, to
  make these more "human readable".

* The experimental calculation of p-values in `equivalence_test()` was replaced
  by a proper calculation p-values. The argument `p_value` was removed and
  p-values are now always included.

* Minor improvements to `print()`, `print_html()` and `print_md()`.

## Bug fixes

* The random effects returned by `model_parameters()` mistakenly displayed the
  residuals standard deviation as square-root of the residual SD.

* Fixed issue with `model_parameters()` for *brmsfit* objects that model
  standard errors (i.e. for meta-analysis).

* Fixed issue in `model_parameters` for `lmerMod` models that, by default,
  returned residual degrees of freedom in the statistic column, but confidence
  intervals were based on `Inf` degrees of freedom instead.

* Fixed issue in `ci_satterthwaite()`, which used `Inf` degrees of freedom
  instead of the Satterthwaite approximation.

* Fixed issue in `model_parameters.mlm()` when model contained interaction
  terms.

* Fixed issue in `model_parameters.rma()` when model contained interaction
  terms.

* Fixed sign error for `model_parameters.htest()` for objects created with
  `t.test.formula()` (issue #552)

* Fixed issue when computing random effect variances in `model_parameters()` for
  mixed models with categorical random slopes.

# parameters 0.14.0

## Breaking changes

* `check_sphericity()` has been renamed into `check_sphericity_bartlett()`.

* Removed deprecated arguments.

* `model_parameters()` for bootstrapped samples used in *emmeans* now treats the
  bootstrap samples as samples from posterior distributions (Bayesian models).

## New supported model classes

* `SemiParBIV` (*GJRM*), `selection` (*sampleSelection*), `htest` from the
  *survey* package, `pgmm` (*plm*).

## General

* Performance improvements for models from package *survey*.

## New functions

* Added a `summary()` method for `model_parameters()`, which is a convenient
  shortcut for `print(..., select = "minimal")`.

## Changes to functions

### `model_parameters()`

* `model_parameters()` gains a `parameters` argument, which takes a regular
  expression as string, to select specific parameters from the returned data
  frame.

* `print()` for `model_parameters()` and `compare_parameters()` gains a `groups`
  argument, to group parameters in the output. Furthermore, `groups` can be used
  directly as argument in `model_parameters()` and `compare_parameters()` and
  will be passed to the `print()` method.

* `model_parameters()` for ANOVAs now saves the type as attribute and prints
  this information as footer in the output as well.

* `model_parameters()` for *htest*-objects now saves the alternative hypothesis
  as attribute and prints this information as footer in the output as well.

* `model_parameters()` passes arguments `type`, `parallel` and `n_cpus` down to
  `bootstrap_model()` when `bootstrap = TRUE`.

### other

* `bootstrap_models()` for *merMod* and *glmmTMB* objects gains further
  arguments to set the type of bootstrapping and to allow parallel computing.

* `bootstrap_parameters()` gains the `ci_method` type `"bci"`, to compute
  bias-corrected and accelerated bootstrapped intervals.

* `ci()` for `svyglm` gains a `method` argument.

## Bug fixes

* Fixed issue in `model_parameters()` for *emmGrid* objects with Bayesian
  models.

* Arguments `digits`, `ci_digits` and `p_digits` were ignored for `print()` and
  only worked when used in the call to `model_parameters()` directly.

# parameters 0.13.0

## General

* Revised and improved the `print()` method for `model_parameters()`.

## New supported model classes

* `blrm` (*rmsb*), `AKP`, `med1way`, `robtab` (*WRS2*), `epi.2by2` (*epiR*),
  `mjoint` (*joineRML*), `mhurdle` (*mhurdle*), `sarlm` (*spatialreg*),
  `model_fit` (*tidymodels*), `BGGM` (*BGGM*), `mvord` (*mvord*)

## Changes to functions

### `model_parameters()`

* `model_parameters()` for `blavaan` models is now fully treated as Bayesian
  model and thus relies on the functions from *bayestestR* (i.e. ROPE, Rhat or
  ESS are reported) .

* The `effects`-argument from `model_parameters()` for mixed models was revised
  and now shows the random effects variances by default (same functionality as
  `random_parameters()`, but mimicking the behaviour from
  `broom.mixed::tidy()`). When the `group_level` argument is set to `TRUE`, the
  conditional modes (BLUPs) of the random effects are shown.

* `model_parameters()` for mixed models now returns an `Effects` column even
  when there is just one type of "effects", to mimic the behaviour from
  `broom.mixed::tidy()`. In conjunction with `standardize_names()` users can get
  the same column names as in `tidy()` for `model_parameters()` objects.

* `model_parameters()` for t-tests now uses the group values as column names.

* `print()` for `model_parameters()` gains a `zap_small` argument, to avoid
  scientific notation for very small numbers. Instead, `zap_small` forces to
  round to the specified number of digits.

* To be internally consistent, the degrees of freedom column for `lqm(m)` and
  `cgam(m)` objects (with *t*-statistic) is called `df_error`.

* `model_parameters()` gains a `summary` argument to add summary information
  about the model to printed outputs.

* Minor improvements for models from *quantreg*.

* `model_parameters` supports rank-biserial, rank epsilon-squared, and Kendall's
  *W* as effect size measures for `wilcox.test()`, `kruskal.test`, and
  `friedman.test`, respectively.

### Other functions

* `describe_distribution()` gets a `quartiles` argument to include 25th and 75th
  quartiles of a variable.

## Bug fixes

* Fixed issue with non-initialized argument `style` in `display()` for
  `compare_parameters()`.

* Make `print()` for `compare_parameters()` work with objects that have "simple"
  column names for confidence intervals with missing CI-level (i.e. when column
  is named `"CI"` instead of, say, `"95% CI"`).

* Fixed issue with `p_adjust` in `model_parameters()`, which did not work for
  adjustment-methods `"BY"` and `"BH"`.

* Fixed issue with `show_sigma` in `print()` for `model_parameters()`.

* Fixed issue in `model_parameters()` with incorrect order of degrees of
  freedom.

# parameters 0.12.0

## General

* Roll-back R dependency to R >= 3.4.

* Bootstrapped estimates (from `bootstrap_model()` or `bootstrap_parameters()`)
  can be passed to `emmeans` to obtain bootstrapped estimates, contrasts, simple
  slopes (etc) and their CIs.

  * These can then be passed to `model_parameters()` and related functions to
    obtain standard errors, p-values, etc.

## Breaking changes

* `model_parameters()` now always returns the confidence level for as additional
  `CI` column.

* The `rule` argument in `equivalenct_test()` defaults to `"classic"`.

## New supported model classes

* `crr` (*cmprsk*), `leveneTest()` (*car*), `varest` (*vars*), `ergm` (*ergm*),
  `btergm` (*btergm*), `Rchoice` (*Rchoice*), `garch` (*tseries*)

## New functions

* `compare_parameters()` (and its alias `compare_models()`) to show / print
  parameters of multiple models in one table.

## Changes to functions

* Estimation of bootstrapped *p*-values has been re-written to be more
  accurate.

* `model_parameters()` for mixed models gains an `effects`-argument, to return
  fixed, random or both fixed and random effects parameters.

* Revised printing for `model_parameters()` for *metafor* models.

* `model_parameters()` for *metafor* models now recognized confidence levels
  specified in the function call (via argument `level`).

* Improved support for effect sizes in `model_parameters()` from *anova*
  objects.

## Bug fixes

* Fixed edge case when formatting parameters from polynomial terms with many
  degrees.

* Fixed issue with random sampling and dropped factor levels in
  `bootstrap_model()`.

# parameters 0.11.0

## New supported model classes

* `coxr` (*coxrobust*), `coeftest` (*lmtest*), `ivfixed` (*ivfixed*), `ivprobit`
  (*ivprobit*), `riskRegression` (*riskRegression*), `fitdistr` (*MASS*),
  `yuen`, `t1way`, `onesampb`, `mcp1` and `mcp2` (*WRS2*), `Anova.mlm` (*car*),
  `rqs` (*quantreg*), `lmodel2` (*lmodel2*), `summary.lm`, `PMCMR`, `osrt` and
  `trendPMCMR` (*PMCMRplus*), `bamlss` (*bamlss*).

## New functions

### Printing and table Formatting

* `print_html()` as an alias for `display(format = "html")`. This allows to
  print tabular outputs from data frames (as returned by most functions in
  _parameters_) into nicely rendered HTML markdown tables.

## Changes to functions

* Added more effect size measures to `model_parameters()` for `htest` objects.

* `model_parameters()` for anova objects gains a `power` argument, to calculate
  the power for each parameter.

* `ci()` for models from *lme4* and *glmmTMB* can now computed profiled
  confidence intervals, using `method = "profile"`. Consequently,
  `model_parameters()` with `df_method = "profile"` also computes profiled
  confidence intervals. For models of class `glmmTMB`, option `"uniroot"` is
  also available.

## Bug fixes

* `model_parameters()` for t-tests when `standardize_d = TRUE`, did not return
  columns for the group-specific means.

* Fixed issue in `p_value()` for `fixest::feols()`.

* Fixed issue in `model_parameters()` for `glmer()` models with p-values that
  were calculated with `df_method = "ml1"` or `df_method = "betwithin"`.

* Fixed issue in `model_parameters()` for multinomial models when response was a
  character vector (and no factor).

* Fixed issue in `print_md()` for model-parameters objects from Bayesian
  models.

* Fixed issues with printing of model parameters for multivariate response
  models from *brms*.

* Fixed issue with paired t-tests and `model_parameters()`.

# parameters 0.10.1

## New functions

* `format_p_adjust()`, to create pretty names for p-adjustment methods.

## Bug fixes

* Fixed breaking code / failing tests due to latest _effectsize_ update.

* Fixed issue with `model_parameters()` for models of class `mlm`.

* Undocumented arguments `digits`, `ci_digits` and `p_digits` worked for
  `print()`, but not when directly called inside `model_parameters()`. Now,
  `model_parameters(model, digits = 5, ci_digits = 8)` works again.

* Fixed some minor printing-issues.

# parameters 0.10.0

## Breaking changes

* The default-method for effect sizes in `model_parameters()` for Anova-models
  (i.e. when arguments `omega_squared`, `eta_squared` or `epsilon_squared` are
  set to `TRUE`) is now `"partial"`, as initially intended.

* Column names for degrees of freedom were revised. `"df_residual"` was replaced
  by the more generic `"df_error"`. Moreover, models of class `htest` now also
  have the column name `"df_error"` and no longer `"df"` (where applicable).

* Some re-exports for functions that were moved to *insight* longer ago, were
  now removed.

## New supported model classes

* `Glm` (*rms*), `mediate` (*mediation*).

* `model_parameters()` supports `Gam` models (*gam*), `ridgelm` (*MASS*),
  `htest` objects from `oneway.test()`, `chisq.test()`, `prop.test()`,
  `mcnemar.test()` and `pairwise.htest` objects, `mcmc.list` (e.g. from
  *bayesGARCH*).

## New functions

### Printing and table Formatting

* `display()`, to format output from package-functions into different formats.

* `print_md()` as an alias for `display(format = "markdown")`. This allows to
  print tabular outputs from data frames (as returned by most functions in
  _parameters_) into nicely rendered markdown tables.

* `format()`, to create a "pretty data frame" with nicer column names and
  formatted values. This is one of the worker-functions behind `print()` or
  `print_md()`.

## Changes to functions

### `model_parameters()`

* `model_parameters()` for Anova-models (of class `aov`, `anova` etc.) gains a
  `ci`-argument, to add confidence intervals to effect size parameters.

* `model_parameters()` for `htest` objects gains a `cramers_v` and `phi`
  argument, to compute effect size parameters for objects from `chisq.test()`,
  and a `standardized_D` argument, to compute effect size parameters for objects
  from `t.test()`.

* `model_parameters()` for `metafor`-models is more stable when called from
  inside functions.

* `model_parameters()` for *metaBMA*-models now includes prior information for
  the meta-parameters.

* `model_parameters()` for meta-analysis-models gains a
  `include_studies`-argument, to include or remove studies from the output.

* `model_parameters()` for gam-models now includes the residual df for smooth
  terms, and no longer the reference df.

* Slightly revised and improved the `print()` method for `model_parameters()`.

### Other functions

* `describe_distribution()` now includes the name of the centrality index in the
  `CI`-column, when `centrality = "all"`.

* `pool_parameters()` gains a `details`-argument. For mixed models, and if
  `details = TRUE`, random effect variances will also be pooled.

## Bug fixes

* Fixed issue in `ci()` for *lme* models with non-positive definite
  variance-covariance.

* Fixed issue in `model_parameters()` for `nnet::multinom()`, `lqmm::lqm()`,
  `mgcv::gam()`, and `margins::margins()` models, and models from package
  *blme*.
## Test environments
* local R installation, R 4.0.5
* ubuntu 16.04 (on github-actions ci), R 4.0.5
* win-builder (on github-actions ci)

## R CMD check results

0 errors | 0 warnings | 0 note

## revdepcheck results

We checked 19 reverse dependencies, comparing R CMD check results across CRAN
and dev versions of this package.

 * We saw 1 new problems: *effectsize*
 * We failed to check 0 packages
 
We have already informed the maintainer of this breaking change and the new version of *effectsize* on GitHub should soon be submitted to CRAN.---
title: "Extracting, Computing and Exploring the Parameters of Statistical Models using R"
authors:
- affiliation: 1
  name: Daniel Lüdecke
  orcid: 0000-0002-8895-3206
- affiliation: 2
  name: Mattan S. Ben-Shachar
  orcid: 0000-0002-4287-4801
- affiliation: 3
  name: Indrajeet Patil
  orcid: 0000-0003-1995-6531
- affiliation: 4
  name: Dominique Makowski
  orcid: 0000-0001-5375-9967

date: "01 July 2020"
output: pdf_document
bibliography: paper.bib
csl: apa.csl
tags:
- R
- easystats
- parameters
- regression
- linear models
- coefficients
affiliations:
- index: 1
  name:  University Medical Center Hamburg-Eppendorf, Germany
- index: 2
  name: Ben-Gurion University of the Negev, Israel
- index: 3
  name: Max Planck Institute for Human Development, Germany
- index: 4
  name: Nanyang Technological University, Singapore
---

# Summary

The recent growth of data science is partly fueled by the ever-growing amount of data and the joint important developments in statistical modeling, with new and powerful models and frameworks becoming accessible to users. Although there exist some generic functions to obtain model summaries and parameters, many package-specific modeling functions do not provide such methods to allow users to access such valuable information. 

# Aims of the Package

**parameters** is an R-package [@rcore] that fills this important gap. Its primary goal is to provide utilities for processing the parameters of various statistical models. Beyond computing p-values, standard errors, confidence intervals (CI), Bayesian indices and other measures for a wide variety of models, this package implements features like parameters bootstrapping and engineering (such as variables reduction and/or selection), as well as tools for data reduction like functions to perform cluster, factor or principal component analysis.

Another important goal of the **parameters** package is to facilitate and streamline the process of reporting results of statistical models, which includes the easy and intuitive calculation of standardized estimates in addition to robust standard errors and p-values. **parameters** therefor offers a simple and unified syntax to process a large variety of (model) objects from many different packages.

**parameters** is part of the [*easystats*](https://github.com/easystats/easystats) ecosystem, a collaborative project created to facilitate the usage of R for statistical analyses.

# Comparison to other Packages

**parameters** functionality is in part comparable to packages like **broom** [@robinson_broom_2020], **finalfit** [@harrison2020finalfit] or **stargazer** [@hlavac_stargazer_2018] (and maybe some more). Yet, there are some notable differences, e.g.:

- **broom** (via `glance()`), **finalfit** (via `ff_metrics()`) and **stargazer** (via `stargazer()`) report fit indices (such as R2 or AIC) by default, while **parameters** does not. However, there is a dedicated package in the *easystats* project for assessing regression model quality and fit indices, **performance** [@luedecke2020performance].
- **parameters** easily allows to compute standardized estimates, robust estimation, small-sample-size corrections for degrees of freedom (like *Satterthwaite* or *Kenward-Roger*), bootstrapping or simulating parameters, and feature reduction. Furthermore, **parameters** provides functions to test for the presence or absence of an effect [_equivalence testing_, see @lakens2020equivalence].
- For most functions, [easy-to-use `plot()`-methods](https://easystats.github.io/see/articles/parameters.html) exist to quickly create nice looking plots (powered by the **see** package [@ludecke2020see]).
- **parameters** is a very lightweight package. Its main functionality only relies on the **insight**, the **bayestestR**, and the **effectsize** packages [@ludecke2019insight; @makowski2019bayetestR; @benshachar2020effecsize] to access and process information contained in models, and these packages in turn only depend on R core packages. However, additional features that do not belong to the core functions of **parameters** require the installation of other packages, such as **sandwich** [@zeileis2006] for robust estimation, **psych** [@revelle_psych_2019] for factor analysis or PCA or **cAIC4** [@saefken_caic4_2018] for parameter selection for mixed models.

# Examples of Features

As stated above, **parameters** creates summary tables of many different statistical models. The workflow is simple: fit a model and pass it to the `model_parameters()` function (or its shortcut, `parameters()`) to obtain information about the model's parameters. 

![](figure1.png)

In the following, we show some brief examples. However, a comprehensive overview including in-depth examples are accessible via the dedicated website (https://easystats.github.io/parameters/).

## Summary of Model Parameters

`model_parameters()` allows you to extract the parameters and their characteristics from various models in a consistent way.

``` r
library(parameters)

model <- lm(Sepal.Length ~ Species, data = iris)
parameters(model)

#> Parameter            | Coefficient |   SE |       95% CI |      p
#> -----------------------------------------------------------------
#> (Intercept)          |        5.01 | 0.07 | [4.86, 5.15] | < .001
#> Species [versicolor] |        0.93 | 0.10 | [0.73, 1.13] | < .001
#> Species [virginica]  |        1.58 | 0.10 | [1.38, 1.79] | < .001
```

Extraction of robust indices is possible for many models, in particular models supported by the **sandwich** [@zeileis2006] and **clubSandwich** [@pustejovsky2020] packages.

``` r
parameters(model, robust = TRUE)

#> Parameter            | Coefficient |   SE |       95% CI |      p
#> -----------------------------------------------------------------
#> (Intercept)          |        5.01 | 0.05 | [4.91, 5.11] | < .001
#> Species [versicolor] |        0.93 | 0.09 | [0.75, 1.11] | < .001
#> Species [virginica]  |        1.58 | 0.10 | [1.38, 1.79] | < .001
```

For linear mixed models, `parameters()` also allows to specify the method for approximating degrees of freedom, which may improve the accurracy for calculated standard errors or p-values.

``` r
library(lme4)
model <- lmer(
  Sepal.Length ~ Sepal.Width * Petal.Length + (1 | Species), 
  data = iris
)

parameters(model, digits = 3)

#> Parameter                  | Coefficient |    SE |         95% CI |      p
#> --------------------------------------------------------------------------
#> (Intercept)                |       0.707 | 0.652 | [-0.57,  1.98] | 0.278 
#> Sepal.Width                |       0.731 | 0.156 | [ 0.43,  1.04] | < .001
#> Petal.Length               |       1.023 | 0.143 | [ 0.74,  1.30] | < .001
#> Sepal.Width * Petal.Length |      -0.084 | 0.040 | [-0.16, -0.01] | 0.035 

parameters(model, digits = 3, df_method = "kenward")

#> Parameter                  | Coefficient |    SE |         95% CI |      p
#> --------------------------------------------------------------------------
#> (Intercept)                |       0.707 | 0.654 | [-0.70,  2.11] | 0.298 
#> Sepal.Width                |       0.731 | 0.157 | [ 0.42,  1.04] | < .001
#> Petal.Length               |       1.023 | 0.145 | [ 0.74,  1.31] | < .001
#> Sepal.Width * Petal.Length |      -0.084 | 0.040 | [-0.16, -0.01] | 0.037 
```

## Visualisation

**parameters** functions also include plotting capabilities via the [**see** package](https://easystats.github.io/see/) [@ludecke2020see]. A complete overview of plotting functions is available at the *see* website (https://easystats.github.io/see/articles/parameters.html).

```r
library(see)

model <- lm(Sepal.Length ~ Petal.Width * Species, data=iris)
plot(parameters(model))
```

![](figure3.png)

# Licensing and Availability

**parameters** is licensed under the GNU General Public License (v3.0), with all source code stored at GitHub (https://github.com/easystats/parameters), and with a corresponding issue tracker for bug reporting and feature enhancements. In the spirit of honest and open science, we encourage requests/tips for fixes, feature updates, as well as general questions and concerns via direct interaction with contributors and developers.

# Acknowledgments

**parameters** is part of the collaborative [*easystats*](https://github.com/easystats/easystats) ecosystem. Thus, we would like to thank the [members of easystats](https://github.com/orgs/easystats/people) as well as the users.

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
reported to the community leaders responsible for enforcement at d.luedecke@uke.de. 
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
# Contributing to parameters

This outlines how to propose a change to **parameters**. 

## Fixing typos

Small typos or grammatical errors in documentation may be edited directly using the GitHub web interface, so long as the changes are made in the _source_ file. If you want to fix typos in the documentation, please edit the related `.R` file in the `R/` folder. Do _not_ edit an `.Rd` file in `man/`.

## Filing an issue

The easiest way to propose a change or new feature is to file an issue. If you've found a
bug, you may also create an associated issue. If possible, try to illustrate your proposal or the bug with a minimal [reproducible example](https://www.tidyverse.org/help/#reprex).

## Pull requests

*  Please create a Git branch for each pull request (PR).
*  Your contributed code should roughly follow the [R style guide](http://style.tidyverse.org), but in particular our [**easystats convention of code-style**](https://github.com/easystats/easystats#convention-of-code-style).
*  parameters uses [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html),
for documentation.
*  parameters uses [testthat](https://cran.r-project.org/package=testthat). Adding tests to the PR makes it easier for me to merge your PR into the code base.
*  If your PR is a user-visible change, you may add a bullet to the top of `NEWS.md` describing the changes made. You may optionally add your GitHub username, and links to relevant issue(s)/PR(s).

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to
abide by its terms.
## revdepcheck results

We checked 23 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 3 new problems
 * We failed to check 9 packages

Issues with CRAN packages are summarised below.

### New problems
(This reports the first line of each new failure)

* pairwiseComparisons
  checking tests ...

* psycho
  checking examples ... ERROR

* psycModel
  checking examples ... ERROR
  checking tests ...

### Failed to check

* correlation      (NA)
* effectsize       (NA)
* ggstatsplot      (NA)
* modelbased       (NA)
* report           (NA)
* see              (NA)
* sjPlot           (NA)
* sjstats          (NA)
* statsExpressions (NA)
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.1.0 (2021-05-18) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |German_Germany.1252          |
|ctype    |German_Germany.1252          |
|tz       |Europe/Berlin                |
|date     |2021-07-22                   |

# Dependencies

|package    |old    |new      |<U+0394>  |
|:----------|:------|:--------|:--|
|parameters |0.14.0 |0.14.0.1 |*  |
|bayestestR |0.10.0 |0.10.0   |   |
|datawizard |NA     |0.1.0    |*  |
|insight    |0.14.2 |0.14.2.1 |*  |

# Revdeps

## Failed to check (9)

|package                                          |version |error  |warning |note |
|:------------------------------------------------|:-------|:------|:-------|:----|
|[correlation](failures.md#correlation)           |0.6.1   |__+1__ |        |     |
|[effectsize](failures.md#effectsize)             |0.4.5   |__+1__ |        |     |
|[ggstatsplot](failures.md#ggstatsplot)           |0.8.0   |__+1__ |        |     |
|[modelbased](failures.md#modelbased)             |0.7.0   |__+1__ |        |1 -1 |
|[report](failures.md#report)                     |0.3.5   |__+1__ |        |     |
|[see](failures.md#see)                           |0.6.4   |__+1__ |        |1    |
|[sjPlot](failures.md#sjplot)                     |2.8.9   |__+1__ |        |1    |
|[sjstats](failures.md#sjstats)                   |0.18.1  |__+1__ |        |     |
|[statsExpressions](failures.md#statsexpressions) |1.1.0   |__+1__ |        |     |

## New problems (3)

|package                                                |version |error  |warning |note |
|:------------------------------------------------------|:-------|:------|:-------|:----|
|[pairwiseComparisons](problems.md#pairwisecomparisons) |3.1.6   |__+1__ |        |     |
|[psycho](problems.md#psycho)                           |0.6.1   |__+1__ |        |     |
|[psycModel](problems.md#psycmodel)                     |0.3.1   |__+2__ |        |1    |

# correlation

<details>

* Version: 0.6.1
* GitHub: https://github.com/easystats/correlation
* Source code: https://github.com/cran/correlation
* Date/Publication: 2021-04-09 06:10:02 UTC
* Number of recursive dependencies: 170

Run `revdep_details(, "correlation")` for more info

</details>

## Newly broken

*   checking whether package 'correlation' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/correlation/new/correlation.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'correlation' ...
** package 'correlation' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'correlation'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/correlation/new/correlation.Rcheck/correlation'


```
### CRAN

```
* installing *source* package 'correlation' ...
** package 'correlation' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (correlation)


```
# effectsize

<details>

* Version: 0.4.5
* GitHub: https://github.com/easystats/effectsize
* Source code: https://github.com/cran/effectsize
* Date/Publication: 2021-05-25 13:00:02 UTC
* Number of recursive dependencies: 266

Run `revdep_details(, "effectsize")` for more info

</details>

## Newly broken

*   checking whether package 'effectsize' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/effectsize/new/effectsize.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'effectsize' ...
** package 'effectsize' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'effectsize'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/effectsize/new/effectsize.Rcheck/effectsize'


```
### CRAN

```
* installing *source* package 'effectsize' ...
** package 'effectsize' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (effectsize)


```
# ggstatsplot

<details>

* Version: 0.8.0
* GitHub: https://github.com/IndrajeetPatil/ggstatsplot
* Source code: https://github.com/cran/ggstatsplot
* Date/Publication: 2021-06-09 11:10:02 UTC
* Number of recursive dependencies: 188

Run `revdep_details(, "ggstatsplot")` for more info

</details>

## Newly broken

*   checking whether package 'ggstatsplot' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/ggstatsplot/new/ggstatsplot.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'ggstatsplot' ...
** package 'ggstatsplot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'ggstatsplot'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/ggstatsplot/new/ggstatsplot.Rcheck/ggstatsplot'


```
### CRAN

```
* installing *source* package 'ggstatsplot' ...
** package 'ggstatsplot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (ggstatsplot)


```
# modelbased

<details>

* Version: 0.7.0
* GitHub: https://github.com/easystats/modelbased
* Source code: https://github.com/cran/modelbased
* Date/Publication: 2021-06-06 05:40:02 UTC
* Number of recursive dependencies: 215

Run `revdep_details(, "modelbased")` for more info

</details>

## Newly broken

*   checking whether package 'modelbased' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/modelbased/new/modelbased.Rcheck/00install.out' for details.
    ```

## Newly fixed

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: 'graphics'
      All declared Imports should be used.
    ```

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: 'glmmTMB'
    ```

## Installation

### Devel

```
* installing *source* package 'modelbased' ...
** package 'modelbased' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'modelbased'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/modelbased/new/modelbased.Rcheck/modelbased'


```
### CRAN

```
* installing *source* package 'modelbased' ...
** package 'modelbased' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (modelbased)


```
# pairwiseComparisons

<details>

* Version: 3.1.6
* GitHub: https://github.com/IndrajeetPatil/pairwiseComparisons
* Source code: https://github.com/cran/pairwiseComparisons
* Date/Publication: 2021-06-01 19:30:02 UTC
* Number of recursive dependencies: 87

Run `revdep_details(, "pairwiseComparisons")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
     ERROR
    Running the tests in 'tests/testthat.R' failed.
    Last 13 lines of output:
        2. | \-rlang::list2(...)
        3. +-pairwiseComparisons::pairwise_comparisons(...)
        4. | \-`%>%`(...)
        5. +-pairwiseComparisons:::tidy_model_parameters(.)
        6. | +-`%>%`(...)
        7. | +-parameters::model_parameters(model, verbose = FALSE, ...)
        8. | \-parameters:::model_parameters.mcp2(model, verbose = FALSE, ...)
        9. |   \-parameters:::.extract_wrs2_mcp12(model)
       10. |     +-base::`$<-`(`*tmp*`, "CI", value = numeric(0))
       11. |     \-base::`$<-.data.frame`(`*tmp*`, "CI", value = numeric(0))
       12. \-parameters::standardize_names(., style = "broom")
      
      [ FAIL 4 | WARN 0 | SKIP 1 | PASS 15 ]
      Error: Test failures
      Execution halted
    ```

# psycho

<details>

* Version: 0.6.1
* GitHub: https://github.com/neuropsychology/psycho.R
* Source code: https://github.com/cran/psycho
* Date/Publication: 2021-01-19 06:40:10 UTC
* Number of recursive dependencies: 75

Run `revdep_details(, "psycho")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in 'psycho-Ex.R' failed
    The error most likely occurred in:
    
    > ### Name: is.standardized
    > ### Title: Check if a dataframe is standardized.
    > ### Aliases: is.standardized
    > 
    > ### ** Examples
    > 
    > library(psycho)
    > library(effectsize)
    Error: package or namespace load failed for 'effectsize':
     object 'check_heterogeneity' is not exported by 'namespace:parameters'
    Execution halted
    ```

# psycModel

<details>

* Version: 0.3.1
* GitHub: NA
* Source code: https://github.com/cran/psycModel
* Date/Publication: 2021-05-21 23:50:03 UTC
* Number of recursive dependencies: 170

Run `revdep_details(, "psycModel")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in 'psycModel-Ex.R' failed
    The error most likely occurred in:
    
    > ### Name: cor_test
    > ### Title: Correlation table
    > ### Aliases: cor_test
    > 
    > ### ** Examples
    > 
    > cor_test(iris, where(is.numeric))
    Error in cor_test(iris, where(is.numeric)) : 
      please install.packages('correlation')
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in 'tests/testthat.R' failed.
    Last 13 lines of output:
        4. |     +-testthat:::.capture(...)
        5. |     | \-base::withCallingHandlers(...)
        6. |     \-rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
        7. +-testthat::expect_warning(...)
        8. | \-testthat:::expect_condition_matching(...)
        9. |   \-testthat:::quasi_capture(...)
       10. |     +-testthat:::.capture(...)
       11. |     | \-base::withCallingHandlers(...)
       12. |     \-rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
       13. \-psycModel::integrated_model_summary(...)
       14.   \-psycModel::model_summary(...)
      
      [ FAIL 1 | WARN 0 | SKIP 4 | PASS 13 ]
      Error: Test failures
      Execution halted
    ```

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: 'lifecycle'
      All declared Imports should be used.
    ```

# report

<details>

* Version: 0.3.5
* GitHub: https://github.com/easystats/report
* Source code: https://github.com/cran/report
* Date/Publication: 2021-06-10 11:50:02 UTC
* Number of recursive dependencies: 212

Run `revdep_details(, "report")` for more info

</details>

## Newly broken

*   checking whether package 'report' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/report/new/report.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'report' ...
** package 'report' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'report'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/report/new/report.Rcheck/report'


```
### CRAN

```
* installing *source* package 'report' ...
** package 'report' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (report)


```
# see

<details>

* Version: 0.6.4
* GitHub: https://github.com/easystats/see
* Source code: https://github.com/cran/see
* Date/Publication: 2021-05-29 09:00:03 UTC
* Number of recursive dependencies: 215

Run `revdep_details(, "see")` for more info

</details>

## Newly broken

*   checking whether package 'see' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/see/new/see.Rcheck/00install.out' for details.
    ```

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: 'glmmTMB'
    ```

## Installation

### Devel

```
* installing *source* package 'see' ...
** package 'see' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'see'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/see/new/see.Rcheck/see'


```
### CRAN

```
* installing *source* package 'see' ...
** package 'see' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (see)


```
# sjPlot

<details>

* Version: 2.8.9
* GitHub: https://github.com/strengejacke/sjPlot
* Source code: https://github.com/cran/sjPlot
* Date/Publication: 2021-07-10 10:30:02 UTC
* Number of recursive dependencies: 217

Run `revdep_details(, "sjPlot")` for more info

</details>

## Newly broken

*   checking whether package 'sjPlot' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjPlot/new/sjPlot.Rcheck/00install.out' for details.
    ```

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: 'glmmTMB'
    ```

## Installation

### Devel

```
* installing *source* package 'sjPlot' ...
** package 'sjPlot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'sjPlot'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjPlot/new/sjPlot.Rcheck/sjPlot'


```
### CRAN

```
* installing *source* package 'sjPlot' ...
** package 'sjPlot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (sjPlot)


```
# sjstats

<details>

* Version: 0.18.1
* GitHub: https://github.com/strengejacke/sjstats
* Source code: https://github.com/cran/sjstats
* Date/Publication: 2021-01-09 13:50:02 UTC
* Number of recursive dependencies: 211

Run `revdep_details(, "sjstats")` for more info

</details>

## Newly broken

*   checking whether package 'sjstats' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjstats/new/sjstats.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'sjstats' ...
** package 'sjstats' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'sjstats'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjstats/new/sjstats.Rcheck/sjstats'


```
### CRAN

```
* installing *source* package 'sjstats' ...
** package 'sjstats' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (sjstats)


```
# statsExpressions

<details>

* Version: 1.1.0
* GitHub: https://github.com/IndrajeetPatil/statsExpressions
* Source code: https://github.com/cran/statsExpressions
* Date/Publication: 2021-05-30 04:30:02 UTC
* Number of recursive dependencies: 166

Run `revdep_details(, "statsExpressions")` for more info

</details>

## Newly broken

*   checking whether package 'statsExpressions' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/statsExpressions/new/statsExpressions.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'statsExpressions' ...
** package 'statsExpressions' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'statsExpressions'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/statsExpressions/new/statsExpressions.Rcheck/statsExpressions'


```
### CRAN

```
* installing *source* package 'statsExpressions' ...
** package 'statsExpressions' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (statsExpressions)


```
# correlation

<details>

* Version: 0.6.1
* GitHub: https://github.com/easystats/correlation
* Source code: https://github.com/cran/correlation
* Date/Publication: 2021-04-09 06:10:02 UTC
* Number of recursive dependencies: 170

Run `revdep_details(, "correlation")` for more info

</details>

## Newly broken

*   checking whether package 'correlation' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/correlation/new/correlation.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'correlation' ...
** package 'correlation' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'correlation'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/correlation/new/correlation.Rcheck/correlation'


```
### CRAN

```
* installing *source* package 'correlation' ...
** package 'correlation' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (correlation)


```
# effectsize

<details>

* Version: 0.4.5
* GitHub: https://github.com/easystats/effectsize
* Source code: https://github.com/cran/effectsize
* Date/Publication: 2021-05-25 13:00:02 UTC
* Number of recursive dependencies: 266

Run `revdep_details(, "effectsize")` for more info

</details>

## Newly broken

*   checking whether package 'effectsize' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/effectsize/new/effectsize.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'effectsize' ...
** package 'effectsize' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'effectsize'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/effectsize/new/effectsize.Rcheck/effectsize'


```
### CRAN

```
* installing *source* package 'effectsize' ...
** package 'effectsize' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (effectsize)


```
# ggstatsplot

<details>

* Version: 0.8.0
* GitHub: https://github.com/IndrajeetPatil/ggstatsplot
* Source code: https://github.com/cran/ggstatsplot
* Date/Publication: 2021-06-09 11:10:02 UTC
* Number of recursive dependencies: 188

Run `revdep_details(, "ggstatsplot")` for more info

</details>

## Newly broken

*   checking whether package 'ggstatsplot' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/ggstatsplot/new/ggstatsplot.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'ggstatsplot' ...
** package 'ggstatsplot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'ggstatsplot'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/ggstatsplot/new/ggstatsplot.Rcheck/ggstatsplot'


```
### CRAN

```
* installing *source* package 'ggstatsplot' ...
** package 'ggstatsplot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (ggstatsplot)


```
# modelbased

<details>

* Version: 0.7.0
* GitHub: https://github.com/easystats/modelbased
* Source code: https://github.com/cran/modelbased
* Date/Publication: 2021-06-06 05:40:02 UTC
* Number of recursive dependencies: 215

Run `revdep_details(, "modelbased")` for more info

</details>

## Newly broken

*   checking whether package 'modelbased' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/modelbased/new/modelbased.Rcheck/00install.out' for details.
    ```

## Newly fixed

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: 'graphics'
      All declared Imports should be used.
    ```

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: 'glmmTMB'
    ```

## Installation

### Devel

```
* installing *source* package 'modelbased' ...
** package 'modelbased' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'modelbased'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/modelbased/new/modelbased.Rcheck/modelbased'


```
### CRAN

```
* installing *source* package 'modelbased' ...
** package 'modelbased' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (modelbased)


```
# report

<details>

* Version: 0.3.5
* GitHub: https://github.com/easystats/report
* Source code: https://github.com/cran/report
* Date/Publication: 2021-06-10 11:50:02 UTC
* Number of recursive dependencies: 212

Run `revdep_details(, "report")` for more info

</details>

## Newly broken

*   checking whether package 'report' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/report/new/report.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'report' ...
** package 'report' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'report'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/report/new/report.Rcheck/report'


```
### CRAN

```
* installing *source* package 'report' ...
** package 'report' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (report)


```
# see

<details>

* Version: 0.6.4
* GitHub: https://github.com/easystats/see
* Source code: https://github.com/cran/see
* Date/Publication: 2021-05-29 09:00:03 UTC
* Number of recursive dependencies: 215

Run `revdep_details(, "see")` for more info

</details>

## Newly broken

*   checking whether package 'see' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/see/new/see.Rcheck/00install.out' for details.
    ```

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: 'glmmTMB'
    ```

## Installation

### Devel

```
* installing *source* package 'see' ...
** package 'see' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'see'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/see/new/see.Rcheck/see'


```
### CRAN

```
* installing *source* package 'see' ...
** package 'see' successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (see)


```
# sjPlot

<details>

* Version: 2.8.9
* GitHub: https://github.com/strengejacke/sjPlot
* Source code: https://github.com/cran/sjPlot
* Date/Publication: 2021-07-10 10:30:02 UTC
* Number of recursive dependencies: 217

Run `revdep_details(, "sjPlot")` for more info

</details>

## Newly broken

*   checking whether package 'sjPlot' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjPlot/new/sjPlot.Rcheck/00install.out' for details.
    ```

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: 'glmmTMB'
    ```

## Installation

### Devel

```
* installing *source* package 'sjPlot' ...
** package 'sjPlot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'sjPlot'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjPlot/new/sjPlot.Rcheck/sjPlot'


```
### CRAN

```
* installing *source* package 'sjPlot' ...
** package 'sjPlot' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (sjPlot)


```
# sjstats

<details>

* Version: 0.18.1
* GitHub: https://github.com/strengejacke/sjstats
* Source code: https://github.com/cran/sjstats
* Date/Publication: 2021-01-09 13:50:02 UTC
* Number of recursive dependencies: 211

Run `revdep_details(, "sjstats")` for more info

</details>

## Newly broken

*   checking whether package 'sjstats' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjstats/new/sjstats.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'sjstats' ...
** package 'sjstats' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'sjstats'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/sjstats/new/sjstats.Rcheck/sjstats'


```
### CRAN

```
* installing *source* package 'sjstats' ...
** package 'sjstats' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (sjstats)


```
# statsExpressions

<details>

* Version: 1.1.0
* GitHub: https://github.com/IndrajeetPatil/statsExpressions
* Source code: https://github.com/cran/statsExpressions
* Date/Publication: 2021-05-30 04:30:02 UTC
* Number of recursive dependencies: 166

Run `revdep_details(, "statsExpressions")` for more info

</details>

## Newly broken

*   checking whether package 'statsExpressions' can be installed ... ERROR
    ```
    Installation failed.
    See 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/statsExpressions/new/statsExpressions.Rcheck/00install.out' for details.
    ```

## Installation

### Devel

```
* installing *source* package 'statsExpressions' ...
** package 'statsExpressions' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error: object 'check_heterogeneity' is not exported by 'namespace:parameters'
Execution halted
ERROR: lazy loading failed for package 'statsExpressions'
* removing 'C:/Users/mail/Documents/R/easystats/parameters/revdep/checks/statsExpressions/new/statsExpressions.Rcheck/statsExpressions'


```
### CRAN

```
* installing *source* package 'statsExpressions' ...
** package 'statsExpressions' successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (statsExpressions)


```
---
output: 
  github_document:
    toc: false
    fig_width: 10.08
    fig_height: 6
tags: [r, reports]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  tidy.opts = list(width.cutoff = 100),
  fig.path = "man/figures/",
  comment = "#>"
)
options(knitr.kable.NA = '',
        digits = 1,
        width = 100)

set.seed(333)
library(parameters)
```

# parameters <img src="man/figures/logo.png" align="right" height="139" />

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02445/status.svg)](https://doi.org/10.21105/joss.02445)
[![downloads](http://cranlogs.r-pkg.org/badges/parameters)](https://cran.r-project.org/package=parameters)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/parameters)](https://cranlogs.r-pkg.org/) 
[![status](https://tinyverse.netlify.com/badge/parameters)](https://CRAN.R-project.org/package=parameters)

------------------------------------------------------------------------

:warning: For Bayesian models, we changed the default the CI width! Please make an [informed
decision](https://easystats.github.io/bayestestR/articles/credible_interval.html)
and set it explicitly (`ci = 0.89`, `ci = 0.95`, or anything else that
you decide) :warning:

------------------------------------------------------------------------

***Describe and understand your model's parameters!***

**parameters**' primary goal is to provide utilities for processing the parameters of various statistical models (see [here](https://easystats.github.io/insight/) for a list of supported models). Beyond computing *p-values*, *CIs*, *Bayesian indices* and other measures for a wide variety of models, this package implements features like *bootstrapping* of parameters and models, *feature reduction* (feature extraction and variable selection), or tools for data reduction like functions to perform cluster, factor or principal component analysis.

Another important goal of the **parameters** package is to facilitate and streamline the process of reporting results of statistical models, which includes the easy and intuitive calculation of standardized estimates or robust standard errors and p-values. **parameters** therefor offers a simple and unified syntax to process a large variety of (model) objects from many different packages.

## Installation

[![CRAN](http://www.r-pkg.org/badges/version/parameters)](https://cran.r-project.org/package=parameters) [![parameters status badge](https://easystats.r-universe.dev/badges/parameters)](https://easystats.r-universe.dev) [![R-check](https://github.com/easystats/parameters/workflows/R-check/badge.svg?branch=main)](https://github.com/easystats/parameters/actions)

Run the following to install the stable release of **parameters** from CRAN:

```{r, warning=FALSE, message=FALSE, eval=FALSE}
install.packages("parameters")
```

Or this one to install the latest development version from R-universe...

```{r, warning=FALSE, message=FALSE, eval=FALSE}
install.packages("parameters", repos = "https://easystats.r-universe.dev")
```

...or from GitHub:

```{r, warning=FALSE, message=FALSE, eval=FALSE}
install.packages("remotes")
remotes::install_github("easystats/parameters")
```

## Documentation

[![Documentation](https://img.shields.io/badge/documentation-parameters-orange.svg?colorB=E91E63)](https://easystats.github.io/parameters/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-parameters-orange.svg?colorB=2196F3)](https://easystats.github.io/parameters/reference/index.html) 

Click on the buttons above to access the package [documentation](https://easystats.github.io/parameters/) and the [easystats blog](https://easystats.github.io/blog/posts/), and check-out these vignettes:

- [Summary of Model Parameters](https://easystats.github.io/parameters/articles/model_parameters.html)
- [Standardized Model Parameters](https://easystats.github.io/parameters/articles/model_parameters_standardized.html)
- [Robust Estimation of Standard Errors, Confidence Intervals and p-values](https://easystats.github.io/parameters/articles/model_parameters_robust.html)
- [Model Parameters and Missing Data](https://easystats.github.io/parameters/articles/model_parameters_mice.html)
- [Feature reduction (PCA, cMDS, ICA...)](https://easystats.github.io/parameters/articles/parameters_reduction.html)
- [Structural models (EFA, CFA, SEM...)](https://easystats.github.io/parameters/articles/efa_cfa.html)
- [Parameters selection](https://easystats.github.io/parameters/articles/parameters_selection.html)
- [A Practical Guide for Panel Data Analysis](https://easystats.github.io/datawizard/articles/demean.html)

## Contributing and Support

In case you want to file an issue or contribute in another way to the package, please follow [this guide](https://github.com/easystats/parameters/blob/master/.github/CONTRIBUTING.md). For questions about the functionality, you may either contact us via email or also file an issue.


# Features
 
## Model's parameters description

```{r echo=FALSE, fig.width=734, fig.align='center', dpi=96}
knitr::include_graphics("man/figures/figure1.png")
```

The [`model_parameters()`](https://easystats.github.io/parameters/articles/model_parameters.html) function (that can be accessed via the `parameters()` shortcut) allows you to extract the parameters and their characteristics from various models in a consistent way. It can be considered as a lightweight alternative to [`broom::tidy()`](https://github.com/tidymodels/broom), with some notable differences: 

- The column names of the returned data frame are *specific* to their content. For instance, the column containing the statistic is named following the statistic name, i.e., *t*, *z*, etc., instead of a generic name such as *statistic* (however, you can get standardized (generic) column names using [`standardize_names()`](https://easystats.github.io/insight/reference/standardize_names.html)).
- It is able to compute or extract indices not available by default, such as *p-values*, *CIs*, etc.
- It includes *feature engineering* capabilities, including parameters [bootstrapping](https://easystats.github.io/parameters/reference/bootstrap_parameters.html).


### Classical Regression Models

```{r, warning=FALSE, message=FALSE}
model <- lm(Sepal.Width ~ Petal.Length * Species + Petal.Width, data = iris)

# regular model parameters
model_parameters(model)

# standardized parameters
model_parameters(model, standardize = "refit")
```

### Mixed Models

```{r, warning=FALSE, message=FALSE}
library(lme4)

model <- lmer(Sepal.Width ~ Petal.Length + (1|Species), data = iris)

# model parameters with CI, df and p-values based on Wald approximation
model_parameters(model, effects = "all")

# model parameters with CI, df and p-values based on Kenward-Roger approximation
model_parameters(model, df_method = "kenward")
```

### Structural Models

Besides many types of regression models and packages, it also works for other types of models, such as [**structural models**](https://easystats.github.io/parameters/articles/efa_cfa.html) (EFA, CFA, SEM...).

```{r, warning=FALSE, message=FALSE}
library(psych)

model <- psych::fa(attitude, nfactors = 3)
model_parameters(model)
```



## Variable and parameters selection

```{r echo=FALSE, fig.width=756, fig.align='center', dpi=96}
knitr::include_graphics("man/figures/figure2.png")
```

[`select_parameters()`](https://easystats.github.io/parameters/articles/parameters_selection.html) can help you quickly select and retain the most relevant predictors using methods tailored for the model type.

```{r, warning=FALSE, message=FALSE}
library(poorman)

lm(disp ~ ., data = mtcars) %>% 
  select_parameters() %>% 
  model_parameters()
```

## Miscellaneous

This packages also contains a lot of [other useful functions](https://easystats.github.io/parameters/reference/index.html):

### Describe a Distribution

```{r, warning=FALSE, message=FALSE}
data(iris)
describe_distribution(iris)
```

### Citation

In order to cite this package, please use the following command:

```{r, comment=""}
citation("parameters")
```
---
title: "Robust Estimation of Standard Errors, Confidence Intervals and p-values"
output: 
  github_document:
    toc: true
    fig_width: 10.08
    fig_height: 6
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, variable selection, feature selection]
vignette: >
  %\VignetteIndexEntry{Robust Estimation of Standard Errors}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = '')
options(digits = 2)
knitr::opts_chunk$set(comment = "#>")

if (!requireNamespace("poorman", quietly = TRUE) ||
    !requireNamespace("clubSandwich", quietly = TRUE) ||
    !requireNamespace("sandwich", quietly = TRUE) ||
    !requireNamespace("lme4", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(poorman)
  library(clubSandwich)
  library(lme4)
}

set.seed(333)
```

The [`model_parameters()`](https://easystats.github.io/parameters/articles/model_parameters.html) function also allows the computation of standard errors, confidence intervals and p-values based on robust covariance matrix estimation from model parameters. Robust estimation is based on the packages **sandwich** and **clubSandwich**, so all models supported by either of these packages work with `model_parameters()` when `robust = TRUE`.

## Classical Regression Models

### Robust Covariance Matrix Estimation from Model Parameters

By default, when `model_parameters(robust = TRUE)`, it internally calls `sandwich::vcovHC(type = "HC3")`. However, there are three arguments that allow for choosing different methods and options of robust estimation: `vcov_estimation`, `vcov_type` and `vcov_args` (see [`?standard_error_robust`](https://easystats.github.io/parameters/reference/standard_error_robust.html) for further details).

Let us start with a simple example, which uses a heteroskedasticity-consistent covariance matrix estimation with estimation-type "HC3" (i.e. `sandwich::vcovHC(type = "HC3")` is called):

```{r}
data(iris)
model <- lm(Petal.Length ~ Sepal.Length * Species + Sepal.Width, data = iris)

# model parameters, where SE, CI and p-values are based on robust estimation
mp <- model_parameters(model, robust = TRUE)
mp

# compare standard errors to result from sandwich-package
mp$SE
unname(sqrt(diag(sandwich::vcovHC(model))))
```

### Cluster-Robust Covariance Matrix Estimation (sandwich)

If another covariance matrix estimation is required, use the `vcov_estimation`-argument. This argument needs the suffix for the related `vcov*()`-functions as value, i.e. `vcov_estimation = "CL"` would call `sandwich::vcovCL()`, or `vcov_estimation = "HAC"` would call `sandwich::vcovHAC()`.

The specific estimation type can be changed with `vcov_type`. E.g., `sandwich::vcovCL()` accepts estimation types HC0 to HC3. In the next example, we use a clustered covariance matrix estimation with HC1-estimation type.

```{r}
# change estimation-type
mp <- model_parameters(model, robust = TRUE, vcov_estimation = "CL", vcov_type = "HC1")
mp

# compare standard errors to result from sandwich-package
mp$SE
unname(sqrt(diag(sandwich::vcovCL(model))))
```

Usually, clustered covariance matrix estimation is used when there is a cluster-structure in the data. The variable indicating the cluster-structure can be defined in `sandwich::vcovCL()` with the `cluster`-argument. In `model_parameters()`, additional arguments that should be passed down to functions from the **sandwich** package can be specified in `vcov_args`:

```{r}
iris$cluster <- factor(rep(LETTERS[1:8], length.out = nrow(iris)))
# change estimation-type, defining additional arguments
mp <- model_parameters(
  model, 
  robust = TRUE, 
  vcov_estimation = "CL", 
  vcov_type = "HC1",
  vcov_args = list(cluster = iris$cluster)
)
mp

# compare standard errors to result from sandwich-package
mp$SE
unname(sqrt(diag(sandwich::vcovCL(model, cluster = iris$cluster))))
```

### Cluster-Robust Covariance Matrix Estimation (clubSandwich)

Cluster-robust estimation of the variance-covariance matrix can also be achieved using `clubSandwich::vcovCR()`. Thus, when `vcov_estimation = "CR"`, the related function from the **clubSandwich** package is called. Note that this function _requires_ the specification of the `cluster`-argument.

```{r}
# create fake-cluster-variable, to demonstrate cluster robust standard errors
iris$cluster <- factor(rep(LETTERS[1:8], length.out = nrow(iris)))

# cluster-robust estimation
mp <- model_parameters(
  model, 
  robust = TRUE, 
  vcov_estimation = "CR", 
  vcov_type = "CR1", 
  vcov_args = list(cluster = iris$cluster)
)
mp

# compare standard errors to result from clubSsandwich-package
mp$SE
unname(sqrt(diag(clubSandwich::vcovCR(model, type = "CR1", cluster = iris$cluster))))
```


### Robust Covariance Matrix Estimation on Standardized Model Parameters

Finally, robust estimation can be combined with standardization. However, robust covariance matrix estimation only works for `standardize = "refit"`.

```{r}
# model parameters, robust estimation on standardized model
model_parameters(model, standardize = "refit", robust = TRUE)
```


## Mixed Models

### Robust Covariance Matrix Estimation for Mixed Models

For linear mixed models, that by definition have a clustered ("hierarchical" or multilevel) structure in the data, it is also possible to estimate a cluster-robust covariance matrix. This is possible due to the **clubSandwich** package, thus we need to define the same arguments as in the above example.

```{r}
library(lme4)
data(iris)
set.seed(1234)
iris$grp <- as.factor(sample(1:3, nrow(iris), replace = TRUE))

# fit example model
model <- lme4::lmer(
  Sepal.Length ~ Species * Sepal.Width + Petal.Length + (1 | grp),
  data = iris
)

# normal model parameters, like from 'summary()'
model_parameters(model)

# model parameters, cluster robust estimation for mixed models
model_parameters(
  model, 
  robust = TRUE, 
  vcov_estimation = "CR", 
  vcov_type = "CR1", 
  vcov_args = list(cluster = iris$grp)
)
```

### Robust Covariance Matrix Estimation on Standardized Mixed Model Parameters

Again, robust estimation can be combined with standardization for linear mixed models as well, which in such cases also only works for `standardize = "refit"`.

```{r}
# model parameters, cluster robust estimation on standardized mixed model
model_parameters(
  model, 
  standardize = "refit",
  robust = TRUE, 
  vcov_estimation = "CR", 
  vcov_type = "CR1", 
  vcov_args = list(cluster = iris$grp)
)
```
---
title: "Bootstrapped parameters"
output: 
  github_document:
    toc: true
    fig_width: 10.08
    fig_height: 6
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, bayesian, gam, smooth]
vignette: >
  %\VignetteIndexEntry{Bootstrapped parameters}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = '')
knitr::opts_chunk$set(comment=">")
options(digits=2)

if (!requireNamespace("ggplot2", quietly = TRUE) ||
  !requireNamespace("parameters", quietly = TRUE) ||
  !requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("tidyr", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
library(ggplot2)
library(poorman)
library(tidyr)
}

set.seed(333)
```


The basic idea of bootstrapping is that inference about a parent population from sample data can be modelled by resampling the sample data. It is often used as an alternative to statistical inference based on the assumption of a parametric model when that assumption is in doubt, or where parametric inference is complicated.

## Comparison between regular and boostrapped estimates

### Data

In order to compare point-estimates with bootstrapped parameters for frequentist models. We generated one large sample (the **parent population**, size `1000000`) of two continuous variables producing a regression coefficient of `0.5`. We then iteratively extracted a subsample of size `30`, computed 3 types of coefficient (regular, bootstrapped median with 1000 and 4000 iterations) that were substracted from the "parent" coefficient. The closer the value is from 0, and the closer it is from the "true" effect.

The data is available on githuband the code to generate it is available [here](https://easystats.github.io/circus/articles/bootstrapped.html).


```{r message=FALSE, warning=FALSE}
library(ggplot2)
library(poorman)
library(tidyr)

df <- read.csv("https://raw.github.com/easystats/circus/master/data/bootstrapped.csv")
```


### Visualisation
```{r message=FALSE, warning=FALSE, fig.height=10, fig.width=8}
library(see)

df_long <- df %>%
  select(Coefficient, Bootstrapped_1000, Bootstrapped_4000) %>%
  gather(Type, Distance) %>%
  mutate(Type = forcats::fct_relevel(Type, c("Coefficient", "Bootstrapped_1000", "Bootstrapped_4000")))

df_long %>% 
  ggplot(aes(y = Distance, x = Type, fill = Type)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#2196F3", "#FF9800", "#f44336")) +
  theme_modern() +
  ylab("Distance (0 is the parent / true effect)")
```


### Testing

#### Bayes factor analysis
```{r message=FALSE, warning=FALSE}
library(BayesFactor)
library(bayestestR)

bayestestR::bayesfactor(BayesFactor::ttestBF(df$Coefficient))
bayestestR::bayesfactor(BayesFactor::ttestBF(df$Bootstrapped_1000))
bayestestR::bayesfactor(BayesFactor::ttestBF(df$Bootstrapped_4000))
```

#### Contrast analysis
```{r message=FALSE, warning=FALSE}
library(emmeans)

lm(Distance ~ Type, data = df_long) %>% 
  emmeans::emmeans(~Type)
```


### Conclusion 

Negligible difference, but bootstrapped (n=1000) seems (very) slightly more accurate.
---
title: "Clustering with easystats"
output: 
  rmarkdown::html_vignette:
vignette: >
  %\VignetteIndexEntry{Clustering with easystats}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  eval = TRUE
)

if (!requireNamespace("NbClust", quietly = TRUE) ||
  !requireNamespace("mclust", quietly = TRUE) ||
  !requireNamespace("cluster", quietly = TRUE) ||
  !requireNamespace("ggplot2", quietly = TRUE) ||
  !requireNamespace("see", quietly = TRUE) ||
  !requireNamespace("fpc", quietly = TRUE) ||
  !requireNamespace("performance", quietly = TRUE) ||
  !requireNamespace("dbscan", quietly = TRUE) ||
  !requireNamespace("factoextra", quietly = TRUE) ||
  !requireNamespace("pvclust", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(NbClust)
  library(performance)
  library(cluster)
  library(fpc)
  library(mclust)
  library(insight)
  library(dbscan)
  library(factoextra)
  library(pvclust)
}
```

This vignette can be referred to by citing the package:

- Lüdecke, D., Ben-Shachar, M. S., Patil, I., \& Makowski, D. (2020). *Extracting, computing and exploring the parameters of statistical models using R*. Journal of Open Source Software, 5(53), 2445. https://doi.org/10.21105/joss.02445



Note that in order to fully use all the methods demonstrated below, you will need to additionally install the packages below:

```{r eval = FALSE}
install.packages(c("NbClust", "mclust", "pvclust", "cluster", "fpc", "dbscan"))
```

## Introduction

Clustering traditionally refers to the identification of groups of observations (i.e., data rows). It differs from methods like [**PCA or Factor Analysis**](https://easystats.github.io/parameters/articles/efa_cfa.html), which are usually applied on variables (i.e., columns). That said, it is possible to *transpose* your data (columns become rows) to apply clustering on variables.

There are many clustering algorithms (see [this for an overview](https://scikit-learn.org/stable/modules/clustering.html)), but they can grouped in two categories: **supervised** and **unsupervised** techniques. In **supervised** techniques, you have to explicitly specify [**how many clusters**](https://easystats.github.io/parameters/reference/n_clusters.html) you want to extract. **Unsupervised** techniques, on the other hand, will estimate this number as part of their algorithm. Note that there are no inherently superior and inferior clustering methods, each come with their sets of limitations and benefits.

As an example in the tutorial below, we will use the **iris** dataset, for which we know that there are 3 "real" clusters (the 3 Species of flowers). Let's first start with visualizing the 3 "real" clusters on a 2D space of the variables created through PCA. 


```{r}
library(ggplot2)
library(parameters)
library(see)

set.seed(33)  # Set random seed

# Select the first 4 numeric columns (drop the Species fator)
data <- iris[1:4]  
head(data)  # Print the 6 first rows

# Run PCA
pca <- principal_components(data, n = 2)
pca_scores <- predict(pca, names = c("PCA_1", "PCA_2"))
pca_scores$True_Clusters <- iris$Species  # Add real clusters

# Visualize
ggplot(pca_scores, aes(x = PCA_1, y = PCA_2, color = True_Clusters)) + 
  geom_point() +
  theme_modern()
```

While the **setosa** species stands out quite clearly in this PCA space, the separation between the two other species appear less clear cut. Let's see how data-driven clustering performs, and if we manage to retrieve these 3 clusters.


## Supervised Clustering Methods

### How Many Clusters to Extract?

There is no easy answer to that important question. The best way is to have strong expectations or hypotheses. If you don't, well, researchers have came up with data-driven solutions to estimate the optimal number of clusters. The problem is that there are now a lot of these numerical methods, and that they don't always agree...

Because there is no clearly better method, we have implemented in *easystats* a consensus-based algorithm that runs many of these methods, and returns the number of clusters that is the most agreed upon.


```{r}
n <- n_clusters(data, package = c("easystats", "NbClust", "mclust"))
n
plot(n)
```

As we can see, most methods suggest the existence of **2 clusters**, followed by a **3-clusters** solution. It seems like the data does not clearly discriminate between the 3 species of flowers. This discrepancy between what is, and what we can recover from real-world data, is a fundamental issue in data science.


### K-Means

We won't go too much into details about the mathematics and intuition behind these clustering methods, as good [resources](https://scikit-learn.org/stable/modules/clustering.html) are available all over the internet. Instead, we'll focus on how to apply them.

K-means is one of the most basic clustering algorithm, available in base R through the `kmeans()` function. However, we provide in easystats a unified function to run different clustering algorithms: [**cluster_analysis()**](https://easystats.github.io/parameters/reference/cluster_analysis.html). *(Note that k-means is a non-deterministic algorithm; running it multiple times will result in different results!)*

Now that we know how many clusters we want to extract (let's say that we have a strong hypothesis on 3, which is partially supported by the consensus method for estimating the optimal number of clusters).


```{r}
rez_kmeans <- cluster_analysis(data, n = 3, method = "kmeans")

rez_kmeans  # Show results
```

Note that we can also visualize the **centers** (i.e., the "average" of each variable for each cluster):


```{r}
plot(summary(rez_kmeans))  # Visualize cluster centers
```

One can extract the cluster assignments to use it as a new variable by using `predict()`.

```{r}
predict(rez_kmeans)  # Get clusters
```



### Hierarchical Clustering

Hierarchical clustering is also a common clustering algorithm, available in base R through the `hclust()` function. This method is a bit different in the sense that is does not straight up return clusters. Instead, in creates a hierarchical structure (a *dendrogram*), a tree from which we can *cut* branches to get a given number of clusters. Note that this "tree" cutting can be done in an unsupervised fashion too using bootstrapping (which we will apply in the next section).

```{r}
rez_hclust <- cluster_analysis(data, n = 3, method = "hclust")

rez_hclust  # Show results

# Visualize
plot(rez_hclust) + theme_modern()  # Visualize 
```

### Hierarchical K-Means

Hierarchical K-Means, as its name suggest, is essentially a combination of K-Means and hierarchical clustering that aims at improving the stability and robustness of the results.

```{r}
rez_hkmeans <- cluster_analysis(data, n = 3, method = "hkmeans")

rez_hkmeans  # Show results

# Visualize
plot(rez_hkmeans) + theme_modern()  # Visualize 
```

### K-Medoids (PAM)

Clustering around "medoids", instead of "centroid", is considered to be a more robust version of K-means. See `cluster::pam()` for more information. 

```{r}
rez_pam <- cluster_analysis(data, n = 3, method = "pam")

rez_pam  # Show results

# Visualize
plot(rez_pam) + theme_modern()  # Visualize 
```


## Unsupervised Clustering Methods

Unsupervised clustering methods estimate the optimal number of clusters themselves (hence, `n = NULL` as we don't pre-specify a given number of clusters). Note that unsupervised methods can sometimes identify observations that do not fit under any clusters (i.e., **"outliers"**). They will be classified as belonging to the cluster "0" (which is not a real cluster, but rather groups all the outliers).

### Bootstrapped Hierarchical Clustering

This method computes p-values for each cluster of the hierarchical cluster structure, and returns the **significant** clusters. This method can return a larger number of smaller clusters and, because it's based on bootstrapping, is quite slow.

```{r}
rez_hclust2 <- cluster_analysis(data, 
                        n = NULL, 
                        method = "hclust", 
                        iterations = 500,
                        ci = 0.90)

rez_hclust2  # Show results
plot(rez_hclust2) + theme_modern()  # Visualize
```


### DBSCAN

Although the DBSCAN method is quite powerful to identify clusters, it is highly dependent on its parameters, namely, `eps` and the `min_size`. Regarding the latter, the minimum size of any cluster is set by default to `0.1` (i.e., 10\% of rows), which is appropriate to avoid having too small clusters. 

The "optimal" **eps** value can be estimated using the [`n_clusters_dbscan()`](https://easystats.github.io/parameters/reference/cluster_analysis.html) function:

```{r}
eps <- n_clusters_dbscan(data, min_size = 0.1) 
eps
plot(eps)
```

It seems like the numeric method to find the elbow of the curve doesn't work well, and returns a value that is too high. Based on visual assessment, the elbow seems to be located around `eps = 1.45`.

```{r}
rez_dbscan <- cluster_analysis(data, method = "dbscan", dbscan_eps = 1.45)

rez_dbscan  # Show results
plot(rez_dbscan) + theme_modern()  # Visualize
```

### Hierarchical K-Means

Hierarchical DBSCAN is a variant that does not require the critical **EPS** argument. It computes the hierarchy of all DBSCAN solutions, and then finds the optimal cuts in the hierarchy using a stability-based extraction method.

```{r}
rez_hdbscan <- cluster_analysis(data, method = "hdbscan")

rez_hdbscan  # Show results

# Visualize
plot(rez_hdbscan) + theme_modern()  # Visualize 
```


### K-Medoids with estimation of number of clusters (pamk)

This is K-Medoids with an integrated estimation of the number of clusters. See `fpc::pamk` for more details.

```{r}
rez_pamk <- cluster_analysis(data, method = "pamk")

rez_pamk  # Show results

# Visualize
plot(rez_pamk) + theme_modern()  # Visualize 
```

### Mixture

Model-based clustering based on finite Gaussian mixture models. Models are estimated by EM algorithm initialized by hierarchical model-based agglomerative clustering. The optimal model is then selected according to BIC.

```{r}
library(mclust)

rez_mixture <- cluster_analysis(data, method = "mixture")

rez_mixture  # Show results

# Visualize
plot(rez_mixture) + theme_modern()  # Visualize 
```

## Metaclustering

One of the core "issue" of statistical clustering is that, in many cases, different methods will give different results. The **metaclustering** approach proposed by *easystats* (that finds echoes in *consensus clustering*; see Monti et al., 2003) consists of treating the unique clustering solutions as a ensemble, from which we can derive a probability matrix. This matrix contains, for each pair of observations, the probability of being in the same cluster. For instance, if the 6th and the 9th row of a dataframe has been assigned to a similar cluster by 5 our of 10 clustering methods, then its probability of being grouped together is 0.5.

Metaclustering is based on the hypothesis that, as each clustering algorithm embodies a different prism by which it sees the data, running an infinite amount of algorithms would result in the emergence of the "true" clusters. As the number of algorithms and parameters is finite, the probabilistic perspective is a useful proxy. This method is interesting where there is no obvious reasons to prefer one over another clustering method, as well as to investigate how robust some clusters are under different algorithms.


```{r}
list_of_results <- list(rez_kmeans, rez_hclust, rez_hkmeans, rez_pam,
                        rez_hclust2, rez_dbscan, rez_hdbscan, rez_mixture)

probability_matrix <- cluster_meta(list_of_results)

# Plot the matrix as a reordered heatmap
heatmap(probability_matrix, scale = "none", 
        col = grDevices::hcl.colors(256, palette = "inferno"))
```



The dendrogram (which is a **hierarchical clustering of the clustering solution**, hence the name of **meta**clustering), as well as the heatmap (in which the darker squares represent a higher probability of belonging to the same cluster) shows that there is one metacluster consisting of the 1-50 first rows (bottom left), and then the rest of the observations are closer to one another. However, two subclusters are still visible, corresponding to the "true" species.

The metaclustering approach confirms our initial hypothesis, *the **setosa** species stands out quite clearly, and the separation between the two other species is less clear cut*.

## Resources

- [Clustering algorithms overview](https://scikit-learn.org/stable/modules/clustering.html)
- [Density-based Clustering](https://www.datanovia.com/en/lessons/dbscan-density-based-clustering-essentials/)---
title: "Summary of Model Parameters"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, variable selection, feature selection]
vignette: >
  %\VignetteIndexEntry{Summary of Model Parameters}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%",
  tidy.opts = list(width.cutoff = 120)
)

if (!requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("effectsize", quietly = TRUE) ||
  !requireNamespace("BayesFactor", quietly = TRUE) ||
  !requireNamespace("lme4", quietly = TRUE) ||
  !requireNamespace("metafor", quietly = TRUE) ||
  !requireNamespace("lavaan", quietly = TRUE) ||
  !requireNamespace("nFactors", quietly = TRUE) ||
  !requireNamespace("EGAnet", quietly = TRUE) ||
  !requireNamespace("brms", quietly = TRUE) ||
  !requireNamespace("psych", quietly = TRUE) ||
  !requireNamespace("rstanarm", quietly = TRUE) ||
  !requireNamespace("glmmTMB", quietly = TRUE) ||
  !requireNamespace("GLMMadaptive", quietly = TRUE) ||
  !requireNamespace("FactoMineR", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(poorman)
  library(effectsize)
  library(EGAnet)
  library(psych)
  library(nFactors)
  library(brms)
  library(GLMMadaptive)
  library(FactoMineR)
  
}

set.seed(333)
```

The `model_parameters()` function (also accessible via the shortcut
`parameters()`) allows you to extract the parameters and their characteristics
from various models in a consistent way. It can be considered as a lightweight
alternative to [`broom::tidy()`](https://github.com/tidymodels/broom), with some
notable differences:

- The names of the returned data frame are **specific** to their content. For
  instance, the column containing the statistic is named following the statistic
  name, i.e., *t*, *z*, etc., instead of a generic name such as *statistic*
  (**however**, you can get standardized (generic) column names using
  [`standardize_names()`](https://easystats.github.io/insight/reference/standardize_names.html)).

- It is able to compute or extract indices not available by default, such as
  **p*-values**, **CIs**, etc.

- It includes **feature engineering** capabilities, including parameters
  [**bootstrapping**](https://easystats.github.io/parameters/reference/bootstrap_model.html).

## Correlations and *t*-tests

### Frequentist

```{r}
cor.test(iris$Sepal.Length, iris$Sepal.Width) %>%
  parameters()
```

```{r}
t.test(mpg ~ vs, data = mtcars) %>%
  parameters()
```


### Bayesian

```{r}
library(BayesFactor)

BayesFactor::correlationBF(iris$Sepal.Length, iris$Sepal.Width) %>%
  parameters()
```

```{r}
BayesFactor::ttestBF(formula = mpg ~ vs, data = mtcars) %>%
  parameters()
```

## ANOVAs

Indices of effect size for ANOVAs, such as partial and non-partial versions of
`eta_squared()`, `epsilon_sqared()` or `omega_squared()` are powered by the
[**effectsize**-package](https://easystats.github.io/effectsize/). However,
**parameters** uses these function to compute such indices for parameters
summaries, including confidence intervals

### Simple

```{r}
aov(Sepal.Length ~ Species, data = iris) %>%
  parameters(
    omega_squared = "partial",
    eta_squared = "partial",
    epsilon_squared = "partial"
  )
```

Let's complicate things further with an interaction term:

```{r}
aov(Sepal.Length ~ Species * Sepal.Width, data = iris) %>%
  parameters(
    omega_squared = "partial",
    eta_squared = "partial",
    ci = .8
  )
```

### Repeated measures

`parameters()` (resp. its alias `model_parameters()`) also works on repeated
measures ANOVAs, whether computed from `aov()` or from a mixed model.

```{r}
aov(mpg ~ am + Error(gear), data = mtcars) %>%
  parameters()
```

## Regressions (GLMs, Mixed Models, GAMs, ...)

`parameters()` (resp. its alias `model_parameters()`) was mainly built with
regression models in mind. It works for many types of models and packages,
including mixed models and Bayesian models.

### GLMs

```{r}
glm(vs ~ poly(mpg, 2) + cyl, data = mtcars, family = binomial()) %>%
  parameters()
```

```{r}
# show Odds Ratios and Wald-method for degrees of freedom
glm(vs ~ poly(mpg, 2) + cyl, data = mtcars, family = binomial()) %>%
  parameters(exponentiate = TRUE, df_method = "wald")
```

```{r}
# show Odds Ratios and include model summary
glm(vs ~ poly(mpg, 2) + cyl, data = mtcars, family = binomial()) %>%
  parameters(exponentiate = TRUE, summary = TRUE)
```

### Mixed Models

```{r}
library(lme4)

lmer(Sepal.Width ~ Petal.Length + (1 | Species), data = iris) %>%
  parameters()
```

### Mixed Models, without Random Effects Variances

```{r}
lmer(Sepal.Width ~ Petal.Length + (1 | Species), data = iris) %>%
  parameters(effects = "fixed")
```

### Mixed Model with Zero-Inflation Model

```{r}
library(GLMMadaptive)
library(glmmTMB)
data("Salamanders")
model <- mixed_model(
  count ~ spp + mined,
  random = ~ 1 | site,
  zi_fixed = ~ spp + mined,
  family = zi.negative.binomial(),
  data = Salamanders
)
parameters(model)
```

### Mixed Models with Dispersion Model

```{r}
library(glmmTMB)

sim1 <- function(nfac = 40, nt = 100, facsd = 0.1, tsd = 0.15, mu = 0, residsd = 1) {
  dat <- expand.grid(fac = factor(letters[1:nfac]), t = 1:nt)
  n <- nrow(dat)
  dat$REfac <- rnorm(nfac, sd = facsd)[dat$fac]
  dat$REt <- rnorm(nt, sd = tsd)[dat$t]
  dat$x <- rnorm(n, mean = mu, sd = residsd) + dat$REfac + dat$REt
  dat
}

set.seed(101)
d1 <- sim1(mu = 100, residsd = 10)
d2 <- sim1(mu = 200, residsd = 5)
d1$sd <- "ten"
d2$sd <- "five"
dat <- rbind(d1, d2)
model <- glmmTMB(x ~ sd + (1 | t), dispformula = ~sd, data = dat)

parameters(model)
```

### Bayesian Models

`model_parameters()` also works with Bayesian models from the **rstanarm**
package:

```{r}
library(rstanarm)

# if you are unfamiliar with the `refresh` argument here, it just avoids
# printing few messages to the console
stan_glm(mpg ~ wt * cyl, data = mtcars, refresh = 0) %>%
  parameters()
```

Additionally, it also works for models from the **brms** package.

For more complex models, specific model components can be printed using the
arguments `effects` and `component` arguments.

```{r}
library(brms)
data(fish)
set.seed(123)

# fitting a model using `brms`
model <- brm(
  bf(
    count ~ persons + child + camper + (1 | persons),
    zi ~ child + camper + (1 | persons)
  ),
  data = fish,
  family = zero_inflated_poisson(),
  refresh = 0
)

parameters(model, component = "conditional")

parameters(model, effects = "all", component = "all")
```

To include information about the random effect parameters (group levels), set
`group_level = TRUE`:

```{r}
parameters(model, effects = "all", component = "conditional", group_level = TRUE)
```

## Structural Models (PCA, EFA, CFA, SEM...)

The **parameters** package extends the support to structural models.

### Principal Component Analysis (PCA) and Exploratory Factor Analysis (EFA) 

```{r}
library(psych)

psych::pca(mtcars, nfactors = 3) %>%
  parameters()
```

We will avoid displaying a graph while carrying out factor analysis:

```{r}
library(FactoMineR)

FactoMineR::FAMD(iris, ncp = 3, graph = FALSE) %>%
  parameters()
```

### Confirmatory Factor Analysis (CFA) and Structural Equation Models (SEM)

#### Frequentist

```{r}
library(lavaan)

model <- lavaan::cfa(" visual  =~ x1 + x2 + x3
                       textual =~ x4 + x5 + x6
                       speed   =~ x7 + x8 + x9 ",
  data = HolzingerSwineford1939
)

model_parameters(model)
```

#### Bayesian

`blavaan` to be done.

## Meta-Analysis

`parameters()` also works for `rma`-objects from the **metafor** package.

```{r}
library(metafor)

mydat <- data.frame(
  effectsize = c(-0.393, 0.675, 0.282, -1.398),
  standarderror = c(0.317, 0.317, 0.13, 0.36)
)

rma(yi = effectsize, sei = standarderror, method = "REML", data = mydat) %>%
  model_parameters()
```

## Plotting Model Parameters

There is a `plot()`-method implemented in the
[**see**-package](https://easystats.github.io/see/). Several examples are shown
[in this vignette](https://easystats.github.io/see/articles/parameters.html).
---
title: "Robust Estimation of Standard Errors, Confidence Intervals, and p-values"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, variable selection, feature selection]
vignette: >
  %\VignetteIndexEntry{Robust Estimation of Standard Errors, Confidence Intervals, and p-values}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%"
)

if (!requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("clubSandwich", quietly = TRUE) ||
  !requireNamespace("sandwich", quietly = TRUE) ||
  !requireNamespace("lme4", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(poorman)
}

set.seed(333)
```

The
[`model_parameters()`](https://easystats.github.io/parameters/articles/model_parameters.html)
function also allows the computation of standard errors, confidence intervals,
and *p*-values based on robust covariance matrix estimation from model
parameters. Robust estimation relies on
the [`sandwich`](https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich-CL.pdf)
and
[`clubSandwich`](https://cran.r-project.org/web/packages/clubSandwich/index.html)
packages. This means that all models supported by either of these packages should work
with `model_parameters()` when `robust = TRUE`.

## Linear Regression Models

### Robust Covariance Matrix Estimation from Model Parameters

By default, when `model_parameters(robust = TRUE)`, it internally calls
`sandwich::vcovHC(type = "HC3")`. However, there are three arguments (see
[`?standard_error_robust`](https://easystats.github.io/parameters/reference/standard_error_robust.html)
for further details) that allow for choosing different methods and options of
robust estimation:
- `vcov_estimation`
- `vcov_type`
- `vcov_args`

Let us start with a simple example, which uses a heteroskedasticity-consistent
covariance matrix estimation with estimation-type "HC3" (i.e.
`sandwich::vcovHC(type = "HC3")`). 

First let's create a simple linear regression model, which we know violates
homoscedasticity assumption, and thus robust estimation methods are to be
considered.

```{r}
data(cars)
model <- lm(dist ~ speed, data = cars)

library(performance)
check_heteroscedasticity(model)
```

We would extract model parameters both with and without robust estimation to
highlight difference it makes to standard errors, confidence intervals, and
*p*-values. Also, note that the coefficient estimate and the *t*-statistic
associated with it remain unchanged.

```{r}
# model parameters, where SE, CI, and p-values are *not* based on robust estimation
model_parameters(model)

# model parameters, where SE, CI, and p-values are based on robust estimation
mp <- model_parameters(model, robust = TRUE)
mp

# compare standard errors to result from sandwich-package
mp$SE
unname(sqrt(diag(sandwich::vcovHC(model))))
```

### Cluster-Robust Covariance Matrix Estimation (sandwich)

If a different type of covariance matrix estimation is required, use the
`vcov_estimation`-argument. This argument needs the suffix for the related
`vcov*()`-functions as value, i.e. `vcov_estimation = "CL"` would call
`sandwich::vcovCL()`, or `vcov_estimation = "HAC"` would call
`sandwich::vcovHAC()`.

The specific estimation type can also be changed with `vcov_type`. For example,
`sandwich::vcovCL()` accepts estimation types `HC0` to `HC3`. In the next
example, we use a clustered covariance matrix estimation with `HC1`-estimation
type.

```{r}
# let's create a more complicated model
data(iris)
model <- lm(Petal.Length ~ Sepal.Length * Species + Sepal.Width, data = iris)

# change estimation-type
mp <- model_parameters(
  model,
  robust = TRUE,
  vcov_estimation = "CL", # type of covariance matrix
  vcov_type = "HC1" # type of robust estimation
)

mp

# compare standard errors to result from sandwich-package
mp$SE
unname(sqrt(diag(sandwich::vcovCL(model))))
```

Usually, clustered covariance matrix estimation is used when there is a
cluster-structure in the data. The variable indicating the cluster-structure can
be defined in `sandwich::vcovCL()` with the `cluster`-argument. In
`model_parameters()`, additional arguments that should be passed down to
functions from the `sandwich` package can be specified in `vcov_args`:

```{r}
iris$cluster <- factor(rep(LETTERS[1:8], length.out = nrow(iris)))

# change estimation-type, defining additional arguments
mp <- model_parameters(
  model,
  robust = TRUE,
  vcov_estimation = "CL",
  vcov_type = "HC1",
  vcov_args = list(cluster = iris$cluster)
)

mp

# compare standard errors to result from sandwich-package
mp$SE
unname(sqrt(diag(sandwich::vcovCL(model, cluster = iris$cluster))))
```

### Cluster-Robust Covariance Matrix Estimation (clubSandwich)

Cluster-robust estimation of the variance-covariance matrix can also be achieved
using `clubSandwich::vcovCR()`. Thus, when `vcov_estimation = "CR"`, the related
function from the `clubSandwich` package is called. Note that this function
_requires_ the specification of the `cluster`-argument.

```{r}
# create fake-cluster-variable, to demonstrate cluster robust standard errors
iris$cluster <- factor(rep(LETTERS[1:8], length.out = nrow(iris)))

# cluster-robust estimation
mp <- model_parameters(
  model,
  robust = TRUE,
  vcov_estimation = "CR",
  vcov_type = "CR1",
  vcov_args = list(cluster = iris$cluster)
)
mp

# compare standard errors to result from clubSsandwich-package
mp$SE

unname(sqrt(diag(clubSandwich::vcovCR(model, type = "CR1", cluster = iris$cluster))))
```

### Robust Covariance Matrix Estimation on Standardized Model Parameters

Finally, robust estimation can be combined with
[standardization](https://easystats.github.io/parameters/articles/model_parameters_standardized.html).
However, robust covariance matrix estimation only works for `standardize = "refit"`.

```{r}
# model parameters, robust estimation on standardized model
model_parameters(model, standardize = "refit", robust = TRUE)
```

## Linear Mixed-Effects Regression Models

### Robust Covariance Matrix Estimation for Mixed Models

For linear mixed-effects models, that by definition have a clustered
(*hierarchical* or *multilevel*) structure in the data, it is also possible to
estimate a cluster-robust covariance matrix. This is possible due to the
`clubSandwich` package, thus we need to define the same arguments as in the
above example.

```{r}
library(lme4)
data(iris)
set.seed(1234)
iris$grp <- as.factor(sample(1:3, nrow(iris), replace = TRUE))

# fit example model
model <- lme4::lmer(
  Sepal.Length ~ Species * Sepal.Width + Petal.Length + (1 | grp),
  data = iris
)

# model parameters without robust estimation
model_parameters(model)

# model parameters with cluster robust estimation
model_parameters(
  model,
  robust = TRUE,
  vcov_estimation = "CR",
  vcov_type = "CR1",
  vcov_args = list(cluster = iris$grp)
)
```

Notice that robust estimation returns different standard errors, confidence
intervals, and *p*-values compared to the standard estimation. Also, note that
the coefficient estimate and the statistic associated with it remain unchanged.

### Robust Covariance Matrix Estimation on Standardized Mixed Model Parameters

Once again, robust estimation can be combined with standardization for linear
mixed-effects models as well and works only with `standardize = "refit"`.

```{r}
# model parameters, cluster robust estimation of standardized mixed model
model_parameters(
  model,
  standardize = "refit",
  robust = TRUE,
  vcov_estimation = "CR",
  vcov_type = "CR1",
  vcov_args = list(cluster = iris$grp)
)
```

Notice how drastically some of the *p*-values change between
robust-unstandardized model and robust-standardized model.

<!-- TO DO: Maybe provide references to read more? -->
---
title: "Feature Reduction (PCA, cMDS, ICA...)"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, variable extraction, feature extraction, dimension extraction]
vignette: >
  %\VignetteIndexEntry{Feature Reduction (PCA, cMDS, ICA, ...)}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r , include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%"
)

if (!requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("nFactors", quietly = TRUE) ||
  !requireNamespace("EGAnet", quietly = TRUE) ||
  !requireNamespace("parameters", quietly = TRUE) ||
  !requireNamespace("insight", quietly = TRUE) ||
  !requireNamespace("psych", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(poorman)
  library(EGAnet)
  library(psych)
  library(nFactors)
}

set.seed(333)
```

Also known as [**feature extraction** or **dimension reduction**](https://en.wikipedia.org/wiki/Feature_extraction) in machine learning, the goal of variable reduction is to **reduce the number of predictors** by deriving  a new set of variables intended to be informative and non-redundant from a set of measured data. This method can be used to **simplify models**, which can benefit model interpretation, shorten fitting time, and improve generalization (by reducing overfitting).

## Quick and Exploratory Method

Let's start by fitting a multiple linear regression model with the `attitude` dataset, available is base R, to predict the overall **rating** by employees of their organization with the remaining variables (handling of employee **complaints**, special **privileges**, opportunity of **learning**, **raises**, a feedback considered too **critical** and opportunity of **advancement**).

```{r}
data("attitude")
model <- lm(rating ~ ., data = attitude)
parameters(model)
```

We can explore a reduction of the number of parameters with the `reduce_parameters()` function.

```{r}
newmodel <- reduce_parameters(model)
parameters(newmodel)
```

This output *hints* at the fact that the model could be represented via **two "latent" dimensions**, one correlated with all the positive things that a company has to offer, and the other one related to the amount of negative critiques received by the employees. These two dimensions have a positive and negative relationship with the company rating, respectively.

> What does `reduce_parameters()` exactly do?

This function performs a reduction in the parameter space (the number of variables). It starts by creating a new set of variables, based on the chosen method (the default method is "**PCA**", but other are available via the `method` argument, such as "**cMDS**", "**DRR**" or "**ICA**"). Then, it names this new dimensions using the original variables that *correlate* the most with it. For instance, in the example above a variable named `raises_0.88/learning_0.82/complaints_0.78/privileges_0.70/advance_0.68` means that the respective variables (`raises`, `learning`, `complaints`, `privileges`, `advance`) correlate maximally (with coefficients of .88, .82, .78, .70, .68, respectively) with this dimension.

```{r}
reduce_parameters(model, method = "cMDS") %>%
  parameters()
```

A different method (**Classical Multidimensional Scaling - cMDS**) suggests that negative critiques do not have a significant impact on the rating, and  that the lack of opportunities of career advancement is a separate dimension with an importance on its own.

Although `reduce_parameters()` function can be useful in exploratory data analysis, it's best to perform the dimension reduction step in a **separate and dedicated stage**, as this is a very important process in the data analysis workflow.

## Principal Component Analysis (PCA)

PCA is a widely used procedure that lies in-between dimension reduction and structural modeling. Indeed, one of the ways of reducing the number of predictors is to extract a new set of uncorrelated variables that will *represent* variance of your initial dataset. But how the original variables relate between themselves can also be a question on its own.

We can apply the `principal_components()` function to do the the predictors of the model:

```{r}
pca <- principal_components(insight::get_predictors(model), n = "auto")
pca
```

The `principal_components()` function automatically selected one component (if the number of components is not specified, this function uses [`n_factors()`](https://easystats.github.io/parameters/reference/n_factors.html) to estimate the optimal number to keep) and returned the **loadings**, i.e., the relationship with all of the original variables.

As we can see here, it seems that our new component captured the essence (more than half of the total variance present in the original dataset) of all our other variables together. We can **extract** the values of this component for each of our observation using the `predict()` method and add in the response variable of our initial dataset.

```{r}
newdata <- predict(pca)
newdata$rating <- attitude$rating
```

We can know update the model with this new component:

```{r}
update(model, rating ~ Component_1, data = newdata) %>%
  parameters()
```

### Using the `psych` package for PCA

You can also use different packages for models, such as [`psych`](https://cran.r-project.org/package=psych) [@revelle2018] or [`FactoMineR`](http://factominer.free.fr/) for PCA or Exploratory Factor Analysis (EFA), as it allows for more flexibility and control when running such procedures. 

The functions from this package are **fully supported** by `parameters` through the `model_parameters()` function. For instance, we can redo the above analysis using the `psych` package as follows:

```{r}
library(psych)

# Fit the PCA
pca <- model_parameters(psych::principal(attitude, nfactors = 1))
pca
```

*Note:* By default, `psych::principal()` uses a **varimax** rotation to extract rotated components, possibly leading to discrepancies in the results.

Finally, refit the model:

```{r eval=FALSE}
df <- cbind(attitude, predict(pca))

update(model, rating ~ PC1, data = df) %>%
  model_parameters()
```

# References
---
title: "Selection of Model Parameters"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, variable selection, feature selection]
vignette: >
  %\VignetteIndexEntry{Selection of Model Parameters}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r , include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%"
)

if (!requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("performance", quietly = TRUE) ||
  !requireNamespace("rstanarm", quietly = TRUE) ||
  # !requireNamespace("projpred", quietly = TRUE) ||
  !requireNamespace("lme4", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(poorman)
  library(rstanarm)
  library(lme4)
}

set.seed(333)
```

Also known as [**feature selection**](https://en.wikipedia.org/wiki/Feature_selection) in machine
learning, the goal of variable selection is to **identify a subset of predictors** to **simplify models**. This can benefit model interpretation,
shorten fitting time, and improve generalization (by reducing overfitting).

There are many different methods. The appropriate method for a given problem
will depend on the model type, the data, the objective, and the theoretical
rationale.

The `parameters` package implements a helper that will **automatically** pick a
method deemed appropriate for the provided model, run the variables selection
and return the **optimal formula**, which you can then re-use to update the
model.

## Simple linear regression

### Fit a powerful model

If you are familiar with R and the formula interface, you know of the
possibility of including a dot (`.`) in the formula, signifying "all the
remaining variables". Curiously, few are aware of the possibility of
additionally easily adding "all the interaction terms". This can be achieved
using the `.*.` notation.

Let's try that with the linear regression predicting **Sepal.Length** with the
[`iris`](https://en.wikipedia.org/wiki/Iris_flower_data_set) dataset, included
by default in R.

```{r}
model <- lm(Sepal.Length ~ . * ., data = iris)
summary(model)
```

***Wow, that's a lot of parameters! And almost none of them are significant!***

Which is ***weird***, considering that **gorgeous $R^2$ of 0.882!** 

*I wish I had that in my research!*

### Too many parameters?

As you might know, having a **model that is too performant is not always a good
thing**. For instance, it can be a marker of
[**overfitting**](https://en.wikipedia.org/wiki/Overfitting): the model
corresponds too closely to a particular set of data, and may therefore fail to
predict future observations reliably. In multiple regressions, in can also fall
under the [**Freedman's paradox**](https://en.wikipedia.org/wiki/Freedman%27s_paradox): some predictors
that have actually no relation to the dependent variable being predicted will be
**spuriously found to be statistically significant**.

Let's run a few checks using the
[**performance**](https://github.com/easystats/performance) package:

```{r}
library(performance)

check_normality(model)
check_heteroscedasticity(model)
check_autocorrelation(model)
check_collinearity(model)
```

The main issue of the model seems to be the high
[multicollinearity](https://en.wikipedia.org/wiki/Multicollinearity). This
suggests that our model might not be able to give valid results about any
individual predictor, nor tell which predictors are redundant with respect to
others.

### Parameters selection

Time to do some variables selection! This can be easily done using the
`select_parameters()` function in `parameters`. It will **automatically** select
the best variables and update the model accordingly. One way of using that is in
a tidy pipeline (using [`%>%`](https://cran.r-project.org/package=magrittr/)),
using this output to update a new model.

```{r}
lm(Sepal.Length ~ . * ., data = iris) %>%
  select_parameters() %>%
  summary()
```

That's still a lot of parameters, but as you can see, almost all of them are
now significant, and the $R^2$ did not change much.

Although appealing, please note that these automated selection methods are
[**quite criticized**](https://towardsdatascience.com/stopping-stepwise-why-stepwise-selection-is-bad-and-what-you-should-use-instead-90818b3f52df),
and should not be used in place of **theoretical** or **hypothetical** reasons
(*i.e.*, you should have *a priori* hypotheses about which parameters of your
model you want to focus on).

## Mixed and Bayesian models

For simple linear regressions as above, the selection is made using the `step()`
function (available in base R). This performs a
[**stepwise**](https://en.wikipedia.org/wiki/Stepwise_regression) selection.
However, this procedures is not available for other types of models, such as
**mixed** or **Bayesian** models.

### Mixed models

For mixed models (of class `merMod`), stepwise selection is based on
`cAIC4::stepcAIC()`. This step function only searches the "best" model based on
the _random effects structure_, i.e. `select_parameters()` adds or excludes
random effects until the `cAIC` can't be improved further.

This is what our initial model looks like.

```{r}
library(lme4)
data("qol_cancer")

# initial model
lmer(
  QoL ~ time + phq4 + age + (1 + time | hospital / ID),
  data = qol_cancer
) %>%
  summary()
```

This is the model selected by `select_parameters()`. Please notice the
differences in the random effects structure between the initial and the selected
models:

```{r}
# multiple models are checked, however, initial models
# already seems to be the best one...
lmer(
  QoL ~ time + phq4 + age + (1 + time | hospital / ID),
  data = qol_cancer
) %>%
  select_parameters() %>%
  summary()
```
---
title: "Model Parameters for Multiply Imputed Repeated Analyses"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, variable selection, feature selection]
vignette: >
  %\VignetteIndexEntry{Model Parameters for Multiply Imputed Repeated Analyses}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)
knitr::opts_chunk$set(comment = "#>")

if (!requireNamespace("mice", quietly = TRUE) ||
  !requireNamespace("GLMMadaptive", quietly = TRUE) ||
  !requireNamespace("lme4", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
}

set.seed(333)
```

## Model Parameters from `mira` objects

`model_parameters()` can be used in combination with the *mice* package to deal
with missing data, in particular to summaries regression models used with
multiple imputed datasets. It computes pooled summaries of multiple imputed
repeated regression analyses, i.e. of objects of class `mira`. Thus,
`model_parameters()` for `mira`-objects is comparable to the `pool()`-function
from *mice*, but only focuses on the final summary of parameters and does not
include the diagnostic statistic per estimate.

```{r message=FALSE, warning=FALSE}
library(mice)
library(parameters)

data("nhanes2")
imp <- mice(nhanes2, printFlag = FALSE)
fit <- with(data = imp, exp = lm(bmi ~ age + hyp + chl))

model_parameters(fit)
```

Not all packages work with `with.mids()` from package *mice*. Thus, for some
modeling packages, it's not possible to perform multiply imputed repeated
analyses, i.e. you cannot work with imputed data for such models. We give an
example for the *GLMMadaptive* package here.

First, we generate a dataset with missing values. We take the data `cbpp` from
*lme4* and randomly assign some missing values into one of the predictors. Then
we impute the data, using `mice()` from package *mice*.

```{r message=FALSE, warning=FALSE}
library(lme4)
library(GLMMadaptive)

data(cbpp)
cbpp$period[sample(1:nrow(cbpp), size = 10)] <- NA

imputed_data <- mice(cbpp, printFlag = FALSE)
```

Using `with` to compute multiple regression analyses for each imputed dataset
fails.

```{r message=FALSE, eval=FALSE}
fit <- with(data = imputed_data, expr = GLMMadaptive::mixed_model(
  cbind(incidence, size - incidence) ~ period,
  random = ~ 1 | herd,
  family = binomial
))
# > Error in as.data.frame(data) :
# >   argument "data" is missing, with no default
```

However, we can use a workaround by using `pool_parameters()`, which works on a
list of model objects. So whenever a model-object is not yet supported by
`mice::with()`, you can instead fit multiple models to the imputed datasets and
pool all parameters with `pool_parameters()`:

The steps would be:

1. Calculate the regression models for each imputed dataset manually (either by
   using `complete()` from package *mice* to get the imputed datasets, or by
   accessing the datasets directly from the `mids` object)

2. Save all model objects in a list.

3. Pass the list to `pool_parameters()`.

```{r message=FALSE}
models <- lapply(1:imputed_data$m, function(i) {
  mixed_model(
    cbind(incidence, size - incidence) ~ period,
    random = ~ 1 | herd,
    data = complete(imputed_data, action = i),
    family = binomial
  )
})
pool_parameters(models)
```

For comparison and to show that the results from `mice:pool()` and
`pool_parameters()` are identical, we take an example that also works with the
_mice_ package:

```{r message=FALSE}
library(mice)
library(parameters)

data("nhanes2")
imp <- mice(nhanes2, printFlag = FALSE)

# approach when model is supported by "mice"
fit <- with(data = imp, exp = lm(bmi ~ age + hyp + chl))
summary(pool(fit))

# approach when model is *not* supported by "mice"
models <- lapply(1:5, function(i) {
  lm(bmi ~ age + hyp + chl, data = complete(imp, action = i))
})
pool_parameters(models)
```

## Model Parameters from `mipo` objects

It is also possible to compute summaries of pooled objects of class `mipo`.

```{r message=FALSE, warning=FALSE}
data("nhanes2")
imp <- mice(nhanes2, printFlag = FALSE)
fit <- with(data = imp, exp = lm(bmi ~ age + hyp + chl))
pooled <- pool(fit)

model_parameters(pooled)
```
---
title: "Structural Models (EFA, CFA, SEM...)"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, efa, cfa, factor analysis, sem, fa, pca, how many factors, n factors]
vignette: >
  %\VignetteIndexEntry{Structural Models (EFA, CFA, SEM, ...)}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r , include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%"
)

if (!requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("see", quietly = TRUE) ||
  !requireNamespace("lavaan", quietly = TRUE) ||
  !requireNamespace("performance", quietly = TRUE) ||
  !requireNamespace("nFactors", quietly = TRUE) ||
  !requireNamespace("datawizard", quietly = TRUE) ||
  !requireNamespace("GPArotation", quietly = TRUE) ||
  !requireNamespace("psych", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(performance)
  library(GPArotation)
  library(psych)
  library(nFactors)
  library(poorman)
  library(lavaan)
}

set.seed(333)
```

# How to perform a Factor Analysis (FA)

The difference between PCA and EFA can be quite hard to intuitively grasp as
their output is very familiar. The idea is that PCA aims at extracting the most
variance possible from all variables of the dataset, whereas EFA aims at
creating consistent factors from the dataset without desperately trying to
represent all the variables.

This is why PCA is popular for feature reduction, as it will try to best
represent the variance contained in the original data, minimizing the loss of
information. On the other hand, EFA is usually in the context of exploring the
latent dimensions that might be hidden in the observed variables, without
necessarily striving to represent the whole dataset.

To illustrate EFA, let us use the [International Personality Item Pool](https://ipip.ori.org/) data available in the
[`psych`](https://www.personality-project.org/r/html/bfi.html) package. It
includes 25 personality self report items. The authors built these items
following the **big 5** personality structure.

## Factor Structure (Sphericity and KMO)

The first step is to test if the dataset is suitable for carrying out factor analysis. There are two

- **Bartlett's Test of Sphericity**: This tests whether a matrix (of correlations)
is significantly different from an identity matrix. The test provides
probability that the correlation matrix has significant correlations among at
least some of the variables in a dataset, a prerequisite for factor analysis to
work. In other words, before starting with factor analysis, one needs to check
whether Bartlett’s test of sphericity is significant.

- **Kaiser, Meyer, Olkin (KMO) Measure of Sampling Adequacy (MSA)**: This test
was introduced by Kaiser (1970) as the Measure of Sampling Adequacy (MSA), later
modified by Kaiser and Rice (1974). The Kaiser-Meyer-Olkin (KMO) statistic,
which can vary from 0 to 1, indicates the degree to which each variable in a set
is predicted without error by the other variables. A value of 0 indicates that
the sum of partial correlations is large relative to the sum correlations,
indicating factor analysis is likely to be inappropriate. A KMO value close to 1
indicates that the sum of partial correlations is not large relative to the sum
of correlations and so factor analysis should yield distinct and reliable
factors.

Both tests can be performed by using the `check_factorstructure()` function.

```{r}
library(parameters)
library(psych)

# Load the data
data <- psych::bfi[, 1:25] # Select only the 25 first columns corresponding to the items
data <- na.omit(data) # remove missing values

# Check factor structure
check_factorstructure(data)
```

## Exploratory Factor Analysis (EFA)

Now that we are confident that our dataset is appropriate, we will explore a
factor structure made of 5 latent variables, corresponding to the items' authors
theory of personality.

```{r}
# Fit an EFA
efa <- psych::fa(data, nfactors = 5) %>%
  model_parameters(sort = TRUE, threshold = "max")

efa
```

As we can see, the 25 items nicely spread on the 5 latent factors, the famous
**big 5**. Based on this model, we can now predict back the scores for each
individual for these new variables:

```{r}
# let's look only at the first five individuals
head(predict(efa, names = c("Neuroticism", "Conscientiousness", "Extraversion", "Agreeableness", "Opennness")), 5)
```

## How many factors to retain in Factor Analysis (FA)

When running a **factor analysis (FA)**, one often needs to specify **how many components** (or latent variables) to retain or to extract. This decision is
often motivated or supported by some statistical indices and procedures aiming
at finding the optimal number of factors.

There are a huge number of methods exist to statistically address this issue,
and they can sometimes give very different results.

> **Unfortunately, there is no consensus on which method to use, or which is the
best.**

### The Method Agreement procedure

The Method Agreement procedure, first implemented in the
[`psycho`](https://neuropsychology.github.io/psycho.R/2018/05/24/n_factors.html)
package [@makowski2018psycho], proposes to rely on the consensus of methods,
rather than on one method in particular.

This procedure can be easily used via the `n_factors()` function, re-implemented
and improved in the [**parameters**](https://github.com/easystats/parameters)
package. One can provide a dataframe, and the function will run a large number
of routines and return the optimal number of factors based on the higher
consensus.

```{r}
n <- n_factors(data)
n
```

Interestingly, the smallest nubmer of factors that most methods suggest is 6,
which is consistent with the newer models of personality (e.g., HEXACO).

More details, as well as a summary table can be obtained as follows:

```{r}
as.data.frame(n)
summary(n)
```

A plot can also be obtained (the `see` package must be loaded):

```{r}
library(see)

plot(n) + theme_modern()
```

## Confirmatory Factor Analysis (CFA)

We've seen above that while an EFA with 5 latent variables works great on our
dataset, a structure with 6 latent factors might in fact be more appropriate.
How can we **statistically test** if that is actually the case? This can be done
using **Confirmatory Factor Analysis (CFA)** (as opposed to **Exploratory** FA),
which bridges factor analysis with Structural Equation Modelling (SEM).

However, in order to do that cleanly, EFA should be **independent** from CFA:
the factor structure should be explored in a **"training" set**, and then tested
(or "confirmed") in a **"testing" set**. 

In other words, the dataset used for exploration and confirmation should not be
the same, a standard widely adopted in the field of machine learning.

### Partition the data

The data can be easily split into two sets with the `data_partition()` function,
through which we will use 70\% of the sample as training and the rest as test.

```{r}
# to have reproducible result, we will also set seed here so that similar
# portions of the data are used each time we run the following code
partitions <- datawizard::data_partition(data, training_proportion = 0.7, seed = 111)
training <- partitions$training
test <- partitions$test
```

### Create CFA structures out of EFA models

In the next step, we will run two EFA models on the training set, specifying 5
and 6 latent factors respectively, that we will then transform into CFA
structures.

```{r}
structure_big5 <- psych::fa(training, nfactors = 5) %>%
  efa_to_cfa()
structure_big6 <- psych::fa(training, nfactors = 6) %>%
  efa_to_cfa()

# Investigate how the models look
structure_big5

structure_big6
```

As we can see, a structure is just a string encoding how the **manifest variables** (the observed variables) are integrated into **latent variables**.

### Fit and Compare models

We can finally apply this structure to the testing dataset using the `lavaan`
package, and compare these models against each other:

```{r}
library(lavaan)
library(performance)

big5 <- lavaan::cfa(structure_big5, data = test)
big6 <- lavaan::cfa(structure_big6, data = test)

performance::compare_performance(big5, big6)
```


```{r eval=FALSE, include=FALSE}
performance::test_likelihoodratio(big5, big6)  # TODO: This doesn't work
```

All in all, it seems that the Big-5 structure remains quite reliable.

# Structural Equation Modeling

<!-- TO DO: needs a better intro -->

The previous example shows one of the enormous amount of modeling possibilities for
structural equation models, in particular an example for mediation analysis,
i.e. a model that estimates indirect effects in partial mediation structures.

```{r}
set.seed(1234)
X <- rnorm(100)
M <- 0.5 * X + rnorm(100)
Y <- 0.7 * M + rnorm(100)
df <- data.frame(X = X, Y = Y, M = M)

model <- " # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         "
fit <- lavaan::sem(model, data = df, test = "Satorra-Bentler")
model_parameters(fit)
```

# References
---
title: "Standardized Model Parameters"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, variable selection, feature selection]
vignette: >
  %\VignetteIndexEntry{Standardized Model Parameters}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%"
)

if (!requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("lme4", quietly = TRUE) ||
  !requireNamespace("effectsize", quietly = TRUE) ||
  !requireNamespace("lm.beta", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(poorman)
  library(effectsize)
  library(lm.beta)
}

set.seed(333)
```

The
[`model_parameters()`](https://easystats.github.io/parameters/articles/model_parameters.html)
function (also accessible via the shortcut `parameters()`) can also be used to
calculate standardized model parameters via the `standardize`-argument. Recall
that standardizing data/variable (*z*-scoring), i.e. centering and scaling,
involves expressing data in terms of standard deviation (i.e., mean = 0, SD =
1). That is, it the process of subtracting the mean and dividing the quantity by
standard deviation. Standardization can help avoid multicollinearity issues when more complex (polynomial, for instance) terms are included in the model.

There are different methods of standardizing model parameters (see also [`?effectsize::standardize_parameters`](https://easystats.github.io/effectsize/reference/standardize_parameters.html)):

- `"refit"`,
- `"posthoc"`
- `"smart"` 
- `"basic"`

If you are interested in more statistical and technical details, and how standardization methods relate to different (standardized) effect size measures, read the following vignette from *effectsize* package, from whence this functionality comes:
<https://easystats.github.io/effectsize/articles/standardize_parameters.html>

## Standardization by re-fitting the model

`standardize = "refit"` is based on a complete model re-fit with a standardized
version of data. Hence, this method is equal to standardizing the variables
*before* fitting the model. It is the most accurate (Neter et al., 1989), but it
is also the most computationally costly and long (especially for heavy models
such as, for instance, Bayesian models). This method is particularly recommended
for complex models that include interactions or transformations (e.g.,
polynomial or spline terms).

When `standardize = "refit"`, `model_parameters()` internally calls
[`effectsize::standardize()`](https://easystats.github.io/effectsize/reference/standardize.html)
to standardize the data that was used to fit the model and updates the model
with the standardized data. Note that `effectsize::standardize()` tries to
detect which variables should be standardized and which not. For instance,
having a `log(x)` in the model formula would exclude `x` from being
standardized, because `x` might get negative values, and thus `log(x)` would no
longer be defined. Factors or dates will also *not* be standardized. Response
variables will be standardized, if appropriate.

```{r}
library(lme4)
data(iris)
set.seed(1234)
iris$grp <- as.factor(sample(1:3, nrow(iris), replace = TRUE))

# fit example model
model <- lme4::lmer(
  Sepal.Length ~ Species * Sepal.Width + Petal.Length + (1 | grp),
  data = iris
)

# classic model parameters
model_parameters(model)

# standardized model parameters
model_parameters(model, standardize = "refit")
```

The second output is identical to following:

```{r}
# standardize continuous variables manually
model2 <- lme4::lmer(
  scale(Sepal.Length) ~ Species * scale(Sepal.Width) + scale(Petal.Length) + (1 | grp),
  data = iris
)

model_parameters(model2)
```

## Post-hoc standardization

`standardize = "posthoc"` aims at emulating the results obtained by `"refit"`
without refitting the model. The coefficients are divided by the standard
deviation of the outcome (which becomes their expression *unit*). Then, the
coefficients related to numeric variables are additionally multiplied by the
standard deviation of the related terms, so that they correspond to changes of 1
SD of the predictor (e.g., "a change in 1 SD of `x` is related to a change of
0.24 of the SD of `y`"). This does not apply to binary variables or factors, so
the coefficients are still related to changes in levels.

This method is not accurate and tends to give aberrant results when interactions
are specified. However, this method of standardization is the "classic" result
obtained by many statistical packages when standardized coefficients are
requested.

When `standardize = "posthoc"`, `model_parameters()` internally calls
[`effectsize::standardize_parameters(method = "posthoc")`](https://easystats.github.io/effectsize/reference/standardize_parameters.html).
Test statistic and p-values are not affected, i.e. they are the same as if no
standardization would be applied.

```{r}
model_parameters(model, standardize = "posthoc")
```

`standardize = "basic"` also applies post-hoc standardization, however, factors
are converted to numeric, which means that it also scales the coefficient by the
standard deviation of model's matrix' parameter of factor levels (transformed to
integers) or binary predictors.

```{r}
model_parameters(model, standardize = "basic")
```

Compare the two outputs above and notice how coefficient estimates, standard
errors, confidence intervals, and *p*-values change for main effect and
interaction effect terms containing `Species` variable, the only factor variable
in our model.

This method is the one implemented by default in other software packages, such as `lm.beta::lm.beta()`:

```{r}
library(lm.beta)
data(iris)
model3 <- lm(Sepal.Length ~ Species * Sepal.Width + Petal.Length, data = iris)
mp <- model_parameters(model3, standardize = "basic")
out <- lm.beta(model3)

data.frame(model_parameters = mp$Std_Coefficient, lm.beta = coef(out))
```


## Smart standardization

`standardize = "smart"` is similar to `standardize = "posthoc"` in that it does
not involve model re-fitting. The difference is that the SD of the response is
computed on the relevant section of the data. For instance, if a factor with 3
levels A (the intercept), B and C is entered as a predictor, the effect
corresponding to B versus A will be scaled by the variance of the response at
the intercept only. As a results, the coefficients for effects of factors are
similar to a Glass' delta.

```{r}
model_parameters(model, standardize = "smart")
```
---
title: "Printing Model Parameters"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, table layout]
vignette: >
  %\VignetteIndexEntry{Printing Model Parameters}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r , include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%",
  tidy.opts = list(width.cutoff = 100)
)

if (!requireNamespace("gt", quietly = TRUE) ||
  !requireNamespace("magrittr", quietly = TRUE) ||
  !requireNamespace("glmmTMB", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(gt)
  library(magrittr)
  library(insight)
}

set.seed(333)
```

`model_parameters()` and `compare_parameters()` are functions that return a data frame of model summaries in a consistent way. The printed table of those summaries is formatted to make the output more readable and removes or collapses redundant columns, to get a compact and yet comprehensive summary table. _(N.B. for developers: the function [`standardize_names()`](https://easystats.github.io/insight/reference/standardize_names.html) standardizes the column names, so column names are consistent and the same for any model object, also in **broom** style, which makes it easy to build your packages on top of the **parameters** package.)_

The default [print-methods](https://easystats.github.io/parameters/reference/print.parameters_model.html) for `model_parameters()` and `compare_parameters()` allows the user to modify the layout and style of the output.

# Summaries for a single model

In the following examples for `model_parameters()`, which returns tabular output for single models, are shown.

## Pretty parameter names formatting

By default, the argument `pretty_names` is `TRUE`, meaning that parameter names are formatted to make them more "human readable", i.e. factor levels are separated from the variable names, interactions are denoted by `*` etc.

```{r}
library(parameters)
data(iris)
model <- lm(Sepal.Length ~ Species * Petal.Length, data = iris)
model_parameters(model)

mp <- model_parameters(model)
print(mp, pretty_names = FALSE)
```

## Splitting model components

Again by default, the argument `split_components` is `TRUE`, which means that models with multiple components like fixed and random effects, count and zero-inflated part etc. are split into separate tables in the output.

```{r}
library(glmmTMB)
data("Salamanders")
model <- glmmTMB(count ~ spp + mined + (1 | site),
                 ziformula = ~spp + mined,
                 family = nbinom2(), 
                 data = Salamanders)
model_parameters(model)
```

Redundant columns are removed. The related model component is shown as table header. However, you can also return a single table:

```{r}
mp <- model_parameters(model)
print(mp, split_component = FALSE)
```

## Adding model summaries

A model summary can be added to the table when `summary = TRUE` in the call to `model_parameters()`:

```{r}
model <- lm(Sepal.Length ~ Species * Petal.Length, data = iris)
model_parameters(model, summary = TRUE)
```

## Changing number of digits

`digits` changes the digits for coefficients, standard errors and statistics. `ci_digits` and `p_digits` are especially for the confidence intervals and p-values.

```{r}
model <- lm(Sepal.Length ~ Species, data = iris)
model_parameters(model, digits = 4)
```

p-values can be displayed in exact, scientific notation if required.

```{r}
model_parameters(model, p_digits = "scientific")
```

## Fixing column widths

By default, the width of table columns is set to the minimum required width. This works well for models that produce just one table. However, for models with multiple components, where each component is shown as separate table, columns are possibly no longer aligned across tables. See the following example from a zero-inflated mixed model that has three components (fixed count, fixed zero-inflated, random effects):

```{r}
data("Salamanders")
# we create very long parameter names for this predictor here
levels(Salamanders$spp) <- paste("long", levels(Salamanders$spp))

model <- glmmTMB(
 count ~ spp + mined + (1 | site),
 ziformula = ~mined,
 family = poisson(),
 data = Salamanders
)

# default printing
model_parameters(model) 
```

The `column_width` argument can be used to either define the width of specific columns, or to fix column widths of the same columns across tables to have the same width. In the latter case, use `column_width = "fixed"` in the `print()` method.

```{r}
mp <- model_parameters(model)
print(mp, column_width = "fixed")
```

If `column_width` is a named vector, names are matched against column names, and those columns gain the specified minimum width.

```{r}
print(mp, column_width = c(SE = 8, `95% CI` = 12, p = 7))
```

## Group parameters

The `groups` argument can be used to group parameters in the table. `groups` must be a named list, where the names of the list elements equal the header of each group, while the values of the list elements equal the parameter names, or the position of the parameters in the table (data frame).

In the following example, we see the names of the parameters in the `Parameter` column, while the rownumbers indicate their position.

```{r}
data(mtcars)
mtcars$cyl <- as.factor(mtcars$cyl)
mtcars$gear <- as.factor(mtcars$gear)
model <- lm(mpg ~ hp + gear * vs + cyl + drat, data = mtcars)

# don't select "Intercept" parameter
mp <- model_parameters(model, drop = "^\\(Intercept")

# inspect data frame
as.data.frame(mp)
```

Now we create a group named `"Engine"`, which encompasses the parameters `"cyl6"`, `"cyl8"`, `"vs"` and `"hp"`. The `"Interactions"` group includes `"gear4:vs"` and `"gear5:vs"`. The group `"controls"` has the parameters from rows 2, 3 and 7.

Note that the parameters in the table summary are re-ordered according to the order specified in `groups`.

```{r}
# group parameters, either by parameter name or position
print(mp, groups = list("Engine" = c("cyl6", "cyl8", "vs", "hp"),
                        "Interactions" = c("gear4:vs", "gear5:vs"), 
                        "Controls" = c(2, 3, 7))) # gear 4 and 5, drat
```

If you prefer tables without vertical borders, use the `sep` argument to define the string that is used as border-separator. This argument is passed down to `insight::export_table()`.

```{r}
# group parameters, either by parameter name or position
print(mp, sep = "  ",
      groups = list("Engine" = c("cyl6", "cyl8", "vs", "hp"),
                    "Interactions" = c("gear4:vs", "gear5:vs"), 
                    "Controls" = c(2, 3, 7)))
```

# Summaries for multiple models

`compare_parameters()` (or its alias `compare_models()`) allows to create tables for multiple models, aligned side by side.

By default, estimates and confidence intervals are shown.

```{r}
data(iris)
lm1 <- lm(Sepal.Length ~ Species, data = iris)
lm2 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm3 <- lm(Sepal.Length ~ Species * Petal.Length, data = iris)
compare_parameters(lm1, lm2, lm3)
```

## Changing style of column output

By default, estimates and confidence intervals are shown. Using `style` allows us to create different output, e.g. standard errors instead of confidence intervals, or including p-values.

```{r}
compare_parameters(lm1, lm2, lm3, style = "se_p")
```

## Defining column names

The column names for the models are by default the objects' names. You can define own names using the `column_names` argument.

```{r}
compare_parameters(
  lm1, lm2, lm3,
  column_names = c("First Model", "Second Model", "Third Model")
)
```

## Group parameters of multiple model tables

Grouping parameters works for `compare_models()` in the same way as shown above for `model_parameters()`.

```{r}
lm1 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm2 <- lm(Sepal.Width ~ Species * Petal.Length, data = iris)

# remove intercept
cp <- compare_parameters(lm1, lm2, drop = "^\\(Intercept")

# look at parameters names, to know their names for "groups" argument
as.data.frame(cp)$Parameter

# create groups. Interactions only present in 2nd model
print(cp, groups = list(Species = c("Species (versicolor)", 
                                    "Species (virginica)"),
                        Interactions = c("Species (versicolor) * Petal Length",
                                         "Species (virginica) * Petal Length"),
                        Controls = "Petal Length"))
```


# Splitting wide tables into multiple table parts

For very wide tables that cannot be displayed properly, you can use the `table_width` argument in the `print()` method to split tables into multiple parts. `table_width` can be a numeric value, or `"auto"`, indicating the width of the complete table. If `table_width = "auto"` and the table is wider than the current available width (i.e. line length) of the console (or any other source for textual output, like markdown files), the table is split into multiple parts. Else, if `table_width` is numeric and table rows are wider than `table_width`, the table is split into multiple parts.

```{r}
data(iris)
lm1 <- lm(Sepal.Length ~ Species, data = iris)
lm2 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm3 <- lm(Sepal.Length ~ Species * Petal.Length, data = iris)
lm4 <- lm(Sepal.Length ~ Species * Petal.Length + Petal.Width, data = iris)

# very wide table
compare_parameters(lm1, lm2, lm3, lm4)

# table split into two parts
tab <- compare_parameters(lm1, lm2, lm3, lm4)
print(tab, table_width = 80)
```

# More advances tables and markdown / HTML formatting

The `print_md()` as well as `print_html()` functions can be used to create markdown (for knitting to PDF or Word) and HTML tables. 

Meanwhile, there are a lot of additional packages that allow users to have even more flexibility regarding table layouts. One package we can recommend is the [*modelsummary* package](https://vincentarelbundock.github.io/modelsummary/).
---
title: "Formatting Model Parameters"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, parameters, table layout]
vignette: >
  %\VignetteIndexEntry{Formatting Model Parameters}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r , include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
options(digits = 2)

knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  out.width = "100%",
  tidy.opts = list(width.cutoff = 100)
)

if (!requireNamespace("broom", quietly = TRUE) ||
  !requireNamespace("gt", quietly = TRUE) ||
  !requireNamespace("magrittr", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(parameters)
  library(broom)
  library(gt)
  library(magrittr)
  library(insight)
}

set.seed(333)
```

The *parameters* package, together with the [*insight*
package](https://easystats.github.io/insight/), provides tools to format the
layout and style of tables from model parameters. When you use the
`model_parameters()` function, you usually don't have to take care about
formatting and layout, at least not for simple purposes like printing to the
console or inside rmarkdown documents. However, sometime you may want to do the
formatting steps manually. This vignette introduces the various functions that
are used for parameters table formatting.

## An Example Model

We start with a model that does not make much sense, but it is useful for
demonstrating the formatting functions.

```{r}
data(iris)
iris$Petlen <- cut(iris$Petal.Length, breaks = c(0, 3, 7))
model <- lm(Sepal.Width ~ poly(Sepal.Length, 2) + Species + Petlen, data = iris)

summary(model)
```

## Formatting Parameter Names

As we can see, in such cases, the standard R output looks a bit cryptic,
although all necessary and important information is included in the summary. The
formatting of coefficients for polynomial transformation is difficult to read,
factors grouped with `cut()` always require a short time of thinking to find out
which of the bound (in this case, `Petlen(3,7]`, 3 and 7) is included in the
range, and names of factor levels are directly concatenated to the name of the
factor variable.

Thus, the first step would be to format the parameter names, which can be done
with `format_parameters()` from the *parameters* package:

```{r}
library(parameters)
format_parameters(model)
```

`format_parameters()` returns a (named) character vector with the original
coefficients as _names_ of each character element, and the formatted names of
the coefficients as values of the character vector. Let's look at the results
again:

```{r}
cat(format_parameters(model), sep = "\n")
```

Now variable names and factor levels, but also polynomial terms or even factors
grouped with `cut()` are much more readable. Factor levels are separated from
the variable name, inside brackets. Same for the coefficients of the different
polynomial degrees. And the exact range for `cut()`-factors is also clearer
now.

## Standardizing Column Names of Parameter Tables

As seen above, the `summary()` returns columns named `Estimate`, `t value` or
`Pr(>|t|)`. While `Estimate` is not specific for certain models, `t value` is.
For logistic regression models, you would get `z value`. Some packages alter the
names, so you get just `t` or `t-value` etc.

`model_parameters()` also uses context-specific column names, where applicable:

```{r}
colnames(model_parameters(model))
```

For Bayesian models, `Coefficient` is usually named `Median` etc. While this
makes sense from a user perspective, because you instantly know which type of
statistic or coefficient you have, it becomes difficult when you need a generic
naming scheme to access model parameters when the input model is unknown. This
is the typical approach from the *broom* package, where you get "standardized"
column names:

```{r}
library(broom)
colnames(tidy(model))
```

To deal with such situations, the *insight* package provides a
`standardize_names()` function, which exactly does that: standardizing the
column names of the input. In the following example, you see that the
statistic-column is no longer named `t`, but `statistic`. `df_error` or
`df_residuals` will be renamed to `df`.

```{r}
library(insight)
library(magrittr)
model %>%
  model_parameters() %>%
  standardize_names() %>%
  colnames()
```

Furthermore, you can request "broom"-style for column names:

```{r}
model %>%
  model_parameters() %>%
  standardize_names(style = "broom") %>%
  colnames()
```

## Formatting Column Names and Columns

Beside formatting parameter names (coefficient names) using
`format_parameters()`, we can do even more to make the output more readable.
Let's look at an example that includes confidence intervals.

```{r}
cbind(summary(model)$coefficients, confint(model))
```

We can get a similar tabular output using *broom*.

```{r}
tidy(model, conf.int = TRUE)
```

Some improvements according to readability could be collapsing and formatting
the confidence intervals, and maybe the p-values. This would require some
effort, for instance, to format the values of the lower and upper confidence
intervals and collapsing them into one column. However, the `format_table()`
function is a convenient function that does all the work for you.

`format_table()` requires a data frame with model parameters as input, however,
there are some requirements to make `format_table()` work. In particular, the
column names must follow a certain pattern to be recognized, and this pattern
may either be the naming convention from *broom* or the [*easystats*
packages](https://easystats.github.io/easystats/).

```{r}
model %>%
  tidy(conf.int = TRUE) %>%
  format_table()
```

When the parameters table also includes degrees of freedom, and the degrees of
freedom are the same for each parameter, then this information is included in
the statistic-column. This is usually the default for `model_parameters()`:

```{r}
model %>%
  model_parameters() %>%
  format_table()
```

## Exporting the Parameters Table

Finally, `export_table()` from *insight* formats the data frame and returns a
character vector that can be printed to the console or inside rmarkdown
documents. The data frame then looks more "table-like".

```{r}
data(mtcars)
export_table(mtcars[1:8, 1:5])
```

Putting all this together allows us to create nice tabular outputs of parameters
tables. This can be done using *broom*:

```{r}
model %>%
  tidy(conf.int = TRUE) %>%
  format_table() %>%
  export_table()
```

Or, in a simpler way and with much more options (like standardizing, robust
standard errors, bootstrapping, ...) using `model_parameters()`, which
`print()`-method does all these steps automatically:

```{r}
model_parameters(model)
```

## Formatting the Parameters Table in Markdown

`export_table()` provides a few options to generate tables in markdown-format.
This allows to easily render nice-looking tables inside markdown-documents.
First of all, use `format = "markdown"` to activate the markdown-formatting.
`caption` can be used to add a table caption. Furthermore, `align` allows to
choose an alignment for all table columns, or to specify the alignment for each
column individually.

The following table has six columns. Using `align = "lcccrr"` would left-align
the first column, center columns two to four, and right-align the last two
columns.

```{r}
model %>%
  tidy(conf.int = TRUE) %>%
  # parenthesis look better in markdown-tables, so we use "brackets" here
  format_table(ci_brackets = c("(", ")")) %>%
  export_table(format = "markdown", caption = "My Table", align = "lcccrr")
```

`print_md()` is a convenient wrapper around `format_table()` and
`export_table(format = "markdown")`, and allows to directly format the output of
functions like `model_parameters()`, `simulate_parameters()` or other
_parameters_ functions in markdown-format.

These tables are also nicely formatted when knitting markdown-documents into
Word or PDF. `print_md()` applies some default settings that have proven to work
well for markdown, PDF or Word tables.

```{r}
model_parameters(model) %>% print_md()
```

A similar option is `print_html()`, which is a convenient wrapper for
`format_table()` and `export_table(format = "html")`. Using HTML in markdown has
the advantage that it will be properly rendered when exporting to PDF.

```{r}
model_parameters(model) %>% print_html()
```

`print_md()` and `print_html()` are considered as main functions for users who
want to generate nicely rendered tables inside markdown-documents. A wrapper
around these both is `display()`, which either calls `print_md()` or
`print_html()`.

```{r}
model_parameters(model) %>% display(format = "html")
```
---
title: "Overview of Vignettes"
output: 
  rmarkdown::html_vignette:
vignette: >
  %\VignetteIndexEntry{Overview of Vignettes}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r message=FALSE, warning=FALSE, include=FALSE}
library(knitr)
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  eval = TRUE
)
```

All package vignettes are available at [https://easystats.github.io/parameters/](https://easystats.github.io/parameters/).

## Function Overview

* [Function Reference](https://easystats.github.io/parameters/reference/index.html)

## Description of Parameters

* [Summary of Model Parameters](https://easystats.github.io/parameters/articles/model_parameters.html)
* [Standardized Model Parameters](https://easystats.github.io/parameters/articles/model_parameters_standardized.html)
* [Robust Estimation of Standard Errors, Confidence Intervals, and p-values](https://easystats.github.io/parameters/articles/model_parameters_robust.html)
* [Model Parameters for Multiply Imputed Repeated Analyses](https://easystats.github.io/parameters/articles/model_parameters_mice.html)

## Formatting and Printing

* [Formatting Model Parameters](https://easystats.github.io/parameters/articles/model_parameters_formatting.html)
* [Printing Model Parameters](https://easystats.github.io/parameters/articles/model_parameters_print.html)

## Dimension Reduction and Clustering

* [Feature Reduction (PCA, cMDS, ICA...)](https://easystats.github.io/parameters/articles/parameters_reduction.html)
* [Structural Models (EFA, CFA, SEM...)](https://easystats.github.io/parameters/articles/efa_cfa.html)
* [Selection of Model Parameters](https://easystats.github.io/parameters/articles/parameters_selection.html)
* [Clustering with easystats](https://easystats.github.io/parameters/articles/clustering.html)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_p_adjust.R
\name{format_p_adjust}
\alias{format_p_adjust}
\title{Format the name of the p-value adjustment methods}
\usage{
format_p_adjust(method)
}
\arguments{
\item{method}{Name of the method.}
}
\value{
A string with the full surname(s) of the author(s), including year of publication, for the adjustment-method.
}
\description{
Format the name of the p-value adjustment methods.
}
\examples{
library(parameters)

format_p_adjust("holm")
format_p_adjust("bonferroni")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape_loadings.R
\name{reshape_loadings}
\alias{reshape_loadings}
\alias{reshape_loadings.parameters_efa}
\alias{reshape_loadings.data.frame}
\title{Reshape loadings between wide/long formats}
\usage{
reshape_loadings(x, ...)

\method{reshape_loadings}{parameters_efa}(x, threshold = NULL, ...)

\method{reshape_loadings}{data.frame}(x, threshold = NULL, loadings_columns = NULL, ...)
}
\arguments{
\item{x}{A data frame or a statistical model.}

\item{...}{Arguments passed to or from other methods.}

\item{threshold}{A value between 0 and 1 indicates which (absolute) values
from the loadings should be removed. An integer higher than 1 indicates the
n strongest loadings to retain. Can also be \code{"max"}, in which case it
will only display the maximum loading per variable (the most simple
structure).}

\item{loadings_columns}{Vector indicating the columns corresponding to loadings.}
}
\description{
Reshape loadings between wide/long formats.
}
\examples{
library(parameters)
library(psych)

pca <- model_parameters(psych::fa(attitude, nfactors = 3))
loadings <- reshape_loadings(pca)

loadings
reshape_loadings(loadings)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_dbscan.R, R/methods_hclust.R,
%   R/methods_kmeans.R, R/methods_mclust.R, R/methods_pam.R
\name{model_parameters.dbscan}
\alias{model_parameters.dbscan}
\alias{model_parameters.hclust}
\alias{model_parameters.pvclust}
\alias{model_parameters.kmeans}
\alias{model_parameters.hkmeans}
\alias{model_parameters.Mclust}
\alias{model_parameters.pam}
\title{Parameters from Cluster Models (k-means, ...)}
\usage{
\method{model_parameters}{dbscan}(model, data = NULL, clusters = NULL, ...)

\method{model_parameters}{hclust}(model, data = NULL, clusters = NULL, ...)

\method{model_parameters}{pvclust}(model, data = NULL, clusters = NULL, ci = 0.95, ...)

\method{model_parameters}{kmeans}(model, ...)

\method{model_parameters}{hkmeans}(model, ...)

\method{model_parameters}{Mclust}(model, data = NULL, clusters = NULL, ...)

\method{model_parameters}{pam}(model, data = NULL, clusters = NULL, ...)
}
\arguments{
\item{model}{Cluster model.}

\item{data}{A data.frame.}

\item{clusters}{A vector with clusters assignments (must be same length as rows in data).}

\item{...}{Arguments passed to or from other methods.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}
}
\description{
Format cluster models obtained for example by \code{\link[=kmeans]{kmeans()}}.
}
\examples{
\donttest{
# DBSCAN ---------------------------
if (require("dbscan", quietly = TRUE)) {
  model <- dbscan::dbscan(iris[1:4], eps = 1.45, minPts = 10)

  rez <- model_parameters(model, iris[1:4])
  rez

  # Get clusters
  predict(rez)

  # Clusters centers in long form
  attributes(rez)$means

  # Between and Total Sum of Squares
  attributes(rez)$Sum_Squares_Total
  attributes(rez)$Sum_Squares_Between

  # HDBSCAN
  model <- dbscan::hdbscan(iris[1:4], minPts = 10)
  model_parameters(model, iris[1:4])
}
}
#
# Hierarchical clustering (hclust) ---------------------------
data <- iris[1:4]
model <- hclust(dist(data))
clusters <- cutree(model, 3)

rez <- model_parameters(model, data, clusters)
rez

# Get clusters
predict(rez)

# Clusters centers in long form
attributes(rez)$means

# Between and Total Sum of Squares
attributes(rez)$Total_Sum_Squares
attributes(rez)$Between_Sum_Squares
\donttest{
#
# pvclust (finds "significant" clusters) ---------------------------
if (require("pvclust", quietly = TRUE)) {
  data <- iris[1:4]
  # NOTE: pvclust works on transposed data
  model <- pvclust::pvclust(datawizard::data_transpose(data),
    method.dist = "euclidean",
    nboot = 50,
    quiet = TRUE
  )

  rez <- model_parameters(model, data, ci = 0.90)
  rez

  # Get clusters
  predict(rez)

  # Clusters centers in long form
  attributes(rez)$means

  # Between and Total Sum of Squares
  attributes(rez)$Sum_Squares_Total
  attributes(rez)$Sum_Squares_Between
}
}
\dontrun{
#
# K-means -------------------------------
model <- kmeans(iris[1:4], centers = 3)
rez <- model_parameters(model)
rez

# Get clusters
predict(rez)

# Clusters centers in long form
attributes(rez)$means

# Between and Total Sum of Squares
attributes(rez)$Sum_Squares_Total
attributes(rez)$Sum_Squares_Between
}
\dontrun{
#
# Hierarchical K-means (factoextra::hkclust) ----------------------
if (require("factoextra", quietly = TRUE)) {
  data <- iris[1:4]
  model <- factoextra::hkmeans(data, k = 3)

  rez <- model_parameters(model)
  rez

  # Get clusters
  predict(rez)

  # Clusters centers in long form
  attributes(rez)$means

  # Between and Total Sum of Squares
  attributes(rez)$Sum_Squares_Total
  attributes(rez)$Sum_Squares_Between
}
}
if (require("mclust", quietly = TRUE)) {
  model <- mclust::Mclust(iris[1:4], verbose = FALSE)
  model_parameters(model)
}
\dontrun{
#
# K-Medoids (PAM and HPAM) ==============
if (require("cluster", quietly = TRUE)) {
  model <- cluster::pam(iris[1:4], k = 3)
  model_parameters(model)
}
if (require("fpc", quietly = TRUE)) {
  model <- fpc::pamk(iris[1:4], criterion = "ch")
  model_parameters(model)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_wrs2.R
\name{model_parameters.t1way}
\alias{model_parameters.t1way}
\title{Parameters from robust statistical objects in \code{WRS2}}
\usage{
\method{model_parameters}{t1way}(model, keep = NULL, verbose = TRUE, ...)
}
\arguments{
\item{model}{Object from \code{WRS2} package.}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from robust statistical objects in \code{WRS2}
}
\examples{
if (require("WRS2") && packageVersion("WRS2") >= "1.1.3") {
  model <- t1way(libido ~ dose, data = viagra)
  model_parameters(model)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_averaging.R, R/methods_betareg.R,
%   R/methods_glmx.R
\name{model_parameters.averaging}
\alias{model_parameters.averaging}
\alias{model_parameters.betareg}
\alias{model_parameters.glmx}
\title{Parameters from special models}
\usage{
\method{model_parameters}{averaging}(
  model,
  ci = 0.95,
  component = c("conditional", "full"),
  exponentiate = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{betareg}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  component = c("conditional", "precision", "all"),
  standardize = NULL,
  exponentiate = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{glmx}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  component = c("all", "conditional", "extra"),
  standardize = NULL,
  exponentiate = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{Model object.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{component}{Model component for which parameters should be shown. May be
one of \code{"conditional"}, \code{"precision"} (\pkg{betareg}),
\code{"scale"} (\pkg{ordinal}), \code{"extra"} (\pkg{glmx}),
\code{"marginal"} (\pkg{mfx}), \code{"conditional"} or \code{"full"} (for
\code{MuMIn::model.avg()}) or \code{"all"}.}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods. For instance, when
\code{bootstrap = TRUE}, arguments like \code{type} or \code{parallel} are
passed down to \code{bootstrap_model()}, and arguments like \code{ci_method}
are passed down to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}.}

\item{bootstrap}{Should estimates be based on bootstrapped model? If
\code{TRUE}, then arguments of \link[=model_parameters.stanreg]{Bayesian regressions} apply (see also
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}).}

\item{iterations}{The number of bootstrap replicates. This only apply in the
case of bootstrapped frequentist models.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from special regression models not listed under one of the previous categories yet.
}
\examples{
library(parameters)
if (require("brglm2", quietly = TRUE)) {
  data("stemcell")
  model <- bracl(
    research ~ as.numeric(religion) + gender,
    weights = frequency,
    data = stemcell,
    type = "ML"
  )
  model_parameters(model)
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to rename
columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_simulate_model.R, R/methods_glmmTMB.R
\name{simulate_model}
\alias{simulate_model}
\alias{simulate_model.glmmTMB}
\title{Simulated draws from model coefficients}
\usage{
simulate_model(model, iterations = 1000, ...)

\method{simulate_model}{glmmTMB}(
  model,
  iterations = 1000,
  component = c("all", "conditional", "zi", "zero_inflated", "dispersion"),
  verbose = FALSE,
  ...
)
}
\arguments{
\item{model}{Statistical model (no Bayesian models).}

\item{iterations}{The number of draws to simulate/bootstrap.}

\item{...}{Arguments passed to or from other methods.}

\item{component}{Should all parameters, parameters for the conditional model,
or for the zero-inflated part of the model be returned? Applies to models
with zero-inflated component. \code{component} may be one of \code{"conditional"},
\code{"zi"}, \code{"zero-inflated"}, \code{"dispersion"} or \code{"all"}
(default). May be abbreviated.}

\item{verbose}{Toggle warnings and messages.}
}
\value{
A data frame.
}
\description{
Simulate draws from a statistical model to return a data frame of estimates.
}
\details{
\subsection{Technical Details}{
\code{simulate_model()} is a computationally faster alternative
to \code{bootstrap_model()}. Simulated draws for coefficients are based
on a multivariate normal distribution (\code{MASS::mvrnorm()}) with mean
\code{mu = coef(model)} and variance \code{Sigma = vcov(model)}.
}
\subsection{Models with Zero-Inflation Component}{
For models from packages \pkg{glmmTMB}, \pkg{pscl}, \pkg{GLMMadaptive} and
\pkg{countreg}, the \code{component} argument can be used to specify
which parameters should be simulated. For all other models, parameters
from the conditional component (fixed effects) are simulated. This may
include smooth terms, but not random effects.
}
}
\examples{
library(parameters)
model <- lm(Sepal.Length ~ Species * Petal.Width + Petal.Length, data = iris)
head(simulate_model(model))
\donttest{
if (require("glmmTMB", quietly = TRUE)) {
  model <- glmmTMB(
    count ~ spp + mined + (1 | site),
    ziformula = ~mined,
    family = poisson(),
    data = Salamanders
  )
  head(simulate_model(model))
  head(simulate_model(model, component = "zero_inflated"))
}
}
}
\seealso{
\code{\link[=simulate_parameters]{simulate_parameters()}},
\code{\link[=bootstrap_model]{bootstrap_model()}},
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_mfx.R
\name{p_value.poissonmfx}
\alias{p_value.poissonmfx}
\alias{p_value.betaor}
\alias{p_value.betamfx}
\title{p-values for Marginal Effects Models}
\usage{
\method{p_value}{poissonmfx}(model, component = c("all", "conditional", "marginal"), ...)

\method{p_value}{betaor}(model, component = c("all", "conditional", "precision"), ...)

\method{p_value}{betamfx}(
  model,
  component = c("all", "conditional", "precision", "marginal"),
  ...
)
}
\arguments{
\item{model}{A statistical model.}

\item{component}{Should all parameters, parameters for the conditional model,
precision-component or marginal effects be returned? \code{component} may be one
of \code{"conditional"}, \code{"precision"}, \code{"marginal"} or \code{"all"} (default).}

\item{...}{Currently not used.}
}
\value{
A data frame with at least two columns: the parameter names and the
p-values. Depending on the model, may also include columns for model
components etc.
}
\description{
This function attempts to return, or compute, p-values of marginal effects
models from package \pkg{mfx}.
}
\examples{
if (require("mfx", quietly = TRUE)) {
  set.seed(12345)
  n <- 1000
  x <- rnorm(n)
  y <- rnegbin(n, mu = exp(1 + 0.5 * x), theta = 0.5)
  d <- data.frame(y, x)
  model <- poissonmfx(y ~ x, data = d)

  p_value(model)
  p_value(model, component = "marginal")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.recode_to_zero}
\alias{.recode_to_zero}
\title{Recode a variable so its lowest value is beginning with zero}
\usage{
.recode_to_zero(x)
}
\description{
Recode a variable so its lowest value is beginning with zero
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_clusterstructure.R
\name{check_clusterstructure}
\alias{check_clusterstructure}
\title{Check suitability of data for clustering}
\usage{
check_clusterstructure(x, standardize = TRUE, distance = "euclidean", ...)
}
\arguments{
\item{x}{A data frame.}

\item{standardize}{Standardize the dataframe before clustering (default).}

\item{distance}{Distance method used. Other methods than "euclidean"
(default) are exploratory in the context of clustering tendency. See
\code{\link[stats:dist]{stats::dist()}} for list of available methods.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
The H statistic (numeric)
}
\description{
This checks whether the data is appropriate for clustering using the Hopkins'
H statistic of given data. If the value of Hopkins statistic is close to 0
(below 0.5), then we can reject the null hypothesis and conclude that the
dataset is significantly clusterable. A value for H lower than 0.25 indicates
a clustering tendency at the \verb{90\%} confidence level. The visual assessment of
cluster tendency (VAT) approach (Bezdek and Hathaway, 2002) consists in
investigating the heatmap of the ordered dissimilarity matrix. Following
this, one can potentially detect the clustering tendency by counting the
number of square shaped blocks along the diagonal.
}
\examples{
\donttest{
library(parameters)
check_clusterstructure(iris[, 1:4])
plot(check_clusterstructure(iris[, 1:4]))
}
}
\references{
\itemize{
\item Lawson, R. G., & Jurs, P. C. (1990). New index for clustering
tendency and its application to chemical problems. Journal of chemical
information and computer sciences, 30(1), 36-41.
\item Bezdek, J. C., & Hathaway, R. J. (2002, May). VAT: A tool for visual
assessment of (cluster) tendency. In Proceedings of the 2002 International
Joint Conference on Neural Networks. IJCNN02 (3), 2225-2230. IEEE.
}
}
\seealso{
\code{\link[=check_kmo]{check_kmo()}}, \code{\link[=check_sphericity_bartlett]{check_sphericity_bartlett()}} and \code{\link[=check_factorstructure]{check_factorstructure()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_factors.R
\name{.n_factors_sescree}
\alias{.n_factors_sescree}
\title{Standard Error Scree and Coefficient of Determination Procedures}
\usage{
.n_factors_sescree(eigen_values = NULL, model = "factors")
}
\description{
Standard Error Scree and Coefficient of Determination Procedures
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reduce_parameters.R
\name{reduce_parameters}
\alias{reduce_parameters}
\alias{reduce_data}
\title{Dimensionality reduction (DR) / Features Reduction}
\usage{
reduce_parameters(x, method = "PCA", n = "max", distance = "euclidean", ...)

reduce_data(x, method = "PCA", n = "max", distance = "euclidean", ...)
}
\arguments{
\item{x}{A data frame or a statistical model.}

\item{method}{The feature reduction method. Can be one of 'PCA', 'cMDS',
'DRR', 'ICA' (see the Details section).}

\item{n}{Number of components to extract. If \code{n="all"}, then \code{n} is
set as the number of variables minus 1 (\code{ncol(x)-1}). If
\code{n="auto"} (default) or \code{n=NULL}, the number of components is
selected through \code{\link[=n_factors]{n_factors()}} resp. \code{\link[=n_components]{n_components()}}.
In \code{\link[=reduce_parameters]{reduce_parameters()}}, can also be \code{"max"}, in which case
it will select all the components that are maximally pseudo-loaded (i.e.,
correlated) by at least one variable.}

\item{distance}{The distance measure to be used. Only applies when
\code{method = "cMDS"}. This must be one of "euclidean", "maximum",
"manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring
can be given.}

\item{...}{Arguments passed to or from other methods.}
}
\description{
This function performs a reduction in the parameter space (the number of
variables). It starts by creating a new set of variables, based on the given
method (the default method is "PCA", but other are available via the
\code{method} argument, such as "cMDS", "DRR" or "ICA"). Then, it names this
new dimensions using the original variables that correlates the most with it.
For instance, a variable named 'V1_0.97/V4_-0.88' means that the V1 and the
V4 variables correlate maximally (with respective coefficients of .97 and
-.88) with this dimension. Although this function can be useful in
exploratory data analysis, it's best to perform the dimension reduction step
in a separate and dedicated stage, as this is a very important process in the
data analysis workflow. \code{reduce_data()} is an alias for
\code{reduce_parameters.data.frame()}.
}
\details{
The different methods available are described below:
\subsection{Supervised Methods}{
\itemize{
\item \strong{PCA}: See \code{\link[=principal_components]{principal_components()}}.

\item \strong{cMDS / PCoA}: Classical Multidimensional Scaling (cMDS) takes a
set of dissimilarities (i.e., a distance matrix) and returns a set of points
such that the distances between the points are approximately equal to the
dissimilarities.

\item \strong{DRR}: Dimensionality Reduction via Regression (DRR) is a very
recent technique extending PCA (Laparra et al., 2015). Starting from a
rotated PCA, it predicts redundant information from the remaining components
using non-linear regression. Some of the most notable advantages of
performing DRR are avoidance of multicollinearity between predictors and
overfitting mitigation. DRR tends to perform well when the first principal
component is enough to explain most of the variation in the predictors.
Requires the \pkg{DRR} package to be installed.

\item \strong{ICA}: Performs an Independent Component Analysis using the
FastICA algorithm. Contrary to PCA, which attempts to find uncorrelated
sources (through least squares minimization), ICA attempts to find
independent sources, i.e., the source space that maximizes the
"non-gaussianity" of all sources. Contrary to PCA, ICA does not rank each
source, which makes it a poor tool for dimensionality reduction. Requires the
\pkg{fastICA} package to be installed.
}
}
See also \href{https://easystats.github.io/parameters/articles/parameters_reduction.html}{package vignette}.
}
\examples{
data(iris)
model <- lm(Sepal.Width ~ Species * Sepal.Length + Petal.Width, data = iris)
model
reduce_parameters(model)

out <- reduce_data(iris, method = "PCA", n = "max")
head(out)
}
\references{
\itemize{
\item Nguyen, L. H., \& Holmes, S. (2019). Ten quick tips for effective
dimensionality reduction. PLOS Computational Biology, 15(6).
\item Laparra, V., Malo, J., & Camps-Valls, G. (2015). Dimensionality
reduction via regression in hyperspectral imagery. IEEE Journal of Selected
Topics in Signal Processing, 9(6), 1026-1036.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_parameters.R
\name{bootstrap_parameters}
\alias{bootstrap_parameters}
\title{Parameters bootstrapping}
\usage{
bootstrap_parameters(
  model,
  iterations = 1000,
  centrality = "median",
  ci = 0.95,
  ci_method = "quantile",
  test = "p-value",
  ...
)
}
\arguments{
\item{model}{Statistical model.}

\item{iterations}{The number of draws to simulate/bootstrap.}

\item{centrality}{The point-estimates (centrality indices) to compute.  Character (vector) or list with one or more of these options: \code{"median"}, \code{"mean"}, \code{"MAP"} or \code{"all"}.}

\item{ci}{Value or vector of probability of the CI (between 0 and 1)
to be estimated. Default to \code{.95} (\verb{95\%}).}

\item{ci_method}{The type of index used for Credible Interval. Can be
\code{"HDI"} (default, see \code{\link[bayestestR:hdi]{hdi()}}), \code{"ETI"}
(see \code{\link[bayestestR:eti]{eti()}}), \code{"BCI"} (see
\code{\link[bayestestR:bci]{bci()}}) or \code{"SI"} (see \code{\link[bayestestR:si]{si()}}).}

\item{test}{The indices to compute. Character (vector) with one or more of these options: \code{"p-value"} (or \code{"p"}), \code{"p_direction"} (or \code{"pd"}), \code{"rope"}, \code{"p_map"}, \code{"equivalence_test"} (or \code{"equitest"}), \code{"bayesfactor"} (or \code{"bf"}) or \code{"all"} to compute all tests. For each "test", the corresponding \pkg{bayestestR} function is called (e.g. \code{\link[bayestestR:rope]{bayestestR::rope()}} or \code{\link[bayestestR:p_direction]{bayestestR::p_direction()}}) and its results included in the summary output.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame summarizing the bootstrapped parameters.
}
\description{
Compute bootstrapped parameters and their related indices such as Confidence Intervals (CI) and p-values.
}
\details{
This function first calls \code{\link[=bootstrap_model]{bootstrap_model()}} to generate
bootstrapped coefficients. The resulting replicated for each coefficient
are treated as "distribution", and is passed to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}
to calculate the related indices defined in the \code{"test"} argument.
\cr\cr
Note that that p-values returned here are estimated under the assumption of
\emph{translation equivariance}: that shape of the sampling distribution is
unaffected by the null being true or not. If this assumption does not hold,
p-values can be biased, and it is suggested to use proper permutation tests
to obtain non-parametric p-values.
}
\section{Using with \strong{emmeans}}{

The output can be passed directly to the various functions from the
\strong{emmeans} package, to obtain bootstrapped estimates, contrasts, simple
slopes, etc. and their confidence intervals. These can then be passed to
\code{model_parameter()} to obtain standard errors, p-values, etc. (see
example).
\cr\cr
Note that that p-values returned here are estimated under the assumption of
\emph{translation equivariance}: that shape of the sampling distribution is
unaffected by the null being true or not. If this assumption does not hold,
p-values can be biased, and it is suggested to use proper permutation tests
to obtain non-parametric p-values.
}

\examples{
\dontrun{
if (require("boot", quietly = TRUE)) {
  set.seed(2)
  model <- lm(Sepal.Length ~ Species * Petal.Width, data = iris)
  b <- bootstrap_parameters(model)
  print(b)

  if (require("emmeans")) {
    est <- emmeans(b, trt.vs.ctrl ~ Species)
    print(model_parameters(est))
  }
}
}
}
\references{
Davison, A. C., & Hinkley, D. V. (1997). Bootstrap methods and their application (Vol. 1). Cambridge university press.
}
\seealso{
\code{\link[=bootstrap_model]{bootstrap_model()}}, \code{\link[=simulate_parameters]{simulate_parameters()}}, \code{\link[=simulate_model]{simulate_model()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R, R/format.R, R/print_html.R,
%   R/print_md.R
\name{display.parameters_model}
\alias{display.parameters_model}
\alias{display.parameters_sem}
\alias{display.parameters_efa_summary}
\alias{display.parameters_efa}
\alias{display.equivalence_test_lm}
\alias{format.parameters_model}
\alias{print_html.parameters_model}
\alias{print_md.parameters_model}
\title{Print tables in different output formats}
\usage{
\method{display}{parameters_model}(
  object,
  format = "markdown",
  pretty_names = TRUE,
  split_components = TRUE,
  select = NULL,
  caption = NULL,
  subtitle = NULL,
  footer = NULL,
  align = NULL,
  digits = 2,
  ci_digits = 2,
  p_digits = 3,
  footer_digits = 3,
  ci_brackets = c("(", ")"),
  show_sigma = FALSE,
  show_formula = FALSE,
  zap_small = FALSE,
  verbose = TRUE,
  ...
)

\method{display}{parameters_sem}(
  object,
  format = "markdown",
  digits = 2,
  ci_digits = 2,
  p_digits = 3,
  ci_brackets = c("(", ")"),
  ...
)

\method{display}{parameters_efa_summary}(object, format = "markdown", digits = 3, ...)

\method{display}{parameters_efa}(
  object,
  format = "markdown",
  digits = 2,
  sort = FALSE,
  threshold = NULL,
  labels = NULL,
  ...
)

\method{display}{equivalence_test_lm}(object, format = "markdown", digits = 2, ...)

\method{format}{parameters_model}(
  x,
  pretty_names = TRUE,
  split_components = TRUE,
  select = NULL,
  digits = 2,
  ci_digits = 2,
  p_digits = 3,
  ci_width = NULL,
  ci_brackets = NULL,
  zap_small = FALSE,
  format = NULL,
  groups = NULL,
  ...
)

\method{print_html}{parameters_model}(
  x,
  pretty_names = TRUE,
  split_components = TRUE,
  select = NULL,
  caption = NULL,
  subtitle = NULL,
  footer = NULL,
  align = NULL,
  digits = 2,
  ci_digits = 2,
  p_digits = 3,
  footer_digits = 3,
  ci_brackets = c("(", ")"),
  show_sigma = FALSE,
  show_formula = FALSE,
  zap_small = FALSE,
  groups = NULL,
  verbose = TRUE,
  ...
)

\method{print_md}{parameters_model}(
  x,
  pretty_names = TRUE,
  split_components = TRUE,
  select = NULL,
  caption = NULL,
  subtitle = NULL,
  footer = NULL,
  align = NULL,
  digits = 2,
  ci_digits = 2,
  p_digits = 3,
  footer_digits = 3,
  ci_brackets = c("(", ")"),
  show_sigma = FALSE,
  show_formula = FALSE,
  zap_small = FALSE,
  groups = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{An object returned by \code{\link[=model_parameters]{model_parameters()}},
\code{\link[=simulate_parameters]{simulate_parameters()}},
\code{\link[=equivalence_test.lm]{equivalence_test()}} or
\code{\link[=principal_components]{principal_components()}}.}

\item{format}{String, indicating the output format. Can be \code{"markdown"}
or \code{"html"}.}

\item{pretty_names}{Return "pretty" (i.e. more human readable) parameter
names.}

\item{split_components}{Logical, if \code{TRUE} (default), For models with
multiple components (zero-inflation, smooth terms, ...), each component is
printed in a separate table. If \code{FALSE}, model parameters are printed
in a single table and a \code{Component} column is added to the output.}

\item{select}{Character vector (or numeric index) of column names that should
be printed. If \code{NULL} (default), all columns are printed. The shortcut
\code{select = "minimal"} prints coefficient, confidence intervals and p-values,
while \code{select = "short"} prints coefficient, standard errors and p-values.}

\item{caption}{Table caption as string. If \code{NULL}, no table caption is printed.}

\item{subtitle}{Table title (same as caption) and subtitle, as strings. If \code{NULL},
no title or subtitle is printed, unless it is stored as attributes (\code{table_title},
or its alias \code{table_caption}, and \code{table_subtitle}). If \code{x} is a list of
data frames, \code{caption} may be a list of table captions, one for each table.}

\item{footer}{Table footer, as string. For markdown-formatted tables, table
footers, due to the limitation in markdown rendering, are actually just a
new text line under the table. If \code{x} is a list of data frames, \code{footer}
may be a list of table captions, one for each table.}

\item{align}{Only applies to HTML tables. May be one of \code{"left"},
\code{"right"} or \code{"center"}.}

\item{digits, ci_digits, p_digits}{Number of digits for rounding or
significant figures. May also be \code{"signif"} to return significant
figures or \code{"scientific"} to return scientific notation. Control the
number of digits by adding the value as suffix, e.g. \code{digits = "scientific4"}
to have scientific notation with 4 decimal places, or \code{digits = "signif5"}
for 5 significant figures (see also \code{\link[=signif]{signif()}}).}

\item{footer_digits}{Number of decimal places for values in the footer summary.}

\item{ci_brackets}{Logical, if \code{TRUE} (default), CI-values are
encompassed in square brackets (else in parentheses).}

\item{show_sigma}{Logical, if \code{TRUE}, adds information about the residual
standard deviation.}

\item{show_formula}{Logical, if \code{TRUE}, adds the model formula to the output.}

\item{zap_small}{Logical, if \code{TRUE}, small values are rounded after
\code{digits} decimal places. If \code{FALSE}, values with more decimal
places than \code{digits} are printed in scientific notation.}

\item{verbose}{Toggle messages and warnings.}

\item{...}{Arguments passed to or from other methods.}

\item{sort}{Sort the loadings.}

\item{threshold}{A value between 0 and 1 indicates which (absolute) values
from the loadings should be removed. An integer higher than 1 indicates the
n strongest loadings to retain. Can also be \code{"max"}, in which case it
will only display the maximum loading per variable (the most simple
structure).}

\item{labels}{A character vector containing labels to be added to the
loadings data. Usually, the question related to the item.}

\item{x}{An object returned by \code{\link[=model_parameters]{model_parameters()}}.}

\item{ci_width}{Minimum width of the returned string for confidence
intervals. If not \code{NULL} and width is larger than the string's length,
leading whitespaces are added to the string. If \code{width="auto"}, width
will be set to the length of the longest string.}

\item{groups}{Named list, can be used to group parameters in the printed output.
List elements may either be character vectors that match the name of those
parameters that belong to one group, or list elements can be row numbers
of those parameter rows that should belong to one group. The names of the
list elements will be used as group names, which will be inserted as "header
row". A possible use case might be to emphasize focal predictors and control
variables, see 'Examples'. Parameters will be re-ordered according to the
order used in \code{groups}, while all non-matching parameters will be added
to the end.}
}
\value{
If \code{format = "markdown"}, the return value will be a character
vector in markdown-table format. If \code{format = "html"}, an object of
class \code{gt_tbl}.
}
\description{
Prints tables (i.e. data frame) in different output formats.
\code{print_md()} is a alias for \code{display(format = "markdown")}.
}
\details{
\code{display()} is useful when the table-output from functions,
which is usually printed as formatted text-table to console, should
be formatted for pretty table-rendering in markdown documents, or if
knitted from rmarkdown to PDF or Word files. See
\href{https://easystats.github.io/parameters/articles/model_parameters_formatting.html}{vignette}
for examples.
}
\examples{
model <- lm(mpg ~ wt + cyl, data = mtcars)
mp <- model_parameters(model)
display(mp)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.find_most_common}
\alias{.find_most_common}
\title{Find most common occurence}
\usage{
.find_most_common(x)
}
\description{
Find most common occurence
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_FactoMineR.R, R/methods_psych.R
\name{model_parameters.PCA}
\alias{model_parameters.PCA}
\alias{model_parameters.principal}
\alias{model_parameters.omega}
\title{Parameters from Structural Models (PCA, EFA, ...)}
\usage{
\method{model_parameters}{PCA}(
  model,
  sort = FALSE,
  threshold = NULL,
  labels = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{principal}(
  model,
  sort = FALSE,
  threshold = NULL,
  labels = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{omega}(model, verbose = TRUE, ...)
}
\arguments{
\item{model}{PCA or FA created by the \pkg{psych} or \pkg{FactoMineR}
packages (e.g. through \code{psych::principal},  \code{psych::fa} or \code{psych::omega}).}

\item{sort}{Sort the loadings.}

\item{threshold}{A value between 0 and 1 indicates which (absolute) values
from the loadings should be removed. An integer higher than 1 indicates the
n strongest loadings to retain. Can also be \code{"max"}, in which case it
will only display the maximum loading per variable (the most simple
structure).}

\item{labels}{A character vector containing labels to be added to the
loadings data. Usually, the question related to the item.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame of loadings.
}
\description{
Format structural models from the \pkg{psych} or \pkg{FactoMineR} packages.
}
\details{
For the structural models obtained with \pkg{psych}, the following indices
are present:
\itemize{
\item \strong{Complexity} (\cite{Hoffman's, 1978; Pettersson and Turkheimer,
2010}) represents the number of latent components needed to account for
the observed variables. Whereas a perfect simple structure solution has a
complexity of 1 in that each item would only load on one factor, a
solution with evenly distributed items has a complexity greater than 1.
\item \strong{Uniqueness} represents the variance that is 'unique' to the
variable and not shared with other variables. It is equal to \verb{1 – communality} (variance that is shared with other variables). A uniqueness
of \code{0.20} suggests that \verb{20\%} or that variable's variance is not shared
with other variables in the overall factor model. The greater 'uniqueness'
the lower the relevance of the variable in the factor model.
\item \strong{MSA} represents the Kaiser-Meyer-Olkin Measure of Sampling
Adequacy (\cite{Kaiser and Rice, 1974}) for each item. It indicates
whether there is enough data for each factor give reliable results for the
PCA. The value should be > 0.6, and desirable values are > 0.8
(\cite{Tabachnick and Fidell, 2013}).
}
}
\examples{
\donttest{
library(parameters)
if (require("psych", quietly = TRUE)) {
  # Principal Component Analysis (PCA) ---------
  pca <- psych::principal(attitude)
  model_parameters(pca)

  pca <- psych::principal(attitude, nfactors = 3, rotate = "none")
  model_parameters(pca, sort = TRUE, threshold = 0.2)

  principal_components(attitude, n = 3, sort = TRUE, threshold = 0.2)


  # Exploratory Factor Analysis (EFA) ---------
  efa <- psych::fa(attitude, nfactors = 3)
  model_parameters(efa, threshold = "max", sort = TRUE, labels = as.character(1:ncol(attitude)))


  # Omega ---------
  omega <- psych::omega(mtcars, nfactors = 3)
  params <- model_parameters(omega)
  params
  summary(params)
}

# FactoMineR ---------
if (require("FactoMineR", quietly = TRUE)) {
  model <- FactoMineR::PCA(iris[, 1:4], ncp = 2)
  model_parameters(model)
  attributes(model_parameters(model))$scores

  model <- FactoMineR::FAMD(iris, ncp = 2)
  model_parameters(model)
}
}

}
\references{
\itemize{
\item Kaiser, H.F. and Rice. J. (1974). Little jiffy, mark iv. Educational and
Psychological Measurement, 34(1):111–117
\item Pettersson, E., \& Turkheimer, E. (2010). Item selection, evaluation, and
simple structure in personality data. Journal of research in personality,
44(4), 407-420.
\item Revelle, W. (2016). How To: Use the psych package for Factor Analysis and
data reduction.
\item Tabachnick, B. G., and Fidell, L. S. (2013). Using multivariate statistics
(6th ed.). Boston: Pearson Education.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_cgam.R, R/methods_gam.R,
%   R/methods_quantreg.R
\name{model_parameters.cgam}
\alias{model_parameters.cgam}
\alias{model_parameters.gam}
\alias{model_parameters.rqss}
\title{Parameters from Generalized Additive (Mixed) Models}
\usage{
\method{model_parameters}{cgam}(
  model,
  ci = 0.95,
  ci_method = "residual",
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)

\method{model_parameters}{gam}(
  model,
  ci = 0.95,
  ci_method = "residual",
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)

\method{model_parameters}{rqss}(
  model,
  ci = 0.95,
  ci_method = "residual",
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{A gam/gamm model.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{ci_method}{Method for computing degrees of freedom for
confidence intervals (CI) and the related p-values. Allowed are following
options (which vary depending on the model class): \code{"residual"},
\code{"normal"}, \code{"likelihood"}, \code{"satterthwaite"}, \code{"kenward"}, \code{"wald"},
\code{"profile"}, \code{"boot"}, \code{"uniroot"}, \code{"ml1"}, \code{"betwithin"}, \code{"hdi"},
\code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, or \code{"bcai"}. See section
\emph{Confidence intervals and approximation of degrees of freedom} in
\code{\link[=model_parameters]{model_parameters()}} for further details. When \code{ci_method=NULL}, in most
cases \code{"wald"} is used then.}

\item{bootstrap}{Should estimates be based on bootstrapped model? If
\code{TRUE}, then arguments of \link[=model_parameters.stanreg]{Bayesian regressions} apply (see also
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}).}

\item{iterations}{The number of bootstrap replicates. This only apply in the
case of bootstrapped frequentist models.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{robust}{Logical, if \code{TRUE}, robust standard errors are calculated
(if possible), and confidence intervals and p-values are based on these
robust standard errors. Additional arguments like \code{vcov_estimation} or
\code{vcov_type} are passed down to other methods, see
\code{\link[=standard_error_robust]{standard_error_robust()}} for details
and \href{https://easystats.github.io/parameters/articles/model_parameters_robust.html}{this vignette}
for working examples.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods. For instance, when
\code{bootstrap = TRUE}, arguments like \code{type} or \code{parallel} are
passed down to \code{bootstrap_model()}, and arguments like \code{ci_method}
are passed down to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Extract and compute indices and measures to describe parameters
of generalized additive models (GAM(M)s).
}
\details{
The reporting of degrees of freedom \emph{for the spline terms}
slightly differs from the output of \code{summary(model)}, for example in the
case of \code{mgcv::gam()}. The \emph{estimated degrees of freedom}, column
\code{edf} in the summary-output, is named \code{df} in the returned data
frame, while the column \code{df_error} in the returned data frame refers to
the residual degrees of freedom that are returned by \code{df.residual()}.
Hence, the values in the the column \code{df_error} differ from the column
\code{Ref.df} from the summary, which is intentional, as these reference
degrees of freedom \dQuote{is not very interpretable}
(\href{https://stat.ethz.ch/pipermail/r-help/2019-March/462135.html}{web}).
}
\examples{
library(parameters)
if (require("mgcv")) {
  dat <- gamSim(1, n = 400, dist = "normal", scale = 2)
  model <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
  model_parameters(model)
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to rename
columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_model_parameters.R, R/methods_mfx.R
\name{model_parameters.default}
\alias{model_parameters.default}
\alias{model_parameters.glm}
\alias{model_parameters.logitor}
\alias{model_parameters.poissonmfx}
\alias{model_parameters.betamfx}
\title{Parameters from (General) Linear Models}
\usage{
\method{model_parameters}{default}(
  model,
  ci = 0.95,
  ci_method = NULL,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  summary = FALSE,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)

\method{model_parameters}{glm}(
  model,
  ci = 0.95,
  ci_method = NULL,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  summary = FALSE,
  verbose = TRUE,
  df_method = ci_method,
  ...
)

\method{model_parameters}{logitor}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = TRUE,
  robust = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{poissonmfx}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  component = c("all", "conditional", "marginal"),
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{betamfx}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  component = c("all", "conditional", "precision", "marginal"),
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{Model object.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{ci_method}{Method for computing degrees of freedom for
confidence intervals (CI) and the related p-values. Allowed are following
options (which vary depending on the model class): \code{"residual"},
\code{"normal"}, \code{"likelihood"}, \code{"satterthwaite"}, \code{"kenward"}, \code{"wald"},
\code{"profile"}, \code{"boot"}, \code{"uniroot"}, \code{"ml1"}, \code{"betwithin"}, \code{"hdi"},
\code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, or \code{"bcai"}. See section
\emph{Confidence intervals and approximation of degrees of freedom} in
\code{\link[=model_parameters]{model_parameters()}} for further details. When \code{ci_method=NULL}, in most
cases \code{"wald"} is used then.}

\item{bootstrap}{Should estimates be based on bootstrapped model? If
\code{TRUE}, then arguments of \link[=model_parameters.stanreg]{Bayesian regressions} apply (see also
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}).}

\item{iterations}{The number of bootstrap replicates. This only apply in the
case of bootstrapped frequentist models.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{robust}{Logical, if \code{TRUE}, robust standard errors are calculated
(if possible), and confidence intervals and p-values are based on these
robust standard errors. Additional arguments like \code{vcov_estimation} or
\code{vcov_type} are passed down to other methods, see
\code{\link[=standard_error_robust]{standard_error_robust()}} for details
and \href{https://easystats.github.io/parameters/articles/model_parameters_robust.html}{this vignette}
for working examples.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}

\item{summary}{Logical, if \code{TRUE}, prints summary information about the
model (model formula, number of observations, residual standard deviation
and more).}

\item{keep, drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods. For instance, when
\code{bootstrap = TRUE}, arguments like \code{type} or \code{parallel} are
passed down to \code{bootstrap_model()}, and arguments like \code{ci_method}
are passed down to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}.}

\item{df_method}{Deprecated. Please use \code{ci_method}.}

\item{component}{Model component for which parameters should be shown. May be
one of \code{"conditional"}, \code{"precision"} (\pkg{betareg}),
\code{"scale"} (\pkg{ordinal}), \code{"extra"} (\pkg{glmx}),
\code{"marginal"} (\pkg{mfx}), \code{"conditional"} or \code{"full"} (for
\code{MuMIn::model.avg()}) or \code{"all"}.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Extract and compute indices and measures to describe parameters of (general)
linear models (GLMs).
}
\section{Confidence intervals and approximation of degrees of freedom}{

There are different ways of approximating the degrees of freedom depending
on different assumptions about the nature of the model and its sampling
distribution. The \code{ci_method} argument modulates the method for computing degrees
of freedom (df) that are used to calculate confidence intervals (CI) and the
related p-values. Following options are allowed, depending on the model
class:

\strong{Classical methods:}

Classical inference is generally based on the \strong{Wald method}.
The Wald approach to inference computes a test statistic by dividing the
parameter estimate by its standard error (Coefficient / SE),
then comparing this statistic against a t- or normal distribution.
This approach can be used to compute CIs and p-values.

\code{"wald"}:
\itemize{
\item Applies to \emph{non-Bayesian models}. For \emph{linear models}, CIs
computed using the Wald method (SE and a \emph{t-distribution with residual df});
p-values computed using the Wald method with a \emph{t-distribution with residual df}.
For other models, CIs computed using the Wald method (SE and a \emph{normal distribution});
p-values computed using the Wald method with a \emph{normal distribution}.
}

\code{"normal"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a normal distribution.
}

\code{"residual"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a \emph{t-distribution with residual df} when possible. If the
residual df for a model cannot be determined, a normal distribution is
used instead.
}

\strong{Methods for mixed models:}

Compared to fixed effects (or single-level) models, determining appropriate
df for Wald-based inference in mixed models is more difficult.
See \href{https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#what-are-the-p-values-listed-by-summaryglmerfit-etc.-are-they-reliable}{the R GLMM FAQ}
for a discussion.

Several approximate methods for computing df are available, but you should
also consider instead using profile likelihood (\code{"profile"}) or bootstrap ("\verb{boot"})
CIs and p-values instead.

\code{"satterthwaite"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the
Wald method (SE and a \emph{t-distribution with Satterthwaite df}); p-values
computed using the Wald method with a \emph{t-distribution with Satterthwaite df}.
}

\code{"kenward"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (\emph{Kenward-Roger SE} and a \emph{t-distribution with Kenward-Roger df});
p-values computed using the Wald method with \emph{Kenward-Roger SE and t-distribution with Kenward-Roger df}.
}

\code{"ml1"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (SE and a \emph{t-distribution with m-l-1 approximated df}); p-values
computed using the Wald method with a \emph{t-distribution with m-l-1 approximated df}.
See \code{\link[=ci_ml1]{ci_ml1()}}.
}

\code{"betwithin"}
\itemize{
\item Applies to \emph{linear mixed models} and \emph{generalized linear mixed models}.
CIs computed using the Wald method (SE and a \emph{t-distribution with between-within df});
p-values computed using the Wald method with a \emph{t-distribution with between-within df}.
See \code{\link[=ci_betwithin]{ci_betwithin()}}.
}

\strong{Likelihood-based methods:}

Likelihood-based inference is based on comparing the likelihood for the
maximum-likelihood estimate to the the likelihood for models with one or more
parameter values changed (e.g., set to zero or a range of alternative values).
Likelihood ratios for the maximum-likelihood and alternative models are compared
to a \eqn{\chi}-squared distribution to compute CIs and p-values.

\code{"profile"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glm}, \code{polr} or \code{glmmTMB}.
CIs computed by \emph{profiling the likelihood curve for a parameter}, using
linear interpolation to find where likelihood ratio equals a critical value;
p-values computed using the Wald method with a \emph{normal-distribution} (note:
this might change in a future update!)
}

\code{"uniroot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glmmTMB}. CIs
computed by \emph{profiling the likelihood curve for a parameter}, using root
finding to find where likelihood ratio equals a critical value; p-values
computed using the Wald method with a \emph{normal-distribution} (note: this
might change in a future update!)
}

\strong{Methods for bootstrapped or Bayesian models:}

Bootstrap-based inference is based on \strong{resampling} and refitting the model
to the resampled datasets. The distribution of parameter estimates across
resampled datasets is used to approximate the parameter's sampling
distribution. Depending on the type of model, several different methods for
bootstrapping and constructing CIs and p-values from the bootstrap
distribution are available.

For Bayesian models, inference is based on drawing samples from the model
posterior distribution.

\code{"quantile"} (or \code{"eti"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{equal tailed intervals} using the quantiles of the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:eti]{bayestestR::eti()}}.
}

\code{"hdi"}
\itemize{
\item Applies to \emph{all models (including Bayesian models)}. For non-Bayesian
models, only applies if \code{bootstrap = TRUE}. CIs computed as \emph{highest density intervals}
for the bootstrap or posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:hdi]{bayestestR::hdi()}}.
}

\code{"bci"} (or \code{"bcai"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{bias corrected and accelerated intervals} for the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:bci]{bayestestR::bci()}}.
}

\code{"si"}
\itemize{
\item Applies to \emph{Bayesian models} with proper priors. CIs computed as
\emph{support intervals} comparing the posterior samples against the prior samples;
p-values are based on the \emph{probability of direction}. See \code{\link[bayestestR:si]{bayestestR::si()}}.
}

\code{"boot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{merMod}. CIs computed
using \emph{parametric bootstrapping} (simulating data from the fitted model);
p-values computed using the Wald method with a \emph{normal-distribution)}
(note: this might change in a future update!).
}

For all iteration-based methods other than \code{"boot"}
(\code{"hdi"}, \code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, \code{"bcai"}),
p-values are based on the probability of direction (\code{\link[bayestestR:p_direction]{bayestestR::p_direction()}}),
which is converted into a p-value using \code{\link[bayestestR:pd_to_p]{bayestestR::pd_to_p()}}.
}

\examples{
library(parameters)
model <- lm(mpg ~ wt + cyl, data = mtcars)

model_parameters(model)

# bootstrapped parameters
model_parameters(model, bootstrap = TRUE)

# standardized parameters
model_parameters(model, standardize = "refit")

# different p-value style in output
model_parameters(model, p_digits = 5)
model_parameters(model, digits = 3, ci_digits = 4, p_digits = "scientific")
\donttest{
# logistic regression model
model <- glm(vs ~ wt + cyl, data = mtcars, family = "binomial")
model_parameters(model)

# show odds ratio / exponentiated coefficients
model_parameters(model, exponentiate = TRUE)
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to
rename columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_DirichletReg.R, R/methods_averaging.R,
%   R/methods_betareg.R, R/methods_cgam.R, R/methods_ordinal.R
\name{p_value.DirichletRegModel}
\alias{p_value.DirichletRegModel}
\alias{p_value.averaging}
\alias{p_value.betareg}
\alias{p_value.cgam}
\alias{p_value.clm2}
\title{p-values for Models with Special Components}
\usage{
\method{p_value}{DirichletRegModel}(model, component = c("all", "conditional", "precision"), ...)

\method{p_value}{averaging}(model, component = c("conditional", "full"), ...)

\method{p_value}{betareg}(model, component = c("all", "conditional", "precision"), ...)

\method{p_value}{cgam}(model, component = c("all", "conditional", "smooth_terms"), ...)

\method{p_value}{clm2}(model, component = c("all", "conditional", "scale"), ...)
}
\arguments{
\item{model}{A statistical model.}

\item{component}{Should all parameters, parameters for the conditional model,
precision- or scale-component or smooth_terms be returned? \code{component}
may be one of \code{"conditional"}, \code{"precision"}, \code{"scale"},
\code{"smooth_terms"}, \code{"full"} or \code{"all"} (default).}

\item{...}{Arguments passed down to \code{standard_error_robust()} when confidence
intervals or p-values based on robust standard errors should be computed.
Only available for models where \code{method = "robust"} is supported.}
}
\value{
The p-values.
}
\description{
This function attempts to return, or compute, p-values of models
with special model components.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters_type.R
\name{parameters_type}
\alias{parameters_type}
\title{Type of model parameters}
\usage{
parameters_type(model, ...)
}
\arguments{
\item{model}{A statistical model.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame.
}
\description{
In a regression model, the parameters do not all have the meaning. For
instance, the intercept has to be interpreted as theoretical outcome value
under some conditions (when predictors are set to 0), whereas other
coefficients are to be interpreted as amounts of change. Others, such as
interactions, represent changes in another of the parameter. The
\code{parameters_type} function attempts to retrieve information and meaning
of parameters. It outputs a dataframe of information for each parameters,
such as the \code{Type} (whether the parameter corresponds to a factor or a
numeric predictor, or whether it is a (regular) interaction or a nested
one), the \code{Link} (whether the parameter can be interpreted as a mean
value, the slope of an association or a difference between two levels) and,
in the case of interactions, which other parameters is impacted by which
parameter.
}
\examples{
library(parameters)

model <- lm(Sepal.Length ~ Petal.Length + Species, data = iris)
parameters_type(model)

model <- lm(Sepal.Length ~ Species + poly(Sepal.Width, 2), data = iris)
parameters_type(model)

model <- lm(Sepal.Length ~ Species + poly(Sepal.Width, 2, raw = TRUE), data = iris)
parameters_type(model)

# Interactions
model <- lm(Sepal.Length ~ Sepal.Width * Species, data = iris)
parameters_type(model)

model <- lm(Sepal.Length ~ Sepal.Width * Species * Petal.Length, data = iris)
parameters_type(model)

model <- lm(Sepal.Length ~ Species * Sepal.Width, data = iris)
parameters_type(model)

model <- lm(Sepal.Length ~ Species / Sepal.Width, data = iris)
parameters_type(model)


# Complex interactions
data <- iris
data$fac2 <- ifelse(data$Sepal.Width > mean(data$Sepal.Width), "A", "B")
model <- lm(Sepal.Length ~ Species / fac2 / Petal.Length, data = data)
parameters_type(model)

model <- lm(Sepal.Length ~ Species / fac2 * Petal.Length, data = data)
parameters_type(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_parameters.R
\name{random_parameters}
\alias{random_parameters}
\title{Summary information from random effects}
\usage{
random_parameters(model, component = "conditional")
}
\arguments{
\item{model}{A mixed effects model (including \code{stanreg} models).}

\item{component}{Should all parameters, parameters for the conditional model,
or for the zero-inflated part of the model be returned? Applies to models
with zero-inflated component. \code{component} may be one of
\code{"conditional"} (default), \code{"zi"} or \code{"zero_inflated"}.
May be abbreviated.}
}
\value{
A data frame with random effects statistics for the variance components,
including number of levels per random effect group, as well as complete
observations in the model.
}
\description{
This function extracts the different variance components of a
mixed model and returns the result as a data frame.
}
\details{
The variance components are obtained from
\code{\link[insight:get_variance]{insight::get_variance()}} and are denoted as following:
\subsection{Within-group (or residual) variance}{
The residual variance, \ifelse{html}{\out{&sigma;<sup>2</sup><sub>&epsilon;</sub>}}{\eqn{\sigma^2_\epsilon}},
is the sum of the distribution-specific variance and the variance due to additive dispersion.
It indicates the \emph{within-group variance}.
}
\subsection{Between-group random intercept variance}{
The random intercept variance, or \emph{between-group} variance
for the intercept (\ifelse{html}{\out{&tau;<sub>00</sub>}}{\eqn{\tau_{00}}}),
is obtained from \code{VarCorr()}. It indicates how much groups
or subjects differ from each other.
}
\subsection{Between-group random slope variance}{
The random slope variance, or \emph{between-group} variance
for the slopes (\ifelse{html}{\out{&tau;<sub>11</sub>}}{\eqn{\tau_{11}}})
is obtained from \code{VarCorr()}. This measure is only available
for mixed models with random slopes. It indicates how much groups
or subjects differ from each other according to their slopes.
}
\subsection{Random slope-intercept correlation}{
The random slope-intercept correlation
(\ifelse{html}{\out{&rho;<sub>01</sub>}}{\eqn{\rho_{01}}})
is obtained from \code{VarCorr()}. This measure is only available
for mixed models with random intercepts and slopes.
}
\strong{Note:} For the within-group and between-group variance, variance
and standard deviations (which are simply the square root of the variance)
are shown.
}
\examples{
if (require("lme4")) {
  data(sleepstudy)
  model <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
  random_parameters(model)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_DirichletReg.R, R/methods_bife.R,
%   R/methods_brglm2.R, R/methods_mlm.R, R/methods_ordinal.R
\name{model_parameters.DirichletRegModel}
\alias{model_parameters.DirichletRegModel}
\alias{model_parameters.bifeAPEs}
\alias{model_parameters.bracl}
\alias{model_parameters.mlm}
\alias{model_parameters.clm2}
\title{Parameters from multinomial or cumulative link models}
\usage{
\method{model_parameters}{DirichletRegModel}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  component = c("all", "conditional", "precision"),
  standardize = NULL,
  exponentiate = FALSE,
  verbose = TRUE,
  ...
)

\method{model_parameters}{bifeAPEs}(model, ...)

\method{model_parameters}{bracl}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{mlm}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)

\method{model_parameters}{clm2}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  component = c("all", "conditional", "scale"),
  standardize = NULL,
  exponentiate = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{A model with multinomial or categorical response value.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{bootstrap}{Should estimates be based on bootstrapped model? If
\code{TRUE}, then arguments of \link[=model_parameters.stanreg]{Bayesian regressions} apply (see also
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}).}

\item{iterations}{The number of bootstrap replicates. This only apply in the
case of bootstrapped frequentist models.}

\item{component}{Model component for which parameters should be shown. May be
one of \code{"conditional"}, \code{"precision"} (\pkg{betareg}),
\code{"scale"} (\pkg{ordinal}), \code{"extra"} (\pkg{glmx}),
\code{"marginal"} (\pkg{mfx}), \code{"conditional"} or \code{"full"} (for
\code{MuMIn::model.avg()}) or \code{"all"}.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods. For instance, when
\code{bootstrap = TRUE}, arguments like \code{type} or \code{parallel} are
passed down to \code{bootstrap_model()}, and arguments like \code{ci_method}
are passed down to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from multinomial or cumulative link models
}
\details{
Multinomial or cumulative link models, i.e. models where the
response value (dependent variable) is categorical and has more than two
levels, usually return coefficients for each response level. Hence, the
output from \code{model_parameters()} will split the coefficient tables
by the different levels of the model's response.
}
\examples{
library(parameters)
if (require("brglm2", quietly = TRUE)) {
  data("stemcell")
  model <- bracl(
    research ~ as.numeric(religion) + gender,
    weights = frequency,
    data = stemcell,
    type = "ML"
  )
  model_parameters(model)
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to rename
columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_clusters.R, R/n_clusters_easystats.R
\name{n_clusters}
\alias{n_clusters}
\alias{n_clusters_elbow}
\alias{n_clusters_gap}
\alias{n_clusters_silhouette}
\alias{n_clusters_dbscan}
\alias{n_clusters_hclust}
\title{Find number of clusters in your data}
\usage{
n_clusters(
  x,
  standardize = TRUE,
  include_factors = FALSE,
  package = c("easystats", "NbClust", "mclust"),
  fast = TRUE,
  nbclust_method = "kmeans",
  n_max = 10,
  ...
)

n_clusters_elbow(
  x,
  standardize = TRUE,
  include_factors = FALSE,
  clustering_function = stats::kmeans,
  n_max = 10,
  ...
)

n_clusters_gap(
  x,
  standardize = TRUE,
  include_factors = FALSE,
  clustering_function = stats::kmeans,
  n_max = 10,
  gap_method = "firstSEmax",
  ...
)

n_clusters_silhouette(
  x,
  standardize = TRUE,
  include_factors = FALSE,
  clustering_function = stats::kmeans,
  n_max = 10,
  ...
)

n_clusters_dbscan(
  x,
  standardize = TRUE,
  include_factors = FALSE,
  method = c("kNN", "SS"),
  min_size = 0.1,
  eps_n = 50,
  eps_range = c(0.1, 3),
  ...
)

n_clusters_hclust(
  x,
  standardize = TRUE,
  include_factors = FALSE,
  distance_method = "correlation",
  hclust_method = "average",
  ci = 0.95,
  iterations = 100,
  ...
)
}
\arguments{
\item{x}{A data frame.}

\item{standardize}{Standardize the dataframe before clustering (default).}

\item{include_factors}{Logical, if \code{TRUE}, factors are converted to numerical
values in order to be included in the data for determining the number of
clusters. By default, factors are removed, because most methods that
determine the number of clusters need numeric input only.}

\item{package}{Package from which methods are to be called to determine the
number of clusters. Can be \code{"all"} or a vector containing
\code{"NbClust"}, \code{"mclust"}, \code{"cluster"} and \code{"M3C"}.}

\item{fast}{If \code{FALSE}, will compute 4 more indices (sets \code{index = "allong"}
in \code{NbClust}). This has been deactivated by default as it is
computationally heavy.}

\item{nbclust_method}{The clustering method (passed to \code{NbClust::NbClust()}
as \code{method}).}

\item{n_max}{Maximal number of clusters to test.}

\item{...}{Arguments passed to or from other methods.}

\item{clustering_function, gap_method}{Other arguments passed to other
functions. \code{clustering_function} is used by \code{fviz_nbclust} and
can be \code{kmeans}, code{cluster::pam}, code{cluster::clara},
code{cluster::fanny}, and more. \code{gap_method} is used by
\code{cluster::maxSE} to extract the optimal numbers of clusters (see its
\code{method} argument).}

\item{method, min_size, eps_n, eps_range}{Arguments for DBSCAN algorithm.}

\item{distance_method}{The distance method (passed to \code{\link[=dist]{dist()}}). Used by
algorithms relying on the distance matrix, such as \code{hclust} or \code{dbscan}.}

\item{hclust_method}{The hierarchical clustering method (passed to \code{\link[=hclust]{hclust()}}).}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{iterations}{The number of bootstrap replicates. This only apply in the
case of bootstrapped frequentist models.}
}
\description{
Similarly to \code{\link[=n_factors]{n_factors()}} for factor / principal component analysis,
\code{n_clusters} is the main function to find out the optimal numbers of clusters
present in the data based on the maximum consensus of a large number of
methods.
\cr
Essentially, there exist many methods to determine the optimal number of
clusters, each with pros and cons, benefits and limitations. The main
\code{n_clusters} function proposes to run all of them, and find out the number of
clusters that is suggested by the majority of methods (in case of ties, it
will select the most parsimonious solution with fewer clusters).
\cr
Note that we also implement some specific, commonly used methods, like the
Elbow or the Gap method, with their own visualization functionalities. See
the examples below for more details.
}
\note{
There is also a \href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
\dontrun{
library(parameters)

# The main 'n_clusters' function ===============================
if (require("mclust", quietly = TRUE) && require("NbClust", quietly = TRUE) &&
  require("cluster", quietly = TRUE) && require("see", quietly = TRUE)) {
  n <- n_clusters(iris[, 1:4], package = c("NbClust", "mclust", "cluster"))
  n
  summary(n)
  as.data.frame(n)
  plot(n)

  # The following runs all the method but it significantly slower
  # n_clusters(iris[1:4], standardize = FALSE, package = "all", fast = FALSE)
}
}
\donttest{
#
# Specific Methods =========================
# Elbow method --------------------
if (require("openxlsx", quietly = TRUE) &&
  require("see", quietly = TRUE) &&
  require("factoextra", quietly = TRUE)) {
  x <- n_clusters_elbow(iris[1:4])
  x
  as.data.frame(x)
  plot(x)
}
}
\donttest{
#
# Gap method --------------------
if (require("see", quietly = TRUE) &&
  require("cluster", quietly = TRUE) &&
  require("factoextra", quietly = TRUE)) {
  x <- n_clusters_gap(iris[1:4])
  x
  as.data.frame(x)
  plot(x)
}
}
\donttest{
#
# Silhouette method --------------------------
if (require("factoextra", quietly = TRUE)) {
  x <- n_clusters_silhouette(iris[1:4])
  x
  as.data.frame(x)
  plot(x)
}
}
\donttest{
#
if (require("dbscan", quietly = TRUE)) {
  # DBSCAN method -------------------------
  # NOTE: This actually primarily estimates the 'eps' parameter, the number of
  # clusters is a side effect (it's the number of clusters corresponding to
  # this 'optimal' EPS parameter).
  x <- n_clusters_dbscan(iris[1:4], method = "kNN", min_size = 0.05) # 5 percent
  x
  head(as.data.frame(x))
  plot(x)

  x <- n_clusters_dbscan(iris[1:4], method = "SS", eps_n = 100, eps_range = c(0.1, 2))
  x
  head(as.data.frame(x))
  plot(x)
}
}
\donttest{
#
# hclust method -------------------------------
if (require("pvclust", quietly = TRUE) &&
  getRversion() >= "3.6.0") {
  # iterations should be higher for real analyses
  x <- n_clusters_hclust(iris[1:4], iterations = 50, ci = 0.90)
  x
  head(as.data.frame(x), n = 10) # Print 10 first rows
  plot(x)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_factorstructure.R
\name{check_factorstructure}
\alias{check_factorstructure}
\title{Check suitability of data for Factor Analysis (FA)}
\usage{
check_factorstructure(x, ...)
}
\arguments{
\item{x}{A dataframe.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A list of lists of indices related to sphericity and KMO.
}
\description{
This checks whether the data is appropriate for Factor Analysis (FA) by
running the \link[=check_sphericity_bartlett]{Bartlett's Test of Sphericity} and
the \link[=check_kmo]{Kaiser, Meyer, Olkin (KMO) Measure of Sampling Adequacy (MSA)}.
}
\examples{
library(parameters)
check_factorstructure(mtcars)
}
\seealso{
\code{\link[=check_kmo]{check_kmo()}}, \code{\link[=check_sphericity_bartlett]{check_sphericity_bartlett()}} and \code{\link[=check_clusterstructure]{check_clusterstructure()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_cplm.R
\name{model_parameters.zcpglm}
\alias{model_parameters.zcpglm}
\title{Parameters from Zero-Inflated Models}
\usage{
\method{model_parameters}{zcpglm}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  component = c("all", "conditional", "zi", "zero_inflated"),
  standardize = NULL,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{A model with zero-inflation component.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{bootstrap}{Should estimates be based on bootstrapped model? If
\code{TRUE}, then arguments of \link[=model_parameters.stanreg]{Bayesian regressions} apply (see also
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}).}

\item{iterations}{The number of bootstrap replicates. This only apply in the
case of bootstrapped frequentist models.}

\item{component}{Model component for which parameters should be shown. May be
one of \code{"conditional"}, \code{"precision"} (\pkg{betareg}),
\code{"scale"} (\pkg{ordinal}), \code{"extra"} (\pkg{glmx}),
\code{"marginal"} (\pkg{mfx}), \code{"conditional"} or \code{"full"} (for
\code{MuMIn::model.avg()}) or \code{"all"}.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{robust}{Logical, if \code{TRUE}, robust standard errors are calculated
(if possible), and confidence intervals and p-values are based on these
robust standard errors. Additional arguments like \code{vcov_estimation} or
\code{vcov_type} are passed down to other methods, see
\code{\link[=standard_error_robust]{standard_error_robust()}} for details
and \href{https://easystats.github.io/parameters/articles/model_parameters_robust.html}{this vignette}
for working examples.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods. For instance, when
\code{bootstrap = TRUE}, arguments like \code{type} or \code{parallel} are
passed down to \code{bootstrap_model()}, and arguments like \code{ci_method}
are passed down to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from zero-inflated models (from packages like \pkg{pscl},
\pkg{cplm} or \pkg{countreg}).
}
\examples{
library(parameters)
if (require("pscl")) {
  data("bioChemists")
  model <- zeroinfl(art ~ fem + mar + kid5 + ment | kid5 + phd, data = bioChemists)
  model_parameters(model)
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to rename
columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robust_estimation.R
\name{standard_error_robust}
\alias{standard_error_robust}
\alias{p_value_robust}
\alias{ci_robust}
\title{Robust estimation}
\usage{
standard_error_robust(
  model,
  vcov_estimation = "HC",
  vcov_type = NULL,
  vcov_args = NULL,
  component = "conditional",
  ...
)

p_value_robust(
  model,
  vcov_estimation = "HC",
  vcov_type = NULL,
  vcov_args = NULL,
  component = "conditional",
  method = NULL,
  ...
)

ci_robust(
  model,
  ci = 0.95,
  method = NULL,
  vcov_estimation = "HC",
  vcov_type = NULL,
  vcov_args = NULL,
  component = "conditional",
  ...
)
}
\arguments{
\item{model}{A model.}

\item{vcov_estimation}{String, indicating the suffix of the
\verb{vcov*()}-function from the \pkg{sandwich} or \pkg{clubSandwich}
package, e.g. \code{vcov_estimation = "CL"} (which calls
\code{\link[sandwich:vcovCL]{sandwich::vcovCL()}} to compute clustered covariance matrix
estimators), or \code{vcov_estimation = "HC"} (which calls
\code{\link[sandwich:vcovHC]{sandwich::vcovHC()}} to compute
heteroskedasticity-consistent covariance matrix estimators).}

\item{vcov_type}{Character vector, specifying the estimation type for the
robust covariance matrix estimation (see
\code{\link[sandwich:vcovHC]{sandwich::vcovHC()}} or \code{clubSandwich::vcovCR()}
for details). Passed down as \code{type} argument to the related \verb{vcov*()}-function
from the  \pkg{sandwich} or \pkg{clubSandwich} package and hence will be
ignored if there is no \code{type} argument (e.g., \code{sandwich::vcovHAC()} will
ignore that argument).}

\item{vcov_args}{List of named vectors, used as additional arguments that are
passed down to the \pkg{sandwich}-function specified in
\code{vcov_estimation}.}

\item{component}{Should all parameters or parameters for specific model
components be returned?}

\item{...}{Arguments passed to or from other methods. For
\code{standard_error()}, if \code{method = "robust"}, arguments
\code{vcov_estimation}, \code{vcov_type} and \code{vcov_args} can be passed
down to \code{standard_error_robust()}.}

\item{method}{Method for computing degrees of freedom for
confidence intervals (CI) and the related p-values. Allowed are following
options (which vary depending on the model class): \code{"residual"},
\code{"normal"}, \code{"likelihood"}, \code{"satterthwaite"}, \code{"kenward"}, \code{"wald"},
\code{"profile"}, \code{"boot"}, \code{"uniroot"}, \code{"ml1"}, \code{"betwithin"}, \code{"hdi"},
\code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, or \code{"bcai"}. See section
\emph{Confidence intervals and approximation of degrees of freedom} in
\code{\link[=model_parameters]{model_parameters()}} for further details.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}
}
\value{
A data frame.
}
\description{
\code{standard_error_robust()}, \code{ci_robust()} and \code{p_value_robust()}
attempt to return indices based on robust estimation of the variance-covariance
matrix, using the packages \pkg{sandwich} and \pkg{clubSandwich}.
}
\note{
These functions rely on the \pkg{sandwich} or \pkg{clubSandwich} package
(the latter if \code{vcov_estimation = "CR"} for cluster-robust standard errors)
and will thus only work for those models supported by those packages.
}
\examples{
if (require("sandwich", quietly = TRUE)) {
  # robust standard errors, calling sandwich::vcovHC(type="HC3") by default
  model <- lm(Petal.Length ~ Sepal.Length * Species, data = iris)
  standard_error_robust(model)
}
\dontrun{
if (require("clubSandwich", quietly = TRUE)) {
  # cluster-robust standard errors, using clubSandwich
  iris$cluster <- factor(rep(LETTERS[1:8], length.out = nrow(iris)))
  standard_error_robust(
    model,
    vcov_type = "CR2",
    vcov_args = list(cluster = iris$cluster)
  )
}
}
}
\seealso{
Working examples cam be found \href{https://easystats.github.io/parameters/articles/model_parameters_robust.html}{in this vignette}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_BayesFM.R
\name{model_parameters.befa}
\alias{model_parameters.befa}
\title{Parameters from Bayesian Exploratory Factor Analysis}
\usage{
\method{model_parameters}{befa}(
  model,
  sort = FALSE,
  centrality = "median",
  dispersion = FALSE,
  ci = 0.95,
  ci_method = "hdi",
  test = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{Bayesian EFA created by the \code{BayesFM::befa}.}

\item{sort}{Sort the loadings.}

\item{centrality}{The point-estimates (centrality indices) to compute.  Character (vector) or list with one or more of these options: \code{"median"}, \code{"mean"}, \code{"MAP"} or \code{"all"}.}

\item{dispersion}{Logical, if \code{TRUE}, computes indices of dispersion related to the estimate(s) (\code{SD} and \code{MAD} for \code{mean} and \code{median}, respectively).}

\item{ci}{Value or vector of probability of the CI (between 0 and 1)
to be estimated. Default to \code{.95} (\verb{95\%}).}

\item{ci_method}{The type of index used for Credible Interval. Can be
\code{"HDI"} (default, see \code{\link[bayestestR:hdi]{hdi()}}), \code{"ETI"}
(see \code{\link[bayestestR:eti]{eti()}}), \code{"BCI"} (see
\code{\link[bayestestR:bci]{bci()}}) or \code{"SI"} (see \code{\link[bayestestR:si]{si()}}).}

\item{test}{The indices of effect existence to compute. Character (vector) or
list with one or more of these options: \code{"p_direction"} (or \code{"pd"}),
\code{"rope"}, \code{"p_map"}, \code{"equivalence_test"} (or \code{"equitest"}),
\code{"bayesfactor"} (or \code{"bf"}) or \code{"all"} to compute all tests.
For each "test", the corresponding \pkg{bayestestR} function is called
(e.g. \code{\link[bayestestR:rope]{rope()}} or \code{\link[bayestestR:p_direction]{p_direction()}}) and its results
included in the summary output.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame of loadings.
}
\description{
Format Bayesian Exploratory Factor Analysis objects from the BayesFM package.
}
\examples{
library(parameters)
\donttest{
if (require("BayesFM")) {
  efa <- BayesFM::befa(mtcars, iter = 1000)
  results <- model_parameters(efa, sort = TRUE)
  results
  efa_to_cfa(results)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equivalence_test.R
\name{equivalence_test.lm}
\alias{equivalence_test.lm}
\alias{equivalence_test.merMod}
\title{Equivalence test}
\usage{
\method{equivalence_test}{lm}(
  x,
  range = "default",
  ci = 0.95,
  rule = "classic",
  verbose = TRUE,
  ...
)

\method{equivalence_test}{merMod}(
  x,
  range = "default",
  ci = 0.95,
  rule = "classic",
  effects = c("fixed", "random"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{A statistical model.}

\item{range}{The range of practical equivalence of an effect. May be
\code{"default"}, to automatically define this range based on properties of the
model's data.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{rule}{Character, indicating the rules when testing for practical
equivalence. Can be \code{"bayes"}, \code{"classic"} or \code{"cet"}. See
'Details'.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods.}

\item{effects}{Should parameters for fixed effects (\code{"fixed"}), random
effects (\code{"random"}), or both (\code{"all"}) be returned? Only applies
to mixed models. May be abbreviated. If the calculation of random effects
parameters takes too long, you may use \code{effects = "fixed"}.}
}
\value{
A data frame.
}
\description{
Compute the (conditional) equivalence test for frequentist models.
}
\details{
In classical null hypothesis significance testing (NHST) within a frequentist
framework, it is not possible to accept the null hypothesis, H0 - unlike
in Bayesian statistics, where such probability statements are possible.
\dQuote{\link{...} one can only reject the null hypothesis if the test
statistics falls into the critical region(s), or fail to reject this
hypothesis. In the latter case, all we can say is that no significant effect
was observed, but one cannot conclude that the null hypothesis is true.}
(\cite{Pernet 2017}). One way to address this issues without Bayesian methods
is \emph{Equivalence Testing}, as implemented in \code{equivalence_test()}.
While you either can reject the null hypothesis or claim an inconclusive result
in NHST, the equivalence test adds a third category, \emph{"accept"}. Roughly
speaking, the idea behind equivalence testing in a frequentist framework is
to check whether an estimate and its uncertainty (i.e. confidence interval)
falls within a region of "practical equivalence". Depending on the rule for
this test (see below), statistical significance does not necessarily indicate
whether the null hypothesis can be rejected or not, i.e. the classical
interpretation of the p-value may differ from the results returned from
the equivalence test.

\subsection{Calculation of equivalence testing}{
\describe{
\item{"bayes" - Bayesian rule (Kruschke 2018)}{
This rule follows the \dQuote{HDI+ROPE decision rule} \cite{(Kruschke,
2014, 2018)} used for the
\code{\link[bayestestR:equivalence_test]{Bayesian counterpart()}}. This
means, if the confidence intervals are completely outside the ROPE, the
"null hypothesis" for this parameter is "rejected". If the ROPE
completely covers the CI, the null hypothesis is accepted. Else, it's
undecided whether to accept or reject the null hypothesis. Desirable
results are low proportions inside the ROPE (the closer to zero the
better).
}
\item{"classic" - The TOST rule (Lakens 2017)}{
This rule follows the \dQuote{TOST rule}, i.e. a two one-sided test
procedure (\cite{Lakens 2017}). Following this rule, practical
equivalence of an effect (i.e. H0) is \emph{rejected}, when the
coefficient is statistically significant \emph{and} the narrow
confidence intervals (i.e. \code{1-2*alpha}) \emph{include} or
\emph{exceed} the ROPE boundaries. Practical equivalence is assumed
(i.e. H0 accepted) when the narrow confidence intervals are completely
inside the ROPE, no matter if the effect is statistically significant
or not. Else, the decision whether to accept or reject H0 is undecided.
}
\item{"cet" - Conditional Equivalence Testing (Campbell/Gustafson 2018)}{
The Conditional Equivalence Testing as described by \cite{Campbell and
Gustafson 2018}. According to this rule, practical equivalence is
rejected when the coefficient is statistically significant. When the
effect is \emph{not} significant and the narrow confidence intervals are
completely inside the ROPE, we accept H0, else it is undecided.
}
}
}
\subsection{Levels of Confidence Intervals used for Equivalence Testing}{
For \code{rule = "classic"}, "narrow" confidence intervals are used for
equivalence testing. "Narrow" means, the the intervals is not 1 - alpha,
but 1 - 2 * alpha. Thus, if \code{ci = .95}, alpha is assumed to be 0.05
and internally a ci-level of 0.90 is used. \code{rule = "cet"} uses
both regular and narrow confidence intervals, while \code{rule = "bayes"}
only uses the regular intervals.
}
\subsection{p-Values}{
The equivalence p-value is the area of the (cumulative) confidence
distribution that is outside of the region of equivalence. It can be
interpreted as p-value for \emph{rejecting} the alternative hypothesis
and \emph{accepting} the null hypothesis.
}
\subsection{Second Generation p-Value (SGPV)}{
Second generation p-values (SGPV) were proposed as a statistic
that represents \dQuote{the proportion of data-supported hypotheses
that are also null hypotheses} \cite{(Blume et al. 2018)}. This statistic
is actually computed in the same way as the percentage inside the ROPE as
returned by \code{equivalence_test()} (see \cite{Lakens and Delacre 2020}
for details on computation of the SGPV). Thus, the \code{"inside ROPE"}
column reflects the SGPV.
}
\subsection{ROPE range}{
Some attention is required for finding suitable values for the ROPE limits
(argument \code{range}). See 'Details' in \code{\link[bayestestR:rope_range]{bayestestR::rope_range()}}
for further information.
}
}
\note{
There is also a \href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
data(qol_cancer)
model <- lm(QoL ~ time + age + education, data = qol_cancer)

# default rule
equivalence_test(model)

# conditional equivalence test
equivalence_test(model, rule = "cet")

# plot method
if (require("see", quietly = TRUE)) {
  result <- equivalence_test(model)
  plot(result)
}
}
\references{
\itemize{
\item Blume, J. D., D'Agostino McGowan, L., Dupont, W. D., & Greevy, R. A.
(2018). Second-generation p-values: Improved rigor, reproducibility, &
transparency in statistical analyses. PLOS ONE, 13(3), e0188299.
https://doi.org/10.1371/journal.pone.0188299

\item Campbell, H., & Gustafson, P. (2018). Conditional equivalence
testing: An alternative remedy for publication bias. PLOS ONE, 13(4),
e0195145. doi: 10.1371/journal.pone.0195145

\item Kruschke, J. K. (2014). Doing Bayesian data analysis: A tutorial with
R, JAGS, and Stan. Academic Press

\item Kruschke, J. K. (2018). Rejecting or accepting parameter values in
Bayesian estimation. Advances in Methods and Practices in Psychological
Science, 1(2), 270-280. doi: 10.1177/2515245918771304

\item Lakens, D. (2017). Equivalence Tests: A Practical Primer for t Tests,
Correlations, and Meta-Analyses. Social Psychological and Personality
Science, 8(4), 355–362. doi: 10.1177/1948550617697177

\item Lakens, D., & Delacre, M. (2020). Equivalence Testing and the Second
Generation P-Value. Meta-Psychology, 4.
https://doi.org/10.15626/MP.2018.933

\item Pernet, C. (2017). Null hypothesis significance testing: A guide to
commonly misunderstood concepts and recommendations for good practice.
F1000Research, 4, 621. doi: 10.12688/f1000research.6963.5
}
}
\seealso{
For more details, see \code{\link[bayestestR:equivalence_test]{bayestestR::equivalence_test()}}.
Further readings can be found in the references.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_analysis.R
\name{cluster_analysis}
\alias{cluster_analysis}
\title{Cluster Analysis}
\usage{
cluster_analysis(
  x,
  n = NULL,
  method = "kmeans",
  include_factors = FALSE,
  standardize = TRUE,
  verbose = TRUE,
  distance_method = "euclidean",
  hclust_method = "complete",
  kmeans_method = "Hartigan-Wong",
  dbscan_eps = 15,
  iterations = 100,
  ...
)
}
\arguments{
\item{x}{A data frame.}

\item{n}{Number of clusters used for supervised cluster methods. If \code{NULL},
the number of clusters to extract is determined by calling \code{\link[=n_clusters]{n_clusters()}}. Note
that this argument does not apply for unsupervised clustering methods like
\code{dbscan}, \code{mixture}, \code{pvclust}, or \code{pamk}.}

\item{method}{Method for computing the cluster analysis. Can be \code{"kmeans"}
(default; k-means using \code{kmeans()}), \code{"hkmeans"} (hierarchical k-means
using \code{factoextra::hkmeans()}), \code{pam} (K-Medoids using \code{cluster::pam()}),
\code{pamk} (K-Medoids that finds out the number of clusters), \code{"hclust"}
(hierarchical clustering using \code{hclust()} or \code{pvclust::pvclust()}),
\code{dbscan} (DBSCAN using \code{dbscan::dbscan()}), \code{hdbscan} (Hierarchical DBSCAN
using \code{dbscan::hdbscan()}), or \code{mixture} (Mixture modelling using
\code{mclust::Mclust()}, which requires the user to run \code{library(mclust)}
before).}

\item{include_factors}{Logical, if \code{TRUE}, factors are converted to numerical
values in order to be included in the data for determining the number of
clusters. By default, factors are removed, because most methods that
determine the number of clusters need numeric input only.}

\item{standardize}{Standardize the dataframe before clustering (default).}

\item{verbose}{Toggle warnings and messages.}

\item{distance_method}{Distance measure to be used for methods based on
distances (e.g., when \code{method = "hclust"} for hierarchical clustering. For
other methods, such as \code{"kmeans"}, this argument will be ignored). Must be
one of \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"}
or \code{"minkowski"}. See \code{\link[=dist]{dist()}} and \code{pvclust::pvclust()} for more
information.}

\item{hclust_method}{Agglomeration method to be used when \code{method = "hclust"}
or \code{method = "hkmeans"} (for hierarchical clustering). This should be one
of \code{"ward"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
\code{"mcquitty"}, \code{"median"} or \code{"centroid"}. Default is \code{"complete"} (see
\code{\link[=hclust]{hclust()}}).}

\item{kmeans_method}{Algorithm used for calculating kmeans cluster. Only applies,
if \code{method = "kmeans"}. May be one of \code{"Hartigan-Wong"} (default),
\code{"Lloyd"} (used by SPSS), or \code{"MacQueen"}. See \code{\link[=kmeans]{kmeans()}} for details on
this argument.}

\item{dbscan_eps}{The 'eps' argument for DBSCAN method. See \code{\link[=n_clusters_dbscan]{n_clusters_dbscan()}}.}

\item{iterations}{The number of replications.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
The group classification for each observation as vector. The
returned vector includes missing values, so it has the same length
as \code{nrow(x)}.
}
\description{
Compute hierarchical or kmeans cluster analysis and return the group
assignment for each observation as vector.
}
\details{
The \code{print()} and \code{plot()} methods show the (standardized) mean value for
each variable within each cluster. Thus, a higher absolute value indicates
that a certain variable characteristic is more pronounced within that
specific cluster (as compared to other cluster groups with lower absolute
mean values).
}
\note{
There is also a \href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
set.seed(33)
# K-Means ====================================================
rez <- cluster_analysis(iris[1:4], n = 3, method = "kmeans")
rez # Show results
predict(rez) # Get clusters
summary(rez) # Extract the centers values (can use 'plot()' on that)
cluster_discrimination(rez) # Perform LDA

# Hierarchical k-means (more robust k-means)
if (require("factoextra", quietly = TRUE)) {
  rez <- cluster_analysis(iris[1:4], n = 3, method = "hkmeans")
  rez # Show results
  predict(rez) # Get clusters
}

# Hierarchical Clustering (hclust) ===========================
rez <- cluster_analysis(iris[1:4], n = 3, method = "hclust")
rez # Show results
predict(rez) # Get clusters

# K-Medoids (pam) ============================================
if (require("cluster", quietly = TRUE)) {
  rez <- cluster_analysis(iris[1:4], n = 3, method = "pam")
  rez # Show results
  predict(rez) # Get clusters
}

# PAM with automated number of clusters
if (require("fpc", quietly = TRUE)) {
  rez <- cluster_analysis(iris[1:4], method = "pamk")
  rez # Show results
  predict(rez) # Get clusters
}

# DBSCAN ====================================================
if (require("dbscan", quietly = TRUE)) {
  # Note that you can assimilate more outliers (cluster 0) to neighbouring
  # clusters by setting borderPoints = TRUE.
  rez <- cluster_analysis(iris[1:4], method = "dbscan", dbscan_eps = 1.45)
  rez # Show results
  predict(rez) # Get clusters
}

# Mixture ====================================================
if (require("mclust", quietly = TRUE)) {
  library(mclust) # Needs the package to be loaded
  rez <- cluster_analysis(iris[1:4], method = "mixture")
  rez # Show results
  predict(rez) # Get clusters
}
}
\references{
\itemize{
\item Maechler M, Rousseeuw P, Struyf A, Hubert M, Hornik K (2014) cluster: Cluster
Analysis Basics and Extensions. R package.
}
}
\seealso{
\itemize{
\item \code{\link[=n_clusters]{n_clusters()}} to determine the number of clusters to extract,
\code{\link[=cluster_discrimination]{cluster_discrimination()}} to determine the accuracy of cluster group
classification via linear discriminant analysis (LDA) and
\code{\link[=check_clusterstructure]{check_clusterstructure()}} to check suitability of data
for clustering.
\item https://www.datanovia.com/en/lessons/
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_factors.R
\name{n_factors}
\alias{n_factors}
\alias{n_components}
\title{Number of components/factors to retain in PCA/FA}
\usage{
n_factors(
  x,
  type = "FA",
  rotation = "varimax",
  algorithm = "default",
  package = c("nFactors", "psych"),
  cor = NULL,
  safe = TRUE,
  n_max = NULL,
  ...
)

n_components(
  x,
  type = "PCA",
  rotation = "varimax",
  algorithm = "default",
  package = c("nFactors", "psych"),
  cor = NULL,
  safe = TRUE,
  ...
)
}
\arguments{
\item{x}{A data frame.}

\item{type}{Can be \code{"FA"} or \code{"PCA"}, depending on what you want to
do.}

\item{rotation}{Only used for VSS (Very Simple Structure criterion, see
\code{\link[psych:VSS]{psych::VSS()}}). The rotation to apply. Can be \code{"none"},
\code{"varimax"}, \code{"quartimax"}, \code{"bentlerT"}, \code{"equamax"},
\code{"varimin"}, \code{"geominT"} and \code{"bifactor"} for orthogonal
rotations, and \code{"promax"}, \code{"oblimin"}, \code{"simplimax"},
\code{"bentlerQ"}, \code{"geominQ"}, \code{"biquartimin"} and
\code{"cluster"} for oblique transformations.}

\item{algorithm}{Factoring method used by VSS. Can be \code{"pa"} for
Principal Axis Factor Analysis, \code{"minres"} for minimum residual (OLS)
factoring, \code{"mle"} for Maximum Likelihood FA and \code{"pc"} for
Principal Components. \code{"default"} will select \code{"minres"} if
\code{type = "FA"} and \code{"pc"} if \code{type = "PCA"}.}

\item{package}{Package from which respective methods are used. Can be
\code{"all"} or a vector containing \code{"nFactors"}, \code{"psych"}, \code{"PCDimension"}, \code{"fit"} or
\code{"EGAnet"}. Note that \code{"fit"} (which actually also relies on the \code{psych}
package) and \code{"EGAnet"} can be very slow for bigger
datasets. Thus, the default is \code{c("nFactors", "psych")}. You must have
the respective packages installed for the methods to be used.}

\item{cor}{An optional correlation matrix that can be used (note that the
data must still be passed as the first argument). If \code{NULL}, will
compute it by running \code{cor()} on the passed data.}

\item{safe}{If \code{TRUE}, the function will run all the procedures in try
blocks, and will only return those that work and silently skip the ones
that may fail.}

\item{n_max}{If set to a value (e.g., \code{10}), will drop from the results all
methods that suggest a higher number of components. The interpretation becomes
'from all the methods that suggested a number lower than n_max, the results
are ...'.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame.
}
\description{
This function runs many existing procedures for determining how many factors
to retain/extract from factor analysis (FA) or dimension reduction (PCA). It
returns the number of factors based on the maximum consensus between methods.
In case of ties, it will keep the simplest model and select the solution
with the fewer factors.
}
\details{
\code{n_components} is actually an alias for \code{n_factors}, with
different defaults for the function arguments.
}
\note{
There is also a
\href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method}
implemented in the
\href{https://easystats.github.io/see/}{\pkg{see}-package}.
\code{n_components()} is a convenient short for \code{n_factors(type = "PCA")}.
}
\examples{
library(parameters)
if (require("nFactors", quietly = TRUE) && require("EGAnet", quietly = TRUE)) {
  n_factors(mtcars, type = "PCA")

  result <- n_factors(mtcars[1:5], type = "FA")
  as.data.frame(result)
  summary(result)
  \dontrun{
  if (require("PCDimension", quietly = TRUE)) {
    # Setting package = 'all' will increase the number of methods (but is slow)
    n_factors(mtcars, type = "PCA", package = "all")
    n_factors(mtcars, type = "FA", algorithm = "mle", package = "all")
  }
  }
}
}
\references{
\itemize{
\item Bartlett, M. S. (1950). Tests of significance in factor analysis.
British Journal of statistical psychology, 3(2), 77-85.

\item Bentler, P. M., & Yuan, K. H. (1996). Test of linear trend in
eigenvalues of a covariance matrix with application to data analysis.
British Journal of Mathematical and Statistical Psychology, 49(2), 299-312.

\item Cattell, R. B. (1966). The scree test for the number of factors.
Multivariate behavioral research, 1(2), 245-276.

\item Finch, W. H. (2019). Using Fit Statistic Differences to Determine the
Optimal Number of Factors to Retain in an Exploratory Factor Analysis.
Educational and Psychological Measurement.

\item Zoski, K. W., & Jurs, S. (1996). An objective counterpart to the
visual scree test for factor analysis: The standard error scree.
Educational and Psychological Measurement, 56(3), 443-451.

\item Zoski, K., & Jurs, S. (1993). Using multiple regression to determine
the number of factors to retain in factor analysis. Multiple Linear
Regression Viewpoints, 20(1), 5-9.

\item Nasser, F., Benson, J., & Wisenbaker, J. (2002). The performance of
regression-based variations of the visual scree for determining the number
of common factors. Educational and psychological measurement, 62(3),
397-419.

\item Golino, H., Shi, D., Garrido, L. E., Christensen, A. P., Nieto, M.
D., Sadana, R., & Thiyagarajan, J. A. (2018). Investigating the performance
of Exploratory Graph Analysis and traditional techniques to identify the
number of latent factors: A simulation and tutorial.

\item Golino, H. F., & Epskamp, S. (2017). Exploratory graph analysis: A
new approach for estimating the number of dimensions in psychological
research. PloS one, 12(6), e0174035.

\item Revelle, W., & Rocklin, T. (1979). Very simple structure: An
alternative procedure for estimating the optimal number of interpretable
factors. Multivariate Behavioral Research, 14(4), 403-414.

\item Velicer, W. F. (1976). Determining the number of components from the
matrix of partial correlations. Psychometrika, 41(3), 321-327.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_kenward.R, R/dof_kenward.R,
%   R/p_value_kenward.R, R/standard_error_kenward.R
\name{ci_kenward}
\alias{ci_kenward}
\alias{dof_kenward}
\alias{p_value_kenward}
\alias{se_kenward}
\title{Kenward-Roger approximation for SEs, CIs and p-values}
\usage{
ci_kenward(model, ci = 0.95)

dof_kenward(model)

p_value_kenward(model, dof = NULL)

se_kenward(model)
}
\arguments{
\item{model}{A statistical model.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{dof}{Degrees of Freedom.}
}
\value{
A data frame.
}
\description{
An approximate F-test based on the Kenward-Roger (1997) approach.
}
\details{
Inferential statistics (like p-values, confidence intervals and
standard errors) may be biased in mixed models when the number of clusters
is small (even if the sample size of level-1 units is high). In such cases
it is recommended to approximate a more accurate number of degrees of freedom
for such inferential statistics. Unlike simpler approximation heuristics
like the "m-l-1" rule (\code{dof_ml1}), the Kenward-Roger approximation is
also applicable in more complex multilevel designs, e.g. with cross-classified
clusters. However, the "m-l-1" heuristic also applies to generalized
mixed models, while approaches like Kenward-Roger or Satterthwaite are limited
to linear mixed models only.
}
\examples{
\donttest{
if (require("lme4", quietly = TRUE)) {
  model <- lmer(Petal.Length ~ Sepal.Length + (1 | Species), data = iris)
  p_value_kenward(model)
}
}
}
\references{
Kenward, M. G., & Roger, J. H. (1997). Small sample inference for
fixed effects from restricted maximum likelihood. Biometrics, 983-997.
}
\seealso{
\code{dof_kenward()} and \code{se_kenward()} are small helper-functions
to calculate approximated degrees of freedom and standard errors for model
parameters, based on the Kenward-Roger (1997) approach.
\cr \cr
\code{\link[=dof_satterthwaite]{dof_satterthwaite()}} and
\code{\link[=dof_ml1]{dof_ml1()}} approximate degrees
of freedom based on Satterthwaite's method or the "m-l-1" rule.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_factorstructure.R
\name{check_sphericity_bartlett}
\alias{check_sphericity_bartlett}
\title{Bartlett's Test of Sphericity}
\usage{
check_sphericity_bartlett(x, ...)
}
\arguments{
\item{x}{A dataframe.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A list of indices related to sphericity.
}
\description{
Bartlett's (1951) test of sphericity tests whether a matrix (of correlations)
is significantly different from an identity matrix. The test provides
probability that the correlation matrix has significant correlations among at
least some of the variables in a dataset, a prerequisite for factor analysis
to work. In other words, before starting with factor analysis, one needs to
check whether Bartlett’s test of sphericity is significant.
}
\details{
This function is strongly inspired by the \code{cortest.bartlett}
function in the \pkg{psych} package (Revelle, 2016). All credit goes to its
author.
}
\examples{
library(parameters)
check_sphericity_bartlett(mtcars)
}
\references{
\itemize{
\item Revelle, W. (2016). How To: Use the psych package for Factor Analysis
and data reduction.
\item Bartlett, M. S. (1951). The effect of standardization on a Chi-square
approximation in factor analysis. Biometrika, 38(3/4), 337-344.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_htest.R
\name{model_parameters.htest}
\alias{model_parameters.htest}
\alias{model_parameters.pairwise.htest}
\title{Parameters from hypothesis tests}
\usage{
\method{model_parameters}{htest}(
  model,
  cramers_v = NULL,
  phi = NULL,
  standardized_d = NULL,
  hedges_g = NULL,
  omega_squared = NULL,
  eta_squared = NULL,
  epsilon_squared = NULL,
  cohens_g = NULL,
  rank_biserial = NULL,
  rank_epsilon_squared = NULL,
  kendalls_w = NULL,
  ci = 0.95,
  alternative = NULL,
  bootstrap = FALSE,
  verbose = TRUE,
  ...
)

\method{model_parameters}{pairwise.htest}(model, verbose = TRUE, ...)
}
\arguments{
\item{model}{Object of class \code{htest} or \code{pairwise.htest}.}

\item{cramers_v, phi}{Compute Cramer's V or phi as index of effect size.
Can be \code{"raw"} or \code{"adjusted"} (effect size will be bias-corrected).
Only applies to objects from \code{chisq.test()}.}

\item{standardized_d}{If \code{TRUE}, compute standardized d as index of
effect size. Only applies to objects from \code{t.test()}. Calculation of
\code{d} is based on the t-value (see \code{\link[effectsize:t_to_r]{effectsize::t_to_d()}})
for details.}

\item{hedges_g}{If \code{TRUE}, compute Hedge's g as index of effect size.
Only applies to objects from \code{t.test()}.}

\item{omega_squared, eta_squared, epsilon_squared}{Logical, if \code{TRUE},
returns the non-partial effect size Omega, Eta or Epsilon squared. Only
applies to objects from \code{oneway.test()}.}

\item{cohens_g}{If \code{TRUE}, compute Cohen's g as index of effect size.
Only applies to objects from \code{mcnemar.test()}.}

\item{rank_biserial}{If \code{TRUE}, compute the rank-biserial correlation as
effect size measure. Only applies to objects from \code{wilcox.test()}.}

\item{rank_epsilon_squared}{If \code{TRUE}, compute the rank epsilon squared
as effect size measure. Only applies to objects from \code{kruskal.test()}.}

\item{kendalls_w}{If \code{TRUE}, compute the Kendall's coefficient of
concordance as effect size measure. Only applies to objects from
\code{friedman.test()}.}

\item{ci}{Level of confidence intervals for effect size statistic. Currently
only applies to objects from \code{chisq.test()} or \code{oneway.test()}.}

\item{alternative}{A character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"two.sided"} (default, two-sided CI),
\code{"greater"} or \code{"less"} (one-sided CI). Partial matching is allowed
(e.g., \code{"g"}, \code{"l"}, \code{"two"}...). See section \emph{One-Sided CIs} in
the \href{https://easystats.github.io/effectsize/}{effectsize_CIs vignette}.}

\item{bootstrap}{Should estimates be bootstrapped?}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters of h-tests (correlations, t-tests, chi-squared, ...).
}
\examples{
model <- cor.test(mtcars$mpg, mtcars$cyl, method = "pearson")
model_parameters(model)

model <- t.test(iris$Sepal.Width, iris$Sepal.Length)
model_parameters(model)

model <- t.test(mtcars$mpg ~ mtcars$vs)
model_parameters(model)

model <- t.test(iris$Sepal.Width, mu = 1)
model_parameters(model)

data(airquality)
airquality$Month <- factor(airquality$Month, labels = month.abb[5:9])
model <- pairwise.t.test(airquality$Ozone, airquality$Month)
model_parameters(model)

smokers <- c(83, 90, 129, 70)
patients <- c(86, 93, 136, 82)
model <- pairwise.prop.test(smokers, patients)
model_parameters(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.data_frame}
\alias{.data_frame}
\title{help-functions}
\usage{
.data_frame(...)
}
\description{
help-functions
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.factor_to_dummy}
\alias{.factor_to_dummy}
\title{Safe transformation from factor/character to numeric}
\usage{
.factor_to_dummy(x)
}
\description{
Safe transformation from factor/character to numeric
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_standard_error.R, R/methods_DirichletReg.R,
%   R/methods_averaging.R, R/methods_base.R, R/methods_betareg.R,
%   R/methods_glmmTMB.R, R/methods_lme4.R, R/methods_mfx.R, R/methods_mixmod.R,
%   R/methods_mixor.R, R/methods_ordinal.R, R/methods_pscl.R,
%   R/methods_survival.R
\name{standard_error}
\alias{standard_error}
\alias{standard_error.default}
\alias{standard_error.DirichletRegModel}
\alias{standard_error.averaging}
\alias{standard_error.factor}
\alias{standard_error.betareg}
\alias{standard_error.glmmTMB}
\alias{standard_error.merMod}
\alias{standard_error.poissonmfx}
\alias{standard_error.betamfx}
\alias{standard_error.MixMod}
\alias{standard_error.mixor}
\alias{standard_error.clm2}
\alias{standard_error.zeroinfl}
\alias{standard_error.coxph}
\title{Standard Errors}
\usage{
standard_error(model, ...)

\method{standard_error}{default}(model, method = NULL, verbose = TRUE, ...)

\method{standard_error}{DirichletRegModel}(model, component = c("all", "conditional", "precision"), ...)

\method{standard_error}{averaging}(model, component = c("conditional", "full"), ...)

\method{standard_error}{factor}(model, force = FALSE, verbose = TRUE, ...)

\method{standard_error}{betareg}(model, component = c("all", "conditional", "precision"), ...)

\method{standard_error}{glmmTMB}(
  model,
  effects = c("fixed", "random"),
  component = c("all", "conditional", "zi", "zero_inflated", "dispersion"),
  verbose = TRUE,
  ...
)

\method{standard_error}{merMod}(model, effects = c("fixed", "random"), method = NULL, ...)

\method{standard_error}{poissonmfx}(model, component = c("all", "conditional", "marginal"), ...)

\method{standard_error}{betamfx}(
  model,
  component = c("all", "conditional", "precision", "marginal"),
  ...
)

\method{standard_error}{MixMod}(
  model,
  effects = c("fixed", "random"),
  component = c("all", "conditional", "zi", "zero_inflated"),
  robust = FALSE,
  verbose = TRUE,
  ...
)

\method{standard_error}{mixor}(model, effects = "all", ...)

\method{standard_error}{clm2}(model, component = c("all", "conditional", "scale"), ...)

\method{standard_error}{zeroinfl}(
  model,
  component = c("all", "conditional", "zi", "zero_inflated"),
  method = NULL,
  verbose = TRUE,
  ...
)

\method{standard_error}{coxph}(model, method = NULL, ...)
}
\arguments{
\item{model}{A model.}

\item{...}{Arguments passed to or from other methods. For
\code{standard_error()}, if \code{method = "robust"}, arguments
\code{vcov_estimation}, \code{vcov_type} and \code{vcov_args} can be passed
down to \code{\link[=standard_error_robust]{standard_error_robust()}}.}

\item{method}{If \code{"robust"}, robust standard errors are computed by
calling \code{\link[=standard_error_robust]{standard_error_robust()}}.
\code{standard_error_robust()}, in turn, calls one of the
\verb{vcov*()}-functions from the \pkg{sandwich} or \pkg{clubSandwich}
package for robust covariance matrix estimators. For linear mixed models,
\code{method} may also be \code{\link[=p_value_kenward]{"kenward"}} or
\code{\link[=p_value_satterthwaite]{"satterthwaite"}}.}

\item{verbose}{Toggle warnings and messages.}

\item{component}{Should all parameters, parameters for the conditional model,
or for the zero-inflated part of the model be returned? Applies to models
with zero-inflated component. \code{component} may be one of \code{"conditional"},
\code{"zi"}, \code{"zero-inflated"}, \code{"dispersion"} or \code{"all"}
(default). May be abbreviated.}

\item{force}{Logical, if \code{TRUE}, factors are converted to numerical
values to calculate the standard error, with the lowest level being the
value \code{1} (unless the factor has numeric levels, which are converted
to the corresponding numeric value). By default, \code{NA} is returned for
factors or character vectors.}

\item{effects}{Should standard errors for fixed effects or random effects be
returned? Only applies to mixed models. May be abbreviated. When standard
errors for random effects are requested, for each grouping factor a list of
standard errors (per group level) for random intercepts and slopes is
returned.}

\item{robust}{Logical, if \code{TRUE}, computes confidence intervals (or p-values)
based on robust standard errors. See \code{\link[=standard_error_robust]{standard_error_robust()}}.}
}
\value{
A data frame with at least two columns: the parameter names and the
standard errors. Depending on the model, may also include columns for model
components etc.
}
\description{
\code{standard_error()} attempts to return standard errors of model
parameters, while \code{standard_error_robust()} attempts to return robust
standard errors.
}
\note{
For Bayesian models (from \pkg{rstanarm} or \pkg{brms}), the standard
error is the SD of the posterior samples.
}
\examples{
model <- lm(Petal.Length ~ Sepal.Length * Species, data = iris)
standard_error(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_factors.R
\name{.n_factors_cng}
\alias{.n_factors_cng}
\title{Cattell-Nelson-Gorsuch CNG Indices}
\usage{
.n_factors_cng(eigen_values = NULL, model = "factors")
}
\description{
Cattell-Nelson-Gorsuch CNG Indices
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_scores.R
\name{get_scores}
\alias{get_scores}
\title{Get Scores from Principal Component Analysis (PCA)}
\usage{
get_scores(x, n_items = NULL)
}
\arguments{
\item{x}{An object returned by \code{\link[=principal_components]{principal_components()}}.}

\item{n_items}{Number of required (i.e. non-missing) items to build the sum
score. If \code{NULL}, the value is chosen to match half of the number of
columns in a data frame.}
}
\value{
A data frame with subscales, which are average sum scores for all
items from each component.
}
\description{
\code{get_scores()} takes \code{n_items} amount of items that load the most
(either by loading cutoff or number) on a component, and then computes their
average.
}
\details{
\code{get_scores()} takes the results from
\code{\link[=principal_components]{principal_components()}} and extracts the variables for each
component found by the PCA. Then, for each of these "subscales", row means
are calculated (which equals adding up the single items and dividing by the
number of items). This results in a sum score for each component from the
PCA, which is on the same scale as the original, single items that were
used to compute the PCA.
}
\examples{
if (require("psych")) {
  pca <- principal_components(mtcars[, 1:7], n = 2, rotation = "varimax")

  # PCA extracted two components
  pca

  # assignment of items to each component
  closest_component(pca)

  # now we want to have sum scores for each component
  get_scores(pca)

  # compare to manually computed sum score for 2nd component, which
  # consists of items "hp" and "qsec"
  (mtcars$hp + mtcars$qsec) / 2
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.compact_character}
\alias{.compact_character}
\title{remove empty string from character}
\usage{
.compact_character(x)
}
\description{
remove empty string from character
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_metafor.R
\name{model_parameters.rma}
\alias{model_parameters.rma}
\title{Parameters from Meta-Analysis}
\usage{
\method{model_parameters}{rma}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  include_studies = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{Model object.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{bootstrap}{Should estimates be based on bootstrapped model? If
\code{TRUE}, then arguments of \link[=model_parameters.stanreg]{Bayesian regressions} apply (see also
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}).}

\item{iterations}{The number of bootstrap replicates. This only apply in the
case of bootstrapped frequentist models.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{include_studies}{Logical, if \code{TRUE} (default), includes parameters
for all studies. Else, only parameters for overall-effects are shown.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods. For instance, when
\code{bootstrap = TRUE}, arguments like \code{type} or \code{parallel} are
passed down to \code{bootstrap_model()}, and arguments like \code{ci_method}
are passed down to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Extract and compute indices and measures to describe parameters of meta-analysis models.
}
\examples{
library(parameters)
mydat <<- data.frame(
  effectsize = c(-0.393, 0.675, 0.282, -1.398),
  stderr = c(0.317, 0.317, 0.13, 0.36)
)
if (require("metafor", quietly = TRUE)) {
  model <- rma(yi = effectsize, sei = stderr, method = "REML", data = mydat)
  model_parameters(model)
}
\dontrun{
# with subgroups
if (require("metafor", quietly = TRUE)) {
  data(dat.bcg)
  dat <- escalc(
    measure = "RR",
    ai = tpos,
    bi = tneg,
    ci = cpos,
    di = cneg,
    data = dat.bcg
  )
  dat$alloc <- ifelse(dat$alloc == "random", "random", "other")
  model <- rma(yi, vi, mods = ~alloc, data = dat, digits = 3, slab = author)
  model_parameters(model)
}

if (require("metaBMA", quietly = TRUE)) {
  data(towels)
  m <- meta_random(logOR, SE, study, data = towels)
  model_parameters(m)
}
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_base.R, R/methods_brms.R,
%   R/methods_rstanarm.R
\name{model_parameters.data.frame}
\alias{model_parameters.data.frame}
\alias{model_parameters.brmsfit}
\alias{model_parameters.stanreg}
\title{Parameters from Bayesian Models}
\usage{
\method{model_parameters}{data.frame}(
  model,
  centrality = "median",
  dispersion = FALSE,
  ci = 0.95,
  ci_method = "hdi",
  test = c("pd", "rope"),
  rope_range = "default",
  rope_ci = 0.95,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)

\method{model_parameters}{brmsfit}(
  model,
  centrality = "median",
  dispersion = FALSE,
  ci = 0.95,
  ci_method = "hdi",
  test = c("pd", "rope"),
  rope_range = "default",
  rope_ci = 0.95,
  bf_prior = NULL,
  diagnostic = c("ESS", "Rhat"),
  priors = FALSE,
  effects = "fixed",
  component = "all",
  exponentiate = FALSE,
  standardize = NULL,
  group_level = FALSE,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)

\method{model_parameters}{stanreg}(
  model,
  centrality = "median",
  dispersion = FALSE,
  ci = 0.95,
  ci_method = "hdi",
  test = c("pd", "rope"),
  rope_range = "default",
  rope_ci = 0.95,
  bf_prior = NULL,
  diagnostic = c("ESS", "Rhat"),
  priors = TRUE,
  effects = "fixed",
  exponentiate = FALSE,
  standardize = NULL,
  group_level = FALSE,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{Bayesian model (including SEM from \pkg{blavaan}. May also be
a data frame with posterior samples.}

\item{centrality}{The point-estimates (centrality indices) to compute.  Character (vector) or list with one or more of these options: \code{"median"}, \code{"mean"}, \code{"MAP"} or \code{"all"}.}

\item{dispersion}{Logical, if \code{TRUE}, computes indices of dispersion related to the estimate(s) (\code{SD} and \code{MAD} for \code{mean} and \code{median}, respectively).}

\item{ci}{Credible Interval (CI) level. Default to \code{0.95} (\verb{95\%}). See
\code{\link[bayestestR:ci]{bayestestR::ci()}} for further details.}

\item{ci_method}{Method for computing degrees of freedom for
confidence intervals (CI) and the related p-values. Allowed are following
options (which vary depending on the model class): \code{"residual"},
\code{"normal"}, \code{"likelihood"}, \code{"satterthwaite"}, \code{"kenward"}, \code{"wald"},
\code{"profile"}, \code{"boot"}, \code{"uniroot"}, \code{"ml1"}, \code{"betwithin"}, \code{"hdi"},
\code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, or \code{"bcai"}. See section
\emph{Confidence intervals and approximation of degrees of freedom} in
\code{\link[=model_parameters]{model_parameters()}} for further details. When \code{ci_method=NULL}, in most
cases \code{"wald"} is used then.}

\item{test}{The indices of effect existence to compute. Character (vector) or
list with one or more of these options: \code{"p_direction"} (or \code{"pd"}),
\code{"rope"}, \code{"p_map"}, \code{"equivalence_test"} (or \code{"equitest"}),
\code{"bayesfactor"} (or \code{"bf"}) or \code{"all"} to compute all tests.
For each "test", the corresponding \pkg{bayestestR} function is called
(e.g. \code{\link[bayestestR:rope]{rope()}} or \code{\link[bayestestR:p_direction]{p_direction()}}) and its results
included in the summary output.}

\item{rope_range}{ROPE's lower and higher bounds. Should be a list of two
values (e.g., \code{c(-0.1, 0.1)}) or \code{"default"}. If \code{"default"},
the bounds are set to \code{x +- 0.1*SD(response)}.}

\item{rope_ci}{The Credible Interval (CI) probability, corresponding to the
proportion of HDI, to use for the percentage in ROPE.}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}

\item{verbose}{Toggle messages and warnings.}

\item{...}{Currently not used.}

\item{bf_prior}{Distribution representing a prior for the computation of
Bayes factors / SI. Used if the input is a posterior, otherwise (in the
case of models) ignored.}

\item{diagnostic}{Diagnostic metrics to compute.  Character (vector) or list
with one or more of these options: \code{"ESS"}, \code{"Rhat"}, \code{"MCSE"} or \code{"all"}.}

\item{priors}{Add the prior used for each parameter.}

\item{effects}{Should results for fixed effects, random effects or both be
returned? Only applies to mixed models. May be abbreviated.}

\item{component}{Which type of parameters to return, such as parameters for the
conditional model, the zero-inflated part of the model, the dispersion
term, or other auxiliary parameters be returned? Applies to models with
zero-inflated and/or dispersion formula, or if parameters such as \code{sigma}
should be included. May be abbreviated. Note that the \emph{conditional}
component is also called \emph{count} or \emph{mean} component, depending on the
model. There are three convenient shortcuts: \code{component = "all"} returns
all possible parameters. If \code{component = "location"}, location parameters
such as \code{conditional}, \code{zero_inflated}, or \code{smooth_terms}, are returned
(everything that are fixed or random effects - depending on the \code{effects}
argument - but no auxiliary parameters). For \code{component = "distributional"}
(or \code{"auxiliary"}), components like \code{sigma}, \code{dispersion}, or \code{beta}
(and other auxiliary parameters) are returned.}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{group_level}{Logical, for multilevel models (i.e. models with random
effects) and when \code{effects = "all"} or \code{effects = "random"},
include the parameters for each group level from random effects. If
\code{group_level = FALSE} (the default), only information on SD and COR
are shown.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from Bayesian models.
}
\note{
When \code{standardize = "refit"}, columns \code{diagnostic},
\code{bf_prior} and \code{priors} refer to the \emph{original}
\code{model}. If \code{model} is a data frame, arguments \code{diagnostic},
\code{bf_prior} and \code{priors} are ignored. \cr \cr There is also a
\href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method}
implemented in the
\href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\section{Confidence intervals and approximation of degrees of freedom}{

There are different ways of approximating the degrees of freedom depending
on different assumptions about the nature of the model and its sampling
distribution. The \code{ci_method} argument modulates the method for computing degrees
of freedom (df) that are used to calculate confidence intervals (CI) and the
related p-values. Following options are allowed, depending on the model
class:

\strong{Classical methods:}

Classical inference is generally based on the \strong{Wald method}.
The Wald approach to inference computes a test statistic by dividing the
parameter estimate by its standard error (Coefficient / SE),
then comparing this statistic against a t- or normal distribution.
This approach can be used to compute CIs and p-values.

\code{"wald"}:
\itemize{
\item Applies to \emph{non-Bayesian models}. For \emph{linear models}, CIs
computed using the Wald method (SE and a \emph{t-distribution with residual df});
p-values computed using the Wald method with a \emph{t-distribution with residual df}.
For other models, CIs computed using the Wald method (SE and a \emph{normal distribution});
p-values computed using the Wald method with a \emph{normal distribution}.
}

\code{"normal"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a normal distribution.
}

\code{"residual"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a \emph{t-distribution with residual df} when possible. If the
residual df for a model cannot be determined, a normal distribution is
used instead.
}

\strong{Methods for mixed models:}

Compared to fixed effects (or single-level) models, determining appropriate
df for Wald-based inference in mixed models is more difficult.
See \href{https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#what-are-the-p-values-listed-by-summaryglmerfit-etc.-are-they-reliable}{the R GLMM FAQ}
for a discussion.

Several approximate methods for computing df are available, but you should
also consider instead using profile likelihood (\code{"profile"}) or bootstrap ("\verb{boot"})
CIs and p-values instead.

\code{"satterthwaite"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the
Wald method (SE and a \emph{t-distribution with Satterthwaite df}); p-values
computed using the Wald method with a \emph{t-distribution with Satterthwaite df}.
}

\code{"kenward"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (\emph{Kenward-Roger SE} and a \emph{t-distribution with Kenward-Roger df});
p-values computed using the Wald method with \emph{Kenward-Roger SE and t-distribution with Kenward-Roger df}.
}

\code{"ml1"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (SE and a \emph{t-distribution with m-l-1 approximated df}); p-values
computed using the Wald method with a \emph{t-distribution with m-l-1 approximated df}.
See \code{\link[=ci_ml1]{ci_ml1()}}.
}

\code{"betwithin"}
\itemize{
\item Applies to \emph{linear mixed models} and \emph{generalized linear mixed models}.
CIs computed using the Wald method (SE and a \emph{t-distribution with between-within df});
p-values computed using the Wald method with a \emph{t-distribution with between-within df}.
See \code{\link[=ci_betwithin]{ci_betwithin()}}.
}

\strong{Likelihood-based methods:}

Likelihood-based inference is based on comparing the likelihood for the
maximum-likelihood estimate to the the likelihood for models with one or more
parameter values changed (e.g., set to zero or a range of alternative values).
Likelihood ratios for the maximum-likelihood and alternative models are compared
to a \eqn{\chi}-squared distribution to compute CIs and p-values.

\code{"profile"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glm}, \code{polr} or \code{glmmTMB}.
CIs computed by \emph{profiling the likelihood curve for a parameter}, using
linear interpolation to find where likelihood ratio equals a critical value;
p-values computed using the Wald method with a \emph{normal-distribution} (note:
this might change in a future update!)
}

\code{"uniroot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glmmTMB}. CIs
computed by \emph{profiling the likelihood curve for a parameter}, using root
finding to find where likelihood ratio equals a critical value; p-values
computed using the Wald method with a \emph{normal-distribution} (note: this
might change in a future update!)
}

\strong{Methods for bootstrapped or Bayesian models:}

Bootstrap-based inference is based on \strong{resampling} and refitting the model
to the resampled datasets. The distribution of parameter estimates across
resampled datasets is used to approximate the parameter's sampling
distribution. Depending on the type of model, several different methods for
bootstrapping and constructing CIs and p-values from the bootstrap
distribution are available.

For Bayesian models, inference is based on drawing samples from the model
posterior distribution.

\code{"quantile"} (or \code{"eti"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{equal tailed intervals} using the quantiles of the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:eti]{bayestestR::eti()}}.
}

\code{"hdi"}
\itemize{
\item Applies to \emph{all models (including Bayesian models)}. For non-Bayesian
models, only applies if \code{bootstrap = TRUE}. CIs computed as \emph{highest density intervals}
for the bootstrap or posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:hdi]{bayestestR::hdi()}}.
}

\code{"bci"} (or \code{"bcai"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{bias corrected and accelerated intervals} for the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:bci]{bayestestR::bci()}}.
}

\code{"si"}
\itemize{
\item Applies to \emph{Bayesian models} with proper priors. CIs computed as
\emph{support intervals} comparing the posterior samples against the prior samples;
p-values are based on the \emph{probability of direction}. See \code{\link[bayestestR:si]{bayestestR::si()}}.
}

\code{"boot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{merMod}. CIs computed
using \emph{parametric bootstrapping} (simulating data from the fitted model);
p-values computed using the Wald method with a \emph{normal-distribution)}
(note: this might change in a future update!).
}

For all iteration-based methods other than \code{"boot"}
(\code{"hdi"}, \code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, \code{"bcai"}),
p-values are based on the probability of direction (\code{\link[bayestestR:p_direction]{bayestestR::p_direction()}}),
which is converted into a p-value using \code{\link[bayestestR:pd_to_p]{bayestestR::pd_to_p()}}.
}

\examples{
\dontrun{
library(parameters)
if (require("rstanarm")) {
  model <- stan_glm(
    Sepal.Length ~ Petal.Length * Species,
    data = iris, iter = 500, refresh = 0
  )
  model_parameters(model)
}
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to
rename columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_PMCMRplus.R, R/methods_multcomp.R
\name{model_parameters.PMCMR}
\alias{model_parameters.PMCMR}
\alias{model_parameters.glht}
\title{Parameters from Hypothesis Testing}
\usage{
\method{model_parameters}{PMCMR}(model, ...)

\method{model_parameters}{glht}(model, ci = 0.95, exponentiate = FALSE, verbose = TRUE, ...)
}
\arguments{
\item{model}{Object of class \code{\link[multcomp:glht]{multcomp::glht()}} (\pkg{multcomp})
or of class \code{PMCMR}, \code{trendPMCMR} or \code{osrt} (\pkg{PMCMRplus}).}

\item{...}{Arguments passed to or from other methods. For instance, when
\code{bootstrap = TRUE}, arguments like \code{type} or \code{parallel} are
passed down to \code{bootstrap_model()}, and arguments like \code{ci_method}
are passed down to \code{\link[bayestestR:describe_posterior]{bayestestR::describe_posterior()}}.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{verbose}{Toggle warnings and messages.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from Hypothesis Testing.
}
\examples{
\donttest{
if (require("multcomp", quietly = TRUE)) {
  # multiple linear model, swiss data
  lmod <- lm(Fertility ~ ., data = swiss)
  mod <- glht(
    model = lmod,
    linfct = c(
      "Agriculture = 0",
      "Examination = 0",
      "Education = 0",
      "Catholic = 0",
      "Infant.Mortality = 0"
    )
  )
  model_parameters(mod)
}
if (require("PMCMRplus", quietly = TRUE)) {
  model <- kwAllPairsConoverTest(count ~ spray, data = InsectSprays)
  model_parameters(model)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_model.R
\name{bootstrap_model}
\alias{bootstrap_model}
\alias{bootstrap_model.default}
\alias{bootstrap_model.merMod}
\title{Model bootstrapping}
\usage{
bootstrap_model(model, iterations = 1000, ...)

\method{bootstrap_model}{default}(
  model,
  iterations = 1000,
  type = "ordinary",
  parallel = c("no", "multicore", "snow"),
  n_cpus = 1,
  verbose = FALSE,
  ...
)

\method{bootstrap_model}{merMod}(
  model,
  iterations = 1000,
  type = "parametric",
  parallel = c("no", "multicore", "snow"),
  n_cpus = 1,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{model}{Statistical model.}

\item{iterations}{The number of draws to simulate/bootstrap.}

\item{...}{Arguments passed to or from other methods.}

\item{type}{Character string specifying the type of bootstrap. For mixed models
of class \code{merMod} or \code{glmmTMB}, may be \code{"parametric"} (default) or
\code{"semiparametric"} (see \code{?lme4::bootMer} for details). For all
other models, see argument \code{sim} in \code{?boot::boot} (defaults to
\code{"ordinary"}).}

\item{parallel}{The type of parallel operation to be used (if any).}

\item{n_cpus}{Number of processes to be used in parallel operation.}

\item{verbose}{Toggle warnings and messages.}
}
\value{
A data frame of bootstrapped estimates.
}
\description{
Bootstrap a statistical model n times to return a data frame of estimates.
}
\details{
By default, \code{boot::boot()} is used to generate bootstraps from
the model data, which are then used to \code{update()} the model, i.e. refit
the model with the bootstrapped samples. For \code{merMod} objects (\strong{lme4})
or models from \strong{glmmTMB}, the \code{lme4::bootMer()} function is used to
obtain bootstrapped samples. \code{bootstrap_parameters()} summarizes the
bootstrapped model estimates.
}
\section{Using with \strong{emmeans}}{

The output can be passed directly to the various functions from the
\strong{emmeans} package, to obtain bootstrapped estimates, contrasts, simple
slopes, etc. and their confidence intervals. These can then be passed to
\code{model_parameter()} to obtain standard errors, p-values, etc. (see
example).
\cr\cr
Note that that p-values returned here are estimated under the assumption of
\emph{translation equivariance}: that shape of the sampling distribution is
unaffected by the null being true or not. If this assumption does not hold,
p-values can be biased, and it is suggested to use proper permutation tests
to obtain non-parametric p-values.
}

\examples{
\dontrun{
if (require("boot", quietly = TRUE)) {
  model <- lm(mpg ~ wt + factor(cyl), data = mtcars)
  b <- bootstrap_model(model)
  print(head(b))

  if (require("emmeans", quietly = TRUE)) {
    est <- emmeans(b, consec ~ cyl)
    print(model_parameters(est))
  }
}
}
}
\seealso{
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}, \code{\link[=simulate_model]{simulate_model()}}, \code{\link[=simulate_parameters]{simulate_parameters()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_parameters.R
\name{compare_parameters}
\alias{compare_parameters}
\alias{compare_models}
\title{Compare model parameters of multiple models}
\usage{
compare_parameters(
  ...,
  ci = 0.95,
  effects = "fixed",
  component = "conditional",
  standardize = NULL,
  exponentiate = FALSE,
  ci_method = "wald",
  p_adjust = NULL,
  style = NULL,
  column_names = NULL,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  df_method = ci_method
)

compare_models(
  ...,
  ci = 0.95,
  effects = "fixed",
  component = "conditional",
  standardize = NULL,
  exponentiate = FALSE,
  ci_method = "wald",
  p_adjust = NULL,
  style = NULL,
  column_names = NULL,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  df_method = ci_method
)
}
\arguments{
\item{...}{One or more regression model objects, or objects returned by
\code{model_parameters()}. Regression models may be of different model
types. Model objects may be passed comma separated, or as a list.
If model objects are passed with names or the list has named elements,
these names will be used as column names.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{effects}{Should parameters for fixed effects (\code{"fixed"}), random
effects (\code{"random"}), or both (\code{"all"}) be returned? Only applies
to mixed models. May be abbreviated. If the calculation of random effects
parameters takes too long, you may use \code{effects = "fixed"}.}

\item{component}{Model component for which parameters should be shown. See
documentation for related model class in \code{\link[=model_parameters]{model_parameters()}}.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{ci_method}{Method for computing degrees of freedom for p values
and confidence intervals (CI). See documentation for related model class
in \code{\link[=model_parameters]{model_parameters()}}.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}

\item{style}{String, indicating which style of output is requested. Following
templates are possible:
\itemize{
\item \code{"ci"}: Estimate and confidence intervals, no asterisks for p-values.
\item \code{"se"}: Estimate and standard errors, no asterisks for p-values.
\item \code{"ci_p"}: Estimate, confidence intervals and asterisks for p-values.
\item \code{"se_p"}: Estimate, standard errors and asterisks for p-values.
\item \code{"ci_p2"}: Estimate, confidence intervals and numeric p-values, in two columns.
\item \code{"se_p2"}: Estimate, standard errors and numeric p-values, in two columns.
}}

\item{column_names}{Character vector with strings that should be used as
column headers. Must be of same length as number of models in \code{...}.}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}

\item{verbose}{Toggle warnings and messages.}

\item{df_method}{Deprecated. Please use \code{ci_method}.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Compute and extract model parameters of multiple regression
models. See \code{\link[=model_parameters]{model_parameters()}} for further details.
}
\details{
This function is in an early stage and does not yet cope with more complex
models, and probably does not yet properly render all model components. It
should also be noted that when including models with interaction terms, not
only do the values of the parameters change, but so does their meaning (from
main effects, to simple slopes), thereby making such comparisons hard.
Therefore, you should not use this function to compare models with
interaction terms with models without interaction terms.
}
\examples{
data(iris)
lm1 <- lm(Sepal.Length ~ Species, data = iris)
lm2 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
compare_parameters(lm1, lm2)

data(mtcars)
m1 <- lm(mpg ~ wt, data = mtcars)
m2 <- glm(vs ~ wt + cyl, data = mtcars, family = "binomial")
compare_parameters(m1, m2)
\dontrun{
# exponentiate coefficients, but not for lm
compare_parameters(m1, m2, exponentiate = "nongaussian")

# change column names
compare_parameters("linear model" = m1, "logistic reg." = m2)
compare_parameters(m1, m2, column_names = c("linear model", "logistic reg."))

# or as list
compare_parameters(list(m1, m2))
compare_parameters(list("linear model" = m1, "logistic reg." = m2))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_model_parameters.R
\name{model_parameters}
\alias{model_parameters}
\alias{parameters}
\title{Model Parameters}
\usage{
model_parameters(model, ...)

parameters(model, ...)
}
\arguments{
\item{model}{Statistical Model.}

\item{...}{Arguments passed to or from other methods. Non-documented
arguments are \code{digits}, \code{p_digits}, \code{ci_digits} and
\code{footer_digits} to set the number of digits for the output.
\code{group} can also be passed to the \code{print()} method. See details
in \code{\link[=print.parameters_model]{print.parameters_model()}} and 'Examples' in
\code{\link[=model_parameters.default]{model_parameters.default()}}.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Compute and extract model parameters. See the documentation for your object's class:
\itemize{
\item{\link[=model_parameters.htest]{Correlations, t-tests, ...} (\code{htest}, \code{pairwise.htest})}
\item{\link[=model_parameters.aov]{ANOVAs} (\code{aov}, \code{anova}, \strong{afex}, ...)}
\item{\link[=model_parameters.default]{Regression models} (\code{lm}, \code{glm}, \strong{survey}, ...)}
\item{\link[=model_parameters.cgam]{Additive models} (\code{gam}, \code{gamm}, ...)}
\item{\link[=model_parameters.zcpglm]{Zero-inflated models} (\code{hurdle}, \code{zeroinfl}, \code{zerocount})}
\item{\link[=model_parameters.mlm]{Multinomial, ordinal and cumulative link models} (\code{bracl}, \code{multinom}, \code{mlm}, ...)}
\item{\link[=model_parameters.averaging]{Other special models} (\code{model.avg}, \code{betareg}, \code{glmx}, ...)}
\item{\link[=model_parameters.merMod]{Mixed models} (\pkg{lme4}, \pkg{nlme}, \pkg{glmmTMB}, \pkg{afex}, ...)}
\item{\link[=model_parameters.BFBayesFactor]{Bayesian tests} (\pkg{BayesFactor})}
\item{\link[=model_parameters.stanreg]{Bayesian models} (\pkg{rstanarm}, \pkg{brms}, \pkg{MCMCglmm}, \pkg{blavaan}, ...)}
\item{\link[=model_parameters.principal]{PCA and FA} (\pkg{psych})}
\item{\link[=model_parameters.lavaan]{CFA and SEM} (\pkg{lavaan})}
\item{\link[=model_parameters.kmeans]{Cluster models} (k-means, ...)}
\item{\link[=model_parameters.rma]{Meta-Analysis via linear (mixed) models} (\code{rma}, \code{metaplus}, \pkg{metaBMA}, ...)}
\item{\link[=model_parameters.glht]{Hypothesis testing} (\code{glht}, \pkg{PMCMRplus})}
\item{\link[=model_parameters.t1way]{Robust statistical tests} (\pkg{WRS2})}
\item{\link[=model_parameters.mira]{Multiply imputed repeated analyses} (\code{mira})}
}
}
\note{
The \code{\link[=print.parameters_model]{print()}} method has several
arguments to tweak the output. There is also a
\href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method}
implemented in the
\href{https://easystats.github.io/see/}{\strong{see}-package}, and a dedicated
method for use inside rmarkdown files,
\code{\link[=print_md.parameters_model]{print_md()}}.
}
\section{Standardization of model coefficients}{

Standardization is based on \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}. In case
of \code{standardize = "refit"}, the data used to fit the model will be
standardized and the model is completely refitted. In such cases, standard
errors and confidence intervals refer to the standardized coefficient. The
default, \code{standardize = "refit"}, never standardizes categorical predictors
(i.e. factors), which may be a different behaviour compared to other R
packages or other software packages (like SPSS). To mimic behaviour of SPSS
or packages such as \pkg{lm.beta}, use \code{standardize = "basic"}.
}

\section{Standardization Methods}{

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
\strong{Note} that \code{standardize_parameters(method = "refit")} may not return
the same results as fitting a model on data that has been standardized with
\code{standardize()}; \code{standardize_parameters()} used the data used by the model
fitting function, which might not be same data if there are missing values.
see the \code{remove_na} argument in \code{standardize()}.
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
within-group variable is found to have access between-group variance.
}
}

\section{Labeling the Degrees of Freedom}{

Throughout the \pkg{parameters} package, we decided to label the residual
degrees of freedom \emph{df_error}. The reason for this is that these degrees
of freedom not always refer to the residuals. For certain models, they refer
to the estimate error - in a linear model these are the same, but in - for
instance - any mixed effects model, this isn't strictly true. Hence, we
think that \code{df_error} is the most generic label for these degrees of
freedom.
}

\section{Confidence intervals and approximation of degrees of freedom}{

There are different ways of approximating the degrees of freedom depending
on different assumptions about the nature of the model and its sampling
distribution. The \code{ci_method} argument modulates the method for computing degrees
of freedom (df) that are used to calculate confidence intervals (CI) and the
related p-values. Following options are allowed, depending on the model
class:

\strong{Classical methods:}

Classical inference is generally based on the \strong{Wald method}.
The Wald approach to inference computes a test statistic by dividing the
parameter estimate by its standard error (Coefficient / SE),
then comparing this statistic against a t- or normal distribution.
This approach can be used to compute CIs and p-values.

\code{"wald"}:
\itemize{
\item Applies to \emph{non-Bayesian models}. For \emph{linear models}, CIs
computed using the Wald method (SE and a \emph{t-distribution with residual df});
p-values computed using the Wald method with a \emph{t-distribution with residual df}.
For other models, CIs computed using the Wald method (SE and a \emph{normal distribution});
p-values computed using the Wald method with a \emph{normal distribution}.
}

\code{"normal"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a normal distribution.
}

\code{"residual"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a \emph{t-distribution with residual df} when possible. If the
residual df for a model cannot be determined, a normal distribution is
used instead.
}

\strong{Methods for mixed models:}

Compared to fixed effects (or single-level) models, determining appropriate
df for Wald-based inference in mixed models is more difficult.
See \href{https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#what-are-the-p-values-listed-by-summaryglmerfit-etc.-are-they-reliable}{the R GLMM FAQ}
for a discussion.

Several approximate methods for computing df are available, but you should
also consider instead using profile likelihood (\code{"profile"}) or bootstrap ("\verb{boot"})
CIs and p-values instead.

\code{"satterthwaite"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the
Wald method (SE and a \emph{t-distribution with Satterthwaite df}); p-values
computed using the Wald method with a \emph{t-distribution with Satterthwaite df}.
}

\code{"kenward"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (\emph{Kenward-Roger SE} and a \emph{t-distribution with Kenward-Roger df});
p-values computed using the Wald method with \emph{Kenward-Roger SE and t-distribution with Kenward-Roger df}.
}

\code{"ml1"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (SE and a \emph{t-distribution with m-l-1 approximated df}); p-values
computed using the Wald method with a \emph{t-distribution with m-l-1 approximated df}.
See \code{\link[=ci_ml1]{ci_ml1()}}.
}

\code{"betwithin"}
\itemize{
\item Applies to \emph{linear mixed models} and \emph{generalized linear mixed models}.
CIs computed using the Wald method (SE and a \emph{t-distribution with between-within df});
p-values computed using the Wald method with a \emph{t-distribution with between-within df}.
See \code{\link[=ci_betwithin]{ci_betwithin()}}.
}

\strong{Likelihood-based methods:}

Likelihood-based inference is based on comparing the likelihood for the
maximum-likelihood estimate to the the likelihood for models with one or more
parameter values changed (e.g., set to zero or a range of alternative values).
Likelihood ratios for the maximum-likelihood and alternative models are compared
to a \eqn{\chi}-squared distribution to compute CIs and p-values.

\code{"profile"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glm}, \code{polr} or \code{glmmTMB}.
CIs computed by \emph{profiling the likelihood curve for a parameter}, using
linear interpolation to find where likelihood ratio equals a critical value;
p-values computed using the Wald method with a \emph{normal-distribution} (note:
this might change in a future update!)
}

\code{"uniroot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glmmTMB}. CIs
computed by \emph{profiling the likelihood curve for a parameter}, using root
finding to find where likelihood ratio equals a critical value; p-values
computed using the Wald method with a \emph{normal-distribution} (note: this
might change in a future update!)
}

\strong{Methods for bootstrapped or Bayesian models:}

Bootstrap-based inference is based on \strong{resampling} and refitting the model
to the resampled datasets. The distribution of parameter estimates across
resampled datasets is used to approximate the parameter's sampling
distribution. Depending on the type of model, several different methods for
bootstrapping and constructing CIs and p-values from the bootstrap
distribution are available.

For Bayesian models, inference is based on drawing samples from the model
posterior distribution.

\code{"quantile"} (or \code{"eti"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{equal tailed intervals} using the quantiles of the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:eti]{bayestestR::eti()}}.
}

\code{"hdi"}
\itemize{
\item Applies to \emph{all models (including Bayesian models)}. For non-Bayesian
models, only applies if \code{bootstrap = TRUE}. CIs computed as \emph{highest density intervals}
for the bootstrap or posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:hdi]{bayestestR::hdi()}}.
}

\code{"bci"} (or \code{"bcai"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{bias corrected and accelerated intervals} for the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:bci]{bayestestR::bci()}}.
}

\code{"si"}
\itemize{
\item Applies to \emph{Bayesian models} with proper priors. CIs computed as
\emph{support intervals} comparing the posterior samples against the prior samples;
p-values are based on the \emph{probability of direction}. See \code{\link[bayestestR:si]{bayestestR::si()}}.
}

\code{"boot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{merMod}. CIs computed
using \emph{parametric bootstrapping} (simulating data from the fitted model);
p-values computed using the Wald method with a \emph{normal-distribution)}
(note: this might change in a future update!).
}

For all iteration-based methods other than \code{"boot"}
(\code{"hdi"}, \code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, \code{"bcai"}),
p-values are based on the probability of direction (\code{\link[bayestestR:p_direction]{bayestestR::p_direction()}}),
which is converted into a p-value using \code{\link[bayestestR:pd_to_p]{bayestestR::pd_to_p()}}.
}

\section{Interpretation of Interaction Terms}{

Note that the \emph{interpretation} of interaction terms depends on many
characteristics of the model. The number of parameters, and overall
performance of the model, can differ \emph{or not} between \code{a * b}
\code{a : b}, and \code{a / b}, suggesting that sometimes interaction terms
give different parameterizations of the same model, but other times it gives
completely different models (depending on \code{a} or \code{b} being factors
of covariates, included as main effects or not, etc.). Their interpretation
depends of the full context of the model, which should not be inferred
from the parameters table alone - rather, we recommend to use packages
that calculate estimated marginal means or marginal effects, such as
\CRANpkg{modelbased}, \CRANpkg{emmeans} or \CRANpkg{ggeffects}. To raise
awareness for this issue, you may use \code{print(...,show_formula=TRUE)}
to add the model-specification to the output of the
\code{\link[=print.parameters_model]{print()}} method for \code{model_parameters()}.
}

\references{
\itemize{
\item Hoffman, L. (2015). Longitudinal analysis: Modeling within-person
fluctuation and change. Routledge.
\item Neter, J., Wasserman, W., & Kutner, M. H. (1989). Applied linear
regression models.
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to
rename columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_factorstructure.R
\name{check_kmo}
\alias{check_kmo}
\title{Kaiser, Meyer, Olkin (KMO) Measure of Sampling Adequacy (MSA) for Factor Analysis}
\usage{
check_kmo(x, ...)
}
\arguments{
\item{x}{A dataframe.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A list of indices related to KMO.
}
\description{
Kaiser (1970) introduced a Measure of Sampling Adequacy (MSA), later modified
by Kaiser and Rice (1974). The Kaiser-Meyer-Olkin (KMO) statistic, which can
vary from 0 to 1, indicates the degree to which each variable in a set is
predicted without error by the other variables.
}
\details{
A value of 0 indicates that the sum of partial correlations is large relative
to the sum correlations, indicating factor analysis is likely to be
inappropriate. A KMO value close to 1 indicates that the sum of partial
correlations is not large relative to the sum of correlations and so factor
analysis should yield distinct and reliable factors.

Kaiser (1975) suggested that KMO > .9 were marvelous, in the .80s,
meritorious, in the .70s, middling, in the .60s, mediocre, in the .50s,
miserable, and less than .5, unacceptable. Hair et al. (2006) suggest
accepting a value > 0.5. Values between 0.5 and 0.7 are mediocre, and values
between 0.7 and 0.8 are good.

This function is strongly inspired by the \code{KMO} function in the
\code{psych} package (Revelle, 2016). All credit goes to its author.
}
\examples{
library(parameters)
check_kmo(mtcars)
}
\references{
\itemize{
\item Revelle, W. (2016). How To: Use the psych package for Factor Analysis
and data reduction.
\item Kaiser, H. F. (1970). A second generation little jiffy.
Psychometrika, 35(4), 401-415.
\item Kaiser, H. F., & Rice, J. (1974). Little jiffy, mark IV. Educational
and psychological measurement, 34(1), 111-117.
\item Kaiser, H. F. (1974). An index of factorial simplicity.
Psychometrika, 39(1), 31-36.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_parameters.R,
%   R/select_parameters.stanreg.R
\name{select_parameters}
\alias{select_parameters}
\alias{select_parameters.lm}
\alias{select_parameters.merMod}
\alias{select_parameters.stanreg}
\title{Automated selection of model parameters}
\usage{
select_parameters(model, ...)

\method{select_parameters}{lm}(model, direction = "both", steps = 1000, k = 2, ...)

\method{select_parameters}{merMod}(model, direction = "backward", steps = 1000, ...)

\method{select_parameters}{stanreg}(model, method = NULL, cross_validation = FALSE, ...)
}
\arguments{
\item{model}{A statistical model (of class \code{lm}, \code{glm},
\code{merMod}, \code{stanreg} or \code{brmsfit}).}

\item{...}{Arguments passed to or from other methods.}

\item{direction}{
    the mode of stepwise search, can be one of \code{"both"},
    \code{"backward"}, or \code{"forward"}, with a default of \code{"both"}.
    If the \code{scope} argument is missing the default for
    \code{direction} is \code{"backward"}.  Values can be abbreviated.
  }

\item{steps}{
    the maximum number of steps to be considered.  The default is 1000
    (essentially as many as required).  It is typically used to stop the
    process early.
  }

\item{k}{
    the multiple of the number of degrees of freedom used for the penalty.
    Only \code{k = 2} gives the genuine AIC: \code{k = log(n)} is sometimes
    referred to as BIC or SBC.
  }

\item{method}{The method used in the variable selection. Can be \code{NULL}
(default), \code{"forward"} or \code{"L1"}. See \code{projpred::varsel}.}

\item{cross_validation}{Select with cross-validation.}
}
\value{
The model refitted with optimal number of parameters.
}
\description{
This function performs an automated selection of the 'best' parameters,
updating and returning the "best" model.
}
\details{
\subsection{Classical lm and glm}{
For frequentist GLMs, \code{select_parameters()} performs an AIC-based
stepwise selection.
}

\subsection{Mixed models}{
For mixed-effects models of class \code{merMod}, stepwise selection is
based on \code{\link[cAIC4:stepcAIC]{cAIC4::stepcAIC()}}. This step function
only searches the "best" model based on the random-effects structure,
i.e. \code{select_parameters()} adds or excludes random-effects until
the cAIC can't be improved further.
}

\subsection{Bayesian models}{
For Bayesian models, it uses the \pkg{projpred} package.
}
}
\examples{
model <- lm(mpg ~ ., data = mtcars)
select_parameters(model)

model <- lm(mpg ~ cyl * disp * hp * wt, data = mtcars)
select_parameters(model)
\donttest{
# lme4 -------------------------------------------
if (require("lme4")) {
  model <- lmer(
    Sepal.Width ~ Sepal.Length * Petal.Width * Petal.Length + (1 | Species),
    data = iris
  )
  select_parameters(model)
}
}

\dontrun{
# rstanarm -------------------------------------------
if (require("rstanarm") && require("projpred")) {
  model <- stan_glm(
    mpg ~ .,
    data = mtcars,
    iter = 500, refresh = 0, verbose = FALSE
  )
  select_parameters(model, cross_validation = TRUE)

  model <- stan_glm(
    mpg ~ cyl * disp * hp,
    data = mtcars,
    iter = 500, refresh = 0, verbose = FALSE
  )
  select_parameters(model, cross_validation = FALSE)
}
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_satterthwaite.R, R/dof_satterthwaite.R,
%   R/p_value_satterthwaite.R, R/standard_error_satterthwaite.R
\name{ci_satterthwaite}
\alias{ci_satterthwaite}
\alias{dof_satterthwaite}
\alias{p_value_satterthwaite}
\alias{se_satterthwaite}
\title{Satterthwaite approximation for SEs, CIs and p-values}
\usage{
ci_satterthwaite(model, ci = 0.95, robust = FALSE, ...)

dof_satterthwaite(model)

p_value_satterthwaite(model, dof = NULL, robust = FALSE, ...)

se_satterthwaite(model)
}
\arguments{
\item{model}{A statistical model.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{robust}{Logical, if \code{TRUE}, computes confidence intervals (or p-values)
based on robust standard errors. See \code{\link[=standard_error_robust]{standard_error_robust()}}.}

\item{...}{Arguments passed down to \code{\link[=standard_error_robust]{standard_error_robust()}}
when confidence intervals or p-values based on robust standard errors
should be computed.}

\item{dof}{Degrees of Freedom.}
}
\value{
A data frame.
}
\description{
An approximate F-test based on the Satterthwaite (1946) approach.
}
\details{
Inferential statistics (like p-values, confidence intervals and
standard errors) may be biased in mixed models when the number of clusters
is small (even if the sample size of level-1 units is high). In such cases
it is recommended to approximate a more accurate number of degrees of freedom
for such inferential statitics. Unlike simpler approximation heuristics
like the "m-l-1" rule (\code{dof_ml1}), the Satterthwaite approximation is
also applicable in more complex multilevel designs. However, the "m-l-1"
heuristic also applies to generalized mixed models, while approaches like
Kenward-Roger or Satterthwaite are limited to linear mixed models only.
}
\examples{
\donttest{
if (require("lme4", quietly = TRUE)) {
  model <- lmer(Petal.Length ~ Sepal.Length + (1 | Species), data = iris)
  p_value_satterthwaite(model)
}
}
}
\references{
Satterthwaite FE (1946) An approximate distribution of estimates of variance components. Biometrics Bulletin 2 (6):110–4.
}
\seealso{
\code{dof_satterthwaite()} and \code{se_satterthwaite()} are small helper-functions
to calculate approximated degrees of freedom and standard errors for model
parameters, based on the Satterthwaite (1946) approach.
\cr \cr
\code{\link[=dof_kenward]{dof_kenward()}} and \code{\link[=dof_ml1]{dof_ml1()}}
approximate degrees of freedom based on Kenward-Roger's method or the "m-l-1" rule.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_cplm.R, R/methods_pscl.R
\name{p_value.zcpglm}
\alias{p_value.zcpglm}
\alias{p_value.zeroinfl}
\title{p-values for Models with Zero-Inflation}
\usage{
\method{p_value}{zcpglm}(model, component = c("all", "conditional", "zi", "zero_inflated"), ...)

\method{p_value}{zeroinfl}(
  model,
  component = c("all", "conditional", "zi", "zero_inflated"),
  method = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{A statistical model.}

\item{component}{Model component for which parameters should be shown. See
the documentation for your object's class in \code{\link[=model_parameters]{model_parameters()}} for
further details.}

\item{...}{Arguments passed down to \code{standard_error_robust()} when confidence
intervals or p-values based on robust standard errors should be computed.
Only available for models where \code{method = "robust"} is supported.}

\item{method}{If \code{"robust"}, and if model is supported by the \pkg{sandwich}
or \pkg{clubSandwich} packages, computes p-values based on robust
covariance matrix estimation.}

\item{verbose}{Toggle warnings and messages.}
}
\value{
A data frame with at least two columns: the parameter names and the p-values.
Depending on the model, may also include columns for model components etc.
}
\description{
This function attempts to return, or compute, p-values of hurdle and
zero-inflated models.
}
\examples{
if (require("pscl", quietly = TRUE)) {
  data("bioChemists")
  model <- zeroinfl(art ~ fem + mar + kid5 | kid5 + phd, data = bioChemists)
  p_value(model)
  p_value(model, component = "zi")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_factors.R
\name{.n_factors_bentler}
\alias{.n_factors_bentler}
\title{Bentler and Yuan's Procedure}
\usage{
.n_factors_bentler(eigen_values = NULL, model = "factors", nobs = NULL)
}
\description{
Bentler and Yuan's Procedure
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_heterogeneity.R
\name{check_heterogeneity}
\alias{check_heterogeneity}
\title{Check model predictor for heterogeneity bias}
\usage{
check_heterogeneity(x, select = NULL, group = NULL)
}
\arguments{
\item{x}{A data frame or a mixed model object.}

\item{select}{Character vector (or formula) with names of variables to select
that should be checked. If \code{x} is a mixed model object, this argument
will be ignored.}

\item{group}{Character vector (or formula) with the name of the variable that
indicates the group- or cluster-ID. If \code{x} is a model object, this
argument will be ignored.}
}
\description{
\code{check_heterogeneity()} checks if model predictors or variables may
cause a heterogeneity bias, i.e. if variables have a within- and/or
between-effect.
}
\note{
This function will be removed in a future update. Please use
\code{performance::check_heterogeneity_bias()}.
}
\seealso{
For further details, see documentation for \code{?datawizard::demean}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_meta.R
\name{cluster_meta}
\alias{cluster_meta}
\title{Metaclustering}
\usage{
cluster_meta(list_of_clusters, rownames = NULL, ...)
}
\arguments{
\item{list_of_clusters}{A list of vectors with the clustering assignments from various methods.}

\item{rownames}{An optional vector of row.names for the matrix.}

\item{...}{Currently not used.}
}
\value{
A matrix containing all the pairwise (between each observation) probabilities of being clustered together by the methods.
}
\description{
One of the core "issue" of statistical clustering is that, in many cases, different methods will give different results. The \strong{metaclustering} approach proposed by \emph{easystats} (that finds echoes in \emph{consensus clustering}; see Monti et al., 2003) consists of treating the unique clustering solutions as a ensemble, from which we can derive a probability matrix. This matrix contains, for each pair of observations, the probability of being in the same cluster. For instance, if the 6th and the 9th row of a dataframe has been assigned to a similar cluster by 5 our of 10 clustering methods, then its probability of being grouped together is 0.5.
\cr\cr
Metaclustering is based on the hypothesis that, as each clustering algorithm embodies a different prism by which it sees the data, running an infinite amount of algorithms would result in the emergence of the "true" clusters. As the number of algorithms and parameters is finite, the probabilistic perspective is a useful proxy. This method is interesting where there is no obvious reasons to prefer one over another clustering method, as well as to investigate how robust some clusters are under different algorithms.
}
\examples{
\dontrun{
data <- iris[1:4]

rez1 <- cluster_analysis(data, n = 2, method = "kmeans")
rez2 <- cluster_analysis(data, n = 3, method = "kmeans")
rez3 <- cluster_analysis(data, n = 6, method = "kmeans")

list_of_clusters <- list(rez1, rez2, rez3)

m <- cluster_meta(list_of_clusters)

# Visualize matrix without reordering
heatmap(m, Rowv = NA, Colv = NA, scale = "none") # Without reordering
# Reordered heatmap
heatmap(m, scale = "none")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_parameters.R
\name{format_parameters}
\alias{format_parameters}
\alias{format_parameters.default}
\title{Parameter names formatting}
\usage{
format_parameters(model, ...)

\method{format_parameters}{default}(model, brackets = c("[", "]"), ...)
}
\arguments{
\item{model}{A statistical model.}

\item{...}{Currently not used.}

\item{brackets}{A character vector of length two, indicating the opening and closing brackets.}
}
\value{
A (names) character vector with formatted parameter names. The value names refer to the original names of the coefficients.
}
\description{
This functions formats the names of model parameters (coefficients)
to make them more human-readable.
}
\section{Interpretation of Interaction Terms}{

Note that the \emph{interpretation} of interaction terms depends on many
characteristics of the model. The number of parameters, and overall
performance of the model, can differ \emph{or not} between \code{a * b}
\code{a : b}, and \code{a / b}, suggesting that sometimes interaction terms
give different parameterizations of the same model, but other times it gives
completely different models (depending on \code{a} or \code{b} being factors
of covariates, included as main effects or not, etc.). Their interpretation
depends of the full context of the model, which should not be inferred
from the parameters table alone - rather, we recommend to use packages
that calculate estimated marginal means or marginal effects, such as
\CRANpkg{modelbased}, \CRANpkg{emmeans} or \CRANpkg{ggeffects}. To raise
awareness for this issue, you may use \code{print(...,show_formula=TRUE)}
to add the model-specification to the output of the
\code{\link[=print.parameters_model]{print()}} method for \code{model_parameters()}.
}

\examples{
model <- lm(Sepal.Length ~ Species * Sepal.Width, data = iris)
format_parameters(model)

model <- lm(Sepal.Length ~ Petal.Length + (Species / Sepal.Width), data = iris)
format_parameters(model)

model <- lm(Sepal.Length ~ Species + poly(Sepal.Width, 2), data = iris)
format_parameters(model)

model <- lm(Sepal.Length ~ Species + poly(Sepal.Width, 2, raw = TRUE), data = iris)
format_parameters(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_discrimination.R
\name{cluster_discrimination}
\alias{cluster_discrimination}
\title{Compute a linear discriminant analysis on classified cluster groups}
\usage{
cluster_discrimination(x, cluster_groups = NULL, ...)
}
\arguments{
\item{x}{A data frame}

\item{cluster_groups}{Group classification of the cluster analysis, which can
be retrieved from the \code{\link[=cluster_analysis]{cluster_analysis()}} function.}

\item{...}{Other arguments to be passed to or from.}
}
\description{
Computes linear discriminant analysis (LDA) on classified cluster groups, and determines the goodness of classification for each cluster group. See \code{MASS::lda()} for details.
}
\examples{
if (requireNamespace("MASS", quietly = TRUE)) {
  # Retrieve group classification from hierarchical cluster analysis
  clustering <- cluster_analysis(iris[, 1:4], n = 3)

  # Goodness of group classification
  cluster_discrimination(clustering)
}
}
\seealso{
\code{\link[=n_clusters]{n_clusters()}} to determine the number of clusters to extract, \code{\link[=cluster_analysis]{cluster_analysis()}} to compute a cluster analysis and \code{\link[=check_clusterstructure]{check_clusterstructure()}} to check suitability of data for clustering.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_order.R
\name{format_order}
\alias{format_order}
\title{Order (first, second, ...) formatting}
\usage{
format_order(order, textual = TRUE, ...)
}
\arguments{
\item{order}{value or vector of orders.}

\item{textual}{Return number as words. If \code{FALSE}, will run \code{\link[insight:format_value]{insight::format_value()}}.}

\item{...}{Arguments to be passed to \code{\link[insight:format_value]{format_value()}} if \code{textual} is \code{FALSE}.}
}
\value{
A formatted string.
}
\description{
Format order.
}
\examples{
format_order(2)
format_order(8)
format_order(25, textual = FALSE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_df_adjust.R
\name{format_df_adjust}
\alias{format_df_adjust}
\title{Format the name of the degrees-of-freedom adjustment methods}
\usage{
format_df_adjust(
  method,
  approx_string = "-approximated",
  dof_string = " degrees of freedom"
)
}
\arguments{
\item{method}{Name of the method.}

\item{approx_string, dof_string}{Suffix added to the name of the method in
the returned string.}
}
\value{
A formatted string.
}
\description{
Format the name of the degrees-of-freedom adjustment methods.
}
\examples{
library(parameters)

format_df_adjust("kenward")
format_df_adjust("kenward", approx_string = "", dof_string = " DoF")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_BayesFactor.R
\name{model_parameters.BFBayesFactor}
\alias{model_parameters.BFBayesFactor}
\title{Parameters from BayesFactor objects}
\usage{
\method{model_parameters}{BFBayesFactor}(
  model,
  centrality = "median",
  dispersion = FALSE,
  ci = 0.95,
  ci_method = "hdi",
  test = c("pd", "rope"),
  rope_range = "default",
  rope_ci = 0.95,
  priors = TRUE,
  cohens_d = NULL,
  cramers_v = NULL,
  include_proportions = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{Object of class \code{BFBayesFactor}.}

\item{centrality}{The point-estimates (centrality indices) to compute.  Character (vector) or list with one or more of these options: \code{"median"}, \code{"mean"}, \code{"MAP"} or \code{"all"}.}

\item{dispersion}{Logical, if \code{TRUE}, computes indices of dispersion related to the estimate(s) (\code{SD} and \code{MAD} for \code{mean} and \code{median}, respectively).}

\item{ci}{Value or vector of probability of the CI (between 0 and 1)
to be estimated. Default to \code{.95} (\verb{95\%}).}

\item{ci_method}{The type of index used for Credible Interval. Can be
\code{"HDI"} (default, see \code{\link[bayestestR:hdi]{hdi()}}), \code{"ETI"}
(see \code{\link[bayestestR:eti]{eti()}}), \code{"BCI"} (see
\code{\link[bayestestR:bci]{bci()}}) or \code{"SI"} (see \code{\link[bayestestR:si]{si()}}).}

\item{test}{The indices of effect existence to compute. Character (vector) or
list with one or more of these options: \code{"p_direction"} (or \code{"pd"}),
\code{"rope"}, \code{"p_map"}, \code{"equivalence_test"} (or \code{"equitest"}),
\code{"bayesfactor"} (or \code{"bf"}) or \code{"all"} to compute all tests.
For each "test", the corresponding \pkg{bayestestR} function is called
(e.g. \code{\link[bayestestR:rope]{rope()}} or \code{\link[bayestestR:p_direction]{p_direction()}}) and its results
included in the summary output.}

\item{rope_range}{ROPE's lower and higher bounds. Should be a list of two
values (e.g., \code{c(-0.1, 0.1)}) or \code{"default"}. If \code{"default"},
the bounds are set to \code{x +- 0.1*SD(response)}.}

\item{rope_ci}{The Credible Interval (CI) probability, corresponding to the
proportion of HDI, to use for the percentage in ROPE.}

\item{priors}{Add the prior used for each parameter.}

\item{cohens_d}{If \code{TRUE}, compute Cohens' \emph{d} as index of effect size. Only
applies to objects from \code{ttestBF()}. See \code{effectsize::cohens_d()} for
details.}

\item{cramers_v}{Compute Cramer's V or phi as index of effect size.
Can be \code{"raw"} or \code{"adjusted"} (effect size will be bias-corrected).
Only applies to objects from \code{chisq.test()}.}

\item{include_proportions}{Logical that decides whether to include posterior
cell proportions/counts for Bayesian contingency table analysis (from
\code{BayesFactor::contingencyTableBF()}). Defaults to \code{FALSE}, as this
information is often redundant.}

\item{verbose}{Toggle off warnings.}

\item{...}{Additional arguments to be passed to or from methods.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from \code{BFBayesFactor} objects from \code{{BayesFactor}} package.
}
\details{
The meaning of the extracted parameters:
\itemize{
\item For \code{\link[BayesFactor:ttestBF]{BayesFactor::ttestBF()}}: \code{Difference} is the raw
difference between the means. \item For
\code{\link[BayesFactor:correlationBF]{BayesFactor::correlationBF()}}: \code{rho} is the linear
correlation estimate (equivalent to Pearson's \emph{r}). \item For
\code{\link[BayesFactor:lmBF]{BayesFactor::lmBF()}} / \code{\link[BayesFactor:generalTestBF]{BayesFactor::generalTestBF()}}
/ \code{\link[BayesFactor:regressionBF]{BayesFactor::regressionBF()}} /
\code{\link[BayesFactor:anovaBF]{BayesFactor::anovaBF()}}: in addition to parameters of the fixed
and random effects, there are: \code{mu} is the (mean-centered) intercept;
\code{sig2} is the model's sigma; \code{g} / \verb{g_*} are the \emph{g}
parameters; See the \emph{Bayes Factors for ANOVAs} paper
(\doi{10.1016/j.jmp.2012.08.001}).
}
}
\examples{
\donttest{
if (require("BayesFactor")) {
  # Bayesian t-test
  model <- ttestBF(x = rnorm(100, 1, 1))
  model_parameters(model)
  model_parameters(model, cohens_d = TRUE, ci = .9)

  # Bayesian contingency table analysis
  data(raceDolls)
  bf <- contingencyTableBF(raceDolls, sampleType = "indepMulti", fixedMargin = "cols")
  model_parameters(bf,
    centrality = "mean",
    dispersion = TRUE,
    verbose = FALSE,
    cramers_v = TRUE
  )
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R, R/equivalence_test.R,
%   R/methods_bayestestR.R, R/n_parameters.R, R/print_md.R, R/reexports.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{display}
\alias{equivalence_test}
\alias{ci}
\alias{n_parameters}
\alias{print_md}
\alias{standardize_names}
\alias{supported_models}
\alias{print_html}
\alias{describe_distribution}
\alias{demean}
\alias{rescale_weights}
\alias{data_to_numeric}
\alias{convert_data_to_numeric}
\alias{skewness}
\alias{kurtosis}
\alias{smoothness}
\alias{center}
\alias{visualisation_recipe}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{bayestestR}{\code{\link[bayestestR]{ci}}, \code{\link[bayestestR]{equivalence_test}}}

  \item{datawizard}{\code{\link[datawizard]{center}}, \code{\link[datawizard]{convert_data_to_numeric}}, \code{\link[datawizard:convert_data_to_numeric]{data_to_numeric}}, \code{\link[datawizard]{demean}}, \code{\link[datawizard]{describe_distribution}}, \code{\link[datawizard:skewness]{kurtosis}}, \code{\link[datawizard]{rescale_weights}}, \code{\link[datawizard]{skewness}}, \code{\link[datawizard]{smoothness}}, \code{\link[datawizard]{visualisation_recipe}}}

  \item{insight}{\code{\link[insight]{display}}, \code{\link[insight]{n_parameters}}, \code{\link[insight:display]{print_html}}, \code{\link[insight:display]{print_md}}, \code{\link[insight]{standardize_names}}, \code{\link[insight:is_model_supported]{supported_models}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.factor_to_numeric}
\alias{.factor_to_numeric}
\title{Safe transformation from factor/character to numeric}
\usage{
.factor_to_numeric(x, lowest = NULL)
}
\description{
Safe transformation from factor/character to numeric
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_lavaan.R
\name{model_parameters.lavaan}
\alias{model_parameters.lavaan}
\title{Parameters from CFA/SEM models}
\usage{
\method{model_parameters}{lavaan}(
  model,
  ci = 0.95,
  standardize = FALSE,
  component = c("regression", "correlation", "loading", "defined"),
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{CFA or SEM created by the \code{lavaan::cfa} or \code{lavaan::sem}
functions.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{standardize}{Return standardized parameters (standardized coefficients).
Can be \code{TRUE} (or \code{"all"} or \code{"std.all"}) for standardized
estimates based on both the variances of observed and latent variables;
\code{"latent"} (or \code{"std.lv"}) for standardized estimates based
on the variances of the latent variables only; or \code{"no_exogenous"}
(or \code{"std.nox"}) for standardized estimates based on both the
variances of observed and latent variables, but not the variances of
exogenous covariates. See \code{lavaan::standardizedsolution} for details.}

\item{component}{What type of links to return. Can be \code{"all"} or some of \code{c("regression", "correlation", "loading", "variance", "mean")}.}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Format CFA/SEM objects from the lavaan package (Rosseel, 2012; Merkle and Rosseel 2018).
}
\note{
There is also a \href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
library(parameters)

# lavaan -------------------------------------
if (require("lavaan", quietly = TRUE)) {

  # Confirmatory Factor Analysis (CFA) ---------

  structure <- " visual  =~ x1 + x2 + x3
                 textual =~ x4 + x5 + x6
                 speed   =~ x7 + x8 + x9 "
  model <- lavaan::cfa(structure, data = HolzingerSwineford1939)
  model_parameters(model)
  model_parameters(model, standardize = TRUE)

  # filter parameters
  model_parameters(
    model,
    parameters = list(
      To = "^(?!visual)",
      From = "^(?!(x7|x8))"
    )
  )

  # Structural Equation Model (SEM) ------------

  structure <- "
    # latent variable definitions
      ind60 =~ x1 + x2 + x3
      dem60 =~ y1 + a*y2 + b*y3 + c*y4
      dem65 =~ y5 + a*y6 + b*y7 + c*y8
    # regressions
      dem60 ~ ind60
      dem65 ~ ind60 + dem60
    # residual correlations
      y1 ~~ y5
      y2 ~~ y4 + y6
      y3 ~~ y7
      y4 ~~ y8
      y6 ~~ y8
  "
  model <- lavaan::sem(structure, data = PoliticalDemocracy)
  model_parameters(model)
  model_parameters(model, standardize = TRUE)
}
}
\references{
\itemize{
\item Rosseel Y (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36.
\item Merkle EC , Rosseel Y (2018). blavaan: Bayesian Structural Equation
Models via Parameter Expansion. Journal of Statistical Software, 85(4),
1-30. http://www.jstatsoft.org/v85/i04/
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dof.R
\name{degrees_of_freedom}
\alias{degrees_of_freedom}
\alias{degrees_of_freedom.default}
\alias{dof}
\title{Degrees of Freedom (DoF)}
\usage{
degrees_of_freedom(model, ...)

\method{degrees_of_freedom}{default}(model, method = "analytical", ...)

dof(model, ...)
}
\arguments{
\item{model}{A statistical model.}

\item{...}{Currently not used.}

\item{method}{Can be \code{"analytical"} (default, DoFs are estimated based
on the model type), \code{"residual"} in which case they are directly taken
from the model if available (for Bayesian models, the goal (looking for
help to make it happen) would be to refit the model as a frequentist one
before extracting the DoFs), \code{"ml1"} (see \code{\link[=dof_ml1]{dof_ml1()}}), \code{"betwithin"}
(see \code{\link[=dof_betwithin]{dof_betwithin()}}), \code{"satterthwaite"} (see \code{\link[=dof_satterthwaite]{dof_satterthwaite()}}),
\code{"kenward"} (see \code{\link[=dof_kenward]{dof_kenward()}}) or \code{"any"}, which tries to extract DoF
by any of those methods, whichever succeeds. See 'Details'.}
}
\description{
Estimate or extract degrees of freedom of models parameters.
}
\details{
Methods for calculating degrees of freedom:
\itemize{
\item \code{"analytical"} for models of class \code{lmerMod}, Kenward-Roger approximated degrees of freedoms are calculated, for other models, \code{n-k} (number of observations minus number of parameters).
\item \code{"residual"} tries to extract residual degrees of freedom, and returns \code{Inf} if residual degrees of freedom could not be extracted.
\item \code{"any"} first tries to extract residual degrees of freedom, and if these are not available, extracts analytical degrees of freedom.
\item \code{"nokr"} same as \code{"analytical"}, but does not Kenward-Roger approximation for models of class \code{lmerMod}. Instead, always uses \code{n-k} to calculate df for any model.
\item \code{"normal"} returns \code{Inf}.
\item \code{"wald"} returns residual df for models with t-statistic, and \code{Inf} for all other models.
\item \code{"kenward"} calls \code{\link[=dof_kenward]{dof_kenward()}}.
\item \code{"satterthwaite"} calls \code{\link[=dof_satterthwaite]{dof_satterthwaite()}}.
\item \code{"ml1"} calls \code{\link[=dof_ml1]{dof_ml1()}}.
\item \code{"betwithin"} calls \code{\link[=dof_betwithin]{dof_betwithin()}}.
}
For models with z-statistic, the returned degrees of freedom for model parameters is \code{Inf} (unless \code{method = "ml1"} or \code{method = "betwithin"}), because there is only one distribution for the related test statistic.
}
\note{
In many cases, \code{degrees_of_freedom()} returns the same as \code{df.residuals()},
or \code{n-k} (number of observations minus number of parameters). However,
\code{degrees_of_freedom()} refers to the model's \emph{parameters} degrees of freedom
of the distribution for the related test statistic. Thus, for models with
z-statistic, results from \code{degrees_of_freedom()} and \code{df.residuals()} differ.
Furthermore, for other approximation methods like \code{"kenward"} or
\code{"satterthwaite"}, each model parameter can have a different degree of
freedom.
}
\examples{
model <- lm(Sepal.Length ~ Petal.Length * Species, data = iris)
dof(model)

model <- glm(vs ~ mpg * cyl, data = mtcars, family = "binomial")
dof(model)
\dontrun{
if (require("lme4", quietly = TRUE)) {
  model <- lmer(Sepal.Length ~ Petal.Length + (1 | Species), data = iris)
  dof(model)
}

if (require("rstanarm", quietly = TRUE)) {
  model <- stan_glm(
    Sepal.Length ~ Petal.Length * Species,
    data = iris,
    chains = 2,
    refresh = 0
  )
  dof(model)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pool_parameters.R
\name{pool_parameters}
\alias{pool_parameters}
\title{Pool Model Parameters}
\usage{
pool_parameters(
  x,
  exponentiate = FALSE,
  effects = "fixed",
  component = "conditional",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{A list of \code{parameters_model} objects, as returned by
\code{\link[=model_parameters]{model_parameters()}}, or a list of model-objects that is
supported by \code{model_parameters()}.}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{effects}{Should parameters for fixed effects (\code{"fixed"}), random
effects (\code{"random"}), or both (\code{"all"}) be returned? Only applies
to mixed models. May be abbreviated. If the calculation of random effects
parameters takes too long, you may use \code{effects = "fixed"}.}

\item{component}{Model component for which parameters should be shown. May be
one of \code{"conditional"}, \code{"precision"} (\pkg{betareg}),
\code{"scale"} (\pkg{ordinal}), \code{"extra"} (\pkg{glmx}),
\code{"marginal"} (\pkg{mfx}), \code{"conditional"} or \code{"full"} (for
\code{MuMIn::model.avg()}) or \code{"all"}.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Currently not used.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
This function "pools" (i.e. combines) model parameters in a similar fashion
as \code{mice::pool()}. However, this function pools parameters from
\code{parameters_model} objects, as returned by
\code{\link[=model_parameters]{model_parameters()}}.
}
\details{
Averaging of parameters follows Rubin's rules (\cite{Rubin, 1987, p. 76}).
The pooled degrees of freedom is based on the Barnard-Rubin adjustment for
small samples (\cite{Barnard and Rubin, 1999}).
}
\note{
Models with multiple components, (for instance, models with zero-inflation,
where predictors appear in the count and zero-inflated part) may fail in
case of identical names for coefficients in the different model components,
since the coefficient table is grouped by coefficient names for pooling. In
such cases, coefficients of count and zero-inflated model parts would be
combined. Therefore, the \code{component} argument defaults to
\code{"conditional"} to avoid this.
}
\examples{
# example for multiple imputed datasets
if (require("mice")) {
  data("nhanes2")
  imp <- mice(nhanes2, printFlag = FALSE)
  models <- lapply(1:5, function(i) {
    lm(bmi ~ age + hyp + chl, data = complete(imp, action = i))
  })
  pool_parameters(models)

  # should be identical to:
  m <- with(data = imp, exp = lm(bmi ~ age + hyp + chl))
  summary(pool(m))
}
}
\references{
Barnard, J. and Rubin, D.B. (1999). Small sample degrees of freedom with
multiple imputation. Biometrika, 86, 948-955. Rubin, D.B. (1987). Multiple
Imputation for Nonresponse in Surveys. New York: John Wiley and Sons.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_factors.R
\name{.n_factors_bartlett}
\alias{.n_factors_bartlett}
\title{Bartlett, Anderson and Lawley Procedures}
\usage{
.n_factors_bartlett(eigen_values = NULL, model = "factors", nobs = NULL)
}
\description{
Bartlett, Anderson and Lawley Procedures
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_ci.R, R/methods_glmmTMB.R, R/methods_lme4.R
\name{ci.default}
\alias{ci.default}
\alias{ci.glmmTMB}
\alias{ci.merMod}
\title{Confidence Intervals (CI)}
\usage{
\method{ci}{default}(x, ci = 0.95, dof = NULL, method = NULL, robust = FALSE, ...)

\method{ci}{glmmTMB}(
  x,
  ci = 0.95,
  dof = NULL,
  method = "wald",
  robust = FALSE,
  component = "all",
  verbose = TRUE,
  ...
)

\method{ci}{merMod}(
  x,
  ci = 0.95,
  dof = NULL,
  method = "wald",
  robust = FALSE,
  iterations = 500,
  ...
)
}
\arguments{
\item{x}{A statistical model.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{dof}{Number of degrees of freedom to be used when calculating
confidence intervals. If \code{NULL} (default), the degrees of freedom are
retrieved by calling \code{\link[=degrees_of_freedom]{degrees_of_freedom()}} with
approximation method defined in \code{method}. If not \code{NULL}, use this argument
to override the default degrees of freedom used to compute confidence
intervals.}

\item{method}{Method for computing degrees of freedom for
confidence intervals (CI) and the related p-values. Allowed are following
options (which vary depending on the model class): \code{"residual"},
\code{"normal"}, \code{"likelihood"}, \code{"satterthwaite"}, \code{"kenward"}, \code{"wald"},
\code{"profile"}, \code{"boot"}, \code{"uniroot"}, \code{"ml1"}, \code{"betwithin"}, \code{"hdi"},
\code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, or \code{"bcai"}. See section
\emph{Confidence intervals and approximation of degrees of freedom} in
\code{\link[=model_parameters]{model_parameters()}} for further details.}

\item{robust}{Logical, if \code{TRUE}, computes confidence intervals (or p-values)
based on robust standard errors. See \code{\link[=standard_error_robust]{standard_error_robust()}}.}

\item{...}{Arguments passed down to \code{\link[=standard_error_robust]{standard_error_robust()}}
when confidence intervals or p-values based on robust standard errors
should be computed.}

\item{component}{Model component for which parameters should be shown. See
the documentation for your object's class in \code{\link[=model_parameters]{model_parameters()}} for
further details.}

\item{verbose}{Toggle warnings and messages.}

\item{iterations}{The number of bootstrap replicates. Only applies to models
of class \code{merMod} when \code{method=boot}.}
}
\value{
A data frame containing the CI bounds.
}
\description{
Compute confidence intervals (CI) for frequentist models.
}
\note{
\code{ci_robust()} resp. \code{ci(robust=TRUE)} rely on the \pkg{sandwich}
or \pkg{clubSandwich} package (the latter if \code{vcov_estimation="CR"} for
cluster-robust standard errors) and will thus only work for those models
supported by those packages.
}
\section{Confidence intervals and approximation of degrees of freedom}{

There are different ways of approximating the degrees of freedom depending
on different assumptions about the nature of the model and its sampling
distribution. The \code{ci_method} argument modulates the method for computing degrees
of freedom (df) that are used to calculate confidence intervals (CI) and the
related p-values. Following options are allowed, depending on the model
class:

\strong{Classical methods:}

Classical inference is generally based on the \strong{Wald method}.
The Wald approach to inference computes a test statistic by dividing the
parameter estimate by its standard error (Coefficient / SE),
then comparing this statistic against a t- or normal distribution.
This approach can be used to compute CIs and p-values.

\code{"wald"}:
\itemize{
\item Applies to \emph{non-Bayesian models}. For \emph{linear models}, CIs
computed using the Wald method (SE and a \emph{t-distribution with residual df});
p-values computed using the Wald method with a \emph{t-distribution with residual df}.
For other models, CIs computed using the Wald method (SE and a \emph{normal distribution});
p-values computed using the Wald method with a \emph{normal distribution}.
}

\code{"normal"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a normal distribution.
}

\code{"residual"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a \emph{t-distribution with residual df} when possible. If the
residual df for a model cannot be determined, a normal distribution is
used instead.
}

\strong{Methods for mixed models:}

Compared to fixed effects (or single-level) models, determining appropriate
df for Wald-based inference in mixed models is more difficult.
See \href{https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#what-are-the-p-values-listed-by-summaryglmerfit-etc.-are-they-reliable}{the R GLMM FAQ}
for a discussion.

Several approximate methods for computing df are available, but you should
also consider instead using profile likelihood (\code{"profile"}) or bootstrap ("\verb{boot"})
CIs and p-values instead.

\code{"satterthwaite"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the
Wald method (SE and a \emph{t-distribution with Satterthwaite df}); p-values
computed using the Wald method with a \emph{t-distribution with Satterthwaite df}.
}

\code{"kenward"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (\emph{Kenward-Roger SE} and a \emph{t-distribution with Kenward-Roger df});
p-values computed using the Wald method with \emph{Kenward-Roger SE and t-distribution with Kenward-Roger df}.
}

\code{"ml1"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (SE and a \emph{t-distribution with m-l-1 approximated df}); p-values
computed using the Wald method with a \emph{t-distribution with m-l-1 approximated df}.
See \code{\link[=ci_ml1]{ci_ml1()}}.
}

\code{"betwithin"}
\itemize{
\item Applies to \emph{linear mixed models} and \emph{generalized linear mixed models}.
CIs computed using the Wald method (SE and a \emph{t-distribution with between-within df});
p-values computed using the Wald method with a \emph{t-distribution with between-within df}.
See \code{\link[=ci_betwithin]{ci_betwithin()}}.
}

\strong{Likelihood-based methods:}

Likelihood-based inference is based on comparing the likelihood for the
maximum-likelihood estimate to the the likelihood for models with one or more
parameter values changed (e.g., set to zero or a range of alternative values).
Likelihood ratios for the maximum-likelihood and alternative models are compared
to a \eqn{\chi}-squared distribution to compute CIs and p-values.

\code{"profile"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glm}, \code{polr} or \code{glmmTMB}.
CIs computed by \emph{profiling the likelihood curve for a parameter}, using
linear interpolation to find where likelihood ratio equals a critical value;
p-values computed using the Wald method with a \emph{normal-distribution} (note:
this might change in a future update!)
}

\code{"uniroot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glmmTMB}. CIs
computed by \emph{profiling the likelihood curve for a parameter}, using root
finding to find where likelihood ratio equals a critical value; p-values
computed using the Wald method with a \emph{normal-distribution} (note: this
might change in a future update!)
}

\strong{Methods for bootstrapped or Bayesian models:}

Bootstrap-based inference is based on \strong{resampling} and refitting the model
to the resampled datasets. The distribution of parameter estimates across
resampled datasets is used to approximate the parameter's sampling
distribution. Depending on the type of model, several different methods for
bootstrapping and constructing CIs and p-values from the bootstrap
distribution are available.

For Bayesian models, inference is based on drawing samples from the model
posterior distribution.

\code{"quantile"} (or \code{"eti"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{equal tailed intervals} using the quantiles of the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:eti]{bayestestR::eti()}}.
}

\code{"hdi"}
\itemize{
\item Applies to \emph{all models (including Bayesian models)}. For non-Bayesian
models, only applies if \code{bootstrap = TRUE}. CIs computed as \emph{highest density intervals}
for the bootstrap or posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:hdi]{bayestestR::hdi()}}.
}

\code{"bci"} (or \code{"bcai"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{bias corrected and accelerated intervals} for the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:bci]{bayestestR::bci()}}.
}

\code{"si"}
\itemize{
\item Applies to \emph{Bayesian models} with proper priors. CIs computed as
\emph{support intervals} comparing the posterior samples against the prior samples;
p-values are based on the \emph{probability of direction}. See \code{\link[bayestestR:si]{bayestestR::si()}}.
}

\code{"boot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{merMod}. CIs computed
using \emph{parametric bootstrapping} (simulating data from the fitted model);
p-values computed using the Wald method with a \emph{normal-distribution)}
(note: this might change in a future update!).
}

For all iteration-based methods other than \code{"boot"}
(\code{"hdi"}, \code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, \code{"bcai"}),
p-values are based on the probability of direction (\code{\link[bayestestR:p_direction]{bayestestR::p_direction()}}),
which is converted into a p-value using \code{\link[bayestestR:pd_to_p]{bayestestR::pd_to_p()}}.
}

\examples{
\donttest{
library(parameters)
if (require("glmmTMB")) {
  model <- glmmTMB(
    count ~ spp + mined + (1 | site),
    ziformula = ~mined,
    family = poisson(),
    data = Salamanders
  )

  ci(model)
  ci(model, component = "zi")
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_betwithin.R, R/dof_betwithin.R,
%   R/p_value_betwithin.R
\name{ci_betwithin}
\alias{ci_betwithin}
\alias{dof_betwithin}
\alias{p_value_betwithin}
\title{Between-within approximation for SEs, CIs and p-values}
\usage{
ci_betwithin(model, ci = 0.95, robust = FALSE, ...)

dof_betwithin(model)

p_value_betwithin(model, dof = NULL, robust = FALSE, ...)
}
\arguments{
\item{model}{A mixed model.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{robust}{Logical, if \code{TRUE}, computes confidence intervals (or p-values)
based on robust standard errors. See \code{\link[=standard_error_robust]{standard_error_robust()}}.}

\item{...}{Arguments passed down to \code{\link[=standard_error_robust]{standard_error_robust()}}
when confidence intervals or p-values based on robust standard errors
should be computed.}

\item{dof}{Degrees of Freedom.}
}
\value{
A data frame.
}
\description{
Approximation of degrees of freedom based on a "between-within" heuristic.
}
\details{
\subsection{Small Sample Cluster corrected Degrees of Freedom}{
Inferential statistics (like p-values, confidence intervals and
standard errors) may be biased in mixed models when the number of clusters
is small (even if the sample size of level-1 units is high). In such cases
it is recommended to approximate a more accurate number of degrees of freedom
for such inferential statistics (see \cite{Li and Redden 2015}). The
\emph{Between-within} denominator degrees of freedom approximation is
recommended in particular for (generalized) linear mixed models with repeated
measurements (longitudinal design). \code{dof_betwithin()} implements a heuristic
based on the between-within approach. \strong{Note} that this implementation
does not return exactly the same results as shown in \cite{Li and Redden 2015},
but similar.
}
\subsection{Degrees of Freedom for Longitudinal Designs (Repeated Measures)}{
In particular for repeated measure designs (longitudinal data analysis),
the \emph{between-within} heuristic is likely to be more accurate than simply
using the residual or infinite degrees of freedom, because \code{dof_betwithin()}
returns different degrees of freedom for within-cluster and between-cluster effects.
}
}
\examples{
\donttest{
if (require("lme4")) {
  data(sleepstudy)
  model <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
  dof_betwithin(model)
  p_value_betwithin(model)
}
}
}
\references{
\itemize{
\item Elff, M.; Heisig, J.P.; Schaeffer, M.; Shikano, S. (2019). Multilevel Analysis with Few Clusters: Improving Likelihood-based Methods to Provide Unbiased Estimates and Accurate Inference, British Journal of Political Science.
\item Li, P., Redden, D. T. (2015). Comparing denominator degrees of freedom approximations for the generalized linear mixed model in analyzing binary outcome in small sample cluster-randomized trials. BMC Medical Research Methodology, 15(1), 38. \doi{10.1186/s12874-015-0026-x}
}
}
\seealso{
\code{dof_betwithin()} is a small helper-function to calculate approximated
degrees of freedom of model parameters, based on the "between-within" heuristic.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_mice.R
\name{model_parameters.mira}
\alias{model_parameters.mira}
\title{Parameters from multiply imputed repeated analyses}
\usage{
\method{model_parameters}{mira}(
  model,
  ci = 0.95,
  exponentiate = FALSE,
  p_adjust = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{An object of class \code{mira}.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods.}
}
\description{
Format models of class \code{mira}, obtained from \code{mice::width.mids()}.
}
\details{
\code{model_parameters()} for objects of class \code{mira} works
similar to \code{summary(mice::pool())}, i.e. it generates the pooled summary
of multiple imputed repeated regression analyses.
}
\examples{
library(parameters)
if (require("mice", quietly = TRUE)) {
  data(nhanes2)
  imp <- mice(nhanes2)
  fit <- with(data = imp, exp = lm(bmi ~ age + hyp + chl))
  model_parameters(fit)
}
\dontrun{
# model_parameters() also works for models that have no "tidy"-method in mice
if (require("mice", quietly = TRUE) && require("gee", quietly = TRUE)) {
  data(warpbreaks)
  set.seed(1234)
  warpbreaks$tension[sample(1:nrow(warpbreaks), size = 10)] <- NA
  imp <- mice(warpbreaks)
  fit <- with(data = imp, expr = gee(breaks ~ tension, id = wool))

  # does not work:
  # summary(pool(fit))

  model_parameters(fit)
}
}



# and it works with pooled results
if (require("mice")) {
  data("nhanes2")
  imp <- mice(nhanes2)
  fit <- with(data = imp, exp = lm(bmi ~ age + hyp + chl))
  pooled <- pool(fit)

  model_parameters(pooled)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_performance.R
\name{cluster_performance}
\alias{cluster_performance}
\alias{cluster_performance.kmeans}
\alias{cluster_performance.hclust}
\alias{cluster_performance.dbscan}
\alias{cluster_performance.parameters_clusters}
\title{Performance of clustering models}
\usage{
cluster_performance(model, ...)

\method{cluster_performance}{kmeans}(model, ...)

\method{cluster_performance}{hclust}(model, data, clusters, ...)

\method{cluster_performance}{dbscan}(model, data, ...)

\method{cluster_performance}{parameters_clusters}(model, ...)
}
\arguments{
\item{model}{Cluster model.}

\item{...}{Arguments passed to or from other methods.}

\item{data}{A data.frame.}

\item{clusters}{A vector with clusters assignments (must be same length as rows in data).}
}
\description{
Compute performance indices for clustering solutions.
}
\examples{
# kmeans
model <- kmeans(iris[1:4], 3)
cluster_performance(model)
# hclust
data <- iris[1:4]
model <- hclust(dist(data))
clusters <- cutree(model, 3)

rez <- cluster_performance(model, data, clusters)
rez
# DBSCAN
if (require("dbscan", quietly = TRUE)) {
  model <- dbscan::dbscan(iris[1:4], eps = 1.45, minPts = 10)

  rez <- cluster_performance(model, iris[1:4])
  rez
}
# Retrieve performance from parameters
params <- model_parameters(kmeans(iris[1:4], 3))
cluster_performance(params)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_factors.R
\name{.n_factors_mreg}
\alias{.n_factors_mreg}
\title{Multiple Regression Procedure}
\usage{
.n_factors_mreg(eigen_values = NULL, model = "factors")
}
\description{
Multiple Regression Procedure
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_BayesFactor.R
\name{p_value.BFBayesFactor}
\alias{p_value.BFBayesFactor}
\title{p-values for Bayesian Models}
\usage{
\method{p_value}{BFBayesFactor}(model, ...)
}
\arguments{
\item{model}{A statistical model.}

\item{...}{Arguments passed down to \code{standard_error_robust()} when confidence
intervals or p-values based on robust standard errors should be computed.
Only available for models where \code{method = "robust"} is supported.}
}
\value{
The p-values.
}
\description{
This function attempts to return, or compute, p-values of Bayesian models.
}
\details{
For Bayesian models, the p-values corresponds to the \emph{probability of
direction} (\code{\link[bayestestR:p_direction]{bayestestR::p_direction()}}), which is converted to a p-value
using \code{bayestestR::convert_pd_to_p()}.
}
\examples{
data(iris)
model <- lm(Petal.Length ~ Sepal.Length + Species, data = iris)
p_value(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_factors.R
\name{.n_factors_scree}
\alias{.n_factors_scree}
\title{Non Graphical Cattell's Scree Test}
\usage{
.n_factors_scree(eigen_values = NULL, model = "factors")
}
\description{
Non Graphical Cattell's Scree Test
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factor_analysis.R, R/principal_components.R,
%   R/utils_pca_efa.R
\name{factor_analysis}
\alias{factor_analysis}
\alias{principal_components}
\alias{rotated_data}
\alias{predict.parameters_efa}
\alias{print.parameters_efa}
\alias{sort.parameters_efa}
\alias{closest_component}
\title{Principal Component Analysis (PCA) and Factor Analysis (FA)}
\usage{
factor_analysis(
  x,
  n = "auto",
  rotation = "none",
  sort = FALSE,
  threshold = NULL,
  standardize = TRUE,
  cor = NULL,
  ...
)

principal_components(
  x,
  n = "auto",
  rotation = "none",
  sort = FALSE,
  threshold = NULL,
  standardize = TRUE,
  ...
)

rotated_data(pca_results)

\method{predict}{parameters_efa}(object, newdata = NULL, names = NULL, keep_na = TRUE, ...)

\method{print}{parameters_efa}(x, digits = 2, sort = FALSE, threshold = NULL, labels = NULL, ...)

\method{sort}{parameters_efa}(x, ...)

closest_component(pca_results)
}
\arguments{
\item{x}{A data frame or a statistical model.}

\item{n}{Number of components to extract. If \code{n="all"}, then \code{n} is
set as the number of variables minus 1 (\code{ncol(x)-1}). If
\code{n="auto"} (default) or \code{n=NULL}, the number of components is
selected through \code{\link[=n_factors]{n_factors()}} resp. \code{\link[=n_components]{n_components()}}.
In \code{\link[=reduce_parameters]{reduce_parameters()}}, can also be \code{"max"}, in which case
it will select all the components that are maximally pseudo-loaded (i.e.,
correlated) by at least one variable.}

\item{rotation}{If not \code{"none"}, the PCA / FA will be computed using the
\pkg{psych} package. Possible options include \code{"varimax"},
\code{"quartimax"}, \code{"promax"}, \code{"oblimin"}, \code{"simplimax"},
or \code{"cluster"} (and more). See \code{\link[psych:fa]{psych::fa()}} for details.}

\item{sort}{Sort the loadings.}

\item{threshold}{A value between 0 and 1 indicates which (absolute) values
from the loadings should be removed. An integer higher than 1 indicates the
n strongest loadings to retain. Can also be \code{"max"}, in which case it
will only display the maximum loading per variable (the most simple
structure).}

\item{standardize}{A logical value indicating whether the variables should be
standardized (centered and scaled) to have unit variance before the
analysis (in general, such scaling is advisable).}

\item{cor}{An optional correlation matrix that can be used (note that the
data must still be passed as the first argument). If \code{NULL}, will
compute it by running \code{cor()} on the passed data.}

\item{...}{Arguments passed to or from other methods.}

\item{pca_results}{The output of the \code{principal_components()} function.}

\item{object}{An object of class \code{parameters_pca} or
\code{parameters_efa}}

\item{newdata}{An optional data frame in which to look for variables with
which to predict. If omitted, the fitted values are used.}

\item{names}{Optional character vector to name columns of the returned data
frame.}

\item{keep_na}{Logical, if \code{TRUE}, predictions also return observations
with missing values from the original data, hence the number of rows of
predicted data and original data is equal.}

\item{digits, labels}{Arguments for \code{print()}.}
}
\value{
A data frame of loadings.
}
\description{
The functions \code{principal_components()} and \code{factor_analysis()} can
be used to perform a principal component analysis (PCA) or a factor analysis
(FA). They return the loadings as a data frame, and various methods and
functions are available to access / display other information (see the
Details section).
}
\details{
\subsection{Methods and Utilities}{
\itemize{
\item \code{\link[=n_components]{n_components()}} and \code{\link[=n_factors]{n_factors()}} automatically
estimates the optimal number of dimensions to retain.

\item \code{\link[=check_factorstructure]{check_factorstructure()}} checks the suitability of the
data for factor analysis using the
\code{\link[=check_sphericity_bartlett]{sphericity()}} and the
\code{\link[=check_kmo]{sphericity()}} KMO measure.

\item{\code{\link[performance:check_itemscale]{performance::check_itemscale()}} computes various measures
of internal consistencies applied to the (sub)scales (i.e., components)
extracted from the PCA.}

\item{Running \code{summary} returns information related to each
component/factor, such as the explained variance and the Eivenvalues.}

\item{Running \code{\link[=get_scores]{get_scores()}} computes scores for each subscale.}

\item{Running \code{\link[=closest_component]{closest_component()}} will return a numeric vector
with the assigned component index for each column from the original data
frame.}

\item{Running \code{\link[=rotated_data]{rotated_data()}} will return the rotated data,
including missing values, so it matches the original data frame.}

\item{Running
\href{https://easystats.github.io/see/articles/parameters.html#principal-component-analysis}{\code{plot()}}
visually displays the loadings (that requires the
\href{https://easystats.github.io/see/}{\pkg{see} package} to work).}
}
}

\subsection{Complexity}{
Complexity represents the number of latent components needed to account
for the observed variables. Whereas a perfect simple structure solution
has a complexity of 1 in that each item would only load on one factor,
a solution with evenly distributed items has a complexity greater than 1
(\cite{Hofman, 1978; Pettersson and Turkheimer, 2010}) .
}

\subsection{Uniqueness}{
Uniqueness represents the variance that is 'unique' to the variable and
not shared with other variables. It is equal to \verb{1 – communality}
(variance that is shared with other variables). A uniqueness of \code{0.20}
suggests that \verb{20\%} or that variable's variance is not shared with other
variables in the overall factor model. The greater 'uniqueness' the lower
the relevance of the variable in the factor model.
}

\subsection{MSA}{
MSA represents the Kaiser-Meyer-Olkin Measure of Sampling Adequacy
(\cite{Kaiser and Rice, 1974}) for each item. It indicates whether there
is enough data for each factor give reliable results for the PCA. The
value should be > 0.6, and desirable values are > 0.8
(\cite{Tabachnick and Fidell, 2013}).
}

\subsection{PCA or FA?}{
There is a simplified rule of thumb that may help do decide whether to run
a factor analysis or a principal component analysis:
\itemize{
\item Run \emph{factor analysis} if you assume or wish to test a
theoretical model of \emph{latent factors} causing observed variables.

\item Run \emph{principal component analysis} If you want to simply
\emph{reduce} your correlated observed variables to a smaller set of
important independent composite variables.
}
(Source: \href{https://stats.stackexchange.com/q/1576/54740}{CrossValidated})
}

\subsection{Computing Item Scores}{
Use \code{\link[=get_scores]{get_scores()}} to compute scores for the "subscales"
represented by the extracted principal components. \code{get_scores()}
takes the results from \code{principal_components()} and extracts the
variables for each component found by the PCA. Then, for each of these
"subscales", raw means are calculated (which equals adding up the single
items and dividing by the number of items). This results in a sum score
for each component from the PCA, which is on the same scale as the
original, single items that were used to compute the PCA.
One can also use \code{predict()} to back-predict scores for each component,
to which one can provide \code{newdata} or a vector of \code{names} for the
components.
}

\subsection{Explained Variance and Eingenvalues}{
Use \code{summary()} to get the Eigenvalues and the explained variance
for each extracted component. The eigenvectors and eigenvalues represent
the "core" of a PCA: The eigenvectors (the principal components)
determine the directions of the new feature space, and the eigenvalues
determine their magnitude. In other words, the eigenvalues explain the
variance of the data along the new feature axes.
}
}
\examples{

library(parameters)

\donttest{
# Principal Component Analysis (PCA) -------------------
if (require("psych")) {
  principal_components(mtcars[, 1:7], n = "all", threshold = 0.2)
  principal_components(mtcars[, 1:7],
    n = 2, rotation = "oblimin",
    threshold = "max", sort = TRUE
  )
  principal_components(mtcars[, 1:7], n = 2, threshold = 2, sort = TRUE)

  pca <- principal_components(mtcars[, 1:5], n = 2, rotation = "varimax")
  pca # Print loadings
  summary(pca) # Print information about the factors
  predict(pca, names = c("Component1", "Component2")) # Back-predict scores

  # which variables from the original data belong to which extracted component?
  closest_component(pca)
  # rotated_data(pca)  # TODO: doesn't work
}

# Automated number of components
principal_components(mtcars[, 1:4], n = "auto")
}



# Factor Analysis (FA) ------------------------
if (require("psych")) {
  factor_analysis(mtcars[, 1:7], n = "all", threshold = 0.2)
  factor_analysis(mtcars[, 1:7], n = 2, rotation = "oblimin", threshold = "max", sort = TRUE)
  factor_analysis(mtcars[, 1:7], n = 2, threshold = 2, sort = TRUE)

  efa <- factor_analysis(mtcars[, 1:5], n = 2)
  summary(efa)
  predict(efa)
\donttest{
  # Automated number of components
  factor_analysis(mtcars[, 1:4], n = "auto")
}
}
}
\references{
\itemize{
\item Kaiser, H.F. and Rice. J. (1974). Little jiffy, mark iv. Educational
and Psychological Measurement, 34(1):111–117

\item Hofmann, R. (1978). Complexity and simplicity as objective indices
descriptive of factor solutions. Multivariate Behavioral Research, 13:2,
247-250, \doi{10.1207/s15327906mbr1302_9}

\item Pettersson, E., & Turkheimer, E. (2010). Item selection, evaluation,
and simple structure in personality data. Journal of research in
personality, 44(4), 407-420, \doi{10.1016/j.jrp.2010.03.002}

\item Tabachnick, B. G., and Fidell, L. S. (2013). Using multivariate
statistics (6th ed.). Boston: Pearson Education.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_efa_to_cfa.R
\name{convert_efa_to_cfa}
\alias{convert_efa_to_cfa}
\alias{convert_efa_to_cfa.fa}
\alias{efa_to_cfa}
\title{Conversion between EFA results and CFA structure}
\usage{
convert_efa_to_cfa(model, ...)

\method{convert_efa_to_cfa}{fa}(model, threshold = "max", names = NULL, ...)

efa_to_cfa(model, ...)
}
\arguments{
\item{model}{An EFA model (e.g., a \code{psych::fa} object).}

\item{...}{Arguments passed to or from other methods.}

\item{threshold}{A value between 0 and 1 indicates which (absolute) values
from the loadings should be removed. An integer higher than 1 indicates the
n strongest loadings to retain. Can also be \code{"max"}, in which case it
will only display the maximum loading per variable (the most simple
structure).}

\item{names}{Vector containing dimension names.}
}
\value{
Converted index.
}
\description{
Enables a conversion between Exploratory Factor Analysis (EFA) and
Confirmatory Factor Analysis (CFA) \code{lavaan}-ready structure.
}
\examples{
\donttest{
library(parameters)
if (require("psych") && require("lavaan")) {
  efa <- psych::fa(attitude, nfactors = 3)

  model1 <- efa_to_cfa(efa)
  model2 <- efa_to_cfa(efa, threshold = 0.3)

  anova(
    lavaan::cfa(model1, data = attitude),
    lavaan::cfa(model2, data = attitude)
  )
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.filter_component}
\alias{.filter_component}
\title{for models with zero-inflation component, return required component of model-summary}
\usage{
.filter_component(dat, component)
}
\description{
for models with zero-inflation component, return required component of model-summary
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_centers.R
\name{cluster_centers}
\alias{cluster_centers}
\title{Find the cluster centers in your data}
\usage{
cluster_centers(data, clusters, fun = mean, ...)
}
\arguments{
\item{data}{A data.frame.}

\item{clusters}{A vector with clusters assignments (must be same length as rows in data).}

\item{fun}{What function to use, \code{mean} by default.}

\item{...}{Other arguments to be passed to or from other functions.}
}
\value{
A dataframe containing the cluster centers. Attributes include performance statistics and distance between each observation and its respective cluster centre.
}
\description{
For each cluster, computes the mean (or other indices) of the variables. Can be used
to retrieve the centers of clusters. Also returns the within Sum of Squares.
}
\examples{
k <- kmeans(iris[1:4], 3)
cluster_centers(iris[1:4], clusters = k$cluster)
cluster_centers(iris[1:4], clusters = k$cluster, fun = median)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.parameters_model.R
\name{print.parameters_model}
\alias{print.parameters_model}
\alias{summary.parameters_model}
\title{Print model parameters}
\usage{
\method{print}{parameters_model}(
  x,
  pretty_names = TRUE,
  split_components = TRUE,
  select = NULL,
  caption = NULL,
  digits = 2,
  ci_digits = 2,
  p_digits = 3,
  footer_digits = 3,
  show_sigma = FALSE,
  show_formula = FALSE,
  zap_small = FALSE,
  groups = NULL,
  column_width = NULL,
  ci_brackets = c("[", "]"),
  ...
)

\method{summary}{parameters_model}(object, ...)
}
\arguments{
\item{x, object}{An object returned by \code{\link[=model_parameters]{model_parameters()}}.}

\item{pretty_names}{Return "pretty" (i.e. more human readable) parameter
names.}

\item{split_components}{Logical, if \code{TRUE} (default), For models with
multiple components (zero-inflation, smooth terms, ...), each component is
printed in a separate table. If \code{FALSE}, model parameters are printed
in a single table and a \code{Component} column is added to the output.}

\item{select}{Character vector (or numeric index) of column names that should
be printed. If \code{NULL} (default), all columns are printed. The shortcut
\code{select = "minimal"} prints coefficient, confidence intervals and p-values,
while \code{select = "short"} prints coefficient, standard errors and p-values.}

\item{caption}{Table caption as string. If \code{NULL}, no table caption is printed.}

\item{digits, ci_digits, p_digits}{Number of digits for rounding or
significant figures. May also be \code{"signif"} to return significant
figures or \code{"scientific"} to return scientific notation. Control the
number of digits by adding the value as suffix, e.g. \code{digits = "scientific4"}
to have scientific notation with 4 decimal places, or \code{digits = "signif5"}
for 5 significant figures (see also \code{\link[=signif]{signif()}}).}

\item{footer_digits}{Number of decimal places for values in the footer summary.}

\item{show_sigma}{Logical, if \code{TRUE}, adds information about the residual
standard deviation.}

\item{show_formula}{Logical, if \code{TRUE}, adds the model formula to the output.}

\item{zap_small}{Logical, if \code{TRUE}, small values are rounded after
\code{digits} decimal places. If \code{FALSE}, values with more decimal
places than \code{digits} are printed in scientific notation.}

\item{groups}{Named list, can be used to group parameters in the printed output.
List elements may either be character vectors that match the name of those
parameters that belong to one group, or list elements can be row numbers
of those parameter rows that should belong to one group. The names of the
list elements will be used as group names, which will be inserted as "header
row". A possible use case might be to emphasize focal predictors and control
variables, see 'Examples'. Parameters will be re-ordered according to the
order used in \code{groups}, while all non-matching parameters will be added
to the end.}

\item{column_width}{Width of table columns. Can be either \code{NULL}, a named
numeric vector, or \code{"fixed"}. If \code{NULL}, the width for each table column is
adjusted to the minimum required width. If a named numeric vector, value
names are matched against column names, and for each match, the specified
width is used. If \code{"fixed"}, and table is split into multiple components,
columns across all table components are adjusted to have the same width.}

\item{ci_brackets}{Logical, if \code{TRUE} (default), CI-values are
encompassed in square brackets (else in parentheses).}

\item{...}{Arguments passed to or from other methods.}
}
\value{
Invisibly returns the original input object.
}
\description{
A \code{print()}-method for objects from \code{\link[=model_parameters]{model_parameters()}}.
}
\details{
\code{summary()} is a convenient shortcut for
\code{print(object, select = "minimal", show_sigma = TRUE, show_formula = TRUE)}.
}
\section{Interpretation of Interaction Terms}{

Note that the \emph{interpretation} of interaction terms depends on many
characteristics of the model. The number of parameters, and overall
performance of the model, can differ \emph{or not} between \code{a * b}
\code{a : b}, and \code{a / b}, suggesting that sometimes interaction terms
give different parameterizations of the same model, but other times it gives
completely different models (depending on \code{a} or \code{b} being factors
of covariates, included as main effects or not, etc.). Their interpretation
depends of the full context of the model, which should not be inferred
from the parameters table alone - rather, we recommend to use packages
that calculate estimated marginal means or marginal effects, such as
\CRANpkg{modelbased}, \CRANpkg{emmeans} or \CRANpkg{ggeffects}. To raise
awareness for this issue, you may use \code{print(...,show_formula=TRUE)}
to add the model-specification to the output of the
\code{\link[=print.parameters_model]{print()}} method for \code{model_parameters()}.
}

\section{Labeling the Degrees of Freedom}{

Throughout the \pkg{parameters} package, we decided to label the residual
degrees of freedom \emph{df_error}. The reason for this is that these degrees
of freedom not always refer to the residuals. For certain models, they refer
to the estimate error - in a linear model these are the same, but in - for
instance - any mixed effects model, this isn't strictly true. Hence, we
think that \code{df_error} is the most generic label for these degrees of
freedom.
}

\examples{
\donttest{
library(parameters)
if (require("glmmTMB", quietly = TRUE)) {
  model <- glmmTMB(
    count ~ spp + mined + (1 | site),
    ziformula = ~mined,
    family = poisson(),
    data = Salamanders
  )
  mp <- model_parameters(model)

  print(mp, pretty_names = FALSE)

  print(mp, split_components = FALSE)

  print(mp, select = c("Parameter", "Coefficient", "SE"))

  print(mp, select = "minimal")
}


# group parameters ------

data(iris)
model <- lm(
  Sepal.Width ~ Sepal.Length + Species + Petal.Length,
  data = iris
)
# don't select "Intercept" parameter
mp <- model_parameters(model, parameters = "^(?!\\\\(Intercept)")
groups <- list(
  "Focal Predictors" = c("Speciesversicolor", "Speciesvirginica"),
  "Controls" = c("Sepal.Length", "Petal.Length")
)
print(mp, groups = groups)

# or use row indices
print(mp, groups = list(
  "Focal Predictors" = c(1, 4),
  "Controls" = c(2, 3)
))

# only show coefficients, CI and p,
# put non-matched parameters to the end

data(mtcars)
mtcars$cyl <- as.factor(mtcars$cyl)
mtcars$gear <- as.factor(mtcars$gear)
model <- lm(mpg ~ hp + gear * vs + cyl + drat, data = mtcars)

# don't select "Intercept" parameter
mp <- model_parameters(model, parameters = "^(?!\\\\(Intercept)")
print(mp, groups = list(
  "Engine" = c("cyl6", "cyl8", "vs", "hp"),
  "Interactions" = c("gear4:vs", "gear5:vs")
))
}
}
\seealso{
There is a dedicated method to use inside rmarkdown files,
\code{\link[=print_md.parameters_model]{print_md()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_cplm.R, R/methods_glmmTMB.R,
%   R/methods_lme4.R, R/methods_mixor.R, R/methods_ordinal.R
\name{model_parameters.cpglmm}
\alias{model_parameters.cpglmm}
\alias{model_parameters.glmmTMB}
\alias{model_parameters.merMod}
\alias{model_parameters.mixor}
\alias{model_parameters.clmm}
\title{Parameters from Mixed Models}
\usage{
\method{model_parameters}{cpglmm}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  effects = "all",
  group_level = FALSE,
  exponentiate = FALSE,
  ci_method = NULL,
  p_adjust = NULL,
  verbose = TRUE,
  df_method = ci_method,
  include_sigma = FALSE,
  ...
)

\method{model_parameters}{glmmTMB}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  effects = "all",
  component = "all",
  group_level = FALSE,
  standardize = NULL,
  exponentiate = FALSE,
  ci_method = "wald",
  robust = FALSE,
  p_adjust = NULL,
  wb_component = TRUE,
  summary = FALSE,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  df_method = ci_method,
  include_sigma = FALSE,
  ...
)

\method{model_parameters}{merMod}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  ci_method = NULL,
  iterations = 1000,
  standardize = NULL,
  effects = "all",
  group_level = FALSE,
  exponentiate = FALSE,
  robust = FALSE,
  p_adjust = NULL,
  wb_component = TRUE,
  summary = FALSE,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  verbose = TRUE,
  df_method = ci_method,
  include_sigma = FALSE,
  ...
)

\method{model_parameters}{mixor}(
  model,
  ci = 0.95,
  effects = "all",
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  exponentiate = FALSE,
  verbose = TRUE,
  include_sigma = FALSE,
  ...
)

\method{model_parameters}{clmm}(
  model,
  ci = 0.95,
  bootstrap = FALSE,
  iterations = 1000,
  standardize = NULL,
  effects = "all",
  group_level = FALSE,
  exponentiate = FALSE,
  ci_method = NULL,
  p_adjust = NULL,
  verbose = TRUE,
  df_method = ci_method,
  include_sigma = FALSE,
  ...
)
}
\arguments{
\item{model}{A mixed model.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{bootstrap}{Should estimates be based on bootstrapped model? If
\code{TRUE}, then arguments of \link[=model_parameters.stanreg]{Bayesian regressions} apply (see also
\code{\link[=bootstrap_parameters]{bootstrap_parameters()}}).}

\item{iterations}{The number of draws to simulate/bootstrap.}

\item{standardize}{The method used for standardizing the parameters. Can be
\code{NULL} (default; no standardization), \code{"refit"} (for re-fitting the model
on standardized data) or one of \code{"basic"}, \code{"posthoc"}, \code{"smart"},
\code{"pseudo"}. See 'Details' in \code{\link[effectsize:standardize_parameters]{effectsize::standardize_parameters()}}.
\strong{Important:}
\itemize{
\item The \code{"refit"} method does \emph{not} standardized categorical predictors (i.e.
factors), which may be a different behaviour compared to other R packages
(such as \pkg{lm.beta}) or other software packages (like SPSS). to mimic
such behaviours, either use \code{standardize="basic"} or standardize the data
with \code{datawizard::standardize(force=TRUE)} \emph{before} fitting the model.
\item For mixed models, when using methods other than \code{"refit"}, only the fixed
effects will be returned.
\item Robust estimation (i.e. \code{robust=TRUE}) of standardized parameters only
works when \code{standardize="refit"}.
}}

\item{effects}{Should parameters for fixed effects (\code{"fixed"}), random
effects (\code{"random"}), or both (\code{"all"}) be returned? Only applies
to mixed models. May be abbreviated. If the calculation of random effects
parameters takes too long, you may use \code{effects = "fixed"}.}

\item{group_level}{Logical, for multilevel models (i.e. models with random
effects) and when \code{effects = "all"} or \code{effects = "random"},
include the parameters for each group level from random effects. If
\code{group_level = FALSE} (the default), only information on SD and COR
are shown.}

\item{exponentiate}{Logical, indicating whether or not to exponentiate the
the coefficients (and related confidence intervals). This is typical for
logistic regression, or more generally speaking, for models with log
or logit links. \strong{Note:} Delta-method standard errors are also
computed (by multiplying the standard errors by the transformed
coefficients). This is to mimic behaviour of other software packages, such
as Stata, but these standard errors poorly estimate uncertainty for the
transformed coefficient. The transformed confidence interval more clearly
captures this uncertainty. For \code{compare_parameters()},
\code{exponentiate = "nongaussian"} will only exponentiate coefficients
from non-Gaussian families.}

\item{ci_method}{Method for computing degrees of freedom for
confidence intervals (CI) and the related p-values. Allowed are following
options (which vary depending on the model class): \code{"residual"},
\code{"normal"}, \code{"likelihood"}, \code{"satterthwaite"}, \code{"kenward"}, \code{"wald"},
\code{"profile"}, \code{"boot"}, \code{"uniroot"}, \code{"ml1"}, \code{"betwithin"}, \code{"hdi"},
\code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, or \code{"bcai"}. See section
\emph{Confidence intervals and approximation of degrees of freedom} in
\code{\link[=model_parameters]{model_parameters()}} for further details. When \code{ci_method=NULL}, in most
cases \code{"wald"} is used then.}

\item{p_adjust}{Character vector, if not \code{NULL}, indicates the method to
adjust p-values. See \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Further
possible adjustment methods are \code{"tukey"}, \code{"scheffe"},
\code{"sidak"} and \code{"none"} to explicitly disable adjustment for
\code{emmGrid} objects (from \pkg{emmeans}).}

\item{verbose}{Toggle warnings and messages.}

\item{df_method}{Deprecated. Please use \code{ci_method}.}

\item{include_sigma}{Logical, if \code{TRUE}, includes the residual standard
deviation. For mixed models, this is defined as the sum of the distribution-specific
variance and the variance for the additive overdispersion term (see
\code{\link[insight:get_variance]{insight::get_variance()}} for details). Defaults to \code{FALSE} for mixed models
due to the longer computation time.}

\item{...}{Arguments passed to or from other methods.}

\item{component}{Should all parameters, parameters for the conditional model,
or for the zero-inflated part of the model be returned? Applies to models
with zero-inflated component. \code{component} may be one of \code{"conditional"},
\code{"zi"}, \code{"zero-inflated"}, \code{"dispersion"} or \code{"all"}
(default). May be abbreviated.}

\item{robust}{Logical, if \code{TRUE}, robust standard errors are calculated
(if possible), and confidence intervals and p-values are based on these
robust standard errors. Additional arguments like \code{vcov_estimation} or
\code{vcov_type} are passed down to other methods, see
\code{\link[=standard_error_robust]{standard_error_robust()}} for details
and \href{https://easystats.github.io/parameters/articles/model_parameters_robust.html}{this vignette}
for working examples.}

\item{wb_component}{Logical, if \code{TRUE} and models contains within- and
between-effects (see \code{datawizard::demean()}), the \code{Component} column
will indicate which variables belong to the within-effects,
between-effects, and cross-level interactions. By default, the
\code{Component} column indicates, which parameters belong to the
conditional or zero-inflated component of the model.}

\item{summary}{Logical, if \code{TRUE}, prints summary information about the
model (model formula, number of observations, residual standard deviation
and more).}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from (linear) mixed models.
}
\note{
If the calculation of random effects parameters takes too long, you may
use \code{effects = "fixed"}. There is also a \href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\section{Confidence intervals for random effect variances}{

For models of class \code{merMod} and \code{glmmTMB}, confidence intervals for random
effect variances can be calculated. For models of class \code{lme4}, when
\code{ci_method} is either \code{"profile"} or \code{"boot"}, and \code{effects} is either
\code{"random"} or \code{"all"}, profiled resp. bootstrapped confidence intervals are
computed for the random effects. For all other options of \code{ci_method},
confidence intervals for random effects will be missing. For models of class
\code{glmmTMB}, confidence intervals for random effect variances always use a
Wald t-distribution approximation.
}

\section{Confidence intervals and approximation of degrees of freedom}{

There are different ways of approximating the degrees of freedom depending
on different assumptions about the nature of the model and its sampling
distribution. The \code{ci_method} argument modulates the method for computing degrees
of freedom (df) that are used to calculate confidence intervals (CI) and the
related p-values. Following options are allowed, depending on the model
class:

\strong{Classical methods:}

Classical inference is generally based on the \strong{Wald method}.
The Wald approach to inference computes a test statistic by dividing the
parameter estimate by its standard error (Coefficient / SE),
then comparing this statistic against a t- or normal distribution.
This approach can be used to compute CIs and p-values.

\code{"wald"}:
\itemize{
\item Applies to \emph{non-Bayesian models}. For \emph{linear models}, CIs
computed using the Wald method (SE and a \emph{t-distribution with residual df});
p-values computed using the Wald method with a \emph{t-distribution with residual df}.
For other models, CIs computed using the Wald method (SE and a \emph{normal distribution});
p-values computed using the Wald method with a \emph{normal distribution}.
}

\code{"normal"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a normal distribution.
}

\code{"residual"}
\itemize{
\item Applies to \emph{non-Bayesian models}. Compute Wald CIs and p-values,
but always use a \emph{t-distribution with residual df} when possible. If the
residual df for a model cannot be determined, a normal distribution is
used instead.
}

\strong{Methods for mixed models:}

Compared to fixed effects (or single-level) models, determining appropriate
df for Wald-based inference in mixed models is more difficult.
See \href{https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#what-are-the-p-values-listed-by-summaryglmerfit-etc.-are-they-reliable}{the R GLMM FAQ}
for a discussion.

Several approximate methods for computing df are available, but you should
also consider instead using profile likelihood (\code{"profile"}) or bootstrap ("\verb{boot"})
CIs and p-values instead.

\code{"satterthwaite"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the
Wald method (SE and a \emph{t-distribution with Satterthwaite df}); p-values
computed using the Wald method with a \emph{t-distribution with Satterthwaite df}.
}

\code{"kenward"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (\emph{Kenward-Roger SE} and a \emph{t-distribution with Kenward-Roger df});
p-values computed using the Wald method with \emph{Kenward-Roger SE and t-distribution with Kenward-Roger df}.
}

\code{"ml1"}
\itemize{
\item Applies to \emph{linear mixed models}. CIs computed using the Wald
method (SE and a \emph{t-distribution with m-l-1 approximated df}); p-values
computed using the Wald method with a \emph{t-distribution with m-l-1 approximated df}.
See \code{\link[=ci_ml1]{ci_ml1()}}.
}

\code{"betwithin"}
\itemize{
\item Applies to \emph{linear mixed models} and \emph{generalized linear mixed models}.
CIs computed using the Wald method (SE and a \emph{t-distribution with between-within df});
p-values computed using the Wald method with a \emph{t-distribution with between-within df}.
See \code{\link[=ci_betwithin]{ci_betwithin()}}.
}

\strong{Likelihood-based methods:}

Likelihood-based inference is based on comparing the likelihood for the
maximum-likelihood estimate to the the likelihood for models with one or more
parameter values changed (e.g., set to zero or a range of alternative values).
Likelihood ratios for the maximum-likelihood and alternative models are compared
to a \eqn{\chi}-squared distribution to compute CIs and p-values.

\code{"profile"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glm}, \code{polr} or \code{glmmTMB}.
CIs computed by \emph{profiling the likelihood curve for a parameter}, using
linear interpolation to find where likelihood ratio equals a critical value;
p-values computed using the Wald method with a \emph{normal-distribution} (note:
this might change in a future update!)
}

\code{"uniroot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{glmmTMB}. CIs
computed by \emph{profiling the likelihood curve for a parameter}, using root
finding to find where likelihood ratio equals a critical value; p-values
computed using the Wald method with a \emph{normal-distribution} (note: this
might change in a future update!)
}

\strong{Methods for bootstrapped or Bayesian models:}

Bootstrap-based inference is based on \strong{resampling} and refitting the model
to the resampled datasets. The distribution of parameter estimates across
resampled datasets is used to approximate the parameter's sampling
distribution. Depending on the type of model, several different methods for
bootstrapping and constructing CIs and p-values from the bootstrap
distribution are available.

For Bayesian models, inference is based on drawing samples from the model
posterior distribution.

\code{"quantile"} (or \code{"eti"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{equal tailed intervals} using the quantiles of the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:eti]{bayestestR::eti()}}.
}

\code{"hdi"}
\itemize{
\item Applies to \emph{all models (including Bayesian models)}. For non-Bayesian
models, only applies if \code{bootstrap = TRUE}. CIs computed as \emph{highest density intervals}
for the bootstrap or posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:hdi]{bayestestR::hdi()}}.
}

\code{"bci"} (or \code{"bcai"})
\itemize{
\item Applies to \emph{all models (including Bayesian models)}.
For non-Bayesian models, only applies if \code{bootstrap = TRUE}. CIs computed
as \emph{bias corrected and accelerated intervals} for the bootstrap or
posterior samples; p-values are based on the \emph{probability of direction}.
See \code{\link[bayestestR:bci]{bayestestR::bci()}}.
}

\code{"si"}
\itemize{
\item Applies to \emph{Bayesian models} with proper priors. CIs computed as
\emph{support intervals} comparing the posterior samples against the prior samples;
p-values are based on the \emph{probability of direction}. See \code{\link[bayestestR:si]{bayestestR::si()}}.
}

\code{"boot"}
\itemize{
\item Applies to \emph{non-Bayesian models} of class \code{merMod}. CIs computed
using \emph{parametric bootstrapping} (simulating data from the fitted model);
p-values computed using the Wald method with a \emph{normal-distribution)}
(note: this might change in a future update!).
}

For all iteration-based methods other than \code{"boot"}
(\code{"hdi"}, \code{"quantile"}, \code{"ci"}, \code{"eti"}, \code{"si"}, \code{"bci"}, \code{"bcai"}),
p-values are based on the probability of direction (\code{\link[bayestestR:p_direction]{bayestestR::p_direction()}}),
which is converted into a p-value using \code{\link[bayestestR:pd_to_p]{bayestestR::pd_to_p()}}.
}

\examples{
library(parameters)
if (require("lme4")) {
  data(mtcars)
  model <- lmer(mpg ~ wt + (1 | gear), data = mtcars)
  model_parameters(model)
}
\donttest{
if (require("glmmTMB")) {
  data(Salamanders)
  model <- glmmTMB(
    count ~ spp + mined + (1 | site),
    ziformula = ~mined,
    family = poisson(),
    data = Salamanders
  )
  model_parameters(model, effects = "all")
}

if (require("lme4")) {
  model <- lmer(mpg ~ wt + (1 | gear), data = mtcars)
  model_parameters(model, bootstrap = TRUE, iterations = 50)
}
}
}
\seealso{
\code{\link[insight:standardize_names]{insight::standardize_names()}} to
rename columns into a consistent, standardized naming scheme.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_aov.R
\name{model_parameters.aov}
\alias{model_parameters.aov}
\title{Parameters from ANOVAs}
\usage{
\method{model_parameters}{aov}(
  model,
  omega_squared = NULL,
  eta_squared = NULL,
  epsilon_squared = NULL,
  df_error = NULL,
  type = NULL,
  ci = NULL,
  alternative = NULL,
  test = NULL,
  power = FALSE,
  keep = NULL,
  drop = NULL,
  parameters = keep,
  table_wide = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{Object of class \code{\link[=aov]{aov()}}, \code{\link[=anova]{anova()}},
\code{aovlist}, \code{Gam}, \code{\link[=manova]{manova()}}, \code{Anova.mlm},
\code{afex_aov} or \code{maov}.}

\item{omega_squared}{Compute omega squared as index of effect size. Can be
\code{"partial"} (the default, adjusted for effect size) or \code{"raw"}.}

\item{eta_squared}{Compute eta squared as index of effect size. Can be
\code{"partial"} (the default, adjusted for effect size), \code{"raw"}  or
\code{"adjusted"} (the latter option only for ANOVA-tables from mixed
models).}

\item{epsilon_squared}{Compute epsilon squared as index of effect size. Can
be \code{"partial"} (the default, adjusted for effect size) or
\code{"raw"}.}

\item{df_error}{Denominator degrees of freedom (or degrees of freedom of the
error estimate, i.e., the residuals). This is used to compute effect sizes
for ANOVA-tables from mixed models. See 'Examples'. (Ignored for
\code{afex_aov}.)}

\item{type}{Numeric, type of sums of squares. May be 1, 2 or 3. If 2 or 3,
ANOVA-tables using \code{car::Anova()} will be returned. (Ignored for
\code{afex_aov}.)}

\item{ci}{Confidence Interval (CI) level for effect sizes
\code{omega_squared}, \code{eta_squared} etc. The default, \code{NULL},
will compute no confidence intervals. \code{ci} should be a scalar between
0 and 1.}

\item{alternative}{A character string specifying the alternative hypothesis;
Controls the type of CI returned: \code{"two.sided"} (default, two-sided CI),
\code{"greater"} or \code{"less"} (one-sided CI). Partial matching is allowed
(e.g., \code{"g"}, \code{"l"}, \code{"two"}...). See section \emph{One-Sided CIs} in
the \href{https://easystats.github.io/effectsize/}{effectsize_CIs vignette}.}

\item{test}{String, indicating the type of test for \code{Anova.mlm} to be
returned. If \code{"multivariate"} (or \code{NULL}), returns the summary of
the multivariate test (that is also given by the \code{print}-method). If
\code{test = "univariate"}, returns the summary of the univariate test.}

\item{power}{Logical, if \code{TRUE}, adds a column with power for each
parameter.}

\item{keep}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{drop}{Character containing a regular expression pattern that
describes the parameters that should be included (for \code{keep}) or excluded
(for \code{drop}) in the returned data frame. \code{keep} may also be a
named list of regular expressions. All non-matching parameters will be
removed from the output. If \code{keep} is a character vector, every parameter
name in the \emph{"Parameter"} column that matches the regular expression in
\code{keep} will be selected from the returned data frame (and vice versa,
all parameter names matching \code{drop} will be excluded). Furthermore, if
\code{keep} has more than one element, these will be merged with an \code{OR}
operator into a regular expression pattern like this: \code{"(one|two|three)"}.
If \code{keep} is a named list of regular expression patterns, the names of the
list-element should equal the column name where selection should be
applied. This is useful for model objects where \code{model_parameters()}
returns multiple columns with parameter components, like in
\code{\link[=model_parameters.lavaan]{model_parameters.lavaan()}}. Note that the regular expression pattern
should match the parameter names as they are stored in the returned data
frame, which can be different from how they are printed. Inspect the
\verb{$Parameter} column of the parameters table to get the exact parameter
names.}

\item{parameters}{Deprecated, alias for \code{keep}.}

\item{table_wide}{Logical that decides whether the ANOVA table should be in
wide format, i.e. should the numerator and denominator degrees of freedom
be in the same row. Default: \code{FALSE}.}

\item{verbose}{Toggle warnings and messages.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame of indices related to the model's parameters.
}
\description{
Parameters from ANOVAs
}
\note{
For ANOVA-tables from mixed models (i.e. \code{anova(lmer())}), only
partial or adjusted effect sizes can be computed. Note that type 3 ANOVAs
with interactions involved only give sensible and informative results when
covariates are mean-centred and factors are coded with orthogonal contrasts
(such as those produced by \code{contr.sum}, \code{contr.poly}, or
\code{contr.helmert}, but \emph{not} by the default \code{contr.treatment}).
}
\examples{
if (requireNamespace("effectsize", quietly = TRUE)) {
  df <- iris
  df$Sepal.Big <- ifelse(df$Sepal.Width >= 3, "Yes", "No")

  model <- aov(Sepal.Length ~ Sepal.Big, data = df)
  model_parameters(
    model,
    omega_squared = "partial",
    eta_squared = "partial",
    epsilon_squared = "partial"
  )

  model_parameters(
    model,
    omega_squared = "partial",
    eta_squared = "partial",
    ci = .9
  )

  model <- anova(lm(Sepal.Length ~ Sepal.Big, data = df))
  model_parameters(model)
  model_parameters(
    model,
    omega_squared = "partial",
    eta_squared = "partial",
    epsilon_squared = "partial"
  )

  model <- aov(Sepal.Length ~ Sepal.Big + Error(Species), data = df)
  model_parameters(model)

  \dontrun{
    if (require("lme4")) {
      mm <- lmer(Sepal.Length ~ Sepal.Big + Petal.Width + (1 | Species),
        data = df
      )
      model <- anova(mm)

      # simple parameters table
      model_parameters(model)

      # parameters table including effect sizes
      model_parameters(
        model,
        eta_squared = "partial",
        ci = .9,
        df_error = dof_satterthwaite(mm)[2:3]
      )
    }
  }
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_glmmTMB.R, R/simulate_parameters.R
\name{simulate_parameters.glmmTMB}
\alias{simulate_parameters.glmmTMB}
\alias{simulate_parameters}
\alias{simulate_parameters.default}
\title{Simulate Model Parameters}
\usage{
\method{simulate_parameters}{glmmTMB}(
  model,
  iterations = 1000,
  centrality = "median",
  ci = 0.95,
  ci_method = "quantile",
  test = "p-value",
  ...
)

simulate_parameters(model, ...)

\method{simulate_parameters}{default}(
  model,
  iterations = 1000,
  centrality = "median",
  ci = 0.95,
  ci_method = "quantile",
  test = "p-value",
  ...
)
}
\arguments{
\item{model}{Statistical model (no Bayesian models).}

\item{iterations}{The number of draws to simulate/bootstrap.}

\item{centrality}{The point-estimates (centrality indices) to compute.  Character (vector) or list with one or more of these options: \code{"median"}, \code{"mean"}, \code{"MAP"} or \code{"all"}.}

\item{ci}{Value or vector of probability of the CI (between 0 and 1)
to be estimated. Default to \code{.95} (\verb{95\%}).}

\item{ci_method}{The type of index used for Credible Interval. Can be
\code{"HDI"} (default, see \code{\link[bayestestR:hdi]{hdi()}}), \code{"ETI"}
(see \code{\link[bayestestR:eti]{eti()}}), \code{"BCI"} (see
\code{\link[bayestestR:bci]{bci()}}) or \code{"SI"} (see \code{\link[bayestestR:si]{si()}}).}

\item{test}{The indices of effect existence to compute. Character (vector) or
list with one or more of these options: \code{"p_direction"} (or \code{"pd"}),
\code{"rope"}, \code{"p_map"}, \code{"equivalence_test"} (or \code{"equitest"}),
\code{"bayesfactor"} (or \code{"bf"}) or \code{"all"} to compute all tests.
For each "test", the corresponding \pkg{bayestestR} function is called
(e.g. \code{\link[bayestestR:rope]{rope()}} or \code{\link[bayestestR:p_direction]{p_direction()}}) and its results
included in the summary output.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame with simulated parameters.
}
\description{
Compute simulated draws of parameters and their related indices such as Confidence Intervals (CI) and p-values. Simulating parameter draws can be seen as a (computationally faster) alternative to bootstrapping.
}
\details{
\subsection{Technical Details}{
\code{simulate_parameters()} is a computationally faster alternative
to \code{bootstrap_parameters()}. Simulated draws for coefficients are based
on a multivariate normal distribution (\code{MASS::mvrnorm()}) with mean
\code{mu = coef(model)} and variance \code{Sigma = vcov(model)}.
}
\subsection{Models with Zero-Inflation Component}{
For models from packages \pkg{glmmTMB}, \pkg{pscl}, \pkg{GLMMadaptive} and
\pkg{countreg}, the \code{component} argument can be used to specify
which parameters should be simulated. For all other models, parameters
from the conditional component (fixed effects) are simulated. This may
include smooth terms, but not random effects.
}
}
\note{
There is also a \href{https://easystats.github.io/see/articles/parameters.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
library(parameters)

model <- lm(Sepal.Length ~ Species * Petal.Width + Petal.Length, data = iris)
simulate_parameters(model)
\dontrun{
if (require("glmmTMB")) {
  model <- glmmTMB(
    count ~ spp + mined + (1 | site),
    ziformula = ~mined,
    family = poisson(),
    data = Salamanders
  )
  simulate_parameters(model, centrality = "mean")
  simulate_parameters(model, ci = c(.8, .95), component = "zero_inflated")
}
}
}
\references{
Gelman A, Hill J. Data analysis using regression and multilevel/hierarchical models. Cambridge; New York: Cambridge University Press 2007: 140-143
}
\seealso{
\code{\link[=bootstrap_model]{bootstrap_model()}}, \code{\link[=bootstrap_parameters]{bootstrap_parameters()}}, \code{\link[=simulate_model]{simulate_model()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3_p_value.R, R/methods_emmeans.R
\name{p_value}
\alias{p_value}
\alias{p_value.default}
\alias{p_value.emmGrid}
\title{p-values}
\usage{
p_value(model, ...)

\method{p_value}{default}(
  model,
  dof = NULL,
  method = NULL,
  robust = FALSE,
  component = "all",
  verbose = TRUE,
  ...
)

\method{p_value}{emmGrid}(model, ci = 0.95, adjust = "none", ...)
}
\arguments{
\item{model}{A statistical model.}

\item{...}{Arguments passed down to \code{standard_error_robust()} when confidence
intervals or p-values based on robust standard errors should be computed.
Only available for models where \code{method = "robust"} is supported.}

\item{dof}{Number of degrees of freedom to be used when calculating
confidence intervals. If \code{NULL} (default), the degrees of freedom are
retrieved by calling \code{\link[=degrees_of_freedom]{degrees_of_freedom()}} with
approximation method defined in \code{method}. If not \code{NULL}, use this argument
to override the default degrees of freedom used to compute confidence
intervals.}

\item{method}{If \code{"robust"}, and if model is supported by the \pkg{sandwich}
or \pkg{clubSandwich} packages, computes p-values based on robust
covariance matrix estimation.}

\item{robust}{Logical, if \code{TRUE}, computes confidence intervals (or p-values)
based on robust standard errors. See \code{\link[=standard_error_robust]{standard_error_robust()}}.}

\item{component}{Model component for which parameters should be shown. See
the documentation for your object's class in \code{\link[=model_parameters]{model_parameters()}} for
further details.}

\item{verbose}{Toggle warnings and messages.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{adjust}{Character value naming the method used to adjust p-values or
confidence intervals. See \code{?emmeans::summary.emmGrid} for details.}
}
\value{
A data frame with at least two columns: the parameter names and the
p-values. Depending on the model, may also include columns for model
components etc.
}
\description{
This function attempts to return, or compute, p-values of a model's
parameters. See the documentation for your object's class:
\itemize{
\item{\link[=p_value.BFBayesFactor]{Bayesian models} (\pkg{rstanarm}, \pkg{brms}, \pkg{MCMCglmm}, ...)}
\item{\link[=p_value.zeroinfl]{Zero-inflated models} (\code{hurdle}, \code{zeroinfl}, \code{zerocount}, ...)}
\item{\link[=p_value.poissonmfx]{Marginal effects models} (\pkg{mfx})}
\item{\link[=p_value.DirichletRegModel]{Models with special components} (\code{DirichletRegModel}, \code{clm2}, \code{cgam}, ...)}
}
}
\note{
\code{p_value_robust()} resp. \code{p_value(robust = TRUE)}
rely on the \pkg{sandwich} or \pkg{clubSandwich} package (the latter if
\code{vcov_estimation = "CR"} for cluster-robust standard errors) and will
thus only work for those models supported by those packages.
}
\examples{
data(iris)
model <- lm(Petal.Length ~ Sepal.Length + Species, data = iris)
p_value(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{qol_cancer}
\alias{qol_cancer}
\title{Sample data set}
\format{
A data frame with 564 rows and 7 variables:
\describe{
\item{ID}{Patient ID}
\item{QoL}{Quality of Life Score}
\item{time}{Timepoint of measurement}
\item{age}{Age in years}
\item{phq4}{Patients' Health Questionnaire, 4-item version}
\item{hospital}{Hospital ID, where patient was treated}
\item{education}{Patients' educational level}
}
}
\description{
A sample data set with longitudinal data, used in the vignette describing the \code{datawizard::demean()} function. Health-related quality of life from cancer-patients was measured at three time points (pre-surgery, 6 and 12 months after surgery).
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci_ml1.R, R/dof_ml1.R, R/p_value_ml1.R
\name{ci_ml1}
\alias{ci_ml1}
\alias{dof_ml1}
\alias{p_value_ml1}
\title{"m-l-1" approximation for SEs, CIs and p-values}
\usage{
ci_ml1(model, ci = 0.95, robust = FALSE, ...)

dof_ml1(model)

p_value_ml1(model, dof = NULL, robust = FALSE, ...)
}
\arguments{
\item{model}{A mixed model.}

\item{ci}{Confidence Interval (CI) level. Default to \code{0.95} (\verb{95\%}).}

\item{robust}{Logical, if \code{TRUE}, computes confidence intervals (or p-values)
based on robust standard errors. See \code{\link[=standard_error_robust]{standard_error_robust()}}.}

\item{...}{Arguments passed down to \code{\link[=standard_error_robust]{standard_error_robust()}}
when confidence intervals or p-values based on robust standard errors
should be computed.}

\item{dof}{Degrees of Freedom.}
}
\value{
A data frame.
}
\description{
Approximation of degrees of freedom based on a "m-l-1" heuristic
as suggested by Elff et al. (2019).
}
\details{
\subsection{Small Sample Cluster corrected Degrees of Freedom}{
Inferential statistics (like p-values, confidence intervals and
standard errors) may be biased in mixed models when the number of clusters
is small (even if the sample size of level-1 units is high). In such cases
it is recommended to approximate a more accurate number of degrees of freedom
for such inferential statistics (see \cite{Li and Redden 2015}). The
\emph{m-l-1} heuristic is such an approach that uses a t-distribution with
fewer degrees of freedom (\code{dof_ml1()}) to calculate p-values
(\code{p_value_ml1()}) and confidence intervals (\code{ci(method = "ml1")}).
}
\subsection{Degrees of Freedom for Longitudinal Designs (Repeated Measures)}{
In particular for repeated measure designs (longitudinal data analysis),
the \emph{m-l-1} heuristic is likely to be more accurate than simply using the
residual or infinite degrees of freedom, because \code{dof_ml1()} returns
different degrees of freedom for within-cluster and between-cluster effects.
}
\subsection{Limitations of the "m-l-1" Heuristic}{
Note that the "m-l-1" heuristic is not applicable (or at least less accurate)
for complex multilevel designs, e.g. with cross-classified clusters. In such cases,
more accurate approaches like the Kenward-Roger approximation (\code{dof_kenward()})
is recommended. However, the "m-l-1" heuristic also applies to generalized
mixed models, while approaches like Kenward-Roger or Satterthwaite are limited
to linear mixed models only.
}
}
\examples{
\donttest{
if (require("lme4")) {
  model <- lmer(Petal.Length ~ Sepal.Length + (1 | Species), data = iris)
  p_value_ml1(model)
}
}
}
\references{
\itemize{
\item Elff, M.; Heisig, J.P.; Schaeffer, M.; Shikano, S. (2019). Multilevel
Analysis with Few Clusters: Improving Likelihood-based Methods to Provide
Unbiased Estimates and Accurate Inference, British Journal of Political
Science.
\item Li, P., Redden, D. T. (2015). Comparing denominator degrees of freedom
approximations for the generalized linear mixed model in analyzing binary
outcome in small sample cluster-randomized trials. BMC Medical Research
Methodology, 15(1), 38. \doi{10.1186/s12874-015-0026-x}
}
}
\seealso{
\code{dof_ml1()} is a small helper-function to calculate approximated
degrees of freedom of model parameters, based on the "m-l-1" heuristic.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{fish}
\alias{fish}
\title{Sample data set}
\description{
A sample data set, used in tests and some examples.
}
\keyword{data}
