
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
