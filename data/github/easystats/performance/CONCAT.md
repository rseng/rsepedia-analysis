
# performance <img src='man/figures/logo.png' align="right" height="139" />

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03139/status.svg)](https://doi.org/10.21105/joss.03139)
[![downloads](http://cranlogs.r-pkg.org/badges/performance)](https://cran.r-project.org/package=performance)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/performance)](https://cranlogs.r-pkg.org/)
[![status](https://tinyverse.netlify.com/badge/performance)](https://CRAN.R-project.org/package=performance)

***Test if your model is a good model!***

A crucial aspect when building regression models is to evaluate the
quality of modelfit. It is important to investigate how well models fit
to the data and which fit indices to report. Functions to create
diagnostic plots or to compute fit measures do exist, however, mostly
spread over different packages. There is no unique and consistent
approach to assess the model quality for different kind of models.

The primary goal of the **performance** package is to fill this gap and
to provide utilities for computing **indices of model quality** and
**goodness of fit**. These include measures like r-squared (R2), root
mean squared error (RMSE) or intraclass correlation coefficient (ICC) ,
but also functions to check (mixed) models for overdispersion,
zero-inflation, convergence or singularity.

## Installation

[![CRAN](http://www.r-pkg.org/badges/version/performance)](https://cran.r-project.org/package=performance)
[![performance status
badge](https://easystats.r-universe.dev/badges/performance)](https://easystats.r-universe.dev)
[![R
check](https://github.com/easystats/performance/workflows/R-check/badge.svg?branch=master)](https://github.com/easystats/performance/actions)

The *performance* package is available on CRAN, while its latest
development version is available on R-universe (from *rOpenSci*).

| Type        | Source     | Command                                                                       |
|-------------|------------|-------------------------------------------------------------------------------|
| Release     | CRAN       | `install.packages("performance")`                                             |
| Development | R-universe | `install.packages("performance", repos = "https://easystats.r-universe.dev")` |

Once you have downloaded the package, you can then load it using:

``` r
library("performance")
```

## Citation

To cite performance in publications use:

``` r
citation("performance")
#> 
#>   Lüdecke et al., (2021). performance: An R Package for Assessment, Comparison and
#>   Testing of Statistical Models. Journal of Open Source Software, 6(60), 3139.
#>   https://doi.org/10.21105/joss.03139
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {{performance}: An {R} Package for Assessment, Comparison and Testing of Statistical Models},
#>     author = {Daniel Lüdecke and Mattan S. Ben-Shachar and Indrajeet Patil and Philip Waggoner and Dominique Makowski},
#>     year = {2021},
#>     journal = {Journal of Open Source Software},
#>     volume = {6},
#>     number = {60},
#>     pages = {3139},
#>     doi = {10.21105/joss.03139},
#>   }
```

## Documentation

[![Documentation](https://img.shields.io/badge/documentation-performance-orange.svg?colorB=E91E63)](https://easystats.github.io/performance/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-performance-orange.svg?colorB=2196F3)](https://easystats.github.io/performance/reference/index.html)

There is a nice introduction into the package on
[youtube](https://www.youtube.com/watch?v=EPIxQ5i5oxs).

## The *performance* workflow

<img src="man/figures/figure_workflow.png" width="75%" />

### Assessing model quality

#### R-squared

**performance** has a generic `r2()` function, which computes the
r-squared for many different models, including mixed effects and
Bayesian regression models.

`r2()` returns a list containing values related to the “most
appropriate” r-squared for the given model.

``` r
model <- lm(mpg ~ wt + cyl, data = mtcars)
r2(model)
#> # R2 for Linear Regression
#>        R2: 0.830
#>   adj. R2: 0.819

model <- glm(am ~ wt + cyl, data = mtcars, family = binomial)
r2(model)
#> # R2 for Logistic Regression
#>   Tjur's R2: 0.705

library(MASS)
data(housing)
model <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
r2(model)
#>   Nagelkerke's R2: 0.108
```

The different R-squared measures can also be accessed directly via
functions like `r2_bayes()`, `r2_coxsnell()` or `r2_nagelkerke()` (see a
full list of functions
[here](https://easystats.github.io/performance/reference/index.html#section-r-functions)).

For mixed models, the *conditional* and *marginal* R-squared are
returned. The *marginal R-squared* considers only the variance of the
fixed effects and indicates how much of the model’s variance is
explained by the fixed effects part only. The *conditional R-squared*
takes both the fixed and random effects into account and indicates how
much of the model’s variance is explained by the “complete” model.

For frequentist mixed models, `r2()` (resp. `r2_nakagawa()`) computes
the *mean* random effect variances, thus `r2()` is also appropriate for
mixed models with more complex random effects structures, like random
slopes or nested random effects (Johnson 2014; Nakagawa, Johnson, and
Schielzeth 2017).

``` r
set.seed(123)
library(rstanarm)

model <- stan_glmer(Petal.Length ~ Petal.Width + (1 | Species), data = iris, cores = 4)

r2(model)
#> # Bayesian R2 with Compatibility Interval
#> 
#>   Conditional R2: 0.953 (95% CI [0.941, 0.963])
#>      Marginal R2: 0.824 (95% CI [0.713, 0.896])

library(lme4)
model <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
r2(model)
#> # R2 for Mixed Models
#> 
#>   Conditional R2: 0.799
#>      Marginal R2: 0.279
```

#### Intraclass Correlation Coefficient (ICC)

Similar to R-squared, the ICC provides information on the explained
variance and can be interpreted as “the proportion of the variance
explained by the grouping structure in the population” (Hox 2010).

`icc()` calculates the ICC for various mixed model objects, including
`stanreg` models.

``` r
library(lme4)
model <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
icc(model)
#> # Intraclass Correlation Coefficient
#> 
#>      Adjusted ICC: 0.722
#>   Conditional ICC: 0.521
```

…and models of class `brmsfit`.

``` r
library(brms)
set.seed(123)
model <- brm(mpg ~ wt + (1 | cyl) + (1 + wt | gear), data = mtcars)
```

``` r
icc(model)
#> # Intraclass Correlation Coefficient
#> 
#>      Adjusted ICC: 0.930
#>   Conditional ICC: 0.771
```

### Model diagnostics

#### Check for overdispersion

Overdispersion occurs when the observed variance in the data is higher
than the expected variance from the model assumption (for Poisson,
variance roughly equals the mean of an outcome).
`check_overdispersion()` checks if a count model (including mixed
models) is overdispersed or not.

``` r
library(glmmTMB)
data(Salamanders)
model <- glm(count ~ spp + mined, family = poisson, data = Salamanders)
check_overdispersion(model)
#> # Overdispersion test
#> 
#>        dispersion ratio =    2.946
#>   Pearson's Chi-Squared = 1873.710
#>                 p-value =  < 0.001
```

Overdispersion can be fixed by either modelling the dispersion parameter
(not possible with all packages), or by choosing a different
distributional family (like Quasi-Poisson, or negative binomial, see
(Gelman and Hill 2007)).

#### Check for zero-inflation

Zero-inflation (in (Quasi-)Poisson models) is indicated when the amount
of observed zeros is larger than the amount of predicted zeros, so the
model is *underfitting* zeros. In such cases, it is recommended to use
negative binomial or zero-inflated models.

Use `check_zeroinflation()` to check if zero-inflation is present in the
fitted model.

``` r
model <- glm(count ~ spp + mined, family = poisson, data = Salamanders)
check_zeroinflation(model)
#> # Check for zero-inflation
#> 
#>    Observed zeros: 387
#>   Predicted zeros: 298
#>             Ratio: 0.77
```

#### Check for singular model fits

A “singular” model fit means that some dimensions of the
variance-covariance matrix have been estimated as exactly zero. This
often occurs for mixed models with overly complex random effects
structures.

`check_singularity()` checks mixed models (of class `lme`, `merMod`,
`glmmTMB` or `MixMod`) for singularity, and returns `TRUE` if the model
fit is singular.

``` r
library(lme4)
data(sleepstudy)

# prepare data
set.seed(123)
sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
sleepstudy$mysubgrp <- NA
for (i in 1:5) {
    filter_group <- sleepstudy$mygrp == i
    sleepstudy$mysubgrp[filter_group] <- sample(1:30, size = sum(filter_group), replace = TRUE)
}

# fit strange model
model <- lmer(Reaction ~ Days + (1 | mygrp/mysubgrp) + (1 | Subject), data = sleepstudy)

check_singularity(model)
#> [1] TRUE
```

Remedies to cure issues with singular fits can be found
[here](https://easystats.github.io/performance/reference/check_singularity.html).

#### Check for heteroskedasticity

Linear models assume constant error variance (homoskedasticity).

The `check_heteroscedasticity()` functions assess if this assumption has
been violated:

``` r
data(cars)
model <- lm(dist ~ speed, data = cars)

check_heteroscedasticity(model)
#> Warning: Heteroscedasticity (non-constant error variance) detected (p = 0.031).
```

#### Comprehensive visualization of model checks

**performance** provides many functions to check model assumptions, like
`check_collinearity()`, `check_normality()` or
`check_heteroscedasticity()`. To get a comprehensive check, use
`check_model()`.

``` r
# defining a model
model <- lm(mpg ~ wt + am + gear + vs * cyl, data = mtcars)

# checking model assumptions
check_model(model)
```

<img src="man/figures/unnamed-chunk-14-1.png" width="60%" />

### Model performance summaries

`model_performance()` computes indices of model performance for
regression models. Depending on the model object, typical indices might
be r-squared, AIC, BIC, RMSE, ICC or LOOIC.

#### Linear model

``` r
m1 <- lm(mpg ~ wt + cyl, data = mtcars)
model_performance(m1)
#> # Indices of model performance
#> 
#> AIC     |     BIC |    R2 | R2 (adj.) |  RMSE | Sigma
#> -----------------------------------------------------
#> 156.010 | 161.873 | 0.830 |     0.819 | 2.444 | 2.568
```

#### Logistic regression

``` r
m2 <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
model_performance(m2)
#> # Indices of model performance
#> 
#> AIC    |    BIC | Tjur's R2 |  RMSE | Sigma | Log_loss | Score_log | Score_spherical |   PCP
#> --------------------------------------------------------------------------------------------
#> 31.298 | 35.695 |     0.478 | 0.359 | 0.934 |    0.395 |   -14.903 |           0.095 | 0.743
```

#### Linear mixed model

``` r
library(lme4)
m3 <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
model_performance(m3)
#> # Indices of model performance
#> 
#> AIC      |      BIC | R2 (cond.) | R2 (marg.) |   ICC |   RMSE |  Sigma
#> -----------------------------------------------------------------------
#> 1755.628 | 1774.786 |      0.799 |      0.279 | 0.722 | 23.438 | 25.592
```

### Models comparison

The `compare_performance()` function can be used to compare the
performance and quality of several models (including models of different
types).

``` r
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
m4 <- glm(counts ~ outcome + treatment, family = poisson())

compare_performance(m1, m2, m3, m4)
#> # Comparison of Model Performance Indices
#> 
#> Name |   Model |      AIC | AIC weights |      BIC | BIC weights |   RMSE |  Sigma | Score_log | Score_spherical |    R2 | R2 (adj.) | Tjur's R2 | Log_loss |   PCP | R2 (cond.) | R2 (marg.) |   ICC | Nagelkerke's R2
#> -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#> m1   |      lm |  156.010 |     < 0.001 |  161.873 |     < 0.001 |  2.444 |  2.568 |           |                 | 0.830 |     0.819 |           |          |       |            |            |       |                
#> m2   |     glm |   31.298 |       1.000 |   35.695 |       1.000 |  0.359 |  0.934 |   -14.903 |           0.095 |       |           |     0.478 |    0.395 | 0.743 |            |            |       |                
#> m3   | lmerMod | 1755.628 |     < 0.001 | 1774.786 |     < 0.001 | 23.438 | 25.592 |           |                 |       |           |           |          |       |      0.799 |      0.279 | 0.722 |                
#> m4   |     glm |   56.761 |     < 0.001 |   57.747 |     < 0.001 |  3.043 |  1.132 |    -2.598 |           0.324 |       |           |           |          |       |            |            |       |           0.657
```

#### General index of model performance

One can also easily compute and a [**composite
index**](https://easystats.github.io/performance/reference/compare_performance.html#details)
of model performance and sort the models from the best one to the worse.

``` r
compare_performance(m1, m2, m3, m4, rank = TRUE)
#> # Comparison of Model Performance Indices
#> 
#> Name |   Model |   RMSE |  Sigma | AIC weights | BIC weights | Performance-Score
#> --------------------------------------------------------------------------------
#> m2   |     glm |  0.359 |  0.934 |       1.000 |       1.000 |           100.00%
#> m4   |     glm |  3.043 |  1.132 |     < 0.001 |     < 0.001 |            46.89%
#> m1   |      lm |  2.444 |  2.568 |     < 0.001 |     < 0.001 |            46.09%
#> m3   | lmerMod | 23.438 | 25.592 |     < 0.001 |     < 0.001 |             0.00%
```

#### Visualisation of indices of models’ performance

Finally, we provide convenient visualisation (the `see` package must be
installed).

``` r
plot(compare_performance(m1, m2, m4, rank = TRUE))
```

<img src="man/figures/unnamed-chunk-20-1.png" width="100%" />

### Testing models

`test_performance()` (and `test_bf`, its Bayesian sister) carries out
the most relevant and appropriate tests based on the input (for
instance, whether the models are nested or not).

``` r
set.seed(123)
data(iris)

lm1 <- lm(Sepal.Length ~ Species, data = iris)
lm2 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm3 <- lm(Sepal.Length ~ Species * Sepal.Width, data = iris)
lm4 <- lm(Sepal.Length ~ Species * Sepal.Width + Petal.Length + Petal.Width, data = iris)

test_performance(lm1, lm2, lm3, lm4)
#> Name | Model |     BF | Omega2 | p (Omega2) |    LR | p (LR)
#> ------------------------------------------------------------
#> lm1  |    lm |        |        |            |       |       
#> lm2  |    lm | > 1000 |   0.69 |     < .001 | -6.25 | < .001
#> lm3  |    lm | > 1000 |   0.36 |     < .001 | -3.44 | < .001
#> lm4  |    lm | > 1000 |   0.73 |     < .001 | -7.77 | < .001
#> Each model is compared to lm1.

test_bf(lm1, lm2, lm3, lm4)
#> Bayes Factors for Model Comparison
#> 
#>       Model                                                    BF
#> [lm2] Species + Petal.Length                             3.45e+26
#> [lm3] Species * Sepal.Width                              4.69e+07
#> [lm4] Species * Sepal.Width + Petal.Length + Petal.Width 7.58e+29
#> 
#> * Against Denominator: [lm1] Species
#> *   Bayes Factor Type: BIC approximation
```

# Code of Conduct

Please note that the performance project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

# Contributing

We are happy to receive bug reports, suggestions, questions, and (most
of all) contributions to fix problems and add features.

Please follow contributing guidelines mentioned here:

<https://easystats.github.io/performance/CONTRIBUTING.html>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gelman_data_2007" class="csl-entry">

Gelman, Andrew, and Jennifer Hill. 2007. *Data Analysis Using Regression
and Multilevel/Hierarchical Models*. Analytical Methods for Social
Research. Cambridge ; New York: Cambridge University Press.

</div>

<div id="ref-hox_multilevel_2010" class="csl-entry">

Hox, J. J. 2010. *Multilevel Analysis: Techniques and Applications*. 2nd
ed. Quantitative Methodology Series. New York: Routledge.

</div>

<div id="ref-johnson_extension_2014" class="csl-entry">

Johnson, Paul C. D. 2014. “Extension of Nakagawa & Schielzeth’s R2 GLMM
to Random Slopes Models.” Edited by Robert B. O’Hara. *Methods in
Ecology and Evolution* 5 (9): 944–46.
<https://doi.org/10.1111/2041-210X.12225>.

</div>

<div id="ref-nakagawa_coefficient_2017" class="csl-entry">

Nakagawa, Shinichi, Paul C. D. Johnson, and Holger Schielzeth. 2017.
“The Coefficient of Determination R2 and Intra-Class Correlation
Coefficient from Generalized Linear Mixed-Effects Models Revisited and
Expanded.” *Journal of The Royal Society Interface* 14 (134): 20170213.
<https://doi.org/10.1098/rsif.2017.0213>.

</div>

</div>
# performance 0.8.0

## Breaking Changes

* The `ci`-level in `r2()` for Bayesian models now defaults to `0.95`, to be in
  line with the latest changes in the *bayestestR* package.

* S3-method dispatch for `pp_check()` was revised, to avoid problems with the
  _bayesplot_ package, where the generic is located.

## General

* Minor revisions to wording for messages from some of the check-functions.

* `posterior_predictive_check()` and `check_predictions()` were added as aliases
  for `pp_check()`.

## New functions

* `check_multimodal()` and `check_heterogeneity_bias()`. These functions will
  be removed from the _parameters_ packages in the future.

## Changes to functions

* `r2()` for linear models can now compute confidence intervals, via the `ci`
  argument.

## Bug fixes

* Fixed issues in `check_model()` for Bayesian models.

* Fixed issue in `pp_check()` for models with transformed response variables,
  so now predictions and observed response values are on the same (transformed)
  scale.

# performance 0.7.3

## Changes to functions

* `check_outliers()` has new `ci` (or `hdi`, `eti`) method to filter based on
  Confidence/Credible intervals.

* `compare_performance()` now also accepts a list of model objects.

* `performance_roc()` now also works for binomial models from other classes 
  than *glm*.

* Several functions, like `icc()` or `r2_nakagawa()`, now have an 
  `as.data.frame()` method.

* `check_collinearity()` now correctly handles objects from forthcoming *afex* 
  update.

# performance 0.7.2

## New functions

* `performance_mae()` to calculate the mean absolute error.

## Bug fixes

* Fixed issue with `"data length differs from size of matrix"` warnings in
  examples in forthcoming R 4.2.

* Fixed issue in `check_normality()` for models with sample size larger than
  5.000 observations.

* Fixed issue in `check_model()` for *glmmTMB* models.

* Fixed issue in `check_collinearity()` for *glmmTMB* models with zero-inflation,
  where the zero-inflated model was an intercept-only model.

# performance 0.7.1

## New supported models

* Add support for `model_fit` (*tidymodels*).

* `model_performance` supports *kmeans* models.

## General

* Give more informative warning when `r2_bayes()` for *BFBayesFactor* objects
  can't be calculated.

* Several `check_*()` functions now return informative messages for invalid
  model types as input.

* `r2()` supports `mhurdle` (*mhurdle*) models.

* Added `print()` methods for more classes of `r2()`.

* The `performance_roc()` and `performance_accuracy()` functions unfortunately 
  had spelling mistakes in the output columns: *Sensitivity* was called 
  *Sensivity* and *Specificity* was called *Specifity*. We think these are 
  understandable mistakes :-)

## Changes to functions

### `check_model()`

* `check_model()` gains more arguments, to customize plot appearance.

* Added option to detrend QQ/PP plots in `check_model()`.

### `model_performance()`

* The `metrics` argument from `model_performance()` and `compare_performance()`
  gains a `"AICc"` option, to also compute the 2nd order AIC.

* `"R2_adj"` is now an explicit option in the `metrics` argument from
  `model_performance()` and `compare_performance()`.

### Other functions

* The default-method for `r2()` now tries to compute an r-squared for all models
  that have no specific `r2()`-method yet, by using following formula:
  `1-sum((y-y_hat)^2)/sum((y-y_bar)^2))`
  
* The column name `Parameter` in `check_collinearity()` is now more
  appropriately named `Term`.

## Bug fixes

* `test_likelihoodratio()` now correctly sorts models with identical fixed
  effects part, but different other model parts (like zero-inflation).

* Fixed incorrect computation of models from inverse-Gaussian families, or
  Gaussian families fitted with `glm()`.

* Fixed issue in `performance_roc()` for models where outcome was not 0/1 coded.

* Fixed issue in `performance_accuracy()` for logistic regression models when
  `method = "boot"`.

* `cronbachs_alpha()` did not work for `matrix`-objects, as stated in the docs.
  It now does.

# performance 0.7.0

## General

* Roll-back R dependency to R >= 3.4.

## Breaking Changes

* `compare_performance()` doesn't return the models' Bayes Factors, now returned
  by `test_performance()` and `test_bf()`.

## New functions to test or compare models

* `test_vuong()`, to compare models using Vuong's (1989) Test.

* `test_bf()`, to compare models using Bayes factors.

* `test_likelihoodratio()` as an alias for `performance_lrt()`.

* `test_wald()`, as a rough approximation for the LRT.

* `test_performance()`, to run the most relevant and appropriate tests based on
  the input.

## Changes to functions

### `performance_lrt()`

* `performance_lrt()` get an alias `test_likelihoodratio()`.

* Does not return AIC/BIC now (as they are not related to LRT *per se* and can
  be easily obtained with other functions).

* Now contains a column with the difference in degrees of freedom between
  models.

* Fixed column names for consistency.

### `model_performance()`

* Added more diagnostics to models of class `ivreg`.

### Other functions

* Revised computation of `performance_mse()`, to ensure that it's always based
  on response residuals.

* `performance_aic()` is now more robust.

## Bug fixes

* Fixed issue in `icc()` and `variance_decomposition()` for multivariate
  response models, where not all model parts contained random effects.

* Fixed issue in `compare_performance()` with duplicated rows.

* `check_collinearity()` no longer breaks for models with rank deficient model
  matrix, but gives a warning instead.

* Fixed issue in `check_homogeneity()` for `method = "auto"`, which wrongly
  tested the response variable, not the residuals.

* Fixed issue in `check_homogeneity()` for edge cases where predictor had
  non-syntactic names.

# performance 0.6.1

## General

* `check_collinearity()` gains a `verbose` argument, to toggle warnings and
  messages.

## Bug fixes

* Fixed examples, now using suggested packages only conditionally.

# performance 0.6.0

## General

* `model_performance()` now supports `margins`, `gamlss`, `stanmvreg` and
  `semLme`.

## New functions

* `r2_somers()`, to compute Somers' Dxy rank-correlation as R2-measure for
  logistic regression models.

* `display()`, to print output from package-functions into different formats.
  `print_md()` is an alias for `display(format = "markdown")`.

## Changes to functions

### `model_performance()`

* `model_performance()` is now more robust and doesn't fail if an index could
  not be computed. Instead, it returns all indices that were possible to
  calculate.

* `model_performance()` gains a default-method that catches all model objects
  not previously supported. If model object is also not supported by the
  default-method, a warning is given.

* `model_performance()` for metafor-models now includes the degrees of freedom
  for Cochran's Q.

### Other functions

* `performance_mse()` and `performance_rmse()` now always try to return the
  (R)MSE on the response scale.

* `performance_accuracy()` now accepts all types of linear or logistic
  regression models, even if these are not of class `lm` or `glm`.

* `performance_roc()` now accepts all types of logistic regression models, even
  if these are not of class `glm`.

* `r2()` for mixed models and `r2_nakagawa()` gain a `tolerance`-argument, to
  set the tolerance level for singularity checks when computing random effect
  variances for the conditional r-squared.

## Bug fixes

* Fixed issue in `icc()` introduced in the last update that make `lme`-models
  fail.

* Fixed issue in `performance_roc()` for models with factors as response.

# performance 0.5.1

## Breaking changes

* Column names for `model_performance()` and `compare_performance()` were
  changed to be in line with the _easystats_ naming convention: `LOGLOSS` is now
  `Log_loss`, `SCORE_LOG` is `Score_log` and `SCORE_SPHERICAL` is now
  `Score_spherical`.

## New functions
* `r2_posterior()` for Bayesian models to obtain posterior distributions of
  R-squared.

## Changes to functions

* `r2_bayes()` works with Bayesian models from `BayesFactor` ( #143 ).

* `model_performance()` works with Bayesian models from `BayesFactor` ( #150 ).

* `model_performance()` now also includes the residual standard deviation.

* Improved formatting for Bayes factors in `compare_performance()`.

* `compare_performance()` with `rank = TRUE` doesn't use the `BF` values when
  `BIC` are present, to prevent "double-dipping" of the BIC values (#144).

* The `method` argument in `check_homogeneity()` gains a `"levene"` option, to
  use Levene's Test for homogeneity.

## Bug fixes

* Fix bug in `compare_performance()` when `...` arguments were function calls to
  regression objects, instead of direct function calls.

# performance 0.5.0

## General

* `r2()` and `icc()` support `semLME` models (package *smicd*).

* `check_heteroscedasticity()` should now also work with zero-inflated mixed
  models from *glmmTMB* and *GLMMadpative*.

* `check_outliers()` now returns a logical vector. Original numerical vector is
  still accessible via `as.numeric()`.

## New functions

* `pp_check()` to compute posterior predictive checks for frequentist models.

## Bug fixes

* Fixed issue with incorrect labeling of groups from `icc()` when `by_group =
  TRUE`.

* Fixed issue in `check_heteroscedasticity()` for mixed models where sigma could
  not be calculated in a straightforward way.

* Fixed issues in `check_zeroinflation()` for `MASS::glm.nb()`.

* Fixed CRAN check issues.

# performance 0.4.8

## General

* Removed suggested packages that have been removed from CRAN.

## Changes to functions

* `icc()` now also computes a "classical" ICC for `brmsfit` models. The former
  way of calculating an "ICC" for `brmsfit` models is now available as new
  function called `variance_decomposition()`.

## Bug fixes

* Fix issue with new version of *bigutilsr* for `check_outliers()`.

* Fix issue with model order in `performance_lrt()`.

# performance 0.4.7

## General

* Support for models from package *mfx*.

## Changes to functions

* `model_performance.rma()` now includes results from heterogeneity test for
  meta-analysis objects.

* `check_normality()` now also works for mixed models (with the limitation that
  studentized residuals are used).

* `check_normality()` gets an `effects`-argument for mixed models, to check
  random effects for normality.

## Bug fixes

* Fixed issue in `performance_accuracy()` for binomial models when response
  variable had non-numeric factor levels.

* Fixed issues in `performance_roc()`, which printed 1 - AUC instead of AUC.

# performance 0.4.6

## General

* Minor revisions to `model_performance()` to meet changes in *mlogit* package.

* Support for `bayesx` models.

## Changes to functions

* `icc()` gains a `by_group` argument, to compute ICCs per different group
  factors in mixed models with multiple levels or cross-classified design.

* `r2_nakagawa()` gains a `by_group` argument, to compute explained variance at
  different levels (following the variance-reduction approach by Hox 2010).

* `performance_lrt()` now works on *lavaan* objects.

## Bug fixes

* Fix issues in some functions for models with logical dependent variable.

* Fix bug in `check_itemscale()`, which caused multiple computations of skewness
  statistics.

* Fix issues in `r2()` for *gam* models.

# performance 0.4.5

## General

* `model_performance()` and `r2()` now support *rma*-objects from package
  *metafor*, *mlm* and *bife* models.

## Changes to functions

* `compare_performance()` gets a `bayesfactor` argument, to include or exclude
  the Bayes factor for model comparisons in the output.

* Added `r2.aov()`.

## Bug fixes

* Fixed issue in `performance_aic()` for models from package *survey*, which
  returned three different AIC values. Now only the AIC value is returned.

* Fixed issue in `check_collinearity()` for *glmmTMB* models when zero-inflated
  formula only had one predictor.

* Fixed issue in `check_model()` for *lme* models.

* Fixed issue in `check_distribution()` for *brmsfit* models.

* Fixed issue in `check_heteroscedasticity()` for *aov* objects.

* Fixed issues for *lmrob* and *glmrob* objects.

# performance 0.4.4

## General

* Removed `logLik.felm()`, because this method is now implemented in the *lfe*
  package.

* Support for `DirichletRegModel` models.

## New functions

* `check_itemscale()` to describe various measures of internal consistencies for
  scales which were built from several items from a PCA, using
  `parameters::principal_components()`.

* `r2_efron()` to compute Efron's pseudo R2.

## Bug fixes

* Fixed issue in documentation of `performance_score()`.

# performance 0.4.3

## General

* Support for `mixor`, `cpglm` and `cpglmm` models.

## New functions

* `performance_aic()` as a small wrapper that returns the AIC. It is a generic
  function that also works for some models that don't have a AIC method (like
  Tweedie models).

* `performance_lrt()` as a small wrapper around `anova()` to perform a
  Likelihood-Ratio-Test for model comparison.

## Bug fixes

* Fix issues with CRAN checks.

## Changes to functions

* `model_performance()` now calculates AIC for Tweedie models.

# performance 0.4.2

## General

* Support for `bracl`, `brmultinom`, `fixest`, `glmx`, `glmmadmb`, `mclogit`,
  `mmclogit`, `vgam` and `vglm` models.

* `model_performance()` now supports *plm* models.

* `r2()` now supports *complmrob* models.

* `compare_performance()` now gets a `plot()`-method (requires package
  **see**).

## Changes to functions

* `compare_performance()` gets a `rank`-argument, to rank models according to
  their overall model performance.

* `compare_performance()` has a nicer `print()`-method now.

* Verbosity for `compare_performance()` was slightly adjusted.

* `model_performance()`-methods for different objects now also have a
  `verbose`-argument.

## Minor changes

* `check_collinearity()` now no longer returns backticks in row- and column
  names.

## Bug fixes

* Fixed issue in `r2()` for `wbm`-models with cross-level interactions.

* `plot()`-methods for `check_heteroscedasticity()` and `check_homogeneity()`
  now work without requiring to load package *see* before.

* Fixed issues with models of class `rlmerMod`.

# performance 0.4.0

## General

* `performance()` is an alias for `model_performance()`.

## Deprecated and Defunct

* `principal_components()` was removed and re-implemented in the
  **parameters**-package. Please use `parameters::principal_components()` now.

## Changes to functions

* `check_outliers()` now also works on data frames.

* Added more methods to `check_outliers()`.

* `performance_score()` now also works on `stan_lmer()` and `stan_glmer()`
  objects.

* `check_singularity()` now works with models of class *clmm*.

* `r2()` now works with models of class *clmm*, *bigglm* and *biglm*.

* `check_overdispersion()` for mixed models now checks that model family is
  Poisson.

## Bug fixes

* Fixed bug in `compare_performance()` that toggled a warning although models
  were fit from same data.

* Fixed bug in `check_model()` for *glmmTMB* models that occurred when checking
  for outliers.

# performance 0.3.0

## General

* Many `check_*()`-methods now get a `plot()`-method. Package **see** is
  required for plotting.

* `model_performance()` gets a preliminary `print()`-method.

## Breaking changes

* The attribute for the standard error of the Bayesian R2 (`r2_bayes()`) was
  renamed from `std.error` to `SE` to be in line with the naming convention of
  other easystats-packages.

* `compare_performance()` now shows the Bayes factor when all compared models
  are fit from the same data. Previous behaviour was that the BF was shown when
  models were of same class.

## Changes to functions

* `model_performance()` now also works for *lavaan*-objects.

* `check_outliers()` gets a `method`-argument to choose the method for detecting
  outliers. Furthermore, two new methods (Mahalanobis Distance and Invariant
  Coordinate Selection) were implemented.

* `check_model()` now performs more checks for GLM(M)s and other model objects.

* `check_model()` gets a `check`-argument to plot selected checks only.

* `r2_nakagawa()` now returns r-squared for models with singular fit, where no
  random effect variances could be computed. The r-squared then does not take
  random effect variances into account. This behaviour was changed to be in line
  with `MuMIn::r.squaredGLMM()`, which returned a value for models with singular
  fit.

* `check_distribution()` now detects negative binomial and zero-inflated
  distributions. Furthermore, attempt to improve accuracy.

* `check_distribution()` now also accepts a numeric vector as input.

* `compare_performance()` warns if models were not fit from same data.

## New check-functions

* `check_homogeneity()` to check models for homogeneity of variances.

## Bug fixes

* Fixed issues with `compare_performance()` and row-ordering.

* Fixed issue in `check_collinearity()` for zero-inflated models, where the
  zero-inflation component had not enough model terms to calculate
  multicollinearity.

* Fixed issue in some `check_*()` and `performance_*()` functions for models
  with binary outcome, when outcome variable was a factor.

# performance 0.2.0

## General

* `r2()` now works for more regression models.

* `r2_bayes()` now works for multivariate response models.

* `model_performance()` now works for more regression models, and also includes
  the log-loss, proper scoring rules and percentage of correct predictions as
  new metric for models with binary outcome.

## New performance-functions

* `performance_accuracy()`, which calculates the predictive accuracy of linear
  or logistic regression models.

* `performance_logloss()` to compute the log-loss of models with binary outcome.
  The log-loss is a proper scoring function comparable to the `rmse()`.

* `performance_score()` to compute the logarithmic, quadratic and spherical
  proper scoring rules.

* `performance_pcp()` to calculate the percentage of correct predictions for
  models with binary outcome.

* `performance_roc()`, to calculate ROC-curves.

* `performance_aicc()`, to calculate the second-order AIC (AICc).

## New check-functions

* `check_collinearity()` to calculate the variance inflation factor and check
  model predictors for multicollinearity.

* `check_outliers()` to check models for influential observations.

* `check_heteroscedasticity()` to check models for (non-)constant error
  variance.

* `check_normality()` to check models for (non-)normality of residuals.

* `check_autocorrelation()` to check models for auto-correlated residuals.

* `check_distribution()` to classify the distribution of a model-family using
  machine learning.

## New indices-functions

* `r2_mckelvey()` to compute McKelvey and Zavoinas R2 value.

* `r2_zeroinflated()` to compute R2 for zero-inflated (non-mixed) models.

* `r2_xu()` as a crude R2 measure for linear (mixed) models.

## Breaking changes

* `model_performance.stanreg()` and `model_performance.brmsfit()` now only
  return one R2-value and its standard error, instead of different (robust) R2
  measures and credible intervals.

* `error_rate()` is now integrated in the `performance_pcp()`-function.

## Changes to functions

* `model_performance.stanreg()` and `model_performance.brmsfit()` now also
  return the _WAIC_ (widely applicable information criterion).

* `r2_nakagawa()` now calculates the full R2 for mixed models with
  zero-inflation.

* `icc()` now returns `NULL` and no longer stops when no mixed model is
  provided.

* `compare_performance()` now shows the Bayes factor when all compared models
  are of same class.

* Some functions get a `verbose`-argument to show or suppress warnings.

## Bug fixes

* Renamed `r2_coxnell()` to `r2_coxsnell()`.

* Fix issues in `r2_bayes()` and `model_performance()` for ordinal models resp.
  models with cumulative link (#48).

* `compare_performance()` did not sort the `name`-column properly, if the
  columns `class` and `name` were not in the same alphabetical order (#51).

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
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

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
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
---
title: "performance: An R Package for Assessment, Comparison and Testing of Statistical Models"
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
  name: Philip Waggoner
  orcid: 0000-0002-7825-7573
- affiliation: 5
  name: Dominique Makowski
  orcid: 0000-0001-5375-9967
date: "23 March 2021"
bibliography: paper.bib
citation_author: Lüdecke
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
  name: Center for Humans and Machines, Max Planck Institute for Human Development, Berlin, Germany
- index: 4
  name: University of Chicago, USA  
- index: 5
  name: Nanyang Technological University, Singapore
output: rticles::joss_article
csl: apa.csl
journal: JOSS
link-citations: yes
---

# Summary

A crucial part of statistical analysis is evaluating a model's quality and fit, or *performance*. During analysis, especially with regression models, investigating the fit of models to data also often involves selecting the best fitting model amongst many competing models. Upon investigation, fit indices should also be reported both visually and numerically to bring readers in on the investigative effort. 

The *performance* R-package [@rcore] provides utilities for computing measures to assess model quality, many of which are not directly provided by R's *base* or *stats* packages. These include measures like $R^2$, intraclass correlation coefficient (ICC), root mean squared error (RMSE), or functions to check for vexing issues like overdispersion, singularity, or zero-inflation. These functions support a large variety of regression models including generalized linear models, (generalized) mixed-effects models, their Bayesian cousins, and many others.

# Statement of Need

While functions to build and produce diagnostic plots or to compute fit statistics exist, these are located across many packages, which results in a lack of a unique and consistent approach to assess the performance of many types of models. The result is a difficult-to-navigate, unorganized ecosystem of individual packages with different syntax, making it onerous for researchers to locate and use fit indices relevant for their unique purposes. The *performance* package in R fills this gap by offering researchers a suite of intuitive functions with consistent syntax for computing, building, and presenting regression model fit statistics and visualizations.

*performance* is part of the [*easystats*](https://github.com/easystats/performance) ecosystem, which is a collaborative project focused on facilitating simple and intuitive usage of R for statistical analysis [@benshachar2020effectsize; @ludecke2020see; @Lüdecke2020parameters; @makowski2019bayetestR; @Makowski2020correlation].

# Comparison to other Packages

Compared to other packages (e.g., *lmtest* [@lmtest], *MuMIn* [@MuMin], *car* [@car], *broom* [@robinson_broom_2020]), the *performance* package offers functions for checking validity *and* model quality systematically and comprehensively for many regression model objects such as (generalized) linear models, mixed-effects models, and Bayesian models. *performance* also offers functions to compare and test multiple models simultaneously to evaluate the best fitting model to the data.

# Features

Beyond validity and quality checks, *performance* also includes plotting functions via the [*see* package](https://easystats.github.io/see/) [@ludecke2020see].^[A complete overview of plotting functions is available at the *see* website (<https://easystats.github.io/see/articles/performance.html>).]

## Checking Model Assumptions

Inferences made from regression models such as significance tests or interpretation of coefficients require meeting several assumptions, which vary based on the type of model. *performance* offers a collection of functions to check if assumptions are met. To demonstrate the efficiency of the package, we provide examples for a few functions, followed by a broader function that runs a comprehensive suite of checks in a single call.

For example, linear (Gaussian) models assume constant error variance (homoscedasticity). We can use `check_heteroscedasticity()` from *performance* to check if this assumption has been violated.

``` {.r}
data(cars)
model <- lm(dist ~ speed, data = cars)

check_heteroscedasticity(model)

#> Warning: Heteroscedasticity (non-constant error variance) detected (p = 0.031).
```

For another example, Poisson regression models assume equidispersion. Violating this assumption leads to *overdispersion*, which occurs when the observed variance in the data is higher than the expected variance from the model. We can call `check_overdispersion()` to check if overdispersion is an issue.

``` {.r}
library(glmmTMB)
data(Salamanders)
model <- glm(count ~ spp + mined, family = poisson, data = Salamanders)

check_overdispersion(model)

#> # Overdispersion test
#> 
#>        dispersion ratio =    2.946
#>   Pearson's Chi-Squared = 1873.710
#>                 p-value =  < 0.001
#>
#> Overdispersion detected.
```

In addition to tests for checking assumptions, *performance* also provides convenience functions to *visually* assess these assumptions of regression models. *performance*'s visual checks detect the type of model passed to the function call, and return the appropriate visual checks for each model type. At present, there are many supported regression models, such as linear models, linear mixed-effects models or their Bayesian equivalents. Inspect the package documentation for a complete listing.

For example, consider the visual checks from a simple linear regression model.

<!-- TO DO: Regenerate plot once feedback from other has been incorporated -->

``` {.r}
library(see)

model <- lm(mpg ~ wt + am + gear + vs * cyl, 
            data = mtcars)

check_model(model)
```

![](figure1.png)

## Computing Quality Indices of Models

*performance* offers a number of indices to assess the goodness of fit of a model. For example, $R^2$, also known as the coefficient of determination, is a popular statistical measure to gauge the amount of the variance in the dependent variable accounted for by the specified model. The `r2()` function from *performance* computes and returns this index for a variety of regression models. Depending on the model, the returned value may be $R^2$, pseudo-$R^2$, or marginal/adjusted $R^2$.

First, consider a simple linear model.

``` {.r}
model <- lm(mpg ~ wt + cyl, data = mtcars)

r2(model)

#> # R2 for Linear Regression
#>        R2: 0.830
#>   adj. R2: 0.819
```

Next, consider a mixed-effects model.

``` {.r}
library(lme4)
model <- lmer(
  Petal.Length ~ Petal.Width + (1 | Species),
  data = iris
)

r2(model)

#> # R2 for Mixed Models
#> 
#>   Conditional R2: 0.933
#>      Marginal R2: 0.303
```

Similar to $R^2$, the Intraclass Correlation Coefficient (ICC) provides information on the explained variance and can be interpreted as the proportion of the variance explained by the grouping structure in the population [@hox2017multilevel]. The `icc()` function from *performance* computes and returns the ICC for various mixed-effects regression models.

``` {.r}
library(brms)
set.seed(123)
model <- brm(mpg ~ wt + (1 | cyl) + (1 + wt | gear), 
  data = mtcars
  )

icc(model)

#> # Intraclass Correlation Coefficient
#> 
#>      Adjusted ICC: 0.874
#>   Conditional ICC: 0.663
```

Instead of computing and returning individual indices, users can obtain *all* indices from the model by simply passing the fitted model object to `model_performance()`. A list of computed indices is returned, which might include $R^2$, AIC, BIC, RMSE, ICC, LOOIC, etc.

For example, consider a simple linear model.

``` {.r}
m1 <- lm(mpg ~ wt + cyl, data = mtcars)

model_performance(m1)

#> # Indices of model performance
#> 
#> AIC     |     BIC |    R2 | R2 (adj.) |  RMSE | Sigma
#> -----------------------------------------------------
#> 156.010 | 161.873 | 0.830 |     0.819 | 2.444 | 2.568
```

Next, consider a mixed-effects model.

``` {.r}
library(lme4)
m3 <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)

model_performance(m3)

#> # Indices of model performance
#> 
#> AIC      |      BIC | R2 (cond.) | R2 (marg.) |   ICC |   RMSE |  Sigma
#> -----------------------------------------------------------------------
#> 1755.628 | 1774.786 |      0.799 |      0.279 | 0.722 | 23.438 | 25.592
```

## Comparing Multiple Models

For multiple models, users can inspect a table of these indices by calling, [`compare_performance()`](https://easystats.github.io/performance/reference/compare_performance.html).

``` {.r}
data(iris)
lm1 <- lm(Sepal.Length ~ Species, data = iris)
lm2 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm3 <- lm(Sepal.Length ~ Species * Sepal.Width, data = iris)
lm4 <- lm(Sepal.Length ~ Species * Sepal.Width + 
          Petal.Length + Petal.Width, data = iris)

compare_performance(lm1, lm2, lm3, lm4)

#> # Comparison of Model Performance Indices
#> 
#> Name | Model |     AIC |     BIC |    R2 | R2 (adj.) |  RMSE | Sigma
#> --------------------------------------------------------------------
#> lm1  |    lm | 231.452 | 243.494 | 0.619 |     0.614 | 0.510 | 0.515
#> lm2  |    lm | 106.233 | 121.286 | 0.837 |     0.833 | 0.333 | 0.338
#> lm3  |    lm | 187.092 | 208.167 | 0.727 |     0.718 | 0.431 | 0.440
#> lm4  |    lm |  78.797 | 105.892 | 0.871 |     0.865 | 0.296 | 0.305
```

As noted previously, in addition to the returning numeric results, *performance* also offers visualizations of model fit indices.

``` {.r}
library(see)

plot(compare_performance(lm1, lm2, lm3, lm4))
```

![](figure2.png)

## Testing Models

While comparing these indices is often useful, making a decision such as whether to keep or drop a model, can often be difficult as some indices can give conflicting suggestions. Additionally, it may be unclear which index to favour in different contexts. This difficulty is one of the reasons why *tests* are often useful as they facilitate decisions via "significance" indices like *p*-values (in a frequentist framework) or [Bayes Factors](https://easystats.github.io/bayestestR/articles/bayes_factors.html) (in a Bayesian framework).

The generic `test_performance()` function computes the appropriate test(s) based on the supplied input. For instance, the following example shows results from *Vuong's Test* [@vuong_likelihood_1989].

``` {.r}
test_performance(lm1, lm2, lm3, lm4)

#> Name | Model |     BF | Omega2 | p (Omega2) |    LR | p (LR)
#> ------------------------------------------------------------
#> lm1  |    lm |        |        |            |       |       
#> lm2  |    lm | > 1000 |   0.69 |     < .001 | -6.25 | < .001
#> lm3  |    lm | > 1000 |   0.36 |     < .001 | -3.44 | < .001
#> lm4  |    lm | > 1000 |   0.73 |     < .001 | -7.77 | < .001
#> Each model is compared to lm1.
```

*performance* also provides `test_bf()` to compare models via Bayes factors in a Bayesian framework.

``` {.r}
test_bf(lm1, lm2, lm3, lm4)

#> Bayes Factors for Model Comparison
#> 
#>       Model                                                  BF
#> [lm2] Species + Petal.Length                             > 1000
#> [lm3] Species * Sepal.Width                              > 1000
#> [lm4] Species * Sepal.Width + Petal.Length + Petal.Width > 1000
#> 
#> * Against Denominator: [lm1] Species
#> *   Bayes Factor Type: BIC approximation
```

# Licensing and Availability

*performance* is licensed under the GNU General Public License (v3.0), with all source code stored at GitHub (<https://github.com/easystats/performance>), and with a corresponding issue tracker for bug reporting and feature enhancements. In the spirit of honest and open science, we encourage requests, tips for fixes, feature updates, as well as general questions and concerns via direct interaction with contributors and developers.

# Acknowledgments

*performance* is part of the collaborative [*easystats*](https://github.com/easystats/easystats) ecosystem. Thus, we thank the [members of easystats](https://github.com/orgs/easystats/people) as well as the users.

# References
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
and orientation.

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
* Focusing on what is best not just for us as individuals, but for the
  overall community

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

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
[INSERT CONTACT METHOD].
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

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

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
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# Contributing to performance

This outlines how to propose a change to **performance**. 

## Fixing typos

Small typos or grammatical errors in documentation may be edited directly using the GitHub web interface, so long as the changes are made in the _source_ file. If you want to fix typos in the documentation, please edit the related `.R` file in the `R/` folder. Do _not_ edit an `.Rd` file in `man/`.

## Filing an issue

The easiest way to propose a change or new feature is to file an issue. If you've found a
bug, you may also create an associated issue. If possible, try to illustrate your proposal or the bug with a minimal [reproducible example](https://www.tidyverse.org/help/#reprex).

## Pull requests

*  Please create a Git branch for each pull request (PR).
*  Your contributed code should roughly follow the [R style guide](http://style.tidyverse.org), but in particular our [**easystats convention of code-style**](https://github.com/easystats/easystats#convention-of-code-style).
*  performance uses [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html),
for documentation.
*  performance uses [testthat](https://cran.r-project.org/package=testthat). Adding tests to the PR makes it easier for me to merge your PR into the code base.
*  If your PR is a user-visible change, you may add a bullet to the top of `NEWS.md` describing the changes made. You may optionally add your GitHub username, and links to relevant issue(s)/PR(s).

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to
abide by its terms.
## revdepcheck results

We checked 18 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 1 new problems
 * We failed to check 0 packages

Issues with CRAN packages are summarised below.

### New problems
(This reports the first line of each new failure)

* see
  checking examples ... ERROR

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.5 (2021-03-31) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |German_Germany.1252          |
|ctype    |German_Germany.1252          |
|tz       |Europe/Berlin                |
|date     |2021-04-09                   |

# Dependencies

|package     |old    |new     |<U+0394>  |
|:-----------|:------|:-------|:--|
|performance |0.7.0  |0.7.1.1 |*  |
|bayestestR  |0.9.0  |0.9.0   |   |
|insight     |0.13.2 |0.13.2  |   |

# Revdeps

## New problems (1)

|package                |version |error  |warning |note |
|:----------------------|:-------|:------|:-------|:----|
|[see](problems.md#see) |0.6.2   |__+1__ |        |     |

# see

<details>

* Version: 0.6.2
* GitHub: https://github.com/easystats/see
* Source code: https://github.com/cran/see
* Date/Publication: 2021-02-04 18:00:02 UTC
* Number of recursive dependencies: 192

Run `revdep_details(, "see")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in 'see-Ex.R' failed
    The error most likely occurred in:
    
    > ### Name: plot.see_performance_roc
    > ### Title: Plot method for ROC curves
    > ### Aliases: plot.see_performance_roc
    > 
    > ### ** Examples
    > 
    > library(performance)
    ...
      6. |     \-ggplot2:::f(l = layers[[i]], d = data[[i]])
      7. |       \-l$compute_aesthetics(d, plot)
      8. |         \-ggplot2:::f(..., self = self)
      9. |           \-base::lapply(aesthetics, eval_tidy, data = data, env = env)
     10. |             \-rlang:::FUN(X[[i]], ...)
     11. +-Specifity
     12. +-rlang:::`$.rlang_data_pronoun`(.data, Specifity)
     13. | \-rlang:::data_pronoun_get(x, nm)
     14. \-rlang:::abort_data_pronoun(x)
    Execution halted
    ```

*Wow, no problems at all. :)*---
output: 
  github_document:
    toc: false
    fig_width: 10
    fig_height: 6
tags: [r, reports]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: paper/paper.bib
---

# performance <img src='man/figures/logo.png' align="right" height="139" />

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  tidy.opts = list(width.cutoff = 100),
  tidy = TRUE,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  dpi = 300,
  fig.path = "man/figures/",
  comment = "#>"
)
options(
  knitr.kable.NA = "",
  digits = 4,
  width = 100
)

library(performance)
```

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03139/status.svg)](https://doi.org/10.21105/joss.03139)
[![downloads](http://cranlogs.r-pkg.org/badges/performance)](https://cran.r-project.org/package=performance) [![total](https://cranlogs.r-pkg.org/badges/grand-total/performance)](https://cranlogs.r-pkg.org/) [![status](https://tinyverse.netlify.com/badge/performance)](https://CRAN.R-project.org/package=performance)

***Test if your model is a good model!***

A crucial aspect when building regression models is to evaluate the quality of modelfit. It is important to investigate how well models fit to the data and which fit indices to report. Functions to create diagnostic plots or to compute fit measures do exist, however, mostly spread over different packages. There is no unique and consistent approach to assess the model quality for different kind of models.

The primary goal of the **performance** package is to fill this gap and to provide utilities for computing **indices of model quality** and **goodness of fit**. These include measures like r-squared (R2), root mean squared error (RMSE) or intraclass correlation coefficient (ICC) , but also functions to check (mixed) models for overdispersion, zero-inflation, convergence or singularity.

## Installation

[![CRAN](http://www.r-pkg.org/badges/version/performance)](https://cran.r-project.org/package=performance) [![performance status badge](https://easystats.r-universe.dev/badges/performance)](https://easystats.r-universe.dev) [![R check](https://github.com/easystats/performance/workflows/R-check/badge.svg?branch=master)](https://github.com/easystats/performance/actions)

The *performance* package is available on CRAN, while its latest development version is available on R-universe (from _rOpenSci_).

Type | Source | Command
---|---|---
Release | CRAN | `install.packages("performance")`
Development | R-universe | `install.packages("performance", repos = "https://easystats.r-universe.dev")`

Once you have downloaded the package, you can then load it using:

```{r, eval=FALSE}
library("performance")
```

## Citation

To cite performance in publications use:

```{r}
citation("performance")
```

## Documentation

[![Documentation](https://img.shields.io/badge/documentation-performance-orange.svg?colorB=E91E63)](https://easystats.github.io/performance/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-performance-orange.svg?colorB=2196F3)](https://easystats.github.io/performance/reference/index.html) 

There is a nice introduction into the package on [youtube](https://www.youtube.com/watch?v=EPIxQ5i5oxs).

## The *performance* workflow

```{r workflow, echo=FALSE, out.width="75%"}
knitr::include_graphics("man/figures/figure_workflow.png")
```

### Assessing model quality

#### R-squared

**performance** has a generic `r2()` function, which computes the r-squared for
many different models, including mixed effects and Bayesian regression models.

`r2()` returns a list containing values related to the "most appropriate"
r-squared for the given model.

```{r}
model <- lm(mpg ~ wt + cyl, data = mtcars)
r2(model)

model <- glm(am ~ wt + cyl, data = mtcars, family = binomial)
r2(model)

library(MASS)
data(housing)
model <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
r2(model)
```

The different R-squared measures can also be accessed directly via functions like `r2_bayes()`, `r2_coxsnell()` or `r2_nagelkerke()` (see a full list of functions [here](https://easystats.github.io/performance/reference/index.html#section-r-functions)).

For mixed models, the _conditional_ and _marginal_ R-squared are returned. The
_marginal R-squared_ considers only the variance of the fixed effects and
indicates how much of the model's variance is explained by the fixed effects
part only. The _conditional R-squared_ takes both the fixed and random effects
into account and indicates how much of the model's variance is explained by the
"complete" model.

For frequentist mixed models, `r2()` (resp. `r2_nakagawa()`) computes the _mean_
random effect variances, thus `r2()` is also appropriate for mixed models with
more complex random effects structures, like random slopes or nested random
effects [@johnson_extension_2014; @nakagawa_coefficient_2017].

```{r}
set.seed(123)
library(rstanarm)

model <- stan_glmer(
  Petal.Length ~ Petal.Width + (1 | Species),
  data = iris,
  cores = 4
)

r2(model)

library(lme4)
model <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
r2(model)
```

#### Intraclass Correlation Coefficient (ICC)

Similar to R-squared, the ICC provides information on the explained variance and
can be interpreted as "the proportion of the variance explained by the grouping
structure in the population" [@hox_multilevel_2010].

`icc()` calculates the ICC for various mixed model objects, including `stanreg`
models.

```{r}
library(lme4)
model <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
icc(model)
```

...and models of class `brmsfit`.
 
```{r, echo=FALSE}
model <- insight::download_model("brms_mixed_1")
```
```{r, eval=FALSE}
library(brms)
set.seed(123)
model <- brm(mpg ~ wt + (1 | cyl) + (1 + wt | gear), data = mtcars)
```
```{r}
icc(model)
```

### Model diagnostics

#### Check for overdispersion

Overdispersion occurs when the observed variance in the data is higher than the
expected variance from the model assumption (for Poisson, variance roughly
equals the mean of an outcome). `check_overdispersion()` checks if a count model
(including mixed models) is overdispersed or not.

```{r}
library(glmmTMB)
data(Salamanders)
model <- glm(count ~ spp + mined, family = poisson, data = Salamanders)
check_overdispersion(model)
```

Overdispersion can be fixed by either modelling the dispersion parameter (not
possible with all packages), or by choosing a different distributional family
(like Quasi-Poisson, or negative binomial, see [@gelman_data_2007]).

#### Check for zero-inflation

Zero-inflation (in (Quasi-)Poisson models) is indicated when the amount of
observed zeros is larger than the amount of predicted zeros, so the model is
_underfitting_ zeros. In such cases, it is recommended to use negative binomial
or zero-inflated models.

Use `check_zeroinflation()` to check if zero-inflation is present in the fitted model.

```{r}
model <- glm(count ~ spp + mined, family = poisson, data = Salamanders)
check_zeroinflation(model)
```

#### Check for singular model fits

A "singular" model fit means that some dimensions of the variance-covariance
matrix have been estimated as exactly zero. This often occurs for mixed models
with overly complex random effects structures.

`check_singularity()` checks mixed models (of class `lme`, `merMod`, `glmmTMB`
or `MixMod`) for singularity, and returns `TRUE` if the model fit is singular.

```{r}
library(lme4)
data(sleepstudy)

# prepare data
set.seed(123)
sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
sleepstudy$mysubgrp <- NA
for (i in 1:5) {
  filter_group <- sleepstudy$mygrp == i
  sleepstudy$mysubgrp[filter_group] <-
    sample(1:30, size = sum(filter_group), replace = TRUE)
}

# fit strange model
model <- lmer(
  Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
  data = sleepstudy
)

check_singularity(model)
```

Remedies to cure issues with singular fits can be found [here](https://easystats.github.io/performance/reference/check_singularity.html). 

#### Check for heteroskedasticity

Linear models assume constant error variance (homoskedasticity).

The `check_heteroscedasticity()` functions assess if this assumption has been
violated:

```{r}
data(cars)
model <- lm(dist ~ speed, data = cars)

check_heteroscedasticity(model)
```

#### Comprehensive visualization of model checks

**performance** provides many functions to check model assumptions, like
`check_collinearity()`, `check_normality()` or `check_heteroscedasticity()`. To
get a comprehensive check, use `check_model()`.

```{r, fig.height=10, fig.width=8, out.width="60%"}
# defining a model
model <- lm(mpg ~ wt + am + gear + vs * cyl, data = mtcars)

# checking model assumptions
check_model(model)
```

### Model performance summaries

`model_performance()` computes indices of model performance for regression
models. Depending on the model object, typical indices might be r-squared, AIC,
BIC, RMSE, ICC or LOOIC.

#### Linear model

```{r}
m1 <- lm(mpg ~ wt + cyl, data = mtcars)
model_performance(m1)
```

#### Logistic regression

```{r}
m2 <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
model_performance(m2)
```

#### Linear mixed model

```{r}
library(lme4)
m3 <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
model_performance(m3)
```

### Models comparison

The `compare_performance()` function can be used to compare the performance and
quality of several models (including models of different types).

```{r}
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
m4 <- glm(counts ~ outcome + treatment, family = poisson())

compare_performance(m1, m2, m3, m4)
```

#### General index of model performance

One can also easily compute and a [**composite index**](https://easystats.github.io/performance/reference/compare_performance.html#details) of model performance and sort the models from the best one to the worse.

```{r}
compare_performance(m1, m2, m3, m4, rank = TRUE)
```

#### Visualisation of indices of models' performance

Finally, we provide convenient visualisation (the `see` package must be
installed).

```{r}
plot(compare_performance(m1, m2, m4, rank = TRUE))
```

### Testing models

`test_performance()` (and `test_bf`, its Bayesian sister) carries out the most
relevant and appropriate tests based on the input (for instance, whether the
models are nested or not).

```{r}
set.seed(123)
data(iris)

lm1 <- lm(Sepal.Length ~ Species, data = iris)
lm2 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm3 <- lm(Sepal.Length ~ Species * Sepal.Width, data = iris)
lm4 <- lm(Sepal.Length ~ Species * Sepal.Width + Petal.Length + Petal.Width, data = iris)

test_performance(lm1, lm2, lm3, lm4)

test_bf(lm1, lm2, lm3, lm4)
```

# Code of Conduct

Please note that the performance project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

# Contributing

We are happy to receive bug reports, suggestions, questions, and (most of all)
contributions to fix problems and add features.

Please follow contributing guidelines mentioned here:

<https://easystats.github.io/performance/CONTRIBUTING.html>

## References
---
title: "Compare, Test, and Select Models"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, performance, r2]
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{Compare, Test and Select Models}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r , include=FALSE}
library(knitr)
library(performance)
options(knitr.kable.NA = "")
knitr::opts_chunk$set(
  comment = ">",
  message = FALSE,
  warning = FALSE,
  out.width = "100%",
  dpi = 450
)
options(digits = 2)

set.seed(333)

library(rstanarm)
```

# Comparing vs. Testing

Let's imagine that we are interested in explaining the variability in the `Sepal.Length` using 3 different predictors. For that, we can build 3 linear models.

```{r}
model1 <- lm(Sepal.Length ~ Petal.Length, data = iris)
model2 <- lm(Sepal.Length ~ Petal.Width, data = iris)
model3 <- lm(Sepal.Length ~ Sepal.Width, data = iris)
```

## Comparing Indices of Model Performance

The eponymous function from the package,
[`performance()`](https://easystats.github.io/performance/reference/model_performance.html),
can be used to compute different indices of performance (an umbrella term for
indices of fit).

```{r}
library(performance)
library(insight)
library(magrittr) # for pipe operator

# we will use `print_md` function to display a well-formatted table
performance(model1) %>%
  print_md()
```

But for multiple models, one can obtain a useful table to compare these indices
at a glance using the
[`compare_performance()`](https://easystats.github.io/performance/reference/compare_performance.html)
function.

```{r}
compare_performance(model1, model2, model3) %>%
  print_md()
```

If you remember your stats lessons, while comparing different model fits, you
would like to choose a model that has a high $R^2$ value (a measure of how much
variance is explained by predictors), low AIC and BIC values, and low root mean
squared error (RMSE). Based on these criteria, we can immediately see that
`model1` has the best fit.

If you don't like looking at tables, you can also plot them using a plotting method supported in `see` package:

```{r}
library(see)

plot(compare_performance(model1, model2, model3))
```

For more, see: <https://easystats.github.io/see/articles/performance.html>

## Testing Models

While **comparing** these indices is often useful, making a decision (for
instance, which model to keep or drop) can often be hard, as the indices can
give conflicting suggestions. Additionally, it is sometimes unclear which index
to favour in the given context.

This is one of the reason why **tests** are useful, as they facilitate decisions
via (infamous) "significance" indices, like *p*-values (in frequentist
framework) or [Bayes Factors](https://easystats.github.io/bayestestR/articles/bayes_factors.html) (in
Bayesian framework).

```{r}
test_performance(model1, model2, model3) %>%
  print_md()
```

However, these tests also have strong limitations and shortcomings, and cannot
be used as the **one criterion to rule them all**!

You can find more information on how these tests [**here**](https://easystats.github.io/performance/reference/test_performance.html).

## Experimenting

Although we have shown here examples only with simple linear models, we will
highly encourage you to try these functions out with models of your choosing.
For example, these functions work with mixed-effects regression models, Bayesian
regression models, etc.

To demonstrate this, we will run Bayesian versions of linear regression models
we just compared:

```{r}
library(rstanarm)

model1 <- stan_glm(Sepal.Length ~ Petal.Length, data = iris, refresh = 0)
model2 <- stan_glm(Sepal.Length ~ Petal.Width, data = iris, refresh = 0)
model3 <- stan_glm(Sepal.Length ~ Sepal.Width, data = iris, refresh = 0)

compare_performance(model1, model2, model3) %>%
  print_md()
```

Note that, since these are Bayesian regression models, the function
automatically picked up the appropriate indices to
compare!

If you are unfamiliar with some of these, explore more [here](https://easystats.github.io/performance/reference/looic.html).

**Now it's your turn to play!** :)---
title: "R-squared (R2)"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, performance, r2]
vignette: >
  %\VignetteIndexEntry{R-squared (R2)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
editor_options: 
  chunk_output_type: console
---

```{r, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")

knitr::opts_chunk$set(
  comment = ">",
  message = FALSE,
  warning = FALSE,
  out.width = "100%",
  collapse = TRUE,
  strip.white = FALSE,
  dpi = 450
)
options(digits = 2)

if (!requireNamespace("BayesFactor", quietly = TRUE) ||
  !requireNamespace("rstanarm", quietly = TRUE) ||
  !requireNamespace("lme4", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(BayesFactor)
  library(rstanarm)
  library(lme4)
}

set.seed(333)
```

# What is the R2?

The **coefficient of determination**, denoted $R^2$ and pronounced "R squared", typically corresponds the proportion of the variance in the dependent variable (the response) that is *explained* (i.e., predicted) by the independent variables (the predictors).

It is an "absolute" index of *goodness-of-fit*, ranging from 0 to 1 (often expressed in percentage), and can be used for model performance assessment or models comparison.

# Different types of R2

As models become more complex, the computation of an $R^2$ becomes increasingly less straightforward.

Currently, depending on the context of the regression model object, one can choose from the following measures supported in `{performance}`:

- Bayesian $R^2$
- Cox & Snell's $R^2$
- Efron's $R^2$
- Kullback-Leibler $R^2$
- LOO-adjusted $R^2$
- McFadden's $R^2$
- McKelvey & Zavoinas $R^2$
- Nagelkerke's $R^2$
- Nakagawa's $R^2$ for mixed models
- Somers' $D_{xy}$ rank correlation for binary outcomes
- Tjur's $R^2$ - coefficient of determination (D)
- Xu' $R^2$ (Omega-squared)
- $R^2$ for models with zero-inflation

```{r , echo=FALSE, include=FALSE, eval=FALSE}
# DONT INCLUDE FOR NOW AS IT's NOT COMPLETE
d <- data.frame(
  "Model_class" = c("lm", "glm"),
  "r2_simple" = c("X", NA),
  "r2_Tjur" = c(NA, "X")
)
knitr::kable(d)
```

*TO BE COMPLETED.*

Before we begin, let's first load the package.

```{r}
library(performance)
```


# R2 for `lm`

```{r}
m_lm <- lm(wt ~ am * cyl, data = mtcars)

r2(m_lm)
```

# R2 for `glm`

In the context of a generalized linear model (e.g., a logistic model which outcome is binary), $R^2$ doesn't measure the percentage of *"explained variance"*, as this concept doesn't apply. However, the $R^2$s that have been adapted for GLMs have retained the name of "R2", mostly because of the similar properties (the range, the sensitivity, and the interpretation as the amount of explanatory power).

# R2 for Mixed Models

## Marginal vs. Conditional R2

For mixed models, `performance` will return two different $R^2$s:

- The **conditional** $R^2$ 
- The **marginal** $R^2$ 

The marginal $R^2$ considers only the variance of the **fixed effects** (without the random effects), while the conditional $R^2$ takes *both* the **fixed and random effects** into account (i.e., the total model).

```{r}
library(lme4)

# defining a linear mixed-effects model
model <- lmer(Petal.Length ~ Petal.Width + (1 | Species), data = iris)

r2(model)
```

Note that `r2` functions only return the $R^2$ values. We would encourage users to instead always use the `model_performance` function to get a more comprehensive set of indices of model fit. 

```{r}
model_performance(model)
```

But, in the current vignette, we would like to exclusively focus on this family of functions and will only talk about this measure.

# R2 for Bayesian Models

```{r}
library(rstanarm)

model <- stan_glm(mpg ~ wt + cyl, data = mtcars, refresh = 0)
r2(model)
```

As discussed above, for mixed-effects models, there will be two components associated with $R^2$.

```{r}
# defining a Bayesian mixed-effects model
model <- stan_lmer(Petal.Length ~ Petal.Width + (1 | Species), data = iris, refresh = 0)

r2(model)
```

Let's look at another regression analysis carried out with `{BayesFactor}` package.

```{r, eval=utils::packageVersion("BayesFactor") >= package_version("0.9.12-4.3")}
library(BayesFactor)
data(puzzles)

m1 <- anovaBF(extra ~ group + ID,
  data = sleep,
  whichRandom = "ID", progress = FALSE
)

r2(m1)

m2 <- generalTestBF(RT ~ shape * color + ID,
  data = puzzles, whichRandom = "ID",
  neverExclude = "ID", progress = FALSE
)

r2(m2)
```


If you want to know more about these indices, you can check out details and references in the functions that compute them [**here**](https://easystats.github.io/performance/reference/index.html#section-r-functions).

# Interpretation

If you want to know about how to *interpret* these $R^2$ values, see these [**interpretation guidelines**](https://easystats.github.io/effectsize/reference/interpret_r2.html).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.kmeans.R
\name{model_performance.kmeans}
\alias{model_performance.kmeans}
\title{Model summary for k-means clustering}
\usage{
\method{model_performance}{kmeans}(model, verbose = TRUE, ...)
}
\arguments{
\item{model}{Object of type \code{kmeans}.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}
}
\description{
Model summary for k-means clustering
}
\examples{
# a 2-dimensional example
x <- rbind(
  matrix(rnorm(100, sd = 0.3), ncol = 2),
  matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2)
)
colnames(x) <- c("x", "y")
model <- kmeans(x, 2)
model_performance(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_performance.R
\name{compare_performance}
\alias{compare_performance}
\title{Compare performance of different models}
\usage{
compare_performance(..., metrics = "all", rank = FALSE, verbose = TRUE)
}
\arguments{
\item{...}{Multiple model objects (also of different classes).}

\item{metrics}{Can be \code{"all"}, \code{"common"} or a character vector of
metrics to be computed. See related
\code{\link[=model_performance]{documentation()}} of object's class for
details.}

\item{rank}{Logical, if \code{TRUE}, models are ranked according to 'best'
overall model performance. See 'Details'.}

\item{verbose}{Toggle off warnings.}
}
\value{
A data frame (with one row per model) and one column per "index" (see
\code{metrics}).
}
\description{
\code{compare_performance()} computes indices of model
performance for different models at once and hence allows comparison of
indices across models.
}
\details{
\subsection{Model Weights}{
When information criteria (IC) are requested in \code{metrics} (i.e., any of \code{"all"},
\code{"common"}, \code{"AIC"}, \code{"AICc"}, \code{"BIC"}, \code{"WAIC"}, or \code{"LOOIC"}), model
weights based on these criteria are also computed. For all IC except LOOIC,
weights are computed as \code{w = exp(-0.5 * delta_ic) / sum(exp(-0.5 * delta_ic))},
where \code{delta_ic} is the difference between the model's IC value and the
smallest IC value in the model set (Burnham & Anderson, 2002).
For LOOIC, weights are computed as "stacking weights" using
\code{\link[loo:loo_model_weights]{loo::stacking_weights()}}.
}

\subsection{Ranking Models}{
When \code{rank = TRUE}, a new column \code{Performance_Score} is returned.
This score ranges from 0\\% to 100\\%, higher values indicating better model
performance. Note that all score value do not necessarily sum up to 100\\%.
Rather, calculation is based on normalizing all indices (i.e. rescaling
them to a range from 0 to 1), and taking the mean value of all indices for
each model. This is a rather quick heuristic, but might be helpful as
exploratory index.
\cr \cr
In particular when models are of different types (e.g. mixed models,
classical linear models, logistic regression, ...), not all indices will be
computed for each model. In case where an index can't be calculated for a
specific model type, this model gets an \code{NA} value. All indices that
have any \code{NA}s are excluded from calculating the performance score.
\cr \cr
There is a \code{plot()}-method for \code{compare_performance()},
which creates a "spiderweb" plot, where the different indices are
normalized and larger values indicate better model performance.
Hence, points closer to the center indicate worse fit indices
(see \href{https://easystats.github.io/see/articles/performance.html}{online-documentation}
for more details).
}
}
\note{
There is also a \href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
data(iris)
lm1 <- lm(Sepal.Length ~ Species, data = iris)
lm2 <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm3 <- lm(Sepal.Length ~ Species * Petal.Length, data = iris)
compare_performance(lm1, lm2, lm3)
compare_performance(lm1, lm2, lm3, rank = TRUE)

if (require("lme4")) {
  m1 <- lm(mpg ~ wt + cyl, data = mtcars)
  m2 <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
  m3 <- lmer(Petal.Length ~ Sepal.Length + (1 | Species), data = iris)
  compare_performance(m1, m2, m3)
}
}
\references{
Burnham, K. P., & Anderson, D. R. (2002).
\emph{Model selection and multimodel inference: A practical information-theoretic approach} (2nd ed.).
Springer-Verlag. \doi{10.1007/b97636}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_coxsnell.R
\name{r2_coxsnell}
\alias{r2_coxsnell}
\title{Cox & Snell's R2}
\usage{
r2_coxsnell(model, ...)
}
\arguments{
\item{model}{Model with binary outcome.}

\item{...}{Currently not used.}
}
\value{
A named vector with the R2 value.
}
\description{
Calculates the pseudo-R2 value based on the proposal from \cite{Cox & Snell
(1989)}.
}
\details{
This index was proposed by \cite{Cox & Snell (1989, pp. 208-9)} and,
apparently independently, by \cite{Magee (1990)}; but had been suggested
earlier for binary response models by \cite{Maddala (1983)}. However, this
index achieves a maximum of less than 1 for discrete models (i.e. models
whose likelihood is a product of probabilities) which have a maximum of 1,
instead of densities, which can become infinite \cite{(Nagelkerke, 1991)}.
}
\examples{
model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
r2_coxsnell(model)
}
\references{
\itemize{
\item Cox, D. R., Snell, E. J. (1989). Analysis of binary data (Vol. 32).
Monographs on Statistics and Applied Probability.

\item Magee, L. (1990). R 2 measures based on Wald and likelihood ratio
joint significance tests. The American Statistician, 44(3), 250-253.

\item Maddala, G. S. (1986). Limited-dependent and qualitative variables in
econometrics (No. 3). Cambridge university press.

\item Nagelkerke, N. J. (1991). A note on a general definition of the
coefficient of determination. Biometrika, 78(3), 691-692.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2.R
\name{r2}
\alias{r2}
\alias{r2.default}
\alias{r2.merMod}
\title{Compute the model's R2}
\usage{
r2(model, ...)

\method{r2}{default}(model, ci = NULL, ci_method = "analytical", verbose = TRUE, ...)

\method{r2}{merMod}(model, tolerance = 1e-05, ...)
}
\arguments{
\item{model}{A statistical model.}

\item{...}{Arguments passed down to the related r2-methods.}

\item{ci}{Confidence Interval (CI) level. Default is \code{NULL}. Confidence
intervals for R2 can be calculated based on different methods, see
\code{ci_method}.}

\item{ci_method}{Method for constructing the R2 confidence interval.
Options are \code{"analytical"} for sampling-theory-based frequentist
intervals and \code{"bootstrap"} for bootstrap intervals. Analytical intervals
are not available for all models. For Bayesian models, \code{\link[=r2_bayes]{r2_bayes()}} is used.}

\item{verbose}{Logical. Should details about R2 and CI methods be given (\code{TRUE}) or not (\code{FALSE})?}

\item{tolerance}{Tolerance for singularity check of random effects, to decide
whether to compute random effect variances for the conditional r-squared
or not. Indicates up to which value the convergence result is accepted. When
\code{r2_nakagawa()} returns a warning, stating that random effect variances
can't be computed (and thus, the conditional r-squared is \code{NA}),
decrease the tolerance-level. See also \code{\link[=check_singularity]{check_singularity()}}.}
}
\value{
Returns a list containing values related to the most appropriate R2
for the given model (or \code{NULL} if no R2 could be extracted). See the
list below:
\itemize{
\item Logistic models: \link[=r2_tjur]{Tjur's R2}
\item General linear models: \link[=r2_nagelkerke]{Nagelkerke's R2}
\item Multinomial Logit: \link[=r2_mcfadden]{McFadden's R2}
\item Models with zero-inflation: \link[=r2_zeroinflated]{R2 for zero-inflated models}
\item Mixed models: \link[=r2_nakagawa]{Nakagawa's R2}
\item Bayesian models: \link[=r2_bayes]{R2 bayes}
}
}
\description{
Calculate the R2, also known as the coefficient of
determination, value for different model objects. Depending on the model,
R2, pseudo-R2, or marginal / adjusted R2 values are returned.
}
\note{
If there is no \code{r2()}-method defined for the given model class,
\code{r2()} tries to return a "generic r2 value, calculated as following:
\verb{1-sum((y-y_hat)^2)/sum((y-y_bar)^2))}
}
\examples{
model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
r2(model)

if (require("lme4")) {
  model <- lmer(Sepal.Length ~ Petal.Length + (1 | Species), data = iris)
  r2(model)
}
}
\seealso{
\code{\link[=r2_bayes]{r2_bayes()}}, \code{\link[=r2_coxsnell]{r2_coxsnell()}}, \code{\link[=r2_kullback]{r2_kullback()}},
\code{\link[=r2_loo]{r2_loo()}}, \code{\link[=r2_mcfadden]{r2_mcfadden()}}, \code{\link[=r2_nagelkerke]{r2_nagelkerke()}},
\code{\link[=r2_nakagawa]{r2_nakagawa()}}, \code{\link[=r2_tjur]{r2_tjur()}}, \code{\link[=r2_xu]{r2_xu()}} and
\code{\link[=r2_zeroinflated]{r2_zeroinflated()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_accuracy.R
\name{performance_accuracy}
\alias{performance_accuracy}
\title{Accuracy of predictions from model fit}
\usage{
performance_accuracy(
  model,
  method = c("cv", "boot"),
  k = 5,
  n = 1000,
  verbose = TRUE
)
}
\arguments{
\item{model}{A linear or logistic regression model. A mixed-effects model is
also accepted.}

\item{method}{Character string, indicating whether cross-validation
(\code{method = "cv"}) or bootstrapping (\code{method = "boot"}) is used to
compute the accuracy values.}

\item{k}{The number of folds for the k-fold cross-validation.}

\item{n}{Number of bootstrap-samples.}

\item{verbose}{Toggle warnings.}
}
\value{
A list with three values: The \code{Accuracy} of the model
predictions, i.e. the proportion of accurately predicted values from the
model, its standard error, \code{SE}, and the \code{Method} used to compute
the accuracy.
}
\description{
This function calculates the predictive accuracy of linear
or logistic regression models.
}
\details{
For linear models, the accuracy is the correlation coefficient
between the actual and the predicted value of the outcome. For
logistic regression models, the accuracy corresponds to the
AUC-value, calculated with the \code{bayestestR::auc()}-function.
\cr \cr
The accuracy is the mean value of multiple correlation resp.
AUC-values, which are either computed with cross-validation
or non-parametric bootstrapping (see argument \code{method}).
The standard error is the standard deviation of the computed
correlation resp. AUC-values.
}
\examples{
model <- lm(mpg ~ wt + cyl, data = mtcars)
performance_accuracy(model)

model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
performance_accuracy(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_mae.R
\name{performance_mae}
\alias{performance_mae}
\alias{mae}
\title{Mean Absolute Error of Models}
\usage{
performance_mae(model, ...)

mae(model, ...)
}
\arguments{
\item{model}{A model.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
Numeric, the mean absolute error of \code{model}.
}
\description{
Compute mean absolute error of models.
}
\examples{
data(mtcars)
m <- lm(mpg ~ hp + gear, data = mtcars)
performance_mae(m)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_model.R
\name{check_model}
\alias{check_model}
\alias{check_model.default}
\title{Visual check of model assumptions}
\usage{
check_model(x, ...)

\method{check_model}{default}(
  x,
  dot_size = 2,
  line_size = 0.8,
  panel = TRUE,
  check = "all",
  alpha = 0.2,
  dot_alpha = 0.8,
  colors = c("#3aaf85", "#1b6ca8", "#cd201f"),
  theme = "see::theme_lucid",
  detrend = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{A model object.}

\item{...}{Currently not used.}

\item{dot_size, line_size}{Size of line and dot-geoms.}

\item{panel}{Logical, if \code{TRUE}, plots are arranged as panels; else,
single plots for each diagnostic are returned.}

\item{check}{Character vector, indicating which checks for should be performed
and plotted. May be one or more of
\verb{"all", "vif", "qq", "normality", "linearity", "ncv", "homogeneity", "outliers", "reqq"}.
\code{"reqq"} is a QQ-plot for random effects and only available for mixed
models.
\code{"ncv"} is an alias for \code{"linearity"}, and checks for non-constant
variance, i.e. for heteroscedasticity, as well as the linear relationship.
By default, all possible checks are performed and plotted.}

\item{alpha, dot_alpha}{The alpha level of the confidence bands and dot-geoms.
Scalar from 0 to 1.}

\item{colors}{Character vector with color codes (hex-format). Must be of
length 3. First color is usually used for reference lines, second color
for dots, and third color for outliers or extreme values.}

\item{theme}{String, indicating the name of the plot-theme. Must be in the
format \code{"package::theme_name"} (e.g. \code{"ggplot2::theme_minimal"}).}

\item{detrend}{Should QQ/PP plots be detrended?}

\item{verbose}{Toggle off warnings.}
}
\value{
The data frame that is used for plotting.
}
\description{
Visual check of model various assumptions (normality of residuals, normality
of random effects, linear relationship, homogeneity of variance,
multicollinearity).
}
\details{
For Bayesian models from packages \strong{rstanarm} or \strong{brms},
models will be "converted" to their frequentist counterpart, using
\href{https://easystats.github.io/bayestestR/reference/convert_bayesian_as_frequentist.html}{\code{bayestestR::bayesian_as_frequentist}}.
A more advanced model-check for Bayesian models will be implemented at a
later stage.
}
\note{
This function just prepares the data for plotting. To create the plots,
\CRANpkg{see} needs to be installed. Furthermore, this function suppresses
all possible warnings. In case you observe suspicious plots, please refer
to the dedicated functions (like \code{check_collinearity()},
\code{check_normality()} etc.) to get informative messages and warnings.
}
\section{Linearity Assumption}{

The plot \strong{Linearity} checks the assumption of linear relationship.
However, the spread of dots also indicate possible heteroscedasticity (i.e.
non-constant variance); hence, the alias \code{"ncv"} for this plot.
\strong{Some caution is needed} when interpreting these plots. Although these
plots are helpful to check model assumptions, they do not necessarily
indicate so-called "lack of fit", e.g. missed non-linear relationships or
interactions. Thus, it is always recommended to also look at
\href{https://strengejacke.github.io/ggeffects/articles/introduction_partial_residuals.html}{effect plots, including partial residuals}.
}

\section{Residuals for (Generalized) Linear Models}{

Plots that check the normality of residuals (QQ-plot) or the homogeneity of
variance use standardized Pearson's residuals for generalized linear models,
and standardized residuals for linear models. The plots for the normality of
residuals (with overlayed normal curve) and for the linearity assumption use
the default residuals for \code{lm} and \code{glm} (which are deviance
residuals for \code{glm}).
}

\examples{
\dontrun{
m <- lm(mpg ~ wt + cyl + gear + disp, data = mtcars)
check_model(m)

if (require("lme4")) {
  m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  check_model(m, panel = FALSE)
}

if (require("rstanarm")) {
  m <- stan_glm(mpg ~ wt + gear, data = mtcars, chains = 2, iter = 200)
  check_model(m)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.R
\name{model_performance}
\alias{model_performance}
\alias{performance}
\title{Model Performance}
\usage{
model_performance(model, ...)

performance(model, ...)
}
\arguments{
\item{model}{Statistical model.}

\item{...}{Arguments passed to or from other methods, resp. for
\code{compare_performance()}, one or multiple model objects (also of
different classes).}
}
\value{
A data frame (with one row) and one column per "index" (see \code{metrics}).
}
\description{
See the documentation for your object's class:
\itemize{
\item \link[=model_performance.lm]{Frequentist Regressions}
\item \link[=model_performance.ivreg]{Instrumental Variables Regressions}
\item \link[=model_performance.merMod]{Mixed models}
\item \link[=model_performance.stanreg]{Bayesian models}
\item \link[=model_performance.lavaan]{CFA / SEM lavaan models}
\item \link[=model_performance.rma]{Meta-analysis models}
}
}
\examples{
model <- lm(mpg ~ wt + cyl, data = mtcars)
model_performance(model)

model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
model_performance(model)
}
\seealso{
\code{\link[=compare_performance]{compare_performance()}} to compare performance of many different models.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_normality.R
\name{check_normality}
\alias{check_normality}
\alias{check_normality.merMod}
\title{Check model for (non-)normality of residuals.}
\usage{
check_normality(x, ...)

\method{check_normality}{merMod}(x, effects = c("fixed", "random"), ...)
}
\arguments{
\item{x}{A model object.}

\item{...}{Currently not used.}

\item{effects}{Should normality for residuals (\code{"fixed"}) or random
effects (\code{"random"}) be tested? Only applies to mixed-effects models.
May be abbreviated.}
}
\value{
Invisibly returns the p-value of the test statistics. A p-value
< 0.05 indicates a significant deviation from normal distribution
}
\description{
Check model for (non-)normality of residuals.
}
\details{
\code{check_normality()} calls \code{stats::shapiro.test}
and checks the standardized residuals (or Studentized residuals for mixed
models) for normal distribution. Note that this formal test almost always
yields significant results for the distribution of residuals and visual
inspection (e.g. Q-Q plots) are preferable.
}
\note{
For mixed-effects models, studentized residuals, and \emph{not}
standardized residuals, are used for the test. There is also a
\href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method}
implemented in the
\href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
m <<- lm(mpg ~ wt + cyl + gear + disp, data = mtcars)
check_normality(m)

# plot results
if (require("see")) {
  x <- check_normality(m)
  plot(x)
}
\dontrun{
# QQ-plot
plot(check_normality(m), type = "qq")

# PP-plot
plot(check_normality(m), type = "pp")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_heteroscedasticity.R
\name{check_heteroscedasticity}
\alias{check_heteroscedasticity}
\alias{check_heteroskedasticity}
\title{Check model for (non-)constant error variance}
\usage{
check_heteroscedasticity(x, ...)

check_heteroskedasticity(x, ...)
}
\arguments{
\item{x}{A model object.}

\item{...}{Currently not used.}
}
\value{
Invisibly returns the p-value of the test statistics. A p-value <
0.05 indicates a non-constant variance (heteroskedasticity).
}
\description{
Significance testing for linear regression models assumes that
the model errors (or residuals) have constant variance. If this assumption
is violated the p-values from the model are no longer reliable.
}
\details{
This test of the hypothesis of (non-)constant error is also called
\emph{Breusch-Pagan test} (\cite{1979}).
}
\note{
There is also a \href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
m <<- lm(mpg ~ wt + cyl + gear + disp, data = mtcars)
check_heteroscedasticity(m)

# plot results
if (require("see")) {
  x <- check_heteroscedasticity(m)
  plot(x)
}
}
\references{
Breusch, T. S., and Pagan, A. R. (1979) A simple test for heteroscedasticity and random coefficient variation. Econometrica 47, 1287-1294.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binned_residuals.R
\name{binned_residuals}
\alias{binned_residuals}
\title{Binned residuals for binomial logistic regression}
\usage{
binned_residuals(model, term = NULL, n_bins = NULL, ...)
}
\arguments{
\item{model}{A \code{glm}-object with \emph{binomial}-family.}

\item{term}{Name of independent variable from \code{x}. If not \code{NULL},
average residuals for the categories of \code{term} are plotted; else,
average residuals for the estimated probabilities of the response are
plotted.}

\item{n_bins}{Numeric, the number of bins to divide the data. If
\code{n_bins = NULL}, the square root of the number of observations is
taken.}

\item{...}{Further argument like \code{size} (for point-size) or
\code{color} (for point-colors).}
}
\value{
A data frame representing the data that is mapped in the accompanying
plot. In case all residuals are inside the error bounds, points are black.
If some of the residuals are outside the error bounds (indicated by the
grey-shaded area), blue points indicate residuals that are OK, while red
points indicate model under- or over-fitting for the relevant range of
estimated probabilities.
}
\description{
Check model quality of binomial logistic regression models.
}
\details{
Binned residual plots are achieved by \dQuote{dividing the data into
categories (bins) based on their fitted values, and then plotting
the average residual versus the average fitted value for each bin.}
\cite{(Gelman, Hill 2007: 97)}. If the model were true, one would
expect about 95\\% of the residuals to fall inside the error bounds.
\cr \cr
If \code{term} is not \code{NULL}, one can compare the residuals in
relation to a specific model predictor. This may be helpful to check if a
term would fit better when transformed, e.g. a rising and falling pattern
of residuals along the x-axis is a signal to consider taking the logarithm
of the predictor (cf. Gelman and Hill 2007, pp. 97-98).
}
\note{
Since \code{binned_residuals()} returns a data frame, the default
action for the result is \emph{printing}. However, the \code{print()}-method for
\code{binned_residuals()} actually creates a plot. For further
modifications of the plot, use \code{print()} and add ggplot-layers to the
return values, e.g. \code{print(binned_residuals(model)) + see::scale_color_pizza()}.
}
\examples{
if (require("see")) {
  # creating a model
  model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")

  # this will automatically plot the results
  (result <- binned_residuals(model))

  # if you assign results to an object, you can also look at the dataframe
  as.data.frame(result)
}
}
\references{
Gelman, A., & Hill, J. (2007). Data analysis using regression and
multilevel/hierarchical models. Cambridge; New York: Cambridge University
Press.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_roc.R
\name{performance_roc}
\alias{performance_roc}
\title{Simple ROC curve}
\usage{
performance_roc(x, ..., predictions, new_data)
}
\arguments{
\item{x}{A numeric vector, representing the outcome (0/1), or a model with
binomial outcome.}

\item{...}{One or more models with binomial outcome. In this case,
\code{new_data} is ignored.}

\item{predictions}{If \code{x} is numeric, a numeric vector of same length
as \code{x}, representing the actual predicted values.}

\item{new_data}{If \code{x} is a model, a data frame that is passed to
\code{predict()} as \code{newdata}-argument. If \code{NULL}, the ROC for
the full model is calculated.}
}
\value{
A data frame with three columns, the x/y-coordinate pairs for the ROC
curve (\code{Sensitivity} and \code{Specificity}), and a column with the
model name.
}
\description{
This function calculates a simple ROC curves of x/y coordinates
based on response and predictions of a binomial model.
}
\note{
There is also a \href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
library(bayestestR)
data(iris)

set.seed(123)
iris$y <- rbinom(nrow(iris), size = 1, .3)
folds <- sample(nrow(iris), size = nrow(iris) / 8, replace = FALSE)
test_data <- iris[folds, ]
train_data <- iris[-folds, ]

model <- glm(y ~ Sepal.Length + Sepal.Width, data = train_data, family = "binomial")
as.data.frame(performance_roc(model, new_data = test_data))

roc <- performance_roc(model, new_data = test_data)
area_under_curve(roc$Specificity, roc$Sensitivity)

m1 <- glm(y ~ Sepal.Length + Sepal.Width, data = iris, family = "binomial")
m2 <- glm(y ~ Sepal.Length + Petal.Width, data = iris, family = "binomial")
m3 <- glm(y ~ Sepal.Length + Species, data = iris, family = "binomial")
performance_roc(m1, m2, m3)

# if you have `see` package installed, you can also plot comparison of
# ROC curves for different models
if (require("see")) plot(performance_roc(m1, m2, m3))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R, R/print_md.R
\name{display.performance_model}
\alias{display.performance_model}
\alias{print_md.performance_model}
\title{Print tables in different output formats}
\usage{
\method{display}{performance_model}(object, format = "markdown", digits = 2, caption = NULL, ...)

\method{print_md}{performance_model}(x, digits = 2, caption = "Indices of model performance", ...)
}
\arguments{
\item{object, x}{An object returned by \code{\link[=model_performance]{model_performance()}}
or \code{\link[=compare_performance]{compare_performance()}}.
or its summary.}

\item{format}{String, indicating the output format. Currently, only
\code{"markdown"} is supported.}

\item{digits}{Number of decimal places.}

\item{caption}{Table caption as string. If \code{NULL}, no table caption is printed.}

\item{...}{Currently not used.}
}
\value{
A character vector. If \code{format = "markdown"}, the return value
will be a character vector in markdown-table format.
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
mp <- model_performance(model)
display(mp)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_homogeneity.R
\name{check_homogeneity}
\alias{check_homogeneity}
\alias{check_homogeneity.afex_aov}
\title{Check model for homogeneity of variances}
\usage{
check_homogeneity(x, method = c("bartlett", "fligner", "levene", "auto"), ...)

\method{check_homogeneity}{afex_aov}(x, method = "levene", ...)
}
\arguments{
\item{x}{A linear model or an ANOVA object.}

\item{method}{Name of the method (underlying test) that should be performed
to check the homogeneity of variances. May either be \code{"levene"} for
Levene's Test for Homogeneity of Variance, \code{"bartlett"} for the
Bartlett test (assuming normal distributed samples or groups),
\code{"fligner"} for the Fligner-Killeen test (rank-based, non-parametric
test), or \code{"auto"}. In the latter case, Bartlett test is used if the
model response is normal distributed, else Fligner-Killeen test is used.}

\item{...}{Arguments passed down to \code{car::leveneTest()}.}
}
\value{
Invisibly returns the p-value of the test statistics. A p-value <
0.05 indicates a significant difference in the variance between the groups.
}
\description{
Check model for homogeneity of variances between groups described
by independent variables in a model.
}
\note{
There is also a \href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
model <<- lm(len ~ supp + dose, data = ToothGrowth)
check_homogeneity(model)

# plot results
if (require("see")) {
  result <- check_homogeneity(model)
  plot(result)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_check.R
\name{check_predictions}
\alias{check_predictions}
\alias{posterior_predictive_check}
\alias{check_posterior_predictions}
\title{Posterior predictive checks}
\usage{
check_predictions(
  object,
  iterations = 50,
  check_range = FALSE,
  re_formula = NULL,
  ...
)

posterior_predictive_check(
  object,
  iterations = 50,
  check_range = FALSE,
  re_formula = NULL,
  ...
)

check_posterior_predictions(
  object,
  iterations = 50,
  check_range = FALSE,
  re_formula = NULL,
  ...
)
}
\arguments{
\item{object}{A statistical model.}

\item{iterations}{The number of draws to simulate/bootstrap.}

\item{check_range}{Logical, if \code{TRUE}, includes a plot with the minimum
value of the original response against the minimum values of the replicated
responses, and the same for the maximum value. This plot helps judging whether
the variation in the original data is captured by the model or not
(\cite{Gelman et al. 2020, pp.163}). The minimum and maximum values of \code{y} should
be inside the range of the related minimum and maximum values of \code{yrep}.}

\item{re_formula}{Formula containing group-level effects (random effects) to
be considered in the simulated data. If \code{NULL} (default), condition
on all random effects. If \code{NA} or \code{~0}, condition on no random
effects. See \code{simulate()} in \strong{lme4}.}

\item{...}{Passed down to \code{simulate()}.}
}
\value{
A data frame of simulated responses and the original response vector.
}
\description{
Posterior predictive checks mean \dQuote{simulating replicated data
under the fitted model and then comparing these to the observed data}
\cite{(Gelman and Hill, 2007, p. 158)}. Posterior predictive checks
can be used to \dQuote{look for systematic discrepancies between real and
simulated data} \cite{(Gelman et al. 2014, p. 169)}.

\pkg{performance} provides posterior predictive check methods for a variety
of frequentist models (e.g., \code{lm}, \code{merMod}, \code{glmmTMB}, ...). For Bayesian
models, the model is passed to \code{\link[bayesplot:pp_check]{bayesplot::pp_check()}}.
}
\details{
An example how posterior predictive checks can also be used for model
comparison is Figure 6 from \cite{Gabry et al. 2019, Figure 6}.
\cr
\if{html}{\cr \figure{pp_check.png}{options: width="90\%" alt="Posterior Predictive Check"} \cr}
The model shown in the right panel (b) can simulate new data that are more
similar to the observed outcome than the model in the left panel (a). Thus,
model (b) is likely to be preferred over model (a).
}
\note{
Every model object that has a \code{simulate()}-method should work with
\code{check_predictions()}. On R 3.6.0 and higher, if \pkg{bayesplot}
(or a package that imports \pkg{bayesplot} such as \pkg{rstanarm} or \pkg{brms})
is loaded, \code{pp_check()} is also available as an alias for \code{check_predictions()}.
}
\examples{
library(performance)
model <- lm(mpg ~ disp, data = mtcars)
if (require("see")) {
  check_predictions(model)
}
}
\references{
\itemize{
\item Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A. (2019). Visualization in Bayesian workflow. Journal of the Royal Statistical Society: Series A (Statistics in Society), 182(2), 389–402. https://doi.org/10.1111/rssa.12378
\item Gelman, A., & Hill, J. (2007). Data analysis using regression and multilevel/hierarchical models. Cambridge; New York: Cambridge University Press.
\item Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2014). Bayesian data analysis. (Third edition). CRC Press.
\item Gelman, A., Hill, J., & Vehtari, A. (2020). Regression and Other Stories. Cambridge University Press.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_singularity.R
\name{check_singularity}
\alias{check_singularity}
\title{Check mixed models for boundary fits}
\usage{
check_singularity(x, tolerance = 1e-05, ...)
}
\arguments{
\item{x}{A mixed model.}

\item{tolerance}{Indicates up to which value the convergence result is
accepted. The larger \code{tolerance} is, the stricter the test
will be.}

\item{...}{Currently not used.}
}
\value{
\code{TRUE} if the model fit is singular.
}
\description{
Check mixed models for boundary fits.
}
\details{
If a model is "singular", this means that some dimensions of the
variance-covariance matrix have been estimated as exactly zero. This
often occurs for mixed models with complex random effects structures.
\cr \cr
\dQuote{While singular models are statistically well defined (it is
theoretically sensible for the true maximum likelihood estimate to
correspond to a singular fit), there are real concerns that (1) singular
fits correspond to overfitted models that may have poor power; (2) chances
of numerical problems and mis-convergence are higher for singular models
(e.g. it may be computationally difficult to compute profile confidence
intervals for such models); (3) standard inferential procedures such as
Wald statistics and likelihood ratio tests may be inappropriate.}
(\cite{lme4 Reference Manual})
\cr \cr
There is no gold-standard about how to deal with singularity and which
random-effects specification to choose. Beside using fully Bayesian methods
(with informative priors), proposals in a frequentist framework are:
\itemize{
\item avoid fitting overly complex models, such that the
variance-covariance matrices can be estimated precisely enough
(\cite{Matuschek et al. 2017})
\item use some form of model selection to choose a model that balances
predictive accuracy and overfitting/type I error (\cite{Bates et al. 2015},
\cite{Matuschek et al. 2017})
\item \dQuote{keep it maximal}, i.e. fit the most complex model consistent
with the experimental design, removing only terms required to allow a
non-singular fit (\cite{Barr et al. 2013})
}
Note the different meaning between singularity and convergence: singularity
indicates an issue with the "true" best estimate, i.e. whether the maximum
likelihood estimation for the variance-covariance matrix of the random
effects is positive definite or only semi-definite. Convergence is a
question of whether we can assume that the numerical optimization has
worked correctly or not.
}
\examples{
if (require("lme4")) {
  data(sleepstudy)
  set.seed(123)
  sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
  sleepstudy$mysubgrp <- NA
  for (i in 1:5) {
    filter_group <- sleepstudy$mygrp == i
    sleepstudy$mysubgrp[filter_group] <-
      sample(1:30, size = sum(filter_group), replace = TRUE)
  }

  model <- lmer(
    Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
    data = sleepstudy
  )

  check_singularity(model)
}
}
\references{
\itemize{
\item Bates D, Kliegl R, Vasishth S, Baayen H. Parsimonious Mixed Models.
arXiv:1506.04967, June 2015.

\item Barr DJ, Levy R, Scheepers C, Tily HJ. Random effects structure for
confirmatory hypothesis testing: Keep it maximal. Journal of Memory and
Language, 68(3):255-278, April 2013.

\item Matuschek H, Kliegl R, Vasishth S, Baayen H, Bates D. Balancing type
I error and power in linear mixed models. Journal of Memory and Language,
94:305-315, 2017.

\item lme4 Reference Manual, \url{https://cran.r-project.org/package=lme4}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.mixed.R
\name{model_performance.merMod}
\alias{model_performance.merMod}
\title{Performance of Mixed Models}
\usage{
\method{model_performance}{merMod}(model, metrics = "all", verbose = TRUE, ...)
}
\arguments{
\item{model}{A mixed effects model.}

\item{metrics}{Can be \code{"all"}, \code{"common"} or a character vector of
metrics to be computed (some of \code{c("AIC", "AICc", "BIC", "R2", "ICC", "RMSE", "SIGMA", "LOGLOSS", "SCORE")}). \code{"common"} will compute AIC,
BIC, R2, ICC and RMSE.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame (with one row) and one column per "index" (see
\code{metrics}).
}
\description{
Compute indices of model performance for mixed models.
}
\details{
This method returns the \emph{adjusted ICC} only, as this is typically of
interest when judging the variance attributed to the random effects part of
the model (see also \code{\link[=icc]{icc()}}).
\cr \cr
Furthermore, see 'Details' in \code{\link[=model_performance.lm]{model_performance.lm()}} for
more details on returned indices.
}
\examples{
if (require("lme4")) {
  model <- lmer(Petal.Length ~ Sepal.Length + (1 | Species), data = iris)
  model_performance(model)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_nakagawa.R
\name{r2_nakagawa}
\alias{r2_nakagawa}
\title{Nakagawa's R2 for mixed models}
\usage{
r2_nakagawa(model, by_group = FALSE, tolerance = 1e-05)
}
\arguments{
\item{model}{A mixed effects model.}

\item{by_group}{Logical, if \code{TRUE}, returns the explained variance
at different levels (if there are multiple levels). This is essentially
similar to the variance reduction approach by \cite{Hox (2010), pp. 69-78}.}

\item{tolerance}{Tolerance for singularity check of random effects, to decide
whether to compute random effect variances for the conditional r-squared
or not. Indicates up to which value the convergence result is accepted. When
\code{r2_nakagawa()} returns a warning, stating that random effect variances
can't be computed (and thus, the conditional r-squared is \code{NA}),
decrease the tolerance-level. See also \code{\link[=check_singularity]{check_singularity()}}.}
}
\value{
A list with the conditional and marginal R2 values.
}
\description{
Compute the marginal and conditional r-squared value for
mixed effects models with complex random effects structures.
}
\details{
Marginal and conditional r-squared values for mixed models are calculated
based on \cite{Nakagawa et al. 2017}. For more details on the computation of
the variances, see \code{?insight::get_variance}.
\cr \cr
The marginal r-squared considers only the variance of the fixed effects,
while the conditional r-squared takes both the fixed and random effects into
account. The random effect variances are actually the mean random effect
variances, thus the r-squared value is also appropriate for mixed models
with random slopes or nested random effects (see \cite{Johnson 2014}).
}
\examples{
if (require("lme4")) {
  model <- lmer(Sepal.Length ~ Petal.Length + (1 | Species), data = iris)
  r2_nakagawa(model)
  r2_nakagawa(model, by_group = TRUE)
}
}
\references{
\itemize{
\item Hox, J. J. (2010). Multilevel analysis: techniques and applications
(2nd ed). New York: Routledge.

\item Johnson, P. C. D. (2014). Extension of Nakagawa & Schielzeth’s R2 GLMM
to random slopes models. Methods in Ecology and Evolution, 5(9), 944–946.
\doi{10.1111/2041-210X.12225}

\item Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
obtaining R2 from generalized linear mixed-effects models. Methods in
Ecology and Evolution, 4(2), 133–142. \doi{10.1111/j.2041-210x.2012.00261.x}

\item Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The
coefficient of determination R2 and intra-class correlation coefficient from
generalized linear mixed-effects models revisited and expanded. Journal of
The Royal Society Interface, 14(134), 20170213. \doi{10.1098/rsif.2017.0213}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_xu.R
\name{r2_xu}
\alias{r2_xu}
\title{Xu' R2 (Omega-squared)}
\usage{
r2_xu(model)
}
\arguments{
\item{model}{A linear (mixed) model.}
}
\value{
The R2 value.
}
\description{
Calculates Xu' Omega-squared value, a simple R2 equivalent for
linear mixed models.
}
\details{
\code{r2_xu()} is a crude measure for the explained variance from
linear (mixed) effects models, which is originally denoted as
\ifelse{html}{\out{&Omega;<sup>2</sup>}}{\eqn{\Omega^2}}.
}
\examples{
model <- lm(Sepal.Length ~ Petal.Length + Species, data = iris)
r2_xu(model)
}
\references{
Xu, R. (2003). Measuring explained variation in linear mixed effects models.
Statistics in Medicine, 22(22), 3527–3541. \doi{10.1002/sim.1572}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_rse.R
\name{performance_rse}
\alias{performance_rse}
\title{Residual Standard Error for Linear Models}
\usage{
performance_rse(model)
}
\arguments{
\item{model}{A model.}
}
\value{
Numeric, the residual standard error of \code{model}.
}
\description{
Compute residual standard error of linear models.
}
\details{
The residual standard error is the square root of the residual
sum of squares divided by the residual degrees of freedom.
}
\examples{
data(mtcars)
m <- lm(mpg ~ hp + gear, data = mtcars)
performance_rse(m)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_itemscale.R
\name{check_itemscale}
\alias{check_itemscale}
\title{Describe Properties of Item Scales}
\usage{
check_itemscale(x)
}
\arguments{
\item{x}{An object of class \code{parameters_pca}, as returned by
\code{parameters::principal_components()}.}
}
\value{
A list of data frames, with related measures of internal
consistencies of each subscale.
}
\description{
Compute various measures of internal consistencies
applied to (sub)scales, which items were extracted using
\code{parameters::principal_components()}.
}
\details{
\code{check_itemscale()} calculates various measures of internal
consistencies, such as Cronbach's alpha, item difficulty or discrimination
etc. on subscales which were built from several items. Subscales are
retrieved from the results of \code{parameters::principal_components()}, i.e.
based on how many components were extracted from the PCA,
\code{check_itemscale()} retrieves those variables that belong to a component
and calculates the above mentioned measures.
}
\note{
\itemize{
\item \emph{Item difficulty} should range between 0.2 and 0.8. Ideal value
is \code{p+(1-p)/2} (which mostly is between 0.5 and 0.8). See
\code{\link[=item_difficulty]{item_difficulty()}} for details.

\item For \emph{item discrimination}, acceptable values are 0.20 or higher;
the closer to 1.00 the better. See \code{\link[=item_reliability]{item_reliability()}} for more
details.

\item In case the total \emph{Cronbach's alpha} value is below the
acceptable cut-off of 0.7 (mostly if an index has few items), the
\emph{mean inter-item-correlation} is an alternative measure to indicate
acceptability. Satisfactory range lies between 0.2 and 0.4. See also
\code{\link[=item_intercor]{item_intercor()}}.
}
}
\examples{
# data generation from '?prcomp', slightly modified
C <- chol(S <- toeplitz(.9^(0:15)))
set.seed(17)
X <- matrix(rnorm(1600), 100, 16)
Z <- X \%*\% C
if (require("parameters") && require("psych")) {
  pca <- principal_components(as.data.frame(Z), rotation = "varimax", n = 3)
  pca
  check_itemscale(pca)
}
}
\references{
\itemize{
\item Briggs SR, Cheek JM (1986) The role of factor analysis in the
development and evaluation of personality scales. Journal of Personality,
54(1), 106-148. doi: 10.1111/j.1467-6494.1986.tb00391.x

\item Trochim WMK (2008) Types of Reliability.
(\href{https://conjointly.com/kb/types-of-reliability/}{web})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_distribution.R
\docType{data}
\name{classify_distribution}
\alias{classify_distribution}
\title{Classify the distribution of a model-family using machine learning}
\format{
An object of class \code{randomForest.formula} (inherits from \code{randomForest}) of length 8.
}
\usage{
classify_distribution
}
\description{
Choosing the right distributional family for regression models is essential
to get more accurate estimates and standard errors. This function may help to
Machine learning model trained to classify distributions
}
\details{
Mean accuracy and Kappa of 0.86 and 0.85, repsectively.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_aicc.R
\name{performance_aicc}
\alias{performance_aicc}
\alias{performance_aic}
\title{Compute the AIC or second-order AIC}
\usage{
performance_aicc(x, ...)

performance_aic(x, ...)
}
\arguments{
\item{x}{A model object.}

\item{...}{Currently not used.}
}
\value{
Numeric, the AIC or AICc value.
}
\description{
Compute the AIC or the second-order Akaike's information criterion (AICc).
\code{performance_aic()} is a small wrapper that returns the AIC. It is a generic
function that also works for some models that don't have a AIC method (like
Tweedie models). \code{performance_aicc()} returns the second-order (or "small
sample") AIC that incorporates a correction for small sample sizes.
}
\examples{
m <- lm(mpg ~ wt + cyl + gear + disp, data = mtcars)
AIC(m)
performance_aicc(m)
}
\references{
\itemize{
\item Akaike, H. (1973) Information theory as an extension of the maximum
likelihood principle. In: Second International Symposium on Information
Theory, pp. 267-281. Petrov, B.N., Csaki, F., Eds, Akademiai Kiado, Budapest.
\item Hurvich, C. M., Tsai, C.-L. (1991) Bias of the corrected AIC criterion
for underfitted regression and time series models. Biometrika 78, 499–509.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_tjur.R
\name{r2_tjur}
\alias{r2_tjur}
\title{Tjur's R2 - coefficient of determination (D)}
\usage{
r2_tjur(model)
}
\arguments{
\item{model}{Binomial Model.}
}
\value{
A named vector with the R2 value.
}
\description{
This method calculates the Coefficient of Discrimination \code{D}
(also known as Tjur's R2; \cite{Tjur, 2009}) for generalized linear (mixed) models
for binary outcomes. It is an alternative to other pseudo-R2 values like
Nagelkerke's R2 or Cox-Snell R2. The Coefficient of Discrimination \code{D}
can be read like any other (pseudo-)R2 value.
}
\examples{
model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
r2_tjur(model)
}
\references{
Tjur, T. (2009). Coefficients of determination in logistic regression models - A new proposal: The coefficient of discrimination. The American Statistician, 63(4), 366-372.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_zeroinflation.R
\name{check_zeroinflation}
\alias{check_zeroinflation}
\title{Check for zero-inflation in count models}
\usage{
check_zeroinflation(x, tolerance = 0.05)
}
\arguments{
\item{x}{Fitted model of class \code{merMod}, \code{glmmTMB}, \code{glm},
or \code{glm.nb} (package \pkg{MASS}).}

\item{tolerance}{The tolerance for the ratio of observed and predicted
zeros to considered as over- or underfitting zeros. A ratio
between 1 +/- \code{tolerance} is considered as OK, while a ratio
beyond or below this threshold would indicate over- or underfitting.}
}
\value{
A list with information about the amount of predicted and observed
zeros in the outcome, as well as the ratio between these two values.
}
\description{
\code{check_zeroinflation()} checks whether count models are
over- or underfitting zeros in the outcome.
}
\details{
If the amount of observed zeros is larger than the amount of
predicted zeros, the model is underfitting zeros, which indicates a
zero-inflation in the data. In such cases, it is recommended to use
negative binomial or zero-inflated models.
}
\examples{
if (require("glmmTMB")) {
  data(Salamanders)
  m <- glm(count ~ spp + mined, family = poisson, data = Salamanders)
  check_zeroinflation(m)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_nagelkerke.R
\name{r2_nagelkerke}
\alias{r2_nagelkerke}
\title{Nagelkerke's R2}
\usage{
r2_nagelkerke(model, ...)
}
\arguments{
\item{model}{A generalized linear model, including cumulative links resp.
multinomial models.}

\item{...}{Currently not used.}
}
\value{
A named vector with the R2 value.
}
\description{
Calculate Nagelkerke's pseudo-R2.
}
\examples{
model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
r2_nagelkerke(model)
}
\references{
Nagelkerke, N. J. (1991). A note on a general definition of the coefficient
of determination. Biometrika, 78(3), 691-692.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_mse.R
\name{performance_mse}
\alias{performance_mse}
\alias{mse}
\title{Mean Square Error of Linear Models}
\usage{
performance_mse(model, ...)

mse(model, ...)
}
\arguments{
\item{model}{A model.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
Numeric, the mean square error of \code{model}.
}
\description{
Compute mean square error of linear models.
}
\details{
The mean square error is the mean of the sum of squared residuals, i.e. it
measures the average of the squares of the errors. Less technically speaking,
the mean square error can be considered as the variance of the residuals,
i.e. the variation in the outcome the model doesn't explain. Lower values
(closer to zero) indicate better fit.
}
\examples{
data(mtcars)
m <- lm(mpg ~ hp + gear, data = mtcars)
performance_mse(m)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_distribution.R
\name{check_distribution}
\alias{check_distribution}
\title{Classify the distribution of a model-family using machine learning}
\usage{
check_distribution(model)
}
\arguments{
\item{model}{Typically, a model (that should response to \code{residuals()}).
May also be a numeric vector.}
}
\description{
Choosing the right distributional family for regression models is essential
to get more accurate estimates and standard errors. This function may help to
check a models' distributional family and see if the model-family probably
should be reconsidered. Since it is difficult to exactly predict the correct
model family, consider this function as somewhat experimental.
}
\details{
This function uses an internal random forest model to classify the
distribution from a model-family. Currently, following distributions are
trained (i.e. results of \code{check_distribution()} may be one of the
following): \code{"bernoulli"}, \code{"beta"}, \code{"beta-binomial"},
\code{"binomial"}, \code{"chi"}, \code{"exponential"}, \code{"F"},
\code{"gamma"}, \code{"lognormal"}, \code{"normal"}, \code{"negative binomial"}, \code{"negative binomial (zero-inflated)"}, \code{"pareto"},
\code{"poisson"}, \code{"poisson (zero-inflated)"}, \code{"uniform"} and
\code{"weibull"}.
\cr \cr
Note the similarity between certain distributions according to shape, skewness,
etc. Thus, the predicted distribution may not be perfectly representing the
distributional family of the underlying fitted model, or the response value.
\cr \cr
There is a \code{plot()} method, which shows the probabilities of all predicted
distributions, however, only if the probability is greater than zero.
}
\note{
This function is somewhat experimental and might be improved in future
releases. The final decision on the model-family should also be based on
theoretical aspects and other information about the data and the model.
\cr \cr
There is also a
\href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method}
implemented in the
\href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
if (require("lme4") && require("parameters") && require("see") && require("patchwork")) {
  data(sleepstudy)

  model <<- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  check_distribution(model)
  plot(check_distribution(model))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_convergence.R
\name{check_convergence}
\alias{check_convergence}
\title{Convergence test for mixed effects models}
\usage{
check_convergence(x, tolerance = 0.001, ...)
}
\arguments{
\item{x}{A \code{merMod}-object.}

\item{tolerance}{Indicates up to which value the convergence result is
accepted. The smaller \code{tolerance} is, the stricter the test
will be.}

\item{...}{Currently not used.}
}
\value{
\code{TRUE} if convergence is fine and \code{FALSE} if convergence
is suspicious. Additionally, the convergence value is returned as attribute.
}
\description{
\code{check_convergence()} provides an alternative convergence
test for \code{merMod}-objects.
}
\details{
\subsection{Convergence and log-likelihood}{
Convergence problems typically arise when the model hasn't converged
to a solution where the log-likelihood has a true maximum. This may result
in unreliable and overly complex (or non-estimable) estimates and standard
errors.
}
\subsection{Inspect model convergence}{
\strong{lme4} performs a convergence-check (see \code{?lme4::convergence}),
however, as as discussed \href{https://github.com/lme4/lme4/issues/120}{here}
and suggested by one of the lme4-authors in
\href{https://github.com/lme4/lme4/issues/120#issuecomment-39920269}{this comment},
this check can be too strict. \code{check_convergence()} thus provides an
alternative convergence test for \code{merMod}-objects.
}
\subsection{Resolving convergence issues}{
Convergence issues are not easy to diagnose. The help page on
\code{?lme4::convergence} provides most of the current advice about
how to resolve convergence issues. Another clue might be large parameter
values, e.g. estimates (on the scale of the linear predictor) larger than
10 in (non-identity link) generalized linear model \emph{might} indicate
\href{https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faqwhat-is-complete-or-quasi-complete-separation-in-logisticprobit-regression-and-how-do-we-deal-with-them/}{complete separation}.
Complete separation can be addressed by regularization, e.g. penalized
regression or Bayesian regression with appropriate priors on the fixed effects.
}
\subsection{Convergence versus Singularity}{
Note the different meaning between singularity and convergence: singularity
indicates an issue with the "true" best estimate, i.e. whether the maximum
likelihood estimation for the variance-covariance matrix of the random effects
is positive definite or only semi-definite. Convergence is a question of
whether we can assume that the numerical optimization has worked correctly
or not.
}
}
\examples{
if (require("lme4")) {
  data(cbpp)
  set.seed(1)
  cbpp$x <- rnorm(nrow(cbpp))
  cbpp$x2 <- runif(nrow(cbpp))

  model <- glmer(
    cbind(incidence, size - incidence) ~ period + x + x2 + (1 + x | herd),
    data = cbpp,
    family = binomial()
  )

  check_convergence(model)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/item_split_half.R
\name{item_split_half}
\alias{item_split_half}
\title{Split-Half Reliability}
\usage{
item_split_half(x, digits = 3)
}
\arguments{
\item{x}{A matrix or a data frame.}

\item{digits}{Amount of digits for returned values.}
}
\value{
A list with two elements: the split-half reliability \code{splithalf}
and the Spearman-Brown corrected split-half reliability
\code{spearmanbrown}.
}
\description{
Compute various measures of internal consistencies
for tests or item-scales of questionnaires.
}
\details{
This function calculates the split-half reliability for items in
\code{x}, including the Spearman-Brown adjustment. Splitting is done by
selecting odd versus even columns in \code{x}. A value closer to 1
indicates greater internal consistency.
}
\examples{
data(mtcars)
x <- mtcars[, c("cyl", "gear", "carb", "hp")]
item_split_half(x)
}
\references{
\itemize{
\item Spearman C. 1910. Correlation calculated from faulty data. British
Journal of Psychology (3): 271-295. \doi{10.1111/j.2044-8295.1910.tb00206.x}
}

-Brown W. 1910. Some experimental results in the correlation of mental
abilities. British Journal of Psychology (3): 296-322.
\doi{10.1111/j.2044-8295.1910.tb00207.x}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/item_intercor.R
\name{item_intercor}
\alias{item_intercor}
\title{Mean Inter-Item-Correlation}
\usage{
item_intercor(x, method = c("pearson", "spearman", "kendall"))
}
\arguments{
\item{x}{A matrix as returned by the \code{cor()}-function,
or a data frame with items (e.g. from a test or questionnaire).}

\item{method}{Correlation computation method. May be one of
\code{"spearman"} (default), \code{"pearson"} or \code{"kendall"}.
You may use initial letter only.}
}
\value{
The mean inter-item-correlation value for \code{x}.
}
\description{
Compute various measures of internal consistencies
for tests or item-scales of questionnaires.
}
\details{
This function calculates a mean inter-item-correlation, i.e.
a correlation matrix of \code{x} will be computed (unless
\code{x} is already a matrix as returned by the \code{cor()}-function)
and the mean of the sum of all item's correlation values is returned.
Requires either a data frame or a computed \code{cor()}-object.
\cr \cr
\dQuote{Ideally, the average inter-item correlation for a set of
items should be between .20 and .40, suggesting that while the
items are reasonably homogenous, they do contain sufficiently
unique variance so as to not be isomorphic with each other.
When values are lower than .20, then the items may not be
representative of the same content domain. If values are higher than
.40, the items may be only capturing a small bandwidth of the construct.}
\cite{(Piedmont 2014)}
}
\examples{
data(mtcars)
x <- mtcars[, c("cyl", "gear", "carb", "hp")]
item_intercor(x)
}
\references{
Piedmont RL. 2014. Inter-item Correlations. In: Michalos AC (eds)
Encyclopedia of Quality of Life and Well-Being Research. Dordrecht: Springer,
3303-3304. \doi{10.1007/978-94-007-0753-5_1493}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_bayes.R
\name{r2_bayes}
\alias{r2_bayes}
\alias{r2_posterior}
\alias{r2_posterior.brmsfit}
\alias{r2_posterior.stanreg}
\alias{r2_posterior.BFBayesFactor}
\title{Bayesian R2}
\usage{
r2_bayes(model, robust = TRUE, ci = 0.95, verbose = TRUE, ...)

r2_posterior(model, ...)

\method{r2_posterior}{brmsfit}(model, verbose = TRUE, ...)

\method{r2_posterior}{stanreg}(model, verbose = TRUE, ...)

\method{r2_posterior}{BFBayesFactor}(model, average = FALSE, prior_odds = NULL, verbose = TRUE, ...)
}
\arguments{
\item{model}{A Bayesian regression model (from \strong{brms},
\strong{rstanarm}, \strong{BayesFactor}, etc).}

\item{robust}{Logical, if \code{TRUE}, the median instead of mean is used to
calculate the central tendency of the variances.}

\item{ci}{Value or vector of probability of the CI (between 0 and 1) to be
estimated.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to \code{r2_posterior()}.}

\item{average}{Compute model-averaged index? See \code{\link[bayestestR:weighted_posteriors]{bayestestR::weighted_posteriors()}}.}

\item{prior_odds}{Optional vector of prior odds for the models compared to
the first model (or the denominator, for \code{BFBayesFactor} objects). For
\code{data.frame}s, this will be used as the basis of weighting.}
}
\value{
A list with the Bayesian R2 value. For mixed models, a list with the
Bayesian R2 value and the marginal Bayesian R2 value. The standard errors
and credible intervals for the R2 values are saved as attributes.
}
\description{
Compute R2 for Bayesian models. For mixed models (including a
random part), it additionally computes the R2 related to the fixed effects
only (marginal R2). While \code{r2_bayes()} returns a single R2 value,
\code{r2_posterior()} returns a posterior sample of Bayesian R2 values.
}
\details{
\code{r2_bayes()} returns an "unadjusted" R2 value. See
\code{\link[=r2_loo]{r2_loo()}} to calculate a LOO-adjusted R2, which comes
conceptually closer to an adjusted R2 measure.
\cr \cr
For mixed models, the conditional and marginal R2 are returned. The marginal
R2 considers only the variance of the fixed effects, while the conditional
R2 takes both the fixed and random effects into account.
\cr \cr
\code{r2_posterior()} is the actual workhorse for \code{r2_bayes()} and
returns a posterior sample of Bayesian R2 values.
}
\examples{
library(performance)
if (require("rstanarm") && require("rstantools")) {
  model <- stan_glm(mpg ~ wt + cyl, data = mtcars, chains = 1, iter = 500, refresh = 0)
  r2_bayes(model)

  model <- stan_lmer(
    Petal.Length ~ Petal.Width + (1 | Species),
    data = iris,
    chains = 1,
    iter = 500,
    refresh = 0
  )
  r2_bayes(model)
}
\dontrun{
if (require("BayesFactor")) {
  data(mtcars)

  BFM <- generalTestBF(mpg ~ qsec + gear, data = mtcars, progress = FALSE)
  FM <- lm(mpg ~ qsec + gear, data = mtcars)

  r2_bayes(FM)
  r2_bayes(BFM[3])
  r2_bayes(BFM, average = TRUE) # across all models


  # with random effects:
  mtcars$gear <- factor(mtcars$gear)
  model <- lmBF(
    mpg ~ hp + cyl + gear + gear:wt,
    mtcars,
    progress = FALSE,
    whichRandom = c("gear", "gear:wt")
  )
  r2_bayes(model)
}

if (require("brms")) {
  model <- brms::brm(mpg ~ wt + cyl, data = mtcars)
  r2_bayes(model)

  model <- brms::brm(Petal.Length ~ Petal.Width + (1 | Species), data = iris)
  r2_bayes(model)
}
}
}
\references{
Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2018). R-squared for Bayesian regression models. The American Statistician, 1–6. \doi{10.1080/00031305.2018.1549100}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.bayesian.R
\name{model_performance.stanreg}
\alias{model_performance.stanreg}
\alias{model_performance.BFBayesFactor}
\title{Performance of Bayesian Models}
\usage{
\method{model_performance}{stanreg}(model, metrics = "all", verbose = TRUE, ...)

\method{model_performance}{BFBayesFactor}(
  model,
  metrics = "all",
  verbose = TRUE,
  average = FALSE,
  prior_odds = NULL,
  ...
)
}
\arguments{
\item{model}{Object of class \code{stanreg} or \code{brmsfit}.}

\item{metrics}{Can be \code{"all"}, \code{"common"} or a character vector of
metrics to be computed (some of \code{c("LOOIC", "WAIC", "R2", "R2_adj", "RMSE", "SIGMA", "LOGLOSS", "SCORE")}). \code{"common"} will compute LOOIC,
WAIC, R2 and RMSE.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}

\item{average}{Compute model-averaged index? See \code{\link[bayestestR:weighted_posteriors]{bayestestR::weighted_posteriors()}}.}

\item{prior_odds}{Optional vector of prior odds for the models compared to
the first model (or the denominator, for \code{BFBayesFactor} objects). For
\code{data.frame}s, this will be used as the basis of weighting.}
}
\value{
A data frame (with one row) and one column per "index" (see
\code{metrics}).
}
\description{
Compute indices of model performance for (general) linear models.
}
\details{
Depending on \code{model}, the following indices are computed:
\itemize{
\item{\strong{ELPD}} {expected log predictive density. Larger ELPD values
mean better fit. See \code{\link[=looic]{looic()}}.}

\item{\strong{LOOIC}} {leave-one-out cross-validation (LOO) information
criterion. Lower LOOIC values mean better fit. See \code{\link[=looic]{looic()}}.}

\item{\strong{WAIC}} {widely applicable information criterion. Lower WAIC
values mean better fit. See \code{?loo::waic}.}

\item{\strong{R2}} {r-squared value, see \code{\link[=r2_bayes]{r2_bayes()}}.}

\item{\strong{R2_adjusted}} {LOO-adjusted r-squared, see
\code{\link[=r2_loo]{r2_loo()}}.}

\item{\strong{RMSE}} {root mean squared error, see
\code{\link[=performance_rmse]{performance_rmse()}}.}

\item{\strong{SIGMA}} {residual standard deviation, see
\code{\link[insight:get_sigma]{insight::get_sigma()}}.}

\item{\strong{LOGLOSS}} {Log-loss, see \code{\link[=performance_logloss]{performance_logloss()}}.}

\item{\strong{SCORE_LOG}} {score of logarithmic proper scoring rule, see
\code{\link[=performance_score]{performance_score()}}.}

\item{\strong{SCORE_SPHERICAL}} {score of spherical proper scoring rule,
see \code{\link[=performance_score]{performance_score()}}.}

\item{\strong{PCP}} {percentage of correct predictions, see
\code{\link[=performance_pcp]{performance_pcp()}}.}
}
}
\examples{
\dontrun{
if (require("rstanarm") && require("rstantools")) {
  model <- stan_glm(mpg ~ wt + cyl, data = mtcars, chains = 1, iter = 500, refresh = 0)
  model_performance(model)

  model <- stan_glmer(
    mpg ~ wt + cyl + (1 | gear),
    data = mtcars,
    chains = 1,
    iter = 500,
    refresh = 0
  )
  model_performance(model)
}

if (require("BayesFactor") && require("rstantools")) {
  model <- generalTestBF(carb ~ am + mpg, mtcars)

  model_performance(model)
  model_performance(model[3])

  model_performance(model, average = TRUE)
}
}
}
\references{
Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2018).
R-squared for Bayesian regression models. The American Statistician, The
American Statistician, 1-6.
}
\seealso{
\link{r2_bayes}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_outliers.R
\name{check_outliers}
\alias{check_outliers}
\alias{check_outliers.default}
\alias{check_outliers.numeric}
\alias{check_outliers.data.frame}
\title{Outliers detection (check for influential observations)}
\usage{
check_outliers(x, ...)

\method{check_outliers}{default}(x, method = c("cook", "pareto"), threshold = NULL, ...)

\method{check_outliers}{numeric}(x, method = "zscore_robust", threshold = NULL, ...)

\method{check_outliers}{data.frame}(x, method = "mahalanobis", threshold = NULL, ...)
}
\arguments{
\item{x}{A model or a data.frame object.}

\item{...}{When \code{method = "ics"}, further arguments in \code{...} are
passed down to \code{ICSOutlier::ics.outlier()}.}

\item{method}{The outlier detection method(s). Can be "all" or some of
c("cook", "pareto", "zscore", "zscore_robust", "iqr", "eti", "hdi", "bci",
"mahalanobis", "mahalanobis_robust", "mcd", "ics", "optics", "lof").}

\item{threshold}{A list containing the threshold values for each method (e.g.
\code{list('mahalanobis' = 7, 'cook' = 1)}), above which an observation is
considered as outlier. If \code{NULL}, default values will be used (see
'Details'). If a numeric value is given, it will be used as the threshold
for any of the method run.}
}
\value{
A logical vector of the detected outliers with a nice printing
method: a check (message) on whether outliers were detected or not. The
information on the distance measure and whether or not an observation is
considered as outlier can be recovered with the \code{as.data.frame}
function.
}
\description{
Checks for and locates influential observations (i.e.,
"outliers") via several distance and/or clustering methods. If several
methods are selected, the returned "Outlier" vector will be a composite
outlier score, made of the average of the binary (0 or 1) results of each
method. It represents the probability of each observation of being
classified as an outlier by at least one method. The decision rule used by
default is to classify as outliers observations which composite outlier
score is superior or equal to 0.5 (i.e., that were classified as outliers
by at least half of the methods). See the \strong{Details} section below
for a description of the methods.
}
\details{
Outliers can be defined as particularly influential observations.
Most methods rely on the computation of some distance metric, and the
observations greater than a certain threshold are considered outliers.
Importantly, outliers detection methods are meant to provide information to
consider for the researcher, rather than to be an automatized procedure
which mindless application is a substitute for thinking.

An \strong{example sentence} for reporting the usage of the composite method
could be:

\emph{"Based on a composite outlier score (see the 'check_outliers' function
in the 'performance' R package; Lüdecke et al., 2021) obtained via the joint
application of multiple outliers detection algorithms (Z-scores, Iglewicz,
1993; Interquartile range (IQR); Mahalanobis distance, Cabana, 2019; Robust
Mahalanobis distance, Gnanadesikan & Kettenring, 1972; Minimum Covariance
Determinant, Leys et al., 2018; Invariant Coordinate Selection, Archimbaud et
al., 2018; OPTICS, Ankerst et al., 1999; Isolation Forest, Liu et al. 2008;
and Local Outlier Factor, Breunig et al., 2000), we excluded n participants
that were classified as outliers by at least half of the methods used."}

\subsection{Model-specific methods}{
\itemize{
\item \strong{Cook's Distance}:
Among outlier detection methods, Cook's distance and leverage are less
common than the basic Mahalanobis distance, but still used. Cook's distance
estimates the variations in regression coefficients after removing each
observation, one by one (Cook, 1977). Since Cook's distance is in the metric
of an F distribution with p and n-p degrees of freedom, the median point of
the quantile distribution can be used as a cut-off (Bollen, 1985). A common
approximation or heuristic is to use 4 divided by the numbers of
observations, which usually corresponds to a lower threshold (i.e., more
outliers are detected). This only works for Frequentist models. For Bayesian
models, see \code{pareto}.

\item \strong{Pareto}:
The reliability and approximate convergence of Bayesian models can be
assessed using the estimates for the shape parameter k of the generalized
Pareto distribution. If the estimated tail shape parameter k exceeds 0.5, the
user should be warned, although in practice the authors of the \code{loo}
package observed good performance for values of k up to 0.7 (the default
threshold used by \code{performance}).
}}

\subsection{Univariate methods}{
\itemize{
\item \strong{Z-scores} \verb{("zscore", "zscore_robust")}:
The Z-score, or standard score, is a way of describing a data point as
deviance from a central value, in terms of standard deviations from the mean
(\code{"zscore"}) or, as it is here the case (\code{"zscore_robust"}) by
default (Iglewicz, 1993), in terms of Median Absolute Deviation (MAD) from
the median (which are robust measures of dispersion and centrality). The
default threshold to classify outliers is 1.959 (\code{threshold = list("zscore" = 1.959)}), corresponding to the 2.5\\% (\code{qnorm(0.975)})
most extreme observations (assuming the data is normally distributed).
Importantly, the Z-score method is univariate: it is computed column by
column. If a dataframe is passed, the Z-score is calculated for each
variable separately, and the maximum (absolute) Z-score is kept for each
observations. Thus, all observations that are extreme on at least one
variable might be detected as outliers. Thus, this method is not suited for
high dimensional data (with many columns), returning too liberal results
(detecting many outliers).

\item \strong{IQR} \code{("iqr")}:
Using the IQR (interquartile range) is a robust method developed by John
Tukey, which often appears in box-and-whisker plots (e.g., in
\code{geom_boxplot}). The interquartile range is the range between the first
and the third quartiles. Tukey considered as outliers any data point that
fell outside of either 1.5 times (the default threshold) the IQR below the
first or above the third quartile. Similar to the Z-score method, this is a
univariate method for outliers detection, returning outliers detected for at
least one column, and might thus not be suited to high dimensional data.

\item \strong{CI} \verb{("ci", "eti", "hdi", "bci")}:
Another univariate method is to compute, for each variable, some sort of
"confidence" interval and consider as outliers values lying beyond the edges
of that interval. By default, \code{"ci"} computes the Equal-Tailed Interval
(\code{"eti"}), but other types of intervals are available, such as Highest
Density Interval (\code{"hdi"}) or the Bias Corrected and Accelerated
Interval (\code{"bci"}). The default threshold is \code{0.95}, considering
as outliers all observations that are outside the 95\\% CI on any of the
variable. See \code{\link[bayestestR:ci]{bayestestR::ci()}} for more details
about the intervals.
}}

\subsection{Multivariate methods}{
\itemize{
\item \strong{Mahalanobis Distance}:
Mahalanobis distance (Mahalanobis, 1930) is often used for multivariate
outliers detection as this distance takes into account the shape of the
observations. The default \code{threshold} is often arbitrarily set to some
deviation (in terms of SD or MAD) from the mean (or median) of the
Mahalanobis distance. However, as the Mahalanobis distance can be
approximated by a Chi squared distribution (Rousseeuw & Van Zomeren, 1990),
we can use the alpha quantile of the chi-square distribution with k degrees
of freedom (k being the number of columns). By default, the alpha threshold
is set to 0.025 (corresponding to the 2.5\\% most extreme observations;
Cabana, 2019). This criterion is a natural extension of the median plus or
minus a coefficient times the MAD method (Leys et al., 2013).

\item \strong{Robust Mahalanobis Distance}:
A robust version of Mahalanobis distance using an Orthogonalized
Gnanadesikan-Kettenring pairwise estimator (Gnanadesikan & Kettenring,
1972). Requires the \pkg{bigutilsr} package. See the
\code{bigutilsr::dist_ogk()} function.

\item \strong{Minimum Covariance Determinant (MCD)}:
Another robust version of Mahalanobis. Leys et al. (2018) argue that
Mahalanobis Distance is not a robust way to determine outliers, as it uses
the means and covariances of all the data - including the outliers - to
determine individual difference scores. Minimum Covariance Determinant
calculates the mean and covariance matrix based on the most central subset of
the data (by default, 66\\%), before computing the Mahalanobis Distance. This
is deemed to be a more robust method of identifying and removing outliers
than regular Mahalanobis distance.

\item \strong{Invariant Coordinate Selection (ICS)}:
The outlier are detected using ICS, which by default uses an alpha threshold
of 0.025 (corresponding to the 2.5\\% most extreme observations) as a cut-off
value for outliers classification. Refer to the help-file of
\code{ICSOutlier::ics.outlier()} to get more details about this procedure.
Note that \code{method = "ics"} requires both \pkg{ICS} and \pkg{ICSOutlier}
to be installed, and that it takes some time to compute the results.

\item \strong{OPTICS}:
The Ordering Points To Identify the Clustering Structure (OPTICS) algorithm
(Ankerst et al., 1999) is using similar concepts to DBSCAN (an unsupervised
clustering technique that can be used for outliers detection). The threshold
argument is passed as \code{minPts}, which corresponds to the minimum size
of a cluster. By default, this size is set at 2 times the number of columns
(Sander et al., 1998). Compared to the others techniques, that will always
detect several outliers (as these are usually defined as a percentage of
extreme values), this algorithm functions in a different manner and won't
always detect outliers. Note that \code{method = "optics"} requires the
\pkg{dbscan} package to be installed, and that it takes some time to compute
the results.

\item \strong{Isolation Forest}:
The outliers are detected using the anomaly score of an isolation forest (a
class of random forest). The default threshold of 0.025 will classify as
outliers the observations located at \verb{qnorm(1-0.025) * MAD)} (a robust
equivalent of SD) of the median (roughly corresponding to the 2.5\\% most
extreme observations). Requires the \pkg{solitude} package.

\item \strong{Local Outlier Factor}:
Based on a K nearest neighbours algorithm, LOF compares the local density of
an point to the local densities of its neighbors instead of computing a
distance from the center (Breunig et al., 2000). Points that have a
substantially lower density than their neighbors are considered outliers. A
LOF score of approximately 1 indicates that density around the point is
comparable to its neighbors. Scores significantly larger than 1 indicate
outliers. The default threshold of 0.025 will classify as outliers the
observations located at \verb{qnorm(1-0.025) * SD)} of the log-transformed
LOF distance. Requires the \pkg{dbscan} package.
}}

\subsection{Threshold specification}{
Default thresholds are currently specified as follows:

\preformatted{
list(
  zscore = stats::qnorm(p = 1 - 0.025),
  iqr = 1.5,
  ci = 0.95,
  cook = stats::qf(0.5, ncol(x), nrow(x) - ncol(x)),
  pareto = 0.7,
  mahalanobis = stats::qchisq(p = 1 - 0.025, df = ncol(x)),
  robust = stats::qchisq(p = 1 - 0.025, df = ncol(x)),
  mcd = stats::qchisq(p = 1 - 0.025, df = ncol(x)),
  ics = 0.025,
  optics = 2 * ncol(x),
  iforest = 0.025,
  lof = 0.025
)
}}
}
\note{
There is also a
\href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method}
implemented in the
\href{https://easystats.github.io/see/}{\pkg{see}-package}. \strong{Please
note} that the range of the distance-values along the y-axis is re-scaled
to range from 0 to 1.
}
\examples{
data <- mtcars # Size nrow(data) = 32

# For single variables ------------------------------------------------------
outliers_list <- check_outliers(data$mpg) # Find outliers
outliers_list # Show the row index of the outliers
as.numeric(outliers_list) # The object is a binary vector...
filtered_data <- data[!outliers_list, ] # And can be used to filter a dataframe
nrow(filtered_data) # New size, 28 (4 outliers removed)

# Find all observations beyond +/- 2 SD
check_outliers(data$mpg, method = "zscore", threshold = 2)

# For dataframes ------------------------------------------------------
check_outliers(data) # It works the same way on dataframes

# You can also use multiple methods at once
outliers_list <- check_outliers(data, method = c(
  "mahalanobis",
  "iqr",
  "zscore"
))
outliers_list

# Using `as.data.frame()`, we can access more details!
outliers_info <- as.data.frame(outliers_list)
head(outliers_info)
outliers_info$Outlier # Including the probability of being an outlier

# And we can be more stringent in our outliers removal process
filtered_data <- data[outliers_info$Outlier < 0.1, ]

# We can run the function stratified by groups using `{dplyr}` package:
if (require("poorman")) {
  iris \%>\%
    group_by(Species) \%>\%
    check_outliers()
}
\dontrun{
# You can also run all the methods
check_outliers(data, method = "all")

# For statistical models ---------------------------------------------
# select only mpg and disp (continuous)
mt1 <- mtcars[, c(1, 3, 4)]
# create some fake outliers and attach outliers to main df
mt2 <- rbind(mt1, data.frame(
  mpg = c(37, 40), disp = c(300, 400),
  hp = c(110, 120)
))
# fit model with outliers
model <- lm(disp ~ mpg + hp, data = mt2)

outliers_list <- check_outliers(model)

if (require("see")) {
  plot(outliers_list)
}

insight::get_data(model)[outliers_list, ] # Show outliers data

if (require("MASS")) {
  check_outliers(model, method = c("mahalabonis", "mcd"))
}
if (require("ICS")) {
  # This one takes some seconds to finish...
  check_outliers(model, method = "ics")
}
}
}
\references{
\itemize{
\item Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2018). ICS for
multivariate outlier detection with application to quality control.
Computational Statistics & Data Analysis, 128, 184-199.
\doi{10.1016/j.csda.2018.06.011}
\item Gnanadesikan, R., & Kettenring, J. R. (1972). Robust estimates, residuals,
and outlier detection with multiresponse data. Biometrics, 81-124.
\item Bollen, K. A., & Jackman, R. W. (1985). Regression diagnostics: An
expository treatment of outliers and influential cases. Sociological Methods
& Research, 13(4), 510-542.
\item Cabana, E., Lillo, R. E., & Laniado, H. (2019). Multivariate outlier
detection based on a robust Mahalanobis distance with shrinkage estimators.
arXiv preprint arXiv:1904.02596.
\item Cook, R. D. (1977). Detection of influential observation in linear
regression. Technometrics, 19(1), 15-18.
\item Iglewicz, B., & Hoaglin, D. C. (1993). How to detect and handle outliers
(Vol. 16). Asq Press.
\item Leys, C., Klein, O., Dominicy, Y., & Ley, C. (2018). Detecting
multivariate outliers: Use a robust variant of Mahalanobis distance. Journal
of Experimental Social Psychology, 74, 150-156.
\item Liu, F. T., Ting, K. M., & Zhou, Z. H. (2008, December). Isolation forest.
In 2008 Eighth IEEE International Conference on Data Mining (pp. 413-422).
IEEE.
\item Lüdecke, D., Ben-Shachar, M. S., Patil, I., Waggoner, P., & Makowski, D.
(2021). performance: An R package for assessment, comparison and testing of
statistical models. Journal of Open Source Software, 6(60), 3139.
\doi{10.21105/joss.03139}
\item Rousseeuw, P. J., & Van Zomeren, B. C. (1990). Unmasking multivariate
outliers and leverage points. Journal of the American Statistical
association, 85(411), 633-639.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_overdispersion.R
\name{check_overdispersion}
\alias{check_overdispersion}
\title{Check overdispersion of GL(M)M's}
\usage{
check_overdispersion(x, ...)
}
\arguments{
\item{x}{Fitted model of class \code{merMod}, \code{glmmTMB}, \code{glm},
or \code{glm.nb} (package \pkg{MASS}).}

\item{...}{Currently not used.}
}
\value{
A list with results from the overdispersion test, like chi-squared
statistics, p-value or dispersion ratio.
}
\description{
\code{check_overdispersion()} checks generalized linear (mixed)
models for overdispersion.
}
\details{
Overdispersion occurs when the observed variance is higher than the
variance of a theoretical model. For Poisson models, variance increases
with the mean and, therefore, variance usually (roughly) equals the mean
value. If the variance is much higher, the data are "overdispersed".

\subsection{Interpretation of the Dispersion Ratio}{
If the dispersion ratio is close to one, a Poisson model fits well to the
data. Dispersion ratios larger than one indicate overdispersion, thus a
negative binomial model or similar might fit better to the data. A p-value <
.05 indicates overdispersion.
}

\subsection{Overdispersion in Poisson Models}{
For Poisson models, the overdispersion test is based on the code from
\cite{Gelman and Hill (2007), page 115}.
}

\subsection{Overdispersion in Mixed Models}{
For \code{merMod}- and \code{glmmTMB}-objects, \code{check_overdispersion()}
is based on the code in the
\href{http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html}{GLMM FAQ},
section \emph{How can I deal with overdispersion in GLMMs?}. Note that this
function only returns an \emph{approximate} estimate of an overdispersion
parameter, and is probably inaccurate for zero-inflated mixed models (fitted
with \code{glmmTMB}).
}

\subsection{How to fix Overdispersion}{
Overdispersion can be fixed by either modeling the dispersion parameter, or
by choosing a different distributional family (like Quasi-Poisson, or
negative binomial, see \cite{Gelman and Hill (2007), pages 115-116}).
}
}
\examples{
if (require("glmmTMB")) {
  data(Salamanders)
  m <- glm(count ~ spp + mined, family = poisson, data = Salamanders)
  check_overdispersion(m)

  m <- glmmTMB(
    count ~ mined + spp + (1 | site),
    family = poisson,
    data = Salamanders
  )
  check_overdispersion(m)
}
}
\references{
\itemize{
\item Bolker B et al. (2017):
\href{http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html}{GLMM FAQ.}
\item Gelman, A., & Hill, J. (2007). Data analysis using regression and
multilevel/hierarchical models. Cambridge; New York: Cambridge University
Press.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_multimodal.R
\name{check_multimodal}
\alias{check_multimodal}
\title{Check if a distribution is unimodal or multimodal}
\usage{
check_multimodal(x, ...)
}
\arguments{
\item{x}{A numeric vector or a data frame.}

\item{...}{Arguments passed to or from other methods.}
}
\description{
For univariate distributions (one-dimensional vectors), this functions
performs a Ameijeiras-Alonso et al. (2018) excess mass test. For multivariate
distributions (dataframes), it uses mixture modelling. However, it seems that
it always returns a significant result (suggesting that the distribution is
multimodal). A better method might be needed here.
}
\examples{
\dontrun{
if (require("multimode")) {
  # Univariate
  x <- rnorm(1000)
  check_multimodal(x)
}

if (require("multimode") && require("mclust")) {
  x <- c(rnorm(1000), rnorm(1000, 2))
  check_multimodal(x)

  # Multivariate
  m <- data.frame(
    x = rnorm(200),
    y = rbeta(200, 2, 1)
  )
  plot(m$x, m$y)
  check_multimodal(m)

  m <- data.frame(
    x = c(rnorm(100), rnorm(100, 4)),
    y = c(rbeta(100, 2, 1), rbeta(100, 1, 4))
  )
  plot(m$x, m$y)
  check_multimodal(m)
}
}
}
\references{
\itemize{
\item Ameijeiras-Alonso, J., Crujeiras, R. M., \& Rodríguez-Casal, A. (2019).
Mode testing, critical bandwidth and excess mass. Test, 28(3), 900-919.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_logloss.R
\name{performance_logloss}
\alias{performance_logloss}
\title{Log Loss}
\usage{
performance_logloss(model, verbose = TRUE, ...)
}
\arguments{
\item{model}{Model with binary outcome.}

\item{verbose}{Toggle off warnings.}

\item{...}{Currently not used.}
}
\value{
Numeric, the log loss of \code{model}.
}
\description{
Compute the log loss for models with binary outcome.
}
\details{
Logistic regression models predict the probability of an outcome of being a
"success" or "failure" (or 1 and 0 etc.). \code{performance_logloss()} evaluates
how good or bad the predicted probabilities are. High values indicate bad
predictions, while low values indicate good predictions. The lower the
log-loss, the better the model predicts the outcome.
}
\examples{
data(mtcars)
m <- glm(formula = vs ~ hp + wt, family = binomial, data = mtcars)
performance_logloss(m)
}
\seealso{
\code{\link[=performance_score]{performance_score()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_zeroinflated.R
\name{r2_zeroinflated}
\alias{r2_zeroinflated}
\title{R2 for models with zero-inflation}
\usage{
r2_zeroinflated(model, method = c("default", "correlation"))
}
\arguments{
\item{model}{A model.}

\item{method}{Indicates the method to calculate R2. See 'Details'. May be
abbreviated.}
}
\value{
For the default-method, a list with the R2 and adjusted R2 values.
For \code{method = "correlation"}, a named numeric vector with the
correlation-based R2 value.
}
\description{
Calculates R2 for models with zero-inflation component, including mixed
effects models.
}
\details{
The default-method calculates an R2 value based on the residual variance
divided by the total variance. For \code{method = "correlation"}, R2 is a
correlation-based measure, which is rather crude. It simply computes the
squared correlation between the model's actual and predicted response.
}
\examples{
\donttest{
if (require("pscl")) {
  data(bioChemists)
  model <- zeroinfl(
    art ~ fem + mar + kid5 + ment | kid5 + phd,
    data = bioChemists
  )

  r2_zeroinflated(model)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.rma.R
\name{model_performance.rma}
\alias{model_performance.rma}
\title{Performance of Meta-Analysis Models}
\usage{
\method{model_performance}{rma}(model, metrics = "all", verbose = TRUE, ...)
}
\arguments{
\item{model}{A \code{rma} object as returned by \code{metafor::rma()}.}

\item{metrics}{Can be \code{"all"} or a character vector of metrics to be
computed (some of \code{c("AIC", "BIC", "I2", "H2", "TAU2", "R2", "CochransQ", "QE", "Omnibus", "QM")}).}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame (with one row) and one column per "index" (see
\code{metrics}).
}
\description{
Compute indices of model performance for meta-analysis model from the
\pkg{metafor} package.
}
\details{
\subsection{Indices of fit}{
\itemize{
\item{\strong{AIC}} {Akaike's Information Criterion, see
\code{?stats::AIC}}

\item{\strong{BIC}} {Bayesian Information Criterion, see
\code{?stats::BIC}}

\item \strong{I2}: For a random effects model, \code{I2} estimates (in
percent) how much of the total variability in the effect size estimates
can be attributed to heterogeneity among the true effects. For a
mixed-effects model, \code{I2} estimates how much of the unaccounted
variability can be attributed to residual heterogeneity.

\item \strong{H2}: For a random-effects model, \code{H2} estimates the
ratio of the total amount of variability in the effect size estimates to
the amount of sampling variability. For a mixed-effects model, \code{H2}
estimates the ratio of the unaccounted variability in the effect size
estimates to the amount of sampling variability.

\item \strong{TAU2}: The amount of (residual) heterogeneity in the random
or mixed effects model.

\item \strong{CochransQ (QE)}: Test for (residual) Heterogeneity. Without
moderators in the model, this is simply Cochran's Q-test.

\item \strong{Omnibus (QM)}: Omnibus test of parameters.

\item \strong{R2}: Pseudo-R2-statistic, which indicates the amount of
heterogeneity accounted for by the moderators included in a fixed-effects
model.
}
See the documentation for \code{?metafor::fitstats}.
}
}
\examples{
if (require("metafor")) {
  data(dat.bcg)
  dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)
  model <- rma(yi, vi, data = dat, method = "REML")
  model_performance(model)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.ivreg.R
\name{model_performance.ivreg}
\alias{model_performance.ivreg}
\title{Performance of instrumental variable regression models}
\usage{
\method{model_performance}{ivreg}(model, metrics = "all", verbose = TRUE, ...)
}
\arguments{
\item{model}{A model.}

\item{metrics}{Can be \code{"all"}, \code{"common"} or a character vector of
metrics to be computed (some of \code{c("AIC", "AICc", "BIC", "R2", "RMSE", "SIGMA", "Sargan", "Wu_Hausman")}). \code{"common"} will compute AIC, BIC,
R2 and RMSE.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}
}
\description{
Performance of instrumental variable regression models
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/item_reliability.R
\name{item_reliability}
\alias{item_reliability}
\title{Reliability Test for Items or Scales}
\usage{
item_reliability(x, standardize = FALSE, digits = 3)
}
\arguments{
\item{x}{A matrix or a data frame.}

\item{standardize}{Logical, if \code{TRUE}, the data frame's vectors will be
standardized. Recommended when the variables have different measures /
scales.}

\item{digits}{Amount of digits for returned values.}
}
\value{
A data frame with the corrected item-total correlations (\emph{item
discrimination}, column \code{item_discrimination}) and Cronbach's Alpha
(if item deleted, column \code{alpha_if_deleted}) for each item
of the scale, or \code{NULL} if data frame had too less columns.
}
\description{
Compute various measures of internal consistencies
for tests or item-scales of questionnaires.
}
\details{
This function calculates the item discriminations (corrected item-total
correlations for each item of \code{x} with the remaining items) and the
Cronbach's alpha for each item, if it was deleted from the scale. The
absolute value of the item discrimination indices should be above 0.1. An
index between 0.1 and 0.3 is considered as "fair", while an index above 0.3
(or below -0.3) is "good". Items with low discrimination indices are often
ambiguously worded and should be examined. Items with negative indices should
be examined to determine why a negative value was obtained (e.g. reversed
answer categories regarding positive and negative poles).
}
\examples{
data(mtcars)
x <- mtcars[, c("cyl", "gear", "carb", "hp")]
item_reliability(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R, R/print_md.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{display}
\alias{print_md}
\alias{print_html}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{insight}{\code{\link[insight]{display}}, \code{\link[insight:display]{print_html}}, \code{\link[insight:display]{print_md}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_heterogeneity_bias.R
\name{check_heterogeneity_bias}
\alias{check_heterogeneity_bias}
\title{Check model predictor for heterogeneity bias}
\usage{
check_heterogeneity_bias(x, select = NULL, group = NULL)
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
\code{check_heterogeneity_bias()} checks if model predictors or variables may
cause a heterogeneity bias, i.e. if variables have a within- and/or
between-effect.
}
\examples{
data(iris)
iris$ID <- sample(1:4, nrow(iris), replace = TRUE) # fake-ID
check_heterogeneity_bias(iris, select = c("Sepal.Length", "Petal.Length"), group = "ID")
}
\seealso{
For further details, read the vignette
\url{https://easystats.github.io/datawizard/articles/demean.html} and also
see documentation for \code{?datawizard::demean}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_loo.R
\name{r2_loo}
\alias{r2_loo}
\alias{r2_loo_posterior}
\alias{r2_loo_posterior.brmsfit}
\alias{r2_loo_posterior.stanreg}
\title{LOO-adjusted R2}
\usage{
r2_loo(model, robust = TRUE, ci = 0.95, verbose = TRUE, ...)

r2_loo_posterior(model, ...)

\method{r2_loo_posterior}{brmsfit}(model, verbose = TRUE, ...)

\method{r2_loo_posterior}{stanreg}(model, verbose = TRUE, ...)
}
\arguments{
\item{model}{A Bayesian regression model (from \strong{brms},
\strong{rstanarm}, \strong{BayesFactor}, etc).}

\item{robust}{Logical, if \code{TRUE}, the median instead of mean is used to
calculate the central tendency of the variances.}

\item{ci}{Value or vector of probability of the CI (between 0 and 1) to be
estimated.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to \code{r2_posterior()}.}
}
\value{
A list with the Bayesian R2 value. For mixed models, a list with the
Bayesian R2 value and the marginal Bayesian R2 value. The standard errors
and credible intervals for the R2 values are saved as attributes.

A list with the LOO-adjusted R2 value. The standard errors
and credible intervals for the R2 values are saved as attributes.
}
\description{
Compute LOO-adjusted R2.
}
\details{
\code{r2_loo()} returns an "adjusted" R2 value computed using a
leave-one-out-adjusted posterior distribution. This is conceptually similar
to an adjusted/unbiased R2 estimate in classical regression modeling. See
\code{\link[=r2_bayes]{r2_bayes()}} for an "unadjusted" R2.
\cr \cr
Mixed models are not currently fully supported.
\cr \cr
\code{r2_loo_posterior()} is the actual workhorse for \code{r2_loo()} and
returns a posterior sample of LOO-adjusted Bayesian R2 values.
}
\examples{
if (require("rstanarm")) {
  model <- stan_glm(mpg ~ wt + cyl, data = mtcars, chains = 1, iter = 500, refresh = 0)
  r2_loo(model)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_somers.R
\name{r2_somers}
\alias{r2_somers}
\title{Somers' Dxy rank correlation for binary outcomes}
\usage{
r2_somers(model)
}
\arguments{
\item{model}{A logistic regression model.}
}
\value{
A named vector with the R2 value.
}
\description{
Calculates the Somers' Dxy rank correlation for logistic regression models.
}
\examples{
\dontrun{
if (require("correlation")) {
  model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
  r2_somers(model)
}
}

}
\references{
Somers, R. H. (1962). A new asymmetric measure of association for
ordinal variables. American Sociological Review. 27 (6).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/icc.R
\name{icc}
\alias{icc}
\alias{variance_decomposition}
\title{Intraclass Correlation Coefficient (ICC)}
\usage{
icc(model, by_group = FALSE, tolerance = 1e-05)

variance_decomposition(model, re_formula = NULL, robust = TRUE, ci = 0.95, ...)
}
\arguments{
\item{model}{A (Bayesian) mixed effects model.}

\item{by_group}{Logical, if \code{TRUE}, \code{icc()} returns the variance
components for each random-effects level (if there are multiple levels).
See 'Details'.}

\item{tolerance}{Tolerance for singularity check of random effects, to decide
whether to compute random effect variances or not. Indicates up to which
value the convergence result is accepted. The larger tolerance is, the
stricter the test will be. See \code{\link[performance:check_singularity]{performance::check_singularity()}}.}

\item{re_formula}{Formula containing group-level effects to be considered in
the prediction. If \code{NULL} (default), include all group-level effects.
Else, for instance for nested models, name a specific group-level effect
to calculate the variance decomposition for this group-level. See 'Details'
and \code{?brms::posterior_predict}.}

\item{robust}{Logical, if \code{TRUE}, the median instead of mean is used to
calculate the central tendency of the variances.}

\item{ci}{Credible interval level.}

\item{...}{Arguments passed down to \code{brms::posterior_predict()}.}
}
\value{
A list with two values, the adjusted and conditional ICC. For
\code{variance_decomposition()}, a list with two values, the decomposed
ICC as well as the credible intervals for this ICC.
}
\description{
This function calculates the intraclass-correlation coefficient (ICC) -
sometimes also called \emph{variance partition coefficient} (VPC) - for mixed
effects models. The ICC can be calculated for all models supported by
\code{insight::get_variance()}. For models fitted with the
\strong{brms}-package, \code{icc()} might fail due to the large variety of
models and families supported by the \strong{brms}-package. In such cases, an
alternative to the ICC is the \code{variance_decomposition()}, which is based
on the posterior predictive distribution (see 'Details').
}
\details{
\subsection{Interpretation}{
The ICC can be interpreted as \dQuote{the proportion of the variance
explained by the grouping structure in the population}. The grouping
structure entails that measurements are organized into groups (e.g., test
scores in a school can be grouped by classroom if there are multiple
classrooms and each classroom was administered the same test) and ICC indexes
how strongly measurements in the same group resemble each other. This index
goes from 0, if the grouping conveys no information, to 1, if all
observations in a group are identical (Gelman \& Hill, 2007, p. 258). In
other word, the ICC \dQuote{can also be interpreted as the expected
correlation between two randomly drawn units that are in the same group}
\cite{(Hox 2010: 15)}, although this definition might not apply to mixed
models with more complex random effects structures.
}
\subsection{Calculation}{
The ICC is calculated by dividing the random effect variance,
\ifelse{html}{\out{&sigma;<sup>2</sup><sub>i</sub>}}{\eqn{\sigma^2_i}}, by
the total variance, i.e. the sum of the random effect variance and the
residual variance, \ifelse{html}{\out{&sigma;<sup>2</sup><sub>&epsilon;</sub>}}{\eqn{\sigma^2_\epsilon}}.
}
\subsection{Adjusted and conditional ICC}{
\code{icc()} calculates an adjusted and conditional ICC, which both take all
sources of uncertainty (i.e. of \emph{all random effects}) into account.
While the \emph{adjusted ICC} only relates to the random effects, the
\emph{conditional ICC} also takes the fixed effects variances into account
(see \cite{Nakagawa et al. 2017}). Typically, the \emph{adjusted} ICC is of
interest when the analysis of random effects is of interest. \code{icc()}
returns a meaningful ICC also for more complex random effects structures,
like models with random slopes or nested design (more than two levels) and
is applicable for models with other distributions than Gaussian. For more
details on the computation of the variances, see
\code{?insight::get_variance}.
}
\subsection{ICC for unconditional and conditional models}{
Usually, the ICC is calculated for the null model ("unconditional model").
However, according to \cite{Raudenbush and Bryk (2002)} or
\cite{Rabe-Hesketh and Skrondal (2012)} it is also feasible to compute the
ICC for full models with covariates ("conditional models") and compare how
much, e.g., a level-2 variable explains the portion of variation in the
grouping structure (random intercept).
}
\subsection{ICC for specific group-levels}{
The proportion of variance for specific levels related to the overall model
can be computed by setting \code{by_group = TRUE}. The reported ICC is
the variance for each (random effect) group compared to the total
variance of the model. For mixed models with a simple random intercept,
this is identical to the classical (adjusted) ICC.
}
\subsection{Variance decomposition for brms-models}{
If \code{model} is of class \code{brmsfit}, \code{icc()} might fail due to
the large variety of models and families supported by the \strong{brms}
package. In such cases, \code{variance_decomposition()} is an alternative
ICC measure. The function calculates a variance decomposition based on the
posterior predictive distribution. In this case, first, the draws from the
posterior predictive distribution \emph{not conditioned} on group-level
terms (\code{posterior_predict(..., re_formula = NA)}) are calculated as
well as draws from this distribution \emph{conditioned} on \emph{all random
effects} (by default, unless specified else in \code{re_formula}) are taken.
Then, second, the variances for each of these draws are calculated. The
"ICC" is then the ratio between these two variances. This is the recommended
way to analyse random-effect-variances for non-Gaussian models. It is then
possible to compare variances across models, also by specifying different
group-level terms via the \code{re_formula}-argument.
\cr \cr
Sometimes, when the variance of the posterior predictive distribution is
very large, the variance ratio in the output makes no sense, e.g. because
it is negative. In such cases, it might help to use \code{robust = TRUE}.
}
}
\examples{
if (require("lme4")) {
  model <- lmer(Sepal.Length ~ Petal.Length + (1 | Species), data = iris)
  icc(model)
}

# ICC for specific group-levels
if (require("lme4")) {
  data(sleepstudy)
  set.seed(12345)
  sleepstudy$grp <- sample(1:5, size = 180, replace = TRUE)
  sleepstudy$subgrp <- NA
  for (i in 1:5) {
    filter_group <- sleepstudy$grp == i
    sleepstudy$subgrp[filter_group] <-
      sample(1:30, size = sum(filter_group), replace = TRUE)
  }
  model <- lmer(
    Reaction ~ Days + (1 | grp / subgrp) + (1 | Subject),
    data = sleepstudy
  )
  icc(model, by_group = TRUE)
}
}
\references{
\itemize{
\item Hox, J. J. (2010). Multilevel analysis: techniques and applications
(2nd ed). New York: Routledge.
\item Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The
coefficient of determination R2 and intra-class correlation coefficient from
generalized linear mixed-effects models revisited and expanded. Journal of
The Royal Society Interface, 14(134), 20170213. \doi{10.1098/rsif.2017.0213}
\item Rabe-Hesketh, S., & Skrondal, A. (2012). Multilevel and longitudinal
modeling using Stata (3rd ed). College Station, Tex: Stata Press
Publication.
\item Raudenbush, S. W., & Bryk, A. S. (2002). Hierarchical linear models:
applications and data analysis methods (2nd ed). Thousand Oaks: Sage
Publications.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.lavaan.R
\name{model_performance.lavaan}
\alias{model_performance.lavaan}
\title{Performance of lavaan SEM / CFA Models}
\usage{
\method{model_performance}{lavaan}(model, metrics = "all", verbose = TRUE, ...)
}
\arguments{
\item{model}{A \pkg{lavaan} model.}

\item{metrics}{Can be \code{"all"} or a character vector of metrics to be
computed (some of \code{c("Chi2", "Chi2_df", "p_Chi2", "Baseline", "Baseline_df", "p_Baseline", "GFI", "AGFI", "NFI", "NNFI", "CFI", "RMSEA", "RMSEA_CI_low", "RMSEA_CI_high", "p_RMSEA", "RMR", "SRMR", "RFI", "PNFI", "IFI", "RNI", "Loglikelihood", "AIC", "BIC", "BIC_adjusted")}).}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame (with one row) and one column per "index" (see
\code{metrics}).
}
\description{
Compute indices of model performance for SEM or CFA models from the
\pkg{lavaan} package.
}
\details{
\subsection{Indices of fit}{
\itemize{
\item \strong{Chisq}: The model Chi-squared assesses overall fit and the
discrepancy between the sample and fitted covariance matrices. Its p-value
should be > .05 (i.e., the hypothesis of a perfect fit cannot be
rejected). However, it is quite sensitive to sample size.

\item \strong{GFI/AGFI}: The (Adjusted) Goodness of Fit is the proportion
of variance accounted for by the estimated population covariance.
Analogous to R2. The GFI and the AGFI should be > .95 and > .90,
respectively.

\item \strong{NFI/NNFI/TLI}: The (Non) Normed Fit Index. An NFI of 0.95,
indicates the model of interest improves the fit by 95\\% relative to the
null model. The NNFI (also called the Tucker Lewis index; TLI) is
preferable for smaller samples. They should be > .90 (Byrne, 1994) or >
.95 (Schumacker & Lomax, 2004).

\item \strong{CFI}: The Comparative Fit Index is a revised form of NFI.
Not very sensitive to sample size (Fan, Thompson, & Wang, 1999). Compares
the fit of a target model to the fit of an independent, or null, model. It
should be > .90.

\item \strong{RMSEA}: The Root Mean Square Error of Approximation is a
parsimony-adjusted index. Values closer to 0 represent a good fit. It
should be < .08 or < .05. The p-value printed with it tests the hypothesis
that RMSEA is less than or equal to .05 (a cutoff sometimes used for good
fit), and thus should be not significant.

\item \strong{RMR/SRMR}: the (Standardized) Root Mean Square Residual
represents the square-root of the difference between the residuals of the
sample covariance matrix and the hypothesized model. As the RMR can be
sometimes hard to interpret, better to use SRMR. Should be < .08.

\item \strong{RFI}: the Relative Fit Index, also known as RHO1, is not
guaranteed to vary from 0 to 1. However, RFI close to 1 indicates a good
fit.

\item \strong{IFI}: the Incremental Fit Index (IFI) adjusts the Normed Fit
Index (NFI) for sample size and degrees of freedom (Bollen's, 1989). Over
0.90 is a good fit, but the index can exceed 1.

\item \strong{PNFI}: the Parsimony-Adjusted Measures Index. There is no
commonly agreed-upon cutoff value for an acceptable model for this index.
Should be > 0.50. }

See the documentation for \code{?lavaan::fitmeasures}.
}

\subsection{What to report}{
Kline (2015) suggests that at a minimum the following indices should be
reported: The model \strong{chi-square}, the \strong{RMSEA}, the \strong{CFI}
and the \strong{SRMR}.
}
}
\examples{
# Confirmatory Factor Analysis (CFA) ---------
if (require("lavaan")) {
  structure <- " visual  =~ x1 + x2 + x3
                 textual =~ x4 + x5 + x6
                 speed   =~ x7 + x8 + x9 "
  model <- lavaan::cfa(structure, data = HolzingerSwineford1939)
  model_performance(model)
}
}
\references{
\itemize{
\item Byrne, B. M. (1994). Structural equation modeling with EQS and
EQS/Windows. Thousand Oaks, CA: Sage Publications.

\item Tucker, L. R., \& Lewis, C. (1973). The reliability coefficient for
maximum likelihood factor analysis. Psychometrika, 38, 1-10.

\item Schumacker, R. E., \& Lomax, R. G. (2004). A beginner's guide to
structural equation modeling, Second edition. Mahwah, NJ: Lawrence Erlbaum
Associates.

\item Fan, X., B. Thompson, \& L. Wang (1999). Effects of sample size,
estimation method, and model specification on structural equation modeling
fit indexes. Structural Equation Modeling, 6, 56-83.

\item Kline, R. B. (2015). Principles and practice of structural equation
modeling. Guilford publications.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_autocorrelation.R
\name{check_autocorrelation}
\alias{check_autocorrelation}
\alias{check_autocorrelation.default}
\title{Check model for independence of residuals.}
\usage{
check_autocorrelation(x, ...)

\method{check_autocorrelation}{default}(x, nsim = 1000, ...)
}
\arguments{
\item{x}{A model object.}

\item{...}{Currently not used.}

\item{nsim}{Number of simulations for the Durbin-Watson-Test.}
}
\value{
Invisibly returns the p-value of the test statistics. A p-value < 0.05
indicates autocorrelated residuals.
}
\description{
Check model for independence of residuals, i.e. for autocorrelation
of error terms.
}
\details{
Performs a Durbin-Watson-Test to check for autocorrelated residuals.
In case of autocorrelation, robust standard errors return more accurate
results for the estimates, or maybe a mixed model with error term for the
cluster groups should be used.
}
\examples{
m <- lm(mpg ~ wt + cyl + gear + disp, data = mtcars)
check_autocorrelation(m)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_mckelvey.R
\name{r2_mckelvey}
\alias{r2_mckelvey}
\title{McKelvey & Zavoinas R2}
\usage{
r2_mckelvey(model)
}
\arguments{
\item{model}{Generalized linear model.}
}
\value{
The R2 value.
}
\description{
Calculates McKelvey & Zavoinas pseudo R2.
}
\details{
McKelvey & Zavoinas R2 is based on the explained variance,
where the variance of the predicted response is divided by the sum
of the variance of the predicted response and residual variance.
For binomial models, the residual variance is either \code{pi^2/3}
for logit-link and 1 for probit-link. For poisson-models, the
residual variance is based on log-normal approximation, similar to
the \emph{distribution-specific variance} as described in
\code{?insight::get_variance}.
}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial:
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12) #
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
model <- glm(counts ~ outcome + treatment, family = poisson())

r2_mckelvey(model)
}
\references{
\itemize{
\item McKelvey, R., Zavoina, W. (1975), "A Statistical Model for the Analysis of Ordinal Level Dependent Variables", Journal of Mathematical Sociology 4, S. 103–120.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_hosmer.R
\name{performance_hosmer}
\alias{performance_hosmer}
\title{Hosmer-Lemeshow goodness-of-fit test}
\usage{
performance_hosmer(model, n_bins = 10)
}
\arguments{
\item{model}{A \code{glm}-object with binomial-family.}

\item{n_bins}{Numeric, the number of bins to divide the data.}
}
\value{
An object of class \code{hoslem_test} with following values:
\code{chisq}, the Hosmer-Lemeshow chi-squared statistic; \code{df}, degrees
of freedom and \code{p.value} the p-value for the goodness-of-fit test.
}
\description{
Check model quality of logistic regression models.
}
\details{
A well-fitting model shows \emph{no} significant difference between
the model and the observed data, i.e. the reported p-value should be
greater than 0.05.
}
\examples{
model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
performance_hosmer(model)
}
\references{
Hosmer, D. W., & Lemeshow, S. (2000). Applied Logistic Regression. Hoboken,
NJ, USA: John Wiley & Sons, Inc. \doi{10.1002/0471722146}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_sphericity.R
\name{check_sphericity}
\alias{check_sphericity}
\title{Check model for violation of sphericity}
\usage{
check_sphericity(x, ...)
}
\arguments{
\item{x}{A model object.}

\item{...}{Arguments passed to \code{car::Anova}.}
}
\value{
Invisibly returns the p-values of the test statistics. A p-value <
0.05 indicates a violation of sphericity.
}
\description{
Check model for violation of sphericity
}
\examples{
if (require("car")) {
  soils.mod <- lm(
    cbind(pH, N, Dens, P, Ca, Mg, K, Na, Conduc) ~ Block + Contour * Depth,
    data = Soils
  )

  check_sphericity(Manova(soils.mod))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_efron.R
\name{r2_efron}
\alias{r2_efron}
\title{Efron's R2}
\usage{
r2_efron(model)
}
\arguments{
\item{model}{Generalized linear model.}
}
\value{
The R2 value.
}
\description{
Calculates Efron's pseudo R2.
}
\details{
Efron's R2 is calculated by taking the sum of the squared model residuals,
divided by the total variability in the dependent variable. This R2 equals
the squared correlation between the predicted values and actual values,
however, note that model residuals from generalized linear models are not
generally comparable to those of OLS.
}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial:
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12) #
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
model <- glm(counts ~ outcome + treatment, family = poisson())

r2_efron(model)
}
\references{
\itemize{
\item Efron, B. (1978). Regression and ANOVA with zero-one data: Measures of
residual variation. Journal of the American Statistical Association, 73,
113-121.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_pcp.R
\name{performance_pcp}
\alias{performance_pcp}
\title{Percentage of Correct Predictions}
\usage{
performance_pcp(model, ci = 0.95, method = "Herron", verbose = TRUE)
}
\arguments{
\item{model}{Model with binary outcome.}

\item{ci}{The level of the confidence interval.}

\item{method}{Name of the method to calculate the PCP (see 'Details').
Default is \code{"Herron"}. May be abbreviated.}

\item{verbose}{Toggle off warnings.}
}
\value{
A list with several elements: the percentage of correct predictions
of the full and the null model, their confidence intervals, as well as the
chi-squared and p-value from the Likelihood-Ratio-Test between the full and
null model.
}
\description{
Percentage of correct predictions (PCP) for models
with binary outcome.
}
\details{
\code{method = "Gelman-Hill"} (or \code{"gelman_hill"}) computes the
PCP based on the proposal from \cite{Gelman and Hill 2017, 99}, which is
defined as the proportion of cases for which the deterministic prediction
is wrong, i.e. the proportion where the predicted probability is above 0.5,
although y=0 (and vice versa) (see also \cite{Herron 1999, 90}).
\cr \cr
\code{method = "Herron"} (or \code{"herron"}) computes a modified version
of the PCP (\cite{Herron 1999, 90-92}), which is the sum of predicted
probabilities, where y=1, plus the sum of 1 - predicted probabilities,
where y=0, divided by the number of observations. This approach is said to
be more accurate.
\cr \cr
The PCP ranges from 0 to 1, where values closer to 1 mean that the model
predicts the outcome better than models with an PCP closer to 0. In general,
the PCP should be above 0.5 (i.e. 50\\%), the closer to one, the better.
Furthermore, the PCP of the full model should be considerably above
the null model's PCP.
\cr \cr
The likelihood-ratio test indicates whether the model has a significantly
better fit than the null-model (in such cases, p < 0.05).
}
\examples{
data(mtcars)
m <- glm(formula = vs ~ hp + wt, family = binomial, data = mtcars)
performance_pcp(m)
performance_pcp(m, method = "Gelman-Hill")
}
\references{
\itemize{
\item Herron, M. (1999). Postestimation Uncertainty in Limited Dependent
Variable Models. Political Analysis, 8, 83–98.
\item Gelman, A., & Hill, J. (2007). Data analysis using regression and
multilevel/hierarchical models. Cambridge; New York: Cambridge University
Press, 99.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cronbachs_alpha.R
\name{cronbachs_alpha}
\alias{cronbachs_alpha}
\title{Cronbach's Alpha for Items or Scales}
\usage{
cronbachs_alpha(x)
}
\arguments{
\item{x}{A matrix or a data frame.}
}
\value{
The Cronbach's Alpha value for \code{x}.
}
\description{
Compute various measures of internal consistencies
for tests or item-scales of questionnaires.
}
\details{
The Cronbach's Alpha value for \code{x}. A value closer to 1
indicates greater internal consistency, where usually following
rule of thumb is applied to interpret the results:
\ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} < 0.5 is unacceptable,
0.5 < \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} < 0.6 is poor,
0.6 < \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} < 0.7 is questionable,
0.7 < \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} < 0.8 is acceptable,
and everything > 0.8 is good or excellent.
}
\examples{
data(mtcars)
x <- mtcars[, c("cyl", "gear", "carb", "hp")]
cronbachs_alpha(x)
}
\references{
Bland, J. M., \& Altman, D. G. Statistics notes: Cronbach's
alpha. BMJ 1997;314:572. 10.1136/bmj.314.7080.572
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/item_difficulty.R
\name{item_difficulty}
\alias{item_difficulty}
\title{Difficulty of Questionnaire Items}
\usage{
item_difficulty(x)
}
\arguments{
\item{x}{Depending on the function, \code{x} may be a \code{matrix} as
returned by the \code{cor()}-function, or a data frame
with items (e.g. from a test or questionnaire).}
}
\value{
A data frame with three columns: The name(s) of the item(s), the item
difficulties for each item, and the ideal item difficulty.
}
\description{
Compute various measures of internal consistencies
for tests or item-scales of questionnaires.
}
\details{
This function calculates the item difficulty, which should
range between 0.2 and 0.8. Lower values are a signal for
more difficult items, while higher values close to one
are a sign for easier items. The ideal value for item difficulty
is \code{p + (1 - p) / 2}, where \code{p = 1 / max(x)}. In most
cases, the ideal item difficulty lies between 0.5 and 0.8.
}
\examples{
data(mtcars)
x <- mtcars[, c("cyl", "gear", "carb", "hp")]
item_difficulty(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_collinearity.R
\name{check_collinearity}
\alias{check_collinearity}
\alias{multicollinearity}
\alias{check_collinearity.default}
\alias{check_collinearity.glmmTMB}
\title{Check for multicollinearity of model terms}
\usage{
check_collinearity(x, ...)

multicollinearity(x, ...)

\method{check_collinearity}{default}(x, verbose = TRUE, ...)

\method{check_collinearity}{glmmTMB}(
  x,
  component = c("all", "conditional", "count", "zi", "zero_inflated"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{A model object (that should at least respond to \code{vcov()},
and if possible, also to \code{model.matrix()} - however, it also should
work without \code{model.matrix()}).}

\item{...}{Currently not used.}

\item{verbose}{Toggle off warnings or messages.}

\item{component}{For models with zero-inflation component, multicollinearity
can be checked for the conditional model (count component,
\code{component = "conditional"} or \code{component = "count"}),
zero-inflation component (\code{component = "zero_inflated"} or
\code{component = "zi"}) or both components (\code{component = "all"}).
Following model-classes are currently supported: \code{hurdle},
\code{zeroinfl}, \code{zerocount}, \code{MixMod} and \code{glmmTMB}.}
}
\value{
A data frame with three columns: The name of the model term, the
variance inflation factor and the factor by which the standard error
is increased due to possible correlation with other terms.
}
\description{
\code{check_collinearity()} checks regression models for
multicollinearity by calculating the variance inflation factor (VIF).
\code{multicollinearity()} is an alias for \code{check_collinearity()}.
(When printed, VIF are also translated to Tolerance values, where
\code{tolerance = 1/vif}.)
}
\details{
\subsection{Multicollinearity}{
Multicollinearity should not be confused with a raw strong correlation
between predictors. What matters is the association between one or more
predictor variables, \emph{conditional on the other variables in the
model}. In a nutshell, multicollinearity means that once you know the
effect of one predictor, the value of knowing the other predictor is rather
low. Thus, one of the predictors doesn't help much in terms of better
understanding the model or predicting the outcome. As a consequence, if
multicollinearity is a problem, the model seems to suggest that the
predictors in question don't seems to be reliably associated with the
outcome (low estimates, high standard errors), although these predictors
actually are strongly associated with the outcome, i.e. indeed might have
strong effect (\cite{McElreath 2020, chapter 6.1}).
\cr \cr
Multicollinearity might arise when a third, unobserved variable has a causal
effect on each of the two predictors that are associated with the outcome.
In such cases, the actual relationship that matters would be the association
between the unobserved variable and the outcome.
\cr \cr
Remember: \dQuote{Pairwise correlations are not the problem. It is the
conditional associations - not correlations - that matter.}
(\cite{McElreath 2020, p. 169})
}

\subsection{Interpretation of the Variance Inflation Factor}{
The variance inflation factor is a measure to analyze the magnitude of
multicollinearity of model terms. A VIF less than 5 indicates a low
correlation of that predictor with other predictors. A value between 5 and
10 indicates a moderate correlation, while VIF values larger than 10 are a
sign for high, not tolerable correlation of model predictors (\cite{James
et al. 2013}). The \emph{Increased SE} column in the output indicates how
much larger the standard error is due to the association with other
predictors conditional on the remaining variables in the model.
}

\subsection{Multicollinearity and Interaction Terms}{
If interaction terms are included in a model, high VIF values are expected.
This portion of multicollinearity among the component terms of an
interaction is also called "inessential ill-conditioning", which leads to
inflated VIF values that are typically seen for models with interaction
terms \cite{(Francoeur 2013)}.
}
}
\note{
There is also a \href{https://easystats.github.io/see/articles/performance.html}{\code{plot()}-method} implemented in the \href{https://easystats.github.io/see/}{\pkg{see}-package}.
}
\examples{
m <- lm(mpg ~ wt + cyl + gear + disp, data = mtcars)
check_collinearity(m)

# plot results
if (require("see")) {
  x <- check_collinearity(m)
  plot(x)
}
}
\references{
\itemize{
\item Francoeur, R. B. (2013). Could Sequential Residual Centering Resolve
Low Sensitivity in Moderated Regression? Simulations and Cancer Symptom
Clusters. Open Journal of Statistics, 03(06), 24-44.

\item James, G., Witten, D., Hastie, T., & Tibshirani, R. (eds.). (2013).
An introduction to statistical learning: with applications in R. New York:
Springer.

\item McElreath, R. (2020). Statistical rethinking: A Bayesian course with
examples in R and Stan. 2nd edition. Chapman and Hall/CRC.

\item Vanhove, J. (2019). Collinearity isn't a disease that needs curing.
\href{https://janhove.github.io/analysis/2019/09/11/collinearity}{webpage}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_score.R
\name{performance_score}
\alias{performance_score}
\title{Proper Scoring Rules}
\usage{
performance_score(model, verbose = TRUE)
}
\arguments{
\item{model}{Model with binary or count outcome.}

\item{verbose}{Toggle off warnings.}
}
\value{
A list with three elements, the logarithmic, quadratic/Brier and spherical score.
}
\description{
Calculates the logarithmic, quadratic/Brier and spherical score
from a model with binary or count outcome.
}
\details{
Proper scoring rules can be used to evaluate the quality of model
predictions and model fit. \code{performance_score()} calculates the logarithmic,
quadratic/Brier and spherical scoring rules. The spherical rule takes values
in the interval \verb{[0, 1]}, with values closer to 1 indicating a more
accurate model, and the logarithmic rule in the interval \verb{[-Inf, 0]},
with values closer to 0 indicating a more accurate model.
\cr \cr
For \code{stan_lmer()} and \code{stan_glmer()} models, the predicted values
are based on \code{posterior_predict()}, instead of \code{predict()}. Thus,
results may differ more than expected from their non-Bayesian counterparts
in \strong{lme4}.
}
\note{
Code is partially based on \href{https://drizopoulos.github.io/GLMMadaptive/reference/scoring_rules.html}{GLMMadaptive::scoring_rules()}.
}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
model <- glm(counts ~ outcome + treatment, family = poisson())

performance_score(model)
\dontrun{
if (require("glmmTMB")) {
  data(Salamanders)
  model <- glmmTMB(
    count ~ spp + mined + (1 | site),
    zi =  ~ spp + mined,
    family = nbinom2(),
    data = Salamanders
  )

  performance_score(model)
}
}
}
\references{
Carvalho, A. (2016). An overview of applications of proper scoring rules. Decision Analysis 13, 223–242. \doi{10.1287/deca.2016.0337}
}
\seealso{
\code{\link[=performance_logloss]{performance_logloss()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_kl.R
\name{r2_kullback}
\alias{r2_kullback}
\title{Kullback-Leibler R2}
\usage{
r2_kullback(model, adjust = TRUE)
}
\arguments{
\item{model}{A generalized linear model.}

\item{adjust}{Logical, if \code{TRUE} (the default), the adjusted R2 value is
returned.}
}
\value{
A named vector with the R2 value.
}
\description{
Calculates the Kullback-Leibler-divergence-based
R2 for generalized linear models.
}
\examples{
model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
r2_kullback(model)
}
\references{
Cameron, A. C. and Windmeijer, A. G. (1997) An R-squared measure of goodness
of fit for some common nonlinear regression models. Journal of Econometrics,
77: 329-342.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_bf.R, R/test_likelihoodratio.R,
%   R/test_performance.R, R/test_vuong.R, R/test_wald.R
\name{test_bf}
\alias{test_bf}
\alias{test_bf.default}
\alias{test_likelihoodratio}
\alias{performance_lrt}
\alias{test_lrt}
\alias{test_performance}
\alias{test_vuong}
\alias{test_wald}
\title{Test if Models are Different}
\usage{
test_bf(...)

\method{test_bf}{default}(..., text_length = NULL)

test_likelihoodratio(..., estimator = "ML")

performance_lrt(..., estimator = "ML")

test_lrt(..., estimator = "ML")

test_performance(..., reference = 1)

test_vuong(...)

test_wald(...)
}
\arguments{
\item{...}{Multiple model objects.}

\item{text_length}{Numeric, length (number of chars) of output lines.
\code{test_bf()} describes models by their formulas, which can lead to
overly long lines in the output. \code{text_length} fixes the length of
lines to a specified limit.}

\item{estimator}{Applied when comparing regression models using
\code{test_likelihoodratio()}. Corresponds to the different estimators for
the standard deviation of the errors. If \code{estimator="OLS"} (default),
then it uses the same method as \code{anova(..., test="LRT")} implemented
in base R, i.e., scaling by n-k (the unbiased OLS estimator) and using this
estimator under the alternative hypothesis. If \code{estimator="ML"}, which
is for instance used by \code{lrtest(...)} in package \pkg{lmtest}, the
scaling is done by n (the biased ML estimator) and the estimator under the
null hypothesis. In moderately large samples, the differences should be
negligible, but it is possible that OLS would perform slightly better in
small samples with Gaussian errors.}

\item{reference}{This only applies when models are non-nested, and determines
which model should be taken as a reference, against which all the other
models are tested.}
}
\value{
A data frame containing the relevant indices.
}
\description{
Testing whether models are "different" in terms of accuracy or explanatory
power is a delicate and often complex procedure, with many limitations and
prerequisites. Moreover, many tests exist, each coming with its own
interpretation, and set of strengths and weaknesses.
\cr \cr
The \code{test_performance()} function runs the most relevant and appropriate
tests based on the type of input (for instance, whether the models are
\emph{nested} or not). However, it still requires the user to understand what the
tests are and what they do in order to prevent their misinterpretation. See
the \strong{details} section for more information regarding the different tests
and their interpretation.
}
\details{
\subsection{Nested vs. Non-nested Models}{
Model's "nesting" is an important concept of models comparison. Indeed, many
tests only make sense when the models are \emph{"nested",} i.e., when their
predictors are nested. This means that all the predictors of a model are
contained within the predictors of a larger model (sometimes referred to as
the encompassing model). For instance, \code{model1 (y ~ x1 + x2)} is
"nested" within \code{model2 (y ~ x1 + x2 + x3)}. Usually, people have a list
of nested models, for instance \code{m1 (y ~ 1)}, \code{m2 (y ~ x1)},
\code{m3 (y ~ x1 + x2)}, \code{m4 (y ~ x1 + x2 + x3)}, and it is conventional
that they are "ordered" from the smallest to largest, but it is up to the
user to reverse the order from largest to smallest. The test then shows
whether a more parsimonious model, or whether adding a predictor, results in
a significant difference in the model's performance. In this case, models are
usually compared \emph{sequentially}: m2 is tested against m1, m3 against m2,
m4 against m3, etc.
\cr\cr
Two models are considered as \emph{"non-nested"} if their predictors are
different. For instance, \code{model1 (y ~ x1 + x2)} and `model2 (y ~ x3
\itemize{
\item x4)\verb{. In the case of non-nested models, all models are usually compared against the same *reference* model (by default, the first of the list). \\cr\\cr Nesting is detected via the }insight::is_nested_models()` function.
Note that, apart from the nesting, in order for the tests to be valid,
other requirements have often to be the fulfilled. For instance, outcome
variables (the response) must be the same. You cannot meaningfully test
whether apples are significantly different from oranges!
}
}

\subsection{Tests Description}{

\itemize{
\item \strong{Bayes factor for Model Comparison} - \code{test_bf()}: If all
models were fit from the same data, the returned \code{BF} shows the Bayes
Factor (see \code{bayestestR::bayesfactor_models()}) for each model against
the reference model (which depends on whether the models are nested or
not). Check out
\href{https://easystats.github.io/bayestestR/articles/bayes_factors.html#bayesfactor_models}{this vignette} for more details.

\item \strong{Wald's F-Test} - \code{test_wald()}: The Wald test is a rough
approximation of the Likelihood Ratio Test. However, it is more applicable
than the LRT: you can often run a Wald test in situations where no other
test can be run. Importantly, this test only makes statistical sense if the
models are nested.\cr Note: this test is also available in base R through the
\code{\link[=anova]{anova()}} function. It returns an \code{F-value} column
as a statistic and its associated \code{p-value}.

\item \strong{Likelihood Ratio Test (LRT)} - \code{test_likelihoodratio()}:
The LRT tests which model is a better (more likely) explanation of the
data. Likelihood-Ratio-Test (LRT) gives usually somewhat close results (if
not equivalent) to the Wald test and, similarly, only makes sense for
nested models. However, Maximum likelihood tests make stronger assumptions
than method of moments tests like the F-test, and in turn are more
efficient. Agresti (1990) suggests that you should use the LRT instead of
the Wald test for small sample sizes (under or about 30) or if the
parameters are large.\cr Note: for regression models, this is similar to
\code{anova(..., test="LRT")} (on models) or \code{lmtest::lrtest(...)},
depending on the \code{estimator} argument. For \code{lavaan} models (SEM,
CFA), the function calls \code{lavaan::lavTestLRT()}.

\item \strong{Vuong's Test} - \code{test_vuong()}: Vuong's (1989) test can
be used both for nested and non-nested models, and actually consists of two
tests.
\itemize{
\item The \strong{Test of Distinguishability} (the \code{Omega2} column and
its associated p-value) indicates whether or not the models can possibly be
distinguished on the basis of the observed data. If its p-value is
significant, it means the models are distinguishable.
\item The \strong{Robust Likelihood Test} (the \code{LR} column and its
associated p-value) indicates whether each model fits better than the
reference model. If the models are nested, then the test works as a robust
LRT. The code for this function is adapted from the \code{nonnest2}
package, and all credit go to their authors.}
}
}
}
\examples{
# Nested Models
# -------------
m1 <- lm(Sepal.Length ~ Petal.Width, data = iris)
m2 <- lm(Sepal.Length ~ Petal.Width + Species, data = iris)
m3 <- lm(Sepal.Length ~ Petal.Width * Species, data = iris)

test_performance(m1, m2, m3)

test_bf(m1, m2, m3)
test_wald(m1, m2, m3) # Equivalent to anova(m1, m2, m3)

# Equivalent to lmtest::lrtest(m1, m2, m3)
test_likelihoodratio(m1, m2, m3, estimator = "ML")

# Equivalent to anova(m1, m2, m3, test='LRT')
test_likelihoodratio(m1, m2, m3, estimator = "OLS")

test_vuong(m1, m2, m3) # nonnest2::vuongtest(m1, m2, nested=TRUE)

# Non-nested Models
# -----------------
m1 <- lm(Sepal.Length ~ Petal.Width, data = iris)
m2 <- lm(Sepal.Length ~ Petal.Length, data = iris)
m3 <- lm(Sepal.Length ~ Species, data = iris)

test_performance(m1, m2, m3)
test_bf(m1, m2, m3)
test_vuong(m1, m2, m3) # nonnest2::vuongtest(m1, m2)

# Tweak the output
# ----------------
test_performance(m1, m2, m3, include_formula = TRUE)


# SEM / CFA (lavaan objects)
# --------------------------
# Lavaan Models
if (require("lavaan")) {
  structure <- " visual  =~ x1 + x2 + x3
                 textual =~ x4 + x5 + x6
                 speed   =~ x7 + x8 + x9

                  visual ~~ textual + speed "
  m1 <- lavaan::cfa(structure, data = HolzingerSwineford1939)

  structure <- " visual  =~ x1 + x2 + x3
                 textual =~ x4 + x5 + x6
                 speed   =~ x7 + x8 + x9

                  visual ~~ 0 * textual + speed "
  m2 <- lavaan::cfa(structure, data = HolzingerSwineford1939)

  structure <- " visual  =~ x1 + x2 + x3
                 textual =~ x4 + x5 + x6
                 speed   =~ x7 + x8 + x9

                  visual ~~ 0 * textual + 0 * speed "
  m3 <- lavaan::cfa(structure, data = HolzingerSwineford1939)

  test_likelihoodratio(m1, m2, m3)

  # Different Model Types
  # ---------------------
  if (require("lme4") && require("mgcv")) {
    m1 <- lm(Sepal.Length ~ Petal.Length + Species, data = iris)
    m2 <- lmer(Sepal.Length ~ Petal.Length + (1 | Species), data = iris)
    m3 <- gam(Sepal.Length ~ s(Petal.Length, by = Species) + Species, data = iris)

    test_performance(m1, m2, m3)
  }
}
}
\references{
\itemize{
\item Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. Econometrica, 57, 307-333.

\item Merkle, E. C., You, D., & Preacher, K. (2016). Testing non-nested
structural equation models. Psychological Methods, 21, 151-163.
}
}
\seealso{
\code{\link[=compare_performance]{compare_performance()}} to compare
the performance indices of many different models.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_performance.lm.R
\name{model_performance.lm}
\alias{model_performance.lm}
\title{Performance of Regression Models}
\usage{
\method{model_performance}{lm}(model, metrics = "all", verbose = TRUE, ...)
}
\arguments{
\item{model}{A model.}

\item{metrics}{Can be \code{"all"}, \code{"common"} or a character vector of
metrics to be computed (some of \code{c("AIC", "AICc", "BIC", "R2", "R2_adj", "RMSE", "SIGMA", "LOGLOSS", "PCP", "SCORE")}). \code{"common"}
will compute AIC, BIC, R2 and RMSE.}

\item{verbose}{Toggle off warnings.}

\item{...}{Arguments passed to or from other methods.}
}
\value{
A data frame (with one row) and one column per "index" (see \code{metrics}).
}
\description{
Compute indices of model performance for regression models.
}
\details{
Depending on \code{model}, following indices are computed:
\itemize{
\item{\strong{AIC}} {Akaike's Information Criterion, see \code{?stats::AIC}}
\item{\strong{AICc}} {Second-order (or small sample) AIC with a correction for small sample sizes}
\item{\strong{BIC}} {Bayesian Information Criterion, see \code{?stats::BIC}}
\item{\strong{R2}} {r-squared value, see \code{\link[=r2]{r2()}}}
\item{\strong{R2_adj}} {adjusted r-squared, see \code{\link[=r2]{r2()}}}
\item{\strong{RMSE}} {root mean squared error, see \code{\link[=performance_rmse]{performance_rmse()}}}
\item{\strong{SIGMA}} {residual standard deviation, see \code{\link[insight:get_sigma]{insight::get_sigma()}}}
\item{\strong{LOGLOSS}} {Log-loss, see \code{\link[=performance_logloss]{performance_logloss()}}}
\item{\strong{SCORE_LOG}} {score of logarithmic proper scoring rule, see \code{\link[=performance_score]{performance_score()}}}
\item{\strong{SCORE_SPHERICAL}} {score of spherical proper scoring rule, see \code{\link[=performance_score]{performance_score()}}}
\item{\strong{PCP}} {percentage of correct predictions, see \code{\link[=performance_pcp]{performance_pcp()}}}
}
}
\examples{
model <- lm(mpg ~ wt + cyl, data = mtcars)
model_performance(model)

model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
model_performance(model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/looic.R
\name{looic}
\alias{looic}
\title{LOO-related Indices for Bayesian regressions.}
\usage{
looic(model, verbose = TRUE)
}
\arguments{
\item{model}{A Bayesian regression model.}

\item{verbose}{Toggle off warnings.}
}
\value{
A list with four elements, the ELPD, LOOIC and their standard errors.
}
\description{
Compute LOOIC (leave-one-out cross-validation (LOO) information
criterion) and ELPD (expected log predictive density) for Bayesian
regressions. For LOOIC and ELPD, smaller and larger values are respectively
indicative of a better fit.
}
\examples{
if (require("rstanarm")) {
  model <- stan_glm(mpg ~ wt + cyl, data = mtcars, chains = 1, iter = 500, refresh = 0)
  looic(model)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_mcfadden.R
\name{r2_mcfadden}
\alias{r2_mcfadden}
\title{McFadden's R2}
\usage{
r2_mcfadden(model, ...)
}
\arguments{
\item{model}{Generalized linear or multinomial logit (\code{mlogit}) model.}

\item{...}{Currently not used.}
}
\value{
For most models, a list with McFadden's R2 and adjusted McFadden's
R2 value. For some models, only McFadden's R2 is available.
}
\description{
Calculates McFadden's pseudo R2.
}
\examples{
if (require("mlogit")) {
  data("Fishing", package = "mlogit")
  Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")

  model <- mlogit(mode ~ price + catch, data = Fish)
  r2_mcfadden(model)
}
}
\references{
\itemize{
\item McFadden, D. (1987). Regression-based specification tests for the
multinomial logit model. Journal of econometrics, 34(1-2), 63-82.

\item McFadden, D. (1973). Conditional logit analysis of qualitative choice
behavior.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_rmse.R
\name{performance_rmse}
\alias{performance_rmse}
\alias{rmse}
\title{Root Mean Squared Error}
\usage{
performance_rmse(model, normalized = FALSE, verbose = TRUE)

rmse(model, normalized = FALSE, verbose = TRUE)
}
\arguments{
\item{model}{A model.}

\item{normalized}{Logical, use \code{TRUE} if normalized rmse should be returned.}

\item{verbose}{Toggle off warnings.}
}
\value{
Numeric, the root mean squared error.
}
\description{
Compute root mean squared error for (mixed effects) models,
including Bayesian regression models.
}
\details{
The RMSE is the square root of the variance of the residuals and indicates
the absolute fit of the model to the data (difference between observed data
to model's predicted values). It can be interpreted as the standard
deviation of the unexplained variance, and is in the same units as the
response variable. Lower values indicate better model fit.
\cr \cr
The normalized RMSE is the proportion of the RMSE related to the
range of the response variable. Hence, lower values indicate
less residual variance.
}
\examples{
if (require("nlme")) {
  m <- lme(distance ~ age, data = Orthodont)

  # RMSE
  performance_rmse(m, normalized = FALSE)

  # normalized RMSE
  performance_rmse(m, normalized = TRUE)
}
}
