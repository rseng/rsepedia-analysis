
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

*Wow, no problems at all. :)*