<!-- README.md is generated from README.Rmd. Please edit that file -->

The mecor Package
=================

This package for R implements measurement error correction methods for
measurement error in a continuous covariate or outcome in a linear model
with a continuous outcome.

Installation
============

The package can be installed via

``` r
devtools::install_github("LindaNab/mecor", build_vignettes = TRUE)
```

Quick demo
==========

``` r
library(mecor)
# load the internal covariate validation study
data("vat", package = "mecor")
head(vat)
# correct the biased exposure-outcome association
mecor(ir_ln ~ MeasError(substitute = wc, reference = vat) + age + sex + tbf, data = vat, method = "standard")
```

More examples
=============

Browse the vignettes of the package for more information.

``` r
browseVignettes(package = "mecor")
```

References
==========

Key reference
-------------

-   Nab L, van Smeden M, Keogh RH, Groenwold RHH. mecor: an R package
    for measurement error correction in linear models with a continuous
    outcome. 2021:208:106238. 
    [doi:10.1016/j.cmpb.2021.106238](https://doi.org/10.1016/j.cmpb.2021.106238)

References to methods implemented in the package
------------------------------------------------

-   Bartlett JW, Stavola DBL, Frost C. Linear mixed models for
    replication data to efficiently allow for covariate measurement
    error. Statistics in Medicine. 2009:28(25):3158–3178.
    [doi:10.1002/sim.3713](https://doi.org/10.1002/sim.3713)

-   Buonaccorsi JP. Measurement error: Models, methods, and
    applications. 2010. Chapman & Hall/CRC, Boca Raton.

-   Carroll RJ, Ruppert D, Stefanski LA, Crainiceanu CM. Measurement
    error in non-linear models: A modern perspective. 2006, 2nd edition.
    Chapman & Hall/CRC, Boca Raton.

-   Keogh RH, Carroll RJ, Tooze JA, Kirkpatrick SI, Freedman LS.
    Statistical issues related to dietary intake as the response
    variable in intervention trials. Statistics in Medicine.
    2016:35(25):4493–4508.
    [doi:10.1002/sim.7011](https://doi.org/10.1002/sim.7011)

-   Keogh RH, White IR. A toolkit for measurement error correction, with
    a focus on nutritional epidemiology. Statistics in Medicine
    2014:33(12):2137–2155.
    [doi:10.1002/sim.6095](https://doi.org/10.1002/sim.6095)

-   Nab L, Groenwold RHH, Welsing PMJ, van Smeden M. Measurement error
    in continuous endpoints in randomised trials: Problems and
    solutions. Statistics in Medicine. 2019:38(27):5182-5196.
    [doi:10.1002/sim.8359](https://doi.org/10.1002/sim.8359)

-   Rosner B, Spiegelman D, Willett WC. Correction of logistic
    regression relative risk estimates and confidence intervals for
    measurement error: The case of multiple covariates measured with
    error. 1990:132(4):734-745.
    [doi:10.1093/oxfordjournals.aje.a115715](https://doi.org/10.1093/oxfordjournals.aje.a115715)

-   Rosner B, Spiegelman D, Willett WC. Correction of logistic
    regression relative risk estimates and confidence intervals for
    random within-person measurement error. American Journal of
    Epidemiology. 1992:136(11):1400-1413.
    [doi:10.1093/oxfordjournals.aje.a116453](https://doi.org/10.1093/oxfordjournals.aje.a116453)

-   Spiegelman D, Carroll RJ, Kipnis V. Efficient regression calibration
    for logistic regression in main study/internal validation study
    designs with an imperfect reference instrument. Statistics in
    Medicine. 2001:20(1):139-160.
    [doi:10.1002/1097-0258(20010115)20:1\<139::AID-SIM644\>3.0.CO;2-K](https://doi.org/10.1002/1097-0258(20010115)20:1%3C139::AID-SIM644%3E3.0.CO;2-K)
# mecor 0.9.0

* First formal release.
* Added a `NEWS.md` file to track changes to the package.

## Resubmission
This is a resubmission. In this version I have:

* Changed the doi in DESCRIPTION

## Test environments
* local OS x86_64-apple-darwin17.0, R 4.0.2
* win-builder (devel)
* win-builder (release)

## R CMD check results

0 errors | 0 warnings | 1 notes

There was 1 NOTE:
Possibly mis-spelled words in DESCRIPTION:
  Buonaccorsi (8:543, 8:1015)
  Crainiceanu (8:144, 8:1114)
  Fieller (8:987)
  JA (8:645)
  JW (8:397)
  Kipnis (8:247)
  RJ (8:114, 8:242, 8:635, 8:1084)
  Rosner (8:810, 8:898)
  Ruppert (8:118, 8:1088)
  Spiegelman (8:220, 8:820, 8:908)
  Stavola (8:401)
  Stefanski (8:129, 8:1099)
  Tooze (8:639)
  Willett (8:835, 8:923)

These are names and not misspelled.

## revdepchecks
No reverse dependencies

---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# The mecor Package
This package for R implements measurement error correction methods for 
measurement error in a continuous covariate or outcome in a linear model 
with a continuous outcome. 

# Installation
The package can be installed via
```r
devtools::install_github("LindaNab/mecor", build_vignettes = TRUE)
```

# Quick demo
```r
library(mecor)
# load the internal covariate validation study
data("icvs", package = "mecor")
head(icvs)
# correct the biased exposure-outcome association
mecor(Y ~ MeasError(X_star, reference = X) + Z, data = icvs, method = "standard")
```

# More examples
Browse the vignettes of the package for more information.
```r
browseVignettes(package = "mecor")
```


# References
## Key reference
- Nab L, van Smeden M, Keogh RH, Groenwold RHH. mecor: an R package for
measurement error correction in linear models with a continuous outcome.

## References to methods implemented in the package
- Bartlett JW, Stavola DBL, Frost C. Linear mixed models for replication data to efficiently allow for covariate measurement error. Statistics in Medicine. 2009:28(25):3158–3178. [doi:10.1002/sim.3713](https://doi.org/10.1002/sim.3713)

- Buonaccorsi JP. Measurement error: Models, methods, and applications. 2010. Chapman & Hall/CRC, Boca Raton.

- Carroll RJ, Ruppert D, Stefanski LA, Crainiceanu CM. Measurement error in non-linear models: A modern perspective. 2006, 2nd edition. Chapman & Hall/CRC, Boca Raton.

- Keogh RH, Carroll RJ, Tooze JA, Kirkpatrick SI, Freedman LS. Statistical issues related to dietary intake as the response variable in intervention trials. Statistics in Medicine. 2016:35(25):4493–4508. [doi:10.1002/sim.7011](https://doi.org/10.1002/sim.7011)

- Keogh RH, White IR. A toolkit for measurement error correction, with a focus on nutritional epidemiology. Statistics in Medicine 2014:33(12):2137–2155. [doi:10.1002/sim.6095](https://doi.org/10.1002/sim.6095)

- Nab L, Groenwold RHH, Welsing PMJ, van Smeden M. Measurement error 
in continuous endpoints in randomised trials: Problems and solutions. Statistics
in Medicine. 2019:38(27):5182-5196. [doi:10.1002/sim.8359](https://doi.org/10.1002/sim.8359)

- Rosner B, Spiegelman D, Willett WC. Correction of logistic regression relative risk estimates and confidence intervals for measurement error: The case of multiple covariates measured with error. 1990:132(4):734-745. [doi:10.1093/oxfordjournals.aje.a115715](https://doi.org/10.1093/oxfordjournals.aje.a115715)

- Rosner B, Spiegelman D, Willett WC. Correction of logistic regression relative risk estimates and confidence intervals for random within-person measurement error. American Journal of Epidemiology. 1992:136(11):1400-1413. [doi:10.1093/oxfordjournals.aje.a116453](https://doi.org/10.1093/oxfordjournals.aje.a116453)

- Spiegelman D, Carroll RJ, Kipnis V. Efficient regression calibration for logistic regression in main study/internal validation study designs with an imperfect reference instrument. Statistics in Medicine. 2001:20(1):139-160. [doi:10.1002/1097-0258(20010115)20:1<139::AID-SIM644>3.0.CO;2-K](https://doi.org/10.1002/1097-0258(20010115)20:1<139::AID-SIM644>3.0.CO;2-K)
---
title: "A user's guide to random measurement error correction in a covariate for R"
author: "Linda Nab"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{random_covme}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)
library(mecor)
```

## Introduction
__mecor__ is an R package for Measurement Error CORrection. __mecor__ implements measurement error correction methods for linear models with continuous outcomes. The measurement error can either occur in a continuous covariate or in the continuous outcome. This vignette discusses how a sensitivity analysis for random covariate measurement error is conducted in __mecor__. 

*Regression calibration* is one of the most popular measurement error correction methods for covariate measurement error. This vignette shows how *regression calibration* is used in __mecor__ to correct for random measurement error in a covariate. Our interest lies in estimating the association between a continuous reference exposure $X$ and a continuous outcome $Y$, given covariates $Z$. Instead of $X$, the substitute error-prone exposure $X^*$ is measured, assumed with random measurement error. It is further assumed that there is no extra information available to quantify the random measurement error in $X^*$. The input for our measurement error correction therefore is constrained to informed guesses about the size of the random measurement error. Literature or expert knowledge could be used to inform these guesses. We refer to the vignettes discussing e.g. *standard regression calibration* for random measurement error correction when validation data is available.

## Random measurement error
We assume that $X^*$ is measured with random measurement error. This means that we assume that $X^* = X + U$, where $U$ has mean 0 and variance $\tau^2$. More specifically, we assume non-differential random measurement error, i.e. $X^*|X$ is independent of $Y$ (our outcome). 

## Random measurement error correction in __mecor__
The object `MeasErrorRandom()` in __mecor__ is used for random measurement error correction in a covariate. We explain the usage of the `MeasErrorRandom()` object in the following. We first introduce the simulated data set `vat`. The simulated data set `vat` is an internal covariate-validation study. We will use this data set to explore random measurement error correction without using the reference measure visceral adipose tissue $VAT$ that is available in the data set. The data set `vat` contains 1000 observations of the outcome insulin resistance $IR_{ln}$, the error-prone exposure waist circumference $WC$ and the covariate age $age$. The reference exposure $VAT$ is observed in approximately 25% of the individuals in the study, but will be ignored. In this example, we assume that there is random measurement error in the exposure $WC$.
```{r load_data, eval = TRUE}
# load internal covariate validation study
data("vat", package = "mecor")
head(vat)
```
When ignoring the measurement error in $WC$, one would naively regress $WC$ and $age$ on $IR_{ln}$. This results in a biased estimation of the exposure-outcome association:
```{r uncorfit, eval = TRUE}
# naive estimate of the exposure-outcome association
data(vat)
lm(ir_ln ~ wc + age, data = vat)
```
Suppose that  $VAT$ is not observed in the internal covariate-validation study `vat`. To correct the bias in the naive association between exposure $WC$ and outcome $IR_{ln}$ given $age$, we need to make an informed guess about the quantity of $\tau^2$. Suppose we assume $\tau^2 = 0.25$. One can proceed as follows using `mecor()`:
```{r extrc, eval = TRUE}
# Use MeasErrorRandom for measurement error correction:
mecor(ir_ln ~ MeasErrorRandom(substitute = wc, variance = 0.25) + age,
      data = vat)
```

## How does __mecor__ do this?
To correct for the random measurement error in $WC$, __mecor__ constructs the calibration model matrix as follows:
```{r extrc2, eval = TRUE}
# First, construct the variance--covariance matrix of X_star and Z: 
# ( Var(X_star)   Cov(X_star, Z)
#   Cov(Z,X_star) Var(Z)       )
# To do so, we design Q, a matrix with 1000 rows (number of observations) and 2 
# columns. The first column of Q contains all 1000 observations of X_star, each 
# minus the mean of X_star. The second column of Q contains all 1000 obervations 
# of Z, each minus the mean of Z. 
Q <- scale(cbind(vat$wc, vat$age), scale = F)
# Subsequently, the variance--covariance matrix of X_star and Z is constructed:
matrix <- t(Q) %*% Q / (length(vat$age) - 1)
# Then, the variance--covariance matrix of X and Z is constructed, by using:
# Var(X) = Var(X_star) - Var(U) <--- Var(U) is the assumed tau^2
# Cov(X, Z) = Cov(X_star, Z)    <--- since U is assumed independent of Z
matrix1 <- matrix
matrix1[1, 1] <- matrix1[1, 1] - 0.25 # tau^2 = 0.25
# Rosner et al. (1992) show that the calibration model matrix can be constructed
# by taking the inverse of the variance--covariance matrix of X and Z and by
# matrix multiplying that matrix with the variance--covariance matrix of X_star
# and Z. 
model_matrix <- solve(matrix1) %*% matrix
model_matrix
matrix1 %*% solve(matrix)
# The resulting matrix is now:
# (1/lambda1        0
#  -lambda2/lambda1 1)
# Where,
# lambda1 = Cov(X,X_star|Z) / Var(X_star|Z)
# lambda2 = Cov(X,Z|X_star) / Var(Z|X_star) 
# Or, more familiar, the calibration model,
# E[X|X_star, Z] = lambda0 + lambda1 * X_star + lambda2 * Z
lambda1 <- 1 / model_matrix[1, 1]
lambda2 <- model_matrix[2,1] * - lambda1
# From standard theory, we have,
# lambda0 = mean(X) - lambda1 * mean(X_star) - lambda2 * mean(Z)
# mean(X) = mean(X_star) since we assume random measurement error
lambda0 <- mean(vat$wc) - lambda1 * mean(vat$wc) - lambda2 * mean(vat$age)
# The calibration model matrix Lambda is defined as:
# (lambda1 lambda0 lambda2
#  0       1       0
#  0       0       1)
model_matrix <- diag(3)
model_matrix[1, 1:3] <- c(lambda1, lambda0, lambda2)
model_matrix
# The calibration model matrix is standard output of mecor, and can be found
# using:
mecor_fit <- mecor(ir_ln ~ MeasErrorRandom(wc, 0.25) + age,
                   data = vat)
mecor_fit$corfit$matrix
```
Subsequently, the naive estimates of the outcome model are multiplied by the inverse of the calibration model matrix to obtain corrected estimates of the outcome model. 
```{r extrc3, eval = TRUE}
# Fit naive outcome model
naive_fit <- lm(ir_ln ~ wc + age, 
                data = vat)
# Save coefficients
beta_star <- naive_fit$coefficients
# To prepare the coefficients for the measurement error correction, exchange the
# intercept and the coefficient for X_star
beta_star[1:2] <- rev(beta_star[1:2]) 
# Perform the measurement error correction:
beta <- beta_star %*% solve(model_matrix)
# Reverse the order 
beta[1:2] <- rev(beta[1:2])
beta # corrected coefficients of the outcome model
```
Which exactly matches the output of `mecor()` above.



---
title: "Measurement error correction in a continuous trial endpoint"
author: "Linda Nab"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{mecor_mece}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)
library(mecor)
```
This vignette shows how one can correct for the bias in the trial's effect estimate using the R package mecor. We make use of external validation data. Suppose an endpoint in a trial is measured with error, i.e., the substitute endpoint $Y^*$ instead of the reference endpoint $Y$ is observed. 

First, we simulate some example data for a trial composed out of two groups. For example, a placebo group ($X = 0$) and an active comparator ($X = 1$). The number of individuals included in the trial is set to 1000, 500 individuals in each group. Suppose the substitute endpoint $Y^*$ is observed instead of $Y$. Further, suppose that an external validation set of sample size 500 is available in which both $Y^*$ and $Y$ are measured. 
```{r example1}
# simulate the trial's data
X <- rep(c(0,1), 500)
Y <- 2 * X + rnorm(1000, 0, 1) # estimand: 2
# introduce measurement error
Y_star <- 1.1 * Y + rnorm(1000, 0, 1)
trial <- cbind.data.frame(X = X, Y_star = Y_star)
# simulate an external validation data set
Y <- rnorm(100, 2, 1)
Y_star <- 1.1 * Y + rnorm(500, 0, 1)
trial_ext <- cbind.data.frame(Y = Y, Y_star = Y_star)
```

When the error is ignored, one would estimate the trial's effect by regressing $X$ on $Y^*$.
```{r example2}
# uncorrected estimate of the trial's effect:
uncor_fit <- lm(Y_star ~ X, data = trial)
uncor_fit$coefficients
```
As you might expect, the trial's effect estimate does not equal 2, to which value the estimand was set when generating the data. To obtain an unbiased trial effect, measurement error correction is needed. First, we estimate the parameters of the measurement error model using our external validation data:
```{r example3}
memod_fit <- lm(Y_star ~ Y, data = trial_ext)
memod_fit$coefficients
```

Then, mecor can be used to correct for the measurement error in the trial's effect estimate as follows:
```{r example4}
cor_fit <- mecor(MeasErrorExt(substitute = Y_star, model = memod_fit) ~ X,
                 data = trial,
                 method = "standard",
                 B = 0 # for bootstrap intervals, set to e.g. 999
                 )
```

Confidence intervals for the corrected estimate can be obtained by using the summary object:
```{r example5}
summary(cor_fit, fieller = TRUE, zerovar = TRUE)
```

When there is no external validation data available. One could conduct a sensitivity analysis by making informed guesses about the parameters values of the measurement error model. Suppose e.g. we guess the following measurement error model: $Y^* = 1.1 Y$. The following code can be used to quantify the impact of the measurement error would be on the trial's effect estimate: 
```{r example6}
sens_fit <- mecor(MeasErrorExt(substitute = Y_star, 
                               model = list(coef = c(0, 1.1))) ~ X, 
                  data = trial,
                  method = "standard"
                  )
sens_fit
```
---
title: "A user's guide to standard covariate measurement error correction for R"
author: "Linda Nab"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{standard_covme}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)
library(mecor)
```
## Introduction
__mecor__ is an R package for Measurement Error CORrection. __mecor__ implements measurement error correction methods for linear models with continuous outcomes. The measurement error can either occur in a continuous covariate or in the continuous outcome. This vignette discusses covariate measurement error correction by means of *standard regression calibration* in __mecor__. 

*Regression calibration* is one of the most popular measurement error correction methods for covariate measurement error. This vignette shows how *standard regression calibration* is applied in an internal validation study, a replicates study, a calibration study and an external validation study. Each of the four studies will be introduced in the subsequent sections, along with examples of how *standard regression calibration* can be applied in each of the four studies using __mecor__'s function `mecor()`. In all four studies, our interest lies in estimating the association between a continuous reference exposure $X$ and a continuous outcome $Y$, given covariates $Z$. Instead of $X$, the substitute error-prone exposure $X^*$ is measured. 

## Internal validation study
The simulated data set `vat` in __mecor__ is an internal covariate-validation study. An internal validation study includes a subset of individuals of whom the reference exposure $X$ is observed. The data set `vat` contains 1000 observations of the outcome insulin resistance $IR_{ln}$, the error-prone exposure waist circumference $WC$ and the covariates sex, age and total body fat, $sex$, $age$ and $TBF$, respectively. The reference exposure visceral adipose tissue $VAT$ is observed in approximately 25% of the individuals in the study.
```{r load_data, eval = TRUE}
# load internal covariate validation study
data("vat", package = "mecor")
head(vat)
```
When ignoring the measurement error in $WC$, one would naively regress $WC$ and $sex$, $age$ and $TBF$ on $IR_{ln}$. This results in a biased estimation of the exposure-outcome association:
```{r uncorfit, eval = TRUE}
# naive estimate of the exposure-outcome association
lm(ir_ln ~ wc + sex + age + tbf, data = vat)
```
Alternatively, one could perform an analysis restricted to the internal validation set: 
```{r ivr, eval = TRUE}
# analysis restricted to the internal validation set
lm(ir_ln ~ vat + sex + age + tbf, data = subset(vat, !is.na(vat)))
```
Although the above would result in an unbiased estimation of the exposure-outcome association, approximately 75% of the data is thrown out. Instead of doing an analysis restricted to the internal validation set, you could use *standard regression calibration* to correct for the measurement error in $X^*$. The following code chunk shows *standard regression calibration* with `mecor()`:
```{r rc, eval = TRUE}
mecor(formula = ir_ln ~ MeasError(substitute = wc, reference = vat) + sex + age + tbf,
      data = vat,
      method = "standard", # defaults to "standard"
      B = 0) # defaults to 0
```
As shown in the above code chunk, the `mecor()` function needs a *formula* argument, a *data* argument, a *method* argument and a *B* argument. Presumably, you are familiar with the structure of a formula in R. The only thing that's different here is the use of a `MeasError()` object in the formula. A `MeasError()` object is used to declare the substitute measure, in our case $WC$, and the reference measure, in our case $VAT$. The *B* argument of `mecor()` is used to calculate bootstrap confidence intervals for the corrected coefficients of the model. Let us construct 95% confidence intervals using the bootstrap with 999 replicates:
```{r rc_se, eval = TRUE, results = 'hide'}
# save corrected fit
rc_fit <- 
  mecor(formula = ir_ln ~ MeasError(substitute = wc, reference = vat) + age + sex + tbf,
        data = vat,
        method = "standard", # defaults to "standard"
        B = 999) # defaults to 0
```
Print the confidence intervals to the console using `summary()`:
```{r rc_se_sum, eval = TRUE}
summary(rc_fit)
```
Two types of 95% confidence intervals are shown in the output of the `summary()` object. Bootstrap confidence intervals and Delta method confidence intervals. The default method to constructing confidence intervals in __mecor__ is the Delta method. Further, Fieller method confidence intervals and zero variance method confidence intervals can be constructed with `summary()`:
```{r rc_se2, eval = TRUE}
# fieller method ci and zero variance method ci and se's for 'rc_fit'
summary(rc_fit, zerovar = TRUE, fieller = TRUE)
```
Fieller method confidence intervals are only constructed for the corrected covariate (in this case $X$). 

## Replicates study
The simulated data set `bloodpressure` in __mecor__ is a replicates study. A replicates study includes a subset of individuals of whom the error-prone substitute exposure is repeatedly measured. The dataset `bloodpressure` contains 1000 observations of the outcome $creatinine$, three replicate measures of the error-prone exposure systolic blood pressure $sbp30, sbp60$ and $sbp120$, and one covariates age $age$. It is assumed that there is 'random' measurement error in the repeatedly measured substitute exposure measure. 
```{r load_data2, eval = TRUE}
# load replicates study
data("bloodpressure", package = "mecor")
head(bloodpressure)
```
When ignoring the measurement error in $sbp30$, one would naively regress $sbp30$, $age$ on $creatinine$. Which results in a biased estimation of the exposure-outcome association:
```{r uncorfit2, eval = TRUE}
# naive estimate of the exposure-outcome association
lm(creatinine ~ sbp30 + age, 
   data = bloodpressure)
```
Or alternatively, one could calculate the mean of each of the three replicate measures. Yet, this would still lead to a biased estimation of the exposure-outcome association:
```{r uncorfit3, eval = TRUE}
## calculate the mean of the three replicate measures
bloodpressure$sbp_123 <- with(bloodpressure, rowMeans(cbind(sbp30, sbp60, sbp120)))
# naive estimate of the exposure-outcome association version 2
lm(creatinine ~ sbp_123 + age,
   data = bloodpressure)
```
For an unbiased estimation of the exposure-outcome association, one could use regression calibration using `mecor()`:
```{r rc2, eval = TRUE}
mecor(formula = creatinine ~ MeasError(sbp30, replicate = cbind(sbp60, sbp120)) + age,
      data = bloodpressure)
```
Instead of using the *reference* argument in the `MeasError()` object, the *replicate* argument is used. Standard errors of the regression calibration estimator and confidence intervals can be constructed similar to what was shown for an internal validation study.

## Calibration study
The simulated data set `sodium` in __mecor__ is a outcome calibration study. In a calibration study, two types of measurements are used to measure the outcome (or exposure). A measurement method prone to 'systematic' error, and a measurement method prone to 'random' error. The measurement prone to 'systematic' error is observed in the full study, the measurement prone to 'classical' error is observed in a subset of the study and repeatedly measured. The dataset `sodium` contains 1000 observations of the systematically error prone outcome $recall$, the randomly error prone outcome $urinary1$ and $urinary2$, and the exposure (in our case a indicator for diet) $diet$. The two replicate measures of the outcome prone to random error are observed in 498 individuals (approximately 50 percent).
```{r load_data3, eval = TRUE}
# load calibration study
data("sodium", package = "mecor")
head(sodium)
```
When ignoring the measurement error in $recall$, one would naively regress $diet$ on $recall$. Which results in a biased estimation of the exposure-outcome association:
```{r uncorfit4, eval = TRUE}
## uncorrected regression
lm(recall ~ diet, 
   data = sodium)
```
Alternatively, one could use the first half of the study population and use the mean of each of the two replicate measures. This would lead to an unbiased estimation of the exposure-outcome association since there is random measurement error int he replicate measures:
```{r uncorfit5, eval = TRUE}
## calculate mean of three replicate measures
sodium$urinary_12 <- with(sodium, rowMeans(cbind(urinary1, urinary2)))
## uncorrected regression version 2
lm(urinary_12 ~ diet,
   data = sodium)
```
For an unbiased estimation of the exposure-outcome association, one could alternatively use *standard regression calibration* using `mecor()`:
```{r rc3, eval = TRUE}
mecor(formula = MeasError(substitute = recall, replicate = cbind(urinary1, urinary2)) ~ diet,
      data = sodium)
```
Standard errors of the regression calibration estimator and confidence intervals can be constructed similar to what was shown for an internal validation study.

## External validation study
The simulated data set `heamoglogin_ext` in __mecor__ is a external outcome-validation study. An external validation study is used when in the main study, no information is available to correct for the measurement error in the outcome $Y$ (or exposure). Suppose for example that venous heamoglobin levels $venous$ are not observed in the internal outcome-validation study `haemoglobin`. An external validation study is a (small) sub study containing observations of the reference measure venous heamoglobin levels $venous$, the error-prone substitute measure $capillary$ (and the covariate(s) $Z$ in case of an covariate-validation study) of the original study. The external validation study is then used to estimate the calibration model, that is subsequently used to correct for the measurement error in the main study. 
```{r load_data4, eval = TRUE}
# load internal covariate validation study
data("haemoglobin_ext", package = "mecor")
head(haemoglobin_ext)
data("haemoglobin", package = "mecor")
```
Suppose reference measure $X$ is not observed in dataset `icvs`. To correct the bias in the naive association between exposure $X^*$ and outcome $Y$ given $Z$, using the external validation study, one can proceed as follows using `mecor()`:
```{r extrc, eval = TRUE}
# Estimate the calibration model in the external validation study
calmod <- lm(capillary ~ venous, 
             data = haemoglobin)
# Use the calibration model for measurement error correction:
mecor(MeasErrorExt(substitute = capillary, model = calmod) ~ supplement,
      data = haemoglobin)
```
In the above, a `MeasErrorExt()` object is used, indicating that external information is used for measurement error correction. The model argument of a `MeasErrorExt()` object takes a linear model of class lm (in the above `calmod`). Alternatively, a named list with the coefficients of the calibration model can be used as follows:
```{r extrc2, eval = TRUE}
# Use coefficients for measurement error correction:
mecor(MeasErrorExt(capillary, model = list(coef = c(-7, 1.1))) ~ supplement,
      data = haemoglobin)
```



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mecor.R
\name{mecor}
\alias{mecor}
\title{mecor: a Measurement Error Correction Package}
\usage{
mecor(formula, data, method = "standard", B = 0)
}
\arguments{
\item{formula}{an object of class \link[stats]{formula} (or one that is
coerced to that class): a symbolic description of the regression model
containing a \link[mecor]{MeasError}, \link[mecor]{MeasErrorExt} or
\link[mecor]{MeasErrorRandom} object in one of the covariates or the outcome.}

\item{data}{a data.frame, list or environment (or object coercible by
as.data.frame to a data frame) containing the variables in the model
specified in \code{formula}.}

\item{method}{a character string indicating the method used to correct for
the measurement error, either "standard" (regression calibration for
covariate measurement error and method of moments for outcome measurement
error), "efficient" (efficient regression calibration for covariate
measurement error and efficient method of moments for outcome measurement
error), "valregcal" (validation regression calibration) or "mle" (maximum
likelihood estimation). Defaults to "standard".}

\item{B}{number of bootstrap samples, defaults to 0.}
}
\value{
\code{mecor} returns an object of \link[base]{class} "mecor".

An object of class \code{mecor} is a list containing the following components:

\item{corfit}{a list containing the corrected fit, including the coefficients
of the corrected fit (\code{coef}) and the variance--covariance matrix of the
coefficients of the corrected fit obtained by the delta method (\code{vcov}),
and more depending on the method used.}
\item{uncorfit}{an \link[stats]{lm.fit} object of the uncorrected fit.}
}
\description{
mecor provides correction methods for measurement error in a continuous
covariate or outcome in linear regression models with a continuous outcome
}
\examples{
## measurement error in a covariate/outcome:
# internal covariate-validation study
data(vat)
out <-
mecor(ir_ln ~ MeasError(wc, reference = vat) + sex + age + tbf,
      data = vat,
      method = "standard",
      B = 999)
# replicates study
data(bloodpressure)
mecor(creatinine ~ MeasError(sbp30, replicate = cbind(sbp60, sbp120)) + age,
      data = bloodpressure,
      method = "mle")
# outcome-calibration study
data(sodium)
mecor(MeasError(recall, replicate = cbind(urinary1, urinary2)) ~ diet,
      data = sodium,
      method = "efficient")
# external outcome-validation study
data(haemoglobin_ext)
calmod_fit <- lm(capillary ~ venous, data = haemoglobin_ext)
data(haemoglobin) # suppose reference venous is not available
mecor(MeasErrorExt(capillary, model = calmod_fit) ~ supplement,
      data = haemoglobin)
# sensitivity analyses
data(vat) # suppose reference vat is not available
# guesstimate the coefficients of the calibration model:
mecor(ir_ln ~ MeasErrorExt(wc, model = list(coef = c(0.2, 0.5, -1.3, 0, 0.6))) + sex + age + tbf,
      data = vat)
# assume random measurement error in wc of magnitude 0.25:
mecor(ir_ln ~ MeasErrorRandom(wc, variance = 0.25) + sex + age + tbf,
      data = vat)
data(bloodpressure) # suppose replicates sbp60 and sbp60 are not available
mecor(creatinine ~ MeasErrorRandom(sbp30, variance = 25) + age,
      data = bloodpressure)

## differential measurement error in the outcome:
# internal outcome-validation study
mecor(MeasError(capillary, reference = venous, differential = supplement) ~ supplement,
      data = haemoglobin,
      method = "standard")
}
\references{
L. Nab, R.H.H. Groenwold, P.M.J. Welsing, and  M. van Smeden.
Measurement error in continuous endpoints in randomised trials: problems and
solutions

L. Nab, M. van Smeden, R.H. Keogh, and R.H.H. Groenwold.
mecor: an R package for measurement error correction in linear models with
continuous outcomes
}
\author{
Linda Nab, \email{l.nab@lumc.nl}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeasError.R
\name{MeasError}
\alias{MeasError}
\title{Create a Measurement Error Object}
\usage{
MeasError(substitute, reference, replicate, differential)
}
\arguments{
\item{substitute}{a vector containing the error-prone measure}

\item{reference}{a vector containing the reference measure assumed without
measurement error}

\item{replicate}{a vector or matrix with replicates of the error-prone
measure with classical measurement error. This can either be
replicates obtained by using the same measurement method as the substitute
measure (replicates study) or replicates using a different measurement method
than the substitute measure (calibration study).}

\item{differential}{a vector containing the variable to which the measurement
error is differential.}
}
\value{
\code{MeasError} returns an object of \link[base]{class} "MeasError".

An object of class \code{MeasError} is a list containing the substitute and
reference (and replicate or differential if applicable) variables and has
attributes input (the name of the substitute and reference or replicate
and differential (if applicable) variables) and call (the matched call).
}
\description{
This function creates a measurement error object, usually used as a covariate
or the outcome in the \code{formula} argument of \link[mecor]{mecor} if one
wants to correct for the measurement error in that variable using a reference
variable or a replicate measure.
}
\examples{
## measurement error in a covariate:
# internal covariate-validation study
data(vat)
with (vat, MeasError(substitute = wc,
                     reference = vat))
# replicates study
data(bloodpressure)
with (bloodpressure, MeasError(substitute = sbp30,
                               replicate = cbind(sbp60, sbp120)))
# outcome-calibration study
data(sodium)
with(sodium, MeasError(substitute = recall,
                       replicate = cbind(urinary1, urinary2)))
## measurement error in the outcome:
# internal outcome-validation study
data(haemoglobin)
with(haemoglobin, MeasError(substitute = capillary,
                            reference = venous))
# internal outcome- validation study with differential measurement error in
# the dependent variable
data(haemoglobin)
with(haemoglobin, MeasError(substitute = capillary,
                            reference = venous,
                            differential = supplement))
}
\author{
Linda Nab, \email{l.nab@lumc.nl}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bloodpressure}
\alias{bloodpressure}
\title{PDAC blood pressure data [replicates study]}
\format{
A data frame with 450 rows and 6 variables:
\describe{
  \item{creatinine}{Serum creatinine (umol/L)}
  \item{age}{Age (years)}
  \item{sbp30}{Systolic blood pressure at 30 minutes (mm Hg)}
  \item{sbp60}{Systolic blood pressure at 60 minutes (mm Hg)}
  \item{sbp90}{Systolic blood pressure at 90 minutes (mm Hg)}
  \item{sbp120}{Systolic blood pressure at 120 minutes (mm Hg)}
}
}
\usage{
bloodpressure
}
\description{
Blood pressure, age and creatinine levels of 450 pregnant women from the
Pregnancy Day Assessment Clinic.
}
\details{
This is a simulated dataset inspired by data that was originally published at the Dryad Digital Repository: <doi:10.5061/dryad.0bq15>
}
\examples{
data("bloodpressure", package = "mecor")
}
\references{
Elizabeth Anne McCarthy, Thomas A Carins, Yolanda Hannigan, Nadia Bardien, Alexis Shub, and Susan P Walker. Data from: Effectiveness and safety of 1 vs 4h blood pressure profile with clinical and laboratory assessment for the exclusion of gestational hypertension and pre-eclampsia: a retrospective study in a university affiliated maternity hospital. Dryad (2015). <doi:10.5061/dryad.0bq15>.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vat_ext}
\alias{vat_ext}
\title{Visceral adipose tissue external data [external covariate-validation study]}
\format{
A data frame with 100 rows and 5 variables:
\describe{
  \item{wc}{Waist circumference (standardised, cm)}
  \item{vat}{Visceral adipose tissue (standardised, cm^2)}
  \item{sex}{Sex (0 = male, 1 = female)}
  \item{age}{Age (years)}
  \item{tbf}{Total body fat (standardised, \%)}
}
}
\usage{
vat_ext
}
\description{
Waist circumference, visceral adipose tissue, sex, age, and total body fat of 100 individuals
}
\details{
This is a simulated data set accompanying the dataset "vat", that is inspired by the NEO data <doi:10.1007/s10654-013-9801-3>. A motivating example using the example data can be found here: <doi:10.1093/aje/kwab114>
}
\examples{
data("vat_ext", package = "mecor")
}
\references{
Renee de Mutsert, Martin den Heijer, Ton J Rabelink, Johannes WA Smit, Johannes A Romijn, Johan W Jukema, Albert de Roos, Christa M Cobbaert, Margreet Kloppenburg, Saskia le Cessie, Saskia Middeldorp, Frits R Rosendaal. The Netherlands epidemiology of obesity (NEO) study: Study design and data collection. European Journal of Epidemiology (2013). <doi:10.1007/s10654-013-9801-3>

Linda Nab, Maarten van Smeden, Renee de Mutsert, Frits R Rosendaal, and Rolf HH Groenwold. Sampling strategies for internal validation samples for exposure measurement error correction: A study of visceral adipose tissue measures replaced by waist circumference measures. American Journal of Epidemiology (2021). <doi:10.1093/aje/kwab114>
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3-mecor.R
\name{summary.mecor}
\alias{summary.mecor}
\title{Summarizing Measurement Error Correction}
\usage{
\method{summary}{mecor}(object, alpha = 0.05, zerovar = FALSE, fieller = FALSE, ...)
}
\arguments{
\item{object}{an object of class "mecor", a result of a call to
\link[mecor]{mecor}.}

\item{alpha}{probability of obtaining a type II error.}

\item{zerovar}{a boolean indicating whether standard errors and confidence
intervals using the zerovariance method must be added to the summary object.}

\item{fieller}{a boolean indicating whether confidence intervals using the
fieller method must be added to the summary object.}

\item{...}{additional arguments affecting the summary produced}
}
\value{
The function \code{summary.mecor} returns a list of summary statistics of the
fitted corrected model and fitted uncorrected model.

\item{call}{the matched call}
\item{c}{summary of the corrected fit}
\item{uc}{summary of the uncorrected fit}
\item{B}{number of bootstrap replicates used}
\item{alpha}{alpha level used}
}
\description{
\code{summary} method for class "mecor"
}
\examples{
## measurement error in a covariate:
# internal covariate-validation study
data(vat)
mecor_fit <- mecor(ir_ln ~ MeasError(wc, reference = vat) + sex + age + tbf,
                   data = vat,
                   method = "standard")
summary(mecor_fit)
summary(mecor_fit, zerovar = TRUE, fieller = TRUE)
summary(mecor_fit, alpha = 0.10)

}
\seealso{
The model fitting function \link[mecor]{mecor}, \link[base]{summary}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sim}
\alias{sim}
\title{Simulated dataset for the \link[mecor]{ipwm} function}
\format{
A data frame with 5000 rows and 14 variables:
\describe{
  \item{L1}{covariate, binary}
  \item{L2}{covariate, continuous}
  \item{L3}{covariate, binary}
  \item{L4}{covariate, continuous}
  \item{L5}{covariate, binary}
  \item{L6}{covariate, binary}
  \item{L7}{covariate, continuous}
  \item{L8}{covariate, binary}
  \item{L9}{covariate, binary}
  \item{L10}{covariate, continuous}
  \item{A}{exposure, binary}
  \item{Y}{outcome, binary}
  \item{B}{misclassified exposure, binary}
  \item{Z}{misclassified outcome, binary}
}
}
\usage{
sim
}
\description{
A simulated dataset containing 5000 observations of the covariates L1-L10,
the true exposure A and true outcome Y, and the misclassified exposure B and
misclassified outcome Z.
}
\examples{
data("sim", package = "mecor")
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ipwm.R
\name{ipwm}
\alias{ipwm}
\title{Weighting for Confounding and Joint Misclassification of Exposure and Outcome}
\usage{
ipwm(
  formulas,
  data,
  outcome_true,
  outcome_mis = NULL,
  exposure_true,
  exposure_mis = NULL,
  nboot = 1000,
  conf_level = 0.95,
  fix_nNAs = FALSE,
  semiparametric = FALSE,
  optim_args = list(method = "BFGS"),
  force_optim = FALSE,
  sp = Inf,
  print = TRUE
)
}
\arguments{
\item{formulas}{a list of objects of \link[base]{class} \code{\link[stats]{formula}} specifying the probability models for the stats::terms of some factorisation of the joint conditional probability function of \code{exposure_true}, \code{exposure_mis}, \code{outcome_true} and \code{outcome_mis}, given covariates}

\item{data}{\code{\link[base]{data.frame}} containing \code{exposure.true}, \code{exposure.mis}, \code{outcome.true}, \code{outcome.mis} and covariates. Missings (\code{NA}s) are allowed on variables \code{exposure_true} and \code{outcome_true}.}

\item{outcome_true}{a character string specifying the name of the true outcome variable that is free of misclassification but possibly unknown (\code{NA}) for some (but not all) subjects}

\item{outcome_mis}{a character string specifying the name of the counterpart of \code{outcome_true} that is available on all subjects but potentially misclassifies subjects' outcomes. The default (\code{outcome_mis = NULL}) indicates absence of outcome misclassification}

\item{exposure_true}{a character string specifying the name of the true exposure variable that is free of misclassification but possibly unknown (\code{NA}) for some (but not all) subjects}

\item{exposure_mis}{a character string specifying the name of the counterpart of \code{exposure_true} that is available on all subjects but potentially misclassifies subjects as exposed or as non-exposed. The default (\code{exposure_mis = NULL}) indicates absence of exposure misclassification}

\item{nboot}{number of bootstrap samples. Setting \code{nboot == 0} results in point estimation only.}

\item{conf_level}{the desired confidence level of the confidence interval}

\item{fix_nNAs}{logical indicator specifying whether or not to fix the joint distribution of \code{is.na(exposure_true)} and \code{is.na(outcome_true)}. If \code{TRUE}, stratified bootstrap sampling is done according to the missing data pattern.}

\item{semiparametric}{logical indicator specifying whether or not to parametrically sample \code{exposure_true}, \code{exposure_mis}, \code{outcome_true} and \code{outcome_mis}. If \code{semiparametric == TRUE}, it is assumed that the missing data pattern is conditionally independent of these variables given covariates. Provided \code{nboot > 0}, the missing data pattern and covariates are sampled nonparametrically. \code{semiparametric} is ignored if \code{nboot == 0}.}

\item{optim_args}{arguments passed onto \code{\link[stats]{optim}} if called. See Details below for more information.}

\item{force_optim}{logical indicator specifying whether or not to force the \code{\link[stats]{optim}} function to be called}

\item{sp}{scalar shrinkage parameter in the interval \code{(0, Inf)}. Values closer to zero result in greater shrinkage of the estimated odds ratio to unity; \code{sp == Inf} results in no shrinkage.}

\item{print}{logical indicator specifying whether or not to print the output.}
}
\value{
\code{ipwm} returns an object of \link[base]{class} \code{ipwm}.
The returned object is a list containing the following elements:

\item{logOR}{the estimated log odds ratio;}
\item{call}{the matched function call.}

If \code{nboot != 0}, the list also contains

\item{SE}{a bootstrap estimate of the standard error for the estimator of the log odds ratio;}
\item{CI}{a bootstrap percentile confidence interval for the log odds ratio.}
}
\description{
\code{ipwm} implements a method for estimating the marginal causal odds ratio by constructing weights (modified inverse probability weights) that address both confounding and joint misclassification of exposure and outcome.
}
\details{
This function is an implementation of the weighting method described by Penning de Vries et al. (2018).
The method defaults to the estimator proposed by Gravel and Platt (2018) in the absence of exposure misclassification.

The function assumes that the exposure or the outcome has a misclassified version. An error is issued when both \code{outcome_mis} and \code{exposure_mis} are set to \code{NULL}.

Provided \code{force_optim = FALSE}, \code{ipwm} is considerably more efficient when the \code{\link[stats]{optim}} function is not invoked; i.e., when (1) \code{exposure_mis = NULL} and the formula for \code{outcome_true} does not contain stats::terms involving \code{outcome_mis} or \code{exposure_true}, (2) \code{outcome_mis = NULL} and the formula for \code{exposure_true} does not contain stats::terms involving \code{exposure_mis} or \code{outcome_true}, or (3) \code{all(is.na(data[, exposure_true]) == is.na(data[, outcome_true]))} and the formulas for \code{exposure_true} and \code{outcome_true} do not contain stats::terms involving \code{exposure_mis} or \code{outcome_mis}. In these cases, \code{ipwm} uses iteratively reweighted least squares via the \code{\link[stats]{glm}} function for maximum likelihood estimation. In all other cases, \code{optim_args} is passed on to \code{\link[stats]{optim}} for optimisation of the joint likelihood of \code{outcome_true}, \code{outcome_mis}, \code{exposure_true} and \code{exposure_mis}.
}
\examples{
data(sim) # simulated data on 10 covariates, exposure A and outcome Y.
formulas <- list(
  Y ~ A + L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10 + B + Z,
  A ~ L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10 + B + Z,
  Z ~ L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10 + B,
  B ~ L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10
)
\dontrun{
ipwm_out <- ipwm(
  formulas = formulas,
  data = sim,
  outcome_true = "Y",
  outcome_mis = "Z",
  exposure_true = "A",
  exposure_mis = "B",
  nboot = 200,
  sp = 1e6
)
ipwm_out
}

}
\references{
Gravel, C. A., & Platt, R. W. (2018). Weighted estimation for confounded binary outcomes subject to misclassification. \emph{Statistics in medicine}, 37(3), 425-436. https://doi.org/10.1002/sim.7522

Penning de Vries, B. B. L., van Smeden, M., & Groenwold, R. H. H. (2020). A weighting method for simultaneous adjustment for confounding and joint exposure-outcome misclassifications. \emph{Statistical Methods in Medical Research}, 0(0), 1-15. https://doi.org/10.1177/0962280220960172
}
\author{
Bas B. L. Penning de Vries, \email{b.b.l.penning_de_vries@lumc.nl}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{haemoglobin}
\alias{haemoglobin}
\title{Low-dose iron supplements haemoglobin data [internal outcome-validation study]}
\format{
A data frame with 400 rows and 3 variables:
\describe{
  \item{capillary}{Haemoglobin levels measured in capillary blood (g/L)}
  \item{supplement}{Low-dose iron supplement (20 mg/d) (0 = no, 1 = yes)}
  \item{venous}{Haemoglobin levels measured in venous blood (g/L)}
}
}
\usage{
haemoglobin
}
\description{
Capillary haemoglobin and venous haemoglobin levels of 400 subjects of a trial investigating the efficacy of low-dose iron supplements during pregnancy.
Venous haemoglobin levels were observed of approximately 25\% of the subjects included in the trial.
}
\details{
This is a simulated data set inspired by a trial investigating low-dose iron supplements <doi:10.1093/ajcn/78.1.145>. A motivating example using the example data can be found here: <doi:10.1002/sim.8359>
}
\examples{
data("haemoglobin", package = "mecor")
}
\references{
Maria Makrides, Caroline A Crowther, Robert A Gibson, Rosalind S Gibson, and C Murray Skeaff. Efficacy and tolerability of low-dose iron supplements during pregnancy: a randomized controlled trial. The American Journal of Clinical Nutrition (2003). <doi:10.1093/ajcn/78.1.145>

Linda Nab, Rolf HH Groenwold, Paco MJ Welsing, and Maarten van Smeden. Measurement error in continuous endpoints in randomised trials: Problems and solutions. Statistics in Medicine (2019). <doi:10.1002/sim.8359>
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeasErrorRandom.R
\name{MeasErrorRandom}
\alias{MeasErrorRandom}
\title{Create a Random Measurement Error Object}
\usage{
MeasErrorRandom(substitute, variance)
}
\arguments{
\item{substitute}{a vector containing the error-prone measure}

\item{variance}{a numeric quantifying the assumed variance of the random measurement error}
}
\value{
\code{MeasErrorRandom} returns an object of \link[base]{class}
"MeasErrorRandom".

An object of class \code{MeasErrorRandom} is a list containing the substitute
variable, the assumed variance of the random measurement error in that variable and, the
attributes input (the name of the substitute variable) and call (the matched
call).
}
\description{
This function creates a random measurement error object, usually used as
a covariate in the \code{formula} argument of \link[mecor]{mecor} if one
wants to correct for random measurement error in that variable
}
\examples{
## random measurement error in a covariate:
# internal covariate-validation study
data(bloodpressure)
with(bloodpressure, MeasErrorRandom(sbp30, variance = 0.25))
}
\author{
Linda Nab, \email{l.nab@lumc.nl}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeasErrorExt.R
\name{MeasErrorExt}
\alias{MeasErrorExt}
\title{Create an External Measurement Error Object}
\usage{
MeasErrorExt(substitute, model)
}
\arguments{
\item{substitute}{a vector containing the error-prone measure}

\item{model}{a fitted linear model of class \link[stats]{lm} or a named
\link[base]{list}. The \link[base]{list} contains a vector named \code{coef}:
the coefficients of the calibration model or measurement error model and an
optional matrix named \code{vcov}: the variance--covariance matrix of the
coefficients}
}
\value{
\code{MeasErrorExt} returns an object of \link[base]{class}
"MeasErrorExt".

An object of class \code{MeasErrorExt} is a list containing the substitute
variable and the fitted calibration model or measurement error model and has
attributes input (the name of the substitute variable) and call (the matched
call).
}
\description{
This function creates an external measurement error object, usually used as
a covariate or the outcome in the \code{formula} argument of
\link[mecor]{mecor} if one wants to correct for the measurement error in that
variable using external data or externally estimated coefficients of the
calibration model (covariate-measurement error) or measurement error model
(outcome-measurement error)
}
\examples{
## measurement error in a outcome:
# external outcome-validation study
data(haemoglobin_ext)
# calibration model
calmod_fit <- lm(capillary ~ venous, data = haemoglobin)
# the external covariate-validation study can be used to correct for the
# measurement error in X_star in the dataset 'icvs', using the fitted
# calibration model
data(haemoglobin)
with (haemoglobin, MeasErrorExt(substitute = capillary,
                                model = calmod_fit))
# identical to:
calmod_coef <- coefficients(calmod_fit)
calmod_vcov <- vcov(calmod_fit)
with (haemoglobin, MeasErrorExt(substitute = capillary,
                                model = list(coef = calmod_coef,
                                             vcov = calmod_vcov)))
# when no external data is available, guesstimations of the coefficients of
# the calibration model can be used instead:
with (haemoglobin, MeasErrorExt(substitute = capillary,
                                model = list(coef = c('(Intercept)' = -7,
                                                      'venous' = 1.1))))
}
\author{
Linda Nab, \email{l.nab@lumc.nl}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{haemoglobin_ext}
\alias{haemoglobin_ext}
\title{Haemoglobin external data [external outcome-validation study]}
\format{
A data frame with 100 rows and 2 variables:
\describe{
  \item{capillary}{Haemoglobin levels measured in capillary blood (g/L)}
  \item{venous}{Haemoglobin levels measured in venous blood (g/L)}
}
}
\usage{
haemoglobin_ext
}
\description{
Capillary haemoglobin and venous haemoglobin levels of 100 individuals.
}
\details{
This is a simulated data set accompanying the dataset "haemoglobin", that is inspired by a trial investigating low-dose iron supplements <doi:10.1093/ajcn/78.1.145>. A motivating example using the example data can be found here: <doi:10.1002/sim.8359>
}
\examples{
data("haemoglobin_ext", package = "mecor")
}
\references{
Maria Makrides, Caroline A Crowther, Robert A Gibson, Rosalind S Gibson, and C Murray Skeaff. Efficacy and tolerability of low-dose iron supplements during pregnancy: a randomized controlled trial. The American Journal of Clinical Nutrition (2003). <doi:10.1093/ajcn/78.1.145>

Linda Nab, Rolf HH Groenwold, Paco MJ Welsing, and Maarten van Smeden. Measurement error in continuous endpoints in randomised trials: Problems and solutions. Statistics in Medicine (2019). <doi:10.1002/sim.8359>
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vat}
\alias{vat}
\title{NEO visceral adipose tissue data [internal covariate-validation study]}
\format{
A data frame with 650 rows and 6 variables:
\describe{
  \item{ir_ln}{Natural logarithm of insulin resistance (fasting glucose (mmol/L) x fasting insulin (mU/L) / 22.5)}
  \item{wc}{Waist circumference (standardised, cm)}
  \item{sex}{Sex (0 = male, 1 = female)}
  \item{age}{Age (years)}
  \item{tbf}{Total body fat (standardised, \%)}
  \item{vat}{Visceral adipose tissue (standardised, cm^2)}
}
}
\usage{
vat
}
\description{
Insulin resistance, waist circumference, sex, age, total body fat and visceral adipose tissue of 650 individuals from the
Netherlands Epidemiology of Obesity (NEO) study. Visceral adipose tissue measurements were taken of approximately 40\% of the individuals, at random.
}
\details{
This is a simulated data set inspired by the NEO data <doi:10.1007/s10654-013-9801-3>. A motivating example using the example data can be found here: <doi:10.1093/aje/kwab114>
}
\examples{
data("vat", package = "mecor")
}
\references{
Renee de Mutsert, Martin den Heijer, Ton J Rabelink, Johannes WA Smit, Johannes A Romijn, Johan W Jukema, Albert de Roos, Christa M Cobbaert, Margreet Kloppenburg, Saskia le Cessie, Saskia Middeldorp, Frits R Rosendaal. The Netherlands epidemiology of obesity (NEO) study: Study design and data collection. European Journal of Epidemiology (2013). <doi:10.1007/s10654-013-9801-3>

Linda Nab, Maarten van Smeden, Renee de Mutsert, Frits R Rosendaal, and Rolf HH Groenwold. Sampling strategies for internal validation samples for exposure measurement error correction: A study of visceral adipose tissue measures replaced by waist circumference measures. American Journal of Epidemiology (2021). <doi:10.1093/aje/kwab114>
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sodium}
\alias{sodium}
\title{TONE sodium data [outcome-calibration study]}
\format{
A data frame with 1000 rows and 4 variables:
\describe{
  \item{recall}{Sodium intake measured by a 24h recall (mg)}
  \item{diet}{Usual diet or sodium-lowering diet (0 = usual, 1 = sodium-lowering)}
  \item{urinary1}{Sodium intake measured in urine (1st measure, mg)}
  \item{urinary2}{Sodium intake measured in urine (2nd measure, mg)}
}
}
\usage{
sodium
}
\description{
Self-reported sodium intake and urinary sodium in the TONE study, a randomized controlled trial designed to
investigate whether a reduction in sodium intake results in satisfactory blood pressure control.
Two replicate urinary sodium measures were available in 50\% of the subjects included in the trial.
}
\details{
This is a simulated data set inspired by the TONE study <doi: 10.1016/1047-2797(94)00056-y>. A motivating example using the example data can be found here: <doi:10.1002/sim.7011>
}
\examples{
data("sodium", package = "mecor")
}
\references{
Lawrence J Appel, Mark Espeland, Paul K Whelton, Therese Dolecek, Shiriki Kumanyika, William B Applegate, Walter H Ettinger, John B Kostis, Alan C Wilson, Clifton Lacy, and Stephen T Miller. Trial of Nonpharmacologic Intervention in the Elderly (TONE). Design and rationale of a blood pressure control trial. Annals of Epidemiology (1995). <doi: 10.1016/1047-2797(94)00056-y>

Ruth H Keogh, Raymond J Carroll, Janet A Tooze, Sharon I Kirkpatrick, Laurence S Freedman. Statistical issues related to dietary intake as the response variable in intervention trials. Statistics in Medicine (2016). <doi:10.1002/sim.7011>
}
\keyword{datasets}
