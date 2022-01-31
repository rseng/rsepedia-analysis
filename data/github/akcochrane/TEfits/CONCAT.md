<!-- README.md is generated from README.Rmd. Please edit that file -->

# TEfits

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/225967950.svg)](https://zenodo.org/badge/latestdoi/225967950)
[![status](https://joss.theoj.org/papers/0d67da372696cc9a817255858d8bb8a7/status.svg)](https://joss.theoj.org/papers/0d67da372696cc9a817255858d8bb8a7)
[![Build
Status](https://travis-ci.com/akcochrane/TEfits.svg?branch=master)](https://travis-ci.com/akcochrane/TEfits)

## Overview to Time-Evolving fits

Behavioral data is described, interpreted, and tested using indices such
as d prime, mean, or psychometric function threshold. The **TEfits**
package serves to allow the same questions to be asked about
time-evolving aspects of these indices, such as the starting level, the
amount of time that the index takes to change, and the asymptotic level
of that index. Nonlinear regression applied to time-evolving functions
is made as intuitive and painless as is feasible, with many extensions
if desired.

The **TEfits** package has a heavy emphasis on interpretability of
parameters. As far as possible, parameters fit by **TEfits** are meant
to reflect human-interpretable representations of time-evolving
processes. Error functions, nonlinear (“change”) functions linking
predicted values to parameters and time, parameter and prediction
boundaries, and goodness-of-fit indices are intended to be clear and
adjustable. An equal emphasis is on ease of use: minimal arguments are
necessary to begin using the primary functions, `TEfit()` and `TEbrm()`,
and many common tasks are fully automated (e.g., optimization starting
points, bootstrapping).

## Installing the package

The R package `devtools` includes a very easy way to install packages
from Github.

    devtools::install_github('akcochrane/TEfits', build_vignettes = TRUE)

Although having vignettes is nice for exploring the functionality of the
package (via `browseVignettes('TEfits')`), building the vignettes takes
a minute or two. Remove the `build_vignettes = TRUE` argument to speed
up installation.

## Simple model of exponential change

A basic maximum-likelihood model nonlinearly relating time to an outcome
variable. The first argument is a data frame, with the first column
being the response variable and the second column being the time
variable. The model is parameterized in terms of the starting value, the
asymptotic value, and the \[base-2\] log of the time taken to change
halfway from the starting to the asymptotic values.

``` r
library(TEfits)

# generate artificial data:
dat_simple  <- data.frame(response=log(2:31)/log(32),trial_number=1:30)

# fit a `TEfit` model
mod_simple <- TEfit(dat_simple[,c('response','trial_number')])

plot(mod_simple,plot_title='Time-evolving fit of artificial data')
```

![](README_files/figure-markdown_github/model_simple-1.png)

``` r
summary(mod_simple)
```

    ## 
    ## >> Formula: response~((pAsym) + ((pStart) - (pAsym)) * 2^((1 - trial_number)/(2^(pRate))))
    ## 
    ## >> Converged: TRUE 
    ## 
    ## >> Fit Values:
    ##        Estimate
    ## pAsym     1.016
    ## pStart    0.251
    ## pRate     2.866
    ## 
    ## >> Goodness-of-fit:
    ##             err  nullErr nPars nObs     Fval Pval  Rsquared       BIC   nullBIC
    ## ols 0.007789065 1.272257     3   30 2191.575    0 0.9938778 -237.4834 -91.41094
    ##      deltaBIC
    ## ols -146.0724
    ## 
    ## >> Test of change in nonindependence:
    ##                          rawSpearman modelConditionalSpearman
    ## response ~ trial_number:          -1               0.03537264
    ##                          proportionalSpearmanChange pValSpearmanChange
    ## response ~ trial_number:                 0.03537264                  0

Alternatively, a similar model can be fit using the Bayesian package
`brms`. This takes a bit longer, but provides much more flexibility and
information about the model.

``` r
# fit a `TEbrm` model
mod_TEbrm <- TEbrm(response ~ trial_number, dat_simple)
```

``` r
conditional_effects(mod_TEbrm)
```

![](README_files/figure-markdown_github/model_simple_TEbrm_output-1.png)

``` r
summary(mod_TEbrm)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: response ~ pAsym + ((pStart) - (pAsym)) * 2^((1 - trial_number)/(2^(pRate))) 
    ##          pStart ~ 1
    ##          pRate ~ 1
    ##          pAsym ~ 1
    ##    Data: attr(rhs_form, "data") (Number of observations: 30) 
    ## Samples: 3 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 3000
    ## 
    ## Population-Level Effects: 
    ##                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## pStart_Intercept     0.25      0.01     0.23     0.28 1.00     1089     1122
    ## pRate_Intercept      2.87      0.08     2.73     3.04 1.00      768     1127
    ## pAsym_Intercept      1.02      0.01     0.99     1.05 1.00      820     1214
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.02      0.00     0.01     0.02 1.00     1444     1292
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

## Bootstrapped model with Bernoulli error function

An example of a maximum-likelihood fit using a Bernoulli response
distribution, with 40 bootstrapped fits.

``` r
# fit a `TEfit` model
mod_boot <- TEfit(dat_simple[,c('response','trial_number')], 
             errFun='bernoulli',
             bootPars=tef_bootList(resamples = 40))
```

``` r
plot(mod_boot,plot_title='Time-evolving fit of artificial data with 95% CI from 40 bootstrapped fits')
```

![](README_files/figure-markdown_github/model_boot_output-1.png)

``` r
summary(mod_boot)
```

    ## 
    ## >> Formula: response~((pAsym) + ((pStart) - (pAsym)) * 2^((1 - trial_number)/(2^(pRate))))
    ## 
    ## >> Converged: TRUE 
    ## 
    ## >> Fit Values:
    ##        Estimate  Q025  Q975 pseudoSE
    ## pAsym     0.998 0.981 1.000    0.005
    ## pRate     2.741 2.522 2.755    0.059
    ## pStart    0.231 0.179 0.268    0.023
    ## 
    ## >> Goodness-of-fit:
    ##                err  nullErr nPars nObs      BIC  nullBIC    deltaBIC
    ## bernoulli 13.42567 16.83409     3   30 37.05493 37.06937 -0.01444199
    ## 
    ## >> Test of change in nonindependence:
    ##                          rawSpearman modelConditionalSpearman
    ## response ~ trial_number:          -1              -0.04694105
    ##                          proportionalSpearmanChange pValSpearmanChange
    ## response ~ trial_number:                 0.04694105                  0
    ## 
    ## >> Percent of resamples predicting an increase in values: 100 
    ## 
    ## >> Timepoint at which resampled estimates diverge from timepoint 1, with Cohen's d>1: 2 
    ## 
    ## >> Bootstrapped parameter correlations:
    ##         pAsym pStart  pRate    err
    ## pAsym   1.000 -0.217 -0.041 -0.143
    ## pStart -0.217  1.000  0.563  0.368
    ## pRate  -0.041  0.563  1.000 -0.115
    ## err    -0.143  0.368 -0.115  1.000

## Fitting multiple models

An example of fitting a given model to subsets of data (e.g., individual
participants within a behavioral study).

``` r
# generate artificial data:
dat <- data.frame(response=rep(dat_simple$response,4)*seq(0,.2,length=120),trial_number=rep(1:30,4),group=rep(letters[1:4],each=30))

# fit a `TEfitAll` model
mod_4group <- TEfitAll(dat[,c('response','trial_number')], 
             groupingVar = dat$group,
             groupingVarName = 'Participant')
```

    ## 
    ## Your rate is very close to the boundary. Consider penalizing the likelihood.. 
    ## Your rate is very close to the boundary. Consider penalizing the likelihood.. 
    ## Your rate is very close to the boundary. Consider penalizing the likelihood.. .

Note the warnings regarding rate parameters; identifiability is a major
concern in nonlinear models, and `TEfits` attempts to notify the user of
potentially problematic situations.

``` r
plot(mod_4group)
```

![](README_files/figure-markdown_github/plot_model_groups-1.png)

``` r
summary(mod_4group)
```

    ## 
    ## >> Formula: response ~ ((pAsym) + ((pStart) - (pAsym)) * 2^((1 - trial_number)/(2^(pRate))))
    ## 
    ## >> Overall effects:
    ##             pAsym     pStart      pRate
    ## mean   0.14922726 0.01639030 3.83366602
    ## stdErr 0.03933407 0.01060451 0.02431497
    ## 
    ##                 err    nullErr nPars nObs      Fval         Pval   Rsquared
    ## mean   3.005041e-04 0.03071614     3   30 1692.5939 1.110223e-16 0.97598962
    ## stdErr 6.864644e-05 0.01187769     0    0  653.4848 1.110223e-16 0.01661145
    ##                BIC    nullBIC   deltaBIC  linkFun errFun changeFun converged
    ## mean   -337.338970 -211.91820 -125.42077 identity    ols      expo         1
    ## stdErr    6.548355   14.35328   19.26152 identity    ols      expo         0
    ##        pValSpearmanChange
    ## mean                    0
    ## stdErr                  0
    ## 
    ## 
    ## >> Max runs: 200  -- Tolerance: 0.05 
    ## 
    ## >> Parameter Pearson product-moment correlations:

    ##         pAsym pStart  pRate
    ## pAsym   1.000  1.000 -0.757
    ## pStart  1.000  1.000 -0.763
    ## pRate  -0.757 -0.763  1.000

An analogous model, this time fitting “participant-level” models as
random effects within a mixed-effects model, can be implemented using
`TEbrm` (the recommended method).

``` r
mod_4group_TEbrm <- TEbrm(response ~
                            tef_change_expo3('trial_number',parForm = ~ (1|group))
                          ,data = dat
)
```

``` r
conditional_effects(mod_4group_TEbrm)
```

![](README_files/figure-markdown_github/model_groups_TEbrm_output-1.png)

``` r
summary(mod_4group_TEbrm)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: response ~ pAsym + ((pStart) - (pAsym)) * 2^((1 - trial_number)/(2^(pRate))) 
    ##          pStart ~ (1 | group)
    ##          pRate ~ (1 | group)
    ##          pAsym ~ (1 | group)
    ##    Data: attr(rhs_form, "data") (Number of observations: 120) 
    ## Samples: 3 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 3000
    ## 
    ## Group-Level Effects: 
    ## ~group (Number of levels: 4) 
    ##                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## sd(pStart_Intercept)     0.04      0.03     0.01     0.11 1.03       96
    ## sd(pRate_Intercept)      1.73      0.78     0.74     3.69 1.02      264
    ## sd(pAsym_Intercept)      0.10      0.10     0.01     0.31 1.14       15
    ##                      Tail_ESS
    ## sd(pStart_Intercept)      846
    ## sd(pRate_Intercept)      1375
    ## sd(pAsym_Intercept)       103
    ## 
    ## Population-Level Effects: 
    ##                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## pStart_Intercept     0.02      0.02    -0.01     0.08 1.02      473      888
    ## pRate_Intercept      4.62      0.69     3.20     5.89 1.03      121     1329
    ## pAsym_Intercept      0.20      0.06     0.06     0.29 1.18       12       25
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.00      0.00     0.00     0.00 1.11     2691     1791
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

## Using a more common linear regression framework

In some cases (such as `mod_simple` above), similar performance can be
attained using a nonlinear transformation of time as a predictor in a
linear model. This method is plotted in green on top of the `mod_simple`
results, with clearly near-identical fits.

``` r
# Fit a `lm` model, first computing the best nonlinear transformation for time:
mod_lm <- TElm(response~trial_number,dat_simple,timeVar = 'trial_number')

plot(mod_simple)

lines(dat_simple$trial_number,fitted(mod_lm),col='green',lty=2,lwd=2)
```

![](README_files/figure-markdown_github/TElm-1.png)

TElm parameter estimates:

| X.Intercept. | trial_number |  rate |
|-------------:|-------------:|------:|
|        1.018 |       -0.766 | 2.878 |

TEfit parameter estimates:

|          | pAsym | pStart | pRate |
|:---------|------:|-------:|------:|
| Estimate | 1.016 |  0.251 | 2.866 |

Note that `TEfit` provides start and asymptote parameters directly,
while `TElm` provides start as an offset from asymptote (ie.,
`Intercept`).

For extensions of this framework see `TEglm`, `TElmem`, and `TEglmem`.

# Testing functionality

`TEfits` includes automatic testing using the `testthat` package and
[Travis-CI](https://travis-ci.com/github/akcochrane/TEfits). If users
wish to run these tests locally, it’s recommended to download/clone the
repo to a local directory `~/TEfits`. Then install and run tests as
follows:

    devtools::install('~/TEfits') # replace '~' with your filepath

    testthat::test_package('TEfits')

# Performance disclaimer

**TEfits** comes with no guarantee of performance. Nonlinear regression
can be very sensitive to small changes in parameterization, optimization
starting values, etc. No universal out-of-the box implementation exists,
and **TEfits** is simply an attempt to create an easy-to-use and robust
framework for behavioral researchers to integrate the dimension of time
into their analyses. **TEfits** may be unstable with poorly-behaved
data, and using the option to bootstrap models or use Bayesian sampling,
and run slight variations to test for robustness, is generally the best
option for assessing fits. In addition, running the same fitting code
multiple times and comparing fit models should provide useful checks.
All of these things take time, and **TEfits** is not built for speed;
please be patient.

# Community guidelines

If you are having technical difficulties, if you would like to report a
bug, or if you want to recommend features, it’s best to open a Github
Issue. Please feel welcome to fork the repository and submit a pull
request as well.
---
title: 'TEfits: Nonlinear regression for time-evolving indices'
authors:
- affiliation: 1
  name: Aaron Cochrane
  orcid: 0000-0001-6691-9149
tags:
- R
- Psychology
- Nonlinear regression
- Generalized nonlinear model
- Perceptul learning
affiliations:
- index: 1
  name: University of Wisconsin - Madison
output: html_document
date: 29 June 2020
bibliography: references.bib
---

# Summary

Within behavioral science, it is common for data to be paradigmatically collected through repeated measurement of behavior (e.g., on each of 400 trials a human presses one of two buttons to indicate which of two possible stimuli they saw). Typical analytic tools used alongside such designs, such as ANOVA, linear regression, or T tests, implicitly assume that the data arises from distributions that are stationary across the repeated individual measurements (i.e., that every trial is independently and identically sampled from the same distribution, or _iid_, conditional on experimentally manipulated or observed variables).  Interestingly, the use of such analytic tools is common even in those areas of behavioral science that are inherently concerned with time-evolving changes in behavior such as learning, memory, priming, adaptation, vigilance, cognitive control. For instance, in learning research it is common for researchers to first divide the repeated measurements into temporal bins (e.g., trials 1-100; 101-200; 201-300). They then calculate means within those bins, before applying the analysis tools above. Such an analytic method is explicitly modelling a process where performance can change between, but not within, bins. That is, conditional stationarity is assumed within the temporal bins. Beyond these fields, research methods in many others (e.g., attention, development, neuroscience, perception) attempt to by-pass the problem of non-stationarity by utilizing practice trials prior to collecting behavioral data. Practice trials are intended to give participants enough practice with the task that they reach a stable level of performance. However, whether this assumption is true is rarely tested. Our previous work has demonstrated, across several different experimental contexts, that by-trial modeling of performance provides estimates of the full timecourse of behavioral change. In doing so, these models of nonlinear monotonic *trend* stationarity provide both better estimates of behavior, as well as allowing for deeper inferences regarding the underlying processes at work, than statistical methods that assume that behavioral data remains unconditionally stationary over the course of a set of measurements [@kattner_trial-dependent_2017].

`TEfits` is a `R` package for fitting and assessing time-evolving models to data common in behavioral science. `TEfits` is designed with behavioral science researchers with a range of interests and expertise in models of time-dependent changes in behavior. Although many excellent nonlinear regression methods exist in `R`, most notably using the powerful and flexible Bayesian package `brms` [@burkner_brms_2017], but also including functions such as `nls` from the core `R` `stats` package, these methods can be difficult to learn and integrate into the workflow of researchers not familiar with nonlinear regression. The user-oriented functions of `TEfits` are designed to be friendly to `R` users with minimal experience implementing nonlinear models. Extensions of this base functionality allow for simple use of various time-evolving indices (e.g., psychometric function threshold or d prime), objective functions, and/or functional forms of time-related change. Default constraints are applied to models for stability and reproducibility, but boundaries on parameters or predicted values are fully user-defineable. `TEfits` is designed to operate with minimal dependencies on other `R` packages. However, certain functions allow for optional simulation from model fits using `MASS` [@venables_ripley], fitting of hierarchical models using `lme4` [@bates_fitting_2015], or re-fitting a model using Bayesian methods using `brms` [@burkner_brms_2017].

`TEfits` is being actively used in learning and memory research, with several manuscripts in preparation or under review, and results using `TEfits` have been presented at academic conferences. Primary areas of use to-date include assessments of the most appropriate learning functional form of visual perceptual improvements, testing for learning and generalization in the field of radiological diagnosis [@johnston_perceptual], and modeling rapid shifts of attentional control in response to environmental statistics. However, the use of `TEfits` could be appropriate in any domain where individuals repeatedly engage with the same task (i.e., most psychological tasks). The ease of use and wide applicability of `TEfits` models should remove barriers from many behavioral scientists' assessment of the assumption of stationarity. Instead of this assumption users are provided a framework for understanding the changes that occur in nonstationary (i.e., _iid_ conditioned on a time-evolving trend) distributions of behavioral data.

# Funding and Support

This work has been supported in part by US Office of Naval Research Grant ONR-N000141712049.

# References
---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# TEfits

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/225967950.svg)](https://zenodo.org/badge/latestdoi/225967950)
[![status](https://joss.theoj.org/papers/0d67da372696cc9a817255858d8bb8a7/status.svg)](https://joss.theoj.org/papers/0d67da372696cc9a817255858d8bb8a7)
[![Build Status](https://travis-ci.com/akcochrane/TEfits.svg?branch=master)](https://travis-ci.com/akcochrane/TEfits)

## Overview to Time-Evolving fits

Behavioral data is described, interpreted, and tested using indices such as d prime, mean, or psychometric function threshold. The **TEfits** package serves to allow the same questions to be asked about time-evolving aspects of these indices, such as the starting level, the amount of time that the index takes to change, and the asymptotic level of that index. Nonlinear regression applied to time-evolving functions is made as intuitive and painless as is feasible, with many extensions if desired. 

The **TEfits** package has a heavy emphasis on interpretability of parameters. As far as possible, parameters fit by **TEfits** are meant to reflect human-interpretable representations of time-evolving processes. Error functions, nonlinear ("change") functions linking predicted values to parameters and time, parameter and prediction boundaries, and goodness-of-fit indices are intended to be clear and adjustable. An equal emphasis is on ease of use: minimal arguments are necessary to begin using the primary functions, `TEfit()` and `TEbrm()`, and many common tasks are fully automated (e.g., optimization starting points, bootstrapping).

## Installing the package

The R package `devtools` includes a very easy way to install packages from Github.

```
devtools::install_github('akcochrane/TEfits', build_vignettes = TRUE)
```

Although having vignettes is nice for exploring the functionality of the package (via `browseVignettes('TEfits')`), building the vignettes takes a minute or two. Remove the `build_vignettes = TRUE` argument to speed up installation.

## Simple model of exponential change

A basic maximum-likelihood model nonlinearly relating time to an outcome variable. The first argument is a data frame, with the first column being the response variable and the second column being the time variable. The model is parameterized in terms of the starting value, the asymptotic value, and the [base-2] log of the time taken to change halfway from the starting to the asymptotic values.

```{r model_simple, message=FALSE, warning=FALSE}

library(TEfits)

# generate artificial data:
dat_simple  <- data.frame(response=log(2:31)/log(32),trial_number=1:30)

# fit a `TEfit` model
mod_simple <- TEfit(dat_simple[,c('response','trial_number')])

plot(mod_simple,plot_title='Time-evolving fit of artificial data')

summary(mod_simple)
```

Alternatively, a similar model can be fit using the Bayesian package `brms`. This takes a bit longer, but provides much more flexibility and information about the model.

```{r model_simple_TEbrm, message=FALSE, warning=FALSE, results='hide'}

# fit a `TEbrm` model
mod_TEbrm <- TEbrm(response ~ trial_number, dat_simple)
```

```{r model_simple_TEbrm_output, message=FALSE, warning=FALSE}
conditional_effects(mod_TEbrm)

summary(mod_TEbrm)
```

## Bootstrapped model with Bernoulli error function

An example of a maximum-likelihood fit using a Bernoulli response distribution, with 40 bootstrapped fits.

```{r model_boot, message=FALSE, warning=FALSE, results='hide'}
# fit a `TEfit` model
mod_boot <- TEfit(dat_simple[,c('response','trial_number')], 
             errFun='bernoulli',
             bootPars=tef_bootList(resamples = 40))
```

```{r model_boot_output, message=FALSE, warning=FALSE}
plot(mod_boot,plot_title='Time-evolving fit of artificial data with 95% CI from 40 bootstrapped fits')

summary(mod_boot)
```

## Fitting multiple models

An example of fitting a given model to subsets of data (e.g., individual participants within a behavioral study).

```{r model_groups, message=FALSE, warning=FALSE}
# generate artificial data:
dat <- data.frame(response=rep(dat_simple$response,4)*seq(0,.2,length=120),trial_number=rep(1:30,4),group=rep(letters[1:4],each=30))

# fit a `TEfitAll` model
mod_4group <- TEfitAll(dat[,c('response','trial_number')], 
             groupingVar = dat$group,
             groupingVarName = 'Participant')

```

Note the warnings regarding rate parameters; identifiability is a major concern in nonlinear models, 
and `TEfits` attempts to notify the user of potentially problematic situations.

```{r plot_model_groups}

plot(mod_4group)

summary(mod_4group)
```

An analogous model, this time fitting "participant-level" models as random effects within a mixed-effects model, can be implemented using `TEbrm` (the recommended method).

```{r model_groups_TEbrm, message=FALSE, warning=FALSE, results='hide'}

mod_4group_TEbrm <- TEbrm(response ~
                            tef_change_expo3('trial_number',parForm = ~ (1|group))
                          ,data = dat
)
```

```{r model_groups_TEbrm_output, message=FALSE, warning=FALSE}
conditional_effects(mod_4group_TEbrm)

summary(mod_4group_TEbrm)
```

## Using a more common linear regression framework

In some cases (such as `mod_simple` above), similar performance can be attained using a nonlinear transformation of time as a predictor in a linear model. This method is plotted in green on top of the `mod_simple` results, with clearly near-identical fits.

```{r TElm}
# Fit a `lm` model, first computing the best nonlinear transformation for time:
mod_lm <- TElm(response~trial_number,dat_simple,timeVar = 'trial_number')

plot(mod_simple)

lines(dat_simple$trial_number,fitted(mod_lm),col='green',lty=2,lwd=2)

```

TElm parameter estimates:

`r  knitr::kable(round(data.frame(t(c(coef(mod_lm),rate=mod_lm$rate))),3))`

TEfit parameter estimates:

`r knitr::kable(round(data.frame(t(coef(mod_simple))),3))`

Note that `TEfit` provides start and asymptote parameters directly, while `TElm` provides 
start as an offset from asymptote (ie., `Intercept`).

For extensions of this framework see `TEglm`, `TElmem`, and `TEglmem`.

# Testing functionality

`TEfits` includes automatic testing using the `testthat` package and [Travis-CI](https://travis-ci.com/github/akcochrane/TEfits). If users wish to run these tests locally, it's recommended to download/clone the repo to a local directory `~/TEfits`. Then install and run tests as follows:

```
devtools::install('~/TEfits') # replace '~' with your filepath

testthat::test_package('TEfits')
```

# Performance disclaimer

**TEfits** comes with no guarantee of performance. Nonlinear regression can be very sensitive to small changes in parameterization, optimization starting values, etc. No universal out-of-the box implementation exists, and **TEfits** is simply an attempt to create an easy-to-use and robust framework for behavioral researchers to integrate the dimension of time into their analyses. **TEfits** may be unstable with poorly-behaved data, and using the option to bootstrap models or use Bayesian sampling, and run slight variations to test for robustness, is generally the best option for assessing fits. In addition, running the same fitting code multiple times and comparing fit models should provide useful checks. All of these things take time, and **TEfits** is not built for speed; please be patient.

# Community guidelines

If you are having technical difficulties, if you would like to report a bug, or if you want to recommend features, it's best to open a Github Issue. Please feel welcome to fork the repository and submit a pull request as well.
---
title: "TEfit tutorial"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{TEfit tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup,echo=F}
library(TEfits)
```

# Example: Memory span

One context in which a time-evolving model would be fit is within cognitive psychology. For instance, within a test of working memory span, participants may complete a set of trials (e.g., 80) each with some set size (e.g., memory load of 6 items). The task is novel to them and difficult, and it is possible that they are learning within the task. This would manifest as an increase in accuracy with increasing trial number.

As a rough simulation of one participants' data, we will generate by-trial accuracies that exponentially increases from a starting point accuracy of .4 to an asymptotic accuracy of .9.

```{r simulate_data}
set.seed(111)

dat <- data.frame(
  accuracy = rbinom(80,6,
           .9-.5*2^(-rnorm(80,1:80,3)/5) # exponential change from .4 to .9 with a rate [half-change time constant] of log2(5+1)=2.58
           )/6,
  trial_number = 1:80
)
```

The typical analytical approach would be to find the average accuracy (`r round(mean(dat$accuracy),3)`). One alternative would be to remove some early set of trials (e.g., 20) as so-called "practice" and average the remaining trials (`r round(mean(dat$accuracy[21:80]),3)`).

## Fitting a time-evolving mean accuracy

This latter method implicitly acknowledges the possibility of nonstationarity in the measure of interest (accuracy), but it addresses the issue by arbitrarily reducing the amount of data being considered. Alternatively, if the asymptotic accuracy is desired, `TEfits` can be used to fit a time-evolving model and find the asymptotic accuracy:

```{r TEfit_simple}
m_simple <- TEfit(dat[,c('accuracy','trial_number')])
```

```{r TEfit_simple_asymptote}
m_simple$model$par['pAsym']
```

This is more accurate than either approach, while utilizing all data. A more general summary can be viewed:

```{r TEfit_simple_summary}
summary(m_simple)
```

`Formula` shows the equation that was fit in the model

`Converged` indicates whether the model converged, as measured by similar parameter values being found by the optimization runs with the lowest errors. *This is a heuristic*, and should ideally be corroborated with parameter distributions from bootstrapped fits.

`Fit Values` prints the fit parameter values. In 3-parameter exponential change, `rate` is log~2~ of the amount of time taken for 50% of change to occur

`Goodness-of-fit` prints several fit indices (e.g., the `err` and `BIC` of the best fit). Also printed are the corresponding fit indices for the null model (i.e., the model that is not time-evolving). The row name indicates the error function (e.g., `ols` or `logcosh`).

`Test of change in nonindependence` prints the rank correlation between the null model residuals and the time variable, the rank correlation between the full model residuals and the time variable, the ratio between the absolute values of these correlations, and the p value as calculated with `psych::r.test()`. If the `tseries` package is available, the p-value from `kpss.test` of stationarity will also be included.

A plot can also be viewed:

```{r TEfit_simple_plot, fig.height=5, fig.width=6}
plot(m_simple)
```

## Adjusting the model to better fit the data: Error function

By default `TEfit` minimizes the sum of squared error (`ols`) between model predictions and the response variable. This allows a fit of the conditional mean of the response variable, and it is directly comparable to "normal" methods of linear regression (e.g., `lm`; see `?TElm`). However, analysts' knowledge about their data may allow for more specific choices. For example, when considering accuracy on a memory span task, this accuracy is bounded at 0 and 1 (i.e., 0% and 100%). A `bernoulli` error function is therefore appropriate whereas an `ols` error function is not (e.g., due to skewed residuals and the ability to evaluate the error at predictions above 1 or below 0). The previous model can be fit using a `bernoulli` error function instead:

```{r TEfit_bernoulli, message=FALSE, warning=FALSE}
m_bernoulli <- TEfit(dat[,c('accuracy','trial_number')],errFun = 'bernoulli')

```
```{r TEfit_bernoulli_summary}
summary(m_bernoulli, printOutput = F)[c("param_vals",'GoF')]
```

Note the different scale of error values here, when using the negative log likelihood of the `bernoulli` distribution, as compared to the `ols` fit.

## Beyond point estimates: Bootstrapped fits

Parameter distributions can be described by resampling the data with replacement and fitting the model to each dataset (see `?tef_bootList`). Here we will use 50 resamples for the sake of speed; using fewer than 500 for final inferences is not recommended, and an order of magnitude higher may be necessary if categorical decisions (e.g., parameter "significance") are near decision threshold values.

```{r TEfit_bootstrap, message=FALSE, warning=FALSE}
m_bootstrap <- TEfit(dat[,c('accuracy','trial_number')],
                     errFun = 'bernoulli',
                     bootPars = tef_bootList(resamples = 50))

```
```{r TEfit_bootstrap_summary, fig.height=5, fig.width=6}
summary(m_bootstrap)
plot(m_bootstrap)
```

Relative to the previous model additional information is provided regarding .025 and .975 quantiles of parameter estimates as well as the pseudo standard error (assuming that CI came from normal distribution; see `?summary.TEfit`). Three new outputs are also included which should be somewhat self-explanatory: `Percent of resamples predicting an increase in values` and `Timepoint at which resampled estimates diverge from timepoint 1, with Cohen's d>1` both aid in interpreting the strength of time-evolving trends. `Bootstrapped parameter correlations` provides a glimpse at the dimensionality of the parameterization.

## Adjusting the model to better fit the data: Limiting parameter values

If we imagine that our working memory task involves 9 possible options with replacement (e.g., the numerals 1 through 9), we know that the guessing rate is 1/9. As with the `bernoulli` distributional information that can be incorporated into a TEfits model, *a priori* bounds to predicted values can be incorporated. If accuracies are theoretically unable to fall below 1/9 or above 1 we can re-fit the previous model with this information.

```{r TEfit_bounded, message=FALSE, warning=FALSE}
m_bounded <- TEfit(dat[,c('accuracy','trial_number')],
                   errFun = 'bernoulli',
                   bootPars = tef_bootList(resamples = 50),
                   control = tef_control(y_lim = c(1/9,1))
                   )

```
```{r TEfit_bounded_summary, fig.height=5, fig.width=6}
summary(m_bounded, printOutput = F)[c("param_vals",'GoF')]
plot(m_bounded)
```

Given the current dataset this should not make much of a difference. However, *a priori* prediction bounds can often provide powerful theory-based constraints on parameters.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_subsetCorrel.R
\name{tef_subsetCorrel}
\alias{tef_subsetCorrel}
\title{Find mean correlations of subsets of vectors}
\usage{
tef_subsetCorrel(
  x,
  y,
  method = "spearman",
  returnVal = "mean",
  subsample = "even_odd",
  iter = 200
)
}
\arguments{
\item{x}{vector to correlate}

\item{y}{vector to correlate}

\item{method}{correlation method (e.g., "spearman" or "pearson")}

\item{returnVal}{return "mean" or (if subsample is random) "quantiles"}

\item{subsample}{should the subsamples be "even_odd" or "random" (half of the first half, half of the second half)}

\item{iter}{if subsample is "random", number of times the vectors are split, each half is subset, those subsets are correlated, and those correlations are averaged}
}
\description{
Selects subsets of x and y vectors,
and calculates the mean correlations between the two subsets.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_change_power4.R
\name{tef_change_power4}
\alias{tef_change_power4}
\title{Construct a 4-parameter power function of change}
\usage{
tef_change_power4(
  timeVar,
  parForm = ~1,
  startForm = ~1,
  rateForm = ~1,
  asymForm = ~1,
  prevTimeForm = ~1,
  rateBase = 2,
  propRemain = 0.25
)
}
\arguments{
\item{timeVar}{String. The name of the variable in the model that corresponds to time. The variable of time should be positive and numeric, and the function of change should be expected to happen with increasing time.}

\item{parForm}{The right-hand side of the formula defining all nonlinear parameters as well as the null [non-time-varying] model.}

\item{startForm}{The right-hand side of the formula defining the start parameter. Overwrites \code{parForm} for this parameter.}

\item{rateForm}{The right-hand side of the formula defining the rate parameter. Overwrites \code{parForm} for this parameter.}

\item{asymForm}{The right-hand side of the formula defining the asymptote parameter. Overwrites \code{parForm} for this parameter.}

\item{prevTimeForm}{The right-hand side of the formula defining the log of the "previous time" parameter. Overwrites \code{parForm} for this parameter. The base of the log is \code{rateBase}.}

\item{rateBase}{Number. The base of the log (e.g., 2 or \code{exp(1)}) of the rate [time constant].}

\item{propRemain}{Change rate is parameterized in terms of the \code{rateBase} log of time to this proportion of change remaining (i.e., \code{1-propRemain} of total change occurs in \code{rateBase^rateParameter} time.}
}
\description{
By defining the model variable associated with time (e.g., trial number), and
formulas defining each of the nonlinear parameters of time-related change,
this function constructs a model that can then be passed to functions for fitting
the model (e.g., \code{\link{TEbrm}}).
}
\details{
Function is \strong{under development}
and is likely to be buggy, and to change frequently.
}
\examples{
equation_to_fit <- tef_change_power4('timeVar',parForm = ~ xvar1*xvar2) # both variables should be numeric for TEfit methods! TEbrm should work with factors as well
}
\seealso{
\code{\link{TEbrm}} for examples of how to use this function in specifying models.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.TEfitAll.R
\name{summary.TEfitAll}
\alias{summary.TEfitAll}
\title{Summarize a collection of time-evolving fits}
\usage{
\method{summary}{TEfitAll}(TEs3s, printOutput = T, printAll = F)
}
\arguments{
\item{TEs3s}{A set of fit TEfit models (output from TEfitAll())}

\item{printOutput}{Print output to console (if T) or return a list of summary items (if F)}

\item{printAll}{Print grouping-level fits?}
}
\description{
Summarize a collection of time-evolving fits
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_runningMean.R
\name{tef_runningMean}
\alias{tef_runningMean}
\title{1-dimensional Gaussian smoothing}
\usage{
tef_runningMean(x, k_hwhm = 2, distr = "gaussian")
}
\arguments{
\item{x}{vector to be smoothed}

\item{k_hwhm}{The half-width half-max of the kernal [when Gaussian = sd*1.17741]. This is the index distance at which an element receives half the weight of the element at the center of the smooothing.}

\item{distr}{The distribition's density function to be used. 'gaussian' or 'cauchy'}
}
\description{
Runs a Gaussian or Cauchy density estimate centered over each
element [numeric, logical, or NA]
of the vector, and calculates
a density-weighted average for that element's index.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted.TEfitAll.R
\name{fitted.TEfitAll}
\alias{fitted.TEfitAll}
\title{Get fitted values and summary statistics from a set of TEfit models}
\usage{
\method{fitted}{TEfitAll}(TEs3s)
}
\arguments{
\item{TEs3s}{A set of models fit by TEfitAll()}
}
\description{
Get fitted values and summary statistics from a set of TEfit models
}
\examples{
\dontrun{
m <- TEfitAll(anstrain[,c('acc','trialNum')],groupingVar = anstrain$subID)
fitted_data <- fitted(m)
plot(fitted_data$meanPred)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEfitAll.R
\name{TEfitAll}
\alias{TEfitAll}
\title{Fit several time-evolving regression models}
\usage{
TEfitAll(
  varIn,
  groupingVar,
  groupingVarName = "grouping_var",
  returnAll = T,
  progressDot = T,
  linkFun = list(link = "identity"),
  errFun = "ols",
  changeFun = "expo",
  bootPars = tef_bootList(),
  blockTimeVar = NULL,
  covarTerms = list(),
  control = tef_control()
)
}
\arguments{
\item{varIn}{Data frame or vector. First column [or vector] must be the time-dependent response variable (left hand side of regression). If available, second column must be the time variable. All other columns are covariates, possibly involved in a link function.}

\item{groupingVar}{Variable (e.g., participant ID) with which to separate TEfit models. Length must be nrows(varIn)}

\item{groupingVarName}{Name of grouping variable}

\item{returnAll}{Logical. Return only a summary (when T), or that summary plus every model, in a list (when F)}

\item{progressDot}{If TRUE, prints a dot after each group fit}

\item{linkFun}{A list defining a link function (i.e., 'identity', 'd_prime', 'weibull', or 'logistic')}

\item{errFun}{A string defining an error function (e.g., 'ols', 'logcosh', 'bernoulli').}

\item{changeFun}{A string defining the functional form of change (e.g., 'expo', 'power', 'weibull')}

\item{bootPars}{A list defining the details for bootstrapped fits. Defaults to no bootstrapping. Necessary for estimates of uncertainty around fits and for covariance between parameters.}

\item{blockTimeVar}{A string identifying which covariate is the time points of sub-scales (e.g., "blocks" of times within the overall timescale of data collection)}

\item{covarTerms}{An optional list of logical vectors indicating whether parameters should vary by covariates. See examples.}

\item{control}{A list of model parameters. Use of tef_control() is highly recommended.}
}
\description{
A wrapper for fitting a \code{\link{TEfit}} model to the data
for every unique value of groupingVar. Defaults to
returning a list including two summaries and all models; returning
only a summary is also an option. Most arguments (except, e.g., \code{groupingVar}, a grouping vector)
are identical to, and are passed directly to, \code{\link{TEfit}}.
}
\examples{
\dontrun{
m <- TEfitAll(anstrain[,c('acc','trialNum')],groupingVar = anstrain$subID,groupingVarName = 'subID',bootPars = tef_bootList(resamples = 20))
summary(m)
}

}
\seealso{
\code{\link{TEfit}} for fitting a single model;
\code{\link{tef_fitAll2brms}} to re-fit the TEfitAll output using \code{\link[brms]{brms-package}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEbrm_advi.R
\name{TEbrm_advi}
\alias{TEbrm_advi}
\title{Run a brm model with ADVI}
\usage{
TEbrm_advi(
  formIn,
  dataIn = data.frame(),
  ...,
  algorithm = "fullrank",
  conv_thresh = 0.5,
  quiet = F
)
}
\arguments{
\item{formIn}{Model formula, as in \code{\link[brms]{brm}}.}

\item{dataIn}{Data, as in \code{\link[brms]{brm}}.}

\item{...}{Any other argument to pass to \code{\link[brms]{brm}}.}

\item{algorithm}{Which ADVI algorithm to use: "meanfield" or "fullrank".}

\item{conv_thresh}{Re-fit models are compared, with the standardized distance (mean_diff / SD) being calculated. Models keep being re-fit until at least 2 models' largest standardized differences are smaller than this value (or until a certain number of models has been fit in total, which scales inversely with this value). The better-fitting of these two models is then returned. Values over 30 will cause an error, which should not be an issue for any normal use of this function.}

\item{quiet}{Progress is printed by default, but can be suppressed with quiet=T.}
}
\description{
Uses Stan's stochastic gradient ascent methods "fullrank" or "meanfield" rather than full
Bayesian sampling. Is likely to be faster than typical sampling for large models, but possibly
less accurate. This function fits the model several times and returns the best model
(initial model selection uses Bayesian R-squared, final model selection uses 10-fold
cross-validation). \strong{Beware: This re-fitting
has a tendency to crash R sometimes.}
}
\details{
Stochastic gradient ascent in Stan uses Automatic Differentiation Variational Inference (ADVI).
}
\note{
Re-fitting a model within this function is not comprehensive. If using ADVI, it is recommended
to use \code{TEbrm_advi} multiple times, and choose the best using comparisons of fit, such as
\code{fit_model$criteria$kfold} (estimated using
 \code{\link[brms]{add_criterion}(fit_model,'kfold')}).
}
\examples{
\dontrun{
## can be used in place of brm
m <- TEbrm_advi(ratio ~ resp, anstrain_s1)
summary(m)
conditional_effects(m)

## Use in the context of TEfits
m1 <- TEbrm(
acc ~ tef_change_expo3('trialNum')
,data = anstrain_s1,
,algorithm = 'meanfield'
)

}
}
\seealso{
\code{\link[rstan]{vb}} and \link{http://mc-stan.org}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_checkPars.R
\name{tef_checkPars}
\alias{tef_checkPars}
\title{Check for bound hitting and other undesireable outcomes within an optim() call}
\usage{
tef_checkPars(
  err,
  guesses,
  curDat,
  pNames,
  evalFun,
  errFun,
  respVar,
  linkFunX = NA,
  y_lim,
  rate_lim,
  shape_lim,
  penalizeRate,
  paramTerms,
  guessGroups = NULL
)
}
\arguments{
\item{err}{Error}

\item{guesses}{Parameter values}

\item{curDat}{Data being fit}

\item{pNames}{Parameter names}

\item{evalFun}{Function being fit}

\item{errFun}{Function to calculate error}

\item{respVar}{Name of the response variable}

\item{linkFunX}{If relevant, the "x" value for a link function (e.g., Weibull, logistic)}

\item{y_lim}{Limits to fit values}

\item{rate_lim}{Limits to rate parameter}

\item{shape_lim}{If using a Weibull change function, limits to Weibull shape parameter}

\item{penalizeRate}{Logical. Should error be penalized if rate is extremely close to the bounds?}

\item{paramTerms}{parameter-level regressions, to be evaluated for checking y_lim and rate_lim}

\item{guessGroups}{deprecated}
}
\description{
TEfits internal.
}
\details{
Sane boundaries for parameters are the only way that many nonlinear regression optimizations
can be identifiable. Fortunately, theory-driven constraints on parameter ranges provide useful
a priori restrictions on the possible ranges for parameters and model predictions.
This function checks the following:
 \itemize{
\item{\code{start and asymptote parameters} -- all models are parameterized in terms of
starting and ending values. This ensures that the starting and ending values comply with
the \emph{y_lim} boundaries; \emph{y_lim} may be user-defined, defined by another model feature (e.g.,
\emph{bernoulli} error function is limited to predicted values of 0 or 1; Weibull link thresholds must be
above 0).}
\item{\code{rate parameter} -- If not user-input, then defined by \code{TEfits::tef_getLinkedFun}.
Defaults, with exponential change, to a minimum that would provide 50% of model change in \code{sd(timeVar[1:10])}
amount of time, and to a maximum value that would provide 80% of model change at \code{max(timeVar)}. Other change
functions have limits that, with their respective parameterizations, are intended to imitate the limits of the
3-parameter exponential (i.e., imitate the overall shape of the curve's extremes) \strong{These are
heuristics.} The default values are intended to be flexible while maintaining a sufficiently constrained curve such that
\emph{both} starting asymptote parameters are interpretable (i.e., if the rate parameter, which is a time constant,
were to be extremely small, the start parameter could become infinitely large or small). If a time-evolving process
occurs on a timescale that cannot be fit by the default boundaries, it is likely that the data is unsufficient to
characterize that process.}
\item{\code{pPrevTime parameter} -- 4-parameter Power change utilizes a so-called "previous learning time" parameter
that assumes that the learning function extents \emph{backward} through time. This parameter must be greater than
0 and smaller than 1*10^5}
}

Users are highly encouraged to use their own boundaries (e.g., \code{y_lim} & \code{rate_lim}), given knowledge of a specific dataset, using \code{\link{tef_control}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.TEfit.R
\name{coef.TEfit}
\alias{coef.TEfit}
\title{Extract coefficients from a TEfit model}
\usage{
\method{coef}{TEfit}(TEs3)
}
\arguments{
\item{TEs3}{A TEfit model, possibly with resampled parameters}
}
\description{
Extract coefficients from a TEfit model
}
\examples{
\dontrun{
m <- TEfit(anstrain_s1[,c('acc','trialNum')], bootPars = tef_bootList(resamples = 50))
coef(m)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEhelp.R
\name{TEhelp}
\alias{TEhelp}
\title{Get help for the TEfits package}
\usage{
TEhelp(topic = "package")
}
\arguments{
\item{topic}{Help topic.}
}
\description{
Many options for fitting using TEfits are specified in a variety of places. To help,
this function provides verbose explanations of how to use various tools. Topics will be
gradually added, with no guarantee that a given topic has an entry yet.
}
\examples{
## Brief package summary:
TEhelp('package')

## List topics with TEhelp entries:
TEhelp('topics')

## Overview of error functions:
TEhelp('errFun')

## The log-cosh error function:
TEhelp('errFun=logcosh')

## The ex-Gaussian error function, with change in the exponential component:
TEhelp('errFun=exGauss_tau')

## The d prime link function:
TEhelp('linkFun=d_prime')

## The summary method for TEfit models:
TEhelp('summary(TEfit)')

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_applyWDM.R
\name{tef_applyWDM}
\alias{tef_applyWDM}
\title{Calculate the Wiener negative log likelihood of a vector of data using a vector of drift rates.}
\usage{
tef_applyWDM(dat, DR, BS, NDT, Bias = 0.5)
}
\arguments{
\item{dat}{Vector of data (e.g., RT). Hitting upper bound must be positively signed; hitting lower bound must be negatively signed.}

\item{DR}{Vector of drift rates. Must be the same length as dat.}

\item{BS}{Boundary separation parameter. Positive scalar.}

\item{NDT}{Non-decision time parameter. Positive scalar.}

\item{Bias}{Bias parameter. Positive scalar.}
}
\description{
\code{\link{TEfit}} internal. Currently not supported.
Within a TEfit call, errFun='wiener_df' is likely to be extremely slow, e.g., 70 seconds per run.
}
\examples{
dat <- c(-1,1,-.5,.5)
driftRate <- c(.2,.3,.4,.5)
tef_applyWDM(dat,DR=driftRate,BS=1,NDT=.3,Bias=.5)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEbrm.R
\name{TEbrm}
\alias{TEbrm}
\title{Fit a time-evolving model with Stan using brms}
\usage{
TEbrm(
  formula,
  data,
  ...,
  chains = 3,
  priorIn = c(),
  algorithm = "sampling",
  link_start_asym = "",
  quiet = F,
  tef_control_list = tef_control()
)
}
\arguments{
\item{formula}{A formula, with the time-varying response variable on the left, followed by \code{~}.  The right side must be either [A] a single variable name corresponding to the dimension of time, or [B] a call to a \code{TEfits} constructor function such as \code{\link{tef_change_expo3}}. See examples.}

\item{data}{Data frame, from which to fit the model.}

\item{...}{Further arguments passed to the \code{\link[brms]{brm}s} model}

\item{chains}{Number of chains to run the model.}

\item{priorIn}{Optional argument to pass priors to the \code{brms} model, alongside the TEfit-default rate prior. If you provide any, you will likely need to provide priors for all nonlinear parameters. \code{brm} error messages tend to be very helpful in this regard. For more explicit and full control of priors, define all desired priors directly with the \code{prior} argument (which passes straight to \code{brm} and overwrites all other defined priors).}

\item{algorithm}{The algorithm to use when fitting the \code{\link[brms]{brm}} model. See \code{\link{TEbrm_advi}} for warnings and implementation of 'meanfield' or 'fullrank.'}

\item{link_start_asym}{Inverse of the link function to use for the start and asymptote parameters. Defaults to what is passed from formIn. Otherwise, the user would most likely to want to use 'exp' or 'inv_logit'. Refer to examples, and to the recommendation to use an "identity" link function rather than a statistical family's default link.}

\item{tef_control_list}{A list of control parameters passed in by \code{tef_control()}}
}
\description{
Formats and runs a \code{\link[brms]{brm}s} model for a time-evolving nonlinear regression.
This is the recommended way to fit models using \code{TEfits}.
Function is \strong{under development}
and is likely to change frequently.
}
\details{
The variable of time should be positive and numeric, with the nonlinear model providing a regression as a
function of that time variable.

The default number of iterations and chains is small, and intended largely for testing model specifications.
For final inferences from a model, it is highly recommended to run a model for many more iterations
(e.g., \code{iter=5000} or \code{iter=10000}).

When specifying statistical families, it is \emph{extremely highly recommended} to specify an "identity" link function,
and then [if appropriate] specifying a link function using the \code{link_start_asym} argument. See example.

Rates [time constants] are estimated on log scales. Within the exponent of the log estimation, the rates have an additive
offset corresponding to the median of the time variable; this allows the rate parameter priors and estimation to be zero-centered.

Currently supported model constructor functions are:
\itemize{
\item{\code{\link{tef_change_expo3}} -- 3-parameter exponential (start, [inverse] rate, and asymptote) -- rate is log of time to some proportion remaining, default is \code{log2} of time to 50 percent remaining}
\item{\code{\link{tef_change_weibull}} -- 4-parameter weibull (start, [inverse] rate, asymptote, and shape) -- Augmented exponential function. Rate is log of time to some proportion remaining, default is \code{log2} of time to 50 percent remaining. Shape indicates acceleration (>0) or deceleration (<0) of the hazard function (which is constant in an exponential function).}
\item{\code{\link{tef_change_power3}} -- 3-parameter power (start, [inverse] rate, and asymptote) -- rate is log of time to some proportion remaining, default is \code{log2} of time to 25 percent remaining, in order to have parameter ranges be similar to the 3-parameter exponential.}
\item{\code{\link{tef_change_power4}} -- 4-parameter power (start, [inverse] rate, asymptote, and "previous learning time") -- Augmented power function, with a "amount of previous learning time" parameter that is estimated on the same scale as \code{rate}. Rate is log of time to some proportion remaining, default is \code{log2} of time to 25 percent remaining, in order to have parameter ranges be similar to the 3-parameter exponential.}
}

Currently supported link functions are:
\itemize{
\item{\code{\link{tef_link_logistic}} -- logistic psychometric function, parameterized in terms of threshold and bias. Takes the output of a \code{tef_change_} function and augments it before passing to \code{TEbrm}}
\item{\code{\link{tef_link_weibull}} -- weibull (AKA Quick) psychometric function, parameterized in terms of threshold and shape. Takes the output of a \code{tef_change_} function and augments it before passing to \code{TEbrm}}
}
}
\note{
Default priors and parameter boundaries are implemented, but all users would benefit from
re-fitting models with various different priors in order to ensure that inferences are not biased
by defaults. The assumptions that guided the creation of the default priors may not be
appropriate for your data. If no \code{link_start_asym} function is used (i.e., the default
'identity') then start and asymptote priors are Gaussian, with the response variable's mean and double
the response variable's SD. If a \code{link_start_asym} function is used (e.g., 'exp' or
'inv_logit') then start and asymptote priors are Gaussian with a mean of zero and a SD of 3 (which
may place quite a bit of the prior density at values more extreme than appropriate for many
users' data). Default [log time constant] rate parameter's prior is Gaussian centered at the
log of the mean of the time variable (with the base of the log defined in \code{tef_control_list}).
The SD of this prior is 1/3 of the mean, and boundaries are implemented at extreme values. All
other priors are \code{\link[brms]{brm}} defaults.
Use \code{\link[brms]{prior_summary}} to examine priors from a fitted model object; see
\code{\link[brms]{set_prior}} for setting priors.

It is \emph{highly recommended} that, if additional customization is desired (e.g., regarding priors),
to run \code{TEbrm} to create a model with a small number of iterations that is "close enough" to the
desired model specification, then use \code{\link[brms]{update.brmsfit}} to "fine-tune" your model directly (see example 10).
}
\examples{
\dontrun{
#-- #-- Example 01: Simple model
#> Default model formula is exponential change, with no covariates or random effects
m1 <- TEbrm(
  acc ~ trialNum   #> equivalent to `acc ~ tef_change_expo3('trialNum')`
  ,data = anstrain_s1
)

prior_summary(m1)
summary(m1)
conditional_effects(m1)
hypothesis(m1,'pAsym_Intercept > pStart_Intercept') #> Test for learning, i.e., whether asymptote was reliably higher than start. With this limited data, that difference is not reliable

#-- #-- Example 02: Random effects
#> using the tef_change_expo3 function to construct the model formula, with priors, and fixed and random effects
m2 <- TEbrm(
  acc ~ tef_change_expo3('trialNum',parForm = ~ sizeRat + (1|subID))
  ,data = anstrain
  ,priorIn = prior(normal(.5,.5),nlpar='pAsym') + prior(normal(.5,.5),nlpar='pStart')   #> for demonstration, also include non-default priors
)

#-- #-- Example 03: Bernoulli family
#> Estimate accuracy using a more appropriate [bernoulli] response function,
#> > and also estimate the start and asymptote parameters using invert-logit links
m3 <- TEbrm(
  acc ~ tef_change_expo3('trialNum')
  ,data = anstrain_s1
  ,link_start_asym = 'inv_logit'
  ,family=bernoulli(link='identity')
)

#-- #-- Example 04: Logistic PF
#> Fit a time-evolving logistic mixed-effects model (see, e.g., Cochrane et al., 2019, AP&P, 10.3758/s13414-018-01636-w).
#> > May take a few minutes to run.
m4 <- TEbrm(
  resp ~ tef_link_logistic(
     tef_change_expo3('trialNum', parForm = ~ (1|subID))
     ,linkX = 'ratio' )
  ,family=bernoulli(link='identity')
  ,iter = 4000   #> most models, in practice, will need more than the default number of iterations in order to converge as well as have a sufficient effective sample size (ESS)
  ,data = anstrain
)

summary(m4) #> note the `exp` inverse link function on pStartXform and pAsymXform(i.e., log link for threshold values)
conditional_effects(m4, 'ratio:trialNum') #> The psychometric function steepens with learning; see ?brms::conditional_effects
cat(attr(m4$right_hand_side,'link_explanation')) #> An explanation of the link function is included

#-- #-- Example 05: Weibull PF
#> Model change in a Weibull psychometric function's threshold
#> > (learning is change in the absolute stimulus strength at which accuracy is 75\%)
d_tmp <- anstrain_s1   #> make temporary data
d_tmp$absRat <- abs(d_tmp$ratio)   #> calculate absolute stimulus strength
m5 <- TEbrm(
  acc ~ tef_link_weibull(
    tef_change_expo3('trialNum'),linkX = 'absRat')
  ,data = d_tmp
)

#-- #-- Example 06: d prime
#> Model change in d-prime
d_tmp <- anstrain_s1 #> make temporary data
d_tmp$dprime <- tef_acc2dprime(d_tmp$acc, d_tmp$ratio > 0)   #> calculate by-trial d-prime; it isn't prototypical data for d-prime because `resp` doesn't categorize `acc` into "hits" and "false alarms", but it's a decent demonstration of the method.
m6 <- TEbrm(
  dprime ~  tef_change_expo3('trialNum')
  ,data = d_tmp
)

#-- #-- Example 07: Power change
#> Rather than a 3-parameter exponential function of change, use a 3-parameter power function of change.
#> > Also, include a covariate for learning rate.
m7 <- TEbrm(
  acc ~ tef_change_power3('trialNum'
      ,rateForm = ~ sizeRat)
  ,data = anstrain_s1
)

#-- #-- Example 08: Weibull change
#> Model learning as a Weibull function of time (3-parameter exponential with one additional "acceleration" or "deceleration" parameter, "pShape")
#> > May take a few minutes to run.
m8 <- TEbrm(
  acc ~ tef_change_weibull('trialNum')
  ,anstrain_s1
)

#> Test whether learning is accelerating or decelerating, relative to a 3-parameter exponential
hypothesis(m8, 'pShape_Intercept = 0')    #> acceleration or deceleration is not evident in this task on this timescale

#-- #-- Example 09: fixing a parameter
#> Fix a parameter to a constant rather than estimating it (here, fix asymptotic accuracy to 90 percent)
m9 <- TEbrm(
  acc ~ tef_change_expo3('trialNum'
                      ,asymForm = .9)
  ,data = anstrain_s1
)

#-- #-- Example 10: updating a model
#> Fit a preliminary model with few iterations, adjust it to a different family, and update it with more iterations and chains
d_tmp <- anstrain_s1 #> make temporary data
d_tmp$acc_smooth <- tef_runningMean(d_tmp$acc)   #> get smoothed accuracy
m10_initial <-  TEbrm(
  acc_smooth ~ tef_change_expo3('trialNum')
  ,data = d_tmp
  ,iter = 200
)

m10_final <- update(
 m10_initial
 ,family=student() #> heavier-tailed (t distribution) regression
 ,iter = 4000
 ,chains=4
)

#-- #-- Example 11: Generalized Additive Models
#> Fit a generalized additive mixed-effects model with accuracy varying by time.
m11 <- TEbrm(acc ~ tef_change_gam('trialNum', groupingVar = 'subID')
,data = anstrain)

#-- #-- Example 12: Using a variational algorithm rather than full sampling
m12 <- TEbrm(acc ~ tef_change_expo3('trialNum')
,algorithm = 'fullrank'
,data = anstrain_s1)

}
}
\seealso{
For additional flexibility, and full explanations of model options, see \code{\link[brms]{brms-package}}. In
particular, the "..." argument that is passed to \code{\link[brms]{brm}} includes many important options, such
as increasing iterations (argument \code{iter})or parallelization (e.g., \code{cores = 3}).

For other approaches to time-evolving models, see \code{\link{TEfits}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_fitAll2brms.R
\name{tef_fitAll2brms}
\alias{tef_fitAll2brms}
\title{Refit a TEfitAll model with brms}
\usage{
tef_fitAll2brms(
  TEs3s,
  fixef = NA,
  nIter = 2000,
  nChains = 3,
  nCores = 2,
  errFun = NA,
  prior_dispersion = 2
)
}
\arguments{
\item{TEs3s}{TEfitAll model}

\item{fixef}{Parameters vary as random effects by the TEs3s grouping variable. However, if you have main effects (e.g., group differences), enter them \emph{as a data frame} here.}

\item{nIter}{number of iterations}

\item{nChains}{number of chains}

\item{nCores}{number of cores}

\item{errFun}{the error function to use. Defaults to the same as the TEfitAll model, if possible.}

\item{prior_dispersion}{This number, multiplied by the SD of each TEfitAll parameter, is used as the prior SD for that parameter.}
}
\value{
A \code{\link[brms]{brms-package}} nonlinear mixed-effects model object.
}
\description{
\emph{This method has been superceded by} \code{\link{TEbrm}}. \emph{Please use
that method instead.}
}
\details{
Passes a \code{\link{TEfitAll}} model to [nonlinear mixed-effects Bayesian] fitting using
\code{\link[brms]{brms-package}}. Note that, due to the extensive time needed to
fit \code{\link[brms]{brms-package}} models,
this function is less tested than most functions in the \code{TEfits} package. Functionality is
\strong{experimental}.

Priors for nonlinear parameters are informed by the distributions of parameters in the \code{TEfitAll} object [models].
However, any fixed effects should be minimally influenced by these priors

\code{TEfitAll} \code{bernoulli} models are fit using either \code{bernoulli} or \code{Beta} response
distributions in \code{brms} depending on whether the \code{TEfitAll} distrIibution is
binary. \code{TEfitAll} \code{logcosh} models are fit using a \code{asym_laplace} response distribution
in brms predicting the .5 quantile.

If sampling issues occur, increased number of iterations are recommended. Also, running one chain at a time
may help; these models should later be merged using \code{brms::combine_models()}.
}
\note{
Under development. Partial functionality.
}
\examples{
\dontrun{
dat <- anstrain
dat$condition <- rep(c('A','B'),each=500)

# Model with time and one categorical fixed effect
mod_tef <- TEfitAll(dat[,c('acc','trialNum')], groupingVar = dat$subID)
mod_brm <- tef_fitAll2brms(mod_tef,nChains=1,fixef=data.frame(condition=dat$condition))

# Model with time, one categorical fixed effect, and one by-groupingVar (subID) random slope
dat$absRat <- scale(abs(dat$ratio))
mod_tef <- TEfitAll(dat[,c('acc','trialNum',"absRat")], groupingVar = dat$subID,covarTerms=list(pRate=c(F)))
mod_brm <- tef_fitAll2brms(mod_tef,nChains=1,fixef=data.frame(condition=dat$condition))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEglm.R
\name{TEglm}
\alias{TEglm}
\title{Generalized linear model with nonlinear time predictor}
\usage{
TEglm(
  formIn,
  dat,
  timeVar,
  family = gaussian,
  startingOffset = T,
  fixRate = NA
)
}
\arguments{
\item{formIn}{model formula, as in glm()}

\item{dat}{model data, as in glm()}

\item{timeVar}{String. Indicates which model predictor is time (i.e., should be transformed)}

\item{family}{passed to glm()}

\item{startingOffset}{By default (if T) time is coded to start at 1 and saturate to 0. If startingOffset is F, time starts at 0 and saturates to 1. May assist in interpreting interactions with other variables, etc.}

\item{fixRate}{If numeric, use this as a rate parameter [binary-log of 50 percent time constant] rather than estimating it (e.g., to improve reproducibility)}
}
\description{
Fit a generalized linear model with time as a covariate,
while estimating the shape of the nonlinear interpolation between starting and ending time.
First resamples data with replacement 200 times, and each time estimates the best-shaped curve to interpolate
between initial time-related offset and asymptotic time (i.e., rate at which effect of time saturates at zero).
Then uses the mean estimated rate to transform the \code{timeVar} predictor into an exponentially decaying variable
interpolating between initial time (time offset magnitude of 1) and arbitrarily large time values (time
offset magnitude 0). Last uses this transformed time variable in a \code{glm} model
(i.e., attempts to answer the question "how different was the start than the end?").
}
\details{
Rate is parameterized as a time constant, or the amount of time it takes for half of change to occur.
The value of rate has a lower bound of
the .0333 quantile of the time variable (i.e., 87.5\% of change happens in the first 10\% of time) and an upper bound of the
.333 quantile of the time variable (i.e., 87.5\% of change takes 100\% of the time to happen). These bounds provide
some robustness in estimates of asympototic effects (i.e., "controlling for time") as well as initial effects
(i.e., "time-related starting offset").

Mean estimated rate is calculated after trimming the upper 25\% and lower 25\% of bootstrapped rate estimates, for robustness to
extremes in resampling.
}
\note{
Although the time variable is transformed to exponentially decay toward zero, this does not necessarily mean
that the model prediction involves an exponential change with time. The nonlinear change in time relates
to the time-associated model coefficients.

The \code{TEglm} approach to including a nonlinear time function in regression is quite different than the
\code{\link{TEfit}} approach. \code{TEglm} utilizes a point estimate for the rate parameter in order to
coerce the model into a generalized linear format; \code{\link{TEfit}} simultaneously finds the best
combination of rate, start, and asympote parameters. In effect, \code{TEglm} treats \emph{magnitude}
of change as being of theoretical interest, while \code{\link{TEfit}} treats the starting value, rate, and
the asymptotic value as each being of theoretical interest.
}
\examples{
dat <- data.frame(trialNum = 1:200, resp = rbinom(200,1,log(11:210)/log(300)))
m_glm <- TEglm(resp ~ trialNum,dat,'trialNum',family=binomial)
summary(m_glm)
m_glm$rate # estimated half-of-change time constant
summary(m_glm$bootRate) # bootstrapped parameter distributions
cor(m_glm$bootRate) # bootstrapped parameter correlations

}
\seealso{
\code{\link{TEglmem}} for mixed-effects extension of \code{TEglm};
\code{\link{TElm}} for a linear model version of \code{TEglm}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEglmem.R
\name{TEglmem}
\alias{TEglmem}
\title{Generalized linear mixed-effects model with nonlinear time random effects}
\usage{
TEglmem(
  formIn,
  dat,
  timeVar,
  groupingVar,
  family = gaussian,
  startingOffset = T,
  nRuns = 5,
  silent = F
)
}
\arguments{
\item{formIn}{model formula, as in \code{glmer}}

\item{dat}{model data, as in \code{glmer}}

\item{timeVar}{String. Indicates which variable in \code{datIn} corresponds to time (i.e., should be transformed). Must be numeric and positive.}

\item{groupingVar}{String. Indicates which variable in \code{datIn} should have a time-related random effect.}

\item{family}{model family, as in \code{glmer}}

\item{startingOffset}{By default (if T) time is coded to start at 1 and saturate to 0. If startingOffset is F, time starts at 0 and saturates to 1. May assist in interpreting interactions with other variables, etc.}

\item{nRuns}{Number of times to run optimization of the rate (i.e., fitting nonlinear transformations of \code{timeVar})}

\item{silent}{Progress is printed by default. silent=T to suppress}
}
\value{
A list including:
\describe{
\item{\code{glmerMod}}{\code{\link[lme4]{glmer}} model fit with transformed time variable}
\item{\code{rates}}{Named vector of rates [\emph{binary log of 50-percent-of-change time constants}]}
\item{\code{timeDat}}{Data frame with original and transformed time variable}
\item{\code{groupMods}}{List of fit \code{\link{TEglm}} models, and the corresponding transformed time variable and named vector of rates}
}
}
\description{
Fits a \code{\link[lme4]{glmer}} generalized linear mixed-effects model with the random effects of
\code{timeVar} for each level of \code{groupingVar}. Provides estimates of time-related change
(i.e., attempts to answer the question "how different was the start than the end?").
}
\details{
First uses \code{\link{TEglm}} to find a rate parameter for each level of \code{groupingVar}, with
the formula extracted using \code{\link{tef_getRanefForm}}. These
rate parameters are used to transform the corresponding \code{timeVar} into a exponentially
saturating variable (see \code{\link{TEglm}}). After finding an initial set of
rate parameters using \code{\link{TEglm}},
\code{TEglmem} attempts to optimize the vector of rate parameters in conjunction with the full
\code{glmer} model.

May be used, with \code{nRuns=0}, to simply
use rate estimates from independent \code{groupingVar}-level \code{\link{TEglm}} models, extracting the corresponding
transformed time variables and using them in a GLMEM.
}
\note{
Random effects and rate estimates may be unstable, and optimization may take
a very long time to run. The primary purpose of this function is to allow for by-\code{groupingVar}
detrending of time-related changes in data (i.e., to estimate and test fixed effects at asymptotic time,
or to estimate and test the magnitude of time-related effects).
If reliable by-\code{groupingVar} parameters are desired, it is highly recommended to use
\code{\link{TEfit}} or \code{\link{TEfitAll}}.

The \code{formIn} must include a random effect of \code{timeVar} by \code{groupingVar}
(e.g., \code{(time_variable | grouping_variable)}).

Although the time variable is transformed to exponentially decay from one toward zero
(or, if \code{startingOffset=F}, from zero toward one), this does not necessarily mean
that the model prediction involves an exponential change with time. The nonlinear change in time relates
to the time-associated model coefficients.
}
\examples{
\dontrun{
m_TEglmem <- TEglmem(resp ~ ratio + trialNum:ratio + (ratio + trialNum:ratio || subID),anstrain, timeVar = 'trialNum',groupingVar = 'subID',family=binomial)
# Typical glmer model:
summary(m_TEglmem$glmerMod)
# Participant-level rate parameters:
m_TEglmem$rates
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_acc2dprime.R
\name{tef_acc2dprime}
\alias{tef_acc2dprime}
\title{Convert an ordered vector of accuracies to d-prime}
\usage{
tef_acc2dprime(
  accuracy,
  stim_present,
  by_index = T,
  trial_hwhm = 3,
  max_dprime = 5
)
}
\arguments{
\item{accuracy}{At each index, what was the accuracy: [bounded at 0 and 1]}

\item{stim_present}{At each index, was the stimulus present or absent: [binary; 0 and 1, or logical]}

\item{by_index}{Should the d-prime be calculated for each index?}

\item{trial_hwhm}{The Gaussian smoother has a half-width-half-max; values are given half weight at this index distance from the center index (as the smoother iterates through each index in turn as the center). An arbitrarily small HWHM will lead to the same behavior as linear interpolation.}

\item{max_dprime}{d-prime becomes infinite as accuracy approaches 0 or 1. This value limits the absolute value of d-prime.}
}
\description{
Given paired information regarding accuracy and stimulus presence,
run a Gaussian-weighted-mean smoother [kernel] over the accuracy vector separately for
stimulus-present and stimulus-absent indices, then compute an index-wise d-prime.
Returning and entire-vector d-prime (i.e., stable across time) is also an option.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.TEfit.R
\name{plot.TEfit}
\alias{plot.TEfit}
\title{Plot a TEfit}
\usage{
\method{plot}{TEfit}(
  TEs3,
  plot_title = "",
  xlabel = "",
  ylabel = "",
  sub_title = "",
  ymin = NA,
  ymax = NA
)
}
\arguments{
\item{TEs3}{TEfit model output}

\item{plot_title}{optional plot title}

\item{xlabel}{optional plot x axis label}

\item{ylabel}{optional plot y axis label}

\item{sub_title}{optional plot caption}

\item{ymin}{optional lower boundary for Y axis}

\item{ymax}{optional upper boundary for X axis}
}
\description{
Plot predicted values of a TEfit model. Predictions are thresholds, if relevant;
otherwise they are overall model predictions
}
\examples{
\dontrun{
m <- TEfit(anstrain_s1[,c('acc','trialNum')])
plot(m)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_control.R
\name{tef_control}
\alias{tef_control}
\title{Get control parameters for a TEfit model}
\usage{
tef_control(
  quietErrs = F,
  suppressWarnings = F,
  nTries = 200,
  y_lim = c(-1e+07, 1e+07),
  rate_lim = c(0, 0),
  shape_lim = c(0, 0),
  expBase = 2,
  rateBase = 2,
  pFix = c(),
  penalizeMean = T,
  penalizeRate = F,
  convergeTol = 0.05,
  stepwise_asym = F,
  explicit = ""
)
}
\arguments{
\item{quietErrs}{logical. Should errors be printed to the Console?}

\item{suppressWarnings}{logical. Should warnings be printed to the Console?}

\item{nTries}{Numeric. What is the maximum number of optimization runs that should be attempted?}

\item{y_lim}{Numeric vector of length 2. Lower and upper bounds of permitted predicted values.}

\item{rate_lim}{Numeric vector of length 2. Lower and upper bounds of permitted rate values [log time constants].}

\item{shape_lim}{Numeric vector of length 2. Lower and upper bounds of permitted shape parameter values (i.e., for Weibull).}

\item{expBase}{For change functions with an exponential component, what should the base of the exponent be?}

\item{rateBase}{What should the base of the rate exponent be?}

\item{pFix}{Named numeric vector allowing specific parameters to be fixed to a constant (i.e., not estimated)}

\item{penalizeMean}{Logical. Should the time-evolving model be penalized if the mean of the time-evolving predicted values diverges from the mean of the null [non-time-evolving] predicted values?}

\item{penalizeRate}{Logical. Should the time-evolving model be penalized if the rate parameter is very near a boundary?}

\item{convergeTol}{Convergence is extremely roughly defined in \code{TEfits} as the SD of the same estimated parameter on different runs with relatively low error. What should this SD be?}

\item{stepwise_asym}{Logical. If a function will saturate by the end of the measurement time, this option allows the asymptote to be estimated from this time period (i.e., as stationary).}

\item{explicit}{Character. Rather than using any of the pre-defined change or link functions, enter the specific function you want to test.}
}
\description{
\code{\link{TEfit}} internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_penalizedErr.R
\name{tef_penalizedErr}
\alias{tef_penalizedErr}
\title{Penalize error for being close to a boundary}
\usage{
tef_penalizedErr(boundedPar, errPar, loBound, upBound, dbeta_penalty = 1.001)
}
\arguments{
\item{boundedPar}{Value of the parameter that has bounds.}

\item{errPar}{Value of the error to be penalized.}

\item{loBound}{Lower parameter boundary.}

\item{upBound}{Upper parameter boundary.}

\item{dbeta_penalty}{The multiplicative error penalty at approximately 5.3 percent away from the boundary. Error increases with increasing proximity to a bound.}
}
\description{
\code{\link{TEfit}} internal.
Formally, `pErr=Err/dbeta(par,dbeta_penalty,dbeta_penalty)` when Err is positive, where
par is the (often rate) parameter normalized between 0 [lower bound] and 1 [upper bound]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_mono_expo.R
\name{tef_mono_expo}
\alias{tef_mono_expo}
\title{Create an integer vector coding a saturating function of time}
\usage{
tef_mono_expo(N, style = "increasing")
}
\arguments{
\item{N}{Number of time values}

\item{style}{Should the integer "blocks" be saturating (increasing in length) ("increasing") or equal in size ("equal")}
}
\description{
Create an integer vector coding a saturating function of time
}
\examples{
\dontrun{
# # Make a temporary data set
d_tmp <- anstrain_s1

# # Get an integer-coded time variable
d_tmp$trialNum_mono <- tef_mono_expo(nrow(d_tmp))

# # Look at the new coding of time
plot(d_tmp$trialNum,d_tmp$trialNum_mono)

# # Fit a model predicting accuracy using a monotonic measure of time; see documentation for brms::mo
library(brms)
m <- brm(acc ~ mo(trialNum_mono), d_tmp)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_change_power3.R
\name{tef_change_power3}
\alias{tef_change_power3}
\title{Construct a 3-parameter power function of change}
\usage{
tef_change_power3(
  timeVar,
  parForm = ~1,
  startForm = ~1,
  rateForm = ~1,
  asymForm = ~1,
  rateBase = 2,
  propRemain = 0.25
)
}
\arguments{
\item{timeVar}{String. The name of the variable in the model that corresponds to time. The variable of time should be positive and numeric, and the function of change should be expected to happen with increasing time.}

\item{parForm}{The right-hand side of the formula defining all nonlinear parameters as well as the null [non-time-varying] model.}

\item{startForm}{The right-hand side of the formula defining the start parameter. Overwrites \code{parForm} for this parameter.}

\item{rateForm}{The right-hand side of the formula defining the rate parameter. Overwrites \code{parForm} for this parameter.}

\item{asymForm}{The right-hand side of the formula defining the asymptote parameter. Overwrites \code{parForm} for this parameter.}

\item{rateBase}{Number. The base of the log (e.g., 2 or \code{exp(1)}) of the rate [time constant].}

\item{propRemain}{Change rate is parameterized in terms of the \code{rateBase} log of time to this proportion of change remaining (i.e., \code{1-propRemain} of total change occurs in \code{rateBase^rateParameter} time.}
}
\description{
By defining the model variable associated with time (e.g., trial number), and
formulas defining each of the nonlinear parameters of time-related change,
this function constructs a model that can then be passed to functions for fitting
the model (e.g., \code{\link{TEbrm}}).
}
\details{
Function is \strong{under development}
and is likely to be buggy, and to change frequently.
}
\examples{
equation_to_fit <- tef_change_power3('timeVar',rateForm = ~ xvar1*xvar2) # both variables should be numeric for TEfit methods! TEbrm should work with factorsas well

equation_to_fit <- tef_change_power3('timeVar'
, parForm = ~ xvar1   # overall parameter formula; overwritten for time-evolving formulas below, leaving this as just the null model
, startForm = ~ xvar2 # start parameter's regression model
, rateForm = ~ (1|participantID)  # rate parameter's [mixed-effects] regression model
, asymForm = ~ xvar3 # asymptote parameter's regression model
)

}
\seealso{
\code{\link{TEbrm}} for examples of how to use this function in specifying models.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_change_gam.R
\name{tef_change_gam}
\alias{tef_change_gam}
\title{Construct a by-time Generalized Additive Model formula}
\usage{
tef_change_gam(
  timeVar,
  timeCovar = c(NULL),
  groupingVar = "",
  bs = c("ts", "cr", "cs", "tp", "cc")
)
}
\arguments{
\item{timeVar}{String. Indicates the name of the time variable}

\item{timeCovar}{Vector of strings, for variables that will be included in the time-evolving model.}

\item{groupingVar}{String. Optional grouping variable (e.g., participant ID, in a behavioral study). If this is non-empty, then all \code{timeCovar} are fit on the level of \code{groupingVar}}

\item{bs}{Two letter character string inidicating the basis for smoothing. See ?mgcv::s}
}
\description{
By defining the model variable associated with time (e.g., trial number), and
formulas defining each of the nonlinear parameters of time-related change,
this function constructs a model that can then be passed to functions for fitting
the model (e.g., \code{\link{TEbrm}}).
}
\details{
Use is primarily for passing to link functions. If no link function is needed, then an equivalent
formula can be passed directly to \code{\link[brms]{brm}}.

Function is \strong{under development}
and is likely to be buggy, and to change frequently.
}
\examples{
equation_to_fit <- tef_change_gam('trialNum')


}
\seealso{
\code{\link{TEbrm}} for examples of how to use this function in specifying models.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEfit.R
\name{TEfit}
\alias{TEfit}
\title{Fit a time-evolving model (nonlinear regression by minimizing error)}
\usage{
TEfit(
  varIn,
  linkFun = list(link = "identity"),
  errFun = "ols",
  changeFun = "expo",
  bootPars = tef_bootList(),
  blockTimeVar = NULL,
  covarTerms = list(),
  control = tef_control()
)
}
\arguments{
\item{varIn}{Data frame or vector. First column [or vector] must be the time-dependent response variable (left hand side of regression). If available, second column must be the time variable. All other columns are covariates, possibly involved in a link function.}

\item{linkFun}{A list defining a link function (i.e., 'identity', 'd_prime', 'weibull', or 'logistic')}

\item{errFun}{A string defining an error function (e.g., 'ols', 'logcosh', 'bernoulli').}

\item{changeFun}{A string defining the functional form of change (e.g., 'expo', 'power', 'weibull')}

\item{bootPars}{A list defining the details for bootstrapped fits. Defaults to no bootstrapping. Necessary for estimates of uncertainty around fits and for covariance between parameters.}

\item{blockTimeVar}{A string identifying which covariate is the time points of sub-scales (e.g., "blocks" of times within the overall timescale of data collection)}

\item{covarTerms}{An optional list of logical vectors indicating whether parameters should vary by covariates. See examples.}

\item{control}{A list of model parameters. Use of tef_control() is highly recommended.}
}
\value{
A \code{TEfit} S3 object including:
\describe{
\item{\code{model}}{The best fit of parameters to the dataset, along with associated items such as tests of goodness-of-fit and nonstationarity}
\item{\code{nullFit}}{The best fit of a model that does not change over time, but otherwise uses the same parameterization as \code{model}}
\item{\code{data}}{Data frame of the input variables fit by the model, with the additiona of fit values}
\item{\code{modList}}{List of model details}
\item{\code{bootList}}{\emph{(if relevant)} List of bootstrapped \code{TEfit} models.}
}
}
\description{
This is the primary function for the \code{TEfits} package (but see \code{\link{TEbrm}} for
a more powerful approach). Fits a
time-evolving regression model. Many options are available for
various error functions, functional forms of change,
nested timescales, bootstrapping/subsampling/cross-validation, and so on.
Various handy S3 methods are available, such as
\code{plot}, \code{summary},
\code{coef}, and \code{simulate}.
}
\details{
TEfit defines a nonlinear regression model and re-fits that model
using \code{optim} numerous times, with parameter values randomly initialized prior to optimization, until
the highest-likelihood fitting runs also have parameters very similar to
one another (i.e., SD less than the convergence criterion). Runs are
implemented in batches of 10. Convergence is a heuristic and should ideally be
corroborated with other measures (e.g., bootstrapping).

Bootstrapping or subsampling is specified as
\code{bootPars=tef_bootList(resamples = 0, bootPercent = 1, bootTries = 20)}.
\code{resamples} refers to the number of times the model is re-fit on resampled data,
\code{bootPercent} is the proportion (between 0 and 1) of the data resampled, and
\code{bootTries} is the number of optimization runs attempted on each subsample.
\code{bootPercent} of 1, the default, implements resampling with replacement (bootstrapping).
\code{bootPercent} less than 1 implements resampling without replacement, fitting
the model to that subsample, and evaluation of the fit values on the
left-out subsample (i.e., cross-validation). \code{bootTries} defaults to a very small number (20).

Currently supported \strong{error functions} are:
\itemize{
\item{\code{ols}, i.e. \code{sum((y-yHat)^2)} -- sum of squared error}
\item{\code{rmse}, i.e. \code{sqrt(mean((y-yHat)^2))} -- root mean squared error}
\item{\code{logcosh}, i.e. \code{sum(log(cosh(y-yHat)))} -- log-hyperbolic-cosine}
\item{\code{bernoulli}, i.e. \code{-sum(y*log(yHat) + (1-y)*log(1-yHat))} -- Bernoulli [binary binomial]}
\item{\code{exGauss_mu}, i.e. \code{-sum(log(retimes::dexgauss(y,mu=yHat,sigma=sigma_param,tau=tau_param)))} --
ex-Gaussian distribution with time-evolving change in the Gaussian mean parameter}
\item{\code{exGauss_tau}, i.e. \code{-sum(log(retimes::dexgauss(y,mu=mu_param,sigma=sigma_param,tau=yHat)))} --
ex-Gaussian distribution with time-evolving change in the tau parameter}
}

Currently supported \strong{link functions} are:
\itemize{
\item{\code{identity} -- implemented by \code{linkFun=list(link='identity')} --
Default. The predicted values of the time-evolving function are the predicted
values of the model.}

\item{\code{logit} -- implemented by
\code{linkFun=list(link='logit',logistX='variableName',
threshChange=T,biasChange=F,fitThresh=.75,lapseRate=.005)}
-- The predicted values of the time-evolving function are threshold values
of \code{logistX} (by default) and/or
the bias value (defaults to constant). \code{link} and \code{logistX} are required.
Other parameters have default values
(i.e., modelled value of the outcome variable [\code{fitThresh}],
offset of the predicted values to prevent a pathological error calculation [\code{lapseRate}])
}
\item{\code{weibull} -- implemented by
\code{linkFun=list(link='weibull',weibullX='variableName',
fitThresh=.75,yIntercept=.5,rhAsymptote=1,lapseRate=.005)}
-- The predicted values of the time-evolving function are threshold values
of \code{weibullX}.
\code{link} and \code{weibullX} are required. Other parameters
have default values
(i.e., modelled value of the outcome variable [\code{fitThresh}],
value of the outcome variable at weibullX==0 [\code{yIntercept}],
value of the outcome variable at weibullX==Inf [\code{rhAsymptote}],
offset of the predicted values to prevent a pathological error calculation [\code{lapseRate}])
.}
\item{\code{d_prime} -- implemented by
\code{linkFun=list(link='d_prime',presence='variableName',
max_d_prime=5,smooth_hwhm=3)}
-- \code{link} and \code{presence} are required, \code{max_d_prime}
and \code{smooth_hwhm} have default values. The pFA and pH are
first calculated using a windowed average of stimulus-present or
stimulus-absent trials (penalized to bound the max d-prime), then calculating the
by-timepoint d-prime, then fitting that d-prime as the response variable using an identity link. See
\code{\link{tef_acc2dprime}} and \code{\link{tef_runningMean}} for details about the intermediate steps.}
}

Currently supported \strong{change functions} are:
\itemize{
\item{\code{expo} -- 3-parameter exponential (start, [inverse] rate, and asymptote) -- rate is log of time to some proportion remaining, default is log2 of time to 50 percent remaining}
\item{\code{expo_block} -- 3-parameter exponential (start, [inverse] rate, and asymptote)
plus 2-paramter multiplicative changes on timescales that are a subset of the whole}
\item{\code{expo_double} -- 4-parameter exponential (start, two equally weighted [inverse] rates, and asymptote)}
\item{\code{power} -- 3-parameter power (start, [inverse] rate, and asymptote) -- rate is log of time to some proportion remaining, defaulting to log2 of time to 25 percent remaining}
\item{\code{power4} -- 4-parameter power (start, [inverse] rate, asymptote, and "previous learning time")}
\item{\code{weibull} -- 4-parameter weibull (start, [inverse] rate, asymptote, and shape) -- rate is same as \code{expo}}
}
}
\note{
By default, the mean of the time-evolving model's fit values should be very similar to the mean of the null fit values.
This is enforced by penalizing the time-evolving model's error multiplicatively by 1 + the square of the difference
between the average of the model prediction and the average of the null [non-time-evolving] prediction. This is intended to
constrain model predictions to a "sane" range. This constraint can be removed with \code{control=tef_control(penalizeMean=F)}.

Currently \strong{known bugs} are:
\itemize{
\item{The logistic linkfun must include threshChange=T. Returns error if only the bias term is allowed to change.}
}
}
\examples{
\dontrun{
## example data:
dat <- data.frame(timeVar = 1:50, respVar = c(seq(.3,.9,length=25),seq(.9,.91,length=25))+rep(c(0,.01),25),covar1=rep(c(1,2),25),covar2 = rep(c(-3,-1,0,2,5),10))

## Default fitting of 'ols' error function and 3-parameter exponential change function:
m <- TEfit(dat[,c('respVar','timeVar')])
summary(m)
# view a plot of the model:
plot(m)

## 'bernoulli' error function:
m <- TEfit(dat[,c('respVar','timeVar')],errFun='bernoulli')
summary(m)

## 3-parameter power change function:
m <- TEfit(dat[,c('respVar','timeVar')],changeFun='power')
summary(m)
# view a plot of the model:
plot(m)

## logistic threshold change:
m <- TEfit(dat[,c('respVar','timeVar','covar2')],errFun='bernoulli',linkFun=list(link='logit',logistX='covar2'))

## include 2 covariates (on all 3 parameters by default):
m <- TEfit(dat[,c('respVar','timeVar','covar1','covar2')])
summary(m) # (likely does not converge due to too many [nonsense] covariates)
plot(m)

## include 2 covariates:
## asymptote and rate are affected by covar1, start and rate are affected by covar2
m <- TEfit(dat[,c('respVar','timeVar','covar1','covar2')],covarTerms=list(pStart=c(F,T),pRate=c(T,T),pAsym=c(T,F)))

## 50 bootstrapped fits:
 m <- TEfit(dat[,c('respVar','timeVar')],bootPars=tef_bootList(resamples=50))
 summary(m)
# view a plot of the model, with CI bands:
plot(m)
# view the predicted values of the model by plotting data simulated from the parameters:
. <- simulate(m,toPlot=T)

 ## 50 random-subsample 80/20 cross-validation fits:
 m <- TEfit(dat[,c('respVar','timeVar')],bootPars=tef_bootList(resamples=50,bootPercent=.8))
 summary(m)

 ## ## ## control parameters:

 ## Increase convergence tolerance to 0.1:
 m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(convergeTol=.1))

 ## Increase the maximum run number to 5000 (defaults to 200):
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(nTries=5000))

 ## If the function will asymptote in the given time period, then one option is to calculate the TE function stepwise: first get a stable fit of last 20\% of timepoints (if there are enough timepoints, average this with a stable fit to the last 10\% of timepoints). Then fit the start and rate of approach to this asymptote:
 m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(stepwise_asym = T))

 ## Put limits on the predicted values:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(y_lim=c(.1,.9)))

 ## Put limits on the rate parameter:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(rate_lim=c(2,4)))

 ## Remove the constraint that the time-evolving fit values should have the same mean as the null fit values:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(penalizeMean=F))

 ## If rate parameter is hitting the boundary, try imposing a slight penalization for extreme rate values:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(penalizeRate=T))

 ## Change the exponential change log base from 2 to e:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(expBase=exp(1)))

 ## Change the rate parameter log base from 2 to e:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(rateBase=exp(1)))

 ## Silence errors:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(quietErrs=T))

 ## Fix a parameter [asymptote] to 0.8:
  m <- TEfit(dat[,c('respVar','timeVar')],control=tef_control(pFix=list(pAsym=.8)))
}
}
\seealso{
For interpreting model outputs: \code{\link{plot.TEfit}}; \code{\link{summary.TEfit}};
\code{\link{coef.TEfit}}; \code{\link{simulate.TEfit}}

\code{\link{TEfitAll}} for fitting a set of \code{TEfit} models.

\code{\link{TEbrm}} for fitting a Bayesian regression model, with many options
including fixed and random effects.

For including a nonlinear time predictor in [generalized] linear regression frameworks:
\code{\link{TElm}}; \code{\link{TEglm}}; \code{\link{TElmem}}; \code{\link{TEglmem}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.TEfit.R
\name{residuals.TEfit}
\alias{residuals.TEfit}
\title{Extract residuals from a TEfit model}
\usage{
\method{residuals}{TEfit}(modelIn, toPlot = F)
}
\arguments{
\item{modelIn}{A TEfit model}
}
\description{
Deprecated
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TElmem.R
\name{TElmem}
\alias{TElmem}
\title{Linear mixed-effects model with nonlinear time random effects}
\usage{
TElmem(
  formIn,
  dat,
  timeVar,
  groupingVar,
  nRuns = 5,
  startingOffset = T,
  silent = F
)
}
\arguments{
\item{formIn}{model formula, as in \code{lmer}}

\item{dat}{model data, as in \code{lmer}}

\item{timeVar}{String. Indicates which variable in \code{datIn} corresponds to time (i.e., should be transformed). Must be numeric and positive.}

\item{groupingVar}{String. Indicates which variable in \code{datIn} should have a time=related random effect.}

\item{nRuns}{Number of times to run optimization of the rate (i.e., fitting nonlinear transformations of \code{timeVar})}

\item{startingOffset}{By default (if T) time is coded to start at 1 and saturate to 0. If startingOffset is F, time starts at 0 and saturates to 1. May assist in interpreting interactions with other variables, etc.}

\item{silent}{Progress is printed by default. silent=T to suppress}
}
\value{
A list including:
\describe{
\item{\code{lmerMod}}{\code{\link[lme4]{lmer}} model fit with transformed time variable}
\item{\code{rates}}{Named vector of rates [\emph{50-percent-of-change time constants}]}
\item{\code{timeDat}}{Data frame with original and transformed time variable}
\item{\code{groupMods}}{List of fit \code{\link{TElm}} models, and the corresponding transformed time variable and named vector of rates}
}
}
\description{
Fits a \code{\link[lme4]{lmer}} linear mixed-effects model with the random effects of
\code{timeVar} for each level of \code{groupingVar}. Provides estimates of time-related change
(i.e., attempts to answer the question "how different was the start than the end?").
}
\details{
First uses \code{\link{TElm}} to find a rate parameter for each level of \code{groupingVar}, with
the formula extracted using \code{\link{tef_getRanefForm}}. These
rate parameters are used to transform the corresponding \code{timeVar} into a exponentially
saturating variable (see \code{\link{TElm}}). After finding an initial set of
rate parameters using \code{\link{TElm}},
\code{TElmem} attempts to optimize the vector of rate parameters in conjunction with the full
\code{lmer} model.

May be used, with \code{nRuns=0}, to simply
use rate estimates from independent \code{groupingVar}-level \code{\link{TElm}} models, extracting the corresponding
transformed time variables, and using them in a LMEM.
}
\note{
Random effects and rate estimates may be unstable, and optimization may take
a very long time to run. The primary purpose of this function is to allow for by-\code{groupingVar}
detrending of time-related changes in data (i.e., to estimate and test fixed effects at asymptotic time,
or to estimate and test the magnitude of time-related effects).
If reliable by-\code{groupingVar} parameters are desired, especially of rate, it is highly recommended to use
\code{\link{TEfit}} or \code{\link{TEfitAll}}.

The \code{formIn} must include a random effect of \code{timeVar} by \code{groupingVar}
(e.g., \code{(time_variable | grouping_variable)})
}
\examples{
\dontrun{
m_TElmem <- TElmem(acc ~ trialNum + (trialNum || subID), anstrain, timeVar = 'trialNum',groupingVar = 'subID')
# Typical lmer model:
summary(m_TElmem$lmerMod) # On average, starting accuracy was .137 worse than asymptotic accuracy
# Participant-level rate parameters:
m_TElmem$rates
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_windowWDM.R
\name{tef_windowWDM}
\alias{tef_windowWDM}
\title{Find a data frame of Wiener Diffusion Model parameters}
\usage{
tef_windowWDM(dat, windowSize = 20, fit_BS = FALSE)
}
\arguments{
\item{dat}{Response time vector. Positive values indicates upper-boundary while negative values indicate lower-boundary.}

\item{windowSize}{Total width of each window for estimating parameters.}

\item{fit_BS}{Logical. If \code{FALSE} only drift rate is estimated, if \code{TRUE} both drift rate and boundary separation are estimated.}
}
\description{
Estimates a set of Wiener Diffusion Model parameters for a vector of formatted
data. Estimation uses \code{\link[RWiener]{wdm}} and data must be formatted as a
1-dimensional vector according the to \code{RWiener} package (i.e., RT * 2 * (boundary=='upper' - 0.5) ).
Assumes that RTs are associated with equally-spaced and sorted timepoints (e.g., trials of a behavioral study).
First fits a \code{\link[RWiener]{wdm}} for the entire RT vector, and applies some of these parameters
(NDT, bias, and possibly BS) to the entire RT vector. Next takes a set of windows, equally spaced and centered at
\code{windowSize/2}, and estimates the DR (and possibly BS) for each window. Averages overlapping windows.
Returns a data frame of estimated parameters, including 95\% CI.
}
\examples{
\dontrun{
dat <- (rnorm(200,.2,.05) + rexp(200,1 / .2))*sample(c(-1,1),200,replace=T)
    m_drY_bsN <- tef_windowWDM(dat,fit_BS = F)
    m_drY_bsY <- tef_windowWDM(dat,fit_BS = T)
    psych::pairs.panels(cbind(trialNum = 1:200, RT = abs(dat),m_drY_bsY[,1:6]))
    }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TElm.R
\name{TElm}
\alias{TElm}
\title{[Robust] linear model with nonlinear time predictor}
\usage{
TElm(
  formIn,
  dat,
  timeVar,
  startingOffset = T,
  robust = F,
  fixRate = NA,
  nBoot = 250
)
}
\arguments{
\item{formIn}{model formula, as in \code{lm()}}

\item{dat}{model data, as in \code{lm()}}

\item{timeVar}{String. Indicates which model predictor is time (i.e., should be transformed). Must be numeric and positive.}

\item{startingOffset}{By default (if T) time is coded to start at 1 and saturate to 0. If startingOffset is F, time starts at 0 and saturates to 1. May assist in interpreting interactions with other variables, etc.}

\item{robust}{Logical. Should \code{\link[MASS]{rlm}} be used?}

\item{fixRate}{If numeric, use this as a rate parameter [binary-log of 50 percent time constant] rather than estimating it (e.g., to improve reproducibility)}

\item{nBoot}{Number of bootstrapped models to fit after rate [time constant] has been estimated (passed to \code{\link{tef_rlm_boot}})}
}
\description{
Fit a linear model or robust linear model with time as a covariate,
while estimating the shape of the nonlinear interpolation between starting and ending time.
First resamples data with replacement 200 times, and each time estimates the best-shaped curve to interpolate
between initial time-related offset and asymptotic time (i.e., rate at which effect of time saturates at zero).
Then uses the mean estimated rate to transform the \code{timeVar} predictor into an exponentially decaying variable
interpolating between initial time (time offset magnitude of 1) and arbitrarily large time values (time
offset magnitude 0). Last uses this transformed time variable in a \code{\link{tef_rlm_boot}} [\code{rlm} or \code{lm}] model
(i.e., attempts to answer the question "how different was the start than the end?").
}
\details{
Rate is parameterized as a time constant, or the amount of time it takes for half of change to occur.
The value of rate has a lower bound of
the .0333 quantile of the time variable (i.e., 87.5\% of change happens in the first 10\% of time) and an upper bound of the
.333 quantile of the time variable (i.e., 87.5\% of change takes 100\% of the time to happen). These bounds provide
some robustness in estimates of asympototic effects (i.e., "controlling for time") as well as initial effects
(i.e., "time-related starting offset"). Uses this transformed time variable in a \code{\link{tef_rlm_boot}} model to estimate
bootstrapped parameter coefficients and out-of-sample prediction.

Mean estimated rate is calculated after trimming the upper 25\% and lower 25\% of bootstrapped rate estimates, for robustness to
extremes in resampling.
}
\note{
The \code{TElm} approach to including a nonlinear time function in regression is quite different than the
\code{\link{TEfit}} approach. \code{TElm} utilizes a point estimate for the rate parameter in order to
coerce the model into a generalized linear format; \code{\link{TEfit}} simultaneously finds the best
combination of rate, start, and asympote parameters. In effect, \code{TElm} treats \emph{magnitude}
of change as being of theoretical interest (and rate as a nuisance parameter to be controlled for),
while \code{\link{TEfit}} treats the starting value, rate, and
the asymptotic value as each being of theoretical interest.
}
\examples{
dat <- data.frame(trialNum = 1:200, resp = log(11:210)+rnorm(200))

# using a Linear Model
m_lm <- TElm(resp ~ trialNum,dat, 'trialNum')
summary(m_lm)
m_lm$bootSummary
m_lm$rate

# using a Robust Linear Model
m_rlm <- TElm(resp ~ trialNum,dat,'trialNum',robust=TRUE)
summary(m_rlm)
m_rlm$bootSummary
m_rlm$rate

# comparing fits
plot(dat[,c('trialNum','resp')])
lines(dat$trialNum,fitted(m_lm),col='blue')
lines(dat$trialNum,fitted(m_rlm),col='red')

# Examining the bootstrapped rates and other parameters together
summary(m_rlm$bootRate)
cor(m_rlm$bootRate)

}
\seealso{
\code{\link{TElmem}} for mixed-effects extension of \code{TElm};
\code{\link{TEglm}} for generalized extension of \code{TElm}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_err.R
\name{tef_err}
\alias{tef_err}
\title{Calculate error for TEFit::tef_fitErr}
\usage{
tef_err(y, yHat, errFun, curDat = NA)
}
\arguments{
\item{y}{Actual response variable values}

\item{yHat}{Predicted response variabel values}

\item{errFun}{Error function (e.g., OLS, logcosh)}

\item{curDat}{Data input. Only really relevant when calculating a full likelihood function and needing info beyond yHat}
}
\description{
\code{\link{TEfit}} internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.TEfitAll.R
\name{plot.TEfitAll}
\alias{plot.TEfitAll}
\title{Plot a collection of models fit by TEfitAll()}
\usage{
\method{plot}{TEfitAll}(TEs3s, ymin = NA, ymax = NA)
}
\arguments{
\item{TEs3s}{A TEfitAll model object (collection of fit models)}

\item{ymin}{minimum y value}

\item{ymax}{maximum y value}
}
\description{
Plot a collection of models fit by TEfitAll()
}
\examples{
\dontrun{
m <- TEfitAll(anstrain[,c('acc','trialNum')],groupingVar = anstrain$subID)
plot(m)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_link_weibull.R
\name{tef_link_weibull}
\alias{tef_link_weibull}
\title{Construct a weibull link function parameterized with threshold and bias}
\usage{
tef_link_weibull(
  changeForm,
  linkX,
  threshVal = 0.75,
  rhAsymptote = 1,
  yIntercept = 0.5,
  lapseRate = 0.005,
  boundScale = 2,
  constantPar_prior = "normal(0,3)"
)
}
\arguments{
\item{changeForm}{The formula describing the change in threshold}

\item{linkX}{Character. The name of the "x" variable in the weibull link function (e.g., stimulus strength in a psychometric function). Should be a positive real numeric variable (e.g., presentation time or number of targets).}

\item{threshVal}{The threshold at which to evaluate the weibull function (i.e., the y-value for which threshold describes the x-value).}

\item{rhAsymptote}{The asymptotic value of the weibull function with large \code{linkX} values (e.g., accuracy at infinitely large stimulus strength).}

\item{yIntercept}{The origin value of the weibull function (with a \code{linkX} value of zero, e.g., accuracy at a stimulus strength of zero; in behavioral data is likely to be "guessing rate").}

\item{lapseRate}{The offset, from rhAsymptote, of the weibull function at arbitrarily large values of \code{linkX}. A small lapse rate improves model fit (see Wichmann and Hill, 2001, P&P).}

\item{boundScale}{Currently not implemented. Upper threshold of threshold estimates, as a multiple of the maximum absolute \code{linkX}.}

\item{constantPar_prior}{The prior to put on the constant component of the weibull function (i.e., shape parameter). Only relevant if passing the model to a \code{\link[brms]{brm}} model (e.g., with \code{\link{TEbrm}}).}
}
\description{
Function is \strong{under development}
and is likely to be buggy, and to change frequently.
}
\details{
\strong{\code{shape}} is a parameter in the weibull function, and \strong{must not} be
used as a name of data variables (e.g., \code{linkX} or within \code{changeForm}).
}
\examples{

equation_to_fit <- tef_link_weibull( tef_change_expo3('trialNum', parForm = ~ (1|subID)) , linkX = 'absoluteRatio' )

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_parseParFormula.R
\name{tef_parseParFormula}
\alias{tef_parseParFormula}
\title{Parse a parameter-level formula for a TEfit model}
\usage{
tef_parseParFormula(parFormula, label = "pPar")
}
\arguments{
\item{parFormula}{Formula for fitting a parameter. Variables should be numeric; non-numeric variables may break. Can also take a string or a scalar, if the parameter is to be fixed to that value.}

\item{label}{Label for the parameters}
}
\description{
Returns the \code{parFormula}, with some additional information that assists TEfits in fitting.
}
\details{
If there are no random-effects terms, parses into a literal equation (e.g., turns 
\code{~ x1 + x2 * x3} into \code{pPar_Intercept + pPar_x1\*x1 + pPar_x2\*x2 + pPar_x3\*x3 + pPar_x2_x3\*x2\*x3}) 
which is returned as the "equation" attribute that has an accompanying "parameters" attribute. Additional attributes
include whether there are random-effects terms (attribute "MEM") and whether there are any parameters at all, or if the
"equation" is actually fixed to a variable or scalar (attributed "is_fixed").
}
\examples{
tef_parseParFormula(~ x1 + x2*x3 + (x2 || subID))
tef_parseParFormula(~ x1 + x2*x3, label = 'thresholdAsymptote')
tef_parseParFormula(2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_bootList.R
\name{tef_bootList}
\alias{tef_bootList}
\title{Make a list of parameters for resampling a TEfit model}
\usage{
tef_bootList(resamples = 0, bootPercent = 1, bootTries = 20)
}
\arguments{
\item{resamples}{Number of resamples. More (e.g., >2000) is better, but requires more time and memory.}

\item{bootPercent}{Proportion of data to resample. Must be 0 < bootPercent <= 1. If 1, data is resampled repeatedly with replacement. If less than 1, that proportion of data is repeatedly randomly selected, the model is fit on the selected data, and fit indices (e.g., error) are calculated for the model predictions on the left-out data (i.e., cross-validation).}

\item{bootTries}{Number of optimization runs for each resample}
}
\description{
Used to assist the user in the specification of a set of parameters for resampling analyses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_link_logistic.R
\name{tef_link_logistic}
\alias{tef_link_logistic}
\title{Construct a logistic link function parameterized with threshold and bias}
\usage{
tef_link_logistic(
  changeForm,
  linkX,
  changePar = c("threshold", "bias"),
  threshVal = 0.75,
  lapseRate = 0.005,
  boundScale = 2,
  constantPar_prior = "normal(0,3)"
)
}
\arguments{
\item{changeForm}{The formula describing the change in either threshold or bias}

\item{linkX}{Character. The name of the "x" variable in the logistic link function (e.g., stimulus strength in a psychometric function)}

\item{changePar}{Character. Which variable, "threshold" or "bias", changes over time. The other one is stable over time; \emph{the stable component inherits its formula from the asymptote parameter of the changing component}.}

\item{threshVal}{The threshold at which to evaluate the logistic function (i.e., the y-value for which threshold describes the x-value).}

\item{lapseRate}{The offset, from 0 or 1, of the logistic function at arbitrarily large (positive or negative) values of \code{linkX}. A small lapse rate improves model fit (see Wichmann and Hill, 2001, P&P).}

\item{boundScale}{Currently not implemented. Upper threshold of threshold estimates, as a multiple of the maximum absolute \code{linkX}.}

\item{constantPar_prior}{The prior to put on the constant component of the logistic function (i.e., either the bias or the log threshold). Only relevant if passing the model to a \code{\link[brms]{brm}} model (e.g., with \code{\link{TEbrm}}).}
}
\description{
Function is \strong{under development}
and is likely to be buggy, and to change frequently.
}
\details{
If bias is changing, threshold inherits its formula from asymptotic bias.
If threshold is changing, bias inherits its formula from asymptotic threshold.
}
\examples{

equation_to_fit <- tef_link_logistic( tef_change_expo3('trialNum', parForm = ~ (1|subID)) , linkX = 'ratio' )

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.TEfit.R
\name{summary.TEfit}
\alias{summary.TEfit}
\title{Summarize a time-evolving fit}
\usage{
\method{summary}{TEfit}(TEs3, printOutput = T)
}
\arguments{
\item{TEs3}{A fit TEfit model}

\item{printOutput}{Print output to console (if T) or return a list of summary items (if F)}
}
\description{
Prints or returns a summary of a TEfit object. This includes parameter values,
convergence, the full formula, goodness-of-fit metrics, measures of change in
conditional independence, and (if applicable) distributional information from resampling.
}
\details{
Pseudo-SE is an approximation to the standard error of the parameter using bootstrapped estimates. It is calculated by
first calculating the .025 and .975 quantile() of the resampled parameters. Then, the absolute
difference is calculated between each of these CI values and the overall fit value. The absolute
differences are divided by qnorm(.975) and averaged in order to get the pseudo-SE, the
pseudo standard deviation of the parameter expected value.
}
\examples{
\dontrun{
m <- TEfit(anstrain_s1[,c('acc','trialNum')])
summary(m)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.TEfit.R
\name{simulate.TEfit}
\alias{simulate.TEfit}
\title{Simulate from a TEfit, assuming multivariate gaussian parameters}
\usage{
\method{simulate}{TEfit}(object, nsim = 100, seed = NULL, newdata = data.frame(), toPlot = F)
}
\arguments{
\item{object}{A TEfit model}

\item{nsim}{number of simulations to generate}

\item{seed}{Null, to retain symmetry with other simulate() methods}

\item{newdata}{If desired, new data from which to generate values}

\item{toPlot}{Logical, for plotting of values.}
}
\description{
Requires the $psych$ and $MASS$ packages.
}
\examples{
\dontrun{simulate(model_fit_by_TEfit)}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_rlm_boot.R
\name{tef_rlm_boot}
\alias{tef_rlm_boot}
\title{Bootstrapped robust linear model}
\usage{
tef_rlm_boot(formIn, datIn, nBoot = 500, useLM = F)
}
\arguments{
\item{formIn}{Model formula, as with \code{lm()} or \code{rlm()}}

\item{datIn}{Data, as with \code{lm()} or \code{rlm()}}

\item{nBoot}{Number of resamples [with replacement]}

\item{useLM}{Override the standard \code{rlm()} implementation to use basic \code{lm()} instead}
}
\value{
An augmented \code{rlm} or \code{lm} object that includes
several new items: \code{$bootSummary}, \code{$boots} (all bootstrapped parameters),
\code{$bootQs} (quantiles
of bootstrapped parameters), \code{$dRsq} (all out-of-sample proportional reduction of error),
\code{$dRsqQs} (quantiles of out-of-sample proportional reduction of error),
and \code{$results} (strings, formatted
for RMarkdown, including whole-sample slope, bootstrapped CI, and median out-of-sample dRSq).

For an explanation specific objects see \code{comment(model$boots)},
\code{comment(model$dRsq)}, or \code{comment(model$bootSummary)}.
}
\description{
Run a \code{\link[MASS]{rlm}} model \code{nBoot} times, then include the bootstrapped parameter estimates
(and summary of quantiles thereof) in the rlm output. Also includes out-of-sample
delta-R-squared. If parameter distributions are *anywhere near* decision thresholds,
using \code{nBoot}>2000 (or even much higher) is recommended.
}
\details{
Wraps \code{rlm()} or \code{lm()} and then bootstraps [fits models to datasets sampled with replacement]
that model \code{nBoot} times. Then fits models to
nBoot random 80 percent of data and tests the delta R-squared of each numeric or logical parameter
when predicting the out-of-sample 20 percent.
}
\examples{
dat <- data.frame(x=rnorm(50))
dat$y <- dat$x + rnorm(50)
dat$z <- dat$y - rnorm(50)
dat$z[2] <- NA # the function is robust to NAs
m <- tef_rlm_boot(y~x*z,dat)
m$bootSummary # to get a summary of the fit[s]
comment(m$bootSummary) # to get an explanation of the summary data frame

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_change_weibull.R
\name{tef_change_weibull}
\alias{tef_change_weibull}
\title{Construct a 4-parameter Weibull function of change}
\usage{
tef_change_weibull(
  timeVar,
  parForm = ~1,
  startForm = ~1,
  rateForm = ~1,
  asymForm = ~1,
  shapeForm = ~1,
  changeBase = 2,
  rateBase = 2
)
}
\arguments{
\item{timeVar}{String. The name of the variable in the model that corresponds to time. The variable of time should be positive and numeric, and the function of change should be expected to happen with increasing time.}

\item{parForm}{The right-hand side of the formula defining all nonlinear parameters as well as the null [non-time-varying] model.}

\item{startForm}{The right-hand side of the formula defining the start parameter. If anything besides \code{~1}, overwrites \code{parForm} for this parameter.}

\item{rateForm}{The right-hand side of the formula defining the rate parameter. If anything besides \code{~1}, overwrites \code{parForm} for this parameter.}

\item{asymForm}{The right-hand side of the formula defining the asymptote parameter. If anything besides \code{~1}, overwrites \code{parForm} for this parameter.}

\item{shapeForm}{The right-hand side of the formula defining the weibull shape (i.e., acceleration or deceleration) parameter. If anything besides \code{~1}, overwrites \code{parForm} for this parameter. Shape is estimated on a log scale. Following this, a shape parameter equal to 0 makes the weibull function equivalent to the 3-parameter exponential function.}

\item{changeBase}{Scalar. The base of the log (e.g., 2 or \code{exp(1)}) of the change function.}

\item{rateBase}{Scalar. The base of the log (e.g., 2 or \code{exp(1)}) of the rate [time constant].}
}
\description{
By defining the model variable associated with time (e.g., trial number), and
formulas defining each of the nonlinear parameters of time-related change,
this function constructs a model that can then be passed to functions for fitting
the model (e.g., \code{\link{TEbrm}}).
}
\details{
Function is \strong{under development}
and is likely to be buggy, and to change frequently.
}
\examples{
equation_to_fit <- tef_change_weibull('timeVar',parForm = ~ xvar1*xvar2) # both variables should be numeric for TEfit methods! TEbrm should work with factors as well

equation_to_fit <- tef_change_weibull('timeVar'
, parForm = ~ xvar1   # overall parameter formula; overwritten for time-evolving formulas below, leaving this as just the null model
, startForm = ~ xvar2 # start parameter's regression model
, rateForm = ~ (1|participantID)  # rate parameter's [mixed-effects] regression model
, asymForm = ~ xvar3 # asymptote parameter's regression model
, shapeForm = ~ (1|participantID)  # shape parameter's [mixed-effects] regression model
)

}
\seealso{
\code{\link{TEbrm}} for examples of how to use this function in specifying models.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_getRanefForm.R
\name{tef_getRanefForm}
\alias{tef_getRanefForm}
\title{Extract group-level regressions from a multilevel model formula}
\usage{
tef_getRanefForm(formIn, groupingVar)
}
\arguments{
\item{formIn}{multilevel model formula}

\item{groupingVar}{grouping variable that will be the left-hand-side of the returned formula}
}
\description{
Returns the groupingVar-level regression formulas. For example, the
formula \code{yVar ~ aVar \* bVar \* cVar + (aVar \* cVar | groupingVar)}
will return the group-level formula \code{groupingVar ~ aVar \* cVar}. Note that
All random effects of \code{groupingVar} must be within a single set of (parentheses).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tef_change_expo3.R
\name{tef_change_expo3}
\alias{tef_change_expo3}
\title{Construct a 3-parameter exponential function of change}
\usage{
tef_change_expo3(
  timeVar,
  parForm = ~1,
  startForm = ~1,
  rateForm = ~1,
  asymForm = ~1,
  changeBase = 2,
  rateBase = 2
)
}
\arguments{
\item{timeVar}{String. The name of the variable in the model that corresponds to time. The variable of time should be positive and numeric, and the function of change should be expected to happen with increasing time. \emph{Extensive testing has only been done on positive integer time variables with a minimum time of 1.}}

\item{parForm}{The right-hand side of the formula defining all nonlinear parameters as well as the null [non-time-varying] model.}

\item{startForm}{The right-hand side of the formula defining the start parameter. If anything besides \code{~1}, overwrites \code{parForm} for this parameter.}

\item{rateForm}{The right-hand side of the formula defining the rate parameter. If anything besides \code{~1}, overwrites \code{parForm} for this parameter.}

\item{asymForm}{The right-hand side of the formula defining the asymptote parameter. If anything besides \code{~1}, overwrites \code{parForm} for this parameter.}

\item{changeBase}{Number. The base of the log (e.g., 2 or exp(1)) of the change function.}

\item{rateBase}{Number. The base of the log (e.g., 2 or exp(1)) of the rate [time constant].}
}
\description{
By defining the model variable associated with time (e.g., trial number), and
formulas defining each of the nonlinear parameters of time-related change,
this function constructs a model that can then be passed to functions for fitting
the model (e.g., \code{\link{TEbrm}}).
}
\details{
Function is \strong{under development}
and is likely to be buggy, and to change frequently.
}
\examples{
equation_to_fit <- tef_change_expo3('timeVar',rateForm = ~ xvar1*xvar2) # both variables should be numeric for TEfit methods! TEbrm should work with factorsas well

equation_to_fit <- tef_change_expo3('timeVar'
, parForm = ~ xvar1   # overall parameter formula; overwritten for time-evolving formulas below, leaving this as just the null model
, startForm = ~ xvar2 # start parameter's regression model
, rateForm = ~ (1|participantID)  # rate parameter's [mixed-effects] regression model
, asymForm = ~ xvar3 # asymptote parameter's regression model
)

}
\seealso{
\code{\link{TEbrm}} for examples of how to use this function in specifying models.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{anstrain_s1}
\alias{anstrain_s1}
\title{Perceptual learning data from one participant}
\format{
An object of class \code{data.frame} with 250 rows and 6 columns.
}
\usage{
anstrain_s1
}
\description{
Subset of data from Cochrane, Cui, Hubbard & Green (2019; DOI: 10.3758/s13414-018-01636).
Data from 1 participant completing 250 trials of perceptual learning in an
Approximate Number System ratio discrimination task. "resp" is participant response to the
"ratio" offset of two interspersed colors of dots, with correctness recorded as "acc".
"sizeRat" is the ratio of sizes of the different colors of dots.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TEfits-package.R
\docType{package}
\name{TEfits}
\alias{TEfits}
\title{TEfits: Time-evolving model fits}
\description{
Data is described, interpreted, and tested using indices such as d prime,
mean, or psychometric function threshold. This package serves to allow the same
questions to be asked about time-evolving aspects of these indices, namely,
the starting level, the amount of time that the index takes to change, and the
asymptotic level of that index. Nonlinear regression applied to time-evolving
functions is made as intuitive and painless as is feasible, with many
extensions if desired.
}
\section{Time-evolving fits}{


The \code{\link{TEfit}} function is the primary user-oriented function of
the \code{TEfits} package. It allows for nonlinear fitting of time-related
change in an outcome variable by estimating that outcome variable's value at
each timepoint. See \code{\link{TEfit}}, \code{\link{TEfitAll}} for fitting
a \code{\link{TEfit}} model to subsets of data (e.g., individual participants),
and \code{vignette('TEfits_tutorial')} for an introduction to the framework.
}

\section{Nonlinear regressors in [generalized] linear models}{


While \code{\link{TEfit}} is intended to interrogate time-evolving trends
themselves within data, \code{TEfits} also includes several extensions to common
regression functions that allow for seamless incorporation of a nonlinear
[exponentially saturating] variable of time. These functions approach time-evolving
dynamics from a common perspective in behavioral research: Behavior "of interest"
occurs after transient initial bias in performance (e.g., initial task learning needs to
occur before the target behavior can be effectively measured). The following functions use a
stepwise method to first estimate the rate of change in the outcome variable, then use these rates to
transform the variable of time into a saturating interpolation between 1 [starting offset] and 0
[asymptotic level], and finally include the time variable in the corresponding
[\code{g}]\code{lm}[\code{er}] model:

\itemize{
\item \code{\link{TElm}} wraps this method for \code{\link[stats]{lm}}
\item \code{\link{TEglm}} wraps this method for \code{\link[stats]{glm}}
\item \code{\link{TElmem}} wraps this method for \code{\link[lme4]{lmer}}
\item \code{\link{TEglmem}} wraps this method for \code{\link[lme4]{glmer}}
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{anstrain}
\alias{anstrain}
\title{Perceptual learning data}
\format{
An object of class \code{data.frame} with 1000 rows and 6 columns.
}
\usage{
anstrain
}
\description{
Subset of data from Cochrane, Cui, Hubbard & Green (2019; DOI: 10.3758/s13414-018-01636).
Data from 4 participants ("subID") completing 250 trials each of perceptual learning in an
Approximate Number System ratio discrimination task. "resp" is participant response to the
"ratio" offset of two interspersed colors of dots, with correctness recorded as "acc".
"sizeRat" is the ratio of sizes of the different colors of dots.
}
\keyword{datasets}
