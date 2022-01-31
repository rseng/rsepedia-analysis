
<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable state and is being activelydeveloped](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test coverage](https://codecov.io/gh/ejikeugba/serp/branch/master/graph/badge.svg)](https://codecov.io/gh/ejikeugba/serp?branch=master)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/serp)](https://CRAN.R-project.org/package=serp)
[![CRAN status](https://www.r-pkg.org/badges/version/serp )](https://CRAN.R-project.org/package=serp)
[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ejikeugba/serp?branch=master&svg=true)](https://ci.appveyor.com/project/ejikeugba/serp)
[![R build status](https://github.com/ejikeugba/serp/workflows/R-CMD-check/badge.svg)](https://github.com/ejikeugba/serp/actions)
[![status](https://joss.theoj.org/papers/6ebd3b75ea792be908f0dadebd7cf81c/status.svg)](https://joss.theoj.org/papers/6ebd3b75ea792be908f0dadebd7cf81c)
<!-- badges: end -->

# Contributing to serp

Thank you for taking the time to contribute to the development of `serp`. You could find the following guidelines useful in making your contributions. 

Before you start:

* It is important to have a valid [GitHub account](https://github.com/signup/free).
* Trivial changes to comments or documentation do not require creating a new issue.

## Did you find a bug?

* Make sure the bug was not already reported in the Github [Issues](https://github.com/ejikeugba/serp/issues).
* [Open an issue](https://github.com/ejikeugba/serp/issues/new) and clearly describe the issue with as much information as possible. A code sample or an executable test case are recommended.
  
## Did you plan to write a patch that fixes a bug?

  * [Open an issue](https://github.com/ejikeugba/serp/issues/new) and clearly describes the problem and discuss how your solution will affect `serp`.
  * Fork the repository on GitHub to work on the patch.
  * Get in touch with the maintainer to refine and prioritize your issue.

## Making changes and Pull requests

* Start your work on your fork of the repository. If you haven't done this before, try using `usethis::create_from_github("ejikeugba/serp", fork = TRUE)`.
* Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
* Create a Git branch for your pull request (PR). You may want to use `usethis::pr_init("brief-description-of-change")`.
* Check for unnecessary whitespace with `git diff --check` and format code.
* Commit messages should be descriptive, mentioning what was changed and why, and also **reference the relevant issue number**. 
* Ensure to add the necessary tests for your changes (testthat preferably).
* Run **all** the tests to assure nothing else was accidentally broken, also keep an eye on the test coverage metric.
* Commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser
    
## Copyright issues
 
* On submission, it is crucial your PR includes the following statement: You own the copyright on the code being contributed, and you hereby grant `serp` repo cph unlimited license to use this code in this version or any future version of `serp`. You reserve all other rights to the code.
* It may not be advisable to contribute third party codes to this project. Useful suggestions are nonetheless welcomed.
* The Pull Requests are thereafter reviewed, with feedbacks communicated as soon as possible.

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

# Additional Resources

* [General GitHub documentation](http://help.github.com/)
* [About pull requests](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# serp <a href="https://ejikeugba.github.io/serp/"><img src='man/figures/hex_logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being
activelydeveloped](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/ejikeugba/serp/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ejikeugba/serp?branch=master)
[![Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/serp)](https://CRAN.R-project.org/package=serp)
[![CRAN
status](https://www.r-pkg.org/badges/version/serp)](https://CRAN.R-project.org/package=serp)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03705/status.svg)](https://doi.org/10.21105/joss.03705)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ejikeugba/serp?branch=master&svg=true)](https://ci.appveyor.com/project/ejikeugba/serp)
[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![R build
status](https://github.com/ejikeugba/serp/workflows/R-CMD-check/badge.svg)](https://github.com/ejikeugba/serp/actions)

<!-- badges: end -->

### Overview

The `serp` R package fits cumulative link models (CLMs) with the
`smooth-effect-on-response penalty (SERP)`. The `cumulative model`
developed by McCullagh (1980) is probably the most frequently used
ordinal model in empirical studies. However, the stochastic ordering
property of the general form of the model poses a very serious challenge
in most empirical applications of the model. For instance, unstable
likelihoods with ill-conditioned parameter space are frequently
encountered during the iterative process. `serp` implements a unique
regularization method for CLMs that provides the means of smoothing the
adjacent categories in the model. At extreme shrinkage, SERP causes all
subject-specific effects associated with each variable in the model to
shrink towards unique global effects. Fitting is done using a modified
Newton’s method. Several standard model performance and descriptive
methods are also available. See [Ugba,
2021](https://doi.org/10.21105/joss.03705), [Ugba et al.,
2021](https://doi.org/10.3390/stats4030037) and [Tutz and Gertheiss,
2016](https://doi.org/10.1177/1471082X16642560) for further details on
the implemented penalty.

### Example

Consider the cumulative logit model of the [wine
dataset](https://ejikeugba.github.io/serp/reference/wine.html), where
the rating of wine bitterness is predicted with the two treatment
factors, temperature and contact.

``` r
## The unpenalized non-proportional odds model returns unbounded estimates, hence,
## not fully identifiable.
f1 <- serp(rating ~ temp + contact, slope = "unparallel",
           reverse = TRUE, link = "logit", data = wine)
coef(f1)
```

``` r
## The penalized non-proportional odds model with a user-supplied lambda gives 
## a fully identified model having bounded estimates. A suitable tuning criterion
## could as well be used to select lambda (e.g., aic or cv) 
f2 <- serp(rating ~ temp + contact, slope = "penalize",
           link = "logit", reverse = TRUE, tuneMethod = "user",
           lambda = 1e1 ,data = wine)
coef(f2)
```

``` r
## A penalized partial proportional odds model with one variable set to 
## global effect is also possible.
f3 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "user",
           lambda = 2e1, globalEff = ~ temp, data = wine)
coef(f3)
```

``` r
## The unpenalized proportional odds model with constrained estimates. 
## Under estreme shrinkage, estimates in f2 equal those in this model.  
f4 <-  serp(rating ~ temp + contact, slope = "parallel",
            reverse = FALSE, link = "logit", data = wine)
summary(f4)
```

### Installation and Use

Before installing `serp`, it is encouraged to have a recent version of
[R](https://cran.r-project.org/bin/windows/base/) installed. The
released version of `serp` can be installed from
[CRAN](https://cran.r-project.org/package=serp) with:

``` r
install.packages("serp")
```

or the development version from
[GitHub](https://github.com/ejikeugba/serp) with:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ejikeugba/serp")
```

Load `serp` into R environment with:

``` r
library(serp)
```

### Community Guidelines

Pull requests are welcomed! Please submit your contributions to `serp`
through the list of `Pull Requests`, following the [contributing
guidelines](https://ejikeugba.github.io/serp/CONTRIBUTING.html). To
report issues and/or seek support, please file a new ticket in the
[issue](https://github.com/ejikeugba/serp/issues) tracker, and expect a
feedback ASAP!

### Code of Conduct

Please note that `serp` is released with a [Contributor Code of
Conduct](https://github.com/ejikeugba/serp/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

### References

McCullagh, P. (1980). Regression Models for Ordinal Data. *Journal of
the Royal Statistical Society. Series B (Methodological)*, 42, 109-142.
<https://doi.org/10.1111/j.2517-6161.1980.tb01109.x>

Randall, J (1989). The analysis of sensory data by generalized linear
model. *Biometrical Journal*, 31, 781–793.
<https://doi.org/10.1002/bimj.4710310703>

Tutz, G. and Gertheiss, J. (2016). Regularized Regression for
Categorical Data (With Discussion and Rejoinder). *Statistical
Modelling*, 16, 161-260. <https://doi.org/10.1177/1471082X16642560>

Ugba, E. R., Mörlein, D. and Gertheiss, J. (2021). Smoothing in Ordinal
Regression: An Application to Sensory Data. *Stats*, 4, 616–633.
<https://doi.org/10.3390/stats4030037>

Ugba, E. R. (2021). serp: An R package for smoothing in ordinal
regression *Journal of Open Source Software*, 6(66), 3705.
<https://doi.org/10.21105/joss.03705>
## serp 0.2.3
- CRAN release
- import from crayon, with colored outputs in returned objects 
- changes made in serp test for improved test coverage report
- bug fix in errorMetrics, with model argument also dropped 
- add print method and class of objects returned by errorMetrics
- message() replaces cat() where appropriate
- add citation for serp package

---
## serp 0.2.2
- JOSS release
- minor changes in serp documentation
- provide coefficient() as an alias for coef()
- update README.md to include community guidelines and contributors code of conduct
- 

---
## serp 0.2.1
- submission to CRAN

---
## serp 0.2.0.9001
- update function included in namespace
- examples included in the different function documentation
- shrinkage parameter upper limit set to 1e10 in serp.control
- bug fix in test function

---
## serp 0.2.0.9000 
- deviance tuning option in serp tuneMethod now replaced by AIC 
- serp output includes residual degrees of freedom (rdf)
- function to compute the effective degrees of freedom and rdf from the trace of the generalized hat matrix provided
- serp.summary documentation has its value segment edited
- changes made in serp test functions
- changes made in README.Rmd and README.md

---
## serp 0.2.0
- serp version 0.2.0 release

---
## serp 0.1.9.9001
- re-submission to CRAN

---
## serp 0.1.9.9000
- fixed all issues spotted out in the last version released on CRAN.  

---
## serp 0.1.9
- Submit serp version 0.1.9 to CRAN

---
## serp 0.1.8.9007
- updates README.md
- description gets additional doi
- references in serp documentation updated

---
## serp 0.1.8.9006
- Bugs in tests fixed

---
## serp 0.1.8.9005
- tests reconstructed to yield improved unit test coverage.
- rewrote some of the error and warning messages in key serp functions.
- Bugs in anova.serp fixed.
- 

---
## serp 0.1.8.9004
- update License in the description 
- Bugs in the example codes fixed
- Bugs in test files corrected

---
## serp 0.1.8.9003
- Column titles of predicted values with reverse TRUE are corrected.
- Bugs in the deviance (dvfun) are corrected.
- serp predict function now handles both single and multiple row input(s).
- serp.fit iteration algorithm improved.
- commented lines in serp main function removed.
- reverse and linkf arguments moved from main function to serpfit.
- loglog and cloglog links now give correct results with reverse='TRUE'.
- long lines of codes in the score function split for easy readability.
- starting values augmented to reflect the different link functions.
- reverse statement removed from cv function.
- reverse.fun is removed from summary.serp function.
- predicted values in errorMetrics get normalized for greater efficiency.

---
## serp 0.1.8.9002
- warning messages in serpfit get updated.
- Bug in serpfit when specifying gridType corrected.
- trainError no longer shows up in serp output.
- reverse.fun is now deprecated.
- penalty.print gets slightly reconstructed.
- Bugs in print.summary.serp fixed.
- Summary.serp drops off trainError with some few more adjustments.
- Predict gets reconstructed, also bugs when supplying 'newdata' corrected.
- Bugs in anova function fixed with slope type now included in the output.
- anova test for the 'penalize' slope currently disabled.
- updates in serp.control warning messages.

---
## serp 0.1.8.9001
- serp gets a new license.

---
## serp 0.1.8.9000
* README.md gets additional badges and a hexagon sticker

---
## serp 0.1.8
* package re-submission to CRAN

- reference to methods used in the package is included in the  description field
- values were added to all exported methods with corresponding explanations
- examples now runs by having all \dontrun{} removed
- all occurrences of <<- in the functions are dropped
- READmE.md is updated
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
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
ejike.ugba@outlook.com.
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
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
## What has changed in the current version 

## serp 0.2.3
* import from crayon, with colored outputs in returned objects 
* minor changes in serp documentation
* provide coefficient() as an alias for coef()
* README.md updated to include community guidelines and contributors code of conduct
* changes made in serp test for improved test coverage report
* bug fix in errorMetrics, with model argument also dropped 
* add print method and class of objects returned by errorMetrics
* message() replaces cat() where appropriate
* add citation for serp package

## local R CMD check results

### Test environments
* local x86_64-w64-mingw32 (64-bit), 4.1.1 (2021-08-10) ;

-- R CMD check results ----------------------------------------- serp 0.2.3 ----
Duration: 3m 25s

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded



## R-hub builder

### Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

R CMD check results
0 errors √ | 0 warnings √ | 0 notes √

Status: OK


## win-builder

### Test environments
- using platform: x86_64-w64-mingw32 (64-bit), R 4.0.3 (2020-10-10)
- using platform: x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2021-11-06 r81149)

R CMD check results

Status: OK


## Downstream dependencies
- There are currently no downstream dependencies for this package
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```


# serp <a href="https://ejikeugba.github.io/serp/"><img src='man/figures/hex_logo.png' align="right" height="139" /></a>


<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being activelydeveloped](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test coverage](https://codecov.io/gh/ejikeugba/serp/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ejikeugba/serp?branch=master)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/serp)](https://CRAN.R-project.org/package=serp)
[![CRAN status](https://www.r-pkg.org/badges/version/serp )](https://CRAN.R-project.org/package=serp)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03705/status.svg)](https://doi.org/10.21105/joss.03705)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ejikeugba/serp?branch=master&svg=true)](https://ci.appveyor.com/project/ejikeugba/serp)
[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![R build status](https://github.com/ejikeugba/serp/workflows/R-CMD-check/badge.svg)](https://github.com/ejikeugba/serp/actions)

<!-- badges: end -->


### Overview
The `serp` R package fits cumulative link models (CLMs) with the `smooth-effect-on-response penalty (SERP)`. The `cumulative model` developed by McCullagh (1980) is probably the most frequently used ordinal model in empirical studies. However, the stochastic ordering property of the general form of the model poses a very serious challenge in most empirical applications of the model. For instance, unstable likelihoods with ill-conditioned parameter space are frequently encountered during the iterative process. `serp` implements a unique regularization method for CLMs that provides the means of smoothing the adjacent categories in the model. At extreme shrinkage, SERP causes all subject-specific effects associated with each variable in the model to shrink towards unique global effects. Fitting is done using a modified Newton's method. Several standard model performance and descriptive methods are also available. See [Ugba, 2021](https://doi.org/10.21105/joss.03705), [Ugba et al., 2021](https://doi.org/10.3390/stats4030037) and [Tutz and Gertheiss, 2016](https://doi.org/10.1177/1471082X16642560) for further details on the implemented penalty.


### Example
Consider the cumulative logit model of the [wine dataset](https://ejikeugba.github.io/serp/reference/wine.html), where the rating of wine bitterness is predicted with the two treatment factors, temperature and contact.

```R
## The unpenalized non-proportional odds model returns unbounded estimates, hence,
## not fully identifiable.
f1 <- serp(rating ~ temp + contact, slope = "unparallel",
           reverse = TRUE, link = "logit", data = wine)
coef(f1)
```

```R
## The penalized non-proportional odds model with a user-supplied lambda gives 
## a fully identified model having bounded estimates. A suitable tuning criterion
## could as well be used to select lambda (e.g., aic or cv) 
f2 <- serp(rating ~ temp + contact, slope = "penalize",
           link = "logit", reverse = TRUE, tuneMethod = "user",
           lambda = 1e1 ,data = wine)
coef(f2)
```

```R
## A penalized partial proportional odds model with one variable set to 
## global effect is also possible.
f3 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "user",
           lambda = 2e1, globalEff = ~ temp, data = wine)
coef(f3)
```

```R
## The unpenalized proportional odds model with constrained estimates. 
## Under estreme shrinkage, estimates in f2 equal those in this model.  
f4 <-  serp(rating ~ temp + contact, slope = "parallel",
            reverse = FALSE, link = "logit", data = wine)
summary(f4)
```

### Installation and Use

Before installing `serp`, it is encouraged to have a recent version of [R](https://cran.r-project.org/bin/windows/base/) installed. The released version of `serp` can be installed from [CRAN](https://cran.r-project.org/package=serp) with:

``` r
install.packages("serp")
```

or the development version from [GitHub](https://github.com/ejikeugba/serp) with:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ejikeugba/serp")
```

Load `serp` into R environment with:
```{r, eval = FALSE}
library(serp)
```


### Community Guidelines

Pull requests are welcomed! Please submit your contributions to `serp` through the list of `Pull Requests`, following the [contributing guidelines](https://ejikeugba.github.io/serp/CONTRIBUTING.html). To report issues and/or seek support, please file a new ticket in the [issue](https://github.com/ejikeugba/serp/issues) tracker, and expect a feedback ASAP! 


### Code of Conduct

Please note that `serp` is released with a [Contributor Code of Conduct](https://github.com/ejikeugba/serp/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.


### References
McCullagh, P. (1980). Regression Models for Ordinal Data. *Journal of the Royal Statistical Society. Series B (Methodological)*, 42, 109-142. https://doi.org/10.1111/j.2517-6161.1980.tb01109.x 

Randall, J (1989). The analysis of sensory data by generalized linear model. *Biometrical Journal*, 31, 781--793. https://doi.org/10.1002/bimj.4710310703

Tutz, G. and Gertheiss, J. (2016). Regularized Regression for Categorical Data (With Discussion and Rejoinder).  *Statistical Modelling*, 16, 161-260. https://doi.org/10.1177/1471082X16642560

Ugba, E. R., Mörlein, D. and Gertheiss, J. (2021). Smoothing in Ordinal Regression: An Application to Sensory Data. *Stats*, 4, 616–633. https://doi.org/10.3390/stats4030037

Ugba, E. R. (2021). serp: An R package for smoothing in ordinal regression *Journal of Open Source Software*, 6(66), 3705. https://doi.org/10.21105/joss.03705
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.R
\name{serp}
\alias{serp}
\title{Smooth Effects on Response Penalty for CLM}
\usage{
serp(
     formula,
     link = c("logit", "probit","loglog", "cloglog", "cauchit"),
     slope = c("penalize", "parallel", "unparallel", "partial"),
     tuneMethod = c("aic", "cv", "finite", "user"),
     reverse = FALSE,
     lambdaGrid = NULL,
     cvMetric = c("brier", "logloss", "misclass"),
     gridType = c("discrete", "fine"),
     globalEff = NULL,
     data,
     subset,
     weights = NULL,
     weight.type = c("analytic", "frequency"),
     na.action = NULL,
     lambda = NULL,
     contrasts = NULL,
     control = list(),
     ...)
}
\arguments{
\item{formula}{regression formula of the form: response ~ predictors. The
response should be a factor (ordered).}

\item{link}{sets the link function for the cumulative link model including:
logit, probit, complementary log-log, cloglog, cauchit.}

\item{slope}{selects the form of coefficients used in the model, with
\code{penalize} denoting the penalized coefficients, \code{unparallel},
\code{parallel} and \code{partial} denoting the unpenalized non-parallel,
parallel and semi-parallel coefficients respectively.}

\item{tuneMethod}{sets the method of choosing an optimal shrinkage
parameter, including: \code{aic}, \code{cv}, \code{finite} and
\code{user}. i.e., the lambda value along parameter shrinkage path at
which the fit's AIC or the k-fold cross-validated test error is
minimal. The finite tuning is used to obtain the model along parameter
shrinkage for which the log-Likelihood exist (is finite). The 'user'
tuning supports a user-supplied lambda value.}

\item{reverse}{false by default, when true the sign of the linear predictor
is reversed.}

\item{lambdaGrid}{optional user-supplied lambda grid for the \code{aic},
and \code{cv} tuning methods, when the discrete \code{gridType}
is chosen. Negative range of values are not allowed. A short lambda grid
could increase computation time assuming large number of predictors and
cases in the model.}

\item{cvMetric}{sets the performance metric for the cv tuning, with the
brier score used by default.}

\item{gridType}{chooses if a discrete or a continuous lambda grid should be
used to select the optimal tuning parameter. The former is used by default
and could be adjusted as desired in \code{serp.control}. The latter
is on the range (0, \code{maxPen}). A user-supplied grid is also possible,
which automatically overrides the internal grid.}

\item{globalEff}{specifies variable(s) to be assigned global effects during
penalization or when \code{slope} is set to \code{partial}. Variables are
specified as a formula with an empty left hand side, for instance,
globalEff = ~predictors.}

\item{data}{optional dataframe explaining the variables used in the formula.}

\item{subset}{specifies which subset of the rows of the data should be used
for fit. All observations are used by default.}

\item{weights}{optional case weights in fitting. Negative weights are not
allowed. Defaults to 1.}

\item{weight.type}{chooses between analytic and frequency weights with the
former used by default. The latter should be used when weights are mere
case counts used to compress the data set.}

\item{na.action}{a function to filter missing data.}

\item{lambda}{a user-supplied single numeric value for the tuning parameter
when using the \code{user} tuning method. Negative values are not
allowed.}

\item{contrasts}{a list of contrasts to be used for some or all of the
factors appearing as variables in the model formula.}

\item{control}{A list of fit control parameters to replace default values
returned by \code{serp.control}. Values not set assume default values.}

\item{...}{additional arguments.}
}
\value{
\item{aic}{the akaike information criterion, with effective degrees of
        freedom obtained from the trace of the generalized hat matrix
        depending on the tuning parameter.}

\item{bic}{the bayesian information criterion, with effective degrees of
        freedom obtained from the trace of the generalized hat matrix
        depending on the tuning parameter.}

\item{call}{the matched call.}

\item{coef}{a vector of coefficients of the fitted model.}

\item{converged}{a character vector of fit convergence status.}

\item{contrasts}{(where relevant) the contrasts used in the model.}

\item{control}{list of control parameters from \code{serp.control}.}

\item{cvMetric}{the performance metric used for cv tuning.}

\item{deviance}{the residual deviance.}

\item{edf}{the (effective) number of degrees of freedom used by the model}

\item{fitted.values}{the fitted probabilities.}

\item{globalEff}{variable(s) in model treated as global effect(s)}

\item{gradient}{a column vector of gradients for the coefficients at the
        model convergence.}

\item{Hessian}{the hessian matrix for the coefficients at the model
        convergence.}

\item{iter}{number of interactions before convergence or non-convergence.}

\item{lambda}{a user-supplied single numeric value for the \code{user}
        tuning tuning method.}

\item{lambdaGrid}{a numeric vector of lambda values used to determine the
        optimum tuning parameter.}

\item{logLik}{the realized log-likelihood at the model convergence.}

\item{link}{character vector indicating the link function of the fit.}

\item{message}{character vector stating the type of convergence obtained}

\item{misc}{a list to hold miscellaneous fit information.}

\item{model}{model.frame having variables from formula.}

\item{na.action}{(where relevant) information on the treatment of NAs.}

\item{nobs}{the number of observations.}

\item{nrFold}{the number of k-fold cross validation for the cv tuning
        method. Default to k = 5.}

\item{rdf}{the residual degrees of freedom}

\item{reverse}{a logical vector indicating the the direction of the
        cumulative probabilities. Default to P(Y<=r).}

\item{slope}{a character vector indicating the type of slope parameters
        fitted. Default to \code{penalize}.}

\item{Terms}{the terms structure describing the model.}

\item{testError}{numeric value of the cross-validated test error at which
        the optimal tuning parameter emerged.}

\item{tuneMethod}{a character vector specifying the method for choosing an
        optimal shrinkage parameter.}

\item{value}{numeric value of AIC or logLik obtained at the optimal tuning
        parameter when using \code{aic} or \code{finite} tuning methods respectively.}

\item{ylev}{the number of the response levels.}
}
\description{
Fits cumulative link models (CLMs) with the
smooth-effect-on-response penalty (SERP) via a modified Newton-Raphson
algorithm. SERP enables the regularization of the parameter space between
the general and the restricted cumulative models, with a resultant shrinkage
of all subject-specific effects to global effects. The Akaike information
critrion (\code{aic}), K-fold cross validation (\code{cv}), among other tuning
aproaches, provide the means of arriving at an optimal tuning parameter in a
in a situation where a user-supplied tuning value is not available.
The \code{slope} argument allows for the selection of a penalized, unparallel,
parallel, or partial slope.
}
\details{
The \code{serp} function fits the cumulative link model (CLM)
with smooth-effect-on-response penalty (SERP). The cumulative
model developed by McCullagh (1980) is probably most frequently
used ordinal model. When motivated by an underlying latent
variable, a simple form of the model is expressed as follows:

\deqn{P(Y\leq r|x) = F(\delta_{0r} + x^T\delta)}

where \eqn{x} is a vector of covariates, \eqn{\delta} a vector
of regression parameters and \eqn{F} a continuous distribution
function. This model assumes that the effect of \eqn{x} does not
depend on the  category. However, with this assumption relaxed,
one obtains the following general cumulative model:

\deqn{P(Y\leq r|x) = F(\delta_{0r} + x^T\delta_{r}),}

where r=1,\dots,k-1. This model, however, has the stochastic ordering
property, which implies that \eqn{P(Y\leq r-1|x) < P(Y\leq r|x)}
holds for all \eqn{x} and all categories \eqn{r}. Such assumption
is often problematic, resulting in unstable likelihoods with
ill-conditioned parameter space during the iterative procedure.

SERP offers a means of arriving at stable estimates of the general model.
It provides a form of regularization that is based on minimizing the
penalized log-likelihood:

\deqn{l_{p}(\delta)=l(\delta)-J_{\lambda}(\delta)}

where \eqn{l(\delta)}, is the log-likelihood of the general cumulative
model and \eqn{J_{\lambda}(\delta)=\lambda J(\delta)} the penalty
function weighted by the turning parameter \eqn{\lambda}. Assuming an
ordered categorical outcome \eqn{Y \in \{1,\dots,k\}}, and considering
that the corresponding parameters \eqn{\delta_{1j},\dots \delta_{k-1,j}}
vary smoothly over the categories, the following penalty
(Tutz and Gertheiss, 2016),

\deqn{J_{\lambda}(\delta)= \sum_{j=1}^{p} \sum_{r=1}^{k-2}
(\delta_{r+1,j}-\delta_{rj})^{2}}

enables the smoothing of response categories such that all
category-specific effects associated with the response turn towards a
common global effect. SERP could also be applied to a semi-parallel model
with only the category-specific part of the model penalized. See,
Ugba (2021), Ugba et al. (2021) for further details and application in
empirical studies.

An object of class \code{serp} with the components listed below,
depending on the type of slope modeled. Other summary methods include:
 \code{summary}, \code{coef}, \code{predict}, \code{vcov},
\code{anova}, etc.
}
\examples{
require(serp)

## The unpenalized non-proportional odds model returns unbounded estimates, hence,
## not fully identifiable.
f1 <- serp(rating ~ temp + contact, slope = "unparallel",
           reverse = TRUE, link = "logit", data = wine)
coef(f1)

## The penalized non-proportional odds model with a user-supplied lambda gives
## a fully identified model with bounded estimates. A suitable tuning criterion
## could as well be used to select lambda (e.g., aic, cv)
f2 <- serp(rating ~ temp + contact, slope = "penalize",
           link = "logit", reverse = TRUE, tuneMethod = "user",
           lambda = 1e1, data = wine)
coef(f2)

## A penalized partial proportional odds model with some variables set to
## global effect is also possible.
f3 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "user",
           lambda = 2e1, globalEff = ~ temp, data = wine)
coef(f3)


## The unpenalized proportional odds model having constrained estimates can
## as well be fit. Under extreme shrinkage, estimates in f2 equal those in
## this model.
f4 <-  serp(rating ~ temp + contact, slope = "parallel",
            reverse = FALSE, link = "logit", data = wine)
summary(f4)

}
\references{
Ugba, E. R. (2021). serp: An R package for smoothing in ordinal regression
    \emph{Journal of Open Source Software}, 6(66), 3705.
    https://doi.org/10.21105/joss.03705

Ugba, E. R., Mörlein, D. and Gertheiss, J. (2021). Smoothing in Ordinal
    Regression: An Application to Sensory Data. \emph{Stats}, 4, 616–633.
    https://doi.org/10.3390/stats4030037

Tutz, G. and Gertheiss, J. (2016). Regularized Regression
    for Categorical Data (With Discussion and Rejoinder).
    \emph{Statistical Modelling}, 16, pp. 161-260.
    https://doi.org/10.1177/1471082X16642560

McCullagh, P. (1980). Regression Models for Ordinal Data.
    \emph{Journal of the Royal Statistical Society. Series B
    (Methodological)}, 42, pp. 109-142.
    https://doi.org/10.1111/j.2517-6161.1980.tb01109.x
}
\seealso{
\code{\link{anova.serp}}, \code{\link{summary.serp}},
\code{\link{predict.serp}}, \code{\link{confint.serp}},
\code{\link{vcov.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wine.R
\docType{data}
\name{wine}
\alias{wine}
\title{Bitterness of wine dataset}
\format{
A data frame with 72 rows and 6 variables:
}
\source{
Taken from Randall (1989).
}
\usage{
wine
}
\value{
\item{\code{response}}{
   scorings of wine bitterness on a 0---100 continuous scale.
 }

\item{\code{rating}}{
   ordered factor with 5 levels; a grouped version of \code{response}.
 }

\item{\code{contact}}{
   factor with two levels (\code{"no"} and \code{"yes"}).
 }

\item{\code{temp}}{
   temperature: factor with two levels.
 }

\item{\code{judge}}{
   factor with nine levels.
 }

\item{\code{bottle}}{
   factor with eight levels.
 }
}
\description{
The \code{wine} dataset adopted from Randall(1989),
represents the outcome of a factorial experiment on factors
determining the bitterness of wine. Two treatment factors
(temperature and contact) with two levels each are provided,
with the rating of wine taken on a continuous scale in the interval
from 0 (none) to 100 (intense). These were subsequently grouped
into five ordered categories ranging from 1 = 'least bitter'
to 5 = 'most bitter'. Altogether, nine different judges assessed
wine from two bottles and out of the four treatment conditions,
making a total of 72 observations.
}
\examples{

\dontrun{
str(wine)
head(wine)
}

}
\references{
Randall, J (1989). The analysis of sensory data by generalized linear
    model. \emph{Biometrical Journal}, 31, 781--793.
    https://doi.org/10.1002/bimj.4710310703
}
\keyword{dataset}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.anova.R
\name{anova.serp}
\alias{anova.serp}
\title{ANOVA method for a fitted serp object}
\usage{
\method{anova}{serp}(object, ..., test = c("Chisq", "none"))
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}

\item{test}{type of test to be conducted.}
}
\value{
\item{model}{the respective model aliases.}

\item{slope}{type of slope fitted, which may be any of, unparallel, parallel,
        or partial slope.}

\item{no.par}{the no of parameters in the model.}

\item{AIC}{the akaike information criterion.}

\item{logLik}{the realized log-likelihood.}

\item{Test}{the different pair(s) of test(s) conducted.}

\item{LR.stat}{the computed Likelihood ratio statistic.}

\item{df}{the degree of freedom.}

\item{Pr(chi)}{the p-value of test statitic.}
}
\description{
Provides a likelihood ratio test for comparing two or more \code{serp} objects. This
does not currently support model(s) with penalized slope.
}
\details{
An ANOVA table with the following components on display:
}
\examples{
library(serp)
m1 <- serp(rating ~ temp + contact, slope = "parallel", link = "logit",
           data = wine)
m2 <- update(m1, ~ contact)
anova(m1, m2)

}
\seealso{
\code{\link{serp}}, \code{\link{confint.serp}}, \code{\link{vcov.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{logLik.serp}
\alias{logLik.serp}
\title{Log-likelihood for a fitted serp object}
\usage{
\method{logLik}{serp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}
}
\value{
A single numeric value of model log-likelihood
}
\description{
Returns the Log-likelihood for a fitted object of class \code{serp}.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "parallel", link = "loglog",
          data = wine)
logLik(m)
}
\seealso{
\code{\link{serp}}, \code{\link{AIC.serp}}, \code{\link{BIC.serp}},
\code{\link{coef.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.vcov.R
\name{vcov.serp}
\alias{vcov.serp}
\title{Variance covariance matrix for a fitted serp object}
\usage{
\method{vcov}{serp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}
}
\value{
A variance covariance matrix of a fitted model.
}
\description{
Provides the Variance covariance matrix of an object of class \code{serp}.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "parallel", link = "logit",
           data = serp::wine)
vcov(m)

}
\seealso{
\code{\link{serp}}

\code{\link{serp}}, \code{\link{anova.serp}}, \code{\link{confint.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{BIC.serp}
\alias{BIC.serp}
\title{BIC for a fitted serp object}
\usage{
\method{BIC}{serp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}
}
\value{
A single numeric value of the model.
}
\description{
Returns the bayesian information criterion of a fitted object of class
\code{serp}. For the penalized slopes, the effective degrees of freedom (edf)
is obtained from the trace of the generalized hat matrix which depends on
the tuning parameter.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "parallel", link = "loglog",
          data = wine)
BIC(m)
}
\seealso{
\code{\link{serp}}, \code{\link{AIC.serp}}, \code{\link{coef.serp}},
\code{\link{logLik.serp}},
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.control.R
\name{serp.control}
\alias{serp.control}
\title{Control parameters for a fitted serp object}
\usage{
serp.control(
             maxits = 5e01,
             eps = 1e-07,
             maxpen = 1e07,
             trace = 0L,
             maxAdjIter = 5e0,
             max.half.iter = 1e01,
             relTol = 1e-03,
             nrFold = 5e0,
             cv.seed = 1e01,
             grid.length = 5e01,
             misclass.thresh = 5e-01,
             minP = .Machine$double.eps,
             ...)
}
\arguments{
\item{maxits}{the maximum number of Newton's iterations. Default to 100.}

\item{eps}{threshold value during optimization at which the iteration
routine terminates. In other words, when the reported change in the
log-likelihood goes below this threshold, convergence is achieved.}

\item{maxpen}{the upper end point of the interval from zero to be searched
for a tuning parameter.}

\item{trace}{prints the Newton's fitting process at each iteration step.If
0 (default) no information is printed, if 1, 2 or 3 different shades of
information are printed.}

\item{maxAdjIter}{the maximum allowable number of Newton step adjustment to
forestall an early optimization failure. Defaults to 5.}

\item{max.half.iter}{the maximum number of iteration step-halfings. Defaults
to 10.}

\item{relTol}{relative convergence tolerance, defaults to 1e-03. checks
relative changes in the parameter estimates between Newton iterations.}

\item{nrFold}{the number of k-fold cross validation for the CV tuning
method. Default to k = 5.}

\item{cv.seed}{single numeric value to change the random seed in CV
tuning.}

\item{grid.length}{the length of the discrete lambda grid for the penalty
method.}

\item{misclass.thresh}{to reset the classification threshold in
\code{errorMetrics} when \code{type} is 'misclass'.}

\item{minP}{A near zero minimum value the fitted probabilities are allowed
to get during iteration to prevent numerical instability .}

\item{...}{additional arguments.}
}
\value{
a list of control parameters.
}
\description{
Default control parameters for 'serp' fit. User-supplied control parameters
could be specified in the main function.
}
\examples{
library(serp)
serp(rating ~ contact, slope = "parallel", link = "logit",
     control = list(maxits = 2e01, eps=1e-05, trace = 2),
     data = wine)

}
\seealso{
\code{\link{serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{print.serp}
\alias{print.serp}
\title{Print method for a fitted serp object}
\usage{
\method{print}{serp}(x, ...)
}
\arguments{
\item{x}{An object of class \code{serp}.}

\item{...}{additional arguments.}
}
\value{
No return value
}
\description{
Prints out a vector of coefficients of the fitted model with some
additional goodness-of-fit measures.
}
\seealso{
\code{\link{serp}}, \code{\link{print.summary.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{AIC.serp}
\alias{AIC.serp}
\title{AIC for a fitted serp object}
\usage{
\method{AIC}{serp}(object, ..., k = 2)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}

\item{k}{fixed value equal to 2.}
}
\value{
A single numeric value of the model AIC.
}
\description{
Returns the akaike information criterion of a fitted object of class
\code{serp}. For the penalized slope, the effective degrees of freedom (edf)
is obtained from the trace of the generalized hat matrix which depends on
the tuning parameter.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "parallel", link = "probit",
          data = wine)
AIC(m)
}
\seealso{
\code{\link{serp}}, \code{\link{BIC.serp}}, \code{\link{coef.serp}},
\code{\link{logLik.serp}},
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{summary.serp}
\alias{summary.serp}
\title{Summary method for a fitted serp object.}
\usage{
\method{summary}{serp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{Not used. Additional summary arguments.}
}
\value{
\item{coefficients}{the matrix of coefficients, standard errors,
        z-values and p-values.}

\item{null.deviance}{the deviance for the intercept only model.}

\item{null.logLik}{the log-likelihood for the intercept only model.}

\item{penalty}{list of penalization information obtained with
        \code{slope} set to "penalize".}

\item{expcoefs}{the exponentiated coefficients.}
}
\description{
This function summarizes the result of a fitted serp object in a dataframe.
}
\details{
an object of class \code{summary.serp}. A list (depending on the type of
\code{slope} used) of all model components defined in the \code{\link{serp}},
function with additional components listed below.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "user",
           lambda = 0, data = wine)
summary(m)
}
\seealso{
\code{\link{anova.serp}}, \code{\link{predict.serp}},
\code{\link{confint.serp}}, \code{\link{vcov.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.confint.R
\name{confint.serp}
\alias{confint.serp}
\title{Confidence interval for a fitted serp object}
\usage{
\method{confint}{serp}(object, ..., parm, level = 0.95)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}

\item{parm}{unused argument.}

\item{level}{significance level.}
}
\value{
A matrix of the the confidence intervals of fitted model.
}
\description{
Provides the confidence interval of estimates for an object of class \code{serp}.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "parallel", link = "logit",
           data = wine)
confint(m)

}
\seealso{
\code{\link{serp}}, \code{\link{anova.serp}}, \code{\link{vcov.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{print.summary.serp}
\alias{print.summary.serp}
\title{Print method for an object of class summary.serp}
\usage{
\method{print}{summary.serp}(x, ...)
}
\arguments{
\item{x}{An object of class \code{summary.serp}.}

\item{...}{additional arguments.}
}
\value{
No return value
}
\description{
Prints out the information supplied via \code{summary.serp} method.
}
\seealso{
\code{\link{serp}}, \code{\link{print.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{predict.serp}
\alias{predict.serp}
\title{Prediction from fitted serp model}
\usage{
\method{predict}{serp}(object, type = c("link", "response", "class"), newdata = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{type}{could be any of these: response, link or terms.}

\item{newdata}{fresh dataset with all relevant variables.}

\item{...}{additional arguments.}
}
\value{
A vector of predicted classes with \code{type} equal to 'class'
or a dataframe of predicted values for \code{type} equal to 'response'
and 'link'.
}
\description{
This function takes a fitted \code{serp} object produced by serp() and
produces predicted values. Type of predictions returned include response,
link and class. Prediction is also possible with new set of values having
the same column names as in the original values used for the model fit.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "user",
           lambda = 1, data = wine)

head(predict(m, type = "link"))
head(predict(m, type = "response"))
predict(m, type = "class")

n.wine <- wine[1:20,]
predict(m, newdata = n.wine, type = "class")

}
\seealso{
\code{\link{anova.serp}}, \code{\link{summary.serp}},
\code{\link{confint.serp}}, \code{\link{vcov.serp}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{coef.serp}
\alias{coef.serp}
\alias{coefficients.serp}
\title{Coefficients for a fitted serp object}
\usage{
\method{coef}{serp}(object, ...)

\method{coefficients}{serp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}
}
\value{
A vector of model coefficients.
}
\description{
Returns the coefficients of a fitted object of class \code{serp}.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "parallel", link = "loglog",
          data = wine)
coef(m)

}
\seealso{
\code{\link{serp}}, \code{\link{AIC.serp}}, \code{\link{BIC.serp}},
\code{\link{logLik.serp}}
}
