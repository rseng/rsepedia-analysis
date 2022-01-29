
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
