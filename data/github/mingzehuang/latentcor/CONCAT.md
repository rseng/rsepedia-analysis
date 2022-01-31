<!-- badges: start -->
[![R-CMD-check](https://github.com/mingzehuang/latentcor/workflows/R-CMD-check/badge.svg)](https://github.com/mingzehuang/latentcor/actions)
[![codecov](https://codecov.io/gh/mingzehuang/latentcor/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mingzehuang/latentcor)
[![CRAN status](https://www.r-pkg.org/badges/version-last-release/latentcor)](https://CRAN.R-project.org/package=latentcor)
[![Launch binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mingzehuang/latentcor/master)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03634/status.svg)](https://doi.org/10.21105/joss.03634)
[![DOI](https://zenodo.org/badge/336304814.svg)](https://zenodo.org/badge/latestdoi/336304814)
<!-- badges: end -->


# latentcor: Latent Correlation for Mixed Types of Data

`latentcor` is an `R` package for estimation of latent correlations with mixed data types (continuous, binary, truncated, and ternary) under the latent Gaussian copula model. For references on the estimation framework, see

  * [Fan, J., Liu, H., Ning, Y., and Zou, H. (2017), “High Dimensional Semiparametric Latent Graphical Model for Mixed Data.” *JRSS B*](https://doi.org/10.1111/rssb.12168). **Continuous/binary** types.

  * [Quan X., Booth J.G. and Wells M.T."Rank-based approach for estimating correlations in mixed ordinal data." *arXiv*](https://arxiv.org/abs/1809.06255) **Ternary** type.

  * [Yoon G., Carroll R.J. and Gaynanova I. (2020). “Sparse semiparametric canonical correlation analysis for data of mixed types”. *Biometrika*](https://doi.org/10.1093/biomet/asaa007). **Truncated** type for zero-inflated data.

  * [Yoon G., Müller C.L. and Gaynanova I. (2021). “Fast computation of latent correlations” *JCGS*](https://doi.org/10.1080/10618600.2021.1882468). **Approximation method of computation**, see [vignette](https://mingzehuang.github.io/latentcor/articles/latentcor.html) for details.

## Statement of Need

No R software package is currently available that allows accurate and fast correlation estimation from mixed variable data in a unifying manner. The R package [`latentcor`](https://CRAN.R-project.org/package=latentcor), introduced here, thus represents the first stand-alone R package for 
computation of latent correlation that takes into account all variable types (continuous/binary/ordinal/zero-inflated), comes with an optimized memory footprint, 
and is computationally efficient, essentially making latent correlation estimation almost as fast as rank-based correlation estimation. 

## Installation

To use `latentcor`, you need to install [`R`](https://cran.r-project.org/). To enhance your user experience, you may use some IDE for it (e.g. [`RStudio`](https://www.rstudio.com/)).

The development version of [`latentcor`](https://github.com/mingzehuang/latentcor) is available on GitHub. You can download it with the help of the `devtools` package in `R` as follow:

```r
install.packages("devtools")
devtools::install_github("https://github.com/mingzehuang/latentcor", build_vignettes = TRUE)
```
The stable release version [`latentcor`](https://CRAN.R-project.org/package=latentcor) is available on CRAN. You can download it in `R` as follow:

```r
install.packages("latentcor")
```

## Example

A simple example estimating latent correlation is shown below.

```r
library(latentcor)

# Generate two variables of sample size 100
# The first variable is ternary (pi0 = 0.3, pi1 = 0.5, pi2 = 1-0.3-0.5 = 0.2) 
# The second variable is continuous. 
# No copula transformation is applied.
X = gen_data(n = 1000, types = c("ter", "con"), XP = list(c(0.3, .5), NA))$X

# Estimate latent correlation matrix with the original method
latentcor(X = X, types = c("ter", "con"), method = "original")$R

# Estimate latent correlation matrix with the approximation method
latentcor(X = X, types = c("ter", "con"))$R

# Speed improvement by approximation method compared with original method
library(microbenchmark)
microbenchmark(latentcor(X, types = c("ter", "con"), method = "original"),
               latentcor(X, types = c("ter", "con")))
# Unit: milliseconds
# min     lq     mean    median     uq     max     neval
# 5.3444 5.8301 7.033555 6.06740 6.74975 20.8878   100
# 1.5049 1.6245 2.009371 1.73805 1.99820  5.0027   100
# This is run on Windows 10 with Intel(R) Core(TM) i5-4570 CPU @ 3.20GHz   3.20 GHz

# Heatmap for latent correlation matrix.
latentcor(X = X, types = c("ter", "con"), showplot = TRUE)$plotR
```
Another example with the `mtcars` dataset.

```r
library(latentcor)
# Use build-in dataset mtcars
X = mtcars
# Check variable types for manual determination
apply(mtcars, 2, table)
# Or use built-in get_types function to get types suggestions
get_types(mtcars)

# Estimate latent correlation matrix with original method
latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"), method = "original")$R
# Estimate latent correlation matrix with approximation method
latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"))$R

# Speed improvement by approximation method compared with original method
library(microbenchmark)
microbenchmark(latentcor(mtcars, types = types, method = "original"),
               latentcor(mtcars, types = types, method = "approx"))
# Unit: milliseconds
#  min       lq        mean      median        uq      max    neval
#  201.9872 215.6438   225.30385 221.5364 226.58330 411.4940   100
#   71.8457  75.1681   82.42531  80.1688  84.77845 238.3793    100
# This is run on Windows 10 with Intel(R) Core(TM) i5-4570 CPU @ 3.20GHz   3.20 GHz

# Heatmap for latent correlation matrix with approximation method.
latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin",
                       "bin", "ter", "con"), showplot = TRUE)$plotR
```

Interactive heatmap see: [interactive heatmap of latent correlations (approx) for mtcars](https://rpubs.com/mingzehuang/797937)

Community Guidelines
--------------------

1.  Contributions and suggestions to the software are always welcome.
    Please consult our [contribution guidelines](https://github.com/mingzehuang/latentcor/blob/master/CONTRIBUTING.md) prior
    to submitting a pull request.
2.  Report issues or problems with the software using github’s [issue
    tracker](https://github.com/mingzehuang/latentcor/issues).
3.  Contributors must adhere to the [Code of Conduct](https://github.com/mingzehuang/latentcor/blob/master/CODE_OF_CONDUCT.md).

Acknowledgments
--------------

We thank Dr. Grace Yoon for providing implementation details of the [`mixedCCA`](https://github.com/irinagain/mixedCCA) R package.
# latentcor 1.1.0

* Rename main function from `estR` to `latentcor`.
* Add help function `get_types` to automatically generate types for data matrix.
* Add more documentation on parameters.
* Add more documentation for speed comparison between original method and approximation.

# latentcor 1.2.0

* Add user-defined `use.nearPD`, so that user can decide if latent correlation matrix should be adjusted to be positive definite automatically.
* Remove redundant code for checking positive definiteness.
* Minor correction on types detection to accommodate `NA` values. 
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
reported to the community leaders responsible for enforcement at mingzehuang@gmail.com.
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

# Contributing to `latentcor`

Contributions to `latentcor` follow the same principles as stated in the guideline for contributors to the `tidyverse` ecosystem of R (see here for more details: [**development contributing guide**](https://rstd.io/tidy-contrib) ).

## Fixing typos

You can refer to [this](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork) for proposing changes with pull request from a fork.
You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file.

If you find any typos in an `.Rd` file, then please make changes in the corresponding `.R` file, as `.Rd` files are automatically generated by [roxygen2](https://roxygen2.r-lib.org/articles/roxygen2.html) and should not be edited by hand.

## Bigger changes

The first step here for you to make any changes is to install `devtools` using `install.packages("devtools")`.
If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it is needed.
If you have found a bug, please file an issue that illustrates the bug with a minimal reproducible example, a
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

To contribute a change to `latentcor`, you follow these steps:

1. Create a branch in git, give it a descriptive name, and make your changes.
2. Push branch to GitHub and issue a pull request (PR).
3. Discuss the pull request.
4. Iterate until either we accept the PR or decide that it's not a good fit for `latentcor`.

If you're not familiar with git or GitHub, please start by reading
<http://r-pkgs.had.co.nz/git.html>

Pull requests will be evaluated against a checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivate the need for change. Please describe the problem your PR addresses and show how your pull request solves it as concisely as possible.

Also, include this motivation in `NEWS` so that when a new release of
`latentcor` comes out it's easy for users to see what has changed. Add your
item at the top of the file and use markdown for formatting. The
news item should end with `(@yourGithubUsername, #the_issue_number)`.

2.  __Only related changes__. Before you submit your pull request, please check to make sure that you haven't accidentally included any unrelated changes. These make it harder to see exactly what's changed, and to evaluate any unexpected side effects.

Each PR corresponds to a git branch, so if you expect to submit multiple changes
make sure to create multiple branches. If you have multiple changes that depend
on each other, start with the first one and don't submit any others until the
first one has been processed.

3. __Use `latentcor` coding style__. We tried to adhere as closely as possible to the [official `tidyverse` style guide](http://style.tidyverse.org) -- please do so as well. Maintaining a consistent style across the whole code base makes it much easier to jump into the code. If you're modifying existing `latentcor` code that doesn't follow the style guide, a separate pull request to fix the style would be greatly appreciated.

4. If you're adding new parameters or a new function, you'll also need to document them with [`roxygen2`](https://github.com/klutometis/roxygen). Make sure to re-run `devtools::document()` on the code before submitting.
---
title: 'latentcor: An R Package for estimating latent correlations from mixed data types'
tags:
- R
- Statistics
- Latent Correlation
date: "15 July 2021"
output:
  html_document:
    df_print: paged
bibliography: paper.bib
authors:
- name: Mingze Huang
  orcid: 0000-0003-3919-1564
  affiliation: 1, 2
- name: Christian L. Müller
  orcid: 0000-0002-3821-7083
  affiliation: 3, 4, 5
- name: Irina Gaynanova
  orcid: 0000-0002-4116-0268
  affiliation: 1
affiliations:
- name: Department of Statistics, Texas A&M University, College Station, TX
  index: 1
- name: Department of Economics, Texas A&M University, College Station, TX
  index: 2
- name: Ludwig-Maximilians-Universität München, Germany
  index: 3
- name: Helmholtz Zentrum München, Germany
  index: 4
- name: Flatiron Institute, New York
  index: 5
---


# Summary

We present `latentcor`, an R package for correlation estimation from data with mixed variable types. Mixed variables types, including continuous, binary, ordinal, zero-inflated, or truncated data are routinely collected in many areas of science. Accurate estimation of correlations among such variables is often the first critical step in statistical analysis workflows. Pearson correlation as the default choice is not well suited for mixed data types as the underlying normality assumption is violated. The concept of semi-parametric latent Gaussian copula models, on the other hand, provides a unifying way to estimate correlations between mixed data types. The R package `latentcor` comprises a comprehensive list of these models, enabling the estimation of correlations between any of continuous/binary/ternary/zero-inflated (truncated) variable types. The underlying implementation takes advantage of a fast multi-linear interpolation scheme with an efficient choice of interpolation grid points, thus giving the package a small memory footprint without compromising estimation accuracy. This makes latent correlation estimation readily available for modern high-throughput data analysis.

# Statement of need
No R software package is currently available that allows accurate and fast correlation estimation from mixed variable data in a unifying manner. 
The popular `cor` function within R package `stats` [@team2013r], for instance, allows to compute Pearson's correlation, Kendall's $\tau$ and Spearman's
$\rho$, and a faster algorithm for calculating Kendall's $\tau$ is implemented in the R package `pcaPP` [@croux2013robust]. Pearson's correlation is not
appropriate for skewed or ordinal data, and its use leads to invalid inference in those cases. While the rank-based Kendall's $\tau$ and Spearman's $\rho$ are
more robust measures of *association*, they cannot directly be used as subsitutes for statistical methods that require Pearson correlation as input (a prominent example is, e.g., graphical model estimation [@xue2012regularized;@yoon2019microbial]). The R package `polycor` [@fox2019poly] is designed for ordinal data and allows to computes
polychoric (ordinal/ordinal) and polyserial (ordinal/continuous) correlations based on the latent Gaussian model. However, the package does not have functionality
for zero-inflated data, nor can it handle skewed continuous measurements as it does not allow for copula transformation. The R package `correlation`
[@makowski2020methods] in the `easystats` collection provides 16 different correlation measures, including polychoric and polyserial correlations. However, 
functionality for correlation estimation from zero-inflated data is lacking. The R package `mixedCCA` [@yoon2020sparse] is based on the latent Gaussian copula 
model and can compute latent correlations between continuous/binary/zero-inflated variable types as an intermediate step for canonical correlation analysis. 
However, `mixedCCA` does not allow for ordinal data types. The R package `latentcor`, introduced here, thus represents the first stand-alone R package for 
computation of latent correlation that takes into account all variable types (continuous/binary/ordinal/zero-inflated), comes with an optimized memory footprint, 
and is computationally efficient, essentially making latent correlation estimation almost as fast as rank-based correlation estimation. 

# Estimation of latent correlations

## The general estimation workflow

The estimation of latent correlations consists of three steps: 

* computing Kendall's $\tau$ between each pair of variables,

* choosing the bridge function $F(\cdot)$ based on the types of variable pairs; the bridge function connects the Kendall's $\tau$ computed from the data, $\widehat \tau$, to the true underlying correlation $\rho$ via moment equation $\mathbb{E}(\widehat \tau) = F(\rho)$;

* estimating latent correlation by calculating $F^{-1}(\widehat \tau)$. 

We summarize the references for the explicit form of $F(\cdot)$ for each variable combination as implemented in `latentcor` below.

+----------------+-----------------------+-----------------+--------------------------+-----------------+
| Type           | continuous            | binary          | ternary                  | zero-inflated\  |
|                |                       |                 |                          | (truncated)     |
+================+=======================+=================+==========================+=================+
| continuous     | @liu2009nonparanormal | \-              | \-                       | \-              |
+----------------+-----------------------+-----------------+--------------------------+-----------------+
| binary         | @fan2017high          | @fan2017high    | \-                       | \-              |
+----------------+-----------------------+-----------------+--------------------------+-----------------+
| ternary        | @quan2018rank         | @quan2018rank   | @quan2018rank            | \-              |
+----------------+-----------------------+-----------------+--------------------------+-----------------+
| zero-inflated\ | @yoon2020sparse       | @yoon2020sparse | See `latentcor`\         | @yoon2020sparse |
| (truncated)    |                       |                 | vignette for derivation\ |                 |
+----------------+-----------------------+-----------------+--------------------------+-----------------+


## Efficient inversion of the bridge function

In `latentcor`, the inversion of the bridge function $F(\cdot)$ can be computed in two ways. The original approach (`method = "original"`) relies on numerical
inversion for each pair of variables based on uni-root optimization [@yoon2020sparse]. Since each pair of variables requires a separate optimization run, the
original approach is computationally expensive when the number of variables is large. The second approach to invert $F(\cdot)$ is through fast multi-linear
interpolation of pre-calculated $F^{-1}$ values at specific sets of interpolation grid points (`method = "approx"`). This construction has been proposed in
[@yoon2021fast] and is available for continuous/binary/truncated pairs in the current version of `mixedCCA`. However, that implementation lacks the ternary
variable case and relies on an interpolation grid with a large memory footprint. `latentcor` includes the ternary case and provides an optimized interpolation 
grid by redefining the bridge functions on a rescaled version of Kendall's $\tau$. Here, the scaling adapts to the smoothness of the underlying type of variables 
by simultaneously controlling the approximation error at the same or lower level. As a result, `latentcor` has significantly smaller memory footprint (see
Table below) and smaller approximation error compared to `mixedCCA`.

\newpage

Memory footprints (in KB):

 | case | mixedCCA | latentcor |
 |-----|----------|----------|
| binary/continuous | 10.08 | 4.22 |
| binary/binary | 303.04 | 69.1 |
| truncated/continuous | 20.99 | 6.16 |
| truncated/binary | 907.95 | 92.25 | 
| truncated/truncated | 687.68 | 84.33 |
| ternary/continuous | - | 125.83 |
| ternary/binary | - | 728.3 |
| ternary/truncated | - | 860.9 |
| ternary/ternary | - | 950.61 |

## Illustrative examples 

To illustrate the excellent performance of latent correlation estimation on mixed data, we consider the simple example of estimating correlations between continuous and ternary variables. 

First, we use `latentcor` to generate synthetic data with two variables of sample size 500, and true latent correlation value of 0.5. We then estimate the correlation using the original method, the approximation method (default), and standard Pearson correlation.
```r
library(latentcor)

# The first variable is ternary 
# The second variable is continuous. 
# No copula transformation is applied.
set.seed(2346)
X = gen_data(n = 500, types = c("ter", "con"), rhos = 0.5)$X
# Estimate correlations
latentcor(X = X, types = c("ter", "con"), method = "original")$R
latentcor(X = X, types = c("ter", "con"))$R
cor(X)
```
The original method estimates the latent correlation to be equal to 0.4766 (and the approximation method is very close with the value 0.4762). 
By contrast, applying Pearson correlation gives an estimate of 0.4224, which is further from the true value 0.5.

To illustrate the bias induced by Pearson correlation estimation, we consider the truncated/continuous case for different values of the true correlation. Figure \ref{fig:R_all}A displays the values obtained by using standard Pearson correlation, revealing a significant estimation bias with respect to the true correlations. Figure \ref{fig:R_all}B displays the estimated latent correlations using the original approach versus the true values of the underlying ternary/continuous correlations. 
The alignment of points around $y=x$ line confirms that the estimation is empirically unbiased. Figure \ref{fig:R_all}C displays the estimated latent correlations using the approximation approach (`method = "approx"`) versus true values of underlying latent correlation. The results are almost indistinguishable from Figure \ref{fig:R_all}B at a fraction of the computational cost.

![Scatter plots of estimated Pearson correlation (panel A) and latent correlations (`original` in panel B, `approx` in panel C) vs. ground truth correlations \label{fig:R_all}](./CombinedCorrelations.pdf)

The script to reproduce the displayed results is available at [latentcor_evaluation](https://github.com/mingzehuang/latentcor_evaluation/blob/master/unbias_check.R).

We next illustrate an application of `latentcor` to the `mtcars` dataset, available in standard R. The `mtcars` dataset comprises eleven variables of continuous, binary, and ternary data type. The function `get_types` can be used to automatically extract these types from the data. After the types are determined, the correlation matrix can be estimated using either the original method or the approximation method.

```r
library(latentcor)
X = mtcars
# Extract variable types
type = get_types(X)
# Estimate correlations
latentcor(mtcars, types = types, method = "original")$R
latentcor(mtcars, types = types)$R
```

Figure \ref{fig:R_cars} shows the $11 \times 11$ matrices with latent correlation estimates (with default `approx` method, left panel), Pearson correlation estimates (middle panel), and their difference in estimation (right panel). Even on this small dataset, we observe absolute differences exceeding $0.2$.    

![Heatmap of latent correlations (`approx`, left panel), Pearson correlation (middle panel), and difference between the two estimators (latent correlation - Pearson correlation) on the mtcars dataset \label{fig:R_cars}](./all_heatmap.pdf)

The script to reproduce Figure \ref{fig:R_cars} is available [here](https://github.com/mingzehuang/latentcor_evaluation/blob/master/all_heatmap.R).
We also provide interactive heatmaps for [estimated latent correlations](https://rpubs.com/mingzehuang/797937), [Pearson correlations](https://rpubs.com/mingzehuang/797945), and [their differences (estimated latent correlations minus Pearson correlations)](https://rpubs.com/mingzehuang/798060) for the `mtcars` data set.

# Basic Usage and Availability

The R package `latentcor` is available on [Github](https://github.com/mingzehuang/latentcor/). A getting started vignette with basic examples is available [here](https://mingzehuang.github.io/latentcor/articles/latentcor.html). A vignette with mathematical background of estimation process as well as effect of optional parameters is available [here](https://mingzehuang.github.io/latentcor/articles/latentcor_math.html).

# Acknowledgments

We thank Dr. Grace Yoon for providing implementation details of the `mixedCCA` R package.

# References
---
title: "Mathematical Framework for latentcor"
author: "Mingze Huang, Christian L. Müller, Irina Gaynanova"
date: "`r Sys.Date()`"
bibliography: latentcor.bib
output: rmarkdown::html_vignette
extra_dependencies: ["amsmath"]
nocite: |
  @croux2013robust
  @filzmoser2021pcapp
  @liu2009nonparanormal
  @fox2019poly
vignette: >
  %\VignetteIndexEntry{Mathematical Framework for latentcor}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(latentcor)
```

# Latent Gaussian Copula Model for Mixed Data

`latentcor` utilizes the powerful semi-parametric latent Gaussian copula models to estimate latent correlations between mixed data types (continuous/binary/ternary/truncated or zero-inflated). Below we review the definitions for each type.

***Definition of continuous model*** [@fan2017high]

A random $X\in\cal{R}^{p}$ satisfies the Gaussian copula (or nonparanormal) model if there exist monotonically increasing $f=(f_{j})_{j=1}^{p}$ with $Z_{j}=f_{j}(X_{j})$ satisfying $Z\sim N_{p}(0, \Sigma)$, $\sigma_{jj}=1$; we denote $X\sim NPN(0, \Sigma, f)$.

```{r continuous}
X = gen_data(n = 6, types = "con")$X
X
```

***Definition of binary model*** [@fan2017high]

A random $X\in\cal{R}^{p}$ satisfies the binary latent Gaussian copula model if there exists $W\sim NPN(0, \Sigma, f)$ such that $X_{j}=I(W_{j}>c_{j})$, where $I(\cdot)$ is the indicator function and $c_{j}$ are constants.

```{r binary}
X = gen_data(n = 6, types = "bin")$X
X
```

***Definition of ternary model*** [@quan2018rank]

A random $X\in\cal{R}^{p}$ satisfies the ternary latent Gaussian copula model if there exists $W\sim NPN(0, \Sigma, f)$ such that $X_{j}=I(W_{j}>c_{j})+I(W_{j}>c'_{j})$, where $I(\cdot)$ is the indicator function and $c_{j}<c'_{j}$ are constants.

```{r ternary}
X = gen_data(n = 6, types = "ter")$X
X
```

***Definition of truncated or zero-inflated model*** [@yoon2020sparse]

A random $X\in\cal{R}^{p}$ satisfies the truncated latent Gaussian copula model if there exists $W\sim NPN(0, \Sigma, f)$ such that $X_{j}=I(W_{j}>c_{j})W_{j}$, where $I(\cdot)$ is the indicator function and $c_{j}$ are constants.

```{r truncated}
X = gen_data(n = 6, types = "tru")$X
X
```


***Mixed latent Gaussian copula model***

The mixed latent Gaussian copula model jointly models $W=(W_{1}, W_{2}, W_{3}, W_{4})\sim NPN(0, \Sigma, f)$ such that $X_{1j}=W_{1j}$, $X_{2j}=I(W_{2j}>c_{2j})$, $X_{3j}=I(W_{3j}>c_{3j})+I(W_{3j}>c'_{3j})$  and $X_{4j}=I(W_{4j}>c_{4j})W_{4j}$.

```{r mixed}
set.seed("234820")
X = gen_data(n = 100, types = c("con", "bin", "ter", "tru"))$X
head(X)
```

# Moment-based estimation of $\Sigma$ based on bridge functions


The estimation of latent correlation matrix $\Sigma$ is achieved via the **bridge function** $F$ which is defined such that $E(\hat{\tau}_{jk})=F(\sigma_{jk})$, where $\sigma_{jk}$ is the latent correlation between variables $j$ and $k$, and $\hat{\tau}_{jk}$ is the corresponding sample Kendall's $\tau$. 


***Kendall's $\tau$ ($\tau_{a}$)***

Given observed $\mathbf{x}_{j}, \mathbf{x}_{k}\in\cal{R}^{n}$,

$$
\hat{\tau}_{jk}=\hat{\tau}(\mathbf{x}_{j}, \mathbf{x}_{k})=\frac{2}{n(n-1)}\sum_{1\le i<i'\le n}sign(x_{ij}-x_{i'j})sign(x_{ik}-x_{i'k}),
$$
where $n$ is the sample size.

`latentcor` calculates pairwise Kendall's $\widehat \tau$ as part of the estimation process

```{r KendallTau}
estimate = latentcor(X, types = c("con", "bin", "ter", "tru"))
K = estimate$K
K
```

Using $F$ and $\widehat \tau_{jk}$, a moment-based estimator is $\hat{\sigma}_{jk}=F^{-1}(\hat{\tau}_{jk})$ with the corresponding  $\hat{\Sigma}$ being consistent for $\Sigma$ [@fan2017high; @quan2018rank; @yoon2020sparse]. 


The explicit form of **bridge function** $F$ has been derived for all combinations of continuous(C)/binary(B)/ternary(N)/truncated(T) variable types, and we summarize the corresponding references. Each of this combinations is implemented in `latentcor`.

|Type | continuous | binary | ternary | zero-inflated (truncated) |
|-----|----------|----------|----------|----------|
|continuous | @liu2009nonparanormal |- | -| - |
|binary | @fan2017high | @fan2017high | - | - |
|ternary | @quan2018rank | @quan2018rank | @quan2018rank | - |
|zero-inflated (truncated) | @yoon2020sparse | @yoon2020sparse | See Appendix | @yoon2020sparse |

Below we provide an explicit form of $F$ for each combination.

**Theorem (explicit form of bridge function)** Let $W_{1}\in\cal{R}^{p_{1}}$, $W_{2}\in\cal{R}^{p_{2}}$, $W_{3}\in\cal{R}^{p_{3}}$, $W_{4}\in\cal{R}^{p_{4}}$ be such that $W=(W_{1}, W_{2}, W_{3}, W_{4})\sim NPN(0, \Sigma, f)$ with $p=p_{1}+p_{2}+p_{3}+p_{4}$. Let $X=(X_{1}, X_{2}, X_{3}, X_{4})\in\cal{R}^{p}$ satisfy $X_{j}=W_{j}$ for $j=1,...,p_{1}$, $X_{j}=I(W_{j}>c_{j})$ for $j=p_{1}+1, ..., p_{1}+p_{2}$, $X_{j}=I(W_{j}>c_{j})+I(W_{j}>c'_{j})$ for $j=p_{1}+p_{2}+1, ..., p_{3}$ and $X_{j}=I(W_{j}>c_{j})W_{j}$ for $j=p_{1}+p_{2}+p_{3}+1, ..., p$ with $\Delta_{j}=f(c_{j})$. The rank-based estimator of $\Sigma$ based on the observed $n$ realizations of $X$ is the matrix $\mathbf{\hat{R}}$ with $\hat{r}_{jj}=1$, $\hat{r}_{jk}=\hat{r}_{kj}=F^{-1}(\hat{\tau}_{jk})$ with block structure

$$
\mathbf{\hat{R}}=\begin{pmatrix}
F_{CC}^{-1}(\hat{\tau}) & F_{CB}^{-1}(\hat{\tau}) & F_{CN}^{-1}(\hat{\tau}) & F_{CT}^{-1}(\hat{\tau})\\
F_{BC}^{-1}(\hat{\tau}) & F_{BB}^{-1}(\hat{\tau}) & F_{BN}^{-1}(\hat{\tau}) & F_{BT}^{-1}(\hat{\tau})\\
F_{NC}^{-1}(\hat{\tau}) & F_{NB}^{-1}(\hat{\tau}) & F_{NN}^{-1}(\hat{\tau}) & F_{NT}^{-1}(\hat{\tau})\\
F_{TC}^{-1}(\hat{\tau}) & F_{TB}^{-1}(\hat{\tau}) & F_{TN}^{-1}(\hat{\tau}) & F_{TT}^{-1}(\hat{\tau})
\end{pmatrix}
$$
$$
F(\cdot)=\begin{cases}
CC:\ 2\sin^{-1}(r)/\pi \\
\\
BC: \ 4\Phi_{2}(\Delta_{j},0;r/\sqrt{2})-2\Phi(\Delta_{j}) \\
\\
BB: \ 2\{\Phi_{2}(\Delta_{j},\Delta_{k};r)-\Phi(\Delta_{j})\Phi(\Delta_{k})\}  \\
\\
NC: \ 4\Phi_{2}(\Delta_{j}^{2},0;r/\sqrt{2})-2\Phi(\Delta_{j}^{2})+4\Phi_{3}(\Delta_{j}^{1},\Delta_{j}^{2},0;\Sigma_{3a}(r))-2\Phi(\Delta_{j}^{1})\Phi(\Delta_{j}^{2})\\
\\
NB: \ 2\Phi_{2}(\Delta_{j}^{2},\Delta_{k},r)\{1-\Phi(\Delta_{j}^{1})\}-2\Phi(\Delta_{j}^{2})\{\Phi(\Delta_{k})-\Phi_{2}(\Delta_{j}^{1},\Delta_{k},r)\} \\
\\
NN: \ 2\Phi_{2}(\Delta_{j}^{2},\Delta_{k}^{2};r)\Phi_{2}(-\Delta_{j}^{1},-\Delta_{k}^{1};r)-2\{\Phi(\Delta_{j}^{2})-\Phi_{2}(\Delta_{j}^{2},\Delta_{k}^{1};r)\}\{\Phi(\Delta_{k}^{2})-\Phi_{2}(\Delta_{j}^{1},\Delta_{k}^{2};r)\} \\
\\
TC: \ -2\Phi_{2}(-\Delta_{j},0;1/\sqrt{2})+4\Phi_{3}(-\Delta_{j},0,0;\Sigma_{3b}(r)) \\
\\
TB: \ 2\{1-\Phi(\Delta_{j})\}\Phi(\Delta_{k})-2\Phi_{3}(-\Delta_{j},\Delta_{k},0;\Sigma_{3c}(r))-2\Phi_{3}(-\Delta_{j},\Delta_{k},0;\Sigma_{3d}(r))  \\
\\
TN: \ -2\Phi(-\Delta_{k}^{1})\Phi(\Delta_{k}^{2}) + 2\Phi_{3}(-\Delta_{k}^{1},\Delta_{k}^{2},\Delta_{j};\Sigma_{3e}(r))+2\Phi_{4}(-\Delta_{k}^{1},\Delta_{k}^{2},-\Delta_{j},0;\Sigma_{4a}(r))+2\Phi_{4}(-\Delta_{k}^{1},\Delta_{k}^{2},-\Delta_{j},0;\Sigma_{4b}(r)) \\
\\
TT: \ -2\Phi_{4}(-\Delta_{j},-\Delta_{k},0,0;\Sigma_{4c}(r))+2\Phi_{4}(-\Delta_{j},-\Delta_{k},0,0;\Sigma_{4d}(r)) \\
\end{cases}
$$

where $\Delta_{j}=\Phi^{-1}(\pi_{0j})$, $\Delta_{k}=\Phi^{-1}(\pi_{0k})$, $\Delta_{j}^{1}=\Phi^{-1}(\pi_{0j})$, $\Delta_{j}^{2}=\Phi^{-1}(\pi_{0j}+\pi_{1j})$, $\Delta_{k}^{1}=\Phi^{-1}(\pi_{0k})$, $\Delta_{k}^{2}=\Phi^{-1}(\pi_{0k}+\pi_{1k})$,

$$
\Sigma_{3a}(r)=
\begin{pmatrix}
1 & 0 & \frac{r}{\sqrt{2}} \\
0 & 1 & -\frac{r}{\sqrt{2}} \\
\frac{r}{\sqrt{2}} & -\frac{r}{\sqrt{2}} & 1
\end{pmatrix}, \;\;\;
\Sigma_{3b}(r)=
\begin{pmatrix}
1 & \frac{1}{\sqrt{2}} & \frac{r}{\sqrt{2}}\\
\frac{1}{\sqrt{2}} & 1 & r \\
\frac{r}{\sqrt{2}} & r & 1
\end{pmatrix}, \;\;\;
\Sigma_{3c}(r)=
\begin{pmatrix}
1 & -r & \frac{1}{\sqrt{2}} \\
-r & 1 & -\frac{r}{\sqrt{2}} \\
\frac{1}{\sqrt{2}} & -\frac{r}{\sqrt{2}} & 1
\end{pmatrix},
$$

$$
\Sigma_{3d}(r)=
\begin{pmatrix}
1 & 0 & -\frac{1}{\sqrt{2}} \\
0 & 1 & -\frac{r}{\sqrt{2}} \\
-\frac{1}{\sqrt{2}} & -\frac{r}{\sqrt{2}} & 1
\end{pmatrix}, \;\;\;
\Sigma_{3e}(r)=
\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & r \\
0 & r & 1
\end{pmatrix},  \;\;\;
\Sigma_{4a}(r)=
\begin{pmatrix}
1 & 0 & 0 & \frac{r}{\sqrt{2}} \\
0 & 1 & -r & \frac{r}{\sqrt{2}} \\
0 & -r & 1 & -\frac{1}{\sqrt{2}} \\
\frac{r}{\sqrt{2}} & \frac{r}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 1
\end{pmatrix},
$$

$$
\Sigma_{4b}(r)=
\begin{pmatrix}
1 & 0 & r & \frac{r}{\sqrt{2}} \\
0 & 1 & 0 & \frac{r}{\sqrt{2}} \\
r & 0 & 1 & \frac{1}{\sqrt{2}} \\
\frac{r}{\sqrt{2}} & \frac{r}{\sqrt{2}} & \frac{1}{\sqrt{2}} & 1
\end{pmatrix}, \;\;\;
\Sigma_{4c}(r)=
\begin{pmatrix}
1 & 0 & \frac{1}{\sqrt{2}} & -\frac{r}{\sqrt{2}} \\
0 & 1 & -\frac{r}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\
\frac{1}{\sqrt{2}} & -\frac{r}{\sqrt{2}} & 1 & -r \\
-\frac{r}{\sqrt{2}} & \frac{1}{\sqrt{2}} & -r & 1
\end{pmatrix}\;\;\text{and}\;\;
\Sigma_{4d}(r)=
\begin{pmatrix}
1 & r & \frac{1}{\sqrt{2}} & \frac{r}{\sqrt{2}} \\
r & 1 & \frac{r}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\
\frac{1}{\sqrt{2}} & \frac{r}{\sqrt{2}} & 1 & r \\
\frac{r}{\sqrt{2}} & \frac{1}{\sqrt{2}} & r & 1
\end{pmatrix}.
$$



# Estimation methods

Given the form of bridge function $F$, obtaining a moment-based estimation $\widehat \sigma_{jk}$ requires inversion of $F$. `latentcor` implements two methods for calculation of the inversion:

  * `method = "original"` [Subsection describing original method and relevant parameter `tol`](#original)
  * `method = "approx"` [Subsection describing approximation method and relevant parameter `ratio`](#approx)
  
Both methods calculate inverse bridge function applied to each element of sample Kendall's $\tau$ matrix. Because the calculation is performed point-wise (separately for each pair of variables), the resulting point-wise estimator of correlation matrix may not be positive semi-definite. `latentcor` performs projection of the pointwise-estimator to the space of positive semi-definite matrices, and allows for shrinkage towards identity matrix using the parameter `nu` (see [Subsection describing adjustment of point-wise estimator and relevant parameter `nu`](#shrinkage)).

## Original method (`method = "original"`) {#original }

Original estimation approach relies on numerical inversion of $F$ based on solving uni-root optimization problem. Given the calculated $\widehat \tau_{jk}$ (sample Kendall's $\tau$ between variables $j$ and $k$), the estimate of latent correlation $\widehat \sigma_{jk}$ is obtained by calling `optimize` function to solve the following optimization problem:
$$
\widehat r_{jk} = \arg\min_{r} \{F(r) - \widehat \tau_{jk}\}^2.
$$
The parameter `tol` controls the desired accuracy of the minimizer and is passed to `optimize`, with the default precision of 1e-8:

```{r callR, warning = FALSE, message = F}
estimate = latentcor(X, types = c("con", "bin", "ter", "tru"), method = "original", tol = 1e-8)
```

***Algorithm for Original method***

**Input**: $F(r)=F(r, \mathbf{\Delta})$ - bridge function based on the type of variables $j$, $k$

   - Step 1. Calculate $\hat{\tau}_{jk}$ using (1).
   
```{r kendall}
estimate$K
```
   
   - Step 2. For binary/truncated variable $j$, set $\hat{\mathbf{\Delta}}_{j}=\hat{\Delta}_{j}=\Phi^{-1}(\pi_{0j})$ with $\pi_{0j}=\sum_{i=1}^{n}\frac{I(x_{ij}=0)}{n}$. For ternary variable $j$, set $\hat{\mathbf{\Delta}}_{j}=(\hat{\Delta}_{j}^{1}, \hat{\Delta}_{j}^{2})$ where $\hat{\Delta}_{j}^{1}=\Phi^{-1}(\pi_{0j})$ and $\hat{\Delta}_{j}^{2}=\Phi^{-1}(\pi_{0j}+\pi_{1j})$ with $\pi_{0j}=\sum_{i=1}^{n}\frac{I(x_{ij}=0)}{n}$ and $\pi_{1j}=\sum_{i=1}^{n}\frac{I(x_{ij}=1)}{n}$.
   
```{r zratios_2}
estimate$zratios
```
   
   - Compute $F^{-1}(\hat{\tau}_{jk})$ as $\hat{r}_{jk}=argmin\{F(r)-\hat{\tau}_{jk}\}^{2}$ solved via `optimize` function in *R* with accuracy `tol`.

```{r estimate2}
estimate$Rpointwise
```


## Approximation method (`method = "approx"`) {#approx}

A faster approximation method is based on multi-linear interpolation of pre-computed inverse bridge function on a fixed grid of points [@yoon2021fast]. This is possible as the inverse bridge function is an analytic function of at most 5 parameters:

  - Kendall's $\tau$
  - Proportion of zeros in the 1st variable 
  - (Possibly) proportion of zeros and ones in the 1st variable
  - (Possibly) proportion of zeros in the 2nd variable
  - (Possibly) proportion of zeros and ones in the 2nd variable


In short, d-dimensional multi-linear interpolation uses a weighted average of $2^{d}$ neighbors to approximate the function values at the points within the d-dimensional cube of the neighbors, and to perform interpolation, `latentcor` takes advantage of the R package `chebpol` [@R-chebpol]. This approximation method has been first described in [@yoon2021fast] for continuous/binary/truncated cases. In `latentcor`, we additionally implement ternary case, and optimize the choice of grid as well as interpolation boundary for faster computations with smaller memory footprint.




```{r callR2, warning = FALSE, message = F}
estimate = latentcor(X, types = c("con", "bin", "ter", "tru"), method = "approx")
```


***Algorithm for Approximation method *** 

**Input**: Let $\check{g}=h(g)$, pre-computed values $F^{-1}(h^{-1}(\check{g}))$ on a fixed grid $\check{g}\in\check{\cal{G}}$ based on the type of variables $j$ and $k$. For binary/continuous case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j})$; for binary/binary case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j}, \check{\Delta}_{k})$; for truncated/continuous case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j})$; for truncated/truncated case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j}, \check{\Delta}_{k})$; for ternary/continuous case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j}^{1}, \check{\Delta}_{j}^{2})$; for ternary/binary case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j}^{1}, \check{\Delta}_{j}^{2}, \check{\Delta}_{k})$; for ternary/truncated case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j}^{1}, \check{\Delta}_{j}^{2}, \check{\Delta}_{k})$; for ternay/ternary case, $\check{g}=(\check{\tau}_{jk}, \check{\Delta}_{j}^{1}, \check{\Delta}_{j}^{2}, \check{\Delta}_{k}^{1}, \check{\Delta}_{k}^{2})$.

  - Step 1 and Step 2 same as Original method.
  
  - Step 3. If $|\hat{\tau}_{jk}|\le \mbox{ratio}\times \bar{\tau}_{jk}(\cdot)$, apply interpolation; otherwise apply Original method.



To avoid interpolation in areas with high approximation errors close to the boundary, we use hybrid scheme in Step 3. The parameter `ratio` controls the size of the region where the interpolation is performed (`ratio = 0` means no interpolation, `ratio = 1` means interpolation is always performed). For the derivation of approximate bound for BC, BB, TC, TB, TT cases see @yoon2021fast. The derivation of approximate bound for NC, NB, NN, NT case is in the Appendix.

$$
\bar{\tau}_{jk}(\cdot)=
\begin{cases}
2\pi_{0j}(1-\pi_{0j})  &   for \; BC \; case\\
2\min(\pi_{0j},\pi_{0k})\{1-\max(\pi_{0j}, \pi_{0k})\}  &   for \; BB \; case\\
2\{\pi_{0j}(1-\pi_{0j})+\pi_{1j}(1-\pi_{0j}-\pi_{1j})\}  &   for \; NC \; case\\
2\min(\pi_{0j}(1-\pi_{0j})+\pi_{1j}(1-\pi_{0j}-\pi_{1j}),\pi_{0k}(1-\pi_{0k}))  &   for \; NB \; case\\
2\min(\pi_{0j}(1-\pi_{0j})+\pi_{1j}(1-\pi_{0j}-\pi_{1j}), \\
\;\;\;\;\;\;\;\;\;\;\pi_{0k}(1-\pi_{0k})+\pi_{1k}(1-\pi_{0k}-\pi_{1k}))  &   for \; NN \; case\\
1-(\pi_{0j})^{2}  &   for \; TC \; case\\
2\max(\pi_{0k},1-\pi_{0k})\{1-\max(\pi_{0k},1-\pi_{0k},\pi_{0j})\}  &   for \; TB \; case\\
1-\{\max(\pi_{0j},\pi_{0k},\pi_{1k},1-\pi_{0k}-\pi_{1k})\}^{2}  &   for \; TN \; case\\
1-\{\max(\pi_{0j},\pi_{0k})\}^{2}  &   for \; TT \; case\\
\end{cases}
$$

By default, `latentcor` uses `ratio = 0.9` as this value was recommended in @yoon2021fast having a good balance of accuracy and computational speed. This value, however, can be modified by the user.

```{r estimate3, warning = FALSE, message = F}
latentcor(X, types = c("con", "bin", "ter", "tru"), method = "approx", ratio = 0.99)$R
latentcor(X, types = c("con", "bin", "ter", "tru"), method = "approx", ratio = 0.4)$R
latentcor(X, types = c("con", "bin", "ter", "tru"), method = "original")$R
```

The lower is the `ratio`, the closer is the approximation method to original method (with `ratio = 0` being equivalent to `method = "original"`), but also the higher is the cost of computations.

```{r speed, warning = FALSE, message = F}
library(microbenchmark)
microbenchmark(latentcor(X, types = c("con", "bin", "ter", "tru"), method = "approx", ratio = 0.99)$R)
microbenchmark(latentcor(X, types = c("con", "bin", "ter", "tru"), method = "approx", ratio = 0.4)$R)
microbenchmark(latentcor(X, types = c("con", "bin", "ter", "tru"), method = "original")$R)
```

**Rescaled Grid for Interpolation**

Since $|\hat{\tau}|\le \bar{\tau}$, the grid does not need to cover the whole domain $\tau\in[-1, 1]$. To optimize memory associated with storing the grid, we rescale $\tau$ as follows:
$\check{\tau}_{jk}=\tau_{jk}/\bar{\tau}_{jk}\in[-1, 1]$, where $\bar{\tau}_{jk}$ is as defined above. 

In addition, for ternary variable $j$, it always holds that $\Delta_{j}^{2}>\Delta_{j}^{1}$ since $\Delta_{j}^{1}=\Phi^{-1}(\pi_{0j})$ and $\Delta_{j}^{2}=\Phi^{-1}(\pi_{0j}+\pi_{1j})$. Thus, the grid should not cover the the area corresponding to $\Delta_{j}^{2}\ge\Delta_{j}^{1}$. We thus rescale as follows: $\check{\Delta}_{j}^{1}=\Delta_{j}^{1}/\Delta_{j}^{2}\in[0, 1]$; $\check{\Delta}_{j}^{2}=\Delta_{j}^{2}\in[0, 1]$.

**Speed Comparison**

To illustrate the speed improvement by `method = "approx"`, we plot the run time scaling behavior of `method = "approx"` and `method = "original"` (setting `types` for `gen_data` by replicating `c("con", "bin", "ter", "tru")` multiple times) with increasing dimensions $p = [20, 40, 100, 200, 400]$ at sample size $n = 100$ using simulation data. Figure below summarizes the observed scaling in a log-log plot. For both methods we observe the expected $O(p^2)$ scaling
behavior with dimension p, i.e., a linear scaling in the log-log plot. However, `method = "approx"` is at least one order of magnitude faster than `method = "original"` independent of the dimension of the problem.

![](./timing_plot.png)

## Adjustment of pointwise-estimator for positive-definiteness {#shrinkage}

Since the estimation is performed point-wise, the resulting matrix of estimated latent correlations is not guaranteed to be positive semi-definite. For example, this could be expected when the sample size is small (and so the estimation error for each pairwise correlation is larger)

```{r, message = FALSE}
set.seed("234820")
X = gen_data(n = 6, types = c("con", "bin", "ter", "tru"))$X
X
out = latentcor(X, types = c("con", "bin", "ter", "tru"))
out$Rpointwise
eigen(out$Rpointwise)$values
```

`latentcor` automatically corrects the pointwise estimator to be positive definite by making two adjustments. First, if `Rpointwise` has smallest eigenvalue less than zero, the `latentcor` projects this matrix to the nearest positive semi-definite matrix. The user is notified of this adjustment through the message (supressed in previous code chunk), e.g.
```{r, message = TRUE}
out = latentcor(X, types = c("con", "bin", "ter", "tru"))
```

Second, `latentcor` shrinks the adjusted matrix of correlations towards identity matrix using the parameter $\nu$ with default value of 0.001 (`nu = 0.001`), so that the resulting `R` is strictly positive definite with the minimal eigenvalue being greater or equal to $\nu$. That is
$$
R = (1 - \nu) \widetilde R + \nu I,
$$
where $\widetilde R$ is the nearest positive semi-definite matrix to `Rpointwise`. 
```{r}
out = latentcor(X, types = c("con", "bin", "ter", "tru"), nu = 0.001)
out$Rpointwise
out$R
```
As a result, `R` and `Rpointwise` could be quite different when sample size $n$ is small. When $n$ is large and $p$ is moderate, the difference is typically driven by parameter `nu`.

```{r}
set.seed("234820")
X = gen_data(n = 100, types = c("con", "bin", "ter", "tru"))$X
out = latentcor(X, types = c("con", "bin", "ter", "tru"), nu = 0.001)
out$Rpointwise
out$R
```

# Appendix

## Derivation of bridge function $F$ for ternary/truncated case

Without loss of generality, let $j=1$ and $k=2$. By the definition of Kendall's $\tau$,
$$
    \tau_{12}=E(\hat{\tau}_{12})=E[\frac{2}{n(n-1)}\sum_{1\leq i\leq i' \leq n} sign\{(X_{i1}-X_{i'1})(X_{i2}-X_{i'2})\}].
$$
Since $X_{1}$ is ternary,
\begin{align}
    &sign(X_{1}-X_{1}') \nonumber\\ =&[I(U_{1}>C_{11},U_{1}'\leq C_{11})+I(U_{1}>C_{12},U_{1}'\leq C_{12})-I(U_{1}>C_{12},U_{1}'\leq C_{11})] \nonumber\\
    &-[I(U_{1}\leq C_{11}, U_{1}'>C_{11})+I(U_{1}\leq C_{12}, U_{1}'>C_{12})-I(U_{1}\leq C_{11}, U_{1}'>C_{12})] \nonumber\\
    =&[I(U_{1}>C_{11})-I(U_{1}>C_{11},U_{1}'>C_{11})+I(U_{1}>C_{12})-I(U_{1}>C_{12},U_{1}'>C_{12}) \nonumber\\
    &-I(U_{1}>C_{12})+I(U_{1}>C_{12},U_{1}'>C_{11})] \nonumber\\
    &-[I(U_{1}'>C_{11})-I(U_{1}>C_{11},U_{1}'>C_{11})+I(U_{1}'>C_{12})-I(U_{1}>C_{12},U_{1}'>C_{12}) \nonumber\\
    &-I(U_{1}'>C_{12})+I(U_{1}>C_{11},U_{1}'>C_{12})] \nonumber\\
    =&I(U_{1}>C_{11})+I(U_{1}>C_{12},U_{1}'>C_{11})-I(U_{1}'>C_{11})-I(U_{1}>C_{11},U_{1}'>C_{12}) \nonumber\\
    =&I(U_{1}>C_{11},U_{1}'\leq C_{12})-I(U_{1}'>C_{11},U_{1}\leq C_{12}).
\end{align}
Since $X_{2}$ is truncated, $C_{1}>0$ and
\begin{align}
    sign(X_{2}-X_{2}')=&-I(X_{2}=0,X_{2}'>0)+I(X_{2}>0,X_{2}'=0) \nonumber\\
    &+I(X_{2}>0,X_{2}'>0)sign(X_{2}-X_{2}') \nonumber\\
    =&-I(X_{2}=0)+I(X_{2}'=0)+I(X_{2}>0,X_{2}'>0)sign(X_{2}-X_{2}').
\end{align}
Since $f$ is monotonically increasing, $sign(X_{2}-X_{2}')=sign(Z_{2}-Z_{2}')$,
\begin{align}
    \tau_{12}=&E[I(U_{1}>C_{11},U_{1}'\leq C_{12}) sign(X_{2}-X_{2}')] \nonumber\\ &-E[I(U_{1}'>C_{11},U_{1}\leq C_{12}) sign(X_{2}-X_{2}')] \nonumber\\
    =&-E[I(U_{1}>C_{11},U_{1}'\leq C_{12}) I(X_{2}=0)] \nonumber\\
    &+E[I(U_{1}>C_{11},U_{1}'\leq C_{12}) I(X_{2}'=0)] \nonumber\\
    &+E[I(U_{1}>C_{11},U_{1}'\leq C_{12})I(X_{2}>0,X_{2}'>0)sign(Z_{2}-Z_{2}')] \nonumber\\
    &+E[I(U_{1}'>C_{11},U_{1}\leq C_{12}) I(X_{2}=0)] \nonumber\\
    &-E[I(U_{1}'>C_{11},U_{1}\leq C_{12}) I(X_{2}'=0)] \nonumber\\
    &-E[I(U_{1}'>C_{11},U_{1}\leq C_{12})I(X_{2}>0,X_{2}'>0)sign(Z_{2}-Z_{2}')]  \nonumber\\
    =&-2E[I(U_{1}>C_{11},U_{1}'\leq C_{12}) I(X_{2}=0)] \nonumber\\
    &+2E[I(U_{1}>C_{11},U_{1}'\leq C_{12}) I(X_{2}'=0)] \nonumber\\
    &+E[I(U_{1}>C_{11},U_{1}'\leq C_{12})I(X_{2}>0,X_{2}'>0)sign(Z_{2}-Z_{2}')] \nonumber\\
    &-E[I(U_{1}'>C_{11},U_{1}\leq C_{12})I(X_{2}>0,X_{2}'>0)sign(Z_{2}-Z_{2}')].
\end{align}
From the definition of $U$, let $Z_{j}=f_{j}(U_{j})$ and $\Delta_{j}=f_{j}(C_{j})$ for $j=1,2$. Using $sign(x)=2I(x>0)-1$, we obtain
\begin{align}
    \tau_{12}=&-2E[I(Z_{1}>\Delta_{11},Z_{1}'\leq \Delta_{12},Z_{2}\leq \Delta_{2})]+2E[I(Z_{1}>\Delta_{11},Z_{1}'\leq \Delta_{12},Z_{2}'\leq \Delta_{2})] \nonumber\\
    &+2E[I(Z_{1}>\Delta_{11},Z_{1}'\leq \Delta_{12})I(Z_{2}>\Delta_{2},Z_{2}'>\Delta_{2},Z_{2}-Z_{2}'>0)] \nonumber\\
    &-2E[I(Z_{1}'>\Delta_{11},Z_{1}\leq \Delta_{12})I(Z_{2}>\Delta_{2},Z_{2}'>\Delta_{2},Z_{2}-Z_{2}'>0)] \nonumber\\
    =&-2E[I(Z_{1}>\Delta_{11},Z_{1}'\leq \Delta_{12}, Z_{2}\leq \Delta_{2})]+2E[I(Z_{1}>\Delta_{11},Z_{1}'\leq \Delta_{12}, Z_{2}'\leq \Delta_{2})] \nonumber\\
    &+2E[I(Z_{1}>\Delta_{11},Z_{1}'\leq\Delta_{12},Z_{2}'>\Delta_{2},Z_{2}>Z_{2}')] \nonumber\\
    &-2E[I(Z_{1}'>\Delta_{11},Z_{1}\leq\Delta_{12},Z_{2}'>\Delta_{2},Z_{2}>Z_{2}')].
\end{align}
Since $\{\frac{Z_{2}'-Z_{2}}{\sqrt{2}}, -Z{1}\}$, $\{\frac{Z_{2}'-Z_{2}}{\sqrt{2}}, Z{1}'\}$ and $\{\frac{Z_{2}'-Z_{2}}{\sqrt{2}}, -Z{2}'\}$ are standard bivariate normally distributed variables with correlation $-\frac{1}{\sqrt{2}}$, $r/\sqrt{2}$ and $-\frac{r}{\sqrt{2}}$, respectively, by the definition of $\Phi_3(\cdot,\cdot, \cdot;\cdot)$ and $\Phi_4(\cdot,\cdot, \cdot,\cdot;\cdot)$ we have
\begin{align}
    F_{NT}(r;\Delta_{j}^{1},\Delta_{j}^{2},\Delta_{k})= & -2\Phi_{3}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},\Delta_{k};\begin{pmatrix}
1 & 0 & -r \\
0 & 1 & 0 \\
-r & 0 & 1
\end{pmatrix} \right\} \nonumber\\
    &+2\Phi_{3}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},\Delta_{k};\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & r \\
0 & r & 1
\end{pmatrix}\right\}\nonumber \\
    & +2\Phi_{4}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k},0;\begin{pmatrix}
1 & 0 & 0 & \frac{r}{\sqrt{2}} \\
0 & 1 & -r & \frac{r}{\sqrt{2}} \\
0 & -r & 1 & -\frac{1}{\sqrt{2}} \\
\frac{r}{\sqrt{2}} & \frac{r}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 1
\end{pmatrix}\right\} \nonumber\\
    &-2\Phi_{4}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k},0;\begin{pmatrix}
1 & 0 & r & -\frac{r}{\sqrt{2}} \\
0 & 1 & 0 & -\frac{r}{\sqrt{2}} \\
r & 0 & 1 & -\frac{1}{\sqrt{2}} \\
-\frac{r}{\sqrt{2}} & -\frac{r}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 1
\end{pmatrix}\right\}.
\end{align}
Using the facts that
\begin{align}
&\Phi_{4}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k},0;\begin{pmatrix}
1 & 0 & r & -\frac{r}{\sqrt{2}} \\
0 & 1 & 0 & -\frac{r}{\sqrt{2}} \\
r & 0 & 1 & -\frac{1}{\sqrt{2}} \\
-\frac{r}{\sqrt{2}} & -\frac{r}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 1
\end{pmatrix}\right\} \nonumber\\ &+\Phi_{4}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k},0;\begin{pmatrix}
1 & 0 & r & \frac{r}{\sqrt{2}} \\
0 & 1 & 0 & \frac{r}{\sqrt{2}} \\
r & 0 & 1 & \frac{1}{\sqrt{2}} \\
\frac{r}{\sqrt{2}} & \frac{r}{\sqrt{2}} & \frac{1}{\sqrt{2}} & 1
\end{pmatrix}\right\} \nonumber\\
=&\Phi_{3}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k};\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & r \\
0 & r & 1
\end{pmatrix}\right\}
\end{align}
and
\begin{align}
&\Phi_{3}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k};\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & r \\
0 & r & 1
\end{pmatrix}\right\}+\Phi_{3}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},\Delta_{k};\begin{pmatrix}
1 & 0 & -r \\
0 & 1 & 0 \\
-r & 0 & 1
\end{pmatrix} \right\} \nonumber\\
=&\Phi_{2}(-\Delta_{j}^{1},\Delta_{j}^{2};0)
=\Phi(-\Delta_{j}^{1})\Phi(\Delta_{j}^{2}).
\end{align}
So that,
\begin{align}
    F_{NT}(r;\Delta_{j}^{1},\Delta_{j}^{2},\Delta_{k})= & -2\Phi(-\Delta_{j}^{1})\Phi(\Delta_{j}^{2}) \nonumber\\
    &+2\Phi_{3}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},\Delta_{k};\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & r \\
0 & r & 1
\end{pmatrix}\right\}\nonumber \\
    & +2\Phi_{4}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k},0;\begin{pmatrix}
1 & 0 & 0 & \frac{r}{\sqrt{2}} \\
0 & 1 & -r & \frac{r}{\sqrt{2}} \\
0 & -r & 1 & -\frac{1}{\sqrt{2}} \\
\frac{r}{\sqrt{2}} & \frac{r}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 1
\end{pmatrix}\right\} \nonumber\\
    &+2\Phi_{4}\left\{-\Delta_{j}^{1},\Delta_{j}^{2},-\Delta_{k},0;\begin{pmatrix}
1 & 0 & r & \frac{r}{\sqrt{2}} \\
0 & 1 & 0 & \frac{r}{\sqrt{2}} \\
r & 0 & 1 & \frac{1}{\sqrt{2}} \\
\frac{r}{\sqrt{2}} & \frac{r}{\sqrt{2}} & \frac{1}{\sqrt{2}} & 1
\end{pmatrix}\right\}.
\end{align}

It is easy to get the bridge function for truncated/ternary case by switching $j$ and $k$.

## Derivation of approximate bound for the ternary/continuous case

Let $n_{0x}=\sum_{i=1}^{n_x}I(x_{i}=0)$, $n_{2x}=\sum_{i=1}^{n_x}I(x_{i}=2)$, $\pi_{0x}=\frac{n_{0x}}{n_{x}}$ and $\pi_{2x}=\frac{n_{2x}}{n_{x}}$, then
\begin{align}
    |\tau(\mathbf{x})|\leq & \frac{n_{0x}(n-n_{0x})+n_{2x}(n-n_{0x}-n_{2x})}{\begin{pmatrix} n \\ 2 \end{pmatrix}} \nonumber\\
    = & 2\{\frac{n_{0x}}{n-1}-(\frac{n_{0x}}{n})(\frac{n_{0x}}{n-1})+\frac{n_{2x}}{n-1}-(\frac{n_{2x}}{n})(\frac{n_{0x}}{n-1})-(\frac{n_{2x}}{n})(\frac{n_{2x}}{n-1})\} \nonumber\\
    \approx & 2\{\frac{n_{0x}}{n}-(\frac{n_{0x}}{n})^2+\frac{n_{2x}}{n}-(\frac{n_{2x}}{n})(\frac{n_{0x}}{n})-(\frac{n_{2x}}{n})^2\} \nonumber\\
    = & 2\{\pi_{0x}(1-\pi_{0x})+\pi_{2x}(1-\pi_{0x}-\pi_{2x})\}
\end{align}

For ternary/binary and ternary/ternary cases, we combine the two individual bounds.


## Derivation of approximate bound for the ternary/truncated case

 Let $\mathbf{x}\in\mathcal{R}^{n}$ and $\mathbf{y}\in\mathcal{R}^{n}$ be the observed $n$ realizations of ternary and truncated variables, respectively. Let $n_{0x}=\sum_{i=0}^{n}I(x_{i}=0)$, $\pi_{0x}=\frac{n_{0x}}{n}$, $n_{1x}=\sum_{i=0}^{n}I(x_{i}=1)$, $\pi_{1x}=\frac{n_{1x}}{n}$, $n_{2x}=\sum_{i=0}^{n}I(x_{i}=2)$, $\pi_{2x}=\frac{n_{2x}}{n}$,
$n_{0y}=\sum_{i=0}^{n}I(y_{i}=0)$, $\pi_{0y}=\frac{n_{0y}}{n}$, $n_{0x0y}=\sum_{i=0}^{n}I(x_{i}=0 \;\& \; y_{i}=0)$, $n_{1x0y}=\sum_{i=0}^{n}I(x_{i}=1 \;\& \; y_{i}=0)$ and
$n_{2x0y}=\sum_{i=0}^{n}I(x_{i}=2 \;\& \; y_{i}=0)$ then
\begin{align}
    |\tau(\mathbf{x}, \mathbf{y})|\leq &
    \frac{\begin{pmatrix}n \\ 2\end{pmatrix}-\begin{pmatrix}n_{0x} \\ 2\end{pmatrix}-\begin{pmatrix}n_{1x} \\ 2\end{pmatrix}-\begin{pmatrix} n_{2x} \\ 2 \end{pmatrix}-\begin{pmatrix}n_{0y} \\ 2\end{pmatrix}+\begin{pmatrix}n_{0x0y} \\ 2 \end{pmatrix}+\begin{pmatrix}n_{1x0y} \\ 2\end{pmatrix}+\begin{pmatrix}n_{2x0y} \\ 2\end{pmatrix}}{\begin{pmatrix}n \\ 2\end{pmatrix}} \nonumber
\end{align}
Since $n_{0x0y}\leq\min(n_{0x},n_{0y})$, $n_{1x0y}\leq\min(n_{1x},n_{0y})$ and $n_{2x0y}\leq\min(n_{2x},n_{0y})$ we obtain
\begin{align}
     |\tau(\mathbf{x}, \mathbf{y})|\leq &
    \frac{\begin{pmatrix}n \\ 2\end{pmatrix}-\begin{pmatrix}n_{0x} \\ 2\end{pmatrix}-\begin{pmatrix}n_{1x} \\ 2\end{pmatrix}-\begin{pmatrix} n_{2x} \\ 2 \end{pmatrix}-\begin{pmatrix}n_{0y} \\ 2\end{pmatrix}}{\begin{pmatrix}n \\ 2\end{pmatrix}} \nonumber\\
    & +  \frac{\begin{pmatrix}\min(n_{0x},n_{0y}) \\ 2 \end{pmatrix}+\begin{pmatrix}\min(n_{1x},n_{0y}) \\ 2\end{pmatrix}+\begin{pmatrix}\min(n_{2x},n_{0y}) \\ 2\end{pmatrix}}{\begin{pmatrix}n \\ 2\end{pmatrix}} \nonumber\\
    \leq & \frac{\begin{pmatrix}n \\ 2\end{pmatrix}-\begin{pmatrix}\max(n_{0x},n_{1x},n_{2x},n_{0y}) \\ 2\end{pmatrix}}{\begin{pmatrix}n \\ 2\end{pmatrix}} \nonumber\\
    \leq & 1-\frac{\max(n_{0x},n_{1x},n_{2x},n_{0y})(\max(n_{0x},n_{1x},n_{2x},n_{0y})-1)}{n(n-1)} \nonumber\\
    \approx & 1-(\frac{\max(n_{0x},n_{1x},n_{2x},n_{0y})}{n})^{2} \nonumber\\
    =& 1-\{\max(\pi_{0x},\pi_{1x},\pi_{2x},\pi_{0y})\}^{2} \nonumber\\
    =& 1-\{\max(\pi_{0x},(1-\pi_{0x}-\pi_{2x}),\pi_{2x},\pi_{0y})\}^{2}
\end{align}

It is easy to get the approximate bound for truncated/ternary case by switching $\mathbf{x}$ and $\mathbf{y}$.

# References


---
title: "latentcor"
author: "Mingze Huang, Christian L. Müller, Irina Gaynanova"
date: "`r Sys.Date()`"
bibliography: latentcor.bib
output: rmarkdown::html_vignette
extra_dependencies: ["amsmath"]
nocite: |
  @croux2013robust
  @filzmoser2021pcapp
  @liu2009nonparanormal
  @fox2019poly
vignette: >
  %\VignetteIndexEntry{latentcor}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(latentcor)
```


# Introduction

R package `latentcor` utilizes the powerful semi-parametric latent Gaussian copula models to estimate latent correlations between mixed data types. The package allows to estimate correlations between any of continuous/binary/ternary/zero-inflated (truncated) variable types. The underlying implementation takes advantage of fast multi-linear interpolation scheme with a clever choice of grid points that give the package a small memory footprint, and allows to use the latent correlations with sub-sampling and bootstrapping.

# Statement of need

No R software package is currently available that allows accurate and fast correlation estimation from mixed variable data in a unifying manner. The R package `latentcor`, introduced here, thus represents the first stand-alone R package for 
computation of latent correlation that takes into account all variable types (continuous/binary/ordinal/zero-inflated), comes with an optimized memory footprint, 
and is computationally efficient, essentially making latent correlation estimation almost as fast as rank-based correlation estimation. 

# Getting started

## A simple example with two variables

First, we will generate a pair of variables with different types using a sample size $n=100$ which will serve as example data. Here first variable will be ternary, and second variable will be continuous.

```{r data_generation}
simdata = gen_data(n = 100, types = c("ter", "con"))
```

The output of `gen_data` is a list with 2 elements:

```{r data_output}
names(simdata)
```

  - `X`: a matrix ($100\times 2$), the first column is the ternary variable; the second column is the continuous variable.
  
```{r data_matrix}
X = simdata$X
head(X, n = 6L)
```

  - `plotX`: NULL (`showplot = FALSE`, can be changed to display the plot of generated data in`gen_data` input).

```{r data_plot}
simdata$plotX
```

Then we can estimate the latent correlation matrix based on these 2 variables using `latentcor` function.

```{r estimation}
estimate = latentcor(X, types = c("ter", "con"))
```

The output of `latentcor` is a list with several elements:

```{r estimation_output}
names(estimate)
```

  - `zratios` is a list has the same length as the number of variables. Here the first element is a ($2\times1$) vector indicating the cumulative proportions for zeros and ones in the ternary variable (e.g. first element in vector is the proportion of zeros, second element in vector is the proportion of zeros and ones.) The second element of the list is NA for continuous variable.

```{r zratios}
estimate$zratios
```

  - `K`: Kendall $\tau$ ($\tau_{a}$) correlation matrix for these 2 variables. 
  
```{r Kendall}
estimate$K
```  

  - `Rpointwise`: matrix of pointwise estimated correlations. Due to pointwise estimation, `Rpointwise` is not guaranteed to be positive semi-definite

```{r latent_correlation_pointwise}
estimate$Rpointwise
``` 

  - `R`: estimated final latent correlation matrix, this matrix is guaranteed to be strictly positive definite (through `nearPD` projection and parameter `nu`, see Mathematical framework for estimation) if `use.nearPD = TRUE`.

```{r latent_correlation}
estimate$R
``` 

  - `plotR`: NULL by default as `showplot = FALSE` in `latentcor`. Otherwise displays a heatmap of latent correlation matrix.

```{r heatmap}
estimate$plotR
``` 

## Example with mtcars dataset

We use the build-in dataset `mtcars`:

```{r mtcars}
head(mtcars, n = 6L)
```

Let's take a look at the unique values for each variable to determine the corresponding data type.

```{r unique}
apply(mtcars, 2, table)
```


Then we can estimate the latent correlation matrix for all variables of `mtcars` by using `latentcor` function.

```{r mtcars_estimation}
estimate_mtcars = latentcor(mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin", "bin", "ter", "con"))
```

Note that the determination of variable types can also be done automatically by `latentcor` package using `get_types` function:
```{r mtcars_types, message = FALSE}
get_types(mtcars)
```
This function is run automatically inside `latentcor` if the `types` are not supplied by the user, however the automatic determination of types takes extra time, so we recommend to specify `types` explicitly if they are known in advance.

The output of `latentcor` for `mtcars`:

```{r mtcars_estimation_output}
names(estimate_mtcars)
```

  - `zratios`: zratios for corresponding variables in `mtcars`.

```{r mtcars_zratios}
estimate_mtcars$zratios
```

  - `K`: Kendall $\tau$ ($\tau_{a}$) correlation matrix for variables in `mtcars`. 
  
```{r mtcars_Kendall}
estimate_mtcars$K
```  

  - `Rpointwise`: matrix of pointwise estimated correlations for `mtcars`.

```{r mtcars_latent_correlation_pointwise}
estimate_mtcars$Rpointwise
``` 

  - `R`: estimated final latent correlation matrix for `mtcars`.

```{r mtcars_latent_correlation}
estimate_mtcars$R
``` 

  - `plotR`: NULL by default as `showplot = FALSE` in `latentcor`. Otherwise displays a heatmap of latent correlation matrix for `mtcars` (See [heatmap of latent correlation (approx) for mtcars](https://rpubs.com/mingzehuang/797937)).

```{r mtcars_heatmap}
estimate_mtcars$plotR
``` 

## Example using latentcor with subsampling

While `latentcor` can determine the types of each variable automatically, it is recommended to call `get_types` first and then supply `types` explicitly to save the computation time, especially when using latentcor with sub-sampling (which we illustrate below).

First, we will generate variables with different types using a sample size $n=100$ which will serve as an example data for subsampling. 

```{r data_generation 2}
simdata2 = gen_data(n = 100, types = c(rep("ter", 3), "con", rep("bin", 3)))
```

To use the data with subsampling, we recommend to first run `get_types` on the full data
```{r types subsampling}
types = get_types(simdata2$X)
types
```

Then, when doing subsampling, we recommend to explicitly supply identified types to `latentcor`. We illustrate using 10 subsamples, each of size 80.
```{r subsampling}
start_time = proc.time()
for (s in 1:10){
  # Select a random subsample of size 80
  subsample = sample(1:100, 80)
  # Estimate latent correlation on subsample specifying the types
  Rs = latentcor(simdata2$X[subsample, ], types = types)
}
proc.time() - start_time
```
Compared with
```{r subsampling 2}
start_time = proc.time()
for (s in 1:10){
  # Select a random subsample of size 80
  subsample = sample(1:100, 80)
  # Estimate latent correlation on subsample specifying the types
  Rs = latentcor(simdata2$X[subsample, ], types = NULL)
}
proc.time() - start_time
```

# References


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_types.R
\name{get_types}
\alias{get_types}
\title{Automatically determine types of each variable (continuous/binary/ternary/truncated) in a data matrix.}
\usage{
get_types(X, tru_prop = 0.05)
}
\arguments{
\item{X}{A numeric data matrix (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed.}

\item{tru_prop}{A scalar between 0 and 1 indicating the minimal proportion of zeros that should be present in a variable to be treated as \code{"tru"} (truncated type or zero-inflated) rather than as \code{"con"} (continuous type). The default value is 0.05 (any variable with more than 5\% of zero values among n samples is treated as truncated or zero-inflated)}
}
\value{
\code{get_types} returns
\itemize{
      \item{types: }{A vector of length p indicating the type of each of the p variables in \code{X}. Each element is one of \code{"con"} (continuous), \code{"bin"} (binary), \code{"ter"} (ternary) or \code{"tru"} (truncated).}
}
}
\description{
Automatically determine types of each variable (continuous/binary/ternary/truncated) in a data matrix.
}
\examples{
X = gen_data(types = c("ter", "con"))$X
get_types(X)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/latentcor.R
\name{latentcor}
\alias{latentcor}
\title{Estimate latent correlation for mixed types.}
\usage{
latentcor(
  X,
  types = NULL,
  method = c("approx", "original"),
  use.nearPD = TRUE,
  nu = 0.001,
  tol = 1e-08,
  ratio = 0.9,
  showplot = FALSE
)
}
\arguments{
\item{X}{A numeric matrix or numeric data frame (n by p), where n is number of samples, and p is number of variables. Missing values (NA) are allowed, in which case the estimation is based on pairwise complete observations.}

\item{types}{A vector of length p indicating the type of each of the p variables in \code{X}. Each element must be one of \code{"con"} (continuous), \code{"bin"} (binary), \code{"ter"} (ternary) or \code{"tru"} (truncated). If the vector has length 1, then all p variables are assumed to be of the same type that is supplied. The default value is \code{NULL}, and the variable types are determined automatically using function \code{\link{get_types}}. As automatic determination of variable types takes extra time, it is recommended to supply the types explicitly when they are known in advance.}

\item{method}{The calculation method for latent correlations. Either \code{"original"} or \code{"approx"}. If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method, see Yoon et al. (2021) for timing comparisons for various variable types. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is \code{"approx"}.}

\item{use.nearPD}{Logical indicator. \code{use.nearPD = TRUE} gets nearest positive definite matrix for the estimated latent correlation matrix with shrinkage adjustment by \code{nu}. Output \code{R} is the same as \code{Rpointwise} if \code{use.nearPD = FALSE}. Default value is \code{TRUE}.}

\item{nu}{Shrinkage parameter for the correlation matrix, must be between 0 and 1. Guarantees that the minimal eigenvalue of returned correlation matrix is greater or equal to \code{nu}. When \code{nu = 0}, no shrinkage is performed, the returned correlation matrix will be semi-positive definite but not necessarily strictly positive definite. When \code{nu = 1}, the identity matrix is returned (not recommended).  The default (recommended) value is 0.001.}

\item{tol}{When \code{method = "original"}, specifies the desired accuracy of the bridge function inversion via uniroot optimization and is passed to \code{\link{optimize}}. The default value is 1e-8. When \code{method = "approx"}, this parameter is ignored.}

\item{ratio}{When \code{method = "approx"}, specifies the boundary value for multilinear interpolation, must be between 0 and 1. When \code{ratio = 0}, no linear interpolation is performed (the slowest execution) which is equivalent to \code{method = "original"}. When \code{ratio = 1}, linear interpolation is always performed (the fastest execution) but may lead to high approximation errors. The default (recommended) value is 0.9 which controls the approximation error and has fast execution, see Yoon et al. (2021) for details. When \code{method = "original"}, this parameter is ignored.}

\item{showplot}{Logical indicator. \code{showplot = TRUE} generates a ggplot object \code{plotR} with the heatmap of latent correlation matrix \code{R}. \code{plotR = NULL} if \code{showplot = FALSE}. Default value is \code{FALSE}.}
}
\value{
\code{latentcor} returns
\itemize{
      \item{zratios: }{A list of of length p corresponding to each variable. Returns NA for continuous variable; proportion of zeros for binary/truncated variables; the cumulative proportions of zeros and ones (e.g. first value is proportion of zeros, second value is proportion of zeros and ones) for ternary variable. }
      \item{K: }{(p x p) Kendall Tau (Tau-a) Matrix for \code{X} }
      \item{R: }{(p x p) Estimated latent correlation matrix for \code{X} }
      \item{Rpointwise: }{(p x p) Point-wise estimates of latent correlations for \code{X}. This matrix is not guaranteed to be semi-positive definite. This is the original estimated latent correlation matrix without adjustment for positive-definiteness.}
      \item{plotR: }{Heatmap plot of latent correlation matrix \code{R}, NULL if \code{showplot = FALSE}}
}
}
\description{
Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/binary/truncated/ternary) based on the latent Gaussian copula model. Missing values (NA) are allowed. The estimation is based on pairwise complete observations.
}
\details{
The function estimates latent correlation by calculating inverse bridge function (Fan et al., 2017) evaluated at the value of sample Kendall's tau (\eqn{\hat \tau}). The bridge function F connects Kendall's tau to latent correlation r so that \eqn{F(r) = E(\hat \tau)}. The form of function F depends on the variable types (continuous/binary/truncated/ternary), but is exact. The exact form of inverse is not available, so has to be evaluated numerically for each pair of variables leading to \code{Rpointwise}.

When \code{method = "original"}, the inversion is done numerically by solving \deqn{minimize_r (F(r) - \hat \tau)^2} using \code{\link{optimize}}. The parameter \code{tol} is used to control the accuracy of the solution.

When \code{method = "approx"}, the inversion is done approximately by interpolating previously calculated and stored values of \eqn{F^{-1}(\hat \tau)}. This is significantly faster than the original method (Yoon et al., 2021) for binary/ternary/truncated cases, however the approximation errors may be non-negligible on some regions of the space. The parameter \code{ratio} controls the region where the interpolation is performed with default recommended value of 0.9 giving a good balance of accuracy and computational speed . Increasing the value of ratio may improve speed (but possibly sacrifice the accuracy), whereas decreasing the value of ratio my improve accuracy (but possibly sacrifice the speed). See Yoon et al. 2021 and vignette for more details.

 In case the pointwise estimator \code{Rpointwise} is has negative eigenvalues, it is projected onto the space of positive semi-definite matrices using \code{\link{nearPD}}. The parameter \code{nu} further allows to perform additional shrinkage towards identity matrix (desirable in cases where the number of variables p is very large) as
 \deqn{R = (1 - \nu) \tilde R + \nu I,}
 where \eqn{\tilde R} is \code{Rpointwise} after projection by \code{\link{nearPD}}.
}
\examples{
# Example 1 - truncated data type, same type for all variables
# Generate data
X = gen_data(n = 300, types = rep("tru", 5))$X
# Estimate latent correlation matrix with original method and check the timing
start_time = proc.time()
R_nc_org = latentcor(X = X, types = "tru", method = "original")$R
proc.time() - start_time
# Estimate latent correlation matrix with approximation method and check the timing
start_time = proc.time()
R_nc_approx = latentcor(X = X, types = "tru", method = "approx")$R
proc.time() - start_time
# Heatmap for latent correlation matrix.
Heatmap_R_nc_approx = latentcor(X = X, types = "tru", method = "approx",
                                showplot = TRUE)$plotR

# Example 2 - ternary/continuous case
X = gen_data()$X
# Estimate latent correlation matrix with original method
R_nc_org = latentcor(X = X, types = c("ter", "con"), method = "original")$R
# Estimate latent correlation matrix with aprroximation method
R_nc_approx = latentcor(X = X, types = c("ter", "con"), method = "approx")$R
}
\references{
Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphical model for mixed data" \doi{10.1111/rssb.12168}.

Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" \doi{10.1093/biomet/asaa007}.

Yoon G., Müller C.L., Gaynanova I. (2021) "Fast computation of latent correlations" \doi{10.1080/10618600.2021.1882468}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_data.R
\name{gen_data}
\alias{gen_data}
\title{Mixed type simulation data generator}
\usage{
gen_data(
  n = 100,
  types = c("ter", "con"),
  rhos = 0.5,
  copulas = "no",
  XP = NULL,
  showplot = FALSE
)
}
\arguments{
\item{n}{A positive integer indicating the sample size. The default value is 100.}

\item{types}{A vector indicating the type of each variable, could be \code{"con"} (continuous), \code{"bin"} (binary), \code{"tru"} (truncated) or \code{"ter"} (ternary). The number of variables is determined based on the length of types, that is \code{p = length(types)}. The default value is \code{c("ter", "con")} which creates two variables: the first one is ternary, the second one is continuous.}

\item{rhos}{A vector with lower-triangular elements of desired correlation matrix, e.g. \code{rhos = c(.3, .5, .7)} means the correlation matrix is \code{matrix(c(1, .3, .5, .3, 1, .7, .5, .7, 1), 3, 3)}. If only a scalar is supplied (\code{length(rhos) = 1}), then equi-correlation matrix is assumed with all pairwise correlations being equal to \code{rhos}. The default value is 0.5 which means correlations between any two variables are 0.5.}

\item{copulas}{A vector indicating the copula tramsformation f for each of the p variables, e.g. U = f(Z). Each element can take value \code{"no"} (f is identity), \code{"expo"} (exponential transformation) or \code{"cube"} (cubic transformation). If the vector has length 1, then the same transformation is applied to all p variables. The default value is \code{"no"}: no copula transformation for any of the variables.}

\item{XP}{A list of length p indicating proportion of zeros (for binary and truncated), and proportions of zeros and ones (for ternary) for each of the variables. For continuous variable, NA should be supplied. If \code{NULL}, the following values are automatically generated as elements of \code{XP} list for the corresponding data types:
For continuous variable, the corresponding value is NA;
for binary or truncated variable, the corresponding value is a number between 0 and 1 representing the proportion of zeros, the default value is 0.5;
for ternary variable, the corresponding value is a pair of numbers between 0 and 1, the first number indicates the the proportion of zeros, the second number indicates the proportion of ones. The sum of a pair of numbers should be between 0 and 1, the default value is \code{c(0.3, 0.5)}.}

\item{showplot}{Logical indicator. If TRUE, generates the plot of the data when number of variables p is no more than 3. The default value is FALSE.}
}
\value{
\code{gen_data} returns a list containing
\itemize{
      \item{X: }{Generated data matrix (n by p) of observed variables.}
      \item{plotX: }{Visualization of the data matrix X.
                    Histogram if \code{p=1}. 2D Scatter plot if \code{p=2}. 3D scatter plot if \code{p=3}. Returns NULL if \code{showplot = FALSE}.}
}
}
\description{
Generates data of mixed types from the latent Gaussian copula model.
}
\examples{
# Generate single continuous variable with exponential transformation (always greater than 0)
# and show histogram.
simdata = gen_data(n = 100, copulas = "expo", types = "con", showplot = FALSE)
X = simdata$X; plotX = simdata$plotX
# Generate a pair of variables (ternary and continuous) with default proportions
# and without copula transformation.
simdata = gen_data()
X = simdata$X
# Generate 3 variables (binary, ternary and truncated)
# corresponding copulas for each variables are "no" (no transformation),
# "cube" (cube transformation) and "cube" (cube transformation).
# binary variable has 30\% of zeros, ternary variable has 20\% of zeros
# and 40\% of ones, truncated variable has 50\% of zeros.
# Then show the 3D scatter plot (data points project on either 0 or 1 on Axis X1;
# on 0, 1 or 2 on Axas X2; on positive domain on Axis X3)
simdata = gen_data(n = 100, rhos = c(.3, .4, .5), copulas = c("no", "cube", "cube"),
          types = c("bin", "ter", "tru"), XP = list(.3, c(.2, .4), .5), showplot = TRUE)
X = simdata$X; plotX = simdata$plotX
# Check the proportion of zeros for the binary variable.
sum(simdata$X[ , 1] == 0)
# Check the proportion of zeros and ones for the ternary variable.
sum(simdata$X[ , 2] == 0); sum(simdata$X[ , 2] == 1)
# Check the proportion of zeros for the truncated variable.
sum(simdata$X[ , 3] == 0)
}
\references{
Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" \doi{10.1111/rssb.12168}.

Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" \doi{10.1093/biomet/asaa007}.
}
