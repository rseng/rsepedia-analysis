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
