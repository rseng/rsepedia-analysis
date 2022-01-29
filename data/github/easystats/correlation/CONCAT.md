
# correlation <img src='man/figures/logo.png' align="right" height="139" />

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02306/status.svg)](https://doi.org/10.21105/joss.02306)
[![downloads](http://cranlogs.r-pkg.org/badges/correlation)](https://cran.r-project.org/package=correlation)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/correlation)](https://cranlogs.r-pkg.org/)
[![status](https://tinyverse.netlify.com/badge/correlation)](https://CRAN.R-project.org/package=correlation)

`correlation` is an
[**easystats**](https://github.com/easystats/easystats) package focused
on correlation analysis. It’s lightweight, easy to use, and allows for
the computation of many different kinds of correlations, such as
**partial** correlations, **Bayesian** correlations, **multilevel**
correlations, **polychoric** correlations, **biweight**, **percentage
bend** or **Sheperd’s Pi** correlations (types of robust correlation),
**distance** correlation (a type of non-linear correlation) and more,
also allowing for combinations between them (for instance, *Bayesian
partial multilevel correlation*).

# Citation

You can reference the package and its documentation as follows:

Makowski, D., Ben-Shachar, M. S., Patil, I., & Lüdecke, D. (2019).
Methods and Algorithms for Correlation Analysis in R. *Journal of Open
Source Software*, *5*(51), 2306. <https://doi.org/10.21105/joss.02306>

# Installation

[![CRAN](http://www.r-pkg.org/badges/version/correlation)](https://cran.r-project.org/package=correlation)
[![Build
Status](https://travis-ci.org/easystats/correlation.svg?branch=master)](https://travis-ci.org/easystats/correlation)
[![codecov](https://codecov.io/gh/easystats/correlation/branch/master/graph/badge.svg)](https://codecov.io/gh/easystats/correlation)

Run the following to install the stable release from CRAN:

``` r
install.packages("correlation")
```

Or this one to install the latest development version:

``` r
install.packages("remotes")
remotes::install_github("easystats/correlation")
```

# Documentation

[![Documentation](https://img.shields.io/badge/documentation-correlation-orange.svg?colorB=E91E63)](https://easystats.github.io/correlation/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-correlation-orange.svg?colorB=2196F3)](https://easystats.github.io/correlation/reference/index.html)

Click on the buttons above to access the package
[documentation](https://easystats.github.io/correlation/) and the
[easystats blog](https://easystats.github.io/blog/posts/), and check-out
these vignettes:

-   [Types of
    Correlation](https://easystats.github.io/correlation/articles/types.html)
-   [Multilevel
    Correlations](https://easystats.github.io/correlation/articles/multilevel.html)

# Features

The *correlation* package can compute many different types of
correlation, including:

✅ **Pearson’s correlation**<br> ✅ **Spearman’s rank correlation**<br>
✅ **Kendall’s rank correlation**<br> ✅ **Biweight midcorrelation**<br>
✅ **Distance correlation**<br> ✅ **Percentage bend correlation**<br>
✅ **Shepherd’s Pi correlation**<br> ✅ **Blomqvist’s coefficient**<br>
✅ **Hoeffding’s D**<br> ✅ **Gamma correlation**<br> ✅ **Gaussian rank
correlation**<br> ✅ **Point-Biserial and biserial correlation**<br> ✅
**Winsorized correlation**<br> ✅ **Polychoric correlation**<br> ✅
**Tetrachoric correlation**<br> ✅ **Multilevel correlation**<br>

An overview and description of these correlations types is [**available
here**](https://easystats.github.io/correlation/articles/types.html).
Moreover, many of these correlation types are available as **partial**
or within a **Bayesian** framework.

# Examples

The main function is
[`correlation()`](https://easystats.github.io/correlation/reference/correlation.html),
which builds on top of
[`cor_test()`](https://easystats.github.io/correlation/reference/cor_test.html)
and comes with a number of possible options.

## Correlation details and matrix

``` r
results <- correlation(iris)
results
## # Correlation Matrix (pearson-method)
## 
## Parameter1   |   Parameter2 |     r |         95% CI | t(148) |         p
## -------------------------------------------------------------------------
## Sepal.Length |  Sepal.Width | -0.12 | [-0.27,  0.04] |  -1.44 | 0.152    
## Sepal.Length | Petal.Length |  0.87 | [ 0.83,  0.91] |  21.65 | < .001***
## Sepal.Length |  Petal.Width |  0.82 | [ 0.76,  0.86] |  17.30 | < .001***
## Sepal.Width  | Petal.Length | -0.43 | [-0.55, -0.29] |  -5.77 | < .001***
## Sepal.Width  |  Petal.Width | -0.37 | [-0.50, -0.22] |  -4.79 | < .001***
## Petal.Length |  Petal.Width |  0.96 | [ 0.95,  0.97] |  43.39 | < .001***
## 
## p-value adjustment method: Holm (1979)
## Observations: 150
```

The output is not a square matrix, but a **(tidy) dataframe with all
correlations tests per row**. One can also obtain a **matrix** using:

``` r
summary(results)
## # Correlation Matrix (pearson-method)
## 
## Parameter    | Petal.Width | Petal.Length | Sepal.Width
## -------------------------------------------------------
## Sepal.Length |     0.82*** |      0.87*** |       -0.12
## Sepal.Width  |    -0.37*** |     -0.43*** |            
## Petal.Length |     0.96*** |              |            
## 
## p-value adjustment method: Holm (1979)
```

Note that one can also obtain the full, **square** and redundant matrix
using:

``` r
summary(results, redundant = TRUE)
## # Correlation Matrix (pearson-method)
## 
## Parameter    | Sepal.Length | Sepal.Width | Petal.Length | Petal.Width
## ----------------------------------------------------------------------
## Sepal.Length |      1.00*** |       -0.12 |      0.87*** |     0.82***
## Sepal.Width  |        -0.12 |     1.00*** |     -0.43*** |    -0.37***
## Petal.Length |      0.87*** |    -0.43*** |      1.00*** |     0.96***
## Petal.Width  |      0.82*** |    -0.37*** |      0.96*** |     1.00***
## 
## p-value adjustment method: Holm (1979)
```

``` r
library(see)

results %>%
  summary(redundant = TRUE) %>%
  plot()
```

![](man/figures/README-7-1.png)<!-- -->

## Correlation tests

The `cor_test()` function, for pairwise correlations, is also very
convenient for making quick scatter plots.

``` r
plot(cor_test(iris, "Sepal.Width", "Sepal.Length"))
```

![](man/figures/README-corr-1.png)<!-- -->

## Grouped dataframes

The `correlation()` function also supports **stratified correlations**,
all within the *tidyverse* workflow!

``` r
iris %>%
  select(Species, Sepal.Length, Sepal.Width, Petal.Width) %>%
  group_by(Species) %>%
  correlation()
## # Correlation Matrix (pearson-method)
## 
## Group      |   Parameter1 |  Parameter2 |    r |        95% CI | t(48) |         p
## ----------------------------------------------------------------------------------
## setosa     | Sepal.Length | Sepal.Width | 0.74 | [ 0.59, 0.85] |  7.68 | < .001***
## setosa     | Sepal.Length | Petal.Width | 0.28 | [ 0.00, 0.52] |  2.01 | 0.101    
## setosa     |  Sepal.Width | Petal.Width | 0.23 | [-0.05, 0.48] |  1.66 | 0.104    
## versicolor | Sepal.Length | Sepal.Width | 0.53 | [ 0.29, 0.70] |  4.28 | < .001***
## versicolor | Sepal.Length | Petal.Width | 0.55 | [ 0.32, 0.72] |  4.52 | < .001***
## versicolor |  Sepal.Width | Petal.Width | 0.66 | [ 0.47, 0.80] |  6.15 | < .001***
## virginica  | Sepal.Length | Sepal.Width | 0.46 | [ 0.20, 0.65] |  3.56 | 0.002**  
## virginica  | Sepal.Length | Petal.Width | 0.28 | [ 0.00, 0.52] |  2.03 | 0.048*   
## virginica  |  Sepal.Width | Petal.Width | 0.54 | [ 0.31, 0.71] |  4.42 | < .001***
## 
## p-value adjustment method: Holm (1979)
## Observations: 50
```

## Bayesian Correlations

It is very easy to switch to a **Bayesian framework**.

``` r
correlation(iris, bayesian = TRUE)
## # Correlation Matrix (pearson-method)
## 
## Parameter1   |   Parameter2 |   rho |         95% CI |      pd | % in ROPE |         Prior |        BF
## ------------------------------------------------------------------------------------------------------
## Sepal.Length |  Sepal.Width | -0.11 | [-0.26,  0.05] |  91.40% |    43.97% | Beta (3 +- 3) |     0.509
## Sepal.Length | Petal.Length |  0.86 | [ 0.82,  0.90] | 100%*** |        0% | Beta (3 +- 3) | > 1000***
## Sepal.Length |  Petal.Width |  0.81 | [ 0.75,  0.86] | 100%*** |        0% | Beta (3 +- 3) | > 1000***
## Sepal.Width  | Petal.Length | -0.41 | [-0.54, -0.28] | 100%*** |        0% | Beta (3 +- 3) | > 1000***
## Sepal.Width  |  Petal.Width | -0.35 | [-0.48, -0.20] | 100%*** |     0.12% | Beta (3 +- 3) | > 1000***
## Petal.Length |  Petal.Width |  0.96 | [ 0.95,  0.97] | 100%*** |        0% | Beta (3 +- 3) | > 1000***
## 
## Observations: 150
```

## Tetrachoric, Polychoric, Biserial, Biweight…

The `correlation` package also supports different types of methods,
which can deal with correlations **between factors**!

``` r
correlation(iris, include_factors = TRUE, method = "auto")
## # Correlation Matrix (auto-method)
## 
## Parameter1         |         Parameter2 |     r |         95% CI | t(148) |         p
## -------------------------------------------------------------------------------------
## Sepal.Length       |        Sepal.Width | -0.12 | [-0.27,  0.04] |  -1.44 | 0.452    
## Sepal.Length       |       Petal.Length |  0.87 | [ 0.83,  0.91] |  21.65 | < .001***
## Sepal.Length       |        Petal.Width |  0.82 | [ 0.76,  0.86] |  17.30 | < .001***
## Sepal.Length       |     Species.setosa | -0.72 | [-0.79, -0.63] | -12.53 | < .001***
## Sepal.Length       | Species.versicolor |  0.08 | [-0.08,  0.24] |   0.97 | 0.452    
## Sepal.Length       |  Species.virginica |  0.64 | [ 0.53,  0.72] |  10.08 | < .001***
## Sepal.Width        |       Petal.Length | -0.43 | [-0.55, -0.29] |  -5.77 | < .001***
## Sepal.Width        |        Petal.Width | -0.37 | [-0.50, -0.22] |  -4.79 | < .001***
## Sepal.Width        |     Species.setosa |  0.60 | [ 0.49,  0.70] |   9.20 | < .001***
## Sepal.Width        | Species.versicolor | -0.47 | [-0.58, -0.33] |  -6.44 | < .001***
## Sepal.Width        |  Species.virginica | -0.14 | [-0.29,  0.03] |  -1.67 | 0.392    
## Petal.Length       |        Petal.Width |  0.96 | [ 0.95,  0.97] |  43.39 | < .001***
## Petal.Length       |     Species.setosa | -0.92 | [-0.94, -0.89] | -29.13 | < .001***
## Petal.Length       | Species.versicolor |  0.20 | [ 0.04,  0.35] |   2.51 | 0.066    
## Petal.Length       |  Species.virginica |  0.72 | [ 0.63,  0.79] |  12.66 | < .001***
## Petal.Width        |     Species.setosa | -0.89 | [-0.92, -0.85] | -23.41 | < .001***
## Petal.Width        | Species.versicolor |  0.12 | [-0.04,  0.27] |   1.44 | 0.452    
## Petal.Width        |  Species.virginica |  0.77 | [ 0.69,  0.83] |  14.66 | < .001***
## Species.setosa     | Species.versicolor | -0.88 | [-0.91, -0.84] | -22.43 | < .001***
## Species.setosa     |  Species.virginica | -0.88 | [-0.91, -0.84] | -22.43 | < .001***
## Species.versicolor |  Species.virginica | -0.88 | [-0.91, -0.84] | -22.43 | < .001***
## 
## p-value adjustment method: Holm (1979)
## Observations: 150
```

## Partial Correlations

It also supports **partial correlations** (as well as Bayesian partial
correlations).

``` r
iris %>%
  correlation(partial = TRUE) %>%
  summary()
## # Correlation Matrix (pearson-method)
## 
## Parameter    | Petal.Width | Petal.Length | Sepal.Width
## -------------------------------------------------------
## Sepal.Length |    -0.34*** |      0.72*** |     0.63***
## Sepal.Width  |     0.35*** |     -0.62*** |            
## Petal.Length |     0.87*** |              |            
## 
## p-value adjustment method: Holm (1979)
```

## Gaussian Graphical Models (GGMs)

Such partial correlations can also be represented as **Gaussian
Graphical Models** (GGM), an increasingly popular tool in psychology. A
GGM traditionally include a set of variables depicted as circles
(“nodes”), and a set of lines that visualize relationships between them,
which thickness represents the strength of association (see [Bhushan et
al.,
2019](https://www.frontiersin.org/articles/10.3389/fpsyg.2019.01050/full)).

``` r
library(see) # for plotting
library(ggraph) # needs to be loaded

mtcars %>%
  correlation(partial = TRUE) %>%
  plot()
```

![](man/figures/README-12-1.png)<!-- -->

## Multilevel Correlations

It also provide some cutting-edge methods, such as Multilevel (partial)
correlations. These are are partial correlations based on linear
mixed-effects models that include the factors as **random effects**.
They can be see as correlations *adjusted* for some group
(*hierarchical*) variability.

``` r
iris %>%
  correlation(partial = TRUE, multilevel = TRUE) %>%
  summary()
## # Correlation Matrix (pearson-method)
## 
## Parameter    | Petal.Width | Petal.Length | Sepal.Width
## -------------------------------------------------------
## Sepal.Length |      -0.17* |      0.71*** |     0.43***
## Sepal.Width  |     0.39*** |       -0.18* |            
## Petal.Length |     0.38*** |              |            
## 
## p-value adjustment method: Holm (1979)
```

However, if the `partial` argument is set to `FALSE`, it will try to
convert the partial coefficient into regular ones.These can be
**converted back** to full correlations:

``` r
iris %>%
  correlation(partial = FALSE, multilevel = TRUE) %>%
  summary()
## Parameter    | Petal.Width | Petal.Length | Sepal.Width
## -------------------------------------------------------
## Sepal.Length |     0.36*** |      0.76*** |     0.53***
## Sepal.Width  |     0.47*** |      0.38*** |            
## Petal.Length |     0.48*** |              |
```

# Contributing and Support

In case you want to file an issue or contribute in another way to the
package, please follow [this
guide](https://easystats.github.io/correlation/CONTRIBUTING.html). For
questions about the functionality, you may either contact us via email
or also file an issue.

# Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://easystats.github.io/correlation/CODE_OF_CONDUCT.html).
By participating in this project you agree to abide by its terms.
# correlation 0.8.0

## Breaking Changes

- `robust` argument, which was deprecated in favour of `ranktransform` in
  `0.6.1` release, has now been removed.

# correlation 0.7.1

## Bug Fixes

- Bug fix in `plot()` methods

# correlation 0.7.0

## Breaking Changes

- Removes `winsorize()` function, which now lives in `datawizard` package.

## New Features

- New `cor_smooth()` function for smoothing non-positive definite matrices.

## Bug Fixes

- When `data2` was specified `correlation()` was over-correcting for all of the
  combinations of variables in the full x and y tables, rather than in just the
  ones specified (#195).

## Minor Changes

- `correlation()` gains a new argument `rename` to rename variables.

- `simualte_simpson()` function is now re-exported from `bayestestR` package.

- `plot()` for `"easycor_test"` objects now produces an annotated scatter plot.

# correlation 0.6.1

## Breaking Changes

- `simualte_simpson()`: The groups are now named after the pattern `"G_"` (can
  be altered with the `group_prefix` argument).

- `robust` argument deprecated in favour of `ranktransform`.

## New Features

- `correlation` gains two new arguments: `select` and `select2` to select
  specific variables from dataframes to compare (#146).

- `as.matrix` method works for grouped correlations (#148).

- New `as.list` method returns a list of various matrices related to correlation
  analysis (correlation, number of observations, *p*-values, etc.).

## Bug Fixes

- The `0.6.0` release introduced a bug in Winsorized Pearson correlation where
  the missing values were removed from the entire data, instead for each pair
  (#151). This is now fixed.

# correlation 0.6.0

## New Features

- Added `verbose` arguments to some functions, to toggle warnings on/off.

- `cor_test()` (and hence, `correlation()`) now default the `winsorize` argument
  to `.1` when it's set to `TRUE`.

- The `Method` column in output dataframe is now more explicit about the
  correlation method used.

## Bug Fixes

- Winsorization doesn't fail when `NA`s are present (#130).

## Minor Changes

- Fixed CRAN check issues due to changes in dependent packages.

# correlation 0.5.0

## Changes

- Added `winsorize()` function.

- Added `winsorize` argument for Winsorized correlations.

- Added `method = "somers"` to `correlation()`, to compute Somers's Dxy rank
  correlation for binary outcomes.

- New function `display()`, to print output into different formats. Currently,
  only markdown is supported. `print_md()` is an alias for `display(format =
  "markdown")`.

## Bug fixes

- Fix bug in `cor_to_p()` that gave slightly different test statistics.

# correlation 0.4.0

## Changes

- Don't error if less than 3 valid observations
  ([#100](https://github.com/easystats/correlation/issues/100)).

- Add "gaussian" rank method.

- Add "gamma" method.

- Add "hoeffding" method.

- Add "blomqvist" method.

## Bug fixes

- Added `Method` column to Bayesian correlations.

- Fix bug when `robust=TRUE`
  ([#87](https://github.com/easystats/effectsize/issues/87)).

# correlation 0.3.0

## Changes

## Bug fixes

# correlation 0.2.1

## Changes

- Added confidence intervals CI support for Spearman and Kendall (#80)

- Improved documentation (#45, #63)

## Bug fixes

- Removed CI threshold column from `distance_mahalanobis()`

- Fixed bug (#76)

# correlation 0.2.0

## Changes

- Some changes were made.

## Bug fixes

- Some bugs were fixed.

# correlation 0.1.0

## Changes

- Initial CRAN release.

- Add `plot()`-method for `summary()`.

## Bug fixes

- Fixed issue in `correlation()` for some edge cases when `include_factors =
  TRUE`.

- Fixed issue in `correlation()` for correlation coefficients with less than
  four complete pairs of observations (in such cases, `cor_test()` now returns
  `NA` for the confidence intervals).

## Test environments
* local R installation, R 4.1.1
* ubuntu 16.04 (on github-actions), R 4.1.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 note

## revdepcheck results

We checked 4 reverse dependencies, comparing R CMD check results across CRAN and
dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
---
title: 'Methods and Algorithms for Correlation Analysis in R'
tags:
- R
- Correlation
- Easystats
authors:
- affiliation: 1
  name: Dominique Makowski
  orcid: 0000-0001-5375-9967
- affiliation: 2
  name: Mattan S. Ben-Shachar
  orcid: 0000-0002-4287-4801
- affiliation: 3
  name: Indrajeet Patil
  orcid: 0000-0003-1995-6531
- affiliation: 4
  name: Daniel Lüdecke
  orcid: 0000-0002-8895-3206
affiliations:
- index: 1
  name: Nanyang Technological University, Singapore
- index: 2
  name: Ben-Gurion University of the Negev, Israel
- index: 3
  name: Max Planck Institute for Human Development, Germany
- index: 4
  name: University Medical Center Hamburg-Eppendorf, Germany
date: "22 March 2020"
bibliography: paper.bib
csl: apa.csl
output: pdf_document
---

# Introduction

Correlations tests are arguably one of the most commonly used statistical procedures, and are used as a basis in many applications such as exploratory data analysis, structural modelling, data engineering etc. In this context, we present **correlation**, a toolbox for the R language [@Rteam] and part of the [**easystats**](https://github.com/easystats/easystats) collection, focused on correlation analysis. Its goal is to be lightweight, easy to use, and allows for the computation of many different kinds of correlations, such as:

- **Pearson's correlation**: This is the most common correlation method. It corresponds to the covariance of the two variables normalized (i.e., divided) by the product of their standard deviations.

$$r_{xy} = \frac{cov(x,y)}{SD_x \times SD_y}$$

- **Spearman's rank correlation**: A non-parametric measure of correlation, the Spearman correlation between two variables is equal to the Pearson correlation between the rank scores of those two variables; while Pearson's correlation assesses linear relationships, Spearman's correlation assesses monotonic relationships (whether linear or not). Confidence Intervals (CI) for Spearman's correlations are computed using the @fieller1957tests correction [see @bishara2017confidence].

$$r_{s_{xy}} = \frac{cov(rank_x, rank_y)}{SD(rank_x) \times SD(rank_y)}$$

- **Kendall's rank correlation**: In the normal case, the Kendall correlation is preferred to the Spearman correlation because of a smaller gross error sensitivity (GES) and a smaller asymptotic variance (AV), making it more robust and more efficient. However, the interpretation of Kendall's tau is less direct compared to that of the Spearman's rho, in the sense that it quantifies the difference between the % of concordant and discordant pairs among all possible pairwise events. Confidence Intervals (CI) for Kendall's correlations are computed using the @fieller1957tests correction [see @bishara2017confidence]. For each pair of observations (i ,j) of two variables (x, y), it is defined as follows:

$$\tau_{xy} = \frac{2}{n(n-1)}\sum_{i<j}^{}sign(x_i - x_j) \times sign(y_i - y_j)$$

- **Biweight midcorrelation**: A measure of similarity that is median-based, instead of the traditional mean-based, thus being less sensitive to outliers. It can be used as a robust alternative to other similarity metrics, such as Pearson correlation [@langfelder2012fast].

- **Distance correlation**: Distance correlation measures both linear and non-linear association between two random variables or random vectors. This is in contrast to Pearson's correlation, which can only detect linear association between two random variables.

- **Percentage bend correlation**: Introduced by Wilcox (1994), it is based on a down-weight of a specified percentage of marginal observations deviating from the median (by default, 20 percent).

- **Shepherd's Pi correlation**: Equivalent to a Spearman's rank correlation after outliers removal (by means of bootstrapped Mahalanobis distance).

- **Point-Biserial and biserial correlation**: Correlation coefficient used when one variable is continuous and the other is dichotomous (binary). Point-Biserial is equivalent to a Pearson's correlation, while Biserial should be used when the binary variable is assumed to have an underlying continuity. For example, anxiety level can be measured on a continuous scale, but can be classified dichotomously as high/low.

- **Polychoric correlation**: Correlation between two theorised normally distributed continuous latent variables, from two observed ordinal variables.

- **Tetrachoric correlation**: Special case of the polychoric correlation applicable when both observed variables are dichotomous.

- **Partial correlation**: Correlation between two variables after adjusting for the (linear) the effect of one or more variables. The correlation test is here run after having partialized the dataset, independently from it. In other words, it considers partialization as an independent step generating a different dataset, rather than belonging to the same model. This is why some discrepancies are to be expected for the *t*- and the *p*-values (but not the correlation coefficient) compared to other implementations such as **ppcor**. Let $e_{x.z}$ be the residuals from the linear prediction of $x$ by $z$ (note that this can be expanded to a multivariate $z$):

$$r_{xy.z} = r_{e_{x.z},e_{y.z}}$$

- **Multilevel correlation**: Multilevel correlations are a special case of partial correlations where the variable to be adjusted for is a factor and is included as a random effect in a mixed model.

These methods allow for different ways of quantifying the link between two variables (see **Figure 1**).

![Illustration of the different correlation estimates (a measure of association, represent by the height of the bars) obtained via different methods for the same data (the scatter plot).](figure1.png)

# Design


It relies on one main function, `correlation()`, which outputs a dataframe containing each pairwise correlation per row. This long format is convenient for further data analysis, but not as much to get a summary, which is usually obtained via a correlation matrix. To address this, we added standard methods, such as `summary()` and `as.matrix()`, to automatically transform the long output to a matrix. Moreover, **correlation** also includes plotting capabilities via the [**see** package](https://easystats.github.io/see/) [@ludecke2019see].

An overview of the features is available on the GitHub page (https://github.com/easystats/correlation). The typical core workflow is as follows:

\footnotesize

``` r
results <- correlation(iris)
results
# Parameter1   |   Parameter2 |     r |         95% CI |     t |  df |      p |  Method | n_Obs
# ---------------------------------------------------------------------------------------------
# Sepal.Length |  Sepal.Width | -0.12 | [-0.27,  0.04] | -1.44 | 148 | 0.152  | Pearson |   150
# Sepal.Length | Petal.Length |  0.87 | [ 0.83,  0.91] | 21.65 | 148 | < .001 | Pearson |   150
# Sepal.Length |  Petal.Width |  0.82 | [ 0.76,  0.86] | 17.30 | 148 | < .001 | Pearson |   150
# Sepal.Width  | Petal.Length | -0.43 | [-0.55, -0.29] | -5.77 | 148 | < .001 | Pearson |   150
# Sepal.Width  |  Petal.Width | -0.37 | [-0.50, -0.22] | -4.79 | 148 | < .001 | Pearson |   150
# Petal.Length |  Petal.Width |  0.96 | [ 0.95,  0.97] | 43.39 | 148 | < .001 | Pearson |   150
```

\normalsize

The output is not a square matrix, but a (tidy) dataframe with all correlations tests per row. One can also obtain a matrix using:

\small

``` r
summary(results)
# Parameter    | Petal.Width | Petal.Length | Sepal.Width
# -------------------------------------------------------
# Sepal.Length |     0.82*** |      0.87*** |       -0.12
# Sepal.Width  |    -0.37*** |     -0.43*** |            
# Petal.Length |     0.96*** |              |
```

\normalsize



# Availability

The **correlation** package can be downloaded and installed from CRAN [1](https://CRAN.R-project.org/package=correlation). It is licensed under the GNU General Public License (v3.0), with all its source code stored at GitHub [2](https://github.com/easystats/correlation), and with a corresponding issue tracker [2](https://github.com/easystats/correlation/issues) for bug reporting and feature enhancements. In the spirit of honest and open science, we encourage requests/tips for fixes, feature updates, as well as general questions and concerns via direct interaction with contributors and developers.

# Acknowledgments

**correlation** is part of the [*easystats*](https://github.com/easystats/easystats) ecosystem [relying on **insight**; @ludecke2019insight and **bayestestR**; @makowski2019bayestestr], a collaborative project created to facilitate the usage of R. Thus, we would like to thank the [council of masters](https://github.com/orgs/easystats/people) of easystats, all other padawan contributors, as well as the users.

# References
# as.data.frame for correlation output

    Code
      as.data.frame(correlation(ggplot2::msleep))
    Output
          Parameter1  Parameter2          r   CI      CI_low     CI_high            t
      1  sleep_total   sleep_rem  0.7517550 0.95  0.61667557  0.84383201     8.756396
      2  sleep_total sleep_cycle -0.4737127 0.95 -0.70581894 -0.14975542    -2.946170
      3  sleep_total       awake -0.9999986 0.95 -0.99999908 -0.99999779 -5328.711772
      4  sleep_total     brainwt -0.3604874 0.95 -0.56942242 -0.10780364    -2.839979
      5  sleep_total      bodywt -0.3120106 0.95 -0.49442632 -0.10327118    -2.955645
      6    sleep_rem sleep_cycle -0.3381235 0.95 -0.61438094  0.01198335    -1.967883
      7    sleep_rem       awake -0.7517713 0.95 -0.84384279 -0.61669876    -8.756832
      8    sleep_rem     brainwt -0.2213348 0.95 -0.47556189  0.06701441    -1.539344
      9    sleep_rem      bodywt -0.3276507 0.95 -0.53530394 -0.08264933    -2.663776
      10 sleep_cycle       awake  0.4737127 0.95  0.14975542  0.70581894     2.946170
      11 sleep_cycle     brainwt  0.8516203 0.95  0.70882870  0.92736294     8.597296
      12 sleep_cycle      bodywt  0.4178029 0.95  0.08089399  0.66902912     2.518773
      13       awake     brainwt  0.3604874 0.95  0.10780364  0.56942242     2.839979
      14       awake      bodywt  0.3119801 0.95  0.10323781  0.49440083     2.955326
      15     brainwt      bodywt  0.9337822 0.95  0.88916423  0.96081138    19.175704
         df_error             p              Method n_Obs
      1        59  3.783810e-11 Pearson correlation    61
      2        30  4.934837e-02 Pearson correlation    32
      3        81 3.627785e-225 Pearson correlation    83
      4        54  4.934837e-02 Pearson correlation    56
      5        81  4.085332e-02 Pearson correlation    83
      6        30  1.167709e-01 Pearson correlation    32
      7        59  3.783810e-11 Pearson correlation    61
      8        46  1.305716e-01 Pearson correlation    48
      9        59  4.934837e-02 Pearson correlation    61
      10       30  4.934837e-02 Pearson correlation    32
      11       28  2.662362e-08 Pearson correlation    30
      12       30  5.202211e-02 Pearson correlation    32
      13       54  4.934837e-02 Pearson correlation    56
      14       81  4.085332e-02 Pearson correlation    83
      15       54  1.281756e-24 Pearson correlation    56

# summary.correlation - target column

    Code
      summary(correlation(ggplot2::msleep), target = "t")
    Output
      # Correlation Matrix (pearson-method)
      
      Parameter   |   bodywt | brainwt |       awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------
      sleep_total |   -2.96* |  -2.84* | -5328.71*** |      -2.95* |   8.76***
      sleep_rem   |   -2.66* |   -1.54 |    -8.76*** |       -1.97 |          
      sleep_cycle |     2.52 | 8.60*** |       2.95* |             |          
      awake       |    2.96* |   2.84* |             |             |          
      brainwt     | 19.18*** |         |             |             |          
      
      p-value adjustment method: Holm (1979)

---

    Code
      summary(correlation(ggplot2::msleep), target = "df_error")
    Output
      # Correlation Matrix (pearson-method)
      
      Parameter   |   bodywt |  brainwt |    awake | sleep_cycle | sleep_rem
      ----------------------------------------------------------------------
      sleep_total |   81.00* |   54.00* | 81.00*** |      30.00* |  59.00***
      sleep_rem   |   59.00* |    46.00 | 59.00*** |       30.00 |          
      sleep_cycle |    30.00 | 28.00*** |   30.00* |             |          
      awake       |   81.00* |   54.00* |          |             |          
      brainwt     | 54.00*** |          |          |             |          
      
      p-value adjustment method: Holm (1979)

---

    Code
      summary(correlation(ggplot2::msleep), target = "p")
    Output
      # Correlation Matrix (pearson-method)
      
      Parameter   |   bodywt |  brainwt |     awake | sleep_cycle | sleep_rem
      -----------------------------------------------------------------------
      sleep_total |     0.04 |     0.05 | 3.63e-225 |        0.05 |  3.78e-11
      sleep_rem   |     0.05 |     0.13 |  3.78e-11 |        0.12 |          
      sleep_cycle |     0.05 | 2.66e-08 |      0.05 |             |          
      awake       |     0.04 |     0.05 |           |             |          
      brainwt     | 1.28e-24 |          |           |             |          
      
      p-value adjustment method: Holm (1979)

# as.list

    Code
      as.list(correlation(mtcars))
    Output
       r 
      ---
      Parameter |  carb |  gear |    am |    vs |  qsec |    wt |  drat |    hp |  disp |   cyl
      -----------------------------------------------------------------------------------------
      mpg       | -0.55 |  0.48 |  0.60 |  0.66 |  0.42 | -0.87 |  0.68 | -0.78 | -0.85 | -0.85
      cyl       |  0.53 | -0.49 | -0.52 | -0.81 | -0.59 |  0.78 | -0.70 |  0.83 |  0.90 |      
      disp      |  0.39 | -0.56 | -0.59 | -0.71 | -0.43 |  0.89 | -0.71 |  0.79 |       |      
      hp        |  0.75 | -0.13 | -0.24 | -0.72 | -0.71 |  0.66 | -0.45 |       |       |      
      drat      | -0.09 |  0.70 |  0.71 |  0.44 |  0.09 | -0.71 |       |       |       |      
      wt        |  0.43 | -0.58 | -0.69 | -0.55 | -0.17 |       |       |       |       |      
      qsec      | -0.66 | -0.21 | -0.23 |  0.74 |       |       |       |       |       |      
      vs        | -0.57 |  0.21 |  0.17 |       |       |       |       |       |       |      
      am        |  0.06 |  0.79 |       |       |       |       |       |       |       |      
      gear      |  0.27 |       |       |       |       |       |       |       |       |      
      
       n_Obs 
      -------
      Parameter |  carb |  gear |    am |    vs |  qsec |    wt |  drat |    hp |  disp |   cyl
      -----------------------------------------------------------------------------------------
      mpg       | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00
      cyl       | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 |      
      disp      | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 |       |      
      hp        | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 |       |       |      
      drat      | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 |       |       |       |      
      wt        | 32.00 | 32.00 | 32.00 | 32.00 | 32.00 |       |       |       |       |      
      qsec      | 32.00 | 32.00 | 32.00 | 32.00 |       |       |       |       |       |      
      vs        | 32.00 | 32.00 | 32.00 |       |       |       |       |       |       |      
      am        | 32.00 | 32.00 |       |       |       |       |       |       |       |      
      gear      | 32.00 |       |       |       |       |       |       |       |       |      
      
       p 
      ---
      Parameter |     carb |     gear |       am |       vs |     qsec |       wt |     drat |       hp |     disp |      cyl
      -----------------------------------------------------------------------------------------------------------------------
      mpg       |     0.02 |     0.10 | 8.27e-03 | 1.09e-03 |     0.22 | 6.86e-09 | 5.86e-04 | 8.05e-06 | 4.78e-08 | 3.18e-08
      cyl       |     0.04 |     0.08 |     0.04 | 9.03e-07 |     0.01 | 5.60e-06 | 2.97e-04 | 1.74e-07 | 9.92e-11 |         
      disp      |     0.30 |     0.02 |     0.01 | 2.04e-04 |     0.20 | 6.60e-10 | 2.04e-04 | 3.36e-06 |          |         
      hp        | 3.44e-05 |     1.00 |     1.00 | 1.24e-04 | 2.13e-04 | 1.29e-03 |     0.17 |          |          |         
      drat      |     1.00 | 2.97e-04 | 1.94e-04 |     0.19 |     1.00 | 1.94e-04 |          |          |          |         
      wt        |     0.20 |     0.01 | 3.83e-04 |     0.02 |     1.00 |          |          |          |          |         
      qsec      | 1.36e-03 |     1.00 |     1.00 | 4.43e-05 |          |          |          |          |          |         
      vs        |     0.02 |     1.00 |     1.00 |          |          |          |          |          |          |         
      am        |     1.00 | 2.80e-06 |          |          |          |          |          |          |          |         
      gear      |     1.00 |          |          |          |          |          |          |          |          |         
      

---

    Code
      suppressWarnings(as.list(msleep %>% group_by(vore) %>% correlation(method = "spearman")))
    Output
      =======
       carni 
      =======
      
       rho 
      -----
      Group |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------
      carni | sleep_total |  -0.48 |   -0.59 | -1.00 |        0.31 |      0.95
      carni |   sleep_rem |  -0.72 |   -0.26 | -0.95 |        0.46 |          
      carni | sleep_cycle |  -0.56 |   -0.80 | -0.31 |             |          
      carni |       awake |   0.48 |    0.59 |       |             |          
      carni |     brainwt |   0.82 |         |       |             |          
      
       n_Obs 
      -------
      Group |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------
      carni | sleep_total |  19.00 |    9.00 | 19.00 |        5.00 |     10.00
      carni |   sleep_rem |  10.00 |    6.00 | 10.00 |        5.00 |          
      carni | sleep_cycle |   5.00 |    4.00 |  5.00 |             |          
      carni |       awake |  19.00 |    9.00 |       |             |          
      carni |     brainwt |   9.00 |         |       |             |          
      
       p 
      ---
      Group |   Parameter | bodywt | brainwt |    awake | sleep_cycle | sleep_rem
      ---------------------------------------------------------------------------
      carni | sleep_total |   0.37 |    0.73 |     0.00 |        1.00 |  3.19e-04
      carni |   sleep_rem |   0.20 |    1.00 | 3.19e-04 |        1.00 |          
      carni | sleep_cycle |   1.00 |    1.00 |     1.00 |             |          
      carni |       awake |   0.37 |    0.73 |          |             |          
      carni |     brainwt |   0.09 |         |          |             |          
      
      
      =======
       herbi 
      =======
      
       rho 
      -----
      Group |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------
      herbi | sleep_total |  -0.77 |   -0.86 | -1.00 |       -0.44 |      0.92
      herbi |   sleep_rem |  -0.72 |   -0.75 | -0.92 |       -0.48 |          
      herbi | sleep_cycle |   0.74 |    0.74 |  0.44 |             |          
      herbi |       awake |   0.77 |    0.86 |       |             |          
      herbi |     brainwt |   0.98 |         |       |             |          
      
       n_Obs 
      -------
      Group |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------
      herbi | sleep_total |  32.00 |   20.00 | 32.00 |       12.00 |     24.00
      herbi |   sleep_rem |  24.00 |   16.00 | 24.00 |       12.00 |          
      herbi | sleep_cycle |  12.00 |   11.00 | 12.00 |             |          
      herbi |       awake |  32.00 |   20.00 |       |             |          
      herbi |     brainwt |  20.00 |         |       |             |          
      
       p 
      ---
      Group |   Parameter |   bodywt |  brainwt |    awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------------
      herbi | sleep_total | 3.42e-06 | 8.71e-06 |     0.00 |        0.35 |  5.00e-09
      herbi |   sleep_rem | 4.65e-04 | 4.43e-03 | 5.00e-09 |        0.35 |          
      herbi | sleep_cycle |     0.03 |     0.04 |     0.35 |             |          
      herbi |       awake | 3.42e-06 | 8.71e-06 |          |             |          
      herbi |     brainwt | 5.17e-13 |          |          |             |          
      
      
      =========
       insecti 
      =========
      
       rho 
      -----
      Group   |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      --------------------------------------------------------------------------
      insecti | sleep_total |  -0.60 |   -0.60 | -1.00 |        0.50 |     -0.40
      insecti |   sleep_rem |   0.80 |    0.80 |  0.40 |       -1.00 |          
      insecti | sleep_cycle |  -0.50 |   -0.50 | -0.50 |             |          
      insecti |       awake |   0.60 |    0.60 |       |             |          
      insecti |     brainwt |   1.00 |         |       |             |          
      
       n_Obs 
      -------
      Group   |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      --------------------------------------------------------------------------
      insecti | sleep_total |   5.00 |    5.00 |  5.00 |        3.00 |      4.00
      insecti |   sleep_rem |   4.00 |    4.00 |  4.00 |        3.00 |          
      insecti | sleep_cycle |   3.00 |    3.00 |  3.00 |             |          
      insecti |       awake |   5.00 |    5.00 |       |             |          
      insecti |     brainwt |   5.00 |         |       |             |          
      
       p 
      ---
      Group   |   Parameter |   bodywt | brainwt |    awake | sleep_cycle | sleep_rem
      -------------------------------------------------------------------------------
      insecti | sleep_total |     1.00 |    1.00 | 1.46e-22 |        1.00 |      1.00
      insecti |   sleep_rem |     1.00 |    1.00 |     1.00 |        0.00 |          
      insecti | sleep_cycle |     1.00 |    1.00 |     1.00 |             |          
      insecti |       awake |     1.00 |    1.00 |          |             |          
      insecti |     brainwt | 5.56e-23 |         |          |             |          
      
      
      ======
       omni 
      ======
      
       rho 
      -----
      Group |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------
      omni  | sleep_total |  -0.10 |   -0.28 | -1.00 |       -0.24 |      0.14
      omni  |   sleep_rem |  -0.20 |   -0.39 | -0.14 |       -0.46 |          
      omni  | sleep_cycle |   0.80 |    0.92 |  0.24 |             |          
      omni  |       awake |   0.10 |    0.28 |       |             |          
      omni  |     brainwt |   0.91 |         |       |             |          
      
       n_Obs 
      -------
      Group |   Parameter | bodywt | brainwt | awake | sleep_cycle | sleep_rem
      ------------------------------------------------------------------------
      omni  | sleep_total |  20.00 |   17.00 | 20.00 |       11.00 |     18.00
      omni  |   sleep_rem |  18.00 |   17.00 | 18.00 |       11.00 |          
      omni  | sleep_cycle |  11.00 |   11.00 | 11.00 |             |          
      omni  |       awake |  20.00 |   17.00 |       |             |          
      omni  |     brainwt |  17.00 |         |       |             |          
      
       p 
      ---
      Group |   Parameter |   bodywt |  brainwt | awake | sleep_cycle | sleep_rem
      ---------------------------------------------------------------------------
      omni  | sleep_total |     1.00 |     1.00 |  0.00 |        1.00 |      1.00
      omni  |   sleep_rem |     1.00 |     1.00 |  1.00 |        1.00 |          
      omni  | sleep_cycle |     0.04 | 7.73e-04 |  1.00 |             |          
      omni  |       awake |     1.00 |     1.00 |       |             |          
      omni  |     brainwt | 7.64e-06 |          |       |             |          
      
      

---

    Code
      suppressWarnings(as.list(mtcars %>% group_by(am) %>% correlation(select = c(
        "cyl", "wt"), select2 = c("hp"), method = "percentage")))
    Output
      ===
       0 
      ===
      
       r 
      ---
      Group | Parameter |   hp
      ------------------------
      0     |       cyl | 0.87
      0     |        wt | 0.83
      
       n_Obs 
      -------
      Group | Parameter |    hp
      -------------------------
      0     |       cyl | 19.00
      0     |        wt | 19.00
      
       p 
      ---
      Group | Parameter |       hp
      ----------------------------
      0     |       cyl | 2.11e-06
      0     |        wt | 1.11e-05
      
      
      ===
       1 
      ===
      
       r 
      ---
      Group | Parameter |   hp
      ------------------------
      1     |       cyl | 0.83
      1     |        wt | 0.80
      
       n_Obs 
      -------
      Group | Parameter |    hp
      -------------------------
      1     |       cyl | 13.00
      1     |        wt | 13.00
      
       p 
      ---
      Group | Parameter |       hp
      ----------------------------
      1     |       cyl | 9.58e-04
      1     |        wt | 1.04e-03
      
      

# display and print method works - markdown

    Code
      print(summary(correlation(iris)))
    Output
      # Correlation Matrix (pearson-method)
      
      Parameter    | Petal.Width | Petal.Length | Sepal.Width
      -------------------------------------------------------
      Sepal.Length |     0.82*** |      0.87*** |       -0.12
      Sepal.Width  |    -0.37*** |     -0.43*** |            
      Petal.Length |     0.96*** |              |            
      
      p-value adjustment method: Holm (1979)

# display and print method works - html

    Code
      print(summary(correlation(iris)), format = "html")
    Output
      Correlation Matrix (pearson-method)
      
      Parameter    | Petal.Width | Petal.Length | Sepal.Width
      -------------------------------------------------------
      Sepal.Length |     0.82*** |      0.87*** |       -0.12
      Sepal.Width  |    -0.37*** |     -0.43*** |            
      Petal.Length |     0.96*** |              |            
      p-value adjustment method: Holm (1979)

# as.matrix works

    Code
      list(mat1, mat2)
    Output
      [[1]]
                 am         wt         hp
      am  1.0000000 -0.6924953 -0.2432043
      wt -0.6924953  1.0000000  0.6587479
      hp -0.2432043  0.6587479  1.0000000
      
      [[2]]
                    wt        hp
      0 - wt 1.0000000 0.6797596
      0 - hp 0.6797596 1.0000000
      1 - wt 1.0000000 0.8145279
      1 - hp 0.8145279 1.0000000
      

# renaming columns

    Code
      correlation(anscombe, select = c("x1", "x2"), rename = c("var1"))
    Warning <simpleWarning>
      Mismatch between number of variables and names.
    Output
      # Correlation Matrix (pearson-method)
      
      Parameter1 | Parameter2 |    r |       95% CI | t(9) |         p
      ----------------------------------------------------------------
      x1         |         x2 | 1.00 | [1.00, 1.00] |  Inf | < .001***
      
      p-value adjustment method: Holm (1979)
      Observations: 11

---

    Code
      correlation(anscombe, select = c("x1", "x2"), rename = c("var1", "var2"))
    Output
      # Correlation Matrix (pearson-method)
      
      Parameter1 | Parameter2 |    r |       95% CI | t(9) |         p
      ----------------------------------------------------------------
      var1       |       var2 | 1.00 | [1.00, 1.00] |  Inf | < .001***
      
      p-value adjustment method: Holm (1979)
      Observations: 11

---

    Code
      correlation(anscombe, select = c("x1", "x2"), select2 = c("y1", "y2"), rename = c(
        "var1", "var2"))
    Output
      # Correlation Matrix (pearson-method)
      
      Parameter1 | Parameter2 |    r |       95% CI | t(9) |       p
      --------------------------------------------------------------
      var1       |         y1 | 0.82 | [0.42, 0.95] | 4.24 | 0.009**
      var1       |         y2 | 0.82 | [0.42, 0.95] | 4.24 | 0.009**
      var2       |         y1 | 0.82 | [0.42, 0.95] | 4.24 | 0.009**
      var2       |         y2 | 0.82 | [0.42, 0.95] | 4.24 | 0.009**
      
      p-value adjustment method: Holm (1979)
      Observations: 11

# display and print method works - markdown

    Code
      display(correlation(iris))
    Output
      
      
      Table: Correlation Matrix (pearson-method)
      
      |Parameter1   |   Parameter2 |     r |         95% CI | t(148) |         p |
      |:------------|:------------:|:-----:|:--------------:|:------:|:---------:|
      |Sepal.Length |  Sepal.Width | -0.12 |  (-0.27, 0.04) |  -1.44 | 0.152     |
      |Sepal.Length | Petal.Length |  0.87 |   (0.83, 0.91) |  21.65 | < .001*** |
      |Sepal.Length |  Petal.Width |  0.82 |   (0.76, 0.86) |  17.30 | < .001*** |
      |Sepal.Width  | Petal.Length | -0.43 | (-0.55, -0.29) |  -5.77 | < .001*** |
      |Sepal.Width  |  Petal.Width | -0.37 | (-0.50, -0.22) |  -4.79 | < .001*** |
      |Petal.Length |  Petal.Width |  0.96 |   (0.95, 0.97) |  43.39 | < .001*** |
      p-value adjustment method: Holm (1979)
      Observations: 150

# display and print method works - HTML

    Code
      print(correlation(subset(mtcars, select = c("wt", "mpg"))), format = "html")
    Output
      Correlation Matrix (pearson-method)
      
      Parameter1 | Parameter2 |     r |         95% CI | t(30) |         p
      --------------------------------------------------------------------
      wt         |        mpg | -0.87 | [-0.93, -0.74] | -9.56 | < .001***
      p-value adjustment method: Holm (1979)Observations: 32

# selecting specific variables works

    Code
      list(df1, df2, df3, df4)
    Output
      [[1]]
      # Correlation Matrix (pearson-method)
      
      Parameter1 | Parameter2 |    r |       95% CI | t(30) |         p
      -----------------------------------------------------------------
      cyl        |         hp | 0.83 | [0.68, 0.92] |  8.23 | < .001***
      wt         |         hp | 0.66 | [0.40, 0.82] |  4.80 | < .001***
      
      p-value adjustment method: Holm (1979)
      Observations: 32
      
      [[2]]
      # Correlation Matrix (pearson-method)
      
      Group | Parameter1 | Parameter2 |    r |       95% CI |    t | df |         p
      -----------------------------------------------------------------------------
      0     |        cyl |         hp | 0.85 | [0.64, 0.94] | 6.53 | 17 | < .001***
      0     |         wt |         hp | 0.68 | [0.33, 0.87] | 3.82 | 17 | 0.001**  
      1     |        cyl |         hp | 0.90 | [0.69, 0.97] | 6.87 | 11 | < .001***
      1     |         wt |         hp | 0.81 | [0.48, 0.94] | 4.66 | 11 | < .001***
      
      p-value adjustment method: Holm (1979)
      Observations: 13-19
      
      [[3]]
      # Correlation Matrix (pearson-method)
      
      Parameter1 | Parameter2 |    r |       95% CI | t(30) |         p
      -----------------------------------------------------------------
      wt         |         hp | 0.66 | [0.40, 0.82] |  4.80 | < .001***
      
      p-value adjustment method: Holm (1979)
      Observations: 32
      
      [[4]]
      # Correlation Matrix (pearson-method)
      
      Parameter1 | Parameter2 |    r |       95% CI | t(30) |         p
      -----------------------------------------------------------------
      wt         |         hp | 0.66 | [0.40, 0.82] |  4.80 | < .001***
      
      p-value adjustment method: Holm (1979)
      Observations: 32
      

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
reported to the community leaders responsible for enforcement at <brenton@wiernik.org>. 
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
# Contribution Guidelines 

<sup>easystats guidelines 0.1.0</sup>

**All people are very much welcome to contribute to code, documentation, testing and suggestions.**

This package aims at being beginner-friendly. Even if you're new to this open-source way of life, new to coding and github stuff, we encourage you to try submitting pull requests (PRs). 

- **"I'd like to help, but I'm not good enough with programming yet"**

It's alright, don't worry! You can always dig in the code, in the documentation or tests. There are always some typos to fix, some docs to improve, some details to add, some code lines to document, some tests to add... **Even the smaller PRs are appreciated**.

- **"I'd like to help, but I don't know where to start"**

You can look around the **issue section** to find some features / ideas / bugs to start working on. You can also open a new issue **just to say that you're there, interested in helping out**. We might have some ideas adapted to your skills.

- **"I'm not sure if my suggestion or idea is worthwile"**

Enough with the impostor syndrom! All suggestions and opinions are good, and even if it's just a thought or so, it's always good to receive feedback.

- **"Why should I waste my time with this? Do I get any credit?"**

Software contributions are getting more and more valued in the academic world, so it is a good time to collaborate with us! Authors of substantial contributions will be added within the **authors** list. We're also very keen on including them to eventual academic publications.


**Anyway, starting is the most important! You will then enter a *whole new world, a new fantastic point of view*... So fork this repo, do some changes and submit them. We will then work together to make the best out of it :)**


## Code

- Please document and comment your code, so that the purpose of each step (or code line) is stated in a clear and understandable way.
- Before submitting a change, please read the [**R style guide**](https://style.tidyverse.org/) and in particular our [**easystats convention of code-style**](https://github.com/easystats/easystats#convention-of-code-style) to keep some consistency in code formatting.
- Regarding the style guide, note this exception: we put readability and clarity before everything. Thus, we like underscores and full names (prefer `model_performance` over `modelperf` and `interpret_odds_logistic` over `intoddslog`).
- Before you start to code, make sure you're on the `dev` branch (the most "advanced"). Then, you can create a new branch named by your feature (e.g., `feature_lightsaber`) and do your changes. Finally, submit your branch to be merged into the `dev` branch. Then, every now and then, the dev branch will merge into `main`, as a new package version.

## Checks to do before submission

- Make sure **documentation** (roxygen) is good
- Make sure to add **tests** for the new functions
- Run:

  - `styler::style_pkg()`: Automatic style formatting
  - `lintr::lint_package()`: Style checks
  - `devtools::check()`: General checks



## Useful Materials

- [Understanding the GitHub flow](https://guides.github.com/introduction/flow/)


## revdepcheck results

We checked 6 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.1.0 (2021-05-18) |
|os       |macOS Mojave 10.14.6         |
|system   |x86_64, darwin17.0           |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |Europe/Berlin                |
|date     |2021-06-23                   |

# Dependencies

|package     |old    |new      |Δ  |
|:-----------|:------|:--------|:--|
|correlation |0.6.1  |0.6.2    |*  |
|bayestestR  |0.10.0 |0.10.5   |*  |
|datawizard  |NA     |0.1.0    |*  |
|effectsize  |0.4.5  |NA       |*  |
|insight     |0.14.2 |0.14.2   |   |
|parameters  |0.14.0 |0.14.0.1 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*