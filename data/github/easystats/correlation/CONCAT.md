
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

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output: github_document
---

# correlation <img src='man/figures/logo.png' align="right" height="139" />

```{r README-1, warning=FALSE, message=FALSE, echo=FALSE}
library(ggplot2)
library(poorman)
library(correlation)

options(digits = 2)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 450,
  message = FALSE,
  warning = FALSE,
  fig.path = "man/figures/"
)
```


[![DOI](https://joss.theoj.org/papers/10.21105/joss.02306/status.svg)](https://doi.org/10.21105/joss.02306)
[![downloads](http://cranlogs.r-pkg.org/badges/correlation)](https://cran.r-project.org/package=correlation)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/correlation)](https://cranlogs.r-pkg.org/) [![status](https://tinyverse.netlify.com/badge/correlation)](https://CRAN.R-project.org/package=correlation) 

`correlation` is an [**easystats**](https://github.com/easystats/easystats) package focused on correlation analysis. It's lightweight, easy to use, and allows for the computation of many different kinds of correlations, such as **partial** correlations, **Bayesian** correlations, **multilevel** correlations, **polychoric** correlations, **biweight**, **percentage bend** or **Sheperd's Pi** correlations (types of robust correlation), **distance** correlation (a type of non-linear correlation) and more, also allowing for combinations between them (for instance, *Bayesian partial multilevel correlation*).

# Citation

You can reference the package and its documentation as follows:

Makowski, D., Ben-Shachar, M. S., Patil, I., \& Lüdecke, D. (2019). Methods and Algorithms for Correlation Analysis in R. _Journal of Open Source Software_,
*5*(51), 2306. https://doi.org/10.21105/joss.02306

# Installation

[![CRAN](http://www.r-pkg.org/badges/version/correlation)](https://cran.r-project.org/package=correlation)
[![Build Status](https://travis-ci.org/easystats/correlation.svg?branch=master)](https://travis-ci.org/easystats/correlation)
[![codecov](https://codecov.io/gh/easystats/correlation/branch/master/graph/badge.svg)](https://codecov.io/gh/easystats/correlation)

Run the following to install the stable release from CRAN:

```{r README-2, eval=FALSE}
install.packages("correlation")
```

Or this one to install the latest development version:

```{r README-3, eval=FALSE}
install.packages("remotes")
remotes::install_github("easystats/correlation")
```

# Documentation

[![Documentation](https://img.shields.io/badge/documentation-correlation-orange.svg?colorB=E91E63)](https://easystats.github.io/correlation/)
[![Blog](https://img.shields.io/badge/blog-easystats-orange.svg?colorB=FF9800)](https://easystats.github.io/blog/posts/)
[![Features](https://img.shields.io/badge/features-correlation-orange.svg?colorB=2196F3)](https://easystats.github.io/correlation/reference/index.html)

Click on the buttons above to access the package [documentation](https://easystats.github.io/correlation/) and the [easystats blog](https://easystats.github.io/blog/posts/), and check-out these vignettes:


- [Types of Correlation](https://easystats.github.io/correlation/articles/types.html)
- [Multilevel Correlations](https://easystats.github.io/correlation/articles/multilevel.html)

# Features

The *correlation* package can compute many different types of correlation,
including:

✅ **Pearson's correlation**<br>
✅ **Spearman's rank correlation**<br>
✅ **Kendall's rank correlation**<br>
✅ **Biweight midcorrelation**<br>
✅ **Distance correlation**<br>
✅ **Percentage bend correlation**<br>
✅ **Shepherd's Pi correlation**<br>
✅ **Blomqvist’s coefficient**<br>
✅ **Hoeffding’s D**<br>
✅ **Gamma correlation**<br>
✅ **Gaussian rank correlation**<br>
✅ **Point-Biserial and biserial correlation**<br>
✅ **Winsorized correlation**<br>
✅ **Polychoric correlation**<br>
✅ **Tetrachoric correlation**<br>
✅ **Multilevel correlation**<br>

An overview and description of these correlations types is [**available here**](https://easystats.github.io/correlation/articles/types.html). Moreover,
many of these correlation types are available as **partial** or within a
**Bayesian** framework.

# Examples

The main function is [`correlation()`](https://easystats.github.io/correlation/reference/correlation.html), which builds on top of [`cor_test()`](https://easystats.github.io/correlation/reference/cor_test.html) and comes with a number of possible options.

## Correlation details and matrix

```{r README-4}
results <- correlation(iris)
results
```

The output is not a square matrix, but a **(tidy) dataframe with all correlations tests per row**. One can also obtain a **matrix** using:

```{r README-5}
summary(results)
```

Note that one can also obtain the full, **square** and redundant matrix using:

```{r README-6}
summary(results, redundant = TRUE)
```


```{r README-7}
library(see)

results %>%
  summary(redundant = TRUE) %>%
  plot()
```

## Correlation tests

The `cor_test()` function, for pairwise correlations, is also very convenient for making quick scatter plots.

```{r README-corr}
plot(cor_test(iris, "Sepal.Width", "Sepal.Length"))
```


## Grouped dataframes

The `correlation()` function also supports **stratified correlations**, all within the
*tidyverse* workflow!

```{r README-8}
iris %>%
  select(Species, Sepal.Length, Sepal.Width, Petal.Width) %>%
  group_by(Species) %>%
  correlation()
```


## Bayesian Correlations

It is very easy to switch to a **Bayesian framework**.

```{r README-9}
correlation(iris, bayesian = TRUE)
```

## Tetrachoric, Polychoric, Biserial, Biweight...

The `correlation` package also supports different types of methods, which can
deal with correlations **between factors**!

```{r README-10}
correlation(iris, include_factors = TRUE, method = "auto")
```


## Partial Correlations

It also supports **partial correlations** (as well as Bayesian partial correlations).

```{r README-11}
iris %>%
  correlation(partial = TRUE) %>%
  summary()
```

## Gaussian Graphical Models (GGMs)

Such partial correlations can also be represented as **Gaussian Graphical
Models** (GGM), an increasingly popular tool in psychology. A GGM traditionally
include a set of variables depicted as circles ("nodes"), and a set of lines
that visualize relationships between them, which thickness represents the
strength of association (see [Bhushan et al., 2019](https://www.frontiersin.org/articles/10.3389/fpsyg.2019.01050/full)).

```{r README-12}
library(see) # for plotting
library(ggraph) # needs to be loaded

mtcars %>%
  correlation(partial = TRUE) %>%
  plot()
```

## Multilevel Correlations

It also provide some cutting-edge methods, such as Multilevel (partial)
correlations. These are are partial correlations based on linear mixed-effects
models that include the factors as **random effects**. They can be see as
correlations *adjusted* for some group (*hierarchical*) variability.

```{r README-13}
iris %>%
  correlation(partial = TRUE, multilevel = TRUE) %>%
  summary()
```

However, if the `partial` argument is set to `FALSE`, it will try to convert the
partial coefficient into regular ones.These can be **converted back** to full
correlations:

```{r README-14}
iris %>%
  correlation(partial = FALSE, multilevel = TRUE) %>%
  summary()
```

# Contributing and Support

In case you want to file an issue or contribute in another way to the package, please follow [this guide](https://easystats.github.io/correlation/CONTRIBUTING.html). For questions about the functionality, you may either contact us via email or also file an issue.

# Code of Conduct

Please note that this project is released with a 
[Contributor Code of Conduct](https://easystats.github.io/correlation/CODE_OF_CONDUCT.html). By participating in this project you agree to abide by its terms.

---
title: "Correlation Types"
output: 
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, correlation, types]
vignette: >
  %\VignetteIndexEntry{Correlation Types}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

```{r, include=FALSE}
library(knitr)
options(knitr.kable.NA = "")
knitr::opts_chunk$set(
  comment = ">",
  out.width = "100%",
  message = FALSE,
  warning = FALSE,
  dpi = 450
)
options(digits = 2)

set.seed(333)


if (!requireNamespace("see", quietly = TRUE) ||
  !requireNamespace("tidyr", quietly = TRUE) ||
  !requireNamespace("poorman", quietly = TRUE) ||
  !requireNamespace("ggplot2", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(see)
  library(tidyr)
  library(poorman)
  library(ggplot2)
}
```

---

This vignette can be cited as:

```{r cite}
citation("correlation")
```

---

## Different Methods for Correlations

Correlations tests are arguably one of the most commonly used statistical
procedures, and are used as a basis in many applications such as exploratory
data analysis, structural modeling, data engineering, etc. In this context, we
present **correlation**, a toolbox for the R language [@Rteam] and part of the
[**easystats**](https://github.com/easystats/easystats) collection, focused on
correlation analysis. Its goal is to be lightweight, easy to use, and allows for
the computation of many different kinds of correlations, such as:

- **Pearson's correlation**: This is the most common correlation method. It
  corresponds to the covariance of the two variables normalized (i.e., divided)
  by the product of their standard deviations.

$$r_{xy} = \frac{cov(x,y)}{SD_x \times SD_y}$$

- **Spearman's rank correlation**: A non-parametric measure of correlation, the
Spearman correlation between two variables is equal to the Pearson correlation
between the rank scores of those two variables; while Pearson's correlation
assesses linear relationships, Spearman's correlation assesses monotonic
relationships (whether linear or not). Confidence Intervals (CI) for Spearman's
correlations are computed using the @fieller1957tests correction [see
@bishara2017confidence].

$$r_{s_{xy}} = \frac{cov(rank_x, rank_y)}{SD(rank_x) \times SD(rank_y)}$$

- **Kendall's rank correlation**: In the normal case, the Kendall correlation is
preferred to the Spearman correlation because of a smaller gross error
sensitivity (GES) and a smaller asymptotic variance (AV), making it more robust
and more efficient. However, the interpretation of Kendall's tau is less direct
compared to that of the Spearman's rho, in the sense that it quantifies the
difference between the % of concordant and discordant pairs among all possible
pairwise events. Confidence Intervals (CI) for Kendall's correlations are
computed using the @fieller1957tests correction [see
@bishara2017confidence]. For each pair of observations (i ,j) of two variables
(x, y), it is defined as follows:

$$\tau_{xy} = \frac{2}{n(n-1)}\sum_{i<j}^{}sign(x_i - x_j) \times sign(y_i - y_j)$$

- **Biweight midcorrelation**: A measure of similarity that is median-based,
instead of the traditional mean-based, thus being less sensitive to outliers. It
can be used as a robust alternative to other similarity metrics, such as Pearson
correlation [@langfelder2012fast].

- **Distance correlation**: Distance correlation measures both linear and
non-linear association between two random variables or random vectors. This is
in contrast to Pearson's correlation, which can only detect linear association
between two random variables.

- **Percentage bend correlation**: Introduced by Wilcox (1994), it is based on a
down-weight of a specified percentage of marginal observations deviating from
the median (by default, 20 percent).

- **Shepherd's Pi correlation**: Equivalent to a Spearman's rank correlation
after outliers removal (by means of bootstrapped Mahalanobis distance).

- **Blomqvist’s coefficient**: The Blomqvist’s coefficient (also referred to as
Blomqvist's Beta or medial correlation; Blomqvist, 1950) is a median-based
non-parametric correlation that has some advantages over measures such as
Spearman's or Kendall's estimates (see Shmid and Schimdt, 2006).

- **Hoeffding’s D**: The Hoeffding’s D statistic is a non-parametric rank based
measure of association that detects more general departures from independence
(Hoeffding 1948), including non-linear associations. Hoeffding’s D varies
between -0.5 and 1 (if there are no tied ranks, otherwise it can have lower
values), with larger values indicating a stronger relationship between the
variables.

- **Gamma correlation**: The Goodman-Kruskal gamma statistic is similar to
Kendall's Tau coefficient. It is relatively robust to outliers and deals well
with data that have many ties.

- **Gaussian rank correlation**: The Gaussian rank correlation estimator is a
simple and well-performing alternative for robust rank correlations (Boudt et
al., 2012). It is based on the Gaussian quantiles of the ranks.

- **Point-Biserial and biserial correlation**: Correlation coefficient used when
one variable is continuous and the other is dichotomous (binary). Point-Biserial
is equivalent to a Pearson's correlation, while Biserial should be used when the
binary variable is assumed to have an underlying continuity. For example,
anxiety level can be measured on a continuous scale, but can be classified
dichotomously as high/low.

- **Winsorized correlation**: Correlation of variables that have been
Winsorized, i.e., transformed by limiting extreme values to reduce the effect of
possibly spurious outliers.

- **Polychoric correlation**: Correlation between two theorised normally
distributed continuous latent variables, from two observed ordinal variables.

- **Tetrachoric correlation**: Special case of the polychoric correlation
applicable when both observed variables are dichotomous.

- **Partial correlation**: Correlation between two variables after adjusting for
the (linear) effect of one or more variables. The correlation test is 
run after having partialized the dataset, independently from it. In other words,
it considers partialization as an independent step generating a different
dataset, rather than belonging to the same model. This is why some discrepancies
are to be expected for the *t*- and the *p*-values (but not the correlation
coefficient) compared to other implementations such as `ppcor`. Let $e_{x.z}$ be
the residuals from the linear prediction of $x$ by $z$ (note that this can be
expanded to a multivariate $z$):

$$r_{xy.z} = r_{e_{x.z},e_{y.z}}$$

- **Multilevel correlation**: Multilevel correlations are a special case of
partial correlations where the variable to be adjusted for is a factor and is
included as a random effect in a mixed-effects model.

## Comparison

We will fit different types of correlations of generated data with different
link strengths and link types.

Let's first load the required libraries for this analysis.

```{r}
library(correlation)
library(bayestestR)
library(see)
library(ggplot2)
library(tidyr)
library(poorman)
```

### Utility functions

```{r}
generate_results <- function(r, n = 100, transformation = "none") {
  data <- bayestestR::simulate_correlation(round(n), r = r)
  
  if (transformation != "none") {
    var <- ifelse(grepl("(", transformation, fixed = TRUE), "data$V2)", "data$V2")
    transformation <- paste0(transformation, var)
    data$V2 <- eval(parse(text = transformation))
  }
  
  out <- data.frame(n = n, transformation = transformation, r = r)

  out$Pearson <- cor_test(data, "V1", "V2", method = "pearson")$r
  out$Spearman <- cor_test(data, "V1", "V2", method = "spearman")$rho
  out$Kendall <- cor_test(data, "V1", "V2", method = "kendall")$tau
  out$Biweight <- cor_test(data, "V1", "V2", method = "biweight")$r
  out$Distance <- cor_test(data, "V1", "V2", method = "distance")$r
  out$Distance <- cor_test(data, "V1", "V2", method = "distance")$r
  
  out
}
```

### Effect of Relationship Type

```{r}
data <- data.frame()
for (r in seq(0, 0.999, length.out = 200)) {
  for (n in c(100)) {
    for (transformation in c(
      "none",
      "exp(",
      "log10(1+max(abs(data$V2))+",
      "1/",
      "tan(",
      "sin(",
      "cos(",
      "cos(2*",
      "abs(",
      "data$V2*",
      "data$V2*data$V2*",
      "ifelse(data$V2>0, 1, 0)*("
    )) {
      data <- rbind(data, generate_results(r, n, transformation = transformation))
    }
  }
}


data %>%
  tidyr::pivot_longer(-c(n, r, transformation),
                      names_to = "Type",
                      values_to = "Estimation") %>% 
  mutate(Type = forcats::fct_relevel(Type, "Pearson", "Spearman", "Kendall", "Biweight", "Distance")) %>%
  ggplot(aes(x = r, y = Estimation, fill = Type)) +
  geom_smooth(aes(color = Type), method = "loess", alpha = 0) +
  geom_vline(aes(xintercept = 0.5), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  see::theme_modern() +
  scale_color_flat_d(palette = "rainbow") +
  scale_fill_flat_d(palette = "rainbow") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~transformation)

model <- data %>%
  tidyr::pivot_longer(-c(n, r, transformation),
                      names_to = "Type",
                      values_to = "Estimation") %>% 
  lm(r ~ Type / Estimation, data = .) %>%
  parameters::parameters()

arrange(model[6:10, ], desc(Coefficient))
```

As we can see, **distance** correlation is able to capture the strength even for
severely non-linear relationships.

# References

---
title: "Multilevel Correlations"
output:
  rmarkdown::html_vignette:
    toc: true
    fig_width: 10.08
    fig_height: 6
tags: [r, correlation, types]
vignette: >
  %\VignetteIndexEntry{Multilevel Correlations}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
bibliography: bibliography.bib
---

---

This vignette can be cited as:

```{r cite}
citation("correlation")
```

---

```{r, include=FALSE}
library(knitr)
options(
  knitr.kable.NA = "",
  digits = 2,
  out.width = "100%",
  message = FALSE,
  warning = FALSE,
  dpi = 450
)

if (!requireNamespace("ggplot2", quietly = TRUE) ||
  !requireNamespace("lme4", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
}
```

## Data

Imagine we have an experiment in which **10 individuals** completed a task with
**100 trials**. For each trial - there will 1000 trials (10 * 1000) in total -
we measured two things, **V1** and **V2**, and we are interested in
**investigating the link between these two variables**.

We will generate data using the
[`simulate_simpson()`](https://easystats.github.io/bayestestR/reference/simulate_simpson.html)
function from this package and look at its summary:

```{r}
library(correlation)

data <- simulate_simpson(n = 100, groups = 10)

summary(data)
```

Now let's visualize the two variables:

```{r}
library(ggplot2)

ggplot(data, aes(x = V1, y = V2)) +
  geom_point() +
  geom_smooth(colour = "black", method = "lm", se = FALSE) +
  theme_classic()
```

That seems pretty straightforward! It seems like there is a **negative
correlation** between V1 and V2. Let's test this.

## Simple correlation

```{r}
correlation(data)
```

Indeed, there is  a **strong, negative and significant correlation** between V1
and V2.

Great, can we go ahead and **publish these results in _PNAS_**?

## The Simpson's Paradox

Not so fast! Ever heard of the [**Simpson's Paradox**](https://en.wikipedia.org/wiki/Simpson%27s_paradox)?

Let's colour our datapoints by group (by individuals):

```{r}
library(ggplot2)

ggplot(data, aes(x = V1, y = V2)) +
  geom_point(aes(colour = Group)) +
  geom_smooth(aes(colour = Group), method = "lm", se = FALSE) +
  geom_smooth(colour = "black", method = "lm", se = FALSE) +
  theme_classic()
```

Mmh, interesting. It seems like, for each subject, the relationship is
different. The (global) negative trend seems to be an artifact of **differences between the groups** and could be spurious!

**Multilevel *(as in multi-group)* ** correlations allow us to account for
**differences between groups**. It is based on a partialization of the group,
entered as a random effect in a mixed linear regression.

You can compute them with the
[**correlations**](https://github.com/easystats/correlation) package by setting
the `multilevel` argument to `TRUE`.

```{r}
correlation(data, multilevel = TRUE)
```

For completeness, let's also see if its Bayesian cousin agrees with it:

```{r}
correlation(data, multilevel = TRUE, bayesian = TRUE)
```

**Dayum!** 
We were too hasty in our conclusions! Taking the group into account
seems to be super important.

_Note_: In this simple case where only two variables are of interest, it would be
of course best to directly proceed using a mixed regression model instead of
correlations. That being said, the latter can be useful for exploratory
analysis, when multiple variables are of interest, or in combination with a
network or structural approach.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R, R/methods_print.R
\name{display.easycormatrix}
\alias{display.easycormatrix}
\alias{print_md.easycorrelation}
\alias{print_html.easycorrelation}
\alias{print_md.easycormatrix}
\alias{print_html.easycormatrix}
\title{Export tables into different output formats}
\usage{
\method{display}{easycormatrix}(
  object,
  format = "markdown",
  digits = 2,
  p_digits = 3,
  stars = TRUE,
  include_significance = NULL,
  ...
)

\method{print_md}{easycorrelation}(x, digits = NULL, p_digits = NULL, stars = NULL, ...)

\method{print_html}{easycorrelation}(x, digits = NULL, p_digits = NULL, stars = NULL, ...)

\method{print_md}{easycormatrix}(
  x,
  digits = NULL,
  p_digits = NULL,
  stars = NULL,
  include_significance = NULL,
  ...
)

\method{print_html}{easycormatrix}(
  x,
  digits = NULL,
  p_digits = NULL,
  stars = NULL,
  include_significance = NULL,
  ...
)
}
\arguments{
\item{object, x}{An object returned by
\code{\link[=correlation]{correlation()}} or its summary.}

\item{format}{String, indicating the output format. Currently, only
\code{"markdown"} is supported.}

\item{digits, p_digits}{To do...}

\item{stars}{To do...}

\item{include_significance}{To do...}

\item{...}{Currently not used.}
}
\value{
A character vector. If \code{format = "markdown"}, the return value
will be a character vector in markdown-table format.
}
\description{
Export tables (i.e. data frame) into different output formats.
\code{print_md()} is a alias for \code{display(format = "markdown")}.
}
\details{
\code{display()} is useful when the table-output from functions,
which is usually printed as formatted text-table to console, should
be formatted for pretty table-rendering in markdown documents, or if
knitted from rmarkdown to PDF or Word files.
}
\examples{
data(iris)
corr <- correlation(iris)
display(corr)

s <- summary(corr)
display(s)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation.R
\name{correlation}
\alias{correlation}
\title{Correlation Analysis}
\usage{
correlation(
  data,
  data2 = NULL,
  select = NULL,
  select2 = NULL,
  rename = NULL,
  method = "pearson",
  p_adjust = "holm",
  ci = 0.95,
  bayesian = FALSE,
  bayesian_prior = "medium",
  bayesian_ci_method = "hdi",
  bayesian_test = c("pd", "rope", "bf"),
  redundant = FALSE,
  include_factors = FALSE,
  partial = FALSE,
  partial_bayesian = FALSE,
  multilevel = FALSE,
  ranktransform = FALSE,
  winsorize = FALSE,
  verbose = TRUE,
  standardize_names = getOption("easystats.standardize_names", FALSE),
  ...
)
}
\arguments{
\item{data}{A data frame.}

\item{data2}{An optional data frame. If specified, all pair-wise correlations
between the variables in \code{data} and \code{data2} will be computed.}

\item{select, select2}{(Ignored if \code{data2} is specified.) Optional names
of variables that should be selected for correlation. Instead of providing
the data frames with those variables that should be correlated, \code{data}
can be a data frame and \code{select} and \code{select2} are (quoted) names
of variables (columns) in \code{data}. \code{correlation()} will then
compute the correlation between \code{data[select]} and
\code{data[select2]}. If only \code{select} is specified, all pairwise
correlations between the \code{select} variables will be computed. This is
a "pipe-friendly" alternative way of using \code{correlation()} (see
'Examples').}

\item{rename}{In case you wish to change the names of the variables in
the output, these arguments can be used to specify these alternative names.
Note that the number of names should be equal to the number of columns
selected. Ignored if \code{data2} is specified.}

\item{method}{A character string indicating which correlation coefficient is
to be used for the test. One of \code{"pearson"} (default),
\code{"kendall"}, \code{"spearman"} (but see also the \code{robust} argument), \code{"biserial"},
\code{"polychoric"}, \code{"tetrachoric"}, \code{"biweight"},
\code{"distance"}, \code{"percentage"} (for percentage bend correlation),
\code{"blomqvist"} (for Blomqvist's coefficient), \code{"hoeffding"} (for
Hoeffding's D), \code{"gamma"}, \code{"gaussian"} (for Gaussian Rank
correlation) or \code{"shepherd"} (for Shepherd's Pi correlation). Setting
\code{"auto"} will attempt at selecting the most relevant method
(polychoric when ordinal factors involved, tetrachoric when dichotomous
factors involved, point-biserial if one dichotomous and one continuous and
pearson otherwise). See below the \strong{details} section for a description of
these indices.}

\item{p_adjust}{Correction method for frequentist correlations. Can be one of
\code{"holm"} (default), \code{"hochberg"}, \code{"hommel"},
\code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"},
\code{"somers"} or \code{"none"}. See
\code{\link[stats:p.adjust]{stats::p.adjust()}} for further details.}

\item{ci}{Confidence/Credible Interval level. If \code{"default"}, then it is
set to \code{0.95} (\verb{95\%} CI).}

\item{bayesian}{If TRUE, will run the correlations under a
Bayesian framework. Note that for partial correlations, you will also need
to set \code{partial_bayesian} to \code{TRUE} to obtain "full" Bayesian
partial correlations. Otherwise, you will obtain pseudo-Bayesian partial
correlations (i.e., Bayesian correlation based on frequentist
partialization).}

\item{bayesian_prior}{For the prior argument, several named values are
recognized: \code{"medium.narrow"}, \code{"medium"}, \code{"wide"}, and
\code{"ultrawide"}. These correspond to scale values of \code{1/sqrt(27)},
\code{1/3}, \code{1/sqrt(3)} and \code{1}, respectively. See the
\code{BayesFactor::correlationBF} function.}

\item{bayesian_ci_method}{See arguments in
\code{\link[=parameters]{model_parameters()}} for \code{BayesFactor} tests.}

\item{bayesian_test}{See arguments in
\code{\link[=parameters]{model_parameters()}} for \code{BayesFactor} tests.}

\item{redundant}{Should the data include redundant rows (where each given
correlation is repeated two times).}

\item{include_factors}{If \code{TRUE}, the factors are kept and eventually
converted to numeric or used as random effects (depending of
\code{multilevel}). If \code{FALSE}, factors are removed upfront.}

\item{partial}{Can be \code{TRUE} or \code{"semi"} for partial and
semi-partial correlations, respectively.}

\item{partial_bayesian}{If TRUE, will run the correlations under a
Bayesian framework. Note that for partial correlations, you will also need
to set \code{partial_bayesian} to \code{TRUE} to obtain "full" Bayesian
partial correlations. Otherwise, you will obtain pseudo-Bayesian partial
correlations (i.e., Bayesian correlation based on frequentist
partialization).}

\item{multilevel}{If \code{TRUE}, the factors are included as random factors.
Else, if \code{FALSE} (default), they are included as fixed effects in the
simple regression model.}

\item{ranktransform}{If \code{TRUE}, will rank-transform the variables prior to
estimating the correlation, which is one way of making the analysis more
resistant to extreme values (outliers). Note that, for instance, a Pearson's
correlation on rank-transformed data is equivalent to a Spearman's rank
correlation. Thus, using \code{robust=TRUE} and \code{method="spearman"} is
redundant. Nonetheless, it is an easy option to increase the robustness of the
correlation as well as flexible way to obtain Bayesian or multilevel
Spearman-like rank correlations.}

\item{winsorize}{Another way of making the correlation more "robust" (i.e.,
limiting the impact of extreme values). Can be either \code{FALSE} or a
number between 0 and 1 (e.g., \code{0.2}) that corresponds to the desired
threshold. See the \code{\link[=winsorize]{winsorize()}} function for more details.}

\item{verbose}{Toggle warnings.}

\item{standardize_names}{This option can be set to \code{TRUE} to run
\code{\link[insight:standardize_names]{insight::standardize_names()}} on the output to get standardized column
names. This option can also be set globally by running
\code{options(easystats.standardize_names = TRUE)}.}

\item{...}{Additional arguments (e.g., \code{alternative}) to be passed to
other methods. See \code{stats::cor.test} for further details.}
}
\value{
A correlation object that can be displayed using the \code{print}, \code{summary} or
\code{table} methods.

\subsection{Multiple tests correction}{
The \code{p_adjust} argument can be used to adjust p-values for multiple
comparisons. All adjustment methods available in \code{p.adjust} function
\code{stats} package are supported.
}
}
\description{
Performs a correlation analysis.
}
\details{
\subsection{Correlation Types}{
\itemize{
\item \strong{Pearson's correlation}: This is the most common correlation
method. It corresponds to the covariance of the two variables normalized
(i.e., divided) by the product of their standard deviations.
\item \strong{Spearman's rank correlation}: A non-parametric measure of rank
correlation (statistical dependence between the rankings of two variables).
The Spearman correlation between two variables is equal to the Pearson
correlation between the rank values of those two variables; while Pearson's
correlation assesses linear relationships, Spearman's correlation assesses
monotonic relationships (whether linear or not). Confidence Intervals (CI)
for Spearman's correlations are computed using the Fieller et al. (1957)
correction (see Bishara and Hittner, 2017).
\item \strong{Kendall's rank correlation}: In the normal case, the Kendall correlation
is preferred than the Spearman correlation because of a smaller gross error
sensitivity (GES) and a smaller asymptotic variance (AV), making it more
robust and more efficient. However, the interpretation of Kendall's tau is
less direct than that of Spearman's rho, in the sense that it quantifies the
difference between the percentage of concordant and discordant pairs among
all possible pairwise events. Confidence Intervals (CI) for Kendall's
correlations are computed using the Fieller et al. (1957) correction (see
Bishara and Hittner, 2017).
\item \strong{Biweight midcorrelation}: A measure of similarity that is
median-based, instead of the traditional mean-based, thus being less
sensitive to outliers. It can be used as a robust alternative to other
similarity metrics, such as Pearson correlation (Langfelder & Horvath,
2012).
\item \strong{Distance correlation}: Distance correlation measures both
linear and non-linear association between two random variables or random
vectors. This is in contrast to Pearson's correlation, which can only detect
linear association between two random variables.
\item \strong{Percentage bend correlation}: Introduced by Wilcox (1994), it
is based on a down-weight of a specified percentage of marginal observations
deviating from the median (by default, \verb{20\%}).
\item \strong{Shepherd's Pi correlation}: Equivalent to a Spearman's rank
correlation after outliers removal (by means of bootstrapped Mahalanobis
distance).
\item \strong{Blomqvist’s coefficient}: The Blomqvist’s coefficient (also
referred to as Blomqvist's Beta or medial correlation; Blomqvist, 1950) is a
median-based non-parametric correlation that has some advantages over
measures such as Spearman's or Kendall's estimates (see Shmid & Schimdt,
2006).
\item \strong{Hoeffding’s D}: The Hoeffding’s D statistics is a
non-parametric rank based measure of association that detects more general
departures from independence (Hoeffding 1948), including non-linear
associations. Hoeffding’s D varies between -0.5 and 1 (if there are no tied
ranks, otherwise it can have lower values), with larger values indicating a
stronger relationship between the variables.
\item \strong{Somers’ D}: The Somers’ D statistics is a non-parametric rank
based measure of association between a binary variable and a continuous
variable, for instance, in the context of logistic regression the binary
outcome and the predicted probabilities for each outcome. Usually, Somers' D
is a measure of ordinal association, however, this implementation it is
limited to the case of a binary outcome.
\item \strong{Point-Biserial and biserial correlation}: Correlation
coefficient used when one variable is continuous and the other is dichotomous
(binary). Point-Biserial is equivalent to a Pearson's correlation, while
Biserial should be used when the binary variable is assumed to have an
underlying continuity. For example, anxiety level can be measured on a
continuous scale, but can be classified dichotomously as high/low.
\item \strong{Gamma correlation}: The Goodman-Kruskal gamma statistic is
similar to Kendall's Tau coefficient. It is relatively robust to outliers and
deals well with data that have many ties.
\item \strong{Winsorized correlation}: Correlation of variables that have
been formerly Winsorized, i.e., transformed by limiting extreme values to
reduce the effect of possibly spurious outliers.
\item \strong{Gaussian rank Correlation}: The Gaussian rank correlation
estimator is a simple and well-performing alternative for robust rank
correlations (Boudt et al., 2012). It is based on the Gaussian quantiles of
the ranks.
\item \strong{Polychoric correlation}: Correlation between two theorized
normally distributed continuous latent variables, from two observed ordinal
variables.
\item \strong{Tetrachoric correlation}: Special case of the polychoric
correlation applicable when both observed variables are dichotomous.
}
}

\subsection{Partial Correlation}{
\strong{Partial correlations} are estimated as the correlation between two
variables after adjusting for the (linear) effect of one or more other
variable. The correlation test is then run after having partialized the
dataset, independently from it. In other words, it considers partialization
as an independent step generating a different dataset, rather than belonging
to the same model. This is why some discrepancies are to be expected for the
t- and p-values, CIs, BFs etc (but \emph{not} the correlation coefficient)
compared to other implementations (e.g., \code{ppcor}). (The size of these
discrepancies depends on the number of covariates partialled-out and the
strength of the linear association between all variables.) Such partial
correlations can be represented as Gaussian Graphical Models (GGM), an
increasingly popular tool in psychology. A GGM traditionally include a set of
variables depicted as circles ("nodes"), and a set of lines that visualize
relationships between them, which thickness represents the strength of
association (see Bhushan et al., 2019).

\strong{Multilevel correlations} are a special case of partial correlations where
the variable to be adjusted for is a factor and is included as a random
effect in a mixed model (note that the remaining continuous variables of the
dataset will still be included as fixed effects, similarly to regular partial
correlations). That said, there is an important difference between using
\code{cor_test()} and \code{correlation()}: If you set \code{multilevel=TRUE} in
\code{correlation()} but \code{partial} is set to \code{FALSE} (as per default), then a
back-transformation from partial to non-partial correlation will be attempted
(through \code{\link[=pcor_to_cor]{pcor_to_cor()}}). However, this is not possible when
using \code{cor_test()} so that if you set \code{multilevel=TRUE} in it, the resulting
correlations are partial one. Note that for Bayesian multilevel correlations,
if \code{partial = FALSE}, the back transformation will also recompute \emph{p}-values
based on the new \emph{r} scores, and will drop the Bayes factors (as they are not
relevant anymore). To keep Bayesian scores, set \code{partial = TRUE}.
}

\subsection{Notes}{
Kendall and Spearman correlations when \code{bayesian=TRUE}: These are technically
Pearson Bayesian correlations of rank transformed data, rather than pure
Bayesian rank correlations (which have different priors).
}
}
\examples{

library(correlation)
results <- correlation(iris)

results
summary(results)
summary(results, redundant = TRUE)

# pipe-friendly usage with  grouped dataframes from {dplyr} package
if (require("poorman")) {
  iris \%>\%
    correlation(select = "Petal.Width", select2 = "Sepal.Length")

  # Grouped dataframe
  # grouped correlations
  iris \%>\%
    group_by(Species) \%>\%
    correlation()

  # selecting specific variables for correlation
  mtcars \%>\%
    group_by(am) \%>\%
    correlation(
      select = c("cyl", "wt"),
      select2 = c("hp")
    )
}

# supplying custom variable names
correlation(anscombe, select = c("x1", "x2"), rename = c("var1", "var2"))

# automatic selection of correlation method
correlation(mtcars[-2], method = "auto")

}
\references{
\itemize{
\item Boudt, K., Cornelissen, J., & Croux, C. (2012). The Gaussian rank
correlation estimator: robustness properties. Statistics and Computing,
22(2), 471-483.
\item Bhushan, N., Mohnert, F., Sloot, D., Jans, L., Albers, C., & Steg, L.
(2019). Using a Gaussian graphical model to explore relationships between
items and variables in environmental psychology research. Frontiers in
psychology, 10, 1050.
\item Bishara, A. J., & Hittner, J. B. (2017). Confidence intervals for
correlations when data are not normal. Behavior research methods, 49(1),
294-309.
\item Fieller, E. C., Hartley, H. O., & Pearson, E. S. (1957). Tests for
rank correlation coefficients. I. Biometrika, 44(3/4), 470-481.
\item Langfelder, P., & Horvath, S. (2012). Fast R functions for robust
correlations and hierarchical clustering. Journal of statistical software,
46(11).
\item Blomqvist, N. (1950). On a measure of dependence between two random
variables,Annals of Mathematical Statistics,21, 593–600
\item Somers, R. H. (1962). A new asymmetric measure of association for
ordinal variables. American Sociological Review. 27 (6).
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_smooth.R
\name{cor_smooth}
\alias{cor_smooth}
\alias{is.positive_definite}
\alias{is_positive_definite}
\title{Smooth a non-positive definite correlation matrix to make it positive definite}
\usage{
cor_smooth(x, method = "psych", verbose = TRUE, ...)

is.positive_definite(x, tol = 10^-12, ...)

is_positive_definite(x, tol = 10^-12, ...)
}
\arguments{
\item{x}{A correlation matrix.}

\item{method}{Smoothing method. Can be \code{psych} (will use
\code{psych::cor.smooth()}), \code{hj} (Jorjani et al., 2003) or \code{lrs} (Schaeffer,
2014). For the two last, will use \code{mbend::bend()} (check its documentation
for details).}

\item{verbose}{Set to \code{FALSE} to silence the function.}

\item{...}{Other arguments to be passed to or from other functions.}

\item{tol}{The minimum eigenvalue to be considered as acceptable.}
}
\description{
Make correlations positive definite using \code{psych::cor.smooth}. If smoothing
is done, inferential statistics (\emph{p}-values, confidence intervals, etc.) are
removed, as they are no longer valid.
}
\examples{
set.seed(123)
data <- as.matrix(mtcars)
# Make missing data so pairwise correlation matrix is non-positive definite
data[sample(seq_len(352), size = 60)] <- NA
data <- as.data.frame(data)
x <- correlation(data)
is.positive_definite(x)

smoothed <- cor_smooth(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_inverse.R
\name{matrix_inverse}
\alias{matrix_inverse}
\title{Matrix Inversion}
\usage{
matrix_inverse(m, tol = .Machine$double.eps^(2/3))
}
\arguments{
\item{m}{Matrix for which the inverse is required.}

\item{tol}{Relative tolerance to detect zero singular values.}
}
\value{
An inversed matrix.
}
\description{
Performs a Moore-Penrose generalized inverse (also called the Pseudoinverse).
}
\examples{
m <- cor(iris[1:4])
matrix_inverse(m)
}
\seealso{
pinv from the pracma package
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_to_cov.R
\name{cor_to_cov}
\alias{cor_to_cov}
\title{Convert a correlation to covariance}
\usage{
cor_to_cov(cor, sd = NULL, variance = NULL, tol = .Machine$double.eps^(2/3))
}
\arguments{
\item{cor}{A correlation matrix, or a partial or a semipartial
correlation matrix.}

\item{sd, variance}{A vector that contains the standard deviations, or the
variance, of the variables in the correlation matrix.}

\item{tol}{Relative tolerance to detect zero singular values.}
}
\value{
A covariance matrix.
}
\description{
Convert a correlation to covariance
}
\examples{
cor <- cor(iris[1:4])
cov(iris[1:4])

cor_to_cov(cor, sd = sapply(iris[1:4], sd))
cor_to_cov(cor, variance = sapply(iris[1:4], var))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{isSquare}
\alias{isSquare}
\title{Check if Square Matrix}
\usage{
isSquare(m)
}
\arguments{
\item{m}{A matrix.}
}
\value{
\code{TRUE} of the matrix is square or \code{FALSE} otherwise.
}
\description{
Check if Square Matrix
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_to_pcor.R, R/cor_to_spcor.R
\name{cor_to_pcor}
\alias{cor_to_pcor}
\alias{pcor_to_cor}
\alias{cor_to_spcor}
\title{Correlation Matrix to (Semi) Partial Correlations}
\usage{
cor_to_pcor(cor, tol = .Machine$double.eps^(2/3))

pcor_to_cor(pcor, tol = .Machine$double.eps^(2/3))

cor_to_spcor(cor = NULL, cov = NULL, tol = .Machine$double.eps^(2/3))
}
\arguments{
\item{cor, pcor}{A correlation matrix, or a partial or a semipartial
correlation matrix.}

\item{tol}{Relative tolerance to detect zero singular values.}

\item{cov}{A covariance matrix (or a vector of the SD of the variables).
Required for semi-partial correlations.}
}
\value{
The (semi) partial correlation matrix.
}
\description{
Convert a correlation matrix to a (semi)partial correlation matrix. Partial
correlations are a measure of the correlation between two variables that
remains after controlling for (i.e., "partialling" out) all the other
relationships. They can be used for graphical Gaussian models, as they
represent the direct interactions between two variables, conditioned on all
remaining variables. This means that the squared partial correlation between
a predictor X1 and a response variable Y can be interpreted as the proportion
of (unique) variance accounted for by X1 relative to the residual or
unexplained variance of Y that cannot be accounted by the other variables.
}
\details{
The semi-partial correlation is similar to the partial correlation statistic.
However, it represents (when squared) the proportion of (unique) variance
accounted for by the predictor X1, relative to the total variance of Y. Thus,
it might be seen as a better indicator of the "practical relevance" of a
predictor, because it is scaled to (i.e., relative to) the total variability
in the response variable.
}
\examples{
cor <- cor(iris[1:4])

# Partialize
cor_to_pcor(cor)
cor_to_spcor(cor, cov = sapply(iris[1:4], sd))

# Inverse
round(pcor_to_cor(cor_to_pcor(cor)) - cor, 2) # Should be 0
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is.cor}
\alias{is.cor}
\title{Check if matrix ressembles a correlation matrix}
\usage{
is.cor(x)
}
\arguments{
\item{x}{A matrix.}
}
\value{
\code{TRUE} of the matrix is a correlation matrix or \code{FALSE} otherwise.
}
\description{
Check if matrix ressembles a correlation matrix
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_sort.R
\name{cor_sort}
\alias{cor_sort}
\title{Sort a correlation matrix to improve readability of groups and clusters}
\usage{
cor_sort(x, distance = "correlation", ...)
}
\arguments{
\item{x}{A correlation matrix.}

\item{distance}{How the distance between each variable should be calculated.
If \code{correlation} (default; suited for correlation matrices), the matrix
will be rescaled to 0-1 (\code{distance = 0} indicating correlation of 1;
\code{distance = 1} indicating correlation of -1). If \code{raw}, then the matrix
will be used as a distance matrix as-is. Can be others (\code{euclidean},
\code{manhattan}, ...), in which case it will be passed to \code{dist()} (see the
arguments for it).}

\item{...}{Other arguments to be passed to or from other functions.}
}
\description{
Sort a correlation matrix based on \code{hclust}.
}
\examples{
x <- correlation(mtcars)

cor_sort(as.matrix(x))
cor_sort(x, hclust_method = "ward.D2") # It can also reorder the long form output
cor_sort(summary(x, redundant = TRUE)) # As well as from the summary
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/z_fisher.R
\name{z_fisher}
\alias{z_fisher}
\title{Fisher z-transformation}
\usage{
z_fisher(r = NULL, z = NULL)
}
\arguments{
\item{r, z}{The r or the z' value to be converted.}
}
\value{
The transformed value.
}
\description{
The Fisher z-transformation converts the standard Pearson's \emph{r} to a normally
distributed variable z'. It is used to compute confidence intervals to
correlations. The z' variable is different from the \emph{z}-statistic.
}
\examples{
z_fisher(r = 0.7)
z_fisher(z = 0.867)

}
\references{
Zar, J.H., (2014). Spearman Rank Correlation: Overview. Wiley StatsRef:
Statistics Reference Online. doi:10.1002/9781118445112.stat05964
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_to_ci.R, R/cor_to_p.R
\name{cor_to_ci}
\alias{cor_to_ci}
\alias{cor_to_p}
\title{Convert correlation to p-values and CIs}
\usage{
cor_to_ci(cor, n, ci = 0.95, method = "pearson", correction = "fieller", ...)

cor_to_p(cor, n, method = "pearson")
}
\arguments{
\item{cor}{A correlation matrix or coefficient.}

\item{n}{The sample size (number of observations).}

\item{ci}{Confidence/Credible Interval level. If \code{"default"}, then it is
set to \code{0.95} (\verb{95\%} CI).}

\item{method}{A character string indicating which correlation coefficient is
to be used for the test. One of \code{"pearson"} (default),
\code{"kendall"}, \code{"spearman"} (but see also the \code{robust} argument), \code{"biserial"},
\code{"polychoric"}, \code{"tetrachoric"}, \code{"biweight"},
\code{"distance"}, \code{"percentage"} (for percentage bend correlation),
\code{"blomqvist"} (for Blomqvist's coefficient), \code{"hoeffding"} (for
Hoeffding's D), \code{"gamma"}, \code{"gaussian"} (for Gaussian Rank
correlation) or \code{"shepherd"} (for Shepherd's Pi correlation). Setting
\code{"auto"} will attempt at selecting the most relevant method
(polychoric when ordinal factors involved, tetrachoric when dichotomous
factors involved, point-biserial if one dichotomous and one continuous and
pearson otherwise). See below the \strong{details} section for a description of
these indices.}

\item{correction}{Only used if method is 'spearman' or 'kendall'. Can be
'fieller' (default; Fieller et al., 1957), 'bw' (only for Spearman) or
'none'. Bonett and Wright (2000) claim their correction ('bw') performs
better, though the Bishara and Hittner (2017) paper favours the Fieller
correction. Both are generally very similar.}

\item{...}{Additional arguments (e.g., \code{alternative}) to be passed to
other methods. See \code{stats::cor.test} for further details.}
}
\value{
A list containing a \emph{p}-value and the statistic or the CI bounds.
}
\description{
Get statistics, \emph{p}-values and confidence intervals (CI) from correlation
coefficients.
}
\examples{
cor.test(iris$Sepal.Length, iris$Sepal.Width)
cor_to_p(-0.1175698, n = 150)
cor_to_p(cor(iris[1:4]), n = 150)
cor_to_ci(-0.1175698, n = 150)
cor_to_ci(cor(iris[1:4]), n = 150)

cor.test(iris$Sepal.Length, iris$Sepal.Width, method = "spearman")
cor_to_p(-0.1667777, n = 150, method = "spearman")
cor_to_ci(-0.1667777, ci = 0.95, n = 150)

cor.test(iris$Sepal.Length, iris$Sepal.Width, method = "kendall")
cor_to_p(-0.07699679, n = 150, method = "kendall")
}
\references{
Bishara, A. J., & Hittner, J. B. (2017). Confidence intervals for
correlations when data are not normal. Behavior research methods, 49(1),
294-309.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation_recipe.cor_test.R,
%   R/visualisation_recipe.easycormatrix.R,
%   R/visualisation_recipe.easycorrelation.R
\name{visualisation_recipe.easycor_test}
\alias{visualisation_recipe.easycor_test}
\alias{visualisation_recipe.easycormatrix}
\alias{visualisation_recipe.easycorrelation}
\title{Visualisation Recipe for 'correlation' Objects}
\usage{
\method{visualisation_recipe}{easycor_test}(
  x,
  show_data = "point",
  show_text = "subtitle",
  smooth = NULL,
  point = NULL,
  text = NULL,
  labs = NULL,
  ...
)

\method{visualisation_recipe}{easycormatrix}(
  x,
  show_data = "tile",
  show_text = "text",
  show_legend = TRUE,
  tile = NULL,
  point = NULL,
  text = NULL,
  scale = NULL,
  scale_fill = NULL,
  labs = NULL,
  type = show_data,
  ...
)

\method{visualisation_recipe}{easycorrelation}(x, ...)
}
\arguments{
\item{x}{A correlation object.}

\item{show_data}{Show data. For correlation matrices, can be \code{"tile"}
(default) or \code{"point"}.}

\item{show_text}{Show labels with matrix values.}

\item{...}{Other arguments passed to other functions.}

\item{show_legend}{Show legend. Can be set to \code{FALSE} to remove the legend.}

\item{tile, point, text, scale, scale_fill, smooth, labs}{Additional aesthetics and
parameters for the geoms (see customization example).}

\item{type}{Alias for \code{show_data}, for backwards compatibility.}
}
\description{
Visualisation Recipe for 'correlation' Objects
}
\examples{
\donttest{
# ==============================================
# Correlation Test
# ==============================================
if (require("see")) {
  rez <- cor_test(mtcars, "mpg", "wt")

  layers <- visualisation_recipe(rez, labs = list(x = "Miles per Gallon (mpg)"))
  layers
  plot(layers)

  plot(rez,
    show_text = "label",
    point = list(color = "#f44336"),
    text = list(fontface = "bold"),
    show_statistic = FALSE, show_ci = FALSE, stars = TRUE
  )
}
}
# ==============================================
# Correlation Matrix
# ==============================================
if (require("see")) {
  rez <- correlation(mtcars)

  x <- cor_sort(as.matrix(rez))
  layers <- visualisation_recipe(x)
  layers
  plot(layers)

  #' Get more details using `summary()`
  x <- summary(rez, redundant = TRUE, digits = 3)
  plot(visualisation_recipe(x))

  # Customize
  x <- summary(rez)
  layers <- visualisation_recipe(x,
    show_data = "points",
    scale = list(range = c(10, 20)),
    scale_fill = list(
      high = "#FF5722",
      low = "#673AB7",
      name = "r"
    ),
    text = list(color = "white"),
    labs = list(title = "My Plot")
  )
  plot(layers) + theme_modern()
}
\donttest{
if (FALSE) {
  # ==============================================
  # Correlation Results (easycorrelation)
  # ==============================================
  if (require("see") && require("tidygraph") && require("ggraph")) {
    rez <- correlation(iris)

    layers <- visualisation_recipe(rez)
    layers
    plot(layers)
  }
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexports.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{simulate_simpson}
\alias{visualisation_recipe}
\alias{print_md}
\alias{print_html}
\alias{display}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{bayestestR}{\code{\link[bayestestR]{simulate_simpson}}}

  \item{datawizard}{\code{\link[datawizard]{visualisation_recipe}}}

  \item{insight}{\code{\link[insight]{display}}, \code{\link[insight:display]{print_html}}, \code{\link[insight:display]{print_md}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_test.R
\name{cor_test}
\alias{cor_test}
\title{Correlation test}
\usage{
cor_test(
  data,
  x,
  y,
  method = "pearson",
  ci = 0.95,
  bayesian = FALSE,
  bayesian_prior = "medium",
  bayesian_ci_method = "hdi",
  bayesian_test = c("pd", "rope", "bf"),
  include_factors = FALSE,
  partial = FALSE,
  partial_bayesian = FALSE,
  multilevel = FALSE,
  ranktransform = FALSE,
  winsorize = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{data}{A data frame.}

\item{x, y}{Names of two variables present in the data.}

\item{method}{A character string indicating which correlation coefficient is
to be used for the test. One of \code{"pearson"} (default),
\code{"kendall"}, \code{"spearman"} (but see also the \code{robust} argument), \code{"biserial"},
\code{"polychoric"}, \code{"tetrachoric"}, \code{"biweight"},
\code{"distance"}, \code{"percentage"} (for percentage bend correlation),
\code{"blomqvist"} (for Blomqvist's coefficient), \code{"hoeffding"} (for
Hoeffding's D), \code{"gamma"}, \code{"gaussian"} (for Gaussian Rank
correlation) or \code{"shepherd"} (for Shepherd's Pi correlation). Setting
\code{"auto"} will attempt at selecting the most relevant method
(polychoric when ordinal factors involved, tetrachoric when dichotomous
factors involved, point-biserial if one dichotomous and one continuous and
pearson otherwise). See below the \strong{details} section for a description of
these indices.}

\item{ci}{Confidence/Credible Interval level. If \code{"default"}, then it is
set to \code{0.95} (\verb{95\%} CI).}

\item{bayesian, partial_bayesian}{If TRUE, will run the correlations under a
Bayesian framework. Note that for partial correlations, you will also need
to set \code{partial_bayesian} to \code{TRUE} to obtain "full" Bayesian
partial correlations. Otherwise, you will obtain pseudo-Bayesian partial
correlations (i.e., Bayesian correlation based on frequentist
partialization).}

\item{bayesian_prior}{For the prior argument, several named values are
recognized: \code{"medium.narrow"}, \code{"medium"}, \code{"wide"}, and
\code{"ultrawide"}. These correspond to scale values of \code{1/sqrt(27)},
\code{1/3}, \code{1/sqrt(3)} and \code{1}, respectively. See the
\code{BayesFactor::correlationBF} function.}

\item{bayesian_ci_method, bayesian_test}{See arguments in
\code{\link[=parameters]{model_parameters()}} for \code{BayesFactor} tests.}

\item{include_factors}{If \code{TRUE}, the factors are kept and eventually
converted to numeric or used as random effects (depending of
\code{multilevel}). If \code{FALSE}, factors are removed upfront.}

\item{partial}{Can be \code{TRUE} or \code{"semi"} for partial and
semi-partial correlations, respectively.}

\item{multilevel}{If \code{TRUE}, the factors are included as random factors.
Else, if \code{FALSE} (default), they are included as fixed effects in the
simple regression model.}

\item{ranktransform}{If \code{TRUE}, will rank-transform the variables prior to
estimating the correlation, which is one way of making the analysis more
resistant to extreme values (outliers). Note that, for instance, a Pearson's
correlation on rank-transformed data is equivalent to a Spearman's rank
correlation. Thus, using \code{robust=TRUE} and \code{method="spearman"} is
redundant. Nonetheless, it is an easy option to increase the robustness of the
correlation as well as flexible way to obtain Bayesian or multilevel
Spearman-like rank correlations.}

\item{winsorize}{Another way of making the correlation more "robust" (i.e.,
limiting the impact of extreme values). Can be either \code{FALSE} or a
number between 0 and 1 (e.g., \code{0.2}) that corresponds to the desired
threshold. See the \code{\link[=winsorize]{winsorize()}} function for more details.}

\item{verbose}{Toggle warnings.}

\item{...}{Additional arguments (e.g., \code{alternative}) to be passed to
other methods. See \code{stats::cor.test} for further details.}
}
\description{
This function performs a correlation test between two variables.
}
\details{
\subsection{Correlation Types}{
\itemize{
\item \strong{Pearson's correlation}: This is the most common correlation
method. It corresponds to the covariance of the two variables normalized
(i.e., divided) by the product of their standard deviations.
\item \strong{Spearman's rank correlation}: A non-parametric measure of rank
correlation (statistical dependence between the rankings of two variables).
The Spearman correlation between two variables is equal to the Pearson
correlation between the rank values of those two variables; while Pearson's
correlation assesses linear relationships, Spearman's correlation assesses
monotonic relationships (whether linear or not). Confidence Intervals (CI)
for Spearman's correlations are computed using the Fieller et al. (1957)
correction (see Bishara and Hittner, 2017).
\item \strong{Kendall's rank correlation}: In the normal case, the Kendall correlation
is preferred than the Spearman correlation because of a smaller gross error
sensitivity (GES) and a smaller asymptotic variance (AV), making it more
robust and more efficient. However, the interpretation of Kendall's tau is
less direct than that of Spearman's rho, in the sense that it quantifies the
difference between the percentage of concordant and discordant pairs among
all possible pairwise events. Confidence Intervals (CI) for Kendall's
correlations are computed using the Fieller et al. (1957) correction (see
Bishara and Hittner, 2017).
\item \strong{Biweight midcorrelation}: A measure of similarity that is
median-based, instead of the traditional mean-based, thus being less
sensitive to outliers. It can be used as a robust alternative to other
similarity metrics, such as Pearson correlation (Langfelder & Horvath,
2012).
\item \strong{Distance correlation}: Distance correlation measures both
linear and non-linear association between two random variables or random
vectors. This is in contrast to Pearson's correlation, which can only detect
linear association between two random variables.
\item \strong{Percentage bend correlation}: Introduced by Wilcox (1994), it
is based on a down-weight of a specified percentage of marginal observations
deviating from the median (by default, \verb{20\%}).
\item \strong{Shepherd's Pi correlation}: Equivalent to a Spearman's rank
correlation after outliers removal (by means of bootstrapped Mahalanobis
distance).
\item \strong{Blomqvist’s coefficient}: The Blomqvist’s coefficient (also
referred to as Blomqvist's Beta or medial correlation; Blomqvist, 1950) is a
median-based non-parametric correlation that has some advantages over
measures such as Spearman's or Kendall's estimates (see Shmid & Schimdt,
2006).
\item \strong{Hoeffding’s D}: The Hoeffding’s D statistics is a
non-parametric rank based measure of association that detects more general
departures from independence (Hoeffding 1948), including non-linear
associations. Hoeffding’s D varies between -0.5 and 1 (if there are no tied
ranks, otherwise it can have lower values), with larger values indicating a
stronger relationship between the variables.
\item \strong{Somers’ D}: The Somers’ D statistics is a non-parametric rank
based measure of association between a binary variable and a continuous
variable, for instance, in the context of logistic regression the binary
outcome and the predicted probabilities for each outcome. Usually, Somers' D
is a measure of ordinal association, however, this implementation it is
limited to the case of a binary outcome.
\item \strong{Point-Biserial and biserial correlation}: Correlation
coefficient used when one variable is continuous and the other is dichotomous
(binary). Point-Biserial is equivalent to a Pearson's correlation, while
Biserial should be used when the binary variable is assumed to have an
underlying continuity. For example, anxiety level can be measured on a
continuous scale, but can be classified dichotomously as high/low.
\item \strong{Gamma correlation}: The Goodman-Kruskal gamma statistic is
similar to Kendall's Tau coefficient. It is relatively robust to outliers and
deals well with data that have many ties.
\item \strong{Winsorized correlation}: Correlation of variables that have
been formerly Winsorized, i.e., transformed by limiting extreme values to
reduce the effect of possibly spurious outliers.
\item \strong{Gaussian rank Correlation}: The Gaussian rank correlation
estimator is a simple and well-performing alternative for robust rank
correlations (Boudt et al., 2012). It is based on the Gaussian quantiles of
the ranks.
\item \strong{Polychoric correlation}: Correlation between two theorized
normally distributed continuous latent variables, from two observed ordinal
variables.
\item \strong{Tetrachoric correlation}: Special case of the polychoric
correlation applicable when both observed variables are dichotomous.
}
}

\subsection{Partial Correlation}{
\strong{Partial correlations} are estimated as the correlation between two
variables after adjusting for the (linear) effect of one or more other
variable. The correlation test is then run after having partialized the
dataset, independently from it. In other words, it considers partialization
as an independent step generating a different dataset, rather than belonging
to the same model. This is why some discrepancies are to be expected for the
t- and p-values, CIs, BFs etc (but \emph{not} the correlation coefficient)
compared to other implementations (e.g., \code{ppcor}). (The size of these
discrepancies depends on the number of covariates partialled-out and the
strength of the linear association between all variables.) Such partial
correlations can be represented as Gaussian Graphical Models (GGM), an
increasingly popular tool in psychology. A GGM traditionally include a set of
variables depicted as circles ("nodes"), and a set of lines that visualize
relationships between them, which thickness represents the strength of
association (see Bhushan et al., 2019).

\strong{Multilevel correlations} are a special case of partial correlations where
the variable to be adjusted for is a factor and is included as a random
effect in a mixed model (note that the remaining continuous variables of the
dataset will still be included as fixed effects, similarly to regular partial
correlations). That said, there is an important difference between using
\code{cor_test()} and \code{correlation()}: If you set \code{multilevel=TRUE} in
\code{correlation()} but \code{partial} is set to \code{FALSE} (as per default), then a
back-transformation from partial to non-partial correlation will be attempted
(through \code{\link[=pcor_to_cor]{pcor_to_cor()}}). However, this is not possible when
using \code{cor_test()} so that if you set \code{multilevel=TRUE} in it, the resulting
correlations are partial one. Note that for Bayesian multilevel correlations,
if \code{partial = FALSE}, the back transformation will also recompute \emph{p}-values
based on the new \emph{r} scores, and will drop the Bayes factors (as they are not
relevant anymore). To keep Bayesian scores, set \code{partial = TRUE}.
}

\subsection{Notes}{
Kendall and Spearman correlations when \code{bayesian=TRUE}: These are technically
Pearson Bayesian correlations of rank transformed data, rather than pure
Bayesian rank correlations (which have different priors).
}
}
\examples{
library(correlation)

cor_test(iris, "Sepal.Length", "Sepal.Width")
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "spearman")
\dontrun{
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "kendall")
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "biweight")
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "distance")
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "percentage")
if (require("wdm", quietly = TRUE)) {
  cor_test(iris, "Sepal.Length", "Sepal.Width", method = "blomqvist")
}
if (require("Hmisc", quietly = TRUE)) {
  cor_test(iris, "Sepal.Length", "Sepal.Width", method = "hoeffding")
}
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "gamma")
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "gaussian")
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "shepherd")
if (require("BayesFactor", quietly = TRUE)) {
  cor_test(iris, "Sepal.Length", "Sepal.Width", bayesian = TRUE)
}

# Robust (these two are equivalent)
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "spearman")
cor_test(iris, "Sepal.Length", "Sepal.Width", method = "pearson", ranktransform = TRUE)

# Winsorized
cor_test(iris, "Sepal.Length", "Sepal.Width", winsorize = 0.2)

# Tetrachoric
if (require("psych", quietly = TRUE)) {
  data <- iris
  data$Sepal.Width_binary <- ifelse(data$Sepal.Width > 3, 1, 0)
  data$Petal.Width_binary <- ifelse(data$Petal.Width > 1.2, 1, 0)
  cor_test(data, "Sepal.Width_binary", "Petal.Width_binary", method = "tetrachoric")

  # Biserial
  cor_test(data, "Sepal.Width", "Petal.Width_binary", method = "biserial")

  # Polychoric
  data$Petal.Width_ordinal <- as.factor(round(data$Petal.Width))
  data$Sepal.Length_ordinal <- as.factor(round(data$Sepal.Length))
  cor_test(data, "Petal.Width_ordinal", "Sepal.Length_ordinal", method = "polychoric")

  # When one variable is continuous, will run 'polyserial' correlation
  cor_test(data, "Sepal.Width", "Sepal.Length_ordinal", method = "polychoric")
}

# Partial
cor_test(iris, "Sepal.Length", "Sepal.Width", partial = TRUE)
cor_test(iris, "Sepal.Length", "Sepal.Width", multilevel = TRUE)
cor_test(iris, "Sepal.Length", "Sepal.Width", partial_bayesian = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mahalanobis.R
\name{distance_mahalanobis}
\alias{distance_mahalanobis}
\title{Mahalanobis distance and confidence interval (CI)}
\usage{
distance_mahalanobis(data, ci = 0.95, iterations = 1000, robust = TRUE, ...)
}
\arguments{
\item{data}{A data frame.}

\item{ci}{Confidence/Credible Interval level. If \code{"default"}, then it is
set to \code{0.95} (\verb{95\%} CI).}

\item{iterations}{The number of draws to simulate/bootstrap (when
\code{robust} is \code{TRUE}).}

\item{robust}{If \code{TRUE}, will run a bootstrapped version of the function
with i iterations.}

\item{...}{Additional arguments (e.g., \code{alternative}) to be passed to
other methods. See \code{stats::cor.test} for further details.}
}
\value{
Description of the Mahalanobis distance.
}
\description{
The Mahalanobis distance (in squared units) measures the distance in
multivariate space taking into account the covariance structure of the data.
Because a few extreme outliers can skew the covariance estimate, the
bootstrapped version is considered as more robust.
}
\examples{
library(correlation)

distance_mahalanobis(iris[, 1:4])
distance_mahalanobis(iris[, 1:4], robust = FALSE)
}
\references{
\itemize{
\item Schwarzkopf, D. S., De Haas, B., & Rees, G. (2012). Better ways to
improve standards in brain-behavior correlation analysis. Frontiers in
human neuroscience, 6, 200.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_text.R
\name{cor_text}
\alias{cor_text}
\title{Correlation text}
\usage{
cor_text(x, show_ci = TRUE, show_statistic = TRUE, show_sig = TRUE, ...)
}
\arguments{
\item{x}{A dataframe with correlation statistics.}

\item{show_ci, show_statistic, show_sig}{Toggle on/off different parts of the text.}

\item{...}{Other arguments to be passed to or from other functions.}
}
\description{
This function returns a formatted character of correlation statistics.
}
\examples{
rez <- cor_test(mtcars, "mpg", "wt")

cor_text(rez)
cor_text(rez, show_statistic = FALSE, show_ci = FALSE, stars = TRUE)

rez <- correlation(mtcars)

cor_text(rez)
}
