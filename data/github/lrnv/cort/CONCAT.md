---
title: 'Empirical and non-parametric copula models with the `cort` R package'
tags:
  - R
  - copula
  - statistics
authors:
  - name: Oskar Laverny^[Corresponding author]
    orcid: 0000-0002-7508-999X
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Institut Camille Jordan, Université Lyon 1, Lyon, France
   index: 1
 - name: SCOR SE, Paris, France
   index: 2
date: 01 Decembre 2020
bibliography: paper.bib
---

# Summary

The R package `cort` implements object-oriented classes and methods to estimate, simulate and visualize certain types of non-parametric copulas.


Copulas are functions that describe dependence structure of a given dataset, or of a given multivariate random event, without describing the univariate events, the marginals. In statistics, it is sometimes useful to separate the marginal distributions (inflation of the money, mortality of the population) from the dependence structure between them, since estimating everything separately is usually easier. Copulas are broadly used in finance, actuarial science, geostatistics, biostatistics, and many other fields, when dealing with dependence.

Copulas are distribution functions on the unit hypercube that have uniform margins (what we call the 'copula constraints'), and hence this package can be classified in 'density estimation software'. Although the estimation of copulas is a widely-treated subject, most performing estimators available in the literature are based on restricted, parametric estimation: vine copulas [@nagler2016evading] and graphical models [@li2018panda] for example are potential solutions but under restrictive assumptions. Classical density estimators such as kernels or wavelets do not satisfy marginal copula constraints. There also exist several tree-structured piece wise constant density estimators, but they do not always lead to proper copulas when applied on pseudo-observations or true copula samples. The new models that are implemented in this package try to solve these issues.

We note that a lot of tools are available in R for copula modeling  through the excellent package `copula` [@cop1;@cop2;@cop3;@cop4]. Most of these tools however focus on parametric estimation. We start to bridge the gap by providing some tools for non-parametric estimation.


# Statement of need 

The Copula recursive tree, or Cort, designed by [@laverny2020dependence] is a flexible, consistent, piece wise linear estimator for a copula [@sklar1959fonctions], leveraging the patchwork copula formalization [@durante2015convergence] and a specific piece wise constant density estimator, the density estimation tree [@ram2011density]. While the patchwork structure imposes the grid, this estimator is data-driven and constructs the grid recursively from the data, minimizing a chosen distance on the copula space. Furthermore, while the addition of the copula constraints makes the available solutions for density estimation unusable, our estimator is only concerned with dependence and guarantees the uniformity of margins. The R package `cort` provides a useful implementation of this model and several potential refinements, allowing for fast computations of Cort trees, and parallel computations of Cort forests.

The main feature implemented in the package is the Cort algorithm, a non-parametric, piece wise constant, copula density estimator. The implementation is recursive and hence quite efficient. The `cort` package is a statistical package that allows to estimate several non-parametric copula models in R. Although the state of the art `copula` package has functions to estimate the empirical copula, we provide a structured set of S4 classes that allows estimation of empirical copulas, checkerboard copulas, Cort copula and bagging of all of these. A specific class exists for bagging Cort models, which implementation runs in parallel, to fasten the computations, using the `future` package [@future]. Most of the underlying machinery and computations are written in `C++`, through the `Rcpp` [@rcpp1;@rcpp2;@rcpp3;@rcpp4] package.

The `cort` package was designed to be used by statisticians who need a non-parametric view of the dependence structure of a given dataset. It features a rich and extensive API to call the statistical fitting procedures, plotting functions and tools to assess the quality of the fits, while complying with the R standards. Examples datasets are included in the package, and the many vignettes give examples of use cases. The package is available on the Comprehensive R Archive Network (CRAN).

# Acknowledgements

We acknowledge contributions from Véronique Maume-deschamps, Didier Rullière and Esterina Masiello, who all gave meaningful insights about performances of the code. We are thankful to the reviewers for the fruitful review of this article and of the package code.

# References
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cort)](https://CRAN.R-project.org/package=cort)
[![CRAN number of downloads](https://cranlogs.r-pkg.org/badges/grand-total/cort)](https://cranlogs.r-pkg.org/badges/grand-total/cort)
[![Codecov test coverage](https://codecov.io/gh/lrnv/cort/branch/master/graph/badge.svg)](https://codecov.io/gh/lrnv/cort?branch=master)
[![DOI](https://zenodo.org/badge/247063359.svg)](https://zenodo.org/badge/latestdoi/247063359)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02653/status.svg)](https://doi.org/10.21105/joss.02653)
[![tic](https://github.com/lrnv/cort/workflows/tic/badge.svg?branch=master)](https://github.com/lrnv/cort/actions)
<!-- badges: end -->

The `cort` package provides S4 classes and methods to fit several copula models: 

* The classic empirical checkerboard copula and the empirical checkerboard copula with known margins, see Cuberos, Masiello and Maume-Deschamps (2019) are proposed. These two models allow to fit copulas in high dimension with a small number of observations, and they are always proper copulas. Some flexibility is added via a possibility to differentiate the checkerboard parameter by dimension. 

* The last model consist of the implementation of the Copula Recursive Tree algorithm, aka. CORT, including the localised dimension reduction, which fits a copula by recursive splitting of the copula domain, see Laverny, Maume-Deschamps, Masiello and Rullière (2020).

* We finally provide an efficient way of mixing copulas, allowing to bag the algorithm into a forest, and a generic way of measuring d-dimensional boxes with a given copula.

## Installation

`cort` is Now on [CRAN](https://CRAN.R-project.org)! You can install the stable version with:

``` r
install.packages("cort")
```

The upstream development version can also be installed with :

``` r
devtools::install_github("lrnv/cort")
```

Note that the installation from github will require the system to have a compiler: 

- Windows: Rtools
- macOS: Xcode CLI
- Linux: r-base-dev (debian)


The vignettes are quite expressive. They give a clear overview of what can be done with this package, how it is coded and why it is useful. Please read them for more details. 

## How to report bugs and get support

To report a bug, feel free to open an issue on the github repository. Support can also be provided through the same chanel if you need it.

## How to contribute

Every contribution is welcome, on the form of pull requests on the github repository. For large modifications, please open an issue for discussions firsts. Concerning the naming convention, the CamelCase functions usually designate classes and constructors of these classes, and all other methods are in snake_case.


## References

Cuberos A, Masiello E, Maume-Deschamps V (2019). “Copulas Checker-Type Approximations: Application to Quantiles Estimation of Sums of Dependent Random Variables.” *Communications in Statistics - Theory and Methods, 1--19. ISSN 0361-0926, 1532-415X.*

Laverny O, Maume-Deschamps V, Masiello E, Rullière D (2020). “Dependence Structure Estimation Using Copula Recursive Trees.” *arXiv preprint arXiv:2005.02912*
# cort (development version)

# cort 0.3.2

## Minor improvements and fixes

* Improved documentation
* Fixed buggy automated tests on CortForest
* Added simulation script to the documentation of the datasets
* Cleared up pkgdown reference page
* Cleaned up some vignettes
* README fixes.

# cort 0.3.1

## Breaking changes

* Infrastructure lightening : the Box class was removed, and some unnecessary generics were also removed.
* cbCopula() now defaults to compute pseudo-observations.
* The forest now defaults to weighting the trees (can be turned off by an option)

## New features

* Core computations have been moved on to Rcpp code for performance. 
* Parallel computations are now possible via furrr in CortForest.
* Solver options and number of bootstrap resamples are now accessible as parameters to the Cort() function.
* Add an option to force the checkerboard grid on trees and inside Cort() and CortForest().
* Adding four example datasets from the paper.

## Minor improvements and fixes

* Removed dependency to magritr.
* Improved documentation
* Cleaned up the code of the Cort algorithm.
* Fixed bug in p-value computations.
* Fixed bug when there is only one leave in the tree.
* Fixed bug in pairs.Cort : dimensions were switched for leaves but not for points.
* Fixed bug in pCopula values for cbCopula objects.


# cort 0.3.0

* First release published on cran !
* Changes to conform to CRAN submissions recomendations
* Added a multiprocess option via furrr

# cort 0.2.0

* Corrections to pass checks on every CI plateform.
* Some corrections to spelling in documentation.

# cort 0.1.0

* Finalised and block the API.

# cort 0.0.0.9001

* The Cort object is now lighter.
* Fastened p-values computations by moving bootstrap to Rcpp.
* Added a vignette with a clayton example.


# cort 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* Making the workspace settings : travis, appveyor, github actions, pkgdowm, ...
* merged the code from former empCop package.
* Implement Cort algorithm
* Implement the CortForest algorithm
* Add many methods to Cort and CortForest class.






## Test environments
* local R installation (Windows 10 x64), R 3.6.0
* Ubuntu 16.04.6 LTS (on travis-ci), R oldrel
* Ubuntu 16.04.6 LTS (on travis-ci), R devel
* Ubuntu 16.04.6 LTS (on travis-ci), R release
* Windows Server x64 (build 17763) (on Appveyor), R release
* Windows Server x64 (build 17763) (on Appveyor), R devel
* macOS-latest (on Github Actions), R devel
* macOS-latest (on Github Actions), R release
* Windows Latest (on Github Actions), R release
* Ubuntu Latest (on Github Actions), R release
* Windows x86_64-w64-mingw32 (on win-builder), R devel
* Windows Server 2008 R2 SP1 (on R-Hub), R devel
* Debian Linux (on R-Hub), R devel
* Ubuntu Linux 16.04 (on R-Hub), R release
* Fedora Linux (on R-Hub), R devel

## R CMD check results
There were no NOTEs, ERRORs or WARNINGs
# cort 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* Trying to make the package run on travis. 
* We will target CRAN with this package. 

# cort

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/lrnv/cort.svg?branch=master)](https://travis-ci.org/lrnv/cort)
[![Codecov test coverage](https://codecov.io/gh/lrnv/cort/branch/master/graph/badge.svg)](https://codecov.io/gh/lrnv/cort?branch=master)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/lrnv/cort?branch=master&svg=true)](https://ci.appveyor.com/project/lrnv/cort)
<!-- badges: end -->

The goal of cort is to ...

## Installation

You can install the released version of cort from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("cort")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(cort)
## basic example code
```

# cort 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* Trying to make the package run on travis. 
* We will target CRAN with this package. 
