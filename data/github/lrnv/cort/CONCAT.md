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
---
title: "1. Empirical Checkerboard Copula"
author: "Oskar Laverny"
date: "`r Sys.Date()`"
bibliography: ../inst/REFERENCES.bib
output:
  rmarkdown::html_vignette:
    fig_caption: yes
vignette: >
  %\VignetteIndexEntry{1. Empirical Checkerboard Copula}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7
)
```
The empirical checkerboard copula is a model for copula defined by [@cuberos2019]. It provides a fast and easy way to model a copula in high-dimensional settings (more columns than rows in the dataset). It only requires one parameter, $m$, the so-called checkerboard parameter. In the first section of this vignette, we will set notations and define the model. The second section will discuss and demonstrate the implementation made in the package.

## Definition of the empirical checkerboard copula

Suppose that we have a dataset with $n$ i.i.d observation of a $d$-dimentional copula (or pseudo-observations). Take the checkerboard parameter $m$ to be an integer dividing $n$. 

Let's consider the ensemble of multi-indexes $\mathcal{I} = \{\mathbf{i} = (i_1,..,i_d) \subset \{1,...,m \}^d\}$ which indexes the boxes : 

$$B_{\mathbf{i}} = \left]\frac{\mathbf{i}-1}{m},\frac{\mathbf{i}}{m}\right]$$

partitioning the space $\mathbb{I}^d = [0,1]^d$.

Let now $\lambda$ be the dimension-unspecific Lebesgue measure on any power of $\mathbb{R}$, that is : 

$$\forall d \in \mathbb{N}, \forall x,y \in \mathbb{R}^p, \lambda(\left[x,y\right]) = \prod\limits_{p=1}^{d} (y_i - x_i)$$

Let furthermore $\mu$ and $\hat{\mu}$ be respectively the true copula measure of the sample at hand and the classical Deheuvels empirical copula, that is : 

- For $n$ i.i.d observation of the copula of dimension $d$, let $\forall i \in \{1,...,d\}, \, R_i^1,...,R_i^d$ be the marginal ranks for the variable $i$. 
- $\forall x \in \mathcal{I}^d, \text{ let } \hat{\mu}([0,x]) = \frac{1}{n} \sum\limits_{k=1}^n \mathbb{1}_{R_1^k\le x_1,...,R_d^k\le x_d}$


We are now ready to define the checkerboard copula, $C$, and the empirical checkerboard copula, $\hat{C}$, by the following : 

$$\forall x \in [0,1]^d, C(x) = \sum\limits_{\mathbf{i}\in\mathcal{I}} {m^d \mu(B_{\mathbf{i}}) \lambda([0,\mathbf{x}]\cap B_{\mathbf{i}})}$$

Where $m^d = \lambda(B_{\mathbf{i}})$. 

This copula is a special form of patchwork copulas (see Durante) and some results are known about it : it is indeed a copula, it converges to the true copula as the mesh (size of boxes) goes to zero, etc..

This package gives a comprehensive implementation of the empirical counterpart of this copula, which has exactly the same expression except that $\mu$, the true copula of the sample, is replace by it's Deheuvel approximation $\hat{\mu}$, that is : 

$$\forall x \in [0,1]^d, \hat{C}(x) = \sum\limits_{\mathbf{i}\in\mathcal{I}} {m^d \hat{\mu}(B_{\mathbf{i}}) \lambda([0,\mathbf{x}]\cap B_{\mathbf{i}})}$$

A known result is that this is a copula if and only if $m$ divides $n$ (see cuberos). In this case, some theoretical assymptotics are avaliables.

The next section discuss the implementation.


## Implementation

In this package, this empirical checkerboard copula is implemented in the `cbCopula` class. This class is a little more general as it allows for a vector $m = (m_1,...,m_d)$ instead of a single value. Each of the $m_i$'s must divide $n$ for this to be a proper copula. With the same train of thoughts as before, we have the following expression for this model : 

$$\forall x \in [0,1]^d, \hat{C}(x) = \sum\limits_{\mathbf{i}\in\mathcal{I}} { \hat{\mu}(B_{\mathbf{i}}) \lambda([0,\mathbf{x}]\cap B_{\mathbf{i}})\prod\limits_{p=1}^{d}m_p}$$

You need to provide a dataset, defining $\hat{\mu}$, to construct the model. For the matter of this vignette, we will use the `LifeCycleSavings` dataset, which has following pairs dependencies plot : 


```{r,fig.cap="Pairs-plot of original peusdo-observation from the data"}
set.seed(1)
data("LifeCycleSavings")
pseudo_data <- (apply(LifeCycleSavings,2,rank,ties.method="max")/(nrow(LifeCycleSavings)+1))

pairs(pseudo_data,lower.panel=NULL)
```


You can see that the variable 2 to 4 have dependencies, while the first and fifth variable seems to be (marginally) independent. The dataset having $n=50$ rows, we will pick a value of $m$ dividing $50$, e.g $m=5$, and use the function `cbCopula` to build our copula model. Since we are already providing the pseudo observations, we will set `pseudo = TRUE`. Providing a single value for the $m$ parameter will set all $m_i$'s equal to that value (the default is the proper checkerboard copula).



```{r}
(cop <- cbCopula(x = pseudo_data,m = 5,pseudo = TRUE))
```

For the moment, only some methods exist for this copula. We can calculate it's values via the `pCopula` method, or simulate from it via the `rCopula` method. the `dCopula` methods gives it's density. Here is an example of simulation from this model : 

```{r, fig.cap = "Pairs-plot of original peusdo-observation from the data (red) with simulated pseudo_observation (black)"}
simu <- rCopula(n = 1000,copula = cop)

pairs(rbind(simu,pseudo_data),
      col=c(rep("black",nrow(simu)),rep("red",nrow(pseudo_data))),
      gap=0,
      lower.panel=NULL,cex=0.5)
```

## About the value of $m$

The value of the checkerboard parameter $m$ condition heavily the copula itself. See the vignette about convex mixtures of m-randomized checkerboards for more details.














---
title: "3. Empirical Checkerboard Copula with known margins"
author: "Oskar Laverny"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_caption: yes
vignette: >
  %\VignetteIndexEntry{3. Empirical Checkerboard Copula with known margins}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7
)
```


The empirical checkerboard copula with known margins is a model for copula that uses priori information on a multidimensional margins to condition a checkerboard construction on the rest of the copula. This vignette describes the model in the first section, and then discuss current implementation in the second section

## Empirical checkerboard copula with known margins

### Preliminary notations


First of all, for a certain dimension of the data at hand $d$ and a certain checkerboard parameter $m$, let's consider the ensemble of multi-indexes $\mathcal{I} = \{\mathbf{i} = (i_1,..,i_d) \subset \{1,...,m \}^d\}$ which indexes the boxes :

$$B_{\mathbf{i}} = \left]\frac{\mathbf{i}-1}{m},\frac{\mathbf{i}}{m}\right]$$

partitioning the space $\mathbb{I}^d = [0,1]^d$. Denote the set of thoose boxes by $\mathcal{B}_\mathcal{I} = \left\{B_{\mathbf{i}}, \mathbf{i} \in \mathcal{I}\right\}$. Furthermore, let's choose $p$ dimensions that would be assigned to a known copula, by setting : $J \subset \{1,...,d\}, |J| = p$ and let's define proper projections for the boxes :

$$B^J_{\mathbf{i}} = \{ \mathbf{x} \in [0,1]^p, x_j \in \left]\frac{i_j-1}{m},\frac{i_j}{m}\right] \forall j \in J \}$$

$$B^{-J}_{\mathbf{i}} = \{\mathbf{x} \in [0,1]^p, x_j \in \left]\frac{i_j-1}{m},\frac{i_j}{m}\right] \forall j \notin J \}$$

such that $B_{\mathbf{i}} = B^J_{\mathbf{i}} \times B^{-J}_{\mathbf{i}}$ for all $\mathbf{i} \in \mathcal{I}$. The tensor product is here understood taking dimensions of $\mathbb{I}^d = [0,1]^d$ in the right order such that $B^J_{\mathbf{i}}$ dimensions end up in place of dimensions indexed by $J$. Think of dimensions as re-ordered such that $J = \{1,...,p\}$.

Let now $\lambda$ be the dimension-unspecified Lebesgue measure on any power of $\mathbb{R}$, that is :

$$\forall d \in \mathbb{N}, \forall \mathbf{x},\mathbf{y} \in \mathbb{R}^p, \lambda(\left[\mathbf{x},\mathbf{y}\right]) = \prod\limits_{p=1}^{d} (y_i - x_i)$$

Let furthermore $\mu^J$ be a copula measure of dimension $p$, corresponding to the known multivariate margin associated to marginals in $J$. Let also $\mu$ and $\hat{\mu}$ be dimension-unspecific  version of the true copula measure of the sample at hand and (respectively) the classical Deheuvels empirical copula, that is :

- For $n$ i.i.d observation of the copula of dimension $d$, let $\forall i \in \{1,...,d\}, \, R_i^1,...,R_i^d$ be the marginal ranks for the variable $i$.
- $\forall \mathbf{x} \in \mathcal{I}^d, \text{ let } \hat{\mu}([0,x]) = \frac{1}{n} \sum\limits_{k=1}^n \mathbb{1}_{R_1^k\le x_1,...,R_d^k\le x_d}$



We are now ready to define the empirical checkerboard copula with known margins.

### Definition, estimation and simulation procedures

The empirical copula with known margins is the copula that correspond to the following simulation procedure.

- Simulate a sample from the known sub-copula $\mu^J$, of dimension $p$, through any avaliable method (depends on the known copula model). Let $B_{\mathbf{i}}^J$ be the (projected) box containing this sample.
- Sample one box $B_{\mathbf{i}}$ among all boxes with projection $B_{\mathbf{i}}^J$ with probability weights :
  -  $\frac{\hat{\mu}(B_{\mathbf{i}})}{\hat{\mu}(B_{\mathbf{i}}^J)}$ if the projected box $B_{\mathbf{i}}^J$ contains one or more (empirical) data point, that is $\hat{\mu}(B_{\mathbf{i}}^J) \neq 0$
	- $\frac{\lambda(B_{\mathbf{i}})}{\lambda(B_{\mathbf{i}}^J)}$ otherwise.
- Simulate uniformly from $B_{\mathbf{i}}^{-J}$

This algorithm simulates first the known part of the model (dimensions in $J$), and then, conditionally, the checkerboard part, ensuring that the known copula is respected. The downside of this behavior is that the checkerboard part may have points outside standard checkerboard boxes, making this part of the copula less sparse than a true checkerboard. But is does become sparser as soon as the data fits the known margins. On the other hand, this algorithm allows for a lot of flexibility, mainly in the following points :


- The "grid" given by $\mathcal{B}_\mathcal{I}$ can be taken more arbitrarily than $m^d$ boxes of same volume, as soon as it's a partition of $\mathbb{I}^d$.
- The known copula is *not* restricted at all and can be chose among all $p$-dimensionnal copulas.
- The estimation of the checkerboard part can be turn into a more flexible *patchwork* construction, by changing the independence copula for an other one inside the boxes. See Durante2013,Durante2015,Durante2015a

We are now going to define properly the measure associated to this simulation procedure, a.k.a the empirical checkerboard copula with known margins. Let $\nu$ be this measure and let $\mathbf{U}$ be a random vector drawn from it. Then $\forall \mathbf{x} \in \mathbb{I}^d$, following the above procedure, we have :

$$\nu([0,\mathbf{x}]) = \sum\limits_{{\mathbf{i}} \in \mathcal{I}} \mathbb{P}(\mathbf{U}^{-J} \in B_{\mathbf{i}}^{-J} \cap [0,\mathbf{x}^{-J}] | \mathbf{U}^{J} \in B_{\mathbf{i}}^{J} \cap [0,\mathbf{x}^{J}]) \mathbb{P}(\mathbf{U}^{J} \in B_{\mathbf{i}}^{J} \cap [0,\mathbf{x}^{J}])$$
	

While the unconditional term is easy to handle since it's the measure associated with the known copula, $\mathbb{P}(\mathbf{U}^{J} \in B_{\mathbf{i}}^{J} \cap [0,\mathbf{x}^{J}]) = \mu^J(B_{\mathbf{i}}^{J} \cap [0,\mathbf{x}^{J}])$, the conditional term can be treated according to the algorithm : it will be $\frac{\lambda(B_{\mathbf{i}}^{-J} \cap [0,\mathbf{x}^{-J}])}{\lambda(B_{\mathbf{i}}^{-J})}$ inside a box chosen with probability conditional on $\hat{\mu}(B_{\mathbf{i}}^J) \neq 0$. We finally get the following definition :



**The empirical checkerboard copula with parameter $m$ and with set of known margins $J$ following the measure $\mu^J$ is the copula corresponding to the measure $\nu$ given by : **

$$
\nu([0,x]) = \sum\limits_{{\mathbf{i}} \in \mathcal{I}} \mu^J(B_{\mathbf{i}}^{J} \cap [0,\mathbf{x}^{J}]) \frac{\lambda(B_{\mathbf{i}}^{-J} \cap [0,\mathbf{x}^{-J}])}{\lambda(B_{\mathbf{i}}^{-J})} \left[ \frac{\hat{\mu}(B_{\mathbf{i}})}{\hat{\mu}(B_{\mathbf{i}}^J)} \mathbb{1}_{\hat{\mu}(B_{\mathbf{i}}^J) \neq 0} + \frac{\lambda(B_{\mathbf{i}})}{\lambda(B_{\mathbf{i}}^J)} \mathbb{1}_{\hat{\mu}(B_{\mathbf{i}}^J) = 0}   \right]
$$



The next section will discuss the current implementation of this copula;

## Current implementation

The package implements the empirical checkerboard copula with known margins through the `cbkmCopula` class. The constructor of the class takes several arguments :

- `x`, the pseudo_data.
- `m=nrow(x)`, repesenting the checkerboard parameter
- `pseudo=FALSE`,  is the pseudo_data already given in a pseudo_observation form ?
- `margins_numbers=NULL`, the margins index that are associated to the known copula, formerly noted $J$
- `known_cop=NULL`, the known copula to be applied to thoose margins : a copula object of right dimension.

for example, let's take the LifeCycleSavings data : 

```{r,fig.cap="Pairs-plot of original peusdo-observations"}
set.seed(1)
data("LifeCycleSavings")
dataset <- (apply(LifeCycleSavings,2,rank,ties.method="max")/(nrow(LifeCycleSavings)+1))
pairs(dataset,col="2",lower.panel=NULL)

```


let's now estimate a checkerboard copula on margins 2 and 3 with a precise $m=25$, and consider it to be known.
```{r}
  known_margins <- c(2,3)
  known_copula <- cbCopula(x = dataset[,known_margins],m = 25,pseudo = TRUE)
```

Then we can construct the ECBC with this known margin :

```{r}
  (cop <- cbkmCopula(x = dataset,m = 5,pseudo = TRUE,margins_numbers = known_margins,known_cop = known_copula))
```

We can then simulate from it :

```{r, fig.cap="Pairs-plot of the original data (red) and simulated data from a good model (black)"}
  simu <- rCopula(1000,cop)
  pairs(rbind(simu,dataset),col=c(rep("black",nrow(simu)),rep("red",nrow(dataset))),gap=0,lower.panel = NULL,cex=0.5)
```



You can see that the known margins were respected, which is the whole point of this model. Let now see an example with a clearly wrongly-specified copula for the known margins :
```{r, fig.cap="Pairs-plot of the original data (red) and simulated data from a wrong model (black)"}
  bad_known_copula <- cbCopula(x = dataset[,known_margins],m = 2,pseudo = TRUE)
  cop <- cbkmCopula(x = dataset,m = 5,pseudo = TRUE,margins_numbers = known_margins,known_cop = bad_known_copula)

  simu <- rCopula(1000,cop)
  pairs(rbind(simu,dataset),col=c(rep("black",nrow(simu)),rep("red",nrow(dataset))),gap=0,lower.panel = NULL,cex=0.5)
```

This example shows the sensitivity of the estimation on other parts of the model to the good behavior of the prior estimate.
The conditioning did suffer for the dependences inside the known multidimensional margin but also for dependencies involving one of the variables from thoose known margins. But the checkerboard construction for the other part of the copula was not harmed.

What now if the 2 parts are clearly independent ?

```{r, fig.cap="Pairs-plot of the original data with independance"}

  true_copula1 <- known_copula
  true_copula2 <- bad_known_copula

  dataset <- cbind(rCopula(1000,true_copula1),rCopula(1000,true_copula2))
  colnames(dataset) <- c("u","v","w","x")
  pairs(dataset,lower.panel=NULL,cex=0.5)

```

```{r, fig.cap="Pairs-plot of the original data (red) and simulated data (black) -- Independance case"}
  cop <- cbkmCopula(x = dataset,m = 5,pseudo = TRUE,margins_numbers = c(1,2),known_cop = true_copula1)
  simu <- rCopula(500,cop)
  pairs(rbind(simu,dataset),col=c(rep("black",nrow(simu)),rep("red",nrow(dataset))),gap=0,lower.panel = NULL,cex=0.5)
```

This fit is quite good.

---
title: "4. Convex mixture of m-randomized checkerboards"
author: "Oskar Laverny"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_caption: yes
vignette: >
  %\VignetteIndexEntry{4. Convex mixture of m-randomized checkerboards}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r setup, include = FALSE}
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7
)
```

# Introduction

In this vignette, we will demonstrate how the checkerboard parameters $m$ can be used to produce an efficient meta-model. We will construct $m$-randomized checkerboard copulas, and based on some fit statistics we will then combine them in a convex mixture. Our model here will use all possibles $m_i$'s values for checkerboard copulas and aggregate all those values by optimizing a quadratic loss under resampling condition. As usual, we will work on the LifeCycleSavings dataset : 

```{r "getting data and parameters"}
set.seed(1)
df <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
d = ncol(df)
n = nrow(df)
nb_replicates = 20 # number of replication of the resampling.
nb_fold = 5 # "k" value for the k-fold resampling method.
nb_cop = 50 # the number of m-randomized checkerboard copulas we will test.
pairs(df,lower.panel = NULL)
```

# Fitting function

The following functions proposes a fit : 

```{r "fitting_functions"}

make_k_fold_samples <- function(data,k = nb_fold,n_repeat = 1){
  
  sapply(1:n_repeat,function(n){
      # shuffle data : 
      data <- data[sample(nrow(data)),]
    
      # create k folds : 
      folds <- cut(seq(1,nrow(data)),breaks=k,labels=FALSE)
      
      sapply(1:k,function(i){
        test_index <- which(folds==i,arr.ind=TRUE)
        return(list(train = data[-test_index,], 
                    test = data[test_index,]))
      },simplify=FALSE)
  })
}

build_random_m <- function(how_much=1,dim = d,nrow_data){
  t(sapply(1:how_much,function(i){
    m_pos = (2:nrow_data)[nrow_data%%(2:nrow_data)==0]
    sample(m_pos,d,replace=TRUE)
  }))
}

build_all_checkerboard <- function(sampled_data,m){
  lapply(sampled_data,function(d){
    apply(m,1,function(m_){
      cbCopula(x = d$train,m = m_,pseudo=TRUE)
    })
  })
}

samples <- make_k_fold_samples(df,k=nb_fold,n_repeat=nb_replicates)
rand_m <- build_random_m(nb_cop,d,nrow(samples[[1]]$train))
cops <- build_all_checkerboard(samples,rand_m)

```



Let's first calculate the empirical copula values : 

```{r}
pEmpCop <- function(points,data=df){
  sapply(1:nrow(points),function(i){
    sum(colSums(t(data) <= points[i,]) == d)
  }) / nrow(data)
}

```


Now, we also need to calculate the errors that our copulas made compared to the empirical copula (our benchmark).
```{r}
error <- function(cop,i,j){
  test <- samples[[i]]$test
  return(sum((pCopula(test,cop) - pEmpCop(test))^2))
  
}

errors <- sapply(1:(nb_replicates*nb_fold),function(i){
  sapply(1:nb_cop,function(j){
   error(cops[[i]][[j]],i,j)
  })
})

rmse_by_model <- sqrt(rowMeans(errors))
plot(rmse_by_model)

```


Each point on the graph correspond to the rmse of a model. We recall the values of $m$ used by those models : 

```{r}
rand_m
```




```{r}
convex_combination <- ConvexCombCopula(unlist(cops,recursive=FALSE),alpha = rep(1/rmse_by_model,nb_replicates*nb_fold))
simu = rCopula(1000,convex_combination)
pairs(simu,lower.panel = NULL,cex=0.5)
```

Which is quite good compared to the dataset we started with. This process is really fast and useful, and can of course be used for high-dimensional datasets.

The parameters `nb_fold`, `nb_repeats`, `nb_cop` could be further optimized. Furthermore, if some multivariate margins are known, the same thing could be done using `cbkmCopula()` instead of `cbCopula()` models.


---
title: "2. The Copula Recursive Tree"
bibliography: ../inst/REFERENCES.bib
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_caption: yes
vignette: >
  %\VignetteIndexEntry{2. The Copula Recursive Tree}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction

```{r setup, include = FALSE}
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
```

In this vignette, we will use the Cort algorithm to fit a copula on a given simulated dataset based on clayton copula simulations. We will show how the algorithm can be used to produce a bona-fide copula, and describe some of the parameters. 

# Dataset

First, let's create and plot the dataset we will work with. For that, we'll use the gamma frailty model for the Clayton copula (but it'll work for any other completely monotonous archimedean generator), as it is done in the `copula` package, see [there](https://github.com/cran/copula/blob/6e0a5f0332776c24b2f53a81b11a8c66cc3f996d/R/claytonCopula.R#L82). The following code is directly taken from the previous link, from the `copula` package : 

```{r clayton-frailty-model}
psi <- function(t,alpha) (1 + sign(alpha)*t) ^ (-1/alpha) # generator
rClayton <- function(n,dim,alpha){
  val <- matrix(runif(n * dim), nrow = n)
  gam <- rgamma(n, shape = 1/alpha, rate = 1)
  gam <- matrix(gam, nrow = n, ncol = dim)
  psi(- log(val) / gam,alpha)
}
```


The following code simulates a dataset and then visualise it : 

```{r dataset} 
set.seed(12)
n = 200 # taken small to reduce runtime of the vignette.
d = 4
n_trees = 5 # taken small to reduce runtime of the vignette.
number_max_dim_forest = 2 # taken small to reduce runtime of the vignette.

data <- matrix(nrow=n,ncol=d)
data[,c(1,4,3)] = rClayton(n=n,dim=d-1,alpha=7)
data[,2] = runif(n)
data[,3] <- 1 - data[,3]

pairs(data,cex=0.6)
  
```

We can clearly see that the second marginal is independent form the rest. In the following we will use this package to fit this dependence structure. 


# Fitting the Cort copula

Now that we have a dataset, we can run the Cort algorithm on it. In the implementation proposed here, this is done via the `cort::Cort()` function, passing first the dataset, and then various parameters. See `?Cort` for a detailed list of parameters. Note that the verbosity level is quite progressive: We will here put it on 4 to see the splitting decisions that the algorithm is making.

```{r run_cort} 
set.seed(12)
(model = Cort(data,verbose_lvl = 1))
```

Looking at the top of the output, we see that the first thing the algorithm did was removing the second dimension due to the independence test. Now that the copula is fitted, we have access to numerous of it's methods. Two plotting functions are exported with this model, the `pairs` function is implemented at a very low level in the class hierarchy and hence is working with almost all copulas of this package, but the `plot` function is only implemented for Cort.  


```{r,fig.cap="Pairs-plot of original data (in black, bottom-left corner) versus a simulation from the model (in red, top-right corner)"}
pairs(model)
```

```{r,fig.cap="Gray boxes representing 2-d projections of the fitted density. In red, the imputed data points."}
plot(model)
```


We see that there are some noise with point were there should not be. A bagged version of the model is accessible via the `CortForest` class, and might be able to correct these problems.






% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{clayton_data}
\alias{clayton_data}
\title{Dataset clayton_data}
\format{
A matrix with 200 rows and 4 columns

The example section below gives the code to re-generate this data if needed.
}
\usage{
clayton_data
}
\description{
This dataset is a simulation of 200 points from a 3-dimensional clayton copula with \eqn{\theta = 7},
hence highly dependent, for the first, third and fourth marginals. The second marginal is added
as independent uniform draws. Lastly, the third marginal is flipped,
inducing a negative dependence structure.
}
\details{
This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
}
\examples{
psi <- function(t,alpha) (1 + sign(alpha)*t) ^ (-1/alpha) # generator
rClayton <- function(n,dim,alpha){
  val <- matrix(runif(n * dim), nrow = n)
  gam <- rgamma(n, shape = 1/alpha, rate = 1)
  gam <- matrix(gam, nrow = n, ncol = dim)
  psi(- log(val) / gam,alpha)
}
set.seed(12,kind = "Mersenne-Twister",normal.kind = "Inversion")
clayton_data <- matrix(nrow=200,ncol=4)
clayton_data[,c(1,4,3)] = rClayton(n=200,dim=3,alpha=7)
clayton_data[,2] = runif(200)
clayton_data[,3] <- 1 - clayton_data[,3]

}
\references{
\insertRef{laverny2020}{cort}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cort.R
\name{Cort-Class}
\alias{Cort-Class}
\alias{Cort}
\title{Cort copulas}
\usage{
Cort(
  x,
  p_value_for_dim_red = 0.75,
  min_node_size = 1,
  pseudo_data = FALSE,
  number_max_dim = NULL,
  verbose_lvl = 1,
  slsqp_options = NULL,
  osqp_options = NULL,
  N = 999,
  force_grid = FALSE
)
}
\arguments{
\item{x}{The data, must be provided as a matrix with each row as an observation.}

\item{p_value_for_dim_red}{a p_value for the localized dimension reduction test}

\item{min_node_size}{The minimum number of observation available in a leaf to initialize a split.}

\item{pseudo_data}{set to True if you are already providing data on the copula space.}

\item{number_max_dim}{The maximum number of dimension a split occurs in. Defaults to be all of the dimensions.}

\item{verbose_lvl}{numeric. set the verbosity. 0 for no output and bigger you set it the most output you get.}

\item{slsqp_options}{options for nloptr::slsqp to find breakpoints : you can change defaults.}

\item{osqp_options}{options for the weights optimization. You can pass a call to osqp::osqpSettings, or NULL for defaults.}

\item{N}{The number of bootstrap samples for p_values computations.}

\item{force_grid}{Set to TRUE to force breakpoints to be on the n-checkerboard grid.}
}
\value{
An instance of the \code{Cort} S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, constraint influence, etc.
Beside returning some inputted parameters, notable slots are :

\itemize{
\item{\code{data} }{Your original data}
\item{\code{dim} }{The dimension of problem, number of columns of your dataset}
\item{\code{f} }{The empirical frequency in the leaves}
\item{\code{p} }{The fitted probabilities of each leaf}
\item{\code{a} }{Minimum points of leaves}
\item{\code{b} }{Maximum points of leaves}
\item{\code{vols} }{Volume of the leaves}
}

More details about these slots can be found in the reference.
}
\description{
Cort class
}
\details{
This class implements the CORT algorithm to a fit a multivariate copula using piece constant density. Given a dataset \code{x}, the function will produce an estimator for the copula of this dataset
that is tree-shaped, by recursive partitioning of the unit hypercube. the \code{min_node_size} parameter controls the stopping conditions for the splitting procedure. Once the space is splitted,
we ran a quadratic solver, which options can be tweaked via the \code{osqp_options} parameter, to ensure that the weights respect the copula conditions.

Once the model is fitted, it can be used through the classical (r/d/p/v)Copula functions to compute, respectively, random number generations, the density, the cdf and the volume function of the copula.

See O. Laverny, E. Masiello, V. Maume-Deschamps and D. Rullière (2020) for the details of this density estimation procedure, and \code{vignettes(package='cort')} for examples of usecases.
}
\examples{
(Cort(LifeCycleSavings[,1:3]))
}
\references{
\insertRef{laverny2020}{cort}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{impossible_data}
\alias{impossible_data}
\title{Dataset impossible_data}
\format{
A matrix with 200 rows and 2 columns

The example section below gives the code to re-generate this data if needed.
}
\usage{
impossible_data
}
\description{
We simulate from a density inside the piecewise linear copula class, by applying the function:
\deqn{h(u) = (u_1,          \frac{u_2}{2} + \frac{1}{2}I_{u_1 \notin (\frac{1}{3}, \frac{2}{3})})}
to a 200x2 uniform sample, and taking ranks.
}
\details{
This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
}
\examples{
set.seed(seed = 12, kind = "Mersenne-Twister", normal.kind = "Inversion")
x = matrix(runif(400),200,2)
x = t(apply(x, 1,function(u){
  if(u[1]< 1/3){
    u[2] = 1/2 + u[2]/2
  } else{ if(u[1]<2/3){
    u[2] = u[2]/2
  } else {
    u[2] = 1/2 + u[2]/2
  }}
  return(u)
}))
impossible_data = apply(x,2,function(x){return(rank(x,ties.method = "max"))})/(201)

}
\references{
\insertRef{laverny2020}{cort}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConvexCombCopula.R
\name{ConvexCombCopula-Class}
\alias{ConvexCombCopula-Class}
\alias{ConvexCombCopula}
\title{Convex Combination of copulas.}
\usage{
ConvexCombCopula(copulas, alpha = rep(1, length(copulas)))
}
\arguments{
\item{copulas}{a list of copulas of same dimension}

\item{alpha}{a vector of (positive) weights}
}
\value{
An instance of the \code{ConvexCombCopula} S4 class. The object represent the copula that results from a convex combinaison of other copulas, and can be used through several methods to query classical (r/d/p/v)Copula methods, etc.
}
\description{
ConvexCombCopula class
}
\details{
The ConvexCombcopula class is used to build convex combinations of copulas,
with given positives weights. The rCopula and pCopula functions works for
those copulas, assuming they work for the given copulas that we combined
in a convex way.

See the corresponding vignette for more details about the implementation.
}
\examples{
dataset <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
copulas <- list(
  cbCopula(dataset[,2:3],m=10),
  cbCopula(dataset[,2:3],m=5)
)
alpha <- c(1,4)
(cop <- ConvexCombCopula(copulas,alpha))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{vCopula}
\alias{vCopula}
\alias{vCopula,matrix,matrix-method}
\alias{vCopula,matrix,matrix,Copula}
\title{Copula volume on hyper-boxes}
\usage{
vCopula(u, v, copula, ...)

\S4method{vCopula}{matrix,matrix}(u, v, copula)
}
\arguments{
\item{u}{numeric matrix : minimum point of the hyper-rectangles, one row per observation.}

\item{v}{numeric matrix : maximum point of the hyper-rectangle, one row per observation.}

\item{copula}{the copula that we compute the measure on the box (u,v)}

\item{...}{other parameter to be passed to methods for this generic.}
}
\value{
the measure of the copula.
}
\description{
u must be piecewise smaller than v, otherwise the function will return an error.
}
\details{
A method is currently implemented for the main virtual class 'Copula', but it assumes
that a pCopula method is avaliable for the given copula. This method could be used with Copulas that are not from this package, assuming that pCopula(u,cop) works.

This function computes the measure of the copula according to the algorithm proposed by Cherubini U, Romagnoli S (2009-oct).
}
\examples{
cop <- cbCopula(LifeCycleSavings,m = 5)
vCopula(rep(0,5),rep(1,5),cop) == 1
vCopula(rep(0,5),rep(0.5,5),cop)

}
\references{
\insertRef{cherubini2009}{cort}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cbkmCopula.R
\name{cbkmCopula-Class}
\alias{cbkmCopula-Class}
\alias{cbkmCopula}
\title{Checkerboards with known margins}
\usage{
cbkmCopula(
  x,
  m = rep(nrow(x), ncol(x)),
  pseudo = FALSE,
  margins_numbers = NULL,
  known_cop = NULL
)
}
\arguments{
\item{x}{the data to be used}

\item{m}{checkerboard parameter}

\item{pseudo}{Boolean, defaults to \code{FALSE}. Set to \code{TRUE} if you are already providing pseudo-data into the \code{x} argument.}

\item{margins_numbers}{numeric integers which determines the margins for the known copula.}

\item{known_cop}{Copula a copula object representing the known copula for the selected margins.}
}
\value{
An instance of the \code{cbkmCopula} S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, etc.
}
\description{
cbkmCopula contructor
}
\details{
Given some empirical data, and given some known copula estimation on a sub-vector of this data,
the checkerboard with known margins construction consist in
a conditional pattern where a checkerboard copula is fitted (similar the the \code{cbCopula} algorithm), but conditionally on some known margins.

See the corresponding vignette for more details.
}
\examples{
dataset <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
known_copula <- cbCopula(dataset[,2:3],m=10)
(cop <- cbkmCopula(x = dataset,
                  m = 5,
                  pseudo = TRUE,
                  margins_numbers = c(2,3),
                  known_cop = known_copula))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R, R/CortForest.R
\name{quad_norm}
\alias{quad_norm}
\alias{quad_norm,Cort-method}
\alias{quad_norm,CortForest-method}
\title{Quadratic norm of the model (if it has one)}
\usage{
quad_norm(object)

\S4method{quad_norm}{Cort}(object)

\S4method{quad_norm}{CortForest}(object)
}
\arguments{
\item{object}{the copula object}
}
\value{
the Integrated square error quad_norm of the model
}
\description{
Currently only implemented for Cort models.
Compute the L2 norm of the model
}
\section{Functions}{
\itemize{
\item \code{quad_norm,Cort-method}: Method for the class Cort

\item \code{quad_norm,CortForest-method}: Method for the class CortForest
}}

\examples{
cop <- Cort(cort::impossible_data)
quad_norm(cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{kendall_func}
\alias{kendall_func}
\alias{kendall_func,Cort-method}
\title{Kendall function of a copula (if it has one)}
\usage{
kendall_func(object, t, ...)

\S4method{kendall_func}{Cort}(object, t, M = 1000)
}
\arguments{
\item{object}{: the tree}

\item{t}{: the value where to compute the kendall function, may be a vector of evaluation values;}

\item{...}{other parameters passed to methods}

\item{M}{the number of simulations}
}
\value{
the quadratic product between the trees
}
\description{
Currently only implemented for Cort models.
Compute the Kendall cdf from the model in a point t
}
\section{Functions}{
\itemize{
\item \code{kendall_func,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(LifeCycleSavings[,1:3])
kendall_func(cop,0.5)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{funcdep_data}
\alias{funcdep_data}
\title{Dataset funcdep_data}
\format{
A matrix with 500 rows and 3 columns

The example section below gives the code to re-generate this data if needed.
}
\usage{
funcdep_data
}
\description{
This dependence structure is constructed by applying the function :
\deqn{h(u_1,u_2,u_3) = (u_{1},\sin(2\pi u_{1})-\frac{u_{2}}{\pi},(1+\frac{u_{3}}{\pi^{2}})(\frac{u_{3}}{2} I_{\frac{1}{4}\ge u_1}-\sin(\pi^{x_{1}}) I_{\frac{1}{4} < u_{1}}))}
to uniformly drawn 3-dimensional random vectors. The dataset is the ranks of \eqn{h(u)}.
}
\details{
This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
}
\examples{
set.seed(seed = 12,kind = "Mersenne-Twister",normal.kind = "Inversion")
x = matrix(runif(1500),500,3)
x[,2] = sin(2*pi*x[,1])-x[,2]/pi
x[,3] = (x[,3]*(x[,1]<1/4)/2 - sin(pi**(x[,1]))*(x[,1]>1/4))*(1+x[,3]/(pi^2))
funcdep_data = apply(x,2,function(x){return(rank(x,ties.method = "max"))})/(501)

}
\references{
\insertRef{laverny2020}{cort}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{biv_rho}
\alias{biv_rho}
\alias{biv_rho,Cort-method}
\title{Spearman's rho matrix of a copula}
\usage{
biv_rho(copula)

\S4method{biv_rho}{Cort}(copula)
}
\arguments{
\item{copula}{the copula object}
}
\value{
the density of the copula on each observation
}
\description{
Computes the bivariate Spearmann's rho matrix for a copula.
}
\section{Functions}{
\itemize{
\item \code{biv_rho,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(LifeCycleSavings[,1:3])
biv_rho(cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{quad_prod_with_data}
\alias{quad_prod_with_data}
\alias{quad_prod_with_data,Cort-method}
\title{Quadratic product with data of the model (if it has one)}
\usage{
quad_prod_with_data(object)

\S4method{quad_prod_with_data}{Cort}(object)
}
\arguments{
\item{object}{the copula object}
}
\value{
the quad_prod_with_data of the model
}
\description{
Currently only implemented for Cort models.
Compute the quadratic product with the empirical density from the data
}
\section{Functions}{
\itemize{
\item \code{quad_prod_with_data,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(LifeCycleSavings[,1:3])
quad_prod_with_data(cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/ConvexCombCopula.R, R/Cort.R,
%   R/CortForest.R, R/cbCopula.R, R/cbkmCopula.R
\name{pCopula}
\alias{pCopula}
\alias{pCopula,matrix,ConvexCombCopula-method}
\alias{pCopula,matrix,Cort-method}
\alias{pCopula,matrix,CortForest-method}
\alias{pCopula,matrix,cbCopula-method}
\alias{pCopula,matrix,cbkmCopula-method}
\title{Copula cdf}
\usage{
pCopula(u, copula, ...)

\S4method{pCopula}{matrix,ConvexCombCopula}(u, copula)

\S4method{pCopula}{matrix,Cort}(u, copula)

\S4method{pCopula}{matrix,CortForest}(u, copula)

\S4method{pCopula}{matrix,cbCopula}(u, copula)

\S4method{pCopula}{matrix,cbkmCopula}(u, copula)
}
\arguments{
\item{u}{numeric matrix : one row per observation}

\item{copula}{the copula object}

\item{...}{other parameter to be passed to methods for this generic.}
}
\value{
The value of the copula on each observation
}
\description{
This function returns the value of the copula itself on given points.
}
\section{Functions}{
\itemize{
\item \code{pCopula,matrix,ConvexCombCopula-method}: Method for the cbCopula

\item \code{pCopula,matrix,Cort-method}: Method for the class Cort

\item \code{pCopula,matrix,CortForest-method}: Method for the class CortForest

\item \code{pCopula,matrix,cbCopula-method}: Method for the cbCopula

\item \code{pCopula,matrix,cbkmCopula-method}: Method for the cbCopula
}}

\examples{
cop <- cbCopula(cort::recoveryourself_data,m = 5)
pCopula(rep(0,2),cop) == 0
pCopula(rep(0.5,2),cop)
pCopula(rep(1,2),cop) == 1

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{biv_tau}
\alias{biv_tau}
\alias{biv_tau,Cort-method}
\title{Kendall's tau matrix of a copula}
\usage{
biv_tau(copula)

\S4method{biv_tau}{Cort}(copula)
}
\arguments{
\item{copula}{the copula object}
}
\value{
the density of the copula on each observation
}
\description{
Computes the bivariate Kendall's tau matrix for a copula.
}
\section{Functions}{
\itemize{
\item \code{biv_tau,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(cort::funcdep_data[1:10,1:3])
biv_tau(cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{recoveryourself_data}
\alias{recoveryourself_data}
\title{Dataset recoveryourself_data}
\format{
A matrix with 500 rows and 2 columns

The example section below gives the code to re-generate this data if needed.
}
\usage{
recoveryourself_data
}
\description{
This dataset is a simple test: we simulate random samples from a density inside the piecewise copula class,
and test whether or not the estimator can recover it. For that, we will use a 2-dimensional sample with 500
observations, uniform on the unit hypercube, and apply the following function:
\deqn{h(u) = (u_1, \frac{u_2 + I_{u_1 \le \frac{1}{4}} + 2I_{u_1 \le \frac{1}{2}} + I_{\frac{3}{4} \le u_1}}{4})}
}
\details{
This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
}
\examples{
set.seed(seed = 12, kind = "Mersenne-Twister", normal.kind = "Inversion")
x = matrix(runif(1000),500,2)
recoveryourself_data = t(apply(x, 1,function(u){
  if(u[1]< 1/4){
    u[2] = 3/4 + u[2]/4
  } else{ if(u[1]<1/2){
    u[2] = 1/2 + u[2]/4
  } else { if(u[1]<3/4){
    u[2] = u[2]/4
  } else {
    u[2] = 1/4 + u[2]/4
  }}}
  return(u)
}))

}
\references{
\insertRef{laverny2020}{cort}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/ConvexCombCopula.R, R/Cort.R,
%   R/CortForest.R, R/cbCopula.R, R/cbkmCopula.R
\name{rCopula}
\alias{rCopula}
\alias{rCopula,numeric,ConvexCombCopula-method}
\alias{rCopula,numeric,Cort-method}
\alias{rCopula,numeric,CortForest-method}
\alias{rCopula,numeric,cbCopula-method}
\alias{rCopula,numeric,cbkmCopula-method}
\title{Copula random generation}
\usage{
rCopula(n, copula, ...)

\S4method{rCopula}{numeric,ConvexCombCopula}(n, copula)

\S4method{rCopula}{numeric,Cort}(n, copula)

\S4method{rCopula}{numeric,CortForest}(n, copula)

\S4method{rCopula}{numeric,cbCopula}(n, copula)

\S4method{rCopula}{numeric,cbkmCopula}(n, copula)
}
\arguments{
\item{n}{the number of simulations}

\item{copula}{the copula object}

\item{...}{other parameter to be passed to methods for this generic.}
}
\value{
A matrix with \code{n} rows, each representing a random vector generated from the provided copula.
}
\description{
Random number generation following the given copula. This function performs the simulation of random vectors following the copula.
}
\section{Functions}{
\itemize{
\item \code{rCopula,numeric,ConvexCombCopula-method}: Method for the cbCopula

\item \code{rCopula,numeric,Cort-method}: Method for the class Cort

\item \code{rCopula,numeric,CortForest-method}: Method for the class CortForest

\item \code{rCopula,numeric,cbCopula-method}: Method for the cbCopula

\item \code{rCopula,numeric,cbkmCopula-method}: Method for the cbCopula
}}

\examples{
cop <- cbCopula(cort::clayton_data,m = 5)
xx <- rCopula(1000,cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CortForest.R
\name{CortForest-Class}
\alias{CortForest-Class}
\alias{CortForest}
\title{Bagged Cort copulas}
\usage{
CortForest(
  x,
  p_value_for_dim_red = 0.75,
  n_trees = 10,
  compte_loo_weights = FALSE,
  min_node_size = 1,
  pseudo_data = FALSE,
  number_max_dim = NULL,
  verbose_lvl = 2,
  force_grid = FALSE,
  oob_weighting = TRUE
)
}
\arguments{
\item{x}{The data, must be provided as a matrix with each row as an observation.}

\item{p_value_for_dim_red}{a p_value for the localised dimension reduction test}

\item{n_trees}{Number of trees}

\item{compte_loo_weights}{Defaults to FALSE. Allows to use an automatic re-weighting of the trees in the forest, based on leave-one-out considerations.}

\item{min_node_size}{The minimum number of observation avaliable in a leaf to initialise a split.}

\item{pseudo_data}{set to True if you are already providing data on the copula space.}

\item{number_max_dim}{The maximum number of dimension a split occurs in. Defaults to be all of the dimensions.}

\item{verbose_lvl}{verbosity level : can be 0 (default) or an integer. bigger the integer bigger the output level.}

\item{force_grid}{boolean (default: FALSE). set to TRUE to force breakpoint to be on the n-checkerboard grid in every tree.}

\item{oob_weighting}{boolean (default : TRUE) option to weight the trees with an oob criterion (otherwise they are equally weighted)}
}
\value{
An instance of the \code{CortForest} S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, constraint influence, etc.
Beside returning some inputted parameters, notable slots are :

\itemize{
\item{\code{trees} }{A list of Cort objects representing each fitted tree in the forest.}
\item{\code{weights} }{The weigths of each tree.}
\item{\code{indexes} }{The indexes of data points that were selected for fitting the trees}
\item{\code{pmf} }{The density of each tree on data points}
\item{\code{norm_matrix} }{The matrix of scalar product between trees}
\item{\code{oob_pmf} }{The density of each tree on data points it did not see during fitting}
\item{\code{oob_kl} }{The out-of-bag Kullback-Leibler divergence of each tree }
\item{\code{oob_ise} }{The out-of-bag Integrated Square Error of each tree}
}

More details about these slots can be found in the reference.
}
\description{
CortForest class
}
\details{
This class implements the bagging of CORT models, with an out-of-bag error minimisation in the weights.

See O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020) for the details of this density estimation procedure, and \code{vignettes(package='cort')} for examples of usecases.
}
\examples{
(CortForest(LifeCycleSavings[,1:3],number_max_dim=2,n_trees=2))
}
\references{
\insertRef{laverny2020}{cort}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cbCopula.R
\name{cbCopula-Class}
\alias{cbCopula-Class}
\alias{cbCopula}
\title{Checkerboard copulas}
\usage{
cbCopula(x, m = rep(nrow(x), ncol(x)), pseudo = FALSE)
}
\arguments{
\item{x}{the data to be used}

\item{m}{checkerboard parameters}

\item{pseudo}{Boolean, defaults to \code{FALSE}. Set to \code{TRUE} if you are already
providing pseudo data into the \code{x} argument.}
}
\value{
An instance of the \code{cbCopula} S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, etc.
}
\description{
cbCopula contructor
}
\details{
The cbCopula class computes a checkerboard copula with a given checkerboard parameter \eqn{m}, as described by A. Cuberos, E. Masiello and V. Maume-Deschamps (2019).
Assymptotics for this model are given by C. Genest, J. Neslehova and R. bruno (2017). The construction of this copula model is as follows :

Start from a dataset with \eqn{n} i.i.d observation of a \eqn{d}-dimensional copula (or pseudo-observations), and a checkerboard parameter \eqn{m},dividing \eqn{n}.

Consider the ensemble of multi-indexes \eqn{I = \{i = (i_1,..,i_d) \subset \{1,...,m \}^d\}} which indexes the boxes :

\deqn{B_{i} = \left]\frac{i-1}{m},\frac{i}{m}\right]}

Let now \eqn{\lambda} be the dimension-unspecific lebesgue measure on any power of \eqn{R}, that is :

\deqn{\forall d \in N, \forall x,y \in R^p, \lambda(\left(x,y\right)) = \prod\limits_{p=1}^{d} (y_i - x_i)}

Let furthermore \eqn{\mu} and \eqn{\hat{\mu}} be respectively the true copula measure of the sample at hand and the classical Deheuvels empirical copula, that is :
\itemize{
\item For \eqn{n} i.i.d observation of the copula of dimension \eqn{d}, let \eqn{\forall i \in \{1,...,d\}, \, R_i^1,...,R_i^d} be the marginal ranks for the variable \eqn{i}.
\item \eqn{\forall x \in I^d} let \eqn{\hat{\mu}((0,x)) = \frac{1}{n} \sum\limits_{k=1}^n I_{R_1^k\le x_1,...,R_d^k\le x_d}}
}

The checkerboard copula, \eqn{C}, and the empirical checkerboard copula, \eqn{\hat{C}}, are then defined by the following :

\deqn{\forall x \in (0,1)^d, C(x) = \sum\limits_{i\in I} {m^d \mu(B_{i}) \lambda((0,x)\cap B_{i})}}

Where \eqn{m^d = \lambda(B_{i})}.

This copula is a special form of patchwork copulas, see F. Durante, J. Fernández Sánchez and C. Sempi (2013) and F. Durante, J. Fernández Sánchez, J. Quesada-Molina and M. Ubeda-Flores (2015).
The estimator has the good property of always being a copula.

The checkerboard copula is a kind of patchwork copula that only uses independent copula as fill-in, only where there are values on the empirical data provided.
To create such a copula, you should provide data and checkerboard parameters (depending on the dimension of the data).
}
\references{
\insertRef{cuberos2019}{cort}

\insertRef{genest2017}{cort}

\insertRef{durante2013}{cort}

\insertRef{durante2015}{cort}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{loss}
\alias{loss}
\alias{loss,Cort-method}
\title{Loss of a copula estimation (if the model has one)}
\usage{
loss(object)

\S4method{loss}{Cort}(object)
}
\arguments{
\item{object}{the copula object}
}
\value{
the Integrated square error loss of the model
}
\description{
Currently only implemented for Cort models.
Compute the loss of the model.
}
\section{Functions}{
\itemize{
\item \code{loss,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(cort::recoveryourself_data[1:10,])
loss(cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R, R/CortForest.R,
%   R/cbCopula.R
\name{dCopula}
\alias{dCopula}
\alias{dCopula,matrix,Cort-method}
\alias{dCopula,matrix,CortForest-method}
\alias{dCopula,matrix,cbCopula-method}
\title{Copula density}
\usage{
dCopula(u, copula, ...)

\S4method{dCopula}{matrix,Cort}(u, copula)

\S4method{dCopula}{matrix,CortForest}(u, copula)

\S4method{dCopula}{matrix,cbCopula}(u, copula)
}
\arguments{
\item{u}{numeric matrix : one row per observation}

\item{copula}{the copula object}

\item{...}{other parameter to be passed to methods for this generic.}
}
\value{
The density of the copula on each observation
}
\description{
This function returns the density of a given copula on given observations.
}
\section{Functions}{
\itemize{
\item \code{dCopula,matrix,Cort-method}: Method for the class Cort

\item \code{dCopula,matrix,CortForest-method}: Method for the class CortForest

\item \code{dCopula,matrix,cbCopula-method}: Method for the cbCopula
}}

\examples{
cop <- cbCopula(cort::funcdep_data[1:10,1:2], m = 5)
dCopula(rep(0,2),cop)
dCopula(rep(0.5,2),cop)
dCopula(rep(1,2),cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{quad_prod}
\alias{quad_prod}
\alias{quad_prod,Cort,Cort-method}
\title{Quadratic product of two copulas (if they have one)}
\usage{
quad_prod(object, other_tree)

\S4method{quad_prod}{Cort,Cort}(object, other_tree)
}
\arguments{
\item{object}{: the tree}

\item{other_tree}{: the other tree}
}
\value{
the quadratic product between the trees
}
\description{
Currently only implemented for Cort models.
Compute the L2 quadratic product of 2 trees
}
\section{Functions}{
\itemize{
\item \code{quad_prod,Cort,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(LifeCycleSavings[,1:3])
all.equal(quad_prod(cop,cop),quad_norm(cop))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R, R/CortForest.R
\name{constraint_infl}
\alias{constraint_infl}
\alias{constraint_infl,Cort-method}
\alias{constraint_infl,CortForest-method}
\title{Constraint influence of the model (if it has one)}
\usage{
constraint_infl(object)

\S4method{constraint_infl}{Cort}(object)

\S4method{constraint_infl}{CortForest}(object)
}
\arguments{
\item{object}{the copula object}
}
\value{
The constraint influence statistic of the model
}
\description{
Currently only implemented for Cort models.
Compute the constraint influence of the model
}
\section{Functions}{
\itemize{
\item \code{constraint_infl,Cort-method}: Method for the class Cort

\item \code{constraint_infl,CortForest-method}: Method for the class CortForest
}}

\examples{
cop <- Cort(cort::recoveryourself_data[1:10,])
constraint_infl(cop)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Cort.R
\name{project_on_dims}
\alias{project_on_dims}
\alias{project_on_dims,Cort-method}
\title{Projection on smaller dimensions of a copula (if implemented)}
\usage{
project_on_dims(object, dims)

\S4method{project_on_dims}{Cort}(object, dims)
}
\arguments{
\item{object}{: the tree}

\item{dims}{the set of dimensions}
}
\value{
other cort object
}
\description{
Currently only implemented for Cort models.
Compute, as a Cort object, the projection on a smaller set of dimensions of a Cort object.
}
\section{Functions}{
\itemize{
\item \code{project_on_dims,Cort-method}: Method for the class Cort
}}

\examples{
cop <- Cort(LifeCycleSavings[,1:3])
projection = project_on_dims(cop,c(1,2))

}
