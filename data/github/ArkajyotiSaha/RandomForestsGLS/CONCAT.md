<!-- badges: start -->
[![R-CMD-check](https://github.com/ArkajyotiSaha/RandomForestsGLS/workflows/R-CMD-check/badge.svg)](https://github.com/ArkajyotiSaha/RandomForestsGLS/actions)
<!-- badges: end -->

# RandomForestsGLS
====

## Overview
The R package `RandomForestsGLS: Random Forests for dependent data` fits non-linear regression models on dependent data with Generalized Least Square (GLS) based Random Forest (RF-GLS). Classical Random forests ignore the correlation structure in the data for purpose of greedy partition, mean estimation and resampling for each regression tree. The package implements a RF-GLS proposed in Saha et al. (2021) which circumvents the aforementioned problems by incorporating a working correlation structure of the data for partition, mean estimation and resampling. In this article, it is shown that the greedy split criterion of classical regression trees can be written as an Ordinary Least Square (OLS) optimization with membership in current leaf nodes forming the design matrix. The article extends RF to RF-GLS in a similar fashion to how OLS is extended to GLS by incorporating the covariance structure of the data in the cost function. This ensures that the node splitting and node representatives involve contribution from points belonging to other nodes, weighed by their respective spatial correlations. In classical Random Forest (RF), data points are resampled/subsampled for each of the regression trees, without accounting for their inherent correlation structure. RF-GLS circumvents this problem by resampling/subsampling uncorrelated contrasts instead of original data points. `RandomForestsGLS` implements a fast version of RF-GLS which approximates the working correlation structure using Nearest Neighbor Gaussian Process (NNGP) (Datta et al., 2016) which makes it suitable for larger datasets.

## Installation
In order to install the development version of the package, please run the following command in R:

```{r }
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ArkajyotiSaha/RandomForestsGLS", ref ="HEAD")
```
For installation of the CRAN version of the package, please use the following command in R:

```{r}
install.packages("RandomForestsGLS")
```

## Example usage: Vignette
The package vignette, available at https://cran.rstudio.com/web/packages/RandomForestsGLS/vignettes/RandomForestsGLS_user_guide.pdf demonstrates with example how the functions available in `RandomForestsGLS` can be used for non-linear regression analysis of dependent data. Specific functions are discussed in much detail in the code documentation of the package. 

## Function documentation

For detailed help on the functions in `RandomForestsGLS` please use the following:
```{r }
?RFGLS_estimate_spatial #(for estimation in spatial data)
?RFGLS_estimate_timeseries #(for estimation in timeseries data)
?RFGLS_predict #(for prediction of mean function)
?RFGLS_predict_spatial #(for prediction of Spatial Response)
```
The function input and outputs are described in detail in the reference manual documentation, available in https://cran.rstudio.com/web/packages/RandomForestsGLS/RandomForestsGLS.pdf .

## Unknown model parameters

### Spatial data
More often than not, the true covariance parameters of the data are not known apriori. Given a choice from a prespecified list of covariance functions (software presently allows for exponential, spherical, Matérn and Gaussian covariances), the software accommodates parameter estimation with the argument `param_estimate = TRUE`. For this, first we fit a classical RF ignoring the correlation structure. Next, we fit a zero mean Gaussian process with the covariance structure of choice on the residuals to obtain the estimate of the model parameters through maximum likelihood estimation. This initial RF fitting and estimation of the spatial parameters from the RF residuals is justified by theoretical results in Saha et al. (2021) proving that even under dependence, naive RF is consistent (although empirically suboptimal). The estimated model parameters are then used to fit RF-GLS. As demonstrated in Saha et al. (2021), this leads to significant improvement over classical RF both in terms of estimation and prediction accuracy. This procedure is again analogous to the linear model setup, where fitting a feasible GLS model  often involves pre-estimating the covariance parameters from residuals from an OLS fit.

It has also been demonstrated (Section 4.4 in Saha et al. (2021)) that even under misspecification in covariance structure (i.e. when the true spatial effects are generated from a different covariance structure than the specified covariance for estimation or arbitrary smooth function), using RF-GLS with exponential covariance (the default choice of covariance function in the software) leads to noticeable improvement over simply using classical RF. 

One additional consideration for model parameter estimation involves the computational complexity of the maximum likelihood estimation of Gaussian process. In its original form, this involves $O(n^3)$ computation and $O(n^2)$ storage space, which may be prohibitive for large $n$. This problem is thoroughly studied in spatial statistics literature (we refer the readers to Datta et al. (2016) for a detailed review). In the present scenario, we perform a fast, linear-time estimation of the model parameters with BRISC (Saha and Datta, 2018), which is directly implemented in the present software through the R package BRISC (https://CRAN.R-project.org/package=BRISC).


### Time-series data
For time-series data, most often the model parameters for the autoregressive process (order of autoregression and the autocorrelation coefficients) are unknown. These are estimated using a strategy similar to the spatial case. With the option `param_estimate = TRUE`, the software estimates the model coefficients by fitting a zero mean autoregressive process of user defined order on the residuals from a naive RF fit. The autoregressive parameter fitting is done using `arima` from base R package `stats`.


## Parallelization

For `RFGLS_estimate_spatial`, `RFGLS_estimate_timeseries`, `RFGLS_predict` and `RFGLS_predict_spatial` one can also take the advantage of parallelization, contingent upon the availability of multiple cores. One aspect of the parallelizatiion is through the multithreaded implementation of derivation of the NNGP (Datta et al., 2016) components following `BRISC` package (which helps when we are dealing with a large dataset). We also run the regression trees on parallel, which can be beneficial when the number of trees is large. With very small dataset and small number of trees, communication overhead between the nodes for parallelization outweighs the benefits of the parallel computing. Hence it is recommended to parallelize only for moderately large dataset and/or large number of trees.


# Package Features
* **Implementation**: The source code of the package are written in [C](https://en.cppreference.com/w/c/language)/[C++](https://isocpp.org/) for sake of optimizing execution time. The functions available to the user are wrappers around the source code, built with `R`'s foreign language interface. For the basic structure of the code, we make use of the open-source code of the regression trees in `R` based implementation of classical RF in `randomForest` package. As the split criterion in RF-GLS involves computationally intensive linear algebra operation in nested loops, we use `Fortran`'s Basic Linear Algebra Subprograms ([BLAS](http://www.netlib.org/blas/)) and Linear Algebra Package ([LAPACK](http://www.netlib.org/lapack/)). This is achieved by storing all matrices in contiguous memory column-major format. We also offer multicore computation by building each regression tree independently.

* **NNGP approximation**: Node splitting in RF-GLS requires optimizing a cost function involving the Cholesky factor of the precision matrix. Use of the full dense precision matrix in spatial processes becomes taxing on typical personal computing resources both in terms of computational cost ($O(n^3)$) and storage cost ($O(n^2)$). In order to circumvent this problem, we use NNGP (Datta et al., 2016) to replace the dense graph among spatial locations with a nearest neighbor graphical model. NNGP components can be combined to obtain a sparse Cholesky factor, which closely approximates the decorrelation performance of the true Cholesky. We implement a convenient nearest neighbor search following [spNNGP](https://CRAN.R-project.org/package=spNNGP) (Finley et al., 2020) and efficient sparse matrix multiplication as in [BRISC](https://CRAN.R-project.org/package=BRISC) (Saha and Datta, 2018). The structure of the loops used in the process facilitates parallelization using `openMP` (Dagum and Menon, 1998) for this stage of calculation. In time series analysis, the sparsity in the precision matrix is inherently induced by AR covariance structure. 

* **Scalable node splitting**: Another aspect of optimization of the proposed algorithm involves resourceful implementation of the cost function optimization. Provided candidate cut direction ($d$), the optimal cutoff point ($c$) is chosen by searching through the "gaps" in the corresponding covariate. Following the implementation of classical Regression Forest in [randomForest](https://CRAN.R-project.org/package=randomForest), we start with a list of ordered covariate values corresponding to the prefixed candidate direction and assign the associated data points to one of the nodes initially. This helps with searching through the "gaps" in the data, as searching the next "gap" is equivalent to switch of membership of the existing smallest member (w.r.t the covariate value in the prefixed direction) of the initially assigned node. In order to determine the cost function corresponding to each of the potential cutoff points, in each iteration, we serially switch the membership of the data points from the initial node. The process is graphically demonstrated in the following figure.

![Serial update of membership and cutoff .\label{fig:example}](figure.png)

Since a serial update only affects one row in two columns of design matrix **Z** (corresponding to the newly formed nodes, an example of changes in **Z** corresponding to the changes in the figure above is demonstrated in the following figure, here we note that the covariate values corresponding to the matrix components are unordered, hence the serial update will not follow the natural ordering), the resulting correlation adjusted effective design matrix (**Q**<sup>1/2</sup> **Z**) used in GLS loss, only experiences changes in the corresponding two columns and the rows corresonding to the points, which have this specific point in their nearest neighbor set. If we consider the specific example in the following figure, location "3" is a nearest neighbor of only locations "4" and "6". Hence the update rule only effects rows 3, 4 and 6. We efficiently implement this in the package which provides efficiency over brute force recomputation of the effective design matrix for each serial update. The process is graphically demonstrated in the following figure.

![Serial update of Z .\label{fig:example2}](figure2.png)

![Changes in correlation adjusted design matrix .\label{fig:example3}](figure3.png)

## Community guidelines

Please report issues, bugs or problem with the software at https://github.com/ArkajyotiSaha/RandomForestsGLS/issues . For contribution to the software and support please get in touch with the maintainer Arkajyoti Saha (arkajyotisaha93@gmail.com).

## Note
Some code blocks are borrowed from the R packages: `spNNGP: Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes` https://CRAN.R-project.org/package=spNNGP and `randomForest: Breiman and Cutler's Random Forests for Classification and Regression` https://CRAN.R-project.org/package=randomForest .
RF-GLS uses nearest neighbor Gaussian process (NNGP) to approximate the covariance structure in the data, tools necessary to implement the NNGP are borrowed from the spNNGP package, which include `util.cpp` and parts of `updateBF_org` and `RFGLS_BFcpp` in `RFGLS.cpp`. The basic building blocks for Random Forest is borrowed from `randomForest` which include parts of `RFGLStree_cpp`, `findBestSplit` `RFGLSpredicttree_cpp` in `RFGLS.cpp`.


## Citation
Please cite the following paper when you use RF-GLS

Saha, Arkajyoti, Sumanta Basu, and Abhirup Datta. "Random forests for spatially dependent data." Journal of the American Statistical Association (2021): 1-19. https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1950003 .


## References

Dagum, Leonardo, and Ramesh Menon. "OpenMP: an industry standard API for shared-memory programming." IEEE computational science and engineering 5, no. 1 (1998): 46-55. 

Datta, Abhirup, Sudipto Banerjee, Andrew O. Finley, and Alan E. Gelfand. "Hierarchical nearest-neighbor Gaussian process models for large geostatistical datasets." Journal of the American Statistical Association 111, no. 514 (2016): 800-812. https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091 .

Saha, Arkajyoti, and Abhirup Datta. "BRISC: bootstrap for rapid inference on spatial covariances." Stat 7, no. 1 (2018): e184. https://onlinelibrary.wiley.com/doi/10.1002/sta4.184 .

Finley, Andrew O., Abhirup Datta, and Sudipto Banerjee. "spNNGP R package for nearest neighbor Gaussian process models." arXiv preprint arXiv:2001.09111 (2020).

Saha, Arkajyoti, Sumanta Basu, and Abhirup Datta. "Random forests for spatially dependent data." Journal of the American Statistical Association (2021): 1-19. https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1950003 .
---
title: 'RandomForestsGLS: An R package for Random Forests for dependent data'
tags:
  - R
  - spatial statistics
  - Gaussian Processes
  - Random forests
  - generalized least-squares
authors:
  - name: Arkajyoti Saha
    affiliation: 1
  - name: Sumanta Basu
    affiliation: 2
  - name: Abhirup Datta
    affiliation: 3
affiliations:
 - name: Departments of Statistics, University of Washington
   index: 1
 - name: Department of Statistics and Data Science, Cornell University
   index: 2
 - name: Department of Biostatistics, Johns Hopkins Bloomberg School of Public Health
   index: 3
date: 23 February 2022
bibliography: paper.bib
---

# Summary
With the modern advances in geographical information systems, remote sensing technologies, and low-cost sensors, we are increasingly encountering datasets where we need to account for spatial or serial dependence. Dependent observations $(y_1, y_2, \cdots, y_n)$ with covariates $(\mathbf x_1,\ldots,\mathbf x_n)$ can be modeled non-parametrically as $y_i = m(\mathbf x_i) + \epsilon_i$, where $m(\mathbf x_i)$ is mean component and $\epsilon_i$ accounts for the dependency in data. We assume that dependence is captured through a covariance function of the correlated stochastic process $\epsilon_i$ (second order dependence). The correlation is typically a function of "spatial distance" or "time-lag" between two observations. 

Unlike linear regression, non-linear Machine Learning (ML) methods for estimating the regression function $m$ can capture complex interactions among the variables. However, they often fail to account for the dependence structure, resulting in sub-optimal estimation. On the other hand, specialized software for spatial/temporal data properly models data correlation but lacks flexibility in modeling the mean function $m$ by only focusing on linear models. `RandomForestsGLS` bridges the gap through a novel rendition of Random Forests (RF) -- namely, RF-GLS -- by explicitly modeling the spatial/serial data correlation in the RF fitting procedure to substantially improve the estimation of the mean function. Additionally, `RandomForestsGLS` leverages kriging to perform predictions at new locations for geo-spatial data.

# Statement of need
`RandomForestsGLS` is a statistical machine learning oriented [R](https://cran.r-project.org) package for fitting RF-GLS on dependent data. The RF-GLS algorithm described in [@saha2021random] involves computationally intensive linear algebra operations in nested loops which are especially slow in an interpreted language like `R`. `RandomForestsGLS` efficiently implements RF-GLS algorithm in a low-level language with a user-friendly interface in `R`, a popular computational software (free under GNU General Public License) in the statistics community. `RandomForestsGLS` focuses on fast, parallelizable implementations of RF-GLS for spatial and time series data, which includes popular choices for covariance functions for both spatial (Matérn GP) and time series (autoregressive) data. The package is primarily designed to be used by researchers associated with the fields of statistical machine learning, spatial statistics, time series analysis, and their scientific applications. A significant part of the code has already been used in @saha2021random. With the combination of speed and ease-of-use that `RandomForestsGLS` brings to the table regarding non-linear regression analysis for dependent data, we hope to see this package being used in a plethora of future scientific and methodological explorations.

# State of the field

Several `R` packages implement classical RF. Most notable of them being [randomForest](https://CRAN.R-project.org/package=randomForest), which implements Breiman's Random Forests [@breiman2001random] for Classification and
Regression using the [Fortran](https://fortran-lang.org/). Some of the other packages are [xgboost](https://CRAN.R-project.org/package=xgboost), [randomForestSRC](https://CRAN.R-project.org/package=randomForestSRC), [ranger](https://CRAN.R-project.org/package=ranger), [Rborist](https://CRAN.R-project.org/package=Rborist). For a detailed overview of it, we refer the reader to [CRAN Task View: Machine Learning & Statistical Learning](https://cran.r-project.org/web/views/MachineLearning.html). To the best of our knowledge, none of these packages explicitly account for spatial and/or temporal correlation.

Classical RF has been used in geo-spatial and temporal applications (see @saha2021random for references) without making methodological adjustments to account for spatial dependencies. ([CRAN Task View: Analysis of Spatial Data](https://cran.r-project.org/web/views/Spatial.html), [CRAN Task View: Time Series Analysis](https://cran.r-project.org/web/views/TimeSeries.html)). Two recent works that attempt to explicitly use spatial information in RF for prediction purposes, are @hengl2018random ( [GeoMLA](https://github.com/thengl/GeoMLA)) and @georganos2019geographical ( [SpatialML](https://CRAN.R-project.org/package=SpatialML)) (see @saha2021random for details). Both approaches try to account for the dependence structure by incorporating additional spatial covariates, which adversely affect the prediction performance in the presence of a dominant covariate effect (@saha2021random). Additionally, unlike RF-GLS, they cannot estimate covariate effects separately from the spatial effect, which can be of independent interest. The RF for temporal data, proposed in @BASAK2019552, suffers from similar shortcomings.

# The RandomForestsGLS package

We provide a brief overview of the functionality of the package. The package [vignette](https://cran.r-project.org/web/packages/RandomForestsGLS/vignettes/RandomForestsGLS_user_guide.pdf) demonstrates with example how to use the package for non-linear regression analysis of dependent data. Specific functions are discussed in detail in the documentation of the package.

## RF to RF-GLS: Accounting for correlation structure

In classical RF, which is an average of many regression trees, each node in a regression tree is split by optimizing the CART split criterion in @breiman1984classification. It can be rewritten in the following way: 
$$
v_{n}^{CART}((d,c)) =  \frac{1}{n} \left( \|\mathbf{Y} - \mathbf Z^{(0)}\boldsymbol{\hat{\beta}}(\mathbf Z^{(0)})\|_2^2 - \|\mathbf Y - \mathbf Z \boldsymbol{\hat{\beta}}(\mathbf Z)\|_2^2 \right).
$$

where $\mathbf Z^{(0)}$ and $\mathbf Z$ are the binary membership matrices for the leaf nodes of the tree before and after the potential node split. $(d,c)$ denotes a potential cut (location of the split), with $d$ and $c$ being the cut direction (choice of the covariate) and cutoff point (value of the covariate) respectively, $\boldsymbol{\hat{\beta}} (\mathbf Z)$ are the leaf node representatives given by OLS estimates corresponding to a design matrix $\mathbf Z$ and can be written as: 

$$\boldsymbol{\hat{\beta}} (\mathbf Z) = \left(\mathbf Z ^\top \mathbf Z \right)^{-1} \mathbf Z ^\top \mathbf y$$

We observe that the split criterion is the difference of OLS loss functions before and after the cut with the respective design matrices of membership of the leaf nodes. We can incorporate the correlation structure of the data in the split criterion by replacing the OLS loss with GLS loss as is traditionally done in linear models. The modified split criterion can be rewritten as:

$$
\begin{aligned}
v_{n,\mathbf Q}^{DART}((d,c)) = 
&\frac{1}{n} \Bigg[\left(\mathbf{Y} - \mathbf{Z}^{(0)}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z^{(0)}) \right)^\top \mathbf Q\left(\mathbf{Y} - \mathbf{Z}^{(0)}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z^{(0)}) \right)\\ &-\left(\mathbf{Y} - \mathbf{Z}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z) \right)^\top \mathbf Q\left(\mathbf{Y} - \mathbf{Z}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z) \right) \Bigg].
\end{aligned}
$$

where $\mathbf Q$ is the inverse of the working covariance matrix that models the spatial/serial dependence and $\boldsymbol{\hat{\beta}}_{GLS} (\mathbf Z)$ are the leaf node representatives given by the GLS estimates corresponding to a design matrix $\mathbf Z$ and can be written as follows:

$$\boldsymbol{\hat{\beta}}_{GLS} (\mathbf Z) = \left(\mathbf Z ^\top \mathbf Q \mathbf Z \right)^{-1} \mathbf Z ^\top \mathbf Q \mathbf y.$$

## Spatial Data
### Model
We consider spatial point-referenced data with the following mixed model:

$$
y_i = m(\mathbf{x}_i) + w(\mathbf{s}_i) + \epsilon_i;
$$

where $y_i, \mathbf{x}_i$ respectively denotes the observed response and the covariate corresponding to the $i^{th}$ observed location $\mathbf{s}_i$. $m(\mathbf{x}_i)$ denotes the covariate effect; spatial random effect, $w (\mathbf{s})$ accounts for spatial dependence beyond covariates modeled using a GP, and $\mathbf{\epsilon}$ accounts for independent and identically distributed random Gaussian noise. 



### Fitting & Prediction
Spatial random effects are modeled using a Gaussian process (GP) as is the practice. We use the computationally convenient Nearest Neighbor GP (NNGP) [@nngp]. Model parameters, if unknown, are estimated from the data (@saha2021random). Alongside predicting covariate effects for a new covariate value, we also offer spatial prediction at new locations with non-linear kriging by combining the non-linear mean estimate and spatial kriging estimate from the [BRISC](https://CRAN.R-project.org/package=BRISC) [@brisc] package.


## Autoregressive (AR) Time Series Data
### Model

RF-GLS can also be used for function estimation in a time series setting under autoregressive (AR) errors. We consider time series data with errors from an AR($q$) (autoregressive process of lag $q$) process as follows:

$$
y_t = m(\mathbf{x}_t) + e_t;\: e_t = \sum_{i = 1}^q\rho_i e_{t-i} + \eta_t
$$

where $y_i, \mathbf{x}_i$ denote the response and the covariate corresponding to the $t^{th}$ time point, respectively, $e_t$ is an AR(q) process, $\eta_t$ denotes i.i.d. white noise and $(\rho_1, \cdots, \rho_q)$ are the model parameters that capture the dependence of $e_t$ on $(e_{t - 1}, \cdots, e_{t-q})$.

### Fitting & Prediction
RF-GLS exploits the sparsity of the closed form precision matrix of the AR process for model fitting and prediction of the mean function $m(.)$. If the AR model parameters (coefficients and order of autoregressive process) are unknown, the code automatically estimates them from AR models with specified lags. The prediction of covariate effect for time series data is like that of spatial data. 

# Discussion

This package provides an efficient, parallel implementation of the RF-GLS method proposed in @saha2021random which accounts for correlated data with modified node splitting criteria and node representative update rule. The package accounts for spatial correlation via Matérn GP or serial autocorrelation. It provides parameter estimation for both covariance structures. Efficient implementation through C/C++ takes advantage of the NNGP approximation and scalable node splitting update rules which reduces the execution time. More details on package features, parameter estimation, and parallelization can be found in the package [README](https://github.com/ArkajyotiSaha/RandomForestsGLS#readme). Since the computational complexity of evaluating potential splits is cubic in the number of leaf nodes for RF-GLS, but constant for standard RF, improving the computational complexity of the algorithm is of independent research interest. We also plan to implement and validate additional forms of time series dependency.

# Acknowledgements

AS and AD were supported by NSF award DMS-1915803. AD was supported by NIEHS award R01ES033739. SB was supported by an NSF award DMS-1812128, and an NIH award R01GM135926.

# References

---
title: "How to use RandomForestsGLS"
output: rmarkdown::pdf_document
vignette: >
  %\VignetteIndexEntry{How to use RandomForestsGLS}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The package `RandomForestsGLS` fits non-linear regression models on dependent data with Generalised Least Square (GLS) based Random Forest (RF-GLS) detailed in Saha, Basu and Datta (2021) https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1950003. We will start by loading the `RandomForestsGLS` R package.

```{r setup}
library(RandomForestsGLS)
```

Next, we discuss how the `RandomForestsGLS` package can be used for estimation and prediction in a non-linear regression setup under correlated errors in different scenarios.

# 1. Spatial Data

We consider spatial point referenced data with the following model:

$$
y_i = m(\mathbf{x}_i) + w(\mathbf{s}_i) + \epsilon_i;
$$
where, $y_i, \mathbf{x}_i$ respectively denotes the observed response and the covariate corresponding to the $i^{th}$ observed location $\mathbf{s}_i$. $m(\mathbf{x}_i)$ denotes the covariate effect, spatial random effect, $w (\mathbf{s})$ accounts for spatial dependence beyond covariates, and $\mathbf{\epsilon}$ accounts for the independent and identically distributed random Gaussian noise. 

In the spatial mixed model setting, the package `RandomForestsGLS` allows for fitting $m(.)$ using RF-GLS. Spatial random effects are modeled using Gaussian Process as is the practice. For model fitting, we use the computationally convenient Nearest Neighbor Gaussian Process (NNGP) (Datta, Banerjee, Finley, and Gelfand (2016)) for the spatial random effects $w(\cdot)$. Along with prediction of the covariate effect (mean function) $m(.)$ we also offer kriging based prediction of spatial responses at new location. 


## Illustration

We simulate a data from the following model:
$$
y_i = 10\sin(\pi x_i) + w (\mathbf{s}_i)+ \epsilon_i; \:\: \epsilon \sim N(\mathbf{0},\: \tau^2 \mathbf{I}), \tau^2 = 0.1; \:\:\: w \sim \textit{exponential GP};\: \sigma^2 = 10; \phi = 1.
$$

Here, the mean function is $E(Y) = 10\sin(\pi X)$; $w$ accounts for the spatial correlation, which is generated as a exponential Gaussian process with spatial variance $\sigma^2 = 10$ and spatial correlation decay $\phi = 1$; and $\epsilon$ is the i.i.d random noise with variance $\tau^2 = 0.1$, which is also called the nugget in spatial literature.

For illustration purposes, we simulate with $n = 200$:

```{r simulation, message=FALSE, warning=FALSE}
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(5)
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))
set.seed(2)
x <- as.matrix(runif(n),n,1)
sigma.sq = 10
phi = 1
tau.sq = 0.1
D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))
```



##  Model fitting

In the package `RandomForestsGLS`, the working precision matrix used in the GLS-loss are NNGP approximations of precision matrices corresponding to Matérn covariance function. 

In order to fit the model, the code requires: 

+ Coordinates (`coords`): an $n \times 2$ matrix of 2-dimensional locations.
+ Response (`y`):  an $n$ length vector of response at the observed coordinates.
+ Covariates (`X`): an $n \times p$ matrix of the covariates in the observation coordinates.
+ Covariates for estimation (`Xtest`): an $ntest \times p$ matrix of the covariates where we want to estimate the function. Must have identical variables as that of `X`. Default is `X`.
+ Minimum size of leaf nodes (`nthsize`): We recommend not setting this value too small, as that will lead to very deep trees that takes a lot of time to be built and can produce unstable estimates. Default value is 20.
+ The parameters corresponding to the covariance function (detailed afterwards).

For the details on choice of other parameters, please refer to the help file of the code `RFGLS_estimate_spatial`, which can be accessed with `?RFGLS_estimate_spatial`. 

### Known Covariance Parameters

If the covariance parameters are known, we set `param_estimate = FALSE` (default value); the code additionally requires the following:

+ Covariance Model (`cov.model`): Supported keywords are: "exponential", "matern", "spherical", and "gaussian" for exponential, Matérn, spherical and Gaussian covariance function respectively. Default value is "exponential".
+ $\sigma^2$ (`sigma.sq`): The spatial variance. Default value is 1.
+ $\tau^2$ (`tau.sq`): The nugget. Default value is 0.01.
+ $\phi$ (`phi`): The spatial correlation decay parameter. Default value is 5.
+ $\nu$ (`nu`): The smoothing parameter corresponding to the Matérn covariance function. Default value is 0.5. 

We can fit the model as follows:

```{r model_fit, message=FALSE, warning=FALSE}
set.seed(1)
est_known <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential",
                                    nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq, 
                                    phi = phi)
```
The estimate of the function at the covariates `Xtest` is given in `estimation_reult$predicted`. For interpretation of the rest of the outputs, please see the help file of the code `RFGLS_estimate_spatial`. Using covariance models other than exponential model are in beta testing stage.

### Unknown Covariance Parameters

If the covariance parameters are not known we set `param_estimate = TRUE`; the code additionally requires the covariance model (`cov.model`) to be used for parameter estimation prior to RF-GLS fitting. We fit the model with unknown covariance parameters as follows.  

```{r model_fit_unknown, message=FALSE, warning=FALSE}
set.seed(1)
est_unknown <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential", 
                                      nthsize = 20, param_estimate = TRUE)
```



## Prediction of mean function

Given a fitted model using `RFGLS_estimate_spatial`, we can estimate the mean function at new covariate values as follows:

```{r model_estimate, message=FALSE, warning=FALSE}
Xtest <- matrix(seq(0,1, by = 1/10000), 10001, 1)
RFGLS_predict_known <- RFGLS_predict(est_known, Xtest)
```

### Performance comparison

We obtain the Mean Integrated Squared Error (MISE) of the estimate $\hat{m}$ from RF-GLS on [0,1] and compare it with that corresponding to the classical Random Forest (RF) obtained using package `randomForest` (with similar minimum nodesize, `nodesize = 20`, as default `nodesize` performs worse). We see that our method has a significantly smaller MISE. Additionally, we show that the MISE obtained with unknown parameters in RF-GLS is comparable to that of the MISE obtained with known covariance parameters.

```{r comparison, message=FALSE, warning=FALSE}
library(randomForest)
set.seed(1)
RF_est <- randomForest(x, y, nodesize = 20)

RF_predict <- predict(RF_est, Xtest)
#RF MISE
mean((RF_predict - 10*sin(pi * Xtest))^2)

#RF-GLS MISE
mean((RFGLS_predict_known$predicted - 10*sin(pi * Xtest))^2)

RFGLS_predict_unknown <- RFGLS_predict(est_unknown, Xtest)
#RF-GLS unknown MISE
mean((RFGLS_predict_unknown$predicted - 10*sin(pi * Xtest))^2)
```

We plot the true $m(x) = 10sin(\pi x)$ along with the loess-smoothed version of estimated $\hat{m}(.)$ obtained from RF-GLS and RF where we show that RF-GLS estimate approximates $m(x)$ better than that corresponding to RF.

```{r plot_comparison, message=FALSE, warning=FALSE}
rfgls_loess_10 <- loess(RFGLS_predict_known$predicted ~ c(1:length(Xtest)), span=0.1)
rfgls_smoothed10 <- predict(rfgls_loess_10)

rf_loess_10 <- loess(RF_predict ~ c(1:length(RF_predict)), span=0.1)
rf_smoothed10 <- predict(rf_loess_10)

xval <- c(10*sin(pi * Xtest), rf_smoothed10, rfgls_smoothed10)
xval_tag <- c(rep("Truth", length(10*sin(pi * Xtest))), rep("RF", length(rf_smoothed10)), 
              rep("RF-GLS",length(rfgls_smoothed10)))
plot_data <- as.data.frame(xval)
plot_data$Methods <- xval_tag
coval <- c(rep(seq(0,1, by = 1/10000), 3))
plot_data$Covariate <- coval

library(ggplot2)
ggplot(plot_data, aes(x=Covariate, y=xval, color=Methods)) +
geom_point() + labs( x = "x") + labs( y = "f(x)")
```


## Prediction of spatial response

Given a fitted model using `RFGLS_estimate_spatial`, we can predict the spatial response/outcome at new locations provided the covariates at that location. This approach performs kriging at a new location using the mean function estimates at the corresponding covariate values. Here we partition the simulated data into training and test sets in 4:1 ratio. Next we perform prediction on the test set using a model fitted on the training set.

```{r prediction_spatial, message=FALSE, warning=FALSE}
est_known_short <- RFGLS_estimate_spatial(coords[1:160,], y[1:160], 
                   matrix(x[1:160,],160,1), ntree = 50,  cov.model = "exponential", 
                   nthsize = 20, param_estimate = TRUE)
RFGLS_predict_spatial <- RFGLS_predict_spatial(est_known_short, coords[161:200,], 
                                               matrix(x[161:200,],40,1))
pred_mat <- as.data.frame(cbind(RFGLS_predict_spatial$prediction, y[161:200]))
colnames(pred_mat) <- c("Predicted", "Observed")
ggplot(pred_mat, aes(x=Observed, y=Predicted)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  ylim(0, 16) + xlim(0, 16)
```



## Misspecification in covariance model

The following example considers a setting when the parameters are estimated from a misspecified covariance model. We simulate the spatial correlation from a Matérn covariance function with smoothing parameter $\nu = 1.5$. While fitting the RF-GLS, we estimate the covariance parameters using an exponential covariance model ($\nu = 0.5$) and show that the obtained MISE can compare favorably to that of classical RF. 

```{r misspec_spatial, message=FALSE, warning=FALSE}
#Data simulation from matern with nu = 1.5
nu = 3/2
R1 <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
diag(R1) <- 1
set.seed(2)
w <- rmvn(1, rep(0,n), sigma.sq*R1)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

#RF-GLS with exponential covariance
set.seed(3)
est_misspec <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential", 
                                      nthsize = 20, param_estimate = TRUE)
RFGLS_predict_misspec <- RFGLS_predict(est_misspec, Xtest)

#RF
set.seed(4)
RF_est <- randomForest(x, y,  nodesize = 20)
RF_predict <- predict(RF_est, Xtest)

#RF-GLS MISE
mean((RFGLS_predict_misspec$predicted - 10*sin(pi * Xtest))^2)
#RF MISE
mean((RF_predict - 10*sin(pi * Xtest))^2)
```

# 2. Autoregressive Time Series Data

RF-GLS can also be used for function estimation in a time series setting under autoregressive errors. We consider time series data with errors from an AR(q) process as follows:

$$
y_t = m(\mathbf{x}_t) + e_t;\: e_t = \sum_{i = 1}^q\rho_i e_{t-i} + \eta_t
$$
where, $y_i, \mathbf{x}_i$ denotes the response and the covariate corresponding to the $t^{th}$ time point, $e_t$ is an AR(q) pprocess, $\eta_t$ denotes the i.i.d. white noise and $(\rho_1, \cdots, \rho_q)$ are the model parameters that captures the dependence of $e_t$ on $(e_{t - 1}, \cdots, e_{t-q})$.


In the AR time series scenario, the package `RandomForestsGLS` allows for fitting $m(.)$ using RF-GLS. RF-GLS exploits the sparsity of the closed form precision matrix of the AR process for model fitting and prediction of mean function $m(.)$. 

## Illustration

Here, we simulate from the AR(1) process as follows:
$$
y = 10\sin(\pi x) + \mathbf{e}; e_t = \rho e_{t-1} + \eta_t; \eta_t \sim N(0,\sigma^2); e_1 = \eta_1; \: \rho = 0.9; \sigma^2 = 10.
$$

Here, $E(Y) = 10\sin(\pi X)$; $\mathbf{e}$ which is an AR(1) process, accounts for the temporal correlation, $\sigma^2$ denotes the variance of white noise part of the AR(1) process and $\rho$ captures the degree of dependence of $e_t$ on $e_{t-1}$. 

For illustration purposes, we simulate with $n = 200$:

```{r simulation_temporal}
rho <- 0.9
set.seed(1)
b <- rho
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(1,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)
y <- c(eps + 10*sin(pi * x))
```


##  Model fitting

In case of time series data, the code requires: 

+ Response (`y`):  an $n$ length vector of response at the observed time points.
+ Covariates (`X`): an $n \times p$ matrix of the covariates in the observation time points.
+ Covariates for estimation (`Xtest`): an $ntest \times p$ matrix of the covariates where we want to estimate the function. Must have identical variables as that of `X`. Default is `X`.
+ Minimum size of leaf nodes (`nthsize`): We recommend not setting this value too small, as that will lead to very deep trees that takes a lot of time to be built and can produce unstable estimates. Default value is 20.
+ The parameters corresponding to the AR process (detailed afterwards).

For the details on choice of other parameters, please refer to the help file of the code `RFGLS_estimate_timeseries`, which can be accessed with `?RFGLS_estimate_timeseries`. 

### Known AR process Parameters

If the AR process parameters are known we set `param_estimate = FALSE` (default value); the code additionally requires `lag_params` = $c(\rho_1, \cdots, \rho_q)$.

We can fit the model as follows:

```{r model_fit_temporal, message=FALSE, warning=FALSE}
set.seed(1)
est_temp_known <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = rho, nthsize = 20)
```


### Unknown AR process Parameters

If the AR process parameters are not known, we set `param_estimate = TRUE`; the code requires the order of the AR process, which is obtained from the length of the `lag_params` input vector. Hence if we want to estimate the parameters from a AR(q) process, `lag_params` should be any vector of length `q`. Here we fit the model with `q = 1`

```{r model_fit_temporal_unknown, message=FALSE, warning=FALSE}
set.seed(1)
est_temp_unknown <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = rho, 
                                              nthsize = 20, param_estimate = TRUE)
```

## Prediction of mean function

This part of time series data analysis is identical to that corresponding to the spatial data.

```{r model_estimate_temporal, message=FALSE, warning=FALSE}
Xtest <- matrix(seq(0,1, by = 1/10000), 10001, 1)
RFGLS_predict_temp_known <- RFGLS_predict(est_temp_known, Xtest)
```

Here also, similar to the spatial data scenario, RF-GLS outperforms classical RF in terms of MISE both with true and estimated AR process parameters.

```{r comparison_temp, message=FALSE, warning=FALSE}
library(randomForest)
set.seed(1)
RF_est_temp <- randomForest(x, y, nodesize = 20)

RF_predict_temp <- predict(RF_est_temp, Xtest)
#RF MISE
mean((RF_predict_temp - 10*sin(pi * Xtest))^2)

#RF-GLS MISE
mean((RFGLS_predict_temp_known$predicted - 10*sin(pi * Xtest))^2)

RFGLS_predict_temp_unknown <- RFGLS_predict(est_temp_unknown, Xtest)
#RF-GLS unknown MISE
mean((RFGLS_predict_temp_unknown$predicted - 10*sin(pi * Xtest))^2)
```

## Misspecification in AR process order

We consider a scenario where the order of autoregression used for RF-GLS model fitting is mis-specified. We simulate the AR errors from an AR(2) process and fit RF-GLS with an AR(1) process.

```{r misspec_temporal, message=FALSE, warning=FALSE}
#Simulation from AR(2) process
rho1 <- 0.7
rho2 <- 0.2
set.seed(2)
b <- c(rho1, rho2)
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(2,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)
y <- c(eps + 10*sin(pi * x))

#RF-GLS with AR(1)
set.seed(3)
est_misspec_temp <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = 0, 
                                              nthsize = 20, param_estimate = TRUE)
RFGLS_predict_misspec_temp <- RFGLS_predict(est_misspec_temp, Xtest)

#RF
set.seed(4)
RF_est_temp <- randomForest(x, y,  nodesize = 20)
RF_predict_temp <- predict(RF_est_temp, Xtest)

#RF-GLS MISE
mean((RFGLS_predict_misspec_temp$predicted - 10*sin(pi * Xtest))^2)
#RF MISE
mean((RF_predict_temp - 10*sin(pi * Xtest))^2)
```


# Parallelization

For `RFGLS_estimate_spatial`, `RFGLS_estimate_timeseries`, `RFGLS_predict` and `RFGLS_predict_spatial` one can also take the advantage of parallelization, contingent upon the availability of multiple cores. The component `h` in all the functions determines the number of cores to be used. Here we demonstrate an example with `h = 2`.

```{r parallel_spatial, message=FALSE, warning=FALSE}
#simulation from exponential distribution
set.seed(5)
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))
set.seed(2)
x <- as.matrix(runif(n),n,1)
sigma.sq = 10
phi = 1
tau.sq = 0.1
nu = 0.5
D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

#RF-GLS model fitting and prediction with parallel computation
set.seed(1)
est_known_pl <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential",
                                    nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq, 
                                    phi = phi, h = 2)
RFGLS_predict_known_pl <- RFGLS_predict(est_known_pl, Xtest, h = 2)

#MISE from single core
mean((RFGLS_predict_known$predicted - 10*sin(pi * Xtest))^2)
#MISE from parallel computation
mean((RFGLS_predict_known_pl$predicted - 10*sin(pi * Xtest))^2)
```

For `RFGLS_estimate_spatial` with very small dataset (`n`) and small number of trees (`ntree`), communication overhead between the nodes for parallelization outweighs the benefits of the parallel computing hence it is recommended to parallelize only for moderately large `n` and/or `ntree`. It is strongly recommended that the max value of `h` is kept strictly less than the number of total cores available. Parallelization for `RFGLS_estimate_timeseries` can be addressed identically. For `RFGLS_predict` and `RFGLS_predict_spatial`, even for large dataset, single core performance is very fast, hence unless `ntest` and `ntree` are very high, we do not recommend using parallelization for `RFGLS_predict` and `RFGLS_predict_spatial`.
---
title: "How to use RandomForestsGLS"
output: rmarkdown::pdf_document
vignette: >
  %\VignetteIndexEntry{How to use RandomForestsGLS}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The package `RandomForestsGLS` fits non-linear regression models on dependent data with Generalised Least Square (GLS) based Random Forest (RF-GLS) detailed in Saha, Basu and Datta (2021) https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1950003. We will start by loading the `RandomForestsGLS` R package.

```{r setup}
library(RandomForestsGLS)
```

Next, we discuss how the `RandomForestsGLS` package can be used for estimation and prediction in a non-linear regression setup under correlated errors in different scenarios.

# 1. Spatial Data

We consider spatial point referenced data with the following model:

$$
y_i = m(\mathbf{x}_i) + w(\mathbf{s}_i) + \epsilon_i;
$$
where, $y_i, \mathbf{x}_i$ respectively denotes the observed response and the covariate corresponding to the $i^{th}$ observed location $\mathbf{s}_i$. $m(\mathbf{x}_i)$ denotes the covariate effect, spatial random effect, $w (\mathbf{s})$ accounts for spatial dependence beyond covariates, and $\mathbf{\epsilon}$ accounts for the independent and identically distributed random Gaussian noise. 

In the spatial mixed model setting, the package `RandomForestsGLS` allows for fitting $m(.)$ using RF-GLS. Spatial random effects are modeled using Gaussian Process as is the practice. For model fitting, we use the computationally convenient Nearest Neighbor Gaussian Process (NNGP) (Datta, Banerjee, Finley, and Gelfand (2016)) for the spatial random effects $w(\cdot)$. Along with prediction of the covariate effect (mean function) $m(.)$ we also offer kriging based prediction of spatial responses at new location. 


## Illustration

We simulate a data from the following model:
$$
y_i = 10\sin(\pi x_i) + w (\mathbf{s}_i)+ \epsilon_i; \:\: \epsilon \sim N(\mathbf{0},\: \tau^2 \mathbf{I}), \tau^2 = 0.1; \:\:\: w \sim \textit{exponential GP};\: \sigma^2 = 10; \phi = 1.
$$

Here, the mean function is $E(Y) = 10\sin(\pi X)$; $w$ accounts for the spatial correlation, which is generated as a exponential Gaussian process with spatial variance $\sigma^2 = 10$ and spatial correlation decay $\phi = 1$; and $\epsilon$ is the i.i.d random noise with variance $\tau^2 = 0.1$, which is also called the nugget in spatial literature.

For illustration purposes, we simulate with $n = 200$:

```{r simulation, message=FALSE, warning=FALSE}
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(5)
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))
set.seed(2)
x <- as.matrix(runif(n),n,1)
sigma.sq = 10
phi = 1
tau.sq = 0.1
D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))
```



##  Model fitting

In the package `RandomForestsGLS`, the working precision matrix used in the GLS-loss are NNGP approximations of precision matrices corresponding to Matérn covariance function. 

In order to fit the model, the code requires: 

+ Coordinates (`coords`): an $n \times 2$ matrix of 2-dimensional locations.
+ Response (`y`):  an $n$ length vector of response at the observed coordinates.
+ Covariates (`X`): an $n \times p$ matrix of the covariates in the observation coordinates.
+ Covariates for estimation (`Xtest`): an $ntest \times p$ matrix of the covariates where we want to estimate the function. Must have identical variables as that of `X`. Default is `X`.
+ Minimum size of leaf nodes (`nthsize`): We recommend not setting this value too small, as that will lead to very deep trees that takes a lot of time to be built and can produce unstable estimates. Default value is 20.
+ The parameters corresponding to the covariance function (detailed afterwards).

For the details on choice of other parameters, please refer to the help file of the code `RFGLS_estimate_spatial`, which can be accessed with `?RFGLS_estimate_spatial`. 

### Known Covariance Parameters

If the covariance parameters are known, we set `param_estimate = FALSE` (default value); the code additionally requires the following:

+ Covariance Model (`cov.model`): Supported keywords are: "exponential", "matern", "spherical", and "gaussian" for exponential, Matérn, spherical and Gaussian covariance function respectively. Default value is "exponential".
+ $\sigma^2$ (`sigma.sq`): The spatial variance. Default value is 1.
+ $\tau^2$ (`tau.sq`): The nugget. Default value is 0.01.
+ $\phi$ (`phi`): The spatial correlation decay parameter. Default value is 5.
+ $\nu$ (`nu`): The smoothing parameter corresponding to the Matérn covariance function. Default value is 0.5. 

We can fit the model as follows:

```{r model_fit, message=FALSE, warning=FALSE}
set.seed(1)
est_known <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential",
                                    nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq, 
                                    phi = phi)
```
The estimate of the function at the covariates `Xtest` is given in `estimation_reult$predicted`. For interpretation of the rest of the outputs, please see the help file of the code `RFGLS_estimate_spatial`. Using covariance models other than exponential model are in beta testing stage.

### Unknown Covariance Parameters

If the covariance parameters are not known we set `param_estimate = TRUE`; the code additionally requires the covariance model (`cov.model`) to be used for parameter estimation prior to RF-GLS fitting. We fit the model with unknown covariance parameters as follows.  

```{r model_fit_unknown, message=FALSE, warning=FALSE}
set.seed(1)
est_unknown <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential", 
                                      nthsize = 20, param_estimate = TRUE)
```



## Prediction of mean function

Given a fitted model using `RFGLS_estimate_spatial`, we can estimate the mean function at new covariate values as follows:

```{r model_estimate, message=FALSE, warning=FALSE}
Xtest <- matrix(seq(0,1, by = 1/10000), 10001, 1)
RFGLS_predict_known <- RFGLS_predict(est_known, Xtest)
```

### Performance comparison

We obtain the Mean Integrated Squared Error (MISE) of the estimate $\hat{m}$ from RF-GLS on [0,1] and compare it with that corresponding to the classical Random Forest (RF) obtained using package `randomForest` (with similar minimum nodesize, `nodesize = 20`, as default `nodesize` performs worse). We see that our method has a significantly smaller MISE. Additionally, we show that the MISE obtained with unknown parameters in RF-GLS is comparable to that of the MISE obtained with known covariance parameters.

```{r comparison, message=FALSE, warning=FALSE}
library(randomForest)
set.seed(1)
RF_est <- randomForest(x, y, nodesize = 20)

RF_predict <- predict(RF_est, Xtest)
#RF MISE
mean((RF_predict - 10*sin(pi * Xtest))^2)

#RF-GLS MISE
mean((RFGLS_predict_known$predicted - 10*sin(pi * Xtest))^2)

RFGLS_predict_unknown <- RFGLS_predict(est_unknown, Xtest)
#RF-GLS unknown MISE
mean((RFGLS_predict_unknown$predicted - 10*sin(pi * Xtest))^2)
```

We plot the true $m(x) = 10sin(\pi x)$ along with the loess-smoothed version of estimated $\hat{m}(.)$ obtained from RF-GLS and RF where we show that RF-GLS estimate approximates $m(x)$ better than that corresponding to RF.

```{r plot_comparison, message=FALSE, warning=FALSE}
rfgls_loess_10 <- loess(RFGLS_predict_known$predicted ~ c(1:length(Xtest)), span=0.1)
rfgls_smoothed10 <- predict(rfgls_loess_10)

rf_loess_10 <- loess(RF_predict ~ c(1:length(RF_predict)), span=0.1)
rf_smoothed10 <- predict(rf_loess_10)

xval <- c(10*sin(pi * Xtest), rf_smoothed10, rfgls_smoothed10)
xval_tag <- c(rep("Truth", length(10*sin(pi * Xtest))), rep("RF", length(rf_smoothed10)), 
              rep("RF-GLS",length(rfgls_smoothed10)))
plot_data <- as.data.frame(xval)
plot_data$Methods <- xval_tag
coval <- c(rep(seq(0,1, by = 1/10000), 3))
plot_data$Covariate <- coval

library(ggplot2)
ggplot(plot_data, aes(x=Covariate, y=xval, color=Methods)) +
geom_point() + labs( x = "x") + labs( y = "f(x)")
```


## Prediction of spatial response

Given a fitted model using `RFGLS_estimate_spatial`, we can predict the spatial response/outcome at new locations provided the covariates at that location. This approach performs kriging at a new location using the mean function estimates at the corresponding covariate values. Here we partition the simulated data into training and test sets in 4:1 ratio. Next we perform prediction on the test set using a model fitted on the training set.

```{r prediction_spatial, message=FALSE, warning=FALSE}
est_known_short <- RFGLS_estimate_spatial(coords[1:160,], y[1:160], 
                   matrix(x[1:160,],160,1), ntree = 50,  cov.model = "exponential", 
                   nthsize = 20, param_estimate = TRUE)
RFGLS_predict_spatial <- RFGLS_predict_spatial(est_known_short, coords[161:200,], 
                                               matrix(x[161:200,],40,1))
pred_mat <- as.data.frame(cbind(RFGLS_predict_spatial$prediction, y[161:200]))
colnames(pred_mat) <- c("Predicted", "Observed")
ggplot(pred_mat, aes(x=Observed, y=Predicted)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  ylim(0, 16) + xlim(0, 16)
```



## Misspecification in covariance model

The following example considers a setting when the parameters are estimated from a misspecified covariance model. We simulate the spatial correlation from a Matérn covariance function with smoothing parameter $\nu = 1.5$. While fitting the RF-GLS, we estimate the covariance parameters using an exponential covariance model ($\nu = 0.5$) and show that the obtained MISE can compare favorably to that of classical RF. 

```{r misspec_spatial, message=FALSE, warning=FALSE}
#Data simulation from matern with nu = 1.5
nu = 3/2
R1 <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
diag(R1) <- 1
set.seed(2)
w <- rmvn(1, rep(0,n), sigma.sq*R1)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

#RF-GLS with exponential covariance
set.seed(3)
est_misspec <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential", 
                                      nthsize = 20, param_estimate = TRUE)
RFGLS_predict_misspec <- RFGLS_predict(est_misspec, Xtest)

#RF
set.seed(4)
RF_est <- randomForest(x, y,  nodesize = 20)
RF_predict <- predict(RF_est, Xtest)

#RF-GLS MISE
mean((RFGLS_predict_misspec$predicted - 10*sin(pi * Xtest))^2)
#RF MISE
mean((RF_predict - 10*sin(pi * Xtest))^2)
```

# 2. Autoregressive Time Series Data

RF-GLS can also be used for function estimation in a time series setting under autoregressive errors. We consider time series data with errors from an AR(q) process as follows:

$$
y_t = m(\mathbf{x}_t) + e_t;\: e_t = \sum_{i = 1}^q\rho_i e_{t-i} + \eta_t
$$
where, $y_i, \mathbf{x}_i$ denotes the response and the covariate corresponding to the $t^{th}$ time point, $e_t$ is an AR(q) pprocess, $\eta_t$ denotes the i.i.d. white noise and $(\rho_1, \cdots, \rho_q)$ are the model parameters that captures the dependence of $e_t$ on $(e_{t - 1}, \cdots, e_{t-q})$.


In the AR time series scenario, the package `RandomForestsGLS` allows for fitting $m(.)$ using RF-GLS. RF-GLS exploits the sparsity of the closed form precision matrix of the AR process for model fitting and prediction of mean function $m(.)$. 

## Illustration

Here, we simulate from the AR(1) process as follows:
$$
y = 10\sin(\pi x) + \mathbf{e}; e_t = \rho e_{t-1} + \eta_t; \eta_t \sim N(0,\sigma^2); e_1 = \eta_1; \: \rho = 0.9; \sigma^2 = 10.
$$

Here, $E(Y) = 10\sin(\pi X)$; $\mathbf{e}$ which is an AR(1) process, accounts for the temporal correlation, $\sigma^2$ denotes the variance of white noise part of the AR(1) process and $\rho$ captures the degree of dependence of $e_t$ on $e_{t-1}$. 

For illustration purposes, we simulate with $n = 200$:

```{r simulation_temporal}
rho <- 0.9
set.seed(1)
b <- rho
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(1,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)
y <- c(eps + 10*sin(pi * x))
```


##  Model fitting

In case of time series data, the code requires: 

+ Response (`y`):  an $n$ length vector of response at the observed time points.
+ Covariates (`X`): an $n \times p$ matrix of the covariates in the observation time points.
+ Covariates for estimation (`Xtest`): an $ntest \times p$ matrix of the covariates where we want to estimate the function. Must have identical variables as that of `X`. Default is `X`.
+ Minimum size of leaf nodes (`nthsize`): We recommend not setting this value too small, as that will lead to very deep trees that takes a lot of time to be built and can produce unstable estimates. Default value is 20.
+ The parameters corresponding to the AR process (detailed afterwards).

For the details on choice of other parameters, please refer to the help file of the code `RFGLS_estimate_timeseries`, which can be accessed with `?RFGLS_estimate_timeseries`. 

### Known AR process Parameters

If the AR process parameters are known we set `param_estimate = FALSE` (default value); the code additionally requires `lag_params` = $c(\rho_1, \cdots, \rho_q)$.

We can fit the model as follows:

```{r model_fit_temporal, message=FALSE, warning=FALSE}
set.seed(1)
est_temp_known <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = rho, nthsize = 20)
```


### Unknown AR process Parameters

If the AR process parameters are not known, we set `param_estimate = TRUE`; the code requires the order of the AR process, which is obtained from the length of the `lag_params` input vector. Hence if we want to estimate the parameters from a AR(q) process, `lag_params` should be any vector of length `q`. Here we fit the model with `q = 1`

```{r model_fit_temporal_unknown, message=FALSE, warning=FALSE}
set.seed(1)
est_temp_unknown <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = rho, 
                                              nthsize = 20, param_estimate = TRUE)
```

## Prediction of mean function

This part of time series data analysis is identical to that corresponding to the spatial data.

```{r model_estimate_temporal, message=FALSE, warning=FALSE}
Xtest <- matrix(seq(0,1, by = 1/10000), 10001, 1)
RFGLS_predict_temp_known <- RFGLS_predict(est_temp_known, Xtest)
```

Here also, similar to the spatial data scenario, RF-GLS outperforms classical RF in terms of MISE both with true and estimated AR process parameters.

```{r comparison_temp, message=FALSE, warning=FALSE}
library(randomForest)
set.seed(1)
RF_est_temp <- randomForest(x, y, nodesize = 20)

RF_predict_temp <- predict(RF_est_temp, Xtest)
#RF MISE
mean((RF_predict_temp - 10*sin(pi * Xtest))^2)

#RF-GLS MISE
mean((RFGLS_predict_temp_known$predicted - 10*sin(pi * Xtest))^2)

RFGLS_predict_temp_unknown <- RFGLS_predict(est_temp_unknown, Xtest)
#RF-GLS unknown MISE
mean((RFGLS_predict_temp_unknown$predicted - 10*sin(pi * Xtest))^2)
```

## Misspecification in AR process order

We consider a scenario where the order of autoregression used for RF-GLS model fitting is mis-specified. We simulate the AR errors from an AR(2) process and fit RF-GLS with an AR(1) process.

```{r misspec_temporal, message=FALSE, warning=FALSE}
#Simulation from AR(2) process
rho1 <- 0.7
rho2 <- 0.2
set.seed(2)
b <- c(rho1, rho2)
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(2,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)
y <- c(eps + 10*sin(pi * x))

#RF-GLS with AR(1)
set.seed(3)
est_misspec_temp <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = 0, 
                                              nthsize = 20, param_estimate = TRUE)
RFGLS_predict_misspec_temp <- RFGLS_predict(est_misspec_temp, Xtest)

#RF
set.seed(4)
RF_est_temp <- randomForest(x, y,  nodesize = 20)
RF_predict_temp <- predict(RF_est_temp, Xtest)

#RF-GLS MISE
mean((RFGLS_predict_misspec_temp$predicted - 10*sin(pi * Xtest))^2)
#RF MISE
mean((RF_predict_temp - 10*sin(pi * Xtest))^2)
```


# Parallelization

For `RFGLS_estimate_spatial`, `RFGLS_estimate_timeseries`, `RFGLS_predict` and `RFGLS_predict_spatial` one can also take the advantage of parallelization, contingent upon the availability of multiple cores. The component `h` in all the functions determines the number of cores to be used. Here we demonstrate an example with `h = 2`.

```{r parallel_spatial, message=FALSE, warning=FALSE}
#simulation from exponential distribution
set.seed(5)
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))
set.seed(2)
x <- as.matrix(runif(n),n,1)
sigma.sq = 10
phi = 1
tau.sq = 0.1
nu = 0.5
D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

#RF-GLS model fitting and prediction with parallel computation
set.seed(1)
est_known_pl <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential",
                                    nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq, 
                                    phi = phi, h = 2)
RFGLS_predict_known_pl <- RFGLS_predict(est_known_pl, Xtest, h = 2)

#MISE from single core
mean((RFGLS_predict_known$predicted - 10*sin(pi * Xtest))^2)
#MISE from parallel computation
mean((RFGLS_predict_known_pl$predicted - 10*sin(pi * Xtest))^2)
```

For `RFGLS_estimate_spatial` with very small dataset (`n`) and small number of trees (`ntree`), communication overhead between the nodes for parallelization outweighs the benefits of the parallel computing hence it is recommended to parallelize only for moderately large `n` and/or `ntree`. It is strongly recommended that the max value of `h` is kept strictly less than the number of total cores available. Parallelization for `RFGLS_estimate_timeseries` can be addressed identically. For `RFGLS_predict` and `RFGLS_predict_spatial`, even for large dataset, single core performance is very fast, hence unless `ntest` and `ntree` are very high, we do not recommend using parallelization for `RFGLS_predict` and `RFGLS_predict_spatial`.
\name{RFGLS_estimate_spatial}
\alias{RFGLS_estimate_spatial}
\title{Function for estimation in spatial data with RF-GLS}

\description{
  The function \code{RFGLS_estimate_spatial} fits univariate non-linear spatial regression models for
  spatial data using RF-GLS in Saha et al. 2020. \code{RFGLS_estimate_spatial} uses the sparse Cholesky representation
  of Vecchia’s likelihood (Vecchia, 1988) developed in Datta et al., 2016 and Saha & Datta, 2018.
  The fitted Random Forest (RF) model is used later for prediction via the \code{RFGLS_predict} and \code{RFGLS_predict_spatial}.

  Some code blocks are borrowed from the R packages: spNNGP:
  Spatial Regression Models for Large Datasets using Nearest
  Neighbor Gaussian Process \cr
  https://CRAN.R-project.org/package=spNNGP and randomForest:
  Breiman and Cutler's Random
  Forests for Classification and Regression \cr
  https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_estimate_spatial(coords, y, X, Xtest = NULL,
                       nrnodes = NULL, nthsize = 20,
                       mtry = 1, pinv_choice = 1,
                       n_omp = 1, ntree = 50, h = 1,
                       sigma.sq = 1, tau.sq = 0.1,
                       phi = 5, nu = 0.5,
                       n.neighbors = 15,
                       cov.model = "exponential",
                       search.type = "tree",
                       param_estimate = FALSE,
                       verbose = FALSE)
}

\arguments{

  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation
  coordinates in \eqn{R^2}{R^2} (e.g., easting and northing). }

  \item{y}{an \eqn{n}{n} length vector of response at the observed coordinates. }

  \item{X}{an \eqn{n \times p}{n x p} matrix of the covariates in the observation coordinates. }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction locations. Its Structure should be
               identical (including intercept) with that of covariates provided for estimation purpose in \code{X}.
               If \code{NULL}, will use \code{X} as \code{Xtest}. Default value is \code{NULL}. }

  \item{nrnodes}{the maximum number of nodes a tree can have. Default choice leads to the deepest tree contigent on \code{nthsize}. For significantly large \eqn{n},
                   one needs to bound it for growing shallow trees which trades off efficiency for computation time. }

  \item{nthsize}{minimum size of leaf nodes. We recommend not setting this value too small, as that will lead to very deep trees
                 that takes a lot of time to be built and can produce unstable estimaes. Default value is 20. }

  \item{mtry}{number of variables randomly sampled at each partition as a candidate split direction. We recommend using
              the value p/3 where p is the number of variables in \code{X}. Default value is 1. }

  \item{pinv_choice}{dictates the choice of method for obtaining the pseudoinverse involved in the cost function and node
                     representative evaluation. if pinv_choice = 0, SVD is used (slower but more stable), if pinv_choice = 1,
                     orthogonal decomposition (faster, may produce unstable results if \code{nthsize} is too low) is used.
                     Default value is 1. }

  \item{n_omp}{number of threads to be used, value can be more than 1 if source code is compiled with OpenMP support.
               Default is 1. }

  \item{ntree}{number of trees to be grown. This value should not be too small. Default value is 50. }

  \item{h}{number of core to be used in parallel computing setup for bootstrap samples. If h = 1, there is no parallelization.
           Default value is 1. }

  \item{sigma.sq}{value of sigma square. Default value is 1. }

  \item{tau.sq}{value of tau square. Default value is 0.1. }

  \item{phi}{value of phi. Default value is 5. }

  \item{nu}{value of nu, only required for matern covariance model. Default value is 0.5. }

  \item{n.neighbors}{number of neighbors used in the NNGP. Default value is 15. }

  \item{cov.model}{keyword that specifies the covariance function to be used in modelling the spatial dependence structure
                   among the observations. Supported keywords are: "exponential", "matern", "spherical", and "gaussian"
                   for exponential, matern, spherical and gaussian covariance function respectively. Default value is "exponential". }

  \item{search.type}{keyword that specifies type of nearest neighbor search algorithm to be used. Supported keywords are:
                     "tree" and "brute". Both of them provide the same result, though "tree" should be faster.
                     Default value is "tree". }

  \item{param_estimate}{if \code{TRUE}, using the residuals obtained from fitting a classical RF with default options and \code{nodesize = nthsize},
                        will estimate the coefficeints corresponding to \code{cov.model} from \code{BRISC_estimate} with the deafult options.
                        Default value is \code{FALSE}. }

  \item{verbose}{if \code{TRUE}, model specifications along with information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is \code{FALSE}. }
}

\value{
  A list comprising:

  \item{P_matrix}{an \eqn{n \times ntree}{n x ntree} matrix of zero indexed resamples. t-th column denote the
                  \eqn{n} resamples used in the t-th tree. }

  \item{predicted_matrix}{an \eqn{ntest \times ntree}{ntest x ntree} matrix of predictions. t-th column denote the
                          predictions at \eqn{ntest} datapoints obtained from the t-th tree. }

  \item{predicted}{preducted values at the \eqn{ntest} prediction points. Average (\code{rowMeans}) of the treewise predctions
                   in \code{predicted_matrix}, }

  \item{X}{the matrix \code{X}. }

  \item{y}{the vector \code{y}. }

  \item{RFGLS_Object}{object required for prediction. }

}

\references{
  Saha, A., Basu, S., & Datta, A. (2020). Random Forests for dependent data. arXiv preprint arXiv:2007.15421.

  Saha, A., & Datta, A. (2018). BRISC: bootstrap for rapid inference on spatial
  covariances. Stat, e184, DOI: 10.1002/sta4.184.

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. Journal of the American Statistical
  Association, 111:800-812.

  Andrew Finley, Abhirup Datta and Sudipto Banerjee (2017). spNNGP: Spatial Regression Models for Large
  Datasets using Nearest Neighbor Gaussian Processes. R package version 0.1.1.
  https://CRAN.R-project.org/package=spNNGP

  Andy Liaw, and Matthew Wiener (2015). randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression.
  R package version 4.6-14. \cr https://CRAN.R-project.org/package=randomForest

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Sumanta Basu \email{sumbose@cornell.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}
}

\examples{

rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(1)
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))

set.seed(2)
x <- as.matrix(rnorm(n),n,1)

sigma.sq = 1
phi = 5
tau.sq = 0.1

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)

y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

estimation_result <- RFGLS_estimate_spatial(coords, y, x, ntree = 10)

}

\keyword{model}

\name{RFGLS_estimate_timeseries}
\alias{RFGLS_estimate_timeseries}
\title{Function for estimation in time-series data with RF-GLS}

\description{
  The function \code{RFGLS_estimate_spatial} fits univariate non-linear regression models for
  time-series data using a RF-GLS in Saha et al. 2020. \code{RFGLS_estimate_spatial} uses the sparse Cholesky representation
  corresponsinding to \code{AR(q)} process. The fitted Random Forest (RF) model is used later for
  prediction via the \code{RFGLS-predict}.

  Some code blocks are borrowed from the R packages: spNNGP:
  Spatial Regression Models for Large Datasets using Nearest Neighbor
  Gaussian Processes \cr https://CRAN.R-project.org/package=spNNGP and
  randomForest: Breiman and Cutler's Random Forests for Classification
  and Regression \cr https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_estimate_timeseries(y, X, Xtest = NULL, nrnodes = NULL,
                          nthsize = 20, mtry = 1,
                          pinv_choice = 1, n_omp = 1,
                          ntree = 50, h = 1, lag_params = 0.5,
                          variance = 1,
                          param_estimate = FALSE,
                          verbose = FALSE)
}

\arguments{

  \item{y}{an \eqn{n}{n} length vector of response at the observed time points. }

  \item{X}{an \eqn{n \times p}{n x p} matrix of the covariates in the observation time points. }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction. Its Structure should be
               identical (including intercept) with that of covariates provided for estimation purpose in \code{X}.
               If \code{NULL}, will use \code{X} as \code{Xtest}. Default value is \code{NULL}. }

  \item{nrnodes}{the maximum number of nodes a tree can have. Default choice leads to the deepest tree contigent on \code{nthsize}. For significantly large \eqn{n},
                 one needs to bound it for growing shallow trees which trades off efficiency for computation time. }

  \item{nthsize}{minimum size of leaf nodes. We recommend not setting this value too small, as that will lead to very deep trees
                 that takes a lot of time to be built and can produce unstable estimaes. Default value is 20. }

  \item{mtry}{number of variables randomly sampled at each partition as a candidate split direction. We recommend using
              the value p/3 where p is the number of variables in \code{X}. Default value is 1. }

  \item{pinv_choice}{dictates the choice of method for obtaining the pseudoinverse involved in the cost function and node
                     representative evaluation. if pinv_choice = 0, SVD is used (slower but more stable), if pinv_choice = 1,
                     orthogonal decomposition (faster, may produce unstable results if \code{nthsize} is too low) is used.
                     Default value is 1. }

  \item{n_omp}{number of threads to be used, value can be more than 1 if source code is compiled with OpenMP support.
               Default is 1. }

  \item{ntree}{number of trees to be grown. This value should not be too small. Default value is 50. }

  \item{h}{number of core to be used in parallel computing setup for bootstrap samples. If h = 1, there is no parallelization.
           Default value is 1. }

  \item{lag_params}{\eqn{q}{q} length vector of AR coefficients. If the parameters need to be estimated from AR(q) process, should be
                    any numeric vector of length q. For notations please see \code{arima}. Default value is 0.5. }

  \item{variance}{variance of the white noise in temporal error. The function estimate is not affected by this. Default value is 1. }

  \item{param_estimate}{if \code{TRUE}, using the residuals obtained from fitting a classical RF default options and \code{nodesize = nthsize},
                        will estimate the coefficeints corresponding to \eqn{AR(q)} from \code{arima} with the option, \code{include.mean = FALSE}.
                        Default value is \code{FALSE}. }

  \item{verbose}{if \code{TRUE}, model specifications along with information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is \code{FALSE}. }
}

\value{
  A list comprising:

  \item{P_matrix}{an \eqn{n \times ntree}{n x ntree} matrix of zero indexed resamples. t-th column denote the
                  \eqn{n} resamples used in the t-th tree. }

  \item{predicted_matrix}{an \eqn{ntest \times ntree}{ntest x ntree} matrix of predictions. t-th column denote the
                          predictions at \eqn{ntest} datapoints obtained from the t-th tree. }

  \item{predicted}{preducted values at the \eqn{ntest} prediction points. Average (\code{rowMeans}) of the treewise predctions
                   in \code{predicted_matrix}, }

  \item{X}{the matrix \code{X}. }

  \item{y}{the vector \code{y}. }

  \item{RFGLS_Object}{object required for prediction. }

}

\references{
  Saha, A., Basu, S., & Datta, A. (2020). Random Forests for dependent data. arXiv preprint arXiv:2007.15421.

  Saha, A., & Datta, A. (2018). BRISC: bootstrap for rapid inference on spatial
  covariances. Stat, e184, DOI: 10.1002/sta4.184.

  Andy Liaw, and Matthew Wiener (2015). randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression. R package version 4.6-14. \cr
  https://CRAN.R-project.org/package=randomForest

  Andrew Finley, Abhirup Datta and Sudipto Banerjee (2017). spNNGP: Spatial Regression Models for Large
  Datasets using Nearest Neighbor Gaussian Processes. R package version 0.1.1.
  https://CRAN.R-project.org/package=spNNGP

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Sumanta Basu \email{sumbose@cornell.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}
}

\examples{

set.seed(2)
n <- 200
x <- as.matrix(rnorm(n),n,1)

sigma.sq <- 1
rho <- 0.5

set.seed(3)
b <- rho
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(1,0,0), ar = b),
                n = n, rand.gen = rnorm, sd = s)

y <- eps + 10*sin(pi * x)

estimation_result <- RFGLS_estimate_timeseries(y, x, ntree = 10)

}

\keyword{model}

\name{RFGLS_predict_spatial}
\alias{RFGLS_predict_spatial}
\title{Spatial response prediction at new location with RF-GLS}

\description{
  The function \code{RFGLS_predict_spatial} performs fast prediction on a set of new locations by combining
  non-linear mean estimate from a fitted RF-GLS model in Saha et al. 2020 with spatial kriging estimate obtained by using Nearest Neighbor Gaussian Processes (NNGP) (Datta et al., 2016).

  Some code blocks are borrowed from the R packages: spNNGP:
  Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes \cr
  https://CRAN.R-project.org/package=spNNGP and randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression \cr https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_predict_spatial(RFGLS_out, coords.0, Xtest,
                      h = 1, verbose = FALSE)
}

\arguments{

  \item{RFGLS_out}{an object obtained from \code{RFGLS_estimate_spatial}. }

  \item{coords.0}{the spatial coordinates corresponding to prediction locations. \cr
                  Its structure should be same as that of coords
                  in \code{BRISC_estimation}. Default covariate value is a column of \eqn{1} to adjust for the mean (intercept). }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction. Its Structure should be
               identical (including intercept) with that of covariates provided for estimation purpose in \code{X}
               in \code{RFGLS_out}. }

  \item{h}{number of core to be used in parallel computing setup for bootstrap samples. If \code{h = 1}, there is no parallelization.
           Default value is 1. }

  \item{verbose}{if \code{TRUE}, model specifications along with information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is \code{FALSE}. }
}

\value{
  A list comprising:

  \item{prediction}{predicted spatial response corresponding to \code{Xtest} and \code{coords.0}. }

}

\references{
  Saha, A., Basu, S., & Datta, A. (2020). Random Forests for dependent data. arXiv preprint arXiv:2007.15421.

  Saha, A., & Datta, A. (2018). BRISC: bootstrap for rapid inference on spatial
  covariances. Stat, e184, DOI: 10.1002/sta4.184.

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. Journal of the American Statistical
  Association, 111:800-812.

  Andrew Finley, Abhirup Datta and Sudipto Banerjee (2017). spNNGP: Spatial Regression Models for Large
  Datasets using Nearest Neighbor Gaussian Processes. R package version 0.1.1.
  https://CRAN.R-project.org/package=spNNGP

  Andy Liaw, and Matthew Wiener (2015). randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression. R package version 4.6-14. \cr
  https://CRAN.R-project.org/package=randomForest

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Sumanta Basu \email{sumbose@cornell.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}
}

\examples{
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(1)
n <- 250
coords <- cbind(runif(n,0,1), runif(n,0,1))

set.seed(2)
x <- as.matrix(rnorm(n),n,1)

sigma.sq = 1
phi = 5
tau.sq = 0.1

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)

y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

estimation_result <- RFGLS_estimate_spatial(coords[1:200,], y[1:200],
                                 matrix(x[1:200,],200,1), ntree = 10)
prediction_result <- RFGLS_predict_spatial(estimation_result,
                           coords[201:250,], matrix(x[201:250,],50,1))

}

\keyword{model}

\name{RFGLS_predict}
\alias{RFGLS_predict}
\title{Prediction of mean function with RF-GLS}

\description{
  The function \code{RFGLS_predict} predicts the mean function at a given set of covariates.
  It uses a fitted RF-GLS model in Saha et al. 2020 to obtain the predictions.

  Some code blocks are borrowed from the R package: randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression \cr https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_predict(RFGLS_out, Xtest, h = 1, verbose = FALSE)
}

\arguments{

  \item{RFGLS_out}{an object obtained from \code{RFGLS_estimate_spatial}
                   or \cr \code{RFGLS_estimate_timeseries}. }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction. Its Structure should be
               identical (including intercept) with that of covariates provided for estimation purpose in \code{X}
               in \code{RFGLS_out}. }

  \item{h}{number of core to be used in parallel computing setup for bootstrap samples. If \code{h = 1}, there is no parallelization.
           Default value is 1. }

  \item{verbose}{if \code{TRUE}, model specifications along with information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is \code{FALSE}. }
}

\value{
  A list comprising:

  \item{predicted_matrix}{an \eqn{ntest \times ntree}{ntest x ntree} matrix of predictions. t-th column denote the
                          predictions at \eqn{ntest} datapoints obtained from the t-th tree. }

  \item{predicted}{preducted values at the \eqn{ntest} prediction points. Average (\code{rowMeans}) of the treewise predctions
                   in \code{predicted_matrix} }

}

\references{
  Saha, A., Basu, S., & Datta, A. (2020). Random Forests for dependent data. arXiv preprint arXiv:2007.15421.

  Andy Liaw, and Matthew Wiener (2015). randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression. R package version 4.6-14. \cr
  https://CRAN.R-project.org/package=randomForest

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Sumanta Basu \email{sumbose@cornell.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}
}

\examples{

rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(2)
n <- 200
x <- as.matrix(rnorm(n),n,1)

sigma.sq <- 1
rho <- 0.5

set.seed(3)
b <- rho
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(1,0,0), ar = b),
                n = n, rand.gen = rnorm, sd = s)

y <- eps + 10*sin(pi * x[,1])

estimation_result <- RFGLS_estimate_timeseries(y, x, ntree = 10)
Xtest <- matrix(seq(0,1, by = 1/1000), 1001, 1)
RFGLS_predict <- RFGLS_predict(estimation_result, Xtest)

}

\keyword{model}

