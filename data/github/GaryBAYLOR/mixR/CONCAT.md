<!-- badges: start -->
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/mixR)](https://CRAN.R-project.org/package=mixR)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/mixR)](https://cran.rstudio.com/web/packages/mixR/index.html)
[![CRAN Monthly Downloads](https://cranlogs.r-pkg.org/badges/mixR)](https://cran.r-project.org/package=mixR)
[![DOI](https://zenodo.org/badge/360674750.svg)](https://zenodo.org/badge/latestdoi/360674750)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04031/status.svg)](https://doi.org/10.21105/joss.04031)

<!-- badges: end -->

# mixR: An R package for finite mixture modeling for both raw and binned data

## Why `mixR`?
R programming language provides a rich collection of packages for building and analyzing finite mixture models which are widely used in unsupervised learning such as model-based clustering and density estimation. For example, 
- [`mclust`](https://cran.r-project.org/web/packages/mclust/index.html) can be used to build Gaussian mixture models with different covariance structures
- [`mixtools`](https://cran.r-project.org/web/packages/mixtools/index.html) implements parametric and non-parametric mixture models as well as mixtures of Gaussian regressions
- [`flexmix`](https://cran.r-project.org/web/packages/flexmix/index.html) provides a general framework for finite mixtures of regression models
- [`mixdist`](https://cran.r-project.org/web/packages/mixdist/index.html) fits mixture models for grouped and conditional data (also called binned data). 

To our knowledge, almost all R packages for finite mixture models are designed to use raw data as the modeling input except `mixdist`. However the popular model selection methods based on information criteria or bootstrapping likelihood ratio test ([McLachlan, 1987](https://doi.org/10.2307/2347790); [Feng & McCulloch, 1996](https://doi.org/10.1111/j.2517-6161.1996.tb02104.x); [Yu & Harvill, 2019](https://doi.org/10.1080/03610926.2018.1494838)) are not implemented in `mixdist`.

`mixR` is a package that aims to bridge this gap and to unify the interface for finite mixture modeling for both raw and binned data.


## Installation

For stable/pre-compiled(for Windows and OS X) version, please install from [CRAN](https://CRAN.R-project.org/package=mixR):

```r
install.packages('mixR')
```

To get the latest development version from Github:
```r
# install.packages('devtools')
devtools::install_github('garybaylor/mixR')
```

## Examples

* Fitting a normal mixture model
```r
library(mixR)

# generate data from a Normal mixture model
set.seed(102)
x1 = rmixnormal(1000, c(0.3, 0.7), c(-2, 3), c(2, 1))

# fit a Normal mixture model
mod1 = mixfit(x1, ncomp = 2)

# plot the fitted model
plot(mod1)

# fit a Normal mixture model (equal variance)
mod1_ev = mixfit(x1, ncomp = 2, ev = TRUE)
```

* Fitting a Weibull mixture model
```r
# generate data from a Weibull mixture model
x2 = rmixweibull(1000, c(0.4, 0.6), c(0.6, 1.3), c(0.1, 0.1))
mod2_weibull = mixfit(x2, family = 'weibull', ncomp = 2)
```
* Fitting a mixture model with binned data
```r
head(Stamp2)
##     lower  upper freq
## 1  0.0595 0.0605    1
## 5  0.0635 0.0645    2
## 6  0.0645 0.0655    1
## 7  0.0655 0.0665    1
## 9  0.0675 0.0685    1
## 10 0.0685 0.0695    7
mod_binned = mixfit(Stamp2, ncomp = 7, family = 'weibull')
plot(mod_binned)

# data binned from numeric data
x1_binned = bin(x1, seq(min(x1), max(x1), length = 30))
mod1_binned = mixfit(x1_binned, ncomp = 2)
```

* Mixture model selection by BIC
```r
# Selecting the best g for Normal mixture model
s_normal = select(x2, ncomp = 2:6)

# Selecting the best g for Weibull mixture model
s_weibull = select(x2, ncomp = 2:6, family = 'weibull')

plot(s_weibull)
plot(s_normal)
```

* Mixture model selection by bootstrap likelihood ratio test (LRT)
```r
b1 = bs.test(x1, ncomp = c(2, 3))
plot(b1, main = 'Bootstrap LRT for Normal Mixture Models (g = 2 vs g = 3)')
b1$pvalue

b2 = bs.test(x2, ncomp = c(2, 4))
plot(b2, main = 'Bootstrap LRT for Normal Mixture Models (g = 2 vs g = 4)')
b2$pvalue
```
For more examples please check the vignette [An Introduction to mixR](https://cran.r-project.org/web/packages/mixR/vignettes/intro-to-mixr.pdf).


## Contributor Code of Conduct
Everyone is welcome to contribute to the project through reporting issues, posting feature requests, updating documentation, submitting pull requests, or contact the project maintainer directly. To maintain a friendly atmosphere and to collaborate in a fun and productive way, we expect contributors to abide by the [Contributor Code of Conduct](https://www.contributor-covenant.org/version/1/0/0/code-of-conduct/).

## Citation

Yu, Y., (2022). mixR: An R package for Finite Mixture Modeling for Both Raw and Binned Data. Journal of Open Source Software, 7(69), 4031, https://doi.org/10.21105/joss.04031

BibTex information
```
@article{Yu2022,
  doi = {10.21105/joss.04031},
  url = {https://doi.org/10.21105/joss.04031},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {69},
  pages = {4031},
  author = {Youjiao Yu},
  title = {mixR: An R package for Finite Mixture Modeling for Both Raw and Binned Data},
  journal = {Journal of Open Source Software}
}
```

## 0.2.0

- Updated the function plot.mixfitEM()
- Updated the function density.mixfitEM()
- Corrected the function to_k_lambda_weibull()
- Solved the NA issue in the output of mixfit() (NaN still happens when EM algorithm fails to converge)
- Added a vignette
- Added tests

## 0.1.1

- Resolved the warnings issues found in [check results](https://www.r-project.org/nosvn/R.check/r-release-windows-ix86+x86_64/mixR-00check.html).


## 0.1.0

- First version
## Test environments
* local OS X installation, R 4.0.5
* ubuntu 16.04 (on travis-ci), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note


---
title: 'mixR: An R package for Finite Mixture Modeling for Both Raw and Binned Data'
tags:
  - R
  - mixture models
  - EM algorithm
  - model selection
authors:
  - name: Youjiao Yu
    affiliation: 1
    orcid: 0000-0003-0519-9605
affiliations:
  - name: Department of Statistical Science, Baylor University
    index: 1
date: "28 December 2021"
bibliography: paper.bib
---

# Statement of need

R [@R] provides a rich collection of packages for building and analyzing finite mixture models, which are widely used in unsupervised learning, such as model-based clustering and density estimation. For example, `mclust` [@mclust] can be used to build Gaussian mixture models with different covariance structures, `mixtools` [@mixtools] implements parametric and non-parametric mixture models as well as mixtures of Gaussian regressions, `flexmix` [@flexmix] provides a general framework for finite mixtures of regression models, `mixdist` [@mixdist] fits mixture models for grouped and conditional data (also called binned data). To our knowledge, almost all R packages for finite mixture models are designed to use raw data as the modeling input except `mixdist`. However, the popular model selection methods based on information criteria or bootstrapping likelihood ratio test (bLRT) [@mclachlan1987; @feng1996; @yu2019] are not implemented in `mixdist`. To bridge this gap and to unify the interface for finite mixture modeling for both raw and binned data, we implement `mixR` package that provides the following primary features.

-   `mixfit()` performs maximum likelihood estimation (MLE) for finite mixture models for Gaussian, Weibull, Gamma, and Log-normal distributions via EM algorithm [@dempster1977]. The model fitting is accelerated via package `Rcpp` [@rcpp].

-   `select()` selects the best model from a series of mixture models with a different number of mixture components by using Bayesian Information Criterion (BIC).

-   `bs.test()` performs bLRT for two mixture models from the same distribution family but with a different number of components.

`mixR` also contains the following additional features.

-   Visualization of the fitted mixture models using `ggplot2` [@ggplot2].
-   Functions to generate random data from mixture models.
-   Functions to convert parameters of Weibull and Gamma mixture models between shape-scale representation used in probability density functions and mean-variance representation which is more intuitive for people to understand the distribution.

# Examples

We demonstrate how to use `mixR` for fitting finite mixture models and selecting mixture models using BIC and bLRT.

## Model fitting

We fit the following four mixture models to a data set that consists of 1000 random data points generated from a Weibull mixture model with two components.

-   Gaussian mixture with two components (`mod1`)
-   Gaussian mixture with two components to the binned data (`mod2`)
-   Gaussian mixture with three components (`mod3`)
-   Weibull mixture with two components (`mod4`)

The fitted coefficients in `mod1` and `mod2` and the top two plots in Figure \ref{fig:plot1} show that binning does not cause much information loss, and we get similar fitted results using either raw data or binned data. This is usually the case when we have at least moderate data size, and the underlying mixture model is not too complex (e.g., too many mixture components). A benefit of binning is that it reduces the computation burden significantly for large data, especially when conducting bLRT, which is computationally intensive. From Figure \ref{fig:plot1} we also observe that Gaussian mixture models can provide a good fit for non-Gaussian data though the number of mixture components tends to be overestimated because more Gaussian components are needed to model the asymmetry and long tails that usually exist in non-Gaussian data.

```{r}
library(mixR)

set.seed(101)
x <- rmixweibull(1000, c(0.4, 0.6), c(0.6, 1.3), c(0.1, 0.1))
x_binned <- bin(x, brks = seq(min(x), max(x), length = 30))

mod1 <- mixfit(x, ncomp = 2)
mod2 <- mixfit(x_binned, ncomp = 2)
mod3 <- mixfit(x, ncomp = 3)
mod4 <- mixfit(x, ncomp = 2, family = 'weibull')

mod1
## Normal mixture model with 2 components
##        comp1     comp2
## pi 0.4210604 0.5789396
## mu 0.6014690 1.3084871
## sd 0.1092375 0.0932826
## 
## EM iterations: 5 AIC: -406.65 BIC: -382.11 log-likelihood: 208.32

mod2
## Normal mixture model with 2 components
##        comp1     comp2
## pi 0.4213019 0.5786981
## mu 0.6018737 1.3091224
## sd 0.1084973 0.0916267
## 
## EM iterations: 9 AIC: 5813.09 BIC: 5837.63 log-likelihood: -2901.54

p1 <- plot(mod1, title = 'Gaussian Mixture (2 components)')
p2 <- plot(mod2, title = 'Gaussian Mixture (binned data 2 components)')
p3 <- plot(mod3, title = 'Gaussian Mixture (3 components)')
p4 <- plot(mod4, title = 'Weibull Mixture (2 components)')
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
```

![(top left) the fitted Gaussian mixture with two components; (top right) the fitted Gaussian mixture with two components to the binned data; (bottom left) the fitted Gaussian mixture with three components; (bottom right) the fitted Weibull mixture with two components \label{fig:plot1}](plot1.png)

## Model selection

Figure \ref{fig:plot2} shows that the best Gaussian mixture model selected by BIC has three components and unequal variances for each component, while the best Weibull mixture model has two components. The bLRT with $H_0: g=2$ versus $H_a: g=3$ for Gaussian mixture models (using the default 100 bootstrap iterations) returns a p-value of zero, showing that Gaussian mixture with three components is significantly better than that with two components. Similarly, the same test for Weibull mixture models returns an insignificant p-value of 0.82, indicating that the Weibull mixture with three components is no better than it with two components.

```{r}
b1 <- select(x, ncomp = 2:4)
b2 <- select(x, ncomp = 2:4, family = 'weibull')
b3 <- bs.test(x, ncomp = c(2, 3))
b4 <- bs.test(x, ncomp = c(2, 3), family = 'weibull')

b3$pvalue
## [1] 0

b4$pvalue
## [1] 0.82

par(mfrow = c(2, 2))
plot(b1)
plot(b2, main = "Weibull Mixture Model Selection by BIC")
plot(b3, main = "Bootstrap LRT for Gaussian Mixture Models\n 
     (g = 2 vs. g = 3)", xlab = 'Bootstrap Test Statistics')
plot(b4, main = "Bootstrap LRT for Weibull Mixture Models\n
     (g = 2 vs. g = 3)", xlab = 'Bootstrap Test Statistics')
```

![(top left) Gaussian mixture model selection using BIC (UV stands for unequal variances for each mixture components and EV stands for equal variance); (top right) Weibull mixture model selection using BIC; (bottom left) bLRT with $H_0: g=2$ versus $H_a: g=3$ for Gaussian mixture models; (bottom right) bLRT with $H_0: g=2$ versus $H_a: g=3$ for Weibull mixture models \label{fig:plot2}](plot2.png)

# Summary

`mixR` unifies the interface for fitting and comparing finite mixture models for both raw data and binned data for distributions including Gaussian, Weibull, Gamma, and Log-normal. The package also provides features for generating random data from mixture models, conversion of parameters for Weibull and Gamma models, and model visualization in `ggplot2`. The heavy computation in `mixR` is completed in C++ using `Rcpp`.

`mixR` is actively used by researchers and practitioners in various fields [@jung2020; @sylvestre2020; @ogana2020; @de2021; @buckland2021; @buchel2021; @yang2021; @yang2021bio].

# References
---
#output:
#  pdf_document:
#    fig_caption: yes
#    toc: yes
#    toc_depth: 3
output: 
  rmarkdown::pdf_document:
    toc: true
    number_sections: true
    fig_caption: true
vignette: >
  %\VignetteIndexEntry{An Introduction to mixR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: master.bib
biblio-style: apsr
title: "An Introduction to mixR"
author: 
  - Youjiao Yu
abstract: "The package **mixR** performs maximum likelihood estimation (MLE) for finite mixture models for families including Normal, Weibull, Gamma and Lognormal via EM algorithm. It also conducts model selection by using Bayesian Information Criterion (BIC) or bootstrap likelihood ratio test (LRT). The data used for mixture model fitting can be raw data or binned data. The model fitting is accelerated by using R package **Rcpp**."
keywords: "mixture models, EM algorithm, model selection"
date: "`r format(Sys.time(), '%B %d, %Y')`"
geometry: margin=1in
fontfamily: mathpazo
# fontsize: 11pt
# spacing: double
endnote: no
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = "#>",
  #cache=TRUE,
  message=FALSE, 
  warning=FALSE,
  #fig.path='figs/',
  #fig.align="center",
  fig.width = 6,
  fig.height = 4)#,
  #cache.path = '_cache/')#,
                      #fig.process = function(x) {
                      #x2 = sub('-\\d+([.][a-z]+)$', '\\1', x)
                      #if (file.rename(x, x2)) x2 else x
                      #}
                      #)
```

# Background
## Mixture models
Finite mixture models can be represented by
$$
f(x; \Phi) = \sum_{j=1}^g \pi_j f_j(x; \theta_j)
$$
where $f(x; \Phi)$ is the probability density function (p.d.f.) or probability mass function
(p.m.f.) of the mixture model, $f_j(x; \theta_j)$ is the p.d.f. or p.m.f. of the $j$th
component of the mixture model, $\pi_j$ is the proportion of the $j$th component,
$\theta_j$ is the parameter of the $j$th component which can be a scalar or a vector,
$\Phi = (\pi_1, \theta_1, \dots, \pi_g, \theta_g)$ is a vector of all the parameters in the mixture model, and $g$ is the total number of components in the mixture model. The MLE of $\Phi$ can be obtained using the EM algorithm [@dempster1977].

## Mixture model selection by BIC
One critical problem for a mixture model is how to estimate $g$ when there is no such *a priori* knowledge. As EM algorithm doesn't estimate $g$ itself, a commonly used approach to estimate $g$ is to fit a series of mixture models with different values of $g$ and then select $g$ using information criteria such as Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), Deviance Information Criterion (DIC), or Integrated Complete-data Likelihood (ICL). Among all information criteria, BIC has shown to outperform other ones in model selection [@steele2010]. BIC is defined as
$$
BIC = k \log(n) - 2 \log(\hat L)
$$
in which $k$ is the total number of parameters in the mixture model, $n$ is the size of data, and $\hat L$ is the estimated maximum likelihood of the model. The model which has the lowest BIC is regarded as the optimal one.

## Mixture model selection by bootstrap LRT
A mixture model with $g = g_1$ components is a nested model of a mixture model with $g = g_2(g_1 < g_2)$ components, as the former model can be regarded as the later one with $\pi_j = 0$ for $g_2 - g_1$ components and $p_j > 0$ for all the remaining $g_1$ components. LRT is a common tool for assessing the goodness of fit of the nested model ($H_0: g = g_1$) compared to the full model ($H_a: g = g_2$). However the regularity condition of the LRT, which requires that the parameter space of the model in the null hypothesis $H_0$ should lie in the interior of the parameter space of the model in the alternative hypothesis $H_a$, doesn't hold for the mixture models [@feng1996], and therefore the test statistic of LRT, denoted as $w(x)$
doesn't follow a known Chi-square distribution under $H_0$. @mclachlan1987 proposed the idea of applying the method of bootstrapping [@efron1994] for approximating the distribution of $w(x)$. The general
 steps of bootstrap LRT are as follows.


1. For the given data $x$, estimate $\Phi$ under both $H_0$ and $H_a$ to get $\hat\Phi_0$ and $\hat\Phi_1$. Calculate the observed log-likelihood $\ell(x; \hat\Phi_0)$ and $\ell(x; \hat\Phi_1)$. The LRT statistic is defined as
$w_0 = -2(\ell(x; \hat\Phi_0) - \ell(x; \hat\Phi_1)).$
2. Generate random data of the same size as the original data $x$ from the model under the $H_0$ using estimated parameter $\hat\Phi_0$, then repeat step 1 using the simulated data. Repeat this process for $B$ times to get a vector of the simulated likelihood ratio test statistics $w_1^{(1)}, \dots, w_1^{(B)}$.
3. Calculate the empirical p-value as
$$
p = \frac{1}{B} \sum_{i=1}^B I(w_1^{(i)} > w_0)
$$
where $I(\cdot)$ is the indicator function.
 

## Fitting mixture models to the binned data
The binned data is present instead of the raw data in some situations, often for the reason of storage convenience or necessity. The binned data is recorded in the form of $(a_i, b_i, n_i)$ where $a_i$ is the left bin value of the $i^{th}$ bin, $b_i$ is the right bin value of the $i^{th}$ bin, and $n_i$ is the number of observations that fall in the $i^{th}$ bin for $i = 1, \dots, r$, where $r$ is the total number of bins.

The MLE for finite mixture models fitted to binned data can also be obtained via EM algorithm by introducing an additional latent variable $x$ that represents the unknown value of the raw data, besides the usually latent variable $z$ that represents the component an observation belongs to. To apply the EM algorithm we first write the complete-data
log-likelihood as
$$
Q(\Phi; \Phi^{(p)}) = \sum_{j = 1}^{g} \sum_{i = 1}^r n_i z^{(p)} [\log f(x^{(p)}; \theta_j)
 + \log \pi_j ]
$$
where $z^{(p)}$ is the expected value of $z$ given $\Phi^{(p)}$ and $x^{(p)}$, the estimated value of $\Phi$ and expected value of $x$ at $p^{th}$ iteration. The estimate of $\Phi$ can be updated alternatively via an E-step, in which we estimate $\Phi$ by maximizing $Q(\Phi; \Phi^{(p)})$, and an M-step, in which we compute $x^{(p)}$ and $z^{(p)}$, until the convergence of the EM algorithm. The M-step may not have a closed-form solution, e.g. in the Weibull mixture model or Gamma mixture model, which if is the case, an iterative approach like Newton's algorithm or bisection method may be used.

## Beyond normality
The Normal distribution is mostly used in a mixture model for continuous data, but there are also circumstances when other distributions fit the data better.  @mclachlan2004 explained that a limitation of the Normal distribution is that when the shapes of the components are skewed, there may not be a one-to-one correspondence between the number of components in the mixture model and that in the data. More than one Normal component is needed to model a skewed component, which may cause overestimation of $g$. For skewed or asymmetric components, other distributions such as Gamma, Log-normal or Weibull might provide better model fitting than the Normal distribution in a mixture model. As an example, @yu2019 demonstrated two examples where Weibull mixture models are preferred.


# mixR package

We present the functions in **mixR** package for (a) fitting finite mixture models for continuous data for families including Normal, Weibull, Gamma and Log-normal via EM algorithm; (b) selecting the optimal number of components for a mixture model using BIC or bootstrap LRT. We also discuss how to fit mixture models with binned data.

## Model fitting

The function `mixfit()` can be used to fit mixture models for four different families -- Normal, Weibull, Gamma, and Log-normal. For Normal distribution, the variances of each component are the same by setting `ev = TRUE`. 
```{r, cache=FALSE, fig.show="hold", out.width="50%", fig.cap="The fitted Normal mixture model with unequal variances (left) and equal variance(right)"}
# generate data from a Normal mixture model
library(mixR)
set.seed(102)
x1 = rmixnormal(1000, c(0.3, 0.7), c(-2, 3), c(2, 1))

# fit a Normal mixture model (unequal variances)
mod1 = mixfit(x1, ncomp = 2); mod1

# fit a Normal mixture model (equal variance)
mod1_ev = mixfit(x1, ncomp = 2, ev = TRUE); mod1_ev

plot(mod1, title = 'Normal Mixture Model (unequal variances)')
plot(mod1_ev, title = 'Normal Mixture Model (equal variance)')
```

The initial values for $\Phi$ are estimated by k-means or hierarchical clustering method if they are not provided by the users. In situations when EM algorithm is stuck in a local minimum and leads to unsatisfactory fitting results, which  happens more likely when the number of components $g$ and/or data size $n$ are large, initial values can be provided manually to get a better fitting.

To illustrate the idea that $g$ tends to be over-estimated when using a Normal mixture model to fit data with asymmetric or skewed components, we simulate data from a Weibull mixture model with $g = 2$, then fit both Normal and Weibull mixture models to the data. First we fit both models with $g = 2$. Weibull distribution provides better fitting than Normal by either visually checking the plots of the fitted results in Figure 2, or the fact that the log-likelihood of the fitted Weibull mixture model (244) is much higher than that of the fitted Normal mixture model (200). Figure 3 shows that the best value of $g$ for Weibull mixture model is two and for Normal mixture model is four, higher than the actual value of $g$. 

```{r, fig.show="hold", out.width="50%", cache=FALSE, eval=TRUE, fig.cap="The fitted Weibull mixture model (left) and Normal mixture model (right) to the same data"}
x2 = rmixweibull(1000, c(0.4, 0.6), c(0.6, 1.3), c(0.1, 0.1))
mod2_weibull = mixfit(x2, family = 'weibull', ncomp = 2); mod2_weibull
mod2_normal = mixfit(x2, ncomp = 2); mod2_normal
plot(mod2_weibull)
plot(mod2_normal)
```

## Model selection by BIC

The function `select()` is used to fit a series of finite mixture models with values of $g$ specified in `ncomp`, and then select the best $g$ by BIC. For Normal mixture models, both equal and unequal variances are considered. Figure 3 shows the value of BIC for Normal and Weibull mixture models with different $g$. For Weibull mixture models, BIC increases monotonically as $g$ increases from two to six, therefore the best value of $g$ is two. The best Normal mixture model is $g = 3$ with unequal variance as its BIC is the lowest. Figure 4 shows the fitted Weibull mixture models and Normal mixture model with the best values of $g$.

```{r, fig.show="hold", out.width="50%", cache=FALSE, eval=TRUE, fig.cap="The value of BIC for Weibull mixture models (left) and Normal mixture models (right) with different values of $g$."}
# Selecting the best g for Weibull mixture model
s_weibull = select(x2, ncomp = 2:6, family = 'weibull')

# Selecting the best g for Normal mixture model
s_normal = select(x2, ncomp = 2:6)
s_normal
plot(s_weibull)
plot(s_normal)
```

```{r, fig.show="hold", out.width="50%", cache=FALSE, eval=TRUE, fig.cap="The fitted Weibull mixture model with $g = 2$ (left) and the Normal mixture model with $g = 4$ and equal variance (right)"}
plot(mod2_weibull)
plot(mixfit(x2, ncomp = 3))
```


## Model selection by bootstrap LRT

The function `bs.test()` performs bootstrap LRT and returns the p-value as well as the test statistics $w_0$ and $w_1$. As an example, the data set `x1` above are generated from a Normal mixture model with $g = 2$. If we conduct bootstrap LRT for $g = 2$ against $g = 3$ for `x1` and set the number of bootstrap iterations `B=100`, we get p-value of 0.48, showing that the Normal mixture model with three components is not any better than the one with two components for data `x1`.

As another example, the data set `x2` above are generated from a Weibull mixture model with $g = 2$. We discussed previously that if we use Normal distribution to fit the mixture model, the best value of $g$ selected by BIC is three A bootstrap LRT of $g = 2$ against $g = 3$ returns zero p-value, indicating that if we fit a Normal mixture model to `x2`, $g = 3$ is a much better fit than $g= 2$, though visually the data shows two modes rather than three. Figure 5 shows the histogram of $w_1$ and the location of $w_0$ (red vertical line) for the above two examples.

```{r, fig.show="hold", out.width="50%", cache=FALSE, eval=TRUE, fig.cap="(left) The bootstrap LRT of $H_0: g = 2$ against $H_1: g = 3$ for fitting Normal mixture models for data `x1`; (right) The bootstrap LRT of $H_0: g = 2$ against $H_1: g = 3$ for fitting Normal mixture models for data `x2`. In each plot the histogram shows the distribution of $w_1$ and the red line shows the value of $w_0$."}
b1 = bs.test(x1, ncomp = c(2, 3))
plot(b1, main = 'Bootstrap LRT for Normal Mixture Models (g = 2 vs g = 3)')
b1$pvalue
b2 = bs.test(x2, ncomp = c(2, 3))
plot(b2, main = 'Bootstrap LRT for Normal Mixture Models (g = 2 vs g = 3)')
b2$pvalue
```


## Mixture model fitting with binned data

The function `mixfit()` can also fit mixture models with binned data, in the form of a three-column matrix each row of which represents a bin with the left bin value, the right bin value, and the total number of data points that fall in each bin (analogous to the data used to create a histogram). The package contains a data set `Stamp2` that is in the binned format. We fit a Weibull mixture model with three components to `Stamp2` as discussed in @yu2019. The `cassie` data set contained in `mixdist` package [@mixdist] is also in binned format and we fit a Gamma mixture model to `cassie` after some reformatting of the data. Figure 6 shows the fitted results of the two mixture models.
```{r, fig.show="hold", out.width="50%", cache=FALSE, eval=TRUE, fig.cap="(left) The Weibull mixture model with three components fitted with `Stamp2` data; (right) The Gamma mixture model with four components fitted with `cassie` data"}
head(Stamp2, 3)
mod_stamp = mixfit(Stamp2, ncomp = 3, family = 'weibull', pi = c(0.3, 0.3, 0.4), 
                   mu = c(0.07, 0.08, 0.1), sd = c(0.002, 0.002, 0.015))

library(mixdist)
data(cassie)
head(cassie, 3)
cassie_binned = data.frame(lower = cassie[1:(nrow(cassie)-1), 1], 
                           upper = cassie[2:nrow(cassie), 1],
                           freq = cassie[1:(nrow(cassie)-1), 2])
cassie_binned = as.matrix(cassie_binned[1:(nrow(cassie_binned)-1), ])

head(cassie_binned, 3)
mod_cassie = mixfit(cassie_binned, ncomp = 4, family = 'gamma')
plot(mod_stamp)
plot(mod_cassie)
```

The function `bin()` is used to create binned data from raw data, and the function `reinstate()` can simulate the raw data from binned data. Figure 7 shows the mixture models fitted with data binned from raw data `x1` and `x2`, with 30 bins for each data set.

```{r, fig.show="hold", out.width="50%", cache=FALSE, eval=TRUE, fig.cap="(left) The Normal mixture model fitted with binned data; (right) The Weibull mixture model fitted with binned data."}
x1_binned = bin(x1, seq(min(x1), max(x1), length = 30))
head(x1_binned, 3)

mod1_binned = mixfit(x1_binned, ncomp = 2)
plot(mod1_binned, xlab = 'x1_binned', 
     title = 'The Normal Mixture Model Fitted With Binned Data')
mod1_binned

x2_binned = bin(x2, seq(min(x2), max(x2), length = 30))
head(x2_binned, 3)

mod2_binned = mixfit(x2_binned, ncomp = 2, family = 'weibull')
plot(mod2_binned, xlab = 'x2_binned', 
     title = 'The Weibull Mixture Model Fitted With Binned Data')
mod2_binned
```

As binning can be considered a way to compress data, binned data can accelerate the fitting of mixture models, especially when the original data set is large. To illustrate, we simulate 100,000 data points from a Normal mixture model with five components, and bin the data with 100 bins. Normal mixture models are fitted on both the simulated raw data and binned data. The results show that model fitting takes 27 seconds on raw data, and less than one second on binned data. Another example shows that fitting a Weibull mixture model with data binned from a data set with one million observations takes just over two seconds \footnote{evaluated on iMac with processor: 3 GHz Quad-Core Intel Core i5, memory: 8 GB 2400 MHz DDR4}.
```{r, fig.show="hold", out.width="50%", cache=FALSE, eval=TRUE, fig.cap="Normal mixture models fitted to the raw data (left) and binned data (right)"}
# a function to generate parameters for a mixture model
generate_params = function(ncomp = 2) {
  pi = runif(ncomp)
  low = runif(1, 0, 0)
  upp = low + runif(1, 0, 10)
  mu = runif(ncomp, low, upp)
  sd = runif(ncomp, (max(mu) - min(mu))/ncomp/10, (max(mu) - min(mu))/ncomp/2)
  list(pi = pi / sum(pi), mu = sort(mu), sd = sd)
}

# simulate data from a Normal mixture model
set.seed(988)
n = 100000
ncomp = 5
params = generate_params(ncomp)
x_large = rmixnormal(n, pi = params$pi, mu = params$mu, sd = params$sd)

# fitting a Normal mixture model with raw data
t1 = Sys.time()
mod_large <- mixfit(x_large, ncomp = ncomp)
t2 = Sys.time()
t2 - t1

plot(mod_large, title = 'Normal Mixture Model Fitted With Raw Data')
mod_large

# fitting a Normal mixture model with binned data
t3 = Sys.time()
x_binned = bin(x_large, seq(min(x_large), max(x_large), length = 100))
mod_binned <- mixfit(x_binned, ncomp = ncomp)
t4 = Sys.time()
t4 - t3

plot(mod_binned, title = 'Normal Mixture Model Fitted With Binned Data')
mod_binned
```

# Citation
Run `citation(package = 'mixR')` to see how to cite package **mixR** in publications.


\newpage
# References
\setlength{\parindent}{-0.2in}
\setlength{\leftskip}{0.2in}
\setlength{\parskip}{8pt}
\vspace*{-0.2in}
\noindent




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibullpar.R
\name{to_k_lambda_weibull}
\alias{to_k_lambda_weibull}
\title{Parameter Conversion for Weibull Distribution}
\usage{
to_k_lambda_weibull(mu, sd)
}
\arguments{
\item{mu}{a numeric vector representing the means of Weibull distributions}

\item{sd}{a numeric vector representing the standard deviations of Weibull distributions.
\code{mu} and \code{sd} should have the same length.}
}
\value{
a list of two items
\item{k}{a vector of the shapes of Weibull distributions}
\item{lambda}{a vector of the scales of Weibull distributions}
}
\description{
The function \code{to_k_lambda_weibull} converts the mean and standard deviation to the shape and scale for
the Weibull distributions.
}
\details{
The purpose of this function is to convert the parameterization of Weibull distribution in the form of
mean and standard deviation to the form of shape and scale. It can be used for specifying the initial
values for the EM algorithm when the first-hand initial values are in the form of mean and standard
deviation from K-means clustering algorithm.
}
\examples{
to_k_lambda_weibull(2, 1)
to_k_lambda_weibull(c(2, 5), c(1, 0.7))

}
\seealso{
\code{\link{to_mu_sd_weibull}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.bootEM.R
\name{plot.bootEM}
\alias{plot.bootEM}
\title{Plot Bootstrap Likelihood Ratio Test}
\usage{
\method{plot}{bootEM}(x, ...)
}
\arguments{
\item{x}{an object of class \code{bootEM}, which is the output of the function \code{\link{bs.test}}.}

\item{...}{the other parameters passed to the function \code{\link{hist}}}
}
\description{
This function is the plot method for the class \code{bootEM}.
}
\details{
The histogram of the bootstrap LRT statistics \eqn{w_1} is plotted, with the
observed LRT statistic imposed in a red vertical line.
}
\examples{
## plotting the bootstrap LRT result
set.seed(100)
x <- rmixnormal(200, c(0.5, 0.5), c(2, 5), c(1, 0.7))
ret <- bs.test(x, ncomp = c(2, 3), B = 30)
plot(ret)

}
\seealso{
\code{\link{bs.test}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_Stamp.R
\docType{data}
\name{Stamp}
\alias{Stamp}
\title{1872 Hidalgo Stamp Data}
\format{
A vector with 485 measurements of the thickness (nm) of the stamps
}
\usage{
Stamp
}
\description{
A vector containing the 1872 Hidalgo stamp data
}
\references{
Izenman, A. J. and Sommer, C. J. Philatelic mixtures and multimodal densities.
\emph{Journal of the American Statistical association}, 83(404):941-953, 1988.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmixnormal.R
\name{rmixnormal}
\alias{rmixnormal}
\title{Generating Random Data From A Normal Mixture Model}
\usage{
rmixnormal(n, pi, mu, sd)
}
\arguments{
\item{n}{a positive integer specifying the number of observations we want to generate from the mixture model}

\item{pi}{a numeric vector for the proportion of each component}

\item{mu}{a numeric vector for the mean of each component}

\item{sd}{a numeric vector for the standard deviation of each component}
}
\value{
The function \code{rmixnormal} returns a numeric vector of random data from the specified normal mixture model.
}
\description{
The function \code{rmixnormal} generates random data from a normal mixture model.
}
\details{
The number of random data from each component \eqn{n_0} (a vector) is generated from a multinomial
distribution Multinom\eqn{(n, pi)}. Then the random data from each component is generated with
the sample sized specified in \eqn{n_0} and parameters of normal distributions specified in
\code{mu} and \code{sd}.
}
\examples{
x <- rmixnormal(1000, c(0.4, 0.6), c(2, 5), c(1, 0.5))
hist(x, breaks = 40)

}
\seealso{
\code{\link{rmixweibull}}, \code{\link{rmixgamma}}, \code{\link{rmixlnorm}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.mixfitEM}
\alias{plot.mixfitEM}
\title{Plotting the Fitted Mixture Models}
\usage{
\method{plot}{mixfitEM}(
  x,
  theme = NULL,
  add_hist = TRUE,
  add_poly = TRUE,
  add_legend = TRUE,
  smoothness = 512,
  trans = 0.5,
  cut = 3.8,
  xlab,
  ylab,
  title,
  breaks,
  plot.title = element_text(hjust = 0.5),
  axis.text.x = element_text(),
  axis.text.y = element_text(),
  axis.title.x = element_text(),
  axis.title.y = element_text(),
  legend.title = element_text(),
  legend.text = element_text(),
  legend.position = "right",
  legend.direction = ifelse(legend.position \%in\% c("top", "bottom"), "horizontal",
    "vertical"),
  ...
)
}
\arguments{
\item{x}{an object of class \code{mixfitEM}, an output from the function \code{\link{mixfit}}}

\item{theme}{a string the specifies the appearance of the plot, which is from the ggplot2 and could
be one of'gray', 'bw' (default), 'linedraw', 'light', 'dark', 'minimal', 'classic', or 'void'.}

\item{add_hist}{a logical value specifying whether a histogram of data should be plotted}

\item{add_poly}{a logical value specifying whether a polygon of each component should be plotted.}

\item{add_legend}{a logical value specifying whether the legend should be plotted.}

\item{smoothness}{a positive integer controlling the smoothness of the density curve in the plot.
The default value is 512 and increasing this value will produce smoother curve.}

\item{trans}{the transparency of the polygons if they are plotted (default 0.5)}

\item{cut}{the number of standard deviations from the center of each component we want to plot 
the density (default 3.8)}

\item{xlab}{the label for x axis}

\item{ylab}{the label for y axis}

\item{title}{the title of the plot}

\item{breaks}{the number of bins used for plotting the histogram}

\item{plot.title}{an object returned by element_text() to specify the appearance of the title}

\item{axis.text.x}{an object returned by element_text() to specify the appearance of the x axis}

\item{axis.text.y}{an object returned by element_text() to specify the appearance of the y axis}

\item{axis.title.x}{an object returned by element_text() to specify the appearance of the label along x axis}

\item{axis.title.y}{an object returned by element_text() to specify the appearance of the label along y axis}

\item{legend.title}{an object returned by element_text() to specify the appearance of the legend title}

\item{legend.text}{an object returned by element_text() to specify the appearance of the legend text}

\item{legend.position}{the position of the legend, could be 'right'(default), 'right', 'top', or 'bottom'}

\item{legend.direction}{the direction of the legend, could be 'vertical' (default) or 'horizontal'}

\item{...}{other arguments}
}
\description{
This is the plot method for the class \code{mixfitEM}. It is used to plot the fitted mixture models
by using base R plotting system or using the package ggplot2.
}
\details{
The function \code{plot.mixfitEM} is used for plotting an object of class \code{mixfitEM}, which is
an output of the function \code{\link{mixfit}}. Users can choose base R plotting system or ggplot2
(the package ggplot2 needs to be installed).
plotting system. The plot is a density plot of the fitted mixture model imposed on top of a histogram.
The parameters that control the appearance of the histogram and the density curve can be changed.
The density curve of each component can be shown or hidden.
}
\examples{
x <- rmixnormal(200, c(0.3, 0.7), c(2, 5), c(1, 0.7))
mod <- mixfit(x, ncomp = 2)
plot(mod)
plot(mod, theme = 'classic') 

}
\seealso{
\code{\link{mixfit}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmixlnorm.R
\name{rmixlnorm}
\alias{rmixlnorm}
\title{Generating Random Data From A Lognormal Mixture Model}
\usage{
rmixlnorm(n, pi, mu, sd)
}
\arguments{
\item{n}{a positive integer specifying the number of observations we want to generate from the mixture model}

\item{pi}{a numeric vector for the proportion of each component}

\item{mu}{a numeric vector for the mean of each component}

\item{sd}{a numeric vector for the standard deviation of each component}
}
\value{
The function \code{rmixlnorm} returns a numeric vector of random data from the specified lognormal mixture model.
}
\description{
The function \code{rmixlnorm} generates random data from a lognormal mixture model.
}
\details{
The number of random data from each component \eqn{n_0} (a vector) is generated from a multinomial
distribution Multinom\eqn{(n, pi)}. Then the random data from each component is generated with
the sample sized specified in \eqn{n_0} and parameters of lognormal distributions specified in
\code{mu} and \code{sd}.
}
\examples{
x <- rmixlnorm(1000, c(0.4, 0.6), c(2, 5), c(1, 0.5))
hist(x, breaks = 40)

}
\seealso{
\code{\link{rmixnormal}}, \code{\link{rmixweibull}}, \code{\link{rmixgamma}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gammapar.R
\name{to_mu_sd_gamma}
\alias{to_mu_sd_gamma}
\title{Parameter Conversion for Gamma Distribution}
\usage{
to_mu_sd_gamma(alpha, lambda)
}
\arguments{
\item{alpha}{a numeric vector representing the shape of one or more than one gamma distributions}

\item{lambda}{a numeric vector representing the rate of one or more than one gamma distributions.
\code{alpha} and \code{lambda} should have the same length.}
}
\value{
a list of two items
\item{mu}{a vector of the means of gamma distributions}
\item{sd}{a vector of the standard deviations of gamma distributions}
}
\description{
The function \code{to_mu_sd_gamma} converts the shape and rate to the mean and standard deviation
}
\details{
The purpose of this function is to convert the parameterization of gamma distribution in the form of
shape and rate to the form of mean and standard deviation.
}
\examples{
to_mu_sd_gamma(2, 1)
to_mu_sd_gamma(c(2, 4), c(1, 1))

}
\seealso{
\code{\link{to_shape_rate_gamma}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin.R
\name{bin}
\alias{bin}
\title{Binning the Raw Data}
\usage{
bin(x, brks)
}
\arguments{
\item{x}{a numeric vector of raw data}

\item{brks}{a numeric vector in increasing order, representing the bin values within each of which 
we want to calculate the frequency of the data}
}
\value{
The function \code{bin} returns a matrix with three columns, representing the value of the left bin, the value of the right bin and the number of observations in \code{x} that falls in each bin.
}
\description{
This function creates a binned data set from the raw data
}
\details{
Given a numeric vector, the function \code{bin} creates a binned data set with bin values provided
by \code{brks}. Fitting mixture models with a large data set may be slow, especially when
we want to fit non-normal mixture models. Binning the data with a relatively
small bin width speeds up the computation of EM algorithm while at the same time keeps the
precision of the estimation result.
}
\examples{
set.seed(99)
x <- rmixnormal(200, c(0.5, 0.5), c(2, 5), c(1, 1))
data <- bin(x, seq(-2, 10, 0.1))
fit1 <- mixfit(x, ncomp = 2)
fit2 <- mixfit(data, ncomp = 2)

}
\seealso{
\code{\link{reinstate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.mixfitEM.R
\name{print.mixfitEM}
\alias{print.mixfitEM}
\title{Print Method for Class \code{mixfitEM}}
\usage{
\method{print}{mixfitEM}(x, digits = getOption("digits"), ...)
}
\arguments{
\item{x}{an object of class \code{mixfitEM}}

\item{digits}{the digits to print for the values in the print output. The default
value is from the global option \code{getOption("digits")}.}

\item{...}{other arguments passed to \code{print}}
}
\description{
This function is the print method for the \code{mixfitEM} class.
}
\details{
\code{print.mixfitEM} prints the value of the parameters of a fitted mixture model, together
with some other information like the number of iterations of the EM algorithm, the loglikelihood,
the value of AIC and BIC.
}
\examples{
x <- rmixnormal(200, c(0.5, 0.5), c(2, 5), c(1, 0.7))
fit <- mixfit(x, ncomp = 2)
print(x)

}
\seealso{
\code{\link{mixfit}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{mixR-package}
\alias{mixR}
\alias{mixR-package}
\title{Finite Mixture Modeling for Raw and Binned Data}
\description{
The package \code{mixR} performs maximum likelihood estimation for finite
mixture models for families including Normal, Weibull, Gamma and Lognormal via EM algorithm.
It also conducts model selection by using information criteria or bootstrap likelihood ratio
test. The data used for mixture model fitting can be raw data or binned data. The model fitting
is accelerated by using R package Rcpp.
}
\details{
Finite mixture models can be represented by
\deqn{f(x; \Phi) = \sum_{j = 1}^g \pi_j f_j(x; \theta_j)}
where \eqn{f(x; \Phi)} is the probability density function (p.d.f.) or probability mass function
(p.m.f.) of the mixture model, \eqn{f_j(x; \theta_j)} is the p.d.f. or p.m.f. of the \eqn{j}th
component of the mixture model, \eqn{\pi_j} is the proportion of the \eqn{j}th component and
\eqn{\theta_j} is the parameter of the \eqn{j}th component, which can be a scalar or a vector,
\eqn{\Phi} is a vector of all the parameters of the mixture model. The maximum likelihood
estimate of the parameter vector \eqn{\Phi} can be obtained by using
the EM algorithm (Dempster \emph{et al}, 1977).
The binned data is present sometimes instead of the raw data, for the reason of storage
convenience or necessity. The binned data is recorded in the form of \eqn{(a_i, b_i, n_i)}
where \eqn{a_i} is the lower bound of the \eqn{i}th bin, \eqn{b_i} is
the upper bound of the \eqn{i}th bin, and \eqn{n_i} is the number of observations that fall
in the \eqn{i}th bin, for \eqn{i = 1, \dots, r}, and \eqn{r} is the total number of bins.

To obtain maximum likelihood estimate of the finite mixture model for binned data, we can
introduce two types of latent variables \eqn{x} and \eqn{z}, where\eqn{x} represents the
value of the unknown raw data, and \eqn{z} is a vector of zeros and one indicating the
component that \eqn{x} belongs to. To use the EM algorithm we first write the complete-data
log-likelihood
\deqn{Q(\Phi; \Phi^{(p)}) = \sum_{j = 1}^{g} \sum_{i = 1}^r n_i z^{(p)} [\log f(x^{(p)}; \theta_j)
 + \log \pi_j ]}
 where \eqn{z^{(p)}} is the expected value of \eqn{z} given the estimated value of \eqn{\Phi}
 and expected value \eqn{x^{(p)}} at \eqn{p}th iteration. The estimated value of \eqn{\Phi}
 can be updated iteratively via the E-step, in which we estimate \eqn{\Phi} by maximizing
 the complete-data loglikelihood, and M-step, in which we calculate the expected value of
 the latent variables \eqn{x} and \eqn{z}. The EM algorithm is terminated by using a stopping
 rule.
 The M-step of the EM algorithm may or may not have closed-form solution (e.g. the Weibull
 mixture model or Gamma mixture model). If not, an iterative approach like Newton's algorithm
 or bisection method may be used.

 For a given data set, when we have no prior information about the number of components
 \eqn{g}, its value should be estimated from the data. Because mixture models don't satisfy
 the regularity condition for the likelihood ratio test (which requires that the true
 parameter under the null hypothesis should be in the interior of the parameter space
 of the full model under the alternative hypothesis), a bootstrap approach is usually
 used in the literature (see McLachlan (1987, 2004), Feng and McCulloch (1996)). The general
 step of bootstrap likelihood ratio test is as follows.
 \enumerate{
 \item For the given data \eqn{x}, estimate \eqn{\Phi} under both the null and the alternative
 hypothesis to get \eqn{\hat\Phi_0} and \eqn{\hat\Phi_1}. Calculate the observed log-likelihood
 \eqn{\ell(x; \hat\Phi_0)} and \eqn{\ell(x; \hat\Phi_1)}. The likelihood ratio test
 statistic is defined as
 \deqn{w_0 = -2(\ell(x; \hat\Phi_0) - \ell(x; \hat\Phi_1)).}
 \item Generate random data of the same size as the original data \eqn{x} from the model
 under the null hypothesis using estimated parameter \eqn{\hat\Phi_0}, then repeat step
 1 using the simulated data. Repeat this process for \eqn{B} times to get a vector of the
 simulated likelihood ratio test statistics \eqn{w_1^{1}, \dots, w_1^{B}}.
 \item Calculate the empirical p-value
 \deqn{p = \frac{1}{B} \sum_{i=1}^B I(w_1^{(i)} > w_0)}
 where \eqn{I} is the indicator function.
 }

 This package does the following three things.
 \enumerate{
 \item Fitting finite mixture models for both raw data and binned data by using
 EM algorithm, together with Newton-Raphson algorithm and bisection method.
 \item Do parametric bootstrap likelihood ratio test for two candidate models.
 \item Do model selection by Bayesian information criterion.
 }

 To speed up computation, the EM algorithm is fulfilled in C++ by using Rcpp
 (Eddelbuettel and Francois (2011)).
}
\references{
Dempster, A. P., Laird, N. M., and Rubin, D. B. Maximum likelihood from incomplete data
via the EM algorithm. \emph{Journal of the royal statistical society. Series B
(methodological)}, pages 1-38, 1977.

Dirk Eddelbuettel and Romain Francois (2011). Rcpp: Seamless R and C++ Integration.
\emph{Journal of Statistical Software}, 40(8), 1-18. URL http://www.jstatsoft.org/v40/i08/.

Efron, B. Bootstrap methods: Another look at the jackknife. \emph{Ann. Statist.},
7(1):1-26, 01 1979.

Feng, Z. D. and McCulloch, C. E. Using bootstrap likelihood ratios in finite mixture
models. \emph{Journal of the Royal Statistical Society. Series B (Methodological)},
pages 609-617, 1996.

Lo, Y., Mendell, N. R., and Rubin, D. B. Testing the number of components in a normal
mixture. \emph{Biometrika}, 88(3):767-778, 2001.

McLachlan, G. J. On bootstrapping the likelihood ratio test statistic for the number
of components in a normal mixture. \emph{Applied statistics}, pages 318-324, 1987.

McLachlan, G. and Jones, P. Fitting mixture models to grouped and truncated data via
the EM algorithm. \emph{Biometrics}, pages 571-578, 1988.

McLachlan, G. and Peel, D. \emph{Finite mixture models}. John Wiley & Sons, 2004.
}
\author{
\strong{Maintainer}: Youjiao Yu \email{jiaoisjiao@gmail.com}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gammapar.R
\name{to_shape_rate_gamma}
\alias{to_shape_rate_gamma}
\title{Parameter Conversion for Gamma Distribution}
\usage{
to_shape_rate_gamma(mu, sd)
}
\arguments{
\item{mu}{a numeric vector representing the means of gamma distributions}

\item{sd}{a numeric vector representing the standard deviations of gamma distributions.
\code{mu} and \code{sd} should have the same length.}
}
\value{
a list of two items
\item{alpha}{a vector of the shapes of gamma distributions}
\item{lambda}{a vector of the rates of gamma distributions}
}
\description{
The function \code{to_shape_rate_gamma} converts the mean and standard deviation to the shape and rate
}
\details{
The purpose of this function is to convert the parameterization of gamma distribution in the form of
mean and standard deviation to the form of shape and rate. It can be used for specifying the initial
values for the EM algorithm when the first-hand initial values are in the form of mean and standard
deviation from K-means clustering algorithm.
}
\examples{
to_shape_rate_gamma(2, 1)
to_shape_rate_gamma(c(2, 4), c(1, 1))

}
\seealso{
\code{\link{to_mu_sd_gamma}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmixgamma.R
\name{rmixgamma}
\alias{rmixgamma}
\title{Generating Random Data From A Gamma Mixture Model}
\usage{
rmixgamma(n, pi, mu, sd)
}
\arguments{
\item{n}{a positive integer specifying the number of observations we want to generate from the mixture model}

\item{pi}{a numeric vector for the proportion of each component}

\item{mu}{a numeric vector for the mean of each component}

\item{sd}{a numeric vector for the standard deviation of each component}
}
\value{
The function \code{rmixgamma} returns a numeric vector of random data from the specified Gamma mixture model.
}
\description{
The function \code{rmixgamma} generates random data from a Gamma mixture model.
}
\details{
The number of random data from each component \eqn{n_0} (a vector) is generated from a multinomial
distribution Multinom\eqn{(n, pi)}. Then the random data from each component is generated with
the sample sized specified in \eqn{n_0} and parameters of Gamma distributions specified in
\code{mu} and \code{sd}.
}
\examples{
x <- rmixgamma(1000, c(0.4, 0.6), c(2, 5), c(1, 0.5))
hist(x, breaks = 40)

}
\seealso{
\code{\link{rmixnormal}}, \code{\link{rmixweibull}}, \code{\link{rmixlnorm}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initz.R
\name{initz}
\alias{initz}
\title{Initialization of the EM Algorithm}
\usage{
initz(x, ncomp, init.method = c("kmeans", "hclust"))
}
\arguments{
\item{x}{a numeric vector of the raw data or a three-column matrix of the binned data}

\item{ncomp}{a positive integer specifying the number of components for a mixture model}

\item{init.method}{the method used for providing initial values, which can be one of
\code{kmeans} or \code{hclust}.}
}
\value{
\code{initz} returns a list with three items
\item{pi }{a numeric vector of component proportions}
\item{mu }{a numeric vector of component means}
\item{sd }{a numeric vector of component standard deviations}
}
\description{
This function returns the mean and standard deviation of each component by using
K-means clustering or hierarchical clustering.
}
\details{
The function \code{initz} returns the mean and standard deviation of each component
of a mixture model by using K-means clustering algorithm, or hierarchical clustering
method. It is used for automatically selecting initial values for the EM algorithm,
so as to enable mixture model selection by bootstrapping likelihood ratio test or
using information criteria.
}
\examples{
x <- rmixnormal(500, c(0.5, 0.5), c(2, 5), c(1, 0.7))
data <- bin(x, seq(-2, 8, 0.25))
par1 <- initz(x, 2)
par2 <- initz(data, 2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lnormpar.R
\name{to_mulog_sdlog_lnorm}
\alias{to_mulog_sdlog_lnorm}
\title{Parameter Conversion for Lognormal Distribution}
\usage{
to_mulog_sdlog_lnorm(mu, sd)
}
\arguments{
\item{mu}{a vector of means of lognormal distributions}

\item{sd}{a vector of standard deviations of lognormal distributions}
}
\value{
a list of two items
\item{mulog}{a vector of lognormal means of lognormal distributions}
\item{sdlog}{a vector of lognormal standard deviations of lognormal distributions}
}
\description{
The function \code{to_mulog_sdlog_lnorm} converts the mean and standard deviation to
the logarithm mean and logarithm standard deviation
}
\details{
The purpose of this function is to convert the parameterization of lognormal distribution in the
form of mean and standard deviation to the form of logarithm mean and logarithm standard deviation.
It can be used for specifying the initial values for the EM algorithm when the first-hand initial values
are in the form of mean and standard deviation from K-means clustering algorithm.
}
\examples{
to_mulog_sdlog_lnorm(2, 1)
to_mulog_sdlog_lnorm(c(2, 4), c(1, 1))

}
\seealso{
\code{\link{to_mu_sd_lnorm}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reinstate.R
\name{reinstate}
\alias{reinstate}
\title{Reinstate the Binned Data to the Raw Data}
\usage{
reinstate(data)
}
\arguments{
\item{data}{a three-column matrix representing the raw data}
}
\value{
The function returns a numeric vector.
}
\description{
This function creates a numeric vector approximating the raw data from binned data
}
\details{
The function \code{reinstate} creates a numeric vector by generating \eqn{n_i}
random data from the Uniform distribution \eqn{U(a_i, b_i)} for \eqn{i = 1, \dots, r}
and then combine all random data together. \eqn{a_i, b_i, n_i}
are the first, second and the third column of the matrix \code{data}
and \eqn{r} is the number of bins.
It is used for enabling parameter initialization for EM algorithm when we fit mixture
models for binned data.
}
\examples{
x <- rnorm(100)
data <- bin(x, seq(-3, 3, 0.25))
y <- reinstate(data)

}
\seealso{
\code{\link{bin}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.selectEM.R
\name{plot.selectEM}
\alias{plot.selectEM}
\title{Plot Method for Class \code{selectEM}}
\usage{
\method{plot}{selectEM}(x, leg.loc = "topright", ...)
}
\arguments{
\item{x}{an object of class \code{selectEM}, which is an output of the function \code{\link{select}}.}

\item{leg.loc}{the location of the legend, which is the same as the first argument of the function}

\item{...}{other arguments passed to \code{plot}
\code{\link{legend}}. The default value is "topright". The user can change its location
(to "topleft", "bottom right" etc.) if the visual plot conflicts with the legend.}
}
\description{
This function plots the result of mixture model selection by BIC.
}
\details{
The function \code{plot.selectEM} is the plot method for the class \code{selectEM}. It plots
the number of components against the corresponding value of BIC. It is used to visually display
the mixture model selection result by BIC.
}
\examples{
x <- rmixnormal(200, c(0.3, 0.7), c(2, 5), c(1, 1))
res <- select(x, ncomp = 1:3)
plot(res)

}
\seealso{
\code{\link{select}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibullpar.R
\name{to_mu_sd_weibull}
\alias{to_mu_sd_weibull}
\title{Parameter Conversion for Weibull Distribution}
\usage{
to_mu_sd_weibull(k, lambda)
}
\arguments{
\item{k}{a numeric vector representing the shape of a series of Weibull distributions}

\item{lambda}{a numeric vector representing the scale of a series of Weibull distributions.
\code{k} and \code{lambda} should have the same length.}
}
\value{
a list of two items
\item{mu}{a vector of the means of Weibull distributions}
\item{sd}{a vector of the standard deviations of Weibull distributions}
}
\description{
The function \code{to_mu_sd_weibull} converts the parameters of shape and scale of weibull distributions to
the parameters of the mean and standard deviation.
}
\details{
The purpose of this function is to convert the parameterization of Weibull distribution in the form of
shape and scale to the form of mean and standard deviation.
}
\examples{
to_mu_sd_weibull(2, 1)
to_mu_sd_weibull(c(2, 4), c(1, 1))

}
\seealso{
\code{\link{to_k_lambda_weibull}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmixweibull.R
\name{rmixweibull}
\alias{rmixweibull}
\title{Generating Random Data From A Weibull Mixture Model}
\usage{
rmixweibull(n, pi, mu, sd)
}
\arguments{
\item{n}{a positive integer specifying the number of observations we want to generate from the mixture model}

\item{pi}{a numeric vector for the proportion of each component}

\item{mu}{a numeric vector for the mean of each component}

\item{sd}{a numeric vector for the standard deviation of each component}
}
\value{
The function \code{rmixweibull} returns a numeric vector of random data from the specified Weibull mixture model.
}
\description{
The function \code{rmixweibull} generates random data from a normal Weibull model.
}
\details{
The number of random data from each component \eqn{n_0} (a vector) is generated from a multinomial
distribution Multinom\eqn{(n, pi)}. Then the random data from each component is generated with
the sample sized specified in \eqn{n_0} and parameters of Weibull distributions specified in
\code{mu} and \code{sd}.
}
\examples{
x <- rmixweibull(1000, c(0.4, 0.6), c(2, 5), c(1, 0.5))
hist(x, breaks = 40)

}
\seealso{
\code{\link{rmixnormal}}, \code{\link{rmixgamma}}, \code{\link{rmixlnorm}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_Stamp2.R
\docType{data}
\name{Stamp2}
\alias{Stamp2}
\title{1872 Hidalgo Stamp Data (Binned)}
\format{
A matrix with 62 rows and 3 columns:
\describe{
 \item{lower}{the lower bin values}
 \item{upper}{the upper bin values}
 \item{freq}{the number of observations in each bin}
}
}
\usage{
Stamp2
}
\description{
A dataset containing the 1872 Hidalgo stamp data in the form of binned data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.selectEM.R
\name{print.selectEM}
\alias{print.selectEM}
\title{Print Method for Class \code{selectEM}}
\usage{
\method{print}{selectEM}(x, ...)
}
\arguments{
\item{x}{an object of class \code{selectEM}}

\item{...}{other arguments passed to \code{print}}
}
\description{
The function prints the result of mixture model selection.
}
\details{
The function \code{print.selectEM} is the print method for the class{selectEM}, which is the output
of the function \code{select}. It prints a data frame which contains the following information of each
candidate mixture models: the number of components, whether the variance is the same for each component
in a mixture model (only for normal), the value of BIC, and an indicator of the best model.
}
\seealso{
\code{\link{select}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{select}
\alias{select}
\title{Finite Mixture Model Selection by Information Criterion}
\usage{
select(
  x,
  ncomp,
  family = c("normal", "weibull", "gamma", "lnorm"),
  mstep.method = c("bisection", "newton"),
  init.method = c("kmeans", "hclust"),
  tol = 1e-06,
  max_iter = 500
)
}
\arguments{
\item{x}{a numeric vector for raw data or a three-column matrix for the binned data}

\item{ncomp}{a vector of positive integers specifying the number of components of the candidate mixture models}

\item{family}{a character string specifying the family of the mixture model. It can only be
one element from \code{normal}, \code{weibull}, \code{gamma} or \code{lnorm}.}

\item{mstep.method}{a character string specifying the method used in M-step of the EM algorithm
when fitting weibull or gamma mixture models. It can be either \code{bisection} or \code{newton}.
The default is \code{bisection}.}

\item{init.method}{a character string specifying the method used for providing initial values
for the parameters for EM algorithm. It can be one of \code{kmeans} or \code{hclust}. The default is
\code{kmeans}}

\item{tol}{the tolerance for the stopping rule of EM algorithm. It is the value to stop EM algorithm when the two
consecutive iterations produces loglikelihood with difference less than \code{tol}. The default value is 1e-6.}

\item{max_iter}{the maximum number of iterations for the EM algorithm (default 500).}
}
\value{
The function returns an object of class \code{selectEM} which contains the following items.
        \item{ncomp}{the specified number of components of the candidate mixture models}
        \item{equal.var}{a logical vector indicating whether the variances of each component in each mixture model
        are constrained to be the same (only for \code{normal} family)}
        \item{bic}{the value of BIC for each mixture model}
        \item{best}{an indicator of the best model}
        \item{family}{the family of the mixture model}
}
\description{
This function selects the best model from a candidate of mixture models based on the
information criterion BIC.
}
\details{
By specifying different number of components, the function \code{select} fits a series
of mixture models for a given family, and a mixture model with minimum value of BIC
is regarded as the best.
}
\examples{
## selecting the optimal normal mixture model by BIC
set.seed(105)
x <- rmixnormal(1000, c(0.3, 0.4, 0.3), c(-4, 0, 4), c(1, 1, 1))
hist(x, breaks = 40)
ret <- select(x, ncomp = 2:5)
## [1] "The final model: normal mixture (equal variance) with 3 components"

## (not run) selecting the optimal Weibull mixture model by BIC
## set.seed(106)
## x <- rmixweibull(1000, c(0.3, 0.4, 0.3), c(2, 5, 8), c(0.7, 0.6, 1))
## ret <- select(x, ncomp = 2:5, family = "weibull")
## [1] "The final model: weibull mixture with 3 components"

## (not run) selecting the optimal Gamma mixture model by BIC
## set.seed(107)
## x <- rmixgamma(1000, c(0.3, 0.7), c(2, 5), c(0.7, 1))
## ret <- select(x, ncomp = 2:5, family = "gamma")
## [1] "The final model: gamma mixture with 2 components"


## (not run) selecting the optimal lognormal mixture model by BIC
## set.seed(108)
## x <- rmixlnorm(1000, c(0.2, 0.3, 0.2, 0.3), c(4, 7, 9, 12), c(1, 0.5, 0.7, 1))
## ret <- select(x, ncomp = 2:6, family = "lnorm")
## [1] "The final model: lnorm mixture with 4 components"

}
\seealso{
\code{\link{plot.selectEM}}, \code{\link{bs.test}}, \code{\link{mixfit}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs.test.R
\name{bs.test}
\alias{bs.test}
\title{Bootstrap Likelihood Ratio Test for Finite Mixture Models}
\usage{
bs.test(
  x,
  ncomp = c(1, 2),
  family = c("normal", "weibull", "gamma", "lnorm"),
  B = 100,
  ev = FALSE,
  mstep.method = c("bisection", "newton"),
  init.method = c("kmeans", "hclust"),
  tol = 1e-06,
  max_iter = 500
)
}
\arguments{
\item{x}{a numeric vector for the raw data or a three-column matrix for the binned data.}

\item{ncomp}{a vector of two positive integers specifying the number of components of the
mixture model under the null and alternative hypothesis.
The first integer should be smaller than the second one. The default value is
\code{c(1, 2)}.}

\item{family}{a character string specifying the family of the mixture model, which can be one
of \code{normal}, \code{weibull}, \code{gamma}, or \code{lnorm} (default \code{normal}).}

\item{B}{the number of bootstrap iterations (default 100).}

\item{ev}{a logical value indicating whether the variance of each component should be the same
or not (default \code{FALSE} for \code{Normal} family and ignored for other family members).}

\item{mstep.method}{the method used in M-step of EM algorithm for \code{weibull} or
\code{gamma} family. It is ignored for \code{normal} or \code{lnorm} family,
which has closed-form solution in the M-step. The default value is \code{bisection}.}

\item{init.method}{a character string specifying the method used for providing the initial values
for the parameters for the EM algorithm. It can be one of \code{kmeans} or \code{hclust}. The default is
\code{kmeans}}

\item{tol}{the tolerance for the stopping rule of EM algorithm. It is the value to stop
EM algorithm when the two consecutive iterations produces log-likelihood with difference
less than \code{tol}. The default value is 1e-6.}

\item{max_iter}{the maximum number of iterations for the EM algorithm (default 500).}
}
\value{
The function \code{bs.test} returns an object of class \code{bootEM} which
contains the following three items.
\item{pvalue}{The p-value of the bootstrap likelihood ratio test}
\item{w0}{the observed likelihood ratio test statistic}
\item{w1}{a vector of simulated likelihood ratio test statistics}
}
\description{
This function performs the likelihood ratio test by parametric bootstrapping for two mixture
models with different number of components.
}
\details{
For the given data \code{x} and the specified family, the function \code{bs.test} conducts
a bootstrap likelihood ratio test for two mixture models with the number of components
under the null and the alternative hypothesis specified in \code{ncomp}.
}
\examples{
## testing normal mixture models with 2 and 3 components
set.seed(100)
x <- rmixnormal(200, c(0.5, 0.5), c(2, 5), c(1, 0.7))
ret <- bs.test(x, ncomp = c(2, 3), B = 30)
ret

## (not run) testing Weibull mixture models with 2 and 3 components
## set.seed(101)
## x <- rmixweibull(200, c(0.3, 0.4, 0.3), c(2, 5, 8), c(1, 0.6, 0.8))
## ret <- bs.test(x, ncomp = c(2, 3), family = "weibull", B = 30)
## ret

## (not run) testing Gamma mixture models with 1 and 2 components
## set.seed(102)
## x <- rgamma(200, 2, 1)
## ret <- bs.test(x, ncomp = c(1, 2), family = "gamma", B = 30)
## ret

}
\seealso{
\code{\link{plot.bootEM}}, \code{\link{mixfit}}, \code{\link{select}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/density.mixfitEM.R
\name{density.mixfitEM}
\alias{density.mixfitEM}
\title{The Density of Finite Mixture Models}
\usage{
\method{density}{mixfitEM}(x, at, smoothness = 512, cut = 3.8, ...)
}
\arguments{
\item{x}{an object of class \code{mixfitEM}}

\item{at}{a scalar or a numeric vector of locations where densities are calculated}

\item{smoothness}{a positive integer controlling the smoothness of the density curve (default 512).
The higher this value is, the more locations of the mixture model the density is calculated.}

\item{cut}{the number of standard deviations away the density is to be computed (default 3.8)}

\item{...}{other arguments}
}
\value{
This function returns a list of class \code{densityEM}, which contains the following
items.
\item{x}{a scalar or a numeric vector of locations where densities are calculated.}
\item{y}{a vector of the densities of the mixture model at the corresponding locations in \code{x}}
\item{comp}{a matrix with columns representing the densities of each component in the mixture model at the corresponding locations in \code{x}}
}
\description{
This function calculates the probability density of a finite mixture model.
}
\details{
The function \code{density.mixfitEM} is the method of the generic function
\code{density} for the class \code{mixfitEM}.
}
\examples{
set.seed(102)
x <- rmixnormal(200, c(0.5, 0.5), c(2, 5), c(1, 0.7))
fit1 <- mixfit(x, ncomp = 2)
d1 = density(fit1)
d2 = density(fit1, at = 0)

}
\seealso{
\code{\link{mixfit}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lnormpar.R
\name{to_mu_sd_lnorm}
\alias{to_mu_sd_lnorm}
\title{Parameter Conversion for Lognormal Distribution}
\usage{
to_mu_sd_lnorm(mulog, sdlog)
}
\arguments{
\item{mulog}{a vector of logarithm means of lognormal distributions}

\item{sdlog}{a vector of logarithm standard deviations of lognormal distributions}
}
\value{
a list of two items
\item{mu}{a vector of the means of lognormal distributions}
\item{sd}{a vector of the standard deviations of lognormal distributions}
}
\description{
The function \code{to_mu_sd_lnorm} converts the logarithm mean and logarithm standard deviation to the
mean and standard deviation
}
\details{
The purpose of this function is to convert the parameterization of lognormal distribution in the
form of logarithm mean and logarithm standard deviation to the form of mean and standard deviation.
}
\examples{
to_mu_sd_lnorm(2, 1)
to_mu_sd_lnorm(c(2, 4), c(1, 1))

}
\seealso{
\code{\link{to_mulog_sdlog_lnorm}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{mixfit}
\alias{mixfit}
\title{Finite Mixture Modeling for Raw Data and Binned Data}
\usage{
mixfit(
  x,
  ncomp = NULL,
  family = c("normal", "weibull", "gamma", "lnorm"),
  pi = NULL,
  mu = NULL,
  sd = NULL,
  ev = FALSE,
  mstep.method = c("bisection", "newton"),
  init.method = c("kmeans", "hclust"),
  tol = 1e-06,
  max_iter = 500
)
}
\arguments{
\item{x}{a numeric vector for the raw data or a three-column matrix for the binned data}

\item{ncomp}{a positive integer specifying the number of components of the mixture model}

\item{family}{a character string specifying the family of the mixture model. It can only be
one element from \code{normal}, \code{weibull}, \code{gamma} or \code{lnorm}.}

\item{pi}{a vector of the initial value for the proportion}

\item{mu}{a vector of the initial value for the mean}

\item{sd}{a vector of the initial value for the standard deviation}

\item{ev}{a logical value controlling whether each component has the same variance when
fitting normal mixture models. It is ignored when fitting other mixture models. The default is \code{FALSE}.}

\item{mstep.method}{a character string specifying the method used in M-step of the EM algorithm
when fitting weibull or gamma mixture models. It can be either \code{bisection} or \code{newton}.
The default is \code{bisection}.}

\item{init.method}{a character string specifying the method used for providing initial values
for the parameters for EM algorithm. It can be one of \code{kmeans} or \code{hclust}. The default is
\code{kmeans}}

\item{tol}{the tolerance for the stopping rule of EM algorithm. It is the value to stop EM algorithm when the two
consecutive iterations produces loglikelihood with difference less than \code{tol}. The default value is 1e-6.}

\item{max_iter}{the maximum number of iterations for the EM algorithm (default 500).}
}
\value{
the function \code{mixfit} return an object of class \code{mixfitEM}, which contains a list of
different number of items when fitting different mixture models. The common items include
\item{pi}{a numeric vector representing the estimated proportion of each component}
\item{mu}{a numeric vector representing the estimated mean of each component}
\item{sd}{a numeric vector representing the estimated standard deviation of each component}
\item{iter}{a positive integer recording the number of EM iteration performed}
\item{loglik}{the loglikelihood of the estimated mixture model for the data \code{x}}
\item{aic}{the value of AIC of the estimated model for the data \code{x}}
\item{bic}{the value of BIC of the estimated model for the data \code{x}}
\item{data}{the data \code{x}}
\item{comp.prob}{the probability that \code{x} belongs to each component}
\item{family}{the family the mixture model belongs to}
For the Weibull mixture model, the following extra items are returned.
\item{k}{a numeric vector representing the estimated shape parameter of each component}
\item{lambda}{a numeric vector representing the estimated scale parameter of each component}
For the Gamma mixture model, the following extra items are returned.
\item{alpha}{a numeric vector representing the estimated shape parameter of each component}
\item{lambda}{a numeric vector representing the estimated rate parameter of each component}
For the lognormal mixture model, the following extra items are returned.
\item{mulog}{a numeric vector representing the estimated logarithm mean of each component}
\item{sdlog}{a numeric vector representing the estimated logarithm standard deviation of
each component}
}
\description{
This function is used to perform the maximum likelihood estimation for
a variety of finite mixture models for both raw and binned data by using
the EM algorithm, together with Newton-Raphson algorithm or bisection method when necessary.
}
\details{
The function \code{mixfit} is the core function in this package. It is used to perform
the maximum likelihood estimation for finite mixture models from the families of normal,
weibull, gamma or lognormal by using the EM algorithm. When the family is \code{weibull}
or \code{gamma}, the M-step of the EM algorithm has no closed-form solution and we can
use Newton algorithm by specifying \code{method = "newton"} or use bisection method by
specifying \code{method = "bisection"}.

The initial values of the EM algorithm can be provided by specifying the proportion of each
component \code{pi}, the mean of each component \code{mu} and the standard deviation of
each component \code{sd}. If one or more of these initial values are not provided, then
their values are estimated by using K-means clustering method or hierarchical clustering
method. If all of \code{pi}, \code{mu}, and \code{sd}
are not provided, then \code{ncomp} should be provided so initial values are automatically
generated. For the normal mixture models, we can
control whether each component has the same variance or not.
}
\examples{
## fitting the normal mixture models
set.seed(103)
x <- rmixnormal(200, c(0.3, 0.7), c(2, 5), c(1, 1))
data <- bin(x, seq(-1, 8, 0.25))
fit1 <- mixfit(x, ncomp = 2)  # raw data
fit2 <- mixfit(data, ncomp = 2)  # binned data
fit3 <- mixfit(x, pi = c(0.5, 0.5), mu = c(1, 4), sd = c(1, 1))  # providing the initial values
fit4 <- mixfit(x, ncomp = 2, ev = TRUE)  # setting the same variance

## (not run) fitting the weibull mixture models
## x <- rmixweibull(200, c(0.3, 0.7), c(2, 5), c(1, 1))
## data <- bin(x, seq(0, 8, 0.25))
## fit5 <- mixfit(x, ncomp = 2, family = "weibull")  # raw data
## fit6 <- mixfit(data, ncomp = 2, family = "weibull")  # binned data

## (not run) fitting the Gamma mixture models
## x <- rmixgamma(200, c(0.3, 0.7), c(2, 5), c(1, 1))
## data <- bin(x, seq(0, 8, 0.25))
## fit7 <- mixfit(x, ncomp = 2, family = "gamma")  # raw data
## fit8 <- mixfit(data, ncomp = 2, family = "gamma")  # binned data

## (not run) fitting the lognormal mixture models
## x <- rmixlnorm(200, c(0.3, 0.7), c(2, 5), c(1, 1))
## data <- bin(x, seq(0, 8, 0.25))
## fit9 <- mixfit(x, ncomp = 2, family = "lnorm")  # raw data
## fit10 <- mixfit(data, ncomp = 2, family = "lnorm")  # binned data

}
\seealso{
\code{\link{plot.mixfitEM}}, \code{\link{density.mixfitEM}},
\code{\link{select}}, \code{\link{bs.test}}
}
