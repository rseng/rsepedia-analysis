---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: 'BGGM: Bayesian Gaussian Graphical Models in R'
tags:
- Gaussian graphical models
- Bayesian
- Bayes factor
- partial correlation
- R
authors:
  - name: Donald R. Williams
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Joris Mulder
    affiliation: 2
affiliations:
 - name: Department of Psychology, University of California, Davis
   index: 1
 - name: Department of Methodology and Statistics, Tilburg University
   index: 2
citation_author: Williams and Mulder
date: 05 May 2020
year: 2020
bibliography: inst/REFERENCES.bib
---

# BGGM: Bayesian Gaussian Graphical Models
The `R` package **BGGM** provides tools for making Bayesian inference in 
Gaussian graphical models (GGM). The methods are organized around two general 
approaches for Bayesian inference: (1) estimation and (2) hypothesis 
testing. The key distinction is that the former focuses on either 
the posterior or posterior predictive distribution [@Gelman1996a; see section 
5 in @rubin1984bayesianly] , whereas the latter focuses on model comparison with the Bayes factor [@Jeffreys1961; @Kass1995].


## What is a Gaussian Graphical Model ?
A Gaussian graphical model captures conditional (in)dependencies among a set 
of variables. These are pairwise relations (partial correlations) controlling for 
the effects of all other variables in the model.

### Applications
The Gaussian graphical model is used across the sciences, including 
(but not limited to) economics [@millington2020partial], climate science 
[@zerenner2014gaussian], genetics [@chu2009graphical], and psychology [@rodriguez2020formalizing]. 

# Overview
The methods in **BGGM** build upon existing algorithms that are well-known in the literature.
The central contribution of **BGGM** is to extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution [@Mulder2018]
  
    + [Estimation](https://github.com/donaldRwilliams/BGGM#bayesian-estimation) [@Williams2019]

2. Bayesian hypothesis testing with the matrix-F prior distribution [@Williams2019_bf]

    + [Exploratory hypothesis testing](https://github.com/donaldRwilliams/BGGM#Exploratory)
  
    + [Confirmatory hypothesis testing](https://github.com/donaldRwilliams/BGGM#Confirmatory)
    
3. Comparing Gaussian graphical models [@Williams2019; @williams2020comparing]
    
    + [Partial correlation differences](https://github.com/donaldRwilliams/BGGM#partial-correlation-differences) 
    
    + [Posterior predictive check](https://github.com/donaldRwilliams/BGGM#posterior-predictive-check)
    
    + [Exploratory hypothesis testing](https://github.com/donaldRwilliams/BGGM#exploratory-groups) 
    
    + [Confirmatory hypothesis testing](https://github.com/donaldRwilliams/BGGM#confirmatory-groups)

4. Extending inference beyond the conditional (in)dependence structure [@Williams2019]

    +  [Predictability](https://github.com/donaldRwilliams/BGGM#Predictability)[e.g., @haslbeck2018well]
    
    +  [Posterior uncertaintyintervals](https://github.com/donaldRwilliams/BGGM#partial-correlation-differences) for the 
       partial correlations
       
    +  [Custom Network Statistics](https://github.com/donaldRwilliams/BGGM#custom-network-statistics)
    
    
## Supported Data Types

* **Continuous**: The continuous method was described in  @Williams2019. Note that 
                  this is based on the customary [Wishartdistribution](https://en.wikipedia.org/wiki/Wishart_distribution).

* **Binary**: The binary method builds directly upon @talhouk2012efficient
  that, in turn, built upon the approaches of @lawrence2008bayesian and
  @webb2008bayesian (to name a few).
  
* **Ordinal**: The ordinal methods require sampling thresholds. There are two approach 
   included in **BGGM**. The customary approach described in @albert1993bayesian 
   (the default) and the 'Cowles' algorithm described in @cowles1996accelerating.
   
* **Mixed**: The mixed data (a combination of discrete and continuous) method was introduced
 in @hoff2007extending. This is a semi-parametric copula model
 (i.e., a copula GGM) based on the ranked likelihood. Note that this can be used for 
 *only* ordinal data (not restricted to "mixed" data).

The computationally intensive tasks are written in `c++` via the `R` package **Rcpp** [@eddelbuettel2011rcpp] and the `c++` library **Armadillo** [@sanderson2016armadillo]. The Bayes factors are computed with the `R` package **BFpack** [@mulder2019bfpack]. Furthermore, there are [plotting](https://github.com/donaldRwilliams/BGGM#example-network-plot) functions
for each method, control variables can be included in the model (e.g., `~ gender`), 
and there is support for missing values (see `bggm_missing`).

## Comparison to Other Software
**BGGM** is the only `R` package to implement all of these algorithms and methods. The `mixed` data approach 
is also implemented in the package **sbgcop** [base `R`, @hoff2007extending]. The `R` package **BDgraph** implements a Gaussian copula graphical model in `c++` [@mohammadi2015bdgraph], but not the binary or ordinal approaches. Furthermore, **BGGM** is the only package for confirmatory testing and comparing graphical models with the methods described in @williams2020comparing.

# Acknowledgements
DRW was supported by a National Science Foundation Graduate Research Fellowship
under Grant No. 1650042 and JM was supported by a ERC Starting Grant (758791).

# References

# Bayesian Gaussian Graphical Models <img src="man/figures/logo.png" align="right" alt="" width="150" />

<!-- badges: start -->

[![CRAN
Version](http://www.r-pkg.org/badges/version/BGGM)](https://cran.r-project.org/package=BGGM)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BGGM)
[![Build
Status](https://travis-ci.org/donaldRwilliams/BGGM.svg?branch=master)](https://travis-ci.org/donaldRwilliams/BGGM)
<!-- badges: end -->

The `R` package **BGGM** provides tools for making Bayesian inference in
Gaussian graphical models (GGM, Williams and Mulder 2020). The methods
are organized around two general approaches for Bayesian inference: (1)
estimation and (2) hypothesis testing. The key distinction is that the
former focuses on either the posterior or posterior predictive
distribution (Gelman, Meng, and Stern 1996; see section 5 in Rubin
1984), whereas the latter focuses on model comparison with the Bayes
factor (Jeffreys 1961; Kass and Raftery 1995).

## <i class="fas fa-cog"></i> Installation

To install the latest release version (`2.0.0`) from CRAN use

``` r
install.packages("BGGM")    
```

The current developmental version can be installed with

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("donaldRwilliams/BGGM")
```

### <i class="fas fa-skull-crossbones"></i> Dealing with Errors

There are automatic checks for
[**BGGM**](https://travis-ci.org/github/donaldRwilliams/BGGM/branches).
However, that only checks for Linux and **BGGM** is built on Windows.
The most common installation errors occur on OSX. An evolving guide to
address these issues is provided in the [Troubleshoot
Section](https://donaldrwilliams.github.io/BGGM/articles/installation.html).

## <i class="fas fa-clipboard-list"></i> Overview

The methods in **BGGM** build upon existing algorithms that are
well-known in the literature. The central contribution of **BGGM** is to
extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution
    (Mulder and Pericchi 2018)
    
      - Estimation (Williams 2018)

2.  Bayesian hypothesis testing with the matrix-F prior distribution
    (Williams and Mulder 2019)
    
      - [Exploratory hypothesis
        testing](https://github.com/donaldRwilliams/BGGM#Exploratory)
    
      - [Confirmatory hypothesis
        testing](https://github.com/donaldRwilliams/BGGM#confirmatory)

3.  Comparing Gaussian graphical models (Williams 2018; Williams et al.
    2020)
    
      - [Partial correlation
        differences](https://github.com/donaldRwilliams/BGGM#partial-correlation-differences)
    
      - [Posterior predictive
        check](https://github.com/donaldRwilliams/BGGM#posterior-predictive-check)
    
      - [Exploratory hypothesis
        testing](https://github.com/donaldRwilliams/BGGM#exploratory-groups)
    
      - [Confirmatory hypothesis
        testing](https://github.com/donaldRwilliams/BGGM#confirmatory-groups)

4.  Extending inference beyond the conditional (in)dependence structure
    (Williams 2018)
    
      - [Predictability](https://github.com/donaldRwilliams/BGGM#predictability)
    
      - [Posterior uncertainty
        intervals](https://github.com/donaldRwilliams/BGGM#posterior-uncertainty)
        for the partial correlations
    
      - [Custom Network
        Statistics](https://github.com/donaldRwilliams/BGGM#custom-network-statistics)

The computationally intensive tasks are written in `c++` via the `R`
package **Rcpp** (Eddelbuettel et al. 2011) and the `c++` library
**Armadillo** (Sanderson and Curtin 2016). The Bayes factors are
computed with the `R` package **BFpack** (Mulder et al. 2019).
Furthermore, there are plotting functions for each method, control
variables can be included in the model (e.g., `~ gender`), and there is
support for missing values (see `bggm_missing`).

## Supported Data Types

  - **Continuous**: The continuous method was described in Williams
    (2018). Note that this is based on the customary [Wishart
    distribution](https://en.wikipedia.org/wiki/Wishart_distribution).

  - **Binary**: The binary method builds directly upon Talhouk, Doucet,
    and Murphy (2012) that, in turn, built upon the approaches of
    Lawrence et al. (2008) and Webb and Forster (2008) (to name a few).

  - **Ordinal**: The ordinal methods require sampling thresholds. There
    are two approach included in **BGGM**. The customary approach
    described in Albert and Chib (1993) (the default) and the ‘Cowles’
    algorithm described in Cowles (1996).

  - **Mixed**: The mixed data (a combination of discrete and continuous)
    method was introduced in Hoff (2007). This is a semi-parametric
    copula model (i.e., a copula GGM) based on the ranked likelihood.
    Note that this can be used for *only* ordinal data (not restricted
    to “mixed” data).

## <i class="fas fa-folder-open"></i> Illustrative Examples

There are several examples in the
[Vignettes](https://donaldrwilliams.github.io/BGGM/articles/) section.

## <i class="fas fa-play-circle"></i> Basic Usage

It is common to have some combination of continuous and discrete (e.g.,
ordinal, binary, etc.) variables. **BGGM** (as of version `2.0.0`) can
readily be used for these kinds of data. In this example, a model is
fitted for the `gss` data in **BGGM**.

### Visualize

The data are first visualized with the **psych** package, which readily
shows the data are “mixed”.

``` r
# dev version
library(BGGM)
library(psych)

# data
Y <- gss

# histogram for each node
psych::multi.hist(Y, density = FALSE)
```

![](man/figures/index_hist.png)

### Fit Model

A Gaussian copula graphical model is estimated as follows

``` r
fit <- estimate(Y, type = "mixed")
```

`type` can be `continuous`, `binary`, `ordinal`, or `mixed`. Note that
`type` is a misnomer, as the data can consist of *only* ordinal
variables (for example).

### Summarize Relations

The estimated relations are summarized with

``` r
summary(fit)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: mixed 
#> Analytic: FALSE 
#> Formula:  
#> Posterior Samples: 5000 
#> Observations (n): 464  
#> Nodes (p): 7 
#> Relations: 21 
#> --- 
#> Call: 
#> estimate(Y = Y, type = "mixed")
#> --- 
#> Estimates:
#>  Relation Post.mean Post.sd Cred.lb Cred.ub
#>  INC--DEG     0.463   0.042   0.377   0.544
#>  INC--CHI     0.148   0.053   0.047   0.251
#>  DEG--CHI    -0.133   0.058  -0.244  -0.018
#>  INC--PIN     0.087   0.054  -0.019   0.196
#>  DEG--PIN    -0.050   0.058  -0.165   0.062
#>  CHI--PIN    -0.045   0.057  -0.155   0.067
#>  INC--PDE     0.061   0.057  -0.050   0.175
#>  DEG--PDE     0.326   0.056   0.221   0.438
#>  CHI--PDE    -0.043   0.062  -0.162   0.078
#>  PIN--PDE     0.345   0.059   0.239   0.468
#>  INC--PCH     0.052   0.052  -0.052   0.150
#>  DEG--PCH    -0.121   0.056  -0.228  -0.012
#>  CHI--PCH     0.113   0.056   0.007   0.224
#>  PIN--PCH    -0.080   0.059  -0.185   0.052
#>  PDE--PCH    -0.200   0.058  -0.305  -0.082
#>  INC--AGE     0.211   0.050   0.107   0.306
#>  DEG--AGE     0.046   0.055  -0.061   0.156
#>  CHI--AGE     0.522   0.039   0.442   0.594
#>  PIN--AGE    -0.020   0.054  -0.122   0.085
#>  PDE--AGE    -0.141   0.057  -0.251  -0.030
#>  PCH--AGE    -0.033   0.051  -0.132   0.063
#> --- 
```

The summary can also be plotted

``` r
plot(summary(fit))
```

![](man/figures/index_summ.png)

### Graph Selection

The graph is selected and plotted with

``` r
E <- select(fit)

plot(E, node_size = 12,
     edge_magnify = 5)
```

![](man/figures/netplot_index.png)

The Bayes factor testing approach is readily implemented by changing
`estimate` to `explore`.

## <i class="fas fa-pen-square"></i> References

<div id="refs" class="references">

<div id="ref-albert1993bayesian">

Albert, James H, and Siddhartha Chib. 1993. “Bayesian Analysis of Binary
and Polychotomous Response Data.” *Journal of the American Statistical
Association* 88 (422): 669–79.

</div>

<div id="ref-cowles1996accelerating">

Cowles, Mary Kathryn. 1996. “Accelerating Monte Carlo Markov Chain
Convergence for Cumulative-Link Generalized Linear Models.” *Statistics
and Computing* 6 (2): 101–11.

</div>

<div id="ref-eddelbuettel2011rcpp">

Eddelbuettel, Dirk, Romain François, J Allaire, Kevin Ushey, Qiang Kou,
N Russel, John Chambers, and D Bates. 2011. “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software* 40 (8): 1–18.

</div>

<div id="ref-Gelman1996a">

Gelman, Andrew, Xiao-Li Meng, and Hal Stern. 1996. “Posterior predictive
assessment of model fitness via realized discrepancies. Vol.6, No.4.”
*Statistica Sinica* 6 (4): 733–807. <https://doi.org/10.1.1.142.9951>.

</div>

<div id="ref-hoff2007extending">

Hoff, Peter D. 2007. “Extending the Rank Likelihood for Semiparametric
Copula Estimation.” *The Annals of Applied Statistics* 1 (1): 265–83.

</div>

<div id="ref-Jeffreys1961">

Jeffreys, Harold. 1961. *The theory of probability*. Oxford: Oxford
University Press.

</div>

<div id="ref-Kass1995">

Kass, Robert E, and Adrian E Raftery. 1995. “Bayes Factors.” *Journal of
the American Statistical Association* 90 (430): 773–95.

</div>

<div id="ref-lawrence2008bayesian">

Lawrence, Earl, Derek Bingham, Chuanhai Liu, and Vijayan N Nair. 2008.
“Bayesian Inference for Multivariate Ordinal Data Using Parameter
Expansion.” *Technometrics* 50 (2): 182–91.

</div>

<div id="ref-mulder2019bfpack">

Mulder, Joris, Xin Gu, Anton Olsson-Collentine, Andrew Tomarken, Florian
Böing-Messing, Herbert Hoijtink, Marlyne Meijerink, et al. 2019.
“BFpack: Flexible Bayes Factor Testing of Scientific Theories in R.”
*arXiv Preprint arXiv:1911.07728*.

</div>

<div id="ref-Mulder2018">

Mulder, Joris, and Luis Pericchi. 2018. “The Matrix-F Prior for
Estimating and Testing Covariance Matrices.” *Bayesian Analysis*, no. 4:
1–22. <https://doi.org/10.1214/17-BA1092>.

</div>

<div id="ref-rubin1984bayesianly">

Rubin, Donald B. 1984. “Bayesianly Justifiable and Relevant Frequency
Calculations for the Applied Statistician.” *The Annals of Statistics*,
1151–72.

</div>

<div id="ref-sanderson2016armadillo">

Sanderson, Conrad, and Ryan Curtin. 2016. “Armadillo: A Template-Based
C++ Library for Linear Algebra.” *Journal of Open Source Software* 1
(2): 26.

</div>

<div id="ref-talhouk2012efficient">

Talhouk, Aline, Arnaud Doucet, and Kevin Murphy. 2012. “Efficient
Bayesian Inference for Multivariate Probit Models with Sparse Inverse
Correlation Matrices.” *Journal of Computational and Graphical
Statistics* 21 (3): 739–57.

</div>

<div id="ref-webb2008bayesian">

Webb, Emily L, and Jonathan J Forster. 2008. “Bayesian Model
Determination for Multivariate Ordinal and Binary Data.” *Computational
Statistics & Data Analysis* 52 (5): 2632–49.

</div>

<div id="ref-Williams2019">

Williams, Donald R. 2018. “Bayesian Estimation for Gaussian Graphical
Models: Structure Learning, Predictability, and Network Comparisons.”
*arXiv*. <https://doi.org/10.31234/OSF.IO/X8DPR>.

</div>

<div id="ref-Williams2019_bf">

Williams, Donald R, and Joris Mulder. 2019. “Bayesian Hypothesis Testing
for Gaussian Graphical Models: Conditional Independence and Order
Constraints.” *PsyArXiv*. <https://doi.org/10.31234/osf.io/ypxd8>.

</div>

<div id="ref-williams2020bggm">

———. 2020. “BGGM: Bayesian Gaussian Graphical Models in R.” *PsyArXiv*.

</div>

<div id="ref-williams2020comparing">

Williams, Donald R, Philippe Rast, Luis R Pericchi, and Joris Mulder.
2020. “Comparing Gaussian Graphical Models with the Posterior Predictive
Distribution and Bayesian Model Selection.” *Psychological Methods*.

</div>

</div>

<img src="readme_models/hex.jpg" width = 250 />

# BGGM: Bayesian Gaussian Graphical Models

[![CRAN
Version](http://www.r-pkg.org/badges/version/BGGM)](https://cran.r-project.org/package=BGGM)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BGGM)
[![Build
Status](https://travis-ci.com/donaldRwilliams/BGGM.svg?branch=master)](https://travis-ci.com/donaldRwilliams/BGGM)
[![status](https://joss.theoj.org/papers/4ecb84c5b3b2a2b5da46be4e0700502f/status.svg)](https://joss.theoj.org/papers/4ecb84c5b3b2a2b5da46be4e0700502f)

The `R` package **BGGM** provides tools for making Bayesian inference in
Gaussian graphical models (GGM). The methods are organized around two
general approaches for Bayesian inference: (1) estimation and (2)
hypothesis testing. The key distinction is that the former focuses on
either the posterior or posterior predictive distribution (Gelman, Meng,
and Stern 1996; see section 5 in Rubin 1984), whereas the latter focuses
on model comparison with the Bayes factor (Jeffreys 1961; Kass and
Raftery 1995).

## What is a Gaussian Graphical Model ?

A Gaussian graphical model captures conditional (in)dependencies among a
set of variables. These are pairwise relations (partial correlations)
controlling for the effects of all other variables in the model.

### Applications

The Gaussian graphical model is used across the sciences, including (but
not limited to) economics (Millington and Niranjan 2020), climate
science (Zerenner et al. 2014), genetics (Chu et al. 2009), and
psychology (Rodriguez et al. 2020).

## Installation

To install the latest release version (`2.0.0`) from CRAN use

``` r
install.packages("BGGM")    
```

The current developmental version can be installed with

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("donaldRwilliams/BGGM")
```

## Overview

The methods in **BGGM** build upon existing algorithms that are
well-known in the literature. The central contribution of **BGGM** is to
extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution
    (Mulder and Pericchi 2018)
    
      - [Estimation](#bayesian-estimation) (Williams 2018)

2.  Bayesian hypothesis testing with the matrix-F prior distribution
    (Williams and Mulder 2019)
    
      - [Exploratory hypothesis testing](#Exploratory)
    
      - [Confirmatory hypothesis testing](#Confirmatory)

3.  Comparing Gaussian graphical models (Williams 2018; Williams et al.
    2020)
    
      - [Partial correlation
        differences](#partial-correlation-differences)
    
      - [Posterior predictive check](#posterior-predictive-check)
    
      - [Exploratory hypothesis testing](#exploratory-groups)
    
      - [Confirmatory hypothesis testing](#confirmatory-groups)

4.  Extending inference beyond the conditional (in)dependence structure
    (Williams 2018)
    
      - [Predictability](#Predictability)
    
      - [Posterior uncertainty
        intervals](#partial-correlation-differences) for the partial
        correlations
    
      - [Custom Network Statistics](#custom-network-statistics)

The computationally intensive tasks are written in `c++` via the `R`
package **Rcpp** (Eddelbuettel et al. 2011) and the `c++` library
**Armadillo** (Sanderson and Curtin 2016). The Bayes factors are
computed with the `R` package **BFpack** (Mulder et al. 2019).
Furthermore, there are [plotting](#example-network-plot) functions for
each method, control variables can be included in the model (e.g., `~
gender`), and there is support for missing values (see `bggm_missing`).

## Supported Data Types

  - **Continuous**: The continuous method was described in Williams
    (2018). Note that this is based on the customary [Wishart
    distribution](https://en.wikipedia.org/wiki/Wishart_distribution).

  - **Binary**: The binary method builds directly upon Talhouk, Doucet,
    and Murphy (2012) that, in turn, built upon the approaches of
    Lawrence et al. (2008) and Webb and Forster (2008) (to name a few).

  - **Ordinal**: The ordinal methods require sampling thresholds. There
    are two approach included in **BGGM**. The customary approach
    described in Albert and Chib (1993) (the default) and the ‘Cowles’
    algorithm described in Cowles (1996).

  - **Mixed**: The mixed data (a combination of discrete and continuous)
    method was introduced in Hoff (2007). This is a semi-parametric
    copula model (i.e., a copula GGM) based on the ranked likelihood.
    Note that this can be used for *only* ordinal data (not restricted
    to “mixed” data).

## Illustrative Examples

The following includes brief examples for *some* of the methods in
**BGGM**.

### Bayesian Estimation

#### Posterior Sampling

An ordinal GGM is estimated with

``` r
library(BGGM)
library(ggplot2)

# data
Y <- ptsd[,1:5] + 1

# ordinal
fit <- estimate(Y, type = "ordinal", 
                analytic = FALSE)
```

Notice the `+ 1`. This is required, because the first category must be
`1` when `type = "ordinal"`. The partial correlations can the be
summarized with

``` r
summary(fit)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Analytic: FALSE 
#> Formula:  
#> Posterior Samples: 250 
#> Observations (n):
#> Nodes (p): 5 
#> Relations: 10 
#> --- 
#> Call: 
#> estimate(Y = Y, type = "ordinal", analytic = FALSE, iter = 250)
#> --- 
#> Estimates:
#>  Relation Post.mean Post.sd Cred.lb Cred.ub
#>    B1--B2     0.258   0.079   0.105   0.418
#>    B1--B3     0.028   0.086  -0.127   0.189
#>    B2--B3     0.517   0.058   0.406   0.616
#>    B1--B4     0.356   0.070   0.210   0.486
#>    B2--B4    -0.076   0.078  -0.219   0.063
#>    B3--B4     0.246   0.077   0.107   0.385
#>    B1--B5     0.131   0.080  -0.020   0.279
#>    B2--B5     0.127   0.083  -0.040   0.284
#>    B3--B5     0.202   0.079   0.063   0.366
#>    B4--B5     0.349   0.070   0.209   0.474
#> --- 
```

The returned object can also be plotted, which allows for visualizing
the posterior uncertainty interval for each relation. An example is
provided below in [Posterior uncertainty
intervals](#posterior-uncertatiny). The partial correlation matrix is
accessed with

``` r
pcor_mat(fit)
```

|    |    B1 |      B2 |    B3 |      B4 |    B5 |
| :- | ----: | ------: | ----: | ------: | ----: |
| B1 | 0.000 |   0.258 | 0.028 |   0.356 | 0.131 |
| B2 | 0.258 |   0.000 | 0.517 | \-0.076 | 0.127 |
| B3 | 0.028 |   0.517 | 0.000 |   0.246 | 0.202 |
| B4 | 0.356 | \-0.076 | 0.246 |   0.000 | 0.349 |
| B5 | 0.131 |   0.127 | 0.202 |   0.349 | 0.000 |

The graph is selected with

``` r
E <- select(fit)
```

and then plotted

``` r
# "communities"
comm <- substring(colnames(Y), 1, 1)

plot(select(fit), 
     groups = comm,
     edge_magnify = 5, 
     palette = "Pastel1", 
     node_size = 12)
```

![](readme_models/plt_est_net.png)

This basic “workflow” can be used with all methods and data types. A
more involved network plot is provided below.

#### Analytic

There is also an analytic solution that is based on the Wishart
distribution. This simple solution provides competitive performance with
“state-of-the-art” methods, assuming that *n* (observations) \> *p*
(variables). The one caveat is that it works only for `type =
"continuous"` (the default).

``` r
# analytic
fit <- estimate(Y, analytic = TRUE)

# network plot
plot(select(fit))
```

This is quite handy when (1) only the conditional dependence structure
is of interest and (2) an immediate solution is desirable. An example of
(2) is provided in [Posterior Predictive
Check](#posterior-predictive-check).

<br>

### Bayesian Hypothesis Testing

The Bayes factor based methods allow for determining the conditional
**in**dependence structure (evidence for the null hypothesis).

#### Exploratory

``` r
# now 10 nodes
Y <- ptsd[,1:10]

# exploratory hypothesis testing
fit<- explore(Y)

# select 
E <- select(fit, alternative = "exhaustive")
```

The option `alternative = "exhaustive"` compares three hypotheses: (1) a
null relation; (2) a positive relation; and (3) a negative relation.

``` r
summary(E)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Alternative: exhaustive 
#> --- 
#> Call:
#> select.explore(object = fit, alternative = "exhaustive")
#> --- 
#> Hypotheses: 
#> H0: rho = 0
#> H1: rho > 0
#> H2: rho < 0 
#> --- 
#> 
#>  Relation Post.mean Post.sd Pr.H0 Pr.H1 Pr.H2
#>  B1--B2    0.263    0.080   0.000 0.999 0.001
#>  B1--B3    0.020    0.081   0.710 0.173 0.116
#>  B2--B3    0.523    0.073   0.000 1.000 0.000
#>  B1--B4    0.362    0.070   0.000 1.000 0.000
#>  B2--B4   -0.082    0.068   0.459 0.061 0.480
#>  B3--B4    0.252    0.073   0.000 1.000 0.000
#>  B1--B5    0.129    0.072   0.120 0.847 0.033
#>  B2--B5    0.118    0.078   0.223 0.726 0.051
#>  B3--B5    0.213    0.077   0.001 0.996 0.003
#>  B4--B5    0.348    0.072   0.000 1.000 0.000
```

The posterior hypothesis probabilities are provided in the last three
columns. When using `plot(E)`, there is a network plot for each
hypothesis.

#### Confirmatory

A central contribution of **BGGM** is confirmatory hypothesis testing of
(in)equality constraints (Hoijtink 2011). By this we are referring to
testing expectations, as opposed to feeding the data to, say,
`estimate`, and seeing what happens to emerge.

In this example, the focus is on suicidal thoughts (`PHQ9`) in a
comorbidity network. Here is an example set of hypotheses

``` r
# data (+ 1)
Y <- depression_anxiety_t1 + 1

# example hypotheses
hyp <- c("PHQ2--PHQ9 > PHQ1--PHQ9 > 0; 
          PHQ2--PHQ9 = PHQ1--PHQ9 = 0")
```

There are two hypotheses separated by (`;`). The first expresses that
the relation `PHQ2--PHQ9` (“feeling down, depressed, or hopeless” and
“suicidal thoughts”) is larger than `PHQ1--PHQ9` (“little interest or
pleasure in doing things” and “suicidal thoughts”). In other words, that
the partial correlation is larger for `PHQ2--PHQ9`. There is an
additional constraint to positive values (`> 0`) for both relations. The
second hypothesis is then a “null” model.

``` r
# (try to) confirm
fit <- confirm(Y = Y, hypothesis = hyp, 
               type = "ordinal")
```

The object `fit` is then printed

``` r
fit

#> BGGM: Bayesian Gaussian Graphical Models 
#> Type: ordinal 
#> --- 
#> Posterior Samples: 250 
#> Observations (n): 403 
#> Variables (p): 16 
#> Delta: 15 
#> --- 
#> Call:
#> confirm(Y = Y + 1, hypothesis = hyp, type = "ordinal", 
#>     iter = 250)
#> --- 
#> Hypotheses: 
#> 
#> H1: PHQ2--PHQ9>PHQ1--PHQ9>0
#> H2: PHQ2--PHQ9=PHQ1--PHQ9=0
#> H3: complement
#> --- 
#> Posterior prob: 
#> 
#> p(H1|data) = 0.895
#> p(H2|data) = 0.002
#> p(H3|data) = 0.103
#> --- 
#> Bayes factor matrix: 
#>       H1      H2    H3
#> H1 1.000 529.910 8.666
#> H2 0.002   1.000 0.016
#> H3 0.115  61.147 1.000
#> --- 
#> note: equal hypothesis prior probabilities
```

The posterior hypothesis probability is `0.895` which provides some
evidence for the order constraint. The Bayes factor matrix then divides
the posterior probabilities. This provide a measure of *relative*
support for which hypothesis the data were more likely under.

Finally, the results can be plotted

``` r
plot(fit) + 
  scale_fill_brewer(palette = "Set2", 
                    name = "Posterior Prob") +
  ggtitle("Confirmatory: Comorbidity Network")
```

![](readme_models/confirm_hyp.png)

This demonstrates that all the `plot()` functions in **BGGM** return
`ggplot` objects that can be further customized. Note that **BGGM** is
not focused on making publication ready plots. Typically the bare
minimum is provided that can then be honed in.

<br>

### Comparing Gaussian Graphical Models

#### Partial Correlation Differences

This method compares groups by computing the difference for each
relation in the model. In other words, there are pairwise contrasts for
each partial correlation, resulting in a posterior distribution for each
difference.

In all examples in this section, personality networks are compared for
males and females.

``` r
# data
Y <- bfi

# males
Ymales <- subset(Y, gender == 1, 
                 select = -c(gender, education))

# females
Yfemales <- subset(Y, gender == 2, 
                 select = -c(gender, education))
```

Fit the model

``` r
fit <- ggm_compare_estimate(Ymales, Yfemales)
```

Then plot the results, in this case the posterior distribution for each
difference

``` r
# plot summary
plot(summary(fit))
```

![](readme_models/ggm_compare_estimate.png)

Note that it is also possible to use `select` for the object `fit` and
then plot the results. This produces a network plot including the
selected differences. Furthermore, it is also possible to plot the
partial correlations (not the differences). This is accomplished by
using `plot` with the summary computed from an `estimate` object ([see
above](#bayesian-estimation)).

#### Posterior Predictive Check

The predictive check method uses Jensen-Shannon divergence (i.e.,
symmetric Kullback-Leibler divergence
[Wikipedia](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence))
and the sum of squared error (for the partial correlation matrices) to
compare groups (Williams et al. 2020).

The following compares the groups

``` r
fit <- ggm_compare_ppc(Ymales, Yfemales)
```

Then print the summary output with

``` r
fit

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Test: Global Predictive Check 
#> Posterior Samples: 500 
#>   Group 1: 805 
#>   Group 2: 1631 
#> Nodes:  25 
#> Relations: 300 
#> --- 
#> Call: 
#> ggm_compare_ppc(Ymales, Yfemales, iter = 500)
#> --- 
#> Symmetric KL divergence (JSD): 
#>  
#>    contrast JSD.obs p_value
#>  Yg1 vs Yg2   0.442       0
#> --- 
#>  
#> Sum of Squared Error: 
#>  
#>    contrast SSE.obs p.value
#>  Yg1 vs Yg2   0.759       0
#> --- 
#> note:
#> JSD is Jensen-Shannon divergence 
```

In this case, there seems to be decisive evidence that the networks are
different (as indicated by the posterior predictive *p*-value). The
predictive distribution can also be plotted

``` r
plot(fit, 
     critical = 0.05)$plot_jsd
```

![](readme_models/ppc_1.png)

where the red region is the “critical” area and the black point is the
observed KL divergence for the networks. This again shows that the
“distance” between the networks is much more than expected, assuming
that the groups were actually the same.

This next example is a new feature in **BGGM** (`2.0.0`), that allows
for comparing GGMs any way the user wants. All that is required is to
(1) decide on a test-statistic and (2) write a custom function.

Here is an example using Hamming distance
([Wikipedia](https://en.wikipedia.org/wiki/Hamming_distance)), which is
essentially the squared error between adjacency matrices (a test for
different structures).

First define the custom function

``` r
f <- function(Yg1, Yg2){

# remove NA
x <- na.omit(Yg1)
y <- na.omit(Yg2)

# nodes
p <- ncol(x)

# identity matrix
I_p <- diag(p)

# estimate graphs
fit1 <-  estimate(x, analytic = TRUE)
fit2 <-  estimate(y, analytic = TRUE)

# select graphs
sel1 <- select(fit1)
sel2 <- select(fit2)

# Hamming distance
sum((sel1$adj[upper.tri(I_p)] - sel2$adj[upper.tri(I_p)])^2)
}
```

Note that (1) `analytic = TRUE` is being used, which is needed in this
case because two graphs are estimated for each iteration (or draw from
the posterior predictive distribution) and (2) `f` requires two datasets
as the input and returns a single number (the chosen test-statistic).
The next step is to compute the observed Hamming distance

``` r
# observed difference
obs <- f(Ymales, Yfemales)
```

then compare the groups

``` r
fit <- ggm_compare_ppc(Ymales, Yfemales,
                             FUN = f,
                             custom_obs  = obs)

# print
fit

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Test: Global Predictive Check 
#> Posterior Samples: 250 
#>   Group 1: 805 
#>   Group 2: 1631 
#> Nodes:  25 
#> Relations: 300 
#> --- 
#> Call: 
#> ggm_compare_ppc(Ymales, Yfemales, iter = 250, FUN = f, custom_obs = obs)
#> --- 
#> Custom: 
#>  
#>    contrast custom.obs p.value
#>  Yg1 vs Yg2         75   0.576
#> --- 
```

In this case, the *p*-value does not indicate that the groups are
different for this test-statistic. This may seem contradictory to the
previous results, but it is important to note that Hamming distance asks
a much different question related to the adjacency matrices (no other
information, such as edge weights, is considered).

#### Exploratory (groups)

The Bayes factor based methods allow for determining the conditional
**in**dependence structure (evidence for the null hypothesis), in this
case for group equality.

Fit the model

``` r
fit <- ggm_compare_explore(Ymales, Yfemales)
```

Then plot the results

``` r
plot(summary(fit)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_blank()) 
```

![](readme_models/plt_ggm_compare_explore.png)

Here the posterior probability for a difference is visualized for each
relation in the GGM. Note that it is also possible to use `select` for
the object `fit` and then plot the results. This produces a network plot
including the selected differences, in addition to a plot depicting the
relations for which there was evidence for the null hypothesis.

#### Confirmatory (groups)

A central contribution of **BGGM** is confirmatory hypothesis testing of
(in)equality constraints (Hoijtink 2011), in this case for comparing
groups. By this we are referring to testing expectations, as opposed to
feeding the data to, say, `estimate`, and seeing what happens to emerge.

In this example, the focus is on agreeableness in a personality network.
Here is a set of hypotheses

``` r
hyp <- c("g1_A2--A4 > g2_A2--A4 > 0 & g1_A4--A5 > g2_A4--A5 > 0;
          g1_A4--A5 = g2_A4--A5 = 0  & g1_A2--A4 = g2_A2--A4 = 0")
```

where the variables are `A2` (“inquire about others’ well being”), `A4`
(“love children”), and `A5` (“make people feel at ease”). The first
hypothesis states that the conditionally dependent effects are larger
for female than males (note the `&`), with the additional constraint to
positive values, whereas the second hypothesis is a “null” model.

The hypothesis is tested with the following

``` r
fit <- ggm_compare_confirm(Yfemales, Ymales, 
                           hypothesis = hyp)

# print
fit

#> BGGM: Bayesian Gaussian Graphical Models
#> Type: continuous
#> ---
#> Posterior Samples: 500
#>   Group 1: 1631
#>   Group 2: 805
#> Variables (p): 25
#> Relations: 300
#> Delta: 15
#> ---
#> Call:
#> ggm_compare_confirm(Yfemales, Ymales, hypothesis = hyp, iter = 500)
#> ---
#> Hypotheses:
#> 
#> H1: g1_A2--A4>g2_A2--A4>0&g1_A4--A5>g2_A4--A5>0
#> H2: g1_A4--A5=g2_A4--A5=0&g1_A2--A4=g2_A2--A4=0
#> H3: complement
#> ---
#> Posterior prob:
#> 
#> p(H1|data) = 0.989
#> p(H2|data) = 0
#> p(H3|data) = 0.011
#> ---
#> Bayes factor matrix:
#>       H1           H2     H3
#> H1 1.000 1.180798e+14 92.115
#> H2 0.000 1.000000e+00  0.000
#> H3 0.011 1.281873e+12  1.000
#> ---
#> note: equal hypothesis prior probabilities
```

The posterior hypothesis probability is 0.989 which provides strong
evidence for the hypothesis that predicted these “agreeableness”
relations would be larger in females than in males. This can also be
plotted, as in [Confirmatory (one group)](#confirmatory). See Rodriguez
et al. (2020) for a full treatment of confirmatory testing in
substantive applications.

### Beyond the Conditional (In)dependence Structure

#### Predictability

In this example, predictability is computed for each node in the network
(see here for rationale Haslbeck and Waldorp 2018). Currently **BGGM**
computes Bayesian variance explained for all data types (Gelman et al.
2019).

The following computes predictability for binary data

``` r
# binary
Y <- women_math

# fit model
fit <- estimate(Y, type = "binary")

# compute r2
r2 <- predictability(fit, iter = 500)

# plot
plot(r2, type = "ridgeline")
```

![](readme_models/predictability.png)

#### Posterior Uncertainty

See [Partial Correlation Differences](#partial-correlation-differences)

#### Custom Network Statistics

A new feature to **BGGM** allows for computing user defined network
statistics, given a partial correlation or weighted adjacency matrix.

Here is an example for bridge centrality (Jones, Ma, and McNally 2019).
The first step is to define the function

``` r
# need this package 
library(networktools)

# custom function
f <- function(x, ...){
 bridge(x, ...)$`Bridge Strength`
}
```

Note that `x` takes the matrix and `f` can return either a single number
or a number for each node. The next step is to fit the model and compute
the network statistic

``` r
# data
Y <- ptsd

# clusters
communities <- substring(colnames(Y), 1, 1)

# estimate the model
fit <- estimate(Y)

# bridge strength
net_stat <- roll_your_own(fit,
                          FUN = f,
                          select = TRUE,
                          communities = communities)
```

The function `f` is provided to `FUN` and `communities` is passed to
`brigde` (inside of `f`) via `...`. The results can be printed

``` r
# print
net_stat

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Network Stats: Roll Your Own
#> Posterior Samples: 100 
#> --- 
#> Estimates: 
#> 
#>  Node Post.mean Post.sd Cred.lb Cred.ub
#>     1     0.340   0.097   0.166   0.546
#>     2     0.319   0.100   0.176   0.513
#>     3     0.000   0.000   0.000   0.000
#>     4     0.337   0.086   0.189   0.489
#>     5     0.559   0.133   0.332   0.791
#>     6     0.188   0.073   0.029   0.320
#>     7     0.505   0.138   0.241   0.781
#>     8     0.153   0.070   0.022   0.286
#>     9     0.175   0.063   0.041   0.281
#>    10     0.000   0.000   0.000   0.000
#>    11     0.365   0.107   0.178   0.627
#>    12     0.479   0.093   0.280   0.637
#>    13     0.155   0.074   0.022   0.301
#>    14     0.000   0.000   0.000   0.000
#>    15     0.374   0.097   0.175   0.550
#>    16     0.174   0.065   0.034   0.295
#>    17     0.000   0.000   0.000   0.000
#>    18     0.491   0.132   0.238   0.745
#>    19     0.613   0.113   0.408   0.825
#>    20     0.144   0.066   0.038   0.289
#> --- 
```

And then plotted

``` r
plot(net_stat)
```

![](readme_models/bridge.png)

There are additional examples in the documentation.

### Example Network Plot

Here is an example of a more involved network plot. In this case, the
graph is estimated with a semi-parametric copula (`type = "mixed"`),
where two control variables are included in the model.

``` r
# personality (includes gender and education)
Y <- bfi

# fit copula GGM
fit <- estimate(Y, type = "mixed")

# select graph
E <- select(fit)
```

The graph is then plotted

``` r
# extract communities
comm <- substring(colnames(Y), 1, 1)

# plot
plot(E, 
     # enlarge edges
     edge_magnify = 5, 
     # cluster nodes
     groups = comm, 
     # change layout
     layout = "circle")$plt +
  # plot title
  ggtitle("Semi-Parametric Copula") +
  # add custom labels
  scale_color_brewer(breaks = c("A", "C", 
                                "E", "N", 
                                "O", "e",  
                                "g"), 
                     labels =   c("A", "C", 
                                 "E", "N", 
                                 "O",  
                                 "Education",   
                                 "Gender"), 
                     palette = "Set2")
```

![](readme_models/plt_net_example.png)

Note that `layout` can be changed to any option provided in the `R`
package **sna** (Butts 2019).

## Additional Features

The primary focus of **BGGM** is Gaussian graphical modeling (the
inverse covariance matrix). The residue is a suite of useful methods not
explicitly for GGMs. For example,

### Bivariate Correlations

Bivariate correlations for `binary` (tetrachoric), `ordinal`
(polychoric), `mixed` (rank based), and `continuous` (Pearson’s) data.

Here is an example for computing tetrachoric correlations:

``` r
# binary data
Y <- women_math[1:500,]

cors <- zero_order_cors(Y, type = "binary", iter = 250)

cors$R
```

|   |       1 |       2 |       3 |       4 |       5 |       6 |
| :- | ------: | ------: | ------: | ------: | ------: | ------: |
| 1 |   1.000 | \-0.198 |   0.506 |   0.122 | \-0.140 |   0.098 |
| 2 | \-0.198 |   1.000 | \-0.482 | \-0.013 | \-0.146 | \-0.146 |
| 3 |   0.506 | \-0.482 |   1.000 |   0.310 | \-0.343 |   0.351 |
| 4 |   0.122 | \-0.013 |   0.310 |   1.000 | \-0.363 |   0.169 |
| 5 | \-0.140 | \-0.146 | \-0.343 | \-0.363 |   1.000 | \-0.194 |
| 6 |   0.098 | \-0.146 |   0.351 |   0.169 | \-0.194 |   1.000 |

The object `cors` also includes the sampled correlation matrices (in
this case 250) in an array.

### Multivariate Regression

Multivariate regression for binary (probit), ordinal (probit), mixed
(rank likelihood), and continuous data.

Here is an example for a multivariate probit model with an ordinal
outcome, where `E5` (“take charge”) and `N5` (“panic easily”) are
predicted by `gender` and `education`:

``` r
# personality data
Y <- bfi

# variables
Y <- subset(Y, select = c("E5", "N5", 
                          "gender", "education"))


mv_probit <- estimate(Y, formula = ~ gender + as.factor(education), 
                      type = "ordinal")
```

Note that **BGGM** does not use the customary `model.matrix`
formulation. This is for good reason, as each variable in the GGM does
not need to be written out. Here we effectively “tricked” **BGGM** to
fit a multivariate probit model (each variable included in `formula` is
removed from `Y`).

``` r
regression_summary(mv_probit)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Formula: ~ gender + as.factor(education) 
#> --- 
#> Coefficients: 
#>  
#> E5 
#>                       Post.mean Post.sd Cred.lb Cred.ub
#> (Intercept)               1.852   0.533   1.049   3.142
#> gender                    0.169   0.066   0.065   0.295
#> as.factor(education)2     0.215   0.109   0.024   0.437
#> as.factor(education)3     0.271   0.104   0.089   0.445
#> as.factor(education)4     0.206   0.103   0.019   0.404
#> as.factor(education)5     0.345   0.128   0.120   0.593
#> --- 
#> N5 
#>                       Post.mean Post.sd Cred.lb Cred.ub
#> (Intercept)               0.210   0.114  -0.012   0.434
#> gender                    0.502   0.140   0.291   0.835
#> as.factor(education)2    -0.127   0.103  -0.345   0.058
#> as.factor(education)3    -0.104   0.081  -0.258   0.034
#> as.factor(education)4    -0.218   0.104  -0.427  -0.024
#> as.factor(education)5    -0.229   0.103  -0.449  -0.038
#> --- 
#> Residual Correlation Matrix: 
#>       E5    N5
#> E5  1.00 -0.18
#> N5 -0.18  1.00
#> ---
```

This basic idea can also be used to fit regression models with a single
outcome.

## Note on Conditional (In)dependence Models for Latent Data

All of the data types (besides continuous) model latent data. That is,
unobserved data that is assumed to be Gaussian distributed. For example,
a tetrachoric correlation (binary data) is a special case of a
polychoric correlation (ordinal data). Both relations are between
"theorized normally distributed continuous *latent* variables
[Wikepedia](https://en.wikipedia.org/wiki/Polychoric_correlation). In
both instances, the corresponding partial correlation between observed
variables is conditioned on the remaining variables in the *latent*
space. This implies that interpretation is similar to continuous data,
but with respect to latent variables. We refer interested users to (see
page 2364, section 2.2, in Webb and Forster 2008).

## High Dimensional Data?

**BGGM** was built specifically for social-behavioral scientists. Of
course, the methods can be used by all researchers. However, there is
currently *not* support for high-dimensional data (i.e., more variables
than observations) that are common place in, say, the genetics
literature. These data are rare in the social-behavioral sciences. In
the future, support for high-dimensional data may be added to **BGGM**.

## Bug Reports, Feature Requests, and Contributing

Bug reports and feature requests can be made by opening an issue on
[Github](https://github.com/donaldRwilliams/BGGM/issues). To contribute
towards the development of **BGGM**, you can start a branch with a pull
request and we can discuss the proposed changes there.

## Comparison to Other Software

**BGGM** is the only `R` package to implement all of these algorithms
and methods. The `mixed` data approach is also implemented in the
package **sbgcop** (base `R`, Hoff 2007). The `R` package **BDgraph**
implements a Gaussian copula graphical model in `c++` (Mohammadi and Wit
2015), but not the binary or ordinal approaches. Furthermore, **BGGM**
is the only package for confirmatory testing and comparing graphical
models with the methods described in Williams et al. (2020).

## References

<div id="refs" class="references">

<div id="ref-albert1993bayesian">

Albert, James H, and Siddhartha Chib. 1993. “Bayesian Analysis of Binary
and Polychotomous Response Data.” *Journal of the American Statistical
Association* 88 (422): 669–79.
<https://doi.org/10.1080/01621459.1993.10476321>.

</div>

<div id="ref-sna">

Butts, Carter T. 2019. *Sna: Tools for Social Network Analysis*.
<https://CRAN.R-project.org/package=sna>.

</div>

<div id="ref-chu2009graphical">

Chu, Jen-hwa, Scott T Weiss, Vincent J Carey, and Benjamin A Raby. 2009.
“A Graphical Model Approach for Inferring Large-Scale Networks
Integrating Gene Expression and Genetic Polymorphism.” *BMC Systems
Biology* 3 (1): 55. <https://doi.org/10.1186/1752-0509-3-55>.

</div>

<div id="ref-cowles1996accelerating">

Cowles, Mary Kathryn. 1996. “Accelerating Monte Carlo Markov Chain
Convergence for Cumulative-Link Generalized Linear Models.” *Statistics
and Computing* 6 (2): 101–11. <https://doi.org/10.1007/bf00162520>.

</div>

<div id="ref-eddelbuettel2011rcpp">

Eddelbuettel, Dirk, Romain François, J Allaire, Kevin Ushey, Qiang Kou,
N Russel, John Chambers, and D Bates. 2011. “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software* 40 (8): 1–18.

</div>

<div id="ref-gelman_r2_2019">

Gelman, Andrew, Ben Goodrich, Jonah Gabry, and Aki Vehtari. 2019.
“R-squared for Bayesian Regression Models.” *American Statistician* 73
(3): 307–9. <https://doi.org/10.1080/00031305.2018.1549100>.

</div>

<div id="ref-Gelman1996a">

Gelman, Andrew, Xiao-Li Meng, and Hal Stern. 1996. “Posterior Predictive
Assessment of Model Fitness via Realized Discrepancies.” *Statistica
Sinica* 6 (4): 733–807.

</div>

<div id="ref-haslbeck2018well">

Haslbeck, Jonas MB, and Lourens J Waldorp. 2018. “How Well Do Network
Models Predict Observations? On the Importance of Predictability in
Network Models.” *Behavior Research Methods* 50 (2): 853–61.
<https://doi.org/10.3758/s13428-017-0910-x>.

</div>

<div id="ref-hoff2007extending">

Hoff, Peter D. 2007. “Extending the Rank Likelihood for Semiparametric
Copula Estimation.” *The Annals of Applied Statistics* 1 (1): 265–83.
<https://doi.org/10.1214/07-AOAS107>.

</div>

<div id="ref-Hoijtink2011">

Hoijtink, H. 2011. *Informative hypotheses: Theory and practice for
behavioral and social scientists*. Chapman; Hall/CRC.

</div>

<div id="ref-Jeffreys1961">

Jeffreys, Harold. 1961. *The theory of probability*. Oxford: Oxford
University Press.

</div>

<div id="ref-jones2019bridge">

Jones, Payton J, Ruofan Ma, and Richard J McNally. 2019. “Bridge
Centrality: A Network Approach to Understanding Comorbidity.”
*Multivariate Behavioral Research*, 1–15.
<https://doi.org/10.1080/00273171.2019.1614898>.

</div>

<div id="ref-Kass1995">

Kass, Robert E, and Adrian E Raftery. 1995. “Bayes Factors.” *Journal of
the American Statistical Association* 90 (430): 773–95.

</div>

<div id="ref-lawrence2008bayesian">

Lawrence, Earl, Derek Bingham, Chuanhai Liu, and Vijayan N Nair. 2008.
“Bayesian Inference for Multivariate Ordinal Data Using Parameter
Expansion.” *Technometrics* 50 (2): 182–91.
<https://doi.org/10.1198/004017008000000064>.

</div>

<div id="ref-millington2020partial">

Millington, Tristan, and Mahesan Niranjan. 2020. “Partial Correlation
Financial Networks.” *Applied Network Science* 5 (1): 11.
<https://doi.org/10.1007/s41109-020-0251-z>.

</div>

<div id="ref-mohammadi2015bdgraph">

Mohammadi, Reza, and Ernst C Wit. 2015. “BDgraph: An R Package for
Bayesian Structure Learning in Graphical Models.” *Journal of
Statistical Software* 89 (3). <https://doi.org/10.18637/jss.v089.i03>.

</div>

<div id="ref-mulder2019bfpack">

Mulder, Joris, Xin Gu, Anton Olsson-Collentine, Andrew Tomarken, Florian
Böing-Messing, Herbert Hoijtink, Marlyne Meijerink, et al. 2019.
“BFpack: Flexible Bayes Factor Testing of Scientific Theories in R.”
*arXiv Preprint arXiv:1911.07728*.

</div>

<div id="ref-Mulder2018">

Mulder, Joris, and Luis Pericchi. 2018. “The Matrix-F Prior for
Estimating and Testing Covariance Matrices.” *Bayesian Analysis*, no. 4:
1–22. <https://doi.org/10.1214/17-BA1092>.

</div>

<div id="ref-rodriguez2020formalizing">

Rodriguez, Josue E, Donald R Williams, Philippe Rast, and Joris Mulder.
2020. “On Formalizing Theoretical Expectations: Bayesian Testing of
Central Structures in Psychological Networks.” *PsyArXiv*.
<https://doi.org/10.31234/osf.io/zw7pf>.

</div>

<div id="ref-rubin1984bayesianly">

Rubin, Donald B. 1984. “Bayesianly Justifiable and Relevant Frequency
Calculations for the Applied Statistician.” *The Annals of Statistics*,
1151–72. <https://doi.org/10.1214/aos/1176346785>.

</div>

<div id="ref-sanderson2016armadillo">

Sanderson, Conrad, and Ryan Curtin. 2016. “Armadillo: A Template-Based
C++ Library for Linear Algebra.” *Journal of Open Source Software* 1
(2): 26. <https://doi.org/10.21105/joss.00026>.

</div>

<div id="ref-talhouk2012efficient">

Talhouk, Aline, Arnaud Doucet, and Kevin Murphy. 2012. “Efficient
Bayesian Inference for Multivariate Probit Models with Sparse Inverse
Correlation Matrices.” *Journal of Computational and Graphical
Statistics* 21 (3): 739–57.
<https://doi.org/10.1080/10618600.2012.679239>.

</div>

<div id="ref-webb2008bayesian">

Webb, Emily L, and Jonathan J Forster. 2008. “Bayesian Model
Determination for Multivariate Ordinal and Binary Data.” *Computational
Statistics & Data Analysis* 52 (5): 2632–49.
<https://doi.org/10.1016/j.csda.2007.09.008>.

</div>

<div id="ref-Williams2019">

Williams, Donald R. 2018. “Bayesian Estimation for Gaussian Graphical
Models: Structure Learning, Predictability, and Network Comparisons.”
*arXiv*. <https://doi.org/10.31234/OSF.IO/X8DPR>.

</div>

<div id="ref-Williams2019_bf">

Williams, Donald R, and Joris Mulder. 2019. “Bayesian Hypothesis Testing
for Gaussian Graphical Models: Conditional Independence and Order
Constraints.” *PsyArXiv*. <https://doi.org/10.31234/osf.io/ypxd8>.

</div>

<div id="ref-williams2020comparing">

Williams, Donald R, Philippe Rast, Luis R Pericchi, and Joris Mulder.
2020. “Comparing Gaussian Graphical Models with the Posterior Predictive
Distribution and Bayesian Model Selection.” *Psychological Methods*.
<https://doi.org/10.1037/met0000254>.

</div>

<div id="ref-zerenner2014gaussian">

Zerenner, Tanja, Petra Friederichs, Klaus Lehnertz, and Andreas Hense.
2014. “A Gaussian Graphical Model Approach to Climate Networks.” *Chaos:
An Interdisciplinary Journal of Nonlinear Science* 24 (2): 023103.
<https://doi.org/10.1063/1.4870402>.

</div>

</div>
# BGGM 2.0.1
This version of BGGM included changes based on the JOSS reviews: see [here](https://github.com/openjournals/joss-reviews/issues/2111) for 
the overview and [here](https://github.com/donaldRwilliams/BGGM/issues?q=is%3Aissue+is%3Aclosed) for specific issues.


# BGGM 2.0.0

**BGGM** was almost completely rewritten for version `2.0.0`. This was due to adding support 
for binary, ordinal, and mixed data, which required that the methods be written in `c ++`. 
Unfortunately, as a result, lots of code from version `1.0.0` is broken.

## Added features

* Full support for binary, ordinal, and mixed data. This is implemented with the argument `type`

* `roll_your_own`: compute custom network statistics from a weighted adjacency matrix or a partial 
correlation matrix

* `pcor_to_cor`: convert the sampled partial correlation matrices into correlation matrices. 

* `zero_order_cors`: compute zero order correlations 

* `convergence`: acf and trace plots

* `posterior_samples`: extract posterior samples

* `regression_summary`: summarize multivariate regression

* `pcor_sum`: Compute and compare partial correlation sums

* `weighted_adj_mat`: Extract the Weighted Adjacency Matrix

* `pcor_mat`: 	Extract the Partial Correlation Matrix

* Five additional data sets were added.

## Extensions
* `ggm_compare_ppc`: added option for custom network statistics

* Added option to control for variables with `formula`

* A progress bar was added to many functions


# BGGM 1.0.0

Initial CRAN release
---
output: github_document
bibliography: inst/REFERENCES.bib
---

```{r, echo = FALSE, message=F}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  dev = "png",
  dpi = 200,
  fig.align = "center",
  knitr::opts_chunk$set(comment = NA)
  )
library(ggplot2)
library(BGGM)
```

<img src="readme_models/hex.jpg" width = 250 />

# BGGM: Bayesian Gaussian Graphical Models

[![CRAN Version](http://www.r-pkg.org/badges/version/BGGM)](https://cran.r-project.org/package=BGGM)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BGGM)
[![Build Status](https://travis-ci.com/donaldRwilliams/BGGM.svg?branch=master)](https://travis-ci.com/donaldRwilliams/BGGM)
[![status](https://joss.theoj.org/papers/4ecb84c5b3b2a2b5da46be4e0700502f/status.svg)](https://joss.theoj.org/papers/4ecb84c5b3b2a2b5da46be4e0700502f)

The `R` package **BGGM** provides tools for making Bayesian inference in 
Gaussian graphical models (GGM). The methods are organized around two general 
approaches for Bayesian inference: (1) estimation and (2) hypothesis 
testing. The key distinction is that the former focuses on either 
the posterior or posterior predictive distribution [@Gelman1996a; see section 5 in @rubin1984bayesianly], whereas the latter focuses on model comparison with the Bayes factor [@Jeffreys1961; @Kass1995].

## What is a Gaussian Graphical Model ?
A Gaussian graphical model captures conditional (in)dependencies among a set 
of variables. These are pairwise relations (partial correlations) controlling for 
the effects of all other variables in the model.

### Applications
The Gaussian graphical model is used across the sciences, including 
(but not limited to) economics [@millington2020partial], climate science 
[@zerenner2014gaussian], genetics [@chu2009graphical], and psychology [@rodriguez2020formalizing]. 


## Installation

To install the latest release version (`2.0.0`) from CRAN use
```{r gh-installation, eval = FALSE}	
install.packages("BGGM")	
```

The current developmental version can be installed with	

```{r, eval = FALSE}	
if (!requireNamespace("remotes")) {	
  install.packages("remotes")	
}	
remotes::install_github("donaldRwilliams/BGGM")
```

## Overview
The methods in **BGGM** build upon existing algorithms that are well-known in the literature.
The central contribution of **BGGM** is to extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution [@Mulder2018]
  
    + [Estimation](#bayesian-estimation) [@Williams2019]

2. Bayesian hypothesis testing with the matrix-F prior distribution [@Williams2019_bf]

    + [Exploratory hypothesis testing](#Exploratory)
  
    + [Confirmatory hypothesis testing](#Confirmatory)
    
3. Comparing Gaussian graphical models [@Williams2019; @williams2020comparing]
    
    + [Partial correlation differences](#partial-correlation-differences) 
    
    + [Posterior predictive check](#posterior-predictive-check)
    
    + [Exploratory hypothesis testing](#exploratory-groups) 
    
    + [Confirmatory hypothesis testing](#confirmatory-groups)

4. Extending inference beyond the conditional (in)dependence structure [@Williams2019]

    +  [Predictability](#Predictability) 
    
    +  [Posterior uncertainty intervals](#partial-correlation-differences) for the 
       partial correlations
       
    +  [Custom Network Statistics](#custom-network-statistics)
    
    
The computationally intensive tasks are written in `c++` via the `R` package **Rcpp** [@eddelbuettel2011rcpp] and the `c++` library **Armadillo** [@sanderson2016armadillo]. The Bayes factors are computed with the `R` package **BFpack** [@mulder2019bfpack]. Furthermore, there are [plotting](#example-network-plot) functions
for each method, control variables can be included in the model (e.g., `~ gender`), 
and there is support for missing values (see `bggm_missing`).

## Supported Data Types

* **Continuous**: The continuous method was described in  @Williams2019. Note that 
                  this is based on the customary [Wishart distribution](https://en.wikipedia.org/wiki/Wishart_distribution).

* **Binary**: The binary method builds directly upon @talhouk2012efficient
  that, in turn, built upon the approaches of @lawrence2008bayesian and
  @webb2008bayesian (to name a few).
  
* **Ordinal**: The ordinal methods require sampling thresholds. There are two approach 
   included in **BGGM**. The customary approach described in @albert1993bayesian 
   (the default) and the 'Cowles' algorithm described in @cowles1996accelerating.
   
* **Mixed**: The mixed data (a combination of discrete and continuous) method was introduced
 in @hoff2007extending. This is a semi-parametric copula model
 (i.e., a copula GGM) based on the ranked likelihood. Note that this can be used for 
 *only* ordinal data (not restricted to "mixed" data).

## Illustrative Examples
The following includes brief examples for *some* of the methods in **BGGM**.

### Bayesian Estimation

#### Posterior Sampling
An ordinal GGM is estimated with

```{r, eval = FALSE}
library(BGGM)
library(ggplot2)

# data
Y <- ptsd[,1:5] + 1

# ordinal
fit <- estimate(Y, type = "ordinal", 
                analytic = FALSE)
```

Notice the `+ 1`. This is required, because the first category must be `1` when `type = "ordinal"`. The partial correlations can the be summarized with

```{r, eval = FALSE}
summary(fit)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Analytic: FALSE 
#> Formula:  
#> Posterior Samples: 250 
#> Observations (n):
#> Nodes (p): 5 
#> Relations: 10 
#> --- 
#> Call: 
#> estimate(Y = Y, type = "ordinal", analytic = FALSE, iter = 250)
#> --- 
#> Estimates:
#>  Relation Post.mean Post.sd Cred.lb Cred.ub
#>    B1--B2     0.258   0.079   0.105   0.418
#>    B1--B3     0.028   0.086  -0.127   0.189
#>    B2--B3     0.517   0.058   0.406   0.616
#>    B1--B4     0.356   0.070   0.210   0.486
#>    B2--B4    -0.076   0.078  -0.219   0.063
#>    B3--B4     0.246   0.077   0.107   0.385
#>    B1--B5     0.131   0.080  -0.020   0.279
#>    B2--B5     0.127   0.083  -0.040   0.284
#>    B3--B5     0.202   0.079   0.063   0.366
#>    B4--B5     0.349   0.070   0.209   0.474
#> --- 
```

The returned object can also be plotted, which allows for visualizing the posterior uncertainty interval for each relation. An example is provided below in [Posterior uncertainty intervals](#posterior-uncertatiny). 
The partial correlation matrix is accessed with 

```{r, eval=FALSE}
pcor_mat(fit)
```

```{r, echo = FALSE, results='asis'}
load(file = "readme_models/pcor_ex.rda")
knitr::kable(pcor_ex, row.names = TRUE)
```

The graph is selected with
```{r, eval=FALSE}
E <- select(fit)
```

and then plotted
```{r, eval = FALSE}
# "communities"
comm <- substring(colnames(Y), 1, 1)

plot(select(fit), 
     groups = comm,
     edge_magnify = 5, 
     palette = "Pastel1", 
     node_size = 12)
```

![](readme_models/plt_est_net.png)

This basic "workflow" can be used with all methods and data types. A more involved network plot
is provided below.

#### Analytic
There is also an analytic solution that is based on the Wishart distribution. This 
simple solution provides competitive performance with "state-of-the-art" methods, 
assuming that *n* (observations) > *p* (variables). The one 
caveat is that it works only for `type = "continuous"` (the default).

```{r, eval=FALSE}
# analytic
fit <- estimate(Y, analytic = TRUE)

# network plot
plot(select(fit))
```

This is quite handy when (1) only the conditional dependence structure is of interest and (2) an immediate solution is desirable. An example of (2) is provided  in [Posterior Predictive Check](#posterior-predictive-check).

<br>

### Bayesian Hypothesis Testing
The Bayes factor based methods allow for determining the conditional 
**in**dependence structure (evidence for the null hypothesis).

#### Exploratory

```{r, eval=FALSE}
# now 10 nodes
Y <- ptsd[,1:10]

# exploratory hypothesis testing
fit<- explore(Y)

# select 
E <- select(fit, alternative = "exhaustive")
```

The option `alternative = "exhaustive"` compares three hypotheses: (1) a null relation; (2) a positive
relation; and (3) a negative relation. 
```{r, eval = FALSE}
summary(E)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Alternative: exhaustive 
#> --- 
#> Call:
#> select.explore(object = fit, alternative = "exhaustive")
#> --- 
#> Hypotheses: 
#> H0: rho = 0
#> H1: rho > 0
#> H2: rho < 0 
#> --- 
#> 
#>  Relation Post.mean Post.sd Pr.H0 Pr.H1 Pr.H2
#>  B1--B2    0.263    0.080   0.000 0.999 0.001
#>  B1--B3    0.020    0.081   0.710 0.173 0.116
#>  B2--B3    0.523    0.073   0.000 1.000 0.000
#>  B1--B4    0.362    0.070   0.000 1.000 0.000
#>  B2--B4   -0.082    0.068   0.459 0.061 0.480
#>  B3--B4    0.252    0.073   0.000 1.000 0.000
#>  B1--B5    0.129    0.072   0.120 0.847 0.033
#>  B2--B5    0.118    0.078   0.223 0.726 0.051
#>  B3--B5    0.213    0.077   0.001 0.996 0.003
#>  B4--B5    0.348    0.072   0.000 1.000 0.000
```

The posterior hypothesis probabilities are provided in the last three columns.
When using `plot(E)`, there is a network plot for each hypothesis.


#### Confirmatory
A central contribution of **BGGM** is confirmatory hypothesis testing of (in)equality constraints [@Hoijtink2011]. By this we are referring to testing expectations, as opposed to feeding the data to, say, `estimate`, and seeing what happens to emerge. 

In this example, the focus is on suicidal thoughts (`PHQ9`) in a comorbidity network. Here is an example set of hypotheses

```{r}
# data (+ 1)
Y <- depression_anxiety_t1 + 1

# example hypotheses
hyp <- c("PHQ2--PHQ9 > PHQ1--PHQ9 > 0; 
          PHQ2--PHQ9 = PHQ1--PHQ9 = 0")
```

There are two hypotheses separated by (`;`). The first expresses that the relation `PHQ2--PHQ9` ("feeling down, depressed, or hopeless" and "suicidal thoughts") is larger than `PHQ1--PHQ9` ("little interest or pleasure in doing things" and "suicidal thoughts"). In other words, that the partial correlation is larger for `PHQ2--PHQ9`. There is an additional constraint to positive values (`> 0`) for both relations. The second hypothesis is then a "null" model.


```{r, echo=FALSE, warning = FALSE}
load(file = "readme_models/fit_hyp1.rda")
fit <- fit_hyp1
```


```{r, eval=FALSE}
# (try to) confirm
fit <- confirm(Y = Y, hypothesis = hyp, 
               type = "ordinal")
```

The object `fit` is then printed
```{r, eval=FALSE}
fit

#> BGGM: Bayesian Gaussian Graphical Models 
#> Type: ordinal 
#> --- 
#> Posterior Samples: 250 
#> Observations (n): 403 
#> Variables (p): 16 
#> Delta: 15 
#> --- 
#> Call:
#> confirm(Y = Y + 1, hypothesis = hyp, type = "ordinal", 
#>     iter = 250)
#> --- 
#> Hypotheses: 
#> 
#> H1: PHQ2--PHQ9>PHQ1--PHQ9>0
#> H2: PHQ2--PHQ9=PHQ1--PHQ9=0
#> H3: complement
#> --- 
#> Posterior prob: 
#> 
#> p(H1|data) = 0.895
#> p(H2|data) = 0.002
#> p(H3|data) = 0.103
#> --- 
#> Bayes factor matrix: 
#>       H1      H2    H3
#> H1 1.000 529.910 8.666
#> H2 0.002   1.000 0.016
#> H3 0.115  61.147 1.000
#> --- 
#> note: equal hypothesis prior probabilities
```

The posterior hypothesis probability is `0.895` which provides some evidence for the order 
constraint. The Bayes factor matrix then divides the posterior probabilities. This provide a measure
of *relative* support for which hypothesis the data were more likely under. 

Finally, the results can be plotted

```{r, eval=FALSE}
plot(fit) + 
  scale_fill_brewer(palette = "Set2", 
                    name = "Posterior Prob") +
  ggtitle("Confirmatory: Comorbidity Network")
```


![](readme_models/confirm_hyp.png)

This demonstrates that all the `plot()` functions in **BGGM** return `ggplot` objects that can be further customized. Note that **BGGM** is not focused on making publication ready plots. Typically the bare minimum is provided that can then be honed in.

<br>

### Comparing Gaussian Graphical Models

#### Partial Correlation Differences

This method compares groups by computing the difference for each relation in the model. In other 
words, there are pairwise contrasts for each partial correlation, resulting in a posterior 
distribution for each difference.

In all examples in this section, personality networks are compared for males and females.

```{r}
# data
Y <- bfi

# males
Ymales <- subset(Y, gender == 1, 
                 select = -c(gender, education))

# females
Yfemales <- subset(Y, gender == 2, 
                 select = -c(gender, education))
```

Fit the model
```{r, eval=FALSE}
fit <- ggm_compare_estimate(Ymales, Yfemales)
```

Then plot the results, in this case the posterior distribution for each difference

```{r, eval=FALSE}
# plot summary
plot(summary(fit))

```

![](readme_models/ggm_compare_estimate.png)

Note that it is also possible to use `select` for the object `fit` and then plot the results. This 
produces a network plot including the selected differences.
Furthermore, it is also possible to plot the partial correlations (not the differences). This
is accomplished by using `plot` with the summary computed from an `estimate` object 
([see above](#bayesian-estimation)).

#### Posterior Predictive Check

The predictive check method uses Jensen-Shannon divergence (i.e., symmetric Kullback-Leibler divergence [Wikipedia](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence)) and the sum of squared error (for the partial correlation matrices) to compare groups [@williams2020comparing].

The following compares the groups
```{r, eval=FALSE}
fit <- ggm_compare_ppc(Ymales, Yfemales)
```

Then print the summary output with
```{r, eval=FALSE}
fit

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Test: Global Predictive Check 
#> Posterior Samples: 500 
#>   Group 1: 805 
#>   Group 2: 1631 
#> Nodes:  25 
#> Relations: 300 
#> --- 
#> Call: 
#> ggm_compare_ppc(Ymales, Yfemales, iter = 500)
#> --- 
#> Symmetric KL divergence (JSD): 
#>  
#>    contrast JSD.obs p_value
#>  Yg1 vs Yg2   0.442       0
#> --- 
#>  
#> Sum of Squared Error: 
#>  
#>    contrast SSE.obs p.value
#>  Yg1 vs Yg2   0.759       0
#> --- 
#> note:
#> JSD is Jensen-Shannon divergence 
```

In this case, there seems to be decisive evidence that the networks are different (as indicated by the
posterior predictive *p*-value). The predictive distribution can also be plotted
```{r, eval=FALSE, out.width = '65%', warning = FALSE, message=FALSE}
plot(fit, 
     critical = 0.05)$plot_jsd
```

![](readme_models/ppc_1.png)

where the red region is the "critical" area and the black point is the observed KL divergence for the networks. This again shows that the "distance" between the networks is much more than expected, assuming that the groups were actually the same.


This next example is a new feature in **BGGM** (`2.0.0`), that allows for comparing GGMs any way the user wants. All that is required is to (1) decide on a test-statistic and (2) write a custom function. 

Here is an example using Hamming distance ([Wikipedia](https://en.wikipedia.org/wiki/Hamming_distance)), which is essentially the squared error between adjacency matrices (a test for different structures).

First define the custom function
```{r}
f <- function(Yg1, Yg2){

# remove NA
x <- na.omit(Yg1)
y <- na.omit(Yg2)

# nodes
p <- ncol(x)

# identity matrix
I_p <- diag(p)

# estimate graphs
fit1 <-  estimate(x, analytic = TRUE)
fit2 <-  estimate(y, analytic = TRUE)

# select graphs
sel1 <- select(fit1)
sel2 <- select(fit2)

# Hamming distance
sum((sel1$adj[upper.tri(I_p)] - sel2$adj[upper.tri(I_p)])^2)
}
```

Note that (1) `analytic = TRUE` is being used, which is needed in this case because two graphs are 
estimated for each iteration (or draw from the posterior predictive distribution) and (2) `f` requires two datasets as the input and returns a single number (the chosen test-statistic). The next step is to compute the observed Hamming distance

```{r, echo=FALSE, warning = FALSE}
load(file = "readme_models/obs.rda")
```

```{r, eval=FALSE}
# observed difference
obs <- f(Ymales, Yfemales)
```

then compare the groups

```{r, eval = FALSE}
fit <- ggm_compare_ppc(Ymales, Yfemales,
                             FUN = f,
                             custom_obs  = obs)

# print
fit

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Test: Global Predictive Check 
#> Posterior Samples: 250 
#>   Group 1: 805 
#>   Group 2: 1631 
#> Nodes:  25 
#> Relations: 300 
#> --- 
#> Call: 
#> ggm_compare_ppc(Ymales, Yfemales, iter = 250, FUN = f, custom_obs = obs)
#> --- 
#> Custom: 
#>  
#>    contrast custom.obs p.value
#>  Yg1 vs Yg2         75   0.576
#> --- 
```

In this case, the *p*-value does not indicate that the groups are different for this test-statistic. This may seem contradictory to the previous results, but it is important to note that
Hamming distance asks a much different question related to the adjacency matrices (no other
information, such as edge weights, is considered).

#### Exploratory (groups)
The Bayes factor based methods allow for determining the conditional 
**in**dependence structure (evidence for the null hypothesis), in this case for group equality.

Fit the model
```{r, eval=FALSE}
fit <- ggm_compare_explore(Ymales, Yfemales)
```

Then plot the results
```{r, eval=FALSE}
plot(summary(fit)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_blank()) 
```

![](readme_models/plt_ggm_compare_explore.png)

Here the posterior probability for a difference is visualized for each relation in the GGM. Note 
that it is also possible to use `select` for the object `fit` and then plot the results. This 
produces a network plot including the selected differences, in addition to a plot depicting the 
relations for which there was evidence for the null hypothesis.

#### Confirmatory (groups)
A central contribution of **BGGM** is confirmatory hypothesis testing of (in)equality constraints [@Hoijtink2011],
in this case for comparing groups. By this we are referring to testing expectations, as opposed to feeding the data to, say, `estimate`, and seeing what happens to emerge. 


In this example, the focus is on agreeableness in a personality network. Here is a set of hypotheses
```{r}
hyp <- c("g1_A2--A4 > g2_A2--A4 > 0 & g1_A4--A5 > g2_A4--A5 > 0;
          g1_A4--A5 = g2_A4--A5 = 0  & g1_A2--A4 = g2_A2--A4 = 0")
```

where the variables are `A2` ("inquire about others' well being"), `A4` ("love children"), 
and `A5` ("make people feel at ease"). The first hypothesis states that the conditionally 
dependent effects are larger for female than males (note the `&`), with the additional 
constraint to positive values, whereas the second hypothesis is a "null" model.

The hypothesis is tested with the following
```{r, eval=FALSE}
fit <- ggm_compare_confirm(Yfemales, Ymales, 
                           hypothesis = hyp)

# print
fit

#> BGGM: Bayesian Gaussian Graphical Models
#> Type: continuous
#> ---
#> Posterior Samples: 500
#>   Group 1: 1631
#>   Group 2: 805
#> Variables (p): 25
#> Relations: 300
#> Delta: 15
#> ---
#> Call:
#> ggm_compare_confirm(Yfemales, Ymales, hypothesis = hyp, iter = 500)
#> ---
#> Hypotheses:
#> 
#> H1: g1_A2--A4>g2_A2--A4>0&g1_A4--A5>g2_A4--A5>0
#> H2: g1_A4--A5=g2_A4--A5=0&g1_A2--A4=g2_A2--A4=0
#> H3: complement
#> ---
#> Posterior prob:
#> 
#> p(H1|data) = 0.989
#> p(H2|data) = 0
#> p(H3|data) = 0.011
#> ---
#> Bayes factor matrix:
#>       H1           H2     H3
#> H1 1.000 1.180798e+14 92.115
#> H2 0.000 1.000000e+00  0.000
#> H3 0.011 1.281873e+12  1.000
#> ---
#> note: equal hypothesis prior probabilities
```

The posterior hypothesis probability is 0.989 which provides strong evidence for the 
hypothesis that predicted these "agreeableness" relations would be larger in females 
than in males. This can also be plotted, as in [Confirmatory (one group)](#confirmatory). 
See @rodriguez2020formalizing for a full treatment of confirmatory testing in substantive applications.

### Beyond the Conditional (In)dependence Structure

#### Predictability

In this example, predictability is computed for each node in the network [see here for rationale @haslbeck2018well]. Currently **BGGM** computes Bayesian variance explained for all data types [@gelman_r2_2019].  

The following computes predictability for binary data

```{r, eval=FALSE}
# binary
Y <- women_math

# fit model
fit <- estimate(Y, type = "binary")

# compute r2
r2 <- predictability(fit, iter = 500)

# plot
plot(r2, type = "ridgeline")
```

![](readme_models/predictability.png)

#### Posterior Uncertainty

See [Partial Correlation Differences](#partial-correlation-differences)

#### Custom Network Statistics
A new feature to **BGGM** allows for computing user defined network statistics, given a partial correlation or
weighted adjacency matrix.

Here is an example for bridge centrality [@jones2019bridge]. The first step is to define the function 

```{r, eval=FALSE}
# need this package 
library(networktools)

# custom function
f <- function(x, ...){
 bridge(x, ...)$`Bridge Strength`
}
```

Note that `x` takes the matrix and `f` can return either a single number or a number for each node. The next step is to fit the model and compute the network statistic

```{r, eval=FALSE}
# data
Y <- ptsd

# clusters
communities <- substring(colnames(Y), 1, 1)

# estimate the model
fit <- estimate(Y)

# bridge strength
net_stat <- roll_your_own(fit,
                          FUN = f,
                          select = TRUE,
                          communities = communities)
```

The function `f` is provided to `FUN` and `communities` is passed to `brigde` (inside of `f`) via `...`. The results can be printed

```{r, eval=FALSE}
# print
net_stat

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Network Stats: Roll Your Own
#> Posterior Samples: 100 
#> --- 
#> Estimates: 
#> 
#>  Node Post.mean Post.sd Cred.lb Cred.ub
#>     1     0.340   0.097   0.166   0.546
#>     2     0.319   0.100   0.176   0.513
#>     3     0.000   0.000   0.000   0.000
#>     4     0.337   0.086   0.189   0.489
#>     5     0.559   0.133   0.332   0.791
#>     6     0.188   0.073   0.029   0.320
#>     7     0.505   0.138   0.241   0.781
#>     8     0.153   0.070   0.022   0.286
#>     9     0.175   0.063   0.041   0.281
#>    10     0.000   0.000   0.000   0.000
#>    11     0.365   0.107   0.178   0.627
#>    12     0.479   0.093   0.280   0.637
#>    13     0.155   0.074   0.022   0.301
#>    14     0.000   0.000   0.000   0.000
#>    15     0.374   0.097   0.175   0.550
#>    16     0.174   0.065   0.034   0.295
#>    17     0.000   0.000   0.000   0.000
#>    18     0.491   0.132   0.238   0.745
#>    19     0.613   0.113   0.408   0.825
#>    20     0.144   0.066   0.038   0.289
#> --- 
```

And then plotted
```{r, eval = FALSE}
plot(net_stat)
```

![](readme_models/bridge.png)

There are additional examples in the documentation.

### Example Network Plot
Here is an example of a more involved network plot. In this case, the graph is estimated with a semi-parametric copula (`type = "mixed"`), where two control variables are included in the model.

```{r, echo=FALSE}
load(file = "readme_models/fit_plt_ex.rda")
fit <- fit_plt_ex
# select graph
E <- select(fit)
```

```{r, eval=FALSE}
# personality (includes gender and education)
Y <- bfi

# fit copula GGM
fit <- estimate(Y, type = "mixed")

# select graph
E <- select(fit)
```

The graph is then plotted
```{r, eval=FALSE }
# extract communities
comm <- substring(colnames(Y), 1, 1)

# plot
plot(E, 
     # enlarge edges
     edge_magnify = 5, 
     # cluster nodes
     groups = comm, 
     # change layout
     layout = "circle")$plt +
  # plot title
  ggtitle("Semi-Parametric Copula") +
  # add custom labels
  scale_color_brewer(breaks = c("A", "C", 
                                "E", "N", 
                                "O", "e",  
                                "g"), 
                     labels =   c("A", "C", 
                                 "E", "N", 
                                 "O",  
                                 "Education",   
                                 "Gender"), 
                     palette = "Set2")
```

![](readme_models/plt_net_example.png)

Note that `layout` can be changed to any option provided in the `R` package **sna** [@sna].


## Additional Features
The primary focus of **BGGM** is Gaussian graphical modeling (the inverse covariance matrix).
The residue is a suite of useful methods not explicitly for GGMs. For example, 

### Bivariate Correlations

Bivariate correlations for `binary` (tetrachoric), `ordinal` (polychoric), `mixed` (rank based),
and `continuous` (Pearson's) data.
  
Here is an example for computing tetrachoric correlations:

```{r, echo=FALSE}
load(file = "readme_models/binary_cors.rda")
```

```{r, eval=FALSE}
# binary data
Y <- women_math[1:500,]

cors <- zero_order_cors(Y, type = "binary", iter = 250)

cors$R
```

```{r, echo = FALSE, results='asis'}
row.names(cors$R_mean) <- 1:6
knitr::kable(round(cors$R_mean, 3),
             col.names = c("1", "2", "3","4","5", "6"), 
             row.names = TRUE)
```

The object `cors` also includes the sampled correlation matrices (in this case 250) in an array. 

### Multivariate Regression

Multivariate regression for binary (probit), ordinal (probit),
mixed (rank likelihood), and continuous data.
  
Here is an example for a multivariate probit model with an ordinal outcome, where 
`E5` ("take charge") and `N5` ("panic easily") are predicted by `gender` and `education`:

```{r,echo=FALSE}
load("readme_models/mv_probit.rda")
```

```{r, eval = F}
# personality data
Y <- bfi

# variables
Y <- subset(Y, select = c("E5", "N5", 
                          "gender", "education"))


mv_probit <- estimate(Y, formula = ~ gender + as.factor(education), 
                      type = "ordinal")

```

Note that **BGGM** does not use the customary `model.matrix` formulation. This is for good reason, as 
each variable in the GGM does not need to be written out. Here we effectively "tricked" **BGGM** to 
fit a multivariate probit model (each variable included in `formula` is removed from `Y`). 

```{r}
regression_summary(mv_probit)
```

This basic idea can also be used to fit regression models with a single outcome.

## Note on Conditional (In)dependence Models for Latent Data

All of the data types (besides continuous) model latent data. That is, unobserved data 
that is assumed to be Gaussian distributed. For example, a  tetrachoric correlation 
(binary data) is a special case of a polychoric correlation (ordinal data). 
Both relations are between "theorized normally distributed continuous *latent* 
variables [Wikepedia](https://en.wikipedia.org/wiki/Polychoric_correlation). 
In both instances, the corresponding partial correlation between observed 
variables is conditioned on the remaining variables in the *latent* space. 
This implies that interpretation is similar to continuous data, but with respect 
to latent variables. We refer interested users to 
[see page 2364, section 2.2, in  @webb2008bayesian].


## High Dimensional Data?

**BGGM** was built specifically for social-behavioral scientists. Of course, the methods
can be used by all researchers. However, there is currently *not* support for high-dimensional data
(i.e., more variables than observations) that are common place in, say, the genetics literature.
These data are rare in the social-behavioral sciences. In the future, support for high-dimensional
data may be added to **BGGM**.

## Bug Reports, Feature Requests, and Contributing
Bug reports and feature requests can be made by opening an issue on [Github](https://github.com/donaldRwilliams/BGGM/issues). To contribute towards
the development of **BGGM**, you can start a branch with a pull request and we can 
discuss the proposed changes there.

## Comparison to Other Software
**BGGM** is the only `R` package to implement all of these algorithms and methods. The `mixed` data approach 
is also implemented in the package **sbgcop** [base `R`, @hoff2007extending]. The `R` package **BDgraph** implements a Gaussian copula graphical model in `c++` [@mohammadi2015bdgraph], but not the binary or ordinal approaches. Furthermore, **BGGM** is the only package for confirmatory testing and comparing graphical models with the methods described in @williams2020comparing.

## References
---
output: github_document
bibliography: inst/REFERENCES.bib
---

```{r, echo = FALSE, message=F}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  dev = "png",
  dpi = 500,
  fig.align = "center",
  knitr::opts_chunk$set(comment = NA)
  )
library(ggplot2)
library(BGGM)
```

# Bayesian Gaussian Graphical Models <img src="man/figures/logo.png" align="right" alt="" width="150" />


<!-- badges: start -->
[![CRAN Version](http://www.r-pkg.org/badges/version/BGGM)](https://cran.r-project.org/package=BGGM)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BGGM)
[![Build Status](https://travis-ci.org/donaldRwilliams/BGGM.svg?branch=master)](https://travis-ci.org/donaldRwilliams/BGGM)
<!-- badges: end -->


The `R` package **BGGM** provides tools for making Bayesian inference in 
Gaussian graphical models [GGM, @williams2020bggm]. The methods are organized around 
two general approaches for Bayesian inference: (1) estimation and (2) hypothesis 
testing. The key distinction is that the former focuses on either the posterior or posterior 
predictive distribution [@Gelman1996a; see section 5 in @rubin1984bayesianly], whereas the 
latter focuses on model comparison with the Bayes factor [@Jeffreys1961; @Kass1995].

## <i class="fas fa-cog"></i> Installation
To install the latest release version (`2.0.0`) from CRAN use
```{r gh-installation, eval = FALSE}	
install.packages("BGGM")	
```

The current developmental version can be installed with	

```{r, eval = FALSE}	
if (!requireNamespace("remotes")) {	
  install.packages("remotes")	
}	
remotes::install_github("donaldRwilliams/BGGM")
```


### <i class="fas fa-skull-crossbones"></i> Dealing with Errors
There are automatic checks for [**BGGM**](https://travis-ci.org/github/donaldRwilliams/BGGM/branches). However, that only checks for Linux and **BGGM** is built on Windows. 
The most common installation errors occur on OSX. An evolving guide to address these
issues is provided in the [Troubleshoot Section](https://donaldrwilliams.github.io/BGGM/articles/installation.html). 

## <i class="fas fa-clipboard-list"></i> Overview
The methods in **BGGM** build upon existing algorithms that are well-known in the literature.
The central contribution of **BGGM** is to extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution [@Mulder2018]
  
    + Estimation [@Williams2019]

2. Bayesian hypothesis testing with the matrix-F prior distribution [@Williams2019_bf]

    + [Exploratory hypothesis testing](https://github.com/donaldRwilliams/BGGM#Exploratory)
  
    + [Confirmatory hypothesis testing](https://github.com/donaldRwilliams/BGGM#confirmatory)
    
3. Comparing Gaussian graphical models [@Williams2019; @williams2020comparing]
    
    + [Partial correlation differences](https://github.com/donaldRwilliams/BGGM#partial-correlation-differences)
    
    + [Posterior predictive check](https://github.com/donaldRwilliams/BGGM#posterior-predictive-check)
    
    + [Exploratory hypothesis testing](https://github.com/donaldRwilliams/BGGM#exploratory-groups)
    
    + [Confirmatory hypothesis testing](https://github.com/donaldRwilliams/BGGM#confirmatory-groups)

4. Extending inference beyond the conditional (in)dependence structure [@Williams2019]

    +  [Predictability](https://github.com/donaldRwilliams/BGGM#predictability)
    
    +  [Posterior uncertainty intervals](https://github.com/donaldRwilliams/BGGM#posterior-uncertainty) for the 
       partial correlations
       
    +  [Custom Network Statistics](https://github.com/donaldRwilliams/BGGM#custom-network-statistics)
    
    
The computationally intensive tasks are written in `c++` via the `R` package **Rcpp** [@eddelbuettel2011rcpp] and the `c++` library **Armadillo** [@sanderson2016armadillo]. The Bayes factors are computed with the `R` package **BFpack** [@mulder2019bfpack]. Furthermore, there are plotting functions
for each method, control variables can be included in the model (e.g., `~ gender`), 
and there is support for missing values (see `bggm_missing`).

## Supported Data Types

* **Continuous**: The continuous method was described in  @Williams2019. Note that 
                  this is based on the customary [Wishart distribution](https://en.wikipedia.org/wiki/Wishart_distribution).

* **Binary**: The binary method builds directly upon @talhouk2012efficient
  that, in turn, built upon the approaches of @lawrence2008bayesian and
  @webb2008bayesian (to name a few).
  
* **Ordinal**: The ordinal methods require sampling thresholds. There are two approach 
   included in **BGGM**. The customary approach described in @albert1993bayesian 
   (the default) and the 'Cowles' algorithm described in @cowles1996accelerating.
   
* **Mixed**: The mixed data (a combination of discrete and continuous) method was introduced
 in @hoff2007extending. This is a semi-parametric copula model
 (i.e., a copula GGM) based on the ranked likelihood. Note that this can be used for 
 *only* ordinal data (not restricted to "mixed" data).


##  <i class="fas fa-folder-open"></i> Illustrative Examples
There are several examples in the [Vignettes](https://donaldrwilliams.github.io/BGGM/articles/) section.

## <i class="fas fa-play-circle"></i> Basic Usage
It is common to have some combination of continuous and discrete (e.g., ordinal, binary, etc.) variables. **BGGM** (as of version `2.0.0`) can readily be used for these kinds of data. In this example, a model is fitted for the `gss` data in **BGGM**. 

### Visualize
The data are first visualized with the **psych** package, which readily shows the data are "mixed". 

```r
# dev version
library(BGGM)
library(psych)

# data
Y <- gss

# histogram for each node
psych::multi.hist(Y, density = FALSE)
```
![](man/figures/index_hist.png)




### Fit Model
A Gaussian copula graphical model is estimated as follows
```r
fit <- estimate(Y, type = "mixed")
```
`type` can be `continuous`, `binary`, `ordinal`, or `mixed`. Note that `type` is a misnomer, as the data can consist of *only* ordinal variables (for example).


### Summarize Relations
The estimated relations are summarized with
```r
summary(fit)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: mixed 
#> Analytic: FALSE 
#> Formula:  
#> Posterior Samples: 5000 
#> Observations (n): 464  
#> Nodes (p): 7 
#> Relations: 21 
#> --- 
#> Call: 
#> estimate(Y = Y, type = "mixed")
#> --- 
#> Estimates:
#>  Relation Post.mean Post.sd Cred.lb Cred.ub
#>  INC--DEG     0.463   0.042   0.377   0.544
#>  INC--CHI     0.148   0.053   0.047   0.251
#>  DEG--CHI    -0.133   0.058  -0.244  -0.018
#>  INC--PIN     0.087   0.054  -0.019   0.196
#>  DEG--PIN    -0.050   0.058  -0.165   0.062
#>  CHI--PIN    -0.045   0.057  -0.155   0.067
#>  INC--PDE     0.061   0.057  -0.050   0.175
#>  DEG--PDE     0.326   0.056   0.221   0.438
#>  CHI--PDE    -0.043   0.062  -0.162   0.078
#>  PIN--PDE     0.345   0.059   0.239   0.468
#>  INC--PCH     0.052   0.052  -0.052   0.150
#>  DEG--PCH    -0.121   0.056  -0.228  -0.012
#>  CHI--PCH     0.113   0.056   0.007   0.224
#>  PIN--PCH    -0.080   0.059  -0.185   0.052
#>  PDE--PCH    -0.200   0.058  -0.305  -0.082
#>  INC--AGE     0.211   0.050   0.107   0.306
#>  DEG--AGE     0.046   0.055  -0.061   0.156
#>  CHI--AGE     0.522   0.039   0.442   0.594
#>  PIN--AGE    -0.020   0.054  -0.122   0.085
#>  PDE--AGE    -0.141   0.057  -0.251  -0.030
#>  PCH--AGE    -0.033   0.051  -0.132   0.063
#> --- 
```

The summary can also be plotted
```r
plot(summary(fit))
```
![](man/figures/index_summ.png)

### Graph Selection
The graph is selected and plotted with
```r
E <- select(fit)

plot(E, node_size = 12,
     edge_magnify = 5)
```

![](man/figures/netplot_index.png)


The Bayes factor testing approach is readily implemented by changing `estimate` to `explore`. 

## <i class="fas fa-pen-square"></i> References
---
title: "In Tandem: Confirmatory and Exploratory Testing"
author: "Donny Williams"
date: "5/23/2020"
bibliography: ../inst/REFERENCES.bib
output: 
  rmarkdown::html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Confirmatory and Exploratory Testing}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```



The blog post, "Tutorial: Bayesian Testing of Central Structures in Psychological Networks," is hosted on a different website.

# <i class="fas fa-external-link-square-alt"></i> [External Link](https://josue.rbind.io/post/tutorial-bayesian-testing/)
---
title: "Graphical VAR"
author: "Donny Williams"
date: "6/04/2020"
bibliography: ../inst/REFERENCES.bib
output:
rmarkdown::html_vignette:
  toc: yes
vignette: >
  %\VignetteIndexEntry{Graphical VAR}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Introduction
This vignette describes the implemention for a "graphical VAR" in `BGGM`. In `BGGM`, this is fitted as a multivariate regression. The  key innovation is a novel prior distribution for the residual covariance matrix. There are a variety of much cooler names than a *mere* "multivariate regression", including "VAR" (vector autoregressive models) and "TSCGM" (time series chain graphical model).

## R package
```
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")
library(BGGM)

# for comparsion
library(vars)

# for plotting
library(qgraph)

# combine plots
library(cowplot)
```

## Data 
I use data from the ifit (fit bit) study. The data were gathered over 100 consecutive days on a variety of variables, including the PANAS scale (positive and negative affect) and the number of steps each day. `BGGM` includes a subset of variables for two individuals.

```r
# data
Y <- subset(ifit, id == 1)[,-1]

# first 3 rows
head(Y, n = 3)

#>  interested disinterested excited upset strong stressed steps
#>         72            10      50     2     50       16  7805
#>         75             6      75     0     76        0 18248
>#         36            58      38     5     45        1 12139
```



# Estimation 
The methods in **BGGM** are organized around Bayesian "estimation" and "hypothesis testing". 
This is to reach a broader audience, as former is more similar to classical 
methods (those more familiar to researchers).

## Fit Model
With the data in hand, the model is fitted as follows

```
# fit model
fit <- var_estimate(Y, beta_sd = 1)
```
Note that `beta_sd` is the prior distribution for the regression coefficients. A smaller value, say, `beta_sd = 0.25`, results in a Bayesian ridge regression. Note also this model, including 5000 draws from the posterior, was estimated in less than 1 second.

The results can then be printed
```r
# print 
fit

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Vector Autoregressive Model (VAR) 
#> --- 
#> Posterior Samples: 5000 
#> Observations (n): 94 
#> Nodes (p): 7 
#> --- 
#> Call: 
#> var_estimate(Y = Y, beta_sd = 10)
#> --- 
#> Partial Correlations: 
#> 
#>               interested disinterested excited  upset strong stressed  steps
#> interested         0.000        -0.170   0.388 -0.217  0.313    0.268  0.089
#> disinterested     -0.170         0.000  -0.172 -0.029  0.094    0.160 -0.078
#> excited            0.388        -0.172   0.000 -0.126  0.500   -0.161 -0.016
#> upset             -0.217        -0.029  -0.126  0.000  0.118    0.350 -0.039
#> strong             0.313         0.094   0.500  0.118  0.000   -0.010  0.176
#> stressed           0.268         0.160  -0.161  0.350 -0.010    0.000 -0.038
#> steps              0.089        -0.078  -0.016 -0.039  0.176   -0.038  0.000
#> --- 
#> Coefficients: 
#> 
#>                  interested disinterested excited  upset strong stressed  steps
#> interested.l1         0.230        -0.009   0.182 -0.102  0.178    0.018  0.113
#> disinterested.l1     -0.051        -0.007   0.056 -0.019  0.049    0.091 -0.023
#> excited.l1           -0.088        -0.196   0.003  0.057 -0.093    0.092  0.106
#> upset.l1             -0.155         0.262  -0.097  0.435  0.057    0.324 -0.091
#> strong.l1             0.026         0.182   0.026  0.048  0.189   -0.073 -0.196
#> stressed.l1          -0.021        -0.014  -0.033 -0.048 -0.079    0.152  0.133
#> steps.l1             -0.157         0.180  -0.211  0.155 -0.092    0.209  0.042
#> --- 
#> Date: Thu Jun 04 08:54:04 2020 
```

Note that the coefficients are comparable, given each variable has been standardized (e.g., the predictors
and the outcome are standardized). `BGGM` does not compute the partial directed correlation (PDC) by default (as in **graphicalVAR**). This is because the standardized effects can readily be tested with the Bayes factor, both across and within each model, whereas this does not seem straightforward for the PDC (which requires a transformation).

### Compare to Classical
Here are the estimates from the `vars` package

```r
t(round(
  vars::Bcoef( 
  vars:::VAR(scale(na.omit(Y)), type = "none")), 
  digits = 3)
)

#>                  interested disinterested excited  upset strong stressed  steps
#> interested.l1         0.229        -0.012   0.184 -0.100  0.180    0.015  0.112
#> disinterested.l1     -0.050        -0.006   0.057 -0.019  0.050    0.092 -0.022
#> excited.l1           -0.088        -0.193   0.002  0.056 -0.091    0.093  0.106
#> upset.l1             -0.155         0.260  -0.096  0.436  0.058    0.321 -0.092
#> strong.l1             0.027         0.182   0.025  0.047  0.188   -0.073 -0.192
#> stressed.l1          -0.021        -0.012  -0.033 -0.046 -0.077    0.152  0.133
#> steps.l1             -0.157         0.183  -0.210  0.153 -0.093    0.207  0.041
```

Recall that the "estimation" methods are similar to, in this case, ordinary least squares. The graphical structure in `BGGM` is determined with credible intervals, which will be quite similar to using confidence
intervals. Hence for those researchers unfamiliar with Bayesian methods the "estimation" methods are perhaps
a nice place to start.

## Summarize Model
The model can also be summarized with

```r
print(
  summary(fit,  cred = 0.95), 
  param = "pcor"
  )


#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Vector Autoregressive Model (VAR) 
#> --- 
#> Partial Correlations: 
#> 
#>                   Relation Post.mean Post.sd Cred.lb Cred.ub
#>  interested--disinterested    -0.170   0.108  -0.382   0.044
#>        interested--excited     0.388   0.085   0.219   0.546
#>     disinterested--excited    -0.172   0.104  -0.369   0.049
#>          interested--upset    -0.217   0.106  -0.417   0.000
#>       disinterested--upset    -0.029   0.101  -0.239   0.161
#>             excited--upset    -0.126   0.098  -0.315   0.066
#>         interested--strong     0.313   0.090   0.135   0.480
#>      disinterested--strong     0.094   0.112  -0.120   0.318
#>            excited--strong     0.500   0.078   0.337   0.645
#>              upset--strong     0.118   0.109  -0.100   0.325
#>       interested--stressed     0.268   0.102   0.058   0.460
#>    disinterested--stressed     0.160   0.100  -0.049   0.351
#>          excited--stressed    -0.161   0.099  -0.358   0.031
#>            upset--stressed     0.350   0.091   0.166   0.519
#>           strong--stressed    -0.010   0.107  -0.212   0.201
#>          interested--steps     0.089   0.108  -0.123   0.297
#>       disinterested--steps    -0.078   0.108  -0.284   0.125
#>             excited--steps    -0.016   0.100  -0.207   0.182
#>               upset--steps    -0.039   0.107  -0.245   0.178
#>              strong--steps     0.176   0.101  -0.024   0.364
#>            stressed--steps    -0.038   0.108  -0.236   0.193
#> --- 
```

The coefficients can also be printed by changing `param` to either `all` or `beta`, The summary can also be plotted. Here are the coefficients

```r
plts <- plot(summary(fit,  cred = 0.95))

cowplot::plot_grid(
  cowplot::plot_grid(
                   plts$beta_plt$interested, 
                   plts$beta_plt$disinterested, 
                   plts$beta_plt$excited, 
                   nrow = 1),
cowplot::plot_grid(  
                   plts$beta_plt$upset,
                   plts$beta_plt$strong,
                   plts$beta_plt$stressed, 
                   nrow = 1
                  ), 
                  nrow = 2)
```
![](../man/figures/var_coef.png)


There is a plot for the partial correlations in the object `plts`.

## Select Graph
The graphs are selected with 

```r
select(fit, cred = 0.95)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Vector Autoregressive Model (VAR) 
#> --- 
#> Posterior Samples: 5000 
#> Credible Interval: 95 % 
#> --- 
#> Call: 
#> var_estimate(Y = Y, beta_sd = 10)
#> --- 
#> Partial Correlations: 
#> 
#>               interested disinterested excited  upset strong stressed steps
#> interested         0.000             0   0.388 -0.217  0.313    0.268     0
#> disinterested      0.000             0   0.000  0.000  0.000    0.000     0
#> excited            0.388             0   0.000  0.000  0.500    0.000     0
#> upset             -0.217             0   0.000  0.000  0.000    0.350     0
#> strong             0.313             0   0.500  0.000  0.000    0.000     0
#> stressed           0.268             0   0.000  0.350  0.000    0.000     0
#> steps              0.000             0   0.000  0.000  0.000    0.000     0
#> --- 
#> Coefficients: 
#> 
#>                  interested disinterested excited upset strong stressed steps
#> interested.l1             0         0.000       0 0.000      0    0.000     0
#> disinterested.l1          0         0.000       0 0.000      0    0.000     0
#> excited.l1                0         0.000       0 0.000      0    0.000     0
#> upset.l1                  0         0.262       0 0.435      0    0.324     0
#> strong.l1                 0         0.000       0 0.000      0    0.000     0
#> stressed.l1               0         0.000       0 0.000      0    0.000     0
#> steps.l1                  0         0.000       0 0.000      0    0.209     0
#> --- 

```


# Plot Graph
For plotting, I use the **qgraph** package.

```r
par(mfrow=c(1,2))
qgraph::qgraph(sel$pcor_weighted_adj, title = "Partials")
qgraph::qgraph(sel$beta_weighted_adj, title = "Coefficients")
```
![](../man/figures/est_net.png)

# Predictability
Finally, it is also possible to compute predictability, in this case Bayesian $R^2$

```r
r2 <- predictability(fit)

# print
r2

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Metric: Bayes R2
#> Type: continuous 
#> --- 
#> Estimates:
#> 
#>           Node Post.mean Post.sd Cred.lb Cred.ub
#>     interested     0.144   0.057   0.050   0.271
#>  disinterested     0.166   0.061   0.060   0.302
#>        excited     0.127   0.054   0.039   0.250
#>          upset     0.220   0.070   0.093   0.368
#>         strong     0.116   0.051   0.035   0.232
#>       stressed     0.227   0.069   0.102   0.373
#>          steps     0.105   0.047   0.032   0.210


```

The object `r2` can also be plotted

```r
plot(r2, type = "ridgeline")
```

![](../man/figures/var_ridgeline.png)


# Explore
Bayesian (exploratory) testing to come...

# Confirm

Bayesian (confirmatory) testing to come...

# Note
---
title: "Three Ways to Test the Same Hypothesis"
author: "Donny Williams"
date: "5/23/2020"
bibliography: ../inst/REFERENCES.bib
output:
rmarkdown::html_vignette:
  toc: yes
vignette: >
  %\VignetteIndexEntry{Three Ways to Test the Same Hypothesis}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Introduction
On a Facebook methods group, there was a question about testing hypotheses in networks. In 
the comments, it was suggested that **BGGM** could be used to test the hypothesis. And it turns
out that **BGGM** really shines for testing expectations [see for example @rodriguez2020formalizing]. 

In this vignette, I demonstrate three ways to go about testing the same hypothesis, which is
essentially testing for a difference in the **sum** of partial correlations between groups.

### R package
```{r, eval = FALSE, message=FALSE}
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")
library(BGGM)
```

### Data
For demonstrative purposes, I use the `bfi` data and test the hypotheses in males and females

```
# data
Y <- bfi

# males
Y_males <- subset(Y, gender == 1, select = -c(education, gender))[,1:5]

# females
Y_females <- subset(Y, gender == 2, select = -c(education, gender))[,1:5]
```

# Approach 1: Posterior Difference
The first approach is rather straightforward, with the caveat that the method needs to be 
implemented by the user. Note that I could certainly implement this in **BGGM**, assuming there 
is enough interest. Please make a feature request [here](https://github.com/donaldRwilliams/BGGM/issues).

## Hypothesis
The hypothesis was that a sum of relations was larger in one group, for example,

$$
\begin{align}
\mathcal{H}_0: (\rho^{male}_{A1--A2}\; + \; \rho^{male}_{A1--A3}) = (\rho^{female}_{A1--A2}\; + \; \rho^{female}_{A1--A3}) \\
\mathcal{H}_1: (\rho^{male}_{A1--A2}\; + \; \rho^{male}_{A1--A3}) > (\rho^{female}_{A1--A2}\; + \; \rho^{female}_{A1--A3})
\end{align}
$$
Note that the hypothesis is related to the sum of relations, which is readily tested in **BGGM**. 

## Fit Models
The first step is to estimate the model for each group

```r
# fit female
fit_female <- estimate(Y_females, seed = 2)

# fit males
fit_male <- estimate(Y_males, seed = 1)
```
For an example, I used the default which is to assume the data is Gaussian. This can be changed with `type = ` either `binary`, `ordinal`, or `mixed`.

## Extract the Samples
The next step is to extract the posterior samples for each relation

```r
post_male <- posterior_samples(fit_male)[,c("A1--A2", "A1--A3")]

post_female <- posterior_samples(fit_female)[,c("A1--A2", "A1--A3")]
```

Note that the column names reflect the upper-triangular elements of the 
partial correlation matrix. Hence, the first name (e.g.,`A1`) must be located before 
the second name (e.g., `A2`) in the data matrix. This can be understood in reference 
to the column numbers: `1--2` is correct whereas `2--1` will result in an error.

## Sum and Compute Difference
The next step is to sum the relations and compute the difference

```r
# sum males
sum_male <- rowSums(post_male) 

# sum females
sum_female <- rowSums(post_female)

# difference
diff <- sum_male - sum_female
```
which can then be plotted
```r
# three column
par(mfrow=c(1,3))

# male sum
hist(sum_male)

# female sum
hist(sum_female)

# difference
hist(diff)
```
![](../man/figures/hyp_3ways_hist.png)




## Posterior Probability
Next compute the posterior probability the sum is larger in males than females

```r
# posterior prob
mean(sum_male > sum_female)

#> 0.737
```
and then the credible interval for the difference

```
quantile(diff, probs = c(0.025, 0.975))

#>        2.5%       97.5% 
#> -0.06498586  0.12481253 
```

# Approach 2: Predictive Check
The next approach is based on a posterior predictive check. The hypothesis is essentially the same as above, but for the predictive distribution, that is,

$$
\begin{align}
\mathcal{H}_0: (\rho^{male^{yrep}}_{A1--A2}\; + \; \rho^{male^{yrep}}_{A1--A3}) = (\rho^{female^{yrep}}_{A1--A2}\; + \; \rho^{female^{yrep}}_{A1--A3}) \\
\mathcal{H}_1: (\rho^{male^{yrep}}_{A1--A2}\; + \; \rho^{male^{yrep}}_{A1--A3}) > (\rho^{female^{yrep}}_{A1--A2}\; + \; \rho^{female^{yrep}}_{A1--A3})
\end{align}
$$
where the only difference is $yrep$. See more details [here](https://donaldrwilliams.github.io/BGGM/articles/ppc_custom.html).

## Define Function
The first step is to define a function to compute the difference in sums
```r
# colnames
cn <- colnames(Y_males)

# function
f <- function(Yg1, Yg2){
  
  # data
  Yg1 <- na.omit(Yg1)  
  Yg2 <- na.omit(Yg2)

  # estimate partials
  fit1 <- pcor_mat(estimate(Yg1, analytic = TRUE))
  fit2 <- pcor_mat(estimate(Yg2, analytic = TRUE))
  
  # names (not needed)
  colnames(fit1) <- cn
  rownames(fit1) <- cn
  colnames(fit2) <- cn
  rownames(fit2) <- cn
  
  # take sum
  sum1 <- fit1["A1", "A2"] + fit1["A1", "A3"]
  sum2 <- fit2["A1", "A2"] + fit2["A1", "A3"]
  
  # difference
  sum1 - sum2

}
```

Note that the function takes two data matrices and then returns a single value. 
Also, the default in **BGGM** does not require a custom function 
(only needs the data from each group).

## Predictive Check
The next step is to compute the observed difference and then perform the check.

```r
# observed
obs <- f(Y_males, Y_females)

# check
ppc <- ggm_compare_ppc(Y_males, Y_females, 
                       iter = 250, 
                       FUN = f, 
                       custom_obs = obs)

# print
ppc

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Test: Global Predictive Check 
#> Posterior Samples: 250 
#>   Group 1: 896 
#>   Group 2: 1813 
#> Nodes:  5 
#> Relations: 10 
#> --- 
#> Call: 
#> ggm_compare_ppc(Y_males, Y_females, iter = 250, FUN = f, custom_obs = obs)
#> --- 
#> Custom: 
#>  
#>    contrast custom.obs p.value
#>  Yg1 vs Yg2      0.029   0.264
#> --- 
```

Note this requires the user to determine $\alpha$. 

## Plot 
The check can also be plotted

```r
plot(ppc)
```
![](../man/figures/hyp_3ways_ppc.png)

where the red is the critical region.

# Approach 3: Bayesian Hypothesis Testing
The above approaches cannot provide evidence that the sum is equal. In other words, just because there was 
not a difference, this does not provide evidence for equality. The Bayes factor methods allow for formally 
assessing the equality model, that is,

$$
\begin{align}
\mathcal{H}_1&: (\rho^{male}_{A1--A2}\; + \; \rho^{male}_{A1--A3}) > (\rho^{female}_{A1--A2}\; + \; \rho^{female}_{A1--A3}) \\
\mathcal{H}_2&: (\rho^{male}_{A1--A2}\; + \; \rho^{male}_{A1--A3}) = (\rho^{female}_{A1--A2}\; + \; \rho^{female}_{A1--A3}) \\
\mathcal{H}_3&: \text{not} \; \mathcal{H}_1 \; \text{or} \; \mathcal{H}_2
\end{align}
$$

where $\mathcal{H}_3$ is the complement and can be understood as neither the first or second hypothesis.

## Test Hypothesis 
The hypothesis is easily translated to `R` code
```r
hyp <- c("g1_A1--A2 + g1_A1--A3 > g2_A1--A2 + g2_A1--A3; 
          g1_A1--A2 + g1_A1--A3 = g2_A1--A2 + g2_A1--A3")
```

Note the `g1` indicates the group and `;` separates the hypotheses. I again assume the data is Gaussian 
(although this can be changed to `type = "ordinal"` or `type = "mixed"`; see [here](https://donaldrwilliams.github.io/BGGM/reference/ggm_compare_confirm.html))

```r
test <- ggm_compare_confirm(Y_males, Y_females, 
                            hypothesis = hyp)


# print
test


#> BGGM: Bayesian Gaussian Graphical Models 
#> Type: continuous 
#> --- 
#> Posterior Samples: 25000 
#>   Group 1: 896 
#>   Group 2: 1813 
#> Variables (p): 5 
#> Relations: 10 
#> Delta: 15 
#> --- 
#> Call:
#> ggm_compare_confirm(Y_males, Y_females, hypothesis = hyp)
#> --- 
#> Hypotheses: 
#> 
#> H1: g1_A1--A2+g1_A1--A3>g2_A1--A2+g2_A1--A3
#> H2: g1_A1--A2+g1_A1--A3=g2_A1--A2+g2_A1--A3
#> H3: complement
#> --- 
#> Posterior prob: 
#> 
#> p(H1|data) = 0.13
#> p(H2|data) = 0.825
#> p(H3|data) = 0.046
#> --- 
#> Bayes factor matrix: 
#>       H1    H2     H3
#> H1 1.000 0.158  2.853
#> H2 6.349 1.000 18.113
#> H3 0.351 0.055  1.000
#> --- 
#> note: equal hypothesis prior probabilities
```


Note the posterior hypothesis probability for the equality model is 0.825. The Bayes factor matrix then divides those values, for example, $BF_{21}$ indicates the data were about 6 times more likely under $\mathcal{H}_2$ than $\mathcal{H}_1$. 

## Plot Hypothesis
The hypothesis can be plotted

```r
plot(test)
```
![](../man/figures/hyp_3ways_pie.png)

### Sensitivity Analysis
It is also important to check the robustness. Here the width of the prior distribution is decreased

```r
test <- ggm_compare_confirm(Y_males, Y_females, 
                            hypothesis = hyp, 
                            prior_sd = 0.15)
# print
test$out_hyp_prob

#> 0.18523406 0.74906147 0.06570447
```
which results in a probability of 0.75 for $\mathcal{H}_2$ ($BF_{21} = 4.04$).


# Conclusion
Three approaches for testing the same hypothesis were demonstrated in this vignette. This highlights that any hypothesis can be tested in **BGGM** and in several ways. 


# References
---
title: "Network Plots"
author: "Donny Williams"
date: "5/20/2020"
bibliography: ../inst/REFERENCES.bib
output: 
  rmarkdown::html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Network Plots}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Introduction
This vignette shows how to make network plots.

### R packages
```{r, eval = FALSE, message=FALSE}
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")
library(BGGM)
library(cowplot)
```
```{r, echo=FALSE, message=FALSE}
library(BGGM)
```

# Estimate
For the estimate methods, it is currently only possible detect non-zero relations and 
the others are set to zero (no connection in the graph). In a future release, it will be possible 
to define a region of equivalence to directly assess null values. Hence, it is important to note those nodes
not connected are not necessarily conditionally independent (absence of evidence is not evidence of absence).

## Fit Model
In this example, I use the `bfi` data which consists of 25 variables measureing different aspects of personality.

```{r, eval=FALSE}
# data
Y <- bfi[,1:25]

# fit model
fit <- estimate(Y)
```

## Select Graph
The next step is to selec the graph or those relations for which the credible excludes zero
```{r, eval=FALSE}
# select the edge set
E <- select(fit, 
            cred = 0.95, 
            alternative = "two.sided")
```

`alternative` can be changed to, say, `"greater"` which would then perform a one-sided hypothesis
test for postive relations. This is ideal for many applications in psychology, because often 
**all** relations are expected to be positive.


## Plot Graph
Here is the basic plot. This works for any object from `select` (e.g., comparing groups).
```{r, eval=FALSE}
plot(E)
```

![](../man/figures/netplot_1.png)

### Customize Plot
The above is `ggplot` that can be futher honed in. Here is an example.
```r
# extract communities
comm <- substring(colnames(Y), 1, 1)

plot(E, 
     # enlarge edges
     edge_magnify = 5, 
     # cluster nodes
     groups = comm, 
     # change layout
     layout = "circle")$plt +
  # add custom labels
  scale_color_brewer(breaks = c("A", 
                                "C", 
                                "E", 
                                "N", 
                                "O"), 
                     labels =   c("Agreeableness", "Conscientiousness", 
                                  "Extraversion", "Neuroticism", 
                                  "Opennness"),
                     palette = "Set2")

```
![](../man/figures/netplot_2.png)

The `edge_magnify` is a value that is multiplied by the edges, `groups` allows for grouping the 
variables (e.g., those thought to belong to the same "community" will be the same color), and the
`scale_color_brewer` is from the package `ggplot2` (`pallete` controls the color of the `groups`). 
By default the edge colors are from a color blind palette. This can be changed in `plot` with 
the arguments `pos_col` (the color for positive edges) and `pos_neg` (the color for negative edges).

This is just scratching the surface of possibilities, as essentially any change 
can be made to the plot. There is lots of support for making nice plots readily available
online.

#### Layout
It is also possible to change the layout. This is done with the **sna** package, which is linked in the documentation for `plot.select` in **BGGM**. Here is an example using `layout = "random"`
```{r, eval=FALSE}
plot(E, 
     # enlarge edges
     edge_magnify = 5, 
     # cluster nodes
     groups = comm, 
     # change layout
     layout = "random")$plt +
  # add custom labels
  scale_color_brewer(breaks = c("A", 
                                "C", 
                                "E", 
                                "N", 
                                "O"), 
                     labels =   c("Agreeableness", "Conscientiousness", 
                                  "Extraversion", "Neuroticism", 
                                  "Opennness"),
                     palette = "Set2")
```
![](../man/figures/netplot_3.png)

# Bayesian Hypothesis Testing
The Bayesian hypothesis testing methods offer several advantages, for example, that 
evidence for the null hypothesis of conditional independence is formally evaluated. 
As a result, the `explore` method in **BGGM** provides plots for both the conditional 
dependence and independence structure, in addition to a plot for which the evidence was 
ambiguous.

To highlight this advantage, `ptsd` data is used that has a relatively small sample size.

```r
# fit model
fit <- explore(Y)

E <- select(fit, BF_cut = 3)
```

Then plot the results. Note that there are three plots, so the package **cowplot** is used
to combine them into one plot.

```r
plts <- plot(E, 
             edge_magnify = 5, 
             groups = comm)

plot_grid(
plts$H1_plt + 
  ggtitle("Conditional Dependence") + 
  theme(legend.position =  "none"),

plts$H0_plt + 
  ggtitle("Conditional Independence") + 
  theme(legend.position =  "none"),
plts$ambiguous_plt + 
  ggtitle("Ambiguous"), 
nrow = 1,
rel_widths = c(1, 1, 1.1)
)
```

![](../man/figures/hyp_plot.png)

As can be seen, there is not evidence for conditional independence for any of the relations. And
the ambiguous network makes clear there is large uncertainty as to what or what might not be the "true" network structure. This basic idea of having three adjacency matrices was proposed in @Williams2019_bf.


# Note
**BGGM** provides a publication ready plot, but it is also limited compared to **qgraph** 
[@epskamp2012qgraph]. The one advantage of **BGGM**  is that all plots are `ggplots` 
which then allows for combining them rather easily. An example is included in another 
vignette that shows how to combine several plots made with various methods in **BGGM**


# References
---
title: "Troubleshoot"
author: "Donny Williams"
date: "5/20/2020"
bibliography: ../inst/REFERENCES.bib
output: 
  rmarkdown::html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Installation}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# OSX

## Error 1: Missing a Fortran Compiler
The most common error seems to be (or similar).

```
E> ld: warning: directory not found for option '-L/usr/local/gfortran/lib'
E> ld: library not found for -lgfortran
E> clang: error: linker command failed with exit code 1 (use -v to see invocation)
```

This indicates that the fortran compiler is missing. This can can be installed [here](https://github.com/fxcoudert/gfortran-for-macOS/releases/tag/8.2?fbclid=IwAR2SyHWB2BzFcY7bpHYW8VzNvlDsy4Gw1QxUMueXB3H0fXicCWoMbE7Ypac) (the `.dmg` file).


## Error 2: Outdated R and/or R-studio
Typically the above has solved the issue. If not, then an additional error could be

```
Error: (converted from warning) Setting LC_CTYPE failed, using "C"
```

This was solved by updating both R and R-studio. More information can be found [here](https://stackoverflow.com/questions/9689104/installing-r-on-mac-warning-messages-setting-lc-ctype-failed-using-c?fbclid=IwAR0DSaPeWOvRyfIsCx4Tjvz9-jZUh2ySXQIHnzqwbqL2_idfPlFF3j6mOe8).

## Error 3: Xcode missing
If that does not work, then perhaps `Xcode` is missing. This can be installed at the "Mac App Store".

## GitHub Issues
The following are links to issues on github for troubleshooting installation of **BGGM** on OSX. 

* [https://github.com/donaldRwilliams/BGGM/issues/26](https://github.com/donaldRwilliams/BGGM/issues/26)(closed)
---
title: "Predictability: Binary, Ordinal, and Continuous"
author: "Donny Williams"
date: "5/20/2020"
bibliography: ../inst/REFERENCES.bib
output: 
  rmarkdown::html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Predictability: Binary, Ordinal, and Continuous}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Background
This vignette describes a new feature to **BGGM** (`2.0.0`) that allows for 
computing network predictability for binary and ordinal data. Currently 
the available option is Bayesian $R^2$ [@gelman_r2_2019].


### R packages
```{r, eval = FALSE, message=FALSE}
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")
library(BGGM)
```

# Binary
The first example looks at Binary data, consisting of 1190 observations and 6 variables. The data are called `women_math` and the variable descriptions are provided in **BGGM**.

The model is estimated with
```{r, eval=FALSE}
# binary data
Y <- women_math

# fit model
fit <- estimate(Y, type = "binary")
```

and then predictability is computed
```{r, eval=FALSE}
r2 <- predictability(fit)

# print
r2

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Metric: Bayes R2
#> Type: binary 
#> --- 
#> Estimates:
#> 
#>  Node Post.mean Post.sd Cred.lb Cred.ub
#>     1     0.016   0.012   0.002   0.046
#>     2     0.103   0.023   0.064   0.150
#>     3     0.155   0.030   0.092   0.210
#>     4     0.160   0.021   0.118   0.201
#>     5     0.162   0.022   0.118   0.202
#>     6     0.157   0.028   0.097   0.208
#> ---
```



There are then two options for plotting. The first is with error bars, denoting the credible interval (i.e., `cred`),
```{r, message=FALSE, eval=FALSE}
plot(r2,
     type = "error_bar",
     size = 4,
     cred = 0.90)
```

![](../man/figures/binary_r2_error.png)

and the second is with a ridgeline plot
```{r, message=FALSE, eval=FALSE}
plot(r2,
     type = "ridgeline",
     cred = 0.50)
```

![](../man/figures/binary_r2_ridge.png)

# Ordinal
In the following, the `ptsd` data is used (5-level Likert). The variable descriptions are provided in **BGGM**. This is based on the polychoric partial correlations, with $R^2$ computed from the corresponding correlations (due to the correspondence between the correlation matrix and multiple regression).

```{r, eval=FALSE}
Y <- ptsd

fit <- estimate(Y + 1, type = "ordinal")
```

The only change is switching type from `"binary` to `ordinal`. One important 
point is the `+ 1`. This is required because for the ordinal approach the first 
category must be 1 (in `ptsd` the first category is coded as 0).

```{r, eval=FALSE}
r2 <- predictability(fit)

# print 
r2 

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Metric: Bayes R2
#> Type: ordinal 
#> --- 
#> Estimates:
#> 
#>  Node Post.mean Post.sd Cred.lb Cred.ub
#>     1     0.487   0.049   0.394   0.585
#>     2     0.497   0.047   0.412   0.592
#>     3     0.509   0.047   0.423   0.605
#>     4     0.524   0.049   0.441   0.633
#>     5     0.495   0.047   0.409   0.583
#>     6     0.297   0.043   0.217   0.379
#>     7     0.395   0.045   0.314   0.491
#>     8     0.250   0.042   0.173   0.336
#>     9     0.440   0.048   0.358   0.545
#>    10     0.417   0.044   0.337   0.508
#>    11     0.549   0.048   0.463   0.648
#>    12     0.508   0.048   0.423   0.607
#>    13     0.504   0.047   0.421   0.600
#>    14     0.485   0.043   0.411   0.568
#>    15     0.442   0.045   0.355   0.528
#>    16     0.332   0.039   0.257   0.414
#>    17     0.331   0.045   0.259   0.436
#>    18     0.423   0.044   0.345   0.510
#>    19     0.438   0.044   0.354   0.525
#>    20     0.362   0.043   0.285   0.454
#> ---
```

Here is the `error_bar` plot.
```{r, eval=FALSE}
plot(r2)
```

![](../man/figures/ordinal_r2_error.png)

Note that the plot object is a `ggplot` which allows for further customization (e.g,. adding the variable names, a title, etc.).

# Continuous
It is quite common to compute predictability assuming that the data are Gaussian. In the context of Bayesian GGMs, this was introduced in [@Williams2019]. This can also be implemented in **BGGM**.

```{r, eval=FALSE}
# fit model
fit <- estimate(Y)

# predictability
r2 <- predictability(fit)
```

`type` is missing which indicates that `continuous` is the default.

# Note
$R^2$ for binary and ordinal data is computed for the underlying latent variables. This is also the case
when `type = "mixed` (a semi-parametric copula). In future releases, there will be support for predicting 
the variables on the observed scale.

# References
---
title: "MCMC Diagnostics"
author: "Donny Williams"
date: "5/20/2020"
bibliography: ../inst/REFERENCES.bib
output: 
  rmarkdown::html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{MCMC Diagnostics}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


# Introduction
The algorithms in **BGGM** are based on Gibbs samplers. In the context of 
covariance matrix estimation, as opposed, to, say, hierarchical models, 
this allows for efficiently sampling the posterior distribution. Furthermore, in all samplers 
the empirical covariance matrix is used as the starting value which reduces 
the length of the burn-in (or warm-up). Still yet it is important to monitor convergence. See [here](http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html) for an introduction to MCMC diagnostics.

### R packages
```{r, eval = FALSE, message=FALSE}
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")
library(BGGM)
```
```{r, echo=FALSE, message=FALSE}
library(BGGM)
```

# ACF plot
This first example includes an "acf" plot that looks at the auto correlation of the samples. In general,
we do not want the samples to be strongly correlated or related to the previous samples (or lags $k$). 
I am not sure there are general guidelines, but typically we do not want "auto correlation...for higher values of k, [because] this indicates a high degree of correlation between our samples and slow mixing " [source](http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html)

Here is an example for ordinal data. 

```{r, eval=FALSE}
# data
Y <- ptsd[,1:10]

# fit model
# + 1 makes first category a 1
fit <- estimate(Y + 1, type = "ordinal")
```

To check the convergence of a partial correlation, we need the parameter name. These are printed as follows
```{r, eval=FALSE}
convergence(fit, print_names = TRUE)

#>  [1] "B1--B2"         "B1--B3"         "B2--B3"         "B1--B4"         "B2--B4"         "B3--B4"         "B1--B5"        
#>  [8] "B2--B5"         "B3--B5"         "B4--B5"         "B1--C1"         "B2--C1"         "B3--C1"         "B4--C1"        
#> [15] "B5--C1"         "B1--C2"         "B2--C2"         "B3--C2"         "B4--C2"         "B5--C2"         "C1--C2"        
#> [22] "B1--D1"         "B2--D1"         "B3--D1"         "B4--D1"         "B5--D1"         "C1--D1"         "C2--D1"        
#> [29] "B1--D2"         "B2--D2"         "B3--D2"         "B4--D2"         "B5--D2"         "C1--D2"         "C2--D2"        
#> [36] "D1--D2"         "B1--D3"         "B2--D3"         "B3--D3"         "B4--D3"         "B5--D3"         "C1--D3"        
#> [43] "C2--D3"         "D1--D3"         "D2--D3"         "B1_(Intercept)" "B2_(Intercept)" "B3_(Intercept)" "B4_(Intercept)"
#> [50] "B5_(Intercept)" "C1_(Intercept)" "C2_(Intercept)" "D1_(Intercept)" "D2_(Intercept)" "D3_(Intercept)"
```


Note the `(Intercept)` which reflect the fact that the ordinal approach is a multivariate probit model with only intercepts.

The next step is to make the plot
```{r, eval=FALSE}
convergence(fit, param = "B1--B2", type = "acf")
```
![](../man/figures/acf_new2.png)



The argument `param` can take any number of parameters and a plot will be made for each (e.g.., `param = c("B1--B2", B1--B3)`). In this case, the auto correlations looks acceptable and actually really good (note the drop to zero). A problematic `acf` plot would have the black lines start at `1.0`
and perhaps never go below `0.20`.  

To make this clear, I simulated time series data taking the code from [here](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/arima.sim.html)

```{r, eval=FALSE}
# sim time series
ts.sim <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 200)

acf(ts.sim)
```

![](../man/figures/acf_problem.png)

This would be considered problematic. If this occurs, one solution could be to thin the samples manually
```{r, eval=FALSE}
# extract samples
samps <- fit$post_samp$pcors

# iterations
iter <- fit$iter

# thinning interval 
thin <-  5

# save every 5th (add 50 which is the burnin)
new_iter <- length(seq(1,to = iter + 50 , by = thin))

# replace (add 50 which is the burnin)
fit$post_samp$pcors <- samps[,,seq(1,to = iter + 50, by = thin)]

# replace iter
fit$iter <- new_iter - 50

# check thinned
convergence(fit, param = "B1--B2", type = "acf")
```

or perhaps just running the model for more iterations (e.g., increasing `iter` in `estimate`). The above is quite convoluted but note convergence should not typically be an issue. And it might come in handy to know that the samples can be replaced and the other functions
in **BGGM** will still work with the object `fit`.


# Trace plot
The next example is a trace plot. Here we are looking for good "mixing".

```{r, eval=FALSE}
convergence(fit, param = "B1--B2", type = "trace")
```

![](../man/figures/trace.png)

Admittedly the term "mixing" is vague. But in general the plot should look like this example,
where there is no place that the chain is "stuck". See [here](https://stats.stackexchange.com/questions/311151/evaluation-of-mcmc-samples) for 
problematic trace plots.
---
title: "Controlling for Variables"
author: "Donny Williams"
date: "5/25/2020"
bibliography: ../inst/REFERENCES.bib
output:
rmarkdown::html_vignette:
  toc: yes
vignette: >
  %\VignetteIndexEntry{Controlling for Variables}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Introduction
This vignette describes how to control for variables. This is a new feature to **BGGM** (version `2.0.0`).

# Example 1: Multivariate Regression
When controlling for variables, a multivariate regression is fitted in **BGGM**. In fact, a GGM can be understood as a multivariate regression with intercepts only models for the predictors.

## Notes about Implementation
**BGGM** does not use the typical approach for multivariate regression in `R`. This avoids having to 
write out each outcome variable, of which there are typically many in a GGM. In **BGGM**, it is assumed
that the data matrix includes only the variables to be included in the GGM and the control variables.

### Correct 
Suppose that we want to control for education level, with five variables included in the 
GGM.

```r
# data
Y <- bfi[,c(1:5, 27)]

# head
head(Y)

#>       A1 A2 A3 A4 A5 education
#> 61617  2  4  3  4  4        NA
#> 61618  2  4  5  2  5        NA
#> 61620  5  4  5  4  4        NA
#> 61621  4  4  6  5  5        NA
#> 61622  2  3  3  4  5        NA
#> 61623  6  6  5  6  5         3
```

Notice that `Y` includes **only** the five variables and `education`.





## Fit Model
This model can then be fitted with

```
fit <- explore(Y, formula = ~ as.factor(education))
```

To show this is indeed a multivariate regression, here are the summarized regression coefficients for the first
outcome.

```
summ_coef <- regression_summary(fit)

# outcome one
summ_coef$reg_summary[[1]]

#>                       Post.mean Post.sd Cred.lb Cred.ub
#> (Intercept)               0.256   0.095   0.072   0.442
#> as.factor(education)2     0.073   0.128  -0.177   0.323
#> as.factor(education)3    -0.202   0.104  -0.405  -0.001
#> as.factor(education)4    -0.462   0.119  -0.691  -0.233
#> as.factor(education)5    -0.578   0.117  -0.815  -0.346
```

And here are the coefficients from `lm` (a univariate regression for `A1`)

```
round(
  cbind(
    # summary: coef and se
    summary( lm(scale(A1, scale = F) ~ as.factor(education), data = Y))$coefficients[,1:2],
    # confidence interval
    confint( lm(scale(A1, scale = F) ~ as.factor(education), data = Y))
), 3)


#>                       Estimate Std. Error  2.5 % 97.5 %
#> (Intercept)              0.256      0.093  0.073  0.438
#> as.factor(education)2    0.072      0.125 -0.172  0.316
#> as.factor(education)3   -0.203      0.101 -0.401 -0.004
#> as.factor(education)4   -0.461      0.116 -0.690 -0.233
#> as.factor(education)5   -0.578      0.115 -0.804 -0.351
```
The estimate are very (very) similar.

## Summary
Note that all the other functions work just the same. For example, the relations controlling for education
are summarized with

```
summary(fit)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: continuous 
#> Analytic: FALSE 
#> Formula: ~ as.factor(education) 
#> Posterior Samples: 5000 
#> Observations (n):
#> Nodes (p): 5 
#> Relations: 10 
#> --- 
#> Call: 
#> estimate(Y = Y, formula = ~as.factor(education))
#> --- 
#> Estimates:
#>  Relation Post.mean Post.sd Cred.lb Cred.ub
#>    A1--A2    -0.239   0.020  -0.278  -0.200
#>    A1--A3    -0.109   0.020  -0.150  -0.070
#>    A2--A3     0.276   0.019   0.239   0.312
#>    A1--A4    -0.013   0.021  -0.055   0.026
#>    A2--A4     0.156   0.020   0.117   0.196
#>    A3--A4     0.173   0.020   0.134   0.214
#>    A1--A5    -0.010   0.020  -0.050   0.029
#>    A2--A5     0.150   0.020   0.111   0.189
#>    A3--A5     0.358   0.018   0.322   0.392
#>    A4--A5     0.121   0.020   0.082   0.159
#> --- 
```


### Incorrect
Now if we wanted to control for education, but also had gender in `Y`, this would be incorrect

```
Y <- bfi[,c(1:5, 26:27)]

head(Y)

#>       A1 A2 A3 A4 A5 gender education
#> 61617  2  4  3  4  4      1        NA
#> 61618  2  4  5  2  5      2        NA
#> 61620  5  4  5  4  4      2        NA
#> 61621  4  4  6  5  5      2        NA
#> 61622  2  3  3  4  5      1        NA
#> 61623  6  6  5  6  5      2         3
```


In this case, with `estimate(Y, formula = as.factor(education))`, the GGM would also include `gender`
(six variables instead of the desired 5). This is because all variables not included in `formula` are included in the GGM. This was adopted in **BGGM** to save the user from having to write out each outcome.  

This differs from `lm`, where each outcome needs to be written out, for example `cbind(A1, A2, A3, A4, A4) ~ as.factor(education)`. This is quite cumbersome for a model that includes many nodes.

# Example 2: Multivariate Probit
The above data is ordinal. In this case, it is possible to fit a multivariate probit model. This is also the approach for binary data in **BGGM**. This is implemented with 

```
fit <- estimate(Y, formula = ~ as.factor(education), 
                type = "ordinal", iter = 1000)
```

Note that the multivariate probit models can also be summarized with `regression_summary`.

# Example 3: Gaussian Copula Graphical Model
This final example fits a Gaussian copula graphical model that can be used for mixed data. In this case,
`formula` is not used and instead all of the variables are included in the GGM. 

## Fit Model

This model is estimated with
```
# data
Y <- na.omit(bfi[,c(1:5, 27)])

# fit type = "mixed"
fit <- estimate(Y, type = "mixed", iter = 1000)

# summary
summary(fit)

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: mixed 
#> Analytic: FALSE 
#> Formula:  
#> Posterior Samples: 1000 
#> Observations (n):
#> Nodes (p): 6 
#> Relations: 15 
#> --- 
#> Call: 
#> estimate(Y = Y, type = "mixed", iter = 1000)
#> --- 
#> Estimates:
#>       Relation Post.mean Post.sd Cred.lb Cred.ub
#>         A1--A2    -0.217   0.048  -0.294  -0.114
#>         A1--A3    -0.063   0.027  -0.113  -0.011
#>         A2--A3     0.364   0.023   0.317   0.410
#>         A1--A4     0.116   0.038   0.048   0.192
#>         A2--A4     0.241   0.031   0.182   0.303
#>         A3--A4     0.228   0.026   0.174   0.275
#>         A1--A5     0.057   0.031   0.003   0.120
#>         A2--A5     0.186   0.027   0.135   0.241
#>         A3--A5     0.438   0.019   0.399   0.474
#>         A4--A5     0.151   0.025   0.103   0.199
#>  A1--education    -0.016   0.069  -0.125   0.119
#>  A2--education     0.063   0.049  -0.016   0.162
#>  A3--education     0.049   0.025   0.002   0.099
#>  A4--education     0.053   0.026   0.005   0.105
#>  A5--education     0.072   0.024   0.024   0.120
#> --- 
```

Here it is clear that education is included in the model, as the relations with the other nodes are included in the output.

## Select Graph
The graph is selected with
```
select(fit)
```


# Note
It is possible to control for variable with all methods in **BGGM**, including when comparing groups, Bayesian hypothesis testing, etc.
---
title: "Custom Network Comparisons"
author: "Donny Williams"
date: "5/19/2020"
bibliography: ../inst/REFERENCES.bib
output: 
  rmarkdown::html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Custom Network Comparisons}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Background
It is quite common to have partial correlation networks (GGMs) for various subgroups, 
say, males and females, a control and treatment group, or perhaps several educational 
levels. In this case, it is important to not only determine whether the groups are 
different, but actually compare the groups in a way that answers a specific question 
of interest. 

To date, most `R` packages provide a few ways to compare groups, including **BGGM** (version `1.0.0`). 
In version `2.0.0`, however, **BGGM** includes a new feature for the function `ggm_compare_ppc` that enables 
users to **compare networks in any way they want**. 

# Basic Idea
The technical details of the approach are described in [@williams2020comparing]. The basic idea is to

1. Draw samples from the posterior distribution, assuming the groups are equal (i.e., the "null" model).

2. Generate the posterior **predictive** distribution for the chosen test-statistic (how the groups 
   are being compared)

    + This can be understood as what we would expect to observe in the future 
    (e.g., in replication), assuming the groups were in fact equal.

3. Compute the test-statistic for the observed groups.

4. Then compare the observed test-statistic to the predictive distribution 
   (what is expected under the "null" model). 

    + If the observed error is larger than the model assuming group equality, this 
      suggests that the groups are different.


In **BGGM**, the default is to compare the groups with respect to (symmetric) Kullback-Leibler 
divergence (i.e., "distance" between multivariate normal distributions)  and the sum of 
squared error (for the partial correlation matrix). This was shown to be quite powerful in @williams2020comparing, while also having a 
low false positive rate. 

In the following, the focus is on defining custom functions
and using them with `ggm_compare_ppc`. In all examples, post-traumatic stress disorder
networks are compared [@fried2018replicability].

### R packages
```{r, eval = FALSE}
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")

```


### Data
Only the correlation matrices are available. Hence, multivariate normal data is generated with that *exact* 
correlation structure via the `R` package **MASS**.
```{r, warning =FALSE, message=FALSE}
# need these packages
library(BGGM)
library(ggplot2)
library(assortnet)
library(networktools)
library(MASS)

# group 1
Yg1 <- MASS::mvrnorm(n = 926, 
                     mu = rep(0, 16), 
                     Sigma = ptsd_cor3, 
                     empirical = TRUE)

# group 2
Yg2 <- MASS::mvrnorm(n = 956, 
                     mu = rep(0, 16), 
                     Sigma = ptsd_cor4, 
                     empirical = TRUE)
```


# Illustrative Examples

## Correlation
This first example looks at the correlation between partial correlations of the two networks. Note that
it could be two networks have what is considered a large correlation. However, the question here is, 
assuming the groups are equal, just how large should the correlation be? This is needed to interpret
the observed test-statistic.


### Step 1: Define Custom Function
The first step is to define a custom function that takes two data matrices and the output 
is the chosen test-statistic (in this case a correlation)
```{r}
f <- function(Yg1, Yg2){
  # number of nodes
  p <- ncol(Yg1)
  
  # index of off-diagonal
  indices <- upper.tri( diag(p))
  
  # group 1:
  # fit model
  g1_fit <- estimate(Yg1, analytic = TRUE)
  # pcors
  g1_pcors <- pcor_mat(g1_fit)[indices]
  
  # group 2
  # fit model
  g2_fit <- estimate(Yg2, analytic = TRUE)
  # pcors
  g2_pcors <- pcor_mat(g2_fit)[indices]
  
  # test-statistic
  cor(g1_pcors, g2_pcors)
  }
```


### Step 2: Compute the Observed Score
The next step is to compute the observed test-statistic, that is, the correlation between the partial correlations.
```{r}
obs <- f(Yg1, Yg2)

# observed
obs
```


### Step 3: Predictive Check
With the function, `f`, and the observed scores, `obs`, in hand, what is left is the predictive check

```{r, message=FALSE, results='hide'}
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000, 
                       loss = FALSE)
```

Note that `loss = FALSE` controls how the p-value is computed. It is an indicator of whether the test-statistic is 
a "loss" (a bad thing). In this case, a large correlation is a good thing so it is set to `FALSE`. The results
can then be printed

```{r}
ppc
```

which shows the posterior predictive p-value is zero. This indicates that the observed correlation is lower than
the entire predictive distribution (the distribution of correlations for future data, assuming group equality)

and finally plot the results
```{r, eval=FALSE}
plot(ppc)
```

```{r, echo=FALSE, message=FALSE, warning=FALSE}
plot(ppc, col_critical = "lightblue", 
     col_noncritical = "lightblue")[[1]] +
  xlab("Predictive Correlation")

```



The density is the predictive distribution for the correlation. Recall that this is the correlation that we would expect, given the groups were actually the same, and the black point is the observed correlation. In this case, it seems quite clear that the "null model" is inadequate--the groups are apparently quite different.

## Hamming Distance
The next example is [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance), which, in this case, is the squared error for the adjacency matrices. It seems reasonable to think of this as a test for 
different network structures or patterns of zeros and ones.


### Step 1: Define Custom Function
The first step is to define a custom function that takes two data matrices and the output 
is the chosen test-statistic (in this case Hamming distance)
```{r}
f <- function(Yg1, Yg2){
  # nodes
  p <- ncol(Yg1)
  
  # index of off-diagonal
  indices <- upper.tri( diag(p))

  # fit models
  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)
  
  # select graphs
  sel1 <- BGGM::select(fit1)
  sel2 <- BGGM::select(fit2)
  
  # hamming distance
  sum((sel1$adj[indices] - sel2$adj[indices]) ^ 2)
}
```


### Step 2: Compute the Observed Score
The next step is to compute the observed test-statistic, that is, the Hamming distance between adjacency matrices

```{r}
obs <- f(Yg1, Yg2)

# observed
obs
```

### Step 3: Predictive Check
With the function, `f`, and the observed scores, `obs`, in hand, what is left is the predictive check
```{r, message=FALSE, results='hide'}
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)
```

The results can then be printed
```{r}
ppc
```

And then plot the results
```{r, message=FALSE, warning=FALSE}
plot(ppc)
```

This result is intriguing. Whereas the correlation looked at the relation between partial correlation, here there seems to be evidence 
that the adjacency matrices are different (perhaps suggesting that the conditional independence structure is different).


## Partial Correlation Matrix Distance
There might also be interest in the so-called correlation matrix distance [@herdin2005correlation]. This is also easily tested, in this case for the partial correlation matrix.

### Step 1: Define Custom Function
```{r}
f <- function(Yg1, Yg2){
  # nodes
  p <- ncol(Yg1)
  
  # index of off-diagonal
  indices <- upper.tri( diag(p))

  # fit models
  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)
  
  pcor1 <- BGGM::pcor_mat(fit1) 
  pcor2 <- BGGM::pcor_mat(fit2)
  
  # CDM for partial correlations
  # note: numerator is the trace; denominator is the Frobenius norm
  1 - (sum(diag(pcor1 %*% pcor2)) / (norm(pcor1, type = "f") * norm(pcor2, type = "f")))
}
```

### Step 2: Compute the Observed Score
The next step is to compute the observed test-statistic, that is, the Partial Correlation Matrix Distance
```{r}
obs <- f(Yg1, Yg2)

# observed
obs
```

### Step 3: Predictive Check
With the function, `f`, and the observed scores, `obs`, in hand, what is left is the predictive check
```{r, message=FALSE, results='hide'}
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)
```

The results can then be printed
```{r}
ppc
```
which again provides a p-value of zero.

Note that the object `ppc` includes the predictive samples that allows for user defined plots (in the event something custom is desired).

```{r}
hist(ppc$predictive_custom, 
     xlim = c(0, obs), 
     main = "Partial Correlation Matrix Distance")
abline(v = obs)
```

Note that the line is the observed which again makes it clear that the distance is quite surprising, 
assuming the null model were true.


## Assortment
This next example is assortment [@newman2003mixing], which is a measure related 
to clustering in a network. Here the test is for a difference in assortment. 
This is computed by taking the difference (absolute value) for each draw 
from the predictive distribution.

### Step 1: Define Custom Function
```{r}
# clusters based on DSM-5
comms <- c(
  rep("A", 4),
  rep("B", 7),
  rep("C", 5)
)

f <- function(Yg1, Yg2){

  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)

  pcor1 <- BGGM::pcor_mat(fit1)
  pcor2 <- BGGM::pcor_mat(fit2)

  assort1 <- assortnet::assortment.discrete(pcor1, types = comms,
                                         weighted = TRUE,
                                         SE = FALSE, M = 1)$r

  assort2 <- assortnet::assortment.discrete(pcor2, types = comms,
                                          weighted = TRUE,
                                          SE = FALSE, M = 1)$r
  (assort1 - assort2)
}
```


### Step 2: Compute the Observed Score
The next step is to compute the observed test-statistic, that is, assortment for the two groups
```{r}
obs <- f(Yg1, Yg2)

# observed
obs
```

### Step 3: Predictive Check
With the function, `f`, and the observed score, `obs`, in hand, the next step is the predictive check

```{r, message=FALSE, results='hide'}
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)
```

The results can then be printed
```{r}
ppc
```

and plotted

```{r}
plot(ppc)
```

which shows that the clustering in the data appears to be different (given the observed value exceeds
the entire predictive distribution).



## Expected Influence
This last example looks at the expected influence for the network [@robinaugh2016identifying]. In this case, the sum of squared error is the test statistic. This is computed from the squared error for each
draw from the predictive distribution.

### Step 1: Define Custom Function
```{r}
f <- function(Yg1, Yg2){

  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)

  pcor1 <- BGGM::pcor_mat(fit1)
  pcor2 <- BGGM::pcor_mat(fit2)

  ei1 <- networktools::expectedInf(pcor1)$step1
 
  ei2 <- networktools::expectedInf(pcor2)$step1
   sum((ei1 - ei2)^2)
}
```


### Step 2: Compute the Observed Score
The next step is to compute the observed test-statistic, that is, the sum of squared error
for expected influence
```{r}
obs <- f(Yg1, Yg2)

# observed
obs
```


### Step 3: Predictive Check
With the function, `f`, and the observed scores, `obs`, in hand, what is left is the predictive check
```{r, message=FALSE, results='hide'}
ppc <- BGGM:::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)
```

The results can then be printed
```{r}
ppc
```

and plotted 

```{r}
hist(ppc$predictive_custom, 
    xlim = c(0, obs),
     main = "Expected Influence\n Sum of Squared Error")
abline(v = obs)
```
which again shows the sum of squared error for expected influence far exceeds what would be expected, assuming
the null model were true.

# Two Notes of Caution

1. Note that only the default in **BGGM** have been shown to have nominal error rates. However, there is a proof that suggests the error rate cannot be larger than $2\alpha$ [@meng1994posterior], and, further, a predictive check is typically below $\alpha$ [i.e., a tendency to be conservative, @gelman2013two]. 

2. Failing to reject the null model does not indicate the groups are the same! To test for 
equality see `ggm_compare_explore` and `ggm_compare_confirm`.

# Conclusion
These example certainly open the door for tailoring network comparison to answer specific research questions.

# References
---
title: "Custom Network Statistics"
author: "Donny Williams"
date: "5/19/2020"
bibliography: ../inst/REFERENCES.bib
output: 
  rmarkdown::html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Custom Network Statistics}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


# Background
This vignette describes a new feature to **BGGM** (`2.0.0`) that allows for 
computing custom network statistics (e.g., centrality). The new function is
called `roll_your_own` and it was suggested by a user of **BGGM**  ([see feature request here](https://github.com/donaldRwilliams/BGGM/issues/12)).

# Basic Idea
The basic idea is to compute the chosen network statistic for each of the sampled partial 
correlation matrices, resulting in a distribution. All that is required is to define a function
that takes either a partial correlation matrix or a weighted adjacency matrix 
(the partial correlation matrix with values set to zero) as the first argument. 
Several examples are provided below.



### R packages
```{r, eval = FALSE}
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")
```


### Data
In all examples, a subset of `ptsd` data is used. The subset includes two of the "communities" of 
symptoms [details for these data can be found in @armour2017network]. The data are ordinal (5-level Likert).
```{r, warning =FALSE, message=FALSE}
# need these packages
library(BGGM)
library(ggplot2)
library(assortnet)
library(networktools)

# data
Y <- ptsd[,1:7]
```

### Fit Model
For these data, the GGM is estimated with a semi-parametric copula [@hoff2007extending]. 
In **BGGM**, this implemented with `type = mixed` which is kind of a misnomer because the data do not 
have to be "mixed" (consisting of continuous and discrete variables). 
Note that the model is fitted only once which highlights that only the posterior samples 
are needed to compute any network statistic.

```{r, message=FALSE, warning=FALSE, eval=FALSE}
library(BGGM)

# copula ggm
fit <- estimate(Y, type = "mixed", iter = 1000)
```


# Examples

## Expected Influence
The first example computes expected influence [@robinaugh2016identifying]. The first step is to define a function
```{r}
# define function
f <- function(x,...){
  networktools::expectedInf(x,...)$step1
}
```

Note that `x` takes the matrix which is then passed to `expectedInf`. The `...` allows for 
passing additional arguments to the `expectedInf` function. An example is provided below. 
With the function defined, the next step is to compute the network statistic.

```{r, eval = FALSE, message=FALSE, results='hide'}
# iter = 250 for demonstrative purposes
# (but note even 1000 iters takes less than 1 second)
# compute
net_stat <- roll_your_own(object = fit,
                          FUN = f,
                          select = FALSE,
                          iter = 250)
# print
net_stat

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Network Stats: Roll Your Own
#> Posterior Samples: 250 
#> --- 
#> Estimates: 
#> 
#>  Node Post.mean Post.sd Cred.lb Cred.ub
#>     1     0.701   0.099   0.508   0.871
#>     2     0.912   0.113   0.722   1.179
#>     3     0.985   0.112   0.742   1.199
#>     4     1.056   0.105   0.851   1.247
#>     5     1.056   0.116   0.862   1.288
#>     6     0.491   0.092   0.329   0.679
#>     7     0.698   0.098   0.521   0.878
#> --- 
```



The option `select = FALSE` indicates to compute the statistics from the partial correlation matrices (nothing set to zero). This can be changed with `select = TRUE`. Internally, each of the sampled
partial correlation matrices is multiplied by the adjacency matrix.
```{r, eval = FALSE, results='hide'}
net_stat <- roll_your_own(object = fit,
                          FUN = f,
                          select = TRUE,
                          iter = 250)

# print
net_stat

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Network Stats: Roll Your Own
#> Posterior Samples: 250 
#> --- 
#> Estimates: 
#> 
#>  Node Post.mean Post.sd Cred.lb Cred.ub
#>     1     0.636   0.136   0.386   0.874
#>     2     0.792   0.113   0.580   0.996
#>     3     0.777   0.122   0.544   1.001
#>     4     0.910   0.121   0.667   1.143
#>     5     0.525   0.104   0.331   0.727
#>     6     0.484   0.110   0.270   0.686
#>     7     0.247   0.081   0.088   0.412
#> --- 
```


The results are then plotted with 
```{r, message=FALSE, eval=FALSE}
plot(net_stat)
```
![](../man/figures/netstat_ridge.png)

## Bridge Strength
The next example computes bridge strength [@jones2019bridge]. This requires the user to define clusters or "communities".

```{r, eval = FALSE, message=FALSE, results='hide'}
# clusters
communities <- substring(colnames(Y), 1, 1)

# function is slow
f <- function(x, ...){
networktools::bridge(x, ...)$`Bridge Strength`
}


# compute
net_stat <- roll_your_own(object = fit,
                          FUN = f, 
                          communities = communities,
                          iter = 250)

# print
net_stat

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Network Stats: Roll Your Own
#> Posterior Samples: 250 
#> --- 
#> Estimates: 
#> 
#>  Node Post.mean Post.sd Cred.lb Cred.ub
#>     1     0.162   0.082   0.035   0.347
#>     2     0.250   0.113   0.061   0.501
#>     3     0.180   0.104   0.049   0.480
#>     4     0.280   0.098   0.090   0.480
#>     5     0.375   0.093   0.196   0.558
#>     6     0.617   0.166   0.339   1.002
#>     7     0.628   0.166   0.400   1.025
#> --- 
```


Notice `communities`. This is passed to `...` in the function `f`, which, in turn, is passed to the function `bridge`. Any number of arguments can be passed this way. Here are the results


This can then be plotted and further customized (the returned object is a `ggplot`)
```{r, message = FALSE, eval=FALSE}
plot(net_stat, 
     fill = "lightblue") + 
  ggtitle("Bridge Strength") + 
  xlab("Score")
```

![](../man/figures/netstat_bridge.png)

## Assortment
The next example computes assortment [@newman2003mixing].


```{r, eval = FALSE, message=FALSE, results='hide'}
# clusters
communities <- substring(colnames(Y), 1, 1)

# define function
f <- function(x,...){
  assortnet::assortment.discrete(x, ...)$r
}

net_stat <- roll_your_own(object = fit,
                          FUN = f,
                          types = communities,
                          weighted = TRUE,
                          SE = FALSE, M = 1, 
                          iter = 250)

# print
net_stat

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Network Stats: Roll Your Own
#> Posterior Samples: 250 
#> --- 
#> Estimates: 
#> 
#>  Post.mean Post.sd Cred.lb Cred.ub
#>      0.261   0.124   -0.01   0.469
#> --- 
```

This example demonstrate that `...` can take several arguments. The results are stored in the `net_stat` object. They can be accessed with 

```{r, eval=FALSE}
hist(net_stat$results, main = "Assortment")
```

![](../man/figures/netstat_assort.png)

# Note
The function `roll_your_own` is expecting the custom function to return either a single number or a number for each node. This ensures all the printing and plotting functions work. However, you could return anything you want and then access the results to plot, summarize, etc.

# References
---
title: "Testing Sums"
author: "Donny Williams"
date: "5/25/2020"
bibliography: ../inst/REFERENCES.bib
output:
rmarkdown::html_vignette:
  toc: yes
vignette: >
  %\VignetteIndexEntry{Testing Sums}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Introduction
This is a follow-up to the vignette ["Three Ways to Test the Same Hypothesis"](https://donaldrwilliams.github.io/BGGM/articles/hyp_3_ways.html). A
new feature, `pcor_sum`, was added to **BGGM** that allows for testing partial correlation sums. 
This differs from the Bayes factor approach ("Approach #3"), in that only the posterior 
distribution is used to determine whether there is a difference in the sums. 

### R package
```{r, eval = FALSE, message=FALSE}
# need the developmental version
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   

# install from github
remotes::install_github("donaldRwilliams/BGGM")
library(BGGM)
```


# One Group
This first example looks at one group, where a sum is tested within the same ptsd network. I focus on the 
relations between the re-experiencing (`B`) and avoidance (`C`) communities. In particular, the sum of relations between the "Intrusion" (5 nodes) community and the "Avoidance" (two nodes) community is tested. 

## Sum to String
For the avoidance symptom "avoidance of thoughts" `C1`, this can be written in `R` code with

```
# ptsd
Y <- ptsd

# paste together sums
paste0(colnames(Y)[1:5],  "--C1", collapse = " + ")

#> "B1--C1 + B2--C1 + B3--C1 + B4--C1 + B5--C1"
```
whereas, for the avoidance symptom "avoidance of reminders" (`C2`), this is written as

```
paste0(colnames(Y)[1:5],  "--C2", collapse = " + ")

#> "B1--C2 + B2--C2 + B3--C2 + B4--C2 + B5--C2"
```
Note that typically this would have to be written out. `paste0` was used in this case to 
avoid typing out all of the relations.

## Fit Model

Here an ordinal GGM is fitted

```
fit <- estimate(Y+1, type = "ordinal", iter = 1000)
```
where the `+1` changes the first category from 0 to 1 (required).

## Test Sums
The next step is to use the `pcor_sum` function. First, I combine the sums into one string separated with `;`.
```
# sum 1
sum1 <- paste0(colnames(Y)[1:5],  "--C1", collapse = " + ")

# sum 2
sum2 <- paste0(colnames(Y)[1:5],  "--C2", collapse = " + ")

# paste together
sums <- paste(sum1, sum2, sep = ";")

# print
sums
#> "B1--C1 + B2--C1 + B3--C1 + B4--C1 + B5--C1;B1--C2 + B2--C2 + B3--C2 + B4--C2 + B5--C2"
```

Next `pcor_sum` is used

```
test_sum <- pcor_sum(fit, relations = sums)

# print
test_sum

# BGGM: Bayesian Gaussian Graphical Models 
# --- 
# Network Stats: Posterior Sum
# Posterior Samples: 1000 
# --- 
# Estimates 
# 
# Sum: 
#                                    Post.mean Post.sd Cred.lb Cred.ub
# B1--C1+B2--C1+B3--C1+B4--C1+B5--C1     0.215   0.096   0.034   0.404
# B1--C2+B2--C2+B3--C2+B4--C2+B5--C2     0.334   0.097   0.145   0.514
# --- 
# 
# Difference:
# B1--C1+B2--C1+B3--C1+B4--C1+B5--C1 - B1--C2+B2--C2+B3--C2+B4--C2+B5--C2 
# 
#  Post.mean Post.sd Cred.lb Cred.ub Prob.greater Prob.less
#     -0.119   0.145  -0.409   0.173        0.205     0.795
# --- 
```

`Prob.greater` is the posterior probability that the first sum is larger than the second sum.

## Plot Results
The object `test_sum` can then be plotted. Note this returns three plots, but only the difference is shown here

```
plot(test_sum)$diff
```

![](../man/figures/test_sum_hist.png)

The histogram is not very smooth in this case because `iter = 1000`, but this of course can be changed.

# Two Groups
This next example is for two groups. The data are called `bfi` and they are in the **BGGM** package. I compare a sum of two relations for questions measuring agreeableness in males and females. The relations tested are as follows


## Sum to String
```r
sums <- c("A3--A4 + A4--A5")
```
where `A1` is "know how to comfort others", `A4` is "love children", and `A5` is "make people feel at ease".

## Fit Models
The next step is to fit the models
```r
# data
Y <- bfi

# males
Y_males <- subset(Y, gender == 1, select = -c(education, gender))[,1:5]

# females
Y_females <- subset(Y, gender == 2, select = -c(education, gender))[,1:5]


fit_female <- estimate(Y_females, seed = 2)

# fit males
fit_male <- estimate(Y_males, seed = 1)
```

## Test Sums
Then test the sum

```r
test_sum <- pcor_sum(fit_female, fit_male, relations = sums)

# print
test_sum

#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Network Stats: Posterior Sum
#> Posterior Samples: 5000 
#> --- 
#> Estimates 
#> 
#> Sum: 
#>                   Post.mean Post.sd Cred.lb Cred.ub
#> g1: A3--A4+A4--A5     0.292   0.026   0.241   0.342
#> g2: A3--A4+A4--A5     0.305   0.036   0.234   0.375
#> --- 
#> 
#> Difference:
#> g1: A3--A4+A4--A5 - g2: A3--A4+A4--A5 
#> 
#>  Post.mean Post.sd Cred.lb Cred.ub Prob.greater Prob.less
#>     -0.014   0.045    -0.1   0.074        0.386     0.614
#> --- 
```
## Sanity Check
For a kind of sanity check, here is the sum for the male group obtained from the point estimates.

```r
pcor_mat(fit_male)["A3", "A4"] + pcor_mat(fit_male)["A4", "A5"] 

#>  0.305
```

This matches the output.

# Notes
By default, the print function for `pcor_sum` provides 95 % credible intervals. This can be changed by
directly using the print function, for example `print(test_sum, cred = 0.99)`, provides
99 % credible intervals.

Currently, this function only supports sums, due to this being of interest for the psychological network
literature in particular. This can be extended to accommodate multiplication, subtraction,
testing values other than zero, etc. Please make a feature request at either
[github](https://github.com/donaldRwilliams/BGGM/issues) or [BGGM-users group](https://groups.google.com/forum/#!forum/bggm-users).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.explore.R
\name{select.explore}
\alias{select.explore}
\title{Graph selection for \code{explore} Objects}
\usage{
\method{select}{explore}(object, BF_cut = 3, alternative = "two.sided", ...)
}
\arguments{
\item{object}{An object of class \code{explore.default}}

\item{BF_cut}{Numeric. Threshold for including an edge (defaults to 3).}

\item{alternative}{A character string specifying the alternative hypothesis. It
must be one of "two.sided" (default), "greater", "less",
or "exhaustive". See note for further details.}

\item{...}{Currently ignored.}
}
\value{
The returned object of class \code{select.explore} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:


\code{alternative = "two.sided"}

 \itemize{

 \item \code{pcor_mat_zero} Selected partial correlation matrix (weighted adjacency).

 \item \code{pcor_mat} Partial correlation matrix (posterior mean).

 \item \code{Adj_10} Adjacency matrix for the selected edges.

 \item \code{Adj_01} Adjacency matrix for which there was
                     evidence for the null hypothesis.
 }

\code{alternative = "greater"} and \code{"less"}

 \itemize{

 \item \code{pcor_mat_zero} Selected partial correlation matrix (weighted adjacency).

 \item \code{pcor_mat} Partial correlation matrix (posterior mean).

 \item \code{Adj_20} Adjacency matrix for the selected edges.

 \item \code{Adj_02} Adjacency matrix for which there was
                     evidence for the null hypothesis (see note).
 }

\code{alternative = "exhaustive"}

\itemize{

\item \code{post_prob} A data frame that included the posterior hypothesis probabilities.

\item \code{neg_mat} Adjacency matrix for which there was evidence for negative edges.

\item \code{pos_mat} Adjacency matrix for which there was evidence for positive edges.

\item \code{neg_mat} Adjacency matrix for which there was
                     evidence for the null hypothesis (see note).

 \item \code{pcor_mat} Partial correlation matrix (posterior mean). The weighted adjacency
 matrices can be computed by multiplying \code{pcor_mat} with an adjacency matrix.

}
}
\description{
Provides the selected graph based on the Bayes factor
\insertCite{Williams2019_bf}{BGGM}.
}
\details{
Exhaustive provides the posterior hypothesis probabilities for
a positive, negative, or null relation \insertCite{@see Table 3 in @Williams2019_bf}{BGGM}.
}
\note{
Care must be taken with the options \code{alternative = "less"} and
      \code{alternative = "greater"}. This is because the full parameter space is not included,
      such, for  \code{alternative = "greater"}, there can be evidence for the "null" when
      the relation is negative. This inference is correct: the null model better predicted
      the data than the positive model. But note this is relative and does \strong{not}
      provide absolute evidence for the null hypothesis.
}
\examples{

\donttest{
#################
### example 1 ###
#################

#  data
Y <- bfi[,1:10]

# fit model
fit <- explore(Y, progress = FALSE)

# edge set
E <- select(fit,
            alternative = "exhaustive")

}
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{explore}} and \code{\link{ggm_compare_explore}} for several examples.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd_cor2}
\alias{ptsd_cor2}
\title{Data: Post-Traumatic Stress Disorder (Sample # 2)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix
was estimated from 365 individuals \insertCite{fried2018replicability}{BGGM}.
}
\details{
\itemize{
  \item Intrusive Thoughts
  \item Nightmares
  \item Flashbacks
  \item Physiological/psychological reactivity
  \item Avoidance of thoughts
  \item Avoidance of situations
  \item Amnesia
  \item Disinterest in activities
  \item Feeling detached
  \item Emotional numbing
  \item Foreshortened future
  \item Sleep problems
  \item Irritability
  \item Concentration problems
  \item Hypervigilance
  \item Startle response
}
}
\examples{
data(ptsd_cor2)
Y <- MASS::mvrnorm(n = 365,
                   mu = rep(0, 16),
                   Sigma = ptsd_cor2,
                   empirical = TRUE)
}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_estimate.default.R
\name{plot.summary.ggm_compare_estimate}
\alias{plot.summary.ggm_compare_estimate}
\title{Plot \code{summary.ggm_compare_estimate} Objects}
\usage{
\method{plot}{summary.ggm_compare_estimate}(x, color = "black", size = 2, width = 0, ...)
}
\arguments{
\item{x}{An object of class \code{ggm_compare_estimate}.}

\item{color}{Character string. The color of the points
(defaults to \code{"black"}).}

\item{size}{Numeric. The size of the points (defaults to 2).}

\item{width}{Numeric. The width of error bar ends (defaults to \code{0}).}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{ggplot}
}
\description{
Visualize the posterior distribution differences.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes
# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
                select = -c(gender,
                            education))[,1:5]

Yfemale <- subset(Y, gender == 2,
                  select = -c(gender,
                              education))[,1:5]

# fit model
fit <- ggm_compare_estimate(Ymale,  Yfemale,
                            type = "ordinal",
                            iter = 250,
                            prior_sd = 0.25,
                            progress = FALSE)

plot(summary(fit))
}

}
\seealso{
\code{\link{ggm_compare_estimate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{depression_anxiety_t1}
\alias{depression_anxiety_t1}
\title{Data: Depression and Anxiety (Time 1)}
\format{
A data frame containing 403 observations (n = 7466) and 16 variables (p = 16) measured on the 4-point
likert scale.
}
\usage{
data("depression_anxiety_t1")
}
\description{
A data frame containing 403 observations (n = 403) and 16 variables (p = 16) measured on the 4-point
likert scale (depression: 9; anxiety: 7).
}
\details{
\strong{Depression}:

\itemize{
  \item \code{PHQ1}  Little interest or pleasure in doing things?
  \item \code{PHQ2}  Feeling down, depressed, or hopeless?
  \item \code{PHQ3}  Trouble falling or staying asleep, or sleeping too much?
  \item \code{PHQ4}  Feeling tired or having little energy?
  \item \code{PHQ5}  Poor appetite or overeating?
  \item \code{PHQ6} Feeling bad about yourself — or that you are a failure or have let
                    yourself or your family down?
  \item \code{PHQ7}  Trouble concentrating on things, such as reading the newspaper or
                     watching television?
  \item \code{PHQ8} Moving or speaking so slowly that other people could have noticed? Or so
                    fidgety or restless that you have been moving a lot more than usual?
  \item \code{PHQ9}  Thoughts that you would be better off dead, or thoughts of hurting yourself
                     in some way?
}

  \strong{Anxiety}
  \itemize{



  \item \code{GAD1} Feeling nervous, anxious, or on edge
  \item \code{GAD2} Not being able to stop or control worrying
  \item \code{GAD3} Worrying too much about different things
  \item \code{GAD4} Trouble relaxing
  \item \code{GAD5} Being so restless that it's hard to sit still
  \item \code{GAD6} Becoming easily annoyed or irritable
  \item \code{GAD7} Feeling afraid as if something awful might happen
}
}
\examples{
data("depression_anxiety_t1")
labels<- c("interest", "down", "sleep",
            "tired", "appetite", "selfest",
           "concen", "psychmtr", "suicid",
           "nervous", "unctrworry", "worrylot",
           "relax", "restless", "irritable", "awful")


}
\references{
Forbes, M. K., Baillie, A. J., & Schniering, C. A. (2016). A structural equation modeling
analysis of the relationships between depression,anxiety, and sexual problems over time.
The Journal of Sex Research, 53(8), 942-954.

Forbes, M. K., Wright, A. G., Markon, K. E., & Krueger, R. F. (2019). Quantifying the reliability and replicability of psychopathology network characteristics.
Multivariate behavioral research, 1-19.

Jones, P. J., Williams, D. R., & McNally, R. J. (2019). Sampling variability is not nonreplication:
a Bayesian reanalysis of Forbes, Wright, Markon, & Krueger.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore.default.R
\name{summary.explore}
\alias{summary.explore}
\title{Summary Method for \code{explore.default} Objects}
\usage{
\method{summary}{explore}(object, col_names = TRUE, ...)
}
\arguments{
\item{object}{An object of class \code{estimate}}

\item{col_names}{Logical. Should the summary include the column names (default is \code{TRUE})?
Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).}

\item{...}{Currently ignored}
}
\value{
A dataframe containing the summarized posterior distributions.
}
\description{
Summarize the posterior distribution for each partial correlation
with the posterior mean and standard deviation.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

Y <- ptsd[,1:5]

fit <- explore(Y, iter = 250,
               progress = FALSE)

summ <- summary(fit)

summ
}
}
\seealso{
\code{\link{select.explore}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore.default.R
\name{explore}
\alias{explore}
\title{GGM: Exploratory Hypothesis Testing}
\usage{
explore(
  Y,
  formula = NULL,
  type = "continuous",
  mixed_type = NULL,
  analytic = FALSE,
  prior_sd = 0.25,
  iter = 5000,
  progress = TRUE,
  impute = FALSE,
  seed = 1,
  ...
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{formula}{An object of class \code{\link[stats]{formula}}. This allows for including
control variables in the model (i.e., \code{~ gender}).}

\item{type}{Character string. Which type of data for \code{Y} ? The options include \code{continuous},
\code{binary}, \code{ordinal}, or \code{mixed} (semi-parametric copula). See the note for further details.}

\item{mixed_type}{Numeric vector. An indicator of length p for which varibles should be treated as ranks.
(1 for rank and 0 to assume normality). The default is to treat all integer variables as ranks
when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{analytic}{Logical. Should the analytic solution be computed (default is \code{FALSE})?
(currently not implemented)}

\item{prior_sd}{Scale of the prior distribution, approximately the standard deviation
of a beta distribution (defaults to 0.25).}

\item{iter}{Number of iterations (posterior samples; defaults to 5000).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{impute}{Logicial. Should the missing values (\code{NA})
be imputed during model fitting (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed.}

\item{...}{Currently ignored (leave empty).}
}
\value{
The returned object of class \code{explore} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:

\itemize{

\item \code{pcor_mat} partial correltion matrix (posterior mean).

\item \code{post_samp} an object containing the posterior samples.

}
}
\description{
Learn the conditional (in)dependence structure with the Bayes factor using the matrix-F
prior distribution \insertCite{Mulder2018}{BGGM}. These methods were introduced in
\insertCite{Williams2019_bf;textual}{BGGM}. The graph is selected with \code{\link{select.explore}} and
then plotted with \code{\link{plot.select}}.
}
\details{
\strong{Controlling for Variables}:

When controlling for variables, it is assumed that \code{Y} includes \emph{only}
the nodes in the GGM and the control variables. Internally, \code{only} the predictors
that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
\code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
should be included in the GGM. An example is provided below.

\strong{Mixed Type}:

 The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
 continuous or \emph{only} discrete variables. This is based on the ranked likelihood which requires sampling
 the ranks for each variable (i.e., the data is not merely transformed to ranks). This is computationally
 expensive when there are many levels. For example, with continuous data, there are as many ranks
 as data points!

 The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
 and the "emprical" distribution is used otherwise. This is accomplished by specifying an indicator
 vector of length \emph{p}. A one indicates to use the ranks, whereas a zero indicates to "ignore"
 that variable. By default all integer variables are handled as ranks.

\strong{Dealing with Errors}:

An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):

\itemize{

\item The first is due to sampling the thresholds, especially when the data is heavily skewed.
      This can result in an ill-defined matrix. If this occurs, we recommend to first try
      decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
      change the data type to \code{type = mixed} which then estimates a copula GGM
      (this method can be used for data containing \strong{only} ordinal variable). This should
      work without a problem.

\item  The second is due to how the ordinal data are categorized. For example, if the error states
       that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
       the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.

}


\strong{Imputing Missing Values}:

Missing values are imputed with the approach described in \insertCite{hoff2009first;textual}{BGGM}.
The basic idea is to impute the missing values with the respective posterior pedictive distribution,
given the observed data, as the model is being estimated. Note that the default is \code{TRUE},
but this ignored when there are no missing values. If set to \code{FALSE}, and there are missing
values, list-wise deletion is performed with \code{na.omit}.
}
\note{
\strong{Posterior Uncertainty}:

A key feature of \bold{BGGM} is that there is a posterior distribution for each partial correlation.
This readily allows for visiualizing uncertainty in the estimates. This feature works
with all data types and is accomplished by plotting the summary of the \code{explore} object
(i.e., \code{plot(summary(fit))}). Note that in contrast to \code{estimate} (credible intervals),
the posterior standard deviation is plotted for \code{explore} objects.


\strong{"Default" Prior}:

 In Bayesian statistics, a default Bayes factor needs to have several properties. I refer
 interested users to \insertCite{@section 2.2 in @dablander2020default;textual}{BGGM}. In
 \insertCite{Williams2019_bf;textual}{BGGM}, some of these propteries were investigated including
 model selection consistency. That said, we would not consider this a "default" (or "automatic")
 Bayes factor and thus we encourage users to perform sensitivity analyses by varying
 the scale of the prior distribution.

 Furthermore, it is important to note there is no "correct" prior and, also, there is no need
 to entertain the possibility of a "true" model. Rather, the Bayes factor can be interpreted as
 which hypothesis best (\strong{relative} to each other) predicts the observed data
 \insertCite{@Section 3.2 in @Kass1995}{BGGM}.

\strong{Interpretation of Conditional (In)dependence Models for Latent Data}:

See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
(i.e, all data types besides \code{"continuous"})
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

###########################
### example 1:  binary ####
###########################
Y <- women_math[1:500,]

# fit model
fit <- explore(Y, type = "binary",
                iter = 250,
                progress = FALSE)

# summarize the partial correlations
summ <- summary(fit)

# plot the summary
plt_summ <- plot(summary(fit))

# select the graph
E <- select(fit)

# plot the selected graph
plt_E <- plot(E)

plt_E$plt_alt
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.estimate.R
\name{coef.estimate}
\alias{coef.estimate}
\title{Compute Regression Parameters for \code{estimate} Objects}
\usage{
\method{coef}{estimate}(object, iter = NULL, progress = TRUE, ...)
}
\arguments{
\item{object}{An Object of class \code{estimate}}

\item{iter}{Number of iterations (posterior samples; defaults to the number in the object).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{coef}, containting two lists.


\itemize{

\item \code{betas} A list of length \emph{p}, each containing a \emph{p} - 1 by \code{iter} matrix of
posterior samples

\item \code{object} An object of class \code{estimate} (the fitted model).
}
}
\description{
There is a direct correspondence between the inverse covariance matrix and
 multiple regression \insertCite{kwan2014regression,Stephens1998}{BGGM}. This readily allows
 for converting the GGM parameters to regression coefficients. All data types are supported.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

#########################
### example 1: binary ###
#########################
# data
Y <- women_math[1:500, ]

# fit model
fit <- estimate(Y, type = "binary",
                iter = 250,
                progress = FALSE)

# summarize the partial correlations
reg <- coef(fit, progress = FALSE)

# summary
summ <- summary(reg)

summ
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_samples.R
\name{posterior_samples}
\alias{posterior_samples}
\title{Extract Posterior Samples}
\usage{
posterior_samples(object, ...)
}
\arguments{
\item{object}{an object of class \code{estimate} or \code{explore}.}

\item{...}{currently ignored.}
}
\value{
A matrix of posterior samples for the partial correlation. Note that if controlling for
        variables (e.g., formula \code{~ age}), the matrix also includes the coefficients from each
        multivariate regression.
}
\description{
Extract posterior samples for all parameters.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

########################################
### example 1: control  with formula ###
########################################
# (the following works with all data types)

# controlling for gender
Y <- bfi

# to control for only gender
# (remove education)
Y <- subset(Y, select = - education)

# fit model
fit <- estimate(Y, formula = ~ gender,
                iter = 250)

# note regression coefficients
samps <- posterior_samples(fit)

hist(samps[,1])
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictability.R
\name{plot.predictability}
\alias{plot.predictability}
\title{Plot \code{predictability} Objects}
\usage{
\method{plot}{predictability}(
  x,
  type = "error_bar",
  cred = 0.95,
  alpha = 0.5,
  scale = 1,
  width = 0,
  size = 1,
  color = "blue",
  ...
)
}
\arguments{
\item{x}{An object of class \code{predictability}}

\item{type}{Character string. Which type of plot ? The options
are \code{"error_bar"} or \code{"ridgeline"} (defaults to \code{"error_bar"}).}

\item{cred}{Numeric. The credible interval width for summarizing the posterior
distributions (defaults to 0.95; must be between 0 and 1).}

\item{alpha}{Numeric. Transparancey of the ridges}

\item{scale}{Numeric. This controls the overlap of densities
for \code{type = "ridgeline"} (defaults to 1).}

\item{width}{Numeric. The width of error bar ends (defaults to \code{0})
for \code{type = "error_bar"}.}

\item{size}{Numeric. The size for the points (defaults to \code{2})
for \code{type = "error_bar"}.}

\item{color}{Character string. What color for the point (\code{type = "error_bar"}) or
tail region (\code{type = "ridgeline"} ) ? Defaults to \code{"blue"}.}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{ggplot}.
}
\description{
Plot \code{predictability} Objects
}
\examples{
\donttest{
Y <- ptsd[,1:5]

fit <- explore(Y, iter = 250,
               progress = FALSE)

r2 <- predictability(fit, iter = 250,
                     progress = FALSE)

plot(r2)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{women_math}
\alias{women_math}
\title{Data: Women and Mathematics}
\format{
A data frame containing 1190 observations (n = 1190) and 6 variables (p = 6) measured on the binary scale
        \insertCite{fowlkes1988evaluating}{BGGM}. These data have been analyzed in \insertCite{tarantola2004mcmc;textual}{BGGM}
        and in \insertCite{madigan1994model}{BGGM}. The variable descriptions were copied from  (section 5.2 )
        \insertCite{@section 5.2, @talhouk2012efficient}{BGGM}
}
\usage{
data("women_math")
}
\description{
A data frame containing 1190 observations (n = 1190) and 6 variables (p = 6) measured on the binary scale.
}
\details{
\itemize{
  \item \code{1}  Lecture attendance (attend/did not attend)
  \item \code{2}  Gender (male/female)
  \item \code{3}  School type (urban/suburban)
  \item \code{4}  “I will be needing Mathematics in my future work” (agree/disagree)
  \item \code{5}  Subject preference (math/science vs. liberal arts)
  \item \code{6} Future plans (college/job)
}
}
\examples{
data("women_math")
}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_estimate.default.R
\name{summary.ggm_compare_estimate}
\alias{summary.ggm_compare_estimate}
\title{Summary method for \code{ggm_compare_estimate} objects}
\usage{
\method{summary}{ggm_compare_estimate}(object, col_names = TRUE, cred = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{ggm_compare_estimate}.}

\item{col_names}{Logical. Should the summary include the column names (default is \code{TRUE})?
Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).}

\item{cred}{Numeric. The credible interval width for summarizing the posterior
distributions (defaults to 0.95; must be between 0 and 1).}

\item{...}{Currently ignored.}
}
\value{
A list containing the summarized posterior distributions.
}
\description{
Summarize the posterior distribution of each partial correlation
difference with the posterior mean and standard deviation.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes
# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
                select = -c(gender,
                            education))[,1:5]

Yfemale <- subset(Y, gender == 2,
                  select = -c(gender,
                              education))[,1:5]

# fit model
fit <- ggm_compare_estimate(Ymale,  Yfemale,
                            type = "ordinal",
                            iter = 250,
                            prior_sd = 0.25,
                            progress = FALSE)

summary(fit)
}
}
\seealso{
\code{\link{ggm_compare_estimate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_prob.R
\name{predicted_probability}
\alias{predicted_probability}
\title{Predicted Probabilities}
\usage{
predicted_probability(object, outcome, Y, ...)
}
\arguments{
\item{object}{An object of class \code{posterior_predict}}

\item{outcome}{Character string. Node for which the probabilities are computed.}

\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
This must include the column names.}

\item{...}{Compute conditional probabilities by specifying a column name in \code{Y}
(besides the \code{outcome}) and a fixed value. This can include
any number of nodes. See example below. Leave this blank to compute
unconditional probabilities for \code{outcome}.}
}
\value{
A list containing a matrix with the computed probabilities
       (a row for each predictive sample and a column for each category).
}
\description{
Compute the predicted probabilities for discrete data, with the possibility
of conditional predictive probabilities (i.e., at fixed values of other nodes)
}
\note{
There are no checks that the conditional probability exists, i.e., suppose
      you wish to condition on, say, B3 = 2 and B4 = 1, yet there is no instance in
      which B3 is 2 AND B4 is 1. This will result in an uninformative error.
}
\examples{
\donttest{
Y <- ptsd
fit <- estimate(as.matrix(Y), iter = 150, type = "mixed")

pred <- posterior_predict(fit, iter = 100)

prob <- predicted_probability(pred,
                              Y = Y,
                              outcome = "B3",
                              B4 = 0,
                              B5 = 0)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{Sachs}
\alias{Sachs}
\title{Data: Sachs Network}
\format{
A data frame containing 7466 cells (n = 7466) and flow cytometry
 measurements of 11 (p = 11) phosphorylated proteins and phospholipids

 @references
 Sachs, K., Gifford, D., Jaakkola, T., Sorger, P., & Lauffenburger, D. A. (2002).
 Bayesian network approach to cell signaling pathway modeling. Sci. STKE, 2002(148), pe38-pe38.
}
\usage{
data("Sachs")
}
\description{
Protein expression in human immune system cells
}
\examples{
data("Sachs")

}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd_cor4}
\alias{ptsd_cor4}
\title{Data: Post-Traumatic Stress Disorder  (Sample # 4)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix
was estimated from 965 individuals \insertCite{fried2018replicability}{BGGM}.
}
\details{
\itemize{
  \item Intrusive Thoughts
  \item Nightmares
  \item Flashbacks
  \item Physiological/psychological reactivity
  \item Avoidance of thoughts
  \item Avoidance of situations
  \item Amnesia
  \item Disinterest in activities
  \item Feeling detached
  \item Emotional numbing
  \item Foreshortened future
  \item Sleep problems
  \item Irritability
  \item Concentration problems
  \item Hypervigilance
  \item Startle response
}
}
\examples{
data(ptsd_cor4)
Y <- MASS::mvrnorm(n = 965,
                   mu = rep(0, 16),
                   Sigma = ptsd_cor4,
                   empirical = TRUE)

}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.select.R
\name{plot.select}
\alias{plot.select}
\title{Network Plot for \code{select} Objects}
\usage{
\method{plot}{select}(
  x,
  layout = "circle",
  pos_col = "#009E73",
  neg_col = "#D55E00",
  node_size = 10,
  edge_magnify = 1,
  groups = NULL,
  palette = "Set3",
  ...
)
}
\arguments{
\item{x}{An object of class \code{select}.}

\item{layout}{Character string. Which graph layout (defaults is \code{circle}) ?
See \link[sna]{gplot.layout}.}

\item{pos_col}{Character string. Color for the positive edges (defaults to \code{green}).}

\item{neg_col}{Character string.  Color for the negative edges (defaults to \code{green}).}

\item{node_size}{Numeric. The size of the nodes (defaults to \code{10}).}

\item{edge_magnify}{Numeric. A value that is multiplied by the edge weights. This increases (> 1) or
decrease (< 1) the line widths (defaults to 1).}

\item{groups}{A character string of length \emph{p} (the number of nodes in the model).
This indicates groups of nodes that should be the same color
(e.g., "clusters" or "communities").}

\item{palette}{A character string sepcifying the palette for the \code{groups}.
(default is \code{Set3}). See \href{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}{palette options here}.}

\item{...}{Additional options passed to \link[GGally]{ggnet2}}
}
\value{
An object (or list of objects) of class \code{ggplot}
that can then be further customized.
}
\description{
Visualize the conditional (in)dependence structure.
}
\note{
A more extensive example of a custom plot is
provided \href{https://donaldrwilliams.github.io/BGGM/articles/netplot.html}{here}
}
\examples{
\donttest{
#########################
### example 1: one ggm ##
#########################

# data
Y <- bfi[,1:25]

# estimate
fit <- estimate(Y, iter = 250,
                progress = FALSE)

# "communities"
comm <- substring(colnames(Y), 1, 1)

# edge set
E <- select(fit)

# plot edge set
plt_E <- plot(E, edge_magnify = 5,
              palette = "Set1",
              groups = comm)


#############################
### example 2: ggm compare ##
#############################
# compare males vs. females

# data
Y <- bfi[,1:26]

Ym <- subset(Y, gender == 1,
             select = -gender)

Yf <- subset(Y, gender == 2,
              select = -gender)

# estimate
fit <- ggm_compare_estimate(Ym, Yf, iter = 250,
                            progress = FALSE)

# "communities"
comm <- substring(colnames(Ym), 1, 1)

# edge set
E <- select(fit)

# plot edge set
plt_E <- plot(E, edge_magnify = 5,
              palette = "Set1",
              groups = comm)


}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convergence.R
\name{convergence}
\alias{convergence}
\title{MCMC Convergence}
\usage{
convergence(object, param = NULL, type = "trace", print_names = FALSE)
}
\arguments{
\item{object}{An object of class \code{estimate} or \code{explore}}

\item{param}{Character string. Names of parameters for which to monitor MCMC convergence.}

\item{type}{Character string. Which type of convergence plot ? The current
options are \code{trace} (default) and \code{acf}.}

\item{print_names}{Logical. Should the parameter names be printed (defaults to \code{FALSE})? This
can be used to first determine the parameter names to specify in \code{type}.}
}
\value{
A list of \code{ggplot} objects.
}
\description{
Monitor convergence of the MCMC algorithms.
}
\note{
An overview of MCMC diagnostics can be found \href{http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html}{here}.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- ptsd[,1:5]

#########################
###### continuous #######
#########################
fit <- estimate(Y, iter = 250,
                progress = FALSE)

# print names first
convergence(fit, print_names = TRUE)

# trace plots
convergence(fit, type = "trace",
            param = c("B1--B2", "B1--B3"))[[1]]

# acf plots
convergence(fit, type = "acf",
            param = c("B1--B2", "B1--B3"))[[1]]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_ppc.default.R
\name{plot.ggm_compare_ppc}
\alias{plot.ggm_compare_ppc}
\title{Plot \code{ggm_compare_ppc} Objects}
\usage{
\method{plot}{ggm_compare_ppc}(
  x,
  critical = 0.05,
  col_noncritical = "#84e184A0",
  col_critical = "red",
  point_size = 2,
  ...
)
}
\arguments{
\item{x}{An object of class \code{ggm_compare_ppc}}

\item{critical}{Numeric. The 'significance' level
(defaults to \code{0.05}).}

\item{col_noncritical}{Character string. Fill color for the non-critical region
(defaults to \code{"#84e184A0"}).}

\item{col_critical}{Character string. Fill color for the critical region
(defaults to \code{"red"}).}

\item{point_size}{Numeric. The point size for the observed score
(defaults to \code{2}).}

\item{...}{Currently ignored.}
}
\value{
An object (or list of objects) of class \code{ggplot}.
}
\description{
Plot the predictive check with \code{\link[ggridges]{ggridges}}
}
\note{
See
\href{https://CRAN.R-project.org/package=ggridges/vignettes/introduction.html}{ggridges} for
many examples.
}
\examples{
\donttest{
# data
Y <- bfi

#############################
######### global ############
#############################
# males
Ym <- subset(Y, gender == 1,
             select = - c(gender, education))

# females

Yf <- subset(Y, gender == 2,
             select = - c(gender, education))


global_test <- ggm_compare_ppc(Ym, Yf,
                               iter = 250,
                               progress = FALSE)

plot(global_test)
}
}
\seealso{
\code{\link{ggm_compare_ppc}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_pcors.R
\name{gen_net}
\alias{gen_net}
\title{Simulate a Partial Correlation Matrix}
\usage{
gen_net(p = 20, edge_prob = 0.3, lb = 0.05, ub = 0.3)
}
\arguments{
\item{p}{number of variables (nodes)}

\item{edge_prob}{connectivity}

\item{lb}{lower bound for the partial correlations}

\item{ub}{upper bound for the partial correlations}
}
\value{
A list containing the following:

\itemize{

\item{\strong{pcor}}: Partial correlation matrix, encoding
the conditional (in)dependence structure.

\item{\strong{cors}}: Correlation matrix.

\item{\strong{adj}}: Adjacency matrix.

\item{\strong{trys}}: Number of attempts to obtain a
positive definite matrix.

}
}
\description{
Simulate a Partial Correlation Matrix
}
\note{
The function checks for a valid matrix (positive definite),
but sometimes this will still fail. For example, for
larger \code{p}, to have large partial correlations this
requires a sparse GGM
(accomplished by setting \code{edge_prob}
to a small value).
}
\examples{

true_net <- gen_net(p = 10)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd_cor3}
\alias{ptsd_cor3}
\title{Data: Post-Traumatic Stress Disorder  (Sample # 3)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix
was estimated from 926 individuals \insertCite{fried2018replicability}{BGGM}.
}
\details{
\itemize{

  \item Intrusive Thoughts
  \item Nightmares
  \item Flashbacks
  \item Physiological/psychological reactivity
  \item Avoidance of thoughts
  \item Avoidance of situations
  \item Amnesia
  \item Disinterest in activities
  \item Feeling detached
  \item Emotional numbing
  \item Foreshortened future
  \item Sleep problems
  \item Irritability
  \item Concentration problems
  \item Hypervigilance
  \item Startle response
}
}
\examples{
data(ptsd_cor3)
Y <- MASS::mvrnorm(n = 926,
                   mu = rep(0, 16),
                   Sigma = ptsd_cor3,
                   empirical = TRUE)

}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_sum.R
\name{pcor_sum}
\alias{pcor_sum}
\title{Partial Correlation Sum}
\usage{
pcor_sum(..., iter = NULL, relations)
}
\arguments{
\item{...}{An object of class \code{estimate}. This can be either one or two fitted objects.}

\item{iter}{Number of iterations (posterior samples; defaults to the number in the object).}

\item{relations}{Character string. Which partial correlations should be summed?}
}
\value{
An object of class \code{posterior_sum}, including the sum and possibly the difference for
two sums.
}
\description{
Compute and test partial correlation sums either within or between GGMs
(e.g., different groups), resulting in a posterior distribution.
}
\details{
Some care must be taken when writing the string for \code{partial_sum}. Below are several examples

\strong{Just a Sum}:
Perhaps a sum is of interest, and not necessarily the difference of two sums. This can be written as

\itemize{
\item \code{partial_sum <-  c("A1--A2 + A1--A3 + A1--A4")}
}

which will sum those relations.

\strong{Comparing Sums}:
When comparing sums, each must be seperated by "\code{;}". For example,

\itemize{
\item \code{partial_sum <-  c("A1--A2 + A1--A3; A1--A2 + A1--A4")}
}

which will sum both and compute the difference. Note that there cannot be more than two sums, such
that \code{c("A1--A2 + A1--A3; A1--A2 + A1--A4; A1--A2 + A1--A5")} will result in an error.

\strong{Comparing Groups}:

When more than one fitted object is suppled to \code{object} it is assumed that the groups
should be compared for the same sum. Hence, in this case, only the sum needs to be written.

\itemize{
\item \code{partial_sum <-  c("A1--A2 + A1--A3 + A1--A4")}
}

The above results in that sum being computed for each group and then compared.
}
\examples{
\donttest{
# data
Y <- bfi

# males
Y_males <- subset(Y, gender == 1, select = -c(education, gender))[,1:5]

# females
Y_females <- subset(Y, gender == 2, select = -c(education, gender))[,1:5]

# males
fit_males <- estimate(Y_males, seed = 1,
                      progress = FALSE)

# fit females
fit_females <- estimate(Y_females, seed = 2,
                        progress = FALSE)


sums <- pcor_sum(fit_males,
                 fit_females,
                 relations = "A1--A2 + A1--A3")
# print
sums

# plot difference
plot(sums)[[3]]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_confirm.R
\name{plot.confirm}
\alias{plot.confirm}
\title{Plot \code{confirm} objects}
\usage{
\method{plot}{confirm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{confirm}}

\item{...}{Currently ignored.}
}
\value{
A \code{ggplot} object.
}
\description{
Plot the posterior hypothesis probabilities as a pie chart, with
each slice corresponding the probability of a given hypothesis.
}
\examples{

\donttest{

#####################################
##### example 1: many relations #####
#####################################

# data
Y <- bfi

hypothesis <- c("g1_A1--A2 > g2_A1--A2 & g1_A1--A3 = g2_A1--A3;
                 g1_A1--A2 = g2_A1--A2 & g1_A1--A3 = g2_A1--A3;
                 g1_A1--A2 = g2_A1--A2 = g1_A1--A3 = g2_A1--A3")

Ymale   <- subset(Y, gender == 1,
                  select = -c(education,
                              gender))[,1:5]


# females
Yfemale <- subset(Y, gender == 2,
                     select = -c(education,
                                 gender))[,1:5]

test <- ggm_compare_confirm(Ymale,
                            Yfemale,
                            hypothesis = hypothesis,
                            iter = 250,
                            progress = FALSE)


# plot
plot(test)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_belief_var.R
\name{prior_belief_var}
\alias{prior_belief_var}
\title{Prior Belief Graphical VAR}
\usage{
prior_belief_var(
  Y,
  prior_temporal = NULL,
  post_odds_cut = 3,
  est_ggm = TRUE,
  prior_ggm = NULL,
  progress = TRUE,
  ...
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n}
(observations) by \emph{p} (variables/nodes).}

\item{prior_temporal}{Matrix of dimensions \emph{p} by \emph{p},
encoding the prior odds for including each relation
in the temporal graph (see '\code{Details}'). If null
a matrix of 1's is used, resulting in equal prior odds.}

\item{post_odds_cut}{Numeric. Threshold for including an edge (defaults to 3).
Note \code{post_odds} refers to posterior odds.}

\item{est_ggm}{Logical. Should the contemporaneous network be estimated
(defaults to \code{TRUE})?}

\item{prior_ggm}{Matrix of dimensions \emph{p} by \emph{p}, encoding the prior
odds for including each relation in the graph
(see '\code{Details}'). If null a matrix of 1's is used,
resulting in equal prior odds.}

\item{progress}{Logical. Should a progress bar be included
(defaults to \code{TRUE}) ?}

\item{...}{Additional arguments passed to \code{\link{explore}}. Ignored
if \code{prior_ggm = FALSE}.}
}
\value{
An object including (\code{est_ggm = FALSE}):

\itemize{

\item{\strong{adj}}: Adjacency matrix

\item{\strong{post_prob}}: Posterior probability for the
                           alternative hypothesis.

}

An object including (\code{est_ggm = TRUE}):
\itemize{

\item{\strong{adj_temporal}}: Adjacency matrix for the temporal network.

\item{\strong{post_prob_temporal}}: Posterior probability for the
                                    alternative hypothesis (temporal edge)

\item{\strong{adj_ggm}}: Adjacency matrix for the contemporaneous
                         network (ggm).

\item{post_prob_ggm}: Posterior probability for the
                      alternative hypothesis (contemporaneous edge)
}
}
\description{
Prior Belief Graphical VAR
}
\details{
Technically, the prior odds is not for including an edge in the graph,
but for (H1)/p(H0), where H1 captures the hypothesized edge size and H0 is the
null model  \insertCite{@see Williams2019_bf}{BGGM}. Accordingly, setting an
entry in \code{prior_ggm} to, say, 10, encodes a prior belief that H1 is 10 times
more likely than H0. Further, setting an entry in \code{prior_ggm} or
\code{prior_var} to 1 results in equal prior odds
(the default in \code{\link{select.explore}}).
}
\note{
The returned matrices are formatted with the rows indicating
the outcome and the columns the predictor. Hence, adj_temporal[1,4] is the temporal
relation of node 4 predicting node 1. This follows the convention of the
\strong{vars} package (i.e., \code{Acoef}).

Further, in order to compute the Bayes factor the data is
standardized (mean = 0 and standard deviation = 1).
}
\examples{
\donttest{
# affect data from 1 person
# (real data)
y <- na.omit(subset(ifit, id == 1)[,2:7])
p <- ncol(y)

# random prior graph
# (dont do this in practice!!)
prior_var = matrix(sample(c(1,10),
                   size = p^2, replace = TRUE),
                   nrow = p, ncol = p)

# fit model
fit <- prior_belief_var(y,
                        prior_temporal = prior_var,
                        post_odds_cut = 3)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.ggm_compare_bf.R
\name{select.ggm_compare_explore}
\alias{select.ggm_compare_explore}
\title{Graph selection for \code{ggm_compare_explore} Objects}
\usage{
\method{select}{ggm_compare_explore}(object, BF_cut = 3, ...)
}
\arguments{
\item{object}{An object of class \code{ggm_compare_explore}.}

\item{BF_cut}{Numeric. Threshold for including an edge (defaults to 3).}

\item{...}{Currently ignored.}
}
\value{
The returned object of class \code{select.ggm_compare_explore} contains
a lot of information that is used for printing and plotting the results.
For users of \strong{BGGM}, the following are the useful objects:


\itemize{

\item \code{adj_10} Adjacency matrix for which there was evidence for a difference.

\item \code{adj_10} Adjacency matrix for which there was evidence for a null relation

\item \code{pcor_mat_10} Selected partial correlation matrix (weighted adjacency; only for two groups).

}
}
\description{
Provides the selected graph (of differences) based on the Bayes factor
\insertCite{williams2020comparing}{BGGM}.
}
\examples{
\donttest{

##################
### example 1: ###
##################
# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
                   select = -c(gender,
                               education))[,1:10]

Yfemale <- subset(Y, gender == 2,
                     select = -c(gender,
                                 education))[,1:10]

# fit model
fit <- ggm_compare_explore(Ymale, Yfemale,
                           iter = 250,
                           type = "continuous",
                           progress = FALSE)


E <- select(fit, post_prob = 0.50)

}

}
\seealso{
\code{\link{explore}} and \code{\link{ggm_compare_explore}} for several examples.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_bf.default.R
\name{ggm_compare_explore}
\alias{ggm_compare_explore}
\title{GGM Compare: Exploratory Hypothesis Testing}
\usage{
ggm_compare_explore(
  ...,
  formula = NULL,
  type = "continuous",
  mixed_type = NULL,
  analytic = FALSE,
  prior_sd = 0.2,
  iter = 5000,
  progress = TRUE,
  seed = 1
)
}
\arguments{
\item{...}{At least two matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{formula}{An object of class \code{\link[stats]{formula}}. This allows for including
control variables in the model (i.e., \code{~ gender}).}

\item{type}{Character string. Which type of data for \code{Y} ? The options include \code{continuous},
\code{binary}, or \code{ordinal}. See the note for further details.}

\item{mixed_type}{Numeric vector. An indicator of length p for which varibles should be treated as ranks.
(1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{analytic}{logical. Should the analytic solution be computed (default is \code{FALSE}) ? See note for details.}

\item{prior_sd}{Numeric. The scale of the prior distribution (centered at zero), in reference to a beta distribtuion.
The `default` is 0.20. See note for further details.}

\item{iter}{number of iterations (posterior samples; defaults to 5000).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed.}
}
\value{
The returned object of class \code{ggm_compare_explore} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:

\itemize{

\item \code{BF_01} A \emph{p} by \emph{p} matrix including
                    the Bayes factor for the null hypothesis.

\item \code{pcor_diff} A \emph{p} by \emph{p} matrix including
                       the difference in partial correlations (only for two groups).

\item \code{samp} A list containing the fitted models (of class \code{explore}) for each group.

}
}
\description{
Compare Gaussian graphical models with exploratory hypothesis testing using the matrix-F prior
distribution \insertCite{Mulder2018}{BGGM}. A test for each partial correlation in the model for any number
of groups. This provides evidence for the null hypothesis of no difference and the alternative hypothesis
of difference. With more than two groups, the test is for \emph{all} groups simultaneously (i.e., the relation
is the same or different in all groups). This method was introduced in \insertCite{williams2020comparing;textual}{BGGM}.
For confirmatory hypothesis testing see \code{confirm_groups}.
}
\details{
\strong{Controlling for Variables}:

When controlling for variables, it is assumed that \code{Y} includes \emph{only}
the nodes in the GGM and the control variables. Internally, \code{only} the predictors
that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
\code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
should be included in the GGM. An example is provided below.

\strong{Mixed Type}:

 The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
 continuous or \emph{only} discrete variables. This is based on the ranked likelihood which requires sampling
 the ranks for each variable (i.e., the data is not merely transformed to ranks). This is computationally
 expensive when there are many levels. For example, with continuous data, there are as many ranks
 as data points!

 The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
 and the "emprical" distribution is used otherwise. This is accomplished by specifying an indicator
 vector of length \emph{p}. A one indicates to use the ranks, whereas a zero indicates to "ignore"
 that variable. By default all integer variables are handled as ranks.

\strong{Dealing with Errors}:

An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):

\itemize{

\item The first is due to sampling the thresholds, especially when the data is heavily skewed.
      This can result in an ill-defined matrix. If this occurs, we recommend to first try
      decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
      change the data type to \code{type = mixed} which then estimates a copula GGM
      (this method can be used for data containing \strong{only} ordinal variable). This should
      work without a problem.

\item  The second is due to how the ordinal data are categorized. For example, if the error states
       that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
       the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.

}
}
\note{
\strong{"Default" Prior}:

 In Bayesian statistics, a default Bayes factor needs to have several properties. I refer
 interested users to \insertCite{@section 2.2 in @dablander2020default;textual}{BGGM}. In
 \insertCite{Williams2019_bf;textual}{BGGM}, some of these propteries were investigated, such
 model selection consistency. That said, we would not consider this a "default" Bayes factor and
 thus we encourage users to perform sensitivity analyses by varying the scale of the prior
 distribution.

 Furthermore, it is important to note there is no "correct" prior and, also, there is no need
 to entertain the possibility of a "true" model. Rather, the Bayes factor can be interpreted as
 which hypothesis best (relative to each other) predicts the observed data
 \insertCite{@Section 3.2 in @Kass1995}{BGGM}.

\strong{Interpretation of Conditional (In)dependence Models for Latent Data}:

See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
(i.e, all data types besides \code{"continuous"})
}
\examples{

\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
                   select = -c(gender,
                               education))[,1:10]

Yfemale <- subset(Y, gender == 2,
                     select = -c(gender,
                                 education))[,1:10]

##########################
### example 1: ordinal ###
##########################

# fit model
fit <- ggm_compare_explore(Ymale,  Yfemale,
                           type = "ordinal",
                           iter = 250,
                           progress = FALSE)
# summary
summ <- summary(fit)

# edge set
E <- select(fit)
}

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcor_2_cor.BGGM.R
\name{pcor_to_cor}
\alias{pcor_to_cor}
\title{Compute Correlations from the Partial Correlations}
\usage{
pcor_to_cor(object, iter = NULL)
}
\arguments{
\item{object}{An object of class \code{estimate} or \code{explore}}

\item{iter}{numeric. How many iterations (i.e., posterior samples) should be used ?
The default uses all of the samples, but note that this can take a long
time with large matrices.}
}
\value{
\itemize{

\item \code{R} An array including the correlation matrices
              (of dimensions \emph{p} by \emph{p} by \emph{iter})

\item \code{R_mean} Posterior mean of the correlations (of dimensions \emph{p} by \emph{p})
}
}
\description{
Convert the partial correlation matrices into correlation matrices. To our knowledge,
this is the only Bayesian
implementation in \code{R} that can estiamte Pearson's,  tetrachoric (binary), polychoric
(ordinal with more than two cateogries), and rank based correlation coefficients.
}
\note{
The 'default' prior distributions are specified for partial correlations in particular. This
means that the implied prior distribution will not be the same for the correlations.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- BGGM::ptsd

#########################
###### continuous #######
#########################

# estimate the model
fit <- estimate(Y, iter = 250,
                progress = FALSE)

# compute correlations
cors <- pcor_to_cor(fit)


#########################
###### ordinal  #########
#########################

# first level must be 1 !
Y <- Y + 1

# estimate the model
fit <- estimate(Y, type =  "ordinal",
                iter = 250,
                progress = FALSE)

# compute correlations
cors <- pcor_to_cor(fit)


#########################
#######   mixed    ######
#########################

# rank based correlations

# estimate the model
fit <- estimate(Y, type =  "mixed",
                iter = 250,
                progress = FALSE)

# compute correlations
cors <- pcor_to_cor(fit)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{rsa}
\alias{rsa}
\title{Data: Resilience Scale of Adults (RSA)}
\format{
A data frame with 28 variables and 1973 observations (5 point Likert scale)
}
\usage{
data("rsa")
}
\description{
A dataset containing items from the Resilience Scale of Adults (RSA). There are 33 items  and
675 observations
}
\details{
\itemize{
  \item \code{1}  My plans for the future are
  \item \code{2}  When something unforeseen happens
  \item \code{3}  My family understanding of what is important in life is
  \item \code{4}  I feel that my future looks
  \item \code{5}  My goals
  \item \code{6}  I can discuss personal issues with
  \item \code{7}  I feel
  \item \code{8}  I enjoy being
  \item \code{9}  Those who are good at encouraging are
  \item \code{10} The bonds among my friends
  \item \code{11} My personal problems
  \item \code{12} When a family member experiences a crisis/emergency
  \item \code{13} My family is characterised by
  \item \code{14} To be flexible in social settings
  \item \code{15} I get support from
  \item \code{16} In difficult periods my family
  \item \code{17} My judgements and decisions
  \item \code{18} New friendships are something
  \item \code{19} When needed, I have
  \item \code{20} I am at my best when I
  \item \code{21} Meeting new people is
  \item \code{22} When I am with others
  \item \code{23} When I start on new things/projects
  \item \code{24} Facing other people, our family acts
  \item \code{25} Belief in myself
  \item \code{26} For me, thinking of good topics of conversation is
  \item \code{27} My close friends/family members
  \item \code{28} I am good at
  \item \code{29} In my family, we like to
  \item \code{30} Rules and regular routines
  \item \code{31} In difficult periods I have a tendency to
  \item \code{32} My goals for the future are
  \item \code{33} Events in my life that I cannot influence
  \item \code{gender} "M" (male) or "F" (female)

}
}
\note{
There are 6 domains

Planned future: items 1, 4, 5, 32

Perception of self: items 2, 11, 17, 25, 31, 33

Family cohesion: items 3, 7, 13, 16, 24, 29

Social resources: items 6, 9, 10, 12, 15, 19, 27

Social Competence: items 8, 14, 18, 21, 22, 26,

Structured style: items 23, 28, 30
}
\examples{
data("rsa")

}
\references{
Briganti, G., & Linkowski, P. (2019). Item and domain network structures of the Resilience
Scale for Adults in 675 university students. Epidemiology and psychiatric sciences, 1-9.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictability.R
\name{predictability}
\alias{predictability}
\title{Predictability: Bayesian Variance Explained (R2)}
\usage{
predictability(
  object,
  select = FALSE,
  cred = 0.95,
  BF_cut = 3,
  iter = NULL,
  progress = TRUE,
  ...
)
}
\arguments{
\item{object}{object of class \code{estimate} or \code{explore}}

\item{select}{logical. Should the graph be selected ? The default is currently \code{FALSE}.}

\item{cred}{numeric. credible interval between 0 and 1  (default is 0.95) that is used for selecting the graph.}

\item{BF_cut}{numeric. evidentiary threshold (default is 3).}

\item{iter}{interger. iterations (posterior samples) used for computing R2.}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{currently ignored.}
}
\value{
An object of classes \code{bayes_R2} and \code{metric}, including

\itemize{

\item \code{scores} A list containing the posterior samples of R2. The is one element

for each node.

}
}
\description{
Compute nodewise predictability or  Bayesian variance explained \insertCite{@R2 @gelman_r2_2019}{BGGM}.
              In the context of GGMs, this method was described in \insertCite{Williams2019;textual}{BGGM}.
}
\note{
\strong{Binary and Ordinal Data}:

R2 is computed from the latent data.


\strong{Mixed Data}:

The mixed data approach is somewhat ad-hoc \insertCite{@see for example p. 277 in  @hoff2007extending;textual}{BGGM}. This
is becaue uncertainty in the ranks is not incorporated, which means that variance explained is computed from
the 'empirical' \emph{CDF}.

\strong{Model Selection}:

Currently the default to include all nodes in the model when computing R2. This can be changed (i.e., \code{select = TRUE}), which
then sets those edges not detected to zero. This is accomplished by subsetting the correlation matrix according to each neighborhood
of relations.
}
\examples{
\donttest{

# data
Y <- ptsd[,1:5]

fit <- estimate(Y, iter = 250, progress = FALSE)

r2 <- predictability(fit, select = TRUE,
                     iter = 250, progress = FALSE)

# summary
r2
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{plot.summary.estimate}
\alias{plot.summary.estimate}
\title{Plot \code{summary.estimate} Objects}
\usage{
\method{plot}{summary.estimate}(x, color = "black", size = 2, width = 0, ...)
}
\arguments{
\item{x}{An object of class \code{summary.estimate}}

\item{color}{Character string. The color for the error bars.
(defaults to \code{"black"}).}

\item{size}{Numeric. The size for the points (defaults to \code{2}).}

\item{width}{Numeric. The width of error bar ends (defaults to \code{0}).}

\item{...}{Currently ignored}
}
\value{
A \code{ggplot} object.
}
\description{
Visualize the posterior distributions for each partial correlation.
}
\examples{
\donttest{
# data
Y <- ptsd[,1:5]

fit <- estimate(Y, iter = 250,
                progress = FALSE)


plot(summary(fit))

}

}
\seealso{
\code{\link{estimate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.ggm_compare_estimate.R
\name{select.ggm_compare_estimate}
\alias{select.ggm_compare_estimate}
\title{Graph Selection for \code{ggm_compare_estimate} Objects}
\usage{
\method{select}{ggm_compare_estimate}(object, cred = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{estimate.default}.}

\item{cred}{Numeric. The credible interval width for selecting the graph
(defaults to 0.95; must be between 0 and 1).}

\item{...}{not currently used}
}
\value{
The returned object of class \code{select.ggm_compare_estimate} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:


\itemize{

\item \code{mean_diff} A list of matrices for each group comparsion (partial correlation differences).

\item \code{pcor_adj} A list of weighted adjacency matrices for each group comparsion.

\item \code{adj} A list of adjacency matrices for each group comparsion.

}
}
\description{
Provides the selected graph (of differences) based on credible intervals for
the partial correlations that did not contain zero
\insertCite{Williams2019}{BGGM}.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes
##################
### example 1: ###
##################
# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
               select = -c(gender,
                           education))

Yfemale <- subset(Y, gender == 2,
                  select = -c(gender,
                              education))

# fit model
fit <- ggm_compare_estimate(Ymale, Yfemale,
                           type = "continuous",
                           iter = 250,
                           progress = FALSE)


E <- select(fit)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fisher_r2z.R
\name{fisher_r_to_z}
\alias{fisher_r_to_z}
\title{Fisher Z Transformation}
\usage{
fisher_r_to_z(r)
}
\arguments{
\item{r}{correlation (can be a vector)}
}
\value{
Fisher Z transformed correlation(s)
}
\description{
Tranform correlations to Fisher's Z
}
\examples{
fisher_r_to_z(0.5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roll_your_own.R
\name{roll_your_own}
\alias{roll_your_own}
\title{Compute Custom Network Statistics}
\usage{
roll_your_own(
  object,
  FUN,
  iter = NULL,
  select = FALSE,
  cred = 0.95,
  progress = TRUE,
  ...
)
}
\arguments{
\item{object}{An object of class \code{estimate}.}

\item{FUN}{A custom function for computing the statistic. The first argument must be
a partial correlation matrix.}

\item{iter}{Number of iterations (posterior samples; defaults to the number in the object).}

\item{select}{Logical. Should the graph be selected ? The default is currently \code{FALSE}.}

\item{cred}{Numeric. Credible interval between 0 and 1  (default is 0.95) that is used for selecting the graph.}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{Arguments passed to the function.}
}
\value{
An object defined by \code{FUN}.
}
\description{
This function allows for computing custom network statistics for
weighted adjacency matrices (partial correlations). The statistics are computed for
each of the sampled matrices, resulting in a distribution.
}
\details{
The user has complete control of this function. Hence, care must be taken as to what \code{FUN}
returns and in what format. The function should return a single number (one for the entire GGM)
or a vector (one for each node). This ensures that the print and \code{\link{plot.roll_your_own}}
will work.


When \code{select = TRUE}, the graph is selected and then the network statistics are computed based on
the weigthed adjacency matrix. This is accomplished internally by multiplying each of the sampled
partial correlation matrices by the adjacency matrix.
}
\examples{
\donttest{
####################################
###### example 1: assortment #######
####################################
# assortment
library(assortnet)

Y <- BGGM::bfi[,1:10]
membership <- c(rep("a", 5), rep("c", 5))

# fit model
fit <- estimate(Y = Y, iter = 250,
                progress = FALSE)

# membership
membership <- c(rep("a", 5), rep("c", 5))

# define function
f <- function(x,...){
 assortment.discrete(x, ...)$r
}


net_stat <- roll_your_own(object = fit,
                          FUN = f,
                          types = membership,
                          weighted = TRUE,
                          SE = FALSE, M = 1,
                          progress = FALSE)

# print
net_stat


############################################
###### example 2: expected influence #######
############################################
# expected influence from this package
library(networktools)

# data
Y <- depression

# fit model
fit <- estimate(Y = Y, iter = 250)

# define function
f <- function(x,...){
     expectedInf(x,...)$step1
}

# compute
net_stat <- roll_your_own(object = fit,
                          FUN = f,
                          progress = FALSE)

#######################################
### example 3: mixed data & bridge ####
#######################################
# bridge from this package
library(networktools)

# data
Y <- ptsd[,1:7]

fit <- estimate(Y,
                type = "mixed",
                iter = 250)

# clusters
communities <- substring(colnames(Y), 1, 1)

# function is slow
f <- function(x, ...){
 bridge(x, ...)$`Bridge Strength`
}

net_stat <- roll_your_own(fit,
                          FUN = f,
                          select = TRUE,
                          communities = communities,
                          progress = FALSE)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_ordinal.R
\name{gen_ordinal}
\alias{gen_ordinal}
\title{Generate Ordinal and Binary data}
\usage{
gen_ordinal(n, p, levels = 2, cor_mat, empirical = FALSE)
}
\arguments{
\item{n}{Number of observations (\emph{n}).}

\item{p}{Number of variables  (\emph{p}).}

\item{levels}{Number of categories (defaults to 2; binary data).}

\item{cor_mat}{A \emph{p} by \emph{p} matrix including the true correlation structure.}

\item{empirical}{Logical. If true, \code{cor_mat} specifies  the empirical not
population covariance matrix.}
}
\value{
A \emph{n} by \emph{p} data matrix.
}
\description{
Generate Multivariate Ordinal and Binary data.
}
\note{
In order to allow users to enjoy the functionality of \bold{BGGM}, we had to make minor changes to the function \code{rmvord_naiv}
from the \code{R} package \bold{orddata} \insertCite{orddata}{BGGM}. All rights to, and credit for, the function \code{rmvord_naiv}
belong to the authors of that package.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
A copy of the GNU General Public License is available online.
}
\examples{
################################
######### example 1 ############
################################

main <-  ptsd_cor1[1:5,1:5]
p <- ncol(main)

pcors <- -(cov2cor(solve(main)) -diag(p))
diag(pcors) <- 1
pcors <- ifelse(abs(pcors) < 0.05, 0, pcors)

inv <-  -pcors
diag(inv) <- 1
cors <- cov2cor( solve(inv))

# example data
Y <- BGGM::gen_ordinal(n = 500, p = 5,
                       levels = 2,
                       cor_mat = cors,
                       empirical = FALSE)



################################
######### example 2 ############
################################
# empirical = TRUE

Y <-  gen_ordinal(n = 500,
                  p = 16,
                  levels = 5,
                  cor_mat = ptsd_cor1,
                  empirical = TRUE)

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_sum.R
\name{plot.pcor_sum}
\alias{plot.pcor_sum}
\title{Plot \code{pcor_sum} Object}
\usage{
\method{plot}{pcor_sum}(x, fill = "#CC79A7", ...)
}
\arguments{
\item{x}{An object of class \code{posterior_sum}}

\item{fill}{Character string. What fill for the histogram
(defaults to colorblind "pink")?}

\item{...}{Currently ignored.}
}
\value{
A list of \code{ggplot} objects
}
\description{
Plot \code{pcor_sum} Object
}
\note{
\strong{Examples}:
}
\seealso{
pcor_sum
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{tas}
\alias{tas}
\title{Data: Toronto Alexithymia Scale (TAS)}
\format{
A data frame with 20 variables and 1925 observations (5 point Likert scale)
}
\usage{
data("tas")
}
\description{
A dataset containing items from the Toronto Alexithymia Scale (TAS). There are 20 variables  and
1925 observations
}
\details{
\itemize{
  \item \code{1} I am often confused about what emotion I am feeling
  \item \code{2}  It is difficult for me to find the right words for my feelings
  \item \code{3} I have physical sensations that even doctors don’t understand
  \item \code{4} I am able to describe my feelings easily
  \item \code{5} I prefer to analyze problems rather than just describe them
  \item \code{6} When I am upset, I don’t know if I am sad, frightened, or angry
  \item \code{7}  I am often puzzled by sensations in my body
  \item \code{8} I prefer just to let things happen rather than to understand why they turned out that way
  \item \code{9}  I have feelings that I can’t quite identify
  \item \code{10} Being in touch with emotions is essential
  \item \code{11}  I find it hard to describe how I feel about people
  \item \code{12} People tell me to describe my feelings more
  \item \code{13} I don’t know what’s going on inside me
  \item \code{14} I often don’t know why I am angry
  \item \code{15} I prefer talking to people about their daily activities rather than their feelings
  \item \code{16}  I prefer to watch “light” entertainment shows rather than psychological dramas
  \item \code{17} It is difficult for me to reveal my innermost feelings, even to close friends
  \item \code{18}  I can feel close to someone, even in moments of silence
  \item \code{19}  I find examination of my feelings useful in solving personal problems
  \item \code{20} Looking for hidden meanings in movies or plays distracts from their enjoyment
  \item \code{gender} "M" (male) or "F" (female)

}
}
\note{
There are three domains

Difficulty identifying feelings: items 1, 3, 6, 7, 9, 13, 14

Difficulty describing feelings: items 2, 4, 11, 12, 17

Externally oriented thinking: items 10, 15, 16, 18, 19
}
\examples{
data("tas")

}
\references{
Briganti, G., & Linkowski, P. (2019). Network approach to items and domains from
the Toronto Alexithymia Scale. Psychological reports.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{bfi}
\alias{bfi}
\title{Data: 25 Personality items representing 5 factors}
\format{
A data frame with 25 variables and 2800 observations (including missing values)
}
\usage{
data("bfi")
}
\description{
This dataset and the corresponding documentation was taken from the \strong{psych} package. We refer users to that
package for further details \insertCite{psych}{BGGM}.
}
\details{
\itemize{
  \item \code{A1} Am indifferent to the feelings of others. (q_146)
  \item \code{A2} Inquire about others' well-being. (q_1162)
  \item \code{A3} Know how to comfort others. (q_1206)
  \item \code{A4} Love children. (q_1364)
  \item \code{A5} Make people feel at ease. (q_1419)
  \item \code{C1} Am exacting in my work. (q_124)
  \item \code{C2} Continue until everything is perfect. (q_530)
  \item \code{C3} Do things according to a plan. (q_619)
  \item \code{C4} Do things in a half-way manner. (q_626)
  \item \code{C5} Waste my time. (q_1949)
  \item \code{E1} Don't talk a lot. (q_712)
  \item \code{E2} Find it difficult to approach others. (q_901)
  \item \code{E3} Know how to captivate people. (q_1205)
  \item \code{E4} Make friends easily. (q_1410)
  \item \code{E5} Take charge. (q_1768)
  \item \code{N1} Get angry easily. (q_952)
  \item \code{N2} Get irritated easily. (q_974)
  \item \code{N3} Have frequent mood swings. (q_1099)
  \item \code{N4} Often feel blue. (q_1479)
  \item \code{N5} Panic easily. (q_1505)
  \item \code{o1} Am full of ideas. (q_128)
  \item \code{o2} Avoid difficult reading material.(q_316)
  \item \code{o3} Carry the conversation to a higher level. (q_492)
  \item \code{o4} Spend time reflecting on things. (q_1738)
  \item \code{o5} Will not probe deeply into a subject. (q_1964)
  \item \code{gender} Males = 1, Females =2
  \item \code{education} 1 = HS, 2 = finished HS, 3 = some college, 4 = college graduate 5 = graduate degree
}
}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvn_imputation.R
\name{impute_data}
\alias{impute_data}
\title{Obtain Imputed Datasets}
\usage{
impute_data(
  Y,
  type = "continuous",
  lambda = NULL,
  mixed_type = NULL,
  iter = 1000,
  progress = TRUE
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{type}{Character string. Which type of data for \code{Y} ? The options include \code{continuous},
\code{binary}, \code{ordinal}, or \code{mixed}. Note that mixed can be used for data with only
ordinal variables. See the note for further details.}

\item{lambda}{Numeric. A regularization parameter, which defaults to p + 2. A larger value results
in more shrinkage.}

\item{mixed_type}{Numeric vector. An indicator of length \emph{p} for which variables should be treated as ranks.
(1 for rank and 0 to assume the observed marginal distribution).
The default is currently to treat all integer variables as ranks when
\code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{iter}{Number of iterations (posterior samples; defaults to 1000).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}
}
\value{
An object of class \code{mvn_imputation}:

\itemize{

\item \code{imputed_datasets} An array including the imputed datasets.

}
}
\description{
Impute missing values, assuming a  multivariate normal distribution, with the posterior
predictive distribution. For binary, ordinal, and mixed (a combination of discrete and continuous)
data, the values are first imputed for the latent data and then converted to the original scale.
}
\details{
Missing values are imputed with the approach described in \insertCite{hoff2009first;textual}{BGGM}.
The basic idea is to impute the missing values with the respective posterior pedictive distribution,
given the observed data, as the model is being estimated. Note that the default is \code{TRUE},
but this ignored when there are no missing values. If set to \code{FALSE}, and there are missing
values, list-wise deletion is performed with \code{na.omit}.
}
\examples{
\donttest{
# obs
n <- 5000

# n missing
n_missing <- 1000

# variables
p <- 16

# data
Y <- MASS::mvrnorm(n, rep(0, p), ptsd_cor1)

# for checking
Ymain <- Y

# all possible indices
indices <- which(matrix(0, n, p) == 0,
                 arr.ind = TRUE)

# random sample of 1000 missing values
na_indices <- indices[sample(5:nrow(indices),
                             size = n_missing,
                             replace = FALSE),]

# fill with NA
Y[na_indices] <- NA

# missing = 1
Y_miss <- ifelse(is.na(Y), 1, 0)

# true values (to check)
true <- unlist(sapply(1:p, function(x)
        Ymain[which(Y_miss[,x] == 1),x] ))

# impute
fit_missing <- impute_data(Y, progress = FALSE, iter = 250)

# impute
fit_missing <- impute_data(Y,
                           progress = TRUE,
                           iter = 250)

}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_confirm.R
\name{ggm_compare_confirm}
\alias{ggm_compare_confirm}
\title{GGM Compare: Confirmatory Hypothesis Testing}
\usage{
ggm_compare_confirm(
  ...,
  hypothesis,
  formula = NULL,
  type = "continuous",
  mixed_type = NULL,
  prior_sd = 0.25,
  iter = 25000,
  impute = TRUE,
  progress = TRUE,
  seed = 1
)
}
\arguments{
\item{...}{At least two matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (nodes).}

\item{hypothesis}{Character string. The hypothesis (or hypotheses) to be tested. See notes for futher details.}

\item{formula}{an object of class \code{\link[stats]{formula}}. This allows for including
control variables in the model (i.e., \code{~ gender}).}

\item{type}{Character string. Which type of data for \code{Y} ? The options include \code{continuous},
\code{binary}, \code{ordinal}, or \code{mixed}. Note that mixed can be used for data with only
ordinal variables. See the note for further details.}

\item{mixed_type}{numeric vector. An indicator of length p for which varibles should be treated as ranks.
(1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{prior_sd}{Numeric. The scale of the prior distribution (centered at zero),
in reference to a beta distribtuion (defaults to 0.25).}

\item{iter}{Number of iterations (posterior samples; defaults to 25,000).}

\item{impute}{Logicial. Should the missing values (\code{NA})
be imputed during model fitting (defaults to \code{TRUE}) ?}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed.}
}
\value{
The returned object of class \code{confirm} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:

\itemize{

\item \code{out_hyp_prob} Posterior hypothesis probabilities.

\item \code{info} An object of class \code{BF} from the R package \strong{BFpack}
                  \insertCite{mulder2019bfpack}{BGGM}

}
}
\description{
Confirmatory hypothesis testing for comparing GGMs. Hypotheses are expressed as equality
and/or ineqaulity contraints on the partial correlations of interest. Here the focus is \emph{not}
on determining the graph (see \code{\link{explore}}) but testing specific hypotheses related to
the conditional (in)dependence structure. These methods were introduced in
\insertCite{Williams2019_bf;textual}{BGGM} and in \insertCite{williams2020comparing;textual}{BGGM}
}
\details{
The hypotheses can be written either with the respective column names or numbers.
For example, \code{g1_1--2} denotes the relation between the variables in column 1 and 2 for group 1.
The \code{g1_} is required and the only difference from \code{\link{confirm}} (one group).
Note that these must correspond to the upper triangular elements of the correlation
matrix. This is accomplished by ensuring that the first number is smaller than the second number.
This also applies when using column names (i.e,, in reference to the column number).


\strong{One Hypothesis}:

To test whether a relation in larger in one group, while both are expected
to be positive,  this can be written as

\itemize{

\item  \code{hyp <-  c(g1_1--2 > g2_1--2 > 0)}
}

This is then compared to the complement.

\strong{More Than One Hypothesis}:

The above hypothesis can also be compared to, say, a null model by using ";"
      to seperate the hypotheses, for example,

\itemize{

\item  \code{hyp <-  c(g1_1--2 > g2_1--2 > 0; g1_1--2 = g2_1--2 = 0)}.

}

Any number of hypotheses can be compared this way.

\strong{Using "&"}

 It is also possible to include \code{&}. This allows for testing one constraint \bold{and}
 another contraint as one hypothesis.

\itemize{

\item \code{hyp <- c("g1_A1--A2 > g2_A1--A2 & g1_A1--A3 = g2_A1--A3")}

}

Of course, it is then possible to include additional hypotheses by separating them with ";".

\strong{Testing Sums}

It might also be interesting to test the sum of partial correlations. For example, that the
sum of specific relations in one group is larger than the sum in another group.

\itemize{

\item \code{hyp <- c("g1_A1--A2 + g1_A1--A3 > g2_A1--A2 + g2_A1--A3;
                      g1_A1--A2 + g1_A1--A3 = g2_A1--A2 + g2_A1--A3")}

}


\strong{Potential Delays}:

There is a chance for a potentially long delay from the time the progress bar finishes
to when the function is done running. This occurs when the hypotheses require further
sampling to be tested, for example, when grouping relations
\code{c("(g1_A1--A2, g2_A2--A3) > (g2_A1--A2, g2_A2--A3)"}.
This is not an error.


\strong{Controlling for Variables}:

When controlling for variables, it is assumed that \code{Y} includes \emph{only}
the nodes in the GGM and the control variables. Internally, \code{only} the predictors
that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
\code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
should be included in the GGM. An example is provided below.

\strong{Mixed Type}:

 The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
 continuous or \emph{only} discrete variables \insertCite{hoff2007extending}{BGGM}. This is based on the
 ranked likelihood which requires sampling the ranks for each variable (i.e., the data is not merely
 transformed to ranks). This is computationally expensive when there are many levels. For example,
 with continuous data, there are as many ranks as data points!

 The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
 and the "emprical" distribution is used otherwise. This is accomplished by specifying an indicator
 vector of length \emph{p}. A one indicates to use the ranks, whereas a zero indicates to "ignore"
 that variable. By default all integer variables are handled as ranks.

\strong{Dealing with Errors}:

An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):

\itemize{

\item The first is due to sampling the thresholds, especially when the data is heavily skewed.
      This can result in an ill-defined matrix. If this occurs, we recommend to first try
      decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
      change the data type to \code{type = mixed} which then estimates a copula GGM
      (this method can be used for data containing \strong{only} ordinal variable). This should
      work without a problem.

\item  The second is due to how the ordinal data are categorized. For example, if the error states
       that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
       the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.

}


\strong{Imputing Missing Values}:

Missing values are imputed with the approach described in \insertCite{hoff2009first;textual}{BGGM}.
The basic idea is to impute the missing values with the respective posterior pedictive distribution,
given the observed data, as the model is being estimated. Note that the default is \code{TRUE},
but this ignored when there are no missing values. If set to \code{FALSE}, and there are missing
values, list-wise deletion is performed with \code{na.omit}.
}
\note{
\strong{"Default" Prior}:

 In Bayesian statistics, a default Bayes factor needs to have several properties. I refer
 interested users to \insertCite{@section 2.2 in @dablander2020default;textual}{BGGM}. In
 \insertCite{Williams2019_bf;textual}{BGGM}, some of these propteries were investigated (e.g.,
 model selection consistency). That said, we would not consider this a "default" or "automatic"
 Bayes factor and thus we encourage users to perform sensitivity analyses by varying the scale of
 the prior distribution (\code{prior_sd}).

 Furthermore, it is important to note there is no "correct" prior and, also, there is no need
 to entertain the possibility of a "true" model. Rather, the Bayes factor can be interpreted as
 which hypothesis best (relative to each other) predicts the observed data
 \insertCite{@Section 3.2 in @Kass1995}{BGGM}.

\strong{Interpretation of Conditional (In)dependence Models for Latent Data}:

 See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
(i.e, all data types besides \code{"continuous"})
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi

###############################
#### example 1: continuous ####
###############################

# males
Ymale   <- subset(Y, gender == 1,
                  select = -c(education,
                              gender))[,1:5]


# females
Yfemale <- subset(Y, gender == 2,
                     select = -c(education,
                                 gender))[,1:5]

 # exhaustive
 hypothesis <- c("g1_A1--A2 >  g2_A1--A2;
                  g1_A1--A2 <  g2_A1--A2;
                  g1_A1--A2 =  g2_A1--A2")

# test hyp
test <- ggm_compare_confirm(Ymale,  Yfemale,
                            hypothesis = hypothesis,
                            iter = 250,
                            progress = FALSE)

# print (evidence not strong)
test

#########################################
#### example 2: sensitivity to prior ####
#########################################
# continued from example 1

# decrease prior SD
test <- ggm_compare_confirm(Ymale,
                            Yfemale,
                            prior_sd = 0.1,
                            hypothesis = hypothesis,
                            iter = 250,
                            progress = FALSE)

# print
test

# indecrease prior SD
test <- ggm_compare_confirm(Ymale,
                            Yfemale,
                            prior_sd = 0.5,
                            hypothesis = hypothesis,
                            iter = 250,
                            progress = FALSE)

# print
test

################################
#### example 3: mixed data #####
################################

hypothesis <- c("g1_A1--A2 >  g2_A1--A2;
                 g1_A1--A2 <  g2_A1--A2;
                 g1_A1--A2 =  g2_A1--A2")

# test (1000 for example)
test <- ggm_compare_confirm(Ymale,
                            Yfemale,
                            type = "mixed",
                            hypothesis = hypothesis,
                            iter = 250,
                            progress = FALSE)

# print
test

##############################
##### example 4: control #####
##############################
# control for education

# data
Y <- bfi

# males
Ymale   <- subset(Y, gender == 1,
                  select = -c(gender))[,c(1:5, 26)]

# females
Yfemale <- subset(Y, gender == 2,
                  select = -c(gender))[,c(1:5, 26)]

# test
test <- ggm_compare_confirm(Ymale,
                             Yfemale,
                             formula = ~ education,
                             hypothesis = hypothesis,
                             iter = 250,
                             progress = FALSE)
# print
test


#####################################
##### example 5: many relations #####
#####################################

# data
Y <- bfi

hypothesis <- c("g1_A1--A2 > g2_A1--A2 & g1_A1--A3 = g2_A1--A3;
                 g1_A1--A2 = g2_A1--A2 & g1_A1--A3 = g2_A1--A3;
                 g1_A1--A2 = g2_A1--A2 = g1_A1--A3 = g2_A1--A3")

Ymale   <- subset(Y, gender == 1,
                  select = -c(education,
                              gender))[,1:5]


# females
Yfemale <- subset(Y, gender == 2,
                     select = -c(education,
                                 gender))[,1:5]

test <- ggm_compare_confirm(Ymale,
                            Yfemale,
                             hypothesis = hypothesis,
                             iter = 250,
                             progress = FALSE)

# print
test
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_BGGM.R
\name{print.BGGM}
\alias{print.BGGM}
\title{Print method for \code{BGGM} objects}
\usage{
\method{print}{BGGM}(x, ...)
}
\arguments{
\item{x}{An object of class \code{BGGM}}

\item{...}{currently ignored}
}
\description{
Mainly used to avoid a plethora of different print
functions that overcrowded the documentation in previous versions
of \strong{BGGM}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd}
\alias{ptsd}
\title{Data: Post-Traumatic Stress Disorder}
\format{
A dataframe with 221 rows and 20 variables
}
\usage{
data("ptsd")
}
\description{
A dataset containing items that measure Post-traumatic stress disorder symptoms \insertCite{armour2017network}{BGGM}.
There are 20 variables (\emph{p}) and  221 observations (\emph{n}).
}
\details{
\itemize{
  \item Intrusive Thoughts
  \item Nightmares
  \item Flashbacks
  \item Emotional cue reactivity
  \item Psychological cue reactivity
  \item Avoidance of thoughts
  \item Avoidance of reminders
  \item Trauma-related amnesia
  \item Negative beliefs
  \item Negative trauma-related emotions
  \item Loss of interest
  \item Detachment
  \item Restricted affect
  \item Irritability/anger
  \item Self-destructive/reckless behavior
  \item Hypervigilance
  \item Exaggerated startle response
  \item Difficulty concentrating
  \item Sleep disturbance
}
}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted_adj_mat.R
\name{weighted_adj_mat}
\alias{weighted_adj_mat}
\title{Extract the Weighted Adjacency Matrix}
\usage{
weighted_adj_mat(object, ...)
}
\arguments{
\item{object}{A model estimated with \strong{BGGM}. All classes are supported, assuming
there is matrix to be extracted.}

\item{...}{Currently ignored.}
}
\value{
The weighted adjacency matrix (partial correlation matrix with zeros).
}
\description{
Extract the weighted adjacency matrix (posterior mean) from
\code{\link{estimate}}, \code{\link{explore}}, \code{\link{ggm_compare_estimate}},
and \code{\link{ggm_compare_explore}} objects.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes
Y <- bfi[,1:5]

# estimate
fit <- estimate(Y, iter = 250,
                progress = FALSE)

# select graph
E <- select(fit)

# extract weighted adj matrix
weighted_adj_mat(E)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictability.R
\name{summary.predictability}
\alias{summary.predictability}
\title{Summary Method for \code{predictability} Objects}
\usage{
\method{summary}{predictability}(object, cred = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{predictability}.}

\item{cred}{Numeric. The credible interval width for summarizing the posterior
distributions (defaults to 0.95; must be between 0 and 1).}

\item{...}{Currently ignored}
}
\description{
Summary Method for \code{predictability} Objects
}
\examples{
\donttest{
Y <- ptsd[,1:5]

fit <- explore(Y, iter = 250,
               progress = FALSE)

r2 <- predictability(fit, iter = 250,
                     progress = FALSE)

summary(r2)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roll_your_own.R
\name{plot.roll_your_own}
\alias{plot.roll_your_own}
\title{Plot \code{roll_your_own} Objects}
\usage{
\method{plot}{roll_your_own}(x, fill = "#CC79A7", alpha = 0.5, ...)
}
\arguments{
\item{x}{An object of class \code{roll_your_own}}

\item{fill}{Character string specifying the color for the ridges.}

\item{alpha}{Numeric. Transparancey of the ridges}

\item{...}{Currently ignored}
}
\value{
An object of class \code{ggplot}
}
\description{
Plot \code{roll_your_own} Objects
}
\examples{
\donttest{
####################################
###### example 1: assortment #######
####################################
# assortment
library(assortnet)

Y <- BGGM::bfi[,1:10]
membership <- c(rep("a", 5), rep("c", 5))

# fit model
fit <- estimate(Y = Y, iter = 250,
                progress = FALSE)

# membership
membership <- c(rep("a", 5), rep("c", 5))

# define function
f <- function(x,...){
 assortment.discrete(x, ...)$r
}

net_stat <- roll_your_own(object = fit,
                          FUN = f,
                          types = membership,
                          weighted = TRUE,
                          SE = FALSE, M = 1,
                          progress = FALSE)

# plot
plot(net_stat)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.VAR_estimate.R
\name{select.var_estimate}
\alias{select.var_estimate}
\title{Graph Selection for \code{var.estimate} Object}
\usage{
\method{select}{var_estimate}(object, cred = 0.95, alternative = "two.sided", ...)
}
\arguments{
\item{object}{An object of class \code{VAR.estimate}.}

\item{cred}{Numeric. The credible interval width for selecting the graph
(defaults to 0.95; must be between 0 and 1).}

\item{alternative}{A character string specifying the alternative hypothesis. It
must be one of "two.sided" (default), "greater"  or "less".
See note for futher details.}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{select.var_estimate}, including

\itemize{

\item {pcor_adj} Adjacency matrix for the partial correlations.

\item {beta_adj} Adjacency matrix for the regression coefficients.

\item {pcor_weighted_adj} Weighted adjacency matrix for the partial correlations.

\item {beta_weighted_adj} Weighted adjacency matrix for the regression coefficients.

\item \code{pcor_mu} Partial correlation matrix (posterior mean).

\item \code{beta_mu} A matrix including the regression coefficients (posterior mean).

}
}
\description{
Graph Selection for \code{var.estimate} Object
}
\examples{
\donttest{
# data
Y <- subset(ifit, id == 1)[,-1]

# fit model with alias (var_estimate also works)
fit <- var_estimate(Y, progress = FALSE)

# select graphs
select(fit, cred = 0.95)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore.default.R
\name{plot.summary.explore}
\alias{plot.summary.explore}
\title{Plot \code{summary.explore} Objects}
\usage{
\method{plot}{summary.explore}(x, color = "black", size = 2, width = 0, ...)
}
\arguments{
\item{x}{An object of class \code{summary.explore}}

\item{color}{Character string. The color for the error bars.
(defaults to \code{"black"}).}

\item{size}{Numeric. The size for the points (defaults to \code{2}).}

\item{width}{Numeric. The width of error bar ends (defaults to \code{0} ).}

\item{...}{Currently ignored}
}
\value{
A \code{ggplot} object
}
\description{
Visualize the posterior distributions for each partial correlation.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

Y <- ptsd[,1:5]

fit <- explore(Y, iter = 250,
               progress = FALSE)

plt <- plot(summary(fit))

plt
}
}
\seealso{
\code{\link{explore}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_estimate.R
\name{plot.summary.var_estimate}
\alias{plot.summary.var_estimate}
\title{Plot \code{summary.var_estimate} Objects}
\usage{
\method{plot}{summary.var_estimate}(x, color = "black", size = 2, width = 0, param = "all", order = TRUE, ...)
}
\arguments{
\item{x}{An object of class \code{summary.var_estimate}}

\item{color}{Character string. The color for the error bars.
(defaults to \code{"black"}).}

\item{size}{Numeric. The size for the points (defaults to \code{2}).}

\item{width}{Numeric. The width of error bar ends (defaults to \code{0}).}

\item{param}{Character string. Which parameters should be plotted ? The options
are \code{pcor}, \code{beta}, or \code{all} (default).}

\item{order}{Logical. Should the relations be ordered by size (defaults to \code{TRUE}) ?}

\item{...}{Currently ignored}
}
\value{
A list of \code{ggplot} objects.
}
\description{
Visualize the posterior distributions of each partial correlation and
regression coefficient.
}
\examples{
\donttest{

# data
Y <- subset(ifit, id == 1)[,-1]

# fit model with alias (var_estimate also works)
fit <- var_estimate(Y, progress = FALSE)

plts <- plot(summary(fit))
plts$pcor_plt
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{asd_ocd}
\alias{asd_ocd}
\title{Data: Autism and Obssesive Compulsive Disorder}
\format{
A correlation matrix including 17 variables. These data were measured on a 4 level likert scale.
}
\usage{
data("asd_ocd")
}
\description{
A correlation matrix with 17 variables in total (autsim: 9; OCD: 8).
The sample size was 213.
}
\details{
\strong{Autism}:

\itemize{

  \item \code{CI}  Circumscribed interests
  \item \code{UP}  Unusual preoccupations
  \item \code{RO}  Repetitive use of objects or interests in parts of objects
  \item \code{CR}  Compulsions and/or rituals
  \item \code{CI}  Unusual sensory interests
  \item \code{SM}  Complex mannerisms or stereotyped body movements
  \item \code{SU}  Stereotyped utterances/delayed echolalia
  \item \code{NIL} Neologisms and/or idiosyncratic language
  \item \code{VR}  Verbal rituals
}

\strong{OCD}

\itemize{
  \item \code{CD} Concern with things touched due to dirt/bacteria
  \item \code{TB} Thoughts of doing something bad around others
  \item \code{CT} Continual thoughts that do not go away
  \item \code{HP} Belief that someone/higher power put reoccurring thoughts in their head
  \item \code{CW} Continual washing
  \item \code{CCh} Continual checking CntCheck
  \item \code{CC} Continual counting/repeating
  \item \code{RD} Repeatedly do things until it feels good or just right

}
}
\examples{
data("asd_ocd")

# generate continuous
Y <- MASS::mvrnorm(n = 213,
                   mu = rep(0, 17),
                   Sigma = asd_ocd,
                   empirical = TRUE)


}
\references{
Jones, P. J., Ma, R., & McNally, R. J. (2019). Bridge centrality:
A network approach
to understanding comorbidity. Multivariate behavioral research, 1-15.

Ruzzano, L., Borsboom, D., & Geurts, H. M. (2015).
Repetitive behaviors in autism and obsessive-compulsive
disorder: New perspectives from a network analysis.
Journal of Autism and Developmental Disorders, 45(1),
192-202. doi:10.1007/s10803-014-2204-9
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_estimate.R
\name{var_estimate}
\alias{var_estimate}
\title{VAR: Estimation}
\usage{
var_estimate(
  Y,
  rho_sd = 0.5,
  beta_sd = 1,
  iter = 5000,
  progress = TRUE,
  seed = 1,
  ...
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{rho_sd}{Numeric. Scale of the prior distribution for the partial correlations,
approximately the standard deviation of a beta distribution
(defaults to 0.50).}

\item{beta_sd}{Numeric. Standard deviation of the prior distribution for the regression coefficients
(defaults to 1). The prior is by default centered at zero and follows a normal distribution
\insertCite{@Equation 9, @sinay2014bayesian}{BGGM}}

\item{iter}{Number of iterations (posterior samples; defaults to 5000).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed (defaults to 1).}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{var_estimate} containing a lot of information that is
used for printing and plotting the results. For users of \strong{BGGM}, the following are the
useful objects:

\itemize{

\item \code{beta_mu} A matrix including the regression coefficients (posterior mean).

\item \code{pcor_mu} Partial correlation matrix (posterior mean).

\item \code{fit} A list including the posterior samples.

}
}
\description{
Estimate VAR(1) models by efficiently sampling from the posterior distribution. This
provides two graphical structures: (1) a network of undirected relations (the GGM, controlling for the
lagged predictors) and (2) a network of directed relations (the lagged coefficients). Note that
in the graphical modeling literature, this model is also known as a time series chain graphical model
\insertCite{abegaz2013sparse}{BGGM}.
}
\details{
Each time series in \code{Y} is standardized (mean  = 0; standard deviation = 1).
}
\note{
\strong{Regularization}:

A Bayesian ridge regression can be fitted by decreasing \code{beta_sd}
(e.g., \code{beta_sd = 0.25}). This could be advantageous for forecasting
(out-of-sample prediction) in particular.
}
\examples{
\donttest{
# data
Y <- subset(ifit, id == 1)[,-1]

# use alias (var_estimate also works)
fit <- var_estimate(Y, progress = FALSE)

fit

}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constrained_post.R
\name{constrained_posterior}
\alias{constrained_posterior}
\title{Constrained Posterior Distribution}
\usage{
constrained_posterior(
  object,
  adj,
  method = "direct",
  iter = 5000,
  progress = TRUE,
  ...
)
}
\arguments{
\item{object}{An object of class \code{estimate} or \code{explore}}

\item{adj}{A \code{p} by \code{p} adjacency matrix. The zero entries denote the
elements that should be constrained to zero.}

\item{method}{Character string. Which method should be used ? Defaults to
the "direct sampler" (i.e., \code{method = "direct"}) described in
\insertCite{@page 122, section 2.4,  @lenkoski2013direct;textual}{BGGM}. The other
option is a Metropolis-Hastings algorithm (\code{MH}).
See details.}

\item{iter}{Number of iterations (posterior samples; defaults to 5000).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{contrained}, including

\itemize{

\item \code{precision_mean} The posterior mean for the precision matrix.

\item \code{pcor_mean} The posterior mean for the precision matrix.

\item \code{precision_samps} A 3d array of dimension \code{p} by \code{p} by \code{iter}
                             including the sampled precision matrices.

 \item \code{pcor_samps} A 3d array of dimension \code{p} by \code{p} by \code{iter}
                             including sampled partial correlations matrices.
}
}
\description{
Compute the posterior distribution
             with off-diagonal elements of the precision matrix constrained
             to zero.
}
\examples{
\donttest{

# data
Y <- bfi[,1:10]

# sample posterior
fit <- estimate(Y, iter = 100)

# select graph
sel <- select(fit)

# constrained posterior
post <- constrained_posterior(object = fit,
                              adj = sel$adj,
                              iter = 100,
                              progress = FALSE)

}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.estimate.R
\name{predict.explore}
\alias{predict.explore}
\title{Model Predictions for \code{explore} Objects}
\usage{
\method{predict}{explore}(
  object,
  newdata = NULL,
  summary = TRUE,
  cred = 0.95,
  iter = NULL,
  progress = TRUE,
  ...
)
}
\arguments{
\item{object}{object of class \code{explore}}

\item{newdata}{an optional data frame for obtaining predictions (e.g., on test data)}

\item{summary}{summarize the posterior samples (defaults to \code{TRUE}).}

\item{cred}{credible interval used for summarizing}

\item{iter}{number of posterior samples (defaults to all in the object).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{currently ignored}
}
\value{
\code{summary = TRUE}: 3D array of dimensions n (observations),
        4 (posterior summary),
        p (number of nodes). \code{summary = FALSE}:
        list containing predictions for each variable
}
\description{
Model Predictions for \code{explore} Objects
}
\examples{
\donttest{
# data
Y <- ptsd

# fit model
fit <- explore(Y, iter = 250,
               progress = FALSE)

# predict
pred <- predict(fit,
                progress = FALSE)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_prior.R
\name{plot_prior}
\alias{plot_prior}
\title{Plot: Prior Distribution}
\usage{
plot_prior(prior_sd = 0.2, iter = 5000)
}
\arguments{
\item{prior_sd}{Scale of the prior distribution, approximately the standard deviation
of a beta distribution (defaults to 0.25).}

\item{iter}{Number of iterations (prior samples; defaults to 5000).}
}
\value{
A \code{ggplot} object.
}
\description{
Visualize the implied prior distribution for the partial correlations. This is
             particularly useful for the Bayesian hypothesis testing methods.
}
\examples{
# note: iter = 250 for demonstrative purposes

plot_prior(prior_sd = 0.25, iter = 250)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_estimate.R
\name{summary.var_estimate}
\alias{summary.var_estimate}
\title{Summary Method for \code{var_estimate} Objects}
\usage{
\method{summary}{var_estimate}(object, cred = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{var_estimate}}

\item{cred}{Numeric. The credible interval width for summarizing the posterior
distributions (defaults to 0.95; must be between 0 and 1).}

\item{...}{Currently ignored.}
}
\value{
A dataframe containing the summarized posterior distributions,
including both the partial correlations and the regression coefficients.

\itemize{

\item \code{pcor_results} A data frame including the summarized partial correlations

\item \code{beta_results} A list containing the summarized regression coefficients (one
data frame for each outcome)
}
}
\description{
Summarize the posterior distribution of each partial correlation
and regression coefficient with the posterior mean, standard deviation, and
credible intervals.
}
\examples{
\donttest{
# data
Y <- subset(ifit, id == 1)[,-1]

# fit model with alias (var_estimate also works)
fit <- var_estimate(Y, progress = FALSE)

# summary ('pcor')
print(
summary(fit, cred = 0.95),
param = "pcor",
)


# summary ('beta')
print(
summary(fit, cred = 0.95),
param = "beta",
)

}
}
\seealso{
\code{\link{var_estimate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{iri}
\alias{iri}
\title{Data: Interpersonal Reactivity Index (IRI)}
\format{
A data frame with 28 variables and 1973 observations (5 point Likert scale)
}
\usage{
data("iri")
}
\description{
A dataset containing items from the Interpersonal Reactivity Index (IRI; an empathy measure). There are 28 variables  and
1973 observations
}
\details{
\itemize{
  \item \code{1} I daydream and fantasize, with some regularity, about things that might happen to me.
  \item \code{2}  I often have tender, concerned feelings for people less fortunate than me.
  \item \code{3} I sometimes find it difficult to see things from the "other guy's" point of view.
  \item \code{4} Sometimes I don't feel very sorry for other people when they are having problems.
  \item \code{5}  I really get involved with the feelings of the characters in a novel.
  \item \code{6} In emergency situations, I feel apprehensive and ill-at-ease.
  \item \code{7}  I am usually objective when I watch a movie or play, and I don't often get completely caught up in it.
  \item \code{8} I try to look at everybody's side of a disagreement before I make a decision.
  \item \code{9}  When I see someone being taken advantage of, I feel kind of protective towards them.
  \item \code{10} I sometimes feel helpless when I am in the middle of a very emotional situation.
  \item \code{11}  I sometimes try to understand my friends better
  by imagining how things look from their perspective
  \item \code{12} Becoming extremely involved in a good book or movie is somewhat rare for me.
  \item \code{13} When I see someone get hurt, I tend to remain calm.
  \item \code{14} Other people's misfortunes do not usually disturb me a great deal.
  \item \code{15}  If I'm sure I'm right about something, I don't waste much
  time listening to other people's arguments.
  \item \code{16}  After seeing a play or movie, I have felt as though I were one of the characters.
  \item \code{17}  Being in a tense emotional situation scares me.
  \item \code{18}  When I see someone being treated unfairly,
  I sometimes don't feel very much pity for them.
  \item \code{19} I am usually pretty effective in dealing with emergencies.
  \item \code{20} I am often quite touched by things that I see happen.
  \item \code{21} I believe that there are two sides to every question and try to look at them both.
  \item \code{22} I would describe myself as a pretty soft-hearted person.
  \item \code{23} When I watch a good movie, I can very easily put myself in
  the place of a leading character
  \item \code{24}  I tend to lose control during emergencies.
  \item \code{25} When I'm upset at someone, I usually try to "put myself in his shoes" for a while.
  \item \code{26} When I am reading an interesting story or novel, I imagine how I would feel if the
  events in the story were happening to me.
  \item \code{27}  When I see someone who badly needs help in an emergency, I go to pieces.
  \item \code{28} Before criticizing somebody, I try to imagine how I would feel if I were in their place.
  \item \code{gender} "M" (male) or "F" (female)

}
}
\note{
There are four domains

Fantasy: items 1, 5, 7, 12, 16, 23, 26

Perspective taking: items 3, 8, 11, 15, 21, 25, 28

Empathic concern: items 2, 4, 9, 14, 18, 20, 22

Personal distress: items 6, 10, 13, 17, 19, 24, 27,
}
\examples{
data("iri")

}
\references{
Briganti, G., Kempenaers, C., Braun, S., Fried, E. I., & Linkowski, P. (2018). Network analysis of
empathy items from the interpersonal reactivity index in 1973
young adults. Psychiatry research, 265, 87-92.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confirm.R
\name{confirm}
\alias{confirm}
\title{GGM: Confirmatory Hypothesis Testing}
\usage{
confirm(
  Y,
  hypothesis,
  prior_sd = 0.25,
  formula = NULL,
  type = "continuous",
  mixed_type = NULL,
  iter = 25000,
  progress = TRUE,
  impute = TRUE,
  seed = 1,
  ...
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{hypothesis}{Character string. The hypothesis (or hypotheses) to be tested. See details.}

\item{prior_sd}{Numeric. Scale of the prior distribution, approximately the standard deviation
of a beta distribution (defaults to 0.25).}

\item{formula}{An object of class \code{\link[stats]{formula}}. This allows for including
control variables in the model (e.g.,, \code{~ gender * education}).}

\item{type}{Character string. Which type of data for \strong{Y} ? The options include \code{continuous},
\code{binary}, \code{ordinal}, or \code{mixed}. See the note for further details.}

\item{mixed_type}{Numeric vector of length \emph{p}. An indicator for which varibles should be treated as ranks.
(1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{iter}{Number of iterations (posterior samples; defaults to 25,000).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{impute}{Logicial. Should the missing values (\code{NA})
be imputed during model fitting (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed.}

\item{...}{Currently ignored.}
}
\value{
The returned object of class \code{confirm} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:

\itemize{

\item \code{out_hyp_prob} Posterior hypothesis probabilities.

\item \code{info} An object of class \code{BF} from the R package \strong{BFpack}.

}
}
\description{
Confirmatory hypothesis testing in GGMs. Hypotheses are expressed as equality
and/or ineqaulity contraints on the partial correlations of interest. Here the focus is \emph{not}
on determining the graph (see \code{\link{explore}}) but testing specific hypotheses related to
the conditional (in)dependence structure. These methods were introduced in
\insertCite{Williams2019_bf;textual}{BGGM}.
}
\details{
The hypotheses can be written either with the respective column names or numbers.
For example, \code{1--2} denotes the relation between the variables in column 1 and 2.
Note that these must correspond to the upper triangular elements of the correlation
matrix. This is accomplished by ensuring that the first number is smaller than the second number.
This also applies when using column names (i.e,, in reference to the column number).

\strong{One Hypothesis}:

To test whether some relations are larger than others, while others
       are expected to be equal, this can be writting as

\itemize{
\item \code{hyp <-  c(1--2 > 1--3  = 1--4 > 0)},
}

where there is an addition additional contraint that all effects are expected to be positive.
This is then compared to the complement.

\strong{More Than One Hypothesis}:

The above hypothesis can also be compared to, say, a null model by using ";"
to seperate the hypotheses, for example,

\itemize{

\item

\code{hyp <-  c(1--2 > 1--3  = 1--4 > 0; 1--2 = 1--3  = 1--4 = 0)}.


}

Any number of hypotheses can be compared this way.

\strong{Using "&"}

 It is also possible to include \code{&}. This allows for testing one constraint \bold{and}
 another contraint as one hypothesis.

\itemize{

\item \code{hyp <- c("A1--A2 > A1--A2 & A1--A3 = A1--A3")}

}

Of course, it is then possible to include additional hypotheses by separating them with ";".
Note also that the column names were used in this example (e.g., \code{A1--A2} is the relation
between those nodes).

\strong{Testing Sums}

It might also be interesting to test the sum of partial correlations. For example, that the
sum of specific relations is larger than the sum of other relations. This can be written as

\itemize{

\item \code{hyp <- c("A1--A2 + A1--A3 > A1--A4 + A1--A5;
                      A1--A2 + A1--A3 = A1--A4 + A1--A5")}

}

\strong{Potential Delays}:

There is a chance for a potentially long delay from the time the progress bar finishes
to when the function is done running. This occurs when the hypotheses require further
sampling to be tested, for example, when grouping relations
\code{c("(A1--A2, A1--A3) > (A1--A4, A1--A5)"}. This is not an error.


\strong{Controlling for Variables}:

When controlling for variables, it is assumed that \code{Y} includes \emph{only}
the nodes in the GGM and the control variables. Internally, \code{only} the predictors
that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
\code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
should be included in the GGM. An example is provided below.

\strong{Mixed Type}:

 The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
 continuous or \emph{only} discrete variables \insertCite{hoff2007extending}{BGGM}. This is based on the
 ranked likelihood which requires sampling the ranks for each variable (i.e., the data is not merely
 transformed to ranks). This is computationally expensive when there are many levels. For example,
 with continuous data, there are as many ranks as data points!

 The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
 and the "emprical" distribution is used otherwise. This is accomplished by specifying an indicator
 vector of length \emph{p}. A one indicates to use the ranks, whereas a zero indicates to "ignore"
 that variable. By default all integer variables are handled as ranks.

\strong{Dealing with Errors}:

An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):

\itemize{

\item The first is due to sampling the thresholds, especially when the data is heavily skewed.
      This can result in an ill-defined matrix. If this occurs, we recommend to first try
      decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
      change the data type to \code{type = mixed} which then estimates a copula GGM
      (this method can be used for data containing \strong{only} ordinal variable). This should
      work without a problem.

\item  The second is due to how the ordinal data are categorized. For example, if the error states
       that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
       the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.

}
}
\note{
\strong{"Default" Prior}:

 In Bayesian statistics, a default Bayes factor needs to have several properties. I refer
 interested users to \insertCite{@section 2.2 in @dablander2020default;textual}{BGGM}. In
 \insertCite{Williams2019_bf;textual}{BGGM}, some of these propteries were investigated (e.g.,
 model selection consistency). That said, we would not consider this a "default" or "automatic"
 Bayes factor and thus we encourage users to perform sensitivity analyses by varying the scale of the prior
 distribution.

 Furthermore, it is important to note there is no "correct" prior and, also, there is no need
 to entertain the possibility of a "true" model. Rather, the Bayes factor can be interpreted as
 which hypothesis best (relative to each other) predicts the observed data
 \insertCite{@Section 3.2 in @Kass1995}{BGGM}.

\strong{Interpretation of Conditional (In)dependence Models for Latent Data}:

 See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
(i.e, all data types besides \code{"continuous"})
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

##########################
### example 1: cheating ##
##########################
# Here a true hypothesis is tested,
# which shows the method works nicely
# (peeked at partials beforehand)

# data
Y <- BGGM::bfi[,1:10]

hypothesis <- c("A1--A2 < A1--A3 < A1--A4 = A1--A5")

# test cheat
test_cheat <-  confirm(Y = Y,
                       type = "continuous",
                       hypothesis  = hypothesis,
                       iter = 250,
                       progress = FALSE)

# print (probabilty of nearly 1 !)
test_cheat
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zero_order.R
\name{zero_order_cors}
\alias{zero_order_cors}
\title{Zero-Order Correlations}
\usage{
zero_order_cors(
  Y,
  type = "continuous",
  iter = 5000,
  mixed_type = NULL,
  progress = TRUE
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{type}{Character string. Which type of data for \code{Y} ? The options include \code{continuous},
\code{binary}, \code{ordinal}, or \code{mixed}. See the note for further details.}

\item{iter}{Number of iterations (posterior samples; defaults to 5000).}

\item{mixed_type}{Numeric vector. An indicator of length p for which varibles should be treated as ranks.
(1 for rank and 0 to assume normality). The default is currently to treat all integer variables as ranks
when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}
}
\value{
\itemize{

\item \code{R} An array including the correlation matrices
              (of dimensions \emph{p} by \emph{p} by \emph{iter})

\item \code{R_mean} Posterior mean of the correlations (of dimensions \emph{p} by \emph{p})
}
}
\description{
Estimate zero-order correlations for any type of data. Note zero-order refers to the fact that
no variables are controlled for (i.e., bivariate correlations). To our knowledge, this is the only Bayesian
implementation in \code{R} that can estiamte Pearson's,  tetrachoric (binary), polychoric
(ordinal with more than two cateogries), and rank based correlation coefficients.
}
\details{
\strong{Mixed Type}:

 The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
 continuous or \emph{only} discrete variables. This is based on the ranked likelihood which requires sampling
 the ranks for each variable (i.e., the data is not merely transformed to ranks). This is computationally
 expensive when there are many levels. For example, with continuous data, there are as many ranks
 as data points!

 The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
 and the "emprical" distribution is used otherwise \insertCite{hoff2007extending}{BGGM}. This is
 accomplished by specifying an indicator vector of length \emph{p}. A one indicates to use the ranks,
 whereas a zero indicates to "ignore" that variable. By default all integer variables are treated as ranks.

\strong{Dealing with Errors}:

An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):

\itemize{

\item The first is due to sampling the thresholds, especially when the data is heavily skewed.
      This can result in an ill-defined matrix. If this occurs, we recommend to first try
      decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
      change the data type to \code{type = mixed} which then estimates a copula GGM
      (this method can be used for data containing \strong{only} ordinal variable). This should
      work without a problem.

\item  The second is due to how the ordinal data are categorized. For example, if the error states
       that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
       the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.

}
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

Y <- ptsd[,1:3]

#################################
####### example 1: Pearson's ####
#################################

fit <- zero_order_cors(Y, type = "continuous",
                       iter = 250,
                       progress = FALSE)


#################################
###### example 2: polychoric ####
#################################

fit <- zero_order_cors(Y+1, type = "ordinal",
                       iter = 250,
                       progress = FALSE)


###########################
##### example 3: rank #####
###########################

fit <- zero_order_cors(Y+1, type = "mixed",
                       iter = 250,
                       progress = FALSE)

############################
## example 4: tetrachoric ##
############################

# binary data
Y <- women_math[,1:3]

fit <- zero_order_cors(Y, type = "binary",
                       iter = 250,
                       progress = FALSE)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.estimate.R
\name{select}
\alias{select}
\title{S3 \code{select} method}
\usage{
select(object, ...)
}
\arguments{
\item{object}{object of class \code{estimate} or\code{explore}}

\item{...}{not currently used}
}
\value{
\code{select} works with the following methods:
\itemize{
\item \code{\link{select.estimate}}
\item \code{\link{select.explore}}
\item \code{\link{select.ggm_compare_estimate}}
}
}
\description{
S3 select method
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precision.R
\name{precision}
\alias{precision}
\title{Precision Matrix Posterior Distribution}
\usage{
precision(object, progress = TRUE)
}
\arguments{
\item{object}{An object of class \code{estimate}.}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}
}
\value{
\itemize{

\item \code{precision_mean} The mean of the precision matrix (\code{p} by \code{p} matrix).

\item \code{precision} 3d array of dimensions \code{p} by \code{p} by \code{iter}
including \strong{unconstrained} (i.e., from th full graph)
precision matrices.

}
}
\description{
Transform the sampled correlation matrices to
precision matrices (i.e., inverse covariance matrices).
}
\note{
The estimated precision matrix is the inverse of the \strong{correlation} matrix.
}
\examples{
\donttest{
# data
Y <- ptsd

# fit model
fit <- estimate(Y)

# precision matrix
Theta <- precision(fit)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{depression_anxiety_t2}
\alias{depression_anxiety_t2}
\title{Data: Depression and Anxiety (Time 2)}
\format{
A data frame containing 403 observations (n = 7466) and 16 variables (p = 16) measured on the 4-point
likert scale.
}
\usage{
data("depression_anxiety_t2")
}
\description{
A data frame containing 403 observations (n = 403) and 16 variables (p = 16) measured on the 4-point
likert scale  (depression: 9; anxiety: 7).
}
\details{
\strong{Depression}:

\itemize{
  \item \code{PHQ1}  Little interest or pleasure in doing things?
  \item \code{PHQ2}  Feeling down, depressed, or hopeless?
  \item \code{PHQ3}  Trouble falling or staying asleep, or sleeping too much?
  \item \code{PHQ4}  Feeling tired or having little energy?
  \item \code{PHQ5}  Poor appetite or overeating?
  \item \code{PHQ6} Feeling bad about yourself — or that you are a failure or have let
                    yourself or your family down?
  \item \code{PHQ7}  Trouble concentrating on things, such as reading the newspaper or
                     watching television?
  \item \code{PHQ8} Moving or speaking so slowly that other people could have noticed? Or so
                    fidgety or restless that you have been moving a lot more than usual?
  \item \code{PHQ9}  Thoughts that you would be better off dead, or thoughts of hurting yourself
                     in some way?
}

  \strong{Anxiety}
  \itemize{



  \item \code{GAD1} Feeling nervous, anxious, or on edge
  \item \code{GAD2} Not being able to stop or control worrying
  \item \code{GAD3} Worrying too much about different things
  \item \code{GAD4} Trouble relaxing
  \item \code{GAD5} Being so restless that it's hard to sit still
  \item \code{GAD6} Becoming easily annoyed or irritable
  \item \code{GAD7} Feeling afraid as if something awful might happen
}
}
\examples{
data("depression_anxiety_t2")
labels<- c("interest", "down", "sleep",
            "tired", "appetite", "selfest",
           "concen", "psychmtr", "suicid",
           "nervous", "unctrworry", "worrylot",
           "relax", "restless", "irritable", "awful")


}
\references{
Forbes, M. K., Baillie, A. J., & Schniering, C. A. (2016). A structural equation modeling
analysis of the relationships between depression,anxiety, and sexual problems over time.
The Journal of Sex Research, 53(8), 942-954.

Forbes, M. K., Wright, A. G., Markon, K. E., & Krueger, R. F. (2019). Quantifying the reliability and replicability of psychopathology network characteristics.
Multivariate behavioral research, 1-19.

Jones, P. J., Williams, D. R., & McNally, R. J. (2019). Sampling variability is not nonreplication:
a Bayesian reanalysis of Forbes, Wright, Markon, & Krueger.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcor_mat.R
\name{pcor_mat}
\alias{pcor_mat}
\title{Extract the Partial Correlation Matrix}
\usage{
pcor_mat(object, difference = FALSE, ...)
}
\arguments{
\item{object}{A model estimated with \strong{BGGM}. All classes are supported, assuming
there is matrix to be extracted.}

\item{difference}{Logical. Should the difference be returned (defaults to \code{FALSE}) ? Note
that this assumes there is a difference (e.g., an object of class \code{ggm_compare_estimate})
and ignored otherwise.}

\item{...}{Currently ignored.}
}
\value{
The estimated partial correlation matrix.
}
\description{
Extract the partial correlation matrix (posterior mean)
from \code{\link{estimate}}, \code{\link{explore}}, \code{\link{ggm_compare_estimate}},
and \code{\link{ggm_compare_explore}} objects. It is also possible to extract the
partial correlation differences for \code{\link{ggm_compare_estimate}} and
\code{\link{ggm_compare_explore}} objects.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- ptsd[,1:5] + 1

# ordinal
fit <- estimate(Y, type = "ordinal",
                iter = 250,
                progress = FALSE)

pcor_mat(fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.estimate.R
\name{predict.estimate}
\alias{predict.estimate}
\title{Model Predictions for \code{estimate} Objects}
\usage{
\method{predict}{estimate}(
  object,
  newdata = NULL,
  summary = TRUE,
  cred = 0.95,
  iter = NULL,
  progress = TRUE,
  ...
)
}
\arguments{
\item{object}{object of class \code{estimate}}

\item{newdata}{an optional data frame for obtaining predictions (e.g., on test data)}

\item{summary}{summarize the posterior samples (defaults to \code{TRUE}).}

\item{cred}{credible interval used for summarizing}

\item{iter}{number of posterior samples (defaults to all in the object).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{currently ignored}
}
\value{
\code{summary = TRUE}: 3D array of dimensions n (observations),
        4 (posterior summary),
        p (number of nodes). \code{summary = FALSE}:
        list containing predictions for each variable
}
\description{
Model Predictions for \code{estimate} Objects
}
\examples{
\donttest{
# # data
Y <- ptsd

fit <- estimate(Y, iter = 250,
                progress = FALSE)

pred <- predict(fit,
                progress = FALSE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bggm_missing.R
\name{bggm_missing}
\alias{bggm_missing}
\title{GGM: Missing Data}
\usage{
bggm_missing(x, iter = 2000, method = "estimate", ...)
}
\arguments{
\item{x}{An object of class \code{mid} \code{\link[mice]{mice}}.}

\item{iter}{Number of iterations for each imputed dataset (posterior samples; defaults to 2000).}

\item{method}{Character string. Which method should be used (default set to \code{estimate})? The current
options are \code{"estimate"} and \code{"explore"}.}

\item{...}{Additional arguments passed to either
\code{\link{estimate}} or \code{\link{explore}}.}
}
\value{
An object of class \code{estimate} or \code{explore}.
}
\description{
Estimation and exploratory hypothesis testing with missing data.
}
\note{
Currently, \strong{BGGM} is compatible with the package \code{\link[mice]{mice}} for handling
      the missing data. This is accomplished by fitting a model for each imputed dataset
      (i.e., more than one to account for uncertainty in the imputation step) and then pooling
      the estimates.

      In a future version, an additional option will be added that allows for
      imputing the missing values during model fitting. This option will be incorporated directly into
      the \code{\link{estimate}} or \code{\link{explore}} functions, such that \code{bggm_missing} will
      always support missing data with \code{\link[mice]{mice}}.


\strong{Support}:

 There is limited support for missing data. As of version \code{2.0.0}, it is possible to
 determine the graphical structure with either  \code{\link{estimate}} or \code{\link{explore}}, in addition
 to plotting the graph with \code{\link{plot.select}}. All data types \emph{are} currently supported.

\strong{Memory Warning}:
 A model is fitted for each imputed dataset. This results in a potentially large object.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# need this package
library(mice, warn.conflicts = FALSE)

# data
Y <- ptsd[,1:5]

# matrix for indices
mat <- matrix(0, nrow = 221, ncol = 5)

# indices
indices <- which(mat == 0, arr.ind = TRUE)

# 50 NAs
Y[indices[sample(1:nrow(indices), 50),]] <- NA

# impute
x <- mice(Y, m = 5, print = FALSE)

#########################
#######   copula    #####
#########################
# rank based parital correlations

# estimate the model
fit_est <-  bggm_missing(x,
                         method = "estimate",
                         type =  "mixed",
                         iter = 250,
                         progress = FALSE)

# select edge set
E <- select(fit_est)

# plot E
plt_E <- plot(E)$plt

plt_E
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_belief_ggm.R
\name{prior_belief_ggm}
\alias{prior_belief_ggm}
\title{Prior Belief Gaussian Graphical Model}
\usage{
prior_belief_ggm(Y, prior_ggm, post_odds_cut = 3, ...)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by
\emph{p} (variables/nodes).}

\item{prior_ggm}{Matrix of dimensions \emph{p} by \emph{p}, encoding the prior
odds for including each relation in the graph (see '\code{Details}')}

\item{post_odds_cut}{Numeric. Threshold for including an edge (defaults to 3).
Note \code{post_odds} refers to posterior odds.}

\item{...}{Additional arguments passed to \code{\link{explore}}.}
}
\value{
An object including:

\itemize{

\item{\strong{adj}}: Adjacency matrix

\item{\strong{post_prob}}: Posterior probability for the
                           alternative hypothesis.

}
}
\description{
Incorporate prior information into the estimation of the
conditional dependence structure. This prior information is expressed as
the prior odds that each relation should be included in the graph.
}
\details{
Technically, the prior odds is not for including an edge in the graph,
but for (H1)/p(H0), where H1 captures the hypothesized edge size and H0 is the
null model  \insertCite{@see Williams2019_bf}{BGGM}. Accordingly, setting an
entry in \code{prior_ggm} to, say, 10, encodes a prior belief that H1 is 10 times
more likely than H0. Further, setting an entry in \code{prior_ggm} to 1 results
in equal prior odds (the default in \code{\link{select.explore}}).
}
\examples{
\donttest{
# Assume perfect prior information
# synthetic ggm
p <- 20
main <- gen_net()

# prior odds 10:1, assuming graph is known
prior_ggm <- ifelse(main$adj == 1, 10, 1)

# generate data
y <- MASS::mvrnorm(n = 200,
                   mu = rep(0, 20),
                   Sigma = main$cors)

# prior est
prior_est <- prior_belief_ggm(Y = y,
                              prior_ggm = prior_ggm,
                              progress = FALSE)

# check scores
BGGM:::performance(Estimate = prior_est$adj,
                   True = main$adj)

# default in BGGM
default_est <- select(explore(y, progress = FALSE))

# check scores
BGGM:::performance(Estimate = default_est$Adj_10,
                   True = main$adj)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.R
\name{map}
\alias{map}
\title{Maximum A Posteriori Precision Matrix}
\usage{
map(Y)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}
}
\value{
An object of class \code{map}, including the precision matrix,
        partial correlation matrix, and regression parameters.
}
\description{
Maximum A Posteriori Precision Matrix
}
\examples{
Y <- BGGM::bfi[, 1:5]

# map
map <- map(Y)
map
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_predict.R
\name{posterior_predict}
\alias{posterior_predict}
\title{Posterior Predictive Distribution}
\usage{
posterior_predict(object, iter = 1000, progress = TRUE)
}
\arguments{
\item{object}{An object of class \code{estimate} or \code{explore}}

\item{iter}{Numeric. Number of samples from the predictive distribution}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE})}
}
\value{
A 3D array containing the predicted datasets
}
\description{
Draw samples from the posterior predictive distribution.
}
\note{
Currently only implemented for \code{type = "mixed"}, \code{type = "ordinal"},
      and \code{type = "binary"}. Note the term mixed is confusing, in that it can
      be used with only, say, ordinal data. In this case, reestimate the model with \code{type = "mixed"}
      until all data types are supported.
}
\examples{
\donttest{
Y <- gss

fit <- estimate(as.matrix(Y),
                impute = TRUE,
               iter = 150, type = "mixed")

yrep <- posterior_predict(fit, iter = 100)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fisher_z2r.R
\name{fisher_z_to_r}
\alias{fisher_z_to_r}
\title{Fisher Z Back Transformation}
\usage{
fisher_z_to_r(z)
}
\arguments{
\item{z}{Fisher Z}
}
\value{
Correlation (s) (backtransformed)
}
\description{
Back tranform Fisher's Z to correlations
}
\examples{
fisher_z_to_r(0.5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd_cor1}
\alias{ptsd_cor1}
\title{Data: Post-Traumatic Stress Disorder (Sample # 1)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix was estimated from 526
individuals \insertCite{fried2018replicability}{BGGM}.
}
\details{
\itemize{
  \item Intrusive Thoughts
  \item Nightmares
  \item Flashbacks
  \item Physiological/psychological reactivity
  \item Avoidance of thoughts
  \item Avoidance of situations
  \item Amnesia
  \item Disinterest in activities
  \item Feeling detached
  \item Emotional numbing
  \item Foreshortened future
  \item Sleep problems
  \item Irritability
  \item Concentration problems
  \item Hypervigilance
  \item Startle response
}
}
\examples{

data(ptsd_cor1)

Y <- MASS::mvrnorm(n = 526,
                   mu = rep(0, 16),
                   Sigma = ptsd_cor1,
                   empirical = TRUE)

}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{estimate}
\alias{estimate}
\title{GGM: Estimation}
\usage{
estimate(
  Y,
  formula = NULL,
  type = "continuous",
  mixed_type = NULL,
  analytic = FALSE,
  prior_sd = 0.25,
  iter = 5000,
  impute = FALSE,
  progress = TRUE,
  seed = 1,
  ...
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{formula}{An object of class \code{\link[stats]{formula}}. This allows for including
control variables in the model (i.e., \code{~ gender}). See the note for further details.}

\item{type}{Character string. Which type of data for \code{Y} ? The options include \code{continuous},
\code{binary}, \code{ordinal}, or \code{mixed}. Note that mixed can be used for data with only
ordinal variables. See the note for further details.}

\item{mixed_type}{Numeric vector. An indicator of length \emph{p} for which variables should be treated as ranks.
(1 for rank and 0 to assume normality). The default is currently to treat all integer variables as ranks
when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{analytic}{Logical. Should the analytic solution be computed (default is \code{FALSE})?}

\item{prior_sd}{Scale of the prior distribution, approximately the standard deviation of a beta distribution
(defaults to 0.50).}

\item{iter}{Number of iterations (posterior samples; defaults to 5000).}

\item{impute}{Logical. Should the missing values (\code{NA})
be imputed during model fitting (defaults to \code{TRUE}) ?}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed.}

\item{...}{Currently ignored.}
}
\value{
The returned object of class \code{estimate} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:

\itemize{

\item \code{pcor_mat} Partial correltion matrix (posterior mean).

\item \code{post_samp} An object containing the posterior samples.

}
}
\description{
Estimate the conditional (in)dependence with either an analytic solution or efficiently
sampling from the posterior distribution. These methods were introduced in \insertCite{Williams2019;textual}{BGGM}.
The graph is selected with \code{\link{select.estimate}} and then plotted with \code{\link{plot.select}}.
}
\details{
The default is to draw samples from the posterior distribution (\code{analytic = FALSE}). The samples are
required for computing edge differences (see \code{\link{ggm_compare_estimate}}), Bayesian R2 introduced in
\insertCite{gelman_r2_2019;textual}{BGGM} (see \code{\link{predictability}}), etc. If the goal is
to *only* determine the non-zero effects, this can be accomplished by setting \code{analytic = TRUE}.
This is particularly useful when a fast solution is needed (see the examples in \code{\link{ggm_compare_ppc}})

\strong{Controlling for Variables}:

When controlling for variables, it is assumed that \code{Y} includes \emph{only}
the nodes in the GGM and the control variables. Internally, \code{only} the predictors
that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
\code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
should be included in the GGM. An example is provided below.

\strong{Mixed Type}:

 The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
 continuous or \emph{only} discrete variables. This is based on the ranked likelihood which requires sampling
 the ranks for each variable (i.e., the data is not merely transformed to ranks). This is computationally
 expensive when there are many levels. For example, with continuous data, there are as many ranks
 as data points!

 The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
 and the "emprical" distribution is used otherwise \insertCite{hoff2007extending}{BGGM}. This is
 accomplished by specifying an indicator vector of length \emph{p}. A one indicates to use the ranks,
 whereas a zero indicates to "ignore" that variable. By default all integer variables are treated as ranks.

\strong{Dealing with Errors}:

An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):

\itemize{

\item The first is due to sampling the thresholds, especially when the data is heavily skewed.
      This can result in an ill-defined matrix. If this occurs, we recommend to first try
      decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
      change the data type to \code{type = mixed} which then estimates a copula GGM
      (this method can be used for data containing \strong{only} ordinal variable). This should
      work without a problem.

\item  The second is due to how the ordinal data are categorized. For example, if the error states
       that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
       the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.

}

\strong{Imputing Missing Values}:

Missing values are imputed with the approach described in \insertCite{hoff2009first;textual}{BGGM}.
The basic idea is to impute the missing values with the respective posterior pedictive distribution,
given the observed data, as the model is being estimated. Note that the default is \code{TRUE},
but this ignored when there are no missing values. If set to \code{FALSE}, and there are missing
values, list-wise deletion is performed with \code{na.omit}.
}
\note{
\strong{Posterior Uncertainty}:

A key feature of \bold{BGGM} is that there is a posterior distribution for each partial correlation.
This readily allows for visiualizing uncertainty in the estimates. This feature works
with all data types and is accomplished by plotting the summary of the \code{estimate} object
(i.e., \code{plot(summary(fit))}). Several examples are provided below.



\strong{Interpretation of Conditional (In)dependence Models for Latent Data}:

See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
(i.e, all data types besides \code{"continuous"})
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

#########################################
### example 1: continuous and ordinal ###
#########################################
# data
Y <- ptsd

# continuous

# fit model
fit <- estimate(Y, type = "continuous",
                iter = 250)

# summarize the partial correlations
summ <- summary(fit)

# plot the summary
plt_summ <- plot(summary(fit))

# select the graph
E <- select(fit)

# plot the selected graph
plt_E <- plot(select(fit))


# ordinal

# fit model (note + 1, due to zeros)
fit <- estimate(Y + 1, type = "ordinal",
                iter = 250)

# summarize the partial correlations
summ <- summary(fit)

# plot the summary
plt <- plot(summary(fit))

# select the graph
E <- select(fit)

# plot the selected graph
plt_E <- plot(select(fit))

##################################
## example 2: analytic solution ##
##################################
# (only continuous)

# data
Y <- ptsd

# fit model
fit <- estimate(Y, analytic = TRUE)

# summarize the partial correlations
summ <- summary(fit)

# plot summary
plt_summ <- plot(summary(fit))

# select graph
E <- select(fit)

# plot the selected graph
plt_E <- plot(select(fit))

}

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{csws}
\alias{csws}
\title{Data: Contingencies of Self-Worth Scale (CSWS)}
\format{
A data frame with 35 variables and 680 observations (7 point Likert scale)
}
\usage{
data("csws")
}
\description{
A dataset containing items from the Contingencies of Self-Worth Scale (CSWS) scale. There are 35 variables  and
680 observations
}
\details{
\itemize{
  \item \code{1} When I think I look attractive, I feel good about myself
  \item \code{2} My self-worth is based on God's love
  \item \code{3} I feel worthwhile when I perform better than others on a task or skill.
  \item \code{4} My self-esteem is unrelated to how I feel about the way my body looks.
  \item \code{5} Doing something I know is wrong makes me lose my self-respect
  \item \code{6} I don't care if other people have a negative opinion about me.
  \item \code{7} Knowing that my family members love me makes me feel good about myself.
  \item \code{8} I feel worthwhile when I have God's love.
  \item \code{9} I can’t respect myself if others don't respect me.
  \item \code{10} My self-worth is not influenced by the quality of my relationships with my family members.
  \item \code{11} Whenever I follow my moral principles, my sense of self-respect gets a boost.
  \item \code{12} Knowing that I am better than others on a task raises my self-esteem.
  \item \code{13} My opinion about myself isn't tied to how well I do in school.
  \item \code{14} I couldn't respect myself if I didn't live up to a moral code.
  \item \code{15} I don't care what other people think of me.
  \item \code{16} When my family members are proud of me, my sense of self-worth increases.
  \item \code{17} My self-esteem is influenced by how attractive I think my face or facial features are.
  \item \code{18} My self-esteem would suffer if I didn’t have God's love.
  \item \code{19} Doing well in school gives me a sense of selfrespect.
  \item \code{20} Doing better than others gives me a sense of self-respect.
  \item \code{21} My sense of self-worth suffers whenever I think I don't look good.
  \item \code{22} I feel better about myself when I know I'm doing well academically.
  \item \code{23} What others think of me has no effect on what I think about myself.
  \item \code{24} When I don’t feel loved by my family, my selfesteem goes down.
  \item \code{25} My self-worth is affected by how well I do when I am competing with others.
  \item \code{26} My self-esteem goes up when I feel that God loves me.
  \item \code{27} My self-esteem is influenced by my academic performance.
  \item \code{28} My self-esteem would suffer if I did something unethical.
  \item \code{29} It is important to my self-respect that I have a family that cares about me.
  \item \code{30} My self-esteem does not depend on whether or not I feel attractive.
  \item \code{31} When I think that I’m disobeying God, I feel bad about myself.
  \item \code{32} My self-worth is influenced by how well I do on competitive tasks.
  \item \code{33} I feel bad about myself whenever my academic performance is lacking.
  \item \code{34} My self-esteem depends on whether or not I follow my moral/ethical principles.
  \item \code{35} My self-esteem depends on the opinions others hold of me.
  \item \code{gender} "M" (male) or "F" (female)

}
}
\note{
There are seven domains

FAMILY SUPPORT: items 7, 10, 16, 24, and 29.

COMPETITION: items 3, 12, 20, 25, and 32.

APPEARANCE: items 1, 4, 17, 21, and 30.

GOD'S LOVE: items 2, 8, 18, 26, and 31.

ACADEMIC COMPETENCE: items 13, 19, 22, 27, and 33.

VIRTUE: items 5, 11, 14, 28, and 34.

APPROVAL FROM OTHERS: items: 6, 9, 15, 23, and 35.
}
\examples{
data("csws")


}
\references{
Briganti, G., Fried, E. I., & Linkowski, P. (2019). Network analysis of Contingencies of Self-Worth
Scale in 680 university students. Psychiatry research, 272, 252-257.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_estimate.default.R
\name{ggm_compare_estimate}
\alias{ggm_compare_estimate}
\title{GGM Compare: Estimate}
\usage{
ggm_compare_estimate(
  ...,
  formula = NULL,
  type = "continuous",
  mixed_type = NULL,
  analytic = FALSE,
  prior_sd = 0.5,
  iter = 5000,
  impute = TRUE,
  progress = TRUE,
  seed = 1
)
}
\arguments{
\item{...}{Matrices (or data frames) of dimensions \emph{n} (observations) by  \emph{p} (variables).
Requires at least two.}

\item{formula}{An object of class \code{\link[stats]{formula}}. This allows for including
control variables in the model (i.e., \code{~ gender}). See the note for further details.}

\item{type}{Character string. Which type of data for \strong{Y} ? The options include \code{continuous},
\code{binary}, \code{ordinal}, or \code{continuous}. See the note for further details.}

\item{mixed_type}{Numeric vector. An indicator of length \emph{p} for which varibles should be treated as ranks.
(1 for rank and 0 to use the 'empirical' or observed distribution). The default is currently to treat all integer variables
as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.}

\item{analytic}{Logical. Should the analytic solution be computed (default is \code{FALSE})? This is only available
for continous data. Note that if \code{type = "mixed"} and \code{analytic = TRUE}, the data will
automatically be treated as continuous.}

\item{prior_sd}{The scale of the prior distribution (centered at zero), in reference to a beta distribtuion
(defaults to 0.50).
See note for further details.}

\item{iter}{Number of iterations (posterior samples; defaults to 5000).}

\item{impute}{Logicial. Should the missing values (\code{NA})
be imputed during model fitting (defaults to \code{TRUE}) ?}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed.}
}
\value{
A list of class \code{ggm_compare_estimate} containing:
 \itemize{
 \item \code{pcor_diffs} partial correlation differences (posterior distribution)
 \item \code{p} number of variable
 \item \code{info} list containing information about each group (e.g., sample size, etc.)
 \item \code{iter} number of posterior samples
 \item \code{call} \code{match.call}
 }
}
\description{
Compare partial correlations that are estimated from any number of groups. This method works for
continuous, binary, ordinal, and mixed data (a combination of categorical and continuous variables).
The approach (i.e., a difference between posterior distributions) was
described in  \insertCite{Williams2019;textual}{BGGM}.
}
\details{
This function can be used to compare the partial correlations for any number of groups.
This is accomplished with pairwise comparisons for each relation. In the case of three groups,
for example, group 1 and group 2 are compared, then group 1 and group 3 are compared, and then
group 2 and group 3 are compared. There is a full distibution for each difference that can be
summarized (i.e., \code{\link{summary.ggm_compare_estimate}}) and then visualized
(i.e., \code{\link{plot.summary.ggm_compare_estimate}}). The graph of difference is selected with
\code{\link{select.ggm_compare_estimate}}).


\strong{Controlling for Variables}:

When controlling for variables, it is assumed that \code{Y} includes \emph{only}
the nodes in the GGM and the control variables. Internally, \code{only} the predictors
that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
\code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
should be included in the GGM. An example is provided below.

\strong{Mixed Type}:

 The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
 continuous or \emph{only} discrete variables. This is based on the ranked likelihood which requires sampling
 the ranks for each variable (i.e., the data is not merely transformed to ranks). This is computationally
 expensive when there are many levels. For example, with continuous data, there are as many ranks
 as data points!

 The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
 and the "emprical" distribution is used otherwise. This is accomplished by specifying an indicator
 vector of length \emph{p}. A one indicates to use the ranks, whereas a zero indicates to "ignore"
 that variable. By default all integer variables are handled as ranks.

\strong{Dealing with Errors}:

An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):

\itemize{

\item The first is due to sampling the thresholds, especially when the data is heavily skewed.
      This can result in an ill-defined matrix. If this occurs, we recommend to first try
      decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
      change the data type to \code{type = mixed} which then estimates a copula GGM
      (this method can be used for data containing \strong{only} ordinal variable). This should
      work without a problem.

\item  The second is due to how the ordinal data are categorized. For example, if the error states
       that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
       the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.

}

\strong{Imputing Missing Values}:

Missing values are imputed with the approach described in \insertCite{hoff2009first;textual}{BGGM}.
The basic idea is to impute the missing values with the respective posterior pedictive distribution,
given the observed data, as the model is being estimated. Note that the default is \code{TRUE},
but this ignored when there are no missing values. If set to \code{FALSE}, and there are missing
values, list-wise deletion is performed with \code{na.omit}.
}
\note{
\strong{Mixed Data}:

The mixed data approach was introduced  \insertCite{@in @hoff2007extending;textual}{BGGM}
(our paper describing an extension to Bayesian hypothesis testing if forthcoming).
This is a semi-paramateric copula model based on the ranked likelihood. This is computationally
expensive when treating continuous data as ranks. The current default is to treat only integer data as ranks.
This should of course be adjusted for continous data that is skewed. This can be accomplished with the
argument \code{mixed_type}. A \code{1} in the numeric vector of length \emph{p}indicates to treat that
respective node as a rank (corresponding to the column number) and a zero indicates to use the observed
(or "emprical") data.


It is also important to note that \code{type = "mixed"} is not restricted to mixed data (containing a combination of
categorical and continuous): all the nodes can be ordinal or continuous (but again this will take some time).


\strong{Interpretation of Conditional (In)dependence Models for Latent Data}:

See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
(i.e, all data types besides \code{"continuous"})



\strong{Additional GGM Compare Methods}

Bayesian hypothesis testing is implemented in \code{\link{ggm_compare_explore}} and
\code{\link{ggm_compare_confirm}} \insertCite{Williams2019_bf}{BGGM}. The latter allows for confirmatory
hypothesis testing.  An approach based on a posterior predictive check is implemented in \code{\link{ggm_compare_ppc}}
\insertCite{williams2020comparing}{BGGM}. This provides  a 'global' test for comparing the entire GGM and a 'nodewise'
test for comparing each variable in the network \insertCite{Williams2019;textual}{BGGM}.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
                   select = -c(gender,
                               education))[,1:10]

Yfemale <- subset(Y, gender == 2,
                     select = -c(gender,
                                 education))[,1:10]

# fit model
fit <- ggm_compare_estimate(Ymale,  Yfemale,
                           type = "ordinal",
                           iter = 250,
                           prior_sd = 0.25,
                           progress = FALSE)

###########################
### example 2: analytic ###
###########################
# only continuous

# fit model
fit <- ggm_compare_estimate(Ymale, Yfemale,
                            analytic = TRUE)

# summary
summ <- summary(fit)

# plot summary
plt_summ <- plot(summary(fit))

# select
E <- select(fit)

# plot select
plt_E <- plot(select(fit))

}

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.estimate.R
\name{select.estimate}
\alias{select.estimate}
\title{Graph Selection for \code{estimate} Objects}
\usage{
\method{select}{estimate}(object, cred = 0.95, alternative = "two.sided", ...)
}
\arguments{
\item{object}{An object of class \code{estimate.default}.}

\item{cred}{Numeric. The credible interval width for selecting the graph
(defaults to 0.95; must be between 0 and 1).}

\item{alternative}{A character string specifying the alternative hypothesis. It
must be one of "two.sided" (default), "greater"  or "less".
See note for futher details.}

\item{...}{Currently ignored.}
}
\value{
The returned object of class \code{select.estimate} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:

\itemize{

\item \code{pcor_adj} Selected partial correlation matrix (weighted adjacency).
\item \code{adj} Adjacency matrix for the selected edges
\item \code{object} An object of class \code{estimate} (the fitted model).

}
}
\description{
Provides the selected graph based on credible intervals for
the partial correlations that did not contain zero
\insertCite{Williams2019}{BGGM}.
}
\details{
This package was built for the social-behavioral sciences in particular. In these applications, there is
strong theory that expects \emph{all} effects to be positive. This is known as a "positive manifold" and
this notion has a rich tradition in psychometrics. Hence, this can be incorporated into the graph with
\code{alternative = "greater"}. This results in the estimated structure including only positive edges.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi[,1:10]

# estimate
fit <- estimate(Y, iter = 250,
                progress = FALSE)


# select edge set
E <- select(fit)

}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{estimate}} and \code{\link{ggm_compare_estimate}} for several examples.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_bf.default.R
\name{summary.ggm_compare_explore}
\alias{summary.ggm_compare_explore}
\title{Summary Method for \code{ggm_compare_explore} Objects}
\usage{
\method{summary}{ggm_compare_explore}(object, col_names = TRUE, ...)
}
\arguments{
\item{object}{An object of class \code{ggm_compare_explore}.}

\item{col_names}{Logical. Should the summary include the column names (default is \code{TRUE})?
Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{summary.ggm_compare_explore}
}
\description{
Summarize the posterior hypothesis probabilities
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
                   select = -c(gender,
                               education))[,1:10]

Yfemale <- subset(Y, gender == 2,
                     select = -c(gender,
                                 education))[,1:10]

##########################
### example 1: ordinal ###
##########################

# fit model
fit <- ggm_compare_explore(Ymale,  Yfemale,
                           type = "ordinal",
                           iter = 250,
                           progress = FALSE)
# summary
summ <- summary(fit)

summ
}
}
\seealso{
\code{\link{ggm_compare_explore}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.explore.R
\name{plot.summary.select.explore}
\alias{plot.summary.select.explore}
\title{Plot \code{summary.select.explore} Objects}
\usage{
\method{plot}{summary.select.explore}(x, size = 2, color = "black", ...)
}
\arguments{
\item{x}{An object of class \code{summary.select.explore}}

\item{size}{Numeric. The size for the points (defaults to 2).}

\item{color}{Character string. The Color for the points}

\item{...}{Currently ignored}
}
\value{
A \code{ggplot} object
}
\description{
Visualize the posterior hypothesis probabilities.
}
\examples{
\donttest{
#  data
Y <- bfi[,1:10]

# fit model
fit <- explore(Y, iter = 250,
               progress = FALSE)

# edge set
E <- select(fit,
            alternative = "exhaustive")

plot(summary(E))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.estimate.R
\name{summary.coef}
\alias{summary.coef}
\title{Summarize \code{coef} Objects}
\usage{
\method{summary}{coef}(object, cred = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{coef}.}

\item{cred}{Numeric. The credible interval width for summarizing the posterior
distributions (defaults to 0.95; must be between 0 and 1).}

\item{...}{Currently ignored}
}
\value{
A list of length \emph{p} including the
        summaries for each multiple regression.
}
\description{
Summarize regression parameters with the posterior mean,
standard deviation, and credible interval.
}
\note{
See \code{\link{coef.estimate}} and \code{\link{coef.explore}} for examples.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.estimate.R
\name{coef.explore}
\alias{coef.explore}
\title{Compute Regression Parameters for \code{explore} Objects}
\usage{
\method{coef}{explore}(object, iter = NULL, progress = TRUE, ...)
}
\arguments{
\item{object}{An Object of class \code{explore}.}

\item{iter}{Number of iterations (posterior samples; defaults to the number in the object).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{coef}, containting two lists.

\itemize{

\item \code{betas} A list of length \emph{p}, each containing a \emph{p} - 1 by \code{iter} matrix of
posterior samples

\item \code{object} An object of class \code{explore} (the fitted model).
}
}
\description{
There is a direct correspondence between the inverse covariance matrix and
 multiple regression \insertCite{kwan2014regression,Stephens1998}{BGGM}. This readily allows
 for converting the GGM parameters to regression coefficients. All data types are supported.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- ptsd[,1:4]

##########################
### example 1: ordinal ###
##########################

# fit model (note + 1, due to zeros)
fit <- explore(Y + 1,
               type = "ordinal",
               iter = 250,
               progress = FALSE)

# summarize the partial correlations
reg <- coef(fit, progress = FALSE)

# summary
summ <- summary(reg)

summ
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{summary.estimate}
\alias{summary.estimate}
\title{Summary method for \code{estimate.default} objects}
\usage{
\method{summary}{estimate}(object, col_names = TRUE, cred = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{estimate}}

\item{col_names}{Logical. Should the summary include the column names (default is \code{TRUE})?
Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).}

\item{cred}{Numeric. The credible interval width for summarizing the posterior
distributions (defaults to 0.95; must be between 0 and 1).}

\item{...}{Currently ignored.}
}
\value{
A dataframe containing the summarized posterior distributions.
}
\description{
Summarize the posterior distribution of each partial correlation
with the posterior mean and standard deviation.
}
\examples{
\donttest{
# data
Y <- ptsd[,1:5]

fit <- estimate(Y, iter = 250,
                progress = FALSE)

summary(fit)

}

}
\seealso{
\code{\link{estimate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{gss}
\alias{gss}
\title{Data: 1994 General Social Survey}
\format{
A data frame containing 1190 observations (n = 1190) and 6 variables (p = 6) measured on the binary scale
        \insertCite{fowlkes1988evaluating}{BGGM}. The variable descriptions were copied from
        \insertCite{@section 4, @hoff2007extending;textual}{BGGM}
}
\usage{
data("gss")
}
\description{
A data frame containing 1002 rows and 7 variables measured on various scales,
including binary and ordered cateogrical (with varying numbers of categories).
There are also missing values in each variable

\itemize{
  \item \code{Inc}  Income of the respondent in 1000s of dollars, binned into 21 ordered categories.
  \item \code{DEG}   Highest degree ever obtained (none, HS, Associates, Bachelors, or Graduate)
  \item \code{CHILD}  Number of children ever had.
  \item \code{PINC}  Financial status of respondent's parents when respondent was 16 (on a 5-point scale).
  \item \code{PDEG}  Maximum of mother's and father's highest degree
  \item \code{PCHILD}  Number of siblings of the respondent plus one
  \item \code{AGE} Age of the respondent in years.
}
}
\examples{
data("gss")
}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.estimate.R
\name{predict.var_estimate}
\alias{predict.var_estimate}
\title{Model Predictions for \code{var_estimate} Objects}
\usage{
\method{predict}{var_estimate}(object, summary = TRUE, cred = 0.95, iter = NULL, progress = TRUE, ...)
}
\arguments{
\item{object}{object of class \code{var_estimate}}

\item{summary}{summarize the posterior samples (defaults to \code{TRUE}).}

\item{cred}{credible interval used for summarizing}

\item{iter}{number of posterior samples (defaults to all in the object).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{...}{Currently ignored}
}
\value{
The predicted values for each regression model.
}
\description{
Model Predictions for \code{var_estimate} Objects
}
\examples{
\donttest{
# data
Y <- subset(ifit, id == 1)[,-1]

# fit model with alias (var_estimate also works)
fit <- var_estimate(Y, progress = FALSE)

# fitted values
pred <- predict(fit, progress = FALSE)

# predicted values (1st outcome)
pred[,,1]

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_bf.default.R
\name{plot.summary.ggm_compare_explore}
\alias{plot.summary.ggm_compare_explore}
\title{Plot \code{summary.ggm_compare_explore} Objects}
\usage{
\method{plot}{summary.ggm_compare_explore}(x, size = 2, color = "black", ...)
}
\arguments{
\item{x}{An object of class \code{summary.ggm_compare_explore}}

\item{size}{Numeric. The size of the points (defaults to 2).}

\item{color}{Character string. The color of the points
(defaults to \code{"black"}).}

\item{...}{Currently ignored.}
}
\value{
A \code{ggplot} object
}
\description{
Visualize the posterior hypothesis probabilities.
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi

# males and females
Ymale <- subset(Y, gender == 1,
                   select = -c(gender,
                               education))[,1:10]

Yfemale <- subset(Y, gender == 2,
                     select = -c(gender,
                                 education))[,1:10]

##########################
### example 1: ordinal ###
##########################

# fit model
fit <- ggm_compare_explore(Ymale,  Yfemale,
                           type = "ordinal",
                           iter = 250,
                           progress = FALSE)
# summary
summ <- summary(fit)

plot(summ)
}
}
\seealso{
\code{\link{ggm_compare_explore}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare_ppc.default.R
\name{ggm_compare_ppc}
\alias{ggm_compare_ppc}
\title{GGM Compare: Posterior Predictive Check}
\usage{
ggm_compare_ppc(
  ...,
  test = "global",
  iter = 5000,
  FUN = NULL,
  custom_obs = NULL,
  loss = TRUE,
  progress = TRUE
)
}
\arguments{
\item{...}{At least two matrices (or data frames) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{test}{Which test should be performed (defaults to \code{"global"}) ? The options include
\code{global} and \code{nodewise}.}

\item{iter}{Number of replicated datasets used to construct the predictivie distribution
(defaults to 5000).}

\item{FUN}{An optional function for comparing GGMs that returns a number. See \strong{Details}.}

\item{custom_obs}{Number corresponding to the observed score for comparing the GGMs. This is
required if a function is provided in \code{FUN}. See \strong{Details}.}

\item{loss}{Logical. If a function is provided, is the measure a "loss function"
(i.e., a large score is bad thing). This determines how the \emph{p}-value
is computed. See \strong{Details}.}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}
}
\value{
The returned object of class \code{ggm_compare_ppc} contains a lot of information that
        is used for printing and plotting the results. For users of \strong{BGGM}, the following
        are the useful objects:

\code{test = "global"}

\itemize{

\item \code{ppp_jsd} posterior predictive \emph{p}-values (JSD).

\item \code{ppp_sse} posterior predictive \emph{p}-values (SSE).

\item \code{predictive_jsd} list containing the posterior predictive distributions (JSD).

\item \code{predictive_sse} list containing the posterior predictive distributions (SSE).

\item \code{obs_jsd} list containing the observed error (JSD).

\item \code{obs_sse} list containing the observed error (SSE).

}


\code{test = "nodewise"}

\itemize{

\item \code{ppp_jsd} posterior predictive \emph{p}-values (JSD).

\item \code{predictive_jsd} list containing the posterior predictive distributions (JSD).

\item \code{obs_jsd} list containing the observed error (JSD).

}

\code{FUN = f()}

\itemize{

\item \code{ppp_custom} posterior predictive \emph{p}-values (custom).

\item \code{predictive_custom} posterior predictive distributions (custom).

\item \code{obs_custom} observed error (custom).

}
}
\description{
Compare GGMs with a posterior predicitve check \insertCite{gelman1996posterior}{BGGM}.
This method was introduced in \insertCite{williams2020comparing;textual}{BGGM}. Currently,
there is a \code{global} (the entire GGM) and a \code{nodewise} test. The default
is to compare GGMs with respect to the posterior predictive distribution of Kullback
Leibler divergence and the sum of squared errors. It is also possible to compare the
GGMs with a user defined test-statistic.
}
\details{
The \code{FUN} argument allows for a user defined test-statisic (the measure used to compare the GGMs).
The function must include only two agruments, each of which corresponds to a dataset. For example,
\code{f <- function(Yg1, Yg2)}, where each Y is dataset of dimensions \emph{n} by \emph{p}. The
groups are then compare within the function, returning a single number. An example is provided below.

Further, when using a custom function care must be taken when specifying the argument \code{loss}.
We recommended to visualize the results with \code{plot} to ensure the \emph{p}-value was computed
in the right direction.
}
\note{
\strong{Interpretation}:

The primary test-statistic is symmetric KL-divergence that is termed Jensen-Shannon divergence (JSD).
This is in essence a likelihood ratio that provides the "distance" between two multivariate normal
distributions. The basic idea is to (1) compute the posterior predictive distribution, assuming group equality
(the null model). This provides the error that we would expect to see under the null model; (2) compute
JSD for the observed groups; and (3) compare the observed JSD to the posterior predictive distribution,
from which a posterior predictive \emph{p}-value is computed.

For the \code{global} check, the sum of squared error is also provided.
This is computed from the partial correlation matrices and it is analagous
to the strength test in \insertCite{van2017comparing;textual}{BGGM}. The \code{nodewise}
test compares the posterior predictive distribution for each node. This is based on the correspondence
between the inverse covariance matrix and multiple regresssion \insertCite{kwan2014regression,Stephens1998}{BGGM}.

If the null model is \code{not} rejected, note that this does \code{not} provide evidence for equality!
Further, if the null model is rejected, this means that the assumption of group equality is not tenable--the
groups are different.

\strong{Alternative Methods}:

There are several methods in \strong{BGGM} for comparing groups. See
\code{\link{ggm_compare_estimate}} (posterior differences for the
partial correlations), \code{\link{ggm_compare_explore}} (exploratory hypothesis testing),
and \code{\link{ggm_compare_confirm}} (confirmatory hypothesis testing).
}
\examples{

\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi

#############################
######### global ############
#############################


# males
Ym <- subset(Y, gender == 1,
             select = - c(gender, education))

# females

Yf <- subset(Y, gender == 2,
             select = - c(gender, education))


global_test <- ggm_compare_ppc(Ym, Yf,
                               iter = 250)

global_test


#############################
###### custom function ######
#############################
# example 1

# maximum difference van Borkulo et al. (2017)

f <- function(Yg1, Yg2){

# remove NA
x <- na.omit(Yg1)
y <- na.omit(Yg2)

# nodes
p <- ncol(Yg1)

# identity matrix
I_p <- diag(p)

# partial correlations

pcor_1 <- -(cov2cor(solve(cor(x))) - I_p)
pcor_2 <- -(cov2cor(solve(cor(y))) - I_p)

# max difference
max(abs((pcor_1[upper.tri(I_p)] - pcor_2[upper.tri(I_p)])))

}

# observed difference
obs <- f(Ym, Yf)

global_max <- ggm_compare_ppc(Ym, Yf,
                              iter = 250,
                              FUN = f,
                              custom_obs = obs,
                              progress = FALSE)

global_max


# example 2
# Hamming distance (squared error for adjacency)

f <- function(Yg1, Yg2){

# remove NA
x <- na.omit(Yg1)
y <- na.omit(Yg2)

# nodes
p <- ncol(x)

# identity matrix
I_p <- diag(p)

fit1 <-  estimate(x, analytic = TRUE)
fit2 <-  estimate(y, analytic = TRUE)

sel1 <- select(fit1)
sel2 <- select(fit2)

sum((sel1$adj[upper.tri(I_p)] - sel2$adj[upper.tri(I_p)])^2)

}

# observed difference
obs <- f(Ym, Yf)

global_hd <- ggm_compare_ppc(Ym, Yf,
                            iter = 250,
                            FUN = f,
                            custom_obs  = obs,
                            progress = FALSE)

global_hd


#############################
########  nodewise ##########
#############################

nodewise <- ggm_compare_ppc(Ym, Yf, iter = 250,
                           test = "nodewise")

nodewise

}

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.explore.R
\name{summary.select.explore}
\alias{summary.select.explore}
\title{Summary Method for \code{select.explore} Objects}
\usage{
\method{summary}{select.explore}(object, col_names = TRUE, ...)
}
\arguments{
\item{object}{object of class \code{select.explore}.}

\item{col_names}{Logical.}

\item{...}{Currently ignored.}
}
\value{
a data frame including the posterior mean, standard deviation,
and posterior hypothesis probabilities for each relation.
}
\description{
Summary Method for \code{select.explore} Objects
}
\examples{
\donttest{
#  data
Y <- bfi[,1:10]

# fit model
fit <- explore(Y, iter = 250,
               progress = FALSE)

# edge set
E <- select(fit,
            alternative = "exhaustive")

summary(E)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression_summary.R
\name{regression_summary}
\alias{regression_summary}
\title{Summarary Method for Multivariate or Univarate Regression}
\usage{
regression_summary(object, cred = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{estimate}}

\item{cred}{Numeric. The credible interval width for summarizing the posterior
distributions (defaults to 0.95; must be between 0 and 1).}

\item{...}{Currently ignored}
}
\value{
A list of length \emph{p} including the
        summaries for each regression.
}
\description{
Summarary Method for Multivariate or Univarate Regression
}
\examples{
\donttest{
# note: iter = 250 for demonstrative purposes

# data
Y <- bfi

Y <- subset(Y, select = c("E5", "N5",
                          "gender", "education"))


fit_mv_ordinal <- estimate(Y, formula = ~ gender + as.factor(education),
                           type = "ordinal",
                           iter = 250,
                           progress = FALSE)

regression_summary(fit_mv_ordinal)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ifit}
\alias{ifit}
\title{Data: ifit Intensive Longitudinal Data}
\format{
A data frame containing 197 observations and 8 variables. The data have been used in
\insertCite{o2020use}{BGGM} and  \insertCite{williams2019bayesian}{BGGM}
}
\usage{
data("ifit")
}
\description{
A data frame containing 8 variables and nearly 200 observations. There are
two subjects, each of which provided data every data for over 90 days. Six variables are from
the PANAS scale (positive and negative affect), the daily number of steps, and the subject id.

\itemize{
  \item \code{id} Subject id
  \item \code{interested}
  \item \code{disinterested}
  \item \code{excited}
  \item \code{upset}
  \item \code{strong}
  \item \code{stressed}
  \item \code{steps} steps recorded by a fit bit
}
}
\examples{
data("ifit")
}
\references{
\insertAllCited{}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BGGM-package.R
\docType{package}
\name{BGGM-package}
\alias{BGGM-package}
\title{BGGM:  Bayesian Gaussian Graphical Models}
\description{
The \code{R} package \strong{BGGM} provides tools for making Bayesian inference in
Gaussian graphical models (GGM). The methods are organized around two general approaches for
Bayesian inference: (1) estimation \insertCite{Williams2019}{BGGM} and (2) hypothesis testing
\insertCite{Williams2019_bf}{BGGM}. The key distinction is that the former focuses on either the
posterior or posterior predictive distribution, whereas the latter focuses on model comparison
with the Bayes factor.

The methods in \strong{BGGM} build upon existing algorithms that are well-known in the literature.
The central contribution of \strong{BGGM} is to extend those approaches:

\enumerate{

\item Bayesian estimation with the novel matrix-F prior distribution \insertCite{Mulder2018}{BGGM}.

 \itemize{

 \item Estimation \code{\link{estimate}}.

 }


\item Bayesian hypothesis testing with the novel matrix-F prior distribution \insertCite{Mulder2018}{BGGM}.

 \itemize{

 \item Exploratory hypothesis testing \code{\link{explore}}.

 \item Confirmatory hypothesis  testing \code{\link{confirm}}.

 }

\item Comparing GGMs \insertCite{williams2020comparing}{BGGM}

 \itemize{

 \item Partial correlation differences \code{\link{ggm_compare_estimate}}.

 \item Posterior predictive check \code{\link{ggm_compare_ppc}}.

 \item Exploratory hypothesis testing \code{\link{ggm_compare_explore}}.

 \item Confirmatory hypothesis testing \code{\link{ggm_compare_confirm}}.


 }


\item Extending inference beyond the conditional (in)dependence structure

\itemize{

\item Predictability with Bayesian variance explained \insertCite{gelman_r2_2019}{BGGM}
      \code{\link{predictability}}.

\item Posterior uncertainty in the partial correlations \code{\link{estimate}}.

\item Custom Network Statistics \code{\link{roll_your_own}}.


}

}

Furthermore, the computationally intensive tasks are written in \code{c++} via the \code{R}
package \strong{Rcpp} \insertCite{eddelbuettel2011rcpp}{BGGM} and the \code{c++}
library \strong{Armadillo} \insertCite{sanderson2016armadillo}{BGGM}, there are plotting functions
for each method, control variables can be included in the model, and there is support for
missing values \code{\link{bggm_missing}}.

\bold{Supported Data Types}:

\itemize{

\item Continuous: The continuous method was described \insertCite{@in  @Williams2019_bf;textual}{BGGM}.

\item Binary: The binary method builds directly upon \insertCite{@in @talhouk2012efficient;textual}{BGGM},
      that, in turn, built upon the approaches of \insertCite{lawrence2008bayesian;textual}{BGGM} and
      \insertCite{webb2008bayesian;textual}{BGGM} (to name a few).

\item Ordinal: Ordinal data requires sampling thresholds. There are two approach included in \bold{BGGM}: (1)
the customary approach described in \insertCite{@in @albert1993bayesian;textual}{BGGM} (the default) and
the 'Cowles' algorithm described in \insertCite{@in @cowles1996accelerating;textual}{BGGM}.


\item Mixed: The mixed data (a combination of discrete and continuous) method was introduced
\insertCite{@in @hoff2007extending;textual}{BGGM}. This is a semi-parametric copula model
(i.e., a copula GGM) based on the ranked likelihood. Note that this can be used for data
consisting entirely of ordinal data.

}

\bold{Additional Features}:

 The primary focus of \code{BGGM} is Gaussian graphical modeling (the inverse covariance matrix).
 The residue is a suite of useful methods not explicitly for GGMs:

 \enumerate{

 \item Bivariate correlations for binary (tetrachoric), ordinal (polychoric), mixed (rank based),
       and continuous (Pearson's) data \code{\link{zero_order_cors}}.

 \item Multivariate regression for binary (probit), ordinal (probit),
       mixed (rank likelihood), and continous data (\code{\link{estimate}}).

 \item Multiple regression for binary (probit), ordinal (probit),
       mixed (rank likelihood), and continuous data (e.g., \code{\link{coef.estimate}}).
 }

\strong{Note on Conditional (In)dependence Models for Latent Data}:

All of the data types (besides continuous) model latent data. That is, unobserved
(latent) data is assumed to be Gaussian. For example, a tetrachoric correlation
(binary data) is a special case of a polychoric correlation (ordinal data).
Both capture relations between "theorized normally distributed continuous
\strong{latent} variables" (\href{https://en.wikipedia.org/wiki/Polychoric_correlation}{Wikipedia}).
In both instances, the corresponding partial correlation between observed variables is conditioned
on the remaining variables in the \emph{latent} space. This implies that interpretation
is similar to continuous data, but with respect to latent variables. We refer interested users
to \insertCite{@page 2364, section 2.2, in  @webb2008bayesian;textual}{BGGM}.


\strong{High Dimensional Data?}

\strong{BGGM} was built specifically for social-behavioral scientists. Of course,
the methods can be used by all researchers. However, there is currently \emph{not} support
for high-dimensional data (i.e., more variables than observations) that are common
place in the genetics literature. These data are rare in the social-behavioral sciences.
In the future, support for high-dimensional data may be added to \strong{BGGM}.
}
\references{
\insertAllCited{}
}
