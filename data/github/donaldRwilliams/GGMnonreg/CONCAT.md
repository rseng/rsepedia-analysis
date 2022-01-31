---
title: 'GGMnonreg: Non-Regularized Gaussian Graphical Models in R'
tags:
- Graphical models
- partial correlations
- Mixed graphical model
- Ising model
authors:
  - name: Donald R. Williams
    affiliation: "1, 2"
affiliations:
 - name: Department of Psychology, University of California, Davis
   index: 1
 - name: NWEA, Portland, USA
   index: 2
citation_author: Williams
date: 08 November 2021
year: 2021
bibliography: inst/REFERENCES.bib
---
  
# Summary
Studying complex relations in multivariate datasets is a common task across the sciences. Cognitive neuroscientists model brain connectivity with the goal of unearthing functional and structural associations between cortical regions [@ortiz_2015]. In clinical psychology, researchers wish to better understand the intricate web of symptom interrelations that underlie mental health disorders 
[@mcnally_2016; @borsboom_small_world]. To this end, graphical modeling has emerged as an oft-used tool in the chest of scientific inquiry. The basic idea is to characterize multivariate relations by learning the conditional dependence structure. The cortical regions or symptoms are *nodes* and the featured connections linking nodes are *edges* that graphically represent the conditional dependence structure.

Graphical modeling is quite common in fields with wide data, that is, when there are more variables ($p$) than observations ($n$). Accordingly, many regularization-based approaches have been developed for those kinds of data. There are key drawbacks of regularization, including, but not limited to, the fact that obtaining a valid measure of parameter uncertainty is very (very) difficult [@Buhlmann2014] and there can
be an inflated false positive rate [see for example, @williams2019nonregularized].

# Statement of Need
More recently, graphical modeling has emerged in psychology (Epskamp et al. 2018), where the data is typically long or low-dimensional ($p < n$; @williams2019nonregularized, @williams_rethinking). The primary purpose of 
**GGMnonreg** is to provide methods that were specifically designed for low-dimensional data (e.g., those common in the social-behavioral sciences).

## Supported Models
* Gaussian graphical model (GGM). The following data types are supported.
  + Gaussian 
  + Ordinal
  + Binary
* Ising model [@marsman_2018]
* Mixed graphical model

## Additional methods
The following are also included

* Expected network replicability [@williams2020learning]
* Compare Gaussian graphical models
* Measure of parameter uncertainty [@williams_2021_conf]
* Edge inclusion "probabilities" [e.g., Figure 6.4 in  @Hastie2015]
* Network visualization
* Constrained precision matrix [the network, given an assumed graph, see p. 631 in @hastie2009elements]
* Predictability [variance explained for each node, @Haslbeck2018]

## Gaussian graphical Model
The following estimates a GGM for 5 post-traumatic stress disorder 
(PTSD) symptoms [@armour2017network]:

```{r}
fit <- ggm_inference(Y = ptsd[,1:5], 
                     boot = FALSE)

fit
#>           1         2         3         4         5
#> 1 0.0000000 0.2262934 0.0000000 0.3335737 0.1547986
#> 2 0.2262934 0.0000000 0.4993419 0.0000000 0.0000000
#> 3 0.0000000 0.4993419 0.0000000 0.2205442 0.1841798
#> 4 0.3335737 0.0000000 0.2205442 0.0000000 0.3407634
#> 5 0.1547986 0.0000000 0.1841798 0.3407634 0.0000000
```

### Predictability
It is common to then estimate "predictability", which corresponds to $R^2$
for each node in the network. In **GGMnonreg**, this is implemented with the 
following code:

```{r}
predictability(fit)

#>   Estimate Est.Error Ci.lb Ci.ub
#> 1     0.45      0.05  0.35  0.54
#> 2     0.50      0.05  0.41  0.59
#> 3     0.55      0.04  0.47  0.64
#> 4     0.50      0.05  0.41  0.59
#> 5     0.46      0.05  0.37  0.55
```


## Ising Model
An Ising model is for binary data. The PTSD symptoms can be binary, indicating
the symptom was either present or absent. This network is estimated with:

```{r}
# make binary
Y <- ifelse(ptsd[,1:5] == 0, 0, 1)

# fit model
fit <- ising_search(Y, IC = "BIC", 
                    progress = FALSE)

fit
#>          1        2        3        4        5
#> 1 0.000000 1.439583 0.000000 1.273379 0.000000
#> 2 1.439583 0.000000 1.616511 0.000000 1.182281
#> 3 0.000000 1.616511 0.000000 1.716747 1.077322
#> 4 1.273379 0.000000 1.716747 0.000000 1.662550
#> 5 0.000000 1.182281 1.077322 1.662550 0.000000
```

## Network Replicability
Recently, the topic of replicability has captivated the network literature. 
To this end, I developed an analytic solution to estimate network replicability [@williams2020learning].

The first step is to define a "true" partial correlation network. As an example, 
I generate a synthetic partial correlation matrix, and then compute expected
network replicability.

```{r}
# edges between 0.05 and 0.25
main <- gen_net(p = 20, 
                lb = 0.05, 
                ub = 0.25)

# enr                
enr(main$pcors, 
    n = 500, 
    replications = 4)

#> Average Replicability: 0.53 
#> Average Number of Edges: 30 (SD = 2.12) 
#> 
#> ----
#> 
#> Cumulative Probability: 
#> 
#>  prop.edges edges Pr(R > prop.edges)
#>         0.0     0               1.00
#>         0.1     6               1.00
#>         0.2    11               1.00
#>         0.3    17               1.00
#>         0.4    23               1.00
#>         0.5    28               0.78
#>         0.6    34               0.02
#>         0.7    40               0.00
#>         0.8    46               0.00
#>         0.9    51               0.00
----
Pr(R > prop.edges):
probability of replicating more than the
correpsonding proportion (and number) of edges
```
On average, we can expect to replicate roughly half of the 
edges in four replication attempts, where replication is defined as detecting 
a given edge in each attempt. Further, the probability of replicating more than 
70% of the edges is zero.


## Network Visualization
A key aspect of graphical modeling is visualizing the conditional dependence structure. To this end, 
**GGMnonreg** makes network plots with **ggplot2** [@ggplotpackage].

```
plot(fit, 
     node_names = colnames(Y), 
     edge_magnify = 2)
```
![Conditional Dependence Structure](man/figures/figure_1.png)

# Acknowledgements
DRW was supported by a National Science Foundation Graduate Research Fellowship
under Grant No. 1650042


# References


<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/nonreg_hex.png" width = 250 />

# GGMnonreg: Non-regularized Gaussian Graphical Models

[![CRAN
Version](http://www.r-pkg.org/badges/version/GGMnonreg)](https://cran.r-project.org/package=GGMnonreg)
[![Downloads](https://cranlogs.r-pkg.org/badges/GGMnonreg)](https://cran.r-project.org/package=GGMnonreg)
[![CircleCI build
status](https://circleci.com/gh/donaldRwilliams/GGMnonreg.svg?style=shield)](https://circleci.com/gh/donaldRwilliams/GGMnonreg)
[![zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.5668161.svg)](https://doi.org/10.5281/zenodo.5668161)

The goal of **GGMnonreg** is to estimate non-regularized graphical
models. Note that the title is a bit of a misnomer, in that Ising and
mixed graphical models are also supported.

Graphical modeling is quite common in fields with *wide* data, that is,
when there are more variables than observations. Accordingly, many
regularization-based approaches have been developed for those kinds of
data. There are key drawbacks of regularization when the goal is
inference, including, but not limited to, the fact that obtaining a
valid measure of parameter uncertainty is very (very) difficult.

More recently, graphical modeling has emerged in psychology (Epskamp et
al. 2018), where the data is typically long or low-dimensional \[*p*
&lt; *n*; Donald R. Williams et al. (2019); Donald R. Williams and Rast
(2019)\]. The primary purpose of **GGMnonreg** is to provide methods
specifically for low-dimensional data (e.g., those common to
psychopathology networks).

## Supported Models

-   Gaussian graphical model. The following data types are supported.
    -   Gaussian
    -   Ordinal
    -   Binary
-   Ising model (Marsman et al. 2017)
-   Mixed graphical model

## Additional methods

The following are also included

-   Expected network replicability (Donald R. Williams 2020)
-   Compare Gaussian graphical models
-   Measure of parameter uncertainty (Donald R. Williams et al. 2019)
-   Edge inclusion “probabilities”
-   Network visualization
-   Constrained precision matrix (the network, given an assumed graph)
-   Predictability (variance explained)

## Installation

To install the latest release version (1.1.0) from CRAN use

``` r
install.packages("GGMnonreg")    
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/GGMnonreg")
```

## Ising

An Ising model is fitted with the following

``` r
library(GGMnonreg)

# make binary
Y <- ifelse(ptsd[,1:5] == 0, 0, 1)

# fit model
fit <- ising_search(Y, IC = "BIC", 
                    progress = FALSE)

fit
#>          1        2        3        4        5
#> 1 0.000000 1.439583 0.000000 1.273379 0.000000
#> 2 1.439583 0.000000 1.616511 0.000000 1.182281
#> 3 0.000000 1.616511 0.000000 1.716747 1.077322
#> 4 1.273379 0.000000 1.716747 0.000000 1.662550
#> 5 0.000000 1.182281 1.077322 1.662550 0.000000
```

Note the same code, more or less, is also used for GGMs and mixed
graphical models.

## Predictability

It is common to compute predictability, or variance explained, for each
node in the network. An advantage of **GGMnonreg** is that a measure of
uncertainty is also provided.

``` r
# data
Y <- na.omit(bfi[,1:5])

# fit model
fit <- ggm_inference(Y, boot = FALSE)

# predictability
predictability(fit)
#>   Estimate Est.Error Ci.lb Ci.ub
#> 1     0.13      0.01  0.11  0.15
#> 2     0.33      0.01  0.30  0.36
#> 3     0.38      0.01  0.35  0.41
#> 4     0.18      0.01  0.15  0.20
#> 5     0.29      0.01  0.26  0.32
```

## Parameter Uncertainty

Confidence intervals for each relation are obtained with

``` r
# data
Y <- na.omit(bfi[,1:5])

# fit model
fit <- ggm_inference(Y, boot = TRUE, 
                     method = "spearman", 
                     B = 100, progress = FALSE)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

confint(fit)
#>              2.5%        97.5%
#>  [1,] -0.29205233 -0.212951642
#>  [2,] -0.14926375 -0.075919715
#>  [3,]  0.25508466  0.326055944
#>  [4,] -0.03778214  0.032287534
#>  [5,]  0.12219960  0.204592760
#>  [6,]  0.13387070  0.206412243
#>  [7,] -0.07641808  0.006226351
#>  [8,]  0.10997005  0.195192621
#>  [9,]  0.34388225  0.418927425
#> [10,]  0.08130533  0.153621543
```

These can then be plotted with, say, **ggplot2** (left to the user).

## Edge Inclusion

When mining data, or performing an automatic search, it is difficult to
make inference on the network parameters (e.g., confidence are not
easily computed). To summarize data mining, **GGMnonreg** provides edge
inclusion “probabilities” (proportion bootstrap samples for which each
relation was detected).

``` r
# data
Y <- na.omit(bfi[,1:5])

# fit model
fit <-  eip(Y, method = "spearman", 
            B  = 100, progress = FALSE)

fit
#>     eip
#> 1  1.00
#> 2  1.00
#> 3  1.00
#> 4  0.06
#> 5  1.00
#> 6  1.00
#> 7  0.31
#> 8  1.00
#> 9  1.00
#> 10 1.00
```

Note in all cases, the provided estimates correspond to the
upper-triangular elements of the network.

## Expected Network Replicability

**GGMnonreg** allows for computing expected network replicability (ENR),
i.e., the number of effects that will be detected in any number of
replications. This is an analytic solution.

The first step is defining a true network

``` r
# first make the true network
cors <- cor(GGMnonreg::ptsd)

# inverse
inv <- solve(cors)

# partials
pcors <-  -cov2cor(inv)

# set values to zero
pcors <- ifelse(abs(pcors) < 0.05, 0, pcors)
```

Then obtain ENR

``` r
fit_enr <- enr(net = pcors, n = 500, replications = 2)

fit_enr
#> Average Replicability: 0.45 
#> Average Number of Edges: 53 (SD = 3.7) 
#> 
#> ----
#> 
#> Cumulative Probability: 
#> 
#>  prop.edges edges Pr(R > prop.edges)
#>         0.0     0               1.00
#>         0.1    12               1.00
#>         0.2    23               1.00
#>         0.3    35               1.00
#>         0.4    47               0.91
#>         0.5    58               0.05
#>         0.6    70               0.00
#>         0.7    82               0.00
#>         0.8    94               0.00
#>         0.9   105               0.00
#> ----
#> Pr(R > prop.edges):
#> probability of replicating more than the
#> correpsonding proportion (and number) of edges
```

Note this is inherently frequentist. As such, over the long run, 45 % of
the edges will be replicated on average. Then we can further infer that,
in hypothetical replication attempts, more than half of the edges will
be replicated only 5 % of the time.

ENR can also be plotted

``` r
plot_enr(fit_enr)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="75%" />

### Intuition

Here is the basic idea of ENR

``` r
# location of edges
index <- which(pcors[upper.tri(diag(20))] != 0)

# convert network into correlation matrix
diag(pcors) <- 1
cors_new <- corpcor::pcor2cor(pcors)

# replicated edges
R <- NA

# increase 1000 to, say, 5,000
for(i in 1:1000){

  # two replications
  Y1 <- MASS::mvrnorm(500, rep(0, 20), cors_new)
  Y2 <- MASS::mvrnorm(500, rep(0, 20), cors_new)

  # estimate network 1
  fit1 <- ggm_inference(Y1, boot = FALSE)

  # estimate network 2
  fit2 <- ggm_inference(Y2, boot = FALSE)

  # number of replicated edges (detected in both networks)
  R[i] <- sum(
    rowSums(
      cbind(fit1$adj[upper.tri(diag(20))][index],
            fit2$adj[upper.tri(diag(20))][index])
    ) == 2)
}
```

Notice that replication of two networks is being assessed over the long
run. In other words, if we draw two random samples, what is the expected
replicability.

Compare analytic to simulation

``` r
# combine simulation and analytic
cbind.data.frame(
  data.frame(simulation = sapply(seq(0, 0.9, 0.1), function(x) {
    mean(R > round(length(index) * x) )
  })),
  data.frame(analytic = round(fit_enr$cdf, 3))
)
#>    simulation analytic
#> 1       1.000    1.000
#> 2       1.000    1.000
#> 3       1.000    1.000
#> 4       1.000    1.000
#> 5       0.897    0.912
#> 6       0.055    0.052
#> 7       0.000    0.000
#> 8       0.000    0.000
#> 9       0.000    0.000
#> 10      0.000    0.000

# average replicability (simulation)
mean(R / length(index))
#> [1] 0.4482051

# average replicability (analytic)
fit_enr$ave_pwr
#> [1] 0.4485122
```

ENR works with any correlation, assuming there is an estimate of the
standard error.

## Network plot

``` r
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# get info for plotting
plot(fit, edge_magnify = 5)
#> Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
#> "none")` instead.
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="75%" />

## Bug Reports, Feature Requests, and Contributing

Bug reports and feature requests can be made by opening an issue on
[Github](https://github.com/donaldRwilliams/GGMnonreg/issues). To
contribute towards the development of **GGMnonreg**, you can start a
branch with a pull request and we can discuss the proposed changes
there.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Epskamp2018ggm" class="csl-entry">

Epskamp, Sacha, Lourens J. Waldorp, Rene Mottus, and Denny Borsboom.
2018. “<span class="nocase">The Gaussian Graphical Model in
Cross-Sectional and Time-Series Data</span>.” *Multivariate Behavioral
Research* 53 (4): 453–80.
<https://doi.org/10.1080/00273171.2018.1454823>.

</div>

<div id="ref-marsman_2018" class="csl-entry">

Marsman, M, D Borsboom, J Kruis, S Epskamp, R van Bork, L J Waldorp, By
Taylor, L Waldorp, H L J van der Maas, and G Maris. 2017. “<span
class="nocase">An Introduction to Network Psychometrics: Relating Ising
Network Models to Item Response Theory Models</span>.” *Taylor &
Francis* 53 (1): 15–35. <https://doi.org/10.1080/00273171.2017.1379379>.

</div>

<div id="ref-williams2020learning" class="csl-entry">

Williams, Donald R. 2020. “Learning to Live with Sampling Variability:
Expected Replicability in Partial Correlation Networks.” *PsyArXiv*.
<https://doi.org/10.31234/osf.io/fb4sa>.

</div>

<div id="ref-williams_rethinking" class="csl-entry">

Williams, Donald R., and Philippe Rast. 2019. “<span class="nocase">Back
to the basics: Rethinking partial correlation network
methodology</span>.” *British Journal of Mathematical and Statistical
Psychology*. <https://doi.org/10.1111/bmsp.12173>.

</div>

<div id="ref-williams2019nonregularized" class="csl-entry">

Williams, Donald R., Mijke Rhemtulla, Anna C Wysocki, and Philippe Rast.
2019. “On Nonregularized Estimation of Psychological Networks.”
*Multivariate Behavioral Research* 54 (5): 719–50.
<https://doi.org/10.1080/00273171.2019.1575716>.

</div>

</div>
---
output: github_document
bibliography: inst/REFERENCES.bib
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "75%"
)
```

<img src="man/figures/nonreg_hex.png" width = 250 />

# GGMnonreg: Non-regularized Gaussian Graphical Models

[![CRAN Version](http://www.r-pkg.org/badges/version/GGMnonreg)](https://cran.r-project.org/package=GGMnonreg)
[![Downloads](https://cranlogs.r-pkg.org/badges/GGMnonreg)](https://cran.r-project.org/package=GGMnonreg)
[![CircleCI build status](https://circleci.com/gh/donaldRwilliams/GGMnonreg.svg?style=shield)](https://circleci.com/gh/donaldRwilliams/GGMnonreg)
[![zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.5668161.svg)](https://doi.org/10.5281/zenodo.5668161)


The goal of **GGMnonreg** is to estimate non-regularized graphical models. Note that 
the title is a bit of a misnomer, in that Ising and mixed graphical models are also supported.

Graphical modeling is quite common in fields with *wide* data, that is, when there are more 
variables than observations. Accordingly, many regularization-based approaches have been developed for those kinds of data. There are key drawbacks of regularization when the goal is inference, 
including, but not limited to, the fact that obtaining a valid measure of parameter uncertainty is very (very) difficult.

More recently, graphical modeling has emerged in psychology [@Epskamp2018ggm], where the data 
is typically long or low-dimensional [*p* < *n*; @williams2019nonregularized; @williams_rethinking]. The primary purpose of  **GGMnonreg** is to provide methods specifically for low-dimensional data 
(e.g., those common to psychopathology networks).

## Supported Models
* Gaussian graphical model. The following data types are supported.
  + Gaussian 
  + Ordinal
  + Binary
* Ising model [@marsman_2018]
* Mixed graphical model

## Additional methods
The following are also included

* Expected network replicability [@williams2020learning]
* Compare Gaussian graphical models
* Measure of parameter uncertainty [@williams2019nonregularized]
* Edge inclusion "probabilities"
* Network visualization
* Constrained precision matrix (the network, given an assumed graph)
* Predictability (variance explained)

## Installation
To install the latest release version (1.1.0) from CRAN use

```r
install.packages("GGMnonreg")    
```

You can install the development version from [GitHub](https://github.com/) with:
``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/GGMnonreg")
```

## Ising
An Ising model is fitted with the following
```{r}
library(GGMnonreg)

# make binary
Y <- ifelse(ptsd[,1:5] == 0, 0, 1)

# fit model
fit <- ising_search(Y, IC = "BIC", 
                    progress = FALSE)

fit
```

Note the same code, more or less, is also used for GGMs and mixed graphical models.

## Predictability
It is common to compute predictability, or variance explained, for each node in the network.
An advantage of **GGMnonreg** is that a measure of uncertainty is also provided.

```{r}
# data
Y <- na.omit(bfi[,1:5])

# fit model
fit <- ggm_inference(Y, boot = FALSE)

# predictability
predictability(fit)
```


## Parameter Uncertainty
Confidence intervals for each relation are obtained with
```{r}
# data
Y <- na.omit(bfi[,1:5])

# fit model
fit <- ggm_inference(Y, boot = TRUE, 
                     method = "spearman", 
                     B = 100, progress = FALSE)

confint(fit)
```
These can then be plotted with, say, **ggplot2** (left to the user).

## Edge Inclusion
When mining data, or performing an automatic search, it is difficult to make inference on the
network parameters (e.g., confidence are not easily computed). To summarize data mining,
**GGMnonreg** provides edge inclusion "probabilities" (proportion bootstrap samples for 
which each relation was detected).

```{r}
# data
Y <- na.omit(bfi[,1:5])

# fit model
fit <-  eip(Y, method = "spearman", 
            B  = 100, progress = FALSE)

fit
```
Note in all cases, the provided estimates correspond to the upper-triangular elements
of the network.

## Expected Network Replicability
**GGMnonreg** allows for computing expected network replicability (ENR), i.e., the number of 
effects that will be detected in any number of replications. This is an analytic solution.

The first step is defining a true network
```{r}
# first make the true network
cors <- cor(GGMnonreg::ptsd)

# inverse
inv <- solve(cors)

# partials
pcors <-  -cov2cor(inv)

# set values to zero
pcors <- ifelse(abs(pcors) < 0.05, 0, pcors)
```

Then obtain ENR
```{r}
fit_enr <- enr(net = pcors, n = 500, replications = 2)

fit_enr
```
Note this is inherently frequentist. As such, over the long run, 45 % of the edges will be replicated on average. Then we can further infer that, in hypothetical replication attempts, more than half of the edges
will be replicated only 5 % of the time.

ENR can also be plotted
```{r}
plot_enr(fit_enr)
```

### Intuition
Here is the basic idea of ENR
```{r}
# location of edges
index <- which(pcors[upper.tri(diag(20))] != 0)

# convert network into correlation matrix
diag(pcors) <- 1
cors_new <- corpcor::pcor2cor(pcors)

# replicated edges
R <- NA

# increase 1000 to, say, 5,000
for(i in 1:1000){

  # two replications
  Y1 <- MASS::mvrnorm(500, rep(0, 20), cors_new)
  Y2 <- MASS::mvrnorm(500, rep(0, 20), cors_new)

  # estimate network 1
  fit1 <- ggm_inference(Y1, boot = FALSE)

  # estimate network 2
  fit2 <- ggm_inference(Y2, boot = FALSE)

  # number of replicated edges (detected in both networks)
  R[i] <- sum(
    rowSums(
      cbind(fit1$adj[upper.tri(diag(20))][index],
            fit2$adj[upper.tri(diag(20))][index])
    ) == 2)
}
```
Notice that replication of two networks is being assessed over the long run. In other words,
if we draw two random samples, what is the expected replicability.

Compare analytic to simulation
```{r}
# combine simulation and analytic
cbind.data.frame(
  data.frame(simulation = sapply(seq(0, 0.9, 0.1), function(x) {
    mean(R > round(length(index) * x) )
  })),
  data.frame(analytic = round(fit_enr$cdf, 3))
)

# average replicability (simulation)
mean(R / length(index))

# average replicability (analytic)
fit_enr$ave_pwr
```

ENR works with any correlation, assuming there is an estimate of the standard error.

## Network plot
```{r, message=FALSE}
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# get info for plotting
plot(fit, edge_magnify = 5)
```

## Bug Reports, Feature Requests, and Contributing
Bug reports and feature requests can be made by opening an issue on [Github](https://github.com/donaldRwilliams/GGMnonreg/issues). To contribute towards
the development of **GGMnonreg**, you can start a branch with a pull request and we can 
discuss the proposed changes there.

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd_cor2}
\alias{ptsd_cor2}
\title{Cor: Post-Traumatic Stress Disorder (Sample # 2)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix
was estimated from 365 individuals \insertCite{fried2018replicability}{GGMnonreg}.
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
% Please edit documentation in R/enr.R
\name{enr}
\alias{enr}
\title{Expected Network Replicability}
\usage{
enr(net, n, alpha = 0.05, replications = 2, type = "pearson")
}
\arguments{
\item{net}{True network of dimensions \emph{p} by \emph{p}.}

\item{n}{Integer. The samples size, assumed equal in the replication
attempts.}

\item{alpha}{The desired significance level (defaults to \code{0.05}). Note that
1 - alpha corresponds to specificity.}

\item{replications}{Integer. The desired number of replications.}

\item{type}{Character string. Which type of correlation coefficients
to be computed. Options include \code{"pearson"} (default)
and \code{"spearman"}.}
}
\value{
An list of class \code{enr} including the following:

\itemize{

\item{\strong{ave_power}}: Average power.

\item{\strong{cdf}}: cumulative distribution function.

\item{\strong{p_s}}: Power for each edge, or the probability
of success for a given trial.

\item{\strong{p}}: Number of nodes.

\item{\strong{n_nonzero}}: Number of edges.

\item{\strong{n}}: Sample size.

\item{\strong{replication}}: Replication attempts.

\item{\strong{var_pwr}}: Variance of power.

\item{\strong{type}}: Type of correlation coefficient.

}
}
\description{
Investigate network replicability for any kind of
             partial correlation, assuming there is an analytic
             solution for the standard error (e.g., Pearson's or Spearman's).
}
\note{
This method was introduced in
\insertCite{williams2020learning;textual}{GGMnonreg}.

The basic idea is to determine the replicability of edges in a
partial correlation network. This requires defining the true
network, which can include edges of various sizes, and then
solving for the proportion of edges that are expected
to be replicated (e.g. in two, three, or four replication attempt).
}
\examples{
\donttest{
# (1) define partial correlation network

# correlations from ptsd symptoms
cors <- cor(GGMnonreg::ptsd)

# inverse
inv <- solve(cors)

# partials
pcors <-  -cov2cor(inv)

# set values to zero
# (this is the partial correlation network)
pcors <- ifelse(abs(pcors) < 0.05, 0, pcors)


# compute ENR in two replication attempts
fit_enr <- enr(net = pcors,
               n = 500,
               replications = 2)


# intuition for the method:
# The above did not require simulation, and here I use simulation
# for the same purpose.

# location of edges
# (where the edges are located in the network)
index <- which(pcors[upper.tri(diag(20))] != 0)

# convert network a into correlation matrix
# (this is needed to simulate data)
diag(pcors) <- 1
cors_new <- corpcor::pcor2cor(pcors)

# replicated edges
# (store the number of edges that were replicated)
R <- NA

# simulate how many edges replicate in two attempts
# (increase 100 to, say, 5,000)
for(i in 1:100){

  # two replications
  Y1 <- MASS::mvrnorm(500, rep(0, 20), cors_new)
  Y2 <- MASS::mvrnorm(500, rep(0, 20), cors_new)

  # estimate network 1
  fit1 <- ggm_inference(Y1, boot = FALSE)

  # estimate network 2
  fit2 <- ggm_inference(Y2, boot = FALSE)

  # number of replicated edges (detected in both networks)
  R[i] <- sum(
    rowSums(
      cbind(fit1$adj[upper.tri(diag(20))][index],
            fit2$adj[upper.tri(diag(20))][index])
    ) == 2)
}


# combine simulation and analytic
cbind.data.frame(
  data.frame(simulation = sapply(seq(0, 0.9, 0.1), function(x) {
    mean(R > round(length(index) * x) )
  })),
  data.frame(analytic = round(fit_enr$cdf, 3))
)

# now compare simulation to the analytic solution
# average replicability (simulation)
mean(R / length(index))

# average replicability (analytic)
fit_enr$ave_pwr
}

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{women_math}
\alias{women_math}
\title{Data: Women and Mathematics}
\format{
A data frame containing 1190 observations (n = 1190) and 6 variables (p = 6) measured on the binary scale
        \insertCite{fowlkes1988evaluating}{GGMnonreg}. These data have been analyzed in \insertCite{tarantola2004mcmc;textual}{GGMnonreg}
        and in \insertCite{madigan1994model}{GGMnonreg}. The variable descriptions were copied from  (section 5.2 )
        \insertCite{@section 5.2, @talhouk2012efficient}{GGMnonreg}
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
\title{Cor: Post-Traumatic Stress Disorder  (Sample # 4)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix
was estimated from 965 individuals \insertCite{fried2018replicability}{GGMnonreg}.
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
% Please edit documentation in R/plot_graph.R
\name{get_graph}
\alias{get_graph}
\title{Get Graph}
\usage{
get_graph(x)
}
\arguments{
\item{x}{An object of class \code{ggmnonreg}}
}
\value{
A list including two matrices (the weighted adjacency and adjacency matrices)
}
\description{
Extract the necessary ingredients to visualize the conditional
             dependence structure.
}
\examples{
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# get info for plotting
get_graph(fit)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd_cor3}
\alias{ptsd_cor3}
\title{Cor: Post-Traumatic Stress Disorder  (Sample # 3)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix
was estimated from 926 individuals \insertCite{fried2018replicability}{GGMnonreg}.
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
% Please edit documentation in R/constrained.R
\name{constrained}
\alias{constrained}
\title{Precision Matrix with Known Graph}
\usage{
constrained(Sigma, adj)
}
\arguments{
\item{Sigma}{Covariance matrix}

\item{adj}{An adjacency matrix that encodes the constraints,
where a zero indicates that element should be zero.}
}
\value{
A list containing the following:

\itemize{

\item{\strong{Theta}}: Inverse of the covariance matrix
(precision matrix), that encodes the conditional
(in)dependence structure.

\item{\strong{Sigma}}: Covariance matrix.

\item{\strong{wadj}}: Weighted adjacency matrix, corresponding
to the partial correlation network.

}
}
\description{
Compute the maximum likelihood estimate of the precision matrix,
given a known graphical structure (i.e., an adjacency matrix).
This approach was originally described in
"The Elements of Statistical Learning"
\insertCite{@see pg. 631, @hastie2009elements}{GGMnonreg}.
}
\note{
The algorithm is written in \code{c++}, and should scale to high dimensions.

Note there are a variety of algorithms for this purpose. Simulation
studies indicated that this approach is both accurate and computationally
efficient \insertCite{@HFT therein, @emmert2019constrained}{GGMnonreg}
}
\examples{
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# constrain to zero
constrained_graph <- constrained(cor(Y), fit$adj)

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eip.R
\name{eip}
\alias{eip}
\title{Edge Inclusion "Probability"}
\usage{
eip(Y, method = "pearson", B = 1000, progress = TRUE)
}
\arguments{
\item{Y}{The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes).}

\item{method}{Character string. Which type of correlation coefficients
to be computed. Options include \code{"pearson"} (default),
\code{"kendall"}, \code{"spearman"}, and \code{"polychoric"}.}

\item{B}{Integer. Number of bootstrap replicates (defaults to \code{1000}).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE})?}
}
\value{
An object of class \code{eip}, including a matrix of edge inclusion
"probabilities".
}
\description{
Compute the proportion of bootstrap samples that each relation was selected,
 corresponding to an edge inclusion "probability".
}
\details{
The order is the upper-triangular.
}
\note{
In the context of regression, this general approach was described in
see Figure 6.4. \insertCite{@see Figure 6.4, @Hastie2015;textual}{GGMnonreg}. In this case,
the selection is based on classical hypothesis testing instead of L1-regularization.
}
\examples{
\donttest{
# data
Y <- ptsd

# eip
fit_eip <- eip(Y, method = "spearman")

# print
fit_eip
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggm_compare.R
\name{ggm_compare}
\alias{ggm_compare}
\title{Compare Gaussian Graphical Models}
\usage{
ggm_compare(Yg1, Yg2, method = "spearman", alpha = 0.05)
}
\arguments{
\item{Yg1}{The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes) for
group one.}

\item{Yg2}{The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes) for
group two.}

\item{method}{Character string. Which type of correlation coefficients
to be computed. Options include \code{"pearson"} (default),
\code{"kendall"}, \code{"spearman"}, and \code{"polychoric"}.}

\item{alpha}{The desired significance level (defaults to \code{0.05}). Note that
1 - alpha corresponds to specificity.}
}
\value{
An object of class \code{ggm_compare} including:

\itemize{

\item{\strong{adj}}: Adjacency matrix, where a 1
indicates a difference.

\item{\strong{wadj}}: Weighted adjacency matrix
(partial correlation differences that were significantly different)

\item{\strong{cis}}: Confidence intervals for the partial correlation
differences.

}
}
\description{
Establish whether each of the corresponding edges
             are significantly different in two groups
}
\examples{
\donttest{
# data

Yg1 <- na.omit(subset(bfi, gender == 1)[,1:10])
Yg2 <- na.omit(subset(bfi, gender == 2)[,1:10])

# compare relations
fit <- ggm_compare(Yg1, Yg2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_graph.R
\name{plot.ggmnonreg}
\alias{plot.ggmnonreg}
\title{Network Plot for \code{graph} Objects}
\usage{
\method{plot}{ggmnonreg}(
  x,
  layout = "circle",
  neg_col = "#D55E00",
  pos_col = "#009E73",
  edge_magnify = 1,
  node_size = 10,
  palette = 2,
  node_names = NULL,
  node_groups = NULL,
  ...
)
}
\arguments{
\item{x}{An object of class \code{graph} obtained from \code{\link[GGMnonreg]{get_graph}}.}

\item{layout}{Character string. Which graph layout (defaults is \code{circle}) ?
See \link[sna]{gplot.layout}.}

\item{neg_col}{Character string. Color for the positive edges
(defaults to a colorblind friendly red).}

\item{pos_col}{Character string.  Color for the negative edges
(defaults to a colorblind friendly green).}

\item{edge_magnify}{Numeric. A value that is multiplied by the edge weights. This increases (> 1) or
decreases (< 1) the line widths (defaults to 1).}

\item{node_size}{Numeric. The size of the nodes (defaults to \code{10}).}

\item{palette}{A character string sepcifying the palette for the \code{groups}.
(default is \code{Set3}). See \href{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}{palette options here}.}

\item{node_names}{Character string. Names for nodes of length \emph{p}.}

\item{node_groups}{A character string of length \emph{p} (the number of nodes in the model).
This indicates groups of nodes that should be the same color
(e.g., "clusters" or "communities").}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{ggplot}
}
\description{
Visualize the conditional (in)dependence structure.
}
\examples{
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# plot graph
plot(fit)
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
\title{Network Predictability (R2)}
\usage{
predictability(x, ci = 0.95)
}
\arguments{
\item{x}{An object of class \code{ggm_inference}}

\item{ci}{Numeric. The confidence interval to be computed (defaults to \code{0.95}).}
}
\value{
An object of class \code{predictability}, including a matrix of R2
for each node.
}
\description{
Network Predictability (R2)
}
\note{
Predictability is variance explained for each node in the
      network \insertCite{Haslbeck2018}{GGMnonreg}.
}
\examples{
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# predictability
r2 <- predictability(fit)

# print
r2

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cis.R
\name{confint.ggm_inference}
\alias{confint.ggm_inference}
\title{Extract Confidence Intervals from \code{ggm_inference} Objects}
\usage{
\method{confint}{ggm_inference}(object, ...)
}
\arguments{
\item{object}{An object of class \code{ggm_inference}.}

\item{...}{Currently ignored.}
}
\value{
A matrix including bootstrap confidence intervals.
}
\description{
Extract Confidence Intervals from \code{ggm_inference} Objects
}
\examples{
\donttest{
#  data
Y <- ptsd

# eip
fit <- ggm_inference(Y, method = "spearman",
boot = TRUE, B = 100)

# cis
confint(fit)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enr.R
\name{plot_enr}
\alias{plot_enr}
\title{Plot \code{enr} Objects}
\usage{
plot_enr(x, iter = 1e+05, fill = "#009E73", alpha = 0.5, ...)
}
\arguments{
\item{x}{An object of class \code{enr}.}

\item{iter}{Integer. How many draws from the
Poisson-binomial distribution (defaults to 1,000)?}

\item{fill}{Which color to fill the density?}

\item{alpha}{Numeric (between 0 and 1). The transparency
for the density.}

\item{...}{Currently ignored}
}
\value{
An object of class \code{ggplot}
}
\description{
Plot the probability mass function for ENR.
}
\examples{
\donttest{
# correlations
cors <- cor(GGMnonreg::ptsd)

# inverse
inv <- solve(cors)

# partials
pcors <-  -cov2cor(inv)

# set values to zero
pcors <- ifelse(abs(pcors) < 0.05, 0, pcors )

est <- enr(net = pcors, n = 500, replications = 2)

# plot
plot_enr(est)
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
package for further details \insertCite{psych}{GGMnonreg}.
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
A dataset containing items that measure Post-traumatic stress disorder symptoms \insertCite{armour2017network}{GGMnonreg}.
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
% Please edit documentation in R/ggm_search.R
\name{ggm_search}
\alias{ggm_search}
\title{Gaussian graphical model: automated search}
\usage{
ggm_search(
  x,
  IC = "BIC",
  type = "neighborhood_selection",
  method = "forward",
  n = NULL
)
}
\arguments{
\item{x}{A data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes) or
a correlation matrix of dimensions \emph{p} by \emph{p}.}

\item{IC}{Character string. The desired information criterion. Options include
\code{"AIC"} and \code{"BIC"} (default).}

\item{type}{Character string. Which search method should be used? The options included
\code{"regression"} and \code{"approx_L0"}. See details.}

\item{method}{Character string. The desired subset selection method
Options includes \code{"forward"} (default), \code{"backward"},
and \code{"exhaustive"}.}

\item{n}{Integer. Sample size. Required if a correlation matrix is provided.}
}
\value{
An object of class \code{ggm_search} including:

\itemize{

\item{\strong{wadj}}: Weighted adjacency matrix, corresponding
to the partial correlation network.

\item{\strong{adj}}: Adjacency matrix (detected effects).

\item{\strong{pcors}}: Partial correlations.

\item{\strong{n}}: Sample size.

\item{\strong{p}}: Number of nodes.

\item{\strong{Y}}: Data.

}
}
\description{
Data mining to learn the graph.
}
\details{
\code{type = "neighborhood_selection"} was described in
\insertCite{williams2019nonregularized;textual}{GGMnonreg}
and \code{type = "approx_L0"} was described in \insertCite{williams2020beyond;textual}{GGMnonreg}.
The penalty for \code{type = "approx_L0"} is called seamless L0 \insertCite{dicker2013variable}{GGMnonreg}
}
\note{
\code{type = "neighborhood_selection"} employs multiple regression to estimate
the graph (requires the data), whereas \code{type = "approx_L0"} directly estimates
the precision matrix (data or a correlation matrix are acceptable). If
data is provided and \code{type = "approx_L0"}, by default Pearson correlations are
used. For another correlation coefficient, provide the desired correlation matrix.

\code{type = "approx_L0"} is a continuous approximation to (non-regularized)
best subset model selection. This is accomplished by using regularization, but
the penalty (approximately) mimics non-regularized estimation.
}
\examples{
\donttest{
# data
Y <- ptsd

# search data
fit <- ggm_search(Y)
}

}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggmnonreg_print.R
\name{print.ggmnonreg}
\alias{print.ggmnonreg}
\title{Print \code{ggmnonreg} Object}
\usage{
\method{print}{ggmnonreg}(x, ...)
}
\arguments{
\item{x}{An object of class \code{ggmnonreg}}

\item{...}{Currently ignored}
}
\value{
No return value.
}
\description{
Print \code{ggmnonreg} Object
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
\title{Cor: Post-Traumatic Stress Disorder (Sample # 1)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix was estimated from 526
individuals \insertCite{fried2018replicability}{GGMnonreg}.
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
% Please edit documentation in R/mixed_search.R
\name{mixed_search}
\alias{mixed_search}
\title{Mixed Graphical Model: automated search}
\usage{
mixed_search(Y, data_type = NULL, IC = "BIC")
}
\arguments{
\item{Y}{A data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes)}

\item{data_type}{Vector of length \emph{p}. The type of data, with options of
"b" (binary), "p" (Poisson), and "g" (Gaussian).}

\item{IC}{Character string. The desired information criterion. Options include
\code{"AIC"} and \code{"BIC"} (default).}
}
\value{
An object of class \code{mixed_search} including

\itemize{

\item{\strong{wadj}}: Weighted adjacency matrix, corresponding
to the partial correlation network.

\item{\strong{adj}}: Adjacency matrix (detected effects).

\item{\strong{pcors}}: Partial correlations.

\item{\strong{n}}: Sample size.

\item{\strong{p}}: Number of nodes.

\item{\strong{Y}}: Data.

}
}
\description{
Data mining to learn the graph.
}
\details{
Only backwards selection is currently implemented.
         Only an adjacency matrix is provided.
}
\examples{
\donttest{
# data
Y <- ifelse( ptsd[,1:5] == 0, 0, 1)

# search data (ising model)
fit <- mixed_search(Y, data_type = rep("b", 5))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GGMnonreg-package.R
\docType{package}
\name{GGMnonreg-package}
\alias{GGMnonreg-package}
\title{GGMnonreg:  Non-Regularized Gaussian Graphical Models}
\description{
The goal of \strong{GGMnonreg} is to estimate non-regularized
graphical models. Note that the title is a bit of a misnomer, in that Ising
and mixed graphical models are also supported. Graphical modeling is quite
common in fields with \emph{wide} data, that is, when there are more variables
than observations. Accordingly, many regularization-based approaches have been
developed for those kinds of data. There are key drawbacks of regularization
when the goal is inference, including, but not limited to, the fact that
obtaining a valid measure of parameter uncertainty is very (very) difficult.

More recently, graphical modeling has emerged in psychology,
where the data are typically long or low-dimensional
\insertCite{williams_rethinking,williams2019nonregularized}{GGMnonreg}.
The primary purpose of  \strong{GGMnonreg} is to provide methods specifically
for low-dimensional data


\strong{Supported Models}

\itemize{
 \item{Gaussian graphical model. The following data types are supported.}
 \itemize{
 \item{Gaussian}
 \item{Ordinal}
 \item{Binary}
 }
 \item{Ising model}
 \item{Mixed graphical model}
}

\strong{Additional Methods}

\itemize{
\item{Expected network replicability} \insertCite{williams2020learning}{GGMnonreg}

\item{Compare Gaussian graphical models}

\item{Measure of uncertainty} \insertCite{williams_2021_conf}{GGMnonreg}

\item{Edge inclusion "probabilities"}

\item{Network visualization}

\item{Constrained precision matrix (the network, given an assumed graph)}

\item{Predictability (variance explained)}
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{gss}
\alias{gss}
\title{Data: 1994 General Social Survey}
\format{
A data frame containing 1190 observations (n = 1190) and 6 variables (p = 6) measured on the binary scale
        \insertCite{fowlkes1988evaluating}{GGMnonreg}. The variable descriptions were copied from
        \insertCite{@section 4, @hoff2007extending;textual}{GGMnonreg}
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
% Please edit documentation in R/ising_search.R
\name{ising_search}
\alias{ising_search}
\title{Ising: automated search}
\usage{
ising_search(Y, IC = "BIC", progress = TRUE)
}
\arguments{
\item{Y}{A data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes).}

\item{IC}{Character string. The desired information criterion. Options include
\code{"AIC"} and \code{"BIC"} (default).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE})?}
}
\value{
An object of class \code{ising_search} including:

\itemize{

\item{\strong{wadj}}: Weighted adjacency matrix, corresponding
to the partial correlation network.

\item{\strong{adj}}: Adjacency matrix (detected effects).

\item{\strong{pcors}}: Partial correlations.

\item{\strong{n}}: Sample size.

\item{\strong{p}}: Number of nodes.

\item{\strong{Y}}: Data.

}
}
\description{
Data mining to learn the graph of binary variables with an Ising model
             \insertCite{lenz1920beitrvsge,ising1925beitrag}{GGMnonreg}.
}
\details{
Currently only backwards selection is currently implemented.
}
\note{
For an excellent overview of the Ising model see \insertCite{marsman2018introduction;textual}{GGMnonreg}.
}
\examples{
\donttest{
# data
Y <- ifelse( ptsd[,1:5] == 0, 0, 1)

# search data
fit <- ising_search(Y)
}
}
\references{
\insertAllCited{}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ifit}
\alias{ifit}
\title{Data: ifit Intensive Longitudinal Data}
\format{
A data frame containing 197 observations and 8 variables. The data have been used in
\insertCite{o2020use}{GGMnonreg} and  \insertCite{williams2019bayesian}{GGMnonreg}
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
% Please edit documentation in R/ggm_inference.R
\name{ggm_inference}
\alias{ggm_inference}
\title{Gaussian graphical model: statistical inference}
\usage{
ggm_inference(
  Y,
  alpha = 0.05,
  control_precision = FALSE,
  boot = TRUE,
  B = 1000,
  cores = 1,
  method = "pearson",
  progress = TRUE
)
}
\arguments{
\item{Y}{The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes).}

\item{alpha}{The desired significance level (defaults to \code{0.05}). Note that
1 - alpha corresponds to specificity.}

\item{control_precision}{Logical. Should precision (i.e., 1 - false discovery rate)
be controlled at the level alpha (defaults to \code{FALSE}) ?}

\item{boot}{Logical. Should a non-parametric bootstrap be employed (defaults to \code{TRUE})?}

\item{B}{Integer. Number of bootstrap replicates (defaults to \code{1000}).}

\item{cores}{Integer. Number of cores to be used when executing in parallel
(defaults to 1).}

\item{method}{Character string. Which type of correlation coefficients
to be computed. Options include \code{"pearson"} (default),
\code{"kendall"}, \code{"spearman"}, and \code{"polychoric"}.}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE})?}
}
\value{
An object of class \code{ggm_inference} including:

\itemize{

\item{\strong{wadj}}: Weighted adjacency matrix, corresponding
to the partial correlation network.

\item{\strong{adj}}: Adjacency matrix (detected effects).

\item{\strong{pcors}}: Partial correlations.

\item{\strong{n}}: Sample size.

\item{\strong{p}}: Number of nodes.

\item{\strong{Y}}: Data.

}
}
\description{
Learn the conditional dependence structure with null hypothesis
             significance testing. This provides a valid measure of parameter
             uncertainty.
}
\examples{
\donttest{

Y <- ptsd

fit <- ggm_inference(Y)

}

}
