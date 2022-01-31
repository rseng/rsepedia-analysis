
<!-- README.md is generated from README.Rmd. Please edit that file -->

# registr <img src="README_files/figures/registr.png" align="right" height = "150" />

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/registr)](https://cran.r-project.org/package=registr)
[![](http://cranlogs.r-pkg.org/badges/grand-total/registr?color=green)](https://cran.r-project.org/package=registr)
[![](https://travis-ci.org/julia-wrobel/registr.svg?branch=master)](https://travis-ci.org/julia-wrobel/registr)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/julia-wrobel/registr?branch=master&svg=true)](https://ci.appveyor.com/project/julia-wrobel/registr)
[![Codecov test
coverage](https://codecov.io/gh/julia-wrobel/registr/branch/master/graph/badge.svg)](https://codecov.io/gh/julia-wrobel/registr/coverage.svg?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02964/status.svg)](https://doi.org/10.21105/joss.02964)
[![R-CMD-check](https://github.com/julia-wrobel/registr/workflows/R-CMD-check/badge.svg)](https://github.com/julia-wrobel/registr/actions)
<!-- badges: end -->

Registration for incomplete exponential family functional data.

-   Authors: [Julia Wrobel](http://juliawrobel.com), [Alexander
    Bauer](https://www.en.stablab.stat.uni-muenchen.de/people/doktoranden/bauer1/index.html),
    [Erin McDonnell](http://eimcdonnell.com/), and [Jeff
    Goldsmith](https://jeffgoldsmith.com/)
-   License: [MIT](https://opensource.org/licenses/MIT). See the
    [LICENSE](LICENSE) file for details
-   Version: 2.1

### What it does

------------------------------------------------------------------------

Functional data analysis is a set of tools for understanding patterns
and variability in data where the basic unit of observation is a curve
measured over some domain such as time or space. An example is an
accelerometer study where intensity of physical activity was measured at
each minute over 24 hours for 50 subjects. The data will contain 50
curves, where each curve is the 24-hour activity profile for a
particular subject.

Classic functional data analysis assumes that each curve is continuous
or comes from a Gaussian distribution. However, applications with
exponential family functional data – curves that arise from any
exponential family distribution, and have a smooth latent mean – are
increasingly common. For example, take the accelerometer data just
mentioned, but assume researchers are interested in *sedentary behavior*
instead of *activity intensity*. At each minute over 24 hours they
collect a binary measurement that indicates whether a subject was active
or inactive (sedentary). Now we have a *binary curve* for each subject –
a trajectory where each time point can take on a value of 0 or 1. We
assume the binary curve has a smooth latent mean, which in this case is
interpreted as the probability of being active at each minute over 24
hours. This is a example of exponential family functional data.

Often in a functional dataset curves have similar underlying patterns
but the main features of each curve, such as the minimum and maximum,
have shifts such that the data appear misaligned. This misalignment can
obscure patterns shared across curves and produce messy summary
statistics. Registration methods reduce variability in functional data
and clarify underlying patterns by aligning curves.

This package implements statistical methods for registering exponential
family functional data. The basic methods are described in more detail
in our [paper](http://juliawrobel.com/Downloads/registration_ef.pdf) and
were further adapted to (potentially) incomplete curve settings where
(some) curves are not observed from the very beginning and/or until the
very end of the common domain. For details on the incomplete curve
methodology and how to use it see the corresponding package vignette.
Instructions for installing the software and using it to register
simulated binary data are provided below.

### Installation

------------------------------------------------------------------------

To install from `CRAN`, please use:

``` r
install.packages("registr")
```

To install the latest version directly from Github, please use:

``` r
install.packages("devtools")
devtools::install_github("julia-wrobel/registr")
```

The `registr` package includes vignettes with more details on package
use and functionality. To install the latest version and pull up the
vignettes please use:

``` r
devtools::install_github("julia-wrobel/registr", build_vignettes = TRUE)
vignette(package = "registr")
```

### How to use it

------------------------------------------------------------------------

This example registers simulated binary data. More details on the use of
the package can be found in the vignettes mentioned above.

The code below uses `registr::simulate_unregistered_curves()` to
simulate curves for 100 subjects with 200 timepoints each, observed over
domain (0,1). All curves have similar structure but the location of the
peak is shifted. On the observed domain *t*<sup>\*</sup> the curves are
unregistered (misaligned). On the domain *t* the curves are registered
(aligned).

``` r
library(registr)

registration_data = simulate_unregistered_curves(I = 100, D = 200, seed = 2018)
```

The plot below shows the unregistered curves and registered curves.

<img src="README_files/figure-gfm/plot_sim_data-1.png" style="display: block; margin: auto;" />

Continuously observed curves are shown above in order to illustrate the
misalignment problem and our simulated data; the simulated dataset also
includes binary values which have been generated by using these
continuous curves as probabilities. The unregistered and registered
binary curves for two subjects are shown below.

<img src="README_files/figure-gfm/plot_2subjs-1.png" style="display: block; margin: auto;" />

Our software registers curves by estimating *t*. For this we use the
function `registration_fpca()`.

``` r
binary_registration = register_fpca(Y = registration_data, family = "binomial", 
                                    Kt = 6, Kh = 4, npc  = 1)
## Running initial registration step
## current iteration: 1
## Running final FPCA step
```

The plot below shows unregistered, true registered, and estimated
registered binary curves for two subjects after fitting our method.

<img src="README_files/figure-gfm/plot_fit-1.png" style="display: block; margin: auto;" />

### Citation

If you like our software, please cite it in your work! To cite the
latest `CRAN` version of the package with `BibTeX`, use

    @Manual{,
        title = {registr: Registration for Exponential Family Functional Data},
        author = {Julia Wrobel and Alexander Bauer and Erin McDonnell and Jeff Goldsmith},
        year = {2022},
        note = {R package version 2.1.0},
        url = {https://CRAN.R-project.org/package=registr},
      }

To cite the 2021 Journal of Open Source Software paper, use

    @article{wrobel2021registr,
      title={registr 2.0: Incomplete Curve Registration for Exponential Family Functional Data},
      author={Wrobel, Julia and Bauer, Alexander},
      journal={Journal of Open Source Software},
      volume={6},
      number={61},
      pages={2964},
      year={2021}
    }

To cite the 2018 Journal of Open Source Software paper, use

    @article{wrobel2018regis,
      title={registr: Registration for Exponential Family Functional Data},
      author={Wrobel, Julia},
      journal={The Journal of Open Source Software},
      volume={3},
      year={2018}
    }

### Contributions

------------------------------------------------------------------------

If you find small bugs, larger issues, or have suggestions, please file
them using the [issue
tracker](https://github.com/julia-wrobel/registr/issues) or email the
maintainer at <julia.wrobel@cuanschutz.edu>. Contributions (via pull
requests or otherwise) are welcome.
# registr 0.1.0

* Added a `NEWS.md` file to track changes to the package.


# registr 1.0.0

* Preparing for first release to CRAN
* Added Erin McDonnell as an author


# registr 1.1.0

* Erin McDonnel added non-parametric updates
* Updated vignette as well
* Use of periodic b-splines
* get her to implement gradient for parametric warping stuff

# registr 2.0.0

* registr is now able to handle incomplete curves (see new vignette)
* Added 'two-step GFPCA' as alternative GFPCA approach
* New plot function for GFPCA results
* Parallelized the registration call
* Added user-specified template functions by argument 'Y_template'
* Added Gamma and Poisson family
* Added Alexander Bauer an an author

# registr 2.1.0

* Changed convergence threshold in register_fpca from 0.00001 to 0.0001 as a more reasonable threshold in many data situations.
* Improved robustness of 'cov_hall' by (i) first centering the curves before taking the cross product and smoothing it and by (ii) ensuring positive (semi-)definiteness of the covariance matrix with 'Matrix::nearPD'.
* Improved speed of 'cov_hall' for large data by using 'mgcv::bam' for smoothing.
* Added argument 'npc_varExplained' to functions 'fpca_gauss' and 'bfpca' to choose the number of FPCs based on the share of explained variation.
* All 'verbose' arguments now take values between 0-4 to better control the level of detail of info messages.
* A lot of minor fixes and refinements.
* Updated potentially invalid URLs


Preparing things below

* mgcv spline choices
* penalization for optimization (not sure how to choose best sigma)
## Test environments
* local R installation, R 3.6.2
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.


## CRAN submission 1.0

Documentation fixes were requested in an email received by Jelena Saf on March 6, 2020. Bullets below detail changes that were made to address these requests.

* Year was added to references in Description field of DESCRIPTION file
* Description field of DESCRIPTION file was lengthened to provide a full paragraph of information
* Author field of DESCRIPTION file was changed to include Erin McDonnell and Jihui Lee was removed (Jihui was accidentally referenced as an author in another file but was not actually a contributor to the package)
* Removed \dontrun{} from examples that can be checked relatively quickly
* Replaced \dontrun{} with \donttest{} for slow examples. The main user-facing functions all have examples that are tested, the \donttest{} examples just provide additional documentation for user on how to use the software for a large real-world data set.
* Have adjusted .Rd files to comply with CRAN standards
* Have added \value section to all .Rd files that are not data files to explain the function results
  


---
title: 'registr 2.0: Incomplete Curve Registration for Exponential Family Functional Data'
tags:
  - R
  - Statistical analysis
  - Functional data
  - Partially observed curves
authors:
  - name: Julia Wrobel
    orcid: 0000-0001-6783-1421
    affiliation: 1
  - name: Alexander Bauer
    orcid: 0000-0003-3495-5131
    affiliation: 2
affiliations:
 - name: Department of Biostatistics and Informatics, University of Colorado Denver, USA
   index: 1
 - name: Department of Statistics, LMU Munich, Germany
   index: 2
date: 15 December 2020
bibliography: paper.bib
---

# Introduction

Functional data are observed in many different fields.
Typical examples are longer-term panel studies where a sequence of measurements
is observed for each subject.
Compared to classical longitudinal studies, functional data analysis focuses
more on the shapes of the (time-dependent) processes by analyzing the observed
curve per subject.
E.g., one can analyze the speed of growth of children until adulthood
in the Berkeley child growth study (see left pane of \autoref{fig:registration}).

Functional data comprise different modes of variation.
In the Berkeley study, not only can growth spurts be more or less pronounced
regarding the actual growth (i.e., _amplitude variation_ along the y-axis), but each spurt
can also be shifted for some months / years for individual subjects (i.e., _phase variation_ along the x-axis).
Observed curves often have to be preprocessed with a _registration method_ in
order to separate phase and amplitude variation before analysis.

Most registration methods can only handle continuous data or data with a Gaussian
structure. However, functional data are often non-Gaussian or even categorical.
E.g., function values could be binary indicators representing physical (in)activity of patients over time [@wrobel2019].
Moreover, most registration approaches are only applicable to completely observed curves that
comprise the underlying process from its very start to its very end.

Basic routines for registering (Gaussian) data are implemented in R package
@R_fda. Performing joint registration and clustering is possible with @R_fdakma.
The popular square-root velocity function (SRVF) framework for curve registration is implemented in @R_fdasrvf
for completely observed curves on a regular grid. Similar to our approach the package
allows for registering all curves to similar shapes which can be well represented
by some low-rank basis.

# Exponential Family-based Registration

The `registr` package is based on the methods outlined in @wrobel2019.
Registration is performed using a likelihood-based approach and estimates
_inverse warping functions_ ${h_i^{-1}: t_i^* \mapsto t}$ that map the observed
time domain $t_i^*$ for subject $i$ to the common time domain $t$.
The overall model is

$$
\begin{aligned}
E\left[Y_i\left(h_i^{-1}(t_i^*)\right) | h_i^{-1}, \alpha(t), \boldsymbol{c}_i, \boldsymbol{\psi}(t) \right] &= \mu_i(t), \\
g\left[\mu_i(t)\right] &= \alpha(t) + \sum_{k = 1}^K c_{ik}\psi_k(t),
\end{aligned}
$$

with $Y_i\left(t_i^*\right)$ and $Y_i\left(h_i^{-1}(t_i^*)\right)$ the unregistered and registered curves, respectively,
and $\mu_i(t)$ the estimated subject-specific means serving as template functions, i.e., the target for the registration.
The assumed distribution with link function $g(\cdot)$ and this conditional expectation allow us to define a log-likelihood $\ell(i)$ for each observed function [see @wrobel2019].
The subject-specific means $\mu_i(t)$ are expressed through a low-rank representation based on
a population-level mean $\alpha(t)$ and a linear combination of population-level basis functions $\psi_k(t)$
and subject-specific scores $\boldsymbol{c}_i$, composed with a fixed link function $g(\cdot)$.
We estimate this representation using a likelihood-based
approach for generalized functional principal component analysis (GFPCA).

The overall model is estimated with function `register_fpca()`, which iterates 
between the estimation of warping
functions (implemented in function `registr()`)
and GFPCA estimation (functions `fpca_gauss()` or `bfpca()` for Gaussian or binomial data, respectively).
This approach is consistent with earlier versions of the `registr` package [compare @wrobel2018].

In version 2.0, the package now includes the _two-step GFPCA_ approach
of @gertheiss2017 to handle further exponential family distributions.
The respective implementation is based on the `gfpca` package of @goldsmith2016.
New distributions are supported both for registration and GFPCA.
Furthermore, for the registration step, the individual template functions (to which each curve is mapped)
can now be flexibly defined by the user with the argument `Y_template` in `registr()` and `register_fpca()`.
This is of relevance since in many settings the overall mean of the unregistered curves
is no reasonable template.

# Incomplete Curve Registration

We extend the approach of @wrobel2019 to
incomplete curves where the underlying process was either not observed
from its very beginning (i.e., _leading incompleteness_) or until its very end
(_trailing incompleteness_), or both (_full incompleteness_).

Since the underlying process is fully contained in the observed interval for complete curves, the first and last value of complete-curve warping functions lie on the diagonal line so that they preserve the overall domain.
For incomplete curves, warping functions are estimated without this
starting point and / or endpoint constraint.

However, fully removing these constraints can lead to extreme distortions
of the time domain.
We include a regularization term $\lambda$ that penalizes the amount of domain dilation
or compression performed by the inverse warping functions.
Mathematically speaking, we add a penalization term to the log likelihood $\ell(i)$
for curve $i$. For a setting with full incompleteness this results in

$$
\begin{aligned}
\ell_{\text{pen}}(i) &= \ell(i) - \lambda \cdot n_i \cdot \text{pen}(i), \\
\text{with} \ \ \ 
\text{pen}(i) &= \left( [\hat{h}_i^{-1}(t_{max,i}^*) - \hat{h}_i^{-1}(t_{min,i}^*)] - [t_{max,i}^* - t_{min,i}^*] \right)^2,
\end{aligned}
$$

where $t^*_{min,i},t^*_{max,i}$ are the minimum / maximum of the observed time domain of curve $i$ and
$\hat{h}^{-1}_i(t^*_{min,i}), \hat{h}^{-1}_i(t^*_{max,i})$ the inverse warping function evaluated at this
minimum / maximum.
For leading incompleteness with $h_i^{-1}(t_{max,i}^*) = t_{max,i}^* \forall i$ this simplifies to
$\text{pen}(i) = \left(\hat{h}_i^{-1}(t_{min,i}^*) - t_{min,i}^*\right)^2$, and for trailing incompleteness with
$h_i^{-1}(t_{min,i}^*) = t_{min,i}^* \forall i$ to
$\text{pen}(i) = \left(\hat{h}_i^{-1}(t_{max,i}^*) - t_{max,i}^*\right)^2$.
The penalization term is scaled by the number of measurements $n_i$ of curve $i$
to ensure a similar impact of the penalization for curves with different numbers
of measurements.
In practical settings, $\lambda$ has to be set manually to specify which kinds of
warpings are deemed unrealistic and should be prevented.
The choice of $\lambda$ should be based on subject knowledge by comparing
the registration results given different $\lambda$ values.

In `registr()` and `register_fpca()` the type of incompleteness can be defined
by argument `incompleteness`.
Further details are given in the package vignette _incomplete_curves_.
When applied to the Berkeley data with simulated full incompleteness,
our approach leads to a reasonable registration as shown in \autoref{fig:registration}.

![Left pane: Berkeley child growth data with simulated incompleteness; center: curves after registration; right: estimated inverse warping functions.\label{fig:registration}](figures/2_registration.png)

# Acknowledgements

We thank Fabian Scheipl and Helmut Küchenhoff for valuable methodological contributions.

# References---
title: 'registr: Registration for Exponential Family Functional Data'
authors:
- affiliation: 1
  name: Julia Wrobel
  orcid: 000-0001-6783-1421
date: '2018-01-18'
output: pdf_document
bibliography: paper.bib
tags:
- R
- Statistical analysis
- Functional data
affiliations:
- index: 1
  name: Columbia University
---

# Summary

Functional data analysis is a set of tools for understanding patterns and variability in data where the basic unit of observation is a curve measured over time, space, or another domain. Classic functional data analysis assumes that each curve is continuous or comes from a Gaussian distribution. However, applications with exponential family functional data -- curves that arise from any exponential family distribution, and have a smooth latent mean -- are increasingly common.  

Often in a functional dataset curves have similar underlying patterns but the main features of each curve, such as the minimum and maximum, have shifts such that the data appear misaligned. This misalignment can obscure patterns shared across curves and produce messy summary statistics. Registration methods reduce variability in functional data and clarify underlying patterns by aligning curves. Our method estimates a map, called a **warping function**, which transforms the domain from so that curves are aligned. The model for registration can be written

$$
E\left[Y_i\left(h_i^{-1}(t_i^*)\right) | c_i, h_i^{-1} \right] = \mu_i(t)
$$

$$
g\left[\mu_i(t)\right]= \alpha(t) + \sum_{k = 1}^K c_{ik}\psi_k(t).
$$

For subject $i$, warping function $h_i^{-1}$ maps the domain on which curves are misaligned, $t_i^*$, to aligned domain $t$ such that $h_i^{-1}(t_i^*) = t$. Then $Y_i\left(t_i^*\right)$ and $Y_i\left(h_i^{-1}(t_i^*)\right)$ are the unregistered and registered functional response curves, respectively. The $\mu_i(t)$ are subject-specific means related to the population-level mean $\alpha(t)$ and a linear combination of population-level basis functions $\psi(t)$ and subject-specific scores $c_i$ through a known link function $g$. 

The `registr` package estimates warping functions and other parameters in this model via a two-step iterative algorithm which is detailed in @wrobel2018. The main function is `register_fpca`, which registers functional data from a specified exponential family distribution. `register_fpca` reads in a long-format functional dataset and outputs an object of class `registration`.

To enhance computational efficency, key algorithm components are implemented in C++ using the R libraries `Rcpp` and `RcppArmadillo` [@rcpp, @rcppArma]. Interactive visualizations are enabled with the `refund.shiny` package [@refund.shiny, @wrobel2016]. 



# References
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# registr <img src="README_files/figures/registr.png" align="right" height = "150" />

<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version/registr)](https://cran.r-project.org/package=registr)
[![](http://cranlogs.r-pkg.org/badges/grand-total/registr?color=green)](https://cran.r-project.org/package=registr)
[![](https://travis-ci.org/julia-wrobel/registr.svg?branch=master)](https://travis-ci.org/julia-wrobel/registr)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/julia-wrobel/registr?branch=master&svg=true)](https://ci.appveyor.com/project/julia-wrobel/registr)
[![Codecov test coverage](https://codecov.io/gh/julia-wrobel/registr/branch/master/graph/badge.svg)](https://codecov.io/gh/julia-wrobel/registr/coverage.svg?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02964/status.svg)](https://doi.org/10.21105/joss.02964)
[![R-CMD-check](https://github.com/julia-wrobel/registr/workflows/R-CMD-check/badge.svg)](https://github.com/julia-wrobel/registr/actions)
<!-- badges: end -->

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE,
											warning = FALSE,
											message = FALSE,
											collapse = TRUE)

```

Registration for incomplete exponential family functional data. 

* Authors: [Julia Wrobel](http://juliawrobel.com), [Alexander Bauer](https://www.en.stablab.stat.uni-muenchen.de/people/doktoranden/bauer1/index.html), [Erin McDonnell](http://eimcdonnell.com/),
and [Jeff Goldsmith](https://jeffgoldsmith.com/)
* License: [MIT](https://opensource.org/licenses/MIT). See the [LICENSE](LICENSE) file for details
* Version: 2.1

### What it does

---------------

Functional data analysis is a set of tools for understanding patterns and variability in data where the basic unit of observation is a curve measured over some domain such as time or space. An example is an accelerometer study where intensity of physical activity was measured at each minute over 24 hours for 50 subjects. The data will contain 50 curves, where each curve is the 24-hour activity profile for a particular subject.

Classic functional data analysis assumes that each curve is continuous or comes from a Gaussian distribution. However, applications with exponential family functional data -- curves that arise from any exponential family distribution, and have a smooth latent mean -- are increasingly common. For example, take the accelerometer data just mentioned, but assume researchers are interested in *sedentary behavior* instead of *activity intensity*. At each minute over 24 hours they collect a binary measurement that indicates whether a subject was active or inactive (sedentary). Now we have a *binary curve* for each subject -- a trajectory where each time point can take on a value of 0 or 1. We assume the binary curve has a smooth latent mean, which in this case is interpreted as the probability of being active at each minute over 24 hours. This is a  example of exponential family functional data. 

Often in a functional dataset curves have similar underlying patterns but the main features of each curve, such as the minimum and maximum, have shifts such that the data appear misaligned. This misalignment can obscure patterns shared across curves and produce messy summary statistics. Registration methods reduce variability in functional data and clarify underlying patterns by aligning curves.

This package implements statistical methods for registering exponential family functional data. The basic methods are described in more detail in our [paper](http://juliawrobel.com/Downloads/registration_ef.pdf) and were further adapted to (potentially) incomplete curve settings where (some) curves
are not observed from the very beginning and/or until the very end of the common domain.
For details on the incomplete curve methodology and how to use it see the corresponding package vignette.
Instructions for installing the software and using it to register simulated binary data are provided below.

### Installation

---------------

To install from `CRAN`, please use:

```{r, eval = FALSE, echo = TRUE}
install.packages("registr")
```


To install the latest version directly from Github, please use:

```{r, eval = FALSE, echo = TRUE}
install.packages("devtools")
devtools::install_github("julia-wrobel/registr")
```


The `registr` package includes vignettes with more details on package use and functionality. To install the latest version and pull up the vignettes please use:

```{r, eval = FALSE, echo = TRUE}
devtools::install_github("julia-wrobel/registr", build_vignettes = TRUE)
vignette(package = "registr")
```


### How to use it

---------------

This example registers simulated binary data. More details on the use of the package can be found in the vignettes mentioned above. 

The code below uses `registr::simulate_unregistered_curves()` to simulate curves for 100 subjects with 200 timepoints each, observed over domain $(0, 1)$. All curves have similar structure but the location of the peak is shifted. On the observed domain $t^*$ the curves are unregistered (misaligned). On the domain $t$ the curves are registered (aligned). 


```{r simulate_data, echo = TRUE}
library(registr)

registration_data = simulate_unregistered_curves(I = 100, D = 200, seed = 2018)
```

The plot below shows the unregistered curves and registered curves.

```{r plot_sim_data, fig.align='center', fig.height=3, fig.width=9}
library(tidyverse)
library(cowplot)

unreg = ggplot(registration_data, aes(x = index, y = boot::inv.logit(latent_mean),
																			group = id)) +
	geom_path(alpha = .25) + theme_bw() + 
	labs(x = "t_star", y = "Prob(Y = 1)")


reg = ggplot(registration_data, aes(x = t, y = boot::inv.logit(latent_mean), 
																		group = id)) +
	geom_path(alpha = .25) + theme_bw() + 
	labs(x = "t", y = "Prob(Y = 1)")

cowplot::plot_grid(unreg, reg, ncol = 2)
```


Continuously observed curves are shown above in order to illustrate the misalignment problem and our simulated data; the simulated dataset also includes binary values which have been generated by using these continuous curves as probabilities. The unregistered and registered binary curves for two subjects are shown below.

```{r plot_2subjs, fig.align='center', fig.height=3, fig.width=9}
IDs = c(63, 85)
sub_data = registration_data %>% filter(id %in% IDs)

unreg = ggplot(sub_data, aes(x = index, y = boot::inv.logit(latent_mean),
														 group = id, color = factor(id))) +
	geom_path() + theme_bw() + theme(legend.position = "none") + 
	geom_point(aes(y = value), alpha = 0.25, size = 0.25) +
	labs(x = "t_star", y = "Prob(Y = 1)")


reg = ggplot(sub_data, aes(x = t, y = boot::inv.logit(latent_mean), 
													 group = id, color = factor(id))) +
	geom_path() + theme_bw() + theme(legend.position = "none") +  
	geom_point(aes(y = value), alpha = 0.25, size = 0.25) +
	labs(x = "t", y = "Prob(Y = 1)")

cowplot::plot_grid(unreg, reg, ncol = 2)
```

Our software registers curves by estimating $t$. For this we use the function `registration_fpca()`.

```{r register_data, echo = TRUE, message = TRUE}
binary_registration = register_fpca(Y = registration_data, family = "binomial", 
                                    Kt = 6, Kh = 4, npc  = 1)
```

The plot below shows unregistered, true registered, and estimated registered binary curves for two subjects after fitting our method.

```{r plot_fit, fig.align='center', fig.height=3, fig.width=9}
sub_data = binary_registration$Y %>% filter(id %in% IDs)

unreg = ggplot(sub_data, aes(x = tstar, y = boot::inv.logit(latent_mean),
														 group = id, color = factor(id))) +
	geom_path() + theme_bw() + theme(legend.position = "none") + 
	geom_point(aes(y = value), alpha = 0.25, size = 0.25) +
	labs(x = "t_star", y = "Prob(Y = 1)")


reg = ggplot(sub_data, aes(x = t, y = boot::inv.logit(latent_mean), 
													 group = id, color = factor(id))) +
	geom_path() + theme_bw() + theme(legend.position = "none") +  
	geom_point(aes(y = value), alpha = 0.25, size = 0.25) +
	labs(x = "t", y = "Prob(Y = 1)")

reg_hat = ggplot(sub_data, aes(x = t_hat, y = boot::inv.logit(latent_mean), 
													 group = id, color = factor(id))) +
	geom_path() + theme_bw() + theme(legend.position = "none") +  
	geom_point(aes(y = value), alpha = 0.25, size = 0.25) +
	labs(x = "t", y = "Prob(Y = 1)")

cowplot::plot_grid(unreg, reg, reg_hat, ncol = 3)

```

### Citation

If you like our software, please cite it in your work! To cite the latest `CRAN` version of the package with `BibTeX`, use


```{}
@Manual{,
    title = {registr: Registration for Exponential Family Functional Data},
    author = {Julia Wrobel and Alexander Bauer and Erin McDonnell and Jeff Goldsmith},
    year = {2022},
    note = {R package version 2.1.0},
    url = {https://CRAN.R-project.org/package=registr},
  }
```



To cite the 2021 Journal of Open Source Software paper, use

```{}
@article{wrobel2021registr,
  title={registr 2.0: Incomplete Curve Registration for Exponential Family Functional Data},
  author={Wrobel, Julia and Bauer, Alexander},
  journal={Journal of Open Source Software},
  volume={6},
  number={61},
  pages={2964},
  year={2021}
}
```



To cite the 2018 Journal of Open Source Software paper, use

```{}
@article{wrobel2018regis,
  title={registr: Registration for Exponential Family Functional Data},
  author={Wrobel, Julia},
  journal={The Journal of Open Source Software},
  volume={3},
  year={2018}
}
```


### Contributions

---------------

If you find small bugs, larger issues, or have suggestions, please file them using the [issue tracker](https://github.com/julia-wrobel/registr/issues) or email the maintainer at <julia.wrobel@cuanschutz.edu>. Contributions (via pull requests or otherwise) are welcome.



---
title: "Registering Incomplete Curves"
author: "Alexander Bauer"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
    number_sections: true
vignette: >
  %\VignetteIndexEntry{Registering Incomplete Curves}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

<style type="text/css">
h1 { /* Header 1 */
  font-size: 26px;
}
h2 { /* Header 2 */
  font-size: 20px;
}
h3 { /* Header 3 */
  font-size: 16px;
}
</style>


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 5
)
```

This vignette outlines the functionalities of the `registr` package with
regard to incomplete curves.

```{r load_libraries, echo = FALSE}
library(registr)
have_ggplot2 = requireNamespace("ggplot2", quietly = TRUE)
if (have_ggplot2) {
  library(ggplot2)
  theme_set(theme_minimal() + theme(plot.title = element_text(hjust = 0.5)))
}
```

# Introduction

Incomplete curves arise in many applications. Incompleteness refers to functional
data where (some) curves were not observed from the very beginning and/or until
the very end of the common domain.
Such a data structure is e.g. observed in the presence of drop-out in panel studies.

We differentiate three different types of (in)completeness:

1. **no incompleteness**,  
where processes where all observed from their very beginning until their very end.
In this case, it is reasonable to assume that both the starting points and the
endpoints of the warping functions in the registration process lie on the diagonal
since the observed process is fully comprised in the observed interval.
2. **leading incompleteness**,  
where processes where not necessarily observed from their very beginning, but
until their very end. In this case it is reasonable to assume that the endpoints
lie on the diagonal since the observed process is observed until its end.
The starting points of the warping functions are able to vary from the diagonal
to handle potential time distortions towards the beginning of the observed domains.
3. **trailing incompleteness**,  
where processes where observed from their very beginning, but not necessarily
until their very end. In this case it is reasonable to assume that the starting points
lie on the diagonal since the observed process is observed from its beginning.
The endpoints of the warping functions are able to vary from the diagonal
to handle potential time distortions towards the end of the observed domains.
4. **full incompleteness**,  
where processes where neither necessarily observed from their very beginning,
nor until their very end. In this case it is not reasonable to assume that either
the starting points or the endpoints lie on the diagonal.
The starting points and the endpoints of the warping functions are able to
vary from the diagonal to handle potential time distortions both towards the
beginning and the end of the observed domains.

Exemplarily, we showcase the following functionalities on data from the Berkeley
Growth Study (see `?growth_incomplete`) where we artificially simulated that
not all children were observed right from the start and that a relevant part of
the children dropped out early of the study at some point in time:

```{r Berkeley data}
dat = registr::growth_incomplete

# sort the data by the amount of trailing incompleteness
ids    = levels(dat$id)
dat$id = factor(dat$id, levels = ids[order(sapply(ids, function(curve_id) {
	max(dat$index[dat$id == curve_id])
}))])

if (have_ggplot2) {
  # spaghetti plot
  ggplot(dat, aes(x = index, y = value, group = id)) +
    geom_line(alpha = 0.2) +
    xlab("t* [observed]") + ylab("Derivative") +
    ggtitle("First derivative of growth curves")
}
```


```{r Berkeley data 2 lasagna, fig.height = 5.5}
if (have_ggplot2) {
  ggplot(dat, aes(x = index, y = id, col = value)) + 
    geom_line(lwd = 2.5) +
    scale_color_continuous("Derivative", high = "midnightblue", low = "lightskyblue1") +
    xlab("t* [observed]") + ylab("curve") +
    ggtitle("First derivative of growth curves") +
    theme(panel.grid  = element_blank(),
          axis.text.y = element_blank())
}
```


# Incomplete curve methodology

We adapt the registration methodology outlined in Wrobel et al. (2019) to
handle incomplete curves. Since each curve potentially has an individual range
of its observed time domain, the spline basis for estimating a curve's warping
function is defined individually for each curve, based on a given number of
basis functions.

It often is a quite strict assumption in incomplete data
settings that all warping functions start and/or end on the diagonal, i.e. that the individual,
observed part of the whole time domain is not (to some extent) distorted.
Therefore, the `registr` package gives the additional option to estimate
warping functions without the constraint that their starting point and/or endpoint
lies on the diagonal.

On the other hand, if we fully remove such constraints, this can result in
very extreme and unrealistic distortions
of the time domain. This problem is further accompanied by the fact that
the assessment of some given warping to be realistic or unrealistic can heavily
vary between different applications.
As of this reason, our method includes a penalization parameter $\lambda$ that
has to be set manually to specify which kinds of distortions are deemed realistic
in the application at hand.

Mathematically speaking, we add a penalization term to the likelihood $\ell(i)$ 
for curve $i$. For a setting with **full incompleteness** (i.e., where both the starting
point and endpoint are free to vary from the diagonal) this results in
$$
\begin{aligned}
\ell_{\text{pen}}(i) &= \ell(i) - \lambda \cdot n_i \cdot \text{pen}(i), \\
\text{with} \ \ \ 
\text{pen}(i) &= \left( \left[\hat{h}_i^{-1}(t_{max,i}^*) - \hat{h}_i^{-1}(t_{min,i}^*)\right] - \left[t_{max,i}^* - t_{min,i}^*\right] \right)^2,
\end{aligned}
$$
where $t^*_{min,i},t^*_{max,i}$ are the minimum / maximum of the observed time domain of curve $i$ and
$\hat{h}^{-1}_i(t^*_{min,i}), \hat{h}^{-1}_i(t^*_{max,i})$ the inverse warping function evaluated at this
minimum / maximum. For leading incompleteness with $h_i^{-1}(t_{max,i}^*) = t_{max,i}^* \forall i$ this simplifies to
$\text{pen}(i) = \left(\hat{h}_i^{-1}(t_{min,i}^*) - t_{min,i}^*\right)^2$, and for trailing incompleteness with
$h_i^{-1}(t_{min,i}^*) = t_{min,i}^* \forall i$ to
$\text{pen}(i) = \left(\hat{h}_i^{-1}(t_{max,i}^*) - t_{max,i}^*\right)^2$.
The penalization term is scaled by the number of measurements $n_i$ of curve $i$
to ensure a similar impact of the penalization for curves with different numbers
of measurements.

The higher the penalization parameter $\lambda$, the more the length of the registered domain
is forced towards the length of the observed domain.
Given a specific application, $\lambda$ should be chosen s.t.
unrealistic distortions of the time domain are prevented.
To do so, the user has to run the registration approach multiple times with
different $\lambda$'s to find an optimal value.


# Application on incomplete growth data

By default, both functions `register_fpca` and `registr` include the argument
`incompleteness = NULL` to constrain all warping functions to start and end on the diagonal.

```{r application 1}
reg1 = registr(Y = dat, family = "gaussian")

if (have_ggplot2) {
  ggplot(reg1$Y, aes(x = tstar, y = index, group = id)) + 
    geom_line(alpha = 0.2) +
    xlab("t* [observed]") + ylab("t [registered]") +
    ggtitle("Estimated warping functions")
}
```

```{r application 1 lasagna, fig.height = 5.5}
if (have_ggplot2) {
  ggplot(reg1$Y, aes(x = index, y = id, col = value)) + 
    geom_line(lwd = 2.5) +
    scale_color_continuous("Derivative", high = "midnightblue", low = "lightskyblue1") +
    xlab("t [registered]") + ylab("curve") +
    ggtitle("Registered curves") +
    theme(panel.grid  = element_blank(),
          axis.text.y = element_blank())
}
```

```{r application 1 spaghetti}
if (have_ggplot2) {
  ggplot(reg1$Y, aes(x = index, y = value, group = id)) +
    geom_line(alpha = 0.3) +
    xlab("t [registered]") + ylab("Derivative") +
    ggtitle("Registered curves")
}
```

The assumption can be dropped by setting `incompleteness` to some other value than NULL and
some nonnegative value for the penalization parameter `lambda_inc`.
The higher `lambda_inc` is chosen, the more the registered domains are forced to have the
same length as the observed domains.

### Small `lambda_inc`

```{r application 2}
reg2 = registr(Y = dat, family = "gaussian",
							 incompleteness = "full", lambda_inc = 0)

if (have_ggplot2) {
  ggplot(reg2$Y, aes(x = tstar, y = index, group = id)) + 
    geom_line(alpha = 0.2) +
    xlab("t* [observed]") + ylab("t [registered]") +
    ggtitle("Estimated warping functions")
}
```

```{r application 2 lasagna, fig.height = 5.5}
if (have_ggplot2) {
  ggplot(reg2$Y, aes(x = index, y = id, col = value)) + 
    geom_line(lwd = 2.5) +
    scale_color_continuous("Derivative", high = "midnightblue", low = "lightskyblue1") +
    xlab("t [registered]") + ylab("curve") +
    ggtitle("Registered curves") +
    theme(panel.grid  = element_blank(),
          axis.text.y = element_blank())
}
```

```{r application 2 spaghetti}
if (have_ggplot2) {
  ggplot(reg2$Y, aes(x = index, y = value, group = id)) +
    geom_line(alpha = 0.3) +
    xlab("t [registered]") + ylab("Derivative") +
    ggtitle("Registered curves")
}
```

### Larger `lambda_inc`

```{r application 3}
reg3 = registr(Y = dat, family = "gaussian",
							 incompleteness = "full", lambda_inc = 5)

if (have_ggplot2) {
  ggplot(reg3$Y, aes(x = tstar, y = index, group = id)) + 
    geom_line(alpha = 0.2) +
    xlab("t* [observed]") + ylab("t [registered]") +
    ggtitle("Estimated warping functions")
}
```

```{r application 3 lasagna, fig.height = 5.5}
if (have_ggplot2) {
  ggplot(reg3$Y, aes(x = index, y = id, col = value)) + 
    geom_line(lwd = 2.5) +
    scale_color_continuous("Derivative", high = "midnightblue", low = "lightskyblue1") +
    xlab("t [registered]") + ylab("curve") +
    ggtitle("Registered curves") +
    theme(panel.grid  = element_blank(),
          axis.text.y = element_blank())
}
```

```{r application 3 spaghetti}
if (have_ggplot2) {
  ggplot(reg3$Y, aes(x = index, y = value, group = id)) +
    geom_line(alpha = 0.3) +
    xlab("t [registered]") + ylab("Derivative") +
    ggtitle("Registered curves")
}
```

### Choosing an optimal `lambda_inc`

As outlined, $\lambda$ should be set to the smallest value that prevents unrealistic
distortions of the time domain. The intuition of what kinds of distortions are _unrealistic_
has to be based on subject knowledge.

In the above example we see that `lambda_inc = 0` leads to extreme compressions of
the time domain, which can clearly be viewed as unrealistic.
These compressions can for example be prevented by setting `lambda_inc = 0.025`.

```{r application 4}
reg4 = registr(Y = dat, family = "gaussian",
							 incompleteness = "full", lambda_inc = .025)

if (have_ggplot2) {
  ggplot(reg4$Y, aes(x = tstar, y = index, group = id)) + 
    geom_line(alpha = 0.2) +
    xlab("t* [observed]") + ylab("t [registered]") +
    ggtitle("Estimated warping functions")
}
```

Underlying functional principal components can be estimated jointly to the
registration by calling `register_fpca()` and visualizing them with
`registr:::plot.fpca()`.

```{r application 4 joint}
reg4_joint = register_fpca(Y = dat, family = "gaussian",
                           incompleteness = "full", lambda_inc = .025,
                           npc = 4)
```

```{r application 4 joint spaghetti}
if (have_ggplot2) {
  ggplot(reg4_joint$Y, aes(x = t_hat, y = value, group = id)) +
    geom_line(alpha = 0.3) +
    xlab("t [registered]") + ylab("Derivative") +
    ggtitle("Registered curves")
}
```

```{r application 4 joint FPC plot, fig.height=6, fig.width=7}
if (have_ggplot2) {
  plot(reg4_joint$fpca_obj)
}
```


# Constraint matrices for the optimization

Warping functions are estimated using the function `constrOptim()`.
For the estimation of the warping function for curve $i$ it uses linear inequality constraints of the form
$$
\boldsymbol{u}_i \cdot \boldsymbol{\beta}_i - \boldsymbol{c}_i \geq \boldsymbol{0},
$$
where $\boldsymbol{\beta}_i$ is the parameter vector and matrix
$\boldsymbol{u}_i$ and vector $\boldsymbol{c}_i$ define the constraints.

For the estimation of a warping function the parameter vector is constrained
s.t. the resulting warping function is monotone and does not exceed the overall
time domain $[t_{min},t_{max}]$.

In the following the constraint matrices are listed for the different settings 
of (in)completeness and assuming a parameter vector of length $p$:
$$
\boldsymbol{\beta}_i =
\left( \begin{array}{c}
\beta_{i1} \\ \beta_{i2} \\ \vdots \\ \beta_{ip}
\end{array} \right) \in \mathbb{R}_{p \times 1}
$$

**Note:**  
All following constraint matrices refer to the estimation of nonparametric inverse
warping functions with `warping = "nonparametric"`.


## Complete curve setting

When all curves were observed completely -- i.e. the underlying processes of
interest were all observed from the beginning until the end -- warping functions
can typically be assumed to start and end on the diagonal, since each process is
completely observed in its observation interval $[t^*_{min,i},t^*_{max,i}] \subset [t_{min},t_{max}]$.

Assuming that both the starting point and the endpoint lie on the diagonal,
we set $\beta_{i1} = t^*_{min,i}$ and $\beta_{ip} = t^*_{max,i}$ and only perform
the estimation for
$$
\left( \begin{array}{c}
\beta_{i2} \\ \beta_{i3} \\ \vdots \\ \beta_{i(p-1)}
\end{array} \right) \in \mathbb{R}_{(p-2) \times 1}
$$

This results in the following constraint matrices, that allow a mapping from the
observed domain $[t^*_{min,i},t^*_{max,i}]$ to the domain itself $[t^*_{min,i},t^*_{max,i}] \subset [t_{min},t_{max}]$:
$$
\begin{aligned}
\boldsymbol{u}_i &=
\left( \begin{array}{cccccccc}
1 & 0 & 0 & 0 & \ldots & 0 & 0 & 0 \\
-1 & 1 & 0 & 0 & \ldots & 0 & 0 & 0 \\
0 & -1 & 1 & 0 & \ldots & 0 & 0 & 0 \\
\vdots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots \\
0 & 0 & 0 & 0 & \ldots & 0 & -1 & 1 \\
0 & 0 & 0 & 0 & \ldots & 0 & 0 & -1
\end{array} \right) \in \mathbb{R}_{(p-1) \times (p-2)} \\
\boldsymbol{c}_i &=
\left( \begin{array}{c}
t^*_{min,i} \\ 0 \\ 0 \\ \vdots \\ 0 \\ -1 \cdot t^*_{max,i}
\end{array} \right) \in \mathbb{R}_{(p-1) \times 1}
\end{aligned}
$$

## Leading incompleteness only

In the case of *leading incompleteness* -- i.e. the underlying processes of interest
were all observed until their very end but not necessarily starting from their beginning -- warping functions
can typically be assumed to end on the diagonal, s.t. one assumes
$\beta_{ip} = t^*_{max,i}$ to let the warping functions end at the last observed
time point $t^*_{max,i}$. The estimation is then performed for the remaining
parameter vector
$$
\left( \begin{array}{c}
\beta_{i1} \\ \beta_{i3} \\ \vdots \\ \beta_{i(p-1)}
\end{array} \right) \in \mathbb{R}_{(p-1) \times 1}
$$

This results in the following constraint matrices, that allow a mapping from the
observed domain $[t^*_{min,i},t^*_{max,i}]$ to the domain $[t_{min},t^*_{max,i}] \subset [t_{min},t_{max}]$:
$$
\begin{aligned}
\boldsymbol{u}_i &=
\left( \begin{array}{cccccccc}
1 & 0 & 0 & 0 & \ldots & 0 & 0 & 0 \\
-1 & 1 & 0 & 0 & \ldots & 0 & 0 & 0 \\
0 & -1 & 1 & 0 & \ldots & 0 & 0 & 0 \\
\vdots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots \\
0 & 0 & 0 & 0 & \ldots & 0 & -1 & 1 \\
0 & 0 & 0 & 0 & \ldots & 0 & 0 & -1
\end{array} \right) \in \mathbb{R}_{p \times (p-1)} \\
\boldsymbol{c}_i &=
\left( \begin{array}{c}
t_{min} \\ 0 \\ 0 \\ \vdots \\ 0 \\ -1 \cdot t^*_{max,i}
\end{array} \right) \in \mathbb{R}_{p \times 1}
\end{aligned}
$$

## Trailing incompleteness only

In the case of *trailing incompleteness* -- i.e. the underlying processes of interest
were all observed from the beginning but not necessarily until their very end -- warping functions
can typically be assumed to start on the diagonal, s.t. one assumes
$\beta_{i1} = t^*_{min,i}$ to let the warping functions start at the first observed
time point $t^*_{min,i}$. The estimation is then performed for the remaining
parameter vector
$$
\left( \begin{array}{c}
\beta_{i2} \\ \beta_{i3} \\ \vdots \\ \beta_{ip}
\end{array} \right) \in \mathbb{R}_{(p-1) \times 1}
$$

This results in the following constraint matrices, that allow a mapping from the
observed domain $[t^*_{min,i},t^*_{max,i}]$ to the domain $[t^*_{min,i},t_{max}] \subset [t_{min},t_{max}]$:
$$
\begin{aligned}
\boldsymbol{u}_i &\text{  identical to the version for leading incompleteness} \\
\boldsymbol{c}_i &=
\left( \begin{array}{c}
t^*_{min,i} \\ 0 \\ 0 \\ \vdots \\ 0 \\ -1 \cdot t_{max}
\end{array} \right) \in \mathbb{R}_{p \times 1}
\end{aligned}
$$

## Leading and trailing incompleteness

In the case of both leading and trailing incompleteness -- i.e. the underlying
processes of interest were neither necessarily observed from their very beginnings nor to their
very ends -- warping functions can typically only be assumed to map the
observed domains $[t^*_{min,i},t^*_{max,i}]$ to the overall domain
$[t_{min},t_{max}]$.

This results in the following constraint matrices:
$$
\begin{aligned}
\boldsymbol{u}_i &=
\left( \begin{array}{cccccccc}
1 & 0 & 0 & 0 & \ldots & 0 & 0 & 0 \\
-1 & 1 & 0 & 0 & \ldots & 0 & 0 & 0 \\
0 & -1 & 1 & 0 & \ldots & 0 & 0 & 0 \\
\vdots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots & \ddots \\
0 & 0 & 0 & 0 & \ldots & 0 & -1 & 1 \\
0 & 0 & 0 & 0 & \ldots & 0 & 0 & -1
\end{array} \right) \in \mathbb{R}_{(p+1) \times p} \\
\boldsymbol{c}_i &=
\left( \begin{array}{c}
t_{min} \\ 0 \\ 0 \\ \vdots \\ 0 \\ -1 \cdot t_{max}
\end{array} \right) \in \mathbb{R}_{(p+1) \times 1}
\end{aligned}
$$


# Help files

Documentation for individual functions gives more information on their arguments and return objects, and can be pulled up via the following:

- `?register_fpca`
- `?registr`
---
title: "registr: a vignette"
author: "Julia Wrobel, Alexander Bauer and Erin McDonnell"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
    number_sections: true
bibliography: references.bib
vignette: >
  %\VignetteEncoding{UTF-8}
  
  %\VignetteIndexEntry{registr: a vignette}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

<style type="text/css">
h1 { /* Header 1 */
  font-size: 26px;
}
h2 { /* Header 2 */
  font-size: 20px;
}
h3 { /* Header 3 */
  font-size: 16px;
}
</style>


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  message   = FALSE,
  warning   = FALSE,
  fig.width = 6
)
```

The `registr` package is for registering, or aligning, exponential family functional data.

This vignette outlines the general functionality of the package.
The package can handle both complete and incomplete functional data,
i.e. curves which were not observed from the very beginning and/or until the very end
of the common domain.
Details on how to handle incomplete curves with the `registr` package can be
found in the separate vignette `"incomplete_curves"`.

```{r load_libraries, echo = FALSE}
library(registr)
library(dplyr)
have_ggplot2 = requireNamespace("ggplot2", quietly = TRUE)
have_cowplot = requireNamespace("cowplot", quietly = TRUE)
if (have_ggplot2 & have_cowplot)  {
  library(ggplot2)
  theme_set(theme_bw())
  library(cowplot)
}
```


# What is exponential family registration?

Functional data analysis is a set of tools for understanding patterns and variability in data where the basic unit of observation is a curve measured over some domain such as time or space. An example is an accelerometer study where intensity of physical activity was measured at each minute over 24 hours for 50 subjects. The data will contain 50 curves, where each curve is the 24-hour activity profile for a particular subject.

Classic functional data analysis assumes that each curve is continuous or comes from a Gaussian distribution. However, applications with exponential family functional data -- curves that arise from any exponential family distribution, and have a smooth latent mean -- are increasingly common. For example, take the accelerometer data just mentioned, but assume researchers are interested in *sedentary behavior* instead of *activity intensity*. At each minute over 24 hours they collect a binary measurement that indicates whether a subject was active or inactive (sedentary). Now we have a *binary curve* for each subject -- a trajectory where each time point can take on a value of 0 or 1. We assume the binary curve has a smooth latent mean, which in this case is interpreted as the probability of being active at each minute over 24 hours. This is a  example of exponential family functional data. 

Often in a functional dataset curves have similar underlying patterns but the main features of each curve, such as the minimum and maximum, have shifts such that the data appear misaligned. This misalignment can obscure patterns shared across curves and produce messy summary statistics. Registration methods reduce variability in functional data and clarify underlying patterns by aligning curves.

At the core of this registration method is generalized functional principal components analysis (GFPCA), a popular technique for extracting patterns shared across curves. 

## The `registr` model and algorithm 

The main model for exponential family registration is

$$
\begin{eqnarray*}
E\left[Y_i\left(h_i^{-1}(t_i^*)\right) | c_i, h_i^{-1} \right] &=& \mu_i(t) \\
g\left[\mu_i(t)\right]&=& \alpha(t) + \sum_{k = 1}^K c_{ik}\psi_k(t).
\end{eqnarray*}
$$
For subject $i$, inverse warping function $h_i^{-1}$ maps unregistered time $t_i^*$ to registered time $t$ such that $h_i^{-1}(t_i^*) = t$. $Y_i\left(t_i^*\right)$ and $Y_i\left(h_i^{-1}(t_i^*)\right)$ are the unregistered and registered response curves, respectively. The subject-specific means $\mu_i(t)$ are related to the population-level mean $\alpha(t)$ and a linear combination of population-level basis functions $\psi(t)$ and subject-specific scores $c_i$ through a known link function $g$. 

The `registr` algorithm is based on this model and iterates between the following steps:

1. Estimate subject-specific means **$\mu_i(t)$** using GFPCA, conditional on current estimate of $h_i^{-1}(t_i^*)$.
1. Estimate inverse warping functions **$h_i^{-1}(t_i^*)$**, conditional on current estimate of $\mu_i(t)$.




The methods implemented in `registr` are described in more detail in this [paper](http://juliawrobel.com/Downloads/registration_ef.pdf).


# The `registr` package

The main function in the package is `register_fpca()`. It calls two sub-functions: a GFPCA function to implement **step 1** of the iterative algorithm, and `registr()`, a function to implement **step 2** of the algorithm. The function that calculates GFPCA can either be based on the variational EM approach outlined in @wrobel_2019 or the two-step method outlined in @gertheiss_2017. In the former case, the called function depends on the family. For `family = "gaussian"` (for continuous data) and `family = "binomial"` (for binary data) the functions `bfpca()` and `fpca_gauss()` are called for the GFPCA step, respectively. In the latter case, function `gfpca_twoStep()` is called, which also supports families `"gamma"` (for strictly positive data where the variance depends on the mean) and `"poisson"` (for nonnegative count data). The `register_fpca()` function iterates between the alignment and template calculation steps until curves are registered. 

## A note on data formatting

Use of this package requires that data be in a specific format: a long-form data frame with variables `id`, `index`, and `value`, where the `value` column contains functional observations for all subjects, the `id` column identifies which observations belong to which subject, and `index` provides the grid (domain) over which the `value`s are observed.  

The variable `id` should be a unique identifier in that each id identifies a single subject. Since we assume there is only one curve per subject for this package, `id` uniquely identifies each curve as well. Other covariates can be included in the data as long as the variables `id`, `index`, and `value` are present.

# Data simulation

There are two functions for simulating data included in the package: `simulate_unregistered_curves()` and `simulate_functional_data()`. Both simulate functional data; the first is intended for demonstrating the registration algorithm and the second is for testing GFPCA sub-functions in the package.

## Simulate data for registration

`simulate_unregistered_curves()` generates curves with both unregistered and registered time grids.The code below generates data with $I = 10$ subjects and $D = 200$ using this function:

```{r sim_data2}
registration_data = simulate_unregistered_curves(I = 50, D = 200, seed = 2018)

head(registration_data)

```

The resulting object,`registration_data`, is a data frame with variables `id`, `value`, `index`, `latent_mean`, and `t`, which is consistent with the format our `registr` software requires. `id` is the identifier for a particular subject, the `value` variable contains binary observations, and `latent_mean` contains continuous observations used to generate the binary observations for the `value` variable. Note that when `family = "binomial"` we will use the binary `value` variable as the observations for each subject and when `family = "gaussian"` we use the `latent_mean` variable as the outcome.

The variables `index` and `t` are both time grids. Evaluated on the grid `index` the data is unregistered, and on the grid `t` the data is registered. Registered and unregistered curves are plotted below.

```{r plot_sim2, echo = FALSE, fig.show='hold'}
if (have_ggplot2 & have_cowplot) {
  gg1 <- registration_data %>%
    ggplot(aes(index, plogis(latent_mean), group = id)) + theme_bw() + 
    geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")
  
  gg2 <- registration_data %>%
    ggplot(aes(t, plogis(latent_mean), group = id)) + theme_bw() + 
    geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")
  
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}
```

Each curve has one main peak, but the location of that peak is shifted. When curves are registered the peaks are aligned.

## Simulate data for GFPCA

`simulate_functional_data()` simulates data with a population-level mean and two orthogonal principal components based on sine and cosine functions. The code below generates data with $I = 100$ subjects and $D = 200$ time points per subject using this function:

```{r sim_data1}
fpca_data = simulate_functional_data(I = 100, D = 200, seed = 2018)

ls(fpca_data)

head(fpca_data$Y)
```

The resulting object,`fpca_data`, is a list that contains the true population-level mean (`alpha`) and principal components (`psi1` and `psi2`), and a dataframe (`Y`). The dataframe `Y` contains variables `id`, `value`, `index` and `latent_mean`. This data is plotted below.

```{r plot1_sim1, fig.show='hold', echo = FALSE}

Y = fpca_data$Y
pc_df = data.frame(pop_mean = fpca_data$alpha, 
									 psi1 = fpca_data$psi1,
									 psi2 = fpca_data$psi2,
									 index = seq(0, 1, length.out = 200),
									 id = 1)

if (have_ggplot2 & have_cowplot) {
  gg1 <- ggplot(Y, aes(index, latent_mean, group = id)) + theme_bw() +
    geom_line(alpha = 0.25) + geom_line(data = pc_df, aes(y = pop_mean), color = "red") 
  
  gg2 <- ggplot(pc_df, aes(index, psi1)) + theme_bw() + geom_line(color = "blue") 
  gg3 <- ggplot(pc_df, aes(index, psi2)) + theme_bw() + geom_line(color = "blue") 
  
  cowplot::plot_grid(gg1, gg2, gg3, nrow = 1)
}
```

The left panel of the figure above shows the latent means for each subject, along with the population-level mean,  $\alpha(t)$, in red. The middle and right panels show the first and second principal components, $\psi_1(t)$ and $\psi_2(t)$, respectively. Using the $logit^{-1}(\cdot)$ function we can convert the subject-specific means to probabilities; these probabilities are used to generate the binary values. Binary values and latent probability curve for one subject in the dataset is shown below.

```{r plot2_sim1, echo = FALSE}
if (have_ggplot2) {
  Y %>%
    filter(id == 7) %>%
    ggplot(aes(index, value)) + theme_bw() +
    geom_point(alpha = 0.75, size = 0.25) + geom_line(aes(y = plogis(latent_mean))) +
    labs(y = "Pr(Y = 1)")
}
```

We can alter the score variance for the principal components using the arguments `lambda1` and `lambda2`. The default setting is for all subjects to have the same number of time points. However, by specifying `vary_D = TRUE`, we can generate data with uneven grid lengths for each subject.

# Joint registration and GFPCA using `register_fpca()`

`register_fpca()` is the main function for the `registr` package. Use the `family` argument to this function to specify what type of exponential family data you would like to align. The package supports `family = "gaussian"` for registering continuous data, `family = "binomial"` for registering binary data, `family = "gamma"` for strictly positive data where the variance depends on the mean and `family = "poisson"` for nonnegative count data.
The type of GFPCA is specified by the argument `fpca_type`, either calling the variational EM approach of @wrobel_2019 (`fpca_type = "variationalEM`; default) or the two-step approach of @gertheiss_2017 (`fpca_type = "two-step"`).

## Analyzing binary data

To register binary data use the following code:

```{r register_binary, message = FALSE}
registr_bin = register_fpca(Y = registration_data, family = "binomial", Kt = 8, Kh = 4, npc = 1, verbose = 2)
```

The argument `Y` specifies the input dataset; this code uses the simulated `registration_data`. `Kt` and `Kh` specify number of B-spline basis functions for the subject-specific means and warping functions, respectively, and `npc` indicates the number of functional principal components to use. The latter can also be chosen based on an explained share of variance, see argument
`npc_varExplained`.

```{r plot_reg_bin, echo = FALSE, fig.show='hold', fig.width=6}
Y = registr_bin$Y
if (have_ggplot2 & have_cowplot) {
  gg1 <- ggplot(Y, aes(tstar, plogis(latent_mean), group = id)) + theme_bw() + 
    geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")
  
  gg2 <- ggplot(Y, aes(t, plogis(latent_mean), group = id)) + theme_bw() + 
    geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")
  
  gg3 <- ggplot(Y, aes(t_hat, plogis(latent_mean), group = id)) + theme_bw() + 
    geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")
  
  cowplot::plot_grid(gg1, gg2, gg3, nrow = 1)
}

```

Underlying probabilities of the binary data are plotted above. At left probabilities on unregistered domain $t^*$, center are probabilities on true registered domain $t$, and at right are probabilities on estimated registered domain $\widehat{t}$. After registration the underlying probabilities are aligned -- though it is important to note that the algorithm registers based on the underlying binary observations, not the true probabilities.

```{r plot_reg_bin_warp, echo = FALSE, fig.show='hold'}
if (have_ggplot2 & have_cowplot) {
  gg1 <- ggplot(Y, aes(tstar, t, group = id)) + theme_bw() + 
    geom_line(alpha = 0.25)
  
  gg2 <- ggplot(Y, aes(tstar, t_hat, group = id)) + theme_bw() + 
    geom_line(alpha = 0.25)
  
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}
```

The true and estimated warping functions are plotted above. 

## Analyzing gaussian data

To register continuous data use the following code:

```{r register_gaussian, message = FALSE}
registration_data$value = registration_data$latent_mean
registr_gauss = register_fpca(Y = registration_data, family = "gaussian", npc = 1, Kt = 10)
```

# Estimating the GFPCA

Approaches for (Generalized) FPCA are implemented in functions specific for
gaussian data (`fpca_gauss()`) and binomial data (`bfpca()`).
Alternatively, `gfpca_twoStep()` allows to apply GFPCA also for other exponential
family distributions.


## Binomial FPCA using `bfpca()`

The `registr` package includes a novel variational EM algorithm for binary functional principal component analysis (bfpca), derived from methods for binary probabilistic PCA [@tipping_1999].

This `bfpca()` function works for data that is sparse and irregular (subjects do not have to be observed on the same grid and do not have to have the same number of grid points), as well as dense, regularly observed data. The following code runs bfpca on the `fpca_data` dataset.

```{r bfpca}
bfpca_object = bfpca(fpca_data$Y, npc = 2, Kt = 8, print.iter = TRUE)
```

The argument `print.iter = TRUE` prints the error after each iteration. The true and estimated population-level mean and FPCs are plotted below.

```{r plot_bfpca, echo = FALSE, fig.show='hold'}
epc_df = data.frame(index     = bfpca_object$t_vec,
                    psi1_est  = bfpca_object$efunctions[,1],
                    psi2_est  = bfpca_object$efunctions[,2],
                    alpha_est = bfpca_object$alpha %>% as.vector())
if (have_ggplot2 & have_cowplot) {
  gg1 <- ggplot() + geom_line(data = pc_df, aes(index, pop_mean), color = "blue") +
    geom_line(data = epc_df, aes(index, alpha_est), linetype = 2, color = "red") +
    theme_bw()
  
  gg2 <- ggplot() + geom_line(data = pc_df, aes(index, psi1), color = "blue") +
    geom_line(data = epc_df, aes(index, psi2_est), linetype = 2, color = "red") +
    theme_bw()
  
  gg3 <- ggplot() + geom_line(data = pc_df, aes(index, psi2), color = "blue") +
    geom_line(data = epc_df, aes(index, psi1_est), linetype = 2, color = "red") +
    theme_bw()
  
  cowplot::plot_grid(gg1, gg2, gg3, nrow = 1)
}
```

The algorithm runs quickly and does a good job recovering the true FPCs. Note that while the truth and estimation are not perfectly aligned, this is to be expected -- the data used to estimate these functions are binary observations that are generated for the truth with some variability, so results are not expected to perfectly align. One would expect results to get better with increasing number of time points per subject.

In `registr`, the estimated FPCs can be easily visualized using the internal function
`plot.fpca`. This function is automatically called when calling the general
`plot` function on an object of class `fpca`.

```{r plot.fpca, fig.show='hold', fig.width=6}
if (have_ggplot2 && requireNamespace("cowplot", quietly = TRUE)) {
  registr:::plot.fpca(bfpca_object)
}
```

## Generalized FPCA using `gfpca_twoStep()`

If `register_fpca()` is called with `fpca_type = "two-step"`, the GFPCA step
is performed with function `gfpca_twoStep()`.
As @gertheiss_2017 outline, in comparison to purely marginal and thus
biased GFPCA approaches, this two-step approach can be seen as a
*"quick-fix" for the marginal approach* of @hall_2008 that *works well in practice*.

Our implementation is based on the codebase accompanying their paper which can
be found at [github.com/jeff-goldsmith/gfpca](https://github.com/jeff-goldsmith/gfpca).
We further adapted the functions to work more efficiently both regarding the 
need for computation time and RAM, especially for large data settings with
thousands of curves.

# Estimating the registration using `registr()`

The registration step of `register_fpca()` calls the `registr` function. Though registration is intended to be performed through the `register_fpca()` function `registr()` can work as a standalone function. `registr()` uses constrained maximization of an exponential family likelihood function to estimate functions that align curves.

The default option `gradient = TRUE` uses an analytic gradient for this optimization problem (available for families `"gaussian"` and `"binomial"`). For families `"gamma"` and `"poisson"`, the gradient is computed numerically and thus less computationally efficient. The difference in computation time between `gradient = TRUE` and `gradient = FALSE` is illustrated in the code below, for `family = "binomial"`. 

```{r registr_function}
data_test_gradient = simulate_unregistered_curves(I = 50, D = 100, seed = 2018)

start_time   = Sys.time()
reg_analytic = registr(Y = data_test_gradient, family = "binomial", gradient = TRUE)
end_time     = Sys.time()

analytic_gradient = as.numeric(round((end_time - start_time), 2))

start_time  = Sys.time()
reg_numeric = registr(Y = data_test_gradient, family = "binomial", gradient = FALSE)
end_time    = Sys.time()

numeric_gradient = as.numeric(round((end_time - start_time), 2))
```

In this example with just 50 subjects and 100 time points per subject, the `registr()` function runs in `r analytic_gradient` seconds with an analytic gradient and `r numeric_gradient` seconds with a numeric gradient. Since the `register_fpca()` algorithm is iterative and calls the `registr()` function several times, using an analytic derivative drastically increases the computational efficiency, especially if the number of subjects in the data is large.

```{r registr_function gradient large data, include=FALSE}
data_test_gradient = simulate_unregistered_curves(I = 1000, D = 500, seed = 2018)

start_time = Sys.time()
reg_analytic = registr(Y = data_test_gradient, family = "binomial", gradient = TRUE)
end_time = Sys.time()

analytic_gradient_large = as.numeric(round((end_time - start_time), 1))

start_time = Sys.time()
reg_numeric = registr(Y = data_test_gradient, family = "binomial", gradient = FALSE)
end_time = Sys.time()

numeric_gradient_large = as.numeric(round((end_time - start_time), 1))
```

Running the above example with 1000 subjects and 500 time points yields computation times of `r analytic_gradient_large` seconds (for the analytic gradient) and `r numeric_gradient_large` (for the numeric gradient).

The registration step can further be parallelized by using the `cores` argument.

# Additional features

## Registering incomplete curves

Incomplete curves arise in many applications. Incompleteness refers to functional data
where (some) curves were not observed from the very beginning and/or until the very end
of the common domain.
Such a data structure is e.g. observed in the presence of drop-out in panel studies.

The `registr` package offers the possibility to flexibly account for different
types of incompleteness structures in the registration using `registr()`
as well as in the joint approach using `register_fpca()`.
All incomplete curve functionalities are outlined in the separate
vignette `"incomplete_curves"`.


## Parametric inverse warping functions

The `registr` package currently supports two types of inverse warping functions: nonparmetric B-spline basis functions (default), or parametric 2-knot piecewise linear functions. With `warping = "piecewise_linear2"`, the registration step estimates the $x$ and $y$ (or $t_i^*$ and $t$) coordinates of each of the two knots to construct an inverse warping function that consists of 3 line segments.

To register data with parametric inverse warping functions, using the following code:

```{r register_parametric}
registration_data = simulate_unregistered_curves(I = 10, D = 50, seed = 2018)

registr_parametric = register_fpca(Y = registration_data, family = "binomial", 
                                   Kt = 8, Kh = 4, npc = 1, gradient = FALSE,
                                   warping = "piecewise_linear2")
```

This argument works for all families. Note that the `gradient` option is currently unavailable for parametric inverse warping functions. Below are the resulting inverse warping functions from both the `nonparametric` (left) and `piecewise_linear2` (right) specifications. The slopes of the 3 line segments that make up one's parametric function can provide some interpretation about how that particular subject's data was warped to align with the population mean.

```{r register_parametric_plots, echo = FALSE, fig.show='hold', eval = have_ggplot2}
if (have_ggplot2 & have_cowplot) {
  gg1 <- ggplot(registr_gauss$Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line() + 
    labs(title = "warping = nonparametric")
  
  gg2 <- ggplot(registr_parametric$Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line() + 
    labs(title = "warping = piecewise_linear2")
  
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}
```

Beyond interpretability, another advantage of parametric inverse warping functions is the ability to specify prior information about how they should look. The `register_fpca()` function can include normally-distributed priors which pull warping functions toward the identity line. Specifically, the priors pull the knot coordinates toward (0.25, 0.25) and (0.75, 0.75). 

To activate the priors, one must specify `priors = TRUE` and choose a value for the argument `prior_sd`, the standard deviation that will be applied to all 4 prior distributions.

```{r register_par_priors}
registr_par_priors = register_fpca(Y = registration_data, family = "binomial", 
                                   Kt = 8, Kh = 4, npc = 1, gradient = FALSE,
                                   warping = "piecewise_linear2",
                                   priors = TRUE, prior_sd = 0.1)
```
As expected, a smaller variance lead to stronger tendency toward the prior means. This is demonstrated in the warping function plots below.

```{r, register_par_priors_plots, echo = FALSE, fig.show='hold'}
registr_par_priors2 = register_fpca(Y = registration_data, family = "binomial", 
                                    Kt = 8, Kh = 4, npc = 1, gradient = FALSE,
																		warping = "piecewise_linear2",
                                    priors = TRUE, prior_sd = 0.01)

if (have_ggplot2 & have_cowplot) {
  gg1 <- ggplot(registr_par_priors$Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line() + 
    labs(title = "sd for all priors = 0.1")
  
  gg2 <- ggplot(registr_par_priors2$Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line() + 
    labs(title = "sd for all priors = 0.01")
  
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}
```

The ability to specify priors is currently only available for binary curve registration, not Gaussian. 

## Periodic B-spline basis functions for GFPCA

In some cases it may be of interest to place periodic boundary conditions on the B-spline basis functions for the population-level and subject-specific mean templates. The periodic conditions ensure that the resulting function starts and ends with the same value. This may be useful when modeling cyclical data, such as daily physical activity patterns. To use periodic B-spline basis functions during registration, use the option `periodic = TRUE`:

```{r fpca_periodic}
registr_periodic = register_fpca(Y = registration_data, family = "binomial", 
                                 Kt = 8, Kh = 4, npc = 1, gradient = FALSE,
                                 periodic = TRUE)
```

Note that the `gradient` option is currently unavailable for `periodic = TRUE`. 

The resulting population mean and principal component from registering the continuous data using `periodic = FALSE` (default) and `periodic = TRUE` are plotted below, with dotted lines to demonstrate the unification of the start and end of the periodic functions.

```{r register_fpca_periodic, echo = FALSE, fig.show='hold'}

registr_non_periodic = register_fpca(Y = registration_data, family = "binomial", 
                                 Kt = 8, Kh = 4, npc = 1, gradient = FALSE,
                                 periodic = FALSE)
if (have_ggplot2 & have_cowplot) {
  gg1 <- tibble(mu = registr_non_periodic$fpca_obj$mu) %>%
    mutate(time = row_number()) %>%
    ggplot(aes(x = time, y = mu)) + 
    theme_bw() + 
    geom_line() + 
    labs(y = "mu (non-periodic)")
  
  gg2 <- tibble(psi1 = registr_non_periodic$fpca_obj$efunctions[,1]) %>%
    mutate(time = row_number()) %>%
    ggplot(aes(x = time, y = psi1)) + 
    theme_bw() + 
    geom_line() + 
    labs(y = "psi1 (non-periodic)")
  
  gg3 <- tibble(mu = registr_periodic$fpca_obj$mu) %>%
    mutate(time = row_number()) %>%
    ggplot(aes(x = time, y = mu)) + 
    theme_bw() + 
    geom_line() + 
    geom_hline(yintercept = registr_periodic$fpca_obj$mu[1], lty = "dotted") + 
    labs(y = "mu (periodic)")
  
  gg4 <- tibble(psi1 = registr_periodic$fpca_obj$efunctions[,1]) %>%
    mutate(time = row_number()) %>%
    ggplot(aes(x = time, y = psi1)) + 
    theme_bw() + 
    geom_line() + 
    geom_hline(yintercept = registr_periodic$fpca_obj$efunctions[1,1], lty = "dotted") + 
    labs(y = "psi1 (periodic)")
  
  cowplot::plot_grid(gg1, gg2, gg3, gg4, nrow = 2)
}
```


## Choosing the template function

By default `registr()` estimates the overall mean of all curves and uses this
mean curve as the template function to which all curves are registered.
In the joint approach (`register_fpca()`) this mean curve is used as template
for an initial registration step, before the main iteration between registration
and GFPCA starts.

In some situations the overall mean is not the most reasonable choice for the
template. See for example the following curves with a *trailing incompleteness*
structure where some processes were not observed until their very end.
The shape of the mean curve does not match any of the observed
curve shapes well.

```{r template data sim, echo=F, fig.width=5}
t = seq(0, 1, length.out = 100)
temp_dat = data.frame(index = rep(t, times = 3),
                      value = c(dnorm(t, mean = 0.5, sd = 0.15),
                                dnorm(t, mean = 0.65, sd = 0.185),
                                dnorm(t, mean = 0.7, sd = 0.18)),
                      id    = factor(rep(1:3, each = length(t))))
if (have_ggplot2) {
  ggplot(temp_dat, aes(x = index, y = value)) + 
    geom_line(aes(col = id)) +
    geom_smooth(se = FALSE, col = "black") +
    ggtitle("Simulated data with mean curve in black")
}
```

Performing the registration with this mean curve as template would lead to the
following result:

```{r temp registration without template, fig.width=5}
reg1 = registr(Y = temp_dat, family = "gaussian", Kh = 4,
               incompleteness = "trailing", lambda_inc = 0)
if (have_ggplot2) {
  ggplot(reg1$Y, aes(x = index, y = value, col = id)) + 
    geom_line() +
    ggtitle("Registration with overall mean (black) as template function")
}
```

Alternatively, the template function can be manually defined using one of the
following options:

1. Define the template as the mean curve only based on some subset of the
observed curves `Y`
2. Define the template as one observed curve
3. Define the template curve independently from the observed curves

The subset of curves (for option 1) or one specific curve (for options 2 and 3)
can be specified using the argument `Y_template`.
In the following example, the red curve with id 1 is used as template:

```{r temp registration with template, fig.width=5}
Y_template = temp_dat %>% filter(id == 1)
reg2 = registr(Y = temp_dat, family = "gaussian", Kh = 4, Y_template = Y_template,
               incompleteness = "trailing", lambda_inc = 0)
if (have_ggplot2) {
  ggplot(reg2$Y, aes(x = index, y = value, col = id)) +
    geom_line() +
    ggtitle("Registration with red curve as template")
}
```


# Help files

Documentation for individual functions gives more information on their arguments and return objects, and can be pulled up via the following:

- `?register_fpca`
- `?registr`
- `?bfpca`
- `?fpca_gauss`
- `?gfpca_twoStep`

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.fpca.R
\name{plot.fpca}
\alias{plot.fpca}
\title{Plot the results of a functional PCA}
\usage{
\method{plot}{fpca}(
  x,
  plot_FPCs = 1:x$npc,
  sd_factor = 2,
  response_function = NULL,
  add_symbols = TRUE,
  subtitle = TRUE,
  xlim = NULL,
  ylim = NULL,
  xlab = "t [registered]",
  ylab = "y",
  ...
)
}
\arguments{
\item{x}{Object of class \code{"fpca"}.}

\item{plot_FPCs}{Optional index vector of the FPCs to be plotted.
Defaults to all FPCs contained in \code{x}.}

\item{sd_factor}{Numeric factor with which the standard deviations of each
FPC's scores are multiplied to display its variation in the plots.
Defaults to 2.}

\item{response_function}{Optional response function to be applied before
plotting the curves. Defaults to \code{NULL}, i.e. the identity function if
\code{x$family} is one of \code{c("gaussian","binomial")} or
\code{exp()} if \code{x$family} is one of \code{c("gamma","poisson")}.}

\item{add_symbols}{Indicator if '+' and '-' symbols should be added to the
plot to highlight the direction of the displayed FPCs. Defaults to TRUE.}

\item{subtitle}{If TRUE (default) the parameter \code{sd_factor}
is displayed in the plot subtitle.}

\item{xlim, ylim}{Optional numeric vectors with limits for the x and y axis.}

\item{xlab, ylab}{Optional titles for the x and y axis.}

\item{...}{Additional arguments passed to \code{\link[ggplot2]{theme}}.}
}
\value{
@return If multiple FPCs are plotted, returns a grid of \code{ggplot}
plots, created with \code{cowplot::plot_grid}. If only one FPC is plotted,
returns a single \code{ggplot} plot.
}
\description{
S3 plot method for class \code{fpca}.
Plot FPCA results by visualizing the variation of the individual FPCs around
the global mean. based on an object created with function
\code{\link{fpca_gauss}}, \code{\link{bfpca}} or \code{\link{gfpca_twoStep}}. \cr \cr
The shares of explained variance are included in the plot titles if
\code{x} contains an element \code{evalues_sum}.
}
\examples{
data(growth_incomplete)

fpca_obj = fpca_gauss(Y = growth_incomplete, npc = 2)
if (requireNamespace("ggplot2", quietly = TRUE) &&
requireNamespace("cowplot", quietly = TRUE)) {
library(ggplot2)
plot(fpca_obj)
}

}
\author{
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{squareTheta}
\alias{squareTheta}
\title{Calculate quadratic form of spline basis functions for the current subject.}
\usage{
squareTheta(xi, theta)
}
\arguments{
\item{xi}{vector of variational parameters for the current subject.}

\item{theta}{spline basis functions for the current subject.}
}
\value{
A matrix of the quadratic form of theta for the current subject.
}
\description{
Calculations quadratic form of theta with diagonalized variational parameter in the center.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_helpers.R
\name{coarsen_index}
\alias{coarsen_index}
\title{Coarsen an index vector to a given resolution}
\usage{
coarsen_index(index, significant_digits)
}
\arguments{
\item{index}{Numeric vector of index values.}

\item{significant_digits}{Positive integer value.}
}
\value{
Numeric vector of rounded index values.
}
\description{
Reduce the resolution of a numeric vector by specifying the number of
\code{significant_digits} to which the numbers should be rounded. \cr \cr
Internal function used to coarsen the index vector before estimating the
two-step GFPCA with \code{\link{gfpca_twoStep}}.
}
\examples{
index_vector = c(0.7892, 0.2984, 0.328)
registr:::coarsen_index(index_vector, 1)
registr:::coarsen_index(index_vector, 3)

index_vector2 = c(2803, -7639, 13)
registr:::coarsen_index(index_vector2, 1)
registr:::coarsen_index(index_vector2, 3)

}
\author{
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/registr.R
\name{registr_oneCurve}
\alias{registr_oneCurve}
\title{Internal function to register one curve}
\usage{
registr_oneCurve(
  obj = NULL,
  Y = NULL,
  Kt = 8,
  Kh = 4,
  family = "gaussian",
  gradient = TRUE,
  incompleteness = NULL,
  lambda_inc = NULL,
  beta = NULL,
  t_min = NULL,
  t_max = NULL,
  periodic = FALSE,
  warping = "nonparametric",
  gamma_scales = NULL,
  global_knots = NULL,
  mean_coefs = NULL,
  ...,
  verbose = 1,
  just_return_list = FALSE
)
}
\arguments{
\item{obj}{Current estimate of FPC object. 
Can be NULL only if Y argument is selected.}

\item{Y}{Dataframe. Should have values id, value, index.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions. Default is 8.}

\item{Kh}{Number of B-spline basis functions used to estimate warping functions \emph{h}. Default is 4.}

\item{family}{One of \code{c("gaussian","binomial","gamma","poisson")}. Defaults to
\code{"gaussian"}.}

\item{gradient}{If \code{TRUE}, uses analytic gradient to calculate derivative. 
If \code{FALSE}, calculates gradient numerically. Not available for families
\code{"gamma","poisson"}.}

\item{incompleteness}{Optional specification of incompleteness structure.
One of \code{c("leading","trailing","full")}, specifying that incompleteness
is present only in the initial measurements, only in the trailing measurements, or
in both, respectively. For details see the accompanying vignette.
Defaults to NULL, i.e. no incompleteness structure.
Can only be set when \code{warping = "nonparametric"}.}

\item{lambda_inc}{Penalization parameter to control the amount of
overall dilation of the domain.
The higher this lambda, the more the registered domains are forced to have the
same length as the observed domains.
Only used if \code{incompleteness} is not NULL.}

\item{beta}{Current estimates for beta for each subject. Default is NULL.}

\item{t_min}{Minimum value to be evaluated on the time domain.
if `NULL`, taken to be minimum observed value.}

\item{t_max}{Maximum value to be evaluated on the time domain. 
if `NULL`, taken to be maximum observed value.}

\item{periodic}{If \code{TRUE}, uses periodic b-spline basis functions. Default is \code{FALSE}.}

\item{warping}{If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.}

\item{gamma_scales}{Only used for \code{family = "gamma"}.
Vector with one entry for each subject, containing the current estimate for the scale parameter of its
gamma distribution. Default is NULL, which sets the starting value for the scale parameter to 1.5.}

\item{global_knots}{knots for the basis/splines, passed to [pbs::pbs()] 
or [stats::bs()]}

\item{mean_coefs}{Mean coefficients for the mean of all curves or 
GFPCA based.  May extract from `obj` object}

\item{...}{additional arguments passed to or from other functions}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{just_return_list}{Do not use.  For developers only}
}
\value{
An list containing:
\item{hinv_innerKnots}{Inner knots for setting up the spline basis
for the inverse warping function.}
\item{hinv_beta}{Estimated B-spline basis coefficients used to construct
subject-specific inverse warping functions.}
\item{t_hat}{Vector of registered time domain.}
\item{loss}{Loss of the optimal solution.}
}
\description{
This internal function is only to be used from within \code{registr}.
It performs the main optimization step with \code{constrOptim} for the
registration of one curve.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Erin McDonnell \email{eim2117@cumc.columbia.edu},
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/register_fpca.R
\name{register_fpca}
\alias{register_fpca}
\title{Register curves using constrained optimization and GFPCA}
\usage{
register_fpca(
  Y,
  Kt = 8,
  Kh = 4,
  family = "gaussian",
  incompleteness = NULL,
  lambda_inc = NULL,
  Y_template = NULL,
  max_iterations = 10,
  npc = NULL,
  npc_criterion = NULL,
  fpca_type = "variationalEM",
  fpca_maxiter = 50,
  fpca_seed = 1988,
  fpca_error_thresh = 1e-04,
  fpca_index_significantDigits = 4L,
  cores = 1L,
  verbose = 1,
  ...
)
}
\arguments{
\item{Y}{Dataframe. Should have values id, value, index.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8. If
\code{fpca_type = "variationalEM"} and \code{npc_criterion} is used,
\code{Kt} is set to 20.}

\item{Kh}{Number of B-spline basis functions used to estimate warping functions \emph{h}. Default is 4.}

\item{family}{One of \code{c("gaussian","binomial","gamma","poisson")}.
Families \code{"gamma"} and \code{"poisson"} are only supported by
\code{fpca_type = "two-step"}. Defaults to \code{"gaussian"}.}

\item{incompleteness}{Optional specification of incompleteness structure.
One of \code{c("leading","trailing","full")}, specifying that incompleteness
is present only in the initial measurements, only in the trailing measurements, or
in both, respectively. For details see the accompanying vignette.
Defaults to NULL, i.e. no incompleteness structure.
Can only be set when \code{warping = "nonparametric"}.}

\item{lambda_inc}{Penalization parameter to control the amount of
overall dilation of the domain.
The higher this lambda, the more the registered domains are forced to have the
same length as the observed domains.
Only used if \code{incompleteness} is not NULL.}

\item{Y_template}{Optional dataframe with the same structure as \code{Y}.
Only used for the initial registration step. If NULL,
curves are registered to the overall mean of all curves in \code{Y} as template function.
If specified, the template function is taken as the mean
of all curves in \code{Y_template}. Defaults to NULL.}

\item{max_iterations}{Number of iterations for overall algorithm. Defaults to 10.}

\item{npc, npc_criterion}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_criterion} has to be set
to a number between 0 and 1. For \code{fpca_type = "two-step"}, it is also
possible to cut off potential tails of subordinate FPCs (see
\code{\link{gfpca_twoStep}} for details).}

\item{fpca_type}{One of \code{c("variationalEM","two-step")}.
Defaults to \code{"variationalEM"}.}

\item{fpca_maxiter}{Only used if \code{fpca_type = "variationalEM"}. Number
to pass to the \code{maxiter} argument of `bfpca()` or `fpca_gauss()`. 
Defaults to 50.}

\item{fpca_seed}{Only used if \code{fpca_type = "variationalEM"}. Number to
pass to the \code{seed} argument of `bfpca()` or `fpca_gauss()`. Defaults to
1988.}

\item{fpca_error_thresh}{Only used if \code{fpca_type = "variationalEM"}.
Number to pass to the \code{error_thresh} argument of `bfpca()` or
`fpca_gauss()`. Defaults to 0.0001.}

\item{fpca_index_significantDigits}{Only used if \code{fpca_type = "two-step"}.
Positive integer \code{>= 2}, stating the number of significant digits to which
the index grid should be rounded in the GFPCA step. Coarsening the index grid
is necessary since otherwise the covariance surface matrix explodes in size
in the presence of too many unique index values (which is the case after some
registration step). Defaults to 4. Set to \code{NULL} to prevent rounding.}

\item{cores}{Number of cores to be used. If \code{cores > 1}, the registration
call is parallelized by using \code{parallel::mclapply} (for Unix-based
systems) or \code{parallel::parLapply} (for Windows). Defaults to 1,
no parallelized call.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{...}{Additional arguments passed to registr and to the gfpca functions
(if \code{fpca_type = "variationalEM"}).}
}
\value{
An object of class \code{registration} containing:
\item{Y}{The observed data plus variables \code{t_star} and \code{t_hat} which are the
unregistered grid and registered grid, respectively.}
\item{fpca_obj}{List of items from FPCA step.}
\item{family}{Used exponential family.}
\item{index_warped}{List of the (warped) index values for each iteration.
Has \code{'convergence$iterations + 2'} elements since the first two elements
contain the original (observed) index and the warped index values from the
preprocessing registration step (see Details), respectively.}
\item{hinv_innerKnots}{List of inner knots for setting up the spline bases
for the inverse warping functions. Only contains \code{NULL} values for
\code{Kh <= 4}.}
\item{hinv_beta}{Matrix of B-spline basis coefficients used to construct the
subject-specific inverse warping functions. From the last performed
registration step. For details see \code{?registr}.}
\item{convergence}{List with information on the convergence of the joint
approach. Containing the following elements: \cr \cr
\emph{converged} \cr
Indicator if the joint algorithm converged or if not
(i.e., \code{max_iterations} was reached) \cr \cr
\emph{iterations} \cr
Number of joint iterations that were performed. \cr \cr
\emph{delta_index} \cr
Vector of mean squared differences between the (warped) index values
(scaled to [0,1] based on the size of the observed domain)
in the current and the previous iteration.
Convergence is reached if this measure drops below 0.0001. \cr \cr
\emph{registration_loss} \cr
Vector of the loss in each iteration of the algorithm.
Calculated in the registration step using the exponential family
likelihood with natural parameter from the FPCA step.
Has \code{'iterations + 1'} elements since the first element contains the
loss of the preprocessing registration step (see Details).
}
}
\description{
Function combines constrained optimization and GFPCA to estimate warping functions for 
exponential family curves. See argument \code{family} for which families are
supported. Warping functions are calculated by the function \code{\link{registr}}.
The GFPCA step can be performed either using the variational EM-based GFPCA
approaches of Wrobel et al. (2019) (\code{fpca_type = "variationalEM"}, default)
or the mixed model-based two-step approach of Gertheiss et al. (2017)
(\code{fpca_type = "two-step"}). \cr \cr
Warping functions by default are forced to start and end on the diagonal to be
domain-preserving. This behavior can be changed by setting
\code{incompleteness} to some other value than NULL and a reasonable \code{lambda_inc} value.
For further details see the accompanying vignette. \cr \cr
The number of functional principal components (FPCs) can either be specified
directly (argument \code{npc}) or chosen based on the explained share of
variance in each iteration (argument \code{npc_criterion}). \cr \cr
By specifying \code{cores > 1} the registration call can be parallelized.
}
\details{
Requires input data \code{Y} to be a dataframe in long format with variables 
\code{id}, \code{index}, and \code{value} to indicate subject IDs, 
observation times on the domain, and observations, respectively.

One joint iteration consists of a GFPCA step and a registration step.
As preprocessing, one initial registration step is performed.
The template function for this registration step is defined by argument
\code{Y_template}.
After convergence or \code{max_iterations} is reached, one final GFPCA step
is performed.
}
\examples{

### complete binomial curves
Y = simulate_unregistered_curves(I = 20, D = 200)

# estimation based on Wrobel et al. (2019)
reg = register_fpca(Y, npc = 2, family = "binomial",
                    fpca_type = "variationalEM", max_iterations = 5)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  ggplot(reg$Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line(alpha = 0.2) + ggtitle("Estimated warping functions")
  
  plot(reg$fpca_obj, response_function = function(x) { 1 / (1 + exp(-x)) })
}


\donttest{

# estimation based on Gertheiss et al. (2017)
reg2 = register_fpca(Y, npc = 2, family = "binomial",
                     fpca_type = "two-step", max_iterations = 5,
                     fpca_index_significantDigits = 4)
                     
# example using accelerometer data from nhanes 2003-2004 study
data(nhanes)
nhanes_short = nhanes[nhanes$id \%in\% unique(nhanes$id)[1:5],]
reg_nhanes   = register_fpca(nhanes_short, npc = 2, family = "binomial", max_iterations = 5)


### incomplete Gaussian curves
data(growth_incomplete)

# Force the warping functions to start and end on the diagonal
reg2a = register_fpca(growth_incomplete, npc = 2, family = "gaussian",
                      incompleteness = NULL, max_iterations = 5)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  
  ggplot(reg2a$Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Estimated warping functions")
  ggplot(reg2a$Y, aes(x = t_hat, y = value, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Registered curves")
}
# Allow the warping functions to not start / end on the diagonal.
# The higher lambda_inc, the more the starting points and endpoints are forced
# towards the diagonal.
reg2b = register_fpca(growth_incomplete, npc = 2, family = "gaussian",
                      incompleteness = "full", lambda_inc = 0.1,
                      max_iterations = 5)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot(reg2b$Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Estimated warping functions")
  ggplot(reg2b$Y, aes(x = t_hat, y = value, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Registered curves")
}

### complete Gamma curves
Y             = simulate_unregistered_curves(I = 20, D = 100)
Y$value       = exp(Y$latent_mean)
registr_gamma = register_fpca(Y, npc = 2, family = "gamma", fpca_type = "two-step",
                              gradient = FALSE, max_iterations = 3)
}

}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu}
Jeff Goldsmith \email{ajg2202@cumc.columbia.edu},
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_unregistered_curves.R
\name{grid_subj_create}
\alias{grid_subj_create}
\title{Generate subject-specific grid (t_star)}
\usage{
grid_subj_create(coefs, D)
}
\arguments{
\item{coefs}{Spline basis coefficients for reconstructing the subject-specific grid.}

\item{D}{Number of grid points per subject.}
}
\value{
A numeric vector.
}
\description{
This function creates subject-specific time grid
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_unregistered_curves.R
\name{mean_curve}
\alias{mean_curve}
\title{Simulate mean curve}
\usage{
mean_curve(grid, period = 2 * pi, spline_based = FALSE)
}
\arguments{
\item{grid}{Grid of x values over which to evaluate the function.}

\item{period}{Controls the period of the mean curve}

\item{spline_based}{If FALSE curve is constructed using sine and cosine functions,
if TRUE, curve is constructed using B-spline basis.}
}
\value{
A numeric vector.
}
\description{
This function generates mean for simulated accelerometer data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpca_gauss.R
\name{fpca_gauss}
\alias{fpca_gauss}
\title{Functional principal components analysis via variational EM}
\usage{
fpca_gauss(
  Y,
  npc = NULL,
  npc_varExplained = NULL,
  Kt = 8,
  maxiter = 20,
  t_min = NULL,
  t_max = NULL,
  print.iter = FALSE,
  row_obj = NULL,
  seed = 1988,
  periodic = FALSE,
  error_thresh = 1e-04,
  subsample = TRUE,
  verbose = 1,
  ...
)
}
\arguments{
\item{Y}{Dataframe. Should have variables id, value, index.}

\item{npc, npc_varExplained}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_varExplained} has to be set
to a share between 0 and 1.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8. If \code{npc_varExplained}
is used, \code{Kt} is set to 20.}

\item{maxiter}{Maximum number of iterations to perform for EM algorithm. Default is 50.}

\item{t_min}{Minimum value to be evaluated on the time domain.}

\item{t_max}{Maximum value to be evaluated on the time domain.}

\item{print.iter}{Prints current error and iteration}

\item{row_obj}{If NULL, the function cleans the data and calculates row indices. 
Keep this NULL if you are using standalone \code{register} function.}

\item{seed}{Set seed for reproducibility. Defaults to 1988.}

\item{periodic}{If TRUE, uses periodic b-spline basis functions. Default is FALSE.}

\item{error_thresh}{Error threshold to end iterations. Defaults to 0.0001.}

\item{subsample}{if the number of rows of the data is greater than 
10 million rows, the `id` values are subsampled to get the mean coefficients.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{...}{Additional arguments passed to or from other functions}
}
\value{
An object of class \code{fpca} containing:
\item{fpca_type}{Information that FPCA was performed with the 'variationEM' approach,
in contrast to registr::gfpca_twoStep.}
\item{t_vec}{Time vector over which the mean \code{mu} and the functional principal
components \code{efunctions} were evaluated.}
\item{knots}{Cutpoints for B-spline basis used to rebuild \code{alpha}.}
\item{efunctions}{\eqn{D \times npc} matrix of estimated FPC basis functions.}
\item{evalues}{Estimated variance of the FPC scores.}
\item{evalues_sum}{Approximation of the overall variance in \code{Y}, based
on an initial run of the FPCA with \code{npc = 20}. Is \code{NULL} if
\code{npc_varExplained} was not specified.}
\item{npc}{number of FPCs.}
\item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
\item{alpha}{Estimated population-level mean.}
\item{mu}{Estimated population-level mean. Same value as \code{alpha} but included for compatibility
with \code{refund.shiny} package.}
\item{subject_coefs}{B-spline basis coefficients used to construct subject-specific means. 
For use in \code{registr()} function.}
\item{Yhat}{FPC approximation of subject-specific means.}
\item{Y}{The observed data.}
\item{family}{\code{gaussian}, for compatibility with \code{refund.shiny} package.}
\item{sigma2}{Estimated error variance}
}
\description{
Function used in the FPCA step for registering functional data,
called by \code{\link{register_fpca}} when \code{family = "gaussian"}. 
Parameters estimated based on probabilistic PCA framework originally 
introduced by Tipping and Bishop in 1999. \cr \cr
The number of functional principal components (FPCs) can either be specified
directly (argument \code{npc}) or chosen based on the explained share of
variance (\code{npc_varExplained}). In the latter case, the explained share of
variance and accordingly the number of FPCs is estimated before the main
estimation step by once running the FPCA with \code{npc = 20} (and
correspondingly \code{Kt = 20}). Doing so, we approximate the overall
variance in the data \code{Y} with the variance represented by the FPC basis
with 20 FPCs.
}
\examples{
data(growth_incomplete)

# estimate 2 FPCs
fpca_obj = fpca_gauss(Y = growth_incomplete, npc = 2)
plot(fpca_obj)

# estimate npc adaptively, to explain 90\% of the overall variation
\donttest{
fpca_obj2 = fpca_gauss(Y = growth_incomplete, npc_varExplained = 0.9)
plot(fpca_obj, plot_FPCs = 1:2)
}

}
\references{
Tipping, M. E. and Bishop, C (1999). Probabilistic Principal Component Analysis.
\emph{Journal of the Royal Statistical Society Series B,}, 592--598.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Jeff Goldsmith \email{ajg2202@cumc.columbia.edu},
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/registr-utils.R
\name{piecewise_linear2_hinv}
\alias{piecewise_linear2_hinv}
\title{Create two-parameter piecewise linear (inverse) warping functions}
\usage{
piecewise_linear2_hinv(grid, knot_locations = c(0.25, 0.3, 0.75, 0.9))
}
\arguments{
\item{grid}{grid of values over which to evaluate the function.}

\item{knot_locations}{controls the x and y locations of the two knots.}
}
\description{
This function uses a 2-knot piecewise linear model to calculate inverse warping 
functions for registration. The parameters \code{knot1_x} and \code{knot1_y}
control the x and y locations of the first knot, and the parameters
\code{knot1_x} and \code{knot1_y} control the x and y locations of the second
knot. The designation (inverse) is intended to communicate that these 
functions take data from the unregistered space to the registered space, 
consistent with functional data literature on registration.
}
\author{
Erin McDonnell \email{eim2117@cumc.columbia.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss_h_gradient.R
\name{loss_h_gradient}
\alias{loss_h_gradient}
\title{Gradient of loss function for registration step}
\usage{
loss_h_gradient(
  Y,
  Theta_h,
  mean_coefs,
  knots,
  beta.inner,
  family = "gaussian",
  incompleteness = NULL,
  lambda_inc = NULL,
  t_min,
  t_max,
  t_min_curve,
  t_max_curve,
  Kt = 8,
  periodic = FALSE,
  warping = "nonparametric"
)
}
\arguments{
\item{Y}{vector of observed points.}

\item{Theta_h}{B-spline basis for inverse warping functions.}

\item{mean_coefs}{spline coefficient vector for mean curve.}

\item{knots}{knot locations for B-spline basis used to estimate mean and FPC basis function.}

\item{beta.inner}{spline coefficient vector to be estimated for warping function h.}

\item{family}{One of \code{c("gaussian","binomial")}. Defaults to \code{"gaussian"}.}

\item{incompleteness}{Optional specification of incompleteness structure.
One of \code{c("leading","trailing","full")}, specifying that incompleteness
is present only in the initial measurements, only in the trailing measurements, or
in both, respectively. For details see the accompanying vignette.
Defaults to NULL, i.e. no incompleteness structure.
Can only be set when \code{warping = "nonparametric"}.}

\item{lambda_inc}{Penalization parameter to control the amount of
overall dilation of the domain.
The higher this lambda, the more the registered domains are forced to have the
same length as the observed domains.
Only used if \code{incompleteness} is not NULL.}

\item{t_min}{minimum and maximum value to be evaluated on the time domain.}

\item{t_max}{minimum and maximum value to be evaluated on the time domain.}

\item{t_min_curve}{minimum and maximum value of the observed time domain of the
(potentially incomplete) curve.}

\item{t_max_curve}{minimum and maximum value of the observed time domain of the
(potentially incomplete) curve.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions. Default is 8.}

\item{periodic}{If \code{TRUE}, uses periodic b-spline basis functions. Default is \code{FALSE}. 
\code{loss_h_gradient()} is currently only available for \code{periodic = FALSE}.}

\item{warping}{If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.
\code{loss_h_gradient()} is currently only available for \code{warping = "nonparametric"}.}
}
\value{
A numeric vector of spline coefficients for the gradient of the loss function.
}
\description{
Gradient of loss function for registration step
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/registr-utils.R
\name{initial_params}
\alias{initial_params}
\title{Create initial parameters for (inverse) warping functions}
\usage{
initial_params(warping = "nonparametric", K, t_vec)
}
\arguments{
\item{warping}{If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.}

\item{K}{Spline basis matrix defined over the interval \code{c(t_min, t_max)}.}

\item{t_vec}{Vector of the observed and potentially irregular time grid.}
}
\description{
Dependent on the specific type of warping functions, this function creates
a vector of initial parameters. For \code{"nonparametric"} warpings that
are based on a given spline basis matrix, the initial parameters are defined
s.t. the resulting (inverse) warping function equals a diagonal line.
For \code{"piecewise_linear2"} warpings a fixed parameter vector is returned.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_functional_data.R
\name{mean_sim}
\alias{mean_sim}
\title{Simulate mean}
\usage{
mean_sim(grid)
}
\arguments{
\item{grid}{Grid of x values over which to evaluate the function.}
}
\description{
This function generates mean for simulated functional data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bfpca.R
\name{bfpca_argPreparation}
\alias{bfpca_argPreparation}
\title{Internal main preparation function for bfpca}
\usage{
bfpca_argPreparation(
  Y,
  Kt,
  time,
  t_min,
  t_max,
  periodic,
  seed,
  subsample,
  verbose
)
}
\arguments{
\item{Y, time, t_min, t_max}{Internal objects created in \code{bfpca}.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8. If \code{npc_varExplained}
is used, \code{Kt} is set to 20.}

\item{periodic}{If TRUE, uses periodic b-spline basis functions. Default is FALSE.}

\item{seed}{Set seed for reproducibility. Defaults to 1988.}

\item{subsample}{if the number of rows of the data is greater than 
10 million rows, the `id` values are subsampled to get the mean coefficients.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}
}
\value{
List with elements \code{knots}, \code{Theta_phi}, \code{xi},
\code{alpha_coefs}.
}
\description{
Internal main preparation function for bfpca
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_covHall.R
\name{cov_hall}
\alias{cov_hall}
\title{Covariance estimation after Hall et al. (2008)}
\usage{
cov_hall(
  Y,
  index_evalGrid,
  Kt = 25,
  Kc = 10,
  family = "gaussian",
  diag_epsilon = 0.01,
  make_pd = TRUE
)
}
\arguments{
\item{Y}{Dataframe. Should have values id, value, index.}

\item{index_evalGrid}{Grid for the evaluation of the covariance structure.}

\item{Kt}{Number of P-spline basis functions for the estimation of the
marginal mean. Defaults to 25.}

\item{Kc}{Number of marginal P-spline basis functions for smoothing the
covariance surface. Defaults to 10.}

\item{family}{One of \code{c("gaussian","binomial","gamma","poisson")}.
Poisson data are rounded before performing
the GFPCA to ensure integer data, see Details section below.
Defaults to \code{"gaussian"}.}

\item{diag_epsilon}{Small constant to which diagonal elements of the
covariance matrix are set if they are smaller. Defaults to 0.01.}

\item{make_pd}{Indicator if positive (semi-)definiteness of the returned
latent covariance should be ensured via \code{Matrix::near_PD()}. Defaults to
TRUE.}
}
\value{
Covariance matrix with dimension \code{time_evalGrid x time_evalGrid}.
}
\description{
Internal function for the estimation of the covariance matrix of the latent
process using the approach of Hall et al. (2008). Used in the
two-step GFPCA approach implemented in \code{\link{gfpca_twoStep}}. \cr \cr
This function is an adaptation of the implementation of Jan
Gertheiss and Ana-Maria Staicu for Gertheiss et al. (2017), with focus on
higher (RAM) efficiency for large data settings.
}
\details{
The implementation deviates from the algorithm described in Hall (2008) in
one crucial step -- we compute the crossproducts of \emph{centered}
observations and smooth the surface of these crossproducts directly instead
of computing and smoothing the surface of crossproducts of uncentered
observations and subsequently subtracting the (crossproducts of the) mean
function. The former seems to yield smoother eigenfunctions and 
fewer non-positive-definite covariance estimates.

If the data \code{Y} or the crossproduct matrix contain more than
\code{100,000} rows or elements, the estimation of the marginal mean or
the smoothing step of the covariance matrix are performed by
using the discretization-based estimation algorithm in \code{\link[mgcv]{bam}}
rather than the \code{\link[mgcv]{gam}} estimation algorithm.
}
\examples{
data(growth_incomplete)

index_grid = c(1.25, seq(from = 2, to = 18, by = 1))
cov_matrix = registr:::cov_hall(growth_incomplete, index_evalGrid = index_grid)

}
\references{
Hall, P., Müller, H. G., & Yao, F. (2008). Modelling sparse
generalized longitudinal observations with latent Gaussian processes.
\emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
70(4), 703--723.

Gertheiss, J., Goldsmith, J., & Staicu, A. M. (2017). A note on
modeling sparse exponential-family functional response curves.
\emph{Computational statistics & data analysis}, 105, 46--52.
}
\author{
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de} and 
Fabian Scheipl, based on work of Jan Gertheiss and Ana-Maria Staicu
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_functional_data.R
\name{psi2_sim}
\alias{psi2_sim}
\title{Simulate PC2}
\usage{
psi2_sim(grid)
}
\arguments{
\item{grid}{Grid of x values over which to evaluate the function.}
}
\description{
This function generates the second principal component for simulated functional data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_unregistered_curves.R
\name{simulate_unregistered_curves}
\alias{simulate_unregistered_curves}
\title{Simulate unregistered curves}
\usage{
simulate_unregistered_curves(
  I = 50,
  D = 100,
  lambda = 15,
  seed = 1988,
  period = 2 * pi,
  spline_based = FALSE,
  phase_variation = TRUE
)
}
\arguments{
\item{I}{Number of subjects. Defaults is 50.}

\item{D}{Number of grid points per subject. Default is 100.}

\item{lambda}{Standard deviation for subject-specific amplitudes.}

\item{seed}{Seed for reproducibility. Default is 1988.}

\item{period}{Controls the period of the mean curve}

\item{spline_based}{If FALSE curve is constructed using sine and cosine functions,
if TRUE, curve is constructed using B-spline basis.}

\item{phase_variation}{If TRUE, creates phase variation 
(registered curves are observed on uneven grid). If FALSE, no phase variation.}
}
\value{
A simulated dataframe with variables id, value, index, latent_mean, and t. Index is the domain
on which curves are unregistered and t is the domain on which curves are registered.
}
\description{
This function simulates unregistered curves, providing the time values for both 
the unregistered curves (t_star) and the registered curves (t). Curves all have one peak, the location
of which is shifted on the unregistered domain, meant to mimic accelerometer data.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Jeff Goldsmith \email{ajg2202@cumc.columbia.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bfpca.R
\name{bfpca}
\alias{bfpca}
\title{Binary functional principal components analysis}
\usage{
bfpca(
  Y,
  npc = NULL,
  npc_varExplained = NULL,
  Kt = 8,
  maxiter = 50,
  t_min = NULL,
  t_max = NULL,
  print.iter = FALSE,
  row_obj = NULL,
  seed = 1988,
  periodic = FALSE,
  error_thresh = 1e-04,
  verbose = 1,
  subsample = TRUE,
  ...
)
}
\arguments{
\item{Y}{Dataframe. Should have variables id, value, index.}

\item{npc}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_varExplained} has to be set
to a share between 0 and 1.}

\item{npc_varExplained}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_varExplained} has to be set
to a share between 0 and 1.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8. If \code{npc_varExplained}
is used, \code{Kt} is set to 20.}

\item{maxiter}{Maximum number of iterations to perform for EM algorithm. Default is 50.}

\item{t_min}{Minimum value to be evaluated on the time domain.}

\item{t_max}{Maximum value to be evaluated on the time domain.}

\item{print.iter}{Prints current error and iteration}

\item{row_obj}{If NULL, the function cleans the data and calculates row indices. 
Keep this NULL if you are using standalone \code{register} function.}

\item{seed}{Set seed for reproducibility. Defaults to 1988.}

\item{periodic}{If TRUE, uses periodic b-spline basis functions. Default is FALSE.}

\item{error_thresh}{Error threshold to end iterations. Defaults to 0.0001.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{subsample}{if the number of rows of the data is greater than 
10 million rows, the `id` values are subsampled to get the mean coefficients.}

\item{...}{Additional arguments passed to or from other functions}
}
\value{
An object of class \code{fpca} containing:
\item{fpca_type}{Information that FPCA was performed with the 'variationEM' approach,
in contrast to registr::gfpca_twoStep.}
\item{t_vec}{Time vector over which the mean \code{mu} and the functional principal
components \code{efunctions} were evaluated.}
\item{knots}{Cutpoints for B-spline basis used to rebuild \code{alpha}.}
\item{efunctions}{\eqn{D \times npc} matrix of estimated FPC basis functions.}
\item{evalues}{Estimated variance of the FPC scores.}
\item{evalues_sum}{Approximation of the overall variance in \code{Y}, based
on an initial run of the FPCA with \code{npc = 20}. Is \code{NULL} if
\code{npc_varExplained} was not specified.}
\item{npc}{number of FPCs.}
\item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
\item{alpha}{Estimated population-level mean.}
\item{mu}{Estimated population-level mean. Same value as \code{alpha} but included for compatibility
with \code{refund.shiny} package.}
\item{subject_coefs}{B-spline basis coefficients used to construct subject-specific means. 
For use in \code{registr()} function.}
\item{Yhat}{FPC approximation of subject-specific means, before applying the
response function.}
\item{Y}{The observed data.}
\item{family}{\code{binomial}, for compatibility with \code{refund.shiny} package.}
\item{error}{vector containing error for each iteration of the algorithm.}
}
\description{
Function used in the FPCA step for registering binary functional data,
called by \code{\link{register_fpca}} when \code{family = "binomial"}. 
This method uses a variational EM algorithm to estimate scores and principal components for 
binary functional data. \cr \cr
The number of functional principal components (FPCs) can either be specified
directly (argument \code{npc}) or chosen based on the explained share of
variance (\code{npc_varExplained}). In the latter case, the explained share of
variance and accordingly the number of FPCs is estimated before the main
estimation step by once running the FPCA with \code{npc = 20} (and
correspondingly \code{Kt = 20}). Doing so, we approximate the overall
variance in the data \code{Y} with the variance represented by the FPC basis
with 20 FPCs.
}
\examples{
Y = simulate_functional_data()$Y

# estimate 2 FPCs
bfpca_obj = bfpca(Y, npc = 2, print.iter = TRUE, maxiter = 25)


\donttest{
plot(bfpca_obj)

# estimate npc adaptively, to explain 90\% of the overall variation
bfpca_obj2 = bfpca(Y, npc_varExplained = 0.9, print.iter = TRUE, maxiter = 30)
plot(bfpca_obj2)
}
}
\references{
Jaakkola, T. S. and Jordan, M. I. (1997).
A variational approach to Bayesian logistic regression models and their extensions. 
\emph{Proceedings of the Sixth International Workshop on Artificial Intelligence 
and Statistics}.

Tipping, M. E. (1999). Probabilistic Visualisation of High-dimensional binary data.
\emph{Advances in neural information processing systems}, 592--598.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Jeff Goldsmith \email{ajg2202@cumc.columbia.edu},
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{growth_incomplete}
\alias{growth_incomplete}
\title{Berkeley Growth Study data with simulated incompleteness}
\format{
A dataframe made up of \describe{
 \item{id}{A unique subject identifier;}
 \item{index}{Observed age of the child's height;}
 \item{value}{First derivative of the height development in the given age.}
}
}
\usage{
data(growth_incomplete)
}
\description{
This dataset from the Berkeley Growth Study comprises the height
development of 39 boys and 54 girls between ages 1 and 18.
It is based on the dataset \code{fda::growth} and focuses not on the observed
heights, but on the first derivatives of the curves. Before taking the
first derivative, the curves were slightly smoothed. \cr \cr
To showcase the functionality of the \code{registr} package regarding the
analysis of incomplete curves, the growth curves were artificially made
incomplete. For each child, leading incompleteness was simulated by drawing
a random initial age in the first quarter of the domain.
Also, trailing incompleteness was simulated by drawing
a random cut-off age in the second half of the domain.
}
\references{
Ramsay, J. O., and Silverman, B. W. (2006),
\emph{Functional Data Analysis, 2nd ed.}, Springer, New York.

Tuddenham, R. D., and Snyder, M. M. (1954).
Physical growth of California boys and girls from birth to age 18.
\emph{University of California Publications in Child Development}, 1, 183-364.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpca_gauss.R
\name{fpca_gauss_argPreparation}
\alias{fpca_gauss_argPreparation}
\title{Internal main preparation function for fpca_gauss}
\usage{
fpca_gauss_argPreparation(
  Y,
  Kt,
  time,
  t_min,
  t_max,
  periodic,
  seed,
  subsample,
  verbose
)
}
\arguments{
\item{Y, time, t_min, t_max}{Internal objects created in \code{fpca_gauss}.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8. If \code{npc_varExplained}
is used, \code{Kt} is set to 20.}

\item{periodic}{If TRUE, uses periodic b-spline basis functions. Default is FALSE.}

\item{seed}{Set seed for reproducibility. Defaults to 1988.}

\item{subsample}{if the number of rows of the data is greater than 
10 million rows, the `id` values are subsampled to get the mean coefficients.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}
}
\value{
List with elements \code{knots}, \code{Theta_phi}, \code{alpha_coefs}.
}
\description{
Internal main preparation function for fpca_gauss
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_twoStep.R
\name{determine_npc}
\alias{determine_npc}
\title{Determine the number of FPCs based on the share of explained variance}
\usage{
determine_npc(evalues, npc_criterion)
}
\arguments{
\item{evalues}{Vector of estimated variances of the FPC scores.}

\item{npc_criterion}{Either (i) a share between 0 and 1, or (ii) a vector with
two elements for the targeted explained share of variance and a cut-off scree
plot criterion, both between 0 and 1. For the latter, e.g.,
\code{npc_criterion = c(0.9,0.02)} tries to choose a number of FPCs that
explains at least 90\% of variation, but only includes FPCs that explain at
least 2\% of variation (even if this means 90\% explained variation is not reached).}
}
\value{
Integer for the number of fucntional principal components
}
\description{
This internal function is called in \code{gfpca_twoStep}, \code{fpca_gauss}
and \code{bfpca} to determine the number of functional principal components
based on their share of explained variance.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_functional_data.R
\name{simulate_functional_data}
\alias{simulate_functional_data}
\title{Simulate functional data}
\usage{
simulate_functional_data(
  lambda1 = 2,
  lambda2 = 1,
  I = 50,
  D = 100,
  seed = 1988,
  vary_D = FALSE
)
}
\arguments{
\item{lambda1}{Standard deviation for PC1 scores.}

\item{lambda2}{Standard deviation for PC2 scores.}

\item{I}{Number of subjects. Defaults is 50.}

\item{D}{Number of grid points per subject. Default is 100.}

\item{seed}{Seed for reproducibility. Default is 1988.}

\item{vary_D}{Indicates if grid length vary by subject. If FALSE all subjects have grid length D.}
}
\value{
A list containing:
\item{Y}{Simulated dataframe with variables id, value, index, and latent_mean.}
\item{psi1}{True values for first principal component.}
\item{psi2}{True values for second principal component.}
\item{alpha}{True values for population-level mean.}

A list containing:
\item{Y}{A dataframe of simulated data.}
\item{psi1}{The first simulated eigenfunction.}
\item{psi2}{The second simulated eigenfunction.}
\item{alpha}{The population mean.}
}
\description{
This function simulates functional data. The data it outputs is generated from a mean function
and two orthogonal principal component basis functions. The mean and principal components are 
based on sine and cosine functions. Subject-specific scores for each PC are drawn from normal 
distributions with standard deviation lambda1 and lambda2.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help.R, R/registr.R
\docType{package}
\name{registr}
\alias{registr}
\title{Register Exponential Family Functional Data}
\usage{
registr(
  obj = NULL,
  Y = NULL,
  Kt = 8,
  Kh = 4,
  family = "gaussian",
  gradient = TRUE,
  incompleteness = NULL,
  lambda_inc = NULL,
  Y_template = NULL,
  beta = NULL,
  t_min = NULL,
  t_max = NULL,
  row_obj = NULL,
  periodic = FALSE,
  warping = "nonparametric",
  gamma_scales = NULL,
  cores = 1L,
  subsample = TRUE,
  verbose = 1,
  ...
)
}
\arguments{
\item{obj}{Current estimate of FPC object. 
Can be NULL only if Y argument is selected.}

\item{Y}{Dataframe. Should have values id, value, index.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions. Default is 8.}

\item{Kh}{Number of B-spline basis functions used to estimate warping functions \emph{h}. Default is 4.}

\item{family}{One of \code{c("gaussian","binomial","gamma","poisson")}. Defaults to
\code{"gaussian"}.}

\item{gradient}{If \code{TRUE}, uses analytic gradient to calculate derivative. 
If \code{FALSE}, calculates gradient numerically. Not available for families
\code{"gamma","poisson"}.}

\item{incompleteness}{Optional specification of incompleteness structure.
One of \code{c("leading","trailing","full")}, specifying that incompleteness
is present only in the initial measurements, only in the trailing measurements, or
in both, respectively. For details see the accompanying vignette.
Defaults to NULL, i.e. no incompleteness structure.
Can only be set when \code{warping = "nonparametric"}.}

\item{lambda_inc}{Penalization parameter to control the amount of
overall dilation of the domain.
The higher this lambda, the more the registered domains are forced to have the
same length as the observed domains.
Only used if \code{incompleteness} is not NULL.}

\item{Y_template}{Optional dataframe with the same structure as \code{Y}.
Only used if \code{obj} is NULL. If \code{Y_template} is NULL,
curves are registered to the overall mean of all curves in \code{Y} as template function.
If \code{Y_template} is specified, the template function is taken as the mean
of all curves in \code{Y_template}. Default is NULL.}

\item{beta}{Current estimates for beta for each subject. Default is NULL.}

\item{t_min}{Minimum value to be evaluated on the time domain.
if `NULL`, taken to be minimum observed value.}

\item{t_max}{Maximum value to be evaluated on the time domain. 
if `NULL`, taken to be maximum observed value.}

\item{row_obj}{If NULL, the function cleans the data and calculates row indices. 
Keep this NULL if you are using standalone \code{registr} function.}

\item{periodic}{If \code{TRUE}, uses periodic b-spline basis functions. Default is \code{FALSE}.}

\item{warping}{If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.}

\item{gamma_scales}{Only used for \code{family = "gamma"}.
Vector with one entry for each subject, containing the current estimate for the scale parameter of its
gamma distribution. Default is NULL, which sets the starting value for the scale parameter to 1.5.}

\item{cores}{Number of cores to be used. If \code{cores > 1}, the registration
call is parallelized by using \code{parallel::mclapply} (for Unix-based
systems) or \code{parallel::parLapply} (for Windows). Defaults to 1,
no parallelized call.}

\item{subsample}{if the number of rows of the data is greater than 
10 million rows, the `id` values are subsampled to get the mean coefficients.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{...}{additional arguments passed to or from other functions}
}
\value{
An list containing:
\item{Y}{The observed data. The variables \code{index} and \code{index_scaled}
contain the new estimated time domain.}
\item{loss}{Value of the loss function after registraton.}
\item{hinv_innerKnots}{List of inner knots for setting up the spline bases
for the inverse warping functions. Only contains \code{NULL} values for
\code{Kh <= 4}.}
\item{hinv_beta}{Matrix of B-spline basis coefficients used to construct
subject-specific inverse warping functions. See examples on how to
reconstruct a warping function based on \code{hinv_innerKnots} and
\code{hinv_beta}.}
}
\description{
Software for registering functional data from the exponential family of distributions.

Function used in the registration step of an FPCA-based approach for 
registering exponential-family, potentially incomplete functional data,
called by \code{\link{register_fpca}}. 
This method uses constrained optimization to estimate spline 
coefficients for warping functions, where the objective function for optimization comes from 
maximizing the EF likelihood subject to monotonicity constraints on the warping functions. 
You have to either specify \code{obj}, which is a fpca 
object from an earlier step, or \code{Y}, a dataframe in long format with variables 
\code{id}, \code{index}, and \code{value} to indicate subject IDs, times, and observations, 
respectively. \cr \cr
Warping functions by default are forced to start and end on the diagonal to be
domain-preserving. This behavior can be changed by setting
\code{incompleteness} to some other value than NULL and a reasonable \code{lambda_inc} value.
For further details see the accompanying vignette. \cr \cr
By specifying \code{cores > 1} the registration call can be parallelized.
}
\details{
The template function for the registration is defined by argument \code{obj}
or \code{Y_template}, depending on if \code{obj} is NULL or not, respectively.
}
\examples{
### complete binomial curves
Y = simulate_unregistered_curves()
register_step = registr(obj = NULL, Y = Y, Kt = 6, Kh = 4, family = "binomial", 
                        gradient = TRUE)
\donttest{
### incomplete Gaussian curves
data(growth_incomplete)

# Force the warping functions to start and end on the diagonal to preserve the domain
register_step2a = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
                          family = "gaussian", gradient = TRUE,
                          incompleteness = NULL)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  ggplot(register_step2a$Y, aes(x = tstar, y = index, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Estimated warping functions")
  ggplot(register_step2a$Y, aes(x = index, y = value, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Registered curves")
}
  
# Example for how to recreate an estimated inverse warping function given
# the output of registr(). Focus on id "boy01".
id         = "boy01"
index_obsRange_i = range(growth_incomplete$index[growth_incomplete$id == id])
index      = seq(min(index_obsRange_i), max(index_obsRange_i), length.out = 100)
# (note that 'index' must contain both the observed min and max in index_obsRange_i)
Theta_h_i  = splines::bs(index, knots = register_step2a$hinv_innerKnots[[id]], intercept = TRUE)
index_reg  = as.vector(Theta_h_i \%*\% register_step2a$hinv_beta[,id])
warp_dat_i = data.frame(index_observed   = index,
                        index_registered = index_reg)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot(warp_dat_i, aes(x = index_observed, y = index_registered)) + geom_line() +
    ggtitle("Extracted warping function for id 'boy01'")
}

# Allow the warping functions to not start / end on the diagonal.
# The higher lambda_inc, the more the starting points and endpoints are
# forced towards the diagonal.
register_step2b = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
                          family = "gaussian", gradient = TRUE,
                          incompleteness = "full", lambda_inc = 1)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot(register_step2b$Y, aes(x = tstar, y = index, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Estimated warping functions")
  ggplot(register_step2b$Y, aes(x = index, y = value, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Registered curves")
}

# Define the template function only over a subset of the curves
# (even though not very reasonable in this example)
template_ids    = c("girl12","girl13","girl14")
Y_template      = growth_incomplete[growth_incomplete$id \%in\% template_ids,]
register_step2c = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
                          family = "gaussian", gradient = TRUE,
                          Y_template = Y_template,
                          incompleteness = "full", lambda_inc = 1)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot(register_step2c$Y, aes(x = index, y = value, group = id)) +
    geom_line(alpha = 0.2) +
    ggtitle("Registered curves")
}

}

}
\author{
Julia Wrobel

Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Erin McDonnell \email{eim2117@cumc.columbia.edu},
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_helpers.R
\name{deriv.inv.logit}
\alias{deriv.inv.logit}
\title{Estimate the derivative of the logit function}
\usage{
\method{deriv}{inv.logit}(x)
}
\arguments{
\item{x}{Value at which the derivative is computed}
}
\description{
Compute the derivative of the logit function for a given point \code{x}.
}
\author{
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_covHall.R
\name{crossprods_regular}
\alias{crossprods_regular}
\title{Crossproduct computation for mostly regular grids}
\usage{
crossprods_regular(Y)
}
\arguments{
\item{Y}{Dataframe with the centered observations.
Should have values id, centered, index.}
}
\description{
Compute the crossproduct in a fast way for mostly regular grids
(index values are mostly *not* unique).
Only used internally in \code{cov_hall()}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bfpca.R
\name{bfpca_optimization}
\alias{bfpca_optimization}
\title{Internal main optimization for bfpca}
\usage{
bfpca_optimization(
  npc,
  npc_varExplained = NULL,
  Kt,
  maxiter,
  print.iter,
  seed,
  periodic,
  error_thresh,
  verbose,
  Y,
  rows,
  I,
  knots,
  Theta_phi,
  xi,
  alpha_coefs
)
}
\arguments{
\item{npc}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_varExplained} has to be set
to a share between 0 and 1.}

\item{npc_varExplained}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_varExplained} has to be set
to a share between 0 and 1.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8. If \code{npc_varExplained}
is used, \code{Kt} is set to 20.}

\item{maxiter}{Maximum number of iterations to perform for EM algorithm. Default is 50.}

\item{print.iter}{Prints current error and iteration}

\item{seed}{Set seed for reproducibility. Defaults to 1988.}

\item{periodic}{If TRUE, uses periodic b-spline basis functions. Default is FALSE.}

\item{error_thresh}{Error threshold to end iterations. Defaults to 0.0001.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{Y, rows, I, knots, Theta_phi, xi, alpha_coefs}{Internal objects created in
\code{bfpca}.}
}
\value{
list with elements \code{t_vec}, \code{Theta_phi_mean}, \code{alpha_coefs},
\code{efunctions}, \code{evalues}, \code{evalues_sum}, \code{scores},
\code{subject_coef}, \code{fittedVals}, \code{error}. See documentation of
\code{\link{fpca_gauss}} for details.
}
\description{
Main optimization function for \code{bfpca}. If \code{npc_varExplained}
is specified, the function simply returns a list with elements \code{npc}
(chosen number of FPCs), \code{evalues} (estimated variances of the first 'npc'
FPCs) and \code{evalues_sum} (sum of the estimated variances of the first 20
FPCs, as approximation of the overall variance).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{expectedScores}
\alias{expectedScores}
\title{Calculate expected score and score variance for the current subject.}
\usage{
expectedScores(Y, mu, psi, theta, theta_quad)
}
\arguments{
\item{Y}{vector of observations for the current subject.}

\item{mu}{vector of spline coefficients for the population mean.}

\item{psi}{matrix of spline coefficients for the principal component basis functions.}

\item{theta}{spline basis functions for the current subject.}

\item{theta_quad}{quadratic form of theta for the current subject.}
}
\value{
A list with expected score mean and variance for the current subject.
}
\description{
Calculations derived using maximum likelihood estimation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpca_gauss.R
\name{fpca_gauss_optimization}
\alias{fpca_gauss_optimization}
\title{Internal main optimization for fpca_gauss}
\usage{
fpca_gauss_optimization(
  npc,
  npc_varExplained = NULL,
  Kt,
  maxiter,
  print.iter,
  seed,
  periodic,
  error_thresh,
  verbose,
  Y,
  rows,
  I,
  knots,
  Theta_phi,
  alpha_coefs
)
}
\arguments{
\item{npc}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_varExplained} has to be set
to a share between 0 and 1.}

\item{npc_varExplained}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_varExplained} has to be set
to a share between 0 and 1.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8. If \code{npc_varExplained}
is used, \code{Kt} is set to 20.}

\item{maxiter}{Maximum number of iterations to perform for EM algorithm. Default is 50.}

\item{print.iter}{Prints current error and iteration}

\item{seed}{Set seed for reproducibility. Defaults to 1988.}

\item{periodic}{If TRUE, uses periodic b-spline basis functions. Default is FALSE.}

\item{error_thresh}{Error threshold to end iterations. Defaults to 0.0001.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{Y, rows, I, knots, Theta_phi, alpha_coefs}{Internal objects created in
\code{fpca_gauss}.}
}
\value{
list with elements \code{t_vec}, \code{Theta_phi_mean}, \code{alpha_coefs},
\code{efunctions}, \code{evalues}, \code{evalues_sum}, \code{scores},
\code{subject_coef}, \code{fittedVals}, \code{sigma2}. See documentation of
\code{\link{fpca_gauss}} for details.
}
\description{
Main optimization function for \code{fpca_gauss}. If \code{npc_varExplained}
is specified, the function simply returns a list with elements \code{npc}
(chosen number of FPCs), \code{evalues} (estimated variances of the first 'npc'
FPCs) and \code{evalues_sum} (sum of the estimated variances of the first 20
FPCs, as approximation of the overall variance).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constrOptim_helpers.R
\name{ensure_proper_beta}
\alias{ensure_proper_beta}
\title{Correct slightly improper parameter vectors}
\usage{
ensure_proper_beta(beta, t_min, t_max)
}
\arguments{
\item{beta}{Parameter vector.}

\item{t_min, t_max}{Minimum and maximum of the underlying time domain in the
registration step.}
}
\value{
A slightly changed parameter vector that ensures a proper solution
in the optimization of the registration step.
}
\description{
Internal function. In the joint iterations between registration and GFPCA,
the optimization with \code{constrOptim()} in the registration step sometimes
leads to slightly improper solutions, which cause the optimization to
throw an error in the following optimization step. This function corrects
the parameter vector if one of the following slight inconsistencies occurs
that can mess with the optimization of \code{constrOptim()}: \cr
- two neighboring values of the parameter vector are too similar \cr
- the initial values of the parameter vector are smaller than \code{t_min},
the minimum of the underlying time domain \cr
- the last values of the parameter vector are greater than \code{t_max},
the maximum of the underlying time domain \cr
- one parameter value is slightly greater than its following value, i.e.
the parameter vector is not monotone.
}
\examples{
beta_improper = c(0.24, 1.000047, 1.000002)
registr:::ensure_proper_beta(beta_improper, t_min = 0, t_max = 1)
}
\author{
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss_h.R
\name{loss_h}
\alias{loss_h}
\title{Loss function for registration step optimization}
\usage{
loss_h(
  Y,
  Theta_h,
  mean_coefs,
  knots,
  beta.inner,
  family,
  t_min,
  t_max,
  t_min_curve,
  t_max_curve,
  incompleteness = NULL,
  lambda_inc = NULL,
  periodic = FALSE,
  Kt = 8,
  warping = "nonparametric",
  priors = FALSE,
  prior_sd = NULL
)
}
\arguments{
\item{Y}{vector of observed points.}

\item{Theta_h}{B-spline basis for inverse warping functions.}

\item{mean_coefs}{spline coefficient vector for mean curve.}

\item{knots}{knot locations for B-spline basis used to estimate mean and FPC basis function.}

\item{beta.inner}{spline coefficient vector to be estimated for warping function h.}

\item{family}{One of \code{c("gaussian","binomial","gamma","poisson")}.
For internal purposes, can also be set to \code{"gamma-varEM"} and
\code{"poisson-varEM"} if the preceding FPCA step in \code{register_fpca} was
performed with \code{fpca_type = "variationalEM"} which uses Gaussian family.}

\item{t_min, t_max}{minimum and maximum value to be evaluated on the time domain.}

\item{t_min_curve, t_max_curve}{minimum and maximum value of the observed time domain of the
(potentially incomplete) curve.}

\item{incompleteness}{Optional specification of incompleteness structure.
One of \code{c("leading","trailing","full")}, specifying that incompleteness
is present only in the initial measurements, only in the trailing measurements, or
in both, respectively. For details see the accompanying vignette.
Defaults to NULL, i.e. no incompleteness structure.
Can only be set when \code{warping = "nonparametric"}.}

\item{lambda_inc}{Penalization parameter to control the amount of
overall dilation of the domain.
The higher this lambda, the more the registered domains are forced to have the
same length as the observed domains.
Only used if \code{incompleteness} is not NULL.}

\item{periodic}{If \code{TRUE} uses periodic b-spline basis functions. Default is \code{FALSE}.}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions. Default is 8.}

\item{warping}{If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.}

\item{priors}{For \code{warping = "piecewise_linear2"} only. Logical indicator of whether to add Normal priors to pull the knots toward the identity line.}

\item{prior_sd}{For \code{warping = "piecewise_linear2"} with \code{priors = TRUE} only. User-specified standard deviation for the Normal priors 
(single value applied to all 4 knot priors).}
}
\value{
The scalar value taken by the loss function.
}
\description{
Loss function for registration step optimization
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Erin McDonnell \email{eim2117@cumc.columbia.edu},
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_twoStep.R
\name{gfpca_twoStep}
\alias{gfpca_twoStep}
\title{Generalized functional principal component analysis}
\usage{
gfpca_twoStep(
  Y,
  family = "gaussian",
  npc = NULL,
  npc_criterion = NULL,
  Kt = 8,
  t_min = NULL,
  t_max = NULL,
  row_obj = NULL,
  index_significantDigits = 4L,
  estimation_accuracy = "high",
  start_params = NULL,
  periodic = FALSE,
  verbose = 1,
  ...
)
}
\arguments{
\item{Y}{Dataframe. Should have values id, value, index.}

\item{family}{One of \code{c("gaussian","binomial","gamma","poisson")}.
Poisson data are rounded before performing
the GFPCA to ensure integer data, see Details section below.
Defaults to \code{"gaussian"}.}

\item{npc, npc_criterion}{The number of functional principal components (FPCs)
has to be specified either directly as \code{npc} or based on their explained
share of variance. In the latter case, \code{npc_criterion} can either be set
to (i) a share between 0 and 1, or (ii) a vector with two elements comprising
the targeted explained share of variance and a cut-off scree plot criterion,
both between 0 and 1. As an example for the latter,
\code{npc_criterion = c(0.9,0.02)} tries to choose a number of FPCs that
explains at least 90\% of variation, but only includes FPCs that explain at
least 2\% of variation (even if this means 90\% explained variation is not reached).}

\item{Kt}{Number of B-spline basis functions used to estimate mean functions
and functional principal components. Default is 8.}

\item{t_min}{Minimum value to be evaluated on the time domain.}

\item{t_max}{Maximum value to be evaluated on the time domain.}

\item{row_obj}{If NULL, the function cleans the data and calculates row indices. 
Keep this NULL if you are using standalone \code{register} function.}

\item{index_significantDigits}{Positive integer \code{>= 2}, stating the number
of significant digits to which the index grid should be rounded. Coarsening the
index grid is necessary since otherwise the covariance surface matrix
explodes in size in the presence of too many unique index values (which is
always the case after some registration step). Defaults to 4. Set to
\code{NULL} to prevent rounding.}

\item{estimation_accuracy}{One of \code{c("high","low")}. When set to \code{"low"},
the mixed model estimation step in \code{lme4} is performed with lower
accuracy, reducing computation time. Defaults to \code{"high"}.}

\item{start_params}{Optional start values for gamm4. Not used if
\code{npc_criterion} is specified.}

\item{periodic}{Only contained for full consistency with \code{fpca_gauss}
and \code{bfpca}. If TRUE, returns the knots vector for periodic b-spline
basis functions. Defaults to FALSE. This parameter does not change the
results of the two-step GFPCA.}

\item{verbose}{Can be set to integers between 0 and 4 to control the level of
detail of the printed diagnostic messages. Higher numbers lead to more detailed
messages. Defaults to 1.}

\item{...}{Additional arguments passed to \code{\link{cov_hall}}.}
}
\value{
An object of class \code{fpca} containing:
\item{fpca_type}{Information that FPCA was performed with the 'two-step' approach,
in contrast to registr::fpca_gauss or registr::bfpca.}
\item{t_vec}{Time vector over which the mean \code{mu} was evaluated.
The resolution is can be specified by setting \code{index_significantDigits}.}
\item{knots}{Cutpoints for B-spline basis used to rebuild \code{alpha}.}
\item{efunctions}{\eqn{D \times npc} matrix of estimated FPC basis functions.}
\item{evalues}{Estimated variance of the FPC scores.}
\item{evalues_sum}{Sum of all (nonnegative) eigenvalues of the smoothed
covariance surface estimated with \code{\link{cov_hall}}. Can be used as an
approximation for the total variance present in \code{Y} to compute the
shares of explained variance of the FPC scores.}
\item{npc}{number of FPCs.}
\item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
\item{alpha}{Estimated population-level mean.}
\item{mu}{Estimated population-level mean. Same value as \code{alpha} but included for compatibility
with \code{refund.shiny} package.}
\item{subject_coefs}{Always \code{NA} but included for full consistency
with \code{fpca_gauss} and \code{bfpca}.} 
\item{Yhat}{FPC approximation of subject-specific means, before applying the
response function.}
\item{Y}{The observed data.}
\item{family}{\code{binomial}, for compatibility with \code{refund.shiny} package.}
\item{gamm4_theta}{Estimated parameters of the mixed model.}
}
\description{
Function for applying FPCA to different exponential family distributions.
Used in the FPCA step for registering functional data,
called by \code{\link{register_fpca}} when \code{fpca_type = "two-step"}. \cr \cr
The method implements the `two-step approach` of Gertheiss et al. (2017)
and is based on the approach of Hall et al. (2008) to estimate functional
principal components. \cr \cr
The number of functional principal components (FPCs) can either be specified
directly (argument \code{npc}) or chosen based on the explained share of
variance (\code{npc_criterion}). Using the latter, we approximate the overall
variance in the data \code{Y} with the variance represented by the smoothed
covariance surface estimated with \code{\link{cov_hall}}.
Note that the Eigenvalue decomposition of this covariance surface
sometimes leads to a long tail of subordinate FPCs with small eigenvalues.
Such subordinate dimensions seem to often represent phase rather than
amplitude variation, and can be cut off by specifying the second element of
argument \code{npc_criterion}. \cr \cr
This function is an adaptation of the implementation of Jan
Gertheiss for Gertheiss et al. (2017), with focus on higher (RAM) efficiency
for large data settings.
}
\details{
For \code{family = "poisson"} the values in \code{Y} are rounded before
performing the GFPCA to ensure integer data. This is done to ensure reasonable
computation times. Computation times tend to explode when estimating the
underlying high-dimensional mixed model with continuous Poisson data based
on the \code{\link{gamm4}} package.

If negative eigenvalues are present, the respective eigenfunctions are dropped
and not considered further.
}
\examples{
data(growth_incomplete)

# estimate 2 FPCs
fpca_obj = gfpca_twoStep(Y = growth_incomplete, npc = 2, family = "gaussian")
plot(fpca_obj)

# estimate npc adaptively, to explain 90\% of the overall variation
fpca_obj2 = gfpca_twoStep(Y = growth_incomplete, npc_criterion = 0.9, family = "gaussian")
plot(fpca_obj2, plot_FPCs = 1:2)

}
\references{
Gertheiss, J., Goldsmith, J., & Staicu, A. M. (2017). A note on
modeling sparse exponential-family functional response curves.
\emph{Computational statistics & data analysis}, 105, 46--52.

Hall, P., Müller, H. G., & Yao, F. (2008). Modelling sparse
generalized longitudinal observations with latent Gaussian processes.
\emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
70(4), 703--723.
}
\author{
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de},
based on work of Jan Gertheiss
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{lambdaF}
\alias{lambdaF}
\title{Apply lambda transformation of variational parameter.}
\usage{
lambdaF(x)
}
\arguments{
\item{x}{The value to which you apply the function}
}
\value{
A numeric value that has been transformed.
}
\description{
Simple function for use within other C++ functions.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{nhanes}
\alias{nhanes}
\title{NHANES activity data}
\format{
A dataframe made up of \describe{
 \item{id}{A unique subject identifier;}
 \item{age}{Age of survey participant;}
 \item{gender}{Gender of survey participant;}
 \item{index}{Observed time of activity measurement. Integers from 1 to 1440, indicating minutes
 from midnight to midnight;}
 \item{value}{Binary value of zero or one indicating inactivity or activity;}
 \item{raw_activity}{Raw activity count.}
}
}
\usage{
data(nhanes)
}
\description{
Subset of 24 hours of activity data for 50 subjects from 2003-2004
National Health and Nutrition Examination Survey (NHANES). 
Each subject is observed over 24 hours on a Sunday and wore the 
activity collection device for a minimum of 10 hours. Activity is measured each minute over 24 hours.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_covHall.R
\name{crossprods_irregular}
\alias{crossprods_irregular}
\title{Crossproduct computation for highly irregular grids}
\usage{
crossprods_irregular(Y)
}
\arguments{
\item{Y}{Dataframe with the centered observations.
Should have values id, centered, index.}
}
\description{
Compute the crossproduct in a fast way for highly irregular grids
(index values are mostly unique).
Only used internally in \code{cov_hall()}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\name{constraints}
\alias{constraints}
\title{Define constraints for optimization of warping functions}
\usage{
constraints(Kh, t_min = 0, t_max = 1, warping = "nonparametric")
}
\arguments{
\item{Kh}{Number of B-spline basis functions used to estimate warping functions \emph{h}.}

\item{t_min}{Minimum value to be evaluated on the time domain.}

\item{t_max}{Maximum value to be evaluated on the time domain.}

\item{warping}{If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.}
}
\value{
An list containing:
\item{ui}{A constraint matrix.}
\item{ci}{A constraint vector.}
}
\description{
Constraints ensure monotonicity of spline coefficients for warping functions 
for use with \code{constrOptim()} function.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu},
Erin McDonnell \email{eim2117@cumc.columbia.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_deriv.R
\name{bs_deriv}
\alias{bs_deriv}
\title{Nth derivative of spline basis}
\usage{
bs_deriv(
  x,
  knots,
  degree = 3L,
  Boundary.knots = range(x),
  derivative = 1,
  intercept = TRUE
)
}
\arguments{
\item{x}{a numeric vector of values at which to evaluate the B-spline functions or derivatives.}

\item{knots}{the internal breakpoints that define the spline.}

\item{degree}{degree of the piecewise polynomial—default is 3 for cubic splines.}

\item{Boundary.knots}{boundary points at which to anchor the B-spline basis. 
Set to [0,1] if you want this to be your domain.}

\item{derivative}{a positive integer value that specifies which derivative to take. Defaults to 1 for 1st derivative.
Value of 0 returns the original set of b-spline basis functions.}

\item{intercept}{if TRUE, an intercept is included in the basis; default is TRUE}
}
\value{
A matrix containing:
\item{basis}{A B-spline basis that can be used to approximate the derivative of a function.}
}
\description{
This function gets derivative of a spline basis. Adapted from \code{bs()} function in \code{splines} package.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_clean.R
\name{data_clean}
\alias{data_clean}
\title{Convert data to a \code{refund} object}
\usage{
data_clean(data)
}
\arguments{
\item{data}{Dataframe. Should have values id, value, index.}
}
\value{
An list containing:
\item{Y}{The original data sorted by id and index.}
\item{Y_rows}{A dataframe containing the first and last row for each subject.}
}
\description{
Function used for data cleaning.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_unregistered_curves.R
\name{amp_curve}
\alias{amp_curve}
\title{Simulate amplitude variance}
\usage{
amp_curve(grid, period = 2 * pi, spline_based = FALSE)
}
\arguments{
\item{grid}{Grid of x values over which to evaluate the function.}

\item{period}{Controls the period of the mean curve}

\item{spline_based}{If FALSE curve is constructed using sine and cosine functions,
if TRUE, curve is constructed using B-spline basis.}
}
\value{
A numeric vector.
}
\description{
This function generates amplitudes for simulated accelerometer data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_functional_data.R
\name{psi1_sim}
\alias{psi1_sim}
\title{Simulate PC1}
\usage{
psi1_sim(grid)
}
\arguments{
\item{grid}{Grid of x values over which to evaluate the function.}
}
\description{
This function generates the first principal component for simulated functional data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{expectedXi}
\alias{expectedXi}
\title{Estimate variational parameter for the current subject.}
\usage{
expectedXi(theta, mu, mi, psi, Ci)
}
\arguments{
\item{theta}{spline basis functions for the current subject.}

\item{mu}{vector of spline coefficients for the population mean.}

\item{mi}{vector of expected mean scores for the current subject.}

\item{psi}{matrix of spline coefficients for the principal component basis functions.}

\item{Ci}{expected covariance matrix of scores for the current subject.}
}
\value{
A vector of variational parameters for the current subject.
}
\description{
Function calculates value of variational parameter using maximum likelihood.
}
