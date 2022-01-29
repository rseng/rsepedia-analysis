
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
