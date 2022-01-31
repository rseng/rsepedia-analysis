
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EFAtools

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/EFAtools)](https://CRAN.R-project.org/package=EFAtools)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02521/status.svg)](https://doi.org/10.21105/joss.02521)
[![R-CMD-check](https://github.com/mdsteiner/EFAdiff/workflows/R-CMD-check/badge.svg)](https://github.com/mdsteiner/EFAdiff/actions)
<!-- badges: end -->

The EFAtools package provides functions to perform exploratory factor
analysis (EFA) procedures and compare their solutions. The goal is to
provide state-of-the-art factor retention methods and a high degree of
flexibility in the EFA procedures. This way, implementations from R
psych and SPSS can be compared. Moreover, functions for Schmid-Leiman
transformation, and computation of omegas are provided. To speed up the
analyses, some of the iterative procedures like principal axis factoring
(PAF) are implemented in C++.

## Installation

You can install the release version from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("EFAtools")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("mdsteiner/EFAtools")
```

To also build the vignette when installing the development version, use:

``` r
install.packages("devtools")
devtools::install_github("mdsteiner/EFAtools", build_vignettes = TRUE)
```

## Example

Here are a few examples on how to perform the analyses with the
different types and how to compare the results using the `COMPARE`
function. For more details, see the vignette by running
`vignette("EFAtools", package = "EFAtools")`. The vignette provides a
high-level introduction into the functionalities of the package.

``` r
# load the package
library(EFAtools)

# Run all possible factor retention methods
N_FACTORS(test_models$baseline$cormat, N = 500, method = "ML")
#> Warning in N_FACTORS(test_models$baseline$cormat, N = 500, method = "ML"): ! 'x' was a correlation matrix but CD needs raw data. Skipping CD.
#>                                                                                                                                                                  ◉ 🏃 ◯ ◯ ◯ ◯ ◯ Running EKC                                                                                                                                                                 ◉ ◉ 🏃 ◯ ◯ ◯ ◯ Running HULL                                                                                                                                                                 ◉ ◉ ◉ 🏃 ◯ ◯ ◯ Running KGC                                                                                                                                                                 ◉ ◉ ◉ ◉ 🏃 ◯ ◯ Running PARALLEL                                                                                                                                                                 ◉ ◉ ◉ ◉ ◉ 🏃 ◯ Running SCREE                                                                                                                                                                 ◉ ◉ ◉ ◉ ◉ ◉ 🏃  Running SMT                                                                                                                                                                 ◉ ◉ ◉ ◉ ◉ ◉ ◉ Done!
#> 
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> Bartlett's test of sphericity
#> 
#> ✓ The Bartlett's test of sphericity was significant at an alpha level of .05.
#>   These data are probably suitable for factor analysis.
#> 
#>   𝜒²(153) = 2173.28, p < .001
#> 
#> Kaiser-Meyer-Olkin criterion (KMO)
#> 
#> ✓ The overall KMO value for your data is marvellous with 0.916.
#>   These data are probably suitable for factor analysis.
#> 
#> ── Number of factors suggested by the different factor retention criteria ──────
#> 
#> ◌ Comparison data: NA
#> ◌ Empirical Kaiser criterion: 2
#> ◌ Hull method with CAF: 3
#> ◌ Hull method with CFI: 1
#> ◌ Hull method with RMSEA: 1
#> ◌ Kaiser-Guttman criterion with PCA: 3
#> ◌ Kaiser-Guttman criterion with SMC: 1
#> ◌ Kaiser-Guttman criterion with EFA: 1
#> ◌ Parallel analysis with PCA: 3
#> ◌ Parallel analysis with SMC: 3
#> ◌ Parallel analysis with EFA: 3
#> ◌ Sequential 𝜒² model tests: 3
#> ◌ Lower bound of RMSEA 90% confidence interval: 2
#> ◌ Akaike Information Criterion: 3
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" />

``` r
# A type SPSS EFA to mimick the SPSS implementation with
# promax rotation
EFA_SPSS <- EFA(test_models$baseline$cormat, n_factors = 3, type = "SPSS",
                  rotation = "promax")

# look at solution
EFA_SPSS
#> 
#> EFA performed with type = 'SPSS', method = 'PAF', and rotation = 'promax'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>       F1      F2      F3  
#> V1   -.048    .035    .613
#> V2   -.001    .067    .482
#> V3    .060    .056    .453
#> V4    .101   -.009    .551
#> V5    .157   -.018    .438
#> V6   -.072   -.049    .704
#> V7    .001    .533    .093
#> V8   -.016    .581    .030
#> V9    .038    .550   -.001
#> V10  -.022    .674   -.071
#> V11   .015    .356    .232
#> V12   .020    .651   -.010
#> V13   .614    .086   -.067
#> V14   .548   -.068    .088
#> V15   .561    .128   -.070
#> V16   .555   -.050    .091
#> V17   .664   -.037   -.027
#> V18   .555    .004    .050
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1      F2      F3  
#> F1    1.000   0.617   0.648
#> F2    0.617   1.000   0.632
#> F3    0.648   0.632   1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                       F1      F2      F3  
#> SS loadings           4.907   0.757   0.643
#> Prop Tot Var          0.273   0.042   0.036
#> Cum Prop Tot Var      0.273   0.315   0.350
#> Prop Comm Var         0.778   0.120   0.102
#> Cum Prop Comm Var     0.778   0.898   1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .50
#> df: 102

# A type psych EFA to mimick the psych::fa() implementation with
# promax rotation
EFA_psych <- EFA(test_models$baseline$cormat, n_factors = 3, type = "psych",
                  rotation = "promax")
```

<img src="man/figures/README-unnamed-chunk-5-3.png" width="100%" />

``` r
# compare the type psych and type SPSS implementations
COMPARE(EFA_SPSS$rot_loadings, EFA_psych$rot_loadings,
        x_labels = c("SPSS", "psych"))
#> Mean [min, max] absolute difference:  0.0090 [ 0.0001,  0.0245]
#> Median absolute difference:  0.0095
#> Max decimals where all numbers are equal: 0
#> Minimum number of decimals provided: 17
#> 
#>        F1      F2      F3  
#> V1    0.0150  0.0142 -0.0195
#> V2    0.0109  0.0109 -0.0138
#> V3    0.0095  0.0103 -0.0119
#> V4    0.0118  0.0131 -0.0154
#> V5    0.0084  0.0105 -0.0109
#> V6    0.0183  0.0169 -0.0245
#> V7   -0.0026 -0.0017  0.0076
#> V8   -0.0043 -0.0035  0.0102
#> V9   -0.0055 -0.0040  0.0117
#> V10  -0.0075 -0.0066  0.0151
#> V11   0.0021  0.0029  0.0001
#> V12  -0.0064 -0.0050  0.0136
#> V13  -0.0109 -0.0019  0.0163
#> V14  -0.0049  0.0028  0.0070
#> V15  -0.0107 -0.0023  0.0161
#> V16  -0.0051  0.0028  0.0074
#> V17  -0.0096 -0.0001  0.0136
#> V18  -0.0066  0.0014  0.0098
```

<img src="man/figures/README-unnamed-chunk-5-4.png" width="100%" />

``` r
# Average solution across many different EFAs with oblique rotations
EFA_AV <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = c("PAF", "ML", "ULS"), rotation = "oblique",
                      show_progress = FALSE)

# look at solution
EFA_AV
#> 
#> Averaging performed with averaging method mean (trim = 0) across 162 EFAs, varying the following settings: method, init_comm, criterion_type, start_method, rotation, k_promax, P_type, and varimax_type.
#> 
#> The error rate is at 0%. Of the solutions that did not result in an error, 100% converged, 0% contained Heywood cases, and 100% were admissible.
#> 
#> 
#> ══ Indicator-to-Factor Correspondences ═════════════════════════════════════════
#> 
#> For each cell, the proportion of solutions including the respective indicator-to-factor correspondence. A salience threshold of 0.3 was used to determine indicator-to-factor correspondences.
#> 
#>       F1      F2      F3 
#> V1    .11     .00    1.00
#> V2    .11     .00    1.00
#> V3    .11     .00     .94
#> V4    .11     .00    1.00
#> V5    .11     .00     .94
#> V6    .11     .00    1.00
#> V7    .11     .94     .00
#> V8    .11    1.00     .00
#> V9    .11     .94     .00
#> V10   .11    1.00     .00
#> V11   .11     .89     .00
#> V12   .11    1.00     .00
#> V13  1.00     .00     .00
#> V14  1.00     .00     .00
#> V15  1.00     .00     .00
#> V16  1.00     .00     .00
#> V17  1.00     .00     .00
#> V18  1.00     .00     .00
#> 
#> 
#> ══ Loadings ════════════════════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>       F1      F2      F3  
#> V1    .025    .048    .576
#> V2    .060    .077    .451
#> V3    .115    .066    .425
#> V4    .157    .007    .518
#> V5    .198   -.002    .412
#> V6    .002   -.028    .658
#> V7    .074    .497    .102
#> V8    .056    .538    .046
#> V9    .100    .510    .018
#> V10   .048    .625   -.046
#> V11   .082    .336    .228
#> V12   .094    .606    .007
#> V13   .597    .083   -.047
#> V14   .531   -.056    .093
#> V15   .548    .122   -.049
#> V16   .540   -.041    .097
#> V17   .633   -.033   -.009
#> V18   .542    .009    .060
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>       F1      F2      F3  
#> V1    0.513   0.086   0.239
#> V2    0.431   0.093   0.186
#> V3    0.394   0.108   0.179
#> V4    0.415   0.110   0.214
#> V5    0.315   0.122   0.177
#> V6    0.514   0.104   0.267
#> V7    0.527   0.255   0.089
#> V8    0.520   0.275   0.078
#> V9    0.470   0.276   0.080
#> V10   0.533   0.313   0.097
#> V11   0.482   0.176   0.102
#> V12   0.548   0.324   0.103
#> V13   0.081   0.289   0.114
#> V14   0.063   0.220   0.117
#> V15   0.091   0.280   0.107
#> V16   0.072   0.230   0.122
#> V17   0.108   0.270   0.124
#> V18   0.081   0.246   0.118
#> 
#> 
#> ══ Factor Intercorrelations from Oblique Solutions ═════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>       F1      F2      F3  
#> F1    1.000   0.431   0.518
#> F2    0.431   1.000   0.454
#> F3    0.518   0.454   1.000
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>       F1      F2      F3  
#> F1    0.000   1.276   0.679
#> F2    1.276   0.000   1.316
#> F3    0.679   1.316   0.000
#> 
#> 
#> ══ Variances Accounted for ═════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>                   F1      F2      F3  
#> SS loadings       2.443   1.929   1.904
#> Prop Tot Var      0.136   0.107   0.106
#> Prop Comm Var     0.389   0.307   0.303
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>                   F1      F2      F3  
#> SS loadings       2.831   1.356   1.291
#> Prop Tot Var      0.157   0.075   0.072
#> Prop Comm Var     0.419   0.215   0.215
#> 
#> 
#> ══ Model Fit ═══════════════════════════════════════════════════════════════════
#> 
#>        M (SD) [Min; Max]
#> 𝜒²: 101.73 (34.62) [53.23; 125.98]
#> df: 102
#> p: .369 (.450) [.054; 1.000]
#> CFI: 1.00 (.00) [1.00; 1.00]
#> RMSEA: .01 (.01) [.00; .02]
#> AIC: -102.27 (34.62) [-150.77; -78.02]
#> BIC: -532.16 (34.62) [-580.66; -507.91]
#> CAF: .50 (.00) [.50; .50]
```

<img src="man/figures/README-unnamed-chunk-5-5.png" width="100%" />

``` r
# Perform a Schmid-Leiman transformation
SL <- SL(EFA_psych)

# Based on a specific salience threshold for the loadings (here: .20):
factor_corres <- SL$sl[, c("F1", "F2", "F3")] >= .2

# Compute omegas from the Schmid-Leiman solution
OMEGA(SL, factor_corres = factor_corres)
#> Omega total, omega hierarchical, and omega subscale for the general factor (top row) and the group factors:
#> 
#>      tot  hier   sub
#> g  0.883 0.750 0.122
#> F1 0.769 0.498 0.272
#> F2 0.764 0.494 0.270
#> F3 0.745 0.543 0.202
```

## Citation

If you use this package in your research, please acknowledge it by
citing:

Steiner, M.D., & Grieder, S.G. (2020). EFAtools: An R package with fast
and flexible implementations of exploratory factor analysis tools.
*Journal of Open Source Software*, *5*(53), 2521.
<https://doi.org/10.21105/joss.02521>

## Contribute or Report Bugs

If you want to contribute or report bugs, please open an issue on GitHub
or email us at <markus.d.steiner@gmail.com> or
<silvia.grieder@gmail.com>.
# EFAtools 0.3.1.9000

## Changes to Functions

* `OMEGA()`: 
    * Added calculation of additional indices of interpretive relevance (H index, explained common variance [ECV], and percent of uncontaminated correlations [PUC]). This is optional and can be avoided by setting `add_ind = FALSE`.
    
    
# EFAtools 0.3.1

## General
* When testing for whether a matrix is singular and thus smoothing should be done, test against .Machine$double.eps^.6 instead of 0, as suggested by Florian Scharf. 

## Changes to Functions

* `EFA()`: 
    * Added warnings if `type = "SPSS"` was used with `method = "ML"` or `method = "ULS"`, or with a rotation other than `none`, `varimax` or `promax`.
    * Avoided smoothing of non-positive definite correlation matrices if `type = "SPSS"` is used.
    * Use Moore-Penrose Pseudo Inverse in computation of SMCs if `type = "psych"` is used, by calling `psych::smc()`.
    * Use `varimax_type = "kaiser"` if `type = "EFAtools"` is used with `varimax` or `promax`.

## Bug Fixes
* `EFA_AVERAGE()`:
    * Added `future.seed = TRUE` to call to `future.apply::future_lapply()` to prevent warnings.
    * Fixed test for Heywood cases from testing whether a communality or loading is greater than .998, to only test whether communalities exceed 1 + .Machine$double.eps
* `print.EFA()`: Fixed test for Heywood cases from testing whether a communality or loading is greater than .998, to only test whether communalities exceed 1 + .Machine$double.eps
* `OMEGA()`: Small bugfix when `lavaan` second-order model is given as input


# EFAtools 0.3.0

## General
* Added examples for `EFA_AVERAGE()` to readme and the EFAtools vignette
* Updated examples in readme and vignettes according to the updated `OMEGA` function

## New Functions

* Added function `EFA_AVERAGE()` and respective print and plot methods, to allow running many EFAs across different implementations to obtain an average solution and test the stability of the results.

## Changes to Functions

* `EFA()`: Defaults that were previously set to `NULL` are now mostly set to `NA`. This was necessary for `EFA_AVERAGE()` to work correctly.
* `PARALLEL()`: Rewrote the generation of random data based eigenvalues to be more stable when SMCs are used.
* `OMEGA()`: Changed expected input for argument `factor_corres` from vector to matrix. Can now be a logical matrix or a numeric matrix with 0's and 1's of the same dimensions as the matrix of group factor loadings. This is more flexible and allows for cross-loadings.

# EFAtools 0.2.0

## General

* Created new vignette *Replicate_SPSS_psych* to show replication of original `psych` and `SPSS` EFA solutions with `EFAtools`.

## New Functions

* Added function `FACTOR_SCORES()` to calculate factor scores from a solution from `EFA()`. This is just a wrapper for the `psych::factor.scores` function.
* Added function `SCREE()` that does a scree plot. Also added respective print and plot
methods.


## Changes to Functions

* `CD()`: Added check for whether entered data is a tibble, and if so, convert to vanilla data.frame to avoid breaking the procedure.
* `EFA()`: 
    * Updated the EFAtools type in PAF and Promax.
    * Added p value for chi square value in output (calculated for ML and ULS fitting methods).
    * Updated the SPSS varimax implementation to fit SPSS results more closely.
    * Created an argument "varimax_type" that is set according to the specified type, but that can also be specified individually. With type R psych and EFAtools, the stats::varimax is called by default (`varimax_type = "svd"`), with type SPSS, the reproduced SPSS varimax implementation is used (`varimax_type = "kaiser"`).
    * Renamed the `kaiser` argument (controls if a Kaiser normalization is done or not) into `normalize` to avoid confusion with the `varimax_type` argument specifications.
* `ML()`: Changed default start method to "psych".
* `N_FACTORS()`:
    * Added option to do a scree plot if "SCREE" is included in the `criteria` argument.
    * Added a progress bar.
* `OMEGA()`: Now also works with a lavaan second-order solution as input. In this case, it does a Schmid-Leiman transformation based on the first- and second-order loadings first and computes omegas based on this Schmid-Leiman solution.
* `SL()`: Now also works with a lavaan second-order solution as input (first- and second-order loadings taken directly from lavaan output).


## Bug Fixes

* `.get_compare_matrix()`: Fixed a bug that occurred when names of data were longer than n_char
* `COMPARE()`: Fixed a bug that occurred when using `reorder = "names"`.
* `EFA()`: RMSEA is now set to 1 if it is > 1.
* `HULL()`: Fixed a bug that occurred when no factors are located on the HULL
* `KMO()`: Fixed a bug that the inverse of the correlation matrix was not taken anew after smoothing was necessary.
* `PARALLEL()`:
    * Fixed a bug that occurred when using `decision_rule = "percentile"`
    * Relocated error messages that were not evaluated if no data were entered (and should be)
* `print.COMPARE()`: Fixed a bug that occurred when using `print_diff = FALSE` in `COMPARE()`.
* `print.KMO()`: Fixed a bug that printed two statements instead of one, when the KMO value was < .6.

## Minor Changes
* `OMEGA()` and `SL()`: Added an error message if the entered term in `g_name` is invalid (i.e., it cannnot be found among the factor names of the entered lavaan solution).


# EFAtools 0.1.1

## Minor Changes

* Added an error message in `PARALLEL()` if no solution has been found after 25 tries.

## Bug Fixes

* Updated different tests

* Deleted no longer used packages from Imports and Suggests in DESCRIPTION

* `PARALLEL()`: fixed a bug in indexing if method `"EFA"` was used.


# EFAtools 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Initial CRAN submission
## Resubmission
This is a resubmission. In this version we have:

* Updated the type "psych" in the EFA() and EFA_AVERAGE() function according to recent changes in the psych package
* Updated the default type "EFAtools" in the EFA() and EFA_AVERAGE() functions
* Small bug fixes

## Test environments
* local OS X R installation (Mojave 10.14.6), R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 0 notes
---
title: 'EFAtools: An R package with fast and flexible implementations of exploratory
  factor analysis tools'
tags:
- R
- exploratory factor analysis
- factor retention methods
- hierarchical factor analysis
- comparison of implementations
date: "19 July 2020"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: Markus D. Steiner
  orcid: 0000-0002-8126-0757
  affiliation: 1
- name: Silvia Grieder
  orcid: 0000-0002-0118-7722
  affiliation: 2
bibliography: paper.bib
affiliations:
- name: Center for Cognitive and Decision Sciences, Department of Psychology, University
    of Basel, Switzerland
  index: 1
- name: Division of Developmental and Personality Psychology, Department of Psychology,
    University of Basel, Switzerland
  index: 2
---

# Summary

In the social sciences, factor analysis is a widely used tool to identify latent constructs underlying task performance or the answers to questionnaire items. Exploratory factor analysis (EFA) is a data-driven approach to factor analysis and is used to extract a smaller number of common factors that represent or explain the common variance of a larger set of manifest variables [see, e.g., @watkins2018exploratory for an overview]. Several decisions have to be made in advance when performing an EFA, including the number of factors to extract, and the extraction and rotation method to be used. After a factor solution has been found, it is useful to subject the resulting factor solution to an orthogonalization procedure to achieve a hierarchical factor solution with one general and several specific factors. This situation especially applies to data structures in the field of intelligence research where usually high, positive factor intercorrelations occur. From this orthogonalized, hierarchical solution, the variance can then be partitioned to estimate the relative importance of the general versus the specific factors using omega reliability coefficients [e.g., @mcdonald1999test].

*EFAtools* is an R package [@R2018R] that enables fast and flexible analyses in an EFA framework, from tests for suitability of the data for factor analysis and factor retention criteria to hierarchical factor analysis with Schmid-Leiman transformation [@schmid1957development] and McDonald's omegas [e.g., @mcdonald1999test]. The package's core functionalities are listed in Table 1. 

# Statement of Need

Compared to other R packages with which EFA can be performed, *EFAtools* has several advantages, including fast implementations using *Rcpp* [@eddelbuettel2017extending; @eddelbuettel2014rcpparmadillo], more flexibility in the adjustment of implementation features, the ability to reproduce the R *psych* [@revelle2019psych] and SPSS [@ibm_spss_2015] implementations of some analyses methods (see vignette *Replicate SPSS and R psych results with EFAtools*), as well as the inclusion of recommended implementations for these methods based on simulation analyses [@grieder_algorithmic_2020]. Finally, the package includes the implementation of the, as of yet, most comprehensive set of factor retention criteria in R, including recently developed criteria such as the Hull method [@lorenzoseva_hull_2011], comparison data [@ruscio_determining_2012], and the empirical Kaiser criterion [@braeken_empirical_2016]. As recommended by @auerswald_how_2019, multiple factor retention criteria should be examined simultaneously to check their convergence, which now is easily possible with a comprehensive function in *EFAtools* incorporating all implemented factor retention criteria for simultaneous application. Minor advantages over and above the existing implementations in R include that when intending to perform a Schmid-Leiman transformation, this can be done on an obliquely rotated solution obtained with functions from the *EFAtools* or the *psych* package instead of being forced to perform the whole EFA procedure again. Moreover, our implementation of McDonald's omegas calculations include the possibility of manual variable-to-factor correspondences (as are needed for variance partitioning for predetermined / theoretical composites) in addition to automatically determined variable-to-factor correspondences (as done, for example, in the *psych* package). Further, the *EFAtools* function to compute McDonald's omegas can easily be applied on *EFAtools* and *psych* Schmid-Leiman solutions as well as on *lavaan* [@rosseel_lavaan_2012] second-order, bifactor, and single factor solutions (including solutions from multiple group analyses).

# Development and Purpose

*EFAtools* was designed for use in the social sciences in general and is especially suitable for research on cognitive abilities or other hierarchically organized constructs as well as for more time-consuming applications such as in simulation analyses. Its development arose from the need for a tool for easy replication and comparison of EFA solutions from different programs, namely R and SPSS [@grieder_algorithmic_2020], and has already been used in another publication [@grieder_exploratory_2019]. The package was then expanded for a broader, easy, fast, and flexible use of EFA tools such that it is now suitable for most projects within the EFA framework.


Table: Core functionalities of *EFAtools*.

| Topic                | Method                     | Function   
|----------------------|----------------------------|---------------|
|Suitability for factor analysis | Bartlett's test of sphericity | `BARTLETT()` |
|                                | Kaiser-Meyer-Olkin criterion | `KMO()` |
|Factor retention criteria | Comparison data                    | `CD()` |
|                          | Empirical Kaiser criterion         | `EKC()` |
|                          | Hull method                        | `HULL()` |
|                          | Kaiser-Guttman criterion           | `KGC()` |
|                          | Parallel analysis                  | `PARALLEL()` |
|                          | Scree plot                         | `SCREE()` |
|                          | Sequential model tests             | `SMT()` |
|                          | RMSEA lower bound criterion        | `SMT()` |
|                          | AIC criterion                      | `SMT()` |
|Factor extraction methods | Principal axis factoring           | `EFA()` |
|                          | Maximum likelihood                 | `EFA()` |
|                          | Unweighted least squares           | `EFA()` |
|Rotation methods | Orthogonal: Varimax, equamax, quartimax, geominT, bentlerT, bifactorT | `EFA()` |
|                 | Oblique: Promax, oblimin, quartimin, simplimax, bentlerQ, geominQ, bifactorQ | `EFA()` |
|Factor scores             | Different methods for calculating factor scores           | `FACTOR_SCORES()` |
|Hierarchical factor analysis | Schmid-Leiman transformation   | `SL()` |
|                          | McDonald's omegas                 | `OMEGA()` |
*Note*. All functions for suitability for factor analysis and factor retention criteria can be called in any desired combination using the `N_FACTORS()` function.


# Installation

The *EFAtools* package can be installed from CRAN using `install.packages("EFAtools")`. Moreover, the development version can be installed from GitHub (https://github.com/mdsteiner/EFAtools) using `devtools::install_github("mdsteiner/EFAtools", build_vignettes = TRUE)`.

# Acknowledgements

We thank Dirk Wulff for helpful suggestions concerning the C++ implementations.

# References
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# EFAtools

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/EFAtools)](https://CRAN.R-project.org/package=EFAtools)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02521/status.svg)](https://doi.org/10.21105/joss.02521)
[![R-CMD-check](https://github.com/mdsteiner/EFAdiff/workflows/R-CMD-check/badge.svg)](https://github.com/mdsteiner/EFAdiff/actions)
<!-- badges: end -->

The EFAtools package provides functions to perform exploratory factor analysis (EFA) procedures and compare their solutions. The goal is to provide state-of-the-art factor retention methods and a high degree of flexibility in the EFA procedures. This way, implementations from R psych and SPSS can be compared. Moreover, functions for Schmid-Leiman transformation, and computation of omegas are provided. To speed up the analyses, some of the iterative procedures like principal axis factoring (PAF) are implemented in C++.

## Installation

You can install the release version from [CRAN](https://cran.r-project.org/) with:
```{r eval=FALSE}
install.packages("EFAtools")
```

You can install the development version from [GitHub](https://github.com/) with:

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("mdsteiner/EFAtools")
```

To also build the vignette when installing the development version, use:

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("mdsteiner/EFAtools", build_vignettes = TRUE)
```

## Example

Here are a few examples on how to perform the analyses with the different types and how to compare the results using the `COMPARE` function. For more details, see the vignette by running `vignette("EFAtools", package = "EFAtools")`. The vignette provides a high-level introduction into the functionalities of the package.

```{r fig.height=3}
# load the package
library(EFAtools)

# Run all possible factor retention methods
N_FACTORS(test_models$baseline$cormat, N = 500, method = "ML")

# A type SPSS EFA to mimick the SPSS implementation with
# promax rotation
EFA_SPSS <- EFA(test_models$baseline$cormat, n_factors = 3, type = "SPSS",
                  rotation = "promax")

# look at solution
EFA_SPSS

# A type psych EFA to mimick the psych::fa() implementation with
# promax rotation
EFA_psych <- EFA(test_models$baseline$cormat, n_factors = 3, type = "psych",
                  rotation = "promax")

# compare the type psych and type SPSS implementations
COMPARE(EFA_SPSS$rot_loadings, EFA_psych$rot_loadings,
        x_labels = c("SPSS", "psych"))

# Average solution across many different EFAs with oblique rotations
EFA_AV <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = c("PAF", "ML", "ULS"), rotation = "oblique",
                      show_progress = FALSE)

# look at solution
EFA_AV

# Perform a Schmid-Leiman transformation
SL <- SL(EFA_psych)

# Based on a specific salience threshold for the loadings (here: .20):
factor_corres <- SL$sl[, c("F1", "F2", "F3")] >= .2

# Compute omegas from the Schmid-Leiman solution
OMEGA(SL, factor_corres = factor_corres)
```

## Citation

If you use this package in your research, please acknowledge it by citing: 

Steiner, M.D., & Grieder, S.G. (2020). EFAtools: An R package with fast and flexible implementations of exploratory factor analysis tools. *Journal of Open Source Software*, *5*(53), 2521. https://doi.org/10.21105/joss.02521

## Contribute or Report Bugs

If you want to contribute or report bugs, please open an issue on GitHub or email us at markus.d.steiner@gmail.com or silvia.grieder@gmail.com.
---
title: "EFAtools"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{EFAtools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.align = "center"
)

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
      stop("Package \"microbenchmark\" needed for this vignette to work. Please install it.",
      call. = FALSE)
}

```



This vignette provides an overview for the functionalities of the EFAtools package. The general aim of the package is to provide flexible implementations of different algorithms for an exploratory factor analyses (EFA) procedure, including factor retention methods, factor extraction and rotation methods, as well as the computation of a Schmid-Leiman solution and McDonald's omega coefficients.

The package was first designed to enable a comparison of EFA (specifically, principal axis factoring with subsequent promax rotation) performed in R using the [**psych**](https://CRAN.R-project.org/package=psych) package and EFA performed in SPSS. That is why some functions allow the specification of a type, including `"psych"` and `"SPSS"`, such that the respective procedure will be executed to match the output of these implementations (which do not always lead to the same results; see separate vignette [**Replicate_SPSS_psych**](Replicate_SPSS_psych.html "Replicate SPSS and R psych results with EFAtools") for a demonstration of the replication of original results). This vignette will go through a complete example, that is, we will first show how to determine the number of factors to retain, then perform different factor extraction methods, run a Schmid-Leiman transformation and compute omegas.

The package can be installed from CRAN using `install.packages("EFAtools")`, or from GitHub using `devtools::install_github("mdsteiner/EFAtools")`, and then loaded using:

```{r}
library(EFAtools)
```

In this vignette, we will use the `DOSPERT_raw` dataset, which contains responses to the Domain Specific Risk Taking Scale (DOSPERT) of 3123 participants. The dataset is contained in the `EFAtools` package, for details, see `?DOSPERT_raw`. Note that this vignette is to provide a general overview and it is beyond its scope to explain all methods and functions in detail. If you want to learn more on the details and methods, please see the respective help functions for explanations and literature references. However, the dataset is rather large, so, just to save time when building the vignette, we will only use the first 500 observations. When you normally do your analyses, you use the full dataset.

```{r}
# only use a subset to make analyses faster
DOSPERT_sub <- DOSPERT_raw[1:500,]
```


## Test Suitability of Data

The first step in an EFA procedure is to test whether your data is suitable for factor analysis. To this end, the `EFAtools` package provides the `BARTLETT()` and the `KMO()` functions. The Bartlett's test of sphericity tests whether a correlation matrix is significantly different from an identity matrix (a correlation matrix with zero correlations between all variables). This test should thus be significant. The Kaiser-Meyer-Olkin criterion (KMO) represents the degree to which each observed variable is predicted by the other variables in the dataset and thus is another indicator for how correlated the different variables are.

We can test whether our `DOSPERT_sub` dataset is suitable for factor analysis as follows.

```{r}
# Bartlett's test of sphericity
BARTLETT(DOSPERT_sub)

# KMO criterion
KMO(DOSPERT_sub)
```

Note that these tests can also be run in the `N_FACTORS()` function.

## Factor Retention Methods

As the goal of EFA is to determine the underlying factors from a set of multiple variables, one of the most important decisions is how many factors can or should be extracted. There exists a plethora of factor retention methods to use for this decision. The problem is that there is no method that consistently outperforms all other methods. Rather, which factor retention method to use depends on the structure of the data: are there few or many indicators, are factors strong or weak, are the factor intercorrelations weak or strong. For rules on which methods to use, see, for example, [Auerswald and Moshagen, (2019)](https://doi.apa.org/doiLanding?doi=10.1037%2Fmet0000200).

There are multiple factor retention methods implemented in the `EFAtools` package. They can either be called with separate functions, or all (or a selection) of them using the `N_FACTORS()` function. 

### Calling Separate Functions

Let's first look at how to determine the number of factors to retain by calling separate functions. For example, if you would like to perform a parallel analysis based on squared multiple correlations (SMC; sometimes also called a parallel analysis with principal factors), you can do the following:

```{r}
# determine the number of factors to retain using parallel analysis
PARALLEL(DOSPERT_sub, eigen_type = "SMC")
```

Generating the plot can also be suppressed if the output is printed explicitly:

```{r}
# determine the number of factors to retain using parallel analysis
print(PARALLEL(DOSPERT_sub, eigen_type = "SMC"), plot = FALSE)
```


Other factor retention methods can be used accordingly. For example, to use the empirical Kaiser criterion, use the `EKC` function:

```{r}
# determine the number of factors to retain using parallel analysis
print(EKC(DOSPERT_sub), plot = FALSE)
```

The following factor retention methods are currently implemented: comparison data (`CD()`), empirical Kaiser criterion (`EKC()`), the hull method (`HULL()`), the Kaiser-Guttman criterion (`KGC()`), parallel analysis (`PARALLEL()`), scree test (`SCREE()`), and sequential model tests (`SMT()`). Many of these functions have multiple versions of the respective factor retention method implemented, for example, the parallel analysis can be done based on eigenvalues found using unity (principal components) or SMCs, or on an EFA procedure. Another example is the hull method, which can be used with different fitting methods (principal axis factoring [PAF], maximum likelihood [ML], or unweighted least squares [ULS]), and different goodness of fit indices. Please see the respective function documentations for details.

### Run Multiple Factor Retention Methods With `N_FACTORS()`

If you want to use multiple factor retention methods, for example, to compare whether different methods suggest the same number of factors, it is easier to use the `N_FACTORS()` function. This is a wrapper around all the implemented factor retention methods. Moreover, it also enables to run the Bartlett's test of sphericity and compute the KMO criterion.

For example, to test the suitability of the data for factor analysis and to determine the number of factors to retain based on parallel analysis (but only using eigen values based on SMCs and PCA), the EKC, and the sequential model test, we can run the following code:

```{r}
N_FACTORS(DOSPERT_sub, criteria = c("PARALLEL", "EKC", "SMT"),
          eigen_type_other = c("SMC", "PCA"))
```


If all possible factor retention methods should be used, it is sufficient to provide the data object (note that this takes a while, as the comparison data is computationally expensive and therefore relatively slow method, especially if larger datasets are used). We additionally specify the method argument to use unweighted least squares (ULS) estimation. This is a bit faster than using principle axis factoring (PAF) and it enables the computation of more goodness of fit indices:

```{r}
N_FACTORS(DOSPERT_sub, method = "ULS")
```

Now, this is not the scenario one is happy about, but it still does happen: There is no obvious convergence between the methods and thus the choice of the number of factors to retain becomes rather difficult (and to some extend arbitrary). We will proceed with 6 factors, as it is what is typically used with DOSPERT data, but this does not mean that other number of factors are not just as plausible.

Note that all factor retention methods, except comparison data (CD), can also be used with correlation matrices. We use `method = "ULS"` and `eigen_type_other = c("SMC", "PCA")` to skip the slower criteria. In this case, the sample size has to be specified:

```{r}
N_FACTORS(test_models$baseline$cormat, N = 500,
          method = "ULS", eigen_type_other = c("SMC", "PCA"))
```

## Exploratory Factor Analysis: Factor Extraction

Multiple algorithms to perform an EFA and to rotate the found solutions are implemented in the `EFAtools` package. All of them can be used using the `EFA()` function. To perform the EFA, you can use one of principal axis factoring (PAF), maximum likelihood estimation (ML), and unweighted least squares (ULS; also sometimes referred to as MINRES). To rotate the solutions, the `EFAtools` package offers varimax and promax rotations, as well as the orthogonal and oblique rotations provided by the `GPArotation` package (i.e., the `GPArotation` functions are called in the `EFA()` function in this case).

You can run an EFA with PAF and no rotation like this:

```{r}
EFA(DOSPERT_sub, n_factors = 6)
```

To rotate the loadings (e.g., using a promax rotation) adapt the `rotation` argument:

```{r}
EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")
```

This now performed PAF with promax rotation with the specification, on average, we found to produce the most accurate results in a simulation analysis (see function documentation). If you want to replicate the implementation of the *psych* R package, you can set the `type` argument to `"psych"`:

```{r}
EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "psych")
```

If you want to use the *SPSS* implementation, you can set the `type` argument to `"SPSS"`: 

```{r}
EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "SPSS")
```

This enables comparisons of different implementations. The `COMPARE()` function provides an easy way to compare how similar two loading (pattern) matrices are:

```{r}
COMPARE(
  EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "psych")$rot_loadings,
  EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "SPSS")$rot_loadings
)

```

*Why would you want to do this?* One of us has had the experience that a reviewer asked whether the results can be reproduced in another statistical program than R. We therefore implemented this possibility in the package for an easy application of large scale, systematic comparisons.

Note that the `type` argument of the `EFA()` function only affects the implementations of principal axis factoring (PAF), varimax and promax rotations. The other procedures are not affected (except the order of the rotated factors for the other rotation methods).

As indicated previously, it is also possible to use different estimation and rotation methods. For example, to perform an EFA with ULS and an oblimin rotation, you can use the following code:

```{r}
EFA(DOSPERT_sub, n_factors = 6, rotation = "oblimin", method = "ULS")
```

Of course, `COMPARE()` can also be used to compare results from different estimation or rotation methods (in fact, to compare any two matrices), not just from different implementations:

```{r}
COMPARE(
  EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")$rot_loadings,
  EFA(DOSPERT_sub, n_factors = 6, rotation = "oblimin", method = "ULS")$rot_loadings,
  x_labels = c("PAF and promax", "ULS and oblimin")
)
```

Finally, if you are interested in factor scores from the EFA solution, these can be obtained with `FACTOR_SCORES()`, a wrapper for `psych::factor.scores()` to be used directly with an output from `EFA()`:

```{r}
EFA_mod <- EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")
fac_scores <- FACTOR_SCORES(DOSPERT_sub, f = EFA_mod)
```


### Performance

To improve performance of the iterative procedures (currently the parallel analysis, and the PAF, ML, and ULS methods) we implemented some of them in C++. For example, the following code compares the EFAtools parallel analysis with the corresponding one implemented in the psych package (the default of `PARALLEL()` is to use 1000 datasets, but 25 is enough to show the difference):

```{r message=FALSE}
microbenchmark::microbenchmark(
  PARALLEL(DOSPERT_sub, eigen_type = "SMC", n_datasets = 25),
  psych::fa.parallel(DOSPERT_sub, SMC = TRUE, plot = FALSE, n.iter = 25)
)
```


Moreover, the following code compares the PAF implementation (of type "psych") of the EFAtools package with the one from the psych package: 

```{r message=FALSE}
microbenchmark::microbenchmark(
  EFA(DOSPERT_raw, 6),
  psych::fa(DOSPERT_raw, 6, rotate = "none", fm = "pa")
)
```

While these differences are not large, they grow larger the more iterations the procedures need, which is usually the case if solutions are more tricky to find. Especially for simulations this might come in handy. For example, in one simulation analysis we ran over 10,000,000 EFAs, thus a difference of about 25 milliseconds per EFA leads to a difference in runtime of almost three days.

### Model Averaging

Instead of relying on one of the many possible implementations of, for example, PAF, and of using just one rotation (e.g., promax), it may be desirable to average different solutions to potentially arrive at a more robust, average solution. The `EFA_AVERAGE()` function provides this possibility. In addition to the average solution it provides the variation across solutions, a matrix indicating the robustness of indicator-to-factor correspondences, and a visualisation of the average solution and the variability across solutions. For example, to average across all available factor extraction methods and across all available oblique rotations, the following code can be run:

```{r}
# Average solution across many different EFAs with oblique rotations
EFA_AV <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = c("PAF", "ML", "ULS"), rotation = "oblique",
                      show_progress = FALSE)

# look at solution
EFA_AV
```


The first matrix of the output tells us that the indicators are mostly allocated to the same factors. However, that some rowsums are larger than one also tells as that there likely are some cross loadings present in some solutions. Moreover, the relatively high percentages of salient pattern coefficients all loading on the first factor may indicate that some rotation methods failed to achieve simple structure and it might be desirable to exclude these from the model averaging procedure. The rest of the output is similar to the normal `EFA()` outputs shown above, only that in addition to the average coefficients their range is also shown. Finally, the plot shows the average pattern coefficients and their ranges. 

**Important disclaimer:** While it is possible that this approach provides more robust results, we are unaware of simulation studies that have investigated and shown this. Therefore, it might make sense to for now use this approach mainly to test the robustness of the results obtained with one single EFA implementation.

## Exploratory Factor Analysis: Schmid-Leiman transformation and McDonald's Omegas

For the Schmid-Leiman transformation and computation of omegas, we will use PAF and promax rotation:

```{r}
efa_dospert <- EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")
efa_dospert
```

The indicator names in the output (i.e., the rownames of the rotated loadings section) tell us which domain (out of ethical, financial, health, recreational, and social risks) an indicator stems from. From the pattern coefficients it can be seen that these theoretical domains are recovered relatively well in the six factor solution, that is, usually, the indicators from the same domain load onto the same factor. When we take a look at the factor intercorrelations, we can see that there are some strong and some weak correlations. It might be worthwhile to explore whether a general factor can be obtained, and which factors load more strongly on it. To this end, we will use a Schmid-Leiman (SL) transformation.

## Schmid-Leiman Transformation

The SL transformation or orthogonalization transforms an oblique solution into a hierarchical, orthogonalized solution. To do this, the `EFAtools` package provides the `SL()` function.

```{r}
sl_dospert <- SL(efa_dospert)
sl_dospert
```

From the output, it can be seen that all, except the social domain indicators substantially load on the general factor. That is, the other domains covary substantially.

## McDonald's Omegas

Finally, we can compute omega estimates and additional indices of interpretive relevance based on the SL solution. To this end, we can either specify the variable-to-factor correspondences, or let them be determined automatically (in which case the highest factor loading will be taken, which might lead to a different solution than what is desired, in the presence of cross-loadings). Given that no cross-loadings are present here, it is easiest to let the function automatically determine the variable-to-factor correspondence. To this end, we will set the `type` argument to `"psych"`.

```{r}
OMEGA(sl_dospert, type = "psych")
```

If we wanted to specify the variable to factor correspondences explicitly (for example, according to theoretical expectations), we could do it in the following way:

```{r}
OMEGA(sl_dospert, factor_corres = matrix(c(rep(0, 18), rep(1, 6), rep(0, 30), 
                                         rep(1, 6), rep(0, 6), 1, 0, 1, 0, 1,
                                         rep(0, 19), rep(1, 6), rep(0, 31), 1, 0,
                                         1, 0, 1, rep(0, 30), rep(1, 6), 
                                         rep(0, 12)), ncol = 6, byrow = FALSE))
```


---
title: "Replicate SPSS and R psych results with EFAtools"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Replicate_SPSS_psych}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.align = "center"
)


if (!requireNamespace("psych", quietly = TRUE)) {
      stop("Package \"psych\" needed for this vignette to work. Please install it.",
      call. = FALSE)
}

```

This vignette demonstrates the replication of exploratory factor analysis (EFA) results (specifically, principal axis factoring [PAF] with subsequent promax rotation) from the SPSS `FACTOR` algorithm and from the `fa()` function from the `psych` R package.
For a general introduction to the `EFAtools` package, please see the [**EFAtools**](EFAtools.html "EFAtools") vignette. Same as in the EFAtools vignette, we will use the DOSPERT data set for this demonstration as well (see `?DOSPERT` for details).

First load the needed packages EFAtools and psych (original SPSS results for some data sets are available in the EFAtools package).

```{r setup}
library(psych)
library(EFAtools)
```

## Principal Axis Factoring

First, we will fit an EFA with PAF and without rotation using the `EFA` function from `EFAtools` using `type = "psych"` and `type = "SPSS"`. These types are intended to mimic the R psych and SPSS results, respectively.

```{r}
# EFAtools::EFA with type = "psych" without rotation
EFA_psych_paf <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                     type = "psych")
# EFAtools::EFA with type = "SPSS" without rotation
EFA_SPSS_paf <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                    type = "SPSS")
```

As a next step, we fit an EFA with the same configurations (PAF and no rotation) using the `fa` function from `psych` with the same data set.

```{r}
# psych::fa without rotation
psych_paf <- psych::fa(DOSPERT$cormat, nfactors = 10, n.obs = DOSPERT$N,
                       fm = "pa", rotate = "none")
```

Now we can compare results from `EFA` with the respective types to the original R `psych` and `SPSS` results using the same data set. This is easily done using the `COMPARE` function available in the `EFAtools` package.

```{r}
# Compare loadings from psych::fa and EFAtools::EFA with type = "psych"
COMPARE(EFA_psych_paf$unrot_loadings, psych_paf$loadings)

# Compare loadings from SPSS and EFAtools::EFA with type = "SPSS"
COMPARE(EFA_SPSS_paf$unrot_loadings, SPSS_27$DOSPERT$paf_load)

```

To see that this close match was not trivial, we can look at the match between the original R `psych` and `SPSS` solutions.

```{r}
# Compare loadings from psych::fa and SPSS
COMPARE(psych_paf$loadings, SPSS_27$DOSPERT$paf_load)
```

We can see that the solutions are slightly different, especially for the 9th and 10th factor. Although the differences are very small here, they can get quite large for other data sets, or get larger after rotation (see below).

## Varimax Rotation

Now we confirmed the replication of PAF results without rotation, we can continue to compare rotated factor solutions. We start by comparing varimax rotated PAF solutions.


```{r}
## Fit the models

# EFAtools::EFA with type = "psych" with varimax rotation
EFA_psych_var <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                     type = "psych", rotation = "varimax")
# EFAtools::EFA with type = "SPSS" with varimax rotation
EFA_SPSS_var <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                    type = "SPSS", rotation = "varimax")
# psych::fa with varimax rotation
psych_var <- psych::fa(DOSPERT$cormat, nfactors = 10, n.obs = DOSPERT$N,
                       fm = "pa", rotate = "varimax")

## Check replication of results

# Compare loadings from psych::fa and EFAtools::EFA with type = "psych"
COMPARE(EFA_psych_var$rot_loadings, psych_var$loadings)

# Compare loadings from SPSS and EFAtools::EFA with type = "SPSS"
COMPARE(EFA_SPSS_var$rot_loadings, SPSS_27$DOSPERT$var_load)

## Compare original results (just to see the difference)

# Compare loadings from psych::fa and SPSS
COMPARE(psych_var$loadings, SPSS_27$DOSPERT$var_load)

```

## Promax Rotation

Finally, we can do the same for promax rotated results as well.

```{r}
## Fit the models

# EFAtools::EFA with type = "psych" with promax rotation
EFA_psych_pro <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                     type = "psych", rotation = "promax")
# EFAtools::EFA with type = "SPSS" with promax rotation
EFA_SPSS_pro <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                    type = "SPSS", rotation = "promax")
# psych::fa with promax rotation
psych_pro <- psych::fa(DOSPERT$cormat, nfactors = 10, n.obs = DOSPERT$N,
                       fm = "pa", rotate = "Promax")

## Check replication of results

# Compare loadings from psych::fa and EFAtools::EFA with type = "psych"
COMPARE(EFA_psych_pro$rot_loadings, psych_pro$loadings)

# Compare loadings from SPSS and EFAtools::EFA with type = "SPSS"
COMPARE(EFA_SPSS_pro$rot_loadings, SPSS_27$DOSPERT$pro_load)

## Compare original results (just to see the difference)

# Compare loadings from psych::fa and SPSS
COMPARE(psych_pro$loadings, SPSS_27$DOSPERT$pro_load)

```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.PARALLEL.R
\name{print.PARALLEL}
\alias{print.PARALLEL}
\title{Print function for PARALLEL objects}
\usage{
\method{print}{PARALLEL}(x, plot = TRUE, ...)
}
\arguments{
\item{x}{a list of class PARALLEL. Output from \link{PARALLEL} function.}

\item{plot}{logical. Whether to plot the results.}

\item{...}{Further arguments for print.}
}
\description{
Print function for PARALLEL objects
}
\examples{
\donttest{
# example without real data
PARALLEL(N = 500, n_vars = 10)

# example with correlation matrix and "ML" estimation
PARALLEL(test_models$case_11b$cormat, N = 500, method = "ML")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KMO.R
\name{KMO}
\alias{KMO}
\title{Kaiser-Meyer-Olkin criterion}
\source{
Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
35, 401-415.

Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational
and Psychological Measurement, 34, 111-117.

Cureton, E. E. & D'Augustino, R. B. (1983). Factor analysis: An
 applied approach. Hillsdale, N.J.: Lawrence Erlbaum Associates, Inc.
}
\usage{
KMO(
  x,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall")
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
data is given as input. Default is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}
}
\value{
A list containing
\item{KMO}{Overall KMO.}
\item{KMO_i}{KMO for each variable.}
\item{settings}{A list of the settings used.}
}
\description{
This function computes the Kaiser-Meyer-Olkin (KMO) criterion overall and for
each variable in a correlation matrix. The KMO represents the degree to
which each observed variable is predicted by the other variables in the
dataset and with this indicates the suitability for factor analysis.
}
\details{
Kaiser (1970) proposed this index, originally called measure of
sampling adequacy (MSA), that indicates how near the inverted correlation
matrix \eqn{R^{-1}} is to a diagonal matrix \eqn{S} to determine a given
correlation matrix's (\eqn{R}) suitability for factor analysis.
The index is
\deqn{KMO = \frac{\sum\limits_{i<j}\sum r_{ij}^2}{\sum\limits_{i<j}\sum r_{ij}^2 + \sum\limits_{i<j}\sum q_{ij}^2}}
with \eqn{Q = SR^{-1}S} and S = \eqn{(diag R^{-1})^{-1/2}} where
\eqn{\sum\limits_{i<j}\sum r_{ij}^2} is the sum of squares of the upper
off-diagonal elements of \eqn{R} and \eqn{\sum\limits_{i<j}\sum q_{ij}^2} is the
sum of squares of the upper off-diagonal elements of \eqn{Q} (see also Cureton & D'Augustino, 1983).

So KMO varies between 0 and 1, with larger values indicating higher suitability
for factor analysis. Kaiser and Rice (1974) suggest that KMO should at least
exceed .50 for a correlation matrix to be suitable for factor analysis.

This function was heavily influenced by the \code{\link[psych:KMO]{psych::KMO}}
function.

See also \code{\link{BARTLETT}} for another test of suitability for factor
analysis.

The \code{KMO} function can also be called together with the
\code{\link{BARTLETT}} function and with factor retention criteria in the
 \code{\link{N_FACTORS}} function.
}
\examples{
KMO(test_models$baseline$cormat)
}
\seealso{
\code{\link{BARTLETT}} for another measure to determine
suitability for factor analysis.

\code{\link{N_FACTORS}} as a wrapper function for this function,
\code{\link{BARTLETT}} and several factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BARTLETT.R
\name{BARTLETT}
\alias{BARTLETT}
\title{Bartlett's test of sphericity}
\source{
Bartlett, M. S. (1951). The effect of standardization on a Chi-square
approximation in factor analysis. Biometrika, 38, 337-344.
}
\usage{
BARTLETT(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall")
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations.}

\item{N}{numeric. The number of observations. Needs only be specified if a
correlation matrix is used.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
is given as input. Default is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}
}
\value{
A list containing
\item{chisq}{The chi square statistic.}
\item{p_value}{The p value of the chi square statistic.}
\item{df}{The degrees of freedom for the chi square statistic.}
\item{settings}{A list of the settings used.}
}
\description{
This function tests whether a correlation matrix is significantly different
from an identity matrix (Bartlett, 1951). If the Bartlett's test is not
significant, the correlation matrix is not suitable for factor analysis
because the variables show too little covariance.
}
\details{
Bartlett (1951) proposed this statistic to determine a correlation
matrix' suitability for factor analysis. The statistic is approximately
chi square distributed with \eqn{df = \frac{p(p - 1)}{2}} and is given by

\deqn{chi^2 = -log(det(R)) (N - 1 - (2 * p + 5)/6)}

where \eqn{det(R)} is the determinant of the correlation matrix, \eqn{N} is
the sample size, and \eqn{p} is the number of variables.

This tests requires multivariate normality. If this condition is not met,
the Kaiser-Meyer-Olkin criterion (\code{\link[EFAtools]{KMO}})
can still be used.

This function was heavily influenced by the \code{\link[psych:cortest.bartlett]{psych::cortest.bartlett}} function from the psych package.

The \code{BARTLETT} function can also be called together with the
 (\code{\link[EFAtools]{KMO}}) function and with factor retention criteria
 in the \code{\link{N_FACTORS}} function.
}
\examples{
BARTLETT(test_models$baseline$cormat, N = 500)

}
\seealso{
\code{\link[EFAtools]{KMO}} for another measure to determine
 suitability for factor analysis.

 \code{\link{N_FACTORS}} as a wrapper function for this function,
 \code{\link[EFAtools]{KMO}} and several factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.LOADINGS.R
\name{print.LOADINGS}
\alias{print.LOADINGS}
\title{Print LOADINGS object}
\usage{
\method{print}{LOADINGS}(x, cutoff = 0.3, digits = 3, ...)
}
\arguments{
\item{x}{class LOADINGS matrix.}

\item{cutoff}{numeric. The number above which to print loadings in bold
default is .3.}

\item{digits}{numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
to round the loadings to (default is 3).}

\item{...}{additional arguments passed to print}
}
\description{
Print LOADINGS object
}
\examples{
EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                    type = "EFAtools", method = "PAF", rotation = "promax")
EFAtools_PAF

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMPARE.R
\name{COMPARE}
\alias{COMPARE}
\title{Compare two vectors or matrices (communalities or loadings)}
\usage{
COMPARE(
  x,
  y,
  reorder = c("congruence", "names", "none"),
  corres = TRUE,
  thresh = 0.3,
  digits = 4,
  m_red = 0.001,
  range_red = 0.001,
  round_red = 3,
  print_diff = TRUE,
  na.rm = FALSE,
  x_labels = c("x", "y"),
  plot = TRUE,
  plot_red = 0.01
)
}
\arguments{
\item{x}{matrix, or vector. Loadings or communalities of a factor
analysis output.}

\item{y}{matrix, or vector. Loadings or communalities of another
factor analysis output to compare to x.}

\item{reorder}{character. Whether and how elements / columns should be
reordered. If "congruence" (default), reordering is done according to Tuckers
correspondence coefficient, if "names", objects according to their names,
if "none", no reordering is done.}

\item{corres}{logical. Whether factor correspondences should be compared if a
matrix is entered.}

\item{thresh}{numeric. The threshold to classify a pattern coefficient as substantial. Default is .3.}

\item{digits}{numeric. Number of decimals to print in the output. Default is 4.}

\item{m_red}{numeric. Number above which the mean and median should be printed
in red (i.e., if .001 is used, the mean will be in red if it is larger than
.001, otherwise it will be displayed in green.) Default is .001.}

\item{range_red}{numeric. Number above which the min and max should be printed
in red (i.e., if .001 is used, min and max will be in red if the max is larger
 than .001, otherwise it will be displayed in green. Default is .001). Note that
 the color of min also depends on max, that is min will be displayed in the
 same color as max.}

\item{round_red}{numeric. Number above which the max decimals to round to where
all corresponding elements of x and y are still equal are displayed in red
(i.e., if 3 is used, the number will be in red if it is smaller than
 3, otherwise it will be displayed in green). Default is 3.}

\item{print_diff}{logical. Whether the difference vector or matrix should be
printed or not. Default is TRUE.}

\item{na.rm}{logical. Whether NAs should be removed in the mean, median, min,
and max functions. Default is FALSE.}

\item{x_labels}{character. A vector of length two containing identifying
labels for the two objects x and y that will be compared. These will be used
as labels on the x-axis of the plot. Default is "x" and "y".}

\item{plot}{logical. If TRUE (default), a plot illustrating the differences
will be shown.}

\item{plot_red}{numeric. Threshold above which to plot the absolute differences
in red. Default is .001.}
}
\value{
A list of class COMPARE containing summary statistics on the differences
 of x and y.

\item{diff}{The vector or matrix containing the differences between x and y.}
\item{mean_abs_diff}{The mean absolute difference between x and y.}
\item{median_abs_diff}{The median absolute difference between x and y.}
\item{min_abs_diff}{The minimum absolute difference between x and y.}
\item{max_abs_diff}{The maximum absolute difference between x and y.}
\item{max_dec}{The maximum number of decimals to which a comparison makes sense.
 For example, if x contains only values up to the third decimals, and y is a
 normal double, max_dec will be three.}
\item{are_equal}{The maximal number of decimals to which all elements of x and y
 are equal.}
\item{diff_corres}{The number of differing variable-to-factor correspondences
 between x and y, when only the highest loading is considered.}
\item{diff_corres_cross}{The number of differing variable-to-factor correspondences
 between x and y when all loadings \code{>= thresh} are considered.}
\item{g}{The root mean squared distance (RMSE) between x and y.}
\item{settings}{List of the settings used.}
}
\description{
The function takes two objects of the same dimensions containing numeric
information (loadings or communalities) and returns a list of class COMPARE
containing summary information of the differences of the objects.
}
\examples{
# A type SPSS EFA to mimick the SPSS implementation
EFA_SPSS_6 <- EFA(test_models$case_11b$cormat, n_factors = 6, type = "SPSS")

# A type psych EFA to mimick the psych::fa() implementation
EFA_psych_6 <- EFA(test_models$case_11b$cormat, n_factors = 6, type = "psych")

# compare the two
COMPARE(EFA_SPSS_6$unrot_loadings, EFA_psych_6$unrot_loadings,
        x_labels = c("SPSS", "psych"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HULL.R
\name{HULL}
\alias{HULL}
\title{Hull method for determining the number of factors to retain}
\source{
Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011).
The Hull method for selecting the number of common factors. Multivariate
Behavioral Research, 46(2), 340-364.
}
\usage{
HULL(
  x,
  N = NA,
  n_fac_theor = NA,
  method = c("PAF", "ULS", "ML"),
  gof = c("CAF", "CFI", "RMSEA"),
  eigen_type = c("SMC", "PCA", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall"),
  n_datasets = 1000,
  percent = 95,
  decision_rule = c("means", "percentile", "crawford"),
  n_factors = 1,
  ...
)
}
\arguments{
\item{x}{matrix or data.frame. Dataframe or matrix of raw data or matrix with
correlations.}

\item{N}{numeric. Number of cases in the data. This is passed to \link{PARALLEL}.
Only has to be specified if x is a correlation matrix, otherwise it is determined
based on the dimensions of x.}

\item{n_fac_theor}{numeric. Theoretical number of factors to retain. The maximum
of this number and the number of factors suggested by \link{PARALLEL} plus
one will be used in the Hull method.}

\item{method}{character. The estimation method to use. One of  \code{"PAF"},
\code{"ULS"}, or  \code{"ML"}, for principal axis factoring, unweighted
least squares, and maximum likelihood, respectively.}

\item{gof}{character. The goodness of fit index to use. Either \code{"CAF"},
\code{"CFI"}, or \code{"RMSEA"}, or any combination of them.
If \code{method = "PAF"} is used, only
the CAF can be used as goodness of fit index. For details on the CAF, see
Lorenzo-Seva, Timmerman, and Kiers (2011).}

\item{eigen_type}{character. On what the eigenvalues should be found in the
parallel analysis. Can be one of \code{"SMC"}, \code{"PCA"}, or \code{"EFA"}.
 If using  \code{"SMC"} (default), the diagonal of the correlation matrices is
  replaced by the squared multiple correlations (SMCs) of the indicators. If
   using  \code{"PCA"}, the diagonal values of the correlation
matrices are left to be 1. If using  \code{"EFA"}, eigenvalues are found on the
correlation  matrices with the final communalities of an EFA solution as
diagonal. This is passed to  \code{\link{PARALLEL}}.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
is given as input. Default is \code{"pairwise.complete.obs"}.}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is  \code{"pearson"}.}

\item{n_datasets}{numeric. The number of datasets to simulate. Default is 1000.
This is passed to \code{\link{PARALLEL}}.}

\item{percent}{numeric. A vector of percentiles to take the simulated eigenvalues from.
Default is 95. This is passed to \code{\link{PARALLEL}}.}

\item{decision_rule}{character. Which rule to use to determine the number of
factors to retain. Default is \code{"means"}, which will use the average
simulated eigenvalues. \code{"percentile"}, uses the percentiles specified
in percent. \code{"crawford"} uses the 95th percentile for the first factor
and the mean afterwards (based on Crawford et al, 2010). This is passed to \code{\link{PARALLEL}}.}

\item{n_factors}{numeric. Number of factors to extract if  \code{"EFA"} is
included in \code{eigen_type}. Default is 1. This is passed to
\code{\link{PARALLEL}}.}

\item{...}{Further arguments passed to \code{\link{EFA}}, also in
\code{\link{PARALLEL}}.}
}
\value{
A list of class HULL containing the following objects
\item{n_fac_CAF}{The number of factors to retain according to the Hull method
with the CAF.}
\item{n_fac_CFI}{The number of factors to retain according to the Hull method
with the CFI.}
\item{n_fac_RMSEA}{The number of factors to retain according to the Hull method
with the RMSEA.}
\item{solutions_CAF}{A matrix containing the CAFs, degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmerman, and Kiers 2011 for details).}
\item{solutions_CFI}{A matrix containing the CFIs, degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmerman, and Kiers 2011 for details).}
\item{solutions_RMSEA}{A matrix containing the RMSEAs, degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmerman, and Kiers 2011 for details).}
\item{n_fac_max}{The upper bound \emph{J} of the number of factors to extract (see details).}
\item{settings}{A list of the settings used.}
}
\description{
Implementation of the Hull method suggested by Lorenzo-Seva, Timmerman,
and Kiers (2011), with an extension to principal axis factoring. See details for
parallelization.
}
\details{
The Hull method aims to find a model with an optimal balance between
 model fit and number of parameters. That is, it aims to retrieve only major
 factors (Lorenzo-Seva, Timmerman, & Kiers, 2011). To this end, it performs
 the following steps (Lorenzo-Seva, Timmerman, & Kiers, 2011, p.351):
 \enumerate{
   \item It performs parallel analysis and adds one to the identified number of factors (this number is denoted \emph{J}). \emph{J} is taken as an upper bound of the number of factors to retain in the hull method. Alternatively, a theoretical number of factors can be entered. In this case \emph{J} will be set to whichever of these two numbers (from parallel analysis or based on theory) is higher.
   \item For all 0 to \emph{J} factors, the goodness-of-fit (one of \emph{CAF}, \emph{RMSEA}, or \emph{CFI}) and the degrees of freedom (\emph{df}) are computed.
   \item The solutions are ordered according to their \emph{df}.
   \item Solutions that are not on the boundary of the convex hull are eliminated (see Lorenzo-Seva, Timmerman, & Kiers, 2011, for details).
   \item All the triplets of adjacent solutions are considered consecutively. The middle solution is excluded if its point is below or on the line connecting its neighbors in a plot of the goodness-of-fit versus the degrees of freedom.
   \item Step 5 is repeated until no solution can be excluded.
   \item The \emph{st} values of the “hull” solutions are determined.
   \item The solution with the highest \emph{st} value is selected.
 }

The \link{PARALLEL} function and the principal axis factoring of the
  different number of factors can be parallelized using the future framework,
  by calling the \link[future:plan]{future::plan} function. The examples
   provide example code on how to enable parallel processing.

  Note that if \code{gof = "RMSEA"} is used, 1 - RMSEA is actually used to
  compare the different solutions. Thus, the threshold of .05 is then .95. This
  is necessary due to how the heuristic to locate the elbow of the hull works.

  The ML estimation method uses the \link[stats:factanal]{stats::factanal}
   starting values. See also the \link{EFA} documentation.

   The \code{HULL} function can also be called together with other factor
   retention criteria in the \code{\link{N_FACTORS}} function.
}
\examples{
\donttest{
# using PAF (this will throw a warning if gof is not specified manually
# and CAF will be used automatically)
HULL(test_models$baseline$cormat, N = 500, gof = "CAF")

# using ML with all available fit indices (CAF, CFI, and RMSEA)
HULL(test_models$baseline$cormat, N = 500, method = "ML")

# using ULS with only RMSEA
HULL(test_models$baseline$cormat, N = 500, method = "ULS", gof = "RMSEA")
}

\dontrun{
# using parallel processing (Note: plans can be adapted, see the future
# package for details)
future::plan(future::multisession)
HULL(test_models$baseline$cormat, N = 500, gof = "CAF")
}
}
\seealso{
Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
\code{\link{KGC}}, \code{\link{PARALLEL}}, \code{\link{SMT}}

\code{\link{N_FACTORS}} as a wrapper function for this and all the
above-mentioned factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RiskDimensions_doc.R
\docType{data}
\name{RiskDimensions}
\alias{RiskDimensions}
\title{RiskDimensions}
\format{
An object of class \code{list} of length 2.
}
\source{
Fischhoff, B, Slovic, P, Lichtenstein, S, Read, S, and Combs, B. (1978). How safe is safe enough? A psychometric study of attitudes towards technological risks and benefits. Policy Sciences, 9, 127-152. doi: 10.1007/BF00143739
}
\usage{
RiskDimensions
}
\description{
A list containing the bivariate correlations (cormat)
of the 9 dimensions on which participants in Fischhoff et al. (1978) rated
different activities and technologies as well as the sample size (N). This was
then analyzed together with ratings of the risks and benefits of these
activities and technologies.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.KMO.R
\name{print.KMO}
\alias{print.KMO}
\title{Print KMO object}
\usage{
\method{print}{KMO}(x, ...)
}
\arguments{
\item{x}{list of class KMO (output from the \link{KMO} function)}

\item{...}{additional arguments passed to print}
}
\description{
Print KMO object
}
\examples{
KMO_base <- KMO(test_models$baseline$cormat)
KMO_base

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.SCREE.R
\name{print.SCREE}
\alias{print.SCREE}
\title{Print function for SCREE objects}
\usage{
\method{print}{SCREE}(x, plot = TRUE, ...)
}
\arguments{
\item{x}{a list of class SCREE Output from \link{SCREE} function.}

\item{plot}{logical. Whether to plot the results.}

\item{...}{Further arguments for print.}
}
\description{
Print function for SCREE objects
}
\examples{
SCREE_base <- SCREE(test_models$baseline$cormat)
SCREE_base

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.HULL.R
\name{print.HULL}
\alias{print.HULL}
\title{Print function for HULL objects}
\usage{
\method{print}{HULL}(x, plot = TRUE, ...)
}
\arguments{
\item{x}{a list of class HULL. Output from the \code{\link{HULL}} function.}

\item{plot}{logical. Whether to plot the results.}

\item{...}{Further arguments for print.}
}
\description{
Print function for HULL objects
}
\examples{
\donttest{
HULL(test_models$baseline$cormat, N = 500, method = "ML")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.KGC.R
\name{print.KGC}
\alias{print.KGC}
\title{Print function for KGC objects}
\usage{
\method{print}{KGC}(x, plot = TRUE, ...)
}
\arguments{
\item{x}{a list of class KGC. Output from \link{KGC} function.}

\item{plot}{logical. Whether to plot the results.}

\item{...}{Further arguments for print.}
}
\description{
Print function for KGC objects
}
\examples{
KGC_base <- KGC(test_models$baseline$cormat)
KGC_base

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UPPS_raw_doc.R
\docType{data}
\name{UPPS_raw}
\alias{UPPS_raw}
\title{UPPS_raw}
\format{
An object of class \code{data.frame} with 645 rows and 45 columns.
}
\source{
Whiteside, S. P., Lynam, D. R., Miller, J. D., & Reynolds, S. K. (2005).
 Validation of the UPPS impulsive behaviour scale: A four-factor model of
 impulsivity. European Journal of Personality, 19 (7), 559–574.

Steiner, M., & Frey, R. (2020). Representative design in psychological assessment: A case study using the Balloon Analogue Risk Task (BART). PsyArXiv Preprint. doi:10.31234/osf.io/dg4ks
}
\usage{
UPPS_raw
}
\description{
A dataframe containing responses to the UPPS personality scale (Whiteside &
Lynam, 2005) of 645 participants of Study 2 of Steiner and Frey (2020). Each
column are the ratings to one of 45 items to assess urgency, premeditation,
perseverance, and sensation seeking. The original data can be accessed via
\url{https://osf.io/kxp8t/}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.factor_corres}
\alias{.factor_corres}
\title{Compute number of non-matching indicator-to-factor correspondences}
\usage{
.factor_corres(x, y, thresh = 0.3)
}
\arguments{
\item{x}{numeric matrix. A matrix of pattern coefficients.}

\item{y}{numeric matrix. A second matrix of coefficients.}

\item{thresh}{numeric. The threshold to classify a pattern coefficient as substantial.}
}
\description{
Compute number of non-matching indicator-to-factor correspondences
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.SL.R
\name{print.SL}
\alias{print.SL}
\title{Print SL object}
\usage{
\method{print}{SL}(x, ...)
}
\arguments{
\item{x}{list. An object of class SL to be printed}

\item{...}{Further arguments for print.}
}
\description{
Print Method showing a summarized output of the \link{SL} function.
}
\examples{
EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
SL(EFA_mod, type = "EFAtools", method = "PAF")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FACTOR_SCORES.R
\name{FACTOR_SCORES}
\alias{FACTOR_SCORES}
\title{Estimate factor scores for an EFA model}
\usage{
FACTOR_SCORES(
  x,
  f,
  Phi = NULL,
  method = c("Thurstone", "tenBerge", "Anderson", "Bartlett", "Harman", "components"),
  impute = c("none", "means", "median")
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data (needed to get
factor scores) or matrix with correlations.}

\item{f}{object of class \code{\link{EFA}} or matrix.}

\item{Phi}{matrix. A matrix of factor intercorrelations. Only needs to be
specified if a factor loadings matrix is entered directly into \code{f}.
Default is \code{NULL}, in which case all intercorrelations are assumed to be zero.}

\item{method}{character. The method used to calculate factor scores. One of
"Thurstone" (regression-based; default), "tenBerge", "Anderson", "Bartlett",
"Harman", or "components".
See \code{\link[psych:factor.scores]{psych::factor.scores}} for details.}

\item{impute}{character. Whether and how missing values in \code{x} should
be imputed. One of "none" (default, only complete cases are scored), "median",
or "mean".}
}
\value{
A list of class FACTOR_SCORES containing the following:

\item{scores}{The factor scores (only if raw data are provided.)}
\item{weights}{The factor weights.}
\item{r.scores}{The correlations of the factor score estimates.}
\item{missing}{A vector of the number of missing observations per subject
(only if raw data are provided.}
\item{R2}{Multiple R2 of the scores with the factors.}
\item{settings}{A list of the settings used.}
}
\description{
This is a wrapper function for
\code{\link[psych:factor.scores]{psych::factor.scores}} to be used directly
with an output from \code{\link{EFA}} or by manually specifying the factor
loadings and intercorrelations. Calculates factor scores according to the
specified methods if raw data are provided, and only factor weights if a
correlation matrix is provided.
}
\examples{
# Example with raw data with method "Bartlett" and no imputation
EFA_raw <- EFA(DOSPERT_raw, n_factors = 10, type = "EFAtools", method = "PAF",
               rotation = "oblimin")
fac_scores_raw <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "Bartlett",
                                impute = "none")

# Example with a correlation matrix (does not return factor scores)
EFA_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
               type = "EFAtools", method = "PAF", rotation = "oblimin")
fac_scores_cor <- FACTOR_SCORES(test_models$baseline$cormat, f = EFA_cor)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SL.R
\name{SL}
\alias{SL}
\title{Schmid-Leiman Transformation}
\source{
Schmid, J. & Leiman, J. M. (1957). The development of hierarchical
factor solutions. Psychometrika, 22(1), 53–61. doi:10.1007/BF02289209

Wolff, H.-G., & Preising, K. (2005). Exploring item and higher order
factor structure with the Schmid-Leiman solution: Syntax codes for SPSS and
SAS. Behavior Research Methods, 37 , 48–58. doi:10.3758/BF03206397
}
\usage{
SL(
  x,
  Phi = NULL,
  type = c("EFAtools", "psych", "SPSS", "none"),
  method = c("PAF", "ML", "ULS"),
  g_name = "g",
  ...
)
}
\arguments{
\item{x}{object of class \code{\link{EFA}}, class \code{\link[psych:fa]{psych::fa}},
class \code{\link[lavaan]{lavaan}} or matrix. If class \code{\link{EFA}} or
class \code{\link[psych:fa]{psych::fa}}, pattern coefficients and factor
intercorrelations are taken from this object. If class \code{\link[lavaan]{lavaan}},
it must be a second-order CFA solution. In this case first-order and second-order
 factor loadings are taken from this object and the \code{g_name} argument has
 to be specified.
x can also be a pattern matrix from an oblique factor solution (see \code{Phi})
or a matrix of first-order factor loadings from a higher-order confirmatory factor
analysis (see \code{L2}).}

\item{Phi}{matrix. A matrix of factor intercorrelations from an oblique factor
solution. Only needs to be specified if a pattern matrix is entered directly
into \code{x}.}

\item{type}{character. One of "EFAtools" (default), "psych", "SPSS", or "none".
This is used to control the procedure of the second-order factor analysis. See
\code{\link{EFA}} for details.}

\item{method}{character. One of "PAF", "ML", or "ULS" to use
principal axis factoring, maximum likelihood, or unweighted least squares
(also called minres), respectively, used in \code{\link{EFA}} to find the second-order
loadings.}

\item{g_name}{character. The name of the general factor. This needs only be
specified if \code{x} is a \code{lavaan} second-order solution. Default is "g".}

\item{...}{Arguments to be passed to \code{\link{EFA}}.}
}
\value{
A list of class SL containing the following
\item{orig_R}{Original correlation matrix.}
\item{sl}{A matrix with general factor loadings, group factor loadings, communalities,
and uniquenesses.}
\item{L2}{Second-order factor loadings.}
\item{vars_accounted}{A matrix of explained variances and sums of squared loadings.}
\item{iter}{The number of iterations needed for convergence in EFA.}
\item{settings}{list. The settings (arguments) used in EFA to get the
second-order loadings.}
}
\description{
This function implements the Schmid-Leiman (SL) transformation
(Schmid & Leiman, 1957). It takes the pattern coefficients and factor
intercorrelations from an oblique factor solution as
input and can reproduce the results from \code{\link[psych:schmid]{psych::schmid}}
and from the SPSS implementation from Wolff & Preising (2005). Other arguments
from \code{\link{EFA}} can be used to control the procedure to find the
second-order loadings more flexibly. The function can also be used on a
second-order confirmatory factor analysis (CFA) solution from lavaan.
}
\details{
The SL transformation (also called SL orthogonalization) is a procedure with
which an oblique factor solution is transformed into a hierarchical,
orthogonalized solution. As a first step, the factor intercorrelations are
again factor analyzed to find second-order factor loadings. If there is only
one higher-order factor, this step of the procedure stops there, resulting in
a second-order factor structure. The first-order factor and the second-order
factor are then orthogonalized, resulting in an orthogonalized factor solution
with proportionality constraints. The procedure thus makes a suggested
hierarchical data structure based on factor intercorrelations explicit. One
major advantage of SL transformation is that it enables variance
partitioning between higher-order and first-order factors, including the
calculation of McDonald's omegas (see \code{\link{OMEGA}}).
}
\examples{
## Use with an output from the EFAtools::EFA function, both with type EFAtools
EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
SL_EFAtools <- SL(EFA_mod, type = "EFAtools", method = "PAF")

\donttest{
## Use with an output from the psych::fa function with type psych in SL
fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                    fm = "pa", rotate = "Promax")
SL_psych <- SL(fa_mod, type = "psych", method = "PAF")
}

## Use more flexibly by entering a pattern matrix and phi directly (useful if
## a factor solution found with another program should be subjected to SL
## transformation)

## For demonstration, take pattern matrix and phi from an EFA output
## This gives the same solution as the first example
EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
SL_flex <- SL(EFA_mod$rot_loadings, Phi = EFA_mod$Phi, type = "EFAtools",
              method = "PAF")

\donttest{
## Use with a lavaan second-order CFA output

# Create and fit model in lavaan (assume all variables have SDs of 1)
mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ F1 + F2 + F3'
fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                   sample.nobs = 500, estimator = "ml")

SL_lav <- SL(fit, g_name = "g")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WJIV_ages_3_5_doc.R
\docType{data}
\name{WJIV_ages_3_5}
\alias{WJIV_ages_3_5}
\title{Woodcock Johnson IV: ages 3 to 5}
\format{
A list of 2 with elements "cormat" (29 x 29 matrix of bivariate correlations)
and "N" (scalar). The correlation matrix contains the following variables:
\describe{
  \item{ORLVOC}{(numeric) - Oral Vocabulary.}
  \item{VRBATN}{(numeric) - Verbal Attention.}
  \item{LETPAT}{(numeric) - Phonological Processing.}
  \item{STYREC}{(numeric) - Story Recall.}
  \item{VISUAL}{(numeric) - Visualization.}
  \item{GENINF}{(numeric) - General Information.}
  \item{CONFRM}{(numeric) - Concept Formation.}
  \item{NUMREV}{(numeric) - Numbers Reversed.}
  \item{NUMPAT}{(numeric) - Number-Pattern Matching.}
  \item{NWDREP}{(numeric) - Nonword Repetition.}
  \item{VAL}{(numeric) - Visual-Auditory Learning.}
  \item{PICREC}{(numeric) - Picture Recognition.}
  \item{MEMWRD}{(numeric) - Memory for Words.}
  \item{PICVOC}{(numeric) - Picture Vocabulary.}
  \item{ORLCMP}{(numeric) - Oral Comprehension.}
  \item{SEGMNT}{(numeric) - Segmentation.}
  \item{RPCNAM}{(numeric) - Rapid Picture Naming.}
  \item{SENREP}{(numeric) - Sentence Repetition.}
  \item{UNDDIR}{(numeric) - Understanding Directions.}
  \item{SNDBLN}{(numeric) - Sound Blending.}
  \item{RETFLU}{(numeric) - Retrieval Fluency.}
  \item{SNDAWR}{(numeric) - Sound Awareness.}
  \item{LWIDNT}{(numeric) - Letter-Word Identification.}
  \item{APPROB}{(numeric) - Applied Problems.}
  \item{SPELL}{(numeric) - Spelling.}
  \item{PSGCMP}{(numeric) - Passage Comprehension.}
  \item{SCI}{(numeric) - Science.}
  \item{SOC}{(numeric) - Social Studies.}
  \item{HUM}{(numeric) - Humanities.}
 }
}
\source{
McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.
}
\usage{
WJIV_ages_3_5
}
\description{
A list containing the bivariate correlations (N = 435) of the 29
cognitive and achievement subtests from the WJ IV for 3- to 5-year-olds from
the standardization sample obtained from
the WJ IV technical Manual (McGrew, LaForte, & Schrank, 2014).
Tables are reproduced with permission from the publisher.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KGC.R
\name{KGC}
\alias{KGC}
\title{Kaiser-Guttman Criterion}
\source{
Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of extraction
methods under realistic conditions. Psychological Methods, 24(4), 468–491.
https://doi.org/10.1037/met0000200

Guttman, L. (1954). Some necessary conditions for common-factor analysis.
Psychometrika, 19, 149 –161. http://dx.doi.org/10.1007/BF02289162

Kaiser, H. F. (1960). The application of electronic computers to factor
analysis. Educational and Psychological Measurement, 20, 141–151.
http://dx.doi.org/10.1177/001316446002000116

Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain. Psychological Bulletin, 99,
432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
}
\usage{
KGC(
  x,
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall"),
  n_factors = 1,
  ...
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations.}

\item{eigen_type}{character. On what the eigenvalues should be found. Can be
either "PCA", "SMC", or "EFA", or some combination of them. If using "PCA",
the diagonal values of the correlation matrices are left to be 1. If using
"SMC", the diagonal of the
correlation matrices is replaced by the squared multiple correlations (SMCs)
of the indicators. If using "EFA", eigenvalues are found on the correlation
matrices with the final communalities of an exploratory factor analysis
solution (default is principal axis factoring extracting 1 factor) as
diagonal.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
data is given as input. Default is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}

\item{n_factors}{numeric. Number of factors to extract if "EFA" is included in
\code{eigen_type}. Default is 1.}

\item{...}{Additional arguments passed to \code{\link{EFA}}. For example,
to change the extraction method (PAF is default).}
}
\value{
A list of class KGC containing

\item{eigen_PCA}{ A vector containing the eigenvalues found with PCA.}
\item{eigen_SMC}{ A vector containing the eigenvalues found with SMCs.}
\item{eigen_EFA}{ A vector containing the eigenvalues found with EFA.}
\item{n_fac_PCA}{ The number of factors to retain according to the Kaiser-
Guttmann criterion with PCA eigenvalues type.}
\item{n_fac_SMC}{ The number of factors to retain according to the Kaiser-
Guttmann criterion with SMC eigenvalues type.}
\item{n_fac_EFA}{ The number of factors to retain according to the Kaiser-
Guttmann criterion with EFA eigenvalues type.}
\item{settings}{A list of the settings used.}
}
\description{
Probably the most popular factor retention criterion. Kaiser and Guttman suggested
to retain as many factors as there are sample eigenvalues greater than 1.
This is why the criterion is also known as eigenvalues-greater-than-one rule.
}
\details{
Originally, the Kaiser-Guttman criterion was intended for the use
with prinicpal components, hence with eigenvalues derived from the original
correlation matrix. This can be done here by setting \code{eigen_type} to
"PCA". However, it is well-known that this criterion is often inaccurate and
that it tends to overestimate the number of factors, especially for unidimensional
or orthogonal factor structures (e.g., Zwick & Velicer, 1986).

The criterion's inaccuracy in these cases is somewhat addressed if it is
applied on the correlation matrix with communalities in the diagonal, either
initial communalities estimated from SMCs (done setting \code{eigen_type} to
"SMC") or final communality estimates from an EFA (done setting \code{eigen_type}
to "EFA"; see Auerswald & Moshagen, 2019). However, although this variant
of the KGC is more accurate in some cases compared to the traditional KGC, it
is at the same time less accurate than the PCA-variant in other cases, and it
is still often less accurate than other factor retention methods, for
example parallel analysis (\code{\link{PARALLEL}}), the Hull method
\code{\link{HULL}}, or sequential \eqn{chi^2} model tests (\code{\link{SMT}};
see Auerswald & Moshagen, 2019).

The \code{KGC} function can also be called together with other factor
retention criteria in the \code{\link{N_FACTORS}} function.
}
\examples{
KGC(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
}
\seealso{
Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
\code{\link{HULL}}, \code{\link{PARALLEL}}, \code{\link{SMT}}

\code{\link{N_FACTORS}} as a wrapper function for this and all the
above-mentioned factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_models_doc.R
\docType{data}
\name{test_models}
\alias{test_models}
\title{Four test models used in Grieder and Steiner (2020)}
\format{
A list of 4 lists "baseline", "case_1a", "case_6b", and"case_11b", each with the following elements.
\describe{
  \item{cormat}{(matrix) - The correlation matrix of the simulated data.}
  \item{n_factors}{(numeric) - The true number of factors.}
  \item{N}{(numeric) - The sample size of the generated data.}
 }
}
\source{
Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
 in R and SPSS. Manuscript in Preparation.
}
\usage{
test_models
}
\description{
Correlation matrices created from simulated data from four of the
\code{\link{population_models}} cases, each with strong factor intercorrelations.
These are used in Grieder & Steiner (2020) to compare the psych and SPSS
implementations in this package with the actual implementations of the programs.
For details on the cases, see \code{\link{population_models}}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{.compute_vars}
\alias{.compute_vars}
\title{Compute explained variances from loadings}
\usage{
.compute_vars(L_unrot, L_rot, Phi = NULL)
}
\arguments{
\item{L_unrot}{matrix. Unrotated factor loadings.}

\item{L_rot}{matrix. Rotated factor loadings.}

\item{Phi}{matrix. Factor intercorrelations. Provide only if oblique rotation
is used.}
}
\value{
A matrix with sum of squared loadings, proportion explained variance
 from total variance per factor, same as previous but cumulative, Proportion
 of explained variance from total explained variance, and same as previous but
 cumulative.
}
\description{
From unrotated loadings compute the communalities and uniquenesses for total
variance. Compute explained variances per factor from rotated loadings (and
factor intercorrelations Phi if oblique rotation was used).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DOSPERT_raw_doc.R
\docType{data}
\name{DOSPERT_raw}
\alias{DOSPERT_raw}
\title{DOSPERT_raw}
\format{
An object of class \code{data.frame} with 3123 rows and 30 columns.
}
\source{
Blais, A.-R., & Weber, E. U. (2002). A domain-specific risk-taking (DOSPERT) scale for adult populations. Judgment and Decision Making, 15(4), 263–290. doi: 10.1002/bdm.414

Frey, R., Duncan, S. M., & Weber, E. U. (2020). Towards a typology of risk preference: Four risk profiles describe two thirds of individuals in a large sample of the U.S. population. PsyArXiv Preprint. doi:10.31234/osf.io/yjwr9
}
\usage{
DOSPERT_raw
}
\description{
A data.frame containing responses to the risk subscale of the Domain Specific
Risk Taking Scale (DOSPERT; Weber, Blais, & Betz, 2002) based on the publicly
available dataset (at \url{https://osf.io/pjt57/}) by Frey, Duncan, and Weber (2020).
The items measure risk-taking propensity on six different domains: social,
recreational, gambling, health/ safety, investment, and ethical.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPSS_23_doc.R
\docType{data}
\name{SPSS_23}
\alias{SPSS_23}
\title{Various outputs from SPSS (version 23) FACTOR}
\format{
A list of 9 containing EFA results for each of the data sets mentioned above. Each of these nine entries is a list of 4 or 8 (see details), of the following structure:
\describe{
  \item{paf_comm}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details, see Grieder and Grob (2019).}
  \item{paf_load}{(matrix) - F1 to FN = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
  \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
  \item{var_load}{(matrix) - F1 to FN = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
  \item{pro_load}{(matrix) - F1 to FN = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
  \item{pro_phi}{(matrix) - F1 to FN = intercorrelations of the promax rotated loadings.}
  \item{sl}{(matrix) - g = General / second order factor of the Schmid-Leiman solution. F1 to FN  = First order factors of the Schmid-Leiman solution. h2 = Communalities of the Schmid-Leiman solution. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
  \item{L2}{(matrix) - Second order loadings used for the Schmid-Leiman transformation. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
 }
}
\source{
Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
 in R and SPSS. Manuscript in Preparation.

Wolff, H.G., & Preising, K. (2005). Exploring item and higher order factor structure with the Schmid-Leiman solution: Syntax codes for SPSS and SAS. Behavior Research Methods, 37, 48–58. doi: 10.3758/BF03206397

Grieder, S., & Grob, A. (2019). Exploratory factor analyses of the intelligence and development scales–2: Implications for theory and practice. Assessment. Advance online publication. doi:10.1177/10731911198450

Grob, A., & Hagmann-von Arx, P. (2018). Intelligence and Development Scales--2 (IDS-2). Intelligenz- und Entwicklungsskalen für Kinder und Jugendliche.
[Intelligence and Development Scales for Children and Adolescents.]. Bern, Switzerland: Hogrefe.

Frey, R., Pedroni, A., Mata, R., Rieskamp, J., & Hertwig, R. (2017). Risk preference shares the psychometric structure of major psychological traits. Science Advances, 3, e1701381.

McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.

Costa, P. T., & McCrae, R. R. (1992). NEO PI-R professional manual. Odessa, FL: Psychological Assessment Resources, Inc.
}
\usage{
SPSS_23
}
\description{
Various outputs from SPSS (version 23) FACTOR for the IDS-2 (Grob & Hagmann-von Arx, 2018), the WJIV (3 to 5 and 20 to 39 years; McGrew, LaForte, & Schrank, 2014), the DOSPERT (Frey et al., 2017; Weber,
Blais, & Betz, 2002), the NEO-PI-R (Costa, & McCrae, 1992), and four simulated datasets (baseline, case_1a, case_6b, and case_11b, see \link{test_models} and \link{population_models}) used in Grieder and Steiner (2020).
}
\details{
The IDS-2, the two WJIV, the DOSPERT, and the NEO-PI-R contain all the above entries, while the four simulated datasets contain only paf_load, var_load, pro_load, and pro_phi.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.SLLOADINGS.R
\name{print.SLLOADINGS}
\alias{print.SLLOADINGS}
\title{Print SLLOADINGS object}
\usage{
\method{print}{SLLOADINGS}(x, cutoff = 0.2, digits = 3, ...)
}
\arguments{
\item{x}{class SLLOADINGS matrix.}

\item{cutoff}{numeric. The number above which to print loadings in bold
(default is .2).}

\item{digits}{numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
to round the loadings to (default is 3).}

\item{...}{additional arguments passed to print}
}
\description{
Print SLLOADINGS object
}
\examples{
EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
SL(EFA_mod, type = "EFAtools", method = "PAF")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.SCREE.R
\name{plot.SCREE}
\alias{plot.SCREE}
\title{Plot SCREE object}
\usage{
\method{plot}{SCREE}(x, ...)
}
\arguments{
\item{x}{a list of class SCREE An output from the \link{SCREE} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \link{SCREE} function
}
\examples{
SCREE_base <- SCREE(test_models$baseline$cormat)
plot(SCREE_base)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.CD.R
\name{print.CD}
\alias{print.CD}
\title{Print function for CD objects}
\usage{
\method{print}{CD}(x, plot = TRUE, ...)
}
\arguments{
\item{x}{a list of class CD. Output from \link{CD} function.}

\item{plot}{logical. Whether to plot the results.}

\item{...}{Further arguments for print.}
}
\description{
Print function for CD objects
}
\examples{
\donttest{
# determine n factors of the GRiPS
CD(GRiPS_raw)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/N_FACTORS.R
\name{N_FACTORS}
\alias{N_FACTORS}
\title{Various Factor Retention Criteria}
\usage{
N_FACTORS(
  x,
  criteria = c("CD", "EKC", "HULL", "KGC", "PARALLEL", "SCREE", "SMT"),
  suitability = TRUE,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall"),
  n_factors_max = NA,
  N_pop = 10000,
  N_samples = 500,
  alpha = 0.3,
  max_iter_CD = 50,
  n_fac_theor = NA,
  method = c("PAF", "ULS", "ML"),
  gof = c("CAF", "CFI", "RMSEA"),
  eigen_type_HULL = c("SMC", "PCA", "EFA"),
  eigen_type_other = c("PCA", "SMC", "EFA"),
  n_factors = 1,
  n_datasets = 1000,
  percent = 95,
  decision_rule = c("means", "percentile", "crawford"),
  show_progress = TRUE,
  ...
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations. If \code{"CD"} is included as a criterion, x must be raw
 data.}

\item{criteria}{character. A vector with the factor retention methods to
perform. Possible inputs are: \code{"CD"}, \code{"EKC"}, \code{"HULL"},
\code{"KGC"}, \code{"PARALLEL"}, \code{"SCREE"}, and \code{"SMT"}
(see details). By default, all factor retention methods are performed.}

\item{suitability}{logical. Whether the data should be checked for suitability
for factor analysis using the Bartlett's test of sphericity and the
Kaiser-Guttmann criterion (see details). Default is \code{TRUE}.}

\item{N}{numeric. The number of observations. Only needed if x is a
correlation matrix.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
data is given as input. Default is \code{"pairwise.complete.obs"}.}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}
Default is  \code{"pearson"}.}

\item{n_factors_max}{numeric. Passed to \code{\link{CD}}.The maximum number
of factors to test against.
Larger numbers will increase the duration the procedure takes, but test more
possible solutions. Maximum possible is number of variables / 2. Default is
NA. If not specified, number of variables / 2 is used.}

\item{N_pop}{numeric. Passed to \code{\link{CD}}. Size of finite populations
of comparison data. Default is 10000.}

\item{N_samples}{numeric. Passed to \code{\link{CD}}. Number of samples drawn
from each population. Default is 500.}

\item{alpha}{numeric. Passed to \code{\link{CD}}. The alpha level used to test
the significance of the improvement added by an additional factor.
Default is .30.}

\item{max_iter_CD}{numeric. Passed to \code{\link{CD}}. The maximum number of
iterations to perform after which the iterative PAF procedure is halted.
 Default is 50.}

\item{n_fac_theor}{numeric. Passed to \code{\link{HULL}}. Theoretical number
of factors to retain. The maximum of this number and the number of factors
suggested by \link{PARALLEL} plus one will be used in the Hull method.}

\item{method}{character. Passed to \code{\link{EFA}} in \code{\link{HULL}},
\code{\link{KGC}}, \code{\link{SCREE}}, and \code{\link{PARALLEL}}. The
estimation method to use. One of  \code{"PAF"}, \code{"ULS"}, or  \code{"ML"},
for principal axis factoring, unweighted least squares, and maximum
likelihood, respectively.}

\item{gof}{character. Passed to \code{\link{HULL}}. The goodness of fit index
to use. Either \code{"CAF"}, \code{"CFI"}, or \code{"RMSEA"}, or any
combination of them. If \code{method = "PAF"} is used, only
the CAF can be used as goodness of fit index. For details on the CAF, see
Lorenzo-Seva, Timmerman, and Kiers (2011).}

\item{eigen_type_HULL}{character. Passed to  \code{\link{PARALLEL}} in
\code{\link{HULL}}. On what the
eigenvalues should be found in the parallel analysis. Can be one of
\code{"SMC"}, \code{"PCA"}, or \code{"EFA"}. If using  \code{"SMC"} (default),
the diagonal of the correlation matrices is
replaced by the squared multiple correlations (SMCs) of the indicators. If
using  \code{"PCA"}, the diagonal values of the correlation
matrices are left to be 1. If using  \code{"EFA"}, eigenvalues are found on the
correlation  matrices with the final communalities of an EFA solution as
diagonal.}

\item{eigen_type_other}{character. Passed to \code{\link{KGC}},
\code{\link{SCREE}}, and \code{\link{PARALLEL}}. The same as eigen_type_HULL,
but multiple inputs
are possible here. Default is to use all inputs, that is, \code{c("PCA",
"SMC", "EFA"})}

\item{n_factors}{numeric. Passed to \code{\link{PARALLEL}} (also within
\code{\link{HULL}}), \code{\link{KGC}}, and \code{\link{SCREE}}. Number of
factors to extract if \code{"EFA"} is included in \code{eigen_type_HULL} or
 \code{eigen_type_other}. Default is 1.}

\item{n_datasets}{numeric. Passed to \code{\link{PARALLEL}} (also within
\code{\link{HULL}}). The number of datasets to simulate. Default is 1000.}

\item{percent}{numeric. Passed to \code{\link{PARALLEL}} (also within
\code{\link{HULL}}). A vector of percentiles to take the simulated eigenvalues
 from. Default is 95.}

\item{decision_rule}{character. Passed to \code{\link{PARALLEL}} (also within
\code{\link{HULL}}). Which rule to use to determine the number of
 factors to retain. Default is \code{"means"}, which will use the average
 simulated eigenvalues. \code{"percentile"}, uses the percentiles specified
 in percent. \code{"crawford"} uses the 95th percentile for the first factor
 and the mean afterwards (based on Crawford et al, 2010).}

\item{show_progress}{logical. Whether a progress bar should be shown in the
console. Default is TRUE.}

\item{...}{Further arguments passed to \code{\link{EFA}} in
\code{\link{PARALLEL}} (also within \code{\link{HULL}}) and \code{\link{KGC}}.}
}
\value{
A list of class N_FACTORS containing
\item{outputs}{A list with the outputs from \code{\link{BARTLETT}} and
 \code{\link[EFAtools]{KMO}} and the factor retention criteria.}
\item{n_factors}{A named vector containing the suggested number of factors
from each factor retention criterion.}
\item{settings}{A list of the settings used.}
}
\description{
Among the most important decisions for an exploratory factor analysis (EFA) is
the choice of the number of factors to retain. Several factor retention
criteria have been developed for this. With this function, various factor
 retention criteria can be performed simultaneously. Additionally, the data
 can be checked for their suitability for factor analysis.
}
\details{
By default, the entered data are checked for suitability for factor analysis
using the following methods (see respective documentations for details):
\itemize{
\item{Bartlett's test of sphericity (see \code{\link{BARTLETT}})}
\item{Kaiser-Meyer-Olkin criterion (see \code{\link[EFAtools]{KMO}})}}

The available factor retention criteria are the following (see respective
 documentations for details):
 \itemize{
\item{Comparison data (see \code{\link{CD}})}
\item{Empirical Kaiser criterion (see \code{\link{EKC}})}
\item{Hull method (see \code{\link{HULL}})}
\item{Kaiser-Guttman criterion (see \code{\link{KGC}})}
\item{Parallel analysis (see \code{\link{PARALLEL}})}
\item{Scree plot (see \code{\link{SCREE}})}
\item{Sequential chi-square model tests, RMSEA lower bound, and AIC
(see \code{\link{SMT}})}
}
}
\examples{
\donttest{
# All criteria, with correlation matrix and fit method "ML" (where needed)
# This will throw a warning for CD, as no raw data were specified
nfac_all <- N_FACTORS(test_models$baseline$cormat, N = 500, method = "ML")

# The same as above, but without "CD"
nfac_wo_CD <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
                        "HULL", "KGC", "PARALLEL", "SCREE", "SMT"), N = 500,
                        method = "ML")

# Use PAF instead of ML (this will take a lot longer). For this, gof has
# to be set to "CAF" for the Hull method.
nfac_PAF <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
                      "HULL", "KGC", "PARALLEL", "SCREE", "SMT"), N = 500,
                      gof = "CAF")

# Do KGC and PARALLEL with only "PCA" type of eigenvalues
nfac_PCA <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
                      "HULL", "KGC", "PARALLEL", "SCREE", "SMT"), N = 500,
                      method = "ML", eigen_type_other = "PCA")

# Use raw data, such that CD can also be performed
nfac_raw <- N_FACTORS(GRiPS_raw, method = "ML")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.CD.R
\name{plot.CD}
\alias{plot.CD}
\title{Plot CD object}
\usage{
\method{plot}{CD}(x, ...)
}
\arguments{
\item{x}{a list of class CD. An output from the \link{CD} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \link{CD} function
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PARALLEL.R
\name{PARALLEL}
\alias{PARALLEL}
\title{Parallel analysis}
\source{
Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
Psychological Methods, 22, 450 – 466. http://dx.doi.org/10.1037/ met0000074

Crawford, A. V., Green, S. B., Levy, R., Lo, W. J., Scott, L.,
Svetina, D., & Thompson, M. S. (2010). Evaluation of parallel analysis methods
for determining the number of factors. Educational and Psychological
Measurement, 70(6), 885-901.

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. Psychometrika, 30(2), 179–185. doi: 10.1007/BF02289447
}
\usage{
PARALLEL(
  x = NULL,
  N = NA,
  n_vars = NA,
  n_datasets = 1000,
  percent = 95,
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall"),
  decision_rule = c("means", "percentile", "crawford"),
  n_factors = 1,
  ...
)
}
\arguments{
\item{x}{matrix or data.frame. The real data to compare the simulated eigenvalues
against. Must not contain variables of classes other than numeric. Can be a
correlation matrix or raw data.}

\item{N}{numeric. The number of cases / observations to simulate. Only has to
be specified if \code{x} is either a correlation matrix or \code{NULL}. If
x contains raw data, \code{N} is found from the dimensions of \code{x}.}

\item{n_vars}{numeric. The number of variables / indicators to simulate.
Only has to be specified if \code{x} is left as \code{NULL} as otherwise the
dimensions are taken from \code{x}.}

\item{n_datasets}{numeric. The number of datasets to simulate. Default is 1000.}

\item{percent}{numeric. The percentile to take from the simulated eigenvalues.
Default is 95.}

\item{eigen_type}{character. On what the eigenvalues should be found. Can be
either "SMC", "PCA", or "EFA". If using "SMC", the diagonal of the correlation
matrix is replaced by the squared multiple correlations (SMCs) of the
indicators. If using "PCA", the diagonal values of the correlation matrices
are left to be 1. If using "EFA", eigenvalues are found on the correlation
matrices with the final communalities of an EFA solution as diagonal.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
is given as input. Default is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}
Default is "pearson".}

\item{decision_rule}{character. Which rule to use to determine the number of
factors to retain. Default is \code{"means"}, which will use the average
simulated eigenvalues. \code{"percentile"}, uses the percentiles specified
in percent. \code{"crawford"} uses the 95th percentile for the first factor
and the mean afterwards (based on Crawford et al, 2010).}

\item{n_factors}{numeric. Number of factors to extract if "EFA" is included in
\code{eigen_type}. Default is 1.}

\item{...}{Additional arguments passed to \code{\link{EFA}}. For example,
the extraction method can be changed here (default is "PAF"). PAF is more
robust, but it will take longer compared to the other estimation methods
available ("ML" and "ULS").}
}
\value{
A list of class PARALLEL containing the following objects
\item{eigenvalues_PCA}{A matrix containing the eigenvalues of the real and the simulated data found with eigen_type = "PCA"}
\item{eigenvalues_SMC}{A matrix containing the eigenvalues of the real and the simulated data found with eigen_type = "SMC"}
\item{eigenvalues_EFA}{A matrix containing the eigenvalues of the real and the simulated data found with eigen_type = "EFA"}
\item{n_fac_PCA}{The number of factors to retain according to the parallel procedure with eigen_type = "PCA".}
\item{n_fac_SMC}{The number of factors to retain according to the parallel procedure with eigen_type = "SMC".}
\item{n_fac_EFA}{The number of factors to retain according to the parallel procedure with eigen_type = "EFA".}
\item{settings}{A list of control settings used in the print function.}
}
\description{
Various methods for performing parallel analysis. This function uses
\link[future.apply]{future_lapply} for which a parallel processing plan can
be selected. To do so, call \code{library(future)} and, for example,
 \code{plan(multisession)}; see examples.
}
\details{
Parallel analysis (Horn, 1965) compares the eigenvalues obtained from
the sample
 correlation matrix against those of null model correlation matrices (i.e.,
 with uncorrelated variables) of the same sample size. This way, it accounts
 for the variation in eigenvalues introduced by sampling error and thus
 eliminates the main problem inherent in the Kaiser-Guttman criterion
 (\code{\link{KGC}}).

 Three different ways of finding the eigenvalues under the factor model are
 implemented, namely "SMC", "PCA", and "EFA". PCA leaves the diagonal elements
 of the correlation matrix as they are and is thus equivalent to what is done
 in PCA. SMC uses squared multiple correlations as communality estimates with
 which the diagonal of the correlation matrix is replaced. Finally, EFA performs
 an \code{\link{EFA}} with one factor (can be adapted to more factors) to estimate
 the communalities and based on the correlation matrix with these as diagonal
 elements, finds the eigenvalues.

 Parallel analysis is often argued to be one of the most accurate factor
 retention criteria. However, for highly correlated
 factor structures it has been shown to underestimate the correct number of
 factors. The reason for this is that a null model (uncorrelated variables)
 is used as reference. However, when factors are highly correlated, the first
 eigenvalue will be much larger compared to the following ones, as
 later eigenvalues are conditional on the earlier ones in the sequence and thus
 the shared variance is already accounted in the first eigenvalue (e.g.,
 Braeken & van Assen, 2017).

 The \code{PARALLEL} function can also be called together with other factor
 retention criteria in the \code{\link{N_FACTORS}} function.
}
\examples{
\donttest{
# example without real data
pa_unreal <- PARALLEL(N = 500, n_vars = 10)

# example with correlation matrix with all eigen_types and PAF estimation
pa_paf <- PARALLEL(test_models$case_11b$cormat, N = 500)

# example with correlation matrix with all eigen_types and ML estimation
# this will be faster than the above with PAF)
pa_ml <- PARALLEL(test_models$case_11b$cormat, N = 500, method = "ML")
}

\dontrun{
# for parallel computation
future::plan(future::multisession)
pa_faster <- PARALLEL(test_models$case_11b$cormat, N = 500)
}
}
\seealso{
Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
\code{\link{HULL}}, \code{\link{KGC}}, \code{\link{SMT}}

\code{\link{N_FACTORS}} as a wrapper function for this and all the
above-mentioned factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMT.R
\name{SMT}
\alias{SMT}
\title{Sequential Chi Square Model Tests, RMSEA lower bound, and AIC}
\source{
Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of extraction
methods under realistic conditions. Psychological Methods, 24(4), 468–491.
https://doi.org/10.1037/met0000200

Browne, M.W., & Cudeck, R. (1992). Alternative ways of assessing model
fit. Sociological Methods and Research, 21, 230–258.

Preacher, K. J., Zhang G., Kim, C., & Mels, G. (2013). Choosing the
Optimal Number of Factors in Exploratory Factor Analysis: A Model Selection
Perspective, Multivariate Behavioral Research, 48(1), 28-56,
doi:10.108/00273171.2012.710386

Steiger, J. H., & Lind, J. C. (1980, May). Statistically based tests
for the number of common factors. Paper presented at the annual meeting of
the Psychometric Society, Iowa City, IA.
}
\usage{
SMT(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall")
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations.}

\item{N}{numeric. The number of observations. Needs only be specified if a
correlation matrix is used.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
data is given as input. Default is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}
}
\value{
A list of class SMT containing
\item{nfac_chi}{The number of factors to retain according to the significance
of the chi square value.}
\item{nfac_RMSEA}{The number of factors to retain according to the RMSEA lower
bound}
\item{nfac_AIC}{The number of factors to retain according to the AIC}
\item{p_null}{The p-value for the null model (zero factors)}
\item{ps_chi}{The p-values for EFA models with increasing numbers of factors,
starting with 1 factor}
\item{RMSEA_LB_null}{The lower bounds of the 90\% confidence interval for the RMSEA
for the null model (zero factors).}
\item{RMSEA_LBs}{The lower bounds of the 90\% confidence interval for the RMSEA
for EFA models with increasing numbers of factors, starting with 1 factor}
\item{AIC_null}{The AICs for the null model (zero factors)}
\item{AICs}{The AICs for EFA models with increasing numbers of factors,
starting with 1 factor}
}
\description{
Sequential Chi Square Model Tests (SMT) are a factor retention method where
multiple
EFAs with increasing numbers of factors are fitted and the number of factors
for which the Chi Square value first becomes non-significant is taken as the
suggested number of factors.
Preacher, Zhang, Kim, & Mels (2013) suggested a similar approach with the
lower bound of the 90\% confidence interval of the Root Mean Square Error of
Approximation (RMSEA; Browne & Cudeck, 1992; Steiger & Lind, 1980), and with
the Akaike Information Criterion (AIC). For the RMSEA, the
number of factors for which this lower bound first falls below .05 is the
suggested number of factors to retain. For the AIC, it is the number of factors
where the AIC is lowest.
}
\details{
As a first step in the procedure, a maximum number of factors to extract is
determined for which the model is still over-identified (df > 0).

Then, EFAs with increasing numbers of factors from 1 to the maximum number are
fitted with maximum likelihood estimation.

For the SMT, first the significance of the chi
square value for a model with 0 factors is determined. If this value is
not significant, 0 factors are suggested to retain. If it is significant,
a model with 1 factor is estimated and the significance of its chi square value
is determined, and so on, until a non-significant result is obtained. The
suggested number of factors is the number of factors for the model where the
chi square value first becomes non-significant.

Regarding the RMSEA, the suggested number of factors is the number of factors
for the model where the lower bound of the 90\% confidence interval of the
RMSEA first falls below the .05 threshold.

Regarding the AIC, the suggested number of factors is the number of factors
for the model with the lowest AIC.

In comparison with other prominent factor retention criteria, SMT performed
well at determining the number of factors to extract in EFA (Auerswald &
Moshagen, 2019). The RMSEA lower bound also performed well at determining the true
number of factors, while the AIC performed well at determining the
most generalizable model (Preacher, Zhang, Kim, & Mels, 2013).

The \code{SMT} function can also be called together with other factor
retention criteria in the \code{\link{N_FACTORS}} function.
}
\examples{
SMT_base <- SMT(test_models$baseline$cormat, N = 500)
SMT_base

}
\seealso{
Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
\code{\link{HULL}}, \code{\link{KGC}}, \code{\link{PARALLEL}}

\code{\link{N_FACTORS}} as a wrapper function for this and all the
above-mentioned factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{.numformat}
\alias{.numformat}
\title{Format numbers for print method}
\usage{
.numformat(x, digits = 2, print_zero = FALSE)
}
\arguments{
\item{x}{numeric. Number to be formatted.}

\item{digits}{numeric. Number of digits after the comma to keep.}

\item{print_zero}{logical. Whether, if a number is between [-1, 1], the
zero should be omitted or printed (default is FALSE, i.e. omit zeros).}
}
\value{
A formated number
}
\description{
Helper function used in the print method for class LOADINGS and SLLOADINGS.
Strips the 0 in front of the decimal point of a number if number < 1, only
keeps the first \code{digits} number of digits, and adds an empty space in
front of the number if the number is positive. This way all returned strings
(except for those > 1, which are exceptions in LOADINGS) have the same number
of characters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.OMEGA.R
\name{print.OMEGA}
\alias{print.OMEGA}
\title{Print OMEGA object}
\usage{
\method{print}{OMEGA}(x, digits = 3, ...)
}
\arguments{
\item{x}{output of class OMEGA (output from the \link{OMEGA} function)}

\item{digits}{numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
to round to (default is 3).}

\item{...}{additional arguments passed to print}
}
\description{
Print OMEGA object
}
\examples{
efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")

OMEGA(sl_mod, type = "EFAtools",
factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.paf_iter}
\alias{.paf_iter}
\title{Perform the iterative PAF procedure}
\usage{
.paf_iter(h2, criterion, R, n_fac, abs_eig, crit_type, max_iter)
}
\arguments{
\item{h2}{numeric. The initial communality estimates.}

\item{criterion}{double. The convergence criterion to use.}

\item{R}{matrix. The correlation matrix with the initial communality estimates in the diagonal.}

\item{n_fac}{numeric. The number of factors to extract.}

\item{abs_eig}{logical. Whether absolute eigenvalues should be used to compute the loadings.}

\item{crit_type}{numeric. Whether maximum absolute differences (crit_type = 1), or sum of differences (crit_type = 2) should be used}

\item{max_iter}{numeric. The number of iterations after which to end the procedure if no convergence has been reached by then.}
}
\description{
Function called from within PAF so usually no call to this is needed by the user.
Provides a C++ implementation of the PAF procedure
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.BARTLETT.R
\name{print.BARTLETT}
\alias{print.BARTLETT}
\title{Print BARTLETT object}
\usage{
\method{print}{BARTLETT}(x, ...)
}
\arguments{
\item{x}{list of class BARTLETT (output from the \link{BARTLETT} function)}

\item{...}{additional arguments passed to print}
}
\description{
Print BARTLETT object
}
\examples{
BARTLETT(test_models$baseline$cormat, N = 500)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IDS2_R_doc.R
\docType{data}
\name{IDS2_R}
\alias{IDS2_R}
\title{Intelligence subtests from the Intelligence and Development Scales--2}
\format{
A 14 x 14 matrix of bivariate correlations
\describe{
  \item{GS}{(numeric) - Geometric shapes.}
  \item{PL}{(numeric) - Plates.}
  \item{TC}{(numeric) - Two characteristics.}
  \item{CB}{(numeric) - Crossing out boxes.}
  \item{NL}{(numeric) - Numbers / letters.}
  \item{NLM}{(numeric) - Numbers / letter mixed.}
  \item{GF}{(numeric) - Geometric figures.}
  \item{RGF}{(numeric) - Rotated geometric figures.}
  \item{CM}{(numeric) - Completing matrices.}
  \item{EP}{(numeric) - Excluding pictures.}
  \item{CA}{(numeric) - Categories.}
  \item{OP}{(numeric) - Opposites.}
  \item{RS}{(numeric) - Retelling a story.}
  \item{DP}{(numeric) - Describing pictures.}
 }
}
\source{
Grieder, S., & Grob, A. (2019). Exploratory factor analyses of the intelligence and development scales--2: Implications for theory and practice. Assessment. Advance online publication. doi:10.1177/10731911198450

Grob, A., & Hagmann-von Arx, P. (2018). Intelligence and Development Scales--2 (IDS-2). Intelligenz- und Entwicklungsskalen für Kinder und Jugendliche.
[Intelligence and Development Scales for Children and Adolescents.]. Bern, Switzerland: Hogrefe.
}
\usage{
IDS2_R
}
\description{
A matrix containing the bivariate correlations of the 14 intelligence subtests from the Intelligence and Development Scales--2 (IDS-2; Grob & Hagmann-von Arx, 2018), an intelligence and development test battery for children and adolescents aged 5 to 20 years, for the standardization and validation sample (N = 1,991). Details can be found in Grieder & Grob (2019).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.N_FACTORS.R
\name{print.N_FACTORS}
\alias{print.N_FACTORS}
\title{Print function for N_FACTORS objects}
\usage{
\method{print}{N_FACTORS}(x, ...)
}
\arguments{
\item{x}{a list of class N_FACTORS. Output from \link{N_FACTORS} function.}

\item{...}{Further arguments for print.}
}
\description{
Print function for N_FACTORS objects
}
\examples{
\donttest{
# All criteria except "CD", with correlation matrix and fit method "ML"
# (where needed)
N_FACTORS(test_models$baseline$cormat, criteria = c("EKC", "HULL", "KGC",
          "PARALLEL", "SCREE", "SMT"), N = 500, method = "ML")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EFAtools-package.R
\docType{package}
\name{EFAtools-package}
\alias{EFAtools}
\alias{EFAtools-package}
\title{EFAtools: Fast and Flexible Implementations of Exploratory Factor Analysis Tools}
\description{
Provides functions to perform exploratory factor analysis (EFA) procedures and compare their solutions. The goal is to provide state-of-the-art factor retention methods and a high degree of flexibility in the EFA procedures. This way, for example, implementations from R 'psych' and 'SPSS' can be compared. Moreover, functions for Schmid-Leiman transformation and the computation of omegas are provided. To speed up the analyses, some of the iterative procedures, like principal axis factoring (PAF), are implemented in C++.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/mdsteiner/EFAtools}
  \item Report bugs at \url{https://github.com/mdsteiner/EFAtools/issues}
}

}
\author{
\strong{Maintainer}: Markus Steiner \email{markus.d.steiner@gmail.com}

Authors:
\itemize{
  \item Silvia Grieder \email{silvia.grieder@gmail.com}
}

Other contributors:
\itemize{
  \item William Revelle [contributor]
  \item Max Auerswald [contributor]
  \item Morten Moshagen [contributor]
  \item John Ruscio [contributor]
  \item Brendan Roche [contributor]
  \item Urbano Lorenzo-Seva [contributor]
  \item David Navarro-Gonzalez [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DOSPERT_doc.R
\docType{data}
\name{DOSPERT}
\alias{DOSPERT}
\title{DOSPERT}
\format{
An object of class \code{list} of length 2.
}
\source{
Weber, E. U., Blais, A.-R., & Betz, N. E. (2002). A domain specific risk-attitude scale: Measuring risk perceptions and risk behaviors. Journal of Behavioral Decision Making, 15(4), 263–290. doi: 10.1002/bdm.414

Frey, R., Pedroni, A., Mata, R., Rieskamp, J., & Hertwig, R. (2017). Risk preference shares the psychometric structure of major psychological traits. Science Advances, 3, e1701381.

\url{https://osf.io/rce7g}
}
\usage{
DOSPERT
}
\description{
A list containing the the bivariate correlations (cormat) of the 40 items of
the Domain Specific Risk Taking Scale (DOSPERT; Weber, Blais, & Betz, 2002)
and the sample size (N) based on the publicly available dataset at
(\url{https://osf.io/rce7g}) of the Basel-Berlin Risk Study (Frey et al., 2017).
The items measure risk-taking propensity on six different domains: social,
recreational, gambling, health/ safety, investment, and ethical.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OMEGA.R
\name{OMEGA}
\alias{OMEGA}
\title{McDonald's omega}
\source{
McDonald, R. P. (1978). Generalizability in factorable domains: ‘‘Domain
validity and generalizability’’. Educational and Psychological Measurement,
38, 75–79.

McDonald, R. P. (1985). Factor analysis and related methods. Hillsdale,
NJ: Erlbaum.

McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah,
NJ: Erlbaum.

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016a). Applying bifactor
statistical indices in the evaluation of psychological measures. Journal of
Personality Assessment, 98, 223-237.

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016b). Evaluating
bifactor models: Calculating and interpreting statistical indices.
Psychological Methods, 21, 137-150.

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct reliability
within latent variable systems. In R. Cudeck, S. du Toit, & D. Sörbom (Eds.),
Structural equation modeling: Present and future—A Festschrift in honor of Karl
Jöreskog (pp. 195–216). Lincolnwood, IL: Scientific Software International.

Sijtsma, K. (2009). On the use, the misuse, and the very limited usefulness
of Cronbach’s alpha. Psychometrika, 74, 107–120.

Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013).
Multidimensionality and structural coefficient bias in structural equation
modeling: A bifactor perspective. Educational and Psychological Measurement,
73, 5–26.

Bonifay, W. E., Reise, S. P., Scheines, R., & Meijer, R. R. (2015).
When are multidimensional data unidimensional enough for structural equation
modeling?: An evaluation of the DETECT multidimensionality index. Structural
Equation Modeling, 22, 504—516.

Gignac, G. E. (2014). On the Inappropriateness of Using Items to
Calculate Total Scale Score Reliability via Coefficient Alpha for Multidimensional
Scales. European Journal of Psychological Assessment, 30, 130-139.
}
\usage{
OMEGA(
  model = NULL,
  type = c("EFAtools", "psych"),
  g_name = "g",
  group_names = NULL,
  add_ind = TRUE,
  factor_corres = NULL,
  var_names = NULL,
  fac_names = NULL,
  g_load = NULL,
  s_load = NULL,
  u2 = NULL,
  cormat = NULL,
  pattern = NULL,
  Phi = NULL,
  variance = c("correlation", "sums_load")
)
}
\arguments{
\item{model}{class \code{\link{SL}}, class \code{\link{schmid}}, or class
\code{lavaan} object. That is, an output object from \code{\link{SL}} or
\code{\link[psych:schmid]{psych::schmid}}, or a \code{lavaan} fit object with a
single factor, second-order, or bifactor solution. If of class \code{lavaan},
only \code{g_name} needs to be specified additionally. If of class
\code{\link{SL}} or \code{\link{schmid}}, only the arguments \code{factor_corres}
and \code{cormat} need to be specified additionally.}

\item{type}{character. Either \code{"EFAtools"} (default) or \code{"psych"}
(see details)}

\item{g_name}{character. The name of the general factor from the lavaan solution.
This needs only be specified if \code{model} is a \code{lavaan} second-order
or bifactor solution. Default is "g".}

\item{group_names}{character. An optional vector of group names. The length
must correspond to the number of groups for which the \code{lavaan} model
was fitted.}

\item{add_ind}{logical. Whether additional indices (H index, ECV, PUC) should
be calculated or not (see details for these indices). If FALSE, only omegas
are returned. Default is \code{TRUE}.}

\item{factor_corres}{matrix. A logical matrix or a numeric matrix containing
0's and 1's that indicates which variable corresponds to which group factor.
Must have the same dimensions as the matrix of group factor loadings from the
SL solution. Cross-loadings are allowed here. See examples for use.}

\item{var_names}{character. A vector with subtest names in the order
of the rows from the SL solution. This needs only be specified if \code{model}
is left \code{NULL}.}

\item{fac_names}{character. An optional vector of group factor names in the
order of the columns of the SL solution. If left \code{NULL}, names of the
group factors from the entered solution are taken.}

\item{g_load}{numeric. A vector of general factor loadings from an SL solution.
This needs only be specified if \code{model} is left \code{NULL}.}

\item{s_load}{matrix. A matrix of group factor loadings from an SL solution.
This needs only be specified if \code{model} is left \code{NULL}.}

\item{u2}{numeric. A vector of uniquenesses from an SL solution. This needs
only be specified if \code{model} is left \code{NULL}.}

\item{cormat}{matrix. A correlation matrix to be used when
\code{variance = "correlation"}. If left \code{NULL} and an \code{\link{SL}}
output is entered in \code{model}, the correlation matrix is taken from the
output. If left \code{NULL} and a \code{\link[psych:schmid]{psych::schmid}}
output is entered, the correlation matrix will be found based on the pattern
matrix and Phi from the \code{\link[psych:schmid]{psych::schmid}} output
using \code{\link[psych:factor.model]{psych::factor.model}}.
If left \code{NULL} and model is also left \code{NULL}, the correlation matrix
is found based on the pattern matrix and Phi entered. However, if the
correlation matrix is available, \code{cormat} should be specified instead
of \code{Phi} and \code{pattern}.}

\item{pattern}{matrix. Pattern coefficients from an oblique factor solution.
This needs only be specified if \code{model} is left \code{NULL},
\code{variance = "correlation"} and \code{cormat} is also left \code{NULL}.}

\item{Phi}{matrix. Factor intercorrelations from an oblique factor solution.
This needs only be specified if \code{model} is left \code{NULL},
\code{variance = "correlation"} and \code{cormat} is also left \code{NULL}.}

\item{variance}{character. If \code{"correlation"} (default), then total
variances for the whole scale as well as for the subscale composites are
calculated based on the correlation
matrix. If \code{"sums_load"}, then total variances are calculated using the
squared sums of general factor loadings and group factor loadings and
the sum of uniquenesses (see details).}
}
\value{
If found for an SL or \code{lavaan} second-order of bifactor solution
without multiple groups:
A matrix with omegas for the whole scale and for the subscales and (only if
\code{add_ind = TRUE}) with the H index, ECV, and PUC.
\item{tot}{Omega total.}
\item{hier}{Omega hierarchical.}
\item{sub}{Omega subscale.}
\item{H}{H index.}
\item{ECV}{Explained common variance.}
\item{PUC}{Percent of uncontaminated correlations.}

If found for a \code{lavaan} single factor solution without multiple groups:
A (named) vector with omega total and (if \code{add_ind = TRUE}) the H index
for the single factor.

If found for a \code{lavaan} output from a multiple group analysis: A list
containing the output described above for each group.
}
\description{
This function finds omega total, hierarchical, and subscale, as well as additional
model-based indices of interpretive relevance (H index, ECV, PUC)
from a Schmid-Leiman (SL) solution or lavaan single factor, second-order (see below),
or bifactor solution. The SL-based omegas can either be found from a
\code{\link[psych:schmid]{psych::schmid}}, \code{\link{SL}}, or,
in a more flexible way, by leaving
\code{model = NULL} and specifying additional arguments. By setting the
\code{type} argument, results from \code{\link[psych:omega]{psych::omega}}
can be reproduced.
}
\details{
## What this function does

This function calculates McDonald's omegas (McDonald, 1978, 1985, 1999),
the H index (Hancock & Mueller, 2001), the explained common variance (ECV;
Sijtsma, 2009), and the percent of uncontaminated correlations (PUC; Bonifay
et al., 2015; Reise et al., 2013).

All types of omegas (total, hierarchical, and subscale) are calculated for
the general factor as well as for the subscales / group factors (see, e.g.,
Gignac, 2014; Rodriguez et al., 2016a, 2016b). Omegas refer to the correlation
between a factor and a unit-weighted composite score and thus the
true score variance in a unit-weighted composite based on the respective
indicators. Omega total is the total true score variance in a composite.
Omega hierarchical is the true score variance in a composite that is attributable
to the general factor, and omega subscale is the true score variance in a
composite attributable to all subscales / group factors (for the whole scale)
or to the specific subscale / group factor (for subscale composites).

The H index (also construct reliability or replicability index) is the
correlation between an optimally-weighted composite score
and a factor (Hancock & Mueller, 2001; Rodriguez et al., 2016a, 2016b). It, too,
can be calculated for the whole scale / general factor as well as for the
subscales / grouup factors. Low values indicate that a latent variable is not well
defined by its indicators.

The ECV (Sijtsma, 2009, Rodriguez et al., 2016a, 2016b) is the ratio of the
variance explained by the general factor and the variance explained by the
general factor and the group factors.

The PUC (Bonifay et al., 2015; Reise et al., 2013, Rodriguez et al., 2016a,
2016b) refers to the proportion
of correlations in the underlying correlation matrix that is not contaminated
by variance of both the general factor and the group factors (i.e., correlations
between indicators from different group factors, which reflect only general
factor variance). The higher the PUC, the more similar a general factor from
a multidimensional model will be to the single factor from a unidimensional
model.

## How to use this function

If \code{model} is a \code{lavaan} second-order or bifactor solution,
only the name of the general factor from the lavaan model needs to be specified
additionally with the \code{g_name} argument. It is then determined whether this
general factor is a second-order factor (second-order model with one second-order
factor assumed) or a breadth factor (bifactor model assumed). Please note that
this function only works for second-order models if they contain no more than
one second-order factor. In case of a second-order solution, a
Schmid-Leiman transformation is performed on the first- and second-order loadings
and omega coefficents are obtained from the transformed (orthogonalized) solution
(see \code{\link{SL}} for more information on Schmid-Leiman transformation).
There is also the possibility to enter a \code{lavaan} single factor solution.
In this case, \code{g_name} is not needed. Finally, if a solution from a
\code{lavaan} multiple group analysis is entered, the indices are computed for
each group.
The type argument is not evaluated if \code{model} is of class
\code{lavaan}.

If \code{model} is of class \code{\link{SL}} or
\code{\link[psych:schmid]{psych::schmid}} only the
\code{type} and, depending on the type (see below), the \code{factor_corres}
arguments need to be specified additionally. If model is of class
\code{\link[psych:schmid]{psych::schmid}} and \code{variance = "correlation"}
(default), it is
recommended to also provide the original correlation matrix in \code{cormat}
to get more accurate results. Otherwise, the correlation matrix will be found
based on the pattern matrix and Phi from the
\code{\link[psych:schmid]{psych::schmid}} output
using the \code{\link[psych:factor.model]{psych::factor.model}} function.

If \code{model = NULL}, the arguments \code{type}, \code{factor_corres}
(depending on the type, see below), \code{var_names}, \code{g_load}, \code{s_load},
and \code{u2} and either \code{cormat} (recommended) or \code{Phi} and
\code{pattern} need to be specified. If \code{Phi} and \code{pattern} are
specified instead of \code{cormat}, the correlation matrix is found using
the \code{\link[psych:factor.model]{psych::factor.model}} function.

The only difference between \code{type = "EFAtools"} and \code{type = "psych"}
is the determination of variable-to-factor correspondences. \code{type = "psych"}
reproduces the \code{\link[psych:omega]{psych::omega}} results, where
variable-to-factor correspondences are found by taking the highest
group factor loading for each variable as the relevant group factor loading.
To do this, \code{factor_corres} must be left \code{NULL}.

The calculation of the total variance (for the whole scale as well as the
subscale composites) can also be controlled in this function using the
\code{variance} argument. For both types---\code{"EFAtools"} and \code{"psych"}
---\code{variance} is set to \code{"correlation"} by default, which means that
total variances are found using the correlation matrix. If
\code{variance = "sums_load"} the total variance is calculated using the
squared sums of general loadings and group factor loadings and the sum of the
uniquenesses. This will only get comparable results to
\code{variance = "correlation"} if no cross-loadings are present and simple
structure is well-achieved in general with the SL solution (i.e., the
uniquenesses should capture almost all of the variance not explained by the
general factor and the variable's allocated group factor).
}
\examples{
\donttest{
## Use with lavaan outputs

# Create and fit bifactor model in lavaan (assume all variables have SDs of 1)
mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
fit_bi <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                      sample.nobs = 500, estimator = "ml", orthogonal = TRUE)

# Compute omegas and additional indices for bifactor solution
OMEGA(fit_bi, g_name = "g")

# Compute only omegas
OMEGA(fit_bi, g_name = "g", add_ind = FALSE)

# Create and fit second-order model in lavaan (assume all variables have SDs of 1)
mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ F1 + F2 + F3'
fit_ho <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                      sample.nobs = 500, estimator = "ml")

# Compute omegas and additional indices for second-order solution
OMEGA(fit_ho, g_name = "g")
}

## Use with an output from the SL function, with type EFAtools
efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")

# Two examples how to specify the indicator-to-factor correspondences:

# Based on a specific salience threshold for the loadings (here: .20):
factor_corres_1 <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2

# Or more flexibly (could also be TRUE and FALSE instead of 0 and 1):
factor_corres_2 <- matrix(c(rep(0, 12), rep(1, 6), rep(0, 6), rep(1, 6),
                         rep(0, 6), rep(1, 6), rep(0, 12)), ncol = 3,
                         byrow = FALSE)

OMEGA(sl_mod, type = "EFAtools", factor_corres = factor_corres_1)

## Use with an output from the psych::schmid function, with type psych for
## OMEGA
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")
# Find correlation matrix from phi and pattern matrix from psych::schmid output
OMEGA(schmid_mod, type = "psych")
# Use specified correlation matrix
OMEGA(schmid_mod, type = "psych", cormat = test_models$baseline$cormat)

## Manually specify components (useful if omegas should be computed for a SL
## or bifactor solution found with another program)
## As an example, we extract the elements from an SL output here. This gives
## the same results as in the second example above.

efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")

factor_corres <- matrix(c(rep(0, 12), rep(1, 6), rep(0, 6), rep(1, 6),
                        rep(0, 6), rep(1, 6), rep(0, 12)), ncol = 3,
                        byrow = FALSE)

OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
      g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
      u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
      factor_corres = factor_corres)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPSS_27_doc.R
\docType{data}
\name{SPSS_27}
\alias{SPSS_27}
\title{Various outputs from SPSS (version 27) FACTOR}
\format{
A list of 9 containing EFA results for each of the data sets mentioned above. Each of these nine entries is a list of 4 or 8 (see details), of the following structure:
\describe{
  \item{paf_comm}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details, see Grieder and Grob (2019).}
  \item{paf_load}{(matrix) - F1 to FN = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
  \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
  \item{var_load}{(matrix) - F1 to FN = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
  \item{pro_load}{(matrix) - F1 to FN = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
  \item{pro_phi}{(matrix) - F1 to FN = intercorrelations of the promax rotated loadings.}
  \item{sl}{(matrix) - g = General / second order factor of the Schmid-Leiman solution. F1 to FN  = First order factors of the Schmid-Leiman solution. h2 = Communalities of the Schmid-Leiman solution. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
  \item{L2}{(matrix) - Second order loadings used for the Schmid-Leiman transformation. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
 }
}
\source{
Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
 in R and SPSS. Manuscript in Preparation.

Wolff, H.G., & Preising, K. (2005). Exploring item and higher order factor structure with the Schmid-Leiman solution: Syntax codes for SPSS and SAS. Behavior Research Methods, 37, 48–58. doi: 10.3758/BF03206397

Grieder, S., & Grob, A. (2019). Exploratory factor analyses of the intelligence and development scales–2: Implications for theory and practice. Assessment. Advance online publication. doi:10.1177/10731911198450

Grob, A., & Hagmann-von Arx, P. (2018). Intelligence and Development Scales--2 (IDS-2). Intelligenz- und Entwicklungsskalen für Kinder und Jugendliche.
[Intelligence and Development Scales for Children and Adolescents.]. Bern, Switzerland: Hogrefe.

Frey, R., Pedroni, A., Mata, R., Rieskamp, J., & Hertwig, R. (2017). Risk preference shares the psychometric structure of major psychological traits. Science Advances, 3, e1701381.

McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.

Costa, P. T., & McCrae, R. R. (1992). NEO PI-R professional manual. Odessa, FL: Psychological Assessment Resources, Inc.
}
\usage{
SPSS_27
}
\description{
Various outputs from SPSS (version 27) FACTOR for the IDS-2 (Grob & Hagmann-von Arx, 2018), the WJIV (3 to 5 and 20 to 39 years; McGrew, LaForte, & Schrank, 2014), the DOSPERT (Frey et al., 2017; Weber,
Blais, & Betz, 2002), the NEO-PI-R (Costa, & McCrae, 1992), and four simulated datasets (baseline, case_1a, case_6b, and case_11b, see \link{test_models} and \link{population_models}) used in Grieder and Steiner (2020).
}
\details{
The IDS-2, the two WJIV, the DOSPERT, and the NEO-PI-R contain all the above entries, while the four simulated datasets contain only paf_load, var_load, pro_load, and pro_phi.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_models_doc.R
\docType{data}
\name{population_models}
\alias{population_models}
\title{population_models}
\format{
A list of 3 lists "loadings", "phis_3", and "phis_6".
\describe{
\code{loadings} contains the following matrices of pattern coefficients:
  \item{baseline}{(matrix) - The pattern coefficients of the baseline model. Three factors with six indicators each, all with pattern coefficients of .6. Same baseline model as used in de Winter and Dodou (2012).}
  \item{case_1a}{(matrix) - Three factors with 2 indicators per factor.}
  \item{case_1b}{(matrix) - Three factors with 3 indicators per factor. Case 5 in de Winter and Dodou (2012).}
  \item{case_1c}{(matrix) - Three factors with 4 indicators per factor.}
  \item{case_1d}{(matrix) - Three factors with 5 indicators per factor.}
  \item{case_2}{(matrix) - Same as baseline model but with low pattern coefficients of .3.}
  \item{case_3}{(matrix) - Same as baseline model but with high pattern coefficients of .9.}
  \item{case_4}{(matrix) - Three factors with different pattern coefficients \emph{between} factors (one factor with .9, one with .6, and one with .3, respectively). Case 7 in de Winter and Dodou (2012).}
  \item{case_5}{(matrix) - Three factors with different pattern coefficients \emph{within} factors (each factor has two pattern coefficients of each .9, .6, and .3). Similar to cases 8/ 9 in de Winter and Dodou (2012).}
  \item{case_6a}{(matrix) - Same as baseline model but with one cross loading of .4. Similar to case 10 in de Winter and Dodou (2012).}
  \item{case_6b}{(matrix) - Same as baseline model but with three cross loading of .4 (One factor with 2 and one with 1 crossloading). Similar to case 10 in de Winter and Dodou (2012).}
  \item{case_7}{(matrix) - Three factors with different number of indicators per factor (2, 4, and 6 respectively). Similar to cases 11/ 12 in de Winter and Dodou (2012).}
  \item{case_8}{(matrix) - Three factors with random variation in pattern coefficients added, drawn from a uniform distribution between [-.2, .2]. Case 13 in de Winter and Dodou (2012).}
  \item{case_9a}{(matrix) - Three factors with 2 indicators per factor, with different pattern coefficients within one of the factors.}
  \item{case_9b}{(matrix) - Three factors with 3 indicators per factor, with different pattern coefficients.}
  \item{case_9c}{(matrix) - Three factors with 4 indicators per factor, with different pattern coefficients.}
  \item{case_9d}{(matrix) - Three factors with 5 indicators per factor, with different pattern coefficients.}
  \item{case_10a}{(matrix) - Six factors with 2 indicators per factor, all with pattern coefficients of .6.}
  \item{case_10b}{(matrix) - Six factors with 3 indicators per factor, all with pattern coefficients of .6.}
  \item{case_10c}{(matrix) - Six factors with 4 indicators per factor, all with pattern coefficients of .6.}
  \item{case_10d}{(matrix) - Six factors with 5 indicators per factor, all with pattern coefficients of .6.}
  \item{case_10e}{(matrix) - Six factors with 6 indicators per factor, all with pattern coefficients of .6.}
  \item{case_11a}{(matrix) - Six factors with 2 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
  \item{case_11b}{(matrix) - Six factors with 3 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
  \item{case_11c}{(matrix) - Six factors with 4 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
  \item{case_11d}{(matrix) - Six factors with 5 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
  \item{case_11e}{(matrix) - Six factors with 6 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
  \item{case_12a}{(matrix) - One factor, with 2 equal pattern coefficients (.6).}
  \item{case_12b}{(matrix) - One factor, with 3 equal pattern coefficients (.6).}
  \item{case_12c}{(matrix) - One factor, with 6 equal pattern coefficients (.6).}
  \item{case_12d}{(matrix) - One factor, with 10 equal pattern coefficients (.6).}
  \item{case_12e}{(matrix) - One factor, with 15 equal pattern coefficients (.6).}
  \item{case_13a}{(matrix) - One factor, with 2 different pattern coefficients (.3, and .6).}
  \item{case_13b}{(matrix) - One factor, with 3 different pattern coefficients (.3, .6, and .9).}
  \item{case_13c}{(matrix) - One factor, with 6 different pattern coefficients (.3, .6, and .9).}
  \item{case_13d}{(matrix) - One factor, with 10 different pattern coefficients (.3, .6, and .9).}
  \item{case_13e}{(matrix) - One factor, with 15 different pattern coefficients (.3, .6, and .9).}
  \item{case_14a}{(matrix) - No factor, 2 variables (0).}
  \item{case_14b}{(matrix) - No factor, 3 variables (0).}
  \item{case_14c}{(matrix) - No factor, 6 variables (0).}
  \item{case_14d}{(matrix) - No factor, 10 variables (0).}
  \item{case_14e}{(matrix) - No factor, 15 variables (0).}
  \code{phis_3} contains the following 3x3 matrices:
  \item{zero}{(matrix) - Matrix of factor intercorrelations of 0. Same intercorrelations as used in de Winter and Dodou (2012).}
  \item{moderate}{(matrix) - Matrix of moderate factor intercorrelations of .3.}
  \item{mixed}{(matrix) - Matrix of mixed (.3, .5, and .7) factor intercorrelations.}
  \item{strong}{(matrix) - Matrix of strong factor intercorrelations of .7. Same intercorrelations as used in de Winter and Dodou (2012).}
  \code{phis_6} contains the following 6x6 matrices:
  \item{zero}{(matrix) - Matrix of factor intercorrelations of 0. Same intercorrelations as used in de Winter and Dodou (2012).}
  \item{moderate}{(matrix) - Matrix of moderate factor intercorrelations of .3.}
  \item{mixed}{(matrix) - Matrix of mixed (around .3, .5, and .7; smoothing was necessary for the matrix to be positive definite) factor intercorrelations.}
  \item{strong}{(matrix) - Matrix of strong factor intercorrelations of .7. Same intercorrelations as used in de Winter and Dodou (2012).}
 }
}
\source{
Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
 in R and SPSS. Manuscript in Preparation.

de Winter, J.C.F., & Dodou, D. (2012). Factor recovery by principal axis factoring and maximum likelihood factor analysis as a function of factor pattern and sample size. Journal of Applied Statistics. 39.
}
\usage{
population_models
}
\description{
Population factor models, some of which (baseline to case_11e) used for the
simulation analyses reported in Grieder and Steiner (2019). All combinations
of the pattern matrices and the factor
intercorrelations were used in the simulations. Many models are based on cases
used in de Winter and Dodou (2012).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WJIV_ages_14_19_doc.R
\docType{data}
\name{WJIV_ages_14_19}
\alias{WJIV_ages_14_19}
\title{Woodcock Johnson IV: ages 14 to 19}
\format{
A list of 2 with elements "cormat" (47 x 47 matrix of bivariate correlations)
and "N" (scalar). The correlation matrix contains the following variables:
\describe{
  \item{ORLVOC}{(numeric) - Oral Vocabulary.}
  \item{NUMSER}{(numeric) - Number Series.}
  \item{VRBATN}{(numeric) - Verbal Attention.}
  \item{LETPAT}{(numeric) - Letter-Pattern Matching.}
  \item{PHNPRO}{(numeric) - Phonological Processing.}
  \item{STYREC}{(numeric) - Story Recall.}
  \item{VISUAL}{(numeric) - Visualization.}
  \item{GENINF}{(numeric) - General Information.}
  \item{CONFRM}{(numeric) - Concept Formation.}
  \item{NUMREV}{(numeric) - Numbers Reversed.}
  \item{NUMPAT}{(numeric) - Number-Pattern Matching.}
  \item{NWDREP}{(numeric) - Nonword Repetition.}
  \item{VAL}{(numeric) - Visual-Auditory Learning.}
  \item{PICREC}{(numeric) - Picture Recognition.}
  \item{ANLSYN}{(numeric) - Analysis-Synthesis.}
  \item{OBJNUM}{(numeric) - Object-Number Sequencing.}
  \item{PAIRCN}{(numeric) - Pair Cancellation.}
  \item{MEMWRD}{(numeric) - Memory for Words.}
  \item{PICVOC}{(numeric) - Picture Vocabulary.}
  \item{ORLCMP}{(numeric) - Oral Comprehension.}
  \item{SEGMNT}{(numeric) - Segmentation.}
  \item{RPCNAM}{(numeric) - Rapid Picture Naming.}
  \item{SENREP}{(numeric) - Sentence Repetition.}
  \item{UNDDIR}{(numeric) - Understanding Directions.}
  \item{SNDBLN}{(numeric) - Sound Blending.}
  \item{RETFLU}{(numeric) - Retrieval Fluency.}
  \item{SNDAWR}{(numeric) - Sound Awareness.}
  \item{LWIDNT}{(numeric) - Letter-Word Identification.}
  \item{APPROB}{(numeric) - Applied Problems.}
  \item{SPELL}{(numeric) - Spelling.}
  \item{PSGCMP}{(numeric) - Passage Comprehension.}
  \item{CALC}{(numeric) - Calculation.}
  \item{WRTSMP}{(numeric) - Writing Samples.}
  \item{WRDATK}{(numeric) - Word Attack.}
  \item{ORLRDG}{(numeric) - Oral Reading.}
  \item{SNRDFL}{(numeric) - Sentence Reading Fluency.}
  \item{MTHFLU}{(numeric) - Math Facts Fluency.}
  \item{SNWRFL}{(numeric) - Sentence Writing Fluency.}
  \item{RDGREC}{(numeric) - Reading Recall.}
  \item{NUMMAT}{(numeric) - Number Matrices.}
  \item{EDIT}{(numeric) - Editing.}
  \item{WRDFLU}{(numeric) - Word Reading Fluency.}
  \item{SPLSND}{(numeric) - Spelling of Sounds.}
  \item{RDGVOC}{(numeric) - Reading Vocabulary.}
  \item{SCI}{(numeric) - Science.}
  \item{SOC}{(numeric) - Social Studies.}
  \item{HUM}{(numeric) - Humanities.}
 }
}
\source{
McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.
}
\usage{
WJIV_ages_14_19
}
\description{
A list containing the bivariate correlations (N = 1,685) of the 47
cognitive and achievement subtests from the WJ IV for 14- to 19-year-olds
from the standardization sample obtained from
the WJ-IV technical manual (McGrew, LaForte, & Schrank, 2014).
Tables are reproduced with permission from the publisher.
}
\keyword{datasets}
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
% Please edit documentation in R/WJIV_ages_9_13_doc.R
\docType{data}
\name{WJIV_ages_9_13}
\alias{WJIV_ages_9_13}
\title{Woodcock Johnson IV: ages 9 to 13}
\format{
A list of 2 with elements "cormat" (47 x 47 matrix of bivariate
correlations) and "N". The correlation matrix contains the following variables:
\describe{
  \item{ORLVOC}{(numeric) - Oral Vocabulary.}
  \item{NUMSER}{(numeric) - Number Series.}
  \item{VRBATN}{(numeric) - Verbal Attention.}
  \item{LETPAT}{(numeric) - Letter-Pattern Matching.}
  \item{PHNPRO}{(numeric) - Phonological Processing.}
  \item{STYREC}{(numeric) - Story Recall.}
  \item{VISUAL}{(numeric) - Visualization.}
  \item{GENINF}{(numeric) - General Information.}
  \item{CONFRM}{(numeric) - Concept Formation.}
  \item{NUMREV}{(numeric) - Numbers Reversed.}
  \item{NUMPAT}{(numeric) - Number-Pattern Matching.}
  \item{NWDREP}{(numeric) - Nonword Repetition.}
  \item{VAL}{(numeric) - Visual-Auditory Learning.}
  \item{PICREC}{(numeric) - Picture Recognition.}
  \item{ANLSYN}{(numeric) - Analysis-Synthesis.}
  \item{OBJNUM}{(numeric) - Object-Number Sequencing.}
  \item{PAIRCN}{(numeric) - Pair Cancellation.}
  \item{MEMWRD}{(numeric) - Memory for Words.}
  \item{PICVOC}{(numeric) - Picture Vocabulary.}
  \item{ORLCMP}{(numeric) - Oral Comprehension.}
  \item{SEGMNT}{(numeric) - Segmentation.}
  \item{RPCNAM}{(numeric) - Rapid Picture Naming.}
  \item{SENREP}{(numeric) - Sentence Repetition.}
  \item{UNDDIR}{(numeric) - Understanding Directions.}
  \item{SNDBLN}{(numeric) - Sound Blending.}
  \item{RETFLU}{(numeric) - Retrieval Fluency.}
  \item{SNDAWR}{(numeric) - Sound Awareness.}
  \item{LWIDNT}{(numeric) - Letter-Word Identification.}
  \item{APPROB}{(numeric) - Applied Problems.}
  \item{SPELL}{(numeric) - Spelling.}
  \item{PSGCMP}{(numeric) - Passage Comprehension.}
  \item{CALC}{(numeric) - Calculation.}
  \item{WRTSMP}{(numeric) - Writing Samples.}
  \item{WRDATK}{(numeric) - Word Attack.}
  \item{ORLRDG}{(numeric) - Oral Reading.}
  \item{SNRDFL}{(numeric) - Sentence Reading Fluency.}
  \item{MTHFLU}{(numeric) - Math Facts Fluency.}
  \item{SNWRFL}{(numeric) - Sentence Writing Fluency.}
  \item{RDGREC}{(numeric) - Reading Recall.}
  \item{NUMMAT}{(numeric) - Number Matrices.}
  \item{EDIT}{(numeric) - Editing.}
  \item{WRDFLU}{(numeric) - Word Reading Fluency.}
  \item{SPLSND}{(numeric) - Spelling of Sounds.}
  \item{RDGVOC}{(numeric) - Reading Vocabulary.}
  \item{SCI}{(numeric) - Science.}
  \item{SOC}{(numeric) - Social Studies.}
  \item{HUM}{(numeric) - Humanities.}
 }
}
\source{
McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.
}
\usage{
WJIV_ages_9_13
}
\description{
A list containing the bivariate correlations (N = 1,572) of the 47
cognitive and achievement subtests from the WJ IV for 9- to 13-year-olds from
the standardization sample obtained from
the WJ-IV technical manual (McGrew, LaForte, & Schrank, 2014).
Tables are reproduced with permission from the publisher.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WJIV_ages_6_8_doc.R
\docType{data}
\name{WJIV_ages_6_8}
\alias{WJIV_ages_6_8}
\title{Woodcock Johnson IV: ages 6 to 8}
\format{
A list of 2 with elements "cormat" (47 x 47 matrix of bivariate
correlations) and "N". The correlation matrix contains the following
variables:
\describe{
  \item{ORLVOC}{(numeric) - Oral Vocabulary.}
  \item{NUMSER}{(numeric) - Number Series.}
  \item{VRBATN}{(numeric) - Verbal Attention.}
  \item{LETPAT}{(numeric) - Letter-Pattern Matching.}
  \item{PHNPRO}{(numeric) - Phonological Processing.}
  \item{STYREC}{(numeric) - Story Recall.}
  \item{VISUAL}{(numeric) - Visualization.}
  \item{GENINF}{(numeric) - General Information.}
  \item{CONFRM}{(numeric) - Concept Formation.}
  \item{NUMREV}{(numeric) - Numbers Reversed.}
  \item{NUMPAT}{(numeric) - Number-Pattern Matching.}
  \item{NWDREP}{(numeric) - Nonword Repetition.}
  \item{VAL}{(numeric) - Visual-Auditory Learning.}
  \item{PICREC}{(numeric) - Picture Recognition.}
  \item{ANLSYN}{(numeric) - Analysis-Synthesis.}
  \item{OBJNUM}{(numeric) - Object-Number Sequencing.}
  \item{PAIRCN}{(numeric) - Pair Cancellation.}
  \item{MEMWRD}{(numeric) - Memory for Words.}
  \item{PICVOC}{(numeric) - Picture Vocabulary.}
  \item{ORLCMP}{(numeric) - Oral Comprehension.}
  \item{SEGMNT}{(numeric) - Segmentation.}
  \item{RPCNAM}{(numeric) - Rapid Picture Naming.}
  \item{SENREP}{(numeric) - Sentence Repetition.}
  \item{UNDDIR}{(numeric) - Understanding Directions.}
  \item{SNDBLN}{(numeric) - Sound Blending.}
  \item{RETFLU}{(numeric) - Retrieval Fluency.}
  \item{SNDAWR}{(numeric) - Sound Awareness.}
  \item{LWIDNT}{(numeric) - Letter-Word Identification.}
  \item{APPROB}{(numeric) - Applied Problems.}
  \item{SPELL}{(numeric) - Spelling.}
  \item{PSGCMP}{(numeric) - Passage Comprehension.}
  \item{CALC}{(numeric) - Calculation.}
  \item{WRTSMP}{(numeric) - Writing Samples.}
  \item{WRDATK}{(numeric) - Word Attack.}
  \item{ORLRDG}{(numeric) - Oral Reading.}
  \item{SNRDFL}{(numeric) - Sentence Reading Fluency.}
  \item{MTHFLU}{(numeric) - Math Facts Fluency.}
  \item{SNWRFL}{(numeric) - Sentence Writing Fluency.}
  \item{RDGREC}{(numeric) - Reading Recall.}
  \item{NUMMAT}{(numeric) - Number Matrices.}
  \item{EDIT}{(numeric) - Editing.}
  \item{WRDFLU}{(numeric) - Word Reading Fluency.}
  \item{SPLSND}{(numeric) - Spelling of Sounds.}
  \item{RDGVOC}{(numeric) - Reading Vocabulary.}
  \item{SCI}{(numeric) - Science.}
  \item{SOC}{(numeric) - Social Studies.}
  \item{HUM}{(numeric) - Humanities.}
 }
}
\source{
McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.
}
\usage{
WJIV_ages_6_8
}
\description{
A list containing the bivariate correlations (N = 825) of the 47
cognitive and achievement subtests from the WJ IV for 6- to 8-year-olds from
the standardization sample obtained from
the WJ-IV technical manual (McGrew, LaForte, & Schrank, 2014).
Tables are reproduced with permission from the publisher.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRiPS_raw_doc.R
\docType{data}
\name{GRiPS_raw}
\alias{GRiPS_raw}
\title{GRiPS_raw}
\format{
An object of class \code{data.frame} with 810 rows and 8 columns.
}
\source{
Zhang, D. C., Highhouse, S., & Nye, C. D. (2018). Development and validation of the general risk propensity scale (GRiPS).Journal of Behavioral Decision Making, 32, 152–167. doi: 10.1002/bdm.2102

Steiner, M., & Frey, R. (2020). Representative design in psychological assessment: A case study using the Balloon Analogue Risk Task (BART). PsyArXiv Preprint. doi:10.31234/osf.io/dg4ks
}
\usage{
GRiPS_raw
}
\description{
A data.frame containing responses to the General Risk Propensity Scale (GRiPS, Zhang,
Highhouse & Nye, 2018) of 810 participants of Study 1 of Steiner and Frey (2020).
The original data can be accessed via \url{https://osf.io/kxp8t/}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.EKC.R
\name{print.EKC}
\alias{print.EKC}
\title{Print function for EKC objects}
\usage{
\method{print}{EKC}(x, plot = TRUE, ...)
}
\arguments{
\item{x}{a list of class EKC. Output from \code{\link{EKC}} function.}

\item{plot}{logical. Whether to plot the results.}

\item{...}{Further arguments for print.}
}
\description{
Print function for EKC objects
}
\examples{
EKC_base <- EKC(test_models$baseline$cormat, N = 500)
EKC_base

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EFA.R
\name{EFA}
\alias{EFA}
\title{Exploratory factor analysis (EFA)}
\source{
Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
 in R and SPSS. Manuscript in Preparation.

Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
rotation to oblique simple structure. British Journal of Statistical Psychology,
17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The
Hull Method for Selecting the Number of Common Factors, Multivariate Behavioral
Research, 46, 340-364, doi: 10.1080/00273171.2011.564527

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233
}
\usage{
EFA(
  x,
  n_factors,
  N = NA,
  method = c("PAF", "ML", "ULS"),
  rotation = c("none", "varimax", "equamax", "quartimax", "geominT", "bentlerT",
    "bifactorT", "promax", "oblimin", "quartimin", "simplimax", "bentlerQ", "geominQ",
    "bifactorQ"),
  type = c("EFAtools", "psych", "SPSS", "none"),
  max_iter = NA,
  init_comm = NA,
  criterion = NA,
  criterion_type = NA,
  abs_eigen = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  varimax_type = NA,
  k = NA,
  normalize = TRUE,
  P_type = NA,
  precision = 1e-05,
  order_type = NA,
  start_method = "psych",
  cor_method = c("pearson", "spearman", "kendall"),
  ...
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations. If raw data is entered, the correlation matrix is found from the
data.}

\item{n_factors}{numeric. Number of factors to extract.}

\item{N}{numeric. The number of observations. Needs only be specified if a
correlation matrix is used. If input is a correlation matrix and \code{N} = NA
(default), not all fit indices can be computed.}

\item{method}{character. One of "PAF", "ML", or "ULS" to use principal axis
factoring, maximum likelihood, or unweighted least squares (also called minres),
respectively, to fit the EFA.}

\item{rotation}{character. Either perform no rotation ("none"; default),
an orthogonal rotation ("varimax", "equamax", "quartimax", "geominT",
"bentlerT", or "bifactorT"), or an oblique rotation ("promax", "oblimin",
"quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ").}

\item{type}{character. If one of "EFAtools" (default), "psych", or "SPSS" is
used, and the following arguments with default NA are left with
NA, these implementations are executed according to the respective program
("psych" and "SPSS") or according to the best solution found in Grieder &
Steiner (2020; "EFAtools"). Individual properties can be adapted using one of
the three types and specifying some of the following arguments. If set to
"none" additional arguments must be specified depending on the \code{method}
and \code{rotation} used (see details).}

\item{max_iter}{numeric. The maximum number of iterations to perform after which
the iterative PAF procedure is halted with a warning. If \code{type} is one of
"EFAtools", "SPSS", or "psych", this is automatically specified if \code{max_iter} is
left to be \code{NA}, but can be overridden by entering a number. Default is
\code{NA}.}

\item{init_comm}{character. The method to estimate the initial communalities
in \code{PAF}. "smc" will use squared multiple correlations, "mac" will use
maximum absolute correlations, "unity" will use 1s (see details).
Default is \code{NA}.}

\item{criterion}{numeric. The convergence criterion used for PAF.
If the change in communalities from one iteration to the next is smaller than
this criterion the solution is accepted and the procedure ends.
Default is \code{NA}.}

\item{criterion_type}{character. Type of convergence criterion used for
PAF. "max_individual" selects the maximum change in any of the
communalities from one iteration to the next and tests it against the
specified criterion. This is also used by SPSS. "sum" takes the difference of
the sum of all communalities in one iteration and the sum of all communalities
in the next iteration and tests this against the criterion. This procedure is
used by the \code{\link[psych:fa]{psych::fa}} function. Default is \code{NA}.}

\item{abs_eigen}{logical. Which algorithm to use in the PAF
iterations. If FALSE, the loadings are computed from the eigenvalues. This is
also used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
loadings are computed with the absolute eigenvalues as done by SPSS.
Default is \code{NA}.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
is given as input. Default is "pairwise.complete.obs".}

\item{varimax_type}{character. The type of the varimax rotation performed.
If "svd", singular value decomposition is used, as \link[stats:varimax]{stats::varimax} does. If "kaiser", the varimax procedure performed in SPSS is used.
This is the original procedure from Kaiser (1958), but with slight alterations
in the varimax criterion (see details, and Grieder & Steiner, 2020). Default is \code{NA}.}

\item{k}{numeric. Either the power used for computing the target matrix P in
the promax rotation or the number of 'close to zero loadings' for the simplimax
rotation (see \code{\link[GPArotation:GPA]{GPArotation::GPFoblq}}). If left to
\code{NA} (default), the value for promax depends on the specified type.
For simplimax, \code{nrow(L)}, where L is the matrix of unrotated loadings,
is used by default.}

\item{normalize}{logical. If \code{TRUE}, a kaiser normalization is
performed before the specified rotation. Default is \code{TRUE}.}

\item{P_type}{character. This specifies how the target
matrix P is computed in promax rotation. If "unnorm" it will use the
unnormalized target matrix as originally done in Hendrickson and White (1964).
This is also used in the psych and stats packages. If "norm" it will use the
normalized target matrix as used in SPSS. Default is \code{NA}.}

\item{precision}{numeric. The tolerance for stopping in the rotation
procedure. Default is 10^-5 for all rotation methods.}

\item{order_type}{character. How to order the factors. "eigen" will reorder
the factors according to the largest to lowest eigenvalues of the matrix of
rotated loadings. "ss_factors" will reorder the factors according to descending
sum of squared factor loadings per factor. Default is \code{NA}.}

\item{start_method}{character. How to specify the starting values for the
optimization procedure for ML. Default is "psych" which takes the
starting values specified in \link[psych:fa]{psych::fa}. "factanal" takes the
starting values specified in the \link[stats:factanal]{stats::factanal} function.
Solutions are very similar.}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}

\item{...}{Additional arguments passed to rotation functions from the \code{GPArotation} package (e.g., \code{maxit} for maximum number of iterations).}
}
\value{
A list of class EFA containing (a subset of) the following:

\item{orig_R}{Original correlation matrix.}
\item{h2_init}{Initial communality estimates from PAF.}
\item{h2}{Final communality estimates from the unrotated solution.}
\item{orig_eigen}{Eigen values of the original correlation matrix.}
\item{init_eigen}{Initial eigenvalues, obtained from the correlation matrix
 with the initial communality estimates as diagonal in PAF.}
\item{final_eigen}{Eigenvalues obtained from the correlation matrix
 with the final communality estimates as diagonal.}
\item{iter}{The number of iterations needed for convergence.}
\item{convergence}{Integer code for convergence as returned by
\code{\link[stats:optim]{stats:optim}} (only for ML and ULS).
0 indicates successful completion.}
\item{unrot_loadings}{Loading matrix containing the final unrotated loadings.}
\item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based on the unrotated loadings.}
\item{fit_indices}{For ML and ULS: Fit indices derived from the unrotated
factor loadings: Chi Square, including significance level, degrees of freedom
(df), Comparative Fit Index (CFI), Root Mean Square Error of Approximation
(RMSEA), including its 90\% confidence interval, Akaike Information Criterion
(AIC), Bayesian Information Criterion (BIC), and the common part accounted
for (CAF) index as proposed by Lorenzo-Seva, Timmerman, & Kiers (2011).
For PAF, only the CAF and dfs are returned.}
\item{rot_loadings}{Loading matrix containing the final rotated loadings
(pattern matrix).}
\item{Phi}{The factor intercorrelations (only for oblique rotations).}
\item{Structure}{The structure matrix (only for oblique rotations).}
\item{rotmat}{The rotation matrix.}
\item{vars_accounted_rot}{Matrix of explained variances and sums of squared
loadings. Based on rotated loadings and, for oblique rotations, the factor
intercorrelations.}
\item{settings}{A list of the settings used.}
}
\description{
This function does an EFA with either \code{PAF}, \code{ML},
or \code{ULS} with or without subsequent rotation.
All arguments with default value \code{NA} can be left to default if \code{type}
is set to one of "EFAtools", "SPSS", or "psych". The respective specifications are
then handled according to the specified type (see details). For all rotations
except varimax and promax, the \code{GPArotation} package is needed.
}
\details{
There are two main ways to use this function. The easiest way is to
use it with a specified \code{type} (see above), which sets most of the other
arguments accordingly. Another way is to use it more flexibly by explicitly
specifying all arguments used and set \code{type} to "none" (see examples).
A mix of the two can also be done by specifying a \code{type} as well as
additional arguments. However, this will throw warnings to avoid unintentional
deviations from the implementations according to the specified \code{type}.

The \code{type} argument is evaluated for PAF and for all rotations (mainly
important for the varimax and promax rotations). The type-specific settings
for these functions are detailed below.

For PAF, the values of \code{init_comm}, \code{criterion}, \code{criterion_type},
and \code{abs_eigen} depend on the \code{type} argument.

\code{type = "EFAtools"} will use the following argument specification:
\code{init_comm = "smc", criterion = .001, criterion_type = "sum",
abs_eigen = TRUE}.

\code{type = "psych"} will use the following argument specification:
\code{init_comm = "smc", criterion = .001, criterion_type = "sum",
abs_eigen = FALSE}.

\code{type = "SPSS"} will use the following argument specification:
\code{init_comm = "smc", criterion = .001, criterion_type = "max_individual",
abs_eigen = TRUE}.

If SMCs fail, SPSS takes "mac". However, as SPSS takes absolute eigenvalues,
this is hardly ever the case. Psych, on the other hand, takes "unity" if SMCs
fail, but uses the Moore-Penrose Psudo Inverse of a matrix, thus, taking "unity"
is only necessary if negative eigenvalues occur afterwards in the iterative
PAF procedure. The EFAtools type setting combination was the best in terms of accuracy
and number of Heywood cases compared to all the
other setting combinations tested in simulation studies in Grieder & Steiner
(2020), which is why this type is used as a default here.

For varimax, the values of \code{varimax_type} and \code{order_type} depend on
the \code{type} argument.

\code{type = "EFAtools"} will use the following argument specification:
\code{varimax_type = "kaiser", order_type = "eigen"}.

\code{type = "psych"} will use the following argument specification:
\code{varimax_type = "svd", order_type = "eigen"}.

\code{type = "SPSS"} will use the following argument specification:
\code{varimax_type = "kaiser", order_type = "ss_factors"}.

For promax, the values of \code{P_type},
\code{order_type}, and \code{k} depend on the \code{type} argument.

\code{type = "EFAtools"} will use the following argument specification:
\code{P_type = "norm", order_type = "eigen", k = 4}.

\code{type = "psych"} will use the following argument specification:
\code{P_type = "unnorm", order_type = "eigen", k = 4}.

\code{type = "SPSS"} will use the following argument specification:
\code{P_type = "norm", order_type = "ss_factors", k = 4}.

The \code{P_type} argument can take two values, "unnorm" and "norm". It controls
which formula is used to compute the target matrix P in the promax rotation.
"unnorm" uses the formula from Hendrickson and White (1964), specifically:
\code{P = abs(A^(k + 1)) / A},
where A is the unnormalized matrix containing varimax rotated loadings.
"SPSS" uses the normalized varimax rotated loadings. Specifically it used the
following formula, which can be found in the SPSS 23 and SPSS 27 Algorithms manuals:
\code{P = abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)}.
As for PAF, the EFAtools type setting combination for promax was the best
compared to the other setting combinations tested in simulation studies in
Grieder & Steiner (2020).

The \code{varimax_type} argument can take two values, "svd", and "kaiser". "svd" uses
singular value decomposition, by calling \link[stats:varimax]{stats::varimax}. "kaiser"
performs the varimax procedure as described in the SPSS 23 Algorithms manual and as described
by Kaiser (1958). However, there is a slight alteration in computing the varimax criterion, which
we found to better align with the results obtain from SPSS. Specifically, the original varimax
criterion as described in the SPSS 23 Algorithms manual is
\code{sum(n*colSums(lambda ^ 4) - colSums(lambda ^ 2) ^ 2) / n ^ 2}, where n is the
number of indicators, and lambda is the rotated loadings matrix. However, we found the following
to produce results more similar to those of SPSS:
\code{sum(n*colSums(abs(lambda)) - colSums(lambda ^ 4) ^ 2) / n^2}.

For all other rotations except varimax and promax, the \code{type} argument
only controls the \code{order_type} argument with the same values as stated
above for the varimax and promax rotations. For these other rotations, the
\code{GPArotation} package is needed. Additional arguments can also be
specified and will be passed to the respective \code{GPArotation} function
(e.g., maxit to change the maximum number of iterations for the rotation procedure).

The \code{type} argument has no effect on ULS and ML. For ULS, no additional
arguments are needed. For ML, an additional argument
\code{start_method} is needed to determine the starting values for the
optimization procedure. Default for this argument is "factanal" which takes
the starting values specified in the \link[stats:factanal]{stats::factanal} function.
}
\examples{
# A type EFAtools (as presented in Steiner and Grieder, 2020) EFA
EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                    type = "EFAtools", method = "PAF", rotation = "none")

# A type SPSS EFA to mimick the SPSS implementation (this will throw a warning,
# see below)
SPSS_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "SPSS", method = "PAF", rotation = "none")

# A type psych EFA to mimick the psych::fa() implementation
psych_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 type = "psych", method = "PAF", rotation = "none")

# Use ML instead of PAF with type EFAtools
EFAtools_ML <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                   type = "EFAtools", method = "ML", rotation = "none")

# Use oblimin rotation instead of no rotation with type EFAtools
EFAtools_oblim <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", method = "PAF", rotation = "oblimin")

# Do a PAF without rotation without specifying a type, so the arguments
# can be flexibly specified (this is only recommended if you know what your
# doing)
PAF_none <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "none", method = "PAF", rotation = "none",
                max_iter = 500, init_comm = "mac", criterion = 1e-4,
                criterion_type = "sum", abs_eigen = FALSE)

# Add a promax rotation
PAF_pro <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
               type = "none", method = "PAF", rotation = "promax",
               max_iter = 500, init_comm = "mac", criterion = 1e-4,
               criterion_type = "sum", abs_eigen = FALSE, k = 3,
               P_type = "unnorm", precision= 1e-5, order_type = "eigen",
               varimax_type = "svd")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.EFA_AVERAGE.R
\name{print.EFA_AVERAGE}
\alias{print.EFA_AVERAGE}
\title{Print EFA_AVERAGE object}
\usage{
\method{print}{EFA_AVERAGE}(x, stat = c("average", "range"), plot = TRUE, ...)
}
\arguments{
\item{x}{list. An object of class EFA_AVERAGE to be printed}

\item{stat}{character. A vector with the statistics to print. Possible inputs
are "average", "sd", "range", "min", and "max". Default is "average" and
"range".}

\item{plot}{logical. Whether a plot of the average and min- max loadings should
be created. Default is TRUE. If more than 10 factors are extracted, no plot is
created.}

\item{...}{Further arguments for print.}
}
\description{
Print Method showing a summarized output of the \link{EFA_AVERAGE} function
}
\examples{

EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
EFA_aver

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCREE.R
\name{SCREE}
\alias{SCREE}
\title{Scree Plot}
\source{
Cattell, R. B. (1966). The scree test for the number of factors.
Multivariate Behavioral Research, 1(2), 245–276.
https://doi.org/10.1207/s15327906mbr0102_10

Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain. Psychological Bulletin, 99,
432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
}
\usage{
SCREE(
  x,
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall"),
  n_factors = 1,
  ...
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations.}

\item{eigen_type}{character. On what the eigenvalues should be found. Can be
either "PCA", "SMC", or "EFA", or some combination of them. If using "PCA",
the diagonal values of the correlation matrices are left to be 1. If using
"SMC", the diagonal of the
correlation matrices is replaced by the squared multiple correlations (SMCs)
of the indicators. If using "EFA", eigenvalues are found on the correlation
matrices with the final communalities of an exploratory factor analysis
solution (default is principal axis factoring extracting 1 factor) as
diagonal.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
data is given as input. Default is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}

\item{n_factors}{numeric. Number of factors to extract if "EFA" is included in
\code{eigen_type}. Default is 1.}

\item{...}{Additional arguments passed to \code{\link{EFA}}. For example,
to change the extraction method (PAF is default).}
}
\value{
A list of class SCREE containing

\item{eigen_PCA}{ A vector containing the eigenvalues found with PCA.}
\item{eigen_SMC}{ A vector containing the eigenvalues found with SMCs.}
\item{eigen_EFA}{ A vector containing the eigenvalues found with EFA.}
\item{settings}{A list of the settings used.}
}
\description{
The scree plot was originally introduced by Cattell (1966) to perform the
scree test. In a scree plot, the eigenvalues of the factors / components are
plotted against the index of the factors / components, ordered from 1 to N
factors components, hence from largest to smallest eigenvalue. According to
the scree test, the number of factors / components to retain is the number of
factors / components to the left of the "elbow" (where the curve starts to
level off) in the scree plot.
}
\details{
As the scree test requires visual examination, the test has been
especially criticized for its subjectivity and with this low inter-rater
reliability. Moreover, a scree plot can be ambiguous if there are either no
clear "elbow" or multiple "elbows", making it difficult to judge just where
the eigenvalues do level off. Finally, the scree test has also been found to
be less accurate than other factor retention criteria. For all these reasons,
the scree test has been recommended against, at least for exclusive use as a
factor retention criterion (Zwick & Velicer, 1986)

The \code{SCREE} function can also be called together with other factor
retention criteria in the \code{\link{N_FACTORS}} function.
}
\examples{
SCREE(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
}
\seealso{
Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
\code{\link{HULL}}, \code{\link{PARALLEL}}, \code{\link{SMT}}

\code{\link{N_FACTORS}} as a wrapper function for this and all the
above-mentioned factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.EFA_AVERAGE.R
\name{plot.EFA_AVERAGE}
\alias{plot.EFA_AVERAGE}
\title{Plot EFA_AVERAGE object}
\usage{
\method{plot}{EFA_AVERAGE}(x, ...)
}
\arguments{
\item{x}{list. An output from the \link{EFA_AVERAGE} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \link{EFA_AVERAGE} function
}
\examples{
EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
EFA_aver

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CD.R
\name{CD}
\alias{CD}
\title{Comparison Data}
\source{
Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of extraction
methods under realistic conditions. Psychological Methods, 24(4), 468–491.
https://doi.org/10.1037/met0000200

Ruscio, J., & Roche, B. (2012). Determining the number of factors to
retain in an exploratory factor analysis using comparison data of known
factorial structure. Psychological Assessment, 24, 282–292.
doi: 10.1037/a0025697
}
\usage{
CD(
  x,
  n_factors_max = NA,
  N_pop = 10000,
  N_samples = 500,
  alpha = 0.3,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall"),
  max_iter = 50
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data.}

\item{n_factors_max}{numeric. The maximum number of factors to test against.
Larger numbers will increase the duration the procedure takes, but test more
possible solutions. If left NA (default) the maximum number of factors for
which the model is still over-identified (df > 0) is used.}

\item{N_pop}{numeric. Size of finite populations of comparison data. Default
is 10000.}

\item{N_samples}{numeric. Number of samples drawn from each population.
Default is 500.}

\item{alpha}{numeric. The alpha level used to test the significance of the
improvement added by an additional factor. Default is .30.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}}. Default
is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}

\item{max_iter}{numeric. The maximum number of iterations to perform after
which the iterative PAF procedure is halted. Default is 50.}
}
\value{
A list of class CD containing

\item{n_factors}{The number of factors to retain according to comparison data results.}
\item{eigenvalues}{A vector containing the eigenvalues of the entered data.}
\item{RMSE_eigenvalues}{A matrix containing the RMSEs between the eigenvalues of the generated data and those of the entered data.}
\item{settings}{A list of the settings used.}
}
\description{
Factor retention method introduced by Ruscio and Roche (2012). The code was
adapted from the CD code by Auerswald and Moshagen (2017) available at
\url{https://osf.io/x5cz2/?view_only=d03efba1fd0f4c849a87db82e6705668}
}
\details{
"Parallel analysis (PA) is an effective stopping rule that compares
the eigenvalues of randomly generated data with those for the actual data.
PA takes into account sampling error, and at present it is widely considered
the best available method. We introduce a variant of PA that goes even further
by reproducing the observed correlation matrix rather than generating random
data. Comparison data (CD) with known factorial structure are first generated
using 1 factor, and then the number of factors is increased until the
reproduction of the observed eigenvalues fails to improve significantly"
(Ruscio & Roche, 2012, p. 282).

The CD implementation here is based on the code by Ruscio and Roche (2012), but
is slightly adapted to increase speed by performing the principal axis factoring
using a C++ based function.

The \code{CD} function can also be called together with other factor retention
criteria in the \code{\link{N_FACTORS}} function.
}
\examples{
\donttest{
# determine n factors of the GRiPS
CD(GRiPS_raw)

# determine n factors of the DOSPERT risk subscale
CD(DOSPERT_raw)
}
}
\seealso{
Other factor retention criteria: \code{\link{EKC}},
 \code{\link{HULL}}, \code{\link{KGC}}, \code{\link{PARALLEL}}, \code{\link{SMT}}

  \code{\link{N_FACTORS}} as a wrapper function for this and all
  the above-mentioned factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WJIV_ages_40_90_doc.R
\docType{data}
\name{WJIV_ages_40_90}
\alias{WJIV_ages_40_90}
\title{Woodcock Johnson IV: ages 40 to 90 plus}
\format{
A list of 2 with elements "cormat" (47 x 47 matrix of bivariate correlations)
and "N". The correlation matrix contains the following variables:
\describe{
  \item{ORLVOC}{(numeric) - Oral Vocabulary.}
  \item{NUMSER}{(numeric) - Number Series.}
  \item{VRBATN}{(numeric) - Verbal Attention.}
  \item{LETPAT}{(numeric) - Letter-Pattern Matching.}
  \item{PHNPRO}{(numeric) - Phonological Processing.}
  \item{STYREC}{(numeric) - Story Recall.}
  \item{VISUAL}{(numeric) - Visualization.}
  \item{GENINF}{(numeric) - General Information.}
  \item{CONFRM}{(numeric) - Concept Formation.}
  \item{NUMREV}{(numeric) - Numbers Reversed.}
  \item{NUMPAT}{(numeric) - Number-Pattern Matching.}
  \item{NWDREP}{(numeric) - Nonword Repetition.}
  \item{VAL}{(numeric) - Visual-Auditory Learning.}
  \item{PICREC}{(numeric) - Picture Recognition.}
  \item{ANLSYN}{(numeric) - Analysis-Synthesis.}
  \item{OBJNUM}{(numeric) - Object-Number Sequencing.}
  \item{PAIRCN}{(numeric) - Pair Cancellation.}
  \item{MEMWRD}{(numeric) - Memory for Words.}
  \item{PICVOC}{(numeric) - Picture Vocabulary.}
  \item{ORLCMP}{(numeric) - Oral Comprehension.}
  \item{SEGMNT}{(numeric) - Segmentation.}
  \item{RPCNAM}{(numeric) - Rapid Picture Naming.}
  \item{SENREP}{(numeric) - Sentence Repetition.}
  \item{UNDDIR}{(numeric) - Understanding Directions.}
  \item{SNDBLN}{(numeric) - Sound Blending.}
  \item{RETFLU}{(numeric) - Retrieval Fluency.}
  \item{SNDAWR}{(numeric) - Sound Awareness.}
  \item{LWIDNT}{(numeric) - Letter-Word Identification.}
  \item{APPROB}{(numeric) - Applied Problems.}
  \item{SPELL}{(numeric) - Spelling.}
  \item{PSGCMP}{(numeric) - Passage Comprehension.}
  \item{CALC}{(numeric) - Calculation.}
  \item{WRTSMP}{(numeric) - Writing Samples.}
  \item{WRDATK}{(numeric) - Word Attack.}
  \item{ORLRDG}{(numeric) - Oral Reading.}
  \item{SNRDFL}{(numeric) - Sentence Reading Fluency.}
  \item{MTHFLU}{(numeric) - Math Facts Fluency.}
  \item{SNWRFL}{(numeric) - Sentence Writing Fluency.}
  \item{RDGREC}{(numeric) - Reading Recall.}
  \item{NUMMAT}{(numeric) - Number Matrices.}
  \item{EDIT}{(numeric) - Editing.}
  \item{WRDFLU}{(numeric) - Word Reading Fluency.}
  \item{SPLSND}{(numeric) - Spelling of Sounds.}
  \item{RDGVOC}{(numeric) - Reading Vocabulary.}
  \item{SCI}{(numeric) - Science.}
  \item{SOC}{(numeric) - Social Studies.}
  \item{HUM}{(numeric) - Humanities.}
 }
}
\source{
McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.
}
\usage{
WJIV_ages_40_90
}
\description{
A list containing the bivariate correlations (N = 1,146) of the 47
cognitive and achievement subtests from the WJ IV for 40- to 90+-year-olds
from the standardization sample obtained from
the WJ-IV technical manual (McGrew, LaForte, & Schrank, 2014).
Tables are reproduced with permission from the publisher.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.PARALLEL.R
\name{plot.PARALLEL}
\alias{plot.PARALLEL}
\title{Plot PARALLEL object}
\usage{
\method{plot}{PARALLEL}(x, ...)
}
\arguments{
\item{x}{list of class PARALLEL. An output from the \link{PARALLEL} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \link{PARALLEL} function
}
\examples{
\donttest{
# example with correlation matrix and "ML" estimation
x <- PARALLEL(test_models$case_11b$cormat, N = 500, method = "ML")
plot(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.parallel_sim}
\alias{.parallel_sim}
\title{Parallel analysis on simulated data.}
\usage{
.parallel_sim(n_datasets, n_vars, N, eigen_type, maxit = 10000L)
}
\arguments{
\item{n_datasets}{numeric. Number of datasets with dimensions (N, n_vars) to simulate.}

\item{n_vars}{numeric. Number of variables / indicators in dataset.}

\item{N}{numeric. Number of cases / observations in dataset.}

\item{eigen_type}{numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of correlation matrix at 1) or PAF (eigen_type = 2; i.e., setting diagonal of correlation matrix to SMCs).}

\item{maxit}{numeric. Maximum iterations to perform after which to abort.}
}
\description{
Function called from within PARALLEL so usually no call to this is needed by the user.
Provides a C++ implementation of the PARALLEL simulation procedure
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.COMPARE.R
\name{print.COMPARE}
\alias{print.COMPARE}
\title{Print COMPARE object}
\usage{
\method{print}{COMPARE}(x, ...)
}
\arguments{
\item{x}{list. An object of class COMPARE to be printed}

\item{...}{Further arguments for print.}
}
\description{
Print Method showing a summarized output of the \code{\link{COMPARE}} function.
}
\examples{
# A type SPSS EFA to mimick the SPSS implementation
EFA_SPSS_5 <- EFA(IDS2_R, n_factors = 5, type = "SPSS")

# A type psych EFA to mimick the psych::fa() implementation
EFA_psych_5 <- EFA(IDS2_R, n_factors = 5, type = "psych")

# compare the two
COMPARE(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
        x_labels = c("SPSS", "psych"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EFA_AVERAGE.R
\name{EFA_AVERAGE}
\alias{EFA_AVERAGE}
\title{Model averaging across different EFA methods and types}
\source{
Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
 in R and SPSS. Manuscript in Preparation.

Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
rotation to oblique simple structure. British Journal of Statistical Psychology,
17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The
Hull Method for Selecting the Number of Common Factors, Multivariate Behavioral
Research, 46, 340-364, doi: 10.1080/00273171.2011.564527

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233
}
\usage{
EFA_AVERAGE(
  x,
  n_factors,
  N = NA,
  method = "PAF",
  rotation = "promax",
  type = "none",
  averaging = c("mean", "median"),
  trim = 0,
  salience_threshold = 0.3,
  max_iter = 10000,
  init_comm = c("smc", "mac", "unity"),
  criterion = c(0.001),
  criterion_type = c("sum", "max_individual"),
  abs_eigen = c(TRUE),
  varimax_type = c("svd", "kaiser"),
  normalize = TRUE,
  k_promax = 2:4,
  k_simplimax = ncol(x),
  P_type = c("norm", "unnorm"),
  precision = 1e-05,
  start_method = c("psych", "factanal"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall"),
  show_progress = TRUE
)
}
\arguments{
\item{x}{data.frame or matrix. Dataframe or matrix of raw data or matrix with
correlations. If raw data is entered, the correlation matrix is found from the
data.}

\item{n_factors}{numeric. Number of factors to extract.}

\item{N}{numeric. The number of observations. Needs only be specified if a
correlation matrix is used. If input is a correlation matrix and \code{N} = NA
(default), not all fit indices can be computed.}

\item{method}{character vector. Any combination of  "PAF", "ML", and "ULS",
to use principal axis factoring, maximum likelihood, or unweighted least
squares (also called minres), respectively, to fit the EFAs. Default is "PAF".}

\item{rotation}{character vector. Either perform no rotation ("none"),
any combination of orthogonal rotations ("varimax", "equamax", "quartimax", "geominT",
"bentlerT", and "bifactorT"; using "orthogonal" runs all of these), or of
oblique rotations ("promax", "oblimin", "quartimin", "simplimax", "bentlerQ",
"geominQ", and "bifactorQ"; using "oblique" runs all of these). Rotation types
(no rotation, orthogonal rotations, and oblique rotations) cannot be mixed.
Default is "promax".}

\item{type}{character vector. Any combination of "none" (default), "EFAtools",
"psych", and "SPSS" can be entered. "none" allows the specification of various
combinations of the arguments controlling both factor extraction methods and
the rotations. The others ("EFAtools", "psych", and "SPSS"), control the execution
of the respective factor extraction method and rotation to be in line with how
it is executed in this package (i.e., the respective default procedure), in the
psych package, and in SPSS. A specific psych implementation exists for PAF, ML, varimax,
and promax. The SPSS implementation exists for PAF, varimax, and promax. For
details, see \code{\link{EFA}}.}

\item{averaging}{character. One of "mean" (default), and "median". Controls
whether the different results should be averaged using the (trimmed) mean,
or the median.}

\item{trim}{numeric. If averaging is set to "mean", this argument controls
the trimming of extremes (for details see \code{\link[base:mean]{base::mean}}).
By default no trimming is done (i.e., trim = 0).}

\item{salience_threshold}{numeric. The threshold to use to classify a pattern
coefficient or loading as salient (i.e., substantial enough to assign it to
a factor). Default is 0.3. Indicator-to-factor correspondences will be inferred
based on this threshold. Note that this may not be meaningful if rotation = "none"
and n_factors > 1 are used, as no simple structure is present there.}

\item{max_iter}{numeric. The maximum number of iterations to perform after which
the iterative PAF procedure is halted with a warning. Default is 10,000. Note
that non-converged procedures are excluded from the averaging procedure.}

\item{init_comm}{character vector. Any combination of "smc", "mac", and "unity".
Controls the methods to estimate the initial communalities in \code{PAF} if
"none" is among the specified types. "smc" will use squared multiple
correlations, "mac" will use maximum absolute correlations, "unity" will use
1s (for details see \code{\link{EFA}}). Default is \code{c("smc", "mac", "unity")}.}

\item{criterion}{numeric vector. The convergence criterion used for PAF if
"none" is among the specified types.
If the change in communalities from one iteration to the next is smaller than
this criterion the solution is accepted and the procedure ends.
Default is \code{0.001}.}

\item{criterion_type}{character vector. Any combination of "max_individual" and
"sum". Type of convergence criterion used for PAF if "none" is among the
specified types. "max_individual" selects the maximum change in any of the
communalities from one iteration to the next and tests it against the
specified criterion. "sum" takes the difference of
the sum of all communalities in one iteration and the sum of all communalities
in the next iteration and tests this against the criterion
(for details see \code{\link{EFA}}). Default is \code{c("sum", "max_individual")}.}

\item{abs_eigen}{logical vector. Any combination of TRUE and FALSE.
Which algorithm to use in the PAF iterations if "none" is among the specified
types. If FALSE, the loadings are computed from the eigenvalues. This is also
used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
loadings are computed with the absolute eigenvalues as done by SPSS
(for details see \code{\link{EFA}}). Default is \code{TRUE}.}

\item{varimax_type}{character vector. Any combination of "svd" and "kaiser".
The type of the varimax rotation performed if "none" is among the specified
types and "varimax", "promax", "orthogonal", or "oblique" is among the specified
rotations. "svd" uses singular value decomposition, as
\link[stats:varimax]{stats::varimax} does, and "kaiser" uses the varimax
procedure performed in SPSS. This is the original procedure from Kaiser (1958),
but with slight alterations in the varimax criterion (for details, see
\code{\link{EFA}} and Grieder & Steiner, 2020).
Default is \code{c("svd", "kaiser")}.}

\item{normalize}{logical vector. Any combination of TRUE and FALSE.
\code{TRUE} performs a kaiser normalization before the specified rotation(s).
Default is \code{TRUE}.}

\item{k_promax}{numeric vector. The power used for computing the target matrix
P in the promax rotation if "none" is among the specified types and "promax"
or "oblique" is among the specified rotations. Default is \code{2:4}.}

\item{k_simplimax}{numeric. The number of 'close to zero loadings' for the
simplimax rotation (see \code{\link[GPArotation:GPA]{GPArotation::GPFoblq}})
if "simplimax" or "oblique" is among the specified rotations. Default
is \code{ncol(x)}, where x is the entered data.}

\item{P_type}{character vector. Any combination of "norm" and "unnorm".
This specifies how the target matrix P is computed in promax rotation if
"none" is among the specified types and "promax" or "oblique" is among the
specified rotations. "unnorm" will use the unnormalized target matrix as
originally done in Hendrickson and White (1964). "norm" will use a
normalized target matrix (for details see \code{\link{EFA}}).
Default is \code{c("norm", "unnorm")}.}

\item{precision}{numeric vector. The tolerance for stopping in the rotation
procedure(s). Default is 10^-5.}

\item{start_method}{character vector. Any combination of "psych" and "factanal".
How to specify the starting values for the optimization procedure for ML.
"psych" takes the starting values specified in \link[psych:fa]{psych::fa}.
"factanal" takes the starting values specified in the
\link[stats:factanal]{stats::factanal} function. Default is
\code{c("psych", "factanal")}.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
is given as input. Default is "pairwise.complete.obs".}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is "pearson".}

\item{show_progress}{logical. Whether a progress bar should be shown in the
console. Default is TRUE.}
}
\value{
A list of class EFA_AVERAGE containing
\item{orig_R}{Original correlation matrix.}
\item{h2}{A list with the average, standard deviation, minimum, maximum, and
range of the final communality estimates across the factor solutions.}
\item{loadings}{A list with the average, standard deviation, minimum, maximum,
and range of the final loadings across the factor solutions. If rotation was
"none", the unrotated loadings, otherwise the rotated loadings (pattern
coefficients).}
\item{Phi}{A list with the average, standard deviation, minimum, maximum, and
range of the factor intercorrelations across factor solutions obtained with
oblique rotations.}
\item{ind_fac_corres}{A matrix with each cell containing the proportion of
the factor solutions in which the respective indicator-to-factor correspondence
occurred, i.e., in which the loading exceeded the specified salience threshold.
Note: Rowsums can exceed 1 due to cross-loadings.}
\item{vars_accounted}{A list with the average, standard deviation, minimum,
maximum, and range of explained variances and sums of squared loadings across
the factor solutions. Based on the unrotated loadings.}
\item{fit_indices}{A matrix containing the average, standard deviation,
minimum, maximum, and range for all applicable fit indices across the respective
factor solutions, and the degrees of freedom (df). If the method argument
contains ML or ULS: Fit indices derived
from the unrotated factor loadings: Chi Square (chisq), including significance
level, Comparative Fit Index (CFI), Root Mean Square Error of Approximation
(RMSEA), Akaike Information Criterion (AIC), Bayesian Information Criterion
(BIC)and the common part accounted for (CAF) index as proposed by
Lorenzo-Seva, Timmerman, & Kiers (2011). For PAF, only the CAF can be
calculated (see details).}
\item{implementations_grid}{A matrix containing, for each performed EFA,
the setting combination, if an error occurred (logical), the error message
(character), an integer code for convergence as returned by
\code{\link[stats:optim]{stats:optim}} (0 indicates successful completion.),
if heywood cases occurred (logical, see details for definition), if the
solution was admissible (logical, see details for definition), and the fit
indices.}
\item{efa_list}{A list containing the outputs of all performed EFAs. The names
correspond to the rownames from the implementations_grid.}
\item{settings}{A list of the settings used.}
}
\description{
Not all EFA procedures always arrive at the same solution. This function allows
you perform a number of EFAs from different methods (e.g., Maximum Likelihood
and Principal Axis Factoring), with different implementations (e.g., the SPSS
and psych implementations of Principal Axis Factoring), and across different
rotations of the same type (e.g., multiple oblique rotations, like promax and
oblimin). EFA_AVERAGE will then run all these EFAs (using the \code{\link{EFA}}
function) and provide a summary across the different solutions.
}
\details{
As a first step in this function, a grid is produced containing the setting
combinations for the to-be-performed EFAs. These settings are then entered as
arguments to the \code{\link{EFA}} function and the EFAs are run in a second
step. After all EFAs are run, the factor solutions are averaged and their
variability determined in a third step.

The grid containing the setting combinations is produced based on the entries
to the respective arguments. To this end, all possible combinations resulting
in unique EFA models are considered. That is, if, for example, the \code{type}
argument was set to \code{c("none", "SPSS")} and one combination of the specific
settings entered was identical to the SPSS combination, this combination
would be included in the grid and run only once. We include here a list
of arguments that are only evaluated under specific conditions:

The arguments \code{init_comm}, \code{criterion}, \code{criterion_type},
\code{abs_eigen} are only evaluated if "PAF" is included in \code{method}
and "none" is included in \code{type}.

The argument \code{varimax_type} is only evaluated if "varimax", "promax",
"oblique", or "orthogonal" is included in \code{rotation} and "none" is
included in \code{type}.

The argument \code{normalize} is only evaluated if \code{rotation} is not
set to "none" and "none" is included in \code{type}.

The argument \code{k_simplimax} is only evaluated if "simplimax" or "oblique"
is included in \code{rotation}.

The arguments \code{k_promax} and \code{P_type} are only evaluated if
"promax" or "oblique" is included in \code{rotation} and "none" is included
in \code{type}.

The argument \code{start_method} is only evaluated if "ML" is included in
\code{method}.

To avoid a bias in the averaged factor solutions from problematic solutions,
these are excluded prior to averaging. A solution is deemed problematic if
at least one of the following is true: an error occurred, the model did not
converge, or there is at least one Heywood case (defined as a loading or communality of >= .998).
Information on errors, convergence, and Heywood cases are returned in the
implementations_grid and a summary of these is given when printing the output.
In addition to these, information on the admissibility of the factor solutions
is also included. A solution was deemed admissible if (1) no error occurred,
(2) the model converged, (3) no Heywood cases are present, and (4) there are
at least two salient loadings (i.e., loadings exceeding the specified
\code{salience_threshold}) for each factor. So, solutions failing one of the
first three of these criteria of admissibility are also deemed problematic and
therefore excluded from averaging. However, solutions failing only
the fourth criterion of admissibility are still included for averaging.
Finally, if all solutions are problematic (e.g., all solutions contain
Heywood cases), no averaging is performed and the respective outputs are NA.
In this case, the implementations_grid should be inspected to see if there
are any error messages, and the separate EFA solutions that are also included
in the output can be inspected as well, for example, to see where Heywood
cases occurred.

A core output of this function includes the average, minimum, and maximum
loadings derived from all non-problematic (see above) factor solutions. Please
note that these are not entire solutions, but the matrices include the average,
minimum, or maximum value for each cell (i.e., each loading separately). This
means that, for example, the matrix with the minimum loadings will contain
the minimum value in any of the factor solutions for each specific loading,
and therefore most likely contains loadings from different factor solutions.
The matrices containing the minimum and maximum factor solutions can
therefore not be interpreted as whole factor solutions.

The output also includes information on the average, minimum, maximum, and
variability of the fit indices across the non-problematic factor solutions.
It is important to note that not all fit indices are computed for all fit
methods: For ML and ULS, all fit indices can be computed, while for PAF, only
the common part accounted for (CAF) index (Lorenzo-Seva, Timmerman, & Kiers, 2011)
can be computed. As a consequence, if only "PAF" is included in the
\code{method} argument, averaging can only be performed for the CAF, and the
other fit indices are NA. If a combination of "PAF" and "ML" and/or "ULS" are
included in the \code{method} argument, the CAF is averaged across all non-
problematic factor solutions, while all other fit indices are only averaged
across the ML and ULS solutions. The user should therefore keep in mind that
the number of EFAs across which the fit indices are averaged can diverge for
the CAF compared to all other fit indices.
}
\examples{
# Averaging across different implementations of PAF and promax rotation (72 EFAs)
Aver_PAF <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)

# Use median instead of mean for averaging (72 EFAs)
Aver_PAF_md <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                           averaging = "median")

# Averaging across different implementations of PAF and promax rotation,
# and across ULS and different versions of ML (108 EFAs)
Aver_meth_ext <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                             method = c("PAF", "ULS", "ML"))

# Averaging across one implementation each of PAF (EFAtools type), ULS, and
# ML with one implementation of promax (EFAtools type) (3 EFAs)
Aver_meth <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                         method = c("PAF", "ULS", "ML"), type = "EFAtools",
                         start_method = "psych")

# Averaging across different oblique rotation methods, using one implementation
# of ML and one implementation of promax (EFAtools type) (7 EFAs)
Aver_rot <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                         method = "ML", rotation = "oblique", type = "EFAtools",
                         start_method = "psych")


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EKC.R
\name{EKC}
\alias{EKC}
\title{Empirical Kaiser Criterion}
\source{
Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of extraction
methods under realistic conditions. Psychological Methods, 24(4), 468–491.
https://doi.org/10.1037/met0000200

Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
Psychological Methods, 22, 450 – 466. http://dx.doi.org/10.1037/ met0000074

Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain. Psychological Bulletin, 99,
432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
}
\usage{
EKC(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall")
)
}
\arguments{
\item{x}{data.frame or matrix. data.frame or matrix of raw data or matrix with
correlations.}

\item{N}{numeric. The number of observations. Only needed if x is a correlation
matrix.}

\item{use}{character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
data is given as input. Default is  \code{"pairwise.complete.obs"}.}

\item{cor_method}{character. Passed to \code{\link[stats:cor]{stats::cor}}.
Default is  \code{"pearson"}.}
}
\value{
A list of class EKC containing

\item{eigenvalues}{A vector containing the eigenvalues found on the correlation matrix of the entered data.}
\item{n_factors}{The number of factors to retain according to the empirical Kaiser criterion.}
\item{references}{The reference eigenvalues.}
\item{settings}{A list with the settings used.}
}
\description{
The empirical Kaiser criterion incorporates random sampling variations of the
eigenvalues from the Kaiser-Guttman criterion (\code{\link{KGC}}; see Auerswald & Moshagen
, 2019; Braeken & van Assen, 2017). The code is based on Auerswald and Moshagen
(2019).
}
\details{
The Kaiser-Guttman criterion was defined with the intend that a factor
 should only be extracted if it explains at least as much variance as a single
 factor (see \code{\link{KGC}}). However, this only applies to population-level
 correlation matrices. Due to sampling variation, the KGC strongly overestimates
 the number of factors to retrieve (e.g., Zwick & Velicer, 1986). To account
 for this and to introduce a factor retention method that performs well with
 small number of indicators and correlated factors (cases where the performance
 of parallel analysis, see \code{\link{PARALLEL}}, is known to deteriorate)
 Braeken and van Assen (2017) introduced the empirical Kaiser criterion in
 which a series of reference eigenvalues is created as a function of the
 variables-to-sample-size ratio and the observed eigenvalues.

 Braeken and van Assen (2017) showed that "(a) EKC performs about as well as
 parallel analysis for data arising from the null, 1-factor, or orthogonal
 factors model; and (b) clearly outperforms parallel analysis for the specific
 case of oblique factors, particularly whenever factor intercorrelation is
 moderate to high and the number of variables per factor is small, which is
 characteristic of many applications these days" (p.463-464).

 The \code{EKC} function can also be called together with other factor
  retention criteria in the \code{\link{N_FACTORS}} function.
}
\examples{
EKC(test_models$baseline$cormat, N = 500)
}
\seealso{
Other factor retention criteria: \code{\link{CD}},
 \code{\link{HULL}}, \code{\link{KGC}}, \code{\link{PARALLEL}},
 \code{\link{SMT}}

  \code{\link{N_FACTORS}} as a wrapper function for this and all
  the above-mentioned factor retention criteria.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.HULL.R
\name{plot.HULL}
\alias{plot.HULL}
\title{Plot HULL object}
\usage{
\method{plot}{HULL}(x, ...)
}
\arguments{
\item{x}{list of class HULL. An output from the \code{\link{HULL}} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \code{\link{HULL}} function
}
\examples{
\donttest{
x <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
plot(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.SMT.R
\name{print.SMT}
\alias{print.SMT}
\title{Print SMT object}
\usage{
\method{print}{SMT}(x, ...)
}
\arguments{
\item{x}{list of class SMT (output from the \link{SMT} function)}

\item{...}{additional arguments passed to print}
}
\description{
Print SMT object
}
\examples{
SMT_base <- SMT(test_models$baseline$cormat, N = 500)
SMT_base

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.EFA.R
\name{print.EFA}
\alias{print.EFA}
\title{Print EFA object}
\usage{
\method{print}{EFA}(x, ...)
}
\arguments{
\item{x}{list. An object of class EFA to be printed}

\item{...}{Further arguments for print.}
}
\description{
Print Method showing a summarized output of the \link{EFA} function
}
\examples{
EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                    type = "EFAtools", method = "PAF", rotation = "promax")
EFAtools_PAF

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.KGC.R
\name{plot.KGC}
\alias{plot.KGC}
\title{Plot KGC object}
\usage{
\method{plot}{KGC}(x, ...)
}
\arguments{
\item{x}{a list of class KGC. An output from the \link{KGC} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \link{KGC} function
}
\examples{
KGC_base <- KGC(test_models$baseline$cormat)
plot(KGC_base)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WJIV_ages_20_39_doc.R
\docType{data}
\name{WJIV_ages_20_39}
\alias{WJIV_ages_20_39}
\title{Woodcock Johnson IV: ages 20 to 39}
\format{
A list of 2 with elements "cormat" (47 x 47 matrix of bivariate correlations)
and "N" (scalar). The correlation matrix contains the following variables:
\describe{
  \item{ORLVOC}{(numeric) - Oral Vocabulary.}
  \item{NUMSER}{(numeric) - Number Series.}
  \item{VRBATN}{(numeric) - Verbal Attention.}
  \item{LETPAT}{(numeric) - Letter-Pattern Matching.}
  \item{PHNPRO}{(numeric) - Phonological Processing.}
  \item{STYREC}{(numeric) - Story Recall.}
  \item{VISUAL}{(numeric) - Visualization.}
  \item{GENINF}{(numeric) - General Information.}
  \item{CONFRM}{(numeric) - Concept Formation.}
  \item{NUMREV}{(numeric) - Numbers Reversed.}
  \item{NUMPAT}{(numeric) - Number-Pattern Matching.}
  \item{NWDREP}{(numeric) - Nonword Repetition.}
  \item{VAL}{(numeric) - Visual-Auditory Learning.}
  \item{PICREC}{(numeric) - Picture Recognition.}
  \item{ANLSYN}{(numeric) - Analysis-Synthesis.}
  \item{OBJNUM}{(numeric) - Object-Number Sequencing.}
  \item{PAIRCN}{(numeric) - Pair Cancellation.}
  \item{MEMWRD}{(numeric) - Memory for Words.}
  \item{PICVOC}{(numeric) - Picture Vocabulary.}
  \item{ORLCMP}{(numeric) - Oral Comprehension.}
  \item{SEGMNT}{(numeric) - Segmentation.}
  \item{RPCNAM}{(numeric) - Rapid Picture Naming.}
  \item{SENREP}{(numeric) - Sentence Repetition.}
  \item{UNDDIR}{(numeric) - Understanding Directions.}
  \item{SNDBLN}{(numeric) - Sound Blending.}
  \item{RETFLU}{(numeric) - Retrieval Fluency.}
  \item{SNDAWR}{(numeric) - Sound Awareness.}
  \item{LWIDNT}{(numeric) - Letter-Word Identification.}
  \item{APPROB}{(numeric) - Applied Problems.}
  \item{SPELL}{(numeric) - Spelling.}
  \item{PSGCMP}{(numeric) - Passage Comprehension.}
  \item{CALC}{(numeric) - Calculation.}
  \item{WRTSMP}{(numeric) - Writing Samples.}
  \item{WRDATK}{(numeric) - Word Attack.}
  \item{ORLRDG}{(numeric) - Oral Reading.}
  \item{SNRDFL}{(numeric) - Sentence Reading Fluency.}
  \item{MTHFLU}{(numeric) - Math Facts Fluency.}
  \item{SNWRFL}{(numeric) - Sentence Writing Fluency.}
  \item{RDGREC}{(numeric) - Reading Recall.}
  \item{NUMMAT}{(numeric) - Number Matrices.}
  \item{EDIT}{(numeric) - Editing.}
  \item{WRDFLU}{(numeric) - Word Reading Fluency.}
  \item{SPLSND}{(numeric) - Spelling of Sounds.}
  \item{RDGVOC}{(numeric) - Reading Vocabulary.}
  \item{SCI}{(numeric) - Science.}
  \item{SOC}{(numeric) - Social Studies.}
  \item{HUM}{(numeric) - Humanities.}
 }
}
\source{
McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
 Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.

Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
Rolling Meadows, IL: Riverside.
}
\usage{
WJIV_ages_20_39
}
\description{
A list containing the bivariate correlations (N = 1,251) of the 47
cognitive and achievement subtests from the WJ IV for the 20- to 39-year-olds
from the standardization sample obtained from
the WJ-IV technical manual (McGrew, LaForte, & Schrank, 2014).
Tables are reproduced with permission from the publisher.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.EKC.R
\name{plot.EKC}
\alias{plot.EKC}
\title{Plot EKC object}
\usage{
\method{plot}{EKC}(x, ...)
}
\arguments{
\item{x}{a list of class EKC. An output from the \link{EKC} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \link{EKC} function
}
\examples{
EKC_base <- EKC(test_models$baseline$cormat, N = 500)
plot(EKC_base)

}
