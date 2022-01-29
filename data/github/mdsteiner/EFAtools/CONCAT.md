
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
