---
affiliations:
- index: 1
  name: 'Comtek Advanced Structures, Ltd.'
authors:
- affiliation: 1
  name: Stefan Kloppenborg
  orcid: '0000-0002-1908-5214'
bibliography: 'paper.bib'
date: '6/24/2020'
output:
  md_document:
    pandoc_args: '--atx-headers'
    preserve_yaml: True
    variant: markdown
  pdf_document: default
tags:
- R
- statistics
- composite materials
- material science
title: |
    cmstatr: An R Package for Statistical Analysis of Composite Material
    Data
---

# Summary

Strength data for composite materials used in aerospace applications,
such as carbon fiber and fiberglass reinforced composites, are normally
analyzed using statistical methods because of the inherent variability
in the constituent materials and in the processing. The design standards
for civil aviation requires that the probability of structural failure
due to this variability must be minimized, and to do so, the designer
must use what are called "Design Values" for each material in stress
analyses and ensure they exceed the actual stresses experienced by those
materials in service. Design Values are set based on the one-sided lower
confidence bound of the material strength. For some types of structure,
the content of this lower confidence bound is $99\%$ with a confidence
level of $95\%$; in this case, the confidence bound is referred to as
A-Basis. For some other types of structure, the content of the lower
confidence bound is instead $90\%$ with a confidence level of $95\%$; in
this case, the confidence bound is referred to as B-Basis. The
statistical methods for calculating these basis values are outlined in
Composite Materials Handbook, Volume 1, Revision G, or CMH-17-1G in
short [@CMH171G]. The use of these methods is widely accepted by
industry and civil aviation regulators.

Design Values are often adjusted to account for anticipated in-service
damage and other factors, however those adjustments are outside the
scope of the present software package.

For a detailed discussion of the theory and applications of tolerance
bounds, the reader is referred to @Meeker_Hahn_Escobar_2017 or
@Krishnamoorthy_Mathew_2008.

Currently, many users use MS Excel spreadsheets to perform these
analyses. The MS Excel spreadsheets typically used, such as `STAT-17`
[@STAT-17], `ASAP` [@Raju_Tomblin_2008] and `CMH17-STATS`
[@CMH17-STATS], use password-protected `VBA` macros to perform the
computations. As such, the code cannot be audited by the user. `cmstatr`
is an R package that addresses this issue by implementing the same
statistical analysis techniques found in CMH-17-1G in an open-source
environment.

# Statement of Need

The purpose of `cmstatr` is to:

-   Provide a consistent user interface for computing A- and B-Basis
    values and performing the related diagnostic tests in the `R`
    programming environment
-   Allow auditing of the code used to compute A- and B-Basis values,
    and performing the related diagnostic tests
-   Enable users to automate computation workflows or to perform
    simulation studies

# Implementation Goals

`cmstatr` aims give a consistent interface for the user. Most functions
are written to work with the `tidyverse` [@tidyverse] and most functions
have similar argument lists. The intent is to make the package easy to
learn and use.

The implementation of `cmstatr` also aims to avoid the use of lookup
tables that are prevalent in calculation spreadsheets and minimize the
use of approximations. While this decision leads to increased
computation time, the typically small data sets (tens to hundreds of
observations) associated with composite material test data and the speed
of modern computers make this practical for interactive programming.

# Example Usage

Normally, to use `cmstatr` the user will load `cmstatr` itself as well
as the `tidyverse` package.

``` {.r}
library(cmstatr)
library(tidyverse)
```

`cmstatr` contains some example data sets, which can be used to
demonstrate the features of the package. One of those data sets ---
`carbon.fabric.2` --- will be used in the following example. This data
set contains results from several mechanical tests of a typical
composite material, and contains the typical measurements obtained from
a test lab. In the following examples, results from tension testing in
the warp fiber direction (`WT`) will be used. Part of this data set is
shown below.

``` {.r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  head(10)
```

    ##    test condition batch thickness nplies strength modulus failure_mode
    ## 1    WT       CTD     A     0.112     14  142.817   9.285          LAT
    ## 2    WT       CTD     A     0.113     14  135.901   9.133          LAT
    ## 3    WT       CTD     A     0.113     14  132.511   9.253          LAT
    ## 4    WT       CTD     A     0.112     14  135.586   9.150          LAB
    ## 5    WT       CTD     A     0.113     14  125.145   9.270          LAB
    ## 6    WT       CTD     A     0.113     14  135.203   9.189          LGM
    ## 7    WT       CTD     A     0.113     14  128.547   9.088          LAB
    ## 8    WT       CTD     B     0.113     14  127.709   9.199          LGM
    ## 9    WT       CTD     B     0.113     14  127.074   9.058          LGM
    ## 10   WT       CTD     B     0.114     14  126.879   9.306          LGM

One common task is to calculate B-Basis values. There are several
statistical methods for doing so, depending on the distribution of the
data. The single-point basis functions automatically perform the
following diagnostic tests:

-   the maximum normed residual test for outliers within a batch
    [@CMH171G],
-   the Anderson--Darling k-Sample test to check if batches are drawn
    from the same (unspecified) distribution [@Scholz_Stephens_1987],
-   the maximum normed residual test for outliers within the data, and
-   the Anderson--Darling test for a particular probability distribution
    [@Lawless_1982].

Assuming that the data from the warp tension (WT) tested at
elevated-temperature/wet condition (ETW) follows a normal distribution,
then this can be done using the function. Note that all of the functions
in `cmstatr` that compute basis values default to computing tolerance
bounds with a content of $p=0.9$ and a confidence of $conf=0.95$, or
B-Basis.

``` {.r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  filter(condition == "ETW") %>%
  basis_normal(strength, batch)
```

    ## Warning: `anderson_darling_normal` failed: Anderson-Darling test rejects
    ## hypothesis that data is drawn from a normal distribution

    ## 
    ## Call:
    ## basis_normal(data = ., x = strength, batch = batch)
    ## 
    ## Distribution:  Normal    ( n = 18 )
    ## The following diagnostic tests failed: 
    ##     `anderson_darling_normal`
    ## B-Basis:   ( p = 0.9 , conf = 0.95 )
    ## 122.9315

All of the various basis functions perform diagnostic tests for each of
the statistical tests mentioned above. If any of the diagnostic tests
failed, a warning is shown and the test failure is also recorded in the
returned object (and shown in that object's `print` method). In the
example above, the output shows that the Anderson--Darling test for
normality [@Lawless_1982] rejects the hypothesis that the data is drawn
from a normal distribution.

Two non-parametric basis calculations, based on @Guenther_1970 and
@Vangel_1994 are also implemented in `cmstatr`. These functions perform
the same diagnostic tests, but omit the Anderson--Darling test for a
particular distribution.

The diagnostic test can be run directly using `cmstatr` as well. For
example, the failed diagnostic test above can be run as follows:

``` {.r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  filter(condition == "ETW") %>%
  anderson_darling_normal(strength)
```

    ## 
    ## Call:
    ## anderson_darling_normal(data = ., x = strength)
    ## 
    ## Distribution:  Normal ( n = 18 ) 
    ## Test statistic:  A = 0.9381665 
    ## OSL (p-value):  0.01103075  (assuming unknown parameters)
    ## Conclusion: Sample is not drawn from a Normal distribution ( alpha = 0.05 )

If the failure of a diagnostic test is decided to be acceptable, the
test result can be overridden to hide the warning in the basis function
output:

``` {.r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  filter(condition == "ETW") %>%
  basis_normal(strength, batch, override = c("anderson_darling_normal"))
```

    ## 
    ## Call:
    ## basis_normal(data = ., x = strength, batch = batch, override = c("anderson_darling_normal"))
    ## 
    ## Distribution:  Normal    ( n = 18 )
    ## The following diagnostic tests were overridden: 
    ##     `anderson_darling_normal`
    ## B-Basis:   ( p = 0.9 , conf = 0.95 )
    ## 122.9315

`cmstatr` also provides functions for calculating basis values from data
pooled across multiple testing environments, as recommended by
[@CMH171G]. There are two methods of pooling: both calculate a measure
of global variation from the data in all tested environmental
conditions, then they calculate a basis value for each condition using
the measure of the global variation and the individual mean values of
each condition. The two methods of pooling have different underlying
assumptions and hence the data must pass a different set diagnostic
tests for each of the two functions.

# Validation and Comparison With Existing Tools

Where possible, `cmstatr` has been verified against the examples given
in the articles in which the statistical methods were published. Unit
tests have been written so that this verification is re-checked
routinely to prevent unintended regressions. Agreement between `cmstatr`
and the examples in the original articles is within the expected numeric
accuracy.

`cmstatr` has also been verified against existing software, such as
`STAT-17` [@STAT-17], `ASAP` [@Raju_Tomblin_2008] and `CMH17-STATS`
[@CMH17-STATS] using several example data sets. Agreement between
`cmstatr` and the other software is generally good, but some results
differ slightly, likely due to various approximations used in the
software. Comparison between `cmstatr` and the other software is
performed within various unit tests to guard against future regressions.

The tests are automatically run each time a change is made to the code
of `cmstatr` using a continuous integration service.

# Reproducibility

It is envisioned that many users of `cmstatr` will use it within an R
Notebook or a Jupyter Notebook. It is further envisioned that this
notebook will be directly converted into the statistical analysis
report. If this is done, the reader of the statistical report will be
able to verify all of the detailed steps used in the statistical
analysis.

# Acknowledgement

The author would like to thank Mr. Billy Cheng for his contributions to
`cmstatr` and this paper. The author would also like to thank Comtek
Advanced Structures Ltd. for its support in developing and releasing the
`cmstatr` package.

# References {#references .unnumbered}
# Logo
The file `logo.svg` can be edited with inkscape.

# Export Sizes
## /man/figures/logo.png
Hide "White Background" layer. Hide "OpenGraphBorder."
Export from inkscape with size 240x278

## /man/figures/logo-wbg-240x278.png
Show "White Background" layer. Hide "OpenGraphBorder."
Export from inkscape with size 240x278

## /man/figures/logo-wbg-1280x640
Show "OpenGraphBorder". Export from inkscape with size 1280x640

# Favico
`favico.ico` can be created by running:

```r
pkgdown::build_favicons()
```

<!-- README.md is generated from README.Rmd. Please edit that file -->

# cmstatr <img src="man/figures/logo.png" align="right" alt="" width="120" />

<!-- badges: start -->

[![R build
status](https://github.com/cmstatr/cmstatr/workflows/R-CMD-check/badge.svg)](https://github.com/cmstatr/cmstatr/actions?workflow=R-CMD-check)
[![`Codecov` test
coverage](https://codecov.io/gh/cmstatr/cmstatr/branch/master/graph/badge.svg)](https://codecov.io/gh/cmstatr/cmstatr?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02265/status.svg)](https://doi.org/10.21105/joss.02265)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/cmstatr)](https://cran.r-project.org/package=cmstatr)
[![](https://cranlogs.r-pkg.org/badges/cmstatr)](https://cran.r-project.org/package=cmstatr)
<!-- badges: end -->

# What It Does

The `cmstatr` package provides functions for performing statistical
analysis of composite material data. The statistical methods implemented
are those described in [CMH-17-1G](https://www.cmh17.org/). This package
focuses on calculating basis values (lower tolerance bounds) for
material strength properties, as well as performing the associated
diagnostic tests. Functions are also provided for testing for
equivalency between alternate samples and the “qualification” or
“baseline” samples.

Additional details about the package are available in the paper by
Kloppenborg (2020, <https://doi.org/10.21105/joss.02265>).

# Installation

To install `cmstatr` from CRAN, simply run:

``` r
install.packages("cmstatr")
```

If you want the latest development version, you can install it from
`github` using `devtools`. This will also install the dependencies
required to build the vignettes. Optionally, change the value of the
argument `ref` to install `cmstatr` from a different branch of the
repository.

``` r
install.packages(c("devtools", "rmarkdown", "dplyr", "tidyr"))
devtools::install_github("cmstatr/cmstatr", build_vignettes = TRUE,
                         ref = "master",
                         build_opts = c("--no-resave-data", "--no-manual"))
```

# Usage

To compute a B-Basis value from an example data set packaged with
`cmstatr` you can do the following:

``` r
library(dplyr)
library(cmstatr)

carbon.fabric.2 %>%
  filter(test == "FC") %>%
  filter(condition == "RTD") %>%
  basis_normal(strength, batch)
#> 
#> Call:
#> basis_normal(data = ., x = strength, batch = batch)
#> 
#> Distribution:  Normal    ( n = 18 )
#> B-Basis:   ( p = 0.9 , conf = 0.95 )
#> 76.88082
```

For more examples of usage of the `cmstatr` package, see the tutorial
vignette, which can be [viewed
online](https://www.cmstatr.net/articles/cmstatr_Tutorial.html), or can
be loaded as follows, once the package is installed:

``` r
vignette("cmstatr_Tutorial")
```

There is also a vignette showing some examples of the types of graphs
that are typically produced when analyzing composite materials. You can
view this [vignette
online](https://www.cmstatr.net/articles/cmstatr_Graphing.html), or you
can load this vignette with:

``` r
vignette("cmstatr_Graphing")
```

# Philosophical Notes

This package expects
[`tidy data`](https://www.jstatsoft.org/article/view/v059i10). That is,
individual observations should be in rows and variables in columns.

Where possible, this package uses general solutions. Look-up tables are
avoided wherever possible.

# Issues

If you’ve found a bug, please open an issue in this repository and
describe the bug. Please include a [reproducible
example](https://reprex.tidyverse.org/) of the bug. If you’re able to
fix the bug, you can do so by submitting a pull request.

If your bug is related to a particular data set, sharing that data set
will help to fix the bug. If you cannot share the data set, please strip
any identifying information and optionally scale the data by an
unspecified factor so that the bug can be reproduced and diagnosed.

# Contributing

Contributions to `cmstatr` are always welcomed. For small changes
(fixing typos or improving the documentation), go ahead and submit a
pull request. For more significant changes, such as new features, please
discuss the proposed change in an issue first.

## Contribution Guidelines

-   Please create a git branch for each pull request (PR)
-   Before submitting a pull request, please make sure that
    `R CMD check` passes with no errors, warnings or notes
-   New and modified code should follow the style guide enforced by the
    [`lintr`](https://cran.r-project.org/package=lintr) package
-   Document all exported functions using
    [`roxygen2`](https://cran.r-project.org/package=roxygen2)
-   Write tests using
    [`testthat`](https://cran.r-project.org/package=testthat). If your
    contribution fixes a bug, then the test(s) that you add should fail
    before your bug-fix patch is applied and should pass after the code
    is patched.
-   For changes that affect the user, add a bullet at the top of
    `NEWS.md` below the current development version

## Development

Testing is performed using `testthat`. Edition 3 of that package is used
and parallel processing enabled. If you wish to use more than two CPUs,
set the environment variable `TESTTHAT_CPUS` to the number of CPUs that
you want to use. One way of doing this is to create the file `.Rprofile`
with the following contents. This file is ignored both by `git` and also
in `.Rbuildingore`.

``` r
Sys.setenv(TESTTHAT_CPUS = 8)
```
# Version 0.9.1
- Updated tests to accommodate upcoming changes to the rlang package.
  No change to test coverage was made.

# Version 0.9.0
- Added the vignette `cmstatr_Validation`
- Updated the expected value of the order statistic of a normally
  distributed variable in the implementation of `hk_ext_z_j_opt`.
  This affects the Basis values computed by `basis_hk_ext` when
  `method="optimum-order"`. Both the new and old implementations appear to
  perform equally well. See the vignette `hk_ext` for more information.
- Added the function `nested_data_plot` for producing nested data plots.
- Added the vignette `hk_ext`
- Updated the vignette `cmstatr_Graphing` to show some examples of the use
  of `nested_data_plot`.
- Added the additional column `batch` to the `carbon.data.2` example data set.
- In `k_factor_normal`, suppress warnings emitted by `qt` when the non-central
  parameter is large.
- Updated the test to use `testthat` edition 3.

# Version 0.8.0
- Updated `basis_anova` so that in cases where the between-batch variance
  is small compared with the within-batch variance, a tolerance factor
  that doesn't consider the structure of the data is used. This matches the
  recommendation of Vangel (1992).
- Added the alias `override="all"` to allow overriding all applicable
  diagnostic tests that are automatically run by the `basis_...` functions.
- Improved documentation of diagnostic tests
- Added `na.rm` argument to `cv` with identical behavior to the `na.rm`
  argument of `mean` and `sd`.
- Fixed bug causing `maximum_normed_residual` to fail with small data sets
  where all but two observations would be considered outliers.
- When diagnostic tests produce an error (when automatically run by the
  `basis_...` functions), the error message now identifies which test
  produced the error.

# Version 0.7.1
- Fixed bug in `glance.equiv_mean_extremum` where it would include empty
  values when a sample was not specified.
- Moved `dplyr` from Suggests to Depends. It is expected that nearly all
  users will use this package in their workflow, and a future version of
  `cmstatr` will also rely on functionality from `dplyr`.
- Changed tests and vignettes such that tests and vignette code
  is not re-run when the necessary packages are not available. Test coverage
  and re-building of vignettes is unchanged when all packages in Depends and
  Suggests are available.

# Version 0.7.0
- Added optional argument to `glance.basis` to add diagnostic test results
  to resulting `data.frame`

# Version 0.6.0
- Improved the documentation for several functions
- Made minor formatting changes to the `print` methods for:
  - `ad_ksample`
  - `anderson_darling`
  - `basis`
  - `equiv_mean_extremum`
  - `equiv_chage_mean`
  - `levene_test`
  - `maximum_normed_residual`
- Added `alpha` into the `mnr` object, and updated `print` and `glance`
  methods to show the value of `alpha` specified by the user

# Version 0.5.2
- Internally use `vapply` instead of `sapply` to improve code safety
- Increased coverage of unit tests

# Version 0.5.1
- Fixed the title of the graphing vignette

# Version 0.5.0
- Renamed `transform_mod_cv_2` to `transform_mod_cv_ad` to better describe
  the purpose of this function.
- Removed the optional argument from `transform_mod_cv`. Now if several
  groups are to be transformed separately, this needs to be done explicitly
  using `dplyr::group_by` or a similar strategy.
- Fixed bug related to the automated diagnostic tests of pooled basis methods
  when `modcv = TRUE`. Previously, the diagnostic tests were performed with
  the unmodified data. After this bug fix, the the data after the modified
  CV transform is used for the diagnostic tests.
- Added `stat` extensions to `ggplot2`:
  - `stat_normal_surv_func` to plot a normal survival function based on
    the data given
  - `stat_esf` to plot an empirical survival function
- Updated cmstatr_Tutorial vignette
- Created cmstatr_Graphing vignette
- Various documentation improvements

# Version 0.4.0
- Added automated diagnostic tests to basis_... methods
- Updated argument names for functions:
  - `transform_mod_cv`
  - `transform_mod_cv_2`
  - `normalize_group_mean`
- Updated cmstatr_Tutorial vignette


# Version 0.3.0
- Added modified CV functionality
- Added glance and augment methods for most objects
- Added function for calculating CV of a sample
- Breaking changes:
  - Renamed function `basis_nonparametric_large_sample` to
    `basis_nonpara_large_sample`
  - Renamed function `nonparametric_binomial_rank` to
    `nonpara_binomial_rank`

# Version 0.2.0
- Added ANOVA basis calculation
- Added non-parametric basis calculations

# Version 0.1.0
- Initial release
This re-submission does not change functionality, but the test suite has
been updated to allow for an upcoming change to rlang when version 1.0.0
of that package is released.

## Test environments
- win-builder (devel, release, oldrelease)
- local Ubuntu 20.04, R 4.1.1
- GitHub Action runners:
  - Windows, R 4.1.1
  - MacOS, R 4.1.1
  - Ubuntu 20.04, R 4.1.1
  - Ubuntu 18.04, R 4.0.5

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
There are no downstream dependencies.
