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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/"
)
```

# cmstatr <img src="man/figures/logo.png" align="right" alt="" width="120" />

<!-- badges: start -->
[![R build status](https://github.com/cmstatr/cmstatr/workflows/R-CMD-check/badge.svg)](https://github.com/cmstatr/cmstatr/actions?workflow=R-CMD-check)
[![`Codecov` test coverage](https://codecov.io/gh/cmstatr/cmstatr/branch/master/graph/badge.svg)](https://codecov.io/gh/cmstatr/cmstatr?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02265/status.svg)](https://doi.org/10.21105/joss.02265)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/cmstatr)](https://cran.r-project.org/package=cmstatr)
[![](https://cranlogs.r-pkg.org/badges/cmstatr)](https://cran.r-project.org/package=cmstatr)
<!-- badges: end -->

# What It Does
The `cmstatr` package provides functions for performing statistical analysis
of composite material data. The statistical methods implemented are those
described in [CMH-17-1G](https://www.cmh17.org/).
This package focuses on calculating basis values (lower tolerance
bounds) for material strength properties, as well as performing the
associated diagnostic tests. Functions are also provided for testing
for equivalency between alternate samples and the "qualification"
or "baseline" samples.

Additional details about the package are available in the paper by
Kloppenborg (2020,
[https://doi.org/10.21105/joss.02265](https://doi.org/10.21105/joss.02265)).

# Installation
To install `cmstatr` from CRAN, simply run:

```r
install.packages("cmstatr")
```

If you want the latest development version, you can install it
from `github` using `devtools`. This will also install the dependencies
required to build the vignettes. Optionally, change the value of the
argument `ref` to install `cmstatr` from a different branch of the
repository.

```r
install.packages(c("devtools", "rmarkdown", "dplyr", "tidyr"))
devtools::install_github("cmstatr/cmstatr", build_vignettes = TRUE,
                         ref = "master",
                         build_opts = c("--no-resave-data", "--no-manual"))
```

# Usage
To compute a B-Basis value from an example data set packaged
with `cmstatr` you can do the following:

```{r message=FALSE}
library(dplyr)
library(cmstatr)

carbon.fabric.2 %>%
  filter(test == "FC") %>%
  filter(condition == "RTD") %>%
  basis_normal(strength, batch)
```

For more examples of usage of the `cmstatr` package,
see the tutorial vignette, which can be
[viewed online](https://www.cmstatr.net/articles/cmstatr_Tutorial.html),
or can be loaded as follows, once the package is installed:

```r
vignette("cmstatr_Tutorial")
```

There is also a vignette showing some examples of the types of graphs
that are typically produced when analyzing composite materials.
You can view this
[vignette online](https://www.cmstatr.net/articles/cmstatr_Graphing.html),
or you can load this vignette with:

```r
vignette("cmstatr_Graphing")
```

# Philosophical Notes
This package expects
[`tidy data`](https://www.jstatsoft.org/article/view/v059i10).
That is, individual observations should be in rows and variables in columns.

Where possible, this package uses general solutions. Look-up tables are avoided
wherever possible.

# Issues
If you've found a bug, please open an issue in this repository and
describe the bug. Please
include a [reproducible example](https://reprex.tidyverse.org/) of the bug.
If you're able to fix the bug, you can do so by submitting a pull request.

If your bug is related to a particular data set, sharing that data set will
help to fix the bug. If you cannot share the data set, please strip any
identifying information and optionally scale the data by an unspecified
factor so that the bug can be reproduced and diagnosed.


# Contributing
Contributions to `cmstatr` are always welcomed. For small changes (fixing typos
or improving the documentation), go ahead and submit a pull request. For more
significant changes, such as new features, please discuss the proposed change
in an issue first.

## Contribution Guidelines
- Please create a git branch for each pull request (PR)
- Before submitting a pull request, please make sure that `R CMD check`
  passes with no errors, warnings or notes
- New and modified code should follow the style guide enforced by the
  [`lintr`](https://cran.r-project.org/package=lintr)
  package
- Document all exported functions using
  [`roxygen2`](https://cran.r-project.org/package=roxygen2)
- Write tests using [`testthat`](https://cran.r-project.org/package=testthat).
  If your contribution fixes a bug, then the test(s) that you add should fail
  before your bug-fix patch is applied and should pass after the code is
  patched.
- For changes that affect the user, add a bullet at the top of `NEWS.md` below
  the current development version
  
## Development
Testing is performed using `testthat`. Edition 3 of that package is used and
parallel processing enabled. If you wish to use more than two CPUs, set the
environment variable `TESTTHAT_CPUS` to the number of CPUs that you want to
use. One way of doing this is to create the file `.Rprofile` with the following
contents. This file is ignored both by `git` and also in `.Rbuildingore`.

```r
Sys.setenv(TESTTHAT_CPUS = 8)
```
---
title: |
  cmstatr: An R Package for Statistical Analysis of Composite Material Data
tags:
  - R
  - statistics
  - composite materials
  - material science
authors:
  - name: Stefan Kloppenborg
    orcid: 0000-0002-1908-5214
    affiliation: 1
affiliations:
 - name: Comtek Advanced Structures, Ltd.
   index: 1
date: "6/24/2020"
bibliography: paper.bib
output:
  md_document:
    preserve_yaml: true
    variant: markdown
    pandoc_args: "--atx-headers"
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Summary
Strength data for composite materials used in aerospace applications,
such as carbon fiber and fiberglass reinforced composites, are normally analyzed 
using statistical methods because of the 
inherent variability in the constituent materials and in the processing.
The design standards for civil aviation requires that the probability
of structural failure due to this variability must be minimized, and to do
so, the designer must use what are called "Design Values" for each material 
in stress analyses and ensure they
exceed the actual stresses experienced by those materials in service.
Design Values are set based on the one-sided lower confidence bound of the
material strength. For some types of structure, the content of this
lower confidence bound is $99\%$ with a confidence level of $95\%$;
in this case, the confidence bound is referred to as A-Basis.
For some other types of structure, the content of the lower
confidence bound is instead $90\%$ with a confidence level of $95\%$;
in this case, the confidence bound is referred to as B-Basis.
The statistical 
methods for calculating these 
basis values are outlined in Composite Materials Handbook, Volume 1, 
Revision G, or CMH-17-1G in short [@CMH171G].
The use of these methods is widely accepted by industry and civil 
aviation regulators.

Design Values are often adjusted to account for anticipated in-service
damage and other factors, however those adjustments are outside the scope
of the present software package.

For a detailed discussion of the theory and applications of tolerance bounds,
the reader is referred to @Meeker_Hahn_Escobar_2017 or
@Krishnamoorthy_Mathew_2008.

Currently, many users use MS Excel spreadsheets to perform these
analyses. The MS Excel spreadsheets typically used, such as 
`STAT-17` [@STAT-17],
`ASAP` [@Raju_Tomblin_2008] and `CMH17-STATS` [@CMH17-STATS], use 
password-protected `VBA` macros to
perform the computations. As such, the code cannot be audited by 
the user. 
`cmstatr` is an R package that addresses this issue by implementing 
the same statistical analysis techniques 
found in CMH-17-1G in an open-source environment.


# Statement of Need
The purpose of `cmstatr` is to:

- Provide a consistent user interface for computing A- and B-Basis
  values and performing the related diagnostic tests in the `R`
  programming environment
- Allow auditing of the code used to compute A- and B-Basis values,
  and performing the related diagnostic tests
- Enable users to automate computation workflows or to perform
  simulation studies


# Implementation Goals
`cmstatr` aims give a consistent interface for the user. Most functions
are written to work with the `tidyverse` [@tidyverse] and most functions
have similar argument lists. The intent is to make the package easy to
learn and use.

The implementation of `cmstatr` also aims to avoid the use of lookup 
tables
that are prevalent in calculation spreadsheets and minimize the use of 
approximations. 
While this decision leads to increased computation time, the typically small 
data sets (tens to hundreds
of observations) associated with composite material test data and the
speed of modern computers make this practical for interactive programming.

# Example Usage
Normally, to use `cmstatr` the user will load `cmstatr` itself as well as
the `tidyverse` package.

```{r message=FALSE}
library(cmstatr)
library(tidyverse)
```

`cmstatr` contains some example data sets, which can be used to demonstrate
the features of the package. One of those data sets --- `carbon.fabric.2` ---
will be used in the following example. This data set contains results from
several mechanical tests of a typical composite material, and contains the
typical measurements obtained from a test lab. In the following
examples, results from tension testing in the warp fiber direction (`WT`) 
will be used.
Part of this data set is shown below.

```{r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  head(10)
```

One common task is to calculate B-Basis values. There are several
statistical methods for doing so, depending on the distribution of
the data.
The single-point basis functions automatically perform the following
diagnostic tests:

- the maximum normed residual test for outliers within a batch [@CMH171G],
- the Anderson--Darling k-Sample test to check if batches are drawn from
  the same (unspecified) distribution [@Scholz_Stephens_1987], 
- the maximum normed residual test for outliers within the data, and 
- the Anderson--Darling test for a particular probability 
  distribution [@Lawless_1982].

Assuming that the data from the warp tension (WT) tested at 
elevated-temperature/wet condition (ETW) follows
a normal distribution, then this can be done using the function.
Note that all of the functions in `cmstatr` that compute
basis values default to computing tolerance bounds with a
content of $p=0.9$ and a confidence of $conf=0.95$, or
B-Basis.

```{r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  filter(condition == "ETW") %>%
  basis_normal(strength, batch)
```

All of the various basis functions perform diagnostic tests for each of the 
statistical tests mentioned above. If any
of the diagnostic tests failed, a warning is shown and the test failure
is also recorded in the returned object (and shown in that object's
`print` method). In the example above, the output shows that the
Anderson--Darling test for normality [@Lawless_1982] rejects the hypothesis
that the data is drawn from a normal distribution. 

Two non-parametric basis calculations, based on @Guenther_1970 and
@Vangel_1994 are also implemented in `cmstatr`. These functions 
perform the same
diagnostic tests, but omit the Anderson--Darling test for a particular
distribution.

The diagnostic test can be run directly using `cmstatr` as well.
For example, the failed diagnostic test above can be run as follows:

```{r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  filter(condition == "ETW") %>%
  anderson_darling_normal(strength)
```

If the failure of a diagnostic test is decided to be acceptable,
the test result can be overridden to hide the warning in the basis function 
output:

```{r}
carbon.fabric.2 %>%
  filter(test == "WT") %>%
  filter(condition == "ETW") %>%
  basis_normal(strength, batch, override = c("anderson_darling_normal"))
```

`cmstatr` also provides functions for calculating basis values from
data pooled across multiple testing environments, as recommended by [@CMH171G]. 
There are two methods of pooling: both calculate a measure of global variation
from the data in all tested environmental conditions, 
then they calculate a basis value for each condition using the measure of
the global variation and the individual mean values of each 
condition. The two methods of pooling have different underlying assumptions
and hence the data must pass a different set diagnostic tests for each of
the two functions.

# Validation and Comparison With Existing Tools
Where possible, `cmstatr` has been verified against the examples given
in the articles in which the statistical methods were published. Unit
tests have been written so that this verification is re-checked routinely
to prevent unintended regressions. Agreement between `cmstatr` and the
examples in the original articles is within the expected numeric accuracy.

`cmstatr` has also been verified against existing software, such as
`STAT-17` [@STAT-17], `ASAP` [@Raju_Tomblin_2008] and
`CMH17-STATS` [@CMH17-STATS] using several example data sets.
Agreement between `cmstatr` and the other software is generally good,
but some results differ slightly, likely due to various approximations used
in the software.
Comparison between `cmstatr` and the other software is performed
within various unit tests to guard against future regressions.

The tests are automatically run each time a change is made to the code
of `cmstatr` using a continuous integration service.


# Reproducibility
It is envisioned that many users of `cmstatr` will use it within an 
R Notebook or a Jupyter Notebook. It is further envisioned that this notebook
will be directly converted into the statistical analysis report. If this is
done, the reader of the statistical report will be able to verify all of the
detailed steps used in the statistical analysis.

# Acknowledgement
The author would like to thank Mr. Billy Cheng for his contributions to 
`cmstatr` and this paper. The author would also like to thank Comtek
Advanced Structures Ltd. for its support in developing and releasing
the `cmstatr` package.

# References
---
title: "Extended Hanson-Koopmans"
author: "Stefan Kloppenborg"
date: "2021-09-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Extended Hanson-Koopmans}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
csl: ieee.csl
references:
- id: Hanson1964
  type: article
  author:
    - given: D. L.
      family: Hanson
    - given: L. H.
      family: Koopmans
  title: Tolerance Limits for the Class of Distributions with Increasing Hazard Rates
  container-title: The Annals of Mathematical Statistics
  volume: "35"
  issue: "4"
  page: 1561-1570
  issued:
    year: "1964"
  DOI: 10.1214/aoms/1177700380
- id: Vangel1994
  type: article
  author:
    - given: Mark
      family: Vangel
  title: One-Sided Nonparametric Tolerance Limits
  container-title: Communications in Statistics - Simulation and Computation
  volume: "23"
  issue: "4"
  page: 1137-1154
  issued:
    year: "1994"
  DOI: 10.1080/03610919408813222
- id: Harter1961
  type: article
  author:
    - given: H. Leon
      family: Harter
  title: Expected values of normal order statistics
  container-title: Biometrika
  volume: "48"
  issue: 1/2
  page: 151-165
  issued:
    year: "1961"
  DOI: https://doi.org/10.2307/2333139
- id: CMH-17-1G
  type: report
  number: CMH-17-1G
  title: Composite Materials Handbook, Volume 1. Polymer Matrix Composites
    Guideline for Characterization of Structural Materials
  publisher: SAE International
  issued:
    year: "2012"
    month: "03"
---



In this vignette, we'll use the following packages:


```r
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
```


The Extended Hanson--Koopmans method is a nonparametric method of
determining tolerance limits (such as A- or B-Basis values).
This method does not assume
any particular distribution, but does require that
$-\log\left(F\right)$ is convex, where $F$ is the cumulative distribution
function (CDF) of the distribution.

The functions `kh_ext_z`, `hk_ext_z_j_opt` and `basis_kh_ext` in `cmstatr`
are based on the Extended Hanson--Koopmans method, developed by
Vangel [@Vangel1994]. This is an extension of the method published in
[@Hanson1964].

Tolerance limits (Basis values) calculated using the Extended
Hanson--Koopmans method are calculated based on two order statistics [^1], $i$
and $j$, and a factor, $z$.
The function `hk_ext_z_j_opt` and the function
`basis_kh_ext` (with `method = "optimum-order"`) set the first of
these order statistics to the first (lowest) order statistic,
and a second order statistic determined by minimizing the following
function:

$$
\left| z E\left(X_{\left(1\right)}\right)
+ \left(1 - z\right) E\left(X_{\left(j\right)}\right)
- \Phi\left(p\right)\right|
$$

Where $E\left(X_{(i)}\right)$ is the expected value of the
$i$`th` order
statistic for a sample drawn from the standard normal distribution,
and $\Phi\left(p\right)$
is the CDF of a standard normal distribution for the content of the tolerance
limit (i.e. $p=0.9$ for B-Basis).

[^1]: The $i$`th` order statistic is the $i$`th`
lowest value in the sample.

The value of $z$ is calculated based on the sample size, $n$, the two order
statistics $i$ and $j$, the content $p$ and the confidence. The calculation
is performed using the method in [@Vangel1994] and
implemented in `kh_ext_z`. The value of $j$ is very sensitive to the way that
the expected value of the order statistics is calculated, and may be sensitive
to numerical precision.

In version 0.8.0 of `cmstatr` and prior, the expected value of an order
statistic for a sample drawn from a standard normal distribution was
determined in a crude way. After version 0.8.0, the method in
[@Harter1961] is used. These method produce different values of $j$ for certain
sample sizes. Additionally, a table of $j$ and $z$ values for various
sample sizes is published in CMH-17-1G[^2] [@CMH-17-1G].
This table gives slightly different
values of $j$ for some sample sizes.

[^2]: Note that CMH-17-1G uses the symbols $r$ and $k$ instead of
      $j$ and $z$.

The values of $j$ and $z$ produced by `cmstatr` in version 0.8.0 and before,
the values produced after version 0.8.0 and the value published in CMH-17-1G
are shown below. All of these values are for B-Basis (90% content,
95% confidence).



```r
factors <- tribble(
  ~n, ~j_pre_080, ~z_pre_080, ~j_post_080, ~z_post_080, ~j_cmh, ~z_cmh,
  2, 2, 35.1768141883907, 2, 35.1768141883907, 2, 35.177,
  3, 3, 7.85866787768029, 3, 7.85866787768029, 3, 7.859,
  4, 4, 4.50522447199018, 4, 4.50522447199018, 4, 4.505,
  5, 4, 4.10074820079326, 4, 4.10074820079326, 4, 4.101,
  6, 5, 3.06444416024793, 5, 3.06444416024793, 5, 3.064,
  7, 5, 2.85751000593839, 5, 2.85751000593839, 5, 2.858,
  8, 6, 2.38240998122575, 6, 2.38240998122575, 6, 2.382,
  9, 6, 2.25292053841772, 6, 2.25292053841772, 6, 2.253,
  10, 7, 1.98762060673102, 6, 2.13665759924781, 6, 2.137,
  11, 7, 1.89699586212496, 7, 1.89699586212496, 7, 1.897,
  12, 7, 1.81410756892749, 7, 1.81410756892749, 7, 1.814,
  13, 8, 1.66223343216608, 7, 1.73773765993598, 7, 1.738,
  14, 8, 1.59916281901889, 8, 1.59916281901889, 8, 1.599,
  15, 8, 1.54040000806181, 8, 1.54040000806181, 8, 1.54,
  16, 9, 1.44512878109546, 8, 1.48539432060546, 8, 1.485,
  17, 9, 1.39799975474842, 9, 1.39799975474842, 8, 1.434,
  18, 9, 1.35353033609361, 9, 1.35353033609361, 9, 1.354,
  19, 10, 1.28991705486727, 9, 1.31146980117942, 9, 1.311,
  20, 10, 1.25290765871981, 9, 1.27163203813793, 10, 1.253,
  21, 10, 1.21771654027026, 10, 1.21771654027026, 10, 1.218,
  22, 11, 1.17330587650406, 10, 1.18418267046374, 10, 1.184,
  23, 11, 1.14324511741536, 10, 1.15218647199938, 11, 1.143,
  24, 11, 1.11442082880151, 10, 1.12153586685854, 11, 1.114,
  25, 11, 1.08682185727661, 11, 1.08682185727661, 11, 1.087,
  26, 11, 1.06032912052507, 11, 1.06032912052507, 11, 1.06,
  27, 12, 1.03307994274081, 11, 1.03485308510789, 11, 1.035,
  28, 12, 1.00982188136729, 11, 1.01034609051393, 12, 1.01
)
```

For the sample sizes where $j$ is the same for each approach, the values
of $z$ are also equal within a small tolerance.







```r
factors %>%
  filter(j_pre_080 == j_post_080 & j_pre_080 == j_cmh)
#> # A tibble: 16 × 7
#>        n j_pre_080 z_pre_080 j_post_080 z_post_080 j_cmh z_cmh
#>    <dbl>     <dbl>     <dbl>      <dbl>      <dbl> <dbl> <dbl>
#>  1     2         2     35.2           2      35.2      2 35.2 
#>  2     3         3      7.86          3       7.86     3  7.86
#>  3     4         4      4.51          4       4.51     4  4.50
#>  4     5         4      4.10          4       4.10     4  4.10
#>  5     6         5      3.06          5       3.06     5  3.06
#>  6     7         5      2.86          5       2.86     5  2.86
#>  7     8         6      2.38          6       2.38     6  2.38
#>  8     9         6      2.25          6       2.25     6  2.25
#>  9    11         7      1.90          7       1.90     7  1.90
#> 10    12         7      1.81          7       1.81     7  1.81
#> 11    14         8      1.60          8       1.60     8  1.60
#> 12    15         8      1.54          8       1.54     8  1.54
#> 13    18         9      1.35          9       1.35     9  1.35
#> 14    21        10      1.22         10       1.22    10  1.22
#> 15    25        11      1.09         11       1.09    11  1.09
#> 16    26        11      1.06         11       1.06    11  1.06
```

The sample sizes where the value of $j$ differs are as follows:


```r
factor_diff <- factors %>%
  filter(j_pre_080 != j_post_080 | j_pre_080 != j_cmh | j_post_080 != j_cmh)
factor_diff
#> # A tibble: 11 × 7
#>        n j_pre_080 z_pre_080 j_post_080 z_post_080 j_cmh z_cmh
#>    <dbl>     <dbl>     <dbl>      <dbl>      <dbl> <dbl> <dbl>
#>  1    10         7      1.99          6       2.14     6  2.14
#>  2    13         8      1.66          7       1.74     7  1.74
#>  3    16         9      1.45          8       1.49     8  1.48
#>  4    17         9      1.40          9       1.40     8  1.43
#>  5    19        10      1.29          9       1.31     9  1.31
#>  6    20        10      1.25          9       1.27    10  1.25
#>  7    22        11      1.17         10       1.18    10  1.18
#>  8    23        11      1.14         10       1.15    11  1.14
#>  9    24        11      1.11         10       1.12    11  1.11
#> 10    27        12      1.03         11       1.03    11  1.03
#> 11    28        12      1.01         11       1.01    12  1.01
```

While there are differences in the three implementations, it's not clear
how much these differences will matter in terms of the tolerance limits
calculated. This can be investigated through simulation.

# Simulation with Normally Distributed Data
First, we'll generate a large number (10,000) of samples of sample size
$n$ from a normal distribution. Since we're generating the samples, we
know the true population parameters, so can calculate the true population
quantiles. We'll use the three sets of $j$ and $z$ values to compute
tolerance limits and compared those tolerance limits to the population
quantiles. The proportion of the calculated tolerance limits below the
population quantiles should be equal to the selected confidence. We'll
restrict the simulation study to the sample sizes where the values of
$j$ and $z$ differ in the three implementations of this method, and we'll
consider B-Basis (90% content, 95% confidence).


```r
mu_normal <- 100
sd_normal <- 6

set.seed(1234567)  # make this reproducible

tolerance_limit <- function(x, j, z) {
  x[j] * (x[1] / x[j]) ^ z
}

sim_normal <- pmap_dfr(factor_diff, function(n, j_pre_080, z_pre_080,
                                             j_post_080, z_post_080,
                                             j_cmh, z_cmh) {
  map_dfr(1:10000, function(i_sim) {
    x <- sort(rnorm(n, mu_normal, sd_normal))
    tibble(
      n = n,
      b_pre_080 = tolerance_limit(x, j_pre_080, z_pre_080),
      b_post_080 = tolerance_limit(x, j_post_080, z_post_080),
      b_cmh = tolerance_limit(x, j_cmh, z_cmh),
      x = list(x)
    )
  }
  )
})
sim_normal
#> # A tibble: 110,000 × 5
#>        n b_pre_080 b_post_080 b_cmh x         
#>    <dbl>     <dbl>      <dbl> <dbl> <list>    
#>  1    10      78.4       77.7  77.7 <dbl [10]>
#>  2    10      82.8       82.0  82.0 <dbl [10]>
#>  3    10      83.3       83.0  83.0 <dbl [10]>
#>  4    10      78.4       77.2  77.2 <dbl [10]>
#>  5    10      87.3       86.6  86.6 <dbl [10]>
#>  6    10      92.3       93.2  93.2 <dbl [10]>
#>  7    10      75.2       77.9  77.9 <dbl [10]>
#>  8    10      75.4       73.9  73.9 <dbl [10]>
#>  9    10      75.5       75.1  75.1 <dbl [10]>
#> 10    10      76.4       78.4  78.4 <dbl [10]>
#> # … with 109,990 more rows
```

One can see that the tolerance limits calculated with each set of
factors for (most) data sets is different. However, this does not necessarily
mean that any set of factors is more or less correct.

The distribution of the tolerance limits for each sample size is as follows:


```r
sim_normal %>%
  pivot_longer(cols = b_pre_080:b_cmh, names_to = "Factors") %>%
  ggplot(aes(x = value, color = Factors)) +
  geom_density() +
  facet_wrap(n ~ .) +
  theme_bw() +
  ggtitle("Distribution of Tolerance Limits for Various Values of n")
```

![plot of chunk distribution-normal](distribution-normal-1.png)

For all samples sizes, the distribution of tolerance limits is actually
very similar between all three sets of factors.

The true population quantile can be calculated as follows:


```r
x_p_normal <- qnorm(0.9, mu_normal, sd_normal, lower.tail = FALSE)
x_p_normal
#> [1] 92.31069
```

The proportion of calculated tolerance limit values that are below the
population quantile can be calculated as follows. We see that the in all
cases the tolerance limits are all conservative, and also that each
set of factors produce similar levels of conservatism.


```r
sim_normal %>%
  mutate(below_pre_080 = b_pre_080 < x_p_normal,
         below_post_080 = b_post_080 < x_p_normal,
         below_cmh = b_cmh < x_p_normal) %>%
  group_by(n) %>%
  summarise(
    prop_below_pre_080 = sum(below_pre_080) / n(),
    prop_below_post_080 = sum(below_post_080) / n(),
    prop_below_cmh = sum(below_cmh) / n()
  )
#> # A tibble: 11 × 4
#>        n prop_below_pre_080 prop_below_post_080 prop_below_cmh
#>    <dbl>              <dbl>               <dbl>          <dbl>
#>  1    10              0.984               0.980          0.980
#>  2    13              0.979               0.975          0.975
#>  3    16              0.969               0.967          0.967
#>  4    17              0.973               0.973          0.971
#>  5    19              0.962               0.961          0.961
#>  6    20              0.964               0.962          0.964
#>  7    22              0.961               0.960          0.960
#>  8    23              0.960               0.959          0.960
#>  9    24              0.962               0.961          0.962
#> 10    27              0.954               0.953          0.954
#> 11    28              0.952               0.952          0.952
```

# Simulation with Weibull Data

Next, we'll do a similar simulation using data drawn from a Weibull
distribution. Again, we'll generate 10,000 samples for each sample size.


```r
shape_weibull <- 60
scale_weibull <- 100

set.seed(234568)  # make this reproducible

sim_weibull <- pmap_dfr(factor_diff, function(n, j_pre_080, z_pre_080,
                                              j_post_080, z_post_080,
                                              j_cmh, z_cmh) {
  map_dfr(1:10000, function(i_sim) {
    x <- sort(rweibull(n, shape_weibull, scale_weibull))
    tibble(
      n = n,
      b_pre_080 = tolerance_limit(x, j_pre_080, z_pre_080),
      b_post_080 = tolerance_limit(x, j_post_080, z_post_080),
      b_cmh = tolerance_limit(x, j_cmh, z_cmh),
      x = list(x)
    )
  }
  )
})
sim_weibull
#> # A tibble: 110,000 × 5
#>        n b_pre_080 b_post_080 b_cmh x         
#>    <dbl>     <dbl>      <dbl> <dbl> <list>    
#>  1    10      95.3       95.1  95.1 <dbl [10]>
#>  2    10      88.5       88.3  88.3 <dbl [10]>
#>  3    10      89.7       89.3  89.3 <dbl [10]>
#>  4    10      94.7       94.4  94.4 <dbl [10]>
#>  5    10      96.9       96.9  96.9 <dbl [10]>
#>  6    10      93.6       93.2  93.2 <dbl [10]>
#>  7    10      86.1       85.5  85.5 <dbl [10]>
#>  8    10      91.9       91.9  91.9 <dbl [10]>
#>  9    10      93.7       93.4  93.4 <dbl [10]>
#> 10    10      90.9       90.4  90.4 <dbl [10]>
#> # … with 109,990 more rows
```


The distribution of the tolerance limits for each sample size is as follows.
Once again, we see that the distribution of tolerance limits is nearly
identical when each of the three sets of factors are used.


```r
sim_weibull %>%
  pivot_longer(cols = b_pre_080:b_cmh, names_to = "Factors") %>%
  ggplot(aes(x = value, color = Factors)) +
  geom_density() +
  facet_wrap(n ~ .) +
  theme_bw() +
  ggtitle("Distribution of Tolerance Limits for Various Values of n")
```

![plot of chunk distribution-Weibull](distribution-Weibull-1.png)

The true population quantile can be calculated as follows:


```r
x_p_weibull <- qweibull(0.9, shape_weibull, scale_weibull, lower.tail = FALSE)
x_p_weibull
#> [1] 96.31885
```

The proportion of calculated tolerance limit values that are below the
population quantile can be calculated as follows. We see that the in all
roughly 95% or more of the tolerance limits calculated for each sample
is below the population quantile. We also see very similar proportions
for each of the three sets of factors considered.


```r
sim_weibull %>%
  mutate(below_pre_080 = b_pre_080 < x_p_weibull,
         below_post_080 = b_post_080 < x_p_weibull,
         below_cmh = b_cmh < x_p_weibull) %>%
  group_by(n) %>%
  summarise(
    prop_below_pre_080 = sum(below_pre_080) / n(),
    prop_below_post_080 = sum(below_post_080) / n(),
    prop_below_cmh = sum(below_cmh) / n()
  )
#> # A tibble: 11 × 4
#>        n prop_below_pre_080 prop_below_post_080 prop_below_cmh
#>    <dbl>              <dbl>               <dbl>          <dbl>
#>  1    10              0.97                0.965          0.965
#>  2    13              0.966               0.964          0.964
#>  3    16              0.959               0.959          0.959
#>  4    17              0.961               0.961          0.96 
#>  5    19              0.957               0.956          0.956
#>  6    20              0.955               0.954          0.955
#>  7    22              0.953               0.952          0.952
#>  8    23              0.950               0.950          0.950
#>  9    24              0.953               0.953          0.953
#> 10    27              0.952               0.951          0.951
#> 11    28              0.950               0.950          0.950
```


# Conclusion
The values of $j$ and $z$ computed by the `kh_Ext_z_j_opt` function differs
for certain samples sizes ($n$) before and after version 0.8.0. Furthermore,
for certain sample sizes, these values differ from those published in
CMH-17-1G. The simulation study presented in this vignette shows that
the tolerance limit (Basis value) might differ for any individual sample
based on which set of $j$ and $z$ are used. However, each set of factors
produces tolerance limit factors that are either correct or conservative.
These three methods have very similar performance, and tolerance limits
produced with any of these three methods are equally valid.

# Session Info
This vignette is computed in advance. A system with the following configuration
was used:


```r
sessionInfo()
#> R version 4.1.1 (2021-08-10)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 20.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C               LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
#>  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8    LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] tidyr_1.1.4   purrr_0.3.4   ggplot2_3.3.5 dplyr_1.0.7  
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.1.1  xfun_0.26         remotes_2.4.0     colorspace_2.0-2  vctrs_0.3.8       generics_0.1.0   
#>  [7] testthat_3.0.4    htmltools_0.5.2   usethis_2.0.1     yaml_2.2.1        utf8_1.2.2        rlang_0.4.11.9001
#> [13] pkgbuild_1.2.0    pillar_1.6.3      glue_1.4.2        withr_2.4.2       DBI_1.1.1         waldo_0.3.1      
#> [19] cmstatr_0.9.0     sessioninfo_1.1.1 lifecycle_1.0.1   stringr_1.4.0     munsell_0.5.0     gtable_0.3.0     
#> [25] devtools_2.4.2    memoise_2.0.0     evaluate_0.14     labeling_0.4.2    knitr_1.35        callr_3.7.0      
#> [31] fastmap_1.1.0     ps_1.6.0          curl_4.3.1        fansi_0.5.0       highr_0.9         scales_1.1.1     
#> [37] kSamples_1.2-9    cachem_1.0.6      desc_1.4.0        pkgload_1.2.2     farver_2.1.0      fs_1.5.0         
#> [43] digest_0.6.28     stringi_1.7.4     processx_3.5.2    SuppDists_1.1-9.5 rprojroot_2.0.2   grid_4.1.1       
#> [49] cli_3.0.1         tools_4.1.1       magrittr_2.0.1    tibble_3.1.4      crayon_1.4.1      pkgconfig_2.0.3  
#> [55] MASS_7.3-54       ellipsis_0.3.2    prettyunits_1.1.1 assertthat_0.2.1  rmarkdown_2.11    rstudioapi_0.13  
#> [61] R6_2.5.1          compiler_4.1.1
```


# References
---
title: "Anderson-Darling k-Sample Test"
author: "Stefan Kloppenborg"
date: "20-Jan-2019"
output: rmarkdown::html_vignette
bibliography: bibliography.json
csl: ieee.csl
vignette: >
  %\VignetteIndexEntry{Anderson-Darling k-Sample Test}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This vignette explores the Anderson--Darling k-Sample test.
CMH-17-1G [@CMH-17-1G] provides a formulation for this test that appears different than the formulation given by Scholz and Stephens in their 1987 paper [@Stephens1987].

Both references use different nomenclature, which is summarized as follows:

Term                                               | CMH-17-1G             | Scholz and Stephens
---------------------------------------------------|-----------------------|---------------------
A sample                                           | $i$                   | $i$
The number of samples                              | $k$                   | $k$
An observation within a sample                     | $j$                   | $j$
The number of observations within the sample $i$   | $n_i$                 | $n_i$
The total number of observations within all samples| $n$                   | $N$
Distinct values in combined data, ordered          | $z_{(1)}$...$z_{(L)}$ | $Z_1^*$...$Z_L^*$
The number of distinct values in the combined data | $L$                   | $L$


Given the possibility of ties in the data, the discrete version of the test must be used
Scholz and Stephens (1987) give the test statistic as:

$$
A_{a k N}^2 = \frac{N - 1}{N}\sum_{i=1}^k \frac{1}{n_i}\sum_{j=1}^{L}\frac{l_j}{N}\frac{\left(N M_{a i j} - n_i B_{a j}\right)^2}{B_{a j}\left(N - B_{a j}\right) - N l_j / 4}
$$


CMH-17-1G gives the test statistic as:

$$
ADK = \frac{n - 1}{n^2\left(k - 1\right)}\sum_{i=1}^k\frac{1}{n_i}\sum_{j=1}^L h_j \frac{\left(n F_{i j} - n_i H_j\right)^2}{H_j \left(n - H_j\right) - n h_j / 4}
$$

By inspection, the CMH-17-1G version of this test statistic contains an extra factor of $\frac{1}{\left(k - 1\right)}$.

Scholz and Stephens indicate that one rejects $H_0$ at a significance level of $\alpha$ when:

$$
\frac{A_{a k N}^2 - \left(k - 1\right)}{\sigma_N} \ge t_{k - 1}\left(\alpha\right)
$$

This can be rearranged to give a critical value:

$$
A_{c r i t}^2 = \left(k - 1\right) + \sigma_N t_{k - 1}\left(\alpha\right)
$$

CHM-17-1G gives the critical value for $ADK$ for $\alpha=0.025$ as:

$$
ADC = 1 + \sigma_n \left(1.96 + \frac{1.149}{\sqrt{k - 1}} - \frac{0.391}{k - 1}\right)
$$

The definition of $\sigma_n$ from the two sources differs by a factor of $\left(k - 1\right)$.

The value in parentheses in the CMH-17-1G critical value corresponds to the interpolation formula for $t_m\left(\alpha\right)$ given in Scholz and Stephen's paper.
It should be noted that this is *not* the student's t-distribution, but rather a distribution referred to as the $T_m$ distribution.

The `cmstatr` package use the package `kSamples` to perform the k-sample Anderson--Darling tests.
This package uses the original formulation from Scholz and Stephens, so the test statistic will differ from that given software based on the CMH-17-1G formulation by a factor of $\left(k-1\right)$.
The conclusions about the null hypothesis drawn, however, will be the same.

# References
---
title: "cmstatr Validation"
author: "Stefan Kloppenborg"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{cmstatr Validation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
csl: ieee.csl
references:
- id: "ASAP2008"
  type: "report"
  number: "ASAP-2008"
  author:
  - given: K.S.
    family: Raju
  - given: J.S.
    family: Tomblin
  title: "AGATE Statistical Analysis Program"
  issued:
    year: 2008
  publisher: Wichita State University
- id: Stephens1987
  type: article
  author:
    - given: F.W.
      family: Scholz
    - given: M.A,
      family: Stephens
  title: K-Sample Anderson-Darling Tests
  container-title: Journal of the American Statistical Association
  volume: "82"
  issue: "399"
  page: 918-924
  issued:
    year: "1987"
    month: "09"
  DOI: 10.1080/01621459.1987.10478517
  URL: https://doi.org/10.1080/01621459.1987.10478517
- id: "CMH-17-1G"
  type: report
  number: "CMH-17-1G"
  title: Composite Materials Handbook, Volume 1. Polymer Matrix Composites
    Guideline for Characterization of Structural Materials
  publisher: SAE International
  issued:
    year: "2012"
    month: "03"
- id: "STAT-17"
  type: report
  number: "STAT-17 Rev 5"
  author:
    - literal: Materials Sciences Corporation
  title: CMH-17 Statistical Analysis for B-Basis and A-Basis Values
  publisher: Materials Sciences Corporation
  publisher-place: Horsham, PA
  issued:
    year: "2008"
    month: "01"
    day: "08"
- id: Vangel1994
  type: article
  author:
    - given: Mark
      family: Vangel
  title: One-Sided Nonparametric Tolerance Limits
  container-title: Communications in Statistics - Simulation and Computation
  volume: "23"
  issue: "4"
  page: 1137-1154
  issued:
    year: "1994"
  DOI: 10.1080/03610919408813222
- id: vangel_lot_2002
  type: article
  author:
    - given: Mark
      family: Vangel
  title: Lot Acceptance and Compliance Testing Using the Sample Mean and an Extremum
  container-title: Technometrics
  volume: "44"
  issue: "3"
  page: 242--249
  issued:
    year: 2002
    month: 8
  DOI: 10.1198/004017002188618428

---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# 1. Introduction

This vignette is intended to contain the same validation that is included
in the test suite within the `cmstatr` package, but in a format that is
easier for a human to read. The intent is that this vignette will include
only those validations that are included in the test suite, but that the
test suite may include more tests than are shown in this vignette.

The following packages will be used in this validation. The version of
each package used is listed at the end of this vignette.

```{r message=FALSE, warning=FALSE}
library(cmstatr)
library(dplyr)
library(purrr)
library(tidyr)
library(testthat)
```

Throughout this vignette, the `testthat` package will be used. Expressions
such as `expect_equal` are used to ensure that the two values are equal
(within some tolerance). If this expectation is not true, the vignette will
fail to build. The tolerance is a relative tolerance: a tolerance of 0.01
means that the two values must be within $1\%$ of each other.
As an example, the following expression checks that the value
`10` is equal to `10.1` within a tolerance of `0.01`. Such an expectation
should be satisfied.


```{r}
expect_equal(10, 10.1, tolerance = 0.01)
```

The `basis_...` functions automatically perform certain diagnostic tests.
When those diagnostic tests are not relevant to the validation, the
diagnostic tests are overridden by passing the argument `override = "all"`.

# 2. Validation Table
The following table provides a cross-reference between the various
functions of the `cmstatr` package and the tests shown within this
vignette. The sections in this vignette are organized by data set.
Not all checks are performed on all data sets.


Function                        | Tests
--------------------------------|---------------------------
`ad_ksample()`                  | [Section 3.1](#cf-ad), [Section 4.1.2](#c11-adk), [Section 6.1](#cl-adk)
`anderson_darling_normal()`     | [Section 4.1.3](#c11-ad), [Section 5.1](#ods-adn)
`anderson_darling_lognormal()`  | [Section 4.1.3](#c11-ad), [Section 5.2](#ods-adl)
`anderson_darling_weibull()`    | [Section 4.1.3](#c11-ad), [Section 5.3](#ods-adw)
`basis_normal()`                | [Section 5.4](#ods-nb)
`basis_lognormal()`             | [Section 5.5](#ods-lb)
`basis_weibull()`               | [Section 5.6](#ods-wb)
`basis_pooled_cv()`             | [Section 4.2.3](#c12-pcv), [Section 4.2.4](#c12-pcvmcv),
`basis_pooled_sd()`             | [Section 4.2.1](#c12-psd), [Section 4.2.2](#c12-psdmcv)
`basis_hk_ext()`                | [Section 4.1.6](#c11-hk-wf), [Section 5.7](#ods-hkb), [Section 5.8](#ods-hkb2)
`basis_nonpara_large_sample()`  | [Section 5.9](#ods-lsnb)
`basis_anova()`                 | [Section 4.1.7](#c11-anova)
`calc_cv_star()`                |
`cv()`                          |
`equiv_change_mean()`           | [Section 5.11](#ods-ecm)
`equiv_mean_extremum()`         | [Section 5.10](#ods-eme)
`hk_ext_z()`                    | [Section 7.3](#pf-hk), [Section 7.4](#pf-hk2)
`hk_ext_z_j_opt()`              | [Section 7.5](#pf-hk-opt)
`k_equiv()`                     | [Section 7.8](#pf-equiv)
`k_factor_normal()`             | [Section 7.1](#pf-kb), [Section 7.2](#pf-ka)
`levene_test()`                 | [Section 4.1.4](#c11-lbc), [Section 4.1.5](#c11-lbb)
`maximum_normed_residual()`     | [Section 4.1.1](#c11-mnr)
`nonpara_binomial_rank()`       | [Section 7.6](#pf-npbinom), [Section 7.7](#pf-npbinom2)
`normalize_group_mean()`        |
`normalize_ply_thickness()`     |
`transform_mod_cv_ad()`         |
`transform_mod_cv()`            |


# 3. `carbon.fabric` Data Set

This data set is example data that is provided with `cmstatr`.
The first few rows of this data are shown below.

```{r}
head(carbon.fabric)
```


## 3.1. Anderson--Darling k-Sample Test {#cf-ad}

This data was entered into ASAP 2008 [@ASAP2008] and the reported
Anderson--Darling k--Sample test statistics were recorded, as were
the conclusions.

The value of the test statistic reported by `cmstatr` and that reported
by ASAP 2008 differ by a factor of $k - 1$, as do the critical values
used. As such, the conclusion of the tests are identical. This is
described in more detail in the
[Anderson--Darling k--Sample Vignette](adktest.html).

When the RTD warp-tension data from this data set is entered into ASAP 2008,
it reports a test statistic of 0.456 and fails to reject the null hypothesis
that the batches are drawn from the same distribution. Adjusting for the
different definition of the test statistic, the results given by `cmstatr`
are very similar.


```{r}
res <- carbon.fabric %>%
  filter(test == "WT") %>%
  filter(condition == "RTD") %>%
  ad_ksample(strength, batch)

expect_equal(res$ad / (res$k - 1), 0.456, tolerance = 0.002)
expect_false(res$reject_same_dist)

res
```

When the ETW warp-tension data from this data set are entered into ASAP 2008,
the reported test statistic is 1.604 and it fails to reject the null
hypothesis that the batches are drawn from the same distribution. Adjusting
for the different definition of the test statistic, `cmstatr` gives nearly
identical results.

```{r}
res <- carbon.fabric %>%
  filter(test == "WT") %>%
  filter(condition == "ETW") %>%
  ad_ksample(strength, batch)

expect_equal(res$ad / (res$k - 1), 1.604, tolerance = 0.002)
expect_false(res$reject_same_dist)

res
```


# 4. Comparison with Examples from CMH-17-1G
## 4.1 Dataset From Section 8.3.11.1.1

CMH-17-1G [@CMH-17-1G] provides an example data set and results from
ASAP [@ASAP2008] and STAT17 [@STAT-17]. This example data set is duplicated
below:

```{r}
dat_8_3_11_1_1 <- tribble(
  ~batch, ~strength, ~condition,
  1, 118.3774604, "CTD", 1, 84.9581364, "RTD", 1, 83.7436035, "ETD",
  1, 123.6035612, "CTD", 1, 92.4891822, "RTD", 1, 84.3831677, "ETD",
  1, 115.2238092, "CTD", 1, 96.8212659, "RTD", 1, 94.8030433, "ETD",
  1, 112.6379744, "CTD", 1, 109.030325, "RTD", 1, 94.3931537, "ETD",
  1, 116.5564277, "CTD", 1, 97.8212659, "RTD", 1, 101.702222, "ETD",
  1, 123.1649896, "CTD", 1, 100.921519, "RTD", 1, 86.5372121, "ETD",
  2, 128.5589027, "CTD", 1, 103.699444, "RTD", 1, 92.3772684, "ETD",
  2, 113.1462103, "CTD", 2, 93.7908212, "RTD", 2, 89.2084024, "ETD",
  2, 121.4248107, "CTD", 2, 107.526709, "RTD", 2, 100.686001, "ETD",
  2, 134.3241906, "CTD", 2, 94.5769704, "RTD", 2, 81.0444192, "ETD",
  2, 129.6405117, "CTD", 2, 93.8831373, "RTD", 2, 91.3398070, "ETD",
  2, 117.9818658, "CTD", 2, 98.2296605, "RTD", 2, 93.1441939, "ETD",
  3, 115.4505226, "CTD", 2, 111.346590, "RTD", 2, 85.8204168, "ETD",
  3, 120.0369467, "CTD", 2, 100.817538, "RTD", 3, 94.8966273, "ETD",
  3, 117.1631088, "CTD", 3, 100.382203, "RTD", 3, 95.8068520, "ETD",
  3, 112.9302797, "CTD", 3, 91.5037811, "RTD", 3, 86.7842252, "ETD",
  3, 117.9114501, "CTD", 3, 100.083233, "RTD", 3, 94.4011973, "ETD",
  3, 120.1900159, "CTD", 3, 95.6393615, "RTD", 3, 96.7231171, "ETD",
  3, 110.7295966, "CTD", 3, 109.304779, "RTD", 3, 89.9010384, "ETD",
  3, 100.078562, "RTD", 3, 99.1205847, "RTD", 3, 89.3672306, "ETD",
  1, 106.357525, "ETW", 1, 99.0239966, "ETW2",
  1, 105.898733, "ETW", 1, 103.341238, "ETW2",
  1, 88.4640082, "ETW", 1, 100.302130, "ETW2",
  1, 103.901744, "ETW", 1, 98.4634133, "ETW2",
  1, 80.2058219, "ETW", 1, 92.2647280, "ETW2",
  1, 109.199597, "ETW", 1, 103.487693, "ETW2",
  1, 61.0139431, "ETW", 1, 113.734763, "ETW2",
  2, 99.3207107, "ETW", 2, 108.172659, "ETW2",
  2, 115.861770, "ETW", 2, 108.426732, "ETW2",
  2, 82.6133082, "ETW", 2, 116.260375, "ETW2",
  2, 85.3690411, "ETW", 2, 121.049610, "ETW2",
  2, 115.801622, "ETW", 2, 111.223082, "ETW2",
  2, 44.3217741, "ETW", 2, 104.574843, "ETW2",
  2, 117.328077, "ETW", 2, 103.222552, "ETW2",
  2, 88.6782903, "ETW", 3, 99.3918538, "ETW2",
  3, 107.676986, "ETW", 3, 87.3421658, "ETW2",
  3, 108.960241, "ETW", 3, 102.730741, "ETW2",
  3, 116.122640, "ETW", 3, 96.3694916, "ETW2",
  3, 80.2334815, "ETW", 3, 99.5946088, "ETW2",
  3, 106.145570, "ETW", 3, 97.0712407, "ETW2",
  3, 104.667866, "ETW",
  3, 104.234953, "ETW"
)
dat_8_3_11_1_1
```

### 4.1.1 Maximum Normed Residual Test {#c11-mnr}
CMH-17-1G Table 8.3.11.1.1(a) provides results of the MNR test from
ASAP for this data set. Batches 2 and 3 of the ETW data is considered
here and the results of `cmstatr` are compared with those published in
CMH-17-1G.

For Batch 2 of the ETW data, the results match those published in the
handbook within a small tolerance. The published test statistic is 2.008.

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW" & batch == 2) %>%
  maximum_normed_residual(strength, alpha = 0.05)

expect_equal(res$mnr, 2.008, tolerance = 0.001)
expect_equal(res$crit, 2.127, tolerance = 0.001)
expect_equal(res$n_outliers, 0)

res
```

Similarly, for Batch 3 of the ETW data, the results of `cmstatr` match
the results published in the handbook within a small tolerance. The published
test statistic is 2.119

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW" & batch == 3) %>%
  maximum_normed_residual(strength, alpha = 0.05)

expect_equal(res$mnr, 2.119, tolerance = 0.001)
expect_equal(res$crit, 2.020, tolerance = 0.001)
expect_equal(res$n_outliers, 1)

res
```


### 4.1.2 Anderson--Darling k--Sample Test {#c11-adk}

For the ETW condition, the ADK test statistic given in [@CMH-17-1G] is
$ADK = 0.793$ and the test concludes that the samples come from the
same distribution. Noting that `cmstatr` uses the definition of the
test statistic given in [@Stephens1987], so the test statistic given
by `cmstatr` differs from that given by ASAP by a factor of $k - 1$,
as described in the [Anderson--Darling k--Sample Vignette](adktest.html).

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW") %>%
  ad_ksample(strength, batch)

expect_equal(res$ad / (res$k - 1), 0.793, tolerance = 0.003)
expect_false(res$reject_same_dist)

res
```

Similarly, for the ETW2 condition, the test statistic given in [@CMH-17-1G]
is $ADK = 3.024$ and the test concludes that the samples come from different
distributions. This matches `cmstatr`

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW2") %>%
  ad_ksample(strength, batch)

expect_equal(res$ad / (res$k - 1), 3.024, tolerance = 0.001)
expect_true(res$reject_same_dist)

res
```

### 4.1.3 Anderson--Darling Tests for Distribution {#c11-ad}

CMH-17-1G Section 8.3.11.2.1 contains results from STAT17 for
the "observed significance level" from the Anderson--Darling test
for various distributions. In this section, the ETW condition from the
present data set is used. The published results are given in the 
following table. The results from `cmstatr` are below and are very
similar to those from STAT17.

Distribution | OSL
-------------|------------
Normal       | 0.006051
Lognormal    | 0.000307
Weibull      | 0.219

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW") %>%
  anderson_darling_normal(strength)
expect_equal(res$osl, 0.006051, tolerance = 0.001)
res
```
```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW") %>%
  anderson_darling_lognormal(strength)
expect_equal(res$osl, 0.000307, tolerance = 0.001)
res
```

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW") %>%
  anderson_darling_weibull(strength)
expect_equal(res$osl, 0.0219, tolerance = 0.002)
res
```


### 4.1.4 Levene's Test (Between Conditions) {#c11-lbc}

CMH-17-1G Section 8.3.11.1.1 provides results from ASAP for Levene's test
for equality of variance between conditions after the ETW and ETW2 conditions
are removed. The handbook shows an F statistic of 0.58, however if this
data is entered into ASAP directly, ASAP gives an F statistic of 0.058,
which matches the result of `cmstatr`.

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition != "ETW" & condition != "ETW2") %>%
  levene_test(strength, condition)
expect_equal(res$f, 0.058, tolerance = 0.01)
res
```

### 4.1.5 Levene's Test (Between Batches) {#c11-lbb}

CMH-17-1G Section 8.3.11.2.2 provides output from STAT17. The
ETW2 condition from the present data set was analyzed by STAT17
and that software reported an F statistic of 0.123 from Levene's
test when comparing the variance of the batches within this condition.
The result from `cmstatr` is similar.

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW2") %>%
  levene_test(strength, batch)
expect_equal(res$f, 0.123, tolerance = 0.005)
res
```

Similarly, the published value of the F statistic for the CTD condition is 
$3.850$. `cmstatr` produces very similar results.

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "CTD") %>%
  levene_test(strength, batch)
expect_equal(res$f, 3.850, tolerance = 0.005)
res
```


### 4.1.6 Nonparametric Basis Values {#c11-hk-wf}
CMH-17-1G Section 8.3.11.2.1 provides STAT17 outputs for the ETW
condition of the present data set. The nonparametric Basis values are listed.
In this case, the Hanson--Koopmans method is used. The published A-Basis
value is 13.0 and the B-Basis is 37.9.


```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW") %>%
  basis_hk_ext(strength, method = "woodward-frawley", p = 0.99, conf = 0.95,
               override = "all")

expect_equal(res$basis, 13.0, tolerance = 0.001)

res
```


```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW") %>%
  basis_hk_ext(strength, method = "optimum-order", p = 0.90, conf = 0.95,
               override = "all")

expect_equal(res$basis, 37.9, tolerance = 0.001)

res
```


### 4.1.7 Single-Point ANOVA Basis Value {#c11-anova}
CMH-17-1G Section 8.3.11.2.2 provides output from STAT17 for
the ETW2 condition from the present data set. STAT17 reports
A- and B-Basis values based on the ANOVA method of 34.6 and
63.2, respectively. The results from `cmstatr` are similar.

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW2") %>%
  basis_anova(strength, batch, override = "number_of_groups",
              p = 0.99, conf = 0.95)
expect_equal(res$basis, 34.6, tolerance = 0.001)
res
```

```{r}
res <- dat_8_3_11_1_1 %>%
  filter(condition == "ETW2") %>%
  basis_anova(strength, batch, override = "number_of_groups")
expect_equal(res$basis, 63.2, tolerance = 0.001)
res
```

## 4.2 Dataset From Section 8.3.11.1.2

[@CMH-17-1G] provides an example data set and results from
ASAP [@ASAP2008]. This example data set is duplicated
below:

```{r}
dat_8_3_11_1_2 <- tribble(
  ~batch, ~strength, ~condition,
  1, 79.04517, "CTD", 1, 103.2006, "RTD", 1, 63.22764, "ETW", 1, 54.09806, "ETW2",
  1, 102.6014, "CTD", 1, 105.1034, "RTD", 1, 70.84454, "ETW", 1, 58.87615, "ETW2",
  1, 97.79372, "CTD", 1, 105.1893, "RTD", 1, 66.43223, "ETW", 1, 61.60167, "ETW2",
  1, 92.86423, "CTD", 1, 100.4189, "RTD", 1, 75.37771, "ETW", 1, 60.23973, "ETW2",
  1, 117.218,  "CTD", 2, 85.32319, "RTD", 1, 72.43773, "ETW", 1, 61.4808,  "ETW2",
  1, 108.7168, "CTD", 2, 92.69923, "RTD", 1, 68.43073, "ETW", 1, 64.55832, "ETW2",
  1, 112.2773, "CTD", 2, 98.45242, "RTD", 1, 69.72524, "ETW", 2, 57.76131, "ETW2",
  1, 114.0129, "CTD", 2, 104.1014, "RTD", 2, 66.20343, "ETW", 2, 49.91463, "ETW2",
  2, 106.8452, "CTD", 2, 91.51841, "RTD", 2, 60.51251, "ETW", 2, 61.49271, "ETW2",
  2, 112.3911, "CTD", 2, 101.3746, "RTD", 2, 65.69334, "ETW", 2, 57.7281,  "ETW2",
  2, 115.5658, "CTD", 2, 101.5828, "RTD", 2, 62.73595, "ETW", 2, 62.11653, "ETW2",
  2, 87.40657, "CTD", 2, 99.57384, "RTD", 2, 59.00798, "ETW", 2, 62.69353, "ETW2",
  2, 102.2785, "CTD", 2, 88.84826, "RTD", 2, 62.37761, "ETW", 3, 61.38523, "ETW2",
  2, 110.6073, "CTD", 3, 92.18703, "RTD", 3, 64.3947,  "ETW", 3, 60.39053, "ETW2",
  3, 105.2762, "CTD", 3, 101.8234, "RTD", 3, 72.8491,  "ETW", 3, 59.17616, "ETW2",
  3, 110.8924, "CTD", 3, 97.68909, "RTD", 3, 66.56226, "ETW", 3, 60.17616, "ETW2",
  3, 108.7638, "CTD", 3, 101.5172, "RTD", 3, 66.56779, "ETW", 3, 46.47396, "ETW2",
  3, 110.9833, "CTD", 3, 100.0481, "RTD", 3, 66.00123, "ETW", 3, 51.16616, "ETW2",
  3, 101.3417, "CTD", 3, 102.0544, "RTD", 3, 59.62108, "ETW",
  3, 100.0251, "CTD",                     3, 60.61167, "ETW",
                                          3, 57.65487, "ETW",
                                          3, 66.51241, "ETW",
                                          3, 64.89347, "ETW",
                                          3, 57.73054, "ETW",
                                          3, 68.94086, "ETW",
                                          3, 61.63177, "ETW"
)
```


### 4.2.1 Pooled SD A- and B-Basis {#c12-psd}
CMH-17-1G Table 8.3.11.2(k) provides outputs from ASAP for the data set
above. ASAP uses the pooled SD method. ASAP produces the following
results, which are quite similar to those produced by `cmstatr`.

Condition | CTD   | RTD   | ETW   | ETW2
----------|-------|-------|-------|------
B-Basis   | 93.64 | 87.30 | 54.33 | 47.12
A-Basis   | 89.19 | 79.86 | 46.84 | 39.69

```{r}

res <- basis_pooled_sd(dat_8_3_11_1_2, strength, condition,
                         override = "all")

expect_equal(res$basis$value[res$basis$group == "CTD"],
           93.64, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
           87.30, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
           54.33, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW2"],
           47.12, tolerance = 0.001)
res
```


```{r}
res <- basis_pooled_sd(dat_8_3_11_1_2, strength, condition,
                       p = 0.99, conf = 0.95,
                       override = "all")
expect_equal(res$basis$value[res$basis$group == "CTD"],
           86.19, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
           79.86, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
           46.84, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW2"],
           39.69, tolerance = 0.001)
res
```

### 4.2.2 Pooled SD A- and B-Basis (Mod CV) {#c12-psdmcv}

After removal of the ETW2 condition, CMH17-STATS reports the pooled A- and 
B-Basis (mod CV) shown in the following table.
`cmstatr` computes very similar values.

Condition | CTD   | RTD   | ETW
----------|-------|-------|------
B-Basis   | 92.25 | 85.91 | 52.97
A-Basis   | 83.81 | 77.48 | 44.47


```{r}
res <- dat_8_3_11_1_2 %>%
  filter(condition != "ETW2") %>%
  basis_pooled_sd(strength, condition, modcv = TRUE, override = "all")
expect_equal(res$basis$value[res$basis$group == "CTD"],
             92.25, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
             85.91, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
             52.97, tolerance = 0.001)
res
```


```{r}
res <- dat_8_3_11_1_2 %>%
  filter(condition != "ETW2") %>%
  basis_pooled_sd(strength, condition,
                  p = 0.99, conf = 0.95, modcv = TRUE, override = "all")
expect_equal(res$basis$value[res$basis$group == "CTD"],
             83.81, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
             77.48, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
             44.47, tolerance = 0.001)
res
```



### 4.2.3 Pooled CV A- and B-Basis {#c12-pcv}
This data set was input into CMH17-STATS and the Pooled CV method was selected.
The results from CMH17-STATS were as follows. `cmstatr` produces very similar
results.

Condition | CTD   | RTD   | ETW   | ETW2
----------|-------|-------|-------|------
B-Basis   | 90.89 | 85.37 | 56.79 | 50.55
A-Basis   | 81.62 | 76.67 | 50.98 | 45.40

```{r}
res <- basis_pooled_cv(dat_8_3_11_1_2, strength, condition, override = "all")
expect_equal(res$basis$value[res$basis$group == "CTD"],
             90.89, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
             85.37, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
             56.79, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW2"],
             50.55, tolerance = 0.001)
res
```

```{r}
res <- basis_pooled_cv(dat_8_3_11_1_2, strength, condition,
                       p = 0.99, conf = 0.95, override = "all")
expect_equal(res$basis$value[res$basis$group == "CTD"],
             81.62, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
             76.67, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
             50.98, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW2"],
             45.40, tolerance = 0.001)
res
```


### 4.2.4 Pooled CV A- and B-Basis (Mod CV) {#c12-pcvmcv}
This data set was input into CMH17-STATS and the Pooled CV method was selected
with the modified CV transform. Additionally, the ETW2 condition was removed.
The results from CMH17-STATS were as follows. `cmstatr` produces very similar
results.

Condition | CTD   | RTD   | ETW   
----------|-------|-------|-------
B-Basis   | 90.31 | 84.83 | 56.43 
A-Basis   | 80.57 | 75.69 | 50.33 


```{r}
res <- dat_8_3_11_1_2 %>%
  filter(condition != "ETW2") %>%
  basis_pooled_cv(strength, condition, modcv = TRUE, override = "all")
expect_equal(res$basis$value[res$basis$group == "CTD"],
             90.31, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
             84.83, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
             56.43, tolerance = 0.001)
res
```

```{r}
res <- dat_8_3_11_1_2 %>%
  filter(condition != "ETW2") %>%
  basis_pooled_cv(strength, condition, modcv = TRUE,
                  p = 0.99, conf = 0.95, override = "all")
expect_equal(res$basis$value[res$basis$group == "CTD"],
             80.57, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "RTD"],
             75.69, tolerance = 0.001)
expect_equal(res$basis$value[res$basis$group == "ETW"],
             50.33, tolerance = 0.001)
res
```


# 5. Other Data Sets
This section contains various small data sets. In most cases, these data sets
were generated randomly for the purpose of comparing `cmstatr` to other
software.

## 5.1 Anderson--Darling Test (Normal) {#ods-adn}
The following data set was randomly generated. When this is
entered into STAT17 [@STAT-17], that software gives the value
$OSL = 0.465$, which matches the result of `cmstatr`
within a small margin.

```{r}
dat <- data.frame(
  strength = c(
    137.4438, 139.5395, 150.89, 141.4474, 141.8203, 151.8821, 143.9245,
    132.9732, 136.6419, 138.1723, 148.7668, 143.283, 143.5429,
    141.7023, 137.4732, 152.338, 144.1589, 128.5218
  )
)
res <- anderson_darling_normal(dat, strength)

expect_equal(res$osl, 0.465, tolerance = 0.001)

res
```

## 5.2 Anderson--Darling Test (Lognormal) {#ods-adl}
The following data set was randomly generated. When this is
entered into STAT17 [@STAT-17], that software gives the value
$OSL = 0.480$, which matches the result of `cmstatr`
within a small margin.

```{r}
dat <- data.frame(
  strength = c(
    137.4438, 139.5395, 150.89, 141.4474, 141.8203, 151.8821, 143.9245,
    132.9732, 136.6419, 138.1723, 148.7668, 143.283, 143.5429,
    141.7023, 137.4732, 152.338, 144.1589, 128.5218
  )
)

res <- anderson_darling_lognormal(dat, strength)

expect_equal(res$osl, 0.480, tolerance = 0.001)

res
```

## 5.3 Anderson--Darling Test (Weibull) {#ods-adw}
The following data set was randomly generated. When this is
entered into STAT17 [@STAT-17], that software gives the value
$OSL = 0.179$, which matches the result of `cmstatr`
within a small margin.

```{r}
dat <- data.frame(
  strength = c(
    137.4438, 139.5395, 150.89, 141.4474, 141.8203, 151.8821, 143.9245,
    132.9732, 136.6419, 138.1723, 148.7668, 143.283, 143.5429,
    141.7023, 137.4732, 152.338, 144.1589, 128.5218
  )
)

res <- anderson_darling_weibull(dat, strength)

expect_equal(res$osl, 0.179, tolerance = 0.002)

res
```

## 5.4 Normal A- and B-Basis {#ods-nb}

The following data was input into STAT17 and the A- and B-Basis values
were computed assuming normally distributed data. The results were 120.336
and 129.287, respectively. `cmstatr` reports very similar values.


```{r}
dat <- c(
  137.4438, 139.5395, 150.8900, 141.4474, 141.8203, 151.8821, 143.9245,
  132.9732, 136.6419, 138.1723, 148.7668, 143.2830, 143.5429, 141.7023,
  137.4732, 152.3380, 144.1589, 128.5218
)

res <- basis_normal(x = dat, p = 0.99, conf = 0.95, override = "all")
expect_equal(res$basis, 120.336, tolerance = 0.0005)
res
```


```{r}
res <- basis_normal(x = dat, p = 0.9, conf = 0.95, override = "all")
expect_equal(res$basis, 129.287, tolerance = 0.0005)
res
```


## 5.5 Lognormal A- and B-Basis {#ods-lb}

The following data was input into STAT17 and the A- and B-Basis values
were computed assuming distributed according to a lognormal distribution.
The results were 121.710
and 129.664, respectively. `cmstatr` reports very similar values.


```{r}
dat <- c(
  137.4438, 139.5395, 150.8900, 141.4474, 141.8203, 151.8821, 143.9245,
  132.9732, 136.6419, 138.1723, 148.7668, 143.2830, 143.5429, 141.7023,
  137.4732, 152.3380, 144.1589, 128.5218
)

res <- basis_lognormal(x = dat, p = 0.99, conf = 0.95, override = "all")
expect_equal(res$basis, 121.710, tolerance = 0.0005)
res
```

```{r}
res <- basis_lognormal(x = dat, p = 0.9, conf = 0.95, override = "all")
expect_equal(res$basis, 129.664, tolerance = 0.0005)
res
```


## 5.6 Weibull A- and B-Basis {#ods-wb}

The following data was input into STAT17 and the A- and B-Basis values
were computed assuming data following the Weibull distribution.
The results were 109.150
and 125.441, respectively. `cmstatr` reports very similar values.


```{r}
dat <- c(
  137.4438, 139.5395, 150.8900, 141.4474, 141.8203, 151.8821, 143.9245,
  132.9732, 136.6419, 138.1723, 148.7668, 143.2830, 143.5429, 141.7023,
  137.4732, 152.3380, 144.1589, 128.5218
)

res <- basis_weibull(x = dat, p = 0.99, conf = 0.95, override = "all")
expect_equal(res$basis, 109.150, tolerance = 0.005)
res
```

```{r}
res <- basis_weibull(x = dat, p = 0.9, conf = 0.95, override = "all")
expect_equal(res$basis, 125.441, tolerance = 0.005)
res
```

## 5.7 Extended Hanson--Koopmans A- and B-Basis {#ods-hkb}

The following data was input into STAT17 and the A- and B-Basis values
were computed using the nonparametric (small sample) method.
The results were 99.651
and 124.156, respectively. `cmstatr` reports very similar values.

```{r}
dat <- c(
  137.4438, 139.5395, 150.8900, 141.4474, 141.8203, 151.8821, 143.9245,
  132.9732, 136.6419, 138.1723, 148.7668, 143.2830, 143.5429, 141.7023,
  137.4732, 152.3380, 144.1589, 128.5218
)

res <- basis_hk_ext(x = dat, p = 0.99, conf = 0.95,
                    method = "woodward-frawley", override = "all")
expect_equal(res$basis, 99.651, tolerance = 0.005)
res
```


```{r}
res <- basis_hk_ext(x = dat, p = 0.9, conf = 0.95,
                    method = "optimum-order", override = "all")
expect_equal(res$basis, 124.156, tolerance = 0.005)
res
```


## 5.8 Extended Hanson--Koopmans B-Basis {#ods-hkb2}
The following random numbers were generated.

```{r}
dat <- c(
  139.6734, 143.0032, 130.4757, 144.8327, 138.7818, 136.7693, 148.636,
  131.0095, 131.4933, 142.8856, 158.0198, 145.2271, 137.5991, 139.8298,
  140.8557, 137.6148, 131.3614, 152.7795, 145.8792, 152.9207, 160.0989,
  145.1920, 128.6383, 141.5992, 122.5297, 159.8209, 151.6720, 159.0156
)
```

All of the numbers above were input into STAT17 and the reported B-Basis
value using the Optimum Order nonparametric method was 122.36798. This
result matches the results of `cmstatr` within a small margin.

```{r}
res <- basis_hk_ext(x = dat, p = 0.9, conf = 0.95,
                    method = "optimum-order", override = "all")
expect_equal(res$basis, 122.36798, tolerance = 0.001)
res
```

The last two observations from the above data set were discarded, leaving
26 observations. This smaller data set was input into STAT17 and that
software calculated a B-Basis value of 121.57073 using the Optimum
Order nonparametric method. `cmstatr` reports a very similar number.

```{r}
res <- basis_hk_ext(x = head(dat, 26), p = 0.9, conf = 0.95,
                    method = "optimum-order", override = "all")
expect_equal(res$basis, 121.57073, tolerance = 0.001)
res
```

The same data set was further reduced such that only the first 22
observations were included. This smaller data set was input into STAT17
and that
software calculated a B-Basis value of 128.82397 using the Optimum
Order nonparametric method. `cmstatr` reports a very similar number.

```{r}
res <- basis_hk_ext(x = head(dat, 22), p = 0.9, conf = 0.95,
                    method = "optimum-order", override = "all")
expect_equal(res$basis, 128.82397, tolerance = 0.001)
res
```

## 5.9 Large Sample Nonparametric B-Basis {#ods-lsnb}

The following data was input into STAT17 and the B-Basis value
was computed using the nonparametric (large sample) method.
The results was 122.738297. `cmstatr` reports very similar values.

```{r}
dat <- c(
  137.3603, 135.6665, 136.6914, 154.7919, 159.2037, 137.3277, 128.821,
  138.6304, 138.9004, 147.4598, 148.6622, 144.4948, 131.0851, 149.0203,
  131.8232, 146.4471, 123.8124, 126.3105, 140.7609, 134.4875, 128.7508,
  117.1854, 129.3088, 141.6789, 138.4073, 136.0295, 128.4164, 141.7733,
  134.455,  122.7383, 136.9171, 136.9232, 138.8402, 152.8294, 135.0633,
  121.052,  131.035,  138.3248, 131.1379, 147.3771, 130.0681, 132.7467,
  137.1444, 141.662,  146.9363, 160.7448, 138.5511, 129.1628, 140.2939,
  144.8167, 156.5918, 132.0099, 129.3551, 136.6066, 134.5095, 128.2081,
  144.0896, 141.8029, 130.0149, 140.8813, 137.7864
)

res <- basis_nonpara_large_sample(x = dat, p = 0.9, conf = 0.95,
                                  override = "all")
expect_equal(res$basis, 122.738297, tolerance = 0.005)
res
```

## 5.10 Acceptance Limits Based on Mean and Extremum {#ods-eme}

Results from `cmstatr`'s `equiv_mean_extremum` function were compared with
results from HYTEQ. The summary statistics for the qualification data
were set as `mean = 141.310` and `sd=6.415`. For a value of `alpha=0.05` and
`n = 9`,
HYTEQ reported thresholds of 123.725 and 137.197 for minimum individual
and mean, respectively. `cmstatr` produces very similar results.

```{r}
res <- equiv_mean_extremum(alpha = 0.05, mean_qual = 141.310, sd_qual = 6.415,
                           n_sample = 9)
expect_equal(res$threshold_min_indiv, 123.725, tolerance = 0.001)
expect_equal(res$threshold_mean, 137.197, tolerance = 0.001)
res
```

Using the same parameters, but using the modified CV method,
HYTEQ produces thresholds of 117.024 and 135.630 for minimum individual
and mean, respectively. `cmstatr` produces very similar results.

```{r}
res <- equiv_mean_extremum(alpha = 0.05, mean_qual = 141.310, sd_qual = 6.415,
                           n_sample = 9, modcv = TRUE)
expect_equal(res$threshold_min_indiv, 117.024, tolerance = 0.001)
expect_equal(res$threshold_mean, 135.630, tolerance = 0.001)
res
```


## 5.11 Acceptance Based on Change in Mean {#ods-ecm}

Results from `cmstatr`'s `equiv_change_mean` function were compared with
results from HYTEQ. The following parameters were used. A value of
`alpha = 0.05` was selected.

Parameter | Qualification | Sample
----------|---------------|---------
Mean      | 9.24          | 9.02
SD        | 0.162         | 0.15785
n         | 28            | 9

HYTEQ gives an acceptance range of 9.115 to 9.365. `cmstatr` produces
similar results.

```{r}
res <- equiv_change_mean(alpha = 0.05, n_sample = 9, mean_sample = 9.02,
                         sd_sample = 0.15785, n_qual = 28, mean_qual = 9.24,
                         sd_qual = 0.162)
expect_equal(res$threshold, c(9.115, 9.365), tolerance = 0.001)
res
```

After selecting the modified CV method, HYTEQ gives an acceptance
range of 8.857 to 9.623. `cmstatr` produces similar results.

```{r}
res <- equiv_change_mean(alpha = 0.05, n_sample = 9, mean_sample = 9.02,
                         sd_sample = 0.15785, n_qual = 28, mean_qual = 9.24,
                         sd_qual = 0.162, modcv = TRUE)
expect_equal(res$threshold, c(8.857, 9.623), tolerance = 0.001)
res
```

# 6. Comparison with Literature
In this section, results from `cmstatr` are compared with values published
in literature.

## 6.1 Anderson--Darling K--Sample Test {#cl-adk}
[@Stephens1987] provides example data that compares measurements obtained
in four labs. Their paper gives values of the ADK test statistic as well as
p-values.

The data in [@Stephens1987] is as follows:

```{r}
dat_ss1987 <- data.frame(
  smoothness = c(
    38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0,
    39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8,
    34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0,
    34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8
  ),
  lab = c(rep("A", 8), rep("B", 8), rep("C", 8), rep("D", 8))
)
dat_ss1987
```

[@Stephens1987] lists the corresponding test statistics
$A_{akN}^2 = 8.3926$ and $\sigma_N = 1.2038$ with the p-value
$p = 0.0022$. These match the result of `cmstatr` within a small margin.

```{r}
res <- ad_ksample(dat_ss1987, smoothness, lab)

expect_equal(res$ad, 8.3926, tolerance = 0.001)
expect_equal(res$sigma, 1.2038, tolerance = 0.001)
expect_equal(res$p, 0.00226, tolerance = 0.01)

res
```


# 7. Comparison with Published Factors
Various factors, such as tolerance limit factors, are published in various
publications. This section compares those published factors with those
computed by `cmstatr`.

## 7.1 Normal kB Factors {#pf-kb}
B-Basis tolerance limit factors assuming a normal distribution are published in
CMH-17-1G. Those factors are reproduced below and are compared with the
results of `cmstatr`. The published factors and those computed by `cmstatr`
are quite similar.

```{r}
tribble(
  ~n, ~kB_published,
  2, 20.581, 36, 1.725, 70, 1.582, 104, 1.522,
  3, 6.157, 37, 1.718, 71, 1.579, 105, 1.521,
  4, 4.163, 38, 1.711, 72, 1.577, 106, 1.519,
  5, 3.408, 39, 1.704, 73, 1.575, 107, 1.518,
  6, 3.007, 40, 1.698, 74, 1.572, 108, 1.517,
  7, 2.756, 41, 1.692, 75, 1.570, 109, 1.516,
  8, 2.583, 42, 1.686, 76, 1.568, 110, 1.515,
  9, 2.454, 43, 1.680, 77, 1.566, 111, 1.513,
  10, 2.355, 44, 1.675, 78, 1.564, 112, 1.512,
  11, 2.276, 45, 1.669, 79, 1.562, 113, 1.511,
  12, 2.211, 46, 1.664, 80, 1.560, 114, 1.510,
  13, 2.156, 47, 1.660, 81, 1.558, 115, 1.509,
  14, 2.109, 48, 1.655, 82, 1.556, 116, 1.508,
  15, 2.069, 49, 1.650, 83, 1.554, 117, 1.507,
  16, 2.034, 50, 1.646, 84, 1.552, 118, 1.506,
  17, 2.002, 51, 1.642, 85, 1.551, 119, 1.505,
  18, 1.974, 52, 1.638, 86, 1.549, 120, 1.504,
  19, 1.949, 53, 1.634, 87, 1.547, 121, 1.503,
  20, 1.927, 54, 1.630, 88, 1.545, 122, 1.502,
  21, 1.906, 55, 1.626, 89, 1.544, 123, 1.501,
  22, 1.887, 56, 1.623, 90, 1.542, 124, 1.500,
  23, 1.870, 57, 1.619, 91, 1.540, 125, 1.499,
  24, 1.854, 58, 1.616, 92, 1.539, 126, 1.498,
  25, 1.839, 59, 1.613, 93, 1.537, 127, 1.497,
  26, 1.825, 60, 1.609, 94, 1.536, 128, 1.496,
  27, 1.812, 61, 1.606, 95, 1.534, 129, 1.495,
  28, 1.800, 62, 1.603, 96, 1.533, 130, 1.494,
  29, 1.789, 63, 1.600, 97, 1.531, 131, 1.493,
  30, 1.778, 64, 1.597, 98, 1.530, 132, 1.492,
  31, 1.768, 65, 1.595, 99, 1.529, 133, 1.492,
  32, 1.758, 66, 1.592, 100, 1.527, 134, 1.491,
  33, 1.749, 67, 1.589, 101, 1.526, 135, 1.490,
  34, 1.741, 68, 1.587, 102, 1.525, 136, 1.489,
  35, 1.733, 69, 1.584, 103, 1.523, 137, 1.488
) %>%
  arrange(n) %>%
  mutate(kB_cmstatr = k_factor_normal(n, p = 0.9, conf = 0.95)) %>%
  rowwise() %>%
  mutate(diff = expect_equal(kB_published, kB_cmstatr, tolerance = 0.001)) %>%
  select(-c(diff))
```

## 7.2 Normal kA Factors {#pf-ka}
A-Basis tolerance limit factors assuming a normal distribution are published in
CMH-17-1G. Those factors are reproduced below and are compared with the
results of `cmstatr`. The published factors and those computed by `cmstatr`
are quite similar.

```{r}
tribble(
  ~n, ~kA_published,
  2, 37.094, 36, 2.983, 70, 2.765, 104, 2.676,
  3, 10.553, 37, 2.972, 71, 2.762, 105, 2.674,
  4, 7.042, 38, 2.961, 72, 2.758, 106, 2.672,
  5, 5.741, 39, 2.951, 73, 2.755, 107, 2.671,
  6, 5.062, 40, 2.941, 74, 2.751, 108, 2.669,
  7, 4.642, 41, 2.932, 75, 2.748, 109, 2.667,
  8, 4.354, 42, 2.923, 76, 2.745, 110, 2.665,
  9, 4.143, 43, 2.914, 77, 2.742, 111, 2.663,
  10, 3.981, 44, 2.906, 78, 2.739, 112, 2.662,
  11, 3.852, 45, 2.898, 79, 2.736, 113, 2.660,
  12, 3.747, 46, 2.890, 80, 2.733, 114, 2.658,
  13, 3.659, 47, 2.883, 81, 2.730, 115, 2.657,
  14, 3.585, 48, 2.876, 82, 2.727, 116, 2.655,
  15, 3.520, 49, 2.869, 83, 2.724, 117, 2.654,
  16, 3.464, 50, 2.862, 84, 2.721, 118, 2.652,
  17, 3.414, 51, 2.856, 85, 2.719, 119, 2.651,
  18, 3.370, 52, 2.850, 86, 2.716, 120, 2.649,
  19, 3.331, 53, 2.844, 87, 2.714, 121, 2.648,
  20, 3.295, 54, 2.838, 88, 2.711, 122, 2.646,
  21, 3.263, 55, 2.833, 89, 2.709, 123, 2.645,
  22, 3.233, 56, 2.827, 90, 2.706, 124, 2.643,
  23, 3.206, 57, 2.822, 91, 2.704, 125, 2.642,
  24, 3.181, 58, 2.817, 92, 2.701, 126, 2.640,
  25, 3.158, 59, 2.812, 93, 2.699, 127, 2.639,
  26, 3.136, 60, 2.807, 94, 2.697, 128, 2.638,
  27, 3.116, 61, 2.802, 95, 2.695, 129, 2.636,
  28, 3.098, 62, 2.798, 96, 2.692, 130, 2.635,
  29, 3.080, 63, 2.793, 97, 2.690, 131, 2.634,
  30, 3.064, 64, 2.789, 98, 2.688, 132, 2.632,
  31, 3.048, 65, 2.785, 99, 2.686, 133, 2.631,
  32, 3.034, 66, 2.781, 100, 2.684, 134, 2.630,
  33, 3.020, 67, 2.777, 101, 2.682, 135, 2.628,
  34, 3.007, 68, 2.773, 102, 2.680, 136, 2.627,
  35, 2.995, 69, 2.769, 103, 2.678, 137, 2.626
) %>%
  arrange(n) %>%
  mutate(kA_cmstatr = k_factor_normal(n, p = 0.99, conf = 0.95)) %>%
  rowwise() %>%
  mutate(diff = expect_equal(kA_published, kA_cmstatr, tolerance = 0.001)) %>%
  select(-c(diff))
```

## 7.3 Nonparametric B-Basis Extended Hanson--Koopmans {#pf-hk}
Vangel [@Vangel1994] provides extensive tables of $z$ for the case where
$i=1$ and $j$ is the median observation. This section checks the results of
`cmstatr`'s function against those tables. Only the odd values of $n$
are checked so that the median is a single observation. The unit tests for
the `cmstatr` package include checks of a variety of values of $p$ and
confidence, but only the factors for B-Basis are checked here.

```{r}
tribble(
  ~n, ~z,
  3,  28.820048,
  5,  6.1981307,
  7,  3.4780112,
  9,  2.5168762,
  11, 2.0312134,
  13, 1.7377374,
  15, 1.5403989,
  17, 1.3979806,
  19, 1.2899172,
  21, 1.2048089,
  23, 1.1358259,
  25, 1.0786237,
  27, 1.0303046,
) %>%
  rowwise() %>%
  mutate(
    z_calc = hk_ext_z(n, 1, ceiling(n / 2), p = 0.90, conf = 0.95)
  ) %>%
  mutate(diff = expect_equal(z, z_calc, tolerance = 0.0001)) %>% 
  select(-c(diff))
```


## 7.4 Nonparametric A-Basis Extended Hanson--Koopmans {#pf-hk2}
CMH-17-1G provides Table 8.5.15, which contains factors for calculating
A-Basis values using the Extended Hanson--Koopmans nonparametric method.
That table is reproduced in part here and the factors are compared with
those computed by `cmstatr`. More extensive checks are performed in the
unit test of the `cmstatr` package. The factors computed by `cmstatr` are
very similar to those published in CMH-17-1G.

```{r}
tribble(
  ~n, ~k,
  2, 80.0038,
  4, 9.49579,
  6, 5.57681,
  8, 4.25011,
  10, 3.57267,
  12, 3.1554,
  14, 2.86924,
  16, 2.65889,
  18, 2.4966,
  20, 2.36683,
  25, 2.131,
  30, 1.96975,
  35, 1.85088,
  40, 1.75868,
  45, 1.68449,
  50, 1.62313,
  60, 1.5267,
  70, 1.45352,
  80, 1.39549,
  90, 1.34796,
  100, 1.30806,
  120, 1.24425,
  140, 1.19491,
  160, 1.15519,
  180, 1.12226,
  200, 1.09434,
  225, 1.06471,
  250, 1.03952,
  275, 1.01773
) %>%
  rowwise() %>%
  mutate(z_calc = hk_ext_z(n, 1, n, 0.99, 0.95)) %>%
  mutate(diff = expect_lt(abs(k - z_calc), 0.0001)) %>% 
  select(-c(diff))
```


## 7.5 Factors for Small Sample Nonparametric B-Basis {#pf-hk-opt}
CMH-17-1G Table 8.5.14 provides ranks orders and factors for computing
nonparametric B-Basis values. This table is reproduced below and
compared with the results of `cmstatr`. The results are similar. In some
cases, the rank order ($r$ in CMH-17-1G or $j$ in `cmstatr`) *and* the 
the factor ($k$) are different. These differences are discussed in detail in
the vignette
[Extended Hanson-Koopmans](hk_ext.html).

```{r}
tribble(
  ~n, ~r, ~k,
  2, 2, 35.177,
  3, 3, 7.859,
  4, 4, 4.505,
  5, 4, 4.101,
  6, 5, 3.064,
  7, 5, 2.858,
  8, 6, 2.382,
  9, 6, 2.253,
  10, 6, 2.137,
  11, 7, 1.897,
  12, 7, 1.814,
  13, 7, 1.738,
  14, 8, 1.599,
  15, 8, 1.540,
  16, 8, 1.485,
  17, 8, 1.434,
  18, 9, 1.354,
  19, 9, 1.311,
  20, 10, 1.253,
  21, 10, 1.218,
  22, 10, 1.184,
  23, 11, 1.143,
  24, 11, 1.114,
  25, 11, 1.087,
  26, 11, 1.060,
  27, 11, 1.035,
  28, 12, 1.010
) %>% 
  rowwise() %>%
  mutate(r_calc = hk_ext_z_j_opt(n, 0.90, 0.95)$j) %>%
  mutate(k_calc = hk_ext_z_j_opt(n, 0.90, 0.95)$z)
```



## 7.6 Nonparametric B-Basis Binomial Rank {#pf-npbinom}
CMH-17-1G Table 8.5.12 provides factors for computing B-Basis values
using the nonparametric binomial rank method. Part of that table is
reproduced below and compared with the results of `cmstatr`.
The results of `cmstatr` are similar to the published values.
A more complete comparison is performed in the units tests of the
`cmstatr` package.

```{r}
tribble(
  ~n, ~rb,
  29, 1,
  46, 2,
  61, 3,
  76, 4,
  89, 5,
  103, 6,
  116, 7,
  129, 8,
  142, 9,
  154, 10,
  167, 11,
  179, 12,
  191, 13,
  203, 14
) %>%
  rowwise() %>%
  mutate(r_calc = nonpara_binomial_rank(n, 0.9, 0.95)) %>%
  mutate(test = expect_equal(rb, r_calc)) %>%
  select(-c(test))
```

## 7.7 Nonparametric A-Basis Binomial Rank {#pf-npbinom2}
CMH-17-1G Table 8.5.13 provides factors for computing B-Basis values
using the nonparametric binomial rank method. Part of that table is
reproduced below and compared with the results of `cmstatr`.
The results of `cmstatr` are similar to the published values.
A more complete comparison is performed in the units tests of the
`cmstatr` package.

```{r}
tribble(
  ~n, ~ra,
  299, 1,
  473, 2,
  628, 3,
  773, 4,
  913, 5
) %>%
  rowwise() %>%
  mutate(r_calc = nonpara_binomial_rank(n, 0.99, 0.95)) %>%
  mutate(test = expect_equal(ra, r_calc)) %>%
  select(-c(test))
```


## 7.8 Factors for Equivalency {#pf-equiv}

Vangel's 2002 paper provides factors for calculating limits for sample
mean and sample extremum for various values of $\alpha$ and sample size ($n$).
A subset of those factors are reproduced below and compared with results
from `cmstatr`. The results are very similar for values of $\alpha$ and $n$
that are common for composite materials.


```{r}
read.csv(system.file("extdata", "k1.vangel.csv",
                     package = "cmstatr")) %>%
  gather(n, k1, X2:X10) %>%
  mutate(n = as.numeric(substring(n, 2))) %>%
  inner_join(
    read.csv(system.file("extdata", "k2.vangel.csv",
                         package = "cmstatr")) %>%
      gather(n, k2, X2:X10) %>%
      mutate(n = as.numeric(substring(n, 2))),
    by = c("n" = "n", "alpha" = "alpha")
  ) %>%
  filter(n >= 5 & (alpha == 0.01 | alpha == 0.05)) %>% 
  group_by(n, alpha) %>%
  nest() %>%
  mutate(equiv = map2(alpha, n, ~k_equiv(.x, .y))) %>%
  mutate(k1_calc = map(equiv, function(e) e[1]),
         k2_calc = map(equiv, function(e) e[2])) %>%
  select(-c(equiv)) %>% 
  unnest(cols = c(data, k1_calc, k2_calc)) %>% 
  mutate(check = expect_equal(k1, k1_calc, tolerance = 0.0001)) %>%
  select(-c(check)) %>% 
  mutate(check = expect_equal(k2, k2_calc, tolerance = 0.0001)) %>%
  select(-c(check))
```


# 8. Session Info
This copy of this vignette was build on the following system.

```{r}
sessionInfo()
```

# 9. References

---
title: "cmstatr Tutorial"
author: "Stefan Kloppenborg"
date: "1-Apr-2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{cmstatr Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# If any of the required packages are unavailable,
# don't re-run the code
required <- c("dplyr", "ggplot2", "tidyr", "cmstatr", "purrr")
if (!all(unlist(lapply(required, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)}
  )))) {
  knitr::opts_chunk$set(eval = FALSE)
}
```

`cmstatr` is an R package for analyzing composite material data for use in the
aerospace industry. The statistical methods are based on those published in
[CMH-17-1G](https://www.cmh17.org/). This package is intended to facilitate
reproducible statistical analysis of composite materials. In this tutorial,
we'll explore the basic functionality of `cmstatr`.

Before we can actually use the package, we'll need to load it. We'll also load
the `dplyr` package, which we'll talk about shortly. There are also a few other
packages that we'll load. These could all be loaded by loading the
`tidyverse` package instead.

```{r message=FALSE}
library(cmstatr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
```

# Input Data
`cmstatr` is built with the assumption that the data is in (so called)
[tidy data](http://vita.had.co.nz/papers/tidy-data.html) format. This means
that the data is in a data frame and that each observation (i.e. test result)
has its own row and that each variable has its own column. Included in this
package is a sample composite material data set (this data set is fictional:
don't use it for anything other than learning this package). The data set
`carbon.fabric.2` has the expected format. We'll just show the first 10 rows
of the data for now.

```{r}
carbon.fabric.2 %>%
  head(10)
```

If your data set is not yet in this type of format (note: that the column
names *do not* need to match the column names in the example), there are
many ways to get it into this format. One of the easier ways of doing so
is to use the [`tidyr`](https://tidyr.tidyverse.org/) package. The use of this
package is outside the scope of this vignette.

# Working With Data
Throughout this vignette, we will be using some of the `tidyverse` tools for
working with data. There are several ways to work with data in R, but in the
opinion of the author of this vignette, the `tidyverse` provides the easiest
way to do so. As such, this is the approach used in this vignette. Feel free
to use whichever approach works best for you.

# Normalizing Data to Cured Ply Thickness
Very often, you'll want to normalize as-measured strength data to a nominal
cured ply thickness for fiber-dominated properties. Very often, this will
reduce the apparent variance in the data. The `normalize_ply_thickness`
function can be used to normalize strength or modulus data to a certain
cured ply thickness. This function takes three arguments: the value to
normalize (i.e.. strength or modulus), the measured thickness and the
nominal thickness. In our case, the nominal cured ply thickness of the
material is $0.0079$. We can then normalize the warp-tension and
fill-compression data as follows:

```{r}
norm_data <- carbon.fabric.2 %>%
  filter(test == "WT" | test == "FC") %>%
  mutate(strength.norm = normalize_ply_thickness(strength,
                                                 thickness / nplies,
                                                 0.0079))

norm_data %>%
  head(10)
```

# Calculating Single-Point Basis Value
The simplest thing that you will likely do is to calculate a basis value based
of a set of numbers that you consider as unstructured data. An example of this
would be calculating the B-Basis of the `RTD` warp tension (`WT`) data.

There are a number of diagnostic tests that we should run before actually
calculating a B-Basis value. We'll talk about those later, but for now, let's
just get right to checking how the data are distributed and calculating the
B-Basis.

We'll use an Anderson--Darling test to check if the data are normally
distributed. The `cmstatr` package provides the function
`anderson_darling_normal` and related functions for other distributions.
We can run an Anderson--Darling test for normality on the warp tension RTD
data as follows. We'll perform this test on the normalized strength.

```{r}
norm_data %>%
  filter(test == "WT" & condition == "RTD") %>%
  anderson_darling_normal(strength.norm)
```

```{r include=FALSE}
# Verify that the AD test always provides the same conclusion
# If this assertion fails, the Vignette needs to be re-written
if (0.05 >= (norm_data %>%
  filter(test == "WT" & condition == "RTD") %>%
  anderson_darling_normal(strength.norm))$osl) {
  stop("Unexpected vale for Anderson-Darling test")
  }
```

Now that we know that this data follows a normal distribution (since the
observed significance level (OSL) of the Anderson--Darling test is
greater than $0.05$), we can
proceed to calculate a basis value based based on the assumption of normally
distributed data. The `cmstatr` package provides the function `basis_normal`
as well as related functions for other distributions. By default, the B-Basis
value is calculated, but other population proportions and confidence bounds
can be specified (for example, specify `p = 0.99, conf = 0.99` for A-Basis).

```{r}
norm_data %>%
  filter(test == "WT" & condition == "RTD") %>%
  basis_normal(strength.norm)
```

We see that the calculated B-Basis is $129.96$. We also see two messages
issued by the `cmstatr` package. These messages relate to the automated
diagnostic tests performed by the basis calculation functions. In this case
we see messages that two of the diagnostic tests were not performed because
we didn't specify the batch of each observation. The batch is not required
for calculating single-point basis values, but it is required for performing
batch-to-batch variability and within-batch outlier diagnostic tests.

The `basis_normal` function performs the following diagnostic tests by default:

- Within batch outliers using `maximum_normed_residual()`
- Between batch variability using `ad_ksample()`
- Outliers using `maximum_normed_residual()`
- Normality of data using `anderson_darling_normal()`

There are two ways that we can deal with the two messages that we see. We can
pass in a column that specifies the batch for each observation, or we can
override those two diagnostic tests so that `cmstatr` doesn't run them.

To override the two diagnostic tests, we set the argument `override` to a list
of the names of the diagnostic tests that we want to skip. The names of the
diagnostic tests that were not run are shown between back-ticks (\`) in the
message. Our call to `basis_normal()` would be updated as follows:

```{r}
norm_data %>%
  filter(test == "WT" & condition == "RTD") %>%
  basis_normal(strength.norm, 
               override = c("outliers_within_batch",
                            "between_batch_variability"))
```

Obviously, you should be cautious about overriding the diagnostic tests.
There are certainly times when it is appropriate to do so, but sound
engineering judgment is required.

The better approach would be to specify the batch. This can be done as
follows:

```{r}
norm_data %>%
  filter(test == "WT" & condition == "RTD") %>%
  basis_normal(strength.norm, batch)
```

Now that batch is specified, we see that one of the diagnostic tests
actually fails: the Anderson--Darling k-Sample test shows that the batches
are not drawn from the same (unspecified) distribution. We can run this
diagnostic test directly to investigate further:

```{r}
norm_data %>%
  filter(test == "WT" & condition == "RTD") %>%
  ad_ksample(strength.norm, batch)
```

For the Anderson--Darling k-Sample test, $\alpha=0.025$ is normally used.
In this case the p-value is $p=0.0026$, so it is no where near $\alpha$
(note the number of decimal places).

We can plot the distribution of this data and make a judgment call about
whether to continue.


```{r}
norm_data %>%
  filter(test == "WT" & condition == "RTD") %>%
  group_by(batch) %>%
  ggplot(aes(x = strength.norm, color = batch)) +
  stat_normal_surv_func() +
  stat_esf() +
  ggtitle("Distribution of Data For Each Batch")
```

We can also run the other diagnostic test by themselves. These are described
in more detail in the following sections.

# Calculating Basis Values by Pooling Across Environments
In this section, we'll use the fill-compression data from the `carbon.fabric.2`
data set.

## Checking for Outliers
After checking that there are a sufficient number of conditions, batches and
specimens and that the failure modes are consistent, we would normally
check if there are outliers within each batch and condition. The maximum
normed residual test can be used for this. The `cmstatr` package provides the
function `maximum_normed_residual` to do this. First, we'll group the data
by condition and batch, then run the test on each group. The
`maximum_normed_residual` function returns an object that contains a number
of values. We'll create a `data.frame` that contains those values.

In order to do this, we need to use the `nest` function from the `tidyr`
package. This is explained in detail
[here](https://tidyr.tidyverse.org/articles/nest.html). Basically,
`nest` allows a column of `list`s or a column of `data.frame`s to be
added to a `data.frame`. Once nested, we can use the `glance` method
to unpack the values returned by `maximum_normed_residual` into a
one-row `data.frame`, and then use `unnest` to flatten this into a
single `data.frame`.


```{r}
norm_data %>%
  filter(test == "FC") %>%
  group_by(condition, batch) %>%
  nest() %>%
  mutate(mnr = map(data,
                   ~maximum_normed_residual(data = .x, x = strength.norm)),
         tidied = map(mnr, glance)) %>%
  select(-c(mnr, data)) %>%  # remove unneeded columns
  unnest(tidied)
```

```{r include=FALSE}
if ((norm_data %>%
  filter(test == "FC") %>%
  group_by(condition, batch) %>%
  summarise(
    n_outliers = maximum_normed_residual(x = strength.norm)$n_outliers
    ) %>%
  ungroup() %>%
  summarise(n_outliers = sum(n_outliers)))[[1]] != 0) {
  stop("Unexpected number of outliers")
  }
```

None of the groups have outliers, so we can continue.

# Batch-to-Batch Distribution
Next, we will use the Anderson--Darling k-Sample test to check that each batch
comes from the same distribution within each condition. We can use the
`ad_ksample` function from `cmstatr` to do so. Once again, we'll use
`nest`/`unnest` and `glance` to do so.

```{r}
norm_data %>%
  filter(test == "FC") %>%
  group_by(condition) %>%
  nest() %>%
  mutate(adk = map(data, ~ad_ksample(data = .x,
                                     x = strength.norm,
                                     groups = batch)),
         tidied = map(adk, glance)) %>%
  select(-c(data, adk)) %>%  # remove unneeded columns
  unnest(tidied)
```

```{r include=FALSE}
if (!all(!(norm_data %>%
  filter(test == "FC") %>%
  group_by(condition) %>%
  summarise(different_dist =
           ad_ksample(x = strength.norm, groups = batch)$reject_same_dist
  ))$different_dist)) {
  stop("Unexpected ADK result")
  }
```

For all conditions, the Anderson--Darling k-Sample test fails to reject the
hypothesis that each batch comes from the same (unspecified) distribution.
We can thus proceed to pooling the data.

## Checking for Outliers Within Each Condition
Just as we did when checking for outlier within each condition and each
batch, we can pool all the batches (within each condition) and check
for outliers within each condition.

```{r}
norm_data %>%
  filter(test == "FC") %>%
  group_by(condition) %>%
  nest() %>%
  mutate(mnr = map(data, ~maximum_normed_residual(data = .x,
                                                  x = strength.norm)),
         tidied = map(mnr, glance)) %>%
  select(-c(mnr, data)) %>%  # remove unneeded columns
  unnest(tidied)
```

```{r include=FALSE}
if ((norm_data %>%
  filter(test == "FC") %>%
  group_by(condition) %>%
  summarise(
    n_outliers = maximum_normed_residual(x = strength.norm)$n_outliers
    ) %>%
  ungroup() %>%
  summarise(n_outliers = sum(n_outliers)))[[1]] != 0) {
  stop("Unexpected number of outliers")
  }
```

We find no outliers, so we can continue.

## Pooling Across Environments
Often it is desirable to pool data across several environments. There
are two methods for doing so: "pooled standard deviation" and
"pooled CV" (CV is an abbreviation for Coefficient of Variation)

First, we will check for
equality of variance among the conditions. We will do so using Levene's test.
The `cmstatr` package provides the function `levene_test` to do so.

```{r}
norm_data %>%
  filter(test == "FC") %>%
  levene_test(strength.norm, condition)
```

```{r include=FALSE}
if (!(norm_data %>%
  filter(test == "FC") %>%
  levene_test(strength.norm, condition))$reject_equal_variance) {
  stop("Unexpected result from Levene's test")
  }
```

The result from Levene's test indicates that the variance for each condition
is not equal. This indicates that the data cannot be pooled using the
"pooled standard deviation" method.

We can check if the data can be pooled using the "pooled CV" method.
We'll start by normalizing the data from each group to the group's mean.
The `cmstatr` package provides the function `normalize_group_mean` for
this purpose.

```{r}
norm_data %>%
  filter(test == "FC") %>%
  mutate(
    strength_norm_group = normalize_group_mean(strength.norm, condition)) %>%
  levene_test(strength_norm_group, condition)
```

```{r include=FALSE}
if ((norm_data %>%
  filter(test == "FC") %>%
  mutate(
    strength_norm_group = normalize_group_mean(strength.norm, condition)) %>%
  levene_test(strength_norm_group, condition))$reject_equal_variance) {
  stop("Unexpected value from Levene's test")
  }
```

The Levene's test thus shows the variances of the pooled data are equal.
We can move on to performing an Anderson--Darling test for normality on
the pooled data.

```{r}
norm_data %>%
  filter(test == "FC") %>%
  mutate(
    strength_norm_group = normalize_group_mean(strength.norm, condition)) %>%
  anderson_darling_normal(strength_norm_group)
```

```{r include=FALSE}
if ((norm_data %>%
  filter(test == "FC") %>%
  mutate(
    strength_norm_group = normalize_group_mean(strength.norm, condition)) %>%
  anderson_darling_normal(strength_norm_group))$osl <= 0.05) {
  stop("Unexpected value from AD test")
  }
```

The Anderson--Darling test indicates that the pooled data is drawn from
a normal distribution, so we can continue with calculating basis values
using the "pooled CV" method.

```{r}
norm_data %>%
  filter(test == "FC") %>%
  basis_pooled_cv(strength.norm, condition, batch)
```

The conditions listed in the output above are in alphabetical order. This
probably isn't what you want. Instead, you probably want the conditions
listed in a certain order. This can be done by ordering the data first as
demonstrated below. You're probably just do this one in at the start of
your analysis.

```{r}
norm_data %>%
  mutate(condition = ordered(condition,
                             c("CTD", "RTD", "ETD", "ETW", "ETW2"))) %>%
  filter(test == "FC") %>%
  basis_pooled_cv(strength.norm, condition, batch)
```

# Equivalency
Eventually, once you've finished calculating all your basis values,
you'll probably want to set specification requirements or evaluate
site/process equivalency. `cmstatr` has functionality to do both.

Let's say that you want to develop specification limits for fill compression
that you're going to put in your material specification. You can do this
as follows:

```{r}
carbon.fabric.2 %>%
  filter(test == "FC" & condition == "RTD") %>%
  equiv_mean_extremum(strength, n_sample = 5, alpha = 0.01)
```

If you're determining equivalency limits for modulus, a different
approach is generally used so that bilateral limits are set. `cmstatr`
can do this as well, using the function `equiv_change_mean`.
---
title: "Plotting Composite Material Data"
author: "Ally Fraser"
date: "2-May-2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Plotting Composite Material Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6
)

# If any of the required packages are unavailable,
# don't re-run the code
required <- c("dplyr", "ggplot2", "tidyr", "cmstatr")
if (!all(unlist(lapply(required, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)}
  )))) {
  knitr::opts_chunk$set(eval = FALSE)
}
```

This vignette demonstrates how to create some of the graphs commonly used
when analyzing composite material data. Here, we rely on the
[`ggplot2`](https://ggplot2.tidyverse.org/) package for graphing. This package
can be loaded either on its own, or through the `tidyverse` meta-package, which
also includes packages such as `dplyr` that we will also use.

We'll need to load a few packages in order to proceed.

```{r message=FALSE}
library(dplyr)
library(ggplot2)
library(tidyr)
library(cmstatr)
```

Throughout this vignette, we'll use one of the example data sets that
comes with `cmstatr` and we'll focus on the warp-tension data as an
example. We'll load the example data in a variable as follows. By default
the condition will be in an arbitrary order, but throughout the visualization,
we'll want the conditions shown in a particular order (from coldest and driest
to hottest and wettest). We can define the order of the conditions using
the `ordered` function. For brevity, only the first few rows of the data set
are displayed below.

```{r}
dat <- carbon.fabric.2 %>%
  filter(test == "WT") %>%
  mutate(condition = ordered(condition, c("CTD", "RTD", "ETW", "ETW2")))

dat %>%
  head(10)
```

We'll then calculate the B-Basis value using the pooling by standard deviation
method. This data set happens to fail some of the diagnostic tests, but for the
purpose of this example, we'll ignore those failures using the
`override` argument.

```{r}
b_basis_pooled <- dat %>%
  basis_pooled_cv(strength, condition, batch,
                  override = c("between_group_variability",
                               "normalized_variance_equal"))

b_basis_pooled
```

The object returned from `basis_pooled_cv` contains a number of values.
One value is a `data.frame` containing the groups (i.e. conditions)
and the corresponding basis values. This looks like the following. We'll
use this in the visualizations.

```{r}
b_basis_pooled$basis
```


# Batch Plots
Batch plots are used to identify differences between batches.
Simple batch plots can be created using box plots and adding
horizontal lines for the basis values as follows. Note that
the heavy line in the box of the box plot is the *median*, not
the mean. The two hinges correspond with the first and third
quantiles and the whiskers extend to the most extreme data point,
or 1.5 times the inner quantile range.

In the code below, we use the function `rename` to rename the
column `group` to `condition`. The `data.frame` produced by
`basis_pooled_cv` uses the columns `value` and `group`, but
to match the data, we need the column with the conditions
to be named `condition`.

```{r}
dat %>%
  ggplot(aes(x = batch, y = strength)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  geom_hline(aes(yintercept = value),
             data = b_basis_pooled$basis %>% rename(condition = group),
             color = "blue") +
  facet_grid(. ~ condition) +
  theme_bw() +
  ggtitle("Batch Plot")
```

# Quantile Plots
A quantile plot provides a graphical summary of sample values. 
This plot displays the sample values and the corresponding quantile.
A quantile  plot can  be used to examine the symmetry and tail sizes of the
underlying distribution. Sharp rises may indicate the presence of outliers. 

```{r}
dat %>%
  ggplot(aes(x = strength, color = condition)) +
  stat_ecdf(geom = "point") +
  coord_flip() +
  theme_bw() +
  ggtitle("Quantile Plot")
```

# Normal Survival Function Plots
An empirical survival function, and the corresponding normal survival
function can be plotted using two `ggplot` "stat" functions provided
by `cmstatr`. In the example below, the empirical survival function
is plotted for each condition, and the survival function for a normal
distribution with the mean and variance from the data is also plotted
(the survival function is 1 minus the cumulative distribution function).
This type of plot can be used to identify how closely the data follows
a normal distribution, and also to compare the distributions of the
various conditions.

```{r}
dat %>%
  ggplot(aes(x = strength, color = condition)) +
  stat_normal_surv_func() +
  stat_esf() +
  theme_bw() +
  ggtitle("Normal Survival Function Plot")
```


# Normal Score Plots
The normal scores plot calculates the normal score and plots it against
the normal score. Normal plots are useful to investigate distributions of
the data.

```{r}
dat %>%
  group_by(condition) %>%
  mutate(norm.score = scale(strength)) %>%
  ggplot(aes(x = norm.score, y = strength, colour = condition)) +
  geom_point() +
  ggtitle("Normal Scores Plot") +
  theme_bw()
```

# Q-Q Plots

A Q-Q plot compares the data against the
theoretical quantiles for a particular distribution.
A line is also plotted showing
the normal distribution with mean and variance from the data. If the data
exactly followed a normal distribution, all points would fall on this line.

```{r}
dat %>%
  ggplot(aes(sample = strength, colour = condition)) +
  geom_qq() +
  geom_qq_line() +
  ggtitle("Q-Q Plot") +
  theme_bw()
```


# Property Plots

Property plots allow for a variety of properties for a group to be 
compared to other properties within the same group, as well as to 
other group properties. The properties included in this plot are A-Basis, 
B-Basis, Pooled A- and B-Basis, Pooled Modified CV (Coefficient of Variation)
A- and B-Basis, Mean, and Min for each group.

The property plots will take a bit of work to construct.

First, the distribution of each group must be determined. Once the
distribution has been determined, the proper basis calculation based 
on that distribution should be filled in below. We also have a
column in the tables below for extra arguments to pass to the `basis`
function, such as overrides required or the method for the `basis_hk_ext`
function to use.

```{r}
b_basis_fcn <- tribble(
  ~condition, ~fcn, ~args,
  "CTD", "basis_normal", list(override = c("between_batch_variability")),
  "RTD", "basis_normal", list(override = c("between_batch_variability")),
  "ETW", "basis_hk_ext", NULL,
  "ETW2", "basis_normal", list(override = c("between_batch_variability"))
)

a_basis_fcn <- tribble(
  ~condition, ~fcn, ~args,
  "CTD", "basis_normal", list(override = c("between_batch_variability")),
  "RTD", "basis_normal", list(override = c("between_batch_variability")),
  "ETW", "basis_hk_ext", list(method = "woodward-frawley"),
  "ETW2", "basis_normal", list(override = c("between_batch_variability"))
)
```

We'll write a function that takes the data and information about
the distribution and computes the single-point
basis value. We'll use this
function for both A- and B-Basis, so we'll add a parameter for the
probability (0.90 or 0.99).

```{r}
single_point_fcn <- function(group_x, group_batch, cond, basis_fcn, p) {
  fcn <- basis_fcn$fcn[basis_fcn$condition == cond[1]]
  extra_args <- basis_fcn$args[basis_fcn$condition == cond[1]]

  args <- c(
    list(x = group_x, batch = group_batch, p = p),
    unlist(extra_args))
  basis <- do.call(fcn, args)
  basis$basis
}

single_point_results <- dat %>%
  group_by(condition) %>%
  summarise(single_point_b_basis = single_point_fcn(
              strength, batch, condition, b_basis_fcn, 0.90),
            single_point_a_basis = single_point_fcn(
              strength, batch, condition, a_basis_fcn, 0.99),
            minimum = min(strength),
            mean = mean(strength)) %>%
  mutate(condition = ordered(condition, c("CTD", "RTD", "ETW", "ETW2")))

single_point_results
```

In the above code, we also ensure that the condition column is still
in the order we expect.

We've already computed the B-Basis of the data using a pooling method.
We'll do the same for A-Basis:

```{r}
a_basis_pooled <- dat %>%
  basis_pooled_cv(strength, condition, batch, p = 0.99,
                  override = c("between_group_variability",
                               "normalized_variance_equal"))

a_basis_pooled
```

As we saw before, the returned object has a property called `basis`,
which is a `data.frame` for the pooling methods.

```{r}
a_basis_pooled$basis
```

We can take this `data.frame` and change the column names
to suit our needs.

```{r}
a_basis_pooled$basis %>%
  rename(condition = group,
         b_basis_pooled = value)
```

We can combine all these steps into one statement. We'll also ensure
that the conditions are listed in the order we want.

```{r}
a_basis_pooled_results <- a_basis_pooled$basis %>%
  rename(condition = group,
         a_basis_pooled = value) %>%
  mutate(condition = ordered(condition, c("CTD", "RTD", "ETW", "ETW2")))

a_basis_pooled_results
```

And the same thing for B-Basis:

```{r}
b_basis_pooled_results <- b_basis_pooled$basis %>%
  rename(condition = group,
         b_basis_pooled = value) %>%
  mutate(condition = ordered(condition, c("CTD", "RTD", "ETW", "ETW2")))

b_basis_pooled_results
```

We can use the function `inner_join` from the `dplyr` package to combine
the three sets of computational results. Each row for each condition
will be concatenated.

```{r}
single_point_results %>%
  inner_join(b_basis_pooled_results, by = "condition") %>%
  inner_join(a_basis_pooled_results, by = "condition")
```

To use this table in the plot we're trying to construct, we want to
"lengthen" the table as follows.

```{r}
single_point_results %>%
  inner_join(b_basis_pooled_results, by = "condition") %>%
  inner_join(a_basis_pooled_results, by = "condition") %>%
  pivot_longer(cols = single_point_b_basis:a_basis_pooled)
```


We can now make a plot based on this:

```{r}
single_point_results %>%
  inner_join(b_basis_pooled_results, by = "condition") %>%
  inner_join(a_basis_pooled_results, by = "condition") %>%
  pivot_longer(cols = single_point_b_basis:a_basis_pooled) %>%
  ggplot(aes(x = condition, y = value)) +
  geom_boxplot(aes(y = strength), data = dat) +
  geom_point(aes(shape = name, color = name)) +
  ggtitle("Property Graph") +
  theme_bw()
```


# Nested Data Plots

`cmstatr` contains the function `nested_data_plot`. This function creates a
plot showing the sources of variation. In the following example, the data
is grouped according to the variables in the `group` argument. The data
is first grouped according to `batch`, then according to `panel`. The
labels located according to the data points that fall under them.
By default, the mean is used, but that `stat` argument can be used to
locate the labels according to `median` or some other statistic.


```{r}
carbon.fabric.2 %>%
  mutate(panel = as.character(panel)) %>%
  filter(test == "WT") %>%
  nested_data_plot(strength,
                   groups = c(batch, panel))
```

Optionally, `fill` or `color` can be set as follows:

```{r}
carbon.fabric.2 %>%
  mutate(panel = as.character(panel)) %>%
  filter(test == "WT" & condition == "RTD") %>%
  nested_data_plot(strength,
                   groups = c(batch, panel),
                   fill = batch,
                   color = panel)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv.R
\name{glance.equiv_change_mean}
\alias{glance.equiv_change_mean}
\title{Glance at a \code{equiv_change_mean} object}
\usage{
\method{glance}{equiv_change_mean}(x, ...)
}
\arguments{
\item{x}{a \code{equiv_change_mean} object returned from
\code{\link[=equiv_change_mean]{equiv_change_mean()}}}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
A one-row \code{\link[tibble:tibble]{tibble::tibble()}} with the following
columns:
\itemize{
\item \code{alpha} the value of alpha passed to this function
\item \code{n_sample} the number of observations in the sample for which
equivalency is being checked. This is either the value \code{n_sample}
passed to this function or the length of the vector \code{data_sample}.
\item \code{mean_sample} the mean of the observations in the sample for
which equivalency is being checked. This is either the value
\code{mean_sample} passed to this function or the mean of the vector
\code{data-sample}.
\item \code{sd_sample} the standard deviation of the observations in the
sample for which equivalency is being checked. This is either the value
\code{mean_sample} passed to this function or the standard deviation of
the vector \code{data-sample}.
\item \code{n_qual} the number of observations in the qualification data
to which the sample is being compared for equivalency. This is either
the value \code{n_qual} passed to this function or the length of the
vector \code{data_qual}.
\item \code{mean_qual} the mean of the qualification data to which the
sample is being compared for equivalency. This is either the value
\code{mean_qual} passed to this function or the mean of the vector
\code{data_qual}.
\item \code{sd_qual} the standard deviation of the qualification data to
which the sample is being compared for equivalency. This is either the
value \code{mean_qual} passed to this function or the standard deviation
of the vector \code{data_qual}.
\item \code{modcv} logical value indicating whether the equivalency
calculations were performed using the modified CV approach
\item \code{sp} the value of the pooled standard deviation. If
\code{modecv = TRUE}, this pooled standard deviation includes the
modification to the qualification CV.
\item \code{t0} the test statistic
\item \code{t_req} the t-value for \eqn{\alpha / 2} and
\eqn{df = n1 + n2 -2}
\item \code{threshold_min} the minimum value of the sample mean that would
result in a pass
\item \code{threshold_max} the maximum value of the sample mean that would
result in a pass
\item \code{result} a character vector of either "PASS" or "FAIL"
indicating the result of the test for change in mean
}
}
\description{
Glance accepts an object of type \code{equiv_change_mean}
and returns a \code{\link[tibble:tibble]{tibble::tibble()}} with
one row of summaries.

Glance does not do any calculations: it just gathers the results in a
tibble.
}
\examples{
x0 <- rnorm(30, 100, 4)
x1 <- rnorm(5, 91, 7)
eq <- equiv_change_mean(data_qual = x0, data_sample = x1, alpha = 0.01)
glance(eq)

## # A tibble: 1 x 14
##   alpha n_sample mean_sample sd_sample n_qual mean_qual sd_qual modcv
##   <dbl>    <int>       <dbl>     <dbl>  <int>     <dbl>   <dbl> <lgl>
## 1  0.01        5        85.8      9.93     30      100.    3.90 FALSE
## # ... with 6 more variables: sp <dbl>, t0 <dbl>, t_req <dbl>,
## #   threshold_min <dbl>, threshold_max <dbl>, result <chr>

}
\seealso{
\code{\link[=equiv_change_mean]{equiv_change_mean()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.R
\name{nonpara_binomial_rank}
\alias{nonpara_binomial_rank}
\title{Rank for distribution-free tolerance bound}
\usage{
nonpara_binomial_rank(n, p, conf)
}
\arguments{
\item{n}{the sample size}

\item{p}{the desired content for the tolerance bound}

\item{conf}{the confidence level for the desired tolerance bound}
}
\value{
The rank corresponding with the desired tolerance bound
}
\description{
Calculates the rank order for finding distribution-free tolerance
bounds for large samples. This function should only be used for
computing B-Basis for samples larger than 28 or A-Basis for samples
larger than 298. This function is used by
\code{\link[=basis_nonpara_large_sample]{basis_nonpara_large_sample()}}.
}
\details{
This function uses the sum of binomial terms to determine the rank
of the ordered statistic that corresponds with the desired tolerance
limit. This approach does not assume any particular distribution. This
approach is described by Guenther (1969) and by CMH-17-1G.

The results of this function have been verified against the tables in
CMH-17-1G and agreement was found for all sample sizes published in
CMH-17-1G for both A- and B-Basis, as well as the sample sizes
\code{n+1} and \code{n-1}, where
\code{n} is the sample size published in CMH-17-1G.

The tables in CMH-17-1G purportedly list the smallest sample sizes
for which a particular rank can be used. That is, for a sample size
one less than the \code{n} published in the table, the next lowest rank
would be used. In some cases, the results of this function disagree by a
rank of one for sample sizes one less than the \code{n} published in the
table. This indicates a disagreement in that sample size at which
the rank should change. This is likely due to numerical
differences in this function and the procedure used to generate the tables.
However, the disagreement is limited to sample 6500 for A-Basis; no
discrepancies have been identified for B-Basis. Since these sample sizes are
uncommon for composite materials
testing, and the difference between subsequent order statistics will be
very small for samples this large, this difference will have no practical
effect on computed tolerance bounds.
}
\examples{
nonpara_binomial_rank(n = 1693, p = 0.99, conf = 0.95)
## [1] 11

# The above example indicates that for a sample of 1693 observations,
# the A-Basis is best approximated as the 11th ordered observation.
# In the example below, the same ordered observation would also be used
# for a sample of size 1702.

nonpara_binomial_rank(n = 1702, p = 0.99, conf = 0.95)
## [1] 11

}
\references{
W. Guenther, “Determination of Sample Size for Distribution-Free
Tolerance Limits,” Jan. 1969.
Available online: \url{https://www.duo.uio.no/handle/10852/48686}

“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
\seealso{
\code{\link[=basis_nonpara_large_sample]{basis_nonpara_large_sample()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv.R
\name{equiv_change_mean}
\alias{equiv_change_mean}
\title{Equivalency based on change in mean value}
\usage{
equiv_change_mean(
  df_qual = NULL,
  data_qual = NULL,
  n_qual = NULL,
  mean_qual = NULL,
  sd_qual = NULL,
  data_sample = NULL,
  n_sample = NULL,
  mean_sample = NULL,
  sd_sample = NULL,
  alpha,
  modcv = FALSE
)
}
\arguments{
\item{df_qual}{(optional) a data.frame containing the qualification data.
Defaults to NULL.}

\item{data_qual}{(optional) a vector of observations from the
"qualification" data to which equivalency is being tested. Or the column of
\code{df_qual} that contains this data. Defaults to NULL}

\item{n_qual}{the number of observations in the qualification data to which
the sample is being compared for equivalency}

\item{mean_qual}{the mean from the qualification data to which the sample
is being compared for equivalency}

\item{sd_qual}{the standard deviation from the qualification data to which
the sample is being compared for equivalency}

\item{data_sample}{a vector of observations from the sample being compared
for equivalency}

\item{n_sample}{the number of observations in the sample being compared for
equivalency}

\item{mean_sample}{the mean of the sample being compared for equivalency}

\item{sd_sample}{the standard deviation of the sample being compared for
equivalency}

\item{alpha}{the acceptable probability of a Type I error}

\item{modcv}{a logical value indicating whether the modified CV approach
should be used. Defaults to \code{FALSE}}
}
\value{
\itemize{
\item \code{call} the expression used to call this function
\item \code{alpha} the value of alpha passed to this function
\item \code{n_sample} the number of observations in the sample for which
equivalency is being checked. This is either the value \code{n_sample}
passed to this function or the length of the vector \code{data_sample}.
\item \code{mean_sample} the mean of the observations in the sample for
which equivalency is being checked. This is either the value
\code{mean_sample} passed to this function or the mean of the vector
\code{data-sample}.
\item \code{sd_sample} the standard deviation of the observations in the
sample for which equivalency is being checked. This is either the value
\code{mean_sample} passed to this function or the standard deviation of
the vector \code{data-sample}.
\item \code{n_qual} the number of observations in the qualification data
to which the sample is being compared for equivalency. This is either
the value \code{n_qual} passed to this function or the length of the
vector \code{data_qual}.
\item \code{mean_qual} the mean of the qualification data to which the
sample is being compared for equivalency. This is either the value
\code{mean_qual} passed to this function or the mean of the vector
\code{data_qual}.
\item \code{sd_qual} the standard deviation of the qualification data to
which the sample is being compared for equivalency. This is either the
value \code{mean_qual} passed to this function or the standard deviation
of the vector \code{data_qual}.
\item \code{modcv} logical value indicating whether the equivalency
calculations were performed using the modified CV approach
\item \code{sp} the value of the pooled standard deviation. If
\code{modecv = TRUE}, this pooled standard deviation includes the
modification to the qualification CV.
\item \code{t0} the test statistic
\item \code{t_req} the t-value for \eqn{\alpha / 2} and
\eqn{df = n1 + n2 -2}
\item \code{threshold} a vector with two elements corresponding to the
minimum and maximum values of the sample mean that would result in a
pass
\item \code{result} a character vector of either "PASS" or "FAIL"
indicating the result of the test for change in mean
}
}
\description{
Checks for change in the mean value between a qualification data set and
a sample. This is normally used to check for properties such as modulus.
This function is a wrapper for a two-sample t--test.
}
\details{
There are several optional arguments to this function. Either (but not both)
\code{data_sample} or all of \code{n_sample}, \code{mean_sample} and
\code{sd_sample} must be supplied. And, either (but not both)
\code{data_qual}
(and also \code{df_qual} if \code{data_qual} is a column name and not a
vector) or all of \code{n_qual}, \code{mean_qual} and \code{sd_qual} must
be supplied. If these requirements are violated, warning(s) or error(s) will
be issued.

This function uses a two-sample t-test to determine if there is a difference
in the mean value of the qualification data and the sample. A pooled
standard deviation is used in the t-test. The procedure is per CMH-17-1G.

If \code{modcv} is TRUE, the standard deviation used to calculate the
thresholds will be replaced with a standard deviation calculated
using the Modified Coefficient of Variation (CV) approach.
The Modified CV approach is a way of adding extra variance to the
qualification data in the case that the qualification data has less
variance than expected, which sometimes occurs when qualification testing
is performed in a short period of time.
Using the Modified CV approach, the standard deviation is calculated by
multiplying \code{CV_star * mean_qual} where \code{mean_qual} is either the
value supplied or the value calculated by \code{mean(data_qual)} and
\eqn{CV*} is determined using \code{\link[=calc_cv_star]{calc_cv_star()}}.

Note that the modified CV option should only be used if that data passes the
Anderson--Darling test.
}
\examples{
equiv_change_mean(alpha = 0.05, n_sample = 9, mean_sample = 9.02,
                  sd_sample = 0.15785, n_qual = 28, mean_qual = 9.24,
                  sd_qual = 0.162, modcv = TRUE)

## Call:
## equiv_change_mean(n_qual = 28, mean_qual = 9.24, sd_qual = 0.162,
##                   n_sample = 9, mean_sample = 9.02, sd_sample = 0.15785,
##                   alpha = 0.05,modcv = TRUE)
##
## For alpha = 0.05
## Modified CV used
##                   Qualification        Sample
##           Number        28               9
##             Mean       9.24             9.02
##               SD      0.162           0.15785
##           Result               PASS
##    Passing Range       8.856695 to 9.623305

}
\references{
“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
\seealso{
\code{\link[=calc_cv_star]{calc_cv_star()}}

\code{\link[stats:t.test]{stats::t.test()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.R
\name{basis}
\alias{basis}
\alias{basis_normal}
\alias{basis_lognormal}
\alias{basis_weibull}
\alias{basis_pooled_cv}
\alias{basis_pooled_sd}
\alias{basis_hk_ext}
\alias{basis_nonpara_large_sample}
\alias{basis_anova}
\title{Calculate basis values}
\usage{
basis_normal(
  data = NULL,
  x,
  batch = NULL,
  p = 0.9,
  conf = 0.95,
  override = c()
)

basis_lognormal(
  data = NULL,
  x,
  batch = NULL,
  p = 0.9,
  conf = 0.95,
  override = c()
)

basis_weibull(
  data = NULL,
  x,
  batch = NULL,
  p = 0.9,
  conf = 0.95,
  override = c()
)

basis_pooled_cv(
  data = NULL,
  x,
  groups,
  batch = NULL,
  p = 0.9,
  conf = 0.95,
  modcv = FALSE,
  override = c()
)

basis_pooled_sd(
  data = NULL,
  x,
  groups,
  batch = NULL,
  p = 0.9,
  conf = 0.95,
  modcv = FALSE,
  override = c()
)

basis_hk_ext(
  data = NULL,
  x,
  batch = NULL,
  p = 0.9,
  conf = 0.95,
  method = c("optimum-order", "woodward-frawley"),
  override = c()
)

basis_nonpara_large_sample(
  data = NULL,
  x,
  batch = NULL,
  p = 0.9,
  conf = 0.95,
  override = c()
)

basis_anova(data = NULL, x, groups, p = 0.9, conf = 0.95, override = c())
}
\arguments{
\item{data}{a data.frame}

\item{x}{the variable in the data.frame for which to find the basis value}

\item{batch}{the variable in the data.frame that contains the batches.}

\item{p}{the content of the tolerance bound. Should be 0.90 for B-Basis
and 0.99 for A-Basis}

\item{conf}{confidence level Should be 0.95 for both A- and B-Basis}

\item{override}{a list of names of diagnostic tests to override,
if desired. Specifying "all" will override all diagnostic
tests applicable to the current method.}

\item{groups}{the variable in the data.frame representing the groups}

\item{modcv}{a logical value indicating whether the modified CV approach
should be used. Only applicable to pooling methods.}

\item{method}{the method for Hanson--Koopmans nonparametric basis values.
should be "optimum-order" for B-Basis and "woodward-frawley"
for A-Basis.}
}
\value{
an object of class \code{basis}
This object has the following fields:
\itemize{
\item \code{call} the expression used to call this function
\item \code{distribution} the distribution used (normal, etc.)
\item \code{p} the value of \eqn{p} supplied
\item \code{conf} the value of \eqn{conf} supplied
\item \code{modcv} a logical value indicating whether the modified
CV approach was used. Only applicable to pooling methods.
\item \code{data} a copy of the data used in the calculation
\item \code{groups} a copy of the groups variable.
Only used for pooling and ANOVA methods.
\item \code{batch} a copy of the batch data used for diagnostic tests
\item \code{modcv_transformed_data} the data after the modified CV transformation
\item \code{override} a vector of the names of diagnostic tests that
were overridden. \code{NULL} if none were overridden
\item \code{diagnostic_results} a named character vector containing the
results of all the diagnostic tests. See the Details section for
additional information
\item \code{diagnostic_failures} a vector containing any diagnostic tests
that produced failures
\item \code{n} the number of observations
\item \code{r} the number of groups, if a pooling method was used.
Otherwise it is NULL.
\item \code{basis} the basis value computed. This is a number
except when pooling methods are used, in which case it is a data.frame.
}
}
\description{
Calculate the basis value for a given data set. There are various functions
to calculate the basis values for different distributions.
The basis value is the lower one-sided tolerance bound of a certain
proportion of the population. For more information on tolerance bounds,
see Meeker, et. al. (2017).
For B-Basis, set the content of tolerance bound to \eqn{p=0.90} and
the confidence level to \eqn{conf=0.95}; for A-Basis, set \eqn{p=0.99} and
\eqn{conf=0.95}. While other tolerance bound
contents and confidence levels may be computed, they are infrequently
needed in practice.

These functions also perform some automated diagnostic
tests of the data prior to calculating the basis values. These diagnostic
tests can be overridden if needed.
}
\details{
\code{data} is an optional argument. If \code{data} is given, it should
be a
\code{data.frame} (or similar object). When \code{data} is specified, the
value of \code{x} is expected to be a variable within \code{data}. If
\code{data} is not specified, \code{x} must be a vector.

When \code{modcv=TRUE} is set, which is only applicable to the
pooling methods,
the data is first modified according to the modified coefficient
of variation (CV)
rules. This modified data is then used when both calculating the
basis values and
also when performing the diagnostic tests. The modified CV approach
is a way of
adding extra variance to datasets with unexpectedly low variance.

\code{basis_normal} calculate the basis value by subtracting \eqn{k} times
the standard deviation from the mean. \eqn{k} is given by
the function \code{\link[=k_factor_normal]{k_factor_normal()}}. The equations in
Krishnamoorthy and Mathew (2008) are used.
\code{basis_normal} also
performs a diagnostic test for outliers (using
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for normality (using
\code{\link[=anderson_darling_normal]{anderson_darling_normal()}}).
If the argument \code{batch} is given, this function also performs
a diagnostic test for outliers within
each batch (using \code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for between batch variability (using
\code{\link[=ad_ksample]{ad_ksample()}}). The argument \code{batch} is only used
for these diagnostic tests.

\code{basis_lognormal} calculates the basis value in the same way
that \code{basis_normal} does, except that the natural logarithm of the
data is taken.

\code{basis_lognormal} function also performs
a diagnostic test for outliers (using
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for normality (using
\code{\link[=anderson_darling_lognormal]{anderson_darling_lognormal()}}).
If the argument \code{batch} is given, this function also performs
a diagnostic test for outliers within
each batch (using \code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for between batch variability (using
\code{\link[=ad_ksample]{ad_ksample()}}). The argument \code{batch} is only used
for these diagnostic tests.

\code{basis_weibull} calculates the basis value for data distributed
according to a Weibull distribution. The confidence level for the
content requested is calculated using the conditional method, as
described in Lawless (1982) Section 4.1.2b. This has good agreement
with tables published in CMH-17-1G. Results differ between this function
and STAT17 by approximately 0.5\\%.

\code{basis_weibull} function also performs
a diagnostic test for outliers (using
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for normality (using
\code{\link[=anderson_darling_weibull]{anderson_darling_weibull()}}).
If the argument \code{batch} is given, this function also performs
a diagnostic test for outliers within
each batch (using \code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for between batch variability (using
\code{\link[=ad_ksample]{ad_ksample()}}). The argument \code{batch} is only used
for these diagnostic tests.

\code{basis_hk_ext} calculates the basis value using the Extended
Hanson--Koopmans method, as described in CMH-17-1G and Vangel (1994).
For nonparametric distributions, this function should be used for samples
up to n=28 for B-Basis and up to \eqn{n=299} for A-Basis.
This method uses a pair of order statistics to determine the basis value.
CMH-17-1G suggests that for A-Basis, the first and last order statistic
is used: this is called the "woodward-frawley" method in this package,
after the paper in which this approach is described (as referenced
by Vangel (1994)). For B-Basis, another approach is used whereby the
first and \code{j-th} order statistic are used to calculate the basis value.
In this approach, the \code{j-th} order statistic is selected to minimize
the difference between the tolerance limit (assuming that the order
statistics are equal to the expected values from a standard normal
distribution) and the population quantile for a standard normal
distribution. This approach is described in Vangel (1994). This second
method (for use when calculating B-Basis values) is called
"optimum-order" in this package.
The results of \code{basis_hk_ext} have been
verified against example results from the program STAT-17. Agreement is
typically well within 0.2\%.

Note that the implementation of \code{hk_ext_z_j_opt} changed after \code{cmstatr}
version 0.8.0. This function is used internally by \code{basis_hk_ext}
when \code{method = "optimum-order"}. This implementation change may mean
that basis values computed using this method may change slightly
after version 0.8.0. However, both implementations seem to be equally
valid. See the included vignette
for a discussion of the differences between the implementation before
and after version 0.8.0, as well as the factors given in CMH-17-1G.
To access this vignette, run: \code{vignette("hk_ext", package = "cmstatr")}

\code{basis_hk_ext} also performs
a diagnostic test for outliers (using
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and performs a pair of tests that the sample size and method selected
follow the guidance described above.
If the argument \code{batch} is given, this function also performs
a diagnostic test for outliers within
each batch (using \code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for between batch variability (using
\code{\link[=ad_ksample]{ad_ksample()}}). The argument \code{batch} is only used
for these diagnostic tests.

\code{basis_nonpara_large_sample} calculates the basis value
using the large sample method described in CMH-17-1G. This method uses
a sum of binomials to determine the rank of the ordered statistic
corresponding with the desired tolerance limit (basis value). Results
of this function have been verified against results of the STAT-17
program.

\code{basis_nonpara_large_sample} also performs
a diagnostic test for outliers (using
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and performs a test that the sample size is sufficiently large.
If the argument \code{batch} is given, this function also performs
a diagnostic test for outliers within
each batch (using \code{\link[=maximum_normed_residual]{maximum_normed_residual()}})
and a diagnostic test for between batch variability (using
\code{\link[=ad_ksample]{ad_ksample()}}). The argument \code{batch} is only used
for these diagnostic tests.

\code{basis_anova} calculates basis values using the ANOVA method.
\code{x} specifies the data (normally strength) and \code{groups}
indicates the group corresponding to each observation. This method is
described in CMH-17-1G, but when the ratio of between-batch mean
square to the within-batch mean square is less than or equal
to one, the tolerance factor is calculated based on pooling the data
from all groups. This approach is recommended by Vangel (1992)
and by Krishnamoorthy and Mathew (2008), and is also implemented
by the software CMH17-STATS and STAT-17.
This function automatically performs a diagnostic
test for outliers within each group
(using \code{\link[=maximum_normed_residual]{maximum_normed_residual()}}) and a test for between
group variability (using \code{\link[=ad_ksample]{ad_ksample()}}) as well as checking
that the data contains at least 5 groups.
This function has been verified against the results of the STAT-17 program.

\code{basis_pooled_sd} calculates basis values by pooling the data from
several groups together. \code{x} specifies the data (normally strength)
and \code{group} indicates the group corresponding to each observation.
This method is described in CMH-17-1G and matches the pooling method
implemented in ASAP 2008.

\code{basis_pooled_cv} calculates basis values by pooling the data from
several groups together. \code{x} specifies the data (normally strength)
and \code{group} indicates the group corresponding to each observation.
This method is described in CMH-17-1G.

\code{basis_pooled_sd} and \code{basis_pooled_cv} both automatically
perform a number of diagnostic tests. Using
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}}, they check that there are no
outliers within each group and batch (provided that \code{batch} is
specified). They check the between batch variability using
\code{\link[=ad_ksample]{ad_ksample()}}. They check that there are no outliers within
each group (pooling all batches) using
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}}. They check for the normality
of the pooled data using \code{\link[=anderson_darling_normal]{anderson_darling_normal()}}.
\code{basis_pooled_sd} checks for equality of variance of all
data using \code{\link[=levene_test]{levene_test()}} and \code{basis_pooled_cv}
checks for equality of variances of all data after transforming it
using \code{\link[=normalize_group_mean]{normalize_group_mean()}}
using \code{\link[=levene_test]{levene_test()}}.

The object returned by these functions includes the named vector
\code{diagnostic_results}. This contains all of the diagnostic tests
performed. The name of each element of the vector corresponds with the
name of the diagnostic test. The contents of each element will be
"P" if the diagnostic test passed, "F" if the diagnostic test failed,
"O" if the diagnostic test was overridden and \code{NA} if the
diagnostic test was skipped (typically because an optional
argument was not supplied).

The following list summarizes the diagnostic tests automatically
performed by each function.
\itemize{
\item \code{basis_normal}
\itemize{
\item \code{outliers_within_batch}
\item \code{between_batch_variability}
\item \code{outliers}
\item \code{anderson_darling_normal}
}
\item \code{basis_lognormal}
\itemize{
\item \code{outliers_within_batch}
\item \code{between_batch_variability}
\item \code{outliers}
\item \code{anderson_darling_lognormal}
}
\item \code{basis_weibull}
\itemize{
\item \code{outliers_within_batch}
\item \code{between_batch_variability}
\item \code{outliers}
\item \code{anderson_darling_weibull}
}
\item \code{basis_pooled_cv}
\itemize{
\item \code{outliers_within_batch}
\item \code{between_group_variability}
\item \code{outliers_within_group}
\item \code{pooled_data_normal}
\item \code{normalized_variance_equal}
}
\item \code{basis_pooled_sd}
\itemize{
\item \code{outliers_within_batch}
\item \code{between_group_variability}
\item \code{outliers_within_group}
\item \code{pooled_data_normal}
\item \code{pooled_variance_equal}
}
\item \code{basis_hk_ext}
\itemize{
\item \code{outliers_within_batch}
\item \code{between_batch_variability}
\item \code{outliers}
\item \code{sample_size}
}
\item \code{basis_nonpara_large_sample}
\itemize{
\item \code{outliers_within_batch}
\item \code{between_batch_variability}
\item \code{outliers}
\item \code{sample_size}
}
\item \code{basis_anova}
\itemize{
\item \code{outliers_within_group}
\item \code{equality_of_variance}
\item \code{number_of_groups}
}
}
}
\examples{
library(dplyr)

# A single-point basis value can be calculated as follows
# in this example, three failed diagnostic tests are
# overridden.

carbon.fabric \%>\%
  filter(test == "FC") \%>\%
  filter(condition == "RTD") \%>\%
  basis_normal(strength, batch,
               override = c("outliers",
                            "outliers_within_batch",
                            "anderson_darling_normal"))

## Call:
## basis_normal(data = ., x = strength, batch = batch,
##     override = c("outliers", "outliers_within_batch",
##    "anderson_darling_normal"))
##
## Distribution:  Normal 	( n = 18 )
## The following diagnostic tests were overridden:
##     `outliers`,
##     `outliers_within_batch`,
##     `anderson_darling_normal`
## B-Basis:   ( p = 0.9 , conf = 0.95 )
## 76.94656

# A set of pooled basis values can also be calculated
# using the pooled standard deviation method, as follows.
# In this example, one failed diagnostic test is overridden.
carbon.fabric \%>\%
  filter(test == "WT") \%>\%
  basis_pooled_sd(strength, condition, batch,
                  override = c("outliers_within_batch"))

## Call:
## basis_pooled_sd(data = ., x = strength, groups = condition,
##                 batch = batch, override = c("outliers_within_batch"))
##
## Distribution:  Normal - Pooled Standard Deviation 	( n = 54, r = 3 )
## The following diagnostic tests were overridden:
##     `outliers_within_batch`
## B-Basis:   ( p = 0.9 , conf = 0.95 )
## CTD  127.6914
## ETW  125.0698
## RTD  132.1457

}
\references{
J. F. Lawless, Statistical Models and Methods for Lifetime Data.
New York: John Wiley & Sons, 1982.

“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.

M. Vangel, “One-Sided Nonparametric Tolerance Limits,”
Communications in Statistics - Simulation and Computation,
vol. 23, no. 4. pp. 1137–1154, 1994.

K. Krishnamoorthy and T. Mathew, Statistical Tolerance Regions: Theory,
Applications, and Computation. Hoboken: John Wiley & Sons, 2008.

W. Meeker, G. Hahn, and L. Escobar, Statistical Intervals: A Guide
for Practitioners and Researchers, Second Edition.
Hoboken: John Wiley & Sons, 2017.

M. Vangel, “New Methods for One-Sided Tolerance Limits for a One-Way
Balanced Random-Effects ANOVA Model,” Technometrics, vol. 34, no. 2.
Taylor & Francis, pp. 176–185, 1992.
}
\seealso{
\code{\link[=hk_ext_z_j_opt]{hk_ext_z_j_opt()}}

\code{\link[=k_factor_normal]{k_factor_normal()}}

\code{\link[=transform_mod_cv]{transform_mod_cv()}}

\code{\link[=maximum_normed_residual]{maximum_normed_residual()}}

\code{\link[=anderson_darling_normal]{anderson_darling_normal()}}

\code{\link[=anderson_darling_lognormal]{anderson_darling_lognormal()}}

\code{\link[=anderson_darling_weibull]{anderson_darling_weibull()}}

\code{\link[=ad_ksample]{ad_ksample()}}

\code{\link[=normalize_group_mean]{normalize_group_mean()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adtest.R
\name{anderson_darling}
\alias{anderson_darling}
\alias{anderson_darling_normal}
\alias{anderson_darling_lognormal}
\alias{anderson_darling_weibull}
\title{Anderson--Darling test for goodness of fit}
\usage{
anderson_darling_normal(data = NULL, x, alpha = 0.05)

anderson_darling_lognormal(data = NULL, x, alpha = 0.05)

anderson_darling_weibull(data = NULL, x, alpha = 0.05)
}
\arguments{
\item{data}{a data.frame-like object (optional)}

\item{x}{a numeric vector or a variable in the data.frame}

\item{alpha}{the required significance level of the test.
Defaults to 0.05.}
}
\value{
an object of class \code{anderson_darling}. This object has the following
fields.
\itemize{
\item \code{call} the expression used to call this function
\item \code{dist} the distribution used
\item \code{data} a copy of the data analyzed
\item \code{n} the number of observations in the sample
\item \code{A} the Anderson--Darling test statistic
\item \code{osl} the observed significance level (p-value),
assuming the
parameters of the distribution are estimated from the data
\item \code{alpha} the required significance level for the test.
This value is given by the user.
\item \code{reject_distribution} a logical value indicating whether
the hypothesis that the data is drawn from the specified distribution
should be rejected
}
}
\description{
Calculates the Anderson--Darling test statistic for a sample given
a particular distribution, and determines whether to reject the
hypothesis that a sample is drawn from that distribution.
}
\details{
The Anderson--Darling test statistic is calculated for the distribution
given by the user.

The observed significance level (OSL), or p-value, is calculated assuming
that the parameters
of the distribution are unknown; these parameters are estimate from the
data.

The function \code{anderson_darling_normal} computes the Anderson--Darling
test statistic given a normal distribution with mean and standard deviation
equal to the sample mean and standard deviation.

The function \code{anderson_darling_lognormal} is the same as
\code{anderson_darling_normal} except that the data is log transformed
first.

The function \code{anderson_darling_weibull} computes the Anderson--Darling
test statistic given a Weibull distribution with shape and scale parameters
estimated from the data using a maximum likelihood estimate.

The test statistic, \code{A}, is modified to account for
the fact that the parameters of the population are not known,
but are instead estimated from the sample. This modification is
a function of the sample size only, and is different for each
distribution (normal/lognormal or Weibull). Several such modifications
have been proposed. This function uses the modification published in
Stephens (1974), Lawless (1982) and CMH-17-1G. Some other implementations
of the Anderson-Darling test, such as the implementation in the
\code{nortest} package, use other modifications, such as the one
published in D'Agostino and Stephens (1986). As such, the p-value
reported by this function may differ from the p-value reported
by implementations of the Anderson--Darling test that use
different modifiers. Only the unmodified
test statistic is reported in the result of this function, but
the modified test statistic is used to compute the OSL (p-value).

This function uses the formulae for observed significance
level (OSL) published in CMH-17-1G. These formulae depend on the particular
distribution used.

The results of this function have been validated against
published values in Lawless (1982).
}
\examples{
library(dplyr)

carbon.fabric \%>\%
  filter(test == "FC") \%>\%
  filter(condition == "RTD") \%>\%
  anderson_darling_normal(strength)
## Call:
## anderson_darling_normal(data = ., x = strength)
##
## Distribution:  Normal ( n = 18 )
## Test statistic:  A = 0.9224776
## OSL (p-value):  0.01212193  (assuming unknown parameters)
## Conclusion: Sample is not drawn from a Normal distribution (alpha = 0.05)

}
\references{
J. F. Lawless, \emph{Statistical models and methods for lifetime data}.
New York: Wiley, 1982.

"Composite Materials Handbook, Volume 1. Polymer Matrix
Composites Guideline for Characterization of Structural
Materials," SAE International, CMH-17-1G, Mar. 2012.

M. A. Stephens, “EDF Statistics for Goodness of Fit and Some
Comparisons,”
Journal of the American Statistical Association, vol. 69, no. 347.
pp. 730–737, 1974.

R. D’Agostino and M. Stephens, Goodness-of-Fit Techniques.
New York: Marcel Dekker, 1986.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnr.R
\name{maximum_normed_residual}
\alias{maximum_normed_residual}
\title{Detect outliers using the maximum normed residual method}
\usage{
maximum_normed_residual(data = NULL, x, alpha = 0.05)
}
\arguments{
\item{data}{a data.frame}

\item{x}{the variable in the data.frame for which to find the MNR
or a vector if \code{data=NULL}}

\item{alpha}{the significance level for the test. Defaults to 0.05}
}
\value{
an object of class \code{mnr}
This object has the following fields:
\itemize{
\item \code{call} the expression used to call this function
\item \code{data} the original data used to compute the MNR
\item \code{alpha} the value of alpha given by the user
\item \code{mnr} the computed MNR test statistic
\item \code{crit} the critical value given the sample size and the
significance level
\item \code{outliers} a data.frame containing the \code{index} and
\code{value} of each of the identified outliers
\item \code{n_outliers} the number of outliers found
}
}
\description{
This function detects outliers using the maximum normed residual
method described in CMH-17-1G. This method identifies a value
as an outlier if the absolute difference between the value and
the sample mean divided by the sample standard deviation
exceeds a critical value.
}
\details{
\code{data} is an optional argument. If \code{data} is given, it
should be a
\code{data.frame} (or similar object). When \code{data} is specified, the
value of \code{x} is expected to be a variable within \code{data}. If
\code{data} is not specified, \code{x} must be a vector.

The maximum normed residual test is a test for outliers. The test statistic
is given in CMH-17-1G. Outliers are identified in the returned object.

The maximum normed residual test statistic is defined as:

\deqn{MNR = max \frac{\left| x_i - \bar{x} \right|}{s} }{
  MNR = max | x_i- x_bar | / s }

When the value of the MNR test statistic exceeds the critical value
defined in Section 8.3.3.1 of CMH-17-1G, the corresponding value
is identified as an outlier. It is then removed from the sample, and
the test statistic is computed again and compared with the critical
value corresponding with the new sample. This process is repeated until
no values are identified as outliers.
}
\examples{
library(dplyr)

carbon.fabric.2 \%>\%
  filter(test=="FC" & condition=="ETW2" & batch=="A") \%>\%
  maximum_normed_residual(strength)

## Call:
## maximum_normed_residual(data = ., x = strength)
##
## MNR =  1.958797  ( critical value = 1.887145 )
##
## Outliers ( alpha = 0.05 ):
##   Index  Value
##       6  44.26

carbon.fabric.2 \%>\%
  filter(test=="FC" & condition=="ETW2" & batch=="B") \%>\%
  maximum_normed_residual(strength)

## Call:
## maximum_normed_residual(data = ., x = strength)
##
## MNR =  1.469517  ( critical value = 1.887145 )
##
## No outliers detected ( alpha = 0.05 )

}
\references{
“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnr.R
\name{augment.mnr}
\alias{augment.mnr}
\title{Augment data with information from an \code{mnr} object}
\usage{
\method{augment}{mnr}(x, data = x$data, ...)
}
\arguments{
\item{x}{an \code{mnr} object created by
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}}}

\item{data}{a \code{data.frame} or
\code{\link[tibble:tibble]{tibble::tibble()}}
containing the original data that was passed to
\code{maximum_normed_residual}}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
When \code{data} is supplied, \code{augment} returns \code{data}, but with
one column appended. When \code{data} is not supplied, \code{augment}
returns a new \code{\link[tibble:tibble]{tibble::tibble()}} with the column
\code{values} containing the original values used by
\code{maximum_normed_residaul} plus one additional column. The additional
column is:
\itemize{
\item \code{.outler} a logical value indicating whether the observation
is an outlier
}
}
\description{
Augment accepts an \code{mnr} object (returned from the function
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}}) and a dataset and adds the column
\code{.outlier} to the dataset. The column \code{.outlier} is a logical
vector indicating whether each observation is an outlier.

When passing data into \code{augment} using the \code{data} argument,
the data must be exactly the data that was passed to
\code{maximum_normed_residual}.
}
\examples{
data <- data.frame(strength = c(80, 98, 96, 97, 98, 120))
m <- maximum_normed_residual(data, strength)

# augment can be called with the original data
augment(m, data)

##   strength .outlier
## 1       80    FALSE
## 2       98    FALSE
## 3       96    FALSE
## 4       97    FALSE
## 5       98    FALSE
## 6      120    FALSE

# or augment can be called without the orignal data and it will be
# reconstructed
augment(m)

## # A tibble: 6 x 2
##   values .outlier
##    <dbl> <lgl>
## 1     80 FALSE
## 2     98 FALSE
## 3     96 FALSE
## 4     97 FALSE
## 5     98 FALSE
## 6    120 FALSE

}
\seealso{
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.R
\name{k_factor_normal}
\alias{k_factor_normal}
\title{Calculate k factor for basis values (\eqn{kB}, \eqn{kA}) with normal
distribution}
\usage{
k_factor_normal(n, p = 0.9, conf = 0.95)
}
\arguments{
\item{n}{the number of observations (i.e. coupons)}

\item{p}{the desired content of the tolerance bound.
Should be 0.90 for B-Basis and 0.99 for A-Basis}

\item{conf}{confidence level. Should be 0.95 for both A- and B-Basis}
}
\value{
the calculated factor
}
\description{
The factors returned by this function are used when calculating basis
values (one-sided confidence bounds) when the data are normally
distributed. The basis value will
be equal to \eqn{\bar{x} - k s}{x_bar - k s},
where \eqn{\bar{x}}{x_bar} is the sample mean,
\eqn{s} is the sample standard deviation and \eqn{k} is the result
of this function.
This function is internally used by \code{\link[=basis_normal]{basis_normal()}} when
computing basis values.
}
\details{
This function calculates the k factors used when determining A- and
B-Basis values for normally distributed data. To get \eqn{kB}, set
the content of the tolerance bound to \code{p = 0.90} and
the confidence level to \code{conf = 0.95}. To get \eqn{kA}, set
\code{p = 0.99} and \code{conf = 0.95}. While other tolerance bound
contents and confidence levels may be computed, they are infrequently
needed in practice.

The k-factor is calculated using equation 2.2.3 of
Krishnamoorthy and Mathew (2008).

This function has been validated against the \eqn{kB} tables in
CMH-17-1G for each value of \eqn{n} from \eqn{n = 2} to \eqn{n = 95}.
It has been validated against the \eqn{kA} tables in CMH-17-1G for each
value of \eqn{n} from \eqn{n = 2} to \eqn{n = 75}. Larger values of \eqn{n}
also match the tables in CMH-17-1G, but R
emits warnings that "full precision may not have been achieved." When
validating the results of this function against the tables in CMH-17-1G,
the maximum allowable difference between the two is 0.002. The tables in
CMH-17-1G give values to three decimal places.

For more information about tolerance bounds in general, see
Meeker, et. al. (2017).
}
\examples{
kb <- k_factor_normal(n = 10, p = 0.9, conf = 0.95)
print(kb)

## [1] 2.35464

# This can be used to caclulate the B-Basis if
# the sample mean and sample standard deviation
# is known, and data is assumed to be normally
# distributed

sample_mean <- 90
sample_sd <- 5.2
print("B-Basis:")
print(sample_mean - sample_sd * kb)

## [1] B-Basis:
## [1] 77.75587

}
\references{
K. Krishnamoorthy and T. Mathew, Statistical Tolerance Regions: Theory,
Applications, and Computation. Hoboken: John Wiley & Sons, 2008.

W. Meeker, G. Hahn, and L. Escobar, Statistical Intervals: A Guide
for Practitioners and Researchers, Second Edition.
Hoboken: John Wiley & Sons, 2017.

“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
\seealso{
\code{\link[=basis_normal]{basis_normal()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.R
\name{hk_ext}
\alias{hk_ext}
\alias{hk_ext_z}
\alias{hk_ext_z_j_opt}
\title{Calculate values related to Extended Hanson--Koopmans tolerance bounds}
\usage{
hk_ext_z(n, i, j, p, conf)

hk_ext_z_j_opt(n, p, conf)
}
\arguments{
\item{n}{the sample size}

\item{i}{the first order statistic (1 <= i < j)}

\item{j}{the second order statistic (i < j <= n)}

\item{p}{the content of the tolerance bound (normally 0.90 or 0.99)}

\item{conf}{the confidence level (normally 0.95)}
}
\value{
For \code{hk_ext_z}, the return value is a numeric value representing
the parameter z (denoted as k in CMH-17-1G).

For \code{hk_ext_z_j_opt}, the return value is named list containing
\code{z} and \code{k}. The former is the value of z, as defined by
Vangel (1994), and the latter is the corresponding order statistic.
}
\description{
Calculates values related to Extended Hanson--Koopmans tolerance bounds
as described by Vangel (1994).
}
\details{
Hanson (1964) presents a nonparametric method for determining
tolerance bounds based on consecutive order statistics.
Vangel (1994) extends this method using non-consecutive order statistics.

The extended Hanson--Koopmans method calculates a tolerance bound
(basis value) based on two order statistics and a weighting value
\code{z}. The value of \code{z} is based on the sample size, which
order statistics are selected, the desired content of the tolerance
bond and the desired confidence level.

The function \code{hk_ext_z} calculates the weighting variable \code{z}
based on selected order statistics \code{i} and \code{j}. Based on this
value \code{z}, the tolerance bound can be calculated as:

\deqn{S = z X_{(i)} + (1 - z) X_{(j)}}{S = z X(i) + (1 - z) X(j)}

Where \eqn{X_{(i)}}{X(i)} and \eqn{X_{(j)}}{X(j)} are the \code{i-th}
and \code{j-th} ordered observation.

The function \code{hk_ext_z_j_opt} determines the value of \code{j} and
the corresponding value of \code{z}, assuming \code{i=1}. The value
of \code{j} is selected such that the computed tolerance limit is
nearest to the desired population quantile for a standard normal
distribution when the order statistics are equal to the expected
value of the order statistics for the standard normal distribution.
}
\examples{
# The factors from Table 1 of Vangel (1994) can be recreated
# using the hk_ext_z function. For the sample size n=21,
# the median is the 11th ordered observation. The factor
# required for calculating the tolerance bound with a content
# of 0.9 and a confidence level of 0.95 based on the median
# and first ordered observation can be calculated as follows.
hk_ext_z(n = 21, i = 1, j = 11, p = 0.9, conf = 0.95)

## [1] 1.204806

# The hk_ext_z_j_opt function can be used to refine this value
# of z by finding an optimum value of j, rather than simply
# using the median. Here, we find that the optimal observation
# to use is the 10th, not the 11th (which is the median).
hk_ext_z_j_opt(n = 21, p = 0.9, conf = 0.95)

## $z
## [1] 1.217717
##
## $j
## [1] 10

}
\references{
M. Vangel, “One-Sided Nonparametric Tolerance Limits,”
Communications in Statistics - Simulation and Computation,
vol. 23, no. 4. pp. 1137–1154, 1994.

D. L. Hanson and L. H. Koopmans,
“Tolerance Limits for the Class of Distributions with Increasing
Hazard Rates,” The Annals of Mathematical Statistics,
vol. 35, no. 4. pp. 1561–1570, 1964.
}
\seealso{
\code{\link[=basis_hk_ext]{basis_hk_ext()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{stat_normal_surv_func}
\alias{stat_normal_surv_func}
\title{Normal Survival Function}
\usage{
stat_normal_surv_func(
  mapping = NULL,
  data = NULL,
  geom = "smooth",
  position = "identity",
  show.legend = NA,
  inherit.aes = TRUE,
  n = 100,
  pad = FALSE,
  ...
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{aes()}.}

\item{data}{The data to be displayed in this layer. This has the
same usage as a \code{ggplot2} \code{stat} function.}

\item{geom}{The geometric object to use to display the data.}

\item{position}{Position argument}

\item{show.legend}{Should this layer be included in the legends?}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetic,
rather than combining with them.}

\item{n}{If \code{NULL}, do not interpolated. Otherwise, the
number of points to interpolate.}

\item{pad}{If \code{TRUE}, pad the ESF with additional points
\verb{(-Inf, 0)} and \verb{(0, Inf)}.}

\item{...}{Other arguments to pass on to \code{layer}.}
}
\description{
The Normal survival function provides a visualization of a
distribution. A normal curve is fit based on the mean and standard
deviation of the data, and the survival function of this normal
curve is plotted. The survival function is simply one minus the
CDF.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/levene.R
\name{levene_test}
\alias{levene_test}
\title{Levene's Test for Equality of Variance}
\usage{
levene_test(data = NULL, x, groups, alpha = 0.05, modcv = FALSE)
}
\arguments{
\item{data}{a data.frame}

\item{x}{the variable in the data.frame or a vector on which to perform the
Levene's test (usually strength)}

\item{groups}{a variable in the data.frame that defines the groups}

\item{alpha}{the significance level (default 0.05)}

\item{modcv}{a logical value indicating whether the modified CV approach
should be used.}
}
\value{
Returns an object of class \code{adk}. This object has the following fields:
\itemize{
\item \code{call} the expression used to call this function
\item \code{data} the original data supplied by the user
\item \code{groups} a vector of the groups used in the computation
\item \code{alpha} the value of alpha specified
\item \code{modcv} a logical value indicating whether the modified
CV approach was used.
\item \code{n} the total number of observations
\item \code{k} the number of groups
\item \code{f} the value of the F test statistic
\item \code{p} the computed p-value
\item \code{reject_equal_variance} a boolean value indicating whether the
null hypothesis that all samples have the same variance is rejected
\item \code{modcv_transformed_data} the data after the modified CV transformation
}
}
\description{
This function performs the Levene's test for equality of variance.
}
\details{
This function performs the Levene's test for equality of variance. The
data is transformed as follows:

\deqn{w_{ij} = \left| x_{ij} - m_i \right|}{wij = | xij - mi |}

Where \eqn{m_i}{mi} is median of the \eqn{ith} group. An F-Test is then
performed on the transformed data.

When \code{modcv=TRUE}, the data from each group is first transformed
according to the modified coefficient of variation (CV) rules before
performing Levene's test.
}
\examples{
library(dplyr)

carbon.fabric.2 \%>\%
  filter(test == "FC") \%>\%
  levene_test(strength, condition)
##
## Call:
## levene_test(data = ., x = strength, groups = condition)
##
## n = 91          k = 5
## F = 3.883818    p-value = 0.00600518
## Conclusion: Samples have unequal variance ( alpha = 0.05 )

}
\references{
“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
\seealso{
\code{\link[=calc_cv_star]{calc_cv_star()}}

\code{\link[=transform_mod_cv]{transform_mod_cv()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adk.R
\name{glance.adk}
\alias{glance.adk}
\title{Glance at a \code{adk} (Anderson--Darling k-Sample) object}
\usage{
\method{glance}{adk}(x, ...)
}
\arguments{
\item{x}{an \code{adk} object}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
A one-row \code{\link[tibble:tibble]{tibble::tibble()}} with the following
columns:
\itemize{
\item \code{alpha} the significance level for the test
\item \code{n} the sample size for the test
\item \code{k} the number of samples
\item \code{sigma} the computed standard deviation of the test statistic
\item \code{ad} the test statistic
\item \code{p} the p-value of the test
\item \code{reject_same_dist} whether the test concludes that the samples
are drawn from different populations
}
}
\description{
Glance accepts an object of type \code{adk} and returns a
\code{\link[tibble:tibble]{tibble::tibble()}} with
one row of summaries.

Glance does not do any calculations: it just gathers the results in a
tibble.
}
\examples{
x <- c(rnorm(20, 100, 5), rnorm(20, 105, 6))
k <- c(rep(1, 20), rep(2, 20))
a <- ad_ksample(x = x, groups = k)
glance(a)

## A tibble: 1 x 7
##   alpha     n     k sigma    ad       p reject_same_dist
##   <dbl> <int> <int> <dbl> <dbl>   <dbl> <lgl>
## 1 0.025    40     2 0.727  4.37 0.00487 TRUE

}
\seealso{
\code{\link[=ad_ksample]{ad_ksample()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-nested.R
\name{nested_data_plot}
\alias{nested_data_plot}
\title{Create a plot of nested sources of variation}
\usage{
nested_data_plot(
  dat,
  x,
  groups = c(),
  stat = "mean",
  ...,
  y_gap = 1,
  divider_color = "grey50",
  point_args = list(),
  dline_args = list(),
  vline_args = list(),
  hline_args = list(),
  label_args = list(),
  connector_args = list()
)
}
\arguments{
\item{dat}{a \code{data.frame} or similar object}

\item{x}{the variable within \code{dat} to plot. Most often this would be a
strength or modulus variable.}

\item{groups}{a vector of variables to group the data by}

\item{stat}{a function for computing the central location for each group.
This is normally "mean" but could be "median" or another
function.}

\item{...}{extra options. See Details.}

\item{y_gap}{the vertical gap between grouping variables}

\item{divider_color}{the color of the lines between grouping variables.
Or \code{NULL} to omit these lines.}

\item{point_args}{arguments to pass to \link[ggplot2:geom_point]{ggplot2::geom_point} when plotting
individual data points.}

\item{dline_args}{arguments to pass to \link[ggplot2:geom_segment]{ggplot2::geom_segment} when plotting
the horizontal lines between data points.}

\item{vline_args}{arguments to pass to \link[ggplot2:geom_segment]{ggplot2::geom_segment} when plotting
vertical lines}

\item{hline_args}{arguments to pass to \link[ggplot2:geom_segment]{ggplot2::geom_segment} when plotting
horizontal lines connecting levels in groups}

\item{label_args}{arguments to pass to \link[ggplot2:geom_text]{ggplot2::geom_label} when plotting
labels}

\item{connector_args}{arguments to pass to \link[ggplot2:geom_point]{ggplot2::geom_point} when
plotting the connection between the vertical lines
and the horizontal lines connecting levels in groups}
}
\description{
Creates a plot showing the breakdown of variation within a sample. This
function uses \link{ggplot2} internally.
}
\details{
Extra options can be included to control aesthetic options. The following
options are supported. Any (or all) can be set to a single variable
in the data set.
\itemize{
\item \code{color}: Controls the color of the data points.
\item \code{fill}: Controls the fill color of the labels. When a particular label
is associated with data points with more than one level of the supplied
variable, the fill is omitted.
}
}
\examples{
library(dplyr)
carbon.fabric.2 \%>\%
  filter(test == "WT" & condition == "RTD") \%>\%
  nested_data_plot(strength,
                   groups = c(batch, panel))

# Labels can be filled too
carbon.fabric.2 \%>\%
  filter(test == "WT" & condition == "RTD") \%>\%
  nested_data_plot(strength,
                   groups = c(batch, panel),
                   fill = batch)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv.R
\name{glance.equiv_mean_extremum}
\alias{glance.equiv_mean_extremum}
\title{Glance at an \code{equiv_mean_extremum} object}
\usage{
\method{glance}{equiv_mean_extremum}(x, ...)
}
\arguments{
\item{x}{an equiv_mean_extremum object returned from
\code{\link[=equiv_mean_extremum]{equiv_mean_extremum()}}}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
A one-row \code{\link[tibble:tibble]{tibble::tibble()}} with the following
columns:
\itemize{
\item \code{alpha} the value of alpha passed to this function
\item \code{n_sample} the number of observations in the sample for which
equivalency is being checked. This is either the value \code{n_sample}
passed to this function or the length of the vector \code{data_sample}.
\item \code{modcv} logical value indicating whether the acceptance
thresholds are calculated using the modified CV approach
\item \code{threshold_min_indiv} The calculated threshold value for
minimum individual
\item \code{threshold_mean} The calculated threshold value for mean
\item \code{result_min_indiv} a character vector of either "PASS" or
"FAIL" indicating whether the data from \code{data_sample} passes the
test for minimum individual. If \code{data_sample} was not supplied,
this value will be \code{NULL}
\item \code{result_mean} a character vector of either "PASS" or
"FAIL" indicating whether the data from \code{data_sample} passes the
test for mean. If \code{data_sample} was not supplied, this value will
be  \code{NULL}
\item \code{min_sample} The minimum value from the vector
\code{data_sample}. if \code{data_sample} was not supplied, this will
have a value of \code{NULL}
\item \code{mean_sample} The mean value from the vector
\code{data_sample}. If \code{data_sample} was not supplied, this will
have a value of \code{NULL}
}
}
\description{
Glance accepts an object of type \code{equiv_mean_extremum} and returns a
\code{\link[tibble:tibble]{tibble::tibble()}} with
one row of summaries.

Glance does not do any calculations: it just gathers the results in a
tibble.
}
\examples{
x0 <- rnorm(30, 100, 4)
x1 <- rnorm(5, 91, 7)
eq <- equiv_mean_extremum(data_qual = x0, data_sample = x1, alpha = 0.01)
glance(eq)

## # A tibble: 1 x 9
##   alpha n_sample modcv threshold_min_indiv threshold_mean
##   <dbl>    <int> <lgl>               <dbl>          <dbl>
## 1  0.01        5 FALSE                86.2           94.9
## # ... with 4 more variables: result_min_indiv <chr>, result_mean <chr>,
## #   min_sample <dbl>, mean_sample <dbl>

}
\seealso{
\code{\link[=equiv_mean_extremum]{equiv_mean_extremum()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{transform_mod_cv}
\alias{transform_mod_cv}
\alias{transform_mod_cv_ad}
\title{Transforms data according to the modified CV rule}
\usage{
transform_mod_cv_ad(x, condition, batch)

transform_mod_cv(x)
}
\arguments{
\item{x}{a vector of data to transform}

\item{condition}{a vector indicating the condition to which each
observation belongs}

\item{batch}{a vector indicating the batch to which each observation
belongs}
}
\value{
A vector of transformed data
}
\description{
Transforms data according to the modified coefficient of variation (CV)
rule. This is used to add additional variance to datasets with
unexpectedly low variance, which is sometimes encountered during
testing of new materials over short periods of time.

Two versions of this transformation are implemented. The first version,
\code{transform_mod_cv()}, transforms the data in a single group (with
no other structure) according to the modified CV rules.

The second
version, \code{transform_mod_cv_ad()}, transforms data that is structured
according to both condition and batch, as is commonly done for
the Anderson--Darling k-Sample and Anderson-Darling tests when pooling
across environments.
}
\details{
The modified CV transformation takes the general form:

\deqn{\frac{S_i^*}{S_i} (x_{ij} - \bar{x_i}) + \bar{x_i}}{
  Si*/Si (xij-x_bar_i) + x_bar_i
}

Where \eqn{S_i^*}{Si*} is the modified standard deviation
(mod CV times mean) for
the \eqn{ith} group; \eqn{S_i}{Si} is the standard deviation
for the \eqn{ith} group, \eqn{\bar{x_i}}{x_bar_i} is
the group mean and \eqn{x_{ij}}{xij} is the observation.

\code{transform_mod_cv()} takes a vector
containing the observations and transforms the data.
The equation above is used, and all observations
are considered to be from the same group.

\code{transform_mod_cv_ad()} takes a vector containing the observations
plus a vector containing the corresponding conditions and a vector
containing the batches. This function first calculates the modified
CV value from the data from each condition (independently). Then,
within each condition, the transformation
above is applied to produce the transformed data \eqn{x'}.
This transformed data is further transformed using the following
equation.

\deqn{x_{ij}'' = C (x'_{ij} - \bar{x_i}) + \bar{x_i}}{
  x_ij'' = C (x'_ij - x_bar_i) + x_bar_i}

Where:

\deqn{C = \sqrt{\frac{SSE^*}{SSE'}}}{C = sqrt(SSE* / SSE')}

\deqn{SSE^* = (n-1) (CV^* \bar{x})^2 - \sum(n_i(\bar{x_i}-\bar{x})^2)}{
  SSE* = (n-1) (CV* x_bar)^2 - sum(n_i(x_bar_i-x_bar)^2)}

\deqn{SSE' = \sum(x'_{ij} - \bar{x_i})^2}{SSE' = sum(x'_ij - x_bar_i)^2}
}
\examples{
# Transform data according to the modified CV transformation
# and report the original and modified CV for each condition

library(dplyr)
carbon.fabric \%>\%
filter(test == "FT") \%>\%
  group_by(condition) \%>\%
  mutate(trans_strength = transform_mod_cv(strength)) \%>\%
  head(10)

## # A tibble: 10 x 6
## # Groups:   condition [1]
##    id         test  condition batch strength trans_strength
##    <chr>      <chr> <chr>     <int>    <dbl>          <dbl>
##  1 FT-RTD-1-1 FT    RTD           1     126.           126.
##  2 FT-RTD-1-2 FT    RTD           1     139.           141.
##  3 FT-RTD-1-3 FT    RTD           1     116.           115.
##  4 FT-RTD-1-4 FT    RTD           1     132.           133.
##  5 FT-RTD-1-5 FT    RTD           1     129.           129.
##  6 FT-RTD-1-6 FT    RTD           1     130.           130.
##  7 FT-RTD-2-1 FT    RTD           2     131.           131.
##  8 FT-RTD-2-2 FT    RTD           2     124.           124.
##  9 FT-RTD-2-3 FT    RTD           2     125.           125.
## 10 FT-RTD-2-4 FT    RTD           2     120.           119.

# The CV of this transformed data can be computed to verify
# that the resulting CV follows the rules for modified CV

carbon.fabric \%>\%
  filter(test == "FT") \%>\%
  group_by(condition) \%>\%
  mutate(trans_strength = transform_mod_cv(strength)) \%>\%
  summarize(cv = sd(strength) / mean(strength),
            mod_cv = sd(trans_strength) / mean(trans_strength))

## # A tibble: 3 x 3
##   condition     cv mod_cv
##   <chr>      <dbl>  <dbl>
## 1 CTD       0.0423 0.0612
## 2 ETW       0.0369 0.0600
## 3 RTD       0.0621 0.0711

}
\seealso{
\code{\link[=calc_cv_star]{calc_cv_star()}}

\code{\link[=cv]{cv()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnr.R
\name{glance.mnr}
\alias{glance.mnr}
\title{Glance at a \code{mnr} (maximum normed residual) object}
\usage{
\method{glance}{mnr}(x, ...)
}
\arguments{
\item{x}{An \code{mnr} object}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
A one-row \code{\link[tibble:tibble]{tibble::tibble()}} with the following
columns:
\itemize{
\item \code{mnr} the computed MNR test statistic
\item \code{alpha} the value of alpha used for the test
\item \code{crit} the critical value given the sample size and the
significance level
\item \code{n_outliers} the number of outliers found
}
}
\description{
Glance accepts an object of type \code{mnr} and returns a
\code{\link[tibble:tibble]{tibble::tibble()}} with
one row of summaries.

Glance does not do any calculations: it just gathers the results in a
tibble.
}
\examples{
x <- c(rnorm(20, 100, 5), 10)
m <- maximum_normed_residual(x = x)
glance(m)

## # A tibble: 1 x 4
##     mnr alpha  crit n_outliers
##   <dbl> <dbl> <dbl>      <dbl>
## 1  4.23  0.05  2.73          1

}
\seealso{
\code{\link[=maximum_normed_residual]{maximum_normed_residual()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.R
\name{cv}
\alias{cv}
\title{Calculate the coefficient of variation}
\usage{
cv(x, na.rm = FALSE)
}
\arguments{
\item{x}{a vector}

\item{na.rm}{logical. Should missing values be removed?}
}
\value{
The calculated CV
}
\description{
The coefficient of variation (CV) is the ratio of the standard
deviation to the mean of a sample. This function takes a vector
of data and calculates the CV.
}
\examples{
set.seed(15)  # make this example reproducible
x <- rnorm(100, mean = 100, sd = 5)
cv(x)
## [1] 0.04944505

# the cv function can also be used within a call to dplyr::summarise
library(dplyr)
carbon.fabric \%>\%
filter(test == "WT") \%>\%
  group_by(condition) \%>\%
  summarise(mean = mean(strength), cv = cv(strength))

## # A tibble: 3 x 3
##   condition  mean     cv
##   <chr>     <dbl>  <dbl>
## 1 CTD        137. 0.0417
## 2 ETW        135. 0.0310
## 3 RTD        142. 0.0451


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{normalize_group_mean}
\alias{normalize_group_mean}
\title{Normalize values to group means}
\usage{
normalize_group_mean(x, group)
}
\arguments{
\item{x}{the variable containing the data to normalized}

\item{group}{the variable containing the groups}
}
\value{
Returns a vector of normalized values
}
\description{
This function computes the mean of each group, then divides each
observation by its corresponding group mean. This is commonly done
when pooling data across environments.
}
\details{
Computes the mean for each group, then divides each value by the mean for
the corresponding group.
}
\examples{
library(dplyr)
carbon.fabric.2 \%>\%
filter(test == "WT") \%>\%
  select(condition, strength) \%>\%
  mutate(condition_norm = normalize_group_mean(strength, condition)) \%>\%
  head(10)

##    condition strength condition_norm
## 1        CTD  142.817      1.0542187
## 2        CTD  135.901      1.0031675
## 3        CTD  132.511      0.9781438
## 4        CTD  135.586      1.0008423
## 5        CTD  125.145      0.9237709
## 6        CTD  135.203      0.9980151
## 7        CTD  128.547      0.9488832
## 8        CTD  127.709      0.9426974
## 9        CTD  127.074      0.9380101
## 10       CTD  126.879      0.9365706

}
\references{
“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adk.R
\name{ad_ksample}
\alias{ad_ksample}
\title{Anderson--Darling K-Sample Test}
\usage{
ad_ksample(data = NULL, x, groups, alpha = 0.025)
}
\arguments{
\item{data}{a data.frame}

\item{x}{the variable in the data.frame on which to perform the
Anderson--Darling k-Sample test (usually strength)}

\item{groups}{a variable in the data.frame that defines the groups}

\item{alpha}{the significance level (default 0.025)}
}
\value{
Returns an object of class \code{adk}. This object has the following fields:
\itemize{
\item \code{call} the expression used to call this function
\item \code{data} the original data used to compute the ADK
\item \code{groups} a vector of the groups used in the computation
\item \code{alpha} the value of alpha specified
\item \code{n} the total number of observations
\item \code{k} the number of groups
\item \code{sigma} the computed standard deviation of the test statistic
\item \code{ad} the value of the Anderson--Darling k-Sample test statistic
\item \code{p} the computed p-value
\item \code{reject_same_dist} a boolean value indicating whether the null
hypothesis that all samples come from the same distribution is rejected
\item \code{raw} the original results returned from
\link[kSamples:ad.test]{ad.test}
}
}
\description{
This function performs an Anderson--Darling k-sample test. This is used to
determine if several samples (groups) share a common (unspecified)
distribution.
}
\details{
This function is a wrapper for the \link[kSamples:ad.test]{ad.test} function from
the package \code{kSamples}. The method "exact" is specified in the call to
\code{ad.test}. Refer to that package's documentation for details.

There is a minor difference in the formulation of the Anderson--Darling
k-Sample test in CMH-17-1G, compared with that in the Scholz and
Stephens (1987). This difference affects the test statistic and the
critical value in the same proportion, and therefore the conclusion of
the test is unaffected. When
comparing the test statistic generated by this function to that generated
by software that uses the CMH-17-1G formulation (such as ASAP, CMH17-STATS,
etc.), the test statistic reported by this function will be greater by
a factor of \eqn{(k - 1)}, with a corresponding change in the critical
value.

For more information about the difference between this function and
the formulation in CMH-17-1G, see the vignette on the subject, which
can be accessed by running \code{vignette("adktest")}
}
\examples{
library(dplyr)

carbon.fabric \%>\%
  filter(test == "WT") \%>\%
  filter(condition == "RTD") \%>\%
  ad_ksample(strength, batch)
##
## Call:
## ad_ksample(data = ., x = strength, groups = batch)
##
## N = 18          k = 3
## ADK = 0.912     p-value = 0.95989
## Conclusion: Samples come from the same distribution ( alpha = 0.025 )

}
\references{
F. W. Scholz and M. Stephens, “K-Sample Anderson--Darling Tests,” Journal
of the American Statistical Association, vol. 82, no. 399. pp. 918–924,
Sep-1987.

“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adtest.R
\name{glance.anderson_darling}
\alias{glance.anderson_darling}
\title{Glance at an \code{anderson_darling} object}
\usage{
\method{glance}{anderson_darling}(x, ...)
}
\arguments{
\item{x}{an \code{anderson_darling} object}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
A one-row \code{\link[tibble:tibble]{tibble::tibble()}} with the following
columns:
\itemize{
\item \code{dist} the distribution used
\item \code{n} the number of observations in the sample
\item \code{A} the Anderson--Darling test statistic
\item \code{osl} the observed significance level (p-value),
assuming the
parameters of the distribution are estimated from the data
\item \code{alpha} the required significance level for the test.
This value is given by the user.
\item \code{reject_distribution} a logical value indicating whether
the hypothesis that the data is drawn from the specified distribution
should be rejected
}
}
\description{
Glance accepts an object of type \code{anderson_darling} and
returns a \code{\link[tibble:tibble]{tibble::tibble()}} with
one row of summaries.

Glance does not do any calculations: it just gathers the results in a
tibble.
}
\examples{
x <- rnorm(100, 100, 4)
ad <- anderson_darling_weibull(x = x)
glance(ad)

## # A tibble: 1 x 6
##   dist        n     A        osl alpha reject_distribution
##   <chr>   <int> <dbl>      <dbl> <dbl> <lgl>
## 1 Weibull   100  2.62 0.00000207  0.05 TRUE

}
\seealso{
\code{\link[=anderson_darling]{anderson_darling()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmstatr.R
\docType{package}
\name{cmstatr-package}
\alias{cmstatr}
\alias{cmstatr-package}
\title{cmstatr: Statistical Methods for Composite Material Data}
\description{
To learn more about \code{cmstatr}, start with the vignettes:
\code{browseVignettes(package = "cmstatr")}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://www.cmstatr.net/}
  \item \url{https://github.com/cmstatr/cmstatr}
  \item Report bugs at \url{https://github.com/cmstatr/cmstatr/issues}
}

}
\author{
\strong{Maintainer}: Stefan Kloppenborg \email{stefan@kloppenborg.ca}

Other contributors:
\itemize{
  \item Billy Cheng \email{bcheng@comtekadvanced.com} [contributor]
  \item Ally Fraser \email{ally.fraser25@gmail.com} [contributor]
  \item Comtek Advanced Structures, Ltd. [funder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{augment}
\alias{tidy}
\alias{glance}
\title{Objects exported from other packages}
\seealso{
\code{\link[generics:augment]{generics::augment()}}

\code{\link[generics:tidy]{generics::tidy()}}

\code{\link[generics:glance]{generics::glance()}}
}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{generics}{\code{\link[generics]{augment}}, \code{\link[generics]{glance}}, \code{\link[generics]{tidy}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv.R
\name{k_equiv}
\alias{k_equiv}
\title{k-factors for determining acceptance based on sample mean and an extremum}
\usage{
k_equiv(alpha, n)
}
\arguments{
\item{alpha}{the acceptable probability of a type I error}

\item{n}{the number of observations in the sample to test}
}
\value{
a vector with elements c(k1, k2). k1 is for testing the sample
extremum. k2 is for testing the sample mean
}
\description{
k-factors for determining acceptance based on sample mean and an extremum
}
\details{
The k-factors returned by this function are used for determining
whether to accept a new dataset.

This function is used as part of the procedure for
determining acceptance limits for a sample mean and sample minimum.
These acceptance limits are often used to set acceptance limits for
material strength for each lot of material, or each new manufacturing
site. When a sample meets the criteria that its mean and its minimum are
both greater than these limits, then one may accept the lot of material
or the new manufacturing site.

This procedure is used to ensure that the strength of material processed
at a second site, or made with a new batch of material are not degraded
relative to the data originally used to determine basis values for the
material. For more information about the use of this procedure, see
CMH-17-1G or PS-ACE 100-2002-006.

According to Vangel (2002), the use of mean and extremum for this purpose
is more powerful than the use of mean and standard deviation.

The results of this function match those published by Vangel within
0.05\\% for \eqn{n > 2} and \eqn{\alpha \le 0.25}. Those results published
by Vangel are identical to those published in CMH-17-1G.

This function uses numerical integration and numerical optimization to
find values of the factors \eqn{k_1} and \eqn{k_2} based on Vangel's
saddle point approximation.

The value \eqn{n} refers to the number of observations in the sample
being compared with the original population (the qualification sample is
usually assumed to be equal to the population statistics).

The value of \eqn{alpha} is the acceptable probability of a type I error.
Normally, this is set to 0.05 for material or process equivalency and 0.01
when setting lot acceptance limits. Though, in principle, this parameter
can be set to any number between 0 and 1. This function, however, has only
been validated in the range of \eqn{1e-5 \le alpha \le 0.5}.
}
\examples{
qual_mean <- 100
qual_sd <- 3.5
k <- k_equiv(0.01, 5)
print("Minimum Individual Acceptance Limit:")
print(qual_mean - qual_sd * k[1])
print("Minimum Average Acceptance Limit:")
print(qual_mean - qual_sd * k[2])

## [1] "Minimum Individual Acceptance Limit:"
## [1] 89.24981
## [1] "Minimum Average Acceptance Limit:"
## [1] 96.00123

}
\references{
M. G. Vangel. Lot Acceptance and Compliance Testing Using the Sample Mean
and an Extremum, Technometrics, vol. 44, no. 3. pp. 242–249. 2002.

“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.

Federal Aviation Administration, “Material Qualification and Equivalency
for Polymer Matrix Composite Material Systems,” PS-ACE 100-2002-006,
Sep. 2003.
}
\seealso{
\code{\link[=equiv_mean_extremum]{equiv_mean_extremum()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{stat_esf}
\alias{stat_esf}
\title{Empirical Survival Function}
\usage{
stat_esf(
  mapping = NULL,
  data = NULL,
  geom = "point",
  position = "identity",
  show.legend = NA,
  inherit.aes = TRUE,
  n = NULL,
  pad = FALSE,
  ...
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{aes()}.}

\item{data}{The data to be displayed in this layer. This has the
same usage as a \code{ggplot2} \code{stat} function.}

\item{geom}{The geometric object to use to display the data.}

\item{position}{Position argument}

\item{show.legend}{Should this layer be included in the legends?}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetic,
rather than combining with them.}

\item{n}{If \code{NULL}, do not interpolated. Otherwise, the
number of points to interpolate.}

\item{pad}{If \code{TRUE}, pad the ESF with additional points
\verb{(-Inf, 0)} and \verb{(0, Inf)}.}

\item{...}{Other arguments to pass on to \code{layer}.}
}
\description{
The empirical survival function (ESF) provides a visualization of a
distribution. This is closely related to the empirical cumulative
distribution function (ECDF). The empirical survival function is
simply ESF = 1 - ECDF.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/levene.R
\name{glance.levene}
\alias{glance.levene}
\title{Glance at a \code{levene} object}
\usage{
\method{glance}{levene}(x, ...)
}
\arguments{
\item{x}{a \code{levene} object returned from \code{\link[=levene_test]{levene_test()}}}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
A one-row \code{\link[tibble:tibble]{tibble::tibble()}} with the following
columns:
\itemize{
\item \code{alpha} the value of alpha specified
\item \code{modcv} a logical value indicating whether the modified
CV approach was used.
\item \code{n} the total number of observations
\item \code{k} the number of groups
\item \code{f} the value of the F test statistic
\item \code{p} the computed p-value
\item \code{reject_equal_variance} a boolean value indicating whether the
null hypothesis that all samples have the same variance is rejected
}
}
\description{
Glance accepts an object of type \code{levene} and returns a
\code{\link[tibble:tibble]{tibble::tibble()}} with
one row of summaries.

Glance does not do any calculations: it just gathers the results in a
tibble.
}
\examples{
df <- data.frame(
  groups = c(rep("A", 5), rep("B", 6)),
  strength = c(rnorm(5, 100, 6), rnorm(6, 105, 7))
)
levene_result <- levene_test(df, strength, groups)
glance(levene_result)

## # A tibble: 1 x 7
##   alpha modcv     n     k      f     p reject_equal_variance
##   <dbl> <lgl> <int> <int>  <dbl> <dbl> <lgl>
## 1  0.05 FALSE    11     2 0.0191 0.893 FALSE

}
\seealso{
\code{\link[=levene_test]{levene_test()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{carbon.fabric}
\alias{carbon.fabric}
\alias{carbon.fabric.2}
\title{Sample data for a generic carbon fabric}
\format{
An object of class \code{data.frame} with 216 rows and 5 columns.

An object of class \code{data.frame} with 177 rows and 9 columns.
}
\usage{
carbon.fabric

carbon.fabric.2
}
\description{
Datasets containing sample data that is typical of a generic carbon
fabric prepreg. This data is used in several examples within the
\code{cmstatr} package. This data is fictional and should
only be used for learning how to use this package.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.R
\name{glance.basis}
\alias{glance.basis}
\title{Glance at a basis object}
\usage{
\method{glance}{basis}(x, include_diagnostics = FALSE, ...)
}
\arguments{
\item{x}{a basis object}

\item{include_diagnostics}{a logical value indicating whether to include
columns for diagnostic tests. Default FALSE.}

\item{...}{Additional arguments. Not used. Included only to match generic
signature.}
}
\value{
A \code{\link[tibble:tibble]{tibble::tibble()}} with the following
columns:
\itemize{
\item \code{p} the the content of the tolerance bound. Normally 0.90 or 0.99
\item \code{conf} the confidence level. Normally 0.95
\item \code{distribution} a string representing the distribution assumed
when calculating the basis value
\item \code{modcv} a logical value indicating whether the modified
CV approach was used. Only applicable to pooling methods.
\item \code{n} the sample size
\item \code{r} the number of groups used in the calculation. This will
be \code{NA} for single-point basis values
\item \code{basis} the basis value
}
}
\description{
Glance accepts an object of type basis and returns a
\code{\link[tibble:tibble]{tibble::tibble()}} with
one row of summaries for each basis value.

Glance does not do any calculations: it just gathers the results in a
tibble.
}
\details{
For the pooled basis methods (\code{basis_pooled_cv} and
\code{basis_pooled_sd}), the \code{\link[tibble:tibble]{tibble::tibble()}}
returned by \code{glance} will have one row for each group included in
the pooling. For all other basis methods, the resulting \code{tibble}
will have a single row.

If \code{include_diagnostics=TRUE}, there will be additional columns
corresponding with the diagnostic tests performed. These column(s) will
be of type character and will contain a "P" if the diagnostic test
passed, a "F" if the diagnostic test failed, an "O" if the diagnostic
test was overridden or \code{NA} if the test was not run (typically
because an optional argument was not passed to the function that
computed the basis value).
}
\examples{
set.seed(10)
x <- rnorm(20, 100, 5)
b <- basis_normal(x = x)
glance(b)

## # A tibble: 1 x 7
##       p  conf distribution modcv     n r     basis
##   <dbl> <dbl> <chr>        <lgl> <int> <lgl> <dbl>
## 1   0.9  0.95 Normal       FALSE    20 NA     92.0


glance(b, include_diagnostics = TRUE)

## # A tibble: 1 x 11
##        p  conf distribution modcv     n r     basis outliers_within…
##    <dbl> <dbl> <chr>        <lgl> <int> <lgl> <dbl> <chr>
##  1   0.9  0.95 Normal       FALSE    20 NA     92.0 NA
## # … with 3 more variables: between_batch_variability <chr>,
## #   outliers <chr>, anderson_darling_normal <chr>

}
\seealso{
\code{\link[=basis]{basis()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv.R
\name{equiv_mean_extremum}
\alias{equiv_mean_extremum}
\title{Test for decrease in mean or minimum individual}
\usage{
equiv_mean_extremum(
  df_qual = NULL,
  data_qual = NULL,
  mean_qual = NULL,
  sd_qual = NULL,
  data_sample = NULL,
  n_sample = NULL,
  alpha,
  modcv = FALSE
)
}
\arguments{
\item{df_qual}{(optional) a data.frame containing the qualification data.
Defaults to NULL.}

\item{data_qual}{(optional) a vector of observations from the
"qualification" data to which equivalency is being tested. Or the column of
\code{df_qual} that contains this data. Defaults to NULL}

\item{mean_qual}{(optional) the mean from the "qualification" data to which
equivalency is being tested. Defaults to NULL}

\item{sd_qual}{(optional) the standard deviation from the "qualification"
data to which equivalency is being tested. Defaults to NULL}

\item{data_sample}{(optional) a vector of observations from the sample for
which equivalency is being tested. Defaults to NULL}

\item{n_sample}{(optional) the number of observations in the sample for
which equivalency will be tested. Defaults to NULL}

\item{alpha}{the acceptable probability of a type I error}

\item{modcv}{(optional) a boolean value indicating whether a modified CV
should be used. Defaults to FALSE, in which case the standard deviation
supplied (or calculated from \code{data_qual}) will be used directly.}
}
\value{
Returns an object of class \code{equiv_mean_extremum}. This object is a list
with the following named elements:
\itemize{
\item \code{call} the expression used to call this function
\item \code{alpha} the value of alpha passed to this function
\item \code{n_sample} the number of observations in the sample for which
equivalency is being checked. This is either the value \code{n_sample}
passed to this function or the length of the vector \code{data_sample}.
\item \code{k1} the factor used to calculate the minimum individual
threshold. The minimum individual threshold is calculated as
\eqn{W_{min} = qual\,mean - k_1 \cdot qual\,sd}{
  Wmin = qual_mean - k1 * qual_sd}
\item \code{k2} the factor used to calculate the threshold for mean. The
threshold for mean is calculated as
\eqn{W_{mean} = qual\,mean - k_2 \cdot qual\,sd}{
  Wmean = qual_mean - k2 * qual_sd}
\item \code{modcv} logical value indicating whether the acceptance
thresholds are calculated using the modified CV approach
\item \code{cv} the coefficient of variation of the qualification data.
This value is not modified, even if \code{modcv=TRUE}
\item \code{cv_star} The modified coefficient of variation. If
\code{modcv=FALSE}, this will be \code{NULL}
\item \code{threshold_min_indiv} The calculated threshold value for
minimum individual
\item \code{threshold_mean} The calculated threshold value for mean
\item \code{result_min_indiv} a character vector of either "PASS" or
"FAIL" indicating whether the data from \code{data_sample} passes the
test for minimum individual. If \code{data_sample} was not supplied,
this value will be \code{NULL}
\item \code{result_mean} a character vector of either "PASS" or
"FAIL" indicating whether the data from \code{data_sample} passes the
test for mean. If \code{data_sample} was not supplied, this value will
be  \code{NULL}
\item \code{min_sample} The minimum value from the vector
\code{data_sample}. if \code{data_sample} was not supplied, this will
have a value of \code{NULL}
\item \code{mean_sample} The mean value from the vector
\code{data_sample}. If \code{data_sample} was not supplied, this will
have a value of \code{NULL}
}
}
\description{
This test is used when determining if a new process or
manufacturing location produces material properties that are
"equivalent" to an existing dataset, and hence the existing
basis values are applicable to the new dataset. This test is also
sometimes used for determining if a new batch of material is acceptable.
This function determines thresholds based on both minimum
individual and mean, and optionally evaluates a sample against those
thresholds. The joint distribution between the sample mean
and sample minimum is used to generate these thresholds.
When there is no true difference between the existing ("qualification")
and the new population from which the sample is obtained, there is a
probability of \eqn{\alpha} of falsely concluding that there is a
difference in mean or variance. It is assumed that both the original
and new populations are normally distributed.
According to Vangel (2002), this test provides improved power compared
with a test of mean and standard deviation.
}
\details{
This function is used to
determine acceptance limits for a sample mean and sample minimum.
These acceptance limits are often used to set acceptance limits for
material strength for each lot of material, or each new manufacturing
site. When a sample meets the criteria that its mean and its minimum are
both greater than these limits, then one may accept the lot of material
or the new manufacturing site.

This procedure is used to ensure that the strength of material processed
at a second site, or made with a new batch of material are not degraded
relative to the data originally used to determine basis values for the
material. For more information about the use of this procedure, see
CMH-17-1G or PS-ACE 100-2002-006.

There are several optional arguments to this function. However, you can't
omit all of the optional arguments. You must supply either
\code{data_sample} or \code{n_sample}, but not both. You must also supply
either \code{data_qual} (and \code{df_qual} if \code{data_qual} is a
variable name and not a vector) or both \code{mean_qual} and \code{sd_qual},
but if you supply \code{data_qual} (and possibly \code{df_qual}) you should
not supply either \code{mean_qual} or \code{sd_qual} (and visa-versa). This
function will issue a warning or error if you violate any of these rules.

If \code{modcv} is TRUE, the standard deviation used to calculate the
thresholds will be replaced with a standard deviation calculated
using the Modified Coefficient of Variation (CV) approach.
The Modified CV approach is a way of adding extra variance to the
qualification data in the case that the qualification data has less
variance than expected, which sometimes occurs when qualification testing
is performed in a short period of time.
Using the Modified CV approach, the standard deviation is calculated by
multiplying \code{CV_star * mean_qual} where \code{mean_qual} is either the
value supplied or the value calculated by \code{mean(data_qual)} and
\eqn{CV*} is the value computed by \code{\link[=calc_cv_star]{calc_cv_star()}}.
}
\examples{
equiv_mean_extremum(alpha = 0.01, n_sample = 6,
                    mean_qual = 100, sd_qual = 5.5, modcv = TRUE)
##
## Call:
## equiv_mean_extremum(mean_qual = 100, sd_qual = 5.5, n_sample = 6,
##     alpha = 0.01, modcv = TRUE)
##
## Modified CV used: CV* = 0.0675 ( CV = 0.055 )
##
## For alpha = 0.01 and n = 6
## ( k1 = 3.128346 and k2 = 1.044342 )
##                   Min Individual   Sample Mean
##      Thresholds:    78.88367        92.95069

}
\references{
M. G. Vangel. Lot Acceptance and Compliance Testing Using the Sample Mean
and an Extremum, Technometrics, vol. 44, no. 3. pp. 242–249. 2002.

“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.

Federal Aviation Administration, “Material Qualification and Equivalency
for Polymer Matrix Composite Material Systems,” PS-ACE 100-2002-006,
Sep. 2003.
}
\seealso{
\code{\link[=k_equiv]{k_equiv()}}

\code{\link[=calc_cv_star]{calc_cv_star()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{normalize_ply_thickness}
\alias{normalize_ply_thickness}
\title{Normalizes strength values to ply thickness}
\usage{
normalize_ply_thickness(strength, measured_thk, nom_thk)
}
\arguments{
\item{strength}{the strength to be normalized. Either a vector or a numeric}

\item{measured_thk}{the measured thickness of the samples. Must be the same
length as strength}

\item{nom_thk}{the nominal thickness. Must be a single numeric value.}
}
\value{
The normalized strength values
}
\description{
This function takes a vector of strength values and a
vector of measured thicknesses, and a nominal thickness
and returns the normalized strength.
}
\details{
It is often necessary to normalize strength values so that variation in
specimen thickness does not unnecessarily increase variation in strength.
See CMH-17-1G, or other references, for information about the cases where
normalization is appropriate.

Either cured ply thickness or laminate thickness may be used for
\code{measured_thk} and \code{nom_thk}, as long as the same decision
made for both values.

The formula applied is:
\deqn{normalized\,value = test\,value \frac{t_{measured}}{t_{nominal}}}{
normalized value = test value * t_measured / t_nominal}

If you need to normalize based on fiber volume fraction (or another method),
you will first need to calculate the nominal cured ply thickness (or laminate
thickness). Those calculations are outside the scope of this documentation.
}
\examples{
library(dplyr)

carbon.fabric.2 \%>\%
select(thickness, strength) \%>\%
  mutate(normalized_strength = normalize_ply_thickness(strength,
                                                       thickness,
                                                       0.105)) \%>\%
  head(10)

##    thickness strength normalized_strength
## 1      0.112  142.817            152.3381
## 2      0.113  135.901            146.2554
## 3      0.113  132.511            142.6071
## 4      0.112  135.586            144.6251
## 5      0.113  125.145            134.6799
## 6      0.113  135.203            145.5042
## 7      0.113  128.547            138.3411
## 8      0.113  127.709            137.4392
## 9      0.113  127.074            136.7558
## 10     0.114  126.879            137.7543


}
\references{
“Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,” SAE International,
CMH-17-1G, Mar. 2012.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{calc_cv_star}
\alias{calc_cv_star}
\title{Calculate the modified CV from the CV}
\usage{
calc_cv_star(cv)
}
\arguments{
\item{cv}{The CV to modify}
}
\value{
The value of the modified CV
}
\description{
This function calculates the modified coefficient of variation (CV)
based on a (unmodified) CV.
The modified CV is calculated based on the rules in CMH-17-1G. Those
rules are:
\itemize{
\item For CV < 4\\\%, CV* = 6\\\%
\item For 4\\\% <= CV < 8\\\%, CV* = CV / 2 + 4\\\%
\item For CV > 8\\\%, CV* = CV
}
}
\examples{
# The modified CV for values of CV smaller than 4\% is 6\%
calc_cv_star(0.01)
## [1] 0.06

# The modified CV for values of CV larger than 8\% is unchanged
calc_cv_star(0.09)
## [1] 0.09

}
\references{
"Composite Materials Handbook, Volume 1. Polymer Matrix Composites
Guideline for Characterization of Structural Materials,"
SAE International, CMH-17-1G, Mar. 2012.
}
\seealso{
\code{\link[=cv]{cv()}}
}
