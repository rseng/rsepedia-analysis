
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `dataquieR`

<!-- badges: start -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.6.0-6666ff.svg)](https://cran.r-project.org/)
[![Pipeline
Status](https://travis-ci.com/libreumg/dataquier.svg?branch=master)](https://app.travis-ci.com/gitlab/libreumg/dataquier)
[![Coverage](https://codecov.io/gl/libreumg/dataquier/branch/master/graph/badge.svg?token=79TK6GQTMG)](https://codecov.io/gl/libreumg/dataquier)
[![CRAN-Version](https://www.r-pkg.org/badges/version/dataquieR)](https://cran.r-project.org/package=dataquieR)
[![CRAN-Downloads](https://cranlogs.r-pkg.org/badges/dataquieR)](https://cran.r-project.org/package=dataquieR)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![`Lifecycle`](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![license](https://img.shields.io/badge/license-BSD_2_clause%20+%20file%20LICENSE-00be00.svg)](https://choosealicense.com/)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03093/status.svg)](https://doi.org/10.21105/joss.03093)

<!-- badges: end -->

The goal of `dataquieR` is to provide functions for assessing data
quality issues in studies, that can be used alone or in a data quality
pipeline. `dataquieR` also implements one generic pipeline producing
`flexdashboard` based HTML5 reports.

See also

[`https://dataquality.ship-med.uni-greifswald.de`](https://dataquality.ship-med.uni-greifswald.de)

------------------------------------------------------------------------

## Installation

You can install the released version of `dataquieR` from
[CRAN](https://CRAN.R-project.org/package=dataquieR) with:

``` r
install.packages("dataquieR")
```

The developer version from
[`GitLab.com`](https://gitlab.com/libreumg/dataquier) can be installed
using:

``` r
if (!requireNamespace("devtools")) {
  install.packages("devtools")
}
devtools::install_gitlab("libreumg/dataquier")
```

For examples and additional documentation, please refer to our
[website](https://dataquality.ship-med.uni-greifswald.de).

## References

-   [Software Paper](https://doi.org/10.21105/joss.03093) [![JOSS
    Article](https://joss.theoj.org/papers/10.21105/joss.03093/status.svg)](https://doi.org/10.21105/joss.03093)
-   [Data Quality Concept
    Paper](https://doi.org/10.1186/s12874-021-01252-7)
-   [Data Quality Concept and Software Web
    Site](https://dataquality.ship-med.uni-greifswald.de)

## Funding

-   German Research Foundation (DFG: `SCHM 2744/3–1`)
-   European Union’s Horizon 2020 research and innovation program (grant
    agreement No 825903.
# dataquieR (development version)
  * Fixed missing distribution plots in `dq_report`
  * Improved captions of ECDF plots
  * Updated citations
  * Better report layout
  * Fixed single variable view if no completeness has been computed
  * loading animation when report is opened in the browser but still being
    rendered

# dataquieR 1.0.9
  * Fixed bug in `sigmagap` and made missing guessing more robust.
  * Fixed checks on missing code detection failing for `logical`.
  * Fixed a damaged check for numeric threshold values in `acc_margins`.
  * Fixed wrongly named `GRADING` columns.
  * Improved parallel execution by automatic detection of cores.
  * Tidy html dependency

# dataquieR (1.0.8)
  * Removed formal arguments from `rbind.ReportSummaryTable` since these are
    not needed anyways and the inherited documentation for those arguments
    `rbind` from `base` contains an invalid URL triggering a `NOTE`.

# dataquieR (1.0.7)
  * ***Fixed bugs in example metadata.***
  * Figures now have size hints as attributes.
  * Added simple type conversion check indicator function of dimension 
    integrity, `int_datatype_matrix`.
  * Corrected some error classifications
  * `prep_study2meta` can now also convert factors to `dataquieR` compatible
    `meta_data`/`study_data`
  * Slightly improved documentation.
  * Bug fix in `com_item_missingness` for textual response variables.
  * Added new output slot with heat-map like tables. Implemented some generics
    for those.

# dataquieR (1.0.6)
  * Robustness: Ensure `DT JS` is always loaded when a dq_report report is 
    rendered
  * Bug fix: More robust handling of DECIMALS variable attribute, if
    this is delivered as a character.
  * Bug Fix: `com_segment_missingness` with 
      `strata_vars` / `group_vars` did not work
  * Bug Fix: If `label_col` was set to something else than `LABEL`, 
      `strata_vars` did not work for `com_unit_missingness`
  * More precise documentation.
  * Fixed a bug in a utility function for the univariate outliers indicator 
    function, which caused many data points flagged as outliers by the sigma-
    gap criterion.
  * Made outlier function aware of too many non-outlier points causing too
    complex graphics (e.g. pdf rendering crashes the PDF reader).
  * Fixes and small improvements in `dq_report`.
  * Switched from `cowplot` to `patchwork` in `acc_margins` yielding figures 
    that can be easier manipulated. Please note, that this change could break
    existing output manipulations, since the structure of the margins plots
    has changed internally. However, output manipulations were hardly
    possible for margins plots before, so it is unlikely, that there
    are pipelines affected.
  * More control about the output of the `acc_loess` function.
  * More robust `prep_create_meta` handling length-0 arguments by ignoring
    these variable attributes at all.
  * Added a classification system for warnings and error messages to
    distinguish errors based on mismatching variables for a function from
    other error messages.
  * https://github.com/openjournals/joss-reviews/issues/3093#issuecomment-840695360
  * Some tidy up and more tests.
  
# dataquieR 1.0.5
  * Fixed two bugs in `con_inadmissible_categorical` (one `resp_var` only and
    value-limits all the same for all `resp_vars`)
  * Changed LICENSE to BSD-2
  * Slightly updated documentation
  * Updated `README`-File

# dataquieR 1.0.4
  * Fixed CITATION, a broken reference in Rd and a problem with the vignette
    on `pandoc`-less systems
  * Improved an inaccurate argument description for multivariate outliers
  * Fixed a problem with error messages, if a `dataquieR` function was called 
    by a generated function `f` that lives in an environment 
    directly inheriting from the empty environment, e.g. 
    `environment(f) <- new.env(parent = emptyenv())`.
  * Marked some examples as `dontrun`, because they sometimes caused `NOTE`s
    on `rhub`.

# dataquieR 1.0.3
  * Addressed all comments by the CRAN reviewers, thank you.

# dataquieR 1.0.2
  * Bug Fix: If an empty data frame was delivered in the `SummaryTable` entry 
    of a result within a `dq_report` output, the `summary` and also 
    `print` generic did not work on the report.
  
# dataquieR 1.0.1
  * Skipping some of the slower tests on CRAN now. On my local system,
   a full `devtools::check(cran = TRUE, env_vars = c(NOT_CRAN = "false"))`
   takes 2:22 minutes now.
  
# dataquieR (1.0.0)
  * Initial CRAN release candidate
## Response to your mail

Dear CRAN team,

thank you for your checking. You asked us to submit an updated version:

> Dear maintainer,
>
> Please see the problems shown on
> <https://cran.r-project.org/web/checks/check_results_dataquieR.html>.
>
> Please correct before 2021-09-08 to safely retain your package on CRAN.
>
> Do remember to look at the 'Additional issues'.
>
> The CRAN Team

Thank you for your mail. We have addressed the issue as described below.

## CRAN had stated the following problems before:

```
Version: 1.0.8
Check: tests
Result: ERROR
     Running ‘testthat.R’ [47s/71s]
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
     * _snaps/int_datatype_matrix/intdtv30000.svg
     * _snaps/int_datatype_matrix/intdtv40000.svg
     * _snaps/int_datatype_matrix/intdtv50000.svg
     * _snaps/int_datatype_matrix/integrity-datatype.svg
     * _snaps/print/app-ex-repsumtab.svg
     * _snaps/print/im-empty-repsumtab.svg
     * _snaps/print/im-ex1-repsumtab.svg
     * _snaps/print/im-ex2-repsumtab.svg
     * _snaps/pro_applicability_matrix/appmatrix-plot-for-segment-v10000-ok.svg
     * _snaps/pro_applicability_matrix/appmatrix-plot-ok.svg
     * _snaps/util_heatmap_1th/util-heatmap-1th-1.svg
     * _snaps/util_heatmap_1th/util-heatmap-1th-2.svg
     * _snaps/util_heatmap_1th/util-heatmap-1th-3.svg
     Error: Test failures
     Execution halted
Flavor: r-release-macos-arm64
```

Also:

```
M1mac noLD
```

Both platforms show the same problem:

```
── Failure (test-prep_study2meta.R:745:3): prep_study2meta works on study_data ──
  `guessed_meta_data` not equal to data.frame(...).
  Component "MISSING_LIST": 1 string mismatch
```

We have addressed these no-long-double related problems and made the heuristics
function `prep_study2meta` more robust in case of lower floating point accuracy.

## Test environments and R CMD check results

### Our developers checked on their local systems:

- local R installation, R 4.1.0, macOS 10.15.7
  `devtools::check(cran = TRUE)`
  
    `0 errors ✓ | 0 warnings ✓ | 0 notes ✓`
    
### We have run checks also on Apple Silicon at https://console.scaleway.com/
  `devtools::check(cran = TRUE)`
  
    `0 errors ✓ | 0 warnings ✓ | 0 notes ✓`

    
### win-builder, CI/CD services:

- Ubuntu on travis-ci, R 4.0.2 --> Ok.
- gitlab our own ci/cd image  --> Ok.
  - debian:bullseye and R 4.0.4 (x86_64-pc-linux-gnu (64-bit))
- gitlab our own ci/cd no-long-double image  --> Ok.
- win-builder (devel, stable, oldstable) --> (Ok, Ok, Ok)
- `devtools::check_rhub` --> All three platforms ("windows-x86_64-devel",
  "ubuntu-gcc-release", and "fedora-clang-devel") Ok.
- r-hub debian-gcc-devel-nold -- Ok.
- r-hub solaris-x86-patched -- Ok.
