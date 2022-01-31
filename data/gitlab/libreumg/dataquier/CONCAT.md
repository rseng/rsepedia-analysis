
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

# `dataquieR`

<!-- badges: start -->
```{r, echo = FALSE}
dep <- as.vector(read.dcf('DESCRIPTION')[, 'Depends'])
m <- regexpr('R *\\(>= \\d+.\\d+.\\d+\\)', dep)
rm <- regmatches(dep, m)
rvers <- gsub('.*(\\d+.\\d+.\\d+).*', '\\1', rm)

minrvers <- paste0("https://img.shields.io/badge/R%3E%3D-", rvers, "-6666ff.svg")

lic <- "https://img.shields.io/badge/license-BSD_2_clause + file LICENSE-00be00.svg"
```
[![minimal R version](`r minrvers`)](https://cran.r-project.org/)
[![Pipeline Status](https://travis-ci.com/libreumg/dataquier.svg?branch=master)](https://app.travis-ci.com/gitlab/libreumg/dataquier)
[![Coverage](https://codecov.io/gl/libreumg/dataquier/branch/master/graph/badge.svg?token=79TK6GQTMG)](https://codecov.io/gl/libreumg/dataquier)
[![CRAN-Version](https://www.r-pkg.org/badges/version/dataquieR)](https://cran.r-project.org/package=dataquieR)
[![CRAN-Downloads](https://cranlogs.r-pkg.org/badges/dataquieR)](https://cran.r-project.org/package=dataquieR)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![`Lifecycle`](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![license](`r lic`)](https://choosealicense.com/)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03093/status.svg)](https://doi.org/10.21105/joss.03093)

<!-- badges: end -->

The goal of `dataquieR` is to provide functions for assessing data quality 
issues in studies, that can be used alone or in a data quality pipeline.
`dataquieR` also implements one generic pipeline producing 
`flexdashboard` based HTML5 reports.

See also

[`https://dataquality.ship-med.uni-greifswald.de`](https://dataquality.ship-med.uni-greifswald.de)

-------------

## Installation

You can install the released version of `dataquieR` from 
[CRAN](https://CRAN.R-project.org/package=dataquieR) with:

```r
install.packages("dataquieR")
```

The developer version from [`GitLab.com`](https://gitlab.com/libreumg/dataquier) 
can be installed using:

```r
if (!requireNamespace("devtools")) {
  install.packages("devtools")
}
devtools::install_gitlab("libreumg/dataquier")
```

For examples and additional documentation, please refer to our
[website](https://dataquality.ship-med.uni-greifswald.de).

## References

- [Software Paper](https://doi.org/10.21105/joss.03093) [![JOSS Article](https://joss.theoj.org/papers/10.21105/joss.03093/status.svg)](https://doi.org/10.21105/joss.03093)
- [Data Quality Concept Paper](https://doi.org/10.1186/s12874-021-01252-7)
- [Data Quality Concept and Software Web Site](https://dataquality.ship-med.uni-greifswald.de)

## Funding

- German Research Foundation (DFG: `SCHM 2744/3–1`)
- European Union’s Horizon 2020 research and innovation program (grant agreement No 825903.
---
title: "dataquieR report"
author: "`r Sys.info()[['user']]`"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  flexdashboard::flex_dashboard:
    storyboard: false
    vertical_layout: "scroll"
    orientation: "columns"
  html_document:
    toc: true
    toc_depth: 2
---

```{css}
/*.chart-wrapper { overflow-x: scroll; }*/
```

```{r echo=FALSE, include=FALSE}
knitr::knit_meta_add(list(rmarkdown::html_dependency_font_awesome()))
```


```{r setup, include=FALSE}
if (!exists("params")) { params <- list() }
params$debug <- FALSE
if (params$debug) {
  template <- "default"
  packageName <- "dataquieR"
  chunk_echo <- FALSE
  chunk_warning <- !FALSE
  chunk_error <- !FALSE
  chunk_message <- FALSE
  progress <- function(...) { invisible(NULL) }
  load("/tmp/report.RData")
}

variables <- report$app_mat$SummaryTable$Variables

df_report <- as.data.frame(report)
if (nrow(df_report) * ncol(df_report) == 0)
  stop("Report is empty, no results at all.")
library(DT)
datatable(data.frame(a = 1, b = Sys.time(), c = "test", d = 1.1, 
                     e = as.factor("test"), stringsAsFactors = FALSE),
          filter = "top")
knitr::opts_chunk$set(echo = chunk_echo, warning = chunk_warning, error = chunk_error, message = chunk_message)
cache <- new.env(parent = emptyenv())
process_result_chunk <- function(iform, outputSlot, title, subreport, variable, by_indicator = FALSE) {
    subreport <- match.arg(subreport, c("", "mv", "table", "plotlist_not_mv", "meta_miss"))
    if (nchar(subreport) > 0)
      subreport <- paste0("_", subreport)
    current_results <- subset(df_report, implementationform == iform, 'results', drop = TRUE)
    resp_vars <- subset(df_report, implementationform == iform, 'resp_vars', drop = TRUE)
    
    child_hash <-
      digest::digest(list(
        iform, outputSlot, title, subreport, variable, by_indicator
      ))
    
    if (exists(child_hash, cache)) {
      return(get(child_hash, cache))
    }
    
    if (by_indicator && subreport == "_table" && outputSlot == "SummaryTable") {
      chunk <- subset(current_results, resp_vars == variable, outputSlot)
      class(chunk) <- union(class(chunk), "SummaryTable_by_indicator")
    } else {
      chunk <- suppressWarnings(suppressMessages(try(
        knitr::knit_child(
          system.file(sprintf("%s_subreport%s.Rmd", template, subreport), package = packageName), 
          envir = environment(), 
          quiet = TRUE,
          options = list( # https://rdrr.io/github/rubenarslan/formr/man/asis_knit_child.html
             fig.path = paste0(knitr::opts_chunk$get("fig.path"), child_hash, "-"), 
             cache.path = paste0(knitr::opts_chunk$get("cache.path"), child_hash, "-")
          )
        ), 
      silent = TRUE)))
    }
    
    
    assign(child_hash, chunk, cache)

    return(chunk)
}
call_plan <- dplyr::tribble(
  ~ iform, ~ outputSlot, ~ title, ~ subreport,
#  "con_inadmissible_categorical", "SummaryTable", "### Inadmissible Categorical Values of %s", "table",
  "con_limit_deviations", "SummaryPlotList", "### Hard Limit Deviations Plots for %s",            "mv",
  "con_detection_limits", "SummaryPlotList", "### Detection Limit Deviations Plots for %s",            "mv",
#  "con_contradictions", "SummaryPlot", "### Contradictions Affecting %s", NULL,
  "acc_distributions",    "SummaryPlotList", "### Distribution Plots for %s",                "mv",
#  "acc_distributions",    "SummaryPlotList", "### Distribution Plots for %s",                "plotlist_not_mv",
  "acc_shape_or_scale",   "SummaryPlot",     "### Distribution Shape-or-Scale Plots for %s", NULL,
  "acc_margins",          "SummaryPlot",     "### Marginal Means Plots for %s",              NULL,
  "acc_loess",            "SummaryPlotList", "### LOESS-smoothed Time Course Plots for %s",  "mv",
#  "com_unit_missingness", "SummaryPlot", "### Unit Missingngess", NULL,
#  "com_segment_missingness", "SummaryPlot", "### Segment Missingness", NULL,
#  "com_item_missingness", "SummaryPlot", "### Item Missingngess", NULL,
  "acc_univariate_outlier", "SummaryPlotList", "### Univariate Outlier of %s", "mv",
  "acc_varcomp", "SummaryTable", "### Variance Components of %s", "table",
  "acc_end_digits", "SummaryPlot", "### End digit preferences of %s", NULL,
#  "acc_multivariate_outlier", "SummaryPlot", "### acc_multivariate_outlier", NULL,
)
if (params$debug) {
  call_plan <- dplyr::tribble(
    ~ iform, ~ outputSlot, ~ title, ~ subreport,
  "acc_varcomp", "SummaryTable", "### Variance Components of %s", "table",
  )
}
invisible(NULL)
```

# Overview {data-navmenu="General information"}

> Analysing a data set with `r nrow(report$study_data)` observations and `r ncol(report$study_data)` variables. 
> The metadata table comprises `r ncol(report$meta_data)` attributes of `r nrow(report$meta_data)` variables.

```{r summary}
s <- summary(report)
sv <- s$StudyVariable
for (cl in colnames(s)) {
  s[[cl]] <- htmltools::htmlEscape(s[[cl]])
}
if (nrow(s) > 0)
  s$StudyVariable <- paste0('<a href="#', sv, '">', s$StudyVariable, '</a>')
s[is.na(s)] <- "N/A"
datatable(s, filter = "top", escape = FALSE, options = list(scrollCollapse = TRUE, scrollY = "75vh"))
```

# Applicability {data-navmenu="General information"}
```{r app_mat, fig.height = 12, fig.width = 8}
report$app_mat$ApplicabilityPlot
```

# Unit missingness {data-navmenu="Completeness"}

```{r unit_miss_tab, fig.height = 12, fig.width = 8}
if (length(report$long_format$com_unit_missingness$results) >= 1)
  datatable(as.data.frame(t(report$long_format$com_unit_missingness$results[[1]]$SummaryData), stringsAsFactors = FALSE), filter = "top",
            options = list(scrollCollapse = TRUE, scrollY = "75vh"))
```

# Segment missingness figure {data-navmenu="Completeness"}

```{r segment_miss_fig, fig.height = 12, fig.width = 8, error=TRUE}
if ((!is.null(report$long_format$com_segment_missingness)) && 
  (length(report$long_format$com_segment_missingness$results) == 1))
  if (inherits(report$long_format$com_segment_missingness$results[[1]], "try-error")) {
    stop(report$long_format$com_segment_missingness$results[[1]])
  } else {
    report$long_format$com_segment_missingness$results[[1]]$SummaryPlot
  }
```

# Segment missingness table {data-navmenu="Completeness"}

```{r segment_miss_tab, fig.height = 12, fig.width = 8, error=TRUE}
if (inherits(report$long_format$com_segment_missingness$results[[1]], "try-error")) {
  stop(report$long_format$com_segment_missingness$results[[1]])
} else {
  datatable(report$long_format$com_segment_missingness$results[[1]]$SummaryData, filter = "top", options = list(scrollCollapse = TRUE, scrollY = "75vh"))
}
```

# Item missingness figure {data-navmenu="Completeness"}

```{r item_miss_fig, fig.height = 12, fig.width = 8, error=TRUE}
if (length(report$long_format$com_item_missingness$results) == 1) {
  if (inherits(report$long_format$com_item_missingness$results[[1]], "try-error")) {
    stop(report$long_format$com_item_missingness$results[[1]])
  } else {
    report$long_format$com_item_missingness$results[[1]]$SummaryPlot
  }
}
```

# Item missingness table {data-navmenu="Completeness"}

```{r item_miss_tab, fig.height = 12, fig.width = 8, error=TRUE}
if (inherits(report$long_format$com_item_missingness$results[[1]], "try-error")) {
  stop(report$long_format$com_item_missingness$results[[1]])
} else {
  datatable(report$long_format$com_item_missingness$results[[1]]$SummaryTable, filter = "top", options = list(scrollCollapse = TRUE, scrollY = "75vh"))
}
```

# Inadmissible categorical values {data-navmenu="Consistency"}

```{r iac_tab, fig.height = 12, fig.width = 8, error=TRUE}
if (inherits(report$long_format$con_inadmissible_categorical$results[[1]], "try-error")) {
  stop(report$long_format$con_inadmissible_categorical$results[[1]])
} else {
  datatable(report$long_format$con_inadmissible_categorical$results[[1]]$SummaryTable, filter = "top", options = list(scrollCollapse = TRUE, scrollY = "75vh"))
}
```

# Contradictions figure {data-navmenu="Consistency"}
```{r contra_fig, fig.height = 12, fig.width = 8}
if (length(report$long_format$con_contradictions$results) == 1)
  if (inherits(report$long_format$con_contradictions$results[[1]], "try-error")) {
    conditionMessage(attr(report$long_format$con_contradictions$results[[1]], "condition"))
  } else {
    report$long_format$con_contradictions$results[[1]]$SummaryPlot
  }
```

# Contradictions table {data-navmenu="Consistency"}

```{r contra_tab, fig.height = 12, fig.width = 8}
progress(30)
if (inherits(report$long_format$con_contradictions$results[[1]], "try-error")) {
  conditionMessage(attr(report$long_format$con_contradictions$results[[1]],
                        "condition"))
} else {
  datatable(report$long_format$con_contradictions$results[[1]]$SummaryData, 
            filter = "top", options = list(scrollCollapse = TRUE, scrollY = "75vh"))
}
```

```{r con_output1, echo=FALSE, results='asis', eval = TRUE}
cp <- call_plan[startsWith(call_plan$iform, "con_"), , FALSE]
by_indicator <- dataquieR::prep_pmap(cp, function(...) {
  single_results <- 
    lapply(setNames(nm = unique(as.character(variables))), function(variable) {
    process_result_chunk(..., variable = variable, by_indicator = TRUE)
  })
  errors <- vapply(single_results, function(x) {
    inherits(x, "try-error") || inherits(x, "error")
  } , FUN.VALUE =  logical(1))
  if (any(errors)) {
    if (chunk_error) {
      invisible(lapply(single_results[errors], function(e) try(stop(e))))
    }
    single_results <- single_results[!errors]
  }
  if (all(vapply(single_results, inherits, "SummaryTable_by_indicator", 
                 FUN.VALUE = logical(1)))
      ) {
    single_results <- 
      single_results[vapply(single_results, length, FUN.VALUE = integer(1)) > 0]
    sr <- unlist(single_results, recursive = FALSE)
    single_errors <- vapply(sr, function(x) {
      w <- attr(x, "warning")
      if (chunk_warning) invisible(lapply(w, function(y) try(warning(y))))
      m <- attr(x, "message")
      if (chunk_message) invisible(lapply(m, function(y) try(message(y))))
      e <- attr(x, "error")
      if (chunk_error) invisible(lapply(e, function(y) try(stop(y))))
      length(e) > 0
    }, FUN.VALUE = logical(1))
    sr <- sr[!single_errors]
    oslo <- c(...)["outputSlot"] # oslo <- "SummaryTable" # for now only this
    if (all(vapply(lapply(sr, names), `%in%`, x = oslo, FUN.VALUE = logical(1)))) {
      sr <- lapply(sr, `[[`, oslo)
      srdf <- try(do.call(rbind.data.frame, sr), silent = TRUE)
      if (is.data.frame(srdf)) {
        single_results <- as.character(htmltools::tagList(DT::datatable(srdf, filter = "top", options = list(scrollCollapse = TRUE, scrollY = "75vh"))))
      }
    }
  }
  iform <- c(...)["iform"]
  title <- c(...)["title"]
  title <- gsub("^#+ ", "", title)
  title <- gsub(" [a-z]+ %s$", "", title)
  r <- c(sprintf('\n# %s {#%s_con_mv data-navmenu="Consistency"}', title, 
                 iform), 
         single_results) 
  
  r[trimws(r) == ""] <- NULL
  progress(40)
  if (length(r) > 1) {
    r
  } else {
    NULL
  }
})
cat(unlist(by_indicator), sep = '\n')
names(by_indicator) <- cp$iform
if (chunk_error) {
  errors <- unlist(unlist(unlist(lapply(by_indicator, lapply, lapply, attr, "error"), recursive = FALSE), recursive = FALSE), recursive = FALSE)
  invisible(mapply(function(e, n) {
    n <- strsplit(n, ".", fixed = TRUE)[[1]]
    if (length(n) < 1) {
      n <- c("f", "v")
    } else if (length(n) < 2) {
      n <- c(n[[1]], n[[1]])
    }
    cond <- simpleError(conditionMessage(e), call(n[[1]], n[[2]]));
    try(stop(cond))
  }, e = errors, n = names(errors)))
}
if (chunk_warning) {
  errors <- unlist(unlist(unlist(lapply(by_indicator, lapply, lapply, attr, "warning"), recursive = FALSE), recursive = FALSE), recursive = FALSE)
  invisible(mapply(function(e, n) {
    n <- strsplit(n, ".", fixed = TRUE)[[1]]
    if (length(n) < 1) {
      n <- c("f", "v")
    } else if (length(n) < 2) {
      n <- c(n[[1]], n[[1]])
    }
    cond <- simpleWarning(conditionMessage(e), call(n[[1]], n[[2]])); 
    try(warning(cond))
  }, e = errors, n = names(errors)))
}
if (chunk_message) {
  errors <- unlist(unlist(unlist(lapply(by_indicator, lapply, lapply, attr, "message"), recursive = FALSE), recursive = FALSE), recursive = FALSE)
  invisible(mapply(function(e, n) {
    n <- strsplit(n, ".", fixed = TRUE)[[1]]
    if (length(n) < 1) {
      n <- c("f", "v")
    } else if (length(n) < 2) {
      n <- c(n[[1]], n[[1]])
    }
    cond <- simpleMessage(conditionMessage(e), call(n[[1]], n[[2]])); 
    try(message(cond))
  }, e = errors, n = names(errors)))
}
```
```{r con_output2, echo=FALSE, results='asis', eval = TRUE}
res <- lapply(unique(as.character(variables)), function(variable) {
    cp <- call_plan[startsWith(call_plan$iform, "con_"), , FALSE]
    cp$variable <- variable
    single_results <- dataquieR::prep_pmap(cp, process_result_chunk)
    errors <- vapply(single_results, function(x) {
      inherits(x, "try-error") || inherits(x, "error")
    } , FUN.VALUE =  logical(1))
    if (any(errors)) {
      if (chunk_error) {
        invisible(lapply(single_results[errors], function(e) try(stop(e))))
      }
      single_results <- single_results[!errors]
    }
    r <- c(sprintf('\n# %s {#%s_con data-navmenu="Consistency"}', dQuote(variable), variable), 
           single_results) 

    r[trimws(r) == ""] <- NULL
    if (length(r) > 1) {
      r
    } else {
      NULL
    }
  })
progress(50)
cat(unlist(res), sep = '\n')
```





```{r acc_output1, echo=FALSE, results='asis', eval = TRUE}
cp <- call_plan[startsWith(call_plan$iform, "acc_"), , FALSE]
by_indicator <- dataquieR::prep_pmap(cp, function(...) {
  single_results <- 
    lapply(setNames(nm = unique(as.character(variables))), function(variable) {
    process_result_chunk(..., variable = variable, by_indicator = TRUE)
  })
  errors <- vapply(single_results, function(x) {
    inherits(x, "try-error") || inherits(x, "error")
  } , FUN.VALUE =  logical(1))
  if (any(errors)) {
    if (chunk_error) {
      invisible(lapply(single_results[errors], function(e) try(stop(e))))
    }
    single_results <- single_results[!errors]
  }
  if (all(vapply(single_results, inherits, "SummaryTable_by_indicator", 
                 FUN.VALUE = logical(1)))
      ) {
    single_results <- 
      single_results[vapply(single_results, length, FUN.VALUE = integer(1)) > 0]
    sr <- unlist(single_results, recursive = FALSE)
    single_errors <- vapply(sr, function(x) {
      w <- attr(x, "warning")
      if (chunk_warning) invisible(lapply(w, function(y) try(warning(y))))
      m <- attr(x, "message")
      if (chunk_message) invisible(lapply(m, function(y) try(message(y))))
      e <- attr(x, "error")
      if (chunk_error) invisible(lapply(e, function(y) try(stop(y))))
      length(e) > 0
    }, FUN.VALUE = logical(1))
    sr <- sr[!single_errors]
    oslo <- c(...)["outputSlot"] # oslo <- "SummaryTable" # for now only this
    if (all(vapply(lapply(sr, names), `%in%`, x = oslo, FUN.VALUE = logical(1)))) {
      sr <- lapply(sr, `[[`, oslo)
      srdf <- try(do.call(rbind.data.frame, sr), silent = TRUE)
      if (is.data.frame(srdf)) {
        single_results <- as.character(htmltools::tagList(DT::datatable(srdf, filter = "top", options = list(scrollCollapse = TRUE, scrollY = "75vh"))))
      }
    }
  }
  iform <- c(...)["iform"]
  title <- c(...)["title"]
  title <- gsub("^#+ ", "", title)
  title <- gsub(" [a-z]+ %s$", "", title)
  r <- c(sprintf('\n# %s {#%s_acc_mv data-navmenu="Accuracy"}', title, iform), 
         single_results) 
  
  r[trimws(r) == ""] <- NULL
  if (length(r) > 1) {
    r
  } else {
    NULL
  }
})
progress(65)
cat(unlist(by_indicator), sep = '\n')
names(by_indicator) <- cp$iform
if (chunk_error) {
  errors <- unlist(unlist(unlist(lapply(by_indicator, lapply, lapply, attr, "error"), recursive = FALSE), recursive = FALSE), recursive = FALSE)
  invisible(mapply(function(e, n) {
    n <- strsplit(n, ".", fixed = TRUE)[[1]]
    if (length(n) < 1) {
      n <- c("f", "v")
    } else if (length(n) < 2) {
      n <- c(n[[1]], n[[1]])
    }
    cond <- simpleError(conditionMessage(e), call(n[[1]], n[[2]])); 
    try(stop(cond))
  }, e = errors, n = names(errors)))
}
if (chunk_warning) {
  errors <- unlist(unlist(unlist(lapply(by_indicator, lapply, lapply, attr, "warning"), recursive = FALSE), recursive = FALSE), recursive = FALSE)
  invisible(mapply(function(e, n) {
    n <- strsplit(n, ".", fixed = TRUE)[[1]]
    if (length(n) < 1) {
      n <- c("f", "v")
    } else if (length(n) < 2) {
      n <- c(n[[1]], n[[1]])
    }
    cond <- simpleWarning(conditionMessage(e), call(n[[1]], n[[2]])); 
    try(warning(cond))
  }, e = errors, n = names(errors)))
}
if (chunk_message) {
  errors <- unlist(unlist(unlist(lapply(by_indicator, lapply, lapply, attr, "message"), recursive = FALSE), recursive = FALSE), recursive = FALSE)
  invisible(mapply(function(e, n) {
    n <- strsplit(n, ".", fixed = TRUE)[[1]]
    if (length(n) < 1) {
      n <- c("f", "v")
    } else if (length(n) < 2) {
      n <- c(n[[1]], n[[1]])
    }
    cond <- simpleMessage(conditionMessage(e), call(n[[1]], n[[2]])); 
    try(message(cond))
  }, e = errors, n = names(errors)))
}
```

```{r acc_output2, echo=FALSE, results='asis', eval = TRUE}
res <- lapply(unique(as.character(variables)), function(variable) {
    cp <- call_plan[startsWith(call_plan$iform, "acc_"), , FALSE]
    cp$variable <- variable
    single_results <- dataquieR::prep_pmap(cp, process_result_chunk)
    errors <- vapply(single_results, function(x) {
      inherits(x, "try-error") || inherits(x, "error")
    } , FUN.VALUE =  logical(1))
    if (any(errors)) {
      if (chunk_error) {
        invisible(lapply(single_results[errors], function(e) try(stop(e))))
      }
      single_results <- single_results[!errors]
    }
    r <- c(sprintf('\n# %s {#%s_acc data-navmenu="Accuracy"}', dQuote(variable), variable), 
           single_results) 

    r[trimws(r) == ""] <- NULL
    if (length(r) > 1) {
      r
    } else {
      NULL
    }
  })
progress(85)
cat(unlist(res), sep = '\n')
```

```{r main_output, echo=FALSE, results='asis', eval = TRUE}
res <- lapply(unique(as.character(variables)), function(variable) {
    cp <- call_plan
    cp$variable <- variable
    single_results <- dataquieR::prep_pmap(cp, process_result_chunk)
    errors <- vapply(single_results, function(x) {
      inherits(x, "try-error") || inherits(x, "error")
    } , FUN.VALUE =  logical(1))
    if (any(errors)) {
      if (chunk_error) {
        invisible(lapply(single_results[errors], function(e) try(stop(e))))
      }
      single_results <- single_results[!errors]
    }
    metamiss <- try(
      process_result_chunk("com_item_missingness", "", "### General Information about %s", "meta_miss", variable),
      silent = TRUE)
    if (inherits(metamiss, "try-error"))
      metamiss <- ""
    r <- c(list(sprintf('\n# %s {#%s data-navmenu="Single Variables"}', dQuote(variable), variable)), 
           metamiss,
           single_results) 

    r[trimws(r) == ""] <- NULL
    if (length(r) > 1) {
      r
    } else {
      NULL
    }
  })
progress(99)
cat(unlist(res), sep = '\n')
```
```{r echo=FALSE}
#if (any(grepl("Limit", title))) {
#  save(resp_vars, outputSlot, current_results, file = "/tmp/con.RData")
#}
if (is.list(resp_vars) && length(resp_vars) == 1 && is.character(resp_vars[[1]])) {
  resp_vars <- resp_vars[[1]]
}
if ((all(!is.na(resp_vars) && (length(resp_vars) == 1))) || outputSlot == "SummaryPlot") {
  stop("Plot result")
}
if (outputSlot == "SummaryTable") {
  stop("Not a table result")
}
all_res <- lapply(current_results, `[[`, outputSlot)
if (chunk_error) {
  errors <- lapply(current_results, attr, "error")
  invisible(lapply(errors, function(e) {
    if (length(e) == 1) if (inherits(e[[1]], "error")) {
      cond <- simpleError(conditionMessage(e[[1]]), call(iform, variable)); 
      try(stop(cond))
    }
  }))
}
if (chunk_warning) {
  errors <- lapply(current_results, attr, "warning")
  invisible(lapply(errors, function(e) {
    if (length(e) == 1) if (inherits(e[[1]], "warning")) {
      cond <- simpleWarning(conditionMessage(e[[1]]), call(iform, variable)); 
      try(warning(cond))
    }
  }))
}
if (chunk_message) {
  errors <- lapply(current_results, attr, "message")
  invisible(lapply(errors, function(e) {
    if (length(e) == 1) if (inherits(e[[1]], "message")) {
      cond <- simpleMessage(conditionMessage(e[[1]]), call(iform, variable)); 
      try(message(cond))
    }
  }))
}
have_a_result <- any(unlist(lapply(all_res, function(r) variable %in% names(r))))
if (!have_a_result) {
  stop("No results here")
}
```
```{r results='asis', eval=have_a_result}
cat(sprintf(title, dQuote(variable)))
cat("\n")
```
```{r echo=FALSE, eval=have_a_result}
# if (variable == "SBP_0" && iform == "acc_univariate_outlier") save(all_res, variable, file = "/tmp/yyy.RData")
invisible(lapply(all_res, function(r) {
  if (!is.null(r[[variable]])) print(r[[variable]])
  }))
```
```{r include=FALSE}
if ((all(is.na(variable) || (length(variable) != 1))) && outputSlot == "SummaryPlotList") {
  stop("PlotList result")
}
if (outputSlot == "SummaryTable") {
  stop("Not a table result")
}
all_res <- current_results[resp_vars == variable]
if (!chunk_error) {
  all_res <- all_res[(vapply(all_res, length, FUN.VALUE = integer(1)) > 0)]
}
have_a_result <- length(all_res) > 0
if (!have_a_result) {
  stop("No results here")
}
```
```{r results='asis', eval=have_a_result}
cat(sprintf(title, dQuote(variable)))
cat("\n")
```
```{r echo=FALSE, eval=have_a_result}
invisible(lapply(all_res, print, slot = outputSlot))
```
```{r include=FALSE}
library(ggplot2)
if ((all(is.na(variable) || (length(variable) != 1))) || outputSlot != "") {
  stop("Not a meta_miss result")
}
all_res <- current_results[resp_vars == variable]
if (!chunk_error) {
  all_res <- all_res[(vapply(all_res, length, FUN.VALUE = integer(1)) > 0)]
}
have_a_result <- length(all_res) > 0
```
::: {style="height: 3em;"}
:::
```{r results='asis', eval=TRUE}
cat(sprintf(title, dQuote(variable)))
cat("\n")
cat("\n#### Metadata\n")
cat("\n")
md <- report$meta_data[report$meta_data[[report$label_col]] == variable, , drop = FALSE]
md <- unlist(md[1, , drop = TRUE])
md <- cbind(`Variable Attribute` = names(md), `Attribute Value` = md)
rownames(md) <- NULL
DT::datatable(md, filter = "top", options = list(pageLength = 4,
                                                 scrollCollapse = TRUE, scrollY = "75vh")) # nrow(md)))
```

#### Completeness

```{r echo=FALSE, eval=TRUE||have_a_result}
DT::datatable(subset(current_results[[1]]$SummaryTable, 
                     Variables == variable,
                     drop = FALSE), filter = "top", 
              options = list(scrollCollapse = TRUE, scrollY = "75vh"))
# invisible(lapply(all_res, print, slot = "SummaryPlot"))
```

```{r echo=FALSE, eval=TRUE||have_a_result}
bp_data <- 
  current_results[[1]]$SummaryPlot$data[
    current_results[[1]]$SummaryPlot$data$Var2 == variable, 
    c("Var1", "Freq")]

n_obs <- subset(current_results[[1]]$SummaryTable, 
                Variables == variable, 
                "Observations N", 
                drop = TRUE)

bp_data$Pct <- bp_data$Freq / n_obs

my_cols <- c("#7f0000", "#b30000", "#d7301f", "#ef6548", "#fc8d59",
             "#fdbb84", "#fdd49e", "#fee8c8", "#2166AC")


ggplot(bp_data, aes(x = Var1, y = Pct, fill = Pct)) + 
  geom_col(show.legend = FALSE) + 
  scale_y_continuous(labels = scales::percent_format()) + 
  scale_fill_gradientn(colors = rev(my_cols)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("Missing Codes") +
  ylab("rel. Frequency") +
  ggtitle(sprintf("Completeness of %s", dQuote(variable)), "")
```

```{r include=FALSE}
if ((all(!is.na(resp_vars) && (length(resp_vars) == 1))) || outputSlot == "SummaryPlot") {
  stop("Plot result")
}
if (outputSlot == "SummaryTable") {
  stop("Not a table result")
}
all_res <- current_results[resp_vars == variable]
have_a_result <- any(vapply(all_res, function(r) {
  length(r) > 0
}, FUN.VALUE = logical(1)))
if (!have_a_result) {
  stop("No results here")
}
```
```{r results='asis', eval=have_a_result}
cat(sprintf(title, dQuote(variable)))
cat("\n")
```
```{r echo=FALSE, eval=have_a_result}
# if (variable == "SBP_0" && iform == "acc_loess") save(all_res, variable, file = "/tmp/yyy.RData")
invisible(lapply(all_res, function(r) {
  lapply(r, function(v) lapply(v, print))
}))
```
```{r include=FALSE}
if (outputSlot != "SummaryTable") {
  stop("Not a table result")
}
all_res <- current_results[resp_vars == variable]
if (!chunk_error) {
  all_res <- all_res[(vapply(all_res, length, FUN.VALUE = integer(1)) > 0)]
}
have_a_result <- length(all_res) > 0
if (!have_a_result) {
  stop("No results here")
}
```
```{r results='asis', eval=have_a_result}
cat(sprintf(title, dQuote(variable)))
cat("\n")
```
```{r echo=FALSE, eval=have_a_result, results='asis'}
for (x in all_res) {
  if (length(attr(x, "message")) > 0) {
    for (m in attr(x, "message"))
      message(m)
  }
  if (length(attr(x, "warning")) > 0) {
    for (w in attr(x, "warning"))
      warning(w)
  }
  error_shown <- FALSE
  if (length(attr(x, "error")) > 0) {
    e <- attr(x, "error")[[1]]
    try(stop(e))
    error_shown <- TRUE
  }
  attr(x, "message") <- NULL
  attr(x, "warning") <- NULL
  attr(x, "error") <- NULL
  class(x) <- setdiff(class(x), 'dataquieR_result')
  if (is.data.frame(x[[outputSlot]]) && nrow(x[[outputSlot]]) > 0 && ncol(x[[outputSlot]]) > 0) {
    print( htmltools::tagList(DT::datatable(x[[outputSlot]], filter = "top",
                                            options = list(scrollCollapse = TRUE, scrollY = "75vh"))) )
  }
}
```
---
title: "dataquieR example report"
author: "Adrian Richter, Stephan Struckmann, Carsten Schmidt"
output:
  rmarkdown::html_vignette:
    css: dfg_qs_style.css
    toc: TRUE
    toc_depth: 3
vignette: >
  %\VignetteIndexEntry{dataquieR example report}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r include=FALSE}
library(knitr)
library(DT)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
if (rmarkdown::pandoc_available(version = "1.12.3")) {
  knit_print.data.frame = function(x, ...) {
    knit_print(DT::datatable(head(x, 10)), ...)
  }
  registerS3method("knit_print", "data.frame", knit_print.data.frame)
}
library(dataquieR)
```

# Preface

This is a brief example report using `dataquieR`'s functions. For a longer and
better elaborated example, please also consider our
[online example with data from SHIP](https://dataquality.ship-med.uni-greifswald.de/tutorials.html#Advanced_Report_SHIP_Example).

# INTEGRITY

## Study data

```{r echo=TRUE, warning=FALSE, message=FALSE}
load(system.file("extdata", "study_data.RData", package = "dataquieR"))
sd1 <- study_data
```

The imported study data consist of:

* N = `r dim(sd1)[1]` observations and
* P = `r dim(sd1)[2]` study variables


## Metadata 


```{r echo=TRUE, warning=FALSE, message=FALSE}
load(system.file("extdata", "meta_data.RData", package = "dataquieR"))
md1 <- meta_data
```

The imported meta data provide information for:

* P = `r dim(md1)[1]` study variables and
* Q = `r dim(md1)[2]` attributes

## Applicability 

The call of this R-function requires two inputs only:

```{r message=FALSE, warning=FALSE}
appmatrix <- pro_applicability_matrix(study_data = sd1, 
                                      meta_data = md1, 
                                      label_col = LABEL)
```

Heatmap-like plot:

```{r message=FALSE, warning=FALSE, fig.height = 10, fig.width = 6}
appmatrix$ApplicabilityPlot
```

# COMPLETENESS

## Unit missingness

```{r message = FALSE, warning = FALSE}
my_unit_missings2 <- com_unit_missingness(study_data  = sd1, 
                                          meta_data   = md1, 
                                          id_vars     = c("CENTER_0", "PSEUDO_ID"), 
                                          strata_vars = "CENTER_0", 
                                          label_col   = "LABEL")
```


```{r}
my_unit_missings2$SummaryData
```

## Segment missingness

```{r message=FALSE, warning=FALSE}
MissSegs <- com_segment_missingness(study_data = sd1, 
                                    meta_data = md1, 
                                    label_col = "LABEL", 
                                    threshold_value = 5, 
                                    direction = "high",
                                    exclude_roles = c("secondary", "process"))
```

```{r message=FALSE, echo=TRUE, warning=FALSE, results = 'hide', fig.keep = 'all',  fig.align="center", fig.height = 3, fig.width = 4}
MissSegs$SummaryPlot
```

### Adding variables for stratification

For some analyses adding new and transformed variable to the study data is necessary.

```{r message=FALSE, warning=FALSE}
# use the month function of the lubridate package to extract month of exam date
require(lubridate)
# apply changes to copy of data
sd2 <- sd1
# indicate first/second half year
sd2$month <- month(sd2$v00013)
```

Static metadata of the variable must be added to the respective metadata.

```{r message=FALSE, warning=FALSE}
MD_TMP <- prep_add_to_meta(VAR_NAMES    = "month",
                           DATA_TYPE    = "integer",
                           LABEL        = "EXAM_MONTH",
                           VALUE_LABELS = "1 = January | 2 = February | 3 = March | 
                                          4 = April | 5 = May | 6 = June | 7 = July |
                                          8 = August | 9 = September | 10 = October |
                                          11 = November | 12 = December",
                           meta_data    = md1)
```

Subsequent call of the R-function may include the new variable.

```{r message=FALSE, warning=FALSE}
MissSegs <- com_segment_missingness(study_data = sd2, 
                                    meta_data = MD_TMP, 
                                    group_vars = "EXAM_MONTH", 
                                    label_col = "LABEL", 
                                    threshold_value = 1, 
                                    direction = "high",
                                    exclude_roles = c("secondary", "process"))
```

```{r message=FALSE, echo=TRUE, warning=FALSE, results = 'hide', fig.keep = 'all',  fig.align="center", fig.height = 6, fig.width = 4}
MissSegs$SummaryPlot
```

## Item missingness

The following implementation considers also labeled missing codes. The use of such a table is optional but recommended. Missing code labels used in the simulated study data are loaded as follows:

```{r message=FALSE, warning=FALSE}
code_labels <- read.csv2(system.file("extdata", 
                                     "Missing-Codes-2020.csv", 
                                     package = "dataquieR"), 
                         stringsAsFactors = FALSE, na.strings = c())
```

```{r message = FALSE, warning = FALSE}
item_miss <- com_item_missingness(study_data      = sd1, 
                                  meta_data       = meta_data, 
                                  label_col       = 'LABEL', 
                                  show_causes     = TRUE, 
                                  cause_label_df  = code_labels,
                                  include_sysmiss = TRUE, 
                                  threshold_value = 80
                                ) 
```

The function call above sets the analyses of causes for missing values to TRUE, includes system missings with an own code, and sets the threshold to 80%.

```{r message=FALSE, echo=TRUE, warning=FALSE}
item_miss$SummaryTable
```

#### Summary plot of item missingness

```{r message=FALSE, echo=TRUE, warning=FALSE, fig.height=5, fig.width = 5}
item_miss$SummaryPlot
```


# CONSISTENCY

## Limit deviations

```{r }
MyValueLimits <- con_limit_deviations(resp_vars  = NULL,
                                      label_col  = "LABEL",
                                      study_data = sd1,
                                      meta_data  = md1,
                                      limits     = "HARD_LIMITS")
```

### Summary table
```{r message=FALSE, echo=TRUE, warning=FALSE}
MyValueLimits$SummaryTable
```


### Summary plot
```{r }
# select variables with deviations
whichdeviate <- as.character(MyValueLimits$SummaryTable$Variables)[MyValueLimits$SummaryTable$GRADING == 1]
```


```{r message=FALSE, echo=TRUE, warning=FALSE, results = 'hide', fig.keep = 'all', fig.align="center", fig.height = 3, fig.width = 4}
ggpubr::ggarrange(plotlist = MyValueLimits$SummaryPlotList[whichdeviate], ncol = 2) 
```

## Inadmissible levels

```{r message=FALSE, warning=FALSE}
IAVCatAll <- con_inadmissible_categorical(study_data = sd1, 
                                          meta_data  = md1, 
                                          label_col  = "LABEL")
```

## Contradictions

```{r message=FALSE, warning=FALSE}
checks <- read.csv(system.file("extdata", 
                               "contradiction_checks.csv",
                               package = "dataquieR"), 
                   header = TRUE, sep = "#")
```

```{r }
AnyContradictions <- con_contradictions(study_data      = sd1,
                                        meta_data       = md1,
                                        label_col       = "LABEL",
                                        check_table     = checks,
                                        threshold_value = 1)
```

```{r message=FALSE, echo=TRUE, warning=FALSE}
AnyContradictions$SummaryTable
```

```{r message=FALSE, echo=TRUE, warning=FALSE, fig.height = 4, fig.width = 6}
AnyContradictions$SummaryPlot 
```

# ACCURACY

```{r echo = TRUE}
ruol <- dataquieR:::acc_robust_univariate_outlier(study_data = sd1, meta_data = md1, label_col = LABEL)

ruol$SummaryPlotList
```


```{r, fig.height = 3, fig.width = 4}
myloess <- dataquieR::acc_loess(resp_vars = "SBP_0",
                                group_vars = "USR_BP_0",
                                time_vars = "EXAM_DT_0",
                                label_col = "LABEL",
                                study_data = sd1,
                                meta_data = md1)

myloess$SummaryPlotList
```




