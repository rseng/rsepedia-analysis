
<!-- README.md is generated from README.Rmd. Please edit that file -->

# excluder

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![lifecycle](man/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/excluder)](https://cran.r-project.org/package=excluder)
<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/excluder)](https://CRAN.R-project.org/package=excluder) -->
<!-- [![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/excluder?color=orange)](https://CRAN.R-project.org/package=excluder) -->

[![R-CMD-check](https://github.com/ropensci/excluder/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/excluder/actions)
[![Codecov test
coverage](https://codecov.io/gh/jeffreyrstevens/excluder/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jeffreyrstevens/excluder?branch=main)

[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/455_status.svg)](https://github.com/ropensci/software-review/issues/455)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03893/status.svg)](https://doi.org/10.21105/joss.03893)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5648202.svg)](https://doi.org/10.5281/zenodo.5648202)
<!-- badges: end -->

The goal of [`{excluder}`](https://docs.ropensci.org/excluder/) is to
facilitate checking for, marking, and excluding rows of data frames for
common exclusion criteria. This package applies to data collected from
[Qualtrics](https://www.qualtrics.com/) surveys, and default column
names come from importing data with the
[`{qualtRics}`](https://docs.ropensci.org/qualtRics/) package.

This may be most useful for [Mechanical Turk](https://www.mturk.com/)
data to screen for duplicate entries from the same location/IP address
or entries from locations outside of the United States. But it can be
used more generally to exclude based on response durations, preview
status, progress, or screen resolution.

More details are available on the package
[website](https://docs.ropensci.org/excluder/) and the [getting started
vignette](https://docs.ropensci.org/excluder/articles/getting_started.html).

## Installation

You can install the stable released version of `{excluder}` from
[CRAN](https://cran.r-project.org/package=excluder) with:

``` r
install.packages("excluder")
```

You can install developmental versions from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/excluder")
```

## Verbs

This package provides three primary verbs:

-   `mark` functions add a new column to the original data frame that
    labels the rows meeting the exclusion criteria. This is useful to
    label the potential exclusions for future processing without
    changing the original data frame.
-   `check` functions search for the exclusion criteria and output a
    message with the number of rows meeting the criteria and a data
    frame of the rows meeting the criteria. This is useful for viewing
    the potential exclusions.
-   `exclude` functions remove rows meeting the exclusion criteria. This
    is safest to do after checking the rows to ensure the exclusions are
    correct.

## Exclusion types

This package provides seven types of exclusions based on Qualtrics
metadata. If you have ideas for other metadata exclusions, please submit
them as [issues](https://github.com/ropensci/excluder/issues). Note, the
intent of this package is not to develop functions for excluding rows
based on survey-specific data but on general, frequently used metadata.

-   `duplicates` works with rows that have duplicate IP addresses and/or
    locations (latitude/longitude).
-   `duration` works with rows whose survey completion time is too short
    and/or too long.
-   `ip` works with rows whose IP addresses are not found in the
    specified country (note: this exclusion type requires an internet
    connection to download the country’s IP ranges).
-   `location` works with rows whose latitude and longitude are not
    found in the United States.
-   `preview` works with rows that are survey previews.
-   `progress` works with rows in which the survey was not complete.
-   `resolution` works with rows whose screen resolution is not
    acceptable.

## Usage

The verbs and exclusion types combine with `_` to create the functions,
such as
[`check_duplicates()`](https://docs.ropensci.org/excluder/reference/check_duplicates.html),
[`exclude_ip()`](https://docs.ropensci.org/excluder/reference/exclude_ip.html),
and
[`mark_duration()`](https://docs.ropensci.org/excluder/reference/mark_duration.html).
Multiple functions can be linked together using the
[`{magrittr}`](https://magrittr.tidyverse.org/) pipe `%>%`. For datasets
downloaded directly from Qualtrics, use
[`remove_label_rows()`](https://docs.ropensci.org/excluder/reference/remove_label_rows.html)
to remove the first two rows of labels and convert date and numeric
columns in the metadata, and use
[`deidentify()`](https://docs.ropensci.org/excluder/reference/deidentify.html)
to remove standard Qualtrics columns with identifiable information
(e.g., IP addresses, geolocation).

### Marking

The `mark_*()` functions output the original data set with a new column
specifying rows that meet the exclusion criteria. These can be piped
together with `%>%` for multiple exclusion types.

``` r
library(excluder)
# Mark preview and short duration rows
df <- qualtrics_text %>%
  mark_preview() %>%
  mark_duration(min_duration = 200)
#> ℹ 2 rows were collected as previews. It is highly recommended to exclude these rows before further processing.
#> ℹ 23 out of 100 rows took less time than 200.
tibble::glimpse(df)
#> Rows: 100
#> Columns: 18
#> $ StartDate               <dttm> 2020-12-11 12:06:52, 2020-12-11 12:06:43, 202…
#> $ EndDate                 <dttm> 2020-12-11 12:10:30, 2020-12-11 12:11:27, 202…
#> $ Status                  <chr> "Survey Preview", "Survey Preview", "IP Addres…
#> $ IPAddress               <chr> NA, NA, "73.23.43.0", "16.140.105.0", "107.57.…
#> $ Progress                <dbl> 100, 100, 100, 100, 100, 100, 100, 100, 100, 1…
#> $ `Duration (in seconds)` <dbl> 465, 545, 651, 409, 140, 213, 177, 662, 296, 2…
#> $ Finished                <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE…
#> $ RecordedDate            <dttm> 2020-12-11 12:10:30, 2020-12-11 12:11:27, 202…
#> $ ResponseId              <chr> "R_xLWiuPaNuURSFXY", "R_Q5lqYw6emJQZx2o", "R_f…
#> $ LocationLatitude        <dbl> 29.73694, 39.74107, 34.03852, 44.96581, 27.980…
#> $ LocationLongitude       <dbl> -94.97599, -121.82490, -118.25739, -93.07187, …
#> $ UserLanguage            <chr> "EN", "EN", "EN", "EN", "EN", "EN", "EN", "EN"…
#> $ Browser                 <chr> "Chrome", "Chrome", "Chrome", "Chrome", "Chrom…
#> $ Version                 <chr> "88.0.4324.41", "88.0.4324.50", "87.0.4280.88"…
#> $ `Operating System`      <chr> "Windows NT 10.0", "Macintosh", "Windows NT 10…
#> $ Resolution              <chr> "1366x768", "1680x1050", "1366x768", "1536x864…
#> $ exclusion_preview       <chr> "preview", "preview", "", "", "", "", "", "", …
#> $ exclusion_duration      <chr> "", "", "", "", "duration_quick", "", "duratio…
```

Use the
[`unite_exclusions()`](https://docs.ropensci.org/excluder/reference/unite_exclusions.html)
function to unite all of the marked columns into a single column.

``` r
# Collapse labels for preview and short duration rows
df <- qualtrics_text %>%
  mark_preview() %>%
  mark_duration(min_duration = 200) %>%
  unite_exclusions()
#> ℹ 2 rows were collected as previews. It is highly recommended to exclude these rows before further processing.
#> ℹ 23 out of 100 rows took less time than 200.
tibble::glimpse(df)
#> Rows: 100
#> Columns: 17
#> $ StartDate               <dttm> 2020-12-11 12:06:52, 2020-12-11 12:06:43, 202…
#> $ EndDate                 <dttm> 2020-12-11 12:10:30, 2020-12-11 12:11:27, 202…
#> $ Status                  <chr> "Survey Preview", "Survey Preview", "IP Addres…
#> $ IPAddress               <chr> NA, NA, "73.23.43.0", "16.140.105.0", "107.57.…
#> $ Progress                <dbl> 100, 100, 100, 100, 100, 100, 100, 100, 100, 1…
#> $ `Duration (in seconds)` <dbl> 465, 545, 651, 409, 140, 213, 177, 662, 296, 2…
#> $ Finished                <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE…
#> $ RecordedDate            <dttm> 2020-12-11 12:10:30, 2020-12-11 12:11:27, 202…
#> $ ResponseId              <chr> "R_xLWiuPaNuURSFXY", "R_Q5lqYw6emJQZx2o", "R_f…
#> $ LocationLatitude        <dbl> 29.73694, 39.74107, 34.03852, 44.96581, 27.980…
#> $ LocationLongitude       <dbl> -94.97599, -121.82490, -118.25739, -93.07187, …
#> $ UserLanguage            <chr> "EN", "EN", "EN", "EN", "EN", "EN", "EN", "EN"…
#> $ Browser                 <chr> "Chrome", "Chrome", "Chrome", "Chrome", "Chrom…
#> $ Version                 <chr> "88.0.4324.41", "88.0.4324.50", "87.0.4280.88"…
#> $ `Operating System`      <chr> "Windows NT 10.0", "Macintosh", "Windows NT 10…
#> $ Resolution              <chr> "1366x768", "1680x1050", "1366x768", "1536x864…
#> $ exclusions              <chr> "preview", "preview", "", "", "duration_quick"…
```

### Checking

The `check_*()` functions output messages about the number of rows that
meet the exclusion criteria. Because checks return only the rows meeting
the criteria, they **should not be connected via pipes** unless you want
to subset the second check criterion within the rows that meet the first
criterion. Thus, in general, `check_*()` functions should be used
individually. If you want to view the potential exclusions for multiple
criteria, use the `mark_*()` functions.

``` r
# Check for preview rows
qualtrics_text %>%
  check_preview()
#> ℹ 2 rows were collected as previews. It is highly recommended to exclude these rows before further processing.
#>             StartDate             EndDate         Status IPAddress Progress
#> 1 2020-12-11 12:06:52 2020-12-11 12:10:30 Survey Preview      <NA>      100
#> 2 2020-12-11 12:06:43 2020-12-11 12:11:27 Survey Preview      <NA>      100
#>   Duration (in seconds) Finished        RecordedDate        ResponseId
#> 1                   465     TRUE 2020-12-11 12:10:30 R_xLWiuPaNuURSFXY
#> 2                   545     TRUE 2020-12-11 12:11:27 R_Q5lqYw6emJQZx2o
#>   LocationLatitude LocationLongitude UserLanguage Browser      Version
#> 1         29.73694         -94.97599           EN  Chrome 88.0.4324.41
#> 2         39.74107        -121.82490           EN  Chrome 88.0.4324.50
#>   Operating System Resolution
#> 1  Windows NT 10.0   1366x768
#> 2        Macintosh  1680x1050
```

### Excluding

The `exclude_*()` functions remove the rows that meet exclusion
criteria. These, too, can be piped together. Since the output of each
function is a subset of the original data with the excluded rows
removed, the order of the functions will influence the reported number
of rows meeting the exclusion criteria.

``` r
# Exclude preview then incomplete progress rows
df <- qualtrics_text %>%
  exclude_duration(min_duration = 100) %>%
  exclude_progress()
#> ℹ 4 out of 100 rows of short and/or long duration were excluded, leaving 96 rows.
#> ℹ 4 out of 96 rows with incomplete progress were excluded, leaving 92 rows.
dim(df)
#> [1] 92 16
```

``` r
# Exclude incomplete progress then preview rows
df <- qualtrics_text %>%
  exclude_progress() %>%
  exclude_duration(min_duration = 100)
#> ℹ 6 out of 100 rows with incomplete progress were excluded, leaving 94 rows.
#> ℹ 2 out of 94 rows of short and/or long duration were excluded, leaving 92 rows.
dim(df)
#> [1] 92 16
```

Though the order of functions should not influence the final data set,
it may speed up processing large files by removing preview and
incomplete progress data first and waiting to check IP addresses and
locations after other exclusions have been performed.

``` r
# Exclude rows
df <- qualtrics_text %>%
  exclude_preview() %>%
  exclude_progress() %>%
  exclude_duplicates() %>%
  exclude_duration(min_duration = 100) %>%
  exclude_resolution() %>%
  exclude_ip() %>%
  exclude_location()
#> ℹ 2 out of 100 preview rows were excluded, leaving 98 rows.
#> ℹ 6 out of 98 rows with incomplete progress were excluded, leaving 92 rows.
#> ℹ 9 out of 92 duplicate rows were excluded, leaving 83 rows.
#> ℹ 2 out of 83 rows of short and/or long duration were excluded, leaving 81 rows.
#> ℹ 4 out of 81 rows with unacceptable screen resolution were excluded, leaving 77 rows.
#> ℹ 2 out of 77 rows with IP addresses outside of US were excluded, leaving 75 rows.
#> ℹ 4 out of 75 rows outside of the US were excluded, leaving 71 rows.
```

## Citing this package

To cite `{excluder}`, use:

> Stevens, J. R. (2021). excluder: An R package that checks for
> exclusion criteria in online data. *Journal of Open Source Software*,
> 6(67), 3893. <https://doi.org/10.21105/joss.03893>

## Contributing to this package

[Contributions](https://docs.ropensci.org/excluder/CONTRIBUTING.html) to
`{excluder}` are most welcome! Feel free to check out [open
issues](https://github.com/ropensci/excluder/issues) for ideas. And
[pull requests](https://github.com/ropensci/excluder/pulls) are
encouraged, but you may want to [raise an
issue](https://github.com/ropensci/excluder/issues/new/choose) or
[contact the maintainer](mailto:jeffrey.r.stevens@gmail.com) first.

Please note that the excluder project is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing
to this project, you agree to abide by its terms.

## Acknowledgments

I thank [Francine Goh](https://orcid.org/0000-0002-7364-4398) and Billy
Lim for comments on an early version of the package, as well as the
insightful feedback from [rOpenSci](https://ropensci.org/) editor [Mauro
Lepore](https://orcid.org/0000-0002-1986-7988) and reviewers [Joseph
O’Brien](https://orcid.org/0000-0001-9851-5077) and [Julia
Silge](https://orcid.org/0000-0002-3671-836X). This work was funded by
US National Science Foundation grant NSF-1658837.
# excluder 0.3.3

## BUG FIXES

* Updating `{dplyr}` to 1.0.8 caught problem with using `across()` inside `is.na()`. Instead, it now uses `if_all()` with `is.na()` as argument. Thanks to [@romainfrancois](https://github.com/romainfrancois) for pull request [\#7](https://github.com/ropensci/excluder/pull/7).
* `remove_label_rows()` now properly converts numeric data in Status and Finished columns.

# excluder 0.3.2

## MINOR IMPROVEMENTS

* The `remove_label_rows()` function can now rename columns to match the default column names used in all of the verb function arguments.
* The `mark_ip_()` function now checks for (1) internet connectivity, (2) whether the IP address data can be downloaded from https://www.iwik.org/ipcountry/, and (3) if the country code is valid. The function fails gracefully if any of these are not met.

### DOCUMENTATION UPDATES

* The README and vignette have clarified that multiple `check_*()` functions should not be piped.
* The `*_ip()` function documentation has clarified the internet connectivity requirements.

### PACKAGE DEVELOPMENT

* The paper for the `{excluder}` package has now been published in [Journal of Open Source Software](https://doi.org/10.21105/joss.03893). The JOSS badge has been added to the README.
* Tests for `*_ip()` functions are now skipped on CRAN to avoid timeout delays.
* Clean up package in preparation for submission to CRAN.


# excluder 0.3.1

### MINOR IMPROVEMENTS

* The messages now use `{cli}` to generate messages that more clearly outputs numbers and text.
* All functions now print the output data frames by default but can all be turned off with `print = FALSE`.
* `*_ip()` functions now include `include_na` as an argument (like `*_duplicates()` and `*_location()` functions), so users can decide whether to include NA values in the data that meet exclusion criteria.
* There are now five new utility functions that help simplify the primary verb functions. `keep_exclusion_column()` allows users to keep the exclusion column in the output from `check_*()` functions and moves the column the first column in the output. `mark_rows()` does the bulk of the work creating new columns for exclusion criteria and marking rows that meet the criteria. `print_data()` controls whether the output is printed to the console. `print_exclusion()` generates the message about how many rows were excluded by the `exclude_*()` functions. `validate_columns()` validates the number, names, and type of columns that are inputted as arguments in the verb functions.
* The NEWS.md file is now based on the rOpenSci template.


### BUG FIXES

* The `unite_exclusions()` function now properly removes multiple separators when multiple exclusion criteria are used.
* The `mark_duplicates()` function now properly counts and includes the correct number of NAs for both IP addresses and locations and properly prints data.


### PACKAGE DEVELOPMENT

* The `{excluder}` package has now been approved by and transferred to [rOpenSci](https://ropensci.org/). The package was peer reviewed by Joseph O'Brien ([@jmobrien](https://github.com/jmobrien)) and Julia Silge ([@juliasilge](https://github.com/juliasilge)), who are now listed as reviewers in the DESCRIPTION file.


# excluder 0.3.0

### NEW FEATURES

* The `mark_durations()` function now marks fast and slow durations separately.
* The primary functionality of the package has moved from the `check_*()` functions to the `mark_*()` functions. Thus, `check_*()` and `exclude_*()` now first call `mark_*()` to mark the rows with exclusion criteria, then filter the excluded rows. The documentation for `check_*()` and `exclude_*()` now inherit the arguments from `mark_*()`. This change has been updated in the README and Getting Started vignette.

### MINOR IMPROVEMENTS

* `exclude__*()` functions now have `print = FALSE` and `quiet = TRUE` set as default argument values.
* Calls to `rbind()` have been replaced with `bind_cols()` and `dplyr::pull()` has been replaced with `[[]]`.
* Calls to `all_of()` and `any_of()` now refer to `{tidyselect}` rather than `{dplyr}`.
* `if()` statements are now more robust by using `identical()` rather than `==` and `&&` instead of `&`.
* The `{stringr}` package is now imported instead of suggested.
* All mark, check, and exclude functions for a particular exclusion type have been combined into a single R file. So now each exclusion type has its own R file. Similarly, data file scripts have been combined into a single file.

### BUG FIXES

* The `*_ip()` functions and documentation have been updated to fix a bug/typo to clarify that they mark, check, and exclude rows with IP addresses outside of the specified country.

### DEPRECATED AND DEFUNCT

* `collapse_exclusions()` has been renamed `unite_exclusions()` to match `{tidyverse}` terminology. `collapse_exclusions()` is now deprecated and will be removed in a future version, use `unite_exlusions()`. `unite_exlusions()` also switched from using NA to "" for rows with no exclusions. Combined columns now no longer have leftover separators.

### DOCUMENTATION FIXES

* Package links are replaced with external URLs.


# excluder 0.2.2

* Lifecycle has been updated to Stable, and the repo status of Active has been added.
* There is now a Getting Started Excluding Data vignette that gives an overview of package usage.
* Users can now specify the separator used between exclusion types in the `collapse_exclusions()` function. They can also opt to not remove the excluded columns.
* A bug in the `collapse_exclusions()` function was fixed, so now a subset of exclusion columns can be collapsed into a single column.
* A bug in the `remove_label_rows()` function was fixed, so now the Finished column is converted to a logical vector.
* A codemeta.json file was created.
* URLs have been replaced or tweaked to address CRAN package check warnings.
* Functions are now organized by topic in the Reference page of the website.

# excluder 0.2.1

* The argument name for the data frame switched from `.data` to `x` to avoid confusion with the `{rlang}` use of `.data`.
* Instead of using `quo()` and `sym()` to create new names for columns used as arguments, `.data[[var]]` is now used.
* The `dupe_count` column was removed from `check_duplicates()` output. Tests were adjusted to account for the new number of columns.
* `check_duplicates()` now specifies the number of NA columns.
* URLs have been replaced or tweaked to address CRAN package check warnings.
* Links function reference pages have been added to the README.

# excluder 0.2.0

* The `deidentify()` function was added, which removes IP address, location, and computer information columns.

# excluder 0.1.0

* The `check_qualtrics` argument was removed from `remove_label_rows()` because the functionality did not make sense. This breaks backwards compatibility.
* `remove_label_rows()` now only filters out label rows if label rows are present and outputs invisibly.
* Tests were added for the `qualtrics_raw` dataset and the `remove_label_rows()` function.
* Package-level documentation was created.

# excluder 0.0.1

* `remove_label_rows()` now converts character columns to dates for multiple date formats, including YYYY-MM-DD HH:MM:SS, YYYY-MM-DD HH:MM, MM:DD:YYYY HH:MM:SS, and MM:DD:YYYY HH:MM (#1).
* Code of Conduction and Contributor Guide are added.
* Citation and Contributor sections are added to README.

# excluder 0.0.0.1

* Initial GitHub release
## R CMD check results

0 errors | 0 warnings | 0 note

Though the first version of this package was just submitted/accepted recently, a recent update to {dplyr} discovered a bug in my usage of {dplyr} functions, and @romainfrancois requested that I release an update as soon as possible (https://github.com/ropensci/excluder/pull/7).
# Contributing to excluder

This outlines how to propose a change to excluder. 

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("jeffreyrstevens/excluder", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the excluder project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Reproducible example**
Please include a [reproducible example](https://reprex.tidyverse.org/) when possible.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
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

# excluder

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![lifecycle](man/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/excluder)](https://cran.r-project.org/package=excluder)
<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/excluder)](https://CRAN.R-project.org/package=excluder) -->
<!-- [![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/excluder?color=orange)](https://CRAN.R-project.org/package=excluder) -->


[![R-CMD-check](https://github.com/ropensci/excluder/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/excluder/actions)
[![Codecov test coverage](https://codecov.io/gh/jeffreyrstevens/excluder/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jeffreyrstevens/excluder?branch=main)

[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/455_status.svg)](https://github.com/ropensci/software-review/issues/455)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03893/status.svg)](https://doi.org/10.21105/joss.03893)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5648202.svg)](https://doi.org/10.5281/zenodo.5648202)
<!-- badges: end -->

The goal of [`{excluder}`](https://docs.ropensci.org/excluder/) is to facilitate checking for, marking, and excluding rows of data frames for common exclusion criteria. This package applies to data collected from [Qualtrics](https://www.qualtrics.com/) surveys, and default column names come from importing data with the  [`{qualtRics}`](https://docs.ropensci.org/qualtRics/) package.

This may be most useful for [Mechanical Turk](https://www.mturk.com/) data to screen for duplicate entries from the same location/IP address or entries from locations outside of the United States. But it can be used more generally to exclude based on response durations, preview status, progress, or screen resolution.

More details are available on the package [website](https://docs.ropensci.org/excluder/) and the [getting started vignette](https://docs.ropensci.org/excluder/articles/getting_started.html).

## Installation

You can install the stable released version of `{excluder}` from [CRAN](https://cran.r-project.org/package=excluder) with:

```{r eval = FALSE}
install.packages("excluder")
```

You can install developmental versions from [GitHub](https://github.com/) with:

```{r eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/excluder")
```


## Verbs
This package provides three primary verbs:

* `mark` functions add a new column to the original data frame that labels the rows meeting the exclusion criteria. This is useful to label the potential exclusions for future processing without changing the original data frame.
* `check` functions search for the exclusion criteria and output a message with the number of rows meeting the criteria and a data frame of the rows meeting the criteria. This is useful for viewing the potential exclusions.
* `exclude` functions remove rows meeting the exclusion criteria. This is safest to do after checking the rows to ensure the exclusions are correct.

## Exclusion types
This package provides seven types of exclusions based on Qualtrics metadata. If you have ideas for other metadata exclusions, please submit them as [issues](https://github.com/ropensci/excluder/issues). Note, the intent of this package is not to develop functions for excluding rows based on survey-specific data but on general, frequently used metadata.

* `duplicates` works with rows that have duplicate IP addresses and/or locations (latitude/longitude).
* `duration` works with rows whose survey completion time is too short and/or too long.
* `ip` works with rows whose IP addresses are not found in the specified country (note: this exclusion type requires an internet connection to download the country's IP ranges).
* `location` works with rows whose latitude and longitude are not found in the United States.
* `preview` works with rows that are survey previews.
* `progress` works with rows in which the survey was not complete.
* `resolution` works with rows whose screen resolution is not acceptable.


## Usage

The verbs and exclusion types combine with `_` to create the functions, such as [`check_duplicates()`](https://docs.ropensci.org/excluder/reference/check_duplicates.html), [`exclude_ip()`](https://docs.ropensci.org/excluder/reference/exclude_ip.html), and [`mark_duration()`](https://docs.ropensci.org/excluder/reference/mark_duration.html). Multiple functions can be linked together using the [`{magrittr}`](https://magrittr.tidyverse.org/) pipe `%>%`. For datasets downloaded directly from Qualtrics, use [`remove_label_rows()`](https://docs.ropensci.org/excluder/reference/remove_label_rows.html) to remove the first two rows of labels and convert date and numeric columns in the metadata, and use [`deidentify()`](https://docs.ropensci.org/excluder/reference/deidentify.html) to remove standard Qualtrics columns with identifiable information (e.g., IP addresses, geolocation).

### Marking
The `mark_*()` functions output the original data set with a new column specifying rows that meet the exclusion criteria. These can be piped together with `%>%` for multiple exclusion types.

```{r mark1}
library(excluder)
# Mark preview and short duration rows
df <- qualtrics_text %>%
  mark_preview() %>%
  mark_duration(min_duration = 200)
tibble::glimpse(df)
```

Use the [`unite_exclusions()`](https://docs.ropensci.org/excluder/reference/unite_exclusions.html) function to unite all of the marked columns into a single column.
```{r mark2}
# Collapse labels for preview and short duration rows
df <- qualtrics_text %>%
  mark_preview() %>%
  mark_duration(min_duration = 200) %>%
  unite_exclusions()
tibble::glimpse(df)
```

### Checking

The `check_*()` functions output messages about the number of rows that meet the exclusion criteria. Because checks return only the rows meeting the criteria, they **should not be connected via pipes** unless you want to subset the second check criterion within the rows that meet the first criterion. Thus, in general, `check_*()` functions should be used individually. If you want to view the potential exclusions for multiple criteria, use the `mark_*()` functions.

```{r check1}
# Check for preview rows
qualtrics_text %>%
  check_preview()
```

### Excluding
The `exclude_*()` functions remove the rows that meet exclusion criteria. These, too, can be piped together. Since the output of each function is a subset of the original data with the excluded rows removed, the order of the functions will influence the reported number of rows meeting the exclusion criteria.

```{r exclude1}
# Exclude preview then incomplete progress rows
df <- qualtrics_text %>%
  exclude_duration(min_duration = 100) %>%
  exclude_progress()
dim(df)
```
```{r exclude2}
# Exclude incomplete progress then preview rows
df <- qualtrics_text %>%
  exclude_progress() %>%
  exclude_duration(min_duration = 100)
dim(df)
```
Though the order of functions should not influence the final data set, it may speed up processing large files by removing preview and incomplete progress data first and waiting to check IP addresses and locations after other exclusions have been performed.

```{r exclude3}
# Exclude rows
df <- qualtrics_text %>%
  exclude_preview() %>%
  exclude_progress() %>%
  exclude_duplicates() %>%
  exclude_duration(min_duration = 100) %>%
  exclude_resolution() %>%
  exclude_ip() %>%
  exclude_location()
```

## Citing this package

To cite `{excluder}`, use:

> Stevens, J. R. (2021). excluder: An R package that checks for exclusion criteria in online data. _Journal of Open Source Software_, 6(67), 3893. https://doi.org/10.21105/joss.03893


## Contributing to this package

[Contributions](https://docs.ropensci.org/excluder/CONTRIBUTING.html) to `{excluder}` are most welcome! Feel free to check out [open issues](https://github.com/ropensci/excluder/issues) for ideas. And [pull requests](https://github.com/ropensci/excluder/pulls) are encouraged, but you may want to [raise an issue](https://github.com/ropensci/excluder/issues/new/choose) or [contact the maintainer](mailto:jeffrey.r.stevens@gmail.com) first.

Please note that the excluder project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

## Acknowledgments

I thank [Francine Goh](https://orcid.org/0000-0002-7364-4398) and Billy Lim for comments on an early version of the package, as well as the insightful feedback from [rOpenSci](https://ropensci.org/) editor [Mauro Lepore](https://orcid.org/0000-0002-1986-7988) and reviewers [Joseph O'Brien](https://orcid.org/0000-0001-9851-5077) and [Julia Silge](https://orcid.org/0000-0002-3671-836X). This work was funded by US National Science Foundation grant NSF-1658837.
---
title: "Get started excluding data"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(excluder)
```

The `{excluder}` package facilitates marking, checking for, and excluding rows of data frames^[Though the functions take data frames input, they produce [tibbles](https://tibble.tidyverse.org/).] for common online participant exclusion criteria. This package applies to online data with a focus on data collected from [Qualtrics](https://www.qualtrics.com) surveys, and default column names come from importing data with the  [`{qualtRics}`](https://docs.ropensci.org/qualtRics/) package. This may be most useful for [Mechanical Turk](https://www.mturk.com/) data to screen for duplicate entries from the same location/IP address or entries from locations outside of the United States.  However, the package can be used for data from other sources that need to be screened for IP address, location, duplicated locations, participant progress, survey completion time, and/or screen resolution.


## Usage

The package has three core verbs:

1. `mark_*()` functions add a new column to the original data frame that labels the rows meeting the exclusion criteria. This is useful to label the potential exclusions for future processing without changing the original data frame.
1. `check_*()` functions search for the exclusion criteria and output a message with the number of rows meeting the criteria and a data frame of the rows meeting the criteria. This is useful for viewing the potential exclusions.
1. `exclude_*()` functions remove rows meeting the exclusion criteria. This is safest to do after checking the rows to ensure the exclusions are correct.

The `check_*()` and `exclude_*()` functions call the `mark_*()` function internally and then filter either the excluded or non-excluded rows. So avoid combining different verbs to sidestep unnecessary passes through the data.

The package provides seven types of exclusions based on Qualtrics metadata:

1. `duplicates` works with rows that have duplicate IP addresses and/or locations (latitude/longitude), using [`janitor::get_dupes()`](https://sfirke.github.io/janitor/reference/get_dupes.html).
1. `duration` works with rows whose survey completion time is too short and/or too long.
1. `ip` works with rows whose IP addresses are not found in the specified country (note: this exclusion type requires an internet connection to download the country's IP ranges), using package [`{iptools}`](https://github.com/hrbrmstr/iptools).
1. `location` works with rows whose latitude and longitude are not found in the United States.
1. `preview` works with rows that are survey previews.
1. `progress` works with rows in which the survey was not complete.
1. `resolution` works with rows whose screen resolution is not acceptable.

The verbs combine with the exclusion types to generate functions. For instance, `mark_duplicates()` will mark duplicate rows and `exclude_preview()` will exclude preview rows.

There are also helper functions:

1. [`unite_exclusions()`](https://docs.ropensci.org//excluder/reference/unite_exclusions.html)  unites all of the columns marked by `mark` functions into a single column (each use of a `mark` function creates a new column).
1. [`deidentify()`](https://docs.ropensci.org//excluder/reference/deidentify.html) removes standard Qualtrics columns with identifiable information.
1. [`remove_label_rows()`](https://docs.ropensci.org//excluder/reference/remove_label_rows.html) removes the first two rows of labels and convert date and numeric columns in the metadata.

## Preparing your data
If you use the [`fetch_survey()`](https://docs.ropensci.org/qualtRics/reference/fetch_survey.html) from the `{qualtRics}` package to import your Qualtrics data, it will automatically remove the first two label rows from the data set. However, if you directly download your data from Qualtrics, it will include two rows in your data set that include label information. This has two side effects: (1) there are non-data rows that need to be removed from your data set, and (2) all of your columns will be imported as character data types.

The [`remove_label_rows()`](https://docs.ropensci.org//excluder/reference/remove_label_rows.html) function will remove these two label rows. Also, by default, it will coerce the Qualtrics metadata columns from character types to the correct formats (e.g., _`StartDate`_ is coerced to a date, _`Progress`_ is coerced to numeric). So if you download your data from Qualtrics, you will need to run this function on your data before proceeding.

```{r}
dplyr::glimpse(qualtrics_raw)
#
# Remove label rows and coerce metadata columns
df <- remove_label_rows(qualtrics_raw) %>% 
  dplyr::glimpse()
```


## Marking observations
The core verbs in this package mark them for future processing. The `mark_*()` suite of functions creates a new column for each mark function used that marks which observations violate the exclusion criterion. They print a message about the number of observations meeting each exclusion criteria. Mark functions return a data frame identical to the original with additional columns marking exclusions.

```{r mark1}
# Mark observations run as preview
df %>% 
  mark_preview() %>% 
  dplyr::glimpse()
```

Notice the new _`exclusion_preview`_ column at the end of the data frame. It has marked the first two observations as `preview`.

Piping multiple mark functions will create multiple rows marking observations for exclusion.

```{r mark2}
# Mark preview and incomplete observations
df %>% 
  mark_preview() %>% 
  mark_progress() %>% 
  dplyr::glimpse()
```

To unite all of the marked columns into a single column, use the `unite_exclusions()` function. This will create a new `exclusions` columns that will unite all exclusions in each observation into a single column. Here we move the combined _`exclusions`_ column to the beginning of the data frame to view it.

```{r mark3}
df %>% 
  mark_preview() %>% 
  mark_duration(min = 500) %>% 
  unite_exclusions() %>% 
  dplyr::relocate(exclusions, .before = StartDate)
```

Multiple exclusions are separated by `,` by default, but the separating character can be controlled by the `separator` argument. By default, the multiple exclusion columns are removed from the final data frame, but this can be turned off by setting the `remove` argument to `FALSE`. 

```{r mark4}
df %>% 
  mark_preview() %>% 
  mark_duration(min = 500) %>% 
  unite_exclusions(separator = ";", remove = FALSE) %>% 
  dplyr::relocate(exclusions, .before = StartDate)
```


## Checking observations

The `check_*()` suite of functions return a data frame that includes only the observations that meet the criterion. Since these functions first run the appropriate `mark_*()` function, they also print a message about the number of observations that meet the exclusion criterion.
```{r check1}
# Check for rows with incomplete data
df %>%
  check_progress()
```
```{r check2}
# Check for rows with durations less than 100 seconds
df %>% 
  check_duration(min_duration = 100)
```

Because checks return only the rows meeting the criteria, they **should not be connected via pipes** unless you want to subset the second check criterion within the rows that meet the first criterion.

```{r check3}
# Check for rows with durations less than 100 seconds in rows that did not complete the survey
df %>%
  check_progress() %>%
  check_duration(min_duration = 100)
```

To check all data for multiple criteria, use the `mark_*()` functions followed by a filter.

```{r mark_check}
# Check for multiple criteria
df %>% 
  mark_preview() %>% 
  mark_duration(min = 500) %>% 
  unite_exclusions() %>% 
  dplyr::filter(exclusions != "")
```

## Excluding observations
The `exclude_*()` suite of function will return a data frame that has removed observations that match the exclusion criteria. Exclude functions print their own messages about the number of observations excluded.

```{r exclude1}
# Exclude survey responses used to preview the survey
df %>% 
  exclude_preview() %>% 
  dplyr::glimpse()
```

Piping will apply subsequent excludes to the data frames with the previous excludes already applied. Therefore, it often makes sense to remove the preview surveys and incomplete surveys before checking other exclusion types to speed processing.

```{r exclude2}
# Exclude preview then incomplete progress rows then duplicate locations and IP addresses
df %>%
  exclude_preview() %>%
  exclude_progress() %>%
  exclude_duplicates(print = FALSE)
```

## Messages and console output
Messages about the number of rows meeting the exclusion criteria are printed to the console by default. These messages are generated by the `mark_*()` functions and carry over to `check_*()` functions. They can be turn off by setting `quiet` to `TRUE`.

```{r quiet}
# Turn off marking/checking messages with quiet = TRUE
df %>%
  check_progress(quiet = TRUE)
```

Note that `exclude_*()` functions have the `mark_*()` messages turned off by default and produce their own messages about exclusions. To silence these messages, set `silent` to `TRUE`.

```{r silent}
# Turn off exclusion messages with silent = TRUE
df %>%
  exclude_preview(silent = TRUE) %>%
  exclude_progress(silent = TRUE) %>%
  exclude_duplicates(silent = TRUE)
```

Though `exclude_*()` functions do not print the data frame to the console, `mark_*()` and `check_*()` do. To avoid printing to the console, set `print` = `FALSE`.

```{r printoff}
# Turn off marking/checking printing data frame with print = FALSE
df %>%
  check_progress(print = FALSE)
```

## Deidentifying data

By default, Qualtrics records participant IP address and location.^[To avoid recording this potentially identifiable information, go to _Survey options_ > _Security_ > _Anonymize responses_.] You can also record properties of the participants' computers such as operating system, web browser type and version, and screen resolution.^[You can turn these on by adding a new question (usually to the beginning of your survey) and choosing _Meta info_ as the question type.] While these pieces of information can be useful for excluding observations, they carry potentially identifiable information. Therefore, you may want to remove them from data frame before saving or processing it. The `deidentify()` function removes potentially identifiable data columns collected by Qualtrics. By default, the function uses a strict rule to remove IP address, location, and computer information (browser type and version, operating system, and screen resolution).

```{r deidentify1}
# Exclude preview then incomplete progress rows
df %>%
  exclude_preview() %>%
  exclude_progress() %>%
  exclude_duplicates() %>%
  deidentify() %>%
  dplyr::glimpse()
```

If the computer information is not considered sensitive, it can be kept by setting the `strict` argument to `FALSE`, thereby only removing IP address and location.

```{r deidentify2}
# Exclude preview then incomplete progress rows
df %>%
  exclude_preview() %>%
  exclude_progress() %>%
  exclude_duplicates() %>%
  deidentify(strict = FALSE) %>%
  dplyr::glimpse()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/location.R
\name{mark_location}
\alias{mark_location}
\title{Mark locations outside of US}
\usage{
mark_location(
  x,
  id_col = "ResponseId",
  location_col = c("LocationLatitude", "LocationLongitude"),
  include_na = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{location_col}{Two element vector specifying columns for latitude
and longitude (in that order).}

\item{include_na}{Logical indicating whether to include rows with NA in
latitude and longitude columns in the output list of potentially excluded
data.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes a column marking rows
that are located outside of the US and (if \code{include_na == FALSE}) rows with
no location information.
For a function that checks for these rows, use \code{\link[=check_location]{check_location()}}.
For a function that excludes these rows, use \code{\link[=exclude_location]{exclude_location()}}.
}
\description{
The \code{mark_location()} function creates a column labeling
rows that have locations outside of the US.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The function only works for the United States.
It uses the #' \code{\link[maps:map.where]{maps::map.where()}} to determine if latitude and longitude
are inside the US.

The function outputs to console a message about the number of rows
with locations outside of the US.
}
\examples{
# Mark locations outside of the US
data(qualtrics_text)
df <- mark_location(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_location()
}
\seealso{
Other location functions: 
\code{\link{check_location}()},
\code{\link{exclude_location}()}

Other mark functions: 
\code{\link{mark_duplicates}()},
\code{\link{mark_duration}()},
\code{\link{mark_ip}()},
\code{\link{mark_preview}()},
\code{\link{mark_progress}()},
\code{\link{mark_resolution}()}
}
\concept{location functions}
\concept{mark functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duplicates.R
\name{exclude_duplicates}
\alias{exclude_duplicates}
\title{Exclude rows with duplicate IP addresses and/or locations}
\usage{
exclude_duplicates(
  x,
  id_col = "ResponseId",
  ip_col = "IPAddress",
  location_col = c("LocationLatitude", "LocationLongitude"),
  dupl_ip = TRUE,
  dupl_location = TRUE,
  include_na = FALSE,
  quiet = TRUE,
  print = TRUE,
  silent = FALSE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{ip_col}{Column name for IP addresses.}

\item{location_col}{Two element vector specifying columns for latitude and
longitude (in that order).}

\item{dupl_ip}{Logical indicating whether to check IP addresses.}

\item{dupl_location}{Logical indicating whether to check latitude and
longitude.}

\item{include_na}{Logical indicating whether to include rows with NAs for
IP address and location as potentially excluded rows.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}

\item{silent}{Logical indicating whether to print message to console. Note
this argument controls the exclude message not the check message.}
}
\value{
An object of the same type as \code{x} that excludes rows
with duplicate IP addresses and/or locations.
For a function that just checks for and returns duplicate rows,
use \code{\link[=check_duplicates]{check_duplicates()}}. For a function that marks these rows,
use \code{\link[=mark_duplicates]{mark_duplicates()}}.
}
\description{
The \code{exclude_duplicates()} function removes
rows of data that have the same IP address and/or same latitude and
longitude. The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
By default, IP address and location are both checked, but they can be
checked separately with the \code{dupl_ip} and \code{dupl_location} arguments.

The function outputs to console separate messages about the number of
rows with duplicate IP addresses and rows with duplicate locations.
These counts are computed independently, so rows may be counted for both
types of duplicates.
}
\examples{
# Exclude duplicate IP addresses and locations
data(qualtrics_text)
df <- exclude_duplicates(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_duplicates()

# Exclude only for duplicate locations
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_duplicates(dupl_location = FALSE)
}
\seealso{
Other duplicates functions: 
\code{\link{check_duplicates}()},
\code{\link{mark_duplicates}()}

Other exclude functions: 
\code{\link{exclude_duration}()},
\code{\link{exclude_ip}()},
\code{\link{exclude_location}()},
\code{\link{exclude_preview}()},
\code{\link{exclude_progress}()},
\code{\link{exclude_resolution}()}
}
\concept{duplicates functions}
\concept{exclude functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{keep_marked_column}
\alias{keep_marked_column}
\title{Keep column with marked rows}
\usage{
keep_marked_column(x, column, keep)
}
\arguments{
\item{x}{Data set.}

\item{column}{Name of exclusion column.}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}
}
\description{
For check_*() functions, keep the column that has marked rows and move to
first column or remove the column depending on \code{keep} flag.
\emph{This function is not exported.}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ip.R
\name{exclude_ip}
\alias{exclude_ip}
\title{Exclude IP addresses from outside of a specified country.}
\usage{
exclude_ip(
  x,
  id_col = "ResponseId",
  ip_col = "IPAddress",
  country = "US",
  include_na = FALSE,
  quiet = TRUE,
  print = TRUE,
  silent = FALSE
)
}
\arguments{
\item{x}{Data frame or tibble (preferably imported from Qualtrics using
\{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{ip_col}{Column name for IP addresses.}

\item{country}{Two-letter abbreviation of country to check (default is "US").}

\item{include_na}{Logical indicating whether to include rows with NA in
IP address column in the output list of potentially excluded data.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}

\item{silent}{Logical indicating whether to print message to console. Note
this argument controls the exclude message not the check message.}
}
\value{
An object of the same type as \code{x} that excludes rows
with IP addresses outside of the specified country.
For a function that checks these rows, use \code{\link[=check_ip]{check_ip()}}.
For a function that marks these rows, use \code{\link[=mark_ip]{mark_ip()}}.
}
\description{
The \code{exclude_ip()} function removes rows of data that have
IP addresses from outside the specified country.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The function uses \code{\link[iptools:country_ranges]{iptools::country_ranges()}} to assign IP addresses to
specific countries using
\href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}{ISO 3166-1 alpha-2 country codes}.

The function outputs to console a message about the number of rows
with IP addresses outside of the specified country. If there are \code{NA}s for IP
addresses (likely due to including preview data---see \code{\link[=check_preview]{check_preview()}}), it
will print a message alerting to the number of rows with \code{NA}s.
}
\note{
This function \strong{requires internet connectivity} as it uses the
\code{\link[iptools:country_ranges]{iptools::country_ranges()}} function, which pulls daily updated data
from \url{http://www.iwik.org/ipcountry/}. It only updates the data once
per session, as it caches the results for future work during the session.
}
\examples{
# Exclude IP addresses outside of the US
data(qualtrics_text)
df <- exclude_ip(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_ip()

# Exclude IP addresses outside of Germany
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_ip(country = "DE")
}
\seealso{
Other ip functions: 
\code{\link{check_ip}()},
\code{\link{mark_ip}()}

Other exclude functions: 
\code{\link{exclude_duplicates}()},
\code{\link{exclude_duration}()},
\code{\link{exclude_location}()},
\code{\link{exclude_preview}()},
\code{\link{exclude_progress}()},
\code{\link{exclude_resolution}()}
}
\concept{exclude functions}
\concept{ip functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resolution.R
\name{mark_resolution}
\alias{mark_resolution}
\title{Mark unacceptable screen resolution}
\usage{
mark_resolution(
  x,
  width_min = 1000,
  height_min = 0,
  id_col = "ResponseId",
  res_col = "Resolution",
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{width_min}{Minimum acceptable screen width.}

\item{height_min}{Minimum acceptable screen height.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{res_col}{Column name for screen resolution (in format widthxheight).}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes a column marking rows
that have unacceptable screen resolutions.
For a function that checks for these rows, use \code{\link[=check_resolution]{check_resolution()}}.
For a function that excludes these rows, use \code{\link[=exclude_resolution]{exclude_resolution()}}.
}
\description{
The \code{mark_resolution()} function creates a column labeling
rows that have unacceptable screen resolution.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.

The function outputs to console a message about the number of rows
with unacceptable screen resolution.
}
\examples{
# Mark low screen resolutions
data(qualtrics_text)
df <- mark_resolution(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_resolution()
}
\seealso{
Other resolution functions: 
\code{\link{check_resolution}()},
\code{\link{exclude_resolution}()}

Other mark functions: 
\code{\link{mark_duplicates}()},
\code{\link{mark_duration}()},
\code{\link{mark_ip}()},
\code{\link{mark_location}()},
\code{\link{mark_preview}()},
\code{\link{mark_progress}()}
}
\concept{mark functions}
\concept{resolution functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preview.R
\name{exclude_preview}
\alias{exclude_preview}
\title{Exclude survey previews}
\usage{
exclude_preview(
  x,
  id_col = "ResponseId",
  preview_col = "Status",
  quiet = TRUE,
  print = TRUE,
  silent = FALSE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{preview_col}{Column name for survey preview.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}

\item{silent}{Logical indicating whether to print message to console. Note
this argument controls the exclude message not the check message.}
}
\value{
An object of the same type as \code{x} that excludes rows
that are survey previews.
For a function that checks for these rows, use \code{\link[=check_preview]{check_preview()}}.
For a function that marks these rows, use \code{\link[=mark_preview]{mark_preview()}}.
}
\description{
The \code{exclude_preview()} function removes
rows that are survey previews.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The preview column in Qualtrics can be a numeric or character vector
depending on whether it is exported as choice text or numeric values.
This function works for both.

The function outputs to console a message about the number of rows
that are survey previews.
}
\examples{
# Exclude survey previews
data(qualtrics_text)
df <- exclude_preview(qualtrics_text)

# Works for Qualtrics data exported as numeric values, too
df <- qualtrics_numeric \%>\%
  exclude_preview()

# Do not print rows to console
df <- qualtrics_text \%>\%
  exclude_preview(print = FALSE)
}
\seealso{
Other preview functions: 
\code{\link{check_preview}()},
\code{\link{mark_preview}()}

Other exclude functions: 
\code{\link{exclude_duplicates}()},
\code{\link{exclude_duration}()},
\code{\link{exclude_ip}()},
\code{\link{exclude_location}()},
\code{\link{exclude_progress}()},
\code{\link{exclude_resolution}()}
}
\concept{exclude functions}
\concept{preview functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ip.R
\name{check_ip}
\alias{check_ip}
\title{Check for IP addresses from outside of a specified country.}
\usage{
check_ip(
  x,
  id_col = "ResponseId",
  ip_col = "IPAddress",
  country = "US",
  include_na = FALSE,
  keep = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame or tibble (preferably imported from Qualtrics using
\{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{ip_col}{Column name for IP addresses.}

\item{country}{Two-letter abbreviation of country to check (default is "US").}

\item{include_na}{Logical indicating whether to include rows with NA in
IP address column in the output list of potentially excluded data.}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes the rows with
IP addresses outside of the specified country.
For a function that marks these rows, use \code{\link[=mark_ip]{mark_ip()}}.
For a function that excludes these rows, use \code{\link[=exclude_ip]{exclude_ip()}}.
}
\description{
The \code{check_ip()} function subsets rows of data, retaining rows
that have IP addresses from outside the specified country.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The function uses \code{\link[iptools:country_ranges]{iptools::country_ranges()}} to assign IP addresses to
specific countries using
\href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}{ISO 3166-1 alpha-2 country codes}.

The function outputs to console a message about the number of rows
with IP addresses outside of the specified country. If there are \code{NA}s for IP
addresses (likely due to including preview data---see \code{\link[=check_preview]{check_preview()}}), it
will print a message alerting to the number of rows with \code{NA}s.
}
\note{
This function \strong{requires internet connectivity} as it uses the
\code{\link[iptools:country_ranges]{iptools::country_ranges()}} function, which pulls daily updated data
from \url{https://www.iwik.org/ipcountry/}. It only updates the data once
per session, as it caches the results for future work during the session.
}
\examples{
# Check for IP addresses outside of the US
data(qualtrics_text)
check_ip(qualtrics_text)

# Remove preview data first
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_ip()

# Check for IP addresses outside of Germany
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_ip(country = "DE")

# Do not print rows to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_ip(print = FALSE)

# Do not print message to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_ip(quiet = TRUE)
}
\seealso{
Other ip functions: 
\code{\link{exclude_ip}()},
\code{\link{mark_ip}()}

Other check functions: 
\code{\link{check_duplicates}()},
\code{\link{check_duration}()},
\code{\link{check_location}()},
\code{\link{check_preview}()},
\code{\link{check_progress}()},
\code{\link{check_resolution}()}
}
\concept{check functions}
\concept{ip functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/location.R
\name{exclude_location}
\alias{exclude_location}
\title{Exclude locations outside of US}
\usage{
exclude_location(
  x,
  id_col = "ResponseId",
  location_col = c("LocationLatitude", "LocationLongitude"),
  include_na = FALSE,
  quiet = TRUE,
  print = TRUE,
  silent = FALSE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{location_col}{Two element vector specifying columns for latitude
and longitude (in that order).}

\item{include_na}{Logical indicating whether to include rows with NA in
latitude and longitude columns in the output list of potentially excluded
data.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}

\item{silent}{Logical indicating whether to print message to console. Note
this argument controls the exclude message not the check message.}
}
\value{
An object of the same type as \code{x} that excludes rows
that are located outside of the US and (if \code{include_na == FALSE}) rows with
no location information.
For a function that checks for these rows, use \code{\link[=check_location]{check_location()}}.
For a function that marks these rows, use \code{\link[=mark_location]{mark_location()}}.
}
\description{
The \code{exclude_location()} function removes
rows that have locations outside of the US.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The function only works for the United States.
It uses the #' \code{\link[maps:map.where]{maps::map.where()}} to determine if latitude and longitude
are inside the US.

The function outputs to console a message about the number of rows
with locations outside of the US.
}
\examples{
# Exclude locations outside of the US
data(qualtrics_text)
df <- exclude_location(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_location()
}
\seealso{
Other location functions: 
\code{\link{check_location}()},
\code{\link{mark_location}()}

Other exclude functions: 
\code{\link{exclude_duplicates}()},
\code{\link{exclude_duration}()},
\code{\link{exclude_ip}()},
\code{\link{exclude_preview}()},
\code{\link{exclude_progress}()},
\code{\link{exclude_resolution}()}
}
\concept{exclude functions}
\concept{location functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duration.R
\name{check_duration}
\alias{check_duration}
\title{Check for minimum or maximum durations}
\usage{
check_duration(
  x,
  min_duration = 10,
  max_duration = NULL,
  id_col = "ResponseId",
  duration_col = "Duration (in seconds)",
  keep = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{min_duration}{Minimum duration that is too fast in seconds.}

\item{max_duration}{Maximum duration that is too slow in seconds.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{duration_col}{Column name for durations.}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes the rows with fast and/or
slow duration.
For a function that marks these rows, use \code{\link[=mark_duration]{mark_duration()}}.
For a function that excludes these rows, use \code{\link[=exclude_duration]{exclude_duration()}}.
}
\description{
The \code{check_duration()} function subsets rows of data, retaining rows
that have durations that are too fast or too slow.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
By default, minimum durations of 10 seconds are checked, but either
minima or maxima can be checked with the \code{min_duration} and
\code{max_duration} arguments. The function outputs to console separate
messages about the number of rows that are too fast or too slow.

This function returns the fast and slow rows.
}
\examples{
# Check for durations faster than 100 seconds
data(qualtrics_text)
check_duration(qualtrics_text, min_duration = 100)

# Remove preview data first
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_duration(min_duration = 100)

# Check only for durations slower than 800 seconds
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_duration(max_duration = 800)

# Do not print rows to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_duration(min_duration = 100, print = FALSE)

# Do not print message to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_duration(min_duration = 100, quiet = TRUE)
}
\seealso{
Other duration functions: 
\code{\link{exclude_duration}()},
\code{\link{mark_duration}()}

Other check functions: 
\code{\link{check_duplicates}()},
\code{\link{check_ip}()},
\code{\link{check_location}()},
\code{\link{check_preview}()},
\code{\link{check_progress}()},
\code{\link{check_resolution}()}
}
\concept{check functions}
\concept{duration functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validate_columns}
\alias{validate_columns}
\title{Check number, names, and type of columns}
\usage{
validate_columns(x, column)
}
\arguments{
\item{x}{Data set.}

\item{column}{Name of column argument to check.}
}
\description{
Determines whether the correct number and names of columns were specified
as arguments to the functions. \emph{This function is not exported.}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_label_rows.R
\name{remove_label_rows}
\alias{remove_label_rows}
\title{Remove two initial rows created in Qualtrics data}
\usage{
remove_label_rows(x, convert = TRUE, rename = FALSE)
}
\arguments{
\item{x}{Data frame (downloaded from Qualtrics).}

\item{convert}{Logical indicating whether to convert/coerce date, logical and
numeric columns from the metadata.}

\item{rename}{Logical indicating whether to rename columns based on first row
of data.}
}
\value{
An object of the same type as \code{x} that excludes Qualtrics label rows and
with date, logical, and numeric metadata columns converted to the correct
data class.
}
\description{
The \code{remove_label_rows()} function filters out the initial label rows from
datasets downloaded from \href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
The function (1) checks if the data set uses Qualtrics column names,
(2) checks if label rows are already used as column names,
(3) removes label rows if present, and (4) converts date, logical, and
numeric metadata columns to proper data type. Datasets imported using
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{qualtRics::fetch_survey()}
should not need this function.

The \code{convert} argument only converts the \emph{StartDate}, \emph{EndDate},
\emph{RecordedDate}, \emph{Progress}, \emph{Finished}, \emph{Duration (in seconds)},
\emph{LocationLatitude}, and \emph{LocationLongitude} columns. To convert other data
columns, see \code{\link[dplyr:mutate]{dplyr::mutate()}}.
}
\examples{
# Remove label rows
data(qualtrics_raw)
df <- remove_label_rows(qualtrics_raw)
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unite_exclusions.R
\name{unite_exclusions}
\alias{unite_exclusions}
\title{Unite multiple exclusion columns into single column}
\usage{
unite_exclusions(
  x,
  exclusion_types = c("duplicates", "duration", "ip", "location", "preview",
    "progress", "resolution"),
  separator = ",",
  remove = TRUE
)
}
\arguments{
\item{x}{Data frame or tibble (preferably exported from Qualtrics).}

\item{exclusion_types}{Vector of types of exclusions to unite.}

\item{separator}{Character string specifying what character to use to
separate multiple exclusion types}

\item{remove}{Logical specifying whether to remove united columns
(default = TRUE) or leave them in the data frame (FALSE)}
}
\value{
An object of the same type as \code{x} that includes the all of the same
rows but with a single \code{exclusion} column replacing all of the specified
\verb{exclusion_*} columns.
}
\description{
Each of the \verb{mark_*()} functions appends a new column to the data.
The \code{unite_exclusions()} function unites all of those columns in a
single column that can be used to filter any or all exclusions downstream.
Rows with multiple exclusions are concatenated with commas.
}
\examples{

# Unite all exclusion types
df <- qualtrics_text \%>\%
  mark_duplicates() \%>\%
  mark_duration(min_duration = 100) \%>\%
  mark_ip() \%>\%
  mark_location() \%>\%
  mark_preview() \%>\%
  mark_progress() \%>\%
  mark_resolution()
df2 <- df \%>\%
  unite_exclusions()

# Unite subset of exclusion types
df2 <- df \%>\%
  unite_exclusions(exclusion_types = c("duplicates", "duration", "ip"))
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duplicates.R
\name{mark_duplicates}
\alias{mark_duplicates}
\title{Mark duplicate IP addresses and/or locations}
\usage{
mark_duplicates(
  x,
  id_col = "ResponseId",
  ip_col = "IPAddress",
  location_col = c("LocationLatitude", "LocationLongitude"),
  dupl_ip = TRUE,
  dupl_location = TRUE,
  include_na = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{ip_col}{Column name for IP addresses.}

\item{location_col}{Two element vector specifying columns for latitude and
longitude (in that order).}

\item{dupl_ip}{Logical indicating whether to check IP addresses.}

\item{dupl_location}{Logical indicating whether to check latitude and
longitude.}

\item{include_na}{Logical indicating whether to include rows with NAs for
IP address and location as potentially excluded rows.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes a column marking rows
with duplicate IP addresses and/or locations.
For a function that just checks for and returns duplicate rows,
use \code{\link[=check_duplicates]{check_duplicates()}}. For a function that excludes these rows,
use \code{\link[=exclude_duplicates]{exclude_duplicates()}}.
}
\description{
The \code{mark_duplicates()} function creates a column labeling
rows of data that have the same IP address and/or same latitude and
longitude. The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
By default, IP address and location are both checked, but they can be
checked separately with the \code{dupl_ip} and \code{dupl_location} arguments.

The function outputs to console separate messages about the number of
rows with duplicate IP addresses and rows with duplicate locations.
These counts are computed independently, so rows may be counted for both
types of duplicates.
}
\examples{
# Mark duplicate IP addresses and locations
data(qualtrics_text)
df <- mark_duplicates(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_duplicates()

# Mark only for duplicate locations
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_duplicates(dupl_location = FALSE)
}
\seealso{
Other duplicates functions: 
\code{\link{check_duplicates}()},
\code{\link{exclude_duplicates}()}

Other mark functions: 
\code{\link{mark_duration}()},
\code{\link{mark_ip}()},
\code{\link{mark_location}()},
\code{\link{mark_preview}()},
\code{\link{mark_progress}()},
\code{\link{mark_resolution}()}
}
\concept{duplicates functions}
\concept{mark functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/progress.R
\name{mark_progress}
\alias{mark_progress}
\title{Mark survey progress}
\usage{
mark_progress(
  x,
  min_progress = 100,
  id_col = "ResponseId",
  finished_col = "Finished",
  progress_col = "Progress",
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{min_progress}{Amount of progress considered acceptable to include.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{finished_col}{Column name for whether survey was completed.}

\item{progress_col}{Column name for percentage of survey completed.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes a column marking rows
that have incomplete progress.
For a function that checks for these rows, use \code{\link[=check_progress]{check_progress()}}.
For a function that excludes these rows, use \code{\link[=exclude_progress]{exclude_progress()}}.
}
\description{
The \code{mark_progress()} function creates a column labeling
rows that have incomplete progress.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The default requires 100\% completion, but lower levels of completion
maybe acceptable and can be allowed by specifying the \code{min_progress}
argument.
The finished column in Qualtrics can be a numeric or character vector
depending on whether it is exported as choice text or numeric values.
This function works for both.

The function outputs to console a message about the number of rows
that have incomplete progress.
}
\examples{
# Mark rows with incomplete progress
data(qualtrics_text)
df <- mark_progress(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_progress()

# Include a lower acceptable completion percentage
df <- qualtrics_numeric \%>\%
  exclude_preview() \%>\%
  mark_progress(min_progress = 98)
}
\seealso{
Other progress functions: 
\code{\link{check_progress}()},
\code{\link{exclude_progress}()}

Other mark functions: 
\code{\link{mark_duplicates}()},
\code{\link{mark_duration}()},
\code{\link{mark_ip}()},
\code{\link{mark_location}()},
\code{\link{mark_preview}()},
\code{\link{mark_resolution}()}
}
\concept{mark functions}
\concept{progress functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unite_exclusions.R
\name{collapse_exclusions}
\alias{collapse_exclusions}
\title{Unite multiple exclusion columns into single column}
\usage{
collapse_exclusions(
  x,
  exclusion_types = c("duplicates", "duration", "ip", "location", "preview",
    "progress", "resolution"),
  separator = ",",
  remove = TRUE
)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

\code{collapse_exclusions()} was renamed to \code{\link[=unite_exclusions]{unite_exclusions()}} to create a more
consistent API with tidyverse's \code{unite()} function---please use
\code{unite_exclusions()}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{print_data}
\alias{print_data}
\title{Print data to console}
\usage{
print_data(x, print)
}
\arguments{
\item{x}{Data set to print or not}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\description{
Prints the data to the console. \emph{This function is not exported.}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duration.R
\name{mark_duration}
\alias{mark_duration}
\title{Mark minimum or maximum durations}
\usage{
mark_duration(
  x,
  min_duration = 10,
  max_duration = NULL,
  id_col = "ResponseId",
  duration_col = "Duration (in seconds)",
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{min_duration}{Minimum duration that is too fast in seconds.}

\item{max_duration}{Maximum duration that is too slow in seconds.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{duration_col}{Column name for durations.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes a column marking rows
with fast and slow duration.
For a function that checks for these rows, use \code{\link[=check_duration]{check_duration()}}.
For a function that excludes these rows, use \code{\link[=exclude_duration]{exclude_duration()}}.
}
\description{
The \code{mark_duration()} function creates a column labeling
rows with fast and/or slow duration.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
By default, minimum durations of 10 seconds are checked, but either
minima or maxima can be checked with the \code{min_duration} and
\code{max_duration} arguments. The function outputs to console separate
messages about the number of rows that are too fast or too slow.

This function returns the fast and slow rows.
}
\examples{
# Mark durations faster than 100 seconds
data(qualtrics_text)
df <- mark_duration(qualtrics_text, min_duration = 100)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_duration()

# Mark only for durations slower than 800 seconds
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_duration(max_duration = 800)
}
\seealso{
Other duration functions: 
\code{\link{check_duration}()},
\code{\link{exclude_duration}()}

Other mark functions: 
\code{\link{mark_duplicates}()},
\code{\link{mark_ip}()},
\code{\link{mark_location}()},
\code{\link{mark_preview}()},
\code{\link{mark_progress}()},
\code{\link{mark_resolution}()}
}
\concept{duration functions}
\concept{mark functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualtrics_data.R
\docType{data}
\name{qualtrics_text}
\alias{qualtrics_text}
\title{Example text-based metadata from simulated Qualtrics study}
\format{
A data frame with 100 rows and 16 variables:
\describe{
\item{StartDate}{date and time data collection started, in ISO 8601 format}
\item{EndDate}{date and time data collection ended, in ISO 8601 format}
\item{Status}{flag for preview (Survey Preview) vs. implemented survey
(IP Address) entries}
\item{IPAddress}{participant IP address (truncated for anonymity)}
\item{Progress}{percentage of survey completed}
\item{Duration (in seconds)}{duration of time required to complete survey,
in seconds}
\item{Finished}{logical for whether survey was completed (TRUE) or progress
was < 100 (FALSE)}
\item{RecordedDate}{date and time survey was recorded, in ISO 8601 format}
\item{ResponseId}{random ID for participants}
\item{LocationLatitude}{latitude geolocated from IP address}
\item{LocationLongitude}{longitude geolocated from IP address}
\item{UserLanguage}{language set in Qualtrics}
\item{Browser}{user web browser type}
\item{Version}{user web browser version}
\item{Operating System}{user operating system}
\item{Resolution}{user screen resolution}
}
}
\usage{
qualtrics_text
}
\description{
A dataset containing the metadata from a standard Qualtrics survey with
browser metadata collected and exported with "Use choice text".
These data were randomly generated using \code{\link[iptools:ip_random]{iptools::ip_random()}} and
\href{https://cran.r-project.org/package=rgeolocate}{rgeolocate::ip2location()} functions.
}
\seealso{
Other data: 
\code{\link{qualtrics_numeric}},
\code{\link{qualtrics_raw}}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{print_exclusion}
\alias{print_exclusion}
\title{Print number of excluded rows}
\usage{
print_exclusion(remaining_data, x, msg)
}
\arguments{
\item{remaining_data}{Data after removing exclusions.}

\item{x}{Original data before removing exclusions.}

\item{msg}{Text to describe what types of rows were excluded.}
}
\description{
Prints a message to the console with the number of excluded rows.
\emph{This function is not exported.}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/location.R
\name{check_location}
\alias{check_location}
\title{Check for locations outside of the US}
\usage{
check_location(
  x,
  id_col = "ResponseId",
  location_col = c("LocationLatitude", "LocationLongitude"),
  include_na = FALSE,
  keep = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{location_col}{Two element vector specifying columns for latitude
and longitude (in that order).}

\item{include_na}{Logical indicating whether to include rows with NA in
latitude and longitude columns in the output list of potentially excluded
data.}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
The output is a data frame of the rows that are located outside of
the US and (if \code{include_na == FALSE}) rows with no location information.
For a function that marks these rows, use \code{\link[=mark_location]{mark_location()}}.
For a function that excludes these rows, use \code{\link[=exclude_location]{exclude_location()}}.
}
\description{
The \code{check_location()} function subsets rows of data, retaining rows
that have locations outside of the US.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The function only works for the United States.
It uses the #' \code{\link[maps:map.where]{maps::map.where()}} to determine if latitude and longitude
are inside the US.

The function outputs to console a message about the number of rows
with locations outside of the US.
}
\examples{
# Check for locations outside of the US
data(qualtrics_text)
check_location(qualtrics_text)

# Remove preview data first
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_location()

# Do not print rows to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_location(print = FALSE)

# Do not print message to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_location(quiet = TRUE)
}
\seealso{
Other location functions: 
\code{\link{exclude_location}()},
\code{\link{mark_location}()}

Other check functions: 
\code{\link{check_duplicates}()},
\code{\link{check_duration}()},
\code{\link{check_ip}()},
\code{\link{check_preview}()},
\code{\link{check_progress}()},
\code{\link{check_resolution}()}
}
\concept{check functions}
\concept{location functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duration.R
\name{exclude_duration}
\alias{exclude_duration}
\title{Exclude rows with minimum or maximum durations}
\usage{
exclude_duration(
  x,
  min_duration = 10,
  max_duration = NULL,
  id_col = "ResponseId",
  duration_col = "Duration (in seconds)",
  quiet = TRUE,
  print = TRUE,
  silent = FALSE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{min_duration}{Minimum duration that is too fast in seconds.}

\item{max_duration}{Maximum duration that is too slow in seconds.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{duration_col}{Column name for durations.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}

\item{silent}{Logical indicating whether to print message to console. Note
this argument controls the exclude message not the check message.}
}
\value{
An object of the same type as \code{x} that excludes rows
with fast and/or slow duration.
For a function that checks for these rows, use \code{\link[=check_duration]{check_duration()}}.
For a function that marks these rows, use \code{\link[=mark_duration]{mark_duration()}}.
}
\description{
The \code{exclude_duration()} function removes
rows of data that have durations that are too fast or too slow.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
By default, minimum durations of 10 seconds are checked, but either
minima or maxima can be checked with the \code{min_duration} and
\code{max_duration} arguments. The function outputs to console separate
messages about the number of rows that are too fast or too slow.

This function returns the fast and slow rows.
}
\examples{
# Exclude durations faster than 100 seconds
data(qualtrics_text)
df <- exclude_duration(qualtrics_text, min_duration = 100)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_duration()

# Exclude only for durations slower than 800 seconds
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_duration(max_duration = 800)
}
\seealso{
Other duration functions: 
\code{\link{check_duration}()},
\code{\link{mark_duration}()}

Other exclude functions: 
\code{\link{exclude_duplicates}()},
\code{\link{exclude_ip}()},
\code{\link{exclude_location}()},
\code{\link{exclude_preview}()},
\code{\link{exclude_progress}()},
\code{\link{exclude_resolution}()}
}
\concept{duration functions}
\concept{exclude functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mark_rows}
\alias{mark_rows}
\title{Return marked rows}
\usage{
mark_rows(x, filtered_data, id_col, exclusion_type)
}
\arguments{
\item{x}{Original data.}

\item{filtered_data}{Data to be excluded.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{exclusion_type}{Column name for exclusion column.}
}
\description{
Create new column marking rows that meet exclusion criteria.
\emph{This function is not exported.}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ip.R
\name{mark_ip}
\alias{mark_ip}
\title{Mark IP addresses from outside of a specified country.}
\usage{
mark_ip(
  x,
  id_col = "ResponseId",
  ip_col = "IPAddress",
  country = "US",
  include_na = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame or tibble (preferably imported from Qualtrics using
\{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{ip_col}{Column name for IP addresses.}

\item{country}{Two-letter abbreviation of country to check (default is "US").}

\item{include_na}{Logical indicating whether to include rows with NA in
IP address column in the output list of potentially excluded data.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes a column marking rows
with IP addresses outside of the specified country.
For a function that checks these rows, use \code{\link[=check_ip]{check_ip()}}.
For a function that excludes these rows, use \code{\link[=exclude_ip]{exclude_ip()}}.
}
\description{
The \code{mark_ip()} function creates a column labeling
rows of data that have IP addresses from outside the specified country.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The function uses \code{\link[iptools:country_ranges]{iptools::country_ranges()}} to assign IP addresses to
specific countries using
\href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}{ISO 3166-1 alpha-2 country codes}.

The function outputs to console a message about the number of rows
with IP addresses outside of the specified country. If there are \code{NA}s for IP
addresses (likely due to including preview data---see \code{\link[=check_preview]{check_preview()}}), it
will print a message alerting to the number of rows with \code{NA}s.
}
\note{
This function \strong{requires internet connectivity} as it uses the
\code{\link[iptools:country_ranges]{iptools::country_ranges()}} function, which pulls daily updated data
from \url{https://www.iwik.org/ipcountry/}. It only updates the data once
per session, as it caches the results for future work during the session.
}
\examples{
# Mark IP addresses outside of the US
data(qualtrics_text)
df <- mark_ip(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_ip()

# Mark IP addresses outside of Germany
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  mark_ip(country = "DE")
}
\seealso{
Other ip functions: 
\code{\link{check_ip}()},
\code{\link{exclude_ip}()}

Other mark functions: 
\code{\link{mark_duplicates}()},
\code{\link{mark_duration}()},
\code{\link{mark_location}()},
\code{\link{mark_preview}()},
\code{\link{mark_progress}()},
\code{\link{mark_resolution}()}
}
\concept{ip functions}
\concept{mark functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/excluder-package.R
\docType{package}
\name{excluder-package}
\alias{excluder}
\alias{excluder-package}
\title{excluder: Checks for Exclusion Criteria in Online Data}
\description{
Data that are collected through online sources such as Mechanical Turk may require excluding rows because of IP address duplication, geolocation, or completion duration. This package facilitates exclusion of these data for Qualtrics datasets.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/excluder/}
  \item \url{https://github.com/ropensci/excluder/}
  \item Report bugs at \url{https://github.com/ropensci/excluder/issues/}
}

}
\author{
\strong{Maintainer}: Jeffrey R. Stevens \email{jeffrey.r.stevens@gmail.com} (\href{https://orcid.org/0000-0003-2375-1360}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Joseph O'Brien (\href{https://orcid.org/0000-0001-9851-5077}{ORCID}) [reviewer]
  \item Julia Silge \email{julia.silge@gmail.com} (\href{https://orcid.org/0000-0002-3671-836X}{ORCID}) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duplicates.R
\name{check_duplicates}
\alias{check_duplicates}
\title{Check for duplicate IP addresses and/or locations}
\usage{
check_duplicates(
  x,
  id_col = "ResponseId",
  ip_col = "IPAddress",
  location_col = c("LocationLatitude", "LocationLongitude"),
  dupl_ip = TRUE,
  dupl_location = TRUE,
  include_na = FALSE,
  keep = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{ip_col}{Column name for IP addresses.}

\item{location_col}{Two element vector specifying columns for latitude and
longitude (in that order).}

\item{dupl_ip}{Logical indicating whether to check IP addresses.}

\item{dupl_location}{Logical indicating whether to check latitude and
longitude.}

\item{include_na}{Logical indicating whether to include rows with NAs for
IP address and location as potentially excluded rows.}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes the rows with
duplicate IP addresses and/or locations. This includes a column
called dupe_count that returns the number of duplicates.
For a function that marks these rows, use \code{\link[=mark_duplicates]{mark_duplicates()}}.
For a function that excludes these rows, use \code{\link[=exclude_duplicates]{exclude_duplicates()}}.
}
\description{
The \code{check_duplicates()} function subsets rows of data, retaining rows
that have the same IP address and/or same latitude and longitude. The
function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
By default, IP address and location are both checked, but they can be
checked separately with the \code{dupl_ip} and \code{dupl_location} arguments.

The function outputs to console separate messages about the number of
rows with duplicate IP addresses and rows with duplicate locations.
These counts are computed independently, so rows may be counted for both
types of duplicates.
}
\examples{
# Check for duplicate IP addresses and locations
data(qualtrics_text)
check_duplicates(qualtrics_text)

# Check only for duplicate locations
qualtrics_text \%>\%
  check_duplicates(dupl_location = FALSE)

# Do not print rows to console
qualtrics_text \%>\%
  check_duplicates(print = FALSE)

# Do not print message to console
qualtrics_text \%>\%
  check_duplicates(quiet = TRUE)
}
\seealso{
Other duplicates functions: 
\code{\link{exclude_duplicates}()},
\code{\link{mark_duplicates}()}

Other check functions: 
\code{\link{check_duration}()},
\code{\link{check_ip}()},
\code{\link{check_location}()},
\code{\link{check_preview}()},
\code{\link{check_progress}()},
\code{\link{check_resolution}()}
}
\concept{check functions}
\concept{duplicates functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualtrics_data.R
\docType{data}
\name{qualtrics_numeric}
\alias{qualtrics_numeric}
\title{Example numeric metadata from simulated Qualtrics study}
\format{
A data frame with 100 rows and 16 variables:
\describe{
\item{StartDate}{date and time data collection started, in ISO 8601 format}
\item{EndDate}{date and time data collection ended, in ISO 8601 format}
\item{Status}{numeric flag for preview (1) vs. implemented survey (0)
entries}
\item{IPAddress}{participant IP address (truncated for anonymity)}
\item{Progress}{percentage of survey completed}
\item{Duration (in seconds)}{duration of time required to complete survey,
in seconds}
\item{Finished}{numeric flag for whether survey was completed (1) or
progress was < 100 (0)}
\item{RecordedDate}{date and time survey was recorded, in ISO 8601 format}
\item{ResponseId}{random ID for participants}
\item{LocationLatitude}{latitude geolocated from IP address}
\item{LocationLongitude}{longitude geolocated from IP address}
\item{UserLanguage}{language set in Qualtrics}
\item{Browser}{user web browser type}
\item{Version}{user web browser version}
\item{Operating System}{user operating system}
\item{Resolution}{user screen resolution}
}
}
\usage{
qualtrics_numeric
}
\description{
A dataset containing the metadata from a standard Qualtrics survey with
browser metadata collected and exported with "Use numeric values".
These data were randomly generated using \code{\link[iptools:ip_random]{iptools::ip_random()}} and
\href{https://cran.r-project.org/package=rgeolocate}{rgeolocate::ip2location()} functions.
}
\seealso{
Other data: 
\code{\link{qualtrics_raw}},
\code{\link{qualtrics_text}}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\arguments{
\item{lhs}{A value or the magrittr placeholder.}

\item{rhs}{A function call using the magrittr semantics.}
}
\value{
The result of calling \code{rhs(lhs)}.
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preview.R
\name{check_preview}
\alias{check_preview}
\title{Check for survey previews}
\usage{
check_preview(
  x,
  id_col = "ResponseId",
  preview_col = "Status",
  keep = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{preview_col}{Column name for survey preview.}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
The output is a data frame of the rows
that are survey previews.
For a function that marks these rows, use \code{\link[=mark_preview]{mark_preview()}}.
For a function that excludes these rows, use \code{\link[=exclude_preview]{exclude_preview()}}.
}
\description{
The \code{check_preview()} function subsets rows of data, retaining rows
that are survey previews.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The preview column in Qualtrics can be a numeric or character vector
depending on whether it is exported as choice text or numeric values.
This function works for both.

The function outputs to console a message about the number of rows
that are survey previews.
}
\examples{
# Check for survey previews
data(qualtrics_text)
check_preview(qualtrics_text)

# Works for Qualtrics data exported as numeric values, too
qualtrics_numeric \%>\%
  check_preview()

# Do not print rows to console
qualtrics_text \%>\%
  check_preview(print = FALSE)

# Do not print message to console
qualtrics_text \%>\%
  check_preview(quiet = TRUE)
}
\seealso{
Other preview functions: 
\code{\link{exclude_preview}()},
\code{\link{mark_preview}()}

Other check functions: 
\code{\link{check_duplicates}()},
\code{\link{check_duration}()},
\code{\link{check_ip}()},
\code{\link{check_location}()},
\code{\link{check_progress}()},
\code{\link{check_resolution}()}
}
\concept{check functions}
\concept{preview functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deidentify.R
\name{deidentify}
\alias{deidentify}
\title{Remove columns that could include identifiable information}
\usage{
deidentify(x, strict = TRUE)
}
\arguments{
\item{x}{Data frame (downloaded from Qualtrics).}

\item{strict}{Logical indicating whether to use strict or non-strict level
of deidentification. Strict removes computer information columns in addition
to IP address and location.}
}
\value{
An object of the same type as \code{x} that excludes Qualtrics columns with
identifiable information.
}
\description{
The \code{deidentify()} function selects out columns from
\href{https://www.qualtrics.com/}{Qualtrics} surveys that may include identifiable
information such as IP address, location, or computer characteristics.
}
\details{
The function offers two levels of deidentification. The default strict level
removes columns associated with IP address and location and computer
information (browser type and version, operating system, and screen
resolution). The non-strict level removes only columns associated with
IP address and location.

Typically, deidentification should be used at the end of a processing pipeline
so that these columns can be used to exclude rows.
}
\examples{
names(qualtrics_numeric)

# Remove IP address, location, and computer information columns
deid <- deidentify(qualtrics_numeric)
names(deid)

# Remove only IP address and location columns
deid2 <- deidentify(qualtrics_numeric, strict = FALSE)
names(deid2)
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preview.R
\name{mark_preview}
\alias{mark_preview}
\title{Mark survey previews}
\usage{
mark_preview(
  x,
  id_col = "ResponseId",
  preview_col = "Status",
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{preview_col}{Column name for survey preview.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
An object of the same type as \code{x} that includes a column marking rows
that are survey previews.
For a function that checks for these rows, use \code{\link[=check_preview]{check_preview()}}.
For a function that excludes these rows, use \code{\link[=exclude_preview]{exclude_preview()}}.
}
\description{
The \code{mark_preview()} function creates a column labeling
rows that are survey previews.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The preview column in Qualtrics can be a numeric or character vector
depending on whether it is exported as choice text or numeric values.
This function works for both.

The function outputs to console a message about the number of rows
that are survey previews.
}
\examples{
# Mark survey previews
data(qualtrics_text)
df <- mark_preview(qualtrics_text)

# Works for Qualtrics data exported as numeric values, too
df <- qualtrics_numeric \%>\%
  mark_preview()
}
\seealso{
Other preview functions: 
\code{\link{check_preview}()},
\code{\link{exclude_preview}()}

Other mark functions: 
\code{\link{mark_duplicates}()},
\code{\link{mark_duration}()},
\code{\link{mark_ip}()},
\code{\link{mark_location}()},
\code{\link{mark_progress}()},
\code{\link{mark_resolution}()}
}
\concept{mark functions}
\concept{preview functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualtrics_data.R
\docType{data}
\name{qualtrics_raw}
\alias{qualtrics_raw}
\title{Example text-based metadata from simulated Qualtrics study}
\format{
A data frame with 102 rows and 16 variables:
\describe{
\item{StartDate}{date and time data collection started, in ISO 8601 format}
\item{EndDate}{date and time data collection ended, in ISO 8601 format}
\item{Status}{flag for preview (Survey Preview) vs. implemented survey
(IP Address) entries}
\item{IPAddress}{participant IP address (truncated for anonymity)}
\item{Progress}{percentage of survey completed}
\item{Duration (in seconds)}{duration of time required to complete survey,
in seconds}
\item{Finished}{logical for whether survey was completed (TRUE) or progress
was < 100 (FALSE)}
\item{RecordedDate}{date and time survey was recorded, in ISO 8601 format}
\item{ResponseId}{random ID for participants}
\item{LocationLatitude}{latitude geolocated from IP address}
\item{LocationLongitude}{longitude geolocated from IP address}
\item{UserLanguage}{language set in Qualtrics}
\item{Browser}{user web browser type}
\item{Version}{user web browser version}
\item{Operating System}{user operating system}
\item{Resolution}{user screen resolution}
}
}
\usage{
qualtrics_raw
}
\description{
A dataset containing the metadata from a standard Qualtrics survey with
browser metadata collected and exported with "Use choice text".
These data were randomly generated using \code{\link[iptools:ip_random]{iptools::ip_random()}} and
\href{https://cran.r-project.org/package=rgeolocate}{rgeolocate::ip2location()} functions.
This dataset includes the two header rows of with column information that is
exported by Qualtrics.
}
\seealso{
Other data: 
\code{\link{qualtrics_numeric}},
\code{\link{qualtrics_text}}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/progress.R
\name{check_progress}
\alias{check_progress}
\title{Check for survey progress}
\usage{
check_progress(
  x,
  min_progress = 100,
  id_col = "ResponseId",
  finished_col = "Finished",
  progress_col = "Progress",
  keep = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{min_progress}{Amount of progress considered acceptable to include.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{finished_col}{Column name for whether survey was completed.}

\item{progress_col}{Column name for percentage of survey completed.}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
The output is a data frame of the rows
that have incomplete progress.
For a function that marks these rows, use \code{\link[=mark_progress]{mark_progress()}}.
For a function that excludes these rows, use \code{\link[=exclude_progress]{exclude_progress()}}.
}
\description{
The \code{check_progress()} function subsets rows of data, retaining rows
that have incomplete progress.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The default requires 100\% completion, but lower levels of completion
maybe acceptable and can be allowed by specifying the \code{min_progress}
argument.
The finished column in Qualtrics can be a numeric or character vector
depending on whether it is exported as choice text or numeric values.
This function works for both.

The function outputs to console a message about the number of rows
that have incomplete progress.
}
\examples{
# Check for rows with incomplete progress
data(qualtrics_text)
check_progress(qualtrics_text)

# Remove preview data first
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_progress()

# Include a lower acceptable completion percentage
qualtrics_numeric \%>\%
  exclude_preview() \%>\%
  check_progress(min_progress = 98)

# Do not print rows to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_progress(print = FALSE)

# Do not print message to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_progress(quiet = TRUE)
}
\seealso{
Other progress functions: 
\code{\link{exclude_progress}()},
\code{\link{mark_progress}()}

Other check functions: 
\code{\link{check_duplicates}()},
\code{\link{check_duration}()},
\code{\link{check_ip}()},
\code{\link{check_location}()},
\code{\link{check_preview}()},
\code{\link{check_resolution}()}
}
\concept{check functions}
\concept{progress functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resolution.R
\name{check_resolution}
\alias{check_resolution}
\title{Check screen resolution}
\usage{
check_resolution(
  x,
  width_min = 1000,
  height_min = 0,
  id_col = "ResponseId",
  res_col = "Resolution",
  keep = FALSE,
  quiet = FALSE,
  print = TRUE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{width_min}{Minimum acceptable screen width.}

\item{height_min}{Minimum acceptable screen height.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{res_col}{Column name for screen resolution (in format widthxheight).}

\item{keep}{Logical indicating whether to keep or remove exclusion column.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}
}
\value{
The output is a data frame of the rows that have unacceptable screen
resolutions. This includes new columns for resolution width and height.
For a function that marks these rows, use \code{\link[=mark_resolution]{mark_resolution()}}.
For a function that excludes these rows, use \code{\link[=exclude_resolution]{exclude_resolution()}}.
}
\description{
The \code{check_resolution()} function subsets rows of data, retaining rows
that have unacceptable screen resolution. This can be used, for example, to
determine data collected via phones when desktop monitors are required.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.

The function outputs to console a message about the number of rows
with unacceptable screen resolution.
}
\examples{
# Check for survey previews
data(qualtrics_text)
check_resolution(qualtrics_text)

# Remove preview data first
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_resolution()

# Do not print rows to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_resolution(print = FALSE)

# Do not print message to console
qualtrics_text \%>\%
  exclude_preview() \%>\%
  check_resolution(quiet = TRUE)
}
\seealso{
Other resolution functions: 
\code{\link{exclude_resolution}()},
\code{\link{mark_resolution}()}

Other check functions: 
\code{\link{check_duplicates}()},
\code{\link{check_duration}()},
\code{\link{check_ip}()},
\code{\link{check_location}()},
\code{\link{check_preview}()},
\code{\link{check_progress}()}
}
\concept{check functions}
\concept{resolution functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resolution.R
\name{exclude_resolution}
\alias{exclude_resolution}
\title{Exclude unacceptable screen resolution}
\usage{
exclude_resolution(
  x,
  width_min = 1000,
  height_min = 0,
  id_col = "ResponseId",
  res_col = "Resolution",
  quiet = TRUE,
  print = TRUE,
  silent = FALSE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{width_min}{Minimum acceptable screen width.}

\item{height_min}{Minimum acceptable screen height.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{res_col}{Column name for screen resolution (in format widthxheight).}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}

\item{silent}{Logical indicating whether to print message to console. Note
this argument controls the exclude message not the check message.}
}
\value{
An object of the same type as \code{x} that excludes rows
that have unacceptable screen resolutions.
For a function that checks for these rows, use \code{\link[=check_resolution]{check_resolution()}}.
For a function that marks these rows, use \code{\link[=mark_resolution]{mark_resolution()}}.
}
\description{
The \code{exclude_resolution()} function removes
rows that have unacceptable screen resolution.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.

The function outputs to console a message about the number of rows
with unacceptable screen resolution.
}
\examples{
# Exclude low screen resolutions
data(qualtrics_text)
df <- exclude_resolution(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_resolution()
}
\seealso{
Other resolution functions: 
\code{\link{check_resolution}()},
\code{\link{mark_resolution}()}

Other exclude functions: 
\code{\link{exclude_duplicates}()},
\code{\link{exclude_duration}()},
\code{\link{exclude_ip}()},
\code{\link{exclude_location}()},
\code{\link{exclude_preview}()},
\code{\link{exclude_progress}()}
}
\concept{exclude functions}
\concept{resolution functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/progress.R
\name{exclude_progress}
\alias{exclude_progress}
\title{Exclude survey progress}
\usage{
exclude_progress(
  x,
  min_progress = 100,
  id_col = "ResponseId",
  finished_col = "Finished",
  progress_col = "Progress",
  quiet = TRUE,
  print = TRUE,
  silent = FALSE
)
}
\arguments{
\item{x}{Data frame (preferably imported from Qualtrics using \{qualtRics\}).}

\item{min_progress}{Amount of progress considered acceptable to include.}

\item{id_col}{Column name for unique row ID (e.g., participant).}

\item{finished_col}{Column name for whether survey was completed.}

\item{progress_col}{Column name for percentage of survey completed.}

\item{quiet}{Logical indicating whether to print message to console.}

\item{print}{Logical indicating whether to print returned tibble to
console.}

\item{silent}{Logical indicating whether to print message to console. Note
this argument controls the exclude message not the check message.}
}
\value{
An object of the same type as \code{x} that excludes rows
that have incomplete progress.
For a function that checks for these rows, use \code{\link[=check_progress]{check_progress()}}.
For a function that marks these rows, use \code{\link[=mark_progress]{mark_progress()}}.
}
\description{
The \code{exclude_progress()} function removes
rows that have incomplete progress.
The function is written to work with data from
\href{https://www.qualtrics.com/}{Qualtrics} surveys.
}
\details{
Default column names are set based on output from the
\href{https://docs.ropensci.org/qualtRics/reference/fetch_survey.html}{\code{qualtRics::fetch_survey()}}.
The default requires 100\% completion, but lower levels of completion
maybe acceptable and can be allowed by specifying the \code{min_progress}
argument.
The finished column in Qualtrics can be a numeric or character vector
depending on whether it is exported as choice text or numeric values.
This function works for both.

The function outputs to console a message about the number of rows
that have incomplete progress.
}
\examples{
# Exclude rows with incomplete progress
data(qualtrics_text)
df <- exclude_progress(qualtrics_text)

# Remove preview data first
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_progress()

# Include a lower acceptable completion percentage
df <- qualtrics_numeric \%>\%
  exclude_preview() \%>\%
  exclude_progress(min_progress = 98)

# Do not print rows to console
df <- qualtrics_text \%>\%
  exclude_preview() \%>\%
  exclude_progress(print = FALSE)
}
\seealso{
Other progress functions: 
\code{\link{check_progress}()},
\code{\link{mark_progress}()}

Other exclude functions: 
\code{\link{exclude_duplicates}()},
\code{\link{exclude_duration}()},
\code{\link{exclude_ip}()},
\code{\link{exclude_location}()},
\code{\link{exclude_preview}()},
\code{\link{exclude_resolution}()}
}
\concept{exclude functions}
\concept{progress functions}
