
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
