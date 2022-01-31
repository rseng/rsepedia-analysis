# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

<!-- README.md is generated from README.Rmd. Please edit that file -->

# visdat <img src="man/figures/visdat-logo.png" align="right" />

<!-- badges: start -->

[![rOpenSci
Badge](https://badges.ropensci.org/87_status.svg)](https://github.com/ropensci/software-review/issues/87)[![JOSS
status](https://joss.theoj.org/papers/10.21105/joss.00355/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00355)[![DOI](https://zenodo.org/badge/50553382.svg)](https://zenodo.org/badge/latestdoi/50553382)[![R-CMD-check](https://github.com/ropensci/visdat/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/visdat/actions)[![Codecov
test
coverage](https://codecov.io/gh/ropensci/visdat/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/visdat?branch=master)[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/visdat)](https://cran.r-project.org/package=visdat)[![CRAN
Logs](http://cranlogs.r-pkg.org/badges/visdat)](https://CRAN.R-project.org/package=visdat)[![Project
Status: Active – The project has reached a stable, usable state and is
being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
<!-- badges: end -->

# How to install

visdat is available on CRAN

``` r
install.packages("visdat")
```

If you would like to use the development version, install from github
with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/visdat")
```

# What does visdat do?

Initially inspired by
[`csv-fingerprint`](https://github.com/setosa/csv-fingerprint),
`vis_dat` helps you visualise a dataframe and “get a look at the data”
by displaying the variable classes in a dataframe as a plot with
`vis_dat`, and getting a brief look into missing data patterns using
`vis_miss`.

`visdat` has 6 functions:

-   `vis_dat()` visualises a dataframe showing you what the classes of
    the columns are, and also displaying the missing data.

-   `vis_miss()` visualises just the missing data, and allows for
    missingness to be clustered and columns rearranged. `vis_miss()` is
    similar to `missing.pattern.plot` from the
    [`mi`](https://CRAN.R-project.org/package=mi) package. Unfortunately
    `missing.pattern.plot` is no longer in the `mi` package (as of
    14/02/2016).

-   `vis_compare()` visualise differences between two dataframes of the
    same dimensions

-   `vis_expect()` visualise where certain conditions hold true in your
    data

-   `vis_cor()` visualise the correlation of variables in a nice heatmap

-   `vis_guess()` visualise the individual class of each value in your
    data

-   `vis_value()` visualise the value class of each cell in your data

-   `vis_binary()` visualise the occurrence of binary values in your
    data

You can read more about visdat in the vignette, [“using
visdat”](http://visdat.njtierney.com/articles/using_visdat.html).

## Code of Conduct

Please note that the visdat project is released with a [Contributor Code
of
Conduct](https://github.com/ropensci/visdat/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

# Examples

## Using `vis_dat()`

Let’s see what’s inside the `airquality` dataset from base R, which
contains information about daily air quality measurements in New York
from May to September 1973. More information about the dataset can be
found with `?airquality`.

``` r
library(visdat)

vis_dat(airquality)
```

![](man/figures/README-vis-dat-aq-1.png)<!-- -->

The plot above tells us that R reads this dataset as having numeric and
integer values, with some missing data in `Ozone` and `Solar.R`. The
classes are represented on the legend, and missing data represented by
grey. The column/variable names are listed on the x axis.

## Using `vis_miss()`

We can explore the missing data further using `vis_miss()`:

``` r
vis_miss(airquality)
```

![](man/figures/README-vis-miss-aq-1.png)<!-- -->

Percentages of missing/complete in `vis_miss` are accurate to 1 decimal
place.

You can cluster the missingness by setting `cluster = TRUE`:

``` r
vis_miss(airquality, 
         cluster = TRUE)
```

![](man/figures/README-vis-miss-aq-cluster-1.png)<!-- -->

Columns can also be arranged by columns with most missingness, by
setting `sort_miss = TRUE`:

``` r
vis_miss(airquality,
         sort_miss = TRUE)
```

![](man/figures/README-vis-miss-aq-sort-miss-1.png)<!-- -->

`vis_miss` indicates when there is a very small amount of missing data
at \<0.1% missingness:

``` r
test_miss_df <- data.frame(x1 = 1:10000,
                           x2 = rep("A", 10000),
                           x3 = c(rep(1L, 9999), NA))

vis_miss(test_miss_df)
```

![](man/figures/README-vis-miss-test-1.png)<!-- -->

`vis_miss` will also indicate when there is no missing data at all:

``` r
vis_miss(mtcars)
```

![](man/figures/README-vis-miss-mtcars-1.png)<!-- -->

To further explore the missingness structure in a dataset, I recommend
the [`naniar`](https://github.com/njtierney/naniar) package, which
provides more general tools for graphical and numerical exploration of
missing values.

## Using `vis_compare()`

Sometimes you want to see what has changed in your data. `vis_compare()`
displays the differences in two dataframes of the same size. Let’s look
at an example.

Let’s make some changes to the `chickwts`, and compare this new dataset:

``` r
set.seed(2019-04-03-1105)
chickwts_diff <- chickwts
chickwts_diff[sample(1:nrow(chickwts), 30),sample(1:ncol(chickwts), 2)] <- NA

vis_compare(chickwts_diff, chickwts)
```

![](man/figures/README-vis-compare-iris-1.png)<!-- -->

Here the differences are marked in blue.

If you try and compare differences when the dimensions are different,
you get an ugly error:

``` r
chickwts_diff_2 <- chickwts
chickwts_diff_2$new_col <- chickwts_diff_2$weight*2

vis_compare(chickwts, chickwts_diff_2)
# Error in vis_compare(chickwts, chickwts_diff_2) : 
#   Dimensions of df1 and df2 are not the same. vis_compare requires dataframes of identical dimensions.
```

## Using `vis_expect()`

`vis_expect` visualises certain conditions or values in your data. For
example, If you are not sure whether to expect values greater than 25 in
your data (airquality), you could write:
`vis_expect(airquality, ~.x>=25)`, and you can see if there are times
where the values in your data are greater than or equal to 25:

``` r
vis_expect(airquality, ~.x >= 25)
```

![](man/figures/README-vis-expect-1.png)<!-- -->

This shows the proportion of times that there are values greater than
25, as well as the missings.

## Using `vis_cor()`

To make it easy to plot correlations of your data, use `vis_cor`:

``` r
vis_cor(airquality)
```

![](man/figures/README-vis-cor-1.png)<!-- -->

## Using `vis_value`

`vis_value()` visualises the values of your data on a 0 to 1 scale.

``` r
vis_value(airquality)
```

![](man/figures/README-vis-value-1.png)<!-- -->

It only works on numeric data, so you might get strange results if you
are using factors:

``` r
library(ggplot2)
vis_value(iris)
```

    data input can only contain numeric values, please subset the data to the numeric values you would like. dplyr::select_if(data, is.numeric) can be helpful here!

So you might need to subset the data beforehand like so:

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

iris %>%
  select_if(is.numeric) %>%
  vis_value()
```

![](man/figures/README-iris-error-fix-1.png)<!-- -->

## Using `vis_binary()`

`vis_binary()` visualises binary values. See below for use with example
data, `dat_bin`

``` r
vis_binary(dat_bin)
```

![](man/figures/README-vis-bin-1.png)<!-- -->

If you don’t have only binary values a warning will be shown.

``` r
vis_binary(airquality)
```

    Error in test_if_all_binary(data) : 
      data input can only contain binary values - this means either 0 or 1, or NA. Please subset the data to be binary values, or see ?vis_value.

## Using `vis_guess()`

`vis_guess()` takes a guess at what each cell is. It’s best illustrated
using some messy data, which we’ll make here:

``` r
messy_vector <- c(TRUE,
                  T,
                  "TRUE",
                  "T",
                  "01/01/01",
                  "01/01/2001",
                  NA,
                  NaN,
                  "NA",
                  "Na",
                  "na",
                  "10",
                  10,
                  "10.1",
                  10.1,
                  "abc",
                  "$%TG")

set.seed(2019-04-03-1106)
messy_df <- data.frame(var1 = messy_vector,
                       var2 = sample(messy_vector),
                       var3 = sample(messy_vector))
```

``` r
vis_guess(messy_df)
vis_dat(messy_df)
```

<img src="man/figures/README-vis-guess-messy-df-1.png" width="50%" /><img src="man/figures/README-vis-guess-messy-df-2.png" width="50%" />

So here we see that there are many different kinds of data in your
dataframe. As an analyst this might be a depressing finding. We can see
this comparison above.

# Thank yous

Thank you to Ivan Hanigan who [first
commented](https://www.njtierney.com/post/2015/11/12/ggplot-missing-data/)
this suggestion after I made a blog post about an initial prototype
`ggplot_missing`, and Jenny Bryan, whose
[tweet](https://twitter.com/JennyBryan/status/679011378414268416) got me
thinking about `vis_dat`, and for her code contributions that removed a
lot of errors.

Thank you to Hadley Wickham for suggesting the use of the internals of
`readr` to make `vis_guess` work. Thank you to Miles McBain for his
suggestions on how to improve `vis_guess`. This resulted in making it at
least 2-3 times faster. Thanks to Carson Sievert for writing the code
that combined `plotly` with `visdat`, and for Noam Ross for suggesting
this in the first place. Thank you also to Earo Wang and Stuart Lee for
their help in getting capturing expressions in `vis_expect`.

Finally thank you to [rOpenSci](https://github.com/ropensci) and it’s
amazing [onboarding
process](https://github.com/ropensci/software-review), this process has
made visdat a much better package, thanks to the editor Noam Ross
(@noamross), and the reviewers Sean Hughes (@seaaan) and Mara Averick
(@batpigandme).

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# visdat 0.6.0.9000 (2021/07/06) "Superman, Lazlo Bane"

## New Feature

* `vis_dat()` `vis_miss()` and `vis_guess()` now render missing values in list-columns (@cregouby #138)
* A new vignette on 

# visdat 0.6.0 (2021/07/05) "Shibuya, Covet, San Holo"

## Bug Fix

* output of plot in `vis_expect` would reorder columns ([#133](https://github.com/ropensci/visdat/issues/133)), fixed in [#143](https://github.com/ropensci/visdat/pull/134) by [@muschellij2](https://github.com/muschellij2).

## New Feature

* `vis_value()` for visualising all values in a dataset. It rescales values to be between 0 and 1. See #100
* `vis_binary()` for visualising datasets with binary values - similar to `vis_value()`, but just for binary data (0, 1, NA). See #112. Thank you to Trish Gilholm for her suggested use case for this.

# visdat 0.5.3 (2019/02/04) "The Legend of LoFi"

## Minor Change

* Update `vis_cor()` to use perceptually uniform colours from `scico` package, using `scico::scico(3, palette = "vik")`.
* Update `vis_cor()` to have fixed legend values from -1 to +1 (#110) using options `breaks` and `limits`. Special thanks to [this SO thread for the answer](https://stackoverflow.com/questions/24265652/label-minimum-and-maximum-of-scale-fill-gradient-legend-with-text-ggplot2)
* Uses `glue` and `glue_collapse()` instead of `paste` and `paste0`
* adds WORDLIST for spelling thanks to `usethis::use_spell_check()`

# visdat 0.5.2 (2018/12/06) "Youth, The Midnight, Kids"

## Minor Change

* Internal error message has been improved by [Nic](https://github.com/thisisnic) in [#102](https://github.com/ropensci/visdat/pull/102)

## Bug Fix

* [Jim Hester](https://github.com/jimhester) fixed recent changes in readr 1.2.0 in PR [#103](https://github.com/ropensci/visdat/pull/103), which changes the default behavior of the `guess_parser`, to not
guess integer types by default. To opt-into the current behavior you
need to pass `guess_integer = TRUE.`

# visdat 0.5.1 (2018/07/02) "The Northern Lights Moonwalker"

## New Feature

* `vis_compare()` for comparing two dataframes of the same dimensions
* `vis_expect()` for visualising where certain values of expectations occur in the data
    * Added NA colours to `vis_expect`
    * Added `show_perc` arg to `vis_expect` to show the percentage of expectations that are TRUE. #73
* `vis_cor` to visualise correlations in a dataframe
* `vis_guess()` for displaying the likely type for each cell in a dataframe
* Added draft `vis_expect` to make it easy to look at certain appearances of numbers in your data.
* visdat is now under the rOpenSci github repository

## Minor Changes

* added CITATION for visdat to cite the JOSS article
* updated options for `vis_cor` to use argument `na_action` not `use_op`.
* cleaned up the organisation of the files and internal functions
* Added appropriate legend and x axis for `vis_miss_ly` - thanks to Stuart Lee
* Updated the `paper.md` for JOSS
* Updated some old links in doco
* Added Sean Hughes and Mara Averick to the DESCRIPTION with `ctb`.
* Minor changes to the paper for JOSS

## Bug Fixes

* Fix bug reported in [#75](https://github.com/ropensci/visdat/issues/75) 
  where `vis_dat(diamonds)` errored `seq_len(nrow(x))` inside internal 
  function `vis_gather_`, used to calculate the row numbers. Using 
  `mutate(rows = dplyr::row_number())` solved the issue.

* Fix bug reported in [#72](https://github.com/ropensci/visdat/issues/72)
  where `vis_miss` errored when one column was given to it. This was an issue
  with using `limits` inside `scale_x_discrete` - which is used to order the
  columns of the data. It is not necessary to order one column of data, so I
  created an if-else to avoid this step and return the plot early.

* Fix visdat x axis alignment when show_perc_col = FALSE - [#82](https://github.com/ropensci/visdat/issues/82)

* fix visdat x axis alignment - [issue 57](https://github.com/ropensci/visdat/issues/57)
* fix bug where the column percentage missing would print to be NA when it was exactly equal to 0.1% missing. - [issue 62](https://github.com/ropensci/visdat/issues/62)
* `vis_cor` didn't gather variables for plotting appropriately - now fixed

# visdat 0.1.0 (2017/07/03) ("JOSS")

- lightweight CRAN submission - will only contain functions `vis_dat` and `vis_miss`

# visdat 0.0.7.9100 (2017/07/03)

## New Features

- `add_vis_dat_pal()` (internal) to add a palette for `vis_dat` and `vis_guess`
- `vis_guess` now gets a palette argument like `vis_dat`
- Added protoype/placeholder functions for `plotly` vis_*_ly interactive graphs:
  - `vis_guess_ly()`
  - `vis_dat_ly()`
  - `vis_compare_ly()`
  These simply wrap `plotly::ggplotly(vis_*(data))`. In the future they will
  be written in `plotly` so that they can be generated much faster

## Minor improvements

- corrected testing for `vis_*` family
- added .svg graphics for correct vdiffr testing
- improved hover print method for plotly.

# visdat 0.0.6.9000 (2017/02/26)

## New Features

- axes in `vis_` family are now flipped by default
- `vis_miss` now shows the % missingness in a column, can be disabled by setting `show_perc_col` argument to FALSE
- removed `flip` argument, as this should be the default 

## Minor Improvements

- added internal functions to improve extensibility and debugging - `vis_create_`, `vis_gather_` and `vis_extract_value_`.
- suppress unneeded warnings arising from compiling factors

# visdat 0.0.5.9000 (2017/01/09)

## Minor Improvements

- Added testing for visualisations with `vdiffr`. Code coverage is now at 99%
- Fixed up suggestions from `goodpractice::gp()`
- Submitted to rOpenSci onboarding
- `paper.md` written and submitted to JOSS

# visdat 0.0.4.9999 (2017/01/08)

## New Feature

- Added feature `flip = TRUE`, to `vis_dat` and `vis_miss`. This flips the x axis and the ordering of the rows. This more closely resembles a dataframe.
- `vis_miss_ly` is a new function that uses plotly to plot missing data, like `vis_miss`, but interactive, without the need to call `plotly::ggplotly` on it. It's fast, but at the moment it needs a bit of love on the legend front to maintain the style and features (clustering, etc) of current `vis_miss`.
- `vis_miss` now gains a `show_perc` argument, which displays the % of missing and complete data. This is switched on by default and addresses issue #19.

## New Feature (under development)

- `vis_compare` is a new function that allows you to compare two dataframes of the same dimension. It gives a fairly ugly warning if they are not of the same dimension.
- `vis_dat` gains a "palette" argument in line with [issue 26](https://github.com/ropensci/visdat/issues/26), drawn from http://colorbrewer2.org/, there are currently three arguments, "default", "qual", and "cb_safe". "default" provides the ggplot defaults, "qual" uses some colour blind **unfriendly** colours, and "cb_safe" provides some colours friendly for colour blindness.

## Minor Improvements

- All lines are < 80 characters long
- removed all instances of `1:rnow(x)` and replaced with `seq_along(nrow(x))`.
- Updated documentation, improved legend and colours for `vis_miss_ly`.
- removed export for `vis_dat_ly`, as it currently does not work.
- Removed a lot of unnecessary @importFrom tags, included magrittr in this, and added magrittr to Imports
- Changes ALL CAPS Headers in news to Title Case
- Made it clear that `vis_guess()` and `vis_compare` are very beta
- updated documentation in README and `vis_dat()`, `vis_miss()`, `vis_compare()`, and `vis_guess()`
- updated pkgdown docs
- updated DESCRIPTION URL and bug report
- Changed the default colours of `vis_compare` to be different to the ggplot2 standards.
- `vis_miss` legend labels are created using the internal function `miss_guide_label`. `miss_guide_label` will check if data is 100% missing or 100% present and display this in the figure. Additionally, if there is less than 0.1% missing data, "<0.1% missingness" will also be displayed. This sort of gets around issue #18 for the moment.
- tests have been added for the `miss_guide_label` legend labels function.
- Changed legend label for `vis_miss`, `vis_dat`, and `vis_guess`. 
- updated README
- Added vignette folder (but not vignettes added yet)
- Added appveyor-CI and travis-CI, addressing issues #22 and #23


## Bug Fixes

- Update `vis_dat()` to use `purrr::dmap(fingerprint)` instead of `mutate_each_()`. This solves issue #3 where `vis_dat` couldn't take variables with spaces in their name.

# visdat 0.0.3.9000
=========================

## New Features

- Interactivity with `plotly::ggplotly`! Funcions `vis_guess()`, `vis_dat()`, and `vis_miss` were updated so that you can make them all interactive using the latest dev version of `plotly` from Carson Sievert.


# visdat 0.0.2.9000
=========================

## New Features

- Introducing `vis_guess()`, a function that uses the unexported function `collectorGuess` from `readr`.


# visdat 0.0.1.9000
=========================

## New Features

- `vis_miss()` and `vis_dat` actually run
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
## Test environments
* local OS X install, R 4.1.0
* github actions testing for devel, release, and ubuntu, windows, and macOSX
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

There were no ERRORs or WARNINGs or NOTEs

## revdepcheck results

We checked 4 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 1 packages

Issues with CRAN packages are summarised below.

### Failed to check

* bsem (R CMD check timed out)
---
title: 'visdat: Visualising Whole Data Frames'
authors:
- affiliation: 1
  name: Nicholas Tierney
  orcid: 0000-0003-1460-8722
date: "01 August 2017"
output:
  html_document:
    keep_md: yes
  pdf_document: default
bibliography: paper.bib
tags:
- visualisation
- R
- exporatory data analysis
affiliations:
- index: 1
  name: Monash University
---


# Summary

When you receive a new dataset you need to look at the data to get a sense of what is in it, and understand potential problems and challenges to get it analysis-ready. "Taking a look at the data" can mean different things. For example: examining statistical summaries (minimum, maximum, mean, inter-quartile range), finding missing values, checking data formatting, creating graphical summaries such as histograms, scatter plots, box plots, and more.

When handling typical real-world data, these preliminary exploratory steps can become difficult to perform when values are not what you expect. For example, income might be a factor instead of numeric, date could be a number not a character string or a date class, or values could be missing when they shouldn't be. Often times, you discover that you had expectations of the data, which are hard to realise until they are a problem. This is similar to how one might not think to buy more light bulbs until one goes out: when you use data in an exploratory scatter plot, or a preliminary model, you often don't realise your data is in the wrong format until that moment. What is needed is a birds eye view of the data, which tells you what classes are in the dataframe, and where the missing data are.

`visdat` is an R [@Rcore] package that provides a tool to "get a look at the data" by creating heatmap-like visualisations of an entire dataframe, which provides information on: classes in the data, missing values, and also comparisons between two datasets. `visdat` takes inspiration from [`csv-fingerprint`](https://github.com/setosa/csv-fingerprint), and is powered by ggplot2 [@ggplot2], which provides a consistent, powerful framework for visualisations that can be extended if needed. 

Plots are presented in an intuitive way, reading top down, just like your data. Below is a plot using `vis_dat()` of some typical data containing missing values and data of a variety of classes.


```r
library(visdat)
vis_dat(typical_data)
```

![](paper_files/figure-html/load-data-1.png)<!-- -->

`visdat` will continue to be improved over time, to improve speed in computation and improve interactive plotting.

# Acknowlegements

I would like to thank the two reviewers, Mara Averick and Sean Hughes, for their helpful suggestions that resulted in a much better package, and rOpenSci for providing the support of the onboarding package review that facilitated these improvements.

# References
# Things to check with reviewers

# @seaaan

> In the "Using visdat" vignette, it says "missing data represented by black", but it shows up as gray on my computer.

- I think I've fixed this now, does it still happen?

> I don't know why, but in the plotly section, a vis_dat() (non-interactive) plot of df_test appears between the first two interactive plots. I can't explain it, hopefully it's just a weird quirk on my computer.

- I can't seem to replicate this, can you maybe take a screenshot of this if this is still happening?


> Did you mean to export guess_type()?

- Yes, as I thought that it might be a useful function for users, although I can't think of a good usecase outside of `vis_guess` right now. Perhaps it might be best for `guess_type` to be unexported?


> vis_compare and vis_guess are indicated as being in beta, which seems also to apply to the plotly versions of the functions. The message that they emit is a helpful indication that they may change in the future. Before submitting to CRAN, however, you might consider moving those functions to the development version of the package and only uploading the functions with a stable API, then adding the beta version later once they stabilize. This is a judgment call; I think I tend towards being conservative on this issue personally.

- Agreed, I need to create a dev branch, where these experimental features can reside.

*task*
- create "dev" branch on github where I put the experimental features of visdat.
# vis_binary sends an error when used with the wrong data

    Code
      vis_binary(iris)
    Error <simpleError>
      data input can only contain binary values - this means either 0 or 1, or NA. Please subset the data to be binary values, or see ?vis_value.

# vis_dat and vis_miss throw warnings when the DF is above size

    Code
      vis_dat(big_df)
    Error <simpleError>
      Data exceeds recommended size for visualisation, please consider
               downsampling your data, or set argument 'warn_large_data' to FALSE.

---

    Code
      vis_miss(big_df)
    Error <simpleError>
      Data exceeds recommended size for visualisation, please consider
               downsampling your data, or set argument 'warn_large_data' to FALSE.

# vis_dat fails when the wrong palette is provided

    Code
      vis_dat(typical_data, palette = "wat")
    Error <simpleError>
      palette arguments need to be either 'qual' 'cb_safe' or 'default'

# vis_dat fails when an object of the wrong class is provided

    Code
      vis_dat(AirPassengers)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

# vis_cor sends an error when used with the wrong data

    Code
      vis_cor(iris)
    Error <simpleError>
      data input can only contain numeric values, please subset the data to the numeric values you would like. dplyr::select_if(data, is.numeric) can be helpful here!

# vis_cor fails when an object of the wrong class is provided

    Code
      vis_cor(AirPassengers)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

# vis_guess fails when the wrong palette is provided

    Code
      vis_guess(typical_data, palette = "wat")
    Error <simpleError>
      palette arguments need to be either 'qual' 'cb_safe' or 'default'

# vis_guess fails when an object of the wrong class is provided

    Code
      vis_guess(AirPassengers)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

# vis_expect fails when an object of the wrong class is provided

    Code
      vis_expect(AirPassengers, ~.x < 20)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

# vis_value sends an error when used with the wrong data

    Code
      vis_value(iris)
    Error <simpleError>
      data input can only contain numeric values, please subset the data to the numeric values you would like. dplyr::select_if(data, is.numeric) can be helpful here!

# vis_miss fails when an object of the wrong class is provided

    Code
      vis_miss(AirPassengers)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

# vis_compare will not accept two dataframes of differing dims

    Code
      vis_compare(iris, iris_add)
    Error <simpleError>
      vis_compare requires identical dimensions of df1 and df2

# vis_compare fails when an object of the wrong class is provided

    Code
      vis_compare(iris, AirPassengers)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

---

    Code
      vis_compare(AirPassengers, iris)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

---

    Code
      vis_compare(AirPassengers, AirPassengers)
    Error <simpleError>
      vis_dat requires a data.frame but the object I see has class/es: ts

Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://stackoverflow.com/> or <https://community.rstudio.com/>.

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](https://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

---

Brief description of the problem

```r
# insert reprex here
```
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.1.0 (2021-05-18) |
|os       |macOS Big Sur 11.4           |
|system   |x86_64, darwin17.0           |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_AU.UTF-8                  |
|ctype    |en_AU.UTF-8                  |
|tz       |Australia/Perth              |
|date     |2021-07-05                   |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|visdat  |0.5.3 |0.6.0 |*  |
|cpp11   |NA    |0.3.1 |*  |

# Revdeps

## Failed to check (1)

|package |version |error |warning |note |
|:-------|:-------|:-----|:-------|:----|
|bsem    |1.0.0   |1     |        |     |

*Wow, no problems at all. :)*---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, echo = FALSE}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)

```

# visdat <img src="man/figures/visdat-logo.png" align="right" />

<!-- badges: start -->
[![rOpenSci Badge](https://badges.ropensci.org/87_status.svg)](https://github.com/ropensci/software-review/issues/87)[![JOSS status](https://joss.theoj.org/papers/10.21105/joss.00355/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00355)[![DOI](https://zenodo.org/badge/50553382.svg)](https://zenodo.org/badge/latestdoi/50553382)[![R-CMD-check](https://github.com/ropensci/visdat/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/visdat/actions)[![Codecov test coverage](https://codecov.io/gh/ropensci/visdat/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/visdat?branch=master)[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/visdat)](https://cran.r-project.org/package=visdat)[![CRAN Logs](http://cranlogs.r-pkg.org/badges/visdat)](https://CRAN.R-project.org/package=visdat)[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
<!-- badges: end -->

# How to install

visdat is available on CRAN

```{r install-cran, eval = FALSE}

install.packages("visdat")

```

If you would like to use the development version, install from github with:

```{r installation, eval = FALSE}

# install.packages("devtools")
devtools::install_github("ropensci/visdat")
```

# What does visdat do?

Initially inspired by
[`csv-fingerprint`](https://github.com/setosa/csv-fingerprint), `vis_dat` helps
you visualise a dataframe and "get a look at the data" by displaying the
variable classes in a dataframe as a plot with `vis_dat`, and getting a brief
look into missing data patterns using `vis_miss`.

`visdat` has 6 functions:

- `vis_dat()` visualises a dataframe showing you what the classes of the columns
are, and also displaying the missing data.

- `vis_miss()` visualises just the missing data, and allows for missingness to
be clustered and columns rearranged. `vis_miss()` is similar to
`missing.pattern.plot` from the
[`mi`](https://CRAN.R-project.org/package=mi) package.
Unfortunately `missing.pattern.plot` is no longer in the `mi` package (as of
14/02/2016).

- `vis_compare()` visualise differences between two dataframes of the same 
  dimensions

- `vis_expect()` visualise where certain conditions hold true in your data

- `vis_cor()` visualise the correlation of variables in a nice heatmap

- `vis_guess()` visualise the individual class of each value in your data
- `vis_value()` visualise the value class of each cell in your data
- `vis_binary()` visualise the occurrence of binary values in your data

You can read more about visdat in the vignette, ["using visdat"](http://visdat.njtierney.com/articles/using_visdat.html).

## Code of Conduct

Please note that the visdat project is released with a [Contributor Code of Conduct](https://github.com/ropensci/visdat/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

# Examples

## Using `vis_dat()`

Let's see what's inside the `airquality` dataset from base R, which contains
information about daily air quality measurements in New York from May to
September 1973. More information about the dataset can be found with
`?airquality`.

```{r vis-dat-aq}

library(visdat)

vis_dat(airquality)

```

The plot above tells us that R reads this dataset as having numeric and integer
values, with some missing data in `Ozone` and `Solar.R`. The classes are
represented on the legend, and missing data represented by grey. The
column/variable names are listed on the x axis.

## Using `vis_miss()`

We can explore the missing data further using `vis_miss()`:

```{r vis-miss-aq}

vis_miss(airquality)

```

Percentages of missing/complete in `vis_miss` are accurate to 1 decimal place.

You can cluster the missingness by setting `cluster = TRUE`:

```{r vis-miss-aq-cluster}

vis_miss(airquality, 
         cluster = TRUE)

```

Columns can also be arranged by columns with most missingness, by setting
`sort_miss = TRUE`:

```{r vis-miss-aq-sort-miss}

vis_miss(airquality,
         sort_miss = TRUE)

```

`vis_miss` indicates when there is a very small amount of missing data at <0.1%
missingness:

```{r vis-miss-test}

test_miss_df <- data.frame(x1 = 1:10000,
                           x2 = rep("A", 10000),
                           x3 = c(rep(1L, 9999), NA))

vis_miss(test_miss_df)

```

`vis_miss` will also indicate when there is no missing data at all:

```{r vis-miss-mtcars}

vis_miss(mtcars)

```

To further explore the missingness structure in a dataset, I recommend the
[`naniar`](https://github.com/njtierney/naniar) package, which provides more
general tools for graphical and numerical exploration of missing values.

## Using `vis_compare()`

Sometimes you want to see what has changed in your data. `vis_compare()` displays the differences in two dataframes of the same size. Let's look at an example.

Let's make some changes to the `chickwts`, and compare this new dataset:

```{r vis-compare-iris}
set.seed(2019-04-03-1105)
chickwts_diff <- chickwts
chickwts_diff[sample(1:nrow(chickwts), 30),sample(1:ncol(chickwts), 2)] <- NA

vis_compare(chickwts_diff, chickwts)

```

Here the differences are marked in blue.

If you try and compare differences when the dimensions are different, you get 
an ugly error:

```{r vis-compare-error, eval = FALSE}

chickwts_diff_2 <- chickwts
chickwts_diff_2$new_col <- chickwts_diff_2$weight*2

vis_compare(chickwts, chickwts_diff_2)
# Error in vis_compare(chickwts, chickwts_diff_2) : 
#   Dimensions of df1 and df2 are not the same. vis_compare requires dataframes of identical dimensions.
```


## Using `vis_expect()`

`vis_expect` visualises certain conditions or values in your data. For example,
If you are not sure whether to expect values greater than 25 in your data
(airquality), you could write: `vis_expect(airquality, ~.x>=25)`, and you can
see if there are times where the  values in your data are greater than or equal
to 25:

```{r vis-expect}

vis_expect(airquality, ~.x >= 25)

```

This shows the proportion of times that there are values greater than 25, as
well as the missings.

## Using `vis_cor()`

To make it easy to plot correlations of your data, use `vis_cor`:

```{r vis-cor}

vis_cor(airquality)

```

## Using `vis_value`

`vis_value()` visualises the values of your data on a 0 to 1 scale.

```{r vis-value}
vis_value(airquality)
```

It only works on numeric data, so you might get strange results if you are using   factors:

```{r iris-error, eval = FALSE}
library(ggplot2)
vis_value(iris)
```

```
data input can only contain numeric values, please subset the data to the numeric values you would like. dplyr::select_if(data, is.numeric) can be helpful here!
```

So you might need to subset the data beforehand like so:

```{r iris-error-fix}
library(dplyr)

iris %>%
  select_if(is.numeric) %>%
  vis_value()
```

## Using `vis_binary()`

`vis_binary()` visualises binary values. See below for use with example data, `dat_bin`

```{r vis-bin}
vis_binary(dat_bin)
```

If you don't have only binary values a warning will be shown.

```{r vis-bin-airq, eval = FALSE}
vis_binary(airquality)
```

```
Error in test_if_all_binary(data) : 
  data input can only contain binary values - this means either 0 or 1, or NA. Please subset the data to be binary values, or see ?vis_value.
```


## Using `vis_guess()`

`vis_guess()` takes a guess at what each cell is. It's best illustrated using
some messy data, which we'll make here:

```{r create-messy-vec}

messy_vector <- c(TRUE,
                  T,
                  "TRUE",
                  "T",
                  "01/01/01",
                  "01/01/2001",
                  NA,
                  NaN,
                  "NA",
                  "Na",
                  "na",
                  "10",
                  10,
                  "10.1",
                  10.1,
                  "abc",
                  "$%TG")

set.seed(2019-04-03-1106)
messy_df <- data.frame(var1 = messy_vector,
                       var2 = sample(messy_vector),
                       var3 = sample(messy_vector))

```


```{r vis-guess-messy-df, fig.show='hold', out.width='50%'}

vis_guess(messy_df)
vis_dat(messy_df)

```

So here we see that there are many different kinds of data in your dataframe. As
an analyst this might be a depressing finding. We can see this comparison above.

# Thank yous

Thank you to Ivan Hanigan who [first
commented](https://www.njtierney.com/post/2015/11/12/ggplot-missing-data/)
this suggestion after I made a blog post about an initial prototype
`ggplot_missing`, and Jenny Bryan, whose
[tweet](https://twitter.com/JennyBryan/status/679011378414268416) got me
thinking about `vis_dat`, and for her code contributions that removed a lot of
errors.

Thank you to Hadley Wickham for suggesting the use of the internals of `readr`
to make `vis_guess` work. Thank you to Miles McBain for his suggestions on how
to improve `vis_guess`. This resulted in making it at least 2-3 times faster.
Thanks to Carson Sievert for writing the code that combined `plotly` with
`visdat`, and for Noam Ross for suggesting this in the first place. Thank you
also to Earo Wang and Stuart Lee for their help in getting capturing expressions
in `vis_expect`.

Finally thank you to [rOpenSci](https://github.com/ropensci) and it's amazing
[onboarding process](https://github.com/ropensci/software-review), this process has
made visdat a much better package, thanks to the editor Noam Ross (@noamross),
and the reviewers Sean Hughes (@seaaan) and Mara Averick (@batpigandme).

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# Response to reviewers

I would like to thank the reviewers for taking the time to review visdat so thoroughly. It is fantastic to hear such encouraging comments from members of the community.

Below I respond to suggested changes and other constructive criticisms.

# @seaaan's Review

1. Documentation: README

> You might consider breaking the "experimental" parts of the package into a separate README like you do with the vignettes.

- Agreed! I would prefer to keep the README as a whole piece, Perhaps instead I could just link to the experimental vignette as in the README. This has been logged in [issue 39](https://github.com/njtierney/visdat/issues/39) and fixed in the indicated commit messages.

> Link to the vignettes from the README.

- Agreed, this has been addressed in [issue 39](https://github.com/njtierney/visdat/issues/39)

> The title of the github project is "A package to assist in visually testing a dataset", which doesn't describe the project that well. I suggest changing it to be the same title as on the DESCRIPTION file, which is effective in quickly summarizing what the package does.

- Agreed, this makes more sense - changed!

> You don't introduce the airquality data set at first usage, but you do introduce it after the second plot. Move the introduction to first use.

- Thank you for this, this has been addressed in [issue 39](https://github.com/njtierney/visdat/issues/39).

> I don't understand this sentence: "When there is <0.1% of missingness, vis_miss indicates that there is >1% missingness."

- Agreed, this has been changed in [issue 39](https://github.com/njtierney/visdat/issues/39).

2. Vignettes:

> Make sure you have permission to reproduce the image from R4DS. The first few paragraphs of the "Using visdat" vignette are a bit unnecessary in my opinion. You could start much more simply, with something like: "When you get a new data set, you need to get a sense of what it contains and potential problems with it." and then continue with the discussion of different ways to approach this challenge (e.g. head, glimpse, ...).

- I agree with you that this is a clearer way to introduce the topic - this is now [changed](https://github.com/njtierney/visdat/commit/3f74141d2e9de4bc0ae870fc874b17a21829618e).

> In the "Using visdat" vignette, it says "missing data represented by black", but it shows up as gray on my computer.

- This is fixed this now.

> I don't know why, but in the plotly section, a vis_dat() (non-interactive) plot of df_test appears between the first two interactive plots. I can't explain it, hopefully it's just a weird quirk on my computer.

- I can't seem to replicate this, can you maybe take a screenshot of this if this is still happening?

> There is a "future work and experimental features" section in the "Using visdat" vignette -- I suggest transferring that to the experimental vignette and just linking to it from the "Using visdat" vignette. Perhaps move the vis_guess() and the interactive plots to the experimental features vignette as well, since they seem to be evolving.

- Great point! Thank you! I have now linked to this additional experimental vignette in the "using visdat" vignette.

> The plot in the experimental features vignette shows up small, I think you need to add knitr::opts_chunk$set(fig.width = 5, fig.height = 4) to the first chunk like you have in the other vignette.

- Good idea! This has been addressed now.

# Function documentation

> Did you mean to export guess_type()?

- Yes, as I thought that it might be a useful function for users, although I can't think of a good usecase outside of `vis_guess` right now. For the moment I will unexport guess_type now.

> Formatting: links, code, and other formatting need to be done with .Rd syntax. For example, for code, use \code{}, not backticks. For bold, use \strong{} instead of asterisks.

- I have now converted this over to use markdown with roxygen, this should be fixed now.

> Return values: The documentation doesn't include the "Value" section which is typically used to say what a function returns. I would add this. E.g. @return A \code{ggplot} object displaying the type of values in the data frame and the position of any missing values.

- Thank you for this - this has been changed for `vis_dat.R`, `vis_compare.R`, `vis_guess.R`, and `vis_miss.R`.

> ?vis_compare: The documentation for this function needs to be updated. It is more a list of ideas for how to implement the function than a description of what it does.

- Thank you for this comment, this has [been updated.](https://github.com/njtierney/visdat/blob/master/R/vis_compare.R)

> ?vis_guess: First sentence of description seems to have an extra word or phrase, not sure exactly what was intended.

- Thank you for this comment, this has [been updated.](https://github.com/njtierney/visdat/blob/master/R/vis_guess.R)

> ?vis_miss_ly: The reference to vis_miss in the "See Also" section should be a link. Like this: \code{\link[visdat]{vis_miss}}

- Thank you for this comment, this has [been updated.](https://github.com/njtierney/visdat/blob/master/R/vis_miss_ly.R)

> ?visdat: It's good that this page exists. However, it doesn't have any content -- add a brief overview of the package with links to the documentation page for the main functions.

- Excellent point, I have addressed this in this [commit](https://github.com/njtierney/visdat/commit/106e24a4e186519f6ede615c3b2822073a03d780)

# Examples

> Don't need to call library for either visdat or packages you use internally (only if you actually call a function from another package in the example code itself).

- Great point, thanks for picking up on this, I have made the final changes [here](https://github.com/njtierney/visdat/commit/d71d2ac8c0708a327d51217f0298b4d3dec2b970)

> example(vis_dat): "palette" is misspelled.

- Thank for you picking up on this, this was fixed [here](https://github.com/njtierney/visdat/commit/31909e241c1e1cc1b640c55edaea9c228f92976e)

> example(vis_compare) is a good example but gives a bunch of warning messages.
> example(vis_guess): gives warning message

- These now both give warning messages like: 

```
vis_compare / (vis_guess) is in BETA! If you have suggestions or errors
post an issue at https://github.com/njtierney/visdat/issues
```

Do you think this is OK?

# Community guidelines

> There are no guidelines for how to contribute, either in the README or in a CONTRIBUTING file.

- This has now been added [here](https://github.com/njtierney/visdat/commit/c7ea5bb6bb12a43719fc4dbdfbccc22153072b62).

# JOSS

> The paper seems more detailed than most other JOSS papers and goes into specific functions. Additionally, does JOSS allow figures? I didn't see any in the papers I looked at or see mentions of them in the author guidelines. But I'm not an expert on JOSS, so maybe someone else can weigh in on that.

- I believe JOSS allows figures - the tidytext paper [has a figure](http://joss.theoj.org/papers/89fd1099620268fe0342ffdcdf66776f).

> The paper doesn't have any references in the "References" section, but does have inline references.

- These references are provided in the yaml, and in a .bib file - hopefully these get compiled with the JOSS paper perhaps @arfon can chime in here?

> If you keep it, make sure you have permission to reproduce the image from R4DS.

- I have removed the figure now, in favour of your earlier suggestion to keep things simpler - thank you for this suggestion! :)

# Functionality

## Common issues

> It is unintuitive to me to have the rows in reverse order (i.e. row 1 at the bottom) and for the columns to be clustered by type, rather than appear in the order they appear in the data frame. I think the default behavior should be for rows and columns to appear in the same order as in the input. Additionally, always putting the titles at the top of the columns makes sense to me. Including sort_type as an option might be useful, but by default it should be off. I'm not quite sure of the use case for the flip argument, but maybe there is one!

> The flip argument should be available for all of the functions if you provide it for any, for the sake of consistency. Possibly also the palette argument.

- I agree that the rows should be in the order that you see the dataframes. I would also prefer for the column names to be on top of the visualisation - I have now implemented this, although I think that the presentation of the column names is perhaps not as nice.

- The flip argument was placed there to experiment with this idea of flipping the plot - as you suggested, to have row 1 at the top. The flip argument has now been removed, as I cannot now really think of any times where I would prefer the old arrangement of rows.

- I think that the default sort_type option should be to put all of the similar columns together - to me this makes it easier to identify what is similar, and what is different; when I start looking at a dataset I usually want to know what variables are what class. I am willing to change this, however.^

- I agree with you re the palette argument being provided for all plots, although I'm not sure how this could work for vis_miss. But perhaps I should allow for users to insert their own palette for this. Personally I think that there should be a standard palette where the common classes (character, factor, integer, numeric, datetime, etc) each have their own distinct colour, and other classes are given different factors.

> A number of the functions emit warning messages ("Warning message: attributes are not identical across measure variables; they will be dropped"), specifically when they are called with a data frame that has >1 factor column with different levels. This arises from the call tidyr::gather_(x, "variables", "value", names(x))$value. This message is not relevant to the end user and should be suppressed.

> There is some code duplication across some of the functions, where you could replace it with a helper function. The code creating the plots is generally similar across functions. In addition, the following code blocks appear in essentially identical form in >1 function:

```{r eval = FALSE}

  dplyr::mutate(rows = seq_len(nrow(.))) %>%
  tidyr::gather_(key_col = "variables",
                 value_col = "valueType",
                 gather_cols = names(.)[-length(.)])
                 
```

> Using a helper function for this next one in particular would allow you to suppress the warning messages just once:

```
  d$value <- tidyr::gather_(x, "variables", "value", names(x))$value
```

> The code for flipping the rows is also duplicated:

```
 if (flip == TRUE){
    suppressMessages({
      vis_dat_plot +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_x_discrete(position = "top",
                                  limits = type_order_index) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(x = "")
    })
  } else if (flip == FALSE){
    vis_dat_plot
  }
  
```

- Thank you for this, I have created a helper function for these code duplications, it makes things much easier to handle! These are now in the [internals file](https://github.com/njtierney/visdat/blob/fix-functionality/R/internals.R)

> vis_compare and vis_guess are indicated as being in beta, which seems also to apply to the plotly versions of the functions. The message that they emit is a helpful indication that they may change in the future. Before submitting to CRAN, however, you might consider moving those functions to the development version of the package and only uploading the functions with a stable API, then adding the beta version later once they stabilize. This is a judgment call; I think I tend towards being conservative on this issue personally.

- Agreed, I will create a dev branch, where these experimental features can reside.

> This is a minor issue, but was so thoroughly drummed into me by CS professors that I have to mention it. In calls like if (sort_type == TRUE), the == TRUE is redundant, so it's equivalent to just write if (sort_type). Similarly, if (x == FALSE) is equivalent to if (!x). There are a number of these in the code for various functions.

- I've actually not heard this before, but I can certainly understand why you would want to avoid called such as if(x == TRUE) over if(x). This has now been changed.

# vis_dat()

> Maximum data frame size: the heatmaps stop working on my computer with very large data sets. E.g. library(nycflights13); vis_dat(flights) hangs for a minute or so and then displays an empty grid. It works with 10000 rows, but not 100000. Similar issue with many columns. Especially because this package is designed for people to get a sense of new data sets, where users may not have realized how big the data are, I think it makes sense to prevent this issue. Possibilities: (1) max_rows = 10000 parameter, which stops if nrow(data) > max_rows and gives a descriptive error message. Then the user can increase the maximum or subset their data frame. Might also need to consider max_cols because 10000 rows and 2 cols is very different from 10000 rows and 500 cols. (2) Downsample the non-missing rows somehow and indicate to the user where omissions were made. (3) ???

- This is an interesting one! It was brought to my attention [here](https://github.com/njtierney/visdat/issues/32), by @MilesMcBain, but I couldn't replicate this issue. Similarly, I can actually display the nycflights data. So it seems that this is to do with the fact that computations power has a bit to do with it. I think that option 1 is the best case here, but option 2 actually gets a little bit difficult, because this leads to another issue of effectively downsampling a dataset whilst maintaining data structure. This could be a problem, say, if I did perform down-sampling and then the user looked at the data with visdat and plotly, they might get a mis-represented sample. Likewise, if there is a very large proportion of missing data, this might be difficult to preserve exactly. Ways around this might include some very explicit warning messages, or perhaps even a note directly on the plot that the data have been downsampled. ^At the moment one solution forward is to have a check first, something like:

```r
vis_dat(really_large_data)

Warning: This dataset contains over 10,000 rows, which may take a while to render. Are you sure that you want to continue?

1: Yes
2. No (You should use a sample of your data)
```

What do you think?

> The palette argument doesn't cause color changes, as far as I can tell. I think you need to do vis_dat_plot <- vis_dat_plot + <palette>, whereas currently the value of vis_dat_plot isn't being updated in the blocks related to color palettes. Additionally, you could do ggplot2::scale_color_brewer(type = "qual", palette = "Set1") rather than explicitly writing out the color names, although the colors would be in a slightly different order than you have now.

- ^I believe that this is now fixed, 

```r
vis_dat(airquality, palette = "qual")
vis_dat(airquality, palette = "cb_safe")
```

Return different plots now, right?

> If the user provides a palette name, you should check that it's valid. If not, throw an error with a message. Currently, if you misspell a palette name, it just goes ahead with the default palette.

- Great point! I have fixed this now with:

```r
     } else  {
       warning("palette arguments need to be either 'qual' 'cb_safe' or 'default'")
     } # close else brace

```

at the end of `vis_dat.R`

# Minor comments

> There are some lines of commented-out code: # x = airquality and # mutate_each_(funs(fingerprint), tbl_vars(.)) %>%

> Inside of the first if-block, there are some comments about arranging columns in order of missingness, but as far as I can tell, that is not done.

- Thank you for picking up on these comments, they are now removed!

## vis_miss()

> It would be cool to indicate the % missing for each column individually (maybe in column labels?) like how you indicate the % missing for the data overall.

- great suggestion, this has been added! And can be turned off in the options with `show_perc_col == FALSE`.

## vis_compare()

> This is a neat function! It would be cool if you could generalize it to also handle data frames of unequal dimensions, perhaps by showing which rows/columns are only in one data frame or the other. Just an idea I had though, not a requirement at all!

- I am glad you like this! It is a difficult problem that I cannot work out, re the rows/cols matching when they are out of the dimension. I'm not quite sure what to do here to make this possible, I think that this should be tackled in a future release of `visdat`, and will be moved to a development branch.

> Warning message: "1: In if (dim(df1) != dim(df2)) { : the condition has length > 1 and only the first element will be used". This is a consequence of calling dim(df1) != dim(df2). The dim function returns a two element vector in this case, but the != only compares the first element of each vector. (I.e. it's only comparing the number of rows). Instead, use !identical(dim(df1), dim(df2)).

> In the error message about data frames of unequal dimensions, you mean to say "vis_compare only handles dataframes of identical dimensions", not "does not handle".

```r

  if (!identical(dim(df1), dim(df2))){
    stop("Dimensions of df1 and df2 are not the same. vis_compare currently only handles dataframes of identical dimensions.")
  }

```

- Thank you for both of these comments! These have been changed in `vis_compare.R`.

## vis_guess

> vis_guess lacks some of the parameters that vis_dat has, but should have the same options (except maybe not sort_type?)

- Agree, this is a good idea, I have added a helper palette function now so that `vis_guess` and `vis_dat` can share the same palette

> Performance: Setting the locale once and then calling collectorGuess instead of guess_parser is about 33% faster on my computer. I think the savings is by avoiding repeatedly calling the locale function and avoiding repeatedly calling stopifnot (inside the readr code).

```
l <- readr::locale()

output[!nas] <- vapply(FUN = function(y) readr:::collectorGuess(y, locale_ = l),
                       X = x[!nas], FUN.VALUE = character(1))
                       
```

> If you're interested in converting to C++ for speed, I played around a little bit with implementing a vectorized version for guess_parser in readr. It's about 20X faster. Basically, I wrote a wrapper around collectorGuess that operates on a vector instead of single elements. This is faster because it does the looping in C++ rather than repeatedly transitioning from R to C++ and because it does some initial code once per vector instead of once per element. It is currently written as a modification of readr because that was more convenient, but I imagine it is possible to implement in your package as well. Let me know if you want to discuss that, but I'm definitely not an expert. To get it up and running, clone readr and then put this in readr/src/collectorGuess.cpp:

```
std::string collectorGuess2(const std::string& input, LocaleInfo locale) {

  // Work from strictest to most flexible
  if (canParse(input, isLogical, &locale))
    return "logical";
  if (canParse(input, isInteger, &locale))
    return "integer";
  if (canParse(input, isDouble, &locale))
    return "double";
  if (canParse(input, isNumber, &locale))
    return "number";
  if (canParse(input, isTime, &locale))
    return "time";
  if (canParse(input, isDate, &locale))
    return "date";
  if (canParse(input, isDateTime, &locale))
    return "datetime";

  // Otherwise can always parse as a character
  return "character";
}

// [[Rcpp::export]]
CharacterVector collectorGuessVector(CharacterVector input, List locale_) {
  LocaleInfo locale(locale_);
  CharacterVector output(input.size());

  if (input.size() == 0 || allMissing(input)) {
    CharacterVector result(1);
    result[0] = "character";
    return result;
  }

  for (int i = 0; i < input.size(); i++) {
    output[i] = collectorGuess2(std::string(input[i]), locale);
  }

  return output;
}

```

> Then put this in readr/R/collectors.R:

```
#' @rdname parse_guess
#' @export
guess_parser_vector <- function(x, locale = default_locale()) {
  stopifnot(is.locale(locale))
  collectorGuessVector(x, locale)
}
```

> Then rebuild readr. Back in visdat, call guess_parser_vector once with the whole vector of unknowns instead of calling guess_parser repeatedly.

- ^OK, so I think that this is a fantastic idea! But I'm not sure on the best way to move this into production of `visdat`, since it sounds like we would need to borrow a chunk of code from readr and then modify it. I wonder if perhaps the solution here is to do a pull request for this specific feature in `readr`?

## vis_miss_ly

> The pop-up window is great! Can you make it the same as for ggplotly(vis_dat())? It's especially helpful to show the row number.

- I can get it close! Getting these the same will 

> Why not have the same arguments for this function as for vis_miss?

- At the moment this is because I do now know the API for plotly that well. Once we are happy with the default arguments and appearance for vis_miss, we can focus on getting the same appearance with plotly

## vis_dat_ly()

> With vis_dat(airquality) %>% ggplotly(), the window that appears on mouseover should say "type" instead of "valueType", "variable" instead of "variables" and "row" instead of "rows"

- I can get this to say "variable" instead of "variables", but I cannot change "valueType", as I think that changing this makes breaking changes. "rows" appears to be autogenerated from the `plotly::ggplotly` code, so I cannot change this.

> If you provide vis_miss_ly, it makes sense to me to also provide vis_dat_ly, and possibly also vis_compare_ly and vis_guess_ly. An initial implementation could just be vis_dat(x) %>% ggplotly(), but then you have the API in place at least.

- I have whipped up `vis_dat_ly`, `vis_compare_ly`, and `vis_guess_ly` using `ggplotly(vis_*)`.

- I would prefer to build it in straight plotly, as plotly is really fast, and the process of converting ggplot plots into plotly is a bit computationally intense. I will see what I can whip up in plotly in the future. I do not know much about it at the moment, but @carson has a free book, so I'll get onto that.

# Tests

> Cool use of vdiffr and generally good coverage. The tests don't all pass on my machine (they are mostly skipped), which I think is because I don't have the .svg reference files on my computer.

- Yes, you require the references.

> I suggest adding a test for each function using typical_data, because that currently emits a warning message. That way once you stop the error message from occurring, you'll have a regression test for the future.

- This is a great suggestion!

> Maybe I don't understand how vdiffr works, but I'm confused about test-vis-dat.R where in the calls to vdiffr you always provide vis_dat_plot, never any of the other plots you created. Additionally, both of the plots you create using non-default palettes are given the same name, so you're overwriting one of them.

- Oops, error on my part, thanks for picking that up! :)

> On my computer, the palette argument doesn't change the appearance of the vis_dat output, but the tests seem to pass, so can you check that you have a test that catches that issue? When I run vdiffr, vis-dat-qualitative-palette.svg and vis-dat-vanilla.svg appear to be the same.

- Yes, done. Thank you for picking up on this!

# @batpigandme's review

> vis_dat's power lies in the fact that it's essentially a visual analogue to one's data, if imagined as a frame/spreadsheet. In order to enhance this, I would suggest that the charts, by default, "read" in the spirit of this analogy by going from top to bottom.

- Absolutely agree, I have made these changes now and it makes a LOT more sense to do it this way!

## Documentation: README

> I have made some minor changes to the wording in README and submitted them as a pull request here, most of which I think are self-explanatory. 

- Thank you! Pull request was merged [here](https://github.com/njtierney/visdat/commit/868140da16fbfc9150d8359731f4e9eb3e817675)

> vis_dat(airquality) is numeric and integer values, when I run it. Thus, I've swapped numeric in for character as it was described.

- Thank you!

> Added canonical link to the mi package, since you do so for wakefield

- Thank you!

> Though the README suggests that when less than one percent is missing, "vis_miss indicates that there is >1%", the chart that results from vis_miss(test_miss_df) reads Missing (<0.1%) in the legend

- Agreed, this has now been fixed.

> When running the first example from README, (vis_dat(airquality)), the user will get a warning re. deprecation of dmap(). If this is not of consequence to the package's functionality, then you might consider adding it to the suppressed messages. If it's something most users supress locally, then you might suggest that somewhere in the README.

```
#> dmap() is deprecated. Please use the new colwise family in dplyr.
#> E.g., summarise_all(), mutate_all(), etc.
```

- Thank you for finding this, I have fixed this by suppressing the error as @seaaan suggested.

> Minor suggestion (which may or may not work comply with the actual vignette formatting) the purpose of demonstrating the difference between vis_guess() and vis_dat, you might consider showing the two side by side (perhaps have the vis_guess() example alone and then side by side). For this you could just add fig.show='hold', out.width='50%' to the chunk options in the Rmd.

- Excellent suggestion, done!

> Based on the README, the status of vis_dat_ly isn't totally clear (since it's not yet a function of its own).

- I have removed the mentions of these, and put them into the experimental features vignette, which is referenced in the README.

> I didn't quite understand the Visualising expectations section. That said, I've also never used expectation in assertr— so, this may fall under the category of things that those who need to know do, and, thus, no more need be said.

- I've constructed a small example to illustrate this at the end of the experimental features vignette.

## Documentation: Vignettes

> Some of the same copy edits I made in README I would also suggest in the vignettes (e.g. the compare this to...Where sentence that sort of goes around the vis_dat figure in the vis_guess example section)

> The experimental features vignette is a bit hard to follow at the moment, which may be in part due to the layout (the warning messages take up a good bit of space). The figure is also quite small. I'd recommend adjusting the settings to match your main vignette to make it more legible.

- Agreed, the warning messages have been fixed, and the plot size changed as well.

> vis_dat_ly(), and vis_miss_ly() aren't in the experimental vignette, though the opening paragraph suggests that they will be

- Thank you for picking that up, they have been added now.

## Documentation: Help

> The vis_guess description is a bit hard to follow. I'd start by stating what it does do, before making the comparison to vis_dat, just so it's a bit more self-contained.

- Done, great suggestion!

> The visdat and visdat-package sections are empty, when, I assume, you'll want the description to be included (right now description just says "visdat").

- Thank you for picking up on this, I've updated the main description file now

## Examples

> Note: the function examples in R all throw warning: dmap() is deprecated.
As noted in @seaaan's comments, there's also the warning re. dataframes of different dimensions
e.g. example("vis_compare") results in

```
  #>  Warning messages:
  #> 1: In if (dim(df1) != dim(df2)) { :
  #>  the condition has length > 1 and only the first element will be used
  #> 2: attributes are not identical across measure variables; they will be dropped
  #> 3: attributes are not identical across measure variables; they will be dropped
  
```

- Agreed, this error message for vis_compare has been addressed.

> vis_compare example also seems to cluster by default (all different data is shown at the bottom, of the two columns)

- I believe this is because I changed those particular values, I'll change the values to be something else now, so that it is clear that `vis_compare()` identifies values...

## Community guidelines

> I don't believe these are currently in the package.

- Thank you for picking up on this, I've added them now, :)

## Paper

> The vis_miss chunk is unnamed. I'd put a name in for you, but I can't knit the doc (see below). Either way, quick fix.

- Done!

> I'm getting an error 83 because pandoc-citeproc is missing. I believe this is something I'm missing on my machine, and, thus you don't need to worry about. But, just in case, I'm letting you know!

- This could be due to the yaml info, does it still happen?

## Functionality

### Installation and Tests

> Installation built successfully (on three different computers, and versions of R at that), as did automated testing. Good work!

- Hooray!

### Functionality

> All of the functions run as expected, and as described. My feedback here primarily concerns messages and warnings, that apply to several functions. 

> Warning messages can be scary, so it's always a relief see them acknowledged or shown in vignette and/or README examples (as is done re. non-identical attributes). vis_compare (a function for which I can imagine plenty of use cases) is currently returning enough text that it may leave the user confused as to what's most important, or an error on their end, etc. Currently, running vis_compare(iris_diff, iris) (as in the given example), returns:

```
#> vis_compare is still in BETA! If you have suggestions or errors,
#>        post an issue at https://github.com/njtierney/visdat/issues
#>  dmap() is deprecated. Please use the new colwise family in dplyr.
#>  E.g., summarise_all(), mutate_all(), etc.
#>  Warning messages:
#>  1: In if (dim(df1) != dim(df2)) { :
#>    the condition has length > 1 and only the first element will be used
#>  2: attributes are not identical across measure variables; they will be dropped
#>  3: attributes are not identical across measure variables; they will be dropped
```

> The function runs fine, and appears to have worked, but it's enough "warning" text, that it may warrant concern for the user.

- Absolutely agreed, this is a problem and barrier, the messages for using `vis_compare` have been reduced to: 

```r
vis_compare(typical_data)
...
```

> vis_guess is where the cell-to-row visual disparity becomes most confusing. I think this is a fantastic function, especially for beginners, and so the ability to essentially "read" the visual, as one would a table (top to bottom), would be invaluable.

- Agreed, this has been addressed by making the "flip" argument become the hard set default. I cannot honestly see a good reason to use the old print method, so I have removed it to make it read top down first.

> vis_dat with plotly / interactive functions are great, and, again, here it would be even more helpful if one could "read" the charts as one would a table.

- This has been implemented for all but `vis_miss_ly` - as this is written in plotly, and I'm not certain how to get it to flip.

# njtierney summary response and remaining tasks to complete

This has been a fantastic process to be a part of. To summarise, the main changes that I made were:

- re-orienting the vis_* family to read in the way that one reads a dataframe
- removing many of the error messages that occur due to some tidy/dplyr collation of factors
- refactoring repeated code into internal functions to more clearly express what was happening, and avoiding repetition. This resulted in the creation of internal functions `vis_gather_`, `vis_extract_value`, and `vis_create_`
- adding percentage missing to the column names in `vis_miss`

- For the moment I have decided to submit visdat to CRAN, with the main functions: 
  - vis_dat
  - vis_miss
  
Then, in a separate branch, there will be all of the features:
  - vis_compare()
  - vis_guess()
  - vis_*_ly() family
  
There is still more work to be done on these functions. Specifically:

- vis_guess. This needs to be made much faster, and will incorporate the C code from @seaaan to do so.
- vis_*_ly() family will be written in pure plotly, for speed.

I now summarise some of the changes that I wanted to confirm with the reviewers, @seaaan and @batpigandme:

1. Is the warning message for vis_guess/vis_compare code OK? `example(vis_guess) / example(vis_compare)`

```
vis_compare / (vis_guess) is in BETA! If you have suggestions or errors
post an issue at https://github.com/njtierney/visdat/issues
```

2. The default sort option. I think that the default sort_type option should be to put all of the similar columns together - to me this makes it easier to identify what is similar, and what is different; when I start looking at a dataset I usually want to know what variables are what class. I am willing to change this, however.

3. Maximum row size issue. This is an interesting one! It was brought to my attention [here](https://github.com/njtierney/visdat/issues/32), by @MilesMcBain, but I couldn't replicate this issue. Similarly, I can actually display the nycflights data. So it seems that this is to do with the fact that computations power has a bit to do with it. I think that option 1 is the best case here, but option 2 actually gets a little bit difficult, because this leads to another issue of effectively downsampling a dataset whilst maintaining data structure. This could be a problem, say, if I did perform down-sampling and then the user looked at the data with visdat and plotly, they might get a mis-represented sample. Likewise, if there is a very large proportion of missing data, this might be difficult to preserve exactly. Ways around this might include some very explicit warning messages, or perhaps even a note directly on the plot that the data have been downsampled. ^At the moment one solution forward is to have a check first, something like:

```r
vis_dat(really_large_data)

Warning: This dataset contains over 10,000 rows, which may take a while to render. Are you sure that you want to continue?

1: Yes
2. No (You should use a sample of your data)
```

4. Speed for vis_guess. @seaaan wrote some nice code to make this work much much faster, and I think that this is a fantastic idea! But I'm not sure on the best way to move this into production of `visdat`, since it sounds like we would need to borrow a chunk of code from readr and then modify it. I wonder if perhaps the solution here is to do a pull request for this specific feature in `readr`?

5. Can the two issues below be replicated?

> I don't know why, but in the plotly section, a vis_dat() (non-interactive) plot of df_test appears between the first two interactive plots. I can't explain it, hopefully it's just a weird quirk on my computer.

> I'm getting an error 83 because pandoc-citeproc is missing. I believe this is something I'm missing on my machine, and, thus you don't need to worry about. But, just in case, I'm letting you know!
---
title: 'visdat: Visualising Whole Data Frames'
authors:
- affiliation: 1
  name: Nicholas Tierney
  orcid: 0000-0003-1460-8722
date: "01 August 2017"
output:
  html_document:
    keep_md: yes
  pdf_document: default
bibliography: paper.bib
tags:
- visualisation
- R
- exporatory data analysis
affiliations:
- index: 1
  name: Monash University
---

# Summary

When you receive a new dataset you need to look at the data to get a sense of what is in it, and understand potential problems and challenges to get it analysis-ready. "Taking a look at the data" can mean different things. For example: examining statistical summaries (minimum, maximum, mean, inter-quartile range), finding missing values, checking data formatting, creating graphical summaries such as histograms, scatter plots, box plots, and more.

When handling typical real-world data, these preliminary exploratory steps can become difficult to perform when values are not what you expect. For example, income might be a factor instead of numeric, date could be a number not a character string or a date class, or values could be missing when they shouldn't be. Often times, you discover that you had expectations of the data, which are hard to realise until they are a problem. This is similar to how one might not think to buy more light bulbs until one goes out: when you use data in an exploratory scatter plot, or a preliminary model, you often don't realise your data is in the wrong format until that moment. What is needed is a birds eye view of the data, which tells you what classes are in the dataframe, and where the missing data are.

`visdat` is an R [@Rcore] package that provides a tool to "get a look at the data" by creating heatmap-like visualisations of an entire dataframe, which provides information on: classes in the data, missing values, and also comparisons between two datasets. `visdat` takes inspiration from [`csv-fingerprint`](https://github.com/setosa/csv-fingerprint), and is powered by ggplot2 [@ggplot2], which provides a consistent, powerful framework for visualisations that can be extended if needed. 

Plots are presented in an intuitive way, reading top down, just like your data. Below is a plot using `vis_dat()` of some typical data containing missing values and data of a variety of classes.

```{r load-data}
library(visdat)
vis_dat(typical_data)

```

`visdat` will continue to be improved over time, to improve speed in computation and improve interactive plotting.

# Acknowlegements

I would like to thank the two reviewers, Mara Averick and Sean Hughes, for their helpful suggestions that resulted in a much better package, and rOpenSci for providing the support of the onboarding package review that facilitated these improvements.

# References
---
title: "Customising colour palettes in visdat"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{customising-colour-palettes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(visdat)
```

# How to provide your own colour palette?

This vignette shoes you how to provide your own colour palette with `visdat`.

A `visdat` plot is a `ggplot` object - so we can use the tools from ggplot to 
tinker with colours. In this case, that is the `scale_fill_manual` function.

A "standard" visdat plot might be like so:

```{r standard}
vis_dat(typical_data)
```

You can name the colours yourself like so (after first loading the `ggplot` package.

```{r custom}
library(ggplot2)
vis_dat(typical_data) +
  scale_fill_manual(
    values = c(
      "character" = "red",
      "factor" = "blue",
      "logical" = "green",
      "numeric" = "purple",
      "NA" = "gray"
  ))
```

This is a pretty, uh, "popping" set of colours? You can also use some hex colours instead.

Say, taken from `palette()`:

```{r show-pal}
palette()
```


```{r pal-hex-visdat}
vis_dat(typical_data) +
  scale_fill_manual(
    values = c(
      "character" = "#61D04F",
      "factor" = "#2297E6",
      "logical" = "#28E2E5",
      "numeric" = "#CD0BBC",
      "NA" = "#F5C710"
  ))
```


How can we get nicer ones?

Well, you can use any of `ggplot`'s `scale_fill_*` functions from inside ggplot2

For example:

```{r scale-fill-brewer}
vis_dat(typical_data) +
  scale_fill_brewer()
```

```{r scale-fill-viridis}
vis_dat(typical_data) +
  scale_fill_viridis_d()
```

Happy colour palette exploring! You might want to take a look at some of the following colour palettes from other packages:

- [scico](https://github.com/thomasp85/scico#ggplot2-support)
- [colorspace](https://cran.r-project.org/web/packages/colorspace/vignettes/colorspace.html#Usage_with_ggplot2)
- [wesanderson](https://github.com/karthik/wesanderson#palettes)
---
title: "Using visdat"
author: "Nicholas Tierney"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using visdat}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, echo = FALSE, include = FALSE}

knitr::opts_chunk$set(fig.width = 5,
                      fig.height = 4)

```

When you get a new data set, you need to look at the data to get a sense of what it contains and potential problems with it. That's a key phrase here "looking at the data" - what does that mean?

On the one hand, you can look at the head of the data:

```{r head-iris}

head(iris)

```

Or you can have a `glimpse` at it through `dplyr::glimpse`

```{r glimpse}
library(dplyr)
glimpse(iris)

```

Here we see that we have doubles, and a factor. We get some insight into the data.

But we don't always have data like the canonical iris dataset. let's take a look at some data that might be a bit more typical of "messy" data using the `typical_data` dataset
from the `visdat` package.

```{r visdat-glimpse}
library(visdat)

glimpse(typical_data)

```

Looking at this, you might then ask:

> Isn't it odd that Income is a factor? And Age is a character? 

And you might start to wonder what else is different, what else changed? 

And it might be a bit unclear where to go from there. Do you plot the data? Why does my plot look weird? What are these other strange features in the data? The `visdat` package provides visualisations of an entire dataframe at once. Initially inspired by [`csv-fingerprint`](https://github.com/setosa/csv-fingerprint), `visdat` provides tools to create heatmap-like visualisations of an entire dataframe. `visdat` provides 2 main functions: `vis_dat` and `vis_miss`.

`vis_dat()` helps explore the data class structure and missingness:

```{r load-data}

vis_dat(typical_data)

```

And the `vis_miss` function provides a custom plot for missing data.

```{r example-vis-miss}

vis_miss(typical_data)

```

The name `visdat` was chosen as it borrows from the idea of [`testdat`](https://github.com/karthik/testdat), which provides unit testing for your data.  In a similar way, `visdat` provides visual tests, the idea being that first you visualise your data (`visdat`), then you run tests from `testdat`, or a package like `assertr`, to fix these errors.

## `vis_dat`

Let's see what's inside the dataset `airquality`, which contains information about daily air quality measurements in New York from May to September 1973. More information about the dataset can be found with `?airquality`.

```{r vis_dat}

vis_dat(airquality)

```
The plot above tells us that R reads this dataset as having numeric and integer values, with some missing data in `Ozone` and `Solar.R`. The classes are represented on the legend, and missing data represented by grey. The column/variable names are listed on the x axis. 

By default, `vis_dat` sorts the columns according to the type of the data in the vectors. You can turn this off by setting `sort_type = FALSE`. This feature is better illustrated using the `typical_data` dataset, created using [wakefield](https://github.com/trinker/wakefield) and contained within `visdat`.

```{r visdat-typical}

vis_dat(typical_data)

vis_dat(typical_data, 
        sort_type = FALSE)

```

## `vis_miss`

We can explore the missing data further using `vis_miss`.

```{r vis_miss}

vis_miss(airquality)

```

Notice that the percentages of missingness are provided in the data. These are accurate to 1 decimal place. `vis_miss` indicates when there is a very small amount of missing data at <0.1% missingness.


```{r vismiss-new-data}

df_test <- data.frame(x1 = 1:10000,
                      x2 = rep("A", 10000),
                      x3 = c(rep(1L, 9999), NA))

vis_miss(df_test)

```

`vis_miss` will also indicate when there is no missing data at all. 

```{r vismiss-mtcars}

df_test <- data.frame(x1 = 1:10000,
                      x2 = rep("tidy", 10000),
                      x3 = rep("data", 10000))

vis_miss(df_test)

```

Columns can be arranged by columns with most missingness, by setting `sort_miss = TRUE`.

```{r vismiss}

vis_miss(airquality,
         sort_miss = TRUE)

```

And missingness can be clustered by setting `cluster = TRUE`

```{r vis_miss-cluster}

vis_miss(airquality, 
         cluster = TRUE)

```

To further explore the missingness structure in a dataset, I recommend the [`naniar`](https://github.com/njtierney/naniar) package, which provides more general tools for graphical and numerical exploration of missing values.


## `vis_compare`

Sometimes you want to see what has changed in your data. `vis_compare()` displays the differences in two dataframes of the same size. Let's look at an example.

Let's make some changes to the `chickwts`, and compare this new dataset.

```{r vis-compare-iris}
set.seed(2019-04-03-1107)
chickwts_diff <- chickwts
chickwts_diff[sample(1:nrow(chickwts), 30),sample(1:ncol(chickwts), 2)] <- NA

vis_compare(chickwts_diff, chickwts)

```

Here the differences are marked in blue.

If you try and compare differences when the dimensions are different, you get an ugly error.

```{r vis-compare-error, eval = FALSE}

chickwts_diff_2 <- chickwts
chickwts_diff_2$new_col <- chickwts_diff_2$weight*2

vis_compare(chickwts, chickwts_diff_2)
# Error in vis_compare(chickwts, chickwts_diff_2) : 
#   Dimensions of df1 and df2 are not the same. vis_compare requires dataframes of identical dimensions.
```


## `vis_expect`

`vis_expect` visualises certain conditions or values in your data. For example,
If you are not sure whether to expect values greater than 25 in your data
(airquality), you could write: `vis_expect(airquality, ~.x >= 25)`, and you can
see if there are times where the  values in your data are greater than or equal
to 25.

```{r vis-expect}

vis_expect(airquality, ~.x >= 25)

```

This shows the proportion of times that there are values greater than 25, as well as the missings.

You could also, for example, explore a 
set of bad strings, or possible NA values and visualise where they are 
using `vis_expect(data, ~.x %in% bad_strings)` where `bad_strings` is a
character vector containing bad strings  like `N A`, `N/A` etc.

```{r vis-expect-bad-strings}

bad_data <- data.frame(x = c(rnorm(100), rep("N/A", 10)),
                       y = c(rep("N A ", 30), rnorm(80)))

vis_expect(bad_data, ~.x %in% c("N/A", "N A "))
```

## `vis_cor`

To make it easy to plot correlations of your data, use `vis_cor`:

```{r vis-cor}

vis_cor(airquality)

```

Under the hood, `vis_cor` is powered by the `cor` function in base R, and takes
a character string indicating which correlation coefficient (or covariance) is
to be computed. One of "pearson" (default), "kendall", or "spearman".

```{r vis-cor-spearman}

vis_cor(airquality, cor_method = "spearman")

```

You can also specify what to do for the missing data using the `na_action`
function, which again borrows from the `cor` methods. This can be "everything",
"all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs"
(default), e.g.:

```{r vis-cor-na-action}

vis_cor(airquality,
        na_action = "complete.obs")

```


## `vis_value`

`vis_value()` visualises the values of your data on a 0 to 1 scale.

```{r vis-value}
vis_value(airquality)
```

It only works on numeric data:

```{r diamonds-error, eval = FALSE}
vis_value(iris)
```


```
data input can only contain numeric values, please subset the data to the numeric values you would like. dplyr::select_if(data, is.numeric) can be helpful here!
``` 

So you might need to subset the data beforehand like so:

```{r diamonds-error-subset}

iris %>%
  select_if(is.numeric) %>%
  vis_value()
```

It can be useful to arrange your data before using `vis_value` to explore possible relationships in the data:

```{r airquality-arrange}
airquality %>%
  arrange(Wind) %>%
  vis_value()
```

## `vis_binary`

`vis_binary()` visualises the occurrence of binary values in your data. It is 
similar to `vis_value()` except it just focusses on values that are NA, 0, and 1.

```{r vis-binary}
vis_binary(dat_bin)
```



## `vis_guess`

`vis_guess()` takes a guess at what each cell is. It's best illustrated using some messy data, which we'll make here.

```{r create-messy-vec}

messy_vector <- c(TRUE,
                  T,
                  "TRUE",
                  "T",
                  "01/01/01",
                  "01/01/2001",
                  NA,
                  NaN,
                  "NA",
                  "Na",
                  "na",
                  "10",
                  10,
                  "10.1",
                  10.1,
                  "abc",
                  "$%TG")

set.seed(1114)
messy_df <- data.frame(var1 = messy_vector,
                       var2 = sample(messy_vector),
                       var3 = sample(messy_vector))

```


```{r vis-guess-messy-df, fig.show='hold', out.width='50%'}

vis_guess(messy_df)
vis_dat(messy_df)

```

So here we see that there are many different kinds of data in your dataframe. As an analyst this might be a depressing finding. We can see this comparison above.

Here, you might just assume your data is weird because it's all factors - or worse, not notice that this is a problem.

At the moment `vis_guess` is very slow. Please take this into consideration when you are using it on data with more than 1000 rows. We're looking into ways of making it faster, potentially using methods from the `parallel` package, or extending the c++ code from `readr:::collectorGuess`.

# Interactivity

You can make the plots in visdat by wrapping them in `plotly::ggplotly`:

```{r intx, eval = FALSE}

library(plotly)
ggplotly(vis_dat(airquality))
ggplotly(vis_miss(airquality))
ggplotly(vis_guess(airquality))

```

In the future these will have their own functions, written in plotly with nice
standardised on-hover behaviour. If you would like to see how these work,
please see the [development version on GitHub](https://github.com/ropensci/visdat).

# Future work

Future work from here is focussed on making `visdat` more stable, improving the 
speed of plotting, and adding interactive versions for each function.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-miss.R
\name{vis_miss}
\alias{vis_miss}
\title{Visualise a data.frame to display missingness.}
\usage{
vis_miss(
  x,
  cluster = FALSE,
  sort_miss = FALSE,
  show_perc = TRUE,
  show_perc_col = TRUE,
  large_data_size = 9e+05,
  warn_large_data = TRUE
)
}
\arguments{
\item{x}{a data.frame}

\item{cluster}{logical. TRUE specifies that you want to use hierarchical
clustering (mcquitty method) to arrange rows according to missingness.
FALSE specifies that you want to leave it as is. Default value is FALSE.}

\item{sort_miss}{logical. TRUE arranges the columns in order of missingness.
Default value is FALSE.}

\item{show_perc}{logical. TRUE now adds in the \\% of missing/complete data
in the whole dataset into the legend. Default value is TRUE.}

\item{show_perc_col}{logical. TRUE adds in the \\% missing data in a given
column into the x axis. Can be disabled with FALSE. Default value is TRUE.}

\item{large_data_size}{integer default is 900000 (given by
`nrow(data.frame) * ncol(data.frame)``). This can be changed. See
note for more details.}

\item{warn_large_data}{logical - warn if there is large data? Default is TRUE
see note for more details}
}
\value{
\code{ggplot2} object displaying the position of missing values in the
dataframe, and the percentage of values missing and present.
}
\description{
\code{vis_miss} provides an at-a-glance ggplot of the missingness inside a
dataframe, colouring cells according to missingness, where black indicates
a missing cell and grey indicates a present cell. As it returns a ggplot
object, it is very easy to customize and change labels.
}
\note{
Some datasets might be too large to plot, sometimes creating a blank
plot - if this happens, I would recommend downsampling the data, either
looking at the first 1,000 rows or by taking a random sample. This means
that you won't get the same "look" at the data, but it is better than
a blank plot! See example code for suggestions on doing this.
}
\examples{

vis_miss(airquality)

vis_miss(airquality, cluster = TRUE)

vis_miss(airquality, sort_miss = TRUE)

\dontrun{
# if you have a large dataset, you might want to try downsampling:
library(nycflights13)
library(dplyr)
flights \%>\%
  sample_n(1000) \%>\%
  vis_miss()

flights \%>\%
  slice(1:1000) \%>\%
  vis_miss()
}

}
\seealso{
\code{\link[=vis_dat]{vis_dat()}} \code{\link[=vis_guess]{vis_guess()}} \code{\link[=vis_expect]{vis_expect()}} \code{\link[=vis_cor]{vis_cor()}} \code{\link[=vis_compare]{vis_compare()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-cor.R
\name{vis_cor}
\alias{vis_cor}
\title{Visualise correlations amongst variables in your data as a heatmap}
\usage{
vis_cor(data, cor_method = "pearson", na_action = "pairwise.complete.obs", ...)
}
\arguments{
\item{data}{data.frame}

\item{cor_method}{correlation method to use, from \code{cor}: "a character
string indicating which correlation coefficient (or covariance) is to be
computed. One of "pearson" (default), "kendall", or "spearman": can be
abbreviated."}

\item{na_action}{The method for computing covariances when there are missing
values present. This can be "everything", "all.obs", "complete.obs",
"na.or.complete", or "pairwise.complete.obs" (default). This option is
taken from the \code{cor} function argument \code{use}.}

\item{...}{extra arguments you may want to pass to \code{cor}}
}
\value{
ggplot2 object
}
\description{
Visualise correlations amongst variables in your data as a heatmap
}
\examples{
vis_cor(airquality)
vis_cor(mtcars)
\dontrun{
# this will error
vis_cor(iris)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{add_vis_dat_pal}
\alias{add_vis_dat_pal}
\title{(Internal) Add a specific palette to a visdat plot}
\usage{
add_vis_dat_pal(vis_plot, palette)
}
\arguments{
\item{vis_plot}{visdat plot created using \code{vis_gather_}, \code{vis_extract_value}
and \code{vis_create_}}

\item{palette}{character "default", "qual" or "cb_safe". "default" (the
default) provides the stock ggplot scale for separating the colours. "qual"
uses an experimental qualitative colour scheme for providing distinct
colours for each Type. "cb_safe" is a set of colours that are appropriate
for those with colourblindness. "qual" and "cb_safe" are drawn from
http://colorbrewer2.org/.}
}
\value{
a visdat plot with a particular palette
}
\description{
(Internal) Add a specific palette to a visdat plot
}
\examples{

\dontrun{
# see internal use inside vis_guess and vis_dat
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-guess.R
\name{guess_type}
\alias{guess_type}
\title{(Internal) Guess the type of each individual cell in a dataframe}
\usage{
guess_type(x)
}
\arguments{
\item{x}{is a vector of values you want to guess}
}
\value{
a character vector that describes the suspected class. e.g., "10" is
an integer, "20.11" is a double, "text" is character, etc.
}
\description{
\code{vis_guess} uses \code{guess_type} to guess cell elements, like \code{fingerprint}.
}
\examples{
\dontrun{
guess_type(1)

guess_type("x")

guess_type(c("1", "0L"))

purrr::map_df(iris, guess_type)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-expect.R
\name{vis_expect}
\alias{vis_expect}
\title{Visualise whether a value is in a data frame}
\usage{
vis_expect(data, expectation, show_perc = TRUE)
}
\arguments{
\item{data}{a data.frame}

\item{expectation}{a formula following the syntax: \verb{~.x \{condition\}}.
For example, writing \code{~.x < 20} would mean "where a variable value is less
than 20, replace with NA", and \code{~.x \%in\% {vector}} would mean "where a
variable has values that are in that vector".}

\item{show_perc}{logical. TRUE now adds in the \\% of expectations are
TRUE or FALSE in the whole dataset into the legend. Default value is TRUE.}
}
\value{
a ggplot2 object
}
\description{
\code{vis_expect} visualises certain conditions or values in your data. For
example, If you are not sure whether to expect -1 in your data, you could
write: \code{vis_expect(data, ~.x == -1)}, and you can see if there are times
where the values in your data are equal to -1. You could also, for example,
explore a set of bad strings, or possible NA values and visualise where
they are using \code{vis_expect(data, ~.x \%in\% bad_strings)} where
\code{bad_strings} is a character vector containing bad strings  like \verb{N A}
\code{N/A} etc.
}
\examples{

dat_test <- tibble::tribble(
            ~x, ~y,
            -1,  "A",
            0,  "B",
            1,  "C",
            NA, NA
            )

vis_expect(dat_test, ~.x == -1)

vis_expect(airquality, ~.x == 5.1)

# explore some common NA strings

common_nas <- c(
"NA",
"N A",
"N/A",
"na",
"n a",
"n/a"
)

dat_ms <- tibble::tribble(~x,  ~y,    ~z,
                         "1",   "A",   -100,
                         "3",   "N/A", -99,
                         "NA",  NA,    -98,
                         "N A", "E",   -101,
                         "na", "F",   -1)

vis_expect(dat_ms, ~.x \%in\% common_nas)


}
\seealso{
\code{\link[=vis_miss]{vis_miss()}} \code{\link[=vis_dat]{vis_dat()}} \code{\link[=vis_guess]{vis_guess()}} \code{\link[=vis_cor]{vis_cor()}} \code{\link[=vis_compare]{vis_compare()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-binary.R
\name{vis_binary}
\alias{vis_binary}
\title{Visualise binary values}
\usage{
vis_binary(
  data,
  col_zero = "salmon",
  col_one = "steelblue2",
  col_na = "grey90",
  order = NULL
)
}
\arguments{
\item{data}{a data.frame}

\item{col_zero}{colour for zeroes, default is "salmon"}

\item{col_one}{colour for ones, default is "steelblue2"}

\item{col_na}{colour for NA, default is "grey90"}

\item{order}{optional character vector of the order of variables}
}
\value{
a ggplot plot of the binary values
}
\description{
Visualise binary values
}
\examples{
vis_binary(dat_bin)

# changing order of variables
# create numeric names
df <-  setNames(dat_bin, c("1.1", "8.9", "10.4"))
df

# not ideal
vis_binary(df)
# good - specify the original order
vis_binary(df, order = names(df))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{label_col_missing_pct}
\alias{label_col_missing_pct}
\title{Create labels for the columns containing the percent missing data}
\usage{
label_col_missing_pct(x, col_order_index)
}
\arguments{
\item{x}{data.frame}

\item{col_order_index}{the order of the columns}
}
\value{
data.frame containing the missingness percent down to 0.1 percent
}
\description{
Create labels for the columns containing the percent missing data
}
\note{
internal
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-expect.R
\name{expect_frame}
\alias{expect_frame}
\title{Create a dataframe to help visualise 'expected' values}
\usage{
expect_frame(data, expectation)
}
\arguments{
\item{data}{data.frame}

\item{expectation}{unquoted conditions or "expectations" to test}
}
\value{
data.frames where expectation are true
}
\description{
Create a dataframe to help visualise 'expected' values
}
\examples{
\dontrun{
dat_test <- tibble::tribble(
            ~x, ~y,
            -1,  "A",
            0,  "B",
            1,  "C"
            )

expect_frame(dat_test,
             ~ .x == -1)
             }
}
\author{
Stuart Lee and Earo Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{scale_01}
\alias{scale_01}
\title{Scale a vector between 0 and one.}
\usage{
scale_01(x)
}
\arguments{
\item{x}{numeric vector}
}
\value{
numeric vector between 0 and 1
}
\description{
Scale a vector between 0 and one.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visdat-package.r
\docType{package}
\name{visdat}
\alias{visdat}
\title{visdat}
\description{
visdat is a package that helps with the preliminary visualisation of data.
visdat makes it easy to visualise your whole dataset so that you can
visually identify problems.
}
\seealso{
It's main functions are:
\itemize{
\item \code{\link[=vis_dat]{vis_dat()}}
\item \code{\link[=vis_miss]{vis_miss()}}
\item \code{\link[=vis_guess]{vis_guess()}}
\item \code{\link[=vis_compare]{vis_compare()}}
\item \code{\link[=vis_expect]{vis_expect()}}
}

Learn more about visdat at \url{http://visdat.njtierney.com/articles/using_visdat.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-compare.R
\name{compare_print}
\alias{compare_print}
\title{(Internal) A utility function for \code{vis_compare}}
\usage{
compare_print(x)
}
\arguments{
\item{x}{a vector}
}
\description{
\code{compare_print} is an internal function that takes creates a dataframe with
information about where there are differences in the dataframe. This
function is used in \code{vis_compare}. It evaluates on the data \code{(df1 == df2)}
and (currently) replaces the "true" (the same) with "Same"
and FALSE with "Different", unless it is missing (coded as NA), in which
case it leaves it as NA.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-binary-data.R
\docType{data}
\name{dat_bin}
\alias{dat_bin}
\title{A small toy dataset of binary data with missings.}
\format{
A data frame with 100 rows and 3 variables:
\describe{
\item{x}{a binary variable with missing values.}
\item{y}{a binary variable with missing values.}
\item{z}{a binary variable with \strong{no} missing values.}
}
}
\usage{
dat_bin
}
\description{
A dataset containing binary values and missing values. It is created to
illustrate the usage of \code{\link[=vis_binary]{vis_binary()}}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{vis_extract_value_}
\alias{vis_extract_value_}
\title{(Internal) Add values of each row as a column}
\usage{
vis_extract_value_(x)
}
\arguments{
\item{x}{dataframe created from \code{vis_gather_}}
}
\value{
the x dataframe with the added column \code{value}.
}
\description{
This adds information about each row, so that when called by plotly, the
values are made visible on hover. Warnings are suppressed because \code{tidyr}
gives a warning about type coercion, which is fine.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-cor.R
\name{gather_cor}
\alias{gather_cor}
\title{(Internal) create a tidy dataframe of correlations suitable for plotting}
\usage{
gather_cor(data, cor_method = "pearson", na_action = "pairwise.complete.obs")
}
\arguments{
\item{data}{data.frame}

\item{cor_method}{correlation method to use, from \code{cor}: "a character
string indicating which correlation coefficient (or covariance) is to be
computed. One of "pearson" (default), "kendall", or "spearman": can be
abbreviated."}

\item{na_action}{The method for computing covariances when there are missing
values present. This can be "everything", "all.obs", "complete.obs",
"na.or.complete", or "pairwise.complete.obs" (default). This option is
taken from the \code{cor} function argument \code{use}.}
}
\value{
tidy dataframe of correlations
}
\description{
(Internal) create a tidy dataframe of correlations suitable for plotting
}
\examples{
gather_cor(airquality)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{test_if_dataframe}
\alias{test_if_dataframe}
\title{Test if input is a data.frame}
\usage{
test_if_dataframe(x)
}
\arguments{
\item{x}{object}
}
\value{
an error if input (x) is not a data.frame
}
\description{
Test if input is a data.frame
}
\examples{
\dontrun{
# success
test_if_dataframe(airquality)
#fail
test_if_dataframe(AirPassengers)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{fingerprint}
\alias{fingerprint}
\title{Take the fingerprint of a data.frame - find the class or return NA}
\usage{
fingerprint(x)
}
\arguments{
\item{x}{a vector}
}
\description{
\code{fingerprint} is an internal function that takes the "fingerprint" of a
dataframe, and currently replaces the contents (x) with the class of a
given object, unless it is missing (coded as \code{NA}), in which case it leaves
it as \code{NA}. The name "fingerprint" is taken from the csv-fingerprint, of
which the package, \code{visdat}, is based upon
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-value.R
\name{vis_value}
\alias{vis_value}
\title{Visualise the value of data values}
\usage{
vis_value(data, na_colour = "grey90", viridis_option = "D")
}
\arguments{
\item{data}{a data.frame}

\item{na_colour}{a character vector of length one describing what colour
you want the NA values to be. Default is "grey90"}

\item{viridis_option}{A character string indicating the colormap option to
use. Four options are available: "magma" (or "A"), "inferno" (or "B"),
"plasma" (or "C"), "viridis" (or "D", the default option) and "cividis"
(or "E").}
}
\value{
a ggplot plot of the values
}
\description{
Visualise all of the values in the data on a 0 to 1 scale. Only works on
numeric data - see examples for how to subset to only numeric data.
}
\examples{

vis_value(airquality)
vis_value(airquality, viridis_option = "A")
vis_value(airquality, viridis_option = "B")
vis_value(airquality, viridis_option = "C")
vis_value(airquality, viridis_option = "E")
\dontrun{
library(dplyr)
diamonds \%>\%
  select_if(is.numeric) \%>\%
  vis_value()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-compare.R
\name{vis_compare}
\alias{vis_compare}
\title{Visually compare two dataframes and see where they are different.}
\usage{
vis_compare(df1, df2)
}
\arguments{
\item{df1}{The first dataframe to compare}

\item{df2}{The second dataframe to compare to the first.}
}
\value{
\code{ggplot2} object displaying which values in each data frame are
present in each other, and which are not.
}
\description{
\code{vis_compare}, like the other \verb{vis_*} families, gives an at-a-glance ggplot
of a dataset, but in this case, hones in on visualising \strong{two} different
dataframes of the same dimension, so it takes two dataframes as arguments.
}
\examples{

# make a new dataset of iris that contains some NA values
aq_diff <- airquality
aq_diff[1:10, 1:2] <- NA
vis_compare(airquality, aq_diff)
}
\seealso{
\code{\link[=vis_miss]{vis_miss()}} \code{\link[=vis_dat]{vis_dat()}} \code{\link[=vis_guess]{vis_guess()}} \code{\link[=vis_expect]{vis_expect()}} \code{\link[=vis_cor]{vis_cor()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-expect.R
\name{expect_guide_label}
\alias{expect_guide_label}
\title{(Internal) Label the legend with the percent of missing data}
\usage{
expect_guide_label(x)
}
\arguments{
\item{x}{is a dataframe passed from \code{vis_miss(x)}.}
}
\value{
a \code{tibble} with two columns \code{p_miss_lab} and \code{p_pres_lab},
containing the labels to use for present and missing. A dataframe is
returned because I think it is a good style habit compared to a list.
}
\description{
\code{miss_guide_label} is an internal function to label the legend of \code{vis_miss}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{vis_gather_}
\alias{vis_gather_}
\title{(Internal) Gather rows into a format appropriate for grid visualisation}
\usage{
vis_gather_(x)
}
\arguments{
\item{x}{a dataframe}
}
\value{
data.frame gathered to have columns "variables", "valueType", and a
row id called "rows".
}
\description{
(Internal) Gather rows into a format appropriate for grid visualisation
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-typical-data.r
\docType{data}
\name{typical_data}
\alias{typical_data}
\title{A small toy dataset of imaginary people}
\format{
A data frame with 5000 rows and 11 variables:
\describe{
\item{ID}{Unique identifier for each individual, a sequential character
vector of zero-padded identification numbers (IDs). see ?wakefield::id}
\item{Race}{Race for each individual, "Black", "White", "Hispanic",
"Asian", "Other", "Bi-Racial", "Native", and "Hawaiin", see
?wakefield::race}
\item{Age}{Age of each individual, see ?wakefield::age}
\item{Sex}{Male or female, see ?wakefield::sex }
\item{Height(cm)}{Height in centimeters, see ?wakefield::height}
\item{IQ}{vector of intelligence quotients (IQ), see ?wakefield::iq}
\item{Smokes}{whether or not this person smokes, see ?wakefield::smokes}
\item{Income}{Yearly income in dollars, see ?wakefield::income}
\item{Died}{Whether or not this person has died yet., see ?wakefield::died}
}
}
\usage{
typical_data
}
\description{
A dataset containing information about some randomly generated people,
created using the excellent \code{wakefield} package. It is created as
deliberately messy dataset.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-typical-data-large.R
\docType{data}
\name{typical_data_large}
\alias{typical_data_large}
\title{A small toy dataset of imaginary people}
\format{
A data frame with 300 rows and 49 variables:
\describe{
\item{Age}{Age of each individual, see ?wakefield::age for more info}
\item{Animal}{A vector of animals, see ?wakefield::animal}
\item{Answer}{A vector of "Yes" or "No"}
\item{Area}{A vector of living areas "Suburban", "Urban", "Rural"}
\item{Car}{names of cars - see ?mtcars}
\item{Children}{vector of number of children - see ?wakefield::children}
\item{Coin}{character vector of "heads" and "tails"}
\item{Color}{vector of vectors from "colors()"}
\item{Date}{vector of "important" dates for an individual}
\item{Death}{TRUE / FALSE for whether this person died}
\item{Dice}{6 sided dice result}
\item{DNA}{vector of GATC nucleobases}
\item{DOB}{birth dates}
\item{Dummy}{a 0/1 dummy var}
\item{Education}{education attainment level}
\item{Employment}{employee status}
\item{Eye}{eye colour}
\item{Grade}{percent grades}
\item{Grade_Level}{favorite school grade}
\item{Group}{control or treatment}
\item{hair}{hair colours - "brown", "black", "blonde", or "red"}
\item{Height}{height in cm}
\item{Income}{yearly income}
\item{Browser}{choice of internet browser}
\item{IQ}{intelligence quotient}
\item{Language}{random language of the world}
\item{Level}{levels between 1 and 4}
\item{Likert}{likert response - "strongly agree", "agree", and so on}
\item{Lorem_Ipsum}{lorem ipsum text}
\item{Marital}{marital status- "married", "divorced", "widowed", "separated", etc}
\item{Military}{miliary branch they are in}
\item{Month}{their favorite month}
\item{Name}{their name}
\item{Normal}{a random normal number}
\item{Political}{their favorite political party}
\item{Race}{their race}
\item{Religion}{their religion}
\item{SAT}{their SAT score}
\item{Sentence}{an uttered sentence}
\item{Sex_1}{sex of their first child}
\item{Sex_2}{sex of their second child}
\item{Smokes}{do they smoke}
\item{Speed}{their median speed travelled in a car}
\item{State}{the last state they visited in the USA}
\item{String}{a random string they smashed out on the keyboard}
\item{Upper}{the last key they hit in upper case}
\item{Valid}{TRUE FALSE answer to a question}
\item{Year}{significant year to that individuals}
\item{Zip}{a zip code they have visited}
}
}
\usage{
typical_data_large
}
\description{
A wider dataset than \code{typical_data} containing information about some
randomly generated people, created using the excellent \code{wakefield}
package. It is created as deliberately odd / eclectic dataset.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-guess.R
\name{vis_guess}
\alias{vis_guess}
\title{Visualise type guess in a data.frame}
\usage{
vis_guess(x, palette = "default")
}
\arguments{
\item{x}{a data.frame}

\item{palette}{character "default", "qual" or "cb_safe". "default" (the
default) provides the stock ggplot scale for separating the colours.
"qual" uses an experimental qualitative colour scheme for providing
distinct colours for each Type. "cb_safe" is a set of colours that are
appropriate for those with colourblindness. "qual" and "cb_safe" are drawn
from http://colorbrewer2.org/.}
}
\value{
\code{ggplot2} object displaying the guess of the type of values in the
data frame and the position of any missing values.
}
\description{
\code{vis_guess} visualises the class of every single individual cell in a
dataframe and displays it as ggplot object, similar to \code{vis_dat}. Cells
are coloured according to what class they are and whether the values are
missing. \code{vis_guess} estimates the class of individual elements using
\code{readr::guess_parser}.  It may be currently slow on larger datasets.
}
\examples{

messy_vector <- c(TRUE,
                 "TRUE",
                 "T",
                 "01/01/01",
                 "01/01/2001",
                 NA,
                 NaN,
                 "NA",
                 "Na",
                 "na",
                 "10",
                 10,
                 "10.1",
                 10.1,
                 "abc",
                 "$\%TG")
set.seed(1114)
messy_df <- data.frame(var1 = messy_vector,
                       var2 = sample(messy_vector),
                       var3 = sample(messy_vector))
vis_guess(messy_df)
}
\seealso{
\code{\link[=vis_miss]{vis_miss()}} \code{\link[=vis_dat]{vis_dat()}} \code{\link[=vis_expect]{vis_expect()}} \code{\link[=vis_cor]{vis_cor()}} \code{\link[=vis_compare]{vis_compare()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-dat.R
\name{vis_dat}
\alias{vis_dat}
\title{Visualises a data.frame to tell you what it contains.}
\usage{
vis_dat(
  x,
  sort_type = TRUE,
  palette = "default",
  warn_large_data = TRUE,
  large_data_size = 9e+05
)
}
\arguments{
\item{x}{a data.frame object}

\item{sort_type}{logical TRUE/FALSE. When TRUE (default), it sorts by the
type in the column to make it easier to see what is in the data}

\item{palette}{character "default", "qual" or "cb_safe". "default" (the
default) provides the stock ggplot scale for separating the colours.
"qual" uses an experimental qualitative colour scheme for providing
distinct colours for each Type. "cb_safe" is a set of colours that are
appropriate for those with colourblindness. "qual" and "cb_safe" are drawn
from http://colorbrewer2.org/.}

\item{warn_large_data}{logical - warn if there is large data? Default is TRUE
see note for more details}

\item{large_data_size}{integer default is 900000 (given by
`nrow(data.frame) * ncol(data.frame)``). This can be changed. See
note for more details.}
}
\value{
\code{ggplot2} object displaying the type of values in the data frame and
the position of any missing values.
}
\description{
\code{vis_dat} gives you an at-a-glance ggplot object of what is inside a
dataframe. Cells are coloured according to what class they are and whether
the values are missing. As \code{vis_dat} returns a ggplot object, it is very
easy to customize and change labels, and customize the plot
}
\note{
Some datasets might be too large to plot, sometimes creating a blank
plot - if this happens, I would recommend downsampling the data, either
looking at the first 1,000 rows or by taking a random sample. This means
that you won't get the same "look" at the data, but it is better than
a blank plot! See example code for suggestions on doing this.
}
\examples{

vis_dat(airquality)

# experimental colourblind safe palette
vis_dat(airquality, palette = "cb_safe")
vis_dat(airquality, palette = "qual")

# if you have a large dataset, you might want to try downsampling:
\dontrun{
library(nycflights13)
library(dplyr)
flights \%>\%
  sample_n(1000) \%>\%
  vis_dat()

flights \%>\%
  slice(1:1000) \%>\%
  vis_dat()
}

}
\seealso{
\code{\link[=vis_miss]{vis_miss()}} \code{\link[=vis_guess]{vis_guess()}} \code{\link[=vis_expect]{vis_expect()}} \code{\link[=vis_cor]{vis_cor()}}
\code{\link[=vis_compare]{vis_compare()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{all_numeric}
\alias{all_numeric}
\title{(Internal) Are they all numeric columns?}
\usage{
all_numeric(x, ...)
}
\arguments{
\item{x}{data.frame}

\item{...}{optional extra inputs}
}
\value{
logical - TRUE means that there is a column with numerics, FALSE means that there is a column that is not numeric
}
\description{
(Internal) Are they all numeric columns?
}
\examples{

\dontrun{
all_numeric(airquality) # TRUE
all_numeric(iris) # FALSE
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{miss_guide_label}
\alias{miss_guide_label}
\title{Label the legend with the percent of missing data}
\usage{
miss_guide_label(x)
}
\arguments{
\item{x}{is a dataframe passed from vis_miss(x).}
}
\value{
a tibble with two columns \code{p_miss_lab} and \code{p_pres_lab},
containing the labels to use for present and missing. A dataframe is
returned because I think it is a good style habit compared to a list.
}
\description{
\code{miss_guide_label} is an internal function for vis_miss to label the legend.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{vis_create_}
\alias{vis_create_}
\title{(Internal) Create a boilerplate for visualisations of the \code{vis_} family}
\usage{
vis_create_(x)
}
\arguments{
\item{x}{a dataframe in longformat as transformed by \code{vis_gather_} and
\code{vis_extract_value}.}
}
\value{
a ggplot object
}
\description{
(Internal) Create a boilerplate for visualisations of the \code{vis_} family
}
