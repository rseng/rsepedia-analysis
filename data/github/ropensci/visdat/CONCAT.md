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

*Wow, no problems at all. :)*