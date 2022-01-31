
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mctq <a href = "https://docs.ropensci.org/mctq/"><img src = "man/figures/logo.png" align="right" height="139" /></a>

<!-- badges: start -->

[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/434_status.svg)](https://github.com/ropensci/software-review/issues/434)
[![CRAN
status](https://www.r-pkg.org/badges/version/mctq)](https://cran.r-project.org/package=mctq)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/mctq)](https://cran.r-project.org/package=mctq)
[![mctq status
badge](https://ropensci.r-universe.dev/badges/mctq)](https://ropensci.r-universe.dev)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![R-CMD-check](https://github.com/ropensci/mctq/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/mctq/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/mctq/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/mctq?branch=main)
[![License:
MIT](https://img.shields.io/badge/license-MIT-green)](https://choosealicense.com/licenses/mit/)
[![Buy Me A Coffee donate
button](https://img.shields.io/badge/buy%20me%20a%20coffee-donate-yellow.svg)](https://ko-fi.com/danielvartan)
<!-- badges: end -->

## Overview

`mctq` is an R package that provides a complete and consistent toolkit
to process the Munich ChronoType Questionnaire (MCTQ), a quantitative
and validated tool to assess chronotypes using peoples’ sleep behavior
presented by Till Roenneberg, Anna Wirz-Justice, and Martha Merrow in
[2003](https://doi.org/10.1177/0748730402239679). The aim of `mctq` is
to facilitate the work of sleep and chronobiology scientists with MCTQ
data while also helping with research reproducibility.

Learn more about the MCTQ questionnaire at
<https://www.thewep.org/documentations/mctq>.

### Wait, an R package for a questionnaire?

Although it may look like a simple questionnaire, MCTQ requires a lot of
date/time manipulation. This poses a challenge for many scientists,
being that most people have difficulties with date/time data, especially
when dealing with an extensive dataset. The `mctq` package comes to
address this issue.

`mctq` can handle the processing tasks for the three MCTQ versions
(standard, micro, and shift) with few dependencies, relying much of its
applications on the [lubridate](https://lubridate.tidyverse.org/) and
[hms](https://hms.tidyverse.org/) packages from
[tidyverse](https://www.tidyverse.org/). We also designed `mctq` with
the user experience in mind, by creating an interface that resembles the
way the questionnaire data is shown in MCTQ publications, and by
providing extensive and detailed documentation about each computation
proposed by the MCTQ authors. The package also includes several utility
tools, along with fictional datasets for testing and learning purposes.

## Prerequisites

You need to have some familiarity with the [R programming
language](https://www.r-project.org/) and with the
[lubridate](https://lubridate.tidyverse.org/) and
[hms](https://hms.tidyverse.org/) packages from
[tidyverse](https://www.tidyverse.org/) to use `mctq` main functions.

In case you don’t feel comfortable with R, we strongly recommend
checking Hadley Wickham and Garrett Grolemund free and online book [R
for Data Science](https://r4ds.had.co.nz/) and the Coursera course from
John Hopkins University [Data Science: Foundations using
R](https://www.coursera.org/specializations/data-science-foundations-r)
(free for audit students).

Please refer to the [lubridate](https://lubridate.tidyverse.org/) and
[hms](https://hms.tidyverse.org/) package documentation to learn more
about them. These two are essential packages to deal with date/time data
in R. We also recommend that you read the [Dates and
times](https://r4ds.had.co.nz/dates-and-times.html) chapter from Wickham
& Grolemund’s book [R for Data Science](https://r4ds.had.co.nz/).

## Installation

You can install the released version of `mctq` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("mctq")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/mctq")
```

## Usage

`mctq` makes use of the [lubridate](https://lubridate.tidyverse.org/)
and [hms](https://hms.tidyverse.org/) packages from
[tidyverse](https://www.tidyverse.org/), which provide special objects
to deal with date/time values in R. If your dataset does not conform to
this structure, you first need to convert your data to it. Please refer
to those package documentations to learn more about them.

Due to the circular nature of time, we strongly recommend that you use
appropriate temporal objects while dealing with date/time in R. That can
help you get rid of several computation mistakes while trying to adapt
your data from a base 10 to a system rooted in a base 12 numerical
system.

### Workdays and work-free days variables

After your data is set to start, just use the `mctq` functions below to
process it.

Note that the `mctq` functions uses a similar naming pattern to that
used in the MCTQ publications. That makes it easy to find and apply any
computation necessary.

-   `fd()`: compute MCTQ work-free days.
-   `so()`: compute MCTQ local time of sleep onset.
-   `gu()`: compute MCTQ local time of getting out of bed.
-   `sdu()`: compute MCTQ sleep duration.
-   `tbt()`: compute MCTQ total time in bed.
-   `msl()`: compute MCTQ local time of mid-sleep.
-   `napd()`: compute MCTQ nap duration (only for MCTQ Shift).
-   `sd24()`: compute MCTQ 24 hours sleep duration (only for MCTQ
    Shift).

Example:

``` r
# Local time of preparing to sleep on workdays
sprep_w <- c(hms::parse_hm("23:45"), hms::parse_hm("02:15"))
# Sleep latency or time to fall asleep after preparing to sleep on workdays
slat_w <- c(lubridate::dminutes(30), lubridate::dminutes(90))
# Local time of sleep onset on workdays
so(sprep_w, slat_w)
#> 00:15:00
#> 03:45:00
```

### Combining workdays and work-free days variables

For computations combining workdays and work-free days, use:

-   `sd_week()`: compute MCTQ average weekly sleep duration.
-   `sd_overall()`: compute MCTQ overall sleep duration (only for MCTQ
    Shift).
-   `sloss_week()`: compute MCTQ weekly sleep loss.
-   `le_week()`: compute MCTQ average weekly light exposure.
-   `msf_sc()`: compute MCTQ chronotype or corrected local time of
    mid-sleep on work-free days.
-   `sjl_rel()` and `sjl()`: compute MCTQ social jet lag.
-   `sjl_weighted()`: compute MCTQ absolute social jetlag across all
    shifts (only for MCTQ Shift).

Example:

``` r
# Local time of mid-sleep on workdays
msw <- c(hms::parse_hm("02:05"), hms::parse_hm("04:05"))
# Local time of mid-sleep on work-free days
msf <- c(hms::parse_hm("23:05"), hms::parse_hm("08:30"))
# Relative social jetlag
sjl_rel(msw, msf)
#> [1] "-10800s (~-3 hours)"  "15900s (~4.42 hours)"
```

See a quick tour of all MCTQ main functions
[here](https://docs.ropensci.org/mctq/articles/mctq.html).

### Utilities

`mctq` is also equipped with many utility functions. The package also
provides fictional datasets of the standard, micro, and shift MCTQ
versions for testing and learning purposes.

All functions are well documented, showing all the guidelines behind the
computations. Click
[here](https://docs.ropensci.org/mctq/reference/index.html) to see a
list of them.

## Citation

If you use `mctq` in your research, please consider citing it. We put a
lot of work to build and maintain a free and open-source R package. You
can find the `mctq` citation below.

``` r
citation("mctq")
#> 
#> To cite {mctq} in publications use:
#> 
#>   Vartanian, D., Benedito-Silva, A. A., & Pedrazzoli, M. (2021).
#>   {mctq}: an R package for the Munich ChronoType Questionnaire.
#>   https://docs.ropensci.org/mctq/
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Unpublished{,
#>     title = {{mctq}: an R package for the Munich ChronoType Questionnaire},
#>     author = {Daniel Vartanian and Ana Amelia Benedito-Silva and Mario Pedrazzoli},
#>     year = {2021},
#>     url = {https://docs.ropensci.org/mctq/},
#>     note = {Lifecycle: maturing},
#>   }
```

## Contributing

`mctq` is a community project, everyone is welcome to contribute. Take a
moment to review our [Guidelines for
Contributing](https://docs.ropensci.org/mctq/CONTRIBUTING.html).

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Support `mctq`

[![ko-fi](https://ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/danielvartan)

Working with science in Brazil is a daily challenge. There are few
funding opportunities available and their value is not enough to live
on. Added to this, every day Brazilian science suffers from deep cuts in
funding, which requires researchers to always look for other sources of
income.

If this package helps you in any way or you simply want to support the
author’s work, please consider donating or even creating a membership
subscription (if you can!). Your support will help with the author’s
scientific pursuit and with the package maintenance.

To make a donation click on the [Ko-fi](https://ko-fi.com/danielvartan)
button above. Please indicate the `mctq` package in your donation
message.

Thank you!
<!--- https://devguide.ropensci.org/releasing.html -->
<!--- https://style.tidyverse.org/news.html -->
<!--- https://semver.org/ -->

# mctq (development version)

* The documentation was updated to add
[CRAN](https://cran.r-project.org/package=mctq) installation instructions.
* Some minor documentation issues were resolved.


# mctq 0.1.0 (2021-11-04)

* Initial [CRAN](https://cran.r-project.org/package=mctq) release. 🎉
* You can now install `mctq` with `install.packages("mctq")`.
* We decided to wait for a little while before releasing a `1.0.0` (stable) 
`mctq` version. We don't intend to make any breaking changes, but we think is
better to wait and see if the user community don't have any issues with the
features.


# mctq 0.0.0.9002 (pre-release)

* `mctq` is now a
[peer-reviewed](https://github.com/ropensci/software-review/issues/434) package
by @ropensci! 🎉
* The package repository was transferred to the @ropensci organization. All
links related to `mctq` have been changed. Old links have a redirect protocol to
point to the new repository and new website.


# mctq 0.0.0.9001 (pre-release)

* @jonkeane was added as a reviewer ('rev').
* @leocadio-miguel was added as a reviewer ('rev').

## Breaking changes

* `assign_date()` now returns only `Interval` objects.
* `convert()` and all `convert_*()` functions were removed. See a dedicated note
about this below.
* `ms()` was renamed to `msl()`. See a dedicated note about this below.
* `sd()` was renamed to `sdu()`. See a dedicated note about this below.
* `shortest_interval()` was renamed to `shorter_interval()`.
* `shorter_interval()` and `longer_interval()` now returns only `Interval`
objects.
* `sum_time()` now have different arguments and was divided in two functions:
`sum_time()` (for non-vectorized sums) and `vct_sum_time()` (for vectorized
sums).
* `sum_time()` now returns only `Duration` objects.
* To avoid any unknown compatibility issues, all packages on imports will now
require the latest version of them at the moment of release.

## New features

* `cycle_time()`, a function to cycle time spans, was introduced.

## Minor improvements and fixes

* `round_time()` is now a S3 generic.
* The user interface has gotten more stylish, thanks to the
[`cli`](https://cli.r-lib.org) package (now on imports).
* A new vignette was introduced, explaining why the `mctq` package use
`Duration` instead of `Period` (objects from the
[lubridate](https://lubridate.tidyverse.org/) package) as the default object for
time spans.
* Several typos were fixed in the documentation.
* Several functions were optimized.

## Note about removing `convert()`

`convert()` was created considering the user experience (sleep and chronobiology
scientists). Since most of them don't have much experience with R and that time
can have different types of representations (e.g., decimal hours, radian),
`convert()` aim was to help transpose those difficulties, posing as an
"universal translator" (🖖).

However, after much thought and consideration, we believe that the `convert()`
feature may be out of the `mctq` scope. It can maybe be part of another package
(a `lubritime` package perhaps? 😄). Other `mctq` tools, like
`shorter_interval()` and `sum_time()`, could also be a part of that package (but
are necessary in `mctq` for the time being). Hence, we decided to remove
`convert()` and to instruct the user to check the
[lubridate](https://lubridate.tidyverse.org/) and
[hms](https://hms.tidyverse.org/) packages for parsing/conversion.

## Note about renaming `sd()` and `ms()`

That was a tough, but necessary, call. Although we tried to preserve the
original author's naming pattern for the MCTQ functions, the name `sd` provokes
a dangerous name collision with the widely used `stats::sd()` function (standard
deviation) and the name `ms` provokes a name collision with `lubridate::ms()`
(a function for parsing minutes and seconds components). That's why we
decided to renamed them. `sdu()` and `msl()` are the only exceptions, all the
other `mctq` functions maintain a strong naming resemblance with the original
author's naming pattern.


# mctq 0.0.0.9000 (pre-release)

* Added a `NEWS.md` file to track changes to the package.
## Resubmission

This is a resubmission. Please see the notes about it below 
(based on the last submission response).

> Thanks,
> 
> Please omit the redundant "An R Package for the" from the title of your
> package

Response: Done.

> If there are references describing the methods in your package, please
> add these in the description field of your DESCRIPTION file in the form
> authors (year) <doi:...>
> authors (year) <arXiv:...>
> > authors (year, ISBN:...)
> or if those are not available: <https:...>
> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
> auto-linking.
> (If you want to add a title as well please put it in quotes: "Title")

Response: I plan to write an article about the package, but, at the moment, 
          there isn't one. I added the reference for the original MCTQ article 
          in the description field. The vignette 'mctq' gives a quick tour 
          of `mctq` main functions.

> Please add \value to .Rd files regarding exported methods and explain
> the functions results in the documentation. Please write about the
> structure of the output (class) and also what the output means. (If a
> function does not return a value, please document that too, e.g.
> \value{No return value, called for side effects} or similar)
> Missing Rd-tags:
>       qplot_walk.Rd: \value

Response: Done.

> It is sufficient to wrap interactive examples only in if (interactive())
> {}. Please remove the additional \dontrun{}.
> e.g. qplot_walk.Rd

Response: Done.

> Please fix and resubmit.
> 
> Best,
> Julia Haider

Response: Thank you! :)

## Test environments

* R-hub windows-x86_64-devel (r-devel)
* R-hub ubuntu-gcc-release (r-release)
* R-hub fedora-clang-devel (r-devel)

## R CMD check results

> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Daniel Vartanian <danvartan@gmail.com>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    ChronoType (2:36, 34:70)
    MCTQ (2:62, 35:20, 36:5)
    Merrow (38:30)
    Roenneberg (37:59)
    Wirz (38:5)
    chronotypes (36:57)
  New submission

0 errors √ | 0 warnings √ | 1 note x

## Notes

* 'ChronoType', 'MCTQ', 'Merrow' (refers to Martha Merrow), 'Roenneberg' (refers
to Till Roenneberg), 'Wirz' (refers to Anna Wirz-Justice), and 'chronotypes' are
not misspelled words.
# Contributing to `mctq`

<!-- This CONTRIBUTING.md was adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

First of all, thanks for considering contributing to `mctq`! 👍 It's people like you that make it rewarding for us - the project maintainers - to work on `mctq`. 😊

`mctq` is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/ropensci/mctq
[issues]: https://github.com/ropensci/mctq/issues
[discussions]: https://github.com/ropensci/mctq/discussions
[new_issue]: https://github.com/ropensci/mctq/issues/new
[new_discussion]: https://github.com/ropensci/mctq/discussions/new
[website]: https://docs.ropensci.org/mctq
[citation]: https://docs.ropensci.org/mctq/authors.html
[email]: mailto:danvartan@gmail.com

## Code of conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

This package generally uses the rOpenSci packaging guidelines for style and structure.

### Share the love ❤️

Think `mctq` is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using `mctq` for a paper you are writing? Consider [citing it][citation].

### Ask a question ⁉️

Using `mctq` and got stuck? Browse the [documentation][website] to see if you can find a solution. Still stuck? Post your question as an [new discussion on GitHub][new_discussion]. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email][email].

### Propose an idea 💡

Have an idea for a new `mctq` feature? Take a look at the [documentation][website] and [discussion list][discussions] to see if it isn't included or suggested yet. If not, suggest your idea as an [discussion on GitHub][new_discussion]. While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using `mctq` and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub][new_issue] so we can fix it. A good bug report makes it easier for us to do so, so please include:

* The content of `utils::sessionInfo()`.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug (tip: use [reprex](https://reprex.tidyverse.org/)).

### Improve the documentation 📖

Noticed a typo on the website? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

#### The website

[This website][website] is generated with [`pkgdown`](http://pkgdown.r-lib.org/). That means we don't have to write any html: content is pulled together from documentation in the code, vignettes, [Markdown](https://guides.github.com/features/mastering-markdown/) files, the package `DESCRIPTION` and `_pkgdown.yml` settings. If you know your way around `pkgdown`, you can [propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to improve documentation. If not, [start a discussion][new_discussion] so that we can point you in the right direction.

#### Function documentation

Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you want to improve a function description:

1. Go to `R/` directory in the [code repository][repo].
2. Look for the file with the function.
3. [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the roxygen comments (starting with `#'`).

### Contribute code 📝

Care to fix bugs or implement new functionality for `mctq`? Awesome! 👏 Have a look at the [issue list][issues] and leave a comment on the things you want to work on. See also the development guidelines below.

## Development guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream main`.
3. Open the RStudio project file (`.Rproj`).
4. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors, warnings and notes.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).

Also note that we use the [rOpenSci packaging guidelines](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md), [tidyverse design guide](https://principles.tidyverse.org/), and [tidyverse style guide](https://style.tidyverse.org/). Your code must conform to this principles and rules.
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
---
name: Issue template
about: rOpenSci issue template.
title: ''
labels: ''
assignees: ''
---

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
name: Task template
about: Use this template to create tasks.
title: ''
labels: ''
assignees: ''
---

> Comments in English, Spanish, or Portuguese are welcome.
> Comentarios en inglés, español o portugués son bienvenidos.
> Comentários em inglês, espanhol ou português são bem-vindos.

* __Due date__: `yyyy/mm/dd`.

---

## Instructions

## Materials

## Checklist

* [ ] `placeholder`.

## Other references
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

library(mctq)
library(lubridate)
library(hms)
```

# mctq <a href = "https://docs.ropensci.org/mctq/"><img src = "man/figures/logo.png" align="right" height="139" /></a>

<!-- badges: start -->
[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/434_status.svg)](https://github.com/ropensci/software-review/issues/434)
[![CRAN status](https://www.r-pkg.org/badges/version/mctq)](https://cran.r-project.org/package=mctq)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/mctq)](https://cran.r-project.org/package=mctq)
[![mctq status badge](https://ropensci.r-universe.dev/badges/mctq)](https://ropensci.r-universe.dev)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![R-CMD-check](https://github.com/ropensci/mctq/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/mctq/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/mctq/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/mctq?branch=main)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](https://choosealicense.com/licenses/mit/)
[![Buy Me A Coffee donate button](https://img.shields.io/badge/buy%20me%20a%20coffee-donate-yellow.svg)](https://ko-fi.com/danielvartan)
<!-- badges: end -->

## Overview

`mctq` is an R package that provides a complete and consistent toolkit to process the Munich ChronoType Questionnaire (MCTQ), a quantitative and validated tool to assess chronotypes using peoples' sleep behavior presented by Till Roenneberg, Anna Wirz-Justice, and Martha Merrow in [2003](https://doi.org/10.1177/0748730402239679). The aim of `mctq` is to facilitate the work of sleep and chronobiology scientists with MCTQ data while also helping with research reproducibility.

Learn more about the MCTQ questionnaire at <https://www.thewep.org/documentations/mctq>.

### Wait, an R package for a questionnaire?

Although it may look like a simple questionnaire, MCTQ requires a lot of date/time manipulation. This poses a challenge for many scientists, being that most people have difficulties with date/time data, especially when dealing with an extensive dataset. The `mctq` package comes to address this issue.

`mctq` can handle the processing tasks for the three MCTQ versions (standard, micro, and shift) with few dependencies, relying much of its applications on the [lubridate](https://lubridate.tidyverse.org/) and [hms](https://hms.tidyverse.org/) packages from [tidyverse](https://www.tidyverse.org/). We also designed `mctq` with the user experience in mind, by creating an interface that resembles the way the questionnaire data is shown in MCTQ publications, and by providing extensive and detailed documentation about each computation proposed by the MCTQ authors. The package also includes several utility tools, along with fictional datasets for testing and learning purposes.

## Prerequisites

You need to have some familiarity with the [R programming language](https://www.r-project.org/) and with the [lubridate](https://lubridate.tidyverse.org/) and [hms](https://hms.tidyverse.org/) packages from [tidyverse](https://www.tidyverse.org/) to use `mctq` main functions.

In case you don't feel comfortable with R, we strongly recommend checking Hadley Wickham and Garrett Grolemund free and online book [R for Data Science](https://r4ds.had.co.nz/) and the Coursera course from John Hopkins University [Data Science: Foundations using R](https://www.coursera.org/specializations/data-science-foundations-r) (free for audit students).

Please refer to the [lubridate](https://lubridate.tidyverse.org/) and [hms](https://hms.tidyverse.org/) package documentation to learn more about them. These two are essential packages to deal with date/time data in R. We also recommend that you read the [Dates and times](https://r4ds.had.co.nz/dates-and-times.html) chapter from Wickham & Grolemund's book [R for Data Science](https://r4ds.had.co.nz/).

## Installation

You can install the released version of `mctq` from [CRAN](https://CRAN.R-project.org) with:

```{r, eval = FALSE}
install.packages("mctq")
```

And the development version from [GitHub](https://github.com/) with:

``` {r, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/mctq")
```

## Usage

`mctq` makes use of the [lubridate](https://lubridate.tidyverse.org/) and [hms](https://hms.tidyverse.org/) packages from [tidyverse](https://www.tidyverse.org/), which provide special objects to deal with date/time values in R. If your dataset does not conform to this structure, you first need to convert your data to it. Please refer to those package documentations to learn more about them.

Due to the circular nature of time, we strongly recommend that you use appropriate temporal objects while dealing with date/time in R. That can help you get rid of several computation mistakes while trying to adapt your data from a base 10 to a system rooted in a base 12 numerical system.

### Workdays and work-free days variables

After your data is set to start, just use the `mctq` functions below to process it.

Note that the `mctq` functions uses a similar naming pattern to that used in the MCTQ publications. That makes it easy to find and apply any computation necessary.

* `fd()`: compute MCTQ work-free days.
* `so()`: compute MCTQ local time of sleep onset.
* `gu()`: compute MCTQ local time of getting out of bed.
* `sdu()`: compute MCTQ sleep duration.
* `tbt()`: compute MCTQ total time in bed.
* `msl()`: compute MCTQ local time of mid-sleep.
* `napd()`: compute MCTQ nap duration (only for MCTQ Shift).
* `sd24()`: compute MCTQ 24 hours sleep duration (only for MCTQ Shift).

Example:

```{r, message = FALSE}
# Local time of preparing to sleep on workdays
sprep_w <- c(hms::parse_hm("23:45"), hms::parse_hm("02:15"))
# Sleep latency or time to fall asleep after preparing to sleep on workdays
slat_w <- c(lubridate::dminutes(30), lubridate::dminutes(90))
# Local time of sleep onset on workdays
so(sprep_w, slat_w)
```

### Combining workdays and work-free days variables

For computations combining workdays and work-free days, use:

* `sd_week()`: compute MCTQ average weekly sleep duration.
* `sd_overall()`: compute MCTQ overall sleep duration (only for MCTQ Shift).
* `sloss_week()`: compute MCTQ weekly sleep loss.
* `le_week()`: compute MCTQ average weekly light exposure.
* `msf_sc()`: compute MCTQ chronotype or corrected local time of mid-sleep on work-free days.
* `sjl_rel()` and `sjl()`: compute MCTQ social jet lag.
* `sjl_weighted()`: compute MCTQ absolute social jetlag across all shifts (only for MCTQ Shift).

Example:

```{r, message = FALSE}
# Local time of mid-sleep on workdays
msw <- c(hms::parse_hm("02:05"), hms::parse_hm("04:05"))
# Local time of mid-sleep on work-free days
msf <- c(hms::parse_hm("23:05"), hms::parse_hm("08:30"))
# Relative social jetlag
sjl_rel(msw, msf)
```

See a quick tour of all MCTQ main functions [here](https://docs.ropensci.org/mctq/articles/mctq.html).

### Utilities

`mctq` is also equipped with many utility functions. The package also provides fictional datasets of the standard, micro, and shift MCTQ versions for testing and learning purposes.

All functions are well documented, showing all the guidelines behind the computations. Click [here](https://docs.ropensci.org/mctq/reference/index.html) to see a list of them.

## Citation

If you use `mctq` in your research, please consider citing it. We put a lot of work to build and maintain a free and open-source R package. You can find the `mctq` citation below.

```{r}
citation("mctq")
```

## Contributing

`mctq` is a community project, everyone is welcome to contribute. Take a moment to review our [Guidelines for Contributing](https://docs.ropensci.org/mctq/CONTRIBUTING.html).

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

## Support `mctq`

[![ko-fi](https://ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/danielvartan)

Working with science in Brazil is a daily challenge. There are few funding opportunities available and their value is not enough to live on. Added to this, every day Brazilian science suffers from deep cuts in funding, which requires researchers to always look for other sources of income.

If this package helps you in any way or you simply want to support the author's work, please consider donating or even creating a membership subscription (if you can!). Your support will help with the author's scientific pursuit and with the package maintenance.

To make a donation click on the [Ko-fi](https://ko-fi.com/danielvartan) button above. Please indicate the `mctq` package in your donation message.

Thank you!
---
title: "Why Duration and not Period?"
output: rmarkdown::html_vignette
description: >
  This article explains why the {mctq} package use `Duration` instead
  of `Period` (objects from the {lubridate} package) as the default object
  for time spans.
vignette: >
  %\VignetteIndexEntry{Why Duration and not Period?}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This article explains why the `mctq` package uses `Duration` instead of `Period` (objects from the [lubridate](https://lubridate.tidyverse.org/) package) as the default object for time spans.

## `Duration` versus `Period` objects

The [`lubridate`](https://lubridate.tidyverse.org/) package offers three types of objects for storing and manipulating time spans: [`Duration`](https://lubridate.tidyverse.org/reference/duration.html), [`Period`](https://lubridate.tidyverse.org/reference/period.html), and [`Interval`](https://lubridate.tidyverse.org/reference/interval.html).

To understand the difference between `Duration` and `Period` objects you must first remember that the timeline is not always consistent, as it can have irregularities caused by, for example, leap years, DST (Daylight Saving Time), or leap seconds. That's when `Period` objects differ from `Duration` objects.

[`Duration`](https://lubridate.tidyverse.org/reference/duration.html) objects represent time spans by their exact number of seconds. That is, a `Duration` object of 1 hour will always represent a 1-hour time span, even with possible timeline irregularities.

```{r echo = TRUE}
start <- lubridate::ymd_hms("2020-01-01 10:00:00", tz = "America/New_York")

start + lubridate::duration(1, units = "hour")

```

[`Period`](https://lubridate.tidyverse.org/reference/period.html) objects work a little bit differently. They are a special type of object developed by the `lubridate` team that represents "human units", ignoring possible timeline irregularities. That is to say that 1 day as `Period` can have different time spans when looking to a timeline after an irregular event.

To illustrate this behavior, take the case of a DST event, starting at 2016-03-13 01:00:00 EST.


```{r echo = TRUE}
start <- lubridate::ymd_hms("2016-03-13 01:00:00", tz = "America/New_York")

start + lubridate::duration(1, units = "hour")
start + lubridate::period(1, units = "hour")
```

You might ask: why the result is `NA` when adding 1 hour as a `Period` object? That's because `Period` objects ignore time irregularities. When the DST starts at `01:00:00` the timeline "jumps" to `03:00:00`, so the period from `02:00:00` to `02:59:59` doesn't exist.

```
base: 2016-03-13 01:00:00, tz = "America/New_York" 

                        DST + 1 hour
-----|---------------|               |---------------|----->
  01:00             NA           03:00           04:00

From the `Duration` perspective: base + 1 hour = 2016-03-13 03:00:00

     |-------------------------------|---------------|
                   1 hour                  1 hour

From the `Period` perspective: base + 1 hour = NA

     |---------------|---------------|---------------|
           1 hour          1 hour          1 hour
```

`Period` objects are useful when you need to consider the human units of time. For example:

```{r echo = TRUE}
start <- lubridate::ymd_hms("2016-03-13 01:00:00", tz = "America/New_York")

start + lubridate::duration(1, units = "day")
start + lubridate::period(1, units = "day")
```

In this case, `1 day`, by human standards, represents the same `time of day` of the next day. But, considering the DST event, that `1 day` has a time span of 23 hours.

You can learn more about [`lubridate`](https://lubridate.tidyverse.org/) time span objects in the [Dates and times](https://r4ds.had.co.nz/dates-and-times.html#time-spans) chapter from Wickham & Grolemund's book "R for Data Science".

## The MCTQ context

At first glance you might think that, since MCTQ was made for human respondents, the best representation for time spans would be the one that better represents "human units", right? That would be fine if we were talking about a time span in a timeline irregularity context, but MCTQ doesn’t deal with this scenario.

When using MCTQ, the interest is to measure the exact time span between one local time to another. By ignoring irregularities in the timeline, `Periods` produce a fluctuating time span, hence `Period` objects are not compatible with other time spans like objects (e.g., `hms`).

```{r echo = TRUE, eval = FALSE}
hms::parse_hm("10:00") + lubridate::period(1, units = "hours")

#> Error: Incompatible classes: <hms> + <Period>
```

In summary, `Period` objects were made considering a very specific context that doesn't apply to MCTQ. That's why `Duration` objects are the default object for time spans.
---
title: "Missing sections"
output: rmarkdown::html_vignette
description: >
  This article shows a possible workaround to deal with missing sections when
  working with a standard or micro Munich Chronotype Questionnaire (MCTQ).
vignette: >
  %\VignetteIndexEntry{Missing sections}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This article shows a possible workaround to a common problem encountered when dealing with a __standard__ Munich Chronotype Questionnaire (MCTQ) and __$\mu$MCTQ__ data.

It's a good idea to have the standard MCTQ questionnaire and the guidelines for the standard MCTQ variable computation open while reading this article/vignette. That way you can have a better understanding of the data objects we are going to deal with. You can download a copy of the MCTQ full standard version  [here](https://www.thewep.org/documentations/mctq/item/english-mctq-full). Click [here](https://www.thewep.org/documentations/mctq/item/mctq-variables) to download a copy of the guidelines for the standard MCTQ variables.

## Working around missing sections

Although the standard and micro versions of the MCTQ asks for respondents to complete the workdays and work-free days sections, even when they do not have a regular work schedule (`wd = 0`) or have a 7 day/week work schedule (`wd = 7`), some of them may still end skipping one of this parts of the questionnaire. In those cases, `sd_week()`, `sloss_week()`, `le_week()`, `msf_sc()`, `sjl_rel()`, and `sjl()` will produce `NA` (Not Available) as output. That's because those computations combine workdays and work-free days variables.

For those special standard and micro MCTQ cases, where one section is missing, a `NA` value is the correct output for the functions mentioned above when `wd` (number of workdays per week) are `wd > 0 & wd < 7`, but it may not be when `wd == 0` or `wd == 7`. While some researchers may just invalidate these latter cases, we propose a different approach.

To illustrate this approach, consider the following.

If a respondent __does not have a regular work schedule__ (`wd == 0`), __only answered the work-free days section__, and __does not use an alarm clock on their free days__ (i.e., `alarm_f == FALSE`), it would be fair to assume that there's no sleep correction (`sc`) to be made, therefore, their chronotype (`msf_sc`) must be equal to their midsleep on work-free days (`msf`).

Following this same line of thought, we can also say that:

* `sd_week` (average weekly sleep duration) must be equal to `sd_f` (sleep duration on work-free days) since the respondent does not have workdays.
* `sloss_week` (weekly sleep loss) must be equal to `0s` since there's no sleep debt.
* `le_week` (average weekly light exposure) must be equal to `le_f` (light exposure on work-free days) since there are no workdays.
* `sjl_rel` (relative social jet lag) and `sjl` (absolute social jet lag) must be equal to `0s` since there's no discrepancy between social and biological time.

Note that the [chronotype computation](https://www.thewep.org/documentations/mctq/item/mctq-variables) follows a similar line of thought.

The opposite scenario, i.e., when the respondent __works 7 days per week__ (`wd == 7`) and __only answered the workdays section__, can also have different outputs. `sloss_week()`, `msf_sc()`, `sjl_rel()`, and `sjl()` should still produce a `NA` as output since there's no way to know the real behavior of the respondent's sleep-wake cycle. But, according to this reasoning, `sd_week` and `le_week` can have different outputs.

* `sd_week` (average weekly sleep duration) must be equal to `sd_w` (sleep duration on workdays) since the respondent does not have work-free days.
* `le_week` (average weekly light exposure) must be equal to `le_w` (light exposure on work-free days) since the respondent does not have work-free days.

If you agree with this line of thought, we recommend creating dummy variables to identify those two situations and then change the case values as mentioned. You can see this procedure in action with the [data wrangling algorithms](https://github.com/ropensci/mctq/blob/main/data-raw/std_mctq.R#L923) made to produce the fictional `std_mctq` dataset, provided by the `mctq` package.

Please note that this workaround is not mentioned or endorsed by the MCTQ authors. If you use it, you must mention this reasoning in your methods section.
---
title: "Introduction to mctq"
output: rmarkdown::html_vignette
description: >
  This article takes a quick tour of all standard Munich Chronotype 
  Questionnaire (MCTQ) main functions from the `mctq` package.
vignette: >
  %\VignetteIndexEntry{Introduction to mctq}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


This article takes a quick tour of the Munich Chronotype Questionnaire (MCTQ) main functions from the `mctq` package. Please see the function documentation and other articles/vignettes for more details.

The same features presented here can also be used with $\mu$MCTQ. To make it easy for newcomers, some MCTQ$^{Shift}$ functions are not shown.

We assume that you already have [R](https://www.r-project.org/) installed and have some familiarity with R and MCTQ data. We also strongly recommend using [RStudio](https://www.rstudio.com/) as your IDE (Integrated Development Environment).

It's a good idea to have the standard MCTQ questionnaire open while reading this introduction. That way you can have a better understanding of the data objects we are going to deal with. You can download a copy of the MCTQ full standard version  [here](https://www.thewep.org/documentations/mctq/item/english-mctq-full).

## First things first

Let's start with the basics. The first thing you must do to use `mctq` is to have some MCTQ data and `mctq` installed and loaded.

Install `mctq` with:

```{r, eval = FALSE}
install.packages("mctq")
```

Great! We now must load the package to memory to start using it. Do this with:

```{r, warning = FALSE, message = FALSE}
library(mctq)
```

Now we just need to get some MCTQ data. For demonstration purposes, we're going to use a small and fictional raw standard MCTQ data provided by the `mctq` package.

This dataset already has valid values. As for any data analysis, you must have clean and valid data before using any analysis tool. If you don't know how to do that, we strongly recommend checking Hadley Wickham and Garrett Grolemund free and online book [R for data Science](https://r4ds.had.co.nz/) and the Coursera course from John Hopkins University [Data Science: Foundations using R](https://www.coursera.org/specializations/data-science-foundations-r) (free for audit students).

Teaching you how to load your data in R is outside the scope of this article. For that, we recommend checking the [readr](https://readr.tidyverse.org/) package from [tidyverse](https://www.tidyverse.org/).

Our fictional MCTQ data will be loaded with the code below. The naming of the variables follows the same naming pattern used in MCTQ publications. You can see the meaning of each variable by running `?std_mctq` in your console (you can also see it in this [link](https://docs.ropensci.org/mctq/reference/std_mctq.html)).

```{r, warning = FALSE, message = FALSE}
library(readr)

data <- readr::read_csv(mctq::raw_data("vignette_mctq.csv"),
                        col_types = readr::cols(.default = "c"))
```

## Converting your data

`mctq` makes use of the [lubridate](https://lubridate.tidyverse.org/) and [hms](https://hms.tidyverse.org/) packages from [tidyverse](https://www.tidyverse.org/), which provide special objects to deal with date/time values in R. If your dataset does not conform to this structure, you first need to convert your data to it.

Due to the circular nature of time, we strongly recommend that you use appropriate temporal objects while dealing with date/time in R. That can help you get rid of several computation mistakes while trying to adapt your data from a base 10 to a system rooted in a base 12 numerical system.

Teaching you how to parse/convert your data is outside the scope of this article. Please refer to the [lubridate](https://lubridate.tidyverse.org/) and [hms](https://hms.tidyverse.org/) package documentation to learn more about them. These two are essential packages to deal with date/time data in R. We also recommend that you read the [Dates and times](https://r4ds.had.co.nz/dates-and-times.html) chapter from Wickham & Grolemund's book "R for Data Science".

Here we are interested in two types of time objects: `Duration` objects, to store time spans, such as sleep latency, and `hms` objects, to store local time values, such as bedtime.

But first, let's take a look at the data, shall we? If you're using [RStudio](https://www.rstudio.com/), you can run the code above and then type `View(data)` in the console to explore it.

```{r}
data
```

As you can see, our data came in different formats. For example, the column `bt_w` (local time of going to bed on workdays) is in `hours:minutes` format, while `slat_w` (sleep latency on workdays) is a duration expressed in minutes.

Our fictional data will be parsed/converted with the code below. Please note that the [lubridate](https://lubridate.tidyverse.org/) and [hms](https://hms.tidyverse.org/) packages are equipped with easy tools and great documentation to help you with this task.

```{r, warnings = FALSE, message = FALSE}
library(dplyr)
library(hms)
library(lubridate)

data <- data %>% dplyr::mutate(
  dplyr::across(c("id", "wd"), as.integer),
  dplyr::across(dplyr::matches("^work$|^alarm_|^wake_|^reasons_f$"),
                as.logical),
  dplyr::across(dplyr::matches("^bt_|^sprep_|^se_"), hms::parse_hm),
  dplyr::across(dplyr::matches("^slat_|^si_"),
                ~ lubridate::dminutes(as.numeric(.x))),
  dplyr::across(dplyr::matches("^le_"), 
                ~ lubridate::as.duration(hms::parse_hm(.x)))
  )
```

Our data is now all set to start. Let's take a look at it.

```{r}
data
```

## Workdays and work-free days variables

`mctq` provides a complete and consistent toolkit to process Munich Chronotype Questionnaire (MCTQ) data. To start this process, we must first compute some MCTQ variables related to each section of the questionnaire.

We're going to use direct assigning while computing the MCTQ variables, just because is more straightforward for the examples. But, we recommend assigning variables to your dataset by using the `mutate()` function, included in the  [dplyr](https://dplyr.tidyverse.org/) package.

### `fd()`: Number of work-free days per week

`fd()` is a simple function that allows you to compute the difference between the number of days in a week (7) and the number of workdays per week (`wd`). It takes only `wd` as argument.

The output must be the total of free days a respondent has in a week.

```{r, warnings = FALSE, message = FALSE}
data$fd <- fd(data$wd)

# Comparing the result
data %>% dplyr::select(wd, fd)
```

### `so()`: Local time of sleep onset

`so()` allows you to compute the local time of sleep onset for workdays (`so_w`) and work-free days (`so_f`). It takes two arguments: `sprep` (local time of preparing to sleep) and `slat` (sleep latency or time to fall asleep after preparing to sleep).

The output must be the sum of `sprep` and `slat` in a circular time frame of 24 hours.

What is a circular time frame of 24 hours? Run `?cycle_time` in your console or click in this [link](https://docs.ropensci.org/mctq/reference/cycle_time.html) for a detail explanation.

```{r, warnings = FALSE, message = FALSE}
data$so_w <- so(data$sprep_w, data$slat_w)
data$so_f <- so(data$sprep_f, data$slat_f)

# Comparing the result
data %>% dplyr::select(sprep_w, slat_w, so_w, sprep_f, slat_f, so_f)
```

### `gu()`: Local time of getting out of bed

`gu()` allows you to compute the local time of getting out of bed for workdays (`gu_w`) and work-free days (`gu_f`). It takes two arguments: `se` (local time of sleep end) and `si` (sleep inertia).

The output must be the sum of `se` and `si` in a circular time frame of 24 hours.

Please note that, despite the name, `si` represents the time that the respondent takes to get up after sleep end. We decided to maintain the original names and abbreviations proposed by the MCTQ authors.

```{r, warnings = FALSE, message = FALSE}
data$gu_w <- gu(data$se_w, data$si_w)
data$gu_f <- gu(data$se_f, data$si_f)

# Comparing the result
data %>% dplyr::select(se_w, si_w, gu_w, se_f, si_f, gu_f)
```

### `sdu()`: Sleep duration

`sdu()` allows you to compute the sleep duration for workdays (`sd_w`) and work-free days (`sd_f`). It takes two arguments: `so` (local time of sleep onset) and `se` (local time of sleep end).

The output must be the difference between `se` and `so` in a circular time frame of 24 hours.

Please note that, although we tried to preserve the original authors' naming pattern for the MCTQ functions, the name `sd` provokes a dangerous name collision with the widely used [stats::sd](https://rdrr.io/r/stats/sd.html) (standard deviation) function. That's why we named it as `sdu`. `sdu()` and `msl()` are the only exceptions, all the other `mctq` functions maintain a strong naming resemblance with the original authors' naming pattern.

```{r, warnings = FALSE, message = FALSE}
data$sd_w <- sdu(data$so_w, data$se_w)
data$sd_f <- sdu(data$so_f, data$se_f)

# Comparing the result
data %>% dplyr::select(so_w, se_w, sd_w, so_f, se_f, sd_f)
```

### `tbt()`: Total time in bed

`tbt()` allows you to compute total time in bed for workdays (`tbt_w`) and work-free days (`tbt_f`). It takes two arguments: `bt` (local time of going to bed) and `gu` (local time of getting out of bed).

The output must be the difference between `gu` and `bt` in a circular time frame of 24 hours.

```{r, warnings = FALSE, message = FALSE}
data$tbt_w <- tbt(data$bt_w, data$gu_w)
data$tbt_f <- tbt(data$bt_f, data$gu_f)

# Comparing the result
data %>% dplyr::select(bt_w, gu_w, tbt_w, bt_f, gu_f, tbt_f)
```

### `msl()`: Local time of mid-sleep

`msl()` allows you to compute the local time of mid-sleep for workdays (`msw`) and work-free days (`msf`). It takes two arguments: `so` (local time of sleep onset) and `sd` (sleep duration).

The output must be the sum of `so` with the half of `sd` duration in a circular time frame of 24 hours.

```{r, warnings = FALSE, message = FALSE}
data$msw <- msl(data$so_w, data$sd_w)
data$msf <- msl(data$so_f, data$sd_f)

# Comparing the result
data %>% dplyr::select(so_w, sd_w, msw, so_f, sd_f, msf)
```

## Combining workdays and work-free days variables

We now have computed all MCTQ variables for each section of the questionnaire. Let's move to some variables that summarize our findings considering workdays and work-free days.

### `sd_week()`: Average weekly sleep duration

`sd_week()` allows you to compute the average weekly sleep duration. It takes three arguments: `sd_w` (sleep duration on workdays), `sd_f` (sleep duration on work-free days), and `wd` (number of workdays per week).

The output must be the weighted mean of `sd_w` and `sd_f`, with `wd` and `fd(wd)` as weights, in a circular time frame of 24 hours.

```{r, warnings = FALSE, message = FALSE}
data$sd_week <- sd_week(data$sd_w, data$sd_f, data$wd)

# Comparing the result
data <- data %>% dplyr::mutate(sd_week_rounded = mctq::round_time(sd_week))
data %>% dplyr::select(wd, sd_w, fd, sd_f, sd_week_rounded)
```

### `sloss_week()`: Weekly sleep loss

`sloss_week()` allows you to compute the weekly sleep loss. It takes three arguments: `sd_w` (sleep duration on workdays), `sd_f` (sleep duration on work-free days), and `wd` (number of workdays per week).

If `sd_week`(average weekly sleep duration) is greater than `sd_w`, the output must be the difference between `sd_week` and `sd_w` times `wd`. Else, it must return the difference between `sd_week` and `sd_f` times `fd(wd`) (number of free days per week). See `?sloss_week` to learn more.

```{r, warnings = FALSE, message = FALSE}
data$sloss_week <- sloss_week(data$sd_w, data$sd_f, data$wd)

# Comparing the result
data <- data %>% dplyr::mutate(
  sloss_week_rounded = mctq::round_time(sloss_week))
data %>% dplyr::select(wd, sd_w, fd, sd_f, sloss_week_rounded)
```

### `le_week()`: Average weekly light exposure

`le_week()` allows you to compute the average weekly light exposure. It takes three arguments: `le_w` (light exposure on workdays), `le_f` (light exposure on work-free days), and `wd` (number of workdays per week).

The output must be the weighted mean of `le_w` and `le_f`, with `wd` and `fd(wd)` as weights, in a circular time frame of 24 hours.

Please note that light exposure is measured only with the full version of the standard MCTQ.

```{r, warnings = FALSE, message = FALSE}
data$le_week <- le_week(data$le_w, data$le_f, data$wd)

# Comparing the result
data <- data %>% dplyr::mutate(le_week_rounded = mctq::round_time(le_week))
data %>% dplyr::select(wd, le_w, fd, le_f, le_week_rounded)
```

### `msf_sc()`: Chronotype or corrected local time of mid-sleep on work-free days

`msf_sc()` allows you to compute the chronotype, or corrected local time of mid-sleep on work-free days. It takes five arguments: `msf` (local time of mid-sleep on work-free days), `sd_w` (sleep duration on workdays), `sd_f` (sleep duration on work-free days), `sd_week`(average weekly sleep duration), and `alarm_f` (a `logical` object indicating if the respondent uses an alarm clock to wake up on work-free days).

If `sd_f` is less or equal than `sd_w`, the output must be `msf`. Else, it must return `msf` minus the difference between `sd_f` and `sd_week` divided by 2. `msf_sc` can only be computed if `alarm_f` is equal to `FALSE` (the function will return `NA` when `alarm_f == TRUE`).

`msf_sc` applies a correction to `msf`, removing an estimation of the effect from accumulated sleep debt on workdays that usually is compensated on work-free days. See `?msf_sc` to learn more.

```{r, warnings = FALSE, message = FALSE}
data$msf_sc <- msf_sc(data$msf, data$sd_w, data$sd_f, data$sd_week, 
                      data$alarm_f)

# Comparing the result
data <- data %>% dplyr::mutate(msf_sc_rounded = mctq::round_time(msf_sc))
data %>% dplyr::select(msf, msf_sc_rounded)
```

### `sjl_rel()`: Relative social jetlag

`sjl_rel()` allows you to compute the relative social jetlag. It takes at least two arguments: `msw` (local time of mid-sleep on workdays) and `msf` (local time of mid-sleep on work-free days).

The output must be the difference between `msf` and `msw` in a circular time frame of 24 hours.

In case you don't know, social jet lag is a concept developed by Wittmann et al. ([2006](https://doi.org/10.1080/07420520500545979)) that represents the discrepancy between social and biological time.

The difference described above may seem trivial or easy to compute, but it's not. See the `vignette("sjl-computation", package = "mctq")` to learn more.

```{r, warnings = FALSE, message = FALSE}
data$sjl_rel <- sjl_rel(data$msw, data$msf)

# Comparing the result
data %>% dplyr::select(msw, msf, sjl_rel)
```

### `sjl()`: Absolute social jetlag

`sjl()` allows you to compute the absolute social jetlag. This function works the same way as `sjl_rel`, but it returns an absolute value. In fact, `sjl_rel()` is just a wrapper function to `sjl()`, but with the `abs` argument set as `FALSE`.

If you already have `sjl_rel` computed, you don't really need to compute it twice, you can just use `abs(sjl_rel)`. That's what we're going to do with our data.

```{r, warnings = FALSE, message = FALSE}
data$sjl <- abs(data$sjl_rel)

# Comparing the result
data %>% dplyr::select(sjl_rel, sjl)
```


## Success!

We have now processed all the MCTQ standard variables proposed by the MCTQ authors.

Before we look at the final data, let's first reorder the columns to a nice logical order and remove some `*_rounded` variables that we created just for show.

```{r}
data <- data %>%
  dplyr::relocate(
            id, work, wd, fd,

            bt_w, sprep_w, slat_w, so_w, se_w, si_w, gu_w, alarm_w,
            wake_before_w, sd_w, tbt_w, le_w, msw,

            bt_f, sprep_f, slat_f, so_f, se_f, si_f, gu_f, alarm_f,
            reasons_f, reasons_why_f, sd_f, tbt_f, le_f, msf,

            sd_week, sloss_week, le_week, msf_sc, sjl_rel, sjl
            ) %>%
  dplyr::select(-dplyr::ends_with("_rounded"))
```

And our final dataset is ...

```{r}
data
```

If you're using [RStudio](https://www.rstudio.com/), you can run all the code showed above and then type `View(data)` in the console to explore the final result.

If you don't feel comfortable with the way `Duration` objects are printed, `mctq` provides a utility function to help you. Just use `pretty_mctq()` to get a better view.

```{r}
pretty_mctq(data, round = FALSE)
```

## Utilities

Before we end, it's important to note that `mctq` also provides several utility tools to help with your MCTQ data. Here's a list of them.

* `assign_date()`: Assign dates to two sequential hours.
* `cycle_time()`: Cycle time objects.
* `pretty_mctq()`: Make a MCTQ dataset more presentable.
* `qplot_walk()`: Walk through distribution plots.
* `random_mctq()`: Build a random MCTQ case.
* `raw_data()`: Get paths to `mctq` raw datasets.
* `round_time()`: Round time values.
* `shorter_interval()`: Find the shorter interval between two hours.
* `sum_time()`: Sum time objects.

`mctq` also provides fictional datasets of the standard, micro, and shift versions for testing and learning purposes.

* `std_mctq`: A fictional standard MCTQ dataset (`?std_mctq`).
* `micro_mctq`: A fictional $\mu$MCTQ dataset (`?micro_mctq`).
* `shift_mctq`: A fictional MCTQ$^{Shift}$ dataset (`?shift_mctq`).

We encouraged you to read the documentation of the features above. You may find it worth your time.
---
title: "Social jetlag computation"
output: rmarkdown::html_vignette
description: >
  This article shows some notes about different approaches that can be used to
  compute the social jetlag for the Munich Chronotype Questionnaire (MCTQ). It
  also explains how the `method` argument from the `sjl()` function  works.
vignette: >
  %\VignetteIndexEntry{Social jetlag computation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This article shows some notes about different approaches that can be used to compute the social jetlag ($SJL$) for the Munich Chronotype Questionnaire (MCTQ). It also explains how the `method` argument from the `sjl()` function works.

It's a good idea to have the standard MCTQ questionnaire and the guidelines for the standard MCTQ variable computation open while reading this article/vignette. That way you can have a better understanding of the data objects we are going to deal with. You can download a copy of the MCTQ full standard version  [here](https://www.thewep.org/documentations/mctq/item/english-mctq-full). Click [here](https://www.thewep.org/documentations/mctq/item/mctq-variables) to download a copy of the guidelines for the standard MCTQ variables.

## The two intervals problem

According to Roenneberg, Allebrandt, Merrow, & Vetter ([2012](https://doi.org/10.1016/j.cub.2012.03.038)) supplemental materials, the relative social jetlag ($SJL_{ rel}$), i.e., the discrepancy between social and biological time, must be computed as the difference between $MSF$ (local time of mid-sleep on work-free days) and $MSW$ (local time of mid-sleep on workdays).

$$SJL_{rel} = MSF - MSW$$

This simple equation may seem trivial until you take into account that you're are dealing with two time values detached from a timeline. In other words, $MSW$ and $MSF$ represent two moments in two different contexts (one on workdays and the other on work-free days).

If you dive into the MCTQ articles, you can see that this computation have two objectives:

1. To represent the distance between $MSW$ and $MSF$ (i.e., the discrepancy).
2. To establish what value comes before or after the other, representing that with a $+/-$ signal. That is, when $MSW$ comes before $MSF$, $SJL_{rel}$ must be positive, and when $MSW$ comes after $MSF$, $SJL_{rel}$ must be negative.

You can find the rationale about the $SJL_{rel}$ signal in Roenneberg, Pilz, Zerbini, & Winnebeck ([2019](https://doi.org/10.3390/biology8030054)) (see item "3.2 Social Jetlag Computation").

Most people have some trouble understanding this. To illustrate what we mean, let's visualize a timeline overlapping an $MSW$ and $MSF$ value:

```
             day 1                        day 2
    MSF                MSW       MSF                MSW
   05:00              21:00     05:00              21:00
-----|------------------|---------|------------------|----->
              16h           8h             16h
          longer int.  shorter int.    longer int.

```

Note that, while doing the representation above, we're dealing with the assumption that $MSW$ and $MSF$ can be represented in a two-day timeline since people don't usually sleep more than 24 hours (basic assumption).

As you can see, by overlapping two time values in a two-day timeline, we need to make a choice of what interval to use. For most people $MSF$ and $MSW$ are close to each other, so, usually, we are looking for the shorter interval between the two. But, in some extreme cases, usually when dealing with shift workers, $MSW$ and $MSF$ distance can surpass 12 hours, making the longer interval the correct answer.

To obtain the $SJL_{rel}$ signal we must check the start value of the interval. If the interval between $MSW$ and $MSF$ starts with $MSW$, that means that $MSW$ comes before $MSF$, hence, the signal must be positive. Else, if the interval between $MSW$ and $MSF$ starts with $MSF$, that means that $MSW$ comes after $MSF$, hence, the signal must be negative.

* Example 1: when $MSF - MSW$ makes a __positive__ $SJL_{rel}$


```
             day 1                        day 2
                       MSW       MSF                
                      21:00     05:00
------------------------|---------|------------------------>

```

* Example 2: when $MSF - MSW$ makes a __negative__ $SJL_{rel}$


```
             day 1                        day 2
                       MSF       MSW                
                      21:00     05:00
------------------------|---------|------------------------>

```

We call this the __two intervals problem__. It represents an unsolvable mathematical scenario, if you deprive it of the respondent context. That can generate minor errors when computing $SJL$, especially if you're dealing with large datasets.

## Methods for computing $SJL$

The `sjl()` function provides an argument called `method` that allows you to choose three different methods to deal with the two intervals problem. Here's how they work.

### `method = "difference"`

By using `method = "difference"`, `sjl()` will do the exact computation proposed by the MCTQ authors, i.e., $SJL$ will be computed as the linear difference between  $MSF$ and $MSW$.

Let's see some examples using this method.

* Example 3: using the `"difference"` method

$MSW = \text{04:00}$

$MSF = \text{06:00}$

$\text{Real difference: + 02:00}$

$MSF - MSW = \text{06:00} - \text{04:00} = \text{+ 02:00}$ (__right__)

* Example 4: using the `"difference"` method

$MSW = \text{23:00}$

$MSF = \text{03:00}$

$\text{Real difference: + 04:00}$

$MSF - MSW = \text{03:00} - \text{23:00} = \text{- 20:00}$ (__wrong__)

As you can see with the second example, the `"difference"` method uses a linear time frame approach, creating problems regarding the circularity of time.

### `method = "shorter"` (default method)

By using `method = "shorter"`, `sjl()` uses the shorter interval between $MSW$ and $MSF$.

This is the most reliable method we found to compute $SJL$, considering the context of the MCTQ data. However, it comes with a limitation: when $MSW$ and $MSF$ values distance themselves by more than 12 hours, `sjl()` can return a wrong output. From our experience with MCTQ data, a $SJL$ greater than 12 hours is highly improbable.

Let's see some examples using this method.

* Example 5: using the `"shorter"` method

$MSW = \text{04:00}$

$MSF = \text{06:00}$

$\text{Real difference: + 02:00}$

```
             day 1                        day 2
    MSF                MSW       MSF                MSW
   06:00              04:00     06:00              04:00
-----|------------------|---------|------------------|----->
             22h            2h             22h
         longer int.   shorter int.    longer int.

```

By using the shorter interval, $MSW$ comes before $MSF$, so $SJL_{rel}$ must be equal to $\text{+ 02:00}$ (__right__).

* Example 6: using the `"shorter"` method

$MSW = \text{23:00}$

$MSF = \text{03:00}$

$\text{Real difference: + 04:00}$

```
             day 1                        day 2
    MSF                MSW       MSF                MSW
   03:00              23:00     03:00              23:00
-----|------------------|---------|------------------|----->
             20h            4h             20h
         longer int.   shorter int.    longer int.

```

By using the shorter interval, $MSW$ comes before $MSF$, so $SJL_{rel}$ must be equal to $\text{+ 04:00}$  (__right__).

* Example 7: when the `"shorter"` method fails

$MSW = \text{12:00}$

$MSF = \text{23:00}$

$\text{Real difference: - 13:00}$

```
             day 1                        day 2
    MSW                     MSF                      MSW
   12:00                   23:00                    12:00
-----|-----------------------|------------------------|----->
               11h                       13h
           shorter int.              longer int.

```

By using the shorter interval, $MSW$ comes before $MSF$, so $SJL_{rel}$ must be equal to $\text{+ 11:00}$ (__wrong__).

You can see example 7 in the `shift_mctq` dataset provided by the `mctq` package (ID 39, on and after night shifts). That's the only MCTQ$^{Shift}$ case in `shift_mctq` where we think that the `"shorter"` method would fail.

### `method = "longer"`

By using `method = "longer"`, `sjl()` uses the longer interval between $MSW$ and $MSF$. It's just the opposite of the `"shorter"` method showed above.

## So, what method should I use?

We recommend that you always use the `"shorter"` method when computing $SJL_{rel}$ or $SJL$ (the default `sjl()` method).

In our tests, the `"shorter"` method demonstrated to be almost fail-safe. You really need to worry about the $SJL$ computation if you are dealing with shift workers.

When dealing with a large MCTQ$^{Shift}$ dataset, it will be very difficult to identify $SJL$ errors, unless you look case by case and check the results with your respondents. That's usually not a viable option. We recommend that you mention which method you use to compute $SJL$ and add it as a possible limitation of your results.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cycle_time.R
\name{cycle_time}
\alias{cycle_time}
\title{Cycle time objects}
\usage{
cycle_time(time, cycle, reverse = TRUE)
}
\arguments{
\item{time}{An object belonging to one of the following classes: \code{numeric},
\code{Duration}, \code{difftime}, or \code{hms}.}

\item{cycle}{A \code{numeric} or \code{Duration} object of length 1, equal or greater
than 0, indicating the cycle length in seconds (see Details to learn more).}

\item{reverse}{(optional) A \code{logical} value indicating if the function must
use a reverse cycle for negative values in \code{time} (see Details to learn
more) (default: \code{TRUE}).}
}
\value{
The same type of object of \code{time} cycled with the \code{cycle} parameter.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{cycle_time()} cycles time span objects in a predetermined cycle length,
converting linear time objects to a circular time frame.
}
\details{
\subsection{Linear versus circular time}{

Time can have different "shapes".

If the objective is to measure the duration (time span) of an event, time is
usually measured considering a linear frame, with a fixed point of
\href{https://en.wikipedia.org/wiki/Origin_(mathematics)}{origin}. In this
context, the time value distance itself to infinity in relation to the
origin.\preformatted{                                   B
                             |----------|
                                        A
                             |---------------------|
 - inf                                                inf +
<----------------------------|----------|----------|------->
 s                           0          5          10     s
                           origin

A + B = 10 + 5 = 15s
}

But that's not the only possible "shape" of time, as it can also be measured
in other contexts.

In a "time of day" context, time will be linked to the rotation of the
earth, "resetting" when a new rotation cycle starts. That brings a different
kind of shape to time: a circular shape. With this shape the time value
encounters the origin at the beginning and end of each cycle.\preformatted{               - <--- h ---> +
                    origin
                . . . 0 . . .
             .                 .
            .                   .
           .                     .
          .                       .
         .                         .
         18                        6
         .                         .
          .                       .
           .                     .
            .                   .
             .                 .
                . . . 12 . . .

18 + 6 = 0h
}

If we transpose this circular time frame to a linear one, it would look like
this:\preformatted{<----|---------------|---------------|---------------|----->
    0h              12h              0h             12h
  origin                           origin
}

Note that now the origin is not fix, but cyclical.

\code{cycle_time()} operates by converting linear time objects using a circular
approach relative to the cycle length (e.g, \code{cycle = 86400} (1 day)).
}

\subsection{Fractional time}{

\code{cycle_time()} uses the \code{\%\%} operator to cycle values. Hence, it can be
subject to catastrophic loss of accuracy if \code{time} is fractional and much
larger than \code{cycle}. A warning is given if this is detected.

\code{\%\%} is a \code{builtin} R function that operates like this:\preformatted{function(a, b) \{
    a - floor(a / b) * b
\}
}
}

\subsection{Negative time cycling}{

If \code{time} have a negative value and \code{reverse == FALSE}, \code{cycle_time()} will
perform the cycle considering the absolute value of \code{time} and return the
result with a negative signal.

However, If \code{time} have a negative value and \code{reverse == TRUE} (default),
\code{cycle_time()} will perform the cycle in reverse, relative to its origin.

Example: If you have a -30h time span in a reversed cycle of 24h, the result
will be 18h. By removing the full cycles of -30h you will get -6h (-30 + 24),
and -6h relative to the origin will be 18h.\preformatted{               - <--- h ---> +
                    origin
                . . . 0 . . .
              .                 .
            .                   .
           .                     .
          .                       .
         .                         .
    (-6) 18                        6 (-18)
         .                         .
          .                       .
           .                     .
            .                   .
             .                 .
                . . . 12 . . .
                    (-12)
}
}

\subsection{\code{Period} objects}{

\code{\link[lubridate:period]{Period}} objects are a special type of object
developed by the \link[lubridate:lubridate-package]{lubridate} team that
represents "human units", ignoring possible timeline irregularities. That is
to say that 1 day as \code{Period} can have different time spans, when looking to
a timeline after a irregularity event.

Since the time span of a \code{Period} object can fluctuate, \code{cycle_time()} don't
accept this kind of object. You can transform it to a \code{Duration} object and
still use the function, but beware that this can produce errors.

Learn more about \code{Period} objects in the \href{https://r4ds.had.co.nz/dates-and-times.html#periods}{Dates and times} chapter of
Wickham & Grolemund (n.d.).
}
}
\examples{
## Scalar example

time <- lubridate::dhours(25)
cycle <- lubridate::ddays(1)
cycle_time(time, cycle)
#> [1] "3600s (~1 hours)" # Expected

time <- lubridate::dhours(-25)
cycle <- lubridate::ddays(1)
reverse <- FALSE
cycle_time(time, cycle, reverse)
#> [1] "-3600s (~-1 hours)" # Expected

time <- lubridate::dhours(-25)
cycle <- lubridate::ddays(1)
reverse <- TRUE
cycle_time(time, cycle, reverse)
#> [1] "82800s (~23 hours)" # Expected

## Vector example

time <- c(lubridate::dmonths(24), lubridate::dmonths(13))
cycle <- lubridate::dyears(1)
cycle_time(time, cycle)
#> [1] "0s"                     "2629800s (~4.35 weeks)" # Expected

time <- c(lubridate::dmonths(24), lubridate::dmonths(-13))
cycle <- lubridate::dyears(1)
reverse <- FALSE
cycle_time(time, cycle, reverse)
#> [1] "0s"                       "-2629800s (~-4.35 weeks)" # Expected

time <- c(lubridate::dmonths(24), lubridate::dmonths(-13))
cycle <- lubridate::dyears(1)
reverse <- TRUE
cycle_time(time, cycle, reverse)
#> [1] "0s"                       "28927800s (~47.83 weeks)" # Expected
}
\references{
Wickham, H., & Grolemund, G. (n.d.). \emph{R for data science}. Sebastopol, CA:
O'Reilly Media. \url{https://r4ds.had.co.nz}
}
\seealso{
Other utility functions: 
\code{\link{assign_date}()},
\code{\link{pretty_mctq}()},
\code{\link{qplot_walk}()},
\code{\link{random_mctq}()},
\code{\link{raw_data}()},
\code{\link{round_time}()},
\code{\link{shorter_interval}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raw_data.R
\name{raw_data}
\alias{raw_data}
\title{Get paths to \code{mctq} raw datasets}
\usage{
raw_data(file = NULL)
}
\arguments{
\item{file}{(optional) a \code{character} object indicating the raw data file
name(s). If \code{NULL}, all raw data file names will be printed (default:
\code{NULL}).}
}
\value{
\itemize{
\item If \code{file = NULL}, a \code{character} object with all file names available.
\item If \code{file != NULL}, a string with the file name path.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{mctq} comes bundled with raw fictional datasets for testing and learning
purposes. \code{raw_data()} makes it easy to access their paths.
}
\examples{
\dontrun{
## To list all raw data file names available

raw_data()

## To get the file path from a specific raw data

raw_data(raw_data()[1])

raw_data("std_mctq.csv")}
}
\seealso{
Other utility functions: 
\code{\link{assign_date}()},
\code{\link{cycle_time}()},
\code{\link{pretty_mctq}()},
\code{\link{qplot_walk}()},
\code{\link{random_mctq}()},
\code{\link{round_time}()},
\code{\link{shorter_interval}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdu.R
\name{sd_overall}
\alias{sd_overall}
\title{Compute MCTQ overall sleep duration (only for MCTQ\eqn{^{Shift}}{ Shift})}
\usage{
sd_overall(sd_w, sd_f, n_w, n_f)
}
\arguments{
\item{sd_w}{A \code{Duration} object corresponding to the \strong{sleep duration between
two days in a particular shift} from a shift version of the MCTQ
questionnaire. You can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{sd_f}{A \code{Duration} object corresponding to the \strong{sleep duration between
two free days after a particular shift} from a shift version of the MCTQ
questionnaire. You can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{n_w}{An \link[checkmate:checkIntegerish]{integerish} \code{numeric} object or
an \code{integer} object corresponding to the \strong{number of days worked in a
particular shift within a shift cycle} from a shift version of the MCTQ
questionnaire.}

\item{n_f}{An \link[checkmate:checkIntegerish]{integerish} \code{numeric} object or
an \code{integer} object corresponding to the \strong{number of free days after a
particular shift within a shift cycle} from a shift version of the MCTQ
questionnaire.}
}
\value{
A \code{Duration} object corresponding to the vectorized weighted mean of
\code{sd_w} and \code{sd_f} with \code{n_w} and \code{n_f} as weights.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sd_overall()} computes the \strong{overall sleep duration in a particular shift}
for the shift version of the Munich Chronotype Questionnaire (MCTQ).

See \code{\link[=sd_week]{sd_week()}} to compute the average weekly sleep duration for the
standard and micro versions of the MCTQ.
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Juda, Vetter, & Roenneberg (2013) and The Worldwide Experimental Platform
(n.d.) guidelines for \code{sd_overall()} (\eqn{\emptyset SD^{M/E/N}}{OSD_M/E/N})
computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each shift section of the
questionnaire. If you're using the three-shift design proposed by the
authors, you need to compute three overall sleep duration (e.g.,
\eqn{\emptyset SD^M}{OSD_M}; \eqn{\emptyset SD^E}{OSD_E}; \eqn{\emptyset
SD^N}{OSD_N}).
\item The overall sleep duration is the weighted average of the shift-specific
mean sleep durations.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{Computation}{

\strong{\deqn{\frac{SD_W^{M/E/N} \times n_W^{M/E/N} + SD_F^{M/E/N} \times
n_F^{M/E/N}}{n_W^{M/E/N} + n_F^{M/E/N}}}{(SD_W_M/E/N * n_W_M/E/N +
SD_F_M/E/N * n_F_M/E/N) / (n_W_M/E/N + n_F_M/E/N)}}

Where:
\itemize{
\item \eqn{SD_W^{M/E/N}}{SD_W_M/E/N} = sleep duration between two days in a
particular shift.
\item \eqn{SD_F^{M/E/N}}{SD_F_M/E/N} = sleep duration between two free days after
a particular shift.
\item \eqn{n_W^{M/E/N}}{n_W_M/E/N} = number of days worked in a particular shift
within a shift cycle.
\item \eqn{n_F^{M/E/N}}{n_F_M/E/N} = number of free days after a particular shift
within a shift cycle.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

sd_w <- lubridate::dhours(5)
sd_f <- lubridate::dhours(9)
n_w <- 2
n_f <- 2
sd_overall(sd_w, sd_f, n_w, n_f)
#> [1] "25200s (~7 hours)" # Expected

sd_w <- lubridate::dhours(3.45)
sd_f <- lubridate::dhours(10)
n_w <- 3
n_f <- 1
sd_overall(sd_w, sd_f, n_w, n_f)
#> [1] "18315s (~5.09 hours)" # Expected

sd_w <- lubridate::as.duration(NA)
sd_f <- lubridate::dhours(12)
n_w <- 4
n_f <- 4
sd_overall(sd_w, sd_f, n_w, n_f)
#> [1] NA # Expected

## Vector example

sd_w <- c(lubridate::dhours(4), lubridate::dhours(7))
sd_f <- c(lubridate::dhours(12), lubridate::dhours(9))
n_w <- c(3, 4)
n_f <- c(2, 4)
sd_overall(sd_w, sd_f, n_w, n_f)
#> [1] "25920s (~7.2 hours)" "28800s (~8 hours)"  # Expected

## Checking second output from vector example

if (requireNamespace("stats", quietly = TRUE)) {
    i <- 2
    x <- c(sd_w[i], sd_f[i])
    w <- c(n_w[i], n_f[i])
    lubridate::as.duration(stats::weighted.mean(x, w))
}
#> [1] "28800s (~8 hours)" # Expected

## Converting the output to `hms`

sd_w <- lubridate::dhours(4.75)
sd_f <- lubridate::dhours(10)
n_w <- 5
n_f <- 2
x <- sd_overall(sd_w, sd_f, n_w, n_f)
x
#> [1] "22500s (~6.25 hours)" # Expected
hms::as_hms(as.numeric(x))
#> 06:15:00 # Expected

## Rounding the output at the seconds level

sd_w <- lubridate::dhours(5.9874)
sd_f <- lubridate::dhours(9.3)
n_w <- 3
n_f <- 2
x <- sd_overall(sd_w, sd_f, n_w, n_f)
x
#> [1] "26324.784s (~7.31 hours)" # Expected
round_time(x)
#> [1] "26325s (~7.31 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/micro_mctq.R
\docType{data}
\name{micro_mctq}
\alias{micro_mctq}
\title{A fictional \eqn{\mu}MCTQ dataset}
\format{
A \code{\link[dplyr:reexports]{tibble}} with 17 columns and 50 rows:

\describe{
\item{id}{
A unique \code{integer} value to identify each respondent in the dataset.
\cr \cr
Type: Control.
\cr \cr
R class: \code{integer}.}

\item{shift_work}{
A \code{logical} value indicating if the respondent has been a shift- or
night-worker in the past three months.
\cr \cr
Statement (\code{EN}): "I have been a shift- or night-worker in the past three
months: Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{wd}{
Number of \strong{workdays} per week.
\cr \cr
Statement (\code{EN}): "Normally, I work ___ days/week".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{integer}.}

\item{fd}{
Number of \strong{work-free days} per week.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{integer}.}

\item{so_w}{
Local time of sleep onset on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "On WORKDAYS ... I normally fall asleep at ___ : ___
AM/PM (this is NOT when you get into bed, but rather when you fall
asleep)".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{se_w}{
Local time of sleep end on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "On WORKDAYS ... I normally wake up at ___ : ___ AM/PM
(this is NOT when you get out of bed, but rather when you wake up)".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{sd_w}{
Sleep duration on \strong{workdays}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{msw}{
Local time of mid-sleep on \strong{workdays}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{so_f}{
Local time of sleep onset on \strong{work-free days} when the respondent
\strong{doesn't} use an alarm clock to wake up.
\cr \cr
Statement (\code{EN}): "On WORK-FREE DAYS when I DON'T use an alarm clock ... I
normally fall asleep at ___ : ___ AM/PM (this is NOT when you get into bed,
but rather when you fall asleep)".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{se_f}{
Local time of sleep end on \strong{work-free days} when the respondent
\strong{doesn't} use an alarm clock to wake up.
\cr \cr
Statement (\code{EN}): "On WORK-FREE DAYS when I DON'T use an alarm clock ... I
normally wake up at ___ : ___ AM/PM (this is NOT when you get out of bed,
but rather when you wake up)".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{sd_f}{
Sleep duration on \strong{work-free days} when the respondent \strong{doesn't} use an
alarm clock to wake up.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{msf}{
Local time of mid-sleep on \strong{work-free days} when the respondent
\strong{doesn't} use an alarm clock to wake up.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{sd_week}{
Average weekly sleep duration.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sloss_week}{
Weekly sleep loss.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{msf_sc}{
Chronotype or corrected local time of mid-sleep on \strong{work-free days}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{sjl_rel}{
Relative social jetlag.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sjl}{
Absolute social jetlag.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}
}
}
\source{
Created by Daniel Vartanian (package author).
}
\usage{
micro_mctq
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

A fictional dataset, \strong{for testing and learning purposes}, composed of
basic/measurable and computed variables of the Munich Chronotype
Questionnaire (MCTQ) micro version.

This data was created following the guidelines in Ghotbi \emph{et.al} (2020), in
addition to the guidelines in Roenneberg, Wirz-Justice, & Merrow (2003),
Roenneberg, Allebrandt, Merrow, & Vetter (2012), and The Worldwide
Experimental Platform (n.d.). See the References and Details sections
to learn more.
}
\details{
\code{micro_mctq} is a tidied, validated, and transformed version of
\code{raw_data("micro_mctq.csv")}.
\subsection{Guidelines}{

To learn more about the Munich Chronotype Questionnaire (MCTQ),
see Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt,
Merrow, & Vetter (2012), Roenneberg et al. (2015), and Roenneberg, Pilz,
Zerbini, & Winnebeck (2019).

To know about different MCTQ versions, see Juda, Vetter, & Roenneberg
(2013) and Ghotbi et.al (2020).

If you're curious about the variable computations and want to have access to
the full questionnaire, see The Worldwide Experimental Platform (n.d.).
}

\subsection{Data building and data wrangling}{

This dataset was created by randomized sampling (see \code{\link[=random_mctq]{random_mctq()}})
and by manual insertions of special cases. Its purpose is to demonstrate
common cases and data issues that researchers may find in their MCTQ data, in
addition to be a suggested data structure for MCTQ data.

You can see the \code{micro_mctq} build and data wrangling processes
\href{https://github.com/ropensci/mctq/blob/main/data-raw/micro_mctq.R}{here}.
}

\subsection{Variable naming}{

The naming of the variables took into account the naming scheme used in MCTQ
publications, in addition to the guidelines of the \href{https://style.tidyverse.org/}{tidyverse style guide}.
}

\subsection{Variable classes}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the \link[hms:hms-package]{hms}
and \link[lubridate:lubridate-package]{lubridate} package.
}

\subsection{\code{Duration} objects}{

If you prefer to view \code{Duration} objects as \code{hms} objects, run
\code{pretty_mctq(micro_mctq)}.
}
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Keller, L. K., Fischer, D., Matera, J. L., Vetter, C., &
Winnebeck, E. C. (2015). Human activity and rest in situ. In A. Sehgal (Ed.),
\emph{Methods in Enzymology} (Vol. 552, pp. 257-283). London, UK: Academic Press.
\doi{10.1016/bs.mie.2014.11.028}.

Roenneberg, T., Pilz, L. K., Zerbini, G., & Winnebeck, E. C. (2019).
Chronotype and social jetlag: a (self-) critical review. \emph{Biology}, \emph{8}(3),
54. \doi{10.3390/biology8030054}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other datasets: 
\code{\link{shift_mctq}},
\code{\link{std_mctq}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sloss_week.R
\name{sloss_week}
\alias{sloss_week}
\title{Compute MCTQ weekly sleep loss}
\usage{
sloss_week(sd_w, sd_f, wd)
}
\arguments{
\item{sd_w}{A \code{Duration} object corresponding to the \strong{sleep duration on
workdays} from a standard or micro version of the MCTQ questionnaire. You
can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{sd_f}{A \code{Duration} object corresponding to the \strong{sleep duration on
work-free days} from a standard or micro version of the MCTQ
questionnaire. You can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{wd}{An \link[checkmate:checkIntegerish]{integerish} \code{numeric} object or
an \code{integer} object corresponding to the \strong{number of workdays per week}
from a standard or micro version of the MCTQ questionnaire.}
}
\value{
A \code{Duration} object corresponding to the weekly sleep loss.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sloss_week()} computes the \strong{weekly sleep loss} for the standard and micro
versions of the Munich Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012) and The Worldwide
Experimental Platform (n.d.) guidelines for \code{sloss_week()}
(\eqn{SLoss_{week}}{SLoss_week}) computation are as follows.
\subsection{Notes}{
\itemize{
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{Computation}{

\strong{\deqn{\textrm{If } SD_{week} > SD_W \; , \; (SD_{week} - SD_W) \times WD}{
If SD_week > SD_W, (SD_week - SD_W) * WD}}
\strong{\deqn{\textrm{If } SD_{week} \leq SD_W \; , \; (SD_{week} - SD_F) \times
FD}{If SD_week <= SD_W, (SD_week - SD_F) * FD}}

Where:
\itemize{
\item \eqn{SD_W} = sleep duration on workdays.
\item \eqn{SD_F} = sleep duration on work-free days.
\item \eqn{SD_{week}}{SD_week} = average weekly sleep duration.
\item \eqn{WD} = number of workdays per week ("I have a regular work schedule and
work ___ days per week").
\item \eqn{FD} = number of work-free days per week.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}
}

\section{Missing sections in standard and micro MCTQ versions}{


Although the standard and micro versions of the MCTQ asks for respondents to
complete the workdays and work-free days sections, even when they do not
have a regular work schedule (\code{wd = 0}) or have a 7 day/week work schedule
(\code{wd = 7}), some of them may still end skipping one of this parts of the
questionnaire. In those cases, \code{sd_week()}, \code{sloss_week()}, \code{le_week()},
\code{msf_sc()}, \code{sjl_rel()}, and \code{sjl()} will produce \code{NA} (Not Available) as
output. That's because those computations combine workdays and work-free days
variables.

For those special standard and micro MCTQ cases, where one section is
missing, a \code{NA} value is the correct output for the functions mentioned above
when \code{wd} (number of workdays per week) are \code{wd > 0 & wd < 7}, but it may not
be when \code{wd == 0} or \code{wd == 7}. There are different approaches to deal with
this issue. See \code{vignette("missing-sections", package = "mctq")} to learn
more.
}

\examples{
## Scalar example

sd_w <- lubridate::dhours(6.5)
sd_f <- lubridate::dhours(7)
wd <- 4
sloss_week(sd_w, sd_f, wd)
#> [1] "3085.71428571429s (~51.43 minutes)" # Expected

sd_w <- lubridate::dhours(7)
sd_f <- lubridate::dhours(8)
wd <- 5
sloss_week(sd_w, sd_f, wd)
#> [1] "5142.85714285714s (~1.43 hours)" # Expected

sd_w <- lubridate::dhours(NA)
sd_f <- lubridate::dhours(9.45)
wd <- 7
sloss_week(sd_w, sd_f, wd)
#> [1] NA # Expected

## Vector example

sd_w <- c(lubridate::dhours(7), lubridate::dhours(8))
sd_f <- c(lubridate::dhours(6.5), lubridate::dhours(8))
wd <- c(2, 0)
sloss_week(sd_w, sd_f, wd)
#> [1] "2571.42857142857s (~42.86 minutes)" "0s" # Expected

## Converting the output to `hms`

sd_w <- lubridate::dhours(4)
sd_f <- lubridate::dhours(5)
wd <- 3
x <- sloss_week(sd_w, sd_f, wd)
x
#> [1] "6171.42857142858s (~1.71 hours)" # Expected
hms::as_hms(as.numeric(x))
#> 01:42:51.428571 # Expected

## Rounding the output at the seconds level

sd_w <- lubridate::dhours(5.8743)
sd_f <- lubridate::dhours(7.4324)
wd <- 6
x <- sloss_week(sd_w, sd_f, wd)
x
#> [1] "4807.85142857144s (~1.34 hours)" # Expected
round_time(x)
#> [1] "4808s (~1.34 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift_mctq.R
\docType{data}
\name{shift_mctq}
\alias{shift_mctq}
\title{A fictional MCTQ\eqn{^{Shift}}{ Shift} dataset}
\format{
A \code{\link[dplyr:reexports]{tibble}} with 128 columns and 50 rows:

\describe{
\item{id}{
A unique \code{integer} value to identify each respondent in the dataset.
\cr \cr
Type: Control.
\cr \cr
R class: \code{integer}.}

\item{n_w_m}{
Number of days \strong{worked in morning shifts} within a shift cycle.
\cr \cr
Type: Basic.
\cr \cr
R class: \code{integer}.}

\item{bt_w_m}{
Local time of going to bed on workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "I go to bed at ___ o'clock'".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{sprep_w_m}{
Local time of preparing to sleep on workdays \strong{between two morning
shifts}.
\cr \cr
Statement (\code{EN}): "I actually get ready to fall asleep at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{slat_w_m}{
Sleep latency or time to fall asleep after preparing to sleep on workdays
\strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "I need ___ minutes to fall asleep".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{so_w_m}{
Local time of sleep onset on workdays \strong{between two morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{se_w_m}{
Local time of sleep end on workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "I wake up at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{tgu_w_m}{
Time to get up on workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "I get up after ___ minutes".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{gu_w_m}{
Local time of getting out of bed on workdays \strong{between two morning
shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{alarm_w_m}{
A \code{logical} value indicating if the respondent uses an alarm clock to wake
up on workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "I wake up at ___ o'clock: ( ___ ) with alarm ( ___ )
without alarm".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{reasons_w_m}{
A \code{logical} value indicating if the respondent has any particular reasons
for why they \strong{cannot} freely choose their sleep times on workdays
\strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "There are particular reasons why I \strong{cannot} freely
choose my sleep times on morning shifts: Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{reasons_why_w_m}{
Particular reasons for why the respondent cannot freely choose their sleep
times on workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "If "Yes": Child(ren)/pet(s) ( ___ ) Hobbies ( ___ )
Others, for example: ___".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{character}.}

\item{sd_w_m}{
Sleep duration on workdays \strong{between two morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{tbt_w_m}{
Total time in bed on workdays \strong{between two morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{msw_m}{
Local time of mid-sleep on workdays \strong{between two morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{nap_w_m}{
A \code{logical} value indicating if the respondent usually takes a nap on
workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "I usually take a nap: Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{napo_w_m}{
Local time of nap onset on workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "If "Yes": I take a nap from ___ o'clock to ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{nape_w_m}{
Local time of nap end on workdays \strong{between two morning shifts}.
\cr \cr
Statement (\code{EN}): "If "Yes": I take a nap from ___ o'clock to ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{napd_w_m}{
Nap duration on workdays \strong{between two morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sd24_w_m}{
24 hours sleep duration (sleep duration + nap duration) on workdays
\strong{between two morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{n_f_m}{
Number of free days \strong{after working in morning shifts} within a shift
cycle.
\cr \cr
Type: Basic.
\cr \cr
R class: \code{integer}.}

\item{bt_f_m}{
Local time of going to bed on work-free days \strong{between two free days after
morning shifts}.
\cr \cr
Statement (\code{EN}): "I go to bed at ___ o'clock'".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{sprep_f_m}{
Local time of preparing to sleep on work-free days \strong{between two free days
after morning shifts}.
\cr \cr
Statement (\code{EN}): "I actually get ready to fall asleep at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{slat_f_m}{
Sleep latency or time to fall asleep after preparing to sleep on work-free
days \strong{between two free days after morning shifts}.
\cr \cr
Statement (\code{EN}): "I need ___ minutes to fall asleep".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{so_f_m}{
Local time of sleep onset on work-free days \strong{between two free days after
morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{se_f_m}{
Local time of sleep end on work-free days \strong{between two free days after
morning shifts}.
\cr \cr
Statement (\code{EN}): "I wake up at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{tgu_f_m}{
Time to get up on work-free days \strong{between two free days after morning
shifts}.
\cr \cr
Statement (\code{EN}): "I get up after ___ minutes".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{gu_f_m}{
Local time of getting out of bed on work-free days \strong{between two free days
after morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{alarm_f_m}{
A \code{logical} value indicating if the respondent uses an alarm clock to wake
up on work-free days \strong{between two free days after morning shifts}.
\cr \cr
Statement (\code{EN}): "I wake up at ___ o'clock: ( ___ ) with alarm ( ___ )
without alarm".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{reasons_f_m}{
A \code{logical} value indicating if the respondent has any particular reasons
for why they \strong{cannot} freely choose their sleep times on work-free days
\strong{between two free days after morning shifts}.
\cr \cr
Statement (\code{EN}): "There are particular reasons why I \strong{cannot} freely
choose my sleep times on morning shifts: Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{reasons_why_f_m}{
Particular reasons for why the respondent cannot freely choose their sleep
times on work-free days \strong{between two free days after morning shifts}.
\cr \cr
Statement (\code{EN}): "If "Yes": Child(ren)/pet(s) ( ___ ) Hobbies ( ___ )
Others, for example: ___".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{character}.}

\item{sd_f_m}{
Sleep duration on work-free days \strong{between two free days after morning
shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{tbt_f_m}{
Total time in bed on work-free days \strong{between two free days after morning
shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{msf_m}{
Local time of mid-sleep on work-free days \strong{between two free days after
morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{nap_f_m}{
A \code{logical} value indicating if the respondent usually takes a nap on
work-free days \strong{between two free days after morning shifts}.
\cr \cr
Statement (\code{EN}): "I usually take a nap: Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{napo_f_m}{
Local time of nap onset on work-free days \strong{between two free days after
morning shifts}.
\cr \cr
Statement (\code{EN}): "If "Yes": I take a nap from ___ o'clock to ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{nape_f_m}{
Local time of nap end on work-free days \strong{between two free days after
morning shifts}.
\cr \cr
Statement (\code{EN}): "If "Yes": I take a nap from ___ o'clock to ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{napd_f_m}{
Nap duration on work-free days \strong{between two free days after morning
shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sd24_f_m}{
24 hours sleep duration (sleep duration + nap duration) on work-free days
\strong{between two free days after morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sd_overall_m}{
Overall sleep duration considering workdays \strong{between two morning shifts}
and work-free days \strong{between two free days after morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{msf_sc_m}{
Corrected local time of mid-sleep on work-free days \strong{between two free days
after morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{sjl_rel_m}{
Relative social jetlag considering workdays \strong{between two morning shifts}
and work-free days \strong{between two free days after morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sjl_m}{
Absolute social jetlag considering workdays \strong{between two morning shifts}
and work-free days \strong{between two free days after morning shifts}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{...}{
For brevity, the subsequent variables, except for \strong{sjl_weighted}
(described below), are not shown here. That's because they have
the same configurations of the variables shown above, differing only
by shift (\strong{evening shift} (\verb{_e}) and \strong{night shift} (\verb{_n})).}

\item{sjl_weighted}{
Absolute social jetlag across all shifts.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}
}
}
\source{
Created by Daniel Vartanian (package author).
}
\usage{
shift_mctq
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

A fictional dataset, for \strong{testing and learning purposes}, composed of
basic/measurable and computed variables of the Munich Chronotype
Questionnaire (MCTQ) shift version.

This data was created following the guidelines in Juda, Vetter, & Roenneberg
(2013), in addition to the guidelines found in Roenneberg, Wirz-Justice, &
Merrow (2003), Roenneberg, Allebrandt, Merrow, & Vetter (2012), and The
Worldwide Experimental Platform (n.d.). See the References and Details
sections to learn more.
}
\details{
\code{shift_mctq} is a tidied, validated, and transformed version of
\code{raw_data("shift_mctq.csv")}.
\subsection{Guidelines}{

To learn more about the Munich Chronotype Questionnaire (MCTQ),
see Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt,
Merrow, & Vetter (2012), Roenneberg et al. (2015), and Roenneberg, Pilz,
Zerbini, & Winnebeck (2019).

To know about different MCTQ versions, see Juda, Vetter, & Roenneberg
(2013) and Ghotbi et al. (2020).

If you're curious about the variable computations and want to have access to
the full questionnaire, see The Worldwide Experimental Platform (n.d.).
}

\subsection{Data building and data wrangling}{

This dataset was created by randomized sampling (see \code{\link[=random_mctq]{random_mctq()}})
and by manual insertions of special cases. Its purpose is to demonstrate
common cases and data issues that researchers may find in their MCTQ data, in
addition to be a suggested data structure for MCTQ data.

You can see the \code{shift_mctq} build and data wrangling processes
\href{https://github.com/ropensci/mctq/blob/main/data-raw/shift_mctq.R}{here}.
}

\subsection{Variable naming}{

The naming of the variables took into account the naming scheme used in MCTQ
publications, in addition to the guidelines of the \href{https://style.tidyverse.org/}{tidyverse style guide}.
}

\subsection{Variable classes}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the \link[hms:hms-package]{hms}
and \link[lubridate:lubridate-package]{lubridate} package.
}

\subsection{\code{Duration} objects}{

If you prefer to view \code{Duration} objects as \code{hms} objects, run
\code{pretty_mctq(shift_mctq)}.
}
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Keller, L. K., Fischer, D., Matera, J. L., Vetter, C., &
Winnebeck, E. C. (2015). Human activity and rest in situ. In A. Sehgal (Ed.),
\emph{Methods in Enzymology} (Vol. 552, pp. 257-283). London, UK: Academic Press.
\doi{10.1016/bs.mie.2014.11.028}.

Roenneberg, T., Pilz, L. K., Zerbini, G., & Winnebeck, E. C. (2019).
Chronotype and social jetlag: a (self-) critical review. \emph{Biology}, \emph{8}(3),
54. \doi{10.3390/biology8030054}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other datasets: 
\code{\link{micro_mctq}},
\code{\link{std_mctq}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sjl.R
\name{sjl_weighted}
\alias{sjl_weighted}
\title{Compute MCTQ absolute social jetlag across all shifts}
\usage{
sjl_weighted(sjl, n_w)
}
\arguments{
\item{sjl}{A \code{list} object with \code{Duration} elements corresponding to the
\strong{absolute social jetlag in each shift} from a shift version of the MCTQ
questionnaire (you can use \code{\link[=sjl]{sjl()}} to compute it). \code{sjl} elements and
values must be paired with \code{n} elements and values.}

\item{n_w}{A \code{list} object with \link[checkmate:checkIntegerish]{integerish}
\code{integer} or \code{double} elements corresponding to the \strong{number of days worked
in each shift within a shift cycle} from a shift version of the MCTQ
questionnaire. \code{n} elements and values must be paired with \code{sjl} elements
and values.}
}
\value{
A \code{Duration} object corresponding to the vectorized weighted mean of
\code{sjl} with \code{n_w} as weights.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sjl_weighted()} computes the \strong{absolute social jetlag across all shifts}
for the shift version of the Munich Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Operation}{


The shift version of the MCTQ was developed for shift-workers rotating
through morning-, evening-, and night-shifts, but it also allows adaptations
to other shift schedules (Juda, Vetter, & Roenneberg, 2013). For this reason,
\code{sjl_weighted()} must operate with any shift combination.

Considering the requirement above, \code{sjl_weighted()} was developed to only
accept list objects as arguments. For this approach to work, both \code{sjl} and
\code{n_w} arguments must be lists with paired elements and values, i.e., the
first element of \code{sjl} (e.g., \code{sjl_m}) must be paired with the first element
of \code{n_w} (e.g., \code{n_w_m}). The function will do the work of combining them
and output a weighted mean.
}

\section{Guidelines}{


Juda, Vetter, & Roenneberg (2013) and The Worldwide Experimental Platform
(n.d.) guidelines for \code{sjl_weighted()} (\eqn{\emptyset
SJL_{weighted}}{OSJL_weighted}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The absolute social jetlag across all shifts (\eqn{\emptyset
SJL_{weighted}}{OSJL_weighted}) is the weighted average of all absolute
social jetlags.
\item The authors describe an equation for a three-shift schedule, but this may
not be your case. That's why this function works a little bit differently
(see the Operation section), allowing you to compute a weighted average with
any shift combination.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{Computation}{

\strong{\deqn{\frac{| SJL^M | \times n_W^M + | SJL^E | \times n_W^E + | SJL^N |
\times n_W^N}{n_W^M + n_W^E + n_W^N}}{(| SJL_M | * n_W_M + | SJL_E | *
n_W_E + | SJL_N | * n_W_N) / (n_W_M + n_W_E + n_W_N)}}

Where:
\itemize{
\item \eqn{SJL^{M/E/N}}{SJL_M/E/N} = absolute social jetlag in each shift.
\item \eqn{n_W^{M/E/N}}{n_W_M/E/N} = number of days worked in each shift within a
shift cycle.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

sjl <- list(sjl_m = lubridate::dhours(1.25),
            sjl_e = lubridate::dhours(0.5),
            sjl_n = lubridate::dhours(3))
n_w <- list(n_w_m = 3, n_w_e = 1, n_w_n = 4)
sjl_weighted(sjl, n_w)
#> [1] "7312.5s (~2.03 hours)" # Expected

sjl <- list(sjl_m = lubridate::dhours(1.25),
            sjl_e = lubridate::as.duration(NA),
            sjl_n = lubridate::dhours(3))
n_w <- list(n_w_m = 3, n_w_e = 1, n_w_n = 4)
sjl_weighted(sjl, n_w)
#> [1] NA # Expected

## Vector example

sjl <- list(sjl_m = c(lubridate::dhours(2), lubridate::dhours(2.45)),
            sjl_e = c(lubridate::dhours(3.21), lubridate::as.duration(NA)),
            sjl_n = c(lubridate::dhours(1.2), lubridate::dhours(5.32)))
n_w <- list(n_w_m = c(1, 3), n_w_e = c(4, 1), n_w_n = c(3, 3))
sjl_weighted(sjl, n_w)
#> [1] "8298s (~2.31 hours)" NA # Expected

## Checking the first output from vector example

if (requireNamespace("stats", quietly = TRUE)) {
    i <- 1
    x <- c(sjl[["sjl_m"]][i], sjl[["sjl_e"]][i], sjl[["sjl_n"]][i])
    w <- c(n_w[["n_w_m"]][i], n_w[["n_w_e"]][i], n_w[["n_w_n"]][i])
    lubridate::as.duration(stats::weighted.mean(x, w))
}
#> [1] "8298s (~2.31 hours)" # Expected

## Converting the output to hms

sjl <- list(sjl_m = lubridate::dhours(0.25),
            sjl_e = lubridate::dhours(1.2),
            sjl_n = lubridate::dhours(4.32))
n_w <- list(n_w_m = 4, n_w_e = 2, n_w_n = 1)
x <- sjl_weighted(sjl, n_w)
x
#> [1] "3970.28571428571s (~1.1 hours)" # Expected
hms::as_hms(as.numeric(x))
#> 01:06:10.285714 # Expected

## Rounding the output at the seconds level

round_time(x)
#> [1] "3970s (~1.1 hours)" # Expected
round_time(hms::as_hms(as.numeric(x)))
#> 01:06:10 # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/le_week.R
\name{le_week}
\alias{le_week}
\title{Compute MCTQ average weekly light exposure}
\usage{
le_week(le_w, le_f, wd)
}
\arguments{
\item{le_w}{A \code{Duration} object corresponding to the \strong{light exposure on
workdays} from a standard version of the MCTQ questionnaire.}

\item{le_f}{A \code{Duration} object corresponding to the \strong{light exposure on
work-free days} from a standard version of the MCTQ questionnaire.}

\item{wd}{An \link[checkmate:checkIntegerish]{integerish} \code{numeric} object or
an \code{integer} object corresponding to the \strong{number of workdays per week}
from a standard version of the MCTQ questionnaire.}
}
\value{
A \code{Duration} object corresponding to the vectorized weighted mean of
\code{le_w} and \code{le_f} with \code{wd} and \code{fd(wd)} as weights.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{le_week()} computes the \strong{average weekly light exposure} for the standard
version of the Munich Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012) and The Worldwide
Experimental Platform (n.d.) guidelines for \code{le_week()}
(\eqn{LE_{week}}{LE_week}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The average weekly light exposure is the weighted average of the light
exposure on work and work-free days in a week.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{Computation}{

\strong{\deqn{\frac{LE_W \times WD + LE_F \times FD}{7}}{
(LE_W * WD + LE_F * FD) / 7}}

Where:
\itemize{
\item \eqn{LE_W} = light exposure on workdays.
\item \eqn{LE_F} = light exposure on work-free days.
\item \eqn{WD} = number of workdays per week ("I have a regular work schedule and
work ___ days per week").
\item \eqn{FD} = number of work-free days per week.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}
}

\section{Missing sections in standard and micro MCTQ versions}{


Although the standard and micro versions of the MCTQ asks for respondents to
complete the workdays and work-free days sections, even when they do not
have a regular work schedule (\code{wd = 0}) or have a 7 day/week work schedule
(\code{wd = 7}), some of them may still end skipping one of this parts of the
questionnaire. In those cases, \code{sd_week()}, \code{sloss_week()}, \code{le_week()},
\code{msf_sc()}, \code{sjl_rel()}, and \code{sjl()} will produce \code{NA} (Not Available) as
output. That's because those computations combine workdays and work-free days
variables.

For those special standard and micro MCTQ cases, where one section is
missing, a \code{NA} value is the correct output for the functions mentioned above
when \code{wd} (number of workdays per week) are \code{wd > 0 & wd < 7}, but it may not
be when \code{wd == 0} or \code{wd == 7}. There are different approaches to deal with
this issue. See \code{vignette("missing-sections", package = "mctq")} to learn
more.
}

\examples{
## Scalar example

le_w <- lubridate::dhours(1.5)
le_f <- lubridate::dhours(3.7)
wd <- 5
le_week(le_w, le_f, wd)
#> [1] "7662.85714285714s (~2.13 hours)" # Expected

le_w <- lubridate::dhours(3)
le_f <- lubridate::dhours(1.5)
wd <- 6
le_week(le_w, le_f, wd)
#> [1] "10028.5714285714s (~2.79 hours)" # Expected

le_w <- lubridate::dhours(5.6)
le_f <- lubridate::as.duration(NA)
wd <- 3
le_week(le_w, le_f, wd)
#> [1] NA # Expected

## Vector example

le_w <- c(lubridate::dhours(3), lubridate::dhours(2.45))
le_f <- c(lubridate::dhours(3), lubridate::dhours(3.75))
wd <- c(4, 5)
le_week(le_w, le_f, wd)
#> [1] "10800s (~3 hours)" # Expected
#> [2] "10157.1428571429s (~2.82 hours)" # Expected

## Checking second output from vector example

if (requireNamespace("stats", quietly = TRUE)) {
    i <- 2
    x <- c(le_w[i], le_f[i])
    w <- c(wd[i], fd(wd[i]))
    lubridate::as.duration(stats::weighted.mean(x, w))
}
#> [1] "10157.1428571429s (~2.82 hours)" # Expected

## Converting the output to `hms`

le_w <- lubridate::dhours(1.25)
le_f <- lubridate::dhours(6.23)
wd <- 3
x <- le_week(le_w, le_f, wd)
x
#> [1] "14744.5714285714s (~4.1 hours)" # Expected
hms::hms(as.numeric(x))
#> 04:05:44.571429 # Expected

## Rounding the output at the seconds level

le_w <- lubridate::dhours(3.4094)
le_f <- lubridate::dhours(6.2345)
wd <- 2
x <- le_week(le_w, le_f, wd)
x
#> [1] "19538.3828571429s (~5.43 hours)" # Expected
round_time(x)
#> [1] "19538s (~5.43 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msl.R
\name{msl}
\alias{msl}
\alias{msw}
\alias{msf}
\title{Compute MCTQ local time of mid-sleep}
\usage{
msl(so, sd)
}
\arguments{
\item{so}{A \code{hms} object corresponding to the \strong{local time of sleep onset}
from a standard, micro, or shift version of the MCTQ questionnaire. You can
use \code{\link[=so]{so()}} to compute it for the standard or shift version.}

\item{sd}{A \code{Duration} object corresponding to the \strong{sleep duration} from a
standard, micro, or shift version of the MCTQ questionnaire. You can use
\code{\link[=sdu]{sdu()}} to compute it for any MCTQ version.}
}
\value{
A \code{hms} object corresponding to the vectorized sum of \code{so} and \code{(sd / 2)} in a circular time frame of 24 hours.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{msl()} computes the \strong{local time of mid-sleep} for standard, micro, and
shift versions of the Munich Chronotype Questionnaire (MCTQ).

Please note that, although we tried to preserve the original authors' naming
pattern for the MCTQ functions, the name \code{ms} provokes a dangerous name
collision with the \code{\link[lubridate:hms]{lubridate::ms()}} function (a function for parsing minutes
and seconds components). That's why we named it \code{msl}. \code{msl()} and \code{\link[=sdu]{sdu()}}
are the only exceptions, all the other \code{mctq} functions maintain a strong
naming resemblance with the original authors' naming pattern.
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Ghotbi et al. (2020), Juda,
Vetter, & Roenneberg (2013), and The Worldwide Experimental Platform (n.d.)
guidelines for \code{msl()} (\eqn{MSW} or \eqn{MSF}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each section of the
questionnaire.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{For standard and micro versions of the MCTQ}{

\strong{\deqn{SO_{W/F} + \frac{SD_{W/F}}{2}}{SO_W/F + (SD_W/F / 2)}}

Where:
\itemize{
\item \eqn{SO_{W/F}}{SO_W/F} = local time of sleep onset on work \strong{or} work-free
days.
\item \eqn{SD_{W/F}}{SD_W/F} = sleep duration on work \strong{or} work-free days.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}

\subsection{For the shift version of the MCTQ}{

\strong{\deqn{SO_{W/F}^{M/E/N} + \frac{SD_{W/F}^{M/E/N}}{2}}{
SO_W/F_M/E/N + (SD_W/F_M/E/N / 2)}}

Where:
\itemize{
\item \eqn{SO_{W/F}^{M/E/N}}{SO_W/F_M/E/N} = local time of sleep onset between
two days in a particular shift \strong{or} between two free days after a
particular shift.
\item \eqn{SD_{W/F}^{M/E/N}}{SD_W/F_M/E/N} = sleep duration between two days in a
particular shift \strong{or} between two free days after a particular shift.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

so <- hms::parse_hm("23:30")
sd <- lubridate::dhours(8)
msl(so, sd)
#> 03:30:00 # Expected

so <- hms::parse_hm("01:00")
sd <- lubridate::dhours(10)
msl(so, sd)
#> 06:00:00 # Expected

so <- hms::as_hms(NA)
sd <- lubridate::dhours(7.5)
msl(so, sd)
#> NA # Expected

## Vector example

so <- c(hms::parse_hm("00:10"), hms::parse_hm("01:15"))
sd <- c(lubridate::dhours(9.25), lubridate::dhours(5.45))
msl(so, sd)
#> [1] 04:47:30 # Expected
#> [1] 03:58:30 # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gu.R
\name{gu}
\alias{gu}
\title{Compute MCTQ local time of getting out of bed}
\usage{
gu(se, si)
}
\arguments{
\item{se}{A \code{hms} object corresponding to the \strong{local time of sleep end}
from a standard or shift version of the MCTQ questionnaire.}

\item{si}{A \code{Duration} object corresponding to the \strong{sleep inertia} or
\strong{time to get up} from a standard or shift version of the MCTQ
questionnaire.}
}
\value{
A \code{hms} object corresponding to the vectorized sum of \code{se} and \code{si}
in a circular time frame of 24 hours.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{gu()} computes the \strong{local time of getting out of bed} for standard and
shift versions of the Munich Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Juda, Vetter, & Roenneberg
(2013), and The Worldwide Experimental Platform (n.d.) guidelines for \code{gu()}
(\eqn{GU}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each section of the
questionnaire.
\item MCTQ\eqn{^{Shift}}{ Shift} uses \eqn{TGU} (time to get up) instead of
\eqn{SI} (sleep inertia). For the purpose of this computation, both represent
the same thing.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{For standard and micro versions of the MCTQ}{

\strong{\deqn{SE_{W/F} + SI_{W/F}}{SE_W/F + SI_W/F}}

Where:
\itemize{
\item \eqn{SE_{W/F}}{SE_W/F} = local time of sleep end on work \strong{or} work-free
days.
\item \eqn{SI_{W/F}}{SI_W/F} = sleep inertia on work \strong{or} work-free days
("after ___ min, I get up").
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}

\subsection{For the shift version of the MCTQ}{

\strong{\deqn{SE_{W/F}^{M/E/N} + TGU_{W/F}^{M/E/N}}{SE_W/F_M/E/N + TGU_W/F_M/E/N}}

Where:
\itemize{
\item \eqn{SE_{W/F}^{M/E/N}}{SE_W/F_M/E/N} = local time of sleep end between two
days in a particular shift \strong{or} between two free days after a particular
shift.
\item \eqn{TGU_{W/F}^{M/E/N}}{TGU_W/F_M/E/N} = time to get up after sleep end
between two days in a particular shift \strong{or} between two free days after a
particular shift ("after ___ min, I get up").
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

gu(hms::parse_hm("08:00"), lubridate::dminutes(10))
#> 08:10:00 # Expected
gu(hms::parse_hm("11:45"), lubridate::dminutes(90))
#> 13:15:00 # Expected
gu(hms::as_hms(NA), lubridate::dminutes(90))
#> NA # Expected

## Vector example

se <- c(hms::parse_hm("12:30"), hms::parse_hm("23:45"))
si <- c(lubridate::dminutes(10), lubridate::dminutes(70))
gu(se, si)
#> 12:40:00 # Expected
#> 00:55:00 # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pretty_mctq.R
\name{pretty_mctq}
\alias{pretty_mctq}
\title{Make an MCTQ dataset more presentable}
\usage{
pretty_mctq(data, round = TRUE, hms = TRUE)
}
\arguments{
\item{data}{A \code{data.frame} object.}

\item{round}{(optional) a \code{logical} value indicating if \code{Duration} and
\code{hms} objects must be rounded at the level of seconds (default: \code{TRUE}).}

\item{hms}{(optional) a \code{logical} value indicating if \code{Duration} and
\code{difftime} objects must be converted to \code{hms} (default: \code{TRUE}).}
}
\value{
A transformed \code{data.frame} object, as indicated in the arguments.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{pretty_mctq()} helps you to transform your Munich Chronotype Questionnaire
(MCTQ) data in many ways. See the Arguments and Details section to learn
more.
}
\details{
\subsection{Rounding}{

Please note that by rounding MCTQ values you discard data. That is to say
that if you need to redo a computation, or do new ones, your values can be
off by a couple of seconds (see
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off error}).

Round your values only if and when you want to present them more clearly,
like in graphical representations. You can also round values to facilitate
data exporting to text formats (like \code{.csv}), but note that this will come
with a precision cost.

Note also that \code{pretty_mctq()} uses \link[=round_time]{round_time()} for
rounding. \code{round_time()} is based on \link[base:Round]{round()}, which uses the
IEC 60559 standard. For more information see \code{?round_time}.
}
}
\examples{
data <- data.frame(
    a = 1,
    b = lubridate::duration(1.12345),
    c = hms::hms(1.12345))

## Rounding time objects from `data`

pretty_mctq(data, round = TRUE, hms = FALSE)

## Converting non-`hms` time objects from `data` to `hms`

pretty_mctq(data, round = FALSE, hms = TRUE)
}
\seealso{
Other utility functions: 
\code{\link{assign_date}()},
\code{\link{cycle_time}()},
\code{\link{qplot_walk}()},
\code{\link{random_mctq}()},
\code{\link{raw_data}()},
\code{\link{round_time}()},
\code{\link{shorter_interval}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdu.R
\name{sdu}
\alias{sdu}
\title{Compute MCTQ sleep duration}
\usage{
sdu(so, se)
}
\arguments{
\item{so}{A \code{hms} object corresponding to the \strong{local time of sleep onset}
from a standard, micro, or shift version of the MCTQ questionnaire. You can
use \code{\link[=so]{so()}} to compute it for the standard or shift version.}

\item{se}{A \code{hms} object corresponding to the \strong{local time of sleep end}
from a standard, micro, or shift version of the MCTQ questionnaire.}
}
\value{
A \code{Duration} object corresponding to the vectorized difference
between \code{se} and \code{so} in a circular time frame of 24 hours.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sdu()} computes the \strong{sleep duration} for standard, micro, and shift
versions of the Munich Chronotype Questionnaire (MCTQ).

Please note that, although we tried to preserve the original authors' naming
pattern for the MCTQ functions, the name \code{sd} provokes a dangerous name
collision with the widely used \code{\link[stats:sd]{stats::sd()}} function (standard deviation).
That's why we named it \code{sdu}. \code{sdu()} and \code{\link[=msl]{msl()}} are the only exceptions,
all the other \code{mctq} functions maintain a strong naming resemblance with the
original authors' naming pattern.
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Ghotbi et al. (2020), Juda,
Vetter, & Roenneberg (2013), and The Worldwide Experimental Platform (n.d.)
guidelines for \code{sdu()} (\eqn{SD}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each section of the
questionnaire.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{For standard and micro versions of the MCTQ}{

\strong{\deqn{SE_{W/F} - SO_{W/F}}{SE_W/F - SO_W/F}}

Where:
\itemize{
\item \eqn{SE_{W/F}}{SE_W/F} = local time of sleep end on work \strong{or} work-free
days.
\item \eqn{SO_{W/F}}{SO_W/F}  = local time of sleep onset on work \strong{or}
work-free days.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}

\subsection{For the shift version of the MCTQ}{

\strong{\deqn{SE_{W/F}^{M/E/N} - SO_{W/F}^{M/E/N}}{SE_W/F_M/E/N - SO_W/F_M/E/N}}

Where:
\itemize{
\item \eqn{SE_{W/F}^{M/E/N}}{SE_W/F_M/E/N} = local time of sleep end between two
days in a particular shift \strong{or} between two free days after a particular
shift.
\item \eqn{SO_{W/F}^{M/E/N}}{SO_W/F_M/E/N}  = local time of sleep onset between
two days in a particular shift \strong{or} between two free days after a
particular shift.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

so <- hms::parse_hm("23:00")
se <- hms::parse_hm("08:00")
sdu(so, se)
#> [1] "32400s (~9 hours)" # Expected

so <- hms::parse_hm("02:00")
se <- hms::parse_hm("12:30")
sdu(so, se)
#> [1] "37800s (~10.5 hours)" # Expected

so <- hms::parse_hm("03:15")
se <- hms::as_hms(NA)
sdu(so, se)
#> [1] NA # Expected

## Vector example

so <- c(hms::parse_hm("04:12"), hms::parse_hm("21:20"))
se <- c(hms::parse_hm("14:30"), hms::parse_hm("03:45"))
sdu(so, se)
#> [1] "37080s (~10.3 hours)" "23100s (~6.42 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_mctq.R
\name{random_mctq}
\alias{random_mctq}
\title{Build a random MCTQ case}
\usage{
random_mctq(model = "standard")
}
\arguments{
\item{model}{A string indicating the data model to return. Valid values are:
\code{"standard"}, "\verb{shift"}, and \code{"micro"} (default: \code{"standard"}).}
}
\value{
A named \code{list} with elements representing each MCTQ basic/measurable
variable of the model indicated in the \code{model} argument.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{random_mctq} builds a fictional Munich Chronotype Questionnaire (MCTQ) case
composed of MCTQ basic/measurable variables.

This function is \strong{for testing and learning purposes only}. Please don't
misuse it.
}
\details{
The case structure (variable names and classes) are the same as the datasets
provided by the \code{mctq} package. See \link{std_mctq}, \link{micro_mctq} and
\link{shift_mctq} to learn more.
\subsection{Requirements}{

This function requires the \code{\link[stats:stats-package]{stats}} package. This
won't be an issue for most people since the package comes with a standard R
installation.

If you don't have the \code{\link[stats:stats-package]{stats}} package, you can
install it with \code{install.packages("stats")}.
}

\subsection{Cases}{

Random standard and micro MCTQ cases were created with the general
population in mind. The data was set to resemble the distribution parameters
shown in Roenneberg, Wirz-Justice, & Merrow (2003).

MCTQ\eqn{^{Shift}}{ Shift} random cases were created based on the shift
configuration from "Study Site 1" shown in Vetter, Juda, & Roenneberg (2012).
The data was set to resemble the distribution parameters shown in Juda,
Vetter, & Roenneberg (2013).

You can see more about the distribution parameters used
\href{https://github.com/ropensci/mctq/blob/main/data-raw/random_mctq.R}{here}.
}
}
\examples{
\dontrun{
random_mctq("standard")
random_mctq("micro")
random_mctq("shift")
}
}
\references{
Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

Vetter, C., Juda, M., & Roenneberg, T. (2012). The influence of internal time,
time awake, and sleep duration on cognitive performance in shiftworkers.
\emph{Chronobiology International}, \emph{29}(8), 1127-1138.
\doi{10.3109/07420528.2012.707999}.
}
\seealso{
Other utility functions: 
\code{\link{assign_date}()},
\code{\link{cycle_time}()},
\code{\link{pretty_mctq}()},
\code{\link{qplot_walk}()},
\code{\link{raw_data}()},
\code{\link{round_time}()},
\code{\link{shorter_interval}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdu.R
\name{napd}
\alias{napd}
\title{Compute MCTQ nap duration (only for MCTQ\eqn{^{Shift}}{ Shift})}
\usage{
napd(napo, nape)
}
\arguments{
\item{napo}{A \code{hms} object corresponding to the \strong{local time of nap onset}
from the shift version of the MCTQ questionnaire.}

\item{nape}{A \code{hms} object corresponding to the \strong{local time of nap end}
from the shift version of the MCTQ questionnaire.}
}
\value{
A \code{Duration} object corresponding to the vectorized difference
between \code{nape} and \code{napo} in a circular time frame of 24 hours.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{napd()} computes the \strong{nap duration} for the shift version of the Munich
Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Juda, Vetter & Roenneberg (2013) and The Worldwide Experimental Platform
(n.d.) guidelines for \code{napd()} (\eqn{NapD}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each shift section of the
questionnaire.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{Computation}{

\strong{\deqn{NapE_{W/F}^{M/E/N} - NapO_{W/F}^{M/E/N}}{
NapE_W/F_M/E/N - NapO_W/F_M/E/N}}

Where:
\itemize{
\item \eqn{NapO_{W/F}^{M/E/N}}{NapO_W/F_M/E/N} = local time of nap onset between
two days in a particular shift \strong{or} between two free days after a
particular shift ("I take a nap from ___ o'clock [...]").
\item \eqn{NapE_{W/F}^{M/E/N}}{NapE_W/F_M/E/N} = local time of nap end between
two days in a particular shift \strong{or} between two free days after a
particular shift ("[...] to ___ o'clock").
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

napo <- hms::parse_hm("12:30")
nape <- hms::parse_hm("14:20")
napd(napo, nape)
#> [1] "6600s (~1.83 hours)"" # Expected

napo <- hms::parse_hm("23:45")
nape <- hms::parse_hm("00:30")
napd(napo, nape)
#> [1] "2700s (~45 minutes)" # Expected

napo <- hms::parse_hm("10:20")
nape <- hms::as_hms(NA)
napd(napo, nape)
#> [1] NA # Expected

## Vector example

napo <- c(hms::parse_hm("01:25"), hms::parse_hm("23:50"))
nape <- c(hms::parse_hm("03:10"), hms::parse_hm("01:10"))
napd(napo, nape)
#> [1] "6300s (~1.75 hours)" "4800s (~1.33 hours)"  # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shorter_interval.R
\name{shorter_interval}
\alias{shorter_interval}
\alias{longer_interval}
\title{Find the shorter interval between two hours}
\usage{
shorter_interval(x, y, inverse = FALSE)

longer_interval(x, y)
}
\arguments{
\item{x, y}{A \code{hms} or \code{POSIXt} object.}

\item{inverse}{(optional) a \code{logical} value indicating if the function must
return an inverse output, i.e., the longer interval between \code{x} and \code{y}.}
}
\value{
\itemize{
\item If \code{inverse = FALSE} (default), an \code{Interval} object with the shorter
interval between \code{x} and \code{y}.
\item If \code{inverse = TRUE}, an \code{Interval} object with the longer interval between
\code{x} and \code{y}.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{shorter_interval()} finds and returns the shorter interval between two
\code{hms} or \code{POSIXt} object hours.

\code{longer_interval()} do the inverse of \code{shorter_interval()}, i.e.,
finds the longer interval between two hours. It's just a wrapper for
\code{shorter_interval(x, y, inverse = TRUE)}.
}
\details{
\subsection{The two intervals problem}{

Given two hours, \code{x} and \code{y}, in a two-day timeline, without date references,
there will be always two possible intervals between them, as illustrated
below.

To figure out what interval is the  shorter or the longer,
\code{shorter_interval()} verify two scenarios: 1. When \code{x} comes before \code{y}; and
2. when \code{x} comes after \code{y}. This only works if \code{x} value is smaller than
\code{y}, therefore, the function will make sure to swap \code{x} and \code{y} values if the
latter assumption is not true.

Because \code{shorter_interval()} objective is to find the shorter interval, if
\code{x} and \code{y} are equal, the shorter interval will have a length of 0 hours,
resulting in an interval from \code{x} to \code{x}. But, if \code{inverse = TRUE} or
\code{longer_interval()} is used instead, the latter condition will return a
interval with 24 hours of length (from \code{x} to \code{x} + 1 day).

In cases when \code{x} and \code{y} distance themselves by 12 hours, there will be no
shorter or longer interval (they will have equal length). In those cases,
\code{shorter_interval()} and \code{longer_interval()} will return the same value
(an interval of 12 hours).\preformatted{             day 1                        day 2
     x                  y         x                  y
   06:00              22:00     06:00              22:00
-----|------------------|---------|------------------|----->
              16h           8h             16h
          longer int.  shorter int.   longer int.

              day 1                      day 2
     y                   x       y                   x
   13:00               08:00   13:00               08:00
-----|-------------------|-------|-------------------|----->
              19h           5h            19h
          longer int.  shorter int.  longer int.

    x,y             x,y             x,y             x,y
     x               y               x               y
   10:00           10:00           10:00           10:00
-----|---------------|---------------|---------------|----->
    0h              0h              0h              0h
            24h             24h             24h

              day 1                      day 2
     y               x               y               x
   12:00           00:00           12:00           00:00
-----|---------------|---------------|---------------|----->
            12h             12h             12h
}
}

\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the \link[hms:hms-package]{hms}
and \link[lubridate:lubridate-package]{lubridate} package.
}

\subsection{Base date and timezone}{

\code{shorter_interval()} uses the
\href{https://en.wikipedia.org/wiki/Unix_time}{Unix epoch} (1970-01-01) date as
the start date for creating intervals.

The output will always have \code{"UTC"} set as timezone. Learn more about
time zones in \link[base:timezones]{base::timezone}.
}

\subsection{\code{POSIXt} objects}{

\code{POSIXt} objects passed as argument to \code{x} or \code{y} will be stripped of their
dates. Only the time will be considered.

Both \code{POSIXct} and \code{POSIXlt} are objects that inherits the class \code{POSIXt}.
Learn more about it in \link[base:DateTimeClasses]{base::DateTimeClasses}.
}

\subsection{\code{NA} values}{

\code{shorter_interval()} will return an \code{Interval} \code{NA}-\code{NA} if \code{x} or \code{y} are
\code{NA}.
}
}
\examples{
## Scalar example

x <- hms::parse_hm("23:00")
y <- hms::parse_hm("01:00")
shorter_interval(x, y)
#> [1] 1970-01-01 23:00:00 UTC--1970-01-02 01:00:00 UTC # Expected

x <- lubridate::as_datetime("1985-01-15 12:00:00")
y <- lubridate::as_datetime("2020-09-10 12:00:00")
shorter_interval(x, y)
#> [1] 1970-01-01 12:00:00 UTC--1970-01-01 12:00:00 UTC # Expected

## Vector example

x <- c(hms::parse_hm("15:30"), hms::parse_hm("21:30"))
y <- c(hms::parse_hm("19:30"), hms::parse_hm("04:00"))
shorter_interval(x, y)
#> [1] 1970-01-01 15:30:00 UTC--1970-01-01 19:30:00 UTC # Expected
#> [2] 1970-01-01 21:30:00 UTC--1970-01-02 04:00:00 UTC # Expected

## Finding the longer interval between two hours

x <- lubridate::parse_date_time("01:10:00", "HMS")
y <- lubridate::parse_date_time("11:45:00", "HMS")
shorter_interval(x, y, inverse = TRUE)
#> [1] 1970-01-01 11:45:00 UTC--1970-01-02 01:10:00 UTC # Expected

x <- lubridate::as_datetime("1915-02-14 05:00:00")
y <- lubridate::as_datetime("1970-07-01 05:00:00")
longer_interval(x, y)
#> [1] 1970-01-01 05:00:00 UTC--1970-01-02 05:00:00 UTC # Expected
}
\seealso{
Other utility functions: 
\code{\link{assign_date}()},
\code{\link{cycle_time}()},
\code{\link{pretty_mctq}()},
\code{\link{qplot_walk}()},
\code{\link{random_mctq}()},
\code{\link{raw_data}()},
\code{\link{round_time}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdu.R
\name{sd24}
\alias{sd24}
\title{Compute MCTQ 24 hours sleep duration (only for MCTQ\eqn{^{Shift}}{ Shift})}
\usage{
sd24(sd, napd, nap)
}
\arguments{
\item{sd}{A \code{Duration} object corresponding to the \strong{sleep duration} from
the shift version of the MCTQ questionnaire. You can use \code{\link[=sdu]{sdu()}} to
compute it.}

\item{napd}{A \code{Duration} object corresponding to the \strong{nap duration} from
the shift version of the MCTQ questionnaire. You can use \code{\link[=napd]{napd()}} to
compute it.}

\item{nap}{A \code{logical} value corresponding to the \strong{"I usually take a nap"}
from the shift version of the MCTQ questionnaire.}
}
\value{
\itemize{
\item If \code{nap == TRUE}, a \code{Duration} object corresponding to the vectorized sum
of \code{sd} and \code{napd} in a circular time frame of 24 hours.
\item If \code{nap == FALSE}, a \code{Duration} object equal to \code{sd}.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sd24()} computes the \strong{24 hours sleep duration} for the shift version of
the Munich Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Juda, Vetter & Roenneberg (2013) and The Worldwide Experimental Platform
(n.d.) guidelines for \code{sd24()} (\eqn{SD24}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each shift section of the
questionnaire.
\item If the respondent don't usually take a nap in a particular shift \strong{or}
between two free days after a particular shift, \code{sd24()} will return only
\eqn{SD_{W/F}^{M/E/N}}{SD_W/F_M/E/N}.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{Computation}{

\strong{\deqn{SD_{W/F}^{M/E/N} + NapD_{W/F}^{M/E/N}}{
SD_W/F_M/E/N + NapD_W/F_M/E/N}}

Where:
\itemize{
\item \eqn{SD_{W/F}^{M/E/N}}{SD_W/F_M/E/N} = sleep duration between two days in a
particular shift \strong{or} between two free days after a particular shift.
\item \eqn{NapD_{W/F}^{M/E/N}}{NapD_W/F_M/E/N} = nap duration between two days in
a particular shift \strong{or} between two free days after a particular shift.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

sd <- lubridate::dhours(6)
napd <- lubridate::dhours(0.5)
nap <- TRUE
sd24(sd, napd, nap)
#> [1] "23400s (~6.5 hours)" # Expected

sd <- lubridate::dhours(9)
napd <- lubridate::dhours(1.5)
nap <- TRUE
sd24(sd, napd, nap)
#> [1] "37800s (~10.5 hours)" # Expected

sd <- lubridate::dhours(6.5)
napd <- lubridate::as.duration(NA)
nap <- FALSE
sd24(sd, napd, nap)
#> [1] "23400s (~6.5 hours)" # Expected

sd <- lubridate::as.duration(NA)
napd <- lubridate::dhours(2.3)
nap <- TRUE
sd24(sd, napd, nap)
#> [1] NA # Expected

## Vector example

sd <- c(lubridate::dhours(7.5), lubridate::dhours(8))
napd <- c(lubridate::dhours(0.75), lubridate::dhours(1))
nap <- c(TRUE, TRUE)
sd24(sd, napd, nap)
#> [1] "29700s (~8.25 hours)" "32400s (~9 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mctq-package.R
\docType{package}
\name{mctq-package}
\alias{mctq}
\alias{mctq-package}
\title{mctq: Tools to Process the Munich ChronoType Questionnaire (MCTQ)}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

A complete and consistent toolkit to process the Munich ChronoType Questionnaire (MCTQ) for its three versions (standard, micro, and shift). MCTQ is a quantitative and validated tool to assess chronotypes using peoples' sleep behavior, originally presented by Till Roenneberg, Anna Wirz-Justice, and Martha Merrow (2003, <doi:10.1177/0748730402239679>).
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/mctq/}
  \item \url{https://github.com/ropensci/mctq/}
  \item Report bugs at \url{https://github.com/ropensci/mctq/issues/}
}

}
\author{
\strong{Maintainer}: Daniel Vartanian \email{danvartan@gmail.com} (\href{https://orcid.org/0000-0001-7782-759X}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Ana Amelia Benedito-Silva \email{aamelia@usp.br} (\href{https://orcid.org/0000-0003-4976-2623}{ORCID}) [scientific advisor]
  \item Mario Pedrazzoli \email{pedrazzo@usp.br} (\href{https://orcid.org/0000-0002-5257-591X}{ORCID}) [scientific advisor]
}

Other contributors:
\itemize{
  \item Jonathan Keane \email{jkeane@gmail.com} (\href{https://orcid.org/0000-0001-7087-9776}{ORCID}) [reviewer]
  \item Mario Andre Leocadio-Miguel \email{miguel.ml.mario@gmail.com} (\href{https://orcid.org/0000-0002-7248-3529}{ORCID}) [reviewer]
  \item Interdisciplinary Sleep Research Group (GIPSO) [funder]
  \item University of Sao Paulo (USP) [funder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_date.R
\name{assign_date}
\alias{assign_date}
\title{Assign dates to two sequential hours}
\usage{
assign_date(start, end, ambiguity = 0)
}
\arguments{
\item{start, end}{A \code{hms} or \code{POSIXt} object indicating the start or end
hour.}

\item{ambiguity}{(optional) a \code{numeric} or \code{NA} value to instruct
\code{assign_date()} on how to deal with ambiguities (see Details) (default:
\code{0}).}
}
\value{
A \code{start}--\code{end} \code{Interval} object.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{assign_date()} assign dates to two sequential hours. It can facilitate
time arithmetic by locating time values without a date reference on a
timeline.
}
\details{
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{\code{ambiguity} argument}{

In cases when \code{start} is equal to \code{end}, there are two possibilities of
intervals between the two hours (ambiguity). That's because \code{start} and \code{end}
can be at the same point in time or they can distance themselves by one day,
considering a two-day timeline.\preformatted{ start,end       start,end       start,end       start,end
   start            end            start            end
   10:00           10:00           10:00           10:00
-----|---------------|---------------|---------------|----->
    0h              0h              0h              0h
            24h             24h             24h
}

You must instruct \code{assign_date()} on how to deal with this problem if it
occurs. There are three options to choose.
\itemize{
\item \code{ambiguity = 0}: to consider the interval between \code{start} and \code{end} as 0
hours, i.e., \code{start} and \code{end} are located at the same point in time
(default).
\item \code{ambiguity = 24}: to consider the interval between \code{start} and \code{end} as 24
hours, i.e., \code{start} and \code{end} distance themselves by one day.
\item \code{ambiguity = NA}: to disregard these cases, assigning \code{NA} as value.
}
}

\subsection{Base date and timezone}{

\code{assign_date()} uses the
\href{https://en.wikipedia.org/wiki/Unix_time}{Unix epoch} (1970-01-01) date as
the start date for creating intervals.

The output will always have \code{"UTC"} set as timezone. Learn more about
time zones in \link[base:timezones]{base::timezone}.
}

\subsection{\code{POSIXt} objects}{

\code{\link[base:as.POSIXlt]{POSIXt}} objects passed as argument to \code{start} or \code{end}
will be stripped of their dates. Only the time will be considered.

Both \code{POSIXct} and \code{POSIXlt} are objects that inherits the class \code{POSIXt}.
Learn more about it in \link[base:DateTimeClasses]{base::DateTimeClasses}.
}

\subsection{\code{NA} values}{

\code{assign_date()} will return an \code{Interval} \code{NA}-\code{NA} if \code{start} or \code{end} are
\code{NA}.
}
}
\examples{
## Scalar example

start <- hms::parse_hms("23:11:00")
end <- hms::parse_hms("05:30:00")
assign_date(start, end)
#> [1] 1970-01-01 23:11:00 UTC--1970-01-02 05:30:00 UTC # Expected

start <- hms::parse_hms("10:15:00")
end <- hms::parse_hms("13:25:00")
assign_date(start, end)
#> [1] 1970-01-01 10:15:00 UTC--1970-01-01 13:25:00 UTC # Expected

start <- hms::parse_hms("05:42:00")
end <- hms::as_hms(NA)
assign_date(start, end)
#> [1] NA--NA # Expected

## Vector example

start <- c(hms::parse_hm("09:45"), hms::parse_hm("20:30"))
end <- c(hms::parse_hm("21:15"), hms::parse_hm("04:30"))
assign_date(start, end)
#> [1] 1970-01-01 09:45:00 UTC--1970-01-01 21:15:00 UTC # Expected
#> [2] 1970-01-01 20:30:00 UTC--1970-01-02 04:30:00 UTC # Expected

## To assign a 24 hours interval to ambiguities

start <- lubridate::as_datetime("1985-01-15 12:00:00")
end <- lubridate::as_datetime("2020-09-10 12:00:00")
assign_date(start, end, ambiguity = 24)
#> [1] 1970-01-01 12:00:00 UTC--1970-01-02 12:00:00 UTC # Expected
}
\seealso{
Other utility functions: 
\code{\link{cycle_time}()},
\code{\link{pretty_mctq}()},
\code{\link{qplot_walk}()},
\code{\link{random_mctq}()},
\code{\link{raw_data}()},
\code{\link{round_time}()},
\code{\link{shorter_interval}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fd.R
\name{fd}
\alias{fd}
\title{Compute MCTQ work-free days}
\usage{
fd(wd)
}
\arguments{
\item{wd}{An \link[checkmate:checkIntegerish]{integerish} \code{numeric} object or
an \code{integer} object corresponding to the \strong{number of workdays per week}
from a standard or micro version of the MCTQ questionnaire.}
}
\value{
An \code{integer} object corresponding to the difference between the
number of days in a week (7) and the number of workdays (\code{wd}).
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{fd()} computes the \strong{number of work-free days per week} for standard
and micro versions of the Munich Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012) and The Worldwide
Experimental Platform (n.d.) guidelines for \code{fd()} (\eqn{FD}) computation are
as follows.

\strong{\deqn{7 - WD}}

Where:
\itemize{
\item \eqn{WD} = number of workdays ("I have a regular work schedule and work ___
days per week").
}
}

\examples{
## Scalar example

fd(5)
#> [1] 2 # Expected
fd(4)
#> [1] 3 # Expected
fd(as.numeric(NA))
#> [1] NA # Expected

## Vector example

fd(0:7)
#> [1] 7 6 5 4 3 2 1 0 # Expected
fd(c(1, NA))
#> [1]  6 NA # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbt.R
\name{tbt}
\alias{tbt}
\title{Compute MCTQ total time in bed}
\usage{
tbt(bt, gu)
}
\arguments{
\item{bt}{A \code{hms} object corresponding to the \strong{local time of going to bed}
from a standard or shift version of the MCTQ questionnaire.}

\item{gu}{A \code{hms} object corresponding to the \strong{local time of getting out of
bed} from a standard or shift version of the MCTQ questionnaire. You can
use \code{\link[=gu]{gu()}} to compute it.}
}
\value{
A \code{Duration} object corresponding to the vectorized difference
between \code{gu} and \code{bt} in a circular time frame of 24 hours.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{tbt()} computes the \strong{total time in bed} for standard and shift versions of
the Munich Chronotype Questionnaire (MCTQ).
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Juda, Vetter, & Roenneberg
(2013), and The Worldwide Experimental Platform (n.d.) guidelines for \code{tbt()}
(\eqn{TBT}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each section of the
questionnaire.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{For standard and micro versions of the MCTQ}{

\strong{\deqn{GU_{W/F} - BT_{W/F}}{GU_W/F - BT_W/F}}

Where:
\itemize{
\item \eqn{BT_{W/F}}{BT_W/F} = local time of going to bed on work \strong{or}
work-free days ("I go to bed at ___ o'clock").
\item \eqn{GU_{W/F}}{GU_W/F} = local time of getting out of bed on work \strong{or}
work-free days.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}

\subsection{For the shift version of the MCTQ}{

\strong{\deqn{GU_{W/F}^{M/E/N} - BT_{W/F}^{M/E/N}}{GU_W/F_M/E/N - BT_W/F_M/E/N}}

Where:
\itemize{
\item \eqn{BT_{W/F}^{M/E/N}}{BT_W/F_M/E/N} = local time of going to bed between
two days in a particular shift \strong{or} between two free days after a
particular shift  ("I go to bed at ___ o'clock").
\item \eqn{GU_{W/F}^{M/E/N}}{GU_W/F_M/E/N} = local time of getting out of bed
between two days in a particular shift \strong{or} between two free days after a
particular shift.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

bt <- hms::parse_hm("22:10")
gu <- hms::parse_hm("06:15")
tbt(bt, gu)
#> [1] "29100s (~8.08 hours)" # Expected

bt <- hms::parse_hm("01:20")
gu <- hms::parse_hm("14:00")
tbt(bt, gu)
#> [1] "45600s (~12.67 hours)" # Expected

bt <- hms::as_hms(NA)
gu <- hms::parse_hm("07:20")
tbt(bt, gu)
#> [1] NA # Expected

## Vector example

bt <- c(hms::parse_hm("23:50"), hms::parse_hm("02:30"))
gu <- c(hms::parse_hm("09:30"), hms::parse_hm("11:25"))
tbt(bt, gu)
#> [1] "34800s (~9.67 hours)" "32100s (~8.92 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sjl.R
\name{sjl}
\alias{sjl}
\alias{sjl_rel}
\title{Compute MCTQ social jetlag}
\usage{
sjl(msw, msf, abs = TRUE, method = "shorter")

sjl_rel(msw, msf, method = "shorter")
}
\arguments{
\item{msw}{A \code{hms} object corresponding to the \strong{local time of mid-sleep on
workdays} from a standard, micro, or shift version of the MCTQ
questionnaire. You can use \code{\link[=ms]{ms()}} to compute it.}

\item{msf}{A \code{hms} object corresponding to the \strong{local time of mid-sleep on
work-free days} from a standard, micro, or shift version of the MCTQ
questionnaire. You can use \code{\link[=msl]{msl()}} to compute it.}

\item{abs}{(optional) a \code{logical} object indicating if the function must
return an absolute social jetlag (default: \code{TRUE}).}

\item{method}{(optional) a string indicating which method the function must
use to compute the social jetlag. See the Methods section to learn
more (default: \code{"shorter"}).}
}
\value{
\itemize{
\item If \code{abs = TRUE}, a \code{Duration} object corresponding to the absolute social
jetlag.
\item If \code{abs = FALSE}, a \code{Duration} object corresponding to the relative social
jetlag.
}

The output may also vary depending on the \code{method} used.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sjl()} computes the \strong{relative or absolute social jetlag} for standard,
micro, and shift versions of the Munich Chronotype Questionnaire (MCTQ).

\code{sjl_rel()} it's just a wrapper for \code{sjl()} with \code{abs = FALSE}.
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Juda, Vetter, & Roenneberg
(2013), and The Worldwide Experimental Platform (n.d.) guidelines for \code{sjl()}
(\eqn{SJL_{rel}}{SJL_rel} and \eqn{SJL}) computation are as follows.
\subsection{Notes}{
\itemize{
\item For MCTQ\eqn{^{Shift}}{ Shift}, the computation below must be applied to
each shift section of the questionnaire.
\item Due to time arithmetic issues, \code{sjl()} does a slightly different
computation by default than those proposed by the authors mentioned above.
See \code{vignette("sjl-computation", package = "mctq")} for more details.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{For standard and micro versions of the MCTQ}{

Relative social jetlag (\eqn{SJL_{rel}}{SJL_rel}):

\strong{\deqn{MSF - MSW}}

Absolute social jetlag (\eqn{SJL}):

\strong{\deqn{| MSF - MSW |}}

Where:
\itemize{
\item \eqn{MSW} = local time of mid-sleep on workdays.
\item \eqn{MSF} = local time of mid-sleep on work-free days.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}

\subsection{For the shift version of the MCTQ}{

Relative social jetlag (\eqn{SJL_{rel}}{SJL_rel}):

\strong{\deqn{MSF^{M/E/N} - MSW^{M/E/N}}{MSF_M/E/N - MSW_M/E/N}}

Absolute social jetlag (\eqn{SJL}):

\strong{\deqn{| MSF^{M/E/N} - MSW^{M/E/N} |}{| MSF_M/E/N - MSW_M/E/N |}}

Where:
\itemize{
\item \eqn{MSW^{M/E/N}}{MSW_M/E/N} = local time of mid-sleep between two days in
a particular shift.
\item \eqn{MSF^{M/E/N}}{MSF_M/E/N} = local time of mid-sleep between two free
days after a particular shift.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\section{Methods for computing the social jetlag}{


There are different approaches to compute the social jetlag (\eqn{SJL}). By
default, \code{sjl()} uses an approach that we call "the shorter interval
approach" (\code{"shorter"}).

The topics below provide a simple explanation of each method supported by
\code{sjl()}. To get a detail understating of this methods, see
\code{vignette("sjl-computation", package = "mctq")}.

Please note that none of the approaches below are related to Jankowski's
(2017) social jetlag sleep-corrected proposal. Since Jankowski's alternative
is still disputed (Roenneberg, Pilz, Zerbini, & Winnebeck, 2019), the \code{mctq}
package currently doesn't provide a function for it. Future versions of the
package may include it.
\itemize{
\item \code{"difference"}
}

By using \code{method = "difference"}, \code{sjl()} will do the exact computation
proposed by the MCTQ authors, i.e., \eqn{SJL} will be computed as the linear
difference between \eqn{MSF} and \eqn{MSW} (see the Guidelines section).

\strong{We do not recommend using this method}, as it has many limitations.
\itemize{
\item \code{"shorter"}
}

This is the default method for \code{sjl()}. It's based on the shorter
interval between \eqn{MSW} and \eqn{MSF}, solving most of the issues
relating to \eqn{SJL} computation.
\itemize{
\item \code{"longer"}
}

The \code{"longer"} method uses the same logic of the \code{"shorter"} method, but,
instead of using the shorter interval between \eqn{MSW} and \eqn{MSF}, it
uses the longer interval between the two, considering a two-day window.

This method may help with special contexts, like when dealing with
shift-workers that have a greater than 12 hours distance between their
mid-sleep hours.
}

\section{Missing sections in standard and micro MCTQ versions}{


Although the standard and micro versions of the MCTQ asks for respondents to
complete the workdays and work-free days sections, even when they do not
have a regular work schedule (\code{wd = 0}) or have a 7 day/week work schedule
(\code{wd = 7}), some of them may still end skipping one of this parts of the
questionnaire. In those cases, \code{sd_week()}, \code{sloss_week()}, \code{le_week()},
\code{msf_sc()}, \code{sjl_rel()}, and \code{sjl()} will produce \code{NA} (Not Available) as
output. That's because those computations combine workdays and work-free days
variables.

For those special standard and micro MCTQ cases, where one section is
missing, a \code{NA} value is the correct output for the functions mentioned above
when \code{wd} (number of workdays per week) are \code{wd > 0 & wd < 7}, but it may not
be when \code{wd == 0} or \code{wd == 7}. There are different approaches to deal with
this issue. See \code{vignette("missing-sections", package = "mctq")} to learn
more.
}

\examples{
## Scalar example

msw <- hms::parse_hm("03:30")
msf <- hms::parse_hm("05:00")
sjl(msw, msf)
#> [1] "5400s (~1.5 hours)" # Expected
sjl(msw, msf, abs = FALSE)
#> [1] "5400s (~1.5 hours)" # Expected
sjl_rel(msw, msf) # Wrapper function
#> [1] "5400s (~1.5 hours)" # Expected

msw <- hms::parse_hm("04:30")
msf <- hms::parse_hm("23:30")
sjl(msw, msf)
#> [1] "18000s (~5 hours)" # Expected
sjl(msw, msf, abs = FALSE)
#> [1] "18000s (~-5 hours)" # Expected
sjl_rel(msw, msf) # Wrapper function
#> [1] "18000s (~-5 hours)" # Expected

msw <- hms::as_hms(NA)
msf <- hms::parse_hm("05:15")
sjl(msw, msf)
#> [1] NA # Expected

## Vector example

msw <- c(hms::parse_hm("02:05"), hms::parse_hm("04:05"))
msf <- c(hms::parse_hm("23:05"), hms::parse_hm("04:05"))
sjl(msw, msf)
#> [1] "10800s (~3 hours)" "0s" # Expected
sjl(msw, msf, abs = FALSE)
#> [1] "-10800s (~-3 hours)" "0s" # Expected
sjl_rel(msw, msf) # Wrapper function
#> [1] "-10800s (~-3 hours)" "0s" # Expected

## Using different methods

msw <- hms::parse_hm("19:15")
msf <- hms::parse_hm("02:30")
sjl(msw, msf, abs = FALSE, method = "difference")
#> [1] "-60300s (~-16.75 hours)" # Expected
sjl(msw, msf, abs = FALSE, method = "shorter") # default method
#> [1] "26100s (~7.25 hours)" # Expected
sjl(msw, msf, abs = FALSE, method = "longer")
#> [1] "-60300s (~-16.75 hours)" # Expected

msw <- hms::parse_hm("02:45")
msf <- hms::parse_hm("04:15")
sjl(msw, msf, abs = FALSE, method = "difference")
#> [1] "5400s (~1.5 hours)" # Expected
sjl(msw, msf, abs = FALSE, method = "shorter") # default method
#> [1] "5400s (~1.5 hours)" # Expected
sjl(msw, msf, abs = FALSE, method = "longer")
#> [1] "-81000s (~-22.5 hours)" # Expected

## Converting the output to `hms`

msw <- hms::parse_hm("01:15")
msf <- hms::parse_hm("03:25")
x <- sjl(msw, msf)
x
#> [1] "7800s (~2.17 hours)" # Expected
hms::as_hms(as.numeric(x))
#> 02:10:00 # Expected

## Rounding the output at the seconds level

msw <- hms::parse_hms("04:19:33.1234")
msf <- hms::parse_hms("02:55:05")
x <- sjl(msw, msf)
x
#> [1] "5068.12339997292s (~1.41 hours)" # Expected
round_time(x)
#> [1] "5068s (~1.41 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Jankowski K. S. (2017). Social jet lag: sleep-corrected formula.
\emph{Chronobiology International}, \emph{34}(4), 531-535.
\doi{10.1080/07420528.2017.1299162}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Pilz, L. K., Zerbini, G., & Winnebeck, E. C. (2019).
Chronotype and social jetlag: a (self-) critical review. \emph{Biology}, \emph{8}(3),
54. \doi{10.3390/biology8030054}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/std_mctq.R
\docType{data}
\name{std_mctq}
\alias{std_mctq}
\title{A fictional standard MCTQ dataset}
\format{
A \code{\link[dplyr:reexports]{tibble}} with 37 columns and 50 rows:

\describe{
\item{id}{
A unique \code{integer} value to identify each respondent in the dataset.
\cr \cr
Type: Control.
\cr \cr
R class: \code{integer}.}

\item{work}{
A \code{logical} value indicating if the respondent has a regular work schedule.
\cr \cr
Statement (\code{EN}): "I have a regular work schedule (this includes being, for
example, a housewife or househusband): Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{wd}{
Number of \strong{workdays} per week.
\cr \cr
Statement (\code{EN}): "I have a regular work schedule and work ___ days per
week".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{integer}.}

\item{fd}{
Number of \strong{work-free days} per week.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{integer}.}

\item{bt_w}{
Local time of going to bed on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "I go to bed at ___ o'clock'".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{sprep_w}{
Local time of preparing to sleep on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "I actually get ready to fall asleep at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{slat_w}{
Sleep latency or time to fall asleep after preparing to sleep on
\strong{workdays}.
\cr \cr
Statement (\code{EN}): "I need ___ minutes to fall asleep".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{so_w}{
Local time of sleep onset on \strong{workdays}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{se_w}{
Local time of sleep end on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "I wake up at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{si_w}{
"Sleep inertia" on \strong{workdays}.
\cr \cr
Despite the name, this variable represents the time the respondent takes to
get up after sleep end.
\cr \cr
Statement (\code{EN}): "After ___ minutes, I get up".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{gu_w}{
Local time of getting out of bed on \strong{workdays}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{alarm_w}{
A \code{logical} value indicating if the respondent uses an alarm clock to wake
up on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "I use an alarm clock on workdays: Yes ( ___ ) No ( ___
)".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{wake_before_w}{
A \code{logical} value indicating if the respondent regularly wakes up
\strong{before} the alarm rings on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "If "Yes": I regularly wake up BEFORE the alarm rings:
Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{sd_w}{
Sleep duration on \strong{workdays}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{tbt_w}{
Total time in bed on \strong{workdays}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{le_w}{
Light exposure on \strong{workdays}.
\cr \cr
Statement (\code{EN}): "On average, I spend the following amount of time
outdoors in daylight (without a roof above my head)".
\cr \cr
Type: Extra.
\cr \cr
R class: \code{Duration}.}

\item{msw}{
Local time of mid-sleep on \strong{workdays}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{bt_f}{
Local time of going to bed on \strong{work-free days}.
\cr \cr
Statement (\code{EN}): "I go to bed at ___ o'clock'".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{sprep_f}{
Local time of preparing to sleep on \strong{work-free days}.
\cr \cr
Statement (\code{EN}): "I actually get ready to fall asleep at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{slat_f}{
Sleep latency or time to fall asleep after preparing to sleep on
\strong{work-free days}.
\cr \cr
Statement (\code{EN}): "I need ___ minutes to fall asleep".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{so_f}{
Local time of sleep onset on \strong{work-free days}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{se_f}{
Local time of sleep end on \strong{work-free days}.
\cr \cr
Statement (\code{EN}): "I wake up at ___ o'clock".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{hms}.}

\item{si_f}{
"Sleep inertia" on \strong{work-free days}.
\cr \cr
Despite the name, this variable represents the time the respondent takes to
get up after sleep end.
\cr \cr
Statement (\code{EN}): "After ___ minutes, I get up".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{Duration}.}

\item{gu_f}{
Local time of getting out of bed on \strong{work-free days}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{alarm_f}{
A \code{logical} value indicating if the respondent uses an alarm clock to wake
up on \strong{work-free days}.
\cr \cr
Statement (\code{EN}): "My wake-up time is due to the use of an alarm
clock: Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{reasons_f}{
A \code{logical} value indicating if the respondent has any particular reasons
for why they \strong{cannot} freely choose their sleep times on \strong{work-free
days}.
\cr \cr
Statement (\code{EN}): "There are particular reasons why I \strong{cannot} freely
choose my sleep times on free days: Yes ( ___ ) No ( ___ )".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{logical}.}

\item{reasons_why_f}{
Particular reasons for why the respondent cannot freely choose their sleep
times on \strong{work-free days}.
\cr \cr
Statement (\code{EN}): "If "Yes": Child(ren)/pet(s) ( ___ ) Hobbies ( ___ )
Others ( ___ ), for example: ___".
\cr \cr
Type: Basic.
\cr \cr
R class: \code{character}.}

\item{sd_f}{
Sleep duration on \strong{work-free days}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{tbt_f}{
Total time in bed on \strong{work-free days}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{le_f}{
Light exposure on \strong{work-free days}.
\cr \cr
Statement (\code{EN}): "On average, I spend the following amount of time
outdoors in daylight (without a roof above my head)".
\cr \cr
Type: Extra.
\cr \cr
R class: \code{Duration}.}

\item{msf}{
Local time of mid-sleep on \strong{work-free days}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{sd_week}{
Average weekly sleep duration.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sloss_week}{
Weekly sleep loss.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{le_week}{
Average weekly light exposure.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{msf_sc}{
Chronotype or corrected local time of mid-sleep on \strong{work-free days}.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{hms}.}

\item{sjl_rel}{
Relative social jetlag.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}

\item{sjl}{
Absolute social jetlag.
\cr \cr
Type: Computed.
\cr \cr
R class: \code{Duration}.}
}
}
\source{
Created by Daniel Vartanian (package author).
}
\usage{
std_mctq
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

A fictional dataset, \strong{for testing and learning purposes}, composed of
basic/measurable and computed variables of the Munich Chronotype
Questionnaire (MCTQ) standard version.

This data was created following the guidelines in Roenneberg, Wirz-Justice, &
Merrow (2003), Roenneberg, Allebrandt, Merrow, & Vetter (2012), and The
Worldwide Experimental Platform (n.d.). See the References and Details
sections to learn more.
}
\details{
\code{std_mctq} is a tidied, validated, and transformed version of
\code{raw_data("std_mctq.csv")}.
\subsection{Guidelines}{

To learn more about the Munich Chronotype Questionnaire (MCTQ),
see Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt,
Merrow, & Vetter (2012), Roenneberg et al. (2015), and Roenneberg, Pilz,
Zerbini, & Winnebeck (2019).

To know about different MCTQ versions, see Juda, Vetter, & Roenneberg
(2013) and Ghotbi et al. (2020).

If you're curious about the variable computations and want to have access to
the full questionnaire, see The Worldwide Experimental Platform (n.d.).
}

\subsection{Data building and data wrangling}{

This dataset was created by randomized sampling (see \code{\link[=random_mctq]{random_mctq()}})
and by manual insertions of special cases. Its purpose is to demonstrate
common cases and data issues that researchers may find in their MCTQ data, in
addition to be a suggested data structure for MCTQ data.

You can see the \code{std_mctq} build and data wrangling processes
\href{https://github.com/ropensci/mctq/blob/main/data-raw/std_mctq.R}{here}.
}

\subsection{Variable naming}{

The naming of the variables took into account the naming scheme used in MCTQ
publications, in addition to the guidelines of the \href{https://style.tidyverse.org/}{tidyverse style guide}.
}

\subsection{Variable classes}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the \link[hms:hms-package]{hms}
and \link[lubridate:lubridate-package]{lubridate} package.
}

\subsection{\code{Duration} objects}{

If you prefer to view \code{Duration} objects as \code{hms} objects, run
\code{pretty_mctq(std_mctq)}.
}
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Keller, L. K., Fischer, D., Matera, J. L., Vetter, C., &
Winnebeck, E. C. (2015). Human activity and rest in situ. In A. Sehgal (Ed.),
\emph{Methods in Enzymology} (Vol. 552, pp. 257-283). London, UK: Academic Press.
\doi{10.1016/bs.mie.2014.11.028}.

Roenneberg, T., Pilz, L. K., Zerbini, G., & Winnebeck, E. C. (2019).
Chronotype and social jetlag: a (self-) critical review. \emph{Biology}, \emph{8}(3),
54. \doi{10.3390/biology8030054}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other datasets: 
\code{\link{micro_mctq}},
\code{\link{shift_mctq}}
}
\concept{datasets}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msl.R
\name{msf_sc}
\alias{msf_sc}
\alias{chronotype}
\title{Compute MCTQ corrected local time of mid-sleep on work-free days}
\usage{
msf_sc(msf, sd_w, sd_f, sd_week, alarm_f)

chronotype(msf, sd_w, sd_f, sd_week, alarm_f)
}
\arguments{
\item{msf}{A \code{hms} object corresponding to the \strong{local time of mid-sleep on
work-free days} from a standard, micro, or shift version of the MCTQ
questionnaire. You can use \code{\link[=msl]{msl()}} to compute it.}

\item{sd_w}{A \code{Duration} object corresponding to the \strong{sleep duration on work
days} from a standard, micro, or shift version of the MCTQ questionnaire.
You can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{sd_f}{A \code{Duration} object corresponding to the \strong{sleep duration on
work-free days} from a standard, micro, or shift version of the MCTQ
questionnaire. You can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{sd_week}{A \code{Duration} object corresponding to the \strong{average weekly
sleep duration} from a standard or micro version of the MCTQ questionnaire
(you can use \code{\link[=sd_week]{sd_week()}} to compute it) \strong{or} the \strong{overall sleep
duration of a particular shift} from a shift version of the MCTQ
questionnaire (you can use \code{\link[=sd_overall]{sd_overall()}} to compute it).}

\item{alarm_f}{A \code{logical} object corresponding to the \strong{alarm clock use on
work-free days} from a standard, micro, or shift version of the MCTQ
questionnaire. Note that, if \code{alarm_f == TRUE}, \code{msf_sc} cannot be
computed, \code{msf_sc()} will return \code{NA} for those cases. For the
\eqn{\mu}MCTQ, this value must be set as \code{FALSE} all times, since the
questionnaire considers only the work-free days when the respondent does
not use an alarm.}
}
\value{
A \code{hms} object corresponding to the MCTQ chronotype or corrected
local time of mid-sleep on work-free days.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{msf_sc()} computes the \strong{chronotype or corrected local time of mid-sleep on
work-free days} for standard, micro, and shift versions of the Munich
Chronotype Questionnaire (MCTQ).

\code{chronotype()} is just a wrapper for \code{msf_sc()}.

When using the shift version of the MCTQ, replace the value of \code{sd_week} to
\code{sd_overall}, as instructed in the Arguments section.
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Ghotbi et al. (2020), Juda,
Vetter, & Roenneberg (2013), and The Worldwide Experimental Platform (n.d.)
guidelines for \code{msf_sc()} (\eqn{MSF_{sc}}{MSF_sc}) computation are as
follows.
\subsection{Notes}{
\itemize{
\item For all cases, \eqn{MSF_{sc}}{MSF_sc} cannot be computed if the participant
wakes up with an alarm clock on work-free days (\eqn{Alarm_F}{Alarm_F}).
\item For MCTQ\eqn{^{Shift}}{ Shift}, the computation below must be applied to
each shift section of the questionnaire.
\item \eqn{MSF_{sc}}{MSF_sc} is a proxy for the participant chronotype in
standard and micro versions of the MCTQ.
\item The basis for estimating chronotype in shift-workers is the mid-sleep on
work-free days after evening shifts (\eqn{MSF^E}{MSF_E}). In case work
schedules do not comprise evening shifts, Juda, Vetter, & Roenneberg (2013)
propose to derive it from the \eqn{MSF_{sc}}{MSF_sc} of other shifts (e.g.,
by using a linear model). Unfortunately, the \code{mctq} package can't help you
with that, as it requires a closer look at your data.
\item \eqn{MSF_{sc}}{MSF_sc} depends on developmental and environmental
conditions (e.g., age, light exposure). For epidemiological and genetic
studies, \eqn{MSF_{sc}}{MSF_sc} must be normalized for age and sex to make
populations of different age and sex compositions comparable (Roenneberg,
Allebrandt, Merrow, & Vetter, 2012).
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{For standard and micro versions of the MCTQ}{

\strong{\deqn{\textrm{If } SD_F \leq SD_W \; , \; MSF}{If SD_F <= SD_W, MSF}}
\strong{\deqn{\textrm{If } SD_F > SD_W \; , \; MSF - \frac{SD_F - SD_{week}}{2}}{
If SD_F > SD_W, MSF - (SD_F - SD_week) / 2}}

Where:
\itemize{
\item \eqn{MSF} = local time of mid-sleep on work-free days.
\item \eqn{SD_W} = sleep duration on workdays.
\item \eqn{SD_F} = sleep duration on work-free days.
\item \eqn{SD_{week}}{SD_week} = average weekly sleep duration.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}

\subsection{For the shift version of the MCTQ}{

\strong{\deqn{\textrm{If } SD_{F}^{M/E/N} \leq SD_{W}^{M/E/N} \; , \; MSF^{M/E/N}}{
If SD_F_M/E/N <= SD_W_M/E/N, MSF}}
\strong{\deqn{\textrm{If } SD_{F}^{M/E/N} > SD_{W}^{M/E/N} \; , \; MSF^{M/E/N} -
\frac{SD_{F}^{M/E/N} - \emptyset SD^{M/E/N}}{2}}{If SD_F_M/E/N > SD_W_M/E/N,
MSF - (SD_F_M/E/N - OSD_M/E/N) / 2}}

Where:
\itemize{
\item \eqn{MSF^{M/E/N}}{MSF_M/E/N} = local time of mid-sleep between two free
days after a particular shift.
\item \eqn{SD_{W}^{M/E/N}}{SD_W_M/E/N} = sleep duration between two days in a
particular shift.
\item \eqn{SD_{F}^{M/E/N}}{SD_F_M/E/N} = sleep duration between two free days
after a particular shift.
\item \eqn{\emptyset SD^{M/E/N}}{OSD_M/E/N} = overall sleep duration of a
particular shift.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\section{Missing sections in standard and micro MCTQ versions}{


Although the standard and micro versions of the MCTQ asks for respondents to
complete the workdays and work-free days sections, even when they do not
have a regular work schedule (\code{wd = 0}) or have a 7 day/week work schedule
(\code{wd = 7}), some of them may still end skipping one of this parts of the
questionnaire. In those cases, \code{sd_week()}, \code{sloss_week()}, \code{le_week()},
\code{msf_sc()}, \code{sjl_rel()}, and \code{sjl()} will produce \code{NA} (Not Available) as
output. That's because those computations combine workdays and work-free days
variables.

For those special standard and micro MCTQ cases, where one section is
missing, a \code{NA} value is the correct output for the functions mentioned above
when \code{wd} (number of workdays per week) are \code{wd > 0 & wd < 7}, but it may not
be when \code{wd == 0} or \code{wd == 7}. There are different approaches to deal with
this issue. See \code{vignette("missing-sections", package = "mctq")} to learn
more.
}

\examples{
## Scalar example

msf <- hms::parse_hms("04:00:00")
sd_w <- lubridate::dhours(6)
sd_f <- lubridate::dhours(7)
sd_week <- lubridate::dhours(6.29)
alarm_f <- FALSE
msf_sc(msf, sd_w, sd_f, sd_week, alarm_f)
#> 03:38:42 # Expected

msf <- hms::parse_hm("01:00:00")
sd_w <- lubridate::dhours(5.5)
sd_f <- lubridate::dhours(9)
sd_week <- lubridate::dhours(6.75)
alarm_f <- FALSE
msf_sc(msf, sd_w, sd_f, sd_week, alarm_f)
#> 23:52:30 # Expected

msf <- hms::parse_hms("05:40:00")
sd_w <- lubridate::dhours(7.5)
sd_f <- lubridate::dhours(10)
sd_week <- lubridate::dhours(8.5)
alarm_f <- TRUE
msf_sc(msf, sd_w, sd_f, sd_week, alarm_f)
#> NA # Expected (`msf_sc` cannot be computed if `alarm_f == TRUE`)

## Vector example

msf <- c(hms::parse_hms("03:45:00"), hms::parse_hm("04:45:00"))
sd_w <- c(lubridate::dhours(9), lubridate::dhours(6.45))
sd_f <- c(lubridate::dhours(5), lubridate::dhours(10))
sd_week <- c(lubridate::dhours(8.5), lubridate::dhours(9.2))
alarm_f <- c(FALSE, FALSE)
msf_sc(msf, sd_w, sd_f, sd_week, alarm_f)
#> 03:45:00 # Expected
#> 04:21:00 # Expected

## A wrapper for msf_sc()

msf <- hms::parse_hms("07:00:00")
sd_w <- lubridate::dhours(6)
sd_f <- lubridate::dhours(12)
sd_week <- lubridate::dhours(9.45)
alarm_f <- FALSE
chronotype(msf, sd_w, sd_f, sd_week, alarm_f)
#> 05:43:30 # Expected

## Rounding the output at the seconds level

msf <- hms::parse_hms("05:40:00")
sd_w <- lubridate::dhours(5.43678)
sd_f <- lubridate::dhours(9.345111)
sd_week <- lubridate::dhours(7.5453)
alarm_f <- FALSE
x <- msf_sc(msf, sd_w, sd_f, sd_week, alarm_f)
x
#> 04:46:00.3402 # Expected
round_time(x)
#> 04:46:00 # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/round_time.R
\name{round_time}
\alias{round_time}
\alias{round_time.Duration}
\alias{round_time.difftime}
\alias{round_time.hms}
\alias{round_time.POSIXct}
\alias{round_time.POSIXlt}
\title{Round time values}
\usage{
round_time(x)

\method{round_time}{Duration}(x)

\method{round_time}{difftime}(x)

\method{round_time}{hms}(x)

\method{round_time}{POSIXct}(x)

\method{round_time}{POSIXlt}(x)
}
\arguments{
\item{x}{An object belonging to one of the following classes: \code{Duration},
\code{Period}, \code{difftime}, \code{hms}, \code{POSIXct}, or \code{POSIXlt}.}
}
\value{
An object of the same class of \code{x} rounded at the level of seconds.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{round_time()} takes a \code{Duration}, \code{difftime}, \code{hms}, \code{POSIXct}, or
\code{POSIXlt} object and round it at the level of seconds.
}
\details{
\subsection{Round standard}{

\code{round_time()} uses \code{\link[base:Round]{base::round()}} for rounding. That is to say that
\code{round_time()} uses the same IEC 60559 standard (\emph{"go to the even digit"})
for rounding off a 5. Therefore, \code{round(0.5)} is equal to 0 and \code{round(-1.5)}
is equal to -2. See \code{?round} to learn more.
}

\subsection{\code{Period} objects}{

\code{\link[lubridate:period]{Period}} objects are a special type of object
developed by the \link[lubridate:lubridate-package]{lubridate} team that
represents "human units", ignoring possible timeline irregularities. That is
to say that 1 day as \code{Period} can have different time spans, when looking to
a timeline after a irregularity event.

Since the time span of a \code{Period} object can fluctuate, \code{round_time()} don't
accept this kind of object. You can transform it to a \code{Duration} object and
still use the function, but beware that this can produce errors.

Learn more about \code{Period} objects in the \href{https://r4ds.had.co.nz/dates-and-times.html#periods}{Dates and times} chapter of
Wickham & Grolemund (n.d.).
}
}
\examples{
## Scalar example

lubridate::dmilliseconds(123456789)
#> [1] "123456.789s (~1.43 days)" # Expected
round_time(lubridate::dmilliseconds(123456789))
#> [1] "123457s (~1.43 days)" # Expected

as.difftime(12345.6789, units = "secs")
#> Time difference of 12345.68 secs # Expected
round_time(as.difftime(12345.6789, units = "secs"))
#> Time difference of 12346 secs # Expected

hms::as_hms(12345.6789)
#> 03:25:45.6789 # Expected
round_time(hms::as_hms(12345.6789))
#> 03:25:46 # Expected

lubridate::as_datetime(12345.6789, tz = "EST")
#> [1] "1969-12-31 22:25:45 EST" # Expected
as.numeric(lubridate::as_datetime(12345.6789, tz = "EST"))
#> [1] 12345.68 # Expected
round_time(lubridate::as_datetime(12345.6789, tz = "EST"))
#> [1] "1969-12-31 22:25:46 EST" # Expected
as.numeric(round_time(lubridate::as_datetime(12345.6789, tz = "EST")))
#> [1] 12346 # Expected

## Vector example

x <- c(lubridate::dhours(5.6987), lubridate::dhours(2.6875154))
x
#> [1] "20515.32s (~5.7 hours)"    "9675.05544s (~2.69 hours)" # Expected
round_time(x)
#> [1] "20515s (~5.7 hours)" "9675s (~2.69 hours)" # Expected
}
\references{
Wickham, H., & Grolemund, G. (n.d.). \emph{R for data science}. Sebastopol, CA:
O'Reilly Media. \url{https://r4ds.had.co.nz}
}
\seealso{
Other date-time rounding functions: \code{\link[hms:round_hms]{hms::round_hms()}}
\code{\link[hms:round_hms]{hms::trunc_hms()}} \code{\link[lubridate:round_date]{lubridate::round_date()}}.

Other utility functions: 
\code{\link{assign_date}()},
\code{\link{cycle_time}()},
\code{\link{pretty_mctq}()},
\code{\link{qplot_walk}()},
\code{\link{random_mctq}()},
\code{\link{raw_data}()},
\code{\link{shorter_interval}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdu.R
\name{sd_week}
\alias{sd_week}
\title{Compute MCTQ average weekly sleep duration}
\usage{
sd_week(sd_w, sd_f, wd)
}
\arguments{
\item{sd_w}{A \code{Duration} object corresponding to the \strong{sleep duration on
workdays} from a standard or micro version of the MCTQ questionnaire. You
can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{sd_f}{A \code{Duration} object corresponding to the \strong{sleep duration on
work-free days} from a standard or micro version of the MCTQ
questionnaire. You can use \code{\link[=sdu]{sdu()}} to compute it.}

\item{wd}{An \link[checkmate:checkIntegerish]{integerish} \code{numeric} object or
an \code{integer} object corresponding to the \strong{number of workdays per week}
from a standard or micro version of the MCTQ questionnaire.}
}
\value{
A \code{Duration} object corresponding to the vectorized weighted mean of
\code{sd_w} and \code{sd_f} with \code{wd} and \code{fd(wd)} as weights.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sd_week()} computes the \strong{average weekly sleep duration} for the standard
and micro versions of the Munich Chronotype Questionnaire (MCTQ).

See \code{\link[=sd_overall]{sd_overall()}} to compute the overall sleep duration of a
particular shift for the shift version of the MCTQ.
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Ghotbi et al. (2020), and
The Worldwide Experimental Platform (n.d.) guidelines for \code{sd_week()}
(\eqn{SD_{week}}{SD_week}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The average weekly sleep duration is the weighted average of the sleep
durations on work and work-free days in a week.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{Computation}{

\strong{\deqn{\frac{SD_W \times WD + SD_F \times FD}{7}}{
(SD_W * WD + SD_F * FD) / 7}}

Where:
\itemize{
\item \eqn{SD_W} = sleep duration on workdays.
\item \eqn{SD_F} = sleep duration on work-free days.
\item \eqn{WD} = number of workdays per week ("I have a regular work schedule and
work ___ days per week").
\item \eqn{FD} = number of work-free days per week.
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}
}

\section{Missing sections in standard and micro MCTQ versions}{


Although the standard and micro versions of the MCTQ asks for respondents to
complete the workdays and work-free days sections, even when they do not
have a regular work schedule (\code{wd = 0}) or have a 7 day/week work schedule
(\code{wd = 7}), some of them may still end skipping one of this parts of the
questionnaire. In those cases, \code{sd_week()}, \code{sloss_week()}, \code{le_week()},
\code{msf_sc()}, \code{sjl_rel()}, and \code{sjl()} will produce \code{NA} (Not Available) as
output. That's because those computations combine workdays and work-free days
variables.

For those special standard and micro MCTQ cases, where one section is
missing, a \code{NA} value is the correct output for the functions mentioned above
when \code{wd} (number of workdays per week) are \code{wd > 0 & wd < 7}, but it may not
be when \code{wd == 0} or \code{wd == 7}. There are different approaches to deal with
this issue. See \code{vignette("missing-sections", package = "mctq")} to learn
more.
}

\examples{
## Scalar example

sd_w <- lubridate::dhours(4)
sd_f <- lubridate::dhours(8)
wd <- 5
sd_week(sd_w, sd_f, wd)
#> [1] "18514.2857142857s (~5.14 hours)" # Expected

sd_w <- lubridate::dhours(7)
sd_f <- lubridate::dhours(7)
wd <- 4
sd_week(sd_w, sd_f, wd)
#> [1] "25200s (~7 hours)" # Expected

sd_w <- lubridate::as.duration(NA)
sd_f <- lubridate::dhours(10)
wd <- 6
sd_week(sd_w, sd_f, wd)
#> [1] NA # Expected

## Vector example

sd_w <- c(lubridate::dhours(4.5), lubridate::dhours(5.45))
sd_f <- c(lubridate::dhours(8), lubridate::dhours(7.3))
wd <- c(3, 7)
sd_week(sd_w, sd_f, wd)
#> [1] "23400s (~6.5 hours)"  "19620s (~5.45 hours)" # Expected

## Checking second output from vector example

if (requireNamespace("stats", quietly = TRUE)) {
    i <- 2
    x <- c(sd_w[i], sd_f[i])
    w <- c(wd[i], fd(wd[i]))
    lubridate::as.duration(stats::weighted.mean(x, w))
}
#> [1] "19620s (~5.45 hours)" # Expected

## Converting the output to `hms`

sd_w <- lubridate::dhours(5.45)
sd_f <- lubridate::dhours(9.5)
wd <- 5
x <- sd_week(sd_w, sd_f, wd)
x
#> [1] "23785.7142857143s (~6.61 hours)" # Expected
hms::as_hms(as.numeric(x))
#> 06:36:25.714286 # Expected

## Rounding the output at the seconds level

sd_w <- lubridate::dhours(4.5)
sd_f <- lubridate::dhours(7.8)
wd <- 3
x <- sd_week(sd_w, sd_f, wd)
x
#> [1] "22988.5714285714s (~6.39 hours)" # Expected
round_time(x)
#> [1] "22989s (~6.39 hours)" # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{so}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/so.R
\name{so}
\alias{so}
\title{Compute MCTQ local time of sleep onset}
\usage{
so(sprep, slat)
}
\arguments{
\item{sprep}{A \code{hms} object corresponding to the \strong{local time of preparing to
sleep} from a standard or shift version of the MCTQ questionnaire.}

\item{slat}{A \code{Duration} object corresponding to the \strong{sleep latency or time
to fall asleep after preparing to sleep} from a standard or shift version
of the MCTQ questionnaire.}
}
\value{
A \code{hms} object corresponding to the vectorized sum of \code{sprep} and
\code{slat} in a circular time frame of 24 hours.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{so()} computes the \strong{local time of sleep onset} for standard and shift
versions of the Munich Chronotype Questionnaire (MCTQ).

Note that this value is collected directly from the questionnaire if you're
using the \eqn{\mu}MCTQ.
}
\details{
\strong{Standard MCTQ} functions were created following the guidelines in
Roenneberg, Wirz-Justice, & Merrow (2003), Roenneberg, Allebrandt, Merrow, &
Vetter (2012), and from The Worldwide Experimental Platform (theWeP, n.d.).

\strong{\eqn{\mu}MCTQ} functions were created following the guidelines in Ghotbi
et al. (2020), in addition to the guidelines used for the standard MCTQ.

\strong{MCTQ\eqn{^{Shift}}{ Shift}} functions were created following the
guidelines in Juda, Vetter, & Roenneberg (2013), in addition to the
guidelines used for the standard MCTQ.

See the References section to learn more.
\subsection{Class requirements}{

The \code{mctq} package works with a set of object classes specially created to
hold time values. These classes can be found in the
\link[lubridate:lubridate-package]{lubridate} and \link[hms:hms-package]{hms}
packages. Please refer to those package documentations to learn more about
them.
}

\subsection{Rounding and fractional time}{

Some operations may produce an output with fractional time (e.g.,
\code{"19538.3828571429s (~5.43 hours)"}, \code{01:15:44.505}). If you want, you
can round it with \link[=round_time]{round_time()}.

Our recommendation is to avoid rounding, but, if you do, make sure that you
only round your values after all computations are done. That way you avoid
\href{https://en.wikipedia.org/wiki/Round-off_error}{round-off errors}.
}
}
\section{Guidelines}{


Roenneberg, Allebrandt, Merrow, & Vetter (2012), Juda, Vetter, & Roenneberg
(2013), and The Worldwide Experimental Platform (n.d.) guidelines for \code{so()}
(\eqn{SO}) computation are as follows.
\subsection{Notes}{
\itemize{
\item The computation below must be applied to each section of the
questionnaire.
\item If you are visualizing this documentation in plain text (\code{ASCII}), you may
have some trouble understanding the equations. If you want a better viewer,
you can see this documentation on the package
\href{https://docs.ropensci.org/mctq/reference/}{website}.
}
}

\subsection{For standard and micro versions of the MCTQ}{

\strong{\deqn{SPrep_{W/F} + SLat_{W/F}}{SPrep_W/F + SLat_W/F}}

Where:
\itemize{
\item \eqn{SPrep_{W/F}}{SPrep_W/F} = local time of preparing to sleep on work
\strong{or} work-free days ("I actually get ready to fall asleep at ___ o'clock").
\item \eqn{SLat_{W/F}}{SLat_W/F} = sleep latency or time to fall asleep after
preparing to sleep on work \strong{or} work-free days ("I need ___ min to fall
asleep").
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days.
}

\subsection{For the shift version of the MCTQ}{

\strong{\deqn{SPrep_{W/F}^{M/E/N} + SLat_{W/F}^{M/E/N}}{SPrep_W/F_M/E/N +
SLat_W/F_M/E/N}}

Where:
\itemize{
\item \eqn{SPrep_{W/F}^{M/E/N}}{SPrep_W/F_M/E/N} = local time of preparing to
sleep between two days in a particular shift \strong{or} between two free days
after a particular shift ("I actually get ready to fall asleep at ___
o'clock").
\item \eqn{SLat_{W/F}^{M/E/N}}{SLat_W/F_M/E/N} = sleep latency or time to fall
asleep after preparing to sleep between two days in a particular shift \strong{or}
between two free days after a particular shift ("I need ___ min to fall
asleep").
}

\strong{*} \eqn{W} = workdays; \eqn{F} = work-free days, \eqn{M} =
morning shift; \eqn{E} = evening shift; \eqn{N} = night shift.
}
}

\examples{
## Scalar example

sprep <- hms::parse_hm("22:00")
slat <- lubridate::dminutes(15)
so(sprep, slat)
#> 22:15:00 # Expected

sprep <- hms::parse_hm("23:30")
slat <- lubridate::dminutes(45)
so(sprep, slat)
#> 00:15:00 # Expected

sprep <- hms::parse_hm("20:45")
slat <- lubridate::as.duration(NA)
so(sprep, slat)
#> NA # Expected

## Vector example

sprep <- c(hms::parse_hm("21:30"), hms::parse_hm("22:15"))
slat <- c(lubridate::dminutes(45), lubridate::dminutes(5))
so(sprep, slat)
#> 22:15:00 # Expected
#> 22:20:00 # Expected
}
\references{
Ghotbi, N., Pilz, L. K., Winnebeck, E. C., Vetter, C., Zerbini, G., Lenssen,
D., Frighetto, G., Salamanca, M., Costa, R., Montagnese, S., & Roenneberg, T.
(2020). The \eqn{\mu}MCTQ: an ultra-short version of the Munich ChronoType
Questionnaire. \emph{Journal of Biological Rhythms}, \emph{35}(1), 98-110.
\doi{10.1177/0748730419886986}.

Juda, M., Vetter, C., & Roenneberg, T. (2013). The Munich ChronoType
Questionnaire for shift-workers (MCTQ\eqn{^{Shift}}{ Shift}). \emph{Journal of
Biological Rhythms}, \emph{28}(2), 130-140. \doi{10.1177/0748730412475041}.

Roenneberg T., Allebrandt K. V., Merrow M., & Vetter C. (2012). Social jetlag
and obesity. \emph{Current Biology}, \emph{22}(10), 939-43.
\doi{10.1016/j.cub.2012.03.038}.

Roenneberg, T., Wirz-Justice, A., & Merrow, M. (2003). Life between clocks:
daily temporal patterns of human chronotypes. \emph{Journal of Biological
Rhythms}, \emph{18}(1), 80-90. \doi{10.1177/0748730402239679}.

The Worldwide Experimental Platform (n.d.). MCTQ.
\url{https://www.thewep.org/documentations/mctq/}
}
\seealso{
Other MCTQ functions: 
\code{\link{fd}()},
\code{\link{gu}()},
\code{\link{le_week}()},
\code{\link{msf_sc}()},
\code{\link{msl}()},
\code{\link{napd}()},
\code{\link{sd24}()},
\code{\link{sd_overall}()},
\code{\link{sd_week}()},
\code{\link{sdu}()},
\code{\link{sjl_weighted}()},
\code{\link{sjl}()},
\code{\link{tbt}()}
}
\concept{MCTQ functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sum_time.R
\name{sum_time}
\alias{sum_time}
\alias{vct_sum_time}
\title{Sum time objects}
\usage{
sum_time(..., cycle = NULL, reverse = TRUE, na_rm = FALSE)

vct_sum_time(..., cycle = NULL, reverse = TRUE, na_rm = FALSE)
}
\arguments{
\item{...}{Objects belonging to one of the following classes: \code{Duration},
\code{difftime}, \code{hms}, \code{POSIXct}, \code{POSIXlt}, or \code{Interval}.}

\item{cycle}{(optional) A \code{numeric} or \code{Duration} object of length 1, equal
or greater than 0, indicating the cycle length in seconds. If \code{NULL} the
function will perform a linear sum (see Details to learn more) (default:
\code{NULL}).}

\item{reverse}{(optional) A \code{logical} value indicating if the function must
use a reverse cycle for negative sums (see Details to learn more) (default:
\code{TRUE}).}

\item{na_rm}{(optional) a \code{logical} value indicating if the function must
remove \code{NA} values while performing the sum (default: \code{FALSE}).}
}
\value{
\itemize{
\item If \code{cycle = NULL}, a \code{Duration} object with a linear sum of the time from
objects in \code{...}.
\item If \code{cycle != NULL}, a \code{Duration} object with a circular sum of the time
from objects in \code{...}.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{sum_time()} returns the sum of the time from different kinds of date/time
objects.

\code{vct_sum_time()} returns the vectorized sum of the time from different kinds
of date/time objects.

Both functions can be set to work with a circular time frame (see Details to
learn more).
}
\details{
\subsection{\code{sum_time()} versus \code{vct_sum_time()}}{

\code{sum_time()} behaves similar to \code{\link[base:sum]{base::sum()}}, in the sense that it
aggregates the time lengths of values in \code{...} into a single data point. For
example, \code{sum_time(c(x, y), z)} will have the same output as \code{sum_time(x, y, z)}.

\code{vct_sum_time()} performs a different type of sum (a vectorized one). Instead
of aggregating the time lengths, the function perform a paired sum between
elements. For example, \code{sum_time(c(x, y), c(w, z))} will return a vector like
\code{c(sum_time(x, w), sum_time(y, z))}. Because of that, \code{vct_sum_time()}
requires that all objects in \code{...} have the same length.
}

\subsection{Linear versus circular time}{

Time can have different "shapes".

If the objective is to measure the duration (time span) of an event, time is
usually measured considering a linear frame, with a fixed point of
\href{https://en.wikipedia.org/wiki/Origin_(mathematics)}{origin}. In this
context, the time value distance itself to infinity in relation to the
origin.\preformatted{                                   B
                             |----------|
                                        A
                             |---------------------|
 - inf                                                inf +
<----------------------------|----------|----------|------->
 s                           0          5          10     s
                           origin

A + B = 10 + 5 = 15s
}

But that's not the only possible "shape" of time, as it can also be measured
in other contexts.

In a "time of day" context, time will be linked to the rotation of the
earth, "resetting" when a new rotation cycle starts. That brings a different
kind of shape to time: a circular shape. With this shape the time value
encounters the origin at the end of each cycle.\preformatted{               - <--- h ---> +
                    origin
                . . . 0 . . .
             .                 .
            .                   .
           .                     .
          .                       .
         .                         .
         18                        6
         .                         .
          .                       .
           .                     .
            .                   .
             .                 .
                . . . 12 . . .

18 + 6 = 0h
}

If we transpose this circular time frame to a linear one, it would look like
this:\preformatted{<----|---------------|---------------|---------------|----->
    0h              12h              0h             12h
  origin                           origin
}

Note that now the origin is not fix, but cyclical.

\code{sum_time()} and \code{vct_sum_time()} can both operate in either a linear or a
circular fashion. If \code{cycle = NULL} (default), the function will use a
linear approach. Else, the function will use a circular approach relative to
the cycle length (e.g, \code{cycle = 86400} (1 day)).
}

\subsection{Fractional time}{

\code{sum_time()} uses the \code{\%\%} operator to cycle values. Hence, it can be subject
to catastrophic loss of accuracy if values in \code{...} are fractional and much
larger than \code{cycle}. A warning is given if this is detected.

\code{\%\%} is a \code{builtin} R function that operates like this:\preformatted{function(a, b) \{
    a - floor(a / b) * b
\}
}
}

\subsection{Negative time cycling}{

If the sum of the time is negative, with a \code{cycle} assigned and
\code{reverse = FALSE}, \code{sum_time()} and \code{vtc_sum_time()} will perform the cycle
considering the absolute value of the sum and return the result with a
negative signal.

However, If the sum of the time have a negative value, with a \code{cycle}
assigned and \code{reverse = TRUE} (default), \code{sum_time()} and \code{vtc_sum_time()}
will perform the cycle in reverse, relative to its origin.

Example: If the sum of the time have a -30h time span in a reversed cycle of
24h, the result will be 18h. By removing the full cycles of -30h you will
get -6h (-30 + 24), and -6h relative to the origin will be 18h.\preformatted{               - <--- h ---> +
                    origin
                . . . 0 . . .
              .                 .
            .                   .
           .                     .
          .                       .
         .                         .
    (-6) 18                        6 (-18)
         .                         .
          .                       .
           .                     .
            .                   .
             .                 .
                . . . 12 . . .
                    (-12)
}
}

\subsection{\code{Period} objects}{

\code{\link[lubridate:period]{Period}} objects are a special type of object
developed by the \link[lubridate:lubridate-package]{lubridate} team that
represents "human units", ignoring possible timeline irregularities. That is
to say that 1 day as \code{Period} can have different time spans, when looking to
a timeline after a irregularity event.

Since the time span of a \code{Period} object can fluctuate, \code{sum_time()}
and \code{vct_sum_time()} don't accept this kind of object. You can transform
it to a \code{Duration} object and still use the functions, but beware that
this can produce errors.

Learn more about \code{Period} objects in the \href{https://r4ds.had.co.nz/dates-and-times.html#periods}{Dates and times} chapter of
Wickham & Grolemund (n.d.).
}

\subsection{\code{POSIXt} objects}{

\code{\link[base:as.POSIXlt]{POSIXt}} objects in \code{...} will be stripped of their
dates. Only the time will be considered.

Both \code{POSIXct} and \code{POSIXlt} are objects that inherits the class \code{POSIXt}.
Learn more about it in \link[base:DateTimeClasses]{base::DateTimeClasses}.
}

\subsection{\code{Interval} objects}{

By using \code{\link[lubridate:interval]{Interval}} objects in \code{...}, \code{sum_time()}
and \code{vct_sum_time()} will consider only their time spans. That is, the
amount of seconds of the intervals.

Learn more about \code{Interval} objects in the \href{https://r4ds.had.co.nz/dates-and-times.html#periods}{Dates and times} chapter of
Wickham & Grolemund (n.d.).
}

\subsection{Timeline irregularities}{

This function does not take into account timeline irregularities (e.g.,
leap years, DST, leap seconds). This may not be an issue for most people, but
it must be considered when doing time arithmetic.
}
}
\examples{
## Non-vectorized sum in an linear time frame

x <- c(as.POSIXct("2020-01-01 15:00:00"), as.POSIXct("1999-05-04 17:30:00"))
y <- lubridate::as.interval(lubridate::dhours(7), as.Date("1970-05-08"))
sum_time(x, y)
#> [1] "142200s (~1.65 days)" # 39:30:00 # Expected

## Non-vectorized sum in a circular time frame of 24 hours

x <- c(lubridate::dhours(25), lubridate::dhours(5), lubridate::dminutes(50))
sum_time(x, cycle = lubridate::ddays())
#> [1] "24600s (~6.83 hours)" # 06:50:00 # Expected

x <- c(hms::parse_hm("00:15"), hms::parse_hm("02:30"), hms::as_hms(NA))
sum_time(x, cycle = lubridate::ddays())
#> NA # Expected
sum_time(x, cycle = lubridate::ddays(), na_rm = TRUE)
#> [1] "9900s (~2.75 hours)" # 02:45:00 # Expected

x <- c(lubridate::dhours(-12), lubridate::dhours(-13))
sum_time(x, cycle = lubridate::ddays(), reverse = FALSE)
#> [1] "-3600s (~-1 hours)" # -01:00:00 # Expected

x <- c(lubridate::dhours(-12), lubridate::dhours(-13))
sum_time(x, cycle = lubridate::ddays(), reverse = TRUE)
#> [1] "82800s (~23 hours)" # 23:00:00 # Expected

## Vectorized sum in an linear time frame

x <- c(lubridate::dhours(6), NA)
y <- c(hms::parse_hm("23:00"), hms::parse_hm("10:00"))
vct_sum_time(x, y)
#> [1] "104400s (~1.21 days)" NA # 29:00:00 NA # Expected
vct_sum_time(x, y, na_rm = TRUE)
#> [1] "104400s (~1.21 days)" "36000s (~10 hours)" # Expected

## Vectorized sum in a circular time frame of 24 hours

x <- c(lubridate::dhours(6), NA)
y <- c(hms::parse_hm("23:00"), hms::parse_hm("10:00"))
vct_sum_time(x, y, cycle = lubridate::ddays())
#> [1] "18000s (~5 hours)" NA  # Expected
vct_sum_time(x, y, cycle = lubridate::ddays(), na_rm = TRUE)
#> [1] "18000s (~5 hours)"  "36000s (~10 hours)" # Expected

x <- c(lubridate::dhours(-49), lubridate::dhours(-24))
y <- c(hms::parse_hm("24:00"), - hms::parse_hm("06:00"))
vct_sum_time(x, y, cycle = lubridate::ddays(), reverse = FALSE)
#> [1] "-3600s (~-1 hours)"  "-21600s (~-6 hours)" # Expected

x <- c(lubridate::dhours(-49), lubridate::dhours(-24))
y <- c(hms::parse_hm("24:00"), - hms::parse_hm("06:00"))
vct_sum_time(x, y, cycle = lubridate::ddays(), reverse = TRUE)
#> [1] "82800s (~23 hours)" "64800s (~18 hours)" # Expected
}
\references{
Wickham, H., & Grolemund, G. (n.d.). \emph{R for data science}. Sebastopol, CA:
O'Reilly Media. \url{https://r4ds.had.co.nz}
}
\seealso{
Other utility functions: 
\code{\link{assign_date}()},
\code{\link{cycle_time}()},
\code{\link{pretty_mctq}()},
\code{\link{qplot_walk}()},
\code{\link{random_mctq}()},
\code{\link{raw_data}()},
\code{\link{round_time}()},
\code{\link{shorter_interval}()}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qplot_walk.R
\name{qplot_walk}
\alias{qplot_walk}
\title{Walk through distribution plots}
\usage{
qplot_walk(
  data,
  ...,
  cols = NULL,
  pattern = NULL,
  ignore = "character",
  remove_id = TRUE,
  midday_change = TRUE
)
}
\arguments{
\item{data}{An \code{atomic} or \code{data.frame} object.}

\item{...}{(optional) additional arguments to be passed to
\code{\link[ggplot2:qplot]{ggplot2::qplot()}}.}

\item{cols}{(optional) (only for data frames) a \code{character} object indicating
column names in \code{data} for plotting. If \code{NULL}, \code{qplot_walk()} will use
all columns in \code{data}. This setting only works if \code{pattern = NULL}
(default: \code{NULL}).}

\item{pattern}{(optional) (only for data frames) a string with a regular
expression to select column names in \code{data} for plotting. This setting
only works if \code{cols = NULL} (default: \code{NULL}).}

\item{ignore}{(optional) (only for data frames) a \code{character} object
indicating which object classes the function must ignore. This setting can
be used with \code{cols} and \code{pattern}. Assign \code{NULL} to disable this behavior
(default: \code{"character"}).}

\item{remove_id}{(optional) (only for data frames) a logical value indicating
if the function must ignore column names in \code{data} that match with the
regular expression \code{"^id$|[\\\\._-]id$"} (default: \code{TRUE}).}

\item{midday_change}{(optional) a logical value indicating if the function
must apply a midday change for \code{hms} variables with values greater than
\code{22:00:00} (see Details section to learn more) (default: \code{TRUE}).}
}
\value{
An invisible \code{NULL}. This function don't aim to return values.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

\code{qplot_walk()} helps you to visually assess the distribution of your data.
It uses \code{\link[ggplot2:qplot]{ggplot2::qplot()}} to walk through each selected variable from a
data frame.
}
\details{
\subsection{Requirements}{

This function requires the \code{\link[grDevices:grDevices-package]{grDevices}} and
\code{\link[ggplot2:ggplot2-package]{ggplot2}} package and can only run in interactive
mode. The \code{grDevices} package comes with a standard R installation and is
typically loaded by default. Most people also run R interactively.

If you don't have any of the two packages mentioned above, you can install
them with \code{install.packages("grDevices")} and \code{install.packages("ggplot2")}.
}

\subsection{Plot recover}{

\code{qplot_walk()} erases all plots after it runs. For that reason, the function
first emits a dialog message warning the user of this behavior before it
runs. If you want to recover a single distribution plot, assign the variable
vector to the \code{data} argument.
}

\subsection{Additional arguments to \code{\link[ggplot2:qplot]{ggplot2::qplot()}}}{

\code{qplot_walk()} uses ggplot2 \code{\link[ggplot2:qplot]{ggplot2::qplot()}} to generate plots. If you are
familiar with \code{\link[ggplot2:qplot]{ggplot2::qplot()}}, you can pass additional arguments to the
function using the ellipsis argument (\code{...}).

Note that \code{x}, \code{y}, and \code{data} arguments are reserved for \code{qplot_walk()}.
}

\subsection{\code{Duration}, \code{Period}, and \code{difftime} objects}{

To help with the visualization, \code{qplot_walk()} automatically
converts \code{Duration}, \code{Period}, and \code{difftime} objects to \code{hms}.
}

\subsection{Midday change}{

Time variables with values greater than \code{22:00:00} will automatically be
converted to \code{POSIXct} and be attached to a two-day timeline using the midday
hour as a cutting point, i.e., all values with 12 hours or more will be
placed on day 1, and all the rest will be placed on day 2.

This is made to better represent time vectors that cross the midnight hour.
You can disable this behavior by using \code{midday_change = FALSE}.

Example: Say you have a vector of time values that cross the midnight hour
(e.g., an \code{hms} vector with \code{22:00}, \code{23:00}, \code{00:00}, \code{01:00} values). If
you use \code{midday_change = FALSE}, your data will be represented linearly.\preformatted{00:00 01:00                                22:00 23:00
  |-----|------------------------------------|-----|------->
}

By using \code{midday_change = TRUE} (default), \code{qplot_walk()} will fit your data
to a circular time frame of 24 hours.\preformatted{             day 1                         day 2
                22:00 23:00 00:00 01:00
------------------|-----|-----|-----|---------------------->
}
}

\subsection{\code{id} variables}{

\code{qplot_walk()} will ignore any variable with the follow name pattern
\code{"^id$|[\\\\._-]id$"}, i.e. any variable named \code{id} or that ends with
\code{.id}, \verb{_id}, or \code{-id}.

You can disable this behavior using \code{remove_id = FALSE}.
}
}
\examples{
if (interactive()) {

## Ploting a single column from `data`

qplot_walk(mctq::std_mctq$bt_w)

## Ploting all columns from `data`

qplot_walk(mctq::std_mctq, ignore = NULL, remove_id = FALSE)

## Ploting selected columns from `data`

qplot_walk(mctq::std_mctq, cols = c("bt_w", "msf_sc"))

## Ploting selected columns from `data` using a name pattern

qplot_walk(mctq::std_mctq, pattern = "_w$")

## Examples using other datasets

if (requireNamespace("datasets", quietly = TRUE)) {
    qplot_walk(datasets::iris)
    qplot_walk(datasets::mtcars)
}
}
}
\seealso{
Other utility functions: 
\code{\link{assign_date}()},
\code{\link{cycle_time}()},
\code{\link{pretty_mctq}()},
\code{\link{random_mctq}()},
\code{\link{raw_data}()},
\code{\link{round_time}()},
\code{\link{shorter_interval}()},
\code{\link{sum_time}()}
}
\concept{utility functions}
