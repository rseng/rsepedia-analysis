
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
