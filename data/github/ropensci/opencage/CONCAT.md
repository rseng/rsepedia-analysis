
<!-- README.md is generated from README.Rmd. Please edit that file -->

# opencage

<!-- badges: start -->

[![CRAN
Status](https://www.r-pkg.org/badges/version/opencage)](https://cran.r-project.org/package=opencage)
[![CRAN
Checks](https://cranchecks.info/badges/worst/opencage)](https://cran.r-project.org/web/checks/check_results_opencage.html)
[![CRAN Downloads per
Month](https://cranlogs.r-pkg.org/badges/opencage)](https://cran.r-project.org/package=opencage)
[![R-CMD-check](https://github.com/ropensci/opencage/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/opencage/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/opencage/coverage.svg?branch=main)](https://codecov.io/github/ropensci/opencage?branch=main)
[![rOpenSci
Peer-Review](https://badges.ropensci.org/36_status.svg)](https://github.com/ropensci/software-review/issues/36)
[![License](https://img.shields.io/cran/l/opencage)](https://opensource.org/licenses/gpl-license)

<!-- badges: end -->

Geocode with the [OpenCage](https://opencagedata.com) API, either from
place name to longitude and latitude (forward geocoding) or from
longitude and latitude to the name and address of the location (reverse
geocoding).

## Installation

Install the package with:

``` r
install.packages("opencage")
```

Or install the development version using
[remotes](https://remotes.r-lib.org) with:

``` r
remotes::install_github("ropensci/opencage")
```

## Quickstart

For the best experience, we recommend that you read through the
“Introduction to opencage” vignette (`vignette("opencage")`), but if you
are in a hurry:

1.  Register at
    [opencagedata.com/users/sign\_up](https://opencagedata.com/users/sign_up).
2.  Generate an API key at the [OpenCage
    dashboard](https://opencagedata.com/dashboard#api-keys).
3.  Save your API key as an [environment
    variable](https://rstats.wtf/r-startup.html#renviron) like
    `OPENCAGE_KEY=yourkey` in `.Renviron`. See `help(oc_config)` for
    alternative ways to set your OpenCage API key.

Now you are ready to turn place names into latitude and longitude
coordinates:

``` r
library(opencage)
oc_forward_df(placename = "Sarzeau")
```

    ## # A tibble: 1 x 4
    ##   placename oc_lat oc_lng oc_formatted         
    ##   <chr>      <dbl>  <dbl> <chr>                
    ## 1 Sarzeau     47.5  -2.76 56370 Sarzeau, France

Or turn a set of coordinates into the name and address of the location:

``` r
oc_reverse_df(latitude = 51.5034070, longitude = -0.1275920)
```

    ## # A tibble: 1 x 3
    ##   latitude longitude oc_formatted                                               
    ##      <dbl>     <dbl> <chr>                                                      
    ## 1     51.5    -0.128 Prime Minister’s Office, 10 Downing Street, Westminster, L~

But remember, the vignettes are really great! We have:

-   “Introduction to opencage” `vignette("opencage")`
-   “Customise your query” `vignette("customise_query")`
-   “Output options” `vignette("output_options")`

## About OpenCage

The [OpenCage](https://opencagedata.com) API supports forward and
reverse geocoding. Sources of OpenCage are open geospatial data
including [OpenStreetMap](https://www.openstreetmap.org),
[DataScienceToolkit](https://github.com/petewarden/dstk),
[GeoPlanet](https://en.wikipedia.org/wiki/GeoPlanet), [Natural Earth
Data](https://www.naturalearthdata.com),
[libpostal](https://github.com/openvenues/libpostal),
[GeoNames](https://www.geonames.org), and [Flickr’s
shapefiles](https://code.flickr.net/2009/05/21/flickr-shapefiles-public-dataset-10/)
plus a whole lot more besides. Refer to the current full [list of
credits](https://opencagedata.com/credits).

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/opencage/issues).
-   License: GPL &gt;= 2
-   Get citation information for `opencage` in R doing
    `citation(package = 'opencage')`
-   Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.
# opencage (development version)

* The geocoding functions will not send a query to the API anymore if no API key is present (#133).

# opencage 0.2.2

* Fixed a test that caused an error on CRAN's Solaris (#131).

# opencage 0.2.1

This is a major rewrite of the {opencage} package. `opencage_forward()` and `opencage_reverse()` have been deprecated and are superseded by `oc_forward()` and `oc_reverse()`, respectively. In addition there are two new functions `oc_forward_df()` and `oc_reverse_df()`, which geocode place names or addresses into geographic coordinates (latitude and longitude) or vice versa, and return a data frame. 

The new features include:

* `oc_forward()` and `oc_reverse()` return either lists of data frames, JSON strings, GeoJSON strings, or URLs to be sent to the API (the latter for debugging purposes).
* `oc_forward_df()` and `oc_reverse_df()` take a data frame or vectors as input and return a data frame with the geocoding results, optionally with the source data frame bound to the results data frame. 
* Almost all arguments of the geocoding functions are vectorised (the exceptions being `output`), so it is possible to serially (reverse) geocode lists of locations or coordinates. The geocoding functions show a progress indicator when more than one `placename` or `latitude`/`longitude` pair is provided.
* The forward geocoding functions now support multiple `countrycode`s in accordance with the OpenCage API (#44). The `countrycode`s can now be provided in upper or lower case (#47).
* A helper function `oc_bbox()` now makes it easier to create a list of bounding boxes from numeric vectors, bbox objects or data frames. 
* `oc_forward` and `oc_forward_df` now support [OpenCage's `proximity` parameter](https://blog.opencagedata.com/post/new-optional-parameter-proximity). The results of the geocoding request will be biased towards that location (#60).
* A helper function `oc_points()` now makes it easier to create a list of point coordinates from numeric vectors or data frames to pass to the `proximity` argument for example. 
* All geocoding functions now support [OpenCage's `roadinfo` parameter](https://blog.opencagedata.com/post/new-optional-parameter-roadinfo) (#65). If set to `TRUE`, OpenCage attempts to match the nearest road (rather than an address) and provides additional road and driving information.
* Language tags passed to the `language` argument are not validated anymore, since the language tags used by OpenStreetMap and hence OpenCage do not always conform with the IETF BCP 47 standard (#90). The `languagecodes`, which were stored in {opencage} as external data, have therefore been omitted from the package. In addition, it is now possible to specify `language = "native"`, so OpenCage will attempt to return the [results in the "official" language](https://blog.opencagedata.com/post/support-for-local-language) of the country. 
* http requests are now handled by {[crul](https://docs.ropensci.org/crul/)}, not {[httr](https://httr.r-lib.org)} (#37).
* API calls are now rate limited (#32). The default limit is set to 1 call per second as per the API limit of the [Free Trial plan](https://opencagedata.com/pricing).
* {opencage} settings like the OpenCage API key or the API rate limit can be configured with `oc_config()`. If you want OpenCage to have no record of the contents of your queries, you can also set the `no_record` parameter for the active R session with `oc_config()` (as opposed to providing the parameter with each function call). All `oc_config()` settings can be set more permanently via `options()` or environment variables, see `help(oc_config)`.

## Breaking changes

* `opencage_forward()`, `opencage_reverse()`, and `opencage_key()` are [soft-deprecated](https://lifecycle.r-lib.org/reference/deprecate_soft.html). 
* `opencage_forward()` and `opencage_reverse()` will always output strings as characters, i.e. they won't coerce to factor depending on the `stringsAsFactor` option.
* `opencage_key()` returns the OpenCage API key invisibly.
* `NA` values are not allowed anymore for the `placename` or `latitude`/`longitude` arguments, because OpenCage throws a HTTP 400 ‘bad query’ error when the query is empty (#98). 

## Minor changes

* The column name for `countrycodes` is now `code`, not `Code`. 
* HTTP error messages are now returned directly from the API and are therefore always up-to-date. The previously used responses in `code_message`, which were stored in {opencage} as external data, have been deleted. For more information on OpenCage's HTTP status codes see https://opencagedata.com/api#codes.
* Fixed two URLs, one of which was rejected on the v0.2.0 submission.

# opencage 0.1.4

* Bug fix: now the `countrycode` argument can be used for Namibia (#24, #25).

# opencage 0.1.3

* Added a `add_request` parameter (for appending original query to results).

# opencage 0.1.2

* Added a `abbrv` parameter, see https://blog.opencagedata.com/post/160294347883/shrtr-pls.

# opencage 0.1.1

* Added a `no_record` parameter, see https://blog.opencagedata.com/post/145602604628/more-privacy-with-norecord-parameter

# opencage 0.1.0

* Added a `NEWS.md` file to track changes to the package.
## Release summary

### v0.2.2
This patch version fixes a test that caused an error on CRAN's Solaris (<https://github.com/ropensci/opencage/pull/131>, <https://www.r-project.org/nosvn/R.check/r-patched-solaris-x86/opencage-00check.html>).

## Test environments
* local x86_64-w64-mingw32/x64, R 4.0.4
* GitHub Actions <https://github.com/ropensci/opencage/actions?query=workflow%3AR-CMD-check>:
  * Ubuntu 20.04, R devel, release and oldrel
  * windows-latest, R release
  * macOS-latest, R release
* R-hub:
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* win-builder (devel)

## R CMD check results (local)

Duration: 2m 10.8s

> checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Daniel Possenriede <possenriede+r@gmail.com>'
  
  Days since last update: 4

0 errors √ | 0 warnings √ | 1 note x

## revdepcheck results

We checked 1 reverse dependency, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
# Contributing to opencage

First of all, thanks for considering contributing to {opencage}! 👍
We welcome bug reports and pull requests that expand and improve the functionality of {opencage} from all contributors.
This document outlines how to propose a change to {opencage}. 

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. 

### Share the love ❤️

Think {opencage} is useful? 
Let others discover it, by telling them in person, via your preferred social medium, or a blog post.
Please also share your use case in our discussion forum at [discuss.ropensci.org](https://discuss.ropensci.org). 

Using {opencage} for a paper you are writing? 
Please consider [citing it](https://docs.ropensci.org/opencage/authors.html).
Get citation information for {opencage} in R with `citation(package = 'opencage')`.

### Ask a question ❓

Using {opencage} and got stuck? 
Browse the [documentation](https://docs.ropensci.org/opencage/) to see if you can find a solution. 
Still stuck? 
Post your question on our [discussion forum](https://discuss.ropensci.org) and tag it with the package name.
While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? 
Email the person listed as maintainer in the `DESCRIPTION` file of this repo.
Keep in mind that private discussions over email don't help others - but of course email is totally warranted if it's a sensitive problem of any kind.

### Improve the documentation ✍

Noticed a typo on the website? 
Think a function could use a better example? 
Good documentation makes all the difference, so your help to improve it is very welcome!

Small typos or grammatical errors in documentation can be edited directly using the GitHub web interface, as long as the changes are made in the _source_ file.

This means you should

* edit a roxygen comment in a `.R` file below `R/`, not the `.Rd` files below `man/`.
* edit the `README.Rmd` file, not the `README.md` file in the package root directory.

Since we use a non-standard workflow to render the vignettes in this package, you should 

* edit the `*.Rmd.src` files in the `vignettes/` directory, not the `*.Rmd` files there.

### Reporting an issue 🐛

Using our_package and discovered a bug? 
That's annoying! 
Don't let others have the same experience and open an [issue report on GitHub](https://github.com/ropensci/opencage/issues/new) so we can fix it. 
Please illustrate the bug with a minimal working example, also known as a [reprex](https://www.tidyverse.org/help/#reprex), i.e. please provide detailed steps to reproduce the bug and any information that might be helpful in troubleshooting. The {[reprex](https://reprex.tidyverse.org/)} 📦 can help you with this. 

### Contribute code  🛠

Care to fix bugs or implement new functionality for {opencage}? 
Awesome! 👏
Before you make a substantial change to the package, it is often preferable to first discuss need and scope for the change with the author(s) of the package in an issue report. 

You should then follow the following process:

* Fork the package and clone onto your computer. 
If you haven't done this before, we recommend using `usethis::create_from_github("ropensci/opencage")`.
See the [Pull Request Helper](https://usethis.r-lib.org/articles/articles/pr-functions.html) vignette for more details on how {[usethis](https://usethis.r-lib.org/)} can assist you with contributing code via pull requests (PR), .
* Install all development dependencies with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
* Create a Git branch for each issue you want to address. 
We recommend using `usethis::pr_init("brief-description-of-change")`.
* Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
The title of your PR should briefly describe the change; the body of your PR should contain "Fixes [#issue-number]".
* Add a bullet point to the top of `NEWS.md` describing the changes made followed by your GitHub username, and links to relevant issue(s)/PR(s).

You should also consider the following:

* Keep the changes in your PR as small and succinct as possible. 
Most importantly only address one issue per PR. 
This makes it easier for us to review and merge your PR. 
* We mostly follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the {[styler](https://styler.r-lib.org/)} package to apply these styles, but please do not restyle code that has nothing to do with your PR. 
* We use {[roxygen2](https://roxygen2.r-lib.org/)}, with [Markdown syntax](https://roxygen2.r-lib.org/articles/rd-formatting.html), for documentation.
* We would prefer it if your PR also included unit tests. 
Contributions with test cases included are easier to accept and unit tests ensure that the functionality you just added will not break in the future.
We use {[testthat](https://testthat.r-lib.org/)} for unit tests and we track test coverage with [covr](https://covr.r-lib.org/) and [Codecov](https://codecov.io/).
* We use {[lintr](https://github.com/jimhester/lintr)} for [static code analysis](https://github.com/jimhester/lintr).
* We use [GitHub Actions](https://docs.github.com/en/actions) for continuous integration. 
Workflows are adapted from [r-lib/actions](https://github.com/r-lib/actions). 
Unfortunately tests requiring an API key will not run on a PR, because neither our nor your API key is available there to prevent it from leaking. 

## rOpenSci discussion forum 👄

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, a use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer discussion.

## License 📜

{opencage} is licensed under the [GPL-2 or later](https://opensource.org/licenses/gpl-license).

## Thanks for contributing! 🙏

For more detailed info about contributing to rOpenSci, please see the [rOpenSci Community Contributing Guide](https://contributing.ropensci.org/). 
<!-- 

Thank you for opening a pull request! 

PLEASE DO NOT SHARE YOUR API KEYS HERE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY! 

Before continuing, please make sure that you 

- have read and followed the guidelines in our [Contributing Guide](https://github.com/ropensci/opencage/blob/master/.github/CONTRIBUTING.md).

Provide a general summary of your changes in the Title above.

-->

## Description

<!--- 
Please describe your changes in detail. 
-->

- [ ] example included
<!--- 
If you are introducing a new feature or changing behaviour of existing
methods/functions, please include a brief example if possible! 
-->
- [ ] related issue mentioned
<!--- 
If this closes an issue make sure include e.g., "fixes #4" or similar - 
or if it just relates to an issue make sure to mention it like "#4" 
-->
- [ ] documentation included
<!--- 
If you are introducing a new feature or changing behaviour of existing
methods/functions, please document the new behaviour and features!
This ideally includes the README and/or vignettes, if appropriate.
-->
- [ ] tests included
<!--- 
Please include tests if you made changes to the package code! 
-->
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''
---

<!-- 

Thank you for writing a bug report! This will help us to improve the package!

REMEMBER: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY 

Before continuing, please make sure that you:

- have checked that there does not exist an issue report for your problem already.
- you are using the latest version of {opencage} and all its dependencies.
- the problem is not already fixed in the development version (install it with `remotes::install_github("ropensci/opencage")`).
- have read the guidelines in our [Contributing Guide](https://github.com/ropensci/opencage/blob/master/.github/CONTRIBUTING.md).

-->

## Description & steps to reproduce
<!--
Please provide a clear and concise description of what the bug is.
Please illustrate the bug with a minimal working example or [reprex](https://www.tidyverse.org/help/#reprex). 
Please provide detailed steps to reproduce the bug and any information that might be helpful in troubleshooting.
Please consider using the {reprex} package for this: https://reprex.tidyverse.org/.
What did you expect to happen instead?
-->

## Session information
<!--
Please include your session information by including the output from `devtools::session_info()` or `sessionInfo()`.
-->
---
name: Feature request
about: Suggest an idea for {opencage}
title: ''
labels: ''
assignees: ''

---

<!-- 

Thank you for writing a feature request! This will help us to improve the package!

REMEMBER: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY 

Before continuing, please make sure that you:

- have checked that there does not exist a similar request already.
- have read the guidelines in our [Contributing Guide](https://github.com/ropensci/opencage/blob/master/.github/CONTRIBUTING.md).

-->

## Description

<!-- 

Please describe what feature you would like to have implemented in {opencage}.

- Describe the problem you are trying to solve and why you cannot solve it yet.
- Describe the solution you'd like to see. 
Try to think about what aspects are important to you, rather than implementation details. 
- Describe alternative solutions or features you have considered.
- What else should we know to understand your problem? 

-->
## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.4 (2021-02-15) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |en                           |
|collate  |German_Germany.1252          |
|ctype    |German_Germany.1252          |
|tz       |Europe/Berlin                |
|date     |2021-02-18                   |

# Dependencies

|package     |old    |new    |<U+0394>  |
|:-----------|:------|:------|:--|
|opencage    |0.2.1  |0.2.2  |*  |
|assertthat  |0.2.1  |0.2.1  |   |
|cachem      |1.0.4  |1.0.4  |   |
|cli         |2.3.0  |2.3.0  |   |
|cpp11       |0.2.6  |0.2.6  |   |
|crayon      |1.4.1  |1.4.1  |   |
|crul        |1.1.0  |1.1.0  |   |
|curl        |4.3    |4.3    |   |
|digest      |0.6.27 |0.6.27 |   |
|dplyr       |1.0.4  |1.0.4  |   |
|ellipsis    |0.3.1  |0.3.1  |   |
|fansi       |0.4.2  |0.4.2  |   |
|fastmap     |1.1.0  |1.1.0  |   |
|generics    |0.1.0  |0.1.0  |   |
|glue        |1.4.2  |1.4.2  |   |
|hms         |1.0.0  |1.0.0  |   |
|httpcode    |0.3.0  |0.3.0  |   |
|jsonlite    |1.7.2  |1.7.2  |   |
|lifecycle   |1.0.0  |1.0.0  |   |
|magrittr    |2.0.1  |2.0.1  |   |
|memoise     |2.0.0  |2.0.0  |   |
|mime        |0.10   |0.10   |   |
|pillar      |1.4.7  |1.4.7  |   |
|pkgconfig   |2.0.3  |2.0.3  |   |
|prettyunits |1.1.1  |1.1.1  |   |
|progress    |1.2.2  |1.2.2  |   |
|purrr       |0.3.4  |0.3.4  |   |
|R6          |2.5.0  |2.5.0  |   |
|ratelimitr  |0.4.1  |0.4.1  |   |
|Rcpp        |1.0.6  |1.0.6  |   |
|rlang       |0.4.10 |0.4.10 |   |
|tibble      |3.0.6  |3.0.6  |   |
|tidyr       |1.1.2  |1.1.2  |   |
|tidyselect  |1.1.0  |1.1.0  |   |
|triebeard   |0.3.0  |0.3.0  |   |
|urltools    |1.7.3  |1.7.3  |   |
|utf8        |1.1.4  |1.1.4  |   |
|vctrs       |0.3.6  |0.3.6  |   |
|withr       |2.4.1  |2.4.1  |   |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output:
  github_document
---
<!-- README.md is generated from README.Rmd. Please edit that file -->
# opencage

<!-- badges: start -->

[![CRAN Status](https://www.r-pkg.org/badges/version/opencage)](https://cran.r-project.org/package=opencage)
[![CRAN Checks](https://cranchecks.info/badges/worst/opencage)](https://cran.r-project.org/web/checks/check_results_opencage.html)
[![CRAN Downloads per Month](https://cranlogs.r-pkg.org/badges/opencage)](https://cran.r-project.org/package=opencage)
[![R-CMD-check](https://github.com/ropensci/opencage/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/opencage/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/opencage/coverage.svg?branch=main)](https://codecov.io/github/ropensci/opencage?branch=main)
[![rOpenSci Peer-Review](https://badges.ropensci.org/36_status.svg)](https://github.com/ropensci/software-review/issues/36)
[![License](https://img.shields.io/cran/l/opencage)](https://opensource.org/licenses/gpl-license)

<!-- badges: end -->

Geocode with the [OpenCage](https://opencagedata.com) API, either from place name to longitude and latitude (forward geocoding) or from longitude and latitude to the name and address of the location (reverse geocoding). 

## Installation

Install the package with:

```{r eval = FALSE}
install.packages("opencage")
```

Or install the development version using [remotes](https://remotes.r-lib.org) with:

```{r, eval = FALSE}
remotes::install_github("ropensci/opencage")

```

## Quickstart

For the best experience, we recommend that you read through the "Introduction to opencage" vignette (`vignette("opencage")`), but if you are in a hurry:

1. Register at [opencagedata.com/users/sign_up](https://opencagedata.com/users/sign_up).
2. Generate an API key at the [OpenCage dashboard](https://opencagedata.com/dashboard#api-keys).
3. Save your API key as an [environment variable](https://rstats.wtf/r-startup.html#renviron) like `OPENCAGE_KEY=yourkey` in `.Renviron`. 
See `help(oc_config)` for alternative ways to set your OpenCage API key.

Now you are ready to turn place names into latitude and longitude coordinates:

```{r forward}
library(opencage)
oc_forward_df(placename = "Sarzeau")
```

Or turn a set of coordinates into the name and address of the location:

```{r reverse}
oc_reverse_df(latitude = 51.5034070, longitude = -0.1275920)
```

But remember, the vignettes are really great! We have:

- "Introduction to opencage" `vignette("opencage")`
- "Customise your query" `vignette("customise_query")`
- "Output options" `vignette("output_options")`

## About OpenCage

The [OpenCage](https://opencagedata.com) API supports forward and reverse geocoding. 
Sources of OpenCage are open geospatial data including 
[OpenStreetMap](https://www.openstreetmap.org), 
[DataScienceToolkit](https://github.com/petewarden/dstk), 
[GeoPlanet](https://en.wikipedia.org/wiki/GeoPlanet), 
[Natural Earth Data](https://www.naturalearthdata.com), 
[libpostal](https://github.com/openvenues/libpostal), 
[GeoNames](https://www.geonames.org), and 
[Flickr's shapefiles](https://code.flickr.net/2009/05/21/flickr-shapefiles-public-dataset-10/) plus a whole lot more besides. 
Refer to the current full [list of credits](https://opencagedata.com/credits).

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/opencage/issues).
* License: GPL >= 2
* Get citation information for `opencage` in R doing `citation(package = 'opencage')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Introduction to opencage"
subtitle: "Forward and Reverse Geocoding"
author: "Daniel Possenriede, Jesse Sadler, Maëlle Salmon"
date: "2021-02-27"
description: >
  "Get started with the opencage R package to geocode with the OpenCage API, either from place name to longitude and latitude (forward geocoding) or from longitude and latitude to the name and address of a location (reverse geocoding)."
output:
  rmarkdown::html_vignette:
    df_print: kable
vignette: >
  %\VignetteIndexEntry{Introduction to opencage}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Geocoding is the process of converting place names or addresses into geographic coordinates – like latitude and longitude – or vice versa.
With {opencage} you can geocode using the OpenCage API, either from place name to longitude and latitude (forward geocoding) or from longitude and latitude to the name and address of a location (reverse geocoding).

This vignette covers the setup process for working with {opencage} and basic workflows for both forward and reverse geocoding. Make sure to also check out "[Customise your query](customise_query.html)" if you want a deeper dive into customizing the various parameters available through the OpenCage API. The "[Output options](output_options.html)" vignette shows additional workflows by modifying the form in which the geocoding results are returned.

## Setup

Before you can use the {opencage} package and query the OpenCage API you need to first register with OpenCage.
Additionally, you may want to set a rate limit (if you have a paid [OpenCage plan](https://opencagedata.com/pricing)), and you might want to prevent OpenCage from storing the content of your queries.
In other words, you need to setup {opencage}, so let's go through the process.

### Authentication

To use the package and authenticate yourself with the OpenCage API, you will need to register at [opencagedata.com/users/sign_up](https://opencagedata.com/users/sign_up) to get an API key.
The "Free Trial" plan provides up to 2,500 API requests a day.
There are paid plans available, if you need to run more API requests.
After you have registered, you can generate an API key with the [OpenCage dashboard](https://opencagedata.com/dashboard#api-keys).

Now we need to ensure that the functions in {opencage} can access your API key.
{opencage} will conveniently retrieve your API key if it is saved in the environment variable `"OPENCAGE_KEY"`.
If it is not, `oc_config()` will help to set that environment variable.

Do not pass the key directly as a parameter to the function.
Doing so risks exposing your API key via your script or your history.
There are three safer ways to set your API key instead:

1. Save your API key as an environment variable in `.Renviron` as described in [What They Forgot to Teach You About R](https://rstats.wtf/r-startup.html#renviron) or [Efficient R Programming](https://csgillespie.github.io/efficientR/set-up.html#renviron).
From there it will be fetched by all functions that call the OpenCage API.
You do not even have to call `oc_config()` to set your key;
you can start geocoding right away.
If you have the {[usethis](https://usethis.r-lib.org)} package installed, you can edit your `.Renviron` with `usethis::edit_r_environ()`.
We strongly recommend storing your API key in the user-level `.Renviron`, as opposed to the project-level.
This makes it less likely you will share sensitive information by mistake.

2. If you use a package like {[keyring](https://github.com/r-lib/keyring)} to store your credentials, you can safely pass your key in a script with a function call like `oc_config(key = keyring::key_get("opencage"))`.

3. If you call `oc_config()` in an `interactive()` session and the `OPENCAGE_KEY` environment variable is not set, it will prompt you to enter the key in the console.

Whatever method you choose, keep your API key secret. OpenCage also features [best practices for keeping your API key safe](https://opencagedata.com/guides/how-to-protect-your-api-key).

### Rate limit

A rate limit is used to control the rate of requests sent, so legitimate requests do not lead to an unintended Denial of Service attack.
The rate limit allowed by the API depends on the OpenCage plan you have and ranges from 1 request/sec for the "Free Trial" plan to 15 requests/sec for the "Medium" or "Large" plans.
See [opencagedata.com/pricing](https://opencagedata.com/pricing) for details and up-to-date information.

If you have a "Free Trial" account with OpenCage, you can skip to the next section, because the rate limit is already set correctly for you at 1 request/sec.

If you have a paid account, you can set the rate limit for the active R session with `oc_config(rate_sec = n)` where `n` is the appropriate rate limit.
You can set the rate limit persistently across sessions by setting an `oc_rate_sec` option in your `.Rprofile`.
If you have the {[usethis](https://usethis.r-lib.org)} package installed, you can edit your `.Rprofile` with `usethis::edit_r_profile()`.

### Privacy

By default, OpenCage will store your queries on its server logs and will cache the forward geocoding requests on their side.
They do this to speed up response times and to be able to debug errors and improve their service.
Logs are automatically deleted after six months according to OpenCage's [page on data protection and GDPR](https://opencagedata.com/gdpr).

If you have concerns about privacy and want OpenCage to have no record of your query, i.e. the place name or latitude and longitude coordinates you want to geocode, you can set a `no_record` parameter to `TRUE`, which tells the API to not log    nor cache the queries.
OpenCage still records that you made a request, but not the specific queries you made.

`oc_config(no_record = TRUE)` sets an `oc_no_record` option for the active R session, so it will be used for all subsequent OpenCage queries.
You can set the `oc_no_record` option persistently across sessions in your `.Rprofile`.

For more information on OpenCage's policies on privacy and data protection see the Legal section in [their FAQs](https://opencagedata.com/faq#legal), their [GDPR page](https://opencagedata.com/gdpr), and, for the `no_record` parameter specifically, see the relevant [blog post](https://blog.opencagedata.com/post/145602604628/more-privacy-with-norecord-parameter).

For increased privacy, {opencage} sets `no_record` to `TRUE`, by default.
Please note, however, that {opencage} always caches the data it receives from the OpenCage API locally, but only for as long as your R session is alive (see below).

### (Don't) show API key

`oc_config()` has another argument, `show_key`. This is only used for debugging and we will explain it in more detail in `vignette("output_options")`.
For now suffice it to say that your OpenCage API key will not be shown in any {opencage} output, unless you change this setting.

### Altogether now

In sum, if you want to set your API key with {keyring}, set the rate limit to 10 (only do this if you have a paid account, please!), and do not want OpenCage to have records of your queries, you would configure {opencage} for the active session like this:


```r
library("opencage")
oc_config(
  key = keyring::key_get("opencage"),
  rate_sec = 10,
  no_record = TRUE
)
```

## Forward geocoding

Now you can start to geocode. Forward geocoding is from location name(s) to latitude and longitude tuple(s).


```r
oc_forward_df(placename = "Sarzeau")
#> # A tibble: 1 x 4
#>   placename oc_lat oc_lng oc_formatted         
#>   <chr>      <dbl>  <dbl> <chr>                
#> 1 Sarzeau     47.5  -2.76 56370 Sarzeau, France
```

All geocoding functions are vectorised, i.e. you can geocode multiple locations with one function call.
Note that behind the scenes the requests are still sent to the API one-by-one.


```r
opera <- c("Palacio de Bellas Artes", "Scala", "Sydney Opera House")
oc_forward_df(placename = opera)
#> # A tibble: 3 x 4
#>   placename               oc_lat oc_lng oc_formatted                                                                     
#>   <chr>                    <dbl>  <dbl> <chr>                                                                            
#> 1 Palacio de Bellas Artes   19.4 -99.1  Palacio de Bellas Artes, Avenida Juárez, Centro Urbano, 06050 Mexico City, Mexico
#> 2 Scala                     45.5   9.19 Teatro alla Scala, Largo Antonio Ghiringhelli, Milan Milan, Italy                
#> 3 Sydney Opera House       -33.9 151.   Sydney Opera House, 2 Macquarie Street, Sydney NSW 2000, Australia
```

By default, `oc_forward_df()` only returns three results columns: `oc_lat` (for latitude), `oc_lon` (for longitude), and  `oc_formatted` (the formatted address).
As you can see, the results columns are all prefixed with `oc_`.
If you specify `oc_forward_df(output = all)`, you will receive all result columns, which are often quite extensive.
Which columns you receive exactly depends on the information OpenCage returns for each specific request.


```r
oc_forward_df(placename = opera, output = "all")
#> # A tibble: 3 x 30
#>   placename  oc_lat oc_lng oc_confidence oc_formatted  oc_northeast_lat oc_northeast_lng oc_southwest_lat oc_southwest_lng oc_iso_3166_1_a~
#>   <chr>       <dbl>  <dbl>         <int> <chr>                    <dbl>            <dbl>            <dbl>            <dbl> <chr>           
#> 1 Palacio d~   19.4 -99.1              9 Palacio de B~             19.4           -99.1              19.4           -99.1  MX              
#> 2 Scala        45.5   9.19             9 Teatro alla ~             45.5             9.19             45.5             9.19 IT              
#> 3 Sydney Op~  -33.9 151.               9 Sydney Opera~            -33.9           151.              -33.9           151.   AU              
#> # ... with 20 more variables: oc_iso_3166_1_alpha_3 <chr>, oc_category <chr>, oc_type <chr>, oc_city <chr>, oc_continent <chr>,
#> #   oc_country <chr>, oc_country_code <chr>, oc_museum <chr>, oc_neighbourhood <chr>, oc_postcode <chr>, oc_road <chr>,
#> #   oc_attraction <chr>, oc_county <chr>, oc_political_union <chr>, oc_state <chr>, oc_suburb <chr>, oc_arts_centre <chr>,
#> #   oc_house_number <chr>, oc_municipality <chr>, oc_state_code <chr>
```

You can also pass a data frame to `oc_forward_df()`.
By default the results columns are added to the input data frame, which is useful for keeping information associated with the place names that are in separate columns.
If you want a data frame with only the geocoding results, set `bind_cols = FALSE`.


```r
concert_df <-
  data.frame(location = c("Elbphilharmonie", "Concertgebouw", "Suntory Hall"))
oc_forward_df(data = concert_df, placename = location)
#> # A tibble: 3 x 4
#>   location        oc_lat oc_lng oc_formatted                                                                 
#>   <chr>            <dbl>  <dbl> <chr>                                                                        
#> 1 Elbphilharmonie   53.5   9.98 Elbe Philharmonic Hall, Platz der Deutschen Einheit 1, 20457 Hamburg, Germany
#> 2 Concertgebouw     52.4   4.88 Concertgebouw, Concertgebouwplein 2, 1071 LN Amsterdam, Netherlands          
#> 3 Suntory Hall      35.7 140.   Suntory Hall, Karayan Plaza, Azabu, Minato, Japan
```

You can use it in a piped workflow as well.


```r
library(dplyr, warn.conflicts = FALSE)
concert_df %>% oc_forward_df(location)
#> # A tibble: 3 x 4
#>   location        oc_lat oc_lng oc_formatted                                                                 
#>   <chr>            <dbl>  <dbl> <chr>                                                                        
#> 1 Elbphilharmonie   53.5   9.98 Elbe Philharmonic Hall, Platz der Deutschen Einheit 1, 20457 Hamburg, Germany
#> 2 Concertgebouw     52.4   4.88 Concertgebouw, Concertgebouwplein 2, 1071 LN Amsterdam, Netherlands          
#> 3 Suntory Hall      35.7 140.   Suntory Hall, Karayan Plaza, Azabu, Minato, Japan
```

## Reverse geocoding

Reverse geocoding works in the opposite direction of forward geocoding: from a pair of coordinates to the name and address most appropriate for the coordinates.


```r
oc_reverse_df(latitude = 51.5034070, longitude = -0.1275920)
#> # A tibble: 1 x 3
#>   latitude longitude oc_formatted                                                                            
#>      <dbl>     <dbl> <chr>                                                                                   
#> 1     51.5    -0.128 Prime Minister’s Office, Westminster, 10 Downing Street, London SW1A 2AA, United Kingdom
```

Note that all coordinates sent to the OpenCage API must adhere to the [WGS 84](https://en.wikipedia.org/wiki/World_Geodetic_System) (also known as [EPSG:4326](https://epsg.io/4326)) [coordinate reference system](https://en.wikipedia.org/wiki/Spatial_reference_system) in decimal format.
This is the coordinate reference system used by the [Global Positioning System](https://en.wikipedia.org/wiki/Global_Positioning_System).
There is usually no reason to send more than six or seven digits past the decimal.
Any further precision gets to the [level of a centimeter](https://en.wikipedia.org/wiki/Decimal_degrees).

Like `oc_forward_df()`, `oc_reverse_df()` is vectorised, can work with numeric vectors and data frames, supports the `output = "all"` argument and can be used with the {magrittr} pipe.

OpenCage only returns at most [one result](https://opencagedata.com/api#ranking) per reverse geocoding request.

## Caching

OpenCage [allows and supports caching](https://opencagedata.com/api#caching).
To minimize the number of requests sent to the API {opencage} uses {[memoise](https://github.com/r-lib/memoise)} to cache results inside the active R session.




```r
system.time(oc_reverse(latitude = 10, longitude = 10))
#>    user  system elapsed 
#>    0.00    0.00    0.86

system.time(oc_reverse(latitude = 10, longitude = 10))
#>    user  system elapsed 
#>    0.01    0.00    0.01
```

To clear the cache of all results either start a new R session or call `oc_clear_cache()`.


```r
oc_clear_cache()
#> [1] TRUE

system.time(oc_reverse(latitude = 10, longitude = 10))
#>    user  system elapsed 
#>    0.03    0.00    0.63
```

As you probably know, cache invalidation is one of the harder things to do in computer science.
Therefore {opencage} only supports invalidating the whole cache and not individual records at the moment.

The underlying data at OpenCage is [updated daily](https://opencagedata.com/faq#general).

## Further information

OpenCage supports a lot of parameters to either target your search area more specifically or to specify what additional information you need.
See the ["Customise your query"](customise_query.html) vignette for details.

Besides `oc_forward_df()` and `oc_reverse_df()`, which always return a single tibble, {opencage} has two sibling functions — `oc_forward()` and `oc_reverse()` — which can be used to return types of output.
Depending on what you specify as the `return` parameter, `oc_forward()` and `oc_reverse()` will return either a list of tibbles (`df_list`, the default), JSON lists (`json_list`), GeoJSON lists (`geojson_list`), or the URL with which the API would be called (`url_only`).
Learn more in the ["Output options"](output_options.html) vignette.

Please report any issues or bugs on [our GitHub repository](https://github.com/ropensci/opencage/issues) and post questions on [discuss.ropensci.org](https://discuss.ropensci.org).
---
title: "Output options"
subtitle: "Get different kinds of output from OpenCage"
author: "Daniel Possenriede, Jesse Sadler, Maëlle Salmon"
date: "2021-02-18"
description: >
  "`oc_forward()`/`oc_reverse()` return lists of various type, namely data frames, JSON, GeoJSON or URLs, depending on the `return` value you specify. The possible `return` values are `df_list`, `json_list`, `geojson_list` and `url_only`."
output:
  rmarkdown::html_vignette:
    df_print: kable
vignette: >
  %\VignetteIndexEntry{Output options}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



{opencage} contains two main types of main functions distinguished by the manner in which the results are returned:

1. The `oc_forward_df()`/`oc_reverse_df()` functions always return a single tibble.
2. The `oc_forward()`/`oc_reverse()` functions return a list of various types, depending on the `return` value you specify.
The possible `return` values are `df_list`, `json_list`, `geojson_list` and `url_only`.

Use of the `oc_forward_df()`/`oc_reverse_df()` functions is demonstrated in the "Introduction to opencage" and the "Customise your query" vignettes (`vignette("opencage")` and `vignette("customise_query")`, respectively) .
The function arguments mentioned in these other vignettes are also generally available with `oc_forward()`/`oc_reverse()`.
Here we will show the different return values available with `oc_forward()`/`oc_reverse()`.

## `df_list`

The default return value is `df_list`.
It returns a list of tibbles.


```r
stations <- c("Casey Station", "McMurdo Station")
oc_forward(stations, return = "df_list")
#> [[1]]
#> # A tibble: 2 x 18
#>   oc_confidence oc_formatted oc_northeast_lat oc_northeast_lng oc_southwest_lat oc_southwest_lng oc_category oc_type oc_continent oc_hamlet
#>           <int> <chr>                   <dbl>            <dbl>            <dbl>            <dbl> <chr>       <chr>   <chr>        <chr>    
#> 1             9 Casey Stati~            -66.3             111.            -66.3             111. commerce    office  Antarctica   Casey St~
#> 2             9 Casey Stati~             NA                NA              NA                NA  travel/tou~ point_~ Antarctica   <NA>     
#> # ... with 8 more variables: oc_office <chr>, oc_iso_3166_1_alpha_2 <chr>, oc_iso_3166_1_alpha_3 <chr>, oc_country <chr>,
#> #   oc_country_code <chr>, oc_point_of_interest <chr>, oc_lat <dbl>, oc_lng <dbl>
#> 
#> [[2]]
#> # A tibble: 1 x 12
#>   oc_confidence oc_formatted oc_northeast_lat oc_northeast_lng oc_southwest_lat oc_southwest_lng oc_category oc_type oc_continent oc_town
#>           <int> <chr>                   <dbl>            <dbl>            <dbl>            <dbl> <chr>       <chr>   <chr>        <chr>  
#> 1             7 McMurdo Sta~            -77.8             167.            -77.9             167. place       city    Antarctica   McMurd~
#> # ... with 2 more variables: oc_lat <dbl>, oc_lng <dbl>
```

The `df_list` type drives the `oc_forward_df()`/`oc_reverse_df()` functions.
You can use the `df_list` output in a `dplyr::mutate()` chain to replicate the functionality of `oc_forward_df()`:


```r
library(dplyr, warn.conflicts = FALSE)

oc_data <-
  tibble(place = stations) %>%
  mutate(oc_result = oc_forward(place))

oc_data
#> # A tibble: 2 x 2
#>   place           oc_result        
#>   <chr>           <list>           
#> 1 Casey Station   <tibble [2 x 18]>
#> 2 McMurdo Station <tibble [1 x 12]>
```

This creates a list column `oc_result`, which can be easily unnested with `tidyr::unnest()`:


```r
library(tidyr, warn.conflicts = FALSE)

oc_data %>% unnest(oc_result)
#> # A tibble: 3 x 20
#>   place oc_confidence oc_formatted oc_northeast_lat oc_northeast_lng oc_southwest_lat oc_southwest_lng oc_category oc_type oc_continent
#>   <chr>         <int> <chr>                   <dbl>            <dbl>            <dbl>            <dbl> <chr>       <chr>   <chr>       
#> 1 Case~             9 Casey Stati~            -66.3             111.            -66.3             111. commerce    office  Antarctica  
#> 2 Case~             9 Casey Stati~             NA                NA              NA                NA  travel/tou~ point_~ Antarctica  
#> 3 McMu~             7 McMurdo Sta~            -77.8             167.            -77.9             167. place       city    Antarctica  
#> # ... with 10 more variables: oc_hamlet <chr>, oc_office <chr>, oc_iso_3166_1_alpha_2 <chr>, oc_iso_3166_1_alpha_3 <chr>,
#> #   oc_country <chr>, oc_country_code <chr>, oc_point_of_interest <chr>, oc_lat <dbl>, oc_lng <dbl>, oc_town <chr>
```

## `json_list`

OpenCage's main output format is JSON. When you specify `json_list` as the return type, you get the JSON as an R `list()`.


```r
oc_forward("Casey Station", return = "json_list")
#> [[1]]
#> [[1]]$documentation
#> [1] "https://opencagedata.com/api"
#> 
#> [[1]]$licenses
#> [[1]]$licenses[[1]]
#> [[1]]$licenses[[1]]$name
#> [1] "see attribution guide"
#> 
#> [[1]]$licenses[[1]]$url
#> [1] "https://opencagedata.com/credits"
#> 
#> 
#> 
#> [[1]]$results
#> [[1]]$results[[1]]
#> [[1]]$results[[1]]$bounds
#> [[1]]$results[[1]]$bounds$northeast
#> [[1]]$results[[1]]$bounds$northeast$lat
#> [1] -66.28255
#> 
#> [[1]]$results[[1]]$bounds$northeast$lng
#> [1] 110.5267
#> 
#> 
#> [[1]]$results[[1]]$bounds$southwest
#> [[1]]$results[[1]]$bounds$southwest$lat
#> [1] -66.28265
#> 
#> [[1]]$results[[1]]$bounds$southwest$lng
#> [1] 110.5266
#> 
#> 
#> 
#> [[1]]$results[[1]]$components
#> [[1]]$results[[1]]$components$`_category`
#> [1] "commerce"
#> 
#> [[1]]$results[[1]]$components$`_type`
#> [1] "office"
#> 
#> [[1]]$results[[1]]$components$continent
#> [1] "Antarctica"
#> 
#> [[1]]$results[[1]]$components$hamlet
#> [1] "Casey Station"
#> 
#> [[1]]$results[[1]]$components$office
#> [1] "Casey Station"
#> 
#> 
#> [[1]]$results[[1]]$confidence
#> [1] 9
#> 
#> [[1]]$results[[1]]$formatted
#> [1] "Casey Station"
#> 
#> [[1]]$results[[1]]$geometry
#> [[1]]$results[[1]]$geometry$lat
#> [1] -66.2826
#> 
#> [[1]]$results[[1]]$geometry$lng
#> [1] 110.5266
#> 
#> 
#> 
#> [[1]]$results[[2]]
#> [[1]]$results[[2]]$components
#> [[1]]$results[[2]]$components$`ISO_3166-1_alpha-2`
#> [1] "AQ"
#> 
#> [[1]]$results[[2]]$components$`ISO_3166-1_alpha-3`
#> [1] "ATA"
#> 
#> [[1]]$results[[2]]$components$`_category`
#> [1] "travel/tourism"
#> 
#> [[1]]$results[[2]]$components$`_type`
#> [1] "point_of_interest"
#> 
#> [[1]]$results[[2]]$components$continent
#> [1] "Antarctica"
#> 
#> [[1]]$results[[2]]$components$country
#> [1] "Antarctica"
#> 
#> [[1]]$results[[2]]$components$country_code
#> [1] "aq"
#> 
#> [[1]]$results[[2]]$components$point_of_interest
#> [1] "Casey Station"
#> 
#> 
#> [[1]]$results[[2]]$confidence
#> [1] 9
#> 
#> [[1]]$results[[2]]$formatted
#> [1] "Casey Station, Antarctica"
#> 
#> [[1]]$results[[2]]$geometry
#> [[1]]$results[[2]]$geometry$lat
#> [1] -66.28225
#> 
#> [[1]]$results[[2]]$geometry$lng
#> [1] 110.5278
#> 
#> 
#> 
#> 
#> [[1]]$status
#> [[1]]$status$code
#> [1] 200
#> 
#> [[1]]$status$message
#> [1] "OK"
#> 
#> 
#> [[1]]$stay_informed
#> [[1]]$stay_informed$blog
#> [1] "https://blog.opencagedata.com"
#> 
#> [[1]]$stay_informed$twitter
#> [1] "https://twitter.com/OpenCage"
#> 
#> 
#> [[1]]$thanks
#> [1] "For using an OpenCage API"
#> 
#> [[1]]$timestamp
#> [[1]]$timestamp$created_http
#> [1] "Thu, 18 Feb 2021 19:02:13 GMT"
#> 
#> [[1]]$timestamp$created_unix
#> [1] 1613674933
#> 
#> 
#> [[1]]$total_results
#> [1] 2
```

## `geojson_list`

When you choose `geojson_list` as the return type, the geocoder response will be returned as GeoJSON specified as an R `list()`.


```r
gjsn_lst <- oc_forward("Casey Station", return = "geojson_list")
gjsn_lst
#> [[1]]
#> $documentation
#> [1] "https://opencagedata.com/api"
#> 
#> $features
#> $features[[1]]
#> $features[[1]]$geometry
#> $features[[1]]$geometry$coordinates
#> $features[[1]]$geometry$coordinates[[1]]
#> [1] 110.5266
#> 
#> $features[[1]]$geometry$coordinates[[2]]
#> [1] -66.2826
#> 
#> 
#> $features[[1]]$geometry$type
#> [1] "Point"
#> 
#> 
#> $features[[1]]$properties
#> $features[[1]]$properties$bounds
#> $features[[1]]$properties$bounds$northeast
#> $features[[1]]$properties$bounds$northeast$lat
#> [1] -66.28255
#> 
#> $features[[1]]$properties$bounds$northeast$lng
#> [1] 110.5267
#> 
#> 
#> $features[[1]]$properties$bounds$southwest
#> $features[[1]]$properties$bounds$southwest$lat
#> [1] -66.28265
#> 
#> $features[[1]]$properties$bounds$southwest$lng
#> [1] 110.5266
#> 
#> 
#> 
#> $features[[1]]$properties$components
#> $features[[1]]$properties$components$`_category`
#> [1] "commerce"
#> 
#> $features[[1]]$properties$components$`_type`
#> [1] "office"
#> 
#> $features[[1]]$properties$components$continent
#> [1] "Antarctica"
#> 
#> $features[[1]]$properties$components$hamlet
#> [1] "Casey Station"
#> 
#> $features[[1]]$properties$components$office
#> [1] "Casey Station"
#> 
#> 
#> $features[[1]]$properties$confidence
#> [1] 9
#> 
#> $features[[1]]$properties$formatted
#> [1] "Casey Station"
#> 
#> 
#> $features[[1]]$type
#> [1] "Feature"
#> 
#> 
#> $features[[2]]
#> $features[[2]]$geometry
#> $features[[2]]$geometry$coordinates
#> $features[[2]]$geometry$coordinates[[1]]
#> [1] 110.5278
#> 
#> $features[[2]]$geometry$coordinates[[2]]
#> [1] -66.28225
#> 
#> 
#> $features[[2]]$geometry$type
#> [1] "Point"
#> 
#> 
#> $features[[2]]$properties
#> $features[[2]]$properties$components
#> $features[[2]]$properties$components$`ISO_3166-1_alpha-2`
#> [1] "AQ"
#> 
#> $features[[2]]$properties$components$`ISO_3166-1_alpha-3`
#> [1] "ATA"
#> 
#> $features[[2]]$properties$components$`_category`
#> [1] "travel/tourism"
#> 
#> $features[[2]]$properties$components$`_type`
#> [1] "point_of_interest"
#> 
#> $features[[2]]$properties$components$continent
#> [1] "Antarctica"
#> 
#> $features[[2]]$properties$components$country
#> [1] "Antarctica"
#> 
#> $features[[2]]$properties$components$country_code
#> [1] "aq"
#> 
#> $features[[2]]$properties$components$point_of_interest
#> [1] "Casey Station"
#> 
#> 
#> $features[[2]]$properties$confidence
#> [1] 9
#> 
#> $features[[2]]$properties$formatted
#> [1] "Casey Station, Antarctica"
#> 
#> 
#> $features[[2]]$type
#> [1] "Feature"
#> 
#> 
#> 
#> $licenses
#> $licenses[[1]]
#> $licenses[[1]]$name
#> [1] "see attribution guide"
#> 
#> $licenses[[1]]$url
#> [1] "https://opencagedata.com/credits"
#> 
#> 
#> 
#> $rate
#> named list()
#> 
#> $status
#> $status$code
#> [1] 200
#> 
#> $status$message
#> [1] "OK"
#> 
#> 
#> $stay_informed
#> $stay_informed$blog
#> [1] "https://blog.opencagedata.com"
#> 
#> $stay_informed$twitter
#> [1] "https://twitter.com/OpenCage"
#> 
#> 
#> $thanks
#> [1] "For using an OpenCage API"
#> 
#> $timestamp
#> $timestamp$created_http
#> [1] "Thu, 18 Feb 2021 19:02:15 GMT"
#> 
#> $timestamp$created_unix
#> [1] 1613674935
#> 
#> 
#> $total_results
#> [1] 2
#> 
#> $type
#> [1] "FeatureCollection"
#> 
#> attr(,"class")
#> [1] "geo_list"
```

In fact, {opencage} returns a list of results in `geo_list` format, which should be compatible with the {[geojsonio](https://docs.ropensci.org/geojsonio/)} package.


```r
class(gjsn_lst[[1]])
#> [1] "geo_list"
```

## `url_only`

`url_only` returns the OpenCage URL for debugging purposes.


```r
oc_forward("Casey Station", return = "url_only")
#> [[1]]
#> [1] "https://api.opencagedata.com/geocode/v1/json?q=Casey%20Station&limit=10&no_annotations=1&roadinfo=0&no_dedupe=0&no_record=1&abbrv=0&add_request=0&key=OPENCAGE_KEY"
```

Your OpenCage API key is masked with the `OPENCAGE_KEY` string, by default.
If you really want {opencage} to display your API key with the URL, set the `show_key` argument in `oc_config()` to `TRUE`.


```r
oc_config(show_key = TRUE)
```

Note that the API key will only be returned with the URL in `base::interactive()` mode.

## `xml`

{opencage} does not support the XML response type at the moment.
Please file an [issue](https://github.com/ropensci/opencage/issues) or a [pull-request](https://github.com/ropensci/opencage/pulls) if you have a use-case that requires this.

## Return query text

`oc_forward()` and `oc_reverse()` have an `add_request` argument, indicating whether the request is returned again with the results.
If the `return` value is a `df_list`, the `placename` or `latitude,longitude` is added as a column to the results without a roundtrip to the API.
`json_list` results will contain all request parameters as returned by the API.
This would normally include your OpenCage API key, but {opencage} masks the key and replaces it with the `OPENCAGE_KEY` string in the output.
`add_request` is currently ignored by OpenCage for GeoJSON results.
---
title: "Customise your query"
subtitle: "Get more and better results from OpenCage"
author: "Daniel Possenriede, Jesse Sadler, Maëlle Salmon"
date: "2021-02-18"
description: >
  "The OpenCage API supports about a dozen parameters to customise a query and here we will explain how to use them."
output:
  rmarkdown::html_vignette:
    df_print: kable
vignette: >
  %\VignetteIndexEntry{Customise your query}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Geocoding is surprisingly hard.
Address formats and spellings differ in and between countries; administrative areas on different levels intersect; names, numbers, and boundaries change over time — you name it.
The OpenCage API, therefore, supports about a dozen parameters to customise queries.
This vignette explains how to use the query parameters with {opencage} to get better geocoding results.

## Multiple results

Forward geocoding typically returns multiple results because many places have the same or similar names.

By default `oc_forward_df()` only returns one result:
the one defined as the best result by the OpenCage API.
To receive more results, modify the `limit` argument, which specifies the maximum number of results that should be returned.
Integer values between 1 and 100 are allowed.


```r
oc_forward_df("Berlin")
#> # A tibble: 1 x 4
#>   placename oc_lat oc_lng oc_formatted   
#>   <chr>      <dbl>  <dbl> <chr>          
#> 1 Berlin      52.5   13.4 Berlin, Germany
```

```r
oc_forward_df("Berlin", limit = 5)
#> # A tibble: 4 x 4
#>   placename oc_lat oc_lng oc_formatted                                               
#>   <chr>      <dbl>  <dbl> <chr>                                                      
#> 1 Berlin      52.5   13.4 Berlin, Germany                                            
#> 2 Berlin      44.5  -71.2 Berlin, New Hampshire, United States of America            
#> 3 Berlin      39.8  -89.9 Berlin, Sangamon County, Illinois, United States of America
#> 4 Berlin      41.6  -72.7 Berlin, Connecticut, United States of America
```


Reverse geocoding only returns [at most one result](https://opencagedata.com/api#ranking).
Therefore, `oc_reverse_df()` does not support the `limit` argument.

OpenCage may sometimes have more than one record of one place.
Duplicated records are not returned by default.
If you set the `no_dedupe` argument to `TRUE`, you will receive duplicated results when available.


```r
oc_forward_df("Berlin", limit = 5, no_dedupe = TRUE)
#> # A tibble: 5 x 4
#>   placename oc_lat oc_lng oc_formatted                                               
#>   <chr>      <dbl>  <dbl> <chr>                                                      
#> 1 Berlin      52.5   13.4 Berlin, Germany                                            
#> 2 Berlin      52.5   13.4 Berlin, Germany                                            
#> 3 Berlin      44.5  -71.2 Berlin, New Hampshire, United States of America            
#> 4 Berlin      39.8  -89.9 Berlin, Sangamon County, Illinois, United States of America
#> 5 Berlin      41.6  -72.7 Berlin, Connecticut, United States of America
```

## Better targeted results

As you can see, place names are often ambiguous.
Happily, the OpenCage API has tools to deal with this problem.
The `countrycode`, `bounds`, and `proximity` arguments can make the query more precise.
`min_confidence` lets you limit the results to those with a specified confidence score (which is not necessarily the "best" or most "relevant" result, though).
These parameters are only relevant and available for forward geocoding.

### `countrycode`

The `countrycode` parameter restricts the results to the given country.
The country code is a two letter code as defined by the [ISO 3166-1 Alpha 2](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2) standard.
E.g. "AR" for Argentina, "FR" for France, and "NZ" for the New Zealand.


```r
oc_forward_df(placename = "Paris", countrycode = "US", limit = 5)
#> # A tibble: 5 x 4
#>   placename oc_lat oc_lng oc_formatted                                         
#>   <chr>      <dbl>  <dbl> <chr>                                                
#> 1 Paris       33.7  -95.6 Paris, TX 75460, United States of America            
#> 2 Paris       38.2  -84.3 Paris, Kentucky, United States of America            
#> 3 Paris       36.3  -88.3 Paris, TN 38242, United States of America            
#> 4 Paris       39.6  -87.7 Paris, IL 61944, United States of America            
#> 5 Paris       44.3  -70.5 Paris, Oxford County, Maine, United States of America
```
Multiple countrycodes per `placename` must be wrapped in a list.
Here is an example with places called "Paris" in Italy and Portugal.


```r
oc_forward_df(placename = "Paris", countrycode = list(c("IT", "PT")), limit = 5)
#> # A tibble: 5 x 4
#>   placename oc_lat oc_lng oc_formatted                                                           
#>   <chr>      <dbl>  <dbl> <chr>                                                                  
#> 1 Paris       44.6   7.28 Brossasco, Cuneo, Italy                                                
#> 2 Paris       46.5  10.4  23032 Valfurva, Italy                                                  
#> 3 Paris       37.4  -8.79 8670-320 Odemira, Portugal                                             
#> 4 Paris       43.5  12.1  Paris, Monterchi, Arezzo, Italy                                        
#> 5 Paris       46.5  11.3  Paris, Via Firenze - Florenzstraße, 56, 39100 Bolzano - Bozen BZ, Italy
```

Despite the name, country codes also exist for territories that are not independent states, e.g. Gibraltar ("GI"), Greenland ("GL"), Guadaloupe ("GP"), or Guam ("GU").
You can look up specific country codes with the {[ISOcodes](https://cran.r-project.org/package=ISOcodes)} or {[countrycodes](https://vincentarelbundock.github.io/countrycode/)} packages or on the [ISO](https://www.iso.org/obp/ui/#search/code/) or [Wikipedia](https://en.wikipedia.org/wiki/ISO_3166-1) webpages. In fact, you can also look up country codes via OpenCage as well. If you were interested in the country code of Curaçao for example, you could run:


```r
oc_forward_df("Curaçao", no_annotations = FALSE)["oc_iso_3166_1_alpha_2"]
#> # A tibble: 1 x 1
#>   oc_iso_3166_1_alpha_2
#>   <chr>                
#> 1 CW
```

### `bounds`

The `bounds` parameter restricts the possible results to a defined [bounding box](https://wiki.openstreetmap.org/wiki/Bounding_Box).
A bounding box is a named numeric vector with four coordinates specifying its south-west and north-east corners: `(xmin, ymin, xmax, ymax)`.
The bounds parameter can most easily be specified with the `oc_bbox()` helper. For example, `bounds = oc_bbox(-0.56, 51.28, 0.27, 51.68)`.
OpenCage provides a '[bounds-finder](https://opencagedata.com/bounds-finder)' to interactively determine bounds values.

Below is an example of the use of `bounds` where the bounding box specifies the the South American continent.


```r
oc_forward_df(placename = "Paris", bounds = oc_bbox(-97, -56, -32, 12), limit = 5)
#> # A tibble: 5 x 4
#>   placename oc_lat oc_lng oc_formatted                                                          
#>   <chr>      <dbl>  <dbl> <chr>                                                                 
#> 1 Paris      -6.71  -69.9 Eirunepé, Região Geográfica Intermediária de Tefé, Brazil             
#> 2 Paris      -3.99  -79.2 110108, Loja, Ecuador                                                 
#> 3 Paris     -13.5   -62.5 Canton Motegua, Municipio Baures, Provincia de Iténez, Bolivia        
#> 4 Paris     -23.5   -47.5 Jardim Zezo Miguel, Sorocaba, Região Metropolitana de Sorocaba, Brazil
#> 5 Paris       6.31  -75.6 París, 051052 Bello, ANT, Colombia
```

Again, you can also use {opencage} to determine a bounding box for subsequent queries.
If you wanted to map the airports on the Hawaiian islands, for example, you could find the appropriate bounding box and then search for the airports:


```r
hi <- oc_forward_df(placename = "Hawaii", no_annotations = FALSE)

hi_bbox <-
  oc_bbox(
    hi$oc_southwest_lng,
    hi$oc_southwest_lat,
    hi$oc_northeast_lng,
    hi$oc_northeast_lat
  )

oc_forward_df(placename = "Airport", bounds = hi_bbox, limit = 10)
#> # A tibble: 10 x 4
#>    placename oc_lat oc_lng oc_formatted                                                                                                    
#>    <chr>      <dbl>  <dbl> <chr>                                                                                                           
#>  1 Airport     21.3  -158. Daniel K. Inouye International Airport, Aolele Street, Honolulu, HI 96820, United States of America             
#>  2 Airport     19.7  -156. Ellison Onizuka Kona International Airport at Keahole, Queen Kaahumanu Highway, Keahole, HI, United States of A~
#>  3 Airport     20.9  -156. Kahului Airport, Airport Access Road, Kahului, HI 96732-3509, United States of America                          
#>  4 Airport     19.7  -155. Hilo International Airport, Banyan Drive, Hilo, HI 96720, United States of America                              
#>  5 Airport     21.4  -158. Kaneohe Bay Marine Corp Air Station/Mokapu Point Airport, 6th Street, Honolulu County, HI 96863, United States ~
#>  6 Airport     19.6  -156. Old Kona Airport State Recreation Area, Laniakea, Kailua-Kona, Hawaii, United States of America                 
#>  7 Airport     21.2  -157. Molokai Airport, Mauna Loa Highway, Maui County, HI 96729, United States of America                             
#>  8 Airport     20.8  -157. Lanai Airport, Lanai Airport Road, Maui County, HI 96763, United States of America                              
#>  9 Airport     21.3  -158. Kalaeloa Airport, Midway Street, Kapolei, HI 96862, United States of America                                    
#> 10 Airport     21.0  -157. Kapalua-West Maui Airport, Honoapiilani Highway, Lahaina, HI 96761, United States of America
```

Note that such a query will only give you airports with "Airport" in their name or address, but not necessarily airfields, airstrips, etc. If you are more interested in these kind of features, you might want to take a look at the {[osmdata](https://docs.ropensci.org/osmdata/)} package.

### `proximity`

The `proximity` parameter provides OpenCage with a hint to bias results in favour of those closer to the specified location.
It is just one of many factors used for ranking results, however, and (some) results may be far away from the location or point passed to the `proximity` parameter.
A point is a named numeric vector of a latitude and longitude coordinate pair in decimal format.
The `proximity` parameter can most easily be specified with the `oc_points()` helper.
For example, `proximity = oc_point(38.0, -84.5)`, if you happen to already know the coordinates.
If not, you can also look them up with {opencage}, of course:


```r
lx <- oc_forward_df("Lexington, Kentucky")

lx_point <- oc_points(lx$oc_lat, lx$oc_lng)

oc_forward_df(placename = "Paris", proximity = lx_point, limit = 5)
#> # A tibble: 4 x 4
#>   placename oc_lat oc_lng oc_formatted                             
#>   <chr>      <dbl>  <dbl> <chr>                                    
#> 1 Paris       38.2 -84.3  Paris, Kentucky, United States of America
#> 2 Paris       48.9   2.35 Paris, France                            
#> 3 Paris       39.6 -87.7  Paris, IL 61944, United States of America
#> 4 Paris       38.8 -85.6  Paris, IN 47230, United States of America
```
Note that the French capital is listed before other places in the US, which are closer to the point provided.
This illustrates how `proximity` is only one of many factors influencing the ranking of results.

### Confidence

`min_confidence` — an integer value between 0 and 10 — indicates the precision of the returned result as defined by its geographical extent, i.e. by the extent of the result's bounding box.
When you specify `min_confidence`, only results with at least the requested confidence will be returned.
Thus, in the following example, the French capital is too large to be returned.


```r
oc_forward_df(placename = "Paris", min_confidence = 7, limit = 5)
#> # A tibble: 1 x 4
#>   placename oc_lat oc_lng oc_formatted                             
#>   <chr>      <dbl>  <dbl> <chr>                                    
#> 1 Paris       38.2  -84.3 Paris, Kentucky, United States of America
```

Note that confidence is not used for the [ranking of results](https://opencagedata.com/api#ranking).
It does not tell you which result is more "correct" or "relevant", nor what type of thing the result is, but rather how small a result is, geographically speaking.
See the [API documentation](https://opencagedata.com/api#confidence) for details.

## Retrieve more information from the API

Besides parameters to target your search better, OpenCage offers parameters to receive more or specific types of information from the API.

### `language`

If you would like to get your results in a specific language, you can pass an [IETF BCP 47 language tag](https://en.wikipedia.org/wiki/IETF_language_tag), such as "tr" for Turkish or "pt-BR" for Brazilian Portuguese, to the `language` parameter.
OpenCage will attempt to return results in that language.


```r
oc_forward_df(placename = "Munich", language = "tr")
#> # A tibble: 1 x 4
#>   placename oc_lat oc_lng oc_formatted           
#>   <chr>      <dbl>  <dbl> <chr>                  
#> 1 Munich      48.1   11.6 Münih, Bavyera, Almanya
```

Alternatively, you can specify the "native" tag, in which case OpenCage will attempt to return the response in the "official" language(s) of the location.
Keep in mind, however, that some countries have more than one official language or that the official language may not be the one actually used day-to-day.


```r
oc_forward_df(placename = "Munich", language = "native")
#> # A tibble: 1 x 4
#>   placename oc_lat oc_lng oc_formatted                
#>   <chr>      <dbl>  <dbl> <chr>                       
#> 1 Munich      48.1   11.6 München, Bayern, Deutschland
```

If the `language` parameter is set to `NULL` (which is the default), the tag is not recognized, or OpenCage does not have a record in that language, the results will be returned in English.


```r
oc_forward_df(placename = "München")
#> # A tibble: 1 x 4
#>   placename oc_lat oc_lng oc_formatted            
#>   <chr>      <dbl>  <dbl> <chr>                   
#> 1 München     48.1   11.6 Munich, Bavaria, Germany
```

To find the correct language tag for your desired language, you can search for the language on the [BCP47 language subtag lookup](https://r12a.github.io/app-subtags/) for example.
Note however, that there are some language tags in use on OpenStreetMap, one of OpenCage's main sources, that do not conform with the IETF BCP 47 standard.
For example, OSM uses [`zh_pinyin`](https://wiki.openstreetmap.org/w/index.php?title=Multilingual_names#China) instead of `zh-Latn-pinyin` for [Hanyu Pinyin](https://en.wikipedia.org/wiki/Pinyin).
It might, therefore, be helpful to consult the details page of the target country on openstreetmap.org to see which language tags are actually used.
In any case, neither the OpenCage API nor the functions in this package will validate the language tags you provide.

For further details, see [OpenCage's API documentation](https://opencagedata.com/api#language).

### Annotations

OpenCage supplies additional information about the result location in what it calls [annotations](https://opencagedata.com/api#annotations).
Annotations include, among a variety of other types of information, country information, time of sunset and sunrise, [UN M49](https://en.wikipedia.org/wiki/UN_M49) codes or the location in different geocoding formats, like [Maidenhead](https://en.wikipedia.org/wiki/Maidenhead_Locator_System), [Mercator projection](https://en.wikipedia.org/wiki/Mercator_projection) ([EPSG:3857](https://epsg.io/3857)), [geohash](https://en.wikipedia.org/wiki/Geohash) or [what3words](https://en.wikipedia.org/wiki/What3words).
Some annotations, like the [Irish Transverse Mercator](https://en.wikipedia.org/wiki/Irish_Transverse_Mercator) (ITM, [EPSG:2157](https://epsg.io/2157)) or the [Federal Information Processing Standards](https://en.wikipedia.org/wiki/Federal_Information_Processing_Standards) (FIPS) code will only be shown when appropriate.

Whether the annotations are shown, is controlled by the `no_annotations` argument.
It is `TRUE` by default, which means that the output will _not_ contain annotations.
(Yes, inverted argument names are confusing, but we just follow OpenCage's lead here.)
When you set `no_annotations` to `FALSE`, all columns are returned (i.e. `output` is implicitly set to `"all"`).
This leads to a results with a lot of columns.


```r
oc_forward_df("Dublin", no_annotations = FALSE)
#> # A tibble: 1 x 71
#>   placename oc_lat oc_lng oc_confidence oc_formatted oc_mgrs oc_maidenhead oc_callingcode oc_flag oc_geohash oc_qibla oc_wikidata
#>   <chr>      <dbl>  <dbl>         <int> <chr>        <chr>   <chr>                  <int> <chr>   <chr>         <dbl> <chr>      
#> 1 Dublin      53.3  -6.26             5 Dublin, Ire~ 29UPV8~ IO63ui83sw               353 "\U000~ gc7x9812v~     114. Q1761      
#> # ... with 59 more variables: oc_dms_lat <chr>, oc_dms_lng <chr>, oc_itm_easting <chr>, oc_itm_northing <chr>, oc_mercator_x <dbl>,
#> #   oc_mercator_y <dbl>, oc_osm_edit_url <chr>, oc_osm_note_url <chr>, oc_osm_url <chr>, oc_un_m49_statistical_groupings <list>,
#> #   oc_un_m49_regions_europe <chr>, oc_un_m49_regions_ie <chr>, oc_un_m49_regions_northern_europe <chr>, oc_un_m49_regions_world <chr>,
#> #   oc_currency_alternate_symbols <list>, oc_currency_decimal_mark <chr>, oc_currency_html_entity <chr>, oc_currency_iso_code <chr>,
#> #   oc_currency_iso_numeric <chr>, oc_currency_name <chr>, oc_currency_smallest_denomination <int>, oc_currency_subunit <chr>,
#> #   oc_currency_subunit_to_unit <int>, oc_currency_symbol <chr>, oc_currency_symbol_first <int>, oc_currency_thousands_separator <chr>,
#> #   oc_roadinfo_drive_on <chr>, oc_roadinfo_speed_in <chr>, oc_sun_rise_apparent <int>, oc_sun_rise_astronomical <int>,
#> #   oc_sun_rise_civil <int>, oc_sun_rise_nautical <int>, oc_sun_set_apparent <int>, oc_sun_set_astronomical <int>, oc_sun_set_civil <int>,
#> #   oc_sun_set_nautical <int>, oc_timezone_name <chr>, oc_timezone_now_in_dst <int>, oc_timezone_offset_sec <int>,
#> #   oc_timezone_offset_string <chr>, oc_timezone_short_name <chr>, oc_what3words_words <chr>, oc_northeast_lat <dbl>,
#> #   oc_northeast_lng <dbl>, oc_southwest_lat <dbl>, oc_southwest_lng <dbl>, oc_iso_3166_1_alpha_2 <chr>, oc_iso_3166_1_alpha_3 <chr>,
#> #   oc_category <chr>, oc_type <chr>, oc_city <chr>, oc_continent <chr>, oc_country <chr>, oc_country_code <chr>, oc_county <chr>,
#> #   oc_county_code <chr>, oc_political_union <chr>, oc_state <chr>, oc_state_code <chr>
```

### `roadinfo`

`roadinfo` indicates whether the geocoder should attempt to match the nearest road (rather than an address) and provide additional road and driving information.
It is `FALSE` by default, which means OpenCage will not attempt to match the nearest road.
Some road and driving information is nevertheless provided as part of the annotations (see above), even when `roadinfo` is set to `FALSE`.


```r
oc_forward_df(placename = c("Europa Advance Rd", "Bovoni Rd"), roadinfo = TRUE)
#> # A tibble: 2 x 29
#>   placename oc_lat oc_lng oc_confidence oc_formatted oc_roadinfo_dri~ oc_roadinfo_one~ oc_roadinfo_road oc_roadinfo_roa~ oc_roadinfo_spe~
#>   <chr>      <dbl>  <dbl>         <int> <chr>        <chr>            <chr>            <chr>            <chr>            <chr>           
#> 1 Europa A~   36.1  -5.34             9 Europa Adva~ right            yes              Europa Advance ~ secondary        km/h            
#> 2 Bovoni Rd   18.3 -64.9              9 Bovoni Hill~ left             <NA>             <NA>             <NA>             mph             
#> # ... with 19 more variables: oc_roadinfo_surface <chr>, oc_northeast_lat <dbl>, oc_northeast_lng <dbl>, oc_southwest_lat <dbl>,
#> #   oc_southwest_lng <dbl>, oc_iso_3166_1_alpha_2 <chr>, oc_iso_3166_1_alpha_3 <chr>, oc_category <chr>, oc_type <chr>,
#> #   oc_continent <chr>, oc_country <chr>, oc_country_code <chr>, oc_postcode <chr>, oc_road <chr>, oc_road_type <chr>, oc_state <chr>,
#> #   oc_county <chr>, oc_peak <chr>, oc_state_code <chr>
```

A [blog post](https://blog.opencagedata.com/post/new-optional-parameter-roadinfo) provides more details.

### Abbreviated addresses

The geocoding functions also have an `abbr` parameter, which is `FALSE` by default.
When it is `TRUE`, the addresses in the `formatted` field of the results are abbreviated (e.g. "Main St." instead of "Main Street").


```r
oc_forward_df("Wall Street")
#> # A tibble: 1 x 4
#>   placename   oc_lat oc_lng oc_formatted                                             
#>   <chr>        <dbl>  <dbl> <chr>                                                    
#> 1 Wall Street   40.7  -74.0 Wall Street, New York, NY 10005, United States of America
oc_forward_df("Wall Street", abbrv = TRUE)
#> # A tibble: 1 x 4
#>   placename   oc_lat oc_lng oc_formatted                    
#>   <chr>        <dbl>  <dbl> <chr>                           
#> 1 Wall Street   40.7  -74.0 Wall St, New York, NY 10005, USA
```

See [this blog post](https://blog.opencagedata.com/post/160294347883/shrtr-pls) for more information.

## Vectorised arguments

All of the function arguments mentioned above are vectorised, so you can send queries like this:


```r
oc_forward_df(
  placename = c("New York", "Rio", "Tokyo"),
  language = c("es", "de", "fr")
)
#> # A tibble: 3 x 4
#>   placename oc_lat oc_lng oc_formatted                                                     
#>   <chr>      <dbl>  <dbl> <chr>                                                            
#> 1 New York    40.7  -74.0 Nueva York, Estados Unidos de América                            
#> 2 Rio        -22.9  -43.2 Rio de Janeiro, Região Metropolitana do Rio de Janeiro, Brasilien
#> 3 Tokyo       35.7  140.  Tokyo, Japon
```

Or geocode place names with country codes in a data frame:


```r
for_df <-
  data.frame(
    location = c("Golden Gate Bridge", "Buckingham Palace", "Eiffel Tower"),
    ccode = c("at", "cg", "be")
  )

oc_forward_df(for_df, placename = location, countrycode = ccode)
#> # A tibble: 3 x 5
#>   location           ccode oc_lat oc_lng oc_formatted                                                                              
#>   <chr>              <chr>  <dbl>  <dbl> <chr>                                                                                     
#> 1 Golden Gate Bridge at     48.0   15.6  Karer, Golden Gate Bridge, 3180 Gemeinde Lilienfeld, Austria                              
#> 2 Buckingham Palace  cg     -4.80  11.8  Buckingham Palace, Boulevard du Général Charles de Gaulle, Pointe-Noire, Congo-Brazzaville
#> 3 Eiffel Tower       be     50.9    4.34 Eiffel Tower, Avenue de Bouchout - Boechoutlaan, 1020 City of Brussels, Belgium
```

This also works with `oc_reverse_df()`, of course.


```r
rev_df <-
  data.frame(
    lat = c(51.952659, 41.401372),
    lon = c(7.632473, 2.128685)
  )

oc_reverse_df(rev_df, lat, lon, language = "native")
#> # A tibble: 2 x 3
#>     lat   lon oc_formatted                                        
#>   <dbl> <dbl> <chr>                                               
#> 1  52.0  7.63 Friedrich-Ebert-Straße 7, 48153 Münster, Deutschland
#> 2  41.4  2.13 Carrer de Calatrava, 68, 08017 Barcelona, España
```

## Further information

For further information about the output and query parameters, see the [OpenCage API docs](https://opencagedata.com/api) and the [OpenCage FAQ](https://opencagedata.com/faq).
When building queries, OpenCage's [best practices](https://opencagedata.com/api#bestpractices) can be very useful, as well as their guide to [geocoding accuracy](https://opencagedata.com/guides/how-to-think-about-geocoding-accuracy).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{opencage_reverse}
\alias{opencage_reverse}
\title{Reverse geocoding}
\usage{
opencage_reverse(
  latitude,
  longitude,
  key = opencage_key(),
  bounds = NULL,
  countrycode = NULL,
  language = NULL,
  limit = 10,
  min_confidence = NULL,
  no_annotations = FALSE,
  no_dedupe = FALSE,
  no_record = FALSE,
  abbrv = FALSE,
  add_request = TRUE
)
}
\arguments{
\item{latitude}{Numeric vectors of latitude and longitude values.}

\item{longitude}{Numeric vectors of latitude and longitude values.}

\item{key}{Your OpenCage API key as a character vector of length one. By
default, \code{\link[=opencage_key]{opencage_key()}} will attempt to retrieve the key from the
environment variable \code{OPENCAGE_KEY}.}

\item{bounds}{Bounding box, ignored for reverse geocoding.}

\item{countrycode}{Country code, ignored for reverse geocoding.}

\item{language}{An \href{https://en.wikipedia.org/wiki/IETF_language_tag}{IETF BCP 47 language tag} (such as "es" for
Spanish or "pt-BR" for Brazilian Portuguese). OpenCage will attempt to
return results in that language. Alternatively you can specify the "native"
tag, in which case OpenCage will attempt to return the response in the
"official" language(s). In case the \code{language} parameter is set to \code{NULL}
(which is the default), the tag is not recognized, or OpenCage does not
have a record in that language, the results will be returned in English.}

\item{limit}{How many results should be returned (1-100), ignored for reverse
geocoding.}

\item{min_confidence}{Numeric vector of integer values between 0 and 10
indicating the precision of the returned result as defined by its
geographical extent, (i.e. by the extent of the result's bounding box). See
the \href{https://opencagedata.com/api#confidence}{API documentation} for
details. Only results with at least the requested confidence will be
returned. Default is \code{NULL}.}

\item{no_annotations}{Logical vector indicating whether additional
information about the result location should be returned. \code{TRUE} by
default, which means that the results will not contain annotations.}

\item{no_dedupe}{Logical vector (default \code{FALSE}), when \code{TRUE} the results
will not be deduplicated.}

\item{no_record}{Logical vector of length one (default \code{FALSE}), when \code{TRUE}
no log entry of the query is created, and the geocoding request is not
cached by OpenCage.}

\item{abbrv}{Logical vector (default \code{FALSE}), when \code{TRUE} addresses in the
\code{formatted} field of the results are abbreviated (e.g. "Main St." instead
of "Main Street").}

\item{add_request}{Logical vector (default \code{FALSE}) indicating whether the
request is returned again with the results. If the \code{return} value is a
\code{df_list}, the query text is added as a column to the results. \code{json_list}
results will contain all request parameters, including the API key used!
This is currently ignored by OpenCage if return value is \code{geojson_list}.}
}
\value{
A list with
\itemize{
\item results as a tibble with one line per result,
\item the number of results as an integer,
\item the timestamp as a POSIXct object,
\item rate_info tibble/data.frame with the maximal number of API calls  per
day for the used key, the number of remaining calls for the day and the time
at which the number of remaining calls will be reset.
}
}
\description{
\ifelse{html}{\figure{lifecycle-soft-deprecated.svg}{options: alt='Soft-deprecated lifecycle'}}{\strong{Soft-deprecated}}

Soft deprecated: use \code{oc_reverse} or \code{oc_reverse_df} for reverse geocoding.
}
\examples{
\dontshow{if (oc_key_present() && oc_api_ok()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

opencage_reverse(
  latitude = 0, longitude = 0,
  limit = 2
)
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_api_ok.R
\name{oc_api_ok}
\alias{oc_api_ok}
\title{Is the OpenCage API available?}
\usage{
oc_api_ok(url = "https://api.opencagedata.com")
}
\arguments{
\item{url}{The URL of the OpenCage API, \url{https://api.opencagedata.com}.}
}
\value{
A single logical value, \code{TRUE} or \code{FALSE}.
}
\description{
Checks whether the OpenCage API can be reached.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_bbox.R
\name{oc_bbox}
\alias{oc_bbox}
\alias{oc_bbox.numeric}
\alias{oc_bbox.data.frame}
\alias{oc_bbox.bbox}
\title{List of bounding boxes for OpenCage queries}
\usage{
oc_bbox(...)

\method{oc_bbox}{numeric}(xmin, ymin, xmax, ymax, ...)

\method{oc_bbox}{data.frame}(data, xmin, ymin, xmax, ymax, ...)

\method{oc_bbox}{bbox}(bbox, ...)
}
\arguments{
\item{...}{Ignored.}

\item{xmin}{Minimum longitude (also known as \code{min_lon}, \code{southwest_lng},
\code{west}, or \code{left}).}

\item{ymin}{Minimum latitude (also known as \code{min_lat}, \code{southwest_lat},
\code{south}, or \code{bottom}).}

\item{xmax}{Maximum longitude (also known as \code{max_lon}, \code{northeast_lng},
\code{east}, or \code{right}).}

\item{ymax}{Maximum latitude (also known as \code{max_lat}, \code{northeast_lat},
\code{north}, or \code{top}).}

\item{data}{A \code{data.frame} containing at least 4 columns with \code{xmin}, \code{ymin},
\code{xmax}, and \code{ymax} values, respectively.}

\item{bbox}{A \code{bbox} object, see \code{sf::st_bbox}.}
}
\value{
A list of bounding boxes, each of class \code{bbox}.
}
\description{
Create a list of bounding boxes for OpenCage queries.
}
\examples{
oc_bbox(-5.63160, 51.280430, 0.278970, 51.683979)

xdf <-
  data.frame(
    place = c("Hamburg", "Hamburg"),
    northeast_lat = c(54.0276817, 42.7397729),
    northeast_lng = c(10.3252805, -78.812825),
    southwest_lat = c(53.3951118, 42.7091669),
    southwest_lng = c(8.1053284, -78.860521)
  )
oc_bbox(
  xdf,
  southwest_lng,
  southwest_lat,
  northeast_lng,
  northeast_lat
)

# create bbox list column with dplyr
library(dplyr)
xdf \%>\%
  mutate(
    bbox =
      oc_bbox(
        southwest_lng,
        southwest_lat,
        northeast_lng,
        northeast_lat
      )
  )

# create bbox list from a simple features bbox
if (requireNamespace("sf", quietly = TRUE)) {
  library(sf)
  bbox <- st_bbox(c(xmin = 16.1, xmax = 16.6, ymax = 48.6, ymin = 47.9),
    crs = 4326
  )
  oc_bbox(bbox)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{opencage_key}
\alias{opencage_key}
\title{Retrieve Opencage API key}
\usage{
opencage_key(quiet = TRUE)
}
\arguments{
\item{quiet}{Logical vector of length one indicating whether the key is
returned quietly or whether a message is printed.}
}
\description{
\ifelse{html}{\figure{lifecycle-soft-deprecated.svg}{options: alt='Soft-deprecated lifecycle'}}{\strong{Soft-deprecated}}

Soft-deprecated and will be removed from the package together with
\code{opencage_forward()} and \code{opencage_reverse()}.

Retrieves the OpenCage API Key from the environment variable \code{OPENCAGE_KEY}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_doc.R
\docType{data}
\name{countrycodes}
\alias{countrycodes}
\title{Country codes}
\format{
All possible ISO 3166-1 Alpha 2 standard country codes.
}
\description{
Country codes
}
\examples{
data("countrycodes")
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/memoise.R
\name{oc_clear_cache}
\alias{oc_clear_cache}
\title{Clear the opencage cache}
\usage{
oc_clear_cache()
}
\description{
Forget past results and reset the \pkg{opencage} cache.
}
\examples{
\dontshow{if (oc_key_present() && oc_api_ok()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

system.time(oc_reverse(latitude = 10, longitude = 10))
system.time(oc_reverse(latitude = 10, longitude = 10))
oc_clear_cache()
system.time(oc_reverse(latitude = 10, longitude = 10))
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{opencage_forward}
\alias{opencage_forward}
\title{Forward geocoding}
\usage{
opencage_forward(
  placename,
  key = opencage_key(),
  bounds = NULL,
  countrycode = NULL,
  language = NULL,
  limit = 10L,
  min_confidence = NULL,
  no_annotations = FALSE,
  no_dedupe = FALSE,
  no_record = FALSE,
  abbrv = FALSE,
  add_request = TRUE
)
}
\arguments{
\item{placename}{A character vector with the location names or addresses to
be geocoded.

If the locations are addresses, see \href{https://github.com/OpenCageData/opencagedata-misc-docs/blob/master/query-formatting.md}{OpenCage's instructions}
on how to format addresses for best forward geocoding results.}

\item{key}{Your OpenCage API key as a character vector of length one. By
default, \code{\link[=opencage_key]{opencage_key()}} will attempt to retrieve the key from the
environment variable \code{OPENCAGE_KEY}.}

\item{bounds}{A list of bounding boxes of length one or \code{length(placename)}.
Bounding boxes are named numeric vectors, each with four coordinates
forming the south-west and north-east corners of the bounding box:
\code{list(c(xmin, ymin, xmax, ymax))}. \code{bounds} restricts the possible results
to the supplied region. It can be specified with the \code{\link[=oc_bbox]{oc_bbox()}} helper.
For example: \code{bounds = oc_bbox(-0.563160, 51.280430, 0.278970, 51.683979)}.
Default is \code{NULL}.}

\item{countrycode}{A two letter code as defined by the \href{https://www.iso.org/obp/ui/#search/code}{ISO 3166-1 Alpha 2} standard that restricts the
results to the given country or countries. E.g. "AR" for Argentina, "FR"
for France, "NZ" for the New Zealand. Multiple countrycodes per \code{placename}
must be wrapped in a list. Default is \code{NULL}.}

\item{language}{An \href{https://en.wikipedia.org/wiki/IETF_language_tag}{IETF BCP 47 language tag} (such as "es" for
Spanish or "pt-BR" for Brazilian Portuguese). OpenCage will attempt to
return results in that language. Alternatively you can specify the "native"
tag, in which case OpenCage will attempt to return the response in the
"official" language(s). In case the \code{language} parameter is set to \code{NULL}
(which is the default), the tag is not recognized, or OpenCage does not
have a record in that language, the results will be returned in English.}

\item{limit}{Numeric vector of integer values to determine the maximum number
of results returned for each \code{placename}. Integer values between 1 and 100
are allowed. Default is 10.}

\item{min_confidence}{Numeric vector of integer values between 0 and 10
indicating the precision of the returned result as defined by its
geographical extent, (i.e. by the extent of the result's bounding box). See
the \href{https://opencagedata.com/api#confidence}{API documentation} for
details. Only results with at least the requested confidence will be
returned. Default is \code{NULL}.}

\item{no_annotations}{Logical vector indicating whether additional
information about the result location should be returned. \code{TRUE} by
default, which means that the results will not contain annotations.}

\item{no_dedupe}{Logical vector (default \code{FALSE}), when \code{TRUE} the results
will not be deduplicated.}

\item{no_record}{Logical vector of length one (default \code{FALSE}), when \code{TRUE}
no log entry of the query is created, and the geocoding request is not
cached by OpenCage.}

\item{abbrv}{Logical vector (default \code{FALSE}), when \code{TRUE} addresses in the
\code{formatted} field of the results are abbreviated (e.g. "Main St." instead
of "Main Street").}

\item{add_request}{Logical vector (default \code{FALSE}) indicating whether the
request is returned again with the results. If the \code{return} value is a
\code{df_list}, the query text is added as a column to the results. \code{json_list}
results will contain all request parameters, including the API key used!
This is currently ignored by OpenCage if return value is \code{geojson_list}.}
}
\value{
A list with
\itemize{
\item results as a tibble with one line per result,
\item the number of results as an integer,
\item the timestamp as a POSIXct object,
\item rate_info tibble/data.frame with the maximal number of API calls  per
day for the used key, the number of remaining calls for the day and the time
at which the number of remaining calls will be reset.
}
}
\description{
\ifelse{html}{\figure{lifecycle-soft-deprecated.svg}{options: alt='Soft-deprecated lifecycle'}}{\strong{Soft-deprecated}}

Soft deprecated: use \code{oc_forward} or \code{oc_forward_df} for forward geocoding.
}
\examples{
\dontshow{if (oc_key_present() && oc_api_ok()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
opencage_forward(placename = "Sarzeau")
opencage_forward(placename = "Islington, London")
opencage_forward(placename = "Triererstr 15,
                              Weimar 99423,
                              Deutschland")
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_reverse.R
\name{oc_reverse_df}
\alias{oc_reverse_df}
\alias{oc_reverse_df.data.frame}
\alias{oc_reverse_df.numeric}
\title{Reverse geocoding with data frames}
\usage{
oc_reverse_df(...)

\method{oc_reverse_df}{data.frame}(
  data,
  latitude,
  longitude,
  bind_cols = TRUE,
  output = c("short", "all"),
  language = NULL,
  min_confidence = NULL,
  roadinfo = FALSE,
  no_annotations = TRUE,
  no_dedupe = FALSE,
  abbrv = FALSE,
  ...
)

\method{oc_reverse_df}{numeric}(
  latitude,
  longitude,
  output = c("short", "all"),
  language = NULL,
  min_confidence = NULL,
  no_annotations = TRUE,
  no_dedupe = FALSE,
  abbrv = FALSE,
  ...
)
}
\arguments{
\item{...}{Ignored.}

\item{data}{A data frame.}

\item{latitude, longitude}{Unquoted variable names of numeric columns or
vectors of latitude and longitude values.}

\item{bind_cols}{When \code{bind_col = TRUE}, the default, the results are column
bound to \code{data}. When \code{FALSE}, the results are returned as a new tibble.}

\item{output}{A character vector of length one indicating whether only the
formatted address (\code{"short"}, the default) or all variables (\code{"all"})
variables should be returned.}

\item{language}{Character vector, or an unquoted variable name of such a
vector, of \href{https://en.wikipedia.org/wiki/IETF_language_tag}{IETF BCP 47 language tags} (such as "es" for
Spanish or "pt-BR" for Brazilian Portuguese). OpenCage will attempt to
return results in that language. Alternatively you can specify the "native"
tag, in which case OpenCage will attempt to return the response in the
"official" language(s). In case the \code{language} parameter is set to \code{NULL}
(which is the default), the tag is not recognized, or OpenCage does not
have a record in that language, the results will be returned in English.}

\item{min_confidence}{Numeric vector of integer values, or an unquoted
variable name of such a vector, between 0 and 10 indicating the precision
of the returned result as defined by its geographical extent, (i.e. by the
extent of the result's bounding box). See the \href{https://opencagedata.com/api#confidence}{API documentation} for details. Only
results with at least the requested confidence will be returned. Default is
\code{NULL}).}

\item{roadinfo}{Logical vector, or an unquoted variable name of such a
vector, indicating whether the geocoder should attempt to match the nearest
road (rather than an address) and provide additional road and driving
information. Default is \code{FALSE}.}

\item{no_annotations}{Logical vector, or an unquoted variable name of such a
vector, indicating whether additional information about the result location
should be returned. \code{TRUE} by default, which means that the results will
not contain annotations.}

\item{no_dedupe}{Logical vector, or an unquoted variable name of such a
vector. Default is \code{FALSE}. When \code{TRUE} the results will not be
deduplicated.}

\item{abbrv}{Logical vector, or an unquoted variable name of such a vector.
Default is \code{FALSE}. When \code{TRUE} addresses in the \code{oc_formatted} variable of
the results are abbreviated (e.g. "Main St." instead of "Main Street").}
}
\value{
A tibble. Column names coming from the OpenCage API are prefixed with
\code{"oc_"}.
}
\description{
Reverse geocoding from latitude and longitude pairs to the names and
addresses of a location.
}
\examples{
\dontshow{if (oc_key_present() && oc_api_ok()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

library(tibble)
df <- tibble(id = 1:4,
             lat = c(-36.85007, 47.21864, 53.55034, 34.05369),
             lng = c(174.7706, -1.554136, 10.000654, -118.242767))

# Return formatted address of lat/lng values
oc_reverse_df(df, latitude = lat, longitude = lng)

# Return more detailed information about the locations
oc_reverse_df(df, latitude = lat, longitude = lng,
              output = "all")

# Return results in a preferred language if possible
oc_reverse_df(df, latitude = lat, longitude = lng,
              language = "fr")

# oc_reverse_df accepts unquoted column names for all
# arguments except bind_cols and output.
# This makes it possible to build up more detailed queries
# through the data frame passed to the data argument.

df2 <- add_column(df,
                  language = c("en", "fr", "de", "en"),
                  confidence = c(8, 10, 10, 10))

# Use language column to specify preferred language of results
# and confidence column to allow different confidence levels
oc_reverse_df(df2, latitude = lat, longitude = lng,
              language = language,
              min_confidence = confidence)
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[=oc_reverse]{oc_reverse()}} for inputs as vectors, or \code{\link[=oc_forward]{oc_forward()}} and
\code{\link[=oc_forward]{oc_forward()}} for forward geocoding. For more information about the API
and the various parameters, see the \href{https://opencagedata.com/api}{OpenCage API documentation}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opencage-package.R
\docType{package}
\name{opencage}
\alias{opencage}
\alias{opencage-package}
\title{opencage: Geocode with the OpenCage API.}
\description{
Geocode with the OpenCage API, either from place name to longitude and
latitude (forward geocoding) or from longitude and latitude to the name and
address of a location (reverse geocoding), see \url{https://opencagedata.com/}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/opencage/}
  \item \url{https://github.com/ropensci/opencage}
  \item Report bugs at \url{https://github.com/ropensci/opencage/issues}
}

}
\author{
\strong{Maintainer}: Daniel Possenriede \email{possenriede+r@gmail.com} (\href{https://orcid.org/0000-0002-6738-9845}{ORCID})

Authors:
\itemize{
  \item Jesse Sadler (\href{https://orcid.org/0000-0001-6081-9681}{ORCID})
  \item Maëlle Salmon (\href{https://orcid.org/0000-0002-2815-0399}{ORCID})
}

Other contributors:
\itemize{
  \item Noam Ross [contributor]
  \item Jake Russ [contributor]
  \item Julia Silge (Julia Silge reviewed the package for rOpenSci, see <https://github.com/ropensci/onboarding/issues/36>.) [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_config.R
\name{oc_config}
\alias{oc_config}
\title{Configure settings}
\usage{
oc_config(
  key = Sys.getenv("OPENCAGE_KEY"),
  rate_sec = getOption("oc_rate_sec", default = 1L),
  no_record = getOption("oc_no_record", default = TRUE),
  show_key = getOption("oc_show_key", default = FALSE),
  ...
)
}
\arguments{
\item{key}{Your OpenCage API key as a character vector of length one. Do not
pass the key directly as a parameter, though. See details.}

\item{rate_sec}{Numeric vector of length one. Sets the maximum number of
requests sent to the OpenCage API per second. Defaults to the value set in
the \code{oc_rate_sec} option, or, in case that does not exist, to 1L.}

\item{no_record}{Logical vector of length one. When \code{TRUE}, OpenCage will not
create log entries of the queries and will not cache the geocoding
requests. Defaults to the value set in the \code{oc_no_record} option, or, in
case that does not exist, to \code{TRUE}.}

\item{show_key}{Logical vector of length one. This is only relevant when
debugging \code{oc_forward()} or \code{oc_reverse()} calls with the \code{return = "url_only"} argument. When \code{TRUE}, the result will show your OpenCage API
key in the URL as stored in the \code{OPENCAGE_KEY} environment variable. When
not \code{TRUE}, the API key will be replaced with the string \code{OPENCAGE_KEY}.
\code{show_key} defaults to the value set in the \code{oc_show_key} option, or, in
case that does not exist, to \code{FALSE}.}

\item{...}{Ignored.}
}
\description{
Configure session settings for \pkg{opencage}.
}
\section{Set your OpenCage API key}{


\pkg{opencage} will conveniently retrieve your API key if it is saved in the
environment variable \code{"OPENCAGE_KEY"}. \code{\link[=oc_config]{oc_config()}} will help to set that
environment variable. Do not pass the key directly as a parameter to the
function, though, as you risk exposing it via your script or your history.
There are three safer ways to set your API key instead:
\enumerate{
\item Save your API key as an environment variable in
\code{\link[base:Startup]{.Renviron}} as described in \href{https://rstats.wtf/r-startup.html#renviron}{What They Forgot to Teach You About R} or \href{https://csgillespie.github.io/efficientR/set-up.html#renviron}{Efficient R Programming}.
From there it will be fetched by all functions that call the OpenCage API.
You do not even have to call \code{oc_config()} to set your key; you can start
geocoding right away. If you have the \pkg{usethis} package installed, you
can edit your \code{\link[base:Startup]{.Renviron}} with \code{usethis::edit_r_environ()}.
We strongly recommend storing your API key in the user-level .Renviron, as
opposed to the project-level. This makes it less likely you will share
sensitive information by mistake.
\item If you use a package like \pkg{keyring} to store your credentials, you can
safely pass your key in a script with a function call like this
\code{oc_config(key = keyring::key_get("opencage"))}.
\item If you call \code{oc_config()} in an \code{\link[base:interactive]{base::interactive()}} session and the
\code{OPENCAGE_KEY} environment variable is not set, it will prompt you to enter
the key in the console.
}
}

\section{Set your OpenCage API rate limit}{


The rate limit allowed by the API depends on the OpenCage plan you purchased
and ranges from 1 request/sec for the "Free Trial" plan to 15 requests/sec
for the "Medium" or "Large" plans, see \url{https://opencagedata.com/pricing} for
details and up-to-date information. You can set the rate limit persistently
across sessions by setting an \code{oc_rate_sec} \link[base:options]{option} in your
\code{\link[base:Startup]{.Rprofile}}. If you have the \pkg{usethis} package
installed, you can edit your \code{\link[base:Startup]{.Rprofile}} with
\code{usethis::edit_r_profile()}.
}

\section{Prevent query logging and caching}{


By default, OpenCage will store your queries in its server logs and will
cache the forward geocoding requests on their side. They do this in order to
speed up response times and to be able to debug errors and improve their
service. Logs are automatically deleted after six months according to
OpenCage's \href{https://opencagedata.com/gdpr}{page on data protection and GDPR}.

If you set \code{no_record} to \code{TRUE}, the query contents are not logged nor
cached. OpenCage still records that you made a request, but not the specific
query you made. \code{oc_config(no_record = TRUE)} sets the \code{oc_no_record}
\link[base:options]{option} for the active R session, so it will be used for all
subsequent OpenCage queries. You can set the \code{oc_no_record}
\link[base:options]{option} persistently across sessions in your
\code{\link[base:Startup]{.Rprofile}}.

For increased privacy \pkg{opencage} sets \code{no_record} to \code{TRUE}, by default.
Please note, however, that \pkg{opencage} always caches the data it receives
from the OpenCage API locally, but only for as long as your R session is
alive.

For more information on OpenCage's policies on privacy and data protection
see \href{https://opencagedata.com/faq#legal}{their FAQs}, their \href{https://opencagedata.com/gdpr}{GDPR page}, and, for the \code{no_record} parameter, see
the relevant \href{https://blog.opencagedata.com/post/145602604628/more-privacy-with-norecord-parameter}{blog post}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_points.R
\name{oc_points}
\alias{oc_points}
\alias{oc_points.numeric}
\alias{oc_points.data.frame}
\title{List of points for OpenCage queries}
\usage{
oc_points(...)

\method{oc_points}{numeric}(latitude, longitude, ...)

\method{oc_points}{data.frame}(data, latitude, longitude, ...)
}
\arguments{
\item{...}{Ignored.}

\item{latitude, longitude}{Numeric vectors of latitude and longitude values.}

\item{data}{A \code{data.frame} containing at least 2 columns with \code{latitude} and
\code{longitude} values.}
}
\value{
A list of points. Each point is a named vector of length 2 containing
a latitude/longitude coordinate pair.
}
\description{
Create a list of points (latitude/longitude coordinate pairs) for OpenCage
queries.
}
\examples{
oc_points(-21.01404, 55.26077)

xdf <-
  data.frame(
    place = c("Hamburg", "Los Angeles"),
    lat = c(53.5503, 34.0536),
    lon = c(10.0006, -118.2427)
  )
oc_points(
  data = xdf,
  latitude = lat,
  longitude = lon
)

# create a list column with points with dplyr
library(dplyr)
xdf \%>\%
  mutate(
    points =
      oc_points(
        lat,
        lon
      )
  )

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_forward.R
\name{oc_forward}
\alias{oc_forward}
\title{Forward geocoding}
\usage{
oc_forward(
  placename,
  return = c("df_list", "json_list", "geojson_list", "url_only"),
  bounds = NULL,
  proximity = NULL,
  countrycode = NULL,
  language = NULL,
  limit = 10L,
  min_confidence = NULL,
  no_annotations = TRUE,
  roadinfo = FALSE,
  no_dedupe = FALSE,
  abbrv = FALSE,
  add_request = FALSE,
  ...
)
}
\arguments{
\item{placename}{A character vector with the location names or addresses to
be geocoded.

If the locations are addresses, see \href{https://github.com/OpenCageData/opencagedata-misc-docs/blob/master/query-formatting.md}{OpenCage's instructions}
on how to format addresses for best forward geocoding results.}

\item{return}{A character vector of length one indicating the return value of
the function, either a list of tibbles (\code{df_list}, the default), a JSON
list (\code{json_list}), a GeoJSON list (\code{geojson_list}), or the URL with which
the API would be called (\code{url_only}).}

\item{bounds}{A list of bounding boxes of length one or \code{length(placename)}.
Bounding boxes are named numeric vectors, each with four coordinates
forming the south-west and north-east corners of the bounding box:
\code{list(c(xmin, ymin, xmax, ymax))}. \code{bounds} restricts the possible results
to the supplied region. It can be specified with the \code{\link[=oc_bbox]{oc_bbox()}} helper.
For example: \code{bounds = oc_bbox(-0.563160, 51.280430, 0.278970, 51.683979)}.
Default is \code{NULL}.}

\item{proximity}{A list of points of length one or \code{length(placename)}. A
point is a named numeric vector of a latitude, longitude coordinate pair in
decimal format. \code{proximity} provides OpenCage with a hint to bias results
in favour of those closer to the specified location. It can be specified
with the \code{\link[=oc_points]{oc_points()}} helper. For example: \code{proximity = oc_points(51.9526, 7.6324)}. Default is \code{NULL}.}

\item{countrycode}{A two letter code as defined by the \href{https://www.iso.org/obp/ui/#search/code}{ISO 3166-1 Alpha 2} standard that restricts the
results to the given country or countries. E.g. "AR" for Argentina, "FR"
for France, "NZ" for the New Zealand. Multiple countrycodes per \code{placename}
must be wrapped in a list. Default is \code{NULL}.}

\item{language}{An \href{https://en.wikipedia.org/wiki/IETF_language_tag}{IETF BCP 47 language tag} (such as "es" for
Spanish or "pt-BR" for Brazilian Portuguese). OpenCage will attempt to
return results in that language. Alternatively you can specify the "native"
tag, in which case OpenCage will attempt to return the response in the
"official" language(s). In case the \code{language} parameter is set to \code{NULL}
(which is the default), the tag is not recognized, or OpenCage does not
have a record in that language, the results will be returned in English.}

\item{limit}{Numeric vector of integer values to determine the maximum number
of results returned for each \code{placename}. Integer values between 1 and 100
are allowed. Default is 10.}

\item{min_confidence}{Numeric vector of integer values between 0 and 10
indicating the precision of the returned result as defined by its
geographical extent, (i.e. by the extent of the result's bounding box). See
the \href{https://opencagedata.com/api#confidence}{API documentation} for
details. Only results with at least the requested confidence will be
returned. Default is \code{NULL}.}

\item{no_annotations}{Logical vector indicating whether additional
information about the result location should be returned. \code{TRUE} by
default, which means that the results will not contain annotations.}

\item{roadinfo}{Logical vector indicating whether the geocoder should attempt
to match the nearest road (rather than an address) and provide additional
road and driving information. Default is \code{FALSE}.}

\item{no_dedupe}{Logical vector (default \code{FALSE}), when \code{TRUE} the results
will not be deduplicated.}

\item{abbrv}{Logical vector (default \code{FALSE}), when \code{TRUE} addresses in the
\code{formatted} field of the results are abbreviated (e.g. "Main St." instead
of "Main Street").}

\item{add_request}{Logical vector (default \code{FALSE}) indicating whether the
request is returned again with the results. If the \code{return} value is a
\code{df_list}, the query text is added as a column to the results. \code{json_list}
results will contain all request parameters, including the API key used!
This is currently ignored by OpenCage if return value is \code{geojson_list}.}

\item{...}{Ignored.}
}
\value{
Depending on the \code{return} argument, \code{oc_forward} returns a list with
either
\itemize{
\item the results as tibbles (\code{"df_list"}, the default),
\item the results as JSON specified as a list (\code{"json_list"}),
\item the results as GeoJSON specified as a list (\code{"geojson_list"}),
or
\item the URL of the OpenCage API call for debugging purposes
(\code{"url_only"}).
}

When the results are returned as (a list of) tibbles, the column names
coming from the OpenCage API are prefixed with \code{"oc_"}.
}
\description{
Forward geocoding from a character vector of location names to latitude and
longitude tuples.
}
\examples{
\dontshow{if (oc_key_present() && oc_api_ok()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

# Geocode a single location, an address in this case
oc_forward(placename = "Triererstr 15, 99432, Weimar, Deutschland")

# Geocode multiple locations
locations <- c("Nantes", "Hamburg", "Los Angeles")
oc_forward(placename = locations)

# Use bounding box to help return accurate results
# for each placename
bounds <- oc_bbox(xmin = c(-2, 9, -119),
                  ymin = c(47, 53, 34),
                  xmax = c(0, 10, -117),
                  ymax = c(48, 54, 35))
oc_forward(placename = locations, bounds = bounds)

# Another way to help specify the desired results
# is with country codes.
oc_forward(placename = locations,
           countrycode = c("ca", "us", "co"))

# With multiple countrycodes per placename
oc_forward(placename = locations,
           countrycode = list(c("fr", "ca") , c("de", "us"), c("us", "co"))
           )

# Return results in a preferred language if possible
oc_forward(placename = c("Brugge", "Mechelen", "Antwerp"),
           language = "fr")

# Limit the number of results per placename and return json_list
oc_forward(placename = locations,
           bounds = bounds,
           limit = 1,
           return = "json_list")
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[=oc_forward_df]{oc_forward_df()}} for inputs as a data frame, or \code{\link[=oc_reverse]{oc_reverse()}} and
\code{\link[=oc_reverse_df]{oc_reverse_df()}} for reverse geocoding. For more information about the API
and the various parameters, see the \href{https://opencagedata.com/api}{OpenCage API documentation}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_reverse.R
\name{oc_reverse}
\alias{oc_reverse}
\title{Reverse geocoding}
\usage{
oc_reverse(
  latitude,
  longitude,
  return = c("df_list", "json_list", "geojson_list", "url_only"),
  language = NULL,
  min_confidence = NULL,
  no_annotations = TRUE,
  roadinfo = FALSE,
  no_dedupe = FALSE,
  abbrv = FALSE,
  add_request = FALSE,
  ...
)
}
\arguments{
\item{latitude, longitude}{Numeric vectors of latitude and longitude values.}

\item{return}{A character vector of length one indicating the return value of
the function, either a list of tibbles (\code{df_list}, the default), a JSON
list (\code{json_list}), a GeoJSON list (\code{geojson_list}), or the URL with which
the API would be called (\code{url_only}).}

\item{language}{An \href{https://en.wikipedia.org/wiki/IETF_language_tag}{IETF BCP 47 language tag} (such as "es" for
Spanish or "pt-BR" for Brazilian Portuguese). OpenCage will attempt to
return results in that language. Alternatively you can specify the "native"
tag, in which case OpenCage will attempt to return the response in the
"official" language(s). In case the \code{language} parameter is set to \code{NULL}
(which is the default), the tag is not recognized, or OpenCage does not
have a record in that language, the results will be returned in English.}

\item{min_confidence}{Numeric vector of integer values between 0 and 10
indicating the precision of the returned result as defined by its
geographical extent, (i.e. by the extent of the result's bounding box). See
the \href{https://opencagedata.com/api#confidence}{API documentation} for
details. Only results with at least the requested confidence will be
returned. Default is \code{NULL}.}

\item{no_annotations}{Logical vector indicating whether additional
information about the result location should be returned. \code{TRUE} by
default, which means that the results will not contain annotations.}

\item{roadinfo}{Logical vector indicating whether the geocoder should attempt
to match the nearest road (rather than an address) and provide additional
road and driving information. Default is \code{FALSE}.}

\item{no_dedupe}{Logical vector (default \code{FALSE}), when \code{TRUE} the results
will not be deduplicated.}

\item{abbrv}{Logical vector (default \code{FALSE}), when \code{TRUE} addresses in the
\code{formatted} field of the results are abbreviated (e.g. "Main St." instead
of "Main Street").}

\item{add_request}{Logical vector (default \code{FALSE}) indicating whether the
request is returned again with the results. If the \code{return} value is a
\code{df_list}, the query text is added as a column to the results. \code{json_list}
results will contain all request parameters, including the API key used!
This is currently ignored by OpenCage if return value is \code{geojson_list}.}

\item{...}{Ignored.}
}
\value{
Depending on the \code{return} argument, \code{oc_reverse} returns a list with
either
\itemize{
\item the results as tibbles (\code{"df_list"}, the default),
\item the results as JSON specified as a list (\code{"json_list"}),
\item the results as GeoJSON specified as a list (\code{"geojson_list"}),
or
\item the URL of the OpenCage API call for debugging purposes
(\code{"url_only"}).
}

When the results are returned as (a list of) tibbles, the column names
coming from the OpenCage API are prefixed with \code{"oc_"}.
}
\description{
Reverse geocoding from numeric vectors of latitude and longitude pairs to
the names and addresses of a location.
}
\examples{
\dontshow{if (oc_key_present() && oc_api_ok()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

# Reverse geocode a single location
oc_reverse(latitude = -36.85007, longitude = 174.7706)

# Reverse geocode multiple locations
lat <- c(47.21864, 53.55034, 34.05369)
lng <- c(-1.554136, 10.000654, -118.242767)

oc_reverse(latitude = lat, longitude = lng)

# Return results in a preferred language if possible
oc_reverse(latitude = lat, longitude = lng,
           language = "fr")

# Return results as a json list
oc_reverse(latitude = lat, longitude = lng,
           return = "json_list")
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[=oc_reverse_df]{oc_reverse_df()}} for inputs as a data frame, or \code{\link[=oc_forward]{oc_forward()}} and
\code{\link[=oc_forward]{oc_forward()}} for forward geocoding. For more information about the API
and the various parameters, see the \href{https://opencagedata.com/api}{OpenCage API documentation}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_key.R
\name{oc_key_present}
\alias{oc_key_present}
\title{Is an OpenCage API key present?}
\usage{
oc_key_present()
}
\value{
A single logical value, \code{TRUE} or \code{FALSE}.
}
\description{
Checks whether a potential OpenCage API key, i.e. a 32 character long,
alphanumeric string, is stored in the environment variable \code{OPENCAGE_KEY}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_forward.R
\name{oc_forward_df}
\alias{oc_forward_df}
\alias{oc_forward_df.data.frame}
\alias{oc_forward_df.character}
\title{Forward geocoding with data frames}
\usage{
oc_forward_df(...)

\method{oc_forward_df}{data.frame}(
  data,
  placename,
  bind_cols = TRUE,
  output = c("short", "all"),
  bounds = NULL,
  proximity = NULL,
  countrycode = NULL,
  language = NULL,
  limit = 1L,
  min_confidence = NULL,
  no_annotations = TRUE,
  roadinfo = FALSE,
  no_dedupe = FALSE,
  abbrv = FALSE,
  ...
)

\method{oc_forward_df}{character}(
  placename,
  output = c("short", "all"),
  bounds = NULL,
  proximity = NULL,
  countrycode = NULL,
  language = NULL,
  limit = 1L,
  min_confidence = NULL,
  no_annotations = TRUE,
  roadinfo = FALSE,
  no_dedupe = FALSE,
  abbrv = FALSE,
  ...
)
}
\arguments{
\item{...}{Ignored.}

\item{data}{A data frame.}

\item{placename}{An unquoted variable name of a character column or vector
with the location names or addresses to be geocoded.

If the locations are addresses, see \href{https://github.com/OpenCageData/opencagedata-misc-docs/blob/master/query-formatting.md}{OpenCage's instructions}
on how to format addresses for best forward geocoding results.}

\item{bind_cols}{When \code{bind_col = TRUE}, the default, the results are column
bound to \code{data}. When \code{FALSE}, the results are returned as a new tibble.}

\item{output}{A character vector of length one indicating whether only
latitude, longitude, and formatted address variables (\code{"short"}, the
default), or all variables (\code{"all"}) variables should be returned.}

\item{bounds}{A list of length one, or an unquoted variable name of a list
column of bounding boxes. Bounding boxes are named numeric vectors, each
with 4 coordinates forming the south-west and north-east corners of the
bounding box: \code{list(c(xmin, ymin, xmax, ymax))}. \code{bounds} restricts the
possible results to the supplied region. It can be specified with the
\code{\link[=oc_bbox]{oc_bbox()}} helper. For example: \code{bounds = oc_bbox(-0.563160, 51.280430, 0.278970, 51.683979)}. Default is \code{NULL}.}

\item{proximity}{A list of length one, or an unquoted variable name of a list
column of points. Points are named numeric vectors with latitude, longitude
coordinate pairs in decimal format. \code{proximity} provides OpenCage with a
hint to bias results in favour of those closer to the specified location.
It can be specified with the \code{\link[=oc_points]{oc_points()}} helper. For example: \code{proximity = oc_points(41.40139, 2.12870)}. Default is \code{NULL}.}

\item{countrycode}{Character vector, or an unquoted variable name of such a
vector, of two-letter codes as defined by the \href{https://www.iso.org/obp/ui/#search/code}{ISO 3166-1 Alpha 2} standard that restricts the
results to the given country or countries. E.g. "AR" for Argentina, "FR"
for France, "NZ" for the New Zealand. Multiple countrycodes per \code{placename}
must be wrapped in a list. Default is \code{NULL}.}

\item{language}{Character vector, or an unquoted variable name of such a
vector, of \href{https://en.wikipedia.org/wiki/IETF_language_tag}{IETF BCP 47 language tags} (such as "es" for
Spanish or "pt-BR" for Brazilian Portuguese). OpenCage will attempt to
return results in that language. Alternatively you can specify the "native"
tag, in which case OpenCage will attempt to return the response in the
"official" language(s). In case the \code{language} parameter is set to \code{NULL}
(which is the default), the tag is not recognized, or OpenCage does not
have a record in that language, the results will be returned in English.}

\item{limit}{Numeric vector of integer values, or an unquoted variable name
of such a vector, to determine the maximum number of results returned for
each \code{placename}. Integer values between 1 and 100 are allowed. Default is
1.}

\item{min_confidence}{Numeric vector of integer values, or an unquoted
variable name of such a vector, between 0 and 10 indicating the precision
of the returned result as defined by its geographical extent, (i.e. by the
extent of the result's bounding box). See the \href{https://opencagedata.com/api#confidence}{API documentation} for details. Only
results with at least the requested confidence will be returned. Default is
\code{NULL}).}

\item{no_annotations}{Logical vector, or an unquoted variable name of such a
vector, indicating whether additional information about the result location
should be returned. \code{TRUE} by default, which means that the results will
not contain annotations.}

\item{roadinfo}{Logical vector, or an unquoted variable name of such a
vector, indicating whether the geocoder should attempt to match the nearest
road (rather than an address) and provide additional road and driving
information. Default is \code{FALSE}.}

\item{no_dedupe}{Logical vector, or an unquoted variable name of such a
vector. Default is \code{FALSE}. When \code{TRUE} the results will not be
deduplicated.}

\item{abbrv}{Logical vector, or an unquoted variable name of such a vector.
Default is \code{FALSE}. When \code{TRUE} addresses in the \code{oc_formatted} variable of
the results are abbreviated (e.g. "Main St." instead of "Main Street").}
}
\value{
A tibble. Column names coming from the OpenCage API are prefixed with
\code{"oc_"}.
}
\description{
Forward geocoding from a column or vector of location names to latitude and
longitude tuples.
}
\examples{
\dontshow{if (oc_key_present() && oc_api_ok()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

library(tibble)
df <- tibble(id = 1:3,
             locations = c("Nantes", "Hamburg", "Los Angeles"))

# Return lat, lng, and formatted address
oc_forward_df(df, placename = locations)

# Return more detailed information about the locations
oc_forward_df(df, placename = locations, output = "all")

# Do not column bind results to input data frame
oc_forward_df(df, placename = locations, bind_cols = FALSE)

# Add more results by changing the limit from the default of 1.
oc_forward_df(df, placename = locations, limit = 5)

# Restrict results to a given bounding box
oc_forward_df(df, placename = locations,
              bounds = oc_bbox(-5, 45, 15, 55))

# oc_forward_df accepts unquoted column names for all
# arguments except bind_cols and output.
# This makes it possible to build up more detailed queries
# through the data frame passed to the data argument.

df2 <- add_column(df,
  bounds = oc_bbox(xmin = c(-2, 9, -119),
                   ymin = c(47, 53, 34),
                   xmax = c(0, 10, -117),
                   ymax = c(48, 54, 35)),
  limit = 1:3,
  countrycode = c("ca", "us", "co"),
  language = c("fr", "de", "en"))

# Use the bounds column to help return accurate results and
# language column to specify preferred language of results
oc_forward_df(df2, placename = locations,
              bounds = bounds,
              language = language)

# Different limit of results for each placename
oc_forward_df(df2, placename = locations,
              limit = limit)

# Specify the desired results by the countrycode column
oc_forward_df(df2, placename = locations,
              countrycode = countrycode)
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[=oc_forward]{oc_forward()}} for inputs as vectors, or \code{\link[=oc_reverse]{oc_reverse()}} and
\code{\link[=oc_reverse_df]{oc_reverse_df()}} for reverse geocoding. For more information about the API
and the various parameters, see the \href{https://opencagedata.com/api}{OpenCage API documentation}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{opencage-deprecated}
\alias{opencage-deprecated}
\title{Deprecated functions in opencage}
\description{
These functions still work but will be removed (defunct) in the next version.
}
\details{
\itemize{
\item \code{\link[=opencage_forward]{opencage_forward()}}
\item \code{\link[=opencage_reverse]{opencage_reverse()}}
\item \code{\link[=opencage_key]{opencage_key()}}
}
}
