
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

*Wow, no problems at all. :)**Wow, no problems at all. :)*