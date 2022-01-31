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
comtradr
========

[![Travis-CI Build Status](https://travis-ci.org/ropensci/comtradr.svg?branch=master)](https://travis-ci.org/ropensci/comtradr) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/comtradr?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/comtradr) [![codecov](https://codecov.io/github/ropensci/comtradr/branch/master/graphs/badge.svg)](https://codecov.io/github/ropensci/comtradr) [![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/comtradr)](https://cran.r-project.org/package=comtradr) [![](https://badges.ropensci.org/141_status.svg)](https://github.com/ropensci/onboarding/issues/141) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/comtradr)](https://www.rpackages.io/package/comtradr)

R package for interacting with the [UN Comtrade Database](https://comtrade.un.org/data/) public API. UN Comtrade provides historical data on the weights and value of specific goods shipped between countries, more info can be found [here](https://comtrade.un.org/). Full API documentation can be found [here](https://comtrade.un.org/data/doc/api/).

This package was inspired by the [R tutorial](https://comtrade.un.org/data/Doc/api/ex/r) posted by Comtrade, and is built using [httr](https://CRAN.R-project.org/package=httr) and [jsonlite](https://CRAN.R-project.org/package=jsonlite).

I've also built a Shiny app for visualizing comtrade shipping data, that's powered by this package. The app can be viewed [here](https://chrismuir.shinyapps.io/comtrade_plot_shinyapp/).

Please [report](https://github.com/ropensci/comtradr/issues) issues, comments, or feature requests.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

For information on citation of this package, use `citation("comtradr")`

Installation
------------

Install from CRAN:

``` r
install.packages("comtradr")
```

Or install from this repo:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/comtradr")
```

Example Usage
-------------

**Example 1**: Return all exports from China to South Korea, United States and Mexico, for all years

``` r
library(comtradr)

# Country names passed to the API query function must be spelled as they appear 
# in the Comtrade DB. Use "ct_country_lookup" to query the country DB and 
#return the exact spelling of specific countries.
ct_country_lookup("korea")
#> [1] "Dem. People's Rep. of Korea" "Rep. of Korea"

# Since we want South Korea, we'll use "Rep. of Korea" within the API query.
example1 <- ct_search(reporters = "China", 
                      partners = c("Rep. of Korea", "USA", "Mexico"), 
                      trade_direction = "exports")

# Inspect the return data
str(example1)
#> 'data.frame':    75 obs. of  35 variables:
#>  $ classification        : chr  "H4" "H4" "H4" "H4" ...
#>  $ year                  : int  2012 2012 2012 2013 2013 2013 2014 2014 2014 2015 ...
#>  $ period                : int  2012 2012 2012 2013 2013 2013 2014 2014 2014 2015 ...
#>  $ period_desc           : chr  "2012" "2012" "2012" "2013" ...
#>  $ aggregate_level       : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ is_leaf_code          : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ trade_flow_code       : int  2 2 2 2 2 2 2 2 2 2 ...
#>  $ trade_flow            : chr  "Export" "Export" "Export" "Export" ...
#>  $ reporter_code         : int  156 156 156 156 156 156 156 156 156 156 ...
#>  $ reporter              : chr  "China" "China" "China" "China" ...
#>  $ reporter_iso          : chr  "CHN" "CHN" "CHN" "CHN" ...
#>  $ partner_code          : int  410 484 842 410 484 842 410 484 842 410 ...
#>  $ partner               : chr  "Rep. of Korea" "Mexico" "USA" "Rep. of Korea" ...
#>  $ partner_iso           : chr  "KOR" "MEX" "USA" "KOR" ...
#>  $ second_partner_code   : logi  NA NA NA NA NA NA ...
#>  $ second_partner        : chr  NA NA NA NA ...
#>  $ second_partner_iso    : chr  NA NA NA NA ...
#>  $ customs_proc_code     : chr  NA NA NA NA ...
#>  $ customs               : chr  NA NA NA NA ...
#>  $ mode_of_transport_code: chr  NA NA NA NA ...
#>  $ mode_of_transport     : chr  NA NA NA NA ...
#>  $ commodity_code        : chr  "TOTAL" "TOTAL" "TOTAL" "TOTAL" ...
#>  $ commodity             : chr  "All Commodities" "All Commodities" "All Commodities" "All Commodities" ...
#>  $ qty_unit_code         : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ qty_unit              : chr  "No Quantity" "No Quantity" "No Quantity" "No Quantity" ...
#>  $ alt_qty_unit_code     : logi  NA NA NA NA NA NA ...
#>  $ alt_qty_unit          : chr  NA NA NA NA ...
#>  $ qty                   : logi  NA NA NA NA NA NA ...
#>  $ alt_qty               : logi  NA NA NA NA NA NA ...
#>  $ netweight_kg          : logi  NA NA NA NA NA NA ...
#>  $ gross_weight_kg       : logi  NA NA NA NA NA NA ...
#>  $ trade_value_usd       : num  8.77e+10 2.75e+10 3.52e+11 9.12e+10 2.90e+10 ...
#>  $ cif_trade_value_usd   : logi  NA NA NA NA NA NA ...
#>  $ fob_trade_value_usd   : logi  NA NA NA NA NA NA ...
#>  $ flag                  : int  0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "url")= chr "https://comtrade.un.org/api/get?max=50000&type=C&freq=A&px=HS&ps=all&r=156&p=410%2C842%2C484&rg=2&cc=TOTAL&fmt=json&head=H"
#>  - attr(*, "time_stamp")= POSIXct, format: "2018-05-04 20:03:23"
#>  - attr(*, "req_duration")= num 1.02
```

**Example 2**: Return all exports related to shrimp from Thailand to all other countries, for years 2007 thru 2011

``` r
library(comtradr)

# Fetch all shrimp related commodity codes from the Comtrade commodities DB. 
# This vector of codes will get passed to the API query.
shrimp_codes <- ct_commodity_lookup("shrimp", return_code = TRUE, return_char = TRUE)

# API query.
example2 <- ct_search(reporters = "Thailand", 
                      partners = "All", 
                      trade_direction = "exports", 
                      start_date = 2007, 
                      end_date = 2011, 
                      commod_codes = shrimp_codes)

# Inspect the output
str(example2)
#> 'data.frame':    1203 obs. of  35 variables:
#>  $ classification        : chr  "H3" "H3" "H3" "H3" ...
#>  $ year                  : int  2007 2007 2007 2007 2007 2007 2007 2007 2007 2007 ...
#>  $ period                : int  2007 2007 2007 2007 2007 2007 2007 2007 2007 2007 ...
#>  $ period_desc           : chr  "2007" "2007" "2007" "2007" ...
#>  $ aggregate_level       : int  6 6 6 6 6 6 6 6 6 6 ...
#>  $ is_leaf_code          : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ trade_flow_code       : int  2 2 2 2 2 2 2 2 2 2 ...
#>  $ trade_flow            : chr  "Export" "Export" "Export" "Export" ...
#>  $ reporter_code         : int  764 764 764 764 764 764 764 764 764 764 ...
#>  $ reporter              : chr  "Thailand" "Thailand" "Thailand" "Thailand" ...
#>  $ reporter_iso          : chr  "THA" "THA" "THA" "THA" ...
#>  $ partner_code          : int  0 36 40 48 56 104 116 124 152 156 ...
#>  $ partner               : chr  "World" "Australia" "Austria" "Bahrain" ...
#>  $ partner_iso           : chr  "WLD" "AUS" "AUT" "BHR" ...
#>  $ second_partner_code   : logi  NA NA NA NA NA NA ...
#>  $ second_partner        : chr  NA NA NA NA ...
#>  $ second_partner_iso    : chr  NA NA NA NA ...
#>  $ customs_proc_code     : chr  NA NA NA NA ...
#>  $ customs               : chr  NA NA NA NA ...
#>  $ mode_of_transport_code: chr  NA NA NA NA ...
#>  $ mode_of_transport     : chr  NA NA NA NA ...
#>  $ commodity_code        : chr  "030613" "030613" "030613" "030613" ...
#>  $ commodity             : chr  "Shrimps & prawns, whether/not in shell, frozen" "Shrimps & prawns, whether/not in shell, frozen" "Shrimps & prawns, whether/not in shell, frozen" "Shrimps & prawns, whether/not in shell, frozen" ...
#>  $ qty_unit_code         : int  8 8 8 8 8 8 8 8 8 8 ...
#>  $ qty_unit              : chr  "Weight in kilograms" "Weight in kilograms" "Weight in kilograms" "Weight in kilograms" ...
#>  $ alt_qty_unit_code     : logi  NA NA NA NA NA NA ...
#>  $ alt_qty_unit          : chr  NA NA NA NA ...
#>  $ qty                   : int  169654441 5545602 1265 29780 2721318 750 8510 13088545 4930 3410678 ...
#>  $ alt_qty               : logi  NA NA NA NA NA NA ...
#>  $ netweight_kg          : int  169654441 5545602 1265 29780 2721318 750 8510 13088545 4930 3410678 ...
#>  $ gross_weight_kg       : logi  NA NA NA NA NA NA ...
#>  $ trade_value_usd       : int  1084677273 36120291 11888 124668 16061545 4521 74842 77292118 64218 18400152 ...
#>  $ cif_trade_value_usd   : logi  NA NA NA NA NA NA ...
#>  $ fob_trade_value_usd   : logi  NA NA NA NA NA NA ...
#>  $ flag                  : int  0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "url")= chr "https://comtrade.un.org/api/get?max=50000&type=C&freq=A&px=HS&ps=2007%2C2008%2C2009%2C2010%2C2011&r=764&p=all&r"| __truncated__
#>  - attr(*, "time_stamp")= POSIXct, format: "2018-05-04 20:03:26"
#>  - attr(*, "req_duration")= num 3.19
```

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
comtradr 0.2.2.09000
====================

## NEW FEATURES

* Modifications to `ct_search()` to add support for commodity code `ag6` ([#30](https://github.com/ropensci/comtradr/pull/30))

## BUG FIXES

* Function `ct_register_token()` now checks if the provided token is recognized by the official API and only grants "premium" credentials if it is ([#34](https://github.com/ropensci/comtradr/issues/34)).

* Passing an API token string to `ct_register_token()` now properly bumps the hourly rate limit up to 10,000
([#21](https://github.com/ropensci/comtradr/issues/21)).

* In func `ct_search()`, passing a character vector of long-form commodity descriptions to arg `commod_codes` will now 
throw an error prior to making an API call, which would fail ([#24](https://github.com/ropensci/comtradr/issues/24)).

* Update the country package data, to stay up to date with the reporter/partner country table that Comtrade is using. This is an update to
the file `inst/extdata/country_table.rda`. ([#29](https://github.com/ropensci/comtradr/issues/29)).

* In func `ct_search()`, improve error messaging when an input country is invalid. ([#31](https://github.com/ropensci/comtradr/issues/31)).

* In func `ct_search()`, fix bug in which running queries using the `SITCrev2` commodity type was returning raw HTML (as opposed to json data). ([#27](https://github.com/ropensci/comtradr/issues/27)).


comtradr 0.2.2
====================

* Remove unused dependency `methods` from `Imports`.

comtradr 0.2.1
====================

## NEW FEATURES

* Modifications to `ct_search()` to allow for pulling all monthly data for an entire year in a single query ([#14](https://github.com/ropensci/comtradr/issues/14))
* For function `ct_search()`, expanded the valid input types for args `start_date` and `end_date` ([#10](https://github.com/ropensci/comtradr/issues/10)).

## BUG FIXES

* `ct_search()` now supports all commodity classifications offered by UN Comtrade ([#16](https://github.com/ropensci/comtradr/issues/16)).
* The updates generated by function `ct_update_databases()` are now properly preserved between R sessions ([#11](https://github.com/ropensci/comtradr/issues/11)).
* Passing `"services"` to arg `type` within function `ct_search()` now uses commodity classification `EB02` by default (previously this would throw an error, fixes [#6](https://github.com/ropensci/comtradr/issues/6)).
* When using commodity classification `EB02` within function `ct_search()`, passing `"TOTAL"` to arg `commod_codes` no longer returns zero results ([#7](https://github.com/ropensci/comtradr/issues/7)).
* `ct_commodity_lookup()` no longer returns zero results when passing all caps input to arg `search_terms` ([#9](https://github.com/ropensci/comtradr/issues/9)).


comtradr 0.1.0
===================

## PKG API CHANGES

* Eliminated functions `ct_commodities_table` and `ct_countries_table`.
* Added new functions `ct_update_databases`, `ct_use_pretty_cols`, `ct_commodity_db_type`, `ct_register_token`, `ct_get_reset_time`, `ct_get_remaining_hourly_queries`.
* Renamed functions: `commodity_lookup` is now `ct_commodity_lookup`, `country_lookup` is now `ct_country_lookup`.
* The commodity and country reference tables are now saved as cached package data, and accessed by `comtradr` functions when necessary. This replaces the need for functions `ct_commodities_table` and `ct_countries_table`.
* Reorder function arguments within function `ct_search`.
* Changed some function argument names to ensure `snake_case` is being used throughout the package.
* `ct_search` now returns a data frame, as opposed to a list.

## MINOR CHANGES

* Added a vignette directory, with an "Intro to comtradr" vignette.
* API requests are now throttled based on the [rate limits](https://comtrade.un.org/data/doc/api/#Limits) imposed by the UN Comtrade.
* Added function for setting a valid API key/token (`ct_register_token`).
* Appending API metadata to each returned data frame as attributes (url of the API call, date-time of the query, duration of the query in seconds).
* Added package level man page.
* Now using native R errors/warnings, as opposed to nesting API status codes in a returned list.
* `Imports` changes: remove `dplyr`, add `magrittr` and `purrr`.
* Expand and improve test coverage via [testthat](https://github.com/hadley/testthat).

## BUG FIXES

* The issues related to type-safety in function `commodity_lookup` have been fixed by importing `purrr` and using `purrr::map` in place of `sapply`. This fixes [issue #2](https://github.com/ropensci/comtradr/issues/2) and [issue #3](https://github.com/ropensci/comtradr/issues/3).


comtradr 0.0.2 (2017-07-03)
===========================

## NEW FEATURES

* commodity_lookup(): Expanded function to accept multiple commodities or commodity codes (as either character vector or numeric vector). Also added argument "return_char" that allows the user to specify list output or char vector output, and argument "return_code" that specifies output of commodity descriptions or commodity codes.

## MINOR IMPROVEMENTS

* Add unit tests via [testthat](https://github.com/hadley/testthat).


comtradr 0.0.1 (2017-04-06)
===========================

## NEW FEATURES

* released to CRAN

## Test environments
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.5.1
* macOS High Sierra 10.13.3 (on travis-ci), R 3.5.0
* CRAN win-builder, R Under development (unstable) (2018-10-04 r75399)

## R CMD check results
0 errors | 0 warnings | 1 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Chris Muir <chrismuirRVA@gmail.com>'

## Downstream dependencies
* Reverse imports: ITNr

revdepcheck results:
0 errors | 0 warnings | 0 notes

---

Package changes implemented in this version:

* Remove unused dependency `methods` from `Imports`.
# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.5.0 (2018-04-23) |
|system   |x86_64, darwin15.6.0         |
|ui       |RStudio (1.1.383)            |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|tz       |America/New_York             |
|date     |2018-10-04                   |

## Packages

|package   |*  |version |date       |source        |
|:---------|:--|:-------|:----------|:-------------|
|comtradr  |*  |0.2.1   |2018-05-05 |cran (@0.2.1) |
|dplyr     |   |0.7.6   |2018-06-29 |cran (@0.7.6) |
|ggplot2   |   |3.0.0   |2018-07-03 |cran (@3.0.0) |
|httr      |   |1.3.1   |2017-08-20 |cran (@1.3.1) |
|jsonlite  |   |1.5     |2017-06-01 |cran (@1.5)   |
|knitr     |   |1.20    |2018-02-20 |cran (@1.20)  |
|magrittr  |   |1.5     |2014-11-22 |cran (@1.5)   |
|purrr     |   |0.2.5   |2018-05-29 |cran (@0.2.5) |
|rmarkdown |   |1.10    |2018-06-11 |cran (@1.10)  |
|testthat  |   |2.0.0   |2017-12-13 |cran (@2.0.0) |

# Check results

1 packages

|package |version | errors| warnings| notes|
|:-------|:-------|------:|--------:|-----:|
|ITNr    |0.2.0   |      0|        0|     0|

0 errors | 0 warnings | 0 notes

---
output:
  github_document: default
  html_document: default
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

comtradr
=======

[![Travis-CI Build Status](https://travis-ci.org/ropensci/comtradr.svg?branch=master)](https://travis-ci.org/ropensci/comtradr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/comtradr?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/comtradr)
[![codecov](https://codecov.io/github/ropensci/comtradr/branch/master/graphs/badge.svg)](https://codecov.io/github/ropensci/comtradr)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/comtradr)](https://cran.r-project.org/package=comtradr)
[![](https://badges.ropensci.org/141_status.svg)](https://github.com/ropensci/onboarding/issues/141)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/comtradr)](https://www.rpackages.io/package/comtradr)

R package for interacting with the [UN Comtrade Database](https://comtrade.un.org/data/) public API. UN Comtrade provides historical data on the weights and value of 
specific goods shipped between countries, more info can be found [here](https://comtrade.un.org/). Full API documentation can be found 
[here](https://comtrade.un.org/data/doc/api/).

This package was inspired by the [R tutorial](https://comtrade.un.org/data/Doc/api/ex/r) posted by Comtrade, and is built using
[httr](https://CRAN.R-project.org/package=httr) and [jsonlite](https://CRAN.R-project.org/package=jsonlite).

I've also built a Shiny app for visualizing comtrade shipping data, that's powered by this package. The app can be viewed [here](https://chrismuir.shinyapps.io/comtrade_plot_shinyapp/).

Please [report](https://github.com/ropensci/comtradr/issues) issues, comments, or feature requests.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

For information on citation of this package, use `citation("comtradr")`

## Installation

Install from CRAN:
```{r eval=FALSE}
install.packages("comtradr")
```

Or install from this repo:
```{r eval=FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/comtradr")
```

## Example Usage

**Example 1**: Return all exports from China to South Korea, United States and Mexico, for all years

```{r}
library(comtradr)

# Country names passed to the API query function must be spelled as they appear 
# in the Comtrade DB. Use "ct_country_lookup" to query the country DB and 
#return the exact spelling of specific countries.
ct_country_lookup("korea")

# Since we want South Korea, we'll use "Rep. of Korea" within the API query.
example1 <- ct_search(reporters = "China", 
                      partners = c("Rep. of Korea", "USA", "Mexico"), 
                      trade_direction = "exports")

# Inspect the return data
str(example1)
```

**Example 2**: Return all exports related to shrimp from Thailand to all other countries, for years 2007 thru 2011

```{r}
library(comtradr)

# Fetch all shrimp related commodity codes from the Comtrade commodities DB. 
# This vector of codes will get passed to the API query.
shrimp_codes <- ct_commodity_lookup("shrimp", return_code = TRUE, return_char = TRUE)

# API query.
example2 <- ct_search(reporters = "Thailand", 
                      partners = "All", 
                      trade_direction = "exports", 
                      start_date = 2007, 
                      end_date = 2011, 
                      commod_codes = shrimp_codes)

# Inspect the output
str(example2)
```

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Intro to comtradr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Intro to comtradr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE}
knitr::opts_chunk$set(comment = "#>", collapse = TRUE, fig.width = 9, fig.height = 6)
```

## Package information

API wrapper for the [UN Comtrade Database](https://comtrade.un.org/data/), which features inter-country trade data dating back to the early 1990's. Full API documentation can be found [here](https://comtrade.un.org/data/doc/api/). This package allows users to interact with the API directly from R, and features functions for making queries and importing data.

## Install and load comtradr

Install from CRAN:
```{r eval = FALSE}
install.packages("comtradr")
```
Or install the development version from GitHub:
```{r eval = FALSE}
devtools::install_github("ChrisMuir/comtradr")
```
Load comtradr
```{r}
library(comtradr)
```

## Making API calls
Lets say we want to get data on all imports into the United States from Germany, France, Japan, and Mexico, for all years.
```{r, echo = FALSE}
v_data_1 <- system.file("extdata", "vignette_data_1.rda", package = "comtradr")
if (!file.exists(v_data_1)) {
  stop("internal vignette data set '~/extdata/vignette_data_1.rda' not found", call. = FALSE)
}
load(v_data_1)
```

```{r, eval = FALSE}
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports")
```
API calls return a tidy data frame.
```{r}
str(q)
```

Here are a few more examples to show the different parameter options:

Limit the search range to shipments between 2010 and 2014.
```{r, eval = FALSE}
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports", 
               start_date = 2010, 
               end_date = 2014)
```

By default, the return data is in yearly amounts. We can pass `"monthly"` to arg `freq` to return data in monthly amounts, however the API limits each "monthly" query to a single year.
```{r, eval = FALSE}
# Get all monthly data for a single year (API max of 12 months per call).
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports", 
               start_date = 2012, 
               end_date = 2012, 
               freq = "monthly")

# Get monthly data for a specific span of months (API max of five months per call).
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports", 
               start_date = "2012-03", 
               end_date = "2012-07", 
               freq = "monthly")
```

Countries passed to parameters `reporters` and `partners` must be spelled as they appear in the Comtrade country reference table. Function `ct_country_lookup` allows us to query the country reference table.
```{r}
ct_country_lookup("korea", "reporter")
ct_country_lookup("bolivia", "partner")
```
```{r, eval = FALSE}
q <- ct_search(reporters = "Rep. of Korea", 
               partners = "Bolivia (Plurinational State of)", 
               trade_direction = "all")
```

Search trade related to specific commodities (say, tomatoes). We can query the Comtrade commodity reference table to see all of the different commodity descriptions available for tomatoes.
```{r}
ct_commodity_lookup("tomato")
```
If we want to search for shipment data on all of the commodity descriptions listed, then we can simply adjust the parameters for `ct_commodity_lookup` so that it will return only the codes, which can then be passed along to `ct_search`.
```{r, eval = FALSE}
tomato_codes <- ct_commodity_lookup("tomato", 
                                    return_code = TRUE, 
                                    return_char = TRUE)

q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Mexico"), 
               trade_direction = "all", 
               commod_codes = tomato_codes)
```
On the other hand, if we wanted to exclude juices and sauces from our search, we can pass a vector of the relevant codes to the API call.
```{r, eval = FALSE}
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Mexico"), 
               trade_direction = "all", 
               commod_codes = c("0702", "070200", "2002", "200210", "200290"))
```

## API search metadata

In addition to the trade data, each API return object contains metadata as attributes.
```{r}
# The url of the API call.
attributes(q)$url
# The date-time of the API call.
attributes(q)$time_stamp

# The total duration of the API call, in seconds.
attributes(q)$req_duration
```

## More on the lookup functions

Functions `ct_country_lookup` and `ct_commodity_lookup` are both able to take multiple search terms as input.
```{r}
ct_country_lookup(c("Belgium", "vietnam", "brazil"), "reporter")

ct_commodity_lookup(c("tomato", "trout"), return_char = TRUE)
```


`ct_commodity_lookup` can return a vector (as seen above) or a named list, using parameter `return_char`
```{r}
ct_commodity_lookup(c("tomato", "trout"), return_char = FALSE)
```


For `ct_commodity_lookup`, if any of the input search terms return zero results and parameter `verbose` is set to `TRUE`, a warning will be printed to console (set `verbose` to `FALSE` to turn off this feature).
```{r}
ct_commodity_lookup(c("tomato", "sldfkjkfdsklsd"), verbose = TRUE)
```

## API rate limits

The Comtrade API imposes rate limits on both guest users and premium users. `comtradr` features automated throttling of API calls to ensure the user stays within the limits defined by Comtrade. Below is a breakdown of those limits, API docs on these details can be found [here](https://comtrade.un.org/data/doc/api/#Authentication).

* Without user token: 1 request per second, 100 requests per hour.
* With valid user token: 1 request per second, 10,000 requests per hour.

In addition to these rate limits, the API imposes some limits on parameter combinations.

* Between args `reporters`, `partners`, and the query date range, only one of these three may use the catch-all input "All".
* For the same group of three (`reporters`, `partners`, date range), if the input is not "All", then the maximum number of input values for each is five. For date range, if not using "All", then the `start_date` and `end_date` must not span more than five months or five years. There is one exception to this rule, if arg `freq` is "monthly", then a single year can be passed to `start_date` and `end_date` and the API will return all of the monthly data for that year.
* For arg `commod_codes`, if not using input "All", then the maximum number of input values is 20 (although "All" is always a valid input).

Additionally, the maximum number of returned records from a single query without a token is 50,000. With a token, that number is 250,000.

`comtradr` features a few functions for working with the API rate limits and tokens.

* `ct_register_token()` allows the user to set an API token within the package environment. 
* `ct_get_remaining_hourly_queries()` will return the number of remaining queries for the current hour.
* `ct_get_reset_time()` will return the date/time in which the current hourly time limit will reset, as a `POSIXct` object.

## Package Data

`comtradr` ships with a few different package data objects, and functions for interacting with and using the package data.

**Country/Commodity Reference Tables**

As explained previously, making API calls with `comtradr` often requires the user to query the country reference table and/or the commodity reference table (this is done using functions `ct_country_lookup` and `ct_commodity_lookup`). Both of these reference tables are generated by the UN Comtrade, and are updated roughly once a year. Since they're updated infrequently, both tables are saved as cached data objects within the `comtradr` package, and are referenced by the package functions when needed.

`comtradr` features a function, `ct_update_databases`, for checking the Comtrade website for updates to either reference table. If updates are found, the function will download the updated table, save it to the package directory, and make it available during the current R session. It will also print a message indicating whether updates were found, like so:
```{r, eval = FALSE}
ct_update_databases()
#> All DB's are up to date, no action required
```
If any updates are found, the message will state which reference table(s) were updated.

The user may force download of both reference tables (regardless of whether updates exist) by using arg `force = TRUE` within function `ct_update_databases`. This is useful in the event that either reference table file is deleted or removed from the package directory. If this is the case, and either reference table file is not found upon package load, then any subsequent `comtradr` functions that require the use of a reference table will result in an error, with the error message prompting the user to run `ct_update_databases(force = TRUE)`.

Additionally, the Comtrade API features a number of different commodity reference tables, based on different trade data classification schemes (for more details, see [this](https://comtrade.un.org/data/doc/api/#DataAvailabilityRequests) page from the API docs). `comtradr` ships with the commodity table for the "Harmonized System", or "HS", scheme. The user may download any of the available commodity tables by specifying arg `commodity_type` within function `ct_update_databases` (e.g., `ct_update_databases(commodity_type = "SITC")` will download the commodity table that follows the "Standard International Trade Classification" scheme). Doing this will replace the commodity table on file with the one specified. To see the classification scheme of the commodity table currently on file, use `ct_commodity_db_type`.
```{r}
ct_commodity_db_type()
```

**"Polished" Column Headers**

`ct_pretty_cols` is a named vector of column header values that provide the option of using column headers that are more polished and human-friendly than those returned by the API function `ct_search`. The polished column headers may be useful when plotting the Comtrade data, or for use in publication tables. The data can be accessed directly by using `data("ct_pretty_cols")`, but there is also a function for applying the polished headers to `comtradr` data frames, `ct_use_pretty_cols`. Below is a quick demonstration.

```{r}
# Column headers returned from function ct_search
colnames(q)
```

```{r}
# Apply polished column headers
q <- ct_use_pretty_cols(q)

# Print new column headers.
colnames(q)
```

## Visualize

Once the data is collected, we can use it to create some basic visualizations.

**Plot 1**: Plot total value (USD) of Chinese exports to Mexico, South Korea and the United States, by year.

```{r, echo = FALSE}
v_data_2 <- system.file("extdata", "vignette_data_2.rda", package = "comtradr")
if (!file.exists(v_data_2)) {
  stop("internal vignette data set '~/extdata/vignette_data_2.rda' not found", call. = FALSE)
}
load(v_data_2)
```
```{r, eval = FALSE}
# Comtrade api query.
df <- ct_search(reporters = "China", 
                partners = c("Rep. of Korea", "USA", "Mexico"), 
                trade_direction = "exports")
```

```{r, warning = FALSE, message = FALSE}
library(ggplot2)

# Apply polished col headers.
df <- ct_use_pretty_cols(df)

# Create plot.
ggplot(df, aes(Year, `Trade Value usd`, color = factor(`Partner Country`), 
               shape = factor(`Partner Country`))) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(min(df$Year), max(df$Year)), 
                     breaks = seq.int(min(df$Year), max(df$Year), 2)) +
  scale_color_manual(values = c("orange", "blue", "red"), 
                     name = "Destination\nCountry") +
  scale_shape_discrete(name = "Destination\nCountry") +
  labs(title = "Total Value (USD) of Chinese Exports, by Year") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

**Plot 2**: Plot the top eight destination countries/areas of Thai shrimp exports, by weight (KG), for 2007 - 2011.

```{r, echo = FALSE}
v_data_3 <- system.file("extdata", "vignette_data_3.rda", package = "comtradr")
if (!file.exists(v_data_3)) {
  stop("internal vignette data set '~/extdata/vignette_data_3.rda' not found", call. = FALSE)
}
load(v_data_3)
```
```{r, eval = FALSE}
# First, collect commodity codes related to shrimp.
shrimp_codes <- ct_commodity_lookup("shrimp", 
                                    return_code = TRUE, 
                                    return_char = TRUE)

# Comtrade api query.
df <- ct_search(reporters = "Thailand", 
                partners = "All", 
                trade_direction = "exports", 
                start_date = 2007, 
                end_date = 2011, 
                commod_codes = shrimp_codes)
```

```{r, warning = FALSE, message = FALSE}
library(ggplot2)
library(dplyr)

# Apply polished col headers.
df <- ct_use_pretty_cols(df)

# Create country specific "total weight per year" dataframe for plotting.
plotdf <- df %>% 
  group_by_(.dots = c("`Partner Country`", "Year")) %>% 
  summarise(kg = as.numeric(sum(`Net Weight kg`, na.rm = TRUE))) %>% 
  as_data_frame()

# Get vector of the top 8 destination countries/areas by total weight shipped 
# across all years, then subset plotdf to only include observations related 
# to those countries/areas.
top8 <- plotdf %>% 
  group_by(`Partner Country`) %>% 
  summarise(kg = as.numeric(sum(kg, na.rm = TRUE))) %>% 
  top_n(8, kg) %>%
  arrange(desc(kg)) %>% 
  .[["Partner Country"]]
plotdf <- plotdf %>% filter(`Partner Country` %in% top8)

# Create plots (y-axis is NOT fixed across panels, this will allow us to ID 
# trends over time within each country/area individually).
qplot(Year, kg, data = plotdf) + 
  geom_line(data = plotdf[plotdf$`Partner Country` %in% names(which(table(plotdf$`Partner Country`) > 1)), ]) + 
  xlim(min(plotdf$Year), max(plotdf$Year)) + 
  labs(title = "Weight (KG) of Thai Shrimp Exports, by Destination Area, 2007 - 2011") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.text = element_text(size = 7)) + 
  facet_wrap(~factor(`Partner Country`, levels = top8), scales = "free", nrow = 2, ncol = 4)
```
---
title: "Intro to comtradr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Intro to comtradr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE}
knitr::opts_chunk$set(comment = "#>", collapse = TRUE, fig.width = 9, fig.height = 6)
```

## Package information

API wrapper for the [UN Comtrade Database](https://comtrade.un.org/data/), which features inter-country trade data dating back to the early 1990's. Full API documentation can be found [here](https://comtrade.un.org/data/doc/api/). This package allows users to interact with the API directly from R, and features functions for making queries and importing data.

## Install and load comtradr

Install from CRAN:
```{r eval = FALSE}
install.packages("comtradr")
```
Or install the development version from GitHub:
```{r eval = FALSE}
devtools::install_github("ChrisMuir/comtradr")
```
Load comtradr
```{r}
library(comtradr)
```

## Making API calls
Lets say we want to get data on all imports into the United States from Germany, France, Japan, and Mexico, for all years.
```{r, echo = FALSE}
v_data_1 <- system.file("extdata", "vignette_data_1.rda", package = "comtradr")
if (!file.exists(v_data_1)) {
  stop("internal vignette data set '~/extdata/vignette_data_1.rda' not found", call. = FALSE)
}
load(v_data_1)
```

```{r, eval = FALSE}
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports")
```
API calls return a tidy data frame.
```{r}
str(q)
```

Here are a few more examples to show the different parameter options:

Limit the search range to shipments between 2010 and 2014.
```{r, eval = FALSE}
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports", 
               start_date = 2010, 
               end_date = 2014)
```

By default, the return data is in yearly amounts. We can pass `"monthly"` to arg `freq` to return data in monthly amounts, however the API limits each "monthly" query to a single year.
```{r, eval = FALSE}
# Get all monthly data for a single year (API max of 12 months per call).
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports", 
               start_date = 2012, 
               end_date = 2012, 
               freq = "monthly")

# Get monthly data for a specific span of months (API max of five months per call).
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Japan", "Mexico"), 
               trade_direction = "imports", 
               start_date = "2012-03", 
               end_date = "2012-07", 
               freq = "monthly")
```

Countries passed to parameters `reporters` and `partners` must be spelled as they appear in the Comtrade country reference table. Function `ct_country_lookup` allows us to query the country reference table.
```{r}
ct_country_lookup("korea", "reporter")
ct_country_lookup("bolivia", "partner")
```
```{r, eval = FALSE}
q <- ct_search(reporters = "Rep. of Korea", 
               partners = "Bolivia (Plurinational State of)", 
               trade_direction = "all")
```

Search trade related to specific commodities (say, tomatoes). We can query the Comtrade commodity reference table to see all of the different commodity descriptions available for tomatoes.
```{r}
ct_commodity_lookup("tomato")
```
If we want to search for shipment data on all of the commodity descriptions listed, then we can simply adjust the parameters for `ct_commodity_lookup` so that it will return only the codes, which can then be passed along to `ct_search`.
```{r, eval = FALSE}
tomato_codes <- ct_commodity_lookup("tomato", 
                                    return_code = TRUE, 
                                    return_char = TRUE)

q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Mexico"), 
               trade_direction = "all", 
               commod_codes = tomato_codes)
```
On the other hand, if we wanted to exclude juices and sauces from our search, we can pass a vector of the relevant codes to the API call.
```{r, eval = FALSE}
q <- ct_search(reporters = "USA", 
               partners = c("Germany", "France", "Mexico"), 
               trade_direction = "all", 
               commod_codes = c("0702", "070200", "2002", "200210", "200290"))
```

## API search metadata

In addition to the trade data, each API return object contains metadata as attributes.
```{r}
# The url of the API call.
attributes(q)$url
# The date-time of the API call.
attributes(q)$time_stamp

# The total duration of the API call, in seconds.
attributes(q)$req_duration
```

## More on the lookup functions

Functions `ct_country_lookup` and `ct_commodity_lookup` are both able to take multiple search terms as input.
```{r}
ct_country_lookup(c("Belgium", "vietnam", "brazil"), "reporter")

ct_commodity_lookup(c("tomato", "trout"), return_char = TRUE)
```


`ct_commodity_lookup` can return a vector (as seen above) or a named list, using parameter `return_char`
```{r}
ct_commodity_lookup(c("tomato", "trout"), return_char = FALSE)
```


For `ct_commodity_lookup`, if any of the input search terms return zero results and parameter `verbose` is set to `TRUE`, a warning will be printed to console (set `verbose` to `FALSE` to turn off this feature).
```{r}
ct_commodity_lookup(c("tomato", "sldfkjkfdsklsd"), verbose = TRUE)
```

## API rate limits

The Comtrade API imposes rate limits on both guest users and premium users. `comtradr` features automated throttling of API calls to ensure the user stays within the limits defined by Comtrade. Below is a breakdown of those limits, API docs on these details can be found [here](https://comtrade.un.org/data/doc/api/#Authentication).

* Without user token: 1 request per second, 100 requests per hour.
* With valid user token: 1 request per second, 10,000 requests per hour.

In addition to these rate limits, the API imposes some limits on parameter combinations.

* Between args `reporters`, `partners`, and the query date range, only one of these three may use the catch-all input "All".
* For the same group of three (`reporters`, `partners`, date range), if the input is not "All", then the maximum number of input values for each is five. For date range, if not using "All", then the `start_date` and `end_date` must not span more than five months or five years. There is one exception to this rule, if arg `freq` is "monthly", then a single year can be passed to `start_date` and `end_date` and the API will return all of the monthly data for that year.
* For arg `commod_codes`, if not using input "All", then the maximum number of input values is 20 (although "All" is always a valid input).

Additionally, the maximum number of returned records from a single query without a token is 50,000. With a token, that number is 250,000.

`comtradr` features a few functions for working with the API rate limits and tokens.

* `ct_register_token()` allows the user to set an API token within the package environment. 
* `ct_get_remaining_hourly_queries()` will return the number of remaining queries for the current hour.
* `ct_get_reset_time()` will return the date/time in which the current hourly time limit will reset, as a `POSIXct` object.

## Package Data

`comtradr` ships with a few different package data objects, and functions for interacting with and using the package data.

**Country/Commodity Reference Tables**

As explained previously, making API calls with `comtradr` often requires the user to query the country reference table and/or the commodity reference table (this is done using functions `ct_country_lookup` and `ct_commodity_lookup`). Both of these reference tables are generated by the UN Comtrade, and are updated roughly once a year. Since they're updated infrequently, both tables are saved as cached data objects within the `comtradr` package, and are referenced by the package functions when needed.

`comtradr` features a function, `ct_update_databases`, for checking the Comtrade website for updates to either reference table. If updates are found, the function will download the updated table, save it to the package directory, and make it available during the current R session. It will also print a message indicating whether updates were found, like so:
```{r, eval = FALSE}
ct_update_databases()
#> All DB's are up to date, no action required
```
If any updates are found, the message will state which reference table(s) were updated.

The user may force download of both reference tables (regardless of whether updates exist) by using arg `force = TRUE` within function `ct_update_databases`. This is useful in the event that either reference table file is deleted or removed from the package directory. If this is the case, and either reference table file is not found upon package load, then any subsequent `comtradr` functions that require the use of a reference table will result in an error, with the error message prompting the user to run `ct_update_databases(force = TRUE)`.

Additionally, the Comtrade API features a number of different commodity reference tables, based on different trade data classification schemes (for more details, see [this](https://comtrade.un.org/data/doc/api/#DataAvailabilityRequests) page from the API docs). `comtradr` ships with the commodity table for the "Harmonized System", or "HS", scheme. The user may download any of the available commodity tables by specifying arg `commodity_type` within function `ct_update_databases` (e.g., `ct_update_databases(commodity_type = "SITC")` will download the commodity table that follows the "Standard International Trade Classification" scheme). Doing this will replace the commodity table on file with the one specified. To see the classification scheme of the commodity table currently on file, use `ct_commodity_db_type`.
```{r}
ct_commodity_db_type()
```

**"Polished" Column Headers**

`ct_pretty_cols` is a named vector of column header values that provide the option of using column headers that are more polished and human-friendly than those returned by the API function `ct_search`. The polished column headers may be useful when plotting the Comtrade data, or for use in publication tables. The data can be accessed directly by using `data("ct_pretty_cols")`, but there is also a function for applying the polished headers to `comtradr` data frames, `ct_use_pretty_cols`. Below is a quick demonstration.

```{r}
# Column headers returned from function ct_search
colnames(q)
```

```{r}
# Apply polished column headers
q <- ct_use_pretty_cols(q)

# Print new column headers.
colnames(q)
```

## Visualize

Once the data is collected, we can use it to create some basic visualizations.

**Plot 1**: Plot total value (USD) of Chinese exports to Mexico, South Korea and the United States, by year.

```{r, echo = FALSE}
v_data_2 <- system.file("extdata", "vignette_data_2.rda", package = "comtradr")
if (!file.exists(v_data_2)) {
  stop("internal vignette data set '~/extdata/vignette_data_2.rda' not found", call. = FALSE)
}
load(v_data_2)
```
```{r, eval = FALSE}
# Comtrade api query.
df <- ct_search(reporters = "China", 
                partners = c("Rep. of Korea", "USA", "Mexico"), 
                trade_direction = "exports")
```

```{r, warning = FALSE, message = FALSE}
library(ggplot2)

# Apply polished col headers.
df <- ct_use_pretty_cols(df)

# Create plot.
ggplot(df, aes(Year, `Trade Value usd`, color = factor(`Partner Country`), 
               shape = factor(`Partner Country`))) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(min(df$Year), max(df$Year)), 
                     breaks = seq.int(min(df$Year), max(df$Year), 2)) +
  scale_color_manual(values = c("orange", "blue", "red"), 
                     name = "Destination\nCountry") +
  scale_shape_discrete(name = "Destination\nCountry") +
  labs(title = "Total Value (USD) of Chinese Exports, by Year") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

**Plot 2**: Plot the top eight destination countries/areas of Thai shrimp exports, by weight (KG), for 2007 - 2011.

```{r, echo = FALSE}
v_data_3 <- system.file("extdata", "vignette_data_3.rda", package = "comtradr")
if (!file.exists(v_data_3)) {
  stop("internal vignette data set '~/extdata/vignette_data_3.rda' not found", call. = FALSE)
}
load(v_data_3)
```
```{r, eval = FALSE}
# First, collect commodity codes related to shrimp.
shrimp_codes <- ct_commodity_lookup("shrimp", 
                                    return_code = TRUE, 
                                    return_char = TRUE)

# Comtrade api query.
df <- ct_search(reporters = "Thailand", 
                partners = "All", 
                trade_direction = "exports", 
                start_date = 2007, 
                end_date = 2011, 
                commod_codes = shrimp_codes)
```

```{r, warning = FALSE, message = FALSE}
library(ggplot2)
library(dplyr)

# Apply polished col headers.
df <- ct_use_pretty_cols(df)

# Create country specific "total weight per year" dataframe for plotting.
plotdf <- df %>% 
  group_by_(.dots = c("`Partner Country`", "Year")) %>% 
  summarise(kg = as.numeric(sum(`Net Weight kg`, na.rm = TRUE))) %>% 
  as_data_frame()

# Get vector of the top 8 destination countries/areas by total weight shipped 
# across all years, then subset plotdf to only include observations related 
# to those countries/areas.
top8 <- plotdf %>% 
  group_by(`Partner Country`) %>% 
  summarise(kg = as.numeric(sum(kg, na.rm = TRUE))) %>% 
  top_n(8, kg) %>%
  arrange(desc(kg)) %>% 
  .[["Partner Country"]]
plotdf <- plotdf %>% filter(`Partner Country` %in% top8)

# Create plots (y-axis is NOT fixed across panels, this will allow us to ID 
# trends over time within each country/area individually).
qplot(Year, kg, data = plotdf) + 
  geom_line(data = plotdf[plotdf$`Partner Country` %in% names(which(table(plotdf$`Partner Country`) > 1)), ]) + 
  xlim(min(plotdf$Year), max(plotdf$Year)) + 
  labs(title = "Weight (KG) of Thai Shrimp Exports, by Destination Area, 2007 - 2011") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.text = element_text(size = 7)) + 
  facet_wrap(~factor(`Partner Country`, levels = top8), scales = "free", nrow = 2, ncol = 4)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ct_search.R
\name{ct_search}
\alias{ct_search}
\title{Get UN Comtrade data}
\usage{
ct_search(
  reporters,
  partners,
  trade_direction = c("all", "imports", "exports", "re_imports", "re_exports"),
  freq = c("annual", "monthly"),
  start_date = "all",
  end_date = "all",
  commod_codes = "TOTAL",
  max_rec = NULL,
  type = c("goods", "services"),
  url = "https://comtrade.un.org/api/get?"
)
}
\arguments{
\item{reporters}{Country(s) of interest, as a character vector. Can either
be a vector of country names, or "All" to represent all countries.}

\item{partners}{Country(s) that have interacted with the reporter
country(s), as a character vector. Can either be a vector of country names,
or "All" to represent all countries.}

\item{trade_direction}{Indication of which trade directions on which to
focus, as a character vector. Must either be "all", or a vector containing
any combination of the following: "imports", "exports", "re_imports",
"re_exports". Default value is "all".}

\item{freq}{Time frequency of the returned results, as a character string.
Must be either "annual" or "monthly". Default value is "annual".}

\item{start_date}{Start date of a time period, or "all". Default value is
"all". See "details" for more info on valid input formats when not using
"all" as input.}

\item{end_date}{End date of a time period, or "all". Default value is
"all". See "details" for more info on valid input formats when not using
"all" as input.}

\item{commod_codes}{Character vector of commodity codes, or "TOTAL". Valid
commodity codes as input will restrict the query to only look for trade
related to those commodities, "TOTAL" as input will return all trade
between the indicated reporter country(s) and partner country(s). Default
value is "TOTAL".}

\item{max_rec}{Max number of records returned from each API call, as an
integer. If max_rec is set to NULL, then value is determined by whether or
not an API token has been registered. API cap without a token is 50000,
cap with a valid token is 250000. Default value is NULL. For details on
how to register a valid token, see \code{\link{ct_register_token}}.}

\item{type}{Type of trade, as a character string. Must be either "goods" or
"services". Default value is "goods".}

\item{url}{Base of the Comtrade url string, as a character string. Default
value is "https://comtrade.un.org/api/get?" and should mot be changed
unless Comtrade changes their endpoint url.}
}
\value{
Data frame of Comtrade shipping data.
}
\description{
Make queries to the UN Comtrade API, data is returned as a tidy data frame.
Comtrade is a DB hosted by the United Nations that houses country-level
shipping data. Full API docs can be found here:
\url{https://comtrade.un.org/data/doc/api/}
}
\details{
Basic rate limit restrictions listed below. For details on how to
 register a valid token, see \code{\link{ct_register_token}}. For API docs
 on rate limits, see \url{https://comtrade.un.org/data/doc/api/#Limits}
 \itemize{
 \item Without authentication token: 1 request per second, 100 requests
   per hour (each per IP address).
 \item With valid authentication token: 1 request per second, 10,000
   requests per hour (each per IP address or authenticated user).
 }

 In addition to these rate limits, the API imposes some limits on
 parameter combinations, they are listed below:
 \itemize{
 \item Between params "reporters", "partners", and the query date range (as
   dictated by the two params "start_date" and "end_date"), only one of
   these three may use the catch-all input "All".
 \item For the same group of three ("reporters", "partners", date range),
   if the input is not "All", then the maximum number of input values
   for each is five. For date range, if not using "All", then the
   "start_date" and "end_date" must not span more than five months or five
   years. There is one exception to this rule, if arg "freq" is "monthly",
   then a single year can be passed to "start_date" and "end_date" and the
   API will return all of the monthly data for that year.
 \item For param "commod_codes", if not using input "All", then the maximum
   number of input values is 20 (although "All" is always a valid input).
 }

 This function returns objects with metadata related to the API call that
 can be accessed via \code{\link{attributes}}. The metadata accessible is:
 \itemize{
 \item url: url of the API call.
 \item time_stamp: date-time of the API call.
 \item req_duration: total duration of the API call, in seconds.
 }

 For args \code{start_date} and \code{end_date}, if inputting a date (as
 opposed to the catch-all input "all"), valid input format is dependent on
 the input passed to arg \code{freq}. If \code{freq} is "annual",
 \code{start_date} and \code{end_date} must be either a string w/ format
 "yyyy" or "yyyy-mm-dd", or a year as an integer (so "2016", "2016-01-01",
 and 2016 would all be valid). If \code{freq} is "monhtly",
 \code{start_date} and \code{end_date} must be a string with format
 "yyyy-mm" or "yyyy-mm-dd" (so "2016-02" and "2016-02-01" would both be
 valid).
}
\examples{
\dontrun{
## Example API call number 1:
# All exports from China to South Korea, United States and Mexico over all
# years.
ex_1 <- ct_search(reporters = "China",
                  partners = c("Rep. of Korea", "USA", "Mexico"),
                  trade_direction = "exports")
nrow(ex_1)

## Example API call number 2:
# All shipments related to shrimp between Canada and all other countries,
# between 2011 and 2015.
# Perform "shrimp" query
shrimp_codes <- ct_commodity_lookup("shrimp",
                                    return_code = TRUE,
                                    return_char = TRUE)

# Make API call
ex_2 <- ct_search(reporters = "Canada",
                  partners = "All",
                  trade_direction = "all",
                  start_date = 2011,
                  end_date = 2015,
                  commod_codes = shrimp_codes)
nrow(ex_2)

# Access metadata
attributes(ex_2)$url
attributes(ex_2)$time_stamp
attributes(ex_2)$req_duration
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ct_use_pretty_cols.R
\name{ct_use_pretty_cols}
\alias{ct_use_pretty_cols}
\title{Use Pretty Column Headers}
\usage{
ct_use_pretty_cols(df)
}
\arguments{
\item{df}{data frame, Comtrade API data frame, returned from function
\code{\link{ct_search}}.}
}
\value{
data frame, input df with polish column headers.
}
\description{
Transform the column headers of return data from function
\code{\link{ct_search}} into a more "polished" set of column headers.
Intended for use with plots, publication tables, etc.
}
\examples{
\dontrun{
# Pull API data
df <- ct_search("Germany", "Canada")

# Use polished column names
df <- ct_use_pretty_cols(df)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comtradr.R
\docType{package}
\name{comtradr}
\alias{comtradr}
\title{Interface to the United Nations Comtrade API}
\description{
Interface with and extract data from the United Nations
  Comtrade API. Comtrade provides country level shipping data for a variety
  of commodities, these functions allow for easy API query and data returned
  as a tidy data frame.
}
\section{Package Vignette}{


\itemize{
  \item \url{../doc/comtradr-vignette.html}
}
}

\section{Documentation for the Comtrade API}{


\itemize{
  \item Main Comtrade Site \url{https://comtrade.un.org/}
  \item Comtrade Data Query Web GUI \url{https://comtrade.un.org/data/}
  \item Full API Documentation \url{https://comtrade.un.org/data/doc/api/}
}
}

\section{Development links}{


\itemize{
  \item \url{https://github.com/ChrisMuir/comtradr}
  \item Report bugs at \url{https://github.com/ChrisMuir/comtradr/issues}
}
}

\section{\code{comtradr} features the following functions}{

\itemize{
  \item \code{\link{ct_commodity_db_type}}
  \item \code{\link{ct_commodity_lookup}}
  \item \code{\link{ct_country_lookup}}
  \item \code{\link{ct_get_remaining_hourly_queries}}
  \item \code{\link{ct_get_reset_time}}
  \item \code{\link{ct_register_token}}
  \item \code{\link{ct_search}}
  \item \code{\link{ct_update_databases}}
  \item \code{\link{ct_use_pretty_cols}}
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_limit.R
\name{ct_register_token}
\alias{ct_register_token}
\title{Comtradr set API token}
\usage{
ct_register_token(token)
}
\arguments{
\item{token}{char string, valid API token.}
}
\value{
Set comtradr API token and update rate limits.
}
\description{
Function to set an API token for the UN Comtrade API. Details on tokens and
  rate limits can be found
  \url{https://comtrade.un.org/data/doc/api/#Authentication}
}
\examples{
\dontrun{
ct_register_token("some_valid_token_str")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_limit.R
\name{ct_get_remaining_hourly_queries}
\alias{ct_get_remaining_hourly_queries}
\title{Comtradr rate limit check}
\usage{
ct_get_remaining_hourly_queries()
}
\value{
numeric value, number of current queries left in the hour.
}
\description{
Get the remaining number of queries left in the current hour.
}
\examples{
ct_get_remaining_hourly_queries()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ct_commodity_lookup.R
\name{ct_commodity_lookup}
\alias{ct_commodity_lookup}
\title{UN Comtrade commodities database query}
\usage{
ct_commodity_lookup(
  search_terms,
  return_code = FALSE,
  return_char = FALSE,
  verbose = TRUE,
  ignore.case = TRUE,
  ...
)
}
\arguments{
\item{search_terms}{Commodity names or commodity codes, as a char or numeric
vector.}

\item{return_code}{Logical, if set to FALSE, the function will return a
set of commodity descriptions along with commodity codes (as a single
string for each match found), if set to TRUE it will return only the
commodity codes. Default value is FALSE.}

\item{return_char}{Logical, if set to FALSE, the function will return the
matches as a named list, if set to TRUE it will return them as a character
vector. Default value is FALSE.}

\item{verbose}{Logical, if set to TRUE, a warning message will print to
console if any of the elements of input "search_terms" returned no matches
(message will indicate which elements returned no data). Default is TRUE.}

\item{ignore.case}{logical, to be passed along to arg ignore.case within
\code{\link{grepl}}. Default value is TRUE.}

\item{...}{additional args to be passed along to \code{\link{grepl}}.}
}
\value{
A list or character vector of commodity descriptions and/or
 commodity codes that are matches with the elements of "search_terms".
}
\description{
The Comtrade API requires that searches for specific commodities be done
using commodity codes. This is a helper function for querying the
Comtrade commodity database. It takes as input a vector of
commodities or commodity codes. Output is a list or vector of commodity
descriptions or codes associated with the input search_terms. For use with
the UN Comtrade API, full API docs can be found at
\url{https://comtrade.un.org/data/doc/api/}
}
\details{
This function uses regular expressions (regex) to find matches
 within the commodity DB. This means it will treat as a match any commodity
 description that contains the input search term. For more on using regex
 within R, see this great tutorial by Gloria Li and Jenny Bryan
 \url{http://stat545.com/block022_regular-expression.html}
}
\examples{
# Look up commodity descriptions related to "halibut"
ct_commodity_lookup("halibut",
                    return_code = FALSE,
                    return_char = FALSE,
                    verbose = TRUE)

# Look up commodity codes related to "shrimp".
ct_commodity_lookup("shrimp",
                    return_code = TRUE,
                    return_char = FALSE,
                    verbose = TRUE)
}
\seealso{
\code{\link{grepl}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_limit.R
\name{ct_get_reset_time}
\alias{ct_get_reset_time}
\title{Comtradr rate limit time check}
\usage{
ct_get_reset_time(set = NULL)
}
\arguments{
\item{set}{logical, if \code{TRUE} and the current reset time is
\code{NULL}, set the reset time to be one hour from the current
\code{Sys.time}.}
}
\value{
date and time in which the hourly query limit will reset. Return is
  a "POSIXct" object (see \code{\link{DateTimeClasses}}).
}
\description{
Get the time in which the hourly limit will reset.
}
\examples{
ct_get_reset_time()

# Get minutes remaining until limit reset, as numeric value.
as.double(ct_get_reset_time() - Sys.time())
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ct_update_databases.R
\name{ct_update_databases}
\alias{ct_update_databases}
\title{Check for updates to country/commodity databases}
\usage{
ct_update_databases(
  force = FALSE,
  verbose = TRUE,
  commodity_type = c("HS", "HS1992", "HS1996", "HS2002", "HS2007", "HS2012", "HS2017",
    "SITC", "SITCrev1", "SITCrev2", "SITCrev3", "SITCrev4", "BEC", "EB02"),
  commodity_url = NULL,
  reporter_url = NULL,
  partner_url = NULL
)
}
\arguments{
\item{force}{logical, if TRUE, both the country and commodity databases
will be downloaded, regardless of the status of the DB's on file. Default
value is FALSE.}

\item{verbose}{logical, if TRUE, an update status message will be printed
to console. Default value is TRUE.}

\item{commodity_type}{Trade data classification scheme to use, see
"details" for a list of the valid inputs. Default value is "HS", which is
the default "type" of the commodity database on file upon install of
\code{comtradr}. Please note that if the value passed to this arg doesn't
match the values in variable "type" of the current commodity DB, then
this function will replace the current commodity DB with that of the
type specified by this arg. If you don't intend to change the type of the
current commodity DB, then no input for this arg is required. To see
the "type" of the current commodity DB, use
\code{\link{ct_commodity_db_type}}.}

\item{commodity_url}{Default value NULL, otherwise this should be the base
url of the Comtrade json data directory. Only necessary if the Comtrade
site changes from "https://comtrade.un.org/data/cache/". This partial
url string will have a commodity extension appended to it to create a
valid url. The commodity extension will be chosen based on the input to
arg \code{commodity_type}.}

\item{reporter_url}{Default value NULL, otherwise this should be a url as a
char string that points to the reporter areas JSON dataset on the Comtrade
website. Only necessary if the Comtrade site changes from
\url{https://comtrade.un.org/data/cache/reporterAreas.json}}

\item{partner_url}{Default value NULL, otherwise this should be a url as a
char string that points to the reporter areas JSON dataset on the Comtrade
website. Only necessary if the Comtrade site changes from
\url{https://comtrade.un.org/data/cache/partnerAreas.json}}
}
\value{
Updated database of commodities and countries.
}
\description{
Use of the Comtrade API requires access to the Comtrade countries database
and commodities database. The \code{comtradr} package keeps each DB saved
as a data frame in the package directory, as Comtrade makes updates to
these DB's infrequently (roughly once per year).
}
\details{
This function will check to see if Comtrade has made any updates to either
database. If an update is found, it will download the updated DB and save
it to the \code{comtradr} package directory, and update the DB for use
within the current R session.

The default for arg \code{commodity_type} is \code{HS}. Below is a
 list of all valid inputs with a very brief description for each, for more
 information on each of these types, see
 \url{https://comtrade.un.org/data/doc/api/#DataAvailabilityRequests}
 \itemize{
 \item \code{HS}: Harmonized System (HS), as reported
 \item \code{HS1992}: HS 1992
 \item \code{HS1996}: HS 1996
 \item \code{HS2002}: HS 2002
 \item \code{HS2007}: HS 2007
 \item \code{HS2012}: HS 2012
 \item \code{HS2017}: HS 2017
 \item \code{SITC}: Standard International Trade Classification (SITC), as
   reported
 \item \code{SITCrev1}: SITC Revision 1
 \item \code{SITCrev2}: SITC Revision 2
 \item \code{SITCrev3}: SITC Revision 3
 \item \code{SITCrev4}: SITC Revision 4
 \item \code{BEC}: Broad Economic Categories
 \item \code{EB02}: Extended Balance of Payments Services Classification
 }
}
\examples{
\dontrun{
ct_update_databases()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ct_commodity_db_type.R
\name{ct_commodity_db_type}
\alias{ct_commodity_db_type}
\title{Get current commodity database type}
\usage{
ct_commodity_db_type()
}
\value{
character vector of the "type" of the current commodity database.
}
\description{
Return the "type" of the current commodity database being used by
\code{comtradr}. For a complete list of the different commodity DB
types, see "details".
}
\details{
Below is a list of all of the commodity database "types", with a
 very brief description for each. For more information on each of these
 types, see
 \url{https://comtrade.un.org/data/doc/api/#DataAvailabilityRequests}
 \itemize{
 \item \code{HS}: Harmonized System (HS), as reported
 \item \code{HS1992}: HS 1992
 \item \code{HS1996}: HS 1996
 \item \code{HS2002}: HS 2002
 \item \code{HS2007}: HS 2007
 \item \code{HS2012}: HS 2012
 \item \code{SITC}: Standard International Trade Classification (SITC), as
   reported
 \item \code{SITCrev1}: SITC Revision 1
 \item \code{SITCrev2}: SITC Revision 2
 \item \code{SITCrev3}: SITC Revision 3
 \item \code{SITCrev4}: SITC Revision 4
 \item \code{BEC}: Broad Economic Categories
 \item \code{EB02}: Extended Balance of Payments Services Classification
 }
}
\examples{
ct_commodity_db_type()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comtradr.R
\docType{data}
\name{ct_pretty_cols}
\alias{ct_pretty_cols}
\title{"pretty" column headers for Comtrade API data.}
\format{
Named vector, with the polished column headers as the names, and
 the machine-readable column headers as the values. Each element is meant
 to be treated as a key-value pair. The function \code{\link{ct_search}}
 returns data with the machine-readable column headers by default.
}
\description{
Named vector of polished column headers, intended for use with plots,
publication tables, etc.
}
\examples{
data(ct_pretty_cols)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ct_country_lookup.R
\name{ct_country_lookup}
\alias{ct_country_lookup}
\title{UN Comtrade country database query}
\usage{
ct_country_lookup(
  search_terms,
  type = c("reporter", "partner"),
  ignore.case = TRUE,
  ...
)
}
\arguments{
\item{search_terms}{Char vector of country names.}

\item{type}{str, the country list to use for the search, valid inputs are
"reporter" and "partner".}

\item{ignore.case}{logical, to be passed along to arg ignore.case within
\code{\link{grepl}}. Default value is TRUE.}

\item{...}{additional args to be passed along to \code{\link{grepl}}.}
}
\value{
A character vector of country names that are complete or partial
 matches with any of the input country names.
}
\description{
Country names passed to the Comtrade API must have precise
spelling/capitalization. This is a helper function for querying the country
names/spelling used by Comtrade.. It takes as input a vector of
country names, output is any country names that contain any of the input
strings, using regex via the base function grepl.
For use with the UN Comtrade API, full API docs can be found at
\url{https://comtrade.un.org/data/doc/api/}
}
\details{
This function uses regular expressions (regex) to find matches
 within the country DB. This means it will treat as a match any country
 string that contains the input search term. For more on using regex
 within R, see this great tutorial by Gloria Li and Jenny Bryan
 \url{http://stat545.com/block022_regular-expression.html}
}
\examples{
# Look up all reporters that contain the terms "korea" and "vietnam"
ct_country_lookup(c("korea", "vietnam"), "reporter")
}
\seealso{
\code{\link{grepl}}
}
