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

