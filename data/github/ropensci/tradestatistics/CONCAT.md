
<!-- README.md is generated from README.Rmd. Please edit that file -->
<p>
<a href="https://www.digitalocean.com/">
<img src="https://opensource.nyc3.cdn.digitaloceanspaces.com/attribution/assets/PoweredByDO/DO_Powered_by_Badge_blue.svg" width="201px">
</a>
</p>

# Open Trade Statistics package <img src="svg/hexicon.svg" width=150 align="right" alt="sticker"/>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![R build
status](https://github.com/ropensci/tradestatistics/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/tradestatistics/actions?workflow=R-CMD-check)
[![CRAN
status](https://www.r-pkg.org/badges/version/tradestatistics)](https://cran.r-project.org/package=tradestatistics)
[![Coverage
status](https://codecov.io/gh/ropensci/tradestatistics/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tradestatistics?branch=master)
[![](https://badges.ropensci.org/274_status.svg)](https://github.com/ropensci/onboarding/issues/274)

[Open Trade Statistics](https://tradestatistics.io) is an effort to open
international trade data. `tradestatistics` provides an easy way to
obtain data from OTS by accessing its API.

This is what the package does:

![Data diagram](svg/data-diagram.svg)

Using `tradestatistics` package is all about efficiency, without this
package you could obtain the same data from the API at the expense of
using additional time and effort for the same results. As an API wrapper
and utility program this package makes data obtaining faster and easier
for you.

## Installation

``` r
# Install stable version from CRAN
install.packages("tradestatistics")

# Install stable version from GitHub
devtools::install_github("ropensci/tradestatistics")
```

## Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://docs.ropensci.org/tradestatistics/CODE_OF_CONDUCT.html).
By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# version 3.0.2

Updates
* Allows to obtain tables in Parquet format from the API, giving a speed-up
 of ~50% for the final user.
  
# version 3.0.1

Updates
* Adds section colors data for visualization, this is taken from the palette
  used in shiny.tradestatistics.io

# Version 3.0

Updates
* Removes all references to tables using communities or short names 
 (both unofficial), reflecting changes in the API
* The functionality remains the same, but now the end user functions don't
 add a 21-colors palette to the data (i.e. see the data section)

Data

* Switches from HS92  to HS12 to reflect product changes with less aggregation
* Drops any data from Harvard (communities and short product names) as these
 depend on using HS92 4 digits, therefore the color palettes were removed as
 these depended on the communities table
* The inflation data was trimmed to a window since the year 2000
* The commodities data now contains information for +5000 products instead of
 +1200 as the aggregation level changed in the API
* Adds RTAs and MFN tariffs for gravity modelling

# Version 2.0

Updates 

* Uses ISO codes as is (affects Aruba, Roumania, Timor-Leste, Antarctica, 
 Saint Barthelemy, Curacao, Sint Maarten and South Sudan)
  
# Version 1.0

Updates

* Reflects API changes with less aggregated data
* Follows UN COMTRADE notation (i.e. commodity instead of product)
* Does not impute data before hand, which is better for most of gravity models use cases
* Provides the data exactly as in the API, returning commodity level data to allow users to do their own aggregation
* Does not drop reference year with inflation adjustment (https://github.com/ropensci/tradestatistics/issues/38)
* Takes max and min available years from the API instead of hardcoded values (https://github.com/ropensci/tradestatistics/pull/39)

# Version 0.4.0

Updates

* Includes `yrpc-ga`, `yrpc-sa`, `yrc-ga` and `yr-sa` tables reflecting API updates
* Simplifies end-user functions a bit (i.e. removes `include_groups` option)
* Optimizes the code a bit, specially at the joins with tables in the package
* Fixes codes duplication when both product and group/community match for a search
* Includes both official and shortened section names

# Version 0.3.1

Updates

* Removes `yrp_short` option reflecting last DB changes

# Version 0.3

Updates

* Much improved coverage to detect almost any possible error
* Fixes case in inflation adjustment when year = reference year

# Version 0.2.8

Updates

* Adds caching (in memory or on disk) option
* Lists Daniela de los Santos and Elio Campitelli as new contributors
* Includes forwards and backwards testing for inflation adjustment
* Testing for in memory caching

# Version 0.2.7

Updates

* Adds feedback provided by Daniela de los Santos
* Now ots_create_tidy_data() has both reporter and partner set to "all" by default

# Version 0.2.5

Updates

* Added dependency on R >= 3.5.0 because serialized objects in serialize/load version 3 cannot be read in older versions of R
* Minimal changes in `ots_create_tidy_data()` to allow multiple countries as arguments, in line with API changes from September 2019

# Version 0.2.4

Updates

* Removes `product_code_length`
* The API was updated with simplified parameters and 2018 data

# Version 0.2.3

Updates

* Fixtures for testthat evaluation

Fixes

* Specific Windows error during check

# Version 0.2.2

Adds

* Inflation data
* Inflation adjustment function
* Minor changes in vignettes

# Version 0.2.1

Fixes

* Consistent use of colour vs color, color is used from now
* Fixed available tables description
* Adds `yrp_short` to available tables
* Adds `use_localhost` option for our own server or users who want to clone the
  database locally, therefore avoid having a separate branh for server installation
  
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# Contributing to tradestatistics

This outlines how to propose a change to tradestatistics. For more detailed
info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  New code should follow the tidyverse [style guide](https://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the tradestatistics project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See tidyverse [development contributing guide](https://rstd.io/tidy-contrib)
for further details.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/"
)
```

<p>
  <a href="https://www.digitalocean.com/">
    <img src="https://opensource.nyc3.cdn.digitaloceanspaces.com/attribution/assets/PoweredByDO/DO_Powered_by_Badge_blue.svg" width="201px">
  </a>
</p>

# Open Trade Statistics package <img src="svg/hexicon.svg" width=150 align="right" alt="sticker"/>

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![R build status](https://github.com/ropensci/tradestatistics/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/tradestatistics/actions?workflow=R-CMD-check)
[![CRAN status](https://www.r-pkg.org/badges/version/tradestatistics)](https://cran.r-project.org/package=tradestatistics)
[![Coverage status](https://codecov.io/gh/ropensci/tradestatistics/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tradestatistics?branch=master)
[![](https://badges.ropensci.org/274_status.svg)](https://github.com/ropensci/onboarding/issues/274)

[Open Trade Statistics](https://tradestatistics.io) is an effort to open international trade data. `tradestatistics` provides an easy way to obtain data from OTS by accessing its API.

This is what the package does:

![Data diagram](svg/data-diagram.svg)

Using `tradestatistics` package is all about efficiency, without this package you could obtain the same data from the API at the expense of using additional time and effort for the same results. As an API wrapper and utility program this package makes data obtaining faster and easier for you.

## Installation

```{r, eval = FALSE}
# Install stable version from CRAN
install.packages("tradestatistics")

# Install stable version from GitHub
devtools::install_github("ropensci/tradestatistics")
```

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/tradestatistics/CODE_OF_CONDUCT.html).
By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Basic usage"
author: "Mauricio Vargas S."
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to use this package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  cache = FALSE,
  collapse = TRUE,
  message = FALSE,
  comment = "#>"
)

datatable <- function(x) {
  DT::datatable(x,
    extensions = "FixedColumns",
    options = list(
      pageLength = 5,
      dom = 'Bfrtip',
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 2, rightColumns = 1)
    )
)}
```

# Introduction

This vignette explains the functions within this package. The idea is to show how this package simplifies obtaining data from (api.tradestatistics.io)[https://api.tradestatistics.io].

To improve the presentation of the tables I shall use `DT` and the `datatable()` function besides `tradestatistics`.
```{r pkgs}
library(tradestatistics)
library(DT)
```

# Package data

## Available tables

Provided that this package obtains data from an API, it is useful to know which tables can be accessed:

```{r tables, eval = T}
datatable(ots_tables)
```

You might notice the tables have a pattern. The letters indicate the presence of columns that account for the level of detail in the data:

* `y`: *y*ear column.
* `r`: *r*eporter column
* `p`: *p*artner column
* `c`: *c*ommodity column

The most aggregated table is `yr` which basically says how many dollars each country exports and imports for a given year.

The less aggregated table is `yrpc` which says how many dollars of each of the 1,242 commodities from the Harmonized System each country exports to other countries and imports from other countries.

For the complete detail you can check [tradestatistics.io](https://tradestatistics.io).

## Country codes

The Package Functions section explains that you don't need to memorize all ISO codes. The functions within this package are designed to match strings (i.e. "United States" or "America") to valid ISO codes (i.e. "USA").

Just as a reference, the table with all valid ISO codes can be accessed by running this:

```{r countries, eval = T}
datatable(ots_countries)
```

## Commodity codes

The Package Functions section explains that you don't need to memorize all HS codes. The functions within this package are designed to match strings (i.e. "apple") to valid HS codes (i.e. "0808").

```{r commodities, eval = T}
datatable(ots_commodities)
```

## Inflation data

This table is provided to be used with `ots_inflation_adjustment()`.

```{r inflation, eval = T}
datatable(ots_inflation)
```

# Package functions

## Country code

The end user can use this function to find an ISO code by providing a country name. This works by implementing partial search.

Basic examples:
```{r country_code}
# Single match with no replacement
datatable(ots_country_code("Chile"))

# Single match with replacement
datatable(ots_country_code("America"))

# Double match with no replacement
datatable(ots_country_code("Germany"))
```

The function `ots_country_code()` is used by `ots_create_tidy_data()` in a way that you can pass parameters like `ots_create_tidy_data(... reporters = "Chile" ...)` and it will automatically replace your input for a valid ISO in case there is a match. This will be covered in detail in the Trade Data section.

## Commodity code

The end user can find a code or a set of codes by looking for keywords for commodities or groups. The function `ots_commodity_code()` allows to search from the official commodities and groups in the Harmonized system:
```{r commodity_code2}
datatable(ots_commodity_code(commodity = " ShEEp ", group = " mEaT "))
```

## Trade data

This function downloads data for a single year and needs (at least) some filter parameters according to the query type.

Here we cover aggregated tables to describe the usage. Note: here you may skip the `use_localhost = FALSE` argument.

### Bilateral trade at commodity level (Year - Reporter - Partner - Commodity Code)

If we want Chile-Argentina bilateral trade at community level in 2019:
```{r yrpc1, eval = T}
yrpc <- ots_create_tidy_data(
  years = 2019,
  reporters = "chl",
  partners = "arg",
  table = "yrpc",
  use_localhost = FALSE
)

datatable(yrpc)
```

We can pass two years or more, several reporters/partners, and filter by commodities with exact codes or code matching based on keywords:
```{r yrpc2, eval = T}
# Note that here I'm passing Peru and not per which is the ISO code for Peru
# The same applies to Brazil
yrpc2 <- ots_create_tidy_data(
  years = 2018:2019,
  reporters = c("chl", "Peru", "bol"),
  partners = c("arg", "Brazil"),
  commodities = c("01", "food"),
  table = "yrpc",
  use_localhost = FALSE
)
datatable(yrpc2)
```

The `yrpc` table returns some fields that deserve an explanation which can be seen at [tradestatistics.io](https://tradestatistics.io). This example is interesting because "01" return a set of commodities (all commodities starting with 01, which is the commodity group "Animals; live"), but "food" return all commodities with a matching description ("1601", "1806", "1904", etc.). In addition, not all the requested commodities are exported from each reporter to each partner, therefore a warning is returned.

### Bilateral trade at aggregated level (Year - Reporter - Partner)

If we want Chile-Argentina bilateral trade at aggregated level in 2018 and 2019:

```{r yrp3, eval = T}
yrp <- ots_create_tidy_data(
  years = 2018:2019,
  reporters = c("chl", "per"),
  partners = "arg",
  table = "yrp",
  use_localhost = FALSE
)

datatable(yrp)
```

This table accepts different years, reporters and partners just like `yrpc`.

### Reporter trade at commodity level (Year - Reporter - Commodity Code) 

If we want Chilean trade at commodity level in 2019 with respect to commodity "010121" which means "Horses; live, pure-bred breeding animals":
```{r yrc2, eval = T}
yrc <- ots_create_tidy_data(
  years = 2019,
  reporters = "chl",
  commodities = "010121",
  table = "yrc",
  use_localhost = FALSE
)

datatable(yrc)
```

This table accepts different years, reporters and commodity codes just like `yrpc`.

All the variables from this table are documented at [tradestatistics.io](https://tradestatistics.io).

### Reporter trade at aggregated level (Year - Reporter)

If we want the aggregated trade of Chile, Argentina and Peru in 2018 and 2019:
```{r yr2, eval = T}
yr <- ots_create_tidy_data(
  years = 2018:2019,
  reporters = c("chl", "arg", "per"),
  table = "yr",
  use_localhost = FALSE
)

datatable(yr)
```

This table accepts different years and reporters just like `yrpc`.

All the variables from this table are documented at [tradestatistics.io](https://tradestatistics.io).

### Commodity trade at aggregated level (Year - Commodity Code)

If we want all commodities traded in 2019:
```{r yc1, eval = T}
yc <- ots_create_tidy_data(
  years = 2019,
  table = "yc",
  use_localhost = FALSE
)

datatable(yc)
```

If we want the traded values of the commodity "010121" which means "Horses; live, pure-bred breeding animals" in 2019:
```{r yc2, eval = T}
yc2 <- ots_create_tidy_data(
  years = 2019,
  commodities = "010121",
  table = "yc",
  use_localhost = FALSE
)

datatable(yc2)
```

This table accepts different years just like `yrpc`.

## Inflation adjustment

Taking the `yr` table from above, we can use `ots_inflation_adjustment()` to convert dollars from 2019 to dollars of 2000:

```{r}
inflation <- ots_inflation_adjustment(yr, reference_year = 2000)
datatable(inflation)
```
\name{ots_countries}
\alias{ots_countries}
\title{A table of official country names, ISO-3 codes and other metadata}
\docType{data}
\description{
Provides official codes taken from the United Nations official sources.
This data is used by the functions provided within this package to validate
user parameters and add both product and country text columns to the data,
therefore reducing the number of API calls and the time to generate the
requested data.
}
\usage{ots_countries}
\format{
  A data frame with 254 observations on the following 6 variables.
  \describe{
    \item{\code{country_iso}}{ISO code of the country (e.g. "chl" means Chile)}
    \item{\code{country_name_english}}{Country name (e.g. Germany)}
    \item{\code{country_fullname_english}}{Country name with indications (e.g. Germany (former Federal Republic of Germany until 1990))}
    \item{\code{continent_id}}{Numeric id of the continent where the country belongs to}
    \item{\code{continent}}{Continent where the country belongs to}
    \item{\code{eu28_member}}{Dummy variable such that 1 means "belongs to EU-28 group" and 0 otherwise}
\  }
}
\examples{
ots_countries
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_cache.R
\name{ots_cache}
\alias{ots_cache}
\title{Caching wrapper to reduce API calls (internal)}
\usage{
ots_cache(use_cache, file, ...)
}
\arguments{
\item{use_cache}{Logical to save and load from cache. If \code{TRUE}, the results will be cached in memory
if \code{file} is \code{NULL} or on disk if `file` is not \code{NULL}.}

\item{file}{Character with the full file path to save the data.}

\item{...}{Additional parameters inherited from \code{ots_create_tidy_data()}.}
}
\description{
Eases saving the data downloaded from \code{api.tradestatistics.io}
and prevents \code{ots_read_from_api()} from downloading the same twice.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tradestatistics-package.R
\docType{package}
\name{tradestatistics-package}
\alias{tradestatistics}
\alias{tradestatistics-package}
\title{tradestatistics: Open Trade Statistics API Wrapper and Utility Program}
\description{
Access 'Open Trade Statistics' API from R to download international trade data.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/tradestatistics/}
  \item Report bugs at \url{https://github.com/ropensci/tradestatistics/issues}
}

}
\author{
\strong{Maintainer}: Mauricio Vargas \email{mvargas@dcc.uchile.cl} (\href{https://orcid.org/0000-0003-1017-7574}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Joshua Kunst (contributed to different parts of the pre-release code) [contributor]
  \item Alexey Kravchenko (reviewed 2021 version of the API) [contributor]
  \item Emma Mendelsohn (updated the functions to take available years from the API instead of hardcoded values) [contributor]
  \item Natalia de los Santos (proposed improvements to default parameters) [contributor]
  \item Elio Campitelli (wrote parts of the client-side caching function) [contributor]
  \item Emily Riederer (reviewed the package for rOpenSci, see https://github.com/ropensci/onboarding/issues/274) [reviewer]
  \item Mark Padgham (reviewed the package for rOpenSci, see https://github.com/ropensci/onboarding/issues/274) [reviewer]
  \item Amanda Dobbyn (reviewed a previous package that evolved into the current package for rOpenSci, see https://github.com/ropensci/onboarding/issues/217) [reviewer]
  \item Jorge Cimentada (reviewed a previous package that evolved into the current package for rOpenSci, see https://github.com/ropensci/onboarding/issues/217) [reviewer]
  \item  UN Comtrade [data contributor]
  \item  The World Bank [data contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_create_tidy_data.R
\name{ots_create_tidy_data_memoised}
\alias{ots_create_tidy_data_memoised}
\title{Downloads and processes the data from the API to return a human-readable tibble (memoised, internal)}
\usage{
ots_create_tidy_data_memoised(
  years = 2018,
  reporters = "usa",
  partners = "all",
  commodities = "all",
  table = "yr",
  max_attempts = 5,
  use_localhost = FALSE
)
}
\description{
A composition of \code{ots_create_tidy_data_unmemoised()} and \code{memoise()} for caching the output
}
\keyword{internal}
\name{ots_inflation}
\alias{ots_inflation}
\docType{data}
\title{
A table with world weigthed mean inflation since 2000
}
\description{
Provides year to year inflations value to be applied as a conversion rate
to express dollars of year Y1 as dollars of year Y2. This dataset is provided
to be used with \code{ots_inflation_adjustment} that converts units forwards and
backwards in time.
}
\usage{data("ots_inflation")}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{\code{from}}{Integer values in the range 2000-2019}
    \item{\code{to}}{Integer values in the range 2001-2020}
    \item{\code{conversion_factor}}{Numeric value expressed as one plus 1-year 
     inflation}
  }
}
\examples{
ots_inflation
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_inflation_adjustment.R
\name{ots_inflation_adjustment}
\alias{ots_inflation_adjustment}
\title{Expresses tidy data from the API in dollars of a reference year}
\usage{
ots_inflation_adjustment(trade_data = NULL, reference_year = NULL)
}
\arguments{
\item{trade_data}{A tibble obtained by using ots_create_tidy_data.
Default set to \code{NULL}.}

\item{reference_year}{Year contained within the years specified in
api.tradestatistics.io/year_range (e.g. \code{2010}).
Default set to \code{NULL}.}
}
\description{
Uses inflation records from The World Bank to
convert trade records and express them in dollars of the same year.
}
\examples{
\dontrun{
# The next example can take more than 5 seconds to compute,
# so this is shown without evaluation according to CRAN rules

# Convert dollars of 2010 to dollars of 2000
d <- ots_create_tidy_data(years = 2010, reporters = "chl", partners = "chn")
ots_inflation_adjustment(trade_data = d, reference_year = 2000)
}
}
\keyword{functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_strings_processing.R
\name{ots_commodity_code}
\alias{ots_commodity_code}
\title{String matching of official commodity/group names and Harmonized System (HS) codes
according to the United Nations nomenclature}
\usage{
ots_commodity_code(commodity = NULL, group = NULL)
}
\arguments{
\item{commodity}{A text string such as "Animals", "COPPER" or "fruits".}

\item{group}{A text string such as "meat", "FISH" or "Dairy".}
}
\value{
A tibble with all possible matches (no uppercase distinction)
showing the commodity name and commodity code
}
\description{
Takes a text string and searches within the
package data for all matching commodity codes in the context of valid API
commodity codes.
}
\examples{
ots_commodity_code(commodity = "ANIMALS ")
ots_commodity_code(group = "  fish")
ots_commodity_code(commodity = "Milk", group = "Dairy")
}
\keyword{functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_read_from_api.R
\name{ots_read_from_api}
\alias{ots_read_from_api}
\title{Reads data from the API (internal function)}
\usage{
ots_read_from_api(
  year = NULL,
  reporter_iso = NULL,
  partner_iso = NULL,
  commodity_code = "all",
  table = "yr",
  max_attempts = 5,
  use_localhost = FALSE
)
}
\arguments{
\item{year}{Year contained within the years specified in
api.tradestatistics.io/year_range (e.g. \code{1980}).
Default set to \code{NULL}.}

\item{reporter_iso}{ISO code for reporter country (e.g. \code{"chl"}). Default set to \code{"all"}.}

\item{partner_iso}{ISO code for partner country (e.g. \code{"chl"}). Default set to \code{"all"}.}

\item{commodity_code}{HS code (e.g. \code{0101} or \code{01}) to filter commodities.
Default set to \code{"all"}.}

\item{table}{Character string to select the table to obtain the data. Default set to \code{yr}
(Year - Reporter).}

\item{max_attempts}{Number of attempts to retry in case of data retrieving failure.
Default set to \code{5}.}

\item{use_localhost}{Logical to determine if the base URL shall be localhost instead
of api.tradestatistics.io. Default set to \code{FALSE}.}
}
\description{
Accesses \code{api.tradestatistics.io} and
performs different API calls to return \code{data.frames} by reading \code{JSON} data
}
\examples{
\dontrun{
# The next examples can take more than 5 seconds to compute,
# so these are shown without evaluation according to CRAN rules

# Run `countries` to display the full table of countries

# What does Chile export to China? (1980)
ots_read_from_api(year = 1980, reporter_iso = "chl", partner_iso = "chn")

# What can we say about chilean Horses export? (1980)
ots_read_from_api(year = 1980, commodity_code = "0101", table = "yc")
ots_read_from_api(year = 1980, reporter_iso = "chl", commodity_code = "0101", table = "yrc")
ots_read_from_api(
  year = 1980, reporter_iso = "chl", partner_iso = "arg", commodity_code = "0101",
  table = "yrpc"
)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_strings_processing.R
\name{ots_country_code}
\alias{ots_country_code}
\title{String matching of official country names and ISO-3 codes according to
the United Nations nomenclature}
\usage{
ots_country_code(countryname = NULL)
}
\arguments{
\item{countryname}{A text string such as "Chile", "CHILE" or "CHL".}
}
\value{
A single character if there is a exact match (e.g.
\code{ots_country_code("Chile")}) or a tibble in case of multiple matches
(e.g. \code{ots_country_code("Germany")})
}
\description{
Takes a text string and searches within the
package data for a country code in the context of valid API country codes.
}
\examples{
ots_country_code("Chile ")
ots_country_code("america")
ots_country_code("UNITED  STATES")
ots_country_code(" united_")
}
\keyword{functions}
\name{ots_commodities}
\alias{ots_commodities}
\title{A table of official commodity names from the Harmonized System rev 2007
 (HS07, also known as H3)}
\docType{data}
\description{
Provides official commodity, group and section codes and names taken from the 
United Nations official sources. This data is used by the functions provided
within this package to complement the data obtained from the API.
}
\usage{ots_commodities}
\format{
  A data frame with 5151 observations on the following 6 variables.
  \describe{
    \item{\code{commodity_code}}{Code of every commodity (e.g. 010110)}
    \item{\code{commodity_fullname_english}}{HS commodity names (e.g. 'Horses, asses, mules and hinnies; live, pure-bred breeding animals')}
    \item{\code{group_code}}{Group code (e.g. 01)}
    \item{\code{group_fullname_english }}{Group name (e.g. 'Animals; live')}
    \item{\code{section_code}}{Section code (e.g. 01)}
    \item{\code{section_fullname_english }}{Section name (e.g. 'Live animals and animal products')}
  }
}
\examples{
ots_commodities
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_create_tidy_data.R
\name{ots_create_tidy_data}
\alias{ots_create_tidy_data}
\title{Downloads and processes the data from the API to return a human-readable tibble}
\usage{
ots_create_tidy_data(
  years = 2019,
  reporters = "all",
  partners = "all",
  commodities = "all",
  table = "yr",
  max_attempts = 5,
  use_cache = FALSE,
  file = NULL,
  use_localhost = FALSE
)
}
\arguments{
\item{years}{Year contained within the years specified in
api.tradestatistics.io/year_range (e.g. \code{c(2002,2004)}, \code{c(2002:2004)} or \code{2002}).
Default set to \code{2019}.}

\item{reporters}{ISO code for reporter country (e.g. \code{"chl"}, \code{"Chile"} or
\code{c("chl", "Peru")}). Default set to \code{"all"}.}

\item{partners}{ISO code for partner country (e.g. \code{"chl"}, \code{"Chile"} or
\code{c("chl", "Peru")}). Default set to \code{"all"}.}

\item{commodities}{HS commodity codes (e.g. \code{"0101"}, \code{"01"} or search
matches for \code{"apple"})
to filter commodities. Default set to \code{"all"}.}

\item{table}{Character string to select the table to obtain the data.
Default set to \code{yr} (Year - Reporter).
Run \code{ots_tables} in case of doubt.}

\item{max_attempts}{How many times to try to download data in case the
API or the internet connection fails when obtaining data. Default set
to \code{5}.}

\item{use_cache}{Logical to save and load from cache. If \code{TRUE}, the results will be cached in memory
if \code{file} is \code{NULL} or on disk if `file` is not \code{NULL}. Default set to \code{FALSE}.}

\item{file}{Optional character with the full file path to save the data. Default set to \code{NULL}.}

\item{use_localhost}{Logical to determine if the base URL shall be localhost instead
of api.tradestatistics.io. Default set to \code{FALSE}.}
}
\value{
A tibble that describes bilateral trade metrics (imports,
exports, trade balance and relevant metrics
such as exports growth w/r to last year) between a \code{reporter}
and \code{partner} country.
}
\description{
Accesses \code{api.tradestatistics.io} and
performs different API calls to transform and return tidy data.
}
\examples{
\dontrun{
# The next examples can take more than 5 seconds to compute,
# so these are just shown without evaluation according to CRAN rules

# Run `ots_countries` to display the full table of countries
# Run `ots_commodities` to display the full table of commodities

# What does Chile export to China? (2002)
ots_create_tidy_data(years = 2002, reporters = "chl", partners = "chn")

# What can we say about Horses export in Chile and the World? (2002)
ots_create_tidy_data(years = 2002, commodities = "010110", table = "yc")
ots_create_tidy_data(years = 2002, reporters = "chl", commodities = "010110", table = "yrc")

# What can we say about the different types of apples exported by Chile? (2002)
ots_create_tidy_data(years = 2002, reporters = "chl", commodities = "apple", table = "yrc")
}
}
\keyword{functions}
\name{ots_sections_colors}
\alias{ots_sections_colors}
\title{A table of official section names from the Harmonized System rev 2007
 (HS07, also known as H3) and unofficial colors to ease visualization}
\docType{data}
\description{
Provides official section names taken from the United Nations official sources
but the colors are absolutely unofficial and based of what I consider a good
palette for 22 sections (21 + 1 for unspecified products). This data is not used
by the functions provided within this package and is provided as reference.
}
\usage{ots_sections_colors}
\format{
  A data frame with 22 observations on the following 2 variables.
  \describe{
    \item{\code{section_fullname_english }}{Section name (e.g. 'Live animals and animal products')}
    \item{\code{section_color}}{Section hex color (e.g. "#74c0e2")}
  }
}
\examples{
ots_sections_colors
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots_create_tidy_data.R
\name{ots_create_tidy_data_unmemoised}
\alias{ots_create_tidy_data_unmemoised}
\title{Downloads and processes the data from the API to return a human-readable tibble (unmemoised, internal)}
\usage{
ots_create_tidy_data_unmemoised(
  years = 2018,
  reporters = "usa",
  partners = "all",
  commodities = "all",
  table = "yr",
  max_attempts = 5,
  use_localhost = FALSE
)
}
\description{
A separation of \code{ots_create_tidy_data()} for making caching optional.
}
\keyword{internal}
\name{ots_tables}
\alias{ots_tables}
\title{Available tables in the API}
\docType{data}
\description{
A table describing existing API tables with both description and source.
This data is used by the functions provided within this package to validate
user parameters.
}
\usage{ots_tables}
\format{
  A data frame with 12 observations on the following 3 variables.
  \describe{
    \item{\code{table}}{Table name}
    \item{\code{description}}{Description of table contents}
    \item{\code{source}}{Source for the data (OTS tables are processed after UN Comtrade raw data)}
  }
}
\examples{
ots_tables
}
\keyword{datasets}
