
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
