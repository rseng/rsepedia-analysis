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
getCRUCLdata: Use and Explore CRU CL v. 2.0 Climatology Elements in R
================

<!-- badges: start -->
[![tic](https://github.com/ropensci/getCRUCLdata/workflows/tic/badge.svg?branch=main)](https://github.com/ropensci/getCRUCLdata/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/getCRUCLdata/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/getCRUCLdata?branch=main)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.466812.svg)](https://doi.org/10.5281/zenodo.466812)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/getCRUCLdata)](https://cran.r-project.org/package=getCRUCLdata)
[![JOSS status](http://joss.theoj.org/papers/10.21105/joss.00230/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00230)
[![](https://badges.ropensci.org/96_status.svg)](https://github.com/ropensci/software-review/issues/96)
<!-- badges: end -->

Author/Maintainer: Adam Sparks

## Introduction to *getCRUCLdata*

The *getCRUCLdata* package provides functions that automate importing CRU CL v. 2.0 climatology data into R, facilitate the calculation of minimum temperature and maximum temperature, and formats the data into a data frame or a [base::list()] of [raster::stack()] objects for use.

CRU CL v. 2.0 data are a gridded climatology of 1961-1990 monthly means released in 2002 and cover all land areas (excluding Antarctica) at 10 arc minutes (0.1666667 degree) resolution.
For more information see the description of the data provided by the University of East Anglia Climate Research Unit (CRU), <https://crudata.uea.ac.uk/cru/data/hrg/tmc/readme.txt>.

## Changes to original CRU CL v. 2.0 data

This package automatically converts elevation values from kilometres to metres.

This package crops all spatial outputs to an extent of ymin = -60, ymax = 85, xmin = -180, xmax = 180.
Note that the original wind data include land area for parts of Antarctica.

# Quick Start

## Install

_getCRUCLdata_ is not available from CRAN.
You can install it from GitHub as follows.

``` r
if (!require("remotes")) {
  install.packages("remotes")
}

install_github("ropensci/getCRUCLdata", build_vignettes = TRUE)
```

Or you can install it from the rOpenSci Universe.

``` r
# Enable the rOpenSci universe
options(repos = c(
    rOpenSci = "https://ropensci.r-universe.dev",
    CRAN = "https://cloud.r-project.org""))
# Install the package
install.packages("getCRUCLdata")
```

-----

# Documentation

For complete documentation see the package website: <https://docs.ropensci.org/getCRUCLdata/>.

# Meta

## CRU CL v. 2.0 reference and abstract

> Mark New (1,\*), David Lister (2), Mike Hulme (3), Ian Makin (4)

> A high-resolution data set of surface climate over global land areas 
> Climate Research, 2000, Vol 21, pg 1-25

> 1)  School of Geography and the Environment, University of Oxford,
>     Mansfield Road, Oxford OX1 3TB, United Kingdom  
> 2)  Climatic Research Unit, and (3) Tyndall Centre for Climate Change
>     Research, both at School of Environmental Sciences, University of
>     East Anglia, Norwich NR4 7TJ, United Kingdom  
> 3)  International Water Management Institute, PO Box 2075, Colombo,
>     Sri Lanka

> **ABSTRACT:** We describe the construction of a 10-minute
> latitude/longitude data set of mean monthly surface climate over
> global land areas, excluding Antarctica. The climatology includes 8
> climate elements - precipitation, wet-day frequency, temperature,
> diurnal temperature range, relative humidity,sunshine duration, ground
> frost frequency and windspeed - and was interpolated from a data set
> of station means for the period centred on 1961 to 1990. Precipitation
> was first defined in terms of the parameters of the Gamma
> distribution, enabling the calculation of monthly precipitation at any
> given return period. The data are compared to an earlier data set at
> 0.5 degrees latitude/longitude resolution and show added value over
> most regions. The data will have many applications in applied
> climatology, biogeochemical modelling, hydrology and agricultural
> meteorology and are available through the School of Geography Oxford
> (<http://www.geog.ox.ac.uk>), the International Water Management
> Institute “World Water and Climate Atlas” (<https://www.iwmi.cgiar.org/>) and
> the Climatic Research Unit (<http://www.cru.uea.ac.uk>).

## Contributors

  - [Adam H. Sparks](https://github.com/adamhsparks)

## Other

  - Please [report any issues or
    bugs](https://github.com/ropensci/getCRUCLdata/issues).

  - License: MIT

  - Get citation information for *getCRUCLdata* in R typing `citation(package = "getCRUCLdata")`

  - Please note that the *getCRUCLdata* project is released with a
  [Contributor Code of Conduct](https://github.com/ropensci/getCRUCLdata/blob/main/CONDUCT.md).
  By participating in the *getCRUCLdata* project you agree to abide by its terms.
# getCRUCLdata (development version)

# getCRUCLdata 0.3.2

## Minor changes

- Correct link redirects in README

- Correct formatting in documentation

- Precompile main vignette

- Add second vignette to illustrate advanced usage

- Remove _pkgdown_ from Suggests

- Move CI to GitHub Actions

# getCRUCLdata 0.3.1

## Bug fixes

- Fix bug in documentation that prevented example from working

## Minor changes

- Update URL in DESCRIPTION file

--------------------------------------------------------------------------------

# getCRUCLdata 0.3.0

## Major changes

- Remove Imports for _dplyr_, _tibble_ and _tidyr_ to lessen dependencies

- Remove Suggests for _readr_ and _sp_

- Enhance documentation

## Bug fixes

- Update tests that spuriously failed on some systems due to tolerances

- Update package to follow CRAN policies

--------------------------------------------------------------------------------

# getCRUCLdata 0.2.5

## Minor changes

- Removes startup message, instead placing information in CITATION file

- Reorganises internal functions consolidating functions all in a single file
and following a standard naming scheme for all internal functions

--------------------------------------------------------------------------------

# getCRUCLdata 0.2.4

## Bug fixes

- Fix bug where `tmp` and `dtr` could not be returned with `tmn` or `tmx` raster
stacks

- Move `rappdirs` to SUGGESTS to fix NOTEs on
https://cran.rstudio.com/web/checks/check_results_getCRUCLdata.html

## Minor changes

- Fix documentation formatting issues

- Enhance `stop` messages for user, just print message, not the function that
called it to clarify

--------------------------------------------------------------------------------

# getCRUCLdata 0.2.3

## Bug fixes

- Fix missing import for `rappdirs`

## Minor changes

- Remove the use of `plyr` in tests

--------------------------------------------------------------------------------

# getCRUCLdata 0.2.2

## Bug fixes

- Fix incorrect ORCID entry author field

--------------------------------------------------------------------------------
# getCRUCLdata 0.2.1

## Minor changes

- Fix ORCID entry in DESCRIPTION per CRAN maintainer's request

- Remove Scott as contributor, the code contributed has been removed

--------------------------------------------------------------------------------

# getCRUCLdata 0.2.0

## Major changes

- Use _hoardr_ for managing cached files

- Fixed a bug where the file cache was not in the proper subdirectory. The file
cache has moved to the proper location in a `R/getCRUCLdata` location rather
than `getCRUCLdata`. You may wish to move files externally to R in order to keep
them in the cache where the package will find them

- Use `lapply` in place of `purrr::map`, _purrr_ is no longer imported

## Minor changes

- Correct documentation where examples pointed to a non-existent list

## Deprecated functions

`CRU_cache_list()` now superseded by`manage_cache$list()`
`CRU_cache_details()` now superseded by `manage_cache$details()`
`CRU_cache_delete()` now superseded by `manage_cache$delete()`
`CRU_cache_delete_all()` now superseded by `manage_cache$delete_all()`

--------------------------------------------------------------------------------

# getCRUCLdata 0.1.10

## Major changes

- Add startup message regarding data source, use and citation

- Include Scott Chamberlain as copyright holder and contributor for file
caching functionality

# getCRUCLdata 0.1.9

## Bug fixes

- Fix issues in cached file management where files were not properly handled

# getCRUCLdata 0.1.8

## Bug Fixes

- Fix bug where `cache` was not specified in internal function, `.set_cache()`,
this caused either of the functions fetching data from CRU to fail

- Fix bug where `cache` directory could not be created on Windows OS machines

- Fix bug where tmx was returned when *either* tmn *or* tmx was requested for
data frame, tmn now returned when requested and tmx now returned when requested.
Raster stacks were not affected by this bug

## Minor Changes

- Replaced `for f in 1:length()` with `for f in seq_along()` for better
programming practices

--------------------------------------------------------------------------------

# getCRUCLdata 0.1.7

## Minor Changes

- Use `file.path` in place of `paste0`

## Bug Fixes

- Fix bug where `rappdirs::user_config_dir()` was incorrectly used in place of
`rappdirs::user_cache_dir()`

--------------------------------------------------------------------------------

# getCRUCLdata 0.1.6

## Minor Changes

- Use _purrr_ in place of _plyr_ functions

- Update DESCRIPTION file to be more complete

- Remove use of "%>%" in functions and remove _magrittr_ import

## Bug Fixes

- Fix bugs in CITATION file

- Format NEWS.md to be more markdown standards compliant

--------------------------------------------------------------------------------

# getCRUCLdata 0.1.5

## Major Changes

- `create_CRU_stack()` and `create_CRU_df()` now only work with locally
 available files. If you need to fetch and create a data frame or raster stack
 of the data, please use the new functions, `get_CRU_stack()` and
 `get_CRU_stack()`

- R >=3.2.0 now required

- Data can be cached using either `get_CRU_stack()` or `get_CRU_df()` for later
 use

## Minor Changes

- Improved documentation with examples on mapping and graphing and more detail
regarding the data itself

- Change the method in which files are downloaded to use `httr::GET()`

- Ingest data using `data.table::fread` to decrease the amount of time necessary
to run the functions

- Functions check to see if data file(s) have already been downloaded during
current R session, if so data file(s) are not requested for download again

- Months are returned as a factor object in the tidy data frame

--------------------------------------------------------------------------------

# getCRUCLdata 0.1.4

## Minor Changes

- Correct fix bug in data frame object generation where elevation was improperly
handled and function would stop

# getCRUCLdata 0.1.3

## Minor Changes

- Correct fix bug in raster object generation where the objects were incorrectly
cropped

- Update documentation with ROxygen 6.0.0

- Minor edits to documentation for clarity

--------------------------------------------------------------------------------

# getCRUCLdata 0.1.2

## Minor Changes

- Correct documentation to read that the data resolution is 10 minute, not 10
seconds

- Correct URLs in DESCRIPTION file

- Add required version for PURRR

- Add required version for R

- Corrected URL pointing to CRU readme.txt file

--------------------------------------------------------------------------------

# getCRUCLdata 0.1.1

## Minor Changes

- Renamed to getCRUdata as suggested by CRAN maintainers

- Revised description file as requested by CRAN maintainers

- Enhanced vignette

--------------------------------------------------------------------------------

## getCRUCL2.0 0.1.0

## Minor Changes

- Initial submission to CRAN
# getCRUCLdata v0.3.2

# Test environments

- local macOS install, R version 4.0.3 (2020-10-10)

- win-builder R Under development (unstable) (2020-10-23 r79366)

- win-builder, R version 4.0.3 (2020-10-10)

# R CMD check results

0 errors | 0 warnings | 1 note

# Resubmission as requested by CRAN maintainers

## Note

- This is a resubmission that fixes a link redirect as requested

# Reverse dependencies

No ERRORs or WARNINGs found.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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
# CONTRIBUTING #

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
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the getCRUCLdata project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
title: 'getCRUCLdata: Use and Explore CRU CL v. 2.0 Climatology Elements in R'
authors:
- affiliation: 1
  name: Adam H Sparks
  orcid: 0000-0002-0061-8359
output: pdf_document
tags:
- climate
- R
- applied climatology
- high resolution surface
- data
affiliations:
  index: 1
  name: University of Southern Queensland, Centre for Crop Health, Toowoomba Queensland 4350, Australia
bibliography: paper.bib
date: "04 April 2017"

---

# Summary

The CRU CL v. 2.0 data are a gridded climatology of 1961-1990 monthly means released in 2002 and cover all land areas (excluding Antarctica) at 10 arcminutes (0.1666667 degree) resolution [@New2002] providing precipitation, cv of precipitation, wet-days, mean temperature, mean diurnal temperature range, relative humidity, sunshine, ground-frost, windspeed and elevation. While these data have a high resolution and are freely available, the data format can be cumbersome for working with. Four functions are provided by _getCRUCLdata_ that automate importing these data into R [@R-base]. All of the functions facilitate the calculation of minimum temperature and maximum temperature, and format the data into a tidy data frame [@Wickham2014] in a _tibble_ [@Wickham2017] object or a list of _raster_ stack objects [@Raster] for use in R or easily exported to a raster format file for use in a geographic information system (GIS). Two functions, `get_CRU_df()` and `get_CRU_stack()` provide the ability to easily download CRU CL v. 2.0 data from the CRU website and import the data into R and allow for caching downloaded data. The other two functions, `create_CRU_df()` and `create_CRU_stack()` allow the user to easily import the data files from a local disk location and transform them into a tidy data frame _tibble_ or _raster_ stack. The data have applications in applied climatology, biogeochemical modelling, hydrology and agricultural meteorology [@New2002].

# References
## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.3 (2020-10-10) |
|os       |macOS Catalina 10.15.7       |
|system   |x86_64, darwin17.0           |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_AU.UTF-8                  |
|ctype    |en_AU.UTF-8                  |
|tz       |Australia/Brisbane           |
|date     |2020-10-25                   |

# Dependencies

|package      |old    |new        |Δ  |
|:------------|:------|:----------|:--|
|getCRUCLdata |0.3.1  |0.3.1.9000 |*  |
|assertthat   |0.2.1  |0.2.1      |   |
|cli          |2.1.0  |2.1.0      |   |
|crayon       |1.3.4  |1.3.4      |   |
|curl         |4.3    |4.3        |   |
|data.table   |1.13.2 |1.13.2     |   |
|digest       |0.6.27 |0.6.27     |   |
|ellipsis     |0.3.1  |0.3.1      |   |
|fansi        |0.4.1  |0.4.1      |   |
|glue         |1.4.2  |1.4.2      |   |
|hoardr       |0.5.2  |0.5.2      |   |
|lifecycle    |0.2.0  |0.2.0      |   |
|magrittr     |1.5    |1.5        |   |
|pillar       |1.4.6  |1.4.6      |   |
|pkgconfig    |2.0.3  |2.0.3      |   |
|R6           |2.4.1  |2.4.1      |   |
|rappdirs     |0.3.1  |0.3.1      |   |
|raster       |3.3-13 |3.3-13     |   |
|Rcpp         |1.0.5  |1.0.5      |   |
|rlang        |0.4.8  |0.4.8      |   |
|sp           |1.4-4  |1.4-4      |   |
|tibble       |3.0.4  |3.0.4      |   |
|utf8         |1.1.4  |1.1.4      |   |
|vctrs        |0.3.4  |0.3.4      |   |

# Revdeps

*Wow, no problems at all. :)*# Check times




*Wow, no problems at all. :)*---
title: "getCRUCLdata"
author: "Adam H. Sparks"
date: "2020-10-26"
output:
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteEngine{knitr::rmarkdown_notangle}
  %\VignetteIndexEntry{getCRUCLdata}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{raster}
  %\VignetteDepends{ggplot2}
  %\VignetteDepends{viridis}
---



## Introduction to _getCRUCLdata_

The _getCRUCLdata_ package provides functions that automate importing CRU CL v. 2.0 climatology data into R, facilitate the calculation of minimum temperature and maximum temperature, and formats the data into a
[tidy data frame](http://vita.had.co.nz/papers/tidy-data.html) as a `tibble::tibble()`
object or a [`list()`](https://www.rdocumentation.org/packages/base/versions/3.4.0/topics/list)
of [`raster::stack()`](https://www.rdocumentation.org/packages/raster/versions/2.5-8/topics/stack)
objects for use in an R session.

CRU CL v. 2.0 data are a gridded climatology of 1961-1990 monthly means released in 2002 and cover all land areas (excluding Antarctica) at 10 arcminutes (0.1666667 degree) resolution.
For more information see the description of the data provided by the University of East Anglia Climate Research Unit (CRU), <https://crudata.uea.ac.uk/cru/data/hrg/tmc/readme.txt>.

## Changes to original CRU CL v. 2.0 data

This package automatically converts elevation values from kilometres to metres.

This package crops all spatial outputs to an extent of ymin = -60, ymax = 85, xmin = -180, xmax = 180. Note that the original wind data include land area for parts of Antarctica.

# Using _getCRUCLdata_

Logical arguments are used to specify the climatology elements to retrieve and parse.
All arguments default to `FALSE`.
The `create_CRU_*()` functions require an additional parameter, `dsn` to be provided that states where the files are locally stored.
The arguments for selecting the climatology elements for importing are:

- **pre** Logical. Fetch precipitation (millimetres/month) from server and return in the data?

- **pre_cv** Logical. Fetch cv of precipitation (percent) from server and return in the data?

- **rd0** Logical. Fetch wet-days (number days with >0.1 millimetres rain per month) and return in the data?

- **dtr** Logical. Fetch mean diurnal temperature range (degrees Celsius) and return it in the data?

- **tmp** Logical. Fetch temperature (degrees Celsius) and return it in the data?

- **tmn** Logical. Calculate minimum temperature values (degrees Celsius) and return it in the data?

- **tmx** Logical. Calculate maximum temperature (degrees Celsius) and return it in the data?

- **reh** Logical. Fetch relative humidity and return it in the data?

- **sunp** Logical. Fetch sunshine, percent of maximum possible (percent of day length) and return it in data?

- **frs** Logical. Fetch ground-frost records (number of days with ground-frost per month) and return it in data?

- **wnd** Logical. Fetch 10m wind speed (metres/second) and return it in the data?

- **elv** Logical. Fetch elevation (and convert to metres from kilometres) and return it in the data?

- **dsn** *For `create_CRU_stack()`* and *`create_CRU_df()`* only.
Local file path where CRU CL v. 2.0 .dat.gz files are located.

### Creating tidy data frames for use in R

The `get_CRU_df()` function automates the download process and creates tidy data frames as a `tibble::tibble()` of the CRU CL v. 2.0 climatology elements.


```r
library(getCRUCLdata)

CRU_data <- get_CRU_df(pre = TRUE,
                       pre_cv = TRUE,
                       rd0 = TRUE,
                       tmp = TRUE,
                       dtr = TRUE,
                       reh = TRUE,
                       tmn = TRUE,
                       tmx = TRUE,
                       sunp = TRUE,
                       frs = TRUE,
                       wnd = TRUE,
                       elv = TRUE)

CRU_data
#> # A tibble: 6,795,150 x 15
#>      lat   lon month   dtr   frs   pre pre_cv   rd0   reh   sun   tmp   wnd   elv   tmx   tmn
#>    <dbl> <dbl> <fct> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1  30.9  35.4 <NA>   NA    NA    NA     NA    NA    NA    NA    NA    NA    -260 NA    NA   
#>  2  31.1  35.4 <NA>   NA    NA    NA     NA    NA    NA    NA    NA    NA    -361 NA    NA   
#>  3  31.2  35.4 <NA>   NA    NA    NA     NA    NA    NA    NA    NA    NA    -336 NA    NA   
#>  4  31.4  35.4 <NA>   NA    NA    NA     NA    NA    NA    NA    NA    NA    -284 NA    NA   
#>  5  31.8  35.6 <NA>   NA    NA    NA     NA    NA    NA    NA    NA    NA    -248 NA    NA   
#>  6  31.9  35.6 <NA>   NA    NA    NA     NA    NA    NA    NA    NA    NA    -210 NA    NA   
#>  7 -59.1 -26.6 jan     2.3  18.7 105.    35.2  17.1  88.6   9.4   0.2   6.4   193  1.35 -0.95
#>  8 -58.4 -26.4 jan     2.5  18.5 107.    36.2  17.2  88.5   9.9   0.4   6.4   239  1.65 -0.85
#>  9 -58.4 -26.2 jan     2.4  18.4 106.    36.2  17.1  88.5  10     0.6   6.4   194  1.8  -0.6 
#> 10 -55.9 -67.2 jan     7.6   8    73.1   44.1  13.3  80.7  34.3   8     5      64 11.8   4.2 
#> # … with 6,795,140 more rows
```

Perhaps you only need one or two elements, it is easy to create a tidy data frame of mean temperature only.


```r
t <- get_CRU_df(tmp = TRUE)

t
#> # A tibble: 6,795,144 x 4
#>      lat   lon month   tmp
#>    <dbl> <dbl> <fct> <dbl>
#>  1 -59.1 -26.6 jan     0.2
#>  2 -58.4 -26.2 jan     0.6
#>  3 -58.4 -26.4 jan     0.4
#>  4 -55.9 -67.2 jan     8  
#>  5 -55.8 -67.2 jan     8.2
#>  6 -55.8 -67.4 jan     8  
#>  7 -55.8 -67.6 jan     8.4
#>  8 -55.6 -67.4 jan     8.3
#>  9 -55.6 -67.6 jan     8.6
#> 10 -55.6 -68.1 jan     8.2
#> # … with 6,795,134 more rows
```

#### Plotting data from the tidy dataframe

Now that we have the data, we can plot it easily using _ggplot2_ and the _viridis_ package for the colour scale.


```r
library(ggplot2)
library(viridis)

ggplot(data = t, aes(x = lon, y = lat, fill = tmp)) +
  geom_tile() +
  scale_fill_viridis(option = "inferno") +
  coord_quickmap() +
  ggtitle("Global Mean Monthly Temperatures 1961-1990") +
  facet_wrap( ~ month, nrow = 4)
```

<img src="plot_t-1.png" title="plot of chunk plot_t" alt="plot of chunk plot_t" style="display: block; margin: auto;" />

We can also generate a violin plot of the same data to visualise how the temperatures change throughout the year.


```r
ggplot(data = t, aes(x = month, y = tmp)) +
  geom_violin() +
  ylab("Temperature (˚C)") +
  labs(title = "Global Monthly Mean Land Surface Temperatures From 1960-1991",
       subtitle = "Excludes Antarctica")
```

<img src="violin_plot-1.png" title="plot of chunk violin_plot" alt="plot of chunk violin_plot" style="display: block; margin: auto;" />

#### Saving the tidy `data.frame` as a CSV (comma separated values file) locally

Save the resulting tidy `data.frame` to local disk as a comma separated (CSV)
file to local disk, using _data.table_'s `fwrite()`.


```r
fwrite(x = t, file = "~/CRU_tmp.csv")
```

### Creating raster stacks for use in R and saving for use in another GIS

For working with spatial data, _getCRUCLdata_ provides a function that create lists of _raster_ stacks of the data.

The `get_CRU_stack()` functions provide similar functionality to `get_CRU_df()`, but rather than returning a tidy data frame, it returns a list of `raster::stack()` objects for use in an R session.

The `get_CRU_stack()` function automates the download process and creates a `raster::stack()` object of the CRU CL v. 2.0 climatology elements.
Illustrated here is creating a `raster::stack()` of all CRU CL v. 2.0 climatology elements available.


```r
CRU_stack <- get_CRU_stack(
  pre = TRUE,
  pre_cv = TRUE,
  rd0 = TRUE,
  tmp = TRUE,
  dtr = TRUE,
  reh = TRUE,
  tmn = TRUE,
  tmx = TRUE,
  sunp = TRUE,
  frs = TRUE,
  wnd = TRUE,
  elv = TRUE
)

CRU_stack
#> $dtr
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      :  jan,  feb,  mar,  apr,  may,  jun,  jul,  aug,  sep,  oct,  nov,  dec 
#> min values :  2.3,  2.1,  2.2,  2.3,  1.8,  2.5,  2.8,  2.4,  2.2,  2.8,  2.6,  2.0 
#> max values : 22.7, 23.1, 23.5, 24.0, 24.0, 25.2, 25.8, 25.6, 25.5, 22.6, 22.9, 21.9 
#> 
#> 
#> $elv
#> class      : RasterLayer 
#> dimensions : 870, 2160, 1879200  (nrow, ncol, ncell)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      : elv 
#> values     : -361, 6486  (min, max)
#> 
#> 
#> $frs
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      :  jan,  feb,  mar,  apr,  may,  jun,  jul,  aug,  sep,  oct,  nov,  dec 
#> min values :    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 
#> max values : 31.0, 28.3, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0 
#> 
#> 
#> $pre
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 24  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      :    jan,    feb,    mar,    apr,    may,    jun,    jul,    aug,    sep,    oct,    nov,    dec, pre_cv_jan, pre_cv_feb, pre_cv_mar, ... 
#> min values :    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,        0.0,      -10.5,        5.8, ... 
#> max values :  910.1,  824.3,  727.3,  741.3, 1100.0, 2512.6, 2505.5, 1799.4,  849.8,  851.6,  843.7,  733.3,      496.2,      495.5,      482.0, ... 
#> 
#> 
#> $rd0
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      :  jan,  feb,  mar,  apr,  may,  jun,  jul,  aug,  sep,  oct,  nov,  dec 
#> min values :    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 
#> max values : 31.0, 28.2, 31.0, 30.0, 30.7, 30.0, 31.0, 31.0, 29.1, 28.4, 28.5, 30.3 
#> 
#> 
#> $reh
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      :   jan,   feb,   mar,   apr,   may,   jun,   jul,   aug,   sep,   oct,   nov,   dec 
#> min values :  18.4,  14.6,  13.5,  13.4,  15.5,  10.2,  10.8,  10.1,  11.0,  14.2,  19.0,  19.7 
#> max values : 100.0, 100.0, 100.0, 100.0,  96.9,  95.1,  96.9,  97.1,  95.5, 100.0, 100.0, 100.0 
#> 
#> 
#> $sun
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      :  jan,  feb,  mar,  apr,  may,  jun,  jul,  aug,  sep,  oct,  nov,  dec 
#> min values :  0.0,  0.0,  3.3,  4.3,  8.1,  6.6,  5.3,  8.4,  4.5,  0.8,  0.0,  0.0 
#> max values : 92.8, 93.0, 90.2, 93.1, 94.0, 98.9, 98.8, 98.8, 99.1, 95.8, 94.6, 93.1 
#> 
#> 
#> $tmp
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      :   jan,   feb,   mar,   apr,   may,   jun,   jul,   aug,   sep,   oct,   nov,   dec 
#> min values : -51.6, -47.6, -45.2, -36.6, -22.2, -16.3, -14.0, -17.3, -26.4, -36.3, -41.6, -49.0 
#> max values :  32.5,  32.1,  32.4,  34.3,  36.0,  38.3,  37.8,  36.8,  34.8,  32.8,  32.4,  32.2 
#> 
#> 
#> $wnd
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : memory
#> names      : jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec 
#> min values : 0.1, 0.1, 0.3, 0.4, 0.3, 0.2, 0.3, 0.4, 0.5, 0.4, 0.2, 0.2 
#> max values : 9.8, 9.6, 9.4, 9.0, 8.7, 8.6, 9.1, 9.3, 9.3, 9.7, 9.6, 9.4 
#> 
#> 
#> $tmn
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : /private/var/folders/1q/045gnpqd7dnfgzshmdn8hshw0000gp/T/Rtmp4VtMtp/raster/r_tmp_2020-10-26_111004_18359_15367.grd 
#> names      : layer.1, layer.2, layer.3, layer.4, layer.5, layer.6, layer.7, layer.8, layer.9, layer.10, layer.11, layer.12 
#> min values :  -55.05,  -52.95,  -48.75,  -41.35,  -28.00,  -21.40,  -18.75,  -22.55,  -31.45,   -40.60,   -45.75,   -52.50 
#> max values :   26.30,   26.25,   27.40,   27.50,   30.00,   30.65,   30.60,   30.40,   28.75,    26.95,    25.90,    26.55 
#> 
#> 
#> $tmx
#> class      : RasterBrick 
#> dimensions : 870, 2160, 1879200, 12  (nrow, ncol, ncell, nlayers)
#> resolution : 0.1666667, 0.1666667  (x, y)
#> extent     : -180, 180, -60, 85  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> source     : /private/var/folders/1q/045gnpqd7dnfgzshmdn8hshw0000gp/T/Rtmp4VtMtp/raster/r_tmp_2020-10-26_111006_18359_67835.grd 
#> names      : layer.1, layer.2, layer.3, layer.4, layer.5, layer.6, layer.7, layer.8, layer.9, layer.10, layer.11, layer.12 
#> min values :  -48.20,  -43.35,  -41.65,  -32.45,  -17.55,  -11.50,  -10.85,  -12.30,  -21.65,   -32.05,   -37.55,   -45.50 
#> max values :   39.70,   38.40,   40.25,   41.85,   43.60,   45.95,   45.70,   44.85,   42.35,    39.50,    39.20,    39.90
```

The `create_CRU_stack()` function works in the same way with only one minor difference.
You must supply the location of the files on the local disk (`dsn`) that you wish to import.


```r
t <- create_CRU_stack(tmp = TRUE, dsn = "~/Downloads")
```

#### Plotting raster stacks of tmin and tmax

Because the stacks are in a `list()`, we need to access each element of the list individually to plot them, that's what the `[[1]]` or `[[2]]` is, the first or second element of the list.
Here using `[[7]]` we will plot the monthly average minimum temperature for all twelve months.


```r
library(raster)

plot(CRU_stack[[7]])
```

To plot only one month from the stack is also possible. Here we plot maximum temperature for July.
Note that we use indexing `[[2]]` as before but append a `$jul` to the object.
This is the name of the layer in the `raster::stack()`.
So, we are telling R to plot the second object in the `CRU_stack` list, which is `tmx` and from that raster stack, plot only the layer for July.


```r
plot(t[[8]]$jul)
```

#### Saving raster objects to local disk

The raster stack objects can be saved to disk as geotiff files (others are available, see help for `raster::writeRaster()` and `raster::writeFormats()` for more options) on the `Data` directory with a tmn or tmx prefix to the month for a file name.


```r
library(raster)

dir.create(file.path("~/Data"), showWarnings = FALSE)
writeRaster(
  t$tmn,
  filename = file.path("~/Data/tmn_", names(t$tmn)),
  bylayer = TRUE,
  format = "GTiff"
)

writeRaster(
  t$tmx,
  filename = file.path("~/Data/tmx_", names(t$tmn)),
  bylayer = TRUE,
  format = "GTiff"
)
```

# CRU CL v. 2.0 reference and abstract

Mark New (1,*), David Lister (2), Mike Hulme (3), Ian Makin (4)
A high-resolution data set of surface climate over global land areas Climate Research, 2000, Vol 21, pg 1-25
(1) School of Geography and the Environment, University of Oxford,
    Mansfield Road, Oxford OX1 3TB, United Kingdom
(2) Climatic Research Unit, and (3) Tyndall Centre for Climate Change Research,
    both at School of Environmental Sciences, University of East Anglia,
    Norwich NR4 7TJ, United Kingdom
(4) International Water Management Institute, PO Box 2075, Colombo, Sri Lanka

**ABSTRACT:** We describe the construction of a 10-minute latitude/longitude
data set of mean monthly surface climate over global land areas, excluding
Antarctica. The climatology includes 8 climate elements - precipitation, wet-day
frequency, temperature, diurnal temperature range, relative humidity,sunshine
duration, ground frost frequency and windspeed - and was interpolated from a
data set of station means for the period centred on 1961 to 1990. Precipitation
was first defined in terms of the parameters of the Gamma distribution, enabling
the calculation of monthly precipitation at any given return period. The data
are compared to an earlier data set at 0.5 degrees latitude/longitude resolution
and show added value over most regions. The data will have many applications in
applied climatology, biogeochemical modelling, hydrology and agricultural
meteorology and are available through the School of Geography Oxford
(http://www.geog.ox.ac.uk), the International Water Management Institute
"World Water and Climate Atlas" (https://www.iwmi.cgiar.org/) and the Climatic
Research Unit (http://www.cru.uea.ac.uk).
---
title: "Advanced usage of getCRUCLdata"
author: "Adam H. Sparks"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Advanced usage of getCRUCLdata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Caching files for later use

When using the `get_CRU_df()` or `get_CRU_stack()` functions, files may be cached in the users' local space for later use (optional) or stored in a temporary directory and deleted when the R session is closed and not saved (this is the default behaviour already illustrated above).
Illustrated here, create a tidy data frame of all CRU CL v. 2.0 climatology elements available and cache
them to save time in the future.
*In order to take advantage of the cached data, you must use the `get_CRU_df()` function again in the future*.
This functionality is somewhat modelled after the `raster::getData()` function that will not download files that already exist in the working directory, however in this case the function is portable and it will work for any working directory.
That is, if you have cached the data and you use `get_CRU_df()` again, it will use the cached data no matter what working directory you are in.
This functionality will be most useful for writing scripts that may be used several times rather than just once off or if you frequently use the data in multiple analyses the data will not be downloaded again if they have been cached.

Create a list of raster stacks of maximum and minimum temperature.
To take advantage of the previously cached files and save time by not downloading files, specify `cache = TRUE`.

```{r, eval=FALSE}
tmn_tmx <- get_CRU_stack(tmn = TRUE,
                         tmx = TRUE,
                         cache = TRUE)
```

## Handling files downloaded outside of R

A second set of functions, `create_CRU_df()` and `create_CRU_stack()`, is provided for users that may have connectivity issues or simply wish to use something other than R to download the data files.
You may also wish to use these if you want to download the data and specify where it is stored rather than using the `cache` functionality of `get_CRU_df()` and `get_CRU_stack()`.

The `create_CRU_df()` and `create_CRU_stack()` functions work in the same way as `get_CRU_df()` and `get_CRU_stack()` functions with only one major difference.
You must supply the location of the files on the local disk (`dsn`) that you wish to import.
That is, the CRU CL v. 2.0 data files *must* be downloaded prior to the use of these functions using a program external to R.

```{r, eval=FALSE}
t <- create_CRU_df(tmp = TRUE, dsn = "~/Downloads")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCRUCLdata-package.R
\docType{package}
\name{getCRUCLdata-package}
\alias{getCRUCLdata}
\alias{getCRUCLdata-package}
\title{getCRUCLdata: 'CRU' 'CL' v. 2.0 Climatology Client}
\description{
Provides functions that automate downloading and importing University of East Anglia Climate Research Unit ('CRU') 'CL' v. 2.0 climatology data, facilitates the calculation of minimum temperature and maximum temperature and formats the data into a data frame or a list of 'raster' 'stack' objects for use. 'CRU' 'CL' v. 2.0 data are a gridded climatology of 1961-1990 monthly means released in 2002 and cover all land areas (excluding Antarctica) at 10 arc minutes (0.1666667 degree) resolution. For more information see the description of the data provided by the University of East Anglia Climate Research Unit, <https://crudata.uea.ac.uk/cru/data/hrg/tmc/readme.txt>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/getCRUCLdata/}
  \item Report bugs at \url{https://github.com/ropensci/getCRUCLdata/issues}
}

}
\author{
\strong{Maintainer}: Adam H. Sparks \email{adamhsparks@gmail.com} (\href{https://orcid.org/0000-0002-0061-8359}{ORCID})

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_CRU_stack.R
\name{create_CRU_stack}
\alias{create_CRU_stack}
\title{Create a list of raster stack objects from local disk files}
\usage{
create_CRU_stack(
  pre = FALSE,
  pre_cv = FALSE,
  rd0 = FALSE,
  tmp = FALSE,
  dtr = FALSE,
  reh = FALSE,
  tmn = FALSE,
  tmx = FALSE,
  sunp = FALSE,
  frs = FALSE,
  wnd = FALSE,
  elv = FALSE,
  dsn = ""
)
}
\arguments{
\item{pre}{Logical. Fetch precipitation (millimetres/month) from server and
return in the data frame? Defaults to \code{FALSE}.}

\item{pre_cv}{Logical. Fetch cv of precipitation (percent) from server and
return in the data frame? Defaults to \code{FALSE}. NOTE. Setting this to
\code{TRUE} will always results in \strong{pre} being set to \code{TRUE} and
returned as well.}

\item{rd0}{Logical. Fetch wet-days (number days with >0.1millimetres rain per
month) and return in the data frame? Defaults to \code{FALSE}.}

\item{tmp}{Logical. Fetch temperature (degrees Celsius) and return it in the
data frame? Defaults to \code{FALSE}.}

\item{dtr}{Logical. Fetch mean diurnal temperature range (degrees Celsius)
and return it in the data frame? Defaults to \code{FALSE}.}

\item{reh}{Logical. Fetch relative humidity and return it in the data frame?
Defaults to \code{FALSE}.}

\item{tmn}{Logical. Calculate minimum temperature values (degrees Celsius)
and return it in the data frame? Defaults to \code{FALSE}.}

\item{tmx}{Logical. Calculate maximum temperature (degrees Celsius) and
return it in the data frame? Defaults to \code{FALSE}.}

\item{sunp}{Logical. Fetch sunshine, percent of maximum possible (percent of
day length) and return it in data frame? Defaults to \code{FALSE}.}

\item{frs}{Logical. Fetch ground-frost records (number of days with ground-
frost per month) and return it in data frame? Defaults to \code{FALSE}.}

\item{wnd}{Logical. Fetch 10m wind speed (metres/second) and return it in the
data frame? Defaults to \code{FALSE}.}

\item{elv}{Logical. Fetch elevation (converted to metres) and return it in
the data frame? Defaults to \code{FALSE}.}

\item{dsn}{Local file path where \acronym{CRU} \acronym{CL} v.2.0 .dat.gz
files are located.}
}
\value{
A \code{\link[base]{list}} of \code{\link{raster}}
\code{\link[raster]{stack}} objects of \acronym{CRU} \acronym{CL} v. 2.0
climatology elements
}
\description{
Automates importing \acronym{CRU} \acronym{CL} v.2.0 climatology
data and creates a \code{\link[raster]{stack}} of the data.  If requested,
minimum and maximum temperature may also be automatically calculated as
described in the data readme.txt file.  This function can be useful if you
have network connection issues that mean automated downloading of the files
using \R does not work properly.

Nomenclature and units from readme.txt:
\describe{
\item{pre}{precipitation (millimetres/month)}
  \describe{
    \item{cv}{cv of precipitation (percent)}
  }
\item{rd0}{wet-days (number days with >0.1mm rain per month)}
\item{tmp}{mean temperature (degrees Celsius)}
\item{dtr}{mean diurnal temperature range (degrees Celsius)}
\item{reh}{relative humidity (percent)}
\item{sunp}{sunshine (percent of maximum possible (percent of day length))}
\item{frs}{ground-frost (number of days with ground-frost per month)}
\item{wnd}{10 metre windspeed (metres/second)}
\item{elv}{elevation (automatically converted to metres)}
}
For more information see the description of the data provided by
\acronym{CRU}, \url{https://crudata.uea.ac.uk/cru/data/hrg/tmc/readme.txt}
}
\note{
This package automatically converts elevation values from kilometres to
metres.

This package crops all spatial outputs to an extent of ymin = -60, ymax = 85,
xmin = -180, xmax = 180. Note that the original wind data include land area
for parts of Antarctica, these data are excluded in the raster stacks
generated by this function.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Create a raster stack of temperature from tmp
# files in the tempdir() directory.

download.file(
  url = "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_tmp.dat.gz",
  destfile = file.path(tempdir(), "grid_10min_tmp.dat.gz")
)

CRU_tmp <- create_CRU_stack(tmp = TRUE, dsn = tempdir())

CRU_tmp
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link{get_CRU_stack}}
}
\author{
Adam H. Sparks, \email{adamhsparks@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manage_cached_files.R
\name{manage_cache}
\alias{manage_cache}
\title{Manage locally cached CRU CL v. 2.0 files}
\description{
Manage cached \CRANpkg{getCRUCLdata} files with \CRANpkg{hoardr}
}
\details{
The default cache directory is
\code{file.path(rappdirs::user_cache_dir(), "R/getCRUCLdata")}, but you can
set your own path using \code{manage_cache$cache_path_set()}

\code{manage_cache$cache_delete} only accepts one file name, while
\code{manage_cache$cache_delete_all}
does not accept any names, but deletes all files. For deleting many specific
files, use \code{manage_cache$cache_delete} in an
\code{\link[base]{lapply}()} type call.
}
\section{Useful user functions}{

\itemize{
 \item \code{manage_cache$cache_path_get()} - get cache path
 \item \code{manage_cache$cache_path_set()} - set cache path
 \item \code{manage_cache$list()} - returns a character vector of full
 path file names
 \item \code{manage_cache$files()} - returns file objects with metadata
 \item \code{manage_cache$details()} - returns files with details
 \item \code{manage_cache$delete()} - delete specific files
 \item \code{manage_cache$delete_all()} - delete all files, returns
 nothing
}
}

\examples{
\dontrun{

# list files in cache
manage_cache$list()

# delete certain database files
manage_cache$delete("file path")
manage_cache$list()

# delete all files in cache
manage_cache$delete_all()
manage_cache$list()

# set a different cache path from the default
manage_cache$cache_path_set("~/tmp")
}
}
\author{
Adam H. Sparks, \email{adamhsparks@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_CRU_df.R
\name{get_CRU_df}
\alias{get_CRU_df}
\title{Download and create a data frame of climatology parameters}
\usage{
get_CRU_df(
  pre = FALSE,
  pre_cv = FALSE,
  rd0 = FALSE,
  tmp = FALSE,
  dtr = FALSE,
  reh = FALSE,
  tmn = FALSE,
  tmx = FALSE,
  sunp = FALSE,
  frs = FALSE,
  wnd = FALSE,
  elv = FALSE,
  cache = FALSE
)
}
\arguments{
\item{pre}{Logical.  Fetch precipitation (millimetres/month) from server and
return in the data frame?  Defaults to \code{FALSE}.}

\item{pre_cv}{Logical.  Fetch cv of precipitation (percent) from server and
return in the data frame?  Defaults to \code{FALSE}.  NOTE.  Setting this to
\code{TRUE} will always results in \strong{pre} being set to \code{TRUE} and
returned as well.}

\item{rd0}{Logical.  Fetch wet-days (number days with >0.1 millimetres rain
per month) and return in the data frame?  Defaults to \code{FALSE}.}

\item{tmp}{Logical.  Fetch temperature (degrees Celsius) and return it in the
data frame?  Defaults to \code{FALSE}.}

\item{dtr}{Logical.  Fetch mean diurnal temperature range (degrees Celsius)
and return it in the data frame?  Defaults to \code{FALSE}.}

\item{reh}{Logical.  Fetch relative humidity and return it in the data frame?
Defaults to FALSE.}

\item{tmn}{Logical.  Calculate minimum temperature values (degrees Celsius)
and return it in the data frame?  Defaults to \code{FALSE}.}

\item{tmx}{Logical.  Calculate maximum temperature (degrees Celsius) and
return it in the data frame?  Defaults to \code{FALSE}.}

\item{sunp}{Logical.  Fetch sunshine, percent of maximum possible (percent of
day length) and return it in data frame?  Defaults to \code{FALSE}.}

\item{frs}{Logical.  Fetch ground-frost records (number of days with ground-
frost per month) and return it in data frame?  Defaults to \code{FALSE}.}

\item{wnd}{Logical.  Fetch 10m wind speed (metres/second) and return it in the
data frame? Defaults to \code{FALSE}.}

\item{elv}{Logical.  Fetch elevation (converted to metres) and return it in
the data frame?  Defaults to \code{FALSE}.}

\item{cache}{Logical.  Store CRU CL v. 2.0 data files locally for later use?
If \code{FALSE}, the downloaded files are removed when R session is closed.
To take advantage of cached files in future sessions, use \code{cache = TRUE}
after the initial download and caching.  Defaults to \code{FALSE}.}
}
\value{
A tidy data frame of \acronym{CRU} \acronym{CL} v. 2.0 climatology
elements as a \code{\link[tibble]{tibble}} object
}
\description{
This function automates downloading and importing \acronym{CRU}
\acronym{CL} v. 2.0 climatology data and creates a data frame of the data.
If requested, minimum and maximum temperature may also be automatically
calculated as described in the data readme.txt file.  Data may be cached for
later use by this function, saving time downloading files in future use of
the function.

Nomenclature and units from readme.txt:
\describe{
\item{pre}{precipitation (millimetres/month)}
  \describe{
   \item{cv}{cv of precipitation (percent)}
  }
\item{rd0}{wet-days (number days with >0.1mm rain per month)}
\item{tmp}{mean temperature (degrees Celsius)}
\item{dtr}{mean diurnal temperature range (degrees Celsius)}
\item{reh}{relative humidity (percent)}
\item{sunp}{sunshine (percent of maximum possible (percent of day length))}
\item{frs}{ground-frost (number of days with ground-frost per month)}
\item{wnd}{10 metre windspeed (metres/second)}
\item{elv}{elevation (automatically converted to metres)}
}
For more information see the description of the data provided by
\acronym{CRU}, \url{https://crudata.uea.ac.uk/cru/data/hrg/tmc/readme.txt}
}
\note{
This package automatically converts elevation values from kilometres to
metres.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Download data and create a data frame of precipitation and temperature
# without caching the data files
CRU_pre_tmp <- get_CRU_df(pre = TRUE, tmp = TRUE)

head(CRU_pre_tmp)
CRU_pre_tmp
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link{create_CRU_stack}}
\code{\link{manage_cache}}
}
\author{
Adam H. Sparks, \email{adamhsparks@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_CRU_df.R
\name{create_CRU_df}
\alias{create_CRU_df}
\title{Create a data of climatology variables from local disk files}
\usage{
create_CRU_df(
  pre = FALSE,
  pre_cv = FALSE,
  rd0 = FALSE,
  tmp = FALSE,
  dtr = FALSE,
  reh = FALSE,
  tmn = FALSE,
  tmx = FALSE,
  sunp = FALSE,
  frs = FALSE,
  wnd = FALSE,
  elv = FALSE,
  dsn = ""
)
}
\arguments{
\item{pre}{Logical. Fetch precipitation (millimetres/month) from server and
return in the data frame? Defaults to \code{FALSE}.}

\item{pre_cv}{Logical. Fetch cv of precipitation (percent) from server and
return in the data frame? Defaults to \code{FALSE}. NOTE. Setting this to
\code{TRUE} will always results in \strong{pre} being set to \code{TRUE} and
returned as well.}

\item{rd0}{Logical. Fetch wet-days (number days with >0.1millimetres rain per
month) and return in the data frame? Defaults to \code{FALSE}.}

\item{tmp}{Logical. Fetch temperature (degrees Celsius) and return it in the
data frame? Defaults to \code{FALSE}.}

\item{dtr}{Logical. Fetch mean diurnal temperature range (degrees Celsius)
and return it in the data frame? Defaults to \code{FALSE}.}

\item{reh}{Logical. Fetch relative humidity and return it in the data frame?
Defaults to \code{FALSE}.}

\item{tmn}{Logical. Calculate minimum temperature values (degrees Celsius)
and return it in the data frame? Defaults to \code{FALSE}.}

\item{tmx}{Logical. Calculate maximum temperature (degrees Celsius) and
return it in the data frame? Defaults to \code{FALSE}.}

\item{sunp}{Logical. Fetch sunshine, percent of maximum possible (percent of
day length) and return it in data frame? Defaults to \code{FALSE}.}

\item{frs}{Logical. Fetch ground-frost records (number of days with ground-
frost per month) and return it in data frame? Defaults to \code{FALSE}.}

\item{wnd}{Logical. Fetch 10m wind speed (metres/second) and return it in the
data frame? Defaults to \code{FALSE}.}

\item{elv}{Logical. Fetch elevation (converted to metres) and return it in
the data frame? Defaults to \code{FALSE}.}

\item{dsn}{Local file path where \acronym{CRU} \acronym{CL} v.2.0 .dat.gz
files are located.}
}
\value{
A tidy data frame of \acronym{CRU} \acronym{CL} v. 2.0 climatology
elements as a \code{\link[base]{data.frame}} object
}
\description{
Automates importing \acronym{CRU} \acronym{CL} v.2.0 climatology
data and creates a tidy data frame of the data.  If requested, minimum and
maximum temperature may also be automatically calculated as described in the
data readme.txt file.  This function can be useful if you have network
connection issues that mean automated downloading of the files using \R
does not work properly.

Nomenclature and units from readme.txt:
\describe{
\item{pre}{precipitation (millimetres/month)}
 \describe{
   \item{cv}{cv of precipitation (percent)}
 }
\item{rd0}{wet-days (number days with >0.1mm rain per month)}
\item{tmp}{mean temperature (degrees Celsius)}
\item{dtr}{mean diurnal temperature range (degrees Celsius)}
\item{reh}{relative humidity (percent)}
\item{sunp}{sunshine (percent of maximum possible (percent of day length))}
\item{frs}{ground-frost (number of days with ground-frost per month)}
\item{wnd}{10 metre windspeed (metres/second)}
\item{elv}{elevation (automatically converted to metres)}
}
For more information see the description of the data provided by
\acronym{CRU}, \url{https://crudata.uea.ac.uk/cru/data/hrg/tmc/readme.txt}
}
\note{
This package automatically converts elevation values from kilometres to
metres.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Create a data frame of temperature from locally available files in the
# tempdir() directory.

download.file(
  url = "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_tmp.dat.gz",
  destfile = file.path(tempdir(), "grid_10min_tmp.dat.gz")
)

CRU_tmp <- create_CRU_df(tmp = TRUE, dsn = tempdir())

CRU_tmp
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link{get_CRU_df}}
}
\author{
Adam H Sparks, \email{adamhsparks@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_CRU_stack.R
\name{get_CRU_stack}
\alias{get_CRU_stack}
\title{Download and create a list of raster stacks of climatology parameters}
\usage{
get_CRU_stack(
  pre = FALSE,
  pre_cv = FALSE,
  rd0 = FALSE,
  tmp = FALSE,
  dtr = FALSE,
  reh = FALSE,
  tmn = FALSE,
  tmx = FALSE,
  sunp = FALSE,
  frs = FALSE,
  wnd = FALSE,
  elv = FALSE,
  cache = FALSE
)
}
\arguments{
\item{pre}{Logical.  Fetch precipitation (millimetres/month) from server and
return in the data frame?  Defaults to \code{FALSE}.}

\item{pre_cv}{Logical.  Fetch cv of precipitation (percent) from server and
return in the data frame?  Defaults to \code{FALSE}.  NOTE.  Setting this to
\code{TRUE} will always results in \strong{pre} being set to \code{TRUE} and
returned as well.}

\item{rd0}{Logical.  Fetch wet-days (number days with >0.1millimetres rain
per month) and return in the data frame? Defaults to \code{FALSE}.}

\item{tmp}{Logical.  Fetch temperature (degrees Celsius) and return it in the
data frame?  Defaults to \code{FALSE}.}

\item{dtr}{Logical.  Fetch mean diurnal temperature range (degrees Celsius)
and return it in the data frame?  Defaults to \code{FALSE}.}

\item{reh}{Logical.  Fetch relative humidity and return it in the data frame?
Defaults to FALSE.}

\item{tmn}{Logical.  Calculate minimum temperature values (degrees Celsius)
and return it in the data frame?  Defaults to \code{FALSE}.}

\item{tmx}{Logical.  Calculate maximum temperature (degrees Celsius) and
return it in the data frame?  Defaults to \code{FALSE}.}

\item{sunp}{Logical.  Fetch sunshine, percent of maximum possible (percent of
day length) and return it in data frame?  Defaults to \code{FALSE}.}

\item{frs}{Logical. Fetch ground-frost records (number of days with ground-
frost per month) and return it in data frame?  Defaults to \code{FALSE}.}

\item{wnd}{Logical.  Fetch 10m wind speed (metres/second) and return it in
the data frame? Defaults to \code{FALSE}.}

\item{elv}{Logical.  Fetch elevation (converted to metres) and return it in
the data frame?  Defaults to \code{FALSE}.}

\item{cache}{Logical.  Store CRU CL v. 2.0 data files locally for later use?
If \code{FALSE}, the downloaded files are removed when R session is closed.
To take advantage of cached files in future sessions, use \code{cache = TRUE}
after the initial download and caching.  Defaults to \code{FALSE}.}
}
\value{
A \code{\link[base]{list}} of \code{\link{raster}}
\code{\link[raster]{stack}} objects of CRU CL v. 2.0 climatology elements
}
\description{
This function automates downloading and importing CRU CL v. 2.0
climatology data into \R and creates a list of raster stacks of the
data.  If requested, minimum and maximum temperature may also be
automatically calculated as described in the data readme.txt file.  Data may
be cached for later use by this function, saving time downloading files in
future use of the function.

Nomenclature and units from readme.txt:
\describe{
\item{pre}{precipitation (millimetres/month)}
  \describe{
    \item{cv}{cv of precipitation (percent)}
  }
\item{rd0}{wet-days (number days with >0.1mm rain per month)}
\item{tmp}{mean temperature (degrees Celsius)}
\item{dtr}{mean diurnal temperature range (degrees Celsius)}
\item{reh}{relative humidity (percent)}
\item{sunp}{sunshine (percent of maximum possible (percent of day length))}
\item{frs}{ground-frost (number of days with ground-frost per month)}
\item{wnd}{10 metre windspeed (metres/second)}
\item{elv}{elevation (automatically converted to metres)}
}
For more information see the description of the data provided by CRU,
\url{https://crudata.uea.ac.uk/cru/data/hrg/tmc/readme.txt}
}
\note{
This package automatically converts elevation values from kilometres to
metres.

This package crops all spatial outputs to an extent of ymin = -60, ymax = 85,
xmin = -180, xmax = 180.  Note that the original wind data include land area
for parts of Antarctica, these data are excluded in the raster stacks
generated by this function.
}
\examples{
\donttest{
# Download data and create a raster stack of precipitation and temperature
# without caching the data files
CRU_pre_tmp <- get_CRU_stack(pre = TRUE, tmp = TRUE)

CRU_pre_tmp
}

}
\seealso{
\code{\link{create_CRU_stack}}
\code{\link{manage_cache}}
}
\author{
Adam H. Sparks, \email{adamhsparks@gmail.com}
}
