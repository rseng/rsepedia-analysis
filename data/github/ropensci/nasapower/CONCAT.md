
# *nasapower*: NASA POWER API Client <img align='right' src='man/figures/logo.png'>

<!-- badges: start -->

[![tic](https://github.com/ropensci/nasapower/workflows/tic/badge.svg?branch=main)](https://github.com/ropensci/nasapower/actions)
[![codecov](https://codecov.io/gh/ropensci/nasapower/branch/master/graph/badge.svg?token=Kq9aea0TQN)](https://codecov.io/gh/ropensci/nasapower)
[![DOI](https://zenodo.org/badge/109224461.svg)](https://zenodo.org/badge/latestdoi/109224461)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![peer-review](https://badges.ropensci.org/155_status.svg)](https://github.com/ropensci/software-review/issues/155)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01035/status.svg)](https://doi.org/10.21105/joss.01035)
[![CRAN](http://www.r-pkg.org/badges/version/nasapower)](https://CRAN.R-project.org/package=nasapower)
<!-- badges: end -->

*nasapower* aims to make it quick and easy to automate downloading
[NASA-POWER](https://power.larc.nasa.gov) global meteorology, surface
solar energy and climatology data in your R session as a tidy data frame
`tibble` object for analysis and use in modelling or other purposes.
POWER (Prediction Of Worldwide Energy Resource) data are freely
available for download with varying spatial resolutions dependent on the
original data and with several temporal resolutions depending on the
POWER parameter and community.

**Note that the data are not static and may be replaced with improved
data.** Please see <https://power.larc.nasa.gov/docs/services/> for
detailed information in this regard.

### Quick start

*nasapower* can easily be installed using the following code.

#### From CRAN

The stable version is available through CRAN.

``` r
install.packages("nasapower")
```

#### From GitHub for the version in-development

A development version is available through GitHub.

``` r
if (!require("remotes")) {
  install.packages("remotes")
}

remotes::install_github("ropensci/nasapower")
```

### Example

Fetch daily “ag” community temperature, relative humidity and
precipitation for January 1, 1985 for Kingsthorpe, Queensland,
Australia.

``` r
library("nasapower")
daily_ag <- get_power(community = "ag",
                      lonlat = c(151.81, -27.48),
                      pars = c("RH2M", "T2M", "PRECTOTCORR"),
                      dates = "1985-01-01",
                      temporal_api = "daily"
                      )
daily_ag
```

    ## NASA/POWER CERES/MERRA2 Native Resolution Daily Data  
    ##  Dates (month/day/year): 01/01/1985 through 01/01/1985  
    ##  Location: Latitude  -27.48   Longitude 151.81  
    ##  Elevation from MERRA-2: Average for 0.5 x 0.625 degree lat/lon region = 442.77 meters 
    ##  Value for missing model data cannot be computed or out of model availability range: NA  
    ##  Parameter(s):  
    ##  
    ##  Parameters: 
    ##  RH2M            MERRA-2 Relative Humidity at 2 Meters (%) ;
    ##  T2M             MERRA-2 Temperature at 2 Meters (C) ;
    ##  PRECTOTCORR     MERRA-2 Precipitation Corrected (mm/day)  
    ##  
    ## # A tibble: 1 × 10
    ##     LON   LAT  YEAR    MM    DD   DOY YYYYMMDD    RH2M   T2M PRECTOTCORR
    ##   <dbl> <dbl> <dbl> <int> <int> <int> <date>     <dbl> <dbl>       <dbl>
    ## 1  152. -27.5  1985     1     1     1 1985-01-01  54.7  24.9         0.9

## Documentation

More documentation is available in the vignette in your R session,
`vignette("nasapower")` or available online,
<https://docs.ropensci.org/nasapower/>.

## Use of POWER Data

While *nasapower* does not redistribute the data or provide it in any
way, we encourage users to follow the requests of the POWER Project
Team.

> When POWER data products are used in a publication, we request the
> following acknowledgement be included: “These data were obtained from
> the NASA Langley Research Center POWER Project funded through the NASA
> Earth Science Directorate Applied Science Program.”

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/nasapower/issues).

-   License: MIT

-   To cite *nasapower*, please use the output from
    `citation(package = "nasapower")`.

-   Please note that the *nasapower* project is released with a
    [Contributor Code of
    Conduct](https://github.com/ropensci/nasapower/blob/main/CODE_OF_CONDUCT.md).
    By participating in the *nasapower* project you agree to abide by
    its terms.

-   The U.S. Earth System Research Laboratory, Physical Science Division
    of the National Atmospheric & Oceanic Administration (NOAA)
    maintains a list of gridded climate data sets that provide different
    data and different resolutions <https://psl.noaa.gov/data/gridded/>.

## References

<https://power.larc.nasa.gov>

<https://power.larc.nasa.gov/docs/methodology/>
# nasapower 4.0.3

## Minor changes

* Fixes tests that should use _vcr_ or be skipped on CRAN.

# nasapower 4.0.2

## Minor changes

* Update checks for number of parameters requested by user, maximum of 15 for hourly and 20 for all other temporal APIs.

* Return API messages to user to assist with troubleshooting when an error occurs server-side, see [Issue 55](https://github.com/ropensci/nasapower/issues/55).

* The list of POWER parameters that can be queried from the API, `parameters`, is now in alphabetical order.

* Add paragraph to vignette describing how to work with possible rate limiting by API endpoints using _ratelimitr_.
This is in place of internally rate-limiting due to the way _ratelimitr_ handles function creation and the fact that the rate limits are extremely generous and may change as the project matures.

# nasapower 4.0.1 (unreleased on CRAN)

## Bug fixes

* Fixes a bug in where `NA` values were improperly handled.
Thanks to [@femiguez](https://github.com/femiguez) for the [Pull Request with the fix](https://github.com/ropensci/nasapower/pull/56).

## Minor changes

* Enforces API limits client-side where the API limits unique queries to 30 per 60 seconds as found and reported by [@camwur](https://github.com/camwur) in [Issue 57](https://github.com/ropensci/nasapower/issues/57).
This can be adjusted in future releases of _nasapower_ if the POWER API changes as has been indicated is possible.

* (Re)enables _vcr_ for better unit testing.

* More comprehensive unit tests.

# nasapower 4.0.0

## Major changes

* Adds support for new NASA POWER API v2.0, which includes new hourly data and other major changes to the API and available data.
See <https://power.larc.nasa.gov/> for fully detailed changes to the data.

* Drops support for the deprecated NASA POWER API V1.0.
Previous versions of _nasapower_ are no longer functional.

* Adds new function, `query_parameters()` to fetch information from the API on individual and all available community/temporal API combination parameters.

* Removes `SSE` community, replaced with `RE`.

* Removes `global` option for geographic coverage as passed along through the `latlon` argument of `get_power()`.

* Directly parse data from API response rather than downloading data to disk and importing.

* The `get_power()` arguments are changed:
  * two new arguments are added,
    * `wind_elevation`, and
    * `wind_surface`,
  * the `temporal_average` argument has been superseded by `temporal_api` to align with the terminology used in the POWER API.
  The `temporal_average` argument will still work, however, a message will be given if a user still uses `temporal_average` to alert the user of the change and ask them to update their scripts.


## Minor changes

* Improved documentation.

* Removes internal references to ICASA format files that are no longer supported in this client.

# nasapower 3.0.1

## Bug fixes

* Fix bug where Solar Radiation, "ALLSKY_SFC_SW_DWN", and perhaps others that were missed, return a numeric `-99.00` value rather than the proper `NA` for missing data.
Thanks to Fernando Miguez, <https://github.com/femiguez>, for the assistance in isolating the issue.

# nasapower 3.0.0

## Major Changes to Functionality

* Due to the removal of the CRAN package _APSIM_ from CRAN, the removal of the `create_met()` function has been implemented sooner than expected to keep _nasapower_ on CRAN.

* Deprecates `create_met()`

## Bug fixes

* Properly deprecates `create_icasa()`

# nasapower 2.0.0

## Bug Fixes

* Correct any missing or redirecting URLs

* Replace deprecated `subclass` with `class` in `new_tibble()`

## Major Changes to Functionality

* Following a UNIX-like philosophy, this release removes functionality to write APSIM .met and DSSAT ICASA files to disk.
_nasapower_ now will only fetch the appropriate data and return a `tibble()` object in-session, please use [apsimx](https://cran.r-project.org/package=apsimx) or the POWER web API data access viewer, <https://power.larc.nasa.gov/data-access-viewer/>, for fetching and/or writing .met or .icasa files, respectively.
Note that  `create_icasa()` ideally should have been deprecated, but the server was not responding properly when queried for some time before the current release of _nasapower_ so the function has been removed.

* Add ability to `get_power()` to accept a user-provided `site_elevation` parameter that is passed to the API.
When this is used it will return a corrected atmospheric pressure value adjusted to the elevation provided.

## Minor and Internal Changes

* Use newest values from POWER team to validate user inputs for API requests, see <https://github.com/ropensci/nasapower/issues/48> for more.

* Replace _raster_ with _terra_ for examples of converting to spatial data in vignettes

* Use _vcr_ for enhanced testing

* Refactor the internal handling of temporary files to allow for more efficient use of the _future_ package


# nasapower 1.1.3

## Bug Fixes

* Corrects bug when querying the SB or SSE communities resulting in an error

* Corrects example in vignette when creating a .met file

## Minor Changes

* Update documentation to use ROxygen 7.0.0

* Add new vignette, "Using nasapower with large geographic areas"

# nasapower 1.1.2

# Minor changes

- Correct URL in BibTeX version of citation

- Suppress output in console from `APSIM::createMetFile()`

- Help file titles are now in sentence case

# nasapower 1.1.1

## Bug fixes

- Fix issues reported at https://cloud.r-project.org//web/checks/check_results_nasapower.html with failing tests.
These tests should be skipped on CRAN but were not.

- Fixes bug where missing values in POWER data were not properly replaced with `NA` in `tibble` and metFile outputs

- Fixes bug in documentation for `create_icasa()` where the parameter for `file_out` was misidentified as just `file`

## Minor changes

- Users are now notified if creating a .met file that has any missing values through a console message and .csv file being written to disk to accompany the resulting .met file describing which values are missing

# nasapower 1.1.0

## Bug fixes

- Fixes bug where .met files were not created properly including where "radn" and "rain" col headers were reversed

- Fix `Warning: Must pass a scalar integer as 'nrow' argument to 'new_tibble()'.`

- Fixes bug where "CLIMATE" could not be requested for a single point

## Major changes

- Change how `GLOBAL` values are requested. This is now specified in `lonlat` in conjunction with `temporal_average = CLIMATOLOGY`.

## Minor changes

- Adds example of fetching climate for a single point

- Refactor code to split internal functions by functionality and add more complete test coverage

# nasapower 1.0.7

## Minor changes

- Removes internal check for data - community agreement, as all data is available for all communities, only the units change

- Update links to latest documentation provided by the POWER team

# nasapower 1.0.6

## Minor changes

- Adds support for WS2M_MIN, WS2M_MAX and WS2M_RANGE in AG community

## Bug fixes

- Fixes bug where previous release did not support WS2M from AG community due to a local typo

# nasapower 1.0.5

## Minor changes

- "Fixes" [Issue 32](https://github.com/ropensci/nasapower/issues/32) where WS2M is not available through `nasapower` until the POWER team can properly address how pre-query validation should be performed

# nasapower 1.0.4

## Minor changes

- Corrects an instance where vignette example executed on CRAN but should not

- Adds link to POWER website in error message when query fails

- Documentation .Rd files are now more readable with better formatting

# nasapower 1.0.3

## Minor changes

- Adds citation information for JOSS paper, http://joss.theoj.org/papers/10.21105/joss.01035

## Documentation changes

- Flesh out examples using `naspower` data with `raster` to create spatial objects for systems with low-RAM where the functionality may not work as expected

- Standardise formatting of vignette subheadings

- Spell check vignette

## Bug fixes

- Fixes tests to not run on CRAN so that errors aren't reported when API is unavailable

# nasapower 1.0.2

## Minor changes

- Updates documentation examples

- Provides nicer method of printing data in R console

- Updates tests for better coverage and removes non-functional tests

- Removes `dplyr` as an Import

## Bug fixes

- Corrects issue where `if()` was called with a vector of length 2 or more

- Corrects logical operators `&&` and `||` where they should be `&` or `|`

- Removes extra code in `create_icasa()` and `create_met()` that performed a duplicated check of `latlon` values

- Removes unnecessary checks for `latlon` in `get_power()`

# nasapower 1.0.1

## Minor changes

- Provides corrections to documentation formatting as per CRAN volunteers' requests

- Provides edits and clarifications in DESCRIPTION's Description and Title about the package's uses and capabilities

# nasapower 1.0.0 (unreleased)

## Major changes

- _nasapower_ is now a part of [rOpenSci](https://ropensci.org/) after [peer-review of the code](https://github.com/ropensci/software-review/issues/155)!

- Provides access to all three communities, AG, SSE and and SB, not just AG

- Uses new 'POWER' 'API' to download new 1/2 x 1/2 degree data

- Adds function `get_power()` to get weather data and optionally metadata as well

- Adds function `create_met()` to create 'APSIM' met objects from 'POWER' data

- Adds function `create_icasa()` to create a text file of weather data for use in 'DSSAT' crop modelling

- Internally, replaces _httr_ package with _crul_

### Deprecated functions

- The `get_cell` and `get_region` functions are deprecated in favour of `get_power()`.
The new POWER interface allows for the specification of single points or regional areas.
Global coverage may be queried for Climatology.
See the help for `?get_power()` for more details.

# nasapower 0.1.4

### Bug Fixes

- Fixes bug related to date columns where `MONTH`, `DAY` and `YYYY-MM-DD` were incorrectly reported in final data frame. This did not affect the weather data, `YEAR` or `DOY` columns.

# nasapower 0.1.3

### Bug fixes

- Fix bug where lon/lat values were improperly assigned internally due to row names not being ordered correctly in `get_region()`

- Fix bug reports link in DESCRIPTION file

- Correct vignette where it had said, "both of which will which will download"

- Correct documentation for `get_region()`, which incorrectly stated that it downloaded data for a 1 x 1 degree cell

### Minor improvements

- Optimise arguments used in `utils::read.table()` to ingest weather data in the `get_cell()` and `get_region()` functions more quickly

- NEWS now formatted more nicely for easier reading

- Add statement about possible performance and memory usage when using
`get_region()` in the vignette

- Add an example of converting the data frame to a spatial object using
_raster_ to create a `raster::brick()`

- Specify in documentation that a range of days to years can be specified for download

## Minor changes

- `get_region()` and `get_cell()` now default to download all weather vars

- Add a check to see if POWER website is responding before making request for data. If not, stop and return error message to user.

# nasapower 0.1.2

### Bug fixes

- Fixes bug where only first date is reported when using `get_region()` with multiple dates. https://github.com/ropensci/nasapower/issues/1

### Minor improvements

- Enhanced documentation

- Superfluous function, `.onLoad()`, removed from zzz.R

- Tidied up startup message

- Clean up vignette

- Build vignette faster

- Remove DATE from DESCRIPTION

# nasapower 0.1.1

### Minor improvements

- Fix issues in documentation, typos, incorrect links, etc.

# nasapower 0.1.0

### New features

* Add new functionality to download regions in addition to single cells

* Add static documentation website, <https://docs.ropensci.org/nasapower/>

* Add startup message

### Minor improvements

* Better documentation

# nasapower 0.0.2

### New features

* Added citation file

# nasapower 0.0.1

* Added a `NEWS.md` file to track changes to the package.

* First release, no changes to report yet
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

Examples of behaviour that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behaviour include:

* The use of sexualised language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behaviour and will take appropriate and fair corrective action in
response to any behaviour that they deem inappropriate, threatening, offensive,
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

Instances of abusive, harassing, or otherwise unacceptable behaviour may be
reported to the community leaders responsible for enforcement at
adamhsparks@gmail.com. All complaints will be reviewed and investigated promptly
and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behaviour deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behaviour was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behaviour. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behaviour.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behaviour, harassment of an
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
# nasapower v4.0.3

This submission supersedes v4.0.2, which should not be released to CRAN due to issues in the tests.
This submission skips all tests requiring API connectivity on CRAN fixing this issue introduced with 4.0.2.
Sorry for the inconvenience, but I thought this was better to submit a new(er) version rather than causing failures on CRAN servers.

## Test environments
* local macOS, Platform: aarch64-apple-darwin20 (64-bit), R 4.1.2
* win-builder, R Under development (unstable) (2021-12-17 r81389 ucrt)
* win-builder, R 4.1.2

## R CMD check results

0 errors | 0 warnings | 1 note

This is a new patch release

## Reverse dependencies

No ERRORs or WARNINGs were found
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

Please note that the nasapower project is released with a
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
---
name: Bug report
about: Create a report to help us improve nasapower
title: ''
labels: ''
assignees: ''

---

**Please do not submit bug reports for issues with the NASA POWER server itself. This package is only an API client, it does not have control over the server hosting the data. If you receive error messages when retrieving data, please check <https://power.larc.nasa.gov/> for any information on the server status first.

**Describe the bug**
A clear and concise description of what the bug is.

**Code to Reproduce**
A clear set of code illustrating the steps to reproduce the behaviour.

**Expected behaviour**
A clear and concise description of what you expected to happen.

**OS and R versions (please complete the following information):**
 - OS: [e.g., Windows10, macOS, Linux]
 - R Version [e.g., 3.6.2]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. 

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
Fetch NASA-POWER Parameters and Include Them as an Internal List
================
Adam H. Sparks
2021-12-05

# Create parameters nested list for internal checks before sending queries to POWER server

These data are used for internal checks to be sure that data requested
from the POWER dataset are valid. The POWER list of parameters that can
be queried is available as an API query from the POWER server.

The list structure will be

-   `parameters`
    -   `HOURLY_AG`
        -   `parameter_1` …
        -   `parameter_n`
    -   `DAILY_AG`
        -   `parameter_1` …
        -   `parameter_n`
    -   `MONTHLY_AG`
        -   `parameter_1` …
        -   `parameter_n`
    -   `CLIMATOLOGY_AG`
        -   `parameter_1` …
        -   `parameter_n`
    -   `HOURLY_RE`
        -   `parameter_1` …
        -   `parameter_n`
    -   `DAILY_RE`
        -   `parameter_1` …
        -   `parameter_n`
    -   `MONTHLY_RE`
        -   `parameter_1` …
        -   `parameter_n`
    -   `CLIMATOLOGY_RE`
        -   `parameter_1` …
        -   `parameter_n`
    -   `HOURLY_SB`
        -   `parameter_1` …
        -   `parameter_n`
    -   `DAILY_SB`
        -   `parameter_1` …
        -   `parameter_n`
    -   `MONTHLY_SB`
        -   `parameter_1` …
        -   `parameter_n`
    -   `CLIMATOLOGY_SB`
        -   `parameter_1` …
        -   `parameter_n`

## POWER JSON file

Using `purrr::map2` and `jsonlite::fromJSON()` read the JSON file into R
creating a single, nested list and reorder it alphabetically by
parameter name.

``` r
library(purrr)
library(jsonlite)
```

    ## 
    ## Attaching package: 'jsonlite'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     flatten

``` r
temporal_api <- c("HOURLY", "DAILY", "MONTHLY", "CLIMATOLOGY")
community <- c("AG", "RE", "SB")

# create all combinations
vals <- expand.grid(temporal_api, community, stringsAsFactors = FALSE)

# create equal length vectors for purrr::map2
temporal_api <- vals[, 1]
community <- vals[, 2]

base_url <- "https://power.larc.nasa.gov/api/system/manager/parameters?"

power_pars <- map2(.x = temporal_api,
                  .y = community,
                  .f = ~ fromJSON(paste0(
                    base_url,
                    "temporal=",
                    .x,
                    "&community=",
                    .y
                  )))

names(power_pars) <- paste(temporal_api, community, sep = "_")

# create a list of vectors for each temporal API/par combination for easier
# checking and validation
parameters <- vector(mode = "list", length = length(power_pars))
names(parameters) <- names(power_pars)
for (i in names(power_pars)) {
  parameters[[i]] <- names(power_pars[[i]])
}

parameters <- map(.x = parameters, .f = sort)
```

## Save list for use in `nasapower` package

Using `usethis::use_data()` save the list as an R data object for use in
the *nasapower* package. These values will not be exposed to the user
and so will not be documented as previously. Users will be be pointed to
functions to interact directly with the POWER APIs to query information
for the temporal API/community combinations or for the parameters
themselves.

``` r
usethis::use_data(parameters, overwrite = TRUE, internal = TRUE)
```

## Session Info

``` r
sessioninfo::session_info()
```

    ## ─ Session info  😂  🧾  🤚🏿   ─────────────────────────────────────────────────
    ##  hash: face with tears of joy, receipt, raised back of hand: dark skin tone
    ## 
    ##  setting  value
    ##  version  R version 4.1.2 (2021-11-01)
    ##  os       macOS Monterey 12.0.1
    ##  system   aarch64, darwin20
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_AU.UTF-8
    ##  ctype    en_AU.UTF-8
    ##  tz       Australia/Perth
    ##  date     2021-12-05
    ##  pandoc   2.14.0.3 @ /Applications/RStudio.app/Contents/MacOS/pandoc/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date (UTC) lib source
    ##  cli           3.1.0   2021-10-27 [1] CRAN (R 4.1.1)
    ##  crayon        1.4.2   2021-10-29 [1] CRAN (R 4.1.1)
    ##  curl          4.3.2   2021-06-23 [1] CRAN (R 4.1.0)
    ##  desc          1.4.0   2021-09-28 [1] CRAN (R 4.1.1)
    ##  digest        0.6.29  2021-12-01 [1] CRAN (R 4.1.2)
    ##  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.1.0)
    ##  evaluate      0.14    2019-05-28 [1] CRAN (R 4.1.0)
    ##  fansi         0.5.0   2021-05-25 [1] CRAN (R 4.1.0)
    ##  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.1.0)
    ##  fs            1.5.1   2021-11-30 [1] CRAN (R 4.1.2)
    ##  glue          1.5.1   2021-11-30 [1] CRAN (R 4.1.2)
    ##  htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.1.1)
    ##  jsonlite    * 1.7.2   2020-12-09 [1] CRAN (R 4.1.0)
    ##  knitr         1.36    2021-09-29 [1] CRAN (R 4.1.1)
    ##  lifecycle     1.0.1   2021-09-24 [1] CRAN (R 4.1.1)
    ##  magrittr      2.0.1   2020-11-17 [1] CRAN (R 4.1.0)
    ##  pillar        1.6.4   2021-10-18 [1] CRAN (R 4.1.1)
    ##  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.1.0)
    ##  purrr       * 0.3.4   2020-04-17 [1] CRAN (R 4.1.2)
    ##  R6            2.5.1   2021-08-19 [1] CRAN (R 4.1.1)
    ##  rlang         0.4.12  2021-10-18 [1] CRAN (R 4.1.1)
    ##  rmarkdown     2.11    2021-09-14 [1] CRAN (R 4.1.1)
    ##  rprojroot     2.0.2   2020-11-15 [1] CRAN (R 4.1.0)
    ##  rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.1.0)
    ##  sessioninfo   1.2.1   2021-11-02 [1] CRAN (R 4.1.2)
    ##  stringi       1.7.6   2021-11-29 [1] CRAN (R 4.1.2)
    ##  stringr       1.4.0   2019-02-10 [1] CRAN (R 4.1.1)
    ##  tibble        3.1.6   2021-11-07 [1] CRAN (R 4.1.1)
    ##  usethis       2.1.3   2021-10-27 [1] CRAN (R 4.1.1)
    ##  utf8          1.2.2   2021-07-24 [1] CRAN (R 4.1.0)
    ##  vctrs         0.3.8   2021-04-29 [1] CRAN (R 4.1.0)
    ##  xfun          0.28    2021-11-04 [1] CRAN (R 4.1.2)
    ##  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.1.0)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
---
title: 'nasapower: A NASA POWER Global Meteorology, Surface Solar Energy and Climatology Data Client for R'
authors:
  - affiliation: 1
    name: Adam H. Sparks
    orcid: 0000-0002-0061-8359
date: 19 Oct 2018
output: pdf_document
bibliography: paper.bib
tags:
  - NASA
  - weather data
  - solar data
  - climate data
  - meteorology
  - agroclimatology
  - climatology
  - alternative energy
  - sustainable buildings
  - R
  - earth science
  - reproducibility
affiliations:
  - name: University of Southern Queensland, Centre for Crop Health, Toowoomba Queensland 4350, Australia
    index: 1

---

# Summary and Statement of Need

_nasapower_ is an R [@RCT2018] package providing functionality to interface with
the NASA POWER API [@StackhouseJr2018] for reproducible data retrieval using R.
Three functions, `get_power()`, `create_met()` and `create_icasa()` are
provided. The `get_power()` function provides complete access to all
functionality that the POWER API provides, which includes three user
communities, AG (agroclimatology), SSE (Surface meteorology and Solar Energy)
and SB (Sustainable Buildings); three temporal averages, Daily, Interannual and
Climatology; three geographic options, single point, regional and global for the
appropriate parameters offered. _nasapower_ uses _lubridate_ [@Grolemund2011]
internally to format and parse dates which are passed along to the the query
constructed using _crul_ [@Chamberlain2018] to interface with the POWER API. The
query returns a json response, which is parsed by _jsonlite_ [@Ooms2014] to
obtain the url of the .csv file that has been requested. The .csv file is
downloaded to local disk using _curl_ [@Ooms2018] and read into R using _readr_
[@Wickham2017]. Data are returned in a tidy data frame [@Wickham2014] as a
_tibble_ [@Mueller2018] with a custom header, which provides POWER metadata. Two
other functions provide functionality to generate weather input files for
agricultural crop modelling. The `create_met()` function is a wrapper for the
`get_power()` function coupled with the `prepareMet()` and `writeMet()`
functions from _APSIM_ [@Fainges2017] to simplify the process of querying the
data and creating text files in the .met format for use in Agricultural
Production Systems sIMulator (APSIM). While the `create_icasa()` function wraps
the `get_power()` into a function that generates and locally saves a text file
in the International Consortium for Agricultural Systems Applications (ICASA)
format for use in the Decision Support System for Agrotechnology Transfer
(DSSAT) framework [@Jones2003; @Hoogenboom2017]. Extended documentation is
provided with examples of converting it to spatial objects using _raster_
[@Hijmans2017].

Integrating this data retrieval and formatting in R will streamline processes
with models such as APSIM [@Keating2003], DSSAT
[@Jones1998; @Jones2003] and EPIRICE [@Savary2012] that can be
linked to or are implemented fully in the R programming language.

# About POWER Data

NASA’s POWER (Prediction Of Worldwide Energy Resource) data [@StackhouseJr2018]
are freely available for download via a
[web interface](https://power.larc.nasa.gov/data-access-viewer/) at a
grid resolution of one-half arc degree longitude by one-half arc degree
latitude. Funded through the NASA Earth Science Directorate Applied Science
Program, the data provide daily global coverage from 1983 until near present for
all parameters except precipitation, which is provided for January 1997 to near
present with a several month delay. The data are widely used in agricultural
modelling for modelling crop yields [@Bai2010; @vanWart2013;
@vanWart2015], other crop simulation exercises [@Ojeda2017], plant disease
modelling [@Savary2012].

While _nasapower_ does not redistribute any of the NASA POWER data, users are
encouraged to please refer to the acknowledgement guidelines available at,
<https://power.larc.nasa.gov/#contact> and properly acknowledge the data as
requested.

> When POWER data products are used in a publication, we request the following
acknowledgment be included: "_These data were obtained from the NASA Langley
Research Center POWER Project funded through the NASA Earth Science Directorate
Applied Science Program._"

# References
## revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                                |
|:--------|:------------------------------------|
|version  |R version 4.1.2 (2021-11-01)         |
|os       |macOS Monterey 12.0.1                |
|system   |x86_64, darwin17.0                   |
|ui       |RStudio                              |
|language |(EN)                                 |
|collate  |en_AU.UTF-8                          |
|ctype    |en_AU.UTF-8                          |
|tz       |Australia/Perth                      |
|date     |2021-12-12                           |
|rstudio  |2021.09.0+351 Ghost Orchid (desktop) |
|pandoc   |2.13 @ /usr/local/bin/pandoc         |

# Dependencies

|package   |old   |new   |Δ  |
|:---------|:-----|:-----|:--|
|nasapower |4.0.1 |4.0.2 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output: github_document
---
# _nasapower_: NASA POWER API Client <img align='right' src='man/figures/logo.png'>

<!-- badges: start -->
[![tic](https://github.com/ropensci/nasapower/workflows/tic/badge.svg?branch=main)](https://github.com/ropensci/nasapower/actions)
[![codecov](https://codecov.io/gh/ropensci/nasapower/branch/master/graph/badge.svg?token=Kq9aea0TQN)](https://codecov.io/gh/ropensci/nasapower)
[![DOI](https://zenodo.org/badge/109224461.svg)](https://zenodo.org/badge/latestdoi/109224461)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![peer-review](https://badges.ropensci.org/155_status.svg)](https://github.com/ropensci/software-review/issues/155)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01035/status.svg)](https://doi.org/10.21105/joss.01035)
[![CRAN](http://www.r-pkg.org/badges/version/nasapower)](https://CRAN.R-project.org/package=nasapower)
<!-- badges: end -->

_nasapower_ aims to make it quick and easy to automate downloading [NASA-POWER](https://power.larc.nasa.gov) global meteorology, surface solar energy and climatology data in your R session as a tidy data frame `tibble` object for analysis and use in modelling or other purposes.
POWER (Prediction Of Worldwide Energy Resource) data are freely available for download with varying spatial resolutions dependent on the original data and with several temporal resolutions depending on the POWER parameter and community.

**Note that the data are not static and may be replaced with improved data.** Please see <https://power.larc.nasa.gov/docs/services/> for detailed information in this regard.

### Quick start

_nasapower_ can easily be installed using the following code.

#### From CRAN

The stable version is available through CRAN.

```{r install-cran, eval=FALSE}
install.packages("nasapower")
```

#### From GitHub for the version in-development

A development version is available through GitHub.

```{r install-dev, eval=FALSE}
if (!require("remotes")) {
  install.packages("remotes")
}

remotes::install_github("ropensci/nasapower")
```

### Example

Fetch daily “ag” community temperature, relative humidity and precipitation for January 1, 1985 for Kingsthorpe, Queensland, Australia.

```{r kingsthorpe}
library("nasapower")
daily_ag <- get_power(community = "ag",
                      lonlat = c(151.81, -27.48),
                      pars = c("RH2M", "T2M", "PRECTOTCORR"),
                      dates = "1985-01-01",
                      temporal_api = "daily"
                      )
daily_ag
```

## Documentation

More documentation is available in the vignette in your R session, `vignette("nasapower")` or available online, <https://docs.ropensci.org/nasapower/>.

## Use of POWER Data

While _nasapower_ does not redistribute the data or provide it in any way, we encourage users to follow the requests of the POWER Project Team.

> When POWER data products are used in a publication, we request the
  following acknowledgement be included: “These data were obtained from
  the NASA Langley Research Center POWER Project funded through the NASA
  Earth Science Directorate Applied Science Program.”

## Meta

  - Please [report any issues or bugs](https://github.com/ropensci/nasapower/issues).

  - License: MIT

  - To cite _nasapower_, please use the output from `citation(package = "nasapower")`.

  - Please note that the _nasapower_ project is released with a [Contributor Code of Conduct](https://github.com/ropensci/nasapower/blob/main/CODE_OF_CONDUCT.md).
    By participating in the _nasapower_ project you agree to abide by its terms.

  - The U.S. Earth System Research Laboratory, Physical Science Division of the National Atmospheric & Oceanic Administration (NOAA) maintains a list of gridded climate data sets that provide different data and different resolutions <https://psl.noaa.gov/data/gridded/>.

## References

<https://power.larc.nasa.gov>

<https://power.larc.nasa.gov/docs/methodology/>
---
title: "Fetch NASA-POWER Parameters and Include Them as an Internal List"
author: "Adam H. Sparks"
date: "`r format(Sys.Date())`"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Create parameters nested list for internal checks before sending queries to POWER server

These data are used for internal checks to be sure that data requested from the POWER dataset are valid.
The POWER list of parameters that can be queried is available as an API query from the POWER server.

The list structure will be

- `parameters`
  - `HOURLY_AG`
    - `parameter_1` ...
    - `parameter_n`
  - `DAILY_AG`
    - `parameter_1` ...
    - `parameter_n`
  - `MONTHLY_AG`
    - `parameter_1` ...
    - `parameter_n`
  - `CLIMATOLOGY_AG`
    - `parameter_1` ...
    - `parameter_n`
  - `HOURLY_RE`
    - `parameter_1` ...
    - `parameter_n`
  - `DAILY_RE`
    - `parameter_1` ...
    - `parameter_n`
  - `MONTHLY_RE`
    - `parameter_1` ...
    - `parameter_n`
  - `CLIMATOLOGY_RE`
    - `parameter_1` ...
    - `parameter_n`
  - `HOURLY_SB`
    - `parameter_1` ...
    - `parameter_n`
  - `DAILY_SB`
    - `parameter_1` ...
    - `parameter_n`
  - `MONTHLY_SB`
    - `parameter_1` ...
    - `parameter_n`
  - `CLIMATOLOGY_SB`
    - `parameter_1` ...
    - `parameter_n`

## POWER JSON file

Using `purrr::map2` and `jsonlite::fromJSON()` read the JSON file into R creating a single, nested list and reorder it alphabetically by parameter name.

```{r fetch-JSON}
library(purrr)
library(jsonlite)

temporal_api <- c("HOURLY", "DAILY", "MONTHLY", "CLIMATOLOGY")
community <- c("AG", "RE", "SB")

# create all combinations
vals <- expand.grid(temporal_api, community, stringsAsFactors = FALSE)

# create equal length vectors for purrr::map2
temporal_api <- vals[, 1]
community <- vals[, 2]

base_url <- "https://power.larc.nasa.gov/api/system/manager/parameters?"

power_pars <- map2(.x = temporal_api,
                  .y = community,
                  .f = ~ fromJSON(paste0(
                    base_url,
                    "temporal=",
                    .x,
                    "&community=",
                    .y
                  )))

names(power_pars) <- paste(temporal_api, community, sep = "_")

# create a list of vectors for each temporal API/par combination for easier
# checking and validation
parameters <- vector(mode = "list", length = length(power_pars))
names(parameters) <- names(power_pars)
for (i in names(power_pars)) {
  parameters[[i]] <- names(power_pars[[i]])
}

parameters <- map(.x = parameters, .f = sort)
```

## Save list for use in `nasapower` package

Using `usethis::use_data()` save the list as an R data object for use in the _nasapower_ package.
These values will not be exposed to the user and so will not be documented as previously.
Users will be be pointed to functions to interact directly with the POWER APIs to query information for the temporal API/community combinations or for the parameters themselves.

```{r save-list, message=FALSE}
usethis::use_data(parameters, overwrite = TRUE, internal = TRUE)
```

## Session Info

```{r session-info}
sessioninfo::session_info()
```
---
title: "nasapower"
author: "Adam H. Sparks"
output:
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{nasapower}
  %\VignetteEngine{knitr::rmarkdown_notangle}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{ratelimitr}
---



## Introduction

_nasapower_ aims to make it quick and easy to automate downloading NASA [POWER](https://power.larc.nasa.gov) global meteorology, surface solar energy and climatology data  data in your R session as a tidy data frame for analysis and use in modelling or other purposes using `get_power()`.
POWER (Prediction Of Worldwide Energy Resource) data are freely available for download through a web interface with a spatial resolution of 0.5 x 0.625 degree latitude and longitude for meteorology and  1 x 1 degree latitude and longitude for solar parameters with various temporal resolutions depending on the POWER parameter and community.

## Using get_power() to fetch POWER data

The `get_power()` function has eight possible arguments and returns a data frame with a metadata header in the current R session.

### Example fetching daily data for a single point

Fetch daily "AG" community temperature, relative humidity and precipitation for January 1985 for Kingsthorpe, Queensland, Australia.


```r
library("nasapower")
daily_single_ag <- get_power(
  community = "ag",
  lonlat = c(151.81, -27.48),
  pars = c("RH2M", "T2M", "PRECTOTCORR"),
  dates = c("1985-01-01", "1985-01-31"),
  temporal_api = "daily"
)

daily_single_ag
#> NASA/POWER CERES/MERRA2 Native Resolution Daily Data  
#>  Dates (month/day/year): 01/01/1985 through 01/31/1985  
#>  Location: Latitude  -27.48   Longitude 151.81  
#>  Elevation from MERRA-2: Average for 0.5 x 0.625 degree lat/lon region = 442.77 meters 
#>  Value for missing model data cannot be computed or out of model availability range: NA  
#>  Parameter(s):  
#>  
#>  Parameters: 
#>  RH2M            MERRA-2 Relative Humidity at 2 Meters (%) ;
#>  T2M             MERRA-2 Temperature at 2 Meters (C) ;
#>  PRECTOTCORR     MERRA-2 Precipitation Corrected (mm/day)  
#>  
#> # A tibble: 31 × 10
#>      LON   LAT  YEAR    MM    DD   DOY YYYYMMDD    RH2M   T2M
#>    <dbl> <dbl> <dbl> <int> <int> <int> <date>     <dbl> <dbl>
#>  1  152. -27.5  1985     1     1     1 1985-01-01  54.7  24.9
#>  2  152. -27.5  1985     1     2     2 1985-01-02  42.1  28.6
#>  3  152. -27.5  1985     1     3     3 1985-01-03  43.4  27.4
#>  4  152. -27.5  1985     1     4     4 1985-01-04  48.9  24.3
#>  5  152. -27.5  1985     1     5     5 1985-01-05  55.3  26.5
#>  6  152. -27.5  1985     1     6     6 1985-01-06  60.2  27.0
#>  7  152. -27.5  1985     1     7     7 1985-01-07  63.1  27.2
#>  8  152. -27.5  1985     1     8     8 1985-01-08  70.6  24.9
#>  9  152. -27.5  1985     1     9     9 1985-01-09  60    26.1
#> 10  152. -27.5  1985     1    10    10 1985-01-10  45.2  27.0
#> # … with 21 more rows, and 1 more variable: PRECTOTCORR <dbl>
```

### Example fetching daily data for an area

Fetch daily "ag" community relative humidity and temperature for south east Queensland region.


```r
daily_region_ag <- get_power(
  community = "ag",
  lonlat = c(150.5, -28.5 , 153.5, -25.5),
  pars = c("RH2M", "T2M"),
  dates = c("1985-01-01", "1985-01-02"),
  temporal_api = "daily"
)

daily_region_ag
#> NASA/POWER CERES/MERRA2 Native Resolution Daily Data  
#>  Dates (month/day/year): 01/01/1985 through 01/02/1985  
#>  Location: Regional  
#>  Elevation from MERRA-2: Average for 0.5 x 0.625 degree lat/lon region = na meters 
#>  Value for missing model data cannot be computed or out of model availability range: NA  
#>  Parameter(s):  
#>  
#>  Parameters: 
#>  RH2M     MERRA-2 Relative Humidity at 2 Meters (%) ;
#>  T2M      MERRA-2 Temperature at 2 Meters (C)  
#>  
#> # A tibble: 72 × 9
#>      LAT   LON  YEAR    MM    DD   DOY YYYYMMDD    RH2M   T2M
#>    <dbl> <dbl> <dbl> <int> <int> <int> <date>     <dbl> <dbl>
#>  1 -28.2  151.  1985     1     1     1 1985-01-01  43.6  26.5
#>  2 -28.2  151.  1985     1     1     1 1985-01-01  44.4  25.8
#>  3 -28.2  152.  1985     1     1     1 1985-01-01  52.6  24.0
#>  4 -28.2  152.  1985     1     1     1 1985-01-01  57.7  23.9
#>  5 -28.2  153.  1985     1     1     1 1985-01-01  61.4  24.9
#>  6 -28.2  153.  1985     1     1     1 1985-01-01  66.1  26.0
#>  7 -27.8  151.  1985     1     1     1 1985-01-01  45.8  26.5
#>  8 -27.8  151.  1985     1     1     1 1985-01-01  47.9  26.0
#>  9 -27.8  152.  1985     1     1     1 1985-01-01  53.4  24.8
#> 10 -27.8  152.  1985     1     1     1 1985-01-01  56.7  25.1
#> # … with 62 more rows
```

### Example fetching interannual data for an area

Fetch interannual solar cooking parameters for south east Queensland region.


```r
interannual_re <- get_power(
  community = "re",
  lonlat = c(150.5, -28.5 , 153.5, -25.5),
  dates = c("1984", "1985"),
  temporal_api = "monthly",
  pars = c("CLRSKY_SFC_SW_DWN",
           "ALLSKY_SFC_SW_DWN")
)

interannual_re
#> NASA/POWER CERES/MERRA2 Native Resolution Monthly and Annual  
#>  Dates (month/day/year): 01/01/1984 through 12/31/1985  
#>  Location: Regional  
#>  Elevation from MERRA-2: Average for 0.5 x 0.625 degree lat/lon region = na meters 
#>  Value for missing model data cannot be computed or out of model availability range: NA  
#>  Parameter(s):  
#>  
#>  Parameters: 
#>  ALLSKY_SFC_SW_DWN     CERES SYN1deg All Sky Surface Shortwave Downward Irradiance (kW-hr/m^2/day) ;
#>  CLRSKY_SFC_SW_DWN     CERES SYN1deg Clear Sky Surface Shortwave Downward Irradiance (kW-hr/m^2/day)  
#>  
#> # A tibble: 144 × 17
#>    PARAMETER    YEAR   LAT   LON   JAN   FEB   MAR   APR   MAY
#>    <chr>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 ALLSKY_SFC…  1984 -25.8  151.  6.01  6.49  5.79  4.67  4.12
#>  2 ALLSKY_SFC…  1984 -25.8  151.  5.92  5.97  5.64  4.37  4.01
#>  3 ALLSKY_SFC…  1984 -25.8  152.  5.92  5.97  5.64  4.37  4.01
#>  4 ALLSKY_SFC…  1984 -25.8  152.  5.96  5.85  5.56  4.26  3.92
#>  5 ALLSKY_SFC…  1984 -25.8  153.  5.96  5.85  5.56  4.26  3.92
#>  6 ALLSKY_SFC…  1984 -25.8  153.  6.23  6.05  5.88  4.26  3.81
#>  7 ALLSKY_SFC…  1984 -26.2  151.  5.97  6.65  6     4.66  4.02
#>  8 ALLSKY_SFC…  1984 -26.2  151.  6     6.38  5.71  4.38  4.01
#>  9 ALLSKY_SFC…  1984 -26.2  152.  6     6.38  5.71  4.38  4.01
#> 10 ALLSKY_SFC…  1984 -26.2  152.  5.75  5.96  5.37  4.13  3.8 
#> # … with 134 more rows, and 8 more variables: JUN <dbl>,
#> #   JUL <dbl>, AUG <dbl>, SEP <dbl>, OCT <dbl>, NOV <dbl>,
#> #   DEC <dbl>, ANN <dbl>
```

### Example fetching climatology data

Climatology data can be retrieved for point or regional areas as demonstrated previously.
Change the `temporal_api` value to "climatology" to get these data.

Fetch "ag" climatology for temperature and relative humidity for Kingsthorpe, Queensland, Australia.


```r
climatology_ag <- get_power(
  community = "ag",
  pars = c("T2M", "RH2M"),
  lonlat = c(151.81, -27.48),
  temporal_api = "climatology"
)

climatology_ag
#> NASA/POWER CERES/MERRA2 Native Resolution Climatology Climatologies  
#>  30-year Meteorological and Solar Monthly & Annual Climatologies (January 1990 - December 2019)  
#>  Location: Latitude  -27.48   Longitude 151.81  
#>  Elevation from MERRA-2: Average for 0.5 x 0.625 degree lat/lon region = 442.77 meters 
#>  Value for missing model data cannot be computed or out of model availability range: NA  
#>  Parameter(s):  
#>  
#>  Parameters: 
#>  T2M      MERRA-2 Temperature at 2 Meters (C) ;
#>  RH2M     MERRA-2 Relative Humidity at 2 Meters (%)  
#>  
#> # A tibble: 2 × 16
#>     LON   LAT PARAMETER   JAN   FEB   MAR   APR   MAY   JUN
#>   <dbl> <dbl> <chr>     <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  152. -27.5 T2M        24.3  23.5  21.8  18.7  15.1  12.2
#> 2  152. -27.5 RH2M       65.5  70.1  70.4  70.1  71.8  75.5
#> # … with 7 more variables: JUL <dbl>, AUG <dbl>, SEP <dbl>,
#> #   OCT <dbl>, NOV <dbl>, DEC <dbl>, ANN <dbl>
```

_Note_ the associated metadata in the data frame header are not saved if the data are exported to a file format other than an R data format, _e.g._, .Rdata, .rda or .rds.

## Interrogating the API for available parameters

The POWER API offers functionality to get detailed information on any parameter offered or all parameters that are offered for a given community and temporal API.
This can be used to find available parameter names and definitions for each community and temporal API.

Fetch the complete available information for the temperature at 2 metres above the Earth's surface, T2M.


```r
query_parameters(par = "T2M")
#> $T2M
#> $T2M$temporal
#> $T2M$temporal$HOURLY
#> $T2M$temporal$HOURLY$name
#> [1] "Temperature at 2 Meters"
#> 
#> $T2M$temporal$HOURLY$definition
#> [1] "The average air (dry bulb) temperature at 2 meters above the surface of the earth."
#> 
#> $T2M$temporal$HOURLY$communities
#> [1] "AG" "RE" "SB"
#> 
#> $T2M$temporal$HOURLY$calculated
#> [1] FALSE
#> 
#> 
#> $T2M$temporal$DAILY
#> $T2M$temporal$DAILY$name
#> [1] "Temperature at 2 Meters"
#> 
#> $T2M$temporal$DAILY$definition
#> [1] "The average air (dry bulb) temperature at 2 meters above the surface of the earth."
#> 
#> $T2M$temporal$DAILY$communities
#> [1] "AG" "RE" "SB"
#> 
#> $T2M$temporal$DAILY$calculated
#> [1] FALSE
#> 
#> 
#> $T2M$temporal$MONTHLY
#> $T2M$temporal$MONTHLY$name
#> [1] "Temperature at 2 Meters"
#> 
#> $T2M$temporal$MONTHLY$definition
#> [1] "The average air (dry bulb) temperature at 2 meters above the surface of the earth."
#> 
#> $T2M$temporal$MONTHLY$communities
#> [1] "AG" "RE" "SB"
#> 
#> $T2M$temporal$MONTHLY$calculated
#> [1] FALSE
#> 
#> 
#> $T2M$temporal$CLIMATOLOGY
#> $T2M$temporal$CLIMATOLOGY$name
#> [1] "Temperature at 2 Meters"
#> 
#> $T2M$temporal$CLIMATOLOGY$definition
#> [1] "The average air (dry bulb) temperature at 2 meters above the surface of the earth."
#> 
#> $T2M$temporal$CLIMATOLOGY$communities
#> [1] "AG" "RE" "SB"
#> 
#> $T2M$temporal$CLIMATOLOGY$calculated
#> [1] FALSE
#> 
#> 
#> 
#> $T2M$type
#> [1] "METEOROLOGY"
```

Fetch complete temporal and community specific attribute information for "T2M" in the "ag" community for the "hourly" temporal API.


```r
query_parameters(par = "T2M",
                 community = "ag",
                 temporal_api = "hourly")
#> $T2M
#> $T2M$type
#> [1] "METEOROLOGY"
#> 
#> $T2M$temporal
#> [1] "HOURLY"
#> 
#> $T2M$source
#> [1] "MERRA2"
#> 
#> $T2M$community
#> [1] "AG"
#> 
#> $T2M$calculated
#> [1] FALSE
#> 
#> $T2M$inputs
#> NULL
#> 
#> $T2M$units
#> [1] "C"
#> 
#> $T2M$name
#> [1] "Temperature at 2 Meters"
#> 
#> $T2M$definition
#> [1] "The average air (dry bulb) temperature at 2 meters above the surface of the earth."
```

Fetch complete temporal and community specific attribute information for all parameters in the "ag" community for the "hourly" temporal API.


```r
query_parameters(community = "ag",
                 temporal_api = "hourly")
#> $PRECSNOLAND
#> $PRECSNOLAND$type
#> [1] "METEOROLOGY"
#> 
#> $PRECSNOLAND$temporal
#> [1] "HOURLY"
#> 
#> $PRECSNOLAND$source
#> [1] "MERRA2"
#> 
#> $PRECSNOLAND$community
#> [1] "AG"
#> 
#> $PRECSNOLAND$calculated
#> [1] FALSE
#> 
#> $PRECSNOLAND$inputs
#> NULL
#> 
#> $PRECSNOLAND$units
#> [1] "mm/hour"
#> 
#> $PRECSNOLAND$name
#> [1] "Snow Precipitation Land"
#> 
#> $PRECSNOLAND$definition
#> [1] "The snow precipitation only over land at the surface of the earth."
#> 
#> 
#> $PRECTOTCORR
#> $PRECTOTCORR$type
#> [1] "METEOROLOGY"
#> 
#> $PRECTOTCORR$temporal
#> [1] "HOURLY"
#> 
#> $PRECTOTCORR$source
#> [1] "MERRA2"
#> 
#> $PRECTOTCORR$community
#> [1] "AG"
#> 
#> $PRECTOTCORR$calculated
#> [1] FALSE
#> 
#> $PRECTOTCORR$inputs
#> NULL
#> 
#> $PRECTOTCORR$units
#> [1] "mm/hour"
#> 
#> $PRECTOTCORR$name
#> [1] "Precipitation Corrected"
#> 
#> $PRECTOTCORR$definition
#> [1] "The bias corrected average of total precipitation at the surface of the earth in water mass (includes water content in snow)."
#> 
#> 
#> $PS
#> $PS$type
#> [1] "METEOROLOGY"
#> 
#> $PS$temporal
#> [1] "HOURLY"
#> 
#> $PS$source
#> [1] "MERRA2"
#> 
#> $PS$community
#> [1] "AG"
#> 
#> $PS$calculated
#> [1] FALSE
#> 
#> $PS$inputs
#> NULL
#> 
#> $PS$units
#> [1] "kPa"
#> 
#> $PS$name
#> [1] "Surface Pressure"
#> 
#> $PS$definition
#> [1] "The average of surface pressure at the surface of the earth."
#> 
#> 
#> $QV10M
#> $QV10M$type
#> [1] "METEOROLOGY"
#> 
#> $QV10M$temporal
#> [1] "HOURLY"
#> 
#> $QV10M$source
#> [1] "MERRA2"
#> 
#> $QV10M$community
#> [1] "AG"
#> 
#> $QV10M$calculated
#> [1] FALSE
#> 
#> $QV10M$inputs
#> NULL
#> 
#> $QV10M$units
#> [1] "g/kg"
#> 
#> $QV10M$name
#> [1] "Specific Humidity at 10 Meters"
#> 
#> $QV10M$definition
#> [1] "The ratio of the mass of water vapor to the total mass of air at 10 meters (kg water/kg total air)."
#> 
#> 
#> $QV2M
#> $QV2M$type
#> [1] "METEOROLOGY"
#> 
#> $QV2M$temporal
#> [1] "HOURLY"
#> 
#> $QV2M$source
#> [1] "MERRA2"
#> 
#> $QV2M$community
#> [1] "AG"
#> 
#> $QV2M$calculated
#> [1] FALSE
#> 
#> $QV2M$inputs
#> NULL
#> 
#> $QV2M$units
#> [1] "g/kg"
#> 
#> $QV2M$name
#> [1] "Specific Humidity at 2 Meters"
#> 
#> $QV2M$definition
#> [1] "The ratio of the mass of water vapor to the total mass of air at 2 meters (kg water/kg total air)."
#> 
#> 
#> $RH2M
#> $RH2M$type
#> [1] "METEOROLOGY"
#> 
#> $RH2M$temporal
#> [1] "HOURLY"
#> 
#> $RH2M$source
#> [1] "MERRA2"
#> 
#> $RH2M$community
#> [1] "AG"
#> 
#> $RH2M$calculated
#> [1] FALSE
#> 
#> $RH2M$inputs
#> [1] "T2M"  "PS"   "QV2M"
#> 
#> $RH2M$units
#> [1] "%"
#> 
#> $RH2M$name
#> [1] "Relative Humidity at 2 Meters"
#> 
#> $RH2M$definition
#> [1] "The ratio of actual partial pressure of water vapor to the partial pressure at saturation, expressed in percent."
#> 
#> 
#> $SNODP
#> $SNODP$type
#> [1] "METEOROLOGY"
#> 
#> $SNODP$temporal
#> [1] "HOURLY"
#> 
#> $SNODP$source
#> [1] "MERRA2"
#> 
#> $SNODP$community
#> [1] "AG"
#> 
#> $SNODP$calculated
#> [1] FALSE
#> 
#> $SNODP$inputs
#> NULL
#> 
#> $SNODP$units
#> [1] "cm"
#> 
#> $SNODP$name
#> [1] "Snow Depth"
#> 
#> $SNODP$definition
#> [1] "The snow depth on land at surface of the earth."
#> 
#> 
#> $T2M
#> $T2M$type
#> [1] "METEOROLOGY"
#> 
#> $T2M$temporal
#> [1] "HOURLY"
#> 
#> $T2M$source
#> [1] "MERRA2"
#> 
#> $T2M$community
#> [1] "AG"
#> 
#> $T2M$calculated
#> [1] FALSE
#> 
#> $T2M$inputs
#> NULL
#> 
#> $T2M$units
#> [1] "C"
#> 
#> $T2M$name
#> [1] "Temperature at 2 Meters"
#> 
#> $T2M$definition
#> [1] "The average air (dry bulb) temperature at 2 meters above the surface of the earth."
#> 
#> 
#> $TS
#> $TS$type
#> [1] "METEOROLOGY"
#> 
#> $TS$temporal
#> [1] "HOURLY"
#> 
#> $TS$source
#> [1] "MERRA2"
#> 
#> $TS$community
#> [1] "AG"
#> 
#> $TS$calculated
#> [1] FALSE
#> 
#> $TS$inputs
#> NULL
#> 
#> $TS$units
#> [1] "C"
#> 
#> $TS$name
#> [1] "Earth Skin Temperature"
#> 
#> $TS$definition
#> [1] "The average temperature at the earth's surface."
#> 
#> 
#> $U10M
#> $U10M$type
#> [1] "METEOROLOGY"
#> 
#> $U10M$temporal
#> [1] "HOURLY"
#> 
#> $U10M$source
#> [1] "MERRA2"
#> 
#> $U10M$community
#> [1] "AG"
#> 
#> $U10M$calculated
#> [1] FALSE
#> 
#> $U10M$inputs
#> NULL
#> 
#> $U10M$units
#> [1] "m/s"
#> 
#> $U10M$name
#> [1] "Eastward Wind at 10 Meters"
#> 
#> $U10M$definition
#> [1] "The estimate of the eastward wind average speed for winds blowing 10 meters above the surface of the earth."
#> 
#> 
#> $U2M
#> $U2M$type
#> [1] "METEOROLOGY"
#> 
#> $U2M$temporal
#> [1] "HOURLY"
#> 
#> $U2M$source
#> [1] "MERRA2"
#> 
#> $U2M$community
#> [1] "AG"
#> 
#> $U2M$calculated
#> [1] FALSE
#> 
#> $U2M$inputs
#> NULL
#> 
#> $U2M$units
#> [1] "m/s"
#> 
#> $U2M$name
#> [1] "Eastward Wind at 2 Meters"
#> 
#> $U2M$definition
#> [1] "The estimate of the eastward wind average speed for winds blowing 2 meters above the surface of the earth."
#> 
#> 
#> $U50M
#> $U50M$type
#> [1] "METEOROLOGY"
#> 
#> $U50M$temporal
#> [1] "HOURLY"
#> 
#> $U50M$source
#> [1] "MERRA2"
#> 
#> $U50M$community
#> [1] "AG"
#> 
#> $U50M$calculated
#> [1] FALSE
#> 
#> $U50M$inputs
#> NULL
#> 
#> $U50M$units
#> [1] "m/s"
#> 
#> $U50M$name
#> [1] "Eastward Wind at 50 Meters"
#> 
#> $U50M$definition
#> [1] "The estimate of the eastward wind average speed for winds blowing 50 meters above the surface of the earth."
#> 
#> 
#> $V10M
#> $V10M$type
#> [1] "METEOROLOGY"
#> 
#> $V10M$temporal
#> [1] "HOURLY"
#> 
#> $V10M$source
#> [1] "MERRA2"
#> 
#> $V10M$community
#> [1] "AG"
#> 
#> $V10M$calculated
#> [1] FALSE
#> 
#> $V10M$inputs
#> NULL
#> 
#> $V10M$units
#> [1] "m/s"
#> 
#> $V10M$name
#> [1] "Northward Wind at 10 Meters"
#> 
#> $V10M$definition
#> [1] "The estimate of the northward wind average speed for winds blowing 10 meters above the surface of the earth."
#> 
#> 
#> $V2M
#> $V2M$type
#> [1] "METEOROLOGY"
#> 
#> $V2M$temporal
#> [1] "HOURLY"
#> 
#> $V2M$source
#> [1] "MERRA2"
#> 
#> $V2M$community
#> [1] "AG"
#> 
#> $V2M$calculated
#> [1] FALSE
#> 
#> $V2M$inputs
#> NULL
#> 
#> $V2M$units
#> [1] "m/s"
#> 
#> $V2M$name
#> [1] "Northward Wind at 2 Meters"
#> 
#> $V2M$definition
#> [1] "The estimate of the northward wind average speed for winds blowing 2 meters above the surface of the earth."
#> 
#> 
#> $V50M
#> $V50M$type
#> [1] "METEOROLOGY"
#> 
#> $V50M$temporal
#> [1] "HOURLY"
#> 
#> $V50M$source
#> [1] "MERRA2"
#> 
#> $V50M$community
#> [1] "AG"
#> 
#> $V50M$calculated
#> [1] FALSE
#> 
#> $V50M$inputs
#> NULL
#> 
#> $V50M$units
#> [1] "m/s"
#> 
#> $V50M$name
#> [1] "Northward Wind at 50 Meters"
#> 
#> $V50M$definition
#> [1] "The estimate of the northward wind average speed for winds blowing 50 meters above the surface of the earth."
#> 
#> 
#> $PSC
#> $PSC$type
#> [1] "METEOROLOGY"
#> 
#> $PSC$temporal
#> [1] "HOURLY"
#> 
#> $PSC$source
#> [1] "POWER"
#> 
#> $PSC$community
#> [1] "AG"
#> 
#> $PSC$calculated
#> [1] TRUE
#> 
#> $PSC$inputs
#> [1] "PS"  "T2M"
#> 
#> $PSC$units
#> [1] "kPa"
#> 
#> $PSC$name
#> [1] "Corrected Atmospheric Pressure (Adjusted For Site Elevation)"
#> 
#> $PSC$definition
#> [1] "Atmospheric pressure associated with the MERRA-2 grid has been adjusted based upon the difference between the elevation of an underlying surface site and the average elevation of the MERRA-2 grid cell."
#> 
#> 
#> $T2MDEW
#> $T2MDEW$type
#> [1] "METEOROLOGY"
#> 
#> $T2MDEW$temporal
#> [1] "HOURLY"
#> 
#> $T2MDEW$source
#> [1] "POWER"
#> 
#> $T2MDEW$community
#> [1] "AG"
#> 
#> $T2MDEW$calculated
#> [1] FALSE
#> 
#> $T2MDEW$inputs
#> [1] "T2M"  "RH2M"
#> 
#> $T2MDEW$units
#> [1] "C"
#> 
#> $T2MDEW$name
#> [1] "Dew/Frost Point at 2 Meters"
#> 
#> $T2MDEW$definition
#> [1] "The dew/frost point temperature at 2 meters above the surface of the earth."
#> 
#> 
#> $T2MWET
#> $T2MWET$type
#> [1] "METEOROLOGY"
#> 
#> $T2MWET$temporal
#> [1] "HOURLY"
#> 
#> $T2MWET$source
#> [1] "POWER"
#> 
#> $T2MWET$community
#> [1] "AG"
#> 
#> $T2MWET$calculated
#> [1] FALSE
#> 
#> $T2MWET$inputs
#> [1] "PS"     "T2M"    "T2MDEW"
#> 
#> $T2MWET$units
#> [1] "C"
#> 
#> $T2MWET$name
#> [1] "Wet Bulb Temperature at 2 Meters"
#> 
#> $T2MWET$definition
#> [1] "The adiabatic saturation temperature which can be measured by a thermometer covered in a water-soaked cloth over which air is passed at 2 meters above the surface of the earth."
#> 
#> 
#> $WD10M
#> $WD10M$type
#> [1] "METEOROLOGY"
#> 
#> $WD10M$temporal
#> [1] "HOURLY"
#> 
#> $WD10M$source
#> [1] "POWER"
#> 
#> $WD10M$community
#> [1] "AG"
#> 
#> $WD10M$calculated
#> [1] TRUE
#> 
#> $WD10M$inputs
#> [1] "U10M" "V10M"
#> 
#> $WD10M$units
#> [1] "Degrees"
#> 
#> $WD10M$name
#> [1] "Wind Direction at 10 Meters"
#> 
#> $WD10M$definition
#> [1] "The average of the wind direction at 10 meters above the surface of the earth."
#> 
#> 
#> $WD2M
#> $WD2M$type
#> [1] "METEOROLOGY"
#> 
#> $WD2M$temporal
#> [1] "HOURLY"
#> 
#> $WD2M$source
#> [1] "POWER"
#> 
#> $WD2M$community
#> [1] "AG"
#> 
#> $WD2M$calculated
#> [1] TRUE
#> 
#> $WD2M$inputs
#> [1] "U2M" "V2M"
#> 
#> $WD2M$units
#> [1] "Degrees"
#> 
#> $WD2M$name
#> [1] "Wind Direction at 2 Meters"
#> 
#> $WD2M$definition
#> [1] "The average of the wind direction at 2 meters above the surface of the earth."
#> 
#> 
#> $WD50M
#> $WD50M$type
#> [1] "METEOROLOGY"
#> 
#> $WD50M$temporal
#> [1] "HOURLY"
#> 
#> $WD50M$source
#> [1] "POWER"
#> 
#> $WD50M$community
#> [1] "AG"
#> 
#> $WD50M$calculated
#> [1] TRUE
#> 
#> $WD50M$inputs
#> [1] "U50M" "V50M"
#> 
#> $WD50M$units
#> [1] "Degrees"
#> 
#> $WD50M$name
#> [1] "Wind Direction at 50 Meters"
#> 
#> $WD50M$definition
#> [1] "The average of the wind direction at 50 meters above the surface of the earth."
#> 
#> 
#> $WS10M
#> $WS10M$type
#> [1] "METEOROLOGY"
#> 
#> $WS10M$temporal
#> [1] "HOURLY"
#> 
#> $WS10M$source
#> [1] "POWER"
#> 
#> $WS10M$community
#> [1] "AG"
#> 
#> $WS10M$calculated
#> [1] TRUE
#> 
#> $WS10M$inputs
#> [1] "U10M" "V10M"
#> 
#> $WS10M$units
#> [1] "m/s"
#> 
#> $WS10M$name
#> [1] "Wind Speed at 10 Meters"
#> 
#> $WS10M$definition
#> [1] "The average of wind speed at 10 meters above the surface of the earth."
#> 
#> 
#> $WS2M
#> $WS2M$type
#> [1] "METEOROLOGY"
#> 
#> $WS2M$temporal
#> [1] "HOURLY"
#> 
#> $WS2M$source
#> [1] "POWER"
#> 
#> $WS2M$community
#> [1] "AG"
#> 
#> $WS2M$calculated
#> [1] TRUE
#> 
#> $WS2M$inputs
#> [1] "U2M" "V2M"
#> 
#> $WS2M$units
#> [1] "m/s"
#> 
#> $WS2M$name
#> [1] "Wind Speed at 2 Meters"
#> 
#> $WS2M$definition
#> [1] "The average of wind speed at 2 meters above the surface of the earth."
#> 
#> 
#> $WS50M
#> $WS50M$type
#> [1] "METEOROLOGY"
#> 
#> $WS50M$temporal
#> [1] "HOURLY"
#> 
#> $WS50M$source
#> [1] "POWER"
#> 
#> $WS50M$community
#> [1] "AG"
#> 
#> $WS50M$calculated
#> [1] TRUE
#> 
#> $WS50M$inputs
#> [1] "U50M" "V50M"
#> 
#> $WS50M$units
#> [1] "m/s"
#> 
#> $WS50M$name
#> [1] "Wind Speed at 50 Meters"
#> 
#> $WS50M$definition
#> [1] "The average of wind speed at 50 meters above the surface of the earth."
#> 
#> 
#> $WSC
#> $WSC$type
#> [1] "METEOROLOGY"
#> 
#> $WSC$temporal
#> [1] "HOURLY"
#> 
#> $WSC$source
#> [1] "POWER"
#> 
#> $WSC$community
#> [1] "AG"
#> 
#> $WSC$calculated
#> [1] TRUE
#> 
#> $WSC$inputs
#> [1] "WS10M" "WS50M"
#> 
#> $WSC$units
#> [1] "m/s"
#> 
#> $WSC$name
#> [1] "Corrected Wind Speed (Adjusted For Elevation)"
#> 
#> $WSC$definition
#> [1] "Wind speed associated with the MERRA-2 grid adjusted using the Gipe Power Law and the elevation difference between the average elevation of the grid cell and the elevation of the underlying surface."
#> 
#> 
#> $ALLSKY_SFC_LW_DWN
#> $ALLSKY_SFC_LW_DWN$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_LW_DWN$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_LW_DWN$source
#> [1] "CERES"
#> 
#> $ALLSKY_SFC_LW_DWN$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_LW_DWN$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_LW_DWN$inputs
#> NULL
#> 
#> $ALLSKY_SFC_LW_DWN$units
#> [1] "W/m^2"
#> 
#> $ALLSKY_SFC_LW_DWN$name
#> [1] "All Sky Surface Longwave Downward Irradiance"
#> 
#> $ALLSKY_SFC_LW_DWN$definition
#> [1] "The downward thermal infrared irradiance under all sky conditions reaching a horizontal plane the surface of the earth. Also known as Horizontal Infrared Radiation Intensity from Sky."
#> 
#> 
#> $ALLSKY_SFC_SW_DIFF
#> $ALLSKY_SFC_SW_DIFF$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_SW_DIFF$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_SW_DIFF$source
#> [1] "CERES"
#> 
#> $ALLSKY_SFC_SW_DIFF$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_SW_DIFF$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_SW_DIFF$inputs
#> NULL
#> 
#> $ALLSKY_SFC_SW_DIFF$units
#> [1] "MJ/hr"
#> 
#> $ALLSKY_SFC_SW_DIFF$name
#> [1] "All Sky Surface Shortwave Diffuse Irradiance"
#> 
#> $ALLSKY_SFC_SW_DIFF$definition
#> [1] "The diffuse (light energy scattered out of the direction of the sun) solar irradiance incident on a horizontal plane at the surface of the earth under all sky conditions."
#> 
#> 
#> $ALLSKY_SFC_SW_DWN
#> $ALLSKY_SFC_SW_DWN$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_SW_DWN$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_SW_DWN$source
#> [1] "CERES"
#> 
#> $ALLSKY_SFC_SW_DWN$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_SW_DWN$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_SW_DWN$inputs
#> NULL
#> 
#> $ALLSKY_SFC_SW_DWN$units
#> [1] "MJ/hr"
#> 
#> $ALLSKY_SFC_SW_DWN$name
#> [1] "All Sky Surface Shortwave Downward Irradiance"
#> 
#> $ALLSKY_SFC_SW_DWN$definition
#> [1] "The total solar irradiance incident (direct plus diffuse) on a horizontal plane at the surface of the earth under all sky conditions. An alternative term for the total solar irradiance is the \"Global Horizontal Irradiance\" or GHI."
#> 
#> 
#> $ALLSKY_SFC_UV_INDEX
#> $ALLSKY_SFC_UV_INDEX$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_UV_INDEX$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_UV_INDEX$source
#> [1] "CERES"
#> 
#> $ALLSKY_SFC_UV_INDEX$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_UV_INDEX$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_UV_INDEX$inputs
#> NULL
#> 
#> $ALLSKY_SFC_UV_INDEX$units
#> [1] "dimensionless"
#> 
#> $ALLSKY_SFC_UV_INDEX$name
#> [1] "All Sky Surface UV Index"
#> 
#> $ALLSKY_SFC_UV_INDEX$definition
#> [1] "The ultraviolet radiation exposure index."
#> 
#> 
#> $ALLSKY_SFC_UVA
#> $ALLSKY_SFC_UVA$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_UVA$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_UVA$source
#> [1] "CERES"
#> 
#> $ALLSKY_SFC_UVA$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_UVA$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_UVA$inputs
#> NULL
#> 
#> $ALLSKY_SFC_UVA$units
#> [1] "W/m^2"
#> 
#> $ALLSKY_SFC_UVA$name
#> [1] "All Sky Surface UVA Irradiance"
#> 
#> $ALLSKY_SFC_UVA$definition
#> [1] "The ultraviolet A (UVA 315nm-400nm) irradiance under all sky conditions."
#> 
#> 
#> $ALLSKY_SFC_UVB
#> $ALLSKY_SFC_UVB$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_UVB$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_UVB$source
#> [1] "CERES"
#> 
#> $ALLSKY_SFC_UVB$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_UVB$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_UVB$inputs
#> NULL
#> 
#> $ALLSKY_SFC_UVB$units
#> [1] "W/m^2"
#> 
#> $ALLSKY_SFC_UVB$name
#> [1] "All Sky Surface UVB Irradiance"
#> 
#> $ALLSKY_SFC_UVB$definition
#> [1] "The ultraviolet B (UVB 280nm-315nm) irradiance under all sky conditions."
#> 
#> 
#> $AOD_55
#> $AOD_55$type
#> [1] "RADIATION"
#> 
#> $AOD_55$temporal
#> [1] "HOURLY"
#> 
#> $AOD_55$source
#> [1] "CERES"
#> 
#> $AOD_55$community
#> [1] "AG"
#> 
#> $AOD_55$calculated
#> [1] FALSE
#> 
#> $AOD_55$inputs
#> NULL
#> 
#> $AOD_55$units
#> [1] "dimensionless"
#> 
#> $AOD_55$name
#> [1] "Aerosol Optical Depth 55"
#> 
#> $AOD_55$definition
#> [1] "The optical thickness at 0.55 um measured vertically; the component of the atmosphere to quantify the removal of radiant energy from an incident beam."
#> 
#> 
#> $AOD_84
#> $AOD_84$type
#> [1] "RADIATION"
#> 
#> $AOD_84$temporal
#> [1] "HOURLY"
#> 
#> $AOD_84$source
#> [1] "CERES"
#> 
#> $AOD_84$community
#> [1] "AG"
#> 
#> $AOD_84$calculated
#> [1] FALSE
#> 
#> $AOD_84$inputs
#> NULL
#> 
#> $AOD_84$units
#> [1] "dimensionless"
#> 
#> $AOD_84$name
#> [1] "Aerosol Optical Depth 84"
#> 
#> $AOD_84$definition
#> [1] "The optical thickness at 0.84 um measured vertically; the component of the atmosphere to quantify the removal of radiant energy from an incident beam."
#> 
#> 
#> $CLOUD_AMT
#> $CLOUD_AMT$type
#> [1] "RADIATION"
#> 
#> $CLOUD_AMT$temporal
#> [1] "HOURLY"
#> 
#> $CLOUD_AMT$source
#> [1] "CERES"
#> 
#> $CLOUD_AMT$community
#> [1] "AG"
#> 
#> $CLOUD_AMT$calculated
#> [1] FALSE
#> 
#> $CLOUD_AMT$inputs
#> NULL
#> 
#> $CLOUD_AMT$units
#> [1] "%"
#> 
#> $CLOUD_AMT$name
#> [1] "Cloud Amount"
#> 
#> $CLOUD_AMT$definition
#> [1] "The average percent of cloud amount during the temporal period."
#> 
#> 
#> $CLOUD_OD
#> $CLOUD_OD$type
#> [1] "RADIATION"
#> 
#> $CLOUD_OD$temporal
#> [1] "HOURLY"
#> 
#> $CLOUD_OD$source
#> [1] "CERES"
#> 
#> $CLOUD_OD$community
#> [1] "AG"
#> 
#> $CLOUD_OD$calculated
#> [1] FALSE
#> 
#> $CLOUD_OD$inputs
#> NULL
#> 
#> $CLOUD_OD$units
#> [1] "dimensionless"
#> 
#> $CLOUD_OD$name
#> [1] "Cloud Optical Visible Depth"
#> 
#> $CLOUD_OD$definition
#> [1] "The vertical optical thickness between the top and bottom of a cloud."
#> 
#> 
#> $CLRSKY_SFC_LW_DWN
#> $CLRSKY_SFC_LW_DWN$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_SFC_LW_DWN$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_SFC_LW_DWN$source
#> [1] "CERES"
#> 
#> $CLRSKY_SFC_LW_DWN$community
#> [1] "AG"
#> 
#> $CLRSKY_SFC_LW_DWN$calculated
#> [1] FALSE
#> 
#> $CLRSKY_SFC_LW_DWN$inputs
#> NULL
#> 
#> $CLRSKY_SFC_LW_DWN$units
#> [1] "W/m^2"
#> 
#> $CLRSKY_SFC_LW_DWN$name
#> [1] "Clear Sky Surface Longwave Downward Irradiance"
#> 
#> $CLRSKY_SFC_LW_DWN$definition
#> [1] "The downward thermal infrared irradiance under clear sky conditions reaching a horizontal plane the surface of the earth. Also known as Horizontal Infrared Radiation Intensity from Sky."
#> 
#> 
#> $CLRSKY_SFC_SW_DIFF
#> $CLRSKY_SFC_SW_DIFF$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_SFC_SW_DIFF$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_SFC_SW_DIFF$source
#> [1] "CERES"
#> 
#> $CLRSKY_SFC_SW_DIFF$community
#> [1] "AG"
#> 
#> $CLRSKY_SFC_SW_DIFF$calculated
#> [1] FALSE
#> 
#> $CLRSKY_SFC_SW_DIFF$inputs
#> NULL
#> 
#> $CLRSKY_SFC_SW_DIFF$units
#> [1] "MJ/hr"
#> 
#> $CLRSKY_SFC_SW_DIFF$name
#> [1] "Clear Sky Surface Shortwave Downward Diffuse Horizontal Irradiance"
#> 
#> $CLRSKY_SFC_SW_DIFF$definition
#> [1] "The diffuse (light energy scattered out of the direction of the sun) solar irradiance incident on a horizontal plane at the surface of the earth under clear sky conditions."
#> 
#> 
#> $CLRSKY_SFC_SW_DIRH
#> $CLRSKY_SFC_SW_DIRH$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_SFC_SW_DIRH$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_SFC_SW_DIRH$source
#> [1] "CERES"
#> 
#> $CLRSKY_SFC_SW_DIRH$community
#> [1] "AG"
#> 
#> $CLRSKY_SFC_SW_DIRH$calculated
#> [1] FALSE
#> 
#> $CLRSKY_SFC_SW_DIRH$inputs
#> NULL
#> 
#> $CLRSKY_SFC_SW_DIRH$units
#> [1] "MJ/m^2/day"
#> 
#> $CLRSKY_SFC_SW_DIRH$name
#> [1] "Clear Sky Surface Shortwave Direct Horizontal Irradiance"
#> 
#> $CLRSKY_SFC_SW_DIRH$definition
#> [1] "The direct solar irradiance incident on a horizontal plane at the surface of the earth under clear sky conditions."
#> 
#> 
#> $CLRSKY_SFC_SW_DWN
#> $CLRSKY_SFC_SW_DWN$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_SFC_SW_DWN$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_SFC_SW_DWN$source
#> [1] "CERES"
#> 
#> $CLRSKY_SFC_SW_DWN$community
#> [1] "AG"
#> 
#> $CLRSKY_SFC_SW_DWN$calculated
#> [1] FALSE
#> 
#> $CLRSKY_SFC_SW_DWN$inputs
#> NULL
#> 
#> $CLRSKY_SFC_SW_DWN$units
#> [1] "MJ/hr"
#> 
#> $CLRSKY_SFC_SW_DWN$name
#> [1] "Clear Sky Surface Shortwave Downward Irradiance"
#> 
#> $CLRSKY_SFC_SW_DWN$definition
#> [1] "The total solar irradiance incident (direct plus diffuse) on a horizontal plane at the surface of the earth under clear sky conditions. An alternative term for the total solar irradiance is the \"Global Horizontal Irradiance\" or GHI."
#> 
#> 
#> $PW
#> $PW$type
#> [1] "RADIATION"
#> 
#> $PW$temporal
#> [1] "HOURLY"
#> 
#> $PW$source
#> [1] "CERES"
#> 
#> $PW$community
#> [1] "AG"
#> 
#> $PW$calculated
#> [1] FALSE
#> 
#> $PW$inputs
#> NULL
#> 
#> $PW$units
#> [1] "cm"
#> 
#> $PW$name
#> [1] "Precipitable Water"
#> 
#> $PW$definition
#> [1] "The total atmospheric water vapor contained in a vertical column of the atmosphere."
#> 
#> 
#> $SZA
#> $SZA$type
#> [1] "RADIATION"
#> 
#> $SZA$temporal
#> [1] "HOURLY"
#> 
#> $SZA$source
#> [1] "CERES"
#> 
#> $SZA$community
#> [1] "AG"
#> 
#> $SZA$calculated
#> [1] FALSE
#> 
#> $SZA$inputs
#> NULL
#> 
#> $SZA$units
#> [1] "Degrees"
#> 
#> $SZA$name
#> [1] "Solar Zenith Angle"
#> 
#> $SZA$definition
#> [1] "The angle between the geodetic zenith vector and a vector from the earth point to the sun integrated over the period."
#> 
#> 
#> $TOA_SW_DWN
#> $TOA_SW_DWN$type
#> [1] "RADIATION"
#> 
#> $TOA_SW_DWN$temporal
#> [1] "HOURLY"
#> 
#> $TOA_SW_DWN$source
#> [1] "CERES"
#> 
#> $TOA_SW_DWN$community
#> [1] "AG"
#> 
#> $TOA_SW_DWN$calculated
#> [1] FALSE
#> 
#> $TOA_SW_DWN$inputs
#> NULL
#> 
#> $TOA_SW_DWN$units
#> [1] "MJ/hr"
#> 
#> $TOA_SW_DWN$name
#> [1] "Top-Of-Atmosphere Shortwave Downward Irradiance"
#> 
#> $TOA_SW_DWN$definition
#> [1] "The total solar irradiance incident (direct plus diffuse) on a horizontal plane at the top of the atmosphere (extraterrestrial radiation)."
#> 
#> 
#> $ALLSKY_KT
#> $ALLSKY_KT$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_KT$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_KT$source
#> [1] "POWER"
#> 
#> $ALLSKY_KT$community
#> [1] "AG"
#> 
#> $ALLSKY_KT$calculated
#> [1] FALSE
#> 
#> $ALLSKY_KT$inputs
#> NULL
#> 
#> $ALLSKY_KT$units
#> [1] "dimensionless"
#> 
#> $ALLSKY_KT$name
#> [1] "All Sky Insolation Clearness Index"
#> 
#> $ALLSKY_KT$definition
#> [1] "A fraction representing clearness of the atmosphere; the all sky insolation that is transmitted through the atmosphere to strike the surface of the earth divided by the average of top of the atmosphere total solar irradiance incident."
#> 
#> 
#> $ALLSKY_NKT
#> $ALLSKY_NKT$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_NKT$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_NKT$source
#> [1] "POWER"
#> 
#> $ALLSKY_NKT$community
#> [1] "AG"
#> 
#> $ALLSKY_NKT$calculated
#> [1] FALSE
#> 
#> $ALLSKY_NKT$inputs
#> NULL
#> 
#> $ALLSKY_NKT$units
#> [1] "dimensionless"
#> 
#> $ALLSKY_NKT$name
#> [1] "All Sky Normalized Insolation Clearness Index"
#> 
#> $ALLSKY_NKT$definition
#> [1] "The average zenith angle-independent expression of the all sky insolation clearness index."
#> 
#> 
#> $ALLSKY_SFC_PAR_TOT
#> $ALLSKY_SFC_PAR_TOT$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_PAR_TOT$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_PAR_TOT$source
#> [1] "POWER"
#> 
#> $ALLSKY_SFC_PAR_TOT$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_PAR_TOT$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_PAR_TOT$inputs
#> NULL
#> 
#> $ALLSKY_SFC_PAR_TOT$units
#> [1] "W/m^2"
#> 
#> $ALLSKY_SFC_PAR_TOT$name
#> [1] "All Sky Surface PAR Total"
#> 
#> $ALLSKY_SFC_PAR_TOT$definition
#> [1] "The total Photosynthetically Active Radiation (PAR) incident on a horizontal plane at the surface of the earth under all sky conditions."
#> 
#> 
#> $ALLSKY_SFC_SW_DNI
#> $ALLSKY_SFC_SW_DNI$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SFC_SW_DNI$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SFC_SW_DNI$source
#> [1] "POWER"
#> 
#> $ALLSKY_SFC_SW_DNI$community
#> [1] "AG"
#> 
#> $ALLSKY_SFC_SW_DNI$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SFC_SW_DNI$inputs
#> NULL
#> 
#> $ALLSKY_SFC_SW_DNI$units
#> [1] "MJ/hr"
#> 
#> $ALLSKY_SFC_SW_DNI$name
#> [1] "All Sky Surface Shortwave Downward Direct Normal Irradiance"
#> 
#> $ALLSKY_SFC_SW_DNI$definition
#> [1] "The direct solar irradiance incident to a horizontal plane normal (perpendicular) to the direction of the sun's position under all sky conditions."
#> 
#> 
#> $ALLSKY_SRF_ALB
#> $ALLSKY_SRF_ALB$type
#> [1] "RADIATION"
#> 
#> $ALLSKY_SRF_ALB$temporal
#> [1] "HOURLY"
#> 
#> $ALLSKY_SRF_ALB$source
#> [1] "POWER"
#> 
#> $ALLSKY_SRF_ALB$community
#> [1] "AG"
#> 
#> $ALLSKY_SRF_ALB$calculated
#> [1] FALSE
#> 
#> $ALLSKY_SRF_ALB$inputs
#> NULL
#> 
#> $ALLSKY_SRF_ALB$units
#> [1] "dimensionless"
#> 
#> $ALLSKY_SRF_ALB$name
#> [1] "All Sky Surface Albedo"
#> 
#> $ALLSKY_SRF_ALB$definition
#> [1] "The all sky rate of reflectivity of the earth's surface; the ratio of the solar energy reflected by the surface of the earth compared to the total solar energy incident reaching the surface of the earth."
#> 
#> 
#> $CLRSKY_KT
#> $CLRSKY_KT$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_KT$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_KT$source
#> [1] "POWER"
#> 
#> $CLRSKY_KT$community
#> [1] "AG"
#> 
#> $CLRSKY_KT$calculated
#> [1] FALSE
#> 
#> $CLRSKY_KT$inputs
#> NULL
#> 
#> $CLRSKY_KT$units
#> [1] "dimensionless"
#> 
#> $CLRSKY_KT$name
#> [1] "Clear Sky Insolation Clearness Index"
#> 
#> $CLRSKY_KT$definition
#> [1] "A fraction representing clearness of the atmosphere; the clear sky insolation that is transmitted through the atmosphere to strike the surface of the earth divided by the average of top of the atmosphere total solar irradiance incident."
#> 
#> 
#> $CLRSKY_NKT
#> $CLRSKY_NKT$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_NKT$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_NKT$source
#> [1] "POWER"
#> 
#> $CLRSKY_NKT$community
#> [1] "AG"
#> 
#> $CLRSKY_NKT$calculated
#> [1] FALSE
#> 
#> $CLRSKY_NKT$inputs
#> NULL
#> 
#> $CLRSKY_NKT$units
#> [1] "dimensionless"
#> 
#> $CLRSKY_NKT$name
#> [1] "Clear Sky Normalized Insolation Clearness Index"
#> 
#> $CLRSKY_NKT$definition
#> [1] "The average zenith angle-independent expression of the clear sky insolation clearness index."
#> 
#> 
#> $CLRSKY_SFC_PAR_TOT
#> $CLRSKY_SFC_PAR_TOT$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_SFC_PAR_TOT$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_SFC_PAR_TOT$source
#> [1] "POWER"
#> 
#> $CLRSKY_SFC_PAR_TOT$community
#> [1] "AG"
#> 
#> $CLRSKY_SFC_PAR_TOT$calculated
#> [1] FALSE
#> 
#> $CLRSKY_SFC_PAR_TOT$inputs
#> NULL
#> 
#> $CLRSKY_SFC_PAR_TOT$units
#> [1] "W/m^2"
#> 
#> $CLRSKY_SFC_PAR_TOT$name
#> [1] "Clear Sky Surface PAR Total"
#> 
#> $CLRSKY_SFC_PAR_TOT$definition
#> [1] "The total Photosynthetically Active Radiation (PAR) incident on a horizontal plane at the surface of the earth under clear sky conditions."
#> 
#> 
#> $CLRSKY_SFC_SW_DNI
#> $CLRSKY_SFC_SW_DNI$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_SFC_SW_DNI$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_SFC_SW_DNI$source
#> [1] "POWER"
#> 
#> $CLRSKY_SFC_SW_DNI$community
#> [1] "AG"
#> 
#> $CLRSKY_SFC_SW_DNI$calculated
#> [1] FALSE
#> 
#> $CLRSKY_SFC_SW_DNI$inputs
#> NULL
#> 
#> $CLRSKY_SFC_SW_DNI$units
#> [1] "MJ/hr"
#> 
#> $CLRSKY_SFC_SW_DNI$name
#> [1] "Clear Sky Surface Shortwave Downward Direct Normal Irradiance"
#> 
#> $CLRSKY_SFC_SW_DNI$definition
#> [1] "The direct solar irradiance incident to a horizontal plane normal (perpendicular) to the direction of the sun's position under clear sky conditions."
#> 
#> 
#> $CLRSKY_SRF_ALB
#> $CLRSKY_SRF_ALB$type
#> [1] "RADIATION"
#> 
#> $CLRSKY_SRF_ALB$temporal
#> [1] "HOURLY"
#> 
#> $CLRSKY_SRF_ALB$source
#> [1] "POWER"
#> 
#> $CLRSKY_SRF_ALB$community
#> [1] "AG"
#> 
#> $CLRSKY_SRF_ALB$calculated
#> [1] FALSE
#> 
#> $CLRSKY_SRF_ALB$inputs
#> NULL
#> 
#> $CLRSKY_SRF_ALB$units
#> [1] "dimensionless"
#> 
#> $CLRSKY_SRF_ALB$name
#> [1] "Clear Sky Surface Albedo"
#> 
#> $CLRSKY_SRF_ALB$definition
#> [1] "The clear sky rate of reflectivity of the earth's surface; the ratio of the solar energy reflected by the surface of the earth compared to the total solar energy incident reaching the surface of the earth."
#> 
#> 
#> $DIFFUSE_ILLUMINANCE
#> $DIFFUSE_ILLUMINANCE$type
#> [1] "RADIATION"
#> 
#> $DIFFUSE_ILLUMINANCE$temporal
#> [1] "HOURLY"
#> 
#> $DIFFUSE_ILLUMINANCE$source
#> [1] "POWER"
#> 
#> $DIFFUSE_ILLUMINANCE$community
#> [1] "AG"
#> 
#> $DIFFUSE_ILLUMINANCE$calculated
#> [1] FALSE
#> 
#> $DIFFUSE_ILLUMINANCE$inputs
#> NULL
#> 
#> $DIFFUSE_ILLUMINANCE$units
#> [1] "lux"
#> 
#> $DIFFUSE_ILLUMINANCE$name
#> [1] "Diffuse Illuminance"
#> 
#> $DIFFUSE_ILLUMINANCE$definition
#> [1] "The average amount of illuminance received directly from the solar disk on a surface perpendicular to the suns rays."
#> 
#> 
#> $DIRECT_ILLUMINANCE
#> $DIRECT_ILLUMINANCE$type
#> [1] "RADIATION"
#> 
#> $DIRECT_ILLUMINANCE$temporal
#> [1] "HOURLY"
#> 
#> $DIRECT_ILLUMINANCE$source
#> [1] "POWER"
#> 
#> $DIRECT_ILLUMINANCE$community
#> [1] "AG"
#> 
#> $DIRECT_ILLUMINANCE$calculated
#> [1] FALSE
#> 
#> $DIRECT_ILLUMINANCE$inputs
#> NULL
#> 
#> $DIRECT_ILLUMINANCE$units
#> [1] "lux"
#> 
#> $DIRECT_ILLUMINANCE$name
#> [1] "Direct Illuminance"
#> 
#> $DIRECT_ILLUMINANCE$definition
#> [1] "The average amount of illuminance received from the sky (excluding the solar disk) on a horizontal plane."
#> 
#> 
#> $GLOBAL_ILLUMINANCE
#> $GLOBAL_ILLUMINANCE$type
#> [1] "RADIATION"
#> 
#> $GLOBAL_ILLUMINANCE$temporal
#> [1] "HOURLY"
#> 
#> $GLOBAL_ILLUMINANCE$source
#> [1] "POWER"
#> 
#> $GLOBAL_ILLUMINANCE$community
#> [1] "AG"
#> 
#> $GLOBAL_ILLUMINANCE$calculated
#> [1] FALSE
#> 
#> $GLOBAL_ILLUMINANCE$inputs
#> NULL
#> 
#> $GLOBAL_ILLUMINANCE$units
#> [1] "lux"
#> 
#> $GLOBAL_ILLUMINANCE$name
#> [1] "Global Illuminance"
#> 
#> $GLOBAL_ILLUMINANCE$definition
#> [1] "The average total amount of direct and diffuse illuminance on a horizontal plane."
#> 
#> 
#> $TOA_SW_DNI
#> $TOA_SW_DNI$type
#> [1] "RADIATION"
#> 
#> $TOA_SW_DNI$temporal
#> [1] "HOURLY"
#> 
#> $TOA_SW_DNI$source
#> [1] "POWER"
#> 
#> $TOA_SW_DNI$community
#> [1] "AG"
#> 
#> $TOA_SW_DNI$calculated
#> [1] FALSE
#> 
#> $TOA_SW_DNI$inputs
#> NULL
#> 
#> $TOA_SW_DNI$units
#> [1] "MJ/hr"
#> 
#> $TOA_SW_DNI$name
#> [1] "Top-Of-Atmosphere Shortwave Direct Normal Radiation"
#> 
#> $TOA_SW_DNI$definition
#> [1] "The total solar irradiance incident (direct plus diffuse) on a horizontal plane where oriented to the sun's position at the top of the atmosphere (extraterrestrial radiation)."
#> 
#> 
#> $ZENITH_LUMINANCE
#> $ZENITH_LUMINANCE$type
#> [1] "RADIATION"
#> 
#> $ZENITH_LUMINANCE$temporal
#> [1] "HOURLY"
#> 
#> $ZENITH_LUMINANCE$source
#> [1] "POWER"
#> 
#> $ZENITH_LUMINANCE$community
#> [1] "AG"
#> 
#> $ZENITH_LUMINANCE$calculated
#> [1] FALSE
#> 
#> $ZENITH_LUMINANCE$inputs
#> NULL
#> 
#> $ZENITH_LUMINANCE$units
#> [1] "cd/m^2"
#> 
#> $ZENITH_LUMINANCE$name
#> [1] "Zenith luminance"
#> 
#> $ZENITH_LUMINANCE$definition
#> [1] "The average amount of luminance at the skys zenith."
```

## A Note on API Throttling

The POWER API endpoints limit queries to prevent overloads due to repetitive and rapid requests.
If you find that the API is throttling your queries, I suggest that you investigate the use of `limit_rate()` from [_ratelimitr_](https://cran.r-project.org/package=ratelimitr) to create self-limiting functions that will respect the rate limits that the API has in place.
It is best to check the [POWER website](https://power.larc.nasa.gov/docs/services/api/#rate-limiting) for the latest rate limits as they differ between temporal APIs and may change over time as the project matures.

## References

<https://power.larc.nasa.gov>

<https://power.larc.nasa.gov/docs/methodology/>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_power.R
\name{get_power}
\alias{get_power}
\title{Get NASA POWER data from the POWER web API}
\usage{
get_power(
  community,
  pars,
  temporal_api = NULL,
  lonlat,
  dates = NULL,
  site_elevation = NULL,
  wind_elevation = NULL,
  wind_surface = NULL,
  temporal_average = NULL
)
}
\arguments{
\item{community}{A character vector providing community name: \dQuote{ag},
\dQuote{re} or \dQuote{sb}.  See argument details for more.}

\item{pars}{A character vector of solar, meteorological or climatology
parameters to download.  When requesting a single point of x, y
coordinates, a maximum of twenty (20) \code{pars} can be specified at one time,
for \dQuote{daily}, \dQuote{monthly} and \dQuote{climatology}
\code{temporal_api}s.  If the \code{temporal_api} is specified as \dQuote{hourly}
only 15 \code{pars} can be specified in a single query.  See \code{temporal_api} for
more.}

\item{temporal_api}{Temporal \acronym{API} end-point for data being queried,
supported values are \dQuote{hourly}, \dQuote{daily}, \dQuote{monthly} or
\dQuote{climatology}.  See argument details for more.}

\item{lonlat}{A numeric vector of geographic coordinates for a cell or region
entered as x, y coordinates.  See argument details for more.}

\item{dates}{A character vector of start and end dates in that order,\cr
\emph{e.g.}, \code{dates = c("1983-01-01", "2017-12-31")}.
Not used when\cr \code{temporal_api} is set to \dQuote{climatology}.
See argument details for more.}

\item{site_elevation}{A user-supplied value for elevation at a single point
in metres.  If provided this will return a corrected atmospheric pressure
value adjusted to the elevation provided.  Only used with \code{lonlat} as a
single point of x, y coordinates, not for use with \dQuote{global} or with
a regional request.}

\item{wind_elevation}{A user-supplied value for elevation at a single point
in metres.  Wind Elevation values in Meters are required to be between 10m
and 300m.  Only used with \code{lonlat} as a single point of x, y coordinates,
not for use with \dQuote{global} or with a regional request.  If this
parameter is provided, the \code{wind-surface} parameter is required with the
request, see
\url{https://power.larc.nasa.gov/docs/methodology/meteorology/wind/}.}

\item{wind_surface}{A user-supplied wind surface for which the corrected
wind-speed is to be supplied.  See \code{wind-surface} section for more detail.}

\item{temporal_average}{Deprecated. This argument has been superseded by
\code{temporal_api} to align with the new \acronym{POWER} \acronym{API}
terminology.}
}
\value{
A data frame as a \code{POWER.Info} class, an extension of the
\link[tibble:tibble]{tibble::tibble}, object of \acronym{POWER} data including location, dates
(not including \dQuote{climatology}) and requested parameters.  A decorative
header of metadata is included in this object.
}
\description{
Get \acronym{POWER} global meteorology and surface solar energy
climatology data and return a tidy data frame \code{\link[tibble:tibble]{tibble::tibble()}}
object.  All options offered by the official \acronym{POWER} \acronym{API}
are supported.  Requests are formed to submit one request per point.  There
is no need to make synchronous requests for multiple parameters for a
single point or regional request.  Requests are limited to 30 unique
requests per 60 seconds.  \CRANpkg{nasapower} attempts to enforce this
client-side.
}
\note{
The associated metadata shown in the decorative header are not saved if
the data are exported to a file format other than a native \R data format,
\emph{e.g.}, .Rdata, .rda or .rds.
}
\section{Argument details for \dQuote{community}}{
 there are three valid
values, one must be supplied. This  will affect the units of the parameter
and the temporal display of time series data.

\describe{
\item{ag}{Provides access to the Agroclimatology Archive, which
contains industry-friendly parameters formatted for input to crop models.}

\item{sb}{Provides access to the Sustainable Buildings Archive, which
contains industry-friendly parameters for the buildings community to include
parameters in multi-year monthly averages.}

\item{re}{Provides access to the Renewable Energy Archive, which contains
parameters specifically tailored to assist in the design of solar and wind
powered renewable energy systems.}
}
}

\section{Argument details for \code{temporal_api}}{
 There are four valid values.
\describe{
\item{hourly}{The hourly average of \code{pars} by hour, day, month and year,
the time zone is UTC.d}
\item{daily}{The daily average of \code{pars} by day, month and year.}
\item{monthly}{The monthly average of \code{pars} by month and year.}
\item{climatology}{Provide parameters as 22-year climatologies (solar)
and 30-year climatologies (meteorology); the period climatology and
monthly average, maximum, and/or minimum values.}
}
}

\section{Argument details for \code{lonlat}}{

\describe{
\item{For a single point}{To get a specific cell, 1/2 x 1/2 degree, supply a
length-two numeric vector giving the decimal degree longitude and latitude
in that order for data to download,\cr
\emph{e.g.}, \code{lonlat = c(-179.5, -89.5)}.}

\item{For regional coverage}{To get a region, supply a length-four numeric
vector as lower left (lon, lat) and upper right (lon, lat) coordinates,
\emph{e.g.}, \code{lonlat = c(xmin, ymin, xmax, ymax)} in that order for a
given region, \emph{e.g.}, a bounding box for the south western corner of
Australia: \code{lonlat = c(112.5, -55.5, 115.5, -50.5)}.  *Maximum area
processed is 4.5 x 4.5 degrees (100 points).}

\item{For global coverage}{To get global coverage for \dQuote{climatology},
supply \dQuote{global} while also specifying \dQuote{climatology} for the
\code{temporal_api}.}
}
}

\section{Argument details for \code{dates}}{
 if one date only is provided, it
will be treated as both the start date and the end date and only a single
day's values will be returned, \emph{e.g.}, \code{dates = "1983-01-01"}.  When
\code{temporal_api} is set to \dQuote{monthly}, use only two year values (YYYY),
\emph{e.g.} \code{dates = c(1983, 2010)}.  This argument should not be used when
\code{temporal_api} is set to \dQuote{climatology} and will be ignored if set.
}

\section{wind-surface}{
 There are 17 surfaces that may be used for corrected
wind-speed values using the following equation:
\deqn{WSC_hgt = WS_10m\times(\frac{hgt}{WS_50m})^\alpha}{WSC_hgt = WS_10m*(hgt/WS_50m)^\alpha}
Valid surface types are described here.

\describe{
\item{vegtype_1}{35-m broadleaf-evergreen trees (70\% coverage)}
\item{vegtype_2}{20-m broadleaf-deciduous trees (75\% coverage)}
\item{vegtype_3}{20-m broadleaf and needleleaf trees (75\% coverage)}
\item{vegtype_4}{17-m needleleaf-evergreen trees (75\% coverage)}
\item{vegtype_5}{14-m needleleaf-deciduous trees (50\% coverage)}
\item{vegtype_6}{Savanna:18-m broadleaf trees (30\%) & groundcover}
\item{vegtype_7}{0.6-m perennial groundcover (100\%)}
\item{vegtype_8}{0.5-m broadleaf shrubs (variable \%) & groundcover}
\item{vegtype_9}{0.5-m broadleaf shrubs (10\%) with bare soil}
\item{vegtype_10}{Tundra: 0.6-m trees/shrubs (variable \%) & groundcover}
\item{vegtype_11}{Rough bare soil}
\item{vegtype_12}{Crop: 20-m broadleaf-deciduous trees (10\%) & wheat}
\item{vegtype_20}{Rough glacial snow/ice}
\item{seaice}{Smooth sea ice}
\item{openwater}{Open water}
\item{airportice}{Airport: flat ice/snow}
\item{airportgrass}{Airport: flat rough grass}
}
}

\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

# Fetch daily "ag" community temperature, relative humidity and precipitation
# for January 1 1985 at Kingsthorpe, Queensland, Australia
ag_d <- get_power(
  community = "ag",
  lonlat = c(151.81, -27.48),
  pars = c("RH2M", "T2M", "PRECTOTCORR"),
  dates = "1985-01-01",
  temporal_api = "daily"
)

ag_d

# Fetch single point climatology for air temperature
ag_c_point <- get_power(
  community = "ag",
  pars = "T2M",
  c(151.81, -27.48),
  temporal_api = "climatology"
)

ag_c_point

# Fetch global ag climatology for air temperature
ag_c_global <- get_power(
  community = "ag",
  pars = "T2M",
  lonlat = "global",
  temporal_api = "climatology"
)

ag_c_global

# Fetch interannual solar cooking parameters for a given region
sse_i <- get_power(
  community = "re",
  lonlat = c(112.5, -55.5, 115.5, -50.5),
  dates = c("1984", "1985"),
  temporal_api = "monthly",
  pars = c("CLRSKY_SFC_SW_DWN", "ALLSKY_SFC_SW_DWN")
)

sse_i
\dontshow{\}) # examplesIf}
}
\references{
\url{https://power.larc.nasa.gov/docs/methodology/}
\url{https://power.larc.nasa.gov}
}
\author{
Adam H. Sparks \email{adamhsparks@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nasapower-package.R
\docType{package}
\name{nasapower-package}
\alias{nasapower}
\alias{nasapower-package}
\title{nasapower: NASA POWER API Client}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Client for 'NASA' 'POWER' global meteorology, surface solar energy and climatology data 'API'. 'POWER' (Prediction Of Worldwide Energy Resource) data are freely available for download with varying spatial resolutions dependent on the original data and with several temporal resolutions depending on the POWER parameter and community. This work is funded through the 'NASA' Earth Science Directorate Applied Science Program. For more on the data themselves, the methodologies used in creating, a web- based data viewer and web access, please see <https://power.larc.nasa.gov/>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/nasapower/}
  \item Report bugs at \url{https://github.com/ropensci/nasapower/issues}
}

}
\author{
\strong{Maintainer}: Adam H. Sparks \email{adamhsparks@gmail.com} (\href{https://orcid.org/0000-0002-0061-8359}{ORCID})

Other contributors:
\itemize{
  \item Scott Chamberlain \email{myrmecocystus@gmail.com} (\href{https://orcid.org/0000-0003-1444-9135}{ORCID}) (Scott Chamberlain reviewed nasapower for rOpenSci, see https://github.com/ropensci/software-review/issues/155) [reviewer]
  \item Hazel Kavili (Hazel Kavili reviewed nasapower for rOpenSci, see https://github.com/ropensci/software-review/issues/155) [reviewer]
  \item Alison Boyer (Alison Boyer reviewed nasapower for rOpenSci, see https://github.com/ropensci/software-review/issues/155) [reviewer]
  \item Fernando Miguez (\href{https://orcid.org/0000-0002-4627-8329}{ORCID}) (Fernando Miguez provided assistance in identifying improper missing value handling in the POWER data, see <https://github.com/femiguez/apsimx/pull/26>) [contributor]
  \item Maëlle Salmon (\href{https://orcid.org/0000-0002-2815-0399}{ORCID}) (Maëlle Salmon contributed a patch to fix issues with using the R package, 'vcr', for testing the 'API' queries, see https://github.com/ropensci/nasapower/pull/64.) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_parameters.R
\name{query_parameters}
\alias{query_parameters}
\title{Query the POWER API for detailed information on parameters}
\usage{
query_parameters(community = NULL, par = NULL, temporal_api = NULL)
}
\arguments{
\item{community}{An optional character vector providing community name:
\dQuote{ag}, \dQuote{sb} or \dQuote{re}.}

\item{par}{An optional character vector of a single solar, meteorological or
climatology parameter to query.  If unsure, omit this argument for for a
full list of all the parameters available for each temporal \acronym{API}
and community.}

\item{temporal_api}{An optional character vector indicating the temporal
\acronym{API} end-point for data being queried, supported values are
\dQuote{hourly}, \dQuote{daily}, \dQuote{monthly} or \dQuote{climatology}.}
}
\value{
A \link{list} object of information for the requested parameter(s) (if
requested), community and temporal \acronym{API}.
}
\description{
Queries the \acronym{POWER} \acronym{API} returning detailed information on
available parameters.
}
\details{
If \code{par} is not provided all possible parameters for the provided
community, \code{community} and temporal \acronym{API}, \code{temporal_api} will be
returned.  If only a single parameter is supplied with no \code{community} or
\code{temporal_api} then the complete attribute information for that parameter
will be returned for all possible communities and temporal \acronym{API}s
combinations.  If all three values are provided, only the information for
that specific combination of parameter, temporal \acronym{API} and community
will be returned.
}
\section{Argument details for \code{temporal_api}}{
 There are four valid values.
\describe{
\item{hourly}{The hourly average of \code{pars} by hour, day, month and year.}
\item{daily}{The daily average of \code{pars} by day, month and year.}
\item{monthly}{The monthly average of \code{pars} by month and year.}
\item{climatology}{Provide parameters as 22-year climatologies (solar)
and 30-year climatologies (meteorology); the period climatology and
monthly average, maximum, and/or minimum values.}
}
}

\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

# fetch the complete set of attribute information for "T2M".
query_parameters(par = "T2M")

# fetch complete temporal and community specific attribute information
# for "T2M" in the "ag" community for the "hourly" temporal API.
query_parameters(par = "T2M",
                 community = "ag",
                 temporal_api = "hourly")

# fetch complete temporal and community specific attribute information
# for all parameters in the "ag" community for the "hourly" temporal API.
query_parameters(community = "ag",
                 temporal_api = "hourly")
\dontshow{\}) # examplesIf}
}
\author{
Adam H. Sparks, \email{adamhsparks@gmail.com}
}
