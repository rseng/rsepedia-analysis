
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

*Wow, no problems at all. :)**Wow, no problems at all. :)*