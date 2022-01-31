# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this 
project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased] yyyy-mm-dd

### Added
  - `nm_to_sm` Convert nautical miles to survey miles. (#99)
  - `tidy_adv` to replace `tidy_fstadv` which will be removed in release 0.2.2 
    (#103)
  - `get_storm_list` returns dataframe of all known cyclones. (#114)

### Changed
  - `gis_download` and `gis_latest` can accept parameters for rgdal::readOGR 
    (#104).
    
  - `ep_prblty_stations` now returns all stations.
  
  - `al_prblty_stations`, `cp_prblty_stations` and `ep_prblty_stations` 
    datasets were modified with additional columns. Changes are documented.

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - `tidy_fcst` now returns all forecast periods; previously only returned hours 12:96 (#107)
  - NHC GIS page went from using POST to GET parameters at some point recently for wind speed probability datasets. Func `gis_wsp` modified accordingly. (#108)

### Security
  - NA

## [0.2.0-6] 2017-07-16

### Added
  - NA

### Changed
  - `get_storms` and `get_storm_data` now use asynchronous http requests to make data collection faster. (#94)
  - Prefaced all built-in data objects with "df"; modified names slightly.

### Removed
  - `load_storm_data` has been removed. Archived data can now be accessed 
    through `rrricanesdata`.

### Deprecated
  - NA

### Fixed
  - NA

### Security
  - NA

## [0.2.0-5.1] 2017-07-10

### Added
  - Instructions for Linux users to install `libgdal1-dev`, `libproj-dev`, and `libxml2-dev`. (#95)

### Added
  - NA

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

### Security
  - NA

## [0.2.0-5] 2017-07-08

### Added
  - NA

### Changed
  - `load_storm_data` now takes `readr::read_csv` parameters.

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

### Security
  - NA

## [0.2.0-4] 2017-06-25

### Added
  - NA

### Changed
  - Added variable `Key` to `discus` dataframes. (#80)
  - Removed variable `Adv` from `posest`. Position estimates do not have advisory numbers. (#81)
  - Fix `scrape_adv_num` to accomodate possible "INTERMEDIATE" text in Public Advisory headers. (#83)
  - Remove variable `Adv` from `update`. Updates do not have advisory numbers. (#84)
  - Added variable `Key` to `get_public` dataframes. (#85)
  - Added variable `Key` to `get_update` dataframes. (#86)
  - Removed non-existent wind radii variables in `get_fstadv`. Hrs 48 and 72 hours only have 34 and 50kt wind fields. Hrs 96 and 120 have none. (#89)

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

### Security
  - NA

## [0.2.0-3] 2017-06-22

### Added
  - Examples for functions `knots_to_mph`, `mb_to_in`, `status_abbr_to_str`, `get_discus`, `get_fstadv`, `tidy_fstadv`, `tidy_wr`, `tidy_fcst` and `tidy_fcst_wr`.

### Changed
  - Added data files to make building vignettes quicker.
  - Added `skip_on_cran` to tests. Additionally, slimmed down some tests. Previous tests exist in branch `tests` and will be redeveloped.
  - Minor documentation updates and corrections.

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

### Security
  - NA

## [0.2.0-2] 2017-06-22

### Added
  - NA

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - Advisories now issued when tropical cyclone development is anticipated, but not yet occurred, and watches and warnings need to be issued. See AL022017, AL032017. 
  - Added additional time zones (HDT, HST, MDT, MST)
  - Appended additional headers (TCE, WTPA, MIATCM, MIAWRKAD1, MIAWRKAP)

### Security
  - NA

## [0.2.0-1] 2017-06-16

### Added
  - GIS functions `gis_advisory`, `gis_breakpoints`, `gis_latest`, `gis_outlook`, `gis_prob_storm_surge`, `gis_windfield` and `gis_wsp` added. These functions return one or more URLs to datasets that can be downloaded with `gis_download`.
  - `shp_to_df` added to convert lines and polygons spatial dataframes to dataframes. Points dataframes can be converted using `tibble::as_dataframe` (target the @data object).

### Changed
  - `load_storm_data` now returns full datasets from the `rrricanesdata` repo including tidied `fstadv` data. See documentation for notes on other products. (#76)

### Removed
  - NA

### Deprecated
  - Not yet deprecated but a warning that `al_prblty_stations`, `cp_prblty_stations` and `ep_prblty_stations` may be removed on a future release. (#46)

### Fixed
  - dplyr 0.6.0 has renamed the .cols parameter of `mutate_at` to .vars. Have modified pkg to accept both dplyr 0.5.0 and >= 0.6.0. This will be removed in future releases. (#74)

### Security
  - NA

## [0.1.3] 2017-06-11

### Added
  - `rrricanes.http_sleep` to control time to sleep between multiple HTTP requests.

### Changed
  - Update documentation for `get_fstadv`, `get_prblty`, `get_wndprb`, `tidy_fstadv`, `tidy_wr`, `tidy_fcst` and `tidy_fcst_wr`.

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - Correct `tidy_fcst` and `tidy_fcst_wr` when all forecast periods do not exist. Previously, was expected that all forecast fields would exist. This may not always be the case. Now works only available forecast periods. (#73)

### Security
  - NA

## [0.1.2] - 2017-06-08

### Added
  - dplyr.progress_bar for all products
  - rrricanes.working_msg option to show current working advisory.
  - `tracking_chart()` for a base world plot. `al_tracking_chart()` for chart centered on Atlantic basin. `ep_tracking_chart()` for chart centered on northeast Pacific.
  - `load_storm_data()` helps get datasets that have already been scraped and processed. Designed to make it more efficient to get data faster.
  - `status_abbr_to_str` converts storm status abbreviations (i.e., TD, TS, HU) to string.
  - `saffir` returns Saffir-Simpson classification of tropical cyclones; abbreviated.
  - `twoal` and `twoep` for Atlantic and east Pacific tropical weather outlooks.
  - Added options `rrricanes.http_timeout` and `rrricanes.http_attempts` to give user more control over failures.

### Changed
  - `get_storm_data` now takes link as first parameter for chaining. Returns a list of dataframes for each product.
  - `tidy_fstadv`, `tidy_wr`, `tidy_fcst` and `tidy_fcst_wr` have been added to replaced now-removed `fstadv_split()`.
  - Completed package top-level documentation.
  

### Removed
  - `fstadv_split`. Dataframe can be split if desired by user. 

### Deprecated
  - NA

### Fixed
  - Fix call to `get_storms` on some Linux distros which generated xpath_element fun error. (#67)
  - Fix call to `get_storm_data`. Issue similar to #67. (#68)
  - Fix call to `gis_wsp`. Call in `rvest::html_nodes` generated "xpath_attrib" error. Add test for `gis_wsp`. (#70)

### Security
  - NA

## [0.1.0] - 2017-05-12

### Added
  - `al_prblty_stations` Get list of locations for wind speed probabilities in Atlantic basin.
  - `cp_prblty_stations` Get list of locations for wind speed probabilities in central Pacific basin.
  - `fstadv_split` split dataframe returned from `fstadv()` to four narrow dataframes.
  - `get_discus` get storm discussions from a storm's archive page.
  - `get_fstadv` get forecast/advisory products from a storm's archive page.
  - `get_nhc_link` returns link to NHC homepage
  - `get_posest` get position estimates from a storm's archive page.
  - `get_prblty` get strike probabilities from a storm's archive page.
  - `get_products` get links to all products for a storm.
  - `get_public` get public advisory statements from a storm's archive page.
  - `get_storms` get a list for storms from a year's archive page.
  - `get_storm_data` get one or multiple products for a storm
  - `get_update` get updates from a storm's archive page.
  - `get_wndprb` get wind speed probabilities from a storm's archive page.
  - `knots_to_mph` Convert values from knots to mph (for wind and gust values).
  - `mb_to_in` convert barometric pressure from millibars to inches.
  - `wndprb` Access a specific wind speed probability for a storm.

### Changed
  - Correct version, CHANGELOG and NEWS from previous "release".

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

### Security
  - NA

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![GitHub
(pre-)release](https://img.shields.io/github/release/ropensci/rrricanes/all.svg)](https://github.com/ropensci/rrricanes/tags)
[![](https://badges.ropensci.org/118_status.svg)](https://github.com/ropensci/onboarding/issues/118)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rrricanes)](https://cran.r-project.org/package=rrricanes)
[![Build
Status](https://img.shields.io/travis/ropensci/rrricanes/master.svg)](https://travis-ci.org/ropensci/rrricanes)
[![AppVeyor Build
Status](https://img.shields.io/appveyor/ci/timtrice/rrricanes-g4dos/master.svg)](https://ci.appveyor.com/project/timtrice/rrricanes-g4dos)
[![codecov](https://codecov.io/gh/ropensci/rrricanes/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rrricanes)

# rrricanes <img src='man/figures/logo.png' align="right" height="138" />

`rrricanes` is a R library that extracts information from [available
archives](http://www.nhc.noaa.gov/archive/1998/1998archive.shtml) on
past and current tropical cyclones. Currently, archives date back to
1998.

Data can be obtained for cyclones in the north Atlantic (considered the
Atlantic Basin) and north-eastern Pacific (the East Pacific Basin from
140°W and eastward.

Central Pacific data (140°W to 180°W) is included if issued by the
[National Hurricane Center](http://www.nhc.noaa.gov/) (generally they’re
issued by the [Central Pacific Hurricane
Center](http://www.prh.noaa.gov/cphc/)).

This library parses the text advisories of all tropical cyclones since
1998. Over the years the formats of the text products have changed and
many are loosely formatted.

I wrote this package with the goal of consolidating messy text data into
well-organized formats that can easily be saved to CSV, SQL and other
data formats.

You may explore some features of the package through the
[shinycanes](https://timtrice.shinyapps.io/shinycanes/) beta web
application (built with R Shiny).

## Advisory Products

Generally speaking, there are five products available for tropical
cyclones issued at 03:00, 09:00, 15:00 and 21:00 UTC;

1.  Storm Discussion - These are technical discussions centered on the
    current structure of the cyclone, satellite presentation, computer
    forecast model tendencies and more.

2.  Forecast/Adivsory - This data-rich product lists the current
    location of the cyclone, its wind structure, forecast and forecast
    wind structure.

3.  Public Advisory - These are general text statements issued for the
    public-at-large. Information in these products is a summary of the
    Forecast/Advisory product along with any watches and warnings
    issued, changed, or cancelled. Public Advisory products are the only
    regularly-scheduled product that may be issued intermittently (every
    three hours and, occasionally, every two hours) when watches and
    warnings are in effect.

4.  Wind Speed Probabilities - These products list the probability of a
    minimum sustained wind speed expected in a given forecast window.
    This product replaces the Strike Probabilities product beginning in
    2006 (see below).

5.  Updates - Tropical Cyclone Updates may be issued at any time if a
    storm is an immediate threat to land or if the cyclone undergoes a
    significant change of strength or structure. The information in this
    product is general.

**Discontinued Products**

These products are included in the package though they have been
discontinued at some point:

1.  Strike Probabilities - List the probability of a tropical cyclone
    passing within 65 nautical miles of a location within a forecast
    window. Replaced in 2006 by the Wind Speed Probabilities product.

2.  Position Estimates - Typically issued as a storm is threatening land
    but generally rare (see Hurricane Ike 2008, Key AL092008). It is
    generally just an update of the current location of the cyclone.
    After the 2011 hurricane season, this product was discontinued;
    Updates are now issued in their place.

## Getting Started

Please view the vignette ‘Getting Started’:

``` r
vignette("getting_started", package = "rrricanes")
```

[Online documentation](https://timtrice.github.io/rrricanes/) is also
available.

### Prerequisites

`rrricanes` requires an active internet connection as data is extracted
from online sources.

Linux users must also have the `libgdal-dev`, `libproj-dev` and
`libxml2-dev` packages installed.

To add `rrricanesdata`, a [package of post-scraped
datasets](https://github.com/ropensci/rrricanesdata),

``` r
install.packages("rrricanesdata", 
                 repos = "https://timtrice.github.io/drat/", 
                 type = "source")
```

To use high resolution tracking maps you will need to install the
`rnaturalearthhires` package.

``` r
install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")
```

### Installing

`rrricanes` is currently only available in GitHub. It can be installed
using the `devtools` package:

``` r
devtools::install_github("ropensci/rrricanes", build_vignettes = TRUE)
```

## Built With

  - [R 3.3.3](https://www.r-project.org/) - The R Project for
    Statistical Computing

## Contributing

Please read
[CONTRIBUTING.md](https://github.com/ropensci/rrricanes/blob/master/.github/CONTRIBUTING.md)
for details on our code of conduct, and the process for submitting pull
requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions
available, see the [tags on this
repository](https://github.com/ropensci/rrricanes/tags).

## Authors

  - **Tim Trice** - *Initial work* -
    [timtrice](https://github.com/timtrice)

See also the list of
[contributors](https://github.com/ropensci/rrricanes/contributors) who
participated in this project.

## License

This project is licensed under the MIT License - see the
[LICENSE.md](LICENSE.md) file for details

## Acknowledgments

  - [Molyneux, James](https://github.com/jimmylovestea)
  - [Padgham, Mark](https://github.com/mpadge)
  - [Robinson, Emily](https://github.com/robinsones)
  - [Rudis, Bob](https://github.com/hrbrmstr)
  - [Salmon, Maëlle](https://github.com/maelle)
  - [Stachelek, Joseph](https://github.com/jsta)

## Known Data Quality Issues

1.  Hurricane Juan (AL152003), Adv 15; no status leads to improper
    `Status` and `Name` values in some datasets.
    ([\#82](https://github.com/ropensci/rrricanes/issues/82))
rrricanes Unreleased (yyyy-mm-dd)
==================================

### NEW FEATURES

* Add `get_storm_list` that returns a list of, it would appear, all known storms for many basins (beyond AL, EP)
  + Currently only returns 2582 cyclones despite more storms present in the file. This is due to a currenty unknown issue in the raw CSV file.

* `get_ftp_storm_data`, given a `stormid` (key) and a product, will return a dataframe of the requested data if it exists for the current storm. This opens up cyclones earlier than 1998 **however** is not yet tested.

### MINOR IMPROVEMENTS

* Add `nm_to_sm` to convert nautical miles to survey miles. (#99)

* Allow parameters for rgdal::readOGR (#104)

* Add `tidy_adv` to replace `tidy_fstadv`; `tidy_fstadv` will be removed in 
  release 0.2.2 (#103)

* `al_prblty_stations` and `cp_prblty_stations` updated 
  to accomodate new dataset. (#46)
  
* `ep_prblty_stations` now returns dataset of east Pacific stations. (#46)

### BUG FIXES

* `tidy_fcst` only returned forecasts for hours 12:96. Now returns 120. (#107)
* `gis_wsp` now uses GET instead of POST to retrieve wind speed probability GIS datasets. (#108)

### DEPRECATED AND DEFUNCT

* NA

rrricanes 0.2.0-6 (2017-07-16)
==================================

### NEW FEATURES

* NA

### MINOR IMPROVEMENTS

* `get_storms` and `get_storm_data` have been rewritten to utilize pkg `crul`'s asynchronous features. This will not make much of a difference in `get_storms` (and may actually be slightly slower; to be explained). But the difference with `get_storm_data` should be very noticeable. There is a limit to hitting the NHC archives; 80 requests per 10 seconds. Both functions send 4 links through every 0.5 seconds to avoid this limit. Timeout issues should no longer occur so options rrricanes.http_attempts and rrricanes.http_timeout have been removed. The primary cause of long processing now is due to scraping, particularly with the `fstadv` products; the amount of data in these products and the unstructured nature of the products require a number of rules. This can probably be simplified in future releases. (#94)

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* `load_storm_data` has been removed in favor of loading `rrricanesdata`. See 
  vignette "installing_rrricanesdata" for more information. 

rrricanes 0.2.0-5 (2017-07-08)
==================================

### NEW FEATURES

* NA

### MINOR IMPROVEMENTS

* `load_storm_data` now takes `readr::read_csv` parameters.

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanes 0.2.0-4 (2017-06-25)
==================================

### NEW FEATURES

* NA

### MINOR IMPROVEMENTS

* `Key` variable added to `discus` dataframes. `Key` will be NA for all cyclones >= 2005. Should not be <= 2006. (#80)
* Removed `Adv` variable from `posest` dataframes. Position estimates do not have advisory numbers. (#81)
* Removed `Adv` variable from `update`. Updates do not have advisory numbers. (#84)
* Added variable `Key` to `get_public` dataframes. (#85)
* Added variable `Key` to `get_update` dataframes. (#86)
* Removed non-existent wind radii variables in `get_fstadv`. Hrs 48 and 72 hours only have 34 and 50kt wind fields. Hrs 96 and 120 have none. (#89)

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanes 0.2.0-3 (2017-06-22)
==================================

### NEW FEATURES

* NA

### MINOR IMPROVEMENTS

* Examples for functions `knots_to_mph`, `mb_to_in`, `status_abbr_to_str`, `get_discus`, `get_fstadv`, `tidy_fstadv`, `tidy_wr`, `tidy_fcst` and `tidy_fcst_wr`.
* Minor documentation updates
* Small enhancements to tests, vignettes

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanes 0.2.0-2 (2017-06-22)
==================================

### NEW FEATURES

* NA

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* Advisories now issued when tropical cyclone development is anticipated, but not yet occurred, and watches and warnings need to be issued. See AL022017, AL032017. 
* Added additional time zones (HDT, HST, MDT, MST)
* Appended additional headers (TCE, WTPA, MIATCM, MIAWRKAD1, MIAWRKAP)

### DEPRECATED AND DEFUNCT

* NA

rrricanes 0.2.0-1 (2017-06-16)
==================================

### NEW FEATURES

* GIS functions now available. Please note some of these products may not exist for every available cyclone/advisory.
 + `gis_advisory` Typically will include current and past track data, forecast track data, forecast cone (margin of error) and wind radius data. 
 + `gis_breakpoints` List of breakpoints typically used for watch/warning areas but is not a requirement.
 + `gis_latest` Retrieves the latest GIS products for all active storms. 
 + `gis_outlook` Retrives the latest tropical weather outlook in shapefile format.
 + `gis_prob_storm_surge` Probabilistic storm surge; a polygon dataset for psurge and esurge products with various criteria.
 + `gis_windfield` Wind radius datasets.
 + `gis_wsp` Wind speed probabilities
 + `gis_download` Use this function to download the URLs returned from the above functions.
 + `shp_to_df` added to convert lines and polygons spatial dataframes to dataframes. Points dataframes can be converted using `tibble::as_dataframe` (target the @data object).

### MINOR IMPROVEMENTS

* [Enhanced documentation](https://ropensci.github.io/rrricanes/) added online using `pkgdown`. 
* `load_storm_data` directly returns dataframes. Additionally, retrieval by basin and years removed in favor of importing complete product datasets. Additionally, documentation has been added to the website on [using data.world](https://ropensci.github.io/rrricanes/articles/articles/data_world.html) as a third option. The difference between these two options is `load_storm_data` will return complete datasets. Using data.world will allow users to write custom queries to retrieve data.  (#76)

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* Not yet deprecated but a warning that `al_prblty_stations`, `cp_prblty_stations` and `ep_prblty_stations` may be removed on a future release. (#46)
* Support for dplyr 0.5.0 will be removed in future releases in favor of dplyr 0.7.0.

rrricanes 0.1.3 (2017-06-11)
============================

### NEW FEATURES

* NA

### MINOR IMPROVEMENTS

* `rrricanes.http_sleep` to control time to sleep between multiple HTTP requests.
* Clarified documentation for `get_fstadv`, `get_prblty`, `get_wndprb`, `tidy_fstadv`, `tidy_wr`, `tidy_fcst` and `tidy_fcst_wr`.

### BUG FIXES

* `tidy_fcst` and `tidy_fcst_wr` would err if all forecast periods were not available for a cyclone. Functions now analyze dataframe to determine what forecast fields exist, then tidies based on the result. (#73)

### DEPRECATED AND DEFUNCT

* NA

rrricanes 0.1.2 (2017-06-08)
============================

### NEW FEATURES

* Changed name from `Hurricanes` to `rrricanes`.

* `get_storm_data` can now be chained to other commands and returns a list of dataframes.

* `load_storm_data` accesses pre-scraped datasets and returns requested products through the github repo `rrricanesdata`. This was done to make it quicker to get data. It should not be relied on to get the most immediate data for current storms. However, it should be fairly up-to-date. Original functions can be used if for some reason immediate data access is needed.

* `saffir` returns Saffir-Simpson classification of tropical cyclones; abbreviated.

* `status_abbr_to_str` converts storm status abbreviations (i.e., TD, TS, HU) to string.

* `twoal` and `twoep` parse tropical weather outlook XML files. Gives current status, if any, of areas of interest in either basin.

### MINOR IMPROVEMENTS

* Modified numerous regex patterns to ensure data quality.
* `tidy_fstadv`, `tidy_wr`, `tidy_fcst` and `tidy_fcst_wr` have been added to replaced now-removed `fstadv_split()`.
* Added loop to make multiple attempts at extracting contents from NHC archives. Options `rrricanes.http_timeout` and `rrricanes.http_attempts` added to give user more control over this. Default is 3 attempts with no more than 5 permitted.

### BUG FIXES

* Too many to recall. Apologies. 
* Call to `get_storms` on some linux distros generated xpath_element error. Corrected. (#67)
* Modify call to `get_storm_data`. Replaced css parameter in `rvest::html_nodes` calls with xpath parameter. Some products (notably, `get_prblty`) do not have a "pre" tag but are text documents (not HTML). Modified `scrape_contents` to return full contents if "pre" tag doesn't exist. Tested `get_discus` and `get_public`; no errors generated. (#68)

### DEPRECATED AND DEFUNCT

* `fstadv_split`. See MINOR IMPROVEMENTS for alternatives.

Hurricanes 0.1.0 (2017-05-12)
================

## Major new features

Retrieve all storm's for a given year (>=1998) and access data from a given storm's history. Can access "current" storm position, structure details, forecast, and discussions.

This release should be considered beta. While I've made every effort to ensure quality there may be an issue here or there. I will work on developing QA/QC scripts as time permits.

Please send any issues or questions to: https://github.com/timtrice/Hurricanes/issues.

### Getting Annual Storm Data

Use `get_storms` to access storm's for a given year.

### Getting Storm Data

Use `get_storm_data` to access one or multiple products for a specific storm.

#### Storm Discussions (discus)

Not parsed but contains technical information on the cyclone, development tendencies and forecast model tendencies.

#### Forecast/Advisory (fstadv)

Contains the meat of data. Current storm information, forecast information, wind and sea data. Can use `fstadv_split()` to break the wide dataframe to multiple, relational dataframes.

#### Position Estimate (posest)

Contains current position estimate for a given storm. Usually issued during threats to land. Not issued for all storms. Not parsed. 

#### Strike Probabilities (prblty)

Strike probabilities for given locations prior to 2006 (See Wind Speed Probabilities for >= 2006). 

#### Public Advisory (public)

General information on a storm. Not parsed.

#### Updates (update)

Quick information given when a storm is threatening or undergoes a significant change. Not issued for all storms.

#### Wind Probabilities (wndprb)

Replaced Strike Probabilities after the 2005 hurricane season. Lists the chances of a location seeing 34kt, 50kt and 64kt winds within a given time frame.
### Error Message

### Reproducible Example

### Traceback

### Session Info# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change.

Please note we have a code of conduct, please follow it in all your interactions with the project.

## Pull Request Process

1. Ensure any install or build dependencies are removed before the end of the layer when doing a
   build.
2. Update the README.md with details of changes to the interface, this includes new environment
   variables, exposed ports, useful file locations and container parameters.
3. Increase the version numbers in any examples files and the README.md to the new version that this
   Pull Request would represent. The versioning scheme we use is [SemVer](http://semver.org/).
4. You may merge the Pull Request in once you have the sign-off of two other developers, or if you
   do not have permission to do that, you may request the second reviewer to merge it for you.

## Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at tim.trice@gmail.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/### Overview

What problem or feature request does this pull request resolve?

### Issues Addressed

 * #

### Changes Addressed/Proposed

 * 

### Reproducible Examples

Demonstrate how the new code performs.

### Notes

Optional. Ancillary topics, caveats, alternative strategies that didn't work out, anything else.

### Tests

Include test cases, and expected output
## Run

```
docker run \
  -dti \
  -e DISABLE_AUTH=true \
  -p 8787:8787 \
  --name rrricanes \
  -v /home/timtrice/Projects/ropensci/rrricanes:/home/rstudio/rrricanes \
  timtrice/rrricanes:release
```

## Shell

```
docker exec -ti rrricanes /bin/bash
```
---
output: github_document
---

[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/) 
[![GitHub (pre-)release](https://img.shields.io/github/release/ropensci/rrricanes/all.svg)](https://github.com/ropensci/rrricanes/tags)
[![](https://badges.ropensci.org/118_status.svg)](https://github.com/ropensci/onboarding/issues/118)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rrricanes)](https://cran.r-project.org/package=rrricanes)
[![Build Status](https://img.shields.io/travis/ropensci/rrricanes/master.svg)](https://travis-ci.org/ropensci/rrricanes)
[![AppVeyor Build Status](https://img.shields.io/appveyor/ci/timtrice/rrricanes-g4dos/master.svg)](https://ci.appveyor.com/project/timtrice/rrricanes-g4dos)
[![codecov](https://codecov.io/gh/ropensci/rrricanes/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rrricanes) 

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# rrricanes <img src='man/figures/logo.png' align="right" height="138" />

`rrricanes` is a R library that extracts information from [available archives](http://www.nhc.noaa.gov/archive/1998/1998archive.shtml) on past and current tropical cyclones. Currently, archives date back to 1998. 

Data can be obtained for cyclones in the north Atlantic (considered the Atlantic Basin) and north-eastern Pacific (the East Pacific Basin from 140&deg;W and eastward. 

Central Pacific data (140&deg;W to 180&deg;W) is included if issued by the [National Hurricane Center](http://www.nhc.noaa.gov/) (generally they're issued by the [Central Pacific Hurricane Center](http://www.prh.noaa.gov/cphc/)).

This library parses the text advisories of all tropical cyclones since 1998. Over the years the formats of the text products have changed and many are loosely formatted. 

I wrote this package with the goal of consolidating messy text data into well-organized formats that can easily be saved to CSV, SQL and other data formats. 

You may explore some features of the package through the [shinycanes](https://timtrice.shinyapps.io/shinycanes/) beta web application (built with R Shiny).

## Advisory Products

Generally speaking, there are five products available for tropical cyclones issued at 03:00, 09:00, 15:00 and 21:00 UTC;  

1. Storm Discussion - These are technical discussions centered on the current structure of the cyclone, satellite presentation, computer forecast model tendencies and more. 

2. Forecast/Adivsory - This data-rich product lists the current location of the cyclone, its wind structure, forecast and forecast wind structure.

3. Public Advisory - These are general text statements issued for the public-at-large. Information in these products is a summary of the Forecast/Advisory product along with any watches and warnings issued, changed, or cancelled. Public Advisory products are the only regularly-scheduled product that may be issued intermittently (every three hours and, occasionally, every two hours) when watches and warnings are in effect.

4. Wind Speed Probabilities - These products list the probability of a minimum sustained wind speed expected in a given forecast window. This product replaces the Strike Probabilities product beginning in 2006 (see below).

5. Updates - Tropical Cyclone Updates may be issued at any time if a storm is an immediate threat to land or if the cyclone undergoes a significant change of strength or structure. The information in this product is general.

__Discontinued Products__

These products are included in the package though they have been discontinued at some point:

1. Strike Probabilities - List the probability of a tropical cyclone passing within 65 nautical miles of a location within a forecast window. Replaced in 2006 by the Wind Speed Probabilities product.

2. Position Estimates - Typically issued as a storm is threatening land but generally rare (see Hurricane Ike 2008, Key AL092008). It is generally just an update of the current location of the cyclone. After the 2011 hurricane season, this product was discontinued; Updates are now issued in their place.

## Getting Started

Please view the vignette 'Getting Started':

```r
vignette("getting_started", package = "rrricanes")
```

[Online documentation](https://timtrice.github.io/rrricanes/) is also available.

### Prerequisites

`rrricanes` requires an active internet connection as data is extracted from online sources.

Linux users must also have the `libgdal-dev`, `libproj-dev` and `libxml2-dev` packages installed.

To add `rrricanesdata`, a [package of post-scraped datasets](https://github.com/ropensci/rrricanesdata), 

```r
install.packages("rrricanesdata", 
                 repos = "https://timtrice.github.io/drat/", 
                 type = "source")
```

To use high resolution tracking maps you will need to install the `rnaturalearthhires` package. 

```r
install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")
```

### Installing

`rrricanes` is currently only available in GitHub. It can be installed using the `devtools` package:

```r
devtools::install_github("ropensci/rrricanes", build_vignettes = TRUE)
```

## Built With

* [R 3.3.3](https://www.r-project.org/) - The R Project for Statistical Computing

## Contributing

Please read [CONTRIBUTING.md](https://github.com/ropensci/rrricanes/blob/master/.github/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/ropensci/rrricanes/tags). 

## Authors

* **Tim Trice** - *Initial work* - [timtrice](https://github.com/timtrice)

See also the list of [contributors](https://github.com/ropensci/rrricanes/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [Molyneux, James](https://github.com/jimmylovestea)
* [Padgham, Mark](https://github.com/mpadge)
* [Robinson, Emily](https://github.com/robinsones)
* [Rudis, Bob](https://github.com/hrbrmstr)
* [Salmon, Maëlle](https://github.com/maelle)
* [Stachelek, Joseph](https://github.com/jsta)

## Known Data Quality Issues

1. Hurricane Juan (AL152003), Adv 15; no status leads to improper `Status` and `Name` values in some datasets. ([#82](https://github.com/ropensci/rrricanes/issues/82))
---
title: "Forecast/Advisory GIS"
author: "Tim Trice"
date: "June 11, 2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
    echo = TRUE, 
    fig.width = 7, 
    fig.asp = 1, 
    fig.align = "center"
)
```

```{r, message = FALSE}
library(dplyr)
library(ggplot2)
library(rrricanes)
library(rrricanesdata)
library(sp)
```

```{r}
key <- "AL092008"
adv <- 42
```

```{r}
fstadv <- fstadv %>% filter(Key == key, Adv <= adv)
```

### GIS Advisory Forecast Track, Cone of Uncertainty, and Watches/Warnings

```{r}
gis_adv <- gis_advisory(key = key, advisory = adv) %>% gis_download()
```

Get bounding box of the forecast polygon.

```{r}
bbox <- bbox(gis_adv$al092008.042_5day_pgn)
```

Generate a base plot of the Atlantic ocean.

```{r}
(bp <- al_tracking_chart(color = "black", fill = "white", size = 0.1, res = 50))
```

I like to add a little cushion for the map inset and forecast cone data. 

```{r}
lat_min <- bbox[2,1] - 5
lat_max <- bbox[2,2] + 5
lon_min <- bbox[1,1] - 10
lon_max <- bbox[1,2] + 10
```

Build a thin tracking map for the inset.

```{r}
bp_inset <- ggplotGrob(bp +
                           geom_rect(mapping = aes(xmin = lon_min, xmax = lon_max,
                                                   ymin = lat_min, ymax = lat_max),
                                     color = "red", alpha = 0) +
                           theme_bw() +
                           theme(axis.title = element_blank(),
                                 axis.ticks = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_blank(),
                                 plot.margin = margin(0, 0, 0, 0, "pt")))
```

Modify original `bp` zoomed in on our area of interest.

```{r}
(bp <- bp +
     coord_equal(xlim = c(lon_min, lon_max),
                 ylim = c(lat_min, lat_max)) +
     scale_x_continuous(expand = c(0, 0)) +
     scale_y_continuous(expand = c(0, 0)) +
     labs(x = "Lon",
          y = "Lat",
          caption = sprintf("rrricanes %s", packageVersion("rrricanes"))))
```

Combine `bp` and `bp_inset` to finalize initial base plot. `bp` will be a base plot without the inset. `bpi` will have the inset.

```{r}
(bpi <- bp + annotation_custom(grob = bp_inset, xmin = lon_max - 5,
                               xmax = lon_max - 1, ymin = -Inf,
                               ymax = lat_min + 5))
```

## Current Advisory Details

Lines and Polygons spatial dataframes can be helpfully converted using `shp_to_df`. The original spatial dataframes can be plotted directly in `ggplot2` but, to my understanding, access to the other variables are not available.

```{r}
# Convert object SpatialLinesDataframe to dataframe
shp_storm_lin <- shp_to_df(gis_adv$al092008.042_5day_lin)
shp_storm_ww <- shp_to_df(gis_adv$al092008.042_ww_wwlin)
# Convert object SpatialPolygonsDataframe to dataframe
shp_storm_pgn <- shp_to_df(gis_adv$al092008.042_5day_pgn)
```

Points dataframes can just be converted with `tibble::as_data_frame`. 

```{r}
# Convert object SpatialPointsDataframe to dataframe
shp_storm_pts <- as_data_frame(gis_adv$al092008.042_5day_pts)
```

Modify `shp_storm_pts$DVLBL` with full strings and ordered factor.

```{r}
shp_storm_pts$DVLBL <- factor(shp_storm_pts$DVLBL, 
                              levels = c("D", "S", "H"), 
                              labels = c("Tropical Depression", 
                                         "Tropical Storm", 
                                         "Hurricane"))
```

Same with `shp_storm_pts$TCWW`:

```{r}
shp_storm_ww$TCWW <- factor(shp_storm_ww$TCWW, 
                            levels = c("TWA", "TWR", "HWA", "HWR"), 
                            labels = c("Tropical Storm Watch", 
                                       "Tropical Storm Warning", 
                                       "Hurricane Watch", 
                                       "Hurricane Warning"))
```

```{r}
bpi + geom_polygon(data = shp_storm_pgn, 
                   aes(x = long, y = lat, group = group),
                   alpha = 0.15, fill = "orange") + 
    geom_path(data = shp_storm_lin, aes(x = long, y = lat, group = group)) + 
    geom_point(data = shp_storm_pts, aes(x = LON, y = LAT, fill = DVLBL,
                                         shape = DVLBL, size = MAXWIND)) + 
    geom_path(data = shp_storm_ww, aes(x = long, y = lat, color = TCWW, 
                                       group = group), size = 1) + 
    scale_shape_manual(values = c(21, 21, 21, 21)) + 
    guides(shape = guide_legend(override.aes = list(size = 3)), 
           size = guide_legend(nrow = 1)) + 
    theme(legend.position = "bottom", 
          legend.box = "vertical")
```

Very often, areas that are under a hurricane watch may also be under a tropical storm warning. The chart above does not show the hurricane watch area.

---
title: "Wind and Pressure"
author: "Tim Trice"
date: "June 14, 2017"
output: html_document
references:
- id: sshws
  title: The Saffir-Simpson Hurricane Wind Scale
  author:
  - family: Schott
    given: Timothy
  - family: Landsea
    given: Chris
  - family: Hafele
    given: Gene
  - family: Lorens
    given: Jeffrey
  - family: Taylor
    given: Arthur
  - family: Thurm
    given: Harvey
  - family: Ward
    given: Bill
  - family: Willis
    given: Mark
  - family: Zaleski
    given: Walt
  URL: 'http://www.nhc.noaa.gov/pdf/sshws.pdf'
  issued:
    year: 2012
    month: 2
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      fig.align = "center")
```

```{r, message = FALSE}
library(dplyr)
library(gganimate)
library(ggplot2)
library(rrricanes)
library(rrricanesdata)
library(tidyr)
```

```{r}
key <- "AL092008"
adv <- 42
```

### Forecast/Advisory Data for Hurricane Ike Adv #42

```{r}
fstadv <- fstadv %>% filter(Key == key, Adv <= adv)
```

## Wind Profile

```{r}
# Plot wind values
fstadv %>% ggplot(aes(x = Date, y = Wind)) + 
  geom_line() + 
  geom_point(aes(color = Status), size = 3) + 
  scale_y_continuous(name = "Wind (kts)") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.box = "vertical") +
  labs(title = "Wind Profile", 
       caption = sprintf("rrricanes %s", packageVersion("rrricanes")))
```

## Pressure Profile

```{r}
# Plot pressure values
fstadv %>% ggplot(aes(x = Date, y = Pressure)) + 
  geom_line() + 
  geom_point(aes(color = Status), size = 3) + 
  scale_y_continuous(name = "Pressure (mb)") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.box = "vertical") +
  labs(title = "Pressure Profile", 
       caption = sprintf("rrricanes %s", packageVersion("rrricanes")))
```

## Wind/Pressure Relational Change

```{r}
fstadv %>%
  mutate(WindDist = (Wind - min(Wind))/(max(Wind) - min(Wind)),
         PressDist = (Pressure - max(Pressure))/(max(Pressure) - min(Pressure))) %>%
  gather(Var, Val, WindDist, PressDist) %>% 
  ggplot(aes(x = Date, y = Val, group = Var, color = Var)) + 
  geom_line(size = 1) + 
  scale_color_discrete(labels = c("Pressure Change", "Wind Change")) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank()) + 
  labs(title = "Wind/Pressure Relational Change", 
       subtitle = "",
       caption = sprintf("rrricanes %s", packageVersion("rrricanes")), 
       y = "")
```

## Wind Radius

```{r}
wr_animate <- 
  fstadv %>% 
  tidy_wr() %>% 
  gather(Quadrant, Radius, NE:NW) %>% 
  ggplot(
    aes(
      x = Quadrant, 
      y = Radius, 
      fill = factor(WindField), 
      frame = Adv
    )
  ) + 
  geom_bar(stat = "identity", position = "identity", width = 1) + 
  guides(fill = guide_legend(title = "Wind Field")) +
  coord_polar() + 
  theme_minimal() + 
  theme(legend.position = "bottom") + 
  labs(
    title = "Wind Radius for Advisory {frame}", 
    subtitle = "Minimum sustained one-minute wind speed in knots",
    caption = sprintf("rrricanes %s", packageVersion("rrricanes")), 
    y = "Radius (nm)"
  ) + 
  transition_time(Date) + 
  ease_aes('linear')

animate(wr_animate, renderer = magick_renderer())
```

One of the interesting things to note about the image above; Hurricane Ike was known for it's very large wind field (relatively speaking) which generated a larger and wider storm surge than normal for it's classification. You can see this very well defined structure expansion between advisories 38 and 42. 

This, in part, led to the modification of the Saffir Simpson Hurricane Scale. 

> the very large Hurricane Ike (with hurricane force winds extending as much as 125 mi from the center) in 2008 made landfall in Texas as a Category 2 hurricane and had peak storm surge values of about 20 ft. In contrast, tiny Hurricane Charley (with hurricane force winds extending at most 25 mi from the center) struck Florida in 2004 as a Category 4 hurricane and produced a peak storm surge of only about 7 ft. These storm surge values were substantially outside of the ranges suggested in the original scale. Thus to help reduce public confusion about the impacts associated with the various hurricane categories as well as to provide a more scientifically defensible scale, the storm surge ranges, flooding impact and central pressure statements are being removed from the scale and only peak winds are employed in this revised version – the Saffir-Simpson Hurricane Wind Scale. [@sshws]

## References
---
title: "GIS data"
author: "Tim Trice"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{GIS data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, message = FALSE}
library(dplyr)
library(ggplot2)
library(rrricanes)
library(sf)
```

Most storms will contain a variation of GIS datasets that can be plotted with `ggplot2`. The helper functions for this have the prefix 'gis'.

**All products are experimental and there maybe fluctuations particularly in current datasets.**

In general, datasets are available for storms dated back to 1998. However, products such as Wind Speed Probabilities only go back to 1999. 

Some datasets require the use of the storm key and an optional advisory. Other products require a datetime value and cannot be isolated by storm key or advisory. The datetime values are not based on the issue time of the advisory, but rather three hours prior. For example, if you are seeking a dataset where the forecast/advisory was issued at 9:00AM UTC, you will want the dataset for 6:00AM UTC. This will be explained a little further below.

## Build a Tracking Chart

There are three functions available to help you plot GIS data; `tracking_chart`, `al_tracking_chart` and `ep_tracking_chart`. `al_tracking_chart` and `ep_tracking_chart` are just helpers centered on the Atlantic and northeast Pacific ocean, respectively. 

```{r}
args(rrricanes::tracking_chart)
```

The `countries` and `states` parameters are TRUE by default. This means a basic call to `tracking_chart` will return a map with country and state borders. The `res` parameter is resolution; one of 110, 50 or 10 nautical miles. Resolutions 110nm and 50nm can be used immediately. To use lower resolution you must install the `rnaturalearthdatahires` package from ropensci:

```{r, eval = FALSE}
install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")
```

`tracking_chart` will print a basic tracking chart (a map of the planet).

```{r, fig.width = 7, fig.align = "center"}
tracking_chart()
```

You can pass typical `aes` parameters to refine the color and fill of the plot; remember the tracking chart is a ggplot object.

```{r, fig.width = 7, fig.align = "center"}
tracking_chart(color = "black", fill = "white", size = 0.1)
```

You may choose to only show coastline data instead. In this case, just set the countries parameter to FALSE.

```{r, fig.width = 7, fig.align = "center"}
tracking_chart(countries = FALSE, res = 50, color = "black", fill = "white", 
               size = 0.1)
```

For the purposes of this vignette we'll focus on Atlantic storms.

```{r, fig.width = 7, fig.align = "center"}
(p <- al_tracking_chart(color = "black", fill = "white", size = 0.1, res = 50))
```

The `res` parameter defines the resolution of the chart presented. Options are in 110nm, 50nm and 10nm. The lower the resolution the longer the chart takes to be built. 

States cannot be drawn on resolution greater than 50nm. 

## GIS Datasets

There are several datasets that are published for active cyclones. The following functions are designed to return the URL to those datasets:

* `gis_advisory`
* `gis_prob_storm_surge`
* `gis_windfield`
* `gis_latest`

## Advisory Package

```{r}
gis_advisory(key = "AL182012", advisory = "18")
```

```{r, eval = FALSE}
df.gis_adv <- gis_advisory(key = "AL182012", advisory = "18") %>% gis_download()
```

```{r}
names(df.gis_adv)
```

For this particular storm and advisory, included are the base line, point and polygon datasets along with a dataset for watches and warnings. The objects returned are spatial dataframes contained within the list of dataframes, `df.gis_adv`.

### Line track

```{r}
str(df.gis_adv$al182012.018_5day_lin)
```

As we're dealing with a SpatialLinesDataFrame which needs to be modified, you can use the helper function `shp_to_df` for plotting.

```{r}
fcst_line <- shp_to_df(df.gis_adv$al182012.018_5day_lin)
```

```{r, fig.width = 7, fig.align = "center"}
p + geom_path(data = fcst_line, aes(long, lat, group = FCSTPRD))
```

There are two groups of data in the set: one for 72-hour forecast period and one for 120-hour. Grouping by `FCSTPRD` will show the forecast track correctly.

There is a pretty wide field of view on the map above. You can use `sp::bbox` to "zoom in" on the map.

```{r}
(bb <- sp::bbox(df.gis_adv$al182012.018_5day_lin))
```

```{r, fig.width = 7, fig.align = "center"}
(p2 <- p + geom_path(data = fcst_line, aes(long, lat, group = FCSTPRD)) + 
    coord_equal(xlim = c(bb[1,1] - 5, bb[1,2] + 5), 
                ylim = c(bb[2,1] - 5, bb[2,2] + 5)))
```

### Point track

Each forecast position is also included in the points dataframe. You can access this in the `df.gis_adv$al182012.018_5day_pts@data` object; no conversion necessary.

```{r, fig.width = 7, fig.align = "center"}
(p3 <- p2 + geom_point(data = df.gis_adv$al182012.018_5day_pts@data, aes(LON, LAT)))
```

## Forecast Cone

Forecast cone data is contained in the polygon dataset. To deal with this dataset you can use the `shp_to_df` function again or take the slightly longer way:

```{r}
fcst_cone <- df.gis_adv$al182012.018_5day_pgn
fcst_cone@data$id <- rownames(fcst_cone@data)
fcst_cone.points <- broom::tidy(fcst_cone, region = "id")
fcst_cone <- dplyr::left_join(fcst_cone.points, fcst_cone@data, by = "id")
```

```{r, fig.width = 7, fig.align = "center"}
p3 + geom_polygon(data = fcst_cone %>% filter(FCSTPRD == 120), 
                 aes(long, lat, group = group, fill = factor(FCSTPRD)), 
                 alpha = 0.5) + 
    geom_polygon(data = fcst_cone %>% filter(FCSTPRD == 72), 
                 aes(long, lat, group = group, fill = factor(FCSTPRD)), 
                 alpha = 0.5)
```

Note that in some GIS packages the forecast cones may be identical (though they shouldn't be). I've noticed it with Hurricanes Ike and Matthew; it's in the raw dataset.

## Watches and Warnings

You can also plot watches and warnings if any are in effect for the package release time. 

```{r}
ww_line <- shp_to_df(df.gis_adv$al182012.018_ww_wwlin)
```


```{r, fig.width = 7, fig.align = "center"}
p3 + geom_path(data = ww_line, 
               aes(x = long, y = lat, group = group, color = TCWW), size = 1)
```

In the example above you can see tropical storm warnings issued for the Bahamas, North Carolina and portions of the South Carolina and Florida coast while tropical storm watches are in effect for northern Florida and South Carolina.

## Probabilistic Storm Surge

The Tropical Cyclone Storm Surge Probabilities data shows the probability, in percent, of a specified storm surge occurring during the forecast period indicated. The product is based upon an ensemble of Sea, Lake, and Overland Surge from Hurricanes (SLOSH) model runs using the National Hurricane Center (NHC) official advisory and accounts for track, size, and intensity errors based on historical errors.

```{r}
gis_prob_storm_surge(key = "AL142016", products = list(psurge = 0), 
                     datetime = "20161006")
```

```{r, eval = FALSE}
df.gis_storm_surge <- gis_prob_storm_surge(key = "AL142016", 
                           products = list(psurge = 0), 
                           datetime = "20161006") %>% 
    last() %>% 
    gis_download()
```

```{r}
prob_surge <- shp_to_df(df.gis_storm_surge$al142016_2016100618_gt0)
bb <- sp::bbox(df.gis_storm_surge$al142016_2016100618_gt0)
```

```{r, fig.width = 7, fig.align = "center"}
p + geom_path(data = prob_surge, 
               aes(x = long, y = lat, group = group, color = PSurge00c), 
              size = 1, alpha = 0.5) + 
    coord_equal(xlim = c(bb[1,1], bb[1,2]), ylim = c(bb[2,1], bb[2,2]))
```

## Current and Forecast Wind Field Radii

Wind radii data may also be available for some cyclone advisory packages. This is the radius to which a minimum sustained wind speed may be felt from the center of circulation.

```{r}
gis_windfield("AL142016", advisory = "33")
```

```{r, eval = FALSE}
df.gis_wind_radii <- gis_windfield("AL142016", advisory = "33") %>% gis_download()
```

```{r}
wf_init <- shp_to_df(df.gis_wind_radii$al142016_2016100606_initialradii)
bb <- sp::bbox(df.gis_wind_radii$al142016_2016100606_forecastradii)
```

```{r, fig.width = 7, fig.align = "center"}
(p4 <- p + geom_polygon(data = wf_init, 
                  aes(x = long, y = lat, fill = factor(RADII)), alpha = 0.5) + 
    coord_equal(xlim = c(bb[1,1], bb[1,2]), ylim = c(bb[2,1], bb[2,2])))
```

Additionally, forecast wind radii data is also generally available in some packages

```{r}
wf_fcst <- shp_to_df(df.gis_wind_radii$al142016_2016100606_forecastradii)
```

```{r, fig.width = 7, fig.align = "center"}
p4 + geom_polygon(data = wf_fcst, 
                  aes(x = long, y = lat, group = group, fill = factor(RADII)), 
                  alpha = 0.5)
```

## Wind Speed Probabilities

Wind Speed probabilities show the chance of experiencing a minimum-sustained winds of 34, 50 and 64 knots with a given period of time (typically, 120 hours). These products are not storm-specific but are global so other active cyclones in other basins may also appear.

```{r}
gis_wsp(datetime = "2016100606", res = 0.5)
```

```{r, eval = FALSE}
df.gis_wsp <- 
  gis_download(
    "https://www.nhc.noaa.gov/gis/forecast/archive/2016100606_wsp_120hrhalfDeg.zip"
  )
```

### Cumulative Probability for >34kt Winds

```{r, fig.width = 7, fig.align = "center", error=T}
bb <- sp::bbox(df.gis_wsp$`2016100606_wsp34knt120hr_halfDeg`)
p + 
  geom_sf(
    data = st_as_sf(df.gis_wsp$`2016100606_wsp34knt120hr_halfDeg`), 
    aes(color = PWIND120)
  ) +
  coord_sf(xlim = c(bb[1,1], bb[1,2]), ylim = c(bb[2,1], bb[2,2]))
```

Cumulative wind speed probability for >50kt winds:

```{r, fig.width = 7, fig.align = "center", error=T}
bb <- sp::bbox(df.gis_wsp$`2016100606_wsp50knt120hr_halfDeg`)
p + 
  geom_sf(
    data = st_as_sf(df.gis_wsp$`2016100606_wsp50knt120hr_halfDeg`), 
    aes(color = PWIND120)
  ) +
  coord_sf(xlim = c(bb[1,1], bb[1,2]), ylim = c(bb[2,1], bb[2,2]))
```

Cumulative probability for >64kt winds:

```{r, fig.width = 7, fig.align = "center", error=T}
bb <- sp::bbox(df.gis_wsp$`2016100606_wsp64knt120hr_halfDeg`)
p + 
  geom_sf(
    data = st_as_sf(df.gis_wsp$`2016100606_wsp64knt120hr_halfDeg`), 
    aes(color = PWIND120)
  ) +
  coord_sf(xlim = c(bb[1,1], bb[1,2]), ylim = c(bb[2,1], bb[2,2]))
```
---
title: "Installing rrricanesdata"
author: "Tim Trice"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Installing rrricanesrdata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo = FALSE, message = FALSE, warning = FALSE}
library(rrricanesdata)
library(rrricanes)
```

With `rrricanes` you can access current and archived advisories as-issued from the National Hurricane Center archives. `rrricanesdata` is a complimentary data package meant to make it faster to get most data avialable from these storms. 

Installing `rrricanesdata` is simple:

```{r, eval = FALSE}
install.packages("rrricanesdata", 
                 repos = "https://timtrice.github.io/drat/", 
                 type = "source")
```

Data in `rrricanesdata` will be updated on the first of every month with a cutoff date of midnight on the last day of the month. So, advisories issued at any time during the current month will not be available; you will need to use any of `rrricanes` `get_*` functions. 

## Datasets

### `adv`
  * Key: Unique identifier of cyclone
  * Adv: Advisory number
  * Date: Date and time of advisory
  * Status: Classification of cyclone
  * Name: Name of cyclone
  * Lat: Latitude of cyclone center
  * Lon: Longitude of cyclone center
  * Wind: Maximum sustained one-minute winds in knots
  * Gust: Maximum sustained one-minute gusts in knots
  * Pressure: Minimum central pressure in millibars
  * PosAcc: Position accuracy of cyclone in nautical miles
  * FwdDir: Compass angle of forward motion
  * FwdSpeed: Forward speed in miles per hour
  * Eye: Size of eye in nautical miles
  * SeasNE: Radius of 12ft seas in northeast quadrant
  * SeasSE: Radius of 12ft seas in southeast quadrant
  * SeasSW: Radius of 12ft seas in southwest quadrant
  * SeasNW: Radius of 12ft seas in northwest quadrant

### `discus`
  * Status: Classification of storm, e.g., Tropical Storm, Hurricane,
    etc.
  * Name: Name of storm
  * Adv: Advisory Number
  * Date: Date of advisory issuance
  * Key: ID of cyclone
  * Contents: Text content of product

### `fcst`
  * Key: Unique identifier of cyclone
  * Adv: Advisory number
  * Date: Date and time of advisory
  * FcstDate: Forecast date and time in UTC
  * Lat: Forecast latitude
  * Lon: Forecast Longitude
  * Wind: Forecast wind in knots
  * Gust: Forecast gust in knots

### `fcst_wr`
  * Key: Unique identifier of cyclone
  * Adv: Advisory number
  * Date: Date and time of advisory
  * FcstDate: Forecast date and time in UTC
  * WindField: Minimum sustained wind field for quadrants
  * NE: Radius in nautical miles for northeast quadrant
  * SE: Radius in nautical miles for southeast quadrant
  * SW: Radius in nautical miles for southwest quadrant
  * NW: Radius in nautical miles for northwest quadrant

### `fstadv`
  * Status: Classification of cyclone
  * Name: Name of cyclone
  * Adv: Advisory number
  * Date: Date and time of advisory
  * Key: Unique identifier of cyclone
  * Lat: Latitude of cyclone center
  * Lon: Longitude of cyclone center
  * Wind: Maximum sustained one-minute winds in knots
  * Gust: Maximum sustained one-minute gusts in knots
  * Pressure: Minimum central pressure in millibars
  * PosAcc: Position accuracy of cyclone in nautical miles
  * FwdDir: Compass angle of forward motion
  * FwdSpeed: Forward speed in miles per hour
  * Eye: Size of eye in nautical miles
  * NE64: Radius of >=64kt winds in northeast quadrant
  * SE64: Radius of >=64kt winds in southeast quadrant
  * SW64: Radius of >=64kt winds in southwest quadrant
  * NW64: Radius of >=64kt winds in northwest quadrant
  * NE50: Radius of >=50kt winds in northeast quadrant
  * SE50: Radius of >=50kt winds in southeast quadrant
  * SW50: Radius of >=50kt winds in southwest quadrant
  * NW50: Radius of >=50kt winds in northwest quadrant
  * NE34: Radius of >=34kt winds in northwest quadrant
  * SE34: Radius of >=34kt winds in southeast quadrant
  * SW34: Radius of >=34kt winds in southwest quadrant
  * NW34: Radius of >=34kt winds in northwest quadrant
  * Hr\{n\}FcstDate: Forecast valid date
  * Hr\{n\}Lat: Forecast latitude in `n` hours
  * Hr\{n\}Lon: Forecast longitude in `n` hours
  * Hr\{n\}Wind: Forecast maximum wind in `n` hours
  * Hr\{n\}Gust: Forecast maximum gust in `n` hours
  * Hr\{n\}NE64: Forecast wind radius in `n` hours
  * Hr\{n\}SE64: Forecast wind radius in `n` hours
  * Hr\{n\}SW64: Forecast wind radius in `n` hours
  * Hr\{n\}NW64: Forecast wind radius in `n` hours
  * Hr\{n\}NE50: Forecast wind radius in `n` hours
  * Hr\{n\}SE50: Forecast wind radius in `n` hours
  * Hr\{n\}SW50: Forecast wind radius in `n` hours
  * Hr\{n\}NW50: Forecast wind radius in `n` hours
  * Hr\{n\}NE34: Forecast wind radius in `n` hours
  * Hr\{n\}SE34: Forecast wind radius in `n` hours
  * Hr\{n\}SW34: Forecast wind radius in `n` hours
  * Hr\{n\}NW34: Forecast wind radius in `n` hours
  * SeasNE: Radius of 12ft seas in northeast quadrant
  * SeasSE: Radius of 12ft seas in southeast quadrant
  * SeasSW: Radius of 12ft seas in southwest quadrant
  * SeasNW: Radius of 12ft seas in northwest quadrant

### `posest`
  * Status: Classification of storm, e.g., Tropical Storm, Hurricane,
    etc.
  * Name: Name of storm
  * Date: Date of advisory issuance
  * Contents: Text content of product

### `prblty`
  * Status: Classification of storm, e.g., Tropical Storm, Hurricane,
    etc.
  * Name: Name of storm
  * Adv: Advisory Number
  * Date: Date of advisory issuance
  * Location: Location for which the probability statistics rely
  * A: Probability of a strike within the next 12 hours
  * B: Probability of a strike between 12 and 24 hours
  * C: Probability of a strike between 24 and 36 hours
  * D: Probability of a strike between 36 and 48 hours
  * E: Probability of a strike between 48 and 72 hours

### `public`
  * Status: Classification of storm, e.g., Tropical Storm, Hurricane,
    etc.
  * Name: Name of storm
  * Adv: Advisory Number
  * Date: Date of advisory issuance
  * Key: Unique ID of the cyclone
  * Contents: Text content of product

### `storms`
  * Key: Storm ID
  * Name: Storm name
  * Wind: Peak wind speed in knots
  * StartDate: Date/time of first advisory
  * EndDate: Date/time of last advisory

### `update`
  * Status: Classification of storm, e.g., Tropical Storm, Hurricane,
    etc.
  * Name: Name of storm
  * Date: Date of advisory issuance
  * Key: Unique ID of cyclone
  * Contents: Text content of product

### `wndprb`
  * Status: Classification of storm, e.g., Tropical Storm, Hurricane,
    etc.
  * Name: Name of storm
  * Adv: Advisory Number
  * Date: Date of advisory issuance
  * Wind: Minimum wind speed for which probabilities reference
  * Wind12: Probability of sustained `Wind` within 12 hours
  * Wind24: Probability of sustained `Wind` within 24 hours
  * Wind24Cum: Cumulative probability through 24 hours
  * Wind36: Probability of sustained `Wind` within 36 hours
  * Wind36Cum: Cumulative probability through 36 hours
  * Wind48: Probability of sustained `Wind` within 48 hours
  * Wind48Cum: Cumulative probability through 48 hours
  * Wind72: Probability of sustained `Wind` within 72 hours
  * Wind72Cum: Cumulative probability through 72 hours
  * Wind96: Probability of sustained `Wind` within 96 hours
  * Wind96Cum: Cumulative probability through 96 hours
  * Wind120: Probability of sustained `Wind` within 120 hours
  * Wind120Cum: Cumulative probability through 120 hours

### `wr`
  * Key: Unique identifier of cyclone
  * Adv: Advisory number
  * Date: Date and time of advisory
  * Windfield: Minimum wind speed expected
  * NE: Radius of `Windfield` in the northeast quadrant
  * SE: Radius of `Windfield` in the southeast quadrant
  * SW: Radius of `Windfield` in the southwest quadrant
  * NW: Radius of `Windfield` in the northwest quadrant
---
title: "Wind Speed Probabilities"
author: "Tim Trice"
date: "June 11, 2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      fig.width = 10, 
                      fig.asp = 1, 
                      fig.align = "center")
```

```{r, message = FALSE}
library(dplyr)
library(ggplot2)
library(rrricanes)
library(rrricanesdata)
library(sp)
```

```{r}
key <- "AL092008"
adv <- 42
```

```{r}
wndprb <- wndprb %>% filter(Key == key, Adv <= adv)
```

### GIS Advisory Forecast Track, Cone of Uncertainty, and Watches/Warnings

```{r}
gis_adv <- gis_advisory(key = key, advisory = adv) %>% gis_download()
```

Get bounding box of the forecast polygon.

```{r}
bbox <- bbox(gis_adv$al092008.042_5day_pgn)
bbox
```

## Build a Tracking Chart

Generate a base plot of the Atlantic ocean.

```{r}
bp <- al_tracking_chart(color = "black", fill = "white", size = 0.1, res = 50)
bp
```

I like to add a little cushion for the map inset and forecast cone data. 

```{r}
lat_min <- bbox[2,1] - 5
lat_max <- bbox[2,2] + 5
lon_min <- bbox[1,1] - 10
lon_max <- bbox[1,2] + 10
```

Build a thin tracking map for the inset.

```{r}
bp_inset <- ggplotGrob(bp +
                           geom_rect(mapping = aes(xmin = lon_min, xmax = lon_max,
                                                   ymin = lat_min, ymax = lat_max),
                                     color = "red", alpha = 0) +
                           theme_bw() +
                           theme(axis.title = element_blank(),
                                 axis.ticks = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_blank(),
                                 plot.margin = margin(0, 0, 0, 0, "pt")))
```

Modify original `bp` zoomed in on our area of interest.

```{r}
bp <- bp +
    coord_equal(xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Lon",
         y = "Lat",
         caption = sprintf("rrricanes %s", packageVersion("rrricanes")))
bp
```

Combine `bp` and `bp_inset` to finalize initial base plot. `bp` will be a base plot without the inset. `bpi` will have the inset.

```{r}
bpi <- bp + annotation_custom(grob = bp_inset, xmin = lon_max - 5,
                              xmax = lon_max - 1, ymin = -Inf,
                              ymax = lat_min + 5)
bpi
```

The `wndprb` will not have coordinates for cities. An option is `al_prblty_stations`. However, please note this function may become [deprecated](https://github.com/ropensci/rrricanes/issues/46). 

```{r}
wndprb <- 
    wndprb %>% 
    left_join(al_prblty_stations(), by = "Location") %>% 
    mutate_at(.vars = c("Lat", "Lon"), .funs = as.numeric)
```

Check `wndprb` for NA values in `Lat`, `Lon`.

```{r}
any(is.na(wndprb$Lat), is.na(wndprb$Lon))
```

```{r}
wndprb_adv42 <- filter(wndprb, Adv == adv, Wind >= 64)

bpi +
    geom_point(
        data = wndprb_adv42, 
        aes(
            x = Lon,
            y = Lat, 
            color = Wind120Cum,
            size = Wind120Cum
        )
    ) + 
    scale_color_gradientn(colors = terrain.colors(10)) + 
    guides(size = FALSE) + 
    theme(
        legend.position = "bottom", 
        legend.box = "vertical"
    ) + 
    labs(title = "Total Probability of Wind >= 64kts within 120 Hours")
```
---
title: "Accumulated Cyclone Energy (ACE)"
author: "Tim Trice"
date: "June 27, 2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      fig.align = "center")
```

```{r, message = FALSE}
library(dplyr)
library(ggplot2)
library(HURDAT)
library(lubridate)
library(readr)
library(rrricanes)
library(rrricanesdata)
```

ACE or Accumulated Cyclone Energy is a method of measuring energy of a cyclone or for an entire season. It is calculated by the formula

$$ \text{ACE} = 10^{-4}\sum{v^2_\text{max}} $$

where $v_\text{max}$ is the wind speed in knots. Values may only be used when a storm is a tropical system with winds of at least 35 knots. Additionally, only six-hour intervals are used.

To calculate ACE you would want to use the `fstadv` dataset and apply the following rules:

* since forecast/advisory products are typically issued at 03:00, 09:00, 15:00 and 21:00 UTC filter out odd hours
* `Status` is Tropical Storm or Hurricane.
* `Wind` is not NA
* group by `Key`
* select `Wind`

```{r}
fstadv <- fstadv %>% 
    filter(hour(Date) %in% c(3, 9, 15, 21), 
           Status %in% c("Tropical Storm", "Hurricane"), 
           !is.na(Wind)) %>% 
    group_by(Key) %>% 
    select(Name, Wind)
```

Now let's summarise our dataset with new variable `ACE`.

```{r}
fstadv %>% 
    summarise(Name = last(Name), 
              ACE = sum(Wind^2) * 1e-04) %>% 
    arrange(desc(ACE)) %>% 
    top_n(10)
```

This matches somewhat well with [Wikipedia](https://en.wikipedia.org/wiki/Accumulated_cyclone_energy#Individual_storms) and other sources. But, you may notice we're missing some storms. `rrricanes` currently only holds data back to 1998; this data is considered "real-time". 

A companion package, [HURDAT](https://github.com/timtrice/HURDAT) is available in CRAN that has data for all cyclones dating back as far as 1851. This package has less data than `rrricanes`. But, as it is based on a post-storm reanalysis project, the data is more accurate. 

Let's revisit the top 10 using `HURDAT`:

```{r}
AL %>% 
    filter(hour(DateTime) %in% c(0, 6, 12, 18), 
           Status %in% c("TS", "HU"), 
           !is.na(Wind)) %>% 
    group_by(Key) %>% 
    summarise(Name = last(Name), 
              ACE = sum(Wind^2) * 1e-04) %>% 
    arrange(desc(ACE)) %>% 
    top_n(10)
```

A couple of things to notice here:

1. in `HURDAT`, the common times used are 00:00, 06:00, 12:00 and 18:00 UTC
2. Our list is more comprehensive than the Wikipedia list as that list only measures storms after 1950. 

`ACE` is slightly higher and that could be for a number of reasons. For example, on re-analysis the Hurricane Research Division may have determined a cyclone was actually tropical (shown in `HURDAT`) when initially it was believed to be extratropical (as shown in `rrricanes`). Or, and more likely, they determined through additional data that a storm was actually stronger than originally though. 

You can also calculate `ACE` for a season. Instead of grouping by `Key` we group by `Year`. I'll stick with `HURDAT` in this example.

```{r}
(df <- AL %>% 
    mutate(Year = year(DateTime)) %>% 
    filter(hour(DateTime) %in% c(0, 6, 12, 18), 
           Status %in% c("TS", "HU"), 
           !is.na(Wind)) %>% 
    group_by(Year) %>% 
    summarise(ACE = sum(Wind^2) * 1e-04) %>% 
    arrange(desc(ACE))) %>% 
    top_n(10)
```

This also matches relatively well with that on [Wikipedia](https://en.wikipedia.org/wiki/Accumulated_cyclone_energy#Atlantic_hurricane_seasons.2C_1950.E2.80.932017) and other sources. 

```{r}
ggplot(df, aes(x = Year, y = ACE)) + 
    geom_bar(stat = "identity") + 
    theme_bw()
```

It would certainly seem that tropical cyclone activity ebbs and flows over time. 
---
title: "Probabilistic Storm Surge"
author: "Tim Trice"
date: "June 11, 2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      fig.width = 10, 
                      fig.asp = 1, 
                      fig.align = "center")
```

```{r, message = FALSE}
library(dplyr)
library(ggplot2)
library(rrricanes)
library(rrricanesdata)
library(sp)
```

## Data

```{r}
key <- "AL092008"
adv <- 42
```

```{r}
fstadv <- fstadv %>% filter(Key == key, Adv <= adv)
```

## Probabilistic Storm Surge

Calculate the last date/time of forecast/advisory product and subtract 3 hrs

```{r}
dt <- (last(fstadv$Date) - (60 * 60 * 3)) %>% 
    strftime(format = "%Y%m%d%H", tz = "UTC")
```

Wrap `gis_download` in `safely`. Not all products will exist for every storm or even every advisory. 

I'm downloading the `psurge` products for 5 feet since there are no other storm surge products available for Ike.

```{r}
dl <- purrr::safely(.f = gis_download)
gis_surge <- gis_prob_storm_surge(key, products = list("psurge" = c(5)), 
                                  datetime = dt) %>% dl()

if (!is.null(gis_surge$error))
    message(gis_surge$error)
```

Generate a base plot of the Atlantic ocean.

```{r}
bp <- al_tracking_chart(color = "black", fill = "white", size = 0.1, res = 50)
```

Since we're dealing with a polygon shapefile, we can get the bounding box of the dataset.

```{r}
bbox <- bbox(gis_surge$result$al092008_2008091112_gt5)
```

Add a little cushion for the map inset.

```{r}
lat_min <- bbox[2,1] - 2
lat_max <- bbox[2,2] + 2
lon_min <- bbox[1,1] - 4
lon_max <- bbox[1,2] + 4
```

Build a map inset.

```{r}
bp_inset <- ggplotGrob(bp +
    geom_rect(mapping = aes(xmin = lon_min, xmax = lon_max,
                            ymin = lat_min, ymax = lat_max),
              color = "red", alpha = 0) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt")))
```

Modify original `bp` zoomed in on our area of interest.

```{r}
bp <- bp +
    coord_equal(xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Lon",
         y = "Lat",
         caption = sprintf("rrricanes %s", packageVersion("rrricanes")))
```

Combine `bp` and `bp_inset` to finalize initial base plot. `bp` will be a base plot without the inset. `bpi` will have the inset.

```{r}
bpi <- bp + annotation_custom(grob = bp_inset, xmin = lon_max - 5,
                             xmax = lon_max - 1, ymin = -Inf,
                             ymax = lat_min + 5)
```

Convert the SpatialPolygonsDataframe to a dataframe.

```{r}
shp_storm_surge <- shp_to_df(gis_surge$result$al092008_2008091112_gt5) 
```

Probability of storm surge greater than five feet.

```{r}
bpi + geom_point(data = shp_storm_surge, 
                aes(x = long, y = lat, color = ProbSurge05), size = 1) + 
    theme(legend.position = "bottom", 
          legend.box = "vertical") +
    labs(title = "Probabilistic Storm Surge > 5ft", 
         caption = sprintf("rrricanes %s", packageVersion("rrricanes")))
```
---
title: "Getting Started"
author: "Tim Trice"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting Started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, message = FALSE}
library(dplyr)
library(rrricanes)
```

## General Introduction

`rrricanes` is intended to give easy access to hurricane archives. It is a web-scraping tool that parses the National Hurricane Center's (NHC) archives to get storm data. Data is available for storms dating back to 1998.

There are two basins which data is available: north Atlantic ("AL") and northeastern Pacific ("EP"). The northeastern Pacific basin typically covers from the west coast of North America to -140&deg; longitude (140&deg;W). 

## Get Storms

By default, `get_storms` will return all storms that have developed for the current year in both basins. If no storms have developed, an error will be generated. For this example, we'll use 2012.

```{r get.storms, eval = FALSE}
df.al_2012 <- get_storms(years = 2012, basins = "AL")
```

## Getting Storm Data

`get_storm_data` can be used to retrieve one or multiple products for one or more cyclones. A list of dataframes is returned.

```{r, eval = FALSE}
df.al_18_2012_fstadv <- df.al_2012 %>% 
    filter(Name == "Hurricane Sandy") %>% 
    .$Link %>% 
    get_storm_data(products = "fstadv")
```

We can get the forecast/advisory data and wind speed probabilities at once:

```{r, eval = FALSE}
df.al_18_2012 <- df.al_2012 %>% 
    filter(Name == "Hurricane Sandy") %>% 
    .$Link %>% 
    get_storm_data(c("fstadv", "wndprb"))
```

`df.al_18_2012` now contains two dataframes for Hurricane Sandy; `fstadv` and `wndprb`.

## Forecast/Advisory Product (`fstadv`)

The core of a storm's dataset is located in the Forecast/Advisory product, `fstadv`. This product contains current location, forecast position, movement and structural details of the cyclone. 

To access only this product, we can use `get_fstadv`:

```{r, eval = FALSE}
df.al_18_2012_fstadv <- df.al_2012 %>% 
    filter(Name == "Hurricane Sandy") %>% 
    .$Link %>% 
    get_fstadv()
```

As you may have noticed above, the dataframe is very wide at 149 variables. There are four groups of variables in this dataset: current details, current wind radii, forecast positions, and forecast wind radii.

### Current Details

Let's look at an example of the current details.

```{r}
str(df.al_18_2012_fstadv %>% select(Status:Eye, SeasNE:SeasNW))
```

The most important variable in this dataset is `Key`. `Key` is a unique identifier for each storm that develops in either basin. It is formatted such as "AABBCCCC" where "AA" is the basin abbreviation (AL or EP), "BB" is the year number of the storm left-padded, and "CC" is the year of the storm.

`Adv` is the second-most important variable here. You'll notice it is in character format. For regularly-scheduled advisories, advisory numbers are always numeric. However, when watches and warnings are in effect, intermediate advisories are issued which are given alpha suffixes; i.e., 1, 2, 3, 3A, 4, 4A, 4B, 5, etc.

Only the Public Advisory (`public`) will be issued more frequently. All other regular products (`discus`, `fstadv`, `prblty`, `wndprb`) are generally issued every six hours.

`Status` lists the current designation of the cyclone, i.e., Tropical Depression, Tropical Storm, etc. A `Name` is given once a storm crosses the threshold of Tropical Storm; that is, winds greater than 33kts. 

`Lat` and `Lon` are the current position of the storm within `PosAcc` nautical miles. All distance measurements are in nautical miles.

`Wind` and `Gust` are current one-minute sustained wind speeds in knots (kts). You can use the function `knots_to_mph` to convert this. All wind speed values are in knots.

`Pressure` is the lowest atmospheric pressure of the cyclone either measured or estimated. It's value is in millibars but you can use `mb_to_in()` to convert to inches.

`FwdDir` and `FwdSpeed` show the compass direction of the forward movement of the cyclone. NA values indicate the storm is stationary or drifting. `FwdSpeed` is measured in knots.

In some cases, where hurricanes have an identifiable `Eye`, it's diameter in nautical miles will also be listed. 

Lastly, the `Seas` variables will exist for a storm of at least tropical storm-strength. This is the distance from the center of circulation that 12ft seas can be found in each quadrant. The measurement is in nautical miles.

Helper function `tidy_fstadv` will subset this data to a narrow dataframe.

```{r}
tidy_fstadv(df.al_18_2012_fstadv)
```

### Wind Radius

Any cyclone of at least tropical storm-strength will have associated wind radius values. This is the distance from the center of circulation that a specified wind speed (34kts, 50kts, 64kts) can be found in each quadrant. Measurement is in nautical miles.

```{r}
str(df.al_18_2012_fstadv %>% select(NE64:NW34))
```

A helper function, `tidy_wr` will reorganize this data into a narrow format and tidied up. Complete wind radius values that are NA are removed for efficiency.

```{r}
tidy_wr(df.al_18_2012_fstadv)
```

### Forecast

Most Forecast/Advisory products will have forecast data associated with it unless the storm has dissipated or is no longer tropical. There may be up to seven forecast positions. These positions are issued by 12-hour intervals through 48 hours where they are then at 24-hour intervals; 12, 24, 36, 48, 72, 96 and 120 hours.

```{r, eval = FALSE}
str(df.al_18_2012_fstadv %>% select(Hr12FcstDate:Hr12Gust))
```

Notice each variable begins with the prefix "Hrn" where n is the forecast period as noted above. Only Date, Lat, Lon, Wind, Gust and wind radius (will discuss shortly) are given for forecast periods.

Use `tidy_fcst` to tidy forecast data.

```{r}
tidy_fcst(df.al_18_2012_fstadv)
```

#### Forecast Dates/Times

A note about forecast times. 

```{r}
df.al_18_2012_fstadv %>% select(Date, Hr12FcstDate) %>% slice(1)
```

Notice the `Date` of this advisory is Oct 22 at 15:00 UTC. The `Hr12FcstDate` is Oct 23, 00:00 UTC. This difference, obviously, is not 12 hours. What gives? Forecast/Advisory products are issued with two "current" positions: one that is current (and provided in the dataset) and a position from three hours prior. So, in this specific advisory the text would contain the position of the storm for Oct 22, 12:00 UTC. It is from this position the forecast points are based. I do not know why.

Therefore, while officially the forecast periods are 12, 24, 36, ... hours, in reality they are 9, 21, 33, ... hours from the issuance time of the product. 

### Forecast Wind Radius

Some forecast positions may also contain wind radius information (only up to 72 hours). 

```{r, eval = FALSE}
str(df.al_18_2012_fstadv %>% select(Hr12NE64:Hr12NW34))
```

Again, these variables are prepended with the prefix prefix "Hrn" where n notes the forecast period. 

`tidy_fcst_wr` will tidy this subset of data.

```{r}
tidy_fcst_wr(df.al_18_2012_fstadv)
```

Please see the National Hurricane Center's website for more information on understanding the [Forecast/Advisory product](http://www.nhc.noaa.gov/help/tcm.shtml?).

## Strike Probabilities (`prblty`)

Strike probabilities were discontinued after the 2005 hurricane season (replaced by Wind Speed Probabilities; `wndprb`). For this example, we'll look at Hurricane Katrina. For this we use the function `get_prblty`.

```{r, eval = FALSE}
df.al_12_2005_prblty <- get_storms(year = 2005, basin = "AL") %>% 
    filter(Name == "Hurricane Katrina") %>% 
    .$Link %>% 
    get_prblty()
```

```{r}
str(df.al_12_2005_prblty)
```

This dataframe contains the possibility of a cyclone passing within 65 nautical miles of `Location`. The variables `A`, `B`, `C`, `D`, and `E` are as they appear in the products and were left as-is to avoid confusion. They're definition is as follows:

* `A` - current through 12 hours.
* `B` - within the next 12-24 hours
* `C` - within the next 24-36 hours
* `D` - within the next 36-48 hours
* `E` - Total probability from current through 48 hours. 

Many values in the text product may be "X" for less than 1% chance of a strike. These values are converted to 0 as the fields are numeric. 

The strike probability products did not contain `Key` which is the unique identifier for every cyclone. So the best way to do any joins will be by `Name`, `Adv` and `Date`.

Strike Probabilities may not exist for most Pacific cyclones.

## Wind Speed Probabilities (`wndprb`)

```{r, eval = FALSE}
df.al_18_2012_wndprb <- df.al_2012 %>% 
    filter(Name == "Hurricane Sandy") %>% 
    .$Link %>% 
    get_wndprb()
```

```{r}
str(df.al_18_2012_wndprb)
```

Wind Speed Probabilities are a bit more advanced than their predecessor. The `Wind` variable is for 34kt, 50kt and 64kt winds expected within a specific time period. 

Each consecutive variable is within a specific time-frame (12, 24, 36, 48, 72, 96 and 120 hours) for both that time frame and cumulative.

For example, `Wind24` is the chance of `Wind` between 12-24 hours. `Wind24Cum` is the cumulative probability from `Date` through 24 hours. 

As with strike probabilities, an "X" in the original text product meant less than 0.5% chance for the specified wind in the specified time period. "X" has been replaced by 0 in this package.

Wind Speed Probabilities may not exist for most Pacific cyclones.

See [Tropical Cyclone Wind Speed Probabilities Products](http://www.nhc.noaa.gov/about/pdf/About_Windspeed_Probabilities.pdf) for more information.

## Other products

Other products are available:

* `get_public` for Public Advisory statements. Think general information for the public audience. May not exist for some Pacific cyclones. Additionally, when watches and warnings are issued, these are issued every 3 hours (and, in some cases, every two).

* `get_discus` for Storm Discussions. These are more technical statements on the structure of a storm, forecast model tendencies and satellite presentation.

* `get_update` These are brief update statements when something considerable has changed in the cyclone or if the cyclone is making landfall.

* `get_posest`. Position estimates are generally issued when a storm is making landfall and may be issued hourly. 

Hurricane Ike, 2008, has both updates and position estimates. 

At this time none of these products are parsed. Only the content of the product is returned.

```{r keys_al}
keys_al <- keys[str_which(keys, "^AL.")]
```

```{r}
if (is_empty(keys_al)) {
  src <- "There are no storms in the Atlantic basin."
} else {
  src <- walk(keys_al, function(x) {
    knit_expand(file = "child_storm.Rmd", 
                arguments = list(key = x))
  })
}
```

`r knit(text = unlist(src))`
---
title: "Practice"
author: "Tim Trice"
date: "June 27, 2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      fig.align = "center")
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r libraries}
library(ggplot2)
library(HURDAT)
library(knitr)
library(purrr)
library(rrricanes)
library(rrricanesdata)
library(sp)
library(stringr)
library(tibble)
```

```{r gis_latest, message = FALSE}
gis <- flatten(gis_latest(verbose = FALSE))
```

```{r}
pts <- as_data_frame(gis$AL112017_pts)
lin <- shp_to_df(gis$AL112017_lin)
windswath <- shp_to_df(gis$AL112017_windswath)
radii <- shp_to_df(gis$AL112017_radii)
initial_radii <- shp_to_df(gis$al112017_2017090409_initialradii)
fcst_radii <- shp_to_df(gis$al112017_2017090409_forecastradii)
fcst_pts <- as_data_frame(gis$al112017_020_5day_pts)
fcst_lin <- shp_to_df(gis$al112017_020_5day_lin)
fcst_cone <- shp_to_df(gis$al112017_020_5day_pgn)
```

```{r}
l <- list("bbox.pts" = bbox(gis$AL112017_pts),
          "bbox.lin" = bbox(gis$AL112017_lin), 
          "bbox.windswath" = bbox(gis$AL112017_windswath), 
          "bbox.radii" = bbox(gis$AL112017_radii), 
          "bbox.initial_radii" = bbox(gis$al112017_2017090309_initialradii), 
          "bbox.fcst_radii" = bbox(gis$al112017_2017090309_forecastradii), 
          "bbox.fcst_pts" = bbox(gis$al112017_016_5day_pts), 
          "bbox.fcst_lin" = bbox(gis$al112017_016_5day_lin), 
          "bbox.fcst_cone" = bbox(gis$al112017_016_5day_pgn))

bb <- matrix(c(map(l, `[[`, 1) %>% flatten_dbl() %>% min(), 
               map(l, `[[`, 2) %>% flatten_dbl() %>% max()), 
             nrow = 2, 
             ncol = 2)
```

```{r}
# Current and past details
plot_lin <- geom_path(data = lin, aes(x = long, y = lat, color = STORMTYPE))

plot_pts <- geom_point(data = pts, 
             aes(x = LON, y = LAT, color = STORMTYPE, size = INTENSITY))

# Forecast details
plot_fcst_lin <- geom_path(data = fcst_lin, 
                           aes(x = long, y = lat, color = STORMTYPE))

plot_fcst_pts <- geom_point(data = fcst_pts, 
                            aes(x = LON, y = LAT, color = STORMTYPE, 
                                size = MAXWIND))
# Forecast cone
plot_fcst_cone <- geom_polygon(data = fcst_cone, 
                               aes(x = long, y = lat, group = group, 
                                   fill = FCSTPRD), alpha = 0.25)

# Wind radii details
plot_initial_radii <- geom_polygon(data = initial_radii,
                                   aes(x = long, y = lat, group = group, 
                                       fill = factor(RADII)), 
                                   alpha = 0.25)
```

```{r}
tracking_chart(color = "black", fill = "white", size = 0.1, res = 50) + 
  plot_lin + 
  plot_pts + 
  plot_fcst_lin + 
  plot_fcst_pts + 
  plot_fcst_cone + 
  coord_equal(xlim = c(-80, -16), 
              ylim = c(11.5, 25)) + 
  theme(legend.position = "bottom", 
        legend.direction = "vertical")
```


```{r keys_ep}
(keys_ep <- keys[str_which(keys, "^EP.")])
```

### {{arguments$key}}

```{r}
# Find names of relevant datasets
ds <- grep(sprintf("^%s.", arguments$key), names(gis), ignore.case = TRUE, 
                   value = TRUE)
```


```{r}
# Get previous track points
points <- as_data_frame(gis[[ds[grep(sprintf("%s_pts", arguments$key), 
                                     ds, ignore.case = TRUE)]]])
# Convert points$STORMTYPE to factor
points$STORMTYPE <- factor(points$STORMTYPE, 
                           levels = c("DB", "LO", "TD", "TS", "HU", "MH"), 
                           labels = c("Disturbance", "Low", 
                                      "Tropical Depression", "Tropical Storm", 
                                      "Hurricane", "Major Hurricane"))

plot_points <- geom_point(data = points, aes(x = LON, y = LAT, size = STORMTYPE))
```

```{r}
# Get forecast track points
fcst_points <- as_data_frame(gis[[ds[grep(sprintf("%s.+_5day_pts", 
                                                  arguments$key), ds, 
                                          ignore.case = TRUE)]]])


# Convert points$STORMTYPE to factor
fcst_points$STORMTYPE <- factor(fcst_points$STORMTYPE, 
                           levels = c("DB", "LO", "TD", "TS", "HU", "MH"), 
                           labels = c("Disturbance", "Low", 
                                      "Tropical Depression", "Tropical Storm", 
                                      "Hurricane", "Major Hurricane"))

plot_fcst_points <- geom_point(data = fcst_points, 
                               aes(x = LON, y = LAT, size = STORMTYPE))
```

```{r}
tracking_chart(color = "black", fill = "white", size = 0.1, res = 50) + 
  plot_points + 
  plot_fcst_points + 
  theme(legend.position = "bottom", 
        legend.box = "vertical")
```

---
title: "Latest Tropical Cyclone Activity"
output: 
  html_document:
    toc: TRUE
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r libraries}
library(ggplot2)
library(HURDAT)
library(knitr)
library(purrr)
library(rrricanes)
library(rrricanesdata)
library(stringr)
library(tibble)
```

```{r gis_latest, message = FALSE}
gis <- flatten(gis_latest(verbose = FALSE))
```

```{r keys}
# Keys of existing storms
keys <- str_extract(names(gis), "(^[[:alpha:]]{2}[[:digit:]]{6})") %>% 
  toupper() %>% 
  unique() %>% 
  .[complete.cases(.)]
```

## Atlantic Basin

```{r al, child = "child_al.Rmd"}
```

## East Pacific Basin

```{r ep, child = "child_ep.Rmd"}
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.al_18_2012}
\alias{df.al_18_2012}
\title{Forecast/Advisory and Wind Speed Probabilities for Hurricane Sandy (AL182012)}
\format{An object of class \code{list} of length 2.}
\source{
\url{http://www.nhc.noaa.gov/archive/2012/SANDY.shtml?}
}
\usage{
df.al_18_2012
}
\description{
Forecast/Advisory and Wind Speed Probabilities for Hurricane Sandy (AL182012)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storm_data.R
\name{extract_product_contents}
\alias{extract_product_contents}
\title{extract_product_contents}
\usage{
extract_product_contents(links, product)
}
\arguments{
\item{links}{URLs to storm products}

\item{product}{specific product to parse}
}
\description{
Get and parse product contents for each links
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{get_nhc_link}
\alias{get_nhc_link}
\title{get_nhc_link}
\usage{
get_nhc_link(withTrailingSlash = TRUE, protocol = "https")
}
\arguments{
\item{withTrailingSlash}{True, by default. False returns URL without
trailing slash.}

\item{protocol}{https or http}
}
\description{
Return root link of NHC archive pages.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.al_12_2005_prblty}
\alias{df.al_12_2005_prblty}
\title{Strike probabilities for Hurricane Katrina (AL122005)}
\format{An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 937 rows and 10 columns.}
\source{
\url{http://www.nhc.noaa.gov/archive/2005/KATRINA.shtml?}
}
\usage{
df.al_12_2005_prblty
}
\description{
Strike probabilities for Hurricane Katrina (AL122005)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrapers.R
\name{scrape_key}
\alias{scrape_key}
\title{scrape_key}
\usage{
scrape_key(header)
}
\arguments{
\item{header}{Header text of product.}
}
\description{
Extract Key from header
}
\seealso{
\code{\link{scrape_header}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tracking_chart.R
\name{al_tracking_chart}
\alias{al_tracking_chart}
\title{al_tracking_chart}
\usage{
al_tracking_chart(...)
}
\arguments{
\item{...}{Additional parameters for \link{tracking_chart} and ggplot2}
}
\value{
ggplot2 object centered on Atlantic basin.
}
\description{
Build tracking chart centered on Atlantic Basin.
}
\examples{
\dontrun{
# Build map with white land areas, thin black borders
al_tracking_chart(color = "black", size = 0.1, fill = "white")

# 50nm resolution, no states
al_tracking_chart(res = 50, states = FALSE, color = "black", size = 0.1,
          fill = "white")

# 50nm resolution, coastlines only
al_tracking_chart(countries = FALSE, res = 50, color = "black", size = 0.1,
          fill = "white")

# Adding and modifying with ggplot functions
al_tracking_chart(color = "black", size = 0.1, fill = "white") +
  ggplot2::labs(x = "Lon", y = "Lat",
  title = "Base Atlantic Tracking Chart")
}
}
\seealso{
\code{\link{tracking_chart}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_latest}
\alias{gis_latest}
\title{gis_latest}
\usage{
gis_latest(basins = c("AL", "EP"), ...)
}
\arguments{
\item{basins}{AL and/or EP.}

\item{...}{additional parameters for rgdal::readOGR}
}
\description{
Latest GIS datasets for \strong{active} cyclones
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_lat_lon}
\alias{fstadv_lat_lon}
\title{fstadv_lat_lon}
\usage{
fstadv_lat_lon(contents)
}
\arguments{
\item{contents}{text contents of FORECAST/ADVISORY}
}
\value{
numeric
}
\description{
Returns numeric for latitude or longitude; negative if in
  southern or eastern hemisphere
}
\details{
Helper function to take character latitude or longitude and,
depending on the value of hemisphere return a positive or negative numeric,
or NA if not found.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wndprb.R
\name{cp_prblty_stations}
\alias{cp_prblty_stations}
\title{cp_prblty_stations}
\usage{
cp_prblty_stations()
}
\description{
Retrieve list of probability stations based in the central
Pacific from the NHC. To be used in tandem with `wndprb` products.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_orig}
\alias{filter_orig}
\title{filter_orig}
\usage{
filter_orig(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
For Katrina, 2005, there are two identical discussions; one of
  which has 'orig' in the URL. Because this link is not captured above it
  will throw an error. This function is an effort to capture it and pass
  it through the validation but, at least for this Katrina it is not
  necessary. That being said, if there is output for any other storms it
  should be reviewed as it is common for the NHC to issue UPDATED statements.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_pos_accuracy}
\alias{fstadv_pos_accuracy}
\title{fstadv_pos_accuracy()}
\usage{
fstadv_pos_accuracy(contents)
}
\arguments{
\item{contents}{text contents of FORECAST/ADVISORY}
}
\value{
numeric
}
\description{
Get position accuracy
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public.R
\name{public}
\alias{public}
\title{public}
\usage{
public(contents)
}
\arguments{
\item{contents}{Link to a storm's specific public advisory product.}
}
\value{
Dataframe
}
\description{
Parse Public Advisory products
}
\details{
Given a direct link to a public advisory product, parse and return
dataframe of values.
}
\seealso{
\code{\link{get_public}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posest.R
\name{posest}
\alias{posest}
\title{posest}
\usage{
posest(contents)
}
\arguments{
\item{contents}{URL of a specific position estimate product}
}
\value{
Dataframe
}
\description{
Extrapolate data from Position Estimate products.
}
\details{
Given a direct link to a position estimate product, parse and return
dataframe of values.
}
\seealso{
\code{\link{get_posest}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{tidy_wr}
\alias{tidy_wr}
\title{tidy_wr}
\usage{
tidy_wr(df)
}
\arguments{
\item{df}{fstadv dataframe object}
}
\description{
Tidy current wind radius of a fstadv dataframe object.
}
\details{
Returns tidy dataframe of current wind radius values for a cyclone.
Returns only complete.cases (based on quadrants).
\describe{
 \item{Key}{Unique identifier of cyclone}
 \item{Adv}{Advisory number}
 \item{Date}{Date and time of advisory}
 \item{Windfield}{Minimum wind speed expected}
 \item{NE}{Radius of `Windfield` in the northeast quadrant}
 \item{SE}{Radius of `Windfield` in the southeast quadrant}
 \item{SW}{Radius of `Windfield` in the southwest quadrant}
 \item{NW}{Radius of `Windfield` in the northwest quadrant}
}
}
\examples{
\dontrun{
get_fstadv("http://www.nhc.noaa.gov/archive/1998/1998ALEXadv.html") \%>\%
  tidy_wr()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{saffir}
\alias{saffir}
\title{saffir}
\usage{
saffir(x)
}
\arguments{
\item{x}{Vector of wind speed values.}
}
\description{
Return category of tropical cyclone based on wind. Saffir-
Simpson Hurricane Scale does not apply to non-tropical cyclones.
}
\examples{
saffir(c(32, 45, 70, 90, 110, 125, 140))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{tidy_fcst_wr}
\alias{tidy_fcst_wr}
\title{tidy_fcst_wr}
\usage{
tidy_fcst_wr(df)
}
\arguments{
\item{df}{fstadv dataframe object}
}
\description{
Tidy forecast wind radii of a fstadv dataframe object
}
\details{
Tidies forecast wind radius for each forecast position. Complete
cases only (by quadrants). Use Key, Adv and Date to join with other tidy
dataframes.

\describe{
 \item{Key}{Unique identifier of cyclone}
 \item{Adv}{Advisory number}
 \item{Date}{Date and time of advisory}
 \item{FcstDate}{Forecast date and time in UTC}
 \item{WindField}{Minimum sustained wind field for quadrants}
 \item{NE}{Radius in nautical miles for northeast quadrant}
 \item{SE}{Radius in nautical miles for southeast quadrant}
 \item{SW}{Radius in nautical miles for southwest quadrant}
 \item{NW}{Radius in nautical miles for northwest quadrant}
}
}
\examples{
\dontrun{
get_fstadv("http://www.nhc.noaa.gov/archive/1998/1998ALEXadv.html") \%>\%
  tidy_fcst_wr()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{get_url_contents}
\alias{get_url_contents}
\title{get_url_contents}
\usage{
get_url_contents(links)
}
\arguments{
\item{link}{URL to download}
}
\description{
Get contents from URL
}
\details{
This function primarily is reserved for extracting the contents of
the individual products \(thought it can be used in other instances\). Often,
there are timeout issues. This is an attempt to try to work around that.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_public}
\alias{filter_public}
\title{filter_public}
\usage{
filter_public(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
Get public advisories
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_fwd_mvmt}
\alias{fstadv_fwd_mvmt}
\title{fstadv_fwd_mvmt}
\usage{
fstadv_fwd_mvmt(contents, what = NULL)
}
\arguments{
\item{contents}{text contents of FORECAST/ADVISORY}

\item{what}{is being retrieved
\itemize{
  \item fwd_dir integer azimuth direction of movement (0 - 360)
  \item fwd_speed integer speed of movement in kts
}}
}
\value{
numeric
}
\description{
Get forward movement direction and speed
}
\details{
If STATIONARY should return NA
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_breakpoints}
\alias{gis_breakpoints}
\title{gis_breakpoints}
\usage{
gis_breakpoints()
}
\description{
Return link to breakpoints shapefile by year
}
\details{
Coastal areas placed under tropical storm and hurricane watches and
  warnings are identified through the use of "breakpoints." A tropical
  cyclone breakpoint is defined as an agreed upon coastal location that can
  be chosen as one of two specific end points or designated places between
  which a tropical storm/hurricane watch/warning is in effect. The U.S.
  National Weather Service designates the locations along the U.S. East,
  Gulf, and California coasts, Puerto Rico, and Hawaii. These points are
  listed in NWS Directive 10-605 (PDF). Individual countries across the
  Caribbean, Central America, and South America provide coastal locations
  for their areas of responsibility to the U.S. National Weather Service for
  the National Hurricane Center's use in tropical cyclone advisories when
  watches/warnings are issued by international partners. The National
  Hurricane Center maintains a list of pre-arranged breakpoints for the U.S.
  Atlantic and Gulf coasts, Mexico, Cuba and the Bahamas. Other sites are
  unofficial and sites not on the list can be selected if conditions warrant.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_wsp}
\alias{gis_wsp}
\title{gis_wsp}
\usage{
gis_wsp(datetime, res = c(5, 0.5, 0.1))
}
\arguments{
\item{datetime}{Datetime in \%Y\%m\%d\%H format. \%m, \%d and \%H are
optional but will return more datasets.}

\item{res}{Resolution as a numeric vector; 5, 0.5, 0.1.}
}
\description{
Wind Speed Probabilities
}
\details{
Probability winds affecting an area within a forecast period.
Datasets contain windfields for 34kt, 50kt and 64kt. Resolution is at 5km,
0.5 degrees and 0.1 degrees. Not all resolutions may be available for all
storms. Not all windfields will be available for all advisories.
}
\examples{
\dontrun{
# Return datasets for January 1, 2016 with resolution of 0.5 degrees
gis_wsp("20160101", res = 0.5)

# Return wsp of 0.1 and 0.5 degree resolution, July, 2015
gis_wsp("201507", res = c(0.5, 0.1))
}
}
\seealso{
\code{\link{gis_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{tidy_fcst}
\alias{tidy_fcst}
\title{tidy_fcst}
\usage{
tidy_fcst(df)
}
\arguments{
\item{df}{fstadv dataframe object}
}
\description{
Tidy forecasts of a fstadv dataframe object.
}
\details{
Gathers all forecast points, tidies dataframe to make one row per
forecast position. Complete cases only. Use Key, Adv and Date to join with
other tidy dataframes.

\describe{
 \item{Key}{Unique identifier of cyclone}
 \item{Adv}{Advisory number}
 \item{Date}{Date and time of advisory}
 \item{FcstDate}{Forecast date and time in UTC}
 \item{Lat}{Forecast latitude}
 \item{Lon}{Forecast Longitude}
 \item{Wind}{Forecast wind in knots}
 \item{Gust}{Forecast gust in knots}
}
}
\examples{
\dontrun{
get_fstadv("http://www.nhc.noaa.gov/archive/1998/1998ALEXadv.html") \%>\%
  tidy_fcst()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discus.R
\name{get_discus}
\alias{get_discus}
\title{get_discus}
\usage{
get_discus(links)
}
\arguments{
\item{links}{URL to storm's archive page.}
}
\description{
Return dataframe of discussion data.
\describe{
  \item{Status}{Classification of storm, e.g., Tropical Storm, Hurricane,
  etc.}
  \item{Name}{Name of storm}
  \item{Adv}{Advisory Number}
  \item{Date}{Date of advisory issuance}
  \item{Key}{ID of cyclone}
  \item{Contents}{Text content of product}
}
}
\examples{
\dontrun{
# Return dataframe of storm discussions for Tropical Storm Alex (AL011998)
get_discus("http://www.nhc.noaa.gov/archive/1998/1998ALEXadv.html")
}
}
\seealso{
\code{\link{get_storms}}, \code{\link{public}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storm_list.R
\name{get_ftp_dirs}
\alias{get_ftp_dirs}
\title{get_ftp_dirs()}
\usage{
get_ftp_dirs(x)
}
\description{
Get a list of the FTP directors in /atcf/archive
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prblty.R
\name{get_prblty}
\alias{get_prblty}
\title{get_prblty}
\usage{
get_prblty(links)
}
\arguments{
\item{links}{URL to storm's archive page.}
}
\description{
Strike probabilities; the chances of the center of a cyclone
passing within 65 nautical miles of a location.
\describe{
  \item{Status}{Classification of storm, e.g., Tropical Storm, Hurricane,
  etc.}
  \item{Name}{Name of storm}
  \item{Adv}{Advisory Number}
  \item{Date}{Date of advisory issuance}
  \item{Location}{Location for which the probability statistics rely}
  \item{A}{Probability of a strike within the next 12 hours}
  \item{B}{Probability of a strike between 12 and 24 hours}
  \item{C}{Probability of a strike between 24 and 36 hours}
  \item{D}{Probability of a strike between 36 and 48 hours}
  \item{E}{Probability of a strike between 48 and 72 hours}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posest.R
\name{get_posest}
\alias{get_posest}
\title{get_posest}
\usage{
get_posest(links)
}
\arguments{
\item{links}{URL to storm's archive page.}
}
\description{
Return dataframe of position estimate data.
}
\details{
This product was discontinued after the 2013 hurricane season and is
now included in the Tropical Cyclone Update product (\code{\link{update}}).
\describe{
  \item{Status}{Classification of storm, e.g., Tropical Storm, Hurricane,
  etc.}
  \item{Name}{Name of storm}
  \item{Date}{Date of advisory issuance}
  \item{Contents}{Text content of product}
}
}
\seealso{
\code{\link{get_storms}}, \code{\link{posest}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_outlook}
\alias{gis_outlook}
\title{gis_outlook}
\usage{
gis_outlook()
}
\description{
Tropical Weather Outlook
}
\seealso{
\code{\link{gis_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_eye}
\alias{fstadv_eye}
\title{fstadv_eye}
\usage{
fstadv_eye(contents)
}
\arguments{
\item{contents}{text contents of FORECAST/ADVISORY}
}
\value{
numeric
}
\description{
Get eye diameter, if available
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storm_data.R
\name{extract_storm_links}
\alias{extract_storm_links}
\title{extract_storm_links}
\usage{
extract_storm_links(links)
}
\arguments{
\item{links}{URLs to a storm's archive page}
}
\description{
Extract product links from a storm's archive page
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.al_2012}
\alias{df.al_2012}
\title{Atlantic cyclones for 2012}
\format{An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 19 rows and 4 columns.}
\source{
\url{http://www.nhc.noaa.gov/archive/2012/}
}
\usage{
df.al_2012
}
\description{
Atlantic cyclones for 2012
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_prob_storm_surge}
\alias{gis_prob_storm_surge}
\title{gis_prob_storm_surge}
\usage{
gis_prob_storm_surge(key, products, datetime = NULL)
}
\arguments{
\item{key}{Key of storm (i.e., AL012008, EP092015)}

\item{products}{list of products and associated n values; psurge (0:20) or
esurge (10, 20, 30, 40, 50).}

\item{datetime}{Datetime in \%Y\%m\%d\%H format.}
}
\description{
Probabilistic Storm Surge
}
\details{
Probabilistic Storm Surge Forecasts
}
\section{Products}{

\describe{
  \item{esurge}{The Tropical Cyclone Storm Surge Exceedances (P-Surge 2.5)
    data shows the probability, in percent, of a specified storm surge,
    including tides, exceeding the specified height, in feet, during the
    forecast period indicated. The 10 percent exceedance height, for example,
    is the storm surge height, including tides, above ground level (AGL) such
    that there is a 10 percent chance of exceeding it. The product is based
    upon an ensemble of Sea, Lake,and Overland Surge from Hurricanes (SLOSH)
    model runs using the National Hurricane Center (NHC) official advisory
    and accounts for track, size, and intensity errors based on historical
    errors and astronomical tide. Valid values are 10, 20, 30, 40 or 50.}
  \item{psurge}{The Tropical Cyclone Storm Surge Probabilities (P-Surge 2.5)
    data shows the probability, in percent, of a specified storm surge
    occurring during the forecast period indicated. The product is based upon
    an ensemble of Sea, Lake, and Overland Surge from Hurricanes (SLOSH)
    model runs using the National Hurricane Center(NHC) official advisory
    and accounts for track, size, and intensity errors based on historical
    errors and astronomical tide. Valid values are 0:20.}
}
}

\examples{
\dontrun{
# Return the last psurge0 product for storm AL092016
gis_prob_storm_surge("AL092016", products = list("psurge" = 0))

# Return the psurge0 and esurge10 products for storm AL092016
gis_prob_storm_surge("AL092016", products = list("psurge" = 0, "esurge" = 10))

# Return all psurge0 products for Sep 2, 2016, storm AL092016
gis_prob_storm_surge("AL092016", products = list("psurge" = 0),
           datetime = "20160902")
}
}
\seealso{
\href{http://www.nhc.noaa.gov/surge/psurge.php}{Tropical Cyclone Storm Surge Probabilities}

\code{\link{gis_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storm_list.R
\name{get_ftp_storm_data}
\alias{get_ftp_storm_data}
\title{get_ftp_storm_data}
\usage{
get_ftp_storm_data(stormid, products = c("discus", "fstadv", "posest",
  "public", "prblty", "update", "wndprb"))
}
\arguments{
\item{stormid}{A six-character alphanumeric string formatted as AABBCCCC
where
\describe{
  \item{AA}{The basin of the storm; AL or EP}
  \item{BB}{Storm number for the year as decimal number
    (e.g., 01, 02, ..., 10, ...)}
  \item{CCCC}{Year with century)}
}}

\item{products}{Products to retrieve; discus, fstadv, posest, public,
prblty, update, and windprb.}
}
\description{
Retrieve text products from the National Hurricane Center's FTP
  server. Not all products may exist for certain storms.
}
\seealso{
\code{\link{get_storm_data}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv}
\alias{fstadv}
\title{fstadv}
\usage{
fstadv(contents)
}
\arguments{
\item{contents}{URL of a specific FORECAST/ADVISORY product}
}
\description{
Extrapolate data from FORECAST/ADVISORY products.
}
\details{
Given a direct link to a forecast/advisory product, parse and
return dataframe of values.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_seas}
\alias{fstadv_seas}
\title{fstadv_seas}
\usage{
fstadv_seas(content)
}
\arguments{
\item{content}{text of product}
}
\value{
boolean
}
\description{
There is only one line of sea data, 12FT seas in each quadrant.
So this should go easier than the wind fields
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{nm_to_sm}
\alias{nm_to_sm}
\title{nm_to_sm}
\usage{
nm_to_sm(x)
}
\arguments{
\item{x}{Nautical miles}
}
\description{
Convert nautical miles to survey miles
}
\examples{
nm_to_sm(c(50, 100, 150))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{get_nhc_ftp_link}
\alias{get_nhc_ftp_link}
\title{get_nhc_ftp_link}
\usage{
get_nhc_ftp_link(withTrailingSlash = TRUE)
}
\arguments{
\item{withTrailingSlash}{True, by default. False returns URL without
trailing slash.}
}
\description{
Return root of NHC FTP server
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tracking_chart.R
\name{tracking_chart}
\alias{tracking_chart}
\title{tracking_chart}
\usage{
tracking_chart(countries = TRUE, states = TRUE, res = 110, ...)
}
\arguments{
\item{countries}{Show country borders. Default TRUE.}

\item{states}{Show state boundaries. Default TRUE. Ignored if `countries` is
FALSE.}

\item{res}{Resolution of charts; 110 (1:110m), 50 (1:50m), 10 (1:10m).
Default is low. The higher the resolution, the longer the plot takes to
appear.}

\item{...}{Additional ggplot2::aes parameters}
}
\value{
Returns ggplot2 object that can be printed directly or have new
  layers added.
}
\description{
Build base tracking chart using ggplot
}
\examples{
\dontrun{
# Build map with white land areas, thin black borders
tracking_chart(color = "black", size = 0.1, fill = "white")

# 50nm resolution, no states
tracking_chart(res = 50, states = FALSE, color = "black", size = 0.1,
       fill = "white")

# 50nm resolution, coastlines only
tracking_chart(countries = FALSE, res = 50, color = "black", size = 0.1,
       fill = "white")

# Adding and modifying with ggplot functions
tracking_chart(color = "black", size = 0.1, fill = "white") +
  ggplot2::labs(x = "Lon", y = "Lat", title = "Base Tracking Chart")
}
}
\seealso{
\code{\link[ggplot2]{aes}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tracking_chart.R
\name{ep_tracking_chart}
\alias{ep_tracking_chart}
\title{ep_tracking_chart}
\usage{
ep_tracking_chart(...)
}
\arguments{
\item{...}{Additional parameters for ggplot2}
}
\value{
ggplot2 object centered on northeast Pacific basin.
}
\description{
Build tracking chart centered on northeast Pacific Basin.
}
\examples{
\dontrun{
# Build map with white land areas, thin black borders
ep_tracking_chart(color = "black", size = 0.1, fill = "white")

# 50nm resolution, no states
ep_tracking_chart(res = 50, states = FALSE, color = "black", size = 0.1,
          fill = "white")

# 50nm resolution, coastlines only
ep_tracking_chart(countries = FALSE, res = 50, color = "black", size = 0.1,
          fill = "white")

# Adding and modifying with ggplot functions
ep_tracking_chart(color = "black", size = 0.1, fill = "white") +
  ggplot2::labs(x = "Lon", y = "Lat",
  title = "Base East Pacific Tracking Chart")
}
}
\seealso{
\code{\link{tracking_chart}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.gis_storm_surge}
\alias{df.gis_storm_surge}
\title{GIS storm surge shapefile dataset for Hurricane Sandy (AL182012)}
\format{An object of class \code{list} of length 1.}
\source{
\url{http://www.nhc.noaa.gov/gis/archive_psurge_results.php?id=al18&year=2012&name=Hurricane\%20SANDY}
}
\usage{
df.gis_storm_surge
}
\description{
GIS storm surge shapefile dataset for Hurricane Sandy (AL182012)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{get_fstadv}
\alias{get_fstadv}
\title{get_fstadv}
\usage{
get_fstadv(links)
}
\arguments{
\item{links}{URL to storms' archive page.}
}
\description{
Return dataframe of forecast/advisory data.
}
\details{
Returns a wide dataframe of most the data available in a cyclones
forecast/advisory product (watches and warnings are not included at this
time).

Overall structure of the dataframe is listed below. Note the following
clarifications:

\enumerate{
  \item The value of `n` in `Hr\{n\}` variables is the forecast period.
    Up to 2002, forecast periods are 12, 24, 36, 48 and 72 hours. After
    2002, forecast periods were extended to 96 and 120 hours. Not all
    forecast periods will be available for every cyclone advisory (e.g.,
    if it is dissipating or expected to dissipate.)
  \item Wind radius data is not included 96 and 120 hour forecast periods.
  \item Forecast dates are not truly 12, 24, ..., 120 hours from the
    date/time of the advisory. The NHC issues two positions in these
    products; one for current and one for three hours prior. It is the
    latter position the forecast date/times are based.
}

\describe{
 \item{Status}{Classification of cyclone}
 \item{Name}{Name of cyclone}
 \item{Adv}{Advisory number}
 \item{Date}{Date and time of advisory}
 \item{Key}{Unique identifier of cyclone}
 \item{Lat}{Latitude of cyclone center}
 \item{Lon}{Longitude of cyclone center}
 \item{Wind}{Maximum sustained one-minute winds in knots}
 \item{Gust}{Maximum sustained one-minute gusts in knots}
 \item{Pressure}{Minimum central pressure in millibars}
 \item{PosAcc}{Position accuracy of cyclone in nautical miles}
 \item{FwdDir}{Compass angle of forward motion}
 \item{FwdSpeed}{Forward speed in miles per hour}
 \item{Eye}{Size of eye in nautical miles}
 \item{NE64}{Radius of >=64kt winds in northeast quadrant}
 \item{SE64}{Radius of >=64kt winds in southeast quadrant}
 \item{SW64}{Radius of >=64kt winds in southwest quadrant}
 \item{NW64}{Radius of >=64kt winds in northwest quadrant}
 \item{NE50}{Radius of >=50kt winds in northeast quadrant}
 \item{SE50}{Radius of >=50kt winds in southeast quadrant}
 \item{SW50}{Radius of >=50kt winds in southwest quadrant}
 \item{NW50}{Radius of >=50kt winds in northwest quadrant}
 \item{NE34}{Radius of >=34kt winds in northwest quadrant}
 \item{SE34}{Radius of >=34kt winds in southeast quadrant}
 \item{SW34}{Radius of >=34kt winds in southwest quadrant}
 \item{NW34}{Radius of >=34kt winds in northwest quadrant}
 \item{Hr\{n\}FcstDate}{Forecast valid date}
 \item{Hr\{n\}Lat}{Forecast latitude in `n` hours}
 \item{Hr\{n\}Lon}{Forecast longitude in `n` hours}
 \item{Hr\{n\}Wind}{Forecast maximum wind in `n` hours}
 \item{Hr\{n\}Gust}{Forecast maximum gust in `n` hours}
 \item{Hr\{n\}NE64}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}SE64}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}SW64}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}NW64}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}NE50}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}SE50}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}SW50}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}NW50}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}NE34}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}SE34}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}SW34}{Forecast wind radius in `n` hours}
 \item{Hr\{n\}NW34}{Forecast wind radius in `n` hours}
 \item{SeasNE}{Radius of 12ft seas in northeast quadrant}
 \item{SeasSE}{Radius of 12ft seas in southeast quadrant}
 \item{SeasSW}{Radius of 12ft seas in southwest quadrant}
 \item{SeasNW}{Radius of 12ft seas in northwest quadrant}
}
}
\examples{
\dontrun{
# Return dataframe of forecast/advisories for Tropical Storm Alex (AL011998)
get_fstadv("http://www.nhc.noaa.gov/archive/1998/1998ALEXadv.html")
}
}
\seealso{
\code{\link{tidy_fstadv}}, \code{\link{tidy_wr}},
\code{\link{tidy_fcst}}, \code{\link{tidy_fcst_wr}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_wind_radius}
\alias{fstadv_wind_radius}
\title{fstadv_wind_radius}
\usage{
fstadv_wind_radius(content)
}
\arguments{
\item{contents}{text of product}
}
\value{
dataframe
}
\description{
Parse wind radius data from product, if exists. This is somewhat
tricky as the wind fields are 64KT, 50KT and 34KT and are listed in
descending order. So the first line will not always be 64KT, 50KT or even
34KT depending on strength of storm. What I do here is just extract the
entire blob and work through it. I'll continue to look for ways to improve
it.

Complimentary to fstadv_get_wind_radius
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_storm_surge_flood}
\alias{gis_storm_surge_flood}
\title{gis_storm_surge_flood}
\usage{
gis_storm_surge_flood(key, advisory = as.numeric(),
  products = c("inundation", "tidalmask"))
}
\arguments{
\item{key}{Key of storm (i.e., AL012008, EP092015)}

\item{advisory}{Advisory number.}

\item{products}{indundation or tidalmask}
}
\description{
Potential Storm Surge Flooding (Inundation)
}
\seealso{
\code{\link{gis_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_fstadv}
\alias{filter_fstadv}
\title{filter_fstadv}
\usage{
filter_fstadv(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
Filter out forecast/advisory links
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_prblty}
\alias{filter_prblty}
\title{filter_prblty}
\usage{
filter_prblty(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
Get strike probabilities.
}
\details{
Strike probability products were terminated at the end of the 2005
season, replaced by wind probabilities.
}
\seealso{
\code{\link{filter_wndprb}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_forecasts}
\alias{fstadv_forecasts}
\title{fstadv_fcst}
\usage{
fstadv_forecasts(content, key, adv, adv_date)
}
\arguments{
\item{content}{text content of FORECAST/ADVISORY}

\item{key}{Storm ID}

\item{adv}{Advisory Number}

\item{adv_date}{Date value of forecast/advisory product.}
}
\value{
boolean
}
\description{
Retrieve forecast data from FORECAST/ADVISORY products. Loads
  into respective dataframes (df_forecasts, df_forecast_winds)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{knots_to_mph}
\alias{knots_to_mph}
\title{knots_to_mph}
\usage{
knots_to_mph(x)
}
\arguments{
\item{x}{wind speed in knots}
}
\value{
x in miles per hour
}
\description{
convert knots (kt) to miles per hour (mph)
}
\examples{
knots_to_mph(65)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wndprb.R
\name{wndprb}
\alias{wndprb}
\title{wndprb}
\usage{
wndprb(contents)
}
\arguments{
\item{contents}{Link to a storm's specific wind probability product.}
}
\description{
Parse wind probability products
}
\details{
Given a direct link to a wind probability product, parse and return
dataframe of values.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_posest}
\alias{filter_posest}
\title{filter_posest}
\usage{
filter_posest(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
Get position estimates
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_update}
\alias{filter_update}
\title{filter_update}
\usage{
filter_update(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
Get updates
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrapers.R
\name{scrape_date}
\alias{scrape_date}
\title{scrape_date}
\usage{
scrape_date(header)
}
\arguments{
\item{header}{Header text of product.}
}
\description{
Scrape date/time of product issuance from header.
}
\seealso{
\code{\link{scrape_header}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{extract_year_archive_link}
\alias{extract_year_archive_link}
\title{extract_year_archive_link}
\usage{
extract_year_archive_link(link)
}
\arguments{
\item{link}{URL of archive page}
}
\value{
year 4-digit numeric
}
\description{
Extracts the year from the archive link.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storms.R
\name{year_archives_link}
\alias{year_archives_link}
\title{year_archives_link}
\usage{
year_archives_link(year)
}
\arguments{
\item{year}{4-digit numeric}
}
\description{
Returns link to a year's archive page
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{get_update}
\alias{get_update}
\title{get_update}
\usage{
get_update(links)
}
\arguments{
\item{links}{URL to storm's archive page.}
}
\description{
Return dataframe of cyclone update data.
\describe{
  \item{Status}{Classification of storm, e.g., Tropical Storm, Hurricane,
  etc.}
  \item{Name}{Name of storm}
  \item{Date}{Date of advisory issuance}
  \item{Key}{Unique ID of cyclone}
  \item{Contents}{Text content of product}
}
}
\seealso{
\code{\link{get_storms}}, \code{\link{update}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.gis_wsp}
\alias{df.gis_wsp}
\title{GIS wind speed probabilities for Hurricane Sandy (AL182012)}
\format{An object of class \code{list} of length 3.}
\source{
\url{http://www.nhc.noaa.gov/gis/archive_wsp.php}
}
\usage{
df.gis_wsp
}
\description{
GIS wind speed probabilities for Hurricane Sandy (AL182012)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storm_data.R
\name{get_storm_data}
\alias{get_storm_data}
\title{get_storm_data}
\usage{
get_storm_data(links, products = c("discus", "fstadv", "posest",
  "public", "prblty", "update", "wndprb"))
}
\arguments{
\item{links}{to storm's archive page.}

\item{products}{Products to retrieve; discus, fstadv, posest, public,
prblty, update, and windprb.}
}
\value{
list of dataframes for each of the products.
}
\description{
Retrieve data from products.
}
\details{
\code{get_storm_data} is a wrapper function to make it more
  convenient to access the various storm products.

Types of products:
\describe{
  \item{discus}{Storm Discussions. This is technical information on the
  cyclone such as satellite presentation, forecast model evaluation, etc.}
  \item{fstadv}{Forecast/Advisory. These products contain the meat of an
  advisory package. Current storm information is available as well as
  structural design and forecast data.}
  \item{posest}{Position Estimate. Issued generally when a storm is
  threatening; provides a brief update on location and winds.}
  \item{public}{Public Advisory. Issued for public knowledge; more often for
  Atlantic than East Pacific storms. Contains general information.}
  \item{prblty}{Strike Probability. Discontinued after the 2005 hurricane
  season, strike probabilities list the chances of x-force winds in a
  particular city.}
  \item{update}{Cyclone Update. Generally issued when a significant change
  occurs in the cyclone.}
  \item{windprb}{Wind Probability. Replace strike probabilities beginning in
  the 2006 season. Nearly identical.}
}

Progress bars are displayed by default. These can be turned off by setting
the dplyr.show_progress to FALSE. Additionally, you can display messages for
each advisory being worked by setting the rrricanes.working_msg to TRUE.
}
\examples{
\dontrun{
## Get public advisories for first storm of 2016 Atlantic season.
get_storms(year = 2016, basin = "AL") \%>\%
  slice(1) \%>\%
  .$Link \%>\%
  get_storm_data(products = "public")
## Get public advisories and storm discussions for first storm of 2017 Atlantic season.
get_storms(year = 2017, basin = "AL") \%>\%
  slice(1) \%>\%
  .$Link \%>\%
  get_storm_data(products = c("discus", "public"))
}
}
\seealso{
\code{\link{get_ftp_storm_data}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{tidy_adv}
\alias{tidy_adv}
\alias{tidy_fstadv}
\title{tidy_adv}
\usage{
tidy_adv(df)

tidy_fstadv(df)
}
\arguments{
\item{df}{fstadv dataframe object}
}
\description{
Tidy current details of a fstadv dataframe object.

\code{tidy_fstadv} will be deprecated in 0.2.2
}
\details{
Returns current data only of a fstadv dataframe. Use Key, Adv and
Date to join with other tidy dataframes.
\describe{
 \item{Key}{Unique identifier of cyclone}
 \item{Adv}{Advisory number}
 \item{Date}{Date and time of advisory}
 \item{Status}{Classification of cyclone}
 \item{Name}{Name of cyclone}
 \item{Lat}{Latitude of cyclone center}
 \item{Lon}{Longitude of cyclone center}
 \item{Wind}{Maximum sustained one-minute winds in knots}
 \item{Gust}{Maximum sustained one-minute gusts in knots}
 \item{Pressure}{Minimum central pressure in millibars}
 \item{PosAcc}{Position accuracy of cyclone in nautical miles}
 \item{FwdDir}{Compass angle of forward motion}
 \item{FwdSpeed}{Forward speed in miles per hour}
 \item{Eye}{Size of eye in nautical miles}
 \item{SeasNE}{Radius of 12ft seas in northeast quadrant}
 \item{SeasSE}{Radius of 12ft seas in southeast quadrant}
 \item{SeasSW}{Radius of 12ft seas in southwest quadrant}
 \item{SeasNW}{Radius of 12ft seas in northwest quadrant}
}
}
\examples{
\dontrun{
get_fstadv("http://www.nhc.noaa.gov/archive/1998/1998ALEXadv.html") \%>\%
  tidy_adv()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{shp_to_df}
\alias{shp_to_df}
\title{shp_to_df}
\usage{
shp_to_df(obj)
}
\arguments{
\item{obj}{Spatial object to convert. See details.}
}
\description{
Convert shapefile object to dataframe
}
\details{
Takes a SpatialLinesDataFrame object or SpatialPolygonsDataFrame
object and converts into a dataframe that can be plotted in ggplot2.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storm_list.R
\name{get_storm_list}
\alias{get_storm_list}
\title{get_storm_list}
\usage{
get_storm_list()
}
\description{
Get storm list
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/products.R
\name{twoep}
\alias{twoep}
\title{twoep}
\usage{
twoep()
}
\description{
East Pacific Tropical Weather Outlook
}
\details{
This function parses the latest xml tropical weather outlook for
the east Pacific. The core data is located in the `channel$item` element
where `title`, `description` and `pubDate` reside. `link` is also
available to point to the NHC website.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storms.R
\name{get_storms}
\alias{get_storms}
\title{get_storms}
\format{A 4xN dataframe
\describe{
  \item{Year}{Numeric, four-digit year of the storm}
  \item{Name}{Character, name of storm mixed-case}
  \item{Basin}{AL (Atlantic) or EP (East Pacific)}
  \item{Link}{URL to storms' product pages}
}

To disable the progress bar set option dplyr.show_progress to FALSE.}
\source{
\url{http://www.nhc.noaa.gov/archive/2016/}
}
\usage{
get_storms(years = format(Sys.Date(), "\%Y"), basins = c("AL", "EP"))
}
\arguments{
\item{years}{numeric or vector, four digits (\%Y format)}

\item{basins}{One or both of c("AL", "EP")}
}
\value{
Dataframe of storms.
}
\description{
Returns storms and product link.
}
\details{
By default returns all storms for the current year. If no storms
have developed will return an empty dataframe.
}
\examples{
# Default. Get all storms, both basins, for last year.
\dontrun{
storms <- get_storms(year = 2016, basin = c("AL", "EP"))

# Get storms for two different years
storms.2010 <- get_storms(c(2010, 2015))

# Get storms for two consecutive years, Atlantic basin only
storms.al.2005 <- get_storms(2005:2007, basin = "AL")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrapers.R
\name{scrape_header}
\alias{scrape_header}
\title{scrape_header}
\usage{
scrape_header(contents, ptn_product_title, advisory_number = TRUE)
}
\arguments{
\item{contents}{Text product}

\item{ptn_product_title}{Pattern of product title to match}

\item{advisory_number}{Default is true; set to false if product does not
have an advisory number.}
}
\description{
Extract status, name, and advisory from products header.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.gis_wind_radii}
\alias{df.gis_wind_radii}
\title{GIS windfield and forecast wind radii for Hurricane Sandy (AL182012)}
\format{An object of class \code{list} of length 2.}
\source{
\url{http://www.nhc.noaa.gov/gis/archive_forecast_info_results.php?id=al18&year=2012&name=Hurricane\%20SANDY}
}
\usage{
df.gis_wind_radii
}
\description{
GIS windfield and forecast wind radii for Hurricane Sandy (AL182012)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_download}
\alias{gis_download}
\title{gis_download}
\usage{
gis_download(url, ...)
}
\arguments{
\item{url}{link to GIS dataset to download.}

\item{...}{additional parameters for rgdal::readOGR}
}
\description{
Get GIS data for storm.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_windfield}
\alias{gis_windfield}
\title{gis_windfield}
\usage{
gis_windfield(key, advisory = as.character())
}
\arguments{
\item{key}{Key of storm (i.e., AL012008, EP092015)}

\item{advisory}{Advisory number. If NULL, all advisories are returned.
Intermediate advisories are acceptable.}
}
\description{
Advisory Wind Field and Forecast Wind Radii
}
\details{
Tropical Cyclone Advisory Wind Field
 http://www.nhc.noaa.gov/gis/archive_forecast_info_results.php?id=al14&year=2016
 http://www.nhc.noaa.gov/gis/forecast/archive/
Example file name: al012017_fcst_001.zip
[basin]{2}[year_num]{2}[year]{4}_fcst_[advisory]{3}.zip
Many storms do not appear to have this data; especially earlier.

Not all advisories will be available for storms. For example,
\href{http://www.nhc.noaa.gov/gis/archive_forecast_info_results.php?id=al14&year=2016}{Hurricane Matthew (AL142016)}
is missing several advisories.
}
\seealso{
\code{\link{gis_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{month_str_to_num}
\alias{month_str_to_num}
\title{month_str_to_num}
\usage{
month_str_to_num(m)
}
\arguments{
\item{m}{Month abbreviated (SEP, OCT, etc.)}
}
\value{
numeric 1-12
}
\description{
Convert three-character month abbreviation to integer
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wndprb.R
\name{ep_prblty_stations}
\alias{ep_prblty_stations}
\title{ep_prblty_stations}
\usage{
ep_prblty_stations()
}
\description{
Retrieve list of probability stations based in the eastern
Pacific from the NHC. To be used in tandem with `wndprb` products.
}
\details{
Originally it was believed this data source would be removed by the
National Hurricane Center but it appears to have been updated. Additional
columns have been added, one up front and three in the back. These columns
all contain the same values each and I am unable to find documentation
describing the values.

Regardless, the data is kept, just in case.
}
\section{Warnings}{


Calling \code{ep_prblty_stations} will generate a warning:

> "Expected 7 pieces. Missing pieces filled with `NA` in 1 rows [41]."

Station SALINA CRUZ actually has six columns.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discus.R
\name{discus}
\alias{discus}
\title{discus}
\usage{
discus(contents)
}
\arguments{
\item{contents}{Link to a storm's specific discussion product.}
}
\value{
Dataframe
}
\description{
Parse storm Discussion products
}
\details{
Given a direct link to a discussion product, parse and return
dataframe of values.
}
\seealso{
\code{\link{get_discus}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\docType{package}
\name{rrricanes}
\alias{rrricanes}
\alias{rrricanes-package}
\title{rrricanes}
\description{
rrricanes is a web-scraping library for R designed to deliver
hurricane data (past and current) into well-organized datasets. With these
datasets you can explore past hurricane tracks, forecasts and structure
elements.

This documentation and additional help articles
\href{https://ropensci.github.io/rrricanes/}{can be found online}.

Text products (Forecast/Advisory, Public Advisory, Discussions and
Probabilities) are only available from 1998 to current. An effort will be
made to add prior data as available.
}
\section{Getting Storms}{

List all storms that have developed by year and basin. Year must be in a
four-digit format (\%Y) and no earlier than 1998. Basin can be one or both
of Atlantic ("AL") or East Pacific ("EP").
\describe{
  \item{\code{\link{get_storms}}}{List all storms by year, basin}
}
}

\section{Getting Storm Data}{


\code{\link{get_storm_data}} can be used to select multiple products,
multiple storms and from multiple basins.

Additional text products are:
\describe{
  \item{\code{\link{get_discus}}}{Storm Discussions}
  \item{\code{\link{get_fstadv}}}{Forecast/Advisory. These products contain a
  bulk of the information for tropical cyclones including current position,
  structure, forecast position and forecast structure.}
  \item{\code{\link{get_posest}}}{Position Estimates. Rare and used generally
  for threatening cyclones. This product was discontinued after the 2013
  season and is now issued as \code{\link{get_update}}.}
  \item{\code{\link{get_prblty}}}{Strike Probabilities. Show the probability
  of the center of a cyclone passing within 65nm of a location for a given
  forecast period. This product was discontinued after 2005, replaced with
  \code{\link{get_wndprb}}.}
  \item{\code{\link{get_public}}}{Public Advisory. General non-structured
  information exists in these products.}
  \item{\code{\link{get_update}}}{Updates. Generally issued when a cyclone
  undergoes a sudden change that requires immediate notice.}
  \item{\code{\link{get_wndprb}}}{Wind Speed Probability. Lists the
  probability of a location experiencing a minimum of 35kt, 50kt or 64kt
  winds for an alotted forecast period or accumulated probability. This
  product replaced \code{\link{get_prblty}} after the 2005 season.}
}

The products above may take some time to load if the NHC website is slow (as
is often the case, unfortunately). For all storm advisories issued outside
of the current month, use the \code{rrricanesdata} package.

To install \code{rrricanesdata}, run

\code{
install.packages("rrricanesdata",
         repos = "https://timtrice.github.io/drat/",
         type = "source")
}

See \code{vignette("installing_rrricanesdata", package = "rrricanes")} for
more information.
}

\section{GIS Data}{


For enhanced plotting of storm data, several GIS datasets are available. The
core GIS functions return URLs to help you refine the data you wish to view.
(Some products will not exist for all storms/advisories). These products are:

\describe{
  \item{\code{\link{gis_advisory}}}{Past track, current position, forecast
    and wind radii}
  \item{\code{\link{gis_breakpoints}}}{Breakpoints for watches and warnings}
  \item{\code{\link{gis_latest}}}{All available GIS products for active
    cyclones}
  \item{\code{\link{gis_outlook}}}{Tropical Weather Outlook}
  \item{\code{\link{gis_prob_storm_surge}}}{Probabilistic Storm Surge}
  \item{\code{\link{gis_windfield}}}{Wind Radii}
  \item{\code{\link{gis_wsp}}}{Wind Speed Probabilities}
}

\code{\link{gis_download}} will download the datasets from the above
functions.

Some GIS datasets will need to be converted to dataframes to plot geoms. Use
\code{\link{shp_to_df}} to convert SpatialLinesDataFrames and
SpatialPolygonsDataFrames. SpatialPointsDataFrames can be converted using
\code{tibble::as_data_frame} targeting the @data object.
}

\section{Package Options}{


\code{dplyr.show_progress} displays the dplyr progress bar when scraping raw
product datasets. In \code{\link{get_storms}}, it is based on the number of
years being requested. In the product functions (i.e.,
\code{\link{get_fstadv}}) it is based on the number of advisories. It can be
misleading when calling \code{\link{get_storm_data}} because it shows the
progress of working through a storm's product advisories but will reset on
new products/storms.



\code{dplyr.show_progress} displays the dplyr progress bar when scraping raw
product datasets. In \code{\link{get_storms}}, it is based on the number of
years being requested. In the product functions (i.e.,
\code{\link{get_fstadv}}) it is based on the number of advisories. It can be
misleading when calling \code{\link{get_storm_data}} because it shows the
progress of working through a storm's product advisories but will reset on
new products/storms.

\code{rrricanes.working_msg} is set to FALSE by default. When TRUE, it will
list the current storm, advisory and date being worked.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.al_18_2012_fstadv}
\alias{df.al_18_2012_fstadv}
\title{Forecast/Advisory for Hurricane Sandy (AL182012)}
\format{An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 31 rows and 117 columns.}
\source{
\url{http://www.nhc.noaa.gov/archive/2012/SANDY.shtml?}
}
\usage{
df.al_18_2012_fstadv
}
\description{
Forecast/Advisory for Hurricane Sandy (AL182012)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{mb_to_in}
\alias{mb_to_in}
\title{mb_to_in}
\usage{
mb_to_in(x)
}
\arguments{
\item{x}{barometric pressure in mb}
}
\value{
x in inches
}
\description{
convert millibars (mb) to inches of mercury (in)
}
\examples{
mb_to_in(999)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/products.R
\name{twoal}
\alias{twoal}
\title{twoal}
\usage{
twoal()
}
\description{
Atlantic Tropical Weather Outlook
}
\details{
This function parses the latest xml tropical weather outlook for
the Atlantic ocean. The core data is located in the `channel$item` element
where `title`, `description` and `pubDate` reside. `link` is also
available to point to the NHC website.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storms.R
\name{extract_storms}
\alias{extract_storms}
\title{extract_storms}
\usage{
extract_storms(basin, contents)
}
\arguments{
\item{basin}{AL or EP}

\item{contents}{contents of basin archive page}
}
\value{
4xN Dataframe
}
\description{
Extract storms for the given basin
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{status_abbr_to_str}
\alias{status_abbr_to_str}
\title{status_abbr_to_str}
\usage{
status_abbr_to_str(x)
}
\arguments{
\item{x}{character vector of status abbreviations}
}
\value{
character vector of strings
}
\description{
Convert Status abbreviation to string
}
\details{
Status abbreviations
\describe{
  \item{DB}{Disturbance (of any intensity)}
  \item{EX}{Extratropical cyclone (of any intensity)}
  \item{HU}{Tropical cyclone of hurricane intensity (> 64 knots)}
  \item{LO}{A low that is neither a tropical cyclone, a subtropical
        cyclone, nor an extratropical cyclone (of any intensity)}
  \item{SD}{Subtropical cyclone of subtropical depression intensity
        (< 34 knots)}
  \item{SS}{Subtropical cyclone of subtropical storm intensity
        (> 34 knots)}
  \item{TD}{Tropical cyclone of tropical depression intensity (< 34 knots)}
  \item{TS}{Tropical cyclone of tropical storm intensity (34-63 knots)}
  \item{WV}{Tropical Wave (of any intensity)}
}
}
\examples{
# Extratropical Cyclone
status_abbr_to_str("EX")

# Hurricane
status_abbr_to_str("HU")
}
\seealso{
\url{http://www.aoml.noaa.gov/hrd/hurdat/newhurdat-format.pdf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prblty.R
\name{prblty}
\alias{prblty}
\title{prblty}
\usage{
prblty(contents)
}
\arguments{
\item{contents}{Link to a storm's specific strike probability advisory product.}
}
\value{
Dataframe
}
\description{
Parse strike probability products
}
\details{
Given a direct link to a strike probability advisory product, parse
and return dataframe of values.
}
\seealso{
\code{\link{get_prblty}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wndprb.R
\name{parse_stations}
\alias{parse_stations}
\title{parse_stations}
\usage{
parse_stations(x)
}
\arguments{
\item{x}{URL of station list}
}
\description{
Parse probability station listings for each basin.
}
\details{
At the moment, documentation on the format is unavailable. The
format changed during the 2017/2018 offseason and now includes a
numeric first column and numeric fifth, sixth and seventh column. All
values are identical per column.

Additionally, as of publication, PATRICK AFB in the Atlantic data source
actually contains eight columns; this is noted in
\code{\link{al_prblty_stations}}. SALINA CRUZ in
\code{\link{ep_prblty_stations}} is short one column.

I see no issues with the extra or missing data as I am unsure the value of
the data to begin with. A warning will be given so the user is aware,
but the important pieces (Location, Lat, Lon) all seem good.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public.R
\name{get_public}
\alias{get_public}
\title{get_public}
\usage{
get_public(links)
}
\arguments{
\item{links}{URL to storm's archive page.}
}
\description{
Return dataframe of public advisory data.
\describe{
  \item{Status}{Classification of storm, e.g., Tropical Storm, Hurricane,
  etc.}
  \item{Name}{Name of storm}
  \item{Adv}{Advisory Number}
  \item{Date}{Date of advisory issuance}
  \item{Key}{Unique ID of the cyclone}
  \item{Contents}{Text content of product}
}
}
\seealso{
\code{\link{get_storms}}, \code{\link{public}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wndprb.R
\name{al_prblty_stations}
\alias{al_prblty_stations}
\title{al_prblty_stations}
\usage{
al_prblty_stations()
}
\description{
Retrieve list of probability stations based in the Atlantic
basin from the NHC. To be used in tandem with `wndprb` products.
}
\details{
Originally it was believed this data source would be removed by the
National Hurricane Center but it appears to have been updated. Additional
columns have been added, one up front and three in the back. These columns
all contain the same values each and I am unable to find documentation
describing the values.

Regardless, the data is kept, just in case.
}
\section{Warnings}{


Calling \code{al_prblty_stations} will generate a warning:

> "Expected 7 pieces. Additional pieces discarded in 1 rows [90]."

Station PATRICK AFB actually has eight columns. The data is kept for
consistency; you decide if you want it or not.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.al_18_2012_wndprb}
\alias{df.al_18_2012_wndprb}
\title{Wind speed probabilities for Hurricane Sandy (AL182012)}
\format{An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 2227 rows and 18 columns.}
\source{
\url{http://www.nhc.noaa.gov/archive/2012/SANDY.shtml?}
}
\usage{
df.al_18_2012_wndprb
}
\description{
Wind speed probabilities for Hurricane Sandy (AL182012)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{update}
\alias{update}
\title{update}
\usage{
update(contents)
}
\arguments{
\item{contents}{Link to a storm's specific update advisory product.}
}
\value{
Dataframe
}
\description{
Parse cyclone update products
}
\details{
Given a direct link to a cyclone update product, parse and return
dataframe of values.
}
\seealso{
\code{\link{get_update}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wndprb.R
\name{get_wndprb}
\alias{get_wndprb}
\title{get_wndprb}
\source{
\url{http://www.nhc.noaa.gov/about/pdf/About_Windspeed_Probabilities.pdf}
}
\usage{
get_wndprb(links)
}
\arguments{
\item{links}{URL to storm's archive page.}
}
\description{
Return dataframe of wind speed probability data.
}
\details{
Wind Speed Probability product replaced Strike Probabilities product
  after the 2005 hurricane season. These products may not be issued for
  every advisory/cyclone.

\describe{
  \item{Status}{Classification of storm, e.g., Tropical Storm, Hurricane,
  etc.}
  \item{Name}{Name of storm}
  \item{Adv}{Advisory Number}
  \item{Date}{Date of advisory issuance}
  \item{Wind}{Minimum wind speed for which probabilities reference}
  \item{Wind12}{Probability of sustained `Wind` within 12 hours}
  \item{Wind24}{Probability of sustained `Wind` within 24 hours}
  \item{Wind24Cum}{Cumulative probability through 24 hours}
  \item{Wind36}{Probability of sustained `Wind` within 36 hours}
  \item{Wind36Cum}{Cumulative probability through 36 hours}
  \item{Wind48}{Probability of sustained `Wind` within 48 hours}
  \item{Wind48Cum}{Cumulative probability through 48 hours}
  \item{Wind72}{Probability of sustained `Wind` within 72 hours}
  \item{Wind72Cum}{Cumulative probability through 72 hours}
  \item{Wind96}{Probability of sustained `Wind` within 96 hours}
  \item{Wind96Cum}{Cumulative probability through 96 hours}
  \item{Wind120}{Probability of sustained `Wind` within 120 hours}
  \item{Wind120Cum}{Cumulative probability through 120 hours}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_prev_pos}
\alias{fstadv_prev_pos}
\title{fstadv_prev_pos}
\usage{
fstadv_prev_pos(contents, adv_date)
}
\description{
Get storm's previous position
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_pressure}
\alias{fstadv_pressure}
\title{fstadv_pressure}
\usage{
fstadv_pressure(contents)
}
\arguments{
\item{contents}{text contents of FORECAST/ADVISORY product}
}
\value{
numeric
}
\description{
Return current minimum central pressure of storm in
  millibars (mb)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gis.R
\name{gis_advisory}
\alias{gis_advisory}
\title{gis_advisory}
\usage{
gis_advisory(key, advisory = as.character())
}
\arguments{
\item{key}{Key of storm (i.e., AL012008, EP092015)}

\item{advisory}{Advisory number. If NULL, all advisories are returned.
Intermediate advisories are acceptable.}
}
\description{
Advisory Forecast Track, Cone of Uncertainty, and
  Watches/Warnings
}
\seealso{
\code{\link{gis_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstadv.R
\name{fstadv_winds_gusts}
\alias{fstadv_winds_gusts}
\title{fstadv_winds_gusts}
\usage{
fstadv_winds_gusts(contents)
}
\arguments{
\item{contents}{text contents of FORECAST/ADVISORY product}
}
\value{
numeric
}
\description{
Get winds or gusts in knots (KT)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_discus}
\alias{filter_discus}
\title{filter_discus}
\usage{
filter_discus(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
Filter out storm discussion links.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{df.gis_adv}
\alias{df.gis_adv}
\title{GIS advisory dataset for Hurricane Sandy Adv 18}
\format{An object of class \code{list} of length 4.}
\source{
\url{http://www.nhc.noaa.gov/gis/archive_forecast_results.php?id=al18&year=2012&name=Hurricane\%20SANDY}
}
\usage{
df.gis_adv
}
\description{
GIS advisory dataset for Hurricane Sandy Adv 18
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter_wndprb}
\alias{filter_wndprb}
\title{filter_wndprb}
\usage{
filter_wndprb(links)
}
\arguments{
\item{links}{Vector of URLs retrieved from storm's archive page.}
}
\description{
Get wind probabilities
}
\details{
Wind probability products replaced the strike probability products
at the beginning of the 2006 season.
}
\seealso{
\code{\link{filter_prblty}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_storm_data.R
\name{get_product}
\alias{get_product}
\title{get_product}
\usage{
get_product(links, product)
}
\description{
This funtion acts as a hub for the individual product extraction
  functions. Given the product and links, it will begin the scraping
  process and return a dataset for that product.
}
\keyword{internal}
