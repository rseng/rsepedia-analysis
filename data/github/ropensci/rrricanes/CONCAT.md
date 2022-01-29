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
