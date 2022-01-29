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
# _GSODR_ <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->
[![tic](https://github.com/ropensci/GSODR/workflows/tic/badge.svg?branch=main)](https://github.com/ropensci/GSODR/actions)
[![codecov](https://codecov.io/gh/ropensci/GSODR/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/GSODR)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.439850.svg)](https://doi.org/10.5281/zenodo.439850) 
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/GSODR)](https://cran.r-project.org/package=GSODR)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) 
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.00177/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00177) 
[![](https://badges.ropensci.org/79_status.svg)](https://github.com/ropensci/software-review/issues/79)
<!-- badges: end -->

--------------------

## A Global Surface Summary of the Day (GSOD) Weather Data Client for R

## Introduction

The GSOD or [Global Surface Summary of the Day (GSOD)](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00516) data provided by the US National Centers for Environmental Information (NCEI) are a valuable source of weather data with global coverage.
_**GSODR**_ aims to make it easy to find, transfer and format the data you need for use in analysis and provides five main functions for facilitating this:

- `get_GSOD()` - this function queries and transfers files from the NCEI's web server, reformats them and returns a data frame.

- `reformat_GSOD()` - this function takes individual station files from the local disk and re-formats them returning a data frame.

- `nearest_stations()` - this function returns a vector of station IDs that fall within the given radius (kilometres) of a point given as latitude and longitude in order from nearest to farthest.

- `update_station_list()` - this function downloads the latest station list from the NCEI's server updates the package's internal database of stations and their metadata.

- `get_inventory()` - this function downloads the latest station inventory information from the NCEI's server and returns the header information about the latest version as a message in the console and a tidy data frame of the stations' inventory for each month that data are reported.

When reformatting data either with `get_GSOD()` or `reformat_GSOD()`, all units are converted to International System of Units (SI), _e.g._, inches to millimetres and Fahrenheit to Celsius.
File output is returned as a `data.table` object, summarising each year by station, which also includes vapour pressure and relative humidity elements calculated from existing data in GSOD.
Additional data are calculated by this R package using the original data and included in the final data.
These include vapour pressure (ea and es) and relative humidity calculated using the improved August-Roche-Magnus approximation (Alduchov and Eskridge 1996).

For more information see the description of the data provided by NCEI, <https://www7.ncdc.noaa.gov/CDO/GSOD_DESC.txt>.

## How to Install

### Stable Version

A stable version of _GSODR_ is available from [CRAN](https://cran.r-project.org/package=GSODR).

```r
install.packages("GSODR")
```

### Development Version

A development version is available from from GitHub.
If you wish to install the development version that may have new features or bug fixes before the CRAN version does (but also may not work properly), please install the [remotes](https://github.com/r-lib/remotes) package, available from CRAN.
We strive to keep the main branch on GitHub functional and working properly.

```r
if (!require("remotes")) {
  install.packages("remotes", repos = "http://cran.rstudio.com/")
  library("remotes")
}

install_github("ropensci/GSODR")
```

## Other Sources of Weather Data in R

There are several other sources of weather data and ways of retrieving them through R.
Several are also [rOpenSci](https://ropensci.org) projects.

The [_**GSODTools**_](https://github.com/environmentalinformatics-marburg/GSODTools) by [Florian Detsch](https://github.com/fdetsch) is an R package that offers similar functionality as _**GSODR**_, but also has the ability to graph the data and working with data for time series analysis.

The [_**gsod**_](https://github.com/databrew/gsod) package from [DataBrew](https://www.databrew.cc/posts/gsod.html) aims to streamline the way that researchers and data scientists interact with and utilise weather data and relies on _**GSODR**_, but provides data in the package rather than downloading so it is faster (though available data may be out of date).

[_**rnoaa**_](https://CRAN.R-project.org/package=rnoaa), from [rOpenSci](https://docs.ropensci.org/rnoaa/) offers tools for interacting with and downloading weather data from the United States National Oceanic and Atmospheric Administration but lacks support for GSOD data.

[_**stationaRy**_](https://cran.r-project.org/package=stationaRy), from Richard Iannone offers hourly meteorological data from stations located all over the world.
There is a wealth of data available, with historic weather data accessible from nearly 30,000 stations.

[_**riem**_](https://CRAN.R-project.org/package=riem) from [rOpenSci](https://docs.ropensci.org/riem/) allows to get weather data from Automated Surface Observing System (ASOS) stations (airports) in the whole world thanks to the Iowa Environment Mesonet website.

[_**weathercan**_](https://CRAN.R-project.org/package=weathercan) from [rOpenSci](https://github.com/ropensci/weathercan) makes it easier to search for and download multiple months/years of historical weather data from Environment and Climate Change Canada (ECCC) website.

[_**clifro**_](https://CRAN.R-project.org/package=clifro) from [rOpenSci](https://docs.ropensci.org/clifro/) is a web portal to the New Zealand National Climate Database and provides public access (via subscription) to around 6,500 various climate stations (see <https://cliflo.niwa.co.nz/> for more information).
Collating and manipulating data from CliFlo (hence clifro) and importing into R for further analysis, exploration and visualisation is now straightforward and coherent.
The user is required to have an Internet connection, and a current CliFlo subscription (free) if data from stations, other than the public Reefton electronic weather station, is sought.

## Notes

### NOAA policy

Users of these data should take into account the following (from the
[NCEI website](https://www7.ncdc.noaa.gov/CDO/cdoselect.cmd?datasetabbv=GSOD&countryabbv=&georegionabbv=)): 

> The following data and products may have conditions placed on their international commercial use. They can be used within the U.S. or for non-commercial international activities without restriction. The non-U.S. data cannot be redistributed for commercial purposes. Re-distribution of these data by others must provide this same notification. A log of IP addresses accessing these data and products will be maintained and may be made available to data providers.  
For details, please consult: [WMO Resolution 40. NOAA Policy](https://community.wmo.int/resolution-40)

## Meta

- Please [report any issues or bugs](https://github.com/ropensci/GSODR/issues).

- License: MIT

- To cite _**GSODR**_, please use: Adam H Sparks, Tomislav Hengl and Andrew Nelson (2017). GSODR: Global Summary Daily Weather Data in R. _The Journal of Open Source Software_, **2(10)**. DOI: 10.21105/joss.00177.

- Please note that the _**GSODR**_ project is released with a [Contributor Code of Conduct](https://github.com/ropensci/GSODR/blob/main/CONDUCT.md) By participating in the _**GSODR**_ project you agree to abide by its terms.

## References

Alduchov, O.A. and Eskridge, R.E., 1996. Improved Magnus form approximation of saturation vapor pressure. Journal of Applied Meteorology and Climatology, 35(4), pp. 601-609 DOI: 10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2.
# GSODR 3.1.4

# GSODR 3.1.2.9000

## Minor changes

* Skip **ANY** and **ALL** tests on CRAN.
This fixes the "problems" with _GSODR_ failing on a Solaris instance when the server failed to respond.

* Update internal `isd-history` database.

* Use `\CRANpkg{}` in place of `\pkg{}` in documentation.

# GSODR 3.1.2

## Bug fixes

* Fix (more) bugs related to `NA` value replacements.

## Minor changes and improvements

* Simplify `NA` value replacement in "PRCP" column.

* The PRCP column values are rounded to two decimal places as in original GSOD data, not one.

* The TEMP_ATTRIBUTES, DEWP_ATTRIBUTES, SLP_ATTRIBUTES, STP_ATTRIBUTES, VISIB_ATTRIBUTES and WDSP_ATTRIBUTES columns are formatted as an integer not character.

* Better tests for the generated weather data `data.table` output checking values and formats.

* Tests are updated for updated data availability in the GSOD data due to continuous improvements to the data set.

* Standardise handling of author/contributor comments.
None have a full stop now in the comment.

* Use `on.exit()` to reset the working directory to the original user-space value after changing the working directory to untar files located in `tempdir()`.

# GSODR 3.1.1

## Bug fixes

* Fixes bug reported in [#84](https://github.com/ropensci/GSODR/issues/84) in the FRSHTT columns where the values were all reported as `NA` even if there were observed values.

* Fixes bug where NA values reported as 99.99, 999.9 or 9999.9 were not replaced with `NA`.

* Fix bug where FRSHTT (Fog, Rain/Drizzle, Snow/Ice, Hail, Tornado, Thunder) column values split into the respective columns only returned `NA`, not the proper values as expected.
Reported in [#84](https://github.com/ropensci/GSODR/issues/84).

## Minor changes

* Examples are no longer wrapped in `\donttest{}` but use `@examplesIf interactive()` instead.

# GSODR 3.1.0

## New features

* Include columns for COUNTRY_NAME (country name in English), ISO2C and ISO3C in the final output from `get_GSOD()` and `reformat_GSOD()`.

* Diffs in the isd_history are now recorded in the `/data-raw/fetch_isd-history.md` file and shipped with GSODR as `isd_history.rda`, which can be viewed by using `load(system.file("extdata", "isd_diff.rda", package = "GSODR"))`.

* Update and improve documentation to reflect country name and ISO code columns.

## Minor improvements

* Update NCEI data usage statement.

# GSODR 3.0.0

## Breaking changes

* Remove parallel processing functionality.
A bug that I was unable to properly debug with `future.apply::future_lapply()` caused the `get_GSOD()` and `reformat_GSOD()` functions to run without completing or responding was fixed by simply using R's base `lapply()` function.
If parallel processing is needed, users should implement their own solutions to download and process in parallel.

## Bug fixes

* Fix bug that caused the package to run without responding.

* Fix test that failed on CRAN's Solaris server for some reason.

* Removes a working DOI link from the reference for the equation used because win-builder checks say it doesn't work (even though it does and there's nothing wrong with the link any time I check).

# GSODR 2.1.2

## Bug fixes

* Fix bug where `nearest_stations()` did not always return the nearest station as the first value in the vector

## Minor changes

* Update internal isd-history database, adding 11 stations

* Fix any links that redirect found in DESCRIPTION, documentation or other materials in the package

# GSODR 2.1.1

## Bug fixes 

* Fix bug where station metadata files could not be updated

## Minor changes

* Update internal station list to latest

* Correct an error in documentation for `update_station_list()`

* Remove spatial vignettes to slim down Suggests and make CI maintenance easier

# GSODR v2.1.0

## Major changes

* Implement new calculations for EA, ES and RH using improved August-Roche-Magnus approximation (Alduchov & Eskridge 1996).
HT Rich Iannone for his use in [stationaRy](https://cran.r-project.org/package=stationaRy).
This will result in different EA, ES and RH calculations from the prior versions of GSODR.
However, this new implementation should be more accurate as discussed in (Alduchov & Eskridge 1996).

> Alduchov, O.A. and Eskridge, R.E., 1996. Improved Magnus form approximation of saturation vapor pressure. Journal of Applied Meteorology and Climatology, 35(4), pp.601-609.

## Minor changes

* Update internal station list to latest

* Enhanced documentation

# GSODR v2.0.1
 
## Bug fixes

* Corrects internal bug that provided a warning message when GSOD files were parsed

* Fixes bug where not all files downloaded were cleaned up on the function exit when fetching station inventories

* Fixes bug where station inventories from `get_inventory()` lacked the location metadata, _i.e._ country and other geographic information

## Minor changes

* Update vignette to use latest functions from tidyr, _i.e._ `tidyr::pivot_longer()`

* Update internal station list to latest

* Tidy up documentation, mainly fix functions' title capitalisation

* `get_GSOD()` checks the number of stations being requested.
If the number is \>10, the entire annual file will be downloaded and requested stations will then
be selected and returned.
This saves time by reducing the number of requests made to the server.
Users should not see any difference other than quicker responses for a large number of requested stations.

## Major changes

* Requires R >= 3.5.0 due to the storage of .Rds files using the latest version

# GSODR 2.0.0

## Bug fixes

* `get_GSOD()` now uses https rather than FTP server, correcting bug where the
data could not be downloaded any longer

## Major changes

* Corrected elevation values are no longer available from GSODR

* Objects are returned as `data.table` objects

## Minor changes

* `get_inventory()` now uses https rather than FTP server

* `update_station_list()` now uses https rather than FTP server

* Implement better error handling when attempting to fetch station inventories

* Reduced package dependencies

* Improved vignettes that are pre-compiled for faster package installation and
updated content with linting and error corrections

* Users may now specify country by FIPS code when using `get_GSOD()`

* Improved test coverage

* Update internal database of station locations

# GSODR 1.3.2

## Bug fixes

* Fixes a bug where extra data could be appended to data frame. See
<https://github.com/ropensci/GSODR/issues/49>.
This also means that when you are retrieving large amounts of data, _e.g._ global data for 20+ years, you won't fill up your hard disk space due to the raw data before processing.

## Minor changes

* Update internal database of station locations

# GSODR 1.3.1

## Bug fixes

* Fix examples that did not run properly

## Minor changes

* Update internal database of station locations

# GSODR 1.3.0

## New Functionality

* Use `future_apply` in processing files after downloading.
This allows for end users to use a parallel process of their choice.

# GSODR 1.2.3

## Bug fixes

* Refactor internal functionality to be more clear and efficient in execution

    * `country-list` is not loaded unless user has specified a country in
      `get_GSOD()`

    * An instance where the FIPS code was determined twice was removed

* Replace `\dontrun{}` with `\donttest{}` in documentation examples

* Ensure that DESCRIPTION file follows CRAN guidelines

## Minor changes

* Format help files, fixing errors and formatting for attractiveness

* Update internal database of station locations

* Store internal database of station locations fields `BEGIN` and `END` as
  integer, not double

* Clarify code of conduct statement in README that it only applies to this,
  GSODR, project

* Prompt user for input with warning about reproducibility if using the
  `update_station_list()` function

* Adds metadata header to the `tibble` returned by `get_inventory()`

* Remove start-up message to conform with rOpenSci guidelines

* Remove extra code, clean up code-chunks and use `hrbrthemes::theme_ipsum()`
  for
  [data-raw/fetch_isd-history.md](https://github.com/ropensci/GSODR/blob/main/data-raw/fetch_isd-history.md)

# GSODR 1.2.2

## Bug fixes

  * Fix bug in creating `isd-history.rda` file where duplicate stations existed
  in the file distributed with `GSODR` but with different corrected elevation
  values

  * Repatch bug reported and fixed previously in version 1.2.0 where Windows
  users could not successfully download files. This somehow snuck back in.

## Minor changes

  * Refactor vignettes for clarity

# GSODR 1.2.1

## Bug fixes

  * Introduce a message if a station ID is requested but files are not found
    on the server.
    This is in response to an inquiry from John Paul Bigouette where a station is reported as having data in the inventory but the files do not exist on the server.

  * Fix bug that removed a few hundred stations from the internal `GSODR` database of stations in the `data-raw` files.

## Minor changes

  * Clean documentation, shortening long lines, fixing formatting,
    incomplete sentences and broken links

  * Clarify the reasons for errors that a user may encounter

  * Update internal databases of station metadata

  * Clean up this file

# GSODR 1.2.0

## Major changes

  * Remove ability to export files from `get_GSOD()` to slim down the package
    dependencies and this functions parameters.
    Examples of how to convert to a spatial object (both _sp_ and _sf_ are shown) and export ESRI Shapefiles and
    GeoPackage files are now included in the vignette.

  * As a result of the previous point, the _sp_ and _rgdal_ packages are no longer Imports but are now in Suggests along with _sf_ for examples in the GSOD vignette.

## Bug fixes

  * Fix a nasty bug where GSOD files downloaded using Windows would not untar properly.
  This caused the `get_GSOD()` function to fail.
  Thanks to Ross Darnell, CSIRO, for reporting this.

  * Correct options in "GSODR use case: Specified years/stations vignette" on line 201 where `file` was incorrectly used in place of `path`.
  Thanks to Ross Darnell, CSIRO, for reporting this.

  * Correct documentation for `reformat_GSOD()`

## Minor changes

  * Update internal databases of station metadata

  * Vignettes contain pre-built figures for faster package installation when building vignettes

# GSODR 1.1.2

## Bug fixes

  * Fix start-up message formatting

  * Correct ORCID comment in author field of DESCRIPTION

  * Update internal databases for country list and isd_history

## Minor changes

  * Add X-schema tags to DESCRIPTION

# GSODR 1.1.1

## Bug fixes

  * `MAX_FLAG` and `MIN_FLAG` columns now report `NA` when there is no flag

## Minor changes

  * Comment for Bob and Hugh in DESCRIPTION now only ORCID url

  * dplyr version set to >= 0.7.0 not 0.7 as before

  * Start-up message statement is more clear in relation to WMO resolution 40, that GSODR does not redistribute any weather data itself

  * Remove unnecessary function, .onLoad(), from zzz.R

  * Function titles in documentation now in title case

  * Correct grammar in documentation

# GSODR 1.1.0

## Bug fixes

  * Fixes bug reported in
    [issue 36](https://github.com/ropensci/GSODR/issues/36)

## Major changes

  * The _data.table_ and _fields_ packages are no longer imported.
  All internal functions now use _dplyr_ or base R functionality, reducing the dependencies of _GSODR_

  * Any data frames returned by _GSODR_ functions are returned as a `tibble()` object

  * The `YEARMODA` column is now returned as `Date` without time, rather than `Character`

  * Add new function, `get_inventory()`, which downloads the NCEI's station inventory document and returns a `tibble()` object of the data

  * Use larger images and provide a table of contents in vignettes

  * Updated and enhanced introductory vignette

  * Update internal stations list

# GSODR 1.0.7

## Bug fixes

  * Fix documentation in vignette where first example would not run due to changes in package data formats

  * Fix bug in GSODR vignette where examples would not run due to libraries not being loaded

  * Fix bug where prior server queries would be pre/appended to subsequent queries

  * Fix bug where invalid stations would return an empty data frame, should stop and return message about checking the `station` value supplied to `get_GSOD()` and check if data are available for the years requested

## Minor changes

  * Update Appendix 2 of GSODR vignette, map of station locations, to be more clear and follow same format as that of `bomrang` package

  * Update example output in GSODR vignette where applicable

## Major changes

  * Update internal stations list


# GSODR 1.0.6

## Bug fixes

  * Fix bug where WSPD (mean wind-speed) conversion was miscalculated

# GSODR 1.0.5

## Major changes

  * Add welcome message on start-up regarding data use and sharing

  * Update internal stations list

## Minor changes

  * Tidy up informative messages that the package returns while running

## Bug fixes

  * Fix bug where "Error in read_connection_(con):" when writing to CSV occurs

  * Fix typo in line 160 of `get_GSOD()` where "Rda" should be "rda" to properly load internal package files

# GSODR 1.0.4

## Major changes

  * Data distributed with GSODR are now internal to the package and not externally exposed to the user

  * Vignettes have been updated and improved with an improved order of information presented and some have been combined for easier use

## Minor changes

  * Clean code using linting

# GSODR 1.0.3

## Major changes

  * Data for station locations and unique identifiers is now provided with the package on installation. Previously this was fetched each time from the ftp server.

  * The station metadata can now be updated if necessary by using `update_station_list()`, this change overwrites the internal data that were originally distributed with the package.
  This operation will fetch the latest list of stations and corresponding information from the NCEI ftp
    server.
    Any changes will be overwritten when the R package is updated, however, the package update should have the same or newer data included, so this should not be an issue.

  * Replace _plyr_ functions with _purrr_, _plyr_ is no longer actively developed

  * _plyr_ is no longer an import

## Minor changes

  * Fix bugs in the vignettes related to formatting and spelling

## Deprecated and defunct

  * `get_station_list()` is no longer supported. Instead use the new

  * `update_station_list()` to update the package's internal station database.

# GSODR 1.0.2.1

## Minor changes

  * Correct references to _GSODRdata_ package where incorrectly referred to as _GSODdata_

# GSODR 1.0.2

## Minor changes

  * Improved documentation (i.e., spelling corrections and more descriptive)

  * More descriptive vignette for "GSODR use case: Specified years/stations vignette"

  * Round MAX/MIN temp to one decimal place, not two

  * Update SRTM elevation data

  * Update country list data

  * Fix missing images in README.html on CRAN

# GSODR 1.0.1

## Minor changes

  * Update documentation for `get_GSOD()` when using `station` parameter

  * Edit paper.md for submission to JOSS

  * Remove extra packages listed as dependencies that are no longer necessary

  * Correct Working_with_spatial_and_climate_data.Rmd where it was missing the first portion of documentation and thus examples did not work

# GSODR 1.0.0

## Major changes

  * The `get_GSOD()` function returns a `data.frame` object in the current R session with the option to save data to local disk

  * Multiple stations can be specified for download rather than just downloading a single station or all stations

  * A new function, `nearest_stations()` is now included to find stations within a user specified radius (in kilometres) of a point given as latitude and longitude in decimal degrees

  * A general use vignette is now included

  * New vignette with a detailed use-case

  * Output files now include fields for State (US only) and Call (International Civil Aviation Organization (ICAO) Airport Code)

  * Use FIPS codes in place of ISO3c for file name and in output files because some stations do not have an ISO country code

  * Spatial file output is now in GeoPackage format (GPKG). This results in a single file output unlike shapefile and allows for long field names

  * Users can specify file name of output

  * R >= 3.2.0 now required

  * Field names in output files use "\_" in place of "."

  * Long field names now used in file outputs

  * Country is specified using FIPS codes in file name and output file contents due to stations occurring in some locales that lack ISO 3166 3 letter country codes

  * The `get_GSOD()` function will retrieve the latest station data from NCDC and automatically merge it with the CGIAR-CSI SRTM elevation values provided by this package. Previously, the package provided it's own list of station information, which was difficult to keep up-to-date

  * A new `reformat_GSOD()` function reformats station files in "WMO-WBAN-YYYY.op.gz" format that have been downloaded from the United States National Climatic Data Center's (NCDC) FTP server.

  * A new function, `get_station_list()` allows for fetching latest station list from the FTP server and querying by the user for a specified station or location.

  * New data layers are provided through a separate package, [`GSODRdata`](https://github.com/adamhsparks/GSODRdata), which provide
    climate data formatted for use with GSODR.

     * CHELSA (climatic surfaces at 1 km resolution),

     * MODCF  * Remotely sensed high-resolution global cloud dynamics for predicting ecosystem and biodiversity distributions,

     * ESACCI  * ESA's CCI-LC snow cover probability and

     * CRU CL2.0 (climatic surfaces at 10 minute resolution).

  * Improved file handling for individual station downloads

  * Missing values are handled as `NA` not -9999

  * Change from GPL >= 3 to MIT licence to bring into line with ropensci packages

  * Now included in ropensci, [ropensci/GSODR](https://github.com/ropensci/GSODR)

## Minor changes

  * `get_GSOD()` function optimised for speed as best possible after FTPing files from NCDC server

  * All files are downloaded from server and then locally processed, previously these were sequentially downloaded by year and then processed

  * A progress bar is now shown when processing files locally after downloading

  * Reduced package dependencies

  * The `get_GSOD()` function now checks stations to see if the years being queried are provided and returns a message alerting user if the station and years requested are not available

  * When stations are specified for retrieval using the `station = ""` parameter, the `get_GSOD()` function now checks to see if the file exists on the server, if it does not, a message is returned and all other stations that have files are processed and returned in output

  * Documentation has been improved throughout package

  * Better testing of internal functions

## Bug Fixes

  * Fixed: Remove redundant code in `get_GSOD()` function

  * Fixed: The stations data frame distributed with the package now include stations that are located above 60 latitude and below -60 latitude

## Deprecated and defunct

  * Missing values are reported as NA for use in R, not -9999 as previously

  * The `path` parameter is now instead called `dsn` to be more inline with other tools like `readOGR()` and `writeOGR()`

  * Shapefile file out is no longer supported. Use GeoPackage (GPKG) instead

  * The option to remove stations with too many missing days is now optional, it now defaults to including all stations, the user must specify how many missing stations to check for an exclude

  * The `max_missing` parameter is now user set, defaults to no check, return all stations regardless of missing days

# GSODR 0.1.9

## Bug Fixes

  * Fix bug in precipitation calculation. Documentation states that PRCP is in mm to hundredths.
  Issues with conversion and missing values meant that this was not the case.
  Thanks to Gwenael Giboire for reporting and help with fixing this

## Minor changes

  * Users can now select to merge output for station queries across multiple years. Previously one year = one file per station.
  Now are set by user, `merge_station_years = TRUE` parameter, only one output file is generated

  * Country list is now included in the package to reduce run time necessary when querying for a specific country.
  However, this means any time that the country-list.txt file is updated, this package needs to be updated as well

  * Updated `stations` list with latest version from NCDC published 12-07-2016

  * Country list is now included in the package to reduce run time necessary when querying for a specific country.
  However, this means any time that the country-list.txt file is updated, this package needs to be updated as well

  * Country level, agroclimatology and global data query conversions and calculations are processed in parallel now to reduce runtime

  * Improved documentation with spelling fixes, clarification and updates

  * Enable `ByteCompile` option upon installation for small increase in speed

  * Use `write.csv.raw` from
    `[iotools]("https://cran.r-project.org/web/packages/iotools/index.html")` to greatly improve runtime by decreasing time used to write CSV files to disk

  * Use `writeOGR()` from `rgdal`, in place of `raster's` `shapefile` to improve runtime by decreasing time used to write shapefiles to disk

  * Country level, agroclimatology and global data query conversions and calculations are processed in parallel now to reduce runtime

  * Improved documentation with spelling fixes, clarification and updates

  * Enable `ByteCompile` option upon installation for small increase in speed

  * Use `write.csv.raw` from `[iotools]("https://cran.r-project.org/web/packages/iotools/index.html")`
     to greatly improve runtime by decreasing time used to write CSV files to disk

# GSODR 0.1.8

## Bug Fixes

  * Fix bug with connection timing out for single station queries commit: [a126641e00dc7acc21844ff0436e5702f8b6e04a](https://github.com/ropensci/GSODR/commit/a126641e00dc7acc21844ff0436e5702f8b6e04a)

  * Somehow the previously working function that checked country names broke with the `toupper()` function.
  A new [function from juba](https://stackoverflow.com/questions/16516593/convert-from-lowercase-to-uppercase-all-values-in-all-character-variables-in-dat) fixes this issue and users can now select country again

  * User entered values for a single station are now checked against actual station values for validity

  * stations.rda is compressed

  * stations.rda now includes a field for "corrected" elevation using hole-filled SRTM data from Jarvis et al. 2008, see [https://github.com/ropensci/GSODR/blob/main/data-raw/fetch_isd-history.md](https://github.com/ropensci/GSODR/blob/devel/data-raw/fetch_isd-history.md) for a description

  * Set NA or missing values in CSV or shapefile to -9999 from -9999.99 to
    align with other data sources such as WorldClim

## Minor changes

  * Documentation is more complete and easier to use

# GSODR 0.1.7

## Bug Fixes

  * Fix issues with MIN/MAX where MIN referred to MAX [(Issue 5)](https://github.com/ropensci/GSODR/issues/5)

  * Fix bug where the `tf` item was incorrectly set as `tf <  * "~/tmp/GSOD-2010.tar`, not `tf <  * tempfile`, in `get_GSOD()` [(Issue 6)](https://github.com/ropensci/GSODR/issues/6)

  * CITATION file is updated and corrected

## Minor changes

  * User now has the ability to generate a shapefile as well as CSV file output [(Issue 3)](https://github.com/ropensci/GSODR/issues/3)

  * Documentation is more complete and easier to use

# GSODR 0.1.6

## Bug Fixes

  * Fix issue when reading .op files into R where temperature was incorrectly read causing negative values where T >= 100F, this issue caused RH values of >100% and incorrect TEMP values [(Issue 1)](https://github.com/ropensci/GSODR/issues/1)

  * Spelling corrections

## Major changes

  * Include MIN/MAX flag column

  * Station data is now included in package rather than downloading from NCDC every time get_GSOD() is run, this data has some corrections where stations with missing LAT/LON values or elevation are omitted, this is **not** the original complete station list provided by NCDC.

# GSODR 0.1.5

## Bug Fixes

  * Fixed bug where YDAY not correctly calculated and reported in CSV file

  * CSV files for station only queries now are names with the Station Identifier.
  Previously named same as global data

  * Likewise, CSV files for agroclimatology now are names with the Station Identifier.
  Previously named same as global data.

## Minor Changes

  * Set values where MIN > MAX to NA

  * Set more MIN/MAX/DEWP values to NA.
  GSOD README indicates that 999 indicates missing values in these columns, this does not appear to always be true.
  There are instances where 99 is the value recorded for missing data.
  While 99F is possible, the vast majority of these recorded values are missing data, thus the function now converts them to NA


# GSODR 0.1.4

## Bug Fixes

  * Fixed bug related to MIN/MAX columns when agroclimatology or all stations are selected where flags were not removed properly from numeric values.

# GSODR 0.1.3

## Bug fixes

  * Bug fix in MIN/MAX with flags. Some columns have differing widths, which caused a flag to be left attached to some values

  * Correct URL in README.md for CRAN to point to CRAN not GitHub

## Minor Changes

  * Set NA to -9999.99

# GSODR 0.1.2

## Bug Fixes

  * Bug fix in importing isd-history.csv file. Previous issues caused all lat/lon/elev values to be >0.

  * Bug fix where WDSP was mistyped as WDPS causing the creation of a new column, rather than the conversion of the existing

  * Bug fix if Agroclimatology selected. Previously this resulted in no records.

  * Set the default encoding to UTF8.

  * Bug fix for country selection. Some countries did not return proper ISO code.

  * Bug fix where WDSP was mistyped as WDPS causing the creation of a new column,   rather than the conversion of the existing

  * Use write.csv, not readr::write_csv due to issue converting double to string: <https://github.com/tidyverse/readr/issues/387>

# GSODR 0.1.1

## Major changes

* Now available on CRAN

* Add single quotes around possibly misspelled words and spell out comma-separated values and geographic information system rather than just using "CSV" or "GIS" in DESCRIPTION.

* Add full name of GSOD (Global Surface Summary of the Day) and URL for GSOD, <https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00516> to DESCRIPTION as requested by CRAN.

* Require user to specify directory for resulting .csv file output so that any files written to disk are interactive and with user's permission

# GSODR 0.1

* Initial submission to CRAN
# GSODR 3.1.2

## Comments for the CRAN team

I realise that the last release was last week, however, I found a few more bugs in the code that could cause some surprises with `NA` values for end users, so I've elected to release a new patch release.

## Test environments
* GitHub Actions (ubuntu-latest): release, devel
* GitHub Actions (windows): release
* Github Actions (macOS): release
* Local macOS M1: release
* win-builder: devel

## R CMD check results
0 ERRORs | 0 WARNINGs | 0 NOTES.

## Reverse dependencies
We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
 
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

Please note that the GSODR project is released with a
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
Fetch and Clean ‘isd_history.csv’ File
================
Adam H. Sparks
2021-10-06

<STYLE type='text/css' scoped>
PRE.fansi SPAN {padding-top: .25em; padding-bottom: .25em};
</STYLE>

# Introduction

The isd_history.csv file details GSOD station metadata. These data
include the start and stop years used by *GSODR* to pre-check requests
before querying the server for download and the country code used by
*GSODR* when sub-setting for requests by country. The following checks
are performed on the raw data file before inclusion in *GSODR*,

-   Check for valid lon and lat values;

    -   isd_history where latitude or longitude are `NA` or both 0 are
        removed leaving only properly georeferenced stations,

    -   isd_history where latitude is \< -90˚ or > 90˚ are removed,

    -   isd_history where longitude is \< -180˚ or > 180˚ are removed.

-   A new field, STNID, a concatenation of the USAF and WBAN fields, is
    added.

# Data Processing

## Set up workspace

``` r
library("sessioninfo")
library("skimr")
library("countrycode")
library("data.table")
```

## Download and clean data

``` r
# download data
new_isd_history <- fread("https://www1.ncdc.noaa.gov/pub/data/noaa/isd-history.csv")
```

## Add/drop columns and save to disk

``` r
# add STNID column
new_isd_history[, STNID := paste(USAF, WBAN, sep = "-")]
setcolorder(new_isd_history, "STNID")
setnames(new_isd_history, "STATION NAME", "NAME")

# drop stations not in GSOD data
new_isd_history[, STNID_len := nchar(STNID)]
new_isd_history <- subset(new_isd_history, STNID_len == 12)

# remove stations where LAT or LON is NA
new_isd_history <- na.omit(new_isd_history, cols = c("LAT", "LON"))

# remove extra columns
new_isd_history[, c("USAF", "WBAN", "ICAO", "ELEV(M)", "STNID_len") := NULL]
```

## Add country names based on FIPS

``` r
new_isd_history <-
  new_isd_history[setDT(countrycode::codelist), on = c("CTRY" = "fips")]

new_isd_history <- new_isd_history[, c(
  "STNID",
  "NAME",
  "LAT",
  "LON",
  "CTRY",
  "STATE",
  "BEGIN",
  "END",
  "country.name.en",
  "iso2c",
  "iso3c"
)]

# clean data
new_isd_history[new_isd_history == -999] <- NA
new_isd_history[new_isd_history == -999.9] <- NA
new_isd_history <-
  new_isd_history[!is.na(new_isd_history$LAT) &
                    !is.na(new_isd_history$LON),]
new_isd_history <-
  new_isd_history[new_isd_history$LAT != 0 &
                    new_isd_history$LON != 0,]
new_isd_history <-
  new_isd_history[new_isd_history$LAT > -90 &
                    new_isd_history$LAT < 90,]
new_isd_history <-
  new_isd_history[new_isd_history$LON > -180 &
                    new_isd_history$LON < 180,]

# set colnames to upper case
names(new_isd_history) <- toupper(names(new_isd_history))
setnames(new_isd_history,
         old = "COUNTRY.NAME.EN",
         new = "COUNTRY_NAME")

# set country names to be upper case for easier internal verifications
new_isd_history[, COUNTRY_NAME := toupper(COUNTRY_NAME)]

# set key for joins when processing CSV files
setkeyv(new_isd_history, "STNID")[]
```

    ##               STNID                         NAME    LAT      LON CTRY STATE
    ##     1: 008268-99999                    WXPOD8278 32.950   65.567   AF      
    ##     2: 010010-99999          JAN MAYEN(NOR-NAVY) 70.933   -8.667   NO      
    ##     3: 010014-99999                   SORSTOKKEN 59.792    5.341   NO      
    ##     4: 010015-99999                   BRINGELAND 61.383    5.867   NO      
    ##     5: 010016-99999                  RORVIK/RYUM 64.850   11.233   NO      
    ##    ---                                                                     
    ## 26562: A00024-53848 CHOCTAW NAVAL OUTLYING FIELD 30.507  -86.960   US    FL
    ## 26563: A00026-94297              COUPEVILLE/NOLF 48.217 -122.633   US    WA
    ## 26564: A00029-63820      EVERETT-STEWART AIRPORT 36.380  -88.985   US    TN
    ## 26565: A00030-93795        CONNELLSVILLE AIRPORT 39.959  -79.657   US    PA
    ## 26566: A00032-25715                 ATKA AIRPORT 52.220 -174.206   US    AK
    ##           BEGIN      END  COUNTRY_NAME ISO2C ISO3C
    ##     1: 20100519 20120323   AFGHANISTAN    AF   AFG
    ##     2: 19310101 20211003        NORWAY    NO   NOR
    ##     3: 19861120 20211003        NORWAY    NO   NOR
    ##     4: 19870117 20081231        NORWAY    NO   NOR
    ##     5: 19870116 19910806        NORWAY    NO   NOR
    ##    ---                                            
    ## 26562: 20070601 20211003 UNITED STATES    US   USA
    ## 26563: 20060324 20150514 UNITED STATES    US   USA
    ## 26564: 20130627 20211004 UNITED STATES    US   USA
    ## 26565: 20210309 20211004 UNITED STATES    US   USA
    ## 26566: 20060101 20211003 UNITED STATES    US   USA

## Show changes from last release

``` r
# ensure we aren't using a locally installed dev version
install.packages("GSODR", repos = "https://cloud.r-project.org/")
```

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/hc/tft3s5bn48gb81cs99mycyf00000gn/T//Rtmp8NcPVB/downloaded_packages

``` r
load(system.file("extdata", "isd_history.rda", package = "GSODR"))

# select only the cols of interest
x <- names(isd_history)
new_isd_history <- new_isd_history[, ..x] 

(isd_diff <- diffobj::diffPrint(new_isd_history, isd_history))
```

<PRE class="fansi fansi-output"><CODE>## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>new_isd_history</span>                                                            
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>isd_history</span>                                                                
## <span style='color: #00BBBB;'>@@ 1,3 / 1,3 @@                                                              </span>
## <span style='color: #BBBB00;'>&lt;</span>               STNID                         NAME    LAT      LON CTRY <span style='color: #BBBB00;'>STATE</span>
## <span style='color: #0000BB;'>&gt;</span>               STNID                            NAME    LAT      LON CTRY   
##       1: 008268-99999                    WXPOD8278 32.950   65.567   AF      
##       2: 010010-99999          JAN MAYEN(NOR-NAVY) 70.933   -8.667   NO      
## <span style='color: #00BBBB;'>@@ 6,19 / 6,19 @@                                                            </span>
##       5: 010016-99999                  RORVIK/RYUM 64.850   11.233   NO      
##      ---                                                                     
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26527:</span> <span style='color: #0000BB;'>A00023-63890</span> <span style='color: #0000BB;'>WHITEHOUSE</span> <span style='color: #0000BB;'>NAVAL</span> <span style='color: #0000BB;'>OUTLYING</span> <span style='color: #0000BB;'>FIELD</span> <span style='color: #0000BB;'>30.350</span>  <span style='color: #0000BB;'>-81.883</span>   <span style='color: #0000BB;'>US</span>   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26562:</span> A00024-53848 CHOCTAW NAVAL OUTLYING FIELD 30.507  -86.960   US    <span style='color: #BBBB00;'>FL</span>
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26528:</span> A00024-53848    CHOCTAW NAVAL OUTLYING FIELD 30.507  -86.960   US   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26563:</span> A00026-94297              COUPEVILLE/NOLF 48.217 -122.633   US    <span style='color: #BBBB00;'>WA</span>
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26529:</span> A00026-94297                 COUPEVILLE/NOLF 48.217 -122.633   US   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26564:</span> A00029-63820      EVERETT-STEWART AIRPORT 36.380  -88.985   <span style='color: #BBBB00;'>US</span>    <span style='color: #BBBB00;'>TN</span>
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26565:</span> <span style='color: #BBBB00;'>A00030-93795</span>        <span style='color: #BBBB00;'>CONNELLSVILLE</span> <span style='color: #BBBB00;'>AIRPORT</span> <span style='color: #BBBB00;'>39.959</span>  <span style='color: #BBBB00;'>-79.657</span>   US    <span style='color: #BBBB00;'>PA</span>
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26530:</span> A00029-63820         EVERETT-STEWART AIRPORT 36.380  -88.985   US   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26566:</span> A00032-25715                 ATKA AIRPORT 52.220 -174.206   US    <span style='color: #BBBB00;'>AK</span>
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26531:</span> A00032-25715                    ATKA AIRPORT 52.220 -174.206   US   
## <span style='color: #BBBB00;'>&lt;</span>           BEGIN      END  COUNTRY_NAME ISO2C ISO3C                         
## <span style='color: #0000BB;'>&gt;</span>        <span style='color: #0000BB;'>STATE</span>    BEGIN      END  COUNTRY_NAME ISO2C ISO3C                   
##       1: 20100519 20120323   AFGHANISTAN    AF   AFG                         
## <span style='color: #BBBB00;'>&lt;</span>     2: 19310101 <span style='color: #BBBB00;'>20211003</span>        NORWAY    NO   NOR                         
## <span style='color: #0000BB;'>&gt;</span>     2:       19310101 <span style='color: #0000BB;'>20210116</span>        NORWAY    NO   NOR                   
## <span style='color: #BBBB00;'>&lt;</span>     3: 19861120 <span style='color: #BBBB00;'>20211003</span>        NORWAY    NO   NOR                         
## <span style='color: #0000BB;'>&gt;</span>     3:       19861120 <span style='color: #0000BB;'>20210116</span>        NORWAY    NO   NOR                   
##       4: 19870117 20081231        NORWAY    NO   NOR                         
##       5: 19870116 19910806        NORWAY    NO   NOR                         
##      ---                                                                     
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26562:</span> 20070601 <span style='color: #BBBB00;'>20211003</span> UNITED STATES    US   USA                         
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26527:</span>    <span style='color: #0000BB;'>FL</span> 20070601 <span style='color: #0000BB;'>20210116</span> UNITED STATES    US   USA                   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26563:</span> <span style='color: #BBBB00;'>20060324</span> <span style='color: #BBBB00;'>20150514</span> UNITED STATES    US   USA                         
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26528:</span>    <span style='color: #0000BB;'>FL</span> <span style='color: #0000BB;'>20070601</span> <span style='color: #0000BB;'>20210116</span> UNITED STATES    US   USA                   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26564:</span> <span style='color: #BBBB00;'>20130627</span> <span style='color: #BBBB00;'>20211004</span> UNITED STATES    US   USA                         
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26529:</span>    <span style='color: #0000BB;'>WA</span> <span style='color: #0000BB;'>20060324</span> <span style='color: #0000BB;'>20150514</span> UNITED STATES    US   USA                   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26565:</span> <span style='color: #BBBB00;'>20210309</span> <span style='color: #BBBB00;'>20211004</span> UNITED STATES    US   USA                         
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26530:</span>    <span style='color: #0000BB;'>TN</span> <span style='color: #0000BB;'>20130627</span> <span style='color: #0000BB;'>20210117</span> UNITED STATES    US   USA                   
## <span style='color: #BBBB00;'>&lt;</span> <span style='color: #BBBB00;'>26566:</span> 20060101 <span style='color: #BBBB00;'>20211003</span> UNITED STATES    US   USA                         
## <span style='color: #0000BB;'>&gt;</span> <span style='color: #0000BB;'>26531:</span>    <span style='color: #0000BB;'>AK</span> 20060101 <span style='color: #0000BB;'>20210117</span> UNITED STATES    US   USA
</CODE></PRE>

## View and save the data

``` r
str(isd_history)
```

    ## Classes 'data.table' and 'data.frame':   26531 obs. of  11 variables:
    ##  $ STNID       : chr  "008268-99999" "010010-99999" "010014-99999" "010015-99999" ...
    ##  $ NAME        : chr  "WXPOD8278" "JAN MAYEN(NOR-NAVY)" "SORSTOKKEN" "BRINGELAND" ...
    ##  $ LAT         : num  33 70.9 59.8 61.4 64.8 ...
    ##  $ LON         : num  65.57 -8.67 5.34 5.87 11.23 ...
    ##  $ CTRY        : chr  "AF" "NO" "NO" "NO" ...
    ##  $ STATE       : chr  "" "" "" "" ...
    ##  $ BEGIN       : int  20100519 19310101 19861120 19870117 19870116 19880320 19861109 19850601 19730101 19310103 ...
    ##  $ END         : int  20120323 20210116 20210116 20081231 19910806 20050228 20210114 20210116 20140523 20041030 ...
    ##  $ COUNTRY_NAME: chr  "AFGHANISTAN" "NORWAY" "NORWAY" "NORWAY" ...
    ##  $ ISO2C       : chr  "AF" "NO" "NO" "NO" ...
    ##  $ ISO3C       : chr  "AFG" "NOR" "NOR" "NOR" ...
    ##  - attr(*, ".internal.selfref")=<externalptr> 
    ##  - attr(*, "sorted")= chr "STNID"

``` r
# write rda file to disk for use with GSODR package
save(isd_history,
     file = "../inst/extdata/isd_history.rda",
     compress = "bzip2")

save(isd_diff,
     file = "../inst/extdata/isd_diff.rda",
     compress = "bzip2")
```

# Notes

## NOAA policy

Users of these data should take into account the following (from the
[NCEI
website](https://www7.ncdc.noaa.gov/CDO/cdoselect.cmd?datasetabbv=GSOD&countryabbv=&georegionabbv=)):

> The following data and products may have conditions placed on their
> international commercial use. They can be used within the U.S. or for
> non-commercial international activities without restriction. The
> non-U.S. data cannot be redistributed for commercial purposes.
> Re-distribution of these data by others must provide this same
> notification. A log of IP addresses accessing these data and products
> will be maintained and may be made available to data providers.  
> For details, please consult: [WMO Resolution 40. NOAA
> Policy](https://community.wmo.int/resolution-40)

## R System Information

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 4.1.1 (2021-08-10)
    ##  os       macOS Big Sur 11.6          
    ##  system   aarch64, darwin20           
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_AU.UTF-8                 
    ##  ctype    en_AU.UTF-8                 
    ##  tz       Australia/Perth             
    ##  date     2021-10-06                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date       lib source        
    ##  assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.1.0)
    ##  base64enc     0.1-3   2015-07-28 [1] CRAN (R 4.1.0)
    ##  cli           3.0.1   2021-07-17 [1] CRAN (R 4.1.0)
    ##  countrycode * 1.3.0   2021-07-15 [1] CRAN (R 4.1.0)
    ##  crayon        1.4.1   2021-02-08 [1] CRAN (R 4.1.0)
    ##  curl          4.3.2   2021-06-23 [1] CRAN (R 4.1.0)
    ##  data.table  * 1.14.2  2021-09-27 [1] CRAN (R 4.1.1)
    ##  DBI           1.1.1   2021-01-15 [1] CRAN (R 4.1.0)
    ##  diffobj       0.3.5   2021-10-05 [1] CRAN (R 4.1.1)
    ##  digest        0.6.28  2021-09-23 [1] CRAN (R 4.1.1)
    ##  dplyr         1.0.7   2021-06-18 [1] CRAN (R 4.1.0)
    ##  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.1.0)
    ##  evaluate      0.14    2019-05-28 [1] CRAN (R 4.1.0)
    ##  fansi         0.5.0   2021-05-25 [1] CRAN (R 4.1.0)
    ##  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.1.0)
    ##  generics      0.1.0   2020-10-31 [1] CRAN (R 4.1.0)
    ##  glue          1.4.2   2020-08-27 [1] CRAN (R 4.1.0)
    ##  htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.1.1)
    ##  jsonlite      1.7.2   2020-12-09 [1] CRAN (R 4.1.0)
    ##  knitr         1.36    2021-09-29 [1] CRAN (R 4.1.1)
    ##  lifecycle     1.0.1   2021-09-24 [1] CRAN (R 4.1.1)
    ##  magrittr      2.0.1   2020-11-17 [1] CRAN (R 4.1.0)
    ##  pacman      * 0.5.1   2019-03-11 [1] CRAN (R 4.1.0)
    ##  pillar        1.6.3   2021-09-26 [1] CRAN (R 4.1.1)
    ##  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.1.0)
    ##  purrr         0.3.4   2020-04-17 [1] CRAN (R 4.1.0)
    ##  R6            2.5.1   2021-08-19 [1] CRAN (R 4.1.1)
    ##  repr          1.1.3   2021-01-21 [1] CRAN (R 4.1.0)
    ##  rlang         0.4.11  2021-04-30 [1] CRAN (R 4.1.0)
    ##  rmarkdown     2.11    2021-09-14 [1] CRAN (R 4.1.1)
    ##  rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.1.0)
    ##  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 4.1.0)
    ##  skimr       * 2.1.3   2021-03-07 [1] CRAN (R 4.1.0)
    ##  stringi       1.7.4   2021-08-25 [1] CRAN (R 4.1.1)
    ##  stringr       1.4.0   2019-02-10 [1] CRAN (R 4.1.0)
    ##  tibble        3.1.5   2021-09-30 [1] CRAN (R 4.1.1)
    ##  tidyselect    1.1.1   2021-04-30 [1] CRAN (R 4.1.0)
    ##  utf8          1.2.2   2021-07-24 [1] CRAN (R 4.1.0)
    ##  vctrs         0.3.8   2021-04-29 [1] CRAN (R 4.1.0)
    ##  withr         2.4.2   2021-04-18 [1] CRAN (R 4.1.0)
    ##  xfun          0.26    2021-09-14 [1] CRAN (R 4.1.1)
    ##  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.1.0)
    ## 
    ## [1] /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library

# GSODR data-raw contents

## Fetch isd-history

This document details how the GSOD station history data file,
["isd-history.csv"](https://www1.ncdc.noaa.gov/pub/data/noaa/isd-history.csv),
is fetched from the NCEI server and saved for inclusion in the _GSODR_ 
package in /data/isd_history.rda. These data are used for determining the years
that a station reported data for filtering user requests before sending them to
the server to reduce failed requests.

[fetch_isd-history.md](fetch_isd-history.md)
---
title: 'GSODR: Global Summary Daily Weather Data in R'
authors:
- affiliation: 1
  name: Adam H Sparks
  orcid: 0000-0002-0061-8359
- affiliation: 2
  name: Tomislav Hengl
  orcid: 0000-0002-9921-5129
- affiliation: 3
  name: Andrew Nelson
  orcid: 0000-0002-7249-3778
date: "27 January 2017"
output: pdf_document
bibliography: paper.bib
tags:
- Global Surface Summary of the Day
- GSOD
- meteorology
- climatology
- weather data
- R
affiliations:
- index: 1
  name: Centre for Crop Health, University of Southern Queensland, Toowoomba Queensland
    4350, Australia
- index: 2
  name: ISRIC - World Soil Information, P.O. Box 353, 6700 AJ Wageningen, The Netherlands
- index: 3
  name: Faculty of Geo-Information and Earth Observation (ITC), University of Twente,
    Enschede 7500 AE, The Netherlands
---

# Summary

The GSODR package [@GSODR] is an R package [@R-base] providing automated
downloading, parsing and cleaning of Global Surface Summary of the
Day (GSOD) [@NCDC] weather data for use in R or saving as local files in either
a Comma Separated Values (CSV) or GeoPackage (GPKG) [@geopackage] file. It
builds on or complements several other scripts and packages. We take advantage
of modern techniques in R to make more efficient use of available computing
resources to complete the process, e.g., data.table [@data.table], plyr [@plyr]
and readr [@readr], which allow the data cleaning, conversions and disk 
input/output processes to function quickly and efficiently. The rnoaa [@rnoaa]
package already offers an excellent suite of tools for interacting with and
downloading weather data from the United States National Oceanic and 
Atmospheric Administration, but lacks options for GSOD data retrieval. Several
other APIs and R packages exist to access weather data, but most are region or
continent specific, whereas GSOD is global. This package was developed to
provide:

  * two functions that simplify downloading GSOD data and formatting it to
  easily be used in research; and

  * a function to help identify stations within a given radius of a point of
interest.

Alternative elevation data based on a 200 meter buffer of 
elevation values derived from the CGIAR-CSI SRTM 90m Database [@Jarvis2008]
are included. These data are useful to help address possible inaccuracies and
in many cases, fill in for missing elevation values in the reported station
elevations.

When using this package, GSOD stations are checked for inaccurate longitude and
latitude values and any stations that have missing or have incorrect values are
omitted from the final data set. Users may set a threshold for station files
with too many missing observations for omission from the final output to help
ensure data quality. All units are converted from the United States Customary
System (USCS) to the International System of Units (SI), e.g., inches to
millimetres and Fahrenheit to Celsius. Wind speed is also converted from knots
to metres per second. Additional useful values, actual vapour pressure,
saturated water vapour pressure, and relative humidity are calculated and
included in the final output. Station metadata are merged with weather data for
the final data set.

# References
## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.1.0 (2021-05-18) |
|os       |macOS Big Sur 11.6           |
|system   |x86_64, darwin17.0           |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_AU.UTF-8                  |
|ctype    |en_AU.UTF-8                  |
|tz       |Australia/Perth              |
|date     |2021-10-06                   |

# Dependencies

|package     |old   |new        |Δ  |
|:-----------|:-----|:----------|:--|
|GSODR       |3.1.2 |3.1.2.9000 |*  |
|countrycode |1.3.0 |1.3.0      |   |
|data.table  |NA    |1.14.2     |*  |
|mime        |NA    |0.12       |*  |
|openssl     |NA    |1.4.5      |*  |
|R.utils     |NA    |2.11.0     |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*