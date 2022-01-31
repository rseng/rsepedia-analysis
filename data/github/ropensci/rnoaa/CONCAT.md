rnoaa
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rnoaa)](https://cranchecks.info/pkgs/rnoaa)
[![R-check](https://github.com/ropensci/rnoaa/workflows/R-check/badge.svg)](https://github.com/ropensci/rnoaa/actions)
[![codecov.io](https://codecov.io/github/ropensci/rnoaa/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rnoaa?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rnoaa?color=C9A115)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rnoaa)](https://cran.r-project.org/package=rnoaa)


`rnoaa` is an R interface to many NOAA data sources. We don't cover all of them, but we include many commonly used sources, and add we are always adding new sources. We focus on easy to use interfaces for getting NOAA data, and giving back data in easy to use formats downstream. We currently don't do much in the way of plots or analysis. To get started see: https://docs.ropensci.org/rnoaa/articles/rnoaa.html

## Data sources in rnoaa

* NOAA NCDC climate data:
    * We are using the NOAA API version 2
    * Docs for the NCDC API are at https://www.ncdc.noaa.gov/cdo-web/webservices/v2
    * GHCN Daily data is available at http://www.ncdc.noaa.gov/ghcn-daily-description via FTP and HTTP
* Severe weather data docs are at https://www.ncdc.noaa.gov/swdiws/
* Sea ice data (ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/shapefiles)
* NOAA buoy data (https://www.ndbc.noaa.gov/)
* ERDDAP data (https://upwell.pfeg.noaa.gov/erddap/index.html)
  * Now in package rerddap (https://github.com/ropensci/rerddap)
* Tornadoes! Data from the NOAA Storm Prediction Center (https://www.spc.noaa.gov/gis/svrgis/)
* HOMR - Historical Observing Metadata Repository (http://www.ncdc.noaa.gov/homr/api)
* GHCND FTP data (ftp://ftp.ncdc.noaa.gov/pub/data/noaa) - NOAA NCDC API has some/all (not sure really) of this data, but FTP allows to get more data more quickly
* Extended Reconstructed Sea Surface Temperature (ERSST) data (https://www.ncdc.noaa.gov/data-access/marineocean-data/extended-reconstructed-sea-surface-temperature-ersst-v4)
* Argo buoys (http://www.argo.ucsd.edu/) - a global array of more than 3,000 free-drifting profiling floats that measures thetemperature and salinity of the upper 2000 m of the ocean
* NOAA CO-OPS - tides and currents data (https://tidesandcurrents.noaa.gov/)
* NOAA Climate Prediction Center (CPC) (http://www.cpc.ncep.noaa.gov/)
* Africa Rainfall Climatology version 2 (ftp://ftp.cpc.ncep.noaa.gov/fews/fewsdata/africa/arc2/ARC2_readme.txt)
* Blended Sea Winds (https://www.ncdc.noaa.gov/data-access/marineocean-data/blended-global/blended-sea-winds)
* Local Climatological Data (https://www.ncdc.noaa.gov/cdo-web/datatools/lcd)
* Storm Events Database (https://www.ncdc.noaa.gov/stormevents/)

## Help/Getting Started

Documentation is at https://docs.ropensci.org/rnoaa/, and there are many vignettes in the package itself, available in your R session, or on CRAN (https://cran.r-project.org/package=rnoaa). The tutorials:

* **Getting started - start here**
* NOAA Buoy vignette
* NOAA National Climatic Data Center (NCDC) vignette (examples)
* NOAA NCDC attributes vignette
* NOAA NCDC workflow vignette
* Sea ice vignette
* Severe Weather Data Inventory (SWDI) vignette
* Historical Observing Metadata Repository (HOMR) vignette

## netcdf data

Some functions use netcdf files, including:

* `ersst`
* `buoy`
* `bsw`
* `argo`
 
You'll need the `ncdf4` package for those functions, and those only. `ncdf4` is in Suggests in this package, meaning you only need `ncdf4` if you are using any of the functions listed above. You'll get an informative error telling you to install `ncdf4` if you don't have it and you try to use the those functions. Installation of `ncdf4` should be straightforward on any system. See https://cran.r-project.org/package=ncdf4

## NOAA NCDC Datasets

There are many NOAA NCDC datasets. All data sources work, except `NEXRAD2` and `NEXRAD3`, for an unknown reason. This relates to `ncdc_*()` functions only.


|Dataset    |Description                 |Start Date |End Date   | Data Coverage|
|:----------|:---------------------------|:----------|:----------|-------------:|
|GHCND      |Daily Summaries             |1763-01-01 |2021-11-26 |          1.00|
|GSOM       |Global Summary of the Month |1763-01-01 |2021-11-01 |          1.00|
|GSOY       |Global Summary of the Year  |1763-01-01 |2021-01-01 |          1.00|
|NEXRAD2    |Weather Radar (Level II)    |1991-06-05 |2021-11-27 |          0.95|
|NEXRAD3    |Weather Radar (Level III)   |1994-05-20 |2021-11-26 |          0.95|
|NORMAL_ANN |Normals Annual/Seasonal     |2010-01-01 |2010-01-01 |          1.00|
|NORMAL_DLY |Normals Daily               |2010-01-01 |2010-12-31 |          1.00|
|NORMAL_HLY |Normals Hourly              |2010-01-01 |2010-12-31 |          1.00|
|NORMAL_MLY |Normals Monthly             |2010-01-01 |2010-12-01 |          1.00|
|PRECIP_15  |Precipitation 15 Minute     |1970-05-12 |2014-01-01 |          0.25|
|PRECIP_HLY |Precipitation Hourly        |1900-01-01 |2014-01-01 |          1.00|


```
#> table updated on 2021-11-29
```

**NOAA NCDC Attributes**

Each NOAA dataset has a different set of attributes that you can potentially get back in your search. See https://www.ncdc.noaa.gov/cdo-web/datasets for detailed info on each dataset. We provide some information on the attributes in this package; see the vignette for attributes (https://docs.ropensci.org/rnoaa/articles/ncdc_attributes.html) to find out more


## Contributors

* Scott Chamberlain (https://github.com/sckott)
* Daniel Hocking (https://github.com/djhocking)
* Brooke Anderson (https://github.com/geanders)
* Maëlle Salmon (https://github.com/maelle)
* Adam Erickson (https://github.com/adam-erickson)
* Nicholas Potter (https://github.com/potterzot)
* Joseph Stachelek (https://github.com/jsta)

## Meta

* Please report any issues or bugs: https://github.com/ropensci/rnoaa/issues
* License: MIT
* Get citation information for `rnoaa` in R doing `citation(package = 'rnoaa')`
* Please note that this package is released with a Contributor Code of Conduct (https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rnoaa 1.3.8
==============

### BUG FIXES

* changed location of temporary cache file writing in the ersst tests to match requirements of CRAN for Mac OS. Missed on v1.3.7 release. 
* removed rappsdir from Imports now that using tools to create cache directories
* removed a couple other internal tests from CRAN with skip_on_cran()


rnoaa 1.3.7
==============

### BUG FIXES

* changed location of temporary file writing to match requirements of CRAN for Mac OS
* removed checks on Ubuntu 16.04. Replaced with checks on latest Ubuntu version


rnoaa 1.3.4
===========

### MINOR IMPROVEMENTS

* update URL for tornadoes data to include data from 2019 (#386) thanks @ryanscharf

### BUG FIXES

* fix for all `ncdc*` functions - response header content type changed - we had a check for proper content type - that check is now more general so that any json content type will be okay (#390)

rnoaa 1.3.2
===========

### MINOR IMPROVEMENTS

* remove ropenaq from Suggests as its been archived on CRAN (#385)
* test fixes (#382)

rnoaa 1.3
=========

### NEW FEATURES

* `ghcnd()` now accepts more than 1 station identifier (#373) PR from @eliocamp
* `ersst()`: use new v5 version of their service - see `?ersst` docs for details (#381) thanks @vonStadarhraun for the tip

### MINOR IMPROVEMENTS

* update `buoy()` docs to state that a special value of `9999` passsed to the `year` parameter will give the most up to date data (aka current data) - and an example added using it (#377)
* update to current dplyr functions from deprecated ones (#375)
* update description of units TMAX and TMIN for dataset GHCND when using function `ncdc()` with `add_units = TRUE` (#378) (#379) PR from @amcdavid

### BUG FIXES

* changed `ghcnd()` handling of unknown/bad/invalid station identifiers: now returns an empty data.frame and gives back empty strings for the two attributes `source` and `file_modified`  (#374)

rnoaa 1.2.0
===========

### DEFUNCT

* The following IBTrACS storm functions are now defunct because they are too cumbersome to maintain: `storm_data()`, `storm_meta()`, `storm_shp()`, and `storm_shp_read()`. associated package datasets `storm_columns` and `storm_names` removed  (#306)

### MINOR IMPROVEMENTS

* vignettes are now all pre-built, and all URLs are now unlinked (#367)
* `lcd()` returns a tibble now instead of a tibble with S3 class `lcd` attached (#369)
* created manual file entry for `stormevents_cache`
* drop `sf` from Suggests, only used in an example

### BUG FIXES

* `arc2()` fix: when using `bbox` parameter, it would not have worked as intended, fixed now (#372)

rnoaa 1.1.0
===========

### MINOR IMPROVEMENTS

* fix coops test (#364)
* remove deprecated parameters in argo and ncdc* functions (#361)

rnoaa 1.0.0
===========

### NEW FEATURES

* the `argo` functions that were not working because of a down API are working again, see `?argo` (#358)
* most of the defunct functions have been removed from the package, but are still referenced in the `?rnoaa-defunct` manual file. `gefs` functions are still in the package as those functions may come back at some point. (#359)
* two things for `arc2()`: 1) now accepts more than 1 date; 2) gains new parameter `box` to accept a bounding box to spatially filter results (uses `dplyr::filter` on the data.frame of spatial data) (#351)

### Documentation

Regarding the documentation site at https://docs.ropensci.org/rnoaa

* the function reference page https://docs.ropensci.org/rnoaa/reference/index.html has been improved; grouping functions by topic area and data source - for easier browsing (#360)
* a getting started vignette has been added, see the "Get Started" tab (#357)

### MINOR IMPROVEMENTS

* `ghncd()` (and all functions that build on `ghcnd()`) can now be altered to use a specific base URL for requests. See the "Base URL" section of the `?ghcnd` docs (#353)
* `ghcnd_splitvars()` speedup, using `data.table` instead of dplyr for manipulation (#352) (#355)
* use `tibble::as_tibble` throughout package instead of `dplyr::tbl_df` (#354)

### BUG FIXES

* fix for `ghcnd_stations()`: internal method `get_inventory()` was not creating a directory first before trying to download a file into that directory (#349) (#350)


rnoaa 0.9.6
===========

### NEW FEATURES

* new function `rnoaa_options()` to toggle package level options; only option for now is `cache_messages`, a boolean to toggle whether the user gets messages about cached files or not. along with this change, messages about cached files and file sizes and locations are now consistently used across all functions that cache files on disk  (#331)
* new manual file `?rnoaa_caching` - with information on how to access and manage cached files for each of the rnoaa functions that caches files on disk (#346)
* `isd()` moved to using `hoardr` caching - see `?isd_cache` for details (#347)

### MINOR IMPROVEMENTS

* remove internal code in many exported functions looking for user input `path` parameter and telling them it's no longer used; been defunct for quite a while
* now able to cache http requests for tests that write to disk (#290) (#345)
* `ghcnd_stations()` now caching data - first time requests should now take just over 1 minute, with subsequent requests (assuming cached data isn't deleted) taking ~ 3 seconds  (#164)
* `autoplot` method `meteo_coverage()` fix to visually display gaps in data (#314) (#333) thanks @philipshirk

### BUG FIXES

* fix for ncdc functions to fail better - NOAA was returning HTML on request failures instead of JSON - catch that better and give proper http status code response (#338)
* `meteo_nearby_stations()` fix: coerce input data.frame to the function to a data.frame before remainder of steps - in case user inputs a tibble (#340)
* `coops_search()` fix: when `product=predictions`, we get no metadata back - so just dont adjust times  (#342)
* `lcd()` changes: gains `lcd_cache` for managing cached files; use a new internal function for safely reading each csv file, with more informative error messages;  (#344)
* `meteo_pull_monitors()` fix: changed internals of `meteo_tidy_ghcnd()` to set -9999 values to NA slightly differently to avoi failing (#348)


rnoaa 0.9.5
===========

### BUG FIXES

* `lcd()` function was unfortunately pulling data from `https://www.ncei.noaa.gov/data/global-hourly/access` - whereas it should have been pulling data from `https://www.ncei.noaa.gov/data/local-climatological-data/access` - fixed now; additionaly, `lcd_cleanup` is defunct because lcd data coming from the appropriate link has all variable names spelled out and data split up (#334) thanks @sayon000 !
* all `gefs*` functions are now defunct - they are being taken out for now until fixed - see the issues for the details (#335) (#336)


rnoaa 0.9.4
===========

### NEW FEATURES

* new gefs function helpers `gefs_dimensions` and `gefs_ensembles` (#327) (#328)

### MINOR IMPROVEMENTS

* `gefs` function fixes: fixed failing test on CRAN having to do with a date mismatch; `gefs` now cleans up temporary files  (#327) (#328)

### BUG FIXES

* Some argo buoy functions use an API and some use an FTP server. The API is down, and no longer exists. The funitons that use the API (`argo_search`, `argo_files`, `argo_qwmo`, `argo_plan`) no longer work, while the functions that use the FTP server still work (`argo_buoy_files`, `argo`)  (#333)


rnoaa 0.9.2
===========

### MINOR IMPROVEMENTS

* gefs gains new parameters `ens` and `time`, that will eventually replace the deprecated parameters `ens_idx` and `time_idx`  (#321) (#324)
* `isd()` now using fetching data using http instead of ftp

### BUG FIXES

* fix to `tornadoes()`: the URL had changed yet again (#322) (#323) thanks @mbjoseph !
* fix to gefs, was failing with some examples (#320) (#321)


rnoaa 0.9.0
===========

### NEW FEATURES

* gains `sea_ice_tabular()` function for fetching tabular .csv sea ice files instead of using the shp files in `sea_ice()` (#194)
* `seaice()` fxn name has changed to `sea_ice()` (#313)
* `sea_ice()` gains option to fetch GeoTIFF format data in addition to shp files (#219) (#313)
* gains new function `lcd_cleanup()` - takes output of call to `lcd()`, parsing additional columns that contain comma separated strings (#283)

### MINOR IMPROVEMENTS

* update README to link to ncdf4 pkg instead of ncdf pkg, and a note about which functions in rnoaa use ncdf4 (b/c ncdf4 is in Suggests) (#299) thanks @denrou
* now using markdown docs (#301)
* update `isd()` docs to highlight that cached files downloaded with the fxn will be used until deleted by the user! See `?isd` docs for details (#205)
* lat/lon param definition in `gefs` only mentioned longitude, now both vars discussed (#317) (#318)
* improve docs for `ncdc()` regarding units, and in readme and vignette as well (#265) (#315) from @amoeba

### BUG FIXES

* fix bug in `cpc_prcp()`: should have allowed dates back to 1948, but only allowed back to 1979 (#300)
* fix to `buoy()` fxn: datasets that did not have lat/lon variables were failing to be parsed by the fxn; now when lat/lon vars missing, we just give. back ncdf4 object for the user to deal with themselves (#303) (#304)
* fix to `gefs()`: longitude on the (-180, 180) scale worked but not on the (0,360) scale (#316) (#318) (#319)
* fix to `tornadoes()`: the URL for the data had changed (#311) (#312) thanks @mbjoseph
* `ncdc()` parameters `startdate`/`enddate` weren't handling dates as input values; now handle date and character inputs  (#307)
* fixed issue in `ghcnd_stations()`; there was an encoding issue with the data returned from NOAA (#305)


rnoaa 0.8.4
===========

### US federal government shutdown

This very long US federal government shutdown has allowed time for building in nicer failure behavior and more documentation for government shutdowns. There's a number of related changes:

* better government shutdown failure behavior for the `swdi()` function (#298)
* better government shutdown failure behavior for the `lcd()` function (#295)
* better failure behavior for all `ncdc*()` functions (#293) (#297)
* added a package level manual file section "Where data comes from and government shutdowns"

### MINOR IMPROVEMENTS

* `swdi()`: changed from downloading data with `download.file` to `crul` (#298)
* fix `arc2()` tests to not have hard-coded dates (#294)


rnoaa 0.8.2
===========

### BUG FIXES

* improvements in failing well when there's a US government shutdown for many functions that work with web REST APIs (#293) (#295)
* fix to `arc2` tests to not be sensitive to the real year that the test is run in, reported in CRAN checks and via email (#294)


rnoaa 0.8.0
===========

### NEW FEATURES

* gains function `bsw()` for Blended Sea Winds data (#246)
* gains function `ghcnd_read()` - to read .dly files directly, e.g., files already downloaded (#223) thanks @shabbychef for the feature request
* gains function `lcd()` for Local Climatological Data (#212)
* gains functions `se_data()` and `se_files()` for the Storm Events Database  (#282)
* `ghcnd()` and `ghcnd_search()` gain a `refresh` parameter to refresh data for the query even if it's already cached locally. in addition, these functions now print messages to tell the user what file path the data is locally cached in, and the min and max dates when using `ghcnd_search()` (#269) thanks @kgmccann
* `ncdc()` gains `add_units` parameter (boolean) to toggle adding units to the output data.frame. default is `add_units=FALSE`. if `add_units=TRUE` we match dataset id and data type id and return units if we have them. do be in touch if you see a problem with these units! `ncdc()` now returns tibbles in the `data` slot (#233) (#266) (#289) (#287)

### MINOR IMPROVEMENTS

* spelling fixes to `coops` docs (#228) thanks @jsta
* fix in `ghcnd_search()` to arrange data by day instead of by month (#247) thanks @asrivas3 for reporting
* move `swdi()` to use `xml2` and `crul` packages instead of `XML` and `httr` (#275)
* added `swdi()` tests (#239) thanks @kevin-ht-ho
* `isd_stations_search()` no longer renames lat and lon column names (#238) thanks @kevin-ht-ho
* add codemeta keywords to description (#287)
* replace `httr` with `crul` throughout package (#186)
* many tests use `vcr` now for caching, more to do waiting on `vcr` being able to handle direct to disk use cases and binary files like pdfs (#284)
* fix links in Code of Conduct

### BUG FIXES

* fixes to `gefs()`. was incorrectly repeating time values within ensembles when it should have repeated time values across ensembles so that each ensemble has a time value for each time period (#230) (#231) thanks for report from @lcsma and fix by @potterzot
* fix bug in `ncdc()` - fix to internal function `parse_ncdc()`, which was failing on `strsplit()` call if attributes was `NULL` (#232) thanks for reporting @andypicke
* fix to `cpc_prcp()`: URLs were changed at some point, fixes for this (#242) 
* fix to `cpc_prcp()`: to read `.gz` files correctly with gzfile instead of file (#248)
* fix in `homr()` - NOAA server gives back a 200 OK response even if that's not the case - but we can check content type to see if there was likely an error (#250)
* fix to `autoplot.meteo_coverage` (#258)
* fix to `buoy_stations()`: was getting the wrong station ids - fixed; and use async HTTP requests to get data faster (#261) thanks @johnharley for the bug report
* fix to `ghcnd_stations()`: was returning an extra empty row  (#267) thanks @joeroe for the bug report
* related to (#267), `ghcnd()` was giving a trailing row of NA's - fixed (#270)
* `swdi()`: `radius` parameter doesn't work, update docs to tell users not to use it (#243)
* `ncdc()` fix: fix to internal function `parse_ncdc()`, errored when more than 1 flag returned (#279) thanks @ghaines3 for the bug report
* buoy fixes (#251)
* `storm_shp()` fix: update to use new URL patterns (#263)
* make `swdi()` fail better (#274) thanks @OrionDarley
* fix `meteo_nearby_stations()` to coerce character to numeric lat and lon values (#257) thanks @mondorescue for the bug report
* fix to `meteo_nearby_stations()` - internally ignored the user supplied column names for lat and lon (#286) thanks @ghaines3 for the bug report
* fixed `tornadoes()` for Windows OS's: `utils::untar()` was failing on windows, changed to using `utils::unzip()` (#203)
* fixed `argo_buoy_files()`: use `fill=TRUE` in the `read.table` call - was erroring on some Windows OS's (#235) thanks @jonmcalder for reporting


rnoaa 0.7.0
===========

Note that some NOAA datasets have changed names:

* `GHCNDMS` is now `GSOM` (Global Summary of the Month)
* `ANNUAL` is now `GSOY` (Global Summary of the Year)

### NEW FEATURES

* `isd()` gains new parameters `additional` to toggle whether the 
non-mandatory ISD fields (additional + remarks) are parsed and 
returned & `force` to toggle whether download new version or use
cached version. `isd_read()` gains new parameter `additional` 
(see description above) (#190)
* New function for Climate Prediction Center data: `cpc_prcp()` (#193)
* New function `arc2()` to get data from Africa Rainfall Climatology 
version 2 (#201)

### MINOR IMPROVEMENTS

* A number of NOAA services now use `https` - changed internal code 
to use `https` from `http` for coops, swdi, ersst, and tornadoes
data sources (#187)
* Changes to sea ice URLs - just internal (#185)
* Fixes to `coops_search()` to handle requests better: only certain
date combinations allowed for certain COOPS products (#213) (#214) 
thanks @tphilippi !
* Now using `hoardr` package to manage caching in some functions. Will
roll out to all functions that cache soon (#191)
* README img location fix requested by CRAN (#207)
* `GHCNDMS` is now `GSOM` and `ANNUAL` is now `GSOY` - added to docs
and examples of using GSOM and GSOY (#189)

### BUG FIXES

* A number of fixes to `isd()` (#168)
* Fixes to `coops_search()` to fix time zone problems (#184) thanks @drf5n
* Fixes to `ghcnd()` - fix some column types that were of inappropriate
type before (#211)
* Fix to `ghcnd()`: we were coercing factors to integers, which caused
nonsense output - first coercing to character now, then integer (#221)
* There were problems in parsing flags (attributes) for some datasets 
via `ncdc()` function. Added metadata to the package to help parse 
flags (#199)


rnoaa 0.6.6
===========

### NEW FEATURES

* `isd()` now using a new package `isdparser` to parse
NOAA ISD files. We still fetch the file within `rnoaa`, but the
file parsing is done by `isdparser` (#176) (#177) (#180) thanks @mrubayet
for the push

### MINOR IMPROVEMENTS

* Fixed precipitation units in docs for `meteo_*` functions (#178)
thanks @mrubayet

### BUG FIXES

* Fixed bug in `ghcnd()` where internal unexported function
was not found (#179)
* Fix to `isd_stations()` and `isd_stations_search()` to work
correctly on Windows (#181) thanks @GuodongZhu
* Changed base URL for all NOAA NCDC functions (those starting with
`ncdc`) to `https` from `http` (#182) thanks @maspotts
* Changed base URL for all NOAA HOMR functions (those starting with
`homr`) to `https` from `http` (#183)


rnoaa 0.6.5
===========

### MINOR IMPROVEMENTS

* Added notes to docs of functions that do file caching - where
to find cached files.
* `meteo_clear_cache` gains parameter `force` to control `force`
parameter in `unlink()`
* Removed `lubridate` usage in `seaiceurls()` function, just using
base R functions.

### BUG FIXES

* Fixed bug which was affecting binary installs only. We accidentally
determined a path on package build, such that the user
of the CRAN binary build machine got inserted into the path.
This is now fixed. (#173)


rnoaa 0.6.4
===========

### NEW FEATURES

* New function `isd_read()` to read ISD output from `isd()` manually
instead of letting `isd()` read in the data. This is useful when you
use `isd()` but need to read the file in later when it's already cached.
(#169)
* Some functions in `rnoaa` cache files that are downloaded from
various NOAA web services. File caching is usually done when data comes
from FTP servers. In some of these functions where we cache data, we used
to write to your home directory, but have now changed all these functions
to write to a proper cache directory in a platform independent way.
We determine the cache directory using `rappdirs::user_cache_dir()`.
Note that this may change your workflow if you'd been depending on
cached files to be a in particular place on your file system. In addition,
the `path` parameter in the changed functions is now defunct, but you
get an informative warning about it (#171)

### MINOR IMPROVEMENTS

* `storm_data()` now returns a tibble/data.frame not inside of a list. We used
to return a list with a single slot `data` with a data.frame, but this was
unnecessary.
* `ghcnd_stations()` now outputs a data.frame (`tbl_df`) by itself,
instead of a data.frame nested in a list. This may change how
you access data from this function. (#163)
* Improved docs on token usage for NCDC functions (with prefix
`ncdc_*()`) (#167)
* Added note to `isd()` docs that when you get an error similar to
`Error: download failed for ftp://ftp.ncdc.noaa.gov/pub/data/noaa/1955/011490-99999-1955.gz`,
the file does not exist on NOAA's ftp servers. If your internet is down,
you'll get a different error saying as much (#170)

rnoaa 0.6.0
===============

### NEW FEATURES

* A large PR was merged with a suite of functions. Most functions added
a prefixed with `meteo_*`, and are meant to find weather monitors near
locations (`meteo_nearby_stations`), find all monitors within a radius
of a location (`meteo_distance`), calculate the distances between a
location and all available stations (`meteo_process_geographic_data`),
calculate the distance between two locations (`meteo_spherical_distance`),
pull GHCND weather data for multiple weather monitors (`meteo_pull_monitors`),
create a tidy GHCND dataset from a single monitor (`meteo_tidy_ghcnd`),
and determine the "coverage" for a station data frame (`meteo_coverage()`).
In addition, `vis_miss()` added to visualize missingness in a data.frame. See
the [PR diff against master](https://github.com/ropensci/rnoaa/pull/159/files)
for all the changes. (#159) Thanks a ton to @geanders _et al_. (@hrbrmstr,
@maelle, @jdunic, @njtierney, @leighseverson, @RyanGan, @mandilin, @jferreri,
@cpatrizio88, @ryan-hicks, @Ewen2015, @mgutilla, @hakessler, @rodlammers)

### MINOR IMPROVEMENTS

* `isd_stations_search()` changed internal structure. We replaced
usage of `geojsonio` and `lawn` for faster `dplyr::filter` for
bbox inputs, and `meteo_distance()` for `lat/long/radius` inputs
. This speeds up this function significantly. Thanks to @lukas-rokka
(#157)
* `isd_stations_search()` and `isd_stations()` now return
tibble's instead of data.frame's
* Removed cached ISD stations dataset within package to reduce
package size. Only change is now that on first use of the function
the user has to download the entire thing, but on subsquent
uses it will pull from the cached version on the users machine.
`isd_stations_search()` now caches using `rappdirs` (#161)
* Convert all `is()` uses to `inherits()`

### BUG FIXES

* Fixed `seaiceeurls()` function that's used to generate urls for
the `seaice()` function - due to change in NOAA urls (#160)
* Fix to function `ghncd_split_vars()` to not fail on `dplyr::contains`
call (#156) thanks @lawinslow !

rnoaa 0.5.6
===============

### MINOR IMPROVEMENTS

* Fixes for new `httr` version to call encoding explicitly (#135)
* Fix to broken link for reference to source code used in `gefs` functions (#121)
* Speed ups implemented for the `isd()` function - it's a time consuming task
as we have to parse a nasty string of characters line by line - more speed
ups to come in future versions (#146)
* Replace `dplyr::rbind_all()` with `dplyr::bind_rows()` as the former is
being deprecated (#152)

### BUG FIXES

* Fix for `isd()` function - was failing on some station names that had
leading zeros. (#136)
* Fix for `ncdc_stations()` - used to allow more than one station id to
be passed in, but internally only handled one. This is a restriction
due to the NOAA NCDC API. Documentation now shows an example of how
to deal with many station ids (#138)
* Fixes to the suite of `ncdc_*()` functions to allow multiple inputs
to those parameters where allowed (#139)
* Fixed bug in `ncdc_plot()` due to new `ggplot2` version (#153)
* Fixed bugs in `argo()` functions: a) with new `httr`, box input of a vector
no longer works, now manually make a character vector; b) errant file param
being passed into the http request, removed (#155)

rnoaa 0.5.2
===============

### NEW FEATURES

* New data source added: ARGO buoy data. See functions starting with `argo()` (#123)
for more, see http://www.argo.ucsd.edu/
* New data source added: CO-OPS tide and current data. See function `coops_search()`
(#111) for idea from @fmichonneau (#124) for implementing @jsta 
also (#126) (#128)

### MINOR IMPROVEMENTS

* `rgdal` moved to Suggests to make usage easier (#125)
* Changes to `ncdc_plot()` - made default brakes to just default to what
`ggplot2` does, but you can still pass in your own breaks (#131)

rnoaa 0.5.0
===============

### NEW FEATURES

* New data source added: NOAA Global Ensemble Forecast System (GEFS) data.
See functions `gefs()`, `gefs_dimension_values()`, `gefs_dimensions()`, `gefs_latitudes()`,
`gefs_longitudes()`, and `gefs_variables()` (#106) (#119)  thanks @potterzot - he's
now an author too
* New data source added: NOAA Extended Reconstructed Sea Surface Temperature
(ERSST) data. See function `ersst()` (#96)
* New function `isd_stations()` to get ISD station data.
* Added code of conduct to code repository

### MINOR IMPROVEMENTS

* Swapped `ncdf` package for `ncdf4` package. Windows binaries weren't
availiable for `ncdf4` prior to now. (#117)
* Proper license info added for javascript modules used inside the
package (#116)
* Improvements to `isd()` function to do transformations of certain
variables to give back data that makes more sense (#115)
* `leaflet`, `geojsonio`, and `lawn` added in Suggests, used in a few
functions.
* Note added to `swdi()` function man page that the `nldn` dataset is
available to military users only (#107)

### BUG FIXES

* Fix to `buoy()` function to accept character class inputs for the
`buoyid` parameter. the error occurred because matching was not
case-insensitive, now works regardless of case (#118)
* Fixes for new `ggplot2` version (#113)
* Built in `GET` request retries for `ghncd` functions as
some URLs fail unpredictably (#110)

rnoaa 0.4.2
===============

### MINOR IMPROVEMENTS

* Explicitly import non-base R pkg functions, so importing from `utils`, `methods`, and `stats` (#103)
* All NCDC legacy API functions are now defunct. See `?rnoaa-defunct` for more information (#104)
* `radius` parameter removed from `ncdc_stations()` function (#102), was already removed internally within the function in the last version, now not in the function definition, see also (#98) and (#99)
* Dropped `plyr` and `data.table` from imports. `plyr::rbind.fill()` and `data.table::rbindlist()` replaced with `dplyr::bind_rows()`.

### BUG FIXES

* Fixed problem with `httr` `v1` where empty list not allowed to pass to
the `query` parameter in `GET` (#101)

rnoaa 0.4.0
===============

### NEW FEATURES

+ Gains a suite of new functions for working with NOAA GHCND data, including
`ghcnd()`, `ghcnd_clear_cache()`, `ghcnd_countries()`, `ghcnd_search()`, `ghcnd_splitvars()`
`ghcnd_states()`, `ghcnd_stations()`, and `ghcnd_version()` (#85) (#86) (#87) (#88) (#94)
+ New contributor Adam Erickson (@DougFirErickson)
+ All NOAA buoy functions put back into the package. They were previously
on a separate branch in the GitHub repository. (#37) (#71) (#100)

### MINOR IMPROVEMENTS

+ Minor adjustments to `isd()` functions, including better man file.
+ Cleaner package imports - importing mostly only functions used in dependencies.
+ Startup message gone.
+ `callopts` parameter changed to `...` in function `swdi()`.
+ More robust test suite.
+ `ncdc()` requires that users do their own paging - previously this was done internally (#77)
+ Many dependencies dropped, simplifying package: `RCurl`, `maptools`, `stringr`, `digest`.
A few new ones added: `dplyr`, `tidyr`.

### DEPRECATED AND DEFUNCT

+ All `erddap` functions now defunct - see the package [rerddap](https://github.com/ropensci/rerddap),
a general purpose R client for ERDDAP servers. (#51) (#73) (#90) (#95)
+ The `extent` function in `noaa_stations()` used to accept either a bounding
box or a point defined by lat/long. The lat/long option dropped as it required
two packages, one of which is a pain to install for many users (#98) (#99)

rnoaa 0.3.3
===============

### NEW FEATURES

+ New data source NOAA legacy API with ISD, daily, and ish data via function
`ncdc_legacy()`. (#54)
+ New function `isd()` to get ISD data from NOAA FTP server. (#76)
+ ERDDAP gridded data sets added. Now tabledap datasets are accessible via
`erddap_table()`, while gridded datasets are available via `erddap_grid()`. Helper
function `erddap_search()` was modified to search for either tabledap or griddap
datasets, and `erddap_info()` gets and prints summary information differently
for tabledap and griddap datasets. (#63)

### MINOR IMPROVEMENTS

+ `erddap_data()` defunct, now as functions `erddap_table()` and `erddap_grid()`, uses new
`store` parameter which takes a function, either `disk(path, overwrite)` to store
on disk or `memory()` to store in R memory.
+ `assertthat` library removed, replaced with `stopifnot()`

rnoaa 0.3.0
===============

### NEW FEATURES

+ New data source added (NOAA torndoes data) via function `tornadoes()`. (#56)
+ New data source added (NOAA storm data from IBTrACS) via functions
`storm_*()`. (#57)
+ New data source added (NOAA weather station metadata from HOMR) via functions
`homr_*()` (#59)
+ New vignettes for storm data and homr data.
+ Some functions in rnoaa now print data.frame outputs as `dplyr`-like outputs
with a summary of the data.frame, as appropriate.

### MINOR IMPROVEMENTS

+ Across all `ncdc_*` functions changed `callopts` parameter to `...`. This parameter
allow you to pass in options to `httr::GET` to modify curl requests. (#61)
+ A new helper function `check_key()` looks for one of two stored keys, as an
environment variable under the name `NOAA_KEY`, or an option variable under the name
`noaakey`. Environment variables can be set during session like `Sys.setenv(VAR = "...")`,
or stored long term in your `.Renviron` file. Option variables can be set during session
like `options(var = "...")`, or stored long term in your `.Rprofile` file.
+ `is.*` and `print.*` functions no longer have public man files, but can be seen via
`rnoaa:::` if needed.

rnoaa 0.2.0
===============

### NEW FEATURES

* New package imports: `sp`, `rgeos`, `assertthat`, `jsonlite`, and `ncdf4`, and new package Suggests: `knitr`, `taxize`
* Most function names changed. All `noaa*()` functions for NCDC data changed to `ncdc*()`. `noaa_buoy()` changed to `buoy()`. `noaa_seaice()` changed to `seaice()`. When you call the old versions an error is thrown, with a message pointing you to the new function name. See ?rnoaa-defunct.
* New vignettes: NCDC attributes, NCDC workflow, Seaice vignette, SWDI vignette, ERDDAP vignette, NOAA buoy vignette.
* New functions to interact with NOAA ERDDAP data: `erddap_info()`, `erddap_data()`, and `erddap_search()`.
* New functions to interact with NOAA buoy data: `buoy()`, including a number of helper functions.
* `ncdc()` now splits apart attributes. Previously, the attributes were returned as a single column, but now there is column for each attribute so data can be easily retrieved. Attribute columns differ for each different `datasetid`.
* `buoy()` function has been removed from the CRAN version of `rnoaa`. Install the version with `buoy()` and associated functions via `devtools::install_github("ropensci/rnoaa", ref="buoy")`

### MINOR IMPROVEMENTS

* `noaa_swdi()` (function changed to `swdi()`) gains new parameter `filepath` to specify path to write a file to if `format=kmz` or `format=shp`. Examples added for using `format=` csv, shp, and kmz.
* Now using internal version of `plyr::compact`.
* Added API response checker/handler to all functions to pass on helpful messages on server errors.
* `ncdc()` gains new parameter `includemetadata`. If TRUE, includes metadata, if not, does not, and response should be faster as does not take time to calculate metadata.
* `noaa_stations()` gains new parameter `radius`. If `extent` is a vector of length 4 (for a bounding box) then radius is ignored, but if you pass in two points to `extent`, it is interpreted as a point, and then `radius` is used as the distance upon which to construct a bounding box. `radius` default is 10 km.

### BUG FIXES

* `datasetid`, `startdate`, and `enddate` are often required parameters, and changes were made to help users with this.


rnoaa 0.1.0
===============

### NEW FEATURES

* Submitted to CRAN.


rnoaa 0.0.8
===============

### NEW FEATURES

* Wrote new functions for NOAA API v2.
* A working vignette now.


rnoaa 0.0.1
===============

### NEW FEATURES

* Wrappers for NOAA API v1 were written, not on CRAN at this point.
## Test environments

* ubuntu 20.04 (local install), R 4.1.2
* macOS-latest (release, on GitHub Actions), R 4.1.2
* macOS-latest (devel, on GitHub Actions), R-dev
* ubuntu-latest (release, on GitHub Actions), R 4.1.2
* win-builder (release, devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

-----

This version updates the writing locations to comply with CRAN policy by switching from using rappdirs::user_cache_dir("rnoaa/ersst") to using tools::R_user_dir("rnoaa/ersst", which = "cache").

The previous release had missed a use of rappdirs::user_cache_dir in the tests.

Thanks!
Daniel Hocking
<!--- Provide a general summary of your changes in the Title above 
Note: Continuous integration checks may fail because secret tokens for the NCDC API are not available in PR's. Don't worry about these failures UNLESS they are clearly due to code in your PR.
-->

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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/rnoaa/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rnoaa.git`
* Make sure to track progress upstream (i.e., on our version of `rnoaa` at `ropensci/rnoaa`) by doing `git remote add upstream https://github.com/ropensci/rnoaa.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/rnoaa`

Note: Continuous integration checks may fail because secret tokens for the NCDC API are not available in PR's. Don't worry about these failures UNLESS they are clearly due to code in your PR.

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Email

Do not email. Open an issue.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<!-- IF YOU DO NOT INCLUDE YOUR SESSION INFO YOUR ISSUE WILL LIKELY BE CLOSED WITHOUT FURTHER CONSIDERATION -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
## revdepcheck results

We checked 2 reverse dependencies (0 from CRAN + 2 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                                |
|:--------|:------------------------------------|
|version  |R version 4.1.2 (2021-11-01)         |
|os       |Ubuntu 20.04.3 LTS                   |
|system   |x86_64, linux-gnu                    |
|ui       |RStudio                              |
|language |(EN)                                 |
|collate  |en_US.UTF-8                          |
|ctype    |en_US.UTF-8                          |
|tz       |America/New_York                     |
|date     |2021-11-27                           |
|rstudio  |2021.09.1+372 Ghost Orchid (desktop) |
|pandoc   |2.5 @ /usr/bin/pandoc                |

# Dependencies

|package      |old      |new      |Δ  |
|:------------|:--------|:--------|:--|
|rnoaa        |1.3.7    |1.3.8    |*  |
|cli          |3.1.0    |3.1.0    |   |
|colorspace   |2.0-2    |2.0-2    |   |
|cpp11        |0.4.1    |0.4.1    |   |
|crayon       |1.4.2    |1.4.2    |   |
|crul         |1.2.0    |1.2.0    |   |
|curl         |4.3.2    |4.3.2    |   |
|data.table   |1.14.2   |1.14.2   |   |
|digest       |0.6.28   |0.6.28   |   |
|dplyr        |1.0.7    |1.0.7    |   |
|ellipsis     |0.3.2    |0.3.2    |   |
|fansi        |0.5.0    |0.5.0    |   |
|farver       |2.1.0    |2.1.0    |   |
|generics     |0.1.1    |0.1.1    |   |
|geonames     |0.999    |0.999    |   |
|ggplot2      |3.3.5    |3.3.5    |   |
|glue         |1.5.0    |1.5.0    |   |
|gridExtra    |2.3      |2.3      |   |
|gtable       |0.3.0    |0.3.0    |   |
|hoardr       |0.5.2    |0.5.2    |   |
|httpcode     |0.3.0    |0.3.0    |   |
|isdparser    |0.4.0    |0.4.0    |   |
|isoband      |0.2.5    |0.2.5    |   |
|jsonlite     |1.7.2    |1.7.2    |   |
|labeling     |0.4.2    |0.4.2    |   |
|lifecycle    |1.0.1    |1.0.1    |   |
|lubridate    |1.8.0    |1.8.0    |   |
|magrittr     |2.0.1    |2.0.1    |   |
|mime         |0.12     |0.12     |   |
|munsell      |0.5.0    |0.5.0    |   |
|pillar       |1.6.4    |1.6.4    |   |
|pkgconfig    |2.0.3    |2.0.3    |   |
|purrr        |0.3.4    |0.3.4    |   |
|R6           |2.5.1    |2.5.1    |   |
|rappdirs     |0.3.3    |0.3.3    |   |
|RColorBrewer |1.1-2    |1.1-2    |   |
|Rcpp         |1.0.7    |1.0.7    |   |
|rjson        |0.2.20   |0.2.20   |   |
|rlang        |0.4.12   |0.4.12   |   |
|scales       |1.1.1    |1.1.1    |   |
|tibble       |3.1.6    |3.1.6    |   |
|tidyr        |1.1.4    |1.1.4    |   |
|tidyselect   |1.1.1    |1.1.1    |   |
|triebeard    |0.3.0    |0.3.0    |   |
|urltools     |1.7.3    |1.7.3    |   |
|utf8         |1.2.2    |1.2.2    |   |
|vctrs        |0.3.8    |0.3.8    |   |
|viridisLite  |0.4.0    |0.4.0    |   |
|withr        |2.4.2    |2.4.2    |   |
|XML          |3.99-0.8 |3.99-0.8 |   |
|xml2         |1.3.2    |1.3.2    |   |

# Revdeps

## Failed to check (2)

|package |version |error |warning |note |
|:-------|:-------|:-----|:-------|:----|
|wildviz |?       |      |        |     |
|Z10     |?       |      |        |     |

*Wow, no problems at all. :)*# wildviz

<details>

* Version: 
* GitHub: https://github.com/ropensci/rnoaa
* Source code: NA
* Number of recursive dependencies: 0

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
# Z10

<details>

* Version: 
* GitHub: https://github.com/ropensci/rnoaa
* Source code: NA
* Number of recursive dependencies: 0

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
rnoaa
=====

```{r echo=FALSE}
library("knitr")
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 10,
  fig.path = "man/figures/",
  cache.path = "inst/cache/"
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rnoaa)](https://cranchecks.info/pkgs/rnoaa)
[![R-check](https://github.com/ropensci/rnoaa/workflows/R-check/badge.svg)](https://github.com/ropensci/rnoaa/actions)
[![codecov.io](https://codecov.io/github/ropensci/rnoaa/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rnoaa?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rnoaa?color=C9A115)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rnoaa)](https://cran.r-project.org/package=rnoaa)


`rnoaa` is an R interface to many NOAA data sources. We don't cover all of them, but we include many commonly used sources, and add we are always adding new sources. We focus on easy to use interfaces for getting NOAA data, and giving back data in easy to use formats downstream. We currently don't do much in the way of plots or analysis. To get started see: https://docs.ropensci.org/rnoaa/articles/rnoaa.html

## Data sources in rnoaa

* NOAA NCDC climate data:
    * We are using the NOAA API version 2
    * Docs for the NCDC API are at https://www.ncdc.noaa.gov/cdo-web/webservices/v2
    * GHCN Daily data is available at http://www.ncdc.noaa.gov/ghcn-daily-description via FTP and HTTP
* Severe weather data docs are at https://www.ncdc.noaa.gov/swdiws/
* Sea ice data (ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/shapefiles)
* NOAA buoy data (https://www.ndbc.noaa.gov/)
* ERDDAP data (https://upwell.pfeg.noaa.gov/erddap/index.html)
  * Now in package rerddap (https://github.com/ropensci/rerddap)
* Tornadoes! Data from the NOAA Storm Prediction Center (https://www.spc.noaa.gov/gis/svrgis/)
* HOMR - Historical Observing Metadata Repository (http://www.ncdc.noaa.gov/homr/api)
* GHCND FTP data (ftp://ftp.ncdc.noaa.gov/pub/data/noaa) - NOAA NCDC API has some/all (not sure really) of this data, but FTP allows to get more data more quickly
* Extended Reconstructed Sea Surface Temperature (ERSST) data (https://www.ncdc.noaa.gov/data-access/marineocean-data/extended-reconstructed-sea-surface-temperature-ersst-v4)
* Argo buoys (http://www.argo.ucsd.edu/) - a global array of more than 3,000 free-drifting profiling floats that measures thetemperature and salinity of the upper 2000 m of the ocean
* NOAA CO-OPS - tides and currents data (https://tidesandcurrents.noaa.gov/)
* NOAA Climate Prediction Center (CPC) (http://www.cpc.ncep.noaa.gov/)
* Africa Rainfall Climatology version 2 (ftp://ftp.cpc.ncep.noaa.gov/fews/fewsdata/africa/arc2/ARC2_readme.txt)
* Blended Sea Winds (https://www.ncdc.noaa.gov/data-access/marineocean-data/blended-global/blended-sea-winds)
* Local Climatological Data (https://www.ncdc.noaa.gov/cdo-web/datatools/lcd)
* Storm Events Database (https://www.ncdc.noaa.gov/stormevents/)

## Help/Getting Started

Documentation is at https://docs.ropensci.org/rnoaa/, and there are many vignettes in the package itself, available in your R session, or on CRAN (https://cran.r-project.org/package=rnoaa). The tutorials:

* **Getting started - start here**
* NOAA Buoy vignette
* NOAA National Climatic Data Center (NCDC) vignette (examples)
* NOAA NCDC attributes vignette
* NOAA NCDC workflow vignette
* Sea ice vignette
* Severe Weather Data Inventory (SWDI) vignette
* Historical Observing Metadata Repository (HOMR) vignette

## netcdf data

Some functions use netcdf files, including:

* `ersst`
* `buoy`
* `bsw`
* `argo`
 
You'll need the `ncdf4` package for those functions, and those only. `ncdf4` is in Suggests in this package, meaning you only need `ncdf4` if you are using any of the functions listed above. You'll get an informative error telling you to install `ncdf4` if you don't have it and you try to use the those functions. Installation of `ncdf4` should be straightforward on any system. See https://cran.r-project.org/package=ncdf4

## NOAA NCDC Datasets

There are many NOAA NCDC datasets. All data sources work, except `NEXRAD2` and `NEXRAD3`, for an unknown reason. This relates to `ncdc_*()` functions only.

```{r echo=FALSE}
library('rnoaa')
dat <- ncdc_datasets()$data
dat <- dat[, !names(dat) %in% 'uid']
dat <- dat[, c('id', 'name', 'mindate', 'maxdate', 'datacoverage')]
names(dat) <- c('Dataset', 'Description', 'Start Date', 'End Date', 'Data Coverage')
knitr::kable(dat)
```

```{r echo=FALSE}
cat(paste0("table updated on ", Sys.Date()))
```

**NOAA NCDC Attributes**

Each NOAA dataset has a different set of attributes that you can potentially get back in your search. See https://www.ncdc.noaa.gov/cdo-web/datasets for detailed info on each dataset. We provide some information on the attributes in this package; see the vignette for attributes (https://docs.ropensci.org/rnoaa/articles/ncdc_attributes.html) to find out more


## Contributors

* Scott Chamberlain (https://github.com/sckott)
* Daniel Hocking (https://github.com/djhocking)
* Brooke Anderson (https://github.com/geanders)
* Maëlle Salmon (https://github.com/maelle)
* Adam Erickson (https://github.com/adam-erickson)
* Nicholas Potter (https://github.com/potterzot)
* Joseph Stachelek (https://github.com/jsta)

## Meta

* Please report any issues or bugs: https://github.com/ropensci/rnoaa/issues
* License: MIT
* Get citation information for `rnoaa` in R doing `citation(package = 'rnoaa')`
* Please note that this package is released with a Contributor Code of Conduct (https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rnoaa sparklines}
-->

rnoaa sparklines example
======

### Plot data from many stations

Load some libraries

```{r}
library(rnoaa)
library(scales)
library(lubridate)
library(maptools)
library(ggplot2)
library(doMC)
library(ggsubplot)
library(maps)
library(plyr)
library(stringr)
```

Find stations to get data from

```{r}
stations <- noaa_stations(datasetid = "GHCND", enddate = "2012-12-01", limit = 120)
res <- stations$data$id
```

Get data from stations. Note in the code below that we are using parallelization. You do not have to do this, you can simply set `.parallel=FALSE`, or equivalently, use `lapply` instead of `llply` (`llply` is from the `plyr` package).

```{r}
noaa_fw <- failwith(NULL, noaa)
registerDoMC(cores = 4)
dat <- llply(res, function(x) noaa_fw(datasetid = "GHCND", datatypeid = 'PRCP', stationid = x, startdate = '2010-06-01', enddate = '2010-09-30'), .parallel = TRUE)
dat <- compact(dat)
length(dat)
```

Make a `data.frame` and fix dates. 

```{r}
df <- ldply(dat, function(x) x$data)
df$date <- ymd(str_replace(as.character(df$date), "T00:00:00\\.000|T00:00:00", ""))
```

Get station lat and long data so that we can put data on a map.

```{r}
latlongs <- llply(res, function(x) 
  noaa_stations(x, datasetid = "GHCND")$data$meta[c("id", "latitude", "longitude")])
latlongs <- ldply(latlongs, function(x) as.data.frame(x))
df2 <- merge(df, latlongs, by.x = "station", by.y = "id")
head(df2)
```

Make a map

```{r}
world_map <- map_data("world")
p <- ggplot() + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "gray40", size = 0.2) + 
  annotate(geom = "text", x = -155, y = -55, 
           label = sprintf("Max value is\n %s mm", max(df2$value)/10))
p + 
  geom_subplot(aes(longitude, latitude, group = station, 
                   subplot = geom_line(aes(date, value)), size = 1), 
               ref = ref_vline(aes(fill = length(value)), thickness = 0.1), 
               width = rel(2), height = rel(5), data = df2) + 
  theme(legend.position = "none")
```---
title: "ropenaq and rnoaa"
subtitle: "Complementing air quality data with weather data using rnoaa"
author: "Maëlle Salmon"
date: "2020-07-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{ropenaq and rnoaa}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



# Introduction: getting air quality data

This vignette aims at explaining how you can complement a data.frame with weather data using rnoaa. In this vignette we shall use air quality data from the OpenAQ platform queried with the ropenaq package, for India. Using ropenaq (https://github.com/ropensci/ropenaq) one can get e.g. PM2.5 values over time in Indianapolis over the last days. For getting all data for march we'll loop over several pages.

First, we need to know how many measures are available for Indianapolis for March 2016.


```r
library("ropenaq")

measurementsIndianapolis <- aq_measurements(city = "Indianapolis", parameter = "pm25",
                                     date_from = as.character(Sys.Date() - 30),
                                     date_to = as.character(Sys.Date()),
                                     limit = 10000 )
save(measurementsIndianapolis, file = "data/measurementsIndianapolis.RData")
```

We filter negative values.


```r

load("data/measurementsIndianapolis.RData")
library("dplyr")
measurementsIndianapolis %>% head() %>% knitr::kable()
measurementsIndianapolis <- filter(measurementsIndianapolis, value > 0)
```

We now transform these data into daily data.


```r
# only keep stations with geographical information
measurementsIndianapolis <- filter(measurementsIndianapolis, !is.na(latitude))
# now transform to daily data
measurementsIndianapolis <- measurementsIndianapolis %>%
  #mutate(day = as.Date(dateLocal)) %>%
  group_by(location, day) %>%
  summarize(value = mean(value),
            longitude = longitude[1],
            latitude = latitude[1]) %>%
  ungroup()
measurementsIndianapolis %>% head() %>% knitr::kable()

```

Air quality and weather are correlated, so one could be interested in getting a time series of say temperature for the same location. The OpenAQ platform itself does not provide weather data but nearly all stations have geographical coordinates. Our goal here will be to use rnoaa to complement this table with precipitation and temperature.

# Find weather stations

For finding the right station(s) we shall use the `meteo_nearby_stations` function. It returns a list with the weather stations nearby each latitude/longitude given as arguments, respecting the other arguments such as maximal radius, first year with data, etc. For finding stations one might have to play a bit with the parameters until there is at least one station for each location.

Here we query stations with a less than 15km distance from the air quality stations, with precipitation and temperature data, and with data starting from 2016. Note that the query takes a while.


```r
library("rnoaa")
station_data <- ghcnd_stations()
lat_lon_df <- select(measurementsIndianapolis,
                     location,
                     latitude,
                     longitude) %>% unique() %>%
  ungroup() %>%
  rename(id = location) %>%
  mutate(id = factor(id))

stationsIndianapolis <- meteo_nearby_stations(lat_lon_df = as.data.frame(lat_lon_df),
                                       station_data = station_data,
                                       radius = 5,
                                       year_min = as.character(format(Sys.Date(), "%Y")),
                                       var = c("TAVG", "PRCP"))
stationsIndianapolis <- unique(bind_rows(stationsIndianapolis) %>% select(- distance))

save(stationsIndianapolis, file = "data/stationsIndianapolis.RData")

```


```r
load("data/stationsIndianapolis.RData")
stationsIndianapolis %>% knitr::kable()
```

Now let us plot the AQ and weather stations on a quick and dirty map with no legend, red for AQ stations, blue for weather stations.

> Not shown


```r
library("ggmap")
map <- get_map(location = "Indianapolis", zoom = 11)
ggmap(map) +
  geom_point(aes(x = longitude, y = latitude),
             data = stationsIndianapolis, col = "blue", size = 4)+
  geom_point(aes(x = longitude, y = latitude),
             data = measurementsIndianapolis, col = "red", size = 4)
```

# Query weather data for these stations

For pulling weather data from these weather monitors, we shall use the `meteo_pull_monitors` function.


```r
library("rnoaa")
monitors <- stationsIndianapolis$id
all_monitors_clean <- meteo_pull_monitors(monitors,
    date_min = as.character(Sys.Date() - 30),
    date_max = as.character(Sys.Date())) %>%
  dplyr::rename(day = date, location = id)
all_monitors_clean %>% head() %>% knitr::kable()
```

Here we notice some values are not available. Therefore, we might need to go back to weather stations searching with, for instance, a larger radius. In this case let's say we're ok with the result of the search.

# Join the two tables, thus complementing the original table

Therefore, in this case we will bind the rows of the air quality table with the weather table.


```r
measurementsIndianapolis <- bind_rows(measurementsIndianapolis, all_monitors_clean)
measurementsIndianapolis %>% head() %>% knitr::kable()
```


 Now some locations are air quality locations and have only missing values in the weather columns, and some locations are weather locations and have only missing values in the air quality columns.

We can plot the data we got.


```r
data_plot <- measurementsIndianapolis %>%
  rename(pm25 = value) %>%
  select(- longitude, - latitude)

data_plot <- tidyr::gather_(data_plot, "parameter", "value",
  names(data_plot)[3:ncol(data_plot)])

library("ggplot2")
ggplot(data_plot) +
  geom_line(aes(x = day, y = value, col = location)) +
  facet_grid(parameter ~ ., scales = "free_y")
```
---
title: "working with buoy data"
author: "Scott Chamberlain"
date: "2020-07-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{working with buoy data}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



This vignette covers NOAA buoy data from the National Buoy Data Center. The
main function to get data is `buoy`, while `buoys` can be used to
get the buoy IDs and web pages for each buoy.



```r
library('rnoaa')
```

## Find out what buoys are available in a dataset


```r
res <- buoys(dataset = "cwind")
```

Inspect the buoy ids, and the urls for them


```r
head(res)
#>      id
#> 1 41001
#> 2 41002
#> 3 41004
#> 4 41006
#> 5 41008
#> 6 41009
#>                                                                        url
#> 1 https://dods.ndbc.noaa.gov/thredds/catalog/data/cwind/41001/catalog.html
#> 2 https://dods.ndbc.noaa.gov/thredds/catalog/data/cwind/41002/catalog.html
#> 3 https://dods.ndbc.noaa.gov/thredds/catalog/data/cwind/41004/catalog.html
#> 4 https://dods.ndbc.noaa.gov/thredds/catalog/data/cwind/41006/catalog.html
#> 5 https://dods.ndbc.noaa.gov/thredds/catalog/data/cwind/41008/catalog.html
#> 6 https://dods.ndbc.noaa.gov/thredds/catalog/data/cwind/41009/catalog.html
```

Or browse them on the web


```r
browseURL(res[1, 2])
```

## Get buoy data

With `buoy` you can get data for a particular dataset, buoy id, year, and datatype. 

Get data for a buoy

> if no year or datatype specified, we get the first file


```r
buoy(dataset = 'cwind', buoyid = 46085)
#> Using c2007.nc
#> Dimensions (rows/cols): [33486 X 5] 
#> 2 variables: [wind_dir, wind_spd] 
#> 
#> # A tibble: 33,486 x 5
#>    time                   lat   lon wind_dir wind_spd
#>    <chr>                <dbl> <dbl>    <int>    <dbl>
#>  1 2007-05-05T02:00:00Z  55.9 -143.      331     2.80
#>  2 2007-05-05T02:10:00Z  55.9 -143.      328     2.60
#>  3 2007-05-05T02:20:00Z  55.9 -143.      329     2.20
#>  4 2007-05-05T02:30:00Z  55.9 -143.      356     2.10
#>  5 2007-05-05T02:40:00Z  55.9 -143.      360     1.5 
#>  6 2007-05-05T02:50:00Z  55.9 -143.       10     1.90
#>  7 2007-05-05T03:00:00Z  55.9 -143.       10     2.20
#>  8 2007-05-05T03:10:00Z  55.9 -143.       14     2.20
#>  9 2007-05-05T03:20:00Z  55.9 -143.       16     2.10
#> 10 2007-05-05T03:30:00Z  55.9 -143.       22     1.60
#> # … with 33,476 more rows
```

Including year


```r
buoy(dataset = 'cwind', buoyid = 41001, year = 1999)
#> Using c1999.nc
#> Dimensions (rows/cols): [52554 X 5] 
#> 2 variables: [wind_dir, wind_spd] 
#> 
#> # A tibble: 52,554 x 5
#>    time                   lat   lon wind_dir wind_spd
#>    <chr>                <dbl> <dbl>    <int>    <dbl>
#>  1 1999-01-01T00:00:00Z  34.7 -72.7      272    11.7 
#>  2 1999-01-01T00:10:00Z  34.7 -72.7      260    11   
#>  3 1999-01-01T00:20:00Z  34.7 -72.7      249     8.70
#>  4 1999-01-01T00:30:00Z  34.7 -72.7      247     8.40
#>  5 1999-01-01T00:40:00Z  34.7 -72.7      240     7.10
#>  6 1999-01-01T00:50:00Z  34.7 -72.7      242     7.90
#>  7 1999-01-01T01:00:00Z  34.7 -72.7      246     8.30
#>  8 1999-01-01T01:10:00Z  34.7 -72.7      297    10.9 
#>  9 1999-01-01T01:20:00Z  34.7 -72.7      299    11.3 
#> 10 1999-01-01T01:30:00Z  34.7 -72.7      299    11.1 
#> # … with 52,544 more rows
```

Including year and datatype


```r
buoy(dataset = 'cwind', buoyid = 45005, year = 2008, datatype = "c")
#> Dimensions (rows/cols): [29688 X 5] 
#> 2 variables: [wind_dir, wind_spd] 
#> 
#> # A tibble: 29,688 x 5
#>    time                   lat   lon wind_dir wind_spd
#>    <chr>                <dbl> <dbl>    <int>    <dbl>
#>  1 2008-04-29T09:00:00Z  41.7 -82.4       10     9   
#>  2 2008-04-29T09:10:00Z  41.7 -82.4        8     9   
#>  3 2008-04-29T09:20:00Z  41.7 -82.4        5     9.30
#>  4 2008-04-29T09:30:00Z  41.7 -82.4       13     9.5 
#>  5 2008-04-29T09:40:00Z  41.7 -82.4       14     9.40
#>  6 2008-04-29T09:50:00Z  41.7 -82.4       12     9.40
#>  7 2008-04-29T14:00:00Z  41.7 -82.4      341     6.5 
#>  8 2008-04-29T14:10:00Z  41.7 -82.4      332     6.80
#>  9 2008-04-29T14:20:00Z  41.7 -82.4      335     6.40
#> 10 2008-04-29T14:30:00Z  41.7 -82.4      332     6.5 
#> # … with 29,678 more rows
```

Including just datatype


```r
buoy(dataset = 'cwind', buoyid = 45005, datatype = "c")
#> Using c1996.nc
#> Dimensions (rows/cols): [26784 X 5] 
#> 2 variables: [wind_dir, wind_spd] 
#> 
#> # A tibble: 26,784 x 5
#>    time                   lat   lon wind_dir wind_spd
#>    <chr>                <dbl> <dbl>    <int>    <dbl>
#>  1 1996-05-15T23:00:00Z  41.7 -82.4      337     2.20
#>  2 1996-05-15T23:10:00Z  41.7 -82.4      282     1   
#>  3 1996-05-15T23:20:00Z  41.7 -82.4      282     2.20
#>  4 1996-05-15T23:30:00Z  41.7 -82.4      258     2.60
#>  5 1996-05-15T23:40:00Z  41.7 -82.4      254     3   
#>  6 1996-05-15T23:50:00Z  41.7 -82.4      252     2.70
#>  7 1996-05-16T00:00:00Z  41.7 -82.4      240     2.10
#>  8 1996-05-16T00:10:00Z  41.7 -82.4      246     2.40
#>  9 1996-05-16T00:20:00Z  41.7 -82.4      251     2.70
#> 10 1996-05-16T00:30:00Z  41.7 -82.4      253     2.90
#> # … with 26,774 more rows
```
---
title: "NOAA NCDC dataset attributes"
author: "Scott Chamberlain"
date: "2020-10-05"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NOAA NCDC dataset attributes}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



The attributes, or "flags", for each row of the output for data may have a flag with it.
Each `datasetid` has it's own set of flags. The following are flag columns, and what they
stand are. `fl_` is the beginning of each flag column name, then one or more characters
to describe the flag, keeping it short to maintain a compact data frame. Some of these
fields are the same across datasetids, but they may have different possible values. See
below details on each dataset.

* fl_c = completeness
* fl_m = measurement
* fl_d =  day
* fl_q = quality
* fl_s = source
* fl_t = time
* fl_cmiss = consecutive missing
* fl_miss = missing
* fl_u = units

__Datasets__

* [GHCND](#ghcnd)
* [GSOM](#gsom)
* [GSOY](#gsoy)
* [NORMAL_ANN](#normalann)
* [NORMAL_DLY](#normaldly)
* [NORMAL_HLY](#normalhly)
* [NORMAL_MLY](#normalmly)
* [PRECIP_HLY](#preciphly)
* [PRECIP_15](#precip15)
* NEXRAD2 - not working yet
* NEXRAD3 - not working yet


### <a href="#ghcnd" name="ghcnd"/>#</a> Dataset: GHCND

#### flm_m (Measurement flag)

* __Blank:__ no measurement information applicable
* __A:__ value in precipitation or snow is a multi-day total, accumulated since last measurement (used on Daily Form pdf file)
* __B:__ precipitation total formed from two 12-hour totals
* __D:__ precipitation total formed from four six-hour totals
* __K:__ converted from knots
* __L:__ temperature appears to be lagged with respect to reported hour of observation
* __O:__ converted from oktas
* __P:__ identified as "missing presumed zero" in DSI 3200 and 3206
* __T:__ trace of precipitation, snowfall, or snow depth
* __W:__ converted from 16-point WBAN code (for wind direction)

#### fl_q (Quality flag)

* __Blank:__ did not fail any quality assurance check
* __D:__ failed duplicate check
* __G:__ failed gap check
* __I:__ failed internal consistency check
* __K:__ failed streak/frequent-value check
* __L:__ failed check on length of multiday period
* __M:__ failed mega-consistency check
* __N:__ failed naught check
* __O:__ failed climatological outlier check
* __R:__ failed lagged range check
* __S:__ failed spatial consistency check
* __T:__ failed temporal consistency check
* __W:__ temperature too warm for snow
* __X:__ failed bounds check
* __Z:__ flagged as a result of an official Datzilla investigation

#### fl_s (Source flag)

* __Blank:__ No source (i.e., data value missing)
* __0:__ U.S. Cooperative Summary of the Day (NCDC DSI-3200)
* __6:__ CDMP Cooperative Summary of the Day (NCDC DSI-3206)
* __7:__ U.S. Cooperative Summary of the Day, Transmitted via WxCoder3 (NCDC DSI-3207)
* __A:__ U.S. Automated Surface Observing System (ASOS) real-time data (since January 1, 2006)
* __a:__ Australian data from the Australian Bureau of Meteorology
* __B:__ U.S. ASOS data for October 2000-December 2005 (NCDC DSI-3211)
* __b:__ Belarus update
* __E:__ European Climate Assessment and Dataset (Klein Tank et al., 2002)
* __F:__ U.S. Fort data
* __G:__ Official Global Climate Observing System (GCOS) or other government-supplied data
* __H:__ High Plains Regional Climate Center real-time data
* __I:__ International collection (non U.S. data received through personal contacts)
* __K:__ U.S. Cooperative Summary of the Day data digitized from paper observer forms (from 2011 to present)
* __M:__ Monthly METAR Extract (additional ASOS data)
* __N:__ Community Collaborative Rain, Hail,and Snow (CoCoRaHS)
* __Q:__ Data from several African countries that had been "quarantined", that is, withheld from public release until permission was granted from the respective meteorological services
* __R:__ NCDC Reference Network Database (Climate Reference Network and Historical Climatology Network-Modernized)
* __r:__ All-Russian Research Institute of Hydrometeorological Information-World Data Center
* __S:__ Global Summary of the Day (NCDC DSI-9618) NOTE: "S" values are derived from hourly synoptic reports exchanged on the Global Telecommunications System (GTS). Daily values derived in this fashion may differ significantly from "true" daily data, particularly for precipitation(i.e., use with caution).
* __u:__ Ukraine update
* __W:__ WBAN/ASOS Summary of the Day from NCDC's Integrated Surface Data (ISD).
* __X:__ U.S. First-Order Summary of the Day (NCDC DSI-3210)
* __Z:__ Datzilla official additions or replacements
* __z:__ Uzbekistan update

#### fl_t (Time of observation flag)

Is the (2 digit hour, 2 digit minute) 24 hour clock time of the observation given as the
local time at the station of record.



### <a href="#gsom" name="gsom"/>#</a> Dataset: GSOM

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/gsom-gsoy.pdf

Observations are synonymous with elements or values, and defined in Table A below. 9999 indicates missing data or data that has not been received.

__flags:__ Missing flag , Consecutive Missing flag

#### fl_miss (Missing flag)

Defined as total number of days observation/element is missing in that month. This can  be taken as a measure of quality or completeness as the higher the number of days sampled in the month, the more representative the value is for the entire month.

#### fl_cmiss (Consecutive missing flag)

Defined as the maximum number of consecutive days in the month that an  observation/element is missing.



### <a href="#gsoy" name="gsoy"/>#</a> Dataset: GSOY

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/gsom-gsoy.pdf

#### fl_m (Measurement flag)

* __A:__ Accumulated amount. This value is a total that may include data from a previous month or months
(TPCP).
* __B:__ Adjusted Total. Monthly value totals based on proportional available data across the entire month.
(CLDD, HTDD)
* __E:__ An estimated monthly or annual total.
* __I:__ Monthly means or totals based on incomplete time series. 1 to 9 days are missing. (MMNT,MMXP,
MMXT, MNTM, TPCP, TSNW)
* __M:__ used to indicate data element missing.
* __S:__ Precipitation for the amount is continuing to be accumulated. Total will be included in a
__subsequent value (TPCP). Example:__ Days 1-20 had 1.35 inches of precipitation, then a period of accumulation began. The element TPCP would then be 00135S and the total accumulated amount value appears in a subsequent monthly value. If TPCP = 0 there was no precipitation measured during the month. flag 1 is set to "S" and the total accumulated amount appears in a subsequent monthly value.
* __T:__ Trace of precipitation, snowfall, or snow depth. The precipitation data value will = "00000".
(EMXP, MXSD, TPCP, TSNW)
* __+:__ The phenomena in question occurred on several days. The date in the DAY field is the last day
of occurrence.
* __Blank:__ No report

#### fl_q (Quality flag)

* __A:__ Accumulated amount
* __E:__ Estimated value
* __+:__ Value occurred on more than one day, last date of occurrence is used

#### fl_d (Number of days flag )

Number of days is given as 00 when all days in the month are considered in
computing data value or otherwise the maximum number of consecutive days in the month considered
in computing the data value.

#### fl_u (Units flag)

* __C:__ Whole degree Celsius
* __D:__ Whole Fahrenheit Degree Day
* __F:__ Whole degree Fahrenheit
* __HI:__ Hundredths of inches
* __I:__ Whole inches
* __M:__ Whole miles
* __MH:__ Miles per hour
* __MM:__ Millimeters
* __NA:__ No units applicable (dimensionless)
* __TC:__ Tenths of degrees Celsius
* __TF:__ Tenths of degrees Fahrenheit
* __TI:__ Tenths of inches
* __TM:__ Tenths of millimeters
* __1:__ Soils, degrees Fahrenheit, soil depths in inches and hundredths
* __2:__ Soils, degrees Celsius, soil depth in whole centimeters
* __3:__ Soils, degrees Celsius, soil, soil depth in inches and hundredths
* __4:__ Soils, degrees Fahrenheit, soil depth in whole centimeters
* __5:__ Soils, If the soil station closed during the current month, '5' indicates the station has closed.



### <a href="#normalann" name="normalann"/>#</a> Dataset: NORMAL_ANN

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/NORMAL_ANN_documentation.pdf

#### Description
The 1981-2010 Normals comprise all climate normals using the thirty year period of temperature,
degree days, precipitation, snowfall, snow depth, wind, etc. Data is organized into hourly, daily,
monthly, seasonal and annual. This document describes the elements and layout of the Seasonal and
Annual Normals which are derived from a composite of climate records from numerous sources that
were merged and then subjected to a suite of quality assurance reviews.

#### fl_c (Completeness flag)

flags accompany every Normals value and indicate the completeness of the data record used to
compute each value, accounting for methodological differences for different product classes. There are
six flag options described generally below. Due to methodological differences, the flags are
applied somewhat differently between the temperature-based normals and the precipitation-based
normals. For the precipitation-based and hourly normals, the following flags were assigned
independently for each normals value reported based on number of years available for that individual
calculation. For temperature-based normals, strong precedence is given to the monthly normals of
maximum and minimum temperature or derived from the flags for these two variables.

* __C:__ complete (all 30 years used)
* __S:__ standard (no more than 5 years missing and no more than 3 consecutive
 years missing among the sufficiently complete years)
* __R:__ representative (observed record utilized incomplete, but value was scaled
 or based on filled values to be representative of the full period of record)
* __P:__ provisional (at least 10 years used, but not sufficiently complete to be
 labeled as standard or representative). Also used for parameter values on
 February 29 as well as for interpolated daily precipitation, snowfall, and
 snow depth percentiles.
* __Q:__ quasi-normal (at least 2 years per month, but not sufficiently complete to
 be labeled as provisional or any other higher flag code. The associated
 value was computed using a pseudonormals approach or derived from monthly
 pseudonormals.
* __Blank:__ the data value is reported as a special value (see section B under III. Additional Information
below).

__Note:__ flags Q and R aren't applicable to average number of days with different precipitation,
snowfall, and snow depth threshold exceedance; precipitation/snowfall/snow
probabilities of occurrence. Further, Q flags are not applicable for standard deviations.

### <a href="#normaldly" name="normaldly"/>#</a> Dataset: NORMAL_DLY

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/NORMAL_DLY_documentation.pdf

#### Description
The 1981-2010 Normals comprise all climate normals using the thirty year period of temperature,
degree days, precipitation, snowfall, snow depth, wind, etc. Data is organized into hourly, daily,
monthly, seasonal and annual. This document describes the elements and layout of the Daily Normals
which are derived from a composite of climate records from numerous sources that were merged and
then subjected to a suite of quality assurance reviews.


#### fl_c (Completeness flag)

Same as NORMAL_ANN, see the description above at [Completeness flag](#completenessflag).


### <a href="#normalhly" name="normalhly"/>#</a> Dataset: NORMAL_HLY

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/NORMAL_HLY_documentation.pdf

#### Description
The 1981-2010 Normals comprise all climate normals using the thirty year period of temperature,
degree days, precipitation, snowfall, snow depth, wind, etc. Data is organized into hourly, daily,
monthly, seasonal and annual normals. This document describes the elements and layout of the Hourly
Normals which are derived from a composite of climate records from numerous sources that were
merged and then subjected to a suite of quality assurance reviews.
The hourly normals provide a suite of descriptive statistics based on hourly observations at a few
hundred stations from across the United States and its Pacific territories. Statistics are provided as 30-
year averages, frequencies of occurrence, and percentiles for each hour and day of the year. These
products are useful in examination of the diurnal change of a particular variable.

For temperature, dew point and mean sea level pressure an average hourly value as well as a 10th and
90th percentile of hourly values is given. For heating and cooling degree hours, an average hourly value is
given using a 65 degree F base. Average hourly values are also given for heat index and wind chill. Cloud
cover statistics include percent frequency of clear, few, scattered, broken and overcast conditions. Wind
statistics include prevailing and secondary wind direction and percent frequency, average wind speed,
percentage of calm winds and mean wind vector direction and magnitude.

The statistics are computed from the ISD-lite dataset. 262 stations were selected from
the ISD-lite data, based on their completeness and membership in a list of what were known as "first
order stations." These are typically airport locations with the needed 24 hours/day observations to
make hourly normals meaningful. All stations had at least 27 of the 30 years represented.

Each hourly normal is computed on the basis of 450 possible values. This is the aggregation of the value
for a particular date and time, plus and minus 7 days, over each of 30 years. If fewer than 350 valid
values are present, the output is given as the special value 9999. No normals are computed for February
29, but data for February 29 is included in the 15 day window for leap years. The original data has been
shifted from Greenwich Mean Time to an end product in local standard time.

#### fl_c (Completeness flag)

Same as NORMAL_ANN, see the description above at [Completeness flag](#completenessflag).

### <a href="#normalmly" name="normalmly"/>#</a> Dataset: NORMAL_MLY

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/NORMAL_MLY_documentation.pdf

#### Description

The 1981-2010 Normals comprise all climate normals using the thirty year period of temperature,
degree days, precipitation, snowfall, snow depth, wind, etc. Data are organized into hourly, daily,
monthly, seasonal and annual. This document describes the elements and layout of the Monthly
Normals which are derived from a composite of climate records from numerous sources that were
merged and then subjected to a suite of quality assurance reviews.

#### fl_c (Completeness flag)

Same as NORMAL_ANN, see the description above at [Completeness flag](#completenessflag).

### <a href="#preciphly" name="preciphly"/>#</a> Dataset: PRECIP_HLY

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/PRECIP_HLY_documentation.pdf

#### Description

Hourly Precipitation Data (labeled Precipitation Hourly in Climate Data Online system) is
a database that gives time-sequenced hourly precipitation amounts for a network of
over 7000 reporting station located primarily in the United States. Data is collected from
a variety of sources including National Weather Service reporting stations, volunteer
cooperative observers, Federal Aviation Administration (FAA), utility companies, etc.

Concerning rain gages/data processing: Data from weighing rain gages, Fischer-Porter
gages, Universal rain gages and in recent years, more modern measuring equipment in
conjunction with automated recording sites, etc. have been used in this dataset over the
period. Precipitation values have been checked and edited as necessary by both
automated and manual methods. Because of some inconsistencies identified with the
earlier data (prior to 1996), historical data were reprocessed in 1997. This rehabilitated
data covered 53 million observations between 1900 and 1995. Similar quality control
checks are in place that maintain consistency between the historical and operationally
received data.

#### fl_m (Measurement flag)

__Note:__ This field is left blank when no flag is needed.

* __a:__ Begin accumulation. A data value of 99999 accompanies this flag. It indicates that the
accumulation has begun at some time during the hour.

* __A:__ End accumulation (an amount is associated with this flag). It indicates that
accumulation has ended sometime during the hour. Accumulated period indicates that
the precipitation amount is correct, but only the exact beginning and ending times are
known. A data value of 99999 occurring on the last day and hour of a month indicates
the accumulation continues into the next month.

* __, (comma):__ Used at the beginning of a data month when an accumulation is in progress
from the previous month. A data value of 99999 always accompanies this flag. This flag
is used prior to 1984.

* __{ :__ Begin deleted period during the hour (inclusive). The original data were received, but
were unreadable or clearly recognized as noise. A value of 99999 accompanies this
flag. Primarily used since 1984. Also used in Alaska for 1976-1978.

* __}:__ End deleted period during the hour (inclusive). The original data were received, but
were unreadable or clearly recognized as noise. A value of 99999 accompanies this
flag. Primarily used since 1984. Also used in Alaska for 1976-1978.

* __[:__ Begin missing period during the hour (inclusive). A value of 99999 accompanies this
flag.

* __]:__ End missing period during the hour (inclusive). A value of 99999 accompanies this
flag. Prior to 1984 if precipitation occurred during the last hour of the missing period, the
ending missing value appears with a non-zero value. Beginning in 1984, the beginning
and ending hours of the missing period are recorded as "99999[" and "99999],"
respectively.. A missing flag indicates that the data were not received. The flag appears
on the first and last day of each month for which data were not received or not
processed by NCDC.

* __E:__ Evaporation may have occurred. Data may or may not be reliable. This flag was used
during the period 1984-1993.

* __g:__ Only used for day 1, hour 0100, when precipitation amount is zero.

* __T:__ Indicates a "trace" amount. Data value with this will be zero. "T" flags appear on
National Weather Service data only since July 1996.

* __M:__ Missing data. No data available for this period.

#### fl_q (Data quality flag)

* __Z:__ Indicates probable amounts as a result of melting frozen precipitation. This flag may
be used to identify those sites that are deficient in the manner the snow shields are  employed.
Used since January 1996.
* __R:__ This data value failed one of the NCDC's quality control tests.
* __Q:__ Pre 1996 usage, Indicates value failed an extreme value test (value will be present).
Data are to be used with caution. Extreme tests used are 1) value was not an accumulated amount and was higher than
the one-hour statewide 100 year return period precipitation amount or 2) if they value
was an accumulated amount and was higher than the 24 hour statewide extreme
precipitation total.
* __Q:__ 1996 to present usage. A single erroneous value (value will be present). Rarely
used since 1996.
* __q:__ An hourly value excludes one or more 15 minute periods. Lowest data resolution is
15 minutes. Used since January 1996.


### <a href="#precip15" name="precip15"/>#</a> Dataset: PRECIP_15

__More info:__ https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/PRECIP_15_documentation.pdf

#### Description

15 Minute Precipitation Data (labeled Precipitation 15 Minute in Climate Data Online
system) is a database that gives time-sequenced quarter-hour precipitation amounts for
a network of over 3600 reporting station located primarily in the United States. Data is
collected from a variety of sources including National Weather Service reporting
stations, volunteer cooperative observers, Federal Aviation Administration (FAA), utility
companies, etc.

Concerning rain gages/data processing: Data from Fischer-Porter gages between May
1971 and December 1983 have been used in this dataset. Precipitation values have
been checked and edited as necessary by both automated and manual methods. Data
processing procedures were updated in January 1984 to produce the element
structured data base files and further enhanced beginning with the January 1996 data
month. Currently, interactive quality control procedures are in place that has added
many checks and features and data are subjected to automated editing procedures that
reduce the manual handling of the data.

#### fl_m (Data measurement flag)

__QPCP__

__Note:__ This field is left blank when no flag is needed.

* __a:__ Begin accumulation. A data value of 99999 accompanies this flag. It indicates that the
accumulation has begun at some time during the 15 minute period.
* __A:__ End accumulation (an amount is associated with this flag). It indicates that
accumulation has ended sometime during the 15 minute period. Accumulated period
indicates that the precipitation amount is correct, but only the exact beginning and
ending times are known. A data value of 99999 occurring on the last day and hour of a
month indicates the accumulation continues into the next month.
* __, (comma):__ Used at the beginning of a data month when an accumulation is in progress.
This flag is used prior to 1984.
* __{ :__ Begin deleted period during the hour (inclusive).
* __}:__ End deleted period during the hour (inclusive).
* __[:__ Begin missing period during the hour (inclusive).
* __]:__ End missing period during the hour (inclusive).
* __E:__ Evaporation may have occurred. Data may or may not be reliable. This flag was used
during the period 1984-1993.
* __g:__ Only used on day 1 when precipitation amount is zero.
* __T:__ Indicates a "trace" amount. Data value with this will be zero. "T" flags appear on
National Weather Service data only.
* __M:__ Missing data. No data available for this period.


__QGAG__

* __a:__ begin accumulation (indicates measurement periods overlapped)
* __A:__ end accumulation
* __[:__ begin missing
* __]:__ end missing
* __{:__ begin delete
* __}:__ end delete
* __S:__ gage reset

#### fl_q (Data quality flag)

__QPCP__

* __X:__ Used for data prior to 1996 as part of a 1997 data rehabilitation effort. Indicates value
failed an extreme value test; data are to be used with caution. Extreme tests were: 1) if
the value was not an accumulated precipitation total, the value failed the one hour
statewide 100 year return period precipitation and 2) if the value was an accumulated
precipitation total, the value failed the 24 hour statewide extreme precipitation total.
* __Z:__ Indicates probable amounts as a result of melting frozen precipitation. This flag may
be used to identify those sites that are deficient in the manner the snow shields are
employed. Used since January 1996.
* __R:__ This data value failed one of NCDC's quality control tests.
* __Q:__ A single erroneous value (value will be present). Used since January 1996.
* __q:__ An hourly value excludes one or more 15 minute periods. Lowest data resolution is
15 minutes. Used since January 1996.
* __A:__ Accumulated period and amount. An accumulated period indicates that the
precipitation amount is correct, but the exact beginning and ending times are only
known to the extent that the precipitation occurred sometime within the accumulation
period.

__QGAG__

* __Q:__ Questionable value. Data not used.
* __P:__ Punched mechanism failure, missing punch assumed. Assumed punch value being
used.
* __V:__ Evaporation likely. Gage value has dropped. Data are being used.


#### fl_u (Units flag)

HI indicates data values (QGAG or QPCP) are in hundredths of inches. HT indicates data values (QGAG or QPCP) are in tenths of inches.
---
title: "NCDC introduction"
author: "Scott Chamberlain"
date: "2020-08-14"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NCDC introduction}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



`rnoaa` is an R wrapper for many NOAA data types, including National Climatic Data Center (NCDC).

## Load rnoaa


```r
library('rnoaa')
```

## Get info on a station by specifying a datasetid, locationid, and stationid


```r
ncdc_stations(datasetid='GHCND', locationid='FIPS:12017', stationid='GHCND:USC00084289')
#> $meta
#> NULL
#> 
#> $data
#>   elevation    mindate    maxdate latitude                  name datacoverage
#> 1      17.7 1899-01-01 2020-08-12 28.80286 INVERNESS 3 SE, FL US            1
#>                  id elevationUnit longitude
#> 1 GHCND:USC00084289        METERS -82.31266
#> 
#> attr(,"class")
#> [1] "ncdc_stations"
```

## Search for data and get a data.frame


```r
out <- ncdc(datasetid='NORMAL_DLY', datatypeid='dly-tmax-normal', startdate = '2010-05-01', enddate = '2010-05-10')
out$data
#> # A tibble: 25 x 5
#>    date                datatype        station           value fl_c 
#>    <chr>               <chr>           <chr>             <int> <chr>
#>  1 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:AQW00061705   869 C    
#>  2 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:CAW00064757   607 Q    
#>  3 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:CQC00914080   840 R    
#>  4 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:CQC00914801   858 R    
#>  5 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:FMC00914395   876 P    
#>  6 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:FMC00914419   885 P    
#>  7 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:FMC00914446   885 P    
#>  8 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:FMC00914482   868 R    
#>  9 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:FMC00914720   899 R    
#> 10 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:FMC00914761   897 P    
#> # … with 15 more rows
```

Note that the `value` column has strangely large numbers for temperature measurements.
By convention, `rnoaa` doesn't do any conversion of values from the APIs and some APIs use seemingly odd units.

You have two options here:

1. Use the `add_units` parameter on `ncdc()` to have `rnoaa` attempt to look up the units. This is a good idea to try first.

2. Consult the documentation for whiechever dataset you're accessing. In this case, `GHCND` has a README (https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt) which indicates `TMAX` is measured in tenths of degrees Celcius.

### See a `data.frame` with units

As mentioned above, you can use the `add_units` parameter with `ncdc()` to ask `rnoaa` to attempt to look up units for whatever data you ask it to return.
Let's ask `rnoaa` to add units to some precipitation (PRCP) data:


```r
with_units <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-05-01', enddate = '2010-10-31', limit=500, add_units = TRUE)
head( with_units$data )
#> # A tibble: 6 x 9
#>   date            datatype station         value fl_m  fl_q  fl_so fl_t  units  
#>   <chr>           <chr>    <chr>           <int> <chr> <chr> <chr> <chr> <chr>  
#> 1 2010-05-01T00:… PRCP     GHCND:USW00014…     0 "T"   ""    0     2400  mm_ten…
#> 2 2010-05-02T00:… PRCP     GHCND:USW00014…    30 ""    ""    0     2400  mm_ten…
#> 3 2010-05-03T00:… PRCP     GHCND:USW00014…    51 ""    ""    0     2400  mm_ten…
#> 4 2010-05-04T00:… PRCP     GHCND:USW00014…     0 "T"   ""    0     2400  mm_ten…
#> 5 2010-05-05T00:… PRCP     GHCND:USW00014…    18 ""    ""    0     2400  mm_ten…
#> 6 2010-05-06T00:… PRCP     GHCND:USW00014…    30 ""    ""    0     2400  mm_ten…
```
From the above output, we can see that the units for `PRCP` values are "mm_tenths" which means tenths of a millimeter.
You won't always be so lucky and sometimes you will have to look up the documentation on your own.

## Plot data, super simple, but it's a start


```r
out <- ncdc(datasetid='NORMAL_DLY', stationid='GHCND:USW00014895', datatypeid='dly-tmax-normal', startdate = '2010-01-01', enddate = '2010-12-10', limit = 300)
ncdc_plot(out)
```

![plot of chunk six](../man/figures/six-1.png)

Note that `PRCP` values are in units of tenths of a millimeter, as we found out above.

## More on plotting

### Example 1

Search for data first, then plot


```r
out <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-05-01', enddate = '2010-10-31', limit=500)
```

Default plot


```r
ncdc_plot(out)
```

![plot of chunk ncdc-plot-zero](../man/figures/ncdc-plot-zero-1.png)

Create 14 day breaks


```r
ncdc_plot(out, breaks="14 days")
```

![plot of chunk ncdc-plot-0](../man/figures/ncdc-plot-0-1.png)

One month breaks


```r
ncdc_plot(out, breaks="1 month", dateformat="%d/%m")
```

![plot of chunk ncdc-plot-1](../man/figures/ncdc-plot-1-1.png)

### Example 2

Search for data


```r
out <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP',
            startdate = '2010-05-01', enddate = '2010-10-31', limit=500)
```

Make a plot, with 6 hour breaks, and date format with only hour


```r
ncdc_plot(out, breaks = "1 month", dateformat = "%d/%m")
```

![plot of chunk ncdc-plot-2](../man/figures/ncdc-plot-2-1.png)

## Combine many calls to noaa function

Search for two sets of data


```r
out1 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-03-01', enddate = '2010-05-31', limit=500)

out2 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-09-01', enddate = '2010-10-31', limit=500)
```

Then combine with a call to `ncdc_combine`


```r
df <- ncdc_combine(out1, out2)
head(df[[1]]); tail(df[[1]])
#> # A tibble: 6 x 8
#>   date                datatype station           value fl_m  fl_q  fl_so fl_t 
#>   <chr>               <chr>    <chr>             <int> <chr> <chr> <chr> <chr>
#> 1 2010-03-01T00:00:00 PRCP     GHCND:USW00014895     0 "T"   ""    0     2400 
#> 2 2010-03-02T00:00:00 PRCP     GHCND:USW00014895     0 "T"   ""    0     2400 
#> 3 2010-03-03T00:00:00 PRCP     GHCND:USW00014895     0 "T"   ""    0     2400 
#> 4 2010-03-04T00:00:00 PRCP     GHCND:USW00014895     0 ""    ""    0     2400 
#> 5 2010-03-05T00:00:00 PRCP     GHCND:USW00014895     0 ""    ""    0     2400 
#> 6 2010-03-06T00:00:00 PRCP     GHCND:USW00014895     0 ""    ""    0     2400
#> # A tibble: 6 x 8
#>   date                datatype station           value fl_m  fl_q  fl_so fl_t 
#>   <chr>               <chr>    <chr>             <int> <chr> <chr> <chr> <chr>
#> 1 2010-10-26T00:00:00 PRCP     GHCND:USW00014895   221 ""    ""    0     2400 
#> 2 2010-10-27T00:00:00 PRCP     GHCND:USW00014895     0 ""    ""    0     2400 
#> 3 2010-10-28T00:00:00 PRCP     GHCND:USW00014895     0 "T"   ""    0     2400 
#> 4 2010-10-29T00:00:00 PRCP     GHCND:USW00014895     0 "T"   ""    0     2400 
#> 5 2010-10-30T00:00:00 PRCP     GHCND:USW00014895     0 ""    ""    0     2400 
#> 6 2010-10-31T00:00:00 PRCP     GHCND:USW00014895     0 ""    ""    0     2400
```

Then plot - the default passing in the combined plot plots the data together. In this case it looks kind of weird since a straight line combines two distant dates.


```r
ncdc_plot(df)
```

![plot of chunk ncdc-plot-line](../man/figures/ncdc-plot-line-1.png)

But we can pass in each separately, which uses `facet_wrap` in `ggplot2` to plot each set of data in its own panel.


```r
ncdc_plot(out1, out2, breaks="45 days")
```

![plot of chunk ncdc-plot-panel](../man/figures/ncdc-plot-panel-1.png)
---
title: "rnoaa introduction"
author: "Scott Chamberlain"
date: "2021-11-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{rnoaa introduction}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



For additional vignettes see https://docs.ropensci.org/rnoaa/

## Installation

__GDAL__

You'll need GDAL (https://gdal.org/) installed first. You may want to use GDAL >= `0.9-1` since that version or later can read TopoJSON format files as well, which aren't required here, but may be useful. Install GDAL:

* OSX - From https://www.kyngchaos.com/software/frameworks/
* Linux - run `sudo apt-get install gdal-bin`
* Windows - From https://trac.osgeo.org/osgeo4w/

Then when you install the R package `rgdal` (`rgeos` also requires GDAL), you'll most likely need to specify where you're `gdal-config` file is on your machine, as well as a few other things. I have an OSX Mavericks machine, and this works for me (there's no binary for Mavericks, so install the source version):


```r
install.packages("https://cran.r-project.org/src/contrib/rgdal_0.9-1.tar.gz", repos = NULL, type="source", configure.args = "--with-gdal-config=/Library/Frameworks/GDAL.framework/Versions/1.10/unix/bin/gdal-config --with-proj-include=/Library/Frameworks/PROJ.framework/unix/include --with-proj-lib=/Library/Frameworks/PROJ.framework/unix/lib")
```

The rest of the installation should be easy. If not, let us know.

__Stable version from CRAN__


```r
install.packages("rnoaa")
```

__or development version from GitHub__


```r
remotes::install_github("ropensci/rnoaa")
```

__Load rnoaa__


```r
library('rnoaa')
```

## NCDC v2 API data

**NCDC Authentication**

You'll need an API key to use the NOAA NCDC functions (those starting with `ncdc*()`) in this package (essentially a password). Go to https://www.ncdc.noaa.gov/cdo-web/token to get one. *You can't use this package without an API key.*

Once you obtain a key, there are two ways to use it.

a) Pass it inline with each function call (somewhat cumbersome)


```r
ncdc(datasetid = 'PRECIP_HLY', locationid = 'ZIP:28801', datatypeid = 'HPCP', startdate = '2013-10-01', enddate = '2013-12-01', limit = 5, token =  "YOUR_TOKEN")
```

b) Alternatively, you might find it easier to set this as an option, either by adding this line to the top of a script or somewhere in your `.rprofile`


```r
options(noaakey = "KEY_EMAILED_TO_YOU")
```

c) You can always store in permamently in your `.Rprofile` file.

###  Fetch list of city locations in descending order


```r
ncdc_locs(locationcategoryid='CITY', sortfield='name', sortorder='desc')
#> $meta
#> $meta$totalCount
#> [1] 1989
#> 
#> $meta$pageCount
#> [1] 25
#> 
#> $meta$offset
#> [1] 1
#> 
#> 
#> $data
#>       mindate    maxdate                  name datacoverage            id
#> 1  1892-08-01 2021-05-31            Zwolle, NL       1.0000 CITY:NL000012
#> 2  1901-01-01 2021-11-23            Zurich, SZ       1.0000 CITY:SZ000007
#> 3  1957-07-01 2021-11-23         Zonguldak, TU       1.0000 CITY:TU000057
#> 4  1906-01-01 2021-11-23            Zinder, NG       0.9025 CITY:NG000004
#> 5  1973-01-01 2021-11-23        Ziguinchor, SG       1.0000 CITY:SG000004
#> 6  1938-01-01 2021-11-23         Zhytomyra, UP       0.9723 CITY:UP000025
#> 7  1948-03-01 2021-11-23        Zhezkazgan, KZ       0.9302 CITY:KZ000017
#> 8  1951-01-01 2021-11-23         Zhengzhou, CH       1.0000 CITY:CH000045
#> 9  1941-01-01 2021-11-23          Zaragoza, SP       1.0000 CITY:SP000021
#> 10 1936-01-01 2009-06-17      Zaporiyhzhya, UP       1.0000 CITY:UP000024
#> 11 1957-01-01 2021-11-23          Zanzibar, TZ       0.8016 CITY:TZ000019
#> 12 1973-01-01 2021-11-23            Zanjan, IR       0.9105 CITY:IR000020
#> 13 1893-01-01 2021-11-26     Zanesville, OH US       1.0000 CITY:US390029
#> 14 1912-01-01 2021-11-23             Zahle, LE       0.9819 CITY:LE000004
#> 15 1951-01-01 2021-11-23           Zahedan, IR       0.9975 CITY:IR000019
#> 16 1860-12-01 2021-11-23            Zagreb, HR       1.0000 CITY:HR000002
#> 17 1929-07-01 2021-10-09         Zacatecas, MX       1.0000 CITY:MX000036
#> 18 1947-01-01 2021-11-23 Yuzhno-Sakhalinsk, RS       1.0000 CITY:RS000081
#> 19 1893-01-01 2021-11-26           Yuma, AZ US       1.0000 CITY:US040015
#> 20 1942-02-01 2021-11-25   Yucca Valley, CA US       1.0000 CITY:US060048
#> 21 1885-01-01 2021-11-26      Yuba City, CA US       1.0000 CITY:US060047
#> 22 1998-02-01 2021-11-23            Yozgat, TU       0.9993 CITY:TU000056
#> 23 1893-01-01 2021-11-26     Youngstown, OH US       1.0000 CITY:US390028
#> 24 1894-01-01 2021-11-26           York, PA US       1.0000 CITY:US420024
#> 25 1869-01-01 2021-11-26        Yonkers, NY US       1.0000 CITY:US360031
#> 
#> attr(,"class")
#> [1] "ncdc_locs"
```

### Get info on a station by specifying a dataset, locationtype, location, and station


```r
ncdc_stations(datasetid='GHCND', locationid='FIPS:12017', stationid='GHCND:USC00084289')
#> $meta
#> NULL
#> 
#> $data
#>   elevation    mindate    maxdate latitude                  name datacoverage
#> 1      17.7 1899-01-01 2021-11-08 28.80286 INVERNESS 3 SE, FL US            1
#>                  id elevationUnit longitude
#> 1 GHCND:USC00084289        METERS -82.31266
#> 
#> attr(,"class")
#> [1] "ncdc_stations"
```


### Search for data


```r
out <- ncdc(datasetid='NORMAL_DLY', stationid='GHCND:USW00014895', datatypeid='dly-tmax-normal', startdate = '2010-05-01', enddate = '2010-05-10')
```

### See a data.frame


```r
head( out$data )
#> # A tibble: 6 × 5
#>   date                datatype        station           value fl_c 
#>   <chr>               <chr>           <chr>             <int> <chr>
#> 1 2010-05-01T00:00:00 DLY-TMAX-NORMAL GHCND:USW00014895   652 S    
#> 2 2010-05-02T00:00:00 DLY-TMAX-NORMAL GHCND:USW00014895   655 S    
#> 3 2010-05-03T00:00:00 DLY-TMAX-NORMAL GHCND:USW00014895   658 S    
#> 4 2010-05-04T00:00:00 DLY-TMAX-NORMAL GHCND:USW00014895   661 S    
#> 5 2010-05-05T00:00:00 DLY-TMAX-NORMAL GHCND:USW00014895   663 S    
#> 6 2010-05-06T00:00:00 DLY-TMAX-NORMAL GHCND:USW00014895   666 S
```

Note that the `value` column has strangely large numbers for temperature measurements.
By convention, `rnoaa` doesn't do any conversion of values from the APIs and some APIs use seemingly odd units.

You have two options here:

1. Use the `add_units` parameter on `ncdc` to have `rnoaa` attempt to look up the units. This is a good idea to try first.

2. Consult the documentation for whiechever dataset you're accessing. In this case, `GHCND` has a README (https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt) which indicates `TMAX` is measured in tenths of degrees Celcius.

### See a `data.frame` with units

As mentioned above, you can use the `add_units` parameter with `ncdc()` to ask `rnoaa` to attempt to look up units for whatever data you ask it to return.
Let's ask `rnoaa` to add units to some precipitation (PRCP) data:


```r
with_units <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-05-01', enddate = '2010-10-31', limit=500, add_units = TRUE)
head( with_units$data )
#> # A tibble: 6 × 9
#>   date                datatype station      value fl_m  fl_q  fl_so fl_t  units 
#>   <chr>               <chr>    <chr>        <int> <chr> <chr> <chr> <chr> <chr> 
#> 1 2010-05-01T00:00:00 PRCP     GHCND:USW00…     0 "T"   ""    0     2400  mm_te…
#> 2 2010-05-02T00:00:00 PRCP     GHCND:USW00…    30 ""    ""    0     2400  mm_te…
#> 3 2010-05-03T00:00:00 PRCP     GHCND:USW00…    51 ""    ""    0     2400  mm_te…
#> 4 2010-05-04T00:00:00 PRCP     GHCND:USW00…     0 "T"   ""    0     2400  mm_te…
#> 5 2010-05-05T00:00:00 PRCP     GHCND:USW00…    18 ""    ""    0     2400  mm_te…
#> 6 2010-05-06T00:00:00 PRCP     GHCND:USW00…    30 ""    ""    0     2400  mm_te…
```
From the above output, we can see that the units for `PRCP` values are "mm_tenths" which means tenths of a millimeter.
You won't always be so lucky and sometimes you will have to look up the documentation on your own.

### Plot data, super simple, but it's a start


```r
out <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-05-01', enddate = '2010-10-31', limit=500)
ncdc_plot(out, breaks="1 month", dateformat="%d/%m")
```

![plot of chunk unnamed-chunk-13](../man/figures/unnamed-chunk-13-1.png)

Note that `PRCP` values are in units of tenths of a millimeter, as we found out above.

### More plotting

You can pass many outputs from calls to the `noaa` function in to the `ncdc_plot` function.


```r
out1 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-03-01', enddate = '2010-05-31', limit=500)
out2 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-09-01', enddate = '2010-10-31', limit=500)
ncdc_plot(out1, out2, breaks="45 days")
```

![plot of chunk unnamed-chunk-14](../man/figures/unnamed-chunk-14-1.png)

### Get table of all datasets


```r
ncdc_datasets()
#> $meta
#> $meta$offset
#> [1] 1
#> 
#> $meta$count
#> [1] 11
#> 
#> $meta$limit
#> [1] 25
#> 
#> 
#> $data
#>                     uid    mindate    maxdate                        name
#> 1  gov.noaa.ncdc:C00861 1763-01-01 2021-11-25             Daily Summaries
#> 2  gov.noaa.ncdc:C00946 1763-01-01 2021-11-01 Global Summary of the Month
#> 3  gov.noaa.ncdc:C00947 1763-01-01 2021-01-01  Global Summary of the Year
#> 4  gov.noaa.ncdc:C00345 1991-06-05 2021-11-25    Weather Radar (Level II)
#> 5  gov.noaa.ncdc:C00708 1994-05-20 2021-11-24   Weather Radar (Level III)
#> 6  gov.noaa.ncdc:C00821 2010-01-01 2010-01-01     Normals Annual/Seasonal
#> 7  gov.noaa.ncdc:C00823 2010-01-01 2010-12-31               Normals Daily
#> 8  gov.noaa.ncdc:C00824 2010-01-01 2010-12-31              Normals Hourly
#> 9  gov.noaa.ncdc:C00822 2010-01-01 2010-12-01             Normals Monthly
#> 10 gov.noaa.ncdc:C00505 1970-05-12 2014-01-01     Precipitation 15 Minute
#> 11 gov.noaa.ncdc:C00313 1900-01-01 2014-01-01        Precipitation Hourly
#>    datacoverage         id
#> 1          1.00      GHCND
#> 2          1.00       GSOM
#> 3          1.00       GSOY
#> 4          0.95    NEXRAD2
#> 5          0.95    NEXRAD3
#> 6          1.00 NORMAL_ANN
#> 7          1.00 NORMAL_DLY
#> 8          1.00 NORMAL_HLY
#> 9          1.00 NORMAL_MLY
#> 10         0.25  PRECIP_15
#> 11         1.00 PRECIP_HLY
#> 
#> attr(,"class")
#> [1] "ncdc_datasets"
```

### Get data category data and metadata


```r
ncdc_datacats(locationid = 'CITY:US390029')
#> $meta
#> $meta$totalCount
#> [1] 39
#> 
#> $meta$pageCount
#> [1] 25
#> 
#> $meta$offset
#> [1] 1
#> 
#> 
#> $data
#>                     name            id
#> 1    Annual Agricultural        ANNAGR
#> 2     Annual Degree Days         ANNDD
#> 3   Annual Precipitation       ANNPRCP
#> 4     Annual Temperature       ANNTEMP
#> 5    Autumn Agricultural         AUAGR
#> 6     Autumn Degree Days          AUDD
#> 7   Autumn Precipitation        AUPRCP
#> 8     Autumn Temperature        AUTEMP
#> 9               Computed          COMP
#> 10 Computed Agricultural       COMPAGR
#> 11           Degree Days            DD
#> 12      Dual-Pol Moments DUALPOLMOMENT
#> 13             Echo Tops       ECHOTOP
#> 14      Hydrometeor Type   HYDROMETEOR
#> 15            Miscellany          MISC
#> 16                 Other         OTHER
#> 17               Overlay       OVERLAY
#> 18         Precipitation          PRCP
#> 19          Reflectivity  REFLECTIVITY
#> 20    Sky cover & clouds           SKY
#> 21   Spring Agricultural         SPAGR
#> 22    Spring Degree Days          SPDD
#> 23  Spring Precipitation        SPPRCP
#> 24    Spring Temperature        SPTEMP
#> 25   Summer Agricultural         SUAGR
#> 
#> attr(,"class")
#> [1] "ncdc_datacats"
```

## Tornado data

The function `tornadoes()` simply gets __all the data__. So the call takes a while, but once done, is fun to play with.


```r
shp <- tornadoes()
#> Error in tornadoes(): could not find function "tornadoes"
library('sp')
plot(shp)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'plot': object 'shp' not found
```

## HOMR metadata

In this example, search for metadata for a single station ID


```r
homr(qid = 'COOP:046742')
```

## Argo buoys data

There are a suite of functions for Argo data, a few egs:


```r
# Spatial search - by bounding box
argo_search("coord", box = c(-40, 35, 3, 2))

# Time based search
argo_search("coord", yearmin = 2007, yearmax = 2009)

# Data quality based search
argo_search("coord", pres_qc = "A", temp_qc = "A")

# Search on partial float id number
argo_qwmo(qwmo = 49)

# Get data
argo(dac = "meds", id = 4900881, cycle = 127, dtype = "D")
```

## CO-OPS data

Get daily mean water level data at Fairport, OH (9063053)


```r
coops_search(station_name = 9063053, begin_date = 20150927, end_date = 20150928,
             product = "daily_mean", datum = "stnd", time_zone = "lst")
#> $metadata
#> $metadata$id
#> [1] "9063053"
#> 
#> $metadata$name
#> [1] "Fairport"
#> 
#> $metadata$lat
#> [1] "41.7597"
#> 
#> $metadata$lon
#> [1] "-81.2811"
#> 
#> 
#> $data
#>            t       v   f
#> 1 2015-09-27 174.430 0,0
#> 2 2015-09-28 174.422 0,0
```

## Additional vignettes

For additional vignettes see https://docs.ropensci.org/rnoaa/
---
title: "Sea ice vignette"
author: "Scott Chamberlain"
date: "2020-07-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Sea ice vignette}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



Get sea ice data at ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/shapefiles

********************



```r
library('rnoaa')
library('dplyr')
library('ggplot2')
```

### Look at data for a series of years for Feb, South pole


```r
res <- sapply(seq(1986, 1990, 1), function(x)
    sea_ice(x, month = 'Feb', pole = 'S'))
lapply(res, head)
#> [[1]]
#>      long     lat order  hole piece id group
#> 1 -125000 2250000     1 FALSE     1  0   0.1
#> 2 -100000 2250000     2 FALSE     1  0   0.1
#> 3 -100000 2200000     3 FALSE     1  0   0.1
#> 4 -125000 2200000     4 FALSE     1  0   0.1
#> 5 -125000 2175000     5 FALSE     1  0   0.1
#> 6 -100000 2175000     6 FALSE     1  0   0.1
#> 
#> [[2]]
#>      long     lat order  hole piece id group
#> 1 -100000 2275000     1 FALSE     1  0   0.1
#> 2  -50000 2275000     2 FALSE     1  0   0.1
#> 3  -50000 2200000     3 FALSE     1  0   0.1
#> 4  -75000 2200000     4 FALSE     1  0   0.1
#> 5  -75000 2175000     5 FALSE     1  0   0.1
#> 6 -100000 2175000     6 FALSE     1  0   0.1
#> 
#> [[3]]
#>       long     lat order  hole piece id group
#> 1 -2300000 3475000     1 FALSE     1  0   0.1
#> 2 -2225000 3475000     2 FALSE     1  0   0.1
#> 3 -2225000 3400000     3 FALSE     1  0   0.1
#> 4 -2250000 3400000     4 FALSE     1  0   0.1
#> 5 -2250000 3425000     5 FALSE     1  0   0.1
#> 6 -2300000 3425000     6 FALSE     1  0   0.1
#> 
#> [[4]]
#>      long     lat order  hole piece id group
#> 1 1225000 2025000     1 FALSE     1  0   0.1
#> 2 1250000 2025000     2 FALSE     1  0   0.1
#> 3 1250000 2000000     3 FALSE     1  0   0.1
#> 4 1275000 2000000     4 FALSE     1  0   0.1
#> 5 1275000 1975000     5 FALSE     1  0   0.1
#> 6 1350000 1975000     6 FALSE     1  0   0.1
#> 
#> [[5]]
#>      long     lat order  hole piece id group
#> 1 -150000 2250000     1 FALSE     1  0   0.1
#> 2 -125000 2250000     2 FALSE     1  0   0.1
#> 3 -125000 2225000     3 FALSE     1  0   0.1
#> 4 -150000 2225000     4 FALSE     1  0   0.1
#> 5 -150000 2250000     5 FALSE     1  0   0.1
#> 6  475000 2375000     1 FALSE     1  1   1.1
```

### Map a single year/month/pole combo


```r
ggplot(res[[1]], aes(long, lat, group=group)) +
    geom_polygon(fill="steelblue") +
    theme_ice()
```

![plot of chunk seaice1](../man/figures/seaice1-1.png)

### Map all years for April only for North pole


```r
dat <- sea_ice(year = 1985:1990, month = 'Apr', pole = 'N')
df <- bind_rows(dat, .id = "x")
ggplot(df, aes(long, lat, group = group)) +
  geom_polygon(fill = "steelblue") +
  theme_ice() +
  facet_wrap(~ x)
```

![plot of chunk seaice2](../man/figures/seaice2-1.png)
---
title: "HOMR metadata"
author: "Scott Chamberlain"
date: "2020-10-06"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{HOMR metadata}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



`HOMR` (Historical Observing Metadata Repository) provides climate station metadata. It's a NOAA service.

Find out more about HOMR at https://www.ncdc.noaa.gov/homr/ and the HOMR API at https://www.ncdc.noaa.gov/homr/api

## Load rnoaa


```r
library('rnoaa')
```

## Search by station identifier

You can do this in various ways. Using the `qid` parameter (stands or qualified ID, as far as I know), you can search by suffix (e.g., `046742`), or both separated by a colon (e.g., `COOP:046742`). 

By station suffix


```r
res <- homr(qid = ':046742')
names(res)
#> [1] "20002078"
names(res[['20002078']])
#>  [1] "id"          "head"        "namez"       "identifiers" "status"     
#>  [6] "platform"    "relocations" "remarks"     "updates"     "elements"   
#> [11] "location"
res$`20002078`[1:3]
#> $id
#> [1] "20002078"
#> 
#> $head
#>                  preferredName latitude_dec longitude_dec precision
#> 1 PASO ROBLES MUNICIPAL AP, CA      35.6697     -120.6283   DDddddd
#>             por.beginDate por.endDate
#> 1 1949-10-05T00:00:00.000     Present
#> 
#> $namez
#>                         name  nameType
#> 1   PASO ROBLES MUNICIPAL AP      COOP
#> 2   PASO ROBLES MUNICIPAL AP PRINCIPAL
#> 3 PASO ROBLES MUNICIPAL ARPT       PUB
```

By both


```r
res <- homr(qid = 'COOP:046742')
names(res)
#> [1] "20002078"
names(res[['20002078']])
#>  [1] "id"          "head"        "namez"       "identifiers" "status"     
#>  [6] "platform"    "relocations" "remarks"     "updates"     "elements"   
#> [11] "location"
res$`20002078`[1:5]
#> $id
#> [1] "20002078"
#> 
#> $head
#>                  preferredName latitude_dec longitude_dec precision
#> 1 PASO ROBLES MUNICIPAL AP, CA      35.6697     -120.6283   DDddddd
#>             por.beginDate por.endDate
#> 1 1949-10-05T00:00:00.000     Present
#> 
#> $namez
#>                         name  nameType
#> 1   PASO ROBLES MUNICIPAL AP      COOP
#> 2   PASO ROBLES MUNICIPAL AP PRINCIPAL
#> 3 PASO ROBLES MUNICIPAL ARPT       PUB
#> 
#> $identifiers
#>      idType          id
#> 1     GHCND USW00093209
#> 2   GHCNMLT USW00093209
#> 3      COOP      046742
#> 4      WBAN       93209
#> 5       FAA         PRB
#> 6      ICAO        KPRB
#> 7     NWSLI         PRB
#> 8 NCDCSTNID    20002078
#> 
#> $status
#> NULL
```

## Search by station parameter

You can also search by station identifier, which is different from the `qid` above. 


```r
res <- homr(station=20002078)
names(res)
#> [1] "20002078"
names(res[['20002078']])
#>  [1] "id"          "head"        "namez"       "identifiers" "status"     
#>  [6] "platform"    "relocations" "remarks"     "updates"     "elements"   
#> [11] "location"
res$`20002078`[4:6]
#> $identifiers
#>      idType          id
#> 1     GHCND USW00093209
#> 2   GHCNMLT USW00093209
#> 3      COOP      046742
#> 4      WBAN       93209
#> 5       FAA         PRB
#> 6      ICAO        KPRB
#> 7     NWSLI         PRB
#> 8 NCDCSTNID    20002078
#> 
#> $status
#> NULL
#> 
#> $platform
#> [1] "COOP"
```

## Search by state and county

By state


```r
res <- homr(state='DE', begindate='2005-01-01', enddate='2005-02-01')
names(res)
#>  [1] "10001871" "10100161" "10100162" "10100164" "10100166" "20004155"
#>  [7] "20004158" "20004160" "20004162" "20004163" "20004167" "20004168"
#> [13] "20004171" "20004176" "20004178" "20004179" "20004180" "20004182"
#> [19] "20004184" "20004185" "30001464" "30001561" "30001831" "30075067"
```

By country


```r
res <- homr(country='GHANA', begindate='2005-01-01', enddate='2005-02-01')
library("dplyr")
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
bind_rows(lapply(res, function(x) x$location$latlon))
#>    latitude_rptd longitude_rptd latitude_dec longitude_dec latitude_dms
#> 1            6.2         -2.333          6.2        -2.333   06,12,00,N
#> 2          6.083           -.25        6.083         -0.25   06,04,59,N
#> 3          5.933          -.983        5.933        -0.983   05,55,59,N
#> 4          9.033         -2.483        9.033        -2.483   09,01,59,N
#> 5           7.75           -2.1         7.75          -2.1   07,45,00,N
#> 6         10.083         -2.508       10.083        -2.508   10,04,59,N
#> 7          9.557          -.863        9.557        -0.863   09,33,25,N
#> 8          5.783           .633        5.783         0.633   05,46,59,N
#> 9          4.867         -2.233        4.867        -2.233   04,52,01,N
#> 10           6.1           .117          6.1         0.117   06,06,00,N
#> 11         7.362         -2.329        7.362        -2.329   07,21,43,N
#> 12         5.617              0        5.617           0.0   05,37,01,N
#> 13          10.9           -1.1         10.9          -1.1   10,54,00,N
#> 14         7.817          -.033        7.817        -0.033   07,49,01,N
#> 15           6.6           .467          6.6         0.467   06,36,00,N
#> 16         5.605          -.167        5.605        -0.167   05,36,18,N
#> 17         6.715         -1.591        6.715        -1.591   06,42,54,N
#> 18         4.896         -1.775        4.896        -1.775   04,53,46,N
#> 19          5.85         -.1833         5.85       -0.1833   05,51,00,N
#> 20          5.55            -.2         5.55          -0.2   05,33,00,N
#> 21           9.4            -.9          9.4          -0.9   09,24,00,N
#> 22             5             -2          5.0          -2.0   05,00,00,N
#> 23          6.47            .33         6.47          0.33   06,28,12,N
#> 24           8.2            .57          8.2          0.57   08,12,00,N
#> 25           9.5           -.85          9.5         -0.85   09,30,00,N
#>    longitude_dms          date.beginDate            date.endDate
#> 1    002,19,59,W                 Unknown                 Present
#> 2    000,15,00,W                 Unknown                 Present
#> 3    000,58,59,W                 Unknown                 Present
#> 4    002,28,59,W                 Unknown                 Present
#> 5    002,06,00,W                 Unknown                 Present
#> 6    002,30,29,W                 Unknown                 Present
#> 7    000,51,47,W                 Unknown                 Present
#> 8    000,37,59,E                 Unknown                 Present
#> 9    002,13,59,W                 Unknown                 Present
#> 10   000,07,01,E                 Unknown                 Present
#> 11   002,19,44,W                 Unknown                 Present
#> 12   000,00,00,E                 Unknown                 Present
#> 13   001,06,00,W                 Unknown                 Present
#> 14   000,01,59,W                 Unknown                 Present
#> 15   000,28,01,E                 Unknown                 Present
#> 16   000,10,01,W                 Unknown                 Present
#> 17   001,35,28,W                 Unknown                 Present
#> 18   001,46,30,W                 Unknown                 Present
#> 19   000,11,00,W                 Unknown                 Present
#> 20   000,12,00,W                 Unknown                 Present
#> 21   000,54,00,W                 Unknown                 Present
#> 22   002,00,00,W                 Unknown                 Present
#> 23   000,19,48,E                 Unknown                 Present
#> 24   000,34,12,E                 Unknown                 Present
#> 25   000,51,00,W 1973-01-01T00:00:00.000 2008-12-31T00:00:00.000
```

By state and county


```r
res <- homr(state='NC', county='BUNCOMBE', headersOnly = TRUE)
head(bind_rows(lapply(res, "[[", "head")))
#>                     preferredName latitude_dec longitude_dec
#> 1           ASHEVILLE 5.6 NNW, NC      35.6534      -82.5709
#> 2 ASHEVILLE HENDERSONVILLE AP, NC     35.43333     -82.48333
#> 3                 WEAVERVILLE, NC         35.7     -82.56667
#> 4       BLACK MOUNTAIN 2.4 SE, NC    35.585695     -82.30557
#> 5                GARREN CREEK, NC     35.51667     -82.33333
#> 6               BILTMORE 2 SE, NC     35.56833       -82.545
#>             por.beginDate             por.endDate precision
#> 1 2007-08-27T00:00:00.000                 Present      <NA>
#> 2 1940-11-01T00:00:00.000 1960-12-31T00:00:00.000      DDMM
#> 3 1946-03-30T00:00:00.000 1992-10-01T00:00:00.000      DDMM
#> 4 2013-10-25T00:00:00.000                 Present      <NA>
#> 5 1936-09-25T00:00:00.000 1962-03-31T00:00:00.000      DDMM
#> 6 1963-10-01T00:00:00.000 2007-11-14T00:00:00.000    DDMMSS
```

## Get header information only


```r
res <- homr(headersOnly=TRUE, state='DE')
head(bind_rows(lapply(res, "[[", "head")))
#>               preferredName     latitude_dec     longitude_dec
#> 1         LEWES 1.5 SSW, DE        38.758827        -75.157875
#> 2      SELBYVILLE 7.1 E, DE        38.461155        -75.089328
#> 3       MILFORD 1.3 WSW, DE 38.9052848815918 -75.4544143676758
#> 4    MIDDLETOWN 4.1 NNW, DE 39.5044784545898 -75.7494506835938
#> 5 WILMINGTON PORTER RES, DE          39.7739          -75.5414
#> 6             BEAR 2 SW, DE          39.5917          -75.7325
#>             por.beginDate             por.endDate precision
#> 1 2014-12-08T00:00:00.000                 Present      <NA>
#> 2 2012-07-14T00:00:00.000                 Present      <NA>
#> 3 2016-03-13T00:00:00.000 2017-12-31T00:00:00.000      <NA>
#> 4 2016-03-12T00:00:00.000                 Present      <NA>
#> 5 1912-07-12T00:00:00.000                 Present   DDddddd
#> 6 2003-02-01T00:00:00.000 2013-04-02T00:00:00.000    DDMMSS
```

## Data definitions

The data returned is the same format for all, so a separate function is provided to get metadata. The function `homr_definitions()` does query the HOMR API, so does get updated metadata - i.e., it's not a static dataset stored locally. 


```r
head( homr_definitions() )
#>   defType  abbr                fullName    displayName
#> 1     ids GHCND        GHCND IDENTIFIER       GHCND ID
#> 2     ids  COOP             COOP NUMBER        COOP ID
#> 3     ids  WBAN             WBAN NUMBER        WBAN ID
#> 4     ids   FAA FAA LOCATION IDENTIFIER         FAA ID
#> 5     ids  ICAO                 ICAO ID        ICAO ID
#> 6     ids TRANS          TRANSMITTAL ID Transmittal ID
#>                                                                                                                                 description
#> 1                                                                          GLOBAL HISTORICAL CLIMATOLOGY NETWORK - DAILY (GHCND) IDENTIFIER
#> 2                                                                                   NATIONAL WEATHER SERVICE COOPERATIVE NETWORK IDENTIFIER
#> 3                                                                                                       WEATHER-BUREAU-ARMY-NAVY IDENTIFIER
#> 4                                                                                                FEDERAL AVIATION ADMINISTRATION IDENTIFIER
#> 5                                                                                      INTERNATIONAL CIVIL AVIATION ORGANIZATION IDENTIFIER
#> 6 MISCELLANEOUS IDENTIFIER THAT DOES NOT FALL INTO AN OFFICIALLY SOURCED CATEGORY AND IS NEEDED IN SUPPORT OF NCEI DATA PRODUCTS AND INGEST
#>   cssaName ghcndName
#> 1     <NA>      <NA>
#> 2     <NA>      <NA>
#> 3     <NA>      <NA>
#> 4     <NA>      <NA>
#> 5     <NA>      <NA>
#> 6     <NA>      <NA>
```
---
title: "NCDC workflow"
author: "Scott Chamberlain"
date: "2020-07-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NCDC workflow}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



This vignette is intended to demonstrate a workflow for using the NOAA NCDC data using the `ncdc*()` functions. It can be confusing to understand how to get at data you want - that's the motivation for this vignette. Other vignettes show more thorough and different examples for specific data sources.


## Load rnoaa


```r
library('rnoaa')
```

## The workflow

* Look for weather stations & get station id(s)
* Find out what type of data is available for those stations
* Search for climate data for stations (optionally specify type of data to get)

### Look for weather stations & get station id(s)



```r
ids <- ncdc_stations(locationid='FIPS:12017')$data$id[1:13]
id <- "GHCND:US1FLCT0002"
```

Just information for one station


```r
ncdc_stations(stationid = id)
```


### Find out what type of data is available for those stations

There are various ways to look for data types available. First, __data categories__:


```r
ncdc_datacats(stationid = id)
```

Another way is looking for __data sets__:


```r
ncdc_datasets(stationid = id)
```

Yet another way is looking for __data types__:


```r
ncdc_datatypes(datasetid = "GHCND", stationid = id)
```

### Search for climate data for stations (optionally specify type of data to get)

Now that you know what kinds of data categories, data sets, and data types are available for your station you can search for data with any of those as filters.

Importantly, note that you have to specify three things in a call to the `ncdc` function:

* `datasetid`
* `startdate`
* `enddate`

Here, we are specifying the `datasetid`, `stationid`, `datatypeid`, `startdate`, and `enddate`


```r
ncdc(datasetid = "GHCND", stationid = id, datatypeid = "PRCP", startdate = "2012-10-01", enddate = "2013-01-01")
```
---
title: "SWDI vignette"
author: "Scott Chamberlain"
date: "2020-10-06"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{SWDI vignette}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



`SWDI` is the Severe Weather Data Inventory. SWDI is an (quoting)

> integrated database of severe weather records for the United States. The records in SWDI come from a variety of sources in the NCDC archive. SWDI provides the ability to search through all of these data to find records covering a particular time period and geographic region, and to download the results of your search in a variety of formats. The formats currently supported are Shapefile (for GIS), KMZ (for Google Earth), CSV (comma-separated), and XML.

Data available in SWDI are:

* Storm Cells from NEXRAD (Level-III Storm Structure Product)
* Hail Signatures from NEXRAD (Level-III Hail Product)
* Mesocyclone Signatures from NEXRAD (Level-III Meso Product)
* Digital Mesocyclone Detection Algorithm from NEXRAD (Level-III MDA Product)
* Tornado Signatures from NEXRAD (Level-III TVS Product)
* Preliminary Local Storm Reports from the NOAA National Weather Service
* Lightning Strikes from Vaisala NLDN

Find out more about SWDI at https://www.ncdc.noaa.gov/ncei-severe-weather-data-inventory


## Load rnoaa


```r
library('rnoaa')
```

## Search for nx3tvs data from 5 May 2006 to 6 May 2006


```r
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506')
#> $meta
#> $meta$totalCount
#> [1] 25
#> 
#> $meta$totalTimeInSeconds
#> [1] 0.015
#> 
#> 
#> $data
#>                   ztime wsr_id cell_id cell_type range azimuth max_shear mxdv
#> 1  2006-05-05T00:05:50Z   KBMX      Q0       TVS     7     217       403  116
#> 2  2006-05-05T00:10:02Z   KBMX      Q0       TVS     5     208       421  120
#> 3  2006-05-05T00:12:34Z   KSJT      P2       TVS    49     106        17   52
#> 4  2006-05-05T00:17:31Z   KSJT      B4       TVS    40     297        25   62
#> 5  2006-05-05T00:29:13Z   KMAF      H4       TVS    53     333        34  111
#> 6  2006-05-05T00:31:25Z   KLBB      N0       TVS    51     241        24   78
#> 7  2006-05-05T00:33:25Z   KMAF      H4       TVS    52     334        46  145
#> 8  2006-05-05T00:37:37Z   KMAF      H4       TVS    50     334        34  107
#> 9  2006-05-05T00:41:51Z   KMAF      H4       TVS    51     335        29   91
#> 10 2006-05-05T00:44:33Z   KLBB      N0       TVS    46     245        35  100
#> 11 2006-05-05T00:46:03Z   KMAF      H4       TVS    49     335        41  127
#> 12 2006-05-05T00:48:55Z   KLBB      N0       TVS    44     246        44  121
#> 13 2006-05-05T00:50:16Z   KMAF      H4       TVS    49     337        33   98
#> 14 2006-05-05T00:54:29Z   KMAF      H4       TVS    47     337        42  126
#> 15 2006-05-05T00:57:42Z   KLBB      N0       TVS    41     251        46  117
#> 16 2006-05-05T00:58:41Z   KMAF      H4       TVS    46     340        29   85
#> 17 2006-05-05T01:02:04Z   KLBB      N0       TVS    39     251        42  102
#> 18 2006-05-05T01:02:53Z   KMAF      H4       TVS    46     339        35  101
#> 19 2006-05-05T01:02:53Z   KMAF      H4       TVS    50     338        27   84
#> 20 2006-05-05T01:06:26Z   KLBB      N0       TVS    36     251        31   70
#> 21 2006-05-05T01:07:06Z   KMAF      F5       TVS    45     342        44  120
#> 22 2006-05-05T01:10:48Z   KLBB      N0       TVS    36     256        37   83
#> 23 2006-05-05T01:11:18Z   KMAF      F5       TVS    45     343        39  108
#> 24 2006-05-05T01:15:30Z   KMAF      F5       TVS    44     344        30   78
#> 25 2006-05-05T01:15:30Z   KMAF      H4       TVS    49     341        26   81
#> 
#> $shape
#>                                         shape
#> 1  POINT (-86.8535716274277 33.0786326913943)
#> 2  POINT (-86.8165772540846 33.0982820681588)
#> 3  POINT (-99.5771091971025 31.1421609654838)
#> 4   POINT (-101.188161700093 31.672392833416)
#> 5  POINT (-102.664426480293 32.7306917937698)
#> 6   POINT (-102.70047613441 33.2380072329615)
#> 7    POINT (-102.6393683028 32.7226656893341)
#> 8  POINT (-102.621904684258 32.6927081076156)
#> 9   POINT (-102.614794815627 32.714139844846)
#> 10 POINT (-102.643380529494 33.3266446067682)
#> 11 POINT (-102.597961935071 32.6839260102062)
#> 12 POINT (-102.613894688178 33.3526192273658)
#> 13 POINT (-102.567153417051 32.6956373348052)
#> 14  POINT (-102.551596970251 32.664939580306)
#> 15 POINT (-102.586119971014 33.4287323151248)
#> 16 POINT (-102.499638479193 32.6644438090742)
#> 17   POINT (-102.5485490063 33.4398330734778)
#> 18  POINT (-102.51446954228 32.6597119240996)
#> 19 POINT (-102.559031583693 32.7166090376869)
#> 20 POINT (-102.492174522228 33.4564626989719)
#> 21 POINT (-102.463540844324 32.6573739036181)
#> 22 POINT (-102.510349454162 33.5066366303981)
#> 23 POINT (-102.448763863447 32.6613484943994)
#> 24   POINT (-102.42842159557 32.649061124799)
#> 25 POINT (-102.504158884526 32.7162751126854)
#> 
#> attr(,"class")
#> [1] "swdi"
```

## Use an id


```r
out <- swdi(dataset='warn', startdate='20060506', enddate='20060507', id=533623)
list(out$meta, head(out$data), head(out$shape))
#> [[1]]
#> [[1]]$totalCount
#> [1] 25
#> 
#> [[1]]$totalTimeInSeconds
#> [1] 0.009
#> 
#> 
#> [[2]]
#>            ztime_start            ztime_end     id         warningtype issuewfo
#> 1 2006-05-05T22:53:00Z 2006-05-06T00:00:00Z 397428 SEVERE THUNDERSTORM     KLCH
#> 2 2006-05-05T22:55:00Z 2006-05-06T00:00:00Z 397429 SEVERE THUNDERSTORM     KLUB
#> 3 2006-05-05T22:55:00Z 2006-05-06T00:00:00Z 397430 SEVERE THUNDERSTORM     KLUB
#> 4 2006-05-05T22:57:00Z 2006-05-06T00:00:00Z 397431 SEVERE THUNDERSTORM     KMAF
#> 5 2006-05-05T23:03:00Z 2006-05-06T00:00:00Z 397434 SEVERE THUNDERSTORM     KMAF
#> 6 2006-05-05T23:14:00Z 2006-05-06T00:15:00Z 397437 SEVERE THUNDERSTORM     KLUB
#>   messageid
#> 1    052252
#> 2    052256
#> 3    052256
#> 4    052258
#> 5    052305
#> 6    052315
#> 
#> [[3]]
#>                                                                                                                                                          shape
#> 1                                                                             POLYGON ((-93.27 30.38, -93.29 30.18, -93.02 30.18, -93.04 30.37, -93.27 30.38))
#> 2                                                                        POLYGON ((-101.93 34.74, -101.96 34.35, -101.48 34.42, -101.49 34.74, -101.93 34.74))
#> 3                POLYGON ((-100.36 33.03, -99.99 33.3, -99.99 33.39, -100.28 33.39, -100.5 33.18, -100.51 33.02, -100.45 32.97, -100.37 33.03, -100.36 33.03))
#> 4                                            POLYGON ((-102.8 30.74, -102.78 30.57, -102.15 30.61, -102.15 30.66, -101.92 30.68, -102.07 30.83, -102.8 30.74))
#> 5                                                                        POLYGON ((-103.02 32.94, -103.03 32.66, -102.21 32.53, -102.22 32.95, -103.02 32.94))
#> 6 POLYGON ((-101.6 33.32, -101.57 33.31, -101.57 33.51, -101.65 33.51, -101.66 33.5, -101.75 33.5, -101.77 33.49, -101.84 33.49, -101.84 33.32, -101.6 33.32))
```

## Get all 'plsr' within the bounding box (-91,30,-90,31)


```r
swdi(dataset='plsr', startdate='20060505', enddate='20060510', bbox=c(-91,30,-90,31))
#> $meta
#> $meta$totalCount
#> [1] 5
#> 
#> $meta$totalTimeInSeconds
#> [1] 0
#> 
#> 
#> $data
#>                  ztime     id        event magnitude            city     county
#> 1 2006-05-09T02:20:00Z 427540         HAIL         1    5 E KENTWOOD TANGIPAHOA
#> 2 2006-05-09T02:40:00Z 427536         HAIL         1    MOUNT HERMAN WASHINGTON
#> 3 2006-05-09T02:40:00Z 427537 TSTM WND DMG     -9999    MOUNT HERMAN WASHINGTON
#> 4 2006-05-09T03:00:00Z 427199         HAIL         0     FRANKLINTON WASHINGTON
#> 5 2006-05-09T03:17:00Z 427200      TORNADO     -9999 5 S FRANKLINTON WASHINGTON
#>   state          source
#> 1    LA TRAINED SPOTTER
#> 2    LA TRAINED SPOTTER
#> 3    LA TRAINED SPOTTER
#> 4    LA   AMATEUR RADIO
#> 5    LA LAW ENFORCEMENT
#> 
#> $shape
#>                  shape
#> 1 POINT (-90.43 30.93)
#> 2  POINT (-90.3 30.96)
#> 3  POINT (-90.3 30.96)
#> 4 POINT (-90.14 30.85)
#> 5 POINT (-90.14 30.78)
#> 
#> attr(,"class")
#> [1] "swdi"
```

## Get all 'nx3tvs' within the tile -102.1/32.6 (-102.15,32.55,-102.25,32.65)


```r
swdi(dataset='nx3tvs', startdate='20060506', enddate='20060507', tile=c(-102.12,32.62))
#> $meta
#> $meta$totalCount
#> [1] 5
#> 
#> $meta$totalTimeInSeconds
#> [1] 0
#> 
#> 
#> $data
#>                  ztime wsr_id cell_id cell_type range azimuth max_shear mxdv
#> 1 2006-05-06T00:41:29Z   KMAF      D9       TVS    37       6        39   85
#> 2 2006-05-06T03:56:18Z   KMAF      N4       TVS    39       3        30   73
#> 3 2006-05-06T03:56:18Z   KMAF      N4       TVS    42       4        20   52
#> 4 2006-05-06T04:00:30Z   KMAF      N4       TVS    38       5        35   86
#> 5 2006-05-06T04:04:44Z   KMAF      N4       TVS    41       8        24   62
#> 
#> $shape
#>                                        shape
#> 1 POINT (-102.112726356403 32.5574494581267)
#> 2  POINT (-102.14873079873 32.5933553250156)
#> 3 POINT (-102.131167022161 32.6426287452898)
#> 4 POINT (-102.123671677514 32.5751241756203)
#> 5 POINT (-102.076389686189 32.6209390786829)
#> 
#> attr(,"class")
#> [1] "swdi"
```

## Counts

Notes:

* stat='count' will only return metadata, nothing in the data or shape slots
* stat='tilesum:...' returns counts in the data slot for each date for that tile, and shape data
* Get number of 'nx3tvs' within 15 miles of latitude = 32.7 and longitude = -102.0

Get daily count nx3tvs features on .1 degree grid centered at `latitude = 32.7` and `longitude = -102.0`


```r
swdi(dataset='nx3tvs', startdate='20060505', enddate='20090516', stat='tilesum:-102.0,32.7')
#> $meta
#> $meta$totalCount
#> [1] 5
#> 
#> $meta$totalTimeInSeconds
#> [1] 0
#> 
#> 
#> $data
#>          day centerlat centerlon fcount
#> 1 2007-03-29      32.7      -102      2
#> 2 2007-09-07      32.7      -102      1
#> 3 2008-05-27      32.7      -102      4
#> 4 2008-06-20      32.7      -102      2
#> 5 2009-04-11      32.7      -102      1
#> 
#> $shape
#>                                                                                   shape
#> 1 POLYGON ((-102.05 32.65, -102.05 32.75, -101.95 32.75, -101.95 32.65, -102.05 32.65))
#> 2 POLYGON ((-102.05 32.65, -102.05 32.75, -101.95 32.75, -101.95 32.65, -102.05 32.65))
#> 3 POLYGON ((-102.05 32.65, -102.05 32.75, -101.95 32.75, -101.95 32.65, -102.05 32.65))
#> 4 POLYGON ((-102.05 32.65, -102.05 32.75, -101.95 32.75, -101.95 32.65, -102.05 32.65))
#> 5 POLYGON ((-102.05 32.65, -102.05 32.75, -101.95 32.75, -101.95 32.65, -102.05 32.65))
#> 
#> attr(,"class")
#> [1] "swdi"
```

## Get data in different formats

### CSV format


```r
head(swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506', format='csv')$data)
#>                  ztime wsr_id cell_id cell_type range azimuth max_shear mxdv
#> 1 2006-05-05T00:05:50Z   KBMX      Q0       TVS     7     217       403  116
#> 2 2006-05-05T00:10:02Z   KBMX      Q0       TVS     5     208       421  120
#> 3 2006-05-05T00:12:34Z   KSJT      P2       TVS    49     106        17   52
#> 4 2006-05-05T00:17:31Z   KSJT      B4       TVS    40     297        25   62
#> 5 2006-05-05T00:29:13Z   KMAF      H4       TVS    53     333        34  111
#> 6 2006-05-05T00:31:25Z   KLBB      N0       TVS    51     241        24   78
#>      lat      lon
#> 1 33.079  -86.854
#> 2 33.098  -86.817
#> 3 31.142  -99.577
#> 4 31.672 -101.188
#> 5 32.731 -102.664
#> 6 33.238 -102.700
```

### SHP format


```r
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506', format='shp', filepath='myfile')
#> shp file downloaded to myfile.zip
```

### KMZ format


```r
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506', format='kmz', radius=15, filepath='myfile.kmz')
#> kmz file downloaded to myfile.kmz
```



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seaice.r
\name{theme_ice}
\alias{theme_ice}
\title{ggplot2 map theme}
\usage{
theme_ice()
}
\description{
ggplot2 map theme
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_utils.r
\name{meteo_coverage}
\alias{meteo_coverage}
\title{Determine the "coverage" for a station data frame}
\usage{
meteo_coverage(
  meteo_df,
  obs_start_date = NULL,
  obs_end_date = NULL,
  verbose = FALSE
)
}
\arguments{
\item{meteo_df}{a \emph{meteo} \code{data.frame}}

\item{obs_start_date}{specify either or both (obs_start_date, obs_end_date)
to constrain coverate tests. These should be \code{Date} objects.}

\item{obs_end_date}{specify either or both (obs_start_date, obs_end_date)
to constrain coverate tests. These should be \code{Date} objects.}

\item{verbose}{if \code{TRUE} will display the coverage summary along
with returning the coverage data.frame}
}
\value{
a \code{list} containing 2 \code{data.frame}s named 'summary' and 'detail'.
The 'summary' \code{data.frame} contains columns: \preformatted{
$ id         (chr)
$ start_date (time)
$ end_date   (time)
$ total_obs  (int)
}
with additional fields (and their coverage percent) depending on which
weather variables were queried and available for the weather station. The
\code{data.frame} named 'detail' contains the same columns as the \code{meteo_df} input
data, but expands the rows to contain \code{NA}s for days without data.
}
\description{
Call this function after pulling down observations for a set of stations
to retrieve the "coverage" (i.e. how complete each field is). If either
or both \code{obs_start_date} or \code{obs_end_date} are specified,
the coverage test will be limited to that date range.
}
\details{
There is an \code{autoplot} method for the output of this function.
}
\examples{
\dontrun{

monitors <- c("ASN00095063", "ASN00024025", "ASN00040112", "ASN00041023",
             "ASN00009998", "ASN00066078", "ASN00003069", "ASN00090162",
             "ASN00040126", "ASN00058161")
obs <- meteo_pull_monitors(monitors)
obs_covr <- meteo_coverage(obs)

if (interactive()) {
  library("ggplot2")
  autoplot(obs_covr)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ghcnd_stations.R
\name{ghcnd_stations}
\alias{ghcnd_stations}
\title{Get information on the GHCND weather stations}
\usage{
ghcnd_stations(refresh = FALSE, ...)
}
\arguments{
\item{refresh}{(logical) If \code{TRUE} force re-download of data.
Default: \code{FALSE}}

\item{...}{In the case of \code{ghcnd()} additional curl options to pass
through to \link[crul:HttpClient]{crul::HttpClient}. In the case of \code{ghcnd_read}
further options passed on to \code{read.csv}}
}
\value{
This function returns a tibble (dataframe) with a weather
station on each row with the following columns:
\itemize{
\item \code{id}: The weather station's ID number. The first two letters
denote the country (using FIPS country codes).
\item \code{latitude}: The station's latitude, in decimal degrees.
Southern latitudes will be negative.
\item \code{longitude}: The station's longitude, in decimal degrees.
Western longitudes will be negative.
\item \code{elevation}: The station's elevation, in meters.
\item \code{name}: The station's name.
\item \code{gsn_flag}: "GSN" if the monitor belongs to the GCOS Surface
Network (GSN). Otherwise either blank or missing.
\item \code{wmo_id}: If the station has a WMO number, this column gives
that number. Otherwise either blank or missing.
\item \code{element}: A weather variable recorded at some point during
that station's history. See the link below in "References" for
definitions of the abbreviations used for this variable.
\item \code{first_year}: The first year of data available at that station
for that weather element.
\item \code{last_year}: The last year of data available at that station
for that weather element.
}

If a weather station has data on more than one weather variable,
it will be represented in multiple rows of this output dataframe.
}
\description{
This function returns an object with a dataframe with meta-information
about all available GHCND weather stations.
}
\note{
Since this function is pulling a large dataset by ftp, it may take
a while to run.
}
\examples{
\dontrun{
# Get stations, ghcnd-stations and ghcnd-inventory merged
(stations <- ghcnd_stations())

library(dplyr)
# filter by state
stations \%>\% filter(state == "IL")
stations \%>\% filter(state == "OR")
# those without state values
stations \%>\% filter(state == "")
# filter by element
stations \%>\% filter(element == "PRCP")
# filter by id prefix
stations \%>\% filter(grepl("^AF", id))
stations \%>\% filter(grepl("^AFM", id))
# filter by station long name
stations \%>\% filter(name == "CALLATHARRA")
}
}
\references{
For more documentation on the returned dataset, see
http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsw.R
\name{bsw}
\alias{bsw}
\title{Blended sea winds (BSW)}
\usage{
bsw(date = NULL, uv_stress = "uv", resolution = "6hrly", ...)
}
\arguments{
\item{date}{(date/character) date, in the form YYYY-MM-DD if resolution
is 6hrly or daily, or in the form YYYY-MM if resolution is monthly.
For resolution=clm can be left NULL. If given, must be in the
range 1987-07-09 to today-1 (yesterday)}

\item{uv_stress}{(character) one of uv or stresss, not sure what these
mean exactly yet. Default: uv}

\item{resolution}{(character) temporal resolution. one of 6hrly,
clm, daily, or monthly. See Details.}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
an object of class \code{ncdf4}
}
\description{
The Blended Sea Winds dataset contains globally gridded, high-resolution
ocean surface vector winds and wind stresses on a global 0.25° grid, and
multiple time resolutions of six-hourly, daily, monthly, and 11-year
(1995–2005) climatological monthlies.
}
\details{
Products are available from July 9th, 1987 - present.

Uses \code{ncdf4} under the hood to read NetCDF files
}
\note{
See \link{bsw_cache} for managing cached files

We only handle the netcdf files for now, we're avoiding the ieee
files, see https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/ieee.html
}
\section{Citing NOAA and BSW data}{

Message from NOAA: "We also ask you to acknowledge us in your use of the
data to help us justify continued service. This may be done by  including
text such as: The wind data are acquired from NOAA's National Climatic
Data Center, via their website We would also  appreciate receiving a
copy of the relevant publication."
}

\section{Temporal resolution}{

\itemize{
\item 6hrly: 6-hourly, 4 global snapshots (u,v) at UTC 00, 06, 12 and 18Z
\item clm: climatological monthlies; also provided is the scalar
mean (u,v,w)
\item daily: averages of the 6hrly time points, thus with a center time
09Z; also provided is the scalar mean, (u,v,w)
\item monthly: averages of daily data; also provided is the scalar
mean (u,v,w)
}
}

\examples{
\dontrun{
# 6hrly data
## uv
x <- bsw(date = "2017-10-01")
## stress
y <- bsw(date = "2011-08-01", uv_stress = "stress")

# daily
z <- bsw(date = "2017-10-01", resolution = "daily")

# monthly
w <- bsw(date = "2011-08", resolution = "monthly")

# clm
# x <- bsw(resolution = "clm")
}
}
\references{
https://www.ncdc.noaa.gov/data-access/marineocean-data/blended-global/blended-sea-winds
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_miss.R
\name{vis_miss}
\alias{vis_miss}
\title{Visualize missingness in a dataframe}
\usage{
vis_miss(x, cluster = FALSE, sort_miss = FALSE)
}
\arguments{
\item{x}{a data.frame}

\item{cluster}{logical \code{TRUE}/\code{FALSE}. \code{TRUE} specifies that you want to use
hierarchical clustering (mcquitty method) to arrange rows according to
missingness. \code{FALSE} specifies that you want to leave it as is.}

\item{sort_miss}{logical \code{TRUE}/\code{FALSE}. \code{TRUE} arranges the columns in
order of missingness.}
}
\description{
Gives you an at-a-glance ggplot of the missingness inside a dataframe,
colouring cells according to missingness, where black indicates a present
cell and grey indicates a missing cell. As it returns a \code{ggplot} object,
it is very easy to customize and change labels, and so on.
}
\details{
\code{vis_miss} visualises a data.frame to display missingness. This is
taken from the visdat package, currently only available on github:
https://github.com/tierneyn/visdat
}
\examples{
\dontrun{
  monitors <- c("ASN00003003", "ASM00094299")
  weather_df <- meteo_pull_monitors(monitors)
  vis_miss(weather_df)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{rnoaa-defunct}
\alias{rnoaa-defunct}
\title{Defunct functions in rnoaa}
\description{
\itemize{
\item \code{noaa}: Function name changed, prefixed with ncdc now
\item \code{noaa_datacats}: Function name changed, prefixed with ncdc now
\item \code{noaa_datasets}: Function name changed, prefixed with ncdc now
\item \code{noaa_datatypes}: Function name changed, prefixed with ncdc now
\item \code{noaa_locs}: Function name changed, prefixed with ncdc now
\item \code{noaa_locs_cats}: Function name changed, prefixed with ncdc now
\item \code{noaa_stations}: Function name changed, prefixed with ncdc now
\item \code{noaa_plot}: Function name changed, prefixed with ncdc now
\item \code{noaa_combine}: Function name changed, prefixed with ncdc now
\item \code{noaa_seaice}: Function name changed to seaice
\item \code{erddap_data}: See package rerddap
\item \code{erddap_clear_cache}: See package rerddap
\item \code{erddap_datasets}: Moved to package rerddap
\item \code{erddap_grid}: Moved to package rerddap
\item \code{erddap_info}: Moved to \code{rerddap::info()}
\item \code{erddap_search}: Moved to \code{rerddap::ed_search}
\item \code{erddap_table}: Moved to \code{rerddap::tabledap}
\item \code{ncdc_leg_variables}: Removed. See \verb{NCDC Legacy} below
\item \code{ncdc_leg_sites}: Removed. See \verb{NCDC Legacy} below
\item \code{ncdc_leg_site_info}: Removed. See \verb{NCDC Legacy} below
\item \code{ncdc_leg_data}: Removed. See \verb{NCDC Legacy} below
\item \code{seaice}: Replaced with \code{\link[=sea_ice]{sea_ice()}}
\item \code{lcd_cleanup}: No longer available. See \code{lcd} docs
\item \code{ghcnd_clear_cache}: No longer available. See \link{rnoaa_caching}
\item \code{storm_shp}: Function defunct.
\item \code{storm_shp_read}: Function defunct.
\item \code{storm_data}: Function defunct.
\item \code{storm_meta}: Function defunct.
}
}
\details{
The functions for working with GEFS ensemble forecast data (prefixed with
"gefs") are defunct, but may come back to rnoaa later:
\itemize{
\item \code{\link[=gefs]{gefs()}}
\item \code{\link[=gefs_dimension_values]{gefs_dimension_values()}}
\item \code{\link[=gefs_dimensions]{gefs_dimensions()}}
\item \code{\link[=gefs_ensembles]{gefs_ensembles()}}
\item \code{\link[=gefs_latitudes]{gefs_latitudes()}}
\item \code{\link[=gefs_longitudes]{gefs_longitudes()}}
\item \code{\link[=gefs_times]{gefs_times()}}
\item \code{\link[=gefs_variables]{gefs_variables()}}
}
}
\section{NCDC Legacy}{

The NCDC legacy API is too unreliable and slow. Use the newer NCDC API via
the functions \code{\link[=ncdc]{ncdc()}}, \code{\link[=ncdc_datacats]{ncdc_datacats()}}, \code{\link[=ncdc_datasets]{ncdc_datasets()}},
\code{\link[=ncdc_datatypes]{ncdc_datatypes()}}, \code{\link[=ncdc_locs]{ncdc_locs()}}, \code{\link[=ncdc_locs_cats]{ncdc_locs_cats()}}, \code{\link[=ncdc_stations]{ncdc_stations()}},
\code{\link[=ncdc_plot]{ncdc_plot()}}, and \code{\link[=ncdc_combine]{ncdc_combine()}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_distance.R
\name{meteo_nearby_stations}
\alias{meteo_nearby_stations}
\title{Find weather monitors near locations}
\usage{
meteo_nearby_stations(
  lat_lon_df,
  lat_colname = "latitude",
  lon_colname = "longitude",
  station_data = ghcnd_stations(),
  var = "all",
  year_min = NULL,
  year_max = NULL,
  radius = NULL,
  limit = NULL
)
}
\arguments{
\item{lat_lon_df}{A dataframe that contains the latitude, longitude, and
a unique identifier for each location (\code{id}). For an example of the
proper format for this dataframe, see the examples below. Latitude and
longitude must both be in units of decimal degrees. Southern latitudes
and Western longitudes should be given as negative values. A tibble
is accepted, but is coerced to a data.frame internally before any usage.}

\item{lat_colname}{A character string giving the name of the latitude column
in the \code{lat_lon_df} dataframe.}

\item{lon_colname}{A character string giving the name of the longitude column
in the \code{lat_lon_df} dataframe.}

\item{station_data}{The output of \code{\link[=ghcnd_stations]{ghcnd_stations()}}, which is
a current list of weather stations available through NOAA for the GHCND
dataset. The format of this is a dataframe
with one row per weather station. Latitude and longitude for the station
locations should be in columns with the names "latitude" and "longitude",
consistent with the output from \code{\link[=ghcnd_stations]{ghcnd_stations()}}. To save time,
run the \code{ghcnd_stations} call and save the output to an object,
rather than rerunning the default every time (see the examples in
\code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}}).}

\item{var}{A character vector specifying either \code{"all"} (pull all
available weather parameters for the site) or the weather parameters to
keep in the final data (e.g., \code{c("TMAX", "TMIN")} to only keep
maximum and minimum temperature). Example choices for this argument
include:
\itemize{
\item \code{PRCP}: Precipitation, in tenths of millimeters
\item \code{TAVG}: Average temperature, in tenths of degrees Celsius
\item \code{TMAX}: Maximum temperature, in tenths of degrees Celsius
\item \code{TMIN}: Minimum temperature, in tenths of degrees Celsius
}

A full list of possible weather variables is available in NOAA's README
file for the GHCND data
(https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt).
Most weather stations will only have a small subset of all the possible
weather variables, so the data generated by this function may not include
all of the variables the user specifies through this argument.}

\item{year_min}{A numeric value giving the earliest year from which you
ultimately want weather data (e.g., 2013, if you only are interested in
data from 2013 and later).}

\item{year_max}{A numeric value giving the latest year from which you
ultimately want weather data.}

\item{radius}{A numeric vector giving the radius (in kilometers) within which
to search for monitors near a location.}

\item{limit}{An integer giving the maximum number of monitors to include for
each location. The \code{x} closest monitors will be kept. Default is NULL
(pull everything available, within the radius if the radius is specified).}
}
\value{
A list containing dataframes with the sets of unique weather stations
within the search radius for each location. Site IDs for the weather
stations given in this dataframe can be used in conjunction with other
functions in the \pkg{rnoaa} package to pull weather data for the
station. The dataframe for each location includes:
\itemize{
\item \code{id}: The weather station ID, which can be used in other
functions to pull weather data from the station;
\item \code{name}: The weather station name;
\item \code{latitude}: The station's latitude, in decimal degrees.
Southern latitudes will be negative;
\item \code{longitude}: The station's longitude, in decimal degrees.
Western longitudes will be negative;
\item \code{distance}: The station's distance, in kilometers, from the
location.
}
}
\description{
This function inputs a dataframe with latitudes and longitudes of locations
and creates a dataframe with monitors within a certain radius of those
locations. The function can also be used, with the \code{limit} argument, to
pull a certain number of the closest weather monitors to each location.
The weather monitor IDs in the output dataframe can be used with other
\pkg{rnoaa} functions to pull data from all available weather stations near
a location (e.g., \code{\link[=meteo_pull_monitors]{meteo_pull_monitors()}}).
}
\details{
Great circle distance is used to determine whether a weather monitor is
within the required radius.
}
\note{
By default, this function will pull the full station list from NOAA
to use to identify nearby locations. If you will be creating lists of
monitors nearby several stations, you can save some time by using the
\code{\link[=ghcnd_stations]{ghcnd_stations()}} function separately to create an object
with all stations and then use the argument \code{station_data} in
this function to reference that object, rather than using this function's
defaults (see examples).
}
\examples{
\dontrun{

station_data <- ghcnd_stations() # Takes a while to run

lat_lon_df <- data.frame(id = c("sydney", "brisbane"),
                         latitude = c(-33.8675, -27.4710),
                         longitude = c(151.2070, 153.0234))
nearby_stations <-  meteo_nearby_stations(lat_lon_df = lat_lon_df,
                    station_data = station_data, radius = 10)

miami <- data.frame(id = "miami", latitude = 25.7617, longitude = -80.1918)

# Get all stations within 50 kilometers
meteo_nearby_stations(lat_lon_df = miami, station_data = station_data,
                      radius = 50, var = c("PRCP", "TMAX"),
                      year_min = 1992, year_max = 1992)
# Get the closest 10 monitors
meteo_nearby_stations(lat_lon_df = miami, station_data = station_data,
                      limit = 10, var = c("PRCP", "TMAX"),
                      year_min = 1992, year_max = 1992)
}
}
\seealso{
The weather monitor IDs generated by this function can be used in
other functions in the \pkg{rnoaa} package, like
\code{\link[=meteo_pull_monitors]{meteo_pull_monitors()}} and \code{\link[=meteo_tidy_ghcnd]{meteo_tidy_ghcnd()}}, to
pull weather data from weather monitors near a location.
}
\author{
Alex Simmons \email{a2.simmons@qut.edu.au},
Brooke Anderson \email{brooke.anderson@colostate.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_datasets.r
\name{ncdc_datasets}
\alias{ncdc_datasets}
\title{Search NOAA datasets}
\usage{
ncdc_datasets(
  datasetid = NULL,
  datatypeid = NULL,
  stationid = NULL,
  locationid = NULL,
  startdate = NULL,
  enddate = NULL,
  sortfield = NULL,
  sortorder = NULL,
  limit = 25,
  offset = NULL,
  token = NULL,
  ...
)
}
\arguments{
\item{datasetid}{(optional) Accepts a single valid dataset id. Data returned
will be from the dataset specified.}

\item{datatypeid}{Accepts a valid data type id or a vector or list of data
type ids. (optional)}

\item{stationid}{Accepts a valid station id or a vector or list of station
ids}

\item{locationid}{Accepts a valid location id or a vector or list of
location ids (optional)}

\item{startdate}{(optional) Accepts valid ISO formated date (yyyy-mm-dd) or date time
(YYYY-MM-DDThh:mm:ss). Data returned will have data after the specified date. The
date range must be less than 1 year.}

\item{enddate}{(optional) Accepts valid ISO formated date (yyyy-mm-dd) or date time
(YYYY-MM-DDThh:mm:ss). Data returned will have data before the specified date. The
date range must be less than 1 year.}

\item{sortfield}{The field to sort results by. Supports id, name, mindate,
maxdate, and datacoverage fields (optional)}

\item{sortorder}{Which order to sort by, asc or desc. Defaults to
asc (optional)}

\item{limit}{Defaults to 25, limits the number of results in the response.
Maximum is 1000 (optional)}

\item{offset}{Defaults to 0, used to offset the resultlist (optional)}

\item{token}{This must be a valid token token supplied to you by NCDC's
Climate Data Online access token generator. (required) See
\strong{Authentication} section below for more details.}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}
(optional)}
}
\value{
A data.frame for all datasets, or a list of length two, each with
a data.frame.
}
\description{
From the NOAA API docs: All of our data are in datasets. To retrieve
any data from us, you must know what dataset it is in.
}
\section{Authentication}{

Get an API key (aka, token) at https://www.ncdc.noaa.gov/cdo-web/token
You can pass your token in as an argument or store it one of two places:

\itemize{
\item your .Rprofile file with the entry
\code{options(noaakey = "your-noaa-token")}
\item your .Renviron file with the entry
\code{NOAA_KEY=your-noaa-token}
}

See \code{\link{Startup}} for information on how to create/find your
.Rrofile and .Renviron files
}

\examples{
\dontrun{
# Get a table of all datasets
ncdc_datasets()

# Get details from a particular dataset
ncdc_datasets(datasetid='ANNUAL')

# Get datasets with Temperature at the time of observation (TOBS) data type
ncdc_datasets(datatypeid='TOBS')
## two datatypeid's
ncdc_datasets(datatypeid=c('TOBS', "ACMH"))

# Get datasets with data for a series of the same parameter arg, in this case
# stationid's
ncdc_datasets(stationid='COOP:310090')
ncdc_datasets(stationid=c('COOP:310090','COOP:310184','COOP:310212'))

# Multiple datatypeid's
ncdc_datasets(datatypeid=c('ACMC','ACMH','ACSC'))
ncdc_datasets(datasetid='ANNUAL', datatypeid=c('ACMC','ACMH','ACSC'))
ncdc_datasets(datasetid='GSOY', datatypeid=c('ACMC','ACMH','ACSC'))

# Multiple locationid's
ncdc_datasets(locationid="FIPS:30091")
ncdc_datasets(locationid=c("FIPS:30103", "FIPS:30091"))
}
}
\references{
https://www.ncdc.noaa.gov/cdo-web/webservices/v2
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc_stations}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_cache.r
\name{meteo_show_cache}
\alias{meteo_show_cache}
\title{Show the \emph{meteo} cache directory}
\usage{
meteo_show_cache()
}
\description{
Displays the full path to the \code{meteo} cache directory
}
\seealso{
Other meteo: 
\code{\link{meteo_clear_cache}()}
}
\concept{meteo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seaice.r
\name{sea_ice}
\alias{sea_ice}
\title{Get sea ice data.}
\usage{
sea_ice(year = NULL, month = NULL, pole = NULL, format = "shp", ...)
}
\arguments{
\item{year}{(numeric) a year}

\item{month}{(character) a month, as character abbrevation of a month}

\item{pole}{(character) one of S (south) or N (north)}

\item{format}{(character) one of shp (default), geotiff-extent (for geotiff
extent data), or geotiff-conc (for geotiff concentration data)}

\item{...}{Further arguments passed on to \code{rgdal::readshpfile()} if
\code{format="shp"} or \code{raster::raster()} if not}
}
\value{
data.frame if \code{format="shp"} (a fortified sp object);
\code{raster::raster()} if not
}
\description{
Get sea ice data.
}
\examples{
\dontrun{
if (requireNamespace("raster")) {

## one year, one moth, one pole
sea_ice(year = 1990, month = "Apr", pole = "N")
sea_ice(year = 1990, month = "Apr", pole = "N", format = "geotiff-extent")
sea_ice(year = 1990, month = "Apr", pole = "N", format = "geotiff-conc")

## one year, one month, many poles
sea_ice(year = 1990, month = "Apr")

## one year, many months, many poles
sea_ice(year = 1990, month = c("Apr", "Jun", "Oct"))

## many years, one month, one pole
sea_ice(year = 1990:1992, month = "Sep", pole = "N")

# get geotiff instead of shp data. 
x <- sea_ice(year = 1990, month = "Apr", format = "geotiff-extent")
y <- sea_ice(year = 1990, month = "Apr", format = "geotiff-conc")
}

}
}
\references{
See the "User Guide" pdf at https://nsidc.org/data/g02135
}
\seealso{
\code{\link[=sea_ice_tabular]{sea_ice_tabular()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs_variables}
\alias{gefs_variables}
\title{This function is defunct.}
\usage{
gefs_variables(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs_ensembles}
\alias{gefs_ensembles}
\title{This function is defunct.}
\usage{
gefs_ensembles(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs}
\alias{gefs}
\title{This function is defunct.}
\usage{
gefs(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs_longitudes}
\alias{gefs_longitudes}
\title{This function is defunct.}
\usage{
gefs_longitudes(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homr_definitions.R
\name{homr_definitions}
\alias{homr_definitions}
\title{Historical Observing Metadata Repository (HOMR) station metadata -
definitions}
\usage{
homr_definitions(...)
}
\arguments{
\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}
optional}
}
\description{
Historical Observing Metadata Repository (HOMR) station metadata -
definitions
}
\examples{
\dontrun{
head( homr_definitions() )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_stations.r
\name{ncdc_stations}
\alias{ncdc_stations}
\title{Get metadata about NOAA NCDC stations.}
\usage{
ncdc_stations(
  stationid = NULL,
  datasetid = NULL,
  datatypeid = NULL,
  locationid = NULL,
  startdate = NULL,
  enddate = NULL,
  sortfield = NULL,
  sortorder = NULL,
  limit = 25,
  offset = NULL,
  datacategoryid = NULL,
  extent = NULL,
  token = NULL,
  dataset = NULL,
  station = NULL,
  location = NULL,
  locationtype = NULL,
  page = NULL,
  ...
)
}
\arguments{
\item{stationid}{A single valid station id, with datasetid namespace,
e.g., GHCND:USW00014895}

\item{datasetid}{(optional) Accepts a valid dataset id or a vector or
list of them. Data returned will be from the dataset specified.}

\item{datatypeid}{Accepts a valid data type id or a vector or list of data
type ids. (optional)}

\item{locationid}{Accepts a valid location id or a vector or list of
location ids (optional)}

\item{startdate}{(optional) Accepts valid ISO formated date (yyyy-mm-dd) or date time
(YYYY-MM-DDThh:mm:ss). Data returned will have data after the specified date. The
date range must be less than 1 year.}

\item{enddate}{(optional) Accepts valid ISO formated date (yyyy-mm-dd) or date time
(YYYY-MM-DDThh:mm:ss). Data returned will have data before the specified date. The
date range must be less than 1 year.}

\item{sortfield}{The field to sort results by. Supports id, name, mindate,
maxdate, and datacoverage fields (optional)}

\item{sortorder}{Which order to sort by, asc or desc. Defaults to
asc (optional)}

\item{limit}{Defaults to 25, limits the number of results in the response.
Maximum is 1000 (optional)}

\item{offset}{Defaults to 0, used to offset the resultlist (optional)}

\item{datacategoryid}{(character, optional) Accepts a valid data category id or a
vector or list of data category ids.}

\item{extent}{(numeric, optional) The geographical extent for which you want to
search. Give four values that defines a bounding box, lat and long for the
southwest corner, then lat and long for the northeast corner. For example:
\code{c(minlat, minlong, maxlat, maxlong)}.}

\item{token}{This must be a valid token token supplied to you by NCDC's
Climate Data Online access token generator. (required) See
\strong{Authentication} section below for more details.}

\item{dataset}{THIS IS A DEPRECATED ARGUMENT. See datasetid.}

\item{station}{THIS IS A DEPRECATED ARGUMENT. See stationid.}

\item{location}{THIS IS A DEPRECATED ARGUMENT. See locationid.}

\item{locationtype}{THIS IS A DEPRECATED ARGUMENT. There is no equivalent argument in v2
of the NOAA API.}

\item{page}{THIS IS A DEPRECATED ARGUMENT. There is no equivalent argument in v2
of the NOAA API.}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}
(optional)}
}
\value{
A list of metadata.
}
\description{
From the NOAA NCDC API docs: Stations are where the data comes from
(for most datasets) and can be considered the smallest granual of location
data. If you know what station you want, you can quickly get all manner of
data from it
}
\section{Authentication}{

Get an API key (aka, token) at https://www.ncdc.noaa.gov/cdo-web/token
You can pass your token in as an argument or store it one of two places:

\itemize{
\item your .Rprofile file with the entry
\code{options(noaakey = "your-noaa-token")}
\item your .Renviron file with the entry
\code{NOAA_KEY=your-noaa-token}
}

See \code{\link{Startup}} for information on how to create/find your
.Rrofile and .Renviron files
}

\examples{
\dontrun{
# Get metadata on all stations
ncdc_stations()
ncdc_stations(limit=5)

# Get metadata on a single station
ncdc_stations(stationid='COOP:010008')

# For many stations use lapply or similar
lapply(c("COOP:010008", "COOP:010063", "COOP:010116"), function(z) {
  ncdc_stations(
   startdate = "2013-01-01",
   enddate = "2014-11-01",
   stationid = z)
}$data)

# Displays all stations within GHCN-Daily (100 Stations per page limit)
ncdc_stations(datasetid = 'GHCND')
ncdc_stations(datasetid = 'ANNUAL')
ncdc_stations(datasetid = 'GSOY')

# Station
ncdc_stations(datasetid='NORMAL_DLY', stationid='GHCND:USW00014895')

# datatypeid
ncdc_stations(datatypeid="ANN-HTDD-NORMAL")
ncdc_stations(datatypeid=c("ANN-HTDD-NORMAL", "ACSC"))

# locationid
ncdc_stations(locationid="CITY:AG000001")
ncdc_stations(locationid="FIPS:30091")
ncdc_stations(locationid=c("FIPS:30103", "FIPS:30091"))

# datacategoryid
ncdc_stations(datacategoryid="ANNPRCP")
ncdc_stations(datacategoryid="AUAGR")
ncdc_stations(datacategoryid=c("ANNPRCP", "AUAGR"))

# Displays all stations within GHCN-Daily (Displaying page 10 of the results)
ncdc_stations(datasetid='GHCND')

# Specify datasetid and locationid
ncdc_stations(datasetid='GHCND', locationid='FIPS:12017')

# Specify datasetid, locationid, and station
ncdc_stations(datasetid='GHCND', locationid='FIPS:12017', stationid='GHCND:USC00084289')

# Specify datasetid, locationidtype, locationid, and station
ncdc_stations(datasetid='GHCND', locationid='FIPS:12017', stationid='GHCND:USC00084289')

# Displays list of stations within the specified county
ncdc_stations(datasetid='GHCND', locationid='FIPS:12017')

# Displays list of Hourly Precipitation locationids between 01/01/1990 and 12/31/1990
ncdc_stations(datasetid='PRECIP_HLY', startdate='19900101', enddate='19901231')

# Search for stations by spatial extent
## Search using a bounding box, w/ lat/long of the SW corner, then of NE corner
ncdc_stations(extent=c(47.5204,-122.2047,47.6139,-122.1065))
}
}
\references{
https://www.ncdc.noaa.gov/cdo-web/webservices/v2
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homr.R
\name{homr}
\alias{homr}
\title{Historical Observing Metadata Repository (HOMR) station metadata}
\usage{
homr(
  qid = NULL,
  qidMod = NULL,
  station = NULL,
  state = NULL,
  county = NULL,
  country = NULL,
  name = NULL,
  nameMod = NULL,
  platform = NULL,
  date = NULL,
  begindate = NULL,
  enddate = NULL,
  headersOnly = FALSE,
  phrData = NULL,
  combine = FALSE,
  ...
)
}
\arguments{
\item{qid}{One of COOP, FAA, GHCND, ICAO, NCDCSTNID, NWSLI, TRANS, WBAN, or
WMO, or any of those plus \code{a-z0-9}, or just \code{a-z0-9}. (qid = qualified ID)}

\item{qidMod}{(character) One of: is, starts, ends, contains. Specifies
how the ID portion of the qid parameter should be applied within the search.
If a qid is passed but the qidMod parameter is not used, qidMod is
assumed to be IS.}

\item{station}{(character) A station id.}

\item{state}{(character) A two-letter state abbreviation. Two-letter code
for US states, Canadian provinces, and other Island areas.}

\item{county}{(character) A two letter county code. US county names, best
used with a state identifier.}

\item{country}{(character) A two letter country code. See here for a list
of valid country names.}

\item{name}{(character) One of \verb{0-9A-Z+}. Searches on any type of
name we have for the station.}

\item{nameMod}{(character) \code{is|starts|ends|contains}. Specifies how the
name parameter should be applied within the search. If a name is passed but
the nameMod parameter is not used, nameMod is assumed to be IS.}

\item{platform}{(character) (aka network) \verb{ASOS|USCRN|USHCN|NEXRAD|AL USRCRN|USRCRN|COOP}. Limit the search to stations of a certain
platform/network type.}

\item{date}{(character) \code{YYYY-MM-DD|all} Limits values to only those that
occurred on a specific date. Alternatively, date=all will return all values
for matched stations. If this field is omitted, the search will return only
the most recent values for each field.}

\item{begindate, enddate}{\code{YYYY-MM-DD}. Limits values to only those that
occurred within a date range.}

\item{headersOnly}{(logical) Returns only minimal information for each
station found (NCDC Station ID, Preferred Name, Station Begin Date, and
Station End Date), but is much quicker than a full query. If you are
performing a search that returns a large number of stations and intend to
choose only one from that list to examine in detail, headersOnly may give
you enough information to find the NCDC Station ID for the station that
you actually want.}

\item{phrData}{(logical) The HOMR web service now includes PHR
(element-level) data when available, in an elements section. Because of
how this data is structured, it can substantially increase the size of any
result which includes it. If you don't need this data you can omit it
by including phrData=false. If the parameter is not set, it will default
to phrData=true.}

\item{combine}{(logical) Combine station metadata or not.}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET} (optional)}
}
\value{
A list, with elements named by the station ids.
}
\description{
Historical Observing Metadata Repository (HOMR) station metadata
}
\details{
Since the definitions for variables are always the same, we don't
include the ability to get description data in this function. Use
\code{\link[=homr_definitions]{homr_definitions()}} to get descriptions information.
}
\examples{
\dontrun{
homr(qid = 'COOP:046742')
homr(qid = ':046742')
homr(qidMod='starts', qid='COOP:0467')
homr(headersOnly=TRUE, state='DE')
homr(headersOnly=TRUE, country='GHANA')
homr(headersOnly=TRUE, state='NC', county='BUNCOMBE')
homr(name='CLAYTON')
res <- homr(state='NC', county='BUNCOMBE', combine=TRUE)
res$id
res$head
res$updates
homr(nameMod='starts', name='CLAY')
homr(headersOnly=TRUE, platform='ASOS')
homr(qid='COOP:046742', date='2011-01-01')
homr(qid='COOP:046742', begindate='2005-01-01', enddate='2011-01-01')
homr(state='DE', headersOnly=TRUE)
homr(station=20002078)
homr(station=20002078, date='all', phrData=FALSE)

# Optionally pass in curl options
homr(headersOnly=TRUE, state='NC', county='BUNCOMBE', verbose = TRUE)
}
}
\references{
https://www.ncdc.noaa.gov/homr/api
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_plot.r
\name{ncdc_plot}
\alias{ncdc_plot}
\title{Plot NOAA climate data.}
\usage{
ncdc_plot(..., breaks = NULL, dateformat = "\%d/\%m/\%y")
}
\arguments{
\item{...}{Input noaa object or objects.}

\item{breaks}{Regularly spaced date breaks for x-axis. See examples for
usage. See \link{date_breaks}. Default: \code{NULL} (uses ggplot2 default break
sformatting)}

\item{dateformat}{Date format using standard POSIX specification for labels
on x-axis. See \code{\link[=date_format]{date_format()}}}
}
\value{
ggplot2 plot
}
\description{
Plot NOAA climate data.
}
\details{
This function accepts directly output from the \code{\link[=ncdc]{ncdc()}} function,
not other functions.

This is a simple wrapper function around some ggplot2 code. There is indeed
a lot you can modify in your plots, so this function just does some basic
stuff. Look at the internals for what the function does.
}
\examples{
\dontrun{
# Search for data first, then plot
out <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP',
   startdate = '2010-05-01', enddate = '2010-10-31', limit=500)
ncdc_plot(out)
ncdc_plot(out, breaks="14 days")
ncdc_plot(out, breaks="1 month", dateformat="\%d/\%m")
ncdc_plot(out, breaks="1 month", dateformat="\%d/\%m")

# Combine many calls to ncdc function
out1 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP',
   startdate = '2010-03-01', enddate = '2010-05-31', limit=500)
out2 <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP',
   startdate = '2010-09-01', enddate = '2010-10-31', limit=500)
df <- ncdc_combine(out1, out2)
ncdc_plot(df)
## or pass in each element separately
ncdc_plot(out1, out2, breaks="45 days")
}
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_stations}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_locs_cats.r
\name{ncdc_locs_cats}
\alias{ncdc_locs_cats}
\title{Get metadata about NOAA location categories.}
\usage{
ncdc_locs_cats(
  datasetid = NULL,
  locationcategoryid = NULL,
  startdate = NULL,
  enddate = NULL,
  sortfield = NULL,
  sortorder = NULL,
  limit = 25,
  offset = NULL,
  token = NULL,
  ...
)
}
\arguments{
\item{datasetid}{A valid dataset id or a vector or list of dataset id's. Data returned will be from the
dataset specified, see datasets() (required)}

\item{locationcategoryid}{A valid location id or a vector or list of location category ids}

\item{startdate}{A valid ISO formatted date (yyyy-mm-dd). Data returned will have
data after the specified date. Paramater can be use independently of enddate (optional)}

\item{enddate}{Accepts valid ISO formatted date (yyyy-mm-dd). Data returned will have data
before the specified date. Paramater can be use independently of startdate (optional)}

\item{sortfield}{The field to sort results by. Supports id, name, mindate, maxdate, and
datacoverage fields (optional)}

\item{sortorder}{Which order to sort by, asc or desc. Defaults to asc (optional)}

\item{limit}{Defaults to 25, limits the number of results in the response. Maximum is
1000 (optional)}

\item{offset}{Defaults to 0, used to offset the resultlist (optional)}

\item{token}{This must be a valid token token supplied to you by NCDC's
Climate Data Online access token generator. (required) See
\strong{Authentication} section below for more details.}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}}
}
\value{
A list containing metadata and the data, or a single data.frame.
}
\description{
Location categories are groupings of similar locations.
}
\details{
Locations can be a specific latitude/longitude point such as a station,
or a label representing a bounding area such as a city.
}
\section{Authentication}{

Get an API key (aka, token) at https://www.ncdc.noaa.gov/cdo-web/token
You can pass your token in as an argument or store it one of two places:

\itemize{
\item your .Rprofile file with the entry
\code{options(noaakey = "your-noaa-token")}
\item your .Renviron file with the entry
\code{NOAA_KEY=your-noaa-token}
}

See \code{\link{Startup}} for information on how to create/find your
.Rrofile and .Renviron files
}

\examples{
\dontrun{
# All location categories, first 25 results
ncdc_locs_cats()

# Find locations with category id of CLIM_REG
ncdc_locs_cats(locationcategoryid='CLIM_REG')

# Displays available location categories within GHCN-Daily dataset
ncdc_locs_cats(datasetid='GHCND')
ncdc_locs_cats(datasetid='GSOY')
ncdc_locs_cats(datasetid='ANNUAL')

# multiple datasetid's
ncdc_locs_cats(datasetid=c('GHCND', 'GSOM'))

# Displays available location categories from start date 1970-01-01
ncdc_locs_cats(startdate='1970-01-01')
}
}
\references{
https://www.ncdc.noaa.gov/cdo-web/webservices/v2
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc_stations}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_locs.r
\name{ncdc_locs}
\alias{ncdc_locs}
\title{Get metadata about NOAA NCDC locations.}
\usage{
ncdc_locs(
  datasetid = NULL,
  locationid = NULL,
  locationcategoryid = NULL,
  startdate = NULL,
  enddate = NULL,
  sortfield = NULL,
  sortorder = NULL,
  limit = 25,
  offset = NULL,
  token = NULL,
  ...
)
}
\arguments{
\item{datasetid}{A valid dataset id or a vector or list of dataset id's. Data returned will be from the
dataset specified, see datasets() (required)}

\item{locationid}{A valid location id or a vector or list of location ids.}

\item{locationcategoryid}{A valid location id or a vector or list of location category ids}

\item{startdate}{A valid ISO formatted date (yyyy-mm-dd). Data returned will have
data after the specified date. Paramater can be use independently of enddate (optional)}

\item{enddate}{Accepts valid ISO formatted date (yyyy-mm-dd). Data returned will have data
before the specified date. Paramater can be use independently of startdate (optional)}

\item{sortfield}{The field to sort results by. Supports id, name, mindate, maxdate, and
datacoverage fields (optional)}

\item{sortorder}{Which order to sort by, asc or desc. Defaults to asc (optional)}

\item{limit}{Defaults to 25, limits the number of results in the response. Maximum is
1000 (optional)}

\item{offset}{Defaults to 0, used to offset the resultlist (optional)}

\item{token}{This must be a valid token token supplied to you by NCDC's
Climate Data Online access token generator. (required) See
\strong{Authentication} section below for more details.}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}}
}
\value{
A list containing metadata and the data, or a single data.frame.
}
\description{
From the NOAA NCDC API docs: Locations can be a specific latitude/longitude
point such as a station, or a label representing a bounding area such as
a city.
}
\section{Authentication}{

Get an API key (aka, token) at https://www.ncdc.noaa.gov/cdo-web/token
You can pass your token in as an argument or store it one of two places:

\itemize{
\item your .Rprofile file with the entry
\code{options(noaakey = "your-noaa-token")}
\item your .Renviron file with the entry
\code{NOAA_KEY=your-noaa-token}
}

See \code{\link{Startup}} for information on how to create/find your
.Rrofile and .Renviron files
}

\examples{
\dontrun{
# All locations, first 25 results
ncdc_locs()

# Fetch more information about location id FIPS:37
ncdc_locs(locationid='FIPS:37')

# Fetch available locations for the GHCND (Daily Summaries) dataset
ncdc_locs(datasetid='GHCND')
ncdc_locs(datasetid=c('GHCND', 'ANNUAL'))
ncdc_locs(datasetid=c('GSOY', 'ANNUAL'))
ncdc_locs(datasetid=c('GHCND', 'GSOM'))

# Fetch all U.S. States
ncdc_locs(locationcategoryid='ST', limit=52)

# Many locationcategoryid's
## this apparently works, but returns nothing often with multiple
## locationcategoryid's
ncdc_locs(locationcategoryid=c('ST', 'ZIP'))

# Fetch list of city locations in descending order
ncdc_locs(locationcategoryid='CITY', sortfield='name', sortorder='desc')
}
}
\references{
https://www.ncdc.noaa.gov/cdo-web/webservices/v2
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc_stations}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd_stations.R
\name{isd_stations}
\alias{isd_stations}
\title{Get NOAA ISD/ISH station data from NOAA FTP server.}
\usage{
isd_stations(refresh = FALSE)
}
\arguments{
\item{refresh}{(logical) Download station data from NOAA ftp server again.
Default: \code{FALSE}}
}
\value{
a tibble (data.frame) with the columns:
\itemize{
\item usaf - USAF number, character
\item wban - WBAN number, character
\item station_name - station name, character
\item ctry - Country, if given, character
\item state - State, if given, character
\item icao - ICAO number, if given, character
\item lat - Latitude, if given, numeric
\item lon - Longitude, if given, numeric
\item elev_m - Elevation, if given, numeric
\item begin - Begin date of data coverage, of form YYYYMMDD, numeric
\item end - End date of data coverage, of form YYYYMMDD, numeric
}
}
\description{
Get NOAA ISD/ISH station data from NOAA FTP server.
}
\details{
The data table is cached, but you can force download of data from
NOAA by setting \code{refresh=TRUE}
}
\note{
See \link{isd_cache} for managing cached files
}
\examples{
\dontrun{
# Get station table
(stations <- isd_stations())

## plot stations
### remove incomplete cases, those at 0,0
df <- stations[complete.cases(stations$lat, stations$lon), ]
df <- df[df$lat != 0, ]
### make plot
library("leaflet")
leaflet(data = df) \%>\%
  addTiles() \%>\%
  addCircles()
}
}
\references{
https://ftp.ncdc.noaa.gov/pub/data/noaa/
}
\seealso{
Other isd: 
\code{\link{isd_read}()},
\code{\link{isd_stations_search}()},
\code{\link{isd}()}
}
\concept{isd}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_distance.R
\name{meteo_process_geographic_data}
\alias{meteo_process_geographic_data}
\title{Calculate the distances between a location and all available stations}
\usage{
meteo_process_geographic_data(station_data, lat, long, units = "deg")
}
\arguments{
\item{station_data}{The output of \code{\link[=ghcnd_stations]{ghcnd_stations()}}, which is
a current list of weather stations available through NOAA for the GHCND
dataset. The format of this is a dataframe
with one row per weather station. Latitude and longitude for the station
locations should be in columns with the names "latitude" and "longitude",
consistent with the output from \code{\link[=ghcnd_stations]{ghcnd_stations()}}. To save time,
run the \code{ghcnd_stations} call and save the output to an object,
rather than rerunning the default every time (see the examples in
\code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}}).}

\item{lat}{Latitude of the location. Southern latitudes should be given
as negative values.}

\item{long}{Longitude of the location. Western longitudes should be given as
negative values.}

\item{units}{Units of the latitude and longitude values. Possible values
are:
\itemize{
\item \code{deg}: Degrees (default);
\item \code{rad}: Radians.
}}
}
\value{
The \code{station_data} dataframe that is input, but with a
\code{distance} column added that gives the distance to the location
(in kilometers), and re-ordered by distance between each station and
the location (closest weather stations first).
}
\description{
This function takes a single location and a dataset of available weather stations
and calculates the distance between the location and each of the stations,
using the great circle method. A new column is added to the dataset of
available weather stations giving the distance between each station and
the input location. The station dataset is then sorted from closest to
furthest distance to the location and returned as the function output.
}
\examples{
\dontrun{
station_data <- ghcnd_stations()
meteo_process_geographic_data(station_data, lat=-33, long=151)
}
}
\author{
Alex Simmons \email{a2.simmons@qut.edu.au},
Brooke Anderson \email{brooke.anderson@colostate.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs_dimensions}
\alias{gefs_dimensions}
\title{This function is defunct.}
\usage{
gefs_dimensions(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seaice_tabular.R
\name{sea_ice_tabular}
\alias{sea_ice_tabular}
\title{Sea ice tabular data}
\usage{
sea_ice_tabular(...)
}
\arguments{
\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET} - beware
that curl options are passed to each http request, for each of
24 requests.}
}
\value{
A data.frame with columns:
\itemize{
\item year (integer)
\item mo (integer)
\item data.type (character)
\item region (character)
\item extent (numeric)
\item area (numeric)
}
}
\description{
Collects \code{.csv} files from NOAA, and binds them together into
a single data.frame. Data across years, with extent and area of ice.
}
\details{
An example file, for January, North pole:
\verb{https://sidads.colorado.edu/DATASETS/NOAA/G02135/north/monthly/data/N_01_extent_v3.0.csv}

a value in any cell of -9999 indicates missing data
}
\examples{
\dontrun{
df <- sea_ice_tabular()
df
}
}
\seealso{
\code{\link[=sea_ice]{sea_ice()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnoaa-package.r
\docType{data}
\name{fipscodes}
\alias{fipscodes}
\title{FIPS codes for US states.}
\format{
A data frame with 3142 rows and 5 variables
}
\description{
A dataset containing the FIPS codes for 51 US states
and territories. The variables are as follows:
}
\details{
\itemize{
\item state. US state name.
\item county. County name.
\item fips_state. Numeric value, from 1 to 51.
\item fips_county. Numeric value, from 1 to 840.
\item fips. Numeric value, from 1001 to 56045.
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lcd.R
\name{lcd_columns}
\alias{lcd_columns}
\title{Specify Variable Types in Local Climatological Data from NOAA}
\usage{
lcd_columns(
  STATION = "character",
  DATE = "POSIXct",
  LATITUDE = "numeric",
  LONGITUDE = "numeric",
  ELEVATION = "numeric",
  NAME = "character",
  REPORT_TYPE = "character",
  SOURCE = "character",
  HourlyAltimeterSetting = "character",
  HourlyDewPointTemperature = "character",
  HourlyDryBulbTemperature = "character",
  HourlyPrecipitation = "character",
  HourlyPresentWeatherType = "character",
  HourlyPressureChange = "character",
  HourlyPressureTendency = "integer",
  HourlyRelativeHumidity = "character",
  HourlySkyConditions = "character",
  HourlySeaLevelPressure = "character",
  HourlyStationPressure = "character",
  HourlyVisibility = "character",
  HourlyWetBulbTemperature = "character",
  HourlyWindDirection = "character",
  HourlyWindGustSpeed = "character",
  HourlyWindSpeed = "character",
  Sunrise = "numeric",
  Sunset = "numeric",
  DailyAverageDewPointTemperature = "character",
  DailyAverageDryBulbTemperature = "character",
  DailyAverageRelativeHumidity = "character",
  DailyAverageSeaLevelPressure = "character",
  DailyAverageStationPressure = "character",
  DailyAverageWetBulbTemperature = "character",
  DailyAverageWindSpeed = "character",
  DailyCoolingDegreeDays = "numeric",
  DailyDepartureFromNormalAverageTemperature = "numeric",
  DailyHeatingDegreeDays = "numeric",
  DailyMaximumDryBulbTemperature = "numeric",
  DailyMinimumDryBulbTemperature = "numeric",
  DailyPeakWindDirection = "numeric",
  DailyPeakWindSpeed = "numeric",
  DailyPrecipitation = "character",
  DailySnowDepth = "character",
  DailySnowfall = "character",
  DailySustainedWindDirection = "numeric",
  DailySustainedWindSpeed = "numeric",
  DailyWeather = "character",
  MonthlyAverageRH = "numeric",
  MonthlyDaysWithGT001Precip = "numeric",
  MonthlyDaysWithGT010Precip = "numeric",
  MonthlyDaysWithGT32Temp = "numeric",
  MonthlyDaysWithGT90Temp = "numeric",
  MonthlyDaysWithLT0Temp = "numeric",
  MonthlyDaysWithLT32Temp = "numeric",
  MonthlyDepartureFromNormalAverageTemperature = "numeric",
  MonthlyDepartureFromNormalCoolingDegreeDays = "numeric",
  MonthlyDepartureFromNormalHeatingDegreeDays = "numeric",
  MonthlyDepartureFromNormalMaximumTemperature = "numeric",
  MonthlyDepartureFromNormalMinimumTemperature = "numeric",
  MonthlyDepartureFromNormalPrecipitation = "numeric",
  MonthlyDewpointTemperature = "numeric",
  MonthlyGreatestPrecip = "numeric",
  MonthlyGreatestPrecipDate = "character",
  MonthlyGreatestSnowDepth = "numeric",
  MonthlyGreatestSnowDepthDate = "character",
  MonthlyGreatestSnowfall = "numeric",
  MonthlyGreatestSnowfallDate = "character",
  MonthlyMaxSeaLevelPressureValue = "numeric",
  MonthlyMaxSeaLevelPressureValueDate = "character",
  MonthlyMaxSeaLevelPressureValueTime = "character",
  MonthlyMaximumTemperature = "numeric",
  MonthlyMeanTemperature = "numeric",
  MonthlyMinSeaLevelPressureValue = "numeric",
  MonthlyMinSeaLevelPressureValueDate = "character",
  MonthlyMinSeaLevelPressureValueTime = "character",
  MonthlyMinimumTemperature = "numeric",
  MonthlySeaLevelPressure = "numeric",
  MonthlyStationPressure = "numeric",
  MonthlyTotalLiquidPrecipitation = "numeric",
  MonthlyTotalSnowfall = "numeric",
  MonthlyWetBulb = "numeric",
  AWND = "numeric",
  CDSD = "numeric",
  CLDD = "numeric",
  DSNW = "numeric",
  HDSD = "numeric",
  HTDD = "numeric",
  NormalsCoolingDegreeDay = "numeric",
  NormalsHeatingDegreeDay = "numeric",
  ShortDurationEndDate005 = "character",
  ShortDurationEndDate010 = "character",
  ShortDurationEndDate015 = "character",
  ShortDurationEndDate020 = "character",
  ShortDurationEndDate030 = "character",
  ShortDurationEndDate045 = "character",
  ShortDurationEndDate060 = "character",
  ShortDurationEndDate080 = "character",
  ShortDurationEndDate100 = "character",
  ShortDurationEndDate120 = "character",
  ShortDurationEndDate150 = "character",
  ShortDurationEndDate180 = "character",
  ShortDurationPrecipitationValue005 = "numeric",
  ShortDurationPrecipitationValue010 = "numeric",
  ShortDurationPrecipitationValue015 = "numeric",
  ShortDurationPrecipitationValue020 = "numeric",
  ShortDurationPrecipitationValue030 = "numeric",
  ShortDurationPrecipitationValue045 = "numeric",
  ShortDurationPrecipitationValue060 = "numeric",
  ShortDurationPrecipitationValue080 = "numeric",
  ShortDurationPrecipitationValue100 = "numeric",
  ShortDurationPrecipitationValue120 = "numeric",
  ShortDurationPrecipitationValue150 = "numeric",
  ShortDurationPrecipitationValue180 = "numeric",
  REM = "character",
  BackupDirection = "character",
  BackupDistance = "character",
  BackupDistanceUnit = "character",
  BackupElements = "character",
  BackupElevation = "character",
  BackupEquipment = "character",
  BackupLatitude = "character",
  BackupLongitude = "character",
  BackupName = "character",
  WindEquipmentChangeDate = "character"
)
}
\arguments{
\item{STATION}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{LATITUDE}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{LONGITUDE}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ELEVATION}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{NAME}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{REPORT_TYPE}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{SOURCE}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyAltimeterSetting}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyDewPointTemperature}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyDryBulbTemperature}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyPrecipitation}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyPresentWeatherType}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyPressureChange}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyPressureTendency}{(character) string indicating variable or column type that is returned, default is "integer". optional}

\item{HourlyRelativeHumidity}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlySkyConditions}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlySeaLevelPressure}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyStationPressure}{(character)string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyVisibility}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyWetBulbTemperature}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyWindDirection}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyWindGustSpeed}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{HourlyWindSpeed}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{Sunrise}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{Sunset}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyAverageDewPointTemperature}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailyAverageDryBulbTemperature}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailyAverageRelativeHumidity}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailyAverageSeaLevelPressure}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailyAverageStationPressure}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailyAverageWetBulbTemperature}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailyAverageWindSpeed}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailyCoolingDegreeDays}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyDepartureFromNormalAverageTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyHeatingDegreeDays}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyMaximumDryBulbTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyMinimumDryBulbTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyPeakWindDirection}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyPeakWindSpeed}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyPrecipitation}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailySnowDepth}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailySnowfall}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{DailySustainedWindDirection}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailySustainedWindSpeed}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DailyWeather}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyAverageRH}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDaysWithGT001Precip}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDaysWithGT010Precip}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDaysWithGT32Temp}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDaysWithGT90Temp}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDaysWithLT0Temp}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDaysWithLT32Temp}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDepartureFromNormalAverageTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDepartureFromNormalCoolingDegreeDays}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDepartureFromNormalHeatingDegreeDays}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDepartureFromNormalMaximumTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDepartureFromNormalMinimumTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDepartureFromNormalPrecipitation}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyDewpointTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyGreatestPrecip}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyGreatestPrecipDate}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyGreatestSnowDepth}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyGreatestSnowDepthDate}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyGreatestSnowfall}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyGreatestSnowfallDate}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyMaxSeaLevelPressureValue}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyMaxSeaLevelPressureValueDate}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyMaxSeaLevelPressureValueTime}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyMaximumTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyMeanTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyMinSeaLevelPressureValue}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyMinSeaLevelPressureValueDate}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyMinSeaLevelPressureValueTime}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{MonthlyMinimumTemperature}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlySeaLevelPressure}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyStationPressure}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyTotalLiquidPrecipitation}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyTotalSnowfall}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{MonthlyWetBulb}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{AWND}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{CDSD}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{CLDD}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{DSNW}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{HDSD}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{HTDD}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{NormalsCoolingDegreeDay}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{NormalsHeatingDegreeDay}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationEndDate005}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate010}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate015}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate020}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate030}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate045}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate060}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate080}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate100}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate120}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate150}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationEndDate180}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{ShortDurationPrecipitationValue005}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue010}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue015}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue020}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue030}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue045}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue060}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue080}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue100}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue120}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue150}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{ShortDurationPrecipitationValue180}{(character) string indicating variable or column type that is returned, default is "numeric". optional}

\item{REM}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupDirection}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupDistance}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupDistanceUnit}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupElements}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupElevation}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupEquipment}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupLatitude}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupLongitude}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{BackupName}{(character) string indicating variable or column type that is returned, default is "character". optional}

\item{WindEquipmentChangeDate}{(character) string indicating variable or column type that is returned, default is "character". optional}
}
\value{
a vector indicating column classes types
}
\description{
Use this function to specify what variable types will be returned by
\link{lcd}. The function returns a named vector with specified column
classes. The defaults are specified in the argument descriptions below.
}
\note{
if the column type is not compatible, \link{lcd} will return a
dataframe with the most appropriate column type and a message indicating
the column was not changed to the specified type.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_ghcnd.R
\name{meteo_tidy_ghcnd}
\alias{meteo_tidy_ghcnd}
\title{Create a tidy GHCND dataset from a single monitor}
\usage{
meteo_tidy_ghcnd(
  stationid,
  keep_flags = FALSE,
  var = "all",
  date_min = NULL,
  date_max = NULL
)
}
\arguments{
\item{stationid}{(character) A character vector giving the identification of
the weather stations for which the user would like to pull data. To get a full
and current list of stations, the user can use the \code{\link[=ghcnd_stations]{ghcnd_stations()}}
function. To identify stations within a certain radius of a location, the
user can use the \code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}} function.}

\item{keep_flags}{TRUE / FALSE for whether the user would like to keep all
the flags for each weather variable. The default is to not keep the
flags (FALSE). See the note below for more information on these flags.}

\item{var}{A character vector specifying either \code{"all"} (pull all
available weather parameters for the site) or the weather parameters to
keep in the final data (e.g., \code{c("TMAX", "TMIN")} to only keep
maximum and minimum temperature). Example choices for this argument
include:
\itemize{
\item \code{PRCP}: Precipitation, in tenths of millimeters
\item \code{TAVG}: Average temperature, in tenths of degrees Celsius
\item \code{TMAX}: Maximum temperature, in tenths of degrees Celsius
\item \code{TMIN}: Minimum temperature, in tenths of degrees Celsius
}

A full list of possible weather variables is available in NOAA's README
file for the GHCND data
(https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt).
Most weather stations will only have a small subset of all the possible
weather variables, so the data generated by this function may not include
all of the variables the user specifies through this argument.}

\item{date_min}{A character string giving the earliest
date of the daily weather time series that the user would
like in the final output. This character string should be formatted as
"yyyy-mm-dd". If not specified, the default is to keep all daily data for
the queried weather site from the earliest available date.}

\item{date_max}{A character string giving the latest
date of the daily weather time series that the user would
like in the final output. This character string should be formatted as
"yyyy-mm-dd". If not specified, the default is to keep all daily data for
the queried weather site through the most current available date.}
}
\value{
A data frame of daily weather data for a single weather monitor,
converted to a tidy format. All weather variables may not exist for all
weather stations. Examples of variables returned are:
\itemize{
\item \code{id}: Character string with the weather station site id
\item \code{date}: Date of the observation
\item \code{prcp}: Precipitation, in tenths of mm
\item \code{tavg}: Average temperature, in degrees Celsius
\item \code{tmax}: Maximum temperature, in degrees Celsius
\item \code{tmin}: Minimum temperature, in degrees Celsius
\item \code{awnd}: Average daily wind speed, in meters / second
\item \code{wsfg}: Peak gust wind speed, in meters / second
}

There are other possible weather variables in the Global Historical
Climatology Network; see
http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt for a full
list. The variables \code{prcp}, \code{tmax}, \code{tmin}, and \code{tavg}
have all been converted from tenths of their metric to the metric (e.g.,
from tenths of degrees Celsius to degrees Celsius). All other variables
are in the units specified in the linked file.
}
\description{
This function inputs an object created by \code{\link{ghcnd}} and cleans up
the data into a tidy form.
}
\note{
The weather flags, which are kept by specifying
\code{keep_flags = TRUE} are:
\itemize{
\item \verb{*_mflag}: Measurement flag, which gives some information on how
the observation was measured.
\item \verb{*_qflag}: Quality flag, which gives quality information on the
measurement, like if it failed to pass certain quality checks.
\item \verb{*_sflag}: Source flag. This gives some information on the
weather collection system (e.g., U.S. Cooperative Summary of the Day,
Australian Bureau of Meteorology) the weather observation comes from.
}

More information on the interpretation of these flags can be found in the
README file for the NCDC's Daily Global Historical Climatology Network's
data at http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
}
\examples{
\dontrun{
# One station in Australia is ASM00094275
meteo_tidy_ghcnd(stationid = "ASN00003003")
meteo_tidy_ghcnd(stationid = "ASN00003003", var = "tavg")
meteo_tidy_ghcnd(stationid = "ASN00003003", date_min = "1989-01-01")
}

}
\seealso{
\code{\link[=meteo_pull_monitors]{meteo_pull_monitors()}}
}
\author{
Brooke Anderson \email{brooke.anderson@colostate.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_datacats.r
\name{ncdc_datacats}
\alias{ncdc_datacats}
\title{Get possible data categories for a particular datasetid, locationid,
stationid, etc.}
\usage{
ncdc_datacats(
  datasetid = NULL,
  datacategoryid = NULL,
  stationid = NULL,
  locationid = NULL,
  startdate = NULL,
  enddate = NULL,
  sortfield = NULL,
  sortorder = NULL,
  limit = 25,
  offset = NULL,
  token = NULL,
  ...
)
}
\arguments{
\item{datasetid}{Accepts a valid dataset id or a vector or list of dataset id's. Data returned
will be from the dataset specified, see datasets() (required)}

\item{datacategoryid}{A valid data category id. Data types returned will be associated
with the data category(ies) specified}

\item{stationid}{Accepts a valid station id or a vector or list of station ids (optional)}

\item{locationid}{Accepts a valid location id or a vector or list of location id's. (optional)}

\item{startdate}{Accepts valid ISO formated date (yyyy-mm-dd). Data returned will have
data after the specified date. Paramater can be use independently of enddate (optional)}

\item{enddate}{Accepts valid ISO formated date (yyyy-mm-dd). Data returned will have data
before the specified date. Paramater can be use independently of startdate (optional)}

\item{sortfield}{The field to sort results by. Supports id, name, mindate, maxdate, and
datacoverage fields (optional)}

\item{sortorder}{Which order to sort by, asc or desc. Defaults to asc (optional)}

\item{limit}{Defaults to 25, limits the number of results in the response. Maximum is
1000 (optional)}

\item{offset}{Defaults to 0, used to offset the resultlist (optional)}

\item{token}{This must be a valid token token supplied to you by NCDC's
Climate Data Online access token generator. (required) See
\strong{Authentication} section below for more details.}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}}
}
\value{
A \code{data.frame} for all datasets, or a list of length two,
each with a data.frame.
}
\description{
Data Categories represent groupings of data types.
}
\details{
Note that calls with both startdate and enddate don't seem to
work, though specifying one or the other mostly works.
}
\section{Authentication}{

Get an API key (aka, token) at https://www.ncdc.noaa.gov/cdo-web/token
You can pass your token in as an argument or store it one of two places:

\itemize{
\item your .Rprofile file with the entry
\code{options(noaakey = "your-noaa-token")}
\item your .Renviron file with the entry
\code{NOAA_KEY=your-noaa-token}
}

See \code{\link{Startup}} for information on how to create/find your
.Rrofile and .Renviron files
}

\examples{
\dontrun{
## Limit to 10 results
ncdc_datacats(limit=10)

## by datasetid
ncdc_datacats(datasetid="ANNUAL")
ncdc_datacats(datasetid=c("ANNUAL", "PRECIP_HLY"))

## Single data category
ncdc_datacats(datacategoryid="ANNAGR")

## Fetch data categories for a given set of locations
ncdc_datacats(locationid='CITY:US390029')
ncdc_datacats(locationid=c('CITY:US390029', 'FIPS:37'))

## Data categories for a given date
ncdc_datacats(startdate = '2013-10-01')

# Get data categories with data for a series of the same parameter arg, in this case
# stationid's
ncdc_datacats(stationid='COOP:310090')
ncdc_datacats(stationid=c('COOP:310090','COOP:310184','COOP:310212'))

## Curl debugging
ncdc_datacats(limit=10, verbose = TRUE)
}
}
\references{
https://www.ncdc.noaa.gov/cdo-web/webservices/v2
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc_stations}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{check_response}
\alias{check_response}
\title{Check response from NOAA, including status codes, server error messages, mime-type, etc.}
\usage{
check_response(x)
}
\description{
Check response from NOAA, including status codes, server error messages, mime-type, etc.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd_read.R
\name{isd_read}
\alias{isd_read}
\title{Read NOAA ISD/ISH local file}
\usage{
isd_read(
  path,
  additional = TRUE,
  parallel = FALSE,
  cores = getOption("cl.cores", 2),
  progress = FALSE
)
}
\arguments{
\item{path}{(character) path to the file. required.}

\item{additional}{(logical) include additional and remarks data sections
in output. Default: \code{TRUE}. Passed on to
\code{\link[isdparser:isd_parse]{isdparser::isd_parse()}}}

\item{parallel}{(logical) do processing in parallel. Default: \code{FALSE}}

\item{cores}{(integer) number of cores to use: Default: 2. We look in
your option "cl.cores", but use default value if not found.}

\item{progress}{(logical) print progress - ignored if \code{parallel=TRUE}.
The default is \code{FALSE} because printing progress adds a small bit of
time, so if processing time is important, then keep as \code{FALSE}}
}
\value{
A tibble (data.frame)
}
\description{
Read NOAA ISD/ISH local file
}
\details{
\code{isd_read} - read a \code{.gz} file as downloaded
from NOAA's website
}
\examples{
\dontrun{
file <- system.file("examples", "011490-99999-1986.gz", package = "rnoaa")
isd_read(file)
isd_read(file, additional = FALSE)
}
}
\references{
https://ftp.ncdc.noaa.gov/pub/data/noaa/
}
\seealso{
\code{\link[=isd]{isd()}}, \code{\link[=isd_stations]{isd_stations()}}, \code{\link[=isd_stations_search]{isd_stations_search()}}

Other isd: 
\code{\link{isd_stations_search}()},
\code{\link{isd_stations}()},
\code{\link{isd}()}
}
\concept{isd}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs_latitudes}
\alias{gefs_latitudes}
\title{This function is defunct.}
\usage{
gefs_latitudes(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{ncdc_theme}
\alias{ncdc_theme}
\title{Theme for plotting NOAA data}
\usage{
ncdc_theme()
}
\description{
Theme for plotting NOAA data
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{check_response_swdi}
\alias{check_response_swdi}
\title{Check response from NOAA SWDI service, including status codes, server error messages,
mime-type, etc.}
\usage{
check_response_swdi(x, format)
}
\description{
Check response from NOAA SWDI service, including status codes, server error messages,
mime-type, etc.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd_stations_search.R
\name{isd_stations_search}
\alias{isd_stations_search}
\title{Search for NOAA ISD/ISH station data from NOAA FTP server.}
\usage{
isd_stations_search(lat = NULL, lon = NULL, radius = NULL, bbox = NULL)
}
\arguments{
\item{lat}{(numeric) Latitude, in decimal degree}

\item{lon}{(numeric) Latitude, in decimal degree}

\item{radius}{(numeric) Radius (in km) to search from the lat,lon
coordinates}

\item{bbox}{(numeric) Bounding box, of the form: min-longitude,
min-latitude, max-longitude, max-latitude}
}
\value{
a data.frame with the columns:
\itemize{
\item usaf - USAF number, character
\item wban - WBAN number, character
\item station_name - station name, character
\item ctry - Country, if given, character
\item state - State, if given, character
\item icao - ICAO number, if given, character
\item lat - Latitude, if given, numeric
\item lon - Longitude, if given, numeric
\item elev_m - Elevation, if given, numeric
\item begin - Begin date of data coverage, of form YYYYMMDD, numeric
\item end - End date of data coverage, of form YYYYMMDD, numeric
\item distance - distance (km) (only present if using lat/lon/radius
parameter combination)
}
}
\description{
Search for NOAA ISD/ISH station data from NOAA FTP server.
}
\details{
We internally call \code{\link[=isd_stations]{isd_stations()}} to get the data.frame
of ISD stations, which is quite fast as long as it's not the first time
called since we cache the table. Before searching, we clean up the
data.frame, removing stations with no lat/long coordinates, those with
impossible lat/long coordinates, and those at 0,0.

When lat/lon/radius input we use \code{\link[=meteo_distance]{meteo_distance()}} to search
for stations, while when bbox is input, we simply use
\code{\link[dplyr:filter]{dplyr::filter()}}
}
\examples{
\dontrun{
## lat, long, radius
isd_stations_search(lat = 38.4, lon = -123, radius = 250)

x <- isd_stations_search(lat = 60, lon = 18, radius = 200)

if (requireNamespace("leaflet")) {
  library("leaflet")
  leaflet() \%>\%
    addTiles() \%>\%
    addCircles(lng = x$lon,
               lat = x$lat,
               popup = x$station_name) \%>\%
    clearBounds()
}

## bounding box
bbox <- c(-125.0, 38.4, -121.8, 40.9)
isd_stations_search(bbox = bbox)
}
}
\references{
https://ftp.ncdc.noaa.gov/pub/data/noaa/
}
\seealso{
Other isd: 
\code{\link{isd_read}()},
\code{\link{isd_stations}()},
\code{\link{isd}()}
}
\concept{isd}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_cache.r
\name{meteo_clear_cache}
\alias{meteo_clear_cache}
\title{Clear \emph{meteo} cached files}
\usage{
meteo_clear_cache(force = FALSE)
}
\arguments{
\item{force}{(logical) force delete. default: \code{FALSE}}
}
\description{
The \emph{meteo} functions use an aplication
}
\note{
This function will clear all cached \emph{meteo} files.
}
\seealso{
Other meteo: 
\code{\link{meteo_show_cache}()}
}
\concept{meteo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_ghcnd.R
\name{meteo_tidy_ghcnd_element}
\alias{meteo_tidy_ghcnd_element}
\title{Restructure element of ghcnd_search list}
\usage{
meteo_tidy_ghcnd_element(x, keep_flags = FALSE)
}
\arguments{
\item{x}{A dataframe with daily observations for a single monitor for a
single weather variable. This dataframe is one of the elements returned
by \code{\link[=ghcnd_search]{ghcnd_search()}}}

\item{keep_flags}{TRUE / FALSE for whether the user would like to keep all
the flags for each weather variable. The default is to not keep the
flags (FALSE). See the note below for more information on these flags.}
}
\value{
A dataframe reformatted to allow easy aggregation of all weather
variables for a single monitor.
}
\description{
This function restructures the output of \code{\link[=ghcnd_search]{ghcnd_search()}}
to add a column giving the variable name (\code{key}) and change the
name of the variable column to \code{value}. These changes facilitate
combining all elements from the list created by \code{\link[=ghcnd_search]{ghcnd_search()}},
to create a tidy dataframe of the weather observations from the station.
}
\author{
Brooke Anderson \email{brooke.anderson@colostate.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs_dimension_values}
\alias{gefs_dimension_values}
\title{This function is defunct.}
\usage{
gefs_dimension_values(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_distance.R
\name{meteo_spherical_distance}
\alias{meteo_spherical_distance}
\title{Calculate the distance between two locations}
\usage{
meteo_spherical_distance(lat1, long1, lat2, long2, units = "deg")
}
\arguments{
\item{lat1}{Latitude of the first location.}

\item{long1}{Longitude of the first location.}

\item{lat2}{Latitude of the second location.}

\item{long2}{Longitude of the second location.}

\item{units}{Units of the latitude and longitude values. Possible values
are:
\itemize{
\item \code{deg}: Degrees (default);
\item \code{rad}: Radians.
}}
}
\value{
A numeric value giving the distance (in kilometers) between the
pair of locations.
}
\description{
This function uses the haversine formula to calculate the great circle
distance between two locations, identified by their latitudes and longitudes.
}
\note{
This function assumes an earth radius of 6,371 km.
}
\examples{

meteo_spherical_distance(lat1 = -27.4667, long1 = 153.0217,
                         lat2 = -27.4710, long2 = 153.0234)
}
\author{
Alex Simmons \email{a2.simmons@qut.edu.au},
Brooke Anderson \email{brooke.anderson@colostate.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_distance.R
\name{meteo_distance}
\alias{meteo_distance}
\title{Find all monitors within a radius of a location}
\usage{
meteo_distance(
  station_data,
  lat,
  long,
  units = "deg",
  radius = NULL,
  limit = NULL
)
}
\arguments{
\item{station_data}{The output of \code{\link[=ghcnd_stations]{ghcnd_stations()}}, which is
a current list of weather stations available through NOAA for the GHCND
dataset. The format of this is a dataframe
with one row per weather station. Latitude and longitude for the station
locations should be in columns with the names "latitude" and "longitude",
consistent with the output from \code{\link[=ghcnd_stations]{ghcnd_stations()}}. To save time,
run the \code{ghcnd_stations} call and save the output to an object,
rather than rerunning the default every time (see the examples in
\code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}}).}

\item{lat}{Latitude of the location. Southern latitudes should be given
as negative values.}

\item{long}{Longitude of the location. Western longitudes should be given as
negative values.}

\item{units}{Units of the latitude and longitude values. Possible values
are:
\itemize{
\item \code{deg}: Degrees (default);
\item \code{rad}: Radians.
}}

\item{radius}{A numeric vector giving the radius (in kilometers) within which
to search for monitors near a location.}

\item{limit}{An integer giving the maximum number of monitors to include for
each location. The \code{x} closest monitors will be kept. Default is NULL
(pull everything available, within the radius if the radius is specified).}
}
\value{
A dataframe of weather stations near the location. This is the
single-location version of the return value for
\code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}}
}
\description{
This function will identify all weather stations with a specified radius of
a location. If no radius is given, the function will return a dataframe
of all available monitors, sorted by distance to the location. The
\code{limit} argument can be used to limit the output dataframe to the \code{x}
closest monitors to the location.
}
\examples{
\dontrun{
station_data <- ghcnd_stations()
meteo_distance(station_data, -33, 151, radius = 10, limit = 10)
meteo_distance(station_data, -33, 151, radius = 10, limit = 3)

# FIXME - units param is ignored
#meteo_distance(station_data, -33, 151, units = 'rad', radius = 10, limit = 3)
}
}
\author{
Alex Simmons \email{a2.simmons@qut.edu.au},
Brooke Anderson \email{brooke.anderson@colostate.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_combine.r
\name{ncdc_combine}
\alias{ncdc_combine}
\title{Coerce multiple outputs to a single data.frame object.}
\usage{
ncdc_combine(...)
}
\arguments{
\item{...}{Objects from another ncdc_* function.}
}
\value{
A data.frame
}
\description{
Coerce multiple outputs to a single data.frame object.
}
\examples{
\dontrun{
# data
out1 <- ncdc(datasetid='GHCND', locationid = 'FIPS:02', startdate = '2010-05-01',
enddate = '2010-05-31', limit=10)
out2 <- ncdc(datasetid='GHCND', locationid = 'FIPS:02', startdate = '2010-07-01',
enddate = '2010-07-31', limit=10)
ncdc_combine(out1, out2)

# data sets
out1 <- ncdc_datasets(datatypeid='TOBS')
out2 <- ncdc_datasets(datatypeid='PRCP')
ncdc_combine(out1, out2)

# data types
out1 <- ncdc_datatypes(datatypeid="ACMH")
out2 <- ncdc_datatypes(datatypeid='PRCP')
ncdc_combine(out1, out2)

# data categories
out1 <- ncdc_datacats(datacategoryid="ANNAGR")
out2 <- ncdc_datacats(datacategoryid='PRCP')
ncdc_combine(out1, out2)

# data locations
out1 <- ncdc_locs(locationcategoryid='ST', limit=52)
out2 <- ncdc_locs(locationcategoryid='CITY', sortfield='name', sortorder='desc')
ncdc_combine(out1, out2)

# data locations
out1 <- ncdc_locs_cats(startdate='1970-01-01')
out2 <- ncdc_locs_cats(locationcategoryid='CLIM_REG')
ncdc_combine(out1, out2)

# stations
out1 <- ncdc_stations(datasetid='GHCND', locationid='FIPS:12017',
stationid='GHCND:USC00084289')
out2 <- ncdc_stations(stationid='COOP:010008')
out3 <- ncdc_stations(datasetid='PRECIP_HLY', startdate='19900101',
enddate='19901231')
out4 <- ncdc_stations(datasetid='GHCND', locationid='FIPS:12017')
ncdc_combine(out1, out2, out3, out4)

# try to combine two different classes
out1 <- ncdc_locs_cats(startdate='1970-01-01')
out2 <- ncdc_stations(stationid='COOP:010008')
out3 <- ncdc_locs_cats(locationcategoryid='CLIM_REG')
ncdc_combine(out1, out2, out3)
}
}
\seealso{
Other ncdc: 
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc_stations}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc_datatypes.r
\name{ncdc_datatypes}
\alias{ncdc_datatypes}
\title{Get possible data types for a particular dataset}
\usage{
ncdc_datatypes(
  datasetid = NULL,
  datatypeid = NULL,
  datacategoryid = NULL,
  stationid = NULL,
  locationid = NULL,
  startdate = NULL,
  enddate = NULL,
  sortfield = NULL,
  sortorder = NULL,
  limit = 25,
  offset = NULL,
  token = NULL,
  ...
)
}
\arguments{
\item{datasetid}{(optional) Accepts a valid dataset id or a vector or list
of them. Data returned will be from the dataset specified.}

\item{datatypeid}{Accepts a valid data type id or a vector or list of data
type ids. (optional)}

\item{datacategoryid}{Optional. Accepts a valid data category id or a vector
or list of data category ids (although it is rare to have a data type with
more than one data category)}

\item{stationid}{Accepts a valid station id or a vector or list of
station ids}

\item{locationid}{Accepts a valid location id or a vector or list of
location ids (optional)}

\item{startdate}{(optional) Accepts valid ISO formated date (yyyy-mm-dd) or date time
(YYYY-MM-DDThh:mm:ss). Data returned will have data after the specified date. The
date range must be less than 1 year.}

\item{enddate}{(optional) Accepts valid ISO formated date (yyyy-mm-dd) or date time
(YYYY-MM-DDThh:mm:ss). Data returned will have data before the specified date. The
date range must be less than 1 year.}

\item{sortfield}{The field to sort results by. Supports id, name, mindate,
maxdate, and datacoverage fields (optional)}

\item{sortorder}{Which order to sort by, asc or desc. Defaults to
asc (optional)}

\item{limit}{Defaults to 25, limits the number of results in the response.
Maximum is 1000 (optional)}

\item{offset}{Defaults to 0, used to offset the resultlist (optional)}

\item{token}{This must be a valid token token supplied to you by NCDC's
Climate Data Online access token generator. (required) See
\strong{Authentication} section below for more details.}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}
(optional)}
}
\value{
A \code{data.frame} for all datasets, or a list of length two,
each with a data.frame
}
\description{
From the NOAA API docs: Describes the type of data, acts as a label.
For example: If it's 64 degrees out right now, then the data type is
Air Temperature and the data is 64.
}
\section{Authentication}{

Get an API key (aka, token) at https://www.ncdc.noaa.gov/cdo-web/token
You can pass your token in as an argument or store it one of two places:

\itemize{
\item your .Rprofile file with the entry
\code{options(noaakey = "your-noaa-token")}
\item your .Renviron file with the entry
\code{NOAA_KEY=your-noaa-token}
}

See \code{\link{Startup}} for information on how to create/find your
.Rrofile and .Renviron files
}

\examples{
\dontrun{
# Fetch available data types
ncdc_datatypes()

# Fetch more information about the ACMH data type id, or the ACSC
ncdc_datatypes(datatypeid="ACMH")
ncdc_datatypes(datatypeid="ACSC")

# datasetid, one or many
## ANNUAL should be replaced by GSOY, but both exist and give
## different answers
ncdc_datatypes(datasetid="ANNUAL")
ncdc_datatypes(datasetid="GSOY")
ncdc_datatypes(datasetid=c("ANNUAL", "PRECIP_HLY"))

# Fetch data types with the air temperature data category
ncdc_datatypes(datacategoryid="TEMP", limit=56)
ncdc_datatypes(datacategoryid=c("TEMP", "AUPRCP"))

# Fetch data types that support a given set of stations
ncdc_datatypes(stationid='COOP:310090')
ncdc_datatypes(stationid=c('COOP:310090','COOP:310184','COOP:310212'))

# Fetch data types that support a given set of loncationids
ncdc_datatypes(locationid='CITY:AG000001')
ncdc_datatypes(locationid=c('CITY:AG000001','CITY:AG000004'))
}
}
\references{
https://www.ncdc.noaa.gov/cdo-web/webservices/v2
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc_stations}()},
\code{\link{ncdc}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncdc.r
\name{ncdc}
\alias{ncdc}
\title{Search for and get NOAA NCDC data}
\usage{
ncdc(
  datasetid = NULL,
  datatypeid = NULL,
  stationid = NULL,
  locationid = NULL,
  startdate = NULL,
  enddate = NULL,
  sortfield = NULL,
  sortorder = NULL,
  limit = 25,
  offset = NULL,
  token = NULL,
  includemetadata = TRUE,
  add_units = FALSE,
  ...
)
}
\arguments{
\item{datasetid}{(required) Accepts a single valid dataset id. Data
returned will be from the dataset specified, see \code{\link{ncdc_datasets}}}

\item{datatypeid}{Accepts a valid data type id or a vector or list of data
type ids. (optional)}

\item{stationid}{Accepts a valid station id or a vector or list of station
ids}

\item{locationid}{Accepts a valid location id or a vector or list of
location ids (optional)}

\item{startdate}{(character/date) Accepts valid ISO formated date (yyyy-mm-dd)
or date time (YYYY-MM-DDThh:mm:ss). Data returned will have data after the
specified date. The date range must be less than 1 year. required.}

\item{enddate}{(character/date) Accepts valid ISO formated date (yyyy-mm-dd) or
date time (YYYY-MM-DDThh:mm:ss). Data returned will have data before the
specified date. The date range must be less than 1 year. required.}

\item{sortfield}{The field to sort results by. Supports id, name, mindate,
maxdate, and datacoverage fields (optional)}

\item{sortorder}{Which order to sort by, asc or desc. Defaults to
asc (optional)}

\item{limit}{Defaults to 25, limits the number of results in the response.
Maximum is 1000 (optional)}

\item{offset}{Defaults to 0, used to offset the resultlist (optional)}

\item{token}{This must be a valid token token supplied to you by NCDC's
Climate Data Online access token generator. (required) See
\strong{Authentication} section below for more details.}

\item{includemetadata}{Used to improve response time by preventing the
calculation of result metadata. Default: TRUE. This does not affect the
return object, in that the named part of the output list called "meta"
is still returned, but is NULL. In practice, I haven't seen response
time's improve, but perhaps they will for you.}

\item{add_units}{(logical) whether to add units information or not.
default: \code{FALSE}. If \code{TRUE}, after getting data from NOAA
we add a new column \code{units}. See "Adding units" in Details
for more}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}
(optional)}
}
\value{
An S3 list of length two, a slot of metadata (meta), and a slot
for data (data). The meta slot is a list of metadata elements, and the
data slot is a data.frame, possibly of length zero if no data is found. Note
that values in the data slot don't indicate their units by default, so you
will want to either use the \code{add_units} parameter (experimental, see Adding
units) or consult the documentation for each dataset to ensure you're using
the correct units.
}
\description{
Search for and get NOAA NCDC data
}
\details{
Note that NOAA NCDC API calls can take a long time depending on the call.
The NOAA API doesn't perform well with very long timespans, and will
time out and make you angry - beware.

Keep in mind that three parameters, datasetid, startdate, and enddate
are required.

Note that the default limit (no. records returned) is 25. Look at the
metadata in \verb{$meta} to see how many records were found. If more were
found than 25, you could set the parameter \code{limit} to something
higher than 25.
}
\section{Authentication}{

Get an API key (aka, token) at https://www.ncdc.noaa.gov/cdo-web/token
You can pass your token in as an argument or store it one of two places:

\itemize{
\item your .Rprofile file with the entry
\code{options(noaakey = "your-noaa-token")}
\item your .Renviron file with the entry
\code{NOAA_KEY=your-noaa-token}
}

See \code{\link{Startup}} for information on how to create/find your
.Rrofile and .Renviron files
}

\section{Flags}{

The attributes, or "flags", for each row of the output for data may have
a flag with it. Each \code{datasetid} has it's own set of flags. The
following are flag columns, and what they stand for. \code{fl_} is the
beginning of each flag column name, then one or more characters to describe
the flag, keeping it short to maintain a compact data frame. Some of
these fields are the same across datasetids. See the vignette
\code{vignette("rnoaa_attributes", "rnoaa")} for description of possible
values for each flag.
\itemize{
\item fl_c completeness
\item fl_d day
\item fl_m measurement
\item fl_q quality
\item fl_s source
\item fl_t time
\item fl_cmiss consecutive missing
\item fl_miss missing
\item fl_u units
}
}

\section{GSOM/GSOY Flags}{

Note that flags are different for GSOM and GSOY datasets. They have their
own set of flags per data class. See
\code{system.file("extdata/gsom.json", package = "rnoaa")} for GSOM
and \code{system.file("extdata/gsom.json", package = "rnoaa")} for GSOY.
Those are JSON files. The \code{\link[=system.file]{system.file()}} call gives you then path,
then read in with \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} which will give a data.frame
of the metadata. For more detailed info but plain text, open
\code{system.file("extdata/gsom_readme.txt", package = "rnoaa")}
and \code{system.file("extdata/gsoy_readme.txt", package = "rnoaa")}
in a text editor.
}

\section{Adding units}{

The \code{add_units} parameter is experimental - USE WITH CAUTION!
If \code{add_units=TRUE} we pull data from curated lists of data
used by matching by datasetid and data type.

We've attempted to gather as much information as possible on the many, many
data types across the many different NOAA data sets. However, we may have
got some things wrong, so make sure to double check data you get if you
do add units.

Get in touch if you find some units that are wrong or missing, and
if you are able to help correct information.
}

\examples{
\dontrun{
# GHCN-Daily (or GHCND) data, for a specific station
ncdc(datasetid='GHCND', stationid='GHCND:USW00014895',
   startdate = '2013-10-01', enddate = '2013-12-01')
### also accepts dates as class Date
ncdc(datasetid='GHCND', stationid='GHCND:USW00014895',
   startdate = as.Date('2013-10-01'), enddate = as.Date('2013-12-01'))

# GHCND data, for a location by FIPS code
ncdc(datasetid='GHCND', locationid = 'FIPS:02', startdate = '2010-05-01',
   enddate = '2010-05-10')

# GHCND data from October 1 2013 to December 1 2013
ncdc(datasetid='GHCND', startdate = '2013-10-01', enddate = '2013-10-05')

# GHCN-Monthly (or GSOM) data from October 1 2013 to December 1 2013
ncdc(datasetid='GSOM', startdate = '2013-10-01', enddate = '2013-12-01')
ncdc(datasetid='GSOM', startdate = '2013-10-01', enddate = '2013-12-01',
   stationid = "GHCND:AE000041196")

# Normals Daily (or NORMAL_DLY) GHCND:USW00014895 dly-tmax-normal data
ncdc(datasetid='NORMAL_DLY', stationid='GHCND:USW00014895',
   startdate = '2010-05-01', enddate = '2010-05-10')

# Dataset, and location in Australia
ncdc(datasetid='GHCND', locationid='FIPS:AS', startdate = '2010-05-01',
    enddate = '2010-05-31')

# Dataset, location and datatype for PRECIP_HLY data
ncdc(datasetid='PRECIP_HLY', locationid='ZIP:28801', datatypeid='HPCP',
   startdate = '2010-05-01', enddate = '2010-05-10')

# multiple datatypeid's
ncdc(datasetid='PRECIP_HLY', datatypeid = 'HPCP',
   startdate = '2010-05-01', enddate = '2010-05-10')

# multiple locationid's
ncdc(datasetid='PRECIP_HLY', locationid=c("FIPS:30103", "FIPS:30091"),
   startdate = '2010-05-01', enddate = '2010-05-10')

# Dataset, location, station and datatype
ncdc(datasetid='PRECIP_HLY', locationid='ZIP:28801',
   stationid='COOP:310301', datatypeid='HPCP',
   startdate = '2010-05-01', enddate = '2010-05-10')

# Dataset, location, and datatype for GHCND
ncdc(datasetid='GHCND', locationid='FIPS:BR', datatypeid='PRCP',
   startdate = '2010-05-01', enddate = '2010-05-10')

# Normals Daily GHCND dly-tmax-normal data
ncdc(datasetid='NORMAL_DLY', datatypeid='dly-tmax-normal',
   startdate = '2010-05-01', enddate = '2010-05-10')

# Normals Daily GHCND:USW00014895 dly-tmax-normal
ncdc(datasetid='NORMAL_DLY', stationid='GHCND:USW00014895',
   datatypeid='dly-tmax-normal',
   startdate = '2010-05-01', enddate = '2010-05-10')

# Hourly Precipitation data for ZIP code 28801
ncdc(datasetid='PRECIP_HLY', locationid='ZIP:28801', datatypeid='HPCP',
   startdate = '2010-05-01', enddate = '2010-05-10')

# 15 min Precipitation data for ZIP code 28801
ncdc(datasetid='PRECIP_15', datatypeid='QPCP',
   startdate = '2010-05-01', enddate = '2010-05-02')

# Search the NORMAL_HLY dataset
ncdc(datasetid='NORMAL_HLY', stationid = 'GHCND:USW00003812',
   startdate = '2010-05-01', enddate = '2010-05-10')

# Search the GSOY dataset
ncdc(datasetid='ANNUAL', locationid='ZIP:28801', startdate = '2010-05-01',
   enddate = '2010-05-10')

# Search the NORMAL_ANN dataset
ncdc(datasetid='NORMAL_ANN', datatypeid='ANN-DUTR-NORMAL',
   startdate = '2010-01-01', enddate = '2010-01-01')

# Include metadata or not
ncdc(datasetid='GHCND', stationid='GHCND:USW00014895',
   startdate = '2013-10-01', enddate = '2013-12-01')
ncdc(datasetid='GHCND', stationid='GHCND:USW00014895',
   startdate = '2013-10-01', enddate = '2013-12-01', includemetadata=FALSE)

# Many stationid's
stat <- ncdc_stations(startdate = "2000-01-01", enddate = "2016-01-01")
## find out what datasets might be available for these stations
ncdc_datasets(stationid = stat$data$id[10])
## get some data
ncdc(datasetid = "GSOY", stationid = stat$data$id[1:10],
   startdate = "2010-01-01", enddate = "2011-01-01")
}

\dontrun{
# NEXRAD2 data
## doesn't work yet
ncdc(datasetid='NEXRAD2', startdate = '2013-10-01', enddate = '2013-12-01')
}
}
\seealso{
Other ncdc: 
\code{\link{ncdc_combine}()},
\code{\link{ncdc_datacats}()},
\code{\link{ncdc_datasets}()},
\code{\link{ncdc_datatypes}()},
\code{\link{ncdc_locs_cats}()},
\code{\link{ncdc_locs}()},
\code{\link{ncdc_plot}()},
\code{\link{ncdc_stations}()}
}
\concept{ncdc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpc.R
\name{cpc_prcp}
\alias{cpc_prcp}
\title{Precipitation data from NOAA Climate Prediction Center (CPC)}
\usage{
cpc_prcp(date, us = FALSE, drop_undefined = FALSE, ...)
}
\arguments{
\item{date}{(date/character) date in YYYY-MM-DD format}

\item{us}{(logical) US data only? default: \code{FALSE}}

\item{drop_undefined}{(logical) drop undefined precipitation
values (values in the \code{precip} column in the output data.frame).
default: \code{FALSE}}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a data.frame, with columns:
\itemize{
\item lon - longitude (0 to 360)
\item lat - latitude (-90 to 90)
\item precip - precipitation (in mm) (see Details for more information)
}
}
\description{
Precipitation data from NOAA Climate Prediction Center (CPC)
}
\details{
Rainfall data for the world (1979-present, resolution 50 km), and
the US (1948-present, resolution 25 km).
}
\note{
See \link{cpc_cache} for managing cached files
}
\section{Data processing in this function}{

Internally we multiply all precipitation measurements by 0.1 as
per the CPC documentation.

Values of -99.0 are classified as "undefined". These values can be
removed by setting \code{drop_undefined = TRUE} in the \code{cpc_prcp}
function call. These undefined values are not dropped by default -
so do remember to set \code{drop_undefined = TRUE} to drop them; or
you can easily do it yourself by e.g., \code{subset(x, precip >= 0)}
}

\examples{
\dontrun{
x = cpc_prcp(date = "2017-01-15")
cpc_prcp(date = "2015-06-05")
cpc_prcp(date = "2017-01-15")
cpc_prcp(date = "2005-07-09")
cpc_prcp(date = "1979-07-19")

# United States data only
cpc_prcp(date = "2005-07-09", us = TRUE)
cpc_prcp(date = "2009-08-03", us = TRUE)
cpc_prcp(date = "1998-04-23", us = TRUE)

# drop undefined values (those given as -99.0)
cpc_prcp(date = "1998-04-23", drop_undefined = TRUE)
}
}
\references{
https://www.cpc.ncep.noaa.gov/
https://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP
https://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/DOCU/PRCP_CU_GAUGE_V1.0CONUS_0.25deg.README
https://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_GLB/DOCU/PRCP_CU_GAUGE_V1.0GLB_0.50deg_README.txt
https://psl.noaa.gov/data/gridded/data.unified.daily.conus.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnoaa-package.r
\docType{package}
\name{rnoaa-package}
\alias{rnoaa-package}
\alias{rnoaa}
\title{rnoaa}
\description{
rnoaa is an R interface to NOAA climate data.
}
\section{Data Sources}{

Many functions in this package interact with the National Climatic Data
Center application programming interface (API) at
https://www.ncdc.noaa.gov/cdo-web/webservices/v2, all of
which functions start with \code{ncdc_}. An access token, or API key, is
required to use all the \code{ncdc_} functions. The key is required by NOAA,
not us. Go to the link given above to get an API key.

More NOAA data sources are being added through time. Data sources and their
function prefixes are:
\itemize{
\item \verb{buoy_*} - NOAA Buoy data from the National Buoy Data Center
\item \verb{ghcnd_*}/\verb{meteo_*} - GHCND daily data from NOAA
\item \verb{isd_*} - ISD/ISH data from NOAA
\item \verb{homr_*} - Historical Observing Metadata Repository (HOMR)
vignette
\item \verb{ncdc_*} - NOAA National Climatic Data Center (NCDC) vignette
(examples)
\item \code{sea_ice} - Sea ice vignette
\item \code{storm_} - Storms (IBTrACS) vignette
\item \code{swdi} - Severe Weather Data Inventory (SWDI) vignette
\item \code{tornadoes} - From the NOAA Storm Prediction Center
\item \verb{argo_*} - Argo buoys
\item \code{coops_search} - NOAA CO-OPS - tides and currents data
\item \code{cpc_prcp} - rainfall data from the NOAA Climate
Prediction Center (CPC)
\item \code{arc2} - rainfall data from Africa Rainfall Climatology
version 2
\item \code{bsw} - Blended sea winds (BSW)
\item \code{ersst} - NOAA Extended Reconstructed Sea Surface
Temperature (ERSST) data
\item \code{lcd} - Local Climitalogical Data from NOAA
}
}

\section{Where data comes from and government shutdowns}{


Government shutdowns can greatly affect data sources in this package.
The following is a breakdown of the functions that fetch data by
HTTP vs. FTP - done this way as we've noticed that during the ealry 2019
border wall shutdown most FTP services were up, while those that were down
were HTTP; though not all HTTP services were down.
\itemize{
\item HTTP info: https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol
\item FTP info: https://en.wikipedia.org/wiki/File_Transfer_Protocol
}

HTTP services (whether service is/was up or down during early 2019 shutdown)
\itemize{
\item \verb{buoy_*} - Up
\item \verb{homr_*} - Up
\item \verb{ncdc_*} - Down
\item \code{swdi} - Down
\item \code{tornadoes} - Down
\item \verb{argo_*} - Up (all HTTP except two fxns, see also FTP below)
\item \code{coops_search} - Up
\item \code{ersst} - Down
\item \code{lcd} - Down
\item \verb{se_*} - Down
}

FTP services (whether service is/was up or down during early 2019 shutdown)
\itemize{
\item \verb{ghcnd_*} - Up
\item \verb{isd_*} - Up
\item \code{sea_ice} - Up
\item \code{storm_} - Up
\item \verb{argo_*} - Up (only two fxns: \code{\link[=argo]{argo()}}, \code{\link[=argo_buoy_files]{argo_buoy_files()}})
\item \code{cpc_prcp} - Up
\item \code{arc2} - Up
\item \code{bsw} - Up
}

We've tried to whenever possible detect whether a service is error
due to a government shutdown and give a message saying so. If you know
a service is down that rnoaa interacts with but we don't fail well
during a shutdown let us know.
}

\section{A note about NCDF data}{


Some functions use netcdf files, including:
\itemize{
\item \code{ersst}
\item \code{buoy}
\item \code{bsw}
\item \code{argo}
}

You'll need the \code{ncdf4} package for those functions, and those only.
\code{ncdf4} is in Suggests in this package, meaning you only need
\code{ncdf4} if you are using any of the functions listed above. You'll get
an informative error telling you to install \code{ncdf4} if you don't have
it and you try to use the those functions. Installation of \code{ncdf4}
should be straightforward on any system.
}

\section{The \code{meteo} family of functions}{


The \code{meteo} family of functions are prefixed with \code{meteo_} and
provide a set of helper functions to:
\itemize{
\item Identify candidate stations from a latitude/longitude pair
\item Retrieve complete data for one or more stations (\code{meteo_coverage()})
}
}

\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swdi.r
\name{swdi}
\alias{swdi}
\title{Get NOAA data for the Severe Weather Data Inventory (SWDI)}
\usage{
swdi(
  dataset = NULL,
  format = "xml",
  startdate = NULL,
  enddate = NULL,
  limit = 25,
  offset = NULL,
  radius = NULL,
  center = NULL,
  bbox = NULL,
  tile = NULL,
  stat = NULL,
  id = NULL,
  filepath = NULL,
  ...
)
}
\arguments{
\item{dataset}{Dataset to query. See below for details.}

\item{format}{File format to download. One of xml, csv, shp, or kmz.}

\item{startdate}{Start date. See details.}

\item{enddate}{End date. See details.}

\item{limit}{Number of results to return. Defaults to 25. Any number
from 1 to 10000000. Time out issues likely to occur at higher limits.}

\item{offset}{Any number from 1 to 10000000. Default is NULL, no offset,
start from 1.}

\item{radius}{Search radius in miles (current limit is 15 miles).
BEWARE: As far as we know, this parameter doesn't do anything, or at
least does not in fact limit the search to the given radius.
DO NOT USE.}

\item{center}{Center coordinate in lon,lat decimal degree format,
e.g.: c(-95.45,36.88)}

\item{bbox}{Bounding box in format of minLon,minLat,maxLon,maxLat,
e.g.: c(-91,30,-90,31)}

\item{tile}{Coordinate in lon,lat decimal degree format,
e.g.: c(-95.45,36.88). The lat/lon values are rounded to the nearest tenth
of degree. For the above example, the matching tile would contain values
from -95.4500 to -95.5499 and 36.8500 to 36.9499}

\item{stat}{One of count or tilesum:$longitude,$latitude. Setting
stat='count' returns number of results only (no actual data).
stat='tilesum:$longitude,$latitude' returns daily feature counts for
a tenth of a degree grid centered at the nearest tenth of a degree to
the supplied values.}

\item{id}{An identifier, e.g., 533623. Not sure how you find these ids?}

\item{filepath}{If kmz or shp chosen the file name and optionally path to
write to. Ignored format=xml or format=csv (optional)}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET} (optional)}
}
\value{
If xml or csv chosen, a list of length three, a slot of metadata
(meta), a slot for data (data), and a slot for shape file data with a
single column 'shape'. The meta slot is a list of metadata elements, and
the data slot is a data.frame, possibly of length zero if no data is found.

If kmz or shp chosen, the file is downloaded to your machine and a message
is printed.
}
\description{
Get NOAA data for the Severe Weather Data Inventory (SWDI)
}
\details{
Options for the dataset parameter. One of (and their data formats):
\itemize{
\item nx3tvs NEXRAD Level-3 Tornado Vortex Signatures (point)
\item nx3meso NEXRAD Level-3 Mesocyclone Signatures (point)
\item nx3hail NEXRAD Level-3 Hail Signatures (point)
\item nx3structure NEXRAD Level-3 Storm Cell Structure Information (point)
\item plsr Preliminary Local Storm Reports (point)
\item warn Severe Thunderstorm, Tornado, Flash Flood and Special Marine
warnings (polygon)
\item nldn Lightning strikes from Vaisala. Available to government and
military users only. If you aren't one of those, you'll get a 400 status
stop message if you request data from this dataset (point)
}

For startdate and enddate, the date range syntax is 'startDate:endDate' or
special option of 'periodOfRecord'. Note that startDate is inclusive and
endDate is exclusive. All dates and times are in GMT. The current limit of
the date range size is one year.

All latitude and longitude values for input parameters and output data are
in the WGS84 datum.
}
\examples{
\dontrun{
# Search for nx3tvs data from 5 May 2006 to 6 May 2006
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506')

# Get all 'nx3tvs' near latitude = 32.7 and longitude = -102.0
swdi(dataset='nx3tvs', startdate='20060506', enddate='20060507',
center=c(-102.0,32.7))

# use an id
swdi(dataset='warn', startdate='20060506', enddate='20060507', id=533623)

# Get all 'plsr' within the bounding box (-91,30,-90,31)
swdi(dataset='plsr', startdate='20060505', enddate='20060510',
bbox=c(-91,30,-90,31))

# Get all 'nx3tvs' within the tile -102.1/32.6 (-102.15,32.55,-102.25,32.65)
swdi(dataset='nx3tvs', startdate='20060506', enddate='20060507',
tile=c(-102.12,32.62))

# Counts
## Note: stat='count' will only return metadata, nothing in the data or shape slots
## Note: stat='tilesum:...' returns counts in the data slot for each date for that tile,
##   and shape data
## Get number of 'nx3tvs' near latitude = 32.7 and longitude = -102.0
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060516',
center=c(-102.0,32.7), stat='count')

## Get daily count nx3tvs features on .1 degree grid centered at latitude = 32.7
## and longitude = -102.0
swdi(dataset='nx3tvs', startdate='20060505', enddate='20090516',
stat='tilesum:-102.0,32.7')

# CSV format
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506', format='csv')

# SHP format
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506', format='shp',
   filepath='myfile')

# KMZ format
swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506', format='kmz',
   filepath='myfile.kmz')

# csv output to SpatialPointsDataFrame
res <- swdi(dataset='nx3tvs', startdate='20060505', enddate='20060506', format="csv")
library('sp')
coordinates(res$data) <- ~lon + lat
res$data
class(res$data)
}
}
\references{
https://www.ncdc.noaa.gov/ncei-severe-weather-data-inventory
https://www.ncdc.noaa.gov/swdiws/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buoy.R
\name{buoy}
\alias{buoy}
\alias{buoys}
\alias{buoy_stations}
\title{Get NOAA buoy data from the National Buoy Data Center}
\usage{
buoy(dataset, buoyid, year = NULL, datatype = NULL, ...)

buoys(dataset)

buoy_stations(refresh = FALSE, ...)
}
\arguments{
\item{dataset}{(character) Dataset name to query. See below for Details.
Required}

\item{buoyid}{Buoy ID, can be numeric/integer/character. Required}

\item{year}{(integer) Year of data collection. Optional. Note there is
a special value \code{9999} that, if found, contains the most up to date
data.}

\item{datatype}{(character) Data type, one of 'c', 'cc', 'p', 'o'. Optional}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}
Optional. A number of different HTTP requests are made internally, but
we only pass this on to the request to get the netcdf file in the internal
function \code{get_ncdf_file()}}

\item{refresh}{(logical) Whether to use cached data (\code{FALSE}) or get
new data (\code{FALSE}). Default: \code{FALSE}}
}
\value{
If netcdf data has lat/lon variables, then we'll parse into a
tidy data.frame. If not, we'll give back the ncdf4 object for the user
to parse (in which case the data.frame will be empty).
}
\description{
Get NOAA buoy data from the National Buoy Data Center
}
\details{
Functions:
\itemize{
\item buoy_stations - Get buoy stations. A cached version of the dataset
is available in the package. Beware, takes a long time to run if you
do \code{refresh = TRUE}
\item buoys - Get available buoys given a dataset name
\item buoy - Get data given some combination of dataset name, buoy ID,
year, and datatype
}

Options for the dataset parameter. One of:
\itemize{
\item adcp - Acoustic Doppler Current Profiler data
\item adcp2 - MMS Acoustic Doppler Current Profiler data
\item cwind - Continuous Winds data
\item dart - Deep-ocean Assessment and Reporting of Tsunamis data
\item mmbcur - Marsh-McBirney Current Measurements data
\item ocean - Oceanographic data
\item pwind - Peak Winds data
\item stdmet - Standard Meteorological data
\item swden - Spectral Wave Density data with Spectral Wave Direction data
\item wlevel - Water Level data
}
}
\examples{
\dontrun{
if (crul::ok("https://dods.ndbc.noaa.gov/thredds", timeout_ms = 1000)) {

# Get buoy station information
x <- buoy_stations()
# refresh stations as needed, takes a while to run
# you shouldn't need to update very often
# x <- buoy_stations(refresh = TRUE)
if (interactive() && requireNamespace("leaflet")){
library("leaflet")
z <- leaflet(data = na.omit(x))
z <- leaflet::addTiles(z)
leaflet::addCircles(z, ~lon, ~lat, opacity = 0.5)
}

# year=9999 to get most current data - not always available
buoy(dataset = "swden", buoyid = 46012, year = 9999)

# Get available buoys
buoys(dataset = 'cwind')

# Get data for a buoy
## if no year or datatype specified, we get the first file
buoy(dataset = 'cwind', buoyid = 46085)

# Including specific year
buoy(dataset = 'cwind', buoyid = 41001, year = 1999)

# Including specific year and datatype
buoy(dataset = 'cwind', buoyid = 45005, year = 2008, datatype = "c")
buoy(dataset = 'cwind', buoyid = 41001, year = 1997, datatype = "c")

# Other datasets
buoy(dataset = 'ocean', buoyid = 41029)

# curl debugging
buoy(dataset = 'cwind', buoyid = 46085, verbose = TRUE)

# some buoy ids are character, case doesn't matter, we'll account for it
buoy(dataset = "stdmet", buoyid = "VCAF1")
buoy(dataset = "stdmet", buoyid = "wplf1")
buoy(dataset = "dart", buoyid = "dartu")

}
}
}
\references{
http://www.ndbc.noaa.gov/, http://dods.ndbc.noaa.gov/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo-autoplot.R
\name{autoplot.meteo_coverage}
\alias{autoplot.meteo_coverage}
\title{autoplot method for meteo_coverage objects}
\usage{
\method{autoplot}{meteo_coverage}(mateo_coverage, old_style = FALSE)
}
\arguments{
\item{mateo_coverage}{the object returned from \code{\link[=meteo_coverage]{meteo_coverage()}}}

\item{old_style}{(logical) create the old style of plots, which is faster, but
does not plot gaps to indicate missing data}
}
\value{
A ggplot2 plot
}
\description{
autoplot method for meteo_coverage objects
}
\details{
see \code{\link[=meteo_coverage]{meteo_coverage()}} for examples
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seaice.r
\name{seaiceeurls}
\alias{seaiceeurls}
\title{Make all urls for sea ice data}
\usage{
seaiceeurls(yr = NULL, mo = NULL, pole = NULL, format = "shp")
}
\arguments{
\item{yr}{(numeric) a year}

\item{mo}{(character) a month, as character abbrevation of a month}

\item{pole}{(character) one of S (south) or N (north)}
}
\value{
A vector of urls (character)
}
\description{
Make all urls for sea ice data
}
\examples{
\dontrun{
# Get all urls
seaiceeurls()

# for some range of years
seaiceeurls(yr = 1980:1983)
seaiceeurls(yr = 1980, mo = c("Jan", "Feb", "Mar"))
seaiceeurls(yr = 1980:1983, mo = c("Jan", "Apr", "Oct"))

# Get urls for Feb of all years, both S and N poles
seaiceeurls(mo='Feb')

# Get urls for Feb of all years, just S pole
seaiceeurls(mo='Feb', pole='S')

# Get urls for Feb of 1980, just S pole
seaiceeurls(yr=1980, mo='Feb', pole='S')

# GeoTIFF
seaiceeurls(yr=1980, mo='Feb', pole='S', format = "geotiff")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/argo.R
\name{argo}
\alias{argo}
\alias{argo_search}
\alias{argo_files}
\alias{argo_qwmo}
\alias{argo_plan}
\alias{argo_buoy_files}
\title{Get Argo buoy data}
\usage{
argo_search(
  func = NULL,
  of = NULL,
  qwmo = NULL,
  wmo = NULL,
  box = NULL,
  area = NULL,
  around = NULL,
  year = NULL,
  yearmin = NULL,
  yearmax = NULL,
  month = NULL,
  monthmin = NULL,
  monthmax = NULL,
  lr = NULL,
  from = NULL,
  to = NULL,
  dmode = NULL,
  pres_qc = NULL,
  temp_qc = NULL,
  psal_qc = NULL,
  doxy_qc = NULL,
  ticket = NULL,
  limit = 10,
  ...
)

argo_files(wmo = NULL, cyc = NULL, ...)

argo_qwmo(qwmo, limit = 10, ...)

argo_plan(...)

argo_buoy_files(dac, id, ...)

argo(dac, id, cycle, dtype, ...)
}
\arguments{
\item{func}{A function, one of n, np, nf, coord, fullcoord, list, ftplist, ticket, version}

\item{of}{of string}

\item{qwmo}{qwmo string}

\item{wmo}{wmo string. mandatory when using \code{argo_files}}

\item{box}{Bounding box, of the form: A) lon, lat for geographical coordinates of
the center of a box, or B) min lon, min lat, width, height, for geographical
coordinates of the center of the box and its width and height, and the longitude must
given between -180W and 180E. Width and height are in degrees of longitude and latitude.}

\item{area}{(integer/character), One of 0, 1, or 2, but can be in combination. e.g. 0, '0,2'
See Details.}

\item{around}{(character) Selects profiles located around a given center point. List of
3 or 4 numerical values depending on how the center point need to be specified:
e.g., '-40,35,100', '6900678,2,200', '6900678,2,200,30'. See Details}

\item{year}{restrict profiles sampled on a single, or a list of given years. One or a
comma separated list of numerical value(s) higher than 0 and lower than 9999.}

\item{yearmin, yearmax}{restrict profiles sampled before (yearmax) and/or after (yearmin)
a given year. A numerical value higher than 0 and lower than 9999. cannot be applied with
the other restriction parameter \code{year}}

\item{month}{restrict profiles sampled on a single, or a list of given month(s). One or a
comma separated list of numerical value(s) higher or equal to 1 and lower or equal to 12.
The month convention is a standard one: January is 1, February is 2, ... December is 12.}

\item{monthmin, monthmax}{restrict profiles sampled before (monthmax) and/or after (monthmin)
a given month. Higher or equal to 1 and lower or equal to 12. The month convention is a
standard one: January is 1, February is 2, ... December is 12. These restrictions cannot be
applied with the other restriction parameter month. At this time, these parameters are not
circular, so that the restriction chain: monthmin=12&monthmax=2 will through an error
and not select December to February profiles. To do so, you need to use a coma separated
list of months using the month restriction parameter.}

\item{lr}{restriction allows you to impose the last report (hence lr) date in days. A
numerical value in days between 1 (profiles sampled yesterday) and 60 (profiles sampled
over the last 60 days). This restriction allows a simple selection of the so-called
'active' floats, ie those which reported a profiles over the last 30 days.}

\item{from, to}{select profiles sampled before (to) and/or after (from) an explicit date
(included). The date is specified following the format: YYYYMMDD, ie the year, month
and day numbers.}

\item{dmode}{(character) imposes a restriction on the Data Mode of profiles. A single value or
a coma separated list of characters defining the Data Mode to select. It can be: R for
"Real Time", A for "Real Time with Adjusted value" and D for "Delayed Mode". See Details.}

\item{pres_qc, temp_qc, psal_qc, doxy_qc}{Quality control. Imposes a restriction on the profile
data quality flag. For a given variable PARAM which can be: pres (pressure),
temp (temperature), psal (salinity) or doxy (oxygen), this restriction selects profiles
having one or a coma separated list of data quality flag. See Details.}

\item{ticket}{(numeric) select profiles with or without a ticket filled in the database. A
value: 0 (no ticket) or 1 (has a ticket). See
http://www.ifremer.fr/lpo/naarc/m/docs/api/database.html for more details.}

\item{limit}{(integer) number to return}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}. Optional}

\item{cyc}{a cycle number}

\item{dac}{(character) Data assembly center code}

\item{id}{(numeric) Buoy identifier}

\item{cycle}{(numeric) Cycle number}

\item{dtype}{(character) Data type, one of \code{D} for delayed, or \code{R} for real-time}
}
\description{
Get Argo buoy data
}
\details{
\code{area} parameter definitions:
\itemize{
\item Value 0 selects profiles located in the North-Atlantic ocean north of 20S
and not in areas 1 and 2.
\item Value 1 selects profiles located in the Mediterranean Sea.
\item Value 2 selects profiles located in the Nordic Seas.
}

\code{around} parameter definitions:
\itemize{
\item Specification 1: The location is specified with specific geographical
coordinates in the following format: around=longitude,latitude,distance - The longitude
must given between -180W and 180E and the distance is in kilometers.
\item Specification 2: The location is the one of an existing profile in the database.
It is thus specified with a float WMO and a cycle number: around=wmo,cyc,distance
This specification can take an optional fourth value specifying the time range in days
around the specified profile.
}

\code{dmode} parameter definitions:
\itemize{
\item Data from Argo floats are transmitted from the float, passed through processing and
automatic quality control procedures. These profiles have a Data Mode called: real-time data.
\item The data are also issued to the Principle Investigators who apply other procedures to
check data quality returned to the global data centre within 6 to 12 months. These profiles
have a Data Mode called: delayed mode data.
\item The adjustments applied to delayed-data may also be applied to real-time data, to
correct sensor drifts for real-time users. These profiles have a Data Mode called: real
time data with adjusted values.
}

\code{*_qc} parameter definitions:
This information was extracted from the netcdf profile variable PROFILE_PARAM_QC. Once
quality control procedures have been applied, a synthetic flag is assigned for each
parameter of each profile under this variable in netcdf files. It indicates the fraction
n of profile levels with good data. It can take one of the following values:
\itemize{
\item A or F: All (n=100\%) or none (n=0\%) of the profile levels contain good data,
\item B,C,D,E: n is in one of the intermediate range: 75-100, 50-75, 25-50 or 0-25
\item empty: No QC was performed.
}
}
\section{File storage}{

We use \pkg{rappdirs} to store files, see
\code{\link[rappdirs]{user_cache_dir}} for how we determine the directory on
your machine to save files to, and run
\code{rappdirs::user_cache_dir("rnoaa/argo")} to get that directory.
}

\section{API Status}{

The API weas down as of 2019-11-07, and probably some time before that. The
following functions were defunct:

\itemize{
\item argo_search
\item argo_files
\item argo_qwmo
\item argo_plan
}

These functions are working again as of 2020-06-12.
}

\examples{
\dontrun{
# Search Argo metadata
## Number of profiles
argo_search("np", limit = 3)
## Number of floats
argo_search("nf", limit = 3)
## Number of both profiles and floats
argo_search("n", limit = 3)
## return the coordinates in time and space of profiles
argo_search("coord", limit = 3)
## return the coordinates in time and space of profiles, plus other metadata
argo_search("fullcoord", limit = 3)

## List various things, e.g,...
### data assembly centers
argo_search("list", "dac")
### data modes
argo_search("list", "dmode", limit = 5)
### World Meteorological Organization unique float ID's
argo_search("list", "wmo", limit = 5)
### Profile years
argo_search("list", "year", limit = 5)

## coord or fullcoord with specific buoy id
argo_search("coord", wmo = 13857, limit = 3)
argo_search("fullcoord", wmo = 13857, limit = 3)

# Spatial search
### search by bounding box (see param def above)
argo_search("coord", box = c(-40, 35, 3, 2))
### search by area
argo_search("coord", area = 0)
### search by around
argo_search("coord", around = '-40,35,100')

# Time based search
### search by year
argo_search("coord", year = 2006)
### search by yearmin and yearmax
argo_search("coord", yearmin = 2007)
argo_search("coord", yearmin = 2007, yearmax = 2009)
### search by month
argo_search("coord", month = '12,1,2')
### search by from or to
argo_search("coord", from = 20090212)
argo_search("coord", to = 20051129)

# Data mode search
argo_search("coord", dmode = "R")
argo_search("coord", dmode = "R,A")

# Data quality based search
argo_search("coord", pres_qc = "A,B")
argo_search("coord", temp_qc = "A")
argo_search("coord", pres_qc = "A", temp_qc = "A")

# Ticket search
argo_search("coord", ticket = 0)

## Search on partial float id number
argo_qwmo(qwmo = 49)
argo_qwmo(qwmo = 49, limit = 2)

## Get files
argo_files(wmo = 13857)
argo_files(wmo = 13857, cyc = 12)
argo_files(wmo = 13857, cyc = 45)

## Get planned buoys data, accepts no parameters
argo_plan()

# Get files for a buoy, must specify data assembly center (dac)
argo_buoy_files(dac = "bodc", id = 1901309)
argo_buoy_files(dac = "kma", id = 2900308)

# Get data
x <- argo_buoy_files(dac = "meds", id = 4900881)
argo(dac = "meds", id = 4900881, cycle = 127, dtype = "D")
}
}
\references{
http://www.ifremer.fr/lpo/naarc/m/docs/api/howto.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seaice.r
\name{readshpfile}
\alias{readshpfile}
\title{Function to read shapefiles}
\usage{
readshpfile(x, storepath = NULL)
}
\arguments{
\item{x}{A url}

\item{storepath}{Path to store data in}
}
\value{
An object of class sp
}
\description{
Function to read shapefiles
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lcd.R
\name{lcd}
\alias{lcd}
\title{Local Climatological Data from NOAA}
\usage{
lcd(station, year, col_types = NULL, ...)
}
\arguments{
\item{station}{(character) station code, e.g., "02413099999". we will allow
integer/numeric passed here, but station ids can have leading zeros, so
it's a good idea to keep stations as character class. required}

\item{year}{(integer) year, e.g., 2017. required}

\item{col_types}{(named character vector) defaults to NULL. Use this argument
to change the returned column type. For example,"character" instead of
"numeric". See or use \link{lcd_columns} to create a named vector with allowed
column names. If the user specified type is not compatible, the function
will choose a type automatically and raise a message. optional}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a data.frame with many columns. the first 10 are metadata:
\itemize{
\item station
\item date
\item latitude
\item longitude
\item elevation
\item name
\item report_type
\item source
}

And the rest should be all data columns. The first part of many column
names is the time period, being one of:
\itemize{
\item hourly
\item daily
\item monthly
\item shortduration
}

So the variable you are looking for may not be the first part of the column
name
}
\description{
Local Climatological Data from NOAA
}
\note{
See \link{lcd_cache} for managing cached files
}
\examples{
\dontrun{
x = lcd(station = "01338099999", year = 2017)
lcd(station = "01338099999", year = 2015)
lcd(station = "02413099999", year = 2009)
lcd(station = "02413099999", year = 2001)

# pass curl options
lcd(station = "02413099999", year = 2002, verbose = TRUE)
}
}
\references{
Docs:
https://www.ncei.noaa.gov/data/local-climatological-data/doc/LCD_documentation.pdf
Data comes from:
https://www.ncei.noaa.gov/data/local-climatological-data/access
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/storm_shp.R
\name{storm_shp}
\alias{storm_shp}
\alias{storm_shp_read}
\title{storm shp files}
\usage{
storm_shp(...)

storm_shp_read(...)
}
\description{
storm shp files

read storm shp files
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd.R
\name{isd}
\alias{isd}
\title{Get and parse NOAA ISD/ISH data}
\usage{
isd(
  usaf,
  wban,
  year,
  overwrite = TRUE,
  cleanup = TRUE,
  additional = TRUE,
  parallel = FALSE,
  cores = getOption("cl.cores", 2),
  progress = FALSE,
  force = FALSE,
  ...
)
}
\arguments{
\item{usaf, wban}{(character) USAF and WBAN code. Required}

\item{year}{(numeric) One of the years from 1901 to the current year.
Required.}

\item{overwrite}{(logical) To overwrite the path to store files in or not,
Default: \code{TRUE}}

\item{cleanup}{(logical) If \code{TRUE}, remove compressed \code{.gz} file
at end of function execution. Processing data takes up a lot of time, so we
cache a cleaned version of the data. Cleaning up will save you on disk
space. Default: \code{TRUE}}

\item{additional}{(logical) include additional and remarks data sections
in output. Default: \code{TRUE}. Passed on to
\code{\link[isdparser:isd_parse]{isdparser::isd_parse()}}}

\item{parallel}{(logical) do processing in parallel. Default: \code{FALSE}}

\item{cores}{(integer) number of cores to use: Default: 2. We look in
your option "cl.cores", but use default value if not found.}

\item{progress}{(logical) print progress - ignored if \code{parallel=TRUE}.
The default is \code{FALSE} because printing progress adds a small bit of
time, so if processing time is important, then keep as \code{FALSE}}

\item{force}{(logical) force download? Default: \code{FALSE}
We use a cached version (an .rds compressed file) if it exists, but
this will override that behavior.}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
A tibble (data.frame).
}
\description{
Get and parse NOAA ISD/ISH data
}
\details{
\code{isd} saves the full set of weather data for the queried
site locally in the directory specified by the \code{path} argument. You
can access the path for the cached file via \code{attr(x, "source")}

We use \pkg{isdparser} internally to parse ISD files. They are
relatively complex to parse, so a separate package takes care of that.

This function first looks for whether the data for your specific
query has already been downloaded previously in the directory given by
the \code{path} parameter. If not found, the data is requested form NOAA's
FTP server. The first time a dataset is pulled down we must a) download the
data, b) process the data, and c) save a compressed .rds file to disk. The
next time the same data is requested, we only have to read back in the
.rds file, and is quite fast. The benfit of writing to .rds files is that
data is compressed, taking up less space on your disk, and data is read
back in quickly, without changing any data classes in your data, whereas
we'd have to jump through hoops to do that with reading in csv. The
processing can take quite a long time since the data is quite messy and
takes a bunch of regex to split apart text strings. We hope to speed
this process up in the future. See examples below for different behavior.
}
\note{
There are now no transformations (scaling, class changes, etc.)
done on the output data. This may change in the future with parameters
to toggle transformations, but none are done for now. See
\code{\link[isdparser:isd_transform]{isdparser::isd_transform()}} for transformation help.
Comprehensive transformations for all variables are not yet available
but should be available in the next version of this package.

See \link{isd_cache} for managing cached files
}
\section{Errors}{

Note that when you get an error similar to \verb{Error: download failed for https://ftp.ncdc.noaa.gov/pub/data/noaa/1955/011490-99999-1955.gz}, the
file does not exist on NOAA's servers. If your internet is down,
you'll get a different error.
}

\examples{
\dontrun{
# Get station table
(stations <- isd_stations())

## plot stations
### remove incomplete cases, those at 0,0
df <- stations[complete.cases(stations$lat, stations$lon), ]
df <- df[df$lat != 0, ]
### make plot
library("leaflet")
leaflet(data = df) \%>\%
  addTiles() \%>\%
  addCircles()

# Get data
(res <- isd(usaf='011490', wban='99999', year=1986))
(res <- isd(usaf='011690', wban='99999', year=1993))
(res <- isd(usaf='109711', wban=99999, year=1970))

# "additional" and "remarks" data sections included by default
# can toggle that parameter to not include those in output, saves time
(res1 <- isd(usaf='011490', wban='99999', year=1986, force = TRUE))
(res2 <- isd(usaf='011490', wban='99999', year=1986, force = TRUE,
  additional = FALSE))

# The first time a dataset is requested takes longer
system.time( isd(usaf='782680', wban='99999', year=2011) )
system.time( isd(usaf='782680', wban='99999', year=2011) )

# Plot data
## get data for multiple stations
res1 <- isd(usaf='011690', wban='99999', year=1993)
res2 <- isd(usaf='782680', wban='99999', year=2011)
res3 <- isd(usaf='008415', wban='99999', year=2016)
res4 <- isd(usaf='109711', wban=99999, year=1970)
## combine data
library(dplyr)
res_all <- bind_rows(res1, res2, res3, res4)
# add date time
library("lubridate")
dd <- sprintf('\%s \%s', as.character(res_all$date), res_all$time)
res_all$date_time <- ymd_hm(dd)
## remove 999's
res_all <- filter(res_all, temperature < 900)

## plot
if (interactive()) {
  library(ggplot2)
  ggplot(res_all, aes(date_time, temperature)) +
    geom_line() +
    facet_wrap(~usaf_station, scales = 'free_x')
}

# print progress
## note: if the file is already on your system, you'll see no progress bar
(res <- isd(usaf='011690', wban='99999', year=1993, progress=TRUE))

# parallelize processing
# (res <- isd(usaf=172007, wban=99999, year=2016, parallel=TRUE))
}
}
\references{
https://ftp.ncdc.noaa.gov/pub/data/noaa/
https://www1.ncdc.noaa.gov/pub/data/noaa
}
\seealso{
Other isd: 
\code{\link{isd_read}()},
\code{\link{isd_stations_search}()},
\code{\link{isd_stations}()}
}
\concept{isd}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ghcnd.R
\name{ghcnd_splitvars}
\alias{ghcnd_splitvars}
\title{Split variables in data returned from \code{ghcnd}}
\usage{
ghcnd_splitvars(x)
}
\arguments{
\item{x}{An object returned from \code{\link[=ghcnd]{ghcnd()}}}
}
\description{
This function is a helper function for \code{\link[=ghcnd_search]{ghcnd_search()}}. It helps
with cleaning up the data returned from \code{\link[=ghcnd]{ghcnd()}}, to get it in a
format that is easier to work with.
}
\note{
See \code{\link[=ghcnd]{ghcnd()}} examples
}
\author{
Scott Chamberlain, Adam Erickson, Elio Campitelli
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ersst.R
\name{ersst}
\alias{ersst}
\title{NOAA Extended Reconstructed Sea Surface Temperature (ERSST) data}
\usage{
ersst(year, month, overwrite = TRUE, version = "v5", ...)
}
\arguments{
\item{year}{(numeric) A year. Must be > 1853. The max value is whatever
the current year is. Required}

\item{month}{A month, character or numeric. If single digit (e.g. 8), we
add a zero in front (e.g., 08). Required}

\item{overwrite}{(logical) To overwrite the path to store files in or not,
Default: \code{TRUE}}

\item{version}{(character) ERSST version. one of "v5" (default) or "v4"}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
An \code{ncdf4} object. See \pkg{ncdf4} for parsing the output
}
\description{
NOAA Extended Reconstructed Sea Surface Temperature (ERSST) data
}
\details{
See \link{ersst_cache} for managing cached files

\code{ersst()} currently defaults to use ERSST v5 - you can set v4 or v5
using the \code{version} parameter

If a request is unsuccesful, the file written to disk is deleted before
the function exits.

If you use this data in your research please cite rnoaa
(\code{citation("rnoaa")}), and cite NOAA ERSST (see citations at link above)
}
\examples{
\dontrun{
# October, 2015
ersst(year = 2015, month = 10)

# May, 2015
ersst(year = 2015, month = 5)
ersst(year = 2015, month = "05")

# February, 1890
ersst(year = 1890, month = 2)

# Process data
library("ncdf4")
res <- ersst(year = 1890, month = 2)
## varibles
names(res$var)
## get a variable
ncdf4::ncvar_get(res, "ssta")
}
}
\references{
https://www.ncdc.noaa.gov/data-access/marineocean-data/extended-reconstructed-sea-surface-temperature-ersst-v5
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tornadoes.R
\name{tornadoes}
\alias{tornadoes}
\title{Get NOAA tornado data.}
\usage{
tornadoes(...)
}
\arguments{
\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET} (optional)}
}
\value{
A Spatial object is returned of class SpatialLinesDataFrame.
}
\description{
This function gets spatial paths of tornadoes from NOAA's National Weather
Service Storm Prediction Center Severe Weather GIS web page.
}
\note{
See \link{torn_cache} for managing cached files
}
\examples{
\dontrun{
shp <- tornadoes()
library('sp')
if (interactive()) {
  # may take 10 sec or so to render
  plot(shp)
}
}
}
\references{
https://www.spc.noaa.gov/gis/svrgis/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arc2.R
\name{arc2}
\alias{arc2}
\title{Arc2 - Africa Rainfall Climatology version 2}
\usage{
arc2(date, box = NULL, ...)
}
\arguments{
\item{date}{(character/date) one or more dates of the form YYYY-MM-DD}

\item{box}{(numeric) vector of length 4, of the form
\verb{xmin, ymin, xmax, ymax}. optional. If not given, no spatial filtering
is done. If given, we use \code{dplyr::filter()} on a combined set of all dates,
then split the output into tibbles by date}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a list of tibbles with columns:
\itemize{
\item date: date (YYYY-MM-DD)
\item lon: longitude
\item lat: latitude
\item precip: precipitation (mm)
}
}
\description{
Arc2 - Africa Rainfall Climatology version 2
}
\note{
See \link{arc2_cache} for managing cached files
}
\section{box parameter}{

The \code{box} parameter filters the arc2 data to a bounding box you supply.
The function that does the cropping to the bounding box is \code{dplyr::filter}.
You can do any filtering you want on your own if you do not supply
\code{box} and then use whatever tools you want to filter the data by
lat/lon, date, precip values.
}

\examples{
\dontrun{
x = arc2(date = "1983-01-01")
arc2(date = "2017-02-14")

# many dates
arc2(date = c("2019-05-27", "2019-05-28"))
arc2(seq(as.Date("2019-04-21"), by = "day", length.out = 5))
## combine outputs 
x <- arc2(seq(as.Date("2019-05-20"), as.Date("2019-05-25"), "days"))
dplyr::bind_rows(x)

# bounding box filter
box <- c(xmin = 9, ymin = 4, xmax = 10, ymax = 5)
arc2(date = "2017-02-14", box = box)
arc2(date = c("2019-05-27", "2019-05-28"), box = box)
arc2(seq(as.Date("2019-05-20"), as.Date("2019-05-25"), "days"), box = box)
}
}
\references{
docs:
https://ftp.cpc.ncep.noaa.gov/fews/fewsdata/africa/arc2/ARC2_readme.txt
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ghcnd.R
\name{ghcnd}
\alias{ghcnd}
\alias{ghcnd_read}
\title{Get all GHCND data from a single weather site}
\usage{
ghcnd(stationid, refresh = FALSE, ...)

ghcnd_read(path, ...)
}
\arguments{
\item{stationid}{(character) A character vector giving the identification of
the weather stations for which the user would like to pull data. To get a full
and current list of stations, the user can use the \code{\link[=ghcnd_stations]{ghcnd_stations()}}
function. To identify stations within a certain radius of a location, the
user can use the \code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}} function.}

\item{refresh}{(logical) If \code{TRUE} force re-download of data.
Default: \code{FALSE}}

\item{...}{In the case of \code{ghcnd()} additional curl options to pass
through to \link[crul:HttpClient]{crul::HttpClient}. In the case of \code{ghcnd_read}
further options passed on to \code{read.csv}}

\item{path}{(character) a path to a file with a \code{.dly} extension - already
downloaded on your computer}
}
\value{
A tibble (data.frame) which contains data pulled from NOAA's FTP
server for the queried weather site. A README file with more information
about the format of this file is available from NOAA
(https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt).
This file is formatted so each line of the file gives the daily weather
observations for a single weather variable for all days of one month of
one year. In addition to measurements, columns are included for certain
flags, which add information on observation sources and quality and are
further explained in NOAA's README file for the data.
}
\description{
This function uses ftp to access the Global Historical Climatology Network
daily weather data from NOAA's FTP server for a single weather site. It
requires the site identification number for that site and will pull the
entire weather dataset for the site.
}
\details{
This function saves the full set of weather data for the queried
site locally in the directory specified by the \code{path} argument.

You can access the path for the cached file via \code{attr(x, "source")}

You can access the last modified time for the cached file via
\code{attr(x, "file_modified")}

Messages are printed to the console about file path and file last modified time
which you can suppress with \code{suppressMessages()}

For those station ids that are not found, we will delete the file locally
so that a bad station id file is not cached. The returned data for a bad
station id will be an empty data.frame and the attributes are empty strings.
}
\note{
See \link{ghcnd_cache} for managing cached files
}
\section{Base URL}{

The base url for data requests can be changed. The allowed urls are:
https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all (default),
https://ncei.noaa.gov/pub/data/ghcn/daily/all

You can set the base url using the \code{RNOAA_GHCND_BASE_URL} environment
variable; see example below.

The reason for this is that sometimes one base url source is temporarily
down, but another base url may work. It doesn't make sense to allow
an arbitrary base URL; open an issue if there is another valid
base URL for GHNCD data that we should add to the allowed set of base urls.
}

\examples{
\dontrun{
# Get data
ghcnd(stationid = "AGE00147704")

stations <- ghcnd_stations()
ghcnd(stations$id[40])

library("dplyr")
ghcnd(stations$id[80300]) \%>\% select(id, element) \%>\% slice(1:3)

# manipulate data
## using built in fxns
dat <- ghcnd(stationid = "AGE00147704")
(alldat <- ghcnd_splitvars(dat))

## using dplyr
library("dplyr")
dat <- ghcnd(stationid = "AGE00147704")
filter(dat, element == "PRCP", year == 1909)

# refresh the cached file
ghcnd(stationid = "AGE00147704", refresh = TRUE)

# Read in a .dly file you've already downloaded
path <- system.file("examples/AGE00147704.dly", package = "rnoaa")
ghcnd_read(path)

# change the base url for data requests
Sys.setenv(RNOAA_GHCND_BASE_URL =
  "https://ncei.noaa.gov/pub/data/ghcn/daily/all")
ghcnd(stations$id[45], verbose = TRUE)
## must be in the allowed set of urls
# Sys.setenv(RNOAA_GHCND_BASE_URL = "https://google.com")
# ghcnd(stations$id[58], verbose = TRUE)
}
}
\seealso{
To generate a weather dataset for a single weather site that has
been cleaned to a tidier weather format, the user should use the
\code{\link[=ghcnd_search]{ghcnd_search()}} function, which calls \code{ghcnd()} and then
processes the output, or \code{\link[=meteo_tidy_ghcnd]{meteo_tidy_ghcnd()}}, which wraps the
\code{\link[=ghcnd_search]{ghcnd_search()}} function to output a tidy dataframe. To pull
GHCND data from multiple monitors, see \code{\link[=meteo_pull_monitors]{meteo_pull_monitors()}}
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com},
Adam Erickson \email{adam.erickson@ubc.ca}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnoaa_options.R
\name{rnoaa_options}
\alias{rnoaa_options}
\title{rnoaa options}
\usage{
rnoaa_options(cache_messages = TRUE)
}
\arguments{
\item{cache_messages}{(logical) whether to emit messages with information
on caching status for function calls that can cache data. default: \code{TRUE}}
}
\description{
rnoaa options
}
\details{
rnoaa package level options; stored in an internal
package environment \code{roenv}
}
\examples{
\dontrun{
rnoaa_options(cache_messages = FALSE)
}
}
\seealso{
\link{rnoaa_caching} for managing cached files
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ghcnd_states.R
\name{ghcnd_states}
\alias{ghcnd_states}
\alias{ghcnd_countries}
\alias{ghcnd_version}
\title{Get meta-data on the GHCND daily data}
\usage{
ghcnd_states(...)

ghcnd_countries(...)

ghcnd_version(...)
}
\arguments{
\item{...}{In the case of \code{ghcnd()} additional curl options to pass
through to \link[crul:HttpClient]{crul::HttpClient}. In the case of \code{ghcnd_read}
further options passed on to \code{read.csv}}
}
\description{
These function allow you to pull the current versions of certain meta-datasets
for the GHCND, including lists of country and state abbreviations used in some
of the weather station IDs and information about the current version of the
data.
}
\details{
Functions:
\itemize{
\item \code{ghcnd_version}: Get current version of GHCND data
\item \code{ghcnd_states}: Get US/Canada state names and 2-letter codes
\item \code{ghcnd_countries}: Get country names and 2-letter codes
}
}
\examples{
\dontrun{
ghcnd_states()
ghcnd_countries()
ghcnd_version()
}
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com},
Adam Erickson \email{adam.erickson@ubc.ca}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/storm_events.R
\name{storm_events}
\alias{storm_events}
\alias{se_data}
\alias{se_files}
\title{NOAA Storm Events data}
\usage{
se_data(year, type, overwrite = TRUE, ...)

se_files(...)
}
\arguments{
\item{year}{(numeric) a four digit year. see output of \code{se_files()}
for available years. required.}

\item{type}{(character) one of details, fatalities, locations, or
legacy. required.}

\item{overwrite}{(logical) To overwrite the path to store files in or not,
Default: \code{TRUE}}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}
(optional)}
}
\value{
A tibble (data.frame)
}
\description{
NOAA Storm Events data
}
\note{
See \link{stormevents_cache} for managing cached files
}
\examples{
\dontrun{
# get list of files and their urls
res <- se_files()
res
tail(res)

# get data
x <- se_data(year = 2013, type = "details")
x

z <- se_data(year = 1988, type = "fatalities")
z

w <- se_data(year = 2003, type = "locations")
w

leg <- se_data(year = 2003, type = "legacy")
leg
}
}
\references{
https://www.ncdc.noaa.gov/stormevents/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gefs_times}
\alias{gefs_times}
\title{This function is defunct.}
\usage{
gefs_times(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnoaa_caching.R
\docType{data}
\name{rnoaa_caching}
\alias{rnoaa_caching}
\alias{isd_cache}
\alias{cpc_cache}
\alias{arc2_cache}
\alias{lcd_cache}
\alias{bsw_cache}
\alias{ersst_cache}
\alias{torn_cache}
\alias{ghcnd_cache}
\alias{stormevents_cache}
\title{rnoaa caching}
\description{
Manage data caches
}
\details{
To get the cache directory for a data source, see the method
\code{x$cache_path_get()}

\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call

Note that cached files will continue to be used until they are deleted.
It's possible to run into problems when changes happen in your R
setup. For example, at least one user reported changing versions
of this package and running into problems because a cached data
file from a previous version of rnoaa did not work with the newer
version of rnoaa. You should occasionally delete all cached files.
}
\section{Useful user functions}{


Assuming x is a \code{HoardClient} class object, e.g., \code{lcd_cache}
\itemize{
\item \code{x$cache_path_get()} get cache path
\item \code{x$cache_path_set()} set cache path
\item \code{x$list()} returns a character vector of full path file names
\item \code{x$files()} returns file objects with metadata
\item \code{x$details()} returns files with details
\item \code{x$delete()} delete specific files
\item \code{x$delete_all()} delete all files, returns nothing
}
}

\section{Caching objects for each data source}{

\itemize{
\item \code{isd()}/\code{isd_stations()}: \code{isd_cache}
\item \code{cpc_prcp()}: \code{cpc_cache}
\item \code{arc2()}: \code{arc2_cache}
\item \code{lcd()}: \code{lcd_cache}
\item \code{bsw()}: \code{bsw_cache}
\item \code{ersst()}: \code{ersst_cache}
\item \code{tornadoes()}: \code{torn_cache}
\item \code{ghcnd()}/\code{ghcnd_search()}: \code{ghcnd_cache}
\item \code{se_data()}/\code{se_files()}: \code{stormevents_cache}
}
}

\seealso{
\code{\link[=rnoaa_options]{rnoaa_options()}} for managing whether you see messages
about cached files when you request data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_ghcnd.R
\name{meteo_pull_monitors}
\alias{meteo_pull_monitors}
\title{Pull GHCND weather data for multiple weather monitors}
\usage{
meteo_pull_monitors(
  monitors,
  keep_flags = FALSE,
  date_min = NULL,
  date_max = NULL,
  var = "all"
)
}
\arguments{
\item{monitors}{A character vector listing the station IDs for all
weather stations the user would like to pull. To get a full and
current list of stations, the user can use the \code{\link[=ghcnd_stations]{ghcnd_stations()}}
function. To identify stations within a certain radius of a location, the
user can use the \code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}} function.}

\item{keep_flags}{TRUE / FALSE for whether the user would like to keep all
the flags for each weather variable. The default is to not keep the
flags (FALSE). See the note below for more information on these flags.}

\item{date_min}{A character string giving the earliest
date of the daily weather time series that the user would
like in the final output. This character string should be formatted as
"yyyy-mm-dd". If not specified, the default is to keep all daily data for
the queried weather site from the earliest available date.}

\item{date_max}{A character string giving the latest
date of the daily weather time series that the user would
like in the final output. This character string should be formatted as
"yyyy-mm-dd". If not specified, the default is to keep all daily data for
the queried weather site through the most current available date.}

\item{var}{A character vector specifying either \code{"all"} (pull all
available weather parameters for the site) or the weather parameters to
keep in the final data (e.g., \code{c("TMAX", "TMIN")} to only keep
maximum and minimum temperature). Example choices for this argument
include:
\itemize{
\item \code{PRCP}: Precipitation, in tenths of millimeters
\item \code{TAVG}: Average temperature, in tenths of degrees Celsius
\item \code{TMAX}: Maximum temperature, in tenths of degrees Celsius
\item \code{TMIN}: Minimum temperature, in tenths of degrees Celsius
}

A full list of possible weather variables is available in NOAA's README
file for the GHCND data
(https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt).
Most weather stations will only have a small subset of all the possible
weather variables, so the data generated by this function may not include
all of the variables the user specifies through this argument.}
}
\value{
A data frame of daily weather data for multiple weather monitors,
converted to a tidy format. All weather variables may not exist for all
weather stations. Examples of variables returned are:
\itemize{
\item \code{id}: Character string with the weather station site id
\item \code{date}: Date of the observation
\item \code{prcp}: Precipitation, in tenths of mm
\item \code{tavg}: Average temperature, in tenths of degrees Celsius
\item \code{tmax}: Maximum temperature, in tenths of degrees Celsius
\item \code{tmin}: Minimum temperature, in tenths of degrees Celsius
\item \code{awnd}: Average daily wind speed, in meters / second
\item \code{wsfg}: Peak gust wind speed, in meters / second
}

There are other possible weather variables in the Global Historical
Climatology Network; see
http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt for a full
list. If the \code{var} argument is something other than "all", then
only variables included in that argument will be included in the output
data frame. All variables are in the units specified in the linked file
(note that, in many cases, measurements are given in tenths of the units
more often used, e.g., tenths of degrees for temperature). All column names
correspond to variable names in the linked file, but with  all uppercase
letters changed to lowercase.
}
\description{
This function takes a vector of one or more weather station IDs. It will pull
the weather data from the Global Historical Climatology Network's daily
data (GHCND) for each of the stations and join them together in a single tidy
dataframe. For any weather stations that the user calls that are not
available by ftp from GHCND, the function will return a warning
giving the station ID.
}
\note{
The weather flags, which are kept by specifying
\code{keep_flags = TRUE} are:
\itemize{
\item \verb{*_mflag}: Measurement flag, which gives some information on how
the observation was measured.
\item \verb{*_qflag}: Quality flag, which gives quality information on the
measurement, like if it failed to pass certain quality checks.
\item \verb{*_sflag}: Source flag. This gives some information on the
weather collection system (e.g., U.S. Cooperative Summary of the Day,
Australian Bureau of Meteorology) the weather observation comes from.
}

More information on the interpretation of these flags can be found in the
README file for the NCDC's Daily Global Historical Climatology Network's
data at http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt

This function converts any value of -9999 to a missing value for the
variables "prcp", "tmax", "tmin", "tavg", "snow", and "snwd". However,
for some weather observations, there still may be missing values coded
using a series of "9"s of some length. You will want to check your final
data to see if there are lurking missing values given with series of "9"s.

This function may take a while to run.
}
\examples{
\dontrun{

monitors <- c("ASN00003003", "ASM00094299", "ASM00094995", "ASM00094998")
all_monitors_clean <- meteo_pull_monitors(monitors)

}

}
\references{
For more information about the data pulled with this function, see:

Menne, M.J., I. Durre, R.S. Vose, B.E. Gleason, and T.G. Houston, 2012:
An overview of the Global Historical Climatology Network-Daily Database.
Journal of Atmospheric and Oceanic Technology, 29, 897-910,
doi:10.1175/JTECH-D-11-00103.1.
}
\author{
Brooke Anderson \email{brooke.anderson@colostate.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coops.R
\name{coops}
\alias{coops}
\alias{coops_search}
\title{Get NOAA co-ops data}
\usage{
coops_search(
  begin_date = NULL,
  end_date = NULL,
  station_name = NULL,
  product,
  datum = NULL,
  units = "metric",
  time_zone = "gmt",
  application = "rnoaa",
  ...
)
}
\arguments{
\item{begin_date}{(numeric) Date in yyyymmdd format. Required}

\item{end_date}{(numeric) Date in yyyymmdd format. Required}

\item{station_name}{(numeric) a station name. Required}

\item{product}{(character) Specify the data type. See below for Details.
Required}

\item{datum}{(character) See below for Details. Required for all water
level products.}

\item{units}{(character) Specify metric or english (imperial) units,
one of 'metric', 'english'.}

\item{time_zone}{(character) Time zone, one of 'gmt', 'lst', 'lst_ldt'.
For GMT, we convert time stamps to GMT. For LST, we look up the time zone
of the station with its lat/lon values, and assign that time zone. When
\code{product="predictions"} we don't adjust times at all.}

\item{application}{(character) If called within an external package, set
to the name of your organization. Optional}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}
Optional}
}
\value{
List, of length one or two.
\itemize{
\item metadata A list of metadata with slots id, name, lat, lon
\item data A data.frame with data
}
}
\description{
Get NOAA co-ops data
}
\details{
Options for the product paramter. One of:
\itemize{
\item water_level - Preliminary or verified water levels, depending on
availability
\item air_temperature - Air temperature as measured at the station
\item water_temperature - Water temperature as measured at the station
\item wind - Wind speed, direction, and gusts as measured at the station
\item air_pressure - Barometric pressure as measured at the station
\item air_gap - Air Gap (distance between a bridge and the water's surface)
at the station
\item conductivity - The water's conductivity as measured at the station
\item visibility - Visibility from the station's visibility sensor. A
measure of atmospheric clarity
\item humidity - Relative humidity as measured at the station
\item salinity - Salinity and specific gravity data for the station
\item one_minute_water_level - One minute water level data for the station
\item predictions - 6 minute predictions water level data for the station
\item hourly_height - Verified hourly height water level data for
the station
\item high_low - Verified high/low water level data for the station
\item daily_mean - Verified daily mean water level data for the station
\item monthly_mean - Verified monthly mean water level data for the station
\item datums - datums data for the stations
\item currents - Currents data for currents stations
}

Maximum Durations in a Single Call:
\itemize{
\item Products water_level through predictions allow requests for up to
\item Products hourly_height and high_low allow requests for up to
\item Products daily_mean and monthly_mean allow requests for up to
}

Options for the datum parameter. One of:
\itemize{
\item MHHW - Mean higher high water
\item MHW - Mean high water
\item MTL - Mean tide level
\item MSL - Mean sea level
\item MLW - Mean low water
\item MLLW - Mean lower low water
\item NAVD - North American Vertical Datum
\item STND - Station datum
}
}
\examples{
\dontrun{
# Get monthly mean sea level data at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20120301,
  end_date = 20141001, datum = "stnd", product = "monthly_mean")

# Get verified water level data at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20140927,
  end_date = 20140928, datum = "stnd", product = "water_level")

# Get daily mean water level data at Fairport, OH (9063053)
coops_search(station_name = 9063053, begin_date = 20150927,
  end_date = 20150928, product = "daily_mean", datum = "stnd",
  time_zone = "lst")

# Get air temperature at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20140927,
  end_date = 20140928, product = "air_temperature")

# Get water temperature at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20140927,
  end_date = 20140928, product = "water_temperature")

# Get air pressure at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20140927,
  end_date = 20140928, product = "air_pressure")

# Get wind at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20140927,
  end_date = 20140928, product = "wind")

# Get hourly water level height at Key West (8724580)
coops_search(station_name = 8724580, begin_date = 20140927,
  end_date = 20140928, product = "hourly_height", datum = "stnd")

# Get high-low water level at Key West (8724580)
coops_search(station_name = 8724580, begin_date = 20140927,
  end_date = 20140928, product = "high_low", datum = "stnd")

# Get currents data at Pascagoula Harbor (ps0401)
coops_search(station_name = "ps0401", begin_date = 20151221,
  end_date = 20151222, product = "currents")

# Get one-minute water level at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20140927,
  end_date = 20140928, datum = "stnd", product = "one_minute_water_level")

# Get datums at Fort Myers, FL (8725520)
coops_search(station_name = 8725520, product = "datums")

# Get water level predictions at Vaca Key (8723970)
coops_search(station_name = 8723970, begin_date = 20140928,
  end_date = 20140929, datum = "stnd", product = "predictions")

}
}
\references{
https://tidesandcurrents.noaa.gov/api/
https://tidesandcurrents.noaa.gov/map/
}
\author{
Scott Chamberlain, Joseph Stachelek, Tom Philippi
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/storms.R, R/storms_meta.R
\name{storm_data}
\alias{storm_data}
\alias{storm_meta}
\title{Get NOAA wind storm tabular data, metadata, or shp files from IBTrACS}
\usage{
storm_data(...)

storm_meta(...)
}
\description{
Get NOAA wind storm tabular data, metadata, or shp files from IBTrACS

storm_meta
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{is.ncdc_data}
\alias{is.ncdc_data}
\alias{is.ncdc_datasets}
\alias{is.ncdc_datatypes}
\alias{is.ncdc_datacats}
\alias{is.ncdc_locs}
\alias{is.ncdc_locs_cats}
\alias{is.ncdc_stations}
\title{Check object class}
\usage{
is.ncdc_data(x)

is.ncdc_datasets(x)

is.ncdc_datatypes(x)

is.ncdc_datacats(x)

is.ncdc_locs(x)

is.ncdc_locs_cats(x)

is.ncdc_stations(x)
}
\arguments{
\item{x}{input}
}
\description{
Check if an object is of class ncdc_data, ncdc_datasets,
ncdc_datatypes, ncdc_datacats, ncdc_locs, ncdc_locs_cats,
or ncdc_stations
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ghcnd_search.R
\name{ghcnd_search}
\alias{ghcnd_search}
\title{Get a cleaned version of GHCND data from a single weather site}
\usage{
ghcnd_search(
  stationid,
  date_min = NULL,
  date_max = NULL,
  var = "all",
  refresh = FALSE,
  ...
)
}
\arguments{
\item{stationid}{(character) A character vector giving the identification of
the weather stations for which the user would like to pull data. To get a full
and current list of stations, the user can use the \code{\link[=ghcnd_stations]{ghcnd_stations()}}
function. To identify stations within a certain radius of a location, the
user can use the \code{\link[=meteo_nearby_stations]{meteo_nearby_stations()}} function.}

\item{date_min}{A character string giving the earliest
date of the daily weather time series that the user would
like in the final output. This character string should be formatted as
"yyyy-mm-dd". If not specified, the default is to keep all daily data for
the queried weather site from the earliest available date.}

\item{date_max}{A character string giving the latest
date of the daily weather time series that the user would
like in the final output. This character string should be formatted as
"yyyy-mm-dd". If not specified, the default is to keep all daily data for
the queried weather site through the most current available date.}

\item{var}{A character vector specifying either \code{"all"} (pull all
available weather parameters for the site) or the weather parameters to
keep in the final data (e.g., \code{c("TMAX", "TMIN")} to only keep
maximum and minimum temperature). Example choices for this argument
include:
\itemize{
\item \code{PRCP}: Precipitation, in tenths of millimeters
\item \code{TAVG}: Average temperature, in tenths of degrees Celsius
\item \code{TMAX}: Maximum temperature, in tenths of degrees Celsius
\item \code{TMIN}: Minimum temperature, in tenths of degrees Celsius
}

A full list of possible weather variables is available in NOAA's README
file for the GHCND data
(https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt).
Most weather stations will only have a small subset of all the possible
weather variables, so the data generated by this function may not include
all of the variables the user specifies through this argument.}

\item{refresh}{(logical) If \code{TRUE} force re-download of data.
Default: \code{FALSE}}

\item{...}{In the case of \code{ghcnd()} additional curl options to pass
through to \link[crul:HttpClient]{crul::HttpClient}. In the case of \code{ghcnd_read}
further options passed on to \code{read.csv}}
}
\value{
A list object with slots for each of the available specified
weather variables. Each element in the list is a separate time series
dataframe with daily observations, as well as flag values, for one of
the weather variables. The flag values give information on the quality
and source of each observation; see the NOAA README file linked above
for more information. Each data.frame is sorted by date, with the
earliest date first.
}
\description{
This function uses ftp to access the Global Historical Climatology Network
daily weather data from NOAA's FTP server for a single weather monitor site.
It requires the site identification number for that site and will pull the
entire weather dataset for the site. It will then clean this data to convert
it to a tidier format and will also, if requested, filter it to a certain
date range and to certain weather variables.
}
\details{
Messages are printed to the console about file path, file last modified time
which you can suppress with \code{suppressMessages()}
}
\note{
This function calls \code{\link[=ghcnd]{ghcnd()}}, which will download and save
data from all available dates and weather variables for the queried
weather station. The step of limiting the dataset to only certain dates
and / or weather variables, using the \code{date_min}, \code{date_max},
and \code{var} arguments, does not occur until after the full data has
been pulled.
}
\examples{
\dontrun{
# Search based on variable and/or date
ghcnd_search("AGE00147704", var = "PRCP")
ghcnd_search("AGE00147704", var = "PRCP", date_min = "1920-01-01")
ghcnd_search("AGE00147704", var = "PRCP", date_max = "1915-01-01")
ghcnd_search("AGE00147704", var = "PRCP", date_min = "1920-01-01",
             date_max = "1925-01-01")
ghcnd_search("AGE00147704", date_min = "1920-01-01", date_max = "1925-01-01")
ghcnd_search("AGE00147704", var = c("PRCP","TMIN"))
ghcnd_search("AGE00147704", var = c("PRCP","TMIN"), date_min = "1920-01-01")
ghcnd_search("AGE00147704", var = "adfdf")

# refresh the cached file
ghcnd_search("AGE00147704", var = "PRCP", refresh = TRUE)
}
}
\seealso{
\code{\link[=meteo_pull_monitors]{meteo_pull_monitors()}}, \code{\link[=meteo_tidy_ghcnd]{meteo_tidy_ghcnd()}}
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com},
Adam Erickson \email{adam.erickson@ubc.ca}
}
