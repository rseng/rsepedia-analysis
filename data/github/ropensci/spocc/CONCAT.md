

# spocc (SPecies OCCurrence) <img src="man/figures/logo.png" align="right" alt="" width="120">

[![R-check](https://github.com/ropensci/spocc/workflows/R-check/badge.svg)](https://github.com/ropensci/spocc/actions?query=workflow%3AR-check)
[![test-sp-sf](https://github.com/ropensci/spocc/workflows/test-sp-sf/badge.svg)](https://github.com/ropensci/spocc/actions?query=workflow%3Atest-sp-sf)
[![codecov.io](https://codecov.io/github/ropensci/spocc/coverage.svg?branch=master)](https://codecov.io/github/ropensci/spocc?branch=master)
[![cran checks](https://cranchecks.info/badges/worst/spocc)](https://cranchecks.info/pkgs/spocc)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/spocc?color=FAB657)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/spocc)](https://cran.r-project.org/package=spocc)

Docs: <https://docs.ropensci.org/spocc/>

At rOpenSci, we have been writing R packages to interact with many sources of species occurrence data, including [GBIF][gbif], [Vertnet][vertnet], [BISON][bison], [iNaturalist][inat], and [eBird][ebird]. Other databases are out there as well, which we can pull in. `spocc` is an R package to query and collect species occurrence data from many sources. The goal is to to create a seamless search experience across data sources, as well as creating unified outputs across data sources.

`spocc` currently interfaces with nine major biodiversity repositories

1. [Global Biodiversity Information Facility (GBIF)][gbif] (via `rgbif`)
GBIF is a government funded open data repository with several partner organizations with the express goal of providing access to data on Earth's biodiversity. The data are made available by a network of member nodes, coordinating information from various participant organizations and government agencies.

2. [iNaturalist][inat]
iNaturalist provides access to crowd sourced citizen science data on species observations.

3. [VertNet][vertnet] (via `rvertnet`)
Similar to `rgbif` and `rbison` (see below), VertNet provides access to more than 80 million vertebrate records spanning a large number of institutions and museums primarly covering four major disciplines (mammology, herpetology, ornithology, and icthyology).

4. [Biodiversity Information Serving Our Nation][bison] (via `rbison`)
Built by the US Geological Survey's core science analytic team, BISON is a portal that provides access to species occurrence data from several participating institutions.

5. [eBird][ebird] (via `rebird`)
ebird is a database developed and maintained by the Cornell Lab of Ornithology and the National Audubon Society. It provides real-time access to checklist data, data on bird abundance and distribution, and communtiy reports from birders.

6. [iDigBio][idigbio] (via `ridigbio`)
iDigBio facilitates the digitization of biological and paleobiological specimens and their associated data, and houses specimen data, as well as providing their specimen data via RESTful web services.

7. [OBIS][obis]
OBIS (Ocean Biogeographic Information System) allows users to search marine species datasets from all of the world's oceans.

8. [Atlas of Living Australia][ala]
ALA (Atlas of Living Australia) contains information on all the known species in Australia aggregated from a wide range of data providers: museums, herbaria, community groups, government departments, individuals and universities; it contains more than 50 million occurrence records.

The inspiration for this comes from users requesting a more seamless experience across data sources, and from our work on a similar package for taxonomy data ([taxize][taxize]).

__BEWARE:__ In cases where you request data from multiple providers, especially when including GBIF, there could be duplicate records since many providers' data eventually ends up with GBIF. See `?spocc_duplicates`, after installation, for more.

## Learn more

spocc documentation: <https://docs.ropensci.org/spocc/>

## Contributing

See [CONTRIBUTING.md](https://github.com/ropensci/spocc/blob/master/.github/CONTRIBUTING.md)

## Installation

Stable version from CRAN


```r
install.packages("spocc", dependencies = TRUE)
```

Or the development version from GitHub


```r
install.packages("remotes")
remotes::install_github("ropensci/spocc")
```


```r
library("spocc")
```

## Clean data

All data cleaning functionality is in a new package [scrubr](https://github.com/ropensci/scrubr). `scrubr` [on CRAN](https://cran.r-project.org/package=scrubr).

## Make maps

All mapping functionality is now in a separate package [mapr](https://github.com/ropensci/mapr) (formerly known as `spoccutils`), to make `spocc` easier to maintain. `mapr` [on CRAN](https://cran.r-project.org/package=mapr).

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/spocc/issues).
* License: MIT
* Get citation information for `spocc` in R doing `citation(package = 'spocc')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
* Sticker: Images come from Phylopic <http://phylopic.org/>


[gbif]: https://www.gbif.org/
[vertnet]: https://github.com/ropensci/rvertnet
[bison]: https://bison.usgs.gov/
[inat]: https://www.inaturalist.org/
[taxize]: https://github.com/ropensci/taxize
[idigbio]: https://www.idigbio.org/
[obis]: https://obis.org/
[ebird]: https://ebird.org/home
[ala]: https://www.ala.org.au/
spocc 1.2.0
===========

### DEFUNCT

* `as.ecoengine()` is now defunct. `ecoengine` R pkg wasn't used, we used internal code, but since `ecoengine` R is archived since Feb 2020, probably best to remove here (#239)

### MINOR IMPROVEMENTS

* using `wellknown` now instead of `wicket` for Well-known text manipulation (#235)

### BUG FIXES

* fix failing checks (#240)
* bug fixed in ecoengine bbox - but now ecoengine removed, see above (#234)


spocc 1.1.0
===========

### DEFUNCT

* `fixnames()` is now defunct. it was deprecated in a previous version. See `scrubr::fix_names()` (#231)

### NEW FEATURES

* `occ()` can now handle sf objects passed to `geometry`. `spocc` itself does not import/suggest sf, but uses some code donated by Michael Sumner to pull out well known text (WKT) needed for spatially defined queries

### MINOR IMPROVEMENTS

* refactor `occ()`: factor out functions already defined inside of `occ`, add assertions for user parameters (#228)
* package logo/sticker added (#188)

### BUG FIXES

* Fix for ALA data source in `occ()`: total records found count was always 0 because ALA changed the records found field to `totalRecords`
* Fix for Vertnet data source in `occ()`: was using an old parameter `query` passed to `rvertnet::searchbyterm()` - changed to `scientificname` instead


spocc 1.0.8
===========

### BUG FIXES

* fix tests that are failing on cran checks (#225)
* fix ecoengine data option: an `if` statement was failing because we were trying to access an element of a list that is not there sometimes, leading to `NULL` which caused the `if` statement to fail


spocc 1.0.2
===========

### BUG FIXES

* preserve exact bytes for some tests that are failing on cran checks, taxize integration tests, and identifier based search tests (#221)


spocc 1.0.0
===========

### NEW FEATURES

* `source = "inat"` can now return photos. do a query as normal for inat data, and index to the `photos` slot in the data.frame, that will give a nested list of data.frames for each record with links/metadata for the photos (#214) (#217)

### MINOR IMPROVEMENTS

* tests now using vcr (#209)
* add notes in section `iNaturalist notes` in the `?occ` manual page for help with iNaturalist pagination and rate limiting (#215)
* vignette title fix for docs pages (#200)
* `tibble::data_frame`/`tibble::as_data_frame` replaced by `tibble::tibble`/`tibble::as_tibble` throughout package

### DEPRECATED

* `fixnames()` is now deprecated; still useable here until the next version released; please move to using `scrubr::fix_names` (#196)

### BUG FIXES

* fix inat data source for `occ()` queries: change from http to https for the inat base url (#213)
* inat fixes: use rbind fill approach for combining rows of data to fill missing columns safely; work with newer version of their API; unnest lat/lon results into tidy column in resulting data.frame (#215)
* OBIS API changed; changed internals for OBIS data in line with the new API; note that pagination for OBIS has changed (see `?occ` for details); `as.obis.numeric` is gone and replaced with `as.obis.character` (#218)
* fix to `fixnames()`, coerce taxon name to character in case the name is factor class (#211)


spocc 0.9.0
===========

### NEW FEATURES

* `occ()` now attempts to collect errors from requests that fail and puts these error messages (character strings) in the `$meta$errors` spot. We can not always collect errors, and some data providers do not error well: they do not provide a meaningful error message other than that there was an error. (#189) (#207)
* `occ()` gains new parameter `throw_warnings` (logical). By default set to `TRUE` (matches previous behavior) and throws warnings about errors that occur and when no results found for a query. We now prefix each warning with the data provider so you can match up an error or warning with a data provider and (hopefully) query. If set to `FALSE`, warnings are suppressed  (#189) (#207)

### DEFUNCT

* AntWeb has been removed from `spocc`. The AntWeb API has been down for a while, and no response from maintainers (#202) (#203)

### BUG FIXES

* fixes to ebird internals - a new version of `rebird` on CRAN requires a few changes in parameters used. Importantly, ebird now wants species codes instead of full scientific names, but we internally attempt to handle this, so users still just pass scientific names (#205)

### DOCUMENTATION

* make pkgdown docs better: organize functions into meaningful sets (#193) (#197) (#199)
* in the `?spocc_duplicates` manual file for duplicate records, refer to `scrubr` and `CoordinateCleaner` packages (#198)
* in `inspect()` manual file, clarify what the function does (#194)
* now we document better when we use one or the other function for BISON data source (#204)
* `occ()` gains a `return` block with detail about what's returned from the function (#208)



spocc 0.8.0
===========

### NEW FEATURES

* `occ()` gains new parameter `date` to do date range based searches across data sources without having to know the vagaries of each data source (#181)

### BUG FIXES

* fix to idigbio geometry queries (#180)
* fix `print.occdatind` so that empty data.frame's don't throw tibble warnings (#184)
* fix to internal method for standardizing dates `stand_dates()` due to ALA giving back a timestamp now (#182) (#185)
* vertnet fixes (#179)
* fix to geometry bounding box queries (#187) thanks @timcdlucas
* fix to output of list names for gbif data source when using taxonomic IDs, was resulting in booleans, should be the taxonomic IDs  (#191)


spocc 0.7.0
===========

### NEW FEATURES

* Removed javascript and V8 package import and using
`wicket` C++ based package instead. So you no longer need V8
which should make installation easier on some platforms. (#172)

### MINOR IMPROVEMENTS

* `httr` replaced with `crul` for HTTP reqeusts (#174)
* Moved to using markdown for docs. The only thing you should 
notice that's different now is doing curl options is slightly
different - it's just `curl::curl_options()` (#176)
* All `as.*()` functions can now pass on curl options to the 
http client (#177)
* Bumped minimum versions for a number of dependencies

### BUG FIXES

* Fix to `foo_ala()` - the internal plugin for `occ()` that 
handles ALA queries: changed query from full text query using 
`q=foo bar` to `q=taxon_name="foo bar"` - in addition, improved
error handling as sometimes `occurrences` slot is returned in 
results but is empty, whereas before it seemd to always be 
absent if no results (#178)


spocc 0.6.0
===========

### NEW FEATURES

* Added a new data source: Atlas of Living Australia (ALA), under
the abbreviation `ala` (#98)
* Added a new data source: Ocean Biogeographic Information System (OBIS), 
under the abbreviation `obis` (#155)

### MINOR IMPROVEMENTS

* Added note to docs and minor tweak to internal methods to account
for max results from iDigBio of 100,000. Now when you request more than 
100K, you should get a warning saying as much (#169)

### BUG FIXES

* Made `occ2df()` more robust to varied inputs - allowing for users
that may on purpose or not have a subset of the data source slots
normally in the `occdat` class object (#171)


spocc 0.5.4
===========

### MINOR IMPROVEMENTS

* `rvertnet`, a dependency dealing with data from Vertnet, was failing 
on certain searches. `rvertnet` was fixed and a new version on CRAN now. 
No changes here other than requiring the new version of `rvertnet` (#168)
* Fix to internal INAT parsers to handle JSON data output instead of 
CSV output. And fix to internal date parsing; INAT changed field for date 
from `datetime` to `observed_on`.
* Move all `is()` to `inherits()`, and namespace all `setNames()` calls
* We are now using `rgbif::occ_data()` instead of `rgbif::occ_search()`
* We are now using `rvertnet::searchbyterm()` instead of 
`rgbif::vertsearch()`

### BUG FIXES

* Fixes to iDigBio internal plugin - we were dropping scientificname
if geometry was passed by the user. Fixed now. (#167)
* Fixed bug in GBIF internal plugin - when more than 1 result given back
(e.g., multiple searches were done, resulting in a list of objects)
we weren't parsing the output correctly. Fixed now. (#166)



spocc 0.5.0
===========

### NEW FEATURES

* `occ()` now allows queries that only pass `from` and one of the data
source opts params (e.g., `gbifopts`) - allows specifying any options
passed down to the internal functions used to do data queries without
having to use the other params in `occ` (#163)

### MINOR IMPROVEMENTS

* Now using `tibble` for representing data.frames (#164)
* Now using explicit `encoding="UTF-8"` in `httr::content()` calls 
to parse raw data from web requests (#160)
* Now using `ridigbio` as its on CRAN - was using 
internal fxns prior to this (#154)

### BUG FIXES

* There was a problem in the ebird parser where it wasn't processing 
results from ebird with no data. A problem with `has_coords` also 
fixed. (#161)


spocc 0.4.5
===========

### MINOR IMPROVEMENTS

* Using `data.table::setDF()` instead of `data.frame()` to set a `data.table` 
style table to a `data.frame`
* Added many more tests to make it less likely errors will occur
* Added `vertnet` as an option to `occ_options()` to get the options for passing
to `vertopts` in `occ()`

### BUG FIXES

* Fix to `print.occdatind()` - which in last version introduced a bug in this
print method - wasn't fatal as only applied to empty slots in the output 
of a call to `occ()`, but nonetheless, not good (#159)


spocc 0.4.4
===========

### MINOR IMPROVEMENTS

* New import `data.table` for fast list to data.frame

### BUG FIXES

* Fix to ecoengine spatial search - internally we were not making the 
bounding box correctly - fixed now (#158)


spocc 0.4.0
===========

### NEW FEATURES

* New function `as.vertnet()` to coerce various inputs (e.g., result from `occ()`, `occ2df()`, or a key itself) to occurrence data objects (#142)
* `occ()` gains two parameters `start` and `page` to facilitate paging 
through results across data sources, instead of having to page 
individually for each data source. Some sources use the `start` parameter, 
while others use the `page` parameter. See __Paging__ section in `?occ` for
details on Paging (#140)

### MINOR IMPROVEMENTS

* Added Code of Conduct

### BUG FIXES

* `wkt_vis()` now works with WKT polygons with multipe polygons, e.g., 
`spocc::wkt_vis("POLYGON((-125 38.4, -121.8 38.4, -121.8 40.9, -125 40.9, -125 38.4), (-115 22.4, -111.8 22.4, -111.8 30.9, -115 30.9, -115 22.4))")` (#147)
* Fix to `print.occdatind()` to print more helpful info when a 
geometry search is used as opposed to a taxonomy based search (#149)
* Fix to `print.occdatind()` to not fail when first element not present; 
proceeds to next slot with data (#143)
* Fixed problem where `occ()` failed when multiple `geometry` elements
passed in along with taxonomic names (#146)
* Fix to `occ2df()` for combining outputs to not fail when AntWeb 
doesn't give back dates (#144) (#145) - thanks @timcdlucas
* Fix to `occ2df()` to not fail when date field missing (#141)


spocc 0.3.2
===========

### NEW FEATURES

* Added iDigBio as a new data source in `spocc` (#136) (#124)

### MINOR IMPROVEMENTS

* Added much more detail on what parameters in child packages are being used inside of the `occ()` function. Each data source is taken care of in a separate package or set of wrapper functions, and the man file now details what API parameters are being queried (#138)

### BUG FIXES

* Fixed bug where when latitude/longitude columns missing, caused problems downstream in printing outputs, etc. Now we put in NA's when those columns missing (#139)
* Fixed bug in inat data source - `Datetime` variable changed to `datetime`
* Fixed bug in vertnet data source - `occurrenceID` variable changed to `occurrenceid`

spocc 0.3.0
===========

### NEW FEATURES

* Mapping functions all gone, and put into a new package `spoccutils` (#132)
* `occ()` gains new parameter `has_coords` - a global parameter (except for ebird and bison) to return only records with lat/long data. (#128)
* `type` (#134) and `rank` (#133) parameters dropped from `occ()` 
* When object returned by `occ()` is printed, we now include a message that total count of records found (not returned) is not completely known __if ebird is included__, because eBird does not include data on records found on their servers with requests to their API (#111)
* New functions `as.*()` (e.g., `as.gbif`) for most data sources. These functions take in occurrence keys or sets of keys, and retrieve detailed occurrence record data for each key (#112)
* New data source: VertNet (#110)
* `occ2df()` now returns more fields. This function collapses all essential fields that are easy to get in all data sources: `name`, `lat`, `long`, `prov`, `date`, `key`. The `key` field is the occurrence key for each record, which you can use to keep track of individual records, get more data on the record, etc. (#103) (#108)
* New function `inspect()` - takes output from `occ()` or individual occurrence keys and gets detailed occurrence data. 

### MINOR IMPROVEMENTS

* Now importing packages: `jsonlite`, `V8`, `utils`, and `methods`. No longer importing: `ggmap`, `maptools`, `rworldmap`, `sp`, `rgeos`, `RColorBrewer`, `rgdal`, and `leafletR`. Pkgs removed mostly due to splitting off some functionality into `spoccutils`. related issues: (#131) (#132)
* Now importing explicitly all non-base R functions that we use: now importing `methods`, `utils` (#120)
* We now attempt to standardize dates across all data sources, and return that in the output of a call to `occ2df()` (#106)
* `wkt_vis()` now only has an option to view a WKT shape in the browser.

### BUG FIXES

* Fixes to being able to pass curl options on to each data source's functions (#107)

spocc 0.2.4
===========

### MINOR IMPROVEMENTS

* Improved documentation for bounding boxes, their expected format, etc. (#96)
* Remove dependency on the following packages: `assertthat`, `plyr`, `data.table`, and `XML` (#102)
* Using package `gistr` now to post interactive geojson maps on GitHub gists (#100)
* `rgbif` now must be `v0.7.7` or greater (the latest version on CRAN).
* Removed the startup message.

### BUG FIXES

* Duplicate, but not working correctly, function `occ2sp()` removed. The function `occ_to_sp()` function is the working version. (#97)
* Fixed bug where some records returned form GBIF did not have lat/long column headers, and we internally rearranged columns, which caused complete stop when that happened. Fixed now. (#101)
* Changed all `\donttest` to `\dontrun` in examples as requested by CRAN maintainers (#99)

spocc 0.2.2
===========

### NEW FEATURES

* Added new function `occ_names()` to search only for taxonomic names. The goal here is to use ths function if there is some question about what names you want to use to search for occurrences with. (#84). Suggested by @jarioksa
* New function `occ_names_options()` to quickly get parameter options to pass to `occ_names()`.
* New `summary()` method for the `occdat` `S3` object that is output from `occ()` (#83)
* In many places in `spocc` (README, vignette, `occ()` documentation file, at package startup), we make it clear that there could be duplicate records returned in certain scenarios. And a new documentation page detailing what to watch out for: `?spocc_duplicates`. (#77)

### MINOR IMPROVEMENTS

* All latitude/longitude column headers are now changed to latitude and longitude, whereas they use to vary from `latitude`, `decimalLatitude`, `Latitude`, `lat`, and `decimal_latitude`. (#91)
* Default is 500 now for the `limit` parameter in `occ()` (#78)
* You can now pass in `limit` to each functions options parameter, and it will work. Each data source can have a different parameter internally from `limit`, but now internally within `spocc`, we allow you to use `limit` so you don't have to know what the data source specific parameter is. (#81)
* There is a now a startup message to give information on the package (#79)
* `occ_options()` gains new parameter `where` to print either in the console or to open man file in the IDE, or prints to console in command line R. 

spocc 0.2.0
===========

### NEW FEATURES

* `occ()` gains new parameter `callopts` to pass on curl debugging options to `httr::GET()` (#35)
* `wkt_vis()` now by default plots a well known text area (WKT) on an interactive mapbox map in your default browser. New parameter `which` allows you to choose the interactive map or a static ggplot2 map. (#70)
* Individual data sources `occ()` gains new class. In the previous version of this package, a `data.frame` was printed. Now the data is assigned the object `occdatind` (short for _occdat individual_).
* `occ()` now uses a print method for the `occdatind` class, adopted from `dplyr` that prints a brief `data.frame`, with columns wrapped to fit the width of your console, and additional columns not printed given at bottom with their class type. Note that the print behavior for the resulting object of an `occ()` call remains the same. (#69) (#74)

### MINOR IMPROVEMENTS

* Added `whisker` as a package import to use in the `wkt_vis()` function. (#70)
* Mapping functions now all accept the same input. Previously `mapggplot()` accepted the output of `occ()`, of class `occdat`, while the other two functions for mapping, `mapleaflet()` and `mapgist()` accepted a `data.frame`. Now all three functions accept  the output of `occ()`, an object of class `occdat`. (#75)
* The `meta` slot in each returned object (indexed by `object$meta`) contains spots for `returned` and `found`, to designate number of records returned, and number of records found. (#64)

### BUG FIXES

* Fixed bug in AntWeb output, where there was supposed to be a column titled `name`. (#71)

spocc 0.1.4
===========

### NEW FEATURES

* Can now do geometry only queries. See examples in `occ()`.
* In addition, you can pass in sp objects of SpatialPolygons or SpatialPolygonsDataFrame classes.

spocc 0.1.2
===========

### NEW FEATURES

* There were quite a few changes in one of the key packages that `spocc` depends on: `rgbif`. A number of input and output parameter names changed. A new version of `rgbif` was pushed to CRAN. (#56)
* New function `clean_spocc()` started (not finished yet) to attempt to clean data. For example, one use case is removing impossible lat/long values (i.e., longitue values greater than absolute 180). Another, not implemented yet, is to remove points that are not in the country or habitat your points are supposed to be in. (#44)
* New function `fixnames()` to trim species names with optional input parameters to make data easier to use for mapping.
* New function `wkt_vis()` to visualize a WKT (well-known text) area on a map. Uses `ggmap` to pull down a Google map so that the visualization has some geographic and natural earth context. We'll soon introduce an interactive version of this function that will bring up a small Shiny app to draw a WKT area, then return those coordinates to your R session. (#34)

### MINOR IMPROVEMENTS

* Added a CONTRIBUTING.md file to the github repo to help guide contributions (#61)
* Packages that require a certain version are forced to be X version or greater. Thes are rinat (>= 0.1.1), rbison (>= 0.3.2), rgbif (>= 0.6.2), ecoengine (>= 1.3), rebird (>= 0.1.1), AntWeb (>= 0.6.1), and leafletR (>= 0.2-0). This should help avoid problems.
* General improvement to function documentation.

spocc 0.1.0
===========

* Initial release to CRAN
## Test environments

* local OS X install, R 4.0.3 Patched
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (release/devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 6 reverse dependencies. No problems related to this package were found. Summary at <https://github.com/ropensci/spocc/tree/master/revdep>

--------

This version fixes some bugs, makes a function defunct, and depends on new version of vcr that fixes failing cran checks.

Thanks!
Scott Chamberlain
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
# CONTRIBUTING 

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/spocc/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/spocc.git`
* Make sure to track progress upstream (i.e., on our version of `spocc` at `ropensci/spocc`) by doing `git remote add upstream https://github.com/ropensci/spocc.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/spocc`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-06-30 r78761) |
|os       |macOS Catalina 10.15.6                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-07-31                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|spocc   |1.0.8 |1.1.0 |*  |
|dplyr   |NA    |1.0.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*```{r echo=FALSE}
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
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

# spocc (SPecies OCCurrence) <img src="man/figures/logo.png" align="right" alt="" width="120">

[![R-check](https://github.com/ropensci/spocc/workflows/R-check/badge.svg)](https://github.com/ropensci/spocc/actions?query=workflow%3AR-check)
[![test-sp-sf](https://github.com/ropensci/spocc/workflows/test-sp-sf/badge.svg)](https://github.com/ropensci/spocc/actions?query=workflow%3Atest-sp-sf)
[![codecov.io](https://codecov.io/github/ropensci/spocc/coverage.svg?branch=master)](https://codecov.io/github/ropensci/spocc?branch=master)
[![cran checks](https://cranchecks.info/badges/worst/spocc)](https://cranchecks.info/pkgs/spocc)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/spocc?color=FAB657)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/spocc)](https://cran.r-project.org/package=spocc)

Docs: <https://docs.ropensci.org/spocc/>

At rOpenSci, we have been writing R packages to interact with many sources of species occurrence data, including [GBIF][gbif], [Vertnet][vertnet], [BISON][bison], [iNaturalist][inat], and [eBird][ebird]. Other databases are out there as well, which we can pull in. `spocc` is an R package to query and collect species occurrence data from many sources. The goal is to to create a seamless search experience across data sources, as well as creating unified outputs across data sources.

`spocc` currently interfaces with nine major biodiversity repositories

1. [Global Biodiversity Information Facility (GBIF)][gbif] (via `rgbif`)
GBIF is a government funded open data repository with several partner organizations with the express goal of providing access to data on Earth's biodiversity. The data are made available by a network of member nodes, coordinating information from various participant organizations and government agencies.

2. [iNaturalist][inat]
iNaturalist provides access to crowd sourced citizen science data on species observations.

3. [VertNet][vertnet] (via `rvertnet`)
Similar to `rgbif` and `rbison` (see below), VertNet provides access to more than 80 million vertebrate records spanning a large number of institutions and museums primarly covering four major disciplines (mammology, herpetology, ornithology, and icthyology).

4. [Biodiversity Information Serving Our Nation][bison] (via `rbison`)
Built by the US Geological Survey's core science analytic team, BISON is a portal that provides access to species occurrence data from several participating institutions.

5. [eBird][ebird] (via `rebird`)
ebird is a database developed and maintained by the Cornell Lab of Ornithology and the National Audubon Society. It provides real-time access to checklist data, data on bird abundance and distribution, and communtiy reports from birders.

6. [iDigBio][idigbio] (via `ridigbio`)
iDigBio facilitates the digitization of biological and paleobiological specimens and their associated data, and houses specimen data, as well as providing their specimen data via RESTful web services.

7. [OBIS][obis]
OBIS (Ocean Biogeographic Information System) allows users to search marine species datasets from all of the world's oceans.

8. [Atlas of Living Australia][ala]
ALA (Atlas of Living Australia) contains information on all the known species in Australia aggregated from a wide range of data providers: museums, herbaria, community groups, government departments, individuals and universities; it contains more than 50 million occurrence records.

The inspiration for this comes from users requesting a more seamless experience across data sources, and from our work on a similar package for taxonomy data ([taxize][taxize]).

__BEWARE:__ In cases where you request data from multiple providers, especially when including GBIF, there could be duplicate records since many providers' data eventually ends up with GBIF. See `?spocc_duplicates`, after installation, for more.

## Learn more

spocc documentation: <https://docs.ropensci.org/spocc/>

## Contributing

See [CONTRIBUTING.md](https://github.com/ropensci/spocc/blob/master/.github/CONTRIBUTING.md)

## Installation

Stable version from CRAN

```{r eval=FALSE}
install.packages("spocc", dependencies = TRUE)
```

Or the development version from GitHub

```{r eval=FALSE}
install.packages("remotes")
remotes::install_github("ropensci/spocc")
```

```{r}
library("spocc")
```

## Clean data

All data cleaning functionality is in a new package [scrubr](https://github.com/ropensci/scrubr). `scrubr` [on CRAN](https://cran.r-project.org/package=scrubr).

## Make maps

All mapping functionality is now in a separate package [mapr](https://github.com/ropensci/mapr) (formerly known as `spoccutils`), to make `spocc` easier to maintain. `mapr` [on CRAN](https://cran.r-project.org/package=mapr).

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/spocc/issues).
* License: MIT
* Get citation information for `spocc` in R doing `citation(package = 'spocc')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
* Sticker: Images come from Phylopic <http://phylopic.org/>


[gbif]: https://www.gbif.org/
[vertnet]: https://github.com/ropensci/rvertnet
[bison]: https://bison.usgs.gov/
[inat]: https://www.inaturalist.org/
[taxize]: https://github.com/ropensci/taxize
[idigbio]: https://www.idigbio.org/
[obis]: https://obis.org/
[ebird]: https://ebird.org/home
[ala]: https://www.ala.org.au/
---
title: spocc introduction
author: Scott Chamberlain
date: "2020-12-18"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{spocc introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



The rOpenSci projects aims to provide programmatic access to scientific data repositories on the web. A vast majority of the packages in our current suite retrieve some form of biodiversity or taxonomic data. Since several of these datasets have been georeferenced, it provides numerous opportunities for visualizing species distributions, building species distribution maps, and for using it analyses such as species distribution models. In an effort to streamline access to these data, we have developed a package called `spocc`, which provides a unified API to all the biodiversity sources that we provide. The obvious advantage is that a user can interact with a common API and not worry about the nuances in syntax that differ between packages. As more data sources come online, users can access even more data without significant changes to their code. However, it is important to note that spocc will never replicate the full functionality that exists within specific packages. Therefore users with a strong interest in one of the specific data sources listed below would benefit from familiarising themselves with the inner working of the appropriate packages.

## Data Sources

`spocc` currently interfaces with nine major biodiversity repositories

1. [Global Biodiversity Information Facility (GBIF)](https://www.gbif.org/) (via `rgbif`)
GBIF is a government funded open data repository with several partner organizations with the express goal of providing access to data on Earth's biodiversity. The data are made available by a network of member nodes, coordinating information from various participant organizations and government agencies.

2. [iNaturalist](https://www.inaturalist.org/)
iNaturalist provides access to crowd sourced citizen science data on species observations.

3. [VertNet](http://vertnet.org/) (via `rvertnet`)
Similar to `rgbif` and `rbison` (see below), VertNet provides access to more than 80 million vertebrate records spanning a large number of institutions and museums primarly covering four major disciplines (mammology, herpetology, ornithology, and icthyology).

4. Biodiversity Information Serving Our Nation (https://bison.usgs.gov/) (via `rbison`)
Built by the US Geological Survey's core science analytic team, BISON is a portal that provides access to species occurrence data from several participating institutions.

5. [eBird](https://ebird.org/home) (via `rebird`)
ebird is a database developed and maintained by the Cornell Lab of Ornithology and the National Audubon Society. It provides real-time access to checklist data, data on bird abundance and distribution, and communtiy reports from birders.

6. [iDigBio](https://www.idigbio.org/) (via `ridigbio`)
iDigBio facilitates the digitization of biological and paleobiological specimens and their associated data, and houses specimen data, as well as providing their specimen data via RESTful web services.

7. [OBIS](https://obis.org/)
OBIS (Ocean Biogeographic Information System) allows users to search marine species datasets from all of the world's oceans.

8. [Atlas of Living Australia](https://www.ala.org.au/)
ALA (Atlas of Living Australia) contains information on all the known species in Australia aggregated from a wide range of data providers: museums, herbaria, community groups, government departments, individuals and universities; it contains more than 50 million occurrence records.

__Important Note:__ It's important to keep in mind that several data providers interface with many of the above mentioned repositories. This means that occurence data obtained from BISON may be duplicates of data that are also available through GBIF. We do not have a way to resolve these duplicates or overlaps at this time but it is an issue we are hoping to resolve in future versions of the package. See `?spocc_duplicates`, after installation, for more.


## Data retrieval

The most significant function in spocc is the `occ` (short for occurrence) function. `occ` takes a query, often a species name, and searches across all data sources specified in the `from` argument. For example, one can search for all occurrences of [Sharp-shinned Hawks](https://www.allaboutbirds.org/guide/sharp-shinned_hawk/id) (_Accipiter striatus_) from the GBIF database with the following R call.


```r
library('spocc')
(df <- occ(query = 'Accipiter striatus', from = 'gbif'))
#> Searched: gbif
#> Occurrences - Found: 1,093,434, Returned: 500
#> Search type: Scientific
#>   gbif: Accipiter striatus (500)
```

The data returned are part of a `S3` class called `occdat`. This class has slots for each of the data sources described above. One can easily switch the source by changing the `from` parameter in the function call above.

Within each data source is the set of species queried. In the above example, we only asked for occurrence data for one species, but we could have asked for any number. Let's say we asked for data for two species: _Accipiter striatus_, and _Pinus contorta_. Then the structure of the response would be

```
response -- |
            | -- gbif ------- |
                              | -- Accipiter_striatus
                              | -- Pinus_contorta

            | -- bison ------ |
                              | -- Accipiter_striatus
                              | -- Pinus_contorta

            ... and so on for each data source

```

If you only request data from gbif, like `from = 'gbif'`, then the other four source slots are present in the response object, but have no data.

You can quickly get just the GBIF data by indexing to it, like


```r
df$gbif
#> Species [Accipiter striatus (500)] 
#> First 10 rows of [Accipiter_striatus]
#> 
#> # A tibble: 500 x 83
#>    name  longitude latitude prov  issues key   scientificName datasetKey
#>    <chr>     <dbl>    <dbl> <chr> <chr>  <chr> <chr>          <chr>     
#>  1 Acci…    -107.      35.1 gbif  cdrou… 2542… Accipiter str… 50c9509d-…
#>  2 Acci…     -90.0     37.1 gbif  cdrou… 2543… Accipiter str… 50c9509d-…
#>  3 Acci…     -99.3     36.5 gbif  cdrou… 2543… Accipiter str… 50c9509d-…
#>  4 Acci…     -76.0     39.6 gbif  cdrou… 2543… Accipiter str… 50c9509d-…
#>  5 Acci…     -73.5     40.7 gbif  gass8… 2543… Accipiter str… 50c9509d-…
#>  6 Acci…    -118.      34.6 gbif  cdrou… 2549… Accipiter str… 50c9509d-…
#>  7 Acci…    -121.      36.6 gbif  cdrou… 2550… Accipiter str… 50c9509d-…
#>  8 Acci…     -97.3     27.6 gbif  cdrou… 2550… Accipiter str… 50c9509d-…
#>  9 Acci…     -88.9     30.5 gbif  cdrou… 2550… Accipiter str… 50c9509d-…
#> 10 Acci…     -96.9     33.1 gbif  cdrou… 2550… Accipiter str… 50c9509d-…
#> # … with 490 more rows, and 75 more variables: publishingOrgKey <chr>,
#> #   installationKey <chr>, publishingCountry <chr>, protocol <chr>,
#> #   lastCrawled <chr>, lastParsed <chr>, crawlId <int>,
#> #   hostingOrganizationKey <chr>, basisOfRecord <chr>, occurrenceStatus <chr>,
#> #   taxonKey <int>, kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   orderKey <int>, familyKey <int>, genusKey <int>, speciesKey <int>,
#> #   acceptedTaxonKey <int>, acceptedScientificName <chr>, kingdom <chr>,
#> #   phylum <chr>, order <chr>, family <chr>, genus <chr>, species <chr>,
#> #   genericName <chr>, specificEpithet <chr>, taxonRank <chr>,
#> #   taxonomicStatus <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <date>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, isInCluster <lgl>,
#> #   geodeticDatum <chr>, class <chr>, countryCode <chr>, country <chr>,
#> #   rightsHolder <chr>, identifier <chr>, `http://unknown.org/nick` <chr>,
#> #   verbatimEventDate <chr>, datasetName <chr>, gbifID <chr>,
#> #   verbatimLocality <chr>, collectionCode <chr>, occurrenceID <chr>,
#> #   taxonID <chr>, catalogNumber <chr>, recordedBy <chr>,
#> #   `http://unknown.org/occurrenceDetails` <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, identifiedBy <chr>, identificationID <chr>,
#> #   informationWithheld <chr>, occurrenceRemarks <chr>,
#> #   identificationRemarks <chr>, infraspecificEpithet <chr>,
#> #   nomenclaturalCode <chr>, locality <chr>, vernacularName <chr>,
#> #   fieldNotes <chr>, verbatimElevation <chr>, behavior <chr>,
#> #   higherClassification <chr>
```

When you get data from multiple providers, the fields returned are slightly different because each data provider uses different formats for their data; different arrangements of data and different variable names for the same thing (e.g., one data provider may call latitude "latitude", while another may call it "lat"). For example:


```r
df <- occ(query = 'Accipiter striatus', from = c('gbif', 'bison'), limit = 25)
df$gbif$data$Accipiter_striatus
#> # A tibble: 25 x 74
#>    name  longitude latitude issues prov  key   scientificName datasetKey
#>    <chr>     <dbl>    <dbl> <chr>  <chr> <chr> <chr>          <chr>     
#>  1 Acci…    -107.      35.1 cdrou… gbif  2542… Accipiter str… 50c9509d-…
#>  2 Acci…     -90.0     37.1 cdrou… gbif  2543… Accipiter str… 50c9509d-…
#>  3 Acci…     -99.3     36.5 cdrou… gbif  2543… Accipiter str… 50c9509d-…
#>  4 Acci…     -76.0     39.6 cdrou… gbif  2543… Accipiter str… 50c9509d-…
#>  5 Acci…     -73.5     40.7 gass8… gbif  2543… Accipiter str… 50c9509d-…
#>  6 Acci…    -118.      34.6 cdrou… gbif  2549… Accipiter str… 50c9509d-…
#>  7 Acci…    -121.      36.6 cdrou… gbif  2550… Accipiter str… 50c9509d-…
#>  8 Acci…     -97.3     27.6 cdrou… gbif  2550… Accipiter str… 50c9509d-…
#>  9 Acci…     -88.9     30.5 cdrou… gbif  2550… Accipiter str… 50c9509d-…
#> 10 Acci…     -96.9     33.1 cdrou… gbif  2550… Accipiter str… 50c9509d-…
#> # … with 15 more rows, and 66 more variables: publishingOrgKey <chr>,
#> #   installationKey <chr>, publishingCountry <chr>, protocol <chr>,
#> #   lastCrawled <chr>, lastParsed <chr>, crawlId <int>,
#> #   hostingOrganizationKey <chr>, basisOfRecord <chr>, occurrenceStatus <chr>,
#> #   taxonKey <int>, kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   orderKey <int>, familyKey <int>, genusKey <int>, speciesKey <int>,
#> #   acceptedTaxonKey <int>, acceptedScientificName <chr>, kingdom <chr>,
#> #   phylum <chr>, order <chr>, family <chr>, genus <chr>, species <chr>,
#> #   genericName <chr>, specificEpithet <chr>, taxonRank <chr>,
#> #   taxonomicStatus <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <date>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, isInCluster <lgl>,
#> #   geodeticDatum <chr>, class <chr>, countryCode <chr>, country <chr>,
#> #   rightsHolder <chr>, identifier <chr>, `http://unknown.org/nick` <chr>,
#> #   verbatimEventDate <chr>, datasetName <chr>, gbifID <chr>,
#> #   verbatimLocality <chr>, collectionCode <chr>, occurrenceID <chr>,
#> #   taxonID <chr>, catalogNumber <chr>, recordedBy <chr>,
#> #   `http://unknown.org/occurrenceDetails` <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, identifiedBy <chr>, identificationID <chr>,
#> #   informationWithheld <chr>, occurrenceRemarks <chr>
df$bison$data$Accipiter_striatus
#> # A tibble: 25 x 35
#>    date       providedScienti…  year countryCode ambiguous latlon
#>    <date>     <chr>            <int> <chr>       <lgl>     <chr> 
#>  1 2001-10-10 Accipiter stria…  2001 US          FALSE     -83.2…
#>  2 1980-09-21 Accipiter stria…  1980 US          FALSE     -75.7…
#>  3 1980-10-04 Accipiter stria…  1980 US          FALSE     -75.7…
#>  4 1966-10-13 Accipiter stria…  1966 US          FALSE     -75.9…
#>  5 1987-10-10 Accipiter stria…  1987 US          FALSE     -74.2…
#>  6 1990-10-03 Accipiter stria…  1990 US          FALSE     -75.7…
#>  7 1994-10-05 Accipiter stria…  1994 US          FALSE     -75.7…
#>  8 1976-09-18 Accipiter stria…  1976 US          FALSE     -75.7…
#>  9 1979-10-07 Accipiter stria…  1979 US          FALSE     -75.7…
#> 10 1980-10-01 Accipiter stria…  1980 US          FALSE     -75.7…
#> # … with 15 more rows, and 29 more variables: computedCountyFips <chr>,
#> #   occurrenceID <chr>, longitude <dbl>, basisOfRecord <chr>,
#> #   providedCommonName <chr>, collectionID <chr>,
#> #   ownerInstitutionCollectionCode <chr>, name <chr>, institutionID <chr>,
#> #   computedStateFips <chr>, license <chr>, TSNs <chr>, providerID <int>,
#> #   stateProvince <chr>, higherGeographyID <chr>, latitude <dbl>, geo <chr>,
#> #   provider <chr>, calculatedCounty <chr>, ITISscientificName <chr>,
#> #   pointPath <chr>, kingdom <chr>, calculatedState <chr>,
#> #   hierarchy_homonym_string <chr>, centroid <chr>, ITIScommonName <chr>,
#> #   resourceID <chr>, ITIStsn <chr>, prov <chr>
```

We provide a function `occ2df` that pulls out a few key columns needed for making maps:


```r
occ2df(df)
#> # A tibble: 50 x 6
#>    name                            longitude latitude prov  date       key      
#>    <chr>                               <dbl>    <dbl> <chr> <date>     <chr>    
#>  1 Accipiter striatus Vieillot, 1…    -107.      35.1 gbif  2020-01-02 25429665…
#>  2 Accipiter striatus Vieillot, 1…     -90.0     37.1 gbif  2020-01-01 25430843…
#>  3 Accipiter striatus Vieillot, 1…     -99.3     36.5 gbif  2020-01-01 25430853…
#>  4 Accipiter striatus Vieillot, 1…     -76.0     39.6 gbif  2020-01-01 25430927…
#>  5 Accipiter striatus Vieillot, 1…     -73.5     40.7 gbif  2020-01-01 25430953…
#>  6 Accipiter striatus Vieillot, 1…    -118.      34.6 gbif  2020-01-03 25499936…
#>  7 Accipiter striatus Vieillot, 1…    -121.      36.6 gbif  2020-01-04 25500018…
#>  8 Accipiter striatus Vieillot, 1…     -97.3     27.6 gbif  2020-01-04 25500046…
#>  9 Accipiter striatus Vieillot, 1…     -88.9     30.5 gbif  2020-01-04 25500173…
#> 10 Accipiter striatus Vieillot, 1…     -96.9     33.1 gbif  2020-01-05 25500177…
#> # … with 40 more rows
```

`occ2df()` not only combines data into a single data.frame, but it also standardizes the key columns (name, longitude, latitude, prov (provider), date, and key (occurrence key)). Note that you can look up the exact occurrence with the data provider using the `key` value.

### Standardized parameters

Each data source has a variety of different ways, or parameters, to use to search its data. Some of the parameters are the same across data sources. In `occ()` we've attempted to surface those similar parameters so you can have a single way to define a parameter and it gets applied to every data source. This way you don't have to know the vagaries of each data source, what formatting they expect, etc.

The standardized parameters in `occ()` are:

- query: a scientific taxon name
- limit: number of records to retrieve
- start: page number to start at
- page: page number to retrieve
- geometry: a "spatial filter" - bounding box, well known text, or an sp or sf polygon or multipolygon
- has_coords: exclude records without latitude/longitude data
- date: a date range

However, not all parameters across data sources are able to be standardized, so you can pass data source specific parameters to their matching parameter name, e.g., pass GBIF parameters to `gbifopts` and ALA parameters to `alaopts`. 


## Clean up taxonomic names

See the vignette [cleaning names](https://docs.ropensci.org/spocc/articles/fixnames)

## Clean data

All data cleaning functionality is in a new package [scrubr](https://github.com/ropensci/scrubr). [On CRAN](https://cran.r-project.org/package=scrubr).

## Make maps

All mapping functionality is now in a separate package [mapr](https://github.com/ropensci/mapr) (formerly known as `spoccutils`), to make `spocc` easier to maintain. [On CRAN](https://cran.r-project.org/package=mapr).
---
title: cleaning names
author: Scott Chamberlain
date: "2020-12-18"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{cleaning names}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



One problem you often run in to is that there can be various names for the same taxon in any one source. For example:


```r
library(spocc)
df <- occ(query = 'Pinus contorta', from = c('gbif', 'bison'), limit = 50)
unique(df$gbif$data$Pinus_contorta$name)
#> [1] "Pinus contorta Douglas ex Loudon"             
#> [2] "Pinus contorta var. contorta"                 
#> [3] "Pinus contorta var. murrayana (Balf.) Engelm."
unique(df$bison$data$Pinus_contorta$name)
#> [1] "Pinus contorta"
```

This is fine, but when trying to make a map in which points are colored for each taxon, you can have many colors for a single taxon, where instead one color per taxon is more appropriate. There is a function in `scrubr` called `fix_names()`, which has a few options in which you can take the shortest names (usually just the plain binomials like _Homo sapiens_), or the original name queried, or a vector of names supplied by the user.


```r
install.packages("scrubr")
```


```r
library(scrubr)
df$gbif$data$Pinus_contorta <- fix_names(df$gbif$data$Pinus_contorta, how = 'shortest')
df$bison$data$Pinus_contorta <- fix_names(df$bison$data$Pinus_contorta, how = 'shortest')
unique(df$gbif$data$Pinus_contorta$name)
#> [1] "Pinus contorta var. contorta"
unique(df$bison$data$Pinus_contorta$name)
#> [1] "Pinus contorta"
df_comb <- occ2df(df)
head(df_comb); tail(df_comb)
#> # A tibble: 6 x 6
#>   name                         longitude latitude prov  date       key       
#>   <chr>                            <dbl>    <dbl> <chr> <date>     <chr>     
#> 1 Pinus contorta var. contorta    -115.      50.9 gbif  2020-01-01 2543085192
#> 2 Pinus contorta var. contorta      17.6     59.8 gbif  2020-01-01 2548826490
#> 3 Pinus contorta var. contorta      19.2     64.0 gbif  2020-01-06 2549045731
#> 4 Pinus contorta var. contorta      19.3     64.0 gbif  2020-01-06 2549053727
#> 5 Pinus contorta var. contorta    -123.      49.3 gbif  2020-01-04 2550016817
#> 6 Pinus contorta var. contorta    -106.      39.8 gbif  2020-01-07 2557738499
#> # A tibble: 6 x 6
#>   name           longitude latitude prov  date       key       
#>   <chr>              <dbl>    <dbl> <chr> <date>     <chr>     
#> 1 Pinus contorta     -115.     45.4 bison 2011-09-19 2097723549
#> 2 Pinus contorta     -115.     45.4 bison 2011-09-19 2097723550
#> 3 Pinus contorta     -115.     45.6 bison 2010-09-04 2097723551
#> 4 Pinus contorta     -116.     45.8 bison 2004-07-27 2097723552
#> 5 Pinus contorta     -116.     45.8 bison 2004-07-27 2097723554
#> 6 Pinus contorta     -116.     45.8 bison 2004-07-27 2097723555
```

Now with one taxon name for each taxon we can more easily make a plot.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.ecoengine.R
\name{as.ecoengine}
\alias{as.ecoengine}
\title{Coerce occurrence keys to ecoenginekey/occkey objects}
\usage{
as.ecoengine(...)
}
\arguments{
\item{...}{ignored}
}
\description{
DEFUNCT
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.gbif.R
\name{as.gbif}
\alias{as.gbif}
\title{Coerce occurrence keys to gbifkey/occkey objects}
\usage{
as.gbif(x, ...)
}
\arguments{
\item{x}{Various inputs, including the output from a call to
\code{\link[=occ]{occ()}} (class occdat), \code{\link[=occ2df]{occ2df()}} (class data.frame),
or a list, numeric, character, gbifkey, or occkey.}

\item{...}{curl options; named parameters passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
One or more in a list of both class gbifkey and occkey
}
\description{
Coerce occurrence keys to gbifkey/occkey objects
}
\details{
Internally, we use \code{\link[rgbif:occ_get]{rgbif::occ_get()}}, whereas
\code{\link[=occ]{occ()}} uses \code{\link[rgbif:occ_data]{rgbif::occ_data()}}. We can use
\code{\link[rgbif:occ_get]{rgbif::occ_get()}} here because we have the occurrence key to
go directly to the occurrence record.
}
\examples{
\dontrun{
spnames <- c('Accipiter striatus', 'Setophaga caerulescens', 
  'Spinus tristis')
out <- occ(query=spnames, from=c('gbif','ebird'), 
  gbifopts=list(hasCoordinate=TRUE), limit=2)
res <- occ2df(out)
(tt <- as.gbif(out))
(uu <- as.gbif(res))
as.gbif(as.numeric(res$key[1]))
as.gbif(res$key[1])
as.gbif(as.list(res$key[1:2]))
as.gbif(tt[[1]])
as.gbif(uu[[1]])
as.gbif(tt[1:2])
}
}
\seealso{
Other coercion: 
\code{\link{as.ala}()},
\code{\link{as.bison}()},
\code{\link{as.idigbio}()},
\code{\link{as.inat}()},
\code{\link{as.obis}()},
\code{\link{as.vertnet}()}
}
\concept{coercion}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_names.R
\name{occ_names}
\alias{occ_names}
\title{Search for species names across many data sources.}
\usage{
occ_names(
  query = NULL,
  from = "gbif",
  limit = 100,
  rank = "species",
  callopts = list(),
  gbifopts = list(),
  bisonopts = list()
)
}
\arguments{
\item{query}{(character) One to many names. Either a scientific name or a
common name. Only scientific names supported right now.}

\item{from}{(character) Data source to get data from, any combination of
gbif or bison}

\item{limit}{(numeric) Number of records to return. This is passed across
all sources. To specify different limits for each source, use the options
for each source (gbifopts, bisonopts). See Details for more.}

\item{rank}{(character) Taxonomic rank to limit search space. Used in GBIF,
but not used in BISON.}

\item{callopts}{Options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}, e.g., for
debugging curl calls, setting timeouts, etc.}

\item{gbifopts}{(list) List of named options to pass on to
\code{\link[rgbif:name_lookup]{rgbif::name_lookup()}}. See also \code{\link[=occ_names_options]{occ_names_options()}}}

\item{bisonopts}{(list) List of named options to pass on to
\code{\link[rbison:bison_tax]{rbison::bison_tax()}}. See also \code{\link[=occ_names_options]{occ_names_options()}}}
}
\description{
Search for species names across many data sources.
}
\details{
Not all 7 data sources available from the \code{\link[=occ]{occ()}} function are
available here, as not all of those sources have functionality to search
for names.

We strongly encourage you to use the \code{taxize} package if you want to
search for taxonomic or common names, convert common to scientific names,
etc. That package was built exactly for that purpose, and we only provide
a bit of name searching here in this function.
}
\examples{
\dontrun{
# Single data sources
## gbif
(res <- occ_names(query = 'Accipiter striatus', from = 'gbif'))
head(res$gbif$data[[1]])

## bison
(res <- occ_names(query = '*bear', from = 'bison'))
res$bison$data
}
}
\seealso{
Other queries: 
\code{\link{occ_names_options}()},
\code{\link{occ_options}()},
\code{\link{occ}()},
\code{\link{spocc_objects}}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixnames.r
\name{fixnames}
\alias{fixnames}
\title{Change names to be the same for each taxon.}
\usage{
fixnames(...)
}
\arguments{
\item{...}{ignored}
}
\description{
DEFUNT: This function has moved to \code{scrubr::fix_names}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{spocc_capwords}
\alias{spocc_capwords}
\title{Capitalize the first letter of a character string.}
\usage{
spocc_capwords(s, strict = FALSE, onlyfirst = FALSE)
}
\arguments{
\item{s}{A character string}

\item{strict}{Should the algorithm be strict about capitalizing.
Default: \code{FALSE}}

\item{onlyfirst}{Capitalize only first word, lowercase all others.
Useful for taxonomic names.}
}
\description{
Capitalize the first letter of a character string.
}
\examples{
\dontrun{
spocc_capwords(c('using AIC for model selection'))
spocc_capwords(c('using AIC for model selection'), strict=TRUE)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect.R
\name{inspect}
\alias{inspect}
\alias{inspect.data.frame}
\alias{inspect.occdat}
\alias{inspect.occkey}
\title{Get more data on individual occurrences}
\usage{
inspect(x, from = "gbif")

\method{inspect}{data.frame}(x, from = "gbif")

\method{inspect}{occdat}(x, from = "gbif")

\method{inspect}{occkey}(x, from = "gbif")
}
\arguments{
\item{x}{The output from \code{\link[=occ]{occ()}} call, output from call to
\code{\link[=occ2df]{occ2df()}}, or an occurrence ID as a occkey class.}

\item{from}{(character) The data provider. One of gbif, bison, inat,
or vertnet}
}
\value{
A list, with each slot named for the data source, and then
within data sources is a slot for each taxon, named by it's occurrence ID.
}
\description{
Fetches the complete record, which may or may not be the same
as requested through \code{\link[=occ]{occ()}}. Some data providers have different ways
to retrieve many occurrence records vs. single occurrence records -
and sometimes the results are more verbose when retrieving a
single occurrence record.
}
\examples{
\dontrun{
spnames <- c('Accipiter striatus', 'Spinus tristis')
out <- occ(query=spnames, from=c('gbif','bison'),
   gbifopts=list(hasCoordinate=TRUE), limit=2)
res <- occ2df(out)
inspect(res)

out <- occ(query=spnames, from='gbif', gbifopts=list(hasCoordinate=TRUE),
  limit=4)
res <- occ2df(out)
inspect(res)

# from occkeys
key <- as.gbif(res$key[1])
inspect(key)

# idigbio
spnames <- c('Accipiter striatus', 'Spinus tristis')
out <- occ(query=spnames, from='idigbio', limit=20)
inspect(out)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_options.r
\name{occ_options}
\alias{occ_options}
\title{Look up options for parameters passed to each source}
\usage{
occ_options(from = "gbif", where = "console")
}
\arguments{
\item{from}{(character) Data source to get data from, any combination of
gbif, bison, ebird, idigibio and/or vertnet. Case doesn't matter.
inat is not included here, see that package's help docs.}

\item{where}{(character) One of console (print to console) or html (opens
help page, if in non-interactive R session, prints help to console).}
}
\value{
Opens up the documentation for the function that is used internally
within the occ function for each source.
}
\description{
Look up options for parameters passed to each source
}
\details{
Any of the parameters passed to e.g. \code{\link[rgbif:occ_data]{rgbif::occ_data()}} from the
\code{rgbif} package can be passed in the associated gbifopts list
in \code{\link[=occ]{occ()}}

Note that the from parameter is lowercased within the function and is
called through match.arg first, so you can match on unique partial
strings too (e.g., 'rv' for 'rvertnet').
}
\examples{
\dontrun{
# opens up documentation for this function
occ_options()

# Open up documentation for the appropriate search function for each source
occ_options('gbif')
occ_options('ebird')
occ_options('bison')
occ_options('idigbio')
occ_options('vertnet')

# Or open in html version
occ_options('bison', 'html')
}
}
\seealso{
Other queries: 
\code{\link{occ_names_options}()},
\code{\link{occ_names}()},
\code{\link{occ}()},
\code{\link{spocc_objects}}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_names_options.R
\name{occ_names_options}
\alias{occ_names_options}
\title{Look up options for parameters passed to each source for occ_names function}
\usage{
occ_names_options(from = "gbif", where = "console")
}
\arguments{
\item{from}{(character) Data source to get data from, any combination of
gbif or bison. Case doesn't matter.}

\item{where}{(character) One of console (print to console) or html
(opens help page, if in non-interactive R session, prints help to console).}
}
\value{
Opens up the documentation for the function that is used internally
within the occ function for each source.
}
\description{
Look up options for parameters passed to each source for occ_names function
}
\details{
Any of the parameters passed to e.g. \code{\link[rgbif:name_lookup]{rgbif::name_lookup()}} from the
\code{rgbif} package can be passed in the associated gbifopts list
in \code{\link[=occ]{occ()}}.

Note that the from parameter is lowercased within the function and is
called through \code{match.arg} first, so you can match on unique partial
strings too (e.g., 'rb' for 'rbison').
}
\examples{
\dontrun{
# opens up documentation for this function
occ_names_options()

# Open up documentation for the appropriate search function for each source
occ_names_options('gbif')
occ_names_options('bison')

# Or open in html version
occ_names_options('bison', 'html')
}
}
\seealso{
Other queries: 
\code{\link{occ_names}()},
\code{\link{occ_options}()},
\code{\link{occ}()},
\code{\link{spocc_objects}}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.idigbio.R
\name{as.idigbio}
\alias{as.idigbio}
\title{Coerce occurrence keys to idigbio objects}
\usage{
as.idigbio(x, ...)
}
\arguments{
\item{x}{Various inputs, including the output from a call to \code{\link[=occ]{occ()}}
(class occdat), \code{\link[=occ2df]{occ2df()}} (class data.frame), or a list, numeric,
character, idigbiokey, or occkey.}

\item{...}{curl options; named parameters passed on to \code{httr::GET()}}
}
\value{
One or more in a list of both class idigbiokey and occkey
}
\description{
Coerce occurrence keys to idigbio objects
}
\details{
Internally, we use \code{idig_view_records}, whereas we use
\code{\link[=idig_search_records]{idig_search_records()}} in the \code{\link[=occ]{occ()}} function.
}
\examples{
\dontrun{
spnames <- c('Accipiter striatus', 'Setophaga caerulescens',
  'Spinus tristis')
out <- occ(query=spnames, from='idigbio', limit=2)
res <- occ2df(out)
(tt <- as.idigbio(out))
(uu <- as.idigbio(res))
as.idigbio(res$key[1])
as.idigbio(as.list(res$key[1:2]))
as.idigbio(tt[[1]])
as.idigbio(uu[[1]])
as.idigbio(tt[1:2])

library("dplyr")
bind_rows(lapply(tt, function(x) data.frame(unclass(x)$data)))
}
}
\seealso{
Other coercion: 
\code{\link{as.ala}()},
\code{\link{as.bison}()},
\code{\link{as.gbif}()},
\code{\link{as.inat}()},
\code{\link{as.obis}()},
\code{\link{as.vertnet}()}
}
\concept{coercion}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wkt_vis.r
\name{wkt_vis}
\alias{wkt_vis}
\title{Visualize well-known text area's on a map.}
\usage{
wkt_vis(x, zoom = 6, maptype = "terrain", browse = TRUE)
}
\arguments{
\item{x}{Input well-known text area (character)}

\item{zoom}{Zoom level, defaults to 6 (numeric)}

\item{maptype}{Map type, default is terrain (character)}

\item{browse}{Open in browser or not. If not, gives back
path to html file. Default: \code{TRUE} (logical)}
}
\description{
This can be helpful in visualizing the area in which you are searching for
occurrences with the \code{\link[=occ]{occ()}} function.
}
\details{
Uses Mapbox's map layers, openes in your default browser
}
\examples{
\dontrun{
poly <- 'POLYGON((-111.06 38.84, -110.80 39.37, -110.20 39.17, -110.20 38.90,
     -110.63 38.67, -111.06 38.84))'
wkt_vis(poly)

poly2 <- 'POLYGON((-125 38.4,-125 40.9,-121.8 40.9,-121.8 38.4,-125 38.4))'
wkt_vis(poly2)

# Multiple polygons
x <- "POLYGON((-125 38.4, -121.8 38.4, -121.8 40.9, -125 40.9, -125 38.4), 
(-115 22.4, -111.8 22.4, -111.8 30.9, -115 30.9, -115 22.4))"
wkt_vis(x)

# don't open in browser
poly2 <- 'POLYGON((-125 38.4,-125 40.9,-121.8 40.9,-121.8 38.4,-125 38.4))'
wkt_vis(poly2, browse = FALSE)
}
}
\seealso{
Other bbox: 
\code{\link{bbox2wkt}()}
}
\concept{bbox}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obis_helpers.R
\name{obis_search}
\alias{obis_search}
\title{OBIS search}
\usage{
obis_search(
  scientificName = NULL,
  size = 500,
  after = NULL,
  taxonid = NULL,
  aphiaid = NULL,
  areaid = NULL,
  datasetid = NULL,
  instituteid = NULL,
  nodeid = NULL,
  startdate = NULL,
  enddate = NULL,
  startdepth = NULL,
  enddepth = NULL,
  geometry = NULL,
  exclude = NULL,
  fields = NULL,
  ...
)
}
\arguments{
\item{size}{(integer) number of results to fetch}

\item{after}{(character) Occurrence UUID up to which to skip.}

\item{taxonid}{(character) Taxon AphiaID.}

\item{areaid}{(character) Area ID.}

\item{datasetid}{(character) Dataset UUID.}

\item{instituteid}{(character) Institute ID.}

\item{nodeid}{(character) Node UUID.}

\item{startdate}{(character) Start date formatted as YYYY-MM-DD.}

\item{enddate}{(character) End date formatted as YYYY-MM-DD.}

\item{startdepth}{(integer) Start depth, in meters.}

\item{enddepth}{(integer) End depth, in meters.}

\item{geometry}{(character) Geometry, formatted as WKT.}

\item{exclude}{(character) set of quality flags to be excluded.
one or more in a vector}

\item{fields}{(character) Field to be included in the result set.
one or more in a vector}

\item{scientificname}{(character) Scientific name. Leave empty to
include all taxa. This is what we pass your name query to}
}
\description{
OBIS search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.bison.R
\name{as.bison}
\alias{as.bison}
\title{Coerce occurrence keys to bisonkey/occkey objects}
\usage{
as.bison(x, ...)
}
\arguments{
\item{x}{Various inputs, including the output from a call to \code{\link[=occ]{occ()}}
(class occdat), \code{\link[=occ2df]{occ2df()}} (class data.frame), or a list, numeric,
character, or bisonkey, or occkey.}

\item{...}{curl options; named parameters passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
One or more in a list of both class bisonkey and occkey
}
\description{
Coerce occurrence keys to bisonkey/occkey objects
}
\details{
Internally, we use \code{\link[rbison:bison_solr]{rbison::bison_solr()}}, same function we use
internally within the \code{\link[=occ]{occ()}} function. Although, we query here with the
\code{occurrenceID} parameter to get the occurrence directly instead of
searching for it.
}
\examples{
\dontrun{
spnames <- c('Accipiter striatus', 'Setophaga caerulescens',
  'Spinus tristis')
out <- occ(query=spnames, from='bison', limit=2)
res <- occ2df(out)
(tt <- as.bison(out))
(uu <- as.bison(res))
as.bison(as.numeric(res$key[1]))
as.bison(res$key[1])
as.bison(as.list(res$key[1:2]))
as.bison(tt[[1]])
as.bison(uu[[1]])
as.bison(tt[1:2])
}
}
\seealso{
Other coercion: 
\code{\link{as.ala}()},
\code{\link{as.gbif}()},
\code{\link{as.idigbio}()},
\code{\link{as.inat}()},
\code{\link{as.obis}()},
\code{\link{as.vertnet}()}
}
\concept{coercion}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.r, R/occ_names.R
\name{spocc_objects}
\alias{spocc_objects}
\alias{print.occdat}
\alias{print.occdatind}
\alias{summary.occdat}
\alias{summary.occdatind}
\alias{print.occnames}
\title{spocc objects and their print, plot, and summary methods}
\usage{
\method{print}{occdat}(x, ...)

\method{print}{occdatind}(x, ...)

\method{summary}{occdat}(object, ...)

\method{summary}{occdatind}(object, ...)

\method{print}{occnames}(x, ...)
}
\arguments{
\item{x}{Input, of class occdatind}

\item{...}{Further args to print, plot or summary methods}

\item{object}{Input to summary methods}

\item{n}{Number of rows to show. If \code{NULL}, the default, will print
all rows if less than option \code{dplyr.print_max}. Otherwise, will
print \code{dplyr.print_min}}
}
\description{
spocc objects and their print, plot, and summary methods
}
\examples{
\dontrun{
# occdat object
res <- occ(query = 'Accipiter striatus', from = 'gbif')
res
print(res)
class(res)

# occdatind object
res$gbif
print(res$gbif)
class(res$gbif)

# print summary of occdat object
summary(res)

# print summary of occdatind object
summary(res$gbif)

# Geometry based searches print slightly differently
bounds <- c(-120, 40, -100, 45)
(res <- occ(from = "idigbio", geometry = bounds, limit = 10))
res$idigbio
## Many bounding boxes/WKT strings
bounds <- list(c(165,-53,180,-29), c(-180,-53,-175,-29))
res <- occ(from = "idigbio", geometry = bounds, limit = 10)
res$idigbio
}
}
\seealso{
Other queries: 
\code{\link{occ_names_options}()},
\code{\link{occ_names}()},
\code{\link{occ_options}()},
\code{\link{occ}()}
}
\concept{queries}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spocc_duplicates.R
\name{spocc_duplicates}
\alias{spocc_duplicates}
\title{A note about duplicate occurrence records}
\description{
BEWARE: spocc provides you a nice interface to many data providers for
species occurrence data. However, in cases where you request data from
GBIF \emph{in addition} to other data sources, there could be duplicate records.
This is because GBIF is, to use an ecology  analogy, a top predator, and
pulls in data from lower nodes in the food chain. For example, iNaturalist
provides data to GBIF, so if you search for occurrence records for
\emph{Pinus contorta} from iNaturalist and GBIF, you could get, for example,
20 of the same records.

We think a single R interface to many occurrence record providers
will provide a consistent way to work with occurrence data, making
analyses and vizualizations more repeatable across providers.

For cleaning data, see packages \code{scrubr}
(\url{https://cran.r-project.org/package=scrubr}) and \code{CoordinateCleaner}
(\url{https://cran.r-project.org/package=CoordinateCleaner})

Do get in touch with us if you have concerns, have ideas for eliminating
duplicates
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wkt_bbox.R
\name{bbox2wkt}
\alias{bbox2wkt}
\alias{wkt2bbox}
\title{Convert a bounding box to a Well Known Text polygon, and a WKT to a
bounding box}
\usage{
bbox2wkt(minx = NA, miny = NA, maxx = NA, maxy = NA, bbox = NULL)

wkt2bbox(wkt)
}
\arguments{
\item{minx}{Minimum x value, or the most western longitude}

\item{miny}{Minimum y value, or the most southern latitude}

\item{maxx}{Maximum x value, or the most eastern longitude}

\item{maxy}{Maximum y value, or the most northern latitude}

\item{bbox}{A vector of length 4, with the elements: minx, miny, maxx, maxy}

\item{wkt}{A Well Known Text string}
}
\value{
bbox2wkt returns an object of class charactere, a Well Known Text
string of the form
'POLYGON((minx miny, maxx miny, maxx maxy, minx maxy, minx miny))'

wkt2bbox returns a numeric vector of length 4, like c(minx, miny,
maxx, maxy).
}
\description{
Convert a bounding box to a Well Known Text polygon, and a WKT to a
bounding box
}
\examples{
# Convert a bounding box to a WKT

## Pass in a vector of length 4 with all values
bbox2wkt(bbox = c(-125.0,38.4,-121.8,40.9))

## Or pass in each value separately
bbox2wkt(-125.0, 38.4, -121.8, 40.9)

# Convert a WKT object to a bounding box
wkt <- "POLYGON((-125 38.4,-125 40.9,-121.8 40.9,-121.8 38.4,-125 38.4))"
wkt2bbox(wkt)

identical(
 bbox2wkt(-125.0, 38.4, -121.8, 40.9),
 "POLYGON((-125 38.4,-121.8 38.4,-121.8 40.9,-125 40.9,-125 38.4))"
)

identical(
 c(-125.0, 38.4, -121.8, 40.9),
 as.numeric(
   wkt2bbox(
     "POLYGON((-125 38.4,-125 40.9,-121.8 40.9,-121.8 38.4,-125 38.4))"
   )
 )
)
}
\seealso{
Other bbox: 
\code{\link{wkt_vis}()}

Other bbox: 
\code{\link{wkt_vis}()}
}
\concept{bbox}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spocc-package.R
\docType{package}
\name{spocc-package}
\alias{spocc-package}
\alias{spocc}
\title{Interface to many species occurrence data sources}
\description{
A programmatic interface to many species occurrence data
sources, including GBIF, USGS's BISON, iNaturalist, Berkeley Ecoinformatics
Engine, eBird, iDigBio, VertNet, OBIS, and ALA. Includes
functionality for retrieving species occurrence data, and
combining that data.
}
\section{Package API}{


The main function to use is \code{\link[=occ]{occ()}} - a single interface to
many species occurrence databases (see below for a list).

Other functions include:
\itemize{
\item \code{\link[=occ2df]{occ2df()}} - Combine results from \code{occ} into a
data.frame
\item \code{\link[=wkt_vis]{wkt_vis()}} - Visualize WKT strings (used to define
geometry based searches for some data sources) in an interactive map
}
}

\section{Currently supported species occurrence data sources}{


\tabular{ll}{
Provider \tab Web \cr
GBIF \tab \url{http://www.gbif.org/} \cr
BISON \tab https://bison.usgs.gov/ \cr
eBird \tab \url{http://ebird.org/content/ebird/} \cr
iNaturalist \tab \url{http://www.inaturalist.org/} \cr
VertNet \tab \url{http://vertnet.org/} \cr
iDigBio \tab \url{https://www.idigbio.org/} \cr
OBIS \tab \url{http://www.iobis.org/} \cr
ALA \tab \url{http://www.ala.org.au/}
}
}

\section{Duplicates}{


See \code{\link[=spocc_duplicates]{spocc_duplicates()}} for more.
}

\section{Clean data}{


All data cleaning functionality is in a new package: \code{scrubr}
(\url{https://github.com/ropensci/scrubr}).
On CRAN: \url{https://cran.r-project.org/package=scrubr}.
See also package
\url{https://cran.r-project.org/package=CoordinateCleaner}
}

\section{Make maps}{


All mapping functionality is now in a separate package: \verb{mapr`` (<https://github.com/ropensci/mapr>) (formerly known as }spoccutils`).
On CRAN: \url{https://cran.r-project.org/package=mapr}
}

\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.ala.R
\name{as.ala}
\alias{as.ala}
\title{Coerce occurrence keys to ALA id objects}
\usage{
as.ala(x, ...)
}
\arguments{
\item{x}{Various inputs, including the output from a call to
\code{\link[=occ]{occ()}} (class occdat), \code{\link[=occ2df]{occ2df()}} (class data.frame),
or a list, numeric, alakey, or occkey.}

\item{...}{curl options; named parameters passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
One or more in a list of both class alakey and occkey
}
\description{
Coerce occurrence keys to ALA id objects
}
\examples{
\dontrun{
spnames <- c('Barnardius zonarius', 'Grus rubicunda', 'Cracticus tibicen')
out <- occ(query=spnames, from='ala', limit=2)
(res <- occ2df(out))
(tt <- as.ala(out))
as.ala(x = res$key[1])
}
}
\seealso{
Other coercion: 
\code{\link{as.bison}()},
\code{\link{as.gbif}()},
\code{\link{as.idigbio}()},
\code{\link{as.inat}()},
\code{\link{as.obis}()},
\code{\link{as.vertnet}()}
}
\concept{coercion}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ.r
\name{occ}
\alias{occ}
\title{Search for species occurrence data across many data sources.}
\usage{
occ(
  query = NULL,
  from = "gbif",
  limit = 500,
  start = NULL,
  page = NULL,
  geometry = NULL,
  has_coords = NULL,
  ids = NULL,
  date = NULL,
  callopts = list(),
  gbifopts = list(),
  bisonopts = list(),
  inatopts = list(),
  ebirdopts = list(),
  vertnetopts = list(),
  idigbioopts = list(),
  obisopts = list(),
  alaopts = list(),
  throw_warnings = TRUE
)
}
\arguments{
\item{query}{(character) One to many scientific names. See Details for what parameter
in each data source we query. Note: ebird now expects species codes instead of
scientific names - we pass you name through \code{\link[rebird:species_code]{rebird::species_code()}} internally}

\item{from}{(character) Data source to get data from, any combination of gbif, bison,
inat, ebird, and/or vertnet}

\item{limit}{(numeric) Number of records to return. This is passed across all sources.
To specify different limits for each source, use the options for each source (gbifopts,
bisonopts, inatopts, and ebirdopts). See Details for more.
Default: 500 for each source. BEWARE: if you have a lot of species to query for (e.g.,
n = 10), that's 10 * 500 = 5000, which can take a while to collect. So, when you first query,
set the limit to something smallish so that you can get a result quickly, then do more as
needed.}

\item{start, page}{(integer) Record to start at or page to start at. See \code{Paging} in
Details for how these parameters are used internally. Optional}

\item{geometry}{(character or nmeric) One of a Well Known Text (WKT) object, a vector of
length 4 specifying a bounding box, or an sf object (sfg, sfc, or sf). This parameter
searches for occurrences inside a
polygon - converted to a polygon from whatever user input is given. A WKT shape written as
\verb{POLYGON((30.1 10.1, 20 40, 40 40, 30.1 10.1))} would be queried as is,
i.e. http://bit.ly/HwUSif. See Details for more examples of WKT objects. The format of a
bounding box is \verb{min-longitude, min-latitude, max-longitude, max-latitude}. Geometry
is not possible with vertnet right now, but should be soon. See Details for more info
on geometry inputs.}

\item{has_coords}{(logical) Only return occurrences that have lat/long data. This works
for gbif, rinat, idigbio, and vertnet, but is ignored for ebird and
bison data sources. You can easily though remove records without lat/long data.}

\item{ids}{Taxonomic identifiers. This can be a list of length 1 to many. See examples for
usage. Currently, identifiers for only 'gbif' and 'bison' for parameter 'from' supported. If
this parameter is used, query parameter can not be used - if it is, a warning is thrown.}

\item{date}{(character/Date) A length 2 vector containing two dates of the form
YYY-MM-DD. These can be character of Date class. These are used to do a date range search.
Of course there are other types of date searches one may want to do but date range
seems like the most common date search use case.}

\item{callopts}{Options passed on to \link[crul:HttpClient]{crul::HttpClient}, e.g.,
for debugging curl calls, setting timeouts, etc.}

\item{gbifopts}{(list) List of named options to pass on to
\code{\link[rgbif:occ_search]{rgbif::occ_search()}}. See also \code{\link[=occ_options]{occ_options()}}}

\item{bisonopts}{(list) List of named options to pass on to \code{\link[rbison:bison]{rbison::bison()}}.
See also \code{\link[=occ_options]{occ_options()}}}

\item{inatopts}{(list) List of named options to pass on to internal function
\code{get_inat_obs}}

\item{ebirdopts}{(list) List of named options to pass on to
\code{\link[rebird:ebirdregion]{rebird::ebirdregion()}} or \code{\link[rebird:ebirdgeo]{rebird::ebirdgeo()}}. See also \code{\link[=occ_options]{occ_options()}}}

\item{vertnetopts}{(list) List of named options to pass on to
\code{\link[rvertnet:searchbyterm]{rvertnet::searchbyterm()}}. See also \code{\link[=occ_options]{occ_options()}}.}

\item{idigbioopts}{(list) List of named options to pass on to
\code{\link[ridigbio:idig_search_records]{ridigbio::idig_search_records()}}. See also \code{\link[=occ_options]{occ_options()}}.}

\item{obisopts}{(list) List of named options to pass on to internal function.
See  https://api.obis.org/#/Occurrence/get_occurrence and \link{obis_search} for
what parameters can be used.}

\item{alaopts}{(list) List of named options to pass on to internal function.
See \verb{Occurrence search} part of the API docs at
\url{http://api.ala.org.au/#ws3} for possible parameters.}

\item{throw_warnings}{(logical) \code{occ()} collects errors returned from each
data provider when they occur, and are accessible in the \verb{$meta$errors} slot
for each data provider. If you set \code{throw_warnings=TRUE}, we give these
request errors as warnings with \code{\link[=warning]{warning()}}. if \code{FALSE}, we don't give warnings,
but you can still access them in the output.}
}
\value{
an object of class \code{occdat}, with a print method to give a brief
summary. The print method only shows results for those that have some
results (those with no results are not shown). The \code{occdat} class is just
a thin wrapper around a named list, wher the top level names are the
data sources:
\itemize{
\item gbif
\item bison
\item inat
\item ebird
\item vertnet
\item idigbio
\item obis
\item ala
}

Note that you only get data back for sources that were specified in the \code{from}
parameter. All others are present, but empty.

Then within each data source is an object of class \code{occdatind} holding another
named list that contains:
\itemize{
\item meta: metadata
\itemize{
\item source: the data source name (e.g., "gbif")
\item time: time the request was sent
\item found: number of records found (number found across all queries)
\item returned: number of records returned (number of rows in all data.frame's
in the \code{data} slot)
\item type: query type, only "sci" for scientific
\item opts: a named list with the options you sent to the data source
\item errors: a character vector of errors returned, if any occurred
}
\item data: named list of data.frame's, named by the queries sent
}
}
\description{
Search on a single species name, or many. And search across a single
or many data sources.
}
\details{
The \code{occ} function is an opinionated wrapper
around the rgbif, rbison, rinat, rebird, rvertnet and
ridigbio packages (as well as internal custom wrappers around some data
sources) to allow data access from a single access point. We take
care of making sure you get useful objects out at the cost of
flexibility/options - although you can still set options for each of the
packages via the gbifopts, bisonopts, inatopts, etc. parameters.
}
\section{Inputs}{

All inputs to \code{occ} are one of:
\itemize{
\item scientific name
\item taxonomic id
\item geometry as bounds, WKT, os Spatial classes
}

To search by common name, first use \code{\link[=occ_names]{occ_names()}} to find scientic names or
taxonomic IDs, then feed those to this function. Or use the \code{taxize} package
to get names and/or IDs to use here.
}

\section{Using the query parameter}{

When you use the \code{query} parameter, we pass your search terms on to parameters
within functions that query data sources you specify. Those parameters are:
\itemize{
\item rgbif - \code{scientificName} in the \code{\link[rgbif:occ_search]{rgbif::occ_search()}} function - API
parameter: same as the \code{occ} parameter
\item rebird - \code{species} in the \code{\link[rebird:ebirdregion]{rebird::ebirdregion()}} or
\code{\link[rebird:ebirdgeo]{rebird::ebirdgeo()}} functions, depending on whether you set
\code{method="ebirdregion"} or \code{method="ebirdgeo"} - API parameters: \code{sci} for both
\code{\link[rebird:ebirdregion]{rebird::ebirdregion()}} and \code{\link[rebird:ebirdgeo]{rebird::ebirdgeo()}}
\item rbison - \code{species} or \code{scientificName} in the \code{\link[rbison:bison]{rbison::bison()}} or
\code{\link[rbison:bison_solr]{rbison::bison_solr()}} functions, respectively. If you don't pass anything to
\code{geometry} parameter we use \code{bison_solr}, and if you do we use \code{bison} - API
parameters: same as \code{occ} parameters
\item rvertnet - \code{taxon} in the \code{\link[rvertnet:vertsearch]{rvertnet::vertsearch()}} function - API
parameter: \code{q}
\item ridigbio - \code{scientificname} in the \code{\link[ridigbio:idig_search_records]{ridigbio::idig_search_records()}}
function - API parameter: \code{scientificname}
\item inat - internal function - API parameter: \code{q}
\item obis - internal function - API parameter: \code{scientificName}
\item ala - internal function - API parameter: \code{q}
}

If you have questions about how each of those parameters behaves with respect to
the terms you pass to it, lookup documentation for those functions, or get in touch
at the development repository \url{https://github.com/ropensci/spocc/issues}
}

\section{iDigBio notes}{

When searching iDigBio note that by deafult we set \code{fields = "all"}, so that we return
a richer suite of fields than the \code{ridigbio} R client gives by default. But you can
changes this by passing in a \code{fields} parameter to \code{idigbioopts} parameter with
the specific fields you want.

Maximum of 100,000 results are allowed to be returned. See
\url{https://github.com/iDigBio/ridigbio/issues/33}
}

\section{BISON notes}{

We use two different functions when you request data from \code{bison}. We use
\code{\link[rbison:bison_solr]{rbison::bison_solr()}} by default as it's more flexible. If you pass a value to the
\code{geometry} parameter we use \code{\link[rbison:bison]{rbison::bison()}}. We'd prefer to just use one function
to simplify things, but \code{\link[rbison:bison_solr]{rbison::bison_solr()}} doesn't support geometry queries.
}

\section{iNaturalist notes}{

We're using the iNaturalist API, docs at
https://api.inaturalist.org/v1/docs/#!/Observations/get_observations

API rate limits: max of 100 requests per minute, though they ask that you try to keep it
to 60 requests per minute or lower. If they notice usage that has serious impact on their
performance they may institute blocks without notification.

There is a hard limit 0f 10,000 observations with the iNaturalist API. We do paging
internally so you may not see this aspect, but for example, if you request 12,000
records, you won't be able to get that many. The API will error at anything more than
10,000. We now error if you request more than 10,000 from iNaturalist. There are
some alternatives:
\itemize{
\item Consider exporting data while logged in
to your iNaturalist account, or the iNaturalist research grade observations within
GBIF - see https://www.gbif.org/dataset/50c9509d-22c7-4a22-a47d-8c48425ef4a7 - at
time of this writing it has 8.5 million observations.
\item Search for iNaturalist data within GBIF. e.g., the following searches for iNaturalist
data within GBIF and allows more than 10,000 records:
``
}
}

\section{limit parameter}{

The \code{limit} parameter is set to a default of 25. This means that you will get \strong{up to}
25 results back for each data source you ask for data from. If there are no results for a
particular source, you'll get zero back; if there are 8 results for a particular source, you'll
get 8 back. If there are 26 results for a particular source, you'll get 25 back. You can always
ask for more or less back by setting the limit parameter to any number. If you want to request
a different number for each source, pass the appropriate parameter to each data source via the
respective options parameter for each data source.
}

\section{WKT}{

WKT objects are strings of pairs of lat/long coordinates that define a shape. Many classes
of shapes are supported, including POLYGON, POINT, and MULTIPOLYGON. Within each defined shape
define all vertices of the shape with a coordinate like 30.1 10.1, the first of which is the
latitude, the second the longitude.

Examples of valid WKT objects:
\itemize{
\item 'POLYGON((30.1 10.1, 10 20, 20 60, 60 60, 30.1 10.1))'
\item 'POINT((30.1 10.1))'
\item 'LINESTRING(3 4,10 50,20 25)'
\item 'MULTIPOINT((3.5 5.6),(4.8 10.5))")'
\item 'MULTILINESTRING((3 4,10 50,20 25),(-5 -8,-10 -8,-15 -4))'
\item 'MULTIPOLYGON(((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2)),((6 3,9 2,9 4,6 3)))'
\item 'GEOMETRYCOLLECTION(POINT(4 6),LINESTRING(4 6,7 10))'
}

Only POLYGON objects are currently supported.

Getting WKT polygons or bounding boxes. We will soon introduce a function to help you select
a bounding box but for now, you can use a few sites on the web.
\itemize{
\item Bounding box - \url{http://boundingbox.klokantech.com/}
\item Well known text - \url{http://arthur-e.github.io/Wicket/sandbox-gmaps3.html}
}
}

\section{geometry parameter}{

The behavior of the \code{occ} function with respect to the \code{geometry} parameter
varies depending on the inputs to the \code{query} parameter. Here are the options:
\itemize{
\item geometry (single), no query - If a single bounding box/WKT string passed in,
and no query, a single query is made against each data source.
\item geometry (many), no query - If many bounding boxes/WKT strings are passed in,
we do a separate query for each bounding box/WKT string against each data source.
\item geometry (single), query - If a single bounding box/WKT string passed in,
and a single query, we do a single query against each data source.
\item geometry (many), query - If many bounding boxes/WKT strings are passed in,
and a single query, we do a separate query for each bounding box/WKT string with the
same queried name against each data source.
\item geometry (single), many query - If a single bounding box/WKT string passed in,
and many names to query, we do a separate query for each name, using the same geometry,
for each data source.
\item geometry (many), many query - If many bounding boxes/WKT strings are passed in,
and many names to query, this poses a problem for all data sources, none of which
accept many bounding boxes of WKT strings. So, in this scenario, we loop over each
name and each geometry query, and then re-combine by queried name, so that you get
back a single group of data for each name.
}
}

\section{Geometry options by data provider}{

\strong{wkt & bbox allowed, see WKT section above}
\itemize{
\item gbif
\item bison
\item obis
\item ala
}

\strong{bbox only}
\itemize{
\item inat
\item idigbio
}

\strong{No spatial search allowed}
\itemize{
\item ebird
\item vertnet
}
}

\section{Notes on the date parameter}{

Date searches with the \code{date} parameter are allowed for all sources
except ebird.

Notes on some special cases
\itemize{
\item idigbio: We search on the \code{datecollected} field. Other date fields can be
searched on, but we chose \code{datecollected} as it seemed most appropriate.
\item vertnet: If you want more flexible date searches, you can pass various
types of date searches to \code{vertnetopts}. See \code{\link[rvertnet:searchbyterm]{rvertnet::searchbyterm()}}
for more information
\item ala: There's some issues with the dates returned from ALA. They are
returned as time stamps, and some seem to be malformed. So do beware
of using ALA dates for important things.
}

Get in touch if you have other date search use cases you think
are widely useful
}

\section{Paging}{

All data sources respond to the \code{limit} parameter passed to \code{occ}.

Data sources, however, vary as to whether they respond to an offset. Here's
the details on which data sources will respond to \code{start} and which
to the \code{page} parameter:
\itemize{
\item gbif - Responds to \code{start}. Default: 0
\item bison - Responds to \code{start}. Default: 0
\item inat - Responds to \code{page}. Default: 1
\item ebird - No paging, both \code{start} and \code{page} ignored.
\item vertnet - No paging implemented here, both \code{start} and \code{page}
ignored. VertNet does have a form of paging, but it uses a cursor, and can't
easily be included  here via parameters. However, \code{rvertnet} does paging
internally for you.  For example, the max records per request for VertNet is
1000; if you request 2000 records, we'll do the first request, and do the
second request for you automatically.
\item idigbio - Responds to \code{start}. Default: 0
\item obis - Does not respond to \code{start}. They only allow a starting occurrence
UUID up to which to skip. So order of results matters a great deal of course.
To paginate with OBIS, do e.g.
\code{obisopts = list(after = "017b7818-5b2c-4c88-9d76-f4471afe5584")}; \code{after} can
be combined with the \code{limit} value you pass in to the main \code{occ()} function
call. See \link{obis_search} for what parameters can be used.
\item ala - Responds to \code{start}. Default: 0
}
}

\section{Photographs}{

The iNaturalist data source provides photographs of the records returned,
if available. For example, the following will give photos from inat:
\code{occ(query = 'Danaus plexippus', from = 'inat')$inat$data$Danaus_plexippus$photos}
}

\section{BEWARE}{

In cases where you request data from multiple providers, especially when
including GBIF, there could be duplicate records since many providers' data eventually
ends up with GBIF. See \code{\link[=spocc_duplicates]{spocc_duplicates()}} for more.
}

\examples{
\dontrun{
# Single data sources
(res <- occ(query = 'Accipiter striatus', from = 'gbif', limit = 5))
res$gbif
(res <- occ(query = 'Accipiter striatus', from = 'ebird', limit = 50))
res$ebird
(res <- occ(query = 'Danaus plexippus', from = 'inat', limit = 50,
  has_coords = TRUE))
res$inat
res$inat$data
data.table::rbindlist(res$inat$data$Danaus_plexippus$photos)
(res <- occ(query = 'Bison bison', from = 'bison', limit = 50))
res$bison
(res <- occ(query = 'Bison bison', from = 'vertnet', limit = 5))
res$vertnet
res$vertnet$data$Bison_bison
occ2df(res)

# Paging
one <- occ(query = 'Accipiter striatus', from = 'gbif', limit = 5)
two <- occ(query = 'Accipiter striatus', from = 'gbif', limit = 5, start = 5)
one$gbif
two$gbif

# iNaturalist limits: they allow at most 10,000; query through GBIF to get
# more than 10,000
# See https://www.gbif.org/dataset/50c9509d-22c7-4a22-a47d-8c48425ef4a7
# x <- occ(query = 'Danaus plexippus', from = 'gbif', limit = 10100, 
#   gbifopts = list(datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7"))
# x$gbif

# Date range searches across data sources
## Not possible for ebird
## bison
occ(query = 'Acer', date = c('2010-08-08', '2010-08-21'), from = 'bison', limit=5)
## ala
occ(date = c('2018-01-01T00:00:00Z', '2018-03-28T00:00:00Z'), from = 'ala', limit = 5)
## gbif
occ(query = 'Accipiter striatus', date = c('2010-08-01', '2010-08-31'), from = 'gbif', limit=5)
## vertnet
occ(query = 'Mustela nigripes', date = c('1990-01-01', '2015-12-31'), from = 'vertnet', limit=5)
## idigbio
occ(query = 'Acer', date = c('2010-01-01', '2015-12-31'), from = 'idigbio', limit=5)
## obis
occ(query = 'Mola mola', date = c('2015-01-01', '2015-12-31'), from = 'obis', limit=5)
## inat
occ(query = 'Danaus plexippus', date = c('2015-01-01', '2015-12-31'), from = 'inat', limit=5)


# Restrict to records with coordinates
occ(query = "Acer", from = "idigbio", limit = 5, has_coords = TRUE)

occ(query = 'Setophaga caerulescens', from = 'ebird', ebirdopts = list(loc='US'))
occ(query = 'Spinus tristis', from = 'ebird', ebirdopts =
   list(method = 'ebirdgeo', lat = 42, lng = -76, dist = 50))

# idigbio data
## scientific name search
occ(query = "Acer", from = "idigbio", limit = 5)
occ(query = "Acer", from = "idigbio", idigbioopts = list(offset = 5, limit  = 3))
## geo search
bounds <- c(-120, 40, -100, 45)
occ(from = "idigbio", geometry = bounds, limit = 10)
## just class arachnida, spiders
occ(idigbioopts = list(rq = list(class = 'arachnida')), from = "idigbio", limit = 10)
## search certain recordsets
sets <- c("1ffce054-8e3e-4209-9ff4-c26fa6c24c2f",
    "8dc14464-57b3-423e-8cb0-950ab8f36b6f", 
    "26f7cbde-fbcb-4500-80a9-a99daa0ead9d")
occ(idigbioopts = list(rq = list(recordset = sets)), from = "idigbio", limit = 10)

# Many data sources
(out <- occ(query = 'Pinus contorta', from=c('gbif','bison','vertnet'), limit=10))

## Select individual elements
out$gbif
out$gbif$data
out$vertnet

## Coerce to combined data.frame, selects minimal set of
## columns (name, lat, long, provider, date, occurrence key)
occ2df(out)

# Pass in limit parameter to all sources. This limits the number of occurrences
# returned to 10, in this example, for all sources, in this case gbif and inat.
occ(query='Pinus contorta', from=c('gbif','inat'), limit=10)

# Geometry
## Pass in geometry parameter to all sources. This constraints the search to the
## specified polygon for all sources, gbif and bison in this example.
## Check out http://arthur-e.github.io/Wicket/sandbox-gmaps3.html to get a WKT string
occ(query='Accipiter', from='gbif',
   geometry='POLYGON((30.1 10.1, 10 20, 20 60, 60 60, 30.1 10.1))')
occ(query='Helianthus annuus', from='bison', limit=50,
   geometry='POLYGON((-111.06 38.84, -110.80 39.37, -110.20 39.17, -110.20 38.90,
                      -110.63 38.67, -111.06 38.84))')

## Or pass in a bounding box, which is automatically converted to WKT (required by GBIF)
## via the bbox2wkt function. The format of a bounding box is
## [min-longitude, min-latitude, max-longitude, max-latitude].
occ(query='Accipiter striatus', from='gbif', geometry=c(-125.0,38.4,-121.8,40.9))

## lots of results, can see how many by indexing to meta
res <- occ(query='Accipiter striatus', from='gbif',
   geometry='POLYGON((-69.9 49.2,-69.9 29.0,-123.3 29.0,-123.3 49.2,-69.9 49.2))')
res$gbif

## You can pass in geometry to each source separately via their opts parameter, at
## least those that support it. Note that if you use rinat, you reverse the order, with
## latitude first, and longitude second, but here it's the reverse for consistency across
## the spocc package
bounds <- c(-125.0,38.4,-121.8,40.9)
occ(query = 'Danaus plexippus', from="inat", geometry=bounds)

## Passing geometry with multiple sources
occ(query = 'Danaus plexippus', from=c("inat","gbif"), geometry=bounds)

## Using geometry only for the query
### A single bounding box
occ(geometry = bounds, from = "gbif", limit=50)
### Many bounding boxes
occ(geometry = list(c(-125.0,38.4,-121.8,40.9), c(-115.0,22.4,-111.8,30.9)), from = "gbif")

## Many geometry and many names
res <- occ(query = c('Danaus plexippus', 'Accipiter striatus'),
   geometry = list(c(-125.0,38.4,-121.8,40.9), c(-115.0,22.4,-111.8,30.9)), from = "bison")
res

## Geometry only with WKT
wkt <- 'POLYGON((-98.9 44.2,-89.1 36.6,-116.7 37.5,-102.5 39.6,-98.9 44.2))'
occ(from = "gbif", geometry = wkt, limit = 10)

# Specify many data sources, another example
ebirdopts = list(loc = 'US'); gbifopts  =  list(country = 'US')
out <- occ(query = 'Setophaga caerulescens', from = c('gbif','inat','bison','ebird'),
    gbifopts = gbifopts, ebirdopts = ebirdopts, limit=20)
occ2df(out)

# Pass in many species names, combine just data to a single data.frame, and
# first six rows
spnames <- c('Accipiter striatus', 'Setophaga caerulescens', 'Spinus tristis')
(out <- occ(query = spnames, from = 'gbif', gbifopts = list(hasCoordinate = TRUE), limit=25))
df <- occ2df(out)
head(df)

# no query, geometry, or ids passed
## many dataset keys to gbif
dsets <- c("14f3151a-e95d-493c-a40d-d9938ef62954", "f934f8e2-32ca-46a7-b2f8-b032a4740454")
occ(limit = 20, from = "gbif", gbifopts = list(datasetKey = dsets))
## class name to idigbio
occ(limit = 20, from = "idigbio", idigbioopts = list(rq = list(class = 'arachnida')))

# taxize integration
## You can pass in taxonomic identifiers
library("taxize")
(ids <- get_ids(c("Chironomus riparius","Pinus contorta"), db = c('itis','gbif')))
occ(ids = ids[[1]], from='bison', limit=20)
occ(ids = ids, from=c('bison','gbif'), limit=20)

(ids <- get_ids("Chironomus riparius", db = 'gbif'))
occ(ids = ids, from='gbif', limit=20)

(ids <- get_gbifid("Chironomus riparius"))
occ(ids = ids, from='gbif', limit=20)

(ids <- get_tsn('Accipiter striatus'))
occ(ids = ids, from='bison', limit=20)

## sf classes
library("sp")
library("sf")
one <- Polygon(cbind(c(91,90,90,91), c(30,30,32,30)))
spone = Polygons(list(one), "s1")
sppoly = SpatialPolygons(list(spone), as.integer(1))

## single polygon in a sf class
x <- st_as_sf(sppoly)
out <- occ(geometry = x, limit=50)
out$gbif$data
mapr::map_leaflet(out)

## single polygon in a sfc class
x <- st_as_sf(sppoly)
out <- occ(geometry = x[[1]], limit=50)
out$gbif$data

## single polygon in a sf POLYGON class
x <- st_as_sf(sppoly)
x <- unclass(x[[1]])[[1]]
class(x)
out <- occ(geometry = x, limit=50)
out$gbif$data

## two polygons in an sf class
one <- Polygon(cbind(c(-121.0,-117.9,-121.0,-121.0), c(39.4, 37.1, 35.1, 39.4)))
two <- Polygon(cbind(c(-123.0,-121.2,-122.3,-124.5,-123.5,-124.1,-123.0),
                     c(44.8,42.9,41.9,42.6,43.3,44.3,44.8)))
spone = Polygons(list(one), "s1")
sptwo = Polygons(list(two), "s2")
sppoly = SpatialPolygons(list(spone, sptwo), 1:2)
sppoly_df <- SpatialPolygonsDataFrame(sppoly, 
   data.frame(a=c(1,2), b=c("a","b"), c=c(TRUE,FALSE),
   row.names=row.names(sppoly)))
x <- st_as_sf(sppoly_df)
out <- occ(geometry = x, limit=50)
out$gbif$data


# curl debugging
occ(query = 'Accipiter striatus', from = 'gbif', limit=10, 
 callopts=list(verbose = TRUE))
occ(query = 'Accipiter striatus', from = 'bison', limit=10, 
 callopts=list(verbose = TRUE))
occ(query = 'Accipiter striatus', from = 'inat', 
 callopts=list(verbose = TRUE))
occ(query = 'Mola mola', from = 'obis', limit = 200, 
 callopts = list(verbose = TRUE))

########## More thorough data source specific examples
# idigbio
## scientific name search
res <- occ(query = "Acer", from = "idigbio", limit = 5)
res$idigbio

## geo search
### bounding box
bounds <- c(-120, 40, -100, 45)
occ(from = "idigbio", geometry = bounds, limit = 10)
### wkt
# wkt <- 'POLYGON((-69.9 49.2,-69.9 29.0,-123.3 29.0,-123.3 49.2,-69.9 49.2))'
wkt <- 'POLYGON((-98.9 44.2,-89.1 36.6,-116.7 37.5,-102.5 39.6,-98.9 44.2))'
occ(from = "idigbio", geometry = wkt, limit = 10)

## limit fields returned
occ(query = "Acer", from = "idigbio", limit = 5,
   idigbioopts = list(fields = "scientificname"))

## offset and max_items
occ(query = "Acer", from = "idigbio", limit = 5,
   idigbioopts = list(offset = 10))

## sort
occ(query = "Acer", from = "idigbio", limit = 5,
   idigbioopts = list(sort = TRUE))$idigbio
occ(query = "Acer", from = "idigbio", limit = 5,
   idigbioopts = list(sort = FALSE))$idigbio

## more complex queries
### parameters passed to "rq", get combined with the name queried
occ(query = "Acer", from = "idigbio", limit = 5,
   idigbioopts = list(rq = list(basisofrecord="fossilspecimen")))$idigbio

#### NOTE: no support for multipolygons yet
## WKT's are more flexible than bounding box's. You can pass in a WKT with multiple
## polygons like so (you can use POLYGON or MULTIPOLYGON) when specifying more than one
## polygon. Note how each polygon is in it's own set of parentheses.
# occ(query='Accipiter striatus', from='gbif',
#    geometry='MULTIPOLYGON((30 10, 10 20, 20 60, 60 60, 30 10),
#                           (30 10, 10 20, 20 60, 60 60, 30 10))')

# OBIS examples
## basic query
(res <- occ(query = 'Mola mola', from = 'obis', limit = 200))
## get to obis data
res$obis
## get obis + gbif data
(res <- occ(query = 'Mola mola', from = c('obis', 'gbif'), limit = 200))
res$gbif
res$obis
## no match found
(res <- occ(query = 'Linguimaera thomsonia', from = 'obis'))
## geometry query
geometry <- "POLYGON((8.98 48.05,15.66 48.05,15.66 45.40,8.98 45.40,8.98 48.05))"
(res <- occ(from = 'obis', geometry = geometry, limit = 50))
res$obis

## Pass in spatial classes
## sp classes no longer supported

## Paging
(res1 <- occ(query = 'Mola mola', from = 'obis', limit = 10))
occ_ids <- res1$obis$data$Mola_mola$id
(res2 <- occ(query = 'Mola mola', from = 'obis',
  limit = 10, obisopts = list(after = occ_ids[length(occ_ids)])))
res1$obis
res2$obis
## Pass in any parameters to obisopts as a list
(res <- occ(query = 'Mola mola', from = 'obis', 
   obisopts = list(startdepth = 40, enddepth = 50)))
min(res$obis$data$Mola_mola$minimumDepthInMeters, na.rm=TRUE)
max(res$obis$data$Mola_mola$maximumDepthInMeters, na.rm=TRUE)


# ALA examples
## basic query
(res <- occ(query = 'Alaba vibex', from = 'ala', limit = 200))
## get to ala data
res$ala
occ2df(res)

# geometry search
(x <- occ(query = "Macropus", from = 'ala',
  geometry = "POLYGON((145 -37,150 -37,150 -30,145 -30,145 -37))"))
x$ala
occ2df(x)
}
}
\seealso{
Other queries: 
\code{\link{occ_names_options}()},
\code{\link{occ_names}()},
\code{\link{occ_options}()},
\code{\link{spocc_objects}}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.obis.R
\name{as.obis}
\alias{as.obis}
\title{Coerce occurrence keys to obis id objects}
\usage{
as.obis(x, ...)
}
\arguments{
\item{x}{Various inputs, including the output from a call to
\code{\link[=occ]{occ()}} (class occdat), \code{\link[=occ2df]{occ2df()}} (class data.frame),
or a list, numeric, obiskey, or occkey.}

\item{...}{curl options; named parameters passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
One or more in a list of both class obiskey and occkey
}
\description{
Coerce occurrence keys to obis id objects
}
\examples{
\dontrun{
spnames <- c('Mola mola', 'Loligo vulgaris', 'Stomias boa')
out <- occ(query=spnames, from='obis', limit=2)
(res <- occ2df(out))
(tt <- as.obis(out))
(uu <- as.obis(res))
as.obis(x = res$key[1])
as.obis(as.list(res$key[1:2]))
as.obis(tt[[1]])
as.obis(uu[[1]])
as.obis(tt[1:2])

library("data.table")
rbindlist(lapply(tt, "[[", "results"),
  use.names = TRUE, fill = TRUE)
}
}
\seealso{
Other coercion: 
\code{\link{as.ala}()},
\code{\link{as.bison}()},
\code{\link{as.gbif}()},
\code{\link{as.idigbio}()},
\code{\link{as.inat}()},
\code{\link{as.vertnet}()}
}
\concept{coercion}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.inat.R
\name{as.inat}
\alias{as.inat}
\title{Coerce occurrence keys to iNaturalist id objects}
\usage{
as.inat(x, ...)
}
\arguments{
\item{x}{Various inputs, including the output from a call to
\code{\link[=occ]{occ()}} (class occdat), \code{\link[=occ2df]{occ2df()}} (class data.frame), or a list, numeric,
character, inatkey, or occkey.}

\item{...}{curl options; named parameters passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
One or more in a list of both class inatkey and occkey
}
\description{
Coerce occurrence keys to iNaturalist id objects
}
\examples{
\dontrun{
spnames <- c('Accipiter striatus', 'Setophaga caerulescens',
  'Spinus tristis')
out <- occ(query=spnames, from='inat', limit=2)
res <- occ2df(out)
(tt <- as.inat(out))
(uu <- as.inat(res))
as.inat(res$key[1])
as.inat(as.list(res$key[1:2]))
as.inat(tt[[1]])
as.inat(uu[[1]])
as.inat(tt[1:2])
}
}
\seealso{
Other coercion: 
\code{\link{as.ala}()},
\code{\link{as.bison}()},
\code{\link{as.gbif}()},
\code{\link{as.idigbio}()},
\code{\link{as.obis}()},
\code{\link{as.vertnet}()}
}
\concept{coercion}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ2df.R
\name{occ2df}
\alias{occ2df}
\title{Combine results from occ calls to a single data.frame}
\usage{
occ2df(obj, what = "data")
}
\arguments{
\item{obj}{Input from occ, an object of class \code{occdat}, or an object
of class \code{occdatind}, the individual objects from each source within the
\code{occdat} class.}

\item{what}{(character) One of data (default) or all (with metadata)}
}
\description{
Combine results from occ calls to a single data.frame
}
\details{
This function combines a subset of data from each data provider to a single
data.frame, or metadata plus data if you request \code{what="all"}. The
single data.frame contains the following columns:
\itemize{
\item name - scientific (or common) name
\item longitude - decimal degree longitude
\item latitude - decimal degree latitude
\item prov - data provider
\item date - occurrence record date
\item key - occurrence record key
}
}
\examples{
\dontrun{
# combine results from output of an occ() call
spnames <- c('Accipiter striatus', 'Setophaga caerulescens',
  'Spinus tristis')
out <- occ(query=spnames, from='gbif', gbifopts=list(hasCoordinate=TRUE),
  limit=10)
occ2df(out)
occ2df(out$gbif)

out <- occ(
  query='Accipiter striatus',
  from=c('gbif','bison','ebird','inat'),
  gbifopts=list(hasCoordinate=TRUE), limit=2)
occ2df(out)
occ2df(out$bison)

# or combine many results from a single data source
spnames <- c('Accipiter striatus', 'Spinus tristis')
out <- occ(query=spnames, from='gbif', limit=2)
occ2df(out$gbif)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_coverage.R
\name{occ_coverage}
\alias{occ_coverage}
\title{Automatically generate coverages for a spocc search}
\usage{
occ_coverage(occObj, coverage = "all")
}
\arguments{
\item{occObj}{an search object returned by occ}

\item{coverage}{a vector of coverage types to generate. These include
'temporal','spatial','taxa', or just 'all'.}
}
\description{
This function will automatically generate metadata for spocc
queries that can then be converted to other standards.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.vertnet.R
\name{as.vertnet}
\alias{as.vertnet}
\title{Coerce occurrence keys to vertnetkey/occkey objects}
\usage{
as.vertnet(x)
}
\arguments{
\item{x}{Various inputs, including the output from a call to \code{\link[=occ]{occ()}}
(class occdat), \code{\link[=occ2df]{occ2df()}} (class data.frame), or a list, numeric,
character, vertnetkey, or occkey.}
}
\value{
One or more in a list of both class vertnetkey and occkey
}
\description{
Coerce occurrence keys to vertnetkey/occkey objects
}
\details{
Internally, we use \code{\link[rvertnet:vert_id]{rvertnet::vert_id()}}, whereas \code{\link[=occ]{occ()}}
uses \code{\link[rvertnet:vertsearch]{rvertnet::vertsearch()}}.
}
\examples{
\dontrun{
# spnames <- c('Accipiter striatus', 'Setophaga caerulescens',
#   'Spinus tristis')
# out <- occ(query=spnames, from='vertnet', has_coords=TRUE, limit=2)
# res <- occ2df(out)
# (tt <- as.vertnet(out))
# (uu <- as.vertnet(res))
# keys <- Filter(Negate(is.na), res$key)
# as.vertnet(keys[1])
# as.vertnet(as.list(keys[1:2]))
# as.vertnet(tt[[1]])
# as.vertnet(uu[[1]])
# as.vertnet(tt[1:2])
}
}
\seealso{
Other coercion: 
\code{\link{as.ala}()},
\code{\link{as.bison}()},
\code{\link{as.gbif}()},
\code{\link{as.idigbio}()},
\code{\link{as.inat}()},
\code{\link{as.obis}()}
}
\concept{coercion}
