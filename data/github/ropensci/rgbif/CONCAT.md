<!-- README.md is generated from README.Rmd. Please edit that file and knit -->



# rgbif <img src="man/figures/logo.png" align="right" alt="" width="120">

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rgbif)](https://cranchecks.info/pkgs/rgbif)
[![R-CMD-check](https://github.com/ropensci/rgbif/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rgbif/actions?query=workflow%3AR-CMD-check)
[![real-requests](https://github.com/ropensci/rgbif/workflows/R-check-real-requests/badge.svg)](https://github.com/ropensci/rgbif/actions?query=workflow%3AR-check-real-requests)
[![codecov.io](https://codecov.io/github/ropensci/rgbif/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rgbif?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rgbif)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rgbif)](https://cran.r-project.org/package=rgbif)
[![DOI](https://zenodo.org/badge/2273724.svg)](https://zenodo.org/badge/latestdoi/2273724)

`rgbif` gives you access to data from [GBIF][] via their REST API. GBIF versions their API - we are currently using `v1` of their API. You can no longer use their old API in this package - see `?rgbif-defunct`.

Please cite rgbif. Run the following to get the appropriate citation for the version you're using:

```r
citation(package = "rgbif")
```

To get started, see:

* rgbif vignette (https://docs.ropensci.org/rgbif/articles/rgbif.html): an introduction to the package's main functionalities.
* Function reference (https://docs.ropensci.org/rgbif/reference/index.html): an overview of all `rgbif` functions.
* Articles (https://docs.ropensci.org/rgbif/articles/index.html): vignettes/tutorials on how to download data, clean data, and work with taxonomic names.

Check out the `rgbif` [paper][] for more information on this package and the sister [Python][pygbif], [Ruby][gbifrb], and [PHP][phpgbif] clients.

Note: Maximum number of records you can get with `occ_search()` and `occ_data()` is 100,000. See https://www.gbif.org/developer/occurrence

## Installation


```r
install.packages("rgbif")
```

Or, install development version


```r
pak::pkg_install("ropensci/rgbif")
# OR
install.packages("rgbif", repos="https://dev.ropensci.org")
```


```r
library("rgbif")
```

## Screencast

<a href="https://vimeo.com/127119010"><img src="man/figures/README-screencast.png" width="400"></a>

## Contributors

This list honors all contributors in alphabetical order. Code contributors are in bold.

[adamdsmith](https://github.com/adamdsmith) - [AgustinCamacho](https://github.com/AgustinCamacho) - [AldoCompagnoni](https://github.com/AldoCompagnoni) - [AlexPeap](https://github.com/AlexPeap) - [andzandz11](https://github.com/andzandz11) - [AshleyWoods](https://github.com/AshleyWoods) - [AugustT](https://github.com/AugustT) - [barthoekstra](https://github.com/barthoekstra) - **[benmarwick](https://github.com/benmarwick)** - [cathynewman](https://github.com/cathynewman) - [cboettig](https://github.com/cboettig) - [coyotree](https://github.com/coyotree) - **[damianooldoni](https://github.com/damianooldoni)** - [dandaman](https://github.com/dandaman) - [djokester](https://github.com/djokester) - [dlebauer](https://github.com/dlebauer) - **[dmcglinn](https://github.com/dmcglinn)** - [dnoesgaard](https://github.com/dnoesgaard) - [DupontCai](https://github.com/DupontCai) - [ecology-data-science](https://github.com/ecology-data-science) - [EDiLD](https://github.com/EDiLD) - [elgabbas](https://github.com/elgabbas) - [emhart](https://github.com/emhart) - [fxi](https://github.com/fxi) - [ghost](https://github.com/ghost) - [gkburada](https://github.com/gkburada) - [hadley](https://github.com/hadley) - [Huasheng12306](https://github.com/Huasheng12306) - [ibartomeus](https://github.com/ibartomeus) - **[JanLauGe](https://github.com/JanLauGe)** - **[jarioksa](https://github.com/jarioksa)** - **[jeroen](https://github.com/jeroen)** - **[jhnwllr](https://github.com/jhnwllr)** - [jhpoelen](https://github.com/jhpoelen) - [jivelasquezt](https://github.com/jivelasquezt) - [jkmccarthy](https://github.com/jkmccarthy) - **[johnbaums](https://github.com/johnbaums)** - [jtgiermakowski](https://github.com/jtgiermakowski) - [jwhalennds](https://github.com/jwhalennds) - **[karthik](https://github.com/karthik)** - [kgturner](https://github.com/kgturner) - [Kim1801](https://github.com/Kim1801) - [ljuliusson](https://github.com/ljuliusson) - [ljvillanueva](https://github.com/ljvillanueva) - [luisDVA](https://github.com/luisDVA) - [martinpfannkuchen](https://github.com/martinpfannkuchen) - [MattBlissett](https://github.com/MattBlissett) - [MattOates](https://github.com/MattOates) - [maxhenschell](https://github.com/maxhenschell) - **[mdsumner](https://github.com/mdsumner)** - [no-la-ngo](https://github.com/no-la-ngo) - [Octoberweather](https://github.com/Octoberweather) - [Pakillo](https://github.com/Pakillo) - **[peterdesmet](https://github.com/peterdesmet)** - [PhillRob](https://github.com/PhillRob) - [poldham](https://github.com/poldham) - [qgroom](https://github.com/qgroom) - [raymondben](https://github.com/raymondben) - [rossmounce](https://github.com/rossmounce) - [sacrevert](https://github.com/sacrevert) - [sagitaninta](https://github.com/sagitaninta) - **[sckott](https://github.com/sckott)** - [scottsfarley93](https://github.com/scottsfarley93) - [simon-tarr](https://github.com/simon-tarr) - **[SriramRamesh](https://github.com/SriramRamesh)** - [stevenpbachman](https://github.com/stevenpbachman) - [stevensotelo](https://github.com/stevensotelo) - **[stevenysw](https://github.com/stevenysw)** - [TomaszSuchan](https://github.com/TomaszSuchan) - [tphilippi](https://github.com/tphilippi) - [vandit15](https://github.com/vandit15) - [vervis](https://github.com/vervis) - **[vijaybarve](https://github.com/vijaybarve)** - [willgearty](https://github.com/willgearty) - [Xuletajr](https://github.com/Xuletajr) - [yvanlebras](https://github.com/yvanlebras) - [zixuan75](https://github.com/zixuan75)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rgbif/issues).
* License: MIT
* Get citation information for `rgbif` in R doing `citation(package = 'rgbif')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

- - -

This package is part of [spocc](https://github.com/ropensci/spocc), along with several other packages, that provide access to occurrence records from multiple data sources.

- - -

[mapr]: https://github.com/ropensci/mapr
[paper]: https://peerj.com/preprints/3304/
[GBIF]: https://www.gbif.org/
[pygbif]: https://github.com/sckott/pygbif
[gbifrb]: https://github.com/sckott/gbifrb
[phpgbif]: https://gitlab.res-telae.cat/restelae/php-gbif
rgbif 3.6.0
===========

### Downloads

* typo in download predicate functions fixed - `mulitpoint` -> `multipoint` (#460) thanks @damianooldoni for catching that
* added three new predicate keys: `stateProvince` (#458), `gadm` (#462), and `occurrenceStatus` (#465)

### MINOR IMPROVEMENTS

* add two new occurrence issues: `FOOTPRINT_SRS_INVALID` and `FOOTPRINT_WKT_INVALID` (#454)
* `occ_download_import()` docs: more information on `data.table::fread` parameters and particular ones that would be useful to sort out data read issues (#461)

### BUG FIXES

* fix `occ_download_get()`: downloaded files used to have a certain content type in response header we checked for, but its changed at least once even in successful responses, so that step has been removed (#464)
* fix `occ_download_import()`: country code for Namibia is `NA` - this was turning into the R missing value `NA` - now fixed (#463)


rgbif 3.5.2
===========

### Download predicates

* in occurrence download predicate builder checks, to better help users, give the name of the key that fails upon failure instead of just the string 'key' (#450)
* occurrence download predicates: new key `coordinateUncertaintyInMeters` added, e.g. usage: `pred_lt("coordinateUncertaintyInMeters",10000)`  (#449)
* `pred_and()` and `pred_or()` slight change: now required that more than one predicate is passed to each of these functions because it doesn't make sense to do an `and` or `or` predicate with only one predicate (#452)
* fix for use of `pred_not(pred_notnull())` (#452)

### MINOR IMPROVEMENTS

* add a new occurrence issue (`TAXON_MATCH_AGGREGATE`) and a new name issue (`BACKBONE_MATCH_AGGREGATE`) (#453)

### BUG FIXES

* remove geoaxe references in man-roxygen template doc files - not using pkg anymore here and that pkg is cran archived too (#448)


rgbif 3.5.0
===========

### MINOR IMPROVEMENTS

* remove package wicket - use package wellknown instead - no user facing changes related to this (#447)
* remove package geoaxe (to be archived on CRAN soon) - use package sf instead (#447)

### BUG FIXES

* fix to download predicate function `pred_not()`: it was not constructing the query correctly, fixed now. user facing change as well: it now expects a predicate to be passed, and only a single predicate as GBIF not predicate only accepts one predicate (#446)


rgbif 3.4.2
===========

### MINOR IMPROVEMENTS

* Add new occurrence issue `DIFFERENT_OWNER_INSTITUTION` (#444)
* re-record all test fixtures

### BUG FIXES

* fix bug in `occ_search()` (#443)


rgbif 3.4
=========

### MINOR IMPROVEMENTS

* Documentation: clarify for `occ_search()` and `occ_data()` what parameters accept many values and which do not; in addition, we clarify which parameters accept multiple values in the same HTTP request, and those that accept multiple values but apply each in separate HTTP requests. See also `?many-values` manual file  (#369)
* `gbif_issues()` gains 9 new occurrence issues (#435)
* for `occ_search()` and `occ_data()`, `basisOfRecord` parameter now supports multiple values, both in one request and in different requests, depending on input format (see "Multiple values passed to a parameter" section in `?occ_search`)  (#437)
* remove vignettes from cran to avoid cran checks - still available on our docs site (#438)
* `occ_download_get()`: GBIF slightly altered download behavior - we now explicitly follow any redirects to get a download (#439)
* `print.occ_download_meta` (used when you run `occ_download_meta()`) was printing `NA` for number of results found if no results were ready yet - now prints `0` instead of `NA` (#440)

### BUG FIXES

* `count_facet()` fixes: fixed internal fxn for `count_facet` for parsing results, was dropping values for facets; added assertions to check parameter types input by user for the fxn; changed so that keys and basisofrecord can be passed together (#436)


rgbif 3.3
=========

### MINOR IMPROVEMENTS

* added two new occurrence issues to `gbif_issues()`: `GEOREFERENCED_DATE_INVALID` and `GEOREFERENCED_DATE_UNLIKELY` (#430)

### BUG FIXES

* fixed an error in `occ_data()` caused by GBIF adding a new field of data to the output of `/occurrence/search/`: gadm. cleaned up internals of `occ_data()` to drop gadm, and other fields that are complex and take time to parse (use `occ_search()` if you want all the data fields)  (#427)
* `gbif_names()` fix: was ending up with invalid URLs to GBIF species pages because we had taxon keys with leading spaces somehow. now all leading and trailing spaces in taxon keys removed before making URLs  (#429)


rgbif 3.2
=========

### MINOR IMPROVEMENTS

* `gbif_issues()` changes: three new occurrence issues added; one name issue removed that's deprecated (#423)
* `gbif_citation()` rights field was empty unless pulling from a downloaded file; now fill in with `license` key; also a fix for when occurrence key passed to the function (#424)
* `establishmentMeans` now supported in `occ_download`/`pred` (#420)

### BUG FIXES

* fix for `occ_download_get()`: response content-type header changed recently, fixed (#422)


rgbif 3.1
=========

### MINOR IMPROVEMENTS

* finally delete code originally extracted from `plyr::rbind.fill` - use `data.table::rbindlist` in all cases (#417)
* fix failing test on cran for `dataset_search()` (#418)
* fix xd refs note on cran (non-file package anchored links) for curl pkg function (#419)

### BUG FIXES

* `occ_download_cancel_staged()` fix: was broken cause we were indexing to a column in a table with `[,"key"]` (#416)


rgbif 3.0
=========

### BREAKING CHANGE

* Many functions (`occ_search`, `occ_get`, `name_usage`, `name_lookup`, `name_suggest`, `name_backbone`, and `dataset_search`) have a `return` parameter to toggle what is returned from the function call. To simplify rgbif maintainence, we've deprecated the `return` parameter. We've left it in each of the functions, but it no longer does anything, other than raising a warning if used. This means that function calls to these functions now always return the same data structure, making it easier to reason about for the user, as well as for us developers trying to make sure the package works as expected under a variety of conditions. If you have been using the `return` parameter, do the same function call as before, but now index to the output you need. This is a breaking change, thus the major version bump  (#413)

### NEW FEATURES

* new function `occ_download_cached()`, which takes the same input as `occ_download()`, but instead of starting a query, it checks if you've recently made the same request (with configureable settings for what "recent" means). This can save time when you're doing occurrence download requests that you may have done in the recent past (#308)

### MINOR IMPROVEMENTS

* configured package to be able to use two different base urls, `api.gbif-uat.org` and `api.gbif.org`. We have only used the latter previously, but now can configure rgbif to use the former, mostly for testing purposes  (#398)
* `occ_download_import()` gains `encoding` parameter that is passed down to `data.table::fread` to make it very clear that encoding can be configured (even though you could have before via `...`) (#414)

### BUG FIXES

* fix tibble construction (#412)


rgbif 2.3
=========

### MINOR IMPROVEMENTS

* max records you can return for `/occurrence/search` route is now 100,000 (used in `occ_data()` and `occ_search()`). updated docs throughout accordingly (#405)
* improved docs in `occ_download_queue()` for how we determine when a job is done. see new section "When is a job done?" (#409)
* print methods `print.occ_download_prep` and `print.occ_download` improved. previously well-known text strings were printed in their entirety. now they are handled to only print so many characters; also applies to any download predicate string that's long (#407)
* `occ_download_get()` now supports using a progress bar by passing in `httr::progress()` (#402)
* `occ_data()` and `occ_search()` gain two new parameters: `recordedByID` and `identifiedByID` (#403)

### BUG FIXES

* fix in `occ_download_queue`: an empty `occ_download_meta()` lead to problems; now removing any `NULL`'s from a list of `occ_download_meta()` outputs before further work (#408)
* fix in `occ_download_queue`: we were not accounting for job status "cancelled" (#409)
* `occ_download_import()` fix: `fill` parameter was set to `TRUE` by default, changed to `FALSE`. improved docs for this fxn on passing down parameters to `data.table::fread` (#404)


rgbif 2.2
=========

### MINOR IMPROVEMENTS

* add a section _Download status_ to the `?downloads` manual file listing all the different download status states a download can have and what they mean (#390)
* fix `gbif_issues`/`gbif_issues_lookup`: added four missing occurrence issues to the package (COORDINATE_PRECISION_INVALID, COORDINATE_UNCERTAINTY_METERS_INVALID, INDIVIDUAL_COUNT_INVALID, and INTERPRETATION_ERROR) (#400)
* doing real tests now for `occ_download()` via vcr (#396)

### BUG FIXES

* fix `name_lookup()`: we were attempting to rearrange columns when no results found, leading to an error (#399)


rgbif 2.1
=========

### DEFUNCT

* the `spellCheck` parameter has been removed from the occurrence routes; thus, the `occ_spellcheck()` function is now defunct - and the parameter `spellCheck` has been removed from `occ_data()` and `occ_search()` (#397)

### MINOR IMPROVEMENTS

* docs fix for `occ_data()`: remove `...` parameter definition as it wasn't used in the function (#394)

### BUG FIXES

* download predicate fxns fix: "within" wasnt being handled properly (#393) thanks @damianooldoni


rgbif 2.0
=========

### NEW FEATURES

* The download query user interface for `occ_download()` has changed in a breaking fashion (thus the major version bump). After installation, see `?download_predicate_dsl`. Much more complex queries are now possible with `occ_download()`. TL;DR: you now construct queries with functions like `pred("taxonKey", 3119195)` rather than passing in strings like `taxonKey = 3119195`, and `pred_gt("elevation", 5000)` instead of `"elevation > 5000"`  (#362)
* gains new function `occ_download_wait()` to re-run `occ_download_meta()` until the download is ready - kinda like `occ_download_queue()` but for a single download (#389)
* `occ_download_dataset_activity()` gains pagination parameters `limit` and `start` to paginate through results (#382)
* `gbif_citation()` now works with the output of `occ_data()` in addition to the other existing inputs it accepts (#392)

### MINOR IMPROVEMENTS

* typo fix in the _geometry_ section of the `occ_download()` manual file (#387)
* vignettes fixes (#391)

### BUG FIXES

* `gbif_citation()` tests needed preserve body bytes for vcr (#384)
* fix to `occ_count()` and `count_facet()`: isGeoreferenced/georeferenced variable needed booleans converted to lowercase before being sent to GBIF (#385) (#386)


rgbif 1.4.0
===========

### NEW FEATURES

* gains new function `mvt_fetch()` for fetching Map Vector Tiles (MVT). mvt used to be an option in `map_fetch()`, but we only returned raw bytes for that option. With `mvt_fetch()` we now leverage the `protolite` package, which parses MVT files, to give back an sf object (#373) thanks to @jeroen for the protolite work to make this work
* associated with above, `map_fetch()` loses the `format = ".mvt"` option; and thus now only returns a `RasterLayer`
* `occ_issues()` and `name_issues()` reworked. Both now use the same underlying internal logic, with occ_issues pulling metadata specfic to occurrence issues and name_issues pulling metadata specific to name issues. name_issues used to only be a data.frame of name issues, but can now be used similarly to occ_issues; you can pass the output of `name_usage()` to name_issues to filter/parse name results by their associated name issues. Associated with this, new function `gbif_issues_lookup` can be used to lookup either occurrence or name issues by their full name or code (#363) (#364)

### MINOR IMPROVEMENTS

* fix examples and tests that had WKT in the wrong winding order (#361)
* parsing GBIF issues in the output of `name_usage()` wasn't working (#328) (#363) (#364)
* `name_lookup()` gains an additional parameters `issue` for filtering name results by name issues (#335) (#363) (#364)
* fixed definitions of `x`, `y`, `z` parameters in `map_fetch()` manual file (#375)
* added examples to `gbif_citation()` manual file for accessing many citations (#379)
* fixed a test for `occ_download_queue()` (#365)
* `name_*` function outpus have changed, so be aware if you're using those functions

### BUG FIXES

* fixed issue with `map_fetch()`: when srs was `EPSG:3857`, the extent we set was incorrectly set as `raster::extent(-180, 180, -85.1, 85.1)`. Now the extent is `raster::extent(-20037508, 20037508, -20037508, 20037508` (#366) (#367) thanks @dmcglinn for reporting and @mdsumner for fixing!
* fix for Windows platforms for `gbif_citation()` for `occ_download_get` objects. we weren't correctly creating the path to a file on windows (#359)
* fix to `print.gbif_data` (#370) (#371)
* `occ_download()` was erroring with a useless error when users try to use the fxn with the same parameter input types as `occ_search`/`occ_data`; when this happens now there is a useful error message (#381)
* fix to `occ_download()`: when `type = "in"` was used, we weren't creating the JSON correctly, fixed now  (#362)


rgbif 1.3.0
===========

### NEW FEATURES

* `occ_download()` and `occ_download_prep()` gain a new parameter `format` for specifying the type of download. options are DWCA (default), SIMPLE_CSV, or SPECIES_LIST. SIMPLE_CSV and SPECIES_LIST are csv formats, while DWCA is the darwin core format (#352)
* now throughout the package you can pass `NA` in addition to `NULL` for a missing parameter - both are removed before being sent to GBIF (#351)

### MINOR IMPROVEMENTS

* replace `tibble::as_data_frame`/`tibble::data_frame` with `tibble::as_tibble` (#350)
* `key` and `gbifID` in the output of `occ_data`/`occ_search`/`occ_get` have been changed so that both are character class (strings) to match how GBIF encodes them (#349)
* fix some test fixtures to use preserve exact bytes so that cran checks on debian clang devel don't fail (#355)

### BUG FIXES

* fix to `occ_download`: fail with useful message when user does not pass in queries as character class (#347)
* fix to `occ_download`: fail with useful message now when user/pwd/email not found or given (#348)


rgbif 1.2.0
===========

### NEW FEATURES

* pkgdown documentation site (#336) (#337) all work done by @peterdesmet
* package gains hex logo (#331) (#332) thanks @peterdesmet
* big change to `elevation()` function: the Google Maps API requires a form of payment up front, and so we've decided to move away from the service. `elevation()` now uses the Geonames service <https://www.geonames.org/>; it does require you to register to get a username, but its a free service. Geonames has a few different data models for elevation and can be chosen in the `elevation_model` parameter (#344) (#345)
* biggish change to `occ_data()`/`occ_search()` output: the data.frame in the `data` slot now always has the first column as the occurrence key (`key`), and the second column is now the scientific name (`scientificName`). the previously used `name` column still exists in the data.frame, so as not to break any user code, but is simply a duplicate of the `scientificName` column. in a future version of this package the `name` column will be dropped (#329)

### MINOR IMPROVEMENTS

* README gains full list of code contributors and any folks involved in github issues (#339) (#343) thanks @peterdesmet
* update pkg citation, include all authors (#338)
* added more to `occ_search()`/`occ_data()`/`occ_download()` documentation on WKT (well-known text) with respect to winding order. GBIF requires counter-clockwise winding order; if you submit clockwise winding order WKT to `occ_search()` or `occ_data()` you should get data back but the WKT is treated as an exclusion, so returns data outside of that shape instead of within it; if you submit clockwise winding order WKT to `occ_download()` you will get no data back (#340)

### BUG FIXES

* fix bug in `occ_download()`, was failing in certain cases because of some bad code in an internal function `catch_err()` (#333)
* `occ_download()` was not returning user name and email in it's print method (#334)
* `occ_issues()` was failing with `occ_data()` or `occ_search()` input when `type="many"` (i.e., when > 1 thing was passed in) (#341)


rgbif 1.1.0
===========

### NEW FEATURES

* tests that make HTTP requests are now cached via the `vcr` package so do not require an internet connection (#306) (#327)
* added name usage issues (similar to occurrence issues) data. in part fixes `name_usage()` problem, more work coming to allow users to use the name issues data like we allow for occurrence issues through `occ_issues()` (#324) 

### MINOR IMPROVEMENTS

* `map_fetch()` changes following changes in GBIF maps API: new parameters `taxonKey`, `datasetkey`, `country`, `publishingOrg`, `publishingCountry` and removed parameters `search` and `id`; note that this changes how queries work with this function (#319)
* added note to `map_fetch()` docs that `style` parameter does not necessarily use the style you give it. not sure why (#302)
* fixed messaging in `occ_download_queue()` to report an accurate number of jobs being processed; before we were just saying "kicking off first 3 requests" even if there were only 1 or 2 (#312)

### BUG FIXES

* fix to `occ_get()` when `verbatim=TRUE` (#318)
* `elevation()` function now fails better. when the API key was invalid the function did not give an informative message; now it does (#322)


rgbif 1.0.2
===========

### MINOR IMPROVEMENTS

* significant change to `occ_download_queue()`: sleep time between successive calls to check on the status of download requests is now 10 seconds or greater. This shouldn't slow down your use of `occ_download_queue()` much because most requests should take more than the 10 seconds to be prepared (#313)
* add tests for download queue method (#315)
* explicitly `@importFrom` fxns used from `lazyeval` package to avoid check note (#316)
* remove `reshape2` and `maps` packages from Suggests (#317)

### BUG FIXES

* fix bug in `name_usage()`: we were screwing up parsing of issues column when single taxon keys passed in (#314)

rgbif 1.0.0
===========

### NEW FEATURES

* `occ_issues()` now works with download data and arbitrary data.frame's (#193)
* New downloads queueing tools: gains functions `occ_download_prep()` for preparing a download request without executing it, and `occ_download_queue()`  for kicking off many download jobs while respecting GBIF's downloads rate limits. See also internal R6 classes for dealing with queuing: `DownReq`, `GifQueue`. See `?occ_download_queue` to get started (#266) (#305) (#311)
* New function `map_fetch()` working with the GBIF maps API <https://www.gbif.org/developer/maps>. See `?map_fetch` to get started (#238) (#269) (#284) thanks to @JanLauGe for the work on this
* `name_lookup()` gains `origin` parameter (#288) (#293) thanks @peterdesmet and @damianooldoni
* `name_lookup()` and `name_usage()` gain internal paging - just as `occ_search()`/`occ_data()` have (#291) (see also #281) thanks @damianooldoni 
* new import `lazyeval`, and new suggests `png` and `raster`
* `occ_search()`/`occ_data()` gain parameter `skip_validate` (boolean) to skip or not stkip WKT validation by the `wicket` package

### MINOR IMPROVEMENTS

* removed warnings about parameters that were removed in previous versions of the package (#189)
* add citation file (#189)
* updated `name_usage()` to check params that now only allow 1 value: name, language, datasetKey, rank (#287)
* `occ_count()` loses `nubKey`, `catalogNumber`, and `hostCountry` as those parameters are no longer accepted by GBIF

### BUG FIXES

* fixed bug in `name_usage()`, was screwing something up internally (#286)
* fixed bug in `occ_data()`: curl options weren't being passed through (#297)
* fixed geometry usage in `occ_search()`/`occ_data()` - skipping the wicket validation and constructing WKT by hand from bounding box (if bounding box given) - the validation that wicket does isn't what GBIF wants (#303)
* add `fill` parameter to  `occ_download_import()` to pass on to `fill` in `data.table::fread`, and set `fill=TRUE` as default.  (#292)
* better failure for `occ_download()` (#300)
* fix bug in `occ_download()` in which a single `taxonKey` passed in was failing (#283)
* `name_usage()` was ignoring `datasetKey` and `uuid` parameters (#290)

### DEFUNCT AND DEPRECATED

* `gbifmap()` has been removed, see the package `mapr` for similar functionality and `map_fetch()` in this package to use the GBIF map API (#298)


rgbif 0.9.9
===========

### NEW FEATURES

* Gains new functions `occ_download_datasets` and `occ_download_dataset_activity` to list datasets for a download,
and list the downloads activity of a dataset (#275) (#276)
* Gains a new vignette covering working with GBIF downloads 
in `rgbif` (#262)

### MINOR IMPROVEMENTS

* Guidance added to docs for downloads functions on length of the
request body (#263)
* Changed authentication details (user name, password, email) for 
downloads to allow any of the options: pass in as arguments,
store as R options, store as environment variables (#187)
* `gbif_citation()` function gains an S3 method for passing the 
output of `occ_download_meta()` to it. In addition, for downloads
`gbif_citation()` now returns a citation for the entire download 
(including) its DOI, in addition to citations for each dataset (#274) 
thanks @dnoesgaard

### BUG FIXES

* Fix documentation bug in `occ_count()`: `georeferenced` had a 
misleading description of what the value `FALSE` did (#265)
* Fixed bug in `gbifmap()` - was failing in some cases - better
error handlingn now (#271) thanks @TomaszSuchan
* Fixed `occ_download_cancel_staged()`: it wasn't passing on authentication
parameters correctly (#280)


rgbif 0.9.8
===========

### NEW FEATURES

* The GBIF API supports passing in many instances of the same
parameter for some parameters on some routes. Previously we
didn't support this feature, but now we do. See the
`?many-values` manual file for details. added docs to
individual functions that support this, and added additional
tests (#200) (#260) (#261)
* We've removed `V8` dependency and replaced with C++ based
WKT parser package `wicket`. We still use `rgeos` for some
WKT parsing. rgbif functions that use wicket: `gbif_bbox2wkt`,
`gbif_wkt2bbox`, `check_wkt` (#243)
* `httr` replaced with `crul` for HTTP reqeusts. As part of
this change, the `...` parameter was replaced in most functions
by `curlopts` which expects a list. Some functions require
a `...` parameter for facet inputs, so `...` is retained
with the addition of `curltops` parameter. A result of this
change is that whereas in the past parameters that were not
defined in a function that also had a `...` parameter
would essentially silently ignore that undefined parameter,
but with functions where `...` was removed a misspelled
or undefined parameter will cause an error with message (#256)

### MINOR IMPROVEMENTS

* moved to markdown docs (#258)
* namespacing calls to base R pkgs instead of importing them

### BUG FIXES

* Fixed problem in `occ_download_import()` to allow import
of csv type download in addition to darwin core archive.
additional change to `occ_download_get` to add `format`
attribute stating which format (#246)
* fix to `occ_download_import` adding `fill=TRUE` to
the `data.table::fread` call (#257)


rgbif 0.9.7
===========

### NEW FEATURES

* `occ_dowload` gains new parameter `body` to allow users to pass in
JSON or a list for the query instead of passing in statements to
`...`. See examples in `?occ_dowload`.

### MINOR IMPROVEMENTS

* Now using `tibble` for compact data.frame output for
`occ_download_import` instead of bespoke internal solution (#240)
* Moved all GBIF API requests to use `https` instead of `http` (#244)
* Improved print method for `occ_download_meta`

### BUG FIXES

* Fix to `occ_download` to structure query correctly when
`type=within` and `geometry` used because the structure is slightly
different than when not using `geometry` (#242)
* Fixed `occ_download` to allow `OR` queries for many values of a
parameter, e.g., `taxonKey=2475470,2480946` will be queried correctly
now as essentially `taxonKey=2475470` or `taxonKey=2480946` (#245)


rgbif 0.9.6
===========

### BUG FIXES

* Fixed a bug in `parsenames()` caused by some slots in the list
being `NULL` (#237)
* Fixed some failing tests: `occ_facet()` tests were failing due to
changes in GBIF API (#239)
* Fixes to `gbif_oai_get_records()` for slight changes in `oai`
dependency pkg (#236)


rgbif 0.9.5
===========

### NEW FEATURES

* `occ_search()` now has faceted search. This feature is not in `occ_data()`
as that function focuses on getting occurrence data quickly, so will not
do get facet data. This means that a new slot is available in the output
object from `occ_search()`, namely `facets`. Note that `rgbif` has had
faceted search for the species search route (`name_lookup()`) and the
registry search route (`dataset_search()`) for quite a while. (#215)
* new function (`occ_facet()`) to facilitate retrieving only
facet data, so no occurrence data is retrieved. (#215) (#229)
* A suite of new parameters added to `occ_search()` and
`occ_data()` following addition the GBIF search API: `subgenusKey`,
`repatriated`, `phylumKey`, `kingdomKey`,
`classKey`, `orderKey`, `familyKey`, `genusKey`, `establishmentMeans`,
`protocol`, `license`, `organismId`, `publishingOrg`, `stateProvince`,
`waterBody`, `locality` (#216) (#224)
* New parameter `spellCheck` added to `occ_search()` and
`occ_data()` that if `TRUE` spell checks anything passed to the `search`
parameter (same as `q` parameter on GBIF API; which is a full text
search) (#227)
* New function `occ_spellcheck` to spell check search terms, returns
`TRUE` if no spelling problems, or a list with info on suggestions
if not.
* Both `occ_search()` and `occ_data()` now have ability to support
queries where `limit=0`, which for one should be possible and not
fail as we did previously, and second, this makes it so that you
can do faceted searches (See above) and not have to wait for occurrence
records to be returned. (#222)
* `MULTIPOLYGON` well known text features now supported in the GBIF
API. Previously, you could not query `geometry` with more than
one polygon (`POLYGON`), but now you can. (#222)

### MINOR IMPROVEMENTS

* Improved docs for `occ_count()`, especially for the set of
allowed parameter options that the GBIF count API supports
* `occ_count()` gains new parameter `typeStatus` to indicate the
specimen type status.
* When no results found, the `data` slot now returns `NULL` instead
of a character string

### BUG FIXES

* Fixes to `gbif_photos()`: 1) Mapbox URLs to their JS and CSS assets
were out of date, and API key needed. 2) In RStudio, the `table` view
was outputting errors due to serving files on `localhost:<port>`
instead of simply opening the file; fixed now by checking platform
and using simple open file command appropriate for the OS. (#228) (#235)


rgbif 0.9.4
===========

### NEW FEATURES

* Now using `tibble` in most of the package when the output is
a data.frame (#204)
* New vignette _Taxonomic Names_ for discussing some common names
problems users may run into, and some strategies for dealing with
taxonomic names when using GBIF (#208) (#209)

### MINOR IMPROVEMENTS

* Replaced `is()` with `inherits()`, no longer importing `methods()` (#219)
* Improved docs for registry functions. Not all options were listed
for the `data` parameter, now they are (#210)
* Fixed documentation error in `gbifmap()` man file (#212) thanks to @rossmounce

### BUG FIXES

* Fixed bug in internal parser within `occ_download()`, in which
strings to parse were not being parsed correctly if spaces weren't in
the right place, should be more robust now, and added tests (#217). Came
from https://discuss.ropensci.org/t/rgbif-using-geometry-in-occ-download/395
* The parameter `type` was being silently ignored in a number of
registry functions. fixed that. (#211)


rgbif 0.9.3
===========

### NEW FEATURES

* `occ_data()` and `occ_search()` gain ability to more flexibly deal with inputs to the
`geometry` parameter. Previously, long WKT strings passed to `occ_search()` or
`occ_data()` would fail because URIs can only be so long. Another option is to use
the download API (see `?downloads`). This version adds the ability to choose what to
do with long WKT strings via the `geom_big` parameter: `asis` (same as previous version),
`bbox` which detects if a WKT sting is likely too long, and creates a bounding box from the
WKT string then once data is retrived, clips the result to the original WKT string; `axe`
uses the `geoaxe` package to chop up the input WKT polygon into many, with toggles in the
new parameters `geom_size` and `geom_n`. (#197) (#199)
* As part of this change, when >1 geometry value passed, or if `geom_big="axe"`, then
named elements of the output get names `geom1`, `geom2`, `geom3`, etc. instead of the
input WKT strings - this is because WKT strings can be very long, and make for very
awkward named access to elements. The original WKT strings can still be accessed via
`attr(result, "args")$geometry`

### MINOR IMPROVEMENTS

* code tidying throughout the package

### BUG FIXES

* Fix parsing bug in `name_usage()` function, see commit [e88cf01cc11cb238d44222346eaeff001c0c637e](https://github.com/ropensci/rgbif/commit/e88cf01cc11cb238d44222346eaeff001c0c637e)
* Fix to tests to use new `testthat` fxn names, e.g., `expect_gt()`
instead of `expect_more_than()`
* Fix to `occ_download()` to parse error correctly when empty body passed from
GBIF (#202)


rgbif 0.9.2
===========

### NEW FEATURES

* New function `occ_data()` - its primary purpose to perform faster data requests. Whereas
`occ_search()` gives you lots of data, including taxonomic hierarchies and media records,
`occ_data()` only gives occurrence data. (#190)

### MINOR IMPROVEMENTS

* Replaced `XML` with `xml2` (#192)
* Speed ups to the following functions due to use of `data.table::rbindlist()` for
fast list to data.frame coercion: `name_lookup()`, `name_backbone()`, `name_suggest()`,
`name_usage()`, and `parsenames()` (#191)
* Changes to `httr` usage to comply with changes in `httr >= v1.1.0`: now setting
encoding explicitly to `UTF-8` and parsing all data manually, using the internal
function `function(x) content(x, "text", encoding = "UTF-8")` (#195)

### BUG FIXES

* Fix to internal function `move_col()` to not fail on fields that don't exist.
Was failing sometimes when no latitude or longitude columns were returned. (#196)

rgbif 0.9.0
===============

### NEW FEATURES

* New set of functions (`gbif_oai_*()`) for working with
GBIF registry OAI-PMH service. Now importing `oai` package to
make working with GBIF's OAI-PMH service easier (#183)
* Added code of conduct (#180)
* Now sending user-agent header with all requests from this
package to GBIF's servers indicating what version of rgbif
and that it's an ropensci package. Looks like
`r-curl/0.9.4 httr/1.0.0 rOpenSci(rgbif/0.9.0)`, with whatever
versions of each package you're using. We also pass a user-agent
string with the header `X-USER-AGENT` in case the `useragent`
header gets stripped somewhere along the line (#185)
* New function `gbif_citation()` helps get citations for datasets
eith using the occurrence search API via `occ_search()` or the
downloads API via `occ_downlad()` (#178) (#179)

### MINOR IMPROVEMENTS

* Using `importFrom` instead of `import` in all cases now.
* Parameter `collectorName` changed to `recordedBy` (#184)

### BUG FIXES

* Fix to `occ_download_meta()` print method to handle 1 or more predicate results (#186)
* Fix to `occ_issues()` to work with `return=data` and `return=all` `occ_search()` output (#188)

rgbif 0.8.9
===============

### MINOR IMPROVEMENTS

* Updated `terraformer.js` javascript code included in the package
along with an update in that codebase (#156)
* The `email` parameter now `NULL` by default in the function
`occ_download()`, so that if not provided or not set in options,
then function fails. (#173)
* Additional explanation added to the `?downloads` help file.
* Added internal checks to `elevation()` to check for coordinates that
are impossible (e.g., latitude > 90), not complete (e.g., lat given,
long not given), or points at `0,0` (just warns, doesn't stop). (#176)
thanks @luisDVA
* General code tidying across package

### BUG FIXES

* A route changed for getting images for a taxon within the `/species`
route, fix to function `name_usage()` (#174)
* Fix to `occ_search()` to remove a block of code to do synonym checking.
This block of code was used if the parameter `scientificName` was passed,
and checked if the name given was a synonym; if yes, we used the accepted
name according to the GBIF backbone taxonomy; if no, we proceeded with the
name given by the user. We removed the block of code because the GBIF
API now essentially does this behind the scenes server side. See
https://github.com/gbif/gbif-api for examples. (#175)

rgbif 0.8.8
===============

### MINOR IMPROVEMENTS

* Additional tests added for `gbif_photos()` and `gbif_names()` (#170)

### BUG FIXES

* Fixed a few tests that were not passing on CRAN.

rgbif 0.8.6
===============

### NEW FEATURES

* New set of functions with names `occ_download*()` for working with the GBIF download API. This is the same service as using the GBIF website, but via an API. See `?downloads`. (#154) (#167)

### MINOR IMPROVEMENTS

* Explicitly import non-base R pkg functions, so importing from `utils`, `methods`, and `stats` (#166)

### BUG FIXES

* Fixed problem with `httr` `v1` where empty list not allowed to pass to
the `query` parameter in `GET` (#163)


rgbif 0.8.4
===============

### NEW FEATURES

* New functions for the `/enumerations` GBIF API route: `enumeration()`
and `enumeration_country()`. Many parts of the GBIF API make use of
enumerations, i.e. controlled vocabularies for specific topics - and are
available via these functions. (#152)

### IMPROVEMENTS

* `elevation()` now requires an API key (#148)
* The `V8` package an Import now, used to do WKT read/create with use of
the Javascript library Terraformer (http://terraformer.io/). Replaces
packages `sp` and `rgeos`, which are no longer imported (#155)
* Changed `occ_search()` parameter `spatialIssues` to
`hasGeospatialIssues` (#151)
* Added note to docs about difference between `/search` and `/count`
services, and how they work. (#150)
* Added tests for habitat parameter in `name_lookup()` (#149)
* Dropped `plyr` from Imports (#159)
* Dropped `stringr` from Imports (#160)
* Dropped `maps` and `grid` packages from Imports (#161)

### BUG FIXES

* Looping over records with `limit` and `start` parameters was in some
cases resulting in duplicate records returned. Problem fixed. (#157)

rgbif 0.8.0
===============

### IMPROVEMENTS

* All example moved to `\dontrun` (#139)
* README fixes for html (#141)
* Fixed documentation in `occ_search()` to give correct values for default and max limit
and start parameters (#145)
* Changed internal `GET` helper function to properly pass on error message (#144)
* Replaced `assertthat::assert_that()` with `stopifnot()` to have one less dependency (#134)
* Fixed `occ_search()` to allow ability to query by only publishingCountry, that is, with no
other parameters if desired (#137)

### BUG FIXES

* Fixed bug in internal `GET()` helper function to just pass `NULL` to the `query` parameter when
the list of length 0 passed, since it caused requests to fail in some cases.
* Fix to `name_lookup()` to force a logical entry for certain parameters - before this fix
if the correct logical param was not passed, the GBIF API went with its default parameter (#135)
* Fixed bug in `name_backbone()` due to change in `namelkupparser()` helper function - fixes
parsing for verbose output (#136)
* Fixed some broken URLs in `occ_search()` documentation (#140)


rgbif 0.7.7
===============

### NEW FEATURES

* New function `occ_issues()` to subset data from `occ_search()` based on GBIF issues. (#) (#122)
* Related to the last bullet, GBIF issues now are returned by default in `occ_search()` results, and are intentionally moved to the beginning of the column order of the data to be more obvious. (#102)
* `occ_search()` now returns all data fields by default. The default setting for the `fields` parameter is `all` - but can be changed. See `?occ_search`
* New function `gbif_names()` to view highlighted terms in name results from a call to `name_lookup()`. (#114)
* New functions: `occ_issues_lookup()` to lookup GBIF issues based on code name or full issue name, and `gbif_issues()` to print the entire issues table.

### IMPROVEMENTS

* Completely replaced `RCurl` with `httr`
* Completely replaced `RJSONIO` with `jsonlite`. Should see slight performance in JSON parsing with `jsonlite`.
* Default number of records in `occ_search()` now 500; was 25. (#113)
* Vignette for old version of GBIF API removed.
* New vignette for cleaning data via GBIF issues added. (#132)
* Functions for working with old GBIF API removed, now defunct. (#116)
* Now better parsing for some functions (`organizations()`, `datasets()`, `networks()`, `nodes()`, `installations()`) to data.frames when possible. (#117)
* Added further help to warn users when searching on ranges in latitude or longitude in `occ_search()` (#123)
* `callopts` parameter changed to `...` throughout all functions. Now pass on options to `httr` as named lists or functions. (#130)
* Beware that GBIF data is becoming Darwin Core compliant - so many parameters throughout this package have changed from sentence_case to camelCase.
* `dataset_search()` and `dataset_suggest()` gain new parameter `publishingOrg`
* Default for `limit` parameter changed to 100 for dataset functions: `dataset_search()`, `dataset_suggest()`, and `datasets()`.
* Default for `limit` parameter changed to 100 for registry functions: `installations()`, `networks()`, `organizations`, and `nodes()`.
* Parameter changes in `networks()`: `name`, `code`, `modifiedsince`, `startindex`, and `maxresults` gone; new parameters `query`, `identifier`, `identifierType`, `limit`, and `start`
* Parameter changes in `nodes()`: new parameters `identifier`, `identifierType`, `limit`, and `start`

### BUG FIXES

* `occ_search()` failed sometimes on species that were not found. Fixed. (#112)
* Added better handling of some server errors to pass on to user. (#115) (#118)
* Fixed incorrect parsing for some cases in `occ_search()` (#119)
* Fixed bad parsing on output from `name_lookup()` (#120)
* Fixed single map option in `gbif_photos()` that caused map with no data. (#121)
* Fixed some parameter names in `name_()` functions according to changes in the GBIF API spec, and fixed documentation to align with GBIF API changes, and added note about maximum limit. (#124) (#127) (#129) Thanks to @willgearty !
* Fixed internals of `occ_search()` so that user can pass in multiple values to the `issue` parameter. (#107)
* Fixed URL to tutorial on ropensci website (#105) Thanks @fxi !

rgbif 0.7.0
===============

### NEW FEATURES

* `occ_search()` now has a `dplyr` like summary output when `return='all'`. See `?occ_search` for examples. You can still easily access all data, by indexing to `meta`, `hierarchy`, `data`, or `media` via e.g., `$data`, `['data']`, or `[['data']]`. (#95)
* Media now returned from the GBIF API. Thus, in `occ_search()`, we now return a media slot in the output list by default.
* New function `gbif_photos()` to view media files (photos in the wild or of museum specimens). Two options are available, `which='map'` creates a single map which presents the image when the user clicks on the point, and `which='table'` in which a table has one row for each image, presenting the image and an interactive map with the single point. (#88)
* Two new packages are imported: `sp` and `whisker`

### IMPROVEMENTS

* GBIF updated their API, now at v1. URL endpoints in `rgbif` changed accordingly. (#92)
* GBIF switched to using 2-letter country codes. Take note. (#90)
* GBIF switched all parameters to `camelCase` from `under_score` style - changed accordingly in `rgbif`.
* Using package custom version of `plyr::compact()` instead of importing from `plyr`.
* In `name_lookup()` removed `facet_only` parameter as it doesn't do anything - use `limit=0` instead. Further, added two new slots of output: `hierarchy` and `names` (for common/vernacular names) (#96). The output can be determined by user via the `return` parameter.
* In `name_suggest()`, if the field `higherClassificationMap` is selected to be returned via the `fields` parameter, a list is returned with a data frame, and a list of the hierarchies separately. If `higherClassificationMap` is not selected, only a data frame is returned.
* `occ_search()` gains new parameters  `mediatype` and `issue` (#93), with detailed list of possible options for the `issue` parameter. Gains new examples for searching for images, examples of calls that will throw errors.
* Updated the vignette.

### BUG FIXES

* Added better error message to `check_wkt()`.
* `facet_only` parameter removed from `dataset_search()` function as it doesn't do anything - use `limit=0` instead.
* Fixed some examples that didn't work correctly.

rgbif 0.6.3
===============

### IMPROVEMENTS

* Added functions `gbif_bbox2wkt()` and `gbif_wkt2bbox()` to convert a bounding box to wkt and a wkt object to a bounding box, respectively. Copied from the `spocc` package. Prefixes to fxn names will avoid conflicts.
* Now spitting out more informative error messages when WKT strings passed in are not properly formed, either from `rgeos::readWKT` or from the returned response from GBIF.

rgbif 0.6.2
===============

### BUG FIXES

* `gbifmap()` was throwing an error because it was looking for two variables `latitude` and `longitude`, which had been changed to `decimalLatitude` and `decimalLongitude`, respectively, in other functions in this package. Fixed. (#81)
* `occ_get()` was updated to include changes in the GBIF API for this endpoint. The fix included fixing the parser for verbatim results, see `rgbif::gbifparser_verbatim`. (#83)
* Fixed bugs in `elevation()` - it was expecting column names to be latitude and longitude, whereas inputs from other `rgbif` functions have changed to decimalLatitude and decimalLongitude.
* Fixed bug in `count_facet()` introduced b/c GBIF no longer accepts hostCountry or nubKey parameters.

### IMPROVEMENTS

* `gist()`, `stylegeojson()`, and `togeojson()` functions now listed as deprecated. Their functionality moved to the `spocc` package (http://cran.r-project.org/web/packages/spocc/index.html). These functions will be removed from this package in a future version. (#82)
* Added a quick sanity test for `gbifmap()`.
* Added tests for `occ_get()` for when `verbatim=TRUE`, which gives back different data than when `verbatim=FALSE`.

rgbif 0.6.0
===============

### BUG FIXES

* A number of variables changed names to better follow the Darwin Core standard. `latitude` is now `decimalLatitude`. `longitude` is now `decimalLongitude`. `clazz` is now `class`. Code in this package changed to accomodate these changes. `date` is now `eventDate`. `georeferenced` is now `hasCoordinate`. Beware of these changes in your own code using `rgbif` - find and replace for these should be easy.
* Changed `altitude` parameter in `occ_search()` to `elevation` - should have been `elevation` the whole time.
* `occ_count()` function with parameter changes: `nubKey` parameter in changed to `taxonKey`. New parameter `protocol`. Parameter `catalogNumber` gone. Parameter `hostCountry` gone. These parameters are still in the function definition, but if called they throw a useful warning telling you the correct parameter names. (#76)
* Fixed bug in `name_lookup()` function that was labeling facet outputs incorrectly. (#77)

### IMPROVEMENTS

* Better checking and parsing of response data from GBIF: Across all functions, we now check that the response content type is `application/json`, then parse JSON ourselves using `RJSONIO::fromJSON` (instead of httr doing it).
* Across all functions, we now return all potenital character class columns as character class (instead of factor), by passing `stringsAsFactors = FALSE` to all `data.frame()` calls.
* Now using assertthat package in various places to give better error messages when the wrong input is passed to a function.
* Four parameters have name changes in the `occ_search()` function. These parameters are still in the function definition, but if called they throw a useful warning telling you the correct parameter names. (#75)
* Updated docs in `name_usage`, `name_backbone`, `name_lookup`, and `name_suggest` functions.
* `sourceId` parameter in `name_usage()` function doesn't work so error message is thrown when used.

### NEW FEATURES

* New function `check_wkt()` to check that well known text string is the right format. (#68)
* New dataset typestatus to look up possible specimen typeStatus values. See #74 for more information.
* GBIF added some new parameters for use in the `occ_search()` function. `scientificName`: search for a species by name (instead of `taxonKey`). `continent`: search by continent. `lastInterpreted`: search by last time GBIF modified the record. `recordNumber`: search by the data collector's specimen record number - this is different from the GBIF record number. `typeStatus`: search by specimen type status. (#74)
* Note that given the new parameters many more options are available for implicit faceted search in which you can pass many values in a vector to do multiple searches like `parameterName = c(x, y, z)`. These parameters are: `taxonKey`, `scientificName`, `datasetKey`, `catalogNumber`, `collectorName`, `geometry`, `country`, `recordNumber`, `search`, `institutionCode`, `collectionCode`, `decimalLatitude`, `decimalLongitude`, `depth`, `year`, `typeStatus`, `lastInterpreted`, and `continent`. This isn't faceted search server side - this is just looping your different values of the parameter against the GBIF API.
* Range queries are a new feature in the GBIF API. Some parameters in `occ_search()` now support range queries: `decimalLatitude`,`decimalLongitude`,`depth`,`elevation`,`eventDate`,`lastInterpreted`,`month`, and `year`. Do a range query for example by `depth=50,100` to ask for occurrences where depth was recorded between 50 and 100 meters. Note that this syntax `depth=c(50,100)` will perform two separate searches, one for `depth=50` and one for `depth=100`. (#71)

rgbif 0.5.0
===============

### IMPROVEMENTS

* Changed name of country_codes() function to gbif_country_codes() to avoid conflicts with other packages.
* Replaced sapply() with vapply() throughout the package as it is more robust and can be faster.
* Added a startup message to the package.
* gbifmap() now plots a map with ggplot2::coord_fixed(ratio=1) so that you don't get wonky maps.
* occ_count() now accepts a call to query publishingCountry with a single parameter (country), to list occurrence counts by publishing country.
* occ_get() and occ_search() lose parameter minimal, and in its place gains parameter fields, in which you can request fields='minimal' to get just name, taxon key, lat and long. Or set to 'all' to get all fields, or selection the fields you want by passing in a vector of field names.

### BUG FIXES

* Updated base url for the GIBF parser function parsenames()
* isocodes dataset now with documentation.

### NEW FEATURES

* New function count_facet() to do facetted count search, as GBIF doesn't allow faceted searches against the count API.
* New function elevation() to get elevation data for a data.frame of lat/long points, or a list of lat/long points. This function uses the Google Elevation API (https://developers.google.com/maps/documentation/elevation/).
* New function installations() to get metadata on installations.

rgbif 0.4.1
===============

### BUG FIXES

* Improved handling of limit parameter in occ_search() so that the correct number of occurrences are returned.
* Fixed various tests that were broken.

### IMPROVEMENTS

* Added missing limit argument in datasets() function man file, also function gains start and callopts parameters.

rgbif 0.4.0
===============

### IMPROVEMENTS

* Data object isocodes gains new column gbif_names, the GBIF specific names for countries.
* Added in deprecation messages throughtout package for functions and arguments that are deprecated.
* tests moved to tests/testthat from inst/tests.
* Vignettes now in vignettes/ directory.

### NEW FEATURES

* New function dataset_suggest(), a quick autocomplete service that returns up to 20 datasets.
* New function name_backbone() looks up names against the GBIF backbone taxonomy.
* New function name_suggest(), a quick autocomplete service that returns up to 20 name usages.
* New function occ_metadata() to search dataset metadata.
* New function parsenames() that parses taxonomic names and returns their components.

rgbif 0.3.9
===============

### IMPROVEMENTS

* Added back in functions, and .Rd files, from old version or rgbif that interacts with the old GBIF API.
* Updated vignette to work with new GBIF API and fxns.

### NEW FEATURES

* Added functions to interact with the new GBIF API, notably: country_codes(), dataset_metrics(), dataset_search(), datasets(), name_lookup(), gbifmap(), gist(), name_lookup(), name_usage(), networks(), nodes(), occ_count(), occ_get(), occ_search(), organizations(), stylegeojson(), togeojson(). See the README for a crosswalk from old functions to new ones.

### BUG FIXES

* test files moved from inst/tests/ to tests/testthat/

rgbif 0.3.2
===============

### BUG FIXES

* Removed georeferencedonly parameter - is deprecated in the GBIF API


rgbif 0.3.0
===============

### IMPROVEMENTS

* Added S3 objects: Output from calls to occurrencelist() and occurrence list_many() now of class gbiflist, and output from calls to densitylist() now of class gbifdens.
* Slight changes to gbifmaps() function.
* url parameter in all functions moved into the function itself as the base GBIF API url doesn't need to be specified by user.
* Vignette added.

### NEW FEATURES

* Added function country_codes() to look up 2 character ISO country codes for use in searches.
* Added function occurrencelist_many() to handle searches of many species.
* Added functions togeojson() and stylegeosjon() to convert a data.frame with lat/long columns to geojson file format, and to add styling to data.frames before using togeojson() .
* occurrencelist() and occurrencelist_many() gain argument fixnames, which lets user change species names in output data.frame according to a variety of scenarios.
* taxonsearch() gains argument accepted_status to accept only those names that have a status of accepted. In addition, this function has significant changes, and examples, to improve performance.


rgbif 0.2.0
===============

### IMPROVEMENTS

* Improved code style, and simplified code in some functions.

### NEW FEATURES

* occurrencelist() now handles scientific notation when maxresults are given in that form.
* occurencelist() now can retrieve any number of records; was previously a max of 1000 records.

### BUG FIXES

* Demo "List" was returning incorrect taxon names - corrected now.
* Removed unused parameter 'latlongdf' in occurencelist().


rgbif 0.1.5
===============

### IMPROVEMENTS

* Changed all functions to use RCurl instead of httr as httr was presenting some problems.
* Two function, capwords and gbifxmlToDataFrame, added with documentation as internal functions.

### NEW FEATURES

* Added function density_spplist to get a species list or data.frame of species and their counts for any degree cell.
* Added function densitylist to access to records showing the density of occurrence records from the GBIF Network by one-degree cell.
* Added function gbifmap to make a simple map to visualize GBIF data.
* Added function occurrencecount to count taxon concept records matching a range of filters.

DEPRECATED

* gbifdatause removed, was just a function to return the data sharing agreement from GBIF.


rgbif 0.1.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.1.0
* ubuntu 16.04 (on GitHub Actions), R 4.1.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 15 reverse dependencies. Reverse dependency checks at <https://github.com/ropensci/rgbif/actions?query=workflow%3Arevdep>. No problems were found related to this package.

--------

This version fixes many bugs and adds some additional search capabilities.

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
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/rgbif/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rgbif.git`
* Make sure to track progress upstream (i.e., on our version of `rgbif` at `ropensci/rgbif`) by doing `git remote add upstream https://github.com/ropensci/rgbif.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rgbif`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Email?

If you do email you'll get a slower response than if you open an issue at https://github.com/ropensci/rgbif/issues. It benefits everyone to have an open discussion in an issue.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<!-- IF YOU DO NOT INCLUDE YOUR SESSION INFO YOUR ISSUE WILL LIKELY BE CLOSED WITHOUT FURTHER CONSIDERATION -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-06-30 r78761) |
|os       |macOS Catalina 10.15.5                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-07-22                                  |

# Dependencies

|package |old   |new      |Δ  |
|:-------|:-----|:--------|:--|
|rgbif   |3.1.0 |3.1.1.91 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*<!-- README.md is generated from README.Rmd. Please edit that file and knit -->

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.path = "inst/assets/img/"
)
```

# rgbif <img src="man/figures/logo.png" align="right" alt="" width="120">

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rgbif)](https://cranchecks.info/pkgs/rgbif)
[![R-CMD-check](https://github.com/ropensci/rgbif/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rgbif/actions?query=workflow%3AR-CMD-check)
[![real-requests](https://github.com/ropensci/rgbif/workflows/R-check-real-requests/badge.svg)](https://github.com/ropensci/rgbif/actions?query=workflow%3AR-check-real-requests)
[![codecov.io](https://codecov.io/github/ropensci/rgbif/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rgbif?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rgbif)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rgbif)](https://cran.r-project.org/package=rgbif)
[![DOI](https://zenodo.org/badge/2273724.svg)](https://zenodo.org/badge/latestdoi/2273724)

`rgbif` gives you access to data from [GBIF][] via their REST API. GBIF versions their API - we are currently using `v1` of their API. You can no longer use their old API in this package - see `?rgbif-defunct`.

Please cite rgbif. Run the following to get the appropriate citation for the version you're using:

```r
citation(package = "rgbif")
```

To get started, see:

* rgbif vignette (https://docs.ropensci.org/rgbif/articles/rgbif.html): an introduction to the package's main functionalities.
* Function reference (https://docs.ropensci.org/rgbif/reference/index.html): an overview of all `rgbif` functions.
* Articles (https://docs.ropensci.org/rgbif/articles/index.html): vignettes/tutorials on how to download data, clean data, and work with taxonomic names.

Check out the `rgbif` [paper][] for more information on this package and the sister [Python][pygbif], [Ruby][gbifrb], and [PHP][phpgbif] clients.

Note: Maximum number of records you can get with `occ_search()` and `occ_data()` is 100,000. See https://www.gbif.org/developer/occurrence

## Installation

```{r eval=FALSE}
install.packages("rgbif")
```

Or, install development version

```{r eval=FALSE}
pak::pkg_install("ropensci/rgbif")
# OR
install.packages("rgbif", repos="https://dev.ropensci.org")
```

```{r}
library("rgbif")
```

## Screencast

<a href="https://vimeo.com/127119010"><img src="man/figures/README-screencast.png" width="400"></a>

## Contributors

This list honors all contributors in alphabetical order. Code contributors are in bold.

[adamdsmith](https://github.com/adamdsmith) - [AgustinCamacho](https://github.com/AgustinCamacho) - [AldoCompagnoni](https://github.com/AldoCompagnoni) - [AlexPeap](https://github.com/AlexPeap) - [andzandz11](https://github.com/andzandz11) - [AshleyWoods](https://github.com/AshleyWoods) - [AugustT](https://github.com/AugustT) - [barthoekstra](https://github.com/barthoekstra) - **[benmarwick](https://github.com/benmarwick)** - [cathynewman](https://github.com/cathynewman) - [cboettig](https://github.com/cboettig) - [coyotree](https://github.com/coyotree) - **[damianooldoni](https://github.com/damianooldoni)** - [dandaman](https://github.com/dandaman) - [djokester](https://github.com/djokester) - [dlebauer](https://github.com/dlebauer) - **[dmcglinn](https://github.com/dmcglinn)** - [dnoesgaard](https://github.com/dnoesgaard) - [DupontCai](https://github.com/DupontCai) - [ecology-data-science](https://github.com/ecology-data-science) - [EDiLD](https://github.com/EDiLD) - [elgabbas](https://github.com/elgabbas) - [emhart](https://github.com/emhart) - [fxi](https://github.com/fxi) - [ghost](https://github.com/ghost) - [gkburada](https://github.com/gkburada) - [hadley](https://github.com/hadley) - [Huasheng12306](https://github.com/Huasheng12306) - [ibartomeus](https://github.com/ibartomeus) - **[JanLauGe](https://github.com/JanLauGe)** - **[jarioksa](https://github.com/jarioksa)** - **[jeroen](https://github.com/jeroen)** - [jhnwllr](https://github.com/jhnwllr) - [jhpoelen](https://github.com/jhpoelen) - [jivelasquezt](https://github.com/jivelasquezt) - [jkmccarthy](https://github.com/jkmccarthy) - **[johnbaums](https://github.com/johnbaums)** - [jtgiermakowski](https://github.com/jtgiermakowski) - [jwhalennds](https://github.com/jwhalennds) - **[karthik](https://github.com/karthik)** - [kgturner](https://github.com/kgturner) - [Kim1801](https://github.com/Kim1801) - [ljuliusson](https://github.com/ljuliusson) - [ljvillanueva](https://github.com/ljvillanueva) - [luisDVA](https://github.com/luisDVA) - [martinpfannkuchen](https://github.com/martinpfannkuchen) - [MattBlissett](https://github.com/MattBlissett) - [MattOates](https://github.com/MattOates) - [maxhenschell](https://github.com/maxhenschell) - **[mdsumner](https://github.com/mdsumner)** - [no-la-ngo](https://github.com/no-la-ngo) - [Octoberweather](https://github.com/Octoberweather) - [Pakillo](https://github.com/Pakillo) - **[peterdesmet](https://github.com/peterdesmet)** - [PhillRob](https://github.com/PhillRob) - [poldham](https://github.com/poldham) - [qgroom](https://github.com/qgroom) - [raymondben](https://github.com/raymondben) - [rossmounce](https://github.com/rossmounce) - [sacrevert](https://github.com/sacrevert) - [sagitaninta](https://github.com/sagitaninta) - **[sckott](https://github.com/sckott)** - [scottsfarley93](https://github.com/scottsfarley93) - [simon-tarr](https://github.com/simon-tarr) - **[SriramRamesh](https://github.com/SriramRamesh)** - [stevenpbachman](https://github.com/stevenpbachman) - [stevensotelo](https://github.com/stevensotelo) - **[stevenysw](https://github.com/stevenysw)** - [TomaszSuchan](https://github.com/TomaszSuchan) - [tphilippi](https://github.com/tphilippi) - [vandit15](https://github.com/vandit15) - [vervis](https://github.com/vervis) - **[vijaybarve](https://github.com/vijaybarve)** - [willgearty](https://github.com/willgearty) - [Xuletajr](https://github.com/Xuletajr) - [yvanlebras](https://github.com/yvanlebras) - [zixuan75](https://github.com/zixuan75)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rgbif/issues).
* License: MIT
* Get citation information for `rgbif` in R doing `citation(package = 'rgbif')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

- - -

This package is part of [spocc](https://github.com/ropensci/spocc), along with several other packages, that provide access to occurrence records from multiple data sources.

- - -

[mapr]: https://github.com/ropensci/mapr
[paper]: https://peerj.com/preprints/3304/
[GBIF]: https://www.gbif.org/
[pygbif]: https://github.com/sckott/pygbif
[gbifrb]: https://github.com/sckott/gbifrb
[phpgbif]: https://gitlab.res-telae.cat/restelae/php-gbif
---
title: Downloading data from GBIF
author: Scott Chamberlain
date: "2021-06-11"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{downloading}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---




GBIF provides two ways to get occurrence data: through the
`/occurrence/search` route (see `occ_search()` and `occ_data()`),
or via the `/occurrence/download` route (many functions, see below).
`occ_search()`/`occ_data()` are more appropriate for smaller data, while
`occ_download*()` functions are more appropriate for larger data requests.
Note that the download service is equivalent to downloading a dataset
from the GBIF website - but doing it here makes it reproducible
(and easier once you learn the ropes)!

The download functions are:

- `occ_download()` - Start a download
- `occ_download_prep()` - Prepare a download request
- `occ_download_queue()` - Start many downloads in a queue
- `occ_download_meta()` - Get metadata progress on a single download
- `occ_download_list()` - List your downloads
- `occ_download_cancel()` - Cancel a download
- `occ_download_cancel_staged()` - Cancels any jobs with status `RUNNING`
or `PREPARING`
- `occ_download_get()` - Retrieve a download
- `occ_download_import()` - Import a download from local file system
- `occ_download_datasets()` - List datasets for a download
- `occ_download_dataset_activity()` - Lists the downloads activity of a dataset

The following functions are the download predicate DSL, used with `occ_download`,
`occ_download_prep()`, and `occ_download_queue()`.

`pred*` functions are named for the 'type' of operation they do, following
the terminology used by GBIF, see
https://www.gbif.org/developer/occurrence#predicates

Function names are given, with the equivalent GBIF type value (e.g.,
`pred_gt` and `greaterThan`)

The following functions take one key and one value:

- `pred`: equals
- `pred_lt`: lessThan
- `pred_lte`: lessThanOrEquals
- `pred_gt`: greaterThan
- `pred_gte`: greaterThanOrEquals
- `pred_not`: not
- `pred_like`: like

The following function is only for geospatial queries, and only
accepts a WKT string:

- `pred_within`: within

The following function is only for stating the you don't want
a key to be null, so only accepts one key:

- `pred_notnull`: isNotNull

The following two functions accept multiple individual predicates,
separating them by either "and" or "or":

- `pred_and`: and
- `pred_or`: or

The following function is special in that it accepts a single key
but many values; stating that you want to search for all the values:

- `pred_in`: in

`occ_download()` is the function to start off with when using the GBIF download
service. With it you can specify what query you want. Unfortunately, the interfaces
to the search vs. download services are different, so we couldn't make the `rgbif`
interface to `occ_search`/`occ_data` the same as `occ_download`.

Be aware that you can only perform 3 downloads simultaneously, so plan wisely.
To help with this limitation, we are working on a queue helper, but it's not
ready yet.

Let's take a look at how to use the download functions.

## Load rgbif


```r
library("rgbif")
```

## Kick off a download

Instead of passing parameters like `taxonkey = 12345` in `occ_search`, for downloads
we construct queries using any of the `pred` functions. This provides for much more
complex queries than can be done with `occ_search/occ_data`. With downloads, you can use operators other than `=` (equal to); and the queries can be very complex.


```r
(res <- occ_download(pred('taxonKey', 7264332), pred('hasCoordinate', TRUE)))
#> <<gbif download>>
#>   Username: sckott
#>   E-mail: myrmecocystus@gmail.com
#>   Download key: 0000796-171109162308116
```

What `occ_download` returns is not the data itself!  When you send the request to
GBIF, they have to prepare it first, then when it's done you can download it.

What `occ_download` returns is some useful metadata that tells you about the
download, and helps us check and know when the download is done.


## Check download status

After running `occ_download`, we can pass the resulting object to
`occ_download_meta` - with primary goal of checking the download status.


```r
occ_download_meta(res)
#> <<gbif download metadata>>
#>   Status: PREPARING
#>   Format: DWCA
#>   Download key: 0000796-171109162308116
#>   Created: 2017-11-10T19:30:32.328+0000
#>   Modified: 2017-11-10T19:30:44.590+0000
#>   Download link: http://api.gbif.org/v1/occurrence/download/request/0000796-171109162308116.zip
#>   Total records: 1425
#>   Request:
#>     type:  and
#>     predicates:
#>       > type: equals, key: TAXON_KEY, value: 7264332
#>       > type: equals, key: HAS_COORDINATE, value: TRUE
```

Continue running `occ_download_meta` until the `Status` value is `SUCCEEDED`
or `KILLED`. If it is `KILLED` that means something went wrong - get in touch
with us. If `SUCCEEDED`, then you can proceed to the next step (downloading
the data with `occ_download_get`).

Before we go to the next step, there's another function to help you out.

With `occ_download_list` you can get an overview of all your download
requests, with


```r
x <- occ_download_list()
x$results <- tibble::as_tibble(x$results)
x
#> $meta
#>   offset limit endofrecords count
#> 1      0    20        FALSE   211
#>
#> $results
#> # A tibble: 20 x 18
#>                        key                    doi
#>  *                   <chr>                  <chr>
#>  1 0000796-171109162308116 doi:10.15468/dl.nv3r5p
#>  2 0000739-171109162308116 doi:10.15468/dl.jmachn
#>  3 0000198-171109162308116 doi:10.15468/dl.t5wjpe
#>  4 0000122-171020152545675 doi:10.15468/dl.yghxj7
#>  5 0000119-171020152545675 doi:10.15468/dl.qiowtc
#>  6 0000115-171020152545675 doi:10.15468/dl.tdbkzn
#>  7 0010067-170714134226665 doi:10.15468/dl.ro6qj1
#>  8 0010066-170714134226665 doi:10.15468/dl.bhekhi
#>  9 0010065-170714134226665 doi:10.15468/dl.xy4nfp
#> 10 0010064-170714134226665 doi:10.15468/dl.hsqp84
#> 11 0010062-170714134226665 doi:10.15468/dl.h2apik
#> 12 0010061-170714134226665 doi:10.15468/dl.1srstq
#> 13 0010059-170714134226665 doi:10.15468/dl.2me5hk
#> 14 0010058-170714134226665 doi:10.15468/dl.sjmxvf
#> 15 0010057-170714134226665 doi:10.15468/dl.f28182
#> 16 0010056-170714134226665 doi:10.15468/dl.4t2qim
#> 17 0010055-170714134226665 doi:10.15468/dl.lumz7s
#> 18 0010054-170714134226665 doi:10.15468/dl.wfkgqm
#> 19 0010053-170714134226665 doi:10.15468/dl.fintow
#> 20 0010050-170714134226665 doi:10.15468/dl.a2h9gu
#> # ... with 16 more variables: license <chr>, created <chr>, modified <chr>,
#> #   status <chr>, downloadLink <chr>, size <dbl>, totalRecords <int>,
#> #   numberDatasets <int>, request.creator <chr>, request.format <chr>,
#> #   request.notificationAddresses <list>, request.sendNotification <lgl>,
#> #   request.predicate.type <chr>, request.predicate.predicates <list>,
#> #   request.predicate.key <chr>, request.predicate.value <chr>
```

## Canceling downloads

If for some reason you need to cancel a download you can do so with
`occ_download_cancel` or `occ_download_cancel_staged`.

`occ_download_cancel` cancels a job by download key, while `occ_download_cancel_staged`
cancels all jobs in `PREPARING` or `RUNNING` stage.

## Fetch data

After you see the `SUCCEEDED` status on calling `occ_download_meta`, you can
then download the data using `occ_download_get`.


```r
(dat <- occ_download_get("0000796-171109162308116"))
#> <<gbif downloaded get>>
#>   Path: ./0000796-171109162308116.zip
#>   File size: 0.35 MB
```

This only download data to your machine - it does not read it into R.
You can now move on to importing into R.

## Import data into R

Pass the output of `occ_download_get` directly to `occ_download_import` -
they can be piped together if you like.


```r
occ_download_get("0000796-171109162308116") %>% occ_download_import()
# OR
dat <- occ_download_get("0000796-171109162308116")
occ_download_import(dat)
#> # A tibble: 1,425 x 235
#>        gbifID abstract accessRights accrualMethod accrualPeriodicity
#>         <int>    <lgl>        <chr>         <lgl>              <lgl>
#>  1 1667184715       NA                         NA                 NA
#>  2 1667182218       NA                         NA                 NA
#>  3 1667179996       NA                         NA                 NA
#>  4 1667179527       NA                         NA                 NA
#>  5 1667171607       NA                         NA                 NA
#>  6 1667165448       NA                         NA                 NA
#>  7 1667163154       NA                         NA                 NA
#>  8 1667162324       NA                         NA                 NA
#>  9 1667162162       NA                         NA                 NA
#> 10 1667161552       NA                         NA                 NA
#> # ... with 1,415 more rows, and 230 more variables: accrualPolicy <lgl>,
#> #   alternative <lgl>, audience <lgl>, available <lgl>,
#> #   bibliographicCitation <lgl>, conformsTo <lgl>, contributor <lgl>,
#> #   coverage <lgl>, created <lgl>, creator <lgl>, date <lgl>,
#> #   dateAccepted <lgl>, dateCopyrighted <lgl>, dateSubmitted <lgl>,
#> #   description <lgl>, educationLevel <lgl>, extent <lgl>, format <lgl>,
#> #   hasFormat <lgl>, hasPart <lgl>, hasVersion <lgl>, identifier <chr>,
#>
#> ... cut for brevity
```

## Citing download data

The nice thing about data retrieved via GBIF's download service is that they
provide DOIs for each download, so that you can give a link that resolves to
the download with metadata on GBIF's website. And it makes for a nice citation.

Using the funciton `gbif_citaiton` we can get citations for our downloads,
with the output from `occ_download_get` or `occ_download_meta`.


```r
occ_download_meta(res) %>% gbif_citation()
#> $download
#> [1] "GBIF Occurrence Download https://doi.org/10.15468/dl.ohjevv Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2017-11-10"
#>
#> $datasets
#> NULL
```

You'll notice that the `datasets` slot is `NULL` - because when using `occ_download_meta`,
we don't yet have any information about which datasets are in the download.

But if you use `occ_download_get` you then have the individual datasets, and we
can get citatations for each idividual dataset in addition to the entire download.


```r
occ_download_get("0000796-171109162308116") %>% gbif_citation()
#> $download
#> [1] "GBIF Occurrence Download https://doi.org/10.15468/dl.nv3r5p Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2017-11-10"
#>
#> $datasets
#> $datasets[[1]]
#> <<rgbif citation>>
#>    Citation: Office of Environment & Heritage (2017). OEH Atlas of NSW
#>         Wildlife. Occurrence Dataset https://doi.org/10.15468/14jd9g accessed
#>         via GBIF.org on 2017-11-10.. Accessed from R via rgbif
#>         (https://github.com/ropensci/rgbif) on 2017-11-10
#>    Rights: This work is licensed under a Creative Commons Attribution (CC-BY)
#>         4.0 License.
#>
#> $datasets[[2]]
#> <<rgbif citation>>
#>    Citation: Creuwels J (2017). Naturalis Biodiversity Center (NL) - Botany.
#>         Naturalis Biodiversity Center. Occurrence Dataset
#>         https://doi.org/10.15468/ib5ypt accessed via GBIF.org on 2017-11-10..
#>         Accessed from R via rgbif (https://github.com/ropensci/rgbif) on
#>         2017-11-10
#>    Rights: To the extent possible under law, the publisher has waived all
#>         rights to these data and has dedicated them to the Public Domain (CC0
#>         1.0). Users may copy, modify, distribute and use the work, including
#>         for commercial purposes, without restriction.
#>
#> ... cutoff for brevity
```

Here, we get the overall citation as well as citations (and data rights) for each dataset.

Please do cite the data you use from GBIF!
---
title: Taxonomic names
author: Scott Chamberlain
date: "2021-06-11"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{taxonomic_names}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



You have probably, or will, run into problems with taxonomic names. For example,
you may think you know how a taxonomic name is spelled, but then GBIF will not
agree with you. Or, perhaps GBIF will have multiple versions of the taxon,
spelled in slightly different ways. Or, the version of the name that they think
is the _right one_ does not match what yo think is the right one.

This isn't really anyone's fault. It's a result of there not being one accepted
taxonomic source of truth across the globe. There are many different taxonomic
databases. GBIF makes their own _backbone taxonomy_ that they use as a source
of internal truth for taxonomic names. The accepted names in the backbone taxonomy
match those in the database of occurrences - so do try to figure out what
backbone taxonomy version of the name you want.

Another source of problems stems from the fact that names are constantly changing.
Sometimes epithets change, sometimes generic names, and sometimes higher names
like family or tribe. These changes can take a while to work their way in to
GBIF's data.

The following are some examples of confusing name bits. We'll update these if
GBIF's name's change. The difference between each pair of names is highlighted
in bold.

## Load rgbif


```r
library("rgbif")
```

## Helper function

To reduce code duplication, we'll use a little helper function to make a call
to `name_backbone()` for each input name, then `rbind` them together:


```r
name_rbind <- function(..., rank = "species") {
  columns <- c('usageKey', 'scientificName', 'canonicalName', 'rank',
    'status', 'confidence', 'matchType', 'synonym')
  df <- lapply(list(...), function(w) {
    rgbif::name_backbone(w, rank = rank)[, columns]
  })
  data.frame(do.call(rbind, df))
}
```

And another function to get the taxonomic data provider


```r
taxon_provider <- function(x) {
  tt <- name_usage(key = x)$data
  datasets(uuid = tt$constituentKey)$data$title
}
```

We use `taxon_provider()` below to get the taxonomy provider in the bulleted list of details
for each taxon (even though you don't see it called, we use it, but the code isn't shown :)).

## Pinus sylvestris vs. P. silvestris


```r
(c1 <- name_rbind("Pinus sylvestris", "Pinus silvestris"))
#>   usageKey      scientificName    canonicalName    rank   status confidence
#> 1  5285637 Pinus sylvestris L. Pinus sylvestris SPECIES ACCEPTED         97
#> 2  5285637 Pinus sylvestris L. Pinus sylvestris SPECIES ACCEPTED         94
#>   matchType synonym
#> 1     EXACT   FALSE
#> 2     FUZZY   FALSE
```

* P. s<b>y</b>lvestris w/ occurrences count from the data provider


```r
occ_count(c1$usageKey[[1]])
#> [1] 642901
taxon_provider(c1$usageKey[[1]])
#> [1] "Catalogue of Life - May 2021"
```

* P. s<b>i</b>lvestris w/ occurrences count from the data provider


```r
occ_count(c1$usageKey[[2]])
#> [1] 642901
taxon_provider(c1$usageKey[[2]])
#> [1] "Catalogue of Life - May 2021"
```

## Macrozamia platyrachis vs. M. platyrhachis


```r
(c2 <- name_rbind("Macrozamia platyrachis", "Macrozamia platyrhachis"))
#>   usageKey                     scientificName           canonicalName    rank
#> 1  2683551 Macrozamia platyrhachis F.M.Bailey Macrozamia platyrhachis SPECIES
#> 2  2683551 Macrozamia platyrhachis F.M.Bailey Macrozamia platyrhachis SPECIES
#>     status confidence matchType synonym
#> 1 ACCEPTED         96     FUZZY   FALSE
#> 2 ACCEPTED         98     EXACT   FALSE
```

* M. platyrachis w/ occurrences count from the data provider


```r
occ_count(c2$usageKey[[1]])
#> [1] 61
taxon_provider(c2$usageKey[[1]])
#> [1] "Catalogue of Life - May 2021"
```

* M. platyr<b>h</b>achis w/ occurrences count from the data provider


```r
occ_count(c2$usageKey[[2]])
#> [1] 61
taxon_provider(c2$usageKey[[2]])
#> [1] "Catalogue of Life - May 2021"
```

## Cycas circinalis vs. C. circinnalis


```r
(c3 <- name_rbind("Cycas circinalis", "Cycas circinnalis"))
#>   usageKey      scientificName    canonicalName    rank   status confidence
#> 1  2683264 Cycas circinalis L. Cycas circinalis SPECIES ACCEPTED         98
#> 2  2683264 Cycas circinalis L. Cycas circinalis SPECIES ACCEPTED         95
#>   matchType synonym
#> 1     EXACT   FALSE
#> 2     FUZZY   FALSE
```

* C. circinalis w/ occurrences count from the data provider


```r
occ_count(c3$usageKey[[1]])
#> [1] 724
taxon_provider(c3$usageKey[[1]])
#> [1] "Catalogue of Life - May 2021"
```

* C. circin<b>n</b>alis w/ occurrences count from the data provider


```r
occ_count(c3$usageKey[[2]])
#> [1] 724
taxon_provider(c3$usageKey[[2]])
#> [1] "Catalogue of Life - May 2021"
```

## Isolona perrieri vs. I. perrierii


```r
(c4 <- name_rbind("Isolona perrieri", "Isolona perrierii"))
#>   usageKey          scientificName     canonicalName    rank   status
#> 1  6308376 Isolona perrierii Diels Isolona perrierii SPECIES ACCEPTED
#> 2  6308376 Isolona perrierii Diels Isolona perrierii SPECIES ACCEPTED
#>   confidence matchType synonym
#> 1         96     FUZZY   FALSE
#> 2         98     EXACT   FALSE
```

* I. perrieri w/ occurrences count from the data provider


```r
occ_count(c4$usageKey[[1]])
#> [1] 92
taxon_provider(c4$usageKey[[1]])
#> [1] "Catalogue of Life - May 2021"
```

* I. perrieri<b>i</b> w/ occurrences count from the data provider


```r
occ_count(c4$usageKey[[2]])
#> [1] 92
taxon_provider(c4$usageKey[[2]])
#> [1] "Catalogue of Life - May 2021"
```

## Wiesneria vs. Wisneria


```r
(c5 <- name_rbind("Wiesneria", "Wisneria", rank = "genus"))
#>   usageKey         scientificName canonicalName  rank   status confidence
#> 1  2864604      Wiesneria Micheli     Wiesneria GENUS ACCEPTED         96
#> 2  7327444 Wisneria Micheli, 1881      Wisneria GENUS  SYNONYM         95
#>   matchType synonym
#> 1     EXACT   FALSE
#> 2     EXACT    TRUE
```

* Wi<b>e</b>sneria w/ occurrences count from the data provider


```r
occ_count(c5$usageKey[[1]])
#> [1] 120
taxon_provider(c5$usageKey[[1]])
#> [1] "Catalogue of Life - May 2021"
```

* Wisneria w/ occurrences count from the data provider


```r
occ_count(c5$usageKey[[2]])
#> [1] 3
taxon_provider(c5$usageKey[[2]])
#> [1] "The Interim Register of Marine and Nonmarine Genera"
```

## The take away messages from this vignette

* Make sure you are using the name you think you're using
* Realize that GBIF's backbone taxonomy is used for occurrence data
* Searching for occurrences by name matches against backbone names, 
not other names (e.g., synonyms)
* GBIF may at some points in time have multiple version of the same name in their own backbone taxonomy - These can usually be separated by data provider (e.g., Catalogue of Life vs. International Plant Names Index)
* There are different ways to search for names - make sure are familiar 
with the four different name search functions, all starting with 
`name_`
---
title: Introduction to rgbif
author: Scott Chamberlain
date: "2021-06-11"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



Seach and retrieve data from the Global Biodiverity Information Facilty (GBIF)

## About the package

`rgbif` is an R package to search and retrieve data from the Global Biodiverity Information Facilty (GBIF). `rgbif` wraps R code around the [GBIF API][gbifapi] to allow you to talk to GBIF from R.


## Get rgbif

Install from CRAN


```r
install.packages("rgbif")
```

Or install the development version from GitHub


```r
remotes::install_github("ropensci/rgbif")
```

Load rgbif


```r
library("rgbif")
```

## Number of occurrences

Search by type of record, all observational in this case


```r
occ_count(basisOfRecord='OBSERVATION')
#> [1] 18483075
```

Records for **Puma concolor** with lat/long data (georeferened) only. Note that `hasCoordinate` in `occ_search()` is the same as `georeferenced` in `occ_count()`.


```r
occ_count(taxonKey=2435099, georeferenced=TRUE)
#> [1] 8093
```

All georeferenced records in GBIF


```r
occ_count(georeferenced=TRUE)
#> [1] 1605270476
```

Records from Denmark


```r
denmark_code <- isocodes[grep("Denmark", isocodes$name), "code"]
occ_count(country=denmark_code)
#> [1] 44757342
```

Number of records in a particular dataset


```r
occ_count(datasetKey='9e7ea106-0bf8-4087-bb61-dfe4f29e0f17')
#> [1] 4591
```

All records from 2012


```r
occ_count(year=2012)
#> [1] 63793855
```

Records for a particular dataset, and only for preserved specimens


```r
occ_count(datasetKey='e707e6da-e143-445d-b41d-529c4a777e8b', basisOfRecord='OBSERVATION')
#> [1] 0
```

## Search for taxon names

Get possible values to be used in taxonomic rank arguments in functions


```r
taxrank()
#> [1] "kingdom"       "phylum"        "class"         "order"        
#> [5] "family"        "genus"         "species"       "subspecies"   
#> [9] "infraspecific"
```

`name_lookup()` does full text search of name usages covering the scientific and vernacular name, the species description, distribution and the entire classification across all name usages of all or some checklists. Results are ordered by relevance as this search usually returns a lot of results.

By default `name_lookup()` returns five slots of information: meta, data, facets, hierarchies, and names. hierarchies and names elements are named by their matching GBIF key in the `data.frame` in the data slot.


```r
out <- name_lookup(query='mammalia')
```


```r
names(out)
#> [1] "meta"        "data"        "facets"      "hierarchies" "names"
```


```r
out$meta
#> # A tibble: 1 x 4
#>   offset limit endOfRecords count
#>    <int> <int> <lgl>        <int>
#> 1      0   100 FALSE         2341
```


```r
head(out$data)
#> # A tibble: 6 x 26
#>       key scientificName datasetKey       nubKey parentKey parent kingdom phylum
#>     <int> <chr>          <chr>             <int>     <int> <chr>  <chr>   <chr> 
#> 1  1.60e8 Mammalia       677a9818-ca48-4…    359 159534408 Chord… Animal… Chord…
#> 2  1.77e8 Mammalia       a364ceb5-1864-4…    359 176575081 Chord… Animal… Chord…
#> 3  1.47e8 Mammalia       a4e57579-9638-4…    359 147439639 Chord… Animal… Chord…
#> 4  1.47e8 Mammalia       be44f76b-73cf-4…    359 147440318 Chord… Animal… Chord…
#> 5  1.64e8 Mammalia       c383b17f-e3bc-4…    359 164285181 Chord… Animal… Chord…
#> 6  1.64e8 Mammalia       00e791be-36ae-4…    359 164302487 Chord… Animal… Chord…
#> # … with 18 more variables: kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   canonicalName <chr>, authorship <chr>, nameType <chr>,
#> #   taxonomicStatus <chr>, rank <chr>, origin <chr>, numDescendants <int>,
#> #   numOccurrences <int>, habitats <chr>, nomenclaturalStatus <lgl>,
#> #   threatStatuses <chr>, synonym <lgl>, class <chr>, constituentKey <chr>,
#> #   extinct <lgl>
```


```r
out$facets
#> NULL
```


```r
out$hierarchies[1:2]
#> $`159534411`
#>     rankkey     name
#> 1 159534378 Animalia
#> 2 159534408 Chordata
#> 
#> $`176575082`
#>     rankkey     name
#> 1 176575080 Animalia
#> 2 176575081 Chordata
```


```r
out$names[2]
#> $<NA>
#> NULL
```

Search for a genus


```r
z <- name_lookup(query='Cnaemidophorus', rank="genus")
z$data
#> # A tibble: 28 x 37
#>        key scientificName datasetKey      nubKey parentKey parent kingdom phylum
#>      <int> <chr>          <chr>            <int>     <int> <chr>  <chr>   <chr> 
#>  1  1.59e8 Cnaemidophorus 23905003-5ee5… 1858636 159439401 Ptero… Animal… Arthr…
#>  2  1.58e8 Cnaemidophorus 4cec8fef-f129… 1858636 157904443 Ptero… Animal… Arthr…
#>  3  1.69e8 Cnaemidophorus 4b3e4a71-704a… 1858636 168525701 Ptero… Animal… Arthr…
#>  4  1.83e8 Cnaemidophorus 848271aa-6a4f… 1858636 182646345 Lepid… Animal… <NA>  
#>  5  1.24e8 Cnaemidophorus fab88965-e69d…      NA 104446806 Ptero… Metazoa Arthr…
#>  6  1.78e8 Cnaemidophorus 6b6b2923-0a10… 1858636 177881660 Ptero… Metazoa Arthr…
#>  7  1.82e8 Cnaemidophorus 16c3f9cb-4b19… 1858636 100557623 Ptero… <NA>    <NA>  
#>  8  1.82e8 Cnaemidophorus cbb6498e-8927… 1858636 182338678 Ptero… Animal… Arthr…
#>  9  1.83e8 Cnaemidophorus 4dd32523-a3a3… 1858636 182545750 Ptero… Animal… Arthr…
#> 10  1.80e8 Cnaemidophorus dbaa27eb-29e7… 1858636    161750 Ptero… Animal… <NA>  
#> # … with 18 more rows, and 29 more variables: order <chr>, family <chr>,
#> #   genus <chr>, kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   orderKey <int>, familyKey <int>, genusKey <int>, canonicalName <chr>,
#> #   authorship <chr>, nameType <chr>, taxonomicStatus <chr>, rank <chr>,
#> #   origin <chr>, numDescendants <int>, numOccurrences <int>, habitats <chr>,
#> #   nomenclaturalStatus <chr>, threatStatuses <chr>, synonym <lgl>,
#> #   class <chr>, taxonID <chr>, acceptedKey <int>, accepted <chr>,
#> #   constituentKey <chr>, publishedIn <chr>, extinct <lgl>, accordingTo <chr>
```

Search for the class mammalia


```r
w <- name_lookup(query='mammalia')
w$data
#> # A tibble: 100 x 26
#>        key scientificName datasetKey      nubKey parentKey parent kingdom phylum
#>      <int> <chr>          <chr>            <int>     <int> <chr>  <chr>   <chr> 
#>  1  1.60e8 Mammalia       677a9818-ca48-…    359 159534408 Chord… Animal… Chord…
#>  2  1.77e8 Mammalia       a364ceb5-1864-…    359 176575081 Chord… Animal… Chord…
#>  3  1.47e8 Mammalia       a4e57579-9638-…    359 147439639 Chord… Animal… Chord…
#>  4  1.47e8 Mammalia       be44f76b-73cf-…    359 147440318 Chord… Animal… Chord…
#>  5  1.64e8 Mammalia       c383b17f-e3bc-…    359 164285181 Chord… Animal… Chord…
#>  6  1.64e8 Mammalia       00e791be-36ae-…    359 164302487 Chord… Animal… Chord…
#>  7  1.68e8 Mammalia       6b64ef7e-82f7-…    359 168337497 Chord… Animal… Chord…
#>  8  1.76e8 Mammalia       cd8ba8a5-8504-…    359 176085130 Chord… Animal… Chord…
#>  9  1.77e8 Mammalia       8901e0e4-c1b9-…    359 176574743 Chord… Animal… Chord…
#> 10  1.52e8 Mammalia       961a3a60-22ec-…    359 152056194 Chord… Animal… Chord…
#> # … with 90 more rows, and 18 more variables: kingdomKey <int>,
#> #   phylumKey <int>, classKey <int>, canonicalName <chr>, authorship <chr>,
#> #   nameType <chr>, taxonomicStatus <chr>, rank <chr>, origin <chr>,
#> #   numDescendants <int>, numOccurrences <int>, habitats <chr>,
#> #   nomenclaturalStatus <lgl>, threatStatuses <chr>, synonym <lgl>,
#> #   class <chr>, constituentKey <chr>, extinct <lgl>
```

Look up the species Helianthus annuus


```r
m <- name_lookup(query = 'Helianthus annuus', rank="species")
m$data
#> # A tibble: 100 x 41
#>        key scientificName  datasetKey     nubKey parentKey parent kingdom phylum
#>      <int> <chr>           <chr>           <int>     <int> <chr>  <chr>   <chr> 
#>  1  1.35e8 Helianthus ann… f82a4f7f-6f8… 9206251 180663995 Aster… Plantae Trach…
#>  2  1.15e8 Helianthus ann… ee2aac07-de9… 9206251 144238801 Helia… Plantae Trach…
#>  3  1.46e8 Helianthus ann… 3f5e930b-52a… 9206251 157140516 Helia… Plantae Angio…
#>  4  1.35e8 Helianthus ann… 29d2d5a6-db2… 9206251 181005136 Aster… Plantae Trach…
#>  5  1.63e8 Helianthus ann… 88217638-274… 9206251 163398972 Aster… Plantae Trach…
#>  6  1.03e8 Helianthus ann… fab88965-e69…      NA 103340270 Helia… Viridi… Strep…
#>  7  1.79e8 Helianthus ann… 6b6b2923-0a1… 9206251 178978795 Helia… Viridi… Strep…
#>  8  1.28e8 Helianthus ann… 41c06f1a-23d… 9206251 146770884 Amara… Plantae <NA>  
#>  9  1.35e8 Helianthus ann… 83ca3188-af1… 9206251 180784910 Aster… Plantae Trach…
#> 10  1.35e8 Helianthus ann… 3cabcf37-db1… 9206251 168335901 Aster… Plantae Trach…
#> # … with 90 more rows, and 33 more variables: order <chr>, family <chr>,
#> #   species <chr>, kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   orderKey <int>, familyKey <int>, speciesKey <int>, canonicalName <chr>,
#> #   nameType <chr>, taxonomicStatus <chr>, rank <chr>, origin <chr>,
#> #   numDescendants <int>, numOccurrences <int>, taxonID <chr>, habitats <chr>,
#> #   nomenclaturalStatus <chr>, threatStatuses <chr>, synonym <lgl>,
#> #   class <chr>, genus <chr>, genusKey <int>, authorship <chr>,
#> #   acceptedKey <int>, accepted <chr>, publishedIn <chr>, accordingTo <chr>,
#> #   basionymKey <int>, basionym <chr>, constituentKey <chr>, extinct <lgl>
```

The function `name_usage()` works with lots of different name endpoints in GBIF, listed at https://www.gbif.org/developer/species#nameUsages


```r
name_usage(key=3119195, language="FRENCH", data='vernacularNames')
#> Records returned [0] 
#> Args [offset=0, limit=100, language=FRENCH] 
#> # A tibble: 0 x 0
```

The function `name_backbone()` is used to search against the GBIF backbone taxonomy


```r
name_backbone(name='Helianthus', rank='genus', kingdom='plants')
#> # A tibble: 1 x 20
#>   usageKey scientificName canonicalName rank  status   confidence matchType
#> *    <int> <chr>          <chr>         <chr> <chr>         <int> <chr>    
#> 1  3119134 Helianthus L.  Helianthus    GENUS ACCEPTED         97 EXACT    
#> # … with 13 more variables: kingdom <chr>, phylum <chr>, order <chr>,
#> #   family <chr>, genus <chr>, kingdomKey <int>, phylumKey <int>,
#> #   classKey <int>, orderKey <int>, familyKey <int>, genusKey <int>,
#> #   synonym <lgl>, class <chr>
```

The function `name_suggest()` is optimized for speed, and gives back suggested names based on query parameters.


```r
head( name_suggest(q='Puma concolor') )
#> $data
#> # A tibble: 32 x 3
#>        key canonicalName                rank      
#>      <int> <chr>                        <chr>     
#>  1 2435099 Puma concolor                SPECIES   
#>  2 8860878 Puma concolor capricornensis SUBSPECIES
#>  3 6164618 Puma concolor browni         SUBSPECIES
#>  4 7193927 Puma concolor concolor       SUBSPECIES
#>  5 6164623 Puma concolor cabrerae       SUBSPECIES
#>  6 8944801 Puma concolor acrocodia      SUBSPECIES
#>  7 9104297 Puma concolor greeni         SUBSPECIES
#>  8 9045222 Puma concolor araucanus      SUBSPECIES
#>  9 6164608 Puma concolor californica    SUBSPECIES
#> 10 6164594 Puma concolor vancouverensis SUBSPECIES
#> # … with 22 more rows
#> 
#> $hierarchy
#> list()
```


## Single occurrence records

Get data for a single occurrence. Note that data is returned as a list, with slots for metadata and data.


```r
occ_get(key=855998194)
#> [[1]]
#> [[1]]$hierarchy
#>               name     key    rank
#> 1         Animalia       1 kingdom
#> 2         Chordata      44  phylum
#> 3         Mammalia     359   class
#> 4         Rodentia    1459   order
#> 5        Sciuridae    9456  family
#> 6          Sciurus 2437489   genus
#> 7 Sciurus vulgaris 8211070 species
#> 
#> [[1]]$media
#> [[1]]$media$`855998194`
#> [[1]]$media$`855998194`[[1]]
#> [[1]]$media$`855998194`[[1]][[1]]
#> [1] "none"
#> 
#> 
#> [[1]]$media$`855998194`$key
#> [1] "855998194"
#> 
#> [[1]]$media$`855998194`$species
#> [1] "Sciurus vulgaris"
#> 
#> [[1]]$media$`855998194`$decimalLatitude
#> [1] 58.40677
#> 
#> [[1]]$media$`855998194`$decimalLongitude
#> [1] 12.04386
#> 
#> [[1]]$media$`855998194`$country
#> [1] "Sweden"
#> 
#> 
#> 
#> [[1]]$data
#>         key                  scientificName decimalLatitude decimalLongitude
#> 1 855998194 Sciurus vulgaris Linnaeus, 1758        58.40677         12.04386
#>           issues
#> 1 cdround,gass84
```

Get many occurrences. `occ_get` is vectorized


```r
occ_get(key=c(855998194, 240713150))
#> [[1]]
#> [[1]]$hierarchy
#>               name     key    rank
#> 1         Animalia       1 kingdom
#> 2         Chordata      44  phylum
#> 3         Mammalia     359   class
#> 4         Rodentia    1459   order
#> 5        Sciuridae    9456  family
#> 6          Sciurus 2437489   genus
#> 7 Sciurus vulgaris 8211070 species
#> 
#> [[1]]$media
#> [[1]]$media$`855998194`
#> [[1]]$media$`855998194`[[1]]
#> [[1]]$media$`855998194`[[1]][[1]]
#> [1] "none"
#> 
#> 
#> [[1]]$media$`855998194`$key
#> [1] "855998194"
#> 
#> [[1]]$media$`855998194`$species
#> [1] "Sciurus vulgaris"
#> 
#> [[1]]$media$`855998194`$decimalLatitude
#> [1] 58.40677
#> 
#> [[1]]$media$`855998194`$decimalLongitude
#> [1] 12.04386
#> 
#> [[1]]$media$`855998194`$country
#> [1] "Sweden"
#> 
#> 
#> 
#> [[1]]$data
#>         key                  scientificName decimalLatitude decimalLongitude
#> 1 855998194 Sciurus vulgaris Linnaeus, 1758        58.40677         12.04386
#>           issues
#> 1 cdround,gass84
#> 
#> 
#> [[2]]
#> [[2]]$hierarchy
#>            name     key    rank
#> 1     Chromista       4 kingdom
#> 2  Foraminifera 8376456  phylum
#> 3  Monothalamea 7882876   class
#> 4  Astrorhizida 8142878   order
#> 5 Astrorhizidae 7747923  family
#> 6      Pelosina 7822114   genus
#> 
#> [[2]]$media
#> [[2]]$media$`240713150`
#> [[2]]$media$`240713150`[[1]]
#> [[2]]$media$`240713150`[[1]][[1]]
#> [1] "none"
#> 
#> 
#> [[2]]$media$`240713150`$key
#> [1] "240713150"
#> 
#> [[2]]$media$`240713150`$decimalLatitude
#> [1] -77.5667
#> 
#> [[2]]$media$`240713150`$decimalLongitude
#> [1] 163.583
#> 
#> [[2]]$media$`240713150`$country
#> [1] "Antarctica"
#> 
#> 
#> 
#> [[2]]$data
#>         key       scientificName decimalLatitude decimalLongitude
#> 1 240713150 Pelosina Brady, 1879        -77.5667          163.583
#>                   issues
#> 1 gass84,ambinst,colmano
```


## Search for occurrences

Note: The maximum number of records you can get with `occ_search()` and `occ_data()` is 100,000. See https://www.gbif.org/developer/occurrence

By default `occ_search()` returns a `dplyr` like output summary in which the data printed expands based on how much data is returned, and the size of your window. You can search by scientific name:


```r
occ_search(scientificName = "Ursus americanus", limit = 20)
#> Records found [20895] 
#> Records returned [20] 
#> No. unique hierarchies [1] 
#> No. media records [20] 
#> No. facets [0] 
#> Args [limit=20, offset=0, scientificName=Ursus americanus, fields=all] 
#> # A tibble: 20 x 79
#>    key    scientificName   decimalLatitude decimalLongitude issues  datasetKey  
#>    <chr>  <chr>                      <dbl>            <dbl> <chr>   <chr>       
#>  1 30179… Ursus americanu…            34.5           -120.  cdroun… 50c9509d-22…
#>  2 30179… Ursus americanu…            41.9            -73.5 cdroun… 50c9509d-22…
#>  3 30179… Ursus americanu…            38.4           -122.  cdroun… 50c9509d-22…
#>  4 30180… Ursus americanu…            37.5           -120.  cdroun… 50c9509d-22…
#>  5 30180… Ursus americanu…            37.5           -120.  gass84  50c9509d-22…
#>  6 30181… Ursus americanu…            42.7            -72.3 cdroun… 50c9509d-22…
#>  7 30181… Ursus americanu…            41.9            -73.6 cdroun… 50c9509d-22…
#>  8 30317… Ursus americanu…            42.2           -123.  cdroun… 50c9509d-22…
#>  9 30317… Ursus americanu…            25.2           -101.  cdroun… 50c9509d-22…
#> 10 30318… Ursus americanu…            42.3            -72.4 cdroun… 50c9509d-22…
#> 11 30318… Ursus americanu…            42.7            -73.2 cdroun… 50c9509d-22…
#> 12 30318… Ursus americanu…            51.9           -120.  cdroun… 50c9509d-22…
#> 13 30319… Ursus americanu…            42.2           -123.  cdroun… 50c9509d-22…
#> 14 30320… Ursus americanu…            43.6            -72.6 cdroun… 50c9509d-22…
#> 15 30320… Ursus americanu…            43.4            -71.9 cdroun… 50c9509d-22…
#> 16 30321… Ursus americanu…            45.3            -84.5 cdroun… 50c9509d-22…
#> 17 30321… Ursus americanu…            34.8           -120.  cdroun… 50c9509d-22…
#> 18 30391… Ursus americanu…            29.2            -81.6 cdroun… 50c9509d-22…
#> 19 30392… Ursus americanu…            48.5           -124.  gass84  50c9509d-22…
#> 20 30393… Ursus americanu…            29.3           -103.  cdroun… 50c9509d-22…
#> # … with 73 more variables: publishingOrgKey <chr>, installationKey <chr>,
#> #   publishingCountry <chr>, protocol <chr>, lastCrawled <chr>,
#> #   lastParsed <chr>, crawlId <int>, hostingOrganizationKey <chr>,
#> #   basisOfRecord <chr>, occurrenceStatus <chr>, taxonKey <int>,
#> #   kingdomKey <int>, phylumKey <int>, classKey <int>, orderKey <int>,
#> #   familyKey <int>, genusKey <int>, speciesKey <int>, acceptedTaxonKey <int>,
#> #   acceptedScientificName <chr>, kingdom <chr>, phylum <chr>, order <chr>,
#> #   family <chr>, genus <chr>, species <chr>, genericName <chr>,
#> #   specificEpithet <chr>, taxonRank <chr>, taxonomicStatus <chr>,
#> #   iucnRedListCategory <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <chr>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, identifiers <chr>,
#> #   facts <chr>, relations <chr>, isInCluster <lgl>, geodeticDatum <chr>,
#> #   class <chr>, countryCode <chr>, recordedByIDs <chr>, identifiedByIDs <chr>,
#> #   country <chr>, rightsHolder <chr>, identifier <chr>,
#> #   http...unknown.org.nick <chr>, informationWithheld <chr>,
#> #   verbatimEventDate <chr>, datasetName <chr>, verbatimLocality <chr>,
#> #   gbifID <chr>, collectionCode <chr>, occurrenceID <chr>, taxonID <chr>,
#> #   catalogNumber <chr>, recordedBy <chr>,
#> #   http...unknown.org.occurrenceDetails <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, occurrenceRemarks <chr>, identifiedBy <chr>,
#> #   identificationID <chr>, name <chr>
```

Or to be more precise, you can search for names first, make sure you have the right name, then pass the GBIF key to the `occ_search()` function:


```r
key <- name_suggest(q='Helianthus annuus', rank='species')$data$key[1]
occ_search(taxonKey=key, limit=20)
#> Records found [87309] 
#> Records returned [20] 
#> No. unique hierarchies [1] 
#> No. media records [20] 
#> No. facets [0] 
#> Args [limit=20, offset=0, taxonKey=9206251, fields=all] 
#> # A tibble: 20 x 86
#>    key    scientificName  decimalLatitude decimalLongitude issues  datasetKey   
#>    <chr>  <chr>                     <dbl>            <dbl> <chr>   <chr>        
#>  1 30179… Helianthus ann…           29.4            -98.5  "cdrou… 50c9509d-22c…
#>  2 30181… Helianthus ann…          -29.6             30.4  "cdrou… 50c9509d-22c…
#>  3 30181… Helianthus ann…           50.8              4.10 "cdrou… 50c9509d-22c…
#>  4 30316… Helianthus ann…           34.3           -118.   "cdrou… 50c9509d-22c…
#>  5 30316… Helianthus ann…            3.04           102.   "cdrou… 50c9509d-22c…
#>  6 30318… Helianthus ann…           35.1           -107.   "cdrou… 50c9509d-22c…
#>  7 30319… Helianthus ann…           33.9           -117.   "cdrou… 50c9509d-22c…
#>  8 30319… Helianthus ann…          -38.0            -59.2  "cdrou… 50c9509d-22c…
#>  9 30320… Helianthus ann…            8.60           -83.4  "cdrou… 50c9509d-22c…
#> 10 30320… Helianthus ann…           44.3            -78.4  "cdrou… 50c9509d-22c…
#> 11 30320… Helianthus ann…           28.3           -105.   "cdrou… 50c9509d-22c…
#> 12 30321… Helianthus ann…          -37.4            145.   "cdrou… 50c9509d-22c…
#> 13 30321… Helianthus ann…           45.4           -114.   "cdrou… 50c9509d-22c…
#> 14 30321… Helianthus ann…           29.6            -95.1  "cdrou… 50c9509d-22c…
#> 15 30321… Helianthus ann…           50.4             33.4  "cdrou… 50c9509d-22c…
#> 16 30321… Helianthus ann…          -32.2            142.   "cdrou… 50c9509d-22c…
#> 17 30391… Helianthus ann…          -25.4            -49.2  "cdrou… 50c9509d-22c…
#> 18 30391… Helianthus ann…           NA               NA    ""      50c9509d-22c…
#> 19 30392… Helianthus ann…           22.8            -98.5  "cdrou… 50c9509d-22c…
#> 20 30393… Helianthus ann…           25.8           -109.   "cdrou… 50c9509d-22c…
#> # … with 80 more variables: publishingOrgKey <chr>, installationKey <chr>,
#> #   publishingCountry <chr>, protocol <chr>, lastCrawled <chr>,
#> #   lastParsed <chr>, crawlId <int>, hostingOrganizationKey <chr>,
#> #   basisOfRecord <chr>, occurrenceStatus <chr>, taxonKey <int>,
#> #   kingdomKey <int>, phylumKey <int>, classKey <int>, orderKey <int>,
#> #   familyKey <int>, genusKey <int>, speciesKey <int>, acceptedTaxonKey <int>,
#> #   acceptedScientificName <chr>, kingdom <chr>, phylum <chr>, order <chr>,
#> #   family <chr>, genus <chr>, species <chr>, genericName <chr>,
#> #   specificEpithet <chr>, taxonRank <chr>, taxonomicStatus <chr>,
#> #   iucnRedListCategory <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <chr>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, identifiers <chr>,
#> #   facts <chr>, relations <chr>, isInCluster <lgl>, geodeticDatum <chr>,
#> #   class <chr>, countryCode <chr>, recordedByIDs <chr>, identifiedByIDs <chr>,
#> #   country <chr>, rightsHolder <chr>, identifier <chr>,
#> #   http...unknown.org.nick <chr>, verbatimEventDate <chr>, datasetName <chr>,
#> #   verbatimLocality <chr>, gbifID <chr>, collectionCode <chr>,
#> #   occurrenceID <chr>, taxonID <chr>, catalogNumber <chr>, recordedBy <chr>,
#> #   http...unknown.org.occurrenceDetails <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, reproductiveCondition <chr>,
#> #   identifiedBy <chr>, identificationID <chr>, name <chr>,
#> #   occurrenceRemarks <chr>, recordedByIDs.type <chr>,
#> #   recordedByIDs.value <chr>, identifiedByIDs.type <chr>,
#> #   identifiedByIDs.value <chr>, identificationRemarks <chr>, gadm <chr>,
#> #   informationWithheld <chr>
```

You can index to different parts of the oupu; here, the metadata:


```r
occ_search(taxonKey=key)$meta
#> $offset
#> [1] 300
#> 
#> $limit
#> [1] 200
#> 
#> $endOfRecords
#> [1] FALSE
#> 
#> $count
#> [1] 87309
```

You can choose what fields to return. This isn't passed on to the API query to GBIF as they don't allow that, but we filter out the columns before we give the data back to you.


```r
occ_search(scientificName = "Ursus americanus", fields=c('name','basisOfRecord','protocol'), limit = 20)
#> Records found [20895] 
#> Records returned [20] 
#> No. unique hierarchies [1] 
#> No. media records [20] 
#> No. facets [0] 
#> Args [limit=20, offset=0, scientificName=Ursus americanus,
#>      fields=name,basisOfRecord,protocol] 
#> # A tibble: 20 x 2
#>    protocol    basisOfRecord    
#>    <chr>       <chr>            
#>  1 DWC_ARCHIVE HUMAN_OBSERVATION
#>  2 DWC_ARCHIVE HUMAN_OBSERVATION
#>  3 DWC_ARCHIVE HUMAN_OBSERVATION
#>  4 DWC_ARCHIVE HUMAN_OBSERVATION
#>  5 DWC_ARCHIVE HUMAN_OBSERVATION
#>  6 DWC_ARCHIVE HUMAN_OBSERVATION
#>  7 DWC_ARCHIVE HUMAN_OBSERVATION
#>  8 DWC_ARCHIVE HUMAN_OBSERVATION
#>  9 DWC_ARCHIVE HUMAN_OBSERVATION
#> 10 DWC_ARCHIVE HUMAN_OBSERVATION
#> 11 DWC_ARCHIVE HUMAN_OBSERVATION
#> 12 DWC_ARCHIVE HUMAN_OBSERVATION
#> 13 DWC_ARCHIVE HUMAN_OBSERVATION
#> 14 DWC_ARCHIVE HUMAN_OBSERVATION
#> 15 DWC_ARCHIVE HUMAN_OBSERVATION
#> 16 DWC_ARCHIVE HUMAN_OBSERVATION
#> 17 DWC_ARCHIVE HUMAN_OBSERVATION
#> 18 DWC_ARCHIVE HUMAN_OBSERVATION
#> 19 DWC_ARCHIVE HUMAN_OBSERVATION
#> 20 DWC_ARCHIVE HUMAN_OBSERVATION
```

Most parameters are vectorized, so you can pass in more than one value:


```r
splist <- c('Cyanocitta stelleri', 'Junco hyemalis', 'Aix sponsa')
keys <- sapply(splist, function(x) name_suggest(x)$data$key[1], USE.NAMES=FALSE)
occ_search(taxonKey=keys, limit=5)
#> Records found [2482598 (1241473), 9362842 (6429075), 2498387 (2228924)] 
#> Records returned [2482598 (5), 9362842 (5), 2498387 (5)] 
#> No. unique hierarchies [2482598 (1), 9362842 (1), 2498387 (1)] 
#> No. media records [2482598 (5), 9362842 (5), 2498387 (5)] 
#> No. facets [2482598 (0), 9362842 (0), 2498387 (0)] 
#> Args [limit=5, offset=0, taxonKey=2482598,9362842,2498387, fields=all] 
#> 3 requests; First 10 rows of data from 2482598
#> 
#> # A tibble: 5 x 80
#>   key    scientificName    decimalLatitude decimalLongitude issues datasetKey   
#>   <chr>  <chr>                       <dbl>            <dbl> <chr>  <chr>        
#> 1 30179… Cyanocitta stell…            49.0            -123. cdrou… 50c9509d-22c…
#> 2 30179… Cyanocitta stell…            47.7            -122. cdrou… 50c9509d-22c…
#> 3 30179… Cyanocitta stell…            37.6            -122. cdrou… 50c9509d-22c…
#> 4 30179… Cyanocitta stell…            36.8            -122. cdrou… 50c9509d-22c…
#> 5 30179… Cyanocitta stell…            50.6            -115. cdrou… 50c9509d-22c…
#> # … with 74 more variables: publishingOrgKey <chr>, installationKey <chr>,
#> #   publishingCountry <chr>, protocol <chr>, lastCrawled <chr>,
#> #   lastParsed <chr>, crawlId <int>, hostingOrganizationKey <chr>,
#> #   basisOfRecord <chr>, occurrenceStatus <chr>, taxonKey <int>,
#> #   kingdomKey <int>, phylumKey <int>, classKey <int>, orderKey <int>,
#> #   familyKey <int>, genusKey <int>, speciesKey <int>, acceptedTaxonKey <int>,
#> #   acceptedScientificName <chr>, kingdom <chr>, phylum <chr>, order <chr>,
#> #   family <chr>, genus <chr>, species <chr>, genericName <chr>,
#> #   specificEpithet <chr>, taxonRank <chr>, taxonomicStatus <chr>,
#> #   iucnRedListCategory <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <chr>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, identifiers <chr>,
#> #   facts <chr>, relations <chr>, isInCluster <lgl>, geodeticDatum <chr>,
#> #   class <chr>, countryCode <chr>, recordedByIDs <chr>, identifiedByIDs <chr>,
#> #   country <chr>, rightsHolder <chr>, identifier <chr>,
#> #   http...unknown.org.nick <chr>, verbatimEventDate <chr>, datasetName <chr>,
#> #   verbatimLocality <chr>, gbifID <chr>, collectionCode <chr>,
#> #   occurrenceID <chr>, taxonID <chr>, catalogNumber <chr>, recordedBy <chr>,
#> #   http...unknown.org.occurrenceDetails <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, identifiedBy <chr>, identificationID <chr>,
#> #   name <chr>, occurrenceRemarks <chr>, gadm <chr>, informationWithheld <chr>
```


********************

## Maps

Using thet GBIF map web tile service, making a raster and visualizing it.


```r
x <- map_fetch(taxonKey = 2480498, year = 2000:2017)
library(raster)
plot(x)
```

![](../man/figures/rgbif_vign_1.png)

[gbifapi]: https://www.gbif.org/developer/summary
---
title: Cleaning data using GBIF issues
author: Scott Chamberlain
date: "2021-06-11"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{cleaning}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`rgbif` now has the ability to clean data retrieved from GBIF based on GBIF issues. These issues are returned in data retrieved from GBIF, e.g., through the `occ_search()` function. Inspired by `magrittr`, we've setup a workflow for cleaning data based on using the operator `%>%`. You don't have to use it, but as we show below, it can make the process quite easy.

Note that you can also query based on issues, e.g., `occ_search(taxonKey=1, issue='DEPTH_UNLIKELY')`. However, we imagine it's more likely that you want to search for occurrences based on a taxonomic name, or geographic area, not based on issues, so it makes sense to pull data down, then clean as needed using the below workflow with `occ_issues()`.

Note that `occ_issues()` only affects the data element in the gbif class that is returned from a call to `occ_search()`. Maybe in a future version we will remove the associated records from the hierarchy and media elements as they are remove from the data element.

`occ_issues()` also works with data from `occ_download()`.

## Get rgbif

Install from CRAN


```r
install.packages("rgbif")
```

Or install the development version from GitHub


```r
remotes::install_github("ropensci/rgbif")
```

Load rgbif


```r
library('rgbif')
```

## Get some data

Get taxon key for _Helianthus annuus_


```r
(key <- name_suggest(q='Helianthus annuus', rank='species')$data$key[1])
#> [1] 9206251
```

Then pass to `occ_search()`


```r
(res <- occ_search(taxonKey=key, limit=100))
#> Records found [87309] 
#> Records returned [100] 
#> No. unique hierarchies [1] 
#> No. media records [100] 
#> No. facets [0] 
#> Args [limit=100, offset=0, taxonKey=9206251, fields=all] 
#> # A tibble: 100 x 101
#>    key    scientificName  decimalLatitude decimalLongitude issues  datasetKey   
#>    <chr>  <chr>                     <dbl>            <dbl> <chr>   <chr>        
#>  1 30179… Helianthus ann…           29.4            -98.5  cdroun… 50c9509d-22c…
#>  2 30181… Helianthus ann…          -29.6             30.4  cdroun… 50c9509d-22c…
#>  3 30181… Helianthus ann…           50.8              4.10 cdroun… 50c9509d-22c…
#>  4 30316… Helianthus ann…           34.3           -118.   cdroun… 50c9509d-22c…
#>  5 30316… Helianthus ann…            3.04           102.   cdroun… 50c9509d-22c…
#>  6 30318… Helianthus ann…           35.1           -107.   cdroun… 50c9509d-22c…
#>  7 30319… Helianthus ann…           33.9           -117.   cdroun… 50c9509d-22c…
#>  8 30319… Helianthus ann…          -38.0            -59.2  cdroun… 50c9509d-22c…
#>  9 30320… Helianthus ann…            8.60           -83.4  cdroun… 50c9509d-22c…
#> 10 30320… Helianthus ann…           44.3            -78.4  cdroun… 50c9509d-22c…
#> # … with 90 more rows, and 95 more variables: publishingOrgKey <chr>,
#> #   installationKey <chr>, publishingCountry <chr>, protocol <chr>,
#> #   lastCrawled <chr>, lastParsed <chr>, crawlId <int>,
#> #   hostingOrganizationKey <chr>, basisOfRecord <chr>, occurrenceStatus <chr>,
#> #   taxonKey <int>, kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   orderKey <int>, familyKey <int>, genusKey <int>, speciesKey <int>,
#> #   acceptedTaxonKey <int>, acceptedScientificName <chr>, kingdom <chr>,
#> #   phylum <chr>, order <chr>, family <chr>, genus <chr>, species <chr>,
#> #   genericName <chr>, specificEpithet <chr>, taxonRank <chr>,
#> #   taxonomicStatus <chr>, iucnRedListCategory <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <chr>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, identifiers <chr>,
#> #   facts <chr>, relations <chr>, isInCluster <lgl>, geodeticDatum <chr>,
#> #   class <chr>, countryCode <chr>, recordedByIDs <chr>, identifiedByIDs <chr>,
#> #   country <chr>, rightsHolder <chr>, identifier <chr>,
#> #   http...unknown.org.nick <chr>, verbatimEventDate <chr>, datasetName <chr>,
#> #   verbatimLocality <chr>, gbifID <chr>, collectionCode <chr>,
#> #   occurrenceID <chr>, taxonID <chr>, catalogNumber <chr>, recordedBy <chr>,
#> #   http...unknown.org.occurrenceDetails <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, reproductiveCondition <chr>,
#> #   identifiedBy <chr>, identificationID <chr>, name <chr>,
#> #   occurrenceRemarks <chr>, recordedByIDs.type <chr>,
#> #   recordedByIDs.value <chr>, identifiedByIDs.type <chr>,
#> #   identifiedByIDs.value <chr>, identificationRemarks <chr>, gadm <chr>,
#> #   informationWithheld <chr>, continent <chr>, municipality <chr>,
#> #   county <chr>, identificationVerificationStatus <chr>, language <chr>,
#> #   type <chr>, vernacularName <chr>, taxonConceptID <chr>, endDayOfYear <chr>,
#> #   locality <chr>, startDayOfYear <chr>, datasetID <chr>, accessRights <chr>,
#> #   bibliographicCitation <chr>, higherClassification <chr>
```

## Examine issues

The dataset `gbifissues` can be retrieved using the function `gbif_issues()`. The dataset's first column `code` is a code that is used by default in the results from `occ_search()`, while the second column `issue` is the full issue name given by GBIF. The third column is a full description of the issue.


```r
head(gbif_issues())
#>    code                              issue
#> 1   bri            BASIS_OF_RECORD_INVALID
#> 2   ccm         CONTINENT_COUNTRY_MISMATCH
#> 3   cdc CONTINENT_DERIVED_FROM_COORDINATES
#> 4 conti                  CONTINENT_INVALID
#> 5  cdiv                 COORDINATE_INVALID
#> 6 cdout            COORDINATE_OUT_OF_RANGE
#>                                                                                                    description
#> 1 The given basis of record is impossible to interpret or seriously different from the recommended vocabulary.
#> 2                                                       The interpreted continent and country do not match up.
#> 3                  The interpreted continent is based on the coordinates, not the verbatim string information.
#> 4                                                                      Uninterpretable continent values found.
#> 5                                      Coordinate value given in some form but GBIF is unable to interpret it.
#> 6                                        Coordinate has invalid lat/lon values out of their decimal max range.
#>         type
#> 1 occurrence
#> 2 occurrence
#> 3 occurrence
#> 4 occurrence
#> 5 occurrence
#> 6 occurrence
```

You can query to get certain issues


```r
gbif_issues()[ gbif_issues()$code %in% c('cdround','cudc','gass84','txmathi'), ]
#>       code                            issue
#> 10 cdround               COORDINATE_ROUNDED
#> 12    cudc COUNTRY_DERIVED_FROM_COORDINATES
#> 23  gass84     GEODETIC_DATUM_ASSUMED_WGS84
#> 39 txmathi           TAXON_MATCH_HIGHERRANK
#>                                                                                                                                 description
#> 10                                                                                  Original coordinate modified by rounding to 5 decimals.
#> 12                                                The interpreted country is based on the coordinates, not the verbatim string information.
#> 23 Indicating that the interpreted coordinates assume they are based on WGS84 datum as the datum was either not indicated or interpretable.
#> 39                                        Matching to the taxonomic backbone can only be done on a higher rank and not the scientific name.
#>          type
#> 10 occurrence
#> 12 occurrence
#> 23 occurrence
#> 39 occurrence
```

The code `cdround` represents the GBIF issue `COORDINATE_ROUNDED`, which means that

> Original coordinate modified by rounding to 5 decimals.

The content for this information comes from https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html

## Parse data based on issues

Now that we know a bit about GBIF issues, you can parse your data based on issues. Using the data generated above, and using the function `%>%` imported from `magrittr`, we can get only data with the issue `gass84`, or `GEODETIC_DATUM_ASSUMED_WGS84` (Note how the records returned goes down to 98 instead of the initial 100).


```r
res %>%
  occ_issues(gass84)
#> Records found [87309] 
#> Records returned [99] 
#> No. unique hierarchies [1] 
#> No. media records [100] 
#> No. facets [0] 
#> Args [limit=100, offset=0, taxonKey=9206251, fields=all] 
#> # A tibble: 99 x 101
#>    key    scientificName  decimalLatitude decimalLongitude issues  datasetKey   
#>    <chr>  <chr>                     <dbl>            <dbl> <chr>   <chr>        
#>  1 30179… Helianthus ann…           29.4            -98.5  cdroun… 50c9509d-22c…
#>  2 30181… Helianthus ann…          -29.6             30.4  cdroun… 50c9509d-22c…
#>  3 30181… Helianthus ann…           50.8              4.10 cdroun… 50c9509d-22c…
#>  4 30316… Helianthus ann…           34.3           -118.   cdroun… 50c9509d-22c…
#>  5 30316… Helianthus ann…            3.04           102.   cdroun… 50c9509d-22c…
#>  6 30318… Helianthus ann…           35.1           -107.   cdroun… 50c9509d-22c…
#>  7 30319… Helianthus ann…           33.9           -117.   cdroun… 50c9509d-22c…
#>  8 30319… Helianthus ann…          -38.0            -59.2  cdroun… 50c9509d-22c…
#>  9 30320… Helianthus ann…            8.60           -83.4  cdroun… 50c9509d-22c…
#> 10 30320… Helianthus ann…           44.3            -78.4  cdroun… 50c9509d-22c…
#> # … with 89 more rows, and 95 more variables: publishingOrgKey <chr>,
#> #   installationKey <chr>, publishingCountry <chr>, protocol <chr>,
#> #   lastCrawled <chr>, lastParsed <chr>, crawlId <int>,
#> #   hostingOrganizationKey <chr>, basisOfRecord <chr>, occurrenceStatus <chr>,
#> #   taxonKey <int>, kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   orderKey <int>, familyKey <int>, genusKey <int>, speciesKey <int>,
#> #   acceptedTaxonKey <int>, acceptedScientificName <chr>, kingdom <chr>,
#> #   phylum <chr>, order <chr>, family <chr>, genus <chr>, species <chr>,
#> #   genericName <chr>, specificEpithet <chr>, taxonRank <chr>,
#> #   taxonomicStatus <chr>, iucnRedListCategory <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <chr>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, identifiers <chr>,
#> #   facts <chr>, relations <chr>, isInCluster <lgl>, geodeticDatum <chr>,
#> #   class <chr>, countryCode <chr>, recordedByIDs <chr>, identifiedByIDs <chr>,
#> #   country <chr>, rightsHolder <chr>, identifier <chr>,
#> #   http...unknown.org.nick <chr>, verbatimEventDate <chr>, datasetName <chr>,
#> #   verbatimLocality <chr>, gbifID <chr>, collectionCode <chr>,
#> #   occurrenceID <chr>, taxonID <chr>, catalogNumber <chr>, recordedBy <chr>,
#> #   http...unknown.org.occurrenceDetails <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, reproductiveCondition <chr>,
#> #   identifiedBy <chr>, identificationID <chr>, name <chr>,
#> #   occurrenceRemarks <chr>, recordedByIDs.type <chr>,
#> #   recordedByIDs.value <chr>, identifiedByIDs.type <chr>,
#> #   identifiedByIDs.value <chr>, identificationRemarks <chr>, gadm <chr>,
#> #   informationWithheld <chr>, continent <chr>, municipality <chr>,
#> #   county <chr>, identificationVerificationStatus <chr>, language <chr>,
#> #   type <chr>, vernacularName <chr>, taxonConceptID <chr>, endDayOfYear <chr>,
#> #   locality <chr>, startDayOfYear <chr>, datasetID <chr>, accessRights <chr>,
#> #   bibliographicCitation <chr>, higherClassification <chr>
```

Note also that we've set up `occ_issues()` so that you can pass in issue names without having to quote them, thereby speeding up data cleaning.

Next, we can remove data with certain issues just as easily by using a `-` sign in front of the variable, like this, removing data with issues `depunl` and `mdatunl`.


```r
res %>%
  occ_issues(-depunl, -mdatunl)
#> Records found [87309] 
#> Records returned [100] 
#> No. unique hierarchies [1] 
#> No. media records [100] 
#> No. facets [0] 
#> Args [limit=100, offset=0, taxonKey=9206251, fields=all] 
#> # A tibble: 100 x 101
#>    key    scientificName  decimalLatitude decimalLongitude issues  datasetKey   
#>    <chr>  <chr>                     <dbl>            <dbl> <chr>   <chr>        
#>  1 30179… Helianthus ann…           29.4            -98.5  cdroun… 50c9509d-22c…
#>  2 30181… Helianthus ann…          -29.6             30.4  cdroun… 50c9509d-22c…
#>  3 30181… Helianthus ann…           50.8              4.10 cdroun… 50c9509d-22c…
#>  4 30316… Helianthus ann…           34.3           -118.   cdroun… 50c9509d-22c…
#>  5 30316… Helianthus ann…            3.04           102.   cdroun… 50c9509d-22c…
#>  6 30318… Helianthus ann…           35.1           -107.   cdroun… 50c9509d-22c…
#>  7 30319… Helianthus ann…           33.9           -117.   cdroun… 50c9509d-22c…
#>  8 30319… Helianthus ann…          -38.0            -59.2  cdroun… 50c9509d-22c…
#>  9 30320… Helianthus ann…            8.60           -83.4  cdroun… 50c9509d-22c…
#> 10 30320… Helianthus ann…           44.3            -78.4  cdroun… 50c9509d-22c…
#> # … with 90 more rows, and 95 more variables: publishingOrgKey <chr>,
#> #   installationKey <chr>, publishingCountry <chr>, protocol <chr>,
#> #   lastCrawled <chr>, lastParsed <chr>, crawlId <int>,
#> #   hostingOrganizationKey <chr>, basisOfRecord <chr>, occurrenceStatus <chr>,
#> #   taxonKey <int>, kingdomKey <int>, phylumKey <int>, classKey <int>,
#> #   orderKey <int>, familyKey <int>, genusKey <int>, speciesKey <int>,
#> #   acceptedTaxonKey <int>, acceptedScientificName <chr>, kingdom <chr>,
#> #   phylum <chr>, order <chr>, family <chr>, genus <chr>, species <chr>,
#> #   genericName <chr>, specificEpithet <chr>, taxonRank <chr>,
#> #   taxonomicStatus <chr>, iucnRedListCategory <chr>, dateIdentified <chr>,
#> #   coordinateUncertaintyInMeters <dbl>, stateProvince <chr>, year <int>,
#> #   month <int>, day <int>, eventDate <chr>, modified <chr>,
#> #   lastInterpreted <chr>, references <chr>, license <chr>, identifiers <chr>,
#> #   facts <chr>, relations <chr>, isInCluster <lgl>, geodeticDatum <chr>,
#> #   class <chr>, countryCode <chr>, recordedByIDs <chr>, identifiedByIDs <chr>,
#> #   country <chr>, rightsHolder <chr>, identifier <chr>,
#> #   http...unknown.org.nick <chr>, verbatimEventDate <chr>, datasetName <chr>,
#> #   verbatimLocality <chr>, gbifID <chr>, collectionCode <chr>,
#> #   occurrenceID <chr>, taxonID <chr>, catalogNumber <chr>, recordedBy <chr>,
#> #   http...unknown.org.occurrenceDetails <chr>, institutionCode <chr>,
#> #   rights <chr>, eventTime <chr>, reproductiveCondition <chr>,
#> #   identifiedBy <chr>, identificationID <chr>, name <chr>,
#> #   occurrenceRemarks <chr>, recordedByIDs.type <chr>,
#> #   recordedByIDs.value <chr>, identifiedByIDs.type <chr>,
#> #   identifiedByIDs.value <chr>, identificationRemarks <chr>, gadm <chr>,
#> #   informationWithheld <chr>, continent <chr>, municipality <chr>,
#> #   county <chr>, identificationVerificationStatus <chr>, language <chr>,
#> #   type <chr>, vernacularName <chr>, taxonConceptID <chr>, endDayOfYear <chr>,
#> #   locality <chr>, startDayOfYear <chr>, datasetID <chr>, accessRights <chr>,
#> #   bibliographicCitation <chr>, higherClassification <chr>
```

## Expand issue codes to full names

Another thing we can do with `occ_issues()` is go from issue codes to full issue names in case you want those in your dataset (here, showing only a few columns to see the data better for this demo):


```r
out <- res %>% occ_issues(mutate = "expand")
head(out$data[,c(1,5)])
#> # A tibble: 6 x 2
#>   key        issues                                         
#>   <chr>      <chr>                                          
#> 1 3017947125 COORDINATE_ROUNDED,GEODETIC_DATUM_ASSUMED_WGS84
#> 2 3018105046 COORDINATE_ROUNDED,GEODETIC_DATUM_ASSUMED_WGS84
#> 3 3018133231 COORDINATE_ROUNDED,GEODETIC_DATUM_ASSUMED_WGS84
#> 4 3031652179 COORDINATE_ROUNDED,GEODETIC_DATUM_ASSUMED_WGS84
#> 5 3031677301 COORDINATE_ROUNDED,GEODETIC_DATUM_ASSUMED_WGS84
#> 6 3031889761 COORDINATE_ROUNDED,GEODETIC_DATUM_ASSUMED_WGS84
```


## Add columns

Sometimes you may want to have each type of issue as a separate column.

Split out each issue type into a separate column, with number of columns equal to number of issue types


```r
out <- res %>% occ_issues(mutate = "split")
head(out$data[,c(1,5:10)])
#> # A tibble: 6 x 7
#>   name      cdround gass84 scientificName   datasetKey      publishingOrgKey    
#>   <chr>     <chr>   <chr>  <chr>            <chr>           <chr>               
#> 1 Helianth… y       y      Helianthus annu… 50c9509d-22c7-… 28eb1a3f-1c15-4a95-…
#> 2 Helianth… y       y      Helianthus annu… 50c9509d-22c7-… 28eb1a3f-1c15-4a95-…
#> 3 Helianth… y       y      Helianthus annu… 50c9509d-22c7-… 28eb1a3f-1c15-4a95-…
#> 4 Helianth… y       y      Helianthus annu… 50c9509d-22c7-… 28eb1a3f-1c15-4a95-…
#> 5 Helianth… y       y      Helianthus annu… 50c9509d-22c7-… 28eb1a3f-1c15-4a95-…
#> 6 Helianth… y       y      Helianthus annu… 50c9509d-22c7-… 28eb1a3f-1c15-4a95-…
#> # … with 1 more variable: installationKey <chr>
```

## Expand and add columns

Or you can expand each issue type into its full name, and split each issue into a separate column.


```r
out <- res %>% occ_issues(mutate = "split_expand")
head(out$data[,c(1,5:10)])
#> # A tibble: 6 x 7
#>   name      COORDINATE_ROUND… GEODETIC_DATUM_AS… scientificName  datasetKey     
#>   <chr>     <chr>             <chr>              <chr>           <chr>          
#> 1 Helianth… y                 y                  Helianthus ann… 50c9509d-22c7-…
#> 2 Helianth… y                 y                  Helianthus ann… 50c9509d-22c7-…
#> 3 Helianth… y                 y                  Helianthus ann… 50c9509d-22c7-…
#> 4 Helianth… y                 y                  Helianthus ann… 50c9509d-22c7-…
#> 5 Helianth… y                 y                  Helianthus ann… 50c9509d-22c7-…
#> 6 Helianth… y                 y                  Helianthus ann… 50c9509d-22c7-…
#> # … with 2 more variables: publishingOrgKey <chr>, installationKey <chr>
```

## Wrap up

We hope this helps users get just the data they want, and nothing more. Let us know if you have feedback on data cleaning functionality in `rgbif` at _info@ropensci.org_ or at [https://github.com/ropensci/rgbif/issues](https://github.com/ropensci/rgbif/issues).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgbif-package.r
\docType{package}
\name{rgbif-package}
\alias{rgbif-package}
\alias{rgbif}
\title{Interface to the Global Biodiversity Information Facility API.}
\description{
rgbif: A programmatic interface to the Web Service methods
provided by the Global Biodiversity Information Facility.
}
\note{
See \link{many-values} for discussion of how functions vary in how
they accept values (single vs. many for the same HTTP request vs. many
for different HTTP requests)
}
\section{About}{


This package gives you access to data from GBIF \url{https://www.gbif.org/}
via their API.
}

\section{A note about the old GBIF API}{


The old GBIF API is now defunct - that is, not available anymore. We used
to have functions that worked with the old API, but those functions are
now not available anymore because GBIF made the old API defunct.
}

\section{Documentation for the GBIF API}{

\itemize{
\item summary \url{https://www.gbif.org/developer/summary} - Summary of
the GBIF API
\item registry \url{https://www.gbif.org/developer/registry} - Metadata
on datasets, and contributing organizations
\item species names \url{https://www.gbif.org/developer/species} - Species
names and metadata
\item occurrences \url{https://www.gbif.org/developer/occurrence} -
Occurrences
\item maps \url{https://www.gbif.org/developer/maps} - Maps - these APIs
are not implemented in \pkg{rgbif}, and are meant more for intergration
with web based maps.
}
}

\author{
Scott Chamberlain

Karthik Ram

Dan Mcglinn

Vijay Barve
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networks.r
\name{networks}
\alias{networks}
\title{Networks metadata.}
\usage{
networks(
  data = "all",
  uuid = NULL,
  query = NULL,
  identifier = NULL,
  identifierType = NULL,
  limit = 100,
  start = NULL,
  curlopts = list()
)
}
\arguments{
\item{data}{The type of data to get. One or more of: 'contact', 'endpoint',
'identifier', 'tag', 'machineTag', 'comment', 'constituents', or the
special 'all'. Default: \code{'all'}}

\item{uuid}{UUID of the data network provider. This must be specified if
data is anything other than 'all'. Only 1 can be passed in}

\item{query}{Query nodes. Only used when \code{data='all'}. Ignored
otherwise.}

\item{identifier}{The value for this parameter can be a simple string or
integer, e.g. \code{identifier=120}. This parameter doesn't seem to work right
now.}

\item{identifierType}{Used in combination with the identifier parameter to
filter identifiers by identifier type. See details. This parameter doesn't
seem to work right now.}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\description{
Networks metadata.
}
\details{
identifierType options:

\itemize{
\item {DOI} No description.
\item {FTP} No description.
\item {GBIF_NODE} Identifies the node (e.g: \code{DK} for Denmark, \code{sp2000}
for Species 2000).
\item {GBIF_PARTICIPANT} Participant identifier from the GBIF IMS
Filemaker system.
\item {GBIF_PORTAL} Indicates the identifier originated from an
auto_increment column in the portal.data_provider or portal.data_resource
table respectively.
\item {HANDLER} No description.
\item {LSID} Reference controlled by a separate system, used for example
by DOI.
\item {SOURCE_ID} No description.
\item {UNKNOWN} No description.
\item {URI} No description.
\item {URL} No description.
\item {UUID} No description.
}
}
\examples{
\dontrun{
networks()
networks(uuid='2b7c7b4f-4d4f-40d3-94de-c28b6fa054a6')

# curl options
networks(curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#networks}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_data.R
\name{occ_data}
\alias{occ_data}
\title{Search for GBIF occurrences - simplified for speed}
\usage{
occ_data(
  taxonKey = NULL,
  scientificName = NULL,
  country = NULL,
  publishingCountry = NULL,
  hasCoordinate = NULL,
  typeStatus = NULL,
  recordNumber = NULL,
  lastInterpreted = NULL,
  continent = NULL,
  geometry = NULL,
  geom_big = "asis",
  geom_size = 40,
  geom_n = 10,
  recordedBy = NULL,
  recordedByID = NULL,
  identifiedByID = NULL,
  basisOfRecord = NULL,
  datasetKey = NULL,
  eventDate = NULL,
  catalogNumber = NULL,
  year = NULL,
  month = NULL,
  decimalLatitude = NULL,
  decimalLongitude = NULL,
  elevation = NULL,
  depth = NULL,
  institutionCode = NULL,
  collectionCode = NULL,
  hasGeospatialIssue = NULL,
  issue = NULL,
  search = NULL,
  mediaType = NULL,
  subgenusKey = NULL,
  repatriated = NULL,
  phylumKey = NULL,
  kingdomKey = NULL,
  classKey = NULL,
  orderKey = NULL,
  familyKey = NULL,
  genusKey = NULL,
  establishmentMeans = NULL,
  protocol = NULL,
  license = NULL,
  organismId = NULL,
  publishingOrg = NULL,
  stateProvince = NULL,
  waterBody = NULL,
  locality = NULL,
  limit = 500,
  start = 0,
  skip_validate = TRUE,
  curlopts = list()
)
}
\arguments{
\item{taxonKey}{(numeric) A taxon key from the GBIF backbone. All included and synonym taxa
are included in the search, so a search for aves with taxononKey=212
(i.e. /occurrence/search?taxonKey=212) will match all birds, no matter which
species. You can pass many keys by passing occ_search in a call to an
lapply-family function (see last example below).}

\item{scientificName}{A scientific name from the GBIF backbone. All included and synonym
taxa are included in the search.}

\item{country}{The 2-letter country code (as per ISO-3166-1) of the country in
which the occurrence was recorded. See here
https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}

\item{publishingCountry}{The 2-letter country code (as per ISO-3166-1) of the
country in which the occurrence was recorded.}

\item{hasCoordinate}{(logical) Return only occurence records with lat/long data (\code{TRUE}) or
all records (\code{FALSE}, default).}

\item{typeStatus}{Type status of the specimen. One of many options. See \code{?typestatus}}

\item{recordNumber}{Number recorded by collector of the data, different from GBIF record
number. See http://rs.tdwg.org/dwc/terms/#recordNumber for more info}

\item{lastInterpreted}{Date the record was last modified in GBIF, in ISO 8601 format:
yyyy, yyyy-MM, yyyy-MM-dd, or MM-dd.  Supports range queries, smaller,larger (e.g.,
'1990,1991', whereas '1991,1990' wouldn't work)}

\item{continent}{Continent. One of africa, antarctica, asia, europe, north_america
(North America includes the Caribbean and reachies down and includes Panama), oceania,
or south_america}

\item{geometry}{Searches for occurrences inside a polygon described in Well Known
Text (WKT) format. A WKT shape written as either POINT, LINESTRING, LINEARRING
POLYGON, or MULTIPOLYGON. Example of a polygon: POLYGON((30.1 10.1, 20, 20 40, 40 40, 30.1 10.1))
would be queried as http://bit.ly/1BzNwDq See also the section \strong{WKT} below.}

\item{geom_big}{(character) One of "axe", "bbox", or "asis" (default). See Details.}

\item{geom_size}{(integer) An integer indicating size of the cell. Default: 40. See Details.}

\item{geom_n}{(integer) An integer indicating number of cells in each dimension. Default: 10.
See Details.}

\item{recordedBy}{The person who recorded the occurrence.}

\item{recordedByID}{(character) Identifier (e.g. ORCID) for the person who
recorded the occurrence}

\item{identifiedByID}{(character) Identifier (e.g. ORCID) for the person who
provided the taxonomic identification of the occurrence.}

\item{basisOfRecord}{Basis of record, as defined in our BasisOfRecord enum here
https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html
Acceptable values are:
\itemize{
\item FOSSIL_SPECIMEN An occurrence record describing a fossilized specimen.
\item HUMAN_OBSERVATION An occurrence record describing an observation made by
one or more people.
\item LITERATURE An occurrence record based on literature alone.
\item LIVING_SPECIMEN An occurrence record describing a living specimen, e.g.
\item MACHINE_OBSERVATION An occurrence record describing an observation made
by a machine.
\item OBSERVATION An occurrence record describing an observation.
\item PRESERVED_SPECIMEN An occurrence record describing a preserved specimen.
\item UNKNOWN Unknown basis for the record.
}}

\item{datasetKey}{The occurrence dataset key (a uuid)}

\item{eventDate}{Occurrence date in ISO 8601 format: yyyy, yyyy-MM, yyyy-MM-dd, or
MM-dd. Supports range queries, smaller,larger (e.g., '1990,1991', whereas '1991,1990'
wouldn't work)}

\item{catalogNumber}{An identifier of any form assigned by the source within a
physical collection or digital dataset for the record which may not unique,
but should be fairly unique in combination with the institution and collection code.}

\item{year}{The 4 digit year. A year of 98 will be interpreted as AD 98. Supports range queries,
smaller,larger (e.g., '1990,1991', whereas '1991,1990' wouldn't work)}

\item{month}{The month of the year, starting with 1 for January. Supports range queries,
smaller,larger (e.g., '1,2', whereas '2,1' wouldn't work)}

\item{decimalLatitude}{Latitude in decimals between -90 and 90 based on WGS 84.
Supports range queries, smaller,larger (e.g., '25,30', whereas '30,25' wouldn't work)}

\item{decimalLongitude}{Longitude in decimals between -180 and 180 based on WGS 84.
Supports range queries (e.g., '-0.4,-0.2', whereas '-0.2,-0.4' wouldn't work).}

\item{elevation}{Elevation in meters above sea level. Supports range queries, smaller,larger
(e.g., '5,30', whereas '30,5' wouldn't work)}

\item{depth}{Depth in meters relative to elevation. For example 10 meters below a
lake surface with given elevation. Supports range queries, smaller,larger (e.g., '5,30',
whereas '30,5' wouldn't work)}

\item{institutionCode}{An identifier of any form assigned by the source to identify
the institution the record belongs to. Not guaranteed to be que.}

\item{collectionCode}{An identifier of any form assigned by the source to identify
the physical collection or digital dataset uniquely within the text of an institution.}

\item{hasGeospatialIssue}{(logical) Includes/excludes occurrence records which contain spatial
issues (as determined in our record interpretation), i.e. \code{hasGeospatialIssue=TRUE}
returns only those records with spatial issues while \code{hasGeospatialIssue=FALSE} includes
only records without spatial issues. The absence of this parameter returns any
record with or without spatial issues.}

\item{issue}{(character) One or more of many possible issues with each occurrence record. See
Details. Issues passed to this parameter filter results by the issue.}

\item{search}{Query terms. The value for this parameter can be a simple word or a phrase.}

\item{mediaType}{Media type. Default is NULL, so no filtering on mediatype. Options:
NULL, 'MovingImage', 'Sound', and 'StillImage'.}

\item{subgenusKey}{(numeric) Subgenus classification key.}

\item{repatriated}{(character) Searches for records whose publishing country
is different to the country where the record was recorded in.}

\item{phylumKey}{(numeric) Phylum classification key.}

\item{kingdomKey}{(numeric) Kingdom classification key.}

\item{classKey}{(numeric) Class classification key.}

\item{orderKey}{(numeric) Order classification key.}

\item{familyKey}{(numeric) Family classification key.}

\item{genusKey}{(numeric) Genus classification key.}

\item{establishmentMeans}{(character) EstablishmentMeans, possible values
include: INTRODUCED, INVASIVE, MANAGED, NATIVE, NATURALISED, UNCERTAIN}

\item{protocol}{(character) Protocol or mechanism used to provide the
occurrence record. See Details for possible values}

\item{license}{(character) The type license applied to the dataset or record.
Possible values: CC0_1_0, CC_BY_4_0, CC_BY_NC_4_0, UNSPECIFIED, and
UNSUPPORTED}

\item{organismId}{(numeric) An identifier for the Organism instance (as
opposed to a particular digital record of the Organism). May be a globally
unique identifier or an identifier specific to the data set.}

\item{publishingOrg}{(character) The publishing organization key (a UUID).}

\item{stateProvince}{(character) The name of the next smaller administrative
region than country (state, province, canton, department, region, etc.) in
which the Location occurs.}

\item{waterBody}{(character) The name of the water body in which the
locations occur}

\item{locality}{(character) The specific description of the place.}

\item{limit}{Number of records to return. Default: 500. Note that the per
request maximum is 300, but since we set it at 500 for the function, we
do two requests to get you the 500 records (if there are that many).
Note that there is a hard maximum of 100,000, which is calculated as the
\code{limit+start}, so \code{start=99,000} and \code{limit=2000} won't work}

\item{start}{Record number to start at. Use in combination with limit to
page through results. Note that we do the paging internally for you, but
you can manually set the \code{start} parameter}

\item{skip_validate}{(logical) whether to skip \code{wellknown::validate_wkt}
call or not. passed down to \code{\link[=check_wkt]{check_wkt()}}. Default: \code{TRUE}}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
An object of class \code{gbif_data}, which is a S3 class list, with
slots for metadata (\code{meta}) and the occurrence data itself
(\code{data}), and with attributes listing the user supplied arguments
and whether it was a "single" or "many" search; that is, if you supply
two values of the \code{datasetKey} parameter to searches are done, and
it's a "many". \code{meta} is a list of length four with offset, limit,
endOfRecords and count fields. \code{data} is a tibble (aka data.frame)
}
\description{
Search for GBIF occurrences - simplified for speed
}
\note{
Maximum number of records you can get with this function is 100,000.
See https://www.gbif.org/developer/occurrence
}
\section{protocol parameter options}{

\itemize{
\item BIOCASE - A BioCASe protocl compliant service.
\item DIGIR - A DiGIR service endpoint.
\item DIGIR_MANIS - A DiGIR service slightly modified for the MANIS
network.
\item DWC_ARCHIVE - A Darwin Core Archive as defined by the Darwin Core
Text Guidelines.
\item EML - A single EML metadata document in any EML version.
\item FEED - Syndication feeds like RSS or ATOM of various flavors.
\item OAI_PMH - The Open Archives Initiative Protocol for Metadata
Harvesting.
\item OTHER - Any other service not covered by this enum so far.
\item TAPIR - A TAPIR service.
\item TCS_RDF - Taxon Concept data given as RDF based on the TDWG ontology.
\item TCS_XML - A Taxon Concept Schema document.
\item WFS - An OGC Web Feature Service.
\item WMS - An OGC Web Map Service.
}
}

\section{Multiple values passed to a parameter}{

There are some parameters you can pass multiple values to in a vector,
each value of which produces a different request (MDR: multiple different
requests). Some parameters allow multiple values to be passed in the same
request (MSR: multiple same request) in a semicolon separated string
(e.g., 'a;b'); if given we'll do a single request with that parameter
repeated for each value given (e.g., \code{foo=a&foo=b} if the parameter
is \code{foo}). Some parameters allow both MDR and MSR.

The following list shows which parameters support MDR and MSR.

\itemize{
\item basisOfRecord: MDR, MSR
\item classKey: MDR, MSR
\item country: MDR, MSR
\item familyKey: MDR, MSR
\item genusKey: MDR, MSR
\item identifiedByID: MDR, MSR
\item kingdomKey: MDR, MSR
\item license: MDR, MSR
\item locality: MDR, MSR
\item catalogNumber: MDR, MSR
\item collectionCode: MDR, MSR
\item continent: MDR, MSR
\item datasetKey: MDR, MSR
\item establishmentMeans: MDR, MSR
\item geometry: MDR, MSR
\item institutionCode: MDR, MSR
\item mediaType: MDR, MSR
\item orderKey: MDR, MSR
\item organismId: MDR, MSR
\item phylumKey: MDR, MSR
\item protocol: MDR, MSR
\item publishingCountry: MDR, MSR
\item publishingOrg: MDR, MSR
\item recordedBy: MDR, MSR
\item recordedByID: MDR, MSR
\item recordNumber: MDR, MSR
\item scientificName: MDR, MSR
\item stateProvince: MDR, MSR
\item subgenusKey: MDR, MSR
\item taxonKey: MDR, MSR
\item typeStatus: MDR, MSR
\item waterBody: MDR, MSR
\item depth: MDR
\item limit: MDR
\item q: MDR
\item year: MDR
\item repatriated: MDR
\item lastInterpreted: MDR
\item decimalLatitude: MDR
\item decimalLongitude: MDR
}

Note that you can not pass a vector > length 1 to more than 1 of the above
MDR parameters at the same time.

see also \code{\link{many-values}}
}

\section{Range queries}{

A range query is as it sounds - you query on a range of values defined by
a lower and upper limit. Do a range query by specifying the lower and upper
limit in a string like \code{depth='50,100'}. It would be more R like to
specify the range in a vector like \code{c(50,100)}, but that sort of syntax
allows you to do many searches, one for each element in the vector -
thus range queries have to differ. The following parameters support
range queries.

\itemize{
\item decimalLatitude
\item decimalLongitude
\item depth
\item elevation
\item eventDate
\item lastInterpreted
\item month
\item year
}

See also above section: semicolon and comma separated strings lead to
different outcomes for some parameters.
}

\section{Hierarchies}{

Hierarchies are returned wih each occurrence object. There is no
option no to return them from the API. However, within the \code{occ_search}
function you can select whether to return just hierarchies, just data, all
of data and hiearchies and metadata, or just metadata. If all hierarchies
are the same we just return one for you.
}

\section{curl debugging}{

You can pass parameters not defined in this function into the call to
the GBIF API to control things about the call itself using \code{curlopts}.
See an example below that passes in the \code{verbose} function to get
details on the http call.
}

\section{Scientific names vs. taxon keys}{

In the previous GBIF API and the version of rgbif that
wrapped that API, you could search the equivalent of this function with a
species name, which was convenient. However, names are messy right. So it
sorta makes sense to sort out the species key numbers you want exactly,
and then get your occurrence data with this function. GBIF has added a
parameter scientificName to allow searches by scientific names in this
function - which includes synonym taxa. \emph{Note:} that if you do use the
scientificName parameter, we will check internally that it's not a
synonym of an accepted name, and if it is, we'll search on the
accepted name. If you want to force searching by a synonym do so by
finding the GBIF identifier first with any \code{name_*} functions,
then pass that ID to the \code{taxonKey} parameter.
}

\section{WKT}{

Examples of valid WKT objects:
\itemize{
\item 'POLYGON((-19.5 34.1, 27.8 34.1, 35.9 68.1, -25.3 68.1, -19.5 34.1))'
\item 'MULTIPOLYGON(((-123 38,-116 38,-116 43,-123 43,-123 38)),((-97 41,-93 41,-93 45,-97 45,-97 41)))'
\item 'POINT(-120 40)'
\item 'LINESTRING(3 4,10 50,20 25)'
\item 'LINEARRING' ???' - Not sure how to specify this. Anyone?
}

Note that GBIF expects counter-clockwise winding order for WKT. You can
supply clockwise WKT, but GBIF treats it as an exclusion, so you get all
data not inside the WKT area. \code{\link[=occ_download]{occ_download()}} behaves differently
in that you should simply get no data back at all with clockwise WKT.
}

\section{Long WKT}{

Options for handling long WKT strings:
Note that long WKT strings are specially handled when using \code{\link{occ_search}} or
\code{\link{occ_data}}. Here are the three options for long WKT strings (> 1500 characters),
set one of these three via the parameter \code{geom_big}:
\itemize{
\item asis - the default setting. This means we don't do anything internally. That is,
we just pass on your WKT string just as we've done before in this package.
\item axe - this option uses the \pkg{sf} package to chop up your WKT string in
to many polygons, which then leads to a separate data request for each polygon piece,
then we combine all dat back together to give to you. Note that if your WKT string
is not of type polygon, we drop back to \code{asis}as there's no way to chop up
linestrings, etc. This option will in most cases be slower than the other two options.
However, this polygon splitting approach won't have the problem of
the disconnect between how many records you want and what you actually get back as
with the bbox option.

This method uses \code{sf::st_make_grid} and \code{sf::st_intersection}, which has
two parameters \code{cellsize} and \code{n}. You can tweak those parameters here by
tweaking \code{geom_size} and \code{geom_n}. \code{geom_size} seems to be more useful in
toggling the number of WKT strings you get back.

See \code{\link{wkt_parse}} to manually break make WKT bounding box from a larger WKT
string, or break a larger WKT string into many smaller ones.

\item bbox - this option checks whether your WKT string is longer than 1500 characters,
and if it is we create a bounding box from the WKT, do the GBIF search with that
bounding box, then prune the resulting data to only those occurrences in your original
WKT string. There is a big caveat however. Because we create a bounding box from the WKT,
and the \code{limit} parameter determines some subset of records to get, then when we
prune the resulting data to the WKT, the number of records you get could be less than
what you set with your \code{limit} parameter. However, you could set the limit to be
high enough so that you get all records back found in that bounding box, then you'll
get all the records available within the WKT.
}
}

\section{issue parameter}{

The options for the issue parameter (from
https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html):
\itemize{
\item BASIS_OF_RECORD_INVALID The given basis of record is impossible to interpret or seriously
different from the recommended vocabulary.
\item CONTINENT_COUNTRY_MISMATCH The interpreted continent and country do not match up.
\item CONTINENT_DERIVED_FROM_COORDINATES The interpreted continent is based on the coordinates,
not the verbatim string information.
\item CONTINENT_INVALID Uninterpretable continent values found.
\item COORDINATE_INVALID Coordinate value given in some form but GBIF is unable to interpret it.
\item COORDINATE_OUT_OF_RANGE Coordinate has invalid lat/lon values out of their decimal max
range.
\item COORDINATE_REPROJECTED The original coordinate was successfully reprojected from a
different geodetic datum to WGS84.
\item COORDINATE_REPROJECTION_FAILED The given decimal latitude and longitude could not be
reprojected to WGS84 based on the provided datum.
\item COORDINATE_REPROJECTION_SUSPICIOUS Indicates successful coordinate reprojection according
to provided datum, but which results in a datum shift larger than 0.1 decimal degrees.
\item COORDINATE_ROUNDED Original coordinate modified by rounding to 5 decimals.
\item COUNTRY_COORDINATE_MISMATCH The interpreted occurrence coordinates fall outside of the
indicated country.
\item COUNTRY_DERIVED_FROM_COORDINATES The interpreted country is based on the coordinates, not
the verbatim string information.
\item COUNTRY_INVALID Uninterpretable country values found.
\item COUNTRY_MISMATCH Interpreted country for dwc:country and dwc:countryCode contradict each
other.
\item DEPTH_MIN_MAX_SWAPPED Set if supplied min>max
\item DEPTH_NON_NUMERIC Set if depth is a non numeric value
\item DEPTH_NOT_METRIC Set if supplied depth is not given in the metric system, for example
using feet instead of meters
\item DEPTH_UNLIKELY Set if depth is larger than 11.000m or negative.
\item ELEVATION_MIN_MAX_SWAPPED Set if supplied min > max elevation
\item ELEVATION_NON_NUMERIC Set if elevation is a non numeric value
\item ELEVATION_NOT_METRIC Set if supplied elevation is not given in the metric system, for
example using feet instead of meters
\item ELEVATION_UNLIKELY Set if elevation is above the troposphere (17km) or below 11km
(Mariana Trench).
\item GEODETIC_DATUM_ASSUMED_WGS84 Indicating that the interpreted coordinates assume they are
based on WGS84 datum as the datum was either not indicated or interpretable.
\item GEODETIC_DATUM_INVALID The geodetic datum given could not be interpreted.
\item IDENTIFIED_DATE_INVALID The date given for dwc:dateIdentified is invalid and cant be
interpreted at all.
\item IDENTIFIED_DATE_UNLIKELY The date given for dwc:dateIdentified is in the future or before
Linnean times (1700).
\item MODIFIED_DATE_INVALID A (partial) invalid date is given for dc:modified, such as a non
existing date, invalid zero month, etc.
\item MODIFIED_DATE_UNLIKELY The date given for dc:modified is in the future or predates unix
time (1970).
\item MULTIMEDIA_DATE_INVALID An invalid date is given for dc:created of a multimedia object.
\item MULTIMEDIA_URI_INVALID An invalid uri is given for a multimedia object.
\item PRESUMED_NEGATED_LATITUDE Latitude appears to be negated, e.g. 32.3 instead of -32.3
\item PRESUMED_NEGATED_LONGITUDE Longitude appears to be negated, e.g. 32.3 instead of -32.3
\item PRESUMED_SWAPPED_COORDINATE Latitude and longitude appear to be swapped.
\item RECORDED_DATE_INVALID A (partial) invalid date is given, such as a non existing date,
invalid zero month, etc.
\item RECORDED_DATE_MISMATCH The recording date specified as the eventDate string and the
individual year, month, day are contradicting.
\item RECORDED_DATE_UNLIKELY The recording date is highly unlikely, falling either into the
future or represents a very old date before 1600 that predates modern taxonomy.
\item REFERENCES_URI_INVALID An invalid uri is given for dc:references.
\item TAXON_MATCH_FUZZY Matching to the taxonomic backbone can only be done using a fuzzy, non
exact match.
\item TAXON_MATCH_HIGHERRANK Matching to the taxonomic backbone can only be done on a higher
rank and not the scientific name.
\item TAXON_MATCH_NONE Matching to the taxonomic backbone cannot be done cause there was no
match at all or several matches with too little information to keep them apart (homonyms).
\item TYPE_STATUS_INVALID The given type status is impossible to interpret or seriously
different from the recommended vocabulary.
\item ZERO_COORDINATE Coordinate is the exact 0/0 coordinate, often indicating a bad null
coordinate.
}
}

\section{Counts}{

There is a slight difference in the way records are counted here vs.
results from \code{\link{occ_count}}. For equivalent outcomes, in this
function use \code{hasCoordinate=TRUE}, and \code{hasGeospatialIssue=FALSE}
to have the same outcome using \code{\link{occ_count}} with
\code{isGeoreferenced=TRUE}
}

\section{occ_data vs. occ_search}{

This does nearly the same thing as \code{\link[=occ_search]{occ_search()}}, but
is simplified for speed, and is for the most common use case where
user just wants occurrence data, and not other information like taxon
hierarchies and media (e.g., images). Alot of time in \code{\link[=occ_search]{occ_search()}}
is used parsing data to be more useable downstream. We do less of that
in this function.

There are a number of data fields GBIF returns that we drop to speed up
processing time within R. These fields take extra time to process
because they are deeply nested and so take extra time to check if
they are empty or not, and if not, figure out how to parse them
into a data.frame. The fields are:
\itemize{
\item \code{gadm}
\item \code{media}
\item \code{facts}
\item \code{relations}
\item \code{extensions}
\item \code{identifiers}
\item \code{recordedByIDs}
\item \code{identifiedByIDs}
}

To get these fields use \code{\link[=occ_search]{occ_search()}} instead.
}

\examples{
\dontrun{
(key <- name_backbone(name='Encelia californica')$speciesKey)
occ_data(taxonKey = key, limit = 4)
(res <- occ_data(taxonKey = key, limit = 400))

# Return 20 results, this is the default by the way
(key <- name_suggest(q='Helianthus annuus', rank='species')$data$key[1])
occ_data(taxonKey=key, limit=20)

# Instead of getting a taxon key first, you can search for a name directly
## However, note that using this approach (with \code{scientificName="..."})
## you are getting synonyms too. The results for using \code{scientifcName}
## and \code{taxonKey} parameters are the same in this case, but I wouldn't
## be surprised if for some names they return different results
occ_data(scientificName = 'Ursus americanus', curlopts=list(verbose=TRUE))
key <- name_backbone(name = 'Ursus americanus', rank='species')$usageKey
occ_data(taxonKey = key)

# Search by dataset key
occ_data(datasetKey='7b5d6a48-f762-11e1-a439-00145eb45e9a', limit=10)

# Search by catalog number
occ_data(catalogNumber="49366", limit=10)
## separate requests: use a vector of strings
occ_data(catalogNumber=c("49366","Bird.27847588"), limit=10)
## one request, many instances of same parameter: use semi-colon sep. string
occ_data(catalogNumber="49366;Bird.27847588", limit=10)

# Use paging parameters (limit and start) to page. Note the different results
# for the two queries below.
occ_data(datasetKey='7b5d6a48-f762-11e1-a439-00145eb45e9a',start=10,limit=5)
occ_data(datasetKey='7b5d6a48-f762-11e1-a439-00145eb45e9a',start=20,limit=5)

# Many dataset keys
## separate requests: use a vector of strings
occ_data(datasetKey=c("50c9509d-22c7-4a22-a47d-8c48425ef4a7",
   "7b5d6a48-f762-11e1-a439-00145eb45e9a"), limit=20)
## one request, many instances of same parameter: use semi-colon sep. string
v="50c9509d-22c7-4a22-a47d-8c48425ef4a7;7b5d6a48-f762-11e1-a439-00145eb45e9a"
occ_data(datasetKey = v, limit=20)

# Search by recorder
occ_data(recordedBy="smith", limit=20)

# Many collector names
## separate requests: use a vector of strings
occ_data(recordedBy=c("smith","BJ Stacey"), limit=10)
## one request, many instances of same parameter: use semi-colon sep. string
occ_data(recordedBy="smith;BJ Stacey", limit=10)

# recordedByID
occ_data(recordedByID="https://orcid.org/0000-0003-1691-239X", limit=20)
## many at once
### separate searches
ids <- c("https://orcid.org/0000-0003-1691-239X",
  "https://orcid.org/0000-0001-7569-1828",
  "https://orcid.org/0000-0002-0596-5376")
res <- occ_data(recordedByID=ids, limit=20)
res[[1]]$data$recordedByIDs[[1]]
res[[2]]$data$recordedByIDs[[1]]
res[[3]]$data$recordedByIDs[[1]]
### all in one search
res <- occ_data(recordedByID=paste0(ids, collapse=";"), limit=20)
unique(vapply(res$data$recordedByIDs, "[[", "", "value"))

# identifiedByID
occ_data(identifiedByID="https://orcid.org/0000-0003-4710-2648", limit=20)

# Pass in curl options for extra fun
occ_data(taxonKey=2433407, limit=20, curlopts=list(verbose=TRUE))
occ_data(taxonKey=2433407, limit=20,
  curlopts = list(
    noprogress = FALSE,
    progressfunction = function(down, up) {
      cat(sprintf("up: \%d | down \%d\n", up, down))
      return(TRUE)
    }
  )
)
# occ_data(taxonKey=2433407, limit=20, curlopts=list(timeout_ms=1))

# Search for many species
splist <- c('Cyanocitta stelleri', 'Junco hyemalis', 'Aix sponsa')
keys <- sapply(splist, function(x) name_suggest(x)$data$key[1], USE.NAMES=FALSE)
## separate requests: use a vector of strings
occ_data(taxonKey = keys, limit=5)
## one request, many instances of same parameter: use semi-colon sep. string
occ_data(taxonKey = paste0(keys, collapse = ";"), limit=5)

# Search using a synonym name
#  Note that you'll see a message printing out that the accepted name will
# be used
occ_data(scientificName = 'Pulsatilla patens', limit=5)

# Search on latitidue and longitude
occ_data(decimalLatitude=40, decimalLongitude=-120, limit = 10)

# Search on a bounding box
## in well known text format
### polygon
occ_data(geometry='POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))',
  limit=20)
### multipolygon
wkt <- 'MULTIPOLYGON(((-123 38,-116 38,-116 43,-123 43,-123 38)),
   ((-97 41,-93 41,-93 45,-97 45,-97 41)))'
occ_data(geometry = gsub("\n\\\\s+", "", wkt), limit = 20)
### polygon and taxonkey
key <- name_suggest(q='Aesculus hippocastanum')$data$key[1]
occ_data(taxonKey=key,
 geometry='POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))',
 limit=20)
## or using bounding box, converted to WKT internally
occ_data(geometry=c(-125.0,38.4,-121.8,40.9), limit=20)

## you can seaerch on many geometry objects
### separate requests: use a vector of strings
wkts <-
c('POLYGON((-102.2 46,-102.2 43.7,-93.9 43.7,-93.9 46,-102.2 46))',
'POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))')
occ_data(geometry = wkts, limit=20)
### one request, many instances of same parameter: use semi-colon sep. string
occ_data(geometry = paste0(wkts, collapse = ";"), limit=20)


# Search on a long WKT string - too long for a GBIF search API request
## By default, a very long WKT string will likely cause a request failure as
## GBIF only handles strings up to about 1500 characters long. You can leave as is, or
##  - Alternatively, you can choose to break up your polygon into many, and do a
##      data request on each piece, and the output is put back together (see below)
##  - Or, 2nd alternatively, you could use the GBIF download API
wkt <- "POLYGON((-9.178796777343678 53.22769021556159,
-12.167078027343678 51.56540789297837,
-12.958093652343678 49.78333685689162,-11.024499902343678 49.21251756301334,
-12.079187402343678 46.68179685941719,-15.067468652343678 45.83103608186854,
-15.770593652343678 43.58271629699817,-15.067468652343678 41.57676278827219,
-11.815515527343678 40.44938999172728,-12.958093652343678 37.72112962230871,
-11.639734277343678 36.52987439429357,-8.299890527343678 34.96062625095747,
-8.739343652343678 32.62357394385735,-5.223718652343678 30.90497915232165,
1.1044063476563224 31.80562077746643,1.1044063476563224 30.754036557416256,
6.905187597656322 32.02942785462211,5.147375097656322 32.99292810780193,
9.629796972656322 34.164474406524725,10.860265722656322 32.91918014319603,
14.551671972656322 33.72700959356651,13.409093847656322 34.888564192275204,
16.748937597656322 35.104560368110114,19.561437597656322 34.81643887792552,
18.594640722656322 36.38849705969625,22.989171972656322 37.162874858929854,
19.825109472656322 39.50651757842751,13.760656347656322 38.89353140585116,
14.112218847656322 42.36091601976124,10.596593847656322 41.11488736647705,
9.366125097656322 43.70991402658437,5.059484472656322 42.62015372417812,
2.3348750976563224 45.21526500321446,-0.7412967773436776 46.80225692528942,
6.114171972656322 47.102229890207894,8.047765722656322 45.52399303437107,
12.881750097656322 48.22681126957933,9.190343847656322 48.693079457106684,
8.750890722656322 50.68283120621287,5.059484472656322 50.40356146487845,
4.268468847656322 52.377558897655156,1.4559688476563224 53.28027243658647,
0.8407344726563224 51.62000971578333,0.5770625976563224 49.32721423860726,
-2.5869999023436776 49.49875947592088,-2.4991092773436776 51.18135535408638,
-2.0596561523436776 52.53822562473851,-4.696374902343678 51.67454591918756,
-5.311609277343678 50.009802108095776,-6.629968652343678 48.75106196817059,
-7.684656152343678 50.12263634382465,-6.190515527343678 51.83776110910459,
-5.047937402343678 54.267098895684235,-6.893640527343678 53.69860705549198,
-8.915124902343678 54.77719740243195,-12.079187402343678 54.52294465763567,
-13.573328027343678 53.437631551347174,
-11.288171777343678 53.48995552517918,
-9.178796777343678 53.22769021556159))"
wkt <- gsub("\n", " ", wkt)

#### Default option with large WKT string fails
# res <- occ_data(geometry = wkt)

#### if WKT too long, with 'geom_big=bbox': makes into bounding box
if (interactive()){
res <- occ_data(geometry = wkt, geom_big = "bbox")
library("rgeos")
library("sp")
wktsp <- readWKT(wkt)
plot(wktsp)
z <- data.frame(res$data)
coordinates(z) <- ~decimalLongitude+decimalLatitude
points(z)
}

#### Or, use 'geom_big=axe'
(res <- occ_data(geometry = wkt, geom_big = "axe"))
##### manipulate essentially number of polygons that result, so number of requests
###### default geom_size is 40
###### fewer calls
(res <- occ_data(geometry = wkt, geom_big = "axe", geom_size=50))
###### more calls
(res <- occ_data(geometry = wkt, geom_big = "axe", geom_size=30))

# Search on country
occ_data(country='US', limit=20)
isocodes[grep("France", isocodes$name),"code"]
occ_data(country='FR', limit=20)
occ_data(country='DE', limit=20)
### separate requests: use a vector of strings
occ_data(country=c('US','DE'), limit=20)
### one request, many instances of same parameter: use semi-colon sep. string
occ_data(country = 'US;DE', limit=20)

# Get only occurrences with lat/long data
occ_data(taxonKey=key, hasCoordinate=TRUE, limit=20)

# Get only occurrences that were recorded as living specimens
occ_data(basisOfRecord="LIVING_SPECIMEN", hasCoordinate=TRUE, limit=20)
## multiple values in a vector = a separate request for each value
occ_data(taxonKey=key,
  basisOfRecord=c("OBSERVATION", "HUMAN_OBSERVATION"), limit=20)
## mutiple values in a single string, ";" separated = one request including all values
occ_data(taxonKey=key,
  basisOfRecord="OBSERVATION;HUMAN_OBSERVATION", limit=20)

# Get occurrences for a particular eventDate
occ_data(taxonKey=key, eventDate="2013", limit=20)
occ_data(taxonKey=key, year="2013", limit=20)
occ_data(taxonKey=key, month="6", limit=20)

# Get occurrences based on depth
key <- name_backbone(name='Salmo salar', kingdom='animals')$speciesKey
occ_data(taxonKey=key, depth=1, limit=20)

# Get occurrences based on elevation
key <- name_backbone(name='Puma concolor', kingdom='animals')$speciesKey
occ_data(taxonKey=key, elevation=50, hasCoordinate=TRUE, limit=20)

# Get occurrences based on institutionCode
occ_data(institutionCode="TLMF", limit=20)
### separate requests: use a vector of strings
occ_data(institutionCode=c("TLMF","ArtDatabanken"), limit=20)
### one request, many instances of same parameter: use semi-colon sep. string
occ_data(institutionCode = "TLMF;ArtDatabanken", limit=20)

# Get occurrences based on collectionCode
occ_data(collectionCode="Floristic Databases MV - Higher Plants", limit=20)
### separate requests: use a vector of strings
occ_data(collectionCode=c("Floristic Databases MV - Higher Plants",
  "Artport"), limit = 20)
### one request, many instances of same parameter: use semi-colon sep. string
occ_data(collectionCode = "Floristic Databases MV - Higher Plants;Artport",
  limit = 20)

# Get only those occurrences with spatial issues
occ_data(taxonKey=key, hasGeospatialIssue=TRUE, limit=20)

# Search using a query string
occ_data(search="kingfisher", limit=20)

# search on repatriated - doesn't work right now
# occ_data(repatriated = "")

# search on phylumKey
occ_data(phylumKey = 7707728, limit = 5)

# search on kingdomKey
occ_data(kingdomKey = 1, limit = 5)

# search on classKey
occ_data(classKey = 216, limit = 5)

# search on orderKey
occ_data(orderKey = 7192402, limit = 5)

# search on familyKey
occ_data(familyKey = 3925, limit = 5)

# search on genusKey
occ_data(genusKey = 1935496, limit = 5)

# search on establishmentMeans
occ_data(establishmentMeans = "INVASIVE", limit = 5)
occ_data(establishmentMeans = "NATIVE", limit = 5)
occ_data(establishmentMeans = "UNCERTAIN", limit = 5)
### separate requests: use a vector of strings
occ_data(establishmentMeans = c("INVASIVE", "NATIVE"), limit = 5)
### one request, many instances of same parameter: use semi-colon sep. string
occ_data(establishmentMeans = "INVASIVE;NATIVE", limit = 5)

# search on protocol
occ_data(protocol = "DIGIR", limit = 5)

# search on license
occ_data(license = "CC_BY_4_0", limit = 5)

# search on organismId
occ_data(organismId = "100", limit = 5)

# search on publishingOrg
occ_data(publishingOrg = "28eb1a3f-1c15-4a95-931a-4af90ecb574d", limit = 5)

# search on stateProvince
occ_data(stateProvince = "California", limit = 5)

# search on waterBody
occ_data(waterBody = "pacific ocean", limit = 5)

# search on locality
occ_data(locality = "Trondheim", limit = 5)
### separate requests: use a vector of strings
res <- occ_data(locality = c("Trondheim", "Hovekilen"), limit = 5)
res$Trondheim$data
res$Hovekilen$data
### one request, many instances of same parameter: use semi-colon sep. string
occ_data(locality = "Trondheim;Hovekilen", limit = 5)


# Range queries
## See Detail for parameters that support range queries
occ_data(depth='50,100', limit = 20)
### this is not a range search, but does two searches for each depth
occ_data(depth=c(50,100), limit = 20)

## Range search with year
occ_data(year='1999,2000', limit=20)

## Range search with latitude
occ_data(decimalLatitude='29.59,29.6', limit = 20)

# Search by specimen type status
## Look for possible values of the typeStatus parameter looking at the typestatus dataset
occ_data(typeStatus = 'allotype', limit = 20)$data[,c('name','typeStatus')]

# Search by specimen record number
## This is the record number of the person/group that submitted the data, not GBIF's numbers
## You can see that many different groups have record number 1, so not super helpful
occ_data(recordNumber = 1, limit = 20)$data[,c('name','recordNumber','recordedBy')]

# Search by last time interpreted: Date the record was last modified in GBIF
## The lastInterpreted parameter accepts ISO 8601 format dates, including
## yyyy, yyyy-MM, yyyy-MM-dd, or MM-dd. Range queries are accepted for lastInterpreted
occ_data(lastInterpreted = '2016-04-02', limit = 20)

# Search for occurrences with images
occ_data(mediaType = 'StillImage', limit = 20)
occ_data(mediaType = 'MovingImage', limit = 20)
occ_data(mediaType = 'Sound', limit = 20)

# Search by continent
## One of africa, antarctica, asia, europe, north_america, oceania, or
## south_america
occ_data(continent = 'south_america', limit = 20)$meta
occ_data(continent = 'africa', limit = 20)$meta
occ_data(continent = 'oceania', limit = 20)$meta
occ_data(continent = 'antarctica', limit = 20)$meta
### separate requests: use a vector of strings
occ_data(continent = c('south_america', 'oceania'), limit = 20)
### one request, many instances of same parameter: use semi-colon sep. string
occ_data(continent = 'south_america;oceania', limit = 20)

# Query based on issues - see Details for options
## one issue
x <- occ_data(taxonKey=1, issue='DEPTH_UNLIKELY', limit = 20)
x$data[,c('name','key','decimalLatitude','decimalLongitude','depth')]
## two issues
occ_data(taxonKey=1, issue=c('DEPTH_UNLIKELY','COORDINATE_ROUNDED'), limit = 20)
# Show all records in the Arizona State Lichen Collection that cant be matched to the GBIF
# backbone properly:
occ_data(datasetKey='84c0e1a0-f762-11e1-a439-00145eb45e9a',
   issue=c('TAXON_MATCH_NONE','TAXON_MATCH_HIGHERRANK'), limit = 20)

# Parsing output by issue
(res <- occ_data(geometry='POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))', limit = 50))
## what do issues mean, can print whole table, or search for matches
head(gbif_issues())
gbif_issues()[ gbif_issues()$code \%in\% c('cdround','cudc','gass84','txmathi'), ]
## or parse issues in various ways
### remove data rows with certain issue classes
library('magrittr')
res \%>\% occ_issues(gass84)
### split issues into separate columns
res \%>\% occ_issues(mutate = "split")
### expand issues to more descriptive names
res \%>\% occ_issues(mutate = "expand")
### split and expand
res \%>\% occ_issues(mutate = "split_expand")
### split, expand, and remove an issue class
res \%>\% occ_issues(-cudc, mutate = "split_expand")
}
}
\references{
https://www.gbif.org/developer/occurrence#search
}
\seealso{
\code{\link[=downloads]{downloads()}}, \code{\link[=occ_search]{occ_search()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_datasets.R
\name{occ_download_datasets}
\alias{occ_download_datasets}
\title{List datasets for a download}
\usage{
occ_download_datasets(key, limit = 20, start = 0, curlopts = list())
}
\arguments{
\item{key}{A key generated from a request, like that from \code{\link[=occ_download]{occ_download()}}}

\item{limit}{(integer/numeric) Number of records to return. Default: 20,
Max: 1000}

\item{start}{(integer/numeric) Record number to start at. Default: 0}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
a list with two slots:
\itemize{
\item meta: a single row data.frame with columns: \code{offset}, \code{limit},
\code{endofrecords}, \code{count}
\item results: a tibble with the results, of three columns: \code{downloadKey},
\code{datasetKey}, \code{numberRecords}
}
}
\description{
List datasets for a download
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\examples{
\dontrun{
occ_download_datasets(key="0003983-140910143529206")
occ_download_datasets(key="0003983-140910143529206", limit = 3)
occ_download_datasets(key="0003983-140910143529206", limit = 3, start = 10)
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{stylegeojson}
\alias{stylegeojson}
\title{Style a data.frame prior to converting to geojson.}
\usage{
stylegeojson(...)
}
\description{
This function is defunct.  See the package spocc for similar functionality.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{occurrencelist_many}
\alias{occurrencelist_many}
\title{occurrencelist_many is the same as occurrencelist, but takes in a vector
of species names.}
\usage{
occurrencelist_many(...)
}
\description{
This function is defunct.
}
\seealso{
occ_search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oai-pmh.R
\name{gbif_oai}
\alias{gbif_oai}
\alias{gbif_oai_identify}
\alias{gbif_oai_list_identifiers}
\alias{gbif_oai_list_records}
\alias{gbif_oai_list_metadataformats}
\alias{gbif_oai_list_sets}
\alias{gbif_oai_get_records}
\title{GBIF registry data via OAI-PMH}
\usage{
gbif_oai_identify(...)

gbif_oai_list_identifiers(
  prefix = "oai_dc",
  from = NULL,
  until = NULL,
  set = NULL,
  token = NULL,
  as = "df",
  ...
)

gbif_oai_list_records(
  prefix = "oai_dc",
  from = NULL,
  until = NULL,
  set = NULL,
  token = NULL,
  as = "df",
  ...
)

gbif_oai_list_metadataformats(id = NULL, ...)

gbif_oai_list_sets(token = NULL, as = "df", ...)

gbif_oai_get_records(ids, prefix = "oai_dc", as = "parsed", ...)
}
\arguments{
\item{...}{Curl options passed on to \code{httr::GET}}

\item{prefix}{(character) A string to specify the metadata format in OAI-PMH
requests issued to the repository. The default (\code{"oai_dc"}) corresponds
to the mandatory OAI unqualified Dublin Core metadata schema.}

\item{from}{(character) string giving datestamp to be used as lower bound
for datestamp-based selective harvesting (i.e., only harvest records with
datestamps in the given range). Dates and times must be encoded using ISO
8601. The trailing Z must be used when including time. OAI-PMH implies
UTC for data/time specifications.}

\item{until}{(character) Datestamp to be used as an upper bound, for
datestamp-based selective harvesting (i.e., only harvest records with
datestamps in the given range).}

\item{set}{(character) A set to be used for selective harvesting (i.e., only
harvest records in the given set).}

\item{token}{(character) a token previously provided by the server to resume
a request where it last left off. 50 is max number of records returned.
We will loop for you internally to get all the records you asked for.}

\item{as}{(character) What to return. One of "df" (for data.frame;
default), "list" (get a list), or "raw" (raw text). For
\code{gbif_oai_get_records}, one of "parsed" or "raw"}

\item{id, ids}{(character) The OAI-PMH identifier for the record. Optional.}
}
\value{
raw text, list or data.frame, depending on requested output via
\code{as} parameter
}
\description{
GBIF registry data via OAI-PMH
}
\details{
These functions only work with GBIF registry data, and do so
via the OAI-PMH protocol
(https://www.openarchives.org/OAI/openarchivesprotocol.html)
}
\examples{
\dontrun{
gbif_oai_identify()

today <- format(Sys.Date(), "\%Y-\%m-\%d")
gbif_oai_list_identifiers(from = today)
gbif_oai_list_identifiers(set = "country:NL")

gbif_oai_list_records(from = today)
gbif_oai_list_records(set = "country:NL")

gbif_oai_list_metadataformats()
gbif_oai_list_metadataformats(id = "9c4e36c1-d3f9-49ce-8ec1-8c434fa9e6eb")

gbif_oai_list_sets()
gbif_oai_list_sets(as = "list")

gbif_oai_get_records("9c4e36c1-d3f9-49ce-8ec1-8c434fa9e6eb")
ids <- c("9c4e36c1-d3f9-49ce-8ec1-8c434fa9e6eb",
         "e0f1bb8a-2d81-4b2a-9194-d92848d3b82e")
gbif_oai_get_records(ids)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elevation.r
\name{elevation}
\alias{elevation}
\title{Get elevation for lat/long points from a data.frame or list of points.}
\usage{
elevation(
  input = NULL,
  latitude = NULL,
  longitude = NULL,
  latlong = NULL,
  elevation_model = "srtm3",
  username = Sys.getenv("GEONAMES_USER"),
  key,
  curlopts,
  ...
)
}
\arguments{
\item{input}{A data.frame of lat/long data. There must be columns
decimalLatitude and decimalLongitude.}

\item{latitude}{A vector of latitude's. Must be the same length as the
longitude vector.}

\item{longitude}{A vector of longitude's. Must be the same length as
the latitude vector.}

\item{latlong}{A vector of lat/long pairs. See examples.}

\item{elevation_model}{(character) one of srtm3 (default), srtm1, astergdem,
or gtopo30. See "Elevation models" below for more}

\item{username}{(character) Required. An GeoNames user name. See Details.}

\item{key, curlopts}{defunct. see docs}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}
see \code{curl::curl_options()} for curl options}
}
\value{
A new column named \code{elevation_geonames} in the supplied data.frame
or a vector with elevation of each location in meters. Note that data from
GBIF can already have a column named \code{elevation}, thus the column we
add is named differently.
}
\description{
Uses the GeoNames web service
}
\section{GeoNames user name}{

To get a GeoNames user name, register for an account at
http://www.geonames.org/login - then you can enable your account for the
GeoNames webservice on your account page
(http://www.geonames.org/manageaccount). Once you are enabled to use
the webservice, you can pass in your username to the \code{username}
parameter. Better yet, store your username in your \code{.Renviron} file, or
similar (e.g., .zshrc or .bash_profile files) and read it in via
\code{Sys.getenv()} as in the examples below. By default we do
\code{Sys.getenv("GEONAMES_USER")} for the \code{username} parameter.
}

\section{Elevation models}{

\itemize{
\item srtm3:
\itemize{
\item sample area: ca 90m x 90m
\item result: a single number giving the elevation in meters according to
srtm3, ocean areas have been masked as "no data" and have been assigned
a value of -32768
}
\item srtm1:
\itemize{
\item sample area: ca 30m x 30m
\item result: a single number giving the elevation in meters according to
srtm1, ocean areas have been masked as "no data" and have been assigned
a value of -32768
}
\item astergdem (Aster Global Digital Elevation Model V2 2011):
\itemize{
\item sample area: ca 30m x 30m, between 83N and 65S latitude
\item result: a single number giving the elevation in meters according to
aster gdem, ocean areas have been masked as "no data" and have been
assigned a value of -32768
}
\item gtopo30:
\itemize{
\item sample area: ca 1km x 1km
\item result: a single number giving the elevation in meters according to
gtopo30, ocean areas have been masked as "no data" and have been
assigned a value of -9999
}
}
}

\examples{
\dontrun{
user <- Sys.getenv("GEONAMES_USER")

occ_key <- name_suggest('Puma concolor')$key[1]
dat <- occ_search(taxonKey = occ_key, limit = 300, hasCoordinate = TRUE)
head( elevation(dat$data, username = user) )

# Pass in a vector of lat's and a vector of long's
elevation(latitude = dat$data$decimalLatitude[1:10],
  longitude = dat$data$decimalLongitude[1:10],
  username = user, verbose = TRUE)

# Pass in lat/long pairs in a single vector
pairs <- list(c(31.8496,-110.576060), c(29.15503,-103.59828))
elevation(latlong=pairs, username = user)

# Pass on curl options
pairs <- list(c(31.8496,-110.576060), c(29.15503,-103.59828))
elevation(latlong=pairs, username = user, verbose = TRUE)

# different elevation models
lats <- dat$data$decimalLatitude[1:5]
lons <- dat$data$decimalLongitude[1:5]
elevation(latitude = lats, longitude = lons, elevation_model = "srtm3")
elevation(latitude = lats, longitude = lons, elevation_model = "srtm1")
elevation(latitude = lats, longitude = lons, elevation_model = "astergdem")
elevation(latitude = lats, longitude = lons, elevation_model = "gtopo30")
}
}
\references{
GeoNames http://www.geonames.org/export/web-services.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gbifmap_dens}
\alias{gbifmap_dens}
\title{Make a simple map to visualize GBIF data density data}
\usage{
gbifmap_dens(...)
}
\description{
This function is defunct.
}
\seealso{
gbifmap
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgbif-package.r
\docType{data}
\name{typestatus}
\alias{typestatus}
\title{Type status options for GBIF searching}
\description{
\itemize{
\item name. Name of type.
\item description. Description of the type.
}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_get.r
\name{occ_download_get}
\alias{occ_download_get}
\title{Get a download from GBIF.}
\usage{
occ_download_get(key, path = ".", overwrite = FALSE, ...)
}
\arguments{
\item{key}{A key generated from a request, like that from \code{occ_download}}

\item{path}{Path to write zip file to. Default: \code{"."}, with a
\code{.zip} appended to the end.}

\item{overwrite}{Will only overwrite existing path if TRUE.}

\item{...}{named curl options passed on to
\link[crul:verb-GET]{crul::verb-GET}. see \code{curl::curl_options()} for curl options}
}
\description{
Get a download from GBIF.
}
\details{
Downloads the zip file to a directory you specify on your machine.
\code{\link[crul:HttpClient]{crul::HttpClient()}} is used internally to write the zip file to
disk. See \link[crul:writing-options]{crul::writing-options}. This function only downloads the file.
See \code{occ_download_import} to open a downloaded file in your R session.
The speed of this function is of course proportional to the size of the
file to download. For example, a 58 MB file on my machine took about
26 seconds.
}
\note{
see \link{downloads} for an overview of GBIF downloads methods

This function used to check for HTTP response content type, but
it has changed enough that we no longer check it. If you run into issues
with this function, open an issue in the GitHub repository.
}
\examples{
\dontrun{
occ_download_get("0000066-140928181241064")
occ_download_get("0003983-140910143529206", overwrite = TRUE)
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbif_names.R
\name{gbif_names}
\alias{gbif_names}
\title{View highlighted terms in name results from GBIF.}
\usage{
gbif_names(input, output = NULL, browse = TRUE)
}
\arguments{
\item{input}{Input output from occ_search}

\item{output}{Output folder path. If not given uses temporary folder.}

\item{browse}{(logical) Browse output (default: \code{TRUE})}
}
\description{
View highlighted terms in name results from GBIF.
}
\examples{
\dontrun{
# browse=FALSE returns path to file
gbif_names(name_lookup(query='snake', hl=TRUE), browse=FALSE)

(out <- name_lookup(query='canada', hl=TRUE, limit=5))
gbif_names(out)
gbif_names(name_lookup(query='snake', hl=TRUE))
gbif_names(name_lookup(query='bird', hl=TRUE))

# or not highlight
gbif_names(name_lookup(query='bird', limit=200))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{occurrencelist_all}
\alias{occurrencelist_all}
\title{Occurrencelist_all carries out an occurrencelist query for a single name and
all its name variants according to GBIF's name matching.}
\usage{
occurrencelist_all(...)
}
\description{
This function is defunct.
}
\seealso{
occ_search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodes.r
\name{nodes}
\alias{nodes}
\title{Nodes metadata.}
\usage{
nodes(
  data = "all",
  uuid = NULL,
  query = NULL,
  identifier = NULL,
  identifierType = NULL,
  limit = 100,
  start = NULL,
  isocode = NULL,
  curlopts = list()
)
}
\arguments{
\item{data}{The type of data to get. One or more of: 'organization',
'endpoint', 'identifier', 'tag', 'machineTag', 'comment',
'pendingEndorsement', 'country', 'dataset', 'installation', or the
special 'all'. Default: \code{'all'}}

\item{uuid}{UUID of the data node provider. This must be specified if data
is anything other than 'all'.}

\item{query}{Query nodes. Only used when \code{data='all'}}

\item{identifier}{The value for this parameter can be a simple string or
integer, e.g. \code{identifier=120}. This parameter doesn't seem to work right
now.}

\item{identifierType}{Used in combination with the identifier parameter to
filter identifiers by identifier type. See details. This parameter doesn't
seem to work right now.}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{isocode}{A 2 letter country code. Only used if data='country'.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\description{
Nodes metadata.
}
\details{
identifierType options:

\itemize{
\item {DOI} No description.
\item {FTP} No description.
\item {GBIF_NODE} Identifies the node (e.g: \code{DK} for Denmark, \code{sp2000}
for Species 2000).
\item {GBIF_PARTICIPANT} Participant identifier from the GBIF IMS
Filemaker system.
\item {GBIF_PORTAL} Indicates the identifier originated from an
auto_increment column in the portal.data_provider or portal.data_resource
table respectively.
\item {HANDLER} No description.
\item {LSID} Reference controlled by a separate system, used for example
by DOI.
\item {SOURCE_ID} No description.
\item {UNKNOWN} No description.
\item {URI} No description.
\item {URL} No description.
\item {UUID} No description.
}
}
\examples{
\dontrun{
nodes(limit=5)
nodes(uuid="1193638d-32d1-43f0-a855-8727c94299d8")
nodes(data='identifier', uuid="03e816b3-8f58-49ae-bc12-4e18b358d6d9")
nodes(data=c('identifier','organization','comment'),
  uuid="03e816b3-8f58-49ae-bc12-4e18b358d6d9")

uuids = c("8cb55387-7802-40e8-86d6-d357a583c596",
  "02c40d2a-1cba-4633-90b7-e36e5e97aba8",
  "7a17efec-0a6a-424c-b743-f715852c3c1f",
  "b797ce0f-47e6-4231-b048-6b62ca3b0f55",
  "1193638d-32d1-43f0-a855-8727c94299d8",
  "d3499f89-5bc0-4454-8cdb-60bead228a6d",
  "cdc9736d-5ff7-4ece-9959-3c744360cdb3",
  "a8b16421-d80b-4ef3-8f22-098b01a89255",
  "8df8d012-8e64-4c8a-886e-521a3bdfa623",
  "b35cf8f1-748d-467a-adca-4f9170f20a4e",
  "03e816b3-8f58-49ae-bc12-4e18b358d6d9",
  "073d1223-70b1-4433-bb21-dd70afe3053b",
  "07dfe2f9-5116-4922-9a8a-3e0912276a72",
  "086f5148-c0a8-469b-84cc-cce5342f9242",
  "0909d601-bda2-42df-9e63-a6d51847ebce",
  "0e0181bf-9c78-4676-bdc3-54765e661bb8",
  "109aea14-c252-4a85-96e2-f5f4d5d088f4",
  "169eb292-376b-4cc6-8e31-9c2c432de0ad",
  "1e789bc9-79fc-4e60-a49e-89dfc45a7188",
  "1f94b3ca-9345-4d65-afe2-4bace93aa0fe")

res <- lapply(uuids, function(x) nodes(x, data='identifier')$data)
res <- res[!sapply(res, NROW)==0]
res[1]

# Pass on curl options
nodes(limit=20, curlopts=list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#nodes}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gist}
\alias{gist}
\title{Post a file as a Github gist}
\usage{
gist(...)
}
\description{
This function is defunct.  See the package gistr for similar functionality.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgbif-package.r
\name{rgbif-defunct}
\alias{rgbif-defunct}
\title{Defunct functions in rgbif}
\description{
\itemize{
\item \code{\link[=density_spplist]{density_spplist()}}: service no longer provided
\item \code{\link[=densitylist]{densitylist()}}: service no longer provided
\item \code{\link[=gbifdata]{gbifdata()}}: service no longer provided
\item \code{\link[=gbifmap_dens]{gbifmap_dens()}}: service no longer provided
\item \code{\link[=gbifmap_list]{gbifmap_list()}}: service no longer provided
\item \code{\link[=occurrencedensity]{occurrencedensity()}}: service no longer provided
\item \code{\link[=providers]{providers()}}: service no longer provided
\item \code{\link[=resources]{resources()}}: service no longer provided
\item \code{\link[=taxoncount]{taxoncount()}}: service no longer provided
\item \code{\link[=taxonget]{taxonget()}}: service no longer provided
\item \code{\link[=taxonsearch]{taxonsearch()}}: service no longer provided
\item \code{\link[=stylegeojson]{stylegeojson()}}: moving this functionality to spocc package, will be
removed soon
\item \code{\link[=togeojson]{togeojson()}}: moving this functionality to spocc package, will be
removed soon
\item \code{\link[=gist]{gist()}}: moving this functionality to spocc package, will be
removed soon
\item \code{\link[=occ_spellcheck]{occ_spellcheck()}}: GBIF has removed the \code{spellCheck} parameter
from their API
}
}
\details{
The above functions have been removed. See
\url{https://github.com/ropensci/rgbif} and poke around the code if you
want to find the old functions in previous versions of the package
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_wkt.r
\name{check_wkt}
\alias{check_wkt}
\title{Check input WKT}
\usage{
check_wkt(wkt = NULL, skip_validate = FALSE)
}
\arguments{
\item{wkt}{(character) one or more Well Known Text objects}

\item{skip_validate}{(logical) whether to skip \code{wellknown::validate_wkt}
call or not. Default: \code{FALSE}}
}
\description{
Check input WKT
}
\examples{
\dontrun{
check_wkt('POLYGON((30.1 10.1, 10 20, 20 60, 60 60, 30.1 10.1))')
check_wkt('POINT(30.1 10.1)')
check_wkt('LINESTRING(3 4,10 50,20 25)')

# check many passed in at once
check_wkt(c('POLYGON((30.1 10.1, 10 20, 20 60, 60 60, 30.1 10.1))',
  'POINT(30.1 10.1)'))

# bad WKT
# wkt <- 'POLYGON((30.1 10.1, 10 20, 20 60, 60 60, 30.1 a))'
# check_wkt(wkt)

# many wkt's, semi-colon separated, for many repeated "geometry" args
wkt <- "POLYGON((-102.2 46.0,-93.9 46.0,-93.9 43.7,-102.2 43.7,-102.2 46.0))
;POLYGON((30.1 10.1, 10 20, 20 40, 40 40, 30.1 10.1))"
check_wkt(gsub("\n", '', wkt))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_search.r
\name{occ_search}
\alias{occ_search}
\title{Search for GBIF occurrences}
\usage{
occ_search(
  taxonKey = NULL,
  scientificName = NULL,
  country = NULL,
  publishingCountry = NULL,
  hasCoordinate = NULL,
  typeStatus = NULL,
  recordNumber = NULL,
  lastInterpreted = NULL,
  continent = NULL,
  geometry = NULL,
  geom_big = "asis",
  geom_size = 40,
  geom_n = 10,
  recordedBy = NULL,
  recordedByID = NULL,
  identifiedByID = NULL,
  basisOfRecord = NULL,
  datasetKey = NULL,
  eventDate = NULL,
  catalogNumber = NULL,
  year = NULL,
  month = NULL,
  decimalLatitude = NULL,
  decimalLongitude = NULL,
  elevation = NULL,
  depth = NULL,
  institutionCode = NULL,
  collectionCode = NULL,
  hasGeospatialIssue = NULL,
  issue = NULL,
  search = NULL,
  mediaType = NULL,
  subgenusKey = NULL,
  repatriated = NULL,
  phylumKey = NULL,
  kingdomKey = NULL,
  classKey = NULL,
  orderKey = NULL,
  familyKey = NULL,
  genusKey = NULL,
  establishmentMeans = NULL,
  protocol = NULL,
  license = NULL,
  organismId = NULL,
  publishingOrg = NULL,
  stateProvince = NULL,
  waterBody = NULL,
  locality = NULL,
  limit = 500,
  start = 0,
  fields = "all",
  return = NULL,
  facet = NULL,
  facetMincount = NULL,
  facetMultiselect = NULL,
  skip_validate = TRUE,
  curlopts = list(),
  ...
)
}
\arguments{
\item{taxonKey}{(numeric) A taxon key from the GBIF backbone. All included and synonym taxa
are included in the search, so a search for aves with taxononKey=212
(i.e. /occurrence/search?taxonKey=212) will match all birds, no matter which
species. You can pass many keys by passing occ_search in a call to an
lapply-family function (see last example below).}

\item{scientificName}{A scientific name from the GBIF backbone. All included and synonym
taxa are included in the search.}

\item{country}{The 2-letter country code (as per ISO-3166-1) of the country in
which the occurrence was recorded. See here
https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}

\item{publishingCountry}{The 2-letter country code (as per ISO-3166-1) of the
country in which the occurrence was recorded.}

\item{hasCoordinate}{(logical) Return only occurence records with lat/long data (\code{TRUE}) or
all records (\code{FALSE}, default).}

\item{typeStatus}{Type status of the specimen. One of many options. See \code{?typestatus}}

\item{recordNumber}{Number recorded by collector of the data, different from GBIF record
number. See http://rs.tdwg.org/dwc/terms/#recordNumber for more info}

\item{lastInterpreted}{Date the record was last modified in GBIF, in ISO 8601 format:
yyyy, yyyy-MM, yyyy-MM-dd, or MM-dd.  Supports range queries, smaller,larger (e.g.,
'1990,1991', whereas '1991,1990' wouldn't work)}

\item{continent}{Continent. One of africa, antarctica, asia, europe, north_america
(North America includes the Caribbean and reachies down and includes Panama), oceania,
or south_america}

\item{geometry}{Searches for occurrences inside a polygon described in Well Known
Text (WKT) format. A WKT shape written as either POINT, LINESTRING, LINEARRING
POLYGON, or MULTIPOLYGON. Example of a polygon: POLYGON((30.1 10.1, 20, 20 40, 40 40, 30.1 10.1))
would be queried as http://bit.ly/1BzNwDq See also the section \strong{WKT} below.}

\item{geom_big}{(character) One of "axe", "bbox", or "asis" (default). See Details.}

\item{geom_size}{(integer) An integer indicating size of the cell. Default: 40. See Details.}

\item{geom_n}{(integer) An integer indicating number of cells in each dimension. Default: 10.
See Details.}

\item{recordedBy}{The person who recorded the occurrence.}

\item{recordedByID}{(character) Identifier (e.g. ORCID) for the person who
recorded the occurrence}

\item{identifiedByID}{(character) Identifier (e.g. ORCID) for the person who
provided the taxonomic identification of the occurrence.}

\item{basisOfRecord}{Basis of record, as defined in our BasisOfRecord enum here
https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html
Acceptable values are:
\itemize{
\item FOSSIL_SPECIMEN An occurrence record describing a fossilized specimen.
\item HUMAN_OBSERVATION An occurrence record describing an observation made by
one or more people.
\item LITERATURE An occurrence record based on literature alone.
\item LIVING_SPECIMEN An occurrence record describing a living specimen, e.g.
\item MACHINE_OBSERVATION An occurrence record describing an observation made
by a machine.
\item OBSERVATION An occurrence record describing an observation.
\item PRESERVED_SPECIMEN An occurrence record describing a preserved specimen.
\item UNKNOWN Unknown basis for the record.
}}

\item{datasetKey}{The occurrence dataset key (a uuid)}

\item{eventDate}{Occurrence date in ISO 8601 format: yyyy, yyyy-MM, yyyy-MM-dd, or
MM-dd. Supports range queries, smaller,larger (e.g., '1990,1991', whereas '1991,1990'
wouldn't work)}

\item{catalogNumber}{An identifier of any form assigned by the source within a
physical collection or digital dataset for the record which may not unique,
but should be fairly unique in combination with the institution and collection code.}

\item{year}{The 4 digit year. A year of 98 will be interpreted as AD 98. Supports range queries,
smaller,larger (e.g., '1990,1991', whereas '1991,1990' wouldn't work)}

\item{month}{The month of the year, starting with 1 for January. Supports range queries,
smaller,larger (e.g., '1,2', whereas '2,1' wouldn't work)}

\item{decimalLatitude}{Latitude in decimals between -90 and 90 based on WGS 84.
Supports range queries, smaller,larger (e.g., '25,30', whereas '30,25' wouldn't work)}

\item{decimalLongitude}{Longitude in decimals between -180 and 180 based on WGS 84.
Supports range queries (e.g., '-0.4,-0.2', whereas '-0.2,-0.4' wouldn't work).}

\item{elevation}{Elevation in meters above sea level. Supports range queries, smaller,larger
(e.g., '5,30', whereas '30,5' wouldn't work)}

\item{depth}{Depth in meters relative to elevation. For example 10 meters below a
lake surface with given elevation. Supports range queries, smaller,larger (e.g., '5,30',
whereas '30,5' wouldn't work)}

\item{institutionCode}{An identifier of any form assigned by the source to identify
the institution the record belongs to. Not guaranteed to be que.}

\item{collectionCode}{An identifier of any form assigned by the source to identify
the physical collection or digital dataset uniquely within the text of an institution.}

\item{hasGeospatialIssue}{(logical) Includes/excludes occurrence records which contain spatial
issues (as determined in our record interpretation), i.e. \code{hasGeospatialIssue=TRUE}
returns only those records with spatial issues while \code{hasGeospatialIssue=FALSE} includes
only records without spatial issues. The absence of this parameter returns any
record with or without spatial issues.}

\item{issue}{(character) One or more of many possible issues with each occurrence record. See
Details. Issues passed to this parameter filter results by the issue.}

\item{search}{Query terms. The value for this parameter can be a simple word or a phrase.}

\item{mediaType}{Media type. Default is NULL, so no filtering on mediatype. Options:
NULL, 'MovingImage', 'Sound', and 'StillImage'.}

\item{subgenusKey}{(numeric) Subgenus classification key.}

\item{repatriated}{(character) Searches for records whose publishing country
is different to the country where the record was recorded in.}

\item{phylumKey}{(numeric) Phylum classification key.}

\item{kingdomKey}{(numeric) Kingdom classification key.}

\item{classKey}{(numeric) Class classification key.}

\item{orderKey}{(numeric) Order classification key.}

\item{familyKey}{(numeric) Family classification key.}

\item{genusKey}{(numeric) Genus classification key.}

\item{establishmentMeans}{(character) EstablishmentMeans, possible values
include: INTRODUCED, INVASIVE, MANAGED, NATIVE, NATURALISED, UNCERTAIN}

\item{protocol}{(character) Protocol or mechanism used to provide the
occurrence record. See Details for possible values}

\item{license}{(character) The type license applied to the dataset or record.
Possible values: CC0_1_0, CC_BY_4_0, CC_BY_NC_4_0, UNSPECIFIED, and
UNSUPPORTED}

\item{organismId}{(numeric) An identifier for the Organism instance (as
opposed to a particular digital record of the Organism). May be a globally
unique identifier or an identifier specific to the data set.}

\item{publishingOrg}{(character) The publishing organization key (a UUID).}

\item{stateProvince}{(character) The name of the next smaller administrative
region than country (state, province, canton, department, region, etc.) in
which the Location occurs.}

\item{waterBody}{(character) The name of the water body in which the
locations occur}

\item{locality}{(character) The specific description of the place.}

\item{limit}{Number of records to return. Default: 500. Note that the per
request maximum is 300, but since we set it at 500 for the function, we
do two requests to get you the 500 records (if there are that many).
Note that there is a hard maximum of 100,000, which is calculated as the
\code{limit+start}, so \code{start=99,000} and \code{limit=2000} won't work}

\item{start}{Record number to start at. Use in combination with limit to
page through results. Note that we do the paging internally for you, but
you can manually set the \code{start} parameter}

\item{fields}{(character) Default ('all') returns all fields. 'minimal'
returns just taxon name, key, latitude, and longitute. Or specify each field
you want returned by name, e.g. fields = c('name','latitude','elevation').}

\item{return}{Defunct. All components (meta, hierarchy, data, media,
facets) are returned now; index to the one(s) you want. See \code{\link[=occ_data]{occ_data()}}
if you just want the data component}

\item{facet}{(character) a character vector of length 1 or greater.
Required.}

\item{facetMincount}{(numeric) minimum number of records to be included
in the faceting results}

\item{facetMultiselect}{(logical) Set to \code{TRUE} to still return counts
for values that are not currently filtered. See examples.
Default: \code{FALSE}

\strong{Faceting}:
All fields can be faceted on except for last "lastInterpreted",
"eventDate", and "geometry"

You can do facet searches alongside searching occurrence data, and
return both, or only return facets, or only occurrence data, etc.}

\item{skip_validate}{(logical) whether to skip \code{wellknown::validate_wkt}
call or not. passed down to \code{\link[=check_wkt]{check_wkt()}}. Default: \code{TRUE}}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}

\item{...}{additional facet parameters}
}
\value{
An object of class \code{gbif}, which is a S3 class list, with
slots for metadata (\code{meta}), the occurrence data itself (\code{data}),
the taxonomic hierarchy data (\code{hier}), and media metadata
(\code{media}).
In addition, the object has attributes listing the user supplied arguments
and whether it was a 'single' or 'many' search; that is, if you supply two
values of the \code{datasetKey} parameter to searches are done, and it's a
'many'. \code{meta} is a list of length four with offset, limit,
endOfRecords and count fields. \code{data} is a tibble (aka data.frame). \code{hier}
is a list of data.frames of the unique set of taxa found, where each
data.frame is its taxonomic classification. \code{media} is a list of media
objects, where each element holds a set of metadata about the media object.
}
\description{
Search for GBIF occurrences
}
\note{
Maximum number of records you can get with this function is 100,000.
See https://www.gbif.org/developer/occurrence
}
\section{protocol parameter options}{

\itemize{
\item BIOCASE - A BioCASe protocl compliant service.
\item DIGIR - A DiGIR service endpoint.
\item DIGIR_MANIS - A DiGIR service slightly modified for the MANIS
network.
\item DWC_ARCHIVE - A Darwin Core Archive as defined by the Darwin Core
Text Guidelines.
\item EML - A single EML metadata document in any EML version.
\item FEED - Syndication feeds like RSS or ATOM of various flavors.
\item OAI_PMH - The Open Archives Initiative Protocol for Metadata
Harvesting.
\item OTHER - Any other service not covered by this enum so far.
\item TAPIR - A TAPIR service.
\item TCS_RDF - Taxon Concept data given as RDF based on the TDWG ontology.
\item TCS_XML - A Taxon Concept Schema document.
\item WFS - An OGC Web Feature Service.
\item WMS - An OGC Web Map Service.
}
}

\section{Multiple values passed to a parameter}{

There are some parameters you can pass multiple values to in a vector,
each value of which produces a different request (MDR: multiple different
requests). Some parameters allow multiple values to be passed in the same
request (MSR: multiple same request) in a semicolon separated string
(e.g., 'a;b'); if given we'll do a single request with that parameter
repeated for each value given (e.g., \code{foo=a&foo=b} if the parameter
is \code{foo}). Some parameters allow both MDR and MSR.

The following list shows which parameters support MDR and MSR.

\itemize{
\item basisOfRecord: MDR, MSR
\item classKey: MDR, MSR
\item country: MDR, MSR
\item familyKey: MDR, MSR
\item genusKey: MDR, MSR
\item identifiedByID: MDR, MSR
\item kingdomKey: MDR, MSR
\item license: MDR, MSR
\item locality: MDR, MSR
\item catalogNumber: MDR, MSR
\item collectionCode: MDR, MSR
\item continent: MDR, MSR
\item datasetKey: MDR, MSR
\item establishmentMeans: MDR, MSR
\item geometry: MDR, MSR
\item institutionCode: MDR, MSR
\item mediaType: MDR, MSR
\item orderKey: MDR, MSR
\item organismId: MDR, MSR
\item phylumKey: MDR, MSR
\item protocol: MDR, MSR
\item publishingCountry: MDR, MSR
\item publishingOrg: MDR, MSR
\item recordedBy: MDR, MSR
\item recordedByID: MDR, MSR
\item recordNumber: MDR, MSR
\item scientificName: MDR, MSR
\item stateProvince: MDR, MSR
\item subgenusKey: MDR, MSR
\item taxonKey: MDR, MSR
\item typeStatus: MDR, MSR
\item waterBody: MDR, MSR
\item depth: MDR
\item limit: MDR
\item q: MDR
\item year: MDR
\item repatriated: MDR
\item lastInterpreted: MDR
\item decimalLatitude: MDR
\item decimalLongitude: MDR
}

Note that you can not pass a vector > length 1 to more than 1 of the above
MDR parameters at the same time.

see also \code{\link{many-values}}
}

\section{Range queries}{

A range query is as it sounds - you query on a range of values defined by
a lower and upper limit. Do a range query by specifying the lower and upper
limit in a string like \code{depth='50,100'}. It would be more R like to
specify the range in a vector like \code{c(50,100)}, but that sort of syntax
allows you to do many searches, one for each element in the vector -
thus range queries have to differ. The following parameters support
range queries.

\itemize{
\item decimalLatitude
\item decimalLongitude
\item depth
\item elevation
\item eventDate
\item lastInterpreted
\item month
\item year
}

See also above section: semicolon and comma separated strings lead to
different outcomes for some parameters.
}

\section{Hierarchies}{

Hierarchies are returned wih each occurrence object. There is no
option no to return them from the API. However, within the \code{occ_search}
function you can select whether to return just hierarchies, just data, all
of data and hiearchies and metadata, or just metadata. If all hierarchies
are the same we just return one for you.
}

\section{curl debugging}{

You can pass parameters not defined in this function into the call to
the GBIF API to control things about the call itself using \code{curlopts}.
See an example below that passes in the \code{verbose} function to get
details on the http call.
}

\section{Scientific names vs. taxon keys}{

In the previous GBIF API and the version of rgbif that
wrapped that API, you could search the equivalent of this function with a
species name, which was convenient. However, names are messy right. So it
sorta makes sense to sort out the species key numbers you want exactly,
and then get your occurrence data with this function. GBIF has added a
parameter scientificName to allow searches by scientific names in this
function - which includes synonym taxa. \emph{Note:} that if you do use the
scientificName parameter, we will check internally that it's not a
synonym of an accepted name, and if it is, we'll search on the
accepted name. If you want to force searching by a synonym do so by
finding the GBIF identifier first with any \code{name_*} functions,
then pass that ID to the \code{taxonKey} parameter.
}

\section{WKT}{

Examples of valid WKT objects:
\itemize{
\item 'POLYGON((-19.5 34.1, 27.8 34.1, 35.9 68.1, -25.3 68.1, -19.5 34.1))'
\item 'MULTIPOLYGON(((-123 38,-116 38,-116 43,-123 43,-123 38)),((-97 41,-93 41,-93 45,-97 45,-97 41)))'
\item 'POINT(-120 40)'
\item 'LINESTRING(3 4,10 50,20 25)'
\item 'LINEARRING' ???' - Not sure how to specify this. Anyone?
}

Note that GBIF expects counter-clockwise winding order for WKT. You can
supply clockwise WKT, but GBIF treats it as an exclusion, so you get all
data not inside the WKT area. \code{\link[=occ_download]{occ_download()}} behaves differently
in that you should simply get no data back at all with clockwise WKT.
}

\section{Long WKT}{

Options for handling long WKT strings:
Note that long WKT strings are specially handled when using \code{\link{occ_search}} or
\code{\link{occ_data}}. Here are the three options for long WKT strings (> 1500 characters),
set one of these three via the parameter \code{geom_big}:
\itemize{
\item asis - the default setting. This means we don't do anything internally. That is,
we just pass on your WKT string just as we've done before in this package.
\item axe - this option uses the \pkg{sf} package to chop up your WKT string in
to many polygons, which then leads to a separate data request for each polygon piece,
then we combine all dat back together to give to you. Note that if your WKT string
is not of type polygon, we drop back to \code{asis}as there's no way to chop up
linestrings, etc. This option will in most cases be slower than the other two options.
However, this polygon splitting approach won't have the problem of
the disconnect between how many records you want and what you actually get back as
with the bbox option.

This method uses \code{sf::st_make_grid} and \code{sf::st_intersection}, which has
two parameters \code{cellsize} and \code{n}. You can tweak those parameters here by
tweaking \code{geom_size} and \code{geom_n}. \code{geom_size} seems to be more useful in
toggling the number of WKT strings you get back.

See \code{\link{wkt_parse}} to manually break make WKT bounding box from a larger WKT
string, or break a larger WKT string into many smaller ones.

\item bbox - this option checks whether your WKT string is longer than 1500 characters,
and if it is we create a bounding box from the WKT, do the GBIF search with that
bounding box, then prune the resulting data to only those occurrences in your original
WKT string. There is a big caveat however. Because we create a bounding box from the WKT,
and the \code{limit} parameter determines some subset of records to get, then when we
prune the resulting data to the WKT, the number of records you get could be less than
what you set with your \code{limit} parameter. However, you could set the limit to be
high enough so that you get all records back found in that bounding box, then you'll
get all the records available within the WKT.
}
}

\section{issue parameter}{

The options for the issue parameter (from
https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html):
\itemize{
\item BASIS_OF_RECORD_INVALID The given basis of record is impossible to interpret or seriously
different from the recommended vocabulary.
\item CONTINENT_COUNTRY_MISMATCH The interpreted continent and country do not match up.
\item CONTINENT_DERIVED_FROM_COORDINATES The interpreted continent is based on the coordinates,
not the verbatim string information.
\item CONTINENT_INVALID Uninterpretable continent values found.
\item COORDINATE_INVALID Coordinate value given in some form but GBIF is unable to interpret it.
\item COORDINATE_OUT_OF_RANGE Coordinate has invalid lat/lon values out of their decimal max
range.
\item COORDINATE_REPROJECTED The original coordinate was successfully reprojected from a
different geodetic datum to WGS84.
\item COORDINATE_REPROJECTION_FAILED The given decimal latitude and longitude could not be
reprojected to WGS84 based on the provided datum.
\item COORDINATE_REPROJECTION_SUSPICIOUS Indicates successful coordinate reprojection according
to provided datum, but which results in a datum shift larger than 0.1 decimal degrees.
\item COORDINATE_ROUNDED Original coordinate modified by rounding to 5 decimals.
\item COUNTRY_COORDINATE_MISMATCH The interpreted occurrence coordinates fall outside of the
indicated country.
\item COUNTRY_DERIVED_FROM_COORDINATES The interpreted country is based on the coordinates, not
the verbatim string information.
\item COUNTRY_INVALID Uninterpretable country values found.
\item COUNTRY_MISMATCH Interpreted country for dwc:country and dwc:countryCode contradict each
other.
\item DEPTH_MIN_MAX_SWAPPED Set if supplied min>max
\item DEPTH_NON_NUMERIC Set if depth is a non numeric value
\item DEPTH_NOT_METRIC Set if supplied depth is not given in the metric system, for example
using feet instead of meters
\item DEPTH_UNLIKELY Set if depth is larger than 11.000m or negative.
\item ELEVATION_MIN_MAX_SWAPPED Set if supplied min > max elevation
\item ELEVATION_NON_NUMERIC Set if elevation is a non numeric value
\item ELEVATION_NOT_METRIC Set if supplied elevation is not given in the metric system, for
example using feet instead of meters
\item ELEVATION_UNLIKELY Set if elevation is above the troposphere (17km) or below 11km
(Mariana Trench).
\item GEODETIC_DATUM_ASSUMED_WGS84 Indicating that the interpreted coordinates assume they are
based on WGS84 datum as the datum was either not indicated or interpretable.
\item GEODETIC_DATUM_INVALID The geodetic datum given could not be interpreted.
\item IDENTIFIED_DATE_INVALID The date given for dwc:dateIdentified is invalid and cant be
interpreted at all.
\item IDENTIFIED_DATE_UNLIKELY The date given for dwc:dateIdentified is in the future or before
Linnean times (1700).
\item MODIFIED_DATE_INVALID A (partial) invalid date is given for dc:modified, such as a non
existing date, invalid zero month, etc.
\item MODIFIED_DATE_UNLIKELY The date given for dc:modified is in the future or predates unix
time (1970).
\item MULTIMEDIA_DATE_INVALID An invalid date is given for dc:created of a multimedia object.
\item MULTIMEDIA_URI_INVALID An invalid uri is given for a multimedia object.
\item PRESUMED_NEGATED_LATITUDE Latitude appears to be negated, e.g. 32.3 instead of -32.3
\item PRESUMED_NEGATED_LONGITUDE Longitude appears to be negated, e.g. 32.3 instead of -32.3
\item PRESUMED_SWAPPED_COORDINATE Latitude and longitude appear to be swapped.
\item RECORDED_DATE_INVALID A (partial) invalid date is given, such as a non existing date,
invalid zero month, etc.
\item RECORDED_DATE_MISMATCH The recording date specified as the eventDate string and the
individual year, month, day are contradicting.
\item RECORDED_DATE_UNLIKELY The recording date is highly unlikely, falling either into the
future or represents a very old date before 1600 that predates modern taxonomy.
\item REFERENCES_URI_INVALID An invalid uri is given for dc:references.
\item TAXON_MATCH_FUZZY Matching to the taxonomic backbone can only be done using a fuzzy, non
exact match.
\item TAXON_MATCH_HIGHERRANK Matching to the taxonomic backbone can only be done on a higher
rank and not the scientific name.
\item TAXON_MATCH_NONE Matching to the taxonomic backbone cannot be done cause there was no
match at all or several matches with too little information to keep them apart (homonyms).
\item TYPE_STATUS_INVALID The given type status is impossible to interpret or seriously
different from the recommended vocabulary.
\item ZERO_COORDINATE Coordinate is the exact 0/0 coordinate, often indicating a bad null
coordinate.
}
}

\section{Counts}{

There is a slight difference in the way records are counted here vs.
results from \code{\link{occ_count}}. For equivalent outcomes, in this
function use \code{hasCoordinate=TRUE}, and \code{hasGeospatialIssue=FALSE}
to have the same outcome using \code{\link{occ_count}} with
\code{isGeoreferenced=TRUE}
}

\examples{
\dontrun{
# Search by species name, using \code{\link{name_backbone}} first to get key
(key <- name_suggest(q='Helianthus annuus', rank='species')$data$key[1])
occ_search(taxonKey=key, limit=2)

# Return 20 results, this is the default by the way
occ_search(taxonKey=key, limit=20)

# Get just metadata
occ_search(taxonKey=key, limit=0)$meta

# Instead of getting a taxon key first, you can search for a name directly
## However, note that using this approach (with \code{scientificName="..."})
## you are getting synonyms too. The results for using \code{scientifcName} and
## \code{taxonKey} parameters are the same in this case, but I wouldn't be surprised if for some
## names they return different results
occ_search(scientificName = 'Ursus americanus')
key <- name_backbone(name = 'Ursus americanus', rank='species')$usageKey
occ_search(taxonKey = key)

# Search by dataset key
occ_search(datasetKey='7b5d6a48-f762-11e1-a439-00145eb45e9a', limit=20)$data

# Search by catalog number
occ_search(catalogNumber="49366", limit=20)
## separate requests: use a vector of strings
occ_search(catalogNumber=c("49366","Bird.27847588"), limit=10)
## one request, many instances of same parameter: use semi-colon sep. string
occ_search(catalogNumber="49366;Bird.27847588", limit=10)

# Get all data, not just lat/long and name
occ_search(taxonKey=key, fields='all', limit=20)

# Or get specific fields. Note that this isn't done on GBIF's side of things. This
# is done in R, but before you get the return object, so other fields are garbage
# collected
occ_search(taxonKey=key, fields=c('name','basisOfRecord','protocol'), limit=20)

# Use paging parameters (limit and start) to page. Note the different results
# for the two queries below.
occ_search(datasetKey='7b5d6a48-f762-11e1-a439-00145eb45e9a',start=10,limit=5)$data
occ_search(datasetKey='7b5d6a48-f762-11e1-a439-00145eb45e9a',start=20,limit=5)$data

# Many dataset keys
## separate requests: use a vector of strings
occ_search(datasetKey=c("50c9509d-22c7-4a22-a47d-8c48425ef4a7",
   "7b5d6a48-f762-11e1-a439-00145eb45e9a"), limit=20)
## one request, many instances of same parameter: use semi-colon sep. string
v="50c9509d-22c7-4a22-a47d-8c48425ef4a7;7b5d6a48-f762-11e1-a439-00145eb45e9a"
occ_search(datasetKey = v, limit=20)

# Occurrence data: lat/long data, and associated metadata with occurrences
## The `data` slot has a data.frame of all data together
## for easy manipulation
occ_search(taxonKey=key, limit=20)$data

# Taxonomic hierarchy data
## In the `hier` slot
occ_search(taxonKey=key, limit=10)$hier

# Search by recorder
occ_search(recordedBy="smith", limit=20)

# Many collector names
occ_search(recordedBy=c("smith","BJ Stacey"), limit=20)

# recordedByID
occ_search(recordedByID="https://orcid.org/0000-0003-1691-239X", limit=20)

# identifiedByID
occ_search(identifiedByID="https://orcid.org/0000-0003-4710-2648", limit=20)

# Pass in curl options for extra fun
occ_search(taxonKey=2433407, limit=20, curlopts=list(verbose=TRUE))$hier
occ_search(taxonKey=2433407, limit=20,
  curlopts = list(
    noprogress = FALSE,
    progressfunction = function(down, up) {
      cat(sprintf("up: \%d | down \%d\n", up, down))
      return(TRUE)
    }
  )
)$hier
# occ_search(taxonKey=2433407, limit=20,
#   curlopts = list(timeout_ms = 1))

# Search for many species
splist <- c('Cyanocitta stelleri', 'Junco hyemalis', 'Aix sponsa')
keys <- sapply(splist, function(x) name_suggest(x)$data$key[1], USE.NAMES=FALSE)
## separate requests: use a vector of strings
occ_search(taxonKey = keys, limit=5)
## one request, many instances of same parameter: use semi-colon sep. string
occ_search(taxonKey = paste0(keys, collapse = ";"), limit=5)

# Search using a synonym name
#  Note that you'll see a message printing out that the accepted name will be used
occ_search(scientificName = 'Pulsatilla patens', fields = c('name','scientificName'), limit=5)

# Search on latitidue and longitude
occ_search(decimalLatitude=48, decimalLongitude=10)

# Search on a bounding box
## in well known text format
### polygon
occ_search(geometry='POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))', limit=20)
### multipolygon
wkt <- 'MULTIPOLYGON(((-123 38,-116 38,-116 43,-123 43,-123 38)),
   ((-97 41,-93 41,-93 45,-97 45,-97 41)))'
occ_search(geometry = gsub("\n\\\\s+", "", wkt), limit = 20)

## taxonKey + WKT
key <- name_suggest(q='Aesculus hippocastanum')$data$key[1]
occ_search(taxonKey=key, geometry='POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))',
   limit=20)
## or using bounding box, converted to WKT internally
occ_search(geometry=c(-125.0,38.4,-121.8,40.9), limit=20)

# Search on a long WKT string - too long for a GBIF search API request
## We internally convert your WKT string to a bounding box
##  then do the query
##  then clip the results down to just those in the original polygon
##  - Alternatively, you can set the parameter `geom_big="bbox"`
##  - An additional alternative is to use the GBIF download API, see ?downloads
wkt <- "POLYGON((-9.178796777343678 53.22769021556159,
-12.167078027343678 51.56540789297837,
-12.958093652343678 49.78333685689162,-11.024499902343678 49.21251756301334,
-12.079187402343678 46.68179685941719,-15.067468652343678 45.83103608186854,
-15.770593652343678 43.58271629699817,-15.067468652343678 41.57676278827219,
-11.815515527343678 40.44938999172728,-12.958093652343678 37.72112962230871,
-11.639734277343678 36.52987439429357,-8.299890527343678 34.96062625095747,
-8.739343652343678 32.62357394385735,-5.223718652343678 30.90497915232165,
1.1044063476563224 31.80562077746643,1.1044063476563224 30.754036557416256,
6.905187597656322 32.02942785462211,5.147375097656322 32.99292810780193,
9.629796972656322 34.164474406524725,10.860265722656322 32.91918014319603,
14.551671972656322 33.72700959356651,13.409093847656322 34.888564192275204,
16.748937597656322 35.104560368110114,19.561437597656322 34.81643887792552,
18.594640722656322 36.38849705969625,22.989171972656322 37.162874858929854,
19.825109472656322 39.50651757842751,13.760656347656322 38.89353140585116,
14.112218847656322 42.36091601976124,10.596593847656322 41.11488736647705,
9.366125097656322 43.70991402658437,5.059484472656322 42.62015372417812,
2.3348750976563224 45.21526500321446,-0.7412967773436776 46.80225692528942,
6.114171972656322 47.102229890207894,8.047765722656322 45.52399303437107,
12.881750097656322 48.22681126957933,9.190343847656322 48.693079457106684,
8.750890722656322 50.68283120621287,5.059484472656322 50.40356146487845,
4.268468847656322 52.377558897655156,1.4559688476563224 53.28027243658647,
0.8407344726563224 51.62000971578333,0.5770625976563224 49.32721423860726,
-2.5869999023436776 49.49875947592088,-2.4991092773436776 51.18135535408638,
-2.0596561523436776 52.53822562473851,-4.696374902343678 51.67454591918756,
-5.311609277343678 50.009802108095776,-6.629968652343678 48.75106196817059,
-7.684656152343678 50.12263634382465,-6.190515527343678 51.83776110910459,
-5.047937402343678 54.267098895684235,-6.893640527343678 53.69860705549198,
-8.915124902343678 54.77719740243195,-12.079187402343678 54.52294465763567,
-13.573328027343678 53.437631551347174,
-11.288171777343678 53.48995552517918,
-9.178796777343678 53.22769021556159))"
wkt <- gsub("\n", " ", wkt)

#### Default option with large WKT string fails
# res <- occ_search(geometry = wkt)

#### if WKT too long, with 'geom_big=bbox': makes into bounding box
res <- occ_search(geometry = wkt, geom_big = "bbox")$data
library("rgeos")
library("sp")
wktsp <- readWKT(wkt)
plot(wktsp)
coordinates(res) <- ~decimalLongitude+decimalLatitude
points(res)

#### Or, use 'geom_big=axe'
(res <- occ_search(geometry = wkt, geom_big = "axe"))
##### manipulate essentially number of polygons that result, so number of requests
###### default geom_size is 40
###### fewer calls
(res <- occ_search(geometry = wkt, geom_big = "axe", geom_size=50))
###### more calls
(res <- occ_search(geometry = wkt, geom_big = "axe", geom_size=30))


# Search on country
occ_search(country='US', fields=c('name','country'), limit=20)
isocodes[grep("France", isocodes$name),"code"]
occ_search(country='FR', fields=c('name','country'), limit=20)
occ_search(country='DE', fields=c('name','country'), limit=20)
### separate requests: use a vector of strings
occ_search(country=c('US','DE'), limit=20)
### one request, many instances of same parameter: use semi-colon sep. string
occ_search(country = 'US;DE', limit=20)

# Get only occurrences with lat/long data
occ_search(taxonKey=key, hasCoordinate=TRUE, limit=20)

# Get only occurrences that were recorded as living specimens
occ_search(taxonKey=key, basisOfRecord="LIVING_SPECIMEN", hasCoordinate=TRUE, limit=20)
## multiple values in a vector = a separate request for each value
occ_search(taxonKey=key,
  basisOfRecord=c("LIVING_SPECIMEN", "HUMAN_OBSERVATION"), limit=20)
## mutiple values in a single string, ";" separated = one request including all values
occ_search(taxonKey=key,
  basisOfRecord="LIVING_SPECIMEN;HUMAN_OBSERVATION", limit=20)

# Get occurrences for a particular eventDate
occ_search(taxonKey=key, eventDate="2013", limit=20)
occ_search(taxonKey=key, year="2013", limit=20)
occ_search(taxonKey=key, month="6", limit=20)

# Get occurrences based on depth
key <- name_backbone(name='Salmo salar', kingdom='animals')$speciesKey
occ_search(taxonKey=key, depth="5", limit=20)

# Get occurrences based on elevation
key <- name_backbone(name='Puma concolor', kingdom='animals')$speciesKey
occ_search(taxonKey=key, elevation=50, hasCoordinate=TRUE, limit=20)

# Get occurrences based on institutionCode
occ_search(institutionCode="TLMF", limit=20)
### separate requests: use a vector of strings
occ_search(institutionCode=c("TLMF","ArtDatabanken"), limit=20)
### one request, many instances of same parameter: use semi-colon sep. string
occ_search(institutionCode = "TLMF;ArtDatabanken", limit=20)

# Get occurrences based on collectionCode
occ_search(collectionCode="Floristic Databases MV - Higher Plants", limit=20)
occ_search(collectionCode=c("Floristic Databases MV - Higher Plants","Artport"))

# Get only those occurrences with spatial issues
occ_search(taxonKey=key, hasGeospatialIssue=TRUE, limit=20)

# Search using a query string
occ_search(search = "kingfisher", limit=20)

# search on repatriated - doesn't work right now
# occ_search(repatriated = "")

# search on phylumKey
occ_search(phylumKey = 7707728, limit = 5)

# search on kingdomKey
occ_search(kingdomKey = 1, limit = 5)

# search on classKey
occ_search(classKey = 216, limit = 5)

# search on orderKey
occ_search(orderKey = 7192402, limit = 5)

# search on familyKey
occ_search(familyKey = 3925, limit = 5)

# search on genusKey
occ_search(genusKey = 1935496, limit = 5)

# search on establishmentMeans
occ_search(establishmentMeans = "INVASIVE", limit = 5)
occ_search(establishmentMeans = "NATIVE", limit = 5)
occ_search(establishmentMeans = "UNCERTAIN", limit = 5)

# search on protocol
occ_search(protocol = "DIGIR", limit = 5)

# search on license
occ_search(license = "CC_BY_4_0", limit = 5)

# search on organismId
occ_search(organismId = "100", limit = 5)

# search on publishingOrg
occ_search(publishingOrg = "28eb1a3f-1c15-4a95-931a-4af90ecb574d", limit = 5)

# search on stateProvince
occ_search(stateProvince = "California", limit = 5)

# search on waterBody
occ_search(waterBody = "AMAZONAS BASIN, RIO JURUA", limit = 5)

# search on locality
res <- occ_search(locality = c("Trondheim", "Hovekilen"), limit = 5)
res$Trondheim$data
res$Hovekilen$data



# Range queries
## See Detail for parameters that support range queries
occ_search(depth='50,100') # this is a range depth, with lower/upper limits in character string
occ_search(depth=c(50,100)) # this is not a range search, but does two searches for each depth

## Range search with year
occ_search(year='1999,2000', limit=20)

## Range search with latitude
occ_search(decimalLatitude='29.59,29.6')

# Search by specimen type status
## Look for possible values of the typeStatus parameter looking at the typestatus dataset
occ_search(typeStatus = 'allotype', fields = c('name','typeStatus'))

# Search by specimen record number
## This is the record number of the person/group that submitted the data, not GBIF's numbers
## You can see that many different groups have record number 1, so not super helpful
occ_search(recordNumber = 1, fields = c('name','recordNumber','recordedBy'))

# Search by last time interpreted: Date the record was last modified in GBIF
## The lastInterpreted parameter accepts ISO 8601 format dates, including
## yyyy, yyyy-MM, yyyy-MM-dd, or MM-dd. Range queries are accepted for lastInterpreted
occ_search(lastInterpreted = '2014-04-02', fields = c('name','lastInterpreted'))

# Search by continent
## One of africa, antarctica, asia, europe, north_america, oceania, or south_america
occ_search(continent = 'south_america')$meta
occ_search(continent = 'africa')$meta
occ_search(continent = 'oceania')$meta
occ_search(continent = 'antarctica')$meta

# Search for occurrences with images
occ_search(mediaType = 'StillImage')$media
occ_search(mediaType = 'MovingImage')$media
occ_search(mediaType = 'Sound')$media

# Query based on issues - see Details for options
## one issue
occ_search(taxonKey=1, issue='DEPTH_UNLIKELY', fields =
   c('name','key','decimalLatitude','decimalLongitude','depth'))
## two issues
occ_search(taxonKey=1, issue=c('DEPTH_UNLIKELY','COORDINATE_ROUNDED'))
# Show all records in the Arizona State Lichen Collection that cant be matched to the GBIF
# backbone properly:
occ_search(datasetKey='84c0e1a0-f762-11e1-a439-00145eb45e9a',
   issue=c('TAXON_MATCH_NONE','TAXON_MATCH_HIGHERRANK'))

# Parsing output by issue
(res <- occ_search(geometry='POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))', limit = 50))
## what do issues mean, can print whole table, or search for matches
head(gbif_issues())
gbif_issues()[ gbif_issues()$code \%in\% c('cdround','cudc','gass84','txmathi'), ]
## or parse issues in various ways
### remove data rows with certain issue classes
library('magrittr')
res \%>\% occ_issues(gass84)
### split issues into separate columns
res \%>\% occ_issues(mutate = "split")
### expand issues to more descriptive names
res \%>\% occ_issues(mutate = "expand")
### split and expand
res \%>\% occ_issues(mutate = "split_expand")
### split, expand, and remove an issue class
res \%>\% occ_issues(-cudc, mutate = "split_expand")

# If you try multiple values for two different parameters you are wacked on the hand
# occ_search(taxonKey=c(2482598,2492010), recordedBy=c("smith","BJ Stacey"))

# Get a lot of data, here 1500 records for Helianthus annuus
# out <- occ_search(taxonKey=key, limit=1500)
# nrow(out$data)

# If you pass in an invalid polygon you get hopefully informative errors

### the WKT string is fine, but GBIF says bad polygon
wkt <- 'POLYGON((-178.59375 64.83258989321493,-165.9375 59.24622380205539,
-147.3046875 59.065977905449806,-130.78125 51.04484764446178,-125.859375 36.70806354647625,
-112.1484375 23.367471303759686,-105.1171875 16.093320185359257,-86.8359375 9.23767076398516,
-82.96875 2.9485268155066175,-82.6171875 -14.812060061226388,-74.8828125 -18.849111862023985,
-77.34375 -47.661687803329166,-84.375 -49.975955187343295,174.7265625 -50.649460483096114,
179.296875 -42.19189902447192,-176.8359375 -35.634976650677295,176.8359375 -31.835565983656227,
163.4765625 -6.528187613695323,152.578125 1.894796132058301,135.703125 4.702353722559447,
127.96875 15.077427674847987,127.96875 23.689804541429606,139.921875 32.06861069132688,
149.4140625 42.65416193033991,159.2578125 48.3160811030533,168.3984375 57.019804336633165,
178.2421875 59.95776046458139,-179.6484375 61.16708631440347,-178.59375 64.83258989321493))'

# occ_search(geometry = gsub("\n", '', wkt))

### unable to parse due to last number pair needing two numbers, not one
# wkt <- 'POLYGON((-178.5 64.8,-165.9 59.2,-147.3 59.0,-130.7 51.0,-125.8))'
# occ_search(geometry = wkt)

### unable to parse due to unclosed string
# wkt <- 'POLYGON((-178.5 64.8,-165.9 59.2,-147.3 59.0,-130.7 51.0))'
# occ_search(geometry = wkt)
### another of the same
# wkt <- 'POLYGON((-178.5 64.8,-165.9 59.2,-147.3 59.0,-130.7 51.0,-125.8 36.7))'
# occ_search(geometry = wkt)

### returns no results
# wkt <- 'LINESTRING(3 4,10 50,20 25)'
# occ_search(geometry = wkt)

### Apparently a point is allowed, but errors
# wkt <- 'POINT(45 -122)'
# occ_search(geometry = wkt)

## Faceting
x <- occ_search(facet = "country", limit = 0)
x$facets
x <- occ_search(facet = "establishmentMeans", limit = 10)
x$facets
x$data
x <- occ_search(facet = c("country", "basisOfRecord"), limit = 10)
x$data
x$facets
x$facets$country
x$facets$basisOfRecord
x$facets$basisOfRecord$count
x <- occ_search(facet = "country", facetMincount = 30000000L, limit = 10)
x$facets
x$data
# paging per each faceted variable
(x <- occ_search(
  facet = c("country", "basisOfRecord", "hasCoordinate"),
  country.facetLimit = 3,
  basisOfRecord.facetLimit = 6,
  limit = 0
))
x$facets


# You can set limit=0 to get number of results found
occ_search(datasetKey = '7b5d6a48-f762-11e1-a439-00145eb45e9a', limit = 0)$meta
occ_search(scientificName = 'Ursus americanus', limit = 0)$meta
occ_search(scientificName = 'Ursus americanus', limit = 0)$meta
}
}
\references{
https://www.gbif.org/developer/occurrence#search
}
\seealso{
\code{\link[=downloads]{downloads()}}, \code{\link[=occ_data]{occ_data()}}, \code{\link[=occ_facet]{occ_facet()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloads.R
\name{downloads}
\alias{downloads}
\title{Downloads interface}
\description{
GBIF provides two ways to get occurrence data: through the
\verb{/occurrence/search} route (see \code{\link[=occ_search]{occ_search()}}),
or via the \verb{/occurrence/download} route (many functions, see below).
\code{\link[=occ_search]{occ_search()}} is more appropriate for smaller data, while
\verb{occ_download*()} functions are more appropriate for larger data requests.
}
\section{Settings}{

You'll use \code{\link[=occ_download]{occ_download()}} to kick off a download. You'll need to
give that function settings from your GBIF profile: your user name, your
password, and your email. These three settings are required to use the
function. You can specify them in one of three ways:
\itemize{
\item Pass them to \code{occ_download} as parameters
\item Use R options: As options either in the current R session using
the \code{\link[=options]{options()}} function, or by setting them in your \code{.Rprofile} file, after
which point they'll be read in automatically
\item Use environment variables: As env vars either in the current R session using
the \code{\link[=Sys.setenv]{Sys.setenv()}} function, or by setting them in your
\code{.Renviron}/\code{.bash_profile} or similar files, after which point they'll be read
in automatically
}
}

\section{BEWARE}{

You can not perform that many downloads, so plan wisely.
See \emph{Rate limiting} below.
}

\section{Rate limiting}{

If you try to launch too many downloads, you will receive an 420
"Enhance Your Calm" response. If there is less then 100 in total
across all GBIF users, then you can have 3 running at a time. If
there are more than that, then each user is limited to 1 only.
These numbers are subject to change.
}

\section{Functions}{

\itemize{
\item \code{\link[=occ_download]{occ_download()}} - Start a download
\item \code{\link[=occ_download_prep]{occ_download_prep()}} - Prepare a download request
\item \code{\link[=occ_download_queue]{occ_download_queue()}} - Start many downloads in a queue
\item \code{\link[=occ_download_cached]{occ_download_cached()}} - Check for downloads already in your GBIF account
\item \code{\link[=occ_download_wait]{occ_download_wait()}} - Re-run \code{occ_download_meta()} until ready
\item \code{\link[=occ_download_meta]{occ_download_meta()}} - Get metadata progress on a single download
\item \code{\link[=occ_download_list]{occ_download_list()}} - List your downloads
\item \code{\link[=occ_download_cancel]{occ_download_cancel()}} - Cancel a download
\item \code{\link[=occ_download_cancel_staged]{occ_download_cancel_staged()}} - Cancels any jobs with status \code{RUNNING}
or \code{PREPARING}
\item \code{\link[=occ_download_get]{occ_download_get()}} - Retrieve a download
\item \code{\link[=occ_download_import]{occ_download_import()}} - Import a download from local file system
\item \code{\link[=occ_download_datasets]{occ_download_datasets()}} - List datasets for a download
\item \code{\link[=occ_download_dataset_activity]{occ_download_dataset_activity()}} - Lists the downloads activity
of a dataset
}

Download query composer methods:

See \link{download_predicate_dsl}
}

\section{Query length}{

GBIF has a limit of 12,000 characters for a download query. This means
that you can have a pretty long query, but at some point it may lead to an
error on GBIF's side and you'll have to split your query into a few.
}

\section{Download status}{

The following statuses can be found with any download:
\itemize{
\item PREPARING: just submitted by user and awaiting processing (typically only
a few seconds)
\item RUNNING: being created (takes typically 1-15 minutes)
\item FAILED: something unexpected went wrong
\item KILLED: user decided to abort the job while it was in PREPARING or RUNNING
phase
\item SUCCEEDED: The download was created and the user was informed
\item FILE_ERASED: The download was deleted according to the retention policy,
see https://www.gbif.org/faq?question=for-how-long-will-does-gbif-store-downloads
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgbif-package.r
\docType{data}
\name{isocodes}
\alias{isocodes}
\title{Table of country two character ISO codes, and GBIF names}
\description{
\itemize{
\item code. Two character ISO country code.
\item name. Name of country.
\item gbif_name. Name of country used by GBIF - this is the name
you want to use when searching by country in this package.
}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_lookup.r
\name{name_lookup}
\alias{name_lookup}
\title{Lookup names in all taxonomies in GBIF.}
\usage{
name_lookup(
  query = NULL,
  rank = NULL,
  higherTaxonKey = NULL,
  status = NULL,
  isExtinct = NULL,
  habitat = NULL,
  nameType = NULL,
  datasetKey = NULL,
  origin = NULL,
  nomenclaturalStatus = NULL,
  limit = 100,
  start = 0,
  facet = NULL,
  facetMincount = NULL,
  facetMultiselect = NULL,
  type = NULL,
  hl = NULL,
  issue = NULL,
  verbose = FALSE,
  return = NULL,
  curlopts = list()
)
}
\arguments{
\item{query}{Query term(s) for full text search.}

\item{rank}{CLASS, CULTIVAR, CULTIVAR_GROUP, DOMAIN, FAMILY, FORM, GENUS,
INFORMAL, INFRAGENERIC_NAME, INFRAORDER, INFRASPECIFIC_NAME,
INFRASUBSPECIFIC_NAME, KINGDOM, ORDER, PHYLUM, SECTION, SERIES, SPECIES,
STRAIN, SUBCLASS, SUBFAMILY, SUBFORM, SUBGENUS, SUBKINGDOM, SUBORDER,
SUBPHYLUM, SUBSECTION, SUBSERIES, SUBSPECIES, SUBTRIBE, SUBVARIETY,
SUPERCLASS, SUPERFAMILY, SUPERORDER, SUPERPHYLUM, SUPRAGENERIC_NAME,
TRIBE, UNRANKED, VARIETY}

\item{higherTaxonKey}{Filters by any of the higher Linnean rank keys. Note
this is within the respective checklist and not searching nub keys
across all checklists. This parameter accepts many inputs in a vector (
passed in the same request).}

\item{status}{Filters by the taxonomic status as one of:
\itemize{
\item ACCEPTED
\item DETERMINATION_SYNONYM Used for unknown child taxa referred to via
spec, ssp, ...
\item DOUBTFUL Treated as accepted, but doubtful whether this is correct.
\item HETEROTYPIC_SYNONYM More specific subclass of SYNONYM.
\item HOMOTYPIC_SYNONYM More specific subclass of SYNONYM.
\item INTERMEDIATE_RANK_SYNONYM Used in nub only.
\item MISAPPLIED More specific subclass of SYNONYM.
\item PROPARTE_SYNONYM More specific subclass of SYNONYM.
\item SYNONYM A general synonym, the exact type is unknown.
}}

\item{isExtinct}{(logical) Filters by extinction status (e.g.
\code{isExtinct=TRUE})}

\item{habitat}{(character) Filters by habitat. One of: marine, freshwater,
or terrestrial}

\item{nameType}{Filters by the name type as one of:
\itemize{
\item BLACKLISTED surely not a scientific name.
\item CANDIDATUS Candidatus is a component of the taxonomic name for a
bacterium that cannot be maintained in a Bacteriology Culture Collection.
\item CULTIVAR a cultivated plant name.
\item DOUBTFUL doubtful whether this is a scientific name at all.
\item HYBRID a hybrid formula (not a hybrid name).
\item INFORMAL a scientific name with some informal addition like "cf." or
indetermined like Abies spec.
\item SCINAME a scientific name which is not well formed.
\item VIRUS a virus name.
\item WELLFORMED a well formed scientific name according to present
nomenclatural rules.
}}

\item{datasetKey}{Filters by the dataset's key (a uuid)}

\item{origin}{(character) Filters by origin. One of:
\itemize{
\item SOURCE
\item DENORMED_CLASSIFICATION
\item VERBATIM_ACCEPTED
\item EX_AUTHOR_SYNONYM
\item AUTONYM
\item BASIONYM_PLACEHOLDER
\item MISSING_ACCEPTED
\item IMPLICIT_NAME
\item PROPARTE
\item VERBATIM_BASIONYM
}}

\item{nomenclaturalStatus}{Not yet implemented, but will eventually allow
for filtering by a nomenclatural status enum.}

\item{limit}{Number of records to return.
Hard maximum limit set by GBIF API: 99999.}

\item{start}{Record number to start at. Default: 0.}

\item{facet}{A vector/list of facet names used to retrieve the 100 most
frequent values for a field. Allowed facets are: datasetKey, higherTaxonKey,
rank, status, isExtinct, habitat, and nameType. Additionally threat and
nomenclaturalStatus are legal values but not yet implemented, so data will
not yet be returned for them.}

\item{facetMincount}{Used in combination with the facet parameter. Set
facetMincount={#} to exclude facets with a count less than {#}, e.g.
http://bit.ly/2osAUQB only shows the type values 'CHECKLIST' and 'OCCURRENCE'
because the other types have counts less than 10000}

\item{facetMultiselect}{(logical) Used in combination with the facet
parameter. Set \code{facetMultiselect=TRUE} to still return counts for
values that are not currently filtered, e.g. http://bit.ly/2JAymaC still
shows all type values even though type is being filtered
by \code{type=CHECKLIST}.}

\item{type}{Type of name. One of occurrence, checklist, or metadata.}

\item{hl}{(logical) Set \code{hl=TRUE} to highlight terms matching the query
when in fulltext search fields. The highlight will be an emphasis tag of
class \code{gbifH1} e.g. \code{query='plant', hl=TRUE}. Fulltext search
fields include: title, keyword, country, publishing country, publishing
organization title, hosting organization title, and description. One
additional full text field is searched which includes information from
metadata documents, but the text of this field is not returned in the
response.}

\item{issue}{Filters by issue. Issue has to be related to names. Type
\code{gbif_issues()} to get complete list of issues.}

\item{verbose}{(logical) If \code{TRUE}, all data is returned as a list for each
element. If \code{FALSE} (default) a subset of the data that is thought to be most
essential is organized into a data.frame.}

\item{return}{Defunct. All components are returned; index to the
one(s) you want}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
An object of class gbif, which is a S3 class list, with slots for
metadata (\code{meta}), the data itself (\code{data}), the taxonomic
hierarchy data (\code{hierarchies}), and vernacular names (\code{names}).
In addition, the object has attributes listing the user supplied arguments
and type of search, which is, differently from occurrence data, always
equals to 'single' even if multiple values for some parameters are given.
\code{meta} is a list of length four with offset, limit, endOfRecords and
count fields. \code{data} is a tibble (aka data.frame) containing all
information about the found taxa. \code{hierarchies} is a list of
data.frame's, one per GBIF key (taxon), containing its taxonomic
classification. Each data.frame contains two columns: \code{rankkey} and
\code{name}. \code{names} returns a list of data.frame's, one per GBIF key
(taxon), containing all vernacular names. Each data.frame contains two
columns: \code{vernacularName} and \code{language}.

A list of length five:
\itemize{
\item \strong{metadata}
\item \strong{data}: either a data.frame (\code{verbose=FALSE}, default) or a list (\code{verbose=TRUE}).
\item \strong{facets}
\item \strong{hierarchies}
\item \strong{names}
}
}
\description{
This service uses fuzzy lookup so that you can put in partial names and
you should get back those things that match. See examples below.

Faceting: If \code{facet=FALSE} or left to the default (NULL), no faceting
is done. And therefore, all parameters with facet in their name are
ignored (facetOnly, facetMincount, facetMultiselect).
}
\section{Repeat parameter inputs}{

Some parameters can take many inputs, and treated as 'OR' (e.g., a or b or
c). The following take many inputs:
\itemize{
\item \strong{rank}
\item \strong{higherTaxonKey}
\item \strong{status}
\item \strong{habitat}
\item \strong{nameType}
\item \strong{datasetKey}
\item \strong{origin}
}

see also \code{\link{many-values}}
}

\examples{
\dontrun{
# Look up names like mammalia
name_lookup(query='mammalia', limit = 20)

# Start with an offset
name_lookup(query='mammalia', limit=1)
name_lookup(query='mammalia', limit=1, start=2)

# large requests (paging is internally implemented).
# hard maximum limit set by GBIF API: 99999
# name_lookup(query = "Carnivora", limit = 10000)

# Get all data and parse it, removing descriptions which can be quite long
out <- name_lookup('Helianthus annuus', rank="species", verbose=TRUE)
lapply(out$data, function(x) {
  x[!names(x) \%in\% c("descriptions","descriptionsSerialized")]
})

# Search for a genus
name_lookup(query="Cnaemidophorus", rank="genus")
# Limit records to certain number
name_lookup('Helianthus annuus', rank="species", limit=2)

# Query by habitat
name_lookup(habitat = "terrestrial", limit=2)
name_lookup(habitat = "marine", limit=2)
name_lookup(habitat = "freshwater", limit=2)

# Using faceting
name_lookup(facet='status', limit=0, facetMincount='70000')
name_lookup(facet=c('status','higherTaxonKey'), limit=0,
  facetMincount='700000')

name_lookup(facet='nameType', limit=0)
name_lookup(facet='habitat', limit=0)
name_lookup(facet='datasetKey', limit=0)
name_lookup(facet='rank', limit=0)
name_lookup(facet='isExtinct', limit=0)

name_lookup(isExtinct=TRUE, limit=0)

# text highlighting
## turn on highlighting
res <- name_lookup(query='canada', hl=TRUE, limit=5)
res$data
name_lookup(query='canada', hl=TRUE, limit=45)
## and you can pass the output to gbif_names() function
res <- name_lookup(query='canada', hl=TRUE, limit=5)
gbif_names(res)

# Lookup by datasetKey (set up sufficient high limit, API maximum: 99999)
# name_lookup(datasetKey='3f8a1297-3259-4700-91fc-acc4170b27ce',
#   limit = 50000)

# Some parameters accept many inputs, treated as OR
name_lookup(rank = c("family", "genus"))
name_lookup(higherTaxonKey = c("119", "120", "121", "204"))
name_lookup(status = c("misapplied", "synonym"))$data
name_lookup(habitat = c("marine", "terrestrial"))
name_lookup(nameType = c("cultivar", "doubtful"))
name_lookup(datasetKey = c("73605f3a-af85-4ade-bbc5-522bfb90d847",
  "d7c60346-44b6-400d-ba27-8d3fbeffc8a5"))
name_lookup(datasetKey = "289244ee-e1c1-49aa-b2d7-d379391ce265",
  origin = c("SOURCE", "DENORMED_CLASSIFICATION"))

# Pass on curl options
name_lookup(query='Cnaemidophorus', rank="genus",
  curlopts = list(verbose = TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/species#searching}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_cancel.R
\name{occ_download_cancel}
\alias{occ_download_cancel}
\alias{occ_download_cancel_staged}
\title{Cancel a download creation process.}
\usage{
occ_download_cancel(key, user = NULL, pwd = NULL, curlopts = list())

occ_download_cancel_staged(
  user = NULL,
  pwd = NULL,
  limit = 20,
  start = 0,
  curlopts = list()
)
}
\arguments{
\item{key}{(character) A key generated from a request, like that from
\code{occ_download}. Required.}

\item{user}{(character) User name within GBIF's website. Required. See
Details.}

\item{pwd}{(character) User password within GBIF's website. Required. See
Details.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}

\item{limit}{Number of records to return. Default: 20}

\item{start}{Record number to start at. Default: 0}
}
\description{
Cancel a download creation process.
}
\details{
Note, these functions only cancel a job in progress. If your
download is already prepared for you, this won't do anything to change
that.

\code{occ_download_cancel} cancels a specific job by download key - returns
success message

\code{occ_download_cancel_staged} cancels all jobs with status \code{RUNNING}
or \code{PREPARING} - if none are found, returns a message saying so -
if some found, they are cancelled, returning message saying so
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\examples{
\dontrun{
# occ_download_cancel(key="0003984-140910143529206")
# occ_download_cancel_staged()
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gbifdata}
\alias{gbifdata}
\title{Get data.frame from occurrencelist, occurrencelist_many, or densitylist.}
\usage{
gbifdata(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{togeojson}
\alias{togeojson}
\title{Convert spatial data files to GeoJSON from various formats.}
\usage{
togeojson(...)
}
\description{
This function is defunct.  See the package togeojson for similar functionality.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{blanktheme}
\alias{blanktheme}
\title{Custom ggplot2 theme}
\usage{
blanktheme()
}
\description{
Custom ggplot2 theme
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{densitylist}
\alias{densitylist}
\title{The density web service provides access to records showing the density
of occurrence records from the GBIF Network by one-degree cell.}
\usage{
densitylist(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxrank.R
\name{taxrank}
\alias{taxrank}
\title{Get the possible values to be used for (taxonomic) rank arguments in GBIF
API methods.}
\usage{
taxrank()
}
\description{
Get the possible values to be used for (taxonomic) rank arguments in GBIF
API methods.
}
\examples{
\dontrun{
taxrank()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_get.r
\name{occ_get}
\alias{occ_get}
\alias{occ_get_verbatim}
\title{Get data for GBIF occurrences by occurrence key}
\usage{
occ_get(
  key,
  fields = "minimal",
  curlopts = list(),
  return = NULL,
  verbatim = NULL
)

occ_get_verbatim(key, fields = "minimal", curlopts = list())
}
\arguments{
\item{key}{(numeric/integer) one or more occurrence keys. required}

\item{fields}{(character) Default ("minimal") will return just taxon name,
key, latitude, and longitute. 'all' returns all fields. Or specify each
field you want returned by name, e.g. fields = c('name',
'decimalLatitude','altitude').}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}

\item{return}{Defunct. All components are returned now; index to the
one(s) you want}

\item{verbatim}{Defunct. verbatim records can now be retrieved using
\code{occ_get_verbatim()}}
}
\value{
For \code{occ_get} a list of lists. For \code{occ_get_verbatim} a data.frame
}
\description{
Get data for GBIF occurrences by occurrence key
}
\examples{
\dontrun{
occ_get(key=855998194)

# many occurrences
occ_get(key=c(101010, 240713150, 855998194))

# Verbatim data
occ_get_verbatim(key=855998194)
occ_get_verbatim(key=855998194, fields='all')
occ_get_verbatim(key=855998194,
 fields=c('scientificName', 'lastCrawled', 'county'))
occ_get_verbatim(key=c(855998194, 620594291))
occ_get_verbatim(key=c(855998194, 620594291), fields='all')
occ_get_verbatim(key=c(855998194, 620594291),
   fields=c('scientificName', 'decimalLatitude', 'basisOfRecord'))

# curl options, pass in a named list
occ_get(key=855998194, curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/occurrence#occurrence}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_backbone_checklist.R
\name{name_backbone_checklist}
\alias{name_backbone_checklist}
\title{Lookup names in the GBIF backbone taxonomy in a checklist.}
\usage{
name_backbone_checklist(name_data = NULL, verbose = FALSE, curlopts = list())
}
\arguments{
\item{name_data}{(data.frame or vector) see details.}

\item{verbose}{(logical) should the matching return non-exact matches}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A \code{data.frame} of matched names.
}
\description{
Lookup names in the GBIF backbone taxonomy in a checklist.
}
\details{
This function is a wrapper for  \code{name_backbone()}, which will work with
a list of names (a vector or a data.frame). The data.frame should have the
following column names, but \strong{only the 'name' column is required}. If only
one column is present, then that column is assumed to be the 'name' column.
\itemize{
\item \strong{name} : (required)
\item \strong{rank} : (optional)
\item \strong{kingdom} : (optional)
\item \strong{phylum} : (optional)
\item \strong{class} : (optional)
\item \strong{order} : (optional)
\item \strong{family} : (optional)
\item \strong{genus} : (optional)
}

The input columns will be returned as "verbatim_name","verbatim_rank",
"verbatim_phylum"...

The following aliases for the 'name' column will work (any case or with '_'
will work) :
\itemize{
\item "scientificName", "ScientificName", "scientific_name" ...
\item "sci_name", "sciname", "SCI_NAME" ...
\item "names", "NAMES" ...
\item "species", "SPECIES" ...
\item "species_name", "speciesname" ...
\item "sp_name", "SP_NAME", "spname" ...
}

If more than one aliases is present and no column is named 'name', then the
left-most column with an acceptable aliased name above is used.

This function can also be used with a character vector of names. In that case
no column names are needed of course.

This function is very similar to the GBIF species-lookup tool.
\href{https://www.gbif.org/tools/species-lookup}{https://www.gbif.org/tools/species-lookup}

If you have 1000s of names to match, it can take some minutes to get back all
of the matches. I have tested it with 30K names, but if you have more,
you might need a different approach. Also you will usually get better matches
if you include the author details.
}
\examples{
\dontrun{

library(rgbif)

name_data <- data.frame(
 scientificName = c(
   "Cirsium arvense (L.) Scop.", # a plant
   "Calopteryx splendens (Harris, 1780)", # an insect
   "Puma concolor (Linnaeus, 1771)", # a big cat
   "Ceylonosticta alwisi (Priyadarshana & Wijewardhane, 2016)", # newly discovered insect 
   "Puma concuolor (Linnaeus, 1771)", # a mis-spelled big cat
   "Fake species (John Waller 2021)", # a fake species
   "Calopteryx" # Just a Genus   
 ), description = c(
   "a plant",
   "an insect",
   "a big cat",
   "newly discovered insect",
   "a mis-spelled big cat",
   "a fake species",
   "just a GENUS"
 ), 
 kingdom = c(
   "Plantae",
   "Animalia",
   "Animalia",
   "Animalia",
   "Animalia",
   "Johnlia",
   "Animalia"
 ))

name_backbone_checklist(name_data)
name_backbone_checklist(name_data,verbose=TRUE) # return non-accepted names too 

# works with just vectors too 
name_list <- c(
"Cirsium arvense (L.) Scop.", 
"Calopteryx splendens (Harris, 1780)", 
"Puma concolor (Linnaeus, 1771)", 
"Ceylonosticta alwisi (Priyadarshana & Wijewardhane, 2016)", 
"Puma concuolor", 
"Fake species (John Waller 2021)", 
"Calopteryx")

name_backbone_checklist(name_list)
name_backbone_checklist(name_list,verbose=TRUE)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{occurrencedensity}
\alias{occurrencedensity}
\title{Returns summary counts of occurrence records by one-degree cell for a single
taxon, country, dataset, data publisher or data network.}
\usage{
occurrencedensity()
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_metadata.r
\name{occ_metadata}
\alias{occ_metadata}
\title{Search for catalog numbers, collection codes, collector names, and
institution codes.}
\usage{
occ_metadata(
  type = "catalogNumber",
  q = NULL,
  limit = 5,
  pretty = TRUE,
  curlopts = list()
)
}
\arguments{
\item{type}{Type of data, one of catalogNumber, collectionCode, recordedBy,
or institutionCode. Unique partial strings work too, like 'cat' for
catalogNumber}

\item{q}{Search term}

\item{limit}{Number of results, default=5}

\item{pretty}{Pretty as true (Default) uses cat to print data, \code{FALSE} gives
character strings.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\description{
Search for catalog numbers, collection codes, collector names, and
institution codes.
}
\examples{
\dontrun{
# catalog number
occ_metadata(type = "catalogNumber", q=122)

# collection code
occ_metadata(type = "collectionCode", q=12)

# institution code
occ_metadata(type = "institutionCode", q='GB')

# recorded by
occ_metadata(type = "recordedBy", q='scott')

# data as character strings
occ_metadata(type = "catalogNumber", q=122, pretty=FALSE)

# Change number of results returned
occ_metadata(type = "catalogNumber", q=122, limit=10)

# Partial unique type strings work too
occ_metadata(type = "cat", q=122)

# Pass on curl options
occ_metadata(type = "cat", q=122, curlopts = list(verbose = TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/occurrence#search}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_suggest.r
\name{name_suggest}
\alias{name_suggest}
\title{Suggest up to 20 name usages.}
\usage{
name_suggest(
  q = NULL,
  datasetKey = NULL,
  rank = NULL,
  fields = NULL,
  start = NULL,
  limit = 100,
  curlopts = list()
)
}
\arguments{
\item{q}{(character, required) Simple search parameter. The value for
this parameter can be a simple word or a phrase. Wildcards can be added to
the simple word parameters only, e.g. q=\emph{puma}}

\item{datasetKey}{(character) Filters by the checklist dataset key (a uuid,
see examples)}

\item{rank}{(character) A taxonomic rank. One of class, cultivar,
cultivar_group, domain, family, form, genus, informal, infrageneric_name,
infraorder, infraspecific_name, infrasubspecific_name, kingdom, order,
phylum, section, series, species, strain, subclass, subfamily, subform,
subgenus, subkingdom, suborder, subphylum, subsection, subseries,
subspecies, subtribe, subvariety, superclass, superfamily, superorder,
superphylum, suprageneric_name, tribe, unranked, or variety.}

\item{fields}{(character) Fields to return in output data.frame (simply
prunes columns off)}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A list, with two elements \code{data} (tibble) and \code{hierarchy} (list of
data.frame's). If 'higherClassificationMap' is one of the \code{fields} requested,
then \code{hierarchy} is a list of data.frame's; if not included, \code{hierarchy}
is an empty list.
}
\description{
A quick and simple autocomplete service that returns up to 20 name
usages by doing prefix matching against the scientific name. Results
are ordered by relevance.
}
\section{Repeat parmeter inputs}{

Some parameters can take many inputs, and treated as 'OR' (e.g., a or b or
c). The following take many inputs:
\itemize{
\item \strong{rank}
\item \strong{datasetKey}
}

see also \link{many-values}
}

\examples{
\dontrun{
name_suggest(q='Puma concolor')
name_suggest(q='Puma')
name_suggest(q='Puma', rank="genus")
name_suggest(q='Puma', rank="subspecies")
name_suggest(q='Puma', rank="species")
name_suggest(q='Puma', rank="infraspecific_name")

name_suggest(q='Puma', limit=2)
name_suggest(q='Puma', fields=c('key','canonicalName'))
name_suggest(q='Puma', fields=c('key','canonicalName',
  'higherClassificationMap'))

# Some parameters accept many inputs, treated as OR
name_suggest(rank = c("family", "genus"))
name_suggest(datasetKey = c("73605f3a-af85-4ade-bbc5-522bfb90d847",
  "d7c60346-44b6-400d-ba27-8d3fbeffc8a5"))

# If 'higherClassificationMap' in fields, a list is returned
name_suggest(q='Puma', fields=c('key','higherClassificationMap'))

# Pass on curl options
name_suggest(q='Puma', limit=200, curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/species#searching}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_fetch.R
\name{map_fetch}
\alias{map_fetch}
\title{Fetch aggregated density maps of GBIF occurrences}
\usage{
map_fetch(
  source = "density",
  x = 0,
  y = 0,
  z = 0,
  format = "@1x.png",
  srs = "EPSG:4326",
  bin = NULL,
  hexPerTile = NULL,
  squareSize = NULL,
  style = "classic.point",
  taxonKey = NULL,
  datasetKey = NULL,
  country = NULL,
  publishingOrg = NULL,
  publishingCountry = NULL,
  year = NULL,
  basisOfRecord = NULL,
  ...
)
}
\arguments{
\item{source}{(character) Either \code{density} for fast, precalculated tiles,
or \code{adhoc} for any search. Default: \code{density}}

\item{x}{(integer) the column. Default: 0}

\item{y}{(integer) the row. Default: 0}

\item{z}{(integer) the zoom. Default: 0}

\item{format}{(character) The data format, one of:
\itemize{
\item \verb{@Hx.png} for a 256px raster tile
\item \verb{@1x.png} for a 512px raster tile (the default)
\item \verb{@2x.png} for a 1024px raster tile
\item \verb{@3x.png} for a 2048px raster tile
\item \verb{@4x.png} for a 4096px raster tile
}}

\item{srs}{(character) Spatial reference system. One of:
\itemize{
\item \code{EPSG:3857} (Web Mercator)
\item \code{EPSG:4326} (WGS84 plate care?)
\item \code{EPSG:3575} (Arctic LAEA on 10 degrees E)
\item \code{EPSG:3031} (Antarctic stereographic)
}}

\item{bin}{(character) \code{square} or \code{hex} to aggregate occurrence counts into
squares or hexagons. Points by default. optional}

\item{hexPerTile}{(integer) sets the size of the hexagons
(the number horizontally across a tile). optional}

\item{squareSize}{(integer) sets the size of the squares. Choose a factor
of 4096 so they tessalate correctly: probably from 8, 16, 32, 64, 128,
256, 512. optional}

\item{style}{(character) for raster tiles, choose from the available styles.
Defaults to classic.point. optional. THESE DON'T WORK YET.}

\item{taxonKey}{(integer/numeric/character) search by taxon key, can only
supply 1. optional}

\item{datasetKey}{(character) search by taxon key, can only supply 1.
optional}

\item{country}{(character) search by taxon key, can only supply 1.
optional}

\item{publishingOrg}{(character) search by taxon key, can only supply 1.
optional}

\item{publishingCountry}{(character) search by taxon key, can only
supply 1. optional}

\item{year}{(integer) integer that limits the search to a certain year or,
if passing a vector of integers, multiple years, for example
\code{1984} or \code{c(2016, 2017, 2018)} or \code{2010:2015} (years 2010 to 2015). optional}

\item{basisOfRecord}{(character) one or more basis of record states to
include records with that basis of record. The full list is: \code{c("OBSERVATION", "HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "MATERIAL_SAMPLE", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN", "LIVING_SPECIMEN", "LITERATURE", "UNKNOWN")}. optional}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
an object of class \code{RasterLayer}
}
\description{
This function is a wrapper for the GBIF mapping api version 2.0.
The mapping API is a web map tile service making it straightforward to
visualize GBIF content on interactive maps, and overlay content from other
sources. It returns tile maps with number of
GBIF records per area unit that can be used in a variety of ways, for example
in interactive leaflet web maps. Map details are specified by a number of
query parameters, some of them optional. Full documentation of the GBIF
mapping api can be found at https://www.gbif.org/developer/maps
}
\details{
This function uses the arguments passed on to generate a query
to the GBIF web map API. The API returns a web tile object as png that is
read and converted into an R raster object. The break values or nbreaks
generate a custom colour palette for the web tile, with each bin
corresponding to one grey value. After retrieval, the raster is reclassified
to the actual break values. This is a somewhat hacky but nonetheless
functional solution in the absence of a GBIF raster API implementation.

We add extent and set the projection for the output. You can reproject
after retrieving the output.
}
\note{
Styles don't work yet, sorry, we'll try to fix it asap.
}
\examples{
\dontrun{
if (
 requireNamespace("png", quietly = TRUE) &&
 requireNamespace("raster", quietly = TRUE)
) {
  x <- map_fetch(taxonKey = 2480498, year = 2007:2011)
  x
  # gives a RasterLayer object
  class(x)
  # visualize
  library(raster)
  plot(x)

  # different srs
  ## 3857
  y <- map_fetch(taxonKey = 2480498, year = 2010, srs = "EPSG:3857")
  plot(y)
  ## 3031
  z <- map_fetch(taxonKey = 2480498, year = 2010, srs = "EPSG:3031", verbose = TRUE)
  plot(z)
  # 3575
  z <- map_fetch(taxonKey = 2480498, year = 2010, srs = "EPSG:3575")
  plot(z)

  # bin
  plot(map_fetch(taxonKey = 212, year = 1998, bin = "hex",
     hexPerTile = 30, style = "classic-noborder.poly"))

  # styles
  plot(map_fetch(taxonKey = 2480498, style = "purpleYellow.point"))

  # query with basisOfRecord
  map_fetch(taxonKey = 2480498, year = 2010,
    basisOfRecord = "HUMAN_OBSERVATION")
  map_fetch(taxonKey = 2480498, year = 2010,
    basisOfRecord = c("HUMAN_OBSERVATION", "LIVING_SPECIMEN"))
 }
}
}
\references{
https://www.gbif.org/developer/maps
}
\seealso{
\code{\link[=mvt_fetch]{mvt_fetch()}}
}
\author{
Laurens Geffert \email{laurensgeffert@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/many-values.R
\name{many-values}
\alias{many-values}
\title{Many value inputs to some parameters}
\description{
Many value inputs to some parameters
}
\details{
There are some differences in how functions across \pkg{rgbif} behave
with respect to many values given to a single parameter (let's call it
\code{foo}).

The following functions originally only iterated over many
values passed to \code{foo} as a vector (e.g., \code{foo = c(1, 2)}) with completely
separate HTTP requests. But now these functions also support passing in many
values to the same HTTP request (e.g., \code{foo = "1;2"}). This is a bit
awkward, but means that we don't break existing code. See
"Multiple values passed to a parameter" in \code{occ_search}/\code{occ_data} for more
information.

\itemize{
\item \code{\link[=occ_search]{occ_search()}}
\item \code{\link[=occ_data]{occ_data()}}
}

The following functions, unlike those above, only support passing in many
values to the same HTTP request, which is done like \code{foo = c("1", "2")}.

\itemize{
\item \code{\link[=dataset_search]{dataset_search()}}
\item \code{\link[=dataset_suggest]{dataset_suggest()}}
\item \code{\link[=name_lookup]{name_lookup()}}
\item \code{\link[=name_suggest]{name_suggest()}}
\item \code{\link[=name_usage]{name_usage()}}
}

Last, some parameters in the functions above don't accept more than one,
and some functions don't have any parameters that accept more than one
value (i.e., none of those listed above).

Each function that has at least some parameters that accept many values
also has documentation on this issue.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wkt_parse.R
\name{wkt_parse}
\alias{wkt_parse}
\title{parse wkt into smaller bits}
\usage{
wkt_parse(wkt, geom_big, geom_size = 40, geom_n = 10)
}
\arguments{
\item{wkt}{(character) A WKT string. Required.}

\item{geom_big}{(character) One of "axe" or "bbox". Required.}

\item{geom_size}{(integer) An integer indicating size of the cell.
Default: 40.}

\item{geom_n}{(integer) An integer indicating number of cells in
each dimension. Default: 10.}
}
\description{
parse wkt into smaller bits
}
\examples{
wkt <- "POLYGON((13.26349675655365 52.53991761181831,18.36115300655365 54.11445544219924,
21.87677800655365 53.80418956368524,24.68927800655365 54.217364774722455,28.20490300655365
54.320018299365124,30.49005925655365 52.85948216284084,34.70880925655365 52.753220564427814,
35.93927800655365 50.46131871049754,39.63068425655365 49.55761261299145,40.86115300655365
46.381388009130845,34.00568425655365 45.279102926537,33.30255925655365 48.636868465271846,
30.13849675655365 49.78513301801265,28.38068425655365 47.2236377039631,29.78693425655365
44.6572866068524,27.67755925655365 42.62220075124676,23.10724675655365 43.77542058000212,
24.51349675655365 47.10412345120368,26.79865300655365 49.55761261299145,23.98615300655365
52.00209943876426,23.63459050655365 49.44345313705238,19.41584050655365 47.580567827212114,
19.59162175655365 44.90682206053508,20.11896550655365 42.36297154876359,22.93146550655365
40.651849782081555,25.56818425655365 39.98171166226459,29.61115300655365 40.78507856230178,
32.95099675655365 40.38459278067577,32.95099675655365 37.37491910393631,26.27130925655365
33.65619609886799,22.05255925655365 36.814081996401605,18.71271550655365 36.1072176729021,
18.53693425655365 39.16878677351903,15.37287175655365 38.346355762190846,15.19709050655365
41.578843777436326,12.56037175655365 41.050735748143424,12.56037175655365 44.02872991212046,
15.19709050655365 45.52594200494078,16.42755925655365 48.05271546733352,17.48224675655365
48.86865641518059,10.62677800655365 47.817178329053135,9.57209050655365 44.154980365192,
8.16584050655365 40.51835445724746,6.05646550655365 36.53210972067291,0.9588092565536499
31.583640057148145,-5.54509699344635 35.68001485298146,-6.77556574344635 40.51835445724746,
-9.41228449344635 38.346355762190846,-12.40056574344635 35.10683619158607,-15.74040949344635
38.07010978950028,-14.68572199344635 41.31532459432774,-11.69744074344635 43.64836179231387,
-8.88494074344635 42.88035509418534,-4.31462824344635 43.52103366008421,-8.35759699344635
47.2236377039631,-8.18181574344635 50.12441989397795,-5.01775324344635 49.55761261299145,
-2.73259699344635 46.25998980446569,-1.67790949344635 44.154980365192,-1.32634699344635
39.30493590580802,2.18927800655365 41.44721797271696,4.47443425655365 43.26556960420879,
2.18927800655365 46.7439668697322,1.83771550655365 50.3492841273576,6.93537175655365
49.671505849335254,5.00177800655365 52.32557322466785,7.81427800655365 51.67627099802223,
7.81427800655365 54.5245591562317,10.97834050655365 51.89375191441792,10.97834050655365
55.43241335888528,13.26349675655365 52.53991761181831))"
wkt <- gsub("\n", " ", wkt)

if (requireNamespace("sf", quietly=TRUE)) {
# to a bounding box in wkt format
wkt_parse(wkt, geom_big = "bbox")

# to many wkt strings, chopped up from input
wkt_parse(wkt, geom_big = "axe")
wkt_parse(wkt, geom_big = "axe", 60)
wkt_parse(wkt, geom_big = "axe", 30)
wkt_parse(wkt, geom_big = "axe", 20)
wkt_parse(wkt, geom_big = "axe", 10)
wkt_parse(wkt, geom_big = "axe", 5)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_queue.R
\name{GbifQueue}
\alias{GbifQueue}
\title{GbifQueue}
\description{
GBIF download queue
}
\examples{
\dontrun{
if (interactive()) { # dont run in automated example runs, too costly
x <- GbifQueue$new(
  occ_download(pred('taxonKey', 3119195), pred("year", 1976)),
  occ_download(pred('taxonKey', 3119195), pred("year", 2001)),
  occ_download(pred('taxonKey', 3119195), pred("year", 2001), pred_lte("month", 8)),
  occ_download(pred('taxonKey', 3119195), pred("year", 2004)),
  occ_download(pred('taxonKey', 3119195), pred("year", 2005))
)
x
x$reqs
x$add_all()
x
x$jobs()
x
x$remove(x$reqs[[1]])
x
x$reqs[[1]]$run()
x$reqs[[1]]$result

# pre-prepared download request
z <- occ_download_prep(
  pred_in("basisOfRecord", c("HUMAN_OBSERVATION","OBSERVATION")),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("year", 1993),
  user = "foo", pwd = "bar", email = "foo@bar.com"
)
out <- GbifQueue$new(.list = list(z))
out
out$reqs
}}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{reqs}}{(list) a named list of objects of class \code{\link[=occ_download]{occ_download()}}}

\item{\code{queue}}{(list) holds the queued jobs}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{GbifQueue$print()}}
\item \href{#method-new}{\code{GbifQueue$new()}}
\item \href{#method-add}{\code{GbifQueue$add()}}
\item \href{#method-add_all}{\code{GbifQueue$add_all()}}
\item \href{#method-remove}{\code{GbifQueue$remove()}}
\item \href{#method-next_}{\code{GbifQueue$next_()}}
\item \href{#method-last_}{\code{GbifQueue$last_()}}
\item \href{#method-jobs}{\code{GbifQueue$jobs()}}
\item \href{#method-clone}{\code{GbifQueue$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{GbifQueue} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{GbifQueue} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$new(..., .list = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{any number of \code{\link[=occ_download]{occ_download()}} requests}

\item{\code{.list}}{any number of \code{\link[=occ_download]{occ_download()}} requests as \code{lazy}
objects, called with e.g., \code{lazyeval::lazy()}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{GbifQueue} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add"></a>}}
\if{latex}{\out{\hypertarget{method-add}{}}}
\subsection{Method \code{add()}}{
Add single jobs to the queue
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$add(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{an \code{\link[=occ_download]{occ_download()}} object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; adds job (\code{x}) to the queue
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add_all"></a>}}
\if{latex}{\out{\hypertarget{method-add_all}{}}}
\subsection{Method \code{add_all()}}{
Add all jobs to the queue
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$add_all()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remove"></a>}}
\if{latex}{\out{\hypertarget{method-remove}{}}}
\subsection{Method \code{remove()}}{
Remove a job from the queue
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$remove(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{an \code{\link[=occ_download]{occ_download()}} object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-next_"></a>}}
\if{latex}{\out{\hypertarget{method-next_}{}}}
\subsection{Method \code{next_()}}{
Get the next job in the \code{queue}. if no more jobs,
returns empty list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$next_()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
next job or empty list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-last_"></a>}}
\if{latex}{\out{\hypertarget{method-last_}{}}}
\subsection{Method \code{last_()}}{
Get the last job in the \code{queue}. if no more jobs,
returns empty list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$last_()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
last job or empty list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-jobs"></a>}}
\if{latex}{\out{\hypertarget{method-jobs}{}}}
\subsection{Method \code{jobs()}}{
Get number of jobs in the \code{queue}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$jobs()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(integer) number of jobs
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GbifQueue$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{taxonsearch}
\alias{taxonsearch}
\title{Search for taxa in GBIF.}
\usage{
taxonsearch(...)
}
\description{
This function is defunct.
}
\seealso{
occ_search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{providers}
\alias{providers}
\title{Get data providers and their unique keys.}
\usage{
providers(...)
}
\description{
This function is defunct.
}
\seealso{
networks organizations datasets
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_spellcheck.R
\name{occ_spellcheck}
\alias{occ_spellcheck}
\title{Spell check search term for occurrence searches}
\usage{
occ_spellcheck(...)
}
\arguments{
\item{...}{ignored}
}
\description{
Spell check search term for occurrence searches
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_dataset_activity.R
\name{occ_download_dataset_activity}
\alias{occ_download_dataset_activity}
\title{Lists the downloads activity of a dataset}
\usage{
occ_download_dataset_activity(
  dataset,
  limit = 20,
  start = 0,
  curlopts = list()
)
}
\arguments{
\item{dataset}{(character) A dataset key}

\item{limit}{(integer/numeric) Number of records to return. Default: 20,
Max: 1000}

\item{start}{(integer/numeric) Record number to start at. Default: 0}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
a list with two slots:
\itemize{
\item meta: a single row data.frame with columns: \code{offset}, \code{limit},
\code{endofrecords}, \code{count}
\item results: a tibble with the nested data flattened, with many
columns with the same \code{download.} or \code{download.request.} prefixes
}
}
\description{
Lists the downloads activity of a dataset
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\examples{
\dontrun{
res <- occ_download_dataset_activity("7f2edc10-f762-11e1-a439-00145eb45e9a")
res
res$meta
res$meta$count

# pagination
occ_download_dataset_activity("7f2edc10-f762-11e1-a439-00145eb45e9a",
limit = 3000)
occ_download_dataset_activity("7f2edc10-f762-11e1-a439-00145eb45e9a",
limit = 3, start = 10)
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_facet.R
\name{occ_facet}
\alias{occ_facet}
\title{Facet GBIF occurrences}
\usage{
occ_facet(facet, facetMincount = NULL, curlopts = list(), ...)
}
\arguments{
\item{facet}{(character) a character vector of length 1 or greater. Required.}

\item{facetMincount}{(numeric) minimum number of records to be included
in the faceting results}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}

\item{...}{Facet parameters, such as for paging based on each facet
variable, e.g., \code{country.facetLimit}}
}
\value{
A list of tibbles (data.frame's) for each facet (each element of
the facet parameter).
}
\description{
Facet GBIF occurrences
}
\details{
All fields can be faceted on except for last "lastInterpreted",
"eventDate", and "geometry"

If a faceted variable is not found, it is silently dropped, returning
nothing for that query
}
\examples{
\dontrun{
occ_facet(facet = "country")

# facetMincount - minimum number of records to be included
#   in the faceting results
occ_facet(facet = "country", facetMincount = 30000000L)
occ_facet(facet = c("country", "basisOfRecord"))

# paging with many facets
occ_facet(
  facet = c("country", "basisOfRecord", "hasCoordinate"),
  country.facetLimit = 3,
  basisOfRecord.facetLimit = 6
)

# paging
## limit
occ_facet(facet = "country", country.facetLimit = 3)
## offset
occ_facet(facet = "country", country.facetLimit = 3,
  country.facetOffset = 3)

# Pass on curl options
occ_facet(facet = "country", country.facetLimit = 3,
  curlopts = list(verbose = TRUE))
}
}
\seealso{
\code{\link[=occ_search]{occ_search()}} also has faceting ability, but
can include occurrence data in addition to facets
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_count.r
\name{occ_count}
\alias{occ_count}
\title{Get number of occurrence records.}
\usage{
occ_count(
  taxonKey = NULL,
  georeferenced = NULL,
  basisOfRecord = NULL,
  datasetKey = NULL,
  date = NULL,
  typeStatus = NULL,
  country = NULL,
  year = NULL,
  from = 2000,
  to = 2012,
  type = "count",
  publishingCountry = "US",
  protocol = NULL,
  curlopts = list()
)
}
\arguments{
\item{taxonKey}{Species key}

\item{georeferenced}{Return only occurrence records with lat/long data
(\code{TRUE}) or those that don't have that data (\code{FALSE}, default). Note that
you can also get record count with \code{\link[=occ_search]{occ_search()}} by setting \code{limit=0}}

\item{basisOfRecord}{Basis of record}

\item{datasetKey}{Dataset key}

\item{date}{Collection date}

\item{typeStatus}{A type status. See \code{\link[=typestatus]{typestatus()}} dataset for
options}

\item{country}{Country data was collected in, two letter abbreviation. See
https://countrycode.org/ for abbreviations.}

\item{year}{Year data were collected in}

\item{from}{Year to start at}

\item{to}{Year to end at}

\item{type}{One of count (default), schema, basisOfRecord, countries, or
year.}

\item{publishingCountry}{Publishing country, two letter ISO country code}

\item{protocol}{Protocol. E.g., 'DWC_ARCHIVE'}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A single numeric value, or a list of numerics.
}
\description{
Get number of occurrence records.
}
\details{
There is a slight difference in the way records are counted here vs.
results from \code{\link[=occ_search]{occ_search()}}. For equivalent outcomes, in the
\code{\link[=occ_search]{occ_search()}} function use \code{hasCoordinate=TRUE}, and
\code{hasGeospatialIssue=FALSE} to have the same outcome for this function
using \code{georeferenced=TRUE}.
}
\section{Supported dimensions}{

That is, there are only a certain set of supported query parameter
combinations that GBIF allows on this API route. They can be found with the
call \code{occ_count(type='schema')}. They are also presented below:
\itemize{
\item basisOfRecord
\item basisOfRecord, country
\item basisOfRecord, country, isGeoreferenced
\item basisOfRecord, country, isGeoreferenced, taxonKey
\item basisOfRecord, country, taxonKey
\item basisOfRecord, datasetKey
\item basisOfRecord, datasetKey, isGeoreferenced
\item basisOfRecord, datasetKey, isGeoreferenced, taxonKey
\item basisOfRecord, datasetKey, taxonKey
\item basisOfRecord, isGeoreferenced, taxonKey
\item basisOfRecord, isGeoreferenced, publishingCountry
\item basisOfRecord, isGeoreferenced, publishingCountry, taxonKey
\item basisOfRecord, publishingCountry
\item basisOfRecord, publishingCountry, taxonKey
\item basisOfRecord, taxonKey
\item country
\item country, datasetKey, isGeoreferenced
\item country, isGeoreferenced
\item country, isGeoreferenced, publishingCountry
\item country, isGeoreferenced, taxonKey
\item country, publishingCountry
\item country, taxonKey
\item country, typeStatus
\item datasetKey
\item datasetKey, isGeoreferenced
\item datasetKey, isGeoreferenced, taxonKey
\item datasetKey, issue
\item datasetKey, taxonKey
\item datasetKey, typeStatus
\item isGeoreferenced
\item isGeoreferenced, publishingCountry
\item isGeoreferenced, publishingCountry, taxonKey
\item isGeoreferenced, taxonKey
\item issue
\item publishingCountry
\item publishingCountry, taxonKey
\item publishingCountry, typeStatus
\item taxonKey
\item taxonKey, typeStatus
\item typeStatus
\item protocol
\item year
}
}

\examples{
\dontrun{
occ_count(basisOfRecord='OBSERVATION')
occ_count(georeferenced=TRUE)
occ_count(country='DE')
occ_count(country='CA', georeferenced=TRUE, basisOfRecord='OBSERVATION')
occ_count(datasetKey='9e7ea106-0bf8-4087-bb61-dfe4f29e0f17')
occ_count(year=2012)
occ_count(taxonKey=2435099)
occ_count(taxonKey=2435099, georeferenced=TRUE)

# Just schema
occ_count(type='schema')

# Counts by basisOfRecord types
occ_count(type='basisOfRecord')

# Counts by basisOfRecord types and taxonkey
occ_count(taxonKey=2435099, basisOfRecord='OBSERVATION')

# Counts by typeStatus
occ_count(typeStatus='ALLOTYPE')
occ_count(typeStatus='HOLOTYPE')

# Counts by countries. publishingCountry must be supplied (default to US)
occ_count(type='countries')

# Counts by year. from and to years have to be supplied, default to 2000
# and 2012
occ_count(type='year', from=2000, to=2012)

# Counts by publishingCountry, must supply a country (default to US)
occ_count(type='publishingCountry')
occ_count(type='publishingCountry', country='BZ')

# Pass on curl options
occ_count(type='year', from=2000, to=2012, curlopts = list(verbose = TRUE))
}
}
\references{
https://www.gbif.org/developer/occurrence#metrics
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_backbone.r
\name{name_backbone}
\alias{name_backbone}
\alias{name_backbone_verbose}
\title{Lookup names in the GBIF backbone taxonomy.}
\usage{
name_backbone(
  name,
  rank = NULL,
  kingdom = NULL,
  phylum = NULL,
  class = NULL,
  order = NULL,
  family = NULL,
  genus = NULL,
  strict = FALSE,
  verbose = FALSE,
  start = NULL,
  limit = 100,
  curlopts = list()
)

name_backbone_verbose(
  name,
  rank = NULL,
  kingdom = NULL,
  phylum = NULL,
  class = NULL,
  order = NULL,
  family = NULL,
  genus = NULL,
  strict = FALSE,
  start = NULL,
  limit = 100,
  curlopts = list()
)
}
\arguments{
\item{name}{(character) Full scientific name potentially with authorship
(required)}

\item{rank}{(character) The rank given as our rank enum. (optional)}

\item{kingdom}{(character) If provided default matching will also try to
match against this if no direct match is found for the name alone.
(optional)}

\item{phylum}{(character) If provided default matching will also try to
match against this if no direct match is found for the name alone.
(optional)}

\item{class}{(character) If provided default matching will also try to
match against this if no direct match is found for the name alone.
(optional)}

\item{order}{(character) If provided default matching will also try to
match against this if no direct match is found for the name alone.
(optional)}

\item{family}{(character) If provided default matching will also try to
match against this if no direct match is found for the name alone.
(optional)}

\item{genus}{(character) If provided default matching will also try to
match against this if no direct match is found for the name alone.
(optional)}

\item{strict}{(logical) If \code{TRUE} it (fuzzy) matches only the given name,
but never a taxon in the upper classification (optional)}

\item{verbose}{(logical) should the function give back more results.
See function \code{name_backbone_verbose()}}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
For \code{name_backbone}, a data.frame for a single taxon with many
columns. For \code{name_backbone_verbose}, a larger number of results in a
data.frame the results of resulting from fuzzy matching.
You will also get back your input name, rank, kingdom, phylum ect. as
columns input_name, input_rank, input_kingdom ect. so you can check the
results.
}
\description{
Lookup names in the GBIF backbone taxonomy.
}
\details{
If you don't get a match, GBIF gives back a data.frame with columns
\code{synonym}, \code{confidence}, and \code{matchType='NONE'}.
}
\examples{
\dontrun{
name_backbone(name='Helianthus annuus', kingdom='plants')
name_backbone(name='Helianthus', rank='genus', kingdom='plants')
name_backbone(name='Poa', rank='genus', family='Poaceae')

# Verbose - gives back alternatives
## Strictness
name_backbone_verbose(name='Poa', kingdom='plants',
  strict=FALSE)
name_backbone_verbose(name='Helianthus annuus', kingdom='plants',
  strict=TRUE)

# Non-existent name - returns list of lenght 3 stating no match
name_backbone(name='Aso')
name_backbone(name='Oenante')

# Pass on curl options
name_backbone(name='Oenante', curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/species#searching}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{occurrenceget}
\alias{occurrenceget}
\title{Get individual records for a given occurrence record.}
\usage{
occurrenceget(...)
}
\description{
This function is defunct.
}
\seealso{
occ_get
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsenames.R
\name{parsenames}
\alias{parsenames}
\title{Parse taxon names using the GBIF name parser.}
\usage{
parsenames(scientificname, curlopts = list())
}
\arguments{
\item{scientificname}{A character vector of scientific names.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A \code{data.frame} containing fields extracted from parsed
taxon names. Fields returned are the union of fields extracted from
all species names in \code{scientificname}.
}
\description{
Parse taxon names using the GBIF name parser.
}
\examples{
\dontrun{
parsenames(scientificname='x Agropogon littoralis')
parsenames(c('Arrhenatherum elatius var. elatius',
             'Secale cereale subsp. cereale', 'Secale cereale ssp. cereale',
             'Vanessa atalanta (Linnaeus, 1758)'))
parsenames("Ajuga pyramidata")
parsenames("Ajuga pyramidata x reptans")

# Pass on curl options
# res <- parsenames(c('Arrhenatherum elatius var. elatius',
#          'Secale cereale subsp. cereale', 'Secale cereale ssp. cereale',
#          'Vanessa atalanta (Linnaeus, 1758)'), curlopts=list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/species#parser}
}
\author{
John Baumgartner (johnbb@student.unimelb.edu.au)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gbifmap}
\alias{gbifmap}
\title{Get Github credentials from use in console}
\usage{
gbifmap(...)
}
\description{
This function is defunct.  See the package gistr for similar functionality.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{occurrencelist}
\alias{occurrencelist}
\title{Occurrencelist searches for taxon concept records matching a range of filters.}
\usage{
occurrencelist(...)
}
\description{
This function is defunct.
}
\seealso{
occ_search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enumeration.R
\name{enumeration}
\alias{enumeration}
\alias{enumeration_country}
\title{Enumerations.}
\usage{
enumeration(x = NULL, curlopts = list())

enumeration_country(curlopts = list())
}
\arguments{
\item{x}{A given enumeration.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
\code{enumeration} returns a character vector, while
\code{enumeration_country} returns a data.frame.
}
\description{
Many parts of the GBIF API make use of enumerations, i.e. controlled
vocabularies for specific topics - and are available via these functions
}
\examples{
\dontrun{
# basic enumeration
enumeration()
enumeration("NameType")
enumeration("MetadataType")
enumeration("TypeStatus")

# country enumeration
enumeration_country()

# curl options
enumeration(curlopts = list(verbose=TRUE))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{density_spplist}
\alias{density_spplist}
\title{The density web service provides access to records showing the density
of occurrence records from the GBIF Network by one-degree cell.}
\usage{
density_spplist(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/installations.r
\name{installations}
\alias{installations}
\title{Installations metadata.}
\usage{
installations(
  data = "all",
  uuid = NULL,
  query = NULL,
  identifier = NULL,
  identifierType = NULL,
  limit = 100,
  start = NULL,
  curlopts = list()
)
}
\arguments{
\item{data}{The type of data to get. One or more of: 'contact', 'endpoint',
'dataset', 'comment', 'deleted', 'nonPublishing', or the special 'all'.
Default: \code{'all'}}

\item{uuid}{UUID of the data node provider. This must be specified if data
is anything other than 'all'.}

\item{query}{Query nodes. Only used when \code{data='all'}. Ignored
otherwise.}

\item{identifier}{The value for this parameter can be a simple string or
integer, e.g. \code{identifier=120}. This parameter doesn't seem to work right
now.}

\item{identifierType}{Used in combination with the identifier parameter to
filter identifiers by identifier type. See details. This parameter doesn't
seem to work right now.}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\description{
Installations metadata.
}
\details{
identifierType options:

\itemize{
\item {DOI} No description.
\item {FTP} No description.
\item {GBIF_NODE} Identifies the node (e.g: \code{DK} for Denmark, \code{sp2000}
for Species 2000).
\item {GBIF_PARTICIPANT} Participant identifier from the GBIF IMS
Filemaker system.
\item {GBIF_PORTAL} Indicates the identifier originated from an
auto_increment column in the portal.data_provider or portal.data_resource
table respectively.
\item {HANDLER} No description.
\item {LSID} Reference controlled by a separate system, used for example
by DOI.
\item {SOURCE_ID} No description.
\item {UNKNOWN} No description.
\item {URI} No description.
\item {URL} No description.
\item {UUID} No description.
}
}
\examples{
\dontrun{
installations(limit=5)
installations(query="france", limit = 25)
installations(uuid="b77901f9-d9b0-47fa-94e0-dd96450aa2b4")
installations(data='contact', uuid="2e029a0c-87af-42e6-87d7-f38a50b78201")
installations(data='endpoint', uuid="b77901f9-d9b0-47fa-94e0-dd96450aa2b4")
installations(data='dataset', uuid="b77901f9-d9b0-47fa-94e0-dd96450aa2b4")
installations(data='deleted', limit = 25)
installations(data='deleted', limit=2)
installations(data=c('deleted','nonPublishing'), limit=2)
installations(identifierType='DOI', limit=2)

# Pass on curl options
installations(data='deleted', curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#installations}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{get_credentials}
\alias{get_credentials}
\title{Get Github credentials from use in console}
\usage{
get_credentials(...)
}
\description{
This function is defunct.  See the package gistr for similar functionality.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_issues.R
\name{name_issues}
\alias{name_issues}
\title{Parse and examine further GBIF name issues on a dataset.}
\usage{
name_issues(.data, ..., mutate = NULL)
}
\arguments{
\item{.data}{Output from a call to \code{\link[=name_usage]{name_usage()}}}

\item{...}{Named parameters to only get back (e.g. bbmn), or to
remove (e.g. -bbmn).}

\item{mutate}{(character) One of:
\itemize{
\item \code{split} Split issues into new columns.
\item \code{expand} Expand issue abbreviated codes into descriptive names.
for downloads datasets, this is not super useful since the
issues come to you as expanded already.
\item \code{split_expand} Split into new columns, and expand issue names.
}

For split and split_expand, values in cells become y ("yes") or n ("no")}
}
\description{
Parse and examine further GBIF name issues on a dataset.
}
\examples{
\dontrun{
# what do issues mean, can print whole table
head(gbif_issues())
# or just name related issues
gbif_issues()[which(gbif_issues()$type \%in\% c("name")),]
# or search for matches
gbif_issues()[gbif_issues()$code \%in\% c('bbmn','clasna','scina'),]
# compare out data to after name_issues use
(aa <- name_usage(name = "Lupus"))
aa \%>\% name_issues("clasna")

## or parse issues in various ways
### remove data rows with certain issue classes
aa \%>\% name_issues(-clasna, -scina)

### expand issues to more descriptive names
aa \%>\% name_issues(mutate = "expand")

### split and expand
aa \%>\% name_issues(mutate = "split_expand")

### split, expand, and remove an issue class
aa \%>\% name_issues(-bbmn, mutate = "split_expand")

## Or you can use name_issues without \%>\%
name_issues(aa, -bbmn, mutate = "split_expand")
}
}
\references{
\url{https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/NameUsageIssue.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgb_country_codes.r
\name{rgb_country_codes}
\alias{rgb_country_codes}
\title{Look up 2 character ISO country codes}
\usage{
rgb_country_codes(country_name, fuzzy = FALSE, ...)
}
\arguments{
\item{country_name}{Name of country to look up}

\item{fuzzy}{If TRUE, uses agrep to do fuzzy search on names.}

\item{...}{Further arguments passed on to agrep or grep}
}
\description{
Look up 2 character ISO country codes
}
\examples{
rgb_country_codes(country_name="United")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{resources}
\alias{resources}
\title{Get data resources and their unique keys.}
\usage{
resources(...)
}
\description{
This function is defunct.
}
\seealso{
networks organizations datasets
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gbifmap_list}
\alias{gbifmap_list}
\title{Make a simple map to visualize GBIF point data.}
\usage{
gbifmap_list(...)
}
\description{
This function is defunct.
}
\seealso{
gbifmap
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgbif-package.r
\docType{data}
\name{occ_fields}
\alias{occ_fields}
\title{Vector of fields in the output for the function \code{\link[=occ_search]{occ_search()}}}
\description{
These fields can be specified in the \code{fields} parameer in the
\code{\link[=occ_search]{occ_search()}} function.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.r
\name{datasets}
\alias{datasets}
\title{Search for datasets and dataset metadata.}
\usage{
datasets(
  data = "all",
  type = NULL,
  uuid = NULL,
  query = NULL,
  id = NULL,
  limit = 100,
  start = NULL,
  curlopts = list()
)
}
\arguments{
\item{data}{The type of data to get. One or more of: 'organization',
'contact', 'endpoint', 'identifier', 'tag', 'machinetag', 'comment',
'constituents', 'document', 'metadata', 'deleted', 'duplicate',
'subDataset', 'withNoEndpoint', or the special 'all'. Default: \code{all}}

\item{type}{Type of dataset. Options: include occurrence, checklist,
metadata, or sampling_event.}

\item{uuid}{UUID of the data node provider. This must be specified if data
is anything other than \code{all}}

\item{query}{Query term(s). Only used when \code{data=all}}

\item{id}{A metadata document id.}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A list.
}
\description{
Search for datasets and dataset metadata.
}
\examples{
\dontrun{
datasets(limit=5)
datasets(type="occurrence", limit=10)
datasets(uuid="a6998220-7e3a-485d-9cd6-73076bd85657")
datasets(data='contact', uuid="a6998220-7e3a-485d-9cd6-73076bd85657")
datasets(data='metadata', uuid="a6998220-7e3a-485d-9cd6-73076bd85657")
datasets(data='metadata', uuid="a6998220-7e3a-485d-9cd6-73076bd85657",
  id=598)
datasets(data=c('deleted','duplicate'))
datasets(data=c('deleted','duplicate'), limit=1)

# curl options
datasets(data=c('deleted','duplicate'), curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#datasets}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
Pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset_metrics.r
\name{dataset_metrics}
\alias{dataset_metrics}
\title{Get details on a GBIF dataset.}
\usage{
dataset_metrics(uuid, curlopts = list())
}
\arguments{
\item{uuid}{(character) One or more dataset UUIDs. See examples.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\description{
Get details on a GBIF dataset.
}
\note{
Dataset metrics are only available for checklist type datasets.
}
\examples{
\dontrun{
dataset_metrics(uuid='863e6d6b-f602-4495-ac30-881482b6f799')
dataset_metrics(uuid='66dd0960-2d7d-46ee-a491-87b9adcfe7b1')
dataset_metrics(uuid=c('863e6d6b-f602-4495-ac30-881482b6f799',
   '66dd0960-2d7d-46ee-a491-87b9adcfe7b1'))
dataset_metrics(uuid='66dd0960-2d7d-46ee-a491-87b9adcfe7b1',
  curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#datasetMetrics}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{create_gist}
\alias{create_gist}
\title{Function that takes a list of files and creates payload for API}
\usage{
create_gist(...)
}
\description{
This function is defunct.  See the package gistr for similar functionality.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organizations.r
\name{organizations}
\alias{organizations}
\title{Organizations metadata.}
\usage{
organizations(
  data = "all",
  uuid = NULL,
  query = NULL,
  limit = 100,
  start = NULL,
  curlopts = list()
)
}
\arguments{
\item{data}{(character) The type of data to get. One or more of:
'organization', 'contact', 'endpoint', 'identifier', 'tag', 'machineTag',
'comment', 'hostedDataset', 'ownedDataset', 'deleted', 'pending',
'nonPublishing', or the special 'all'. Default: \code{'all'}}

\item{uuid}{(character) UUID of the data node provider. This must be
specified if data is anything other than 'all'.}

\item{query}{(character) Query nodes. Only used when \code{data='all'}}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A list of length one or two. If \code{uuid} is NULL, then a
data.frame with call metadata, and a data.frame, but if \code{uuid} given,
then a list.
}
\description{
Organizations metadata.
}
\examples{
\dontrun{
organizations(limit=5)
organizations(query="france", limit=5)
organizations(uuid="4b4b2111-ee51-45f5-bf5e-f535f4a1c9dc")
organizations(data='contact', uuid="4b4b2111-ee51-45f5-bf5e-f535f4a1c9dc")
organizations(data='pending')
organizations(data=c('contact','endpoint'),
  uuid="4b4b2111-ee51-45f5-bf5e-f535f4a1c9dc")

# Pass on curl options
organizations(query="spain", curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#organizations}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbif_citation.R
\name{gbif_citation}
\alias{gbif_citation}
\title{Get citation for datasets used}
\usage{
gbif_citation(x)
}
\arguments{
\item{x}{(character) Result of call to \code{\link[=occ_search]{occ_search()}}, \code{\link[=occ_data]{occ_data()}},
\code{\link[=occ_download_get]{occ_download_get()}}, \code{\link[=occ_download_meta]{occ_download_meta()}}, a dataset key, or occurrence key
(character or numeric)}
}
\value{
list with S3 class assigned, used by a print method to pretty print
citation information. Though you can unclass the output or just index to the
named items as needed.
}
\description{
Get citation for datasets used
}
\details{
Returns a set of citations, one for each dataset. We pull out
unique dataset keys and get citations, so the length of citations may not
be equal to the number of records you pass in.

Currently, this function gives back citations at the dataset level, not
at the individual occurrence level. If occurrence keys are passed in, then
we track down the dataset the key is from, and get the citation for
the dataset.
}
\examples{
\dontrun{
res1 <- occ_search(taxonKey=9206251, limit=2)
(xx <- gbif_citation(res1))

# each individual citation object is a list
## rights and/or citation may be NULL
xx[[1]]
xx[[1]]$rights
xx[[1]]$citation
xx[[1]]$citation$title
xx[[1]]$citation$text
xx[[1]]$citation$accessed
xx[[1]]$citation$citation

## access many citations
unlist(lapply(xx, "[[", c("citation", "citation")))

res2 <- occ_search(datasetKey='7b5d6a48-f762-11e1-a439-00145eb45e9a',
limit=20)
(xx <- gbif_citation(res2))

# if no datasetKey field included, we attempt to identify the dataset
## key field included - still works
res3 <- occ_search(taxonKey=9206251, fields=c('name','basisOfRecord','key'),
  limit=20)
(xx <- gbif_citation(res3))
## key field not included - errors
# res3 <- occ_search(taxonKey=9206251, fields=c('name','basisOfRecord','
#    protocol'), limit=20)
# (xx <- gbif_citation(res3))

# occ_data
res1 <- occ_data(taxonKey=9206251, limit=2)
(xx <- gbif_citation(res1))

# character class inputs
## pass in a dataset key
gbif_citation(x='0ec3229f-2b53-484e-817a-de8ceb1fce2b')
## pass in an occurrence key
# gbif_citation(x='1101144669')

# pass in an occurrence key as a numeric (won't work for a dataset key)
# gbif_citation(x=1101144669)

# Downloads
## occ_download_get()
# d1 <- occ_download(pred("country", "BG"), pred_gte("year", 2020))
# occ_download_meta(d1) # wait until status = succeeded
# d1 <- occ_download_get(d1, overwrite = TRUE)
# gbif_citation(d1)

## occ_download_meta()
# key <- "0000122-171020152545675"
# res <- occ_download_meta(key)
# gbif_citation(res)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derived_dataset.R
\name{derived_dataset}
\alias{derived_dataset}
\alias{derived_dataset_prep}
\title{Register a derived dataset for citation.}
\usage{
derived_dataset(
  citation_data = NULL,
  title = NULL,
  description = NULL,
  source_url = NULL,
  gbif_download_doi = NULL,
  user = NULL,
  pwd = NULL,
  curlopts = list()
)

derived_dataset_prep(
  citation_data = NULL,
  title = NULL,
  description = NULL,
  source_url = NULL,
  gbif_download_doi = NULL,
  user = NULL,
  pwd = NULL,
  curlopts = list()
)
}
\arguments{
\item{citation_data}{(required) A data.frame with \strong{two columns}. The first
column should be GBIF \strong{datasetkey uuids} and the second column should be
\strong{occurrence counts} from each of your datasets, representing the
contribution of each GBIF dataset to your final derived dataset.}

\item{title}{(required) The title for your derived dataset.}

\item{description}{(required) A description of the dataset. Perhaps
describing how it was created.}

\item{source_url}{(required) A link to where the dataset is stored.}

\item{gbif_download_doi}{(optional) A DOI from an original GBIF download.}

\item{user}{(required) Your GBIF username.}

\item{pwd}{(required) Your GBIF password.}

\item{curlopts}{a list of arguments to pass to curl.}
}
\value{
A list.
}
\description{
Register a derived dataset for citation.
}
\section{Usage}{

Create a \strong{citable DOI} for a dataset derived from GBIF mediated
occurrences.

\strong{Use-case (1)} your dataset was obtained with \code{occ_search} and
never returned a \strong{citable DOI}, but you want to cite the data in a
research paper.

\strong{Use-case (2)} your dataset was obtained using \code{occ_download} and you
got a DOI, but the data underwent extensive filtering using
\code{CoordinateCleaner} or some other cleaning pipeline. In this case be sure
to fill in your original \code{gbif_download_doi}.

\strong{Use-case (3)} your dataset was generated using a GBIF cloud export but
you want a DOI to cite in your research paper.

Use \code{derived_dataset} to create a custom citable meta-data description and
most importantly a DOI link between an external archive (e.g. Zenodo) and the
datasets involved in your research or analysis.

All fields (except \code{gbif_download_doi}) are required for the registration to
work.

We recommend that you run \code{derived_dataset_prep} to check registration
details before making it final with \code{derived_dataset}.
}

\section{Authentication}{

Some \code{rgbif} functions require your \strong{GBIF credentials}.

For the \code{user} and \code{pwd} parameters, you can set them in one of
three ways:
\enumerate{
\item Set them in your \code{.Renviron}/\code{.bash_profile} (or similar) file with the
names \code{GBIF_USER}, \code{GBIF_PWD}, and \code{GBIF_EMAIL}
\item Set them in your \code{.Rprofile} file with the names \code{gbif_user} and
\code{gbif_pwd}.
\item Simply pass strings to each of the parameters in the function
call.
}

We strongly recommend the \strong{first option} - storing your details as
environment variables - as it's the most widely used way to store secrets.

You can edit your \code{.Renviron} with \code{usethis::edit_r_environ()}.

After editing, your \code{.Renviron} file should look something like this...

GBIF_USER="jwaller"\cr
GBIF_PWD="fakepassword123"\cr
GBIF_EMAIL="jwaller@gbif.org"\cr

See \code{?Startup} for help.
}

\examples{
\dontrun{
#  data <- data.frame(
#  datasetKey = c(
#  "3ea36590-9b79-46a8-9300-c9ef0bfed7b8",
#  "630eb55d-5169-4473-99d6-a93396aeae38",
#  "806bf7d4-f762-11e1-a439-00145eb45e9a"),
#  count = c(3, 1, 2781)
#  )

## If output looks ok, run derived_dataset to register the dataset
#  derived_dataset_prep(
#  citation_data = data,
#  title = "Test for derived dataset",
#  description = "This data was filtered using a fake protocol",
#  source_url = "https://zenodo.org/record/4246090#.YPGS2OgzZPY"
#  )

#  derived_dataset(
#  citation_data = data,
#  title = "Test for derived dataset",
#  description = "This data was filtered using a fake protocol",
#  source_url = "https://zenodo.org/record/4246090#.YPGS2OgzZPY"
#  )

## Example with occ_search and dplyr
# library(dplyr)

# citation_data <- occ_search(taxonKey=212, limit=20)$data \%>\%
#   group_by(datasetKey) \%>\% 
#   count()

# # You would still need to upload your data to Zenodo or something similar 
# derived_dataset_prep(
#   citation_data = citation_data,
#   title="Bird data downloaded for test",
#   description="This data was downloaded using rgbif::occ_search and was 
#   later uploaded to Zenodo.",
#   source_url="https://zenodo.org/record/4246090#.YPGS2OgzZPY",
#   gbif_download_doi = NULL,
# )


}

}
\references{
\url{https://data-blog.gbif.org/post/derived-datasets/}
\url{https://www.gbif.org/derived-dataset/about}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{taxoncount}
\alias{taxoncount}
\title{Search by taxon to retrieve number of records in GBIF.}
\usage{
taxoncount(...)
}
\description{
This function is defunct.
}
\seealso{
occ_count
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_queue.R
\name{DownReq}
\alias{DownReq}
\title{DownReq}
\description{
handles single requests for \link{GbifQueue}
}
\examples{
\dontrun{
res <- DownReq$new(occ_download_prep(pred("basisOfRecord", "LITERATURE"), 
  pred("year", "1956")
))
res
# res$run()
# res
# res$status()
# res$result
}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{req}}{(list) internal holder for the request}

\item{\code{type}}{(list) type, one of "lazy" (to be lazy evaluated) or "pre"
(run with \code{occ_download_exec} internal fxn)}

\item{\code{result}}{(list) holds the result of the http request}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{DownReq$new()}}
\item \href{#method-print}{\code{DownReq$print()}}
\item \href{#method-run}{\code{DownReq$run()}}
\item \href{#method-status}{\code{DownReq$status()}}
\item \href{#method-clone}{\code{DownReq$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{DownReq} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DownReq$new(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{either a lazy object with an object of class \code{occ_download}, or an
object of class \code{occ_download_prep}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{DownReq} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{DownReq} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DownReq$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-run"></a>}}
\if{latex}{\out{\hypertarget{method-run}{}}}
\subsection{Method \code{run()}}{
execute http request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DownReq$run()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing, puts the http response in \verb{$result}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status"></a>}}
\if{latex}{\out{\hypertarget{method-status}{}}}
\subsection{Method \code{status()}}{
check http request status
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DownReq$status()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
output of \code{\link[=occ_download_meta]{occ_download_meta()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DownReq$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbox.R
\name{gbif_bbox2wkt}
\alias{gbif_bbox2wkt}
\alias{gbif_wkt2bbox}
\title{Convert a bounding box to a Well Known Text polygon, and a WKT to a
bounding box}
\usage{
gbif_bbox2wkt(minx = NA, miny = NA, maxx = NA, maxy = NA, bbox = NULL)

gbif_wkt2bbox(wkt = NULL)
}
\arguments{
\item{minx}{(numeric) Minimum x value, or the most western longitude}

\item{miny}{(numeric) Minimum y value, or the most southern latitude}

\item{maxx}{(numeric) Maximum x value, or the most eastern longitude}

\item{maxy}{(numeric) Maximum y value, or the most northern latitude}

\item{bbox}{(numeric) A vector of length 4, with the elements: minx, miny,
maxx, maxy}

\item{wkt}{(character) A Well Known Text object.}
}
\value{
gbif_bbox2wkt returns an object of class charactere, a Well
Known Text string of the form
'POLYGON((minx miny, maxx miny, maxx maxy, minx maxy, minx miny))'.

gbif_wkt2bbox returns a numeric vector of length 4, like
c(minx, miny, maxx, maxy)
}
\description{
Convert a bounding box to a Well Known Text polygon, and a WKT to a
bounding box
}
\examples{
\dontrun{
# Convert a bounding box to a WKT
## Pass in a vector of length 4 with all values
gbif_bbox2wkt(bbox=c(-125.0,38.4,-121.8,40.9))

## Or pass in each value separately
gbif_bbox2wkt(minx=-125.0, miny=38.4, maxx=-121.8, maxy=40.9)

# Convert a WKT object to a bounding box
wkt <- "POLYGON((-125 38.4,-125 40.9,-121.8 40.9,-121.8 38.4,-125 38.4))"
gbif_wkt2bbox(wkt)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_queue.R
\name{occ_download_queue}
\alias{occ_download_queue}
\title{Download requests in a queue}
\usage{
occ_download_queue(..., .list = list(), status_ping = 10)
}
\arguments{
\item{...}{any number of \code{\link[=occ_download]{occ_download()}} requests}

\item{.list}{any number of \code{\link[=occ_download_prep]{occ_download_prep()}} requests}

\item{status_ping}{(integer) seconds between pings checking status of
the download request. generally larger numbers for larger requests.
default: 10 (i.e., 10 seconds). must be 10 or greater}
}
\value{
a list of \code{occ_download} class objects, see \code{\link[=occ_download_get]{occ_download_get()}}
to fetch data
}
\description{
Download requests in a queue
}
\details{
This function is a convenience wrapper around \code{\link[=occ_download]{occ_download()}},
allowing the user to kick off any number of requests, while abiding by
GBIF rules of 3 concurrent requests per user.
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\section{How it works}{

It works by using lazy evaluation to collect your requests into a queue
(but does not use lazy evaluation if use the \code{.list} parameter).
Then it kicks of the first 3 requests. Then in a while loop, we check
status of those requests, and when any request finishes (see
\verb{When is a job done?} below), we kick off the
next, and so on. So in theory, there may not always strictly be 3 running
concurrently, but the function will usually provide for 3 running
concurrently.
}

\section{When is a job done?}{

We mark a job as done by checking the \verb{/occurrence/download/} API route
with our \code{\link[=occ_download_meta]{occ_download_meta()}} function. If the status of the job is
any of "succeeded", "killed", or "cancelled", then we mark the job as done
and move on to other jobs in the queue.
}

\section{Beware}{

This function is still in development. There's a lot of complexity
to this problem. We'll be rolling out fixes and improvements in future
versions of the package, so expect to have to adjust your code
with new versions.
}

\examples{
\dontrun{
if (interactive()) { # dont run in automated example runs, too costly
# passing occ_download() requests via ...
out <- occ_download_queue(
  occ_download(pred('taxonKey', 3119195), pred("year", 1976)),
  occ_download(pred('taxonKey', 3119195), pred("year", 2001)),
  occ_download(pred('taxonKey', 3119195), pred("year", 2001),
    pred_lte("month", 8)),
  occ_download(pred('taxonKey', 5229208), pred("year", 2011)),
  occ_download(pred('taxonKey', 2480946), pred("year", 2015)),
  occ_download(pred("country", "NZ"), pred("year", 1999),
    pred("month", 3)),
  occ_download(pred("catalogNumber", "Bird.27847588"),
    pred("year", 1998), pred("month", 2))
)

# supports <= 3 requests too
out <- occ_download_queue(
  occ_download(pred("country", "NZ"), pred("year", 1999), pred("month", 3)),
  occ_download(pred("catalogNumber", "Bird.27847588"), pred("year", 1998),
    pred("month", 2))
)

# using pre-prepared requests via .list
keys <- c(7905507, 5384395, 8911082)
queries <- list()
for (i in seq_along(keys)) {
  queries[[i]] <- occ_download_prep(
    pred("taxonKey", keys[i]),
    pred_in("basisOfRecord", c("HUMAN_OBSERVATION","OBSERVATION")),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("year", 1993)
  )
}
out <- occ_download_queue(.list = queries)
out

# another pre-prepared example
yrs <- 1930:1934
queries <- list()
for (i in seq_along(yrs)) {
  queries[[i]] <- occ_download_prep(
    pred("taxonKey", 2877951),
    pred_in("basisOfRecord", c("HUMAN_OBSERVATION","OBSERVATION")),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("year", yrs[i])
  )
}
out <- occ_download_queue(.list = queries)
out
}}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset_suggest.r
\name{dataset_suggest}
\alias{dataset_suggest}
\title{Suggest datasets in GBIF.}
\usage{
dataset_suggest(
  query = NULL,
  country = NULL,
  type = NULL,
  subtype = NULL,
  keyword = NULL,
  publishingOrg = NULL,
  hostingOrg = NULL,
  publishingCountry = NULL,
  decade = NULL,
  continent = NULL,
  limit = 100,
  start = NULL,
  pretty = FALSE,
  description = FALSE,
  curlopts = list()
)
}
\arguments{
\item{query}{Query term(s) for full text search.  The value for this
parameter can be a simple word or a phrase. Wildcards can be added to the
simple word parameters only, e.g. \code{q=*puma*}}

\item{country}{NOT YET IMPLEMENTED. Filters by country as given in
isocodes$gbif_name, e.g. \code{country=CANADA}}

\item{type}{Type of dataset, options include occurrene, metadata, checklist,
sampling_event
(http://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/DatasetType.html)}

\item{subtype}{NOT YET IMPLEMENTED. Will allow filtering of datasets by
their dataset subtypes, DC or EML.}

\item{keyword}{Keyword to search by. Datasets can be tagged by keywords,
which you can search on. The search is done on the merged collection of
tags, the dataset keywordCollections and temporalCoverages.}

\item{publishingOrg}{Publishing organization. A uuid string. See
\code{\link{organizations}}}

\item{hostingOrg}{Hosting organization. A uuid string. See
\code{\link{organizations}}}

\item{publishingCountry}{Publishing country. See options at
isocodes$gbif_name}

\item{decade}{Decade, e.g., 1980. Filters datasets by their temporal coverage
broken down to decades. Decades are given as a full year, e.g. 1880, 1960,
2000, etc, and will return datasets wholly contained in the decade as well
as those that cover the entire decade or more. Facet by decade to get the
break down, e.g. /search?facet=DECADE&facet_only=true (see example below)}

\item{continent}{Not yet implemented, but will eventually allow filtering
datasets by their continent(s) as given in our Continent enum.}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{pretty}{Print informative metadata using \code{\link{cat}}. Not easy
to manipulate output though.}

\item{description}{Return descriptions only (TRUE) or all data (FALSE,
default)}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A data.frame, list, or message printed to console (using
\code{pretty=TRUE}).
}
\description{
Suggest datasets in GBIF.
}
\section{Repeat parmeter inputs}{

Some parameters can tak emany inputs, and treated as 'OR' (e.g., a or b or
c). The following take many inputs:
\itemize{
\item \code{type}
\item \code{keyword}
\item \code{publishingOrg}
\item \code{hostingOrg}
\item \code{publishingCountry}
\item \code{decade}
}
}

\examples{
\dontrun{
# Suggest datasets of type "OCCURRENCE".
# dataset_suggest(query="Amazon", type="OCCURRENCE")

# Suggest datasets tagged with keyword "france".
# dataset_suggest(keyword="france")

# Fulltext search for all datasets having the word "amsterdam" somewhere in
# its metadata (title, description, etc).
# dataset_suggest(query="amsterdam")

# Limited search
# dataset_suggest(type="OCCURRENCE", limit=2)
# dataset_suggest(type="OCCURRENCE", limit=2, start=10)

# Return just descriptions
# dataset_suggest(type="OCCURRENCE", limit = 5, description=TRUE)

# Return metadata in a more human readable way (hard to manipulate though)
# dataset_suggest(type="OCCURRENCE", limit = 5, pretty=TRUE)

# Search by country code. Lookup isocodes first, and use US for United States
isocodes[agrep("UNITED", isocodes$gbif_name),]
# dataset_suggest(country="US", limit = 25)

# Search by decade
# dataset_suggest(decade=1980, limit = 30)

# Some parameters accept many inputs, treated as OR
# dataset_suggest(type = c("metadata", "checklist"))
# dataset_suggest(keyword = c("fern", "algae"))
# dataset_suggest(publishingOrg = c("e2e717bf-551a-4917-bdc9-4fa0f342c530",
#   "90fd6680-349f-11d8-aa2d-b8a03c50a862"))
# dataset_suggest(hostingOrg = c("c5f7ef70-e233-11d9-a4d6-b8a03c50a862",
#   "c5e4331-7f2f-4a8d-aa56-81ece7014fc8"))
# dataset_suggest(publishingCountry = c("DE", "NZ"))
# dataset_suggest(decade = c(1910, 1930))

# curl options
# dataset_suggest(type="OCCURRENCE", limit = 2, curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#datasetSearch}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_usage.r
\name{name_usage}
\alias{name_usage}
\title{Lookup details for specific names in all taxonomies in GBIF.}
\usage{
name_usage(
  key = NULL,
  name = NULL,
  data = "all",
  language = NULL,
  datasetKey = NULL,
  uuid = NULL,
  rank = NULL,
  shortname = NULL,
  start = 0,
  limit = 100,
  return = NULL,
  curlopts = list()
)
}
\arguments{
\item{key}{(numeric or character) A GBIF key for a taxon}

\item{name}{(character) Filters by a case insensitive, canonical namestring,
e.g. 'Puma concolor'}

\item{data}{(character) Specify an option to select what data is returned.
See Description below.}

\item{language}{(character) Language, default is english}

\item{datasetKey}{(character) Filters by the dataset's key (a uuid). Must
be length=1}

\item{uuid}{(character) A dataset key}

\item{rank}{(character) Taxonomic rank. Filters by taxonomic rank as
one of: CLASS, CULTIVAR, CULTIVAR_GROUP, DOMAIN, FAMILY, FORM, GENUS,
INFORMAL, INFRAGENERIC_NAME, INFRAORDER, INFRASPECIFIC_NAME,
INFRASUBSPECIFIC_NAME, KINGDOM, ORDER, PHYLUM, SECTION, SERIES, SPECIES,
STRAIN, SUBCLASS, SUBFAMILY, SUBFORM, SUBGENUS, SUBKINGDOM, SUBORDER,
SUBPHYLUM, SUBSECTION, SUBSERIES, SUBSPECIES, SUBTRIBE, SUBVARIETY,
SUPERCLASS, SUPERFAMILY, SUPERORDER, SUPERPHYLUM, SUPRAGENERIC_NAME,
TRIBE, UNRANKED, VARIETY}

\item{shortname}{(character) A short name for a dataset - it may
not do anything}

\item{start}{Record number to start at. Default: 0.}

\item{limit}{Number of records to return. Default: 100.}

\item{return}{Defunct. All components are returned; index to the
one(s) you want}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
An object of class gbif, which is a S3 class list, with slots for
metadata (\code{meta}) and the data itself (\code{data}). In addition, the
object has attributes listing the user supplied arguments and type of
search, which is, differently from occurrence data, always equals to
'single' even if multiple values for some parameters are given. \code{meta}
is a list of length four with offset, limit, endOfRecords and count fields.
\code{data} is a tibble (aka data.frame) containing all information about
the found taxa.
}
\description{
Lookup details for specific names in all taxonomies in GBIF.
}
\details{
This service uses fuzzy lookup so that you can put in partial names and
you should get back those things that match. See examples below.

This function is different from \code{\link[=name_lookup]{name_lookup()}} in that that function
searches for names. This function encompasses a bunch of API endpoints,
most of which require that you already have a taxon key, but there is one
endpoint that allows name searches (see examples below).

Note that \code{data="verbatim"} hasn't been working.

Options for the data parameter are: 'all', 'verbatim', 'name', 'parents',
'children', 'related', 'synonyms', 'descriptions','distributions', 'media',
'references', 'speciesProfiles', 'vernacularNames', 'typeSpecimens', 'root'

This function used to be vectorized with respect to the \code{data}
parameter, where you could pass in multiple values and the function
internally loops over each option making separate requests. This has been
removed. You can still loop over many options for the \code{data} parameter,
just use an \code{lapply} family function, or a for loop, etc.

See \code{\link[=name_issues]{name_issues()}} for more information about issues in \code{issues} column.
}
\section{Repeat parameter inputs}{

These parameters used to accept many inputs, but no longer do:

\itemize{
\item \strong{rank}
\item \strong{name}
\item \strong{langugae}
\item \strong{datasetKey}
}

see also \code{\link{many-values}}
}

\examples{
\dontrun{
# A single name usage
name_usage(key=1)

# Name usage for a taxonomic name
name_usage(name='Puma', rank="GENUS")

# Name usage for all taxa in a dataset
# (set sufficient high limit, but less than 100000)
# name_usage(datasetKey = "9ff7d317-609b-4c08-bd86-3bc404b77c42", 
#  limit = 10000)
# All name usages
name_usage()

# References for a name usage
name_usage(key=2435099, data='references')

# Species profiles, descriptions
name_usage(key=3119195, data='speciesProfiles')
name_usage(key=3119195, data='descriptions')
name_usage(key=2435099, data='children')

# Vernacular names for a name usage
name_usage(key=3119195, data='vernacularNames')

# Limit number of results returned
name_usage(key=3119195, data='vernacularNames', limit=3)

# Search for names by dataset with datasetKey parameter
name_usage(datasetKey="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c")

# Search for a particular language
name_usage(key=3119195, language="FRENCH", data='vernacularNames')

# get root usage with a uuid
name_usage(data = "root", uuid = "73605f3a-af85-4ade-bbc5-522bfb90d847")

# search by language
name_usage(language = "spanish")

# Pass on curl options
name_usage(name='Puma concolor', limit=300, curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/species#nameUsages}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download.R
\name{occ_download}
\alias{occ_download}
\alias{occ_download_prep}
\title{Spin up a download request for GBIF occurrence data.}
\usage{
occ_download(
  ...,
  body = NULL,
  type = "and",
  format = "DWCA",
  user = NULL,
  pwd = NULL,
  email = NULL,
  curlopts = list()
)

occ_download_prep(
  ...,
  body = NULL,
  type = "and",
  format = "DWCA",
  user = NULL,
  pwd = NULL,
  email = NULL,
  curlopts = list()
)
}
\arguments{
\item{...}{For \code{occ_download()} and \code{occ_download_prep()}, one or more
objects of class \code{occ_predicate} or \code{occ_predicate_list}, created by
\verb{pred*} functions (see \link{download_predicate_dsl}). If you use this, don't
use \code{body} parameter.}

\item{body}{if you prefer to pass in the payload yourself, use this
parameter. if use this, don't pass anythig to the dots. accepts
either an R list, or JSON. JSON is likely easier, since the JSON
library \pkg{jsonlite} requires that you unbox strings that shouldn't
be auto-converted to arrays, which is a bit tedious for large queries.
optional}

\item{type}{(character) One of equals (=), and (&), or (|), lessThan (<),
lessThanOrEquals (<=), greaterThan (>), greaterThanOrEquals (>=), in,
within, not (!), like, isNotNull}

\item{format}{(character) The download format. One of 'DWCA' (default),
'SIMPLE_CSV', or 'SPECIES_LIST'}

\item{user}{(character) User name within GBIF's website. Required. See
"Authentication" below}

\item{pwd}{(character) User password within GBIF's website. Required. See
"Authentication" below}

\item{email}{(character) Email address to recieve download notice done
email. Required. See "Authentication" below}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\description{
Spin up a download request for GBIF occurrence data.
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\section{geometry}{

When using the geometry parameter, make sure that your well known text
(WKT) is formatted as GBIF expects it. They expect WKT to have a
counter-clockwise winding order. For example, the following is clockwise
\verb{POLYGON((-19.5 34.1, -25.3 68.1, 35.9 68.1, 27.8 34.1, -19.5 34.1))},
whereas they expect the other order:
\verb{POLYGON((-19.5 34.1, 27.8 34.1, 35.9 68.1, -25.3 68.1, -19.5 34.1))}

note that coordinate pairs are \verb{longitude latitude}, longitude first, then
latitude

you should not get any results if you supply WKT that has clockwise
winding order.

also note that \code{\link[=occ_search]{occ_search()}}/\code{\link[=occ_data]{occ_data()}} behave differently with
respect to WKT in that you can supply clockwise WKT to those
functions but they treat it as an exclusion, so get all data not
inside the WKT area.
}

\section{Methods}{

\itemize{
\item \code{occ_download_prep}: prepares a download request, but DOES NOT execute it.
meant for use with \code{\link[=occ_download_queue]{occ_download_queue()}}
\item \code{occ_download}: prepares a download request and DOES execute it
}
}

\section{Authentication}{

For \code{user}, \code{pwd}, and \code{email} parameters, you can set them in one of
three ways:
\itemize{
\item Set them in your \code{.Rprofile} file with the names \code{gbif_user},
\code{gbif_pwd}, and \code{gbif_email}
\item Set them in your \code{.Renviron}/\code{.bash_profile} (or similar) file with the
names \code{GBIF_USER}, \code{GBIF_PWD}, and \code{GBIF_EMAIL}
\item Simply pass strings to each of the parameters in the function
call
}

We strongly recommend the second option - storing your details as
environment variables as it's the most widely used way to store secrets.

See \code{?Startup} for help.
}

\section{Query length}{

GBIF has a limit of 12,000 characters for a download query. This means
that you can have a pretty long query, but at some point it may lead to an
error on GBIF's side and you'll have to split your query into a few.
}

\examples{
\dontrun{
# occ_download(pred("basisOfRecord", "LITERATURE"))
# occ_download(pred("taxonKey", 3119195), pred_gt("elevation", 5000))
# occ_download(pred_gt("decimalLatitude", 50))
# occ_download(pred_gte("elevation", 9000))
# occ_download(pred_gte('decimalLatitude", 65))
# occ_download(pred("country", "US"))
# occ_download(pred("institutionCode", "TLMF"))
# occ_download(pred("catalogNumber", 217880))

# download format
# z <- occ_download(pred_gte("decimalLatitude", 75),
#  format = "SPECIES_LIST")

# res <- occ_download(pred("taxonKey", 7264332), pred("hasCoordinate", TRUE))

# pass output directly, or later, to occ_download_meta for more information
# occ_download(pred_gt('decimalLatitude', 75)) \%>\% occ_download_meta

# Multiple queries
# occ_download(pred_gte("decimalLatitude", 65),
#  pred_lte("decimalLatitude", -65), type="or")
# gg <- occ_download(pred("depth", 80), pred("taxonKey", 2343454),
#  type="or")
# x <- occ_download(pred_and(pred_within("POLYGON((-14 42, 9 38, -7 26, -14 42))"),
#  pred_gte("elevation", 5000)))

# complex example with many predicates
# shows example of how to do date ranges for both year and month
# res <- occ_download(
#  pred_gt("elevation", 5000),
#  pred_in("basisOfRecord", c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
#  pred("country", "US"),
#  pred("hasCoordinate", TRUE),
#  pred("hasGeospatialIssue", FALSE),
#  pred_gte("year", 1999),
#  pred_lte("year", 2011),
#  pred_gte("month", 3),
#  pred_lte("month", 8)
# )

# Using body parameter - pass in your own complete query
## as JSON
query1 <- '{"creator":"sckott",
  "notification_address":["stuff1@gmail.com"],
  "predicate":{"type":"and","predicates":[
    {"type":"equals","key":"TAXON_KEY","value":"7264332"},
    {"type":"equals","key":"HAS_COORDINATE","value":"TRUE"}]}
 }'
# res <- occ_download(body = query1, curlopts=list(verbose=TRUE))

## as a list
library(jsonlite)
query <- list(
  creator = unbox("sckott"),
  notification_address = "stuff1@gmail.com",
  predicate = list(
    type = unbox("and"),
    predicates = list(
      list(type = unbox("equals"), key = unbox("TAXON_KEY"),
        value = unbox("7264332")),
      list(type = unbox("equals"), key = unbox("HAS_COORDINATE"),
        value = unbox("TRUE"))
    )
  )
)
# res <- occ_download(body = query, curlopts = list(verbose = TRUE))

# Prepared query
occ_download_prep(pred("basisOfRecord", "LITERATURE"))
occ_download_prep(pred("basisOfRecord", "LITERATURE"), format = "SIMPLE_CSV")
occ_download_prep(pred("basisOfRecord", "LITERATURE"), format = "SPECIES_LIST")
occ_download_prep(pred_in("taxonKey", c(2977832, 2977901, 2977966, 2977835)))
occ_download_prep(pred_within("POLYGON((-14 42, 9 38, -7 26, -14 42))"))

## a complicated example
occ_download_prep(
  pred_in("basisOfRecord", c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION")),
  pred_in("taxonKey", c(2498343, 2481776, 2481890)),
  pred_in("country", c("GB", "IE")),
  pred_or(pred_lte("year", 1989), pred("year", 2000))
)

# x = occ_download(
#   pred_in("basisOfRecord", c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION")),
#   pred_in("taxonKey", c(9206251, 3112648)),
#   pred_in("country", c("US", "MX")),
#   pred_and(pred_gte("year", 1989), pred_lte("year", 1991))
# )
# occ_download_meta(x)
# z <- occ_download_get(x)
# df <- occ_download_import(z)
# str(df)
# library(dplyr)
# unique(df$basisOfRecord)
# unique(df$taxonKey)
# unique(df$countryCode)
# sort(unique(df$year))
}
}
\references{
See the API docs
\url{https://www.gbif.org/developer/occurrence#download} for more info,
and the predicates docs
\url{https://www.gbif.org/developer/occurrence#predicates}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_cached.R
\name{occ_download_cached}
\alias{occ_download_cached}
\title{Check for downloads already in your GBIF account}
\usage{
occ_download_cached(
  ...,
  body = NULL,
  type = "and",
  format = "DWCA",
  user = NULL,
  pwd = NULL,
  email = NULL,
  refresh = FALSE,
  age = 30,
  curlopts = list()
)
}
\arguments{
\item{...}{For \code{occ_download()} and \code{occ_download_prep()}, one or more
objects of class \code{occ_predicate} or \code{occ_predicate_list}, created by
\verb{pred*} functions (see \link{download_predicate_dsl}). If you use this, don't
use \code{body} parameter.}

\item{body}{if you prefer to pass in the payload yourself, use this
parameter. if use this, don't pass anythig to the dots. accepts
either an R list, or JSON. JSON is likely easier, since the JSON
library \pkg{jsonlite} requires that you unbox strings that shouldn't
be auto-converted to arrays, which is a bit tedious for large queries.
optional}

\item{type}{(character) One of equals (=), and (&), or (|), lessThan (<),
lessThanOrEquals (<=), greaterThan (>), greaterThanOrEquals (>=), in,
within, not (!), like, isNotNull}

\item{format}{(character) The download format. One of 'DWCA' (default),
'SIMPLE_CSV', or 'SPECIES_LIST'}

\item{user}{(character) User name within GBIF's website. Required. See
"Authentication" below}

\item{pwd}{(character) User password within GBIF's website. Required. See
"Authentication" below}

\item{email}{(character) Email address to recieve download notice done
email. Required. See "Authentication" below}

\item{refresh}{(logical) refresh your list of downloads. on the first
request of each R session we'll cache your stored GBIF occurrence
downloads locally. you can refresh this list by setting \code{refresh=TRUE};
if you're in the same R session, and you've done many download requests,
then refreshing may be a good idea if you're using this function}

\item{age}{(integer) number of days after which you want a new
download. default: 30}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\description{
Check for downloads already in your GBIF account
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\examples{
\dontrun{
# these are examples from the package maintainer's account;
# outcomes will vary by user
occ_download_cached(pred_gte("elevation", 12000L))
occ_download_cached(pred("catalogNumber", 217880))
occ_download_cached(pred_gte("decimalLatitude", 65),
  pred_lte("decimalLatitude", -65), type="or")
occ_download_cached(pred_gte("elevation", 12000L))
occ_download_cached(pred_gte("elevation", 12000L), refresh = TRUE)
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbif_photos.r
\name{gbif_photos}
\alias{gbif_photos}
\title{View photos from GBIF.}
\usage{
gbif_photos(input, output = NULL, which = "table", browse = TRUE)
}
\arguments{
\item{input}{Input output from occ_search}

\item{output}{Output folder path. If not given uses temporary folder.}

\item{which}{One of map or table (default).}

\item{browse}{(logical) Browse output (default: \code{TRUE})}
}
\description{
View photos from GBIF.
}
\details{
The max number of photos you can see when which="map" is ~160,
so cycle through if you have more than that.
}
\section{BEWARE}{
 The maps in the table view may not show up correctly if
you are using RStudio
}

\examples{
\dontrun{
res <- occ_search(mediaType = 'StillImage', limit = 100)
gbif_photos(res)
gbif_photos(res, which='map')

res <- occ_search(scientificName = "Aves", mediaType = 'StillImage',
  limit=150)
gbif_photos(res)
gbif_photos(res, output = '~/barfoo')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_import.R
\name{occ_download_import}
\alias{occ_download_import}
\alias{as.download}
\alias{as.download.character}
\alias{as.download.download}
\title{Import a downloaded file from GBIF.}
\usage{
occ_download_import(
  x = NULL,
  key = NULL,
  path = ".",
  fill = FALSE,
  encoding = "UTF-8",
  ...
)

as.download(path = ".", key = NULL)

\method{as.download}{character}(path = ".", key = NULL)

\method{as.download}{download}(path = ".", key = NULL)
}
\arguments{
\item{x}{The output of a call to \code{occ_download_get}}

\item{key}{A key generated from a request, like that from
\code{occ_download}}

\item{path}{Path to unzip file to. Default: \code{"."} Writes to
folder matching zip file name}

\item{fill}{(logical) (default: \code{FALSE}). If \code{TRUE} then in case
the rows have unequal length, blank fields are implicitly filled.
passed on to \code{fill} parameter in \link[data.table:fread]{data.table::fread}.}

\item{encoding}{(character) encoding to read in data; passed to
\code{\link[data.table:fread]{data.table::fread()}}. default: "UTF-8". other allowed options:
"Latin-1" and "unknown". see \code{?data.table::fread} docs}

\item{...}{parameters passed on to \code{\link[data.table:fread]{data.table::fread()}}. See \code{fread}
docs for details. Some \code{fread} parameters that may be particular useful
here are: \code{select} (select which columns to read in; others are dropped),
\code{nrows} (only read in a certain number of rows)}
}
\value{
a tibble (data.frame)
}
\description{
Import a downloaded file from GBIF.
}
\details{
You can provide either x as input, or both key and path. We use
\code{\link[data.table:fread]{data.table::fread()}} internally to read data.
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\section{Problems reading data}{

You may run into errors when using \code{occ_download_import()}; most often
these are due to \code{\link[data.table:fread]{data.table::fread()}} not being able to parse the
\code{occurrence.txt} file correctly. The \code{fill} parameter passes down to
\code{\link[data.table:fread]{data.table::fread()}} and the \code{...} allows you to pass on any other
parameters that \code{\link[data.table:fread]{data.table::fread()}} accepts. Read the docs for \code{fread}
for help.
}

\section{countryCode result column and Namibia}{

The country code for Namibia is \code{"NA"}. Unfortunately in R an \code{"NA"} string
will be read in to R as an NA/missing. To avoid this, in this function
we read in the data, then convert an NA/missing values to the character
string \code{"NA"}. When a country code is truly missing it will be an empty
string.
}

\examples{
\dontrun{
# First, kick off at least 1 download, then wait for the job to be complete
# Then use your download keys
res <- occ_download_get(key="0000066-140928181241064", overwrite=TRUE)
occ_download_import(res)

occ_download_get(key="0000066-140928181241064", overwrite = TRUE) \%>\%
  occ_download_import

# coerce a file path to the right class to feed to occ_download_import
# as.download("0000066-140928181241064.zip")
# as.download(key = "0000066-140928181241064")
# occ_download_import(as.download("0000066-140928181241064.zip"))

# download a dump that has a CSV file
# res <- occ_download_get(key = "0001369-160509122628363", overwrite=TRUE)
# occ_download_import(res)
# occ_download_import(key = "0001369-160509122628363")

# download and import a species list (in csv format)
# x <- occ_download_get("0000172-190415153152247")
# occ_download_import(x)
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_suggest.r
\name{suggestfields}
\alias{suggestfields}
\title{Fields available in gbif_suggest function}
\usage{
suggestfields()
}
\description{
Fields available in gbif_suggest function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_list.R
\name{occ_download_list}
\alias{occ_download_list}
\title{Lists the downloads created by a user.}
\usage{
occ_download_list(
  user = NULL,
  pwd = NULL,
  limit = 20,
  start = 0,
  curlopts = list()
)
}
\arguments{
\item{user}{(character) User name within GBIF's website. Required. See
Details.}

\item{pwd}{(character) User password within GBIF's website. Required. See
Details.}

\item{limit}{(integer/numeric) Number of records to return. Default: 20,
Max: 1000}

\item{start}{(integer/numeric) Record number to start at. Default: 0}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
a list with two slots:
\itemize{
\item meta: a single row data.frame with columns: \code{offset}, \code{limit},
\code{endofrecords}, \code{count}
\item results: a tibble with the nested data flattened, with many
columns with the same \code{request.} prefix
}
}
\description{
Lists the downloads created by a user.
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\examples{
\dontrun{
occ_download_list(user="sckott")
occ_download_list(user="sckott", limit = 5)
occ_download_list(user="sckott", start = 21)
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_predicate_dsl.R
\name{download_predicate_dsl}
\alias{download_predicate_dsl}
\alias{pred}
\alias{pred_gt}
\alias{pred_gte}
\alias{pred_lt}
\alias{pred_lte}
\alias{pred_not}
\alias{pred_like}
\alias{pred_within}
\alias{pred_notnull}
\alias{pred_or}
\alias{pred_and}
\alias{pred_in}
\title{Download predicate DSL (domain specific language)}
\usage{
pred(key, value)

pred_gt(key, value)

pred_gte(key, value)

pred_lt(key, value)

pred_lte(key, value)

pred_not(...)

pred_like(key, value)

pred_within(value)

pred_notnull(key)

pred_or(..., .list = list())

pred_and(..., .list = list())

pred_in(key, value)
}
\arguments{
\item{key}{(character) the key for the predicate. See "Keys" below}

\item{value}{(various) the value for the predicate}

\item{..., .list}{For \code{pred_or()} or \code{pred_and()}, one or more objects of
class \code{occ_predicate}, created by any \verb{pred*} function}
}
\description{
Download predicate DSL (domain specific language)
}
\section{predicate methods and their equivalent types}{


\verb{pred*} functions are named for the 'type' of operation they do, following
the terminology used by GBIF, see
https://www.gbif.org/developer/occurrence#predicates

Function names are given, with the equivalent GBIF type value (e.g.,
\code{pred_gt} and \code{greaterThan})

The following functions take one key and one value:
\itemize{
\item \code{pred}: equals
\item \code{pred_lt}: lessThan
\item \code{pred_lte}: lessThanOrEquals
\item \code{pred_gt}: greaterThan
\item \code{pred_gte}: greaterThanOrEquals
\item \code{pred_like}: like
}

The following function is only for geospatial queries, and only
accepts a WKT string:
\itemize{
\item \code{pred_within}: within
}

The following function is only for stating the you don't want
a key to be null, so only accepts one key:
\itemize{
\item \code{pred_notnull}: isNotNull
}

The following two functions accept multiple individual predicates,
separating them by either "and" or "or":
\itemize{
\item \code{pred_and}: and
\item \code{pred_or}: or
}

The not predicate accepts one predicate; that is, this negates whatever
predicate is passed in, e.g., not the taxonKey of 12345:
\itemize{
\item \code{pred_not}: not
}

The following function is special in that it accepts a single key
but many values; stating that you want to search for all the values:
\itemize{
\item \code{pred_in}: in
}
}

\section{What happens internally}{

Internally, the input to \verb{pred*} functions turns into JSON to be sent to
GBIF. For example ...

\code{pred_in("taxonKey", c(2480946, 5229208))} gives:\preformatted{\{
   "type": "in",
   "key": "TAXON_KEY",
   "values": ["2480946", "5229208"]
 \}
}

\code{pred_gt("elevation", 5000)} gives:\preformatted{\{
   "type": "greaterThan",
   "key": "ELEVATION",
   "value": "5000"
\}
}

\code{pred_or(pred("taxonKey", 2977832), pred("taxonKey", 2977901))} gives:\preformatted{\{
  "type": "or",
  "predicates": [
     \{
       "type": "equals",
       "key": "TAXON_KEY",
       "value": "2977832"
     \},
     \{
       "type": "equals",
       "key": "TAXON_KEY",
       "value": "2977901"
     \}
  ]
\}
}
}

\section{Keys}{


Acceptable arguments to the \code{key} parameter are (with the version of
the key in parens that must be sent if you pass the query via the \code{body}
parameter; see below for examples). Open an issue in the GitHub
repository for this package if you know of a key that should
be supported that is not yet.
\itemize{
\item taxonKey (TAXON_KEY)
\item scientificName (SCIENTIFIC_NAME)
\item country (COUNTRY)
\item publishingCountry (PUBLISHING_COUNTRY)
\item hasCoordinate (HAS_COORDINATE)
\item hasGeospatialIssue (HAS_GEOSPATIAL_ISSUE)
\item typeStatus (TYPE_STATUS)
\item recordNumber (RECORD_NUMBER)
\item lastInterpreted (LAST_INTERPRETED)
\item continent (CONTINENT)
\item geometry (GEOMETRY)
\item basisOfRecord (BASIS_OF_RECORD)
\item datasetKey (DATASET_KEY)
\item eventDate (EVENT_DATE)
\item catalogNumber (CATALOG_NUMBER)
\item year (YEAR)
\item month (MONTH)
\item decimalLatitude (DECIMAL_LATITUDE)
\item decimalLongitude (DECIMAL_LONGITUDE)
\item elevation (ELEVATION)
\item depth (DEPTH)
\item institutionCode (INSTITUTION_CODE)
\item collectionCode (COLLECTION_CODE)
\item issue (ISSUE)
\item mediatype (MEDIA_TYPE)
\item recordedBy (RECORDED_BY)
\item establishmentMeans (ESTABLISHMENT_MEANS)
\item coordinateUncertaintyInMeters (COORDINATE_UNCERTAINTY_IN_METERS)
\item gadm (GADM_GID) (for the Database of Global Administrative Areas)
\item stateProvince (STATE_PROVINCE)
\item occurrenceStatus (OCCURRENCE_STATUS)
}
}

\examples{
pred("taxonKey", 3119195)
pred_gt("elevation", 5000)
pred_gte("elevation", 5000)
pred_lt("elevation", 1000)
pred_lte("elevation", 1000)
pred_within("POLYGON((-14 42, 9 38, -7 26, -14 42))")
pred_and(pred_within("POLYGON((-14 42, 9 38, -7 26, -14 42))"),
  pred_gte("elevation", 5000))
pred_or(pred_lte("year", 1989), pred("year", 2000))
pred_and(pred_lte("year", 1989), pred("year", 2000))
pred_in("taxonKey", c(2977832, 2977901, 2977966, 2977835))
pred_in("basisOfRecord", c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION"))
pred_not(pred("taxonKey", 729))
pred_like("catalogNumber", "PAPS5-560\%")
pred_notnull("issue")
pred("basisOfRecord", "LITERATURE")
pred("hasCoordinate", TRUE)
pred("stateProvince", "California")
pred("hasGeospatialIssue", FALSE)
pred_within("POLYGON((-14 42, 9 38, -7 26, -14 42))")
pred_or(pred("taxonKey", 2977832), pred("taxonKey", 2977901),
  pred("taxonKey", 2977966))
pred_in("taxonKey", c(2977832, 2977901, 2977966, 2977835))
}
\references{
Download predicates docs:
\url{https://www.gbif.org/developer/occurrence#predicates}
}
\seealso{
Other downloads: 
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_wait.R
\name{occ_download_wait}
\alias{occ_download_wait}
\title{Wait for an occurrence download to be done}
\usage{
occ_download_wait(x, status_ping = 5, curlopts = list(), quiet = FALSE)
}
\arguments{
\item{x}{and object of class \code{occ_download} or downloadkey}

\item{status_ping}{(integer) seconds between each \code{\link[=occ_download_meta]{occ_download_meta()}}
request. default is 5, and cannot be < 3}

\item{curlopts}{(list) curl options, as named list, passed on to
\code{\link[=occ_download_meta]{occ_download_meta()}}}

\item{quiet}{(logical) suppress messages. default: \code{FALSE}}
}
\value{
an object of class \code{occ_download_meta}, see \code{\link[=occ_download_meta]{occ_download_meta()}}
for details
}
\description{
Wait for an occurrence download to be done
}
\note{
\code{\link[=occ_download_queue]{occ_download_queue()}} is similar, but handles many requests
at once; \code{occ_download_wait} handles one request at a time
}
\examples{
\dontrun{
x <- occ_download(
  pred("taxonKey", 9206251),
  pred_in("country", c("US", "MX")),
  pred_gte("year", 1971)
)
res <- occ_download_wait(x)
occ_download_meta(x)

# works also with a downloadkey
occ_download_wait("0000066-140928181241064") 

}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_meta}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{taxonget}
\alias{taxonget}
\title{Get taxonomic information on a specific taxon or taxa in GBIF by their taxon
concept keys.}
\usage{
taxonget(...)
}
\description{
This function is defunct.
}
\seealso{
name_usage
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{occurrencecount}
\alias{occurrencecount}
\title{Counts taxon concept records matching a range of filters.}
\usage{
occurrencecount(...)
}
\description{
This function is defunct.
}
\seealso{
occ_count
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvt_fetch.R
\name{mvt_fetch}
\alias{mvt_fetch}
\title{Fetch Map Vector Tiles (MVT)}
\usage{
mvt_fetch(
  source = "density",
  x = 0,
  y = 0,
  z = 0,
  srs = "EPSG:4326",
  bin = NULL,
  hexPerTile = NULL,
  squareSize = NULL,
  style = "classic.point",
  taxonKey = NULL,
  datasetKey = NULL,
  country = NULL,
  publishingOrg = NULL,
  publishingCountry = NULL,
  year = NULL,
  basisOfRecord = NULL,
  ...
)
}
\arguments{
\item{source}{(character) Either \code{density} for fast, precalculated tiles,
or \code{adhoc} for any search. Default: \code{density}}

\item{x}{(integer) the column. Default: 0}

\item{y}{(integer) the row. Default: 0}

\item{z}{(integer) the zoom. Default: 0}

\item{srs}{(character) Spatial reference system for the output (input srs for mvt
from GBIF is always \code{EPSG:3857}). One of:
\itemize{
\item \code{EPSG:3857} (Web Mercator)
\item \code{EPSG:4326} (WGS84 plate care?)
\item \code{EPSG:3575} (Arctic LAEA on 10 degrees E)
\item \code{EPSG:3031} (Antarctic stereographic)
}}

\item{bin}{(character) \code{square} or \code{hex} to aggregate occurrence counts into
squares or hexagons. Points by default. optional}

\item{hexPerTile}{(integer) sets the size of the hexagons
(the number horizontally across a tile). optional}

\item{squareSize}{(integer) sets the size of the squares. Choose a factor
of 4096 so they tessalate correctly: probably from 8, 16, 32, 64, 128,
256, 512. optional}

\item{style}{(character) for raster tiles, choose from the available styles.
Defaults to classic.point. optional. THESE DON'T WORK YET.}

\item{taxonKey}{(integer/numeric/character) search by taxon key, can only
supply 1. optional}

\item{datasetKey}{(character) search by taxon key, can only supply 1.
optional}

\item{country}{(character) search by taxon key, can only supply 1.
optional}

\item{publishingOrg}{(character) search by taxon key, can only supply 1.
optional}

\item{publishingCountry}{(character) search by taxon key, can only
supply 1. optional}

\item{year}{(integer) integer that limits the search to a certain year or,
if passing a vector of integers, multiple years, for example
\code{1984} or \code{c(2016, 2017, 2018)} or \code{2010:2015} (years 2010 to 2015). optional}

\item{basisOfRecord}{(character) one or more basis of record states to
include records with that basis of record. The full list is: \code{c("OBSERVATION", "HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "MATERIAL_SAMPLE", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN", "LIVING_SPECIMEN", "LITERATURE", "UNKNOWN")}. optional}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
an sf object
}
\description{
This function is a wrapper for the GBIF mapping api version 2.0.
The mapping API is a web map tile service making it straightforward to
visualize GBIF content on interactive maps, and overlay content from other
sources. It returns maps vector tiles with number of
GBIF records per area unit that can be used in a variety of ways, for example
in interactive leaflet web maps. Map details are specified by a number of
query parameters, some of them optional. Full documentation of the GBIF
mapping api can be found at https://www.gbif.org/developer/maps
}
\details{
This function uses the arguments passed on to generate a query
to the GBIF web map API. The API returns a web tile object as png that is
read and converted into an R raster object. The break values or nbreaks
generate a custom colour palette for the web tile, with each bin
corresponding to one grey value. After retrieval, the raster is reclassified
to the actual break values. This is a somewhat hacky but nonetheless
functional solution in the absence of a GBIF raster API implementation.

We add extent and set the projection for the output. You can reproject
after retrieving the output.
}
\examples{
\dontrun{
if (
 requireNamespace("sf", quietly = TRUE) &&
 requireNamespace("protolite", quietly = TRUE)
) {
  x <- mvt_fetch(taxonKey = 2480498, year = 2007:2011)
  x
  
  # gives an sf object
  class(x)
  
  # different srs
  ## 3857
  y <- mvt_fetch(taxonKey = 2480498, year = 2010, srs = "EPSG:3857")
  y
  ## 3031
  z <- mvt_fetch(taxonKey = 2480498, year = 2010, srs = "EPSG:3031", verbose = TRUE)
  z
  # 3575
  z <- mvt_fetch(taxonKey = 2480498, year = 2010, srs = "EPSG:3575")
  z

  # bin
  x <- mvt_fetch(taxonKey = 212, year = 1998, bin = "hex",
     hexPerTile = 30, style = "classic-noborder.poly")
  x

  # query with basisOfRecord
  mvt_fetch(taxonKey = 2480498, year = 2010,
    basisOfRecord = "HUMAN_OBSERVATION")
  mvt_fetch(taxonKey = 2480498, year = 2010,
    basisOfRecord = c("HUMAN_OBSERVATION", "LIVING_SPECIMEN"))
 }
}
}
\references{
https://www.gbif.org/developer/maps
}
\seealso{
\code{\link[=map_fetch]{map_fetch()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_download_meta.R
\name{occ_download_meta}
\alias{occ_download_meta}
\title{Retrieves the occurrence download metadata by its unique key.}
\usage{
occ_download_meta(key, curlopts = list())
}
\arguments{
\item{key}{A key generated from a request, like that from
\code{occ_download}}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
an object of class \code{occ_download_meta}, a list with slots for
the download key, the DOI assigned to the download, license link,
the request details you sent in the \code{occ_download()} request,
and metadata about the size and date/time of the request
}
\description{
Retrieves the occurrence download metadata by its unique key.
}
\note{
see \link{downloads} for an overview of GBIF downloads methods
}
\examples{
\dontrun{
occ_download_meta(key="0003983-140910143529206")
occ_download_meta("0000066-140928181241064")
}
}
\seealso{
Other downloads: 
\code{\link{download_predicate_dsl}},
\code{\link{occ_download_cached}()},
\code{\link{occ_download_cancel}()},
\code{\link{occ_download_dataset_activity}()},
\code{\link{occ_download_datasets}()},
\code{\link{occ_download_get}()},
\code{\link{occ_download_import}()},
\code{\link{occ_download_list}()},
\code{\link{occ_download_queue}()},
\code{\link{occ_download_wait}()},
\code{\link{occ_download}()}
}
\concept{downloads}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset_search.r
\name{dataset_search}
\alias{dataset_search}
\title{Search datasets in GBIF.}
\usage{
dataset_search(
  query = NULL,
  country = NULL,
  type = NULL,
  keyword = NULL,
  publishingOrg = NULL,
  hostingOrg = NULL,
  publishingCountry = NULL,
  decade = NULL,
  facet = NULL,
  facetMincount = NULL,
  facetMultiselect = NULL,
  limit = 100,
  start = NULL,
  pretty = FALSE,
  return = NULL,
  curlopts = list()
)
}
\arguments{
\item{query}{Query term(s) for full text search.  The value for this
parameter can be a simple word or a phrase. Wildcards can be added to the
simple word parameters only, e.g. \code{q=*puma*}}

\item{country}{NOT YET IMPLEMENTED. Filters by country as given in
isocodes$gbif_name, e.g. \code{country=CANADA}}

\item{type}{Type of dataset, options include occurrene, metadata, checklist,
sampling_event
(http://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/DatasetType.html)}

\item{keyword}{Keyword to search by. Datasets can be tagged by keywords,
which you can search on. The search is done on the merged collection of
tags, the dataset keywordCollections and temporalCoverages.}

\item{publishingOrg}{Publishing organization. A uuid string. See
\code{\link{organizations}}}

\item{hostingOrg}{Hosting organization. A uuid string. See
\code{\link{organizations}}}

\item{publishingCountry}{Publishing country. See options at
isocodes$gbif_name}

\item{decade}{Decade, e.g., 1980. Filters datasets by their temporal coverage
broken down to decades. Decades are given as a full year, e.g. 1880, 1960,
2000, etc, and will return datasets wholly contained in the decade as well
as those that cover the entire decade or more. Facet by decade to get the
break down, e.g. /search?facet=DECADE&facet_only=true (see example below)}

\item{facet}{A list of facet names used to retrieve the 100 most frequent
values for a field. Allowed facets are: datasetKey, highertaxonKey, rank,
status, extinct, habitat, and nameType. Additionally threat and
nomenclaturalStatus are legal values but not yet implemented, so data will
not yet be returned for them.}

\item{facetMincount}{Used in combination with the facet parameter. Set
\code{facetMincount={#}} to exclude facets with a count less than {#}}

\item{facetMultiselect}{Used in combination with the facet parameter. Set
facetMultiselect=true to still return counts for values that are not
currently filtered}

\item{limit}{Number of records to return. Default: 100. Maximum: 1000.}

\item{start}{Record number to start at. Default: 0. Use in combination
with \code{limit} to page through results.}

\item{pretty}{Print informative metadata using \code{\link{cat}}. Not easy
to manipulate output though.}

\item{return}{Defunct. All components are returned; index to the
one(s) you want}

\item{curlopts}{list of named curl options passed on to
\code{\link[crul]{HttpClient}}. see \code{curl::curl_options}
for curl options}
}
\value{
A data.frame, list, or message printed to console (using
\code{pretty=TRUE}).
}
\description{
This function does not search occurrence data, only metadata on the datasets
that contain occurrence data.
}
\section{Repeat parmeter inputs}{

Some parameters can tak emany inputs, and treated as 'OR' (e.g., a or b or
c). The following take many inputs:
\itemize{
\item \code{type}
\item \code{keyword}
\item \code{publishingOrg}
\item \code{hostingOrg}
\item \code{publishingCountry}
\item \code{decade}
}
}

\examples{
\dontrun{
# Gets all datasets of type "OCCURRENCE".
dataset_search(type="OCCURRENCE", limit = 10)

# Fulltext search for all datasets having the word "amsterdam" somewhere in
# its metadata (title, description, etc).
dataset_search(query="amsterdam", limit = 10)

# Limited search
dataset_search(type="OCCURRENCE", limit=2)
dataset_search(type="OCCURRENCE", limit=2, start=10)

# Return metadata in a more human readable way (hard to manipulate though)
dataset_search(type="OCCURRENCE", pretty=TRUE, limit = 10)

# Search by country code. Lookup isocodes first, and use US for United States
isocodes[agrep("UNITED", isocodes$gbif_name),]
dataset_search(country="US", limit = 10)

# Search by decade
dataset_search(decade=1980, limit = 10)

# Faceting
## just facets
dataset_search(facet="decade", facetMincount="10", limit=0)

## data and facets
dataset_search(facet="decade", facetMincount="10", limit=2)

# Some parameters accept many inputs, treated as OR
dataset_search(type = c("metadata", "checklist"))$data
dataset_search(keyword = c("fern", "algae"))$data
dataset_search(publishingOrg = c("e2e717bf-551a-4917-bdc9-4fa0f342c530",
  "90fd6680-349f-11d8-aa2d-b8a03c50a862"))$data
dataset_search(hostingOrg = c("c5f7ef70-e233-11d9-a4d6-b8a03c50a862",
  "c5e4331-7f2f-4a8d-aa56-81ece7014fc8"))$data
dataset_search(publishingCountry = c("DE", "NZ"))$data
dataset_search(decade = c(1910, 1930))$data

## curl options
dataset_search(facet="decade", facetMincount="10", limit=2,
  curlopts = list(verbose=TRUE))
}
}
\references{
\url{https://www.gbif.org/developer/registry#datasetSearch}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_issues.r
\name{occ_issues}
\alias{occ_issues}
\title{Parse and examine further GBIF occurrence issues on a dataset.}
\usage{
occ_issues(.data, ..., mutate = NULL)
}
\arguments{
\item{.data}{Output from a call to \code{\link[=occ_search]{occ_search()}}, \code{\link[=occ_data]{occ_data()}}, or
\code{\link[=occ_download_import]{occ_download_import()}}. The data from \code{occ_download_import}
is just a regular data.frame so you can pass in a data.frame to this
function, but if it doesn't have certain columns it will fail.}

\item{...}{Named parameters to only get back (e.g. cdround), or to
remove (e.g. -cdround).}

\item{mutate}{(character) One of:
\itemize{
\item \code{split} Split issues into new columns.
\item \code{expand} Expand issue abbreviated codes into descriptive names.
for downloads datasets, this is not super useful since the
issues come to you as expanded already.
\item \code{split_expand} Split into new columns, and expand issue names.
}

For split and split_expand, values in cells become y ("yes") or n ("no")}
}
\description{
Parse and examine further GBIF occurrence issues on a dataset.
}
\details{
See also the vignette \strong{Cleaning data using GBIF issues}

Note that you can also query based on issues, e.g.,
\code{occ_search(taxonKey=1, issue='DEPTH_UNLIKELY')}. However, I imagine
it's more likely that you want to search for occurrences based on a
taxonomic name, or geographic area, not based on issues, so it makes sense
to pull data down, then clean as needed using this function.

This function only affects the \code{data} element in the \code{gbif} class that is
returned from a call to \code{\link[=occ_search]{occ_search()}}. Maybe in a future version
we will remove the associated records from the \code{hierarchy} and \code{media}
elements as they are removed from the \code{data} element.

You'll notice that we sort columns to make it easier to glimpse the important
parts of your data, namely taxonomic name, taxon key, latitude and longitude,
and the issues. The columns are unchanged otherwise.
}
\examples{
\dontrun{
# what do issues mean, can print whole table
head(gbif_issues())
# or just occurrence related issues
gbif_issues()[which(gbif_issues()$type \%in\% c("occurrence")),]
# or search for matches
iss <- c('cdround','cudc','gass84','txmathi')
gbif_issues()[ gbif_issues()$code \%in\% iss, ]

# compare out data to after occ_issues use
(out <- occ_search(limit=100))
out \%>\% occ_issues(cdround)

# occ_data
(out <- occ_data(limit=100))
out \%>\% occ_issues(cdround)

# Parsing output by issue
(res <- occ_data(
  geometry='POLYGON((30.1 10.1,40 40,20 40,10 20,30.1 10.1))',
  limit = 600))

## or parse issues in various ways
### include only rows with cdround issue
gg <- res \%>\% occ_issues(cdround)
NROW(res$data)
NROW(gg$data)
head(res$data)[,c(1:5)]
head(gg$data)[,c(1:5)]

### remove data rows with certain issue classes
res \%>\% occ_issues(-cdround, -cudc)

### split issues into separate columns
res \%>\% occ_issues(mutate = "split")
res \%>\% occ_issues(-cudc, -mdatunl, mutate = "split")
res \%>\% occ_issues(gass84, mutate = "split")

### expand issues to more descriptive names
res \%>\% occ_issues(mutate = "expand")

### split and expand
res \%>\% occ_issues(mutate = "split_expand")

### split, expand, and remove an issue class
res \%>\% occ_issues(-cdround, mutate = "split_expand")

## Or you can use occ_issues without \%>\%
occ_issues(res, -cdround, mutate = "split_expand")

# from GBIF downloaded data via occ_download_* functions
res <- occ_download_get(key="0000066-140928181241064", overwrite=TRUE)
x <- occ_download_import(res)
occ_issues(x, -txmathi)
occ_issues(x, txmathi)
occ_issues(x, gass84)
occ_issues(x, zerocd)
occ_issues(x, gass84, txmathi)
occ_issues(x, mutate = "split")
occ_issues(x, -gass84, mutate = "split")
occ_issues(x, mutate = "expand")
occ_issues(x, mutate = "split_expand")

# occ_search/occ_data with many inputs - give slightly different output
# format than normal 2482598, 2498387
xyz <- occ_data(taxonKey = c(9362842, 2492483, 2435099), limit = 300)
xyz
length(xyz) # length 3
names(xyz) # matches taxonKey values passed in
occ_issues(xyz, -gass84)
occ_issues(xyz, -cdround)
occ_issues(xyz, -cdround, -gass84)
}
}
\references{
https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbif_issues_lookup.r
\name{gbif_issues_lookup}
\alias{gbif_issues_lookup}
\title{Lookup issue definitions and short codes}
\usage{
gbif_issues_lookup(issue = NULL, code = NULL)
}
\arguments{
\item{issue}{Full name of issue, e.g, CONTINENT_COUNTRY_MISMATCH}

\item{code}{An issue short code, e.g. 'ccm'}
}
\description{
Lookup issue definitions and short codes
}
\examples{
gbif_issues_lookup(issue = 'CONTINENT_COUNTRY_MISMATCH')
gbif_issues_lookup(code = 'ccm')
gbif_issues_lookup(issue = 'COORDINATE_INVALID')
gbif_issues_lookup(code = 'cdiv')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbif_issues.R
\name{gbif_issues}
\alias{gbif_issues}
\title{List all GBIF issues and their codes.}
\source{
https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html
https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/NameUsageIssue.html
}
\usage{
gbif_issues()
}
\description{
Returns a data.frame of all GBIF issues with the following columns:
\itemize{
\item \code{code}: issue short code, e.g. \code{gass84}
\item \code{code}: issue full name, e.g. \code{GEODETIC_DATUM_ASSUMED_WGS84}
\item \code{description}: issue description
\item \code{type}: issue type, either related to \code{occurrence} or \code{name}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_facet.r
\name{count_facet}
\alias{count_facet}
\title{Facetted count occurrence search.}
\usage{
count_facet(keys = NULL, by = "country", countries = 10, removezeros = FALSE)
}
\arguments{
\item{keys}{(numeric) GBIF keys, a vector. optional}

\item{by}{(character) One of georeferenced, basisOfRecord, country, or
publishingCountry. default: country}

\item{countries}{(numeric) Number of countries to facet on, or a vector of
country names. default: 10}

\item{removezeros}{(logical) remove zeros or not? default: \code{FALSE}}
}
\description{
Facetted count occurrence search.
}
\examples{
\dontrun{
# Select number of countries to facet on
count_facet(by='country', countries=3, removezeros = TRUE)
# Or, pass in country names
count_facet(by='country', countries='AR', removezeros = TRUE)

spplist <- c('Geothlypis trichas','Tiaris olivacea','Pterodroma axillaris',
             'Calidris ferruginea','Pterodroma macroptera',
             'Gallirallus australis',
             'Falco cenchroides','Telespiza cantans','Oreomystis bairdi',
             'Cistothorus palustris')
keys <- sapply(spplist,
  function(x) name_backbone(x, rank="species")$usageKey)
count_facet(keys, by='country', countries=3, removezeros = TRUE)
count_facet(keys, by='country', countries=3, removezeros = FALSE)
count_facet(by='country', countries=20, removezeros = TRUE)
count_facet(keys, by='basisOfRecord', countries=5, removezeros = TRUE)

# Pass in country names instead
countries <- isocodes$code[1:10]
count_facet(by='country', countries=countries, removezeros = TRUE)

# get occurrences by georeferenced state
## across all records
count_facet(by='georeferenced')

## by keys
count_facet(keys, by='georeferenced')

# by basisOfRecord
count_facet(by="basisOfRecord")
}
}
