ckanr
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/ckanr)](https://cranchecks.info/pkgs/ckanr)
[![R-check](https://github.com/ropensci/ckanr/workflows/R-check/badge.svg)](https://github.com/ropensci/ckanr/actions/)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ckanr?color=FAB657)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/ckanr)](https://cran.r-project.org/package=ckanr)

`ckanr` is an R client for the CKAN API.

## Description

CKAN is an open source set of tools for hosting and providing data on the web. (CKAN users could include non-profits, museums, local city/county governments, etc.).

`ckanr` allows users to interact with those CKAN websites to create, modify, and manage datasets, as well as search and download pre-existing data, and then to proceed using in R for data analysis (stats/plotting/etc.). It is meant to be as general as possible, allowing you to work with any CKAN instance.

Get started: <https://docs.ropensci.org/ckanr/>

## Installation

Stable CRAN version


```r
install.packages("ckanr")
```

Development version


```r
install.packages("remotes")
remotes::install_github("ropensci/ckanr")
```


```r
library('ckanr')
```

Note: the default base CKAN URL is set to <https://data.ontario.ca/>
Functions requiring write permissions in CKAN additionally require a privileged
CKAN API key.
You can change this using `ckanr_setup()`, or change the URL using the `url`
parameter in each function call.
To set one or both, run:


```r
ckanr_setup() # restores default CKAN url to https://data.ontario.ca/
ckanr_setup(url = "https://data.ontario.ca/")
ckanr_setup(url = "https://data.ontario.ca/", key = "my-ckan-api-key")
```

## ckanr package API

There are a suite of CKAN things (package, resource, etc.) that each have a set of functions in this package. The functions for each CKAN thing have an S3 class that is returned from most functions, and can be passed to most other functions (this also facilitates piping). The following is a list of the function groups for certain CKAN things, with the prefix for the functions that work with that thing, and the name of the S3 class:

+ Packages (aka packages) - `package_*()` - `ckan_package`
+ Resources - `resource_*()` - `ckan_resource`
+ Related - `related_*()` - `ckan_related`
+ Users - `user_*()` - `ckan_user`
+ Groups - `group_*()` - `ckan_group`
+ Tags - `tag_*()` - `ckan_tag`
+ Organizations  - `organization_*()` - `ckan_organization`
+ Groups - `group_*()` - `ckan_group`
+ Users - `user_*()` - `ckan_user`
+ Related items - `related_*()` - `ckan_related`

The S3 class objects all look very similar; for example:

```r
<CKAN Resource> 8abc92ad-7379-4fb8-bba0-549f38a26ddb
  Name: Data From Digital Portal
  Description:
  Creator/Modified: 2015-08-18T19:20:59.732601 / 2015-08-18T19:20:59.657943
  Size:
  Format: CSV
```

All classes state the type of object, have the ID to the right of the type, then have a varying set of key-value fields deemed important. This printed object is just a summary of an R list, so you can index to specific values (e.g., `result$description`). If you feel there are important fields left out of these printed summaries, let us know.

> note: Many examples are given in brief for readme brevity


## Contributors

(alphebetical)

* Scott Chamberlain
* Imanuel Costigan
* Sharla Gelfand
* Florian Mayer
* Wush Wu

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/ckanr/issues).
* License: MIT
* Get citation information for `ckanr` in R doing `citation(package = 'ckanr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
ckanr 0.6.0
===========

### NEW FEATURES

* parameter `http_method` gained in `resource_create()`, `package_update()`, and `package_patch()`; it's passed to `as.ckan_package()` internally, but does not affect the HTTP request for the main point of the function (#163) thanks @hannaboe
* gains new function `organization_purge()` to purge an organization (which requires sysadmin) (#166) thanks @nicholsn

### MINOR IMPROVEMENTS

* update URLs for known CKAN instances behind the `servers()` function (#162) (#167) (#170)
* .Rbuildignore README.md and vignettes (#171)
* `extras` now passed in HTTP request in `package_create()` as a top level part of the request body rather than as a named `extras` element (#158) thanks @galaH
* change in `ckan_fetch()`: now when a zip file has a subdirectory an `NA` is returned rather than `character(0)` (#164)
* change in the `...` parameter in `ckan_fetch()`: was used to pass through curl options to the http request but now is used to pass through additional parameters to either `read.csv`, `xml2::read_xml`, `jsonlite::fromJSON`, `sf::st_read()`, `read.table()`, or `readxl::read_excel` (#165)
* updated docs for `package_create()` and `group_create()` - change `groups` parameter description to explain what kind of input is expected (#168)


ckanr 0.5.0
===========

### NEW FEATURES

* `package_create()` gains parameter `private` (boolean) (#145)
* add support for resource extras. `resource_create()` and `resource_update()` gain new parameter `extras`, while `resource_patch()` function doesn't change but gains an example of adding an extra (#149) (#150) thanks @nicholsn
* `package_patch()` gains `extras` parameter (#94) see also (#147)

### MINOR IMPROVEMENTS

* replace httr with crul throughout package (#86) (#151)
* use markdown for docs (#148)
* `package_patch()`, `package_show()`, `package_activity_list()`, `package_delete()`, `package_update()`, and `related_create()` now pass on `key` parameter value to `as.ckan_package` internally; same for `resource_create()` and `resource_show()`, but passed to `as.ckan_resource()`  (#145) (#146)
* `servers()` gains two additional CKAN urls (#155)
* `package_show()` called `as.ckan_package()` within it, which itself calls `package_show()` - fixed now  (#127)

### BUG FIXES

* fix for `resource_search()` and `tag_search`: both were not allowing a query to be more than length 1 (#153)
* fix for `print.ckan_package`: wasn't handling well results from `package_search()` that had a named list of locale specific results (#152)


ckanr 0.4.0
===========

### NEW FEATURES

* `ckan_fetch()` gains parameter `key` for a CKAN API key; if given the API key is now included in the request headers (#133) see also (#122) by @sharlagelfand
* `ckan_fetch()` gains ability to read xls/xlsx files with multiple sheets (#135) by @sharlagelfand

### MINOR IMPROVEMENTS

* `ckan_fetch()` now sets `stringsAsFactors = FALSE` when reading data (#141) (#142) thanks @LVG77 @sharlagelfand
* in `ckan_fetch()`, use `basename(x)` instead of `gsub(paste0(tempdir(), "/"), "", x)`, to get file path (#140) by @sharlagelfand
* in `package_search()` handle better cases where the CKAN version can not be determined (#139) && fix logic for when `default_schema` and `include_private` parameters are included based on the CKAN version (#137) by @sharlagelfand
* improve `ckan_fetch()`: old behavior of the fxn with zip files was that it only worked if the zip file contained shp files; works more generally now, e.g., a zip file containing a csv file (#132) by @sharlagelfand
* fix `ckan_fetch()` examples that weren't working (#134) by @sharlagelfand
* fix to parsing CKAN version numbers, new internal fxn `parse_version_number()` - now properly parses CKAN version numbers that include patch and dev versions (#136) by @sharlagelfand


ckanr 0.3.0
===========

### NEW FEATURES

* new package author Sharla Gelfand !!!
* new functions for users: `user_create()` and `user_delete()` (#82)
* `package_show()` gains `key` parameter to pass an API key (#97)
* `package_search()` gains new parameters: `include_drafts`, `include_private`, `use_default_schema`, and `facet.mincount` (#107)
* function `fetch()` changed to `ckan_fetch()`
* gains function `organization_delet()` to delete an organization (#83)
* gains function `ckan_version()` to get version info for a CKAN instance
* gains methods for creating a CKAN remote instance as a dplyr backend: gains `src_ckan()` and it's s3 methods `tbl` and `src_tbls`, `sql_translate_env`. in addition gains the S3 methods `db_begin`, `db_explain`, `db_has_table`, `db_insert_into`, `db_query_fields`, `db_query_rows`

### MINOR IMPROVEMENTS

* fix some tests (#62)
* fix to `ds_create()` to properly format body with json data (#85) thanks @mattfullerton
* tests added for `ckan_fetch()` (#118) thanks @sharlagelfand
* `ckan_fetch()` gains `format` parameter if the user knows the file format (useful when the file format can not be guessed) (#117) thanks @sharlagelfand
* `ckan_fetch` gain support for handling geojson (#123) thanks @sharlagelfand
* `ckan_fetch` was writing to current working directory in some cases - fixed to writing to temp files and cleaning up (#125) (#128) (#129) thanks @sharlagelfand
* add USDA CKAN instance to the `servers()` function (#68)
* `ds_create_dataset()` marked as deprecated; see `recourse_create()` instead (#80) (via #79)
* removed the internal `stop()` call in `tag_create()`: now can be used, though haven't been able to test this function as you need to be a sysadmin to use it (#81)
* `ds_search()`: code spacing fixes (#69)
* `resource_update()` gains more examples and tests (#66)
* CKAN API key standardization: `key` parameter now in all fxns that make http requests - and reordering of `url` and `key` params in that order across all functions (#122) (#124)
* repair ORCID links in DESCRIPTION file (#124) by Florian

### BUG FIXES

* fix to `resource_create()`: `upload` param was inappropriately a required param (#75) thanks @mingbogo
* fixes to `resource_update()`: date sent in `last_modified` in request body needed to be converted to character (#96) (thanks @jasonajones73); and the date format needed fixing (#119) (thanks @florianm)
* fix to `ckan_fetch()` - use `sf` instead of `maptools`; in addition `ckan_fetch` can now parse xlsx files in addition to xls files;  (#114) (#115) thanks @sharlagelfand
* fix to `package_search()`: this route fails if parameters that did not exist in the CKAN instance are given; internally remove parameters as needed from query params by pinging the CKAN instance for its version (#120)


ckanr 0.1.0
===========

### NEW FEATURES

* Releasd to CRAN.
## Test environments

* local macOS install, R 4.0.3 Patched
* ubuntu 14.04 (on github actions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 2 reverse dependencies. Summary at <https://github.com/ropensci/rgbif/blob/master/revdep/README.md>. No problems were found related to this package.

--------

This version adds a new function, a few functions gain a new parameter, and many minor fixes.

Thanks very much,
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

* Submit an issue on the [Issues page](https://github.com/ropensci/ckanr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/ckanr.git`
* Make sure to track progress upstream (i.e., on our version of `ckanr` at `ropensci/ckanr`) by doing `git remote add upstream https://github.com/ropensci/ckanr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/ckanr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Email

Don't send email. Open an issue instead.
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If ckanr threw an error, please also post the traceback() of the error. Otherwise, delete all this and proceed :) -->

<details>
<summary><strong>Session Info and Traceback</strong></summary>

<!-- SessionInfo captures the exact circumstances under which the error happened: OS, installed packages. -->
```r
# devtools::session_info()  # or
utils::sessionInfo()

```

<!-- If ckanr threw an error, run traceback() immediately after the error and paste the results here. -->
```r
traceback()

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.3 Patched (2021-02-02 r79929) |
|os       |macOS Catalina 10.15.7                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-02-03                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|ckanr   |0.5.0 |0.6.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*