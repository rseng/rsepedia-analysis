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

*Wow, no problems at all. :)**Wow, no problems at all. :)*ckanr
=====

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("ckanr")
```

Development version

```{r eval=FALSE}
install.packages("remotes")
remotes::install_github("ropensci/ckanr")
```

```{r}
library('ckanr')
```

Note: the default base CKAN URL is set to <https://data.ontario.ca/>
Functions requiring write permissions in CKAN additionally require a privileged
CKAN API key.
You can change this using `ckanr_setup()`, or change the URL using the `url`
parameter in each function call.
To set one or both, run:

```{r}
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
---
title: "ckanr vignette"
author: "Scott Chamberlain"
date: "2020-07-29"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{ckanr vignette}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



## Install

Stable version from CRAN


```r
install.packages("ckanr")
```

Development version from GitHub


```r
remotes::install_github("ropensci/ckanr")
```


```r
library("ckanr")
```

Note: the default base CKAN URL is set to [http://data.techno-science.ca/](http://data.techno-science.ca/). You can change this using `ckanr_setup()`, or change the URL using the `url` 
parameter in each function call.

To set one or both, run:


```r
# restores default CKAN url to http://data.techno-science.ca/
ckanr_setup() 
# Just set url
ckanr_setup(url = "http://data.techno-science.ca/")
# set url and key
ckanr_setup(url = "http://data.techno-science.ca/", key = "my-ckan-api-key")
```

## Changes


```r
changes(limit = 2, as = "table")[, 1:4]
#>                                user_id                  timestamp
#> 1 27778230-2e90-4818-9f00-bbf778c8fa09 2016-06-14T21:31:28.306231
#> 2 27778230-2e90-4818-9f00-bbf778c8fa09 2016-06-14T21:30:26.594125
#>                              object_id                          revision_id
#> 1 99f457c9-ea24-48a1-87be-b52385825b6a e2d9463d-e97c-48f5-a816-7fe26ee60dcd
#> 2 99f457c9-ea24-48a1-87be-b52385825b6a 9d846213-1389-4dab-bfe5-77dd3256995a
```

## List datasets


```r
package_list(as = "table")
#>  [1] "artifact-data-agriculture"                                  
#>  [2] "artifact-data-aviation"                                     
#>  [3] "artifact-data-bookbinding"                                  
#>  [4] "artifact-data-chemistry"                                    
#>  [5] "artifact-data-communications"                               
#>  [6] "artifact-data-computing-technology"                         
#>  [7] "artifact-data-domestic-technology"                          
#>  [8] "artifact-data-energy-electric"                              
#>  [9] "artifact-data-exploration-and-survey"                       
#> [10] "artifact-data-fisheries"                                    
...
```

## List tags


```r
tag_list('aviation', as = 'table')
#>   vocabulary_id                     display_name
#> 1            NA                         Aviation
#> 2            NA Canada Aviation and Space Museum
#>                                     id                             name
#> 1 cc1db2db-b08b-4888-897f-a17eade2461b                         Aviation
#> 2 8d05a650-bc7b-4b89-bcc8-c10177e60119 Canada Aviation and Space Museum
```

## Show tags

Subset for readme brevity


```r
tag_show('Aviation')
#> <CKAN Tag> cc1db2db-b08b-4888-897f-a17eade2461b 
#>   Name: Aviation
#>   Display name: Aviation
#>   Vocabulary id: 
#>   No. Packages: 2
#>   Packages (up to 5): artifact-data-aviation, cstmc-smstc-artifacts-artefact
```

## List groups


```r
group_list(as = 'table')[, 1:3]
#>                         display_name description
#> 1                     Communications            
#> 2 Domestic and Industrial Technology            
#> 3                         Everything            
#> 4                           Location            
#> 5                          Resources            
#> 6         Scientific Instrumentation            
#> 7                     Transportation            
#>                                title
#> 1                     Communications
#> 2 Domestic and Industrial Technology
#> 3                         Everything
#> 4                           Location
#> 5                          Resources
#> 6         Scientific Instrumentation
#> 7                     Transportation
```

## Show groups

Subset for readme brevity


```r
group_show('communications', as = 'table')$users
#>   openid about capacity     name                    created
#> 1     NA  <NA>    admin     marc 2014-10-24T14:44:29.885262
#> 2     NA          admin sepandar 2014-10-23T19:40:42.056418
#>                         email_hash sysadmin
#> 1 a32002c960476614370a16e9fb81f436    FALSE
#> 2 10b930a228afd1da2647d62e70b71bf8     TRUE
#>   activity_streams_email_notifications  state number_of_edits
#> 1                                FALSE active             516
#> 2                                FALSE active              44
#>   number_administered_packages display_name fullname
#> 1                           40         marc     <NA>
#> 2                            1     sepandar         
#>                                     id
#> 1 27778230-2e90-4818-9f00-bbf778c8fa09
#> 2 b50449ea-1dcc-4d52-b620-fc95bf56034b
```

## Show a package


```r
package_show('34d60b13-1fd5-430e-b0ec-c8bc7f4841cf', as = 'table')$resources[, 1:10]
#>                      resource_group_id cache_last_updated
#> 1 ea8533d9-cdc6-4e0e-97b9-894e06d50b92                 NA
#> 2 ea8533d9-cdc6-4e0e-97b9-894e06d50b92                 NA
#> 3 ea8533d9-cdc6-4e0e-97b9-894e06d50b92                 NA
#> 4 ea8533d9-cdc6-4e0e-97b9-894e06d50b92                 NA
#> 5 ea8533d9-cdc6-4e0e-97b9-894e06d50b92                 NA
#>           revision_timestamp webstore_last_updated
#> 1 2016-06-13T20:05:16.818800                    NA
#> 2 2014-11-04T02:59:50.567068                    NA
#> 3 2014-11-05T21:23:58.533397                    NA
#> 4 2014-11-05T21:25:16.848423                    NA
#> 5 2016-06-13T20:06:50.013746                    NA
#>                                     id size  state hash
#> 1 be2b0af8-24a8-4a55-8b30-89f5459b713a   NA active     
#> 2 7d65910e-4bdc-4f06-a213-e24e36762767   NA active     
#> 3 97622ad7-1507-4f6a-8acb-14e826447389   NA active     
#> 4 7a72498a-c49c-4e84-8b10-58991de10df6   NA active     
#> 5 7e2cb5de-550d-41a8-ab9d-b2ec35b6671a   NA active     
#>                                    description format
#> 1                                  XML Dataset    XML
#> 2 Data dictionary for CSTMC artifact datasets.    XLS
#> 3       Tips for using the artifacts datasets.   .php
#> 4       Tips for using the artifacts datasets.   .php
#> 5                          Jeux de données XML    XML
```

## Search for packages


```r
out <- package_search(q = '*:*', rows = 2, as = "table")$results
out[, !names(out) %in% 'resources'][, 1:10]
#>                      license_title maintainer relationships_as_object private
#> 1 Open Government Licence - Canada                               NULL   FALSE
#> 2 Open Government Licence - Canada                               NULL   FALSE
#>   maintainer_email         revision_timestamp
#> 1                  2014-10-28T21:27:57.475091
#> 2                  2014-10-28T20:40:55.803602
#>                                     id           metadata_created
#> 1 99f457c9-ea24-48a1-87be-b52385825b6a 2014-10-24T17:39:06.411039
#> 2 443cb020-f2ae-48b1-be67-90df1abd298e 2014-10-28T20:39:23.561940
#>            metadata_modified author
#> 1 2016-06-14T21:31:27.983485       
#> 2 2016-06-14T18:59:17.786219
```

## Search for resources


```r
resource_search(q = 'name:data', limit = 2, as = 'table')
#> $count
#> [1] 74
#> 
#> $results
#>                      resource_group_id cache_last_updated webstore_last_updated
#> 1 01a82e52-01bf-4a9c-9b45-c4f9b92529fa                 NA                    NA
#> 2 01a82e52-01bf-4a9c-9b45-c4f9b92529fa                 NA                    NA
#>                                     id size  state last_modified hash
#> 1 e179e910-27fb-44f4-a627-99822af49ffa   NA active            NA     
#> 2 ba84e8b7-b388-4d2a-873a-7b107eb7f135   NA active            NA     
#>                                    description format mimetype_inner url_type
#> 1                                  XML Dataset    XML             NA       NA
#> 2 Data dictionary for CSTMC artifact datasets.    XLS             NA       NA
#>   mimetype cache_url                                         name
#> 1       NA        NA Artifact Data - Exploration and Survey (XML)
#> 2       NA        NA                              Data Dictionary
#>                      created
#> 1 2014-10-28T15:50:35.374303
#> 2 2014-11-03T18:01:02.094210
#>                                                                                                                                                    url
#> 1              http://source.techno-science.ca/datasets-donn%C3%A9es/artifacts-artefacts/groups-groupes/exploration-and-survey-exploration-et-leve.xml
#> 2 http://source.techno-science.ca/datasets-donn%C3%A9es/artifacts-artefacts/cstmc-artifact-data-dictionary-dictionnaire-de-donnees-artefacts-smstc.xls
#>   webstore_url position                          revision_id resource_type
#> 1           NA        0 a22e6741-3e89-4db0-a802-ba594b1c1fad            NA
#> 2           NA        1 da1f8585-521d-47ef-8ead-7832474a3421            NA
```

## ckanr's dplyr interface
`ckanr` implements a `dplyr` SQL interface to CKAN's datastore. 
You can access any resource in the datastore directly using only the CKAN
resource ID.

Note: this will only work for resources which were uploaded successfully to the 
datastore - they will show the green "Data API" button in CKAN.


```r
ckan <- ckanr::src_ckan("https://my.ckan.org/")
res_id <- "my-ckan-resource-id"
dplyr::tbl(src = ckan$con, from = res_id) %>% as_tibble(.)
```

## Example of using a different CKAN API

See `ckanr::servers()` for a list of CKAN servers. Ther are 127 as of 2020-07-29.

### The UK Natural History Museum

Website: <https://data.nhm.ac.uk/>

List datasets


```r
ckanr_setup(url = "https://data.nhm.ac.uk")
package_list(as = "table")
#>  [1] "3d-cetacean-scanning"                                                                           
#>  [2] "3d-laser-scan-of-nhm-pv-r-9372-palaeosauropus-sp"                                               
#>  [3] "abyssline"                                                                                      
#>  [4] "african-spiny-solanum"                                                                          
#>  [5] "aleyrodoidea-slide-collection"                                                                  
#>  [6] "alice-test-images"                                                                              
#>  [7] "alignments-of-co1-nd1-and-16s-rrna-for-the-land-snail-corilla"                                  
#>  [8] "al-sabouni-et-al-reproducibility"                                                               
#>  [9] "american-phlebotominae-nhm"                                                                     
#> [10] "an-influence-of-environmental-variability-on-insects-wing-shape-a-case-study-of-british-odonata"
...
```

Tags

_list_


```r
head(tag_list(as = "table"))
#>   vocabulary_id display_name                                   id         name
#> 1            NA          16S a2fabd81-89df-4520-96d2-d5cd7bdcf094          16S
#> 2            NA          18S 1f9670fe-c4c5-40bb-a02a-8889a11788e7          18S
#> 3            NA          28S 71e10e81-c643-4b67-af02-4f9516b4238b          28S
#> 4            NA           3d 098ccfee-a0fe-451d-b748-b617086d146c           3d
#> 5            NA           3D b26c5cce-dc41-40c6-9cbf-f966aa75d458           3D
#> 6            NA 3D modelling 119868d7-1753-4c35-8641-14681dc472a7 3D modelling
```

_show_


```r
tag_show('arthropods', as = 'table')
#> $vocabulary_id
#> NULL
#> 
#> $display_name
#> [1] "arthropods"
#> 
#> $id
#> [1] "f9245868-f4cb-4c85-a59d-11692db19e86"
#> 
#> $name
#> [1] "arthropods"
```

Packages

_search_


```r
out <- package_search(q = '*:*', rows = 2, as = 'table')
out$results[, 1:10]
#>                  license_title maintainer contributors relationships_as_object
#> 1 Creative Commons Attribution         NA                                 NULL
#> 2                      CC0-1.0         NA         <NA>                    NULL
#>   private maintainer_email num_tags            affiliation update_frequency
#> 1   FALSE               NA        1 Natural History Museum                 
#> 2   FALSE               NA        1                   <NA>           weekly
#>                                     id
#> 1 d68e20f4-a56d-4a8a-a8d7-dc478ba64c76
#> 2 56e711e6-c847-4f99-915a-6894bb5c5dea
```

_show_


```r
package_show(id = "56e711e6-c847-4f99-915a-6894bb5c5dea", as = "table")
#> $domain
#> [1] "data.nhm.ac.uk"
#> 
#> $license_title
#> [1] "CC0-1.0"
#> 
#> $maintainer
#> NULL
#> 
#> $relationships_as_object
...
```

### The National Geothermal Data System

Website: <http://geothermaldata.org/>


```r
ckanr_setup("http://search.geothermaldata.org")
x <- package_search(q = '*:*', rows = 1)
x$results
#> [[1]]
#> <CKAN Package> 787986f3-fd65-4e99-b28e-077c69933c76 
#>   Title: Hawthorne Nevada Deep Direct-Use Feasibility Study - Data Used for Geothermal Resource Conceptual Modeling and Power Capacity Estimates Hawthorne_fracture_data_HWAAD-2A_HWAAD-3.xlsx
#>   Creator/Modified: 2020-05-07T23:35:47.237011 / 2020-05-07T23:35:47.328456
#>   Resources (up to 5): Hawthorne_fracture_data_HWAAD-2A_HWAAD-3.xlsx
#>   Tags (up to 5): 2-meter temperatures, DDU, HWAAD-2, HWAAD-2A, HWAAD-3
#>   Groups (up to 5):
NA
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_activity_list.R
\name{user_activity_list}
\alias{user_activity_list}
\title{Return a list of a user's activities}
\usage{
user_activity_list(
  id,
  offset = 0,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) User identifier.}

\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a list of a user's activities
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/",
key = getOption("ckan_demo_key"))

# list package activity
user_activity_list('sckottie')

# input a ckan_user object
(x <- user_show('sckottie'))
user_activity_list(x)

# output different data formats
user_activity_list(x, as = "table")
user_activity_list(x, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ckan_group.R
\name{as.ckan_group}
\alias{as.ckan_group}
\alias{is.ckan_group}
\title{ckan_group class helpers}
\usage{
as.ckan_group(x, ...)

is.ckan_group(x)
}
\arguments{
\item{x}{Variety of things, character, list, or ckan_group class object}

\item{...}{Further args passed on to \code{\link[=group_show]{group_show()}} if character given}
}
\description{
ckan_group class helpers
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

(grps <- group_list())
grps[[3]]

# create item class from only an item ID
as.ckan_group(grps[[3]]$id)

# gives back itself
(x <- as.ckan_group(grps[[3]]$id))
as.ckan_group(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_create.R
\name{package_create}
\alias{package_create}
\title{Create a package}
\usage{
package_create(
  name = NULL,
  title = NULL,
  private = FALSE,
  author = NULL,
  author_email = NULL,
  maintainer = NULL,
  maintainer_email = NULL,
  license_id = NULL,
  notes = NULL,
  package_url = NULL,
  version = NULL,
  state = "active",
  type = NULL,
  resources = NULL,
  tags = NULL,
  extras = NULL,
  relationships_as_object = NULL,
  relationships_as_subject = NULL,
  groups = NULL,
  owner_org = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{name}{(character) the name of the new dataset, must be between 2 and
100 characters long and contain only lowercase alphanumeric characters, -
and _, e.g. 'warandpeace'}

\item{title}{(character) the title of the dataset (optional, default: same
as name)}

\item{private}{(logical) whether the dataset should be private (optional,
default: \code{FALSE}), requires a value for \code{owner_org} if \code{TRUE}}

\item{author}{(character) the name of the dataset's author (optional)}

\item{author_email}{(character) the email address of the dataset's author
(optional)}

\item{maintainer}{(character) the name of the dataset's maintainer
(optional)}

\item{maintainer_email}{(character) the email address of the dataset's
maintainer (optional)}

\item{license_id}{(license id string) - the id of the dataset's license,
see license_list() for available values (optional)}

\item{notes}{(character) a description of the dataset (optional)}

\item{package_url}{(character) a URL for the dataset's source (optional)}

\item{version}{(string, no longer than 100 characters) - (optional)}

\item{state}{(character) the current state of the dataset, e.g. 'active'
or 'deleted', only active datasets show up in search results and other
lists of datasets, this parameter will be ignored if you are not authorized
to change the state of the dataset (optional, default: 'active')}

\item{type}{(character) the type of the dataset (optional), IDatasetForm
plugins associate themselves with different dataset types and provide custom
dataset handling behaviour for these types}

\item{resources}{(list of resource dictionaries) - the dataset's resources,
see \code{\link[=resource_create]{resource_create()}} for the format of resource dictionaries (optional)}

\item{tags}{(list of tag dictionaries) - the dataset's tags, see
\code{\link[=tag_create]{tag_create()}} for the format of tag dictionaries (optional)}

\item{extras}{(list of dataset extra dictionaries) - the dataset's extras
(optional), extras are arbitrary (key: value) metadata items that can be
added to datasets, each extra dictionary should have keys 'key' (a string),
'value' (a string)}

\item{relationships_as_object}{(list of relationship dictionaries) - see
\code{package_relationship_create} for the format of relationship dictionaries
(optional)}

\item{relationships_as_subject}{(list of relationship dictionaries) - see
\code{package_relationship_create} for the format of relationship dictionaries
(optional)}

\item{groups}{(data.frame) the groups to which the dataset
belongs, each row should have one or more of the
following columns which identify an existing group: 'id' (the id of the group,
string), or 'name' (the name of the group, string), to see which groups
exist call \code{\link[=group_list]{group_list()}}. see example (optional)}

\item{owner_org}{(character) the id of the dataset's owning organization,
see \code{\link[=organization_list]{organization_list()}} or \code{organization_list_for_user} for available
values (optional)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Create a package
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create a package
## Example 1
(res <- package_create("foobar4", author="Jane Doe"))
res$author

## Example 2 - create package, add a resource
(res <- package_create("helloworld", author="Jane DOe"))

# include a group
# package_create("brownbear", groups = data.frame(id = "some-id"))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/servers.R
\name{servers}
\alias{servers}
\title{CKAN server URLS and other info}
\usage{
servers()
}
\description{
CKAN server URLS and other info
}
\details{
Comes from the links at https://ckan.org/about/instances/

There were a number of other URLs for CKAN instances in the CKAN URL
above, but some sites are now gone completely, or if they do exist,
I can't figure out how to get access to the CKAN API on their instance.
}
\examples{
\dontrun{
servers()
ckan_info(servers()[5])

# what version is each CKAN server running
out <- lapply(servers()[1:6], function(w) {
  cat(w, sep='\n') 
  ckan_info(w)
})
vapply(out, "[[", "", "ckan_version")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ckan_related.R
\name{as.ckan_related}
\alias{as.ckan_related}
\alias{is.ckan_related}
\title{ckan_related class helpers}
\usage{
as.ckan_related(x, ...)

is.ckan_related(x)
}
\arguments{
\item{x}{Variety of things, character, list, or ckan_related class object}

\item{...}{Further args passed on to \code{\link[=related_show]{related_show()}} if character given}
}
\description{
ckan_related class helpers
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
key = getOption("ckan_demo_key"))

(x <- package_create("foobbbbbarrrrr") \%>\%
   related_create(title = "my resource",
                  type = "visualization"))

# create item class from only an item ID
as.ckan_related(x$id)

# gives back itself
(x <- as.ckan_related(x$id))
as.ckan_related(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/related_show.R
\name{related_show}
\alias{related_show}
\title{Show a related item}
\usage{
related_show(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Related item identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Show a related item
}
\details{
By default the help and success slots are dropped, and only the
result slot is returned. You can request raw json with \code{as = 'json'}
then parse yourself to get the help slot.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create a package and a related item
res <- package_create("hello-pluto2") \%>\%
   related_create(title = "my resource",
                  type = "visualization")

# show the related item
related_show(res)
related_show(res$id)

# get data back in different formats
related_show(res, as = 'json')
related_show(res, as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_show.R
\name{user_show}
\alias{user_show}
\title{Show a user.}
\usage{
user_show(
  id,
  user_obj = NULL,
  include_datasets = FALSE,
  include_num_followers = FALSE,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Package identifier.}

\item{user_obj}{(user dictionary) The user dictionary of the user (optional)}

\item{include_datasets}{(logical) Include a list of datasets the user has
created. If it is the same user or a sysadmin requesting, it includes
datasets that are draft or private. (optional, default:False, limit:50)}

\item{include_num_followers}{(logical) Include the number of followers
the user has (optional, default:False)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Show a user.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# show user
user_show('sckottie')

# include datasets
user_show('sckottie', include_datasets = TRUE)

# include datasets
user_show('sckottie', include_num_followers = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/related_create.R
\name{related_create}
\alias{related_create}
\title{Create a related item}
\usage{
related_create(
  id,
  title,
  type,
  description = NULL,
  related_id = NULL,
  related_url = NULL,
  image_url = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) id of package that the related item should be added to.
This should be an alphanumeric string. Required.}

\item{title}{(character) Title of the related item. Required.}

\item{type}{(character) The type of the related item. One of API,
application, idea, news article, paper, post or visualization. Required.}

\item{description}{(character) description (optional). Optional}

\item{related_id}{(character) An id to assign to the related item. If blank,
an ID will be assigned for you. Optional}

\item{related_url}{(character) A url to associated with the related item.
Optional}

\item{image_url}{(character) A url to associated image. Optional}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Create a related item
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("hello-mars"))

# create a related item
related_create(res, title = "asdfdaf", type = "idea")

# pipe operations together
package_create("foobbbbbarrrr") \%>\%
   related_create(title = "my resource",
                  type = "visualization")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/changes.R
\name{changes}
\alias{changes}
\title{Get an activity stream of recently changed datasets on a site.}
\usage{
changes(
  offset = 0,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Get an activity stream of recently changed datasets on a site.
}
\examples{
\dontrun{
changes()
changes(as = 'json')
changes(as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_update.R
\name{package_update}
\alias{package_update}
\title{Update a package}
\usage{
package_update(
  x,
  id,
  http_method = "GET",
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{x}{(list) A list with key-value pairs}

\item{id}{(character) Package identifier}

\item{http_method}{(character) which HTTP method (verb) to use; one of
"GET" or "POST". Default: "GET"}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Update a package
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# Create a package
(pkg <- package_create("hello-world11", author="Jane Doe"))

# Next show the package to see the fields
(res <- package_show(pkg$id))

## update just chosen things
# Make some changes
x <- list(maintainer_email = "heythere2@things.com")

# Then update the packge
package_update(x, pkg$id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/related_update.R
\name{related_update}
\alias{related_update}
\title{Update a related item}
\usage{
related_update(
  id,
  title,
  type,
  description = NULL,
  related_id = NULL,
  related_url = NULL,
  image_url = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) id of related item to update. This should be an
alphanumeric string. Required.}

\item{title}{(character) Title of the related item. Required.}

\item{type}{(character) The type of the related item. One of API,
application, idea, news article, paper, post or visualization. Required.}

\item{description}{(character) description (optional). Optional}

\item{related_id}{(character) An id to assign to the related item. If blank,
an ID will be assigned for you. Optional}

\item{related_url}{(character) A url to associated with the related item.
Optional}

\item{image_url}{(character) A url to associated image. Optional}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Update a related item
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create a package and related item
res <- package_create("hello-saturn2") \%>\%
   related_create(title = "my resource",
                  type = "visualization")

# update the related item
related_update(res, title = "her resource", type = "idea")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag_show.R
\name{tag_show}
\alias{tag_show}
\title{Show a tag.}
\usage{
tag_show(
  id,
  include_datasets = FALSE,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) The name or id of the tag}

\item{include_datasets}{include a list of up to 1000 of the tag's datasets.
Limit 1000 datasets, use \code{\link[=package_search]{package_search()}} for more.
(optional, default: \code{FALSE})}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Show a tag.
}
\examples{
\dontrun{
# get tags with tag_list()
tags <- tag_list()
tags[[30]]$id

# show a tag
(x <- tag_show(tags[[30]]$id))

# give back different data formats
tag_show(tags[[30]]$id, as = 'json')
tag_show(tags[[30]]$id, as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_revision_list.R
\name{package_revision_list}
\alias{package_revision_list}
\title{Return a dataset (package's) revisions as a list of dictionaries.}
\usage{
package_revision_list(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Package identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a dataset (package's) revisions as a list of dictionaries.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("dolphins"))

# list package revisions
package_revision_list(res$id)

# Make change to the package
x <- list(title = "dolphins and things")
package_patch(x, id = res$id)

# list package revisions
package_revision_list(res$id)

# Output different formats
package_revision_list(res$id, as = "table")
package_revision_list(res$id, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_list.R
\name{group_list}
\alias{group_list}
\title{List groups.}
\usage{
group_list(
  offset = 0,
  limit = 31,
  sort = NULL,
  groups = NULL,
  all_fields = FALSE,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{sort}{Field to sort on. You can specify ascending (e.g., score desc) or
descending (e.g., score asc), sort by two fields (e.g., score desc, price asc),
or sort by a function (e.g., sum(x_f, y_f) desc, which sorts by the sum of
x_f and y_f in a descending order).}

\item{groups}{(character) A list of names of the groups to return, if given
only groups whose names are in this list will be returned}

\item{all_fields}{(logical) Return full group dictionaries instead of just
names. Default: \code{FALSE}}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
List groups.
}
\examples{
\dontrun{
group_list(limit = 3)
group_list(limit = 3, as = 'json')
group_list(limit = 3, as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_patch.R
\name{resource_patch}
\alias{resource_patch}
\title{Update a resource's metadata}
\usage{
resource_patch(
  x,
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{x}{(list) A list with key-value pairs}

\item{id}{(character) Resource ID to update (required)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Update a resource's metadata
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org", key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("twist", author="Alexandria"))

# then create a resource
file <- system.file("examples", "actinidiaceae.csv", package = "ckanr")
(xx <- resource_create(package_id = res$id, description = "my resource"))

# Get a resource
res <- resource_show(xx$id)
res$description

# Make some changes
x <- list(description = "My newer description")
z <- resource_patch(x, id = res)
z$description

# Add an extra key:value pair
extra <- list("extra_key" = "my special value")
zz <- resource_patch(extra, id = res)
zz$extra_key
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ckan_package.R
\name{as.ckan_package}
\alias{as.ckan_package}
\alias{is.ckan_package}
\title{ckan_package class helpers}
\usage{
as.ckan_package(x, ...)

is.ckan_package(x)
}
\arguments{
\item{x}{Variety of things, character, list, or ckan_package class object}

\item{...}{Further args passed on to \code{\link[=package_show]{package_show()}} if
character given. In particular, if GET is not supported you can
try the \code{http_method} parameter to set a different HTTP verb}
}
\description{
ckan_package class helpers
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
  key = getOption("ckan_demo_key"))

(pkgs <- package_search())
pkgs$results
pkgs$results[[3]]

# create item class from only an item ID
as.ckan_package(pkgs$results[[3]]$id)

# gives back itself
(x <- as.ckan_package(pkgs$results[[3]]$id))
as.ckan_package(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds_create_dataset.R
\name{ds_create_dataset}
\alias{ds_create_dataset}
\title{Datastore - create a new resource on an existing dataset}
\usage{
ds_create_dataset(
  package_id,
  name,
  path,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{package_id}{(character) Existing package ID (required)}

\item{name}{(character) Name of the new resource (required)}

\item{path}{(character) Path of the file to add (required)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Datastore - create a new resource on an existing dataset
}
\details{
This function is deprecated - will be defunct in the next version
of this package
}
\examples{
\dontrun{
path <- system.file("examples", "actinidiaceae.csv", package = "ckanr")
ckanr_setup(url = "https://demo.ckan.org/", key = "my-demo-ckan-org-api-key")
ds_create_dataset(package_id='testingagain', name="mydata", path = path)

# Testing: see ?ckanr_setup to set test settings
ckanr_setup(test_url = "http://my-ckan.org/",
            test_key = "my-ckan-api-key",
            test_did="an-existing-package-id",
            test_rid="an-existing-resource-id")
ds_create_dataset(package_id=get_test_pid(), name="mydata",
                  path=system.file("examples",
                                   "actinidiaceae.csv",
                                   package = "ckanr"),
                  key = get_test_key(),
                  url = get_test_url())
}
}
\references{
http://docs.ckan.org/en/latest/api/index.html#ckan.logic.action.create.resource_create
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckanr-package.R
\name{ckanr-deprecated}
\alias{ckanr-deprecated}
\title{Deprecated functions in \pkg{ckanr}}
\description{
These functions still work but will be removed (defunct) in the next version.
}
\details{
\itemize{
\item \code{\link[=ds_create_dataset]{ds_create_dataset()}}: The functionality of this function is already in
another function in this package. See function \code{\link[=resource_create]{resource_create()}}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_show.R
\name{package_show}
\alias{package_show}
\title{Show a package.}
\usage{
package_show(
  id,
  use_default_schema = FALSE,
  http_method = "GET",
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character/ckan_package) Package identifier, or a \code{ckan_package}
object}

\item{use_default_schema}{(logical) Use default package schema instead of a
custom schema defined with an IDatasetForm plugin. Default: \code{FALSE}}

\item{http_method}{(character) which HTTP method (verb) to use; one of
"GET" or "POST". Default: "GET"}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Show a package.
}
\details{
By default the help and success slots are dropped, and only the
result slot is returned. You can request raw json with \code{as = 'json'}
then parse yourself to get the help slot.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/",
  key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("purposeful55"))

# show package
## From the output of package_create
package_show(res)
## Or, from the ID
package_show(res$id)

# get data back in different formats
package_show(res$id, as = 'json')
package_show(res$id, as = 'table')

# use default schema or not
package_show(res$id, TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_patch.R
\name{group_patch}
\alias{group_patch}
\title{Update a group's metadata}
\usage{
group_patch(
  x,
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{x}{(list) A list with key-value pairs}

\item{id}{(character) Resource ID to update (required)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Update a group's metadata
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org", key = getOption("ckan_demo_key"))

# Create a package
(res <- group_create("hello-my-world2"))

# Get a resource
grp <- group_show(res$id)
grp$title
grp$author_email

# Make some changes
x <- list(title = "!hello world!", maintainer_email = "hello@world.com")
group_patch(x, id = grp)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ping.R
\name{ping}
\alias{ping}
\title{Ping a CKAN server to test that it's up or down.}
\usage{
ping(url = get_default_url(), key = get_default_key(), as = "logical", ...)
}
\arguments{
\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Ping a CKAN server to test that it's up or down.
}
\examples{
\dontrun{
ping()
ping(as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds_create.R
\name{ds_create}
\alias{ds_create}
\title{Add a new table to a datastore}
\usage{
ds_create(
  resource_id = NULL,
  resource = NULL,
  force = FALSE,
  aliases = NULL,
  fields = NULL,
  records = NULL,
  primary_key = NULL,
  indexes = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{resource_id}{(string) Resource id that the data is going to be stored
against.}

\item{resource}{(dictionary) Resource dictionary that is passed to
\code{\link[=resource_create]{resource_create()}}. Use instead of \code{resource_id} (optional)}

\item{force}{(logical) Set to \code{TRUE} to edit a read-only resource.
Default: \code{FALSE}}

\item{aliases}{(character) Names for read only aliases of the resource.
(optional)}

\item{fields}{(list) Fields/columns and their extra metadata. (optional)}

\item{records}{(list) The data, eg: \verb{[\{"dob": "2005", "some_stuff": ["a", "b"]\}]} (optional)}

\item{primary_key}{(character) Fields that represent a unique key (optional)}

\item{indexes}{(character) Indexes on table (optional)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
BEWARE: This function still doesn't quite work yet.
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
  key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("foobarrrrr", author="Jane Doe"))

# then create a resource
file <- system.file("examples", "actinidiaceae.csv", package = "ckanr")
(xx <- resource_create(package_id = res$id,
                       description = "my resource",
                       name = "bears",
                       upload = file,
                       rcurl = "http://google.com"
))
ds_create(resource_id = xx$id, records = iris, force = TRUE)
resource_show(xx$id)
}
}
\references{
http://bit.ly/ds_create
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckan_info.R
\name{ckan_info}
\alias{ckan_info}
\alias{ckan_version}
\title{Get information on a CKAN server}
\usage{
ckan_info(url = get_default_url(), ...)

ckan_version(url = get_default_url(), ...)
}
\arguments{
\item{url}{Base url to use. Default: \url{https://data.ontario.ca}. See
also \code{\link[=ckanr_setup]{ckanr_setup()}} and \code{\link[=get_default_url]{get_default_url()}}. (required)}

\item{...}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET} (optional)}
}
\value{
for \code{ckan_info} a list with many slots with various info.
for \code{ckan_version}, list of length two, with actual version as character,
and another with version converted to numeric (any dots or letters removed)
}
\description{
Get information on a CKAN server
}
\examples{
\dontrun{
ckan_info()
ckan_info(servers()[5])

ckan_version(servers()[5])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr.R
\name{src_ckan}
\alias{src_ckan}
\alias{dplyr-interface}
\title{Connect to CKAN with dplyr}
\usage{
src_ckan(url)
}
\arguments{
\item{url, }{the url of the CKAN instance}
}
\description{
Use \code{src_ckan} to connect to an existing CKAN instance and \code{tbl} to
connect to tables within that CKAN based on the DataStore Data API.
}
\examples{
\dontrun{
library("dplyr")

# To connect to a CKAN instance first create a src:
my_ckan <- src_ckan("http://demo.ckan.org")

# List all tables in the CKAN instance
db_list_tables(my_ckan$con)

# Then reference a tbl within that src
my_tbl <- tbl(src = my_ckan, name = "44d7de5f-7029-4f3a-a812-d7a70895da7d")

# You can use the dplyr verbs with my_tbl. For example:
dplyr::filter(my_tbl, GABARITO == "C")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/revision_list.R
\name{revision_list}
\alias{revision_list}
\title{Return a list of the IDs of the site's revisions.}
\usage{
revision_list(
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a list of the IDs of the site's revisions.
}
\examples{
\dontrun{
revision_list()
revision_list(as = "table")
revision_list(as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_list.R
\name{user_list}
\alias{user_list}
\title{Return a list of the site's user accounts.}
\usage{
user_list(
  q = NULL,
  order_by = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{q}{(character) Restrict the users returned to those whose names
contain a string}

\item{order_by}{(character) Which field to sort the list by
(optional, default: 'name')}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a list of the site's user accounts.
}
\examples{
\dontrun{
# all users
user_list()

# search for a user
user_list(q = "j")

# different data formats
user_list(q = "j", as = "table")
user_list(q = "j", as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_list.R
\name{package_list}
\alias{package_list}
\title{List datasets.}
\usage{
package_list(
  offset = 0,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
List datasets.
}
\examples{
\dontrun{
package_list()
package_list(as = 'json')
package_list(as = 'table')

package_list(url = 'https://data.nhm.ac.uk')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_create.R
\name{user_create}
\alias{user_create}
\title{Create a user.}
\usage{
user_create(
  name,
  email,
  password,
  id = NULL,
  fullname = NULL,
  about = NULL,
  openid = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{name}{(character) the name of the new user, a string between 2 and 100
characters in length, containing only lowercase alphanumeric
characters, - and _ (required)}

\item{email}{(character) the email address for the new user (required)}

\item{password}{(character) the password of the new user, a string of at
least 4 characters (required)}

\item{id}{(character) the id of the new user (optional)}

\item{fullname}{(character) user full name}

\item{about}{(character) a description of the new user (optional)}

\item{openid}{(character) an openid (optional)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Create a user.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://data-demo.dpaw.wa.gov.au",
  key = "824e7c50-9577-4bfa-bf32-246ebed1a8a2")

# create a user
user_create(name = 'stacy', email = "stacy@aaaaa.com",
password = "helloworld")
}
}
\references{
http://docs.ckan.org/en/latest/api/index.html#ckan.logic.action.create.user_create
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_activity_list.R
\name{package_activity_list}
\alias{package_activity_list}
\title{Return a list of the package's activity}
\usage{
package_activity_list(
  id,
  offset = 0,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Package identifier.}

\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a list of the package's activity
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("owls64"))

# list package activity
package_activity_list(res$id)

# make a change
x <- list(maintainer = "Jane Forest")
package_update(x, res)

# list activity again
package_activity_list(res)

# output different data formats
package_activity_list(res$id, as = "table")
package_activity_list(res$id, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ckan_organization.R
\name{as.ckan_organization}
\alias{as.ckan_organization}
\alias{is.ckan_organization}
\title{ckan_organization class helpers}
\usage{
as.ckan_organization(x, ...)

is.ckan_organization(x)
}
\arguments{
\item{x}{Variety of things, character, list, or ckan_organization class
object}

\item{...}{Further args passed on to \code{\link[=organization_show]{organization_show()}} if character
given}
}
\description{
ckan_organization class helpers
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
key = getOption("ckan_demo_key"))

(orgs <- organization_list(limit = 3))
orgs[[3]]

# create item class from only an item ID
as.ckan_organization(orgs[[3]]$id)

# gives back itself
(x <- as.ckan_organization(orgs[[3]]$id))
as.ckan_organization(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organization_list.R
\name{organization_list}
\alias{organization_list}
\title{List organization}
\usage{
organization_list(
  order_by = c("name", "package"),
  decreasing = FALSE,
  organizations = NULL,
  all_fields = TRUE,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{order_by}{(character, only the first element is used).
The field to sort the list by, must be \code{name} or \code{packages}.}

\item{decreasing}{(logical). Is the sort-order is decreasing or not.}

\item{organizations}{(character or \code{NULL}). A list of names of the
organizations to return. \code{NULL} returns all organizations.}

\item{all_fields}{(logical). Return the name or all fields of the object.}

\item{limit}{(numeric) The maximum number of organizations to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
List organization
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/")

# list organizations
res <- organization_list()
res[1:2]

# Different data formats
organization_list(as = 'json')
organization_list(as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/related_list.R
\name{related_list}
\alias{related_list}
\title{List related items}
\usage{
related_list(
  offset = 0,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
List related items
}
\examples{
\dontrun{
related_list()
related_list(as = 'json')
related_list(as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_follower_count.R
\name{user_follower_count}
\alias{user_follower_count}
\title{Return a a user's follower count}
\usage{
user_follower_count(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) User identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a a user's follower count
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# list package activity
user_follower_count('sckottie')

# input a ckan_user object
(x <- user_show('sckottie'))
user_follower_count(x)

# output different data formats
user_follower_count(x, as = "table")
user_follower_count(x, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ckan_user.R
\name{as.ckan_user}
\alias{as.ckan_user}
\alias{is.ckan_user}
\title{ckan_user class helpers}
\usage{
as.ckan_user(x, ...)

is.ckan_user(x)
}
\arguments{
\item{x}{Variety of things, character, list, or ckan_user class object}

\item{...}{Further args passed on to \code{\link[=user_show]{user_show()}} if character given}
}
\description{
ckan_user class helpers
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
  key = getOption("ckan_demo_key"))

(usrs <- user_list())
usrs[1:3]
usrs[[3]]

# create item class from only an item ID
as.ckan_user(usrs[[3]]$id)

# gives back itself
(x <- as.ckan_user(usrs[[3]]$id))
as.ckan_user(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_follower_list.R
\name{user_follower_list}
\alias{user_follower_list}
\title{Return a a user's follower count}
\usage{
user_follower_list(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) User identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a a user's follower count
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# list package activity
user_follower_list('sckottie')

# input a ckan_user object
(x <- user_show('sckottie'))
user_follower_list(x)

# output different data formats
user_follower_list(x, as = "table")
user_follower_list(x, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag_search.R
\name{tag_search}
\alias{tag_search}
\title{Searcn tags.}
\usage{
tag_search(
  query = NULL,
  vocabulary_id = NULL,
  offset = 0,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{query}{(character) A tag name query to search for, if given only tags
whose names contain this string will be returned; one or more search
strings}

\item{vocabulary_id}{(character) The id or name of a vocabulary,
if give only tags that belong to this vocabulary will be returned}

\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Searcn tags.
}
\examples{
\dontrun{
tag_search(query = 'ta')
tag_search(query = c('ta', 'al'))

# different formats back
tag_search(query = 'ta', as = 'json')
tag_search(query = 'ta', as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organization_delete.R
\name{organization_delete}
\alias{organization_delete}
\title{Delete an organization}
\usage{
organization_delete(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) name or id of the organization}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\value{
an empty list on success
}
\description{
Delete an organization
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org", key=getOption("ckan_demo_key"))

# create an organization
(res <- organization_create("foobar", title = "Foo bars",
  description = "love foo bars"))

# delete the organization just created
res$id
organization_delete(id = res$id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_show.R
\name{group_show}
\alias{group_show}
\title{Show a package}
\usage{
group_show(
  id,
  include_datasets = TRUE,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Package identifier.}

\item{include_datasets}{(logical) Include a list of the group's datasets.
Default: \code{TRUE}}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Show a package
}
\details{
By default the help and success slots are dropped, and only the
result slot is returned. You can request raw json with \code{as = 'json'}
then parse yourself to get the help slot.
}
\examples{
\dontrun{
res <- group_list()

# via a group name/id
group_show(res[[1]]$name)

# or via an object of class ckan_group
group_show(res[[1]])

# return different data formats
group_show(res[[1]]$name, as = 'json')
group_show(res[[1]]$name, as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_create.R
\name{resource_create}
\alias{resource_create}
\title{Create a resource}
\usage{
resource_create(
  package_id = NULL,
  rcurl = NULL,
  revision_id = NULL,
  description = NULL,
  format = NULL,
  hash = NULL,
  name = NULL,
  resource_type = NULL,
  mimetype = NULL,
  mimetype_inner = NULL,
  webstore_url = NULL,
  cache_url = NULL,
  size = NULL,
  created = NULL,
  last_modified = NULL,
  cache_last_updated = NULL,
  webstore_last_updated = NULL,
  upload = NULL,
  extras = NULL,
  http_method = "GET",
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{package_id}{(character) id of package that the resource should be
added to. This should be an alphanumeric string. Required.}

\item{rcurl}{(character) url of resource. Required.}

\item{revision_id}{(character) revision id (optional)}

\item{description}{(character) description (optional). Required.}

\item{format}{(character) format (optional)}

\item{hash}{(character) hash (optional)}

\item{name}{(character) name (optional). Required.}

\item{resource_type}{(character) resource type (optional)}

\item{mimetype}{(character) mime type (optional)}

\item{mimetype_inner}{(character) mime type inner (optional)}

\item{webstore_url}{(character) webstore url (optional)}

\item{cache_url}{(character) cache url(optional)}

\item{size}{(integer) size (optional)}

\item{created}{(character) iso date string (optional)}

\item{last_modified}{(character) iso date string (optional)}

\item{cache_last_updated}{(character) iso date string (optional)}

\item{webstore_last_updated}{(character) iso date string (optional)}

\item{upload}{(character) A path to a local file (optional)}

\item{extras}{(list) - the resources' extra metadata fields (optional)}

\item{http_method}{(character) which HTTP method (verb) to use; one of
"GET" or "POST". Default: "GET"}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Create a resource
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/",
 key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("foobarrrr", author="Jane Doe"))

# then create a resource
file <- system.file("examples", "actinidiaceae.csv", package = "ckanr")
(xx <- resource_create(package_id = res$id,
                       description = "my resource",
                       name = "bears",
                       upload = file,
                       extras = list(species = "grizzly"),
                       rcurl = "http://google.com"
))

package_create("foobbbbbarrrr") \%>\%
   resource_create(description = "my resource",
                   name = "bearsareus",
                   upload = file,
                   extras = list(my_extra = "some value"),
                   rcurl = "http://google.com")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organization_create.R
\name{organization_create}
\alias{organization_create}
\title{Create an organization}
\usage{
organization_create(
  name = NULL,
  id = NULL,
  title = NULL,
  description = NULL,
  image_url = NULL,
  state = "active",
  approval_status = NULL,
  extras = NULL,
  packages = NULL,
  users = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{name}{(character) the name of the organization, a string between 2 and
100 characters long, containing only lowercase alphanumeric characters,
\itemize{
\item and _
}}

\item{id}{the id of the organization (optional)}

\item{title}{(character) the title of the organization (optional)}

\item{description}{(character) the description of the organization (optional)}

\item{image_url}{(character) the URL to an image to be displayed on the
organization's page (optional)}

\item{state}{(character) the current state of the organization, e.g.
'active' or 'deleted', only active organization show up in search results
and other lists of organization, this parameter will be ignored if you are
not authorized to change the state of the organization (optional).
Default: 'active'}

\item{approval_status}{(character) Approval status}

\item{extras}{The organization's extras (optional), extras are arbitrary
(key: value) metadata items that can be added to organizations, each extra
dictionary should have keys 'key' (a string), 'value' (a string)
\code{package_relationship_create} for the format of relationship dictionaries
(optional)}

\item{packages}{(list of dictionaries) the datasets (packages) that belong
to the organization, a list of dictionaries each with keys 'name' (string,
the id or name of the dataset) and optionally 'title' (string, the title
of the dataset)}

\item{users}{(character) the users that belong to the organization, a list
of dictionaries each with key 'name' (string, the id or name of the user)
and optionally 'capacity' (string, the capacity in which the user is a
member of the organization)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Create an organization
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create an organization
(res <- organization_create("foobar", title = "Foo bars",
  description = "love foo bars"))
res$name
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckanr_settings.R
\name{ckanr_setup}
\alias{ckanr_setup}
\title{Configure default CKAN settings}
\usage{
ckanr_setup(
  url = "https://data.ontario.ca/",
  key = NULL,
  test_url = NULL,
  test_key = NULL,
  test_did = NULL,
  test_rid = NULL,
  test_gid = NULL,
  test_oid = NULL,
  test_behaviour = NULL,
  proxy = NULL
)
}
\arguments{
\item{url}{A CKAN URL (optional), default: https://data.ontario.ca}

\item{key}{A CKAN API key (optional, character)}

\item{test_url}{(optional, character) A valid CKAN URL for testing purposes}

\item{test_key}{(optional, character) A valid CKAN API key privileged to
create datasets at \code{test_url}}

\item{test_did}{(optional, character) A valid CKAN dataset ID, existing at
\code{test_url}}

\item{test_rid}{(optional, character) A valid CKAN resource ID, attached to
\code{did}}

\item{test_gid}{(optional, character) A valid CKAN group name at \code{test_url}}

\item{test_oid}{(optional, character) A valid CKAN organization name at
\code{test_url}}

\item{test_behaviour}{(optional, character) Whether to fail ("FAIL") or skip
("SKIP") writing tests in case of problems with the configured test CKAN.}

\item{proxy}{an object of class \code{request} from a call to
\code{\link[crul:proxies]{crul::proxy()}}}
}
\description{
Configure default CKAN settings
}
\details{
\code{\link[=ckanr_setup]{ckanr_setup()}} sets CKAN connection details. ckanr's functions
default to use the default URL and API key unless specified explicitly.

ckanr's automated tests require a valid CKAN URL, a privileged API key
for that URL, plus the IDs of an existing dataset and an existing resource,
repectively.

The writing tests (create, update, delete) can fail for two reasons:
failures in ckanr's code which the tests aim to detect,
or failures in the configured CKAN, which are not necessarily a problem
with ckanr's code but prevent the tests to prove otherwise.

Setting \code{test_behaviour} to \code{"SKIP"} will allow writing tests to skip
if the configured test CKAN fails. This is desirable to e.g. test the other
functions even if the tester has no write access to a CKAN instance.

Setting \code{test_behaviour} to \code{"FAIL"} will let the tester find any
problems with both the configured test CKAN and the writing functions.
}
\examples{
# CKAN users without admin/editor privileges could run:
ckanr_setup(url = "https://data.ontario.ca/")

# Privileged CKAN editor/admin users can run:
ckanr_setup(url = "https://data.ontario.ca/", key = "some-CKAN-API-key")

# ckanR developers/testers can run:
ckanr_setup(url = "https://data.ontario.ca/", key = "some-CKAN-API-key",
           test_url = "http://test-ckan.gov/",test_key = "test-ckan-API-key",
           test_did = "test-ckan-dataset-id",test_rid = "test-ckan-resource-id",
           test_gid = "test-group-name", test_oid = "test-organzation-name",
           test_behaviour = "FAIL")

# Not specifying the default CKAN URL will reset the CKAN URL to its default
# "https://data.ontario.ca/":
ckanr_setup()

# set a proxy
ckanr_setup(proxy = crul::proxy("64.251.21.73:8080"))
ckanr_settings()
## run without setting proxy to reset to no proxy
ckanr_setup()
ckanr_settings()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_search.R
\name{package_search}
\alias{package_search}
\title{Search for packages.}
\usage{
package_search(
  q = "*:*",
  fq = NULL,
  sort = NULL,
  rows = NULL,
  start = NULL,
  facet = FALSE,
  facet.limit = NULL,
  facet.field = NULL,
  facet.mincount = NULL,
  include_drafts = FALSE,
  include_private = FALSE,
  use_default_schema = FALSE,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{q}{Query terms, defaults to '\emph{:}', or everything.}

\item{fq}{Filter query, this does not affect the search, only what gets
returned}

\item{sort}{Field to sort on. You can specify ascending (e.g., score desc)
or descending (e.g., score asc), sort by two fields (e.g., score desc,
price asc), or sort by a function (e.g., sum(x_f, y_f) desc, which sorts
by the sum of x_f and y_f in a descending order).}

\item{rows}{Number of records to return. Defaults to 10.}

\item{start}{Record to start at, default to beginning.}

\item{facet}{(logical) Whether to return facet results or not.
Default: \code{FALSE}}

\item{facet.limit}{(numeric) This param indicates the maximum number of
constraint counts that should be returned for the facet fields.
A negative value means unlimited. Default: 100.
Can be specified on a per field basis.}

\item{facet.field}{(charcter) This param allows you to specify a field which
should be treated as a facet. It will iterate over each Term in the field
and generate a facet count using that Term as the constraint. This parameter
can be specified multiple times to indicate multiple facet fields. None of
the other params in this section will have any effect without specifying at
least one field name using this param.}

\item{facet.mincount}{(integer) the minimum counts for facet fields should
be included in the results}

\item{include_drafts}{(logical) if \code{TRUE} draft datasets will be
included. A user will only be returned their own draft datasets, and a
sysadmin will be returned all draft datasets. default: \code{FALSE}.
first CKAN version: 2.6.1; dropped from request if CKAN version is older
or if CKAN version isn't available via \code{\link[=ckan_version]{ckan_version()}}}

\item{include_private}{(logical) if \code{TRUE} private datasets will be
included. Only private datasets from the user’s organizations will be
returned and sysadmins will be returned all private datasets.
default: \code{FALSE}
first CKAN version: 2.6.1; dropped from request if CKAN version is older
or if CKAN version isn't available via \code{\link[=ckan_version]{ckan_version()}}}

\item{use_default_schema}{(logical) use default package schema instead of a
custom schema defined with an IDatasetForm plugin. default: \code{FALSE}
first CKAN version: 2.3.5; dropped from request if CKAN version is older
or if CKAN version isn't available via \code{\link[=ckan_version]{ckan_version()}}}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Search for packages.
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org", key=getOption("ckan_demo_key"))

package_search(q = '*:*')
package_search(q = '*:*', rows = 2, as = 'json')
package_search(q = '*:*', rows = 2, as = 'table')

package_search(q = '*:*', sort = 'score asc')
package_search(q = '*:*', fq = 'num_tags:[3 TO *]')$count
package_search(q = '*:*', fq = 'num_tags:[2 TO *]')$count
package_search(q = '*:*', fq = 'num_tags:[1 TO *]')$count
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_update.R
\name{resource_update}
\alias{resource_update}
\title{Update a resource}
\usage{
resource_update(
  id,
  path = NULL,
  extras = list(),
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Resource ID to update (required)}

\item{path}{(character) Local path of the file to upload (optional)}

\item{extras}{(list) - the resources' extra metadata fields (optional)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\value{
The HTTP response from CKAN, formatted as list (default), table,
or JSON.
}
\description{
This function can be used to update a resource's file attachment
and "extra" metadata fields. Any update will also set the metadata key
"last_updated". Other metadata, such as name or description, are not updated.

The new file must exist on a local path. R objects have to be written to a
file, e.g. using \code{tempfile()} - see example.

For convenience, CKAN base url and API key default to the global options,
which are set by \code{ckanr_setup}.
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# Get file
path <- system.file("examples", "actinidiaceae.csv", package = "ckanr")

# Create package, then a resource within that package
(res <- package_create("newpackage10"))
(xx <- resource_create(package_id = res$id,
                       description = "my resource",
                       name = "bears",
                       upload = path,
                       rcurl = "http://google.com"
))

# Modify dataset, here lowercase strings in one column
dat <- read.csv(path, stringsAsFactors = FALSE)
dat$Family <- tolower(dat$Family)
newpath <- tempfile(fileext = ".csv")
write.csv(dat, file = newpath, row.names = FALSE)

# Upload modified dataset
## Directly from output of resource_create
resource_update(xx, path=newpath)

## or from the resource id
resource_update(xx$id, path=newpath)

## optionally include extra tags
resource_update(xx$id, path=newpath,
                extras = list(some="metadata"))
                
# Update a resource's extra tags
## add extra tags without uploading a new file
resource_update(id,
                extras = list(some="metadata"))

## or remove all extra tags
resource_update(id, extras = list())                 

#######
# Using default settings
ckanr_setup(url = "http://demo.ckan.org/", key = "my-demo-ckan-org-api-key")
path <- system.file("examples", "actinidiaceae.csv", package = "ckanr")
resource_update(id="an-existing-resource-id", path = path)

# Using an R object written to a tempfile, and implicit CKAN URL and API key
write.csv(data <- installed.packages(), path <- tempfile(fileext = ".csv"))
ckanr_setup(url = "http://demo.ckan.org/", key = "my-demo-ckan-org-api-key")
resource_update(id="an-existing-resource-id", path = path)

# Testing: see ?ckanr_setup to set default test CKAN url, key, package id
ckanr_setup(test_url = "http://my-ckan.org/",
            test_key = "my-ckan-api-key",
            test_did = "an-existing-package-id",
            test_rid = "an-existing-resource-id")
resource_update(id = get_test_rid(),
                path = system.file("examples",
                                   "actinidiaceae.csv",
                                   package = "ckanr"),
                key = get_test_key(),
                url = get_test_url())

# other file formats
## html
path <- system.file("examples", "mapbox.html", package = "ckanr")

# Create package, then a resource within that package
(res <- package_create("mappkg"))
(xx <- resource_create(package_id = res$id,
                       description = "a map, yay",
                       name = "mapyay",
                       upload = path,
                       rcurl = "http://google.com"
))
browseURL(xx$url)

# Modify dataset, here lowercase strings in one column
dat <- readLines(path)
dat <- sub("-111.06", "-115.06", dat)
newpath <- tempfile(fileext = ".html")
cat(dat, file = newpath, sep = "\n")

# Upload modified dataset
## Directly from output of resource_create
(xxx <- resource_update(xx, path=newpath))
browseURL(xxx$url)
}
}
\references{
http://docs.ckan.org/en/latest/api/index.html#ckan.logic.action.create.resource_create
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organization_show.R
\name{organization_show}
\alias{organization_show}
\title{Show an organization}
\usage{
organization_show(
  id,
  include_datasets = FALSE,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Organization id or name.}

\item{include_datasets}{(logical). Whether to include a list of the
organization datasets}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Show an organization
}
\details{
By default the help and success slots are dropped, and only the
result slot is returned. You can request raw json with \code{as = 'json'}
then parse yourself to get the help slot.
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

res <- organization_create("stuffthings2")
organization_show(res$id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckanr-package.R
\docType{package}
\name{ckanr-package}
\alias{ckanr-package}
\alias{ckanr}
\title{R client for the CKAN API}
\description{
ckanr is a full client for the CKAN API, wrapping all
APIs, including for reading and writing data. Please get in touch
(\url{https://github.com/ropensci/ckanr/issues} or \url{https://discuss.ropensci.org/})
if you have problems, or have use cases that we don't cover yet.
}
\section{CKAN API}{


Document for the CKAN API is at \url{https://docs.ckan.org/en/latest/api/index.html}.
We'll always be following the lastest version of the API.
}

\section{ckanr package API}{


The functions can be grouped into those for setup, packages,
resources, tags, organizations, groups, and users.
\itemize{
\item Setup - The main one is \code{\link[=ckanr_setup]{ckanr_setup()}} - and many related
functions, e.g., \code{\link[=get_default_key]{get_default_key()}}
\item Packages - Create a package with \code{\link[=package_create]{package_create()}}, and see
other functions starting with \verb{package_*}
\item Resources - Create a package with \code{\link[=resource_create]{resource_create()}}, and see
other functions starting with \verb{resource_*}
\item Tags - List tags with \code{\link[=tag_list]{tag_list()}}, and see
other functions starting with \verb{tag_*}
\item Organizations - List organizations with \code{\link[=organization_list]{organization_list()}},
show a specific organization with \code{\link[=organization_show]{organization_show()}}, and
create with \code{\link[=organization_create]{organization_create()}}
\item Groups - List groups with \code{\link[=group_list]{group_list()}}, and see
other functions starting with \verb{group_*}
\item Users - List users with \code{\link[=user_list]{user_list()}}, and see
other functions starting with \verb{user_*}
\item Related items - See functions starting with \verb{related_*}
}
}

\section{Datastore}{


We are also working on supporting the Datastore extension
(\url{https://docs.ckan.org/en/latest/maintaining/datastore.html}).
We currently have these functions:
\itemize{
\item \code{\link[=ds_create]{ds_create()}}
\item \code{\link[=ds_create_dataset]{ds_create_dataset()}}
\item \code{\link[=ds_search]{ds_search()}}
\item \code{\link[=ds_search_sql]{ds_search_sql()}}
}
}

\section{Fetch}{


Data can come back in a huge variety of formats. We've attempted a function to
help you fetch not just metadata but the actual data for a link to a file on
a CKAN instance. Though if you know what you're doing, you can easily use
whatever is your preferred tool for the job (e.g., maybe you like
\code{\link[=read.csv]{read.csv()}} for reading csv files).
}

\section{CKAN Instances}{


We have a helper function (\code{\link[=servers]{servers()}}) that spits out the current
CKAN instances we know about, with URLs to their base URLs that should work
using this package. That is, not necessarily landing pages of each instance,
although, the URL may be the landing page and the base API URL.
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}

Florian Mayer \email{florian.wendelin.mayer@gmail.com}

Wush Wu

Imanuel Costigan \email{i.costigan@me.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag_list.R
\name{tag_list}
\alias{tag_list}
\title{List tags.}
\usage{
tag_list(
  query = NULL,
  vocabulary_id = NULL,
  all_fields = FALSE,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{query}{(character) A tag name query to search for, if given only tags
whose names contain this string will be returned}

\item{vocabulary_id}{(character) The id or name of a vocabulary, if given,
only tags that belong to this vocabulary will be returned}

\item{all_fields}{(logical) Return full tag dictionaries instead of
just names. Default: \code{FALSE}}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
List tags.
}
\examples{
\dontrun{
# list all tags
tag_list()

# search for a specific tag
tag_list(query = 'aviation')

# all fields
tag_list(all_fields = TRUE)

# give back different data formats
tag_list('aviation', as = 'json')
tag_list('aviation', as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckan_classes.R
\name{ckan_classes}
\alias{ckan_classes}
\title{ckanr S3 classes}
\description{
ckanr S3 classes
}
\section{The classes}{

\itemize{
\item ckan_package - CKAN package
\item ckan_resource - CKAN resource
\item ckan_related - CKAN related item
}
}

\section{Coercion}{

The functions \verb{as.ckan_*()} for each CKAN object type coerce something
to a S3 class of that type. For example, you can coerce a package ID as a
character string into an \code{ckan_package} object by calling
\verb{as.ckan_package(<id>}.
}

\section{Testing for classes}{

To test whether an object is of a particular \verb{ckan_*} class, there is a
\verb{is._ckan_*()} function for all of the classes listed above. You can use
one of those functions to get a logical back, \code{TRUE} or \code{FALSE}.
}

\section{Manipulation}{

These are simple S3 classes, basically an R list with an attached class
so we can know what to do with the object and have flexible inputs and
outputs from functions. You can edit one of these classes yourself
by simply changing values in the list.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ckan_tag.R
\name{as.ckan_tag}
\alias{as.ckan_tag}
\alias{is.ckan_tag}
\title{ckan_tag class helpers}
\usage{
as.ckan_tag(x, ...)

is.ckan_tag(x)
}
\arguments{
\item{x}{Variety of things, character, list, or ckan_tag class object}

\item{...}{Further args passed on to \code{\link[=tag_show]{tag_show()}} if character given}
}
\description{
ckan_tag class helpers
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
key = getOption("ckan_demo_key"))

(tags <- tag_search(query = 'ta'))
tags[[3]]

# create item class from only an item ID
as.ckan_tag(tags[[3]]$id)

# gives back itself
(x <- as.ckan_tag(tags[[3]]$id))
as.ckan_tag(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_search.R
\name{resource_search}
\alias{resource_search}
\title{Search for resources.}
\usage{
resource_search(
  q,
  sort = NULL,
  offset = NULL,
  limit = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{q}{Query terms. It is a string of the form \code{field:term} or a
vector/list of strings, each of the same form.  Within each string, \code{field}
is a field or extra field on the Resource domain object. If \code{field} is
hash, then an attempt is made to match the \code{term} as a \emph{prefix} of the
Resource.hash field. If \code{field} is an extra field, then an attempt is
made to match against the extra fields stored against the Resource.}

\item{sort}{Field to sort on. You can specify ascending (e.g., score desc)
or descending (e.g., score asc), sort by two fields (e.g., score desc,
price asc), or sort by a function (e.g., sum(x_f, y_f) desc, which sorts
by the sum of x_f and y_f in a descending order).}

\item{offset}{Record to start at, default to beginning.}

\item{limit}{Number of records to return.}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Search for resources.
}
\examples{
\dontrun{
resource_search(q = 'name:data')
resource_search(q = 'name:data', as = 'json')
resource_search(q = 'name:data', as = 'table')
resource_search(q = 'name:data', limit = 2, as = 'table')
resource_search(q=c("description:encoded", "name:No.2"),url='demo.ckan.org')
}
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
% Please edit documentation in R/ckan_fetch.R
\name{ckan_fetch}
\alias{ckan_fetch}
\title{Download a file}
\usage{
ckan_fetch(
  x,
  store = "session",
  path = "file",
  format = NULL,
  key = get_default_key(),
  ...
)
}
\arguments{
\item{x}{URL for the file}

\item{store}{One of session (default) or disk. session stores in R session,
and disk saves the file to disk.}

\item{path}{if \code{store="disk"}, you must give a path to store file to}

\item{format}{Format of the file. Required if format is not detectable
through file URL.}

\item{key}{A CKAN API key (optional, character)}

\item{...}{Arguments passed on to the function used to read the file, if \code{store="session"}. See details.}
}
\description{
Download a file
}
\details{
The \code{...} argument can be used to pass additional arguments to the
function that actually reads the file. It is only used when \code{store="session"}
and if the file is \emph{not} a ZIP file.

This list shows what function is used to read what file format, so that you can
see what additional arguments are available:
\itemize{
\item csv: \link[utils]{read.csv}
\item xls, xlsx: \link[readxl]{read_excel}
\item xml: \link[xml2]{read_xml}
\item html: \link[xml2]{read_html}
\item json: \link[jsonlite]{fromJSON}
\item shp, geojson: \link[sf]{st_read}
\item txt: \link[utils]{read.table}
}
}
\examples{
\dontrun{
# CSV file
ckanr_setup("http://datamx.io")
res <- resource_show(id = "6145a539-cbde-4b0d-a3d3-d1a5eb013f5c",
as = "table")
head(ckan_fetch(res$url))
ckan_fetch(res$url, "disk", "myfile.csv")

# CSV file, format not available
ckanr_setup("https://ckan0.cf.opendata.inter.prod-toronto.ca")
res <- resource_show(id = "c57c3e1c-20e2-470f-bc82-e39a0264be31",
as = "table")
res$url
res$format
head(ckan_fetch(res$url, format = res$format))

# Excel file - requires readxl package
ckanr_setup("http://datamx.io")
res <- resource_show(id = "e883510e-a082-435c-872a-c5b915857ae1",
as = "table")
head(ckan_fetch(res$url))

# Excel file, multiple sheets - requires readxl package
ckanr_setup()
res <- resource_show(id = "ce02a1cf-35f1-41df-91d9-11ed1fdd4186",
as = "table")
x <- ckan_fetch(res$url)
names(x)
head(x[["Mayor - Maire"]])

# XML file - requires xml2 package
# ckanr_setup("http://data.ottawa.ca")
# res <- resource_show(id = "380061c1-6c46-4da6-a01b-7ab0f49a881e",
# as = "table")
# ckan_fetch(res$url)

# HTML file - requires xml2 package
ckanr_setup("http://open.canada.ca/data/en")
res <- resource_show(id = "80321bac-4283-487c-93bd-c65acaa660f5",
as = "table")
ckan_fetch(res$url)
library("xml2")
xml_text(xml_find_first(xml_children(ckan_fetch(res$url))[[1]], "title"))

# JSON file, by default reads in to a data.frame for ease of use
ckanr_setup("http://data.surrey.ca")
res <- resource_show(id = "8d07c662-800d-4977-9e3e-5a3d2d1e99ab",
as = "table")
head(ckan_fetch(res$url))

# SHP file (spatial data, ESRI format) - requires sf package
ckanr_setup("https://ckan0.cf.opendata.inter.prod-toronto.ca")
res <- resource_show(id = "27362290-8bbf-434b-a9de-325a6c2ef923",
as = "table")
x <- ckan_fetch(res$url)
class(x)
plot(x[, c("AREA_NAME", "geometry")])

# GeoJSON file - requires sf package
ckanr_setup("http://datamx.io")
res <- resource_show(id = "b1cd35b7-479e-4fa0-86e9-e897d3c617e6",
as = "table")
x <- ckan_fetch(res$url)
class(x)
plot(x[, c("mun_name", "geometry")])

# ZIP file - packages required depends on contents
ckanr_setup("https://ckan0.cf.opendata.inter.prod-toronto.ca")
res <- resource_show(id = "bb21e1b8-a466-41c6-8bc3-3c362cb1ed55",
as = "table")
x <- ckan_fetch(res$url)
names(x)
head(x[["ChickenpoxAgegroups2017.csv"]])

# TXT file
ckanr_setup("https://ckan0.cf.opendata.inter.prod-toronto.ca")
res <- resource_show(id = "e4211f49-611f-438c-a444-aaa7f3f84117",
as = "table")
x <- ckan_fetch(res$url)
head(x)

# TXT file, semicolon used as separator
ckanr_setup("https://data.coat.no")
res <- resource_show(id = "384fe537-e0bd-4e57-8a0d-420b7a745196",
as = "table")
x <- ckan_fetch(res$url, sep = ";")
head(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_followee_count.R
\name{user_followee_count}
\alias{user_followee_count}
\title{Return a a user's follower count}
\usage{
user_followee_count(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) User identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Return a a user's follower count
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# list package activity
user_followee_count('sckottie')

# input a ckan_user object
(x <- user_show('sckottie'))
user_followee_count(x)

# output different data formats
user_followee_count(x, as = "table")
user_followee_count(x, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_delete.R
\name{package_delete}
\alias{package_delete}
\title{Delete a package}
\usage{
package_delete(id, url = get_default_url(), key = get_default_key(), ...)
}
\arguments{
\item{id}{(character) The id of the package. Required.}

\item{url}{Base url to use. Default: https://data.ontario.ca
See also \code{\link[=ckanr_setup]{ckanr_setup()}} and \code{\link[=get_default_url]{get_default_url()}}}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST} (optional)}
}
\description{
Delete a package
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org", key = getOption("ckan_demo_key"))

# create a package
(res <- package_create("lions-bears-tigers"))

# show the package
package_show(res)

# delete the package
package_delete(res)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_delete.R
\name{resource_delete}
\alias{resource_delete}
\title{Delete a resource.}
\usage{
resource_delete(id, url = get_default_url(), key = get_default_key(), ...)
}
\arguments{
\item{id}{(character) Resource identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca
See also \code{\link[=ckanr_setup]{ckanr_setup()}} and \code{\link[=get_default_url]{get_default_url()}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST} (optional)}
}
\description{
Delete a resource.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = Sys.getenv("CKAN_DEMO_KEY"))

# create a package
(res <- package_create("yellow9"))

# then create a resource
file <- system.file("examples", "actinidiaceae.csv", package = "ckanr")
(xx <- resource_create(res,
                       description = "my resource",
                       name = "bears",
                       upload = file,
                       rcurl = "http://google.com"
))

# delete the resource
resource_delete(xx)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dashboard_count.R
\name{dashboard_count}
\alias{dashboard_count}
\title{Number of new activities of an authorized user}
\usage{
dashboard_count(
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Number of new activities of an authorized user
}
\details{
Important: Activities from the user herself are not counted by this
function even though they appear in the dashboard (users don't want to be
notified about things they did themselves).
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# count
dashboard_count()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dashboard_activity_list.R
\name{dashboard_activity_list}
\alias{dashboard_activity_list}
\title{Authorized user's dashboard activity stream}
\usage{
dashboard_activity_list(
  limit = 31,
  offset = 0,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{limit}{(integer) The maximum number of activities to return (optional).
Default: 31}

\item{offset}{(integer) Where to start getting activity items from (optional).
Default: 0}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Authorized user's dashboard activity stream
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# get activity
(res <- dashboard_activity_list())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_create.R
\name{group_create}
\alias{group_create}
\title{Create a group}
\usage{
group_create(
  name = NULL,
  id = NULL,
  title = NULL,
  description = NULL,
  image_url = NULL,
  type = NULL,
  state = "active",
  approval_status = NULL,
  extras = NULL,
  packages = NULL,
  groups = NULL,
  users = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{name}{(character) the name of the new dataset, must be between 2 and
100 characters long and contain only lowercase alphanumeric characters,
\itemize{
\item and _, e.g. 'warandpeace'
}}

\item{id}{(character) The id of the group (optional)}

\item{title}{(character) The title of the dataset (optional, default:
same as name)}

\item{description}{(character) The description of the group (optional)}

\item{image_url}{(character) The URL to an image to be displayed on the
group's page (optional)}

\item{type}{(character) The type of the dataset (optional), IDatasetForm
plugins associate themselves with different dataset types and provide custom
dataset handling behaviour for these types}

\item{state}{(character) The current state of the dataset, e.g. 'active' or
'deleted', only active datasets show up in search results and other lists
of datasets, this parameter will be ignored if you are not authorized to
change the state of the dataset (optional, default: 'active')}

\item{approval_status}{(character) Approval status (optional)}

\item{extras}{(list of dataset extra dictionaries) The dataset's extras
(optional), extras are arbitrary (key: value) metadata items that can be
added to datasets, each extra dictionary should have keys 'key' (a string),
'value' (a string)}

\item{packages}{(data.frame) The datasets (packages) that belong
to the group, a data.frame, where each row has column 'name' (string, the id
or name of the dataset) and optionally 'title' (string, the title of
the dataset)}

\item{groups}{(data.frame) The groups to which the dataset
belongs (optional), each data.frame row should have one or more of the
following columns which identify an existing group: 'id' (the id of the group,
string), or 'name' (the name of the group, string), to see which groups
exist call \code{\link[=group_list]{group_list()}}}

\item{users}{(list of dictionaries) The users that belong to the group,
a list of dictionaries each with key 'name' (string, the id or name of the
user) and optionally 'capacity' (string, the capacity in which the user is
a member of the group)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Create a group
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org", key = getOption("ckan_demo_key"))

# create a group
(res <- group_create("fruitloops2", description="A group about fruitloops"))
res$users
res$num_followers
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckanr_settings.R
\name{ckanr_settings}
\alias{ckanr_settings}
\alias{get_default_url}
\alias{get_default_key}
\alias{get_test_url}
\alias{get_test_key}
\alias{get_test_did}
\alias{get_test_rid}
\alias{get_test_gid}
\alias{get_test_oid}
\alias{get_test_behaviour}
\title{Get or set ckanr CKAN settings}
\usage{
ckanr_settings()

get_default_url()

get_default_key()

get_test_url()

get_test_key()

get_test_did()

get_test_rid()

get_test_gid()

get_test_oid()

get_test_behaviour()
}
\value{
\code{ckanr_settings} prints your base url, API key (if used), and
optional test server settings (URL, API key, a dataset ID and a resource ID).
\code{ckanr_setup} sets your production and test settings, while
\verb{get_test_*} get each of those respective settings.
\code{test_behaviour} indicates whether the CKANR test suite will skip
("SKIP") or fail ("FAIL") writing tests in case the configured test
CKAN settings don't work.
}
\description{
Get or set ckanr CKAN settings
}
\examples{
ckanr_settings()
}
\seealso{
\code{\link[=ckanr_setup]{ckanr_setup()}},
\code{\link[=get_default_url]{get_default_url()}}, \code{\link[=get_default_key]{get_default_key()}}, \code{\link[=get_test_url]{get_test_url()}},
\code{\link[=get_test_key]{get_test_key()}}, \code{\link[=get_test_did]{get_test_did()}}, \code{\link[=get_test_rid]{get_test_rid()}},
\code{\link[=get_test_gid]{get_test_gid()}}, \code{get_test_oid}, \code{get_test_behaviour}
}
\concept{ckanr settings}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_update.R
\name{group_update}
\alias{group_update}
\title{Update a group}
\usage{
group_update(
  x,
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{x}{(list) A list with key-value pairs}

\item{id}{(character) Package identifier}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Update a group
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# First, create a group
grp <- group_create("water-bears2")
group_show(grp)

## update just chosen things
# Make some changes
x <- list(description = "A group about water bears and people that love them")

# Then update the packge
group_update(x, id = grp)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_delete.R
\name{group_delete}
\alias{group_delete}
\title{Delete a group}
\usage{
group_delete(id, url = get_default_url(), key = get_default_key(), ...)
}
\arguments{
\item{id}{(character) The id of the group. Required.}

\item{url}{Base url to use. Default: https://data.ontario.ca
See also \code{\link[=ckanr_setup]{ckanr_setup()}} and \code{\link[=get_default_url]{get_default_url()}}}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST} (optional)}
}
\description{
Delete a group
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org", key = getOption("ckan_demo_key"))

# create a group
(res <- group_create("lions", description="A group about lions"))

# show the group
group_show(res$id)

# delete the group
group_delete(res)
## or with it's id
# group_delete(res$id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds_search.R
\name{ds_search}
\alias{ds_search}
\title{Datastore - search or get a dataset from CKRAN datastore}
\usage{
ds_search(
  resource_id = NULL,
  filters = NULL,
  q = NULL,
  plain = NULL,
  language = NULL,
  fields = NULL,
  offset = NULL,
  limit = NULL,
  sort = NULL,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{resource_id}{(character) id or alias of the resource to be searched
against}

\item{filters}{(character) matching conditions to select, e.g
\verb{\{"key1": "a", "key2": "b"\}} (optional)}

\item{q}{(character) full text query (optional)}

\item{plain}{(character) treat as plain text query (optional, default:
\code{TRUE})}

\item{language}{(character) language of the full text query (optional,
default: english)}

\item{fields}{(character) fields to return (optional, default: all fields
in original order)}

\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 100)}

\item{sort}{Field to sort on. You can specify ascending (e.g., score desc) or
descending (e.g., score asc), sort by two fields (e.g., score desc,
price asc), or sort by a function (e.g., sum(x_f, y_f) desc, which sorts
by the sum of x_f and y_f in a descending order). (optional)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Datastore - search or get a dataset from CKRAN datastore
}
\details{
From the help for this method "The datastore_search action allows
you to search data in a resource. DataStore resources that belong to
private CKAN resource can only be read by you if you have access to the
CKAN resource and send the appropriate authorization."

Setting \code{plain=FALSE} enables the entire PostgreSQL \emph{full text search query
language}. A listing of all available resources can be found at the alias
\emph{table_metadata} full text search query language:
http://www.postgresql.org/docs/9.1/static/datatype-textsearch.html#DATATYPE-TSQUERY
}
\examples{
\dontrun{
ckanr_setup(url = 'https://data.nhm.ac.uk/')

ds_search(resource_id = '8f0784a6-82dd-44e7-b105-6194e046eb8d')
ds_search(resource_id = '8f0784a6-82dd-44e7-b105-6194e046eb8d',
  as = "table")
ds_search(resource_id = '8f0784a6-82dd-44e7-b105-6194e046eb8d',
  as = "json")

ds_search(resource_id = '8f0784a6-82dd-44e7-b105-6194e046eb8d', limit = 1,
  as = "table")
ds_search(resource_id = '8f0784a6-82dd-44e7-b105-6194e046eb8d', q = "a*")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tag_create.R
\name{tag_create}
\alias{tag_create}
\title{Create a tag}
\usage{
tag_create(
  name,
  vocabulary_id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{name}{(character) The name for the new tag, a string between 2 and 100
characters long containing only alphanumeric characters and -, _ and .,
e.g. 'Jazz'}

\item{vocabulary_id}{(character) The id of the vocabulary that the new
tag should be added to, e.g. the id of vocabulary 'Genre'}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
IMPORTANT: You must be a sysadmin to create vocabulary tags.
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
  key = Sys.getenv("CKAN_DEMO_KEY"))
tag_create(name = "TestTag1", vocabulary_id = "Testing1")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_show.R
\name{resource_show}
\alias{resource_show}
\title{Show a resource.}
\usage{
resource_show(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) Resource identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Show a resource.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/",
key = Sys.getenv("CKAN_DEMO_KEY"))

# create a package
(res <- package_create("yellow7"))

# then create a resource
file <- system.file("examples", "actinidiaceae.csv", package = "ckanr")
(xx <- resource_create(package_id = res$id,
                       description = "my resource",
                       name = "bears",
                       upload = file,
                       rcurl = "http://google.com"
))

# show the resource
resource_show(xx$id)


# eg. from the NHM CKAN store
resource_show(id = "05ff2255-c38a-40c9-b657-4ccb55ab2feb",
              url = "http://data.nhm.ac.uk")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_patch.R
\name{package_patch}
\alias{package_patch}
\title{Update a package's metadata}
\usage{
package_patch(
  x,
  id = NULL,
  extras = NULL,
  http_method = "GET",
  key = get_default_key(),
  url = get_default_url(),
  as = "list",
  ...
)
}
\arguments{
\item{x}{(list) A list with key-value pairs}

\item{id}{(character) Resource ID to update (optional, required if
x does not have an "id" field)}

\item{extras}{(character vector) - the dataset's extras
(optional), extras are arbitrary (key: value) metadata items that can be
added to datasets, each extra dictionary should have keys 'key' (a string),
'value' (a string)}

\item{http_method}{(character) which HTTP method (verb) to use; one of
"GET" or "POST". Default: "GET"}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Update a package's metadata
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org", key = getOption("ckan_demo_key"))

# Create a package
(res <- package_create("hello-world13", author="Jane Doe"))

# Get a resource
res <- package_show(res$id)
res$title

# patch
package_patch(res, extras = list(list(key = "foo", value = "bar")))
unclass(package_show(res))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds_search_sql.R
\name{ds_search_sql}
\alias{ds_search_sql}
\title{Datastore - search or get a dataset from CKRAN datastore}
\usage{
ds_search_sql(
  sql,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{sql}{(character) A single SQL select statement. (required)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Datastore - search or get a dataset from CKRAN datastore
}
\examples{
\dontrun{
url <- 'https://demo.ckan.org/'
sql <- 'SELECT * from "f4129802-22aa-4437-b9f9-8a8f3b7b2a53" LIMIT 2'
ds_search_sql(sql, url = url, as = "table")
sql2 <- 'SELECT "Species","Genus","Family" from "f4129802-22aa-4437-b9f9-8a8f3b7b2a53" LIMIT 2'
ds_search_sql(sql2, url = url, as = "table")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_list_current.R
\name{package_list_current}
\alias{package_list_current}
\title{List current packages with resources.}
\usage{
package_list_current(
  offset = 0,
  limit = 31,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{offset}{(numeric) Where to start getting activity items from
(optional, default: 0)}

\item{limit}{(numeric) The maximum number of activities to return
(optional, default: 31)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
List current packages with resources.
}
\examples{
\dontrun{
package_list_current()
package_list_current(as = 'json')
package_list_current(as = 'table')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/related_delete.R
\name{related_delete}
\alias{related_delete}
\title{Delete a related item.}
\usage{
related_delete(id, url = get_default_url(), key = get_default_key(), ...)
}
\arguments{
\item{id}{(character) Resource identifier.}

\item{url}{Base url to use. Default: https://data.ontario.ca See
also \code{\link[=ckanr_setup]{ckanr_setup()}} and \code{\link[=get_default_url]{get_default_url()}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST} (optional)}
}
\description{
Delete a related item.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://demo.ckan.org/", key = getOption("ckan_demo_key"))

# create a package and a related item
res <- package_create("hello-venus2") \%>\%
   related_create(title = "my resource",
                  type = "visualization")

# show the related item
related_delete(res)
## or with id itself:
## related_delete(res$id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organization_purge.R
\name{organization_purge}
\alias{organization_purge}
\title{Purge an organization}
\usage{
organization_purge(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) name or id of the organization}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\value{
an empty list on success
}
\description{
IMPORTANT: You must be a sysadmin to purge an organization. Once an
organization is purged, it is permanently removed from the system.
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org", key=getOption("ckan_demo_key"))

# create an organization
(res <- organization_create("foobar", title = "Foo bars",
  description = "love foo bars"))

# delete the organization just created
res$id
organization_delete(id = res$id)

# purge the organization just deleted
res$id
organization_purge(id = res$id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_delete.R
\name{user_delete}
\alias{user_delete}
\title{Delete a user.}
\usage{
user_delete(
  id,
  url = get_default_url(),
  key = get_default_key(),
  as = "list",
  ...
)
}
\arguments{
\item{id}{(character) the id of the new user (required)}

\item{url}{Base url to use. Default: https://data.ontario.ca/ See
also \code{\link{ckanr_setup}} and \code{\link{get_default_url}}.}

\item{key}{A privileged CKAN API key, Default: your key set with
\code{\link{ckanr_setup}}}

\item{as}{(character) One of list (default), table, or json. Parsing with
table option uses \code{jsonlite::fromJSON(..., simplifyDataFrame = TRUE)},
which attempts to parse data to data.frame's when possible, so the result
can vary from a vector, list or data.frame. (required)}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}} (optional)}
}
\description{
Delete a user.
}
\examples{
\dontrun{
# Setup
ckanr_setup(url = "https://data-demo.dpaw.wa.gov.au",
key = "824e7c50-9577-4bfa-bf32-246ebed1a8a2")

# create a user
res <- user_delete(name = 'stacy', email = "stacy@aaaaa.com",
password = "helloworld")

# then, delete a user
user_delete(id = "stacy")
}
}
\references{
http://docs.ckan.org/en/latest/api/index.html#ckan.logic.action.delete.user_delete
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ckan_resource.R
\name{as.ckan_resource}
\alias{as.ckan_resource}
\alias{is.ckan_resource}
\title{ckan_resource class helpers}
\usage{
as.ckan_resource(x, ...)

is.ckan_resource(x)
}
\arguments{
\item{x}{Variety of things, character, list, or ckan_package class object}

\item{...}{Further args passed on to \code{\link[=resource_show]{resource_show()}} if character given}
}
\description{
ckan_resource class helpers
}
\examples{
\dontrun{
ckanr_setup(url = "https://demo.ckan.org/",
key = getOption("ckan_demo_key"))

(resrcs <- resource_search(q = 'name:data'))
resrcs$results
resrcs$results[[3]]

# create item class from only an item ID
as.ckan_resource(resrcs$results[[3]]$id)

# gives back itself
(x <- as.ckan_resource(resrcs$results[[3]]$id))
as.ckan_resource(x)
}
}
